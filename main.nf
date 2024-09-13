#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
import java.text.SimpleDateFormat

nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/ingress'

process checkSampleSheet {
    label "artic"
    cpus 1
    input:
        file "sample_sheet.txt"
    output:
        file "samples.txt"
    """
    workflow-glue check_sample_sheet sample_sheet.txt samples.txt
    """
}


process runArtic {
    label "artic"
    cpus params.artic_threads
    errorStrategy 'retry'
    maxRetries 3

    input:
        tuple val(meta), path(fastq_file), path(fastq_stats)
        path bed
        path ref
    output:
        path "${meta.alias}.consensus.fasta", emit: consensus
        path "${meta.alias}.pass.named.stats", emit: vcf_stats
        path "${meta.alias}.artic.log.txt", emit: artic_log
        path "${meta.alias}.amplicon_depths.tsv", emit: amplicon_depths
        path "${meta.alias}.sorted.bam*", emit: raw_bam
        tuple(
            val(meta.alias),
            path("${meta.alias}.normalised.named.vcf.gz"),
            path("${meta.alias}.normalised.named.vcf.gz.tbi"),
            emit: pass_vcf)
        tuple(
            val(meta.alias),
            path("${meta.alias}.primertrimmed.rg.sorted.bam"),
            path("${meta.alias}.primertrimmed.rg.sorted.bam.bai"),
            emit: primertrimmed_bam)
        tuple(
            val(meta.alias),
            path("${meta.alias}.trimmed.rg.sorted.bam"),
            path("${meta.alias}.trimmed.rg.sorted.bam.bai"),
            emit: trimmed_bam)
    script:
    // we use `params.override_basecaller_cfg` if present; otherwise use
    // `meta.basecall_models[0]` (there should only be one value in the list because
    // we're running ingress with `allow_multiple_basecall_models: false`; note that
    // `[0]` on an empty list returns `null`)
    String basecall_model = params.override_basecaller_cfg ?: meta.basecall_models[0]
    if (!basecall_model) {
        error "Found no basecall model information in the input data for " + \
            "sample '$meta.alias'. Please provide it with the " + \
            "`--override_basecaller_cfg` parameter."
    }

    """
    run_artic.sh \
        ${meta.alias} ${fastq_file} ${params._min_len} ${params._max_len} \
        ${basecall_model}:consensus  ${bed} ${ref} \
        ${task.cpus} ${params._max_softclip_length} ${params.normalise} ${params.min_reads} \
        > ${meta.alias}.artic.log.txt 2>&1
    bcftools stats ${meta.alias}.normalised.named.vcf.gz > ${meta.alias}.pass.named.stats
    """
}


process combineDepth {
  label "artic"
  cpus 1
  input:
    path "depth_stats/*"
  output:
    file "all_depth.txt"
  script:
  """
    header_file=`ls depth_stats/* | head -1`
    head -1 \${header_file} > all_depth.txt
    cat depth_stats/* | grep -v depth_fwd >> all_depth.txt
  """
}


process genotypeSummary {
    // Produce a genotype summary spreadsheet
    label "artic"
    cpus 1
    input:
        tuple val(alias), file(vcf), file(tbi), file(bam), file(bam_index)
        file "reference.vcf"
    output:
        file "*genotype.csv"
    script:
        def lab_id = params.lab_id ? "--lab_id ${params.lab_id}" : ""
        def testkit = params.testkit ? "--testkit ${params.testkit}" : ""
    """
    workflow-glue genotype_summary \
        -b $bam \
        -v $vcf \
        -d reference.vcf \
        --sample $alias \
        $lab_id \
        $testkit \
        -o ${csvName}.genotype.csv
    """
}


process combineGenotypeSummaries {
    label "artic"
    cpus 1
    input:
        file "summary_*.csv"
    output:
        file "genotype_summary.csv"
    """
    workflow-glue combine_genotype_summaries -g *.csv -o genotype_summary.csv
    """
}


process getVersions {
    label "artic"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    medaka --version | sed 's/ /,/' >> versions.txt
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    bcftools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    artic --version | sed 's/ /,/' >> versions.txt
    """
}


process getParams {
    label "wf_common"
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}


process report_no_data {
    label "artic"
    cpus 1
    input:
        path "versions/*"
        val error
        path "params.json"
    output:
        path "wf-artic-*.html"
        path "*.json", optional: true
    script:
    // when genotype_variants is false the channel contains a mock file
    def report_name = "wf-artic-report.html"
    def error_message = error
    """
    workflow-glue report_error \
        --output $report_name \
        --revision $workflow.revision --params params.json --commit $workflow.commitId \
        --versions versions --error_message \"$error_message\"
    """
}


process allConsensus {
    label "artic"
    cpus 1
    input:
        file "*"
    output:
        file "all_consensus.fasta"
        file "consensus_status.txt"
    """
    ls *.consensus.fasta | xargs cat > all_consensus.fasta
    grep "^>" all_consensus.fasta \
        | awk 'BEGIN{OFS="\\t"; print "sample\\tpass"}{print substr(\$1, 2), \$2!="Artic-Fail"}' \
        >> consensus_status.txt
    """
}


process allVariants {
    label "artic"
    cpus 1
    input:
        tuple val(alias), file(vcfs), file(tbis)
        file reference
    output:
        tuple file("all_variants.vcf.gz"), file("all_variants.vcf.gz.tbi")
    """
    for vcf in \$(ls *.vcf.gz)
    do
        bcftools norm -c s -O z --fasta-ref $reference \$vcf > norm.\$vcf
        bcftools index -t norm.\$vcf
    done
    if [[ \$(ls norm.*.vcf.gz | wc -l) == "1" ]]; then
        mv norm.*.vcf.gz all_variants.vcf.gz
        mv norm.*.vcf.gz.tbi all_variants.vcf.gz.tbi
    else
        bcftools merge -o all_variants.vcf.gz -O z norm.*.vcf.gz
        bcftools index -t all_variants.vcf.gz
    fi
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "artic"

    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        file fname
    output:
        file fname
    """
    echo "Writing output files"
    """
}

process get_bed_ref {
    label "artic"
    cpus 1
    input:
        path scheme_dir
        val scheme_name
        val scheme_version
    output:
        path "scheme.bed", emit: bed
        path "reference.fasta", emit: ref

    """
    cp ${scheme_name}/${scheme_version}/primer.bed scheme.bed
    cp ${scheme_name}/${scheme_version}/reference.fasta reference.fasta
    """
}

// workflow module
workflow pipeline {
    take:
        samples
        // scheme_directory
        scheme_dir
        scheme_name
        scheme_version
        reference
        primers
    main:
        software_versions = getVersions()
        // workflow_params = getParams()
        combined_genotype_summary = Channel.empty()

        if ((samples.getClass() == String) && (samples.startsWith("Error"))){
            samples = channel.of(samples)
            html_doc = report_no_data(
                software_versions.collect(),
                samples,
                workflow_params)
            results = html_doc[0].concat(html_doc[1])
        } else {
            // remove samples that only appeared in the sample sheet but didn't have any
            // reads
            samples = samples
            | map { meta, reads, stats ->
                if (!reads) {
                    log.warn "No input data found for sample '$meta.alias'; skipping..."
                } else {
                    [meta, reads, stats]
                }
            }
            artic = runArtic(samples, primers, reference)
            // all_depth = combineDepth(artic.depth_stats.collect())
            // collate consensus and variants
            all_consensus = allConsensus(artic.consensus.collect())
            all_variants = allVariants(
                artic.pass_vcf.toList().transpose().toList(), reference)
            // genotype summary
            genotype_summary = Channel.fromPath("$projectDir/data/OPTIONAL_FILE")

            results = all_consensus[0].concat(
                all_consensus[1],
                all_variants[0].flatten(),
                artic.primertrimmed_bam.flatMap { it -> [ it[1], it[2] ] },
                artic.pass_vcf.flatMap { it -> [ it[1], it[2] ] },
                artic.artic_log,
                artic.consensus,
                artic.amplicon_depths,
                artic.raw_bam
                )
            }
    emit:
        results            
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)

// here we should check if the scheme exists, if not, list schemes and exit




workflow {

    Pinguscript.ping_start(nextflow, workflow, params)

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";

    if (!params.custom_scheme){

      schemes = file(projectDir.resolve("./data/primer_schemes/**bed"), type: 'file', maxdepth: 10)

      valid_scheme_versions = []

      log.info """
      ------------------------------------
      Available Primer Schemes:
      ------------------------------------
      """
      log.info """  Name\t\tVersion"""
      for (scheme in schemes){
        main = scheme.toString().split("primer_schemes/")[1]
        name = main.split("/")[0]
        version = """${main.split("/")[1]}/${main.split("/")[2]}"""
        valid_scheme_versions.add(version)
        log.info """${c_green}  ${name}\t${version}\t${c_reset}"""
      }

      log.info """
      ------------------------------------
      """

      if (params.list_schemes) {
        exit 1
      }



      if (!valid_scheme_versions.any { it == params.scheme_version}) {
          println("`--scheme_version` should be one of: $valid_scheme_versions, for `--scheme_name`: $params.scheme_name")
          exit 1
      }

      if (params.sample && params.detect_samples) {
          println("Select either `--sample` or `--detect_samples`, not both")
          exit 1
      }

      if (!params.min_len) {
          params.remove('min_len')
          if (params.scheme_version.startsWith("yale-mpox") || params.scheme_version.startsWith("rigshospitalet") || params.scheme_version.startsWith("artic-mpox") || params.scheme_version.startsWith("bccdc-mpox")) {
              params._min_len = 500
          } else {
              params._min_len = 500
          }
      } else {
          params._min_len = params.min_len
          params.remove('min_len')
      }
      if (!params.max_len) {
          params.remove('max_len')
            if (params.scheme_version.startsWith("yale-mpox") || params.scheme_version.startsWith("rigshospitalet") || params.scheme_version.startsWith("artic-mpox") || params.scheme_version.startsWith("bccdc-mpox")) {
              params._max_len = 3000
          } else {
                params._max_len = 2500
          }
      } else {
          params._max_len = params.max_len
          params.remove('max_len')
      }
    
      
      scheme_dir_name = "primer_schemes"
      schemes = """./data/${scheme_dir_name}/${params.scheme_name}"""
      scheme_dir = file(projectDir.resolve(schemes), type:'file', checkIfExists:true)
    
      primers_path = """./data/${scheme_dir_name}/${params.scheme_name}/${params.scheme_version}/primer.bed"""
      primers = file(projectDir.resolve(primers_path), type:'file', checkIfExists:true)

      reference_path = """./data/${scheme_dir_name}/${params.scheme_name}/${params.scheme_version}/reference.fasta"""
      reference = file(projectDir.resolve(reference_path),type:'file', checkIfExists:true)

      params._scheme_version = params.scheme_version
      params._scheme_name = params.scheme_name

    } else {
      //custom scheme path defined
      log.info """${c_purple}Custom primer scheme selected: ${params.custom_scheme} (WARNING: We do not validate your scheme - use at your own risk!)${c_reset}"""
      //check path for required files
      primers = file("""${params.custom_scheme}/primer.bed""", type:'file', checkIfExists:true)
      reference = file("""${params.custom_scheme}/reference.fasta""", type:'file', checkIfExists:true)

      // check to make sure min and max length have been set
      if (!params.max_len || !params.min_len) {
          log.info """${c_purple}EXITING: --min_len and --max_len parameters must be specified when using custom schemes.${c_reset}"""
          exit 1
      }

      params._max_len = params.max_len
      params.remove('max_len')

      params._min_len = params.min_len
      params.remove('min_len')

      params._scheme_version = 'None'
      params._scheme_name = params.scheme_name

      scheme_dir =  params.custom_scheme
    }

    if (!params.max_softclip_length) {
        params.remove('max_softclip_length')
        params._max_softclip_length = 0
    }
    else{
        params._max_softclip_length = params.max_softclip_length
        params.remove('max_softclip_length')
    }


    // check fastq dataset and run workflow
    samples = fastq_ingress([
        "input":params.fastq,
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "stats": true,
        "per_read_stats": true,
        "analyse_unclassified":params.analyse_unclassified,
        "allow_multiple_basecall_models":false,
    ])
    
    results = pipeline(samples, scheme_dir, params._scheme_name, params._scheme_version, reference,
        primers)
        | output
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
