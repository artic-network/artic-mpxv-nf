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
    
    shell:
    """
    workflow-glue check_sample_sheet sample_sheet.txt samples.txt
    """
}

process runArtic {
    label "artic"
    cpus params.artic_threads
    errorStrategy { task.exitStatus in ((130..145) + 104) ? "retry" : "terminate" }
    maxRetries 3

    input:
        tuple val(meta), path(fastq_file), path(fastq_stats), path(scheme_dir)
        val models_ok
    output:
        path "${meta.alias}.consensus.fasta", emit: consensus
        path "${meta.alias}.minion.log.txt", emit: artic_log
        path "${meta.alias}.amplicon_depths.tsv", emit: amplicon_depths, optional: true
        path "${meta.alias}.sorted.bam*", emit: raw_bam, optional: true
        tuple(
            val(meta.alias),
            path("${meta.alias}.normalised.named.vcf.gz"),
            path("${meta.alias}.normalised.named.vcf.gz.tbi"),
            emit: pass_vcf, optional: true)
        tuple(
            val(meta.alias),
            path("${meta.alias}.primertrimmed.rg.sorted.bam"),
            path("${meta.alias}.primertrimmed.rg.sorted.bam.bai"),
            emit: primertrimmed_bam, optional: true)
        tuple(
            val(meta.alias),
            path("${meta.alias}.trimmed.rg.sorted.bam"),
            path("${meta.alias}.trimmed.rg.sorted.bam.bai"),
            emit: trimmed_bam, optional: true)
    script:
    // we use `params.override_basecaller_cfg` if present; otherwise use
    // `meta.basecall_models[0]` (there should only be one value in the list because
    // we're running ingress with `allow_multiple_basecall_models: false`; note that
    // `[0]` on an empty list returns `null`)
    if (!meta.basecall_models[0] && !params.override_model) {
        error "Found no basecall model information in the input data for " + \
            "sample '$meta.alias' and no model provided with the `--override_model` " + \
            "parameter. Please rerun with the `--override_model` parameter."
    }

    if (params.override_model) {
        model_str = "--model ${params.override_model}"
    } else {
        model_str = ""
    }

    if (!params.custom_scheme) {
        """
        artic guppyplex --skip-quality-check \
            --min-length ${params._min_len} --max-length ${params._max_len} \
            --directory . --prefix ${meta.alias}

        artic minion --normalise ${params.normalise} --threads ${task.cpus} \
            --read-file ${meta.alias}_..fastq \
            --scheme-name ${params.parsed_scheme_name} \
            --scheme-length ${params.parsed_scheme_length} \
            --scheme-version ${params.parsed_scheme_version} \
            --scheme-directory ${params.store_dir}/primer-schemes/ \
            ${model_str} \
            ${meta.alias}

        zcat ${meta.alias}.normalised.vcf.gz | sed 's/SAMPLE/${meta.alias}/' | bgzip > ${meta.alias}.normalised.named.vcf.gz
        bcftools index -t ${meta.alias}.normalised.named.vcf.gz
        """        
    } else {
        """
        artic guppyplex --skip-quality-check \
            --min-length ${params._min_len} --max-length ${params._max_len} \
            --directory . --prefix ${meta.alias}

        artic minion --normalise ${params.normalise} --threads ${task.cpus} \
            --read-file ${meta.alias}_..fastq \
            --bed ${scheme_dir}/primer.bed \
            --ref ${scheme_dir}/reference.fasta \
            ${model_str} \
            ${meta.alias}

        zcat ${meta.alias}.normalised.vcf.gz | sed 's/SAMPLE/${meta.alias}/' | bgzip > ${meta.alias}.normalised.named.vcf.gz
        bcftools index -t ${meta.alias}.normalised.named.vcf.gz
        """
    }
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
    run_clair3.sh --version | sed 's/ /,/' >> versions.txt
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

process squirrel {
    label "squirrel"
    cpus params.squirrel_threads

    input:
        path "all_consensus.fasta"
    output:
        path "squirrel/all_consensus*.aln.fasta", emit: alignment
        path "squirrel", emit: all
        path "squirrel.version", emit: version

    script:
    """
    export XDG_CACHE_HOME=\$PWD/.cache
    squirrel --version 2>&1 | sed 's/: /,/' > squirrel.version
    squirrel "all_consensus.fasta" -o squirrel --no-mask --seq-qc --outfile all_consensus.aln.fasta --tempdir squirrel_tmp -t ${task.cpus} --clade ${params.clade} $params._squirrel_options
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output_results {
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

process get_models {
    label "artic"
    cpus 1
    output:
        val "models_ok"
    
    script:

    """
    artic_get_models
    """
}

// workflow module
workflow pipeline {
    take:
        samples
        scheme_dir
    main:
        software_versions = getVersions()
        // workflow_params = getParams()
        // combined_genotype_summary = Channel.empty()

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
            if (workflow.session.config.conda) {
                get_models()
                ch_models_ok = get_models.out
            } else {
                ch_models_ok = Channel.value("models_ok")
            }

            ch_artic_in = samples.map { meta, reads, stats ->
                [
                    meta,
                    reads,
                    stats,
                    scheme_dir
                ]
            }

            artic = runArtic(ch_artic_in, ch_models_ok)
            // all_depth = combineDepth(artic.depth_stats.collect())
            // collate consensus and variants
            all_consensus = allConsensus(artic.consensus.collect())
            // genotype summary
            genotype_summary = Channel.fromPath("$projectDir/data/OPTIONAL_FILE")

            // squirrel
            if (!params.skip_squirrel) {
                squirrel(all_consensus[0])
                software_versions = software_versions.mix(squirrel.out.version)
                
                results = all_consensus[0].concat(
                    all_consensus[1],
                    artic.primertrimmed_bam.flatMap { it -> [ it[1], it[2] ] },
                    artic.pass_vcf.flatMap { it -> [ it[1], it[2] ] },
                    artic.artic_log,
                    artic.consensus,
                    artic.amplicon_depths,
                    artic.raw_bam.flatMap { it -> [ it[0], it[1] ] },
                    squirrel.out.alignment.flatMap { it -> [ it ] },
                    squirrel.out.all.flatMap { it -> [ it ] },
                    )
            } else {
                results = all_consensus[0].concat(
                    all_consensus[1],
                    artic.primertrimmed_bam.flatMap { it -> [ it[1], it[2] ] },
                    artic.pass_vcf.flatMap { it -> [ it[1], it[2] ] },
                    artic.artic_log,
                    artic.consensus,
                    artic.amplicon_depths,
                    artic.raw_bam.flatMap { it -> [ it[0], it[1] ] }
                )
            }
            
            }
    emit:
        results            
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)


workflow {

    Pinguscript.ping_start(nextflow, workflow, params)

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";

    if (!params.custom_scheme){
      
      params._bed = false
      params._ref = false
      scheme_dir = []

      if (!params.scheme_version) {
          error "${c_purple}EXITING: --scheme_version parameter must be specified.${c_reset}"
      }

      if (!params.min_len) {
          params.remove('min_len')
          if (params.scheme_version.startsWith("yale-mpox") || params.scheme_version.startsWith("rigshospitalet") || params.scheme_version.startsWith("artic-mpox") || params.scheme_version.startsWith("bccdc-mpox") || params.scheme_version.startsWith("artic-inrb-mpox")) {
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
            if (params.scheme_version.startsWith("yale-mpox") || params.scheme_version.startsWith("rigshospitalet") || params.scheme_version.startsWith("artic-mpox") || params.scheme_version.startsWith("bccdc-mpox") || params.scheme_version.startsWith("artic-inrb-mpox")) {
              params._max_len = 3000
          } else {
                params._max_len = 2500
          }
      } else {
          params._max_len = params.max_len
          params.remove('max_len')
      }
      
      scheme_splits = params.scheme_version.split("/")
      if (scheme_splits.size() != 3) {
        error "${c_purple}EXITING: --scheme_version parameter must be in the format 'scheme_name/scheme_length/scheme_version' (e.g. 'artic-inrb-mpox/2500/v1.0.0').${c_reset}"
      }

      params.parsed_scheme_name = scheme_splits[0]
      params.parsed_scheme_length = scheme_splits[1]
      params.parsed_scheme_version = scheme_splits[2]
      println """${c_purple}Using primer scheme: ${params.parsed_scheme_name} (${params.parsed_scheme_length} bp, version ${params.parsed_scheme_version})${c_reset}"""

    } else {
      //custom scheme path defined
      log.info """${c_purple}Custom primer scheme selected: ${params.custom_scheme} (WARNING: We do not validate your scheme - use at your own risk!)${c_reset}"""
      //check path for required files
      params._bed = file("""${params.custom_scheme}/primer.bed""", type:'file', checkIfExists:true)
      params._ref = file("""${params.custom_scheme}/reference.fasta""", type:'file', checkIfExists:true)

      scheme_dir = file(params.custom_scheme, type: 'dir', checkIfExists: true)

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

    }

    if (!params.max_softclip_length) {
        params.remove('max_softclip_length')
        params._max_softclip_length = 0
    }
    else{
        params._max_softclip_length = params.max_softclip_length
        params.remove('max_softclip_length')
    }
    
    // Squirrel options
    if (params.squirrel_options == false){
        params.remove('squirrel_options')
        params._squirrel_options = ''
    } else {
        params._squirrel_options = params.squirrel_options
        params.remove('squirrel_options')
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
    
    pipeline(samples, scheme_dir)
    output_results(pipeline.out.results.toList())
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
