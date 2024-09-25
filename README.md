# Artic Network MPXV Analysis

Run the ARTIC fieldbioinformatics workflow on multiplexed MPXV ONT data


## Introduction

The artic-mpxv-nf workflow implements an ARTIC FieldBioinformatics
workflow for the purpose of preparing consensus sequences from MPXV
genomes that have been DNA sequenced using a pooled tiling amplicon strategy.

The workflow consumes a folder containing demultiplexed sequence reads as
prepared by either MinKNOW, Guppy, or Dorado. The workflow needs to know the primer
scheme that has been used during genome amplification and library preparation
e.g. `yale-mpox/v1.0.1` or `erasmus/v1.0.0`. Other parameters can be specified too e.g.
assign sample names to the barcodes or to adjust the length distribution of
acceptable amplicon sequences.


# Credits / Acknowledgements

This pipeline only works due to the ongoing efforts of many people performing the often thankless
job of developing and maintaining bioinformatics software, including but not limited to:
* Minimap2 - Heng Li et al, citation: `Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. doi:10.1093/bioinformatics/bty191`
* Samtools - Heng Li et al, citation: `Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li. GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008`
* Bcftools - Heng Li et al, citation: `Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li. GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008`
* Bwa - Heng Li, et al, citation: `Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v2 [q-bio.GN]. (if you use the BWA-MEM algorithm or the fastmap command, or want to cite the whole BWA package)`
* Muscle - Robert Edgar, et al, citation: `Edgar RC., Muscle5: High-accuracy alignment ensembles enable unbiased assessments of sequence homology and phylogeny. Nature Communications 13.1 (2022): 6968.`
* Longshot - Peter Edge et al, citation: `Edge, P. and Bansal, V., 2019. Longshot enables accurate variant calling in diploid genomes from single-molecule long read sequencing. Nature communications, 10(1), pp.1-10.`
* cyvcf2 - Brent Pederson et al, citation: `Brent S Pedersen, Aaron R Quinlan, cyvcf2: fast, flexible variant analysis with Python, Bioinformatics, Volume 33, Issue 12, June 2017, Pages 1867–1869, https://doi.org/10.1093/bioinformatics/btx057`
* Pysam - Anreas Heger et al, `https://github.com/pysam-developers/pysam`
* Clair3 - Zhenxian Zheng et al, citation: `Zheng, Z., Li, S., Su, J. et al. Symphonizing pileup and full-alignment for deep learning-based long-read variant calling. Nat Comput Sci 2, 797–803 (2022). https://doi.org/10.1038/s43588-022-00387-x`
* Medaka - Oxford Nanopore Technologies, Chris Wright et al, `https://github.com/nanoporetech/medaka`


## Compute requirements

Recommended requirements:

+ CPUs = 4
+ Memory = 8GB

Minimum requirements:

+ CPUs = 2
+ Memory = 4GB

Approximate run time: 5 minutes per sample

ARM processor support: False




## Install and run


These are instructions to install and run the workflow on command line.
You can also access the workflow via the
[EPI2ME Desktop application](https://labs.epi2me.io/downloads/).

The workflow uses [Nextflow](https://www.nextflow.io/) to manage
compute and software resources,
therefore Nextflow will need to be
installed before attempting to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop
or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html)
to provide isolation of the required software.
Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.
This is controlled by the
[`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles)
parameter as exemplified below.

It is not required to clone or download the git repository
in order to run the workflow.
More information on running EPI2ME workflows can
be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow.
This will pull the repository in to the assets folder of
Nextflow and provide a list of all parameters
available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-artic --help
```
To update a workflow to the latest version on the command line use
the following command:
```
nextflow pull epi2me-labs/wf-artic
```

A demo dataset is provided for testing of the workflow.
It can be downloaded and unpacked using the following commands:
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-artic/wf-artic-demo.tar.gz
tar -xzvf wf-artic-demo.tar.gz
```
The workflow can then be run with the downloaded demo data using:
```
nextflow run epi2me-labs/wf-artic \
	--fastq 'wf-artic-demo/fastq' \
	--sample_sheet 'wf-artic-demo/sample_sheet.csv' \
	--scheme_name 'SARS-CoV-2' \
	--scheme_version 'Midnight-ONT/V3' \
	-profile standard
```

For further information about running a workflow on
the command line see https://labs.epi2me.io/wfquickstart/




## Related protocols

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.

The Midnight protocol for sample preparation and sequencing can be found in the [Nanopore community](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/pcr-tiling-of-sars-cov-2-virus-rbk114-and-midnight-rt/v/mrt_9186_v114_revd_19apr2023).




## Input example

<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
This workflow accepts FASTQ files as input.

The FASTQ input parameters for this workflow accept one of three cases: (i) the path to a single FASTQ; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second cases (i and ii), a sample name can be supplied with `--sample`. In the last case (iii), the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`.

```
(i)                     (ii)                 (iii)    
input_reads.fastq   ─── input_directory  ─── input_directory
                        ├── reads0.fastq     ├── barcode01
                        └── reads1.fastq     │   ├── reads0.fastq
                                             │   └── reads1.fastq
                                             ├── barcode02
                                             │   ├── reads0.fastq
                                             │   ├── reads1.fastq
                                             │   └── reads2.fastq
                                             └── barcode03
                                              └── reads0.fastq
```



## Input parameters

### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |


### Primer Scheme Selection

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| scheme_name | string | Primer scheme name. | This should be set to `SARS-CoV-2`, or `spike-seq` or your custom scheme name. This affects the choice of scheme versions you can use. The only scheme versions compatible with `spike-seq` are `ONT/V1` and `ONT/V4.1` | SARS-CoV-2 |
| scheme_version | string | Primer scheme version. | This is the version of the primer scheme to use, more details about primer shemes can be found [here](https://labs.epi2me.io/ont-midnight-scheme-update/). | ARTIC/V3 |
| custom_scheme | string | Path to a custom scheme. | If you have a custom primer scheme you can enter the details here. This must be the full path to the directory containing your appropriately named scheme bed and fasta files; <SCHEME_NAME>.bed and <SCHEME_NAME>.fasta. More details [here](https://labs.epi2me.io/ont-midnight-scheme-update/). |  |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files. |  |  |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | output |


### Reporting Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| report_depth | integer | Min. depth for percentage coverage. (e.g. 89% genome covered at > `report_depth`) |  | 100 |
| report_clade | boolean | Show results of Nextclade analysis in report. |  | True |
| report_coverage | boolean | Show genome coverage traces in report. |  | True |
| report_lineage | boolean | Show results of Pangolin analysis in report. |  | True |
| report_variant_summary | boolean | Show variant information in report. |  | True |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| artic_threads | number | Number of CPU threads to use per artic task. | The total CPU resource used by the workflow is constrained by the executor configuration. | 4 |
| pangolin_threads | number | Number of CPU threads to use per pangolin task. | The total CPU resource used by the workflow is constrained by the executor configuration. | 4 |
| genotype_variants | string | Report genotyping information for scheme's known variants of interest, optionally provide file path as argument. |  |  |
| list_schemes | boolean | List primer schemes and exit without running analysis. |  | False |
| min_len | number | Minimum read length (default: set by scheme). |  |  |
| max_len | number | Maximum read length (default: set by scheme). |  |  |
| max_softclip_length | integer | Remove reads with alignments showing large soft clipping |  |  |
| update_data | boolean | Update Pangolin and Nextclade data at runtime. |  | True |
| pangolin_options | string | Pass options to Pangolin, for example "--analysis-mode fast --min-length 26000". |  |  |
| nextclade_data_tag | string | The tag of the nextclade data packet |  |  |
| normalise | integer | Depth ceiling for depth of coverage normalisation |  | 200 |
| override_basecaller_cfg | string | Override auto-detected basecaller model that processed the signal data; used to select an appropriate Medaka model. | Per default, the workflow tries to determine the basecall model from the input data. This parameter can be used to override the detected value (or to provide a model name if none was found in the inputs). However, users should only do this if they know for certain which model was used as selecting the wrong option might give sub-optimal results. A list of recent models can be found here: https://github.com/nanoporetech/dorado#DNA-models. |  |


### Miscellaneous Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| lab_id | string | Laboratory identifier, used in reporting. |  |  |
| testkit | string | Test kit identifier, used in reporting. |  |  |






## Outputs

Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Workflow report | ./wf-artic-report.html | Report for all samples. | aggregated |
| Consensus sequences | ./all_consensus.fasta | Final consensus sequences for all samples in the analysis. | aggregated |
| Pangolin results | ./lineage_report.csv | Pangolin results for each of the samples in the analysis. | aggregated |
| Nextclade results | ./nextclade.json | Nextclade results for each of the samples in the analysis. | aggregated |
| Coverage data | ./all_depth.txt | Coverage of the reference genome in 20 base windows in all the samples in the analysis. | aggregated |
| Variants | ./{{ alias }}.pass.named.vcf.gz | A VCF file containing high confidence variants in the sample when compared to the reference. | per-sample |
| Variants index | ./{{ alias }}.pass.named.vcf.gz.tbi | An index file for the variants. | per-sample |
| Alignments | ./{{ alias }}.primertrimmed.rg.sorted.bam | A BAM file containing the reads for the sample aligned to the reference. | per-sample |
| Alignments index | ./{{ alias }}.primertrimmed.rg.sorted.bam.bai | An index file for the alignments. | per-sample |




## Pipeline overview

The pipeline is largely a wrapper around the [Artic Network](https://artic.network/) [Field Bioinformatics](https://github.com/artic-network/fieldbioinformatics) analysis package.

### 1. Concatenates input files and generate per read stats.

The [fastcat/bamstats](https://github.com/epi2me-labs/fastcat) tool is used to concatenate multifile samples to be processed by the workflow. It will also output per read stats including average read lengths and qualities. Reads are additionally filtered for sequence length and quality characteristics.

### 2. Mapping and primer trimming (Artic)

Concatenated reads are mapped to the reference SARS-CoV-2 genome using [minimap2](https://github.com/lh3/minimap2). A primer scheme-specific BED file is used to identify the regions of
the mapped sequences that correspond to synthetic sequences (primers) - these regions are clipped to ensure that sequences are entirely of biological origin.

### 3. Variant calling and consensus generation (Artic)

The retained sequences are used to prepare a consensus sequence that is then polished using [Medaka](https://github.com/nanoporetech/medaka) and variant calling is performed to produce a VCF file of genetic differences relative to the reference genome.

### 4. Lineage/clade assignment

The consensus sequence is annotated for virus clade information using [NextClade](https://clades.nextstrain.org/), and strain assignment is performed using [Pangolin](https://github.com/cov-lineages/pangolin).



## Troubleshooting

+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).



## FAQ's

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-template/issues) page or start a discussion on the [community](https://community.nanoporetech.com/). 



## Related blog posts

+ [SARS-CoV-2 Midnight Analysis](https://labs.epi2me.io/sarscov2-midnight-analysis/)
+ [Midnight Scheme Update](https://labs.epi2me.io/ont-midnight-scheme-update/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.


## Problems and Solutions

1. Error during `fastcat` step
  ```
  Unable to find image 'ontresearch/wf-common:sha8b5843d549bb210558cbb676fe537a153ce771d6' locally
  docker: Error response from daemon: Get "https://registry-1.docker.io/v2/ontresearch/wf-common/manifests/sha256:de260cd5ba0131cb3c255ada8dfb3ca8a6f920f707eda4124601f92c95c71dd3": EOF.
  See 'docker run --help'.
  ```

  This is an internet problem and this particular container is a very large file - the docker pull request timed out while running. This file only needs to be pulled once. Either just try again, or in a terminal run
  ```
  docker pull ontresearch/wf-common:sha8b5843d549bb210558cbb676fe537a153ce771d6
  ```
  and resume the run once it has downloaded

