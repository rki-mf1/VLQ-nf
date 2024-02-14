#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Nextflow -- Analysis pipeline to estimate sc2 lineage abundances from RT-qPCR sequences from wastewater samples
Authors: Eva, Martin
Based on a research project by Jasmijn A. Baaijens et al.
*/



/**************************
* META & HELP MESSAGES
**************************/

/**************************
* Help messages, user inputs & checks
**************************/
// help message
if (params.help) { exit 0, helpMSG() }

// Log infos based on user inputs
defaultMSG()

// error codes
if (params.profile) {
  exit 1, "--profile is WRONG use -profile" }

if (workflow.profile == 'standard') {
  "NO EXECUTION PROFILE SELECTED, using [-profile local,docker]" }

if (!params.reference){
  if (!params.gisaid) { exit 1, "No GISAID input data provided, please pass input using --gisaid parameter!" }
}

if (!params.query) {
  exit 1, "No query data specified, please provide folder containing fastq(.gz) query files!"
 }



/**************************
* INPUT channels
**************************/
if (!params.reference) {
  gisaid_meta_ch = channel.fromPath("${params.gisaid}/*metadata*", checkIfExists: true)
  gisaid_seq_ch = channel.fromPath("${params.gisaid}/*fa*", checkIfExists: true)
  gisaid_meta_ch.view()
  gisaid_seq_ch.view()
}

wildtype = channel.fromPath("${projectDir}/assets/NC_045512.2.fasta", checkIfExists:true)

if (params.single_end) {
  QUERY = channel.fromPath("${params.query}/*.fastq*", checkIfExists: true)
}
else {
  QUERY = channel.fromFilePairs("${params.query}/*{1,2}.fastq*", checkIfExists: true)
}
QUERY.view()



/**************************
* PROCESSES
**************************/
// include processes that should be used outside of a sub-workflow logic

include { process_gisaid } from './process/process_gisaid'
include { filter_by_metadata } from './process/filter_by_metadata'
include { filter_sequences; filter_sequences as filter_sequences_by_aaf } from './process/filter_sequences'
include { variant_call } from './process/variant_call'
include { merge_vcf } from './process/merge_vcf'
include { filter_by_aaf } from './process/filter_by_aaf'
include { build_index } from './process/build_index'
include { kallisto_prediction_single } from './process/kallisto_prediction_single'
include { kallisto_prediction_paired } from './process/kallisto_prediction_paired'



/**************************
* SUB-WORKFLOWS
**************************/
workflow process_input_data {

  take:
    gisaid_meta_ch
    gisaid_seq_ch

  main:

    gisaid_seq_ch.splitFasta(by: 2000, file:"sample.fasta").map{ file -> tuple("gisaid_${file.baseName}", file) }.set{ seq_ch }
    process_gisaid(gisaid_meta_ch)
    gisaid_in_ch = process_gisaid.out.gisaid_processed

  emit:
    gisaid_in_ch
    seq_ch

}


workflow build_reference_db {

  take:
    gisaid_in_ch
    seq_ch

  main:
    filter_by_metadata(gisaid_in_ch)
    selection_df = filter_by_metadata.out.selection_df
    selection_df.splitCsv(header: true, sep: '\t').map{ row -> row.record_id}.collectFile(newLine: true).set{ selected_ids }
    chunk_collector1 = seq_ch.combine(selected_ids)

    filter_sequences(chunk_collector1)
    // filter out empty files: covers case that a set of sequences to be filtered out happens to be comprised by one fasta chunk
    filtered_fasta = filter_sequences.out.filtered_fasta.filter{ it[1].size()>0 }
    chunk_collector2 = filtered_fasta.combine(wildtype).combine(selection_df)

    variant_call(chunk_collector2)
    variant_call_log = variant_call.out.log.collectFile(name: "variant_call.log", storeDir:"${params.runinfo}")
    chunk_lineages = variant_call.out.chunk_lineages
    chunk_lineages.collect{ it.splitCsv(header: false) }.flatten().unique().set{ lineage_ch }

    merge_vcf(lineage_ch)
    merge_log = merge_vcf.out.log.collectFile(name: "merge_vcf.log", storeDir:"${params.runinfo}")
    merged_lineage_vcf = merge_vcf.out.lineage
    merged_lineage_vcf.collectFile(newLine: true).set{ lineage_collector }

    filter_by_aaf(lineage_collector, selection_df)
    final_selection_df = filter_by_aaf.out.final_selection_df
    final_selection_df.splitCsv(header: true, sep: '\t').map{ row -> row.record_id}.collectFile(newLine: true).set{ final_ids }
    chunk_collector3 = filtered_fasta.combine( final_ids )

    filter_sequences_by_aaf(chunk_collector3)
    final_fasta_chunk = filter_sequences_by_aaf.out.filtered_fasta.filter{ it[1].size()>0 }
    final_fasta_chunk.map{ it -> it[1] }.set{ final_fasta }
    final_fasta.collectFile(newLine: true, name: "reference.fasta", storeDir: "${params.databases}/build_reference/").set{ reference_ch }

  emit:
    reference_ch
    final_selection_df

 }


workflow predict_abundances {

  take:
    reference_ch
    QUERY
    final_selection

  main:
    build_index(reference_ch)
    kallisto_idx = build_index.out.kallisto_idx

    if (params.single_end) {
      kallisto_in = QUERY.combine(kallisto_idx).combine(final_selection)
      kallisto_prediction_single(kallisto_in)
      kallisto_log = kallisto_prediction_single.out.log.collectFile(name: "predict_abundances.log", storeDir: "${params.runinfo}")
      prediction_ch = kallisto_prediction_single.out.prediction_ch
    }
    else {
      QUERY_data = QUERY.map{ it -> it[1] }
      QUERY_name = QUERY.map{ it -> it[0]}
      kallisto_in = QUERY_data.combine(kallisto_idx).combine(final_selection)
      kallisto_prediction_paired(kallisto_in, QUERY_name)
      kallisto_log = kallisto_prediction_paired.out.log.collectFile(name: "predict_abundances.log", storeDir: "${params.runinfo}")
      prediction_ch = kallisto_prediction_paired.out.prediction_ch
    }

  emit:
    prediction_ch
 }



/**************************
* MAIN WORKFLOW ENTRY POINT
**************************/
workflow {

  if (!params.reference) {
    process_input_data(gisaid_meta_ch, gisaid_seq_ch)
    gisaid_in_ch = process_input_data.out.gisaid_in_ch
    seq_ch  = process_input_data.out.seq_ch
    gisaid_in_ch.view()

    build_reference_db( gisaid_in_ch, seq_ch )
    reference_ch = build_reference_db.out.reference_ch
    final_selection_df = build_reference_db.out.final_selection_df
    final_selection_df.view()
    reference_ch.view()
  }
  else {
    reference_ch = channel.fromPath("${params.databases}/build_reference/*.fasta", checkIfExists: true)
    final_selection_df = channel.fromPath("${params.databases}/build_reference/*.csv", checkIfExists: true)
  }

  predict_abundances(reference_ch, QUERY, final_selection_df)
  prediction_ch = predict_abundances.out.prediction_ch
  prediction_ch.view()


  /**********************************************************
  * Coming soon: visualization of reference and output data
  ***********************************************************/

}



/*************
* OUTPUT
*************/
workflow.onComplete {

summary = """---Right here you will find a nice summary soon---"""

/*log_file = file("${params.runinfo}/summary.log")
logReader = log_file.newReader()
String line
while (line=logReader.readLine()) {
if (line.startsWith('---')) { summary = summary + line + "\n" }
}
logReader.close()
*/

log.info """
Execution status: ${ workflow.success ? 'OK' : 'failed' }
______________________________________
\u001B[36mExecution summary\033[0m
______________________________________
$summary
Summary report:                 ${params.runinfo}/
Lineage abundance predictions:  ${params.output}/NAME_OF_QUERY/kallisto_out/predictions.tsv

______________________________________
""".stripIndent()
}



/*************
* --help
*************/
def helpMSG() {
    log.info """
    ____________________________________________________________________________________________
    Workflow: VLQ-nf

    Usage example:
    nextflow run main.nf --gisaid PATH_TO_GISAID_DATA/ --query PATH_TO_QUERY_FASTQs/ --country Germany,England --startdate YYYY-MM-DD --enddate YYY-MM-DD


    Mandatory arguments:
    --gisaid                    Full path to folder storing metadata and sequence files from gisaid for reference building
                                Currently, metadata has to be tab-delimited csv or tsv file and sequence info has to be fasta or fasta.gz
                                Required columns: ['Virus name', 'Accession ID', 'Collection date', 'Location', 'Sequence length', 'Host', 'Pango lineage', 'N-Content']
    or
    --reference                 If true, the prebuilt reference in --databases is used to build an index and analyse the query
                                [default: false]
    --query                     Full path to folder storing query fastq files
                                (fastq files are allowed to be .gz compressed)

    Reference building:
    --continent                 List of comma-separated continents to consider when selecting reference samples based
                                on geography. Inform yourself how continents are captured by your input data resource.
                                Whitespaces in country names need to be replaced by underscore.
                                Example: "North_America"
                                [default: empty str]
    --country                   List of comma-separated countries to consider when selecting reference samples based
                                on geography. Inform yourself how countries are captured by your input data resource.
                                Whitespaces in country names need to be replaced by underscore.
                                Example: "Bosnia_Herzegovina"
                                [default: empty str]
    --startdate                 Earliest sampling date to consider when selecting reference samples based on the timepoint
                                of sample drawing. [default: empty str, format: YYYY-MM-DD]
    --enddate                   Latest sampling date to consider when selecting reference samples based on the timepoint
                                of sample drawing. [default: empty str, format: YYYY-MM-DD]
    --min_len                   Don't select sequences with less than the specified minimal number of non-ambiguous nucleotides.
                                [default: 29500]
    --k                         Specify how many samples to randomly select per lineage when filtering based on metadata.
                                [default: 1000]
    --seed                      Random seed for sampling steps during sequence selection (see --k). [default: 0]
    --min_aaf                   Minimal alternative allele frequency (AAF) to consider for representative lineage variation.
                                [default: 0.5]
    --max_per_lineage           Maximum number of lineages to select to represent a lineages genomic variation.
                                [default: 0]
    Kallisto:
    --single_end                If true, run kallisto parametrized for single-end reads, else for paired-end reads (assuming "1" and "2" as filename suffices)
                                [default: true]
    Currently, only the following kallisto parameters are adjustable via command line parameters:
    --fragment_length           Estimated average fragment length. [default: 200]
    --fragment_length_sd        Estimated standard deviation of fragment length. [default: 20]
    --kallisto_threads          Number of threads to use. [default: 20]
                                For more information on kallisto, e.g. regarding single and paired-end fastq queries, see
                                https://pachterlab.github.io/kallisto/manual
    --bootstrap                 Determine whether kallisto should be run in bootstrap mode and how many bootstraps should be used.
                                [default: 0]
    For more information on kallisto, see https://pachterlab.github.io/kallisto/manual

    Output:
    --min_ab                    Summarize output for all lineages whose estimated abundance is above this minimum threshold.
                                [default: 0]

    Nextflow options:
    --output                    Full path to output folder
    --runinfo                   Full path to folder storing log files
    --databases                 Full path to folder storing intermediate files and reference data

    Note: VLQ-nf requires full file paths

    ____________________________________________________________________________________________
    """.stripIndent()
}


def defaultMSG() {
    log.info """
    ______________________________________
    Workflow: VLQ-nf

    Profile:                $workflow.profile
    Current User:           $workflow.userName
    Nextflow-version:       $nextflow.version
    Starting time:          $nextflow.timestamp
    Workflow hash:          $workflow.commitId
        --workdir           $params.workdir
        --databases         $params.databases
        --output            $params.output
        --max_cores         $params.max_cores
        --memory            $params.memory
        --cachedir          $params.cachedir
    ______________________________________
    """.stripIndent()
}
