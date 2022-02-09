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

if (!params.gisaid) { exit 1, "No GISAID input data provided, please pass input using --gisaid parameter!" }
if (params.desh && !params.gisaid_desh_map) { exit 1, "If using both GISAID and DESH data, a mapping file for EPI_ISL ids is required!"}
else if (!params.gisaid && params.desh) { exit 1, "EPI_ISL mapping file provided without DESH data! (DESH duplicates would be removed from GISAID dataset leading to data loss)" }

if (!params.query) {
  exit 1, "No query data specified, please provide folder containing fastq(.gz) query files!"
 }



/**************************
* INPUT channels
**************************/

gisaid_meta_ch = channel.fromPath("${params.gisaid}/*metadata*", checkIfExists: true)
gisaid_seq_ch = channel.fromPath("${params.gisaid}/*fasta*", checkIfExists: true)
gisaid_meta_ch.view()
gisaid_seq_ch.view()

wildtype = channel.fromPath("${projectDir}/assets/NC_045512.2.fasta", checkIfExists:true)

QUERY = channel.fromPath("${params.query}/*.fastq*", checkIfExists: true)



/**************************
* PROCESSES
**************************/
// include processes that should be used outside of a sub-workflow logic

//include { download_desh } from './process/download_desh'
include { process_desh } from './process/process_desh'
include { process_gisaid } from './process/process_gisaid'
include { filter_by_metadata } from './process/filter_by_metadata'
include { gzip_fasta; gzip_fasta as gzip_fasta2 } from './process/gzip_fasta'
include { filter_sequences; filter_sequences as filter_sequences_by_aaf } from './process/filter_sequences'
include { variant_call } from './process/variant_call'
include { merge_vcf } from './process/merge_vcf'
include { filter_by_aaf } from './process/filter_by_aaf'
include { build_index } from './process/build_index'
include { kallisto_prediction } from './process/kallisto_prediction'



/**************************
* DATABASES
*************************/

workflow get_desh {

  main:
    if (params.desh) {
      desh_seq_ch = channel.fromPath("${params.desh_data}/*fasta*", checkIfExists: true)
      desh_meta_ch = channel.fromPath("${params.desh_data}/SARS-CoV-2-Sequenzdaten_Deutschland.csv.xz", checkIfExists: true)
      desh_lineage_ch = channel.fromPath("${params.desh_data}/SARS-CoV-2-Entwicklungslinien_Deutschland.csv.xz", checkIfExists: true)

      desh_meta_ch.concat(desh_lineage_ch).collect().set{ desh_csv_ch }
      gisaid_desh_map = channel.fromPath("${params.gisaid_desh_map}", checkIfExists: true)
    }

    else {
      gisaid_desh_map = channel.empty()
      desh_csv_ch = channel.empty()
      desh_seq_ch = channel.empty()
     }

  emit:
    desh_csv_ch
    desh_seq_ch
    gisaid_desh_map
}



/**************************
* SUB-WORKFLOWS
**************************/

workflow process_input_data {

  take:
    desh_csv_ch
    desh_seq_ch
    gisaid_desh_map

  main:
    //gzip_fasta(gisaid_seq_ch)
    //gisaid_seq_gz = gzip_fasta.out.gz_fasta
    gisaid_seq_gz = gisaid_seq_ch

    if (params.desh) {
      //gzip_fasta2(desh_seq_ch)
      //desh_seq_gz = gzip_fasta2.out.gz_fasta
      desh_seq_gz = desh_seq_ch

      gisaid_seq_gz.splitFasta(by: 2000, file:"sample.fasta").map{ file -> tuple("gisaid_${file.baseName}", file) }.set{ gisaid_fasta_ch }
      desh_seq_gz.splitFasta(by: 2000, file:"sample.fasta").map{ file -> tuple("desh_${file.baseName}", file) }.set{ desh_fasta_ch }
      gisaid_fasta_ch.concat(desh_fasta_ch).set{ seq_ch }

      process_desh(desh_csv_ch, desh_seq_gz)
      desh_in_ch = process_desh.out.desh_processed
      desh_log = process_desh.out.log

    }
    else {
      gisaid_seq_gz.splitFasta(by: 2000, file:"sample.fasta").map{ file -> tuple("gisaid_${file.baseName}", file) }.set{ seq_ch }
      desh_in_ch = channel.empty()
      desh_log = channel.empty()
     }

    process_gisaid(gisaid_meta_ch, gisaid_desh_map.ifEmpty(file("${projectDir}/assets/DUMMY")))
    gisaid_in_ch = process_gisaid.out.gisaid_processed
    gisaid_log = process_gisaid.out.log

    gisaid_log.concat(desh_log).collectFile(name: "process_input_data.log", storeDir: "${params.runinfo}").set{ process_input_data_log }

  emit:
    gisaid_in_ch
    desh_in_ch
    seq_ch
    process_input_data_log

}


workflow build_reference_db {

  take:
    desh_in_ch
    gisaid_in_ch
    seq_ch

  main:
    filter_by_metadata(desh_in_ch, gisaid_in_ch)
    selection_df = filter_by_metadata.out.selection_df
    selection_log = filter_by_metadata.out.log
    selection_df.splitCsv(header: true, sep: '\t').map{ row -> row.fasta_id}.collectFile(newLine: true).set{ selected_ids }
    chunk_collector1 = seq_ch.combine(selected_ids)
    selection_df.splitCsv(header: true, sep: '\t').map{ row -> row.fasta_id}.flatten().unique().count().view()
    selection_df.splitCsv(header: true, sep: '\t').map{ row -> row.lineage}.flatten().unique().count().view()

    filter_sequences(chunk_collector1)
    filter_log = filter_sequences.out.log.collectFile()
    // filter out empty files: covers case that a set of sequences to be filtered out happens to  be copmrised by one fasta chunk
    filtered_fasta = filter_sequences.out.filtered_fasta.filter{ it[1].size()>0 }
    filtered_fasta_counter = filter_sequences.out.sample_counter.filter{ it.size()>0 }
    filtered_fasta_count = filtered_fasta_counter.collect{ it.splitCsv(header: false) }.flatten().map{ it -> it.replaceAll( />/, "" )}
    filtered_fasta_count.filter{ item -> item=='EPI_ISL_1195379'}.view()
    filtered_fasta_count.count().view()
    filtered_lineages = filtered_fasta_count.join(selection_df.splitCsv(header: true, sep: '\t').map{ row ->  tuple(row.fasta_id, row.lineage) }).map{ it -> it[1] }.unique()
    filtered_lineages.count().view()
    chunk_collector2 = filtered_fasta.combine(wildtype).combine(selection_df)
    chunk_collector2.filter{ item -> item=='EPI_ISL_1195379'}.view()

    variant_call(chunk_collector2)
    variant_call_log = variant_call.out.log.collectFile()
    chunk_lineages = variant_call.out.chunk_lineages
    chunk_samples = variant_call.out.chunk_samples
    chunk_lineages.collect{ it.splitCsv(header: false) }.flatten().unique().set{ lineage_ch }
    chunk_samples.collect{ it.splitCsv(header: false) }.flatten().unique().set{ sample_ch }
    sample_ch.filter{ item -> item=='EPI_ISL_1195379'}.view()
    lineage_ch.count().view()
    sample_ch.count().view()

    merge_vcf(lineage_ch)
    merge_log = merge_vcf.out.log.collectFile()
    merged_lineage_vcf = merge_vcf.out.lineage
    merged_lineage_vcf.collectFile(newLine: true).set{ lineage_collector }
    vcf_dir = file("${params.databases}/vcf")
    tmp_file = vcf_dir.deleteDir()
    println tmp_file ? "Delete $tmp_file" : "Cannot delete: $tmp_file"


    filter_by_aaf(lineage_collector, selection_df)
    final_selection_log = filter_by_aaf.out.log
    final_selection_df = filter_by_aaf.out.final_selection_df
    final_selection_df.splitCsv(header: true, sep: '\t').map{ row -> row.fasta_id}.collectFile(newLine: true).set{ final_ids }
    final_ids.collect{ it.splitCsv(header:false) }.flatten().count().view()
    chunk_collector3 = filtered_fasta.combine( final_ids )

    filter_sequences_by_aaf(chunk_collector3)
    final_filter_log = filter_sequences_by_aaf.out.log.collectFile()
    final_fasta_headers = filter_sequences_by_aaf.out.sample_counter.filter{ it.size() > 0 }
    final_fasta_ids = final_fasta_headers.collect{ it.splitCsv(header: false) }.flatten().map{ it -> it.replaceAll( />/, "" )}
    final_lineages = final_fasta_ids.join(final_selection_df.splitCsv(header: true, sep: '\t').map{ row ->  tuple(row.fasta_id, row.lineage) }).map{ it -> it[1] }.unique()
    final_fasta_ids.count().view()
    final_lineages.count().view()
    final_fasta_chunk = filter_sequences_by_aaf.out.filtered_fasta.filter{ it[1].size()>0 }
    final_fasta_chunk.map{ it -> it[1] }.set{ final_fasta }
    final_fasta.view()
    final_fasta.collectFile(newLine: true, name: "reference.fasta", storeDir: "${params.databases}/build_reference/").set{ reference_ch }

    final_filter_log.concat(final_selection_log, merge_log, variant_call_log, filter_log,  selection_log).collectFile(name: "build_reference_db.log", storeDir:"${params.runinfo}").set{ build_reference_db_log }

  emit:
    reference_ch
    final_selection_df
    build_reference_db_log

 }


workflow predict_abundances {

   take:
     reference_ch
     QUERY
     final_selection

   main:
     build_index(reference_ch)
     build_index_log = build_index.out.log
     kallisto_idx = build_index.out.kallisto_idx
     kallisto_in = QUERY.combine(kallisto_idx).combine(final_selection)

     kallisto_prediction(kallisto_in)
     kallisto_log = kallisto_prediction.out.log.collectFile()
     prediction_ch = kallisto_prediction.out.prediction_ch

     kallisto_log.concat(build_index_log).collectFile(name: "predict_abundances.log", storeDir: "${params.runinfo}").set{ predict_abundances_log }

   emit:
     prediction_ch
     predict_abundances_log
 }



/**************************
* MAIN WORKFLOW ENTRY POINT
**************************/
workflow {

  get_desh()
  desh_csv_ch = get_desh.out.desh_csv_ch
  desh_seq_ch = get_desh.out.desh_seq_ch
  gisaid_desh_map = get_desh.out.gisaid_desh_map
  desh_csv_ch.view()
  desh_seq_ch.view()

  process_input_data(desh_csv_ch, desh_seq_ch, gisaid_desh_map)
  desh_in_ch = process_input_data.out.desh_in_ch
  gisaid_in_ch = process_input_data.out.gisaid_in_ch
  seq_ch  = process_input_data.out.seq_ch
  process_input_data_log = process_input_data.out.process_input_data_log
  desh_in_ch.view()
  gisaid_in_ch.view()
  process_input_data_log.view()

  build_reference_db( desh_in_ch.ifEmpty(file("${projectDir}/assets/DUMMY")), gisaid_in_ch, seq_ch )
  reference_ch = build_reference_db.out.reference_ch
  final_selection_df = build_reference_db.out.final_selection_df
  build_reference_db_log = build_reference_db.out.build_reference_db_log
  final_selection_df.view()
  reference_ch.view()
  //build_reference_db_log.view()

  /*
  predict_abundances(reference_ch, QUERY, final_selection_df)
  prediction_ch = predict_abundances.out.prediction_ch
  predict_abundances_log = predict_abundances.out.predict_abundances_log
  prediction_ch.view()
  predict_abundances_log.view()

  process_input_data_log.concat(build_reference_db_log, predict_abundances_log).collectFile(name: "summary.log", storeDir: "${params.runinfo}").set{ summary_log }
  */
  build_reference_db_log.concat(process_input_data_log).collectFile(name: "summary.log", storeDir: "${params.runinfo}").set{ summary_log }

  summary_log.view()

  /**********************************************************
  * Coming soon: visualization of reference and output data
  ***********************************************************/

}



/*************
* OUTPUT
*************/


workflow.onComplete {

summary = """"""
log_file = file("${params.runinfo}/summary.log")
logReader = log_file.newReader()
String line
while (line=logReader.readLine()) {
if (line.startsWith('---')) { summary = summary + line + "\n" }
}
logReader.close()

log.info """
Execution status: ${ workflow.success ? 'OK' : 'failed' }
______________________________________
\u001B[36mExecution summary\033[0m
______________________________________
$summary
Summary report:                 ${params.runinfo}summary.log
Lineage abundance predictions:  ${params.output}QUERY/kallisto_out/predictions.tsv

______________________________________
""".stripIndent()
}


/*************
* --help
*************/
def helpMSG() {
    log.info """
    ____________________________________________________________________________________________
    Workflow: LinACov (working title)

    Usage example:
    nextflow run main.nf --gisaid PATH/TO/GISAID_DATA/ --query PATH/TO/QUERY_FASTQs/ --country Germany England --startdate YYYY-MM-DD --enddate YYY-MM-DD --vocs VOC.1 VOC.2 VOC.2.1


    Mandatory arguments:
    --gisaid                    Path to folder storing metadata and sequence files from gisaid
                                (file names: "*metadata*" (tsv) and "*fasta*" (fasta)).
    --query                     Path to folder storing query fastq files
                                (fastq files are allowed to be .gz compressed)

    Optional arguments:
    Input:
    --desh                      Boolean parameter, if true, program includes DESH data into the reference building
                                [default: false]
    --desh_data                 Path to folder storing metadata, lineage data and sequence data from gisaid
                                metadata filenames:     SARS-CoV-2-Sequenzdaten_Deutschland.csv.xz,
                                                        SARS-CoV-2-Entwicklungslinien_Deutschland.csv.xz,
                                sequence data filename: "*fasta*" (fasta))
                                [default: empty str]
    --gisaid_desh_map           Path to .csv file mapping accession ids of overlapping samples between desh and gisaid.
                                ATTENTION: If desh input is provided, it is recommended to also input such a mapping
                                file to avoid duplicates in the reference set! [default: empty str]
    Reference building:
    --continent                 List of comma-separated continents to consider when selecting reference samples based
                                on geography. [default: empty str]
    --country                   List of comma-separated countries to consider when selecting reference samples based
                                on geography. [default: empty str]
    --startdate                 Earliest sampling date to consider when selecting reference samples based on the timepoint
                                of sample drawing. [default: empty str, format: YYYY-MM-DD]
    --enddate                   Latest sampling date to consider when selecting reference samples based on the timepoint
                                of sample drawing. [default: empty str, format: YYYY-MM-DD]
    --min_len                   Don't select sequences with less than the specified minimal number of non-ambiguous nucleotides.
                                [default: 29500]
    --k                         Specify how many samples to randomly select per lineage when filtering based on metadata.
                                [default: 1000]
    --seed                      Random seed for sequence selection (see --k). [default: 0]
    --min_aaf                   Minimal alternative allele frequency (AAF) to consider for representative lineage variation.
                                [default: 0.5]
    --max_per_lineage           Maximum number of lineages to select to represent a lineages genomic variation.
                                [default: 100]
    Kallisto:
    --fragment_length           Estimated average fragment length. [default: 200]
    --fragment_length_sd        Estimated standard deviation of fragment length. [default: 20]
    --kallisto_threads          Number of threads to use. [default: 20]
                                For more information on kallisto, e.g. regarding single and paired-end fastq queries, see
                                https://pachterlab.github.io/kallisto/manual
    Output:
    --min_ab                    Summarize output for all lineages whose estimated abundance is above this minimum threshold.
                                [default: 0]
    --vocs                      Comma-separated list of variants of interest. If empty, output is generated for all lineages
                                comprised in the reference. [default: empty list]


    Nextflow options:

    LSF computing:
    For execution of the workflow on a HPC with LSF adjust the following parameters:

    Profiles:

    ____________________________________________________________________________________________
    """.stripIndent()
}

def defaultMSG() {
    log.info """
    ______________________________________
    Workflow: LinACov (working title)

    Profile:                $workflow.profile
    Current User:           $workflow.userName
    Nextflow-version:       $nextflow.version
    Starting time:          $nextflow.timestamp
    Workflow hash:          $workflow.commitId
        --workdir           $params.workdir
        --databases         $params.databases
        --output            $params.output
        --cores             $params.cores
        --max_cores         $params.max_cores
        --memory            $params.memory
        --cachedir          $params.cachedir
    ______________________________________
    """.stripIndent()
}
