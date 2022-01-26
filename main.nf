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

/* Allows the workflow to be started only if at least gisaid data input is provided.
   If both gisaid and desh data are provided, the epi_isl mapping file is required */
if (!params.gisaid) { exit 1, "No GISAID input data provided, please pass input using --gisaid parameter!" }
if (params.desh && !params.gisaid_desh_map) { exit 1, "If using both GISAID and DESH data, a mapping file for epi_isl ids is required!"}
else if (!params.gisaid && params.desh) { exit 1, "Epi_isl mapping file provided without DESH data! (DESH duplicates would be removed from GISAID dataset leading to data loss)" }

if (!params.query) {
  exit 1, "No query data specified, please provide folder containing fastq query files!"
 }



/**************************
* INPUT channels
**************************/
gisaid_meta_ch = channel.fromPath("${params.gisaid}/metadata_tsv_*.tar.xz", checkIfExists: true)
gisaid_meta_ch.view()
gisaid_seq_ch = channel.fromPath("${params.gisaid}/sequences_fasta_*.tar.xz", checkIfExists: true)


wildtype = file("${projectDir}/assets/NC_045512.2.fasta", checkIfExists:true)

/*  all fastq sample files for a country have to be merged into one fastq file
    all country fastq files are stored in $params.query/
    fastq files can be compressed
*/
QUERY = channel.fromPath("${params.query}/*.fastq*", checkIfExists: true)




/**************************
* PROCESSES
**************************/
// include processes that should be used outside of a sub-workflow logic

include { download_desh } from './process/download_desh'
include { process_desh } from './process/process_desh'
include { process_gisaid } from './process/process_gisaid'
include { filter_by_metadata } from './process/filter_by_metadata'
include { unxz} from './process/unxz'
include { untar} from './process/untar'
include { variant_call } from './process/variant_call'
include { merge_vcf } from './process/merge_vcf'
include { select_by_aaf } from './process/select_by_aaf'
include { build_index } from './process/build_index'
include { kallisto_prediction } from './process/kallisto_prediction'



/**************************
* DATABASES
*************************/

workflow get_desh {
  main:
    // Note: There might be sth not as I want it: even if desh data is already present, download_desh was tried to run
    if (params.desh == true) {
      desh_meta_preload = file("${params.databases}/DESH/SARS-CoV-2-Sequenzdaten_Deutschland.csv.xz")
      desh_lineage_preload = file("${params.databases}/DESH/SARS-CoV-2-Entwicklungslinien_Deutschland.csv.xz")
      desh_seq_preload = file("${params.databases}/DESH/SARS-CoV-2-Sequenzdaten_Deutschland.fasta.xz")

      if (desh_meta_preload.exists() && desh_lineage_preload.exists() && desh_seq_preload.exists()) {
        desh_meta_ch = desh_meta_preload
        desh_lineage_ch = desh_lineage_preload
        desh_seq_ch = desh_seq_preload
        }
      else {
        download_desh()
        desh_meta_ch = download_desh.out.meta
        desh_lineage_ch = download_desh.out.lines
        desh_seq_ch = download_desh.out.seq
      }

      desh_meta_ch.concat(desh_lineage_ch).collect().set{ desh_csv_ch }
      desh_csv_ch.view()
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
    untar(gisaid_seq_ch)
    gisaid_seq_gz = untar.out.untar_file

    if (params.desh) {
      unxz(desh_seq_ch)
      desh_seq_gz = unxz.out.unxz_file

      gisaid_seq_gz.splitFasta(record: [id: true, seqString: true]).set{ gisaid_fasta_ch }
      desh_seq_gz.splitFasta(record: [id: true, seqString: true]).set{ desh_fasta_ch }
      gisaid_fasta_ch.concat(desh_fasta_ch).set{ seq_ch }

      process_desh(desh_csv_ch, desh_seq_gz)
      desh_in_ch = process_desh.out.desh_processed
    }
    else {
      // MARK
      gisaid_seq_gz.splitFasta(record: [id: true, seqString: true]).set{ seq_ch }
      // Note: limited to 100
      //gisaid_seq_ch.splitFasta(limit: 100, record: [id: true, seqString: true]).set{ seq_ch }
      desh_in_ch = channel.empty()
     }

    process_gisaid(gisaid_meta_ch, gisaid_desh_map.ifEmpty(file("${projectDir}/assets/DUMMY")))
    gisaid_in_ch = process_gisaid.out.gisaid_processed

  emit:
    gisaid_in_ch
    desh_in_ch
    seq_ch

}


workflow build_reference_db {

  take:
    desh_in_ch
    gisaid_in_ch
    seq_ch

  main:
    filter_by_metadata(desh_in_ch, gisaid_in_ch)
    selection_df = filter_by_metadata.out.selection_df // output: path to df
    selection_df.view()
    // Get sequence data for selected samples:
    selection_df.splitCsv(header: true, sep: '\t').map{ row -> tuple(row.fasta_id, row.record_id, row.lineage) }.set{ lineage_sample_map }
    lineage_sample_map.join(seq_ch.map{ fasta -> tuple(fasta.id, fasta.seqString) }).set{ lineage_fasta_map } // lineage_fasta_map = [fasta_id, record_id, lineage, fasta_seqString]

    variant_call(lineage_fasta_map, wildtype)
    lineage_collector = variant_call.out.lineage
    // Get unique lineages and merge vcf files for corresponding samples
    lineage_collector.map{ it[0] }.unique().set{ lineage_set }
    merge_vcf(lineage_set)
    vcf_ch = merge_vcf.out.merged_vcf

    select_by_aaf(vcf_ch, selection_df)
    final_selection = select_by_aaf.out.final_selection
    final_selection.view()
    final_selection.splitCsv(header: true, sep:'\t').map{ row -> tuple(row.fasta_id) }.set{ selected_ids }
    selected_ids.join(seq_ch.map{ fasta -> tuple(fasta.id, fasta.seqString) }).collectFile(name: "reference.fasta", storeDir: "${params.databases}/build_reference/"){ '>'+it[0]+'\n'+it[1]+'\n' }.set{ reference_ch }
    reference_ch.view()

  emit:
    reference_ch
    final_selection
 }


workflow predict_abundances {

  take:
    reference_ch
    QUERY
    final_selection

  main:
    build_index(reference_ch)
    kallisto_idx = build_index.out.kallisto_idx

    kallisto_prediction(QUERY, kallisto_idx, final_selection)
    prediction_ch = kallisto_prediction.out.prediction_ch

  emit:
    prediction_ch
}



/**************************
* MAIN WORKFLOW ENTRY POINT
**************************/
workflow {

  get_desh()
  desh_csv_ch = get_desh.out.desh_csv_ch
  desh_seq_ch = get_desh.out.desh_seq_ch
  gisaid_desh_map = get_desh.out.gisaid_desh_map

  process_input_data(desh_csv_ch, desh_seq_ch, gisaid_desh_map)
  desh_in_ch = process_input_data.out.desh_in_ch
  gisaid_in_ch = process_input_data.out.gisaid_in_ch
  seq_ch  = process_input_data.out.seq_ch
  desh_in_ch.view()
  gisaid_in_ch.view()

  build_reference_db( desh_in_ch.ifEmpty(file("${projectDir}/assets/DUMMY")), gisaid_in_ch.ifEmpty(file("${projectDir}/assets/DUMMY")), seq_ch )
  reference_ch = build_reference_db.out.reference_ch
  final_selection = build_reference_db.out.final_selection

  predict_abundances(reference_ch, QUERY, final_selection)
  prediction_ch = predict_abundances.out.prediction_ch
  prediction_ch.view()

}



/*************
* OUTPUT
*************/






/*************
* --help
*************/
def helpMSG() {
    log.info """
    ____________________________________________________________________________________________
    Workflow: LinACov (working title)

    Usage example:
    nextflow run main.nf --gisaid PATH/TO/GISAID_DATA/ --query PATH/TO/QUERY_FASTQs/ --country Germany England --startdate YYYY-MM-DD --enddate YYY-MM-DD --vocs B.1.177 P.1



    Mandatory arguments:
    --gisaid                    Path to folder storing metadata and sequence files from gisaid
                                (file names: "metadata_tsv_YYY_MM_DD.tar.xz" and "sequences_fasta_YYYY_MM_DD.tar.xz")
    --query                     Path to folder storing query fastq files
                                (fastq files are allowed to be .gz compressed)


    Optional arguments:
    Input:
    --desh                      Path to folder storing metadata and sequence files from gisaid
                                (file names: "metadata_tsv_YYY_MM_DD.tar.xz" and "sequences_fasta_YYYY_MM_DD.tar.xz").
                                [default: empty str]
    --gisaid_desh_map           Path to .csv file mapping accession ids of overlapping samples between desh and gisaid.
                                ATTENTION: If desh input is provided, it is recommended to also input such a mapping
                                file to avoid duplicates in the reference set! [default: empty str]
    Reference building:
    --country                   List of blank separated countries to consider when selecting reference samples based
                                on geography. [default: empty str]
    --startdate                 Earliest sampling date to consider when selecting reference samples based on the timepoint
                                of sample drawing. [default: empty str, format: YYY-MM-DD]
    --enddate                   Latest sampling date to consider when selecting reference samples based on the timepoint
                                of sample drawing. [default: empty str, format: YYY-MM-DD]
    --min_len                   Don't select sequences with less than the specified minimal number of non-ambiguous nucleotides.
                                [default: 29500]
    --k                         Specify how many samples to randomly select per lineage when filtering based on metadata.
                                [default: 1000]
    --seed                      Random seed for sequence selection (see --k). [default: 0]
    --min_aaf                   Minimal alternative allele frequency (AAF) to consider for representative lineage variation.
                                [default: 0.5]
    --max_per_lineage           Maximum number of lineages to select to represent a lineages genomic variation.
                                [default: 100]
    Output:
    --min_ab                    Summarize output for all lineages whose estimated abundance is above this minimum threshold.
                                [default: 0]
    --vocs                      Whitespace-separated list of variants of interest. If empty, output is generated for all lineages
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
