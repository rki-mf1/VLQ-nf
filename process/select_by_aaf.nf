process select_by_aaf {
  publishDir "${params.runinfo}/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "select_by_aaf.log"}
  publishDir "${params.databases}/build_reference/", mode: 'copy', pattern: "aaf_selection.csv", saveAs: {filename -> "reference.csv"}

  input:
  tuple val(lineage), path(merged_vcf), path(merged_freq)
  path selection_df

  output:
  path "aaf_selection.csv", emit: final_selection
  path ".command.log"

  script:
  """
  aaf_select.py --lineage $lineage --metadata $selection_df --vcf $merged_vcf --freq $merged_freq --min_aaf $params.min_aaf --max_per_lineage $params.max_per_lineage

  """

}
