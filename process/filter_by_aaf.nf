process filter_by_aaf {
  //publishDir "${params.runinfo}/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "filter_by_aaf.log"}
  publishDir "${params.databases}/build_reference/", mode: 'copy', pattern: "aaf_selection.csv", saveAs: {filename -> "reference.csv"}

  input:
  path merged_vcf_lineages
  path selection_df

  output:
  path "aaf_selection.csv", emit: final_selection_df
  path ".command.log", emit: log

  script:
  """
  #!/bin/bash
  echo "_______________________________________________________________"
  echo "--- Filter lineage-representative samples based on alternate allele frequencies"

  filter_by_aaf.py --metadata $selection_df --vcf ${params.databases}/build_reference/merged_vcf/*_merged.vcf.gz --freq ${params.databases}/build_reference/merged_vcf/*_merged.frq --min_aaf $params.min_aaf --max_per_lineage $params.max_per_lineage
  """

}
