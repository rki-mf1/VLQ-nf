process merge_vcf {
  //publishDir "${params.runinfo}/merge_vcf/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "${lineage}_merge_vcf.log"}
  publishDir "${params.databases}/build_reference/merged_vcf/", mode: 'copy', pattern: "${lineage}_merged.*"

  input:
  val lineage

  output:
  path "${lineage}_merged.vcf.gz"
  path "${lineage}_merged.frq"
  val lineage, emit: lineage
  path ".command.log", emit: log

  script:
  """
  #!/bin/bash
  ##############################################
  # Source: https://github.com/baymlab/wastewater_analysis
  #############################################
  echo "# Merge all vcf files for lineage ${lineage}"

  sample_count=\$(ls ${params.databases}/build_reference/vcf/${lineage}_*_merged.vcf.gz | wc -l)

  if [[ \$sample_count -eq 1 ]]; then
    cp ${params.databases}/build_reference/vcf/${lineage}_*_merged.vcf.gz ${lineage}_merged.vcf.gz
  elif [[ \$sample_count -gt 1 ]]; then
    for file in ${params.databases}/build_reference/vcf/${lineage}_*_merged.vcf.gz; do
      bcftools index -f \$file
    done
    bcftools merge -o ${lineage}_merged.vcf.gz -O z ${params.databases}/build_reference/vcf/${lineage}_*_merged.vcf.gz
  fi
  # TODO: site-pi information needed?
  vcftools --gzvcf ${lineage}_merged.vcf.gz --out ${lineage}_merged --site-pi
  vcftools --gzvcf ${lineage}_merged.vcf.gz --out ${lineage}_merged --freq

  """

}
