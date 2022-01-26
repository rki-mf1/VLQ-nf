process merge_vcf {
  publishDir "${params.runinfo}/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "merge_vcf.log"}
  publishDir "${params.databases}/build_reference/vcf/", mode: 'copy', pattern: "${lineage}_merged.*"

  input:
  val lineage

  output:
  tuple val(lineage), path("${lineage}_merged.vcf.gz"), path("${lineage}_merged.frq"), emit: merged_vcf
  path ".command.log"

  script:
  """
  #!/bin/bash
  #SBATCH -c 20
  #SBATCH -t 0-2:00
  #SBATCH -p short
  #SBATCH --mem=10G
  #SBATCH -o call_variants.out
  #SBATCH -e call_variants.err
  ##############################################
  # Source: https://github.com/baymlab/wastewater_analysis
  #############################################

  sample_count=\$(ls ${projectDir}/${params.databases}/build_reference/vcf/${lineage}/*.vcf.gz | wc -l);
  if [[ \$sample_count -eq 1 ]]; then \
      cp ${projectDir}/${params.databases}/build_reference/vcf/${lineage}/*.vcf.gz ${lineage}_merged.vcf.gz;
  else \
      bcftools merge -o ${lineage}_merged.vcf.gz -O z ${projectDir}/${params.databases}/build_reference/vcf/${lineage}/*.vcf.gz;
  fi;
  # TODO: site-pi information needed?
  vcftools --gzvcf ${lineage}_merged.vcf.gz --out ${lineage}_merged --site-pi;
  vcftools --gzvcf ${lineage}_merged.vcf.gz --out ${lineage}_merged --freq;

  rm -rf ${projectDir}/${params.databases}/build_reference/vcf/${lineage}/*

  """

}
