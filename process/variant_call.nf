process variant_call {
  maxForks 4
  publishDir "${params.runinfo}/variant_call/", mode: 'copy', pattern: "*.paftools.log"
  publishDir "${params.databases}/build_reference/vcf/${lineage}/", mode: 'copy', pattern: "${record_id}.vcf.gz*"

  input:
  tuple val(fasta_id), val(record_id), val(lineage), val(seq)
  path wildtype

  output:
  tuple val(lineage), path("${record_id}.vcf.gz*"), emit: lineage
  path "*.paftools.log"

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

  echo ">${fasta_id}" > ${record_id}.fasta
  echo $seq >> ${record_id}.fasta

  PAF=\$(which paftools.js)

  echo \${PAF}
  minimap2 -c -x asm20 --end-bonus 100 -t 20 --cs $wildtype ${record_id}.fasta 2>${record_id}.paftools.log | sort -k6,6 -k8,8n > ${record_id}.paf && \${PAF} call -s ${record_id} -L 100 -f $wildtype ${record_id}.paf > ${record_id}.vcf 2>>${record_id}.paftools.log;
  bgzip -f ${record_id}.vcf;
  bcftools index -f ${record_id}.vcf.gz;

  """

}
