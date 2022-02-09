process variant_call {
  maxForks 4
  publishDir "${params.runinfo}/build_reference/", mode: 'copy', pattern: "${chunk_id}_samples.txt"
  publishDir "${params.databases}/build_reference/vcf", mode: 'copy', pattern: "*_${chunk_id}_merged.vcf.gz"

  input:
  tuple val(chunk_id), path(multifasta), path(wildtype), path(selection_df)

  output:
  path "${chunk_id}_lineages.txt", emit: chunk_lineages
  path "${chunk_id}_samples.txt", emit: chunk_samples
  path "*_${chunk_id}_merged.vcf.gz"
  path "${chunk_id}_paftools.log", emit: log

  script:
  """
  #!/bin/env bash
  #SBATCH -c 20
  #SBATCH -t 0-2:00
  #SBATCH -p short
  #SBATCH --mem=10G
  #SBATCH -o call_variants.out
  #SBATCH -e call_variants.err
  ##############################################
  # Source: https://github.com/baymlab/wastewater_analysis
  #############################################
  echo "Call variants for fasta chunk ${chunk_id}"

  PAF=\$(which paftools.js)
  mkdir -p fasta/

  echo Split multifasta and collect lineage annotation
  split_fasta_collect_lineage.py -meta $selection_df -multifasta $multifasta

  for file in fasta/*.fasta; do
    basename=\${file##*/}
    fasta_id=\${basename%.fasta}
    lineage=\${fasta_id%%_*}

    minimap2 -c -x asm20 --end-bonus 100 -t 20 --cs $wildtype \$file 2>\${fasta_id}.paftools.log | sort -k6,6 -k8,8n > \${fasta_id}.paf && \${PAF} call -s \$fasta_id -L 100 -f $wildtype \${fasta_id}.paf > \${fasta_id}.vcf 2>>\${fasta_id}.paftools.log;
    bgzip -f \${fasta_id}.vcf;
    bcftools index -f \${fasta_id}.vcf.gz;
    mv  \${fasta_id}.vcf*  \${fasta_id}.paf* fasta/\$lineage
  done

  touch ${chunk_id}_lineages.txt
  for dir in fasta/*; do
    if [ -d \${dir} ]; then
      lineage=\${dir##*/}
      echo merge vcf files for lineage \${lineage}
      echo \$lineage >> ${chunk_id}_lineages.txt

      sample_count=\$(ls \${dir}/*.vcf.gz | wc -l);
      if [[ \$sample_count -eq 1 ]]; then
        cp \${dir}/*.vcf.gz \${lineage}_${chunk_id}_merged.vcf.gz
      else
        bcftools merge -o \${lineage}_${chunk_id}_merged.vcf.gz -O z \$dir/*.vcf.gz
      fi
    fi
  done

  touch ${chunk_id}_paftools.log
  touch ${chunk_id}_samples.txt
  cat .command.log >> ${chunk_id}_paftools.log
  for file in fasta/*/*paftools.log; do
    cat \$file >> ${chunk_id}_paftools.log
    filename=\${file##*/}
    tmp_filename=\${filename#*_}
    echo \${tmp_filename%.paftools.log}>> ${chunk_id}_samples.txt
  done

  """
}
