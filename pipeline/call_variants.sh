#!/bin/bash
#SBATCH -c 20
#SBATCH -t 0-2:00
#SBATCH -p short
#SBATCH --mem=10G
#SBATCH -o call_variants.out
#SBATCH -e call_variants.err

ref_dir=$1
cd ${ref_dir}

# for every lineage go to dir with sampeld fasta files
while read lineage; do \
    cd $lineage;
    # for every fasta file
    for fasta in *.fa; do \
        # align and sort, zip and index resulting vcf
        # install minimap2, k8 and paftools.js all via https://github.com/lh3/minimap2/tree/master/misc
        # TODO: add minimap2, k8 temporarily to $PATH via export PATH="$PATH:`pwd`:`pwd`/misc"
        minimap2 -c -x asm20 --end-bonus 100 -t 20 --cs $HOME/Dokumente/RKI/Internship/wildtype/NC_045512.2.fasta $fasta 2>${fasta%.fa}.paftools.log | sort -k6,6 -k8,8n > ${fasta%.fa}.paf && k8 $HOME/Dokumente/RKI/sc2_sewage/minimap2/misc/paftools.js call -s ${fasta%.fa} -L 100 -f /home/eva/Dokumente/RKI/Internship/wildtype/NC_045512.2.fasta ${fasta%.fa}.paf > ${fasta%.fa}.vcf 2>>${fasta%.fa}.paftools.log;
        bgzip -f ${fasta%.fa}.vcf;
        bcftools index -f ${fasta%.fa}.vcf.gz;
    done;
    cd ..;
    sample_count=$(ls ${lineage}/*.vcf.gz | wc -l);
    if [[ ${sample_count} -eq 1 ]]; then \
        cp ${lineage}/*.vcf.gz ${lineage}_merged.vcf.gz;
    else \
        bcftools merge -o ${lineage}_merged.vcf.gz -O z ${lineage}/*.vcf.gz;
    fi;
    vcftools --gzvcf ${lineage}_merged.vcf.gz --out ${lineage}_merged --site-pi;
    vcftools --gzvcf ${lineage}_merged.vcf.gz --out ${lineage}_merged --freq;
done < lineages.txt

cd ..
