#!/bin/bash


echo "Step 1: Build a reference set" >
# pre-select based on nonN counts and country
#python pipeline/preprocess_references.py -m metadata.tsv -f sequences_germany_2021_10_01.fasta --country Germany
# call variants
#pipeline/call_variants.sh seqs_per_lineage
# select sequences per lineage such that each typical mutation with at least 50% frequency is captured at least once
#python pipeline/select_samples.py --vcf seqs_per_lineage/*_merged.vcf.gz --freq seqs_per_lineage/*_merged.frq -m metadata.tsv -f sequences_germany_2021_10_01.fasta -l seqs_per_lineage/lineages.txt -o reference_set


# Step 2: Preprocess seq data???


# Step 3: Predict variant abundances
#kallisto index -i reference_set/sequences.kallisto_idx reference_set/sequences.fasta
#for sample in GR NR SL;
#do
#  kallisto quant -i reference_set/sequences.kallisto_idx -o kallisto_out/${sample} --single -l 200 -s 20 -t 20 wastewater_reads/sra_${sample}.fastq.gz
  ## summarize lineage abundances
#  python pipeline/output_abundances.py --metadata reference_set/metadata.tsv kallisto_out/${sample}/abundance.tsv
#done
