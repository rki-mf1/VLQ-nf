#!/bin/bash

SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")
PAFTOOLS="/home/eva/Downloads/minimap2-master/misc/paftools.js"
G_META=$1
G_FASTA=$2
D_META=$3
D_LINES=$4
D_FASTA=$5
D_EPI=$6
REFERENCE=$7
QUERY=$8
OUTDIR=$9
LOG_FILE="${OUTDIR}/sc2_sewage.log"

mkdir -p $OUTDIR
touch $LOG_FILE
exec &> >(tee "$LOG_FILE")

echo "Step 1: Build a reference set"
echo ">>> pre-select based on nonN counts and country"
python ${SCRIPTPATH}/pipeline_baymlab/preprocess_references.py -gisaid ${G_META} ${G_FASTA} -desh ${D_META} ${D_LINES} ${D_FASTA} -epi ${D_EPI} --country Germany France Denmark England --startdate 2020-12-01 --enddate 2020-12-31 -o ${OUTDIR}/seqs_per_lineage --log ${LOG_FILE}
echo ">>> call variants: find log info in ${OUTDIR}/seqs_per_lineage/<FASTA>.paftools.log"
${SCRIPTPATH}/pipeline_baymlab/call_variants.sh ${OUTDIR}/seqs_per_lineage ${REFERENCE} ${PAFTOOLS}
echo ">>> select sequences per lineage such that each typical mutation with at least 50% frequency is captured at least once"
mkdir -p ${OUTDIR}/reference_set
python ${SCRIPTPATH}/pipeline_baymlab/select_samples.py -f ${G_FASTA} ${D_FASTA} --vcf ${OUTDIR}/seqs_per_lineage/*_merged.vcf.gz --freq ${OUTDIR}/seqs_per_lineage/*_merged.frq -o ${OUTDIR} --log ${LOG_FILE}


# Step 2: Preprocess seq data???


echo "Step 2: Predict variant abundances"
mkdir -p ${OUTDIR}/kallisto_out
kallisto index -i ${OUTDIR}/reference_set/sequences.kallisto_idx ${OUTDIR}/reference_set/sequences.fasta
for sample in SL;
do
  kallisto quant -i ${OUTDIR}/reference_set/sequences.kallisto_idx -o ${OUTDIR}/kallisto_out/${sample} --single -l 200 -s 20 -t 20 ${QUERY}/sra_${sample}.fastq.gz
  ## summarize lineage abundances
  python ${SCRIPTPATH}/pipeline_baymlab/output_abundances.py --metadata ${OUTDIR}/seqs_per_lineage/metadata.tsv ${OUTDIR}/kallisto_out/${sample}/abundance.tsv --log ${LOG_FILE}
done
