#!/bin/bash

SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")
################################################
# Adjust local file path for paftools.js
################################################
PAFTOOLS="/home/eva/Downloads/minimap2-master/misc/paftools.js"
# GISAID input
G_META=$1
G_FASTA=$2
# DESH input
D_META=$3
D_LINES=$4
D_FASTA=$5
D_EPI=$6
# reference and query sequences
REFERENCE=$7
QUERY=$8
# output dir and log file path
OUTDIR=$9
LOG_FILE="${OUTDIR}/sc2_sewage.log"
#LOG_FILE="${OUTDIR}/sc2_sewage_kallisto.log"

# Check output structure
mkdir -p $OUTDIR
touch $LOG_FILE
exec &> >(tee "$LOG_FILE")



echo "Step 1: Build a reference set"
echo ">>> pre-select based on nonN counts and country"
#python ${SCRIPTPATH}/pipeline_baymlab/preprocess_references.py -gisaid ${G_META} ${G_FASTA} --country "Germany,England" --startdate 2021-02-01 --enddate 2021-04-30 -o ${OUTDIR}/seqs_per_lineage --log ${LOG_FILE}

echo ">>> call variants: find log info in ${OUTDIR}/seqs_per_lineage/<FASTA>.paftools.log"
#${SCRIPTPATH}/pipeline_baymlab/call_variants.sh ${OUTDIR}/seqs_per_lineage ${REFERENCE} ${PAFTOOLS}

echo ">>> select sequences per lineage such that each typical mutation with at least 50% frequency is captured at least once"
mkdir -p ${OUTDIR}/reference_set
python ${SCRIPTPATH}/pipeline_baymlab/select_samples.py -f ${G_FASTA} --vcf ${OUTDIR}/seqs_per_lineage/*_merged.vcf.gz --freq ${OUTDIR}/seqs_per_lineage/*_merged.frq -o ${OUTDIR} --log ${LOG_FILE}



# Step 2: Preprocess fastq data???



#echo "Step 2: Predict variant abundances"
#mkdir -p ${OUTDIR}/kallisto_out
#kallisto index -i ${OUTDIR}/reference_set/sequences.kallisto_idx ${OUTDIR}/reference_set/sequences.fasta
#for sample in ${QUERY}/*;
#do
#  f=${sample##*/}
#  f=${f%.*}
#  kallisto quant -i ${OUTDIR}/reference_set/sequences.kallisto_idx -o ${OUTDIR}/kallisto_out/${f} --single -l 200 -s 20 -t 20 ${sample}
  # summarize lineage abundances
#  python ${SCRIPTPATH}/pipeline_baymlab/output_abundances.py --metadata ${OUTDIR}/seqs_per_lineage/metadata.tsv ${OUTDIR}/kallisto_out/${f}/abundance.tsv --log ${LOG_FILE}
#done
exit 0
