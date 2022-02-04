#!/usr/bin/env python3

##############################################
# Source: https://github.com/baymlab/wastewater_analysis
#############################################

#############################################
# IMPORTS
##############################################
import sys, os, argparse
import pandas as pd
from Bio import SeqIO
import datetime as dt
import gzip
import numpy as np

################################################
# FUNCTIONS
################################################
def main():
    parser = argparse.ArgumentParser(description="Preprocess DESH reference collection: format and clean data.")
    parser.add_argument('-meta, --meta', dest='meta', nargs=1, type=str, help="Specify path to desh sample metadata csv file")
    parser.add_argument('-lineage, --lineage', dest='lineage', nargs=1, type=str, help="Specify path to lineage metadata csv file")
    parser.add_argument('-fasta, --fasta', dest='fasta', nargs=1, type=str, help="Specify path to desh sequence fasta file")
    args = parser.parse_args()

    print('read data')
    sample_metadata = args.meta[0]
    lineage_metadata = args.lineage[0]
    desh_fasta = args.fasta[0]

    print('process desh data')
    metadata_df = read_filter_desh(sample_metadata, lineage_metadata, desh_fasta)
    metadata_df.to_csv('processed_desh_metadata.tsv', sep="\t", index=False)

    return None



def read_filter_desh(sample_metadata, lineage_metadata, desh_fasta):
    """
    Read desh sample and lineage metadata from csv into dataframe, cleanse data frame.
    """
    sample_df = pd.read_csv(sample_metadata, header=0, dtype=str)
    lineage_df = pd.read_csv(lineage_metadata, header=0, dtype=str)
    desh_df = pd.merge(sample_df, lineage_df, on='IMS_ID', how='inner')

    # ensure correct data types for calculating number of non-ambiguous bases
    print('collect non-N counts per sample')
    fasta = SeqIO.to_dict(SeqIO.parse(gzip.open(desh_fasta, 'rt'),'fasta'))
    nonN = []
    for seq_id in desh_df['IMS_ID']:
        if seq_id in fasta:
            seq = fasta[seq_id].seq
            nonN.append(len(seq)-seq.count('N'))
        else:
            # Introduced this when working with data where sequence records were
            # pre-filtered and thus contained less samples than the input metadata tables
            print(f"ATTENTION: {seq_id} is not contained in {desh_fasta}")
            nonN.append(0)
    desh_df['nonN'] = nonN

    print('only select samples with lineage annotation')
    # remove samples wich have no pangolin lineage assigned (NaN or None)
    desh_df = desh_df[desh_df.lineage != 'None']
    desh_df = desh_df[~desh_df.lineage.isna()]

    print('format dates to datetime')
    # adjust date representation in dataframe
    desh_df['date'] = pd.to_datetime(desh_df['DATE_DRAW'], yearfirst=True)

    print('remove duplicate entries based on IMS_ID')
    desh_df.drop_duplicates(subset=['IMS_ID'], inplace=True, keep=False, ignore_index=True)

    print('reformat identification columns (add IMS_ID as fasta id)')
    # rename IMS_ID, add fasta_id column
    desh_df['fasta_id'] = desh_df['IMS_ID']
    desh_df.rename(columns={'IMS_ID':'record_id'}, inplace=True)

    print('add country information')
    desh_df['country'] = ['Germany']*desh_df.shape[0]

    return desh_df[['record_id', 'fasta_id', 'nonN','lineage','date','country']]


############################################
# MAIN
##########################################

if __name__ == "__main__":
    sys.exit(main())
