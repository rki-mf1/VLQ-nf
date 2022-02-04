#!/usr/bin/env python3

##############################################
# Source: https://github.com/baymlab/wastewater_analysis
#############################################

#############################################
# IMPORTS
##############################################
import sys, os, argparse
import pandas as pd
import datetime as dt

################################################
# FUNCTIONS
################################################
def main():
    parser = argparse.ArgumentParser(description="Preprocess GISAID reference collection: clean and format data, if DESH data is provided, remove duplicate EPI_ISL entries.")
    parser.add_argument('-meta, --meta', dest='meta', nargs=1, type=str, help="Specify path to gisaid metadata tsv file")
    parser.add_argument('-epi_map, --epi_map', dest='epi_map', type=str, help="Specify path to csv file containing EPI ISL ids for the desh data set")
    args = parser.parse_args()

    print('read metadata')
    gisaid_metadata = args.meta[0]

    print('process gisaid data')
    gisaid_df = read_filter_gisaid(gisaid_metadata, args.epi_map)
    gisaid_df.to_csv('processed_gisaid_metadata.tsv', sep="\t", index=False)

    return None



def read_filter_gisaid(metadata_file, epi_isl_file):
    """
    Read gisaid metadata from tsv into dataframe, preprocess data
    """
    df = pd.read_csv(metadata_file, sep='\t', header=0, dtype=str)

    print('consider only human samples')
    df = df.loc[df.Host == "Human"]

    if epi_isl_file != None:
        gisaid_ids = pd.read_csv(epi_isl_file, header=None, names=['EPI_ISL'])
        print('remove sample records that are included in desh data set')
        df = df.loc[~df['Accession ID'].isin(gisaid_ids['EPI_ISL'])]

    print("ensure correct data types for calculating number of non-ambiguous bases (Sequence length and N-Content)")
    df['Sequence length'] = df['Sequence length'].astype(float)
    print('remove samples wich have no information on N content')
    df = df[df["N-Content"].notna()]
    df['N-Content'] = df['N-Content'].astype(float)
    df['nonN'] = (1-df['N-Content'])*df['Sequence length']
    df['nonN'] = df['nonN'].astype(int)

    print('remove samples wich have no pangolin lineage assigned (NaN or None)')
    df = df[df["Pango lineage"].notna()]
    df = df[df["Pango lineage"] != "None"]
    df.rename(columns={'Pango lineage':'lineage'}, inplace=True)

    print('adjust date representation in dataframe')
    df["date"] = pd.to_datetime(df["Collection date"], yearfirst=True)

    #print('remove duplicate sequences')
    # TODO: EPI id unique, but fasta header not...that's crazy. mighty keep the sample with lower N content but how to distinguish in fasta later? => tmp solution: drop both duplicates
    #df.drop_duplicates(subset=["Virus name","date","Submission date"],inplace=True,ignore_index=True, keep=False)
    #print('add fasta sequence header for mapping')
    #df['fasta_id'] = df[['Virus name', 'Collection date', 'Submission date']].apply(lambda x: '|'.join(x), axis=1)
    ####### NOTE: removed since currently working on GISAID data from RKI API where fasta id == EPI_ISL id
    df['fasta_id'] = df['Accession ID']
    df.rename(columns={'Accession ID':'record_id'}, inplace=True)

    print('format country')
    df['country'] = df["Location"].apply(lambda x: x.split('/')[1].strip())

    return df[['record_id', 'fasta_id','nonN','lineage','date','country']]



############################################
# MAIN
##########################################

if __name__ == "__main__":
    sys.exit(main())
