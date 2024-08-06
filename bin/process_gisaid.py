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
    args = parser.parse_args()

    print('read data')
    gisaid_metadata = args.meta[0]

    print('process gisaid metadata')
    gisaid_df = read_filter_gisaid(gisaid_metadata)
    gisaid_df.to_csv('processed_gisaid_metadata.tsv', sep="\t", index=False)

    return None



def read_filter_gisaid(metadata_file):
    """
    Read gisaid metadata from tsv into dataframe, preprocess data
    """
    df = pd.read_csv(metadata_file, sep='\t', header=0, dtype=str)
    x = df.shape[0]

    print('consider only human samples')
    df = df.loc[df.Host == "Human"]

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

    print('remove samples wich have no location information (NaN or None)')
    df = df[df["Location"].notna()]
    df = df[df["Location"] != "None"]
    print("replace whitespace with underscore")
    df['continent'] = df["Location"].apply(lambda x: x.split('/')[0].strip().replace(' ','_'))

    print('adjust date representation in dataframe')
    df["date"] = pd.to_datetime(df["Collection date"], yearfirst=True)

    print('remove duplicate sequences')
    # TODO: EPI id unique, but fasta header not...that's crazy. mighty keep the sample with lower N content but how to distinguish in fasta later? => tmp solution: drop both duplicates
    df = df.sort_values(by='nonN', ascending=False)
    df = df.drop_duplicates(subset=["Accession ID"],ignore_index=True)
    #print('add fasta sequence header for mapping')
    #df['fasta_id'] = df[['Virus name', 'Collection date', 'Submission date']].apply(lambda x: '|'.join(x), axis=1)
    ####### NOTE: removed since currently working on GISAID data from RKI API where fasta id == EPI_ISL id
    #df['fasta_id'] = df['Accession ID']
    print("use accession id as record id")
    df.rename(columns={'Accession ID':'record_id'}, inplace=True)

    print('retrieve country information, replace whitespace by underscore')
    df['country'] = df["Location"].apply(lambda x: x.split('/')[1].strip().replace(' ','_'))

    print(f'--- Remaining samples: {df.shape[0]}/{x}')

    return df[['record_id', 'nonN','lineage','date','country', 'continent']]



############################################
# MAIN
##########################################

if __name__ == "__main__":
    sys.exit(main())
