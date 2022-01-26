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
    parser = argparse.ArgumentParser(description="Preprocess reference collection: filter samples by temporal and geographical relevance, randomly select samples from remaining dataset.")
    parser.add_argument('-gisaid, --gisaid', dest='gisaid', type=str, help="Specify paths to gisaid metadata tsv file and sequence fasta file")
    parser.add_argument('-desh, --desh', dest='desh', type=str, help="Specify paths to desh sample and lineage metadata csv files and desh sequence fasta file")
    parser.add_argument('--country', dest='country', nargs='*', type=str, help="only consider sequences found in specified list of comme-separated countries")
    parser.add_argument('--startdate', dest='startdate', nargs='*', type=str, help="only consider sequences found on or after this date; input should be ISO format")
    parser.add_argument('--enddate', dest='enddate', nargs='*', type=str, help="only consider sequences found on or before this date; input should be ISO format")
    parser.add_argument('--min_len', dest='min_len', type=int, default=29500, help="Don't select sequences with less than the specified minimal number of non-ambiguous nucleotides.")
    parser.add_argument('-k', dest='select_k', type=int, default=1000, help="Specify how many samples to randomly select per lineage")
    parser.add_argument('--seed', dest='seed', default=0, type=int, help="random seed for sequence selection")
    args = parser.parse_args()

    country_filter = args.country
    startdate = args.startdate
    enddate = args.enddate

    print('read metadata')
    if args.gisaid != None:
        metadata_df = pd.read_csv(args.gisaid, sep='\t', header=0, dtype={'record_id':str, 'fasta_id':str,'nonN':int,'lineage':str,'date':str,'country':str})
    elif args.desh != None:
        metadata_df = pd.read_csv(args.desh, sep='\t', header=0, dtype={'record_id':str, 'fasta_id':str,'nonN':int,'lineage':str,'date':str,'country':str})
    else:
        gisaid_df = pd.read_csv(args.gisaid, sep='\t', header=0, dtype={'record_id':str, 'fasta_id':str,'nonN':int,'lineage':str,'date':str,'country':str})
        desh_df = pd.read_csv(args.desh, sep='\t', header=0, dtype={'record_id':str, 'fasta_id':str,'nonN':int,'lineage':str,'date':str,'country':str})
        metadata_df = pd.concat([gisaid_df, desh_df], axis=0)
    lineages = metadata_df["lineage"].unique()

    print('filter data set by location, date of sample collection, maximum allowed count of unknown bases')
    print(f'Before filtering the input data comprises {metadata_df.shape[0]} samples and {lineages.shape[0]} lineages')
    if len(country_filter) != 0:
        print(f'filter by countries: {country_filter}')
        metadata_df = metadata_df.loc[metadata_df["country"].isin(country_filter)]
    if len(startdate) != 0:
        startdate = startdate[0]
        print(f'filter by earliest sampling date: {startdate}')
        metadata_df = metadata_df.loc[metadata_df["date"] >= pd.to_datetime(args.startdate)]
    if len(enddate) != 0:
        enddate = enddate[0]
        print(f'filter by latest sampling date: {enddate}')
        metadata_df = metadata_df.loc[metadata_df["date"] <= pd.to_datetime(args.enddate)]
    if args.min_len:
        print(f"filter by minimum non 'N' count: {args.min_len}")
        metadata_df = metadata_df.loc[metadata_df['nonN'] >= args.min_len]

############################### TODO
    print(f"Randomly select {args.select_k} samples per lineage")
    # select sequences for each lineage from filtered metadata
    selection_df = pd.DataFrame(columns=metadata_df.columns)
    lineages_with_sequence = []
    for lin_id in lineages:
        samples = metadata_df.loc[metadata_df["lineage"] == lin_id]
        select_n = min(len(samples), args.select_k)
        if select_n == 0:
            continue
        else:
            # randomly select sequences
            selection = samples.sample(n=select_n, random_state=args.seed)
            selection_df = pd.concat([selection_df, selection], axis=0)
            lineages_with_sequence.append(lin_id)
    print("Sequences selected for {} lineages".format(len(lineages_with_sequence)))

    # TODO: do we want to log how many samples were selected from GISAID and DESH, respectively?

    print("Total number of selected sequences: {} ".format(selection_df.shape[0]))

    # write lineages
    with open("lineages.csv", 'w') as f:
        for lin_id in sorted(lineages_with_sequence):
            f.write("{}\n".format(lin_id))

    # store filtered metadata
    selection_df.to_csv('selection_by_metadata.csv', sep="\t", index=False)

    return None



############################################
# MAIN
##########################################

if __name__ == "__main__":
    sys.exit(main())
