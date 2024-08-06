#!/usr/bin/env python3

###############################################################
### Source: https://github.com/baymlab/wastewater_analysis
###############################################################

################################################################
### IMPORTS
################################################################
import sys, os, argparse
import pandas as pd
import datetime as dt

###########################################
### FUNCTIONS
###########################################
def main():
    parser = argparse.ArgumentParser(description="Preprocess reference collection: filter samples by temporal and geographical relevance, randomly select samples from remaining dataset.")
    parser.add_argument('-gisaid, --gisaid', dest='gisaid', nargs=1, type=str, help="Specify paths to gisaid metadata tsv file and sequence fasta file")
    parser.add_argument('--country', dest='country', nargs='*', type=str, help="only consider sequences found in specified list of comma separated countries")
    parser.add_argument('--continent', dest='continent', nargs='*', type=str, help="only consider sequences found in specified list of comma separated continents")
    parser.add_argument('--startdate', dest='startdate', nargs='*', type=str, help="only consider sequences found on or after this date; input should be ISO format")
    parser.add_argument('--enddate', dest='enddate', nargs='*', type=str, help="only consider sequences found on or before this date; input should be ISO format")
    parser.add_argument('--min_len', dest='min_len', type=int, default=29500, help="Don't select sequences with less than the specified minimal number of non#ambiguous nucleotides.")
    parser.add_argument('-k', dest='select_k', type=int, default=1000, help="Specify how many samples to randomly select per lineage")
    parser.add_argument('--seed', dest='seed', default=0, type=int, help="random seed for sequence selection")
    args = parser.parse_args()

    print('read data')
    metadata_df = pd.read_csv(args.gisaid[0], sep='\t', header=0, dtype={'record_id':str, 'nonN':int,'lineage':str,'date':str,'country':str, 'continent':str})

    lineages = metadata_df["lineage"].unique()
    metadata_df['date'] = pd.to_datetime(metadata_df['date'])

    print('filter samples by geographical origin, date of collection, maximum allowed count of unknown bases in genomic sequence')
    print(f'--- Before filtering the input data comprises {metadata_df.shape[0]} samples and {lineages.shape[0]} lineages')
    if len(args.continent) != 0:
        continent_filter = args.continent[0].split(',')
        print(f'--- filter by continent: {continent_filter}')
        subset_df_continent = metadata_df.loc[metadata_df["continent"].isin(continent_filter)]
    if len(args.country) != 0:
        country_filter = args.country[0].split(',')
        print(f'--- filter by countries: {country_filter}')
        subset_df_country = metadata_df.loc[metadata_df["country"].isin(country_filter)]
    if len(args.continent) != 0:
        if len(args.country) != 0:
            metadata_df = pd.concat([subset_df_continent, subset_df_country], axis=0)
            metadata_df.drop_duplicates(subset=['record_id'], inplace=True, ignore_index=True, keep='first')
        else:
            metadata_df = subset_df_continent
    elif len(args.country) != 0:
        metadata_df = subset_df_country

    if len(args.startdate) != 0:
        startdate = args.startdate[0]
        print(f'--- filter by earliest sampling date: {startdate}')
        metadata_df = metadata_df.loc[metadata_df["date"] >= pd.to_datetime(startdate)]
    if len(args.enddate) != 0:
        enddate = args.enddate[0]
        print(f'--- filter by latest sampling date: {enddate}')
        metadata_df = metadata_df.loc[metadata_df["date"] <= pd.to_datetime(enddate)]
    if args.min_len:
        print(f"--- filter by minimum non 'N' count: {args.min_len}")
        metadata_df = metadata_df.loc[metadata_df['nonN'] >= args.min_len]
    print(f"--- Remaining {metadata_df.shape[0]} samples and {len(metadata_df.lineage.unique())} lineages")

    print(f"--- Randomly select {args.select_k} samples per lineage")
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
    print("--- Sequences selected for {} lineages".format(len(lineages_with_sequence)))
    print("--- Total number of selected sequences: {} ".format(selection_df.shape[0]))

    # store filtered metadata
    selection_df.to_csv('selection_by_metadata.csv', sep="\t", index=False)

    return None



#############################################################
### MAIN
#############################################################
if __name__ == "__main__":
    sys.exit(main())
