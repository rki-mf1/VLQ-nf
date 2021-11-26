#!/usr/bin/env python3

##############################################
# Source: https://github.com/baymlab/wastewater_analysis
#############################################

#############################################
# IMPORTS
##############################################
import sys, os, argparse
import glob
import logging
import pandas as pd
from Bio import SeqIO
import datetime as dt

################################################
# FUNCTIONS
################################################
def main():
    parser = argparse.ArgumentParser(description="Preprocess reference collection: randomly select samples and write into individual files in lineage-specific directories.")
    parser.add_argument('-gisaid, --gisaid', dest='gisaid', nargs='*', type=str, help="Specify paths to gisaid metadata tsv file and sequence fasta file")
    parser.add_argument('-desh, --desh', dest='desh', nargs='*', type=str, help="Specify paths to desh sample and lineage metadata csv files and desh sequence fasta file")
    parser.add_argument('-epi, --desh_epi', dest='desh_epi', type=str, help="Specify path to csv file containing EPI ISL ids for the desh data set")
    parser.add_argument('--country', dest='country', type=str, help="only consider sequences found in specified list of comme-separated countries")
    parser.add_argument('--startdate', dest='startdate', type=str, help="only consider sequences found on or after this date; input should be ISO format")
    parser.add_argument('--enddate', dest='enddate', type=str, help="only consider sequences found on or before this date; input should be ISO format")
    parser.add_argument('--min_len', dest='min_len', type=int, default=29500, help="Don't select sequences with less than the specified minimal number of non-ambiguous nucleotides.")
    parser.add_argument('-k', dest='select_k', type=int, default=1000, help="Specify how many samples to randomly select per lineage")
    parser.add_argument('--seed', dest='seed', default=0, type=int, help="random seed for sequence selection")
    parser.add_argument('-o, --outdir', dest='outdir', type=str, default="seqs_per_lineage", help="output directory")
    parser.add_argument('--log', dest='log', type=str, help="path to log file")
    args = parser.parse_args()

    country_filter = [c.strip() for c in args.country.split(',')]

    global logger
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    output_file_handler = logging.FileHandler(args.log, mode='a')
    stdout_handler = logging.StreamHandler(sys.stdout)
    logger.addHandler(output_file_handler)
    logger.addHandler(stdout_handler)

    # create output directory
    try:
        os.mkdir(args.outdir)
    except FileExistsError:
        pass

    # read metadata
    logger.debug('read metadata')
    if args.gisaid !=None and args.desh != None:
        gisaid_metadata = args.gisaid[0]
        gisaid_fasta = args.gisaid[1]
        gisaid_df = read_filter_gisaid(gisaid_metadata, args.desh_epi)
        sample_metadata = args.desh[0]
        lineage_metadata = args.desh[1]
        desh_fasta = args.desh[2]
        desh_df = read_filter_desh(sample_metadata, lineage_metadata, desh_fasta)
        metadata_df = pd.concat([gisaid_df, desh_df], axis=0)
    else:
        if args.gisaid!=None:
            gisaid_metadata = args.gisaid[0]
            gisaid_fasta = args.gisaid[1]
            metadata_df = read_filter_gisaid(metadata_file = gisaid_metadata, epi_isl_file=None)
        if args.desh!=None:
            sample_metadata = args.desh[0]
            lineage_metadata = args.desh[1]
            desh_fasta = args.desh[2]
            metadata_df = read_filter_desh(sample_metadata, lineage_metadata, desh_fasta)
    lineages = metadata_df["lineage"].unique()

    # filter data set by location, date of sample collection, maximum allowed count of unknown bases
    if args.country:
        logger.debug(f'filter by countries: {country_filter}')
        metadata_df = metadata_df.loc[metadata_df["country"].isin(country_filter)]
    if args.startdate:
        logger.debug(f'filter by earliest sampling date: {args.startdate}')
        metadata_df = metadata_df.loc[metadata_df["date"] >= pd.to_datetime(args.startdate)]
    if args.enddate:
        logger.debug(f'filter by latest sampling date: {args.enddate}')
        metadata_df = metadata_df.loc[metadata_df["date"] <= pd.to_datetime(args.enddate)]
    if args.min_len:
        logger.debug(f"filter by minimum non 'N' count: {args.min_len}")
        metadata_df = metadata_df.loc[metadata_df['nonN'] >= args.min_len]
    logger.debug(f"Randomly select {args.select_k} samples per lineage")

    # select sequences for each lineage from filtered metadata
    logger.debug("select sequences")
    selection_dict = {}
    lineages_with_sequence = []
    for lin_id in lineages:
        samples = metadata_df.loc[metadata_df["lineage"] == lin_id]
        select_n = min(len(samples), args.select_k)
        if select_n == 0:
            continue
        else:
            # create lineage directory
            try:
                os.mkdir("{}/{}".format(args.outdir, lin_id))
            except FileExistsError:
                # empty existing directory
                old_files = glob.glob("{}/{}/*".format(args.outdir, lin_id))
                for f_trash in old_files:
                    os.remove(f_trash)
            # randomly select sequences
            selection = samples.sample(n=select_n, random_state=args.seed)
            if select_n == 1:
                seq_id = selection["record_id"].item()
                seq_name =  selection["fasta_id"].item()
                selection_dict[seq_name] = (lin_id, seq_id)
            else:
                seq_ids = list(selection["record_id"])
                seq_names = list(selection["fasta_id"])
                for i, seq_name in enumerate(seq_names):
                    seq_id =  seq_ids[i]
                    selection_dict[seq_name] = (lin_id, seq_id)
            lineages_with_sequence.append(lin_id)
    logger.debug("Sequences selected for {} lineages".format(len(lineages_with_sequence)))

    # write sequences to separate files (structured by lineage)
    logger.debug("searching fasta and writing sequences to output directory...")
    logger.debug("... GISAID data ...")
    n_gisaid = get_sequences(gisaid_fasta, selection_dict, args.outdir)
    if args.desh != None:
        logger.debug("... DESH data ...")
        n_desh = get_sequences(desh_fasta, selection_dict, args.outdir)
    else:
        n_desh = 0
    logger.debug("Total number of selected sequences: {} ".format(n_gisaid+n_desh))

    # write lineages
    with open("{}/lineages.txt".format(args.outdir), 'w') as f:
        for lin_id in sorted(lineages_with_sequence):
            f.write("{}\n".format(lin_id))

    # store filtered metadata
    metadata_df.to_csv(args.outdir+'/metadata.tsv', sep="\t", index=False)

    return None



def read_filter_gisaid(metadata_file, epi_isl_file):
    """
    Read gisaid metadata from tsv into dataframe, preprocess data
    """
    df = pd.read_csv(metadata_file, sep='\t', header=0, dtype=str)

    # consider only human samples
    df = df.loc[df.Host == "Human"]
    if epi_isl_file != None:
        gisaid_ids = pd.read_csv(epi_isl_file, header=None, names=['EPI_ISL'])
        # Remove sample records that are included in desh data set
        df = df.loc[~df['Accession ID'].isin(gisaid_ids['EPI_ISL'])]

    # ensure correct data types for calculating number of non-ambiguous bases
    df['Sequence length'] = df['Sequence length'].astype(float)
    # remove samples wich have no information on N content
    df = df[df["N-Content"].notna()]
    df['N-Content'] = df['N-Content'].astype(float)
    df['nonN'] = (1-df['N-Content'])*df['Sequence length']

    # remove samples wich have no pangolin lineage assigned (NaN or None)
    df = df[df["Pango lineage"].notna()]
    df = df[df["Pango lineage"] != "None"]
    df.rename(columns={'Pango lineage':'lineage'}, inplace=True)

    # adjust date representation in dataframe
    df["date"] = pd.to_datetime(df["Collection date"], yearfirst=True)

    # remove duplicate sequences
    # TODO: EPI id unique, but fasta header not...that's crazy. mighty keep the sample with lower N content but how to distinguish in fasta later? => tmp solution: drop both duplicates
    df.drop_duplicates(subset=["Virus name","date","Submission date"],inplace=True,ignore_index=True, keep=False)
    # add fasta sequence header
    df['fasta_id'] = df[['Virus name', 'Collection date', 'Submission date']].apply(lambda x: '|'.join(x), axis=1)
    df.rename(columns={'Accession ID':'record_id'}, inplace=True)

    # filter by country
    #df = df.loc[df["Location"].apply(lambda x: x.split('/')[1].strip()).isin(country_filter)]
    df['country'] = df["Location"].apply(lambda x: x.split('/')[1].strip())

    return df[['record_id', 'fasta_id','nonN','lineage','date','country']]



def read_filter_desh(sample_metadata, lineage_metadata, desh_fasta):
    """
    Read desh sample and lineage metadata from csv into dataframe, preprocess data
    """
    sample_df = pd.read_csv(sample_metadata, header=0, dtype=str)
    lineage_df = pd.read_csv(lineage_metadata, header=0, dtype=str)
    desh_df = pd.merge(sample_df, lineage_df, on='IMS_ID', how='inner')

    # ensure correct data types for calculating number of non-ambiguous bases
    fasta = SeqIO.to_dict(SeqIO.parse(open(desh_fasta),'fasta'))
    nonN = []
    for seq_id in desh_df['IMS_ID']:
        seq = fasta[seq_id].seq
        nonN.append(len(seq)-seq.count('N'))
    desh_df['nonN'] = nonN

    # remove samples wich have no pangolin lineage assigned (NaN or None)
    desh_df = desh_df[desh_df.lineage != 'None']
    desh_df = desh_df[~desh_df.lineage.isna()]

    # adjust date representation in dataframe
    desh_df['date'] = pd.to_datetime(desh_df['DATE_DRAW'], yearfirst=True)

    # remove duplicate sequences
    desh_df.drop_duplicates(subset=['IMS_ID'], inplace=True, keep=False, ignore_index=True)

    # rename IMS_ID, add fasta_id column
    desh_df['fasta_id'] = desh_df['IMS_ID']
    desh_df.rename(columns={'IMS_ID':'record_id'}, inplace=True)

    desh_df['country'] = ['Germany']*desh_df.shape[0]

    return desh_df[['record_id', 'fasta_id', 'nonN','lineage','date','country']]


# Idea :
# input list of fasta files and iterate or apply function to each fasta file and update selection_dict by removing already visited header id
# let's wait until transforming pipeline into nextflow workflow
def get_sequences(fasta, selection_dict, outdir):
    """
    Write selected sequences from fasta file to separate fasta files in outdir based on
    headers stored in selection_dict
    """
    with open(fasta, 'r') as f_in:
        keep_line = False
        line_idx = 0
        selection_idx = 0
        for line in f_in:
            if line[0] == '>':
                # sequence identifier
                line_idx += 1
                if line_idx % 100000 == 0:
                    logger.debug("{} sequences from input fasta processed".format(line_idx))
                    logger.debug("{} sequences from selection found".format(selection_idx))
                header = line.rstrip('\n').lstrip('>')
                if header in selection_dict.keys():
                    lin_id, seq_id = selection_dict[header]
                    keep_line = True
                    selection_idx += 1
                    outfile = "{}/{}/{}.fa".format(outdir, lin_id, seq_id)
                    f_out = open(outfile, 'w')
                    f_out.write(">{}\n".format(header))
                    continue
                else:
                    # item not found as sequence was not selected
                    keep_line = False
            elif keep_line:
                # write nucleotide sequence
                f_out.write(line.rstrip('\n'))
        logger.debug("{} sequences from input fasta processed".format(line_idx))
        logger.debug("{} sequences from selection found".format(selection_idx))

    return selection_idx

############################################
# MAIN
##########################################

if __name__ == "__main__":
    sys.exit(main())
