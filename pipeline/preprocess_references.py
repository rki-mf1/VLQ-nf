#!/usr/bin/env python3

import sys
import os
import argparse
import pandas as pd
import glob
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description="Preprocess reference collection: randomly select samples and write into individual files in lineage-specific directories.")
    parser.add_argument('-m, --metadata', dest='metadata', type=str, help="metadata tsv file for full sequence database")
    parser.add_argument('-f, --fasta', dest='fasta_in', type=str, help="fasta file representing full sequence database")
    parser.add_argument('-k', dest='select_k', type=int, default=1000, help="randomly select 1000 sequences per lineage")
    parser.add_argument('--min_len', dest='min_len', type=int, default=29500, help="Don't select sequences with less than the specified minimal number of non-ambiguous nucleotides.")
    parser.add_argument('--country', dest='country', type=str, help="only consider sequences found in specified country")
    parser.add_argument('--state', dest='state', type=str, help="only consider sequences found in specified state")
    parser.add_argument('--seed', dest='seed', default=0, type=int, help="random seed for sequence selection")
    parser.add_argument('-o, --outdir', dest='outdir', type=str, default="seqs_per_lineage", help="output directory")
    args = parser.parse_args()

    # create output directory
    try:
        os.mkdir(args.outdir)
    except FileExistsError:
        pass

    # read metadata
    print('read metadata')
    metadata_df = read_metadata(args.metadata)
    lineages = metadata_df["Pango lineage"].unique()

    # select sequences
    print("select sequences")
    selection_dict = {}
    lineages_with_sequence = []
    for lin_id in lineages:
        # filter for lineage, country and length
        samples = metadata_df.loc[metadata_df["Pango lineage"] == lin_id]
        if args.country:
            samples = samples.loc[samples["Virus name"].str.contains(args.country)]
        if args.state:
            samples = samples.loc[samples["Location"].str.contains(args.state)]
        if args.min_len:
            samples = samples.loc[(1-samples["N-Content"])*samples["Sequence length"] >= args.min_len]
        select_n = min(len(samples), args.select_k)
        if select_n == 0:
            #print("WARNING: no sequences satisfying country and length restrictions for lineage {}".format(lin_id))
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
                gisaid_id = selection["Accession ID"].item()
                seq_name =  selection["seq_name"].item()
                selection_dict[seq_name] = (lin_id, gisaid_id)
            else:
                gisaid_ids = list(selection["Accession ID"])
                seq_names = list(selection["seq_name"])
                for i, seq_name in enumerate(seq_names):
                    gisaid_id =  gisaid_ids[i]
                    selection_dict[seq_name] = (lin_id, gisaid_id)
            lineages_with_sequence.append(lin_id)

    print("Sequences selected for {} lineages".format(len(lineages_with_sequence)))

    # write sequences to separate files
    with open(args.fasta_in, 'r') as f_in:
        keep_line = False
        line_idx = 0
        selection_idx = 0
        for line in f_in:
            if line[0] == '>':
                # sequence identifier
                line_idx += 1
                if line_idx % 100000 == 0:
                    print("{} sequences from input fasta processed".format(line_idx))
                    print("{} sequences from selection found".format(selection_idx))
                seq_id = line.rstrip('\n').lstrip('>')
                if seq_id in selection_dict.keys():
                    lin_id, gisaid_id = selection_dict[seq_id]
                    keep_line = True
                    selection_idx += 1
                    outfile = "{}/{}/{}.fa".format(args.outdir, lin_id, gisaid_id)
                    f_out = open(outfile, 'w')
                    f_out.write(">{}\n".format(seq_id))
                    continue
                    # print("keeping sequence {}".format(seq_id))
                else:
                    # item not found as sequence was not selected
                    # print("not keeping sequence {}".format(seq_id))
                    keep_line = False
            elif keep_line:
                # write nucleotide sequence
                f_out.write(line)
        print("{} sequences from input fasta processed".format(line_idx))
        print("{} sequences from selection found".format(selection_idx))

    # write lineages
    with open("{}/lineages.txt".format(args.outdir), 'w') as f:
        for lin_id in sorted(lineages_with_sequence):
            f.write("{}\n".format(lin_id))

    return


def read_metadata(metadata_file):
    """
    Read metadata from tsv into dataframe
    """
    df = pd.read_csv(metadata_file, sep='\t', header=0, dtype=str)
    # remove samples wich have no pangolin lineage assigned (NaN or None)
    df = df[df["Pango lineage"].notna()]
    df = df[df["Pango lineage"] != "None"]
    # add fasta sequence header
    df['seq_name'] = df[['Virus name', 'Collection date', 'Submission date']].apply(lambda x: '|'.join(x), axis=1)
    # ensure correct data types for calculating number of non-ambiguous bases
    df['Sequence length'] = df['Sequence length'].astype(float)
    df['N-Content'] = df['N-Content'].astype(float)
    return df



if __name__ == "__main__":
    sys.exit(main())
