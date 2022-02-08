#!/usr/bin/env python3

##############################################
# Source: https://github.com/baymlab/wastewater_analysis
#############################################

###########################
# IMPORTS
###########################
import sys, os, argparse, logging
import pandas as pd

###########################
# FUNCTIONS
###########################
def main():
    parser = argparse.ArgumentParser(description="Plot abundances from file.")
    parser.add_argument('abundances', type=str, help="abundance file")
    parser.add_argument('--metadata', type=str, help="metadata file")
    parser.add_argument('-m', dest='min_ab', type=float, default=0, help="minimal frequency (%) to output variant")
    parser.add_argument('--voc', dest='voc', nargs='*', type=str, help="comma-separated list of strains of interest, output abundance for these only")
    parser.add_argument('-o', dest='outfile', type=str, help="write output to tsv file")
    args = parser.parse_args()
    outfile = args.outfile

    if len(args.voc) != 0:
        vocs = args.voc[0].split(',')
        print(f"--- Predicted abundances summarized for {vocs}")
    else:
        vocs = args.voc
        print(f"--- Predicted abundances summarized for all lineages")

    print("read abundance predictions and metadata")
    if args.metadata:
        df = pd.read_csv(args.metadata, sep='\t', header=0, dtype={'record_id':str, 'fasta_id':str,'nonN':int,'lineage':str,'date':str,'country':str, 'continent':str})

    print("get abundances per reference sequence, group by lineage and summarize quantification")
    abundance_dict = {}
    abundance_format = ""
    with open(args.abundances, 'r') as f:
        # for each matched reference sequence, add the quantified abundance of
        # mapped reads to the overall abundance of the correpsonding lineage
        for line in f:
            line = line.rstrip('\n').split('\t')
            if line[-1] == "tpm":
                # kallisto header
                abundance_format = "kallisto"
                continue
            elif line[-1] == "NumReads":
                # salmon header
                abundance_format = "salmon"
                continue
            if abundance_format == "":
                print("ERROR: abundance file format not recognized as kallisto or salmon")
                sys.exit(1)
            # previuosly replaced blank spaces in country name since kallisto doesn't accept spaces in fasta headers
            seqname = line[0]
            if '/' in seqname:
                header = seqname.split('/')
                header[1] = header[1].replace('_',' ')
                seqname = '/'.join(header)
            if args.metadata:
                variant = df.loc[df["fasta_id"] == seqname]["lineage"]
                variant = variant.iloc[0]
            else:
                variant = seqname
            if abundance_format == "kallisto":
                tpm = float(line[-1])
            else:
                tpm = float(line[-2])
            abundance = tpm / 10**6
            if variant in abundance_dict:
                abundance_dict[variant][0] += tpm
                abundance_dict[variant][1] += abundance
            else:
                abundance_dict[variant] = [tpm, abundance]

    print("compute corrected abundances and filter for input VOCs")
    total_ab = sum([v[1] for v in abundance_dict.values()])
    with open(outfile, 'w') as f:
        f.write("** evaluating {}\n".format(args.abundances))
        f.write("** {}\n".format(' '.join(sys.argv)))
        f.write("* variant\ttpm\tfreq(%)\tadj_freq(%)\n")
        for variant, values in abundance_dict.items():
            tpm, ab = values
            corrected_ab = ab / total_ab
            if ab >= args.min_ab / 100:
                if len(vocs) == 0 or variant in vocs:
                    f.write("{}\t{:.0f}\t{:.2f}\t{:.2f}\n".format(
                            variant, tpm, ab * 100, corrected_ab * 100))

    return


###########################
# MAIN
###########################
if __name__ == "__main__":
    sys.exit(main())
