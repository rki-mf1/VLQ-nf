#!/usr/bin/env python3

#############################################
# IMPORTS
##############################################
import sys, os, argparse
import pandas as pd
from Bio import SeqIO

################################################
# FUNCTIONS
################################################
def main():
    parser = argparse.ArgumentParser(description="Split multifasta and collect lineage annotation")
    parser.add_argument('-meta, --meta', dest='meta', nargs=1, type=str, help="Specify path to desh sample metadata csv file")
    parser.add_argument('-multifasta, --multifasta', dest='multifasta', nargs=1, type=str, help="Specify path to multifasta file")
    args = parser.parse_args()
    sample_metadata = args.meta[0]
    multifasta = args.multifasta[0]

    print('read data')
    metadata_df = pd.read_csv(sample_metadata, sep='\t', header=0, dtype={'record_id':str, 'fasta_id':str,'nonN':int,'lineage':str,'date':str,'country':str})

    print('parse multifasta records, get lineage and write to fasta')
    for record in SeqIO.parse(open(multifasta),'fasta'):
        lineage = metadata_df[metadata_df['fasta_id'] == record.id].lineage.item()
        dir_name = "fasta/"+lineage
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
        filename = lineage+"_"+str(record.id)+".fasta"
        with open("fasta/"+filename, "w") as f:
            SeqIO.write(record, f , "fasta")

    return None


############################################
# MAIN
##########################################

if __name__ == "__main__":
    sys.exit(main())
