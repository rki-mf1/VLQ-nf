#!/usr/bin/env python3

##############################################
# Source: https://github.com/baymlab/wastewater_analysis
#############################################

#############################################
# IMPORTS
##############################################
import sys, os, argparse
import pandas as pd
import vcf

################################################
# FUNCTIONS
################################################
def main():
    parser = argparse.ArgumentParser(description="Build reference set consisting of a selection of samples per pangolin lineage.")
    parser.add_argument('--lineage', required=True, type=str, nargs='+', help="lineage")
    parser.add_argument('--metadata', required=True, type=str, nargs='+', help="sample metadata")
    parser.add_argument('--vcf', required=True, type=str, nargs='+', help="merged vcf file per lineage")
    parser.add_argument('--freq', required=True, type=str, nargs='+', help="merged allele frequency files per lineage")
    parser.add_argument('--min_aaf', default=0.5, type=float, help="minimal alternative allele frequency (AAF) to consider variation")
    parser.add_argument('--max_per_lineage', default=100, type=int, help="select at most k sequences per lineage")
    args = parser.parse_args()

    # Read meta
    # Select reference sequences per pango lineage and write to separate fasta files
    metadata_df = pd.read_csv(args.metadata[0], header=0, sep='\t',dtype={'record_id':str, 'fasta_id':str,'nonN':int,'lineage':str,'date':str,'country':str})
    metadata_df = metadata_df[metadata_df["lineage"]==args.lineage[0]]
    selection_df = select_ref_genomes(args.lineage[0], metadata_df, args.max_per_lineage, args.vcf[0],
                                      args.freq[0], args.min_aaf)
    selection_df.to_csv('aaf_selection.csv', sep='\t', index=False)

    return



def select_ref_genomes(lineage, meta_df, max_per_lineage, vcf_file, freq_file, min_aaf):
    """
    For every pangolin lineage, select sample records such that every characteristic variant with >= min_aaf
    frequency is captured at least once.
    """
    selection_ids = []
    # sort by descending nonN count and actuality
    samples = meta_df.sort_values(by=["nonN", "date"], ascending=False)
    # read allele frequencies and extract sites with AAF >= minimal alt allele frequency
    variant_positions = []
    with open(freq_file, 'r') as f:
        for line in f:
            line = line.split('\t')
            if line[0] == "CHROM":
                continue
            ref_info = line[4]
            ref_allele, freq = ref_info.split(':')
            ref_allele_freq = float(freq)
            alt_allele_freq = 1 - ref_allele_freq
            if alt_allele_freq > min_aaf:
                variant_positions.append(int(line[1]))
    print("{} total # sites with alt allele frequency > {} = {}".format(lineage, min_aaf, len(variant_positions)))

    # if for current lineage no variant position has > min_aaf, keep only one sample for lineage
    selection_count = 0
    if len(variant_positions) == 0:
        selection_ids.append(samples[0])
        selection_count += 1
    else:
        # read vcf and process samples: keep sequences for current lineage such that
        # all mutations with more than min_aaf are captured at least once
        vcf_reader = vcf.Reader(open(vcf_file, 'rb'))
        vcf_samples = vcf_reader.samples
        sample_patterns = {sample : [] for sample in vcf_samples}
        for record in vcf_reader:
            # check if POS shows variant frequency > min_aaf
            if record.POS in variant_positions:
                # get all alleles at POS among the samples in lineage (collect variation)
                alleles = [record.REF] + [str(x) for x in record.ALT]
                for sample in vcf_samples:
                    genotype = record.genotype(sample)['GT']
                    if genotype[0] != '.':
                        allele_idx = int(genotype[0])
                        allele = alleles[allele_idx]
                        sample_patterns[sample].append(allele)
                    else:
                        sample_patterns[sample].append(genotype[0])
        variation_seen = {pos : [] for pos in variant_positions}

        # for each sample in lineage, check if it carries an allele for some relevant variant position that has not been captured yet
        # found one/some? => mark observed allele(s) for that variant position(s) and select sample as representative for genetic variation in lineage
        for sample in vcf_samples:
            select = False
            variation = sample_patterns[sample]
            for i, pos in enumerate(variant_positions):
                allele = variation[i]
                if allele != '.':
                    if allele not in variation_seen[pos]:
                        select = True
                        variation_seen[pos].append(allele)
            if select:
                selection_ids.append(sample)
                selection_count += 1
                if selection_count == max_per_lineage:
                    break
    print("{} sequences selected for lineage {}".format(selection_count,lineage))
    if selection_count == 0:
        print("ERROR: no sequences selected for lineage {}".format(lineage))
        sys.exit(1)

    print("{} sequences selected in total".format(len(selection_ids)))
    selection_df = meta_df.loc[meta_df["record_id"].isin(selection_ids)]

    return selection_df



############################################
# MAIN
##########################################
if __name__ == "__main__":
    sys.exit(main())
