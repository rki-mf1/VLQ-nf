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
from glob import glob

################################################
# FUNCTIONS
################################################
def main():
    parser = argparse.ArgumentParser(description="Build reference set consisting of a selection of samples per pangolin lineage.")
    parser.add_argument('--metadata', required=True, type=str, nargs='+', help="sample metadata")
    parser.add_argument('--vcf', required=True, type=str, nargs='+', help="merged vcf file per lineage")
    parser.add_argument('--freq', required=True, type=str, nargs='+', help="merged allele frequency files per lineage")
    parser.add_argument('--min_aaf', default=0.5, type=float, help="minimal alternative allele frequency (AAF) to consider variation")
    parser.add_argument('--max_per_lineage', default=100, type=int, help="select at most k sequences per lineage")
    args = parser.parse_args()


    print("read data")
    vcf_list = [ glob(file)[0] for file in args.vcf]
    freq_list = [ glob(file)[0] for file in args.freq]
    metadata_df = pd.read_csv(args.metadata[0], header=0, sep='\t',dtype={'record_id':str, 'nonN':int,'lineage':str,'date':str,'country':str, 'continent':str})

    print(f"select genomic sequences for every lineage such that every characteristic variant with frequency>={args.min_aaf} is captured at least once")
    selection_df = select_ref_genomes(metadata_df, args.max_per_lineage, vcf_list, freq_list, args.min_aaf)
    selection_df.to_csv('aaf_selection.csv', sep='\t', index=False)

    return



def select_ref_genomes(meta_df, max_per_lineage, vcf_list, freq_list, min_aaf):
    """
    For every pangolin lineage, select sample records such that every characteristic variant with >= min_aaf
    frequency is captured at least once.
    """
    # check which lineages are present
    lineages = meta_df["lineage"].unique()
    print(f"{len(lineages)} lineages")
    # assign vcfs and allele frequency files to lineages, assuming vcfs are in current directory and named after the corresponding lineage
    vcf_dict = {vcf.split('/')[-1].split('_')[0] : vcf for vcf in vcf_list}
    freq_dict = {fname.split('/')[-1].split('_')[0] : fname for fname in freq_list}
    selection_ids = []
    for lin_id in lineages:
        print(f"Begin aaf filtering for lineage {lin_id}")
        lin_samples = meta_df.loc[meta_df["lineage"] == lin_id]
        # sort by descending nonN count and actuality
        lin_samples = lin_samples.sort_values(by=["nonN", "date", "record_id"], ascending=False)
        print(lin_samples)
        # read allele frequencies and extract sites with AAF >= minimal alt allele frequency
        try:
            allele_freq_file = freq_dict[lin_id]
        except KeyError as e:
            print("WARNING: skipping lineage {}, allele frequency info missing".format(lin_id))
            continue

        # collect variant positions with more than min_aaf allele frequency
        # read allele frequencies and extract sites with AAF >= minimal alt allele frequency
        variant_positions = []
        with open(allele_freq_file, 'r') as f:
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
        print("{} total # sites with alt allele frequency > {} = {}".format(lin_id, min_aaf, len(variant_positions)))

        # if for current lineage no variant position has > min_aaf, keep only one sample for lineage
        selection_count = 0
        if len(variant_positions) == 0:
            selection_ids.append(samples[0])
            selection_count += 1
        else:
            # read vcf and process samples: keep sequences for current lineage such that
            # all mutations with more than min_aaf are captured at least once
            try:
                vcf_file = vcf_dict[lin_id]
            except KeyError as e:
                print("WARNING: skipping lineage {}, VCF info missing".format(lin_id))
                continue
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
            name_map = {}
            for sample in sample_patterns.keys():
                name_map['_'.join(sample.split('_')[1:])] = sample
            lin_samples = lin_samples.loc[lin_samples["record_id"].isin(name_map.keys())]
            for sample in lin_samples["record_id"]:
                select = False
                variation = sample_patterns[name_map[sample]]
                print(f"sample: {sample}")
                print(f"variation pattern:{variation}")
                for i, pos in enumerate(variant_positions):
                    allele = variation[i]
                    if allele != '.':
                        if allele not in variation_seen[pos]:
                            select = True
                            variation_seen[pos].append(allele)
                if select:
                    selection_ids.append(sample)
                    selection_count += 1
                    print(f"Selected sample {sample}, current count of representative sequences for {lin_id}: {selection_count}")
                    if selection_count == max_per_lineage:
                        break
                else:
                    print(f"Did not select sample {sample}")
        print("--- {} sequences selected for lineage {}".format(selection_count,lin_id))
        if selection_count == 0:
            print("ERROR: no sequences selected for lineage {}".format(lin_id))
            #sys.exit(1)

    print("--- {} sequences selected in total".format(len(selection_ids)))
    selection_df = meta_df.loc[meta_df["record_id"].isin(selection_ids)]

    return selection_df



############################################
# MAIN
##########################################
if __name__ == "__main__":
    sys.exit(main())
