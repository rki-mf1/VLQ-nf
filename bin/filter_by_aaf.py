#!/usr/bin/env python3

##############################################
# Source: https://github.com/baymlab/wastewater_analysis
# This is the final aaf filtering version we use for the benchmarking study
# Note: - Filter by population aaf
#       - generate minimal sequence set to capture every filtered variant at least once
#       - increase reference index for each lineage by adding filtered sequences until the maximum allowe number of sequences per lineage is met
#############################################

#############################################
# IMPORTS
##############################################
import sys, os, argparse
import pandas as pd
import vcf
import random
from glob import glob

################################################
# FUNCTIONS
################################################
def main():
    parser = argparse.ArgumentParser(description="Build reference set consisting of a selection of characteristic sequences per pangolin lineage.")
    parser.add_argument('--metadata', required=True, type=str, nargs='+', help="sample metadata")
    parser.add_argument('--vcf', required=True, type=str, nargs='+', help="merged vcf file per lineage")
    parser.add_argument('--max_per_lineage', default=0, type=int, help="maximum number of sequences per lineage")
    parser.add_argument('--min_aaf', default=0.5, type=float, help="minimal alternative allele frequency (AAF) to consider relevant variation")
    args = parser.parse_args()
    random.seed(0)
    print("read data")
    vcf_list = [ glob(file)[0] for file in args.vcf]
    metadata_df = pd.read_csv(args.metadata[0], header=0, sep='\t',dtype={'record_id':str, 'nonN':int,'lineage':str,'date':str,'country':str, 'continent':str})
    random.seed(0)
    selection_df = select_ref_genomes(metadata_df, vcf_list, args.min_aaf, args.max_per_lineage)
    selection_df.to_csv('aaf_selection.csv', sep='\t', index=False)

    return



def select_ref_genomes(meta_df, vcf_list, min_aaf, max_per_lineage):
    """
    For every pangolin lineage, select samples such that every characteristic variant with
    population aaf > min_aaf is captured at least once.
    """
    lineages = meta_df.lineage.unique()
    vcf_dict = {vcf.split('/')[-1].split('_')[0] : vcf for vcf in vcf_list}
    seqs_per_lineage = {}
    minimal_set = []
    repertoire_per_lineage = {}
    lineage_without_passing_variants = {}
    for lin_id in lineages:
        print(f"Begin aaf filtering for lineage {lin_id}")
        try:
            vcf_file = vcf_dict[lin_id]
        except KeyError as e:
            print("WARNING: skipping lineage {}, VCF info missing".format(lin_id))
            continue

        lin_samples = meta_df.loc[meta_df["lineage"] == lin_id]
        lin_samples = lin_samples.sort_values(by=["nonN", "date", "record_id"], ascending=False)
        vcf_reader = vcf.Reader(open(vcf_file, 'rb'))
        vcf_samples = vcf_reader.samples
        N = len(vcf_samples)
        variant_sample_dict = {'POS:GT':[], 'record_id': []}
        seqs_per_lineage[lin_id] = []
        characteristic_alleles = 0

        # for every variant, read the sequences containing it, store in dict
        # variant_sample_dict = {'POS:GT': [..], 'record_id': [...]}
        for record in vcf_reader:
            for sample in record.samples:
                if sample['GT'] != './.':
                    variant_sample_dict['POS:GT'].append(str(record.POS)+':'+sample['GT'][0])
                    variant_sample_dict['record_id'].append('_'.join(sample.sample.split('_')[1:]))
        # compressed_variant_sample_dict = {'12:1': [id1, id2, ..., idn], '140:2': [id2, ..., idm], ... }
        compressed_variant_sample_dict = {v: [] for v in set(variant_sample_dict['POS:GT'])}
        for var,seq_id in list(zip(variant_sample_dict['POS:GT'],variant_sample_dict['record_id'])):
            # only keep records that passed previous filtering steps
            if seq_id in lin_samples['record_id'].values:
                compressed_variant_sample_dict[var].append(seq_id)
        # Filter variants by population aaf across sequences
        filtered_variant_dict = {k: v for k,v in compressed_variant_sample_dict.items() if len(v)/float(N)>min_aaf}
        filtered_variant_dict_seqs = list(set([seq for l in filtered_variant_dict.values() for seq in l]))
        characteristic_alleles = len(filtered_variant_dict)
        print(f"{lin_id} total # sites with allele frequency > {min_aaf} = {characteristic_alleles}")
        print(f"# sequences available = {len(filtered_variant_dict_seqs)}")

        # select minimal reference
        if characteristic_alleles != 0:
            variant_sample_df = pd.DataFrame.from_dict(variant_sample_dict)
            variant_sample_df = variant_sample_df.loc[variant_sample_df['record_id'].isin(filtered_variant_dict_seqs)]
            variant_sample_df = variant_sample_df.loc[variant_sample_df['POS:GT'].isin(filtered_variant_dict.keys())]
            if len(variant_sample_df) != 0:
                df_for_minimal_selection = variant_sample_df
                captured_vars = []
                while len(captured_vars) < characteristic_alleles:
                    seq_order = df_for_minimal_selection['record_id'].value_counts().index.values
                    seq_selected = seq_order[0]
                    seqs_per_lineage[lin_id].append(seq_selected)
                    print(f"Added {seq_selected} to {lin_id} reference")
                    seq_vars = df_for_minimal_selection.loc[df_for_minimal_selection['record_id'] == seq_selected]['POS:GT'].values
                    captured_vars.extend(seq_vars)
                    print(f"Captured {len(seq_vars)} variants")
                    df_for_minimal_selection = df_for_minimal_selection.loc[df_for_minimal_selection['record_id']!=seq_selected]
                    df_for_minimal_selection = df_for_minimal_selection[~df_for_minimal_selection['POS:GT'].isin(captured_vars)]
                    print(f"Remaining variants to cover: {df_for_minimal_selection['POS:GT'].unique()} ")
                min = len(seqs_per_lineage[lin_id])
                minimal_set.append(min)
                print(f"{lin_id} minimum # sequences necessary to capture every important variant at least once = {min}")
                repertoire_per_lineage[lin_id] = variant_sample_df.loc[~variant_sample_df['record_id'].isin(seqs_per_lineage[lin_id])]
                print(f"Added {len(repertoire_per_lineage[lin_id]['record_id'].unique())} seq to {lin_id} repertoire")
            else:
                print(f"ATTENTION: No sequences selected for {lin_id} due to missing overlap between characteristic variants and prefiltered sequence data")
        else:
            lineage_without_passing_variants[lin_id] = lin_samples['record_id'].values
            print(f"ATTENTION: No sequences selected for {lin_id} due to no variants passing the AAF filter")
    # increase reference up to the largest minimal reference or the pre-determined maximum #seqs
    max_no_seqs = max(max(minimal_set), max_per_lineage)
    for k in lineage_without_passing_variants.keys():
        v = list(lineage_without_passing_variants[k])
        seqs_per_lineage[k] = random.sample(v,max_no_seqs)
    print(f"Select up to {max_no_seqs} samples per lineage")
    for lin_id,var_dict in repertoire_per_lineage.items():
        diff = max_no_seqs - len(seqs_per_lineage[lin_id])
        print(f"Select {diff} more sequences for {lin_id}")
        if diff > 0:
            available_seqs = list(var_dict['record_id'].unique())
            if len(available_seqs) >= diff:
                backup = random.sample(available_seqs, diff)
            else:
                backup = available_seqs
            seqs_per_lineage[lin_id].extend(backup)
            print(f"Successsfully added {len(backup)} more seqs for {lin_id}")
        print(f"{lin_id} # selected sequences = {len(seqs_per_lineage[lin_id])}")

    selected_ids = list(set([i for id_list in seqs_per_lineage.values() for i in id_list]))
    selection_df = meta_df.loc[meta_df["record_id"].isin(selected_ids)]

    return selection_df



############################################
# MAIN
##########################################
if __name__ == "__main__":
    sys.exit(main())
