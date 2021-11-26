#!/usr/bin/env python3

##############################################
# Source: https://github.com/baymlab/wastewater_analysis
#############################################

#############################################
# IMPORTS
##############################################
import sys, os, argparse
import logging
import pandas as pd
import vcf

################################################
# FUNCTIONS
################################################
def main():
    parser = argparse.ArgumentParser(description="Build reference set consisting of a selection of samples per pangolin lineage.")
    parser.add_argument('-f, --fasta', dest='fasta_in', nargs='+', type=str, help="fasta files representing full gisaid and desh sequence databases")
    parser.add_argument('--vcf', required=True, type=str, nargs='+', help="vcf files per lineage")
    parser.add_argument('--freq', required=True, type=str, nargs='+', help="allele frequency files per lineage")
    parser.add_argument('--min_aaf', default=0.5, type=float, help="minimal alternative allele frequency (AAF) to consider variation")
    parser.add_argument('--max_per_lineage', default=100, type=int, help="select at most k sequences per lineage")
    parser.add_argument('-o, --outdir', dest='outdir', type=str, default='.', help="output directory")
    parser.add_argument('--log', dest='log', type=str, help="path to log file")
    args = parser.parse_args()

    global logger
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    output_file_handler = logging.FileHandler(args.log, mode='a')
    stdout_handler = logging.StreamHandler(sys.stdout)
    logger.addHandler(output_file_handler)
    logger.addHandler(stdout_handler)

    # Create output directory
    try:
        os.mkdir(args.outdir)
    except FileExistsError:
        pass

    # Read meta
    # Select reference sequences per pango lineage and write to separate fasta files
    metadata_df = pd.read_csv(args.outdir+'/seqs_per_lineage/metadata.tsv', sep='\t')
    selection_df = select_ref_genomes(metadata_df, args.max_per_lineage, args.vcf,
                                      args.freq, args.min_aaf)

    # Filter collectino of input fasta sequences according to selection and write new fasta
    fasta_out = args.outdir + "/reference_set/sequences.fasta"
    filter_fasta(args.fasta_in, fasta_out, selection_df)
    logger.debug("Selected sequences written to {}".format(fasta_out))

    return



def select_ref_genomes(metadata_df, max_per_lineage, vcf_list, freq_list, min_aaf):
    """
    For every pangolin lineage, select sample records such that every characteristic variant with >= 50%
    frequency is captured at least once.
    """
    # check which lineages are present
    lineages = metadata_df["lineage"].unique()
    logger.debug("# lineages = {}".format(len(lineages)))
    # assign vcfs and allele frequency files to lineages, assuming vcfs are in current directory and named after the corresponding lineage
    vcf_dict = {vcf.split('/')[-1].split('_')[0] : vcf for vcf in vcf_list}
    freq_dict = {fname.split('/')[-1].split('_')[0] : fname for fname in freq_list}
    # select samples for every lineage
    selection_ids = []
    for lin_id in lineages:
        samples = metadata_df.loc[metadata_df["lineage"] == lin_id]
        # sort by descending nonN count and actuality
        samples = samples.sort_values(by=["nonN", "date"], ascending=False)
        # read allele frequencies and extract sites with AAF >= minimal alt allele frequency
        try:
            allele_freq_file = freq_dict[lin_id]
        except KeyError as e:
            logger.debug("WARNING: skipping lineage {}, allele frequency info missing".format(lin_id))
            continue
        # collect variant positions with more than min_aaf allele frequency
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
        logger.debug("{} total # sites with alt allele frequency > {} = {}".format(
                lin_id, min_aaf, len(variant_positions)))

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
                logger.debug("WARNING: skipping lineage {}, VCF info missing".format(lin_id))
                continue
            vcf_reader = vcf.Reader(open(vcf_file, 'rb'))
            samples = vcf_reader.samples
            sample_patterns = {sample : [] for sample in samples}
            for record in vcf_reader:
                # check if POS shows variant frequency > min_aaf
                if record.POS in variant_positions:
                    # get all alleles at POS among the samples in lineage (collect variation)
                    alleles = [record.REF] + [str(x) for x in record.ALT]
                    for sample in samples:
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
            for sample in samples:
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
        logger.debug("{} sequences selected for lineage {}".format(selection_count,
                                                            lin_id))
        if selection_count == 0:
            logger.debug("ERROR: no sequences selected for lineage {}".format(lin_id))
            sys.exit(1)

    logger.debug("{} sequences selected in total".format(len(selection_ids)))
    selection_df = metadata_df.loc[
                        metadata_df["record_id"].isin(selection_ids)]

    return selection_df


# Idea:
# store selected sequences in a desh and gisaid fasta respectively (kallisto should take multiple input fastas)
# Let's wait until transforming pipeline into nextflow workflow
def filter_fasta(fasta_in, fasta_out, selection_df):
    """
    Filter fasta according to selected metadata
    """
    keep_line = False
    selection_identifiers = selection_df["fasta_id"].unique()
    f_out = open(fasta_out, 'w')
    for fasta in fasta_in:
        with open(fasta, 'r') as f_in:
            for line in f_in:
                if line[0] == '>':
                    # sequence identifier
                    seq_id = line.rstrip('\n').lstrip('>')
                    if seq_id in selection_identifiers:
                        # kallisto doesn't accept spaces in fasta headers
                        if '/' in seq_id:
                            header = line.split('/')
                            header[1] = header[1].replace(' ','_')
                            line = '/'.join(header)
                        f_out.write(line)
                        keep_line = True
                    else:
                        keep_line = False
                elif keep_line:
                    # nucleotide sequence
                    f_out.write(line)
    f_out.close()

    return



############################################
# MAIN
##########################################
if __name__ == "__main__":
    sys.exit(main())
