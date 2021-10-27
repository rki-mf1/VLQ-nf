# sc2_sewage
Estiamte sc2 lineage abundances from wastewater samples

# Not included in repository
* GISAID data (fasta) and metadata (tsv)
* Wastewater reads (fastq): Download [sample reads](https://journals.asm.org/doi/full/10.1128/MRA.00280-21) from [SRA](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=search_seq_name)
* minimap2 and k8 for using paftools.js for alignment and variant calling: [minimap2, k8, paftools.js](https://github.com/lh3/minimap2) need to be installed and added to PATH (see [call_variants.sh](https://github.com/EvaFriederike/sc2_sewage/tree/main/pipeline_baymlab/call_variants.sh))


# First results
[SL](https://github.com/EvaFriederike/sc2_sewage/tree/main/SL_lineages) contains the lineage and lineage abundance predictions for Sindlingen using the pipeline from baymlab and GISAID sequences (download 2020-10-01) from Germany and Germany in December 2020, respectively.

## Obsolete
[Frankfurt_lineages](https://github.com/EvaFriederike/sc2_sewage/tree/main/Frankfurt_lineages) contains the lineage and lineage abundance predictions for wastewater reads from Sindlingen und Griesheim using the pipeline from baymlab and GISAID sequences (download 2020-10-01) from Germany.

# Usage
*Note:*Adjust path to installed paftools.js in script.sh ($PAFTOOLS)
```
./script.sh METADATA.TSV FASTA.TSV WILDTYPE.FASTA PATH_TO_QUERY_FASTQ_FILES OUTDIR
```

# References
* [Variant abundance estimation for SARS-CoV-2 in 1 wastewater using RNA-Seq quantification](https://www.medrxiv.org/content/10.1101/2021.08.31.21262938v1.full.pdf)
* [GitHub Repository: SARS-CoV-2 variant abundance estimation in wastewater](https://github.com/baymlab/wastewater_analysis)
