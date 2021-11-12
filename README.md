# sc2_sewage
Estiamte sc2 lineage abundances from wastewater samples

# Requirements
* [GISAID](https://www.epicov.org/epi3/frontend#2a39e0) data (fasta) and metadata (tsv)
* [DESH](https://github.com/robert-koch-institut/SARS-CoV-2-Sequenzdaten_aus_Deutschland) metadata (csv), lineage data (csv) and sequence data (fasta)
* [Mapping file](https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/DESH/Abfrage-GISAID.pdf?__blob=publicationFile) for DESH and GISAID ids (csv)
* [SC2 reference sequence](https://www.ncbi.nlm.nih.gov/sars-cov-2/) (fasta)
* Wastewater reads (fastq)
* minimap2 and k8 for using paftools.js for alignment and variant calling: [minimap2, k8, paftools.js](https://github.com/lh3/minimap2) need to be installed and added to PATH (see [call_variants.sh](https://github.com/EvaFriederike/sc2_sewage/tree/main/pipeline_baymlab/call_variants.sh))

# Usage
*Test data:*
- Download sample reads from [Agrawal et al.](https://journals.asm.org/doi/full/10.1128/MRA.00280-21) from [SRA](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=search_seq_name)

*Note:*
- Adjust path to installed paftools.js in [script.sh](https://github.com/EvaFriederike/sc2_sewage/tree/main/script.sh) ($PAFTOOLS variable)
- Use full path for ``WILDTYPE.FASTA``
- more granular parameter setting for building the kallisto reference db currently have to be adjusted in the script
```
./script.sh GISAID_METADATA.TSV GISAID_SEQUENCES.FASTA DESH_METADATA.CSV DESH_LINEAGES.CSV DESH_SEQUENCES.FASTA DESH_EPI_ISL.CSV WILDTYPE.FASTA PATH_TO_QUERY_FASTQ_FILES OUTDIR
```

## Obsolete
* [Frankfurt_lineages](https://github.com/EvaFriederike/sc2_sewage/tree/main/Frankfurt_lineages) contains the lineage and lineage abundance predictions for wastewater reads from Sindlingen und Griesheim using the pipeline from baymlab and GISAID sequences (download 2020-10-01) from Germany.
* [SL](https://github.com/EvaFriederike/sc2_sewage/tree/main/SL_lineages) contains the lineage and lineage abundance predictions for Sindlingen using the pipeline from baymlab and GISAID sequences (download 2020-10-01) from Germany and Germany in December 2020, respectively.
* tbd: lineage abundance predictions for Sindlingen using the pipeline with GISAID and DESH data from decemeber 2020 considering (1) german sequences only and (2) sequences from multiple countries

# References
* [Variant abundance estimation for SARS-CoV-2 in 1 wastewater using RNA-Seq quantification](https://www.medrxiv.org/content/10.1101/2021.08.31.21262938v1.full.pdf)
* [GitHub Repository: SARS-CoV-2 variant abundance estimation in wastewater](https://github.com/baymlab/wastewater_analysis)
