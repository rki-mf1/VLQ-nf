# sc2_sewage
Estiamte sc2 lineage abundances from wastewater samples.
The implemented pipeline is based on the [manuscript](https://www.medrxiv.org/content/10.1101/2021.08.31.21262938v1.full.pdf) and [code](https://github.com/baymlab/wastewater_analysis) by Baaijens, Zulli, Ott et al..
The code in this branch represents the pipeline at an early development stage before we started transforming it into a Nextflow workflow. This branch does not implement the same pipeline as the Nextflow branch.

# Requirements
* [GISAID](https://www.epicov.org/epi3/frontend#2a39e0) data (fasta) and metadata (tsv)
* [DESH](https://github.com/robert-koch-institut/SARS-CoV-2-Sequenzdaten_aus_Deutschland) metadata (csv), lineage data (csv) and sequence data (fasta)
* [Mapping file](https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/DESH/Abfrage-GISAID.pdf?__blob=publicationFile) for DESH and GISAID ids (csv containing the EPI\_ISL ids of the DESH submissions in GISAID)
* [SC2 reference sequence](https://www.ncbi.nlm.nih.gov/sars-cov-2/) (fasta)
* Wastewater reads (fastq)
* minimap2 and k8 for using paftools.js for alignment and variant calling: [minimap2, k8, paftools.js](https://github.com/lh3/minimap2) need to be installed and added to PATH: 
 
add minimap2, k8 temporarily to $PATH via 

```export PATH="$PATH:`pwd`:`pwd`/misc" ```

**Note:** Despite $PATH update, I had to specify the complete file path to paftools.js additionally in line 8 of [script.sh](https://github.com/EvaFriederike/sc2_sewage/tree/main/script.sh) ($PAFTOOLS variable)

# Usage
**Test data:**
- Download sample reads from communities around Frankfurt by [Agrawal et al.](https://journals.asm.org/doi/full/10.1128/MRA.00280-21) from [SRA](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=search_seq_name) (see section _Data Availability_)
- Required data structure: all fastq files lie in the same folder ``PATH_TO_QUERY_FASTQ_FILES``

**Note:**
- Use full path for ``WILDTYPE.FASTA``
- For input paths to directories (e.g. ``OUTDIR``) loose the last slash (e.g. _example/path/to/outdir_)
- Even if one would like to only use GISAID data to build the reference, with the current implementation still an input for all parameter positions of _script.sh_ are needed (annoying, will change soon)
- More granular parameter setting for building the kallisto reference db currently have to be adjusted in the script (e.g. if only GISAID data should be used, just loose the ``-gisaid`` parameter  for _preprocess\_references.py_)

```
./script.sh GISAID_METADATA.TSV GISAID_SEQUENCES.FASTA DESH_METADATA.CSV DESH_LINEAGES.CSV DESH_SEQUENCES.FASTA DESH_EPI_ISL.CSV WILDTYPE.FASTA PATH_TO_QUERY_FASTQ_FILES OUTDIR
```

# Analysis Steps
1. Reference set building
- Cleanse input data sets
  - GISAID:
    - select human samples only
    - if DESH data provided, filter out overlapping samples between the two data sets
    - drop samples with no 'N' count or lineage information available
    - drop duplicates by virus name, submission date and collection date (need to construct fasta header to access sequences and there would exist duplicates otherwise)
    - create fasta headers from 'Virus name', 'Collection date' and 'Submission date'
    - create consistent scheme for GISAID and DESH by renaming collection date ('date') and 'Accession ID' ('record\_id')
    - filter by country
  - DESH:
    - get non 'N' count from sequence data
    - filter out samples with no lineage information available
    - drop duplicates by 'IMS\_ID'
    - create consistent scheme for GISAID and DESH by renaming collection date ('date') and 'IMS\_ID' ('record\_id')
- Filter records by sample collection date
- Remove sequences with less than X non-ambiguous calls
- Sample k sequences to represent each lineage

2. Call variants
3. Select final reference set from variant data
- Filter samples per lineage such that each variant with freq >50% is represented at least once
3. Apply kallisto to prepared reference set and input fastq reads
4. Summarize predicted lineage abundances

# Output
All output is written to ``OUTDIR``:
* sc2\_sewage.log
* seqs\_per\_lineage contains information on all sequences and lineages that were selected in the first and second analysis steps:
  * a folder for each lineage that samples from the input data sources were selected for. Each lineage folder contains selected samples assigned with that lineage (fasta) and the individual variant call output files (paftools and vcf files)
  * for each lineage one file with merged variant call results across all selected samples plus frequency files for the called alternate alleles
  * _lineages.txt_ containing all covered lineages
* reference\_set stores the final selected kallisto reference after the third analysis step as well as the kallisto index
* kallisto\_out stores the final output. Each analysed fastq query gets a subfolder containig the kallisto output (run info, h5 and _abundance.tsv_) and the derived lineage abundances (_predictions.tsv_)


## Obsolete
* [Frankfurt_lineages](https://github.com/EvaFriederike/sc2_sewage/tree/main/Frankfurt_lineages) contains the lineage and lineage abundance predictions for wastewater reads from Sindlingen und Griesheim using the pipeline from baymlab and GISAID sequences (download 2020-10-01) from Germany.
* [SL](https://github.com/EvaFriederike/sc2_sewage/tree/main/SL_lineages) contains the lineage and lineage abundance predictions for Sindlingen using the pipeline from baymlab and GISAID sequences (download 2020-10-01) from Germany and Germany in December 2020, respectively.
* tbd: lineage abundance predictions for Sindlingen using the pipeline with GISAID and DESH data from decemeber 2020 considering (1) german sequences only and (2) sequences from multiple countries

# References
* [Variant abundance estimation for SARS-CoV-2 in 1 wastewater using RNA-Seq quantification](https://www.medrxiv.org/content/10.1101/2021.08.31.21262938v1.full.pdf)
* [GitHub Repository: SARS-CoV-2 variant abundance estimation in wastewater](https://github.com/baymlab/wastewater_analysis)
