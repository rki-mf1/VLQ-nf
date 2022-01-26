
# sc2_sewage
Estiamte sc2 lineage abundances from wastewater samples


# Requirements
* [GISAID](https://www.epicov.org/epi3/frontend#2a39e0) sequence data (fasta) and metadata (tsv)
* [Mapping file](https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/DESH/Abfrage-GISAID.pdf?__blob=publicationFile) for DESH and GISAID ids (csv containing the EPI\_ISL ids of the DESH submissions in GISAID)
* Wastewater reads (fastq)


# Usage
**Test data:**
- Download sample reads from communities around Frankfurt by [Agrawal et al.](https://journals.asm.org/doi/full/10.1128/MRA.00280-21) from [SRA](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=search_seq_name) (see section _Data Availability_). Required data structure: all fastq files lie in the same folder ``PATH/TO/QUERY_FASTQ_FILES/``
- ``nf_test/`` contains small test files to show that the workflow functionally works. Both metadata and sequence files contain the first 100 entries of a larger gisaid dataset. For application with larger files and valdiation the issue with insufficient RM/memory has to be fixed.

First, create a conda environment ensuring all required packages and versions. Activate and run workflow in the conda env:
```
conda create -n sc2-sewage -c bioconda -c anaconda bcftools=1.3.1 biopython=1.78 htslib=1.3.1 kallisto=0.46.0 minimap2=2.17 pandas=1.1.3 pyvcf=0.68
conda activate sc2-sewage
```
```
nextflow run main.nf --gisaid nf_test/ --query PATH/TO/QUERY_FASTQ_FILES/ --gisaid_desh_map PATH/TO/MAPPING_CSV
```

# Parameter Handling
```
nextflow run main.nf --help
```
* **TODO**: give table

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
## Results
All output is written to ``results/``:
* ``results/NAME_OF_QUERIED_FASTQ/kallisto_out/`` stores the kallisto output (run info, h5 and _abundance.tsv_) and final output table with the relative  lineage abundances (_predictions.tsv_)
* ``nextflow-autodownload-databases/build_reference/ ``stores variant call files as well as the final reference set (fasta and metadata csv)

## Meta
* Process log files are currently written to ``nextflow-run-infos/``.
* NF specific output (e.g. pipeline dag and execution timeline) are also stored in ``nextflow-run-infos/``.

# References
* [Variant abundance estimation for SARS-CoV-2 in 1 wastewater using RNA-Seq quantification](https://www.medrxiv.org/content/10.1101/2021.08.31.21262938v1.full.pdf)
* [GitHub Repository: SARS-CoV-2 variant abundance estimation in wastewater](https://github.com/baymlab/wastewater_analysis)
