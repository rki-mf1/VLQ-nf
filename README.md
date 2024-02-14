# VLQ-nf

Estimates SARS-CoV-2 lineage abundances from wastewater samples. Here, we provide a Nextflow pipeline implementation of the original idea by Baaijens, Zulli, Ott _et al_. Please check their [manuscript](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02805-9) and [code](https://github.com/baymlab/wastewater_analysis) and acknowledge their work and cite their paper when using anything from this repository!

## Summary

The goal is to detect SARS-CoV-2 lineages from wastewater samples and estimate the abundance of each lineage within a sample. The implemented method repurposes the [Kallisto](https://pachterlab.github.io/kallisto) tool like first presented by [Baaijens, Zulli, Ott _et al_](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02805-9). In a first step, a reference data set containing multiple sequences for multiple SARS-CoV-2 lineages is reconstructed from [GISAID](https://gisaid.org/) data. In the second step, a Kallisto index is built and FASTQ samples are pseudo-aligned to estimate the abundance of each reference lineage in the samples.

## Program Composition

### 1. Reference building

| Process | Description |
| --- |--- |
| process_gisaid.nf | Cleanse input data set and create consistent scheme: human samples containing information on sample origin, date and lineage annotation.|
| filter_by_metadata.nf | Filter records by geographic origin, sampling date, number of ambiguous bases. Sample X genomic sequences per lineage. |
| filter_sequences.nf| Extract sequence data for selected samples from input FASTA. |
| variant_call.nf| Call variants. |
| merge_vcf.nf| Merge VCF data per lineage. |
| filter_by_aaf.nf| Filter samples for each lineage, so that those variants that were observed in more than Y% of the variant called sequences are represented at least once. |
| filter_sequences_by_aaf.nf| Extact sequence data for samples selected in AAF-based filtering. |

### 2. Abundance estimatation

| Process | Description |
|--- |--- |
| build_index.nf | Build kallisto index structure from the final reference data set |
| kallisto_prediction_(single,paired).nf | Quantify lineage abundances per input sample |

## Requirements
This tool was developed and tested on Nextflow (v.21.10.6) and Conda (v.4.12.0).

## Usage
First, create a conda environment ensuring all required packages and versions. Activate and run workflow in the conda env:

```bash
conda create -n sc2-sewage -c bioconda -c anaconda -c conda-forge bcftools=1.3.1 biopython=1.78 htslib=1.3.1 kallisto=0.46.0 minimap2=2.17 pandas=1.1.3 pyvcf seqtk=1.3 gzip glob2=0.7
conda activate sc2-sewage
```

Execute program like this to reconstruct a reference data set and using kallisto to analyse your query samples:

- be sure to always use full file paths
- check which parameters (e.g. output locations) might have default values set in nextflow.config
```bash
nextflow run main.nf --gisaid PATH_TO_GISAID_DATA/ --query PATH_TO_SEQ_DATA/ --country Germany,England --startdate 2021-02-01 --enddate 2021-04-30
```

If you want to reuse an already processed reference data set, use like this:

- be aware, that existing reference csv and fasta files are looked up in the directory defined by `--databases`
```bash
nextflow run main.nf --reference true --query PATH_TO_SEQ_DATA/
```

**Default settings for all parameters can be found in nextflow.config**

### Example
The ``example/`` directory contains a small example from https://github.com/baymlab/wastewater_analysis to test the pipeline.
The examplary reference metadata and sequence files are stored in ``example/resource/``, the wastewater sequencing data in `example/samples/`. In order to work with this implementation, we adjusted the metdata scheme and fasta headers according to the required data structure.
```
nextflow run main.nf --gisaid PATH_TO/example/resource --query PATH_TO/example/samples --k 10 --single_end false --output PATH_TO/example/results --runinfo PATH_TO/example/runinfo --databases PATH_TO/example/databases
```
The predictions can be found in ``example/results/input/kallisto_out/predictions.tsv`` and should show that 100% of this sample is SARS-CoV-2.

## Data
* [**GISAID**](https://www.epicov.org/epi3/frontend#2a39e0): download data set with 
	* sequence data (fasta, fasta.gz) and 
	* metadata (tab-delimited csv) required columns: ['Virus name', 'Accession ID', 'Collection date', 'Location', 'Sequence length', 'Host', 'Pango lineage', 'N-Content'])
* **Wastewater samples** (fastq, fastq.gz)

## Program Handling
```
nextflow run main.nf --help
```
| Parameter | Type | Description | Default |
|---    |---    |---    |---    |
| **Mandatory:** 					||||
| - -gisaid | String | Full system path to folder storing metadata and sequence files from gisaid for reference building. <br> Currently, metadata has to be tab-delimited csv or tsv file and sequence info has to be fastaor fasta.gz. <br>**Required columns:** ['Virus name', 'Accession ID', 'Collection date', 'Location', 'Sequence length', 'Host', 'Pango lineage', 'N-Content'] | None |
| or - -reference | String      |  If true, the prebuilt reference in --databases is used to build an index and analyse the query | false |
| - -query | String | Path to folder storing query fastq or fastq.gz files | None |
| **Reference building:** |       |       |       |
| - -continent  | List of Strings       | List of comma-separated continents to consider when selecting reference samples based on geography. Inform yourself how continents are captured by your input data resource. Whitespaces in country names need to be replaced by underscore. <br>Example: "North_America" | " " |
| - -country | List of Strings | List of comma-separated countries to consider when selecting reference samples based on geography. Inform yourself how country names are captured by GISAID. Whitespaces in country names need to be replaced by underscore.<br> Example: Bosnia_Herzegovina | None |
| - -startdate | String | Earliest sampling date YYYY-MM-DD to consider when selecting reference samples based on the timepoint of sample drawing. | " " |
| - -enddate | String | Latest sampling date YYYY-MM-DD to consider when selecting reference samples based on the timepoint of sample drawing. | " " |
| - -min_len | Int      | Don't select sequences with less than the specified minimal number of non-ambiguous nucleotides. | 29500 |
| - -k | Int | Specify how many samples to randomly select per lineage when filtering based on metadata. | 1000 |
| - -seed | Int | Random seed for sampling steps during sequence selection  | 0 |
| - -min_aaf | Float | Minimal alternative allele frequency (AAF) to consider for representative lineage variation. | 0.5 |
| - -max_per_lineage | Int | Maximum number of lineages to select to represent a lineages genomic variation. | 0 |
|  **Kallisto**   |       |   |   |
|- -single_end| boolean | If true, run kallisto parametrized for single-end reads, else for paired-end reads | true |
| - -fragment_length    | Int | Estimated average fragment length. | 200 |
| - -fragment_length_sd | Int | Estimated standard deviation of fragment length. | 20 |
| - -kallisto_threads | Int | Number of threads to use. | 20 |
| - -bootstrap | int | Determine whether kallisto should be run in bootstrap mode and how many bootstraps should be used. | 0 |
|  **Output**     |       |       |       |
| - -min_ab | Float | Summarize output for all lineages whose estimated abundance is above this minimum threshold. | 0 |
|  **Nextflow options**   |       |  Paths of listed folders need to be relative to location of main.nf |         |
| -output | String      | path to output folder | sc2-sewage-results/ |
| -runinfo | String     | path to folder storing log files | sc2-sewage-runinfo/        |
| -databases | String | path to folder storing intermediate files and reference data | sc2-sewage-databases/    |

For more information on kallisto, e.g. regarding single and paired-end fastq queries, see https://pachterlab.github.io/kallisto/manual

## Output
### Results
``--output``:
* ``NAME_OF_QUERY_SAMPLE/kallisto_out/`` stores the kallisto output fiels and the final output table containing the lineage abundances (``predictions.tsv``)

``--databases``:
* ``build_reference/`` stores variant call files as well as the final reference set (fasta and csv)

``--runinfo``:
* Process-specific log files
* Nextflow tracking files like the pipeline DAG, execution report and timeline.


## References
* [Variant abundance estimation for SARS-CoV-2 in wastewater using RNA-Seq quantification](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02805-9)
* [Original GitHub repository: SARS-CoV-2 variant abundance estimation in wastewater](https://github.com/baymlab/wastewater_analysis)
