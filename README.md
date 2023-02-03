
# VLQ-nf

Estimates SARS-CoV-2 lineage abundances from wastewater samples. Here, we provide a Nextflow pipeline implementation of the original idea by Baaijens, Zulli, Ott _et al_. Please check their [manuscript](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02805-9) and [code](https://github.com/baymlab/wastewater_analysis) and acknowledge their work and cite their paper when using anything from this repository! 

## Summary

The goal is to detect SARS-CoV-2 lineages from wastewater samples and estimate the abundance of each lineage within a sample. The implemented method repurposes the [Kallisto](https://pachterlab.github.io/kallisto) tool like first presented by [Baaijens, Zulli, Ott _et al_](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02805-9). In a first step, a reference index containing multiple sequences for multiple SARS-CoV-2 lineages is reconstructed from [GISAID](https://gisaid.org/) data. In the second step, a Kallisto index is built and FASTQ samples are pseudo-aligned to estimate the abundance of each reference lineage in the samples.

## Program Composition
* main.nf: main script called by user
* Reference building:
    | Process | Description |
    |--- |--- |
   | process_desh.nf | Cleanse input data set and create consistent scheme: human samples containing information origin, samplin date and lineage annotation |
    | process_gisaid.nf | Cleanse input data set and create consistent scheme: human samples containing information origin, samplin date and lineage annotation. Samples already contained in provided DESH data are ignored |
    | filter_by_metadata.nf | Filter records by country or continent, sampling date, number of ambiguous bases. Sample X samples per lineage. |
    | filter_sequences.nf| Extract sequence data for selected samples from input FASTA. |
    | variant_call.nf| Perform variant calling on selected sequences. |
    | merge_vcf.nf| Merge VCF data per lineage. |
    | filter_by_aaf.nf| Filter samples for each lineage, so that those variants that were observed in more than X% of the variant called sequences are represented at least once. Sample selection is controlled by a cutoff. |
    | filter_sequences_by_aaf.nf| Extact sequence data for samples selected in AAF filtering. |
* Abundance estimatation:
*   | Process | Description |
    |--- |--- |
    | build_index.nf | Build kallisto index structure from the final sequence data set |
    | kallisto_prediction_(single,paired).nf | Quantify lineage abundances per input sample |




## Requirements
This tool was developed and tested on Nextflow (v.21.10.6) and Conda (v.4.12.0).

## Usage
First, create a conda environment ensuring all required packages and versions. Activate and run workflow in the conda env:

```bash
conda create -n sc2-sewage -c bioconda -c anaconda -c conda-forge bcftools=1.3.1 biopython=1.78 htslib=1.3.1 kallisto=0.46.0 minimap2=2.17 pandas=1.1.3 pyvcf seqtk=1.3 gzip glob2=0.7
conda activate sc2-sewage
```

Execute program like this for using only GISAID data:

```bash
nextflow run main.nf --gisaid /PATH/TO/GISAID_DATA/ --query test/ --country Germany,England --startdate 2021-02-01 --enddate 2021-04-30
```

**UNDER CONSTRUCTION:**
*Execute like this for using both GISAID and DESH as data sources:*

```bash
nextflow run main.nf --desh /PATH/TO/GISAID_DATA/ --desh true --desh_data /PATH/TO/DESH_DATA/ --gisaid_desh_map /PATH/TO/MAPPING/FILE --query test/ --country Germany,England --startdate 2021-02-01 --enddate 2021-04-30
```

## Data:
* [**GISAID**](https://www.epicov.org/epi3/frontend#2a39e0): download data set with sequence data (FASTA, FASTA.GZ) and metadata (tab-delimited CSV, Required columns: ['Virus name', 'Accession ID', 'Collection date', 'Location', 'Sequence length', 'Host', 'Pango lineage', 'N-Content'])
**ATTENTION**: The FASTA headers nee to be the EPI_ISL Accession IDs.
* **UNDER CONSTRUCTION:** *[DESH](https://github.com/robert-koch-institut/SARS-CoV-2-Sequenzdaten_aus_Deutschland) lineage data, metadata and sequence data (TSV and FASTA OR FASTA.GZ)*
* **UNDER CONSTRUCTION:** *[Mapping file](https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/DESH/Abfrage-GISAID.pdf?__blob=publicationFile) for DESH and GISAID ids (TSV file containing the EPI\_ISL ids of the DESH submissions in GISAID)*
* **Wastewater samples** (FASTQ, FASTQ.GZ)

## Program Handling
```
nextflow run main.nf --help
```
| Parameter | Type | Description | Default |  
|---	|---	|---	|---	|
| *Mandatory:* |  	|  	|  	|
| - -gisaid | String | Full system path to folder storing metadata and sequence files from gisaid for reference building. <br> Currently, metadata has to be tab-delimited CSV or TSV file and sequence info has to be FASTA or FASTA.GZ. <br>**Required columns:** ['Virus name', 'Accession ID', 'Collection date', 'Location', 'Sequence length', 'Host', 'Pango lineage', 'N-Content'] | None |  
| or - -reference | String	|  If true, the prebuilt reference in --databases is used to build an index and analyse the query | false |  
| - -query | String | Path to folder storing query FASTQ or FASTQ.GZ files | None |  
| *Optional:* |  	|  	|  	|  
| - -desh | Boolean| !!! UNDER CONSTRUCTION - NOT GUARANTEED TO WORK WITH CURRENT VERSION !!! <br>If true, program includes DESH data into reference reconstruction | false |
| - -desh_data | String | !!! UNDER CONSTRUCTION - NOT GUARANTEED TO WORK WITH CURRENT VERSION !!! <br>Full system path to folder storing metadata, lineage data and sequence data from gisaid<br> **metadata filenames:** SARS-CoV-2-Sequenzdaten_Deutschland.csv, SARS-CoV-2-Entwicklungslinien_Deutschland.csv, <br>**sequence data filename:** "\*fasta*" (FASTA or FASTA.GZ) | None |
| - -gisaid_desh_map | String | !! UNDER CONSTRUCTION - NOT GUARANTEED TO WORK WITH CURRENT VERSION !!! <br>Full system path to CSV file mapping accession ids of overlapping samples between desh and gisaid. If desh input is provided, it is recommended to also input such a mapping file to avoid duplicates in the reference set! | None |  
| *Reference building:* |  	|  	|  	|  	
| - -continent	| List of Strings	| List of comma-separated continents to consider when selecting reference samples based on geography. Inform yourself how continents are captured by your input data resource. Whitespaces in country names need to be replaced by underscore. <br>Example: "North_America" | " " |  
| - -country | List of Strings | List of comma-separated countries to consider when selecting reference samples based on geography. Inform yourself how country names are captured by GISAID. Whitespaces in country names need to be replaced by underscore.<br> Example: Bosnia_Herzegovina | None |  
| - -startdate | String | Earliest sampling date YYYY-MM-DD to consider when selecting reference samples based on the timepoint of sample drawing. | " " |  
| - -enddate | String | Latest sampling date YYYY-MM-DD to consider when selecting reference samples based on the timepoint of sample drawing. | " " |  
| - -min_len | Int	| Don't select sequences with less than the specified minimal number of non-ambiguous nucleotides. | 29500 |  
| - -k | Int | Specify how many samples to randomly select per lineage when filtering based on metadata. | 1000 |  
| - -seed | Int | Random seed for sampling steps during sequence selection  | 0 |  
| - -min_aaf | Float | Minimal alternative allele frequency (AAF) to consider for representative lineage variation. | 0.5 |  
| - -max_per_lineage | Int | Maximum number of lineages to select to represent a lineages genomic variation. | 0 |  
|  *Kallisto*	|  	|   |  	|  	
|- -single_end| boolean | If true, run kallisto parametrized for single-end reads, else for paired-end reads | true |
| - -fragment_length	| Int | Estimated average fragment length. | 200 |
| - -fragment_length_sd | Int | Estimated standard deviation of fragment length. | 20 |  
| - -kallisto_threads | Int | Number of threads to use. | 20 |  
| - -bootstrap | int | Determine whether kallisto should be run in bootstrap mode and how many bootstraps should be used. | 0 |
|  *Output*	|  	|  	|  	|  	
| - -min_ab | Float | Summarize output for all lineages whose estimated abundance is above this minimum threshold. | 0 |
|  *Nextflow options*	|  	|  Paths of listed folders need to be relative to location of main.nf |  	|  	
| -output | String	| path to output folder | sc2-sewage-results/ |
| -runinfo | String	| path to folder storing log files | sc2-sewage-runinfo/	|
| -databases | String | path to folder storing intermediate files and reference data | sc2-sewage-databases/	|

For more information on kallisto, e.g. regarding single and paired-end fastq queries, see https://pachterlab.github.io/kallisto/manual

## Output
### Results
``--output``:
* ``NAME_OF_QUERY/kallisto_out/`` stores the kallisto output fiels and the final output table containing the lineage abundances (``predictions.tsv``)

``--databases``:
* ``build_reference/ ``stores variant call files as well as the final reference set (FASTA and CSV)

``--runinfo``:
* Process-specific log files
* Nextflow tracking files like the pipeline DAG, execution report and timeline.

## Example
The ``example`` directory contains a small example to test the pipeline.
A toy example of reference sequences is given in ``resource/sequences.fa`` with metadata in
``resource/metadata.tsv``. This example is taken from https://github.com/baymlab/wastewater_analysis. In order to work with this implementation, we adjusted the metdata scheme and fasta headers according to the required data structure.
```
nextflow run main.nf --gisaid example/resource --query example/samples --k 10 --single_end false --output example/results --runinfo example/runinfo --databases example/databases
```
The predictions can be found in ``kallisto/predictions.tsv`` and should show that 100% of this sample is SARS-CoV-2.

## References
* [Variant abundance estimation for SARS-CoV-2 in wastewater using RNA-Seq quantification](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02805-9)
* [Original GitHub repository: SARS-CoV-2 variant abundance estimation in wastewater](https://github.com/baymlab/wastewater_analysis)
