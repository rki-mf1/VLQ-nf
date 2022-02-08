# sc2_sewage
Estimate sc2 lineage abundances from wastewater samples

## Summary
The goal is to detect sc2 lineages from wastewater samples and estimate the *abundance* of each lineage within a sample. The method employs the kallisto tool:  sc2 reads are mapped against representative genomes of multiple lineages and assigned one or multiple lineages based on the computed pseudo-alignment. The proportion of reads that were assigned to a certain lineage *X* is interpreted as the abundance of lienage *X*  in the sample.

## Program Composition
tbd
### Analysis Steps
1. Reference set building
- Cleanse input data sets
  - GISAID:
    - select human samples only
    - if DESH data provided, filter out overlapping samples between the two data sets
    - drop samples with no 'N' count or lineage information available
    - drop duplicates by virus name, submission date and collection date (need to construct fasta header to access sequences and there would exist duplicates otherwise)
    - create fasta headers from 'Virus name', 'Collection date' and 'Submission date'
    - create consistent scheme for GISAID and DESH by renaming collection date ('date') and 'Accession ID' ('record\_id')
    - filter by country or continent
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

## Requirements
**Tools:**
* Nextflow (v.?)
* Conda (v.?)

**Data:**
* [GISAID](https://www.epicov.org/epi3/frontend#2a39e0) sequence data (*fasta*) and metadata (*tsv*), all *.xz*-compressed
* [DESH](https://github.com/robert-koch-institut/SARS-CoV-2-Sequenzdaten_aus_Deutschland) lineage data, metadata and sequence data (*csv* and *fasta*, all *.xz*-compressed)
* [Mapping file](https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/DESH/Abfrage-GISAID.pdf?__blob=publicationFile) for DESH and GISAID ids (*tsv* file containing the EPI\_ISL ids of the DESH submissions in GISAID)
* Wastewater reads (*fastq*, *fastq*.gz)
**IMPORTANT: CURRENTLY THE TOOL ONLY RUNS FOR GISAID DATA THAT WAS RETRIEVED FROM THE RKI GIASID API**
This means, in order to run the program on manually downloaded data, the fasta ids in sequence file need to be mapped with the metadata file and rpelaced by the corresponding Accession id which is the EPI\_ISL id.

## Usage
First, create a conda environment ensuring all required packages and versions. Activate and run workflow in the conda env:
```
conda create -n sc2-sewage -c bioconda -c anaconda -c conda-forge bcftools=1.3.1 biopython=1.78 htslib=1.3.1 kallisto=0.46.0 minimap2=2.17 pandas=1.1.3 pyvcf seqtk=1.3 gzip glob2=0.7
conda activate sc2-sewage
```
Execute program like this for using only GISAID data:
```
nextflow run main.nf --gisaid /PATH/TO/GISAID_DATA/ --query PATH/TO/QUERY_FASTQ_FILES/ --country Germany,England --startdate 2021-02-01 --enddate 2021-04-30
```
Execute like this for using both GISAID and DESH as data sources:
```
nextflow run main.nf --desh /PATH/TO/GISAID_DATA/ --desh true --desh_data /PATH/TO/DESH_DATA/ --gisaid_desh_map /PATH/TO/MAPPING/FILE --query PATH/TO/QUERY_FASTQ_FILES/ --country Germany,England --startdate 2021-02-01 --enddate 2021-04-30
```
### Notes on input and query data:
- ``--gisaid`` input should be the path to a folder containing two files:
  1. The sequence fasta file having the substring "fasta" somewhere in its filename
  2. The metadata file having the substring "metadata" somewhere in its filename. The sequence file currently has to be '.gz'.
- ``--desh`` input is handled by a boolean and a path parameter ``--desh_data``:
  - If ``--desh`` is _true_, lineage data, sequence data and metadata are retrieved from the folder located at ``--desh_data``. The  sequence data file should contain the substring "fasta" in its name and currently has to already be decompressed from *.xz* format  and *.gz*-compressed.
  - If ``--desh`` is _false_, only GISAID data is used.
- ``--query`` input should be the path to a folder containing *.fastq* or *fastq.gz* files that should be analysed on the same reference data set lie in the same folder ``PATH/TO/QUERY_FASTQ_FILES/``

### Test data
tbd


## Program Handling
```
nextflow run main.nf --help
```
| Parameter | Type | Description | Default |  
|---	|---	|---	|---	|
| *Mandatory:* |  	|  	|  	|
| --argument |  	|  	|  	|  
| *Optional:* |  	|  	|  	|  
| --argument 	|  	|  	|  	|  
|  *Reference building:*	|  	|  	|  	|  	
| --argument 	|  	|  	|  	|  
|  *Kallisto*	|  	|  	|  	|  	
| --argument 	|  	|  	|  	|  
**Always use full file paths**
## Output
### Results
All output is written to ``results/``:
* ``results/NAME_OF_QUERIED_FASTQ/kallisto_out/`` stores the kallisto output (run info, h5 and _abundance.tsv_) and final output table with the relative  lineage abundances (_predictions.tsv_)
* ``nextflow-autodownload-databases/build_reference/ ``stores variant call files as well as the final reference set (*fasta* and metadata *csv*)

### Execution Metadata
* Process log files are currently written to ``nextflow-run-infos/``.
* NF specific output (e.g. pipeline dag and execution timeline) are also stored in ``nextflow-run-infos/``.

## References
* [Variant abundance estimation for SARS-CoV-2 in 1 wastewater using RNA-Seq quantification](https://www.medrxiv.org/content/10.1101/2021.08.31.21262938v1.full.pdf)
* [GitHub Repository: SARS-CoV-2 variant abundance estimation in wastewater](https://github.com/baymlab/wastewater_analysis)
