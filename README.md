# **PaliDIS 2.0** - **Pali**ndromic **D**etection of **I**nsertion **S**equences

PaliDIS is a Nextflow pipeline that predicts insertion sequence annotations of paired-end, short-read metagenomic data.

The tool is based upon identifying inverted terminal repeats (ITRs) (figure below) of insertion sequences using an efficient maximal exact matching algorithm between reads of large metagenomic datasets, which are then resolved on metagenomic assemblies.

<img src="img/insertion_sequence.png" alt="insertion sequence" width="400"/>

## Installation
- Install [Nextflow](https://www.nextflow.io/)
- Install [Docker](https://www.docker.com/) ([Singularity](https://sylabs.io/singularity/) if using a HPC)
- Git clone this repo
```bash
git clone https://github.com/blue-moon22/Palidis.git
cd Palidis
```

## Usage

The pipeline is divided into two workflows: **1) Get Candidate ITRs** and **2) Get IS Annotations**.

If you are running this on an HPC, you will need to specify `-profile <executor>` in the command. Currently, the pipeline only supports `LSF` (using `-profile lsf`). It is possible to add another profile to the [nextflow config](https://www.nextflow.io/docs/latest/config.html) to make this pipeline compatible with other HPC executors. If you do so, you are welcome to folk this repo and make a pull request to include your new profile for others to use.

### 1. Get Candidate ITRs
#### Usage
This workflow generates candidate ITRs for each sample.
```bash
nextflow palidis.nf --get_candidate_itrs --manifest <manifest_file> --batch_name <batch_name> -resume
```
The output is stored in a directory in the current run directory specified with `--batch_name`.

A tab-delimited manifest must be specified for `--manifest` containing the absolute paths with headers `lane_id`, `read1`, `read2`, `sample_id` and `assembly_path`, e.g. this manifest contains three samples (the first having two lanes and the other two having one lane):

lane_id | read1 | read2 | sample_id | assembly_path
:---: | :---: | :---: | :---: | :---:
lane1 | /path/to/file/lane1_1.fq.gz | /path/to/file/lane1_2.fq.gz | my_sample | /path/to/file/contigs.fasta
lane2 | /path/to/file/lane2_1.fq.gz | /path/to/file/lane2_2.fq.gz | my_sample1 | /path/to/file/my_sample1_contigs.fasta
lane3 | /path/to/file/lane3_1.fq.gz | /path/to/file/lane3_2.fq.gz | my_sample2 | /path/to/file/my_sample2_contigs.fasta
lane4 | /path/to/file/lane4_1.fq.gz | /path/to/file/lane4_2.fq.gz | my_sample3 | /path/to/file/my_sample3_contigs.fasta

#### Output
This workflow generates FASTA files with clipped reads representing candidate ITRs and tab-delimited files containing the ITR-containing (non-clipped) reads' positional information for each sample within the `batch_name` directory, e.g.:

sample_id | contig | read | position
:---: | :---: | :---: | :---:
my_sample1 | contig_name1 | read_name_with_itr1 | 14
my_sample1 | contig_name2 | read_name_with_itr2 | 15
my_sample1 | contig_name3 | read_name_with_itr3 | 88
my_sample1 | contig_name4 | read_name_with_itr5 | 24


### 2. Get IS Annotations
#### Usage
Once you are satisfied all samples have been processed through the first workflow, you can run the second. This workflow pools, clusters and filters the candidate ITRs to create a non-redundant catalogue of ITRs. The output file prefix must be specified by `--output_prefix`, alongside the `--batch_name` directory.

```bash
nextflow palidis.nf --get_IS_annotations --batch_name <batch_name> --output_prefix <output_file_prefix> -resume
```

You can include another ITR catalogue into this workflow by specifying `--include_IR_db /path/to/other_clipped_reads_db.fasta`. This will create a more comprehensive ITR catalogue, which will generate a better representation of rarer ITR clusters.

#### Output
A non-redundant catalogue of ITRs called `<output_file_prefix>_clipped_reads_db.fasta` and a final tab-delimited file of insertion sequence annotations called `<output_file_prefix>_insertion_sequence_annotations.tab` combining all samples are generated. The annotation file consists of the sample_id, contig name, start and end positions of the first ITR, start and end positions of the second ITR and the cluster(s) they belong to, e.g.:

sample_id | contig | itr1_start_position | itr1_end_position | itr2_start_position | itr2_end_position | itr_clusters
:---: | :---: | :---: | :---: | :---: | :---: | :---:
sample_id1 | contig_name1 | 29 | 43 | NA | NA | 1217817
sample_id2 | contig_name2 | 23 | 43 | 2769 | 2822 | 1217817;656079
sample_id3 | contig_name3 | 29 | 43 | NA | NA | 1217817

Although two flanking ITRs may be found, it is possible that positions could not be predicted (represented by `NA`). (This happens when a read maps to a contig, but its paired read containing the ITR does not.)
