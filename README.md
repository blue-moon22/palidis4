[![GitHub release (latest by date)](https://img.shields.io/github/v/release/blue-moon22/palidis)](https://github.com/blue-moon22/palidis/releases)
![GitHub Workflow Status](https://img.shields.io/github/workflow/status/blue-moon22/palidis/test)

<img src="img/logo.png" alt="logo" width="400"/>

# **PaliDIS** - **Pali**ndromic **D**etection of **I**nsertion **S**equences
## Introduction

PaliDIS is a Nextflow pipeline that quickly discovers novel insertion sequences.

The tool is based upon identifying inverted terminal repeats (ITRs) (figure below) using paired-end, short-read metagenomic data.

For each sample, the pipeline produces two output files: **1. FASTA file of insertion sequences** and **2. Information for each insertions sequence**

<img src="img/insertion_sequence.png" alt="insertion sequence" width="400"/>

## Pipeline summary
**Steps:**
1. Pre-process FASTQ.GZ reads [`convertToFasta`]
2. Efficient maximal exact matching to get repeat sequences using [pal-MEM](https://github.com/blue-moon22/pal-MEM) [`palmem`]
3. Map reads against assemblies using Bowtie2 [`filterContigs` `buildDB` `mapreads`]
4. Get candidate ITRs by distance filters [`getCandidateITRs`]
5. Cluster candidate ITRs using CD-HIT-EST [`clusterReads`]
6. Get putative ITRs by cluster concordance and output Insertion Sequences [`getITRs`]
7. Search against ISfinder [`buildBLASTDB` `buildBLASTDB` `searchISfinder`]
8. _Optional:_ Search against a COB index to predict IS origin [`searchCOBSIndex`]
7. Combine ISfinder and optional COB index search results [`getISInfoWithCOBS` `getISInfoWithoutCOBS`]

## Installation on HPC
- Install [Nextflow](https://www.nextflow.io/)
- Install [Docker](https://www.docker.com/) if using own machine or install [Singularity](https://sylabs.io/singularity/)/load a singularity module if using a shared HPC
- Clone this repo:
```bash
git clone --recursive -j8 https://github.com/blue-moon22/Palidis.git
cd palidis
```
If you have already cloned this repo with `git clone https://github.com/blue-moon22/Palidis.git`, you also need to get the submodules:
```bash
cd palidis
git submodule update --init --recursive
```

## Usage

### Without COBS Index Search
```bash
nextflow palidis.nf --manifest <manifest_file> --batch_name <batch_name> -c configs/conf/<name_of_config>.config
```

### With COBS Index Search
Download a COBS index database of all genomes. This is a very large file of just under 1 Terabyte so you will need a good internet connection and storage.
```
wget http://ftp.ebi.ac.uk/pub/databases/ENA2018-bacteria-661k/661k.cobs_compact
```

Run command with `--cobs_index` option
```bash
nextflow palidis.nf --manifest <manifest_file> --batch_name <batch_name> --cobs_index 661k.cobs_compact -c configs/conf/<name_of_config>.config
```

### Mandatory arguments
#### `<batch_name>`

`<batch_name>` must be the directory that the output is stored in.

#### `<manifest_file>`

A tab-delimited manifest must be specified for `--manifest` containing the absolute paths with headers `lane_id`, `read1`, `read2`, `sample_id` and `contigs_path`, e.g. this manifest contains three samples (the first having two lanes and the other two having one lane):

lane_id | read1 | read2 | sample_id | contigs_path
:---: | :---: | :---: | :---: | :---:
lane1 | /path/to/file/lane1_1.fq.gz | /path/to/file/lane1_2.fq.gz | my_sample | /path/to/file/contigs.fasta
lane2 | /path/to/file/lane2_1.fq.gz | /path/to/file/lane2_2.fq.gz | my_sample1 | /path/to/file/my_sample1_contigs.fasta
lane3 | /path/to/file/lane3_1.fq.gz | /path/to/file/lane3_2.fq.gz | my_sample2 | /path/to/file/my_sample2_contigs.fasta
lane4 | /path/to/file/lane4_1.fq.gz | /path/to/file/lane4_2.fq.gz | my_sample3 | /path/to/file/my_sample3_contigs.fasta

#### `<name_of_config>`

This represents the institution or HPC name. You can find your institutional HPC's config in `configs/conf` (which is linked to the configs directory in [nf-core](https://github.com/nf-core). For example, running on Sanger's HPC: `-c configs/conf/sanger.conf`

### Optional arguments
```
  --min_itr_length    Minimum length of ITR. (Default: 25)
  --max_itr_length    Maximum length of ITR. (Default: 50)
  --kmer_length       k-mer length for maximal exact matching. (Default: 15)
  --min_is_len        Minimum length of insertion sequence. (Default: 500)
  --max_is_len        Maximum length of insertion sequence. (Default: 3000)
  --cd_hit_G          -G option for CD-HIT-EST. (Default: 0)
  --cd_hit_aL         -aL option for CD-HIT-EST. (Default: 0.0)
  --cd_hit_aS         -aS option for CD-HIT-EST. (Default: 0.9)
  --cd_hit_c          -c option for CD-HIT-EST. (Default: 0.9)
  --e_value           -evalue option for BLASTn against ISfinder database. (Default: 1e-50)
  --cobs_index        Location of COBS index file for optional COBS index search of predicted IS origin. (Default: "")
  --cobs_threshold    K-mer threshold for identifying sequences in COBS index. (Default: 1)
  -resume             Resume the pipeline
```

## Output
There are two output files stored in a directory specified with `--batch_name`:

**1. FASTA file of insertion sequences**

**2. Information for each insertions sequence**

e.g. (includes information from optional COB index search)

IS_name | sample_id | contig | itr1_start_position | itr1_end_position | itr2_start_position | itr2_end_position | itr_cluster | ISfinder_name | ISfinder_origin | predicted_IS_family | COB_index_biosample_id | COB_index_origin
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
IS_name1 | sample_id1 | contig_name1 | 29 | 53 | 1004 | 1028 | 12 | ISBvu4 | Bacteroides vulgatus | | SAMN00627906 | Bacteroides vulgatus CL09T03C04 |
IS_name2 | sample_id1 | contig_name2 | 23 | 53 | 2769 | 2832 | 65 | | | ISLre2 | | | |
