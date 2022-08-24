

[![GitHub release (latest by date)](https://img.shields.io/github/v/release/blue-moon22/palidis)](https://github.com/blue-moon22/palidis/releases)
![GitHub Workflow Status](https://img.shields.io/github/workflow/status/blue-moon22/palidis/test)

<img src="img/logo.png" alt="logo" width="400"/>

# **PaliDIS** - **Pali**ndromic **D**etection of **I**nsertion **S**equences

## Contents
- [ Introduction ](#introduction)
- [ Pipeline description ](#description)
- [ Installation ](#installation)
- [ Usage ](#usage)
- [ Output ](#output)

<a name="introduction"></a>
## Introduction

PaliDIS is a Nextflow pipeline that quickly discovers novel insertion sequences.

The tool is based upon identifying inverted terminal repeats (ITRs) (figure below) using paired-end, short-read metagenomic data.

For each sample, the pipeline produces two output files: **1. FASTA file of insertion sequences** and **2. Information for each insertions sequence**

<img src="img/insertion_sequence.png" alt="insertion sequence" width="400"/>

<a name="description"></a>
## Pipeline description
**Steps:**
1. Pre-process FASTQ.GZ reads [`convertToFasta`]
2. Efficient maximal exact matching to get repeat sequences using [pal-MEM](https://github.com/blue-moon22/pal-MEM) [`palmem`]
3. Map reads against assemblies using Bowtie2 [`filterContigs` `buildDB` `mapreads`]
4. Get candidate ITRs by distance filters [`getCandidateITRs`]
5. Cluster candidate ITRs using CD-HIT-EST [`clusterReads`]
6. Get putative ITRs by cluster concordance and output Insertion Sequences [`getITRs`]
7. Find transposase [`runProdigal`] [`installInterproscan` `runInterproscan`]
8. _Optional:_ Search against a COB index to predict IS origin [`searchCOBSIndex`]
9. Combine insertion sequence information [`getISInfoWithCOBS` `getISInfoWithoutCOBS`]

<a name="installation"></a>
## Installation
- Install [Nextflow](https://www.nextflow.io/)
- Install [Docker](https://www.docker.com/) if using own machine or install [Singularity](https://sylabs.io/singularity/)/load a singularity module if using a shared HPC
- Clone this repo:
```bash
git clone --recursive -j8 https://github.com/blue-moon22/palidis.git
cd palidis
```
_Note: You may be warned to first call `git config --global --add safe.directory`._

    If you have already cloned this repo with `git clone https://github.com/blue-moon22/palidis.git`, you also need to get the submodules:
```bash
cd palidis
git submodule update --init --recursive
```

<a name="usage"></a>
## Usage

### Without COBS Index Search
```bash
nextflow palidis.nf --manifest <manifest_file> --batch_name <batch_name> -c configs/conf/<name_of_config>.config
```
**If you are running this on an LSF scheduler, also include `--lsf true`.**

### With COBS Index Search
Download a COBS index database of all genomes. This is a very large file of just under 1 Terabyte so you will need a good internet connection and storage.
```
wget http://ftp.ebi.ac.uk/pub/databases/ENA2018-bacteria-661k/661k.cobs_compact
```

Run command with `--cobs_index` option
```bash
nextflow palidis.nf --manifest <manifest_file> --batch_name <batch_name> --cobs_index 661k.cobs_compact -c configs/conf/<name_of_config>.config
```
**If you are running this on an LSF schedular, also include `--lsf true`.**

### Mandatory arguments
#### `<batch_name>`

`<batch_name>` must be the directory that the output is stored in.

#### `<manifest_file>`

A tab-delimited manifest must be specified for `--manifest` containing the absolute paths with headers `lane_id`, `read1`, `read2`, `sample_id` and `contigs_path`, e.g. this manifest contains three samples (the first having two lanes and the other two having one lane):

lane_id | read1 | read2 | sample_id | contigs_path
:---: | :---: | :---: | :---: | :---:
lane1 | /path/to/file/lane1_1.fq.gz | /path/to/file/lane1_2.fq.gz | my_sample1 | /path/to/file/contigs.fasta
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
  --cobs_index        Location of COBS index file for optional COBS index search of predicted IS origin. (Default: "")
  --cobs_threshold    K-mer threshold for identifying sequences in COBS index. (Default: 1)
  -resume             Resume the pipeline
```

### Testing
If you would like to test whether this pipeline produces the expected output on your system, run this command. If successful, it should print `Test passed. All outputs expected.`.
```
./tests/regression_tests.sh
```

<a name="output"></a>
## Output
There are two output files stored in a directory specified with `--batch_name`:

**1. FASTA file of insertion sequences**

**2. Information for each insertions sequence**

e.g. (includes information from optional COB index search)

IS_name | sample_id | contig | itr1_start_position | itr1_end_position | itr2_start_position | itr2_end_position | itr_cluster | COBS_index_biosample_id | COBS_index_origin
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
IS_name1 | sample_id1 | contig_name1 | 29 | 53 | 1004 | 1028 | 12 | SAMN00627906 | Bacteroides vulgatus CL09T03C04 |
IS_name2 | sample_id1 | contig_name2 | 23 | 53 | 2769 | 2832 | 65 | | |

### Interpretation
Header | Description
:--- | :---
**IS_name** | Name assigned by PaliDIS which contains the ITR cluster (see below) and length e.g. `IS_cluster_0_length_1072`
**sample_id** | Sample ID that was given in manifest
**contig** | Name of the contig that was given by the header in the contig file provided by the manifest
**itr1_start_position** | The position of the first nucleotide of the left-hand Inverted Terminal Repeat (ITR) sequence
**itr1_end_position** | The position of the last nucleotide of the left-hand ITR sequence
**itr2_start_position** | The position of the first nucleotide of the right-hand ITR sequence
**itr2_end_position** | The position of the last nucleotide of the right-hand ITR sequence
**itr_cluster** | The ITR cluster that was assigned to both ITRs (in Step 5)
**COBS_index_biosample_id** | The NCBI Biosample ID of a sequenced sample in the COBS index
**COBS_index_origin** | The taxonomy of the sequenced sample in the COBS index
