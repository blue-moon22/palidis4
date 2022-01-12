<img src="img/logo.png" alt="logo" width="400"/>

# **PaliDIS v2.6** - **Pali**ndromic **D**etection of **I**nsertion **S**equences

PaliDIS is a Nextflow pipeline that predicts insertion sequence annotations of paired-end, short-read metagenomic data.

The tool is based upon identifying inverted terminal repeats (ITRs) (figure below) of insertion sequences using an efficient maximal exact matching algorithm between reads of large metagenomic datasets, which are then resolved on metagenomic contigs.

<img src="img/insertion_sequence.png" alt="insertion sequence" width="400"/>

## Contents
- [ Installation ](#installation)
- [ Pipeline summary ](#summary)
- [ Workflow 1: Get IS annotations ](#workflow1)
    - [ Usage ](#usage1)
    - [ Output ](#output1)
    - [ Options ](#options1)
- [ Workflow 2: Create catalog ](#workflow2)
    - [ Usage ](#usage2)
    - [ Output ](#output2)
    - [ Options ](#options2)

<a name="installation"></a>
## Installation
- Install [Nextflow](https://www.nextflow.io/)
- Install [Docker](https://www.docker.com/) ([Singularity](https://sylabs.io/singularity/) if using a HPC)
- Git clone this repo
```bash
git clone https://github.com/blue-moon22/Palidis.git
cd Palidis
```

<a name="summary"></a>
## Pipeline summary
There are two workflows: **1) Get IS annotations** that gets the positional information of ITRs on metagenomic assemblies and **2) Create catalog** that creates a catalog of IS annotations from multiple samples and clusters ITRs into distinct ITR clusters.

**Get IS annotations steps:**
1. Index and label FASTQ.GZ reads
2. Efficient maximal exact matching to get inverted repeats using [pal-MEM](https://github.com/blue-moon22/pal-MEM)
3. Map reads against contigs using Bowtie2
4. Get candidate ITRs by paired read concordance
5. Cluster candidate ITRs using CD-HIT-EST
6. Get putative ITRs by cluster concordance and pairs are found between 700 and 3000 bp (typical IS length) of each other
7. Collect Insertion Sequence annotations on contigs

**Create catalog:**
1. Clusters all ITRs from samples using CD-HIT-EST
2. Creates a catalog of Insertion Sequence annotations for all samples with assigned ITR Clusters

<a name="workflow1"></a>
## Workflow 1: Get IS annotations
<a name="usage1"></a>
### Usage

Workflow `get_IS_annotations` generates Insertion Sequence annotations for each sample.
```bash
nextflow palidis.nf --get_IS_annotations --manifest <manifest_file> --batch_name <batch_name> --min_itr_length <min_itr_length> --kmer_length <kmer_length> --resume -profile <executor>
```
The output is stored in a directory in the current run directory specified with `--batch_name`.

A tab-delimited manifest must be specified for `--manifest` containing the absolute paths with headers `lane_id`, `read1`, `read2`, `sample_id` and `contigs_path`, e.g. this manifest contains three samples (the first having two lanes and the other two having one lane):

lane_id | read1 | read2 | sample_id | contigs_path
:---: | :---: | :---: | :---: | :---:
lane1 | /path/to/file/lane1_1.fq.gz | /path/to/file/lane1_2.fq.gz | my_sample | /path/to/file/contigs.fasta
lane2 | /path/to/file/lane2_1.fq.gz | /path/to/file/lane2_2.fq.gz | my_sample1 | /path/to/file/my_sample1_contigs.fasta
lane3 | /path/to/file/lane3_1.fq.gz | /path/to/file/lane3_2.fq.gz | my_sample2 | /path/to/file/my_sample2_contigs.fasta
lane4 | /path/to/file/lane4_1.fq.gz | /path/to/file/lane4_2.fq.gz | my_sample3 | /path/to/file/my_sample3_contigs.fasta

If you are running this on an HPC, you will need to specify `-profile <executor>` in the command. Currently, the pipeline only supports `LSF` (using `-profile lsf`) and `rosalind` (using `-profile rosalind`). If you use `rosalind`, the default partition is `brc`. If you want to use a different partition, include option `--partition <name>`. It is possible to add another profile to the [nextflow config](https://www.nextflow.io/docs/latest/config.html) to make this pipeline compatible with other HPC executors. If you do so, you are welcome to fork this repo and make a pull request to include your new profile for others to use.

<a name="output1"></a>
### Output
Three files for each sample are generated in a directory specified by `batch_name`:

**1. A non-redundant catalogue of ITRs**: `<sample_id>_ITRs.fasta`

**2. Insertion sequence annotations**: `<sample_id>_insertion_sequence_annotations.tab`
The annotation file is in a tab-delimited format consisting of the sample_id, contig name, start and end positions of the first ITR, start and end positions of the second ITR and the cluster(s) they belong to (`itr_clusters`), e.g.:

sample_id | contig | itr1_start_position | itr1_end_position | itr2_start_position | itr2_end_position | itr_cluster
:---: | :---: | :---: | :---: | :---: | :---: | :---:
sample_id1 | contig_name1 | 29 | 43 | NA | NA | 1217817
sample_id1 | contig_name2 | 23 | 43 | 2769 | 2822 | 656079

Although two flanking ITRs may be found, it is possible that positions could not be predicted (represented by `NA`). (This happens when a read maps to a contig, but its paired read containing the ITR does not.)

**3. ITR clusters and their reads of origin:** `<sample>_reads_itr_clusters.txt`
This file is in a tab-delimited format containing the names of the ITR clusters and the reads that the ITRs of those ITR clusters originate from.

<a name="options1"></a>
### Options
```
  min_itr_length      Minimum length of ITR. (Default: 14)
  kmer_length         k-mer length for maximal exact matching. (Default: 10)
  split               Split reference in pal-MEM by this number. (Default: 20)
  cd_hit_G            -G option for CD-HIT-EST. (Default: 0)
  cd_hit_aL           -aL option for CD-HIT-EST. (Default: 0.0)
  cd_hit_aS           -aS option for CD-HIT-EST. (Default: 1.0)
```

<a name="workflow2"></a>
## Workflow 2: Create catalog
<a name="usage2"></a>
### Usage
Workflow `create_catalog` generates a catalog of Insertion Sequence annotations for all samples.
```bash
nextflow run palidis.nf --create_catalog --batch_name <batch_name> --min_itr_length <min_itr_length> -profile <executor>
```
The `batch_name` is the directory that contains outputs from the previous workflow for all samples.

<a name="output2"></a>
### Output
One tab-delimited catalog is created called `<batch_name>_insertion_sequence_annotations_catalog.tab` that contains the sample_id, contig name, start and end positions of the first ITR, start and end positions of the second ITR and the newly assigned cluster(s) they belong to in the catalog (`itr_cluster_catalog`), e.g.:

sample_id | contig | itr1_start_position | itr1_end_position | itr2_start_position | itr2_end_position | itr_cluster_catalog
:---: | :---: | :---: | :---: | :---: | :---: | :---:
sample_id1 | contig_name1 | 29 | 43 | NA | NA | 101
sample_id1 | contig_name2 | 23 | 43 | 2769 | 2822 | 102
sample_id2 | contig_name3 | 25 | 55 | 5738 | 5768 | 101

<a name="options2"></a>
### Options
```
  min_itr_length      Minimum length of ITR. (Default: 14)
  cd_hit_G            -G option for CD-HIT-EST. (Default: 0)
  cd_hit_aL           -aL option for CD-HIT-EST. (Default: 0.0)
  cd_hit_aS           -aS option for CD-HIT-EST. (Default: 1.0)
```
