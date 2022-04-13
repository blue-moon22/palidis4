<img src="img/logo.png" alt="logo" width="400"/>

# **PaliDIS v2.8.0** - **Pali**ndromic **D**etection of **I**nsertion **S**equences

PaliDIS is a Nextflow pipeline that quickly discovers novel insertion sequences.

The tool is based upon identifying inverted terminal repeats (ITRs) (figure below) using paired-end, short-read mixed microbial genomic data (e.g. metagenomes).

<img src="img/insertion_sequence.png" alt="insertion sequence" width="400"/>

## Installation
- Install [Nextflow](https://www.nextflow.io/)
- Install [Docker](https://www.docker.com/) if using own machine or install [Singularity](https://sylabs.io/singularity/)/load a singularity module if using a shared HPC
- Git clone this repo
```bash
git clone https://github.com/blue-moon22/Palidis.git
cd Palidis
```

## Pipeline summary
**Steps:**
1. Pre-process FASTQ.GZ reads [`convertToFasta`]
2. Efficient maximal exact matching to get repeat sequences using [pal-MEM](https://github.com/blue-moon22/pal-MEM) [`palmem`]
3. Map reads against assemblies using Bowtie2 [`filterContigs` `buildDB` `mapreads`]
4. Get candidate ITRs by distance filters [`getCandidateITRs`]
5. Cluster candidate ITRs using CD-HIT-EST [`clusterReads`]
6. Get putative ITRs by cluster concordance and output Insertion Sequences [`getITRs`]

## Usage
```bash
nextflow palidis.nf --manifest <manifest_file> --batch_name <batch_name> -profile <executor>
```
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

#### `<executor>`

If you are running this on your own machine, then you should specify `-profile standard`.

If you are running this on an HPC, you will need to specify `-profile <executor>` in the command. Currently, the pipeline only supports `farm` (using `-profile farm`) and `rosalind` (using `-profile rosalind`). If you use `rosalind`, the default partition is `brc`. If you want to use a different partition, include option `--partition <name>`.

It is possible to add another profile to the [nextflow config](https://www.nextflow.io/docs/latest/config.html) to make this pipeline compatible with other HPC executors. If you do so, you are welcome to fork this repo and make a pull request to include your new profile for others to use. You may be able to find a basic config for your HPC [here](https://github.com/nf-core/configs/tree/master/conf).

### More Options
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
  -resume             Resume the pipeline
```

## Output
There are two output files stored in a directory specified with `--batch_name`:

**1. FASTA file of insertion sequences**

**2. Information for each insertions sequence** e.g.

IS_name | sample_id | contig | itr1_start_position | itr1_end_position | itr2_start_position | itr2_end_position | itr_cluster
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
IS_name1 | sample_id1 | contig_name1 | 29 | 53 | 1004 | 1028 | 12
IS_name2 | sample_id1 | contig_name2 | 23 | 53 | 2769 | 2832 | 65
