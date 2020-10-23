# PaliDIS v1.0.0 - Palindromic Detection of Insertion Sequences
Identifies insertion sequences from paired-end short-read metagenomic reads using Nextflow.

## DESCRIPTION
PaliDIS is a Nextflow pipeline consisting of five steps:
1. Assembling metagenomic reads to make contigs using [SPAdes](https://github.com/ablab/spades)
2. Clustering reads to make representative sequences using [MMSeq2](https://github.com/soedinglab/MMseqs2)
3. Using [pal-MEM](https://github.com/blue-moon22/pal-MEM) to identify reads with inverted repeats
4. Identifying contigs associated with inverted repeats mapping reads from 3) using Bowtie2
5. Identifying inverted terminal repeats of insertion sequences by transposase annotation using HMMER3

NOTE: Reads need to be quality controlled and filtered (i.e. removing human DNA) before using PaliDIS

### INSTALLATION INSTRUCTIONS
PaliDIS relies on Nextflow and Docker images of pipeline software.
Download:
1. [Nextflow](https://www.nextflow.io/).
2. [Docker](https://www.docker.com/).
3. Docker images:
```
docker pull staphb/spades:3.14.0 && docker pull soedinglab/mmseqs2:latest && docker pull bluemoon222/pal-mem:latest && docker pull python && docker pull comics/bowtie2:2.3.4.1 && docker pull pegi3s/samtools_bcftools:1.9 && docker pull ubuntu:xenial
```
4. PaliDIS:
```
git clone https://github.com/blue-moon22/PaliDIS.git
```

### USAGE
```
nextflow run palidis.nf --reads '<sample>_{1,2}.fastq.gz'
```

As an example/to test PaliDIS, within the directory run:
```
nextflow run palidis.nf --reads 'data/ERR589346_head_{1,2}.fastq.gz'
```
This should take approx. 5 minutes. (Make sure to include the quotes '')
You should get two empty files ERR589346_head_contig_is_sites.fasta and ERR589346_head_contig_is_sites.tab

If you want to run all samples in the directory:
```
nextflow run palidis.nf --reads '*_{1,2}.fastq.gz'
```

If you are running this on large files (above 2 GB), it is recommended to use a cluster or HPC. If you use a particular batch software on your cluster, you can specify the type (e.g. SGE/SLURM), number of cpus and other parameters in the nextflow.config file. This [page](https://www.nextflow.io/docs/latest/config.html) gives you an overview of the configuration file and this [page](https://www.nextflow.io/docs/latest/executor.html?highlight=sge) specifies what to include for particular batch software.

### OPTIONS
PaliDIS comes with two more options: --kmer_length and --mem_length, which represent the options -k and -l in [pal-MEM](https://github.com/blue-moon22/pal-MEM). These select the k-mer lengths to be queried and the minimum maximal exact match (MEM) length of an inverted terminal repeat.
