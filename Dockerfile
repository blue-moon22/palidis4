#######################################
# Dependencies docker image for PaliDIS
#######################################

# syntax=docker/dockerfile:1
FROM ubuntu:20.04 AS builder
ARG DEBIAN_FRONTEND=noninteractive
WORKDIR /opt
RUN apt-get update -y -qq && apt-get install -y -qq \
        wget \
        python3 \
        python3-pip \
        gzip \
        git \
        build-essential \
        autoconf \
        automake \
        cmake \
        libboost-all-dev \
        unzip \
        perl \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-gnutls-dev \
        libssl-dev \
        libncurses5-dev \
        r-base \
        openjdk-11-jre \
      && ln -s /usr/bin/python3 /usr/bin/python \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*
ENV JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
ENV CLASSPATH=/usr/lib/jvm/java-11-openjdk-amd64/bin

# Install pal-MEM
FROM builder AS build1
WORKDIR /opt
ARG PALMEM_VERSION=2.3.4
RUN git clone --branch v${PALMEM_VERSION} https://github.com/blue-moon22/pal-MEM.git \
  && cd pal-MEM \
  && rm -rf .git \
  && mkdir build \
  && cd build \
  && cmake .. \
  && make \
  && mv pal-mem /opt

# Install CD-HIT
FROM builder AS build2
WORKDIR /opt
ARG CDHIT_VERSION=4.8.1
RUN wget -q -O- https://github.com/weizhongli/cdhit/releases/download/V${CDHIT_VERSION}/cd-hit-v${CDHIT_VERSION}-2019-0228.tar.gz | tar -xzf - \
  && cd cd-hit-v${CDHIT_VERSION}-2019-0228 \
  && make \
  && mv cd-hit-est /opt

# Install bowtie2
FROM builder AS build3
WORKDIR /opt
ARG BOWTIE2_VERSION=2.4.2
RUN wget -q https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-linux-x86_64.zip \
  && unzip bowtie2-2.4.2-linux-x86_64.zip \
  && rm -f bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip \
  && mv bowtie2-${BOWTIE2_VERSION}-linux-x86_64/bowtie* /opt

# Install samtools
FROM builder AS build4
WORKDIR /opt
ARG SAMTOOLS_VERSION=1.11
RUN wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && tar -xf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && rm -f samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && cd samtools-${SAMTOOLS_VERSION} \
  && ./configure --prefix=/opt/samtools-${SAMTOOLS_VERSION} \
  && make \
  && make install \
  && mv samtools /opt

# Install BLAST
FROM builder AS build5
WORKDIR /opt
ARG BLAST_VERSION=2.12.0
RUN wget -q https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz \
 && tar -xf ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz \
 && mv ncbi-blast-${BLAST_VERSION}+/bin/blastn /opt \
 && rm ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz

FROM builder
COPY --from=build1 /opt/pal-mem /bin
COPY --from=build2 /opt/cd-hit-est /bin
COPY --from=build3 /opt/bowtie* /bin
COPY --from=build4 /opt/samtools /bin
COPY --from=build5 /opt/blastn /bin
