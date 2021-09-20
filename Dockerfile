#######################################################
# Dependencies docker image for the GBS Typer pipeline.
#######################################################

FROM ubuntu:20.10
ARG DEBIAN_FRONTEND=noninteractive

WORKDIR /opt

# Base system
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
      && ln -s /usr/bin/python3 /usr/bin/python \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

# Install pal-MEM
ARG PALMEM_VERSION=2.2.1
RUN git clone --branch v${PALMEM_VERSION} https://github.com/blue-moon22/pal-MEM.git \
  && cd pal-MEM \
  && rm -rf .git \
  && mkdir build \
  && cd build \
  && cmake .. \
  && make

# Install Python3 packages
RUN pip3 install Bio

# Install CD-HIT
ARG CDHIT_VERSION=4.8.1
RUN wget -q -O- https://github.com/weizhongli/cdhit/releases/download/V${CDHIT_VERSION}/cd-hit-v${CDHIT_VERSION}-2019-0228.tar.gz | tar -xzf - \
  && cd cd-hit-v${CDHIT_VERSION}-2019-0228 \
  && make

# Install bowtie2
ARG BOWTIE2_VERSION=2.4.2
RUN wget -q https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-linux-x86_64.zip \
  && unzip bowtie2-2.4.2-linux-x86_64.zip \
  && rm -f bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip

# Install samtools
ARG SAMTOOLS_VERSION=1.11
RUN wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && tar -xf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && rm -f samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && cd samtools-${SAMTOOLS_VERSION} \
  && ./configure --prefix=/opt/samtools-${SAMTOOLS_VERSION} \
  && make \
  && make install

# Install R packages
RUN Rscript -e "install.packages(pkgs = c('optparse', 'dplyr'))"

# Install seqtk
ARG SEQTK_VERSION=1.3
RUN git clone --branch v${SEQTK_VERSION} https://github.com/lh3/seqtk.git \
  && cd seqtk \
  && make

# Add paths
ENV PATH="/opt/seqtk:/opt/pal-MEM/build:/opt/cd-hit-v${CDHIT_VERSION}-2019-0228:/opt/bowtie2-${BOWTIE2_VERSION}-linux-x86_64:/opt/samtools-${SAMTOOLS_VERSION}:/opt/prodigal-${PRODIGAL_VERSION}:/opt/hmmer-${HMMER_VERSION}/bin:${PATH}"
