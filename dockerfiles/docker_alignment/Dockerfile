FROM ubuntu:16.04

MAINTAINER Chris Dean <cdean11@rams.colostate.edu>

RUN apt-get update && apt-get install -y \
        aufs-tools \
        automake \
        build-essential \
	cmake \
        wget \
        git \
        libbz2-dev \
        unzip \
        python \
	tar \
	zlib1g-dev \
	libncurses5-dev \
	libdatetime-perl \
	libxml-simple-perl \
	libdigest-md5-perl \
	bioperl \
	ncbi-blast+ \
	openjdk-8-jre \
        openjdk-8-jdk \
	tcsh \
	xvfb \
	python-numpy \
        python-qt4 \
        python-six \
        python-pip \
        python-lxml \
&& rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/BenLangmead/bowtie2.git && \
        cd bowtie2 && \
        make NO_TBB=1 && \
        make install

RUN git clone https://github.com/samtools/htslib.git && \
	cd htslib && \
	autoheader && \
	autoconf && \
	./configure && \
	make && \
	make install

RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
	tar -xvf samtools-1.3.1.tar.bz2 && \
	rm samtools-1.3.1.tar.bz2 && \
	cd samtools-1.3.1 && \
	make && \
	make prefix=/usr/local install

RUN wget https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2 && \
	tar -xvf bcftools-1.3.1.tar.bz2 && \
	rm bcftools-1.3.1.tar.bz2 && \
	cd bcftools-1.3.1 && \
	make && \
	make install

RUN git clone --recursive git://github.com/ekg/freebayes.git && \
	cd freebayes && \
	make && \
	make install

RUN git clone https://github.com/tseemann/prokka.git && \
	prokka/bin/prokka --setupdb

RUN wget -N ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz && \
        gunzip linux64.tbl2asn.gz && \
        chmod +x linux64.tbl2asn && \
        mv linux64.tbl2asn prokka/binaries/linux/tbl2asn

RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip && \
        unzip Trimmomatic-0.36.zip && \
        rm Trimmomatic-0.36.zip && \
        cd Trimmomatic-0.36 && \
        chmod u+x trimmomatic-0.36.jar

RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.6/jellyfish-2.2.6.tar.gz && \
        tar -xvf jellyfish-2.2.6.tar.gz && \
        cd jellyfish-2.2.6 && \
        ./configure && \
        make -j 4 && \
        mv bin/jellyfish /usr/bin

RUN git clone https://github.com/cdeanj/kSNP3.git && \
	mv kSNP3 /usr/local

RUN pip install --upgrade ete3

RUN git clone https://github.com/cdeanj/newickphylo.git && \
	chmod u+x newickphylo/phylo && \
	mv newickphylo/phylo /usr/local/bin && \
	rm -rf newickphylo 

ENV PATH="/prokka/bin:${PATH}"
