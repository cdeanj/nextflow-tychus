FROM ubuntu:14.04

MAINTAINER Chris Dean <cdean11@rams.colostate.edu>

RUN apt-get update && apt-get install -y \
        aufs-tools \
        automake \
        build-essential \
        cmake \
        zlib1g-dev \
        wget \
        git \
        libbz2-dev \
        openjdk-7-jre-headless \
        unzip \
        python \
&& rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/BenLangmead/bowtie2.git && \
        cd bowtie2 && \
        make && \
        cp bowtie2* /usr/local/bin

RUN git clone https://github.com/cdeanj/coverage_sampler.git && \
	cd coverage_sampler && \
	make && \
	cp csa /usr/local/bin
