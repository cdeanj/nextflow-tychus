FROM ubuntu:15.10

MAINTAINER Chris Dean <cdean11@rams.colostate.edu>

RUN apt-get update && apt-get install -y \
        aufs-tools \
        automake \
        build-essential \
        wget \
        git \
        libbz2-dev \
        unzip \
        python \
	tar \
	zlib1g-dev \
	mummer \
	ncbi-blast+ \
	r-base \
	r-base-dev \
	libdatetime-perl \
	libxml-simple-perl \
	libdigest-md5-perl \
	bioperl \
	libsparsehash-dev \
	libboost1.58-dev \
	openjdk-8-jre \
	openjdk-8-jdk \
&& rm -rf /var/lib/apt/lists/*

RUN wget http://kmergenie.bx.psu.edu/kmergenie-1.7016.tar.gz && \
	tar -xvf kmergenie-1.7016.tar.gz && \
	rm kmergenie-1.7016.tar.gz && \
	cd kmergenie-1.7016 && \
	make && \
	make install

RUN git clone https://github.com/cdeanj/idba.git && \
	cd idba && \
	./build.sh && \
	mv bin/idba_ud /usr/local/bin && \
	mv bin/fq2fa /usr/local/bin

RUN git clone https://github.com/dzerbino/velvet.git && \
	cd velvet && \
	make MAXKMERLENGTH=100 && \
	mv velvetg velveth /usr/local/bin

RUN wget http://spades.bioinf.spbau.ru/release3.9.0/SPAdes-3.9.0-Linux.tar.gz && \
    	tar -xzf SPAdes-3.9.0-Linux.tar.gz && \
    	mv SPAdes-3.9.0-Linux/bin/* /usr/local/bin && \
	mv SPAdes-3.9.0-Linux/share/* /usr/local/share && \
	rm ../../SPAdes-3.9.0-Linux.tar.gz

RUN wget http://sb.nhri.org.tw/CISA/upload/en/2012/9/CISA_20120904-14032356.tar && \
	tar -xvf CISA_20120904-14032356.tar && \
	rm CISA_20120904-14032356.tar && \
	chmod 755 -R CISA1.2 && \
	cd CISA1.2

RUN git clone https://github.com/tseemann/prokka.git && \
	prokka/bin/prokka --setupdb

RUN git clone https://github.com/bcgsc/abyss.git && \
	cd abyss && \
	./autogen.sh && \
	./configure --enable-maxk=160 && \
	make && \
	make install

RUN wget -N ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz && \
	gunzip linux64.tbl2asn.gz && \
	chmod +x linux64.tbl2asn && \
	mv linux64.tbl2asn prokka/binaries/linux/tbl2asn

ENV PATH="/CISA1.2:${PATH}"
ENV PATH="/prokka/bin:${PATH}"
