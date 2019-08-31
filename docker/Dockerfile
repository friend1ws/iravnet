FROM ubuntu:18.04
MAINTAINER Yuichi Shiraishi <friend1ws@gmail.com> 

RUN apt-get update && apt-get install -y \
    git \
    wget \
    bzip2 \
    make \
    gcc \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    python \
    python3 \
    python3-pip


RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar jxvf samtools-1.9.tar.bz2 && \
    cd samtools-1.9/htslib-1.9 && ./configure && make && make install && \
    cd ../ && ./configure --without-curses && make && make install

RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
    tar jxvf bcftools-1.9.tar.bz2 && \
    cd bcftools-1.9 && ./configure && make && make install

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.27.0/bedtools-2.27.0.tar.gz && \
    tar -zxvf bedtools-2.27.0.tar.gz && \
    cd bedtools2 && make && make install

RUN wget https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/archive/v1.1.tar.gz && \
    tar zxvf v1.1.tar.gz && \
    cd Complete-Striped-Smith-Waterman-Library-1.1/src && \
    gcc -Wall -O3 -pipe -fPIC -shared -rdynamic -o libssw.so ssw.c ssw.h

ENV LD_LIBRARY_PATH /Complete-Striped-Smith-Waterman-Library-1.1/src:$LD_LIBRARY_PATH

RUN pip3 install --upgrade pip
RUN pip3 install --upgrade setuptools

ENV HTSLIB_LIBRARY_DIR /usr/local/lib
ENV HTSLIB_INCLUDE_DIR /usr/local/include

RUN pip3 install cython==0.29.11
RUN pip3 install annot_utils==0.3.0
RUN pip3 install intron_retention_utils==0.6.1

RUN wget https://github.com/pysam-developers/pysam/archive/v0.15.2.tar.gz && \
    tar -zxvf v0.15.2.tar.gz && \
    cd pysam-0.15.2 && python3 setup.py build install

ENV LD_LIBRARY_PATH /usr/local/lib:$LD_LIBRARY_PATH

RUN pwd

RUN git clone https://github.com/friend1ws/iravnet.git && \
    cd iravnet && pip3 install .

