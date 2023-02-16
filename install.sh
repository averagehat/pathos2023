#!/bin/bash
set -e
MINICONDA=$1
# wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh 
# bash Miniconda2-latest-Linux-x86_64.sh -b -p $MINICONDA
BIN=$MINICONDA/bin/

# == PriceSeqFilter
PRICESOURCE=PriceSource140408.tar.gz 
wget http://derisilab.ucsf.edu/software/price/$PRICESOURCE
tar -xvf $PRICESOURCE
cd PriceSource140408 && make && ln -s $PWD/PriceSeqFilter $BIN/PriceSeqFilter && cd ..

# == GSNAPL # bin/gmap
conda install -y gmap

# == RepeatMasker  # bin/RepeatMasker
# conda install repeatmasker

# == STAR # bin/STAR
conda install -y star

# == cd-hit-dup
# conda install cd-hit  # doesn'twork

# == cd-hit-dup # cd-hit-dup
git clone https://github.com/weizhongli/cdhit
cd cdhit/cd-hit-auxtools/ && make && ln -s $PWD/cdhit/cd-hit-auxtools/cd-hit-dup  $BIN/cd-hit-dup

# == bowtie2 (already installed)
# the below version is screwy I think.
# conda install -y bowtie2

#  ??? 
# I think this samtools verison works right.
conda install samtools
conda install ete2

# == RapSearch2 # RAPSearch2/bin/rapsearch
# git clone https://github.com/zhaoyanswill/RAPSearch2
# cd RAPSearch2 && ./install && cd .. && ln -s $PWD/RAPSearch2/bin/rapsearch $BIN/rapsearch


