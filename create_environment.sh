#!/bin/bash

conda create -y -n makebigwigfiles python=2.7
source activate makebigwigfiles
conda install -y -c bioconda \
ucsc-bedgraphtobigwig=357 \
ucsc-bedsort=357 \
samtools=1.6 \
bedtools=2.26 \
pybedtools=0.7.10 \
pysam=0.15.1

# assumes setup.py is in this directory! #
python setup.py install