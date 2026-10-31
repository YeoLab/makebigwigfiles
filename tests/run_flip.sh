#!/bin/bash

# making sure the dont_flip flag is actually doing that

python /tscc/projects/ps-yeolab3/bay001/codebase/makebigwigfiles/makebigwigfiles/make_bigwig_files.py \
--bam inputs/eCLIP/flipped/204_01_RBFOX2.merged.r2.bam \
--genome inputs/hg19.chrom.sizes \
--bw_pos outputs/obtained/flipped/204_01_RBFOX2.merged.r2.norm.pos.bw \
--bw_neg outputs/obtained/flipped/204_01_RBFOX2.merged.r2.norm.neg.bw \
--direction f
