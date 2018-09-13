#!/usr/bin/env python

# __author__ = 'gpratt'

import os
import subprocess
from subprocess import call
import argparse
import pysam

def flip_strands(fn, out_fn):
    """
    Creates a bam file whose strands are flipped from the original file.
    """
    with pysam.AlignmentFile(fn, "rb") as aligned:
        with pysam.AlignmentFile(out_fn, "wb", template=aligned) as aligned_out:
            for read in aligned:
                read.is_reverse = not read.is_reverse
                aligned_out.write(read)

def genome_coverage_bed(in_bam=None, in_bed=None, out_bed_graph=None, genome=None, strand=None, split=True,
                        dont_flip=False, out_flipped_bam=None, five=False):
    with open(out_bed_graph, 'w') as out_bed_graph:
        if in_bam is not None and in_bed is not None:
            raise Exception("can't pass both bam and bed file to this function")

        if dont_flip and in_bam:
            priming_call = "genomeCoverageBed -ibam {}".format(in_bam)
        elif in_bam is not None:
            flip_strands(in_bam, out_flipped_bam)
            check_for_index(out_flipped_bam)
            priming_call = "genomeCoverageBed -ibam {}".format(out_flipped_bam)
            # priming_call = "samtools view -h " + in_bam + " | awk 'BEGIN {OFS=\"\\t\"} {if(!!and($2,0x0080)) {if(!!and($2, 0x0004)) {$2 = $2 - 16} else {$2 = $2 + 16}}; print $0}' | samtools view -bS - | genomeCoverageBed -ibam stdin "
            #priming_call = "samtools view -h " + in_bam + " | samtools view -bS - | genomeCoverageBed -ibam stdin "
            #priming_call = "genomeCoverageBed -ibam " + in_bam

        if in_bed is not None:
            priming_call = "genomeCoverageBed -i {}".format(in_bed)

        priming_call += " -bg "
        if strand:
            priming_call += " -strand {} ".format(strand)

        if split:
            priming_call += " -split "

        priming_call += " -g {} ".format(genome)


        print("genome_coverage_bed priming_call:", priming_call, out_bed_graph)
        subprocess.check_call(priming_call, shell=True, stdout=out_bed_graph)


def normalize_bed_graph(in_bed_graph, in_bam, out_bed_graph):
    with open(out_bed_graph, 'w') as out_bed_graph:
        priming_call = "normalize_bedGraph.py "
        priming_call += " --bg {} ".format(in_bed_graph)
        priming_call += " --bam {}".format(in_bam)
        print("normalize_bed_graph priming_call:", priming_call, out_bed_graph)
        subprocess.call(priming_call, shell=True, stdout=out_bed_graph)


def bed_graph_to_big_wig(in_bed_graph, genome, out_big_wig):
    # BUG: getting error: 204_01_RBFOX2.merged.r2.norm.pos.bg is not case-sensitive sorted at line 1118962.
    #      Please use "sort -k1,1 -k2,2n" with LC_COLLATE=C,  or bedSort and try again.
    # FIX as seen here : http://seqanswers.com/forums/showthread.php?t=63932
    # also see: https://github.com/daler/pybedtools/issues/178

    # TODO temporarily going back to no sorting , to work around error message:
    # TODO needLargeMem: trying to allocate 0 bytes (limit: 100000000000)
    priming_call = "bedGraphSortToBigWig {} {} {}".format(in_bed_graph, genome, out_big_wig)
    print("bed_graph_to_big_wig priming_call:", priming_call)
    with open(os.devnull, 'w') as fnull:
        subprocess.call(priming_call, shell=True, stdout=fnull)


def neg_bed_graph(in_bed_graph, out_bed_graph):
    priming_call = "negBedGraph.py "
    priming_call += " --bg {}".format(in_bed_graph)
    print("neg_bed_graph priming_call:", priming_call, out_bed_graph)
    with open(out_bed_graph, 'w') as out_bed_graph:
        subprocess.call(priming_call, shell=True, stdout=out_bed_graph)


def check_for_index(bamfile):
    """

    Checks to make sure a BAM file has an index, if the index does not exist it is created

    Usage undefined if file does not exist (check is made earlier in program)
    bamfile - a path to a bam file

    """

    if not os.path.exists(bamfile):
        raise NameError("file %s does not exist" % (bamfile))

    if os.path.exists(bamfile + ".bai"):
        return

    if not bamfile.endswith(".bam"):
        raise NameError("file %s not of correct type" % (bamfile))
    else:
        process = call(["samtools", "index", str(bamfile)])

        if process == -11:
            raise NameError("file %s not of correct type" % (bamfile))

def main():
    parser = argparse.ArgumentParser(description="Makes Pretty bed Graph Files!")
    parser.add_argument("--bam", help="bam file to make bedgraphs from", required=True)
    parser.add_argument("--genome", help="chromosome sizes because some things need it", required=True)
    parser.add_argument("--dont_flip", help="by default assumes trueseq reversed strand, this disables that assumption (use flip for ",
                        action="store_true", default=False)    #False)
    parser.add_argument("--bw_pos", help="positive bw file name", required=False)
    parser.add_argument("--bw_neg", help="negative bw file name", required=False)
    parser.add_argument("--bw", help="(if unsplit) bw file name", required=False)
    parser.add_argument("--five_prime", help="plot only 5' end of fragment", action="store_true", default=False)
    parser.add_argument("--no_strand", help="don't split the bigwig file into positive and negative strands", required=False, default=False, action='store_true')

    args = parser.parse_args()
    bamFile = args.bam
    genome = args.genome
    no_strand = args.no_strand
    dont_flip = args.dont_flip
    five = args.five_prime

    flipped_bam = bamFile.replace(".bam", ".flip.bam")

    check_for_index(bamFile)

    if not no_strand:
        bedGraphFilePos = bamFile.replace(".bam", ".pos.bg")
        bedGraphFilePosNorm = bedGraphFilePos.replace(".pos.bg", ".norm.pos.bg")

        bedGraphFileNeg = bamFile.replace(".bam", ".neg.bg")
        bedGraphFileNegNorm = bedGraphFileNeg.replace(".neg.bg", ".norm.neg.bg")
        bedGraphFileNegInverted = bedGraphFileNegNorm.replace(".bg", ".t.bg")

        genome_coverage_bed(in_bam=bamFile, out_bed_graph=bedGraphFilePos, strand="+", genome=genome,
                            dont_flip=dont_flip, out_flipped_bam=flipped_bam, five=five)
        genome_coverage_bed(in_bam=bamFile, out_bed_graph=bedGraphFileNeg, strand="-", genome=genome,
                            dont_flip=dont_flip, out_flipped_bam=flipped_bam, five=five)
        if dont_flip:
            normalize_bed_graph(in_bed_graph=bedGraphFilePos, in_bam=bamFile, out_bed_graph=bedGraphFilePosNorm)
            normalize_bed_graph(in_bed_graph=bedGraphFileNeg, in_bam=bamFile, out_bed_graph=bedGraphFileNegNorm)
        else:
            normalize_bed_graph(in_bed_graph=bedGraphFilePos, in_bam=flipped_bam, out_bed_graph=bedGraphFilePosNorm)
            normalize_bed_graph(in_bed_graph=bedGraphFileNeg, in_bam=flipped_bam, out_bed_graph=bedGraphFileNegNorm)



        neg_bed_graph(in_bed_graph=bedGraphFileNegNorm, out_bed_graph=bedGraphFileNegInverted)

        bed_graph_to_big_wig(bedGraphFilePosNorm, genome, args.bw_pos)
        bed_graph_to_big_wig(bedGraphFileNegInverted, genome, args.bw_neg)
    else:
        bedGraphFile = bamFile.replace(".bam", ".bg")
        bedGraphFileNorm = bedGraphFile.replace(".bg", ".norm.bg")

        genome_coverage_bed(in_bam=bamFile, out_bed_graph=bedGraphFile, strand=False, genome=genome,
                            dont_flip=dont_flip, out_flipped_bam=flipped_bam)
        if dont_flip:
            normalize_bed_graph(in_bed_graph=bedGraphFile, in_bam=flipped_bam, out_bed_graph=bedGraphFileNorm)
        else:
            normalize_bed_graph(in_bed_graph=bedGraphFile, in_bam=bamFile, out_bed_graph=bedGraphFileNorm)
        bed_graph_to_big_wig(bedGraphFileNorm, genome, args.bw)

if __name__ == "__main__":
    main()
