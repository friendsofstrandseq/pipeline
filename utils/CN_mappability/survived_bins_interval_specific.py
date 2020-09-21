#!/usr/bin/python
# -*- coding: utf-8 -*-

# Hufsah 31.8.2020
# script to give a list of bins that survived all the filtraion i.e. the ones
# that would be used for normalization for each segment # changing it to
# take a track with 5 columns chr, start, end, correctly_mapped_reads,
# incorrctly_mapped_reads...this track is just filtered for mapq so further
# filtration for incorrect/correct reads >0.1 bins needs to be done here.
# Plus the usual rule to keep only the bins with mapping_reads >=75

import argparse
import sys
import pysam
from pathlib import Path
import glob
import pdb

from collections import defaultdict


def counts(input_bed, mapping_counts, output_file):
    dictionary = defaultdict(lambda : defaultdict(tuple))
    survived_bins = open(output_file, 'w')
    mapping_counts_file = open(mapping_counts, 'r')
    for lines in mapping_counts_file:
        line = lines.strip().split(' ')
        chrom = line[0]
        try:
            interval_start = int(line[1])
        except:
            pdb.set_trace()
        interval_end = int(line[2])
        reads_originated = 100
        reads_mapped_correctly = int(line[3])
        reads_mapped_incorrectly = int(line[4])
        dictionary[chrom][(interval_start, interval_end)] = \
            (reads_originated, reads_mapped_correctly,
             reads_mapped_incorrectly)
    segment = 0
    with open(input_bed, 'r') as bed_file:

        # next(bed_file) #skipping the header in bed file

        for line in bed_file:
            segment += 1  # next inversion
            line_r = line.strip().split('\t')
            chromosome = line_r[0]
            sub_dictionary = dictionary[chromosome]
            seg_start = int(line_r[1])
            seg_end = int(line_r[2])
            seg_start_bin = seg_start
            for m in sub_dictionary:
                correctly_mapped_count = sub_dictionary[m][1]
                incorrectly_mapped_count = sub_dictionary[m][2]

                # if segment start is towards the right of this bin_end, keep moving

                if seg_start_bin >= m[1]:
                    continue
                else:

                    # if segment ends before the current bin_ends
                    # then we are done with this interval

                    if seg_end < m[1]:
                        break
                    elif seg_start_bin > m[0]:

                    # if the segment starts somewhere inside this bin, skip to the next

                        continue
                    elif seg_start_bin <= m[0] and seg_end > m[1]:

                    # otherwise start

                        seg_start_bin = m[1]
                        if correctly_mapped_count >= 75 \
                            and incorrectly_mapped_count \
                            / correctly_mapped_count <= 0.1:
                            survived_bins.write(str(chromosome) + '\t'
                                    + str(m[0]) + '\t' + str(m[1])
                                    + '\t' + line_r[3] + '\t'
                                    + str(correctly_mapped_count) + '\n'
                                    )
                        else:
                            pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-b', '--input_bed', type=str,
                        help='The bed file with arbitrary segments e.g. inversion callset'
                        , required=True)
    parser.add_argument('-o', '--output', type=str,
                        help='The output file with survived bins for each segment'
                        , required=True)
    parser.add_argument('-i', '--mapping_counts', type=str,
                        help='The mapping track from simulated data already filtered for mapq<10'
                        , required=True)
    args = parser.parse_args()
    counts(args.input_bed, args.mapping_counts, args.output)
