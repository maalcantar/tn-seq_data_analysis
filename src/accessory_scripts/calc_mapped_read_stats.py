#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on July 12 23:15:42 2022
@authors: alcantar and english
run in ngs_processing_source_fastp environment
"""
import argparse

import pandas as pd

import time

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='path to num_read_pairs.txt')

    args = parser.parse_args()

    path_to_read_counts=args.i

    # extrapolate path to deduplication stats from the text file containing
    # sample names and number of reads
    # deduplication stats file is important because it contains information
    # on number of mapped reads pre- and post-deduplication
    path_to_dedup_stats = path_to_read_counts.replace('bowtie2_stats/num_read_pairs.txt',
                                                      'bowtie2_output_dedup/deduplication_stats.csv')
    read_name_to_counts_dict = dict()

    # go through num_read_pairs text file and extract number of reads per sample
    with open(path_to_read_counts) as read_counts_file:
        read_counts_lines = read_counts_file.readlines()
        read_counts_lines = [read.strip('\n') for read in read_counts_lines]
        # skip first line (could also do read_counts_file.readlines()[1:])
        for read in read_counts_lines:
            if read == 'Number of read pairs':
                continue

            read_split = read.split(' : ')
            # current sample name
            read_name_tmp = read_split[0].replace('_out_R1.fastq.gz','').replace('_out_R2.fastq.gz','')
            # reads in current sample
            num_counts_tmp = read_split[1]

            # update dictionary with sample name and number of reads
            read_name_to_counts_dict.update({read_name_tmp: num_counts_tmp})

    # extract number of mapped reads from deduplication stats csv file
    dedup_stats_df = pd.read_csv(path_to_dedup_stats, index_col=0)
    sample_names = list(dedup_stats_df['sample_name'])
    mapped_reads_orig = list(dedup_stats_df['original_read_num'])
    read_name_to_mapped_reads_dict = dict(zip(sample_names, mapped_reads_orig))

    # initialize lists that will contain important alignment stats
    # these will be used to build final dataframe
    sample_name = []
    orig_num_reads = []
    num_mapped_reads = []
    read_to_percent_mapped = []

    # calculate alignment stats for each sample
    for sample in read_name_to_counts_dict:
        orig_num_reads_tmp = int(read_name_to_counts_dict[sample])
        mapped_reads_tmp = int(read_name_to_mapped_reads_dict[sample])

        percentage_of_mapped_reads = mapped_reads_tmp/orig_num_reads_tmp * 100

        sample_name.append(sample)
        orig_num_reads.append(orig_num_reads_tmp)
        num_mapped_reads.append(mapped_reads_tmp)
        read_to_percent_mapped.append(percentage_of_mapped_reads)

    sample_to_mapped_read = pd.DataFrame()
    sample_to_mapped_read['sample_name'] = sample_name
    sample_to_mapped_read['orig_num_reads'] = orig_num_reads
    sample_to_mapped_read['num_mapped_reads'] = num_mapped_reads
    sample_to_mapped_read['percent_mapped_reads'] = read_to_percent_mapped

    out_path = path_to_read_counts.replace('num_read_pairs.txt',
                                           'mapped_stats.csv')
    sample_to_mapped_read.to_csv(out_path)

main()
