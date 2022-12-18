#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Feb 21 12:10:30 2022
@authors: alcantar and english
"""
import argparse

import pandas as pd
import subprocess
import glob

def deduplication_stats(path_to_dedup_bams):

    '''
    produces a dataframe summarizing deduplication statistics: original number of reads, reads after
    deduplication, reads removed, and duplication rate (reads removed / original number of reads). These
    stats are taken from umi_tools output files.

    PARAMETERS
    --------------------
    path_to_dedup_bams: str
        path to directory containing deduplicated bam files and umi_tools output

    RETURNS
    --------------------
    dedup_log_df: pandas dataframe
        dataframe summarizing deduplication statistics
    '''

    dedup_log_files = glob.glob(path_to_dedup_bams + '*_dedup_log')
    dedup_log_files.sort()
    dedup_log_dict_list = []
    for dedup_log_path in dedup_log_files:

        sample_name = dedup_log_path.split('/')[-1].replace('_dedup_log','')

        # use command line to extract lines with read information
        # need to be done sequentially
        ps1 = subprocess.Popen(('cat', dedup_log_path), stdout=subprocess.PIPE)
        reads_input_line= subprocess.check_output(('grep', 'Input Reads:'), stdin=ps1.stdout)
        ps1.wait()
        ps2 = subprocess.Popen(('cat', dedup_log_path), stdout=subprocess.PIPE)
        reads_output_line = subprocess.check_output(('grep', 'reads out: '), stdin=ps2.stdout)
        ps2.wait()

        ## code for jupyter notebook ##
        #  reads_input_line = !cat $dedup_log_path | grep 'Input Reads:'
        # reads_output_line = !cat $dedup_log_path | grep 'reads out: '
        #int(reads_input_line[0].split(':')[4].split(',')[0].rstrip().lstrip())
        #int(reads_output_line[0].split(':')[-1].rstrip().lstrip())

        input_read_num = int(str(reads_input_line).split(':')[4].split(',')[0].rstrip().lstrip())
        output_read_num = int(str(reads_output_line).split(':')[-1].replace('\\n\'','').rstrip().lstrip())
        reads_removed = input_read_num - output_read_num
        deduplication_rate = (reads_removed/input_read_num) * 100

        dedup_log_dict_list.append({'sample_name': sample_name,
                                    'original_read_num': input_read_num,
                                   'dedup_read_num': output_read_num,
                                   'reads_removed': reads_removed,
                                   'duplication_rate': deduplication_rate})

    dedup_log_df = pd.DataFrame(dedup_log_dict_list)

    return(dedup_log_df)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='directory containing deduplicated bamfiles')
    args = parser.parse_args()

    path_to_dedup_bams=args.i
    dedup_log_df = deduplication_stats(path_to_dedup_bams)
    dedup_log_df.to_csv(path_to_dedup_bams + 'deduplication_stats.csv')

if __name__ == "__main__":
    main()
