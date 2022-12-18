#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on June 7 11:59:39 2022
@authors: alcantar and english
run in ngs_processing_source_fastp environment
"""
import argparse

import pandas as pd
import numpy as np
import os.path
import time
import subprocess
import glob
import os
import sys
import tqdm

def add_mapped_end_points(barcode_pos, output_path_prefix):

    '''
    deduces end points of mapped reads from left-most mapped read position and insert size.
    appends end points to original dataframe

    PARAMETERS
    --------------------
    barcode_pos_path: str
        path to dataframe containing Rb-Tnseq information
        at the bare minimum, should contain: read_name_bc, barcode, flag, start, and insert_size
    output_path_prefix: str
        prefix for output path
    RETURNS
    --------------------
    None -- just saves the following two files:
    BC_positions_df: pandas dataframe
        dataframe containing Rb-Tnseq information and read end points.
        output suffix is _barcode_positions_final.csv
    BC_positions_bed_df: pandas dataframe
        same as BC_positions_df except in bed format
        col1: chromosome name; col2: start point, col3: end point
        col4: unique name
        output suffix is _barcode_positions_bed_final.bed
    '''

    # read in barcode lookup file
    BC_positions_df = pd.read_csv(barcode_pos)

    # obtain unique read names and initialize dictionary that maps
    # read names to end points
    read_names_unique = pd.unique(BC_positions_df['read_name_bc'])
    read_to_end_dict = dict()

    # set read_name_bc as index to speed up indexing
    BC_positions_df_idxed = BC_positions_df.set_index('read_name_bc')

    ########## START DATAFRAME MERGING PROCESS ##########
    new_df = BC_positions_df_idxed.reset_index()[['read_name_bc','barcode','flag','start', 'insert_size']]
    new_df['start'] = new_df['start'].astype('str')
    new_df['insert_size'] = new_df['insert_size'].astype('str')
    combined_new_df = new_df.groupby(['read_name_bc']).agg({'insert_size': '#'.join,'start':'#'.join})
    ########## END DATAFRAME MERGING PROCESS ##########

    # creating mappings between read names and other dataframe attributes
    read_names_list = list(combined_new_df.index)
    insert_size_list = list(combined_new_df['insert_size'])
    start_points_list = list(combined_new_df['start'])

    # create dictionaries with mappings
    read_name_to_insert_size_dict = dict(zip(read_names_list, insert_size_list))
    read_name_to_start_points_dict = dict(zip(read_names_list, start_points_list))

    # loop through all unique read names
    for read_name in tqdm.tqdm(read_names_unique):

        # extract insert size and starting point for current read name
        insert_size_tmp = read_name_to_insert_size_dict[read_name]
        start_point_tmp = read_name_to_start_points_dict[read_name]

        # split strings based on delimiters
        insert_size_tmp_split = insert_size_tmp.split('#')
        start_point_tmp_split = start_point_tmp.split('#')
        insert_size_tmp_split = list(map(int, insert_size_tmp_split))
        start_point_tmp_split = list(map(int, start_point_tmp_split))

        insert_size_tmp_split_sorted, start_point_tmp_split_sorted = (list(t) for t in zip(*sorted(zip(insert_size_tmp_split, start_point_tmp_split))))

        # estimate end points based on different number of reads in a pair mapping
        if len(insert_size_tmp_split_sorted) == 2:
            # case where both reads map. use insert size and left-end position from
            # forward read to estimate end point
            end_site = start_point_tmp_split_sorted[1] + insert_size_tmp_split_sorted[1]

        elif len(insert_size_tmp_split_sorted) == 1:
            # case where only 1 read mapped. use only available left-end points
            # and insert size to estimate end site
            end_site = start_point_tmp_split_sorted[0] + insert_size_tmp_split_sorted[0]

        else:
            # case where nothing maps. should not occur unless there's an error with
            # extracting unique read names, or if something mapped twice
            print('Error! No matching reads found or read mapped to >1 location')
        read_to_end_dict.update({read_name: end_site})

    # initialize list with end points and read directions
    # read directions will help obtain unique read names for each read
    end_sites = []
    read_directions = []

    # loop through all reads
    for read in tqdm.tqdm(BC_positions_df.itertuples()):

        ## read tuple details
        # 0 - Index
        # 1 - _1
        # 2 - read_name_bc
        # 3 - barcode
        # 4 - flag
        # 5 - start
        # 6 - insert_size

        read_name = read[2]
        insert_size = read[6]
        barcode = read[3]

        # assign read as being forward or reverse read.
        if insert_size > 0:
            read_direction = '_F'
        else:
            read_direction = '_R'

        read_with_strand = read_name + read_direction
        read_directions.append(read_with_strand)

        end_site_tmp = read_to_end_dict[read_name]
        end_sites.append(end_site_tmp)
    BC_positions_df['end'] = end_sites
    BC_positions_df['read_name_with_dir'] = read_directions
    BC_positions_df['chrom'] = 'AP012306.1'

    # save most informative columns
    BC_positions_df = BC_positions_df[['read_name_bc', 'read_name_with_dir',
                                   'barcode', 'flag', 'start', 'end',
                                   'insert_size','chrom']]

    # sort by starting position (which is required by bedtools closest)
    BC_positions_df = BC_positions_df.sort_values(by='start', ignore_index=True)
    BC_positions_bed_df = BC_positions_df[['chrom', 'start', 'end', 'read_name_with_dir']]


    BC_positions_final_path = output_path_prefix+'_barcode_positions_final.csv'
    BC_positions_bed_final_path = output_path_prefix+'_barcode_positions_bed_final.bed'

    BC_positions_df.to_csv(BC_positions_final_path)
    BC_positions_bed_df.to_csv(BC_positions_bed_final_path, header=False, index=False, sep='\t')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='barcode_positions_lookup file')
    parser.add_argument('-o', help='output_directory')
    args = parser.parse_args()

    barcode_pos_lookup=args.i
    output_path_prefix = args.o

    # make sure output path does NOT end in a forward slash
    if output_path_prefix[-1] == '/':
        output_path_prefix = output_path_prefix[:-1]

    add_mapped_end_points(barcode_pos_lookup, output_path_prefix)

if __name__ == "__main__":
    main()
