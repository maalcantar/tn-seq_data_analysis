#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on August 1 14:11:52 2022
@authors: alcantar and english
example usage
python 09C_insert_site_plot_generator.py --i ../../tn-seq_data/tn-seq_outputs/mae-011_rb-tn-seq/lineages/mae-011_lineages_parsed_annotated.csv
"""

import argparse
import pandas as pd
import numpy as np
import os
import os.path
import glob
import sys
import csv

import pickle
from matplotlib import pyplot as plt
import matplotlib.patches as patches


def calc_reads_per_pos(read_annotation_df, sample_name,
                      ecoli_genome_len=3976195, init_size=4000000):

    '''
    for each sample, calculate number of reads mapping to each position
    in the genome. results are stored as a numpy arrays inside of a dictionary

    PARAMETERS
    --------------------
    read_annotation_df: str
        dataframe with annotated barcode cluster file
        that contains the following columns:
        0) read_name 1) start 2) end 3) barcode 4) flag 5) insert size 6) gene1
        7)gene2 8)gene1_barcode 9) gene2_barcode 10) cluster_id
        11) clustered_barcode 12) gene1_clust_id 13) gene2_clust_id 14) gene1_clust_bc
        15) gene2_clust_bc
    sample_name: str
        current sample name ['LB', 'M9', 'selec.']
    ecoli_genome_len: int
        size of e coli genome in bp. this will change depending on which strain
        is being used. default is genome size of E. coli MDS42
    init_size: int
        size of initialized array. this should be larger than ecoli_genome_len
        to account for any soft-clipping that may have occurred with read mapper.
        default is 4000000

    RETURNS
    --------------------
    sample_to_read_pos_dict_tmp: dictionary
        dictionary containing mapped reads per genome position for each sample.
        this dictionary is also saved as a pickle file so that this script does
        not have to be run completely every time a new plot is created.
    '''

    print(f"starting condition: {sample_name}")

    print('creating array containing counts per position')
    mapped_reads_per_pos_F = np.zeros(init_size, dtype=int)
    mapped_reads_per_pos_R = np.zeros(init_size, dtype=int)

    for read in read_annotation_df.itertuples():
        read_name = read[1]
        start = read[2]-1
        # end = read[3] [start:end]

        if read_name[-1] == 'F':
            mapped_reads_per_pos_F[start] = mapped_reads_per_pos_F[start] + 1
        elif read_name[-1] == 'R':
            mapped_reads_per_pos_R[start] = mapped_reads_per_pos_R[start] + 1

    mapped_reads_per_pos_F = mapped_reads_per_pos_F[0:ecoli_genome_len]
    mapped_reads_per_pos_R = mapped_reads_per_pos_R[0:ecoli_genome_len]

    print('finished condition!')
    print('_____________new condition____________________________')

    return(mapped_reads_per_pos_F, mapped_reads_per_pos_R)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='path to directory lineage csv file, when running program from scratch')
    parser.add_argument('-p', action='store_true', help='run script with pickle files')
    args = parser.parse_args()

    # path to csv file defining lineages
    path_to_lineages = args.i

    # initialize common variables
    ecoli_MDS42_genome_len = 3976195
    genome_positions = list(range(1, ecoli_MDS42_genome_len + 1))
    possible_conditions_list = ['LB', 'M9', 'selec.']

    # make output folder path
    if path_to_lineages is not None:
        output_dir = '/'.join(path_to_lineages.split('/')[:-2]) + '/reads_per_position_dedup_insert_site_plots/'
    else:
        sys.exit('Error! Need to supply either a valid path to sam files or valid path to pickle file. Aborting.')

    # create ouput folders if needed
    CHECK_FOLDER_MAIN = os.path.isdir(output_dir)
    if not CHECK_FOLDER_MAIN:
        os.makedirs(output_dir)
        print(f"created folder : {output_dir}")

    # read dataframe with all metadata
    lineages_df = pd.read_csv(path_to_lineages, index_col=0)

    # hard-coded; consider determining paths programatically
    # if another dataset comes in
    annotated_reads_path_prefix = output_dir.replace('reads_per_position_dedup_insert_site_plots',
                                                     'barcode_clustering')
    #'../../tn-seq_data/tn-seq_outputs/mae-011_rb-tn-seq/barcode_clustering/'
    annotated_reads_path_suffix = '_clustered_barcode_annotated_reads.csv'

    # loop through each row (i.e., lineage) in dataframe
    for lineage in lineages_df.itertuples():

        # extract lineage features from dataframe
        lineage_tmp = lineage[1]
        min_count_thresh = lineage[9]
        lineage_name = lineage[5]
        plasmid = lineage[2]

        # based on plasmid, set how many barcodes should be considered
        if plasmid == 'pMP-001_2':
            max_num_barcodes = 14
        elif plasmid == 'pMP-001_8':
            max_num_barcodes = 9
        else:
            max_num_barcodes = 9

        plasmid_output_dir = output_dir + plasmid
        # create ouput folders each plasmid
        CHECK_FOLDER_plas_out = os.path.isdir(plasmid_output_dir)
        if not CHECK_FOLDER_plas_out:
            os.makedirs(plasmid_output_dir)
            print(f"created folder : {plasmid_output_dir}")
        condition_output_dir = plasmid_output_dir + '/' + lineage_name

        # create ouput folder for each condition
        CHECK_FOLDER_cond_out = os.path.isdir(condition_output_dir)
        if not CHECK_FOLDER_cond_out:
            os.makedirs(condition_output_dir)
            print(f"created folder : {condition_output_dir}")

        # create folder where lineages that will contain data for lineages
        # that grew
        condition_growth_output_dir = condition_output_dir + '/growth'
        CHECK_FOLDER_cond_growth_out = os.path.isdir(condition_growth_output_dir)
        if not CHECK_FOLDER_cond_growth_out:
            os.makedirs(condition_growth_output_dir)
            print(f"created folder : {condition_growth_output_dir}")

        # create folder where lineages that will contain data for lineages
        # that did NOT grew
        condition_no_growth_output_dir = condition_output_dir + '/no_growth'
        CHECK_FOLDER_cond_no_growth_out = os.path.isdir(condition_no_growth_output_dir)
        if not CHECK_FOLDER_cond_no_growth_out:
            os.makedirs(condition_no_growth_output_dir)
            print(f"created folder : {condition_no_growth_output_dir}")

        LB_time_point_sample_name = '_'.join(lineage[6].split('_')[0:3])
        M9_time_point_sample_name = '_'.join(lineage[7].split('_')[0:3])
        selection_time_point_sample_name = '_'.join(lineage[8].split('_')[0:3])

        # if pickle; skip this

        # read in relevant annotated read maps dataframe
        # contains alignments, coordinates, and their corresponding barcodes
        LB_annotated_reads_path = annotated_reads_path_prefix + \
        LB_time_point_sample_name + \
        annotated_reads_path_suffix

        M9_annotated_reads_path = annotated_reads_path_prefix + \
        M9_time_point_sample_name + \
        annotated_reads_path_suffix

        selection_annotated_reads_path = annotated_reads_path_prefix + \
        selection_time_point_sample_name + \
        annotated_reads_path_suffix

        LB_annotated_reads_df = pd.read_csv(LB_annotated_reads_path,
                                            index_col=0,
                                            low_memory=False)
        M9_annotated_reads_df = pd.read_csv(M9_annotated_reads_path,
                                            index_col=0,
                                            low_memory=False)
        selection_annotated_reads_df = pd.read_csv(selection_annotated_reads_path,
                                                   index_col=0,
                                                   low_memory=False)

        # if pickle; extract barcodes from path; use glob.glob
        # would need to sort by barcode number: https://stackoverflow.com/questions/17555218/python-how-to-sort-a-list-of-lists-by-the-fourth-element-in-each-list
        # barcodes from founding colonies
        LB_barcodes = list(LB_annotated_reads_df['clustered_barcode'].value_counts().index[0:max_num_barcodes])
        # barcodes from M9 growth -- manually checked and should always be the same as LB_barcodes
        M9_barcodes = list(M9_annotated_reads_df['clustered_barcode'].value_counts().index[0:max_num_barcodes])
        # barcodes after selection -- only take those above defined threshold
        selection_barcodes = selection_annotated_reads_df['clustered_barcode'].value_counts()#.index[0:max_num_barcodes]
        selection_barcodes = list(selection_barcodes[selection_barcodes >= min_count_thresh-1].index)

        for BC_number, barcode in enumerate(LB_barcodes):

            print(f'starting barcode: {barcode}')
            # if pickle; skip this
            LB_annotated_reads_BC_tmp_df = LB_annotated_reads_df[LB_annotated_reads_df['clustered_barcode'] == barcode]
            M9_annotated_reads_BC_tmp_df = M9_annotated_reads_df[M9_annotated_reads_df['clustered_barcode'] == barcode]

            LB_mapped_pos_F, LB_mapped_pos_R = calc_reads_per_pos(LB_annotated_reads_BC_tmp_df, 'LB')
            M9_mapped_pos_F, M9_mapped_pos_R  = calc_reads_per_pos(M9_annotated_reads_BC_tmp_df, 'M9')

            # mapped_pos_tmp = [LB_mapped_pos, M9_mapped_pos]

            # check if selection group does not contain current barcode
            # this would indicate NO growth
            # else, the barcode is present, which indicates growth
            if barcode in selection_barcodes:
                selection_annotated_reads_BC_tmp_df = selection_annotated_reads_df[selection_annotated_reads_df['clustered_barcode'] == barcode]
                selec_mapped_pos_F, selec_mapped_pos_R = calc_reads_per_pos(selection_annotated_reads_BC_tmp_df, 'SELEC.')
                BC_tmp_output_dir = condition_growth_output_dir + '/BC_' + str(BC_number) + '_' + barcode
                CHECK_FOLDER_BC_tmp_out = os.path.isdir(BC_tmp_output_dir)

                if not CHECK_FOLDER_BC_tmp_out:
                    os.makedirs(BC_tmp_output_dir)
                    print(f"created folder : {BC_tmp_output_dir}")
            else:
                continue
                # BC_tmp_output_dir = condition_no_growth_output_dir + '/BC_' + str(BC_number) + '_' + barcode
                # CHECK_FOLDER_BC_tmp_out = os.path.isdir(BC_tmp_output_dir)
                # if not CHECK_FOLDER_BC_tmp_out:
                #     os.makedirs(BC_tmp_output_dir)
                #     print(f"created folder : {BC_tmp_output_dir}")


            print('creating and saving insert_site_plot' )
            M9_out_path =  BC_tmp_output_dir + '/' + lineage_name + '_' + 'BC_' + str(BC_number) + '_M9.insert_site_plot'
            selec_out_path =  M9_out_path.replace('M9.insert_site_plot', 'selec.insert_site_plot')

            with open(M9_out_path, 'w') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerows(zip(M9_mapped_pos_F,M9_mapped_pos_R))

            with open(selec_out_path, 'w') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerows(zip(selec_mapped_pos_F,selec_mapped_pos_R))


if __name__ == "__main__":
    main()
