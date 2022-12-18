#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on July 29 18:27:51 2022
@authors: alcantar and english
example usage
python 08C_lineage_peak_drawer.py --i ../../tn-seq_data/tn-seq_outputs/mae-013_rb-tn-seq/lineages/mae-013_lineages_parsed_annotated.csv -p
"""

import argparse
import pandas as pd
import numpy as np
import os
import os.path
import glob
import sys

import pickle
from matplotlib import pyplot as plt
import matplotlib.patches as patches

def normalize_array(initial_array, normalization='min_max'):

    '''
    normalizes array

    PARAMETERS
    --------------------
    initial_array: numpy array
        1-dimensional array indicating number of mapped reads at each position
        in the genome
    normalization_type: str
        type of normalization to apply to array
        [min_max: normalize between 0 and 1,
        fraction: divide by number of reads,
        RPKM: calculate RPKM for each position]

    RETURNS
    --------------------
    normalized_array: numpy array
        normalized array
    '''

    if normalization=='min_max':
        normalized_array = (initial_array - np.min(initial_array)) / np.max(initial_array)
    elif normalization=='fraction':
        normalized_array = initial_array/np.sum(initial_array)*100
    elif normalization=='RPKM':
        normalized_array = (testinitial_array_arr/(np.sum(initial_array)/1e6))/1e3
    else:
        print('specified normalization not supported. defaulting to min_max normalization.')
        normalized_array = (initial_array - np.min(initial_array)) / np.max(initial_array)
    return(normalized_array)

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
    mapped_reads_per_pos = np.zeros(init_size, dtype=int)

    for read in read_annotation_df.itertuples():
        start = read[2]-1
        end = read[3]

        mapped_reads_per_pos[start:end] = mapped_reads_per_pos[start:end] + 1

    mapped_reads_per_pos = mapped_reads_per_pos[0:ecoli_genome_len]

    print('finished condition!')
    print('_____________new condition____________________________')

    return(mapped_reads_per_pos)

def extract_barcodes_from_directories(path_to_barcode_directories):
    '''
    Extracts and sorts all barcodes

    PARAMETERS
    --------------------
    path_to_barcode_directories: list of strs
    directories containing barcodes

    RETURNS
    --------------------
    sample_to_read_pos_dict_tmp: dictionary
      dictionary containing mapped reads per genome position for each sample.
      this dictionary is also saved as a pickle file so that this script does
      not have to be run completely every time a new plot is created.
    '''
    curr_barcodes = glob.glob(path_to_barcode_directories + '/*')
    curr_barcodes = [BC.split('/')[-1] for BC in curr_barcodes]

    return(curr_barcodes)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='path to directory lineage csv file, when running program from scratch')
    parser.add_argument('-p', action='store_true', help='run script with pickle files')
    args = parser.parse_args()

    # path to csv file defining lineages
    path_to_lineages = args.i
    pickle_run = args.p

    # initialize common variables
    ecoli_MDS42_genome_len = 3976195
    genome_positions = list(range(1, ecoli_MDS42_genome_len + 1))
    possible_conditions_list = ['LB', 'M9', 'selec.']

    # make output folder path
    if path_to_lineages is not None:
        output_dir = '/'.join(path_to_lineages.split('/')[:-2]) + '/reads_per_position_dedup/'
    else:
        sys.exit('Error! Need to supply either a valid path to sam files or valid path to pickle file. Aborting.')

    # create ouput folders if needed
    CHECK_FOLDER_MAIN = os.path.isdir(output_dir)
    if not CHECK_FOLDER_MAIN:
        os.makedirs(output_dir)
        print(f"created folder : {output_dir}")

    # read dataframe with all metadata
    lineages_df = pd.read_csv(path_to_lineages, index_col=0)

    # if another dataset comes in
    annotated_reads_path_prefix = output_dir.replace('reads_per_position_dedup',
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
        evolution_type = lineage[10]
        extra_peak = lineage[11]

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

        if not pickle_run:
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
        else:
            growth_curr_barcodes = extract_barcodes_from_directories(condition_growth_output_dir)
            no_growth_curr_barcodes = extract_barcodes_from_directories(condition_no_growth_output_dir)

            # handle all barcodes
            LB_barcodes = growth_curr_barcodes + no_growth_curr_barcodes
            LB_barcodes.sort(key=lambda x: int(x.split('_')[1]))
            LB_barcodes = [BC.split('_')[-1] for BC in LB_barcodes]
            # handle no growth barcodes
            selection_barcodes = [BC.split('_')[-1] for BC in growth_curr_barcodes]
        for BC_number, barcode in enumerate(LB_barcodes):

            # if pickle; skip this

            if not pickle_run:
                LB_annotated_reads_BC_tmp_df = LB_annotated_reads_df[LB_annotated_reads_df['clustered_barcode'] == barcode]
                M9_annotated_reads_BC_tmp_df = M9_annotated_reads_df[M9_annotated_reads_df['clustered_barcode'] == barcode]

                LB_mapped_pos = calc_reads_per_pos(LB_annotated_reads_BC_tmp_df, 'LB')
                M9_mapped_pos = calc_reads_per_pos(M9_annotated_reads_BC_tmp_df, 'M9')

                mapped_pos_tmp = [LB_mapped_pos, M9_mapped_pos]

            # check if selection group does not contain current barcode
            # this would indicate NO growth
            # else, the barcode is present, which indicates growth
            if barcode in selection_barcodes:
                if not pickle_run:
                    selection_annotated_reads_BC_tmp_df = selection_annotated_reads_df[selection_annotated_reads_df['clustered_barcode'] == barcode]
                    selection_mapped_pos = calc_reads_per_pos(selection_annotated_reads_BC_tmp_df, 'SELEC.')
                    mapped_pos_tmp.append(selection_mapped_pos)
                    num_conditions = len(mapped_pos_tmp)
                else:
                    num_conditions = 3 # growth case; 3 plots

                BC_tmp_output_dir = condition_growth_output_dir + '/BC_' + str(BC_number) + '_' + barcode
                CHECK_FOLDER_BC_tmp_out = os.path.isdir(BC_tmp_output_dir)
                if not CHECK_FOLDER_BC_tmp_out:
                    os.makedirs(BC_tmp_output_dir)
                    print(f"created folder : {BC_tmp_output_dir}")
            else:
                BC_tmp_output_dir = condition_no_growth_output_dir + '/BC_' + str(BC_number) + '_' + barcode
                CHECK_FOLDER_BC_tmp_out = os.path.isdir(BC_tmp_output_dir)
                if not CHECK_FOLDER_BC_tmp_out:
                    os.makedirs(BC_tmp_output_dir)
                    print(f"created folder : {BC_tmp_output_dir}")

                num_conditions = 2 # no growth case; 2  plots

            pkl_out_path =  BC_tmp_output_dir + '/mapped_reads_per_pos_arrays.pkl'

            if not pickle_run:
                reads_per_position_dict = dict(zip(possible_conditions_list[0:num_conditions],
                                                   mapped_pos_tmp))

           # if pickle present, load pickle file here instead
           # would also need a flag to not read in dataframes

            if not pickle_run:
                print('creating and saving pickle file')
                with open(pkl_out_path, 'wb') as handle:
                    pickle.dump(reads_per_position_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
            else:
                curr_pkl = '/'.join(pkl_out_path.split('/')[6:])
                print(f'loading pickle file: {curr_pkl}')
                with open(pkl_out_path, 'rb') as handle:
                    reads_per_position_dict = pickle.load(handle)
                    print('pickle read')

            max_val = 0
            mapped_pos_norm_tmp = []

            for cond in reads_per_position_dict:
                tmp_array = reads_per_position_dict[cond]
                tmp_array = normalize_array(tmp_array,normalization='fraction' )
                if np.max(tmp_array) > max_val:
                    max_val = np.max(tmp_array)

                # update dictionary with normalized values. this new array is not saved
                # to disk
                reads_per_position_dict.update({cond:tmp_array})

            #`create output folder for figures`
            BC_figs_output_dir = BC_tmp_output_dir + '/figs'
            CHECK_FOLDER_BC_figs_out = os.path.isdir(BC_figs_output_dir)
            if not CHECK_FOLDER_BC_figs_out:
                os.makedirs(BC_figs_output_dir)
                print(f"created folder : {BC_figs_output_dir}")
            BC_figs_output_path = BC_figs_output_dir + '/lineage_peaks'

            title_name = plasmid + ' | ' + lineage_name + ' | ' + 'BC_' + str(BC_number) + '_' + barcode
            output_fname = plasmid + '_' + lineage_name + '_' + 'BC_' + str(BC_number)
            # create line plots
            num_plots = num_conditions
            fig, ax = plt.subplots(num_plots,1,figsize = (6.0,4.0), dpi=400)
            # colors = ['#000000', '#E877AE', '#426BB4','#1C4847', '#409091','#442D85',
            #           '#9375B5', '#851A19','#894E21', '#9A9A9B','#93C945']

            ########## create lineplots ##########
            for plot_num, condition_tmp in enumerate(reads_per_position_dict):
                ax[plot_num].plot(genome_positions, reads_per_position_dict[condition_tmp],
                            color='k', alpha=0.35)
                ax[plot_num].spines['right'].set_visible(False)
                ax[plot_num].spines['top'].set_visible(False)
                ax[plot_num].set_ylabel(condition_tmp, fontsize=12)

            plt.xticks(np.arange(0, 4000001, 1000000))
            ax[0].set_ylim([-0.01, max_val])
            ax[1].set_ylim([-0.01, max_val])

            ax[0].tick_params(
            axis='x',
            which='both',
            bottom=False,
            top=False,
            labelbottom=False)

            # if there are 3 plots, then remove x axis from second plot
            if num_plots == 3:
                ax[1].tick_params(
                axis='x',
                which='both',
                bottom=False,
                top=False,
                labelbottom=False)
                ax[2].set_ylim([-0.01, max_val])
            ax[0].set_title(title_name)

            fig.tight_layout()
            plt.savefig(BC_figs_output_path + '_' + output_fname + '.png')
            plt.savefig(BC_figs_output_path + '_' + output_fname + '.pdf')
            plt.savefig(BC_figs_output_path + '_' + output_fname + '.svg')
            plt.close()

            ########## create insets ###########
            if num_plots == 3:
                print('PLOTTING INSET')
                if extra_peak=='lacI':
                    fig, ax = plt.subplots(figsize = (6.0,4.0), dpi=400)
                    ax.plot(genome_positions, reads_per_position_dict['selec.'],
                                color='k', alpha=0.35)
                    ax.spines['right'].set_visible(False)
                    ax.spines['top'].set_visible(False)
                    ax.set_ylabel('selec. inset', fontsize=12)
                    plt.xlabel('Genome position (bp)', fontsize=12)
                    fig.text(-0.02, 0.5, 'Normalized peak intensity', va='center',
                              rotation='vertical', fontsize=12)
                    max_val_selec_tmp = max(reads_per_position_dict['selec.']) + 0.05

                    plt.xticks(np.arange(289000, 294001, 1000))
                    ax.set_ylim([-0.01, max_val_selec_tmp])
                    ax.set_xlim([292000, 294000])
                    # lacZ coordinates
                    rect1 = patches.Rectangle((289420, 0.02), 3075, 0.04, linewidth=1,
                                                edgecolor='r', facecolor='none')
                    ax.add_patch(rect1)

                    # lacI coordinates
                    rect2 = patches.Rectangle((292617, 0.02), 1083, 0.04, linewidth=1,
                                                edgecolor='r', facecolor='none')
                    ax.add_patch(rect2)
                    fig.tight_layout()
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_lacI_inset' + '.png')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_lacI_inset' + '.pdf')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_lacI_inset' + '.svg')
                    plt.close()
                elif extra_peak=='cpxA':
                    fig, ax = plt.subplots(figsize = (6.0,4.0), dpi=400)
                    ax.plot(genome_positions, reads_per_position_dict['selec.'],
                                color='k', alpha=0.35)
                    ax.spines['right'].set_visible(False)
                    ax.spines['top'].set_visible(False)
                    ax.set_ylabel('selec. inset', fontsize=12)
                    plt.xlabel('Genome position (bp)', fontsize=12)
                    fig.text(-0.02, 0.5, 'Normalized peak intensity', va='center',
                              rotation='vertical', fontsize=12)
                    max_val_selec_tmp = max(reads_per_position_dict['selec.']) + 0.05

                    plt.xticks(np.arange(3532800, 3534200, 1000))
                    ax.set_ylim([-0.01, max_val_selec_tmp])
                    ax.set_xlim([3532515, 3534500])
                    # cpxA coordinates
                    rect1 = patches.Rectangle((3532815, 0.02), 1374, 0.04, linewidth=1,
                                                edgecolor='r', facecolor='none')
                    ax.add_patch(rect1)

                    # cpxR coordinates
                    rect2 = patches.Rectangle((3534185, 0.02), 699, 0.04, linewidth=1,
                                                edgecolor='r', facecolor='none')
                    ax.add_patch(rect2)

                    fig.tight_layout()
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_cpxA_inset' + '.png')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_cpxA_inset' + '.pdf')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_cpxA_inset' + '.svg')
                    plt.close()
                elif extra_peak=='nudH':
                    fig, ax = plt.subplots(figsize = (6.0,4.0), dpi=400)
                    ax.plot(genome_positions, reads_per_position_dict['selec.'],
                                color='k', alpha=0.35)
                    ax.spines['right'].set_visible(False)
                    ax.spines['top'].set_visible(False)
                    ax.set_ylabel('selec. inset', fontsize=12)
                    plt.xlabel('Genome position (bp)', fontsize=12)
                    fig.text(-0.02, 0.5, 'Normalized peak intensity', va='center',
                              rotation='vertical', fontsize=12)
                    max_val_selec_tmp = max(reads_per_position_dict['selec.']) + 0.05

                    plt.xticks(np.arange(2476500, 2481000, 1000))
                    ax.set_ylim([-0.01, max_val_selec_tmp])
                    ax.set_xlim([2476500, 2481000])

                    # ptsP coordinates
                    rect1 = patches.Rectangle((2476723, 0.02), 2247, 0.04, linewidth=1,
                                                edgecolor='r', facecolor='none')
                    ax.add_patch(rect1)

                    # nudH coordinates
                    rect2 = patches.Rectangle((2478982, 0.02), 531, 0.04, linewidth=1,
                                                edgecolor='r', facecolor='none')
                    ax.add_patch(rect2)

                    # mutH coordinates
                    rect3 = patches.Rectangle((2480197, 0.02), 690, 0.04, linewidth=1,
                                                edgecolor='r', facecolor='none')
                    ax.add_patch(rect3)

                    fig.tight_layout()
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_nudH_inset' + '.png')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_nudH_inset' + '.pdf')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_nudH_inset' + '.svg')
                    plt.close()
                else:
                    pass

                fig, ax = plt.subplots(figsize = (6.0,4.0), dpi=400)
                ax.plot(genome_positions, reads_per_position_dict['selec.'],
                            color='k', alpha=0.35)
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.set_ylabel('selec. inset', fontsize=12)
                plt.xlabel('Genome position (bp)', fontsize=12)
                fig.text(-0.02, 0.5, 'Normalized peak intensity', va='center',
                          rotation='vertical', fontsize=12)
                max_val_selec_tmp = max(reads_per_position_dict['selec.']) + 0.05
                if evolution_type == 'arbutin':
                    plt.xticks(np.arange(3332000, 3336001, 1000))
                    ax.set_ylim([-0.01, max_val_selec_tmp])
                    ax.set_xlim([3332900, 3336000])
                    # bglG coordinates
                    rect1 = patches.Rectangle((3334944, 0.02), 836, 0.04, linewidth=1,
                                                edgecolor='r', facecolor='none')
                    ax.add_patch(rect1)

                    # bglF coordinates
                    rect2 = patches.Rectangle((3332933, 0.02), 1877, 0.04, linewidth=1,
                                                edgecolor='r', facecolor='none')
                    ax.add_patch(rect2)

                    fig.tight_layout()
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_inset' + '.png')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_inset' + '.pdf')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_inset' + '.svg')
                    plt.close()

                elif evolution_type == 'serine':
                    plt.xticks(np.arange(1574000, 1579001, 1000))
                    ax.set_ylim([-0.01, max_val_selec_tmp])
                    # ax.set_xlim([1570500, 1583000])
                    ax.set_xlim([1574500, 1579000])

                    # yeaB coordinates
                    rect1 = patches.Rectangle((1576316, 0.02), 578, 0.04, linewidth=1,
                                               edgecolor='r', facecolor='none')
                    ax.add_patch(rect1)

                    # sdaA coordinates
                    rect2 = patches.Rectangle((1577078, 0.02), 1364, 0.04, linewidth=1,
                                               edgecolor='r', facecolor='none')
                    ax.add_patch(rect2)

                    # pabB coordinates
                    rect3 = patches.Rectangle((1574951, 0.02), 1361, 0.04, linewidth=1,
                                               edgecolor='r', facecolor='none')
                    ax.add_patch(rect3)

                    fig.tight_layout()
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_inset' + '.png')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_inset' + '.pdf')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_inset' + '.svg')
                    plt.close()

                elif evolution_type == 'bmdg':
                    # include coordinate for genes of interest of second carbon source
                    # plt.xticks(np.arange(2544800, 2555801, 500))
                    plt.xticks(np.arange(2549800, 2551801, 500)) # true range
                    ax.set_ylim([-0.01, max_val_selec_tmp])
                    # ax.set_xlim([2544800, 2555801])
                    ax.set_xlim([2549800, 2551800])

                    # yqfB coordinates
                    rect1 = patches.Rectangle((2549970, 0.02), 312, 0.04, linewidth=1,
                                               edgecolor='r', facecolor='none')
                    ax.add_patch(rect1)

                    # bglA coordinates
                    rect2 = patches.Rectangle((2550320 , 0.02), 1440, 0.04, linewidth=1,
                                               edgecolor='r', facecolor='none')
                    ax.add_patch(rect2)

                    fig.tight_layout()
                    plt.savefig(BC_figs_output_path + '_' + output_fname + 'yqfB_inset' + '.png')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + 'yqfB_inset' + '.pdf')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + 'yqfB_inset' + '.svg')
                    plt.close()


                    # frdD locus
                    fig, ax = plt.subplots(figsize = (6.0,4.0), dpi=400)
                    ax.plot(genome_positions, reads_per_position_dict['selec.'],
                                color='k', alpha=0.35)
                    ax.spines['right'].set_visible(False)
                    ax.spines['top'].set_visible(False)
                    ax.set_ylabel('selec. inset', fontsize=12)
                    plt.xlabel('Genome position (bp)', fontsize=12)
                    fig.text(-0.02, 0.5, 'Normalized peak intensity', va='center',
                              rotation='vertical', fontsize=12)
                    max_val_selec_tmp = max(reads_per_position_dict['selec.']) + 0.05

                    # include coordinate for genes of interest of second carbon source
                    # plt.xticks(np.arange(2544800, 2555801, 500))
                    plt.xticks(np.arange(3806800, 3809001, 1000)) # true range
                    ax.set_ylim([-0.01, max_val_selec_tmp])
                    # ax.set_xlim([2544800, 2555801])
                    ax.set_xlim([3806800, 3809000])

                    # ampC coordinates
                    rect1 = patches.Rectangle((3806913, 0.02), 1134, 0.04, linewidth=1,
                                               edgecolor='r', facecolor='none')
                    ax.add_patch(rect1)

                    # frdD coordinates
                    rect2 = patches.Rectangle((3808109 , 0.02), 360, 0.04, linewidth=1,
                                               edgecolor='r', facecolor='none')
                    ax.add_patch(rect2)
                    # frdC coordinates
                    rect3 = patches.Rectangle((3808479 , 0.02), 396, 0.04, linewidth=1,
                                               edgecolor='r', facecolor='none')
                    ax.add_patch(rect3)

                    fig.tight_layout()
                    plt.savefig(BC_figs_output_path + '_' + output_fname + 'frdD_inset' + '.png')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + 'frdD_inset' + '.pdf')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + 'frdD_inset' + '.svg')
                    plt.close()

                elif evolution_type == 'gly-glu':

                    plt.xticks(np.arange(2844100, 2846201, 500))
                    ax.set_ylim([-0.01, max_val_selec_tmp])
                    # ax.set_xlim([1570500, 1583000])
                    ax.set_xlim([2844100, 2846200])

                    # sspB coordinates
                    rect1 = patches.Rectangle((2844131, 0.02), 498, 0.04, linewidth=1,
                                               edgecolor='r', facecolor='none')
                    ax.add_patch(rect1)

                    # sspA coordinates
                    rect2 = patches.Rectangle((2844634, 0.02), 639, 0.04, linewidth=1,
                                               edgecolor='r', facecolor='none')
                    ax.add_patch(rect2)

                    # rpsI coordinates
                    rect3 = patches.Rectangle((2845667, 0.02), 393, 0.04, linewidth=1,
                                               edgecolor='r', facecolor='none')
                    ax.add_patch(rect3)

                    fig.tight_layout()
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_inset' + '.png')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_inset' + '.pdf')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_inset' + '.svg')
                    plt.close()

                elif evolution_type == 'serine_bmdg':
                    # L-serine
                    plt.xticks(np.arange(1574000, 1579001, 1000))
                    ax.set_ylim([-0.01, max_val_selec_tmp])
                    # ax.set_xlim([1570500, 1583000])
                    ax.set_xlim([1574500, 1579000])

                    # yeaB coordinates
                    rect1 = patches.Rectangle((1576316, 0.02), 578, 0.04, linewidth=1,
                                               edgecolor='r', facecolor='none')
                    ax.add_patch(rect1)

                    # sdaA coordinates
                    rect2 = patches.Rectangle((1577078, 0.02), 1364, 0.04, linewidth=1,
                                               edgecolor='r', facecolor='none')
                    ax.add_patch(rect2)

                    # pabB coordinates
                    rect3 = patches.Rectangle((1574951, 0.02), 1361, 0.04, linewidth=1,
                                               edgecolor='r', facecolor='none')
                    ax.add_patch(rect3)

                    fig.tight_layout()
                    plt.savefig(BC_figs_output_path + '_' + output_fname + 'serine_inset' + '.png')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + 'serine_inset' + '.pdf')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + 'serine_inset' + '.svg')
                    plt.close()
                    # bmdg
                    fig, ax = plt.subplots(figsize = (6.0,4.0), dpi=400)
                    ax.plot(genome_positions, reads_per_position_dict['selec.'],
                                color='k', alpha=0.35)
                    ax.spines['right'].set_visible(False)
                    ax.spines['top'].set_visible(False)
                    ax.set_ylabel('selec. inset', fontsize=12)
                    plt.xlabel('Genome position (bp)', fontsize=12)
                    fig.text(-0.02, 0.5, 'Normalized peak intensity', va='center',
                              rotation='vertical', fontsize=12)
                    max_val_selec_tmp = max(reads_per_position_dict['selec.']) + 0.05

                    # include coordinate for genes of interest of second carbon source
                    # plt.xticks(np.arange(2544800, 2555801, 500))
                    plt.xticks(np.arange(2549800, 2551801, 500)) # true range
                    ax.set_ylim([-0.01, max_val_selec_tmp])
                    # ax.set_xlim([2544800, 2555801])
                    ax.set_xlim([2549800, 2551800])

                    # yqfB coordinates
                    rect1 = patches.Rectangle((2549970, 0.02), 312, 0.04, linewidth=1,
                                               edgecolor='r', facecolor='none')
                    ax.add_patch(rect1)

                    # bglA coordinates
                    rect2 = patches.Rectangle((2550320 , 0.02), 1440, 0.04, linewidth=1,
                                               edgecolor='r', facecolor='none')
                    ax.add_patch(rect2)

                    fig.tight_layout()
                    plt.savefig(BC_figs_output_path + '_' + output_fname + 'bmdg_inset' + '.png')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + 'bmdg_inset' + '.pdf')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + 'bmdg_inset' + '.svg')
                    plt.close()

                    # frdD locus
                    fig, ax = plt.subplots(figsize = (6.0,4.0), dpi=400)
                    ax.plot(genome_positions, reads_per_position_dict['selec.'],
                                color='k', alpha=0.35)
                    ax.spines['right'].set_visible(False)
                    ax.spines['top'].set_visible(False)
                    ax.set_ylabel('selec. inset', fontsize=12)
                    plt.xlabel('Genome position (bp)', fontsize=12)
                    fig.text(-0.02, 0.5, 'Normalized peak intensity', va='center',
                              rotation='vertical', fontsize=12)
                    max_val_selec_tmp = max(reads_per_position_dict['selec.']) + 0.05

                    # include coordinate for genes of interest of second carbon source
                    # plt.xticks(np.arange(2544800, 2555801, 500))
                    plt.xticks(np.arange(3806800, 3809001, 1000)) # true range
                    ax.set_ylim([-0.01, max_val_selec_tmp])
                    # ax.set_xlim([2544800, 2555801])
                    ax.set_xlim([3806800, 3809000])

                    # ampC coordinates
                    rect1 = patches.Rectangle((3806913, 0.02), 1134, 0.04, linewidth=1,
                                               edgecolor='r', facecolor='none')
                    ax.add_patch(rect1)

                    # frdD coordinates
                    rect2 = patches.Rectangle((3808109 , 0.02), 360, 0.04, linewidth=1,
                                               edgecolor='r', facecolor='none')
                    ax.add_patch(rect2)
                    # frdC coordinates
                    rect3 = patches.Rectangle((3808479 , 0.02), 396, 0.04, linewidth=1,
                                               edgecolor='r', facecolor='none')
                    ax.add_patch(rect3)

                    fig.tight_layout()
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_frdD_inset' + '.png')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_frdD_inset' + '.pdf')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_frdD_inset' + '.svg')
                    plt.close()

                    # cpxA
                    fig, ax = plt.subplots(figsize = (6.0,4.0), dpi=400)
                    ax.plot(genome_positions, reads_per_position_dict['selec.'],
                                color='k', alpha=0.35)
                    ax.spines['right'].set_visible(False)
                    ax.spines['top'].set_visible(False)
                    ax.set_ylabel('selec. inset', fontsize=12)
                    plt.xlabel('Genome position (bp)', fontsize=12)
                    fig.text(-0.02, 0.5, 'Normalized peak intensity', va='center',
                              rotation='vertical', fontsize=12)
                    max_val_selec_tmp = max(reads_per_position_dict['selec.']) + 0.05

                    plt.xticks(np.arange(3532800, 3534200, 1000))
                    ax.set_ylim([-0.01, max_val_selec_tmp])
                    ax.set_xlim([3532515, 3534500])
                    # cpxA coordinates
                    rect1 = patches.Rectangle((3532815, 0.02), 1374, 0.04, linewidth=1,
                                                edgecolor='r', facecolor='none')
                    ax.add_patch(rect1)

                    # cpxR coordinates
                    rect2 = patches.Rectangle((3534185, 0.02), 699, 0.04, linewidth=1,
                                                edgecolor='r', facecolor='none')
                    ax.add_patch(rect2)

                    fig.tight_layout()
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_cpxA_inset' + '.png')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_cpxA_inset' + '.pdf')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_cpxA_inset' + '.svg')
                    plt.close()

                    # lrp
                    fig, ax = plt.subplots(figsize = (6.0,4.0), dpi=400)
                    ax.plot(genome_positions, reads_per_position_dict['selec.'],
                                color='k', alpha=0.35)
                    ax.spines['right'].set_visible(False)
                    ax.spines['top'].set_visible(False)
                    ax.set_ylabel('selec. inset', fontsize=12)
                    plt.xlabel('Genome position (bp)', fontsize=12)
                    fig.text(-0.02, 0.5, 'Normalized peak intensity', va='center',
                              rotation='vertical', fontsize=12)
                    max_val_selec_tmp = max(reads_per_position_dict['selec.']) + 0.05

                    plt.xticks(np.arange(700500, 703501, 1000))
                    ax.set_ylim([-0.01, max_val_selec_tmp])
                    ax.set_xlim([700500, 703500])
                    # trxB coordinates
                    rect1 = patches.Rectangle((700923, 0.02), 267, 0.04, linewidth=1,
                                                edgecolor='r', facecolor='none')
                    ax.add_patch(rect1)

                    # # cpxA coordinates
                    # rect2 = patches.Rectangle((795328, 0.02), 495, 0.04, linewidth=1,
                    #                             edgecolor='r', facecolor='none')
                    # ax.add_patch(rect2)

                    fig.tight_layout()
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_lrp_inset' + '.png')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_lrp_inset' + '.pdf')
                    plt.savefig(BC_figs_output_path + '_' + output_fname + '_lrp_inset' + '.svg')
                    plt.close()

                else:
                    print(f'evolution type {evolution_type} not recognized. No inset will be plotted.')
                    pass

            ########## create heatmaps ##########
            BC_heatmaps_output_path = BC_figs_output_dir + '/heatmaps'

            fig, ax = plt.subplots(num_plots,1,figsize = (8.0,2.0), dpi=400)
            # colors = ['#000000', '#E877AE', '#426BB4','#1C4847', '#409091','#442D85',
            #           '#9375B5', '#851A19','#894E21', '#9A9A9B','#93C945']

            for plot_num, condition_tmp in enumerate(reads_per_position_dict):
                n = 20000 #10000 #bp -- cluster size
                tmp_avg =  np.array([max(reads_per_position_dict[condition_tmp][i:i+n]) for i in range(0,len(reads_per_position_dict[condition_tmp]),n)])

                h=ax[plot_num].imshow(tmp_avg[np.newaxis,:], aspect="auto",
                vmin = 0,vmax = 0.09, cmap = 'Greys') #cmap = 'Greys', 'gray'
                ax[plot_num].set_yticks([])
                x_label_list = ['0','1e6', '2e6', '3e6'] #, '4e6'
                ax[plot_num].set_xticks([0, 50,100,150]) #,400[100,200,300 # this needs to be changed depending on cluster size
                ax[plot_num].set_xticklabels(x_label_list)
                ax[plot_num].set_ylabel(condition_tmp, fontsize=12)
            fig.colorbar(h)

            ax[0].tick_params(
            axis='x',
            which='both',
            bottom=True,
            top=False,
            labelbottom=False)
            ax[0].set_title(title_name)

            if num_plots == 3:
                ax[1].tick_params(
                axis='x',
                which='both',
                bottom=True,
                top=False,
                labelbottom=False)
                ax[2].set_ylim([-0.01, max_val])
            ax[-1].set_xlabel('Genome position (bp)', fontsize=12)

            fig.tight_layout()
            plt.savefig(BC_heatmaps_output_path + '_' + output_fname + '.png')
            plt.savefig(BC_heatmaps_output_path + '_' + output_fname + '.pdf')
            plt.savefig(BC_heatmaps_output_path + '_' + output_fname + '.svg')
            plt.close()

if __name__ == "__main__":
    main()
