#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on June 21 09:40:37 2022
@authors: alcantar and english
"""
import argparse

import pandas as pd
import numpy as np
import os
import os.path
import glob
import tqdm
import sys

from Bio import SeqIO

import pickle
from matplotlib import pyplot as plt
import matplotlib.patches as patches

'''
example usage (terminal):
# running directly from pickle file
python 08_peak_drawer.py --p ../../tn-seq_data/tn-seq_outputs/mae-005/reads_per_position_dedup/mapped_reads_per_pos_arrays.pkl --startpoints ../../tn-seq_data/tn-seq_outputs/mae-005/reads_per_position_dedup/startpoints.txt --endpoints ../../tn-seq_data/tn-seq_outputs/mae-005/reads_per_position_dedup/endpoints.txt -e arbutin -l 3 4
# running from mapped sam files
python 08_peak_drawer.py --s ../../tn-seq_data/tn-seq_outputs/mae-007/bowtie2_output_dedup/ --startpoints ../../tn-seq_data/tn-seq_outputs/mae-007/reads_per_position_dedup/controls_pENG135.txt --endpoints ../../tn-seq_data/tn-seq_outputs/mae-007/reads_per_position_dedup/conditions_pENG135_noAtc.txt -o pENG135_noAtc
'''
def sam_to_df(mapped_reads_sam_path):

    '''
    converts sam file to a pandas dataframe. this function is adapted adapted from
    our map_barcodes_positions.py script

    PARAMETERS
    --------------------
    mapped_reads_sam: str
        path to sam file. this sam file should only contain mapped reads (can be
        obtained by running: samtools view -F 4)

     RETURNS
    --------------------
    sam_mapped_df: pandas dataframe
        dataframe containing information from sam file with only mapped reads
    '''


    print('converting SAM file to dataframe')

    # loop through each line (i.e., mapped read) in sam file
    with open(mapped_reads_sam_path, 'r') as temp_f:
        # get number of columns in each line
        col_count = [ len(l.split("\t")) for l in temp_f.readlines() ]

    # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
    column_names = [i for i in range(0, max(col_count))]

    # Read csv
    sam_mapped_df = pd.read_csv(mapped_reads_sam_path, sep="\t", header=None,\
                                 names=column_names, low_memory=False)

    # extract relevant columns
    sam_mapped_df = sam_mapped_df[[0,1,3,8,9]]
    sam_mapped_df.columns = ['read_name', 'flag','start',
                             'insert_size','sequence']

    return(sam_mapped_df)

def add_mapped_end_points(sam_df):

    '''
    deduces end points of mapped reads from left-most mapped read position
    and insert size. appends end points to original dataframe

    PARAMETERS
    --------------------
    sam_df: str
        path to dataframe containing Rb-Tnseq information
        at the bare minimum, should contain: read_name, barcode, flag, start,
        and insert_size
    output_path_prefix: str
        prefix for output path
    RETURNS
    --------------------
    None -- just saves the following two files:
    sam_df: pandas dataframe
        dataframe containing Rb-Tnseq information and read end points.
        output suffix is _barcode_positions_final.csv
    '''

    print('calculating mapped location')

    # obtain unique read names and initialize dictionary that maps
    # read names to end points
    read_names_unique = pd.unique(sam_df['read_name'])
    read_to_end_dict = dict()

    # set read_name as index to speed up indexing
    sam_df_idxed = sam_df.set_index('read_name')

    ########## START DATAFRAME MERGING PROCESS ##########
    new_df = sam_df_idxed.reset_index()[['read_name','flag','start', 'insert_size']]
    new_df['start'] = new_df['start'].astype('str')
    new_df['insert_size'] = new_df['insert_size'].astype('str')
    combined_new_df = new_df.groupby(['read_name']).agg({'insert_size': '#'.join,'start':'#'.join})
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

        insert_size_tmp_split_sorted, start_point_tmp_split_sorted = (list(t)
        for t in zip(*sorted(zip(insert_size_tmp_split, start_point_tmp_split))))

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
    for read in tqdm.tqdm(sam_df.itertuples()):
        # read tuple details
        # 0 - Index
        # 1 - read_name
        # 2 - flag
        # 3 - start
        # 4 - insert_size
        # 5 - sequence

        read_name = read[1]
        insert_size = read[4]

        # assign read as being forward or reverse read.
        if insert_size > 0:
            read_direction = '_F'
        else:
            read_direction = '_R'

        read_with_strand = read_name + read_direction
        read_directions.append(read_with_strand)

        end_site_tmp = read_to_end_dict[read_name]
        end_sites.append(end_site_tmp)

    sam_df['end'] = end_sites
    sam_df['read_name_with_dir'] = read_directions
    sam_df['chrom'] = 'AP012306.1'

    # save most informative columns
    sam_df = sam_df[['read_name', 'read_name_with_dir',
                                      'flag', 'start', 'end',
                                   'insert_size','chrom']]

    # sort by starting position
    sam_final_df = sam_df.sort_values(by='start', ignore_index=True)
#     BC_positions_bed_df = sam_df[['chrom', 'start', 'end', 'read_name_with_dir']]

    return(sam_final_df)

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

def calc_reads_per_pos(mapped_sam_files, output_path,
                      ecoli_genome_len=3976195, init_size=4000000):

    '''
    for each sample, calculate number of reads mapping to each position
    in the genome. results are stored as a numpy arrays inside of a dictionary

    PARAMETERS
    --------------------
    mapped_sam_files: str
        list containing paths to sam files with mapped reads
    output_path: str
        output path for dictionary containing results
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

    sample_to_read_pos_dict = dict()
    for sam_file_path in mapped_sam_files:

        current_sample = sam_file_path.split('/')[-1]
        print(f"starting sample: {current_sample}")

        sam_df = sam_to_df(sam_file_path)
        sam_final_df = add_mapped_end_points(sam_df)

        print('creating array containing counts per position')
        mapped_reads_per_pos = np.zeros(init_size, dtype=int)

        for read in sam_final_df.itertuples():
            start = read[4]-1 # -1 because SAM files are not 0-indexed
            end = read[5]

            mapped_reads_per_pos[start:end] = mapped_reads_per_pos[start:end] + 1

        mapped_reads_per_pos = mapped_reads_per_pos[0:ecoli_genome_len]

        sample_to_read_pos_dict.update({current_sample: mapped_reads_per_pos})

        print('finished sample!')
        print('_____________new sample____________________________')

    with open(output_path, 'wb') as handle:
        pickle.dump(sample_to_read_pos_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return(sample_to_read_pos_dict)

def load_count_dict(path_to_pickle):

    '''
    load in pickle file containing the sample_to_read_pos_dict dictionary
    (reads per position for each sample).this will only work if the peak drawing
    pipeline has been run completely before.

    PARAMETERS
    --------------------
    path_to_pickle: str
        path to pickle file containing the sample_to_read_pos_dict dictionary

    RETURNS
    --------------------
    sample_to_read_pos_dict_tmp: dictionary
        dictionary containing mapped reads per genome position for each sample

    '''

    with open(path_to_pickle, 'rb') as handle:
        sample_to_read_pos_dict_tmp = pickle.load(handle)

    return(sample_to_read_pos_dict_tmp)

def convert_sample_txt_to_list(path_to_sample_names):

    '''
    stores sample names from a text file into a list. this will be useful for
    indicating which samples to plot on the peak map

    PARAMETERS
    --------------------
    path_to_sample_names: str
        path to text file containing relevant sample names

    RETURNS
    --------------------
    samples_list: list of str
        list containing sample names

    '''
    samples_list = []
    with open(path_to_sample_names, 'r') as samples_file:
        for sample in samples_file:
            sample = sample.rstrip('\n')
            if '.deduplicated.mapped.sam' not in sample:
                sample = sample + '.deduplicated.mapped.sam'
            samples_list.append(sample)

    return(samples_list)

def find_potential_insertion_sites(sample_to_read_pos_dict, evolution_type,
                                   ecoli_genome_sequence):

    '''
    searches for TA site in a 11bp surrounding the position with the most
    read counts

    PARAMETERS
    --------------------
    sample_to_read_pos_dict_tmp: dictionary
        dictionary containing mapped reads per genome position for each sample
    evolution_type: str
        carbon source or condition in which cells were evolved [arbutin, serine]
    ecoli_genome_sequence:str
        reference sequence of E. coli genome reads were mapped onto


    RETURNS
    --------------------
    insert_sites_df: pandas dataframe
        dataframe with potential insert sites and flanking regions

    '''

    # define regions to check based on evolution type
    if evolution_type == 'arbutin':
        begin_check_pos = 3334000
        end_check_pos = 3336000
    elif evolution_type == 'serine':
        begin_check_pos = 1574500
        end_check_pos = 1579000

    insert_site_dict_list = []
    for sample in sample_to_read_pos_dict:

        # find location of peak with maximum intensity
        tmp_array = sample_to_read_pos_dict[sample][begin_check_pos:end_check_pos]
        m = max(tmp_array)
        max_index = [i + begin_check_pos for i, j in enumerate(tmp_array) if j == m]

        # find sequence and flanking sequences around max peak
        insert_site = max_index[0] #max_indices_per_sample[sample][0]
        insert_site_region = ecoli_genome_sequence[insert_site-1: insert_site+2]
        insert_site_upstream = ecoli_genome_sequence[insert_site-5: insert_site-1]
        insert_site_downstream = ecoli_genome_sequence[insert_site+2: insert_site+6]
        insert_site_all = insert_site_upstream + insert_site_region + insert_site_downstream
        insert_range_pos = str(insert_site-3) + ':' + str(insert_site+6)

        # check if TA site(s) is/are present in predicted insertion region
        if 'TA' in insert_site_all or 'AT' in insert_site_all:
            TA_present = 1
        else:
            TA_present = 0

        # store results in dictionary
        insert_site_dict_list.append({'sample': sample,
                                      'insert_site_upstream': insert_site_upstream,
                                     'insert_site_region': insert_site_region,
                                     'insert_site_downstream': insert_site_downstream,
                                     'TA_present': TA_present,
                                      'insert_range_pos': insert_range_pos
                                     })
    # convert results to dataframe for exporting
    insert_sites_df = pd.DataFrame(insert_site_dict_list)

    return(insert_sites_df)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--s', help='path to directory with mapped sam files')
    parser.add_argument('--p', help='path to pickle file with reads per position')
    parser.add_argument('--startpoints', help='file with sample names for startpoints')
    parser.add_argument('--endpoints', help='file with sample names for endpoints')
    parser.add_argument('-l', nargs='*', help='list of samples to plot')
    parser.add_argument('-e', help='evolution type (e.g., arbutin, serine)')
    parser.add_argument('-o', help='output figure name', required=True)

    args = parser.parse_args()

    path_to_mapped_sam_files=args.s
    path_to_pickle=args.p
    startpoints_path=args.startpoints
    endpoints_path=args.endpoints
    samples_to_extract=args.l
    evolution_type= args.e

    output_fig_name=args.o

    # define constant parameters
    ecoli_MDS42_genome_len = 3976195
    genome_positions = list(range(1, ecoli_MDS42_genome_len + 1))

    # check if user provided path to sam file (run whole pipeline) or path
    # to a pickle file (skip read end prediction and alignment counting)
    # note that presence of path to sam files takes precedence over pickle
    if path_to_mapped_sam_files is not None:
        output_dir = '/'.join(path_to_mapped_sam_files.split('/')[:-2]) + '/reads_per_position_dedup/'
    elif path_to_pickle is not None:
        output_dir = '/'.join(path_to_pickle.split('/')[:-2]) + '/reads_per_position_dedup/'
    else:
        sys.exit('Error! Need to supply either a valid path to sam files' \
         ' or valid path to pickle file. Aborting.')

    output_dir_figs = output_dir + 'figs'
    output_path = output_dir + 'mapped_reads_per_pos_arrays.pkl'

    # create ouput folders if needed
    CHECK_FOLDER = os.path.isdir(output_dir)
    CHECK_FOLDER_FIGS = os.path.isdir(output_dir + 'figs/')
    if not CHECK_FOLDER:
        os.makedirs(output_dir)
        print(f"created folder : {output_dir}")

    if not CHECK_FOLDER_FIGS:
        os.makedirs(output_dir_figs)
        print(f"created folder : {output_dir_figs}")

    if path_to_mapped_sam_files is not None:

        mapped_sam_files = glob.glob(path_to_mapped_sam_files + '*deduplicated.mapped.sam')
        mapped_sam_files.sort()

        sample_to_read_pos_dict = calc_reads_per_pos(mapped_sam_files, output_path)

    elif path_to_pickle is not None:
        print(f'Grabbing pickle: {path_to_pickle}')
        sample_to_read_pos_dict = load_count_dict(path_to_pickle)

    else:
         sys.exit('Error! Need to supply either a valid path to sam files' \
          ' or valid path to pickle file. Aborting.')

    # normalize alignment arrays using desired method
    # also, keep track of max value to help define y-axis range when plotting
    print('normalizing arrays')
    max_val = 0
    for sample in sample_to_read_pos_dict:
        tmp_array = sample_to_read_pos_dict[sample]
        tmp_array = normalize_array(tmp_array,normalization='fraction' )
        if np.max(tmp_array) > max_val:
            max_val = np.max(tmp_array)
        # update dictionary with normalized values. this new array is not saved
        # to disk
        sample_to_read_pos_dict.update({sample: tmp_array})

    startpoints_list = convert_sample_txt_to_list(startpoints_path)
    endpoints_list = convert_sample_txt_to_list(endpoints_path)

    if samples_to_extract is None:
        samples_to_extract = list(range(1, len(startpoints_list)+1))
    else:
        samples_to_extract = [int(i) for i in samples_to_extract]

    # extract whichever samples user wanted to plot.
    # NOTE: argument passed to '-l' is not zero-indexed.
    startpoints_to_plot_list = [startpoints_list[i-1]  for i in samples_to_extract]
    endpoints_to_plot_list = [endpoints_list[i-1]  for i in samples_to_extract]

    # do plotting
    fig, ax = plt.subplots(2,1,figsize = (6.0,4.0), dpi=400)
    colors = ['#000000', '#E877AE', '#426BB4','#1C4847', '#409091','#442D85',
              '#9375B5', '#851A19','#894E21', '#9A9A9B','#93C945']
    num_samples_to_plot = len(endpoints_to_plot_list)
    for end_point, start_point, color in zip(endpoints_to_plot_list, startpoints_to_plot_list, colors[0:num_samples_to_plot]):

        ax[0].plot(genome_positions, sample_to_read_pos_dict[start_point],
                    color=color, alpha=0.35)
        ax[0].spines['right'].set_visible(False)
        ax[0].spines['top'].set_visible(False)

        ax[1].plot(genome_positions, sample_to_read_pos_dict[end_point],
                    color=color, alpha=0.35)
        ax[1].spines['right'].set_visible(False)
        ax[1].spines['top'].set_visible(False)

    plt.xlabel('Genome position (bp)', fontsize=12)
    fig.text(-0.02, 0.5, 'Normalized peak intensity', va='center',
            rotation='vertical', fontsize=12)
    plt.xticks(np.arange(0, 4000001, 1000000))
    ax[0].set_ylim([-0.01, max_val])
    ax[1].set_ylim([-0.01, max_val])

    ax[0].tick_params(
    axis='x',
    which='both',
    bottom=False,
    top=False,
    labelbottom=False)

    fig.tight_layout()
    plt.savefig(output_dir_figs + '/' + output_fig_name + '.png')
    plt.savefig(output_dir_figs + '/' + output_fig_name + '.pdf')
    plt.savefig(output_dir_figs + '/' + output_fig_name + '.svg')

    fig, ax = plt.subplots(figsize = (6.0,4.0), dpi=400)

    # plot insets
    for end_point, color in zip(endpoints_to_plot_list, colors[0:num_samples_to_plot]):

        ax.plot(genome_positions, sample_to_read_pos_dict[end_point],
                color=color, alpha=0.35)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    plt.xlabel('Genome position (bp)', fontsize=12)
    fig.text(-0.02, 0.5, 'Normalized peak intensity', va='center',
             rotation='vertical', fontsize=12)

    if evolution_type == 'arbutin':
        plt.xticks(np.arange(3332000, 3336001, 1000))
        ax.set_ylim([-0.01, max_val])
        ax.set_xlim([3332900, 3336000])

        # bglG coordinates
        rect1 = patches.Rectangle((3334944, 0.02), 836, 0.04, linewidth=1,
                                    edgecolor='r', facecolor='none')
        ax.add_patch(rect1)

        # bglF coordinates
        rect2 = patches.Rectangle((3332933, 0.02), 1877, 0.04, linewidth=1,
                                    edgecolor='r', facecolor='none')
        ax.add_patch(rect2)

    elif evolution_type == 'serine':
        plt.xticks(np.arange(1574000, 1579001, 1000))
        ax.set_ylim([-0.01, max_val])
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
    elif evolution_type == 'gly-glu':
        # include coordinate for genes of interest of second carbon source
        pass
    else:
        # print check statement and return None to stop rest of code from running
        print(f'evolution type {evolution_type} not recognized. No inset will be plotted.')
        return

    fig.tight_layout()
    plt.savefig(output_dir_figs + '/' +output_fig_name + '_peak_inset.png')
    plt.savefig(output_dir_figs + '/' +output_fig_name + '_peak_inset.pdf')
    plt.savefig(output_dir_figs + '/' +output_fig_name + '_peak_inset.svg')

    # parse e coli genome so we can find regions flanking insertion site
    MDS42_genome = '../../tn-seq_data/NGS_ref_genomes/MDS42/AP012306.fa'
    MDS42_genome_parsable = SeqIO.parse(open(MDS42_genome),'fasta')

    for fasta in MDS42_genome_parsable:
        MDS42_genome_sequence = str(fasta.seq)

    insert_sites_df = find_potential_insertion_sites(sample_to_read_pos_dict,
                                                    evolution_type,
                                                    MDS42_genome_sequence)
    insert_sites_df.to_csv(output_dir + 'predicted_insertion_site.csv')

if __name__ == "__main__":
    main()
