

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on May 27 12:13:45 2022
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

def add_mapped_end_points(barcode_pos_path):
    '''
    deduces end points of mapped reads from left-most mapped read position and insert size.
    appends end points to original dataframe

    PARAMETERS
    --------------------
    barcode_pos_path: str
        path to dataframe containing Rb-Tnseq information
        at the bare minimum, should contain: read_name_bc, barcode, flag, start, and insert_size
    RETURNS
    --------------------
    BC_positions_final_path: str
        path to finalized dataframe containing Rb-Tnseq information and read end points
    BC_positions_bed_final_path: str
        path to bedfile version of BC_positions_final_path. this is eventually
        fed into bedtools closest
    output_path_prefix: str
        prefix to output paths for all future files

    Additionally, saves the following two files:
    BC_positions_df: pandas dataframe
        dataframe containing Rb-Tnseq information and read end points
    BC_positions_bed_df: pandas dataframe
        same as BC_positions_df except in bed format
        col1: chromosome name; col2: start point, col3: end point
        col4: unique name

    '''

    BC_positions_df = pd.read_csv(barcode_pos_path)

    # obtain unique read names and initialize dictionary that maps
    # read names to end points
    read_names_unique = pd.unique(BC_positions_df['read_name_bc'])
    read_to_end_dict = dict()

    # loop through all unique read names
    for read_name in read_names_unique:

        # look at entries that match read name (should be 1 or 2, depending on whether one or
        # both reads in a pair mapped)
        tmp_df = BC_positions_df[BC_positions_df['read_name_bc'] == read_name]

        if len(tmp_df) == 2:
            # case where both reads map. use insert size and left-end position from
            # forward read to estimate end point
            tmp_df = tmp_df[tmp_df['insert_size'] > 0]
            end_site = tmp_df['start'] + tmp_df['insert_size']

        elif len(tmp_df) == 1:
            # case where only 1 read mapped. use only available left-end points
            # and insert size to estimate end site
            end_site = tmp_df['start'] + abs(tmp_df['insert_size'])

        else:
            # case where nothing maps. should not occur unless there's an error with
            # extracting unique read names, or if something mapped twice
            print('Error! No matching reads found or read mapped to >1 location')
        end_site = list(end_site)[0]
        read_to_end_dict.update({read_name: end_site})

    # initialize list with end points and read directions
    # read directions will help obtain unique read names for each read
    end_sites = []
    read_directions = []

    # loop through all reads
    for idx, read in BC_positions_df.iterrows():
        read_name = read['read_name_bc']
        insert_size = read['insert_size']
        barcode = read['read_name_bc']

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

    BC_positions_df = BC_positions_df[['read_name_bc', 'read_name_with_dir',
                                   'barcode', 'flag', 'start', 'end',
                                   'insert_size','chrom']]
    BC_positions_df = BC_positions_df.sort_values(by='start', ignore_index=True)
    BC_positions_bed_df = BC_positions_df[['chrom', 'start', 'end', 'read_name_with_dir']]

    current_sample = barcode_pos_path.split('/')[-1][:-30]
    output_path_prefix = '/'.join(barcode_pos_path.split('/')[:-2]) +'/RbTnseq_preprocessing/'+ current_sample

    BC_positions_final_path = output_path_prefix+'_barcode_positions_final.csv'
    BC_positions_bed_final_path = output_path_prefix+'_barcode_positions_bed_final.bed'

    BC_positions_df.to_csv(BC_positions_final_path)
    BC_positions_bed_df.to_csv(BC_positions_bed_final_path, header=False, index=False, sep='\t')

    return(BC_positions_final_path, BC_positions_bed_final_path, output_path_prefix)


def index_containing_substring(list_to_search, substring):
    '''
    adapted from https://stackoverflow.com/questions/2170900/get-first-list-index-containing-sub-string
    search a list for a specific substring, and return the index of the first
    element containing the substring.
    PARAMETERS
    --------------------
    list_to_search: list
        list of strings that will be searched for the presence of a specific substring.
    substring: str
        substring to search for.
    RETURNS
    --------------------
    index: int
        index indicating which list element contains the substring. only the first
        index is returned.
    '''

    for index, string in enumerate(list_to_search):
        if substring in string:
              return(index)
    return -1

def extract_gene_name(attributes_str, gene_id_prefix):
    '''
    searches attributes string of a genetic locus annotation and attempts to return the
    associated gene name
    PARAMETERS
    --------------------
    attributes_str: str
        string containing genetic locus annotation, with different attributes separated by
        the ';' character
    gene_id_prefix: str
        expected prefix before the gene name. most common prefixes are: ['gene=', 'Name=', 'protein_id=']
    RETURNS
    --------------------
    gene_name: str
        name of gene (e.g., rpoB)
    '''
    attribute_str_split = attributes_str.split(';')
    gene_idx = index_containing_substring(attribute_str_split,gene_id_prefix) # find index using helper function
    gene_name = attribute_str_split[gene_idx].split('=')[-1]

    return(gene_name)


def annotate_reads(BC_positions_path, closest_feats_orig_path, output_path_prefix):

    '''
    annotate reads with the two most proximal genes

    PARAMETERS
    --------------------
    BC_positions_path: str
        path to dataframe containing finalized Rb-Tnseq information. this path
        should be output by add_mapped_end_points
    closest_feats_orig_path: str
        bedfile containing closest gene features to each read
    output_path_prefix: str
        output path prefix to use

    RETURNS
    --------------------
    NONE: just creates a csv file with gene features closest to each read

    '''

    # read in closest peaks file and tidy column names
    # this closest peaks file indicates which genomic features a specific peak is nearest
    # if not closest_feats_orig_path:
    #     closest_feats_orig_path = BC_positions_path[:-16].replace('macs3_output','bedtools_output') + '_bedtools_closest_peaks.bed'

    while not os.path.exists(closest_feats_orig_path):
        time.sleep(1)

    closest_feats_orig_df = pd.read_csv(closest_feats_orig_path, sep='\t', header=None)
    closest_feats_columns_orig = list(closest_feats_orig_df.columns)
    closest_feats_columns_new = ['read_chr', 'read_start', 'read_end', 'read_name',
                                      'feature_chr', 'database',
                                     'feature', 'feature_start', 'feature_end',
                                    'score', 'strand', 'frame', 'attribute','dist_read_feature']

    closest_feats_column_convert_dict = dict(zip(closest_feats_columns_orig,closest_feats_columns_new))
    closest_feats_final_df = closest_feats_orig_df.rename(columns=closest_feats_column_convert_dict)

    # read in peak file from macs3 output
    BC_positions_fname = BC_positions_path #r"01_S95_L002_macs3_peaks.xls"
    read_file_macs_df = pd.read_csv(BC_positions_fname).drop('Unnamed: 0', axis=1)

    reads_annotated_list = []
    for peak, read_attributes in read_file_macs_df.iterrows():

        # find all genes that are 'near' peak
        read_name = read_attributes['read_name_with_dir']
        closest_tmp_df = closest_feats_final_df.copy()[\
                               (closest_feats_final_df['read_name'] == read_name) &\
                               (closest_feats_final_df['feature'] == 'gene')]

        # sort genes by proximity to peak, and take the two most proximal genes
        closest_tmp_df['dist_read_feature'] = closest_tmp_df['dist_read_feature'].abs()
        closest_tmp_sort_df = closest_tmp_df.sort_values(by='dist_read_feature', ascending=False)
        closest_tmp_sort_top_two_df = closest_tmp_sort_df[:2]

        gene_hits_tmp = [] # list for storing two closest peaks

        # grab gene names for genes most proximal to peaks
        for read_tmp, read_attributes_tmp in closest_tmp_sort_top_two_df.iterrows():
            attribute_str = read_attributes_tmp['attribute']
            if 'gene=' in attribute_str:
                gene_hits_tmp.append(extract_gene_name(attribute_str,'gene='))
            elif 'Name=' in attribute_str:
                gene_hits_tmp.append(extract_gene_name(attribute_str,'Name='))
            elif 'protein_id=' in attribute_str:
                gene_hits_tmp.append(extract_gene_name(attribute_str,'protein_id='))
            else:
                gene_hits_tmp.append('NO_GENE')

        # use empty string to indicate no proximal gene(s)
        if len(gene_hits_tmp) < 2:
            gene_hits_tmp.extend(['']*(2-len(gene_hits_tmp)))

        # storage peak annotations in list
        reads_annotated_list.append({'read_name': read_name,
                               'start': read_attributes['start'],
                               'end': read_attributes['end'],
                                'barcode': read_attributes['barcode'],
                                'flag': read_attributes['flag'],
                               'insert size': read_attributes['insert_size'],
                               'gene1': gene_hits_tmp[0],
                               'gene2': gene_hits_tmp[1]})

    # peak annotations in dataframe
    reads_annotated_df = pd.DataFrame(reads_annotated_list)#.sort_values(by='fold-enrichment', ascending=False)
    reads_annotated_df['gene1_barcode'] = reads_annotated_df['gene1'] + '_' + reads_annotated_df['barcode']
    reads_annotated_df['gene2_barcode'] = reads_annotated_df['gene2'] + '_' + reads_annotated_df['barcode']

    reads_annotated_output_path = output_path_prefix + '_annotated_reads.csv'
    reads_annotated_df.to_csv(reads_annotated_output_path)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='directory containing barcode positions csv files')
    args = parser.parse_args()

    path_to_barcode_pos=args.i
    barcodes_positions_files = glob.glob(path_to_barcode_pos + '*_barcodes_positions_lookup.csv')

    export_dir = '/'.join(barcodes_positions_files[0].split('/')[:-2]) + '/RbTnseq_preprocessing'
    
    # adapted from https://djangocentral.com/check-if-a-directory-exists-if-not-create-it/
    CHECK_FOLDER = os.path.isdir(export_dir)
    if not CHECK_FOLDER:
        os.makedirs(export_dir)
        print("created folder : ", export_dir)

    for sample_file in barcodes_positions_files:
        print('Processing sample: ' + sample_file)
        BC_positions_final_path, BC_positions_bed_final_path, output_path_prefix = add_mapped_end_points(sample_file)
        bedtools_closest_output_path = output_path_prefix + '_bedtools_closest_reads.bed'
        Process=subprocess.Popen('bedtools closest -D a -k 4 -t'\
        ' all -a %s -b ../../tn-seq_data/NGS_ref_genomes/MDS42/AP012306_sorted.gff3 >'\
        '%s' % (BC_positions_bed_final_path,bedtools_closest_output_path,), shell=True)

        annotate_reads(BC_positions_final_path, bedtools_closest_output_path, output_path_prefix)

main()
