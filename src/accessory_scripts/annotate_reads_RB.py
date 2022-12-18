#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on June 7 12:00:04 2022
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

    closest_feats_orig_df = pd.read_csv(closest_feats_orig_path, sep='\t', header=None)
    closest_feats_columns_orig = list(closest_feats_orig_df.columns)
    closest_feats_columns_new = ['read_chr', 'read_start', 'read_end', 'read_name',
                                      'feature_chr', 'database',
                                     'feature', 'feature_start', 'feature_end',
                                    'score', 'strand', 'frame', 'attribute','dist_read_feature']

    closest_feats_column_convert_dict = dict(zip(closest_feats_columns_orig,closest_feats_columns_new))
    closest_feats_final_df = closest_feats_orig_df.rename(columns=closest_feats_column_convert_dict)
    closest_feats_final_df = closest_feats_final_df[closest_feats_final_df['feature'] == 'gene']
    closest_feats_final_df_idxed = closest_feats_final_df.set_index('read_name')

    # read in peak file from macs3 output
    BC_positions_fname = BC_positions_path
    read_file_macs_df = pd.read_csv(BC_positions_fname).drop('Unnamed: 0', axis=1)

    reads_annotated_list = []

    ########## START DATAFRAME MERGING PROCESS ##########
    new_df = closest_feats_final_df_idxed.reset_index()[['read_name','attribute','dist_read_feature']]
    new_df['dist_read_feature'] = new_df['dist_read_feature'].astype('str')
    combined_new_df = new_df.groupby(['read_name']).agg({'attribute': '#'.join,'dist_read_feature':'#'.join})
    # combined_new_df.loc[read_file_macs_df['read_name_with_dir']] # validate successful merge
    ########## END DATAFRAME MERGING PROCESS ##########

    read_names_list = list(combined_new_df.index)
    attributes_list = list(combined_new_df['attribute'])
    dist_read_feature_list = list(combined_new_df['dist_read_feature'])
    read_name_to_attributes_dict = dict(zip(read_names_list, attributes_list))
    read_name_to_dist_read_dict = dict(zip(read_names_list, dist_read_feature_list))

    for read_attributes in tqdm.tqdm(read_file_macs_df.itertuples()):
        ## read_attributes tuple details
        # 0 - Index
        # 1 - read_name_bc
        # 2 - read_name_with_dir
        # 3 - barcode
        # 4 - flag
        # 5 - start
        # 6 - end
        # 7 - insert_size
        # 8 - chrom

        read_name = read_attributes[2]
        att_tmp = read_name_to_attributes_dict[read_name]
        dist_tmp = read_name_to_dist_read_dict[read_name]

        att_tmp_split = att_tmp.split('#')
        dist_tmp_split = dist_tmp.split('#')
        dist_tmp_split = list(map(int, dist_tmp_split))
        dist_tmp_split = list(map(abs, dist_tmp_split))

        dist_tmp_split_sorted, att_tmp_split_sorted = (list(t) for t in zip(*sorted(zip(dist_tmp_split, att_tmp_split))))

        dist_tmp_split_sorted = dist_tmp_split_sorted[:2]
        att_tmp_split_sorted = att_tmp_split_sorted[:2]

        gene_hits_tmp = [] # list for storing two closest peaks

        # grab gene names for genes most proximal to peaks
        for attribute_str in att_tmp_split_sorted:
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
                                   'start': read_attributes[5],
                                   'end': read_attributes[6],
                                    'barcode': read_attributes[3],
                                    'flag': read_attributes[4],
                                   'insert size': read_attributes[7],
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
    parser.add_argument('-i', help='path to _barcode_positions_final.csv')
    parser.add_argument('-c', help='path to _bedtools_closest_reads.bed')
    parser.add_argument('-o', help='output path prefix')
    args = parser.parse_args()

    BC_positions_final_path=args.i
    bedtools_closest_output_path=args.c
    output_path_prefix=args.o

    annotate_reads(BC_positions_final_path, bedtools_closest_output_path, output_path_prefix)

if __name__ == "__main__":
    main()
