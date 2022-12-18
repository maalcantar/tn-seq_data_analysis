#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Feb 21 18:20:44 2022
@authors: alcantar and english
"""
import argparse

import pandas as pd
import numpy as np
import glob
import os
import sys


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

def annotate_peaks(peak_file_macs_path, export_name, closest_peak_feats_orig_path=None):

    '''
    annotate peaks with the two most proximal genes

    PARAMETERS
    --------------------
    peak_file_macs_path: str
        path to file containing peak location and attributes (e.g. fold-enrichment and q-value).
        this is output by the 'macs3 callpeak' function in script '01_fastp_align_fastq.sh'
    export_name: str
        name of exported csv with peak annotations
    closest_peak_feats_orig_path: str (optional)
        path to file containing peaks most proximal to genomic features.
        this is output by the 'bedtools closest' function run in script '01_fastp_align_fastq.sh'
        if filepath is not specified, it is interpreted from peak_file_macs_path

    RETURNS
    --------------------
    peaks_annotated_df: pandas dataframe
        dataframe with peak annotations: peak id, peak start, peak end, peak length,
        fold-enrichment, -log10(pvalue), -log10(qvalue), most proximal gene, and
        second most proximal gene
    '''

    # read in closest peaks file and tidy column names
    # this closest peaks file indicates which genomic features a specific peak is nearest
    if not closest_peak_feats_orig_path:
        closest_peak_feats_orig_path = peak_file_macs_path[:-16].replace('macs3_output','bedtools_output') + '_bedtools_closest_peaks.bed'
    closest_peak_feats_orig_df = pd.read_csv(closest_peak_feats_orig_path, sep='\t', header=None)
    closest_peak_feats_columns_orig = list(closest_peak_feats_orig_df.columns)
    closest_peak_feats_columns_new = ['peak_chr', 'peak_summit_start', 'peak_summit_end', 'peak_id',
                                     'log_q_val', 'feature_chr', 'database',
                                     'feature', 'feature_start', 'feature_end',
                                    'score', 'strand', 'frame', 'attribute','dist_peak_feature']

    closest_peak_feats_column_convert_dict = dict(zip(closest_peak_feats_columns_orig,closest_peak_feats_columns_new))
    closest_peak_feats_final_df = closest_peak_feats_orig_df.rename(columns=closest_peak_feats_column_convert_dict)

    # read in peak file from macs3 output
    peak_file_macs_fname = peak_file_macs_path #r"01_S95_L002_macs3_peaks.xls"
    peak_file_macs_df = pd.read_table(peak_file_macs_fname,skiprows=24, header=0)

    peaks_annotated_list = []
    for peak, peak_attributes in peak_file_macs_df.iterrows():

        # find all genes that are 'near' peak
        peak_id = peak_attributes['name']
        closest_tmp_df = closest_peak_feats_final_df.copy()[\
                               (closest_peak_feats_final_df['peak_id'] == peak_id) &\
                               (closest_peak_feats_final_df['feature'] == 'gene')]

        # sort genes by proximity to peak, and take the two most proximal genes
        closest_tmp_df['dist_peak_feature'] = closest_tmp_df['dist_peak_feature'].abs()
        closest_tmp_sort_df = closest_tmp_df.sort_values(by='dist_peak_feature', ascending=False)
        closest_tmp_sort_top_two_df = closest_tmp_sort_df[:2]

        gene_hits_tmp = [] # list for storing two closest peaks

        # grab gene names for genes most proximal to peaks
        for peak_tmp, peak_attributes_tmp in closest_tmp_sort_top_two_df.iterrows():
            attribute_str = peak_attributes_tmp['attribute']
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
        peaks_annotated_list.append({'peak_id': 'peak' + str(peak),
                               'start': peak_attributes['start'],
                               'end': peak_attributes['end'],
                               'length': peak_attributes['length'],
                               'fold-enrichment': peak_attributes['fold_enrichment'],
                               '-log10(pvalue)': peak_attributes['-log10(pvalue)'],
                               '-log10(qvalue)': peak_attributes['-log10(qvalue)'],
                               'gene1': gene_hits_tmp[0],
                               'gene2': gene_hits_tmp[1]})

    # peak annotations in dataframe
    peaks_annotated_df = pd.DataFrame(peaks_annotated_list)#.sort_values(by='fold-enrichment', ascending=False)
    peaks_annotated_df.to_csv(export_name)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='directory containing macs3 peak data')
    args = parser.parse_args()

    path_to_macs3=args.i
    peak_files = glob.glob(path_to_macs3 + '*macs3_peaks.xls')

    for sample_file in peak_files:
        peak_file_macs_path = sample_file
        basename = peak_file_macs_path[:-16]
        closest_peak_feats_orig_path = basename.replace('macs3_output','bedtools_output') + \
        '_bedtools_closest_peaks.bed'


        export_name = basename.replace('macs3_output','annotated_peaks') + '_annotated_peak.csv'
        export_name_prefix = '/'.join(basename.replace('macs3_output','annotated_peaks').split('/')[:-1])+'/'

        # adapted from https://djangocentral.com/check-if-a-directory-exists-if-not-create-it/
        CHECK_FOLDER = os.path.isdir(export_name_prefix)
        if not CHECK_FOLDER:
            os.makedirs(export_name_prefix)
            print("created folder : ", export_name_prefix)

        annotate_peaks(peak_file_macs_path, export_name)
        
if __name__ == "__main__":
    main()
