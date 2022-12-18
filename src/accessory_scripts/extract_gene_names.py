#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mar 1 09:56:32 2022
@authors: alcantar and english
"""

import argparse
import pandas as pd

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

def extract_gene_names(bed_dir):
    '''
    extracts all gene names and corresponding locus tags from a bedfile.

    PARAMETERS
    --------------------
    bed_dir: str
        location of bedfile you want to parse

    RETURNS
    --------------------
    gene_annotation_df: pandas dataframe
        dataframe with gene SYMBOLS (first column) and corresponding locus tags (second column)
    '''

    # give headers descriptive names
    bedtools_file_MDS42_df = pd.read_csv(bed_dir, sep='\t', header=None )
    new_col_names = {0:'genome_name',
                    1: 'start',
                    2: 'end',
                    3: 'locus_tag',
                    4: 'score',
                    5: 'strand',
                    6: 'database',
                    7: 'feature_annotation',
                    8: 'score_2',
                    9: 'attributes'}
    bedtools_file_MDS42_df = bedtools_file_MDS42_df.rename(columns=new_col_names)

    # only look at genes
    bedtools_file_MDS42_df = bedtools_file_MDS42_df[bedtools_file_MDS42_df['feature_annotation'] == 'gene'].reset_index(drop=True)

    # initialize lisds for storing gene SYMBOLS and locus tags
    gene_name_dict_list = []
    gene_hits = []
    locus_tags = []
    gene_start = []
    gene_end = []

    # extract all gene SYMBOLS and locus tags
    for gene, gene_attributes_tmp in bedtools_file_MDS42_df.iterrows():
                attribute_str = gene_attributes_tmp['attributes']
                gene_start = gene_attributes_tmp['start']
                gene_end = gene_attributes_tmp['end']

                if 'gene=' in attribute_str:
                    gene_hits.append(extract_gene_name(attribute_str,'gene='))
                elif 'Name=' in attribute_str:
                    gene_hits.append(extract_gene_name(attribute_str,'Name='))
                elif 'protein_id=' in attribute_str:
                    gene_hits.append(extract_gene_name(attribute_str,'protein_id='))
                else:
                    gene_hits.append('NO_GENE')

                if 'locus_tag=' in attribute_str:
                    locus_tags.append(extract_gene_name(attribute_str,'locus_tag='))
                elif 'Parent=' in attribute_str:
                    locus_tags.append(extract_gene_name(attribute_str,'Parent='))
                else:
                    locus_tags.append('NO_LOCUS_TAG')

                gene_name_dict_list.append({'gene_name': gene_hits[-1],
                                           'locus_tag': locus_tags[-1],
                                           'start': gene_start,
                                           'end': gene_end})

    # create dataframe
    gene_annotation_df = pd.DataFrame(gene_name_dict_list)

    # output csv file with gene SYMBOLS and locus tags
    if 'AP012306' in bed_dir.split('/')[-1]:
        out_dir = '/'.join(bed_dir.split('/')[:-1]) + '/AP012306_gene_list.csv'
    elif 'U00096' in bed_dir.split('/')[-1]:
        out_dir = '/'.join(bed_dir.split('/')[:-1]) + '/U00096_gene_list.csv'
    else:
        dir_loc_tmp = '/' + bed_dir.split('/')[-1].split('_')[0] + '_gene_list.csv'
        out_dir = '/'.join(bed_dir.split('/')[:-1]) + dir_loc_tmp

    gene_annotation_df.to_csv(out_dir)

    return(gene_annotation_df)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='path to reference genome bedfile')
    args = parser.parse_args()

    bed_dir=args.i

    gene_annotation_df = extract_gene_names(bed_dir)

if __name__ == "__main__":
    main()
