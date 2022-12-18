
"""
Created on Mar 4 15:28:00 2022
@authors: alcantar and english
"""

import argparse
import pandas as pd
import numpy as np
import glob
import os
import sys

'''
example usage (terminal):
accessory_scripts % python3 make_PGT_bedfile.py -b ../../../tn-seq_data/NGS_ref_genomes/MG1655/U00096_sorted.bed
'''

# Names of the columns for a bed9 format file + the attributes column from a gff3 file
bedfile_heads = ['chrom',\
                'chromStart',\
                'chromEnd',\
                'name',\
                'score',\
                'strand',\
                'thickStart',\
                'thickEnd',\
                'itemRgb',\
                'attributes']

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

def extract_gene_name(attribute_str_split, gene_id_prefix):
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
    gene_idx = index_containing_substring(attribute_str_split,gene_id_prefix) # find index using helper function
    gene_name = attribute_str_split[gene_idx].split('=')[-1]

    return(gene_name)

# Step 1: convert scores to zero
# Convert "." to zero in score to 0
def zeroscore_genomefile(filepath):

    '''
    Takes a bed format file and replaces the score value of '.' with '0' in order for pyGenomTracks to read past the 4th column without throwing an error

    PARAMETERS
    --------------------
    filepath: path
        path to the sorted reference genome bedfile.

    RETURNS
    --------------------
    filepath_zeroscore.bed: .bed file
        a tab separated text file with a .bed suffix, equivalent to the input sorted bedfile but the score column has been replaced with zeros.
    '''
    # Read in the bed file
    bedfile = pd.read_csv(filepath, sep='\t',header=None,names=['chrom','chromStart','chromEnd','name','score','strand','source','type','phase','attributes'], \
                   skiprows=0)
    bedfile_zeroscore = bedfile.copy(deep=True)
    index = 0

    # iterate through and replace the '.' scores with 0
    for row in bedfile['score']:
        if bedfile.loc[index,'score'] == '.':
            bedfile_zeroscore.loc[index,'score'] = 0
        index +=1

    bedfile_zeroscore.to_csv(filepath[:-4]+'_zeroscore.bed', sep='\t', header=False, index = False)
    return bedfile_zeroscore



def create_PGT_bedfile(filepath,zeroscore_bedfile):

    '''
    Takes a zero-scored bed format file and performs the following:
        1. Makes the thickStart and thickEnd columns equivalent to the chromStart and chromEnd columns
        2. Gets a human-readable name entry for plotting in pyGenomTracks based on the different types of gbkeys (add if needed)
        3. Adds a user-defined color to the itemRgb column for plotting on pyGenomTracks
        4. Removes the final attribute column to make this a correct bed9 file format

    PARAMETERS
    --------------------
    filepath: path
        path to the sorted reference genome bedfile.
    zeroscore_bedfile: dataframe
        bedfile dataframe from the function zeroscore_genomefile, which is the genome bed file but with the score column set to zero

    RETURNS
    --------------------
    filepath_PGT_format.bed: .bed file
        a tab separated text file with a .bed suffix, for use with pyGenomTracks plotting
    '''

    outfilepath = filepath[:-10]
        # Find all the unique types of gbkey
    zeroscore_bed_copy = zeroscore_bedfile.copy()
    zeroscore_bed_copy.columns = bedfile_heads

    # Replace the gff3 format columns with the thickStart and thickEnd
    # To do this, uses chromStart and chromEnd as default
    zeroscore_bed_copy['thickStart'] = zeroscore_bed_copy['chromStart']
    zeroscore_bed_copy['thickEnd'] = zeroscore_bed_copy['chromEnd']

    #colors =
    #rgb(185,144,149)", <-> #b99095
    #"rgb(252,181,172)",<-> #fcb5ac
    #"rgb(181,229,207)",<-> #b5e5cf
    #"rgb(61,91,89)",<-> #3d5b59
    #"rgb(41,195,234)",<-> #29c3ea
    #"rgb(254,217,108)",<-> #fed96c
    #"rgb(236,220,235)",<-> #ecdceb
    #"rgb(214,214,214)", <-> #d6d6d6

    ## Iterate through the bedfile and replace the name with a human-readable one for plotting
    ## To do this, uses a series of user-defined cases based on the unique gbkeys found in the genome

    # For each gbkey
    # Get a name based on the attributes column of the gff3 files
    # Assign a colour to the itemRGB column that is user defined

    for x in range(len(zeroscore_bed_copy['attributes'])):
        attributes_str_split = zeroscore_bed_copy['attributes'][x].split(';')
        for attribute in attributes_str_split:

            # Handle CDS features
            if 'gbkey=CDS' in attribute:
                #print('found cds')
                gene_id_prefix = 'Name='
                name = extract_gene_name(attributes_str_split, gene_id_prefix)
                color = '181,229,207'
                zeroscore_bed_copy.loc[x, 'name'] = name
                zeroscore_bed_copy.loc[x, 'itemRgb'] = color

            # Handle gene features
            elif 'gbkey=Gene' in attribute:
                #print('found gene')
                gene_id_prefix = 'gene='
                name = extract_gene_name(attributes_str_split, gene_id_prefix)
                color = '252,181,172'
                zeroscore_bed_copy.loc[x, 'name'] = name
                zeroscore_bed_copy.loc[x, 'itemRgb'] = color

            # Handle genome features
            elif 'gbkey=Src' in attribute:
                #print('found genome')
                gene_id_prefix = 'substrain='
                name = extract_gene_name(attributes_str_split, gene_id_prefix)
                color = '214,214,214'
                zeroscore_bed_copy.loc[x, 'name'] = name
                zeroscore_bed_copy.loc[x, 'itemRgb'] = color

            # Handle misc_RNA features
            elif 'gbkey=misc_RNA' in attribute:
                #print('found misc_RNA')
                gene_id_prefix =  'locus_tag='
                name = extract_gene_name(attributes_str_split, gene_id_prefix)
                color = '185,144,149'
                zeroscore_bed_copy.loc[x, 'name'] = name
                zeroscore_bed_copy.loc[x, 'itemRgb'] = color

            # Handle misc_feature features
            elif 'gbkey=misc_feature' in attribute:
                #print('found misc_feature')
                gene_id_prefix = 'Note='
                name = extract_gene_name(attributes_str_split, gene_id_prefix)
                color = '214,214,214'
                zeroscore_bed_copy.loc[x, 'name'] = name
                zeroscore_bed_copy.loc[x, 'itemRgb'] = color

            # Handle ncRNA features
            elif 'gbkey=ncRNA' in attribute:
                #print('found ncRNA')
                gene_id_prefix = 'locus_tag='
                name = extract_gene_name(attributes_str_split, gene_id_prefix)
                color = '61,91,89'
                zeroscore_bed_copy.loc[x, 'name'] = name
                zeroscore_bed_copy.loc[x, 'itemRgb'] = color

            # Handle rRNA features
            elif 'gbkey=rRNA' in attribute:
                #print('found rRNA')
                gene_id_prefix = 'locus_tag='
                name = extract_gene_name(attributes_str_split, gene_id_prefix)
                color = '185,144,149'
                zeroscore_bed_copy.loc[x, 'name'] = name
                zeroscore_bed_copy.loc[x, 'itemRgb'] = color

            # Handle tRNA features
            elif 'gbkey=tRNA' in attribute:
                #print('found tRNA')
                gene_id_prefix = 'locus_tag='
                name = extract_gene_name(attributes_str_split, gene_id_prefix)
                color = '236,220,235'
                zeroscore_bed_copy.loc[x, 'name'] = name
                zeroscore_bed_copy.loc[x, 'itemRgb'] = color

            # Handle tmRNA features
            elif 'gbkey=tmRNA' in attribute:
                #print('found tmRNA')
                gene_id_prefix = 'locus_tag='
                name = extract_gene_name(attributes_str_split, gene_id_prefix)
                color = '236,220,235'
                zeroscore_bed_copy.loc[x, 'name'] = name
                zeroscore_bed_copy.loc[x, 'itemRgb'] = color

            # Handle mobile_element features
            elif 'gbkey=mobile_element' in attribute:
                #print('found mobile_element')
                gene_id_prefix = 'mobile_element_type='
                name = extract_gene_name(attributes_str_split, gene_id_prefix)
                color = '41,195,234'
                zeroscore_bed_copy.loc[x, 'name'] = name
                zeroscore_bed_copy.loc[x, 'itemRgb'] = color

            # Handle rep_origin features
            elif 'gbkey=rep_origin' in attribute:
                #print('found rep_origin')
                gene_id_prefix = 'Note='
                name = extract_gene_name(attributes_str_split, gene_id_prefix)
                color = '254,217,108'
                zeroscore_bed_copy.loc[x, 'name'] = name
                zeroscore_bed_copy.loc[x, 'itemRgb'] = color

    # Remove column 10 to make this a bed9 file format
    zeroscore_bed_copy.drop(labels='attributes', axis=1,inplace=True)
    # Remove the genome feature to make the pyGenomTracks plots look cleaner
    zeroscore_bed_copy.drop(labels=0, axis=0,inplace=True)
    zeroscore_bed_copy.to_csv(outfilepath+'PGT_format.bed', sep='\t', header=False, index = False)
    return zeroscore_bed_copy

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', help='path to the sorted genome bedfile')
    args = parser.parse_args()

    path_to_bedfile=args.b

    bedfile_zeroscore = zeroscore_genomefile(path_to_bedfile)
    final_pgt_bedfile = create_PGT_bedfile(path_to_bedfile,bedfile_zeroscore)

main()
