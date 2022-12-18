
"""
Created on May 26 17:08:00 2022
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
accessory_scripts % python3 map_barcodes_positions.py -barcode ../BC-1_S1_L001_barcodes.txt -align ../BC-4_S4_L001_paired_mapped.sam
'''

def barcode_to_position_mapper(barcode_file_txt, mapped_reads_sam):
    
    
    barcode_data = pd.read_csv(barcode_file_txt, sep="\t", header=None)
    barcode_data = barcode_data[[0,6]]
    barcode_data.columns = ['read_name', 'barcode']

    ### Loop the data lines
    with open(mapped_reads_sam, 'r') as temp_f:
        # get No of columns in each line
        col_count = [ len(l.split("\t")) for l in temp_f.readlines() ]

    ### Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
    column_names = [i for i in range(0, max(col_count))]

    ### Read csv
    alignment_data = pd.read_csv(mapped_reads_sam, sep="\t", header=None,\
                                 names=column_names)
    
    alignment_data = alignment_data[[0,1,3,8,9]]
    alignment_data.columns = ['read_name', 'flag','start','insert_size','sequence']


    bc_location_cross_ref = \
    pd.merge(barcode_data, alignment_data, \
                 how='inner', left_on=barcode_data['read_name'], right_on=alignment_data['read_name'], sort=False, suffixes=('_bc', '_bowtie'), \
                 copy=True, indicator=True)
    
    bc_location_cross_ref = bc_location_cross_ref[['read_name_bc','barcode','flag','start','insert_size']]
    
    
    output_path = barcode_file_txt[:-12]+'barcodes_positions_lookup.csv'
    bc_location_cross_ref.to_csv(output_path, sep=',')
    
    return bc_location_cross_ref

    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-barcode', help='path to the .txt file with the barcode sequences')
    parser.add_argument('-align', help='path to the filtered, aligned SAM file')

    args = parser.parse_args()

    barcode_file_txt=args.barcode
    mapped_reads_sam=args.align

    barcode_to_position_mapper(barcode_file_txt, mapped_reads_sam)    

main()
