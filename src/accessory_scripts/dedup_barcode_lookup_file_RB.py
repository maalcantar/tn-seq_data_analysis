#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on June 7 11:11:47 2022
@authors: alcantar and english
run in ngs_processing_source_fastp environment
"""
import argparse

import pandas as pd

def deduplicate_barcode_positions_lookup(deduplicated_fastq_fnames_path):

    '''
    uses deduplicated fastq read names to deduplicate barcode positions lookup file

    PARAMETERS
    --------------------
    deduplicated_fastq_fnames_path: str
        path to original barcode lookup file that was created by 01C script
    RETURNS
    --------------------
    NONE: just saves deduplicates barcode lookup file
    '''

    # read in file with deduplicated fastq read names
    sample_name = deduplicated_fastq_fnames_path.split('/')[-1].replace('_dedup_read_names.txt','')
    dedup_read_names_df = pd.read_csv(deduplicated_fastq_fnames_path, sep=" ", header=None)
    dedup_read_names_df.columns = ['read_name']

## code is useful if you want to parse down barcode .txt file directly
    barcode_file_txt_path = '/'.join(deduplicated_fastq_fnames_path.replace('fastq_dedup', 'seqkit_filtering_RB').split('/')[:-1]) + '/' + sample_name + '_barcodes.txt'
    barcode_data_df = pd.read_csv(barcode_file_txt_path, sep="\t", header=None)
    barcode_data_df = barcode_data_df.rename(columns = {0:'read_name', 6: 'barcode'})

    # read in non-deduplicated barcode lookup file
    barcode_positions_lookup_path = barcode_file_txt_path.replace('barcodes.txt', 'barcodes_positions_lookup.csv')
    barcode_positions_lookup_df = pd.read_csv(barcode_positions_lookup_path).drop('Unnamed: 0', axis=1)
    dedup_read_names = set(list(dedup_read_names_df['read_name']))

    # parse down file barcode lookup file so it only contains deduplicated reads
    barcode_positions_lookup_dedup_df = barcode_positions_lookup_df[barcode_positions_lookup_df['read_name_bc'].isin(dedup_read_names)].reset_index(drop=True)

    # save new deduplicated lookup file
    barcode_positions_lookup_dedup_out_path = barcode_positions_lookup_path.replace('seqkit_filtering_RB', 'seqkit_filtering_RB_dedup')
    barcode_positions_lookup_dedup_df.to_csv(barcode_positions_lookup_dedup_out_path)

    barcode_data_dedup_df = barcode_data_df[barcode_data_df['read_name'].isin(dedup_read_names)].reset_index(drop=True)
    barcode_data_dedup_out_path = barcode_file_txt_path.replace('seqkit_filtering_RB', 'seqkit_filtering_RB_dedup')
    barcode_data_dedup_df.to_csv(barcode_data_dedup_out_path, index=False, header=False, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='path to .txt file with deduplicated fastq read names')

    args = parser.parse_args()

    deduplicated_fastq_fnames_path=args.i

    deduplicate_barcode_positions_lookup(deduplicated_fastq_fnames_path)

if __name__ == "__main__":
    main()
