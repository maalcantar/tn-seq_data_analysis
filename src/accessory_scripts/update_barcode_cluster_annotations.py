#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on May 30 14:17:22 2022
@authors: alcantar and english
"""
import argparse
import pandas as pd

def annotate_reads_with_clustered_bc(barcode_clustering_path):
    '''
    Update the annotated reads dataframe with clustered barcodes

    PARAMETERS
    --------------------
    barcode_clustering_path: str
        path prefix for barcode data for current sample

    RETURNS
    --------------------
    NONE: just creates a csv file with updated annotated reads

    '''

    # read in dataframes with barcode cluster information
    unique_bc_to_clusterid_path = barcode_clustering_path + 'baracode_clustering_output_barcode.csv'
    clusterid_to_seed_bc_path = barcode_clustering_path + 'baracode_clustering_output_cluster.csv'
    unique_bc_to_clusterid_df = pd.read_csv(unique_bc_to_clusterid_path, low_memory=False)
    clusterid_to_seed_bc_df = pd.read_csv(clusterid_to_seed_bc_path,low_memory=False)

    # make a dictionary with unique barcodes are keys and cluster ids as values
    unique_bcs = unique_bc_to_clusterid_df['Unique.reads']
    cluster_ids = unique_bc_to_clusterid_df['Cluster.ID']
    unique_bc_to_cluster_id_dict = dict(zip(unique_bcs, cluster_ids))

    # make a dictionary with cluster ids as keys and barcode centers (i.e., true barcodes) as values
    cluster_ids_centers = clusterid_to_seed_bc_df['Cluster.ID']
    bc_centers = clusterid_to_seed_bc_df['Center']
    cluster_id_to_centers_dict = dict(zip(cluster_ids_centers, bc_centers))

    # recreate path to RbTnseq_preprocessing data -- which contains the original annotated reads dataframe
    RbTnseq_preprocessing_path = barcode_clustering_path.replace('barcode_clustering', 'RbTn-seq_preprocessing') + 'annotated_reads.csv'
    annotated_reads_df = pd.read_csv(RbTnseq_preprocessing_path, index_col=0, low_memory=False)

    # initialize lists that will contain cluster ids and centers for each read
    cluster_ids_list = []
    center_bc_list = []

    # loop through all reads and assign cluster id and barcode center
    for barcode in annotated_reads_df['barcode']:
        cluster_id_tmp = unique_bc_to_cluster_id_dict[barcode]
        center_bc_tmp = cluster_id_to_centers_dict[cluster_id_tmp]
        cluster_ids_list.append(cluster_id_tmp)
        center_bc_list.append(center_bc_tmp)

    annotated_reads_df['cluster_id'] = cluster_ids_list
    annotated_reads_df['clustered_barcode'] = center_bc_list

    # create gene-bc id mapping using clustered barcodes
    annotated_reads_df['gene1_clust_id'] = annotated_reads_df['gene1'] + '.' + annotated_reads_df['cluster_id'].apply(str)
    annotated_reads_df['gene2_clust_id'] = annotated_reads_df['gene2'] + '.' + annotated_reads_df['cluster_id'].apply(str)
    # create gene-bc sequence mapping using clustered barcodes
    annotated_reads_df['gene1_clust_bc'] = annotated_reads_df['gene1'] + '.' + annotated_reads_df['clustered_barcode'].apply(str)
    annotated_reads_df['gene2_clust_bc'] = annotated_reads_df['gene2'] + '.' + annotated_reads_df['clustered_barcode'].apply(str)

    # only return barcodes that match the expected length
    annotated_reads_df = annotated_reads_df[annotated_reads_df['clustered_barcode'].apply(lambda x: len(x) == 20)]

    # filter for barcodes above a defined read count (helps remove spurious barcodes)
    # should be changed based on approximate number of starting colonies
    annotated_reads_df = annotated_reads_df[annotated_reads_df.groupby('clustered_barcode').clustered_barcode.transform('count') >= 1000].reset_index(drop=True)

    # create path to output final dataframe
    output_path = barcode_clustering_path + 'clustered_barcode_annotated_reads.csv'

    annotated_reads_df.to_csv(output_path)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='path prefix for barcode data for current sample')
    args = parser.parse_args()

    barcode_clustering_path=args.i

    annotate_reads_with_clustered_bc(barcode_clustering_path)

if __name__ == "__main__":
    main()
