#!/bin/bash

######################################################

# This script performs the following:
# 1) clusters barcodes using Bartender (https://github.com/LaoZZZZZ/bartender-1.1)
# 2) assigns clustered barcodes to original reads

### PARAMETERS
## $1: list of .txt files with barcodes extracted from SEQKIT (scipt 01C)
##
# Example: sh 07C_cluster_barcodes.sh \
# ../../tn-seq_data/tn-seq_outputs/mae-009/seqkit_filtering/*_barcodes.txt

### OUTPUTS
# ${barcode_cluster_out_path}_baracode_clustering_output_barcode.csv: unique read to cluster id
# ${barcode_cluster_out_path}_baracode_clustering_output_cluster.csv: cluster id to centered barcode
# ${barcode_cluster_out_path}_baracode_clustering_output_quality.csv: base at each barcode position
# ${barcode_cluster_out_path}_baracode_read_mapping.csv: barcode to read name mapping
# ${barcode_cluster_out_path}_clustered_barcode_annotated_reads.csv: updated read annotations with clustered barcodes
######################################################

# Activate virtual enviroment
source activate ngs_processing_source_fastp

for barcode_file in $@
  do
    fname=$(basename $barcode_file)
  	dname=$(dirname $barcode_file)
  	fpath=$dname/${fname%S*}

    barcode_cluster_out_path=${dname%seqkit_filtering*}barcode_clustering/${fname%_barcodes.txt*}

    if [[ "$barcode_file" == "$1" ]]
      then
        if [ ! -d "${dname%seqkit_filtering*}barcode_clustering/" ]
          then
            echo "making barcode clustering directory"
    		    mkdir ${dname%seqkit_filtering*}barcode_clustering/
	      fi
    fi

    echo "Beginning $fname"

    awk '{print $7 "," $1}' $barcode_file > ${barcode_cluster_out_path}_baracode_read_mapping.csv

    bartender_single_com -f ${barcode_cluster_out_path}_baracode_read_mapping.csv\
     -o ${barcode_cluster_out_path}_baracode_clustering_output -d 2

     echo ${barcode_cluster_out_path}

     python accessory_scripts/update_barcode_cluster_annotations.py -i ${barcode_cluster_out_path}_

done
