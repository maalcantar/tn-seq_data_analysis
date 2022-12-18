#!/bin/bash

######################################################

# This script performs the following:
# 1) adds end points for all mapped reads
# 2) runs bedtools closest function
# 3) does peak annotation

### PARAMETERS
## $1: list of barcode lookup files

# Example: sh 05C_peak_annotation.sh \
# ../../tn-seq_data/tn-seq_outputs/mae-010/seqkit_filtering_RB/*_barcodes_positions_lookup.csv

### OUTPUTS
# ${barcode_cluster_out_path}_barcode_positions_final.csv: csv with mapping and
# barcode information
# ${barcode_cluster_out_path}_barcode_positions_bed_final.bed: same as above but
#  in .bed format -- compatible with bectools functions
# ${barcode_cluster_out_path}_bedtools_closest_reads.bed: genes closest to each read
# ${barcode_cluster_out_path}_annotated_reads.csv: reads annotated with closest genes
# and barcode information

######################################################

# Activate virtual enviroment
source activate ngs_processing_source_fastp

for barcode_lookup in $@
  do
    fname=$(basename $barcode_lookup)
  	dname=$(dirname $barcode_lookup)
  	fpath=$dname/${fname%S*}


    barcode_cluster_out_path=${dname%seqkit_filtering*}RbTn-seq_preprocessing/${fname%_barcodes_positions_lookup.csv*}

    if [[ "$barcode_lookup" == "$1" ]]
      then
        if [ ! -d "${dname%seqkit_filtering*}RbTn-seq_preprocessing/" ]
          then
            echo "making RbTn-seq_preprocessing directory"
    		    mkdir ${dname%seqkit_filtering*}RbTn-seq_preprocessing/
	      fi
    fi

  # add end point to all mapped reads
  python accessory_scripts/add_mapped_end_points_RB.py -i $barcode_lookup \
                                                 -o $barcode_cluster_out_path

  # find closest genes
  bedtools closest -D a -k 2 -t all \
  -a ${barcode_cluster_out_path}_barcode_positions_bed_final.bed \
  -b ../../tn-seq_data/NGS_ref_genomes/MDS42/AP012306_sorted.gff3 > ${barcode_cluster_out_path}_bedtools_closest_reads.bed

  # annotate each read with closest gene-barcode combination
  python accessory_scripts/annotate_reads_RB.py -i ${barcode_cluster_out_path}_barcode_positions_final.csv \
                                                 -c ${barcode_cluster_out_path}_bedtools_closest_reads.bed \
                                                 -o $barcode_cluster_out_path

done
