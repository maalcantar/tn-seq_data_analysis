#!/bin/bash

######################################################

# This script performs the following:
# 1) converts .bam file to .sam file (samtools view)
# 2) extracts only mapped reads from .sam file (samtools view)

### PARAMETERS
# $1: List of .bam files convert to mapped .sam files
# Example: sh bam_to_sam.sh
# ../../../tn-seq_data/tn-seq_outputs/mae-005/bowtie2_output_dedup/*.deduplicated.bam

### OUTPUTS
## SAMTOOLS
## ${bowtie_out_path}.deduplicated.mapped.sam: sam file with only mapped reads
# bowtie_out_path is the same path as input bam file

######################################################

# Activate virtual enviroment
source activate ngs_processing_source_fastp

nproc=6

for bam in "$@"
do
  fname=$(basename $bam .bam) # assumes we only ever get a bam file
  dname=$(dirname $bam)

  bowtie_out_path=$dname/$fname

  echo $bam

  samtools view -@ 6 -h ${bam} | samtools view -F 4 -o ${bowtie_out_path}.mapped.sam

  echo _____________new sample____________________________


done
