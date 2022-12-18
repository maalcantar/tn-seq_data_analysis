#!/bin/bash

######################################################

# This script performs the following:
# 1) sorts bam files (samtools sort)
# 2) indexes bam files (samtools index)
# 3) deduplicates aligned reads based on alignment coordinate and UMI
# 4) converts deduplicated bam file to a deduplicated fastq file

### PARAMETERS
# $1: metadata in .txt format (tab separated, remember to remove weird extra characters).
# headers should include:
# sample_number (e.g. 4), experiment_name (e.g. mae-001),
# sample_name (e.g. S1_L001_R1_001.fastq.gz) - should match .fastq.gz file names,
# sample_type (e.g. endpoint, pooled, startpoint), strain_name (e.g. MG1655 or MDS42),
# helper_plasmid (e.g. pENG133), outgrowth_media (e.g. lb), and inducer (e.g. atc)
## $2: List of .bam files to be deduplicated
# Example: sh deduplicate_bam_files_RB.sh
# ../../../tn-seq_data/Quintara_NGS_Raw/mae_011_test/mae-011_metadata.txt
# ../../../tn-seq_data/tn-seq_outputs/mae-011/bowtie2_output_RB/*_L001.bam

### OUTPUTS
## BOWTIE + SAMTOOLS
## ${bowtie_out_path}.deduplicated.bam: bam file deduplicated by UMI and alignment

## MACS3
## ${fpath}_macs3_treat_pileup.bdg: pileup signals
# (normalized according to --scale-to option)
## ${fpath}_macs3_control_lambda.bdg:  local biases estimated for each genomic
# location from the control sample
## ${fpath}_macs3_peaks.xls: excel file with information about called peaks
## ${fpath}_macs3_peaks.narrowPeak: contains the peak locations together with peak
# summit, p-value, and q-value
## ${fpath}_macs3_summits.bed: peak summits locations for every peak
## ${fpath}_macs3_bedcomp.bdg: bedcomp file compares two signal tracks

## BEDTOOLS
## ${fpath}_bedtools_closest_peaks.bed: closest genes in the genome to the summits
## ${fpath}_bedtools_sloppy_peaks.bed: peaks with extra bases added
## ${fpath}_peak_intersect_genome.bed: genomic features annotated within the sloppy window

## ${fastq_dedup_out_path}_out_R1.deduplicated.fastq.gz: compressed deduplicated fastq file with forward reads
## ${fastq_dedup_out_path}_out_R2.deduplicated.fastq.gz: compressed deduplicated fastq file with reverse reads
######################################################

# Activate virtual enviroment
source activate ngs_processing_source_fastp

nproc=7
METADATA=$1
shift

for bam in "$@"
do
  fname=$(basename $bam .bam) # assumes we only ever get a bam file
  dname=$(dirname $bam)

  experiment_name=$(cat $METADATA | grep $fname | cut -f2)
	strain=$(cat $METADATA | grep $fname | cut -f5)
	plasmid=$(cat $METADATA | grep $fname | cut -f6)

  bowtie_out_path=$dname/$fname


  echo Experiment name $experiment_name
  echo Strain $strain
  echo $bam

  bowtie_out_path=../../../tn-seq_data/tn-seq_outputs/$experiment_name/bowtie2_output_dedup/${fname}
  macs_out_path=../../../tn-seq_data/tn-seq_outputs/$experiment_name/macs3_output_dedup/${fname}
  bedtools_out_path=../../../tn-seq_data/tn-seq_outputs/$experiment_name/bedtools_output_dedup/${fname}
  fastq_dedup_out_path=../../../tn-seq_data/tn-seq_outputs/$experiment_name/fastq_dedup/${fname}
  seqkit_filtering_RB_dedup_out_path=../../../tn-seq_data/tn-seq_outputs/$experiment_name/seqkit_filtering_RB_dedup/${fname}


  echo Beginning $fname

  if [[ "$bam" == "$1" ]]
    then
      if [ ! -d "../../../tn-seq_data/tn-seq_outputs/$experiment_name/bowtie2_output_dedup/" ]
      then
        mkdir ../../../tn-seq_data/tn-seq_outputs/$experiment_name/bowtie2_output_dedup/
      fi

      if [ ! -d "../../../tn-seq_data/tn-seq_outputs/$experiment_name/macs3_output_dedup/" ]
      then
        mkdir ../../../tn-seq_data/tn-seq_outputs/$experiment_name/macs3_output_dedup/
      fi

      if [ ! -d "../../../tn-seq_data/tn-seq_outputs/$experiment_name/bedtools_output_dedup/" ]
      then
        mkdir ../../../tn-seq_data/tn-seq_outputs/$experiment_name/bedtools_output_dedup/
            fi
      if [ ! -d "../../../tn-seq_data/tn-seq_outputs/$experiment_name/fastq_dedup/" ]
      then
        mkdir ../../../tn-seq_data/tn-seq_outputs/$experiment_name/fastq_dedup/
            fi
      if [ ! -d "../../../tn-seq_data/tn-seq_outputs/$experiment_name/seqkit_filtering_RB_dedup/" ]
      then
        mkdir ../../../tn-seq_data/tn-seq_outputs/$experiment_name/seqkit_filtering_RB_dedup/
            fi
  fi

  # Load the reference genome files from the folder NGS_ref_genomes
  if [[ "$strain" == "MG1655" ]]
    then
      echo "Strain is MG1655."
      genome="../../../tn-seq_data/NGS_ref_genomes/MG1655/U00096"
      sorted_gene_list_gff="../../../tn-seq_data/NGS_ref_genomes/MG1655/U00096_sorted.gff3"
      gene_bedfile="../../../tn-seq_data/NGS_ref_genomes/MG1655/U00096_sorted.bed"
      genome_file_sloppy="../../../tn-seq_data/NGS_ref_genomes/MG1655/U00096_genome_file.txt"

  elif [[ "$strain" == "MDS42" ]]
    then
      echo "Strain is MDS42"
      genome="../../../tn-seq_data/NGS_ref_genomes/MDS42/AP012306"
      sorted_gene_list_gff="../../../tn-seq_data/NGS_ref_genomes/MDS42/AP012306_sorted.gff3"
      gene_bedfile="../../../tn-seq_data/NGS_ref_genomes/MDS42/AP012306_sorted.bed"
      genome_file_sloppy="../../../tn-seq_data/NGS_ref_genomes/MDS42/AP012306_genome_file.txt"
  else
    echo $strain "strain is not currently supported -- aborting."
    exit 1
  fi

  if [[ "$plasmid" == "pENG133" ]]
    then
      echo "Plasmid is pENG133 (tetR-based)."
      plasmid_sequence="../../../tn-seq_data/NGS_ref_genomes/plasmids/pENG133"
      #add additional plasmid-specific files
  elif [[ "$plasmid" == "pENG134" ]]
    then
      echo "Plasmid is pENG134 (araC-based)."
      plasmid_sequence="../../../tn-seq_data/NGS_ref_genomes/plasmids/pENG134"
      # add additional plasmid-specific files
  else
    echo $plasmid "plasmid is not currently supported."
    # exit 1
  fi

  # sort and index bam files to make compatible with umi_tools
  samtools sort $bam -@ $nproc -o ${bowtie_out_path}.sorted.bam
  samtools index -@ $nproc ${bowtie_out_path}.sorted.bam

  # run umi_tools deduplication by UMI and mapping coordinate
  umi_tools dedup  --umi-separator=: --paired --log=${bowtie_out_path}_dedup_log \
  -I ${bowtie_out_path}.sorted.bam --output-stats=${bowtie_out_path}_deduplicated -S \
  ${bowtie_out_path}.deduplicated.bam

  # # remove sorted bam file and index
  rm ${bowtie_out_path}.sorted.bam
  rm ${bowtie_out_path}.sorted.bam.bai

  # # sort bam file by read name to make compatible with bedtools
  samtools sort -@ $nproc -n -o ${bowtie_out_path}.deduplicated.sorted.bam \
  ${bowtie_out_path}.deduplicated.bam

  # # convert bam file to deduplicated fastq files
  bedtools bamtofastq -i ${bowtie_out_path}.deduplicated.sorted.bam \
  -fq ${fastq_dedup_out_path}_out_tmp_R1.deduplicated.fastq \
  -fq2 ${fastq_dedup_out_path}_out_tmp_R2.deduplicated.fastq

  # bioawk -c fastx '{id=substr($name, 0, length($name) - 11); print ">"id"\n"$seq}'
  # {print "@"$name"\n"$seq"\n+\n"$qual}'
  # remove UMI tag from read name -- this will make read names compatible with
  # previous processing steps and files (such as the barcode lookup file )
  bioawk -c fastx '{id=substr($name, 0, length($name) - 11); print ">"id"\n"$seq"\n+\n"$qual}'   \
  ${fastq_dedup_out_path}_out_tmp_R1.deduplicated.fastq > \
  ${fastq_dedup_out_path}_out_R1.deduplicated.fastq
  bioawk -c fastx '{id=substr($name, 0, length($name) - 11); print ">"id"\n"$seq"\n+\n"$qual}'  \
  ${fastq_dedup_out_path}_out_tmp_R2.deduplicated.fastq > \
  ${fastq_dedup_out_path}_out_R2.deduplicated.fastq

  # extract read names from forward fastq
  # adapted from https://gist.github.com/poeli/34aa3919f74d20f28eb24d5e00f42e89
  # use 'NR%4==1 {print substr($1,2,length($1)-12)}' if you want to remove
  # UMI from read name here instead of through bioawk
  awk 'NR%4==1 {print substr($1,2,length($1))}' \
  ${fastq_dedup_out_path}_out_R1.deduplicated.fastq > \
  ${fastq_dedup_out_path}_dedup_read_names.txt

  # remove temp files that still contain the UMI tag in read name
  rm ${fastq_dedup_out_path}_out_tmp_R1.deduplicated.fastq
  rm ${fastq_dedup_out_path}_out_tmp_R2.deduplicated.fastq

  # compress fastq files
  gzip ${fastq_dedup_out_path}_out_R1.deduplicated.fastq
  gzip ${fastq_dedup_out_path}_out_R2.deduplicated.fastq

  # script to do barcode matching
  python dedup_barcode_lookup_file_RB.py \
  -i ${fastq_dedup_out_path}_dedup_read_names.txt


  # use macs3 to call peaks for insertions, with no control.
  # genome size is set to 5e6 for e. coli
  #macs3 callpeak -f BAMPE -t ${bowtie_out_path}.deduplicated.bam -g 5e6 -n ${macs_out_path}_macs3 -B -q 0.01

  # make a bedcomp file using macs3 by comparing the pileup data with the generated control
  #macs3 bdgcmp -t ${macs_out_path}_macs3_treat_pileup.bdg -c ${macs_out_path}_macs3_control_lambda.bdg -o ${macs_out_path}_macs3_bedcomp.bdg -m subtract
  # use bedtools to find the closest genes in the genome to the summits identified by macs3 as statistically significant
  #bedtools closest -D a -k 4 -t all -a ${macs_out_path}_macs3_summits.bed -b $sorted_gene_list_gff > ${bedtools_out_path}_bedtools_closest_peaks.bed

  # add some "slop"/"fuzz" arond the peaks to increase the number of interesting/relevant genome featues that it finds and also make a window for plotting later
  # the window of interest has been manually set here to 1500bp
  #bedtools slop -i ${macs_out_path}_macs3_summits.bed -g $genome_file_sloppy -b 1500 > ${bedtools_out_path}_bedtools_sloppy_peaks.bed
  # find the genomic features annotated within the sloppy window
  #bedtools intersect -wa -a $gene_bedfile -b ${bedtools_out_path}_bedtools_sloppy_peaks.bed > ${bedtools_out_path}_peak_intersect_genome.bed

  echo _____________new sample____________________________


done

python deduplication_stats.py -i ../../../tn-seq_data/tn-seq_outputs/$experiment_name/bowtie2_output_dedup/
