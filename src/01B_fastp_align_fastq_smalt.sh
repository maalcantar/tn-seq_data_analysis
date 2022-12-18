#!/bin/bash

######################################################

# This script performs the following:
# 1) adapter removal and quality trimming with fastp
# (https://github.com/OpenGene/fastp )
# 2) aligns reads to a reference genome using smalt
# 3) (optional) call peaks using macs3
# 4) (optional) finds genomic features overlapping with peaks
# Run the script from within its own directory (i.e. the src folder directory of tn-seq_data_analysis)

### PARAMETERS
## $1: metadata in .txt format (tab separated, remember to remove weird extra characters).
# headers should include:
# sample_number (e.g. 4), experiment_name (e.g. mae-001),
# sample_name (e.g. S1_L001_R1_001.fastq.gz) - should match .fastq.gz file names,
# sample_type (e.g. endpoint, pooled, startpoint), strain_name (e.g. MG1655 or MDS42),
# helper_plasmid (e.g. pENG133), outgrowth_media (e.g. lb), and inducer (e.g. atc)
## $2: list of fastq or fastq.gz files to be aligned.
# Only specify one of each read pair; the other
# filename is assumed (R2 if R1, or R1 if R2.)
# Example: sh 01B_fastp_align_fastq_smalt.sh
#          ../../tn-seq_data/Quintara_NGS_Raw/2022-04-18_MAE-006_L-serine_growth/Fastq/mae-006_metadata.txt
#          ../../tn-seq_data/Quintara_NGS_Raw/2022-04-18_MAE-006_L-serine_growth/Fastq/*R1_001.fastq.gz

### Required additional files
## Reference genomes (in ../../tn-seq_data/NGS_ref_genomes/strain_name/)
# These should be .gff3 files downloaded from NCBI
# Sort .gff3 files if necessary (e.g. bedtools sort -i AP012306.gff3 > AP012306_sorted.gff3)
# Use smalt-index to build the reference genomes in the same folder (outputs .bt2, .rev.bt2)
# Make a genome file for bedtools: chomosome_name \t number_bases (e.g. "AP012306.1 3976195"), save as .txt

### OUTPUTS
## FASTP:
## ${fpath}.json: fastp report in json format
## ${fpath}.html: fastp report in html formats
## ${fpath}_out_R1.fastq.gz: quality filtered forward read
## ${fpath}_out_R1.fastq.gz: quality filtered reverse read

## SMALT + SAMTOOLS
## ${fpath}_out_R1.sam: alignment file in sam format
## ${fpath}_out_R1.bam: alignment file in bam format
## ${fpath}.sorted.bam: alignment file in sorted bam format

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

######################################################

# Activate virtual enviroment
source activate ngs_processing_source_fastp

nproc=7

# Maually enter information about the different UMIs
UMI_LEN=9;
ADAPT1_length=27 # trimmed from end of read1
ADAPT2_length=17 # trimmed from end of read2

# Maximum insert size for the alignment (i.e. max distance between paired reads)
max_ins=1000

# Save metadata to variable and shift inputs
METADATA=$1
shift

# Make a list to record all the file names in
declare -a list_filenames=()

# Check to make sure files have the expected formats
for fastq in $@
do
	if [ -z $(basename $fastq | grep -i .fastq) ]
	then
		echo $(basename $fastq) "does not have .fastq suffix - aborting"
		exit 1
	fi
done

for fastq in "$@"
do
	fname=$(basename $fastq)
	dname=$(dirname $fastq)
	fpath=$dname/${fname%_R[12]_001.fastq*}
	experiment_name=$(cat $METADATA | grep $fname | cut -f2)
	strain=$(cat $METADATA | grep $fname | cut -f5)
	plasmid=$(cat $METADATA | grep $fname | cut -f6)

    echo Experiment name $experiment_name
    echo Strain $strain

	fastp_out_path=../../tn-seq_data/tn-seq_outputs/$experiment_name/fastp_output/${fname%_R[12]_001.fastq*}
	smalt_out_path=../../tn-seq_data/tn-seq_outputs/$experiment_name/smalt_output/${fname%_R[12]_001.fastq*}
	macs_out_path=../../tn-seq_data/tn-seq_outputs/$experiment_name/macs3_output/${fname%_R[12]_001.fastq*}
	bedtools_out_path=../../tn-seq_data/tn-seq_outputs/$experiment_name/bedtools_output/${fname%_R[12]_001.fastq*}

	echo Beginning $fname

	if [[ "$fastq" == "$1" ]]
		then
			if [ ! -d "../tn-seq_data/tn-seq_outputs/$experiment_name/fastp_output/" ]
		  then
		    mkdir ../../tn-seq_data/tn-seq_outputs/$experiment_name/fastp_output/
			fi

			if [ ! -d "../tn-seq_data/tn-seq_outputs/$experiment_name/smalt_output/" ]
		  then
		    mkdir ../../tn-seq_data/tn-seq_outputs/$experiment_name/smalt_output/
			fi

			if [ ! -d "../tn-seq_data/tn-seq_outputs/$experiment_name/macs3_output/" ]
		  then
		    mkdir ../../tn-seq_data/tn-seq_outputs/$experiment_name/macs3_output/
			fi

			if [ ! -d "../tn-seq_data/tn-seq_outputs/$experiment_name/bedtools_output/" ]
		  then
		    mkdir ../../tn-seq_data/tn-seq_outputs/$experiment_name/bedtools_output/
            fi
	fi

	# Load the reference genome files from the folder NGS_ref_genomes
	if [[ "$strain" == "MG1655" ]]
		then
			echo "Strain is MG1655."
			genome="../../tn-seq_data/NGS_ref_genomes/MG1655/U00096_smalt_idx"
			sorted_gene_list_gff="../../tn-seq_data/NGS_ref_genomes/MG1655/U00096_sorted.gff3"
			gene_bedfile="../../tn-seq_data/NGS_ref_genomes/MG1655/U00096_sorted.bed"
			genome_file_sloppy="../../tn-seq_data/NGS_ref_genomes/MG1655/U00096_genome_file.txt"

	elif [[ "$strain" == "MDS42" ]]
		then
			echo "Strain is MDS42"
			genome="../../tn-seq_data/NGS_ref_genomes/MDS42/AP012306_smalt_idx"
			sorted_gene_list_gff="../../tn-seq_data/NGS_ref_genomes/MDS42/AP012306_sorted.gff3"
			gene_bedfile="../../tn-seq_data/NGS_ref_genomes/MDS42/AP012306_sorted.bed"
			genome_file_sloppy="../../tn-seq_data/NGS_ref_genomes/MDS42/AP012306_genome_file.txt"
	else
		echo $strain "strain is not currently supported -- aborting."
		exit 1
	fi

	if [[ "$plasmid" == "pENG133" ]]
		then
			echo "Plasmid is pENG133 (tetR-based)."
			plasmid_sequence="../../tn-seq_data/NGS_ref_genomes/plasmids/pENG133"
			#add additional plasmid-specific files
	elif [[ "$plasmid" == "pENG134" ]]
		then
			echo "Plasmid is pENG134 (araC-based)."
			plasmid_sequence="../../tn-seq_data/NGS_ref_genomes/plasmids/pENG134"
			# add additional plasmid-specific files

	elif [[ "$plasmid" == "pENG131" ]]
		then
			echo "Plasmid is pENG131 (tetR-based)."
			plasmid_sequence="../../tn-seq_data/NGS_ref_genomes/plasmids/pENG131"
			# add additional plasmid-specific files
	else
		echo $plasmid "plasmid is not currently supported -- aborting."
		exit 1
	fi

	# Append the file name to the accumulating list of file names to be used later
	list_filenames+=("${fpath}_bedtools_sloppy_peaks.bed")

  if [ ! -e ${fastp_out_path}_out_R1.fastq.gz ]
  then
	  # First pass to merge reads and trim adapters
    fastp \
    -i ${fpath}_R1_001.fastq.gz \
    -I ${fpath}_R2_001.fastq.gz \
    -o ${fastp_out_path}_out_R1.fastq.gz \
    -O ${fastp_out_path}_out_R2.fastq.gz \
    -U --umi_loc=read2 --umi_len=$UMI_LEN \
    --trim_front1=$ADAPT1_length \
    --trim_front2=$ADAPT2_length \
    -j ${fastp_out_path}.json \
    -h ${fastp_out_path}.html
	fi

  # if [ ! -e ${fpath}.sam ]
  # then
    # align to genome, output .bam files for each sample
		echo "mapping with smalt"
		smalt map \
		-n $nproc \
		-x \
		-y 0.96 \
		-r -1 \
		$genome \
		${fastp_out_path}_out_R1.fastq.gz \
		${fastp_out_path}_out_R2.fastq.gz | samtools view -bS - > ${smalt_out_path}.bam

  # **** everything below is optional. to run, just uncomment desired blocks ****

    # **** add code here to align to the plasmid as well ****
  # fi
	#
	# if [ ! -e ${fpath}_plasmid.sam ]
  # then
  #   # align to genome, output .sam files for each sample
  #   bowtie2 \
  #     --sensitive-local \
  #     -x $plasmid_sequence \
  #     -1 ${fastp_out_path}_out_R1.fastq.gz \
  #     -2 ${fastp_out_path}_out_R2.fastq.gz \
  #     --maxins $max_ins \
  #     --no-mixed \
  #     --no-discordant \
  #     --no-unal \
  #     --quiet \
  #     -p $nproc | samtools view -bS - > ${bowtie_out_path}_plasmid.bam
	#
  # fi

  # convert .sam format to .bam format
  # samtools view -bS ${bowtie_out_path}.sam > ${bowtie_out_path}.bam

  # sort the .bam file relative to genome base index
  #samtools sort ${bowtie_out_path}.bam -o ${bowtie_out_path}.sorted.bam

  # use macs3 to call peaks for insertions, with no control.
  # genome size is set to 5e6 for e. coli
#   macs3 callpeak -f BAMPE -t ${bowtie_out_path}.bam -g 5e6 -n ${macs_out_path}_macs3 -B -q 0.01
#
# #   # make a bedcomp file using macs3 by comparing the pileup data with the generated control
#   macs3 bdgcmp -t ${macs_out_path}_macs3_treat_pileup.bdg -c ${macs_out_path}_macs3_control_lambda.bdg -o ${macs_out_path}_macs3_bedcomp.bdg -m subtract
# #
# #   # use bedtools to find the closest genes in the genome to the summits identified by macs3 as statistically significant
#   bedtools closest -D a -k 4 -t all -a ${macs_out_path}_macs3_summits.bed -b $sorted_gene_list_gff > ${bedtools_out_path}_bedtools_closest_peaks.bed
#
#   # add some "slop"/"fuzz" arond the peaks to increase the number of interesting/relevant genome featues that it finds and also make a window for plotting later
#   # the window of interest has been manually set here to 1500bp
#   bedtools slop -i ${macs_out_path}_macs3_summits.bed -g $genome_file_sloppy -b 1500 > ${bedtools_out_path}_bedtools_sloppy_peaks.bed
#   # find the genomic features annotated within the sloppy window
#   bedtools intersect -wa -a $gene_bedfile -b ${bedtools_out_path}_bedtools_sloppy_peaks.bed > ${bedtools_out_path}_peak_intersect_genome.bed

  echo _____________new sample____________________________

done
