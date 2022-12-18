#!/bin/bash

######################################################

# This script performs the following (using pyGenomeTracks: https://pygenometracks.readthedocs.io/en/latest/)
# Takes in pyGenomeTracks .ini files and then creates a pyGenomeTracks figure .png file

### PARAMETERS
## $1: metadata in .txt format (tab separated, remember to remove weird extra characters).
# headers should include:
# sample_number (e.g. 4), experiment_name (e.g. mae-001),
# sample_name (e.g. S1_L001_R1_001.fastq.gz) - should match .fastq.gz file names,
# sample_type (e.g. endpoint, pooled, startpoint), strain_name (e.g. MG1655 or MDS42),
# helper_plasmid (e.g. pENG133), outgrowth_media (e.g. lb), and inducer (e.g. atc)
## $2: list of fastq or fastq.gz files to be aligned
# Only specify one of each read pair; the other
# filename is assumed (R2 if R1, or R1 if R2.)
# Example: ./04_pygenometracks_viz_tracks.sh
#          ../../tn-seq_data/Quintara_NGS_Raw/2021-11-15_tn107-109_p146_arbutin_MiSeq/Fastq/mae-001_metadata.txt
#          ../../tn-seq_data/Quintara_NGS_Raw/2021-11-15_tn107-109_p146_arbutin_MiSeq/Fastq/*R1_001.fastq.gz

### REQUIRED ADDITIONAL FILES
## pyGenomeTracks .ini files for each sample
## A pyGenomeTracks .ini file for all of the samples combined

### OUTPUTS

## PYGENOMETRACKS
## ${fpath}_track.png: a figure of a single track for each sample
## ${fpath}_regional_images/*.png: a folder for each sample .ini file that contains all of the regional figures for each insertion
## ${fpath}_all_samples_track.png: a figure of all the peaks from each sample combined

######################################################

# Activate virtual enviroment
source activate pygenometracks

nproc=6

# Save metadata to variable and shift inputs
METADATA=$1
shift

# Check to make sure files have the expected formats
for fastq in $@
do
	if [ -z $(basename $fastq | grep -i .fastq) ]
	then
		echo $(basename $fastq) "does not have .fastq suffix - aborting"
		exit 1
	fi
done

# Iterate through the samples and for each one, make a pyGenomeTracks .ini file
for fastq in "$@"
do
	fname=$(basename $fastq)
	dname=$(dirname $fastq)
	fpath=$dname/${fname%_R[12]_001.fastq*}
	experiment_name=$(cat $METADATA | grep $fname | cut -f2)
	strain=$(cat $METADATA | grep $fname | cut -f5)
	plasmid=$(cat $METADATA | grep $fname | cut -f6)

    echo Experiment name $experiment_name

    # Make a new subfolder to store the pyGenomeTracks .ini files
    pygenometracks_out_path=../../tn-seq_data/tn-seq_outputs/$experiment_name/pygenometracks_output/${fname%_R[12]_001.fastq*}
    bedtools_out_path=../../tn-seq_data/tn-seq_outputs/$experiment_name/bedtools_output/${fname%_R[12]_001.fastq*}


    [ -d "../../tn-seq_data/tn-seq_outputs/$experiment_name/pygenometracks_output/${fname%_R[12]_001.fastq*}_regional_images" ] && echo "Directory exists."
    mkdir ../../tn-seq_data/tn-seq_outputs/$experiment_name/pygenometracks_output/${fname%_R[12]_001.fastq*}_regional_images

	# Load the reference genome files from the folder NGS_ref_genomes
	if [[ "$strain" == "MG1655" ]]
		then
			echo "Strain is MG1655."
            genome_size=4641652
            chromosome_name=U00096.3

	elif [[ "$strain" == "MDS42" ]]
		then
			echo "Strain is MDS42"
            genome_size=3976195
            chromosome_name=AP012306.1
	else
		echo $strain "strain is not currently supported -- aborting."
		exit 1
	fi

	echo Beginning $fname

    # Make a separate pyGenomeTracks figure for each sample
    pyGenomeTracks --tracks ${pygenometracks_out_path}_track_bedcomp.ini --BED ${bedtools_out_path}_bedtools_sloppy_peaks.bed -o ../../tn-seq_data/tn-seq_outputs/$experiment_name/pygenometracks_output/${fname%_R[12]_001.fastq*}_regional_images/peak.png
    #pyGenomeTracks --tracks ${pygenometracks_out_path}_track_bedcomp.ini --region $chromosome_name:0-$genome_size -o ${pygenometracks_out_path}_global_bedcomp.png


done

# Make a single pyGenomeTracks figure for all samples
#pyGenomeTracks --tracks ../../tn-seq_data/tn-seq_outputs/$experiment_name/pygenometracks_output/all_samples_tracks.ini --region $chromosome_name:0-$genome_size -o ../../tn-seq_data/tn-seq_outputs/$experiment_name/pygenometracks_output/all_samples_tracks.png
