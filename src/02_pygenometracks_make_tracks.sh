#!/bin/bash

######################################################

# This script performs the following (using pyGenomeTracks: https://pygenometracks.readthedocs.io/en/latest/)
# Makes a pyGenomTracks .ini file for each sample that can then be edited by the user before visualization using pyGenomeTracks
# The .ini file encodes a figure that contains the reference genome with features, and then the peaks for each index sample below it
# This script also compiles all of the peaks from all of the samples into a single track

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
# Example: ./02_pygenometracks_make_tracks.sh
#          ../../tn-seq_data/Quintara_NGS_Raw/2021-11-15_tn107-109_p146_arbutin_MiSeq/Fastq/mae-001_metadata.txt
#          ../../tn-seq_data/Quintara_NGS_Raw/2021-11-15_tn107-109_p146_arbutin_MiSeq/Fastq/*R1_001.fastq.gz

### REQUIRED ADDITIONAL FILES
## Reference genomes (in ../../tn-seq_data/NGS_ref_genomes/strain_name/)
## Pileup files from MACS3 representing peaks across the genome (${macs_out_path}_macs3_treat_pileup.bdg)
## Intersect files for MACS3 peaks with genome features (${bedtools_out_path}_peak_intersect_genome.bed)

### OUTPUTS

## PYGENOMETRACKS
## ${fpath}_tracks.ini: closest genes in the genome to the summits
## ${fpath}/all_tracks.ini: peaks with extra bases added

######################################################

# Activate virtual enviroment
source activate pygenometracks

nproc=6

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
    echo Strain $strain

	if [[ "$fastq" == "$1" ]]
		then
			if [ ! -d "../tn-seq_data/tn-seq_outputs/$experiment_name/pygenometracks_output/" ]
		  then
		    mkdir ../../tn-seq_data/tn-seq_outputs/$experiment_name/pygenometracks_output/
			fi
	fi

    # Make a new subfolder to store the pyGenomeTracks .ini files
    pygenometracks_out_path=../../tn-seq_data/tn-seq_outputs/$experiment_name/pygenometracks_output/${fname%_R[12]_001.fastq*}

	bowtie_out_path=../../tn-seq_data/tn-seq_outputs/$experiment_name/bowtie2_output/${fname%_R[12]_001.fastq*}
	macs_out_path=../../tn-seq_data/tn-seq_outputs/$experiment_name/macs3_output/${fname%_R[12]_001.fastq*}
	bedtools_out_path=../../tn-seq_data/tn-seq_outputs/$experiment_name/bedtools_output/${fname%_R[12]_001.fastq*}



	# Load the gene bedfiles files from the folder NGS_ref_genomes
	if [[ "$strain" == "MG1655" ]]
		then
			echo "Strain is MG1655."
			gene_bedfile="../../tn-seq_data/NGS_ref_genomes/MG1655/U00096_PGT_format.bed"

	elif [[ "$strain" == "MDS42" ]]
		then
			echo "Strain is MDS42"
			gene_bedfile="../../tn-seq_data/NGS_ref_genomes/MDS42/AP012306_PGT_format.bed"

	else
		echo $strain "strain is not currently supported -- aborting."
		exit 1
	fi

	echo Beginning $fname
    # Make the tracks.ini file for a single sample
    #make_tracks_file --trackFiles $gene_bedfile ${macs_out_path}_macs3_treat_pileup.bdg  -o ${pygenometracks_out_path}_track.ini
		make_tracks_file --trackFiles $gene_bedfile ${macs_out_path}_macs3_bedcomp.bdg  -o ${pygenometracks_out_path}_track_bedcomp.ini
    # Appened the sample name to the growing list of files, so that they can all be combined into a single track
    list_filenames+=("${bedtools_out_path}_peak_intersect_genome.bed" "${macs_out_path}_macs3_treat_pileup.bdg")

done

# Join arguments with delimiter
# @Params
# $1: The delimiter string
# ${@:2}: The arguments to join
# @Output
# >&1: The arguments separated by the delimiter string
array::join() {
  (($#)) || return 1 # At least delimiter required
  local -- delim="$1" str IFS=
  shift
  str="${*/#/$delim}" # Expand arguments with prefixed delimiter (Empty IFS)
  echo "${str:${#delim}}" # Echo without first delimiter
}

array::join ' ' "${list_filenames[@]}"

# Make a single track file that is the sum of all the peaks of all of the samples
make_tracks_file --trackFiles ${list_filenames[@]} -o ../../tn-seq_data/tn-seq_outputs/$experiment_name/pygenometracks_output/all_samples_tracks.ini
