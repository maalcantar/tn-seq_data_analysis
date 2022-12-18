#!/bin/bash

######################################################

# This script performs the following:
# 1) uses samtools flagstat to obtain mapping statistics
# 2) uses a custom python script to calculate percent of reads that mapped (this
# is implemented as described in: https://www.biostars.org/p/335579/ )


### PARAMETERS
## $1: list of bam files

# Example: sh calculate_percent_mapped_reads.sh \
# ../../../tn-seq_data/tn-seq_outputs/mae-005/fastp_output/*_out_R1.fastq.gz

### OUTPUTS
# ${bowtie_stats_outpath}_num_read_pairs: number of read pairs (from filtered
# fastq) for each sample
# ${bowtie_stats_outpath}_mapped_stats.csv: summary of mapping statistics

######################################################

# Activate virtual enviroment
source activate ngs_processing_source_fastp

for fastq_file in $@
  do
    fname=$(basename $fastq_file)
  	dname=$(dirname $fastq_file)
  	fpath=$dname/${fname%.fastq.gz*}

    bowtie_stats_outpath=${dname%fastp_output*}bowtie2_stats

    if [[ "$fastq_file" == "$1" ]]
      then
        if [ ! -d "$bowtie_stats_outpath" ]
          then
            echo "making bowtie2_stats directory"
    		    mkdir $bowtie_stats_outpath
	      fi

        # create / overwrite text file that will contain number of QC'ed reads
        # for each sample
      echo "Number of read pairs" > $bowtie_stats_outpath/num_read_pairs.txt
    fi

    echo "calculating stats for:" $fastq_file
    # calculate number of reads in current sample
    num_reads_tmp=$(echo $(gzcat $fastq_file|wc -l)/4|bc)
    # write sample name and number of reads to output text file
    echo $fname : $num_reads_tmp >> $bowtie_stats_outpath/num_read_pairs.txt

    echo __________new sample_____________________
done

# pass text file to python script that summarizes alignment stats
echo $bowtie_stats_outpath/num_read_pairs.txt

python calc_mapped_read_stats.py -i $bowtie_stats_outpath/num_read_pairs.txt
