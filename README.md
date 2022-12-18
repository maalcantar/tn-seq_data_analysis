# A self-propagating, barcoded transposon platform for dynamic gene network re-wiring: Tn_Seq data analysis

This repository contains all code needed to reproduce (RB-)Tn-Seq data processing and analyses described in:
>Max A. English, Miguel A. Alcantar, and James J. Collins. A self-propagating, barcoded transposon platform for dynamic gene network re-wiring. <i>In Revision</i>.

Code authors: Max A. English and Miguel A. Alcantar.

# Installation & requirements  

This repository, including all code needed to reproduce analyses, can be installed using:

~~~
git clone https://github.com/maalcantar/tn-seq_data_analysis.git
cd tn-seq_data_analysis
conda create --name ngs_processing_source_fastp python=3.6
source activate ngs_processing_source_fastp
pip install -r requirements_python.txt #requirements for python scripts and Jupyter Notebooks
~~~

Additional requirements: 
* Bartender v1.1 (https://github.com/LaoZZZZZ/bartender-1.1)
* BEDTools v2.30.0
* Bowtie2 v2.4.4
* Fastp v0.23.2 (https://github.com/OpenGene/fastp)
* MACS v3.0.0a7 (https://github.com/macs3-project/MACS)
* Samtools v1.6 
* Seqkit v2.2.0 (https://bioinf.shenwei.me/seqkit/)
* UMI-tools v1.1.0 (https://github.com/CGATOxford/UMI-tools)

# Directory structure

### source code

All code is in  <code>src/</code>, which contains a combination of bash and python scripts. The numbering at the beginning of each file name indicates the order in which that script should be run. Scripts containing the suffix "C" (e.g., 01C_fastp_align_fastq_RB-tn.sh) are specifically designed to handle sequencing reads from RB-Tn-Seq experiments. All other scripts are specific for sequencing reads from regular Tn-Seq experiments.

### data

Due to the large size of output files, outputs are generally written to a separate directory: <code>tn-seq_data/tn-seq_outputs/</code>. 
Additionally, scripts for read alignment require the user to have reference genomes for *E. coli* MDS42 ([AP012306.1](https://www.ncbi.nlm.nih.gov/nuccore/AP012306)) and *E. coli* MG1655 ([U00096.3](https://www.ncbi.nlm.nih.gov/nuccore/545778205)) in <code>tn-seq_data/NGS_ref_genomes/MDS42/</code> and <code>tn-seq_data/NGS_ref_genomes/MG1655/</code>, respectively. 

Raw sequencing data are publicly available under NCBI Bioproject PRJNAXXXXXX.
