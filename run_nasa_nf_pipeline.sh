#!/bin/env bash

#SBATCH --mem=40GB
#SBATCH --mail-type=FAIL,END
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --job-name=nextflow_chip

source $HOME/.bashrc_rj_test.sh

conda activate nextflow_three

# --PE parameter for single end reads
# --SE parameter for pair end reads
# when useing SE do --single_end_reads and give the path to your single end reads with a glob pattern if you have other reads in that directory you dont want
# if you have an adapter sequence use --ada_seq, then specify the sequence with --adapter_seq_str which will take the string sequence
# use --genome and give the path including the file of your genome of choice
# --BL parameter : to tell the pipeline use the process that filters for black list regions
# --blacklist_path : give the path to the blacklist bed file you have and include the file in the path. the defualt used is a path to the hg19 v2 black list. so if using a different species or a different human genome use the correct blacklist and not the default.


########### for PE data ##############
# I will not put an option to specify adapters and put your own sequence for the Pair End part of this pipeline
# The reason being i specified in fastp that we will look adapters for PE and just trim them. read the parameters used in the fastp_PE process

######################################


#nextflow run nasa_pipeline.nf -profile 'nasa_pipeline ' \
#-resume \
#--SE \
#--single_end_reads '/rugpfs/fs0/risc_lab/store/hcanaj/HC_ENDseq_Novaseq_010925/read1_fastqs/*_1.fastq.gz' \
#--ada_seq --adapter_seq_str 'AGATCGGAAGAGC' \
#--BL \



nextflow run nasa_pipeline.nf -profile 'nasa_pipeline' \
-resume \
--PE \
--BL 

#-with-dag fastq2bam_nf_pipeline_flowchart.png
 

 # NOTE: If you want to make your own nextflow diagram to see how the pipeline works run this command

# -with-dag fastq2bam_nf_pipeline_flowchart.png, it seems like this alone overwrites the resume option
# add -preview to render it without having to run the pipeline 