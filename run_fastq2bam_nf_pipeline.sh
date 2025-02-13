#!/bin/env bash

#SBATCH --mem=40GB
#SBATCH --mail-type=FAIL,END
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --job-name=nextflow_chip

source $HOME/.bashrc_rj_test.sh

conda activate nextflow_three

########## for SE data ###############
# --SE parameter for pair end reads
# when using SE do --single_end_reads and give the path to your single end reads with a glob pattern if you have other files in that directory you dont want ex: path/to/single_end_reads/*file*.fastq
# if you have an adapter sequence use --ada_seq, then specify the sequence with --adapter_seq_str which will take the string sequence
# use --genome and give the path including the file of your genome of choice
# --BL  : this parameter is to tell the pipeline use the process that filters for black list regions
# --blacklist_path : give the path to the blacklist bed file you have and include the file in the path. the defualt used is a path to the hg19 v2 black list. so if using a different species or a different human genome use the correct blacklist and not the default.
########################################

########### for PE data ##############
# I will not put an option to specify adapters and put your own sequence for the Pair End part of this pipeline
# The reason being i specified in fastp that we will look adapters for PE and just trim them. read the parameters used for the fastp tool in the fastp_PE process

# --genome : give the path to the reference genome file you want to use. if you dont specify then the defualt hg19 genome will be used
# --PE : lets the pipeline know you have pair end data
# --paired_end_reads : you need to then also use this parameter to specify the path and the glob pattern to get your forward and reverse reads together in the same input channel. EX: path/to/pair_end_reads/*my_pair_end_file*_{R1,R2}*.fastq 
# --BL : you need to specify BL(black list) if you want the pipeline to do blacklist region filtering.
# --blacklist_path : once you choose --BL, use this parameter to specify the blacklist bed file if you changed the genome from the defualt genome. you dont have to do this if you didnt use the --genome parameter to choose a different genome.

######################################


#nextflow run fastq2bam_nextflow_pipeline.nf -profile 'nasa_pipeline ' \
#-resume \
#--SE \
#--single_end_reads '/rugpfs/fs0/risc_lab/store/hcanaj/HC_ENDseq_Novaseq_010925/read1_fastqs/*_1.fastq.gz' \
#--ada_seq --adapter_seq_str 'AGATCGGAAGAGC' \
#--BL 




nextflow run fastq2bam_nextflow_pipeline.nf -profile 'nasa_pipeline' \
-resume \
--PE \
--BL 

 

 # NOTE: If you want to make your own nextflow diagram to see how the pipeline works run this command

# -with-dag fastq2bam_nf_pipeline_flowchart.png, it seems like this alone overwrites the resume option
# add -preview to render it without having to run the pipeline 