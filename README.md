# NEXDEP-Nextflow_DNA_Epigenomic_Pipeline
Welcome to NEXDEP our Nextflow DNA Epigenomic Pipeline. You can take your DNA based reads from fastq files to bam files. This pipeline will also include processes and workflows for any DNA data that is created from ChIP-seq, Cut&Run, Cut&Tag, ATAC-seq, Bisulfate Methylation and more assays that rely on some form of DNA alignment. 


### This is the pipeline so far, it will take your dna based reads and preprocess them from the fastq to bam level. The pipeline handles ATAC data and shifts the genomic coordinates in preprocessing. There are many parameters so read through them and choose the best one that works for your project. This pipeline so far includes workflows specific for hera's project with NASA data, so those parameters will only be used for that project. I will add workflows that will take the dna based bam files and call peaks and other normal downstream analysis methods for this type of data. MORE COMMING SOON

# How to prepare your data (naming of files/fastqs )

```
This is the standard naming scheme that will allow the pipeline to take your fastq files, and for workflows I make in the future that take bam files; you will have to follow this naming scheme

I will call each area between the underscore as a field and only use underscores (_) to separate fields not dashes (-)

# Field 1: condition label. This is like the control and treatment names. In the H1 project for example, the two conditions were scrm and H1low

# Field 2: experiment type. This can practically be anything, but it should in theory group the conditions so you can compare the control vs the treatment conditions. In the H1 project the conditions were histone marks like H3k27me3 and H3k9me3 and others. So you can compare the scrm vs H1low conditions in the H3k27me3 experiment type, or the scrm vs H1low conditions in the H3k9me3 experiment type. This can also be time points if in a different experiment like comparing treatment vs control conditions in the timepoint 1 experiment type.

# Field 3: replicate label. This can be labeled the way you want but and example will be "r1" for replicate 1 or "r2" for replicate 2. Or you can say "replicate1" and "replicate2" as long as an underscore isnt separating anything in the naming scheme of replicates label. Now this is better to use biological replicates here and merge all technical replicates if you have both bio and tech reps. But if you only have tech reps just keep them separate and treat them as bio reps

# Field 4...: the pipeline does not care about any of the fields after field 3. You can put what ever you like here to remind you of the type of files you have.

Example:
scrm_H3k27me3_r1.fastq
H1low_H3k27me3_r1.fastq
scrm_H3k9me3_r1.fastq
H1low_H3k9me3_r1.fastq

or bam example:
scrm_H3k27me3_r1.bam
H1low_H3k27me3_r1.bam
scrm_H3k9me3_r1.bam
H1low_H3k9me3_r1.bam



```


# How to run the pipeline

**Make an sbatch script if you're on an HPC with SLURM**

```
#SBATCH --mail-type=FAIL,END
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=nextflow_chip
```
**You need to include the bash file that contains all the conda environments we will use**
```
source /lustre/fs4/home/rjohnson/.bashrc_rj_test.sh
```

**Now, the only environment you need to activate in your scrip is nextflow**

```
conda activate nextflow_three
```

**Next, read the parameters section and choose which ones you will add for your analysis**
```
nextflow run fastq2bam_nextflow_pipeline.nf -profile 'fastq2bam2_pipeline' \
-resume \
--genome '/lustre/fs4/risc_lab/store/risc_data/downloaded/hg38/genome/Sequence/WholeGenomeFasta/genome.fa' \
--PE \
--BL \
--blacklist_path '/lustre/fs4/risc_lab/store/risc_data/downloaded/hg38/blacklist/hg38-blacklist.v2.bed' \
--paired_end_reads '/ru-auth/local/home/rjohnson/pipelines/peak_calling_analysis_pipeline/test_published_data/sra_data/*{_1,_2}*' \
--use_effectiveGenomeSize \
--num_effectiveGenomeSize '2864785220' 
```

# Parameters Section

**Parameters for Pair-End Data**

```
I will not put an option to specify adapters and put your own sequence for the Pair End part of this pipeline
The reason being i specified in fastp that we will look adapters for PE and just trim them. read the parameters used for the fastp tool in the fastp_PE process

--test: parameter will make the pipeline only take 3 of your fastq files(or fastq pairs in pair end) in your directory of many fastq files. without this the pipeline will run and process all of your fastq files

--ATAC : if you have atac-seq data, please specify this parameter

--bisulfate_methylation : This parameter will let the pipeline know that you have a bisulfate experiment and fastq files that need to be processed/ aligned with that in mind. This parameter will also make deeptools use a bin size of 10 when creating bigwig and bedgraph files for visualization in IGV, the bedgraph file will also have the regions of zero coverage not reported for a cleaner view. This parameter will also make use of the tool MethylDackel to give back exact locations of the methylated Cytosines in a CpG context.

--genome : give the path to the reference genome file you want to use. if you dont specify then the defualt hg19 genome will be used

--PE : lets the pipeline know you have pair end data
--paired_end_reads : you need to then also use this parameter to specify the path and the glob pattern to get your forward and reverse reads together in the same input channel. EX: path/to/pair_end_reads/*my_pair_end_file*_{R1,R2}*.fastq 

--BL : you need to specify BL(black list) if you want the pipeline to do blacklist region filtering.

--blacklist_path : once you choose --BL, use this parameter to specify the blacklist bed file if you changed the genome from the defualt genome. you dont have to do this if you didnt use the --genome parameter to choose a different genome.

--use_effectiveGenomeSize : this should be called if you want the pipeline to use the path where deeptools bamcoverage will take the effective genome size. this parameter does not take the number see next parameter
--num_effectiveGenomeSize : if you used the parameter --use_effectiveGenomeSize then you need to use this one also. this one takes the number as a str and you can go to deeptools website to find the correct effective genome size number to use that here: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html

--spike_in : lets the pipeline know it should run the workflow for spike ins
--t7, --lambda : choose one or both and the pair end spike in workflow for these will be executed. Have to use with parameter --spike_in
--yeast : choose the yeast spike in if you have ricc_seq data or data that wants yeast spike in, and the pair end spike end will run. must use with parameter --spike_in

--give_peak_files : this is a parameter for the nasa project. Put all peaks in a directory and give the glob pattern '*.narrowPeak' for the pipeline to find your narrow peak files and use them. give the absolute path followed by the glob pattern. if this is not specified then the pipeline will used the peak files that I choose from the data I have access to.

--depth_intersection : this new parameter should be for anyone that has alignment bam files and want to check it's depth by seeing how many reads intersect with a given set of already created peak files

--end_seq or --gloe_seq : the end_seq data and the gloe_seq data names are ordered differently have the user specify if end seq or gloe seq so i can use the correct order for geting name metadata


```









**Parameters For Single-End Data**
```
--test: parameter will make the pipeline only take 3 of your fastq files(or fastq pairs in pair end) in your directory of many fastq files. without this the pipeline will run and process all of your fastq files

--ATAC : if you have atac-seq data, please specify this parameter

--bisulfate_methylation : This parameter will let the pipeline know that you have a bisulfate experiment and fastq files that need to be processed/ aligned with that in mind. This parameter will also make deeptools use a bin size of 10 when creating bigwig and bedgraph files for visualization in IGV, the bedgraph file will also have the regions of zero coverage not reported for a cleaner view. This parameter will also make use of the tool MethylDackel to give back exact locations of the methylated Cytosines in a CpG context.

--SE parameter for pair end reads
when using SE do
--single_end_reads and give the path to your single end reads with a glob pattern if you have other files in that directory you dont want ex: path/to/single_end_reads/*file*.fastq

if you have an adapter sequence use --ada_seq, then specify the sequence with --adapter_seq_str which will take the string sequence

--genome: and give the path including the file of your genome of choice

--BL  : this parameter is to tell the pipeline use the process that filters for black list regions

--blacklist_path : give the path to the blacklist bed file you have and include the file in the path. the defualt used is a path to the hg19 v2 black list. so if using a different species or a different human genome use the correct blacklist and not the default.

--use_effectiveGenomeSize : this should be called if you want the pipeline to use the path where deeptools bamcoverage will take the effective genome size. this parameter does not take the number see next parameter

--num_effectiveGenomeSize : if you used the parameter --use_effectiveGenomeSize then you need to use this one also. this one takes the number and you can go to deeptools website to find the correct effective genome size number to use here: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html

--spike_in : lets the pipeline know it should run the workflow for spike ins

--t7, --lambda : choose one or both and the single end spike in workflow for these will be executed. Have to use with parameter --spike_in
--yeast : choose the yeast spike in if you have ricc_seq data or data that wants yeast spike in, and the single end spike end will run. must use with parameter --spike_in

--give_peak_files : this is a parameter for the nasa project. Put all peaks in a directory and give the glob pattern '*.narrowPeak' for the pipeline to find your narrow peak files and use them. give the absolute path followed by the glob pattern. if this is not specified then the pipeline will used the peak files that I choose from the data I have access to.

--depth_intersection : this new parameter should be for anyone that has alignment bam files and want to check it's depth by seeing how many reads intersect with a given set of already created peak files

--end_seq or --gloe_seq : the end_seq data and the gloe_seq data names are ordered differently have the user specify if end seq or gloe seq so i can use the correct order for geting name metadata
```



