// first i will activate dsl2
nextflow.enable.dsl=2

// NOTE: i would prefer if we could just keep all the reads in one directory and just use glob pattern to grab what we want

// first process is to download genomes that will be used or just get the path to the genome
//process download_genomes{}

// this is the process for if the adapter sequence is known


process fastp_SE_adapter_known {
    // using this conda yml file that I created. Nextflow will now make its own conda environment from the dependencies found within.
    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/fastp_rj_env.yml'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/fastp_rj'

    publishDir './results_SE/fastp_qc_single_end', mode: 'copy', pattern:'*_fp_filt.fastq'
    publishDir './results_SE/fastp_qc_single_end/html_reports', mode: 'copy', pattern:'*.html'

    input:
    // input names dont have to be the exact same as what is seen in the workflow section of this script
    // as long and they are in the same order, you will have the correct input
    path(fastq_files) // the files
    val(fastq_names) // the names of the files as a value input channel
    val(adapter_seq)


    output:
    path("${out_name}"), emit: filtered_fastqs
    path("${fastq_names}*.html"), emit: fastp_html_reports


    script:

    // getting the output name
    out_name = "${fastq_names}_fp_filt.fastq"


    """
    #!/usr/bin/env bash

    # ASK HERA ABOUT THE ADAPTERS USED AND FIND THEM FROM THE ILLUMINA KIT

    # remember this process is for the single end reads so only one read in and one filtered qc read out
    # adapter trimming is enabled by defualt for single end but use --detect_adapter_for_pe to enable it in paired end data
    # you can specify the adapter sequence if you have them. look at documentation here: https://open.bioqueue.org/home/knowledge/showKnowledge/sig/fastp
    
    # --adapter_sequence: string sequence that represents the adapter used. in the illumina NEBNext libraries adapter read 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA. Adapter read 2 is AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
    # --dedup enables the deduplication to drop the duplicated reads
    # --dup_calc_accuracy: defualt is 3 but can go from 1~6. higher levels use more memory
    # --trim_poly_g: novaSeq data polyG can happen in read tails since G means no signal in the illumina two-color system. fastp can detect and trim them.
    # --qualified_quality_phred: the quality value that a base is qualified defualt is 15
    # --unqualified_percent_limit: how many percent of bases are allowed to be unqualified defualt is 40
    # --umi: not using this yet have to use --umi_loc also and i dont know that, fastp can extract the unique molecular identifiers and append them to the first part of the read names, so the umi's will be present in SAM/BAM records If the UMI is in the reads, then it will be shifted from read so that the read will become shorter. If the UMI is in the index, it will be kept.
    # --overrepresentation_analysis and --overrepresentation_sampling see fastp documentation
    # --html: this is how the report file will be made

    # their fastq2bam single end pipeline used trim_galore tool for trimming adapters.
    # I want to implement a way to choose to let fastp do this by default or use a known sequence that the user can input.
    ###### can add or remove more options when needed #####
    
    fastp \
    --in1 "${fastq_files}" \
    --out1 "${out_name}" \
    --adapter_sequence "${adapter_seq}" \
    --dedup \
    --dup_calc_accuracy 5 \
    --trim_poly_g \
    --qualified_quality_phred 15 \
    --unqualified_percent_limit 40 \
    --overrepresentation_analysis \
    --overrepresentation_sampling 20 \
    --html "${fastq_names}_fastp.html"






    """

}







// This next process will run some qc to look at the fastq files and trim adapters from the single end reads
process fastp_SE {
    // using this conda yml file that I created. Nextflow will now make its own conda environment from the dependencies found within.
    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/fastp_rj_env.yml'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/fastp_rj'

    publishDir './results_SE/fastp_qc_single_end', mode: 'copy', pattern:'*_fp_filt.fastq'
    publishDir './results_SE/fastp_qc_single_end/html_reports', mode: 'copy', pattern:'*.html'

    input:
    // input names dont have to be the exact same as what is seen in the workflow section of this script
    // as long and they are in the same order, you will have the correct input
    path(fastq_files) // the files
    val(fastq_names) // the names of the files as a value input channel



    output:
    path("${out_name}"), emit: filtered_fastqs
    path("${fastq_names}*.html"), emit: fastp_html_reports


    script:

    // getting the output name
    out_name = "${fastq_names}_fp_filt.fastq"


    """
    #!/usr/bin/env bash

    # ASK HERA ABOUT THE ADAPTERS USED AND FIND THEM FROM THE ILLUMINA KIT

    # remember this process is for the single end reads so only one read in and one filtered qc read out
    # adapter trimming is enabled by defualt for single end but use --detect_adapter_for_pe to enable it in paired end data
    # you can specify the adapter sequence if you have them. look at documentation here: https://open.bioqueue.org/home/knowledge/showKnowledge/sig/fastp
    
    # --adapter_sequence: string sequence that represents the adapter used. in the illumina NEBNext libraries adapter read 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA. Adapter read 2 is AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
    # --dedup enables the deduplication to drop the duplicated reads
    # --dup_calc_accuracy: defualt is 3 but can go from 1~6. higher levels use more memory
    # --trim_poly_g: novaSeq data polyG can happen in read tails since G means no signal in the illumina two-color system. fastp can detect and trim them.
    # --qualified_quality_phred: the quality value that a base is qualified defualt is 15
    # --unqualified_percent_limit: how many percent of bases are allowed to be unqualified defualt is 40
    # --umi: not using this yet have to use --umi_loc also and i dont know that, fastp can extract the unique molecular identifiers and append them to the first part of the read names, so the umi's will be present in SAM/BAM records If the UMI is in the reads, then it will be shifted from read so that the read will become shorter. If the UMI is in the index, it will be kept.
    # --overrepresentation_analysis and --overrepresentation_sampling see fastp documentation
    # --html: this is how the report file will be made

    # their fastq2bam single end pipeline used trim_galore tool for trimming adapters.
    # I want to implement a way to choose to let fastp do this by default or use a known sequence that the user can input.
    ###### can add or remove more options when needed #####
    
    fastp \
    --in1 "${fastq_files}" \
    --out1 "${out_name}" \
    --dedup \
    --dup_calc_accuracy 5 \
    --trim_poly_g \
    --qualified_quality_phred 15 \
    --unqualified_percent_limit 40 \
    --overrepresentation_analysis \
    --overrepresentation_sampling 20 \
    --html "${fastq_names}_fastp.html"






    """

}

// this next process is for fastqc tool
// dont forget about multiqc m
process fastqc_SE {
    // using the conda environment 
    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/fastqc_rj_env.yml'
    publishDir './results_SE/fastqc_htmls', mode: 'copy', pattern: '*.html'


    input:

    path(fastq_filt_files)
    val(fastq_filt_names)


    output:
    path("*.html"), emit: fastqc_htmls
    path("*.zip"), emit: fastqc_zip_files

    script:
    out_name = "${fastq_filt_files}"

    """
    #!/usr/bin/env bash

    # I need to add the adapter sequences later, ask hera.

    fastqc "${fastq_filt_files}"
    



    """
}

process multiqc_SE {
    // this yml file doesnt work
    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/multiqc_rj_env.yml'

    //conda '/ru-auth/local/home/rjohnson/miniconda3/envs/multiqc_rj'
    
    publishDir './results_SE/multiQC_collection', mode: 'copy', pattern: '*.html'

    input:
    path(fastp_filt_html)



    output:

    path("*.html"), emit: multiqc_html_collection

    script:


    """
    #!/usr/bin/env bash

    # I think i jsut have to pass the html files to multiqc for it to complie them
    multiqc . \
    --interactive \
    --profile-runtime \
    --title "Single-end QC"




    """

}


// Creating two processes that will index the reference genome

process bwa_index_genome {
    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/bwa_rj_env.yml'

    publishDir './genome_index_bwa', mode: 'copy', pattern: '*'

    input:
    path(ref_genome)
    


    output:

    path("*"), emit: genome_index_files


    script:
    // not getting the basename because bwa expects the index files will have the exact file name of the genome but with an .{ext} on the end. 
    // example genome.fa, genome.fa.amb, genome.fa.ann
    // not genome.fa, genome.amb, genome.ann
    //genome_file_name = "${ref_genome.baseName}"

    """
    #!/usr/bin/env bash

    ############### parameters used ###############
    # -p: a string representing the prefix of the output database [same as the db filename] so i think just the base name of the genome file
    # -a: a string, choosing the algorithm to construct the BWT index. either is or bwtsw. I'll use bwtsw


    bwa index \
    -p "${ref_genome}" \
    -a bwtsw \
    "${ref_genome}"


    """

}

// creating a process that will align the reads to the genome. i will take in the reference genome, the index files, the filtered fastq's and their names

process bwa_align_SE {
    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/bwa_rj_env.yml'

    publishDir './results_SE/bwa_outputs_singleEnd_SAM', mode: 'copy', pattern: '*.sam'
    publishDir './results_SE/sai_alignment_files', mode: 'copy', pattern: '*.sai'


    input:
    path(ref_genome)
    path(genome_index_files)
    path(fastq_filt_files)
    val(fastq_filt_names)


    output:

    path("*.sam"), emit: sam_se_files
    path("*.sai"), emit: sai_align_files


    script:
    sai_output_file = "${fastq_filt_names}_out.sai"
    sam_name = "${fastq_filt_names}.sam"

    """
    #!/usr/bin/env bash

    ############# Parameters used ############
    # first i need to get the sai file by using bwa aln. This gives the SA coordinates of the input reads
    # -t (nThrds) number of threads for multi threading mode. defualt is 1

    # using bwa samse : this will generate alignments in the SAM format given single-end reads
    # the two parameters may be used in the future.
    # -n: takes an integer. max number of alignments to output in the XA tag for reads paired properly
    # -r: takes a string. specify the read group
    #
    ##########################################

    ls .

    bwa aln \
    -t 20 \
    "${ref_genome}"  \
    "${fastq_filt_files}" \
    > "${sai_output_file}"


    bwa samse \
    "${ref_genome}" \
    "${sai_output_file}" \
    "${fastq_filt_files}" \
    > "${sam_name}"


    """
}

process samtools_sort {
    // using the conda yml file for samtools
    // it doesnt work
    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/samtools_rj_env.yml'

    //conda '/ru-auth/local/home/rjohnson/miniconda3/envs/samtools_rj' // this was samtools version 1.3 which doesnt have samtools fixmate option -m

    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/samtools-1.21_spec_env_rj.txt' // that is the explicit file but if that doesnt work try the yml file samtools-1.21_env_rj.yml; and if that doesnt work use the path to the 1.21 environment

    if (params.PE) {

        publishDir './results_PE/sorted_bam_files', mode: 'copy', pattern: '*_sorted.bam'
        publishDir './results_PE/indexed_bam_files', mode: 'copy', pattern: '*.{bai, csi}'
        publishDir './results_PE/flag_stat_log', mode: 'copy', pattern: '*stats.log'

    }

    else if(params.SE) {

        publishDir './results_SE/sorted_bam_files', mode: 'copy', pattern: '*_sorted.bam'
        publishDir './results_SE/indexed_bam_files', mode: 'copy', pattern: '*.{bai, csi}'
        publishDir './results_SE/flag_stat_log', mode: 'copy', pattern: '*stats.log'

    }
    //publishDir './sorted_bam_files', mode: 'copy', pattern: '*_sorted.bam'
    //publishDir './indexed_bam_files', mode: 'copy', pattern: '*.{bai, csi}'
    //publishDir './flag_stat_log', mode: 'copy', pattern: '*stat.log'

    input:
    path(sam_files)


    output:

    path("*_sorted.bam"), emit: sorted_bams
    //tuple path("*.{bai,csi}"), emit: indexed_bams
    //tuple path("*.bai"), path("*.csi"), emit: indexed_bams
    path("*.bai"), emit: indexed_bams
    tuple path("*_sorted.bam"), path("*.bai"), emit: bam_index_tuple

    path("*stats.log"), emit: flag_stats_log
    path("*stats.txt"), emit: norm_stats_txt
    path("*stats.tsv"), emit: tsv_SN_stats

    script:

    // i will start using baseName inside the process since its easier to keep track of different names an uses less inputs into a process
    out_bam_filt = "${sam_files.baseName}_bam_filt.bam"
    out_bam_name_sort = "${sam_files.baseName}_name_ordered.bam"
    out_bam_coor_sort = "${sam_files.baseName}_filt_coor_sorted.bam"
    out_bam_fixmate = "${sam_files.baseName}_fixmate.bam"
    out_bam_final = "${sam_files.baseName}_markdup_filt_coor_sorted.bam"
    flagstats_log_pe = "${sam_files.baseName}_pair_end_stats.log"
    samtools_stats_log_pe = "${sam_files.baseName}_pe_stats.txt"
    tsv_file_with_stats_pe = "${sam_files.baseName}_pe_SN_stats.tsv"

    out_bam_se_filt = "${sam_files.baseName}_se_file.bam"
    out_bam_coor_sort_se = "${sam_files.baseName}_se_file_coor_sorted.bam"
    flagstats_log_se = "${sam_files.baseName}_single_end_stats.log"
    samtools_stats_log_se = "${sam_files.baseName}_se_stats.txt"
    tsv_file_with_stats_se = "${sam_files.baseName}_se_SN_stats.tsv"

    if (params.PE) {

        """
        #!/usr/bin/env bash

        ################# samtools parameters used ################
        # for samtools view
        # --min-MQ or -q : takes an INT and will skip alignments with a MAPQ smaller than INT
        # --bam or -b : output in the bam format
        # this version of bwa didnt recognize --bam or --min-MQ so i just used -b and -q respetively.

        # for samtools sort
        # -o : takes a file. it writes the final sorted output to file rather than standard output
        # -O : write the final output as sam, bam, or cram

        # samtools fixmate : preparing for finding the duplicates
        # -O : specify the format and i choose bam
        # -m : add ms(mate score) tags. these are used by markdup to select the best reads to keep 
        # other possible option to look at is -r : remove secondary and unmapped reads ?
        
        # for samtools markdup : can only be done on coordinate sorted bam files and run it through samtools fixmate first
        # -r : remove duplicate reads
        # --mode or -m : duplicate decision method for paired reads. values are "t" or "s". read documentation but i choose s becasue it tends to return more results. just incase i will remove this option in the pair end mode by adding an if else statement in this process.

        # now for samtools index to get index files
        # -b, --bai: create a bai index; this version of samtools does not support --bai --csi just use -b -c
        # -c, --csi: create a csi index
        # -o, --output: write the output index to a file specified  only when one alignment file is being indexed

        # using samtools flagstat: generate log files so i can use multiqc to get stats of all files into one html file
        # no parameters needed. just need to give the final bam file that went through all the processing
        ###########################################################

        # I should add a samtools filtering. looking to only get mapq scores higher than 30

        samtools view \
        -q 30 \
        -b \
        "${sam_files}" \
        > "${out_bam_filt}" 
        
        
        # first i have to name sort to use fixmate
        samtools sort \
        -o "${out_bam_name_sort}" \
        -n \
        -O bam \
        "${out_bam_filt}"


        samtools fixmate \
        -O bam \
        -m \
        "${out_bam_name_sort}"\
        "${out_bam_fixmate}"


        # now i will coordinate sort here 
        samtools sort \
        -o "${out_bam_coor_sort}" \
        -O bam \
        "${out_bam_fixmate}"


        # i might need to put coordinate sorted bam into markdup
        # this works but removing the duplicates results in the file being very small meaning too many reads were removed that were considered duplicates.
        # this results in the next process deeptools not being able to create a normalized bedgraph file

        #samtools markdup \
        #"\${out_bam_coor_sort}" \
        #"\${out_bam_final}"

        # so i will just use the out file from the coordinate sort samtools sort section instead of using out_bam_final
        samtools index \
        -b \
        "${out_bam_coor_sort}"

        samtools flagstat \
        "${out_bam_coor_sort}" \
        > "${flagstats_log_pe}"

        # adding another way to get stats from each bam file
        samtools stats \
        "${out_bam_coor_sort}" \
        > "${samtools_stats_log_pe}"

        # now only putting the stats into a tsv file
        less "${samtools_stats_log_pe}" | grep ^SN | cut -f 2-3 >  "${tsv_file_with_stats_pe}"
        """



    }
    else if (params.SE) {

        """
        #!/usr/bin/env bash

        ################# samtools parameters used ################
        # for samtools view
        # --min-MQ or -q : takes an INT and will skip alignments with a MAPQ smaller than INT
        # --bam or -b : output in the bam format
        # this version of bwa didnt recognize --bam or --min-MQ so i just used -b and -q respetively.

        # for samtools sort
        # -o : takes a file. it writes the final sorted output to file rather than standard output
        # -O : write the final output as sam, bam, or cram

        # samtools fixmate : preparing for finding the duplicates
        # -O : specify the format and i choose bam
        # -m : add ms(mate score) tags. these are used by markdup to select the best reads to keep 
        # other possible option to look at is -r : remove secondary and unmapped reads ?
        
        # for samtools markdup : can only be done on coordinate sorted bam files and run it through samtools fixmate first
        # -r : remove duplicate reads
        # --mode or -m : duplicate decision method for paired reads. values are "t" or "s". read documentation but i choose s becasue it tends to return more results. just incase i will remove this option in the pair end mode by adding an if else statement in this process.

        # now for samtools index to get index files
        # -b, --bai: create a bai index; this version of samtools does not support --bai --csi just use -b -c
        # -c, --csi: create a csi index
        # -o, --output: write the output index to a file specified  only when one alignment file is being indexed

        # using samtools flagstat: generate log files so i can use multiqc to get stats of all files into one html file
        # no parameters needed. just need to give the final bam file that went through all the processing


        ###########################################################

        # I should add a samtools filtering. looking to only get mapq scores higher than 30

        samtools view \
        -q 30 \
        -b \
        "${sam_files}" \
        > "${out_bam_se_filt}" 
        



        # now i will coordinate sort here 
        samtools sort \
        -o "${out_bam_coor_sort_se}" \
        -O bam \
        "${out_bam_se_filt}"


        # so i will just use the out file from the coordinate sort samtools sort section instead of using out_bam_final
        samtools index \
        -b \
        "${out_bam_coor_sort_se}"

        samtools flagstat \
        "${out_bam_coor_sort_se}" \
        > "${flagstats_log_se}"

        # adding another way to get stats from each bam file
        samtools stats \
        "${out_bam_coor_sort}" \
        > "${samtools_stats_log_se}"

        # now only putting the stats into a tsv file
        less "${samtools_stats_log_se}" | grep ^SN | cut -f 2-3 >  "${tsv_file_with_stats_se}"
        """

    }

    
}

process deeptools_make_bed {
    // this conda env yml file didnt work. have to use the actual env
    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/deeptools_rj_env.yml'
    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/deeptools_rj'

    // this section is just a simple if else statement controlling the directories that are created and when the files end up
    // will copy and paste in other processes that need it.
    if (params.PE) {

        if (params.BL) {

            publishDir './results_PE/bl_filt_bed/bed_graphs_deeptools/', mode: 'copy', pattern: '*'

        }
        else {
            publishDir './results_PE/no_bl_filt/bed_graphs_deeptools/', mode: 'copy', pattern: '*'
        }
    
    }
    else {

        if (params.BL) {

            publishDir './results_SE/bl_filt_bed/bed_graphs_deeptools/', mode: 'copy', pattern: '*'

        }
        else {
            publishDir './results_SE/no_bl_filt/bed_graphs_deeptools/', mode: 'copy', pattern: '*'
        }
    }


    input:

    tuple path(bams), path(index)


    output:
    path("${out_bed_name}*"), emit: bed_files_normalized

    script:

    out_bed_name="${bams.baseName}_normalized_cpm.bed"

    if (params.use_effectiveGenomeSize) {

        """

        ###### Using deeptools parameters ###############

        # first converting the bam file to a bed file using bamCoverage. I can also make a bigwig file if needed, it stores data better but is binary and cannot be opened in text editor
        # -b or --bam: takes the bam file that will be processed
        # -o or --outFileName: is the name you want the output file to have
        # -of or --outFileFormat: is the type of output file you want; either "bigwig" or "bedgraph"
        # --scaleFactor: the computed scaling factor (or 1, if not applicable) will be multiplied by this.
        # -bs or --binSize: are the size of the bins in bases, for output of bigwig or bedgraph. default is 50
        # -p or --numberOfProcessors: this is the number of processers you want to use. Not using this option yet but if needed I will use it.
        # --normalizeUsing: choose the type of normalization
        # bamCoverage offers normalization by scaling factor, Reads Per Kilobase per Million mapped reads (RPKM), counts per million (CPM), bins per million mapped reads (BPM) and 1x depth (reads per genome coverage, RPGC).
        # --effectiveGenomeSize: choose the mappable genome size for your organism of choice used as reference. find length here: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
        # not using effectiveGenomeSize since multiple users will use this pipeline and might not be using the same organism.
        # actually i decided to use effectiveGenomeSize afterall since i can split this process into using it or not.

        # NOTE: since all the files will be processed using this tool and parameters, they will all be directly comparable in UCSC or IGV without needing to edit track heights.
        #################################################


        bamCoverage \
        --bam "${bams}" \
        --outFileName "${out_bed_name}" \
        --outFileFormat "bedgraph" \
        --scaleFactor 1 \
        --binSize 50 \
        --normalizeUsing CPM \
        --effectiveGenomeSize "${params.num_effectiveGenomeSize}"


        """


    }
    else {

        """

        bamCoverage \
        --bam "${bams}" \
        --outFileName "${out_bed_name}" \
        --outFileFormat "bedgraph" \
        --scaleFactor 1 \
        --binSize 50 \
        --normalizeUsing CPM

        """

    }

    /*"""
    ###### Using deeptools parameters ###############

    # first converting the bam file to a bed file using bamCoverage. I can also make a bigwig file if needed, it stores data better but is binary and cannot be opened in text editor
    # -b or --bam: takes the bam file that will be processed
    # -o or --outFileName: is the name you want the output file to have
    # -of or --outFileFormat: is the type of output file you want; either "bigwig" or "bedgraph"
    # --scaleFactor: the computed scaling factor (or 1, if not applicable) will be multiplied by this.
    # -bs or --binSize: are the size of the bins in bases, for output of bigwig or bedgraph. default is 50
    # -p or --numberOfProcessors: this is the number of processers you want to use. Not using this option yet but if needed I will use it.
    # --normalizeUsing: choose the type of normalization
    # bamCoverage offers normalization by scaling factor, Reads Per Kilobase per Million mapped reads (RPKM), counts per million (CPM), bins per million mapped reads (BPM) and 1x depth (reads per genome coverage, RPGC).
    # --effectiveGenomeSize: choose the mappable genome size for your organism of choice used as reference. find length here: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
    # not using effectiveGenomeSize since multiple users will use this pipeline and might not be using the same organism.

    # NOTE: since all the files will be processed using this tool and parameters, they will all be directly comparable in UCSC or IGV without needing to edit track heights.
    #################################################

    bamCoverage \
    --bam "${bams}" \
    --outFileName "${out_bed_name}" \
    --outFileFormat "bedgraph" \
    --scaleFactor 1 \
    --binSize 50 \
    --normalizeUsing CPM


    """*/
}

process bedtools_filt_blacklist {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/bedtools_rj'

    //publishDir './blacklist_filt_bam', mode: 'copy', pattern: '*.bam'

    if (params.PE) {

        publishDir './results_PE/blacklist_filt_bam', mode: 'copy', pattern: '*.bam'
    
    }
    else {

        publishDir './results_SE/blacklist_filt_bam', mode: 'copy', pattern: '*.bam'    
    }

    input:

    tuple path(bam), path(index)
    path(blacklist_file)

    output:
    path("${out_bl_filtered_bam}"), emit: bl_filtered_bams


    script:
    out_bl_filtered_bam = "${bam.baseName}_BL_filt.bam"

    """
    ############ using bedtools intersect parameters ##############
    # -a : takes a bam/bed/gff/vcf file A
    # -b : takes one or more of the files bam/bed/gff/vcf and calls it B
    # -v : if any regions in A do not have a overlap in B then keep only those regions. that will remove the blacklist regions found in the blacklist B file

    #####################################################

    bedtools intersect \
    -a "${bam}" \
    -b "${blacklist_file}" \
    -v \
    > "${out_bl_filtered_bam}"


    """

}

process samtools_bl_index {
    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/samtools_rj'

    

    //publishDir './blacklist_filt_bam/bl_filt_index', mode: 'copy', pattern:'*.bai'
    if (params.PE) {

        if (params.ATAC) {

            publishDir './results_PE/ATAC_blacklist_filt_bam/bl_filt_index', mode: 'copy', pattern: '*.bai'
            publishDir './results_PE/ATAC_blacklist_filt_bam', mode: 'copy', pattern: '*_sort2.bam'

        }
        else {

            publishDir './results_PE/blacklist_filt_bam/bl_filt_index', mode: 'copy', pattern: '*.bai'
            publishDir './results_PE/blacklist_filt_bam', mode: 'copy', pattern: '*_sort2.bam'
        }

       
    
    }
    else {

        if (params.ATAC) {

            publishDir './results_SE/ATAC_blacklist_filt_bam/bl_filt_index', mode: 'copy', pattern: '*.bai'
            publishDir './results_SE/ATAC_blacklist_filt_bam', mode: 'copy', pattern: '*_sort2.bam'

        }
        else {

            publishDir './results_SE/blacklist_filt_bam/bl_filt_index', mode: 'copy', pattern: '*.bai'
            publishDir './results_SE/blacklist_filt_bam', mode: 'copy', pattern: '*_sort2.bam'
        }

          
    }

    input:
    path(bl_filt_bam)


    output:

    tuple path("${out_bam_name_sort}"), path("*.bai"), emit: bl_filt_bam_index_tuple


    script:
    
    out_bam_name_sort = "${bl_filt_bam.baseName}_sort2.bam"

    """
    ####### parameters for indexing bam ######
    # -b : will create a bai file

    ##########################################

    # just do some sorting 

    samtools sort \
    -o "${out_bam_name_sort}" \
    -O bam \
    "${bl_filt_bam}"



    samtools index \
    -b \
    "${out_bam_name_sort}" 

    """
}

process fastp_PE {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/fastp_rj'

    publishDir './results_PE/fastp_pe_results/filt_fastqs', mode: 'copy', pattern: '*_filt_{R1,R2}*'

    publishDir './results_PE/fastp_pe_results/merged_filt_fastqs', mode: 'copy', pattern: '*_merged*'

    publishDir './results_PE/fastp_pe_results/failed_qc_reads', mode: 'copy', pattern: '*_failed_filter*'

    publishDir './results_PE/fastp_pe_results/htmls', mode: 'copy', pattern: '*fastp.html'


    input:

    tuple val(fastq_name), path(fastq)


    output:

    tuple val(fastq_name), path("${out_name_1}"), path("${out_name_2}"), emit: filt_PE_tuple
    //path("${merged_reads_file}"), emit: merged_filt_reads
    path("${html_file_name}"), emit: html_fastp_out
    path("${failed_reads_file}"), emit: failed_reads_out


    script:

    out_name_1 = "${fastq_name}_filt_R1_0.fastq"
    out_name_2 = "${fastq_name}_filt_R2_0.fastq"
    failed_reads_file = "${fastq_name}_failed_filter_reads.fastq"
    merged_reads_file = "${fastq_name}_merged_file_reads.fastq"
    html_file_name = "${fastq_name}_R1_R2_fastp.html"

    """
    #!/usr/bin/env bash

    # this below is just debugging
    echo "this is the file name: ${fastq_name}; this is the forward read (r1): ${fastq[0]}; this is the reverse read(r2): ${fastq[1]}"


    ############# Now the PE fastp parameters ################
    # --in1 or -i : takes the first pair end read
    # --in2 or -I : takes the second pair end read
    # --out1 or -o : is the output file name of the fastp first read
    # --out2 or -O : is the output file name of the fastp second read
    # --failed_out : specify the file to store reads that cannot pass the filters
    # --merge or -m : for paired end input, merge each pair of reads into a single read if they are overlapped. the merged reads are written to --merged_out, the unmerged reads will be written to the --out1 and --out2
    # --detect_adapter_for_pe : to enable auto detection for PE data. The auto detection for adapter is for SE only, so have to turn on for PE
    # --dedup : enable deduplication to drop the duplicated reads pairs
    # --dup_calc_accuracy : accuracy level to calculate duplication (1~6) higher level uses more memory
    # --trim_poly_g : force polyG tail trimming. polyG can happen if there is no signal in the illumina two-color systems
    # --trim_poly_x : enable polyX trimming in 3' ends. useful because polyX(polyA) can be found in the tails of mRNA seq reads. DONT NEED THIS IN THIS PIPELINE
    # --qualified_quality_phred : this is the quality value that a base is qualified. default 15 but i choose 20 meaning phred quality >= 20
    # --unqualified_percent_limit : how many percent of bases are allowed to be unqualified. defualt is 40 meaning 40%
    # --n_base_limit : if one read's number of N base is > n ase imit  , then this read/pair is discarded default is 5. see fastp documentation
    # --average_qual :  if one read's average quality score < avg ual , then this read/pair is discarded. Default 0 means no requirement
    # --correction or -c : enable base correction in overlapped regions (only for PE data), default is disabled. fastp performs overlap analysis for PE data, which try to find an overlap of each pair of reads. When this option is enabled, and if an proper overlap is found, it can correct mismatched base pairs in overlapped regions of paired end reads, if one base is with high quality while the other is with ultra low quality. If a base is corrected, the quality of its paired base will be assigned to it so that they will share the same quality.
    # --overlap_len_require : the min length to detect oveerlapped region of PE reads. this will affect overlap analysis based PE merge, adapter trimming. default is 30 but i put 20 since the risca_lab snakemake pipeline looks for 20 
    # --overlap_diff_limit : the maximum number of mismatched bases to detect overlapped region of PE reads. it will affect any trimming or merge parameters. the default is 5 but i think the risc_lab snakemake pipeline allowed only for 1 mismatched bases. the snake make pipeline does only allow for 1 mismatches. But the author said "for now" so I assumed they will look to change it in the future. maybe i can change it back to the default which is 5
    # --overrepresentation_analysis : enable overrepresented sequence analysis
    # --overrepresentation_sampling : the number of reads computed for overrepresentation analysis (1~10000). default 20 i used 30
    # --html or -h : the html format report file name. default is fastp.html but i made my own name using the base name string (key) of the two fastq files that were input in this channel
    # --thread or -w : worker thread number default is 2; i used 15
    ###########################################################

    # NOTE: I will remove the merged reads and only keep all the reads that pass the filtering in their corresponding forward and reverse reads
    #--merge \
    #--merged_out "\${merged_reads_file}" 

    fastp \
    --in1 "${fastq[0]}" \
    --in2 "${fastq[1]}" \
    --out1 "${out_name_1}" \
    --out2 "${out_name_2}" \
    --failed_out "${failed_reads_file}" \
    --detect_adapter_for_pe \
    --dedup \
    --dup_calc_accuracy 5 \
    --trim_poly_g \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 40 \
    --n_base_limit 5 \
    --average_qual 0 \
    --correction \
    --overlap_len_require 20 \
    --overlap_diff_limit 1 \
    --overrepresentation_analysis \
    --overrepresentation_sampling 30 \
    --html "${html_file_name}" \
    --thread 15

    """

}

process fastqc_PE {

    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/fastqc_rj_env.yml'

    publishDir './results_PE/fastqc_pe_files', mode: 'copy', pattern: '*'

    input:
   
    tuple val(fastq_name), path(filt_r1), path(filt_r2)


    output:

    path ("*.html"), emit: fastqc_htmls
    path("*.zip"), emit: fastqc_zip_files



    script:


    """
    #!/usr/bin/env bash

    ########### parameters for fastqc pe reads ###############
    # no real parameters 

    ###########################################################

    fastqc \
    "${filt_r1}" 

    fastqc \
    "${filt_r2}" 




    """


}


process multiqc_PE {

    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/multiqc_rj_env.yml'

    publishDir './results_PE/multiqc_PE_output', mode: 'copy', pattern: '*'


    input:
    path(fastqc_zip_files)


    output:

    path("*"), emit: summary_of_PE_filt
    


    script:


    """
    #!/usr/bin/env bash

    ######### parameters for making the html file with all the fastqc output data #############

    
    # --ai : generate an AI summary of the report
    # --ai-summary-full : generate a detailed AI summary of the report
    # --interactive : use only interactive plots
    # --profile-runtime : add an analysis of how long multiqc takes to run to the report
    # --ai-provider : choose ai provider defualt sequra . openai or anthropic
    ############################################################################

    multiqc . \
    --ai-summary-full \
    --ai-provider openai \
    --interactive \
    --profile-runtime \
    --title "Pair end QC"


    """
}

process bwa_PE_aln {

    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/bwa_rj_env.yml'

    publishDir './results_PE/pe_bwa_files/pe_sam_files', mode: 'copy', pattern: '*.sam'
    publishDir './results_PE/pe_bwa_files/pe_sai_index_files', mode: 'copy', pattern: '*.sai'
    //cache false 


    input:
    tuple val(filt_fastq_name), path(fastq_r1), path(fastq_r2)
    path(genome)
    path(genome_index)


    output:

    path("*.sam"), emit: pe_sam_files
    path("*.sai"), emit: pe_sai_files



    script:

    sai_out_file_r1 = "${filt_fastq_name}_filt_r1.sai"
    sai_out_file_r2 = "${filt_fastq_name}_filt_r2.sai"

    out_sam_file = "${filt_fastq_name}_filt_r1_r2.sam"


    """
    #!/usr/bin/env bash

    ######## bwa aln parameters / bwa sampe params #########
    # -t : allows for the amout of threads you want this process to use

    #



    #########################################################


    bwa aln \
    -t 20 \
    "${genome}" \
    "${fastq_r1}" \
    > "${sai_out_file_r1}"

    bwa aln \
    -t 20 \
    "${genome}" \
    "${fastq_r2}" \
    > "${sai_out_file_r2}"


    bwa sampe \
    "${genome}" \
    "${sai_out_file_r1}" \
    "${sai_out_file_r2}" \
    "${fastq_r1}" \
    "${fastq_r2}" \
    > "${out_sam_file}"

    """
}

process multiqc_bam_stats {

    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/multiqc_rj_env.yml'

    if (params.PE) {

        publishDir './results_PE/flag_stat_log/complete_log', mode: 'copy', pattern: '*.html'
    
    }
    else {

        publishDir './results_SE/flag_sta_log/complete_log', mode: 'copy', pattern: '*.html'    
    }

    input:

    path(stats_log_files)
    path(norm_stats_files)


    output:
    path("*.html"), emit: bams_multiqc_html


    script:

    if (params.PE) {

        """

        #### parameters for multiqc ###

        # no real parameters
        ###############################
        
        multiqc . \
        --title "Pair end bams QC"



        """

    }
    else if (params.SE) {

        """

        #### parameters for multiqc ###

        # no real parameters
        ###############################
        
        multiqc . \
        --title "Single end bams QC"



        """

    }
    
}

process deeptools_aln_shift {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/deeptools_rj'

    if (params.PE) {

        publishDir './results_PE/atac_shift_bam', mode: 'copy', pattern: '*shift.bam'
    }
    else {

        publishDir './results_SE/atac_shift_bam', mode: 'copy', pattern: '*shift.bam'
    }


    input:
    tuple path(bam), path(index)


    output:
    tuple path("*shift.bam"), emit: atac_shifted_bam


    script:
    out_file = "${bam.baseName}_ATAC_shift.bam"

    """
    #!/usr/bin/env bash 

    ##### deeptools alignmentSieve params #####
    # --bam : this takes your input bam file
    # --numberofProcessors : number of threads 
    # --ATACshift : Shift the produced BAM file or BEDPE regions as commonly done for ATAC-seq
    ###########################################

    alignmentSieve \
    --bam "${bam}" \
    --numberOfProcessors 12 \
    --ATACshift \
    -o "${out_file}"




    """
}


// finally using modules

include {samtools_index_sort} from './modules/fastq2bam_dna_modules.nf'



// i want to include another workflow script

include {breakDensityWrapper_workflow} from './workflows/breakDensityWrapper_workflow.nf'

workflow {

    // this is the end seq alignment steps first


    // i will use a path already in the hpc as the defualt human genome but the user can change the genome by using -genome parameter and putting the path to a new genome in the command line when calling nextflow run
    params.genome = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg19/genome/Sequence/Bowtie2Index/genome.fa')

    // putting the human genome in a channel
    // keeping the human genome in a value channel so i can have other processes run more than once.
    genome_ch = Channel.value(params.genome)

    params.blacklist_path = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg19/blacklist/hg19-blacklist.v2.bed')
                
    blacklist_ch = Channel.value(params.blacklist_path)

    // i want to add an if then logic to the pipeline so i know which type of reads are comming in paired end or single end

    //if ( params.PE )
         // lets get the channel for the reads first

      //   fastp_PE()
         
         //align_PE_reads(genome_ch)
    
    
    
    
    if ( params.SE ) {

        
            

            

            // lets get the channel for the single end reads first
            // only use the single end read 1 data from the end seq which are already stored here: /rugpfs/fs0/risc_lab/store/hcanaj/HC_ENDseq_Novaseq_010925/read1_fastqs

        

            params.single_end_reads = file('/rugpfs/fs0/risc_lab/store/hcanaj/HC_ENDseq_Novaseq_010925/read1_fastqs/*_1.fastq.gz')
            se_reads_files = Channel.fromPath(params.single_end_reads)
            
            // now let's get the basename of the single end reads
            // removing both the .gz and the .fastq
            // I would normally use file.baseName here but it had the .gz and the .fastq
            se_reads_files.flatten()
                            .map{ file -> file.name.replace('.fastq.gz','')}
                            .set{se_reads_name}
            
            // let's view both the files and the names to make sure they match in order
            //se_reads_files.view()
            //se_reads_name.view()
            // this is where i send both the input file and their corresponding basenames to the fastp_SE process
            

            // if the adapter sequence is known then input it as a string if not dont use the parameter
            if ( params.ada_seq ) {

                params.adapter_seq_str = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA' // this is just a place holder value for the adapter sequence
                adapter_ch = Channel.value(params.adapter_seq_str)

                fastp_SE_adapter_known(se_reads_files.take(3), se_reads_name.take(3), adapter_ch) // will have to make a new process for if the adapter sequence is known

                fastq_filts = fastp_SE_adapter_known.out.filtered_fastqs
                //fastp_SE.out.view()
                fastq_filts.map{file -> file.baseName}
                            .set{fastq_filts_name}

                // now getting the html files since i think fastqc combines them into one, that might be multiqc
                fastp_filt_html = fastp_SE_adapter_known.out.fastp_html_reports

            }    
            else {

                fastp_SE(se_reads_files.take(3), se_reads_name.take(3))

                    // take all of the filtered fastq files and put them in a channel name
                // since the fastq files might be in a different order, if i need to get their base names I will have to do it from this new channel below
                fastq_filts = fastp_SE.out.filtered_fastqs
                //fastp_SE.out.view()
                fastq_filts.map{file -> file.baseName}
                            .set{fastq_filts_name}

                // now getting the html files since i think fastqc combines them into one, that might be multiqc
                fastp_filt_html = fastp_SE.out.fastp_html_reports



            }
            

            //fastp_SE(se_reads_files.take(3), se_reads_name.take(3)) // REMEMBER TO REMOVE THIS TESTING FEATURE WHERE IT WILL ONLY TAKE THE FIRST 3

            

            
            //fastp_SE.out.view()

            //fastq_filts.view()
            //fastp_filt_html.view()
            //fastp_filt_html.collect().view()


            // now creating a fastqc process
            fastqc_SE(fastq_filts, fastq_filts_name)

            fastqc_html_files = fastqc_SE.out.fastqc_htmls
            fastqc_zips = fastqc_SE.out.fastqc_zip_files

            // now using multiqc to combine all of the se zip files. Multiqc takes the zip files generated by fastqc and puts them in a single html file
            multiqc_SE(fastqc_zips.collect())

            // first have a seprate process that indexes the reference genome using bwa or bwa mem so this part doesnt have to be done again and will be cached
            bwa_index_genome(genome_ch)

            // collecting the genome index files from the last process 
            // not sure if i should keep track of the order the files are in first
            // it looks like they are in the same order that they appeared in the published dir using ll
            //bwa_index_genome.out.genome_index_files.view()
            genome_index_files_ch = bwa_index_genome.out.genome_index_files

            // Now I need to pass the human genome file to the process to index the genome file. Also I will add the filtered fastq files from fastp into this process that will be aligned to the genome. 
            // each run of this only takes 20-30 min to run but since the hpc only is allowing 2-3 to run at one time it takes 3 hours
            bwa_align_SE(genome_ch, genome_index_files_ch, fastq_filts, fastq_filts_name )

            //bwa_align_SE.out.sam_se_files.view()

            // making a channel for the sam files generated
            sam_files = bwa_align_SE.out.sam_se_files

            // now I want to take any sam files generated by the bwa and use samtools to order them and convert them into bam files
            // I will hopefully be able to do this outside of the if else statement so the sam file from both conditions can be passed to the same samtools process
            // since i am just testing the pipeline i should find a way to do this on only a few files (about 3-4)
            
            samtools_sort(sam_files) // using take 3 should only take the first 3 files from the sam_files channel which should have 64 sam files. This is just for production and testing. will remove when running pipeline for real.

            //samtools_sort.out.sorted_bams.view()
            //samtools_sort.out.indexed_bams.view()
            //samtools_sort.out.bam_index_tuple.view()

            sorted_bams_ch = samtools_sort.out.sorted_bams
            indexed_bams_ch = samtools_sort.out.indexed_bams
            both_bam_index_ch = samtools_sort.out.bam_index_tuple
            flagstat_log_ch = samtools_sort.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
            norm_stats_txt_ch = samtools_sort.out.norm_stats_txt.collect()
            tsv_SN_stats_ch = samtools_sort.out.tsv_SN_stats.collect()
            // if you want to filter black list use the param --BL in the command line when calling nextflow
            if ( params.BL ) {

                // using bedtools to filter black list but first giving the user an option to put the correct black list for an organism
                // by defualt it will use the hg19 v2 blacklist, but if you used a different organism or human genome use the appropriate blacklist
                //params.blacklist_path = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg19/blacklist/hg19-blacklist.v2.bed')
                
                //blacklist_ch = Channel.value(params.blacklist_path)

                bedtools_filt_blacklist(both_bam_index_ch, blacklist_ch)

                bl_filt_bams_ch = bedtools_filt_blacklist.out.bl_filtered_bams
                // i will need to index the black list filtered bam again so i have to create a different samtools process for this
                samtools_bl_index(bl_filt_bams_ch)

                bl_filt_bam_tuple_ch = samtools_bl_index.out.bl_filt_bam_index_tuple

                if ( params.ATAC ) {


                    // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

                    deeptools_aln_shift(bl_filt_bam_tuple_ch)

                    atac_shift_bam_ch = deeptools_aln_shift.out.atac_shifted_bam
                    atac_shift_bam_ch.view()

                    // now i have to re index this new atac shifted bam. dispite the name of the process I can just pass any future created bam to this channel to be indexed

                    samtools_index_sort(atac_shift_bam_ch )

                    // so now name the tuple channel output appropriately 
                    atac_shift_bam_index_ch = samtools_index_sort.out.bl_filt_bam_index_tuple

                    // now making the bed files for atac seq

                    deeptools_make_bed(atac_shift_bam_index_ch)

                    deeptools_make_bed.out.bed_files_normalized.view()

                    bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized



                }
                else {



                    // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

                    deeptools_make_bed(bl_filt_bam_tuple_ch)

                    deeptools_make_bed.out.bed_files_normalized.view()

                    bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized          


                }

                // then i need to pass the indexed_bl_bam and the bam to the deeptools process

                //deeptools_make_bed(bl_filt_bam_tuple_ch)
                //bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized

            }

            else {


                if (params.ATAC) {


                    // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

                    deeptools_aln_shift(both_bam_index_ch)

                    atac_shift_bam_ch = deeptools_aln_shift.out.atac_shifted_bam

                    // now i have to re index this new atac shifted bam. dispite the name of the process I can just pass any future created bam to this channel to be indexed

                    samtools_index_sort(atac_shift_bam_ch)

                    // so now name the tuple channel output appropriately 
                    atac_shift_bam_index_ch = samtools_index_sort.out.bl_filt_bam_index_tuple

                    // now making the bed files for atac seq

                    deeptools_make_bed(atac_shift_bam_index_ch)

                    deeptools_make_bed.out.bed_files_normalized.view()

                    bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized



                }
                else {



                    // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

                    deeptools_make_bed(both_bam_index_ch)

                    deeptools_make_bed.out.bed_files_normalized.view()

                    bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized          


                }

                // Now i want to pass the tuple that has the bam and its corresponding index file into a process that will create a bigwig file for visulization, created from the bam file. This will show read coverage in the genome without looking for significant areas
                //deeptools_make_bed(both_bam_index_ch)

                //deeptools_make_bed.out.bed_files_normalized.view()
                //bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized
            }
    

            // i will either use the bams or the bed files for any future processes depending on what tool needs what.
    }

    // now I need to make the paired-end part of this pipeline.
    else if (params.PE) {

        // this will take the paired end reads and keep them together
        params.paired_end_reads = '/rugpfs/fs0/risc_lab/store/hcanaj/HC_GLOEseq_Novaseq_010925/fastqs_read1_read2/*_{R1,R2}*'

        pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads)

        //pe_fastqs_ch.view()

        fastp_PE(pe_fastqs_ch.take(3))

        // checking the channels to see if everything works
        //fastp_PE.out.filt_PE_tuple.view()

        //fastp_PE.out.html_fastp_out.view()
        //fastp_PE.out.failed_reads_out.view()
        //fastp_PE.out.

        pe_filt_tuple_ch = fastp_PE.out.filt_PE_tuple

        fastqc_PE(pe_filt_tuple_ch)

        //fastqc_PE.out.fastqc_zip_files.view()


        // then collect them so i can view in one file using multiqc.
        collection_fastqc_ch =fastqc_PE.out.fastqc_zip_files.collect()

        multiqc_PE(collection_fastqc_ch)

        //multiqc_PE.out.summary_of_PE_filt.view()

        // first use the process to index the reference genome since the process exists already for the se

        bwa_index_genome(genome_ch)

        //bwa_index_genome.out.genome_index_files.view()

        genome_index_ch = bwa_index_genome.out.genome_index_files

        // now to use bwa aln and bwa sampe to align the filtered pair end reads to the reference genome and then to create the sam file respectively

        //pe_filt_tuple_ch.view()
        //genome_ch.view()
        //genome_index_ch.view()

        bwa_PE_aln(pe_filt_tuple_ch, genome_ch, genome_index_ch)

        // now check to see if the output channels are good
        //bwa_PE_aln.out.pe_sam_files.view()
        
        
        // now i need to make the parameters for  if the bam file will be blacklist filtered or not
        
        // using this channel for both if Blacklist or not
        sam_files_pe_ch = bwa_PE_aln.out.pe_sam_files

        if (params.BL) {

            
            // i have to make a bam file to then use bedtools intersect to get the blacklist
            // using the same samtools sort process found in the SE part of the pipeline
            samtools_sort(sam_files_pe_ch)

            //samtools_sort.out.bam_index_tuple.view()

            bam_index_tuple_ch = samtools_sort.out.bam_index_tuple
            flagstat_log_ch = samtools_sort.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
            norm_stats_txt_ch = samtools_sort.out.norm_stats_txt.collect()
            tsv_SN_stats_ch = samtools_sort.out.tsv_SN_stats.collect()

            // this will give a blacklist filtered bam but i need to index it again
            bedtools_filt_blacklist(bam_index_tuple_ch, blacklist_ch)
            //bedtools_filt_blacklist.out.bl_filtered_bams.view()

            bl_filt_bams_ch = bedtools_filt_blacklist.out.bl_filtered_bams
            // so using the process to only index which means it will take the blacklist bam file
            
            samtools_bl_index(bl_filt_bams_ch) 

            bl_filt_bam_tuple_ch = samtools_bl_index.out.bl_filt_bam_index_tuple

            // so this would give a bam that is bl filtered and has an index

            if ( params.ATAC ) {


                // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

                deeptools_aln_shift(bl_filt_bam_tuple_ch)

                atac_shift_bam_ch = deeptools_aln_shift.out.atac_shifted_bam
                //atac_shift_bam_ch.view()

                // now i have to re index this new atac shifted bam. dispite the name of the process I can just pass any future created bam to this channel to be indexed

                samtools_index_sort(atac_shift_bam_ch )

                // so now name the tuple channel output appropriately 
                atac_shift_bam_index_ch = samtools_index_sort.out.bl_filt_bam_index_tuple

                // now making the bed files for atac seq

                deeptools_make_bed(atac_shift_bam_index_ch)

                //deeptools_make_bed.out.bed_files_normalized.view()

                bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized



            }
            else {



                // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

                deeptools_make_bed(bl_filt_bam_tuple_ch)

                //deeptools_make_bed.out.bed_files_normalized.view()

                bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized          


            }
            
            
            // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

            //deeptools_make_bed(bl_filt_bam_tuple_ch)

            //deeptools_make_bed.out.bed_files_normalized.view()

            //bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized

            

        }
        else {
            samtools_sort(sam_files_pe_ch)

            //samtools_sort.out.bam_index_tuple.view()

            flagstat_log_ch = samtools_sort.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
            norm_stats_txt_ch = samtools_sort.out.norm_stats_txt.collect()
            tsv_SN_stats_ch = samtools_sort.out.tsv_SN_stats.collect()

            bam_index_tuple_ch = samtools_sort.out.bam_index_tuple

            // add the if logic for ATAC here

            if (params.ATAC) {


                // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

                deeptools_aln_shift(bam_index_tuple_ch)

                atac_shift_bam_ch = deeptools_aln_shift.out.atac_shifted_bam

                // now i have to re index this new atac shifted bam. dispite the name of the process I can just pass any future created bam to this channel to be indexed

                samtools_index_sort(atac_shift_bam_ch)

                // so now name the tuple channel output appropriately 
                atac_shift_bam_index_ch = samtools_index_sort.out.bl_filt_bam_index_tuple

                // now making the bed files for atac seq

                deeptools_make_bed(atac_shift_bam_index_ch)

                //deeptools_make_bed.out.bed_files_normalized.view()

                bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized



            }
            else {



                // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

                deeptools_make_bed(bam_index_tuple_ch)

                //deeptools_make_bed.out.bed_files_normalized.view()

                bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized          


            }






            // just using the original sorted and indexed bam
            
            // i will comment this out below since it was the original but the above if logic should work if there is ATAC data or if there is not

            //deeptools_make_bed(bam_index_tuple_ch)

            //deeptools_make_bed.out.bed_files_normalized.view()

            //bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized


        }

        

 
    }

    // making a multiqc process for the samtools flagstat log files. this should be able to take the flagstat_log_ch from any part of the choosen paths
    multiqc_bam_stats(flagstat_log_ch, norm_stats_txt_ch)


    
    // I want to have a parameter that takes peakfiles. The default will be the IMR90 narrowPeak files
    params.peak_files = files('/lustre/fs4/home/ascortea/store/ascortea/beds/IMR90/*.narrowPeak')
    //now making the channel for the files
    peak_files_ch = Channel.value(params.peak_files)


    if (params.calc_break_density){
    // i want to call the workflow breakDensityWrapper_workflow and pass the bam_index_tuple_ch as an input from either path where the user chose to do blacklist filter or not. Then also pass the peak files that already exists or are created later as input
    
        if (params.PE) {

            if (params.BL) {

                breakDensityWrapper_workflow(bl_filt_bam_tuple_ch, peak_files_ch)

            }
            else {
                breakDensityWrapper_workflow(bam_index_tuple_ch, peak_files_ch)

            }
            

        }

        if (params.SE) {

            if (params.BL) {

                breakDensityWrapper_workflow(bl_filt_bam_tuple_ch, peak_files_ch)

            }
            else {

                breakDensityWrapper_workflow(both_bam_index_ch, peak_files_ch)

            }

            
        }
    }

    // I want to make a log file with all the stats from using samtools stats on each bam file

    if (params.PE) {

        //(tsv_SN_stats_ch)
    }
    if (params.SE) {

        //py_calc_stats_log(tsv_SN_stats_ch)

    }
    
}