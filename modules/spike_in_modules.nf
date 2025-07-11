
process fastp_SE_adapter_known_spike_in {

    label 'spike_in_big_resoruces'

    // using this conda yml file that I created. Nextflow will now make its own conda environment from the dependencies found within.
    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/fastp_rj_env.yml'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/fastp_rj'

    
    // publishDir "./results_SE/fastp_qc_single_end/spike_in${spike_name}", mode: 'copy', pattern:'*_fp_filt.fastq'
    // publishDir "./results_SE/fastp_qc_single_end/html_reports/spike_in${spike_name}", mode: 'copy', pattern:'*.html'

    input:
    // input names dont have to be the exact same as what is seen in the workflow section of this script
    // as long and they are in the same order, you will have the correct input
    path(fastq_files) // the files
    val(fastq_names) // the names of the files as a value input channel
    val(adapter_seq)
    tuple val(spike_name), path(spike_genome)

    publishDir "${params.base_out_dir}/fastp_qc_single_end/spike_in${spike_name}", mode: 'copy', pattern:'*_fp_filt.fastq'
    publishDir "${params.base_out_dir}/fastp_qc_single_end/html_reports/spike_in${spike_name}", mode: 'copy', pattern:'*.html'

    output:
    path("${out_name}"), emit: filtered_fastqs
    path("${fastq_names}*.html"), emit: fastp_html_reports


    script:

    // getting the output name
    out_name = "${fastq_names}${spike_name}_fp_filt.fastq"


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
process fastp_SE_spike_in {

    label 'spike_in_big_resoruces'

    // using this conda yml file that I created. Nextflow will now make its own conda environment from the dependencies found within.
    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/fastp_rj_env.yml'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/fastp_rj'

    // publishDir "./results_SE/fastp_qc_single_end/spike_in${spike_name}", mode: 'copy', pattern:'*_fp_filt.fastq'
    // publishDir "./results_SE/fastp_qc_single_end/html_reports/spike_in${spike_name}", mode: 'copy', pattern:'*.html'

    input:
    // input names dont have to be the exact same as what is seen in the workflow section of this script
    // as long and they are in the same order, you will have the correct input
    path(fastq_files) // the files
    val(fastq_names) // the names of the files as a value input channel
    tuple val(spike_name), path(spike_genome)

    publishDir "${params.base_out_dir}/fastp_qc_single_end/spike_in${spike_name}", mode: 'copy', pattern:'*_fp_filt.fastq'
    publishDir "${params.base_out_dir}/fastp_qc_single_end/html_reports/spike_in${spike_name}", mode: 'copy', pattern:'*.html'

    output:
    path("${out_name}"), emit: filtered_fastqs
    path("${fastq_names}*.html"), emit: fastp_html_reports


    script:

    // getting the output name
    out_name = "${fastq_names}${spike_name}_fp_filt.fastq"


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
process fastqc_SE_spike_in {

    label 'spike_in_small_resources'

    // using the conda environment 
    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/fastqc_rj_env.yml'
    
    // publishDir "./results_SE/fastqc_htmls/spike_in${spike_name}", mode: 'copy', pattern: '*.html'


    input:

    path(fastq_filt_files)
    val(fastq_filt_names)
    tuple val(spike_name), path(spike_genome)

    publishDir "${params.base_out_dir}/fastqc_htmls/spike_in${spike_name}", mode: 'copy', pattern: '*.html'

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

process multiqc_SE_spike_in {

    label 'spike_in_small_resources'

    // this yml file doesnt work
    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/multiqc_rj_env.yml'

    //conda '/ru-auth/local/home/rjohnson/miniconda3/envs/multiqc_rj'
    
    // publishDir "./results_SE/multiQC_collection/spike_in${spike_name}", mode: 'copy', pattern: '*.html'

    input:
    path(fastp_filt_html)
    tuple val(spike_name), path(spike_genome)

    publishDir "${params.base_out_dir}/multiQC_collection/spike_in${spike_name}", mode: 'copy', pattern: '*.html'


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



process bwa_index_genome_spike_in {

    label 'spike_in_big_resoruces'

    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/bwa_rj_env.yml'

    // publishDir "./genome_index_bwa/spike_in/genome_index${spike_name}", mode: 'copy', pattern: '*'

    // dont need the if then statement anymore
    
    


    input:
    //path(ref_genome)
    tuple val(spike_name), path(ref_genome)

    publishDir "./genome_index_bwa/spike_in/genome_index${spike_name}", mode: 'copy', pattern: '*'


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

process bwa_align_SE_spike_in {

    label 'spike_in_big_resoruces'

    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/bwa_rj_env.yml'

    // // i think i can do the same here now without having too many if conditions
    // publishDir "./results_SE/bwa_outputs_singleEnd_SAM/spike_in/sam${spike_name}", mode: 'copy', pattern: '*.sam'
    // publishDir "./results_SE/sai_alignment_files/spike_in/sai${spike_name}", mode: 'copy', pattern: '*.sai'

    //publishDir './results_SE/bwa_outputs_singleEnd_SAM/spike_in', mode: 'copy', pattern: '*.sam'
    //publishDir './results_SE/sai_alignment_files/spike_in', mode: 'copy', pattern: '*.sai'

    // a user running the single end path will only be using end_seq spike ins, not sure if ricc-seq will be pair end or single end
    
    

    input:
    // path(ref_genome)
    path(genome_index_files)
    path(fastq_filt_files)
    val(fastq_filt_names)
    tuple val(spike_name), path(spike_genome)

    // i think i can do the same here now without having too many if conditions
    publishDir "${params.base_out_dir}/bwa_outputs_singleEnd_SAM/spike_in/sam${spike_name}", mode: 'copy', pattern: '*.sam'
    publishDir "${params.base_out_dir}/sai_alignment_files/spike_in/sai${spike_name}", mode: 'copy', pattern: '*.sai'

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
    "${spike_genome}"  \
    "${fastq_filt_files}" \
    > "${sai_output_file}"


    bwa samse \
    "${spike_genome}" \
    "${sai_output_file}" \
    "${fastq_filt_files}" \
    > "${sam_name}"


    """
}




process samtools_sort_spike_in {

    label 'spike_in_big_resoruces'
    // using the conda yml file for samtools
    // it doesnt work
    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/samtools_rj_env.yml'

    //conda '/ru-auth/local/home/rjohnson/miniconda3/envs/samtools_rj' // this was samtools version 1.3 which doesnt have samtools fixmate option -m

    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/samtools-1.21_spec_env_rj.txt' // that is the explicit file but if that doesnt work try the yml file samtools-1.21_env_rj.yml; and if that doesnt work use the path to the 1.21 environment

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/samtools-1.21_rj'

    input:
    path(sam_files)
    tuple val(spike_name), path(spike_genome)

    publishDir "${params.base_out_dir}/sorted_bam_files/spike_in${spike_name}", mode: 'copy', pattern: '*_sorted.bam'
    publishDir "${params.base_out_dir}/sorted_bam_files/spike_in${spike_name}", mode: 'copy', pattern: '*.{bai, csi}'
    // publishDir "${params.base_out_dir}/flag_stat_log/spike_in${spike_name}", mode: 'copy', pattern: '*stats.log'
    // publishDir "${params.base_out_dir}/stats_tsv_files/spike_in${spike_name}", mode: 'copy', pattern: '*stats.tsv'

    // Determine the base directory (PE/SE) first
    // if (params.PE) {
    //     // Handle PE with spike-ins
        
    //     publishDir "./results_PE/sorted_bam_files/spike_in${spike_name}", mode: 'copy', pattern: '*_sorted.bam'
    //     publishDir "./results_PE/sorted_bam_files/spike_in${spike_name}", mode: 'copy', pattern: '*.{bai, csi}'
    //     publishDir "./results_PE/flag_stat_log/spike_in${spike_name}", mode: 'copy', pattern: '*stats.log'
    //     publishDir "./results_PE/stats_tsv_files/spike_in${spike_name}", mode: 'copy', pattern: '*stats.tsv'
    
        
    // } 
    // if (params.SE) {

    //     // Handle SE with spike-ins
        
    //     publishDir "./results_SE/sorted_bam_files/spike_in${spike_name}", mode: 'copy', pattern: '*_sorted.bam'
    //     publishDir "./results_SE/sorted_bam_files/spike_in${spike_name}", mode: 'copy', pattern: '*.{bai, csi}'
    //     publishDir "./results_SE/flag_stat_log/spike_in${spike_name}", mode: 'copy', pattern: '*stats.log'
    //     publishDir "./results_SE/stats_tsv_files/spike_in${spike_name}", mode: 'copy', pattern: '*stats.tsv'
        
    // }
    // //publishDir './sorted_bam_files', mode: 'copy', pattern: '*_sorted.bam'
    //publishDir './indexed_bam_files', mode: 'copy', pattern: '*.{bai, csi}'
    //publishDir './flag_stat_log', mode: 'copy', pattern: '*stat.log'

    // input:
    // path(sam_files)
    // tuple val(spike_name), path(spike_genome)

    output:

    path("*_sorted.bam"), emit: sorted_bams
    //tuple path("*.{bai,csi}"), emit: indexed_bams
    //tuple path("*.bai"), path("*.csi"), emit: indexed_bams
    path("*.bai"), emit: indexed_bams
    tuple path("*_sorted.bam"), path("*.bai"), emit: bam_index_tuple

    // path("*stats.log"), emit: flag_stats_log
    // path("*stats.txt"), emit: norm_stats_txt
    // path("*stats.tsv"), emit: tsv_SN_stats

    script:

    // i will start using baseName inside the process since its easier to keep track of different names an uses less inputs into a process
    out_bam_filt = "${sam_files.baseName}_bam_filt.bam"
    out_bam_name_sort = "${sam_files.baseName}_name_ordered.bam"
    out_bam_coor_sort = "${sam_files.baseName}_filt_coor_sorted.bam"
    out_bam_fixmate = "${sam_files.baseName}_fixmate.bam"
    // out_bam_final = "${sam_files.baseName}_markdup_filt_coor_sorted.bam"
    // flagstats_log = "${sam_files.baseName}_flag_stats.log"
    // samtools_stats_log = "${sam_files.baseName}_stats.txt"
    // tsv_file_with_stats = "${sam_files.baseName}_SN_stats.tsv"


    

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

    #samtools flagstat \
    #"\${out_bam_coor_sort}" \
    #> "\${flagstats_log}"

    # adding another way to get stats from each bam file
    #samtools stats \
    #"\${out_bam_coor_sort}" \
    #> "\${samtools_stats_log}"

    # now only putting the stats into a tsv file
    #less "\${samtools_stats_log}" | grep ^SN | cut -f 2-3 >  "\${tsv_file_with_stats}"
    """





    
}


// might not need this next process here
/*
process bam_log_calc_spike_in {

    label 'normal_small_resources'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/samtools-1.21_rj'

    publishDir "${params.base_out_dir}/flag_stat_log/spike_in${spike_name}", mode: 'copy', pattern: '*stats.log'
    publishDir "${params.base_out_dir}/stats_tsv_files/spike_in${spike_name}", mode: 'copy', pattern: '*stats.tsv'


    input:
    tuple path(bam), path(index_bam)


    output:
    path("*stats.log"), emit: flag_stats_log
    path("*stats.txt"), emit: norm_stats_txt
    path("*stats.tsv"), emit: tsv_SN_stats



    script:
    flagstats_log = "${bam.baseName}_flag_stats.log"
    samtools_stats_log = "${bam.baseName}_stats.txt"
    tsv_file_with_stats = "${bam.baseName}_SN_stats.tsv"


    """
    #!/usr/bin/env bash

    # using samtools flagstat: generate log files so i can use multiqc to get stats of all files into one html file
    # no parameters needed. just need to give the final bam file that went through all the processing


    samtools flagstat \
    "${bam}" \
    > "${flagstats_log}"

    # adding another way to get stats from each bam file
    samtools stats \
    "${bam}" \
    > "${samtools_stats_log}"

    # now only putting the stats into a tsv file
    less "${samtools_stats_log}" | grep ^SN | cut -f 2-3 >  "${tsv_file_with_stats}"




    """
}
*/



process deeptools_make_bed_spike_in {
    
    label 'spike_in_small_resources'

    // this conda env yml file didnt work. have to use the actual env
    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/deeptools_rj_env.yml'
    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/deeptools_rj'

    // this section is just a simple if else statement controlling the directories that are created and when the files end up
    // will copy and paste in other processes that need it.

    // the logic should be: if reads are pair end, go to the results_PE dir and in the spike_in dir make it the spike_name based on what spike in type called this process
    




    input:

    tuple path(bams), path(index)
    tuple val(spike_name), path(spike_genome)

    
    publishDir "${params.base_out_dir}/no_bl_filt/bed_graphs_deeptools/spike_in${spike_name}", mode: 'copy', pattern: '*'

    // if (params.PE) {
    //     publishDir "./results_PE/no_bl_filt/bed_graphs_deeptools/spike_in${spike_name}", mode: 'copy', pattern: '*'
    // }
    // if (params.SE) {
    //     publishDir "./results_SE/no_bl_filt/bed_graphs_deeptools/spike_in${spike_name}", mode: 'copy', pattern: '*'
    // }

    output:
    path("${out_bed_name}*"), emit: bed_files_normalized

    script:

    out_bed_name="${bams.baseName}${spike_name}_normalized_cpm.bed"

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

// not doing black list filter for spike in. cant find blacklist for t7, yeast, or lambda
/*
process bedtools_filt_blacklist_spike_in {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/bedtools_rj'

    //publishDir './blacklist_filt_bam', mode: 'copy', pattern: '*.bam'

    // if (params.PE) {

    //     publishDir './results_PE/blacklist_filt_bam', mode: 'copy', pattern: '*.bam'
    
    // }
    // else {

    //     publishDir './results_SE/blacklist_filt_bam', mode: 'copy', pattern: '*.bam'    
    // }

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

}*/

/*
process samtools_bl_index_spike_in {
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
}*/

process fastp_PE_spike_in {

    label 'spike_in_big_resoruces'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/fastp_rj'

    input:

    tuple val(fastq_name), path(fastq)
    tuple val(spike_name), path(spike_genome)

    publishDir "${params.base_out_dir}/fastp_pe_results/filt_fastqs/spike_in${spike_name}", mode: 'copy', pattern: '*_filt_{R1,R2}*'

    publishDir "${params.base_out_dir}/fastp_pe_results/merged_filt_fastqs/spike_in${spike_name}", mode: 'copy', pattern: '*_merged*'

    publishDir "${params.base_out_dir}/fastp_pe_results/failed_qc_reads/spike_in${spike_name}", mode: 'copy', pattern: '*_failed_filter*'

    publishDir "${params.base_out_dir}/fastp_pe_results/htmls/spike_in${spike_name}", mode: 'copy', pattern: '*fastp.html'


    // input:

    // tuple val(fastq_name), path(fastq)
    // tuple val(spike_name), path(spike_genome)

    output:

    tuple val(fastq_name), path("${out_name_1}"), path("${out_name_2}"), emit: filt_PE_tuple
    //path("${merged_reads_file}"), emit: merged_filt_reads
    path("${html_file_name}"), emit: html_fastp_out
    path("${failed_reads_file}"), emit: failed_reads_out


    script:

    out_name_1 = "${fastq_name}${spike_name}_filt_R1_0.fastq"
    out_name_2 = "${fastq_name}${spike_name}_filt_R2_0.fastq"
    failed_reads_file = "${fastq_name}${spike_name}_failed_filter_reads.fastq"
    merged_reads_file = "${fastq_name}${spike_name}_merged_file_reads.fastq"
    html_file_name = "${fastq_name}${spike_name}_R1_R2_fastp.html"

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

process fastqc_PE_spike_in {

    label 'spike_in_small_resources'

    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/fastqc_rj_env.yml'

    // publishDir "./results_PE/fastqc_pe_files/spike_in${spike_name}", mode: 'copy', pattern: '*'

    input:
   
    tuple val(fastq_name), path(filt_r1), path(filt_r2)
    tuple val(spike_name), path(spike_genome)

    publishDir "${params.base_out_dir}/fastqc_pe_files/spike_in${spike_name}", mode: 'copy', pattern: '*'

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


process multiqc_PE_spike_in {

    label 'spike_in_small_resources'

    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/multiqc_rj_env.yml'

    // publishDir "./results_PE/multiqc_PE_output/spike_in${spike_name}", mode: 'copy', pattern: '*'


    input:
    path(fastqc_zip_files)
    tuple val(spike_name), path(spike_genome)

    publishDir "${params.base_out_dir}/multiqc_PE_output/spike_in${spike_name}", mode: 'copy', pattern: '*'

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

process bwa_PE_aln_spike_in {

    label 'spike_in_big_resoruces'

    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/bwa_rj_env.yml'

    //publishDir './results_PE/pe_bwa_files/pe_sam_files/spike_in', mode: 'copy', pattern: '*.sam'
    //publishDir './results_PE/pe_bwa_files/pe_sai_index_files/spike_in', mode: 'copy', pattern: '*.sai'
    //cache false 

    // publishDir "./results_PE/pe_bwa_files/pe_sam_files/spike_in/sam${spike_name}", mode: 'copy', pattern: '*.sam'
    // publishDir "./results_PE/pe_bwa_files/pe_sai_index_files/spike_in/sai${spike_name}", mode: 'copy', pattern: '*.sai'

    

    input:
    tuple val(filt_fastq_name), path(fastq_r1), path(fastq_r2)
    //path(genome)
    path(genome_index)
    tuple val(spike_name), path(genome)

    publishDir "${params.base_out_dir}/pe_bwa_files/pe_sam_files/spike_in/sam${spike_name}", mode: 'copy', pattern: '*.sam'
    publishDir "${params.base_out_dir}/pe_bwa_files/pe_sai_index_files/spike_in/sai${spike_name}", mode: 'copy', pattern: '*.sai'

    output:

    path("*.sam"), emit: pe_sam_files
    path("*.sai"), emit: pe_sai_files



    script:

    sai_out_file_r1 = "${filt_fastq_name}${spike_name}_filt_r1.sai"
    sai_out_file_r2 = "${filt_fastq_name}${spike_name}_filt_r2.sai"

    out_sam_file = "${filt_fastq_name}${spike_name}_filt_r1_r2.sam"


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

/*
process multiqc_bam_stats_spike_in {

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
    
}*/




// samtools_index_sort_spike_in copied from the original but adding a sub dir for gloe_seq, end_seq and ricc_seq

process samtools_index_sort_spike_in {

    label 'spike_in_big_resoruces'
    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/samtools_rj'

    //publishDir './blacklist_filt_bam/bl_filt_index', mode: 'copy', pattern:'*.bai'

    // not doing blacklist filter for spike_ins

    

    input:
    path(bam) // changed this to be just bam
    tuple val(spike_name), path(spike_genome)


     
    publishDir "${params.base_out_dir}/sorted_bam_files/ATAC_filt_bam/spike_in${spike_name}_bam", mode: 'copy', pattern: '*.{bam,bai}'

    // if (params.PE) {

    //     publishDir "./results_PE/ATAC_filt_bam/spike_in${spike_name}_bam", mode: 'copy', pattern: '*.{bam,bai}'

    // }
    // if (params.SE) {
    //     publishDir "./results_SE/ATAC_filt_bam/spike_in${spike_name}_bam", mode: 'copy', pattern: '*.{bam,bai}'
    // }

    output:

    tuple path("${out_bam_name_sort}"), path("*.bai"), emit: bam_index_tuple
    
    

    script:
    
    out_bam_name_sort = "${bam.baseName}${spike_name}_sort2.bam"

    """
    ####### parameters for indexing bam ######
    # -b : will create a bai file

    ##########################################

    # just do some sorting 

    samtools sort \
    -o "${out_bam_name_sort}" \
    -O bam \
    "${bam}"



    samtools index \
    -b \
    "${out_bam_name_sort}" 

    """
}

process deeptools_aln_shift_spike_in {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/deeptools_rj'

    // commenting the publish dir here out since the samtools_index_sort_spike_in will have the bam and it will be sorted anyway
    /*if (params.PE) {

        publishDir './results_PE/atac_shift_bam', mode: 'copy', pattern: '*shift.bam'
    }
    else {

        publishDir './results_SE/atac_shift_bam', mode: 'copy', pattern: '*shift.bam'
    }*/


    input:
    tuple path(bam), path(index)
    tuple val(spike_name), path(spike_genome)

    // not using this dir because a different dir will have the bam and it will be sorted.
    //publishDir = "${params.base_out_dir}/atac_shift_bam/spike_in${spike_name}", mode: 'copy', pattern: '*shift.bam'

    output:
    path("*shift.bam"), emit: atac_shifted_bam


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

// process overlap_window_spike_in {

//     // this process might have 945 instances  17 bedfiles times the number of peakfiles

//     label 'super_small_resources'

//     publishDir "${params.base_out_dir}/alignment_peak_overlap_qc/spike_in/${spike_in}", mode: 'copy', pattern:"*${condition}*.tsv"

//     //debug true

//     conda '/ru-auth/local/home/rjohnson/miniconda3/envs/bedtools_rj'

//     input:

//     // the order of the condition and everything else is 0GyP, PLC, 0GyC, Cells

//     tuple val(grouping_name), val(condition), val(spike_in), val(basename), path(bedfiles), path(peakfile)

//     //path(peaks)

//     // path(zero_gy)

//     // path(plc)

//     // path(cells)

//     // path(all_peaks)



//     output:
//     path("*_t7_*.tsv"), emit: t7_qc_files, optional: true
//     path("*_lambda_*.tsv"), emit: lambda_qc_files, optional: true
//     path("*_yeast_*.tsv"), emit: yeast_qc_files, optional: true



//     script:

//     out_intersect = "${basename}_intersect_${peakfile}.bed"
//     // out_0GyP = "${basename[0]}_intersect_${peakfile}.bed"

//     // out_0GyC = "${basename[2]}_intersect_${peakfile}.bed"

//     // out_plc = "${basename[1]}_intersect_${peakfile}.bed"

//     // out_cells = "${basename[3]}_intersect_${peakfile}.bed"

//     tsv_qc_file = "spikein_qc_${spike_in}_overlap_${peakfile}.tsv"

    

//     """
//     #!/usr/bin/env bash

//     # I need to only have the first 3 fields of both the peak files and the bed files
//     # using -i inplace with awk only works if you have gawk version, and this hpc does so I am fine with editing the file without changing the name.

//     # bed
//     awk -i inplace '{print \$1"\t"\$2"\t"\$3}' ${bedfiles}  

    

//     # peakfile
//     awk -i inplace '{print \$1"\t"\$2"\t"\$3}' ${peakfile}

//     # first get all of the counts for reads the a bam file in each condition that intersect with the peaks in each peak files

//     # lets do debugging

    

//     # now I just need to run bedtools on each of the 4 files in each process instance (17 total instances) but multiplied by now adding the peak files through the combine operator
    
//     # first 0GyP
//     bedtools window -a ${peakfile} -b ${bedfiles} -w 150 -bed > ${out_intersect} 

//     intersect_counts=\$(less ${out_intersect} | wc -l)
//     total_bed_counts=\$(less ${bedfiles} | wc -l)
//     percent_intersect=\$(awk "BEGIN {print (\$intersect_counts/(\$total_bed_counts+1))*100}" )


//     # now I need to get the word count of each of the conditions, that will represent how many reads are in a this instanced peak file for each condition

//     echo -e "File_base_name\tpeak_file_name\t${condition}\ttotal_${condition}_reads\tpercent_${condition}" > ${tsv_qc_file} # this is the header
//     echo -e "${grouping_name}\t${peakfile}\t\${intersect_counts}\t\${total_bed_counts}\t\${percent_intersect}" >> ${tsv_qc_file} 




//     """
// }


process overlap_window_spike_in {

    // this process might have 945 instances  17 bedfiles times the number of peakfiles

    label 'super_small_resources'

    //publishDir "${params.base_out_dir}/alignment_peak_overlap_qcspike_in/${spike_in[0]}", mode: 'copy', pattern:'*'

    //debug true

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/bedtools_rj'

    input:

    // the order of the condition and everything else is 0GyP, PLC, 0GyC, Cells

    tuple val(grouping_name), val(condition), val(spike_in), val(basename), path(bedfiles), path(peakfile)

    //path(peaks)

    // path(zero_gy)

    // path(plc)

    // path(cells)

    // path(all_peaks)



    output:
    path("${tsv_qc_file}"), emit: tsv_qc_files

    



    script:

    out_0GyP = "${basename[0]}_intersect_${peakfile}.bed"

    out_0GyC = "${basename[2]}_intersect_${peakfile}.bed"

    out_plc = "${basename[1]}_intersect_${peakfile}.bed"

    out_cells = "${basename[3]}_intersect_${peakfile}.bed"

    tsv_qc_file = "qc_${grouping_name}_overlap_${peakfile}.tsv"


    """
    #!/usr/bin/env bash

    # I need to only have the first 3 fields of both the peak files and the bed files
    # using -i inplace with awk only works if you have gawk version, and this hpc does so I am fine with editing the file without changing the name.

    # 0Gyp
    awk -i inplace '{print \$1"\t"\$2"\t"\$3}' ${bedfiles[0]}  

    # 0Gyc
    awk -i inplace '{print \$1"\t"\$2"\t"\$3}' ${bedfiles[2]}

    # plc 
    awk -i inplace '{print \$1"\t"\$2"\t"\$3}' ${bedfiles[1]}

    # cells
    awk -i inplace '{print \$1"\t"\$2"\t"\$3}' ${bedfiles[3]}

    # peakfile
    awk -i inplace '{print \$1"\t"\$2"\t"\$3}' ${peakfile}

    # first get all of the counts for reads the a bam file in each condition that intersect with the peaks in each peak files

    # lets do debugging

    echo "this is the 0GyP: ${bedfiles[0]}, the 0GyC: ${bedfiles[2]} this is plc: ${bedfiles[1]}, this is cells: ${bedfiles[3]}, and this is the peak file: ${peakfile}"

    # now I just need to run bedtools on each of the 4 files in each process instance (17 total instances) but multiplied by now adding the peak files through the combine operator
    
    # first 0GyP
    bedtools window -a ${peakfile} -b ${bedfiles[0]} -w 150 -bed > ${out_0GyP} 

    gyp_counts=\$(less ${out_0GyP} | wc -l)
    total_gyp_counts=\$(less ${bedfiles[0]} | wc -l)
    percent_gyp=\$(awk "BEGIN {print (\$gyp_counts/(\$total_gyp_counts+1))*100}" )

    # second 0GyC
    bedtools window -a ${peakfile} -b ${bedfiles[2]} -w 150 -bed > ${out_0GyC}
     
    gyc_counts=\$(less ${out_0GyC} | wc -l)
    total_gyc_counts=\$(less ${bedfiles[2]} | wc -l)
    percent_gyc=\$(awk "BEGIN {print (\$gyc_counts/(\$total_gyc_counts+1))*100}" )

    # third plc
    bedtools window -a ${peakfile} -b ${bedfiles[1]} -w 150 -bed > ${out_plc} 

    plc_counts=\$(less ${out_plc} | wc -l)
    total_plc_counts=\$(less ${bedfiles[1]} | wc -l)
    percent_plc=\$(awk "BEGIN {print (\$plc_counts/(\$total_plc_counts+1))*100}" )


    # fourth cells
    bedtools window -a ${peakfile} -b ${bedfiles[3]} -w 150 -bed > ${out_cells} 

    cells_counts=\$(less ${out_cells} | wc -l)
    total_cell_counts=\$(less ${bedfiles[3]} | wc -l)
    percent_cell=\$(awk "BEGIN {print (\$cells_counts/(\$total_cell_counts+1))*100}" )

    # now I need to get the word count of each of the conditions, that will represent how many reads are in a this instanced peak file for each condition

    echo -e "File_base_name\tpeak_file_name\t0GyP_in_peak\t0GyC_in_peak\tPLC_in_peak\tcells(200)_in_peak\ttotal_0GyP\ttotal_0GyC\ttotal_PlC\ttotal_cells\tpercent_GyP\tpercent_GyC\tpercent_PLC\tpercent_cells" > ${tsv_qc_file} # this is the header
    echo -e "${grouping_name}\t${peakfile}\t\${gyp_counts}\t\${gyc_counts}\t\${plc_counts}\t\${cells_counts}\t\${total_gyp_counts}\t\${total_gyc_counts}\t\${total_plc_counts}\t\${total_cell_counts}\t\${percent_gyp}\t\${percent_gyc}\t\${percent_plc}\t\${percent_cell}" >> ${tsv_qc_file} 




    """
}


process breakDensityWrapper_spike_in_process {

    //debug true

    conda '/ru-auth/local/home/risc_soft/miniconda3/envs/fastq2bam'

    publishDir "${params.base_out_dir}/break_density_calc/spike_in", mode: 'copy', pattern: '*'
    

    input:
    // input has to be files that the breakDensityWrapper.sh script takes
    
    // it takes the bam files that have reads aligned to the reference genome. So I think it is best to collect all of the generated bams and ensure PLC is one of the bams.
    path(bams) 
    
    path(index)

    // now putting the sorted break files
    path(sorted_breaks)
    
    // then it takes the peak files found in the directory  /lustre/fs4/home/ascortea/store/ascortea/beds
    path(peak_files)


    output:

    // output will be an AdjustedEnrichment.tsv
    path("adjustedEnrichment*.tsv"), emit: adjusted_E_tsv
    path("Adjusted_Enrichment_of_*_Plot.pdf"), emit: break_plot_pdf
    path("densityCalculations*.log"), emit: density_calc_log


    script:

    // if i get the peak file name, I can use the base name in the adjustment output file
    unique_adj_enrich_file_name = "adjustedEnrichment_spike_in_${peak_files}.tsv"
    unique_adj_plot_name = "Adjusted_Enrichment_of_${peak_files}_Break_Density_Within_Peaks_Plot.pdf"
    unique_density_calc_name = "densityCalculations_spike_in_${peak_files}.log"

    //bam_with_sort_name = "${bams.baseName}.sorted.bam"
    //index_with_sort_name = "${index.baseName}.sorted.bam.bai"


    """
    #!/usr/bin/env bash

    # nextflow can find stuff in the bin dir but not recursively, so i have to specify the sub dir

    # see if i need to rename the index files here or if doing it in the nextflow portion works
    
    # also debug here
    
    echo "this is the bam file: "${bams}". And this is the index: "${index}" "

    # one way to do this is to copy the files that breakDensityWrapper.sh calls here to this directory. instead of going in to them and going back 3 directories.    

    breakDensityWrapper.sh "${bams}" "${peak_files}"
    
    # this works but for some reason it is not seeing the output so it can be put in the published dir and also in the emit channels.

    # now making unique files that represent the names of the peak files 

    mv adjustedEnrichment.tsv "${unique_adj_enrich_file_name}" 
    mv Adjusted_Enrichment_of_Break_Density_Within_Peaks_Plot.pdf "${unique_adj_plot_name}"
    mv densityCalculations.log "${unique_density_calc_name}"


    """



}

process mk_break_points_spike_in {

    label 'normal_big_resources'
    // this is my attempt at creating the break density script, but i will try another process where I just call the break density wrapper.

    // all the wrapper does is make each peak file work with each bam file. so just make a process that takes the bam file and run all the break density scripts with each peak file. each bam will work with all the peak files in their own process instance.

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/bedtools_rj'

    publishDir "${params.base_out_dir}/break_point_bed/spike_in", mode: 'copy', pattern: '*bed'
    //publishDir './results_PE/'
    publishDir "${params.base_out_dir}/break_point_bed/bampe_unmapped_counts/spike_in", mode: 'copy', pattern: '*bampe.log'


    input:
    // take the bam files either bl filtered or not bl filtered from only the pair end path
    path(bams)
    // i think i need to have a narrowpeak file or a bed file that has peak information as input also for the next part,
    // have to call peaks and make a narrowpeak or bed peak file
    //path(bed_peaks)

    output:
    path("*breaks.bed"), emit: break_files
    path("*breaks.sorted.bed"), emit: sorted_break_bed
    path("*unmapped.bampe.log"), emit: unmapped_bampe_log


    script:

    out_bampe_name = "${bams.baseName}.bampe.bed"
    break_point_name = "${bams.baseName}.breaks.bed"
    sorted_break_point = "${bams.baseName}.breaks.sorted.bed"

    // making a new file name for the filtered breaks.bed file removing the coordinates < 0

    filt_break_point_name = "${bams.baseName}.bampe.filt.bed"

    // making a log file that records how many read ends I removed that were unmapped
    unmapped_bampe_log = "${bams.baseName}.unmapped.bampe.log"

    // making sure the bed file is sorted so just doing it again
    //sorted_bed_file = 

    """
    #!/usr/bin/env bash

    #### parameters bedtools ######
    # -bedpe : write bam alignments in bedpe format will have second field as start coordinates for forward read and 6th field as start coordinates for reverse reads
    # -i : the input bam file


    ###############################

    # this will create a bampe bed file so we can get the first field the second field and the sixth field
    # can't use bedpe since these files are not sorted by RG and we cant since the fastq files dont have that information
    # this breaks becasue when aligning to the MT genome some read ends are unmapped the chromosome and strand will be "." and the start and end coordinates will be set to -1 and sorting will not like this.
    # I can fix this by removing those lines using awk
    # if that doesnt work I can just not use bedpe option

    bedtools bamtobed \
    -bedpe \
    -i "${bams}" \
    > "${out_bampe_name}"

    

    awk '{print \$1"\t"\$2"\t"\$2}' "${out_bampe_name}" > "${break_point_name}"
    awk '{print \$1"\t"\$6"\t"\$6}' "${out_bampe_name}" >> "${break_point_name}"

    # this is saying if field 2 is greater then 0 print the entire line and append it to the new file
    awk '\$1 !="." {print \$0 }' "${break_point_name}" >> "${filt_break_point_name}"

    echo -e "unmapped read ends in "${break_point_name}":\t\$(awk '\$2 < 0 {print \$0}' "${break_point_name}" | wc -l)" >> "${unmapped_bampe_log}"

     

    # sorting the break bed files now
    bedtools sort \
    -i "${filt_break_point_name}" \
    > "${sorted_break_point}"


    numBreaks=\$(wc -l "${break_point_name}")


    #bedtools sort \
    #-i "\${bed}" \
    #> "\${sorted_break_point}"




    """
}


process break_concat_results_spike_in {

    label 'normal_small_resources'

    publishDir "${params.base_out_dir}/break_density_calc/complete_break_density_calc/spike_in", mode: 'copy', pattern: '*'

    // only using the conda for the second part
    //conda '/ru-auth/local/home/risc_soft/miniconda3/etc/profile.d/conda.sh'

    //conda '/ru-auth/local/home/risc_soft/miniconda3/envs/fastq2bam'
    
    input:
    path(adj_enrich_tsvs)

    path(density_calc_logs)



    output:

    path("${adj_enrich_tsv_out}"), emit: complete_adj_enrichment

    path("${density_calc_log_out}"), emit: complete_density_calc



    script:

    adj_enrich_tsv_out = "spike_in_full_adj_enrichment.tsv"

    density_calc_log_out = "spike_in_full_density_calc.log"

    """
    #!/usr/bin/env bash

    # I think I can separate one of the files to get the header

    head -n 1 "${adj_enrich_tsvs[0]}" > "${adj_enrich_tsv_out}"
    tail -n +2 ${adj_enrich_tsvs} >> "${adj_enrich_tsv_out}"

    # do the same for density calculations

    head -n 1 "${density_calc_logs[0]}" > "${density_calc_log_out}"
    tail -n +2 ${density_calc_logs} >> "${density_calc_log_out}"


    # I just need the full density calculations log and then I can make another full adj enrichment tsv and it will make the pdf plot
    # i can compare my full adj enrichment tsv with the one created below using parts of andrews scripts.

    source /ru-auth/local/home/risc_soft/miniconda3/etc/profile.d/conda.sh
    conda activate rstudio
    Rscript ../../../bin/calculateAdjustedEnrichmentV2.R full_density_calc.log
    echo "All done!"






    
    """
}


