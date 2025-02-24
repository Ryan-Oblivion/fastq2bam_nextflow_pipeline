// i need to put this before the modules are loaded 
// try this here
//params.spike_name = params.gloe_seq ? '_t7' : (params.end_seq ? '_lambda' : (params.ricc_seq ? '_yeast': ''))

// I need access to the processes so I put them in the module file





// used the t7 genome for the gloe_seq run of spike in
params.t7_genome = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/T7_Phage/genome/Sequence/WholeGenomeFasta/genome.fa')
//t7_phage_genome_ch = Channel.value(params.t7_genome)
t7_genome_tuple =Channel.value(params.t7_genome).map{file -> tuple '_t7', file}
//.set{t7_genome_tuple}


// next have to do ricc_seq runs with the yeast genome
params.yeast_genome = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/S_pombe_EF2/genome/Sequence/WholeGenomeFasta/genome.fa')
//yeast_genome_ch = Channel.value(params.yeast_genome)
yeast_genome_tuple = Channel.value(params.yeast_genome).map{file -> tuple '_yeast', file}
//.set{yeast_genome_tuple}


// now for the lambda with end_seq
params.lambda_genome = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/Lambda_cI857ind_1_Sam_7/genome/Sequence/Bowtie2Index/genome.fa')
//lambda_genome_ch = Channel.value(params.lambda_genome)

lambda_genome_tuple =Channel.value(params.lambda_genome).map{file -> tuple '_lambda', file }
//.set{lambda_genome_tuple}


// maybe this has to be outside of the workflow to work
//params.spike_name = params.gloe_seq ? '_t7' : (params.end_seq ? '_lambda' : (params.ricc_seq ? '_yeast': ''))

// if i pair the genome with the spike in name I will be able to pass the proper spike in name along the process




include {
    fastp_SE_adapter_known_spike_in;
    fastp_SE_spike_in;
    fastqc_SE_spike_in;
    multiqc_SE_spike_in;
    bwa_index_genome_spike_in;
    bwa_align_SE_spike_in;
    samtools_sort_spike_in;
    deeptools_make_bed_spike_in;
    fastp_PE_spike_in;
    fastqc_PE_spike_in;
    multiqc_PE_spike_in;
    bwa_PE_aln_spike_in;
    samtools_index_sort_spike_in;
    deeptools_aln_shift_spike_in
    
} from '../modules/spike_in_modules.nf'


workflow gloe_seq_pe_workflow {

    emit:
    spike_in_bam_index_tuple_ch
    //spike_in_bed_files_norm_ch



    main:
    
    
    // testing something where I make a variable based on the spike in being ran
    //spike_name = '_t7' // this will be placed dynamically in the process as part of the output file names and dir

    // this will take the paired end reads and keep them together
    params.paired_end_reads = '/rugpfs/fs0/risc_lab/store/hcanaj/HC_GLOEseq_Novaseq_010925/fastqs_read1_read2/*_{R1,R2}*'

    pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads)

    //pe_fastqs_ch.view()

    fastp_PE_spike_in(pe_fastqs_ch.take(3),t7_genome_tuple)

    // checking the channels to see if everything works
    //fastp_PE_spike_in.out.filt_PE_tuple.view()

    //fastp_PE_spike_in.out.html_fastp_out.view()
    //fastp_PE_spike_in.out.failed_reads_out.view()
    //fastp_PE_spike_in.out.

    pe_filt_tuple_ch = fastp_PE_spike_in.out.filt_PE_tuple

    fastqc_PE_spike_in(pe_filt_tuple_ch,t7_genome_tuple)

    //fastqc_PE_spike_in.out.fastqc_zip_files.view()


    // then collect them so i can view in one file using multiqc.
    collection_fastqc_ch =fastqc_PE_spike_in.out.fastqc_zip_files.collect()

    multiqc_PE_spike_in(collection_fastqc_ch,t7_genome_tuple)

    //multiqc_PE.out.summary_of_PE_filt.view()

    // first use the process to index the reference genome since the process exists already for the se

    bwa_index_genome_spike_in(t7_genome_tuple)

    //bwa_index_genome_spike_in.out.genome_index_files.view()

    genome_index_ch = bwa_index_genome_spike_in.out.genome_index_files

    // now to use bwa aln and bwa sampe to align the filtered pair end reads to the reference genome and then to create the sam file respectively

    //pe_filt_tuple_ch.view()
    //t7_phage_genome_ch.view()
    //genome_index_ch.view()

    bwa_PE_aln_spike_in(pe_filt_tuple_ch, genome_index_ch, t7_genome_tuple)

    // now check to see if the output channels are good
    //bwa_PE_aln_spike_in.out.pe_sam_files.view()
    
    
    // now i need to make the parameters for  if the bam file will be blacklist filtered or not
    
    // using this channel for both if Blacklist or not
    sam_files_pe_ch = bwa_PE_aln_spike_in.out.pe_sam_files

    /*if (params.BL) {

        
        // i have to make a bam file to then use bedtools intersect to get the blacklist
        // using the same samtools sort process found in the SE part of the pipeline
        samtools_sort_spike_in(sam_files_pe_ch)

        //samtools_sort_spike_in.out.bam_index_tuple.view()

        spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple
        flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
        norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
        tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()

        // this will give a blacklist filtered bam but i need to index it again
        bedtools_filt_blacklist(spike_in_bam_index_tuple_ch, blacklist_ch)
        //bedtools_filt_blacklist.out.bl_filtered_bams.view()

        bl_filt_bams_ch = bedtools_filt_blacklist.out.bl_filtered_bams
        // so using the process to only index which means it will take the blacklist bam file
        
        samtools_bl_index(bl_filt_bams_ch) 

        spike_in_bam_index_tuple_ch = samtools_bl_index.out.bl_filt_bam_index_tuple

        

    }
    else {
        samtools_sort_spike_in(sam_files_pe_ch)

        //samtools_sort_spike_in.out.bam_index_tuple.view()

        flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
        norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
        tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()

        spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple

        
    }*/

    // not running blacklist filtering on spike ins
    samtools_sort_spike_in(sam_files_pe_ch, t7_genome_tuple)

    //samtools_sort_spike_in.out.bam_index_tuple.view()

    flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
    norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
    tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()

    spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple

        




}


workflow end_seq_pe_workflow {

    emit: 
    spike_in_bam_index_tuple_ch


    main:

    // testing something where I make a variable based on the spike in being ran
    //spike_name = '_t7' // this will be placed dynamically in the process as part of the output file names and dir

    // this will take the paired end reads and keep them together
    params.paired_end_reads = '/rugpfs/fs0/risc_lab/store/hcanaj/HC_GLOEseq_Novaseq_010925/fastqs_read1_read2/*_{R1,R2}*'

    pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads)

    //pe_fastqs_ch.view()

    fastp_PE_spike_in(pe_fastqs_ch.take(3),lambda_genome_tuple)

    // checking the channels to see if everything works
    //fastp_PE_spike_in.out.filt_PE_tuple.view()

    //fastp_PE_spike_in.out.html_fastp_out.view()
    //fastp_PE_spike_in.out.failed_reads_out.view()
    //fastp_PE_spike_in.out.

    pe_filt_tuple_ch = fastp_PE_spike_in.out.filt_PE_tuple

    fastqc_PE_spike_in(pe_filt_tuple_ch,lambda_genome_tuple)

    //fastqc_PE_spike_in.out.fastqc_zip_files.view()


    // then collect them so i can view in one file using multiqc.
    collection_fastqc_ch =fastqc_PE_spike_in.out.fastqc_zip_files.collect()

    multiqc_PE_spike_in(collection_fastqc_ch,lambda_genome_tuple)

    //multiqc_PE.out.summary_of_PE_filt.view()

    // first use the process to index the reference genome since the process exists already for the se

    bwa_index_genome_spike_in(lambda_genome_tuple)

    //bwa_index_genome_spike_in.out.genome_index_files.view()

    genome_index_ch = bwa_index_genome_spike_in.out.genome_index_files

    // now to use bwa aln and bwa sampe to align the filtered pair end reads to the reference genome and then to create the sam file respectively

    //pe_filt_tuple_ch.view()
    //t7_phage_genome_ch.view()
    //genome_index_ch.view()

    bwa_PE_aln_spike_in(pe_filt_tuple_ch,  genome_index_ch, lambda_genome_tuple)

    // now check to see if the output channels are good
    //bwa_PE_aln_spike_in.out.pe_sam_files.view()
    
    
    // now i need to make the parameters for  if the bam file will be blacklist filtered or not
    
    // using this channel for both if Blacklist or not
    sam_files_pe_ch = bwa_PE_aln_spike_in.out.pe_sam_files

    /*if (params.BL) {

        
        // i have to make a bam file to then use bedtools intersect to get the blacklist
        // using the same samtools sort process found in the SE part of the pipeline
        samtools_sort_spike_in(sam_files_pe_ch)

        //samtools_sort_spike_in.out.bam_index_tuple.view()

        spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple
        flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
        norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
        tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()

        // this will give a blacklist filtered bam but i need to index it again
        bedtools_filt_blacklist(spike_in_bam_index_tuple_ch, blacklist_ch)
        //bedtools_filt_blacklist.out.bl_filtered_bams.view()

        bl_filt_bams_ch = bedtools_filt_blacklist.out.bl_filtered_bams
        // so using the process to only index which means it will take the blacklist bam file
        
        samtools_bl_index(bl_filt_bams_ch) 

        spike_in_bam_index_tuple_ch = samtools_bl_index.out.bl_filt_bam_index_tuple

        

    }
    else {
        samtools_sort_spike_in(sam_files_pe_ch)

        //samtools_sort_spike_in.out.bam_index_tuple.view()

        flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
        norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
        tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()

        spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple

        
    }*/

    // not running blacklist filtering on spike ins
    samtools_sort_spike_in(sam_files_pe_ch, lambda_genome_tuple)

    //samtools_sort_spike_in.out.bam_index_tuple.view()

    flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
    norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
    tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()

    spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple

    



}