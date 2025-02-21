
// I need access to the processes so I put them in the module file

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
    
} from '../modules/spike_in_modules.nf'



// used the t7 genome for the gloe_seq run of spike in
params.t7_genome = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/T7_Phage/genome/Sequence/WholeGenomeFasta/genome.fa')
t7_phage_genome_ch = Channel.value(params.t7_genome)

// next have to do ricc_seq runs with the yeast genome
params.yeast_genome = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/S_pombe_EF2/genome/Sequence/WholeGenomeFasta/genome.fa')
yeast_genome_ch = Channel.value(params.yeast_genome)

// now for the lambda with end_seq
params.lambda_genome = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/Lambda_cI857ind_1_Sam_7/genome/Sequence/Bowtie2Index/genome.fa')
lambda_genome_ch = Channel.value(params.lambda_genome)

workflow spike_in_runs_workflow {

    emit:
    bam_index_tuple_ch



    main:

    // gloe_seq is pair end so if I copy and paste the pe path just say params.gloe_seq

    if (params.gloe_seq) {

        // this will take the paired end reads and keep them together
        params.paired_end_reads = '/rugpfs/fs0/risc_lab/store/hcanaj/HC_GLOEseq_Novaseq_010925/fastqs_read1_read2/*_{R1,R2}*'

        pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads)

        //pe_fastqs_ch.view()

        fastp_PE_spike_in(pe_fastqs_ch.take(3))

        // checking the channels to see if everything works
        //fastp_PE_spike_in.out.filt_PE_tuple.view()

        //fastp_PE_spike_in.out.html_fastp_out.view()
        //fastp_PE_spike_in.out.failed_reads_out.view()
        //fastp_PE_spike_in.out.

        pe_filt_tuple_ch = fastp_PE_spike_in.out.filt_PE_tuple

        fastqc_PE_spike_in(pe_filt_tuple_ch)

        //fastqc_PE_spike_in.out.fastqc_zip_files.view()


        // then collect them so i can view in one file using multiqc.
        collection_fastqc_ch =fastqc_PE_spike_in.out.fastqc_zip_files.collect()

        multiqc_PE_spike_in(collection_fastqc_ch)

        //multiqc_PE.out.summary_of_PE_filt.view()

        // first use the process to index the reference genome since the process exists already for the se

        bwa_index_genome_spike_in(t7_phage_genome_ch)

        //bwa_index_genome_spike_in.out.genome_index_files.view()

        genome_index_ch = bwa_index_genome_spike_in.out.genome_index_files

        // now to use bwa aln and bwa sampe to align the filtered pair end reads to the reference genome and then to create the sam file respectively

        //pe_filt_tuple_ch.view()
        //t7_phage_genome_ch.view()
        //genome_index_ch.view()

        bwa_PE_aln_spike_in(pe_filt_tuple_ch, t7_phage_genome_ch, genome_index_ch)

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

            bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple
            flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
            norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
            tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()

            // this will give a blacklist filtered bam but i need to index it again
            bedtools_filt_blacklist(bam_index_tuple_ch, blacklist_ch)
            //bedtools_filt_blacklist.out.bl_filtered_bams.view()

            bl_filt_bams_ch = bedtools_filt_blacklist.out.bl_filtered_bams
            // so using the process to only index which means it will take the blacklist bam file
            
            samtools_bl_index(bl_filt_bams_ch) 

            bam_index_tuple_ch = samtools_bl_index.out.bl_filt_bam_index_tuple

            

        }
        else {
            samtools_sort_spike_in(sam_files_pe_ch)

            //samtools_sort_spike_in.out.bam_index_tuple.view()

            flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
            norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
            tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()

            bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple

            
        }*/

        // not running blacklist filtering on spike ins
        samtools_sort_spike_in(sam_files_pe_ch)

        //samtools_sort_spike_in.out.bam_index_tuple.view()

        flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
        norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
        tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()

        bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple

 
    }
    

    


}