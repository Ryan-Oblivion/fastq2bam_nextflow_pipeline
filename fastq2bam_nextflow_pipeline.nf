// first i will activate dsl2
nextflow.enable.dsl=2


// params.pe_out_dir = './results_PE'
// params.se_out_dir = './results_SE'

params.base_out_dir = params.PE ? './results_PE' : (params.SE ? './results_SE' : '')
// finally using modules



include {
    fastp_SE_adapter_known;
    fastp_SE;
    fastqc_SE;
    multiqc_SE;
    bwa_index_genome;
    bwa_align_SE;
    samtools_sort;
    deeptools_make_bed;
    bedtools_filt_blacklist;
    samtools_bl_index;
    fastp_PE;
    fastqc_PE;
    multiqc_PE;
    bwa_PE_aln;
    multiqc_bam_stats;
    deeptools_aln_shift;
    samtools_index_sort;
    mk_break_points;
    breakDensityWrapper_process;
    py_calc_stats_log;





} from './modules/fastq2bam_dna_modules.nf'



//include {samtools_index_sort} from './modules/fastq2bam_dna_modules.nf'



// i want to include another workflow script

include {breakDensityWrapper_workflow} from './workflows/breakDensityWrapper_workflow.nf'

//include {py_calc_stats_log} from './modules/fastq2bam_dna_modules.nf'

//include {spike_in_runs_workflow} from './workflows/spike_in_workflow.nf'
include {
    pe_t7_spike_in_workflow;
    pe_lambda_spike_in_workflow;
    pe_yeast_spike_in_workflow;
    se_t7_spike_in_workflow;
    se_lambda_spike_in_workflow;
    se_yeast_spike_in_workflow  
            
            
} from './workflows/2nd_spike_in_workflow.nf'



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

    

    // params.single_end_reads = file('/rugpfs/fs0/risc_lab/store/hcanaj/HC_ENDseq_Novaseq_010925/read1_fastqs/*_1.fastq.gz')
    // params.paired_end_reads = '/rugpfs/fs0/risc_lab/store/hcanaj/HC_GLOEseq_Novaseq_010925/fastqs_read1_read2/*_{R1,R2}*'

        

    // // this will give the user to run in test mode where the pipeline will only take 3 of the fastq files in the directory full of fastq files
    // // or it will run in normal mode where you want all your data processed
    // if (params.test) {
    //     se_reads_files = Channel.fromPath(params.single_end_reads).take(3)
    //     pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads).take(3)
    
    // }else {

    //     se_reads_files = Channel.fromPath(params.single_end_reads)
    //     pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads)

    // }
    
    
    if ( params.SE ) {

        
            

            

            // lets get the channel for the single end reads first
            // only use the single end read 1 data from the end seq which are already stored here: /rugpfs/fs0/risc_lab/store/hcanaj/HC_ENDseq_Novaseq_010925/read1_fastqs

            params.single_end_reads = file('/rugpfs/fs0/risc_lab/store/hcanaj/HC_ENDseq_Novaseq_010925/read1_fastqs/*_1.fastq.gz')
    

        

            // this will give the user to run in test mode where the pipeline will only take 3 of the fastq files in the directory full of fastq files
            // or it will run in normal mode where you want all your data processed
            if (params.test) {
                se_reads_files = Channel.fromPath(params.single_end_reads).take(3)
                
            
            }else {

                se_reads_files = Channel.fromPath(params.single_end_reads)
               

            }



            //params.single_end_reads = file('/rugpfs/fs0/risc_lab/store/hcanaj/HC_ENDseq_Novaseq_010925/read1_fastqs/*_1.fastq.gz')
            //se_reads_files = Channel.fromPath(params.single_end_reads)
            
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

                fastp_SE_adapter_known(se_reads_files, se_reads_name, adapter_ch) // will have to make a new process for if the adapter sequence is known

                fastq_filts = fastp_SE_adapter_known.out.filtered_fastqs
                //fastp_SE.out.view()
                fastq_filts.map{file -> file.baseName}
                            .set{fastq_filts_name}

                // now getting the html files since i think fastqc combines them into one, that might be multiqc
                fastp_filt_html = fastp_SE_adapter_known.out.fastp_html_reports

            }    
            else {

                fastp_SE(se_reads_files, se_reads_name)

                    // take all of the filtered fastq files and put them in a channel name
                // since the fastq files might be in a different order, if i need to get their base names I will have to do it from this new channel below
                fastq_filts = fastp_SE.out.filtered_fastqs
                //fastp_SE.out.view()
                fastq_filts.map{file -> file.baseName}
                            .set{fastq_filts_name}

                // now getting the html files since i think fastqc combines them into one, that might be multiqc
                fastp_filt_html = fastp_SE.out.fastp_html_reports



            }
            

            //fastp_SE(se_reads_files, se_reads_name) // REMEMBER TO REMOVE THIS TESTING FEATURE WHERE IT WILL ONLY TAKE THE FIRST 3

            

            
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
            bam_index_tuple_ch = samtools_sort.out.bam_index_tuple
            flagstat_log_ch = samtools_sort.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
            norm_stats_txt_ch = samtools_sort.out.norm_stats_txt.collect()
            tsv_SN_stats_ch = samtools_sort.out.tsv_SN_stats.collect()
            // if you want to filter black list use the param --BL in the command line when calling nextflow
            if ( params.BL ) {

                // using bedtools to filter black list but first giving the user an option to put the correct black list for an organism
                // by defualt it will use the hg19 v2 blacklist, but if you used a different organism or human genome use the appropriate blacklist
                //params.blacklist_path = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg19/blacklist/hg19-blacklist.v2.bed')
                
                //blacklist_ch = Channel.value(params.blacklist_path)

                bedtools_filt_blacklist(bam_index_tuple_ch, blacklist_ch)

                bl_filt_bams_ch = bedtools_filt_blacklist.out.bl_filtered_bams
                // i will need to index the black list filtered bam again so i have to create a different samtools process for this
                samtools_bl_index(bl_filt_bams_ch)

                bam_index_tuple_ch = samtools_bl_index.out.bl_filt_bam_index_tuple

                /*if ( params.ATAC ) {


                    // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

                    deeptools_aln_shift(bam_index_tuple_ch)

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

                    deeptools_make_bed(bam_index_tuple_ch)

                    deeptools_make_bed.out.bed_files_normalized.view()

                    bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized          


                }*/

                // then i need to pass the indexed_bl_bam and the bam to the deeptools process

                //deeptools_make_bed(bam_index_tuple_ch)
                //bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized

            }

            /*
            else {


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

                    deeptools_make_bed.out.bed_files_normalized.view()

                    bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized



                }
                else {



                    // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

                    deeptools_make_bed(bam_index_tuple_ch)

                    deeptools_make_bed.out.bed_files_normalized.view()

                    bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized          


                }

                // Now i want to pass the tuple that has the bam and its corresponding index file into a process that will create a bigwig file for visulization, created from the bam file. This will show read coverage in the genome without looking for significant areas
                //deeptools_make_bed(bam_index_tuple_ch)

                //deeptools_make_bed.out.bed_files_normalized.view()
                //bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized
            }*/
            
    

            // i will either use the bams or the bed files for any future processes depending on what tool needs what.
    }

    // now I need to make the paired-end part of this pipeline.
    else if (params.PE) {

        // this will take the paired end reads and keep them together
        // params.paired_end_reads = '/rugpfs/fs0/risc_lab/store/hcanaj/HC_GLOEseq_Novaseq_010925/fastqs_read1_read2/*_{R1,R2}*'

        // pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads)

        
        params.paired_end_reads = '/rugpfs/fs0/risc_lab/store/hcanaj/HC_GLOEseq_Novaseq_010925/fastqs_read1_read2/*_{R1,R2}*'

            

        // this will give the user to run in test mode where the pipeline will only take 3 of the fastq files in the directory full of fastq files
        // or it will run in normal mode where you want all your data processed
        if (params.test) {
            
            pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads).take(3)
        
        }else {

            
            pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads)

        }


        //pe_fastqs_ch.view()

        fastp_PE(pe_fastqs_ch)

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

            bam_index_tuple_ch = samtools_bl_index.out.bl_filt_bam_index_tuple

            // so this would give a bam that is bl filtered and has an index

            /*if ( params.ATAC ) {


                // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

                deeptools_aln_shift(bam_index_tuple_ch)

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

                deeptools_make_bed(bam_index_tuple_ch)

                //deeptools_make_bed.out.bed_files_normalized.view()

                bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized          


            }*/
            
            
            // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

            //deeptools_make_bed(bam_index_tuple_ch)

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

            /*if (params.ATAC) {


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


            }*/






            // just using the original sorted and indexed bam
            
            // i will comment this out below since it was the original but the above if logic should work if there is ATAC data or if there is not

            //deeptools_make_bed(bam_index_tuple_ch)

            //deeptools_make_bed.out.bed_files_normalized.view()

            //bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized


        }

        

 
    }



    // I think i can just put the if ATAC script here separately and just grab all of the bam tuple channels

    if (params.ATAC) {


        // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

        deeptools_aln_shift(bam_index_tuple_ch)

        atac_shift_bam_ch = deeptools_aln_shift.out.atac_shifted_bam

        // now i have to re index this new atac shifted bam. dispite the name of the process I can just pass any future created bam to this channel to be indexed

        samtools_index_sort(atac_shift_bam_ch)

        // so now name the tuple channel output appropriately 
        atac_shift_bam_index_ch = samtools_index_sort.out.bam_index_tuple // i changed the emit ch to just be bam_index_tuple

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


    // making a multiqc process for the samtools flagstat log files. this should be able to take the flagstat_log_ch from any part of the choosen paths
    multiqc_bam_stats(flagstat_log_ch, norm_stats_txt_ch)


    
    // I want to have a parameter that takes peakfiles. The default will be the IMR90 narrowPeak files
    params.peak_files = files('/lustre/fs4/home/ascortea/store/ascortea/beds/IMR90/*.narrowPeak')
    //now making the channel for the files
    peak_files_ch = Channel.value(params.peak_files)


    if (params.calc_break_density){
    // i want to call the workflow breakDensityWrapper_workflow and pass the bam_index_tuple_ch as an input from either path where the user chose to do blacklist filter or not. Then also pass the peak files that already exists or are created later as input
    
        if (params.PE) {
            
            breakDensityWrapper_workflow(bam_index_tuple_ch, peak_files_ch)

        }

        if (params.SE) {
         
            breakDensityWrapper_workflow(bam_index_tuple_ch, peak_files_ch)
  
        }
    }

    // I want to make a log file with all the stats from using samtools stats on each bam file

    if (params.PE) {

        //py_calc_stats_log(tsv_SN_stats_ch)
    }
    if (params.SE) {

        //py_calc_stats_log(tsv_SN_stats_ch)

    }

    
    // looking to run the spike_in workflow
    

    if (params.spike_in) {

        // all this is doing is running the normal fastq2bam2 pipeline but with the specified genomes for spike in
       
        if (params.PE) {

            if (params.t7 ) {

                pe_t7_spike_in_workflow()
                // now getting the output channels from the workflows just incase i need to use them in a downstream analysis
                pe_t7_bam_index_tuple_ch = pe_t7_spike_in_workflow.out.spike_in_bam_index_tuple_ch
                pe_t7_bed_files = pe_t7_spike_in_workflow.out.spike_in_bed_files_norm_ch
            }
            if (params.lambda) {

                pe_lambda_spike_in_workflow()

                pe_lambda_bam_index_tuple_ch = pe_lambda_spike_in_workflow.out.spike_in_bam_index_tuple_ch
                pe_lambda_bed_files = pe_lambda_spike_in_workflow.out.spike_in_bed_files_norm_ch
            }
            if (params.yeast) {

                // so i would put the yeast spike in workflow here, for example.
                pe_yeast_spike_in_workflow()

                pe_yeast_bam_index_tuple_ch = pe_yeast_spike_in_workflow.out.spike_in_bam_index_tuple_ch
                pe_yeast_bed_files = pe_yeast_spike_in_workflow.out.spike_in_bed_files_norm_ch
            }
        }
        if (params.SE) {

            //for this i need to make the single end workflow for spike ins
            if (params.t7) {

                se_t7_spike_in_workflow()

                se_t7_bam_index_tuple_ch = se_t7_spike_in_workflow.out.spike_in_bam_index_tuple_ch
                se_t7_bed_files = se_t7_spike_in_workflow.out.spike_in_bed_files_norm_ch
            }
            if (params.lambda) {

                se_lambda_spike_in_workflow()
                se_lambda_bam_index_tuple_ch = se_lambda_spike_in_workflow.out.spike_in_bam_index_tuple_ch
                se_lambda_bed_files = se_lambda_spike_in_workflow.out.spike_in_bed_files_norm_ch
            }           
            if (params.yeast) {

                // so i would put the yeast spike in workflow here, for example.
                se_yeast_spike_in_workflow()
                se_yeast_bam_index_tuple_ch = se_yeast_spike_in_workflow.out.spike_in_bam_index_tuple_ch
                se_yeast_bed_files = se_yeast_spike_in_workflow.out.spike_in_bed_files_norm_ch
            }

        }

        // not sure if i would need to make this channel the normal 'bam_index_tuple_ch', or do what i did with atac-seq and give it a unique name calling it 'spike_in_bam_index_ch'
        // since spike in can occur at the same time as either pe or se, i need to give this channel a unique name. I dont think spike in bam files need to have break density calculated
        //spike_in_bam_index_ch = spike_in_runs_workflow.out.spike_in_bam_index_tuple_ch
        //spike_in_bed_files_norm_ch = spike_in_runs_workflow.out.spike_in_bed_files_norm_ch

    }

    
    
}