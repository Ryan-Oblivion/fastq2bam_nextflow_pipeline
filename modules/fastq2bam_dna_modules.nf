process samtools_index_sort {
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


process mk_break_points {
    // this is my attempt at creating the break density script, but i will try another process where I just call the break density wrapper.

    // all the wrapper does is make each peak file work with each bam file. so just make a process that takes the bam file and run all the break density scripts with each peak file. each bam will work with all the peak files in their own process instance.

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/bedtools_rj'

    publishDir './results_PE/break_point_bed', mode: 'copy', pattern: '*_breaks.bed'
    publishDir './results_PE/'

    input:
    // take the bam files either bl filtered or not bl filtered from only the pair end path
    path(bams)
    // i think i need to have a narrowpeak file or a bed file that has peak information as input also for the next part,
    // have to call peaks and make a narrowpeak or bed peak file
    //path(bed_peaks)

    output:


    script:

    out_bampe_name = "${bams.baseName}_bampe.bed"
    break_point_name = "${bams.baseName}_breaks.bed"
    sorted_break_point = "${bams.baseName}_sorted.bed"

    // making sure the bed file is sorted so just doing it again
    sorted_bed_file = 

    """
    #!/usr/bin/env bash

    #### parameters bedtools ######
    # -bedpe : write bam alignments in bedpe format will have second field as start coordinates for forward read and 6th field as start coordinates for reverse reads
    # -i : the input bam file


    ###############################

    # this will create a bampe bed file so we can get the first field the second field and the sixth field

    bedtools bamtobed \
    -bedpe \
    -i "${bams}" \
    > "${out_bampe_name}"

    awk '{print \$1"\t"\$2"\t"\$2}' "${out_bampe_name}" > "${break_point_name}"
    awk '{print \$1"\t"\$6"\t"\$6}' "${out_bampe_name}" >> "${break_point_name}"

    # sorting the break bed files now
    bedtools sort \
    -i "${break_point_name}" \
    > "${sorted_break_point}"


    numBreaks=\$(wc -l "${break_point_name}")


    bedtools sort \
    -i "${bed}" \
    > 




    """
}




process breakDensityWrapper_process {

    publishDir './break_density_calc', mode: 'copy', pattern: '*'
    

    input:
    // input has to be files that the breakDensityWrapper.sh script takes
    
    // it takes the bam files that have reads aligned to the reference genome. So I think it is best to collect all of the generated bams and ensure PLC is one of the bams.
    path(bams)
    
    // then it takes the peak files found in the directory  /lustre/fs4/home/ascortea/store/ascortea/beds
    path(peak_files)


    output:

    // output will be an AdjustedEnrichment.tsv
    path("adjustedEnrichment.tsv"), emit: adjusted_E_tsv
    path("Adjusted_Enrichment_of_*_Plot.pdf"), emit: break_plot_pdf
    path("densityCalculations.log"), emit: density_calc_log


    script:


    """
    #!/usr/bin/env bash

    # nextflow can find stuff in the bin dir but not recursively, so i have to specify the sub dir

    breakDensityWrapper.sh "${bams}" "${peak_files}"
    
    # this works but for some reason it is not seeing the output so it can be put in the published dir and also in the emit channels.




    """



}
