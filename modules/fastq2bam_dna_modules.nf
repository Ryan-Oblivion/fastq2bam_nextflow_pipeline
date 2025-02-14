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



/*
process breakDensityWrapper {


    input:


    output:



    script:


    """
    #!/usr/bin/env bash


    source /ru-auth/local/home/risc_soft/miniconda3/etc/profile.d/conda.sh

    #Main
    # Iterate through arguments looking for bams, and for each bam iterate through arguments again looking for beds or peaks to construct command for getting break density within the bed.
    #breakDensities=()
    echo "Making densityCalculations.log\n"
    echo "bam	bed	Expected_Density	Observed_Density	Enrichment" > densityCalculations.log
    for peak in $@
        do 
            echo "iterating"
            if [[ \$peak == *.bed ]] || [[ \$peak == *Peak ]]
            then
                echo "Working on \$peak"
                for bam in $@
                    do
                        if [[ \$bam == *.bam ]]
                            then
                                echo "Working on \$bam"
                                echo "expectedDensity"
                                conda activate rstudio
                                expectedDensity=$(Rscript /lustre/fs4/home/ascortea/Risc_scratch/ascortea/scripts/BreakDensity/getUniformBreakDensity.R \$peak hg19)
                                echo \$expectedDensity
                                conda activate fastq2bam
                                bash /lustre/fs4/home/ascortea/Risc_scratch/ascortea/scripts/BreakDensity/getBreakDensityInPeaksV3.sh \$bam \$peak
                                observedDensity=$(awk '{ total += \$4 } END { print total }' \${bam##}.\${peak##}.numBreaksInPeaksNormalized.bed)
                                echo \$observedDensity
                                enrichment=$(awk -v obs="\$observedDensity" -v expt="\$expectedDensity" 'BEGIN {print obs/expt}')
                                echo \$enrichment
                                echo "saving results into densityCalculations.log"
                                echo "\$bam	\$peak	\$expectedDensity	\$observedDensity	\$enrichment" >> densityCalculations.log
                        fi
                    done
            fi
        done
    echo "Making plot and .tsv of adjustedEnrichments"
    conda activate rstudio
    Rscript /lustre/fs4/home/ascortea/Risc_scratch/ascortea/scripts/BreakDensity/calculateAdjustedEnrichmentV2.R densityCalculations.log
    echo "All done!"





    """



}*/
