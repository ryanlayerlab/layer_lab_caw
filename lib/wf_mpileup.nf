tools = params.globals.tools

workflow wf_mpileup{
    take: _int_bam_recal
     take: _fasta
    take: _fasta_fai
    main:
        Mpileup(_int_bam_recal,
                 _fasta,
                _fasta_fai 
        )
        
        MergeMpileup(
            Mpileup.out.groupTuple(by:[0, 1]) 
        )
} // end of wf_mpileup

// STEP MPILEUP.1

process Mpileup {
    label 'memory_singleCPU_2_task'
    tag {idSample + "-" + intervalBed.baseName}
    
    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: { it == "${idSample}.pileup.gz" ? "VariantCalling/${idSample}/mpileup/${it}" : '' }

    input:
        tuple idPatient, idSample, file(bam), file(bai), file(intervalBed)
        file(fasta)
        file(fastaFai)

    output:
        tuple idPatient, idSample, file("${prefix}${idSample}.pileup.gz")

    when: 'controlfreec' in tools || 'mpileup' in tools

    script:
    prefix = params.no_intervals ? "" : "${intervalBed.baseName}_"
    intervalsOptions = params.no_intervals ? "" : "-l ${intervalBed}"
    """
    init.sh
    samtools mpileup \
        -f ${fasta} ${bam} \
        ${intervalsOptions} \
    | bgzip --threads ${task.cpus} -c > ${prefix}${idSample}.pileup.gz
    """
}


// STEP MPILEUP.2 - MERGE

process MergeMpileup {
    tag {idSample}

    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: { it == "${idSample}.pileup.gz" ? "VariantCalling/${idSample}/mpileup/${it}" : '' }

    input:
        tuple idPatient, idSample, file(mpileup)

    output:
        tuple idPatient, idSample, file("${idSample}.pileup.gz")

    when: !(params.no_intervals) && 'controlfreec' in tools || 'mpileup' in tools

    script:
    """
    init.sh
    for i in `ls -1v *.pileup.gz`;
        do zcat \$i >> ${idSample}.pileup
    done

    bgzip --threads ${task.cpus} -c ${idSample}.pileup > ${idSample}.pileup.gz

    rm ${idSample}.pileup
    """
}