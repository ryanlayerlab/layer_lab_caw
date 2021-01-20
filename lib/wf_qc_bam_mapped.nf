skip_qc = params.globals.skip_qc

workflow wf_qc_bam_mapped{
    take: _bams_mapped // tuple idPatient, idSample, file(bam), file(bai)
    take: _target_bed
    main:
        // For the qc, we do not need the bai's
        _bams_mapped_qc = _bams_mapped.map
                                { idPatient, idSample, bam, bai -> 
                                [idPatient, idSample, bam]
                                }
        _bam_qc_out = Channel.empty()
        if (!('bamqc' in skip_qc)){
            BamQC(
            _bams_mapped_qc,
            _target_bed
            )
            _bam_qc_out = BamQC.out
        }
    emit: bam_qc = _bam_qc_out
}

process BamQC {
    // label 'memory_max'
    label 'container_llab'
    label 'cpus_16'
    // cache false

    tag {idPatient + "-" + idSample}

    // publishDir "${params.outdir}/Reports/${idSample}/bamQC", mode: params.publish_dir_mode
    publishDir "${params.outdir}/Reports/${idSample}/bamQC/", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, file(bam) 
        file(targetBED)

    output:
        file("${bam.baseName}")

    // when: !('bamqc' in skipQC)

    script:
    use_bed = params.target_bed ? "-gff ${targetBED}" : ''
    """
    init.sh
    qualimap --java-mem-size=${task.memory.toGiga()}G \
        bamqc \
        -bam ${bam} \
        --paint-chromosome-limits \
        --genome-gc-distr HUMAN \
        $use_bed \
        -nt ${task.cpus} \
        -skip-duplicated \
        --skip-dup-mode 0 \
        -outdir ${bam.baseName} \
        -outformat HTML
    """
}