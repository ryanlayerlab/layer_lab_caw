_skip_qc = params.globals.skip_qc
workflow wf_qc_bam_recal{

    // take: _bam_raw_qc // tuple idPatient, idSample, file(bam), file(bam // take the raw bams
    take: _bams_recal // tuple idPatient, idSample, file(bam), file(bai)
    // take: _bams_recal_on_target // tuple idPatient, idSample, file(bam), file(bai)
    take: _target_bed
    take: _bait_bed
    take: _fasta
    take: _fasta_fai
    take: _dict

    main:
     // For the qc, we do not need the bai's
        _bams_recal_qc = _bams_recal.map
                                { idPatient, idSample, bam, bai -> 
                                [idPatient, idSample, bam]
                                }
        
        _ch_samtools_stats = Channel.empty()
        if(!('samtools' in _skip_qc)){
            SamtoolsStats(_bams_recal_qc)
            _ch_samtools_stats = SamtoolsStats.out
        }

        _ch_alignment_summary = Channel.empty()
        if(! ('alignment_summary' in _skip_qc)) {
            CollectAlignmentSummaryMetrics(
                _bams_recal_qc,
                _fasta,
                _fasta_fai,
                _dict
            )
            _ch_alignment_summary = CollectAlignmentSummaryMetrics.out
        }
        _ch_insert_size_metrics = Channel.empty()
        _ch_insert_size_pdfs = Channel.empty()
        if(!('insert_size_metrics' in _skip_qc)) {
            CollectInsertSizeMetrics(
                _bams_recal_qc
            )
            _ch_insert_size_metrics = CollectInsertSizeMetrics.out[0]
            _ch_insert_size_pdfs = CollectInsertSizeMetrics.out[1]
        }

       _ch_hs_metrics = Channel.empty()
       if(!('hs_metrics' in _skip_qc) && params.bait_bed){
            CollectHsMetrics(
                _bams_recal_qc,
                _fasta,
                _fasta_fai,
                _dict,
                _target_bed,
                _bait_bed
                // '' // this is the suffix to the output file
            )
            _ch_hs_metrics = CollectHsMetrics.out
       }
       
    emit:
        samtools_stats =  _ch_samtools_stats
        alignment_summary_metrics = _ch_alignment_summary
        insert_size_metrics = _ch_insert_size_metrics
        insert_size_metrics_pdf = _ch_insert_size_pdfs
        hs_metrics = _ch_hs_metrics
        // bam_qc = BamQC.out
} // end of wf_qc_bams_recal


// STEP 5: QC

process SamtoolsStats {
    label 'container_llab'
    label 'cpus_2'

    tag {idPatient + "-" + idSample}

    publishDir "${params.outdir}/Reports/${idSample}/SamToolsStats", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, file(bam)

    output:
        file ("${bam}.samtools.stats.out")

    // when: !('samtools' in _skip_qc)

    script:
    """
    init.sh
    samtools stats ${bam} > ${bam}.samtools.stats.out
    """
}



process CollectAlignmentSummaryMetrics{
    label 'container_llab'
    label 'cpus_16'
    tag {idPatient + "-" + idSample}
    
    publishDir "${params.outdir}/Reports/${idSample}/alignment_summary/", mode: params.publish_dir_mode
    
    input:
    tuple idPatient, idSample, file(bam) 
    file(fasta) 
    file(fastaFai)
    file(dict)

    output:
    file("${bam.baseName}_alignment_metrics.txt")
    
    // when: ! ('alignment_summary' in _skip_qc)
    
    script:
    """
    init.sh
    gatk --java-options -Xmx32G CollectAlignmentSummaryMetrics --VALIDATION_STRINGENCY LENIENT \
    -I ${bam} \
    -O ${bam.baseName}_alignment_metrics.txt \
    -R ${fasta}
    """
}

process CollectInsertSizeMetrics{
    label 'container_llab'
    label 'cpus_16'
    tag {idPatient + "-" + idSample}
    
    publishDir "${params.outdir}/Reports/${idSample}/insert_size_metrics/", mode: params.publish_dir_mode
    
    input:
    tuple idPatient, idSample, file(bam)

    output:
    file("${bam.baseName}_insert_size_metrics.txt")
    file("${bam.baseName}_insert_size_histogram.pdf")
    
    
    // when: !('insert_size_metrics' in _skip_qc)

    script:
    """
    init.sh
    gatk --java-options -Xmx32G CollectInsertSizeMetrics --VALIDATION_STRINGENCY LENIENT \
    -I ${bam} \
    -O ${bam.baseName}_insert_size_metrics.txt \
    -H ${bam.baseName}_insert_size_histogram.pdf 
    """
}

process CollectHsMetrics{
    label 'container_llab'
    label 'cpus_16'
    tag {idPatient + "-" + idSample}
    
    publishDir "${params.outdir}/Reports/${idSample}/hs_metrics/", mode: params.publish_dir_mode
    
    input:
    // tuple idPatient, idSample, file(bam), file(bai)
    tuple idPatient, idSample, file(bam)
    file(fasta) 
    file(fastaFai)
    file(dict)
    file(targetBED)
    file(baitBED)
    // val (output_suffix)

    output:
    file("${bam.baseName}.txt")
    // file("${bam.baseName}_${output_suffix}.txt")
    
    
    // when: !('hs_metrics' in _skip_qc) && params.bait_bed
    script:
    """
    init.sh
    gatk BedToIntervalList -I ${targetBED} -O target.interval_list -SD ${dict}
    gatk BedToIntervalList -I ${baitBED} -O bait.interval_list -SD ${dict}

    gatk --java-options -Xmx32G CollectHsMetrics --VALIDATION_STRINGENCY LENIENT \
    -I ${bam} \
    -O ${bam.baseName}.txt \
    -TI target.interval_list \
    -BI bait.interval_list \
    -R ${fasta}
    """
}
