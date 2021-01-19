step = params.step
tools =  params.globals.tools
status_map =  params.globals.status_map
gender_map =  params.globals.gender_map
workflow wf_recal_bam{
    take: _md_bams
    take: _bed_intervals
    take: _fasta
    take: _fasta_fai
    take: _dict
    take: _dbsnp
    take: _dbsnp_index
    take: _known_indels
    take: _known_indels_index
          
    main:
        _bams_recal = Channel.empty()
        _run_wf = params.known_indels  && step != 'variantcalling' &&
            ('haplotypecaller' in tools || 
            'mutect2' in tools ||
            'mutect2_single' in tools
            )
         if (_run_wf){
            // Create a scattered pattern 
            _int_md_bams = _md_bams.combine(_bed_intervals)
            BaseRecalibrator(
                _int_md_bams,
                _fasta,
                _fasta_fai,
                _dict,
                _dbsnp,
                _dbsnp_index,
                _known_indels,
                _known_indels_index
            )

            table_gather_bqsr_reports = 
                !params.no_intervals ? BaseRecalibrator.out.groupTuple(by:[0, 1]) : BaseRecalibrator.out
            GatherBQSRReports(table_gather_bqsr_reports)
            
            bam_apply_bqsr = _md_bams
                            .join(GatherBQSRReports.out.recal_table, by:[0,1])
            // bam_apply_bqsr = Sambamba_MD.out
            //                 .join(GatherBQSRReports.out.recal_table, by:[0,1])

            bam_apply_bqsr = bam_apply_bqsr.combine(_bed_intervals)
            ApplyBQSR(
                bam_apply_bqsr,
                _fasta,
                _fasta_fai,
                _dict
            )

            bam_merge_bam_recal = ApplyBQSR.out.groupTuple(by:[0, 1])
            /* When using intervals, merge (and in the same process index bam files)
                Which one of the MergeBamRecal, or IndexBamRecal runs is controlled by the 'when' clause
                in the processes. We finally mix them together as we know they are mutulaly exclusive.
            */
            MergeBamRecal(bam_merge_bam_recal)
            // When not using intervals, just index the bam coming from ApplyBQSR
            IndexBamRecal(bam_merge_bam_recal)
        
            _bams_recal = MergeBamRecal.out.bam_recal.mix(IndexBamRecal.out.bam_recal)
            // bams_recal_qc = MergeBamRecal.out.bam_recal_qc.mix(IndexBamRecal.out.bam_recal_qc)
        
            
            // Creating a TSV file to restart from this step
            _bams_recal.map { idPatient, idSample, bamFile, baiFile ->
                gender = gender_map[idPatient]
                status = status_map[idPatient, idSample]
                bam = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam"
                bai = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bai"
                "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"
            }.collectFile(
                name: 'recalibrated.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
            )
            // Store tsv's for individual samples as well, so we can run a single sample too
            _bams_recal
                .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") {
                    idPatient, idSample, bamFile, baiFile ->
                    status = status_map[idPatient, idSample]
                    gender = gender_map[idPatient]
                    bam = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam"
                    bai = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bai"
                    ["recalibrated_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"]
            }
        } // end of if
       
    emit:
        bams_recal = _bams_recal
}
process BaseRecalibrator {
    // label 'cpus_1'
    label 'cpus_8'
    // cache false
    // label 'memory_max'
    tag {idPatient + "-" + idSample + "-" + intervalBed.baseName}
    // tag {idPatient + "-" + idSample}

    input:
        tuple idPatient, idSample, file(bam), file(bai), file(intervalBed)
        file(fasta) 
        file(fastaFai)
        file(dict)
        file(dbsnp)
        file(dbsnpIndex) 
        file(knownIndels)
        file(knownIndelsIndex)

    output:
        tuple idPatient, idSample, file("${prefix}${idSample}.recal.table")
        // set idPatient, idSample into recalTableTSVnoInt

    script:
    dbsnpOptions = params.dbsnp ? "--known-sites ${dbsnp}" : ""
    knownOptions = params.known_indels ? knownIndels.collect{"--known-sites ${it}"}.join(' ') : ""
    prefix = params.no_intervals ? "" : "${intervalBed.baseName}_"
    intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
    // intervalsOptions = ""
    // TODO: --use-original-qualities ???
    """
    init.sh
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        BaseRecalibrator \
        -I ${bam} \
        -O ${prefix}${idSample}.recal.table \
        -R ${fasta} \
        ${intervalsOptions} \
        ${dbsnpOptions} \
        ${knownOptions} \
        --verbosity INFO
    """
}

// STEP 3.5: MERGING RECALIBRATION TABLES

process GatherBQSRReports {
    label 'memory_singleCPU_2_task'
    label 'cpus_8'
    echo true
    tag {idPatient + "-" + idSample}

    publishDir "${params.outdir}/Preprocessing/${idSample}/DuplicateMarked", mode: params.publish_dir_mode, overwrite: false

    input:
        tuple idPatient, idSample, file(recal)

    output:
        tuple idPatient, idSample, file("${idSample}.recal.table"), emit: recal_table
        // set idPatient, idSample into recalTableTSV

    script:
    input = recal.collect{"-I ${it}"}.join(' ')
    """
    init.sh
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        GatherBQSRReports \
        ${input} \
        -O ${idSample}.recal.table \
    """
}

// STEP 4: RECALIBRATING

process ApplyBQSR {
    label 'memory_singleCPU_2_task'
    label 'cpus_8'
    // label 'cpus_32'
    // label 'memory_max'
    tag {idPatient + "-" + idSample + "-" + intervalBed.baseName}
    // tag {idPatient + "-" + idSample }

    input:
        tuple idPatient, idSample, file(bam), file(bai), file(recalibrationReport), file(intervalBed)
        file(fasta)
        file(fastaFai) 
        file(dict)

    output:
        tuple idPatient, idSample, file("${prefix}${idSample}.recal.bam")

    script:
    prefix = params.no_intervals ? "" : "${intervalBed.baseName}_"
    intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
    """
    init.sh
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        ApplyBQSR \
        -R ${fasta} \
        --input ${bam} \
        --output ${prefix}${idSample}.recal.bam \
        ${intervalsOptions} \
        --bqsr-recal-file ${recalibrationReport}
    """
}

// STEP 4.5: MERGING THE RECALIBRATED BAM FILES

process MergeBamRecal {
    label 'cpus_8'

    tag {idPatient + "-" + idSample}

    publishDir "${params.outdir}/Preprocessing/${idSample}/Recalibrated", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, file(bam)

    output:
        tuple idPatient, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bai"), emit: bam_recal
        tuple idPatient, idSample, file("${idSample}.recal.bam"), emit: bam_recal_qc
        // set idPatient, idSample into bamRecalTSV

    // when: !(params.no_intervals)

    script:
    """
    init.sh
    samtools merge --threads ${task.cpus} ${idSample}.recal.bam ${bam}
    samtools index ${idSample}.recal.bam
    mv ${idSample}.recal.bam.bai ${idSample}.recal.bai
    """
}

process IndexBamRecal {
    label 'cpus_8'

    tag {idPatient + "-" + idSample}

    publishDir "${params.outdir}/Preprocessing/${idSample}/Recalibrated", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, file("${idSample}.recal.bam")

    output:
        tuple idPatient, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bam.bai"), emit: bam_recal
        tuple idPatient, idSample, file("${idSample}.recal.bam"), emit: bam_recal_qc
        
    when: params.no_intervals

    script:
    """
    init.sh
    samtools index ${idSample}.recal.bam
    """
}
