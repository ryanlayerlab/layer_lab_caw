tools = params.globals.tools
skip_qc = params.globals.skip_qc
workflow wf_fastqc_fq{
    take: _input_samples
    main:
        out_ch = Channel.empty()
        if (params.step == 'mapping' && !('fastqc' in skip_qc)){
            FastQCFQ(_input_samples)
            out_ch = FastQCFQ.out
        }
    
    emit:
        fastqc_reports = out_ch
} // end of wf_fastqc_fq

process FastQCFQ {
    label 'FastQC'
    label 'cpus_2'
    label 'container_llab'

    tag {idPatient + "-" + idRun}

    publishDir "${params.outdir}/Reports/${idSample}/FastQC/${idSample}_${idRun}", 
    mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, idRun, file("${idSample}_${idRun}_R1.fastq.gz"), 
        file("${idSample}_${idRun}_R2.fastq.gz")

    output:
        file("*.{html,zip}")

    // when: !('fastqc' in skipQC) && (step == 'mapping')
    
    script:
    """
    init.sh
    fastqc -t 2 -q ${idSample}_${idRun}_R1.fastq.gz ${idSample}_${idRun}_R2.fastq.gz
    """
}
