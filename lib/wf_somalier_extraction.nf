tools = params.globals.tools

workflow wf_somalier_extraction{
    take: _dm_bam
    take: _fasta
    take: _fasta_fai
    take: _somalier_sites
    main:
        SomalierExtraction(_dm_bam,
                 _fasta,
                _fasta_fai,
                _somalier_sites
        )
} // end of wf_mpileup

// STEP MPILEUP.1

process SomalierExtraction {
    label 'container_llab'
    label 'cpus_8'
    tag {idSample}
    
    publishDir "${params.outdir}/Somalier_extracted/", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, file(bam), file(bai)
        file(fasta)
        file(fasta_fai)
        file(somalier_sites)

    output:
        tuple idPatient, idSample, file("${idSample}.somalier")

    when: params.somalier_sites

    script:
    """
    init.sh
    somalier extract  --sites ${somalier_sites} -f ${fasta}  ${bam}
    """
}
