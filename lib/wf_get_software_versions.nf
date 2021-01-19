
workflow wf_get_software_versions{
    main:
        _ch_sv = Channel.empty()
        if(!('versions' in params.skip_qc)){
            GetSoftwareVersions()
        }
    emit:
        GetSoftwareVersions.out
} // end of wf_get_software_versions

/*
 * Parse software version numbers
 */
process GetSoftwareVersions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    output:
        // file 'software_versions_mqc.yaml', emit: yamlSoftwareVersion
        file 'software_versions_mqc.yaml'

    // when: !('versions' in skipQC)

    script:
    """
    init.sh
    bcftools version > v_bcftools.txt 2>&1 || true
    bwa &> v_bwa.txt 2>&1 || true
    configManta.py --version > v_manta.txt 2>&1 || true
    configureStrelkaGermlineWorkflow.py --version > v_strelka.txt 2>&1 || true
    echo "${workflow.manifest.version}" &> v_pipeline.txt 2>&1 || true
    echo "${workflow.nextflow.version}" &> v_nextflow.txt 2>&1 || true
    echo "SNPEFF version"\$(snpEff -h 2>&1) > v_snpeff.txt
    fastqc --version > v_fastqc.txt 2>&1 || true
    freebayes --version > v_freebayes.txt 2>&1 || true
    gatk ApplyBQSR --help 2>&1 | grep Version: > v_gatk.txt 2>&1 || true
    multiqc --version &> v_multiqc.txt 2>&1 || true
    qualimap --version &> v_qualimap.txt 2>&1 || true
    R --version &> v_r.txt  || true
    samtools --version &> v_samtools.txt 2>&1 || true
    tiddit &> v_tiddit.txt 2>&1 || true
    vcftools --version &> v_vcftools.txt 2>&1 || true
    vep --help &> v_vep.txt 2>&1 || true

    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}
