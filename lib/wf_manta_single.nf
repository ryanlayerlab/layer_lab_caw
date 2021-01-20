
tools = params.globals.tools

workflow wf_germline_cnv{
    // deepvariant output
    //tuple val('DeepVariant'), idSample, file("${idSample}.vcf.gz")
    take: _raw_bam
    take: _target_bed
    take: _fasta
    take: _fasta_fai
    take: _cnvkit_ref
    
    main:
    //     StrelkaSingle(
    //         _bam_recal,
    //         _fasta,
    //         _fasta_fai,
    //         _target_bed
    // )
        MarkDuplicates(
            _raw_bam
        )
        // CNVKitSingle(
        //      MarkDuplicates.out,
        //     _fasta,
        //     _fasta_fai,
        //     _cnvkit_ref
        // )
        MantaSingle(
            MarkDuplicates.out,
            _fasta,
            _fasta_fai,
            _target_bed
        )
} // end of wf_germline_cnv


// STEP MANTA.1 - SINGLE MODE

process MantaSingle {
    label 'container_sarek'
    label 'cpus_32'
    // label 'memory_max'

    tag {idSample}

    publishDir "${params.outdir}/VariantCalling/${idSample}/Manta", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, file(bam), file(bai)
        file(fasta)
        file(fastaFai)
        file(targetBED)

    output:
        tuple val("Manta"), idPatient, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi")

    when: 'manta' in tools

    script:
    beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
    options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
    status = status_map[idPatient, idSample]
    inputbam = status == 0 ? "--bam" : "--tumorBam"
    vcftype = status == 0 ? "diploid" : "tumor"
    """
    ${beforeScript}
    configManta.py \
        ${inputbam} ${bam} \
        --reference ${fasta} \
        ${options} \
        --runDir Manta

    python Manta/runWorkflow.py -m local -j ${task.cpus}

    mv Manta/results/variants/candidateSmallIndels.vcf.gz \
        Manta_${idSample}.candidateSmallIndels.vcf.gz
    mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        Manta_${idSample}.candidateSmallIndels.vcf.gz.tbi
    mv Manta/results/variants/candidateSV.vcf.gz \
        Manta_${idSample}.candidateSV.vcf.gz
    mv Manta/results/variants/candidateSV.vcf.gz.tbi \
        Manta_${idSample}.candidateSV.vcf.gz.tbi
    mv Manta/results/variants/${vcftype}SV.vcf.gz \
        Manta_${idSample}.${vcftype}SV.vcf.gz
    mv Manta/results/variants/${vcftype}SV.vcf.gz.tbi \
        Manta_${idSample}.${vcftype}SV.vcf.gz.tbi
    """
}