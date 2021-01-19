// STEP STRELKA.1 - SINGLE MODE
tools =  params.globals.tools
process StrelkaSingle {
    label 'cpus_max'
    // label 'memory_max'

    tag {idSample}

    publishDir "${params.outdir}/VariantCalling/${idSample}/Strelka", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, file(bam), file(bai)
        file(fasta)
        file(fastaFai)
        file(targetBED)

    output:
        set val("Strelka"), idPatient, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi")

    when: 'strelka' in tools

    script:
    beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
    options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
    """
    init.sh
    ${beforeScript}
    configureStrelkaGermlineWorkflow.py \
        --bam ${bam} \
        --referenceFasta ${fasta} \
        ${options} \
        --runDir Strelka

    python Strelka/runWorkflow.py -m local -j ${task.cpus}

    mv Strelka/results/variants/genome.*.vcf.gz \
        Strelka_${idSample}_genome.vcf.gz
    mv Strelka/results/variants/genome.*.vcf.gz.tbi \
        Strelka_${idSample}_genome.vcf.gz.tbi
    mv Strelka/results/variants/variants.vcf.gz \
        Strelka_${idSample}_variants.vcf.gz
    mv Strelka/results/variants/variants.vcf.gz.tbi \
        Strelka_${idSample}_variants.vcf.gz.tbi
    """
}
