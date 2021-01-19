tools =  params.globals.tools
workflow wf_tiddit_gm_sv{
    take: _md_bam
    take: _fasta
    take: _fasta_fai
    main:
        TIDDIT(
            _md_bam,
            _fasta,
            _fasta_fai
        )
} // end of wf_germline_cnv

// STEP TIDDIT

process TIDDIT {
    tag {idSample}

    publishDir "${params.outdir}/VariantCalling/${idSample}/TIDDIT", mode: params.publish_dir_mode

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (it == "TIDDIT_${idSample}.vcf") "VariantCalling/${idSample}/TIDDIT/${it}"
            else "Reports/${idSample}/TIDDIT/${it}"
        }

    input:
        tuple idPatient, idSample, file(bam), file(bai)
        file(fasta) 
        file(fastaFai)

    output:
        tuple val("TIDDIT"), idPatient, idSample, file("*.vcf.gz"), file("*.tbi"), emit: vcfTIDDIT
        tuple file("TIDDIT_${idSample}.old.vcf"), file("TIDDIT_${idSample}.ploidy.tab"), file("TIDDIT_${idSample}.signals.tab"), file("TIDDIT_${idSample}.wig"), file("TIDDIT_${idSample}.gc.wig"), emit: tidditOut

    when: 'tiddit' in tools

    script:
    """
    tiddit --sv -o TIDDIT_${idSample} --bam ${bam} --ref ${fasta}

    mv TIDDIT_${idSample}.vcf TIDDIT_${idSample}.old.vcf

    grep -E "#|PASS" TIDDIT_${idSample}.old.vcf > TIDDIT_${idSample}.vcf

    bgzip --threads ${task.cpus} -c TIDDIT_${idSample}.vcf > TIDDIT_${idSample}.vcf.gz

    tabix TIDDIT_${idSample}.vcf.gz
    """
}