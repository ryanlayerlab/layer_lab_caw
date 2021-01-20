tools =  params.globals.tools
/*
================================================================================
                                   ANNOTATION
================================================================================
*/
workflow wf_annotate{
    // to do
    take: _md_bam
    take: _fasta
    take: _fasta_fai
    main:
        /* Annotations */
        ch_vcfs_to_annotate = Channel.empty()
        if (step == 'annotate') {
            // ch_vcfs_to_annotate = getVCFsToAnnotate(params.outdir, annotate_tools)
            ch_vcfs_to_annotate = Channel.empty()
            // vcf_no_annotate = Channel.empty()

            if (tsv_path == [] || !tsv_path) {
            // Sarek, by default, annotates all available vcfs that it can find in the VariantCalling directory
            // Excluding vcfs from FreeBayes, and g.vcf from HaplotypeCaller
            // Basically it's: results/VariantCalling/*/{HaplotypeCaller,Manta,Mutect2,SentieonDNAseq,SentieonDNAscope,SentieonTNscope,Strelka,TIDDIT}/*.vcf.gz
            // Without *SmallIndels.vcf.gz from Manta, and *.genome.vcf.gz from Strelka
            // The small snippet `vcf.minus(vcf.fileName)[-2]` catches idSample
            // This field is used to output final annotated VCFs in the correct directory
                // log.info ("annotate_tools: ${annotate_tools}")
                ch_vcfs_to_annotate = 
                Channel.empty().mix(
                Channel.fromPath("${params.outdir}/VariantCalling/*/HaplotypeCaller/*.vcf.gz")
                    .flatten().map{vcf -> ['HaplotypeCaller', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
                // Channel.fromPath("${params.outdir}/VariantCalling/*/Manta/*[!candidate]SV.vcf.gz")
                Channel.fromPath("${params.outdir}/VariantCalling/*/Manta/*diploidSV.vcf.gz")
                    .flatten().map{vcf -> ['Manta', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
                Channel.fromPath("${params.outdir}/VariantCalling/*/Mutect2/*.vcf.gz")
                    .flatten().map{vcf -> ['Mutect2', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
                Channel.fromPath("${params.outdir}/VariantCalling/*/Strelka/*{somatic,variant}*.vcf.gz")
                    .flatten().map{vcf -> ['Strelka', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
                Channel.fromPath("${params.outdir}/VariantCalling/*/TIDDIT/*.vcf.gz")
                    .flatten().map{vcf -> ['TIDDIT', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
                )
                .filter {
                    annotate_tools == [] || (annotate_tools != [] && it[0] in annotate_tools)
                }
            } else if (annotate_tools == []) {
            // Annotate user-submitted VCFs
            // If user-submitted, Sarek assume that the idSample should be assumed automatically
            ch_vcfs_to_annotate = Channel.fromPath(tsv_path)
                .map{vcf -> ['userspecified', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
            } else exit 1, "specify only tools or files to annotate, not both"
        }
        
        // log.info "annotate_tools: ${annotate_tools}"
        // ch_vcfs_to_annotate.dump(tag: 'ch_vcf_to_annotate')
        ch_vcfs_to_annotate = ch_vcfs_to_annotate.mix(
                            ConcatVCF.out.vcf_concatenated_to_annotate,
                            StrelkaSingle.out.map{
                                variantcaller, idPatient, idSample, vcf, tbi ->
                                    [variantcaller, idSample, vcf[1]]
                            },
                            MantaSingle.out.map {
                                variantcaller, idPatient, idSample, vcf, tbi ->
                                [variantcaller, idSample, vcf[2]]
                            },
                            TIDDIT.out.vcfTIDDIT.map {
                                variantcaller, idPatient, idSample, vcf, tbi ->
                                [variantcaller, idSample, vcf]
                                }
                            )

        // ch_vcf_snpEff = ch_vcfs_to_annotate.mix(ConcatVCF.out.vcf_concatenated_to_annotate)
        ch_vcf_snpEff = ch_vcfs_to_annotate
        // ch_vcf_snpEff = Channel.empty()

    ch_vcf_vep = ch_vcf_snpEff.map {
                variantCaller, idSample, vcf ->
                [variantCaller, idSample, vcf, null]
        }
    //     // PrintCh(ch_vcf_snpeff)
        SnpEff( ch_vcf_snpEff,
                ch_snpEff_cache,
                ch_snpEff_db
                )
        CompressVCFsnpEff(SnpEff.out.snpEff_vcf)
        
        VEP(ch_vcf_vep,
            ch_vep_cache,
            ch_vep_cache_version,
            ch_cadd_InDels,
            ch_cadd_InDels_tbi,
            ch_cadd_WG_SNVs,
            ch_cadd_WG_SNVs_tbi,
            _fasta,
            _fasta_fai
        )

        VEPmerge(CompressVCFsnpEff.out.compressVCFsnpEffOut,
            ch_vep_cache,
            ch_vep_cache_version,
            ch_cadd_InDels,
            ch_cadd_InDels_tbi,
            ch_cadd_WG_SNVs,
            ch_cadd_WG_SNVs_tbi,
            _fasta,
            _fasta_fai
        )
        CompressVCFvep(VEP.out.vep_vcf.mix(
                        VEPmerge.out.vep_vcf_merge)
                        )
} // end of wf_germline_cnv
// STEP SNPEFF

process SnpEff {
    label 'container_llab'
    tag {"${idSample} - ${variantCaller} - ${vcf}"}
    // cache false
    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
        if (it == "${reducedVCF}_snpEff.ann.vcf") null
        else "Reports/${idSample}/snpEff/${variantCaller}/${it}"
    }

    input:
        tuple variantCaller, idSample, file(vcf) 
        file(dataDir)
        // path(dataDir)
        val snpeffDb

    output:
        tuple file("${reducedVCF}_snpEff.txt"), file("${reducedVCF}_snpEff.html"), file("${reducedVCF}_snpEff.csv"), emit:snpEff_report
        tuple variantCaller, idSample, file("${reducedVCF}_snpEff.ann.vcf"), emit: snpEff_vcf

    when: 'snpeff' in tools || 'merge' in tools

    script:
    reducedVCF = reduceVCF(vcf.fileName)
    cache = (params.snpEff_cache && params.annotation_cache) ? "-dataDir \${PWD}/${dataDir}" : ""
    // cache = (params.snpeff_cache && params.annotation_cache) ? "-dataDir ${dataDir}" : ""
    """
    init.sh
    snpEff -Xmx${task.memory.toGiga()}g \
        ${snpeffDb} \
        -csvStats ${reducedVCF}_snpEff.csv \
        -nodownload \
        ${cache} \
        -canon \
        -v \
        ${vcf} \
        > ${reducedVCF}_snpEff.ann.vcf

    mv snpEff_summary.html ${reducedVCF}_snpEff.html
    mv ${reducedVCF}_snpEff.genes.txt ${reducedVCF}_snpEff.txt
    """
}

// snpeffReport = snpeffReport.dump(tag:'snpEff report')

// STEP COMPRESS AND INDEX VCF.1 - SNPEFF

process CompressVCFsnpEff {
    label 'container_llab'
    tag {"${idSample} - ${vcf}"}

    publishDir "${params.outdir}/Annotation/${idSample}/snpEff/${variantCaller}", mode: params.publish_dir_mode

    input:
        tuple variantCaller, idSample, file(vcf)

    output:
        tuple variantCaller, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi"), emit: compressVCFsnpEffOut

    script:
    """
    init.sh
    bgzip < ${vcf} > ${vcf}.gz
    tabix ${vcf}.gz
    """
}