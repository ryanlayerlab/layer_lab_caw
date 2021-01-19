tools =  params.globals.tools
skip_qc =  params.globals.skip_qc
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
// STEP VEP.1

process VEP {
    label 'VEP'
    label 'cpus_4'

    tag {"${idSample} - ${variantCaller} - ${vcf}"}

    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
        if (it == "${reducedVCF}_VEP.summary.html") "Reports/${idSample}/VEP/${variantCaller}/${it}"
        else null
    }

    input:
        tuple variantCaller, idSample, file(vcf), file(idx) 
        file(dataDir) 
        val cache_version 
        file(cadd_InDels) 
        file(cadd_InDels_tbi) 
        file(cadd_WG_SNVs) 
        file(cadd_WG_SNVs_tbi) 
        file(fasta)
        file(fasta_fai)
    output:
        tuple variantCaller, idSample, file("${reducedVCF}_VEP.ann.vcf"), emit: vep_vcf
        // file("${reducedVCF}_VEP.summary.html") , emit: vep_report
        file("${reducedVCF}_VEP.summary.html") 

    when: ('vep' in tools) && params.cadd_InDels && params.cadd_WG_SNVs

    script:
    reducedVCF = reduceVCF(vcf.fileName)
    genome = params.genome == 'smallGRCh37' ? 'GRCh37' : params.genome

    dir_cache = (params.vep_cache && params.annotation_cache) ? " \${PWD}/${dataDir}" : "/.vep"
    // cadd = (params.cadd_cache && params.cadd_WG_SNVs && params.cadd_InDels) ? "--plugin CADD,whole_genome_SNVs.tsv.gz,InDels.tsv.gz" : ""
    cadd = (params.cadd_WG_SNVs && params.cadd_InDels) ? "--plugin CADD,whole_genome_SNVs.tsv.gz,InDels.tsv.gz" : ""
    genesplicer = params.genesplicer ? "--plugin GeneSplicer,/opt/miniconda/envs/layer_lab_dna_seq/bin/genesplicer,/opt/miniconda/envs/layer_lab_dna_seq/share/genesplicer-1.0-1/human,context=200,tmpdir=\$PWD/${reducedVCF}" : "--offline"
    """
    init.sh
    mkdir ${reducedVCF}

    vep \
        -i ${vcf} \
        -o ${reducedVCF}_VEP.ann.vcf \
        --assembly ${genome} \
        --species ${params.species} \
        ${cadd} \
        ${genesplicer} \
        --cache \
        --cache_version ${cache_version} \
        --dir_cache ${dir_cache} \
        --everything \
        --filter_common \
        --fork ${task.cpus} \
        --format vcf \
        --per_gene \
        --stats_file ${reducedVCF}_VEP.summary.html \
        --total_length \
        --vcf \
        --offline
    rm -rf ${reducedVCF}
    """
}

// vepReport = vepReport.dump(tag:'VEP')

// STEP VEP.2 - VEP AFTER SNPEFF

process VEPmerge {
    label 'VEP'
    label 'cpus_4'

    tag {"${idSample} - ${variantCaller} - ${vcf}"}

    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
        if (it == "${reducedVCF}_VEP.summary.html") "Reports/${idSample}/VEP/${variantCaller}/${it}"
        else null
    }

    input:
        tuple variantCaller, idSample, file(vcf), file(idx) 
        file(dataDir) 
        val cache_version 
        file(cadd_InDels) 
        file(cadd_InDels_tbi) 
        file(cadd_WG_SNVs) 
        file(cadd_WG_SNVs_tbi) 
        file(fasta)
        file(fasta_fai)
    output:
        tuple variantCaller, idSample, file("${reducedVCF}_VEP.ann.vcf"), emit: vep_vcf_merge
        // file("${reducedVCF}_VEP.summary.html") , emit: vep_report_merge
        file("${reducedVCF}_VEP.summary.html") 

    when: ('merge' in tools) && params.cadd_InDels && params.cadd_WG_SNVs

    script:
    reducedVCF = reduceVCF(vcf.fileName)
    genome = params.genome == 'smallGRCh37' ? 'GRCh37' : params.genome
    dir_cache = (params.vep_cache && params.annotation_cache) ? " \${PWD}/${dataDir}" : "/.vep"
    // cadd = (params.cadd_cache && params.cadd_WG_SNVs && params.cadd_InDels) ? "--plugin CADD,whole_genome_SNVs.tsv.gz,InDels.tsv.gz" : ""
    cadd = (params.cadd_WG_SNVs && params.cadd_InDels) ? "--plugin CADD,whole_genome_SNVs.tsv.gz,InDels.tsv.gz" : ""
    genesplicer = params.genesplicer ? "--plugin GeneSplicer,/opt/miniconda/envs/layer_lab_dna_seq/bin/genesplicer,/opt/miniconda/envs/layer_lab_dna_seq/share/genesplicer-1.0-1/human,context=200,tmpdir=\$PWD/${reducedVCF}" : "--offline"
    """
    init.sh
    mkdir ${reducedVCF}

    vep \
        -i ${vcf} \
        -o ${reducedVCF}_VEP.ann.vcf \
        --assembly ${genome} \
        --species ${params.species} \
        ${cadd} \
        ${genesplicer} \
        --cache \
        --cache_version ${cache_version} \
        --dir_cache ${dir_cache} \
        --everything \
        --filter_common \
        --fork ${task.cpus} \
        --format vcf \
        --per_gene \
        --stats_file ${reducedVCF}_VEP.summary.html \
        --total_length \
        --vcf \
        --offline

    rm -rf ${reducedVCF}
    """
}


// STEP COMPRESS AND INDEX VCF.2 - VEP

process CompressVCFvep {
    tag {"${idSample} - ${vcf}"}

    publishDir "${params.outdir}/Annotation/${idSample}/VEP/${variantCaller}", 
    mode: params.publish_dir_mode

    input:
        tuple variantCaller, idSample, file(vcf) 

    output:
        tuple variantCaller, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") 

    script:
    """
    init.sh
    bgzip < ${vcf} > ${vcf}.gz
    tabix ${vcf}.gz
    """
}