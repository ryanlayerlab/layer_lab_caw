tools = params.globals.tools

workflow wf_hap_py{
    take: deepvariant_vcfs
    take: haplotypecaller_vcfs
    take: _target_bed
    take: _bait_bed
    take: _fasta
    take: _fasta_fai
    main:
    /* Take all vcf's from deepvariant and haplotypecaller, but filter out all but
        those that has 'GIAB' in the sample name.
        The truth set is going to be the GIAB vcfs obtained from their website 
    
    vcfs_giab structure: ['DeepVariant', 'GIABO1', 'GIABO1_S53', vcf.gz, vcf.gz.tbi]    
        */
    vcfs_giab = Channel.empty()
                    .mix(
                            deepvariant_vcfs,
                            haplotypecaller_vcfs
                        )
                    .filter{
                        // The sampleID is the 3rd component of the channel elements
                        "${it[2]}".contains('GIAB')
                    }
                    .dump(tag: "vcfs_giab")

    /*
    The vcfs_giab structure:
    ['Variant Caller', idPatiend, idSample, vcf.gz, vcf.gz.tbi]
    */

    /*
    hap_py_against_giab structure
    [benchmarked_against, 'source caller', idPatient, idSample, src.vcf.gz, src.vcf.gz.tbi, 
    ['GIAB', 'DeepVariant', 'GIABO1', 'GIABO1_S53', GIABO1_S53.vcf.gz, \
    highconf_PGandRTGphasetransfer.vcf.gz, \
    highconf_PGandRTGphasetransfer.vcf.gz.tbi, \
    highconf_nosomaticdel.bed]
    */
    hap_py_against_giab = Channel.from('GIAB') // this will be what the vcfs are being validated against
                            .combine(vcfs_giab) // above filtered vcfs
                            .combine(ch_giab_highconf_vcf) // truth set vcf
                            .combine(ch_giab_highconf_tbi) // truth set index
                            .combine(ch_giab_highconf_regions) // high quality regions from GIAB
                            .dump(tag: 'hap_py_against_giab')

    vcfs_dv_against_hc = deepvariant_vcfs
                      .join(haplotypecaller_vcfs, by:[1,2])                      
                      .map{
                            idPatient, idSample, val_dv, dv_vcf, dv_tbi, val_hc, hc_vcf, hc_tbi ->
                            [ val_dv, idPatient, idSample, dv_vcf, dv_tbi, hc_vcf, hc_tbi]
                            }
                    //   .dump(tag: 'vcfs_against_hc')
    hap_py_dv_against_hc =  Channel.from('hc') // this will be what the vcfs are being validated against
                            .combine(vcfs_dv_against_hc) // above filtered vcfs
                            .combine(ch_giab_highconf_regions) // high quality regions from GIAB
                            .dump(tag: 'hap_py_dv_against_hc')

    // For this channel to go ahead, both the haplotypecaller and the deepvariant should be in tools
    if (!(('benchmark_dv_against_hc' in tools) && ('benchmark_dv_and_hc_against_giab' in tools))) {
        hap_py_dv_against_hc = Channel.empty()
    }

    hap_py_combined = hap_py_against_giab.mix(hap_py_dv_against_hc)
                       .dump(tag:'hap_py_combined')

    HapPy(hap_py_combined,
          _target_bed,
          _bait_bed,
         _fasta,
         _fasta_fai
        )
    emit:
        HapPy.out
} // end of wf_hap_py

process HapPy {
    label 'cpus_32'

    tag {idSample}

    publishDir "${params.outdir}/Validation/Against_${truth_set_type}/${idSample}/${variantCaller}", mode: params.publish_dir_mode

    input:
        tuple truth_set_type, variantCaller, idPatient, idSample, file(vcf), file(tbi), file (truth_set_vcf), file (truth_set_tbi), file(highqual_regions)         
        file(target_bed)
        file(bait_bed)
        file(fasta)
        file(fastaFai)

    output:
        file("Against_${truth_set_type}.${vcf.baseName}.*")

    when: ( 
            ('benchmark_dv_and_hc_against_giab' in tools || 'benchmark_dv_against_hc' in tools )
            && ('haplotypecaller' in tools || 'deepvariant' in tools) 
            && params.giab_highconf_vcf 
            && params.giab_highconf_tbi 
            && params.giab_highconf_regions
        )

    script:
    // bn = "{vcf.baseName}"
    """
    init.sh
    export HGREF=$fasta
    mkdir scratch
    hap.py  \
        ${truth_set_vcf} \
        ${vcf} \
        -f ${highqual_regions} \
        --scratch-prefix scratch \
        --engine vcfeval \
        -T ${target_bed} \
        --threads ${task.cpus} \
        -o Against_${truth_set_type}.${vcf.baseName}
    """
}