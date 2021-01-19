
tools =  params.globals.tools
workflow wf_somatic_pon{
    take:  _vcfs_normal
    take: _fasta
    take: _fasta_fai
    take: _dict
    take: _germline_resource
    take: _germline_resource_index
    take: _target_bed

    main:
        SomaticPonGenomicsDBImport(
            _vcfs_normal,
            _target_bed
        )

        CreateSomaticPON(
            SomaticPonGenomicsDBImport.out,
            _fasta,
            _fasta_fai,
            _dict,
            _germline_resource,
            _germline_resource_index
            )
    emit:
        
        somatic_pon_gdb = SomaticPonGenomicsDBImport.out
        somatic_pon = CreateSomaticPON.out
        
} // end of wf_somatic_pon
// STEP GATK GenomicsDBImport
process SomaticPonGenomicsDBImport {
    label 'cpus_32'

    publishDir "${params.outdir}/Preprocessing/Somatic_pon_db", mode: params.publish_dir_mode

    input:
    file("vcfs/*")
    file(targetBED)

    output:
    file("somatic_pon.gdb")

    when: 'gen_somatic_pon' in tools

    script:
    sample_map="cohort_samples.map"
    
    // gDB = chr
    """
    init.sh
    vcfs=' '
    for x in `ls vcfs/*.vcf.gz`
    do
        base_name=`basename \${x}`
        without_ext=\${base_name%.vcf.gz}
        sample_name=\${without_ext##*_}
        echo "\${sample_name}\t\$x" >> $sample_map 
    done

    gatk --java-options -Xmx${task.memory.toGiga()}g \
    GenomicsDBImport  \
    --genomicsdb-workspace-path somatic_pon.gdb \
    -L ${targetBED} \
    --sample-name-map $sample_map \
    --merge-input-intervals \
    --reader-threads ${task.cpus}
    """
}

process CreateSomaticPON{
    label 'cpus_max'
    // label 'memory_max'
     publishDir "${params.outdir}/Preprocessing/Somatic_pon", mode: params.publish_dir_mode

    input:
    file(pon) 
    file(fasta)
    file(fastaFai)
    file(dict)
    file(germlineResource)
    file(germlineResourceIndex)
    
    output:
    tuple file(out_file), file ("${out_file}.tbi")

    when: 'gen_somatic_pon' in tools

    script:
    args_file = "normals_for_pon_vcf.args"
    out_file = "somatic_pon.vcf.gz" 
    pon_db = "gendb://${pon}"
    
    """
    init.sh
     gatk --java-options -Xmx${task.memory.toGiga()}g \
     CreateSomaticPanelOfNormals -R ${fasta} \
     --germline-resource ${germlineResource} \
    -V ${pon_db} \
    -O ${out_file}
    """
}