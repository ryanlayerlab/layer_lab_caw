tools = params.globals.tools

include {ConcatVCF} from './utility'
workflow wf_jointly_genotype_gvcf{
    // input _scattered_gvcf_HC structure 
    // tuple idPatient, idSample, file(intervalBed), file("${intervalBed.baseName}_${idSample}.g.vcf")
    take: _scattered_gvcf_HC
    take: _target_bed
     take: _fasta
    take: _fasta_fai
    take: _dict
    take: _dbsnp
    take: _dbsnp_index
    // take: _bed_intervls
    main:
    // Transform the haplotypecaller output in the above format to a format where gvcf's are gathered per interval, so carry
    // multiple gvcf (for multiple samples)
    // When grouping gvcf's for per interval basis, we only keep one interval bed as the list carry same intervals for each group
    // This interval then become the key, and a disnguishing element for this joint genotyping
    // patientSampleIdMap will have all the samples names that has a gvcf in a particular group

    mapped_gvcf_GenotypeGVCFs = _scattered_gvcf_HC
                                .map{idPatient, idSample, intervalBed, gvcf ->
                                patientSampleIdMap = [:]
                                patientSampleIdMap['idPatient'] = idPatient
                                patientSampleIdMap['idSample'] = idSample
                                [intervalBed.baseName, intervalBed, patientSampleIdMap , gvcf]}
                                .groupTuple(by:[0])
                                .map{interval_name, liIntervalBed, patientSampleIdMap, liGvcf -> 
                                    [interval_name, liIntervalBed.first(), patientSampleIdMap, liGvcf]
                                    }
                                // .dump(tag: 'Collected HaplotypeCaller output')
    
    
    GenomicsDBImport(mapped_gvcf_GenotypeGVCFs)                              
    GenotypeGVCFs(
        GenomicsDBImport.out,
        _fasta,
        _fasta_fai,
        _dict,
        _dbsnp,
        _dbsnp_index
    )
    // GenotypeGVCFs output
    //tuple val("HaplotypeCaller"),  val(patientSampleIdMap), file(interval_bed), file("vcf"), file("vcf.idx")
    
    // A Cohort vcf
    vcf_cohort_concatenate_vcf = GenotypeGVCFs.out.vcf_GenotypeGVCFs
    .map{caller, li_patient_sample_id_map, interval_bed,  vcf, vcf_idx ->
        [vcf]
    }
    .collect()
    // .dump(tag: 'cohor_vcf: ')

    CohortConcatVCF(
        vcf_cohort_concatenate_vcf,
        _fasta_fai,
        _target_bed,
    )

    // Per Sample vcf (but jointly genotyped)
     ch_select_variants = GenotypeGVCFs.out.vcf_GenotypeGVCFs
    .flatMap{ caller, li_patient_sample_id_map, interval_bed,  vcf, vcf_idx ->
        per_sample_list=[]
        li_patient_sample_id_map.each { entry ->
            per_sample_list.add([caller, entry.idPatient, entry.idSample, interval_bed, vcf, vcf_idx])
        }
        per_sample_list
    }
    // .flatten()
    // .dump(tag: 'ch_select_variants')

    SelectVariants(
        ch_select_variants,
        _fasta,
        _fasta_fai,
        _dict
    )
    // SelectVariants output
    //tuple val("HaplotypeCaller_Jointly_Genotyped"), id_patient, id_sample, file("vcf")
    vcf_ConcatenateVCFs = SelectVariants.out.vcf_SelectVariants.groupTuple(by:[0, 1, 2])
    // Now add the individually called vcfs too
    // vcf_ConcatenateVCFs = vcf_ConcatenateVCFs.mix(vcf_HaplotypeCaller)
    // if (!params.no_gvcf){ // if user specified, noGVCF, skip saving the GVCFs from HaplotypeCaller
    //     vcf_ConcatenateVCFs = vcf_ConcatenateVCFs.mix(gvcf_HaplotypeCaller)
    // }
    // vcf_ConcatenateVCFs.dump('concat_vcf: ')

    ConcatVCF(
        vcf_ConcatenateVCFs,
        _fasta_fai,
        _target_bed,
        'HC_jointly_genotyped', // prefix for output files
        'vcf', // extension for the output files
        'HC_jointly_genotyped_gvcf' // output directory name
        )
    
    // Create a channel to hold GIAB samples for validation using hap.py
    // hc_jointly_genotyped_vcfs = ConcatVCF.out.concatenated_vcf_with_index
    //                       .filter{  "${it[0]}" == 'HaplotypeCaller_Jointly_Genotyped'}
    //                     //   .dump(tag: 'hc_jointly_genotyped_vcfs')

    emit:
    vcfs_with_indexes = ConcatVCF.out.concatenated_vcf_with_index
    vcfs_without_indexes = ConcatVCF.out.concatenated_vcf_without_index
    // cohort_vcf_with_index = CohortConcatVCF.out.cohort_vcf_with_index
    // cohort_vcf_without_index = CohortConcatVCF.out[1]
} // end of wf_haplotypecaller

process GenomicsDBImport {
    label 'container_llab'
    label 'cpus_16'
    // echo true
    tag{interval_name}
    // publishDir "${OUT_DIR}/misc/genomicsdb/", mode: 'copy', overwrite: false

    input:
    // tuple val(interval_name), file(interval_bed), val(list_id_patient), val(list_id_sample), file(gvcfs)
    tuple val(interval_name), file(interval_bed), val(patientSampleIdMap), file(gvcfs)
    
    output:
    tuple val(interval_name), file(interval_bed), val(patientSampleIdMap), file ("${interval_name}.gdb")

    when: 'haplotypecaller' in tools

    script:
    sample_map="cohort_samples.map"
    interval_name_with_underscore="${interval_name}_"
    // gDB = chr
    """
    init.sh
    for x in *.g.vcf
    do
        bgzip \$x
        tabix \${x}.gz
    done

    for x in *.g.vcf.gz
    do
        
        base_name=`basename \$x .g.vcf.gz`
        sample=\${base_name#$interval_name_with_underscore}
        echo "\${sample}\t\${x}" >> ${sample_map}
    done
    
    gatk --java-options -Xmx${task.memory.toGiga()}g \
    GenomicsDBImport \
    --genomicsdb-workspace-path ${interval_name}.gdb \
    -L $interval_bed \
    --sample-name-map ${sample_map} \
    --reader-threads ${task.cpus}

    """
}

// STEP GATK HAPLOTYPECALLER.2


process GenotypeGVCFs {
    label 'container_llab'
    label 'cpus_8'
    tag {interval_bed.baseName}
    input:
        // tuple val(interval_name), file(interval_bed), val(list_id_patient), val(list_id_sample), file(gdb)
        tuple val(interval_name), file(interval_bed), val(patientSampleIdMap), file(gdb)
        file(fasta)
        file(fastaFai)
        file(dict)
        file(dbsnp)
        file(dbsnpIndex)

    output:
    // tuple val("HaplotypeCaller"), list_id_patient, list_id_sample, file("${interval_name}.vcf"), file("${interval_name}.vcf.idx"), emit: vcf_GenotypeGVCFs
    tuple val("HaplotypeCaller"),  val(patientSampleIdMap), file(interval_bed), file("${interval_name}.vcf"), file("${interval_name}.vcf.idx"), emit: vcf_GenotypeGVCFs
    // tuple val("HaplotypeCaller_Jointly_Genotyped"), val(interval_name), file(interval_bed), val(list_id_patient), val(list_id_sample), file ("${interval_name}.vcf"), file ("${interval_name}.vcf.idx"), emit: vcf_GenotypeGVCFs
    
    when: 'haplotypecaller' in tools

    script:
    // Using -L is important for speed and we have to index the interval files also
    """
    init.sh
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        GenotypeGVCFs \
        -R ${fasta} \
        -L ${interval_bed} \
        -D ${dbsnp} \
        -V gendb://${gdb} \
        --create-output-variant-index \
        -O "${interval_name}.vcf"
    """
}

process SelectVariants {
    label 'container_llab'
    label 'cpus_8'
    tag {interval_bed.baseName}
    input:
        // tuple val(caller), val(id_patient), val(id_sample), val(interval_name), file(interval_bed), file (vcf), file (vcf_idx)
        tuple val(caller), val(id_patient), val(id_sample), file(interval_bed), file (vcf), file (vcf_idx)
        file(fasta)
        file(fastaFai)
        file(dict)

    output:
    // tuple val("HaplotypeCaller"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.vcf"), emit: vcf_GenotypeGVCFs
    // tuple val("HaplotypeCaller"), val(interval_name), file(interval_bed), val(list_id_patient), val(list_id_sample), file ("${interval_name}.vcf"), file ("${interval_name}.vcf.idx"), emit: vcf_GenotypeGVCFs
    tuple val("HaplotypeCaller_Jointly_Genotyped"), id_patient, id_sample, file("${interval_bed.baseName}_${id_sample}.vcf"), emit: vcf_SelectVariants
    
    when: 'haplotypecaller' in tools

    script:
    // Using -L is important for speed and we have to index the interval files also
    """
    init.sh
    gatk --java-options -Xmx${task.memory.toGiga()}g \
            SelectVariants \
            -R ${fasta} \
            -L ${interval_bed} \
            -V ${vcf} \
            -O ${interval_bed.baseName}_${id_sample}.vcf \
            -sn ${id_sample}
    """
}



process CohortConcatVCF {
    label 'container_llab'
    label 'cpus_8'

    tag {'CohortConcatVCF'}

    publishDir "${params.outdir}/VariantCalling/HC_cohort_vcf", mode: params.publish_dir_mode

    input:
        file(vcFiles)
        file(fastaFai)
        file(targetBED)

    output:
    // we have this funny *_* pattern to avoid copying the raw calls to publishdir
        tuple file("HC_cohort.vcf.gz"), file("HC_cohort.vcf.gz.tbi"), emit: cohort_vcf_with_index
        // tuple file("HC_cohort.vcf.gz"), emit: cohort_vcf_without_index
        file("HC_cohort.vcf.gz")
        // file("HC_cohort.vcf.gz"), file("HC_cohort.vcf.gz.tbi"), emit: cohortvcfwithindex
        // file("HC_cohort.vcf.gz"), emit: cohortvcfwithoutindex

    when: ('haplotypecaller' in tools || 'mutect2' in tools || 'freebayes' in tools)

    script:
    options = params.target_bed ? "-t ${targetBED}" : ""
    """
    init.sh
    concatenateVCFs.sh -i ${fastaFai} -c ${task.cpus} -o HC_cohort.vcf ${options}
    """
}

