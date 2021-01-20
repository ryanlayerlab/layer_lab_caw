tools = params.globals.tools
model = params.model
workflow wf_deepvariant{
    take: _dm_bams
    take: _target_bed
    take: _fasta
    take: _fasta_fai
    take: _fasta_gz
    take: _fasta_gz_fai
    take: _fasta_gzi

    main:
        DV_MakeExamples(_dm_bams,
                _fasta,
                _fasta_fai,
                _fasta_gz,
                _fasta_gz_fai,
                _fasta_gzi,
                _target_bed)
        DV_CallVariants(DV_MakeExamples.out.shared_examples)
        DV_PostprocessVariants(
            DV_CallVariants.out.variant_tf_records,
            _fasta_gz,
            _fasta_gz_fai,
            _fasta_gzi
        )

    emit:
        vcf = DV_PostprocessVariants.out.vcf
        vcf_to_annotate = DV_PostprocessVariants.out.vcf_to_annotate
} // end of wf_deepvariant

/*
================================================================================
                               Google Deepvariant  for Germline variant calling
================================================================================
*/

/********************************************************************
  process make_examples
  Getting bam files and converting them to images ( named examples )
********************************************************************/

process DV_Combined{
    label "cpus_max"
    tag "${bam}"
    publishDir "${params.outdir}/VariantCalling/${idSample}/DeepVariant", mode: params.publish_dir_mode
    
    input:
    tuple idPatient, idSample, file(bam), file(bai)
    file (fasta)
    file (fai )
    file (fastagz)
    file (gzfai)
    file (gzi)
    file (target_bed)
    
    output:
        tuple idPatient, idSample, file("${bam.baseName}.vcf")
        file("*.html")

    when: 'deepvariant' in tools
    
    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
    --model_type ${model} \
    --ref ${fasta} \
    --reads ${bam} \
     --regions ${target_bed} \
    --output_vcf "${bam.baseName}.vcf" \
    --num_shards ${task.cpus}
    """
}

process DV_MakeExamples{
    label 'cpus_max'
    tag "${bam}"
    publishDir "${params.outdir}/Preprocessing/${idSample}/DV_MakeExamples/", mode: params.publish_dir_mode,
    saveAs: {filename -> "logs/log"}

    input:
    tuple idPatient, idSample, file(bam), file(bai)
    file (fasta)
    file (fai )
    file (fastagz)
    file (gzfai)
    file (gzi)
    file (target_bed)

    output:
    tuple idPatient, idSample, file(bam), file('*_shardedExamples'), emit: shared_examples
    when: 'deepvariant' in tools
    script:
    // sharded_folder = "${bam}.tfrecord@${task.cpus}.gz"
    use_regions = params.target_bed ? "--regions ${target_bed}" :  ''
    """
    mkdir logs
    mkdir ${bam.baseName}_shardedExamples
    dv_make_examples.py \
    --cores ${task.cpus} \
    --sample ${bam} \
    --ref ${fastagz} \
    --reads ${bam} \
    ${use_regions} \
    --logdir logs \
    --examples ${bam.baseName}_shardedExamples

    """
}

/********************************************************************
  process call_variants
  Doing the variant calling based on the ML trained model.
********************************************************************/

process DV_CallVariants{
   label 'cpus_max'
  tag "${bam}"

  input:
  tuple idPatient, idSample, file(bam), file(shardedExamples)

  output:
  tuple idPatient, idSample, file(bam), file('*_call_variants_output.tfrecord'), emit: variant_tf_records
  

  when: 'deepvariant' in tools
  script:
  """
    dv_call_variants.py \
    --cores ${task.cpus} \
    --sample ${bam} \
    --outfile ${bam.baseName}_call_variants_output.tfrecord \
    --examples $shardedExamples \
    --model ${model}

  """
}



/********************************************************************
  process postprocess_variants
  Trasforming the variant calling output (tfrecord file) into a standard vcf file.
********************************************************************/

process DV_PostprocessVariants{

  label 'cpus_32'
  tag "${bam}"
  publishDir "${params.outdir}/VariantCalling/${idSample}/DeepVariant", mode: params.publish_dir_mode

  input:
  tuple idPatient, idSample, file(bam), file(call_variants_tfrecord)
//   set file(bam),file('call_variants_output.tfrecord') from called_variants
  file fastagz
  file gzfai
  file gzi

  output:
   tuple val('DeepVariant'), idPatient, idSample, file("${idSample}.vcf.gz"), file("${idSample}.vcf.gz.tbi") , emit: vcf
   tuple val('DeepVariant'), idSample, file("${idSample}.vcf.gz") , emit: vcf_to_annotate
//    file("${idSample}.*.html"), emit: html_report
   file('*.html')

  when: 'deepvariant' in tools
  script:
  """
  dv_postprocess_variants.py \
  --ref ${fastagz} \
  --infile ${call_variants_tfrecord} \
  --outfile "${idSample}.vcf" \

  bgzip "${idSample}.vcf"
  tabix -p vcf "${idSample}.vcf.gz"
  """
}
