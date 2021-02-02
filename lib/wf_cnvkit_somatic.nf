/* CNVKit related processes */
// CNVKIT related
tools = params.globals.tools
status_map = params.globals.status_map

workflow wf_cnvkit_somatic{
    take: _bams_normal
    take: _bais_normal
    take: _bams_tumor
    take: _bais_tumor
    take: _fasta
    take: _fasta_fai
    take: _target_bed
    
    main:
   CNVKitSomatic(_bams_normal,
            _bais_normal,
            _bams_tumor,
            _bais_tumor,
            _fasta,
            _fasta_fai,
            _target_bed)
} // end of wf_germline_cnv

process CNVKitSomatic{
    label 'container_llab'
    label 'cpus_32'
    publishDir "${params.outdir}/VariantCalling/", mode: params.publish_dir_mode
  
  input:
    file(bams_normal)
    file(bais_normal)
    file(bams_tumor)
    file(bais_tumor)
    file(fasta)
    file(fastaFai)
    file(targetBED)
  
  output:
  file('CNVKit_somatic_batch')

  when:'cnvkit_somatic' in tools

  script:
  """
  init.sh
  cnvkit.py batch -p32 ${bams_tumor} \
      --normal ${bams_normal} \
      --targets ${targetBED} \
      --fasta ${fasta}  \
      --output-reference my_reference.cnn \
      --output-dir CNVKit_somatic_batch \
      --scatter \
      --diagram
  """
}

// process CNVKitSingle{
//     label 'container_llab'
//     label 'cpus_8'
//     publishDir "${params.outdir}/VariantCalling/${idSample}/CNVKit", mode: params.publish_dir_mode
    
//     input:
//         tuple idPatient, idSample, file(bam), file(bai)
//         file(fasta)
//         file(fastaFai)
//         file(cnvkit_ref)
//         // file(targetBED)
    
//     output:
//     tuple val("cnvkit_single"), idPatient, idSample, file("*")

//     when:'cnvkit_somatic' in tools

//     script:
    
//     """
//     init.sh
//     cnvkit.py batch ${bam} \
//         --reference ${cnvkit_ref} \
//         --scatter \
//         --diagram
//     """
// }
