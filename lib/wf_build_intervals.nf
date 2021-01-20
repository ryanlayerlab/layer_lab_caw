include {hasExtension} from './utility'
workflow wf_build_intervals{
    take: fasta_fai
    main:
        // _ch_bed_intervals = Channel
        BuildIntervals(fasta_fai)

        _intervals = params.no_intervals ? "null" : \
                params.intervals && !('annotate' in params.step) ? \
                Channel.value(file(params.intervals)) : BuildIntervals.out
    
        CreateIntervalBeds(_intervals)
    emit:
        bed_intervals = CreateIntervalBeds.out.flatten()
}


// STEP 0: CREATING INTERVALS FOR PARALLELIZATION (PREPROCESSING AND VARIANT CALLING)
process BuildIntervals {
  label 'container_llab'
  tag {fastaFai}

  publishDir params.outdir

  input:
    file(fastaFai)

  output:
    file("${fastaFai.baseName}.bed")

  when: !(params.intervals) && !('annotate' in step) && !(params.no_intervals)

  script:
  """
  init.sh
  awk -v FS='\t' -v OFS='\t' '{ print \$1, \"0\", \$2 }' ${fastaFai} > ${fastaFai.baseName}.bed
  """
}

process CreateIntervalBeds {
  label 'container_llab'
    tag {intervals.fileName}

    input:
        file(intervals)

    output:
        file '*.bed'

    when: (!params.no_intervals) && step != 'annotate'

    script:
    // If the interval file is BED format, the fifth column is interpreted to
    // contain runtime estimates, which is then used to combine short-running jobs
    if (hasExtension(intervals, "bed"))
        """
        init.sh
        awk -vFS="\t" '{
          t = \$5  # runtime estimate
          if (t == "") {
            # no runtime estimate in this row, assume default value
            t = (\$3 - \$2) / ${params.nucleotides_per_second}
          }
          if (name == "" || (chunk > 600 && (chunk + t) > longest * 1.05)) {
            # start a new chunk
            name = sprintf("%s_%d-%d.bed", \$1, \$2+1, \$3)
            chunk = 0
            longest = 0
          }
          if (t > longest)
            longest = t
          chunk += t
          print \$0 > name
        }' ${intervals}
        """
    else if (hasExtension(intervals, "interval_list"))
        """
        init.sh
        grep -v '^@' ${intervals} | awk -vFS="\t" '{
          name = sprintf("%s_%d-%d", \$1, \$2, \$3);
          printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
        }'
        """
    else
        """
        init.sh
        awk -vFS="[:-]" '{
          name = sprintf("%s_%d-%d", \$1, \$2, \$3);
          printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
        }' ${intervals}
        """
}