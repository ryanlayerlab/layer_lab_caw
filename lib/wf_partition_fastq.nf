step = tools = params.globals.step
workflow wf_partition_fastq{
    main:
    _pair_reads = Channel.empty()
    // Close the input_pair_read if the starting step is not mapping
    if (step == 'mapping'){           
        if (params.split_fastq){
            PartitionFastQ(ch_input_sample)
            _pair_reads = _pair_reads.mix(PartitionFastQ.out)
            .flatMap{            
                idPatient, idSample, idRun, reads_1, reads_2 ->
                myList= []
                reads_1.each { read_1 -> 
                                split_index = read_1.fileName.toString().minus("r1_split_").minus(".fastq.gz")
                                parent = read_1.parent
                                read_2_fn = read_1.fileName.toString().replace("r1_split_", "r2_split_")
                                read_2 = "${parent}/${read_2_fn}"
                                new_id_run = "${idRun}_${split_index}"
                                myList.add([idPatient, idSample, new_id_run, read_1, file(read_2)])
                            }
                myList     
            }
        }else{
            _pair_reads = ch_input_sample
        }
    }
    emit:
        pair_reads = _pair_reads
} // end of wf_partition_fastq

process PartitionFastQ {
    // label 'PartitionFastQ'
    label 'cpus_32'

    tag {idPatient + "-" + idRun}

    // publishDir "${params.outdir}/Reports/${idSample}/FastQC/${idSample}_${idRun}", 
    // mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, idRun, file("${idSample}_${idRun}_R1.fastq.gz"), 
        file("${idSample}_${idRun}_R2.fastq.gz")

    output:
        tuple idPatient, idSample, idRun, file("r1_split_*"), file("r2_split_*") 
        // val (tuple_list)

    when: params.split_fastq && step == 'mapping'
    
    script:
    """
    init.sh
    partition.sh \
                in=${idSample}_${idRun}_R1.fastq.gz \
                in2=${idSample}_${idRun}_R2.fastq.gz  \
                out=r1_split_%.fastq.gz \
                out2=r2_split_%.fastq.gz \
                ways=${params.split_fastq}
    """
}