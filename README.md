# Layer Lab Cancer Analysis Workflow
 `PLEASE NOTE THAT THIS IS VERY MUCH A WORK IN PROGRESS AND THE CONFIGURATION AND THE CODE MIGHT NEED TWEAKING BEFORE RUNNING THE PIPELINE`

The pipeline is based on [nextflow](https://www.nextflow.io) and uses the following SV callers:
1. [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035535892-Somatic-copy-number-variant-discovery-CNVs-) (Mutect2 based workflow with small panel of normals)
2. [SavvyCNV](https://github.com/rdemolgen/SavvySuite) (based on off-read signals)
3. [CNVKit](https://cnvkit.readthedocs.io/en/stable/)

## Pipeline Architecture
Below are the major components of the pipeline
- Nextflow is the main glue that wraps the individual commands and shell scripts (bwa mem for example), and orchestrates the workflow
- Configuration files (under the `conf` directory) carries the information that nextflow uses to work out which [executor](https://www.nextflow.io/docs/latest/executor.html) to use to run the pipeline
- Helper utilities (under `lib`) are some helper functions written to aid the main pipeline
- Individual scripts wrapped inside `nextflow processes` to carry out the actual tasks such as the *alignment* or *marking duplicate reads*, or running a *structural variant caller*
- A Singularity container is our preferred way to run the pipeline. Currently we include the singularity configuration in the subdirectory *containers/dna_seq*

## How to run the pipeline
In the current configuration on our Colorado University Boulder Fiji Cluster, we are using a singularity container with all the required software installed in it. So if you plan to use this, you need to have `singularity` (current version we are using is 3.1.1, but should run with others too) on your Unix *PATH* .

For a run on your machine (or any other infrastructure such as an HPC, or AWS), you need to add corresponding configuration under the `conf`. The *configuration* defines a *profile* in Nextflow lingo and needs to be passed at the command-line when running the pipeline. See the `Makefile` in the top-level directory for an example run. At the command-prompt you will need to pass a *samples.tsv* (tab delimited) carrying at least three columns. The first column specifying the *sample name*, and the next two the paths to the *first* and the *second reads* respectively.
### Important parameters
A typical execution of the pipeline can be started as

    nextflow run main.nf \
		--sample-tsv=/Shares/layer_shared/projects/cancer_center_tiny/samples.tsv \
		--conditions-tsv /Shares/layer_shared/projects/cancer_center_tiny/conditions.tsv \
		--run-name='cancer_center_tiny' -profile fiji -resume

#### Passing parameters to the pipeline
Nextflow (and this pipeline) allows the following three ways to pass params to a pipeline run:
1. At the *commandline*, using syntex like *- -sample-tsv*, or *- -results-dir* (note the double dashes before param names)
2. Still at the commandline but putting params in a *yaml* file and by passing that file as *- -params my_params.yaml*. See below  an excerpt from a params.yaml file:

        sampleTsv: "/Shares/layer_shared/projects/cancer_center_tiny/samples.tsv"
        runName: "cancer_center_tiny"
        conditionsTsv: "/Shares/layer_shared/projects/cancer_center_tiny/conditions.tsv"
        sendEmail: False
        logFilesDur:    "logs"
        nfOptions: "-with-reports -with-trace -with-timeline -with-dag flowchart.png"
3. As part of the *profile config* file. See examples in the *conf* subdirectory.

#### *samples.tsv*
This file contains the sample names and the corresponding base calls (fastq files, the *read1* and *read2*)
The first column is the sample name, and the next two columns represents the first and the second read respectively. 

`See a sample tsv below`

    KU1919GC	KU1919GC/KU1919GC_TD180603264_HT3LKCCXY_L3_1.fq.gz	KU1919GC/KU1919GC_TD180603264_HT3LKCCXY_L3_2.fq.gz
    KU1919P	KU1919P/KU1919P_TD180603261_HT3LKCCXY_L2_1.fq.gz	KU1919P/KU1919P_TD180603261_HT3LKCCXY_L2_2.fq.gz
`The pipeline at the moment is not configured to run for GATK somatic2 workflow without a Panel of Normals.`

#### *conditions.tsv*
This file shows the *conditions* of the experiment (*parental* and the *conditions*, or *controls* and the *cases* etc). Each row in the file is treated as a record where the first column names the sample in the *parental* (or *control*), and the rest of the columns are treated as the corresponding *conditions* (or *cases*). The *parental* samples are used to create the so called *panel of normals* or *pon*, which is supposed to capture the sequencing artifacts or the *normal* samples. Please see [GATK Panel of normal](https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON-) for some more information regarding the *pons*

`See a conditions tsv below`

    KU1919P         KU1919C KU1919G KU1919GC
    RT112P          RT112C  RT112G  RT112GC

## Extras
#### Building singularity container
The `containers/dna_seq` subdirectory contains the required *environment.yaml* and the singularity recipe file to generate the singularity image that we use. Our current working version of singularity is 3.1.1. Depending upon your singularity configuration you will need something like the following to generate the image:

    sudo singularity build -F layer_lab_dna_seq_gatk_4.1.4.sif layer_lab_dna_seq.def

