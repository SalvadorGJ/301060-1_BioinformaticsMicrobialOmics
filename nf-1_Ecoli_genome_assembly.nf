#!/usr/bin/env nextflow

// nf-1_Ecoli_genome_assembly.nf

/*
 * Defines the pipeline inputs parameters
 */

// Pipeline root directory
params.rootDir                  = "/lisc/user/gonzalez/MicroBioInfo/Homework_1/"

// Input data
params.rawDataDir               = "${params.rootDir}/data/raw_data/"

// Processed data
params.processedDataDir         = "${params.rootDir}/data/processed_data/"

// Results directory
params.resultsDir               = "${params.rootDir}/results/"

//* Input data
params.illumminaReads           = "SRR15361790"
params.nanoporeReads            = "SRR17645346"


/**********************************************************************************/
/******************************* W O R K F L O W S ********************************/
/**********************************************************************************/

workflow {

    // Download the FASTQ reads by the SRA IDs
    Channel
        .of(["${params.illumminaReads}", "${params.nanoporeReads}"])
        .set { sra_ids_ch }
    raw_reads_ch = download_data(sra_ids_ch)

    // Compute raw FASTQC statistics workflow
    raw_fastqc(raw_reads_ch.illumina)

    // Trim raw reads
    trimmed_reads_ch = trimmomatic(raw_reads_ch.illumina)

    // Compute processed FASTQC statistics workflow
    processed_fastqc(trimmed_reads_ch)

    // Illumina-only assembly
    assembly_illumina_only(trimmed_reads_ch) | quast_illumina_only

    // Hybrid assembly
    assembly_hybrid(trimmed_reads_ch, raw_reads_ch.nanopore) | quast_hybrid

    // Nanopore-only assembly
    assembly_nanopore_only(raw_reads_ch.nanopore) | quast_nanopore_only

}

workflow raw_fastqc {

    take:
    raw_fastqc_input_ch

    main:

    // Run fastqc report
    Channel
        .from('raw')
        .set { sample_type_ch }
    fastqc(raw_fastqc_input_ch.combine(sample_type_ch))

}

workflow processed_fastqc {

    take:
    processed_fastqc_input_ch

    main:

    // Run fastqc report
    Channel
        .from('processed')
        .set { sample_type_ch }
    fastqc(processed_fastqc_input_ch.combine(sample_type_ch))

}

/**********************************************************************************/
/******************************* P R O C E S S E S ********************************/
/**********************************************************************************/

/**********************************************************************************/
/****************************** Download FASTQ reads ******************************/
process download_data {
    executor = 'slurm'
    cpus = 6
    time = { 10.min * task.attempt }
    memory = { 2.GB * task.attempt }
    errorStrategy = { 'retry' }
    maxRetries = 1
    module 'sratoolkit/3.1.1'

    input:
    tuple val(illumina), val(nanopore)

    output:
    tuple path("input/${illumina}*_1.fastq"), path("input/${illumina}*_2.fastq"), emit: illumina
    path("input/${nanopore}*.fastq"), emit: nanopore

    script:
    """
    # Download Illumina reads
    fasterq-dump ${illumina} -O input -t /tmp

    # Download Nanopore reads
    fasterq-dump ${nanopore} -O input -t /tmp -e ${task.cpus} 
    """
}

/**********************************************************************************/
/******************************* FASTQC statistics ********************************/
process fastqc {
    executor = 'slurm'
    cpus = 2
    time = { 5.min * task.attempt }
    memory = { 2.GB * task.attempt }
    errorStrategy = { 'retry' }
    maxRetries = 1
    module 'fastqc/0.12.1'

    // FASTQC results
    publishDir "${params.resultsDir}/fastqc/${sample_type}/", mode: 'copy', overwrite: true

    input:
    tuple path(read1), path(read2), val(sample_type)

    output:
    path("*html")

    script:
    """
    # Quality report of reads
    fastqc -t ${task.cpus} ${read1} ${read2} -o .
    """
}

/**********************************************************************************/
/********************************* Trim raw reads *********************************/
process trimmomatic {
    executor = 'slurm'
    cpus = 4
    time = { 10.min * task.attempt }
    memory = { 2.GB * task.attempt }
    errorStrategy = { 'retry' }
    maxRetries = 1
    module 'trimmomatic/0.39'

    input:
    tuple path(read1), path(read2)

    output:
    tuple path("*_1.trimmomatic_paired.fq"), path("*_2.trimmomatic_paired.fq")

    script:
    """
    # Trim raw reads and remove adapters
    java -jar /lisc/app/trimmomatic/0.39/trimmomatic-0.39.jar PE \
        ${read1} ${read2} \
        ecoli_s788309_1.trimmomatic_paired.fq ecoli_s788309_1.trimmomatic_unpaired.fq \
        ecoli_s788309_2.trimmomatic_paired.fq ecoli_s788309_2.trimmomatic_unpaired.fq \
        ILLUMINACLIP:/lisc/app/trimmomatic/0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:30:20 MINLEN:90 -threads ${task.cpus}
    """
}

/**********************************************************************************/
/***************************** Illumina-only assembly *****************************/
process assembly_illumina_only {
    executor = 'slurm'
    cpus = 4
    time = { 30.min * task.attempt }
    memory = { 5.GB * task.attempt }
    errorStrategy = { 'retry' }
    maxRetries = 1
    module 'spades/4.0.0'

    // Genome assembly
    publishDir "${params.processedDataDir}/", pattern: "spades_filtered_illumina-only/scaffolds.fasta", mode: 'copy', saveAs: { filename -> "spades_filtered_illumina-only_scaffolds.fasta" }, overwrite: true

    input:
    tuple path(read1), path(read2)

    output:
    tuple path(read1), path(read2), path("spades_filtered_illumina-only/scaffolds.fasta")

    script:
    """
    # Illumina-only assembly
    spades.py -1 ${read1} -2 ${read2} -o spades_filtered_illumina-only --isolate -t ${task.cpus}
    """
}

/**********************************************************************************/
/*********************** Illumina-only assembly statistics ************************/
process quast_illumina_only {
    executor = 'slurm'
    cpus = 4
    time = { 15.min * task.attempt }
    memory = { 5.GB * task.attempt }
    errorStrategy = { 'retry' }
    maxRetries = 1
    module 'conda/'

    // QUAST report
    publishDir "${params.resultsDir}/quast/", mode: 'copy', saveAs: { filename -> "spades_filtered_illumina-only.report.txt" }, overwrite: true

    input:
    tuple path(read1), path(read2), path(assembly)

    output:
    path("quast_results/latest/report.txt")

    script:
    """
    # Statistics and evaluation of Illumina-only assembly
    conda activate quast
    quast.py -1 ${read1} -2 ${read2} -t ${task.cpus} ${assembly}
    """
}

/**********************************************************************************/
/******************************** Hybrid assembly *********************************/
process assembly_hybrid {
    executor = 'slurm'
    cpus = 4
    time = { 30.min * task.attempt }
    memory = { 5.GB * task.attempt }
    errorStrategy = { 'retry' }
    maxRetries = 1
    module 'spades/4.0.0'

    // Genome assembly
    publishDir "${params.processedDataDir}/", pattern: "spades_filtered_hybrid/scaffolds.fasta", mode: 'copy', saveAs: { filename -> "spades_filtered_hybrid_scaffolds.fasta" }, overwrite: true

    input:
    tuple path(read1), path(read2)
    path(nanopore)

    output:
    tuple path(read1), path(read2), path(nanopore), path("spades_filtered_hybrid/scaffolds.fasta")

    script:
    """
    # Hybrid assembly
    spades.py -1 ${read1} -2 ${read2} --nanopore ${nanopore} -o spades_filtered_hybrid --isolate -t ${task.cpus}
    """
}

/**********************************************************************************/
/*************************** Hybrid assembly statistics ***************************/
process quast_hybrid {
    executor = 'slurm'
    cpus = 4
    time = { 30.min * task.attempt }
    memory = { 5.GB * task.attempt }
    errorStrategy = { 'retry' }
    maxRetries = 1
    module 'conda/'

    // QUAST report
    publishDir "${params.resultsDir}/quast/", mode: 'copy', saveAs: { filename -> "spades_filtered_hybrid.report.txt" }, overwrite: true

    input:
    tuple path(read1), path(read2), path(nanopore), path(assembly)

    output:
    path("quast_results/latest/report.txt")

    script:
    """
    # Statistics and evaluation of hybrid assembly
    conda activate quast
    quast.py -1 ${read1} -2 ${read2} --nanopore ${nanopore} -t ${task.cpus} ${assembly}
    """
}

/**********************************************************************************/
/***************************** Nanopore-only assembly *****************************/
process assembly_nanopore_only {
    executor = 'slurm'
    cpus = 30
    time = { 15.min * task.attempt }
    memory = { 15.GB * task.attempt }
    errorStrategy = { 'retry' }
    maxRetries = 1
    module 'conda/'

    // Genome assembly
    publishDir "${params.processedDataDir}/", pattern: "flye_nanopore-only/assembly.fasta", mode: 'copy', saveAs: { filename -> "flye_nanopore-only_scaffolds.fasta" }, overwrite: true

    input:
    path(nanopore)

    output:
    tuple path(nanopore), path("flye_nanopore-only/assembly.fasta")

    script:
    """
    # Nanopore-only assembly
    conda activate flye
    flye --nano-raw ${nanopore} -o flye_nanopore-only -t ${task.cpus} -i 5 -g 5m
    """
}

/**********************************************************************************/
/*********************** Nanopore-only assembly statistics ************************/
process quast_nanopore_only {
    executor = 'slurm'
    cpus = 4
    time = { 15.min * task.attempt }
    memory = { 5.GB * task.attempt }
    errorStrategy = { 'retry' }
    maxRetries = 1
    module 'conda/'

    // QUAST report
    publishDir "${params.resultsDir}/quast/", mode: 'copy', saveAs: { filename -> "flye_nanopore-only.report.txt" }, overwrite: true

    input:
    tuple path(nanopore), path(assembly)

    output:
    path("quast_results/latest/report.txt")

    script:
    """
    # Statistics and evaluation of hybrid assembly
    conda activate quast
    quast.py --nanopore ${nanopore} -t ${task.cpus} ${assembly}
    """
}
