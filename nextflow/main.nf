#!/usr/bin/env nextflow

// Enable DSL2 syntax 
nextflow.enable.dsl=2

// ==============================================================================
// Define Processes (Analysis Modules)
// ==============================================================================

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/FastQC", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{zip,html}", emit: qc_reports

    script:
    """
    fastqc -t ${task.cpus} ${reads}
    """
}

process CUTADAPT {
    tag "$sample_id"
    publishDir "${params.outdir}/Trimmed", mode: 'copy', pattern: '*.log'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_trimmed_R{1,2}.fastq"), emit: trimmed_reads
    path "*.log", emit: logs

    script:
    """
    cutadapt \\
        -a ${params.adapter_r1} \\
        -A ${params.adapter_r2} \\
        -o ${sample_id}_trimmed_R1.fastq \\
        -p ${sample_id}_trimmed_R2.fastq \\
        -m ${params.min_len} \\
        -q ${params.quality} \\
        -j ${task.cpus} \\
        ${reads[0]} ${reads[1]} > ${sample_id}_cutadapt.log
    """
}

process ALIGNMENT_HISAT2 {
    tag "$sample_id"
    publishDir "${params.outdir}/Aligned", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path index

    output:
    path "${sample_id}_sorted.bam", emit: bam
    path "${sample_id}_sorted.bam.bai", emit: bai

    script:
    """
    # Run Hisat2 alignment and pipe directly to samtools to save I/O
    hisat2 -p ${task.cpus} -x ${index} -1 ${reads[0]} -2 ${reads[1]} | \\
    samtools view -bS - | \\
    samtools sort -@ ${task.cpus} -o ${sample_id}_sorted.bam
    
    # Generate BAM index
    samtools index -@ ${task.cpus} ${sample_id}_sorted.bam
    """
}

process QUANTIFICATION {
    publishDir "${params.outdir}/Quantification", mode: 'copy'

    input:
    path bams
    path gtf

    output:
    path "counts.txt", emit: count_matrix
    path "counts.txt.summary"

    script:
    """
    featureCounts \\
        -T ${task.cpus} \\
        -p -B -C \\
        -a ${gtf} \\
        -o counts.txt \\
        ${bams}
    """
}

// ==============================================================================
// Main Workflow definition
// ==============================================================================

workflow {
    // 1. Generate input channels
    ch_reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
    ch_index = file(params.index)
    ch_gtf   = file(params.gtf)

    // 2. Quality Control
    FASTQC(ch_reads)

    // 3. Adapter Trimming
    CUTADAPT(ch_reads)

    // 4. Read Alignment (Takes trimmed reads from CUTADAPT)
    ALIGNMENT_HISAT2(CUTADAPT.out.trimmed_reads, ch_index)

    // 5. Gene Quantification (Collects all BAMs to run a single featureCounts instance)
    QUANTIFICATION(ALIGNMENT_HISAT2.out.bam.collect(), ch_gtf)
}
