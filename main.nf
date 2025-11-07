#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {
    // Create channel from samplesheet: [ sample_id, fq1, fq2, index_file ]
    ch_samplesheet = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header:true, sep:',')
        .map { row -> tuple(row.sample_id, file(row.fastq_1), file(row.fastq_2), file(row.index)) }

    // 1. Pre-trimming QC
    QC_PRE(ch_samplesheet)

    // 2. Adapter and Quality Trimming
    TRIMMING(ch_samplesheet)

    // 3. Diversity Trimming
    DIVERSITY_TRIMMING(TRIMMING.out.trimmed)

    // 4. Post-trimming QC
    QC_POST(DIVERSITY_TRIMMING.out.div_trimmed)

    // 5. Alignment
    ALIGNMENT(DIVERSITY_TRIMMING.out.div_trimmed)

    // 6. Deduplication (Corrected 3-step process)
    DEDUP_PREP(ALIGNMENT.out.aligned_bam)
    STRIP_SAM(DEDUP_PREP.out.sorted_sam)
    DEDUP(STRIP_SAM.out.stripped_sam)


    // 7. Methylation Extraction
    METHYLATION_EXTRACTION(DEDUP.out.dedup_bam)


    // 8. FINAL MULTIQC STEP
    qc_pre_zips = QC_PRE.out.zip_files.map { id, zip1, zip2 -> [zip1, zip2] }
    qc_post_zips = QC_POST.out.zip_files.map { id, zip1, zip2 -> [zip1, zip2] }
    // Mix all reports and logs from all previous steps into one channel
    ch_multiqc_reports = Channel.empty()
        .mix(qc_pre_zips)
        .mix(qc_post_zips)
        .mix(TRIMMING.out.report)
        .mix(ALIGNMENT.out.report)
        .mix(METHYLATION_EXTRACTION.out.report)
        .flatten()
        .collect() // Wait for all processes to finish and collect all files

    ch_multiqc_config = file('multiqc_config.yaml')

    MULTIQC(ch_multiqc_reports, ch_multiqc_config)
    //MULTIQC(ch_multiqc_reports)
}

// --- PROCESS DEFINITIONS ---

process QC_PRE {
    tag "FastQC (before): ${sample_id}"
    label 'fastqc'
    publishDir "${params.outdir}/qc_pre", mode: 'copy', pattern: "*.{zip,html}"

    input:
    tuple val(sample_id), path(fastq_1), path(fastq_2), path(index)

    output:
    tuple val(sample_id), path("${sample_id}_R1.fastqc.html"), path("${sample_id}_R2.fastqc.html"), emit: html_files
    tuple val(sample_id), path("${sample_id}_R1_fastqc.zip"), path("${sample_id}_R2_fastqc.zip"), emit: zip_files

    script:
    """
    fastqc --outdir . --threads ${task.cpus} ${fastq_1} ${fastq_2}

    mv \$(basename ${fastq_1} .fastq.gz)_fastqc.html "${sample_id}_R1.fastqc.html"
    mv \$(basename ${fastq_2} .fastq.gz)_fastqc.html "${sample_id}_R2.fastqc.html"
    mv \$(basename ${fastq_1} .fastq.gz)_fastqc.zip  "${sample_id}_R1_fastqc.zip"
    mv \$(basename ${fastq_2} .fastq.gz)_fastqc.zip  "${sample_id}_R2_fastqc.zip"
    """
}

process TRIMMING {
    tag "TrimGalore: ${sample_id}"
    label 'trimming'
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy', pattern: '{*.trim.fq.gz,*.trim_report.txt}'

    input:
    tuple val(sample_id), path(fastq_1), path(fastq_2), path(index)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trim.fq.gz"), path("${sample_id}_R2.trim.fq.gz"), path(index), emit: trimmed
    path("*trimming_report.txt"), emit: report

    script:
    """
    trim_galore --paired --cores ${task.cpus} --output_dir . ${fastq_1} ${fastq_2}

    # Find the output files
    R1_TRIMMED=\$(find . -name "*_val_1.fq.gz")
    R2_TRIMMED=\$(find . -name "*_val_2.fq.gz")

    # Rename them to the clean, standardized format 
    mv "\$R1_TRIMMED" "${sample_id}_R1.trim.fq.gz"
    mv "\$R2_TRIMMED" "${sample_id}_R2.trim.fq.gz"
    """
}

process DIVERSITY_TRIMMING {
    label 'diversity_trimming'
    tag "Diversity Trim: ${sample_id}"
    publishDir "${params.outdir}/divtrimmed", mode: 'copy', pattern: '{*.divtrim.fq.gz}'
    stageInMode = 'copy'

    input:
    tuple val(sample_id), path(trimmed_1), path(trimmed_2), path(index)

    output:
    tuple val(sample_id), path("${sample_id}_R1.divtrim.fq.gz"), path("${sample_id}_R2.divtrim.fq.gz"), path(index), emit: div_trimmed
    
    script:
    """
    python ${params.scripts_dir}/trim_div.py -1 ${trimmed_1} -2 ${trimmed_2}
    R1_DIVTRIMMED=\$(find . -name "*_R1*trimmed.fq.gz")
    R2_DIVTRIMMED=\$(find . -name "*_R2*trimmed.fq.gz")

    mv "\$R1_DIVTRIMMED" "${sample_id}_R1.divtrim.fq.gz"
    mv "\$R2_DIVTRIMMED" "${sample_id}_R2.divtrim.fq.gz"

    ln -s ${index} index
    """
}

process QC_POST {
    tag "FastQC (after): ${sample_id}"
    label 'fastqc'
    publishDir "${params.outdir}/qc_post", mode: 'copy', pattern: "*.{zip,html}"

    input:
    tuple val(sample_id), path(div_trimmed_1), path(div_trimmed_2), path(index)

    output:
    tuple val(sample_id), path("${sample_id}_R1.fastqc_post.html"), path("${sample_id}_R2.fastqc_post.html"), emit: html_files
    tuple val(sample_id), path("${sample_id}_R1_post_fastqc.zip"), path("${sample_id}_R2_post_fastqc.zip"), emit: zip_files

    script:
    """
    fastqc --outdir . --threads ${task.cpus} ${div_trimmed_1} ${div_trimmed_2}

    mv \$(basename ${div_trimmed_1} .fq.gz)_fastqc.html "${sample_id}_R1.fastqc_post.html"
    mv \$(basename ${div_trimmed_2} .fq.gz)_fastqc.html "${sample_id}_R2.fastqc_post.html"
    mv \$(basename ${div_trimmed_1} .fq.gz)_fastqc.zip  "${sample_id}_R1_post_fastqc.zip"
    mv \$(basename ${div_trimmed_2} .fq.gz)_fastqc.zip  "${sample_id}_R2_post_fastqc.zip"
    """
}


process ALIGNMENT {
    tag "Bismark Align: ${sample_id}"
    label 'bismark'
    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: '{*.align.bam,*.align_report.txt}'

    input:
    tuple val(sample_id), path(div_trimmed_1), path(div_trimmed_2), path(index)

    output:
    tuple val(sample_id), path("${sample_id}.align.bam"), path(index), emit: aligned_bam
    path("*_PE_report.txt"), emit: report

    script:
    """
    bismark --bowtie2 \\
        --genome_folder "${params.genome_dir}" \\
        --parallel ${task.cpus} \\
        -1 ${div_trimmed_1} \\
        -2 ${div_trimmed_2} \\
        -o .

    DEFAULT_BAM=\$(find . -name "*_pe.bam")
    DEFAULT_REPORT=\$(find . -name "*PE_report.txt")

    mv "\$DEFAULT_BAM" "${sample_id}.align.bam"
    mv "\$DEFAULT_REPORT" "${sample_id}_PE_report.txt"
    """
}

process DEDUP_PREP {
    label 'samtools'
    tag "Prep: ${sample_id}"
    stageInMode = 'copy'
    
    input:
    tuple val(sample_id), path(aligned_bam), path(index)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.sam"), path(index), emit: sorted_sam

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -T tmp_sort ${aligned_bam}
    samtools view -h -o ${sample_id}.sorted.sam ${sample_id}.sorted.bam
    """
}

process STRIP_SAM {
    label 'dedup' // Assign to the python container
    tag "Strip SAM: ${sample_id}"
    publishDir "${params.outdir}/dedup_prepped", mode: 'copy', pattern: '*.stripped.sam'
    stageInMode = 'copy'

    input:
    tuple val(sample_id), path(sorted_sam), path(index)

    output:
    tuple val(sample_id), path("${sample_id}.stripped.sam"), path(index), emit: stripped_sam

    script:
    """
    bash ${params.scripts_dir}/strip_bs.sh ${sorted_sam}
    
    # Check if stripping produced a file, exit if it failed
    if [ ! -s "${sorted_sam}_stripped.sam" ]; then
        echo "ERROR: strip_bs.sh failed to produce output for ${sorted_sam}." >&2
        exit 1
    fi

    mv "${sorted_sam}_stripped.sam" "${sample_id}.stripped.sam"
    """
}


process DEDUP {
    label 'no_container' 
    tag "Dedup: ${sample_id}"
    publishDir "${params.outdir}/deduplicated_sam", mode: 'copy', pattern: '*.rmdup.sam'
    stageInMode = 'copy'
    
    input:
    tuple val(sample_id), path(stripped_sam), path(index)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), emit: dedup_bam

    script:
    """
    TEMP_DIR="${sample_id}_temp"
    mkdir -p "\${TEMP_DIR}"
    module load gcccore/11.3.0 python/2.7.18

    # Run nudup.py on the correctly prepared stripped_sam file
    python2.7 ${params.scripts_dir}/nudup.py -2 -f "${index}" -o "${sample_id}" -T "\${TEMP_DIR}" --rmdup-only "${stripped_sam}"

    find . -name "${sample_id}*.bam" | xargs -I {} mv {} "${sample_id}.dedup.bam"

    rm -rf "\${TEMP_DIR}"
    """
}



process METHYLATION_EXTRACTION {
    label 'bismark'
    tag "Methyl Extract: ${sample_id}"
    publishDir "${params.outdir}/methyl_extract", mode: 'copy'

    input:
    tuple val(sample_id), path(dedup_bam)

    output:
    tuple val(sample_id), path("${sample_id}.bismark.cov.gz"), emit: cov
    path("*_splitting_report.txt"), emit: report
    
    script:
    """
    samtools sort -n -@ ${task.cpus} -o "${sample_id}.sorted.bam" ${dedup_bam}

    bismark_methylation_extractor \\
        --paired-end \\
        --comprehensive \\
        --bedGraph \\
        --counts \\
        --gzip \\
        --parallel ${task.cpus} \\
        --output . \\
        "${sample_id}.sorted.bam"

    # Find the output coverage file and rename it to the desired clean format
    COV_FILE=\$(find . -name "*.cov.gz")
    mv "\$COV_FILE" "${sample_id}.bismark.cov.gz"
    """
}

process MULTIQC {
    tag "MultiQC (Final Report)"
    label 'multiqc'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    val(reports) // Takes the collected list of all report files
    path(config)

    output:
    path("multiqc_report.html")
    path("multiqc_report_data", type: 'dir')

    script:
    """
    printf "%s\n" ${reports.join(' ')} > custom_file_list.txt

    multiqc --file-list custom_file_list.txt -c ${config} -n multiqc_report.html
    """
}

