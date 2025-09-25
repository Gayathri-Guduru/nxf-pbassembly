process BAM2FASTQ {
    //conda "./environment.yaml"
    //container "quay.io/staphb/samtools"
    publishDir "$params.publish_dir/fastq", mode: 'copy'

    input:
        tuple val(sample), path(input), path(genome), path(gff)

    output: 
        tuple val(sample), path("*.fastq"), emit: reads

    script:
    """
    samtools bam2fq ${input} > ${sample}.fastq
    """
}
