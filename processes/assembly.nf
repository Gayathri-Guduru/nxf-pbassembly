process ASSEMBLY {
    conda "./environment.yaml"
    publishDir "$params.publish_dir/assembly", mode: 'copy'

    input:
        tuple val(sample), path(reads)

    output: 
    tuple val(sample), path("*.r_utg.gfa")       , emit: raw_unitigs
    tuple val(sample), path("*.ec.bin")          , emit: corrected_reads
    tuple val(sample), path("*.ovlp.source.bin") , emit: source_overlaps
    tuple val(sample), path("*.ovlp.reverse.bin"), emit: reverse_overlaps
    tuple val(sample), path("*.p_utg.gfa")       , emit: processed_unitigs, optional: true
    tuple val(sample), path("*.asm.p_ctg.gfa")   , emit: primary_contigs  , optional: true
    tuple val(sample), path("*.asm.a_ctg.gfa")   , emit: alternate_contigs, optional: true
    tuple val(sample), path("*.hap1.p_ctg.gfa")  , emit: paternal_contigs , optional: true
    tuple val(sample), path("*.hap2.p_ctg.gfa")  , emit: maternal_contigs , optional: true

    script:
    """
    hifiasm -o ${sample}.asm -t $params.cpus ${reads}

    """
}