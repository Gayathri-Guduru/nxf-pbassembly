process ASSEMBLY_HIFIASM {
    publishDir "${params.publish_dir}/assembly/hifiasm/${sample}", mode: 'copy'
    tag "${sample}"
    errorStrategy 'ignore'  // Continue workflow even if this process fails

    input:
        tuple val(sample), path(reads)
    output:
        tuple val(sample), path("*.asm.bp.p_ctg.gfa"), emit: primary_contigs, optional: true
        tuple val(sample), path("*.asm.bp.a_ctg.gfa"), emit: alternate_contigs, optional: true
        tuple val(sample), path("*.bp.r_utg.gfa"), emit: raw_unitigs, optional: true
        tuple val(sample), path("*.ec.bin"), emit: corrected_reads, optional: true
        tuple val(sample), path("*.ovlp.source.bin"), emit: source_overlaps, optional: true
        tuple val(sample), path("*.ovlp.reverse.bin"), emit: reverse_overlaps, optional: true
        tuple val(sample), path("*.bp.p_utg.gfa"), emit: processed_unitigs, optional: true
        tuple val(sample), path("*.asm.bp.hap1.p_ctg.gfa"), emit: paternal_contigs, optional: true
        tuple val(sample), path("*.asm.bp.hap2.p_ctg.gfa"), emit: maternal_contigs, optional: true
        tuple val(sample), path("${sample}.assembly_status"), emit: assembly_status

    script:
    """
    #!/bin/bash

    # Create a status file to track assembly results
    echo "RUNNING" > ${sample}.assembly_status

    # Run hifiasm
    hifiasm -o ${sample}.asm -t ${params.cpus} ${reads}
    RESULT=\$?

    if [ \$RESULT -eq 0 ] && [ -f "${sample}.asm.bp.p_ctg.gfa" ]; then
        echo "SUCCESS" > ${sample}.assembly_status
        echo "HiFiasm completed successfully for ${sample}"
    else
        echo "HiFiasm failed to produce expected outputs for ${sample}"
        echo "FAILED" > ${sample}.assembly_status
        # Don't exit with error to allow the workflow to continue
    fi
    """
}
