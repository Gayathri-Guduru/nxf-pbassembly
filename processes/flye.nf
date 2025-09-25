process ASSEMBLY_FLYE {
    publishDir "${params.publish_dir}/assembly/flye/${sample}", mode: 'copy'
    tag "${sample}"
    errorStrategy 'ignore'  // Continue workflow even if this process fails

    input:
        tuple val(sample), path(reads), val(genome_size)
    output:
        tuple val(sample), path("assembly.fasta"), emit: primary_contigs, optional: true
        tuple val(sample), path("assembly_graph.gfa"), emit: assembly_graph, optional: true
        tuple val(sample), path("assembly_info.txt"), emit: assembly_info, optional: true
        tuple val(sample), path("params.json"), emit: params_json, optional: true
        tuple val(sample), path("flye.log"), emit: flye_log, optional: true
        tuple val(sample), path("${sample}.assembly_status"), emit: assembly_status

    script:
    def genome_size_param = genome_size ? "--genome-size ${genome_size}m" : ""
    """
    #!/bin/bash

    # Create a status file to track assembly results
    echo "RUNNING" > ${sample}.assembly_status

    # First try with standard parameters (with or without genome size)
    flye --pacbio-hifi ${reads} ${genome_size_param} --out-dir ./flye_tmp --threads ${params.cpus}
    RESULT=\$?

    # If standard assembly failed for any reason, try meta mode
    if [ \$RESULT -ne 0 ]; then
        echo "Standard assembly failed with exit code \$RESULT, trying with --meta option"
        # Save the original log before removing directory
        if [ -f "./flye_tmp/flye.log" ]; then
            mkdir -p logs
            cp ./flye_tmp/flye.log ./logs/standard_failed.log
        fi

        # Remove any partial output from the failed run
        rm -rf ./flye_tmp

        # Try again with --meta option (keeping genome_size if specified)
        flye --pacbio-hifi ${reads} ${genome_size_param} --meta --out-dir ./flye_tmp --threads ${params.cpus}
        META_RESULT=\$?

        if [ \$META_RESULT -eq 0 ] && [ -f "./flye_tmp/assembly.fasta" ]; then
            echo "Assembly completed successfully using --meta option as fallback" >> ./flye_tmp/flye.log
            echo "SUCCESS" > ${sample}.assembly_status
            
            # Copy results since the assembly was successful
            cp ./flye_tmp/assembly.fasta ./
            cp ./flye_tmp/assembly_graph.gfa ./
            cp ./flye_tmp/assembly_info.txt ./
            cp ./flye_tmp/params.json ./
            cp ./flye_tmp/flye.log ./
        else
            echo "Both standard and meta assembly approaches failed for ${sample}"
            echo "FAILED" > ${sample}.assembly_status
            # Exit with 0 to not trigger a failure in the workflow
            exit 0
        fi
    else
        if [ -f "./flye_tmp/assembly.fasta" ]; then
            echo "SUCCESS" > ${sample}.assembly_status
            
            # Copy results since the assembly was successful
            cp ./flye_tmp/assembly.fasta ./
            cp ./flye_tmp/assembly_graph.gfa ./
            cp ./flye_tmp/assembly_info.txt ./
            cp ./flye_tmp/params.json ./
            cp ./flye_tmp/flye.log ./
        else
            echo "Flye process completed with code 0 but did not produce assembly.fasta"
            echo "FAILED" > ${sample}.assembly_status
        fi
    fi
    """
}
