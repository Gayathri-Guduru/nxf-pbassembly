process GFA_TO_FASTA {
    tag "${sample}"
    publishDir "${params.publish_dir}/gfa_2_fasta", mode: 'copy'

    input:
        tuple val(sample), path(gfa_file), val(assembler)

    output:
        tuple val(sample), path("${sample}.fasta"), val(assembler)

    script:
    """
    # Find gfatools in a more robust way
    if command -v gfatools >/dev/null 2>&1; then
        GFATOOLS="gfatools"
    else
        # Try to find gfatools anywhere in the container
        GFATOOLS_PATH=\$(find / -name gfatools -type f -executable 2>/dev/null | head -n 1 || echo "")
        
        if [ -n "\$GFATOOLS_PATH" ]; then
            GFATOOLS="\$GFATOOLS_PATH"
        else
            echo "ERROR: Could not find gfatools executable in container"
            exit 1
        fi
    fi
    
    echo "Using gfatools at: \$GFATOOLS"
    
    # Convert GFA to FASTA
    \$GFATOOLS gfa2fa ${gfa_file} > ${sample}.fasta

    # Print some debug information
    echo "Converted GFA to FASTA for sample: ${sample}"
    echo "Assembler: ${assembler}"
    echo "Output file size: \$(ls -lh ${sample}.fasta | awk '{print \$5}')"
    """
}
