process PROKKA {
    publishDir "${params.publish_dir}/prokka/${assembler}", mode: 'copy'
    tag "${sample}_${assembler}"

    input:
        tuple val(sample), path(fasta)
        val(assembler)

    output:
        path "${sample}/*"
        path("${sample}/${sample}.gff"), emit: gff, optional: true
        path("${sample}/${sample}.gbk"), emit: gbk, optional: true
        path("${sample}/${sample}.fna"), emit: fna, optional: true
        path("${sample}/${sample}.faa"), emit: faa, optional: true
        path("${sample}/${sample}.ffn"), emit: ffn, optional: true
        tuple val(sample), path("${sample}/${sample}.ffn"), emit: ffn_tuple, optional: true
        path("${sample}/${sample}.sqn"), emit: sqn, optional: true
        path("${sample}/${sample}.fsa"), emit: fsa, optional: true
        path("${sample}/${sample}.tbl"), emit: tbl, optional: true
        path("${sample}/${sample}.err"), emit: err, optional: true
        path("${sample}/${sample}.log"), emit: log, optional: true
        path("${sample}/${sample}.txt"), emit: txt, optional: true
        tuple val(sample), path("${sample}/${sample}.tsv"), emit: tsv, optional: true
        // Add new output that includes assembler information
        tuple val(sample), path("${sample}/${sample}.tsv"), path("${sample}/${sample}.ffn"), val(assembler), emit: for_16s, optional: true
        path "versions.yml", emit: versions

    script:
    def prefix = "${sample}"
    """
    # Find the prokka executable
    echo "Searching for prokka executable..."

    # Common locations for prokka in containers
    LOCATIONS=(
        "/opt/conda/bin/prokka"
        "/usr/local/bin/prokka"
        "/usr/bin/prokka"
        "/bin/prokka"
    )

    PROKKA_PATH=""
    for loc in "\${LOCATIONS[@]}"; do
        if [ -f "\$loc" ] && [ -x "\$loc" ]; then
            PROKKA_PATH="\$loc"
            break
        fi
    done

    # If not found in common locations, search the entire container
    if [ -z "\$PROKKA_PATH" ]; then
        FOUND_PATHS=\$(find / -name "prokka" -type f -executable 2>/dev/null || true)
        if [ -n "\$FOUND_PATHS" ]; then
            PROKKA_PATH=\$(echo "\$FOUND_PATHS" | head -n 1)
        fi
    fi

    if [ -z "\$PROKKA_PATH" ]; then
        # Try using which as a last resort
        PROKKA_PATH=\$(which prokka 2>/dev/null || true)
    fi

    # Check if we found Prokka
    if [ -z "\$PROKKA_PATH" ]; then
        echo "ERROR: Could not find prokka executable in container"
        echo "Contents of /usr/local/bin:"
        ls -la /usr/local/bin || true
        echo "Contents of /opt/conda/bin:"
        ls -la /opt/conda/bin || true
        exit 1
    fi

    echo "Using Prokka at: \$PROKKA_PATH"

    # Create versions file
    PROKKA_VERSION=\$(\$PROKKA_PATH --version 2>&1 | sed 's/^.*prokka //' || echo "unknown")
    cat <<END_VERSIONS > versions.yml
"${task.process}":
  prokka: \$PROKKA_VERSION
END_VERSIONS

    # Run Prokka directly on the FASTA file with the --compliant flag to work around issues
    # The --compliant flag makes Prokka more conservative and helps avoid issues with dependencies
    \$PROKKA_PATH --outdir ${prefix} --prefix ${prefix} --cpus ${params.cpus} --kingdom Bacteria --compliant ${fasta}
    """
}
