process QUAST {
    // Maintain original directory structure
    publishDir path: { "${params.publish_dir}/quast/${assembler}" }, mode: 'copy'
    tag "${sample}_${assembler}"

    input:
        tuple val(sample), path(assembly), val(assembler)
        tuple val(sample), path(bam), path(reference), path(gff)

    output:
        path "${sample}", emit: quast_dir
        path '*.tsv', optional: true, emit: tsv

    script:
    """
    # Find Python without using 'which'
    if [ -f /usr/bin/python3 ]; then
        PYTHON="/usr/bin/python3"
    elif [ -f /usr/bin/python ]; then
        PYTHON="/usr/bin/python"
    elif [ -f /opt/conda/bin/python ]; then
        PYTHON="/opt/conda/bin/python"
    else
        # Default to python3 and hope it's in PATH
        PYTHON="python3"
    fi

    echo "Using Python: \$PYTHON"
    echo "Sample: ${sample}"
    echo "Assembly file: ${assembly}"
    echo "Assembler: ${assembler}"

    # Set PYTHONPATH with careful variable handling
    # Initialize PYTHONPATH to empty if not set
    export PYTHONPATH=""

    # Run QUAST with specific paths (no min-contig parameter)
    if [ -f /quast-5.3.0/quast.py ]; then
        echo "Using QUAST at: /quast-5.3.0/quast.py"
        export PYTHONPATH="/quast-5.3.0:\$PYTHONPATH"
        \$PYTHON /quast-5.3.0/quast.py -o ${sample} -r ${reference} -g ${gff} ${assembly}
    elif [ -f /quast-5.3.0/build/scripts-2.7/quast.py ]; then
        echo "Using QUAST at: /quast-5.3.0/build/scripts-2.7/quast.py"
        export PYTHONPATH="/quast-5.3.0:\$PYTHONPATH"
        \$PYTHON /quast-5.3.0/build/scripts-2.7/quast.py -o ${sample} -r ${reference} -g ${gff} ${assembly}
    else
        # Find quast.py
        echo "Searching for quast.py..."
        QUAST_PATHS=\$(find / -name "quast.py" -type f 2>/dev/null || echo "")

        if [ -n "\$QUAST_PATHS" ]; then
            QUAST_PATH=\$(echo "\$QUAST_PATHS" | head -n 1)
            echo "Using QUAST at: \$QUAST_PATH"
            QUAST_DIR=\$(dirname "\$QUAST_PATH")
            export PYTHONPATH="\$QUAST_DIR:\$PYTHONPATH"
            \$PYTHON "\$QUAST_PATH" -o ${sample} -r ${reference} -g ${gff} ${assembly}
        else
            echo "ERROR: QUAST not found in any expected location."
            exit 1
        fi
    fi
    """
}
