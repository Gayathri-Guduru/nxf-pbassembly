process EXTRACT_16S {
    publishDir path: "${params.publish_dir}/16S_sequences/${assembler}", mode: 'copy'
    tag "${sample}_${assembler}"

    input:
        tuple val(sample), path(tsv), path(ffn), val(assembler)

    output:
        tuple val(sample), path("${sample}/${sample}_16S.csv"), path("${sample}/${sample}_16S.fasta"), val(assembler), emit: results

    script:
    """
    # Create output directory
    mkdir -p ${sample}
    
    # Run the external Python script
    python3 ${projectDir}/bin/extract_16S.py ${tsv} ${ffn} ${sample} ${sample}
    
    # Print debug info
    echo "Extracted 16S sequences for ${sample} using ${assembler} assembler"
    echo "Saving to: ${params.publish_dir}/16S_sequences/${assembler}/${sample}/"
    """
}
