process METAEUK {
    container = 'gayathriguduru/metaeuk:latest'
    // Publish output inside a folder named after the sample
    publishDir "${params.publish_dir}/metaeuk/${assembler}/${sample}", mode: 'copy'
    tag "${sample}_${assembler}"
    
    // Add errorStrategy to continue workflow even if this process fails
    errorStrategy 'ignore'

    input:
        // Expect tuple: (sample, assembly FASTA, assembler, protein_reference)
        tuple val(sample), path(fasta), val(assembler), val(protein_reference)

    output:
        // Expect a directory named "metaeuk_output" containing the output files
        path "metaeuk_output", emit: results

    script:
        """
        echo "============================================================"
        echo "MetaEuk PROCESS STARTING FOR: ${sample} (${assembler})"
        echo "============================================================"
        echo "Sample ID: ${sample}"
        echo "Assembler: ${assembler}"
        echo "Assembly file: ${fasta} (size: \$(ls -lh ${fasta} | awk '{print \$5}'))"
        echo "Protein reference: ${protein_reference}"
        echo "Current directory: \$(pwd)"
        echo "Directory contents:"
        ls -la
        echo "============================================================"

        # Ensure we are in the current working directory
        cd \$(pwd)

        # 1) Check if protein_reference is "NA" or empty, then create empty output dir and exit
        if [ "${protein_reference}" = "NA" ] || [ "${protein_reference}" = "" ]; then
            echo "No valid protein reference provided. Skipping MetaEuk for sample: ${sample}"
            mkdir -p metaeuk_output
            touch metaeuk_output/.skipped
            echo "Sample: ${sample} - MetaEuk skipped (bacterial sample)" > metaeuk_output/status.txt
            exit 0
        fi

        echo "PROCEEDING WITH METAEUK FOR FUNGAL SAMPLE: ${sample}"

        # 2) If protein_reference is a URL, download it; otherwise use the provided file.
        if [[ "${protein_reference}" =~ ^https?:// || "${protein_reference}" =~ ^ftp:// ]]; then
            echo "Downloading protein reference from URL: ${protein_reference}"
            wget -O reference.faa.gz "${protein_reference}" || {
                echo "ERROR: Failed to download ${protein_reference}"
                mkdir -p metaeuk_output
                echo "Sample: ${sample} - MetaEuk failed (download error)" > metaeuk_output/status.txt
                exit 1
            }
            PROTEIN_REF=reference.faa.gz
        else
            echo "Using provided protein reference file: ${protein_reference}"
            PROTEIN_REF=${protein_reference}
            
            # Verify the protein reference file exists
            if [ ! -f "\${PROTEIN_REF}" ]; then
                echo "ERROR: Protein reference file \${PROTEIN_REF} does not exist"
                mkdir -p metaeuk_output
                echo "Sample: ${sample} - MetaEuk failed (protein reference not found)" > metaeuk_output/status.txt
                exit 1
            fi
        fi

        # Use fixed names for the databases (you can adjust these as needed)
        PROTEIN_DB="Sac_${sample}_referenceDB"
        ASSEMBLY_DB="Sac_${sample}_${assembler}DB"

        echo "Creating protein database"
        /usr/local/bin/metaeuk_avx2 createdb \$PROTEIN_REF \$PROTEIN_DB || {
            echo "ERROR: Failed to create protein database"
            mkdir -p metaeuk_output
            echo "Sample: ${sample} - MetaEuk failed (protein database creation error)" > metaeuk_output/status.txt
            exit 1
        }

        echo "Creating assembly database"
        /usr/local/bin/metaeuk_avx2 createdb ${fasta} \$ASSEMBLY_DB || {
            echo "ERROR: Failed to create assembly database"
            mkdir -p metaeuk_output
            echo "Sample: ${sample} - MetaEuk failed (assembly database creation error)" > metaeuk_output/status.txt
            exit 1
        }

        echo "Running MetaEuk easy-predict"
        # Force the output files to be prefixed with 'metaeuk_output'
        /usr/local/bin/metaeuk_avx2 easy-predict \$ASSEMBLY_DB \$PROTEIN_DB metaeuk_output tmp || {
            echo "ERROR: MetaEuk easy-predict failed"
            mkdir -p metaeuk_output
            echo "Sample: ${sample} - MetaEuk failed (easy-predict error)" > metaeuk_output/status.txt
            exit 1
        }

        echo "MetaEuk annotation complete for sample: ${sample}"

        # Ensure all output files are organized in the 'metaeuk_output' folder
        mkdir -p metaeuk_output
        mv metaeuk_output.fas metaeuk_output.codon.fas metaeuk_output.headersMap.tsv metaeuk_output.gff metaeuk_output/ 2>/dev/null || true
        echo "Sample: ${sample} - MetaEuk completed successfully" > metaeuk_output/status.txt
        echo "============================================================"
        echo "MetaEuk PROCESS COMPLETED FOR: ${sample} (${assembler})"
        echo "============================================================"
        """
}
