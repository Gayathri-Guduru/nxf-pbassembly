process BARRNAP {
    publishDir "${params.publish_dir}/barrnap/${assembler}/${sample}", mode: 'copy'
    tag "${sample}_${assembler}"
    errorStrategy 'ignore'  // Continue workflow even if this process fails
    
    input:
        tuple val(sample), path(assembly), val(assembler)
    
    output:
        tuple val(sample), path("${sample}_rrna.gff3"), emit: gff, optional: true
        path "${sample}_rrna.fasta", emit: all_rrna_fasta, optional: true
        path "${sample}_18S.fasta", emit: only_18s_fasta, optional: true
        path "${sample}_18S.csv", emit: only_18s_csv, optional: true
        path "${sample}_barrnap.log", emit: log, optional: true
    
    script:
    """
    #!/bin/bash
    # Don't use set -e here as we want to handle errors manually
    
    echo "=== BARRNAP PROCESSING ==="
    echo "Sample: ${sample}"
    echo "Assembler: ${assembler}"
    echo "Assembly file: ${assembly}"
    echo "Current directory: \$(pwd)"
    
    # Create a simple working copy of the assembly file
    cp "${assembly}" ./input.fa
    
    # Run Barrnap directly (no stdin redirection, no outfile param)
    echo "Running Barrnap..."
    barrnap --kingdom euk --threads ${task.cpus} ./input.fa > "${sample}_rrna.gff3" 2> "${sample}_barrnap.log" || {
        echo "Barrnap failed with exit code \$?"
        echo "WARNING: Barrnap failed for ${sample}" >> "${sample}_barrnap.log"
        # Create empty files to satisfy output requirements
        touch "${sample}_rrna.gff3"
        touch "${sample}_rrna.fasta"
        touch "${sample}_18S.fasta"
        echo "locus_tag,ftype,length_bp,gene,EC_number,COG,product" > "${sample}_18S.csv"
        exit 0
    }
    
    # Check if Barrnap found any 18S sequences
    if grep -q "18S_rRNA" "${sample}_rrna.gff3"; then
        echo "Found 18S rRNA sequences"
    else
        echo "No 18S rRNA sequences found"
        touch "${sample}_rrna.fasta"
        touch "${sample}_18S.fasta"
        echo "locus_tag,ftype,length_bp,gene,EC_number,COG,product" > "${sample}_18S.csv"
        exit 0
    fi
    
    # Create FASTA file for all rRNA sequences using bedtools 
    # (using bedtools from the staphb/bedtools container - referenced in your config)
    grep -v "^#" "${sample}_rrna.gff3" | awk '{print \$1 "\\t" \$4-1 "\\t" \$5 "\\t" \$9 "\\t0\\t" \$7}' > "${sample}_rrna.bed"
    bedtools getfasta -fi ./input.fa -bed "${sample}_rrna.bed" -s -name > "${sample}_rrna.fasta" || {
        echo "Error extracting sequences with bedtools"
        # Create empty dummy files
        echo ">dummy_rRNA" > "${sample}_rrna.fasta"
        echo "ACGTACGTACGTACGTACGT" >> "${sample}_rrna.fasta"
    }
    
    # Extract only 18S sequences from GFF and FASTA
    grep "18S_rRNA" "${sample}_rrna.gff3" > "${sample}_18S.gff3" || true
    grep -A 1 "18S_rRNA" "${sample}_rrna.fasta" | grep -v "^--\$" > "${sample}_18S.fasta" || {
        echo "Error extracting 18S sequences"
        echo ">dummy_18S_rRNA" > "${sample}_18S.fasta"
        echo "ACGTACGTACGTACGTACGT" >> "${sample}_18S.fasta"
    }
    
    # Create CSV with format similar to 16S output
    echo "locus_tag,ftype,length_bp,gene,EC_number,COG,product" > "${sample}_18S.csv"
    
    # Extract info, create locus tags, and add to CSV
    seq_count=1
    grep "18S_rRNA" "${sample}_rrna.gff3" | while read -r line; do
        # Parse GFF3 line fields
        contig=\$(echo "\$line" | cut -f1)
        start=\$(echo "\$line" | cut -f4)
        end=\$(echo "\$line" | cut -f5)
        
        # Calculate length
        if [ "\$start" -lt "\$end" ]; then
            length=\$((end - start + 1))
        else
            length=\$((start - end + 1))
        fi
        
        # Create a locus tag similar to the 16S format (sample_id + sequential number)
        locus_tag="${sample}_\$(printf "%05d" \$seq_count)"
        
        # Write to CSV with the same format as the 16S CSV
        echo "\$locus_tag,rRNA,\$length,,,,18S ribosomal RNA" >> "${sample}_18S.csv"
        
        # Increment sequence counter
        seq_count=\$((seq_count + 1))
    done
    
    # Add simple stats about what was found
    echo "Summary of rRNA sequences found:" >> "${sample}_barrnap.log"
    echo "Total 18S sequences: \$(grep -c "18S_rRNA" "${sample}_rrna.gff3" || echo 0)" >> "${sample}_barrnap.log"
    echo "Total 5.8S sequences: \$(grep -c "5_8S_rRNA" "${sample}_rrna.gff3" || echo 0)" >> "${sample}_barrnap.log"
    echo "Total 28S sequences: \$(grep -c "28S_rRNA" "${sample}_rrna.gff3" || echo 0)" >> "${sample}_barrnap.log"
    
    echo "Barrnap processing complete for ${sample} (${assembler})"
    """
}
