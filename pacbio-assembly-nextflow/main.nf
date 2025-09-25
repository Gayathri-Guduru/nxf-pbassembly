#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters
params.input = params.input ?: "./data/samplesheet2.csv"
params.input_type = params.input_type ?: "bam"
params.cpus = params.cpus ?: 12
params.publish_dir = params.publish_dir ?: "./results"
params.run_hifiasm = params.run_hifiasm ?: true
params.run_flye = params.run_flye ?: true
params.run_prokka = params.run_prokka ?: true
params.run_metaeuk = params.run_metaeuk ?: true
params.run_barrnap = params.run_barrnap ?: true
params.protein_reference = params.protein_reference ?: "NA"

// Add explicit logging of parameters
log.info ""
log.info "==================================="
log.info "  PacBio Assembly Pipeline         "
log.info "==================================="
log.info "Input file       : ${params.input}"
log.info "Input type       : ${params.input_type}"
log.info "CPUs             : ${params.cpus}"
log.info "Results dir      : ${params.publish_dir}"
log.info "Run HiFiAsm      : ${params.run_hifiasm}"
log.info "Run Flye         : ${params.run_flye}"
log.info "Run Prokka       : ${params.run_prokka}"
log.info "Run MetaEuk      : ${params.run_metaeuk}"
log.info "Run Barrnap      : ${params.run_barrnap}"
log.info "Protein Ref      : ${params.protein_reference}"
log.info "Flye meta fallback: Enabled"
log.info "==================================="
log.info ""

// Convert input_type to lowercase for case-insensitive comparison
def input_type = params.input_type.toLowerCase()

if (!params.input) {
    exit 1, "--input sample sheet is a required input!"
}
if (input_type != "bam" && input_type != "fastq") {
    exit 1, "Invalid --input_type parameter. Please use 'bam' or 'fastq'. You specified: ${params.input_type}"
}
if (!params.run_hifiasm && !params.run_flye) {
    exit 1, "At least one assembler (hifiasm or flye) must be enabled!"
}

// Include processes
include { BAM2FASTQ } from './processes/bam2fastq'
include { FASTQC } from './processes/fastqc'
include { ASSEMBLY_HIFIASM } from './processes/hifiasm'
include { ASSEMBLY_FLYE } from './processes/flye'
include { QUAST as QUAST_HIFIASM; QUAST as QUAST_FLYE } from './processes/quast'  // Create distinct aliases
include { GFA_TO_FASTA } from './processes/gfa_to_fasta'
include { PROKKA } from './processes/prokka'
include { EXTRACT_16S } from './processes/extract_16S_sequences'
include { METAEUK } from './processes/metaeuk.nf'
include { BARRNAP } from './processes/barrnap'

workflow {

    log.info "Starting workflow with input_type: ${input_type}"

    // Read complete metadata including protein_reference with sanitized sample names
    ch_full_metadata = Channel
        .fromPath(file(params.input))
        .splitCsv(header: true)
        .map { row ->
            def protein_ref = row.containsKey('protein_reference') ? row.protein_reference : "NA"
            // Clean sample ID by removing trailing spaces and replacing spaces with underscores
            def cleanSampleID = "$row.sampleID".trim().replaceAll("\\s+", "_")

            tuple(
                cleanSampleID,  // use cleaned sample ID
                file("$row.bam"),
                file("$row.reference"),
                file("$row.gff"),
                protein_ref
            )
        }

    if (input_type == "bam") {
        log.info "Processing BAM input files"
        ch_metadata = ch_full_metadata.map { sample, bam, reference, gff, protein_ref ->
            tuple(sample, bam, reference, gff)
        }
        BAM2FASTQ(ch_metadata)
        ch_reads = BAM2FASTQ.out.reads
            .join(ch_metadata, by: 0)
            .map { sample, reads, metadata_sample, bam, reference, gff ->
                def genome_size = null
                tuple(sample, reads, genome_size)
            }
    } else {
        log.info "Processing FASTQ input files"
        ch_fastq = Channel
            .fromPath(file(params.input))
            .splitCsv(header: true)
            .map { row ->
                log.info "Processing row: $row"
                def genome_size = row.containsKey('genome_size') ? row.genome_size : null
                // Clean sample ID by removing trailing spaces and replacing spaces with underscores
                def cleanSampleID = "$row.sampleID".trim().replaceAll("\\s+", "_")

                tuple(
                    cleanSampleID,  // use cleaned sample ID
                    file("$row.fastq"),
                    genome_size
                )
            }
        ch_reads = ch_fastq
        ch_metadata = ch_full_metadata.map { sample, bam, reference, gff, protein_ref ->
            tuple(sample, bam, reference, gff)
        }
    }

    ch_reads_fastqc = ch_reads.map { sample, reads, genome_size -> tuple(sample, reads) }
    ch_reads_hifiasm = ch_reads.map { sample, reads, genome_size -> tuple(sample, reads) }
    ch_reads_flye = ch_reads

    FASTQC(ch_reads_fastqc)

    ch_hifiasm_contigs = Channel.empty()
    ch_flye_contigs = Channel.empty()

    // Create a channel to collect failed samples
    ch_failed_samples = Channel.empty()

    // --- HiFiasm Assembly ---
    if (params.run_hifiasm) {
        ASSEMBLY_HIFIASM(ch_reads_hifiasm)

        // Create a status channel from the assembly_status output
        ch_hifiasm_status = ASSEMBLY_HIFIASM.out.assembly_status
            .map { sample, status_file ->
                def status_text = status_file.text.trim()
                tuple(sample, status_text)
            }

        // Only use contigs that exist and are valid
        ch_hifiasm_contigs = ASSEMBLY_HIFIASM.out.primary_contigs
            .filter { sample, contigs -> 
                // Verify file exists and is not empty
                def valid = contigs.exists() && contigs.size() > 0
                if (!valid) {
                    log.warn "HiFiasm did not produce valid contigs for sample: ${sample}, skipping downstream processes"
                }
                return valid
            }
            .map { sample, contigs ->
                log.info "HiFiasm assembly succeeded for sample: ${sample}"
                tuple(sample, contigs, 'hifiasm')
            }

        // Use conditional process execution without isEmpty()
        ch_hifiasm_contigs_verified = ch_hifiasm_contigs.map { it -> 
            log.info "Valid HiFiasm contig found for: ${it[0]}"
            return it 
        }
        
        // Run GFA_TO_FASTA on verified contigs
        // We don't need to check if the channel is empty - Nextflow will handle this
        GFA_TO_FASTA(ch_hifiasm_contigs_verified)
        ch_hifiasm_contigs = GFA_TO_FASTA.out

        // More robust QUAST input preparation with explicit filtering
        ch_hifiasm_quast_input_with_metadata = ch_hifiasm_contigs
            .join(ch_metadata, remainder: true)
            .filter { data -> 
                // Check if we have all 6 elements and none are null
                if (data.size() != 6 || data.any { it == null }) {
                    log.warn "Incomplete data for QUAST for sample: ${data[0]}, skipping"
                    return false
                }
                return true
            }
         
        // Separate assembly from metadata to match QUAST process input expectations   
        ch_hifiasm_quast_assemblies = ch_hifiasm_quast_input_with_metadata.map { sample, assembly, assembler, bam, reference, gff ->
            println "DEBUG - QUAST HiFiasm input for: ${sample}"
            tuple(sample, assembly, assembler)
        }
        
        ch_hifiasm_quast_metadata = ch_hifiasm_quast_input_with_metadata.map { sample, assembly, assembler, bam, reference, gff ->
            tuple(sample, bam, reference, gff)
        }
        
        // No need to check if empty - Nextflow handles this
        QUAST_HIFIASM(ch_hifiasm_quast_assemblies, ch_hifiasm_quast_metadata)
    }

    // --- Flye Assembly ---
    if (params.run_flye) {
        ASSEMBLY_FLYE(ch_reads_flye)

        // Create a status channel from the assembly_status output
        ch_flye_status = ASSEMBLY_FLYE.out.assembly_status
            .map { sample, status_file ->
                def status_text = status_file.text.trim()
                tuple(sample, status_text)
            }

        // Only use contigs that exist and are valid
        ch_flye_contigs = ASSEMBLY_FLYE.out.primary_contigs
            .filter { sample, contigs -> 
                // Verify file exists and is not empty
                def valid = contigs.exists() && contigs.size() > 0
                if (!valid) {
                    log.warn "Flye did not produce valid contigs for sample: ${sample}, skipping downstream processes"
                }
                return valid
            }
            .map { sample, contigs ->
                log.info "Flye assembly succeeded for sample: ${sample}"
                tuple(sample, contigs, 'flye')
            }

        // More robust QUAST input preparation with explicit filtering
        ch_flye_quast_input_with_metadata = ch_flye_contigs
            .join(ch_metadata, remainder: true)
            .filter { data -> 
                // Check if we have all 6 elements and none are null
                if (data.size() != 6 || data.any { it == null }) {
                    log.warn "Incomplete data for QUAST for sample: ${data[0]}, skipping"
                    return false
                }
                return true
            }
            
        // Separate assembly from metadata to match QUAST process input expectations
        ch_flye_quast_assemblies = ch_flye_quast_input_with_metadata.map { sample, assembly, assembler, bam, reference, gff ->
            println "DEBUG - QUAST Flye input for: ${sample}"
            tuple(sample, assembly, assembler)
        }
        
        ch_flye_quast_metadata = ch_flye_quast_input_with_metadata.map { sample, assembly, assembler, bam, reference, gff ->
            tuple(sample, bam, reference, gff)
        }
        
        // Count valid samples for QUAST
        ch_flye_quast_count = ch_flye_quast_assemblies.count()
        
        // Log the count but don't use it to determine whether to run processes
        ch_flye_quast_count.view { count -> 
            if (count == 0) {
                log.warn "No valid Flye assemblies for QUAST analysis" 
            } else {
                log.info "Found ${count} valid Flye assemblies for QUAST analysis"
            }
            return count
        }
        
        // Nextflow will automatically handle empty channels - no need to check
        QUAST_FLYE(ch_flye_quast_assemblies, ch_flye_quast_metadata)
    }

    // Debug output of Flye contigs to see what's available
    ch_flye_contigs.view { sample, contigs, assembler ->
        "DEBUG - Flye contigs for sample: ${sample}, assembler: ${assembler}, file: ${contigs.getName()}"
    }

    // Debug output of HiFiasm contigs to see what's available
    ch_hifiasm_contigs.view { sample, contigs, assembler ->
        "DEBUG - HiFiasm contigs for sample: ${sample}, assembler: ${assembler}, file: ${contigs.getName()}"
    }

    // --- Identify Failed Samples but Don't Exit Pipeline ---
    if (params.run_hifiasm && params.run_flye) {
        // Join the status channels from both assemblers
        ch_sample_statuses = ch_hifiasm_status.join(ch_flye_status, remainder: true, by: 0)
            .map { sample, hifiasm_status, flye_status ->
                def hifiasm_ok = hifiasm_status == "SUCCESS"
                def flye_ok = flye_status == "SUCCESS"

                if (!hifiasm_ok && !flye_ok) {
                    log.error "Sample ${sample} failed with both assemblers"
                    return tuple(sample, "FAILED")
                } else {
                    log.info "Sample ${sample} succeeded with at least one assembler"
                    if (hifiasm_ok) log.info "  - HiFiasm: SUCCESS"
                    else log.info "  - HiFiasm: FAILED"
                    if (flye_ok) log.info "  - Flye: SUCCESS"
                    else log.info "  - Flye: FAILED"
                    return tuple(sample, "SUCCESS")
                }
            }

        // Write the list of failed samples to a file but don't exit
        ch_sample_statuses
            .filter { sample, status -> status == "FAILED" }
            .map { sample, status -> sample }
            .collectFile(name: 'failed_samples.txt', newLine: true, storeDir: params.publish_dir)
            .subscribe {
                if (file(it).text.trim()) {
                    def failedSamples = file(it).text.trim().split('\n')
                    if (failedSamples.size() > 0) {
                        log.error "Both assemblers failed for the following samples: ${failedSamples.join(', ')}"
                        log.warn "Pipeline continuing without these samples: ${failedSamples.join(', ')}"
                    }
                }
            }
    }
    else if (params.run_hifiasm) {
        ch_hifiasm_status
            .filter { sample, status -> status == "FAILED" }
            .map { sample, status -> sample }
            .collectFile(name: 'failed_samples.txt', newLine: true, storeDir: params.publish_dir)
            .subscribe {
                if (file(it).text.trim()) {
                    def failedSamples = file(it).text.trim().split('\n')
                    if (failedSamples.size() > 0) {
                        log.error "HiFiasm failed for the following samples: ${failedSamples.join(', ')}"
                        log.warn "Pipeline continuing without these samples: ${failedSamples.join(', ')}"
                    }
                }
            }
    }
    else if (params.run_flye) {
        ch_flye_status
            .filter { sample, status -> status == "FAILED" }
            .map { sample, status -> sample }
            .collectFile(name: 'failed_samples.txt', newLine: true, storeDir: params.publish_dir)
            .subscribe {
                if (file(it).text.trim()) {
                    def failedSamples = file(it).text.trim().split('\n')
                    if (failedSamples.size() > 0) {
                        log.error "Flye failed for the following samples: ${failedSamples.join(', ')}"
                        log.warn "Pipeline continuing without these samples: ${failedSamples.join(', ')}"
                    }
                }
            }
    }

    // --- Mix All Contigs for PROKKA and other tools ---
    // Only mixing valid contigs from both assemblers
    ch_all_contigs = ch_hifiasm_contigs.mix(ch_flye_contigs)

    // --- Run Prokka on All Assemblies if Enabled ---
    if (params.run_prokka) {
        // Track the count for logging only
        ch_all_contigs.count().view { count -> 
            if (count == 0) {
                log.warn "No valid assemblies available for Prokka annotation" 
            } else {
                log.info "Found ${count} valid assemblies for Prokka annotation"
            }
            return count
        }
        
        // No need to check empty channel - just run PROKKA and let Nextflow handle empty channels
        ch_prokka_input = ch_all_contigs
        PROKKA(
            ch_prokka_input.map { sample, assembly, assembler -> tuple(sample, assembly) },
            ch_prokka_input.map { sample, assembly, assembler -> assembler }
        )
        
        // Only proceed with 16S extraction if Prokka outputs exist
        ch_16s_input = PROKKA.out.for_16s
            .filter { sample, tsv, ffn, assembler -> 
                tsv.exists() && ffn.exists() 
            }
        
        // Nextflow will handle empty channel correctly - no need to explicitly check
        EXTRACT_16S(ch_16s_input)
        
        PROKKA.out.for_16s
            .filter { sample, tsv, ffn, assembler -> tsv.exists() && ffn.exists() }
            .view { sample, tsv, ffn, assembler ->
                "DEBUG - PROKKA output for 16S extraction: sample=${sample}, tsv=${tsv.getName()}, ffn=${ffn.getName()}, assembler=${assembler}"
            }
    }

    // Create a channel with the protein reference info for sample filtering
    ch_protein_refs = ch_full_metadata
        .map { sample, bam, reference, gff, protein_ref ->
            log.info "Sample ${sample} has protein reference: '${protein_ref}'"
            tuple(sample, protein_ref)
        }

    // --- Run Barrnap only on fungal samples if enabled ---
    if (params.run_barrnap) {
        // Filter assemblies to only include fungal samples with valid assemblies
        ch_barrnap_input = ch_all_contigs
            .combine(ch_protein_refs, by: 0)
            .filter { data ->
                // Skip if any data is missing
                if (data.size() != 4 || data.any { it == null }) {
                    log.warn "Incomplete data for Barrnap for sample: ${data[0]}, skipping"
                    return false
                }
                
                def sample = data[0]
                def assembly = data[1]
                def assembler = data[2]
                def protein_ref = data[3]
                
                // Check if it's a fungal sample
                def is_fungal = protein_ref != "NA" && protein_ref != ""
                
                // Check if assembly file exists and is not empty
                def has_valid_assembly = assembly.exists() && assembly.size() > 0
                
                if (!has_valid_assembly) {
                    log.warn "Invalid assembly file for Barrnap for sample: ${sample}, skipping"
                    return false
                }
                
                if (is_fungal) {
                    log.info "FUNGAL SAMPLE: Running Barrnap on ${sample} with ${assembler} assembler"
                } else {
                    log.info "BACTERIAL SAMPLE: Skipping Barrnap for ${sample} with ${assembler} assembler"
                }
                
                return is_fungal && has_valid_assembly
            }
            .map { sample, assembly, assembler, protein_ref ->
                tuple(sample, assembly, assembler)
            }

        // Debug view of what's being sent to Barrnap
        ch_barrnap_input.view { sample, assembly, assembler ->
            "BARRNAP INPUT: ${sample}_${assembler} - ${assembly.getName()}"
        }

        // Just log the count for information purposes
        ch_barrnap_input.count().view { count -> 
            if (count == 0) {
                log.warn "No valid fungal assemblies for Barrnap analysis" 
            } else {
                log.info "Found ${count} valid fungal assemblies for Barrnap analysis"
            }
            return count
        }
        
        // Run Barrnap (Nextflow handles empty channels)
        BARRNAP(ch_barrnap_input)

        // Debug output for monitoring - FIXED THE LAMBDA SYNTAX HERE
        BARRNAP.out.only_18s_fasta
            .filter { fasta -> fasta.exists() }
            .view { fasta ->
                "DEBUG - Extracted 18S sequences: ${fasta}"
            }

        BARRNAP.out.only_18s_csv
            .filter { csv -> csv.exists() }
            .view { csv ->
                "DEBUG - 18S sequences CSV summary: ${csv}"
            }
    }

    // --- Run MetaEuk if Enabled (UPDATED TO PROCESS BOTH ASSEMBLY TYPES) ---
    if (params.run_metaeuk) {
        // Debug view of protein references
        ch_protein_refs.view { sample, protein_ref ->
            "PROTEIN REF: Sample ${sample}, protein_ref: ${protein_ref}"
        }

        // Debug view of available contigs for MetaEuk
        ch_all_contigs.view { sample, assembly, assembler ->
            "AVAILABLE FOR METAEUK: Sample ${sample}, assembler ${assembler}, file: ${assembly.getName()}"
        }

        // Process contigs for MetaEuk using combine instead of join with improved filtering
        ch_metaeuk_input = ch_all_contigs
            .combine(ch_protein_refs, by: 0)
            .filter { data ->
                // Skip if we're missing data or have null values
                if (data.size() != 4 || data.any { it == null }) {
                    log.warn "Incomplete data for MetaEuk for sample: ${data[0]}, skipping"
                    return false
                }
                
                def sample = data[0]
                def assembly = data[1]
                def assembler = data[2]
                def protein_ref = data[3]
                
                // Check if it's a fungal sample
                def run_metaeuk = protein_ref != "NA" && protein_ref != ""
                
                // Check if assembly file exists and is not empty
                def has_valid_assembly = assembly.exists() && assembly.size() > 0
                
                if (!has_valid_assembly) {
                    log.warn "Invalid assembly file for MetaEuk for sample: ${sample}, skipping"
                    return false
                }
                
                if (run_metaeuk) {
                    log.info "ACCEPTED: Running MetaEuk for fungal sample: ${sample} with ${assembler} assembler"
                } else {
                    log.info "FILTERED OUT: Skipping MetaEuk for bacterial sample: ${sample} with ${assembler} assembler"
                }
                
                return run_metaeuk && has_valid_assembly
            }
            .map { sample, assembly, assembler, protein_ref ->
                log.info "SENDING TO METAEUK: Sample ${sample}, assembler ${assembler}, assembly: ${assembly.getName()}, protein_ref: ${protein_ref}"
                tuple(sample, assembly, assembler, protein_ref)
            }

        // Debug view of what's being sent to MetaEuk
        ch_metaeuk_input.view { sample, contigs, assembler, protein_ref ->
            "METAEUK INPUT: ${sample}_${assembler} - ${contigs.getName()}"
        }

        // Track counts for logging purposes
        ch_metaeuk_input.ifEmpty { 
            log.warn "No fungal samples with valid protein references found for MetaEuk processing" 
        }.count().view { count ->
            if (count > 0) {
                log.info "Found ${count} valid inputs for MetaEuk"
            }
            return count
        }
        
        // Run MetaEuk (Nextflow handles empty channels)
        METAEUK(ch_metaeuk_input)
    }
}
