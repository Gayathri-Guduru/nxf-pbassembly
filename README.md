# nxf-pbassembly

## Introduction
The nxf-pbassembly pipeline processes PacBio HiFi reads to perform genome assembly and annotation for both bacterial and fungal samples.

## Quick start
```
git clone https://github.com/Gayathri-Guduru/pacbio-assembly-nextflow.git
cd nxf-pbassembly
```

- If running locally or VM:
```nextflow run main.nf -profile slurm --input_type fastq ```

- If running on HPC:
``` sbatch main.sh ```

## Input

The input file is ```samplesheet2.csv``` present in the ```data``` folder. It consists of 9 columns 
- ***sampleID*** - The species names (eg:Streptococcus_mutans)
- ***species_type*** - The species type (Bacteria/Fungi/Protozoa)
- ***barcode*** - ID labelled for the species while sequencing
- ***bam*** - Bam files of the species
- ***reference***: Links to the genome sequence (.fna.gz file) downloaded from NCBI.(ex: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/738/105/GCF_009738105.1_ASM973810v1/GCF_009738105.1_ASM973810v1_genomic.fna.gz))
- ***gff***: Links to the corresponding annotation file (.gff.gz file) from NCBI. (ex: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/738/105/GCF_009738105.1_ASM973810v1/GCF_009738105.1_ASM973810v1_genomic.gff.gz)) 
- ***protein_reference***: Links to protein FASTA (.faa.gz) files from NCBI. If no protein reference is used, the value is NA. 
  **Note: Protein references are used only for fungal species, not bacterial.**
- ***fastq***: Links to raw sequencing reads (.fastq.gz files). In this case, files are stored in an S3 bucket.
- ***genome_size***: Estimated genome size (in Mb).

## Pipeline Workflow

<img width="519" height="701" alt="pacbio_assembly_workflow" src="https://github.com/user-attachments/assets/12697d41-77a1-4318-ad5c-a188f2f5cd6e" />

The pipeline is composed of multiple processing steps, each handled by a separate Nextflow (.nf) module present in the ```processes``` folder. The workflow typically proceeds in the following order:

- ***fastqc*** – Performs quality control on raw sequencing reads.
- ***longqc*** – Additional QC tailored for long-read sequencing data.
- ***bam2fastq*** – Converts BAM files to FASTQ format. If the input files are already fastq, then use ```--input_type fastq```. The default takes ```bam``` as input.
  
- ***assembly (flye / hifiasm)*** – Assembles genomes from PACBIO data using assembly tools (e.g., Flye and hifiasm).
  - **Inputs**: FASTQ files
  - **HIFIASM Outputs**: Primary contigs (*.asm.bp.p_ctg.gfa), alternate contigs, unitigs, and an assembly status file (${sample}.assembly_status).
  - **FLYE Outputs**: Primary contigs (assembly.fasta), assembly graph (assembly_graph.gfa), and other metadata files, plus an assembly status file (${sample}.assembly_status).
 
- ***gfa_to_fasta*** – Unlike Flye, HIFIASM only outputs gfa and not fasta file, conversion of Primary contigs (*.asm.bp.p_ctg.gfa) into FASTA format is needed for further downstream annotation.
- ***quast*** – Evaluates the quality of genome assemblies (N50, completeness, etc.). Uses reference and GFF files from the samplesheet.
   -  **Inputs**: FASTA files
   -  **Outputs**: QUAST output directory (${sample}) and TSV reports.
     
- ***prokka*** – Performs genome annotation to identify coding sequences, functional elements and 16S rRNA sequences (primarily for bacteria).
   -  **Inputs**: FASTA files
   -  **Outputs**: Annotation files (GFF, GBK, FNA, FAA, FFN, TSV, etc.).
     
- ***extract_16S_sequences*** – Runs only on ```bacterial``` samples. Extracts 16S rRNA sequences for downstream phylogenetic or taxonomic analysis.
  - **Inputs**: TSV and FFN files from PROKKA.
  - **Outputs**: 16S_csv, 16S_fasta files containing extracted 16S sequences and metadata.
   
- ***metaeuk*** – Performs protein-level gene prediction for fungal samples using MetaEuk with a provided protein reference.
  - **Inputs**: FASTA files from GFA_TO_FASTA and protein reference from the samplesheet.
  - **Outputs**: Directory (metaeuk_output) containing gene prediction results (FASTA, GFF, etc.).

- ***barrnap*** – Annotates rRNA sequences (e.g., 18S, 5.8S, 28S) in fungal assemblies.
    - **Inputs**: FASTA files from GFA_TO_FASTA.
    - **Outputs**: GFF3 file (${sample}_rrna.gff3), rRNA FASTA, 18S-specific FASTA and CSV, and a log file.


