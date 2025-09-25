# nxf-pbassembly

## Introduction
The nxf-pbassembly pipeline processes PacBio HiFi reads to perform genome assembly and annotation for both bacterial and fungal samples.

## Quick start
```
git clone https://github.com/Gayathri-Guduru/pacbio-assembly-nextflow.git
cd nxf-pbassembly
```

### If running locally or VM:

```nextflow run main.nf -profile slurm --input_type fastq ```

### If running on HPC:
``` sbatch main.sh ```

## Input

The input file is ```samplesheet2.csv``` present in the ```data``` folder. It consists of 9 columns 
- sampleID - The species names (eg:Streptococcus_mutans)
- species_type - The species type (Bacteria/Fungi/Protozoa)
- barcode - ID labelled for the species while sequencing
- bam - Bam files of the species
- reference: Links to the genome sequence (.fna.gz file) downloaded from NCBI.(ex: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/738/105/GCF_009738105.1_ASM973810v1/GCF_009738105.1_ASM973810v1_genomic.fna.gz))
- gff: Links to the corresponding annotation file (.gff.gz file) from NCBI. (ex: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/738/105/GCF_009738105.1_ASM973810v1/GCF_009738105.1_ASM973810v1_genomic.gff.gz)) 
- protein_reference: Links to protein FASTA (.faa.gz) files from NCBI. If no protein reference is used, the value is NA. 
  **Note: Protein references are used only for fungal species, not bacterial.**
- fastq: Links to raw sequencing reads (.fastq.gz files). In this case, files are stored in an S3 bucket.
- genome_size: Estimated genome size (in Mb).

## Pipeline Workflow

The pipeline is composed of multiple processing steps, each handled by a separate Nextflow (.nf) module present in the ```processes``` folder. The workflow typically proceeds in the following order:

- fastqc – Performs quality control on raw sequencing reads.
- longqc – Additional QC tailored for long-read sequencing data.
- bam2fastq – Converts BAM files to FASTQ format. If the input files are already fastq, then use ```--input_type fastq```. The default takes ```bam``` as input.
  
- assembly (flye / hifiasm) – Assembles genomes from PACBIO data using assembly tools (e.g., Flye and hifiasm).
  - Inputs: FASTQ files
  - HIFIASM Outputs: Primary contigs (*.asm.bp.p_ctg.gfa), alternate contigs, unitigs, and an assembly status file (${sample}.assembly_status).
  - FLYE Outputs: Primary contigs (assembly.fasta), assembly graph (assembly_graph.gfa), and other metadata files, plus an assembly status file (${sample}.assembly_status).
 
- gfa_to_fasta – Unlike Flye, HIFIASM only outputs gfa and not fasta file, conversion of Primary contigs (*.asm.bp.p_ctg.gfa) into FASTA format is needed for further downstream annotation.
- quast – Evaluates the quality of genome assemblies (N50, completeness, etc.). Uses reference and GFF files from the samplesheet.
   -  Inputs: FASTA files
   -  Outputs: QUAST output directory (${sample}) and TSV reports.
- prokka – Performs genome annotation to identify coding sequences, functional elements and 16S rRNA sequences (primarily for bacteria).
   -  Inputs: FASTA files
   -  Outputs: Annotation files (GFF, GBK, FNA, FAA, FFN, TSV, etc.).
- barrnap – Identifies rRNA sequences (e.g., 16S, 23S, 5S).
- extract_16S_sequences – Extracts 16S rRNA sequences for downstream phylogenetic or taxonomic analysis.
- metaeuk – Performs protein-level annotation, primarily for fungal species.

