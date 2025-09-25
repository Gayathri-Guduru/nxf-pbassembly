# nxf-pbassembly

## Introduction
A basic pipeline for assembling genomes of isolates sequenced with PacBio HiFi. 

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



