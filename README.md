# nxf-pbassembly

## Introduction
A basic pipeline for assembling genomes of isolates sequenced with PacBio HiFi. 

## Quick start
```
git clone https://github.com/Gayathri-Guduru/pacbio-assembly-nextflow.git
cd nxf-pbassembly

If running locally or VM:
nextflow run main.nf -profile slurm --input_type fastq

Use the **--resume ** option to continue the pipeline from where it left off.

If running on HPC:
sbatch main.sh
```


