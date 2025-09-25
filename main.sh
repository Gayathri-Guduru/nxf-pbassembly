#!/bin/bash
#SBATCH --job-name=pbassembly
#SBATCH --output=pbassembly_%j.out
#SBATCH --error=pbassembly_%j.err
#SBATCH --time=60:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=5G
#SBATCH --partition=devel

set -euo pipefail

echo "=== Job started on $(date) ==="
echo "Running on node: $(hostname)"

# module load apptainer (if required)
if ! command -v apptainer &> /dev/null; then
    echo "ERROR: Apptainer is not available in the environment."
    exit 1
fi

echo "Apptainer version: $(apptainer --version)"
echo "Starting Nextflow..."

nextflow run main.nf -profile slurm --input_type fastq -resume

echo "=== Job finished on $(date) ==="
