#!/bin/bash
#BSUB -J cleaned_fastqc[1-30]
#BSUB -o logs/cleaned_fastqc_%J_%I.out
#BSUB -e logs/cleaned_fastqc_%J_%I.err
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -M 40000
#BSUB -W 24:00
#BSUB -q gage

# Purpose: Run FastQC on paired-end reads after fastp trimming and filtering.
# Input: sample_list.txt with one sample name per line.
# Expected read layout: 00_trimmedData/<sample>_trimmed_R1.fq.gz and
#                       00_trimmedData/<sample>_trimmed_R2.fq.gz
# Output: FastQC reports in FastQC_after_trim/
# Requirement: The logs/ directory must exist before submitting this LSF job.
# Replace /path/to/ASE with the path to the required Conda environment.

source ~/.bashrc
conda activate /path/to/ASE

sample_list="$PWD/sample_list.txt"
trimmed_reads_dir="$PWD/00_trimmedData"
fastqc_out_dir="$PWD/FastQC_after_trim"

sample=$(sed -n "${LSB_JOBINDEX}p" "$sample_list")

if [[ -z "$sample" ]]; then
    echo "ERROR: No sample found for LSF array index ${LSB_JOBINDEX}." >&2
    exit 1
fi

r1="${trimmed_reads_dir}/${sample}_trimmed_R1.fq.gz"
r2="${trimmed_reads_dir}/${sample}_trimmed_R2.fq.gz"

if [[ ! -f "$r1" || ! -f "$r2" ]]; then
    echo "ERROR: Missing trimmed FASTQ file for sample ${sample}." >&2
    echo "Expected: $r1" >&2
    echo "Expected: $r2" >&2
    exit 1
fi

mkdir -p "$fastqc_out_dir"

fastqc -o "$fastqc_out_dir" "$r1" "$r2"
