#!/bin/bash
#BSUB -J raw_fastqc[1-30]
#BSUB -o logs/raw_fastqc_%J_%I.out
#BSUB -e logs/raw_fastqc_%J_%I.err
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -M 40000
#BSUB -W 24:00
#BSUB -q gage

# Purpose: Run FastQC on raw paired-end reads.
# Input: sample_list.txt with one sample name per line.
# Expected read layout: 00_rawData/<sample>/<sample>_1.fq.gz and <sample>_2.fq.gz
# Output: FastQC reports in FastQC_before_trim/
# Replace /path/to/ASE with the path to the required Conda environment.

source ~/.bashrc
conda activate /path/to/ASE

sample_list="$PWD/sample_list.txt"
raw_reads_dir="$PWD/00_rawData"
fastqc_out_dir="$PWD/FastQC_before_trim"

sample=$(sed -n "${LSB_JOBINDEX}p" "$sample_list")

if [[ -z "$sample" ]]; then
    echo "ERROR: No sample found for LSF array index ${LSB_JOBINDEX}." >&2
    exit 1
fi

r1="${raw_reads_dir}/${sample}/${sample}_1.fq.gz"
r2="${raw_reads_dir}/${sample}/${sample}_2.fq.gz"

if [[ ! -f "$r1" || ! -f "$r2" ]]; then
    echo "ERROR: Missing raw FASTQ file for sample ${sample}." >&2
    echo "Expected: $r1" >&2
    echo "Expected: $r2" >&2
    exit 1
fi

mkdir -p "$fastqc_out_dir"

fastqc -o "$fastqc_out_dir" "$r1" "$r2"
