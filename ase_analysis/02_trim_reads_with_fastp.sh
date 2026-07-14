#!/bin/bash
#BSUB -J trim_reads[1-30]
#BSUB -o logs/fastp_%J_%I.out
#BSUB -e logs/fastp_%J_%I.err
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -M 40000
#BSUB -W 24:00
#BSUB -q gage

# Purpose: Trim and filter raw paired-end reads with fastp.
# Input: sample_list.txt with one sample name per line.
# Expected read layout: 00_rawData/<sample>/<sample>_1.fq.gz and <sample>_2.fq.gz
# Outputs: Cleaned reads in 00_trimmedData/ and fastp reports in Fastp/
# Requirement: The logs/ directory must exist before submitting this LSF job.
# Replace /path/to/ASE with the path to the required Conda environment.

source ~/.bashrc
conda activate /path/to/ASE

sample_list="$PWD/sample_list.txt"
raw_reads_dir="$PWD/00_rawData"
trimmed_reads_dir="$PWD/00_trimmedData"
fastp_report_dir="$PWD/Fastp"
threads=4

sample=$(sed -n "${LSB_JOBINDEX}p" "$sample_list")

if [[ -z "$sample" ]]; then
    echo "ERROR: No sample found for LSF array index ${LSB_JOBINDEX}." >&2
    exit 1
fi

r1="${raw_reads_dir}/${sample}/${sample}_1.fq.gz"
r2="${raw_reads_dir}/${sample}/${sample}_2.fq.gz"
r1_out="${trimmed_reads_dir}/${sample}_trimmed_R1.fq.gz"
r2_out="${trimmed_reads_dir}/${sample}_trimmed_R2.fq.gz"
json_report="${fastp_report_dir}/${sample}.json"
html_report="${fastp_report_dir}/${sample}.html"

if [[ ! -f "$r1" || ! -f "$r2" ]]; then
    echo "ERROR: Missing raw FASTQ file for sample ${sample}." >&2
    echo "Expected: $r1" >&2
    echo "Expected: $r2" >&2
    exit 1
fi

mkdir -p "$trimmed_reads_dir" "$fastp_report_dir"

fastp \
    -j "$json_report" \
    -h "$html_report" \
    -w "$threads" \
    -i "$r1" \
    -I "$r2" \
    -o "$r1_out" \
    -O "$r2_out"
