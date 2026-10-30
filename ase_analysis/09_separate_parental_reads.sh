#!/bin/bash
#BSUB -J separate_parental_reads[1-30]
#BSUB -o logs/separate_parental_reads_%J_%I.out
#BSUB -e logs/separate_parental_reads_%J_%I.err
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -M 40000
#BSUB -W 48:00

# Purpose: Separate Parent1- and Parent2-associated reads using STAR vA tags.
# Input: sample_list.txt with one sample name per line.
# Expected BAM: 07_starMappingWasp/pass2/<sample>_WASPfilterPASS_RG_sort.bam
# Outputs: Parent-specific BAM files in 07_starMappingWasp/pass2/
# Requirement: The logs/ directory must exist before submitting this LSF job.
# Replace /path/to/ASE with the path to the required Conda environment.
# samtools must be available in the environment or command PATH.

source ~/.bashrc
conda activate /path/to/ASE

sample_list="$PWD/sample_list.txt"
star_wasp_dir="$PWD/07_starMappingWasp/pass2"

sample=$(sed -n "${LSB_JOBINDEX}p" "$sample_list")

if [[ -z "$sample" ]]; then
    echo "ERROR: No sample found for LSF array index ${LSB_JOBINDEX}." >&2
    exit 1
fi

input_bam="${star_wasp_dir}/${sample}_WASPfilterPASS_RG_sort.bam"
parent1_bam="${star_wasp_dir}/${sample}_WASPfilterPASS_RG_sort.vA1.bam"
parent2_bam="${star_wasp_dir}/${sample}_WASPfilterPASS_RG_sort.vA2.bam"

if [[ ! -f "$input_bam" ]]; then
    echo "ERROR: Input BAM not found: $input_bam" >&2
    exit 1
fi

echo "[${sample}] Extracting Parent1-associated reads."

# The vA tag and field position were verified for these STAR outputs.
samtools view -h "$input_bam" |
    grep -E 'vA:B:c,1|^@' |
    awk '$1 ~ /^@/ || $12 !~ /,2|,3/' |
    samtools view -Shb - \
    > "$parent1_bam"

if [[ ! -s "$parent1_bam" ]]; then
    echo "ERROR: Parent1 BAM was not created: $parent1_bam" >&2
    exit 1
fi

echo "[${sample}] Extracting Parent2-associated reads."

# The vA tag and field position were verified for these STAR outputs.
samtools view -h "$input_bam" |
    grep -E 'vA:B:c,2|^@' |
    awk '$1 ~ /^@/ || $12 !~ /,1|,3/' |
    samtools view -Shb - \
    > "$parent2_bam"

if [[ ! -s "$parent2_bam" ]]; then
    echo "ERROR: Parent2 BAM was not created: $parent2_bam" >&2
    exit 1
fi

echo "[${sample}] Completed:"
echo "  Parent1: $parent1_bam"
echo "  Parent2: $parent2_bam"
