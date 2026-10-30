#!/bin/bash
#BSUB -J maf_to_paf
#BSUB -o logs/maf_to_paf_%J.out
#BSUB -e logs/maf_to_paf_%J.err
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -M 300000
#BSUB -W 24:00

# Purpose: Convert the B73–Mo17 AnchorWave MAF alignment to PAF for SyRI.
# Input: AnchorWave MAF file from step 2.
# Output: PAF file in 02_output/03_syri_input/
# Requirement: The logs/ directory must exist before submitting this LSF job.
# Replace /path/to/wgatools with the path to the required Conda environment.

source ~/.bashrc
conda activate /path/to/wgatools

alignment_dir="$PWD/02_output/02_anchorwave_alignment"
syri_input_dir="$PWD/02_output/03_syri_input"

input_maf="${alignment_dir}/Mo17_to_B73v5.maf"
output_paf="${syri_input_dir}/Mo17_to_B73v5.paf"

mkdir -p "$syri_input_dir"

wgatools maf2paf \
    "$input_maf" \
    -o "$output_paf"
