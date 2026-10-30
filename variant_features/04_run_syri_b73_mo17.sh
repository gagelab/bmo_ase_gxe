#!/bin/bash
#BSUB -J syri_b73_mo17
#BSUB -o logs/syri_b73_mo17_%J.out
#BSUB -e logs/syri_b73_mo17_%J.err
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -M 300000
#BSUB -W 24:00

# Purpose: Identify genomic differences between B73 and Mo17 with SyRI.
# Input: B73–Mo17 PAF alignment from step 3.
# Outputs: SyRI result files in 02_output/04_syri_output/
# Requirement: The logs/ directory must exist before submitting this LSF job.
# Replace /path/to/syri with the path to the required Conda environment.

source ~/.bashrc
conda activate /path/to/syri

input_dir="$PWD/01_input"
syri_input_dir="$PWD/02_output/03_syri_input"
syri_output_dir="$PWD/02_output/04_syri_output"

b73_fasta="${input_dir}/Zm-B73-REFERENCE-NAM-5.0.fa"
mo17_fasta="${input_dir}/Zm-Mo17-REFERENCE-CAU-2.0.fa"
input_paf="${syri_input_dir}/Mo17_to_B73v5.paf"
output_prefix="Mo17_toB73v5_paf_"

mkdir -p "$syri_output_dir"

syri \
    -c "$input_paf" \
    -F P \
    --dir "$syri_output_dir" \
    --prefix "$output_prefix" \
    --cigar \
    -f \
    --log DEBUG \
    -r "$b73_fasta" \
    -q "$mo17_fasta"
