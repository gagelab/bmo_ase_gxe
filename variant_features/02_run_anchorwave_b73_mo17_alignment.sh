#!/bin/bash
#BSUB -J anchorwave_alignment
#BSUB -o logs/anchorwave_alignment_%J.out
#BSUB -e logs/anchorwave_alignment_%J.err
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -M 300000
#BSUB -W 24:00

# Purpose: Generate a whole-genome alignment between B73 and Mo17 with AnchorWave.
# Inputs: B73 and Mo17 FASTA files and AnchorWave input files from step 1.
# Outputs: Anchor, MAF, and filtered MAF files in 02_output/02_anchorwave_alignment/
# Requirement: The logs/ directory must exist before submitting this LSF job.
# Replace /path/to/anchorwave with the path to the required Conda environment.

source ~/.bashrc
conda activate /path/to/anchorwave

input_dir="$PWD/01_input"
anchorwave_input_dir="$PWD/02_output/01_anchorwave_input"
alignment_output_dir="$PWD/02_output/02_anchorwave_alignment"

b73_gff="${input_dir}/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"
b73_fasta="${input_dir}/Zm-B73-REFERENCE-NAM-5.0.fa"
mo17_fasta="${input_dir}/Zm-Mo17-REFERENCE-CAU-2.0.fa"

cds_fasta="${anchorwave_input_dir}/B73_cds.fa"
b73_sam="${anchorwave_input_dir}/B73.sam"
mo17_sam="${anchorwave_input_dir}/Mo17.sam"

output_prefix="${alignment_output_dir}/Mo17_to_B73v5"

mkdir -p "$alignment_output_dir"

anchorwave genoAli \
    -i "$b73_gff" \
    -as "$cds_fasta" \
    -r "$b73_fasta" \
    -a "$mo17_sam" \
    -ar "$b73_sam" \
    -s "$mo17_fasta" \
    -n "${output_prefix}.anchors" \
    -o "${output_prefix}.maf" \
    -f "${output_prefix}.f.maf" \
    -w 38000 \
    -fa3 200000 \
    -B -6 \
    -O1 -8 \
    -E1 -2 \
    -O2 -75 \
    -E2 -1 \
    -IV
