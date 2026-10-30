#!/bin/bash
#BSUB -J anchorwave_input
#BSUB -o logs/anchorwave_input_%J.out
#BSUB -e logs/anchorwave_input_%J.err
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -M 40000
#BSUB -W 48:00

# Purpose: Extract B73 CDS anchors and map them to the B73 and Mo17 genomes.
# Inputs: B73 GFF3, B73 v5 FASTA, and Mo17 CAU 2.0 FASTA in 01_input/
# Outputs: CDS FASTA and B73/Mo17 SAM files in 02_output/01_anchorwave_input/
# Requirement: The logs/ directory must exist before submitting this LSF job.
# Replace /path/to/anchorwave with the path to the required Conda environment.

source ~/.bashrc
conda activate /path/to/anchorwave

input_dir="$PWD/01_input"
output_dir="$PWD/02_output/01_anchorwave_input"
threads=10

b73_gff="${input_dir}/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"
b73_fasta="${input_dir}/Zm-B73-REFERENCE-NAM-5.0.fa"
mo17_fasta="${input_dir}/Zm-Mo17-REFERENCE-CAU-2.0.fa"

cds_fasta="${output_dir}/B73_cds.fa"
b73_sam="${output_dir}/B73.sam"
mo17_sam="${output_dir}/Mo17.sam"

mkdir -p "$output_dir"

anchorwave gff2seq \
    -i "$b73_gff" \
    -r "$b73_fasta" \
    -o "$cds_fasta"

minimap2 \
    -x splice \
    -t "$threads" \
    -k 12 \
    -a \
    -p 0.4 \
    -N 20 \
    "$b73_fasta" \
    "$cds_fasta" \
    > "$b73_sam"

minimap2 \
    -x splice \
    -t "$threads" \
    -k 12 \
    -a \
    -p 0.4 \
    -N 20 \
    "$mo17_fasta" \
    "$cds_fasta" \
    > "$mo17_sam"

echo "AnchorWave input preparation completed:"
echo "  CDS anchors: $cds_fasta"
echo "  B73 mapping: $b73_sam"
echo "  Mo17 mapping: $mo17_sam"
