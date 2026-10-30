#!/bin/bash
#BSUB -J reference_indexes
#BSUB -o logs/reference_indexes_%J.out
#BSUB -e logs/reference_indexes_%J.err
#BSUB -n 16
#BSUB -W 24:00
#BSUB -R "span[hosts=1]"

# Purpose: Prepare the reference indexes required by STAR and GATK.
# Inputs: Reference FASTA and GTF files in 01_genomeFiles/
# Outputs: STAR genome index in 02_genomeIndex/ and FASTA .fai and .dict files
# Requirement: The logs/ directory must exist before submitting this LSF job.
# Replace /path/to/ASE with the path to the required Conda environment.
# Replace /path/to/picard.jar with the path to the Picard JAR file.
# STAR, samtools, and Java must be available in the environment or command PATH.

source ~/.bashrc
conda activate /path/to/ASE

genome_files_dir="$PWD/01_genomeFiles"
genome_index_dir="$PWD/02_genomeIndex"
picard_jar="/path/to/picard.jar"
threads=16

fasta="${genome_files_dir}/Zm-B73-REFERENCE-NAM-5.0.fa"
gtf="${genome_files_dir}/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gtf"
fasta_index="${fasta}.fai"
sequence_dict="${genome_files_dir}/Zm-B73-REFERENCE-NAM-5.0.dict"

if [[ ! -f "$fasta" ]]; then
    echo "ERROR: Reference FASTA not found: $fasta" >&2
    exit 1
fi

if [[ ! -f "$gtf" ]]; then
    echo "ERROR: Reference GTF not found: $gtf" >&2
    exit 1
fi

if [[ ! -f "$picard_jar" ]]; then
    echo "ERROR: Picard JAR not found: $picard_jar" >&2
    exit 1
fi

mkdir -p "$genome_index_dir"

if [[ ! -f "$fasta_index" ]]; then
    echo "Creating the FASTA index."

    samtools faidx "$fasta"
else
    echo "FASTA index already exists: $fasta_index"
fi

if [[ ! -f "$sequence_dict" ]]; then
    echo "Creating the reference sequence dictionary."

    java -jar "$picard_jar" CreateSequenceDictionary \
        R="$fasta" \
        O="$sequence_dict"
else
    echo "Reference sequence dictionary already exists: $sequence_dict"
fi

if [[ ! -f "${genome_index_dir}/Genome" ]]; then
    echo "Creating the STAR genome index."

    STAR \
        --runThreadN "$threads" \
        --runMode genomeGenerate \
        --genomeDir "$genome_index_dir" \
        --genomeFastaFiles "$fasta" \
        --sjdbGTFfile "$gtf" \
        --sjdbOverhang 100
else
    echo "STAR genome index already exists: $genome_index_dir"
fi

echo "Reference index preparation completed."
