#!/bin/bash
#BSUB -J star_pre_wasp[1-30]
#BSUB -o logs/star_pre_wasp_%J_%I.out
#BSUB -e logs/star_pre_wasp_%J_%I.err
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -M 40000
#BSUB -W 48:00

# Purpose: Perform two-pass STAR alignment before WASP filtering.
# Input: sample_list.txt with one sample name per line.
# Expected read layout: 00_trimmedData/<sample>_trimmed_R1.fq.gz and
#                       00_trimmedData/<sample>_trimmed_R2.fq.gz
# Output: MAPQ-filtered, read-grouped BAM files in 03_starMapBeforeWasp/pass2/
# Requirement: The logs/ directory must exist before submitting this LSF job.
# Replace /path/to/ASE with the path to the required Conda environment.
# Replace /path/to/picard.jar with the path to the Picard JAR file.
# STAR, samtools, and Java must be available in the environment or command PATH.

source ~/.bashrc
conda activate /path/to/ASE

sample_list="$PWD/sample_list.txt"
trimmed_reads_dir="$PWD/00_trimmedData"
genome_files_dir="$PWD/01_genomeFiles"
genome_index_dir="$PWD/02_genomeIndex"
star_output_dir="$PWD/03_starMapBeforeWasp"

star_pass1_dir="$star_output_dir/pass1"
star_pass2_dir="$star_output_dir/pass2"
second_pass_index_root="$star_output_dir/passAll"

fasta="${genome_files_dir}/Zm-B73-REFERENCE-NAM-5.0.fa"
picard_jar="/path/to/picard.jar"
threads=10

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

if [[ ! -f "$fasta" ]]; then
    echo "ERROR: Reference FASTA not found: $fasta" >&2
    exit 1
fi

if [[ ! -d "$genome_index_dir" ]]; then
    echo "ERROR: STAR genome index directory not found: $genome_index_dir" >&2
    exit 1
fi

if [[ ! -f "$picard_jar" ]]; then
    echo "ERROR: Picard JAR not found: $picard_jar" >&2
    exit 1
fi

mkdir -p \
    "$star_pass1_dir" \
    "$star_pass2_dir" \
    "$second_pass_index_root"

pass1_prefix="${star_pass1_dir}/${sample}."
pass2_prefix="${star_pass2_dir}/${sample}."
second_pass_index="${second_pass_index_root}/${sample}"

sorted_bam="${star_pass2_dir}/${sample}.sorted.bam"
filtered_bam="${star_pass2_dir}/${sample}_filtered_MQ40.bam"
final_bam="${star_pass2_dir}/${sample}_MQ40_rgadded_sorted.bam"

echo "[${sample}] Step 1: Run first-pass STAR alignment."

STAR \
    --genomeDir "$genome_index_dir" \
    --readFilesCommand zcat \
    --readFilesIn "$r1" "$r2" \
    --outFileNamePrefix "$pass1_prefix" \
    --runThreadN "$threads"

echo "[${sample}] Step 2: Build the sample-specific second-pass index."

mkdir -p "$second_pass_index"

STAR \
    --runMode genomeGenerate \
    --genomeDir "$second_pass_index" \
    --genomeFastaFiles "$fasta" \
    --sjdbFileChrStartEnd "${pass1_prefix}SJ.out.tab" \
    --sjdbOverhang 100 \
    --runThreadN "$threads" \
    --limitGenomeGenerateRAM 40000000000

echo "[${sample}] Step 3: Run second-pass STAR alignment."

STAR \
    --genomeDir "$second_pass_index" \
    --readFilesCommand zcat \
    --readFilesIn "$r1" "$r2" \
    --outFileNamePrefix "$pass2_prefix" \
    --outReadsUnmapped Fastx \
    --outSAMmapqUnique 60 \
    --runThreadN "$threads" \
    --outSAMtype BAM SortedByCoordinate

rm -rf "$second_pass_index"

echo "[${sample}] Step 4: Retain alignments with MAPQ >= 40."

mv "${pass2_prefix}Aligned.sortedByCoord.out.bam" "$sorted_bam"

samtools view \
    -b \
    -q 40 \
    "$sorted_bam" \
    > "$filtered_bam"

samtools index "$filtered_bam"

echo "[${sample}] Step 5: Add read-group information."

java -jar "$picard_jar" AddOrReplaceReadGroups \
    I="$filtered_bam" \
    O="$final_bam" \
    SO=coordinate \
    RGID="$sample" \
    RGLB="$sample" \
    RGPL=Illumina \
    RGPU="$sample" \
    RGSM="$sample"

if [[ ! -s "$final_bam" ]]; then
    echo "ERROR: Final BAM file was not created: $final_bam" >&2
    exit 1
fi

echo "[${sample}] Step 6: Remove intermediate files."

rm -f \
    "${pass1_prefix}Aligned.out.sam" \
    "${pass1_prefix}Log.final.out" \
    "${pass1_prefix}Log.out" \
    "${pass1_prefix}Log.progress.out" \
    "${pass1_prefix}SJ.out.tab" \
    "$sorted_bam" \
    "$filtered_bam" \
    "${filtered_bam}.bai" \
    "${pass2_prefix}SJ.out.tab"

rm -rf "${pass2_prefix}_STARtmp"

echo "[${sample}] Completed: $final_bam"
