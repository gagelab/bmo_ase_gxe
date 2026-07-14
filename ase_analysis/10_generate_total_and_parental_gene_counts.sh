#!/bin/bash
#BSUB -J featurecounts
#BSUB -o logs/featurecounts_%J.out
#BSUB -e logs/featurecounts_%J.err
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -M 40000
#BSUB -W 48:00
#BSUB -q gage

# Purpose: Generate total and Parent1/Parent2 gene-count tables with featureCounts.
# Inputs:
#   03_starMapBeforeWasp/pass2/*_MQ40_rgadded_sorted.bam
#   07_starMappingWasp/pass2/*WASPfilterPASS_RG_sort.vA*.bam
#   Reference GTF in 01_genomeFiles/
# Outputs:
#   Total gene counts from pre-WASP BAM files
#   Parent1 and Parent2 gene counts from parent-separated BAM files
# Requirement: The logs/ directory must exist before submitting this LSF job.
# Replace /path/to/ASE with the path to the required Conda environment.
# featureCounts must be available in the environment or command PATH.

source ~/.bashrc
conda activate /path/to/ASE

pre_wasp_bam_dir="$PWD/03_starMapBeforeWasp/pass2"
parental_bam_dir="$PWD/07_starMappingWasp/pass2"
genome_files_dir="$PWD/01_genomeFiles"
count_dir="$PWD/08_featureCount"

gtf="${genome_files_dir}/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gtf"
total_count_output="${count_dir}/total_featurecounts_0.results"
parental_count_output="${count_dir}/allele_featurecounts_0.results"
threads=10

if [[ ! -f "$gtf" ]]; then
    echo "ERROR: Reference GTF not found: $gtf" >&2
    exit 1
fi

mkdir -p "$count_dir"

shopt -s nullglob

total_bams=(
    "$pre_wasp_bam_dir"/*_MQ40_rgadded_sorted.bam
)

parent1_bams=(
    "$parental_bam_dir"/*WASPfilterPASS_RG_sort.vA1.bam
)

parent2_bams=(
    "$parental_bam_dir"/*WASPfilterPASS_RG_sort.vA2.bam
)

if [[ ${#total_bams[@]} -eq 0 ]]; then
    echo "ERROR: No pre-WASP BAM files found in $pre_wasp_bam_dir" >&2
    exit 1
fi

if [[ ${#parent1_bams[@]} -eq 0 ]]; then
    echo "ERROR: No Parent1 BAM files found in $parental_bam_dir" >&2
    exit 1
fi

if [[ ${#parent2_bams[@]} -eq 0 ]]; then
    echo "ERROR: No Parent2 BAM files found in $parental_bam_dir" >&2
    exit 1
fi

# Preserve the original wildcard ordering so that vA1 and vA2 BAM files
# for each sample occur next to one another in the featureCounts output.
parental_bams=(
    "$parental_bam_dir"/*WASPfilterPASS_RG_sort.vA*.bam
)

echo "Step 1: Count total reads assigned to genes."

featureCounts \
    -T "$threads" \
    -s 0 \
    -p \
    -g gene_id \
    -t exon \
    -a "$gtf" \
    -o "$total_count_output" \
    "${total_bams[@]}"

if [[ ! -s "$total_count_output" ]]; then
    echo "ERROR: Total gene-count output was not created: $total_count_output" >&2
    exit 1
fi

echo "Step 2: Count Parent1- and Parent2-associated reads assigned to genes."

featureCounts \
    -T "$threads" \
    -s 0 \
    -p \
    -g gene_id \
    -t exon \
    -a "$gtf" \
    -o "$parental_count_output" \
    "${parental_bams[@]}"

if [[ ! -s "$parental_count_output" ]]; then
    echo "ERROR: Parental gene-count output was not created: $parental_count_output" >&2
    exit 1
fi

echo "Completed:"
echo "  Total gene counts: $total_count_output"
echo "  Parent1/Parent2 gene counts: $parental_count_output"
