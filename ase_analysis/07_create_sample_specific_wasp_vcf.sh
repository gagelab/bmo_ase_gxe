#!/bin/bash
#BSUB -J wasp_vcf[1-30]
#BSUB -o logs/wasp_vcf_%J_%I.out
#BSUB -e logs/wasp_vcf_%J_%I.err
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -M 40000
#BSUB -W 48:00

# Purpose: Select informative SNPs and create a sample-specific VCF for STAR-WASP.
# Input: sample_list.txt with one sample name per line.
# Expected ASE table: 05_ASEReadCounter/<sample>_ASE_read_counts.table
# Expected VCF: 04_vcfFiles/Mo17.diploid.Final.vcf.gz
# Outputs: Sample-specific VCFs in 06_SampleSpecificVCF/ and filtered
#          ASEReadCounter tables in 05_ASEReadCounter/
# Requirement: The logs/ directory must exist before submitting this LSF job.
# Replace /path/to/ASE with the path to the required Conda environment.

source ~/.bashrc
conda activate /path/to/ASE

sample_list="$PWD/sample_list.txt"
vcf_dir="$PWD/04_vcfFiles"
read_count_dir="$PWD/05_ASEReadCounter"
sample_vcf_dir="$PWD/06_SampleSpecificVCF"

input_vcf="${vcf_dir}/Mo17.diploid.Final.vcf.gz"
minimum_depth=10
maximum_log2_ratio=2

sample=$(sed -n "${LSB_JOBINDEX}p" "$sample_list")

if [[ -z "$sample" ]]; then
    echo "ERROR: No sample found for LSF array index ${LSB_JOBINDEX}." >&2
    exit 1
fi

ase_table="${read_count_dir}/${sample}_ASE_read_counts.table"
selected_sites="${read_count_dir}/${sample}_chr_Position_Selected_Thr10_ASE2.txt"
plotting_table="${read_count_dir}/${sample}_ASE_read_counts_2.table"
sample_vcf="${sample_vcf_dir}/${sample}_phasedGT.vcf"

if [[ ! -f "$ase_table" ]]; then
    echo "ERROR: ASEReadCounter table not found: $ase_table" >&2
    exit 1
fi

if [[ ! -f "$input_vcf" ]]; then
    echo "ERROR: Input VCF not found: $input_vcf" >&2
    exit 1
fi

mkdir -p "$sample_vcf_dir"

# Retain SNPs with both alleles observed, total depth >= 10,
# and an absolute log2 reference/alternate count ratio <= 2.
tail -n +2 "$ase_table" |
    awk -v minimum_depth="$minimum_depth" \
        'BEGIN { OFS="\t" }
         $6 > 0 && $7 > 0 && $8 >= minimum_depth { print }' |
    awk -v maximum_log2_ratio="$maximum_log2_ratio" \
        'BEGIN { OFS="\t" }
         {
             log2_ratio = log($6 / $7) / log(2)
             if (log2_ratio <= maximum_log2_ratio &&
                 log2_ratio >= -maximum_log2_ratio) {
                 print $1, $2
             }
         }' |
    awk 'BEGIN { OFS="\t" } { print $1 "." $2 }' \
    > "$selected_sites"

if [[ ! -s "$selected_sites" ]]; then
    echo "ERROR: No SNPs passed the selection criteria for sample ${sample}." >&2
    exit 1
fi

# Subset the WGS-derived VCF to SNPs selected for this RNA-seq sample.
zcat "$input_vcf" |
    awk 'BEGIN { OFS="\t" } { $11=$1 "." $2; print }' |
    grep -F -w -f "$selected_sites" |
    cut -f1-10 |
    sed '1i #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPhasedGenotype' \
    > "$sample_vcf"

if [[ ! -s "$sample_vcf" ]]; then
    echo "ERROR: Sample-specific VCF was not created: $sample_vcf" >&2
    exit 1
fi

# Retain sites with reads supporting both alleles for ASE-ratio plotting.
awk 'BEGIN { OFS="\t" } $6 > 0 && $7 > 0 { print }' \
    "$ase_table" \
    > "$plotting_table"

echo "[${sample}] Completed: $sample_vcf"
