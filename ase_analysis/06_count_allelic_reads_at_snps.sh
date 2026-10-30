#!/bin/bash
#BSUB -J ase_read_counter[1-30]
#BSUB -o logs/ase_read_counter_%J_%I.out
#BSUB -e logs/ase_read_counter_%J_%I.err
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -M 40000
#BSUB -W 48:00

# Purpose: Count reference- and alternate-allele reads at WGS-derived SNPs.
# Input: sample_list.txt with one sample name per line.
# Expected BAM: 03_starMapBeforeWasp/pass2/<sample>_MQ40_rgadded_sorted.bam
# Expected VCF: 04_vcfFiles/Mo17.diploid.Final.vcf.gz
# Output: ASEReadCounter tables in 05_ASEReadCounter/
# Requirement: The logs/ directory must exist before submitting this LSF job.
# Replace /path/to/ASE with the path to the required Conda environment.
# Replace /path/to/gatk with the path to the locally installed GATK launcher.

source ~/.bashrc
conda activate /path/to/ASE

sample_list="$PWD/sample_list.txt"
genome_files_dir="$PWD/01_genomeFiles"
bam_dir="$PWD/03_starMapBeforeWasp/pass2"
vcf_dir="$PWD/04_vcfFiles"
read_count_dir="$PWD/05_ASEReadCounter"

gatk="/path/to/gatk"

fasta="${genome_files_dir}/Zm-B73-REFERENCE-NAM-5.0.fa"
fasta_index="${fasta}.fai"
sequence_dict="${genome_files_dir}/Zm-B73-REFERENCE-NAM-5.0.dict"
vcf="${vcf_dir}/Mo17.diploid.Final.vcf.gz"

sample=$(sed -n "${LSB_JOBINDEX}p" "$sample_list")

if [[ -z "$sample" ]]; then
    echo "ERROR: No sample found for LSF array index ${LSB_JOBINDEX}." >&2
    exit 1
fi

bam="${bam_dir}/${sample}_MQ40_rgadded_sorted.bam"
bam_index="${bam}.bai"
output_table="${read_count_dir}/${sample}_ASE_read_counts.table"

if [[ ! -f "$bam" ]]; then
    echo "ERROR: Input BAM not found: $bam" >&2
    exit 1
fi

if [[ ! -f "$fasta" ]]; then
    echo "ERROR: Reference FASTA not found: $fasta" >&2
    exit 1
fi

if [[ ! -f "$fasta_index" ]]; then
    echo "ERROR: Reference FASTA index not found: $fasta_index" >&2
    exit 1
fi

if [[ ! -f "$sequence_dict" ]]; then
    echo "ERROR: Reference sequence dictionary not found: $sequence_dict" >&2
    exit 1
fi

if [[ ! -f "$vcf" ]]; then
    echo "ERROR: Input VCF not found: $vcf" >&2
    exit 1
fi

if [[ ! -x "$gatk" ]]; then
    echo "ERROR: GATK launcher not found or not executable: $gatk" >&2
    exit 1
fi

mkdir -p "$read_count_dir"

if [[ ! -f "$bam_index" ]]; then
    samtools index "$bam"
fi

"$gatk" ASEReadCounter \
    -R "$fasta" \
    -I "$bam" \
    -V "$vcf" \
    -DF NotDuplicateReadFilter \
    --output-format TABLE \
    -O "$output_table"

if [[ ! -s "$output_table" ]]; then
    echo "ERROR: ASEReadCounter output was not created: $output_table" >&2
    exit 1
fi

echo "[${sample}] Completed: $output_table"
