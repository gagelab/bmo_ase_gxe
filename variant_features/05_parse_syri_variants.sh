#!/bin/bash

# Purpose: Parse the B73–Mo17 SyRI VCF into SNP, small-indel, and large-indel files.
# Input: SyRI VCF in 02_output/04_syri_output/
# Outputs:
#   SNPs
#   Indels smaller than 50 bp
#   Indels 50 bp or larger

syri_output_dir="$PWD/02_output/04_syri_output"

input_vcf="${syri_output_dir}/Mo17_toB73v5_paf_syri.vcf"
filtered_vcf="${syri_output_dir}/Mo17_toB73v5_paf_syri_noStartPOS0.vcf"

snp_vcf="${syri_output_dir}/Mo17_toB73v5_paf_syri_noStartPOS0_SNPs.vcf"
small_indel_vcf="${syri_output_dir}/Mo17_toB73v5_paf_syri_noStartPOS0_INDELs_less50bp.vcf"
large_indel_vcf="${syri_output_dir}/Mo17_toB73v5_paf_syri_noStartPOS0_INDELs_more50bp.vcf"

if [[ ! -f "$input_vcf" ]]; then
    echo "ERROR: SyRI VCF not found: $input_vcf" >&2
    exit 1
fi

# Remove records with a start position of zero while retaining the VCF header.
awk 'BEGIN { OFS="\t" } /^#/ || $2 > 0 { print }' \
    "$input_vcf" \
    > "$filtered_vcf"

# Parse variants and preserve the original final ordering:
# VCF header, followed by deletions, then insertions.
awk \
    -v snp_vcf="$snp_vcf" \
    -v small_indel_vcf="$small_indel_vcf" \
    -v large_indel_vcf="$large_indel_vcf" \
    '
    /^#/ {
        headers[++header_count] = $0
        next
    }

    $3 ~ /^SNP/ {
        snps[++snp_count] = $0
        next
    }

    $3 ~ /^DEL/ {
        if (length($4) < 50) {
            small_deletions[++small_del_count] = $0
        } else {
            large_deletions[++large_del_count] = $0
        }
        next
    }

    $3 ~ /^INS/ {
        if (length($5) < 50) {
            small_insertions[++small_ins_count] = $0
        } else {
            large_insertions[++large_ins_count] = $0
        }
    }

    END {
        for (i = 1; i <= header_count; i++) {
            print headers[i] > snp_vcf
            print headers[i] > small_indel_vcf
            print headers[i] > large_indel_vcf
        }

        for (i = 1; i <= snp_count; i++) {
            print snps[i] > snp_vcf
        }

        for (i = 1; i <= small_del_count; i++) {
            print small_deletions[i] > small_indel_vcf
        }
        for (i = 1; i <= small_ins_count; i++) {
            print small_insertions[i] > small_indel_vcf
        }

        for (i = 1; i <= large_del_count; i++) {
            print large_deletions[i] > large_indel_vcf
        }
        for (i = 1; i <= large_ins_count; i++) {
            print large_insertions[i] > large_indel_vcf
        }
    }
    ' "$filtered_vcf"

echo "SyRI variant parsing completed:"
echo "  SNPs: $snp_vcf"
echo "  Small indels: $small_indel_vcf"
echo "  Large indels: $large_indel_vcf"
