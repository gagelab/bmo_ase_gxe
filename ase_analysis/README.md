# Allele-Specific Expression Analysis

These scripts were adapted and modified from the ASEPipeline developed for the Hu et al. (2022) study, *Allele-specific Expression Reveals Multiple Paths to Highland Adaptation in Maize*. The original scripts are available in the [Maize Highland Adaptation allele-specific expression repository](https://github.com/hh622/Maize_Highland_Adaptation_allele_specific_expression).

## Scripts

1. **`01_run_fastqc_raw_reads.sh`**  
   Runs FastQC on the raw paired-end RNA-seq reads before trimming and cleaning.

2. **`02_trim_reads_with_fastp.sh`**  
   Trims and filters the raw paired-end reads with fastp and generates HTML and JSON quality reports.

3. **`03_run_fastqc_cleaned_reads.sh`**  
   Runs FastQC on the paired-end reads after fastp trimming and filtering.

4. **`04_prepare_reference_indexes.sh`**  
   Creates the STAR genome index, FASTA index, and sequence dictionary required for read mapping and GATK ASEReadCounter.

5. **`05_star_pre_wasp_mapping.sh`**  
   Performs two-pass STAR mapping, retains alignments with a mapping quality of at least 40, and adds read-group information before WASP filtering.

6. **`06_count_allelic_reads_at_snps.sh`**  
   Uses GATK ASEReadCounter to count reference- and alternate-allele reads at WGS-derived SNP positions for each sample.

7. **`07_create_sample_specific_wasp_vcf.sh`**  
   Selects informative SNPs based on allelic read depth and balance and creates a sample-specific phased VCF for STAR-WASP mapping.

8. **`08_star_wasp_mapping.sh`**  
   Performs two-pass STAR mapping with sample-specific variants and retains reads that pass WASP remapping-bias filtering.

9. **`09_separate_parental_reads.sh`**  
   Separates WASP-passing reads into Parent1- and Parent2-associated BAM files using the STAR allele tags.

10. **`10_generate_total_and_parental_gene_counts.sh`**  
    Uses featureCounts to generate total gene-count tables and separate Parent1 and Parent2 gene-count columns.
