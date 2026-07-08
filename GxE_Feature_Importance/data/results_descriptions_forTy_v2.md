### `sv_per_gene.tsv` — Variant counts per gene
**Used in: Figure 2** (`scripts/figures/make_variants_fig.R`)

One row per gene (11,553 genes total). Contains counts of sequence variants found within each gene's promoter region (+/-500bp from TSS).

* `gene` : Gene identifier (B73v5) 
* `chr`, `strand`: Genomic coordinates
* `is_gxe`: TRUE if gene is a GxE gene (padj < 0.1 in DESeq2 GxE test)
* `is_bg`: TRUE if gene is a background gene (no significant GxE or G effect)
* `n_small_indel`: Count of small indels (< 50 bp) in promoter
* `n_large_indel`: Count of large indels (≥ 50 bp) in promoter
* `n_SNP`: Count of B73/Mo17 SNPs in promoter
* `n_sv_total`: Total variant count (SNP + small + large)
* `has_any_sv`: Boolean: does this gene have at least one variant?

### `dap_seq_per_gene.tsv` — DAP-seq allele-specific peak counts per gene
**Used in: Figure 3, panels A–B** (`scripts/figures/make_dap_figure.R`)

One row per gene (11,553 genes). Contains counts of DAP-seq TF binding peaks that are genotype-specific (B73-only or Mo17-only) vs shared, within each gene's promoter. Also includes per-variant-type breakdowns of how many variants of each class overlap peaks of each specificity category (formerly a separate file; now consolidated here).

* `gene`: Gene identifier
* `chr`, `prom_s`, `prom_e`, `strand`: Promoter coordinates
* `is_gxe`: GxE gene flag
* `is_bg`: Background gene flag
* `n_b73_spec`: Count of B73-specific DAP peaks in promoter
* `n_mo17_spec`: Count of Mo17-specific DAP peaks in promoter
* `n_shared`: Count of peaks present in both genotypes
* `n_diff`: Total genotype-specific peaks (B73-specific + Mo17-specific)
* `n_total`: All peaks (specific + shared)
* `frac_diff`: Fraction of peaks that are genotype-specific (`n_diff / n_total`)
* `gene_class`: Categorical: "GxE", "Background", or other grouping

**Small indel (< 50 bp) × peak-type overlap counts:**
* `n_small_indel_b73spec`: Small indels overlapping B73-specific peaks
* `n_small_indel_mo17spec`: Small indels overlapping Mo17-specific peaks
* `n_small_indel_shared`: Small indels overlapping shared peaks
* `n_small_indel_no_peak`: Small indels in promoter with no peak overlap
* `n_small_indel_genotype_spec`: Total small indels overlapping any genotype-specific peak (`n_small_indel_b73spec + n_small_indel_mo17spec`)
* `n_small_indel_total`: All small indels in promoter

**Large indel (≥ 50 bp) × peak-type overlap counts:**
* `n_large_indel_b73spec`: Large indels overlapping B73-specific peaks
* `n_large_indel_mo17spec`: Large indels overlapping Mo17-specific peaks
* `n_large_indel_shared`: Large indels overlapping shared peaks
* `n_large_indel_no_peak`: Large indels in promoter with no peak overlap
* `n_large_indel_genotype_spec`: Total large indels overlapping any genotype-specific peak
* `n_large_indel_total`: All large indels in promoter

**SNP × peak-type overlap counts:**
* `n_snp_b73spec`: SNPs overlapping B73-specific peaks
* `n_snp_mo17spec`: SNPs overlapping Mo17-specific peaks
* `n_snp_shared`: SNPs overlapping shared peaks
* `n_snp_no_peak`: SNPs in promoter with no peak overlap
* `n_snp_genotype_spec`: Total SNPs overlapping any genotype-specific peak
* `n_snp_total`: All SNPs in promoter

### `indel_dap_overlap_per_gene.tsv` — Indels disrupting allele-specific peaks, per gene
**Used in: Figure 3, panels D–E** (`scripts/figures/make_dap_figure.R`)

**NOTE:** This file is a backward-compatible summary; the full per-variant-type breakdown is now also available in `dap_seq_per_gene.tsv`. This information is redundant to that file.
One row per gene (11,553 genes). Counts how many indels in each gene's promoter physically overlap a genotype-specific DAP-seq peak. 

* `gene`: Gene identifier
* `is_gxe`: GxE gene flag
* `is_bg`: Background gene flag
* `n_b73spec`: Indels overlapping B73-specific DAP peaks
* `n_mo17spec`: Indels overlapping Mo17-specific DAP peaks
* `n_shared`: Indels overlapping shared peaks
* `n_no_peak`: Indels in promoter with no peak overlap
* `n_genotype_spec`: Total indels overlapping any genotype-specific peak (`n_b73spec + n_mo17spec`)
* `n_total`: Total indels in promoter
