# GxE-ASE Analysis Rerun Workflow

All commands run from the project root:
```
cd /path/to/gxe_ase
```

---

## Step 0.1 - Download publicly available resources and prepare genome comparisons

Get B73 reference genome and annotation:
```bash
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -O data/annotation.gff
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz -P data/
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz.fai -P data/
```

Get genotype specific and shared DAPseq peaks, and motifs, from Galli et al 2025 https://doi.org/10.1038/s41477-025-02007-8
```bash
mkdir -p data/epigenome/dap_seq/normalized_specific_and_shared_peaks/
wget https://zenodo.org/records/14991916/files/normalized_specific_and_shared_peaks.zip?download=1 && \
unzip "normalized_specific_and_shared_peaks.zip?download=1" -d data/epigenome/dap_seq/normalized_specific_and_shared_peaks/ && \
rm "normalized_specific_and_shared_peaks.zip?download=1"

wget https://zenodo.org/records/14991916/files/motifs.zip?download=1 && \
unzip "motifs.zip?download=1" -d data/epigenome/dap_seq/ && \
rm "motifs.zip?download=1"
```

To repeat the phenotypic GxE mapping and overlap with GxE genes, download these files into ./data/:
 * HapMap3.1 GBS SNPs from "https://cornell.box.com/s/o7wtp1ewuqlw3dalr1920lungxnomnrg". The file name should be ZeaGBSv27_publicSamples_raw_AGPv4-181023.vcf.gz
 * Sample info from "https://cornell.app.box.com/s/v5rsmdtdg0g5ecjtawfonvavuzuffp6z"

---

## Step 0.2 — Generate ASE read counts *(ase_analysis/)*

Ten LSF cluster scripts (`bsub`, hard-coded paths — not run like the numbered scripts below), adapted from the ASEPipeline in Hu et al. (2022), that take raw paired-end RNA-seq reads through to per-gene allele-specific counts. This is what ultimately produces `data/counts_for_DESeq2.txt` and `data/meta_for_DESeq2.txt`, the inputs to Step 1.

1. `01_run_fastqc_raw_reads.sh` — FastQC on raw reads
2. `02_trim_reads_with_fastp.sh` — trim/filter reads with fastp
3. `03_run_fastqc_cleaned_reads.sh` — FastQC on trimmed reads
4. `04_prepare_reference_indexes.sh` — build the STAR index, FASTA index, and sequence dictionary
5. `05_star_pre_wasp_mapping.sh` — two-pass STAR mapping (MQ ≥ 40), add read groups
6. `06_count_allelic_reads_at_snps.sh` — GATK ASEReadCounter at WGS-derived SNP positions
7. `07_create_sample_specific_wasp_vcf.sh` — build a sample-specific phased VCF for WASP
8. `08_star_wasp_mapping.sh` — two-pass STAR mapping with WASP remapping-bias filtering
9. `09_separate_parental_reads.sh` — split WASP-passing reads into Parent1/Parent2 BAMs
10. `10_generate_total_and_parental_gene_counts.sh` — featureCounts total + parental gene counts

See `ase_analysis/README.md` for more detail.

---

## Step 0.3 — Call B73/Mo17 variants *(variant_features/)*

Five scripts that align the B73 and Mo17 genomes and call variants between them amd produces the three VCFs read throughout Steps 3–5.

1. `01_prepare_anchorwave_inputs.sh` — extract B73 CDS anchors, map to B73 and Mo17 with minimap2
2. `02_run_anchorwave_b73_mo17_alignment.sh` — whole-genome AnchorWave alignment (MAF)
3. `03_convert_anchorwave_maf_to_paf.sh` — convert the MAF alignment to PAF for SyRI
4. `04_run_syri_b73_mo17.sh` — call structural differences between B73 and Mo17 with SyRI
5. `05_parse_syri_variants.sh` — split the SyRI VCF into `..._SNPs.vcf`, `..._INDELs_less50bp.vcf`, and `..._INDELs_more50bp.vcf`

Also LSF cluster scripts (except script 5, which is plain bash) with hard-coded paths.

---

## Step 1 — DESeq2 GxE-ASE test *(R)* — `scripts/1_deseq_analysis.R`

Fits the DESeq2 allele-specific-expression model and defines the gene sets everything downstream uses: GxE-ASE genes (padj < 0.1 on the environment x allele interaction), a pre-filtered background set, and G genes (constitutive ASE in the same direction in both environments).

Reads `data/counts_for_DESeq2.txt` and `data/meta_for_DESeq2.txt`. Produces `data/GxE_gene_IDs.txt`, `data/background_gene_IDs.txt`, `data/G_gene_IDs.txt`, `data/DEG_GxE_results.txt`, `data/GxE_allele_specific_test_results.txt`, `data/DEG_Env_results.txt`, and `data/DESeq2_results.rdata`.

```bash
Rscript scripts/1_deseq_analysis.R
```

---

## Step 2 — Sample-size / power check *(R)* — `scripts/2_subsample_test_GxE.R`

Repeatedly subsamples the DESeq2 dataset down to 3–14 replicates per environment and reruns the model, showing how the number of detected GxE genes scales with sample size. Can be run any time after Step 1.

Reads `data/DESeq2_results.rdata` (Step 1). Produces `results/n_siginificant_subsampling.txt`.

```bash
Rscript scripts/2_subsample_test_GxE.R
```

---

## Step 3 — DAP-seq peak analysis *(Python)* — `scripts/3_dap_seq_analysis.py`

Counts genotype-specific vs. shared DAP-seq peaks in each gene's promoter (±500 bp of TSS), plus how many indels/SNPs overlap those peaks.

Reads `data/annotation.gff`, the Step 1 gene sets, `data/tf_annotation.tsv`, the DAP-seq peak files, and the three Mo17/B73 VCFs. Produces `results/dap_seq_per_gene.tsv`, `results/dap_seq_per_tf.tsv`, `results/indel_dap_overlap_per_gene.tsv`.

```bash
python3 scripts/3_dap_seq_analysis.py
```

---

## Step 4 — SV enrichment *(Python)* — `scripts/4_sv_enrichment.py`

Tests whether GxE-ASE genes carry more B73/Mo17 SNPs and indels in their promoters than background genes (Fisher's exact + Mann-Whitney), including a distance-from-TSS breakdown.

Reads `results/dap_seq_per_gene.tsv` (Step 3), the Step 1 gene sets, and the three Mo17/B73 VCFs. Produces `results/sv_enrichment_results.tsv`, `results/sv_per_gene.tsv`, `results/sv_distance_results.tsv`, and two figures.

```bash
python3 scripts/4_sv_enrichment.py
```

---

## Step 5 — DAP-seq peak / variant proximity *(Python)* — `scripts/5_peak_variant_proximity.py`

Tests whether genotype-specific DAP-seq peaks sit closer to B73/Mo17 variants than shared peaks do, split by GxE vs. background genes. Prints results to stdout — writes no files.

Reads `results/dap_seq_per_gene.tsv` (Step 3), the Step 1 gene sets, the three Mo17/B73 VCFs, and the DAP-seq peak files.

```bash
python3 scripts/5_peak_variant_proximity.py
```

---

## Step 6 — Constitutive (G) vs. GxE enrichment *(Python)* — `scripts/6_g_gene_enrichment.py`

Splits genes into GxE-only, G-only, GxE+G, and Neither, and re-checks the Step 3/4 enrichment signals within each group — a robustness check that the signals reflect GxE specifically, not general B73/Mo17 divergence.

Reads the Step 1 gene sets plus `results/sv_per_gene.tsv` (Step 4), `results/dap_seq_per_gene.tsv`, and `results/indel_dap_overlap_per_gene.tsv` (Step 3). Produces `results/g_gene_enrichment.tsv` and a figure.

```bash
python3 scripts/6_g_gene_enrichment.py
```

---

## Step 7 — Filter GBS genotypes to IBM samples *(bash)* — `scripts/7_filter_ZeaGBS_to_IBM.sh`

Filters the full public HapMap v3 GBS VCF (from Step 0) down to the IBM population + Mo17 samples using the sample metadata spreadsheet.

Reads `data/AllZeaGBSv2.7_publicSamples_metadata20140411.xlsx` and `data/ZeaGBSv27_publicSamples_raw_AGPv4-181023.vcf.gz` (Step 0). Produces `data/ZeaGBSv27_IBM_raw_AGPv4.vcf.gz`.

```bash
bash scripts/7_filter_ZeaGBS_to_IBM.sh
```

---

## Step 8 — Build IBM recombination bins *(R)* — `scripts/8_format_gbs_genotypes.R`

Builds recombination bins for the IBM (B73 x Mo17) RIL population from the filtered GBS genotypes — lifts AGPv4 coordinates to NAM v5 and imputes with a qtl2 HMM.

Reads `data/ZeaGBSv27_IBM_raw_AGPv4.vcf.gz` (Step 7) and `data/AllZeaGBSv2.7_publicSamples_metadata20140411.xlsx`. Produces `data/IBM_recomb_bins_fromGBS.tsv` and `data/qtl2_input/`.

```bash
Rscript scripts/8_format_gbs_genotypes.R
```

---

## Step 9 — IBM phenotypic GxE QTL scan *(R)* — `scripts/9_map_gxe_qtl.R`

Two-stage GxE QTL scan of 19 NAM phenotypic traits across the IBM population's recombination bins (block-adjust, then an incremental F-test for the genotype x environment interaction).

Reads `data/IBM_recomb_bins_fromGBS.tsv` (Step 8), `data/IBM_Name_M00_Z017_decoder.txt`, and `data/NAM_all_traits.txt`. Produces `results/IBM_GxE_results.tsv`.

```bash
Rscript scripts/9_map_gxe_qtl.R
```

---

## Step 10 — GxE-ASE x QTL enrichment *(Python)* — `scripts/10_gxe_qtl_enrichment.py`

Tests whether GxE-ASE genes co-localize with the phenotypic GxE QTL bins from Step 9 (Stouffer-combined p-values across traits, binary + continuous enrichment, chromosome-stratified permutation).

Reads `results/IBM_GxE_results.tsv` (Step 9), `data/IBM_recomb_bins_fromGBS.tsv` (Step 8), `data/annotation.gff`, and the Step 1 gene sets/results. Produces `results/gbs_bin_stouffer_stats.tsv`, `results/gbs_gene_bin_assignments.tsv`, `results/gbs_enrichment_results.tsv`, `results/gbs_enrichment_by_trait.tsv`, and a figure.

```bash
python3 scripts/10_gxe_qtl_enrichment.py
```

---

## Manuscript figures

Four R scripts in `scripts/figures/` turn the tables above into the paper figures, each reading one or two result files:

- `make_GxE_and_power_figure.R` — GxE scatter + subsampling power (Steps 1, 2) → `figures/1_GxE_genes_and_subsample.{png,pdf}`
- `make_variants_fig.R` — SV enrichment (Step 4) → `figures/2_variants_fig.{pdf,png}`
- `make_dap_figure.R` — DAP-seq peak and indel/SNP-disruption enrichment (Step 3) → `figures/3_dap_fig.{pdf,png}`
- `make_G_vs_GxE_figure.R` — G vs. GxE enrichment comparison (Step 6) → `figures/4_compare_G_GxE.{pdf,png}`

```bash
for f in scripts/figures/*.R; do Rscript "$f"; done
```

---

## Modeling feature importance in predicting transcriptional GxE

`GxE_Feature_Importance/scripts/Feature_Importance.Rmd` (helpers in `Functions_Feature_Importance.R`) fits LASSO, SVM, random forest, and logistic models to rank which SV/DAP-seq features best predict GxE-ASE status. It reads its own snapshot of the per-gene tables (`GxE_Feature_Importance/data/sv_per_gene.tsv`, `dap_seq_per_gene.tsv`) rather than the live `results/` files, so re-copy those if Steps 3–4 are rerun. Outputs are written outside this repo.

---

## Suggested run order

`ase_analysis/` → Step 1, and `variant_features/` → Steps 3–5, are separate upstream pipelines (cluster jobs, not part of the numbered sequence) that must complete first. From there, the script numbers match dependency order, so running 1 through 10 in sequence works:

```
ase_analysis/      → Step 1 (produces data/counts_for_DESeq2.txt, data/meta_for_DESeq2.txt)
variant_features/  → Steps 3, 4, 5 (produces the three Mo17/B73 VCFs)

1 → 2                 (2 only needs 1)
1 → 3 → 4 → 5         (DAP-seq / SV enrichment track; 5 needs 1 and 3)
      └→ 6            (needs 1, 3, and 4)
7 → 8 → 9 → 10        (phenotypic QTL track; 10 also needs 1 and data/annotation.gff)
figures                (need results/ and data/ files from the tracks above)
```
