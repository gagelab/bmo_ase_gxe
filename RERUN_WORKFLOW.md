# GxE-ASE Analysis Rerun Workflow

TODO:
 - Get B73 Mo17 vcf code from Sontosh
 - Results of script 2 (DAP analysis) don't agree with expected (printed at end)
 - Remove mechanistic section 6? Null results, low sample size(7 DEG TFs)
 - Remove marginal sensitivity section 7. Tests results with pvalue near .05, against some other parameters

All commands run from the project root:
```
cd /path/to/gxe_ase
```

---

## Step 0 - Download publicly available resources and prepare genome comparisons

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

Compute SNPs, small Indels (<50bp) and large Indels (>50bp) between B73 and Mo17:
TODO: Fill this in

---

## Step 1 — DESeq2 GxE-ASE test  *(R)*

Pre-requisites: `data/counts_for_DESeq2.txt`, `data/meta_for_DESeq2.txt`

Produces: `data/DEG_GxE_results.txt`, `data/GxE_allele_specific_test_results.txt`,
`data/GxE_gene_IDs.txt`, `data/background_gene_IDs.txt`,
`data/DEG_Env_results.txt`, `data/DESeq2_results.rdata`

```r
Rscript scripts/1_deseq_analysis.R
```

**Expected:** 248 GxE genes (padj < 0.1), 11,389 background genes.
Check `data/GxE_gene_IDs.txt`: `wc -l data/GxE_gene_IDs.txt` → 248

---

## Step 2 — DAP-seq analysis *(~5–10 min)*

Generates all DAP-seq derived result files directly from the raw supplementary
peak files from the O'Malley et al. DAP-seq dataset, the B73v5 gene annotation,
and the small-indel VCF. No intermediate files are required beyond the gene lists
from Step 1.

Pre-requisites (computed from previous steps):
- `data/GxE_gene_IDs.txt`, `data/background_gene_IDs.txt` — from Step 1
- `data/annotation.gff` — B73v5 gene annotation (TSS coordinates), from Step 0
- `data/Mo17_toB73v5_paf_syri_noStartPOS0_INDELs_less50bp.vcf` — small indels < 50 bp, from Step 0
- `data/epigenome/dap_seq/normalized_specific_and_shared_peaks/B73v5_MATCH_Mo17-B73v5_specific_B73v5_shared/*_withcoords.tsv` — B73-specific and shared peaks (B73v5 coordinates in cols 0–2), from Step 0
- `data/epigenome/dap_seq/normalized_specific_and_shared_peaks/Mo17_MATCH_B73v5-Mo17_specific_Mo17_shared/*_withcoords.tsv` — Mo17-specific peaks (B73 locus coordinates in cols 10–12, assigned by bedtools closest), from Step 0
Pre-requisites (in git repo):
- `data/tf_annotation.tsv` — stable TF → gene ID lookup (58 expressed TFs)

Produces: `results/dap_seq_per_gene.tsv`, `results/dap_seq_per_tf.tsv`,
`results/indel_dap_overlap_per_gene.tsv`

```bash
python3 scripts/2b_dap_seq_analysis.py
```

**Expected:**
- `dap_seq_per_gene.tsv`: 11,637 genes
- DAP genotype-specific enrichment: GxE ~46.0% vs BG ~37.2%, OR ≈ 1.44, p ≈ 0.005
- Indel×DAP enrichment: GxE ~21.8% vs BG ~14.6%, OR ≈ 1.62, p ≈ 0.003

---

## Step 3 — SV (small indel) enrichment *(~20 min — distance figure is slow)*

Reads: `data/annotation.gff`, `data/GxE_gene_IDs.txt`, `data/background_gene_IDs.txt`,
`data/Mo17_toB73v5_paf_syri_noStartPOS0_INDELs_less50bp.vcf`,
`data/Mo17_toB73v5_paf_syri_noStartPOS0_INDELs_more50bp.vcf`,
`data/Mo17_toB73v5_paf_syri_noStartPOS0_SNPs.vcf`

Produces: `results/sv_enrichment_results.tsv`, `results/sv_per_gene.tsv`,
`figures/sv_enrichment_figure.png`

> ⚠️ The distance-from-TSS section of this script is slow (~20 min).
> Core enrichment results are available quickly; the full run is needed for
> the distance figure only.

```bash
python3 scripts/3_sv_enrichment.py
```

**Expected core results:**
- Small indel: GxE 86.3% vs BG 80.0%, OR ≈ 1.57, Fisher p ≈ 0.015, MWU p ≈ 6.8e-7

---

## Step 4 — Threshold robustness

Reads: `data/GxE_gene_IDs.txt`, `data/background_gene_IDs.txt`,
`results/dap_seq_per_gene.tsv` (Step 2), `results/sv_per_gene.tsv` (Step 3),
`results/indel_dap_overlap_per_gene.tsv` (Step 2)

Produces: `results/threshold_sensitivity_results.tsv`,
`figures/threshold_sensitivity_figure.png`

```bash
python3 scripts/4_threshold_sensitivity.py
```

**Expected:** 6 thresholds × 3 tests; all ORs > 1.1.
Canonical threshold (padj < 0.1): DAP OR ≈ 1.44, SV OR ≈ 1.57, Indel×DAP OR ≈ 1.62.

---

## Step 5 — Indel PWM scoring  *(slow — ~15 min)*

Reads:
- `data/GxE_gene_IDs.txt`, `data/background_gene_IDs.txt` — from Step 1
- `data/GxE_allele_specific_test_results.txt` — from Step 1
- `data/annotation.gff` — B73v5 gene annotation
- `data/Mo17_toB73v5_paf_syri_noStartPOS0_INDELs_less50bp.vcf` — small indels
- `data/epigenome/dap_seq/normalized_specific_and_shared_peaks/.../*_withcoords.tsv` — peak files (same as Step 2)
- `data/epigenome/dap_seq/gem02_rep_memechip00_onefile/gem02_rep_memechip00_m1.txt` — MEME motif file
- `data/Zm-B73-REFERENCE-NAM-5.0.fa.gz` — B73 reference FASTA (indexed)

Produces: `results/indel_pwm_scored_records.tsv`, `results/indel_pwm_stats.tsv`,
`figures/indel_pwm_scoring_figure.png`

```bash
python3 scripts/5_indel_pwm_scoring.py
```

After this runs, re-label the records with the current gene lists:

```bash
python3 - <<'EOF'
import pandas as pd
gxe = set(open("data/GxE_gene_IDs.txt").read().split())
bg  = set(open("data/background_gene_IDs.txt").read().split())
pwm = pd.read_csv("results/indel_pwm_scored_records.tsv", sep="\t")
pwm["is_gxe"] = pwm["gene"].isin(gxe)
pwm["is_bg"]  = pwm["gene"].isin(bg)
pwm.to_csv("results/indel_pwm_scored_records.tsv", sep="\t", index=False)
b73 = pwm[pwm["peak_geno"]=="B73"]
pct = 100*(b73["delta"]>0).mean()
print(f"B73-spec Δ>0: {pct:.1f}%  (expect ~60.7%)")
EOF
```

**Expected:** ~60.7% of B73-specific peak records have Δ > 0 (binom p ≈ 1.5e-36).
GxE-specific rate ≈ 61.6% — not significantly enriched vs background.

---

## Step 6 — TF mechanism analysis

Reads:
- `data/GxE_gene_IDs.txt`, `data/background_gene_IDs.txt` — from Step 1
- `data/GxE_allele_specific_test_results.txt`, `data/DEG_Env_results.txt` — from Step 1
- `data/tf_annotation.tsv` — stable TF → gene ID lookup
- `data/annotation.gff` — B73v5 gene annotation
- `results/dap_seq_per_tf.tsv` — **from Step 2** (TF DEG status is derived from this)
- `data/epigenome/dap_seq/normalized_specific_and_shared_peaks/.../*_withcoords.tsv` — peak files (same as Step 2)

Produces: `results/mechanism_per_gene.tsv`, `results/mechanism_per_tf.tsv`,
`figures/mechanism_figure.png`

```bash
python3 scripts/6_tf_gxe_mechanism.py
```

---

## Step 7 — Marginal sensitivity analysis

Reads: `data/DEG_GxE_results.txt` (Step 1), `results/gene_positional_features.tsv`,
`results/mechanism_per_gene.tsv` (Step 6), `data/pangene_table.tsv`

Produces: `results/marginal_sensitivity_*.tsv`,
`figures/marginal_sensitivity_figure.png`

```bash
python3 scripts/7_marginal_sensitivity.py
```

---

## Step 8 — GO enrichment

Reads: `data/GxE_gene_IDs.txt`, `data/background_gene_IDs.txt` (Step 1),
`data/annotation.gff`, GO annotation file

Produces: `results/go_enrichment_full.tsv`, `results/go_proxy_results.tsv`,
`figures/go_proxy_enrichment_figure.png`

```bash
python3 scripts/8_go_enrichment.py
```

---

## Step 9 — Manuscript figures  *(fast — reads pre-computed results)*

Reads: `results/sv_per_gene.tsv` (Step 3), `results/dap_seq_per_gene.tsv` (Step 2),
`results/indel_dap_overlap_per_gene.tsv` (Step 2),
`results/indel_pwm_stats.tsv` (Step 5),
`results/threshold_sensitivity_results.tsv` (Step 4)

Produces: `figures/figure1_gxe_sv_enrichment.{pdf,png}`,
`figures/figure2_indel_dap_pwm.{pdf,png}`,
`figures/figure3_threshold_robustness.{pdf,png}`

```bash
python3 scripts/9_manuscript_figures.py
```

**Quick verification of key figure values:**
```bash
python3 - <<'EOF'
import pandas as pd
from scipy.stats import fisher_exact

gxe = set(open("data/GxE_gene_IDs.txt").read().split())
bg  = set(open("data/background_gene_IDs.txt").read().split())
sv  = pd.read_csv("results/sv_per_gene.tsv", sep="\t")
g = sv[sv["is_gxe"]]; b = sv[sv["is_bg"]]
OR, p = fisher_exact([[g["has_sv"].sum(), len(g)-g["has_sv"].sum()],
                      [b["has_sv"].sum(), len(b)-b["has_sv"].sum()]])
print(f"GxE n={len(g)}, BG n={len(b)}")
print(f"SV enrichment: OR={OR:.3f}, Fisher p={p:.4g}")
print(f"GxE: {100*g['has_sv'].mean():.1f}%  BG: {100*b['has_sv'].mean():.1f}%")
print("Expected: OR≈1.570, p≈0.015, 86.3% vs 80.0%")
EOF
```

---

## Null result scripts (run independently, not needed for main figures)

```bash
# Phenotypic GxE vs ASE GxE overlap (interesting null)
python3 scripts/null_results/gbs_enrichment.py

# Directional DAP-seq test (null — no directional enrichment in GxE genes)
python3 scripts/null_results/directional_dap.py

# IBM phenotypic scan (R — generates IBM_GxE_results_v3.tsv used by gbs_enrichment.py)
# Note: requires IBM phenotype and genotype data files
Rscript scripts/null_results/ibm_phenotypic_scan.R
```

---

## Directory structure

```
gxe_ase/
├── scripts/
│   ├── 1_deseq_analysis.R        # DESeq2 GxE-ASE model + gene set definitions
│   ├── 2b_dap_seq_analysis.py    # DAP-seq per-gene counts and enrichment
│   ├── 3_sv_enrichment.py        # SV enrichment test
│   ├── 4_threshold_sensitivity.py
│   ├── 5_indel_pwm_scoring.py    # Motif disruption at DAP peaks
│   ├── 6_tf_gxe_mechanism.py     # TF binding mechanism test
│   ├── 7_marginal_sensitivity.py
│   ├── 8_go_enrichment.py
│   ├── 9_manuscript_figures.py   # Main text figures (Figs 1–3)
│   └── null_results/
│       ├── directional_dap.py    # Null: no directional enrichment
│       ├── gbs_enrichment.py     # Null: no phenotypic GxE overlap
│       └── ibm_phenotypic_scan.R # IBM QTL scan (context for above)
├── data/
│   ├── counts_for_DESeq2.txt     # Raw allele-specific counts (input to Step 1)
│   ├── meta_for_DESeq2.txt       # Sample metadata (input to Step 1)
│   ├── annotation.gff            # B73v5 gene annotation
│   ├── tf_annotation.tsv         # Stable TF→gene ID lookup (58 expressed TFs)
│   ├── pangene_table.tsv         # Pan-genome gene classifications
│   ├── Mo17_toB73v5_*.vcf        # Small indels and large SVs (Mo17 vs B73v5)
│   ├── Zm-B73-REFERENCE-NAM-5.0.fa.gz  # B73 reference FASTA (Step 5)
│   └── epigenome/dap_seq/
│       ├── normalized_specific_and_shared_peaks/  # Raw DAP-seq peak files (Steps 2, 5, 6)
│       │   ├── B73v5_MATCH_Mo17-B73v5_specific_B73v5_shared/
│       │   └── Mo17_MATCH_B73v5-Mo17_specific_Mo17_shared/
│       └── gem02_rep_memechip00_onefile/  # MEME motif file (Step 5)
├── results/                      # Output TSVs from analyses
└── figures/                      # PNG and PDF figures
```

---

## Notes on reproducibility

**Dependency chain:**
```
Step 1 (DESeq2)
  └─→ Step 2 (DAP-seq)  ──────────────────────┐
       └─→ Step 4 (threshold sensitivity)      │
       └─→ Step 9 (manuscript figures)         │
  └─→ Step 3 (SV enrichment)                  │
       └─→ Step 4                              │
       └─→ Step 9                              │
  └─→ Step 5 (PWM scoring)                    │
       └─→ Step 9                              │
  └─→ Step 6 (TF mechanism) ←── Step 2 output ┘
       └─→ Step 7 (marginal sensitivity)
  └─→ Step 7
  └─→ Step 8 (GO enrichment)
```

**Files generated from raw data by scripts:**
- Step 1 → `data/GxE_gene_IDs.txt`, `data/background_gene_IDs.txt`,
  `data/DEG_GxE_results.txt`, `data/GxE_allele_specific_test_results.txt`,
  `data/DEG_Env_results.txt`
- Step 2 → `results/dap_seq_per_gene.tsv`, `results/dap_seq_per_tf.tsv`,
  `results/indel_dap_overlap_per_gene.tsv`
- Steps 3–9 → all remaining `results/` and `figures/` files

**Stable reference files (not regenerated by pipeline):**
- `data/tf_annotation.tsv` — maps TF names to B73v5 gene IDs and TF families for
  the 58 TFs with expression data; derived once from the DAP-seq metadata.
- `data/IBM_GxE_results_v3.tsv` — IBM phenotypic QTL scan results (from
  `scripts/null_results/ibm_phenotypic_scan.R`; requires IBM genotype/phenotype data).
- VCFs, GFF, reference FASTA, epigenome peak files — external inputs.
