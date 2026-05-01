"""
threshold_sensitivity.py
------------------------
Re-runs key enrichment tests (DAP-seq, SV, indel × DAP overlap) across
a range of GxE significance thresholds to evaluate how sensitive the
results are to the original loose p < 0.1 cutoff.

Thresholds tested:
  raw pvalue : 0.1 (original), 0.05, 0.01, 0.001
  FDR (padj) : 0.1, 0.05

For each threshold:
  1. DAP enrichment  — Fisher's exact: proportion of genes with ≥1
                       genotype-specific DAP peak in promoter
  2. SV enrichment   — Fisher's exact: proportion of genes with ≥1
                       small indel in promoter; MWU on count
  3. Indel × DAP     — Fisher's exact: proportion of genes where ≥1
                       promoter indel overlaps a genotype-specific peak

Per-gene summary files are loaded once; only gene-set labels change.
"""

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import gzip
from pathlib import Path

BASE = "./"
DATA = f"{BASE}/data"
RESULTS = f"{BASE}/results"
FIGURES = f"{BASE}/figures"

# ── 1. Load GxE test results ──────────────────────────────────────────────────
print("Loading test results …")
gxe_res = pd.read_csv(f"{DATA}/DEG_GxE_results.txt", sep="\t")
# NA padj = failed DESeq2 independent filtering; treat as padj=1 (non-significant)
n_na_padj = gxe_res['padj'].isna().sum()
gxe_res['padj'] = gxe_res['padj'].fillna(1.0)
print(f"  Filled {n_na_padj:,} NA padj values with 1.0 (failed independent filtering)")
# Keep only genes with a valid pvalue AND in the pre-filtered universe.
# Pre-filtered universe = genes in dap_seq_per_gene.tsv (>10 counts in >5 columns).
# Using this as the background denominator prevents inflating ORs by including
# low-expressed, untestable genes.
_dap_preview = pd.read_csv(f"{RESULTS}/dap_seq_per_gene.tsv", sep="\t", usecols=["gene"])
pre_filtered_universe = set(_dap_preview["gene"].tolist())
print(f"  Pre-filtered universe (dap_seq_per_gene): {len(pre_filtered_universe):,} genes")
gxe_valid = gxe_res.dropna(subset=["pvalue"]).copy()
gxe_valid = gxe_valid[gxe_valid["GeneID"].isin(pre_filtered_universe)].copy()
print(f"  Genes with pvalue AND in pre-filtered universe: {len(gxe_valid):,}")

# ── 2. Define thresholds ──────────────────────────────────────────────────────
THRESHOLDS = [
    ("p < 0.100", "pvalue",  0.100),
    ("p < 0.050", "pvalue",  0.050),
    ("p < 0.010", "pvalue",  0.010),
    ("p < 0.001", "pvalue",  0.001),
    ("padj < 0.100", "padj", 0.100),
    ("padj < 0.050", "padj", 0.050),
]

print("\nGene counts per threshold:")
for label, col, thresh in THRESHOLDS:
    sub = gxe_valid.dropna(subset=[col])
    n_gxe = (sub[col] < thresh).sum()
    n_bg  = len(gxe_valid) - n_gxe   # all tested-with-pvalue genes not in GxE set
    print(f"  {label:18s}  GxE={n_gxe:5,}  BG={n_bg:6,}")

# ── 3. Load per-gene data files ───────────────────────────────────────────────
print("\nLoading per-gene summaries …")

dap = pd.read_csv(f"{RESULTS}/dap_seq_per_gene.tsv", sep="\t")
# Rename for consistency
dap = dap.rename(columns={"gene": "GeneID"})
print(f"  dap_seq_per_gene: {len(dap):,} genes")

indel_dap = pd.read_csv(f"{RESULTS}/indel_dap_overlap_per_gene.tsv", sep="\t")
indel_dap = indel_dap.rename(columns={"gene": "GeneID"})
print(f"  indel_dap_overlap_per_gene: {len(indel_dap):,} genes")

# ── 4. Compute SV counts per gene (from VCF + promoter coords) ───────────────
# Check if already computed
sv_path = f"{RESULTS}/sv_per_gene.tsv"
if Path(sv_path).exists():
    print(f"  sv_per_gene.tsv found — loading")
    sv = pd.read_csv(sv_path, sep="\t").rename(columns={"gene": "GeneID"})
else:
    print("\nComputing SV counts per gene (one-time, ~2 min) …")

    VCF_SMALL = f"{BASE}/Mo17_toB73v5_paf_syri_noStartPOS0_INDELs_less50bp.vcf"

    def load_vcf_positions(vcf_path):
        rows = []
        opener = gzip.open if str(vcf_path).endswith(".gz") else open
        with opener(vcf_path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.split("\t", 5)
                chrom = parts[0].strip()
                if not chrom.startswith("chr"):
                    chrom = "chr" + chrom
                try:
                    rows.append((chrom, int(parts[1])))
                except (ValueError, IndexError):
                    continue
        return pd.DataFrame(rows, columns=["chr", "pos"])

    vcf_df = load_vcf_positions(VCF_SMALL)
    print(f"  Loaded {len(vcf_df):,} small indel positions")

    # Vectorised count per promoter using searchsorted
    sv_rows = []
    for chrom, v_grp in vcf_df.groupby("chr"):
        var_pos = np.sort(v_grp["pos"].values)
        mask = dap["chr"] == chrom
        for _, row in dap[mask].iterrows():
            s = int(row["prom_s"])
            e = int(row["prom_e"])
            lo = np.searchsorted(var_pos, s, side="left")
            hi = np.searchsorted(var_pos, e, side="right")
            sv_rows.append({"GeneID": row["GeneID"], "n_small_indel": hi - lo})

    sv = pd.DataFrame(sv_rows)
    # Genes on chromosomes with 0 indels got no rows — fill them
    sv_all = dap[["GeneID"]].copy()
    sv_all = sv_all.merge(sv, on="GeneID", how="left").fillna(0)
    sv_all["n_small_indel"] = sv_all["n_small_indel"].astype(int)
    sv_all["has_sv"] = sv_all["n_small_indel"] > 0
    sv = sv_all
    sv.to_csv(sv_path, sep="\t", index=False)
    print(f"  Saved sv_per_gene.tsv  ({sv['has_sv'].sum():,} genes with ≥1 SV)")

sv["has_sv"] = sv["n_small_indel"] > 0

# ── 5. Merge all per-gene data ────────────────────────────────────────────────
print("\nMerging per-gene data …")

per_gene = dap[["GeneID", "chr", "strand"]].copy()
per_gene["has_geno_spec_dap"] = (
    dap["n_b73_spec"].fillna(0) + dap["n_mo17_spec"].fillna(0) > 0
).values
per_gene["n_geno_spec_dap"] = (
    dap["n_b73_spec"].fillna(0) + dap["n_mo17_spec"].fillna(0)
).values

per_gene = per_gene.merge(
    sv[["GeneID", "n_small_indel", "has_sv"]], on="GeneID", how="left"
)
per_gene = per_gene.merge(
    indel_dap[["GeneID", "n_genotype_spec"]].rename(
        columns={"n_genotype_spec": "n_indel_at_geno_spec"}),
    on="GeneID", how="left"
)
per_gene["has_indel_at_geno_spec"] = per_gene["n_indel_at_geno_spec"].fillna(0) > 0
per_gene["n_small_indel"]          = per_gene["n_small_indel"].fillna(0).astype(int)
per_gene["has_sv"]                 = per_gene["has_sv"].fillna(False).astype(bool)

print(f"  Merged gene table: {len(per_gene):,} genes")

# ── 6. Run tests across thresholds ───────────────────────────────────────────
print("\nRunning sensitivity analysis …")
print("=" * 70)

results = []

for label, col, thresh in THRESHOLDS:
    sub = gxe_valid.dropna(subset=[col])
    gxe_ids = set(sub.loc[sub[col] < thresh, "GeneID"])
    bg_ids  = set(sub.loc[sub[col] >= thresh, "GeneID"])
    n_gxe   = len(gxe_ids)
    n_bg    = len(bg_ids)

    # Subset per-gene table
    g = per_gene[per_gene["GeneID"].isin(gxe_ids)]
    b = per_gene[per_gene["GeneID"].isin(bg_ids)]

    row = dict(threshold=label, n_gxe=n_gxe, n_bg=n_bg)

    print(f"\n── {label}  (GxE n={n_gxe:,}  BG n={n_bg:,}) ──")

    # ── A. DAP-seq enrichment ─────────────────────────────────────────────────
    has_g = g["has_geno_spec_dap"].sum()
    has_b = b["has_geno_spec_dap"].sum()
    ct = [[has_g, n_gxe - has_g],
          [has_b, n_bg  - has_b]]
    or_dap, p_dap = fisher_exact(ct)
    mwu_dap, mwu_p_dap = mannwhitneyu(
        g["n_geno_spec_dap"].fillna(0),
        b["n_geno_spec_dap"].fillna(0),
        alternative="greater")
    print(f"  DAP enrichment:   GxE {has_g}/{n_gxe} ({100*has_g/n_gxe:.1f}%)  "
          f"BG {has_b}/{n_bg} ({100*has_b/n_bg:.1f}%)  "
          f"OR={or_dap:.3f}  Fisher p={p_dap:.3g}  MWU p={mwu_p_dap:.3g}")
    row.update(dap_or=or_dap, dap_p=p_dap, dap_mwu_p=mwu_p_dap,
               pct_gxe_dap=100*has_g/n_gxe, pct_bg_dap=100*has_b/n_bg)

    # ── B. SV enrichment (small indels in promoter) ───────────────────────────
    has_g_sv = g["has_sv"].sum()
    has_b_sv = b["has_sv"].sum()
    ct_sv = [[has_g_sv, n_gxe - has_g_sv],
             [has_b_sv, n_bg  - has_b_sv]]
    or_sv, p_sv = fisher_exact(ct_sv)
    mwu_sv, mwu_p_sv = mannwhitneyu(
        g["n_small_indel"].fillna(0),
        b["n_small_indel"].fillna(0),
        alternative="greater")
    print(f"  SV enrichment:    GxE {has_g_sv}/{n_gxe} ({100*has_g_sv/n_gxe:.1f}%)  "
          f"BG {has_b_sv}/{n_bg} ({100*has_b_sv/n_bg:.1f}%)  "
          f"OR={or_sv:.3f}  Fisher p={p_sv:.3g}  MWU p={mwu_p_sv:.3g}")
    row.update(sv_or=or_sv, sv_p=p_sv, sv_mwu_p=mwu_p_sv,
               pct_gxe_sv=100*has_g_sv/n_gxe, pct_bg_sv=100*has_b_sv/n_bg)

    # ── C. Indel × DAP overlap ────────────────────────────────────────────────
    has_g_id = g["has_indel_at_geno_spec"].sum()
    has_b_id = b["has_indel_at_geno_spec"].sum()
    ct_id = [[has_g_id, n_gxe - has_g_id],
             [has_b_id, n_bg  - has_b_id]]
    or_id, p_id = fisher_exact(ct_id)
    print(f"  Indel × DAP:      GxE {has_g_id}/{n_gxe} ({100*has_g_id/n_gxe:.1f}%)  "
          f"BG {has_b_id}/{n_bg} ({100*has_b_id/n_bg:.1f}%)  "
          f"OR={or_id:.3f}  Fisher p={p_id:.3g}")
    row.update(indel_dap_or=or_id, indel_dap_p=p_id,
               pct_gxe_indap=100*has_g_id/n_gxe, pct_bg_indap=100*has_b_id/n_bg)

    results.append(row)

res_df = pd.DataFrame(results)
print()

# ── 7. Figure ─────────────────────────────────────────────────────────────────
print("Generating figure …")

labels      = res_df["threshold"].tolist()
x           = np.arange(len(labels))
bar_width   = 0.32
orig_idx    = 4   # index of the chosen threshold (padj<0.1)

fig = plt.figure(figsize=(16, 11))
fig.subplots_adjust(hspace=0.55, wspace=0.38)

gs = fig.add_gridspec(3, 3)

ANALYSES = [
    ("DAP enrichment\n(≥1 genotype-specific peak)",
     "dap_or", "dap_p", "pct_gxe_dap", "pct_bg_dap"),
    ("SV enrichment\n(≥1 small indel in promoter)",
     "sv_or", "sv_p", "pct_gxe_sv", "pct_bg_sv"),
    ("Indel × DAP overlap\n(indel at genotype-specific peak)",
     "indel_dap_or", "indel_dap_p", "pct_gxe_indap", "pct_bg_indap"),
]

colors_gxe = "#2196F3"
colors_bg  = "#9E9E9E"
or_color   = "#D32F2F"

for col_i, (title, or_col, p_col, pct_gxe_col, pct_bg_col) in enumerate(ANALYSES):
    ax_bar = fig.add_subplot(gs[0, col_i])
    ax_or  = fig.add_subplot(gs[1, col_i])
    ax_p   = fig.add_subplot(gs[2, col_i])

    pct_gxe = res_df[pct_gxe_col].values
    pct_bg  = res_df[pct_bg_col].values
    ors     = res_df[or_col].values
    ps      = res_df[p_col].values

    # Row 0: % bars GxE vs BG
    ax_bar.bar(x - bar_width/2, pct_gxe, bar_width,
               color=colors_gxe, label="GxE-ASE", edgecolor="black", lw=0.5)
    ax_bar.bar(x + bar_width/2, pct_bg,  bar_width,
               color=colors_bg,  label="Background", edgecolor="black", lw=0.5)
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(labels, rotation=35, ha="right", fontsize=7)
    ax_bar.set_ylabel("% genes", fontsize=8)
    ax_bar.set_title(title, fontsize=8.5, fontweight="bold")
    ax_bar.legend(fontsize=6.5)
    # Mark original threshold
    ax_bar.axvline(orig_idx, color="orange", lw=1.2, ls="--", alpha=0.8, label="chosen (padj<0.1)")

    # Row 1: Odds ratio
    ax_or.plot(x, ors, "o-", color=or_color, lw=2, ms=6)
    ax_or.axhline(1.0, color="black", lw=0.8, ls="--")
    ax_or.axvline(orig_idx, color="orange", lw=1.2, ls="--", alpha=0.8)
    ax_or.set_xticks(x)
    ax_or.set_xticklabels(labels, rotation=35, ha="right", fontsize=7)
    ax_or.set_ylabel("Odds ratio (Fisher)", fontsize=8)
    ax_or.set_title("Odds ratio across thresholds", fontsize=8)

    # Row 2: -log10(p)
    neg_log_p = -np.log10(np.clip(ps, 1e-15, 1))
    ax_p.bar(x, neg_log_p, color=or_color, alpha=0.7, edgecolor="black", lw=0.5)
    ax_p.axhline(-np.log10(0.05), color="gray",  lw=1, ls="--", label="p=0.05")
    ax_p.axhline(-np.log10(0.01), color="black", lw=1, ls=":",  label="p=0.01")
    ax_p.axvline(orig_idx, color="orange", lw=1.2, ls="--", alpha=0.8, label="chosen (padj<0.1)")
    ax_p.set_xticks(x)
    ax_p.set_xticklabels(labels, rotation=35, ha="right", fontsize=7)
    ax_p.set_ylabel("-log₁₀(Fisher p)", fontsize=8)
    ax_p.set_title("Significance across thresholds", fontsize=8)
    ax_p.legend(fontsize=6.5)

fig.suptitle(
    "Sensitivity of enrichment results to GxE-ASE significance threshold\n"
    "Orange dashed line = chosen threshold (padj < 0.1)",
    fontsize=11, fontweight="bold"
)

out_fig = f"{FIGURES}/threshold_sensitivity_figure.png"
fig.savefig(out_fig, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved → {out_fig}")

# ── 8. Summary table ──────────────────────────────────────────────────────────
print("\nSummary table:")
display_cols = ["threshold", "n_gxe",
                "dap_or", "dap_p",
                "sv_or",  "sv_p",
                "indel_dap_or", "indel_dap_p"]

print(res_df[display_cols].to_string(index=False,
    float_format=lambda x: f"{x:.4g}" if abs(x) < 100 else f"{x:.1f}"))

out_tsv = f"{RESULTS}/threshold_sensitivity_results.tsv"
res_df.to_csv(out_tsv, sep="\t", index=False)
print(f"\nSaved → {out_tsv}")
print("Done.")
