"""
sv_enrichment.py
----------------
Tests whether GxE-ASE genes have more B73/Mo17 structural variants (small
indels <50bp, large indels >50bp) in their promoter regions compared to
background genes.

Hypothesis: cis-regulatory sequence differences between B73 and Mo17 in the
promoter region are a mechanistic driver of GxE-ASE, so GxE-ASE genes should
be enriched for any class of SV in their promoters.

Tests:
  1. Binary Fisher's exact: proportion of genes with ≥1 variant (any class)
  2. Count comparison: Mann-Whitney U on variant count per gene
  3. Per-class breakdown (small indels, large indels)
  4. Effect of variant class and window distance from TSS
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
import gzip

BASE_DIR  = "./"
DATA_DIR  = f"{BASE_DIR}/data"
RESULTS   = f"{BASE_DIR}/results"
FIGURES   = f"{BASE_DIR}/figures"
PROMO_WIN = 500  # bp from TSS

VCF_FILES = {
    "small_indel": f"{DATA_DIR}/Mo17_toB73v5_paf_syri_noStartPOS0_INDELs_less50bp.vcf",
    "large_indel": f"{DATA_DIR}/Mo17_toB73v5_paf_syri_noStartPOS0_INDELs_more50bp.vcf",
    "SNP": f"{DATA_DIR}/Mo17_toB73v5_paf_syri_noStartPOS0_SNPs.vcf"
}

# ── Load promoter regions ─────────────────────────────────────────────────────
print("Loading promoter regions from dap_seq_per_gene.tsv …")
dap = pd.read_csv(f"{RESULTS}/dap_seq_per_gene.tsv", sep="\t")

gxe_ids = set(pd.read_csv(f"{DATA_DIR}/GxE_gene_IDs.txt",        header=None)[0])
bg_ids  = set(pd.read_csv(f"{DATA_DIR}/background_gene_IDs.txt", header=None)[0])

dap["is_gxe"] = dap["gene"].isin(gxe_ids)
dap["is_bg"]  = dap["gene"].isin(bg_ids)
dap = dap.reset_index(drop=True)

print(f"  GxE-ASE genes:  {dap['is_gxe'].sum():,}")
print(f"  Background genes: {dap['is_bg'].sum():,}")

# ── VCF loader — returns dict: chrom → sorted pos array ──────────────────────
def load_vcf_by_chr(vcf_path):
    """
    Load VCF positions into a dict: chrom → sorted int64 numpy array.
    Pre-sorting once per VCF is the key to making per-gene lookups fast.
    """
    buf = {}
    print(f"  Loading {Path(vcf_path).name} …")
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 5)
            chrom = parts[0].strip()
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            buf.setdefault(chrom, []).append(int(parts[1]))
    vcf = {chrom: np.array(sorted(pos), dtype=np.int64) for chrom, pos in buf.items()}
    total = sum(len(v) for v in vcf.values())
    print(f"    {total:,} variants across {len(vcf)} chromosomes")
    return vcf

# ── Count variants per gene using pre-built chrom arrays ─────────────────────
def count_variants_in_promoters(vcf_by_chr, promoters_df):
    """
    For each gene's promoter [prom_s, prom_e), count overlapping variants.
    vcf_by_chr: dict chrom → sorted int64 pos array (built once per VCF).
    Returns int array aligned with promoters_df row order.
    """
    counts = np.zeros(len(promoters_df), dtype=np.int32)
    for chrom, pos_arr in vcf_by_chr.items():
        mask = (promoters_df["chr"] == chrom).values
        if not mask.any():
            continue
        idx   = np.where(mask)[0]
        prom_s = promoters_df["prom_s"].values[idx]
        prom_e = promoters_df["prom_e"].values[idx]
        lo = np.searchsorted(pos_arr, prom_s, side="left")
        hi = np.searchsorted(pos_arr, prom_e, side="right")
        counts[idx] = hi - lo
    return counts

# ── Count variants per gene for each variant class ────────────────────────────
print("\nCounting variants per gene …")
vcf_cache = {}   # keep loaded VCFs for the distance analysis below

for vcf_class, vcf_path in VCF_FILES.items():
    print(f"\n  Processing {vcf_class} …")
    vcf_by_chr = load_vcf_by_chr(vcf_path)
    vcf_cache[vcf_class] = vcf_by_chr
    col = f"n_{vcf_class}"
    dap[col] = count_variants_in_promoters(vcf_by_chr, dap)
    print(f"    Genes with ≥1 variant: {(dap[col] > 0).sum():,}")
    print(f"    Mean per gene: {dap[col].mean():.3f}")

sv_cols = [f"n_{c}" for c in VCF_FILES.keys()]
dap["n_sv_total"] = dap[sv_cols].sum(axis=1)
dap["has_any_sv"] = dap["n_sv_total"] > 0
print(f"\n  Total genes with ≥1 SV in promoter: {dap['has_any_sv'].sum():,}")

# ── Statistical tests ─────────────────────────────────────────────────────────
print("\n── Statistical Tests ──────────────────────────────────────────────")

results = []
for vcf_class in list(VCF_FILES.keys()) + ["sv_total"]:
    col = f"n_{vcf_class}"
    g = dap[dap["is_gxe"]]
    b = dap[dap["is_bg"]]

    has_g = (g[col] > 0).sum()
    has_b = (b[col] > 0).sum()
    ct = np.array([[has_g, len(g) - has_g],
                   [has_b, len(b) - has_b]])
    or_fish, p_fish = stats.fisher_exact(ct)
    mwu_stat, mwu_p = stats.mannwhitneyu(g[col].fillna(0), b[col].fillna(0),
                                          alternative="greater")

    print(f"\n  {vcf_class}:")
    print(f"    GxE: {has_g}/{len(g)} have ≥1 ({100*has_g/len(g):.1f}%),  "
          f"mean={g[col].mean():.3f}")
    print(f"    BG:  {has_b}/{len(b)} have ≥1 ({100*has_b/len(b):.1f}%),  "
          f"mean={b[col].mean():.3f}")
    print(f"    Fisher OR={or_fish:.3f}, p={p_fish:.4g}")
    print(f"    MWU (count, one-sided GxE>BG): p={mwu_p:.4g}")

    results.append(dict(variant_class=vcf_class, n_gxe=len(g), n_bg=len(b),
                        has_gxe=has_g, has_bg=has_b,
                        pct_gxe=100*has_g/len(g), pct_bg=100*has_b/len(b),
                        mean_gxe=g[col].mean(), mean_bg=b[col].mean(),
                        fisher_or=or_fish, fisher_p=p_fish,
                        mwu_p=mwu_p))

res_df = pd.DataFrame(results)

# ── Distance-from-TSS analysis ─────────────────────────────────────────────────
print("\n── Distance-from-TSS analysis (small indels, binned) ──")

# Build TSS position arrays per chromosome — computed once, outside all loops.
# tss_by_chr: chrom → (tss_positions array, is_gxe bool array, is_bg bool array,
#                       upstream_sign array)
# upstream_sign: +1 strand genes look upstream (lower coords), -1 strand genes
# look upstream (higher coords). For a ± window we just use absolute distance,
# so we store the raw TSS and compute the window symmetrically.
#
# For each gene and bin [b_start, b_end):
#   + strand: upstream window = [tss - b_end, tss - b_start)
#   - strand: upstream window = [tss + b_start, tss + b_end)

# Pre-extract columns as numpy arrays for speed
chroms   = dap["chr"].values
prom_s   = dap["prom_s"].values
prom_e   = dap["prom_e"].values
strands  = dap["strand"].values
is_gxe   = dap["is_gxe"].values
is_bg    = dap["is_bg"].values

# TSS: prom_s + PROMO_WIN for + strand; prom_e - PROMO_WIN for - strand
tss = np.where(strands == "+", prom_s + PROMO_WIN, prom_e - PROMO_WIN)

# Use the already-loaded small-indel VCF
vcf_small = vcf_cache["small_indel"]

bin_edges  = np.arange(0, PROMO_WIN + 1, 100)
bin_labels = [f"{b}–{b+100}" for b in bin_edges[:-1]]

# Pre-group gene indices by chromosome so the outer bin loop doesn't re-scan
genes_by_chr = {}
for chrom in set(chroms):
    genes_by_chr[chrom] = np.where(chroms == chrom)[0]

bin_results = []
for b_start, b_end, label in zip(bin_edges[:-1], bin_edges[1:], bin_labels):
    # Per-gene binary flag: has ≥1 small indel in this distance band?
    has_indel = np.zeros(len(dap), dtype=bool)

    for chrom, gene_idx in genes_by_chr.items():
        if chrom not in vcf_small:
            continue
        pos_arr = vcf_small[chrom]   # pre-sorted array, built once above

        g_tss     = tss[gene_idx]
        g_strands = strands[gene_idx]

        # Upstream window per gene (vectorised)
        win_s = np.where(g_strands == "+", g_tss - b_end,   g_tss + b_start)
        win_e = np.where(g_strands == "+", g_tss - b_start, g_tss + b_end)

        lo = np.searchsorted(pos_arr, win_s, side="left")
        hi = np.searchsorted(pos_arr, win_e, side="right")
        has_indel[gene_idx] = (hi - lo) > 0

    pct_gxe = 100 * has_indel[is_gxe].mean() if is_gxe.any() else 0.0
    pct_bg  = 100 * has_indel[is_bg].mean()  if is_bg.any()  else 0.0
    bin_results.append(dict(distance_band=label, d_start=b_start,
                            pct_gxe=pct_gxe, pct_bg=pct_bg))
    print(f"  {label}: GxE {pct_gxe:.1f}%  BG {pct_bg:.1f}%")

bin_df = pd.DataFrame(bin_results)

# ── Main enrichment figure ────────────────────────────────────────────────────
print("\nGenerating figures …")

fig, axes = plt.subplots(2, 4, figsize=(16, 10))
fig.subplots_adjust(hspace=0.4, wspace=0.4)

class_labels = {
    "small_indel": "Small indels (<50 bp)",
    "large_indel": "Large indels (>50 bp)",
    "SNP": "SNPs",
    "sv_total":    "All SVs combined",
}
colors = {"GxE-ASE": "#2196F3", "Background": "#9E9E9E"}

for ax_i, vcf_class in enumerate(["small_indel", "large_indel", "SNP", "sv_total"]):
    col = f"n_{vcf_class}"
    r = res_df[res_df["variant_class"] == vcf_class].iloc[0]

    # Top row: % genes with ≥1 SV
    ax = axes[0, ax_i]
    bars = ax.bar(["GxE-ASE", "Background"],
                  [r["pct_gxe"], r["pct_bg"]],
                  color=[colors["GxE-ASE"], colors["Background"]],
                  edgecolor="black", linewidth=0.5)
    ax.set_ylabel("% genes with ≥1 variant")
    ax.set_title(f"{chr(65+ax_i)}  {class_labels[vcf_class]}\n"
                 f"Fisher OR={r['fisher_or']:.3f}, p={r['fisher_p']:.3g}\n"
                 f"MWU p={r['mwu_p']:.3g}", fontsize=8)
    for bar, val in zip(bars, [r["pct_gxe"], r["pct_bg"]]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                f"{val:.1f}%", ha="center", va="bottom", fontsize=8)

    # Bottom row: count distribution
    ax2 = axes[1, ax_i]
    g_counts = dap[dap["is_gxe"]][col].fillna(0)
    b_counts = dap[dap["is_bg"]][col].fillna(0)
    cap = min(max(g_counts.max(), b_counts.max()),
              int(np.percentile(np.concatenate([g_counts, b_counts]), 98)) + 1)
    bins = np.arange(0, cap + 2) - 0.5
    ax2.hist(g_counts.clip(upper=cap), bins=bins, color=colors["GxE-ASE"],
             alpha=0.6, density=True, label=f"GxE (n={len(g_counts)})",
             histtype="stepfilled")
    ax2.hist(b_counts.clip(upper=cap), bins=bins, color=colors["Background"],
             alpha=0.4, density=True, label=f"BG (n={len(b_counts)})",
             histtype="step", lw=1.5)
    ax2.set_xlabel(f"Variants per promoter (±{PROMO_WIN} bp)")
    ax2.set_ylabel("Density")
    ax2.legend(fontsize=7)
    ax2.set_title(f"{chr(68+ax_i)}  Mean: GxE={r['mean_gxe']:.2f}, "
                  f"BG={r['mean_bg']:.2f}", fontsize=8)

fig.suptitle(f"Structural Variant Enrichment in Promoters of GxE-ASE Genes\n"
             f"±{PROMO_WIN} bp from TSS | B73/Mo17 SV calls (syri)",
             fontsize=11, fontweight="bold")
out_fig = f"{FIGURES}/sv_enrichment_figure.png"
fig.savefig(out_fig, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved → {out_fig}")

# ── Distance-from-TSS figure ──────────────────────────────────────────────────
fig2, ax = plt.subplots(figsize=(8, 4))
x = np.arange(len(bin_df))
ax.plot(x, bin_df["pct_gxe"], "o-", color=colors["GxE-ASE"],  lw=2, label="GxE-ASE")
ax.plot(x, bin_df["pct_bg"],  "s--", color=colors["Background"], lw=1.5, label="Background")
ax.set_xticks(x)
ax.set_xticklabels(bin_df["distance_band"], rotation=30, ha="right")
ax.set_xlabel("Distance from TSS (bp, upstream)")
ax.set_ylabel("% genes with small indel")
ax.set_title("Small Indel Enrichment by Distance from TSS")
ax.legend()
ax2_r = ax.twinx()
ax2_r.plot(x, bin_df["pct_gxe"] - bin_df["pct_bg"], "^-",
           color="darkred", lw=1.5, alpha=0.7, label="GxE − BG")
ax2_r.set_ylabel("Difference (GxE − BG, %)")
ax2_r.legend(loc="upper right")
fig2.tight_layout()
out_fig2 = f"{FIGURES}/sv_distance_figure.png"
fig2.savefig(out_fig2, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved → {out_fig2}")

# ── Save results ──────────────────────────────────────────────────────────────
res_df.to_csv(f"{RESULTS}/sv_enrichment_results.tsv", sep="\t", index=False)
sv_out_cols = ["gene", "chr", "strand", "is_gxe", "is_bg",
               "n_small_indel", "n_large_indel", "n_SNP", "n_sv_total", "has_any_sv"]
dap[sv_out_cols].to_csv(f"{RESULTS}/sv_per_gene.tsv", sep="\t", index=False)
print(f"  Saved → {RESULTS}/sv_enrichment_results.tsv")
print(f"  Saved → {RESULTS}/sv_per_gene.tsv")
print("\nDone.")
