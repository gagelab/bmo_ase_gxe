"""
8_gxe_qtl_enrichment.py
------------------------
Test whether GxE-ASE genes from the allele-specific expression analysis
co-localize with GxE phenotypic QTL bins from the IBM population scan.

Analysis pipeline:
  1. Load QTL scan results (IBM_GxE_results.tsv); combine p-values across
     14 non-redundant traits with Stouffer's method. BH-FDR is computed for
     reference but a naive raw p-value threshold (SIG_P_THRESH, default
     0.05) is used to flag "signal" bins -- BH-FDR across ~15k bins/trait
     is far too conservative for this screening-level overlap test.
  2. Assign each background gene to a recombination bin by midpoint.
  3. Binary enrichment: Fisher's exact + chromosome-stratified permutation
     (are GxE-ASE genes more likely to fall in a bin with raw
     stouffer_GEp < SIG_P_THRESH?).
  4. Continuous enrichment: Spearman and partial Spearman (controlling for
     baseMean) of bin QTL signal vs ASE -log10(p).
  5. Trait-stratified binary enrichment (per trait, raw GEp < SIG_P_THRESH).

Reads:
  results/IBM_GxE_results.tsv         (output of 7_map_gxe_qtl.R)
  data/IBM_recomb_bins_fromGBS.tsv    (output of 6_format_gbs_genotypes.R)
  data/annotation.gff
  data/GxE_allele_specific_test_results.txt
  data/GxE_gene_IDs.txt
  data/background_gene_IDs.txt

Writes:
  results/gbs_bin_stouffer_stats.tsv
  results/gbs_gene_bin_assignments.tsv
  results/gbs_enrichment_results.tsv
  results/gbs_enrichment_by_trait.tsv
  figures/gbs_enrichment_figure.png
"""

import re
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

BASE     = "./"
RESULTS  = f"{BASE}/results"
FIGS     = f"{BASE}/figures"

# Gene lists, ASE results, and annotation live in the shared reference data dir
# (these files are not duplicated in bmo_ase_gxe/data/)
DATA = f"{BASE}/data"

QTL_FILE  = f"{RESULTS}/IBM_GxE_results.tsv"
BINS_FILE = f"{DATA}/IBM_recomb_bins_fromGBS.tsv"
GFF_FILE  = f"{DATA}/annotation.gff"
ASE_FILE  = f"{DATA}/GxE_allele_specific_test_results.txt"
GXE_IDS   = f"{DATA}/GxE_gene_IDs.txt"
BG_IDS    = f"{DATA}/background_gene_IDs.txt"

SEED    = 2126
N_PERM  = 10_000

# Naive (uncorrected) p-value threshold used to flag a recombination bin as
# "showing GxE QTL signal". This analysis is a screening-level overlap test,
# not a genome-wide significance scan, so BH-FDR across ~15k bins/trait is
# far too conservative (it would require raw p <~ 3e-6 for the strongest
# bin to survive). A flat threshold gives a consistent, interpretable
# definition of "signal" for both the combined (Stouffer) bins and the
# per-trait bins.
SIG_P_THRESH = 0.01

# 14 non-redundant traits (drop redundant phenology and size composites)
KEEP_TRAITS = [
    "cob_diameter", "cob_length", "cob_mass", "days_to_silk",
    "ear_row_num", "kernels_per_row", "leaf_length", "leaf_width",
    "ph", "tassel_length", "tassel_primary_branch_num",
    "total_kernel_weight", "twenty_kernel_weight", "upper_leaf_angle",
]

# ── Helpers ───────────────────────────────────────────────────────────────────

def stouffer_p(pvec):
    """Combine p-values with Stouffer's Z-score method."""
    pv = np.clip(np.asarray(pvec, float), 1e-10, 1 - 1e-10)
    Z  = stats.norm.ppf(1 - pv).sum() / np.sqrt(len(pv))
    return float(stats.norm.sf(Z))

def partial_spearman(x, y, z):
    """Spearman correlation between x and y after partialling out z (rank-based)."""
    rx, ry, rz = [stats.rankdata(v) for v in (x, y, z)]
    def resid(rv):
        rz_c = rz - rz.mean(); rv_c = rv - rv.mean()
        b = (rz_c @ rv_c) / (rz_c @ rz_c) if (rz_c @ rz_c) > 0 else 0.0
        return rv_c - b * rz_c
    ex, ey = resid(rx), resid(ry)
    r, _   = stats.spearmanr(ex, ey)
    return float(r), ex, ey

# ── 1. Load QTL results and compute per-bin Stouffer p ───────────────────────
print("1. Loading QTL scan results ...", flush=True)
qtl = pd.read_csv(QTL_FILE, sep="\t")
qtl = qtl[qtl["trait"].isin(KEEP_TRAITS)].copy()
print(f"   {len(qtl):,} rows  ({qtl['trait'].nunique()} traits  x "
      f"{qtl.groupby(['chr','start','end']).ngroups:,} unique bins)")

qtl_pivot = qtl.pivot_table(
    index=["chr", "start", "end"], columns="trait", values="GEp", aggfunc="first"
)

print("   Combining traits with Stouffer's method ...", flush=True)
stouffer_rows = []
for (chr_, start, end), row in qtl_pivot.iterrows():
    pvec = row.dropna().values
    if len(pvec) < 2:
        continue
    stouffer_rows.append({
        "Chr": int(chr_), "Start": int(start), "End": int(end),
        "n_traits": len(pvec), "stouffer_GEp": stouffer_p(pvec),
    })

bin_stouffer = pd.DataFrame(stouffer_rows)
# BH-FDR is computed and retained in the output for reference, but "signal"
# bins are defined below via a naive raw p-value threshold (SIG_P_THRESH) --
# BH-FDR across ~15k bins/trait is far too conservative for this
# screening-level overlap test.
_, fdr, _, _ = multipletests(bin_stouffer["stouffer_GEp"].values, method="fdr_bh")
bin_stouffer["stouffer_FDR"] = fdr
n_sig = (bin_stouffer["stouffer_GEp"] < SIG_P_THRESH).sum()
print(f"   {len(bin_stouffer):,} testable bins; {n_sig:,} with raw stouffer_GEp < {SIG_P_THRESH} "
      f"(naive threshold, not BH-FDR)")
print(f"   Stouffer GEp median: {bin_stouffer['stouffer_GEp'].median():.4f}")

# ── 2. Parse gene coordinates from GFF ───────────────────────────────────────
print("\n2. Parsing gene coordinates ...", flush=True)
gene_rows = []
with open(GFF_FILE) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.strip().split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        m = re.match(r"(?:chr)?(\d+)$", f[0], re.I)
        if not m or not (1 <= int(m.group(1)) <= 10):
            continue
        try:
            s, e = int(f[3]), int(f[4])
        except ValueError:
            continue
        gid_m = re.search(r"ID=([^;]+)", f[8])
        if not gid_m:
            continue
        gid = gid_m.group(1).split(":")[0]
        gene_rows.append({"GeneID": gid, "chr": int(m.group(1)),
                          "gene_start": s, "gene_end": e,
                          "gene_mid": (s + e) / 2})

genes = pd.DataFrame(gene_rows).drop_duplicates("GeneID").reset_index(drop=True)
print(f"   {len(genes):,} genes")

# ── 3. Load ASE results and gene sets ─────────────────────────────────────────
print("\n3. Loading ASE results ...", flush=True)
ase = pd.read_csv(ASE_FILE, sep="\t").dropna(subset=["pvalue"]).copy()
ase["neg_log10_p"] = -np.log10(ase["pvalue"].clip(1e-300, 1))
bg  = set(pd.read_csv(BG_IDS,  header=None)[0])
sig = set(pd.read_csv(GXE_IDS, header=None)[0])
print(f"   Background: {len(bg):,} | GxE-ASE: {len(sig):,}")

ase_full = ase.merge(genes[["GeneID","chr","gene_start","gene_end","gene_mid"]],
                     on="GeneID", how="inner")
print(f"   {len(ase_full):,} ASE genes with physical coordinates")

# ── 4. Assign genes to bins by midpoint ──────────────────────────────────────
print("\n4. Assigning genes to recombination bins ...", flush=True)

n_genes = len(ase_full)
# GEp and FDR initialised as NaN (float); bin coordinates as -1 (int)
assigned = {
    "stouffer_GEp": np.full(n_genes, np.nan, dtype=float),
    "stouffer_FDR": np.full(n_genes, np.nan, dtype=float),
    "bin_start":    np.full(n_genes, -1,     dtype=np.int64),
    "bin_end":      np.full(n_genes, -1,     dtype=np.int64),
}
in_bin    = np.zeros(n_genes, dtype=bool)

for chr_num in range(1, 11):
    g_idx = np.where(ase_full["chr"].values == chr_num)[0]
    b_sub = bin_stouffer[bin_stouffer["Chr"] == chr_num].sort_values("Start").reset_index(drop=True)
    if len(g_idx) == 0 or len(b_sub) == 0:
        continue
    mids   = ase_full["gene_mid"].values[g_idx]
    starts = b_sub["Start"].values
    ends   = b_sub["End"].values
    for ii, gi in enumerate(g_idx):
        mid  = mids[ii]
        mask = (starts <= mid) & (ends >= mid)
        hits = np.where(mask)[0]
        if hits.size == 0:
            continue
        b = b_sub.iloc[hits[0]]
        in_bin[gi]                    = True
        assigned["stouffer_GEp"][gi]  = b["stouffer_GEp"]
        assigned["stouffer_FDR"][gi]  = b["stouffer_FDR"]
        assigned["bin_start"][gi]     = int(b["Start"])
        assigned["bin_end"][gi]       = int(b["End"])

ase_full["in_bin"]       = in_bin
ase_full["stouffer_GEp"] = assigned["stouffer_GEp"]
ase_full["stouffer_FDR"] = assigned["stouffer_FDR"]
ase_full["bin_start"]    = assigned["bin_start"].astype(int)
ase_full["bin_end"]      = assigned["bin_end"].astype(int)
ase_full["bin_sig"]      = (assigned["stouffer_GEp"] < SIG_P_THRESH) & in_bin

# BUG FIX: include both bg and sig genes (they are mutually exclusive sets,
# so isin(bg) alone would exclude all 248 GxE genes)
ase_bg = ase_full[ase_full["GeneID"].isin(bg | sig) & ase_full["in_bin"]].copy()
ase_bg["is_gxe"] = ase_bg["GeneID"].isin(sig)
print(f"   {len(ase_bg):,} universe genes assigned to a bin")
print(f"   GxE-ASE genes in background+bin: {ase_bg['is_gxe'].sum():,}")
print(f"   In bins with raw p < {SIG_P_THRESH}:   {ase_bg['bin_sig'].sum():,}")

# ── 5. Binary enrichment ─────────────────────────────────────────────────────
print("\n" + "="*65, flush=True)
print("5. BINARY ENRICHMENT", flush=True)

in_sig  = ase_bg["bin_sig"].values
is_gxe  = ase_bg["is_gxe"].values
a = int(( is_gxe &  in_sig).sum())
b = int(( is_gxe & ~in_sig).sum())
c = int((~is_gxe &  in_sig).sum())
d = int((~is_gxe & ~in_sig).sum())

print(f"\n   Contingency table (raw stouffer_GEp < {SIG_P_THRESH} QTL bins):")
print(f"   {'':20}  {'In sig bin':>10}  {'Not in sig bin':>14}")
print(f"   {'GxE-ASE':20}  {a:>10,}  {b:>14,}")
print(f"   {'Non-GxE':20}  {c:>10,}  {d:>14,}")

OR_obs = (a * d) / (b * c) if b * c > 0 else np.inf
OR_fisher, p_fisher = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
print(f"\n   OR (obs) = {OR_obs:.4f}")
print(f"   Fisher exact: OR={OR_fisher:.4f}  p={p_fisher:.4f}")

print(f"\n   Chromosome-stratified permutation (N={N_PERM:,}) ...", flush=True)
rng       = np.random.default_rng(SEED)
chr_arr   = ase_bg["chr"].values
chr_grps  = {c: np.where(chr_arr == c)[0] for c in np.unique(chr_arr)}
gxe_float = is_gxe.astype(float)
sig_int   = in_sig.astype(int)

perm_ors = np.zeros(N_PERM)
for i in range(N_PERM):
    pl = gxe_float.copy()
    for idx in chr_grps.values():
        tmp = pl[idx].copy(); rng.shuffle(tmp); pl[idx] = tmp
    pa = int((pl.astype(bool) &  sig_int.astype(bool)).sum())
    pb = int((pl.astype(bool) & ~sig_int.astype(bool)).sum())
    pc = int((~pl.astype(bool) &  sig_int.astype(bool)).sum())
    pd_ = int((~pl.astype(bool) & ~sig_int.astype(bool)).sum())
    perm_ors[i] = (pa * pd_) / (pb * pc) if pb * pc > 0 else np.nan

perm_ors = perm_ors[~np.isnan(perm_ors)]
perm_p   = (np.sum(perm_ors >= OR_obs) + 1) / (len(perm_ors) + 1)
print(f"   Permutation p (one-sided >= obs OR): {perm_p:.4f}")
print(f"   Null OR: mean={perm_ors.mean():.4f} sd={perm_ors.std():.4f}")

# ── 6. Continuous enrichment ─────────────────────────────────────────────────
print("\n" + "="*65, flush=True)
print("6. CONTINUOUS ENRICHMENT", flush=True)

cont_df = ase_bg[ase_bg["stouffer_GEp"].notna()].copy()
cont_df["qtl_signal"] = -np.log10(cont_df["stouffer_GEp"].clip(1e-300, 1))
cont_df["ase_signal"] = cont_df["neg_log10_p"]
cont_df["log10_bm"]   = np.log10(cont_df["baseMean"].clip(0.01))

r_sp, p_sp = stats.spearmanr(cont_df["qtl_signal"], cont_df["ase_signal"])
print(f"\n   Standard Spearman:  r={r_sp:.4f}  p={p_sp:.4f}")

r_partial, ex, ey = partial_spearman(
    cont_df["qtl_signal"].values,
    cont_df["ase_signal"].values,
    cont_df["log10_bm"].values,
)
n_c      = len(cont_df)
t_part   = r_partial * np.sqrt((n_c - 2) / (1 - r_partial**2 + 1e-30))
p_partial = 2 * stats.t.sf(abs(t_part), df=n_c - 2)
print(f"   Partial Spearman (controlling baseMean): r={r_partial:.4f}  p={p_partial:.4f}")

gxe_vals = cont_df.loc[ cont_df["bin_sig"], "ase_signal"].values
non_vals  = cont_df.loc[~cont_df["bin_sig"], "ase_signal"].values
mw_p = rbc = np.nan
if len(gxe_vals) > 0 and len(non_vals) > 0:
    U, mw_p = stats.mannwhitneyu(gxe_vals, non_vals, alternative="greater")
    rbc = (2 * U / (len(gxe_vals) * len(non_vals))) - 1
    print(f"   Mann-Whitney (sig-bin > not): rbc={rbc:.4f}  p={mw_p:.4f}")

print(f"\n   Permutation for Spearman (N={N_PERM:,}) ...", flush=True)
rng2      = np.random.default_rng(SEED + 1)
chr_arr2  = cont_df["chr"].values
chr_grps2 = {c: np.where(chr_arr2 == c)[0] for c in np.unique(chr_arr2)}
qtl_vals2 = cont_df["qtl_signal"].values
ase_vals2 = cont_df["ase_signal"].values

perm_rs = np.zeros(N_PERM)
for i in range(N_PERM):
    pq = qtl_vals2.copy()
    for idx in chr_grps2.values():
        tmp = pq[idx].copy(); rng2.shuffle(tmp); pq[idx] = tmp
    perm_rs[i], _ = stats.spearmanr(pq, ase_vals2)

perm_p_sp = (np.sum(perm_rs >= r_sp) + 1) / (N_PERM + 1)
print(f"   Permutation Spearman p: {perm_p_sp:.4f}")

# ── 7. Trait-stratified binary enrichment ────────────────────────────────────
print("\n" + "="*65, flush=True)
print("7. TRAIT-STRATIFIED BINARY ENRICHMENT", flush=True)

trait_results = []
for trait in KEEP_TRAITS:
    tq = qtl[qtl["trait"] == trait][["chr","start","end","GEp"]].dropna(subset=["GEp"])
    # Naive raw p-value threshold (see SIG_P_THRESH) -- BH-FDR across
    # ~15k bins/trait is far too conservative for this screening test.
    sig_set = set(
        zip(tq.loc[tq["GEp"] < SIG_P_THRESH, "chr"].astype(int),
            tq.loc[tq["GEp"] < SIG_P_THRESH, "start"].astype(int),
            tq.loc[tq["GEp"] < SIG_P_THRESH, "end"].astype(int))
    )
    if not sig_set:
        trait_results.append({"trait": trait, "n_sig_bins": 0,
                               "OR": np.nan, "fisher_p": np.nan, "perm_p": np.nan})
        continue

    ase_t = ase_bg.copy()
    ase_t["bin_sig_t"] = ase_t.apply(
        lambda r: (int(r["chr"]), int(r["bin_start"]), int(r["bin_end"])) in sig_set
        if r["bin_start"] >= 0 else False, axis=1
    )
    at = int(( ase_t["is_gxe"] &  ase_t["bin_sig_t"]).sum())
    bt = int(( ase_t["is_gxe"] & ~ase_t["bin_sig_t"]).sum())
    ct = int((~ase_t["is_gxe"] &  ase_t["bin_sig_t"]).sum())
    dt = int((~ase_t["is_gxe"] & ~ase_t["bin_sig_t"]).sum())
    OR_t, pf_t = stats.fisher_exact([[at, bt], [ct, dt]], alternative="two-sided")
    OR_obs_t   = (at * dt) / (bt * ct) if bt * ct > 0 else np.nan

    rng_t  = np.random.default_rng(SEED + 100)
    ig_t   = ase_t["is_gxe"].values.astype(float)
    bs_t   = ase_t["bin_sig_t"].values.astype(int)
    chr_t  = ase_t["chr"].values
    cg_t   = {c: np.where(chr_t == c)[0] for c in np.unique(chr_t)}
    p_ors  = np.zeros(1000)
    for pi in range(1000):
        pl = ig_t.copy()
        for idx in cg_t.values():
            tmp = pl[idx].copy(); rng_t.shuffle(tmp); pl[idx] = tmp
        pa = int((pl.astype(bool) &  bs_t.astype(bool)).sum())
        pb = int((pl.astype(bool) & ~bs_t.astype(bool)).sum())
        pc = int((~pl.astype(bool) &  bs_t.astype(bool)).sum())
        pd_ = int((~pl.astype(bool) & ~bs_t.astype(bool)).sum())
        p_ors[pi] = (pa * pd_) / (pb * pc) if pb * pc > 0 else np.nan
    p_ors  = p_ors[~np.isnan(p_ors)]
    pp_t   = (np.sum(p_ors >= OR_obs_t) + 1) / (len(p_ors) + 1) if not np.isnan(OR_obs_t) else np.nan

    trait_results.append({
        "trait": trait, "n_sig_bins": len(sig_set),
        "n_in_sig_gxe": at, "n_in_sig_nongxe": ct,
        "OR": OR_obs_t, "fisher_p": pf_t, "perm_p": pp_t,
    })
    print(f"   {trait:<35}  OR={OR_obs_t:6.3f}  fisher_p={pf_t:.4f}  "
          f"perm_p={pp_t:.4f}  n_sig_bins={len(sig_set):,}")

trait_df = pd.DataFrame(trait_results)

# ── 8. Save results ───────────────────────────────────────────────────────────
print("\nSaving results ...", flush=True)
import os; os.makedirs(RESULTS, exist_ok=True); os.makedirs(FIGS, exist_ok=True)

bin_stouffer.to_csv(f"{RESULTS}/gbs_bin_stouffer_stats.tsv", sep="\t", index=False)

ase_bg[["GeneID","chr","gene_start","gene_end","is_gxe",
        "baseMean","neg_log10_p","stouffer_GEp","stouffer_FDR",
        "bin_start","bin_end","bin_sig"]].to_csv(
    f"{RESULTS}/gbs_gene_bin_assignments.tsv", sep="\t", index=False)

pd.DataFrame([{
    "binary_OR_obs":      OR_obs,   "binary_OR_fisher":   OR_fisher,
    "binary_fisher_p":    p_fisher, "binary_perm_p":      perm_p,
    "spearman_r":         r_sp,     "spearman_p":         p_sp,
    "spearman_perm_p":    perm_p_sp,
    "partial_spearman_r": r_partial, "partial_spearman_p": p_partial,
    "mw_rbc":             rbc,      "mw_p":               mw_p,
    "n_genes":            len(ase_bg), "n_gxe":            int(ase_bg["is_gxe"].sum()),
    "n_sig_bins":         n_sig,    "n_testable_bins":    len(bin_stouffer),
}]).to_csv(f"{RESULTS}/gbs_enrichment_results.tsv", sep="\t", index=False)

trait_df.to_csv(f"{RESULTS}/gbs_enrichment_by_trait.tsv", sep="\t", index=False)
print("   Saved gbs_bin_stouffer_stats.tsv")
print("   Saved gbs_gene_bin_assignments.tsv")
print("   Saved gbs_enrichment_results.tsv")
print("   Saved gbs_enrichment_by_trait.tsv")

# ── 9. Figures ────────────────────────────────────────────────────────────────
print("\nGenerating figures ...", flush=True)
COLORS = {"gxe": "#e74c3c", "bg": "#3498db", "neutral": "#95a5a6"}

fig, axes = plt.subplots(2, 3, figsize=(18, 11))
fig.suptitle("GBS Recombination Bin Enrichment\n"
             "(IBM GxE-ASE genes vs GxE phenotypic QTL)",
             fontsize=13, fontweight="bold")

# A — Stouffer p-value distribution
ax = axes[0, 0]
ax.hist(bin_stouffer["stouffer_GEp"], bins=50, color=COLORS["bg"], edgecolor="black", lw=0.4)
ax.axvline(SIG_P_THRESH, color="red", lw=1.5, ls="--", label=f"p = {SIG_P_THRESH}")
ax.set_xlabel("Stouffer combined GEp"); ax.set_ylabel("Bins")
ax.set_title(f"(A) QTL signal per bin\n({len(bin_stouffer):,} bins, {len(KEEP_TRAITS)} traits)")
ax.text(0.97, 0.97, f"raw p < {SIG_P_THRESH}: {n_sig:,}", transform=ax.transAxes,
        ha="right", va="top", fontsize=9)
ax.legend(fontsize=8)

# B — Gene-bin assignments by chromosome
ax = axes[0, 1]
chr_c = ase_bg.groupby("chr").agg(total=("GeneID","count"),
                                   sig=("bin_sig","sum")).reset_index()
x = np.arange(len(chr_c)); w = 0.4
ax.bar(x - w/2, chr_c["total"], w, label="All bg genes", color=COLORS["neutral"],
       edgecolor="black", lw=0.4)
ax.bar(x + w/2, chr_c["sig"],   w, label="In sig QTL bin", color=COLORS["gxe"],
       edgecolor="black", lw=0.4)
ax.set_xticks(x); ax.set_xticklabels([f"Chr{c}" for c in chr_c["chr"]],
                                       rotation=45, fontsize=8)
ax.set_ylabel("Genes"); ax.set_title("(B) Gene-bin assignments by chromosome")
ax.legend(fontsize=8)

# C — Binary enrichment bar
ax = axes[0, 2]
pct_gxe = a / (a + b) * 100 if (a + b) > 0 else 0
pct_non = c / (c + d) * 100 if (c + d) > 0 else 0
bars = ax.bar(["GxE-ASE", "Non-GxE"], [pct_gxe, pct_non],
              color=[COLORS["gxe"], COLORS["bg"]], edgecolor="black", lw=0.5, width=0.5)
for bar, pct, n in zip(bars, [pct_gxe, pct_non], [a+b, c+d]):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
            f"{pct:.1f}%\n(n={n:,})", ha="center", va="bottom", fontsize=9)
ax.set_ylabel("% genes in significant QTL bin")
ax.set_title(f"(C) Binary enrichment\nOR={OR_obs:.3f}  perm_p={perm_p:.4f}")

# D — ASE signal by QTL bin quartile
ax = axes[1, 0]
cont_plot = cont_df.copy()
cont_plot["qtl_q"] = pd.qcut(cont_plot["qtl_signal"], q=4,
                               labels=["Q1\n(weak)", "Q2", "Q3", "Q4\n(strong)"],
                               duplicates="drop")
bp = ax.boxplot([cont_plot.loc[cont_plot["qtl_q"] == q, "ase_signal"].values
                 for q in cont_plot["qtl_q"].cat.categories],
                labels=cont_plot["qtl_q"].cat.categories,
                patch_artist=True, notch=False, showfliers=False,
                medianprops=dict(color="red", lw=2))
for patch, col in zip(bp["boxes"], ["#d6eaf8","#85c1e9","#2e86c1","#1a5276"]):
    patch.set_facecolor(col)
ax.set_xlabel("QTL signal quartile"); ax.set_ylabel("-log10(ASE p-value)")
ax.set_title(f"(D) ASE signal by QTL quartile\nSpearman r={r_sp:.4f}  perm_p={perm_p_sp:.4f}")

# E — Partial Spearman residual scatter
ax = axes[1, 1]
ax.hexbin(ex, ey, gridsize=50, cmap="Blues", mincnt=1)
ax.set_xlabel("QTL signal residual\n(baseMean partialled out)")
ax.set_ylabel("ASE signal residual")
ax.set_title(f"(E) Partial Spearman\nr={r_partial:.4f}  p={p_partial:.4f}")
if len(ex) > 10:
    z = np.polyfit(ex, ey, 1)
    xl = np.linspace(ex.min(), ex.max(), 100)
    ax.plot(xl, np.polyval(z, xl), "r-", lw=1.5)

# F — Trait-stratified ORs
ax = axes[1, 2]
td = trait_df.dropna(subset=["OR"]).sort_values("OR", ascending=True)
cols_t = [COLORS["gxe"] if p < 0.05 else COLORS["neutral"]
          for p in td["perm_p"].fillna(1)]
ax.barh(range(len(td)), td["OR"], color=cols_t, edgecolor="black", lw=0.4)
ax.axvline(1, color="black", lw=1, ls="--")
ax.set_yticks(range(len(td)))
ax.set_yticklabels(td["trait"].str.replace("_", "\n", regex=False), fontsize=7)
ax.set_xlabel("Odds Ratio")
ax.set_title("(F) Trait-stratified ORs\n(red = perm_p < 0.05)")

fig.tight_layout()
fig_path = f"{FIGS}/gbs_enrichment_figure.png"
fig.savefig(fig_path, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"   Saved {fig_path}")

# ── 10. Summary ───────────────────────────────────────────────────────────────
print("\n" + "="*65, flush=True)
print("FINAL SUMMARY", flush=True)
print(f"  Testable bins:        {len(bin_stouffer):,}")
print(f"  QTL bins with raw p < {SIG_P_THRESH}: {n_sig:,} (naive threshold, not BH-FDR)")
print(f"  Background genes in bin: {len(ase_bg):,}  (GxE-ASE: {int(ase_bg['is_gxe'].sum())})")
print(f"\n  BINARY:     OR={OR_obs:.4f}  Fisher p={p_fisher:.4f}  perm_p={perm_p:.4f}")
print(f"  CONTINUOUS: Spearman r={r_sp:.4f} (p={p_sp:.4f}, perm_p={perm_p_sp:.4f})")
print(f"              Partial Spearman r={r_partial:.4f} (p={p_partial:.4f})")
print(f"              Mann-Whitney rbc={rbc:.4f}  p={mw_p:.4f}")
sig_traits = trait_df[trait_df["perm_p"].fillna(1) < 0.05]
print(f"\n  TRAIT STRATIFIED — {len(sig_traits)} trait(s) with perm_p < 0.05:")
for _, row in sig_traits.iterrows():
    print(f"    {row['trait']:<35}  OR={row['OR']:.3f}  perm_p={row['perm_p']:.4f}")
if len(sig_traits) == 0:
    print("    None")
print("\nDONE", flush=True)
