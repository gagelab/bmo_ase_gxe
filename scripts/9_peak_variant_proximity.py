"""
10_peak_variant_proximity.py
-----------------------------
Test whether genotype-specific DAP-seq peaks sit in more polymorphic
sequence than shared peaks, within gene promoter regions.

Hypothesis (from prior work): allele-specific TF binding is enriched
for nearby sequence variation. The prediction is that B73-specific and
Mo17-specific peaks have more Mo17/B73 variants within a flanking window
than shared peaks do.

For each peak that overlaps a gene's ±500 bp promoter window, count all
Mo17 variants (SNPs + small indels + large indels) within ±FLANK bp of
the peak boundaries. Compare counts between peak types using Mann-Whitney U.

Tests performed:
  1. Overall (all genes): genotype-specific vs shared, B73-specific vs
     shared, Mo17-specific vs shared.
  2. Split by gene type (GxE-ASE vs background).
  3. Interaction: is the genotype-specific enrichment larger in GxE genes
     than in background genes?

Reads:
  results/dap_seq_per_gene.tsv
  data/GxE_gene_IDs.txt
  data/background_gene_IDs.txt
  data/Mo17_toB73v5_paf_syri_noStartPOS0_SNPs.vcf
  data/Mo17_toB73v5_paf_syri_noStartPOS0_INDELs_less50bp.vcf
  data/Mo17_toB73v5_paf_syri_noStartPOS0_INDELs_more50bp.vcf
  data/epigenome/dap_seq/.../B73v5_.../*specific.narrowPeak
  data/epigenome/dap_seq/.../B73v5_.../*shared.narrowPeak
  data/epigenome/dap_seq/.../Mo17_.../*specific.narrowPeak

Writes: nothing (results printed to stdout only)

Runtime: ~3-4 min (VCF loading dominates)
"""

import glob, os, re
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

BASE    = "./"
DATA    = f"{BASE}/data"
RESULTS = f"{BASE}/results"

B73_DIR  = (f"{DATA}/epigenome/dap_seq/normalized_specific_and_shared_peaks/"
            "B73v5_MATCH_Mo17-B73v5_specific_B73v5_shared")
MO17_DIR = (f"{DATA}/epigenome/dap_seq/normalized_specific_and_shared_peaks/"
            "Mo17_MATCH_B73v5-Mo17_specific_Mo17_shared")

DAP_TSV  = f"{RESULTS}/dap_seq_per_gene.tsv"
GXE_FILE = f"{DATA}/GxE_gene_IDs.txt"
BG_FILE  = f"{DATA}/background_gene_IDs.txt"

VCF_FILES = {
    "snp":         f"{DATA}/Mo17_toB73v5_paf_syri_noStartPOS0_SNPs.vcf",
    "small_indel": f"{DATA}/Mo17_toB73v5_paf_syri_noStartPOS0_INDELs_less50bp.vcf",
    "large_indel": f"{DATA}/Mo17_toB73v5_paf_syri_noStartPOS0_INDELs_more50bp.vcf",
}

FLANK = 50    # bp on each side of the peak

# ── 1. Load all variants → per-chrom sorted position arrays ──────────────────
print("Loading variants ...", flush=True)
var_buf = defaultdict(list)
for vtype, path in VCF_FILES.items():
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 3)
            if len(parts) < 2:
                continue
            chrom = parts[0].strip()
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            try:
                var_buf[chrom].append(int(parts[1]) - 1)   # convert to 0-based
                n += 1
            except (ValueError, IndexError):
                continue
    print(f"  {vtype}: {n:,}")

var_by_chr = {c: np.array(sorted(v), dtype=np.int64) for c, v in var_buf.items()}
print(f"  Total: {sum(len(v) for v in var_by_chr.values()):,}")

def count_variants_near(chrom, peak_s, peak_e):
    """Count variants in [peak_s - FLANK, peak_e + FLANK]."""
    if chrom not in var_by_chr:
        return 0
    arr = var_by_chr[chrom]
    return int(np.searchsorted(arr, peak_e + FLANK, side="right") -
               np.searchsorted(arr, peak_s - FLANK, side="left"))

# ── 2. Load peaks ─────────────────────────────────────────────────────────────
print("\nLoading DAP-seq peaks ...", flush=True)

def _tf_b73(fname):
    m = re.search(r"^(.+?)_B73_", os.path.basename(fname))
    return m.group(1) if m else None

def _tf_mo17(fname):
    m = re.search(r"^(.+?)_Mo17_", os.path.basename(fname))
    return m.group(1) if m else None

def load_peaks(peak_dir, tf_extractor, pattern):
    """Return list of (chrom, start, end, tf)."""
    peaks = []
    for fpath in glob.glob(os.path.join(peak_dir, pattern)):
        tf = tf_extractor(fpath)
        if not tf:
            continue
        with open(fpath) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                try:
                    peaks.append((parts[0], int(parts[1]), int(parts[2]), tf))
                except (ValueError, IndexError):
                    continue
    return peaks

b73spec_peaks  = load_peaks(B73_DIR,  _tf_b73,  "*specific.narrowPeak")
mo17spec_peaks = load_peaks(MO17_DIR, _tf_mo17, "*specific.narrowPeak")
shared_peaks   = load_peaks(B73_DIR,  _tf_b73,  "*shared.narrowPeak")

print(f"  B73-specific peaks:  {len(b73spec_peaks):,}")
print(f"  Mo17-specific peaks: {len(mo17spec_peaks):,}")
print(f"  Shared peaks:        {len(shared_peaks):,}")

# ── 3. Load gene info and build promoter index ────────────────────────────────
print("\nLoading gene info ...", flush=True)
dap     = pd.read_csv(DAP_TSV, sep="\t")
gxe_ids = set(open(GXE_FILE).read().split())
bg_ids  = set(open(BG_FILE).read().split()) - gxe_ids

# Per-chromosome sorted list of (prom_s, prom_e, gene, is_gxe)
prom_by_chr = defaultdict(list)
for _, row in dap[dap["is_gxe"] | dap["is_bg"]].iterrows():
    prom_by_chr[row["chr"]].append(
        (int(row["prom_s"]), int(row["prom_e"]),
         row["gene"], row["gene"] in gxe_ids))
for chrom in prom_by_chr:
    prom_by_chr[chrom].sort()

print(f"  GxE genes: {len(gxe_ids):,}   Background genes: {len(bg_ids):,}")

# ── 4. Count nearby variants for each peak-in-promoter instance ───────────────
print("\nScoring peaks ...", flush=True)

def score_peaks(peak_list, peak_type):
    """
    For each peak overlapping a gene promoter, record the number of
    variants within ±FLANK bp of the peak.
    """
    records = []
    by_chr = defaultdict(list)
    for chrom, ps, pe, _ in peak_list:
        by_chr[chrom].append((ps, pe))

    for chrom, peaks in by_chr.items():
        if chrom not in prom_by_chr:
            continue
        proms = prom_by_chr[chrom]
        for ps, pe in peaks:
            for prom_s, prom_e, gene, is_gxe in proms:
                if prom_s >= pe:
                    break    # promoters are sorted; no further overlaps possible
                if prom_e <= ps:
                    continue
                records.append({
                    "peak_type": peak_type,
                    "n_var":     count_variants_near(chrom, ps, pe),
                    "gene":      gene,
                    "is_gxe":    is_gxe,
                })
    return records

recs = []
for peak_list, ptype in [(b73spec_peaks,  "B73-specific"),
                         (mo17spec_peaks, "Mo17-specific"),
                         (shared_peaks,   "Shared")]:
    r = score_peaks(peak_list, ptype)
    print(f"  {ptype}: {len(r):,} peak-in-promoter instances")
    recs.extend(r)

df = pd.DataFrame(recs)

# ── 5. Statistical tests and reporting ───────────────────────────────────────
print("\n" + "=" * 65)
print(f"VARIANTS WITHIN ±{FLANK} bp OF PEAK  |  peak-in-promoter instances")
print("=" * 65)

GENO_TYPES = ["B73-specific", "Mo17-specific"]

def report_group(label, sub):
    geno = sub[sub["peak_type"].isin(GENO_TYPES)]
    shar = sub[sub["peak_type"] == "Shared"]
    b73s = sub[sub["peak_type"] == "B73-specific"]
    m17s = sub[sub["peak_type"] == "Mo17-specific"]

    print(f"\n  {label}  "
          f"(geno-spec n={len(geno):,}, shared n={len(shar):,})")

    for a_lbl, a, b_lbl, b in [
        ("Genotype-specific", geno, "Shared", shar),
        ("B73-specific",      b73s, "Shared", shar),
        ("Mo17-specific",     m17s, "Shared", shar),
    ]:
        if len(a) < 2 or len(b) < 2:
            continue
        _, p_gt = mannwhitneyu(a["n_var"], b["n_var"], alternative="greater")
        _, p_2s = mannwhitneyu(a["n_var"], b["n_var"], alternative="two-sided")
        print(f"    {a_lbl:20s} vs {b_lbl:8s}  "
              f"mean={a['n_var'].mean():.3f} vs {b['n_var'].mean():.3f}  "
              f"MWU p(>)={p_gt:.4g}  p(2-sided)={p_2s:.4g}")

report_group("ALL GENES", df)
report_group("GxE GENES", df[df["is_gxe"]])
report_group("BG  GENES", df[~df["is_gxe"]])

# Interaction: is the geno-spec enrichment larger in GxE vs BG genes?
print(f"\n  INTERACTION — is the enrichment at genotype-specific peaks "
      f"larger in GxE genes than BG genes?")
for ptype in ["B73-specific", "Mo17-specific", "Genotype-specific"]:
    mask = (df["peak_type"].isin(GENO_TYPES)
            if ptype == "Genotype-specific"
            else (df["peak_type"] == ptype))
    gxe_g = df[df["is_gxe"]  & mask]
    bg_g  = df[~df["is_gxe"] & mask]
    gxe_s = df[df["is_gxe"]  & (df["peak_type"] == "Shared")]
    bg_s  = df[~df["is_gxe"] & (df["peak_type"] == "Shared")]

    if len(gxe_g) < 2 or len(bg_g) < 2:
        continue

    d_gxe = gxe_g["n_var"].mean() - gxe_s["n_var"].mean()
    d_bg  = bg_g["n_var"].mean()  - bg_s["n_var"].mean()
    _, p  = mannwhitneyu(gxe_g["n_var"], bg_g["n_var"], alternative="greater")

    print(f"    {ptype:22s}  "
          f"GxE delta={d_gxe:+.3f}  BG delta={d_bg:+.3f}  "
          f"MWU GxE>BG p={p:.4g}")

print("\nDone.", flush=True)
