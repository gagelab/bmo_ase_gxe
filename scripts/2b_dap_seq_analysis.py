#!/usr/bin/env python3
"""
2b_dap_seq_analysis.py

Generates three DAP-seq derived result files from raw peak files + annotation:

  1. results/dap_seq_per_gene.tsv
       Per-gene peak counts in ±500 bp promoter windows (one row per universe gene).
       Columns: chr, prom_s, prom_e, strand, gene, is_gxe, is_bg,
                n_b73_spec, n_shared, n_mo17_spec, n_diff, n_total, frac_diff, gene_class

  2. results/dap_seq_per_tf.tsv
       Per-TF Fisher enrichment of GxE vs background genes.
       Columns: TF, n_peaks_b73spec, n_peaks_mo17spec, n_peaks_shared,
                n_gxe_diff, n_bg_diff, n_gxe_any, n_bg_any, n_gxe_b73spec, n_gxe_mo17spec,
                OR_diff, p_diff, OR_any, p_any, padj_diff, padj_any,
                log2OR_diff, log2OR_any, [key, v5_geneID, tf_family for annotated TFs]

  3. results/indel_dap_overlap_per_gene.tsv
       Per-gene counts of small indels (< 50 bp) overlapping DAP peaks.
       Columns: gene, is_gxe, is_bg, n_b73spec, n_mo17spec, n_shared,
                n_no_peak, n_genotype_spec, n_total

  This script replaces the inline Step 3 flag-update code in RERUN_WORKFLOW.md.
  Run after 1_deseq_analysis.R (which produces GxE_gene_IDs.txt / background_gene_IDs.txt).

  Runtime: ~5-10 min (peak-file loading dominates).

Reads:
  data/annotation.gff
  data/GxE_gene_IDs.txt
  data/background_gene_IDs.txt
  data/tf_annotation.tsv   (stable TF → v5_geneID, tf_family lookup for expressed TFs)
  data/epigenome/dap_seq/normalized_specific_and_shared_peaks/.../*_withcoords.tsv
  data/Mo17_toB73v5_paf_syri_noStartPOS0_INDELs_less50bp.vcf
"""

import os
import re
import time
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

t0 = time.time()

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE      = "./"
DATA      = f"{BASE}/data"
RESULTS   = f"{BASE}/results"
B73_DIR   = (f"{DATA}/epigenome/dap_seq/normalized_specific_and_shared_peaks/"
             "B73v5_MATCH_Mo17-B73v5_specific_B73v5_shared")
MO17_DIR  = (f"{DATA}/epigenome/dap_seq/normalized_specific_and_shared_peaks/"
             "Mo17_MATCH_B73v5-Mo17_specific_Mo17_shared")
GFF        = f"{DATA}/annotation.gff"
VCF_SMALL  = f"{DATA}/Mo17_toB73v5_paf_syri_noStartPOS0_INDELs_less50bp.vcf"
PROMO_WIN  = 500   # bp on each side of TSS (0-based coords)

# ── 1. Load gene sets ─────────────────────────────────────────────────────────
print("1. Loading gene sets...", flush=True)
gxe_genes = set(open(f"{DATA}/GxE_gene_IDs.txt").read().split())
bg_genes  = set(open(f"{DATA}/background_gene_IDs.txt").read().split())
universe  = gxe_genes | bg_genes
n_gxe, n_bg = len(gxe_genes), len(bg_genes)
print(f"   GxE: {n_gxe:,}  Background: {n_bg:,}  Universe: {len(universe):,}")

# ── 2. Parse GFF → promoter windows ──────────────────────────────────────────
print(f"\n2. Parsing promoter windows from GFF (±{PROMO_WIN} bp)...", flush=True)
rows = []
with open(GFF) as fh:
    for line in fh:
        if line.startswith('#') or '\tgene\t' not in line:
            continue
        p = line.rstrip().split('\t')
        if len(p) < 9 or p[2] != 'gene':
            continue
        m = re.search(r'ID=(Zm\d+eb\d+)', p[8])
        if not m:
            continue
        gid    = m.group(1)
        if gid not in universe:
            continue
        strand = p[6]
        s = int(p[3]) - 1   # GFF is 1-based; convert to 0-based
        e = int(p[4])        # 0-based exclusive
        tss    = s if strand == '+' else e
        prom_s = max(0, tss - PROMO_WIN)
        prom_e = tss + PROMO_WIN
        rows.append({'chr': p[0], 'prom_s': prom_s, 'prom_e': prom_e,
                     'strand': strand, 'gene': gid})

proms = pd.DataFrame(rows).drop_duplicates('gene').reset_index(drop=True)
print(f"   {len(proms):,} promoters parsed  "
      f"({proms['gene'].isin(gxe_genes).sum()} GxE, "
      f"{proms['gene'].isin(bg_genes).sum()} BG)")

# Index promoters by chromosome for fast overlap lookup
_prom_idx = {}
for chrom, grp in proms.groupby('chr'):
    _prom_idx[chrom] = {
        'row':    grp.index.values,
        'prom_s': grp['prom_s'].values,
        'prom_e': grp['prom_e'].values,
        'is_gxe': grp['gene'].isin(gxe_genes).values,
        'is_bg':  grp['gene'].isin(bg_genes).values,
    }

# ── 3. Load all TF peak files ─────────────────────────────────────────────────
def _parse_tf(fname, genotype_token):
    """Extract TF name from a peak filename."""
    m = re.match(rf'^(.+?)_{genotype_token}_\d+_MATCH', fname)
    return m.group(1) if m else None

def _load_narrowpeak(fpath, col_chr=0, col_start=1, col_end=2):
    """
    Load a _withcoords.tsv peak file.  Returns dict: chr → (sorted start array, end array).

    B73-specific and shared peak files:  use cols 0,1,2 (B73 peak coordinates).
    Mo17-specific peak files:            use cols 10,11,12 (B73 gene coordinates assigned
                                         by bedtools closest; Mo17 peak coords in cols 0-2
                                         are in Mo17 genome space and cannot be intersected
                                         directly with B73 promoter windows).
    This matches the approach in scripts/6_tf_gxe_mechanism.py (load_peaks_from_dir).
    """
    peaks_by_chr = {}
    buf = {}
    min_cols = max(col_chr, col_start, col_end) + 1
    with open(fpath) as fh:
        for line in fh:
            p = line.split('\t')
            if len(p) < min_cols:
                continue
            try:
                chrom = p[col_chr].strip()
                s = int(p[col_start])
                e = int(p[col_end])
            except (ValueError, IndexError):
                continue
            if e <= s:
                continue
            buf.setdefault(chrom, []).append((s, e))
    for chrom, pairs in buf.items():
        arr = np.array(sorted(pairs), dtype=np.int64)
        peaks_by_chr[chrom] = (arr[:, 0], arr[:, 1])
    return peaks_by_chr

def _count_overlaps(peaks_by_chr, prom_idx, n_proms):
    """
    For each gene promoter, count the number of peaks overlapping it.
    Uses binary search on sorted peak starts for efficiency.
    Returns an int32 array of length n_proms.
    """
    counts = np.zeros(n_proms, dtype=np.int32)
    for chrom, pc in prom_idx.items():
        if chrom not in peaks_by_chr:
            continue
        pa_s, pa_e = peaks_by_chr[chrom]
        for i in range(len(pc['row'])):
            ps, pe, row = pc['prom_s'][i], pc['prom_e'][i], pc['row'][i]
            # peaks with start < pe (upper bound)
            right = int(np.searchsorted(pa_s, pe, side='left'))
            if right > 0:
                counts[row] += int(np.sum(pa_e[:right] > ps))
    return counts

def _has_hit(peaks_by_chr, prom_idx, n_proms):
    """
    Boolean array: True if gene's promoter overlaps ≥1 peak.
    """
    hit = np.zeros(n_proms, dtype=bool)
    for chrom, pc in prom_idx.items():
        if chrom not in peaks_by_chr:
            continue
        pa_s, pa_e = peaks_by_chr[chrom]
        for i in range(len(pc['row'])):
            ps, pe, row = pc['prom_s'][i], pc['prom_e'][i], pc['row'][i]
            right = int(np.searchsorted(pa_s, pe, side='left'))
            if right > 0 and np.any(pa_e[:right] > ps):
                hit[row] = True
    return hit

# Accumulate counts across all TF files
print(f"\n3. Loading DAP-seq peak files and counting promoter overlaps...", flush=True)

# Accumulators for per-gene totals
n_b73_spec  = np.zeros(len(proms), dtype=np.int32)
n_mo17_spec = np.zeros(len(proms), dtype=np.int32)
n_shared    = np.zeros(len(proms), dtype=np.int32)

# Per-TF records
tf_records = {}   # TF → dict of arrays

b73_files  = [f for f in os.listdir(B73_DIR)  if f.endswith('_specific_withcoords.tsv')]
mo17_files = [f for f in os.listdir(MO17_DIR) if f.endswith('_specific_withcoords.tsv')]
b73_shared_files = [f for f in os.listdir(B73_DIR) if f.endswith('_shared_withcoords.tsv')]

# Build TF → file mapping for all three types
b73_by_tf   = {}
mo17_by_tf  = {}
shared_by_tf = {}

for fname in b73_files:
    tf = _parse_tf(fname, 'B73')
    if tf:
        b73_by_tf[tf] = os.path.join(B73_DIR, fname)

for fname in mo17_files:
    tf = _parse_tf(fname, 'Mo17')
    if tf:
        mo17_by_tf[tf] = os.path.join(MO17_DIR, fname)

for fname in b73_shared_files:
    tf = _parse_tf(fname, 'B73')
    if tf:
        shared_by_tf[tf] = os.path.join(B73_DIR, fname)

all_tfs = sorted(set(b73_by_tf) | set(mo17_by_tf))
print(f"   {len(all_tfs)} TFs with peak data", flush=True)

n_p = len(proms)
is_gxe_arr = proms['gene'].isin(gxe_genes).values
is_bg_arr  = proms['gene'].isin(bg_genes).values

# Storage for indel-overlap: union peak sets across all TFs
_all_b73spec_by_chr  = {}
_all_mo17spec_by_chr = {}
_all_shared_by_chr   = {}

for i_tf, tf in enumerate(all_tfs):
    if (i_tf + 1) % 50 == 0:
        print(f"   {i_tf+1}/{len(all_tfs)} TFs  ({time.time()-t0:.0f}s)", flush=True)

    # B73-specific and shared peaks: cols 0,1,2 are B73 genome coordinates.
    # Mo17-specific peaks: cols 10,11,12 are bedtools-assigned B73 gene locus
    # coordinates (cols 0-2 are Mo17 genome coordinates, not usable with B73 GFF).
    peaks_b73  = _load_narrowpeak(b73_by_tf[tf])   if tf in b73_by_tf   else {}
    peaks_mo17 = _load_narrowpeak(mo17_by_tf[tf], col_chr=10, col_start=11, col_end=12) \
                 if tf in mo17_by_tf  else {}
    peaks_sh   = _load_narrowpeak(shared_by_tf[tf]) if tf in shared_by_tf else {}

    # Per-gene totals
    cnt_b73  = _count_overlaps(peaks_b73,  _prom_idx, n_p)
    cnt_mo17 = _count_overlaps(peaks_mo17, _prom_idx, n_p)
    cnt_sh   = _count_overlaps(peaks_sh,   _prom_idx, n_p)

    n_b73_spec  += cnt_b73
    n_mo17_spec += cnt_mo17
    n_shared    += cnt_sh

    # Per-TF enrichment booleans
    hit_b73spec  = (cnt_b73 > 0)
    hit_mo17spec = (cnt_mo17 > 0)
    hit_diff = hit_b73spec | hit_mo17spec
    hit_any  = hit_diff | (cnt_sh > 0)

    n_gxe_b73spec  = int(hit_b73spec[is_gxe_arr].sum())
    n_gxe_mo17spec = int(hit_mo17spec[is_gxe_arr].sum())
    n_gxe_diff = int(hit_diff[is_gxe_arr].sum())
    n_bg_diff  = int(hit_diff[is_bg_arr].sum())
    n_gxe_any  = int(hit_any[is_gxe_arr].sum())
    n_bg_any   = int(hit_any[is_bg_arr].sum())

    # n_peaks_* = total peaks in this file (across all chr)
    def _total_peaks(d):
        return sum(len(v[0]) for v in d.values())

    tf_records[tf] = {
        'n_peaks_b73spec':  _total_peaks(peaks_b73),
        'n_peaks_mo17spec': _total_peaks(peaks_mo17),
        'n_peaks_shared':   _total_peaks(peaks_sh),
        'n_gxe_diff':   n_gxe_diff,
        'n_bg_diff':    n_bg_diff,
        'n_gxe_any':    n_gxe_any,
        'n_bg_any':     n_bg_any,
        'n_gxe_b73spec': n_gxe_b73spec,
        'n_gxe_mo17spec': n_gxe_mo17spec,
    }

    # Accumulate union peak sets for indel-overlap computation
    for chrom, (ps, pe) in peaks_b73.items():
        if chrom not in _all_b73spec_by_chr:
            _all_b73spec_by_chr[chrom] = ([], [])
        _all_b73spec_by_chr[chrom][0].extend(ps.tolist())
        _all_b73spec_by_chr[chrom][1].extend(pe.tolist())

    for chrom, (ps, pe) in peaks_mo17.items():
        if chrom not in _all_mo17spec_by_chr:
            _all_mo17spec_by_chr[chrom] = ([], [])
        _all_mo17spec_by_chr[chrom][0].extend(ps.tolist())
        _all_mo17spec_by_chr[chrom][1].extend(pe.tolist())

    for chrom, (ps, pe) in peaks_sh.items():
        if chrom not in _all_shared_by_chr:
            _all_shared_by_chr[chrom] = ([], [])
        _all_shared_by_chr[chrom][0].extend(ps.tolist())
        _all_shared_by_chr[chrom][1].extend(pe.tolist())

print(f"   All TF files loaded  ({time.time()-t0:.0f}s)", flush=True)

# Sort union peak arrays
def _sort_peaks(d):
    out = {}
    for chrom, (sl, el) in d.items():
        arr = np.array(sorted(zip(sl, el)), dtype=np.int64)
        out[chrom] = (arr[:, 0], arr[:, 1])
    return out

_all_b73spec_by_chr  = _sort_peaks(_all_b73spec_by_chr)
_all_mo17spec_by_chr = _sort_peaks(_all_mo17spec_by_chr)
_all_shared_by_chr   = _sort_peaks(_all_shared_by_chr)

# ── 4. Build and write dap_seq_per_gene.tsv ───────────────────────────────────
print("\n4. Writing results/dap_seq_per_gene.tsv...", flush=True)

n_diff  = n_b73_spec + n_mo17_spec
n_total = n_diff + n_shared
frac_diff = np.where(n_total > 0, n_diff / n_total.astype(float), np.nan)

proms['is_gxe'] = is_gxe_arr
proms['is_bg']  = is_bg_arr
proms['n_b73_spec']  = n_b73_spec
proms['n_shared']    = n_shared
proms['n_mo17_spec'] = n_mo17_spec
proms['n_diff']      = n_diff
proms['n_total']     = n_total
proms['frac_diff']   = frac_diff
proms['gene_class']  = 'Other'
proms.loc[proms['is_bg'],  'gene_class'] = 'Background'
proms.loc[proms['is_gxe'], 'gene_class'] = 'GxE-ASE'

col_order = ['chr', 'prom_s', 'prom_e', 'strand', 'gene', 'is_gxe', 'is_bg',
             'n_b73_spec', 'n_shared', 'n_mo17_spec', 'n_diff', 'n_total',
             'frac_diff', 'gene_class']
proms[col_order].to_csv(f"{RESULTS}/dap_seq_per_gene.tsv", sep='\t', index=False)
print(f"   Saved dap_seq_per_gene.tsv ({len(proms):,} genes)")
print(f"   GxE genes with any DAP peak: "
      f"{((n_total > 0) & is_gxe_arr).sum()}/{n_gxe} "
      f"({100*((n_total > 0) & is_gxe_arr).sum()/n_gxe:.1f}%)")
print(f"   BG genes with any DAP peak: "
      f"{((n_total > 0) & is_bg_arr).sum()}/{n_bg} "
      f"({100*((n_total > 0) & is_bg_arr).sum()/n_bg:.1f}%)")

# ── 5. Build and write dap_seq_per_tf.tsv ────────────────────────────────────
print("\n5. Writing results/dap_seq_per_tf.tsv...", flush=True)

def _fisher_or(a, b, n_a, n_b):
    """OR and p-value from 2×2 table. a and b are successes, n_a/n_b are totals."""
    table = [[a, n_a - a], [b, n_b - b]]
    or_, p = fisher_exact(table, alternative='greater')
    return or_, p

tf_rows = []
for tf, rec in tf_records.items():
    or_diff, p_diff = _fisher_or(rec['n_gxe_diff'], rec['n_bg_diff'], n_gxe, n_bg)
    or_any,  p_any  = _fisher_or(rec['n_gxe_any'],  rec['n_bg_any'],  n_gxe, n_bg)
    tf_rows.append({
        'TF':              tf,
        'n_peaks_b73spec': rec['n_peaks_b73spec'],
        'n_peaks_mo17spec':rec['n_peaks_mo17spec'],
        'n_peaks_shared':  rec['n_peaks_shared'],
        'n_gxe_diff':   rec['n_gxe_diff'],
        'n_bg_diff':    rec['n_bg_diff'],
        'n_gxe_any':    rec['n_gxe_any'],
        'n_bg_any':     rec['n_bg_any'],
        'n_gxe_b73spec': rec['n_gxe_b73spec'],
        'n_gxe_mo17spec':rec['n_gxe_mo17spec'],
        'OR_diff':  or_diff,
        'p_diff':   p_diff,
        'OR_any':   or_any,
        'p_any':    p_any,
    })

tf_df = pd.DataFrame(tf_rows).sort_values('p_diff').reset_index(drop=True)

# BH correction
_, padj_diff, _, _ = multipletests(tf_df['p_diff'], method='fdr_bh')
_, padj_any,  _, _ = multipletests(tf_df['p_any'],  method='fdr_bh')
tf_df['padj_diff'] = padj_diff
tf_df['padj_any']  = padj_any
tf_df['log2OR_diff'] = np.log2(tf_df['OR_diff'].replace(0, np.nan))
tf_df['log2OR_any']  = np.log2(tf_df['OR_any'].replace(0, np.nan))

# Join with TF annotation for expressed TFs
tf_annot = pd.read_csv(f"{DATA}/tf_annotation.tsv", sep='\t')
tf_df = tf_df.merge(tf_annot, on='TF', how='left')

tf_df.to_csv(f"{RESULTS}/dap_seq_per_tf.tsv", sep='\t', index=False)
print(f"   Saved dap_seq_per_tf.tsv ({len(tf_df)} TFs)")
print(f"   Top TF by OR_diff (genotype-specific binding enrichment in GxE):")
top = tf_df.nsmallest(3, 'p_diff')[['TF', 'OR_diff', 'p_diff', 'padj_diff']]
for _, r in top.iterrows():
    print(f"     {r['TF']:20s}  OR={r['OR_diff']:.3f}  p={r['p_diff']:.4g}  padj={r['padj_diff']:.4g}")

# ── 6. Indel-overlap per gene ─────────────────────────────────────────────────
print(f"\n6. Computing indel × DAP-peak overlaps...", flush=True)
print(f"   Loading VCF: {Path(VCF_SMALL).name} ...", flush=True)

indels_by_chr = {}
with open(VCF_SMALL) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        p = line.split('\t', 5)
        chrom = p[0]
        if not chrom.startswith('chr'):
            chrom = 'chr' + chrom
        pos = int(p[1])
        indels_by_chr.setdefault(chrom, []).append(pos)

# Sort and convert to arrays
for chrom in indels_by_chr:
    indels_by_chr[chrom] = np.array(sorted(indels_by_chr[chrom]), dtype=np.int64)

n_total_indels = sum(len(v) for v in indels_by_chr.values())
print(f"   Loaded {n_total_indels:,} small indels across {len(indels_by_chr)} chromosomes")

def _count_peak_hits_per_position(positions, peaks_s, peaks_e):
    """
    For each position in `positions`, count how many peaks (pa_s[i], pa_e[i]) contain it.
    Position pos is contained if pa_s[i] <= pos < pa_e[i].
    Returns array of counts, same length as positions.
    """
    out = np.zeros(len(positions), dtype=np.int32)
    for j, pos in enumerate(positions):
        right = int(np.searchsorted(peaks_s, pos + 1, side='left'))
        if right > 0:
            out[j] = int(np.sum(peaks_e[:right] > pos))
    return out

# Per-gene accumulators
ig_n_b73spec  = np.zeros(n_p, dtype=np.int32)
ig_n_mo17spec = np.zeros(n_p, dtype=np.int32)
ig_n_shared   = np.zeros(n_p, dtype=np.int32)
ig_n_no_peak  = np.zeros(n_p, dtype=np.int32)

for chrom, pc in _prom_idx.items():
    if chrom not in indels_by_chr:
        continue
    all_pos = indels_by_chr[chrom]

    pb73  = _all_b73spec_by_chr.get(chrom,  (np.array([], np.int64), np.array([], np.int64)))
    pmo17 = _all_mo17spec_by_chr.get(chrom, (np.array([], np.int64), np.array([], np.int64)))
    psh   = _all_shared_by_chr.get(chrom,   (np.array([], np.int64), np.array([], np.int64)))

    for i in range(len(pc['row'])):
        ps, pe, row = pc['prom_s'][i], pc['prom_e'][i], pc['row'][i]

        # Indels in this promoter
        lo = int(np.searchsorted(all_pos, ps,   side='left'))
        hi = int(np.searchsorted(all_pos, pe,   side='left'))
        if lo >= hi:
            continue
        pos_in_prom = all_pos[lo:hi]

        # Count overlapping peaks per type per indel position
        hits_b73  = _count_peak_hits_per_position(pos_in_prom, pb73[0],  pb73[1])
        hits_mo17 = _count_peak_hits_per_position(pos_in_prom, pmo17[0], pmo17[1])
        hits_sh   = _count_peak_hits_per_position(pos_in_prom, psh[0],   psh[1])

        total_hits = hits_b73 + hits_mo17 + hits_sh
        ig_n_b73spec[row]  += int(hits_b73.sum())
        ig_n_mo17spec[row] += int(hits_mo17.sum())
        ig_n_shared[row]   += int(hits_sh.sum())
        ig_n_no_peak[row]  += int(np.sum(total_hits == 0))

ig_n_genotype_spec = ig_n_b73spec + ig_n_mo17spec
ig_n_total         = ig_n_b73spec + ig_n_mo17spec + ig_n_shared

olap_df = pd.DataFrame({
    'gene':           proms['gene'].values,
    'is_gxe':         is_gxe_arr,
    'is_bg':          is_bg_arr,
    'n_b73spec':      ig_n_b73spec,
    'n_mo17spec':     ig_n_mo17spec,
    'n_shared':       ig_n_shared,
    'n_no_peak':      ig_n_no_peak,
    'n_genotype_spec':ig_n_genotype_spec,
    'n_total':        ig_n_total,
})
olap_df.to_csv(f"{RESULTS}/indel_dap_overlap_per_gene.tsv", sep='\t', index=False)
print(f"   Saved indel_dap_overlap_per_gene.tsv ({len(olap_df):,} genes)")

# Quick enrichment check
g = olap_df[olap_df['is_gxe']]
b = olap_df[olap_df['is_bg']]
n_g_hit = (g['n_genotype_spec'] > 0).sum()
n_b_hit = (b['n_genotype_spec'] > 0).sum()
from scipy.stats import fisher_exact as fe
or_id, p_id = fe([[n_g_hit, len(g)-n_g_hit], [n_b_hit, len(b)-n_b_hit]])
print(f"   Indel×DAP enrichment: GxE {n_g_hit}/{len(g)} "
      f"({100*n_g_hit/len(g):.1f}%) vs BG {n_b_hit}/{len(b)} "
      f"({100*n_b_hit/len(b):.1f}%)  OR={or_id:.3f}  Fisher p={p_id:.4g}")

# ── Summary ───────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print(f"Done in {(time.time()-t0)/60:.1f} min")
print(f"Outputs:")
print(f"  {RESULTS}/dap_seq_per_gene.tsv")
print(f"  {RESULTS}/dap_seq_per_tf.tsv")
print(f"  {RESULTS}/indel_dap_overlap_per_gene.tsv")
print(f"\nVerification (expected from original analysis):")
print(f"  dap_seq_per_gene.tsv: ~11,637 rows")
print(f"  dap_seq_per_tf.tsv: 197-198 TFs")
print(f"  Indel×DAP OR ≈ 1.52, Fisher p ≈ 0.007")
