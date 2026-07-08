"""
indel_pwm_scoring.py
--------------------
Mechanistic test: do small indels at genotype-specific DAP peaks actually
disrupt / create the TF-binding motif for that specific TF?

Approach
--------
For each indel in a gene's promoter that overlaps a B73-specific or Mo17-specific
DAP peak, we:
  1. Fetch B73 genome context around the indel.
  2. Build REF sequence (B73) and ALT sequence (Mo17 with indel applied).
  3. Score REF with the B73 log-odds PWM and ALT with the Mo17 log-odds PWM
     (max over all positions, both strands).
  4. Δscore = max_score(REF, B73_PWM) − max_score(ALT, Mo17_PWM)
     • B73-specific peaks → expect Δ > 0   (B73 has better motif match)
     • Mo17-specific peaks → expect Δ < 0  (Mo17 has better motif match)

Tests
-----
  • Binomial: fraction of B73-spec records with Δ>0 (vs 0.5 null)
  • Binomial: fraction of Mo17-spec records with Δ<0 (vs 0.5 null)
  • Mann-Whitney: |Δscore| GxE-ASE vs background (are GxE genes more disrupted?)
  • Spearman: |lFC| vs |Δscore| among GxE-ASE genes
"""

import os, re, glob, sys, time
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import binomtest, mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

t0 = time.time()

BASE   = "./"
DATA_DIR = f"{BASE}/data"
RESULTS  = f"{BASE}/results"
FIGURES  = f"{BASE}/figures"
MEME   = f"{DATA_DIR}/epigenome/dap_seq/gem02_rep_memechip00_onefile/gem02_rep_memechip00_m1.txt"
B73FA  = f"{DATA_DIR}/Zm-B73-REFERENCE-NAM-5.0.fa.gz"
PEAKS_DIR = f"{DATA_DIR}/epigenome/dap_seq/normalized_specific_and_shared_peaks"
B73_PEAK_DIR = f"{PEAKS_DIR}/B73v5_MATCH_Mo17-B73v5_specific_B73v5_shared"
MO17_PEAK_DIR = f"{PEAKS_DIR}/Mo17_MATCH_B73v5-Mo17_specific_Mo17_shared"
# VCF_SMALL  = f"{DATA_DIR}/Mo17_toB73v5_paf_syri_noStartPOS0_INDELs_less50bp.vcf"
VCF_SMALL  = f"{DATA_DIR}/Mo17_toB73v5_paf_syri_noStartPOS0_SNPs.vcf"
VCF_LARGE  = f"{DATA_DIR}/Mo17_toB73v5_paf_syri_noStartPOS0_INDELs_more50bp.vcf"
DAP_TSV    = f"{RESULTS}/dap_seq_per_gene.tsv"
GXE_FILE   = f"{DATA_DIR}/GxE_gene_IDs.txt"
BG_FILE    = f"{DATA_DIR}/background_gene_IDs.txt"
GXE_RES    = f"{DATA_DIR}/GxE_allele_specific_test_results.txt"
OUT_DIR    = RESULTS
FLANK      = 50   # bp flanking the indel on each side for PWM scanning

PSEUDO = 0.01    # pseudocount before log-odds to avoid log(0)
BG_FREQ = 0.25   # uniform background

print("=" * 60)
print("Mechanistic Indel × PWM Scoring")
print("=" * 60)

# ── 1. Parse MEME file ────────────────────────────────────────────────────────
print("\n[1] Parsing MEME motif file …")

def rc_dna(seq):
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]

def score_sequence(seq, log_odds, width):
    """Return max log-odds score over all positions and both strands."""
    seq = seq.upper()
    if len(seq) < width:
        return -np.inf
    nmap = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    best = -np.inf
    for strand_seq in [seq, rc_dna(seq)]:
        for i in range(len(strand_seq) - width + 1):
            kmer = strand_seq[i:i+width]
            s = sum(log_odds[j, nmap.get(nuc, -1)]
                    for j, nuc in enumerate(kmer)
                    if nmap.get(nuc, -1) >= 0)
            if s > best:
                best = s
    return best

tf_pwm = {}   # {tf_name: {'B73': {'log_odds': array, 'width': int},
              #             'Mo17': {'log_odds': array, 'width': int}}}

current_name = None
current_geno = None
current_tf   = None
current_width = None
reading_matrix = False
matrix_rows    = []

with open(MEME) as fh:
    for line in fh:
        line = line.rstrip()
        if line.startswith("MOTIF"):
            # Save previous if any
            if current_tf and matrix_rows:
                mat = np.array(matrix_rows, dtype=float)
                # Add pseudocount and convert to log-odds
                mat = (mat + PSEUDO) / (1 + 4 * PSEUDO)
                log_odds = np.log2(mat / BG_FREQ)
                if current_tf not in tf_pwm:
                    tf_pwm[current_tf] = {}
                tf_pwm[current_tf][current_geno] = {
                    'log_odds': log_odds,
                    'width': current_width if current_width else len(mat)
                }
            # Parse new MOTIF line
            pts = line.split()
            current_name = pts[1]  # e.g. ABI19_B73_1_m1 or ABI19_Mo17_1
            # Extract TF name and genotype from the motif name
            m = re.search(r'^(.+?)_(B73|Mo17)_', current_name)
            if m:
                current_tf   = m.group(1)
                current_geno = m.group(2)
            else:
                current_tf   = None
                current_geno = None
            current_width  = None
            matrix_rows    = []
            reading_matrix = False

        elif line.startswith("letter-probability matrix"):
            reading_matrix = True
            m = re.search(r'w=\s*(\d+)', line)
            if m:
                current_width = int(m.group(1))

        elif reading_matrix:
            parts = line.split()
            if len(parts) == 4:
                try:
                    matrix_rows.append([float(x) for x in parts])
                except ValueError:
                    reading_matrix = False
            else:
                reading_matrix = False

# Save last motif
if current_tf and matrix_rows:
    mat = np.array(matrix_rows, dtype=float)
    mat = (mat + PSEUDO) / (1 + 4 * PSEUDO)
    log_odds = np.log2(mat / BG_FREQ)
    if current_tf not in tf_pwm:
        tf_pwm[current_tf] = {}
    tf_pwm[current_tf][current_geno] = {
        'log_odds': log_odds,
        'width': current_width if current_width else len(mat)
    }

n_b73_pwm  = sum(1 for v in tf_pwm.values() if 'B73'  in v)
n_mo17_pwm = sum(1 for v in tf_pwm.values() if 'Mo17' in v)
n_both     = sum(1 for v in tf_pwm.values() if 'B73' in v and 'Mo17' in v)
print(f"  TFs with B73 PWM:  {n_b73_pwm}")
print(f"  TFs with Mo17 PWM: {n_mo17_pwm}")
print(f"  TFs with both:     {n_both}")

# ── 2. Load genotype-specific DAP peaks, record TF per peak ──────────────────
print("\n[2] Loading genotype-specific DAP peaks …")

def extract_tf_from_peak_filename(fname, geno):
    """Extract TF name from peak file basename.
    B73 example: ABI19_B73_1_MATCH_ABI19_Mo17_1-B73v5_specific.narrowPeak
    Mo17 example: ABI19_Mo17_1_MATCH_ABI19_B73_1-Mo17_specific.narrowPeak
    """
    bn = os.path.basename(fname)
    # Get first field before _MATCH_
    first_part = bn.split("_MATCH_")[0]
    m = re.search(r'^(.+?)_(B73|Mo17)_', first_part)
    return m.group(1) if m else None

def load_peaks_with_tf(peak_dir, geno, pattern="*specific.narrowPeak"):
    """Load all specific peak files; return list of (chr, start, end, tf)."""
    peaks = []
    files = glob.glob(os.path.join(peak_dir, pattern))
    for fpath in files:
        tf = extract_tf_from_peak_filename(fpath, geno)
        if tf is None:
            continue
        with open(fpath) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                try:
                    chrom, s, e = parts[0], int(parts[1]), int(parts[2])
                    peaks.append((chrom, s, e, tf))
                except (ValueError, IndexError):
                    continue
    return peaks

b73_peaks  = load_peaks_with_tf(B73_PEAK_DIR, "B73",  "*specific.narrowPeak")
mo17_peaks = load_peaks_with_tf(MO17_PEAK_DIR, "Mo17", "*specific.narrowPeak")
print(f"  B73-specific peaks:  {len(b73_peaks):,}")
print(f"  Mo17-specific peaks: {len(mo17_peaks):,}")

# Build per-chromosome interval trees (sorted arrays for searchsorted)
def build_chr_index(peaks):
    """peaks: list of (chr, start, end, tf)
    Returns: {chrom: {'starts': array, 'ends': array, 'tfs': list}}
    """
    from collections import defaultdict
    by_chr = defaultdict(list)
    for chrom, s, e, tf in peaks:
        by_chr[chrom].append((s, e, tf))
    idx = {}
    for chrom, items in by_chr.items():
        items.sort(key=lambda x: x[0])
        idx[chrom] = {
            'starts': np.array([x[0] for x in items], dtype=np.int64),
            'ends':   np.array([x[1] for x in items], dtype=np.int64),
            'tfs':    [x[2] for x in items]
        }
    return idx

b73_idx  = build_chr_index(b73_peaks)
mo17_idx = build_chr_index(mo17_peaks)

def find_overlapping_peaks(chrom, pos, idx):
    """Return list of (tf,) for peaks overlapping pos."""
    if chrom not in idx:
        return []
    d = idx[chrom]
    # Peaks where start <= pos < end
    lo = np.searchsorted(d['starts'], pos, side='right') - 1
    # Check from lo downward (starts sorted, but end can vary)
    results = []
    # Use a safe range search
    hi = np.searchsorted(d['starts'], pos + 1, side='left')
    for i in range(max(0, lo), min(len(d['starts']), hi + 1)):
        if d['starts'][i] <= pos < d['ends'][i]:
            results.append(d['tfs'][i])
    return results

# ── 3. Load VCF indels ────────────────────────────────────────────────────────
print("\n[3] Loading VCF indels …")

def load_vcf(vcf_path):
    """Returns DataFrame with chr, pos(0-based), ref, alt columns."""
    rows = []
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 5)
            if len(parts) < 5:
                continue
            chrom = parts[0].strip()
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            try:
                pos1 = int(parts[1])   # 1-based VCF
                ref  = parts[3].strip().upper()
                alt  = parts[4].strip().upper()
                rows.append((chrom, pos1 - 1, ref, alt))   # convert to 0-based
            except (ValueError, IndexError):
                continue
    return pd.DataFrame(rows, columns=["chr", "pos", "ref", "alt"])

vcf_small = load_vcf(VCF_SMALL)
vcf_large = load_vcf(VCF_LARGE)
# vcf_all   = pd.concat([vcf_small, vcf_large], ignore_index=True)
vcf_all   = load_vcf(VCF_SMALL)
print(f"  Small indels: {len(vcf_small):,}")
print(f"  Large indels: {len(vcf_large):,}")
print(f"  Total:        {len(vcf_all):,}")

# Build per-chromosome sorted arrays for fast position lookup
vcf_by_chr = {}
for chrom, grp in vcf_all.groupby("chr"):
    grp_sorted = grp.sort_values("pos").reset_index(drop=True)
    vcf_by_chr[chrom] = grp_sorted

# ── 4. Load gene info ─────────────────────────────────────────────────────────
print("\n[4] Loading gene info …")

dap = pd.read_csv(DAP_TSV, sep="\t")
gxe_ids = set(open(GXE_FILE).read().split())
bg_ids  = set(open(BG_FILE).read().split()) - gxe_ids

# Load GxE test results for lFC
res = pd.read_csv(GXE_RES, sep="\t").dropna(subset=["log2FoldChange"])
res_map = dict(zip(res["GeneID"], res["log2FoldChange"]))

print(f"  Genes in dap_seq_per_gene.tsv: {len(dap):,}")
print(f"  GxE-ASE genes: {dap['is_gxe'].sum():,}  BG genes: {dap['is_bg'].sum():,}")

# ── 5. Open genome FASTA ──────────────────────────────────────────────────────
print("\n[5] Opening genome FASTA …")
fa = pysam.FastaFile(B73FA)
print(f"  References: {fa.nreferences}")

# ── 6. Score indels at genotype-specific peaks ────────────────────────────────
print("\n[6] Scoring indels at genotype-specific peaks …")
print("    (this may take a few minutes)")

results = []
n_no_tf_pwm = 0
n_scored    = 0
n_skipped_seq = 0

# Only process genes with at least one genotype-specific indel in the promoter
genes_to_process = dap[(dap["is_gxe"] | dap["is_bg"]) &
                       (dap["n_b73_spec"] + dap["n_mo17_spec"] > 0)]
print(f"  Genes with ≥1 genotype-specific peak: {len(genes_to_process):,}")

for _, gene_row in genes_to_process.iterrows():
    gene_id = gene_row["gene"]
    chrom   = gene_row["chr"]
    prom_s  = int(gene_row["prom_s"])
    prom_e  = int(gene_row["prom_e"])
    is_gxe  = bool(gene_row["is_gxe"])
    is_bg   = bool(gene_row["is_bg"])
    lfc     = res_map.get(gene_id, np.nan)

    # Find indels in this promoter
    if chrom not in vcf_by_chr:
        continue
    chrom_vcf = vcf_by_chr[chrom]
    lo = np.searchsorted(chrom_vcf["pos"].values, prom_s, side="left")
    hi = np.searchsorted(chrom_vcf["pos"].values, prom_e, side="right")
    if lo >= hi:
        continue

    promo_indels = chrom_vcf.iloc[lo:hi]

    for _, indel in promo_indels.iterrows():
        pos0 = int(indel["pos"])    # 0-based
        ref  = indel["ref"]
        alt  = indel["alt"]

        # Skip non-indels (same length = SNP or MNP)
        # if len(ref) == len(alt):
        # continue

        # Check which genotype-specific peaks this indel overlaps
        for peak_geno, peak_idx in [("B73", b73_idx), ("Mo17", mo17_idx)]:
            overlapping_tfs = find_overlapping_peaks(chrom, pos0, peak_idx)
            if not overlapping_tfs:
                continue

            for tf in set(overlapping_tfs):
                # Get PWMs for this TF
                if tf not in tf_pwm:
                    n_no_tf_pwm += 1
                    continue
                tf_data = tf_pwm[tf]

                # For B73-specific peaks: score REF with B73 PWM, ALT with Mo17 PWM
                # For Mo17-specific peaks: score REF with B73 PWM, ALT with Mo17 PWM
                # (same operation; sign of Δ tells the story)
                if "B73" not in tf_data or "Mo17" not in tf_data:
                    n_no_tf_pwm += 1
                    continue

                b73_lo  = tf_data["B73"]["log_odds"]
                b73_w   = tf_data["B73"]["width"]
                mo17_lo = tf_data["Mo17"]["log_odds"]
                mo17_w  = tf_data["Mo17"]["width"]

                # Fetch genome context
                ctx_start = max(0, pos0 - FLANK)
                ctx_end   = pos0 + len(ref) + FLANK
                try:
                    genome_ctx = fa.fetch(chrom, ctx_start, ctx_end).upper()
                except Exception:
                    n_skipped_seq += 1
                    continue

                offset = pos0 - ctx_start

                # Build REF context (B73 genome — already has the ref allele)
                ref_ctx = genome_ctx

                # Build ALT context (Mo17 — replace ref with alt)
                alt_ctx = genome_ctx[:offset] + alt + genome_ctx[offset + len(ref):]

                # Verify REF allele matches genome
                genome_ref = genome_ctx[offset: offset + len(ref)]
                if genome_ref.upper() != ref.upper():
                    # REF mismatch — skip (stale coordinate or chr naming issue)
                    n_skipped_seq += 1
                    continue

                # Score both sequences with both PWMs
                score_ref_b73  = score_sequence(ref_ctx,  b73_lo,  b73_w)
                score_alt_mo17 = score_sequence(alt_ctx, mo17_lo, mo17_w)
                delta = score_ref_b73 - score_alt_mo17

                # Also record the reciprocal scores for diagnostics
                score_ref_mo17 = score_sequence(ref_ctx,  mo17_lo, mo17_w)
                score_alt_b73  = score_sequence(alt_ctx,  b73_lo,  b73_w)

                n_scored += 1
                results.append(dict(
                    gene=gene_id,
                    is_gxe=is_gxe,
                    is_bg=is_bg,
                    lfc=lfc,
                    chrom=chrom,
                    pos=pos0,
                    ref=ref,
                    alt=alt,
                    tf=tf,
                    peak_geno=peak_geno,  # "B73" or "Mo17"
                    score_ref_b73=score_ref_b73,
                    score_alt_mo17=score_alt_mo17,
                    score_ref_mo17=score_ref_mo17,
                    score_alt_b73=score_alt_b73,
                    delta=delta,         # score_ref_b73 - score_alt_mo17
                    abs_delta=abs(delta)
                ))

elapsed = time.time() - t0
print(f"\n  Scored:               {n_scored:,} indel × peak × TF records")
print(f"  Skipped (no PWM):     {n_no_tf_pwm:,}")
print(f"  Skipped (seq issue):  {n_skipped_seq:,}")
print(f"  Elapsed: {elapsed:.1f}s")

if len(results) == 0:
    print("\n  ERROR: No records scored. Check MEME parsing and peak loading.")
    sys.exit(1)

df = pd.DataFrame(results)
print(f"\n  Records per category:")
print(f"    GxE  × B73-specific peak:  {(df['is_gxe'] & (df['peak_geno']=='B73')).sum():,}")
print(f"    GxE  × Mo17-specific peak: {(df['is_gxe'] & (df['peak_geno']=='Mo17')).sum():,}")
print(f"    BG   × B73-specific peak:  {(df['is_bg']  & (df['peak_geno']=='B73')).sum():,}")
print(f"    BG   × Mo17-specific peak: {(df['is_bg']  & (df['peak_geno']=='Mo17')).sum():,}")

# ── 7. Statistical tests ──────────────────────────────────────────────────────
print("\n[7] Statistical Tests")
print("─" * 50)

stat_rows = []

for pg in ["B73", "Mo17"]:
    sub = df[df["peak_geno"] == pg]
    gxe_sub = sub[sub["is_gxe"]]
    bg_sub  = sub[sub["is_bg"]]

    if len(sub) == 0:
        print(f"  {pg}-specific peaks: NO DATA")
        continue

    print(f"\n── {pg}-specific peaks (n={len(sub):,}) ──")

    # Test 1: directional sign test
    # B73-specific: expect Δ > 0; Mo17-specific: expect Δ < 0
    expected_positive = (pg == "B73")
    n_correct = (sub["delta"] > 0).sum() if expected_positive else (sub["delta"] < 0).sum()
    n_total   = len(sub)
    binom_p   = binomtest(n_correct, n_total, p=0.5, alternative="greater").pvalue
    pct        = 100 * n_correct / n_total
    direction  = "Δ>0" if expected_positive else "Δ<0"
    print(f"  Sign test ({direction}): {n_correct}/{n_total} ({pct:.1f}%)  binomial p={binom_p:.4g}")

    # Separately for GxE vs BG
    for label, grp in [("GxE", gxe_sub), ("BG", bg_sub)]:
        if len(grp) == 0:
            continue
        n_c = (grp["delta"] > 0).sum() if expected_positive else (grp["delta"] < 0).sum()
        pct_g = 100 * n_c / len(grp)
        bp    = binomtest(n_c, len(grp), p=0.5, alternative="greater").pvalue
        print(f"    {label}  ({len(grp):,} records): {n_c}/{len(grp)} ({pct_g:.1f}%)  binomial p={bp:.4g}")

    # Test 2: |Δ| GxE vs BG (Mann-Whitney)
    if len(gxe_sub) > 0 and len(bg_sub) > 0:
        mwu_stat, mwu_p = mannwhitneyu(
            gxe_sub["abs_delta"], bg_sub["abs_delta"], alternative="greater")
        print(f"  |Δscore| GxE vs BG: MWU p={mwu_p:.4g}  "
              f"(GxE mean={gxe_sub['abs_delta'].mean():.3f},  BG mean={bg_sub['abs_delta'].mean():.3f})")

    # Test 3: Spearman |lFC| vs |Δ| (only GxE genes)
    gxe_both = gxe_sub.dropna(subset=["lfc"])
    # Collapse to gene-level (mean |Δ| per gene)
    gene_level = gxe_both.groupby("gene").agg(
        abs_lfc=("lfc", lambda x: x.abs().mean()),
        mean_abs_delta=("abs_delta", "mean")
    ).reset_index()
    if len(gene_level) >= 10:
        rho, sp_p = spearmanr(gene_level["abs_lfc"], gene_level["mean_abs_delta"])
        print(f"  Spearman |lFC| ~ |Δscore| (gene-level, GxE, n={len(gene_level)}): "
              f"ρ={rho:.3f}  p={sp_p:.4g}")

    stat_rows.append(dict(
        peak_geno=pg,
        n_records=len(sub),
        n_gxe=len(gxe_sub),
        n_bg=len(bg_sub),
        pct_correct_sign=pct,
        binom_p=binom_p,
        mwu_p=mwu_p if (len(gxe_sub)>0 and len(bg_sub)>0) else np.nan,
        mean_abs_delta_gxe=gxe_sub["abs_delta"].mean() if len(gxe_sub)>0 else np.nan,
        mean_abs_delta_bg=bg_sub["abs_delta"].mean() if len(bg_sub)>0 else np.nan,
    ))

stat_df = pd.DataFrame(stat_rows)

# Per-TF summary (collapse to gene × TF level for cleaner stats)
print("\n── Top TFs by Δscore magnitude (≥20 GxE records) ──")
tf_summary = (
    df[df["is_gxe"]]
      .groupby(["tf", "peak_geno"])
      .agg(n=("delta", "size"),
           mean_delta=("delta", "mean"),
           mean_abs_delta=("abs_delta", "mean"),
           pct_correct=("delta", lambda x:
               100 * (x > 0).sum() / len(x) if x.name else 0))
      .reset_index()
)
# Recalculate pct_correct properly
rows_top = []
for (tf, pg), grp in df[df["is_gxe"]].groupby(["tf", "peak_geno"]):
    if len(grp) < 20:
        continue
    expected_pos = (pg == "B73")
    n_c = (grp["delta"] > 0).sum() if expected_pos else (grp["delta"] < 0).sum()
    rows_top.append(dict(
        tf=tf, peak_geno=pg, n=len(grp),
        mean_delta=grp["delta"].mean(),
        mean_abs_delta=grp["abs_delta"].mean(),
        pct_correct_sign=100*n_c/len(grp)
    ))
if rows_top:
    top_df = pd.DataFrame(rows_top).sort_values("mean_abs_delta", ascending=False)
    print(top_df.head(20).to_string(index=False))

# ── 8. Figure ─────────────────────────────────────────────────────────────────
print("\n[8] Generating figures …")

fig, axes = plt.subplots(2, 2, figsize=(12, 9))
fig.subplots_adjust(hspace=0.45, wspace=0.35)
colors = {"GxE-ASE": "#2196F3", "Background": "#9E9E9E"}

for ax_col, pg in enumerate(["B73", "Mo17"]):
    sub     = df[df["peak_geno"] == pg]
    gxe_sub = sub[sub["is_gxe"]]
    bg_sub  = sub[sub["is_bg"]]

    # Top: Δscore distribution
    ax = axes[0, ax_col]
    for grp, label, color in [(gxe_sub, "GxE-ASE", colors["GxE-ASE"]),
                               (bg_sub,  "Background", colors["Background"])]:
        if len(grp) == 0:
            continue
        vals = grp["delta"].clip(-10, 10)
        ax.hist(vals, bins=40, density=True, alpha=0.5,
                color=color, label=f"{label} (n={len(grp):,})", histtype="stepfilled")
    ax.axvline(0, color="black", lw=1, ls="--")
    expected_sign = ">" if pg == "B73" else "<"
    ax.set_xlabel(f"Δscore (B73_REF PWM − Mo17_ALT PWM)")
    ax.set_ylabel("Density")
    ax.set_title(f"{pg}-specific peaks\n"
                 f"Expected Δ{expected_sign}0 if {pg} has better motif",
                 fontsize=9)
    ax.legend(fontsize=7)

    # Bottom: |Δscore| GxE vs BG violin-style (histogram overlay)
    ax2 = axes[1, ax_col]
    for grp, label, color in [(gxe_sub, "GxE-ASE", colors["GxE-ASE"]),
                               (bg_sub,  "Background", colors["Background"])]:
        if len(grp) == 0:
            continue
        vals = grp["abs_delta"].clip(0, 15)
        ax2.hist(vals, bins=30, density=True, alpha=0.5,
                 color=color, label=f"{label} (mean={grp['abs_delta'].mean():.2f})",
                 histtype="stepfilled")
    ax2.set_xlabel("|Δscore|")
    ax2.set_ylabel("Density")
    ax2.set_title(f"|Δscore| distribution — {pg}-specific peaks", fontsize=9)
    ax2.legend(fontsize=7)

fig.suptitle("Mechanistic Indel × PWM Scoring\n"
             "Δscore = score(B73_REF, B73_PWM) − score(Mo17_ALT, Mo17_PWM)",
             fontsize=11, fontweight="bold")

out_fig = f"{FIGURES}/indel_pwm_scoring_figure.png"
fig.savefig(out_fig, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved → {out_fig}")

# ── 9. Save outputs ───────────────────────────────────────────────────────────
df.to_csv(f"{OUT_DIR}/indel_pwm_scored_records.tsv", sep="\t", index=False)
stat_df.to_csv(f"{OUT_DIR}/indel_pwm_stats.tsv", sep="\t", index=False)
print(f"  Saved → {OUT_DIR}/indel_pwm_scored_records.tsv")
print(f"  Saved → {OUT_DIR}/indel_pwm_stats.tsv")

print(f"\n✓ Done in {time.time()-t0:.1f}s")
