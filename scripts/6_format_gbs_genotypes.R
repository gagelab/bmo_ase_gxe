################################################################################
# 6_format_gbs_genotypes.R
#
# Build IBM recombination bins from HapMap v3 GBS data.
#
# Pipeline:
#   1. Load IBM sample subset of HapMap v3 GBS VCF (AGPv4 coordinates)
#   2. Filter to Mo17-segregating SNPs; resolve duplicate samples
#   3. Lift over coordinates from AGPv4 -> NAM v5
#   4. Impute missing genotypes with qtl2 HMM (riself model)
#   5. Remove individuals with excessive crossover counts (top 10%)
#   6. Collapse consecutive identical calls into recombination bins
#   7. Disjoin bins across all individuals -> consensus bin x sample matrix
#
# Reads (external, not in repo — see RERUN_WORKFLOW.md for download sources):
#   data/ZeaGBSv27_IBM_raw_AGPv4.vcf.gz   (IBM sample subset of HapMap v3 GBS)
#   data/AllZeaGBSv2.7_publicSamples_metadata20140411.xlsx
#
# Writes:
#   data/IBM_recomb_bins_fromGBS.tsv
#
# Runtime: ~30–60 min (HMM imputation dominates)
################################################################################

library(tidyverse)
library(rtracklayer)
library(qtl2)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE          <- "./"
DATA      <- file.path(BASE, "data")   # large shared data files
# Raw GBS VCF and metadata must be downloaded — see RERUN_WORKFLOW.md
GBS_VCF       <- file.path(DATA, "ZeaGBSv27_IBM_raw_AGPv4.vcf.gz")
META_XLSX     <- file.path(DATA, "AllZeaGBSv2.7_publicSamples_metadata20140411.xlsx")
CHAIN_URL     <- "https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/chain_files/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain"
OUT_BINS      <- file.path(DATA, "IBM_recomb_bins_fromGBS.tsv")

Sys.setenv("VROOM_CONNECTION_SIZE" = 10000000)

# ── 1. Load and encode genotypes ───────────────────────────────────────────────
message("1. Loading GBS VCF ...")
geno <- read_tsv(GBS_VCF, comment = "##")
genomat <- as.matrix(geno[, 10:ncol(geno)])
genomat[genomat == "./."] <- NA
genomat[genomat == "0/0"] <- 0
genomat[genomat == "1/1"] <- 2
genomat[genomat %in% c("0/1", "1/0")] <- NA   # hets -> missing
genomat[!is.na(genomat) & !(genomat %in% c(0, 2))] <- NA
mode(genomat) <- "numeric"

# ── 2. Filter to Mo17-informative SNPs ────────────────────────────────────────
message("2. Filtering to Mo17-segregating SNPs ...")
mo17_cols    <- grepl("Mo17", colnames(genomat))
mo17_summary <- tibble(
  SNP  = seq_len(nrow(genomat)),
  mean = rowMeans(genomat[, mo17_cols], na.rm = TRUE),
  n    = rowSums(!is.na(genomat[, mo17_cols]))
)
# Keep SNPs where all typed Mo17 samples are homozygous Mo17 (mean == 2, n >= 2)
seg_snps <- mo17_summary %>% filter(mean == 2, n >= 2) %>% pull(SNP)
genomat   <- genomat[seg_snps, ]
message(sprintf("   %d Mo17-segregating SNPs retained", nrow(genomat)))

# ── 3. Merge duplicate IBM samples ────────────────────────────────────────────
message("3. Resolving duplicate samples ...")
meta <- readxl::read_xlsx(META_XLSX) %>%
  filter(GermplasmSet == "IBM" | DNASample == "Mo17")
ids <- grep("M[0-9]{3}", unique(meta$DNASample), value = TRUE)

# Drop samples where duplicates show <95% concordance (likely mix-ups)
agree_rates <- vapply(ids, function(id) {
  s <- genomat[, grepl(id, colnames(genomat)), drop = FALSE]
  if (ncol(s) < 2) return(NA_real_)
  sum(s[, 1] == s[, 2], na.rm = TRUE) /
    sum(!is.na(s[, 1]) & !is.na(s[, 2]))
}, numeric(1))

drop_ids <- names(agree_rates[!is.na(agree_rates) & agree_rates < 0.95])
drop_ids  <- c(drop_ids, "M0132(mix)")

keep_ids <- ids[!ids %in% drop_ids]
merged   <- matrix(NA_real_, nrow = nrow(genomat), ncol = length(keep_ids),
                   dimnames = list(NULL, keep_ids))

for (id in keep_ids) {
  s <- genomat[, grepl(id, colnames(genomat)), drop = FALSE]
  merged[, id] <- if (ncol(s) > 1) rowMeans(s, na.rm = TRUE) else s[, 1]
}

# ── 4. Quality filter: missingness and allele frequency ───────────────────────
message("4. Applying SNP and sample quality filters ...")
miss_ind <- colMeans(is.na(merged))
miss_snp <- rowMeans(is.na(merged))
af       <- rowMeans(merged, na.rm = TRUE)

keep_ind <- miss_ind < 0.8
keep_snp <- miss_snp < 0.8 & af > 0.5 & af < 1.5

merged      <- merged[keep_snp, keep_ind]
geno_meta_v4 <- geno[seg_snps, ][keep_snp, 1:9]
message(sprintf("   %d SNPs x %d samples after filtering",
                nrow(merged), ncol(merged)))

# ── 5. Lift over v4 -> v5 coordinates ────────────────────────────────────────
message("5. Lifting over AGPv4 -> NAM v5 coordinates ...")
chain_file <- file.path(DATA, basename(CHAIN_URL))
if (!file.exists(chain_file)) {
  download.file(CHAIN_URL, chain_file)
}
chain    <- import.chain(chain_file)
coords_v4 <- makeGRangesFromDataFrame(
  geno_meta_v4,
  seqnames.field = "#CHROM",
  start.field    = "POS",
  end.field      = "POS"
)
coords_v5  <- liftOver(coords_v4, chain)
keep_sites <- which(lengths(coords_v5) == 1)
coords_v5  <- unlist(coords_v5[keep_sites])

merged       <- merged[keep_sites, ]
geno_meta_v5 <- geno_meta_v4[keep_sites, ]
geno_meta_v5$`#CHROM` <- as.numeric(seqnames(coords_v5))
geno_meta_v5$POS      <- as.numeric(start(coords_v5))
geno_meta_v5          <- arrange(geno_meta_v5, `#CHROM`, POS)
merged                <- merged[order(as.numeric(seqnames(coords_v5)),
                                      as.numeric(start(coords_v5))), ]
message(sprintf("   %d SNPs retained after liftover", nrow(merged)))

# ── 6. Impute with qtl2 HMM ───────────────────────────────────────────────────
message("6. Imputing genotypes with qtl2 HMM ...")

# Encode as character for qtl2 (A = B73, B = Mo17, - = missing)
geno_char <- matrix("-", nrow = nrow(merged), ncol = ncol(merged),
                    dimnames = dimnames(merged))
geno_char[!is.na(merged) & merged == 0] <- "A"
geno_char[!is.na(merged) & merged == 2] <- "B"

marker_ids <- paste0("S", geno_meta_v5$`#CHROM`, "_", geno_meta_v5$POS)
rownames(geno_char) <- marker_ids
geno_for_qtl2 <- t(geno_char)   # samples x markers

qtl2_dir <- file.path(DATA, "qtl2_input")
dir.create(qtl2_dir, showWarnings = FALSE)

write.csv(
  data.frame(id = rownames(geno_for_qtl2), geno_for_qtl2, check.names = FALSE),
  file.path(qtl2_dir, "geno.csv"), row.names = FALSE, quote = FALSE
)

# Use physical position (Mb * 4 expansion) as a proxy genetic map for IBM RILs
gmap <- data.frame(
  marker = marker_ids,
  chr    = geno_meta_v5$`#CHROM`,
  pos    = geno_meta_v5$POS / 1e6 * 4
)
write.csv(gmap, file.path(qtl2_dir, "gmap.csv"), row.names = FALSE, quote = FALSE)
write.csv(gmap, file.path(qtl2_dir, "pmap.csv"), row.names = FALSE, quote = FALSE)

qtl2::write_control_file(
  file.path(qtl2_dir, "control.yaml"),
  crosstype  = "riself",
  geno_file  = "geno.csv",
  gmap_file  = "gmap.csv",
  pmap_file  = "pmap.csv",
  geno_codes = c(A = 1, B = 2),
  alleles    = c("A", "B"),
  overwrite  = TRUE
)

cross   <- read_cross2(file.path(qtl2_dir, "control.yaml"))
pr      <- calc_genoprob(cross, error_prob = 0.05, map_function = "haldane")
imputed <- maxmarg(pr, minprob = 0.95)   # 1 = B73, 2 = Mo17; NA where uncertain

# ── 7. Remove individuals with excessive recombination ────────────────────────
message("7. Filtering high-recombination individuals ...")
xo_counts <- count_xo(imputed)
total_xo  <- rowSums(xo_counts)
cutoff    <- quantile(total_xo, 0.9)
excessive <- names(total_xo)[total_xo > cutoff]
message(sprintf("   Removing %d individuals (>%.0f total crossovers)",
                length(excessive), cutoff))
for (chr in seq_along(imputed)) {
  imputed[[chr]] <- imputed[[chr]][!rownames(imputed[[chr]]) %in% excessive, ]
}

# ── 8. Collapse runs -> recombination bins ────────────────────────────────────
message("8. Calling recombination bins ...")
call_long <- imap_dfr(imputed, function(mat, chr) {
  as_tibble(mat, rownames = "Sample") %>%
    pivot_longer(-Sample, names_to = "marker", values_to = "call") %>%
    filter(!is.na(call)) %>%
    mutate(
      chr = as.numeric(sub("S([0-9]{1,2})_[0-9]+",  "\\1", marker)),
      pos = as.numeric(sub("S[0-9]{1,2}_([0-9]+)", "\\1", marker))
    )
})

bins_by_sample <- call_long %>%
  arrange(Sample, chr, pos) %>%
  group_by(Sample, chr) %>%
  mutate(bin_id = cumsum(c(1L, diff(call) != 0L))) %>%
  group_by(Sample, chr, bin_id, call) %>%
  summarise(start = min(pos), end = max(pos), n_markers = n(), .groups = "drop") %>%
  split(.$Sample) %>%
  lapply(function(df) {
    GRanges(
      seqnames  = df$chr,
      ranges    = IRanges(df$start, df$end),
      genotype  = c("B73", "Mo17")[df$call],
      n_markers = df$n_markers
    )
  }) %>%
  GRangesList()

# ── 9. Disjoin to consensus bins across all individuals ───────────────────────
message("9. Constructing consensus bin x sample matrix ...")
all_ranges    <- unlist(bins_by_sample)
consensus_bins <- disjoin(all_ranges)

geno_matrix <- vapply(names(bins_by_sample), function(samp) {
  gr   <- bins_by_sample[[samp]]
  hits <- findOverlaps(consensus_bins, gr)
  calls <- rep(NA_character_, length(consensus_bins))
  calls[queryHits(hits)] <- gr$genotype[subjectHits(hits)]
  calls
}, character(length(consensus_bins)))

consensus_df <- as.data.frame(consensus_bins) %>%
  dplyr::select(seqnames, start, end) %>%
  dplyr::rename(chr = seqnames) %>%
  bind_cols(as_tibble(geno_matrix))

write_tsv(consensus_df, OUT_BINS)
message(sprintf("Done. %d consensus bins written to %s", nrow(consensus_df), OUT_BINS))
