################################################################################
# 0_deseq_padj01.R
#
# GxE allele-specific expression analysis using DESeq2.
# This is a corrected copy of 0_deseq.R with the following changes:
#
#   1. FIXED BUG — toPlot/bg_genes forward reference:
#      In the original script, `toPlot` (used to define sig_ge and bg_genes at
#      lines 62-94) was not constructed until lines 105-112. The code only ran
#      correctly in interactive mode. This version defines toPlot and bg_genes
#      immediately after DESeq2 results are extracted.
#
#   2. ADDED — padj column in output:
#      The original toPlot join omitted padj from the `ge` results object.
#      DEG_GxE_results.txt now includes padj.
#
#   3. CHANGED — GxE threshold: p < 0.1 -> padj < 0.1:
#      FDR < 0.1 (Benjamini-Hochberg) is more principled than raw p < 0.1.
#      Threshold sensitivity analysis confirms all enrichment signals
#      strengthen at this threshold (DAP OR 1.89, SV OR 3.14, indel×DAP OR
#      1.88) with a clean background of comparable size (~11,389 genes).
#
#   4. ADDED — explicit pre-filter for background gene set:
#      Background genes must have >10 counts in >5 sample×allele columns.
#      This was already the implicit filter in the original analysis (the
#      intersection of tested genes and this count filter happened to produce
#      exactly 11,637 genes), but is now explicitly coded and documented.
#      The pre-filter ensures all background genes had sufficient expression
#      to be realistically testable for allele-specific effects, providing a
#      biologically meaningful comparison universe.
#
################################################################################

# setwd("/sessions/dazzling-dreamy-archimedes/mnt/bmo_ase/gxe_ase")

library(tidyverse)
library(magrittr)
library(DESeq2)

padj_thresh <- 0.1   # FDR threshold for GxE gene set definition

# Directory layout (relative to project root gxe_ase/)
DATA_DIR    <- "data"
RESULTS_DIR <- "results"
FIGURES_DIR <- "figures"

# ------------------------------------------------------------------------------
# DATA LOADING
# ------------------------------------------------------------------------------

meta <- read_tsv(file.path(DATA_DIR, "meta_for_DESeq2.txt")) %>%
  rename("Condition" = "Environment")
counts_raw <- read_tsv(file.path(DATA_DIR, "counts_for_DESeq2.txt"))
geneids <- counts_raw$Geneid[, drop = TRUE]
counts <- as.matrix(counts_raw[, -1])
rownames(counts) <- geneids

# Drop two NC observations to keep the dataset balanced.
# A10 and A11 were dropped in the original analysis; models fail with them
# included due to imbalance in the Environment:sampleID design term.
drops <- c("A10", "A11")
meta2   <- meta %>% filter(!LibraryID %in% drops)
counts2 <- counts[, !grepl(paste0(drops, collapse = "|"), colnames(counts))]

# Add within-environment sample ID (runs 1:n_samples within each env)
sampleID <- meta2 %>%
  distinct(Environment, LibraryID) %>%
  group_by(Environment) %>%
  mutate(sampleID = 1:n())

meta2 <- left_join(meta2, sampleID) %>%
  mutate(
    sampleID    = factor(sampleID),
    Allele      = factor(Allele, levels = c("Ref", "Alt")),
    Environment = factor(Environment, levels = c("NC", "MO"))
  )

# Long-format counts for downstream pre-filter
counts2_long <- counts2 %>%
  as_tibble(rownames = "GeneID") %>%
  pivot_longer(-GeneID, names_to = "RowID", values_to = "Count")

# ------------------------------------------------------------------------------
# DESEQ2 MODEL
#
# Design follows the allele-specific expression approach described at:
# https://rpubs.com/mikelove/ase
#
# The model tests for changes in the Alt/Ref (Mo17/B73) ratio between
# environments (MO vs NC). Key design terms:
#   ~ Environment + Environment:sampleID + Environment:Allele
#
# Environment:sampleID accounts for sample-to-sample variation in total
# expression within each environment (acts as a sample normalization per env).
# Environment:Allele estimates the allelic ratio within each environment.
#
# Size factors are set to 1 because we only compare alleles within the same
# sample — inter-sample normalization is not appropriate here.
# ------------------------------------------------------------------------------

dds2 <- DESeqDataSetFromMatrix(
  countData = counts2,
  colData   = meta2,
  design    = ~ Environment + Environment:sampleID + Environment:Allele
)
sizeFactors(dds2) <- rep(1, nrow(meta2))
dds2 <- DESeq(dds2)

# Environment-specific allelic ratios (log2 Alt/Ref in each environment)
mo_ratio <- results(dds2, name = "EnvironmentMO.AlleleAlt")
nc_ratio <- results(dds2, name = "EnvironmentNC.AlleleAlt")

# GxE interaction: tests whether the Alt/Ref ratio DIFFERS between environments
ge <- results(dds2, contrast = list("EnvironmentMO.AlleleAlt",
                                     "EnvironmentNC.AlleleAlt"))

# Save DESeq2 objects for downstream analyses (co-expression VST extraction, etc.)
save(mo_ratio, nc_ratio, ge, dds2, meta2, counts2,
     file = file.path(DATA_DIR, "data/DESeq2_results.rdata"))
message("DESeq2 results saved to data/DESeq2_results.rdata")

# ------------------------------------------------------------------------------
# BUILD RESULTS TABLE (toPlot)
#
# ORIGINAL BUG: toPlot was constructed at the bottom of the script but used
# near the top to define sig_ge and bg_genes. Fixed here by constructing
# toPlot immediately after results extraction.
# Also fixed: padj was missing from the original join; now included.
# ------------------------------------------------------------------------------

toPlot <- full_join(
  full_join(
    mo_ratio %>%
      as_tibble(rownames = "GeneID") %>%
      dplyr::select(GeneID, MO_Log2FC = log2FoldChange, MO_padj = padj),
    nc_ratio %>%
      as_tibble(rownames = "GeneID") %>%
      dplyr::select(GeneID, NC_Log2FC = log2FoldChange, NC_padj = padj)
  ),
  ge %>%
    as_tibble(rownames = "GeneID") %>%
    dplyr::select(GeneID, pvalue, padj)
)

# Write full GxE results (MO_Log2FC, MO_padj, NC_Log2FC, NC_padj, pvalue, padj)
# MO_padj / NC_padj are the per-environment allele-ratio FDR values; used to
# define the G (constitutive ASE) gene set below.
write_tsv(toPlot, file.path(DATA_DIR, "DEG_GxE_results.txt"))
message("GxE results written to data/DEG_GxE_results.txt")

# Write full DESeq2 GxE output (baseMean, log2FoldChange, lfcSE, stat, pvalue, padj).
# Used by scripts 5, 6, 9, and null_results scripts.  Complements DEG_GxE_results.txt
# (which has MO/NC per-environment LFCs) with baseMean and the interaction stat.
write_tsv(as_tibble(ge, rownames = "GeneID"),
          file.path(DATA_DIR, "GxE_allele_specific_test_results.txt"))
message("Full GxE DESeq2 results written to data/GxE_allele_specific_test_results.txt")

# ------------------------------------------------------------------------------
# GENE SET DEFINITIONS
#
# GxE gene set: padj < padj_thresh AND in pre-filtered universe
#
# Background gene set (pre-filter):
#   Genes that (a) were tested in the DESeq2 model (have a pvalue) AND
#   (b) have >10 counts in >5 sample×allele columns.
#   This ensures background genes had sufficient expression to be realistically
#   testable for allele-specific effects.
#   This filter was implicitly applied in the original analysis; it is now
#   explicitly documented.
#
# NOTE: DESeq2's independent filtering assigns NA padj to low-expressed genes
# that fail its adaptive count threshold. These are treated as padj = 1 for
# threshold purposes (they are genuinely non-significant; median pvalue ~0.84).
# Their presence does not affect the GxE set; they join the background.
# ------------------------------------------------------------------------------

# Genes that were tested (have a valid pvalue)
tested_genes <- toPlot %>%
  drop_na(pvalue) %>%
  pull(GeneID)

# Explicit pre-filter: >10 counts in >5 sample×allele columns
enough_counts <- counts2_long %>%
  group_by(GeneID) %>%
  summarise(n_above_10 = sum(Count > 10)) %>%
  filter(n_above_10 > 5) %>%
  pull(GeneID)

# Universe = tested AND passes count filter
universe <- intersect(tested_genes, enough_counts)
message(sprintf("Pre-filtered universe: %d genes (%d tested, %d pass count filter)",
                length(universe), length(tested_genes), length(enough_counts)))

# GxE genes: padj < threshold within universe
# NA padj → treat as 1.0 (failed independent filtering = non-significant)
padj_vec <- toPlot$padj
padj_vec[is.na(padj_vec)] <- 1.0
names(padj_vec) <- toPlot$GeneID

sig_ge  <- universe[padj_vec[universe] < padj_thresh]
bg_genes <- setdiff(universe, sig_ge)

message(sprintf("GxE genes (padj < %.2f, pre-filtered): %d",
                padj_thresh, length(sig_ge)))
message(sprintf("Background genes (pre-filtered, not GxE): %d",
                length(bg_genes)))

# Sanity check: all GxE genes should be in universe
stopifnot(all(sig_ge %in% universe))

# Write gene sets
write.table(sig_ge,   file.path(DATA_DIR, "GxE_gene_IDs.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(bg_genes, file.path(DATA_DIR, "background_gene_IDs.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
message("Gene sets written.")

# ------------------------------------------------------------------------------
# ENVIRONMENT GENE SET
# ------------------------------------------------------------------------------

# Environment main effect (summed allele expression, ignores ASE)
dds_e <- DESeqDataSetFromMatrix(
  countData = counts2,
  colData   = meta2,
  design    = ~ Environment
)
dds_e  <- DESeq(dds_e)
e_res  <- results(dds_e)
sig_e  <- rownames(e_res)[!is.na(e_res$pvalue) & e_res$pvalue < 0.05]
write_tsv(as_tibble(e_res, rownames = "GeneID"), file.path(DATA_DIR, "DEG_Env_results.txt"))

# ------------------------------------------------------------------------------
# G GENE SET (constitutive allele-specific expression)
#
# G genes have a consistent allelic imbalance in BOTH environments — Mo17/B73
# ratio is shifted the same direction in MO and NC.
#
# Definition (stricter than the original raw-p < 0.05 approach):
#   (1) FDR < padj_thresh (0.1) in BOTH the MO-specific and NC-specific allele
#       ratio tests (MO_padj and NC_padj columns in toPlot).
#   (2) Same sign of log2FC in both environments — guards against genes where
#       independent filtering happens to yield low padj in both environments
#       even when allele ratios point in opposite directions, which would
#       actually be a GxE gene, not a G gene.
#   (3) Restricted to the pre-filtered universe for comparability with GxE.
#
# NA padj → treated as 1.0, consistent with the GxE definition above.
# ------------------------------------------------------------------------------

mo_padj_vec <- toPlot$MO_padj; mo_padj_vec[is.na(mo_padj_vec)] <- 1.0
nc_padj_vec <- toPlot$NC_padj; nc_padj_vec[is.na(nc_padj_vec)] <- 1.0
mo_lfc_vec  <- toPlot$MO_Log2FC
nc_lfc_vec  <- toPlot$NC_Log2FC
names(mo_padj_vec) <- names(nc_padj_vec) <-
  names(mo_lfc_vec) <- names(nc_lfc_vec) <- toPlot$GeneID

sig_g <- universe[
  mo_padj_vec[universe] < padj_thresh &
  nc_padj_vec[universe] < padj_thresh &
  !is.na(mo_lfc_vec[universe]) &
  !is.na(nc_lfc_vec[universe]) &
  sign(mo_lfc_vec[universe]) == sign(nc_lfc_vec[universe])
]

message(sprintf("G genes (padj < %.2f in both envs, same direction, pre-filtered): %d",
                padj_thresh, length(sig_g)))
message(sprintf("  GxE ∩ G overlap: %d genes (%.1f%% of GxE, %.1f%% of G)",
                length(intersect(sig_ge, sig_g)),
                100 * length(intersect(sig_ge, sig_g)) / max(length(sig_ge), 1),
                100 * length(intersect(sig_ge, sig_g)) / max(length(sig_g),  1)))

write.table(sig_g, file.path(DATA_DIR, "G_gene_IDs.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
message("G gene IDs written to data/G_gene_IDs.txt")

# ------------------------------------------------------------------------------
# FIGURE: GxE scatter (MO_Log2FC vs NC_Log2FC)
# ------------------------------------------------------------------------------

pthresh_plot <- padj_thresh  # colour by padj < 0.1

p_scatter <- ggplot(mapping = aes(MO_Log2FC, NC_Log2FC)) +
  geom_point(
    data = toPlot %>% filter(is.na(padj) | padj >= pthresh_plot),
    size = 0.5, color = "gray"
  ) +
  geom_point(
    data = toPlot %>%
      filter(!is.na(padj) & padj < pthresh_plot) %>%
      arrange(desc(padj)) %>%
      mutate(neg_log_padj = pmin(-log10(padj), 3)),
    size = 0.5,
    mapping = aes(color = neg_log_padj)
  ) +
  scale_color_continuous(
    breaks = seq(1, 3, 0.5),
    labels = c("1.0", "1.5", "2.0", "2.5", "\u2265 3.0"),
    name   = expression(-log[10](p[adj]))
  ) +
  theme_classic() +
  labs(
    x        = expression(paste(Log[2](Alt/Ref), " in MO")),
    y        = expression(paste(Log[2](Alt/Ref), " in NC")),
    title    = sprintf("GxE-ASE genes (padj < %.1f, n = %d)",
                       pthresh_plot, length(sig_ge)),
    subtitle = sprintf("Background: %d pre-filtered genes", length(bg_genes))
  ) +
  coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5))

ggsave(file.path(FIGURES_DIR, "gxe_scatter_padj01.pdf"), p_scatter,
       width = 5, height = 4.5)
message("Scatter plot saved.")
