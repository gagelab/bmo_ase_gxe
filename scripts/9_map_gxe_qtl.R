################################################################################
# 7_map_gxe_qtl.R
#
# Two-stage GxE QTL scan across the IBM (Intermated B73 x Mo17) population
# for 19 NAM phenotypic traits.
#
# MODEL RATIONALE:
#   A naive lm(trait ~ Genotype * env) applied to RIL data confounds the
#   polygenic genetic background (each line is genetically unique) with the
#   residual, inflating the GxE F-test denominator ~2x and driving p-values
#   toward 1 (median ~0.82 vs expected 0.5 under null).
#
#   Stage 1 — Block adjustment (once per trait):
#     Within each environment, regress out incomplete-block effects to remove
#     within-environment spatial variation. Preserves grand mean.
#
#   Stage 2 — Locus scan:
#     G main effect:   lm(block-adj line mean ~ Genotype)
#                      Error = within-genotype-class line variance.
#     G×E interaction: F-test of Genotype:env columns added to a null model
#                      that already includes Entry_ID (absorbs polygenic main
#                      effect) and env. Uses incremental QR decomposition so
#                      the null model matrix is factored once per trait rather
#                      than once per marker.
#
# Reads:
#   data/IBM_recomb_bins_fromGBS.tsv   (output of 6_format_gbs_genotypes.R)
#   data/IBM_Name_M00_Z017_decoder.txt
#   data/NAM_all_traits.txt
#
# Writes:
#   results/IBM_GxE_results.tsv
#
# Runtime: ~20–40 min with parallel workers
################################################################################

library(tidyverse)
library(future.apply)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE        <- "./"
RESULTS     <- file.path(BASE, "results")
DATA    <- file.path(BASE, "data")   # large shared data files
BINS_FILE   <- file.path(DATA, "IBM_recomb_bins_fromGBS.tsv")
DECODER     <- file.path(DATA, "IBM_Name_M00_Z017_decoder.txt")
PHENO_FILE  <- file.path(DATA, "NAM_all_traits.txt")
OUT_FILE    <- file.path(RESULTS, "IBM_GxE_results.tsv")

dir.create(RESULTS, showWarnings = FALSE)

# ── Load genotype bins ─────────────────────────────────────────────────────────
message("Loading recombination bins ...")
geno <- read_tsv(BINS_FILE, show_col_types = FALSE)
message(sprintf("  %d bins x %d samples", nrow(geno), ncol(geno) - 3))

# ── Load and join phenotypic data ─────────────────────────────────────────────
message("Loading phenotypic data ...")
name_decoder <- read_tsv(DECODER, show_col_types = FALSE) %>%
  dplyr::select(name, synonyms) %>%
  mutate(
    Entry_ID  = sub(".*(M[0-9]{4})", "\\1", name),
    ZE_ID     = sub(".*(Z[0-9]{3}E[0-9]{4}).*", "\\1", synonyms),
    pop       = as.numeric(sub("Z([0-9]{3})E[0-9]{4}", "\\1", ZE_ID)),
    entry_num = as.numeric(sub("Z[0-9]{3}E([0-9]{4})", "\\1", ZE_ID))
  ) %>%
  drop_na() %>%
  dplyr::select(Entry_ID, pop, entry_num)

pheno <- read_tsv(PHENO_FILE, show_col_types = FALSE) %>%
  filter(pedigree == "IBM") %>%
  left_join(name_decoder, by = c("pop", "entry_num"))

trait_cols <- colnames(pheno)[14:32]   # 19 phenotypic traits

message(sprintf("  %d phenotype rows | %d traits | %d environments",
                nrow(pheno), length(trait_cols), n_distinct(pheno$env)))
message(sprintf("  IBM inbreds matched to genotypes: %d",
                sum(unique(pheno$Entry_ID) %in% colnames(geno))))

# ── Stage 1: Block-adjust within environments ─────────────────────────────────
message("\nStage 1: Block-adjusting phenotypes ...")
pheno_adj <- pheno %>%
  mutate(block_id = paste(block, rep, field, sep = "-")) %>%
  group_by(env) %>%
  group_modify(function(env_data, key) {
    for (tr in trait_cols) {
      y <- env_data[[tr]]
      adj_name <- paste0(tr, "_adj")
      if (sum(!is.na(y)) < 10) {
        env_data[[adj_name]] <- y
        next
      }
      m <- tryCatch(
        lm(as.formula(paste(tr, "~ block_id")), data = env_data,
           na.action = na.exclude),
        error = function(e) NULL
      )
      env_data[[adj_name]] <- if (is.null(m)) y else
        residuals(m) + mean(y, na.rm = TRUE)
    }
    env_data
  }) %>%
  ungroup()

adj_cols <- paste0(trait_cols, "_adj")
message("Stage 1 complete.")

# ── Stage 2: GxE locus scan ───────────────────────────────────────────────────
message("\nStage 2: GxE locus scan ...")
plan(multisession, workers = max(1L, parallelly::availableCores() - 1L))

# Pre-extract genotype matrix (markers x samples)
geno_ids  <- colnames(geno)[-(1:3)]   # all sample columns
geno_mat  <- as.matrix(geno[, geno_ids])   # "B73"/"Mo17"/NA
rownames(geno_mat) <- seq_len(nrow(geno))
n_markers <- nrow(geno)

# Incremental F-test via QR decomposition.
# Given the QR of the null model (Q0) and its RSS + df, tests whether adding
# X1 (the Genotype:env columns) significantly reduces RSS.
f_test_incremental <- function(Q0, rss_null, df_null, y, X1) {
  X1r   <- X1 - Q0 %*% crossprod(Q0, X1)   # project out null column space
  qr1   <- qr(X1r)
  rank1 <- qr1$rank
  if (rank1 == 0L) return(list(GEF = NA_real_, GEp = NA_real_))
  Q1    <- qr.Q(qr1)[, seq_len(rank1), drop = FALSE]
  rss_reduction <- sum(crossprod(Q1, y)^2)
  df_alt  <- df_null - rank1
  rss_alt <- rss_null - rss_reduction
  if (df_alt <= 0 || rss_alt <= 0) return(list(GEF = NA_real_, GEp = NA_real_))
  GEF <- (rss_reduction / rank1) / (rss_alt / df_alt)
  list(GEF = GEF, GEp = pf(GEF, rank1, df_alt, lower.tail = FALSE))
}

scan_one_trait <- function(trait) {
  adj_trait  <- paste0(trait, "_adj")
  test_pheno <- pheno_adj %>%
    dplyr::select(Entry_ID, env, all_of(adj_trait)) %>%
    drop_na() %>%
    mutate(env = factor(env))

  envs  <- levels(test_pheno$env)
  n_env <- length(envs)

  # Factor the null model ONCE for this trait
  X0_full   <- model.matrix(~ Entry_ID + env, data = test_pheno)
  qr0_full  <- qr(X0_full)
  Q0_full   <- qr.Q(qr0_full)
  y_full    <- test_pheno[[adj_trait]]
  resid0    <- y_full - Q0_full %*% crossprod(Q0_full, y_full)
  rss_null0 <- sum(resid0^2)
  df_null0  <- length(y_full) - qr0_full$rank

  entry_row_map <- split(seq_len(nrow(test_pheno)), test_pheno$Entry_ID)
  line_means <- test_pheno %>%
    group_by(Entry_ID) %>%
    summarise(y = mean(.data[[adj_trait]], na.rm = TRUE), .groups = "drop")

  out <- vector("list", n_markers)
  for (g in seq_len(n_markers)) {
    marker_geno <- geno_mat[g, ]
    shared_ids  <- intersect(names(marker_geno)[!is.na(marker_geno)],
                             names(entry_row_map))
    geno_vec    <- marker_geno[shared_ids]
    if (length(unique(geno_vec)) < 2) next

    # G main effect (two-sample F on per-line means)
    lm_data <- line_means %>%
      filter(Entry_ID %in% shared_ids) %>%
      mutate(Genotype = marker_geno[Entry_ID])
    grp <- split(lm_data$y, lm_data$Genotype)
    if (length(grp) < 2 || any(lengths(grp) < 2)) next

    G_result <- tryCatch({
      n1 <- length(grp[["B73"]]);   n2 <- length(grp[["Mo17"]])
      m1 <- mean(grp[["B73"]]);     m2 <- mean(grp[["Mo17"]])
      gm <- (n1 * m1 + n2 * m2) / (n1 + n2)
      ss_b <- n1 * (m1 - gm)^2 + n2 * (m2 - gm)^2
      ss_w <- sum((grp[["B73"]] - m1)^2) + sum((grp[["Mo17"]] - m2)^2)
      df_w <- n1 + n2 - 2
      GF   <- (ss_b / 1) / (ss_w / df_w)
      list(GF = GF, Gp = pf(GF, 1, df_w, lower.tail = FALSE))
    }, error = function(e) list(GF = NA_real_, Gp = NA_real_))

    # G×E interaction
    rows_keep   <- unlist(entry_row_map[shared_ids])
    geno_per_row <- marker_geno[test_pheno$Entry_ID[rows_keep]]
    env_per_row  <- test_pheno$env[rows_keep]
    ge_tab       <- table(geno_per_row, env_per_row)
    if (nrow(ge_tab) < 2 || ncol(ge_tab) < n_env || any(ge_tab < 5)) next

    GxE_result <- tryCatch({
      if (length(rows_keep) == length(y_full)) {
        Q0_use      <- Q0_full
        y_use       <- y_full
        rss_null_use <- rss_null0
        df_null_use  <- df_null0
      } else {
        X0_sub   <- X0_full[rows_keep, , drop = FALSE]
        X0_sub   <- X0_sub[, colSums(abs(X0_sub)) > 0, drop = FALSE]
        qr0_sub  <- qr(X0_sub)
        Q0_use   <- qr.Q(qr0_sub)
        y_use    <- y_full[rows_keep]
        resid_s  <- y_use - Q0_use %*% crossprod(Q0_use, y_use)
        rss_null_use <- sum(resid_s^2)
        df_null_use  <- length(y_use) - qr0_sub$rank
      }
      is_mo17 <- as.numeric(geno_per_row == "Mo17")
      X1      <- model.matrix(~ 0 + env_per_row) * is_mo17
      f_test_incremental(Q0_use, rss_null_use, df_null_use, y_use, X1)
    }, error = function(e) list(GEF = NA_real_, GEp = NA_real_))

    out[[g]] <- tibble(
      chr   = geno$chr[g],
      start = geno$start[g],
      end   = geno$end[g],
      trait = trait,
      GF    = G_result$GF,
      Gp    = G_result$Gp,
      GEF   = GxE_result$GEF,
      GEp   = GxE_result$GEp
    )
  }
  bind_rows(out)
}

results <- future_lapply(trait_cols, scan_one_trait, future.seed = TRUE) %>%
  bind_rows()

message(sprintf("Scan complete. %d locus-trait tests.", nrow(results)))

# ── Save ───────────────────────────────────────────────────────────────────────
write_tsv(results, OUT_FILE)
message(sprintf("Results written to %s", OUT_FILE))
