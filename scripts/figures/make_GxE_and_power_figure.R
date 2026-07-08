library(tidyverse)
library(patchwork)

subsetting_results = read_tsv("results/n_siginificant_subsampling.txt")
deseq_results = read_tsv("data/DEG_GxE_results.txt")

# GxE genes figures
scatter_fig = deseq_results %>%
  mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
  arrange(-padj) %>%
  mutate(padj = ifelse(padj > .1, NA, padj),
         alpha = ifelse(is.na(padj), .2, 1)) %>%
ggplot(aes(MO_Log2FC, NC_Log2FC, color=-log10(padj), alpha=alpha)) +
  geom_vline(xintercept=0, lty=2) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(size=.75) +
  theme_classic() +
  labs(x=expression(paste(Log[2](Mo17/B73), " in MO")),
       y=expression(paste(Log[2](Mo17/B73), " in NC")),
       color=expression(-log[10](p[adj]))) +
  lims(x=c(-2,2),
       y=c(-2,2)) +
  scale_color_viridis_c(option="plasma", end=0.9) +
  coord_fixed() +
  guides(alpha = "none", color="none")

deseq_results %>%
  filter(padj < 0.1) %>%
  mutate(single_env = case_when(MO_padj < 0.1 & NC_padj > 0.1 ~ TRUE,
                                MO_padj > 0.1 & NC_padj < 0.1 ~ TRUE,
                                TRUE ~ FALSE)) %>%
  dplyr::select(GeneID, single_env, padj, MO_Log2FC, NC_Log2FC) %>%
  pivot_longer(4:5, names_to = "Loc", values_to = "Log2FC") %>%
  arrange(-padj) %>%
  ggplot(aes(Loc, Log2FC, group=GeneID, color=-log10(padj))) +
  # ggplot(aes(Log2FC)) +
  # geom_histogram()
  geom_line(alpha=1) +
  geom_jitter(height = 0, width=0.05) +
  scale_color_viridis_c(option="plasma", end=0.9)

sig_genes <- deseq_results %>%
  filter(padj < 0.1) %>%
  mutate(delta = NC_Log2FC - MO_Log2FC,
         avg_log2fc = (MO_Log2FC + NC_Log2FC) / 2,
         sign_change = sign(MO_Log2FC) != sign(NC_Log2FC))

signchange_fig = sig_genes %>%
  mutate(category = ifelse(sign_change, "Sign change", "Magnitude change")) %>%
  dplyr::select(GeneID, category, padj, MO_Log2FC, NC_Log2FC) %>%
  pivot_longer(c(MO_Log2FC, NC_Log2FC),
               names_to = "Loc", values_to = "Log2FC") %>%
  mutate(Loc = ifelse(Loc == "MO_Log2FC", "MO", "NC")) %>%
  ggplot(aes(x = Loc, y = Log2FC, group = GeneID,
             color = -log10(padj))) +
  geom_line(alpha = 0.5) +
  geom_point(size = 1, alpha=0.5) +
  facet_wrap(~category) +
  scale_color_viridis_c(option = "plasma", end = 0.9) +
  theme_classic() +
  labs(x = "Environment",
       y = expression(Log[2](Mo17/B73)),
       color = expression(-log[10](p[adj])))

# Test for difference between categories in abs magnitude of location change:
sig_genes %>%
  mutate(category = ifelse(sign_change, "Sign change", "Magnitude change")) %>%
  dplyr::select(GeneID, category, padj, MO_Log2FC, NC_Log2FC) %>%
  mutate(diff = abs(MO_Log2FC - NC_Log2FC)) %>% wilcox.test(diff ~ category, data=.)

# Subsetting figure
subset_fig = ggplot(subsetting_results, aes(N, nSig, group=N)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(height=0, width=0.1, alpha=0.25) +
  labs(x="Number of samples from each environment",
       y="Significant GxE genes") +
  theme_classic()

scatter_fig + guide_area() + signchange_fig + plot_spacer() + subset_fig +
  patchwork::plot_layout(guides = "collect",
                         design = "111123333
                                   111123333
                                   111145555
                                   111145555") +
  plot_annotation(tag_levels = "A") &
  theme(legend.title.position = "left", legend.title = element_text(angle=90, hjust=0.5)) &
  guides(color = guide_colorbar(barwidth = unit(0.1, "in"), 
                               barheight = unit(1.5, "in")))
ggsave("figures/1_GxE_genes_and_subsample.png", width=10, height=4.5)
ggsave("figures/1_GxE_genes_and_subsample.pdf", width=10, height=4.5)
