library(tidyverse)
library(patchwork)

dap_gene = read_tsv("results/dap_seq_per_gene.tsv")
dap_tf = read_tsv("results/dap_seq_per_tf.tsv")
# dap_indel = read_tsv("results/indel_dap_overlap_per_gene.tsv")

test_fisher = dap_gene %>%
  mutate(any_diff = n_diff > 0) %>%
  dplyr::select(is_gxe, any_diff) %>%
  table() %>%
  fisher.test(alternative = "greater")
test_mw = wilcox.test(dap_gene$n_diff[dap_gene$is_gxe],
                      dap_gene$n_diff[!dap_gene$is_gxe],
                      alternative = "greater")
test_results = tibble(
  Fisher_p = test_fisher$p.value,
  OddsRatio = test_fisher$estimate,
  MW_p = test_mw$p.value
)

dap_binary = dap_gene %>%
  mutate(any_diff = n_diff > 0) %>%
  group_by(is_gxe) %>%
  summarise(prop_diff = sum(any_diff)/n())


(dap_fisher = ggplot(dap_binary,
                     aes(is_gxe, prop_diff, fill=is_gxe)) +
    geom_col(alpha=0.3) +
    geom_text(x=-Inf, y=Inf, hjust=-0.1, vjust=1.5,
              mapping=aes(label=sprintf("Fisher p=%0.1g\nOR: %1.2f", Fisher_p, OddsRatio)), 
              inherit.aes = FALSE, data=test_results,
              size=3) +
    theme_classic() +
    # Remove all the X axis text etc - not needed
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none") +
    labs(y="Proportion genes with\nallele-specific DAP peaks",
         title="A") +
    # Make extra space for the p-value text
    scale_y_continuous(expand=expansion(mult=c(0,0.25)))
  )

(dap_mw = dap_gene %>%
  mutate(n_diff = ifelse(n_diff > 5, 5, n_diff)) %>%
  ggplot(aes(n_diff, after_stat(density), fill=is_gxe)) +
    geom_histogram(alpha=0.3, position="identity", binwidth = 1) +
    geom_text(x=-Inf, y=Inf, hjust=-0.1, vjust=1.5,
              mapping=aes(label=sprintf("MW p=%0.1g", MW_p)), 
              inherit.aes = FALSE, data=test_results,
              size=3) +
    theme_classic() +
    theme(legend.position="none") +
    scale_x_continuous(breaks=seq(0, 5), labels = c(as.character(0:4), "\u22655")) +
    labs(y="Density", x="Number of genotype-specific\nDAP peaks",
         title="B") +
    # Make extra space for the p-value text
    scale_y_continuous(expand=expansion(mult=c(0,0.2)))
)

or_fig = dap_tf %>%
  arrange(log2OR_diff) %>%
  drop_na(log2OR_diff) %>%
  mutate(idx = 1:n()) %>%
ggplot(aes(idx, log2OR_diff, fill=is_gxe)) +
  geom_col(color="white", fill="black", linewidth = .1) +
  theme_classic() +
  labs(x="Ordered transcription factors",
       y=expression(log[2]("Odds Ratio")),
       title="C")


#### Figures for indel overlapping DAP peaks ####
dap_indel_table = dap_gene %>%
  dplyr::select(is_gxe, n_small_indel_genotype_spec) %>%
  mutate(n_small_indel_genotype_spec = n_small_indel_genotype_spec > 0) %>%
  table()

test_fisher = fisher.test(dap_indel_table, alternative = "greater")
test_mw = wilcox.test(dap_gene$n_small_indel_genotype_spec[dap_gene$is_gxe],
                      dap_gene$n_small_indel_genotype_spec[!dap_gene$is_gxe],
                      alternative = "greater")
indel_test_results = tibble(
  Fisher_p = test_fisher$p.value,
  OddsRatio = test_fisher$estimate,
  MW_p = test_mw$p.value
)

(dap_indel_or = dap_gene %>%
  group_by(is_gxe) %>%
  summarise(prop_specific_indel_peak = sum(n_small_indel_genotype_spec > 0)/n()) %>%
  ggplot(aes(is_gxe, prop_specific_indel_peak, fill=is_gxe)) +
    geom_col(alpha=0.3) +
    geom_text(x=-Inf, y=Inf, hjust=-0.1, vjust=1.5,
              mapping=aes(label=sprintf("Fisher p=%0.1g\nOR: %1.2f", Fisher_p, OddsRatio)), 
              inherit.aes = FALSE, data=indel_test_results,
              size=3) +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
    labs(y="Proportion genes with\nindels disrupting DAP peaks",
         title="D") +
    scale_y_continuous(expand=expansion(mult=c(0,0.25)))
)

(dap_indel_dist = dap_gene %>%
    mutate(n_genotype_spec = ifelse(n_small_indel_genotype_spec > 5, 5, n_small_indel_genotype_spec)) %>%
    ggplot(aes(n_genotype_spec, after_stat(density), fill=is_gxe)) +
    geom_histogram(alpha=0.3, binwidth = 1, position="identity") +
    geom_text(x=-Inf, y=Inf, hjust=-0.1, vjust=1.5,
              mapping=aes(label=sprintf("MW p=%0.1g", MW_p)), 
              inherit.aes = FALSE, data=indel_test_results,
              size=3) +
    theme_classic() +
    labs(x="Number of small indels overlapping\ngenotype-specific DAP peaks",
         y="Density",
         title="E") +
    scale_x_continuous(breaks=seq(0, 5), labels = c(as.character(0:4), "\u22655")) +
    scale_y_continuous(expand=expansion(mult=c(0,0.25))) 
)

#### Figures for SNPs overlapping DAP peaks ####
dap_snp_table = dap_gene %>%
  dplyr::select(is_gxe, n_snp_genotype_spec) %>%
  mutate(n_snp_genotype_spec = n_snp_genotype_spec > 0) %>%
  table()

test_fisher = fisher.test(dap_snp_table, alternative = "greater")
test_mw = wilcox.test(dap_gene$n_snp_genotype_spec[dap_gene$is_gxe],
                      dap_gene$n_snp_genotype_spec[!dap_gene$is_gxe],
                      alternative = "greater")
snp_test_results = tibble(
  Fisher_p = test_fisher$p.value,
  OddsRatio = test_fisher$estimate,
  MW_p = test_mw$p.value
)

(dap_snp_or = dap_gene %>%
    group_by(is_gxe) %>%
    summarise(prop_specific_snp_peak = sum(n_snp_genotype_spec > 0)/n()) %>%
    ggplot(aes(is_gxe, prop_specific_snp_peak, fill=is_gxe)) +
    geom_col(alpha=0.3) +
    geom_text(x=-Inf, y=Inf, hjust=-0.1, vjust=1.5,
              mapping=aes(label=sprintf("Fisher p=%0.1g\nOR: %1.2f", Fisher_p, OddsRatio)), 
              inherit.aes = FALSE, data=snp_test_results,
              size=3) +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
    labs(y="Proportion genes with\nSNPs disrupting DAP peaks",
         title="F") +
    scale_y_continuous(expand=expansion(mult=c(0,0.25)))
)

(dap_snp_dist = dap_gene %>%
    mutate(n_genotype_spec = ifelse(n_snp_genotype_spec > 5, 5, n_snp_genotype_spec)) %>%
    ggplot(aes(n_genotype_spec, after_stat(density), fill=is_gxe)) +
    geom_histogram(alpha=0.3, binwidth = 1, position="identity") +
    geom_text(x=-Inf, y=Inf, hjust=-0.1, vjust=1.5,
              mapping=aes(label=sprintf("MW p=%0.1g", MW_p)), 
              inherit.aes = FALSE, data=snp_test_results,
              size=3) +
    theme_classic() +
    labs(x="Number of SNPs overlapping\ngenotype-specific DAP peaks",
         y="Density",
         title="G") +
    scale_x_continuous(breaks=seq(0, 5), labels = c(as.character(0:4), "\u22655")) +
    scale_y_continuous(expand=expansion(mult=c(0,0.25)))
)


# Make row titles 
r1 = ggplot() +
  annotate("text",
           label="Genotype-specific DAP-seq peaks:",
           x=-0.13,
           y=0.5,
           size=5.5,
           hjust=0,
           vjust=0.5) +
  # scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  # scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void()  +
  coord_cartesian(xlim=c(0,1), ylim=c(0,1), clip="off")
  # theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

r2 = ggplot(aes(fill=is_gxe), data=dap_gene) +
  annotate("text",
           label="Small indels overlapping genotype-specific DAP-seq peaks:",
           x=-0.13,
           y=0.5,
           size=5.5,
           hjust=0,
           vjust=0.5) +
  # scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  # scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin=margin(6,0,0,0)) +
  coord_cartesian(xlim=c(0,1), ylim=c(0,1), clip="off")

r3 = ggplot(aes(fill=is_gxe), data=dap_gene) +
  annotate("text",
           label="SNPs overlapping genotype-specific DAP-seq peaks:",
           x=-0.13,
           y=0.5,
           size=5.5,
           hjust=0,
           vjust=0.5) +
  # scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  # scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin=margin(6,0,0,0)) +
  coord_cartesian(xlim=c(0,1), ylim=c(0,1), clip="off")


#### Compile figures together ####
layout = c("AAA
            BCD
            BCD
            BCD
            BCD
            BCD
            EEE
            FGH
            FGH
            FGH
            FGH
            FGH
            III
            JKL
            JKL
            JKL
            JKL
            JKL")  

( r1 +
  dap_fisher + dap_mw + or_fig +
  r2 + 
  dap_indel_or + dap_indel_dist + guide_area() +
  r3 +
  dap_snp_or + dap_snp_dist + plot_spacer()) +
  patchwork::plot_layout(guides = "collect", #ncol=3, nrow=2,
                         # heights = c(1,4,1,4),
                         design = layout) &
  # plot_annotation(tag_levels = "A") & 
  scale_fill_manual(labels=c("Background Genes (n=11,334)", "AxE Genes (n=219)"), 
                    values=c("black", "blue")) &
  theme(legend.position = "right",
        legend.title = element_blank(),
        plot.title.position = "plot")
ggsave("figures/3_dap_fig.pdf", width=10, height=9)
ggsave("figures/3_dap_fig.png", width=10, height=9)

