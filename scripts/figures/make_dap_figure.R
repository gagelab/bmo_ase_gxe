library(tidyverse)
library(patchwork)

dap_gene = read_tsv("results/dap_seq_per_gene.tsv")
dap_tf = read_tsv("results/dap_seq_per_tf.tsv")

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
          axis.title.x = element_blank()) +
    labs(y="Proportion genes with \u22651\ngenotype-specific DAP peaks ") +
    # Make extra space for the p-value text
    scale_y_continuous(expand=expansion(mult=c(0,0.25)))
  )

(dap_mw = dap_gene %>%
  mutate(n_diff = ifelse(n_diff > 10, 10, n_diff)) %>%
  ggplot(aes(n_diff, after_stat(density), fill=is_gxe)) +
    geom_histogram(alpha=0.3, position="identity", binwidth = 1) +
    geom_text(x=-Inf, y=Inf, hjust=-0.1, vjust=1.5,
              mapping=aes(label=sprintf("MW p=%0.1g", MW_p)), 
              inherit.aes = FALSE, data=test_results,
              size=3) +
    theme_classic() +
    scale_x_continuous(breaks=seq(0, 10)) +
    labs(y="Density", x="Number of genotype-specific DAP peaks") +
    # Make extra space for the p-value text
    scale_y_continuous(expand=expansion(mult=c(0,0.2)))
)

or_fig = dap_tf %>%
  arrange(log2OR_diff) %>%
  drop_na(log2OR_diff) %>%
  mutate(idx = 1:n()) %>%
ggplot(aes(idx, log2OR_diff)) +
  geom_col(color="white", fill="black", linewidth = .1) +
  theme_classic() +
  labs(x="Ordered transcription factors",
       y=expression(log[2]("Odds Ratio")))


(dap_fisher + dap_mw + guide_area() + or_fig ) + 
  patchwork::plot_layout(guides="collect", heights=c(1, .5), nrow = 2) +
  plot_annotation(tag_levels = "A") & 
  scale_fill_manual(labels=c("Background Genes (n=11,389)", "GxE ASE Genes (n=248)"), 
                    values=c("black", "blue")) &
  theme(legend.position = "right",
        legend.title = element_blank())
ggsave("figures/3_dap_fig.pdf", width=6, height=4.5)
ggsave("figures/3_dap_fig.png", width=6, height=4.5)

