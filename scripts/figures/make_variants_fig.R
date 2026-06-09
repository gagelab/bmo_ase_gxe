library(tidyverse)
library(patchwork)

sv = read_tsv("results/sv_per_gene.tsv")

# Make clean tibble for plotting distribution and testing MW
sv_continuous = sv %>%
  pivot_longer(contains("n_"), names_to = "Category", values_to = "Count",
               names_prefix = "n_") %>%
  filter(Category != "sv_total") %>%
  mutate(Category = case_when(Category == "small_indel" ~ "Small Indels",
                              Category == "large_indel" ~ "Large Indels",
                              Category == "SNP" ~ "SNPs")) %>%
  mutate(Count = case_when(Category == "Large Indels" & Count > 4 ~ 4,
                           Category == "Small Indels" & Count > 20 ~ 20,
                           Category == "SNPs" & Count > 50 ~ 50,
                           TRUE ~ Count))

# Make clean tibble for processing binary proportion and fisher test
sv_binary = sv_continuous %>%
  group_by(Category, is_gxe) %>%
  summarise(Proportion = sum(Count > 0) / n()) %>%
  ungroup()

# Put test statistics into a tibble
test_results = tibble()
for(type in unique(sv_continuous$Category)){
  # Subset to just the SV type we're focusing on
  #  and then make a table of counts for Fisher Test
  test_counts = sv_continuous %>%
    filter(Category == type) %>%
    mutate(has_var = Count > 0) %>%
    dplyr::select(is_gxe, has_var) %>%
    table()
  test_fisher = fisher.test(test_counts, alternative = "greater")

  # Subset to the SV we're focusing on and run MW test
  test_dist = sv_continuous %>%
    filter(Category == type) %>%
    dplyr::select(is_gxe, Count)
  test_mw = wilcox.test(test_dist$Count[test_dist$is_gxe],
                     test_dist$Count[!test_dist$is_gxe], alternative="greater")
  
  # Put results into a tibble
  test_results = bind_rows(
    test_results,
    tibble(Category = type, 
           OddsRatio = test_fisher$estimate, 
           Fisher_p=test_fisher$p.value,
           MW_p=test_mw$p.value)
  )
}


# Bar plot showing relative proportion of variant in GxE vs Background
fisher_bars = ggplot(sv_binary, aes(x=is_gxe, y=Proportion, fill=is_gxe)) +
  geom_col(alpha=0.3, position="dodge") +
  geom_text(x=-Inf, y=Inf, hjust=-0.1, vjust=1.5,
            mapping=aes(label=sprintf("Fisher p=%0.3f\nOR: %1.2f", Fisher_p, OddsRatio)), 
            inherit.aes = FALSE, data=test_results,
            size=3) +
  facet_wrap(~Category) +
  theme_classic() +
  # Remove all the X axis text etc - not needed
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  # Make extra space for the p-value text
  scale_y_continuous(expand=expansion(mult=c(0,0.25))) +
  labs(y="Proportion of genes\nwith \u22651 variant")
mw_dists = ggplot(sv_continuous,
                  aes(Count, after_stat(density), fill=is_gxe)) +
  geom_histogram(alpha=0.3, position="identity", binwidth = 1) +
  geom_text(x=-Inf, y=Inf, hjust=-0.1, vjust=1.5,
            mapping=aes(label=sprintf("MW p=%0.1g", MW_p)), 
            inherit.aes = FALSE, data=test_results,
            size=3) +
  facet_wrap(~Category, scales="free") +
  theme_classic() +
  labs(y="Density", x="Number of variants") +
  # Remove facet label strip since this will be on the bottom
  theme(strip.text = element_blank()) +
  # Make extra space for the p-value text
  scale_y_continuous(expand=expansion(mult=c(0,0.25)))

(fisher_bars / mw_dists) + 
  patchwork::plot_layout(guides="collect") & 
  scale_fill_manual(labels=c("Background Genes (n=11,389)", "GxE ASE Genes (n=248)"), 
                    values=c("black", "blue")) &
  theme(legend.position = "bottom",
        legend.title = element_blank())
ggsave("figures/2_variants_fig.pdf", width=6, height=8)
ggsave("figures/2_variants_fig.png", width=6, height=8)

