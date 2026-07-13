library(tidyverse)
library(patchwork)

results = read_tsv("results/g_gene_enrichment.tsv")
nG = unique(results$n_group[results$group == "G-only"])
nGE = unique(results$n_group[results$group == "GxE-only"])
nGGE = unique(results$n_group[results$group == "GxE+G"])

OR_fig = ggplot(results %>% 
         filter(metric != "Large indels in promoter") %>%
           mutate(metric = ifelse(metric == "Small indel × DAP peak", "Indels disrupt DAP peak", metric),
                  metric = ifelse(metric == "SNP × DAP peak", "SNPs disrupt DAP peak", metric)) %>%
           mutate(metric = factor(metric, levels=unique(.$metric)[c(1:5)])),
       aes(group, fisher_OR, fill=group)) +
  geom_col(alpha=0.3) + 
  facet_wrap(~metric, nrow=1) +
  theme_classic() +
  scale_fill_manual(labels=c(sprintf("Genotype effect alone (n=%i)", nG), 
                             sprintf("GxE effect alone (n=%i)", nGE),
                             sprintf("Genotype and GxE effect (n=%i)", nGGE)), 
                    values=c("red", "blue", "purple"),
                    name="Genes with significant:") +
  labs(y="Odds Ratio") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  guides(fill = guide_legend(title.position="top"))

ggsave("figures/4_compare_G_GxE.pdf", width=10, height=4)
ggsave("figures/4_compare_G_GxE.png", width=10, height=4)

# pvalue_fig = ggplot(results %>% 
#                       filter(metric != "Large indels in promoter") %>%
#                       mutate(metric = ifelse(metric == "Indel × DAP peak", "Indels disrupt DAP peak", metric)) %>%
#                       mutate(metric = factor(metric, levels=unique(.$metric)[c(1,2,3)])),
#                 aes(group, -log10(mwu_p), fill=group)) +
#   geom_col(alpha=0.3) + 
#   facet_wrap(~metric) +
#   theme_classic() +
#   scale_fill_manual(labels=c("Genotype effect alone (n=1,586)", 
#                              "GxE effect alone (n=132)",
#                              "Genotype and GxE effect (n=116)"), 
#                     values=c("red", "blue", "purple"),
#                     name="Genes with significant:") +
#   labs(y=expression(-log[10]("p value"))) +
#   theme(axis.text.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         legend.position = "bottom") +
#   guides(fill = guide_legend(title.position="top"))
