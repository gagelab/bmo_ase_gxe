setwd("~/projects/bmo_ase/gxe_ase/")

library(tidyverse)
library(magrittr)
library(DESeq2)

pthresh = 0.1
niter=100

# TODO: Need to make sure sampling doesn't replicate

load("data/DESeq2_results.rdata")


draw_unique_sample <- function(population, k) {
  # This function draws a sample of size k from population, creates a key that
  #  is just the sorted sample set, pasted together. Then it checks an 
  #  environment named 'seen' for that key. If the key doesn't exist, it returns
  #  this sample, otherwise it contineus to try sampling until it gets an unseen
  #  sample.
  # Seen can be initialized outside the function with:
  #    seen <- new.env(hash = TRUE, parent = emptyenv())
  repeat {
    samp <- sample(population, k)
    key <- paste(sort(samp), collapse = ",")
    
    if (!exists(key, envir = seen, inherits = FALSE)) {
      seen[[key]] <- TRUE
      return(samp)
    }
  }
}

NC_samples = meta2 %>% filter(Environment == "NC") %>% pull(LibraryID) %>% unique()
MO_samples = meta2 %>% filter(Environment == "MO") %>% pull(LibraryID) %>% unique()
subsetting_results = tibble(
  N=numeric(),
  iter=numeric(),
  nSig=numeric(),
  median_p = numeric()
)


set.seed(2126)
for(N in 3:14){
  seen <- new.env(hash = TRUE, parent = emptyenv())
  for(iter in 1:niter){
    print(sprintf("Sample size: %i   Iteration: %i", N, iter))
    
    # Only run once for n=14, since no subsampling possible
    if(N == 14 & iter > 1){ next }
    if(N == 13 & iter > 14){ next } # 14 choose 13 = 14
    if(N == 12 & iter > 91){ next } # 14 choose 12 = 91
    
    # Sample and subset data to N samples
    NC_subset = draw_unique_sample(NC_samples, N)
    # Use setequal() to make sure NC_subset doesn't occur in previous samples
    
    MO_subset = draw_unique_sample(MO_samples, N)
    
    # Create metadata and counts tables
    meta_subset = meta2 %>% 
      filter(LibraryID %in% c(NC_subset, MO_subset)) %>%
      dplyr::select(-c(sampleID))
    sampleID = meta_subset %>%
      distinct(Environment, LibraryID) %>%
      group_by(Environment) %>%
      mutate(sampleID = 1:n(),
             sampleID = factor(sampleID))
    meta_subset = left_join(meta_subset, sampleID,
                            by=c("Environment", "LibraryID"))
    
    keep_counts = paste(c(NC_subset, MO_subset), collapse="|")
    count_subset = counts2[,grepl(keep_counts, colnames(counts2))]
    
    # Check the data match up
    if(! all(meta_subset$RowID == colnames(count_subset))){
      stop("Meta rows don't match Count colnames")
    }
    
    # Make DESeq object and test
    suppressMessages(
      dds_sub <- DESeqDataSetFromMatrix(count_subset, meta_subset,
                                       design = ~ Environment + Environment:sampleID + Environment:Allele)
    )
    sizeFactors(dds_sub) <- rep(1, nrow(meta_subset))
    dds_sub = DESeq(dds_sub,
                    parallel = TRUE,
                    quiet=TRUE)
    ge = results(dds_sub, contrast=list("EnvironmentMO.AlleleAlt", "EnvironmentNC.AlleleAlt"))
    n_sig = sum(ge$padj < pthresh, na.rm = TRUE)
    ntest = sum(!is.na(ge$pvalue))
    
    # Add results to results table
    subsetting_results = bind_rows(
      subsetting_results,
      tibble(N=N, iter=iter, nSig = n_sig, nTest = ntest,
             median_p = median(ge$pvalue[ge$padj < pthresh], na.rm=TRUE))
    )
  }
}
print("All Done.")
write_tsv(subsetting_results, "results/n_siginificant_subsampling.txt")

# Figure showing results
ggplot(subsetting_results, aes(N, nSig, group=N)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(height=0, width=0.1, alpha=0.25) +
  labs(x="Number of samples from each environment",
       y="Number of significant GxE genes") +
  theme_classic()
ggsave("figures/1_subsample_test_GxE.png", width=8, height=5)
# Check number of genes as a proportion of all genes
message("Summary of number of significant genes. Median per sample size (3-14)")
subsetting_results %>% mutate(prop = nSig/nTest) %>% 
  group_by(N) %>%
  summarise(prop = median(prop)) %>%
  pull(prop) %>% summary()
message("Summary of percent of significant genes. Median per sample size (3-14)")
subsetting_results %>%  
  group_by(N) %>%
  summarise(nSig = median(nSig)) %>%
  pull(nSig) %>% summary()

# Figure similar to the one above, but showing proportion significant genes
# ggplot(subsetting_results, aes(N, nSig/nTest, group=N)) +
#   geom_boxplot(outliers = FALSE) +
#   geom_jitter(height=0, width=0.1, alpha=0.25) +
#   labs(x="Number of samples from each environment",
#        y="Proportion of significant GxE genes") +
#   theme_classic()
