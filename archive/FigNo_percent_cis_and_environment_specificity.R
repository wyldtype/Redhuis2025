sapply(c("tidyverse"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Barkai_data_analysis/")

# Script to construct Figure 2: Networks in both species, bipartite graph, and divergent gene depletion in conserved blocks (probably)

# load single gene model dataset
load(file = "data_files/single_gene_models.RData")

# re-call sig if you want to change thresholds
p_thresh <- 0.05/(length(unique(spaldf$gene_name))*5)
eff_thresh <- 1
spaldf$sig <- spaldf$pvalue < p_thresh & abs(spaldf$effect_size) > eff_thresh & spaldf$isexpressed
sum(spaldf$sig)
spaldf %>% group_by(gene_name) %>% summarise(sig_gene = any(sig)) %>% select(sig_gene) %>% pull() %>% sum()
# re-creating genedf, where each row is one gene
genedf <-  spaldf %>% 
  pivot_wider(names_from = c("experiment", "coefficient"), values_from = c("effect_size", "pvalue", "sig"), id_cols = c("gene_name"))
table(genedf$gene_name) %>% table() # checking that every gene is only represented once

# Quantifying %cis
# For now we'll define %cis as the proportion of the effect size of each species-divergent gene that can be explained by hybrid divergence in each experiment
expdf <- spaldf %>% filter(experiment != "all") %>% 
  pivot_wider(id_cols = c(gene_name, experiment), 
              names_from = coefficient, 
              values_from = c(effect_size, sig)) %>% 
  filter(sig_species) %>% 
  mutate(pct_cis = effect_size_allele/effect_size_species)
ggplot(expdf, aes(x = effect_size_species, y = effect_size_allele)) + geom_point(aes(color = sig_allele))
# %cis for now is simply (allele effect size)/(species effect size). 
# It'll be positive if they're in the same direction
# It'll be > 1 if the hybrid showed a stronger effect than the parents (compensated, aka hybrid misexpressed)
ggplot(expdf, aes(x = pct_cis)) + geom_histogram() + geom_vline(xintercept = 1, color = "red") + ggtitle("% cis among species-divergent genes") + theme_classic()

# QC: %cis should be similar for the same gene in all 4 environments (but doesn't need to be, if a trans regulator is only expressed in certain environments or something)
# Null comparison: 4 random genes in the 4 environments
pct_cis_sd_same_gene <- expdf %>% group_by(gene_name) %>% summarise(pct_cis_sd = sd(pct_cis, na.rm = TRUE)) # many will be NA if they only had one experiment where they were divergent
pct_cis_sd_null_set <- expdf %>% mutate(random_group = sample(c(1:floor(nrow(expdf)/4)), nrow(expdf), replace = TRUE)) %>% group_by(random_group) %>% summarise(pct_cis_sd = sd(pct_cis, na.rm = TRUE))
qcdf <- tibble(sd = c(pct_cis_sd_same_gene$pct_cis_sd, pct_cis_sd_null_set$pct_cis_sd),
               type = c(rep("same_gene", length(pct_cis_sd_same_gene$pct_cis_sd)),
                        rep("null_set", length(pct_cis_sd_null_set$pct_cis_sd))))
ggplot(qcdf, aes(x = type, y = sd)) + geom_boxplot() # null should have higher sd

# Are genes with higher pct cis more likely to be divergent in multiple environments?
# the dumbest way I can think to check this is to see if the number of measurements for each gene correlates with pct_cis (or abs(pct_cis - 1)?)
plotdf <- expdf %>% group_by(c(gene_name)) %>% summarise(avg_pct_cis = mean(pct_cis), n_exp_sig = length(unique(experiment)))
plotdf %>% ggplot(aes(x = n_exp_sig, y = avg_pct_cis)) + geom_boxplot(aes(group = n_exp_sig)) + 
  geom_hline(yintercept = 1, color = "gold") + theme_classic() + ylim(c(-2, 2.5)) +
  xlab("Number of experiments gene \nis divergent in") +
  ylab("Average % cis")
plotdf %>% ggplot(aes(x = n_exp_sig, y = abs(1 - avg_pct_cis))) + geom_boxplot(aes(group = n_exp_sig)) + 
  geom_hline(yintercept = 0, color = "gold") + theme_classic() + ylim(c(0, 2.5)) +
  xlab("Number of experiments gene \nis divergent in") +
  ylab("|1 - Average % cis|")

# %cis vs block size

# load networks
load(file = "data_files/Networks_mergemodules25.RData")

expdf$module_color_cer <- sapply(expdf$gene_name, \(g) {
  g_idx <- which(colnames(counts_top3000$cer) == g)
  color <- NA
  if (length(g_idx) == 1) {
    color <- colors$cer[g_idx] %>% unlist() %>% as.character()
  }
  return(color)
}) %>% as.character()
expdf$module_color_par <- sapply(expdf$gene_name, \(g) {
  g_idx <- which(colnames(counts_top3000$par) == g)
  color <- NA
  if (length(g_idx) == 1) {
    color <- colors$par[g_idx] %>% unlist() %>% as.character()
  }
  return(color)
}) %>% as.character()
expdf$block_name <- paste(expdf$module_color_cer, expdf$module_color_par, sep = "_")
expdf$block_size <- sapply(expdf$block_name, \(b) {
  if (b == "NA_NA") {
    return(NA)
  }
  return(sum(expdf$block_name == b))
})
plotdf <- expdf
plotdf$random_mod_color <- map2(plotdf$module_color_cer, plotdf$module_color_par, \(x, y) {
  return(sample(c(x, y), 1))
}) %>% unlist() 
plotdf$block_size_rank <- sapply(plotdf$block_size, \(x) {
  return(rank(unique(plotdf$block_size), ties.method = "random")[which(unique(plotdf$block_size) == x)])
}) %>% as.numeric()
ggplot(plotdf, aes(x = block_size_rank, y = pct_cis)) + geom_point(color = plotdf$random_mod_color) + theme_classic()

# how about off-on genes? What module blocks are they in?
expdf %>% filter(gene_name %in% offongenes_parentsandhybrids) %>% select(block_name) %>% table() # visually does seem to be in smaller blocks
plotdf <- expdf %>% filter(block_name != "NA_NA") %>% 
  mutate(isoffoninparents = gene_name %in% c(offongenes_parentsandhybrids, offongenes_parentsonly),
                           islargeblock = block_size > 22) %>% 
  group_by(block_name) %>% summarise(n_offon = sum(isoffoninparents),
                                     block_size = unique(block_size))
ggplot(plotdf, aes(x = block_size, y = n_offon)) + geom_point()
