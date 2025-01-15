sapply(c("tidyverse","DESeq2"), require, character.only=TRUE)

# load barkai data as DESeq2 results tables (TODO: finish that DEcide_cis_trans.R script first)
load("Barkai_DESeq2_Results.RData") # this is a list called res_list containing 10 DESeq results (5 parental, 5 hybrid)

# following DESeq2 manual: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# following this time course example from Love et al.'s RNAseq workflow: https://doi.org/10.12688/f1000research.7035.2

# list of duplicate and singleton genes
marchant <- read.table("marchant_paralogs.txt", header = TRUE)
duplicates <- c(marchant$P1[marchant$Duplication != "S"], marchant$P2[marchant$Duplication != "S"]) %>% unique()
singletons <- c(marchant$P1[marchant$Duplication == "S"], marchant$P2[marchant$Duplication == "S"]) %>% unique()

# TODO: Use pvalues in results tables to ID cis (DE in parent and hybrid results) trans (DE in parent only) gene set and look for duplicate enrichment in each set. 
# Ultimate goal is 2 line plots (cis set and trans set) where X is number of experiments DE in (1-4) and Y is %DE, with asterisks for statistical significance in Fisher Exact tests

# Function to perform Fisher's Exact test for duplicate gene enrichment
# input: DESeqResults object r, where rownames are genes and there is a column that gives pvalues for the DE test,
# boolean mask DE_mask of length ngenes (nrows of r) saying whether or not each gene meets our DE criteria
# output: a row of a tibble that contains the pvalue, odds ratio, and percent duplicates
duplicate_enrichment_exact_test <- function(r, mask) {
  margin_d <- sum(row.names(r) %in% duplicates, na.rm = TRUE)
  margin_s <- ngenes - margin_d
  margin_DE <- nrow(r[mask,])
  margin_nonDE <- ngenes - margin_DE
  d_DE <- sum((row.names(r) %in% duplicates) & mask, na.rm = TRUE)
  s_DE <- sum((row.names(r) %in% singletons) & mask, na.rm = TRUE)
  d_nonDE <- sum((row.names(r) %in% duplicates) & mask, na.rm = TRUE)
  s_nonDE <- sum((row.names(r) %in% singletons) & mask, na.rm = TRUE)
  contingency_table <- matrix(c(d_DE, s_DE, d_nonDE, s_nonDE), nrow = 2)
  test_output <- fisher.test(contingency_table)
  return(tibble("pvalue" = test_output$p.value,
                "odds_ratio" = test_output$estimate,
                "percent_duplicates" = sum(row.names(geneset) %in% duplicates)/nrow(geneset)))
}

# Fisher's Exact Test unbinned - DE in all 4 experiments
DE_mask <- res$pvalue < p_threshold
DE_mask[is.na(DE_mask)] <- FALSE
duplicate_enrichment_exact_test(res, DE_mask)

# DE in Cell Cycle
DE_mask_CC <- res_CC$pvalue < p_threshold
DE_mask_CC[is.na(DE_mask_CC)] <- FALSE
duplicate_enrichment_exact_test(res_CC, DE_mask_CC)

# DE in Low Nitrogen
DE_mask_LowN <- res_LowN$pvalue < p_threshold
DE_mask_LowN[is.na(DE_mask_LowN)] <- FALSE
duplicate_enrichment_exact_test(res_LowN, DE_mask_LowN)

# DE in Low Pi
DE_mask_LowPi <- res_LowPi$pvalue < p_threshold
DE_mask_LowPi[is.na(DE_mask_LowPi)] <- FALSE
duplicate_enrichment_exact_test(res_LowPi, DE_mask_LowPi)

# DE in Growth Curve
DE_mask_OD <- res_OD$pvalue < p_threshold
DE_mask_OD[is.na(DE_mask_OD)] <- FALSE
duplicate_enrichment_exact_test(res_OD, DE_mask_OD)

# TODO: exact tests for genes DE in 1, 2, or 3 experiments


############## Appendix: Exact tests binned by expression level ###########

# binning into 10, 30, or 50 bins (equal partition of expression level)
avg_log_expr_vec <- log(apply(counts(dds), 1, mean, na.rm=TRUE))
avg_log_expr_df <- tibble("gene" = row.names(dds), "avg_log_expr" = avg_log_expr_vec)
avg_log_expr_df$is_duplicate <- avg_log_expr_df$gene %in% duplicates
avg_log_expr_df$is_singleton <- avg_log_expr_df$gene %in% singletons
bin_10 <- cut(avg_log_expr_vec, breaks = 10, labels = FALSE)
bin_30 <- cut(avg_log_expr_vec, breaks = 30, labels = FALSE)
bin_50 <- cut(avg_log_expr_vec, breaks = 50, labels = FALSE)

# binning into 10, 30, 50 bins (equal partition of genes per bin)
bin_10 <- ntile(avg_log_expr_vec, 10)
bin_30 <- ntile(avg_log_expr_vec, 30)
bin_50 <- ntile(avg_log_expr_vec, 50)

# plotting to make sure the genes are well distributed between bins (Note: dplyr's ntile would make them perfectly equal, if this seems like a better way to go)
par(mfrow = c(1,3))
hist(bin_10, breaks = 10)
hist(bin_30, breaks = 30)
hist(bin_50, breaks = 50)
dev.off()

# Fisher's Exact Tests for duplicates vs singletons for each bin
fishsets <- list(tibble("pvalue_DE" = 1, "pvalue_DEnotAll" = 1, "OR_DE" = 1, "OR_DEnotAll" = 1, "percent_duplicates" = 0, .rows = 10),
                 tibble("pvalue_DE" = 1, "pvalue_DEnotAll" = 1, "OR_DE" = 1, "OR_DEnotAll" = 1, "percent_duplicates" = 0, .rows = 30),
                 tibble("pvalue_DE" = 1, "pvalue_DEnotAll" = 1, "OR_DE" = 1, "OR_DEnotAll" = 1, "percent_duplicates" = 0, .rows = 50))
binsets <- list(bin_10, bin_30, bin_50)
for (i in c(1,2,3)) {
  for (j in as.numeric(rownames(table(binsets[[i]])))) {
    # retaining only known duplicates and singeltons in the current bin of interest
    mask <- (binsets[[i]] == j) & ((row.names(res) %in% singletons) |
                                     (row.names(res) %in% duplicates))
    geneset <- res[mask,]
    geneset_CC <- res_CC[mask,]
    geneset_LowN <- res_LowN[mask,]
    geneset_HAP4 <- res_HAP4[mask,]
    geneset_LowPi <- res_LowPi[mask,]
    ngenes <- nrow(geneset)
    percent_duplicates <- sum(row.names(geneset) %in% duplicates)/nrow(geneset)
    
    # marginal totals for Fisher's exact (i.e. a+c, b+d, a+b, c+d in 2x2 contingency table)
    margin_d <- sum(row.names(geneset) %in% duplicates, na.rm = TRUE)
    margin_s <- ngenes - margin_d
    margin_DE <- sum(geneset$pvalue < p_threshold, na.rm = TRUE)
    margin_nonDE <- ngenes - margin_DE
    # sum of genes with significant pvalues for at least one but not all of the 4 experiments:
    notAll_mask <- (geneset_CC$pvalue < p_threshold | geneset_LowN$pvalue < p_threshold | geneset_LowPi$pvalue < p_threshold | geneset_HAP4$pvalue < p_threshold) &
      !(geneset_CC$pvalue < p_threshold & geneset_LowN$pvalue < p_threshold & geneset_LowPi$pvalue < p_threshold & geneset_HAP4$pvalue < p_threshold)
    margin_DEnotAll <- sum(notAll_mask, na.rm = TRUE)
    margin_nonDEnotAll <- ngenes - margin_DEnotAll
    
    # box totals for Fisher's exact (i.e. a, b, c, d in 2x2 contingency table)
    # DE overall
    d_DE <- sum((row.names(geneset) %in% duplicates) & (geneset$pvalue < p_threshold), na.rm = TRUE)
    s_DE <- sum((row.names(geneset) %in% singletons) & (geneset$pvalue < p_threshold), na.rm = TRUE)
    d_nonDE <- sum((row.names(geneset) %in% duplicates) & (geneset$pvalue >= p_threshold), na.rm = TRUE)
    s_nonDE <- sum((row.names(geneset) %in% singletons) & (geneset$pvalue >= p_threshold), na.rm = TRUE)
    contingency_table_DE <- matrix(c(d_DE, s_DE, d_nonDE, s_nonDE), nrow = 2)
    # DE not all
    d_DEnotAll <- sum((row.names(geneset) %in% duplicates) & notAll_mask, na.rm = TRUE)
    s_DEnotAll <- margin_DEnotAll - d_DEnotAll
    d_nonDEnotAll <- sum((row.names(geneset) %in% duplicates) & !notAll_mask, na.rm = TRUE)
    s_nonDEnotAll <- margin_nonDEnotAll - d_nonDEnotAll
    contingency_table_DEnotAll <- matrix(c(d_DEnotAll, s_DEnotAll, d_nonDEnotAll, s_nonDEnotAll), nrow = 2)
    
    test_output_DE <- fisher.test(contingency_table_DE)
    test_output_DEnotAll <- fisher.test(contingency_table_DEnotAll)
    
    # The way I have it round, odds ratio is > 1 if duplicates are enriched for DEs and < 1 if singletons are enriched for DEs
    fishsets[[i]][j, "pvalue_DE"] <- test_output_DE$p.value
    fishsets[[i]][j, "pvalue_DEnotAll"] <- test_output_DEnotAll$p.value
    fishsets[[i]][j, "lod_DE"] <- log(test_output_DE$estimate)
    fishsets[[i]][j, "lod_DEnotAll"] <- log(test_output_DEnotAll$estimate)
    fishsets[[i]][j, "percent_duplicates"] <- percent_duplicates
  }
}

# Plotting
# DE in all 4 experiments
p_10_DE <- ggplot(data = fishsets[[1]], aes(x = c(1:10), y = lod_DE)) + 
  geom_point(aes(color = pvalue_DE)) + 
  geom_line(aes(x = c(1:10), y = percent_duplicates), color = "orange") +
  geom_hline(yintercept = 0, color = "red", linetype = "dotted") +
  xlab("Bin Number") + scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(name = "Log(Odds Ratio) \n% Duplicate/% Singleton",
                     limits = c(-1.5, 1.5),
                     sec.axis = sec_axis(trans = ~.*1, name = "Percent Duplicates")) +
  labs(color = "p-value") +
  theme(axis.title.y.right = element_text(colour = "orange"))
p_30_DE <- ggplot(data = fishsets[[2]], aes(x = c(1:30), y = lod_DE)) + 
  geom_point(aes(color = pvalue_DE)) + 
  geom_line(aes(x = c(1:30), y = percent_duplicates), color = "orange") +
  geom_hline(yintercept = 0, color = "red", linetype = "dotted") +
  xlab("Bin Number") + scale_x_continuous(breaks = seq(from = 0, to = 30, by = 5)) +
  scale_y_continuous(name = "Log(Odds Ratio) \n% Duplicate/% Singleton",
                     limits = c(-1.5, 1.5),
                     sec.axis = sec_axis(trans = ~.*1, name = "Percent Duplicates")) +
  labs(color = "p-value") +
  theme(axis.title.y.right = element_text(colour = "orange"))
p_50_DE <- ggplot(data = fishsets[[3]], aes(x = c(1:50), y = lod_DE)) + 
  geom_point(aes(color = pvalue_DE)) + 
  geom_line(aes(x = c(1:50), y = percent_duplicates), color = "orange") +
  geom_hline(yintercept = 0, color = "red", linetype = "dotted") +
  xlab("Bin Number") + scale_x_continuous(breaks = seq(from = 0, to = 50, by = 5)) +
  scale_y_continuous(name = "Log(Odds Ratio) \n% Duplicate/% Singleton",
                     limits = c(-1.5, 1.5),
                     sec.axis = sec_axis(trans = ~.*1, name = "Percent Duplicates")) +
  labs(color = "p-value") +
  theme(axis.title.y.right = element_text(colour = "orange"))

# DE in at least one but not all 4 experiments
p_10_DEnotAll <- ggplot(data = fishsets[[1]], aes(x = c(1:10), y = lod_DEnotAll)) + 
  geom_point(aes(color = pvalue_DEnotAll)) + 
  geom_line(aes(x = c(1:10), y = percent_duplicates), color = "orange") +
  geom_hline(yintercept = 0, color = "red", linetype = "dotted") +
  xlab("Bin Number") + scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(name = "Log(Odds Ratio) \n% Duplicate/% Singleton",
                     limits = c(-1.5, 1.5),
                     sec.axis = sec_axis(trans = ~.*1, name = "Percent Duplicates")) +
  labs(color = "p-value") +
  theme(axis.title.y.right = element_text(colour = "orange"))
p_30_DEnotAll <- ggplot(data = fishsets[[2]], aes(x = c(1:30), y = lod_DEnotAll)) + 
  geom_point(aes(color = pvalue_DEnotAll)) + 
  geom_line(aes(x = c(1:30), y = percent_duplicates), color = "orange") +
  geom_hline(yintercept = 0, color = "red", linetype = "dotted") +
  xlab("Bin Number") + scale_x_continuous(breaks = seq(from = 0, to = 30, by = 5)) +
  scale_y_continuous(name = "Log(Odds Ratio) \n% Duplicate/% Singleton",
                     limits = c(-1.5, 1.5),
                     sec.axis = sec_axis(trans = ~.*1, name = "Percent Duplicates")) +
  labs(color = "p-value") +
  theme(axis.title.y.right = element_text(colour = "orange"))
p_50_DEnotAll <- ggplot(data = fishsets[[3]], aes(x = c(1:50), y = lod_DEnotAll)) + 
  geom_point(aes(color = pvalue_DEnotAll)) + 
  geom_line(aes(x = c(1:50), y = percent_duplicates), color = "orange") +
  geom_hline(yintercept = 0, color = "red", linetype = "dotted") +
  xlab("Bin Number") + scale_x_continuous(breaks = seq(from = 0, to = 50, by = 5)) +
  scale_y_continuous(name = "Log(Odds Ratio) \n% Duplicate/% Singleton",
                     limits = c(-1.5, 1.5),
                     sec.axis = sec_axis(trans = ~.*1, name = "Percent Duplicates")) +
  labs(color = "p-value") +
  theme(axis.title.y.right = element_text(colour = "orange"))

library(gridExtra)
library(cowplot)
leg <- get_legend(p_10_DE)

p <- list(p_10_DE + theme(legend.position="none"), p_30_DE + theme(legend.position="none"), p_50_DE + theme(legend.position="none"), leg,
          p_10_DEnotAll + theme(legend.position="none"), p_30_DEnotAll + theme(legend.position="none"), p_50_DEnotAll + theme(legend.position="none"))
grid.arrange(arrangeGrob(p[[1]], top="10 Bins", left="DE in all 4 experiments"),arrangeGrob(p[[2]],top="30 Bins"),arrangeGrob(p[[3]],top="50 Bins"),p[[4]], arrangeGrob(p[[5]], left="DE in at least one\nbut not all 4 experiments"),p[[6]],p[[7]], nrow=2)

