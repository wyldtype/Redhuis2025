sapply(c("tidyverse","DESeq2"), require, character.only=TRUE)

# load barkai data as DESeq2 datasets
load("Barkai_DESeq2_Data_Sets.RData")

# following DESeq2 manual: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# differential expression across experiments
library("BiocParallel")
register(MulticoreParam(4))
dds <- DESeq(dds, parallel = TRUE) # takes a hot sec
res <- results(dds)

# differential expression per experiment
# CellCycle
library("BiocParallel")
register(MulticoreParam(4))
dds_CC <- DESeq(dds_CC, parallel = TRUE)
res_CC <- results(dds_CC)
sum(res_CC$padj < 0.05, na.rm=TRUE) # so significant!

# troubleshooting low p-values
# library("IHW")
# resIHW_CC <- results(dds_CC, filterFun = ihw)
# sum(resIHW_CC$padj < 0.05, na.rm=TRUE) # yeahhhh this didn't help

# visualizing significant and non-significant genes
# plotCounts(dds_CC, gene=which.min(res_CC$padj), intgroup = "allele") # well yeah I see why that one's significant
# sig_CC <- as.data.frame(res_CC) %>% filter(padj < 10e-8, na.rm = TRUE)
# nonsig_CC <- as.data.frame(res_CC) %>% filter(padj >= 10e-8, na.rm = TRUE)
# plotCounts(dds_CC, gene=which(res_CC$padj == sample(sig_CC$padj, 1)), intgroup = "allele")
# plotCounts(dds_CC, gene=which(res_CC$padj == sample(nonsig_CC$padj, 1)), intgroup = "allele")

# conclusion: they look good as long as our p-value threshold is 10e-8, not 0.05

# LowN
dds_LowN <- DESeq(dds_LowN, parallel = TRUE)
res_LowN <- results(dds_LowN)

# LowPi
dds_LowPi <- DESeq(dds_LowPi, parallel = TRUE)
res_LowPi <- results(dds_LowPi)

# HAP4
dds_HAP4 <- DESeq(dds_HAP4, parallel = TRUE)
res_HAP4 <- results(dds_HAP4)

# finding sets of genes that are significant in each experiment and significant overall
# then checking for enrichment for duplicates
marchant <- read.table("marchant_paralogs.txt", header = TRUE)
marchant_duplicates_df <- marchant[marchant$Duplication != "S",]
marchant_duplicates <- c(marchant_duplicates_df$P1, marchant_duplicates_df$P2) %>% unique()

diverged_overall <- rownames(res)[res$padj < 10e-8]
diverged_CC <- rownames(res_CC)[res_CC$padj < 10e-8]
diverged_LowN <- rownames(res_LowN)[res_LowN$padj < 10e-8]
diverged_LowPi <- rownames(res_LowPi)[res_LowPi$padj < 10e-8]
diverged_HAP4 <- rownames(res_HAP4)[res_HAP4$padj < 10e-8]

# weird genes (ones that diverged in overall comparison but not in any of the experiment spceific ones)
weirdos <- diverged_overall[!(diverged_overall %in% diverged_CC) &
                              !(diverged_overall %in% diverged_LowN) &
                              !(diverged_overall %in% diverged_LowPi) &
                              !(diverged_overall %in% diverged_HAP4)]

diverged_CConly <- diverged_CC[!(diverged_CC %in% diverged_overall) &
                                 !(diverged_CC %in% diverged_LowN) &
                                 !(diverged_CC %in% diverged_LowPi) &
                                 !(diverged_CC %in% diverged_HAP4)]

diverged_LowNonly <- diverged_LowN[!(diverged_LowN %in% diverged_overall) &
                                 !(diverged_LowN %in% diverged_CC) &
                                 !(diverged_LowN %in% diverged_LowPi) &
                                 !(diverged_LowN %in% diverged_HAP4)]

diverged_LowPionly <- diverged_LowPi[!(diverged_LowPi %in% diverged_overall) &
                                     !(diverged_LowPi %in% diverged_CC) &
                                     !(diverged_LowPi %in% diverged_LowN) &
                                     !(diverged_LowPi %in% diverged_HAP4)]

diverged_HAP4only <- diverged_HAP4[!(diverged_HAP4 %in% diverged_overall) &
                                       !(diverged_HAP4 %in% diverged_CC) &
                                       !(diverged_HAP4 %in% diverged_LowN) &
                                       !(diverged_HAP4 %in% diverged_LowPi)]

diverged_in_all4 <- diverged_CC[diverged_CC %in% diverged_LowN &
                                 diverged_CC %in% diverged_LowPi &
                                 diverged_CC %in% diverged_HAP4]

diverged_in_only_1_exp <- c(diverged_CConly, diverged_LowNonly, diverged_LowPionly, diverged_HAP4only)
diverged_in_only_1_exp <- unique(diverged_in_only_1_exp) # this should do nothing

# fisher's exact for just genes that diverged in one experiment
n_diverging_dupes <- sum(diverged_in_only_1_exp %in% marchant_duplicates)
diverging_dupes_enrichment_counts <- matrix(nrow = 2, c(n_diverging_dupes, length(diverged_in_only_1_exp)-n_diverging_dupes, length(marchant_duplicates), nrow(dds)-length(marchant_duplicates)))
fisher.test(diverging_dupes_enrichment_counts)
cell_11 <- diverging_dupes_enrichment_counts[1,1]
cell_21 <- diverging_dupes_enrichment_counts[2,1]
cell_22 <- diverging_dupes_enrichment_counts[2,2]
cell_12 <- diverging_dupes_enrichment_counts[1,2]

# plotting for env-specific divergence
exact_test_plotting_data <- tibble("fraction_of_genes" = c(cell_11/(cell_11 + cell_21), cell_12/(cell_12 + cell_22)),
                                   "group" = c("...in genes that diverged\n in only one experiment", "...in whole genome"))
ggplot(data = exact_test_plotting_data, aes(x = group, y = fraction_of_genes, fill = group)) + geom_bar(stat = "identity") +
  xlab("") + ylab("Percent of genes that are duplicates")

# fisher's exact for genes that diverged in all 4 experiments
n_diverging_dupes <- sum(diverged_in_all4 %in% marchant_duplicates)
diverging_dupes_enrichment_counts <- matrix(nrow = 2, c(n_diverging_dupes, length(diverged_in_all4)-n_diverging_dupes, length(marchant_duplicates), nrow(dds)-length(marchant_duplicates)))
fisher.test(diverging_dupes_enrichment_counts)
cell_11 <- diverging_dupes_enrichment_counts[1,1]
cell_21 <- diverging_dupes_enrichment_counts[2,1]
cell_22 <- diverging_dupes_enrichment_counts[2,2]
cell_12 <- diverging_dupes_enrichment_counts[1,2]

# plotting for non-env-specific divergence
exact_test_plotting_data <- tibble("fraction_of_genes" = c(cell_11/(cell_11 + cell_21), cell_12/(cell_12 + cell_22)),
                                   "group" = c("...in genes that diverged\n in all 4 experiments", "...in whole genome"))
ggplot(data = exact_test_plotting_data, aes(x = group, y = fraction_of_genes, fill = group)) + geom_bar(stat = "identity") +
  xlab("") + ylab("Percent of genes that are duplicates")
