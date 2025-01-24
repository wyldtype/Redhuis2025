sapply(c("tidyr", "dplyr", "purrr", "ggplot2", "ggpubr", "DESeq2"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Redhuis2025/")
load("data_files/Cleaned_TFdel_Unnormalized_Counts.RData")
load("data_files/FinalDataframe3Disp.RData")
eff_thresh <- 1
p_thresh <- 0.05

#### Data organizing ####
counts <- list(cerTP1 = counts_TFdel$cer[,infos_TFdel$cer$time_point_num == 0],
               parTP1 = counts_TFdel$par[,infos_TFdel$par$time_point_num == 0],
               hycTP1 = counts_TFdel_allele$cer[,infos_TFdel_allele$cer$time_point_num == 0],
               hypTP1 = counts_TFdel_allele$par[,infos_TFdel_allele$par$time_point_num == 0],
               cerTP2 = counts_TFdel$cer[,infos_TFdel$cer$time_point_num == 60],
               parTP2 = counts_TFdel$par[,infos_TFdel$par$time_point_num == 60],
               hycTP2 = counts_TFdel_allele$cer[,infos_TFdel_allele$cer$time_point_num == 60],
               hypTP2 = counts_TFdel_allele$par[,infos_TFdel_allele$par$time_point_num == 60],
               cerTP3 = counts_TFdel$cer[,infos_TFdel$cer$time_point_num == 960],
               parTP3 = counts_TFdel$par[,infos_TFdel$par$time_point_num == 960],
               hycTP3 = counts_TFdel_allele$cer[,infos_TFdel_allele$cer$time_point_num == 960],
               hypTP3 = counts_TFdel_allele$par[,infos_TFdel_allele$par$time_point_num == 960])
infos <- list(cerTP1 = infos_TFdel$cer[infos_TFdel$cer$time_point_num == 0,],
              parTP1 = infos_TFdel$par[infos_TFdel$par$time_point_num == 0,],
              hycTP1 = infos_TFdel_allele$cer[infos_TFdel_allele$cer$time_point_num == 0,],
              hypTP1 = infos_TFdel_allele$par[infos_TFdel_allele$par$time_point_num == 0,],
              cerTP2 = infos_TFdel$cer[infos_TFdel$cer$time_point_num == 60,],
              parTP2 = infos_TFdel$par[infos_TFdel$par$time_point_num == 60,],
              hycTP2 = infos_TFdel_allele$cer[infos_TFdel_allele$cer$time_point_num == 60,],
              hypTP2 = infos_TFdel_allele$par[infos_TFdel_allele$par$time_point_num == 60,],
              cerTP3 = infos_TFdel$cer[infos_TFdel$cer$time_point_num == 960,],
              parTP3 = infos_TFdel$par[infos_TFdel$par$time_point_num == 960,],
              hycTP3 = infos_TFdel_allele$cer[infos_TFdel_allele$cer$time_point_num == 960,],
              hypTP3 = infos_TFdel_allele$par[infos_TFdel_allele$par$time_point_num == 960,])

### removing samples missing replicates
# DESeq2 still tries to estimate a log2 fold change when a genotype has only one replicate
# # Example: YAP1 in hybrids (any timepoint)
# # in parents it has 2 replicates:
# infos$cerTP1 |> select(genotype, time_point_str) |> table()
# # in hybrids, one:
# infos$hycTP1 |> select(genotype, time_point_str) |> table()
# test <- DESeq(DESeqDataSetFromMatrix(countData = counts$hycTP1,
#                                      colData = infos$hycTP1,
#                                      design = ~ genotype))
# # while most padj are NA, not all are:
# results(test, contrast = c("genotype", "YAP1delete", "WT"),
#         alpha = 0.05)[,"padj"] |> table(useNA = "always")

# looping through counts/info to
# remove genotypes that don't have 
# 2 replicates at each timepoint
filterGenotypes <- function(.info, .counts) {
  tab <- .info |> select(genotype, time_point_str) |> table()
  no_rep_genotypes <- rownames(tab)[which(tab < 2)]
  return(filter(.info, !(genotype %in% no_rep_genotypes)))
}
infos <- map(infos, filterGenotypes)
# also updating cols of counts to match columns
updateCountsCols <- function(.counts, .info) {
  return(.counts[,.info$sample_name])
}
counts <- map2(counts, infos, updateCountsCols)

# now YAP1 isn't in hybrid
infos$cerTP1 |> select(genotype, time_point_str) |> table()
infos$hycTP1 |> select(genotype, time_point_str) |> table()

# converting genotype to a factor, after removing genotypes,
# so different number of levels per dataset
infos <- map(infos, mutate, 
             genotype = genotype |> 
               as.factor() |> 
               relevel(ref = "WT"))

# pairing counts and sample info
dds <- map2(counts, infos, \(x, y) {
  output <- DESeqDataSetFromMatrix(countData = x,
                                   colData = y,
                                   design = ~ genotype)
  return(output)
})

#### Model fitting ####
library("BiocParallel")
register(MulticoreParam(4))
dds <- map(dds, DESeq, parallel = TRUE)

# checking positive control: TDH3 decreases in GCR2delete
# TP1
test <- results(dds$cerTP1, contrast = c("genotype", "GCR2delete", "WT"),
                alpha = 0.05)
test["YGR192C",] # down in cer
test <- results(dds$parTP1, contrast = c("genotype", "GCR2delete", "WT"),
                alpha = 0.05)
test["YGR192C",] # down in par
test <- results(dds$hycTP1, contrast = c("genotype", "GCR2delete", "WT"),
                alpha = 0.05)
test["YGR192C",] # down in hyc
test <- results(dds$hypTP1, contrast = c("genotype", "GCR2delete", "WT"),
                alpha = 0.05)
test["YGR192C",] # down in hyp
# TP2
test <- results(dds$cerTP2, contrast = c("genotype", "GCR2delete", "WT"),
                alpha = 0.05)
test["YGR192C",] # down in cer
test <- results(dds$parTP2, contrast = c("genotype", "GCR2delete", "WT"),
                alpha = 0.05)
test["YGR192C",] # down in par
test <- results(dds$hycTP2, contrast = c("genotype", "GCR2delete", "WT"),
                alpha = 0.05)
test["YGR192C",]
test <- results(dds$hypTP2, contrast = c("genotype", "GCR2delete", "WT"),
                alpha = 0.05)
test["YGR192C",] # down in both hybrid alleles
# TP3
test <- results(dds$cerTP3, contrast = c("genotype", "GCR2delete", "WT"),
                alpha = 0.05)
test["YGR192C",] # down in cer
test <- results(dds$parTP3, contrast = c("genotype", "GCR2delete", "WT"),
                alpha = 0.05)
test["YGR192C",] # not really down in par
test <- results(dds$hycTP3, contrast = c("genotype", "GCR2delete", "WT"),
                alpha = 0.05)
test["YGR192C",]
test <- results(dds$hypTP3, contrast = c("genotype", "GCR2delete", "WT"),
                alpha = 0.05)
test["YGR192C",] # down in both hybrid alleles

#### Formatting data ####
# extracting LFCs, SEs, and pvalues from each TFdel-WT comparison in each 
# species and output as TFdeldfs
griddfTP1 <- bind_rows(infos$cerTP1,
                    infos$parTP1,
                    infos$hycTP1,
                    infos$hypTP1) |> 
  select(genotype, allele, organism) |> 
  filter(genotype != "WT") |> 
  mutate(timepoint = "TP1",
         orgallele = if_else(organism == "hyb",
                             false = organism,
                             true = if_else(allele == "cer",
                                            true = "hyc",
                                            false = "hyp"))) |> 
  unique()
griddfTP2 <- bind_rows(infos$cerTP2,
                       infos$parTP2,
                       infos$hycTP2,
                       infos$hypTP2) |> 
  select(genotype, allele, organism) |> 
  filter(genotype != "WT") |> 
  mutate(timepoint = "TP2",
         orgallele = if_else(organism == "hyb",
                             false = organism,
                             true = if_else(allele == "cer",
                                            true = "hyc",
                                            false = "hyp"))) |> 
  unique()
griddfTP3 <- bind_rows(infos$cerTP3,
                       infos$parTP3,
                       infos$hycTP3,
                       infos$hypTP3) |> 
  select(genotype, allele, organism) |> 
  filter(genotype != "WT") |> 
  mutate(timepoint = "TP3",
         orgallele = if_else(organism == "hyb",
                             false = organism,
                             true = if_else(allele == "cer",
                                            true = "hyc",
                                            false = "hyp"))) |> 
  unique()
griddf <- bind_rows(griddfTP1, griddfTP2, griddfTP3)

TFdeldf <- map(c(1:nrow(griddf)), \(i) {
  del <- griddf$genotype[i] |> as.character()
  org <- griddf$orgallele[i]
  tp <- griddf$timepoint[i]
  cat(i, del, org, tp, "\n")
  res <- results(dds[[paste0(org, tp)]], contrast = c("genotype", del, "WT"),
                 alpha = 0.05, cooksCutoff = FALSE, independentFiltering = FALSE)
  return(tibble(gene_name = rownames(counts_TFdel$cer),
                deletion = gsub("delete", "", del),
                organism = org,
                timepoint = tp,
                basemean = res$baseMean,
                lfc = res$log2FoldChange,
                lfcSE = res$lfcSE,
                pval = res$pvalue,
                padj = res$padj))
}) |> bind_rows()

# A few gene/TF/tp combos didn't converge in DESeq2 or were missing replicates:
length(unique(TFdeldf$deletion))
sum(is.na(TFdeldf))
sum(is.na(TFdeldf$padj))
table(is.na(TFdeldf$padj), TFdeldf$timepoint)
sum(is.na(TFdeldf$lfc)) + sum(is.na(TFdeldf$lfcSE)) + sum(is.na(TFdeldf$pval)) + sum(is.na(TFdeldf$padj))

TFdeldf <- drop_na(TFdeldf)
length(unique(TFdeldf$deletion)) # we shouldn't lose any TFs entirely

sum(TFdeldf$padj < p_thresh)
sum(TFdeldf$padj >= p_thresh) # should have way fewer significant effects

#### Negative control: WT samples 28 vs 2 replicates, "leave one out" as sham TFdels ####
set.seed(23)

# looping through counts, change counts for each TFdel to be 2 randomly sampled WT samples (or however many replicates that TF has)
counts_sham1 <- map2(counts, infos, \(x, y) {
  for (del in setdiff(unique(y$genotype), "WT")) {
    samps <- filter(y, genotype == del) |> select(sample_name) |> pull()
    wt_samps <- filter(y, genotype == "WT") |> select(sample_name) |> pull()
    x[, samps] <- x[, sample(wt_samps, length(samps), replace = FALSE)] # because we loop through this sampling for every deletion it is effectively sampling with replacement. We just want to make sure the same sample isn't used twice for the same genotype, as variation for every gene would be 0
  }
  return(x)
})

set.seed(42)

# same thing twice, just randomizing which WT samples go to which TF
counts_sham2 <- map2(counts, infos, \(x, y) {
  for (del in setdiff(unique(y$genotype), "WT")) {
    samps <- filter(y, genotype == del) |> select(sample_name) |> pull()
    wt_samps <- filter(y, genotype == "WT") |> select(sample_name) |> pull()
    x[, samps] <- x[, sample(wt_samps, length(samps), replace = FALSE)] # because we loop through this sampling for every deletion it is effectively sampling with replacement. We just want to make sure the same sample isn't used twice for the same genotype, as variation for every gene would be 0
  }
  return(x)
})

# pairing counts and sample info
dds_sham1 <- map2(counts_sham1, infos, \(x, y) {
  output <- DESeqDataSetFromMatrix(countData = x,
                                   colData = y,
                                   design = ~ genotype)
  return(output)
})

dds_sham2 <- map2(counts_sham2, infos, \(x, y) {
  output <- DESeqDataSetFromMatrix(countData = x,
                                   colData = y,
                                   design = ~ genotype)
  return(output)
})

### Model fitting
library("BiocParallel")
register(MulticoreParam(4))
dds_sham1 <- map(dds_sham1, DESeq, parallel = TRUE)
dds_sham2 <- map(dds_sham2, DESeq, parallel = TRUE)

# checking control: TDH3 decreases in GCR2delete
# TP1
test <- results(dds_sham1$cerTP1, contrast = c("genotype", "GCR2delete", "WT"),
                alpha = 0.05)
test["YGR192C",] # shouldn't be DE anymore, WT vs WT

### Formatting data
# extracting LFCs, SEs, and pvalues from each TFdel-WT comparison in each 
# species and output as TFdeldfs
TFdeldf_sham1 <- map(c(1:nrow(griddf)), \(i) {
  del <- griddf$genotype[i] |> as.character()
  org <- griddf$orgallele[i]
  tp <- griddf$timepoint[i]
  cat(i, del, org, tp, "\n")
  res <- results(dds_sham1[[paste0(org, tp)]], contrast = c("genotype", del, "WT"),
                 alpha = 0.05, cooksCutoff = FALSE, independentFiltering = FALSE)
  return(tibble(gene_name = rownames(counts_TFdel$cer),
                deletion = gsub("delete", "", del),
                organism = org,
                timepoint = tp,
                basemean = res$baseMean,
                lfc = res$log2FoldChange,
                lfcSE = res$lfcSE,
                pval = res$pvalue,
                padj = res$padj))
}) |> bind_rows()

TFdeldf_sham2 <- map(c(1:nrow(griddf)), \(i) {
  del <- griddf$genotype[i] |> as.character()
  org <- griddf$orgallele[i]
  tp <- griddf$timepoint[i]
  cat(i, del, org, tp, "\n")
  res <- results(dds_sham2[[paste0(org, tp)]], contrast = c("genotype", del, "WT"),
                 alpha = 0.05, cooksCutoff = FALSE, independentFiltering = FALSE)
  return(tibble(gene_name = rownames(counts_TFdel$cer),
                deletion = gsub("delete", "", del),
                organism = org,
                timepoint = tp,
                basemean = res$baseMean,
                lfc = res$log2FoldChange,
                lfcSE = res$lfcSE,
                pval = res$pvalue,
                padj = res$padj))
}) |> bind_rows()

# A few gene/TF/tp combos didn't converge in DESeq2 or were missing replicates:
TFdeldf_sham1 <- drop_na(TFdeldf_sham1)
TFdeldf_sham2 <- drop_na(TFdeldf_sham2)

sum(TFdeldf_sham1$padj < p_thresh)
sum(TFdeldf_sham1$padj >= p_thresh) # should have way fewer significant effects
sum(TFdeldf_sham2$padj < p_thresh)
sum(TFdeldf_sham2$padj >= p_thresh)

#### Saving ####
save(TFdeldf, TFdeldf_sham1, TFdeldf_sham2, file = "data_files/TFdel_DESeq2.RData")

############## Probably Archive - Data Exploration moved to Fig script ########################## 

#### How many genes up or down in each deletion? ####
# for genes in each divergence category:
# 1-1, 2-2, 1-2, 2-1
# Note some genes are missing from finaldf b/c they are lowly
# expressed in LowN:
sum(!unique(TFdeldf$gene_name) %in% finaldf$gene_name[finaldf$experiment == "LowN"])
TFdeldf <- finaldf |> filter(experiment == "LowN") |> 
  select(gene_name, cer, par) |> 
  left_join(TFdeldf, by = "gene_name")
# change to different clusters to explore
plotdf <- TFdeldf |> filter(cer == 2 & par == 1 & timepoint == "TP3") |> 
  mutate(lfc_dir = if_else(pval < p_thresh & abs(lfc) > eff_thresh,
                           true = sign(lfc),
                           false = NA)) |> 
  drop_na() |> 
  group_by(deletion, organism) |> 
  summarise(n_up = sum(lfc_dir == 1),
            n_down = sum(lfc_dir == -1))
ggplot(filter(plotdf, organism %in% c("cer", "par")), 
       aes(x = n_up, y = n_down)) +
  geom_point(aes(color = deletion, shape = organism)) +
  theme(legend.position = "none") +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic()
# TODO: it doesn't super matter how many genes are up or down,
# That's a statement about how powerful each TF is at regulating stuff
# (which we know is variable), plus if there are more
# genes in one category (i.e. 700 in 1-1 vs 50 in 2-1), there will
# invariably be more genes affected
# More important to see how all the ups are pushing towards
# the species that is up and each TF only affects one allele
# or the other
table(TFdeldf$cer, TFdeldf$par)
# Instead three points to convey, which might need two separate plots:
# 1) deletions tend to push one species towards the other's expression,
# rather than away from
# 2) deletions tend to only affect one species or the other
# 3) divergent groups have a higher proportion of stronger effects than conserved

# Exploring example
# HAP1 del in 1-1 at TP3 interesting example b/c clearly has cer
# average lower that WT and par not changed
plotdf <- TFdeldf |> filter(cer == 1 & par == 1 & 
                              timepoint == "TP3" &
                              deletion == "HAP1" &
                              organism %in% c("cer", "par")) |> 
  select(gene_name, deletion, lfc, pval, organism) |> 
  mutate(sig = pval < p_thresh & abs(lfc) > 1) |> 
  pivot_wider(id_cols = c("gene_name", "deletion"),
              names_from = "organism", values_from = c("lfc", "pval", "sig"))

ggplot(plotdf, aes(x = lfc_cer, y = lfc_par)) +
  geom_point(aes(color = paste(sig_cer, sig_par))) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlim(c(-5, 5)) +
  ylim(c(-5, 5))

# MSN2, same deal
plotdf <- TFdeldf |> filter(cer == 1 & par == 1 & 
                              timepoint == "TP3" &
                              deletion == "MSN2" &
                              organism %in% c("cer", "par")) |> 
  select(gene_name, deletion, lfc, pval, organism) |> 
  mutate(sig = pval < p_thresh & abs(lfc) > 2) |> 
  pivot_wider(id_cols = c("gene_name", "deletion"),
              names_from = "organism", values_from = c("lfc", "pval", "sig"))

ggplot(plotdf, aes(x = lfc_cer, y = lfc_par)) +
  geom_point(aes(color = paste(sig_cer, sig_par))) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(slope = 1, yintercept = 0) +
  xlim(c(-5, 5)) +
  ylim(c(-5, 5)) # But this is much more believably having the
# net effect of shifting cer down

#### Exploration: LFC WT vs TFdel ####
# Selecting highly expressed, significant results
sig_gene_dels <- filter(TFdeldf, mean > 100 & abs(lfc) > 1 & pval < 1e-5) |> 
  select(gene_name, deletion) |> 
  unique()
TFdeldf_sig <- map2(sig_gene_dels$gene_name, sig_gene_dels$deletion, \(.g, .del) {
  TFdeldf |> filter(gene_name == .g & deletion == .del)
}) |> bind_rows()

# cer LFC(WT, TFdel) vs par LFC(WT, TFdel), separate plots for parents and hybrids
# and cer/par LFC(WT, TFdel) vs hyc/hyp LFC(WT, TFdel)
plotlist_hyb <- vector(mode = "list", length = length(unique(griddf$genotype)))
plotlist_parents <- vector(mode = "list", length = length(unique(griddf$genotype)))
for (i in c(1:length(unique(TFdeldf_sig$deletion)))) {
  del <- unique(TFdeldf_sig$deletion)[i]
  # hybrid
  plotdf <- left_join(filter(TFdeldf_sig, deletion == del & organism == "hyc"),
                      filter(TFdeldf_sig, deletion == del & organism == "hyp"),
                      by = c("gene_name", "deletion"),
                      suffix = c("_cer", "_par"))
  max_lfc <- c(plotdf$lfc_cer, plotdf$lfc_par) |> abs() |> max()
  p <- ggplot(data = plotdf, aes(x = lfc_cer, y = lfc_par)) + 
    geom_point(alpha = 0.5) +
    geom_text(aes(label = gene_name), check_overlap = FALSE) +
    ggtitle(del) +
    geom_abline(intercept = 0, slope = 1, color = "gold") +
    xlab(paste("log2(mean WT expr) - log2(mean", del, "delete expr)\n Hybrid Scer allele")) +
    ylab(paste("log2(mean WT expr) - log2(mean", del, "delete expr)\n Hybrid Spar allele")) +
    theme_classic() +
    xlim(c(-max_lfc - 1, max_lfc + 1)) +
    ylim(c(-max_lfc - 1, max_lfc + 1))
  plotlist_hyb[[i]] <- p
  # parents
  plotdf <- left_join(filter(TFdeldf_sig, deletion == del & organism == "cer"),
                      filter(TFdeldf_sig, deletion == del & organism == "par"),
                      by = c("gene_name", "deletion"),
                      suffix = c("_cer", "_par"))
  max_lfc <- c(plotdf$lfc_cer, plotdf$lfc_par) |> abs() |> max()
  p <- ggplot(data = plotdf, aes(x = lfc_cer, y = lfc_par)) + 
    geom_point(alpha = 0.5) +
    geom_text(aes(label = gene_name), check_overlap = FALSE) +
    ggtitle(del) +
    geom_abline(intercept = 0, slope = 1, color = "gold") +
    xlab(paste("log2(mean WT expr) - log2(mean", del, "delete expr)\n S. cerevisiae")) +
    ylab(paste("log2(mean WT expr) - log2(mean", del, "delete expr)\n S. paradoxus")) +
    theme_classic() +
    xlim(c(-max_lfc - 1, max_lfc + 1)) +
    ylim(c(-max_lfc - 1, max_lfc + 1))
  plotlist_parents[[i]] <- p
}
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/TFdelLFCratios_hyb.pdf",
    width = 12, height = ceiling(length(plotlist_hyb)/3)*4)
ggarrange(plotlist = plotlist_hyb, 
          nrow = ceiling(length(plotlist_hyb)/3), ncol = 3)
dev.off()
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/TFdelLFCratios_parents.pdf",
    width = 12, height = ceiling(length(plotlist_parents)/3)*4)
ggarrange(plotlist = plotlist_parents, 
          nrow = ceiling(length(plotlist_parents)/3), ncol = 3)
dev.off()

#### Visualizing example genes in cer par hyc hyp ####
# normalizing counts for visualization
countsPerMillion <- function(.cts) {
  librarySizes <- colSums(.cts, na.rm = TRUE)
  output <- apply(.cts, 1, function(x) {
    normalized <- (x/librarySizes)*1e6
    return(round(normalized))
  })
  return(t(output)) # For some unhinged reason, a vector output of apply across ROWS forms the COLUMNS of a new matrix
}
# tests for countsPerMillion
test_cts <- counts_TFdel$cer # change to counts_allele to change which set you're testing
test_rowIdx <- sample(c(1:nrow(test_cts)), 1)
test_colIdx <- sample(c(1:ncol(test_cts)), 1)
test_count <- test_cts[test_rowIdx, test_colIdx]
test_cpm <- countsPerMillion(test_cts)
test_cpm_cpm <- countsPerMillion(test_cpm)
((test_count/colSums(test_cts, na.rm = TRUE)[test_colIdx])*1e6) %>% round() # what it should be
test_cts[test_rowIdx, test_colIdx] # what it is before using our function
test_cpm[test_rowIdx, test_colIdx] # what it is using our function
test_cpm_cpm[test_rowIdx, test_colIdx] # what it is if you run cpm too many times (should be the same as test_cpm)

counts_norm <- map(counts_TFdel, countsPerMillion)
counts_norm_allele <- map(counts_TFdel_allele, countsPerMillion)

compareQuartet <- function(.g, .del) {
  plotdf <- bind_rows(tibble(expr = counts_norm$cer[.g, infos_TFdel$cer$genotype %in% c(paste0(.del, "delete"), "WT")],
                             timepoint = infos_TFdel$cer$time_point_str[infos_TFdel$cer$genotype %in% c(paste0(.del, "delete"), "WT")],
                             genotype = infos_TFdel$cer$genotype[infos_TFdel$cer$genotype %in% c(paste0(.del, "delete"), "WT")],
                             type = "cer"),
                      tibble(expr = counts_norm$par[.g, infos_TFdel$par$genotype %in% c(paste0(.del, "delete"), "WT")],
                             timepoint = infos_TFdel$par$time_point_str[infos_TFdel$par$genotype %in% c(paste0(.del, "delete"), "WT")],
                             genotype = infos_TFdel$par$genotype[infos_TFdel$par$genotype %in% c(paste0(.del, "delete"), "WT")],
                             type = "par"),
                      tibble(expr = counts_norm_allele$cer[.g, infos_TFdel_allele$cer$genotype %in% c(paste0(.del, "delete"), "WT")],
                             timepoint = infos_TFdel_allele$cer$time_point_str[infos_TFdel_allele$cer$genotype %in% c(paste0(.del, "delete"), "WT")],
                             genotype = infos_TFdel_allele$cer$genotype[infos_TFdel_allele$cer$genotype %in% c(paste0(.del, "delete"), "WT")],
                             type = "hyc"),
                      tibble(expr = counts_norm_allele$par[.g, infos_TFdel_allele$par$genotype %in% c(paste0(.del, "delete"), "WT")],
                             timepoint = infos_TFdel_allele$par$time_point_str[infos_TFdel_allele$par$genotype %in% c(paste0(.del, "delete"), "WT")],
                             genotype = infos_TFdel_allele$par$genotype[infos_TFdel_allele$par$genotype %in% c(paste0(.del, "delete"), "WT")],
                             type = "hyp"))
  p <- ggplot(plotdf, aes(x = type, y = log2(expr + 1))) +
    geom_violin(data = filter(plotdf, genotype == "WT"), color = "grey") + 
    geom_jitter(data = filter(plotdf, genotype == paste0(.del, "delete")), aes(color = type)) +
    facet_wrap(~ timepoint) +
    theme_classic() +
    xlab("") +
    ylab("expression (log2)") +
    ggtitle(paste(.g, .del, sep = ", ")) +
    theme(legend.title = element_blank())
  return(p)
}
compareQuartet(.g = "YGR192C", .del = "GCR2")

# TODO: Use compare quartet to look at more examples

# less-specific to single TFs
compareQuartet(.g = "YPL201C", .del = "TEC1") # maybe hybrid cer allele is unusually high, but this just looks like WT cis-divergence with cer higher
compareQuartet(.g = "YDL200C", .del = "ARG81") # same deal as above, but maybe hyc is a little low
compareQuartet(.g = "YDL200C", .del = "AFT1")
compareQuartet(.g = "YDR033W", .del = "MBP1") # par has lowe rexpression, not allele-specific

compareQuartet(.g = "YPR194C", .del = "GLN3") # yes.
compareQuartet(.g = "YPR194C", .del = "SOK2") # yes?
compareQuartet(.g = "YPR194C", .del = "MIG1") # no?
compareQuartet(.g = "YPR194C", .del = "MBP1") # no.

# more specific
compareQuartet(.g = "YDR439W", .del = "AFT1") # deletion equalizes expression in hybrid alleles (by bringing hyc down at 1st 2 timepoints)
compareQuartet(.g = "YMR058W", .del = "AFT1") # deletion shuts this gene off, but par and hyp recover better
compareQuartet(.g = "YLR194C", .del = "INO4") # well that's funky. Deletion causes this gene to increase overtime, but especially in par and hyp
compareQuartet(.g = "YFL051C", .del = "SOK2") # barely expressed in cer/hyc, highest in par and hyp at TP1, deletion attenuates TP1 expression
compareQuartet(.g = "YAL034W-A", .del = "ROX1") # expressed higher in par/hyp, deletion equalizes things
compareQuartet(.g = "YAL034W-A", .del = "MSN2") # same deal here

#### Regress lfc diff by baseline allele difference ####

plotlist_hyb <- vector(mode = "list", length = length(unique(griddf$genotype)))
plotlist_parents <- vector(mode = "list", length = length(unique(griddf$genotype)))
for (i in c(1:length(unique(TFdeldf_sig$deletion)))) {
  del <- unique(TFdeldf_sig$deletion)[i]
  # hybrid
  plotdf <- left_join(filter(TFdeldf, deletion == del & organism == "hyc"),
                      filter(TFdeldf, deletion == del & organism == "hyp"),
                      by = c("gene_name", "deletion"),
                      suffix = c("_cer", "_par")) |> 
    mutate(lfc_diff = lfc_cer - lfc_par)
  plotdf <- left_join(plotdf, tibble(gene_name = rownames(TFdeldf_hybrid),
                                     lfc_allele = TFdeldf_hybrid$log2FoldChange),
                      by = "gene_name", relationship = "many-to-one")
  max_lfc <- c(plotdf$lfc_diff, plotdf$lfc_allele) |> abs() |> max()
  p <- ggplot(data = plotdf, aes(x = lfc_allele, y = lfc_diff)) + 
    geom_point(alpha = 0.5) +
    geom_text(aes(label = gene_name), check_overlap = TRUE) +
    ggtitle(del) +
    geom_abline(intercept = 0, slope = 1, color = "gold") +
    xlab(paste("Overall")) +
    ylab(paste("TFdel")) +
    theme_classic() +
    xlim(c(-max_lfc - 1, max_lfc + 1)) +
    ylim(c(-max_lfc - 1, max_lfc + 1))
  plotlist_hyb[[i]] <- p
  # parents
  plotdf <- left_join(filter(TFdeldf, deletion == del & organism == "cer"),
                      filter(TFdeldf, deletion == del & organism == "par"),
                      by = c("gene_name", "deletion"),
                      suffix = c("_cer", "_par")) |> 
    mutate(lfc_diff = lfc_cer - lfc_par)
  plotdf <- left_join(plotdf, tibble(gene_name = rownames(TFdeldf_parents),
                                     lfc_allele = TFdeldf_parents$log2FoldChange),
                      by = "gene_name", relationship = "many-to-one")
  max_lfc <- c(plotdf$lfc_diff, plotdf$lfc_allele) |> abs() |> max()
  p <- ggplot(data = plotdf, aes(x = lfc_allele, y = lfc_diff)) + 
    geom_point(alpha = 0.5) +
    geom_text(aes(label = gene_name), check_overlap = TRUE) +
    ggtitle(del) +
    geom_abline(intercept = 0, slope = 1, color = "gold") +
    xlab(paste("Overall")) +
    ylab(paste("TFdel")) +
    theme_classic() +
    xlim(c(-max_lfc - 1, max_lfc + 1)) +
    ylim(c(-max_lfc - 1, max_lfc + 1))
  plotlist_parents[[i]] <- p
}
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/TFdelLFCdiffs_hyb.pdf",
    width = 12, height = ceiling(length(plotlist_hyb)/3)*4)
ggarrange(plotlist = plotlist_hyb, 
          nrow = ceiling(length(plotlist_hyb)/3), ncol = 3)
dev.off()
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/TFdelLFCdiffs_parents.pdf",
    width = 12, height = ceiling(length(plotlist_parents)/3)*4)
ggarrange(plotlist = plotlist_parents, 
          nrow = ceiling(length(plotlist_parents)/3), ncol = 3)
dev.off()


                      
