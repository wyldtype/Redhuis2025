sapply(c("tidyr", "dplyr", "purrr", "ggplot2", "bipartite", "ggpubr"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Redhuis2024/")
source("functions_for_figure_scripts.R")
load("data_files/FinalDataframe.RData")
load("data_files/TFdel.RData")

TFnames <- setdiff(unique(gsub("delete", "", sample_info_tfdel$genotype)), "WT")
#### Comparing cluster expression in different TF deletions ####
checkTF <- function(.tf, .group, .parents_or_hybrid = "parents",
                    .show_wt = TRUE, .show_tfdel = TRUE, 
                    .plotlims = NULL, .normalization = "log2") {
  gene_idxs <- finaldf |> filter(experiment == "LowN" &
                                   group == .group) |> 
    select(gene_name) |> pull()
  if (.parents_or_hybrid == "parents") {
    cer_wt_cols <- which(sample_info_tfdel$genotype == "WT" & sample_info_tfdel$allele == "cer")
    par_wt_cols <- which(sample_info_tfdel$genotype == "WT" & sample_info_tfdel$allele == "par")
    cer_tfdel_cols <- which(sample_info_tfdel$genotype == paste0(.tf, "delete") & sample_info_tfdel$allele == "cer")
    par_tfdel_cols <- which(sample_info_tfdel$genotype == paste0(.tf, "delete") & sample_info_tfdel$allele == "par")
    p <- plotExpressionProfileTFdel(.cts1 = counts_tfdel[gene_idxs, cer_wt_cols],
                                    .cts2 = counts_tfdel[gene_idxs, par_wt_cols],
                                    .cts3 = counts_tfdel[gene_idxs, cer_tfdel_cols],
                                    .cts4 = counts_tfdel[gene_idxs, par_tfdel_cols],
                                    .info1 = sample_info_tfdel[cer_wt_cols,],
                                    .info2 = sample_info_tfdel[par_wt_cols,],
                                    .info3 = sample_info_tfdel[cer_tfdel_cols,],
                                    .info4 = sample_info_tfdel[par_tfdel_cols,],
                                    .normalization = .normalization,
                                    .show_points_wt = TRUE,
                                    .show_points_tfdel = TRUE,
                                    .show_wt = .show_wt,
                                    .show_tfdel = .show_tfdel,
                                    .show_lines_tfdel = TRUE,
                                    .plotlims = .plotlims)
  }
  if (.parents_or_hybrid == "hybrid") {
    cer_wt_cols <- which(sample_info_tfdel_allele$genotype == "WT" & sample_info_tfdel_allele$allele == "cer")
    par_wt_cols <- which(sample_info_tfdel_allele$genotype == "WT" & sample_info_tfdel_allele$allele == "par")
    cer_tfdel_cols <- which(sample_info_tfdel_allele$genotype == paste0(.tf, "delete") & sample_info_tfdel_allele$allele == "cer")
    par_tfdel_cols <- which(sample_info_tfdel_allele$genotype == paste0(.tf, "delete") & sample_info_tfdel_allele$allele == "par")
    p <- plotExpressionProfileTFdel(.cts1 = counts_tfdel_allele[gene_idxs, cer_wt_cols],
                                    .cts2 = counts_tfdel_allele[gene_idxs, par_wt_cols],
                                    .cts3 = counts_tfdel_allele[gene_idxs, cer_tfdel_cols],
                                    .cts4 = counts_tfdel_allele[gene_idxs, par_tfdel_cols],
                                    .info1 = sample_info_tfdel_allele[cer_wt_cols,],
                                    .info2 = sample_info_tfdel_allele[par_wt_cols,],
                                    .info3 = sample_info_tfdel_allele[cer_tfdel_cols,],
                                    .info4 = sample_info_tfdel[par_tfdel_cols,],
                                    .normalization = .normalization,
                                    .show_points_wt = TRUE,
                                    .show_points_tfdel = TRUE,
                                    .show_wt = .show_wt,
                                    .show_tfdel = .show_tfdel,
                                    .show_lines_tfdel = TRUE,
                                    .plotlims = .plotlims)
  }
  return(p)
}
# tests for checkTF
random_tf <- sample(TFnames, size = 1)
annotate_figure(checkTF(.tf = random_tf, .group = "dyn21"), top = random_tf)
annotate_figure(checkTF(.tf = random_tf, .group = "dyn21", .parents_or_hybrid = "hybrid"), top = random_tf)

#### plotting TFdels vs WT in each divergence group ####
filter(finaldf, experiment == "LowN") |> select(group) |> table() |> sort(decreasing = TRUE)
### dynamics divergers 2-1
# what are the 2-1 genes doing in other environments?
groupdf <- finaldf |> pivot_longer(cols = c("cer", "par"),
                                   names_to = "allele",
                                   values_to = "cluster") |> 
  pivot_wider(id_cols = c("gene_name", "allele"),
              values_from = "cluster",
              names_from = c("experiment"))
groupdf |> filter((allele == "cer" & LowN == 2) |
                    (allele == "par" & LowN == 1)) |> 
  group_by(allele, CC, HAP4, LowPi, Heat, Cold) |> 
  summarise(count = n()) |> arrange(desc(count))
# collecting plots
plotlist <- vector(mode = "list", length = length(TFnames)*2)
names(plotlist) <- c(paste(TFnames, "parents", sep = "_"),
                     paste(TFnames, "hybrid", sep = "_")) |> sort(decreasing = TRUE)
for (tf in TFnames) {
  plotlist[[paste(tf, "parents", sep = "_")]] <- annotate_figure(checkTF(.tf = tf, .group = "dyn21", .parents_or_hybrid = "parents"), top = tf)
  plotlist[[paste(tf, "hybrid", sep = "_")]] <- annotate_figure(checkTF(.tf = tf, .group = "dyn21", .parents_or_hybrid = "hybrid"), top = tf)
}
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/TFdel/dyn21.pdf",
    width = 9, height = ceiling(length(plotlist)/2)*4)
ggarrange(plotlist = plotlist, ncol = 2, nrow = ceiling(length(plotlist)/2),
          common.legend = TRUE)
dev.off()
# dyn 1-2
plotlist <- vector(mode = "list", length = length(TFnames)*2)
names(plotlist) <- c(paste(TFnames, "parents", sep = "_"),
                     paste(TFnames, "hybrid", sep = "_")) |> sort(decreasing = TRUE)
for (tf in TFnames) {
  plotlist[[paste(tf, "parents", sep = "_")]] <- annotate_figure(checkTF(.tf = tf, .group = "dyn12", .parents_or_hybrid = "parents"), top = tf)
  plotlist[[paste(tf, "hybrid", sep = "_")]] <- annotate_figure(checkTF(.tf = tf, .group = "dyn12", .parents_or_hybrid = "hybrid"), top = tf)
}
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/TFdel/dyn12.pdf",
    width = 9, height = ceiling(length(plotlist)/2)*4)
ggarrange(plotlist = plotlist, ncol = 2, nrow = ceiling(length(plotlist)/2),
          common.legend = TRUE)
dev.off()
### level plots
# levuppar1
plotlist <- vector(mode = "list", length = length(TFnames)*2)
names(plotlist) <- c(paste(TFnames, "parents", sep = "_"),
                     paste(TFnames, "hybrid", sep = "_")) |> sort(decreasing = TRUE)
for (tf in TFnames) {
  plotlist[[paste(tf, "parents", sep = "_")]] <- annotate_figure(checkTF(.tf = tf, .group = "levuppar1", .parents_or_hybrid = "parents"), top = tf)
  plotlist[[paste(tf, "hybrid", sep = "_")]] <- annotate_figure(checkTF(.tf = tf, .group = "levuppar1", .parents_or_hybrid = "hybrid"), top = tf)
}
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/TFdel/levuppar1.pdf",
    width = 9, height = ceiling(length(plotlist)/2)*4)
ggarrange(plotlist = plotlist, ncol = 2, nrow = ceiling(length(plotlist)/2),
          common.legend = TRUE)
dev.off()
# levuppar2
plotlist <- vector(mode = "list", length = length(TFnames)*2)
names(plotlist) <- c(paste(TFnames, "parents", sep = "_"),
                     paste(TFnames, "hybrid", sep = "_")) |> sort(decreasing = TRUE)
for (tf in TFnames) {
  plotlist[[paste(tf, "parents", sep = "_")]] <- annotate_figure(checkTF(.tf = tf, .group = "levuppar2", .parents_or_hybrid = "parents"), top = tf)
  plotlist[[paste(tf, "hybrid", sep = "_")]] <- annotate_figure(checkTF(.tf = tf, .group = "levuppar2", .parents_or_hybrid = "hybrid"), top = tf)
}
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/TFdel/levuppar2.pdf",
    width = 9, height = ceiling(length(plotlist)/2)*4)
ggarrange(plotlist = plotlist, ncol = 2, nrow = ceiling(length(plotlist)/2),
          common.legend = TRUE)
dev.off()
# levupcer1
plotlist <- vector(mode = "list", length = length(TFnames)*2)
names(plotlist) <- c(paste(TFnames, "parents", sep = "_"),
                     paste(TFnames, "hybrid", sep = "_")) |> sort(decreasing = TRUE)
for (tf in TFnames) {
  plotlist[[paste(tf, "parents", sep = "_")]] <- annotate_figure(checkTF(.tf = tf, .group = "levupcer1", .parents_or_hybrid = "parents"), top = tf)
  plotlist[[paste(tf, "hybrid", sep = "_")]] <- annotate_figure(checkTF(.tf = tf, .group = "levupcer1", .parents_or_hybrid = "hybrid"), top = tf)
}
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/TFdel/levupcer1.pdf",
    width = 9, height = ceiling(length(plotlist)/2)*4)
ggarrange(plotlist = plotlist, ncol = 2, nrow = ceiling(length(plotlist)/2),
          common.legend = TRUE)
dev.off()
# levupcer2
plotlist <- vector(mode = "list", length = length(TFnames)*2)
names(plotlist) <- c(paste(TFnames, "parents", sep = "_"),
                     paste(TFnames, "hybrid", sep = "_")) |> sort(decreasing = TRUE)
for (tf in TFnames) {
  plotlist[[paste(tf, "parents", sep = "_")]] <- annotate_figure(checkTF(.tf = tf, .group = "levupcer2", .parents_or_hybrid = "parents"), top = tf)
  plotlist[[paste(tf, "hybrid", sep = "_")]] <- annotate_figure(checkTF(.tf = tf, .group = "levupcer2", .parents_or_hybrid = "hybrid"), top = tf)
}
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/TFdel/levupcer2.pdf",
    width = 9, height = ceiling(length(plotlist)/2)*4)
ggarrange(plotlist = plotlist, ncol = 2, nrow = ceiling(length(plotlist)/2),
          common.legend = TRUE)
dev.off()

#### GLN3 deletion ####
# From previous plots, GLN3delete was identified to cause genes
# diverging in both 2-1 and 1-2 to 
# decrease over time in parents and hybrid (appear like 2 cluster)
# in both alleles
# BUT the same change to GLN3 binding/expression *can't* be 
# causing both the divergence in 2-1 and 1-2 genes b/c the deletion
# has the same effect on both sets of genes
# 1-2
checkTF("GLN3", "dyn12", .show_wt = TRUE, .show_tfdel = FALSE, .plotlims = c(7.5, 8.5))
checkTF("GLN3", "dyn12", .show_wt = TRUE, .show_tfdel = TRUE, .plotlims = c(7.5, 8.5))
checkTF("GLN3", "dyn12", .show_wt = TRUE, .show_tfdel = FALSE, .parents_or_hybrid = "hybrid", .plotlims = c(7.5, 8.5))
checkTF("GLN3", "dyn12", .show_wt = TRUE, .show_tfdel = TRUE, .parents_or_hybrid = "hybrid", .plotlims = c(7.5, 8.5))

# 2-1
checkTF("GLN3", "dyn21", .show_wt = TRUE, .show_tfdel = FALSE, .plotlims = c(7.5, 8.5))
checkTF("GLN3", "dyn21", .show_wt = TRUE, .show_tfdel = TRUE, .plotlims = c(7.5, 8.5))
checkTF("GLN3", "dyn21", .show_wt = TRUE, .show_tfdel = FALSE, .parents_or_hybrid = "hybrid", .plotlims = c(7.5, 8.5))
checkTF("GLN3", "dyn21", .show_wt = TRUE, .show_tfdel = TRUE, .parents_or_hybrid = "hybrid", .plotlims = c(7.5, 8.5))

#### Early/late cer/par-ification poster children ####
checkTF("NRG1", "dyn12", .show_wt = TRUE, .show_tfdel = FALSE, .plotlims = c(7.5, 8.4)) # WT
checkTF("NRG1", "dyn12", .show_wt = TRUE, .show_tfdel = TRUE, .plotlims = c(7.5, 8.4)) # no effect
checkTF("URE2", "dyn12", .show_wt = TRUE, .show_tfdel = TRUE, .plotlims = c(7.5, 8.4)) # both timepoints
checkTF("MSN2", "dyn12", .show_wt = TRUE, .show_tfdel = TRUE, .plotlims = c(7.5, 8.4)) # late Spar-ify

checkTF("INO4", "dyn12", .show_wt = TRUE, .show_tfdel = TRUE, .plotlims = c(7.5, 8.4)) # late Scer-ify
checkTF("GLN3", "dyn12", .show_wt = TRUE, .show_tfdel = TRUE, .plotlims = c(7.5, 8.4)) # early Spar-ify
checkTF("AFT1", "dyn12", .show_wt = TRUE, .show_tfdel = TRUE, .plotlims = c(7.5, 8.4)) # early Scer-ify

# 21
checkTF("NRG1", "dyn21", .show_wt = TRUE, .show_tfdel = FALSE, .plotlims = c(7.5, 8.4)) # WT
checkTF("NRG1", "dyn21", .show_wt = TRUE, .show_tfdel = TRUE, .plotlims = c(7.5, 8.4)) # no effect
checkTF("URE2", "dyn21", .show_wt = TRUE, .show_tfdel = TRUE, .plotlims = c(7.5, 8.4)) # both timepoints
checkTF("MSN2", "dyn21", .show_wt = TRUE, .show_tfdel = TRUE, .plotlims = c(7.5, 8.4)) # late Spar-ify

checkTF("DAL80", "dyn21", .show_wt = TRUE, .show_tfdel = TRUE, .plotlims = c(7.5, 8.4))
checkTF("AFT1", "dyn21", .show_wt = TRUE, .show_tfdel = TRUE, .plotlims = c(7.5, 8.4)) # early Scer-ify
checkTF("RTG1", "dyn21", .show_wt = TRUE, .show_tfdel = TRUE, .plotlims = c(7.5, 8.4)) # early Scer-ify

#### Additional filter for lowly expressed genes ####
# even though we already did this to the full dataset, now we need to filter 
# additional genes that aren't expressed in the LowN environment
# for each gene, plot mean counts (lcpm) across all samples versus
# var counts (lcpm) across all samples (4 plots: cer, par, hyc, hyp)
plotdf <- bind_rows(tibble(mean_expr = colMeans(log2(counts_TFdel$cer + 1), na.rm = TRUE),
                           sd_expr = apply(log2(counts_TFdel$cer), 2, sd, na.rm = TRUE),
                           type = "cer"),
                    tibble(mean_expr = colMeans(log2(counts_TFdel$par + 1), na.rm = TRUE),
                           sd_expr = apply(log2(counts_TFdel$par + 1), 2, sd, na.rm = TRUE),
                           type = "par"),
                    tibble(mean_expr = colMeans(log2(counts_TFdel_allele$cer + 1), na.rm = TRUE),
                           sd_expr = apply(log2(counts_TFdel_allele$cer + 1), 2, sd, na.rm = TRUE),
                           type = "hyc"),
                    tibble(mean_expr = colMeans(log2(counts_TFdel_allele$par + 1), na.rm = TRUE),
                           sd_expr = apply(log2(counts_TFdel_allele$par + 1), 2, sd, na.rm = TRUE),
                           type = "hyp"))
cutoffExpr <- 5
ggplot(plotdf, aes(x = mean_expr, y = sd_expr)) + 
  geom_hex() +
  facet_wrap(~type) +
  geom_vline(xintercept = cutoffExpr, color = "red") 

# keeping genes that are above the expression threshold in any of the 4 counts sets
keep <- (colMeans(log2(counts_TFdel$cer + 1), na.rm = TRUE) > cutoffExpr) |
  (colMeans(log2(counts_TFdel$par + 1), na.rm = TRUE) > cutoffExpr) |
  (colMeans(log2(counts_TFdel_allele$cer + 1), na.rm = TRUE) > cutoffExpr) |
  (colMeans(log2(counts_TFdel_allele$par + 1), na.rm = TRUE) > cutoffExpr)
sum(keep)
sum(!keep)
keep_genes <- colnames(counts_TFdel$cer)[keep]
# filtering
TFdeldf_cer <- filter(TFdeldf_cer, gene_name %in% keep_genes)
TFdeldf_par <- filter(TFdeldf_par, gene_name %in% keep_genes)
TFdeldf_hyc <- filter(TFdeldf_hyc, gene_name %in% keep_genes)
TFdeldf_hyp <- filter(TFdeldf_hyp, gene_name %in% keep_genes)

#### Individual TF WT vs TFdel LFC(hyc vs hyp) scatterplots ####
# motivation: genes off y=x line have the most divergent cis responses between hyc and hyp
# (just hybrid for now, to keep things as controlled as possible)
# LFC will be positive when m2 > m1

# TODO: exploring taking mean across timepoints for more consistent effects
# change means and sds back to vectors if it obscures timepoint-specific effects too much

calculateLogFoldChange <- function(.m1, .m2, .sd1, .sd2, .threshold = 100) {
  .m1 <- mean(.m1)
  .m2 <- mean(.m2)
  .sd1 <- mean(.sd1)
  .sd2 <- mean(.sd2)
  isLowlyExpressed <- .m1 < .threshold & .m2 < .threshold
  m1_less <- ((.m1 < .m2 - .sd2) & (.m1 < .m2 + .sd2))
  m1_greater <- ((.m1 > .m2 - .sd2) & (.m1 > .m2 + .sd2))
  m2_less <- ((.m2 < .m1 - .sd1) & (.m2 < .m1 + .sd1))
  m2_greater <- ((.m2 > .m1 - .sd1) & (.m2 > .m1 + .sd1))
  isOutsideCI <- (m1_less | m1_greater) & (m2_less | m2_greater)
  lfc <- log2(.m2 + 1) - log2(.m1 + 1)
  #lfc[!isOutsideCI] <- NA
  lfc[isLowlyExpressed] <- NA
  return(lfc)
}
# TODO: incorporate uncertainty for replicates with a lot of variation
# for now omitted b/c it leads to too many adverse effects 
# (example: hypothetical gene which straddles the threshold for ommission)
# tests for calculateLogFoldChange
calculateLogFoldChange(1, 2, 1, 1) # too noisy and lowly expressed
calculateLogFoldChange(1, 2, .1, .1) # too lowly expressed
calculateLogFoldChange(101, 2, 150, 4) # too noisy
calculateLogFoldChange(100, 2, .9, .9) # significant negative LFC
calculateLogFoldChange(.m1 = c(1, 1, 101, 100), .m2 = c(2, 2, 2, 2),
                       .sd1 = c(1, .1, 150, .9), .sd2 = c(1, .1, 4, .9)) # same examples, vectorized

# real example, YOR348C in GLN3 delete
# lfc_del (should be positive at TP3, cause it's higher in cer allele)
calculateLogFoldChange(.m1 = TFdeldf_hyp |> filter(deletion == "GLN3" & gene_name == "YOR348C") |> 
                         select(mean_del) |> pull(),
                       .m2 = TFdeldf_hyc |> filter(deletion == "GLN3" & gene_name == "YOR348C") |> 
                         select(mean_del) |> pull(),
                       .sd1 = TFdeldf_hyp |> filter(deletion == "GLN3" & gene_name == "YOR348C") |> 
                         select(sd_del) |> pull(),
                       .sd2 = TFdeldf_hyc |> filter(deletion == "GLN3" & gene_name == "YOR348C") |> 
                         select(sd_del) |> pull())
# should be even stronger in parents TP3
calculateLogFoldChange(.m1 = TFdeldf_par |> filter(deletion == "GLN3" & gene_name == "YOR348C") |> 
                         select(mean_del) |> pull(),
                       .m2 = TFdeldf_cer |> filter(deletion == "GLN3" & gene_name == "YOR348C") |> 
                         select(mean_del) |> pull(),
                       .sd1 = TFdeldf_par |> filter(deletion == "GLN3" & gene_name == "YOR348C") |> 
                         select(sd_del) |> pull(),
                       .sd2 = TFdeldf_cer |> filter(deletion == "GLN3" & gene_name == "YOR348C") |> 
                         select(sd_del) |> pull())

griddf <- expand_grid(deletion = unique(TFdeldf_hyc$deletion),
                      timepoint = unique(TFdeldf_hyc$time_point_str))
# TFdeldf_hyb <- map2(griddf$deletion, griddf$timepoint, \(.del, .tp) {
#   mean_hyc_wt <- TFdeldf_hyc |> filter(deletion == .del & time_point_str == .tp) |> 
#     select(mean_wt) |> pull()
#   mean_hyp_wt <- TFdeldf_hyp |> filter(deletion == .del & time_point_str == .tp) |> 
#     select(mean_wt) |> pull()
#   mean_hyc_del <- TFdeldf_hyc |> filter(deletion == .del & time_point_str == .tp) |> 
#     select(mean_del) |> pull()
#   mean_hyp_del <- TFdeldf_hyp |> filter(deletion == .del & time_point_str == .tp) |> 
#     select(mean_del) |> pull()
#   sd_hyc_wt <- TFdeldf_hyc |> filter(deletion == .del & time_point_str == .tp) |> 
#     select(sd_wt) |> pull()
#   sd_hyp_wt <- TFdeldf_hyp |> filter(deletion == .del & time_point_str == .tp) |> 
#     select(sd_wt) |> pull()
#   sd_hyc_del <- TFdeldf_hyc |> filter(deletion == .del & time_point_str == .tp) |> 
#     select(sd_del) |> pull()
#   sd_hyp_del <- TFdeldf_hyp |> filter(deletion == .del & time_point_str == .tp) |> 
#     select(sd_del) |> pull()
#   gene_name_hyc <- TFdeldf_hyc |> filter(deletion == .del & time_point_str == .tp) |> 
#     select(gene_name) |> pull()
#   gene_name_hyp <- TFdeldf_hyp |> filter(deletion == .del & time_point_str == .tp) |> 
#     select(gene_name) |> pull()
#   if (sum(gene_name_hyc == gene_name_hyp) != length(gene_name_hyc)) {
#     stop("gene names are not in the same order, mean expr ratios will not be ratios of the same gene!")
#   }
#   # log(cer/par), wt samples
#   lfc_wt <- calculateLogFoldChange(.m1 = mean_hyp_wt, .m2 = mean_hyc_wt, .sd1 = sd_hyp_wt, .sd2 = sd_hyc_wt)
#   # log(cer/par), TFdel samples
#   lfc_del <- calculateLogFoldChange(.m1 = mean_hyp_del, .m2 = mean_hyc_del, .sd1 = sd_hyp_del, .sd2 = sd_hyc_del)
#   # log(tfdel/wt), cer samples (not filtered for low expression)
#   lfc_c <- calculateLogFoldChange(.m1 = mean_hyc_wt, .m2 = mean_hyc_del, .sd1 = sd_hyc_wt, .sd2 = sd_hyc_del, .threshold = 0)
#   # log(tfdel/wt), par samples (not filtered for low expression)
#   lfc_p <- calculateLogFoldChange(.m1 = mean_hyp_wt, .m2 = mean_hyp_del, .sd1 = sd_hyp_wt, .sd2 = sd_hyp_del, .threshold = 0)
#   return(tibble(gene_name = gene_name_hyc,
#                 deletion = .del,
#                 time_point_str = .tp,
#                 lfc_wt = lfc_wt,
#                 lfc_del = lfc_del,
#                 lfc_c = lfc_c,
#                 lfc_p = lfc_p,
#                 mean_c_wt = mean_hyc_wt,
#                 mean_p_wt = mean_hyp_wt,
#                 mean_c_del = mean_hyc_del,
#                 mean_p_del = mean_hyp_del))
# }) |> reduce(.f = bind_rows)

TFdeldf_hyb <- full_join(x = TFdeldf_hyc |> 
                           select(sd_wt, mean_wt, sd_del, mean_del, 
                                  deletion, time_point_str, gene_name) |> 
                           rename(c("sd_wt_c"="sd_wt",
                                      "sd_del_c"="sd_del",
                                      "mean_wt_c"="mean_wt",
                                      "mean_del_c"="mean_del")),
                         y = TFdeldf_hyp |> 
                           select(sd_wt, mean_wt, sd_del, mean_del, 
                                  deletion, time_point_str, gene_name) |> 
                           rename(c("sd_wt_p"="sd_wt",
                                      "sd_del_p"="sd_del",
                                      "mean_wt_p"="mean_wt",
                                      "mean_del_p"="mean_del")),
                                  by = c("gene_name", "deletion", "time_point_str")) |> 
  group_by(gene_name, deletion) |> 
  summarise(lfc_wt = calculateLogFoldChange(.m1 = mean_wt_p, .m2 = mean_wt_c, .sd1 = sd_wt_p, .sd2 = sd_wt_c),
            lfc_del = calculateLogFoldChange(.m1 = mean_del_p, .m2 = mean_del_c, .sd1 = sd_del_p, .sd2 = sd_del_c))

TFdeldf_hyb # TODO: why the NAs (it's exactly half of values, so probably something to do with an allele)


plotlist <- vector(mode = "list", length = nrow(griddf))
for (i in c(1:nrow(griddf))) {
  del <- griddf$deletion[i]
  tp <- griddf$timepoint[i]
  plotdf <- TFdeldf_hyb |> filter(deletion == del & time_point_str == tp) |> 
    select(lfc_wt, lfc_del, gene_name) |> 
    drop_na() # removing genes that are lowly expressed in both hyc and hyp of the same genotype (note this might remove some genes that are really off in a certain TFdel, but only when it's true of both alleles, so not divergent)
  max_lfc <- c(plotdf$lfc_wt, plotdf$lfc_del) |> abs() |> max()
  p <- ggplot(data = plotdf, aes(x = lfc_wt, y = lfc_del)) + 
    geom_point(alpha = 0.5) +
    geom_text(aes(label = gene_name), check_overlap = TRUE) +
    ggtitle(paste(del, tp)) +
    geom_abline(intercept = 0, slope = 1, color = "gold") +
    xlab("log2(mean Hc expr) - log2(mean Hp expr)\nWT") +
    ylab("log2(mean Hc expr) - log2(mean Hp expr)\nTF deletion") +
    theme_classic() +
    xlim(c(-max_lfc, max_lfc)) +
    ylim(c(-max_lfc, max_lfc))
  plotlist[[i]] <- p
}
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/TFdelLFCratios_hyb.pdf",
    width = 12, height = 95)
ggarrange(plotlist = plotlist, nrow = 19, ncol = 3)
dev.off()

# repeat in parents to see if same bias exists in certain TFs/timepoints
TFdeldf_parents <- map2(griddf$deletion, griddf$timepoint, \(.del, .tp) {
  mean_cer_wt <- TFdeldf_cer |> filter(deletion == .del & time_point_str == .tp) |> 
    select(mean_wt) |> pull()
  mean_par_wt <- TFdeldf_par |> filter(deletion == .del & time_point_str == .tp) |> 
    select(mean_wt) |> pull()
  mean_cer_del <- TFdeldf_cer |> filter(deletion == .del & time_point_str == .tp) |> 
    select(mean_del) |> pull()
  mean_par_del <- TFdeldf_par |> filter(deletion == .del & time_point_str == .tp) |> 
    select(mean_del) |> pull()
  sd_cer_wt <- TFdeldf_cer |> filter(deletion == .del & time_point_str == .tp) |> 
    select(sd_wt) |> pull()
  sd_par_wt <- TFdeldf_par |> filter(deletion == .del & time_point_str == .tp) |> 
    select(sd_wt) |> pull()
  sd_cer_del <- TFdeldf_cer |> filter(deletion == .del & time_point_str == .tp) |> 
    select(sd_del) |> pull()
  sd_par_del <- TFdeldf_par |> filter(deletion == .del & time_point_str == .tp) |> 
    select(sd_del) |> pull()
  gene_name_cer <- TFdeldf_cer |> filter(deletion == .del & time_point_str == .tp) |> 
    select(gene_name) |> pull()
  gene_name_par <- TFdeldf_par |> filter(deletion == .del & time_point_str == .tp) |> 
    select(gene_name) |> pull()
  if (sum(gene_name_cer == gene_name_par) != length(gene_name_cer)) {
    stop("gene names are not in the same order, mean expr ratios will not be ratios of the same gene!")
  }
  lfc_wt <- calculateLogFoldChange(.m1 = mean_par_wt, .m2 = mean_cer_wt, .sd1 = sd_par_wt, .sd2 = sd_cer_wt)
  lfc_del <- calculateLogFoldChange(.m1 = mean_par_del, .m2 = mean_cer_del, .sd1 = sd_par_del, .sd2 = sd_cer_del)
  lfc_c <- calculateLogFoldChange(.m1 = mean_cer_wt, .m2 = mean_cer_del, .sd1 = sd_cer_wt, .sd2 = sd_cer_del)
  lfc_p <- calculateLogFoldChange(.m1 = mean_par_wt, .m2 = mean_par_del, .sd1 = sd_par_wt, .sd2 = sd_par_del)
  return(tibble(gene_name = gene_name_cer,
                deletion = .del,
                time_point_str = .tp,
                lfc_wt = lfc_wt,
                lfc_del = lfc_del,
                lfc_c = lfc_c,
                lfc_p = lfc_p,
                mean_c_wt = mean_cer_wt,
                mean_p_wt = mean_par_wt,
                mean_c_del = mean_cer_del,
                mean_p_del = mean_par_del))
}) |> reduce(.f = bind_rows)

plotlist <- vector(mode = "list", length = nrow(griddf))
for (i in c(1:nrow(griddf))) {
  del <- griddf$deletion[i]
  tp <- griddf$timepoint[i]
  plotdf <- TFdeldf_parents |> filter(deletion == del & time_point_str == tp) |> 
    select(lfc_wt, lfc_del, gene_name) |> 
    drop_na()
  max_lfc <- c(plotdf$lfc_wt, plotdf$lfc_del) |> abs() |> max()
  p <- ggplot(data = plotdf, aes(x = lfc_wt, y = lfc_del)) + 
    geom_point(alpha = 0.5) +
    geom_text(aes(label = gene_name), check_overlap = TRUE) +
    ggtitle(paste(del, tp)) +
    geom_abline(intercept = 0, slope = 1, color = "gold") +
    xlab("log2(mean Pc expr) - log2(mean Pp expr)\nWT") +
    ylab("log2(mean Pc expr) - log2(mean Pp expr)\nTF deletion") +
    theme_classic() +
    xlim(c(-max_lfc, max_lfc)) +
    ylim(c(-max_lfc, max_lfc))
  plotlist[[i]] <- p
}
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/TFdelLFCratios_parents.pdf",
    width = 12, height = 95)
ggarrange(plotlist = plotlist, nrow = 19, ncol = 3)
dev.off()

#### How are LFC values distributed for parents and hybrid? ####
library(ggridges)
# hybrid
plotdf <- pivot_longer(TFdeldf_hyb, 
                       cols = c("lfc_wt", 
                                "lfc_del", 
                                "lfc_c", 
                                "lfc_p"),
                       names_to = "lfc_type", values_to = "lfc_value")
ggplot(plotdf, aes(x = lfc_value)) + 
  geom_density_ridges(aes(y = lfc_type, fill = lfc_type))
# parents
plotdf <- pivot_longer(TFdeldf_parents, 
                       cols = c("lfc_wt", 
                                "lfc_del", 
                                "lfc_c", 
                                "lfc_p"),
                       names_to = "lfc_type", values_to = "lfc_value")
ggplot(plotdf, aes(x = lfc_value)) + 
  geom_density_ridges(aes(y = lfc_type, fill = lfc_type))

# are the genes with divergent responses to TF deletion in parents shared with hybrids?
plotdf <- bind_rows(TFdeldf_hyb |> 
                      mutate(type = "hybrid",
                             lfc_diff = lfc_wt - lfc_del),
                    TFdeldf_parents |>
                      mutate(type = "parents",
                             lfc_diff = lfc_wt - lfc_del))
plotdf <- pivot_wider(plotdf, id_cols = c("gene_name", "deletion", "time_point_str"), 
                      names_from = "type", values_from = "lfc_diff",
                      names_prefix = "lfc_diff_")
griddf <- expand_grid(deletion = unique(plotdf$deletion),
                      time_point_str = unique(plotdf$time_point_str))
plotlist <- vector(mode = "list", length = nrow(griddf))
radius <- 2
circledf <- tibble(x = radius*cos(seq(from = 0, to = 2*pi, length.out = 100)),
                   y = radius*sin(seq(from = 0, to = 2*pi, length.out = 100)))
for (i in c(1:nrow(griddf))) {
  del <- griddf$deletion[i]
  tp <- griddf$time_point_str[i]
  p_df <- plotdf |> filter(deletion == del & time_point_str == tp)
  # setting genes too lowly expressed in parents or hybrid to 0
  p_df$lfc_diff_hybrid[is.na(p_df$lfc_diff_hybrid)] <- 0
  p_df$lfc_diff_parents[is.na(p_df$lfc_diff_parents)] <- 0
  max_expr <- max(abs(c(p_df$lfc_diff_hybrid, p_df$lfc_diff_parents)), na.rm = TRUE)
  p <- ggplot(p_df, aes(x = lfc_diff_parents, y = lfc_diff_hybrid)) +
    geom_point() +
    geom_text(aes(label = gene_name), check_overlap = TRUE) +
    geom_abline(slope = 1, intercept = 0, color = "gold") +
    geom_vline(xintercept = 0, color = "red") +
    geom_hline(yintercept = 0, color = "blue") +
    geom_path(data = circledf, aes(x = x, y = y), color = "forestgreen") +
    xlim(c(-max_expr, max_expr)) +
    ylim(c(-max_expr, max_expr)) +
    theme_classic() +
    ggtitle(paste(del, tp))
  plotlist[[i]] <- p
}

pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/TFdelLFCdiff.pdf",
    width = 15, height = 95)
ggarrange(plotlist = plotlist, nrow = 19, ncol = 3)
dev.off()


#### examining interesting looking examples ####
compareQuartet <- function(.g, .del) {
  max_expr <- c(counts_TFdel_allele$cer[.g, infos_TFdel_allele$cer$genotype == .del, drop = FALSE],
                counts_TFdel_allele$par[.g, infos_TFdel_allele$par$genotype == .del, drop = FALSE],
                counts_TFdel_allele$cer[.g, infos_TFdel_allele$cer$genotype == "WT", drop = FALSE],
                counts_TFdel_allele$par[.g, infos_TFdel_allele$par$genotype == "WT", drop = FALSE],
                counts_TFdel$cer[.g, infos_TFdel$cer$genotype == .del, drop = FALSE],
                counts_TFdel$par[.g, infos_TFdel$par$genotype == .del, drop = FALSE],
                counts_TFdel$cer[.g, infos_TFdel$cer$genotype == "WT", drop = FALSE],
                counts_TFdel$par[.g, infos_TFdel$par$genotype == "WT", drop = FALSE]) |> 
    unlist() |> 
    max(na.rm = TRUE)
  p_hybrid <- plotExpressionProfileQuartet(.cts1 = counts_TFdel_allele$cer[.g, infos_TFdel_allele$cer$genotype == .del, drop = FALSE],
                                           .cts2 = counts_TFdel_allele$par[.g, infos_TFdel_allele$par$genotype == .del, drop = FALSE],
                                           .cts3 = counts_TFdel_allele$cer[.g, infos_TFdel_allele$cer$genotype == "WT", drop = FALSE],
                                           .cts4 = counts_TFdel_allele$par[.g, infos_TFdel_allele$par$genotype == "WT", drop = FALSE],
                                           .info1 = filter(infos_TFdel_allele$cer, genotype == .del),
                                           .info2 = filter(infos_TFdel_allele$par, genotype == .del),
                                           .info3 = filter(infos_TFdel_allele$cer, genotype == "WT"),
                                           .info4 = filter(infos_TFdel_allele$par, genotype == "WT"),
                                           .color1 = "orange1",
                                           .color2 = "blue2",
                                           .color3 = "orange4",
                                           .color4 = "blue4",
                                           .name1 = paste("hyc", .del),
                                           .name2 = paste("hyp", .del),
                                           .name3 = "hyc WT",
                                           .name4 = "hyp WT",
                                           .normalization = "log2",
                                           .method = "line",
                                           .show_points = TRUE,
                                           .show_confidence_intervals = FALSE) 
  p_parents <- plotExpressionProfileQuartet(.cts1 = counts_TFdel$cer[.g, infos_TFdel$cer$genotype == .del, drop = FALSE],
                                            .cts2 = counts_TFdel$par[.g, infos_TFdel$par$genotype == .del, drop = FALSE],
                                            .cts3 = counts_TFdel$cer[.g, infos_TFdel$cer$genotype == "WT", drop = FALSE],
                                            .cts4 = counts_TFdel$par[.g, infos_TFdel$par$genotype == "WT", drop = FALSE],
                                            .info1 = filter(infos_TFdel$cer, genotype == .del),
                                            .info2 = filter(infos_TFdel$par, genotype == .del),
                                            .info3 = filter(infos_TFdel$cer, genotype == "WT"),
                                            .info4 = filter(infos_TFdel$par, genotype == "WT"),
                                            .color1 = "orange1",
                                            .color2 = "blue2",
                                            .color3 = "orange4",
                                            .color4 = "blue4",
                                            .name1 = paste("cer", .del),
                                            .name2 = paste("par", .del),
                                            .name3 = "cer WT",
                                            .name4 = "par WT",
                                            .normalization = "log2",
                                            .method = "line",
                                            .show_points = TRUE,
                                            .show_confidence_intervals = FALSE)
  return(annotate_figure(ggarrange(p_parents + ylim(c(0, max_expr)),
                   p_hybrid + ylim(c(0, max_expr)), nrow = 1, ncol = 2),
                   top = .g))
} # TODO: this function no longer works since I've made the TFdel data unnormalized with replicates

# DE in parents and hybrid (TFdel LFC goes to 0)
compareQuartet(.del = "SOK2delete", .g = "YPR194C") # that's real. paradoxus allele is restored upon SOK2 deletion. cerevisiae allele is unaffected

# divergent in parents, kind of divergent in hybrid
compareQuartet(.del = "GLN3delete", .g = "YOR348C")

# divergent in parents but not hybrid
compareQuartet(.del = "GLN3delete", .g = "YBR291C")
compareQuartet(.del = "DAL80delete", .g = "YLR164W")

# This is probably real but looks a lot like a divergent gene with no TFdel effect:
# (what makes me think it's real is that it has the same effect in the parent and hybrid)
# YOR329C in AFT1 TP1 parents and hybrid (TFdel LFC goes to 0)
compareQuartet(.del = "AFT1delete", .g = "YOR329C") 
# gene is expressed higher in par and hyp
# AFT1 deletion drops par/hyp (and slightly raises cer/hyc) YPD expression
# then over N starvation, the gene returns to its WT divergent expr

# YDL200C in many hybrid TFs but not parents (TFdel LFC goes to 0)
compareQuartet(.del = "AFT1delete", .g = "YDL200C")
compareQuartet(.del = "GLN3delete", .g = "YDL200C")
compareQuartet(.del = "INO4delete", .g = "YDL200C")
compareQuartet(.del = "MBP1delete", .g = "YDL200C")
compareQuartet(.del = "MIG1delete", .g = "YDL200C")
compareQuartet(.del = "PHD1delete", .g = "YDL200C")
compareQuartet(.del = "ROX1delete", .g = "YDL200C")
compareQuartet(.del = "SOK2delete", .g = "YDL200C") # these are basically all driven by a spike in TP2 Hp

# YPL201C in many parent and hybrid TFs (much higher expr in cer/hyc)
compareQuartet(.del = "AFT1delete", .g = "YPL201C") # Hc TP2 spike
compareQuartet(.del = "GLN3delete", .g = "YPL201C") # Hc TP2 spike
compareQuartet(.del = "INO4delete", .g = "YPL201C") # Hc TP2 spike, Pp also goes up TP3 though
compareQuartet(.del = "MBP1delete", .g = "YPL201C") # Hc TP2 spike
compareQuartet(.del = "MIG1delete", .g = "YPL201C") # Hc and Hp increase in deletion
compareQuartet(.del = "PHD1delete", .g = "YPL201C") # Hc overall increase
compareQuartet(.del = "ROX1delete", .g = "YPL201C") # Hc TP2 spike
compareQuartet(.del = "SOK2delete", .g = "YPL201C") # Hc TP2 spike

# YLL057C and YBR054W in GCR2 TP3 hybrids not parents (TFdel reveals divergence) (YLL057C also in MIG1 TP2 and MSN2 TP2 and NRG1 TP2 parents)
compareQuartet(.del = "GCR2delete", .g = "YBR054W")

# TEC1 YOL052C-A is higher in parents TP2, and maybe us unnamed point in hybrid TP2?
compareQuartet(.del = "TEC1delete", .g = "YOL052C-A") # seems like mostly just driven by that 0 count in Pp TP2

# looking at LFCdiff plot outliers
# positive hyb-specific divergence (misexpression)
compareQuartet(.del = "GZF3delete", .g = "YGR199W")
compareQuartet(.del = "CHA4delete", .g = "YGR199W")
compareQuartet(.del = "DAL80delete", .g = "YGR199W") # same situation in all of these: this gene down at TP1 in deletion, but not specific to one TFdel
# diverged in both parents and hybrid same direction
compareQuartet(.del = "DAL80delete", .g = "YOL058W") # TP1 is outside LFCdiff forestgreen circle, but mainly an artifact of this gene not being strongly expressed until TP2
# diverged in parents, not in hybrid
compareQuartet(.del = "GCR2delete", .g = "YNL160W") # GCR2 deletion sends cer consistently up. In the hybrid, the effect is the same as cer and not allele specific

#### Choosing DE criteria ####
# In order to decide how many effects are shared between parent and hybrid,
# we need to decide on the criteria (discrete or continuous) for DE upon TF deletion
# and DE between cer and par

# here's a random gene/deletion to get a feel for
random_gene <- sample(intersect(pull(select(drop_na(TFdeldf_parents), gene_name)), 
                                pull(select(drop_na(TFdeldf_hyb), gene_name))), 1)
random_del <- sample(unique(TFdeldf_parents$deletion), 1)
compareQuartet(.del = paste0(random_del, "delete"), .g = random_gene)
bind_rows(
  x = TFdeldf_hyb |> filter(gene_name == random_gene & deletion == random_del) |> 
  mutate(lfc_cer = log2(mean_hyc_del + 1) - log2(mean_hyc_wt + 1),
         lfc_par = log2(mean_hyp_del + 1) - log2(mean_hyp_wt + 1)) |> 
  select(time_point_str, lfc_wt, lfc_del, lfc_cer, lfc_par),
  y = TFdeldf_parents |> filter(gene_name == random_gene & deletion == random_del) |> 
  mutate(lfc_cer = log2(mean_cer_del + 1) - log2(mean_cer_wt + 1),
         lfc_par = log2(mean_par_del + 1) - log2(mean_par_wt + 1)) |> 
  select(time_point_str, lfc_wt, lfc_del, lfc_cer, lfc_par))
# random gene/del observations: 
# TFdel effect usually not strong, TFdel genotype close to WT of same species/allele
# plenty of genes have different cer/par level, usually less severe in hybrid


# here's a good example of a gene that should be considered DE and divergent,
# but the effect is stronger in parent than hybrid (still present in hybrid)
compareQuartet(.del = "GLN3delete", .g = "YOR348C")
# here's an example of one where nothing is divergent in hybrid, but expr is
# affected by deletion:
compareQuartet(.del = "GLN3delete", .g = "YBR291C")
# example of what we don't want to call divergent (but miiiight be DE?)
compareQuartet(.del = "TEC1delete", .g = "YOL052C-A") # this gene was flagged for TP2, but it's just b/c the par TP2 is lowly expressed and ends up at 0




#### Upset plots: how many DE are shared between parents and hybrid ###

# TODO: use upset plots to see how many gene/deletions/timepoints
# with LFC of sufficient magnitude are shared between parents and hybrid
# (there's probably a better way to set up the dataframe than what I'm doing
# below it's just 7pm and I'm tired)
# TODO: I'll probably need to select the genes farthest from y=x in
# the wt vs tfdel lfc(c/p) plots (red and green genes in tilt troubleshoot plots)
cutoff_lfc <- 2
DEup <- bind_rows(mutate(filter(TFdeldf_hyb, lfc_wt > cutoff_lfc),
                                 type = "wt_hyb"),
                          mutate(filter(TFdeldf_parents, lfc_wt > cutoff_lfc),
                                 type = "wt_parents"),
                          mutate(filter(TFdeldf_hyb, lfc_del > cutoff_lfc),
                                 type = "del_hyb"),
                          mutate(filter(TFdeldf_parents, lfc_del > cutoff_lfc),
                                 type = "del_parents")) |> 
  select(gene_name, deletion, time_point_str, lfc_wt, lfc_del, lfc_c, lfc_p, type)
DEdown <- bind_rows(mutate(filter(TFdeldf_hyb, lfc_wt < -cutoff_lfc),
                         type = "wt_hyb"),
                  mutate(filter(TFdeldf_parents, lfc_wt < -cutoff_lfc),
                         type = "wt_parents"),
                  mutate(filter(TFdeldf_hyb, lfc_del < -cutoff_lfc),
                         type = "del_hyb"),
                  mutate(filter(TFdeldf_parents, lfc_del < -cutoff_lfc),
                         type = "del_parents")) |> 
  select(gene_name, deletion, time_point_str, lfc_wt, lfc_del, lfc_c, lfc_p, type)

# or doing one value for hybrid and one for parent per TF/deletion/timepoint:
plotdf <- bind_rows(TFdeldf_hyb |> 
                      drop_na() |> 
                      mutate(type = "hybrid",
                             lfc_diff = lfc_wt - lfc_del),
                    TFdeldf_parents |> 
                      drop_na() |> 
                      mutate(type = "parents",
                             lfc_diff = lfc_wt - lfc_del)) |> 
  filter(abs(lfc_diff) > 0.5)
lt <- list(hyb = filter(plotdf, type == "hybrid"),
           parents = filter(plotdf, type == "parents")) |> 
  map(.f = mutate, gene_tf_tp = paste0(gene_name, deletion, time_point_str)) |> 
  map(.f = select, gene_tf_tp) |> 
  map(.f = pull)
plotdf <- make_comb_mat(lt)
p <- UpSet(plotdf, set_order = c("hybrid", "parents"),
           comb_order = order(comb_size(plotdf)), # this, I will say, makes ggplot look like a freakin' dream. I literally just wanted to rename the histogram y axis and add counts to the top of the bars
           top_annotation = HeatmapAnnotation( 
             "number of genes" = anno_barplot(comb_size(plotdf), 
                                              ylim = c(0, max(comb_size(plotdf))*1.1),
                                              border = FALSE, 
                                              gp = gpar(fill = "black"), 
                                              height = unit(4, "cm")), 
             annotation_name_side = "left", 
             annotation_name_rot = 90))
draw(p)
decorate_annotation("number of genes", {
  grid.text(comb_size(plotdf)[column_order(p)], x = seq_along(comb_size(plotdf)), y = unit(comb_size(plotdf)[column_order(p)], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
})

# upset plot of how many genes have LFC > 0.5 in cer/par/hyc/hyp/combinations of those 4
library(ComplexHeatmap)

lt <- list(wt_hyb = filter(DEup, type == "wt_hyb"),
           wt_parents = filter(DEup, type == "wt_parents"),
           del_hyb = filter(DEup, type == "del_hyb"),
           del_parents = filter(DEup, type == "del_parents")) |> 
  map(.f = mutate, gene_tf = paste0(gene_name, deletion)) |> 
  map(.f = select, gene_tf) |> 
  map(.f = pull)
           
# upset
plotdf <- make_comb_mat(lt)
p <- UpSet(plotdf, set_order = c("wt_hyb", "wt_parents", "del_hyb", "del_parents"),
           comb_order = order(comb_size(plotdf)), # this, I will say, makes ggplot look like a freakin' dream. I literally just wanted to rename the histogram y axis and add counts to the top of the bars
           top_annotation = HeatmapAnnotation( 
             "number of genes" = anno_barplot(comb_size(plotdf), 
                                              ylim = c(0, max(comb_size(plotdf))*1.1),
                                              border = FALSE, 
                                              gp = gpar(fill = "black"), 
                                              height = unit(4, "cm")), 
             annotation_name_side = "left", 
             annotation_name_rot = 90))
draw(p)
decorate_annotation("number of genes", {
  grid.text(comb_size(plotdf)[column_order(p)], x = seq_along(comb_size(plotdf)), y = unit(comb_size(plotdf)[column_order(p)], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
})
# repeat for down hits
lt <- list(wt_hyb = filter(DEdown, type == "wt_hyb"),
           wt_parents = filter(DEdown, type == "wt_parents"),
           del_hyb = filter(DEdown, type == "del_hyb"),
           del_parents = filter(DEdown, type == "del_parents")) |> 
  map(.f = mutate, gene_tf = paste0(gene_name, deletion)) |> 
  map(.f = select, gene_tf) |> 
  map(.f = pull)

# upset
plotdf <- make_comb_mat(lt)
p <- UpSet(plotdf, set_order = c("wt_hyb", "wt_parents", "del_hyb", "del_parents"),
           comb_order = order(comb_size(plotdf)), # this, I will say, makes ggplot look like a freakin' dream. I literally just wanted to rename the histogram y axis and add counts to the top of the bars
           top_annotation = HeatmapAnnotation( 
             "number of genes" = anno_barplot(comb_size(plotdf), 
                                              ylim = c(0, max(comb_size(plotdf))*1.1),
                                              border = FALSE, 
                                              gp = gpar(fill = "black"), 
                                              height = unit(4, "cm")), 
             annotation_name_side = "left", 
             annotation_name_rot = 90))
draw(p)
decorate_annotation("number of genes", {
  grid.text(comb_size(plotdf)[column_order(p)], x = seq_along(comb_size(plotdf)), y = unit(comb_size(plotdf)[column_order(p)], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
})




#### Hybrid high LFC troubleshooting continued: Visualizing single genes ####
TFdeldf_hybhigh <- TFdeldf_hyc[(TFdeldf_hyc$lfc > 0.5) | 
                                 (TFdeldf_hyp$lfc > 0.5),] |> 
  select(gene_name, deletion) |> 
  drop_na()

# selecting a random hyb high gene
random_gene_tf <- TFdeldf_hybhigh[sample(c(1:nrow(TFdeldf_hybhigh)), 1),]
random_gene_tf
random_gene_tf$genotype <- paste0(random_gene_tf$deletion, "delete")

# problem gene/TFs to troubleshoot:
# # YLR053C NRG1 does not look right for hyb. Supposedly doubled expr in deletion, clearly did not (unless 0 counts are causing problems)
# random_gene_tf <- tibble(gene_name = "YLR053C", deletion = "NRG1", genotype = "NRG1delete")
# # YLL061W   MIG1 shouldn't be DE for any of them (only DE for hyp)
# random_gene_tf <- tibble(gene_name = "YLL061W", deletion = "MIG1", genotype = "MIG1delete")

# How DE is it? Is it also DE in parents?
# means
list(mutate(TFdeldf_cer, type = "cer"), 
               mutate(TFdeldf_par, type = "par"), 
               mutate(TFdeldf_hyc, type = "hyc"), 
               mutate(TFdeldf_hyp, type = "hyp")) |> 
  map(.f = filter, gene_name == random_gene_tf$gene_name &
        deletion == random_gene_tf$deletion) |> 
  reduce(.f = bind_rows) |> 
  select(gene_name, deletion, lfc, type)

# jitter plot of this gene's raw counts in WT versus deletion 
# for cer/par/hyc/hyp (8 categories total)
plotdf <- list(bind_cols(tibble(expr = counts_TFdel$cer[,random_gene_tf$gene_name]),
               infos_TFdel$cer),
     bind_cols(tibble(expr = counts_TFdel$par[,random_gene_tf$gene_name]),
               infos_TFdel$par),
     bind_cols(tibble(expr = counts_TFdel_allele$cer[,random_gene_tf$gene_name]),
               infos_TFdel_allele$cer),
     bind_cols(tibble(expr = counts_TFdel_allele$par[,random_gene_tf$gene_name]),
               infos_TFdel_allele$par)) |> 
  map(.f = select, expr, organism, allele, genotype, time_point_str, well_flask_ID) |> 
  map(.f = filter, genotype %in% c("WT", random_gene_tf$genotype)) |> 
  reduce(.f = bind_rows)
plotdf$group <- paste(plotdf$organism, plotdf$allele, plotdf$genotype)
ggplot(plotdf, aes(x = group, y = expr)) + 
  geom_jitter(aes(color = group, shape = time_point_str)) +
  theme(axis.text.x = element_text(angle = 90))

# problem examples:
# YLR053C NRG1 does not look right for hyb. Supposedly doubled expr in deletion, clearly did not (unless 0 counts are causing problems)
# random_gene_tf
# plotdf |> filter(organism == "hyb" & allele == "cer") |> 
#   group_by(genotype) |> summarise(mean_expr = mean(expr))
# log2(187 + 1) - log2(108 + 1) # that's not 2.5
# 
# # YLL061W   MIG1 shouldn't be DE for any of them (only DE for hyp)
# random_gene_tf
# plotdf |> filter(organism == "hyb" & allele == "par" & time_point_str == "16 h, low N") |>
#   group_by(genotype) |> summarise(mean_expr = mean(expr), sd_expr = sd(expr))
#
# YAL062W   DAL80 only has 16h timepoint DE, but it's extremely significant in all 4

# does this still happen when we only accept genes DE at all 3 timepoints?
sum(TFdeldf_cer$DE)
sum(TFdeldf_par$DE)
sum(TFdeldf_hyc$DE)
sum(TFdeldf_hyp$DE)
sum(TFdeldf_cer$DE & TFdeldf_par$DE)
sum(TFdeldf_hyc$DE & TFdeldf_hyp$DE)
sum(TFdeldf_cer$DE & TFdeldf_hyc$DE)
sum(TFdeldf_par$DE & TFdeldf_hyp$DE)

# distribution of LFCs for all genes
plotlim <- list(TFdeldf_cer, TFdeldf_par, TFdeldf_hyc, TFdeldf_hyp) |> 
  map(.f = filter, DE) |> 
  map(.f = select, mean_lfc) |> 
  map(.f = pull) |> 
  reduce(.f = c) |> 
  quantile(probs = c(0.01, 0.99)) |> 
  abs() |> 
  max()

# cer vs par
plotdf <- left_join(TFdeldf_cer, TFdeldf_par, 
                    by = c("deletion", "gene_name"),
                    suffix = c("_cer", "_par"))
plotdf$DE_class <- if_else(plotdf$DE_cer, 
                           true = if_else(plotdf$DE_par,
                                          true = "both",
                                          false = "cerevisiae only"),
                           false = if_else(plotdf$DE_par,
                                           true = "paradoxus only",
                                           false = "neither"))
p_cerpar <- ggplot(filter(plotdf, DE_class != "neither"),
                   aes(x = mean_lfc_cer, y = mean_lfc_par)) + 
  geom_point(aes(color = DE_class), alpha = 0.5) +
  xlim(c(-plotlim, plotlim)) +
  ylim(c(-plotlim, plotlim)) +
  theme_classic() +
  geom_abline(slope = 1, intercept = 0) +
  xlab("S. cerevisiae") +
  ylab("S. paradoxus")
# hyc vs hyp
plotdf <- left_join(TFdeldf_hyc, TFdeldf_hyp, 
                    by = c("deletion", "gene_name"),
                    suffix = c("_cer", "_par"))
plotdf$DE_class <- if_else(plotdf$DE_cer, 
                           true = if_else(plotdf$DE_par,
                                          true = "both",
                                          false = "cerevisiae only"),
                           false = if_else(plotdf$DE_par,
                                           true = "paradoxus only",
                                           false = "neither"))
p_hychyp <- ggplot(filter(plotdf, DE_class != "neither"),
       aes(x = mean_lfc_cer, y = mean_lfc_par)) + 
  geom_point(aes(color = DE_class), alpha = 0.5) +
  xlim(c(-plotlim, plotlim)) +
  ylim(c(-plotlim, plotlim)) +
  theme_classic() +
  geom_abline(slope = 1, intercept = 0) +
  xlab("F1 hybrid, cerevisiae allele") +
  ylab("F1 hybrid, paradoxus allele")
# hyc vs cer
plotdf <- left_join(TFdeldf_cer, TFdeldf_hyc, 
                    by = c("deletion", "gene_name"),
                    suffix = c("_cer", "_hyc"))
plotdf$DE_class <- if_else(plotdf$DE_cer, 
                           true = if_else(plotdf$DE_hyc,
                                          true = "both",
                                          false = "parent only"),
                           false = if_else(plotdf$DE_hyc,
                                           true = "hybrid only",
                                           false = "neither"))
p_cerhyc <- ggplot(filter(plotdf, DE_class != "neither"), 
       aes(x = mean_lfc_cer, y = mean_lfc_hyc)) + 
  geom_point(aes(color = DE_class), alpha = 0.5) +
  xlim(c(-plotlim, plotlim)) +
  ylim(c(-plotlim, plotlim)) +
  theme_classic() +
  geom_abline(slope = 1, intercept = 0) +
  xlab("S. cerevisiae") +
  ylab("F1 hybrid, cerevisiae allele")
# hyp vs par
plotdf <- left_join(TFdeldf_par, TFdeldf_hyp, 
                    by = c("deletion", "gene_name"),
                    suffix = c("_par", "_hyp"))
plotdf$DE_class <- if_else(plotdf$DE_par, 
                           true = if_else(plotdf$DE_hyp,
                                          true = "both",
                                          false = "parent only"),
                           false = if_else(plotdf$DE_hyp,
                                           true = "hybrid only",
                                           false = "neither"))
p_parhyp <- ggplot(filter(plotdf, DE_class != "neither"),
       aes(x = mean_lfc_par, y = mean_lfc_hyp)) + 
  geom_point(aes(color = DE_class), alpha = 0.5) +
  xlim(c(-plotlim, plotlim)) +
  ylim(c(-plotlim, plotlim)) +
  theme_classic() +
  geom_abline(slope = 1, intercept = 0) +
  xlab("S. cerevisiae") +
  ylab("F1 hybrid, cerevisiae allele")

annotate_figure(ggarrange(p_cerpar, p_hychyp,
                          p_cerhyc, p_parhyp, nrow = 2, ncol = 2),
                top = "Expression change upon TF deletion\nall TFs, 300 random genes")


# Let's visualize some of those genes that are DE in hyc and hyp
# versus the ones DE in hyc or hyp
# do they really look DE? What about in the parents?

# hyp not hyc
random_gene_tf <- left_join(TFdeldf_hyc, TFdeldf_hyp, 
                      by = c("deletion", "gene_name"),
                      suffix = c("_hyc", "_hyp")) |> 
  filter(DE_hyp & abs(mean_lfc_hyp) > 2 & 
           !DE_hyc & abs(mean_lfc_hyc) <= 2) |> 
  ungroup() |> 
  slice_sample(n = 1) |> 
  select(gene_name, deletion)
random_gene_tf$genotype <- paste0(random_gene_tf$deletion, "delete")

# or hard code an example you're interested in
random_gene_tf <- tibble(gene_name = "YDL200C",
                         deletion = "AFT1")
random_gene_tf$genotype <- paste0(random_gene_tf$deletion, "delete")

# How DE is it? Is it also DE in parents?
TFdeldf_hyc |> filter(gene_name == random_gene_tf$gene_name &
                        deletion == random_gene_tf$deletion)
TFdeldf_hyp |> filter(gene_name == random_gene_tf$gene_name &
                        deletion == random_gene_tf$deletion)
TFdeldf_cer |> filter(gene_name == random_gene_tf$gene_name &
                        deletion == random_gene_tf$deletion)
TFdeldf_par |> filter(gene_name == random_gene_tf$gene_name &
                        deletion == random_gene_tf$deletion)

# visualizing
# Between hybrid alleles
# .cts1 = hyc del
# .cts2 = hyp del
# .cts3 = hyc WT
# .cts4 = hyp WT
p_hybrid <- plotExpressionProfileQuartet(.cts1 = counts_TFdel_allele$cer[infos_TFdel_allele$cer$genotype == random_gene_tf$genotype, random_gene_tf$gene_name, drop = FALSE],
                             .cts2 = counts_TFdel_allele$par[infos_TFdel_allele$par$genotype == random_gene_tf$genotype, random_gene_tf$gene_name, drop = FALSE],
                             .cts3 = counts_TFdel_allele$cer[infos_TFdel_allele$cer$genotype == "WT", random_gene_tf$gene_name, drop = FALSE],
                             .cts4 = counts_TFdel_allele$par[infos_TFdel_allele$par$genotype == "WT", random_gene_tf$gene_name, drop = FALSE],
                             .info1 = filter(infos_TFdel_allele$cer, genotype == random_gene_tf$genotype),
                             .info2 = filter(infos_TFdel_allele$par, genotype == random_gene_tf$genotype),
                             .info3 = filter(infos_TFdel_allele$cer, genotype == "WT"),
                             .info4 = filter(infos_TFdel_allele$par, genotype == "WT"),
                             .color1 = "orange1",
                             .color2 = "blue2",
                             .color3 = "orange4",
                             .color4 = "blue4",
                             .name1 = paste("hyc", random_gene_tf$deletion),
                             .name2 = paste("hyp", random_gene_tf$deletion),
                             .name3 = "hyc WT",
                             .name4 = "hyp WT",
                             .normalization = "log2",
                             .method = "line",
                             .show_points = TRUE,
                             .show_confidence_intervals = FALSE)

# Between parents
# .cts1 = cer del
# .cts2 = par del
# .cts3 = cer WT
# .cts4 = par WT
p_parents <- plotExpressionProfileQuartet(.cts1 = counts_TFdel$cer[infos_TFdel$cer$genotype == random_gene_tf$genotype, random_gene_tf$gene_name, drop = FALSE],
                             .cts2 = counts_TFdel$par[infos_TFdel$par$genotype == random_gene_tf$genotype, random_gene_tf$gene_name, drop = FALSE],
                             .cts3 = counts_TFdel$cer[infos_TFdel$cer$genotype == "WT", random_gene_tf$gene_name, drop = FALSE],
                             .cts4 = counts_TFdel$par[infos_TFdel$par$genotype == "WT", random_gene_tf$gene_name, drop = FALSE],
                             .info1 = filter(infos_TFdel$cer, genotype == random_gene_tf$genotype),
                             .info2 = filter(infos_TFdel$par, genotype == random_gene_tf$genotype),
                             .info3 = filter(infos_TFdel$cer, genotype == "WT"),
                             .info4 = filter(infos_TFdel$par, genotype == "WT"),
                             .color1 = "orange1",
                             .color2 = "blue2",
                             .color3 = "orange4",
                             .color4 = "blue4",
                             .name1 = paste("cer", random_gene_tf$deletion),
                             .name2 = paste("par", random_gene_tf$deletion),
                             .name3 = "cer WT",
                             .name4 = "par WT",
                             .normalization = "log2",
                             .method = "line",
                             .show_points = TRUE,
                             .show_confidence_intervals = FALSE)

annotate_figure(ggarrange(p_parents, p_hybrid),
                top = paste(random_gene_tf$gene_name,
                            random_gene_tf$deletion,
                            "Parents (left) Hybrids (right)"))
# Example of hyc not hyp, not seen in parents: YPR173C HAP1 
# looks legit. It's def expressed higher in cer than par,
# but not much response to the deletion in either species
# in hybrid, hyc allele is expressed higher than hyp in WT
# and in the deletion, it jumps even higher
# same deal with YDR289C in SOK2 mutant

# Example of hyc not hyp, reflected in parents: YKL212W (SAC1) in MSN2
# MSN2 WT hyc allele is much lower than MSN2 delete. Almost like
# MSN2 is repressing the hyc SAC1 allele
# In the cer parent, MSN2 has the same effect (repressing SAC1), not seen in par
# interestingly, the cer WT is at the same level as par WT and MSN2delete
# but the hyc WT is below the level of hyp WT/MSN2delete and hyc MSN2delete
# Why would this be? If cer WT has the same SAC1 expression as par,
# Why would hyc have abnormally low expression versus hyp?

# Example we've already looked at: TDH3 in GCR2 (positive control)
# Parents and hybrids are indistinguishable:
# par/hyp is expressed lower than cer/hyc
# all have lower expression in GCR2 mutant versus WT

# hyp vs cer (should be least similar, or comprable to cer par)
plotdf <- left_join(TFdeldf_cer, TFdeldf_hyp, 
                    by = c("deletion", "gene_name"),
                    suffix = c("_cer", "_hyp"))
plotdf$DE_class <- if_else(plotdf$DE_cer, 
                           true = if_else(plotdf$DE_hyp,
                                          true = "both",
                                          false = "parent only"),
                           false = if_else(plotdf$DE_hyp,
                                           true = "hybrid only",
                                           false = "neither"))
p_cerhyp <- ggplot(filter(plotdf, DE_class != "neither"), 
                   aes(x = mean_lfc_cer, y = mean_lfc_hyp)) + 
  geom_point(aes(color = DE_class), alpha = 0.5) +
  xlim(c(-plotlim, plotlim)) +
  ylim(c(-plotlim, plotlim)) +
  theme_classic() +
  geom_abline(slope = 1, intercept = 0) +
  xlab("S. cerevisiae") +
  ylab("F1 hybrid, paradoxus allele")
# hyc vs par
plotdf <- left_join(TFdeldf_par, TFdeldf_hyc, 
                    by = c("deletion", "gene_name"),
                    suffix = c("_par", "_hyc"))
plotdf$DE_class <- if_else(plotdf$DE_par, 
                           true = if_else(plotdf$DE_hyc,
                                          true = "both",
                                          false = "parent only"),
                           false = if_else(plotdf$DE_hyc,
                                           true = "hybrid only",
                                           false = "neither"))
p_parhyc <- ggplot(filter(plotdf, DE_class != "neither"),
                   aes(x = mean_lfc_par, y = mean_lfc_hyc)) + 
  geom_point(aes(color = DE_class), alpha = 0.5) +
  xlim(c(-plotlim, plotlim)) +
  ylim(c(-plotlim, plotlim)) +
  theme_classic() +
  geom_abline(slope = 1, intercept = 0) +
  xlab("S. paradoxus") +
  ylab("F1 hybrid, cerevisiae allele")
ggarrange(p_cerhyp, p_parhyc)

# checking hybrid for genes that are only DE in hyc or hyp but not both
# also imposing a min smd cutoff of 2
# hyc not hyp:
sum(!TFdeldf_hyc$DE & 
      TFdeldf_hyp$DE & 
      !(abs(TFdeldf_hyc$mean_lfc) > 2) & 
      (abs(TFdeldf_hyp$mean_lfc) > 2))
# hyp not hyc
sum(!TFdeldf_hyc$DE & 
      TFdeldf_hyp$DE & 
      !(abs(TFdeldf_hyc$mean_lfc) > 2) & 
      (abs(TFdeldf_hyp$mean_lfc) > 2))
# both
sum(TFdeldf_hyc$DE & 
      TFdeldf_hyp$DE & 
      (abs(TFdeldf_hyc$mean_lfc) > 2) & 
      (abs(TFdeldf_hyp$mean_lfc) > 2))
# they're relatively even categories

# are hyc not hyp genes more likely to be cer not par?
sum(TFdeldf_hyc$DE & 
      !TFdeldf_hyp$DE & 
      (abs(TFdeldf_hyc$mean_lfc) > 2) & 
      !(abs(TFdeldf_hyp$mean_lfc) > 2) &
      (abs(TFdeldf_cer$mean_lfc) > 2) &
      TFdeldf_cer$DE)
sum(TFdeldf_hyc$DE & 
      !TFdeldf_hyp$DE & 
      (abs(TFdeldf_hyc$mean_lfc) > 2) & 
      !(abs(TFdeldf_hyp$mean_lfc) > 2) &
      (abs(TFdeldf_par$mean_lfc) > 2) &
      TFdeldf_par$DE) # no, even more paradoxus hits in common than cer, but not many in either category
# what about hyp not hyc?
sum(!TFdeldf_hyc$DE & 
      TFdeldf_hyp$DE & 
      !(abs(TFdeldf_hyc$mean_lfc) > 2) & 
      (abs(TFdeldf_hyp$mean_lfc) > 2) &
      (abs(TFdeldf_par$mean_lfc) > 2) &
      TFdeldf_par$DE)
sum(!TFdeldf_hyc$DE & 
      TFdeldf_hyp$DE & 
      !(abs(TFdeldf_hyc$mean_lfc) > 2) & 
      (abs(TFdeldf_hyp$mean_lfc) > 2) &
      (abs(TFdeldf_cer$mean_lfc) > 2) &
      TFdeldf_cer$DE)
# the main conclusion is that most hybrid DE genes simply aren't shared with either parent

# hybrids do have a class of genes with quite high smd not seen in parents
# are these all in one TF or anything like that?
hyc_not_cer_table <- left_join(TFdeldf_cer, TFdeldf_hyc, 
          by = c("deletion", "gene_name"),
          suffix = c("_cer", "_hyc")) |> 
  filter(DE_cer | DE_hyc) |> 
  mutate(hyc_not_cer = if_else(DE_cer, true = if_else(DE_hyc, true = "both",
                                                      false = "cer only"),
                               false = "hyc only")) |> 
  select(deletion, hyc_not_cer) |> table()

plot(hyc_not_cer_table[,"cer only"] + hyc_not_cer_table[,"both"],
     hyc_not_cer_table[,"hyc only"], pch = NA,
     xlab = "cer only or both cer and hyc", ylab = "hyc only")
text(hyc_not_cer_table[,"cer only"] + hyc_not_cer_table[,"both"],
     hyc_not_cer_table[,"hyc only"], labels = rownames(hyc_not_cer_table))
abline(a = 0, b = 1, col = "red")

hyp_not_par_table <- left_join(TFdeldf_par, TFdeldf_hyp, 
                               by = c("deletion", "gene_name"),
                               suffix = c("_par", "_hyp")) |> 
  filter(DE_par | DE_hyp) |> 
  mutate(hyp_not_par = if_else(DE_par, true = if_else(DE_hyp, true = "both",
                                                      false = "par only"),
                               false = "hyp only")) |> 
  select(deletion, hyp_not_par) |> table()
hyp_not_par_table
plot(hyp_not_par_table[,"par only"] + hyp_not_par_table[,"both"],
     hyp_not_par_table[,"hyp only"], pch = NA,
     xlab = "par only or both par and hyp", ylab = "hyp only")
text(hyp_not_par_table[,"par only"] + hyp_not_par_table[,"both"],
     hyp_not_par_table[,"hyp only"], labels = rownames(hyp_not_par_table))
abline(a = 0, b = 1, col = "red")

# nothing obvious, might become obvious with full data
# NRG1 is the only TF with abnormally high effects in both hybrid alleles, interestingly most are different TFs for each allele
# GCR2 has abnormally low effects in hybrid alleles



#### Troubleshooting the tilt: why are LFC ratios larger magnitude in hyb TFdel than WT? ####

# can this be replicated in mean expr rather than ratios? single allele/TF/timepoint?
# example: CHA4 16h
# hyc
ggplot(filter(TFdeldf_hyc, time_point_str == "16 h, low N" & deletion == "CHA4"), 
       aes(x = log2(mean_wt + 1), y = log2(mean_del + 1))) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "gold")
# hyp
ggplot(filter(TFdeldf_hyp, time_point_str == "16 h, low N" & deletion == "CHA4"), 
       aes(x = log2(mean_wt + 1), y = log2(mean_del + 1))) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "gold")
# nope it's only apparent in the ratio:
plotdf <- bind_rows(mutate(filter(TFdeldf_hyc, time_point_str == "16 h, low N" & deletion == "CHA4"),
                           type = "hyc"),
                    mutate(filter(TFdeldf_hyp, time_point_str == "16 h, low N" & deletion == "CHA4"),
                           type = "hyp")) |> 
  pivot_wider(id_cols = "gene_name", names_from = "type", values_from = c("mean_wt", "mean_del")) |> 
  mutate(lfc_wt = log2(mean_wt_hyc + 1) - log2(mean_wt_hyp + 1),
         lfc_del = log2(mean_del_hyc + 1) - log2(mean_del_hyp + 1)) |> 
  filter((mean_del_hyc > 100 |  mean_del_hyp > 100) & 
           (mean_wt_hyc > 100 |  mean_wt_hyp > 100))
# which genes are the largest offenders? What changes about their expression from TP1?
plotdf$lfc_diff <- plotdf$lfc_del - plotdf$lfc_wt
arrange(plotdf, desc(lfc_diff))
highest_tilt <- arrange(plotdf, desc(lfc_diff)) |> slice(1:10) |> select(gene_name) |> pull()
arrange(plotdf, lfc_diff)
lowest_tilt <- arrange(plotdf, lfc_diff) |> slice(1:10) |> select(gene_name) |> pull()
plotdf$tilt <- if_else(plotdf$gene_name %in% highest_tilt, true = "high",
                       false = if_else(plotdf$gene_name %in% lowest_tilt, true = "low",
                                       false = "none"))
ggplot(plotdf, 
       aes(x = lfc_wt, 
           y = lfc_del)) + 
  geom_point(aes(color = tilt)) +
  geom_abline(slope = 1, intercept = 0, color = "gold")
# replace "16 h, low N" with "0 h, YPD" to see how this tilt isn't present at TP1

# most dramatic example:
compareQuartet(.g = "YBR085W", .del = "CHA4delete")
# in the parents, both species have lower expression by 16h in CHA4 delete
# in the hybrid, it seems that the hyc allele "takes over" all expression by 16 h

# is this a common pattern? Are these instances where one allele "takes over" expression
# for the other?
compareQuartet(highest_tilt[2], .del = "CHA4delete")
compareQuartet(highest_tilt[3], .del = "CHA4delete") # here's one where the difference is between WT alleles and CHA4 deletion causes both alleles to no longer be expressed (i.e. their difference is erased)
compareQuartet(highest_tilt[4], .del = "CHA4delete")
compareQuartet(highest_tilt[5], .del = "CHA4delete")
compareQuartet(highest_tilt[6], .del = "CHA4delete")
compareQuartet(highest_tilt[7], .del = "CHA4delete")
compareQuartet(highest_tilt[8], .del = "CHA4delete") # same as highest_tilt[3] (YDL037C): it's actually the WT par that's persisting while the rest are tanked
compareQuartet(highest_tilt[9], .del = "CHA4delete")
compareQuartet(highest_tilt[10], .del = "CHA4delete") # criss cross
# most have hyp TFdel allele tanked

# is the other end of the tilt the opposite? Is the cer allele tanked in TFdel?
compareQuartet(lowest_tilt[1], .del = "CHA4delete")
compareQuartet(lowest_tilt[2], .del = "CHA4delete") 
compareQuartet(lowest_tilt[3], .del = "CHA4delete") 
compareQuartet(lowest_tilt[4], .del = "CHA4delete")
compareQuartet(lowest_tilt[5], .del = "CHA4delete")
compareQuartet(lowest_tilt[6], .del = "CHA4delete") # this one is actually seen in the parents too
compareQuartet(lowest_tilt[7], .del = "CHA4delete")
compareQuartet(lowest_tilt[8], .del = "CHA4delete")
compareQuartet(lowest_tilt[9], .del = "CHA4delete")
compareQuartet(lowest_tilt[10], .del = "CHA4delete") # yes! all 10

# TODO: repeat this with GZF3 or TEC1 0/1h (more 1h) to make sure this same thing is also happening
# sometimes at a middle timepoint and therefore not somehow the result of the 16h WT libraries being too large
# TEC1 tilt at 1h
plotdf <- bind_rows(mutate(filter(TFdeldf_hyc, time_point_str == "1 h, low N" & deletion == "TEC1"),
                           type = "hyc"),
                    mutate(filter(TFdeldf_hyp, time_point_str == "1 h, low N" & deletion == "TEC1"),
                           type = "hyp")) |> 
  pivot_wider(id_cols = "gene_name", names_from = "type", values_from = c("mean_wt", "mean_del")) |> 
  mutate(lfc_wt = log2(mean_wt_hyc + 1) - log2(mean_wt_hyp + 1),
         lfc_del = log2(mean_del_hyc + 1) - log2(mean_del_hyp + 1))
# which genes are the largest offenders? What changes about their expression from TP1?
plotdf$lfc_diff <- plotdf$lfc_del - plotdf$lfc_wt
arrange(plotdf, desc(lfc_diff))
highest_tilt <- arrange(plotdf, desc(lfc_diff)) |> slice(1:10) |> select(gene_name) |> pull()
arrange(plotdf, lfc_diff)
lowest_tilt <- arrange(plotdf, lfc_diff) |> slice(1:10) |> select(gene_name) |> pull()
plotdf$tilt <- if_else(plotdf$gene_name %in% highest_tilt, true = "high",
                       false = if_else(plotdf$gene_name %in% lowest_tilt, true = "low",
                                       false = "none"))
ggplot(plotdf, 
       aes(x = lfc_wt, 
           y = lfc_del)) + 
  geom_point(aes(color = tilt)) +
  geom_abline(slope = 1, intercept = 0, color = "gold")
compareQuartet(highest_tilt[1], .del = "TEC1delete")
compareQuartet(highest_tilt[2], .del = "TEC1delete")
compareQuartet(highest_tilt[3], .del = "TEC1delete") # parent does this one too
compareQuartet(highest_tilt[4], .del = "TEC1delete")
compareQuartet(highest_tilt[5], .del = "TEC1delete")
compareQuartet(highest_tilt[6], .del = "TEC1delete")
compareQuartet(highest_tilt[7], .del = "TEC1delete")
compareQuartet(highest_tilt[8], .del = "TEC1delete")
compareQuartet(highest_tilt[9], .del = "TEC1delete")
compareQuartet(highest_tilt[10], .del = "TEC1delete")

compareQuartet(lowest_tilt[1], .del = "TEC1delete")
compareQuartet(lowest_tilt[2], .del = "TEC1delete")
compareQuartet(lowest_tilt[3], .del = "TEC1delete") 
compareQuartet(lowest_tilt[4], .del = "TEC1delete") # this effect isn't TF del specific and is present in parents
compareQuartet(lowest_tilt[5], .del = "TEC1delete")
compareQuartet(lowest_tilt[6], .del = "TEC1delete")
compareQuartet(lowest_tilt[7], .del = "TEC1delete")
compareQuartet(lowest_tilt[8], .del = "TEC1delete")
compareQuartet(lowest_tilt[9], .del = "TEC1delete")
compareQuartet(lowest_tilt[10], .del = "TEC1delete")

# YBR085W shows up in multiple TF deletions
compareQuartet("YBR085W", .del = "INO4delete")
compareQuartet("YBR085W", .del = "CHA4delete")
compareQuartet("YBR085W", .del = "ARG81delete")
compareQuartet("YBR085W", .del = "AFT1delete")
compareQuartet("YBR085W", .del = "GCR2delete")

# TODO: it seems to be the case that for certain genes,
# YBR085W is one,
# hyc WT is biased to be lower than all TFdels and
# hyp WT is biased to be higher than all TFdels
# (so WTs "meet in the middle")
# is this true across the board? Most genes won't
# have this problem, but when they do, is this the
# direction of bias?
# The way we tell this apart from specific TFdel responses 
# is the key that it's ALL TFdels, so treating genotype as
# 2 groups, WT and TFdel, should still capture the difference

# TODO: plot LFCs of WT vs mean across all TFdel in hyc vs hyp
# to see if A) the hyc ratios are biased to be more negative
# than hyp (or if it's an even number in each group, but 
# gene-specific), and B) if there are obvious genes that have this
# that can be separated from genes that don't
plotdf <- TFdeldf_hyc |> group_by(gene_name) |> 
  summarise(mean_lfc = mean(lfc, na.rm = TRUE),
            mean_wt_all_tp_del = mean(mean_wt, na.rm = TRUE),
            mean_del_all_tp_del = mean(mean_del, na.rm = TRUE)) |> 
  mutate(allele = "cer")
plotdf <- TFdeldf_hyp |> group_by(gene_name) |> 
  summarise(mean_lfc = mean(lfc, na.rm = TRUE),
            mean_wt_all_tp_del = mean(mean_wt, na.rm = TRUE),
            mean_del_all_tp_del = mean(mean_del, na.rm = TRUE)) |> 
  mutate(allele = "par") |> 
  bind_rows(plotdf)
plotdf$lfc_all_tp_del <- log2(plotdf$mean_del_all_tp_del + 1) -
  log2(plotdf$mean_wt_all_tp_del + 1)
plotdf <- plotdf |> pivot_wider(id_cols = "gene_name",
                                names_from = allele,
                                values_from = "lfc_all_tp_del",
                                names_prefix = "lfc_")
p <- ggplot(plotdf, aes(x = lfc_cer, y = lfc_par)) +
  geom_point(aes(color = gene_name == "YBR085W")) +
  theme(legend.position = "none")
library(ggExtra)
ggMarginal(p)
# Where I left off: YBR085W is indeed an outlier on this plot,
# it has a consistently positive LFC in cer upon any TF deletion,
# not seen in par. But I worry I've taken too many means of means
# to see true patterns. What I think is happening is that the ratio
# of counts for certain genes shifts to one allele and that this is seen
# across all TF deletions and not in WT. So a plot that better
# pairs cer and par counts that are from that same sample (so the
# ratio might be cer/par instead of wt/tfdel) but still allows
# for the fact that the phenomenon is that it's seen across
# ALL TF deletions for certain genes

#### "Drop earring" plots of single hybrid samples in single genes versus parental WT expr ratio ####

# TODO: plotdf where 
#     1) each row is a hybrid sample
#     2) there's cer/par expression ratios (LFCs) in WT parents and that hybrid sample
# then plot points parental WT LFC versus hybrid single sample LFC
# with vertical lines connecting samples from same condition (i.e. pairs of replicates)
# color = hybrid genotype == "WT"
# shape = timepoint
plotdf <- expand_grid(gene_name = colnames(counts_TFdel_allele$cer),
                      sample_name = gsub("_hyc_", "_hyb_", infos_TFdel_allele$cer$sample_name)) |> 
  left_join(y = select(infos_TFdel_allele$cer, sample_name, organism, genotype, time_point_str,
                       well_flask_ID, condition) |> 
              mutate(sample_name = gsub("_hyc_", "_hyb_", sample_name)),
            by = "sample_name", relationship = "many-to-one")
# adding reference parental WT LFC between cer and par
# (same number for each gene now)
plotdf <- map(colnames(counts_TFdel_allele$cer), \(.g) {
  cer_mean <- counts_TFdel$cer[infos_TFdel$cer$genotype == "WT", .g] |> mean(na.rm = TRUE)
  par_mean <- counts_TFdel$par[infos_TFdel$par$genotype == "WT", .g] |> mean(na.rm = TRUE)
  if (cer_mean < 30 & par_mean < 30) { # (example: log2(1) - log2(4) = -2, which seems like a massive LFC but it's a difference of 3 counts)
    return(tibble(parents_LFC = NA, 
                  gene_name = .g))
  }
  else {
    return(tibble(parents_LFC = log2(cer_mean + 1) - log2(par_mean + 1), 
                  gene_name = .g))
  }
}) |> reduce(.f = bind_rows) |> 
  right_join(y = plotdf, by = "gene_name", relationship = "one-to-many")

# just the WT samples that were in the same wellplate as the TFdel first, cause it's faster
WT_wellplate_names <- select(infos_TFdel_allele$cer, sample_name) |> 
  pull() |> 
  grep(pattern = "_hyc_WTs_", value = TRUE) |> 
  gsub(pattern = "_hyc_", replacement = "_hyb_") 
WT_withTFdels_name <- plotdf |> 
  filter(genotype == "WT") |> 
  select(sample_name) |> 
  pull() |> 
  setdiff(y = WT_wellplate_names)
plotdf_WT <- filter(plotdf, sample_name %in% WT_withTFdels_name)
plotdf_WT$hyb_LFC <- map2(plotdf_WT$sample_name, plotdf_WT$gene_name, \(.s, .g) {
  cat("working on hybrid WT", max(which(plotdf_WT$gene_name == .g)), "/", nrow(plotdf), "\n")
  hyc_count <- counts_TFdel_allele$cer[gsub("_hyb_", "_hyc_", .s), .g]
  hyp_count <- counts_TFdel_allele$par[gsub("_hyb_", "_hyp_", .s), .g]
  if (hyc_count < 30 & hyp_count < 30) { # (example: log2(1) - log2(4) = -2, which seems like a massive LFC but it's a difference of 3 counts)
    return(NA)
  }
  else {
    return(log2(hyc_count + 1) - log2(hyp_count + 1))
  }
}) |> unlist()
ggplot(filter(plotdf_WT, time_point_str == "16 h, low N"), aes(x = parents_LFC, y = hyb_LFC)) + 
  geom_point(aes(color = time_point_str))

# just one gene which had hyc counts plummet in CHA4 deletion
plotdf_YAR050W <- filter(eardropdf, gene_name == "YAR050W")
plotdf_YAR050W$hyb_LFC <- map2(plotdf_YAR050W$sample_name, plotdf_YAR050W$gene_name, \(.s, .g) {
  hyc_count <- counts_TFdel_allele$cer[gsub("_hyb_", "_hyc_", .s), .g]
  hyp_count <- counts_TFdel_allele$par[gsub("_hyb_", "_hyp_", .s), .g]
  if (hyc_count < 30 & hyp_count < 30) { # (example: log2(1) - log2(4) = -2, which seems like a massive LFC but it's a difference of 3 counts)
    return(NA)
  }
  else {
    return(log2(hyc_count + 1) - log2(hyp_count + 1))
  }
}) |> unlist()
ggplot(plotdf_YAR050W, aes(x = parents_LFC, y = hyb_LFC)) + 
  geom_jitter(aes(color = genotype == "CHA4delete")) + 
  xlim(c(-5, 5)) +
  ylim(c(-5, 5)) +
  theme(legend.position = "none")

# adding hybrid cer/par sample LFC
# (one number per sample = much longer computation time)
plotdf$hyb_LFC <- map2(plotdf$sample_name, plotdf$gene_name, \(.s, .g) {
  cat("working on hybrid", max(which(plotdf$gene_name == .g)), "/", nrow(plotdf), "\n")
  hyc_count <- counts_TFdel_allele$cer[gsub("_hyb_", "_hyc_", .s), .g]
  hyp_count <- counts_TFdel_allele$par[gsub("_hyb_", "_hyp_", .s), .g]
  if (hyc_count < 30 & hyp_count < 30) { # (example: log2(1) - log2(4) = -2, which seems like a massive LFC but it's a difference of 3 counts)
    return(NA)
  }
  else {
    return(log2(hyc_count + 1) - log2(hyp_count + 1))
  }
}) |> unlist()

# if it gets too crazy to have 5000 genes x 19 TFs x 3 timepoints drop earrings, we can
# first try plotting a random sample of drop earrings
# or collapse each gene/TF/TP combo into one earring for each gene?
# this would require taking mean cer/par WT parent across timepoints
# color within it which TFs/TPs are which?
# earrings get longer, should be same number of points though




#### TF deletions ####
# QC: do TF deletion genotypes have their corresponding TF downregulated?
# in cer
plotdf1 <- TFdel_lookup |> 
  filter(common %in% common_TFs) |> 
  left_join(bind_cols(tibble(systematic = colnames(counts_TFdel$cer)),
                                               t(counts_TFdel$cer[infos_TFdel$cer$genotype == "WT",])), by = "systematic")
plotdf1 <- pivot_longer(plotdf1, cols = colnames(plotdf1)[3:ncol(plotdf1)], names_to = "sample_name", values_to = "expr") |> 
  left_join(select(infos_TFdel$cer, sample_name, well_flask_ID), by = "sample_name")
plotdf1$genotype <- "WT"
plotdf2 <- bind_cols(tibble(genotype = infos_TFdel$cer$genotype[infos_TFdel$cer$genotype != "WT"],
                            well_flask_ID = infos_TFdel$cer$well_flask_ID[infos_TFdel$cer$genotype != "WT"],
                            sample_name = rownames(counts_TFdel$cer[infos_TFdel$cer$genotype != "WT",])),
                     counts_TFdel$cer[infos_TFdel$cer$genotype != "WT",])
plotdf2 <- select(plotdf2, c("genotype", "sample_name", "well_flask_ID", intersect(TFdel_lookup$systematic, colnames(plotdf2))))
plotdf2 <- plotdf2 |> pivot_longer(cols = colnames(plotdf2)[4:ncol(plotdf2)],
                                   names_to = "systematic", values_to = "expr") |> 
  left_join(TFdel_lookup, by = "systematic")
plotdf2 <- apply(plotdf2, 1, \(x) {
  if (x["genotype"] == "WT") {
    return(x)
  }
  if (x["genotype"] != "WT") {
    if (x["common"] != gsub("delete", "", x["genotype"])) {
      return()
    }
    if (x["common"] == gsub("delete", "", x["genotype"])) {
      return(x)
    }
  } 
}) |> bind_rows()
plotdf2$expr <- as.numeric(plotdf2$expr)
plotdf <- bind_rows(plotdf1, plotdf2)
plotdf$plot_order <- apply(plotdf, 1, \(x) {
  if (x["genotype"] == "WT") {
    return(paste(x["common"], "WT", sep = "_"))
  }
  else {
    return(x["common"])
  }
})

### supplementary figure: TF deletions repress the detection of their gene (i.e. SUM1delete has lower SUM1 expression than WT)
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/TFdel_QC_cer.pdf",
    width = 12, height = 3)
ggplot(plotdf, aes(x = plot_order, y = expr)) + 
  geom_point(aes(color = common, shape = well_flask_ID == "rep1")) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
  xlab("Transcription factor deletion") +
  ylab("Expression of deleted gene \nin WT vs deletion mutant") +
  ggtitle("S. cerevisiae")
dev.off()

# repeat for par
plotdf1 <- TFdel_lookup |> 
  filter(common %in% common_TFs) |> 
  left_join(bind_cols(tibble(systematic = colnames(counts_TFdel$par)),
                                               t(counts_TFdel$par[infos_TFdel$par$genotype == "WT",])), by = "systematic")
plotdf1 <- pivot_longer(plotdf1, cols = colnames(plotdf1)[3:ncol(plotdf1)], names_to = "sample_name", values_to = "expr") |> 
  left_join(select(infos_TFdel$par, sample_name, well_flask_ID), by = "sample_name")
plotdf1$genotype <- "WT"
plotdf2 <- bind_cols(tibble(genotype = infos_TFdel$par$genotype[infos_TFdel$par$genotype != "WT"],
                            well_flask_ID = infos_TFdel$par$well_flask_ID[infos_TFdel$par$genotype != "WT"],
                            sample_name = rownames(counts_TFdel$par[infos_TFdel$par$genotype != "WT",])),
                     counts_TFdel$par[infos_TFdel$par$genotype != "WT",])
plotdf2 <- select(plotdf2, c("genotype", "sample_name", "well_flask_ID", intersect(TFdel_lookup$systematic, colnames(plotdf2))))
plotdf2 <- plotdf2 |> pivot_longer(cols = colnames(plotdf2)[4:ncol(plotdf2)],
                                   names_to = "systematic", values_to = "expr") |> 
  left_join(TFdel_lookup, by = "systematic")
plotdf2 <- apply(plotdf2, 1, \(x) {
  if (x["genotype"] == "WT") {
    return(x)
  }
  if (x["genotype"] != "WT") {
    if (x["common"] != gsub("delete", "", x["genotype"])) {
      return()
    }
    if (x["common"] == gsub("delete", "", x["genotype"])) {
      return(x)
    }
  } 
}) |> bind_rows()
plotdf2$expr <- as.numeric(plotdf2$expr)
plotdf <- bind_rows(plotdf1, plotdf2)
plotdf$plot_order <- apply(plotdf, 1, \(x) {
  if (x["genotype"] == "WT") {
    return(paste(x["common"], "WT", sep = "_"))
  }
  else {
    return(x["common"])
  }
})
### supplementary figure: paradoxus version
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/TFdel_QC_par.pdf",
    width = 12, height = 3)
ggplot(plotdf, aes(x = plot_order, y = expr)) + 
  geom_point(aes(color = common, shape = well_flask_ID == "rep1")) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
  xlab("Transcription factor deletion") +
  ylab("Expression of deleted gene \nin WT vs deletion mutant") +
  ggtitle("S. paradoxus")
dev.off()
# Based on this, we should remove GCN4, as it seems like one replicate got contaminated with non-deletion cells in paradoxus
# (note in some TFdel iterations, we've previously removed GCN4 due to not having
# 3 timepoints x 2 replicates in cer, par, *and* hyb)
common_TFs <- setdiff(common_TFs, "GCN4")
TFdeldf_cer <- TFdeldf_cer |> filter(deletion %in% common_TFs)
TFdeldf_par <- TFdeldf_par |> filter(deletion %in% common_TFs)
TFdeldf_hyc <- TFdeldf_hyc |> filter(deletion %in% common_TFs)
TFdeldf_hyp <- TFdeldf_hyp |> filter(deletion %in% common_TFs)
TFdf <- TFdf |> filter(common %in% common_TFs)
TFenrichdf00 <- TFenrichdf00 |> filter(deletion %in% common_TFs)
random_TFenrichdf00 <- random_TFenrichdf00 |> filter(deletion %in% common_TFs)
TFenrichdf10 <- TFenrichdf10 |> filter(deletion %in% common_TFs)
random_TFenrichdf10 <- random_TFenrichdf10 |> filter(deletion %in% common_TFs)
TFenrichdf25 <- TFenrichdf25 |> filter(deletion %in% common_TFs)
random_TFenrichdf25 <- random_TFenrichdf25 |> filter(deletion %in% common_TFs)
TFenrichdf35 <- TFenrichdf35 |> filter(deletion %in% common_TFs)
random_TFenrichdf35 <- random_TFenrichdf35 |> filter(deletion %in% common_TFs)

# in addition to above visual check, are the TFs coefficients all low pvalue and negative?
# cer
for (tf in common_TFs) {
  orf <- TFdel_lookup |> filter(common == tf) |> select(systematic) |> pull()
  coeff_cer <- TFdeldf_cer |> filter(gene_name == orf & deletion == tf) |> select(coef) |> pull()
  pval_cer <- TFdeldf_cer |> filter(gene_name == orf & deletion == tf) |> select(pval) |> pull()
  cat(tf, orf, "coefficient:", coeff_cer, "pvalue:", pval_cer, "\n")
} # The empties are genes that don't have expression data (CFB1 and HAP1)
# all strongly negative

# par
for (tf in common_TFs) {
  orf <- TFdel_lookup |> filter(common == tf) |> select(systematic) |> pull()
  coeff_par <- TFdeldf_par |> filter(gene_name == orf & deletion == tf) |> select(coef) |> pull()
  pval_par <- TFdeldf_par |> filter(gene_name == orf & deletion == tf) |> select(pval) |> pull()
  cat(tf, orf, "coefficient:", coeff_par, "pvalue:", pval_par, "\n")
} # all strongly negative except INO4, it's lowly expressed, so deletion can't be distinguished from WT low expression

# visualizing how many samples per genotype. 
# Every genotype should have at least 2 replicates x 3 timepoints
table(infos_TFdel$cer$genotype)
names(table(infos_TFdel$cer$genotype))[table(infos_TFdel$cer$genotype) < 6]
table(infos_TFdel$par$genotype)
names(table(infos_TFdel$par$genotype))[table(infos_TFdel$par$genotype) < 6]
table(infos_TFdel$cer$genotype)
names(table(infos_TFdel$cer$genotype))[table(infos_TFdel$cer$genotype) < 6]
table(infos_TFdel$par$genotype)
names(table(infos_TFdel$par$genotype))[table(infos_TFdel$par$genotype) < 6]

# how many are DE in each species?
sum(abs(TFdeldf_cer$coef) > coef_thresh & TFdeldf_cer$pval < p_thresh, na.rm = TRUE) # all cer
sum(abs(TFdeldf_par$coef) > coef_thresh & TFdeldf_par$pval < p_thresh, na.rm = TRUE) # all par
sum(abs(TFdeldf_cer$coef) > coef_thresh & TFdeldf_cer$pval < p_thresh &
      abs(TFdeldf_par$coef) > coef_thresh & TFdeldf_par$pval < p_thresh, na.rm = TRUE) # both cer and par
# more of medium effect:
sum(abs(TFdeldf_cer$coef) > 1 & TFdeldf_cer$pval < p_thresh, na.rm = TRUE) # all cer
sum(abs(TFdeldf_par$coef) > 1 & TFdeldf_par$pval < p_thresh, na.rm = TRUE) # all par
sum(abs(TFdeldf_cer$coef) > 1 & TFdeldf_cer$pval < p_thresh &
      abs(TFdeldf_par$coef) > 1 & TFdeldf_par$pval < p_thresh, na.rm = TRUE) # both cer and par
# a ton of small effect:
sum(abs(TFdeldf_cer$coef) > 0.25 & TFdeldf_cer$pval < p_thresh, na.rm = TRUE) # all cer
sum(abs(TFdeldf_par$coef) > 0.25 & TFdeldf_par$pval < p_thresh, na.rm = TRUE) # all par
sum(abs(TFdeldf_cer$coef) > 0.25 & TFdeldf_cer$pval < p_thresh &
    abs(TFdeldf_par$coef) > 0.25 & TFdeldf_par$pval < p_thresh, na.rm = TRUE) # both cer and par

#### Why are so many genes only DE in one species? ####

# Basically QC section, we just want to see if
# the actual genes that are DE upon each TF deletion in each species is
# as unique as appears to be
# (and not just the result of my significance cutoff being too stringent)

# The way to see this is to check if most genes have lfc in the same direction
# regardless of significance (or maybe 0.05 significance if the lower confidence 
# lfc's are just absolutely bonkers)

# first, are there any genes with a very nonsig pvalue and high coefficient? 
sum(abs(TFdeldf_cer$coef) > 10 & TFdeldf_cer$pval >= 0.05)
sum(abs(TFdeldf_par$coef) > 10 & TFdeldf_par$pval >= 0.05) # yes there are. So we should filter them out
# the counts over the red line are our problem genes
ggplot(filter(TFdeldf_cer, pval >= 0.05), aes(x = abs(coef))) + 
  geom_histogram() + 
  geom_vline(xintercept = coef_thresh, color = "red") # there's a couple way up high, but it's hard to see them
# log scale
ggplot(filter(TFdeldf_cer, pval >= 0.05), aes(x = log(abs(coef) + 1e-9))) + 
  geom_histogram() + 
  geom_vline(xintercept = log(coef_thresh + 1e-9), color = "red")
# the counts over the red line are true significant
ggplot(filter(TFdeldf_cer, pval < 0.05), aes(x = abs(coef))) + 
  geom_histogram() + 
  geom_vline(xintercept = coef_thresh, color = "red")
# log scale
ggplot(filter(TFdeldf_cer, pval < 0.05), aes(x = log(abs(coef) + 1e-9))) + 
  geom_histogram() + 
  geom_vline(xintercept = log(coef_thresh + 1e-9), color = "red")

# plotting lfc of each gene in response to same TF deletion in cer vs par
plotdf <- bind_rows(mutate(filter(TFdeldf_cer, pval < 0.05), species = "cer"),
                    mutate(filter(TFdeldf_par, pval < 0.05), species = "par")) |> 
  pivot_wider(id_cols = c("gene_name", "deletion"), values_from = "coef",
              names_from = "species") |> drop_na() # we have to drop na for all the genes that are only pval < 0.05 in one species or the other
coef_thresh <- 2
ggplot(plotdf, aes(x = cer, y = par)) + geom_hex() + 
  geom_hline(yintercept = coef_thresh, color = "red") +
  geom_hline(yintercept = -coef_thresh, color = "red") +
  geom_vline(xintercept = coef_thresh, color = "red") +
  geom_vline(xintercept = -coef_thresh, color = "red") +
  geom_vline(xintercept = 0, color = "green3") +
  geom_hline(yintercept = 0, color = "green3")
  
# effects are decently correlated between species---most
sum(sign(plotdf$cer) == sign(plotdf$par)) # number of genes in top right and bottom left green-demarcated quadrants
sum(sign(plotdf$cer) != sign(plotdf$par)) # number of genes in top left and bottom right green-demarcated quadrants
# there are more genes with effects in the same direction for a given tf, but there are still plenty that don't
coef_thresh <- 1
sum(plotdf$cer > coef_thresh & plotdf$par > coef_thresh) + sum(plotdf$cer < -coef_thresh & plotdf$par < -coef_thresh)  # number of genes in top right and bottom left red-demarcated boxes
sum(plotdf$cer > coef_thresh & plotdf$par < -coef_thresh) + sum(plotdf$cer < -coef_thresh & plotdf$par > coef_thresh) # number of genes in top right and bottom left red-demarcated boxes
# the fraction that have directions in opposite directions drops off when you use strong thresholds (from 40% to 22% at l2fc of 1, to 10% at l2fc of 2)

# in conclusion: the genes that are DE in both species 
# are mostly in the same direction
# But, most of the DE genes are only DE in one species, 
# because the strength of most effects have changed

# for the genes that are only DE in one species, 
# is the coef (regardless of pvalue) in the same direction?
plotdf <- bind_rows(mutate(TFdeldf_cer, species = "cer"),
                    mutate(TFdeldf_par, species = "par")) |> 
  pivot_wider(id_cols = c("gene_name", "deletion"), 
              values_from = c("coef", "pval"),
              names_from = "species") |> 
  filter((pval_cer < p_thresh & abs(coef_cer) > coef_thresh) | 
           pval_par < p_thresh & abs(coef_par) > coef_thresh) # filtering for significance in at least one species
plotdf$direction <- map2(plotdf$coef_cer, plotdf$coef_par, \(x, y) {
  if ((x > coef_thresh & y > coef_thresh) | (x < -coef_thresh & y < -coef_thresh)) {
    return("DE same direction")
  }
  if ((x > coef_thresh & y < -coef_thresh) | (x < -coef_thresh & y > coef_thresh)) {
    return("DE opposite direction")
  }
  if ((x > coef_thresh) | (x < -coef_thresh)) {
    return("DE S. cerevisiae only")
  }
  if ((y > coef_thresh) | (y < -coef_thresh)) {
    return("DE S. paradoxus only")
  }
}) |> unlist()
ggplot(plotdf, aes(x = coef_cer, y = coef_par)) + 
  geom_point(aes(color = direction), alpha = 0.5) +
  scale_color_discrete(type = c("gold",
                                "mediumseagreen",
                                "orange1",
                                "blue2"), 
                       limits = c("DE same direction",
                                  "DE opposite direction",
                                  "DE S. cerevisiae only",
                                  "DE S. paradoxus only"),
                       labels = c("DE same direction",
                                  "DE opposite direction",
                                  "DE S. cerevisiae only",
                                  "DE S. paradoxus only"))
# those very crazy outliers in one species or the other have very 
# nonsig pvalues in that species---they just happen to be sig in the other species:
plotdf |> filter(abs(coef_cer) > 20) |> select(pval_cer) |> min() # huge, nonsig
plotdf |> filter(abs(coef_cer) > 20) |> select(pval_par) |> min() # tiny, sig
plotdf |> filter(abs(coef_par) > 20) |> select(pval_par) |> min() # huge, nonsig
plotdf |> filter(abs(coef_par) > 20) |> select(pval_cer) |> min() # tiny, sig
# without crazy outliers
p_point <- ggplot(filter(plotdf, abs(coef_cer) < 20 & abs(coef_par) < 20), 
       aes(x = coef_cer, y = coef_par)) + 
  geom_point(aes(color = direction), alpha = 0.5) +
  scale_color_discrete(type = c("gold",
                                "mediumseagreen",
                                "orange1",
                                "blue2"), 
                       limits = c("DE same direction",
                                  "DE opposite direction",
                                  "DE S. cerevisiae only",
                                  "DE S. paradoxus only"),
                       labels = c("DE same direction",
                                  "DE opposite direction",
                                  "DE S. cerevisiae only",
                                  "DE S. paradoxus only")) +
  theme_classic() +
  geom_hline(yintercept = coef_thresh, color = "red") +
  geom_hline(yintercept = -coef_thresh, color = "red") +
  geom_vline(xintercept = coef_thresh, color = "red") +
  geom_vline(xintercept = -coef_thresh, color = "red") +
  xlab("S. cerevisiae log2 fold change") +
  ylab("S. paradoxus log2 fold change")
p_count <- ggplot(filter(plotdf, abs(coef_cer) < 20 & abs(coef_par) < 20), 
                  aes(x = direction)) + 
  geom_bar(aes(fill = direction)) +
  scale_fill_discrete(type = c("gold",
                                "mediumseagreen",
                                "orange1",
                                "blue2"), 
                       limits = c("DE same direction",
                                  "DE opposite direction",
                                  "DE S. cerevisiae only",
                                  "DE S. paradoxus only"),
                       labels = c("DE same direction",
                                  "DE opposite direction",
                                  "DE S. cerevisiae only",
                                  "DE S. paradoxus only")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
p <- ggarrange(p_point, p_count, nrow = 1, ncol = 2, common.legend = TRUE,
          legend = "left", widths = c(2, 1))
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/TFdel_lfc_correlation_parents.pdf",
    width = 8, height = 5)
annotate_figure(p, top = "Log2 fold change upon TF deletion, S. paradoxus vs. S. cerevisiae\nall genes, all TF deletions")
dev.off()
table(plotdf$direction)

# are the lfc's correlated?
test <- filter(plotdf, abs(coef_cer) < 20 & 
                 abs(coef_par) < 20)
plot(test$coef_cer, test$coef_par)
cor(test$coef_cer, test$coef_par)
cor(test$coef_cer, test$coef_par, method = "spearman") # no

# let's visualize a random example just to make sure nothing weird's going on
# in cer
random_cerDE <- plotdf |> 
  filter(direction == "DE S. cerevisiae only" & abs(coef_cer) > 2 & pval_cer < 1e-5) |> 
  select(gene_name, deletion) |> 
  slice_sample(n = 1)
plotdf |> filter(gene_name == random_cerDE$gene_name & deletion == random_cerDE$deletion)
gdf <- bind_rows(mutate(bind_cols(infos_TFdel$cer, tibble(expr = counts_TFdel$cer[,random_cerDE$gene_name])), species = "cer"),
                    mutate(bind_cols(infos_TFdel$par, tibble(expr = counts_TFdel$par[,random_cerDE$gene_name])), species = "par"))
ggplot(filter(gdf, species == "cer"), aes(x = genotype, y = expr)) + 
  geom_point(aes(color = genotype == paste0(random_cerDE$deletion, "delete"),
                 shape = time_point_str)) +
  theme(legend.position = "bottom") +
  ggtitle(paste(random_cerDE$gene_name, "expression in cer\n for all TFdels,\n focus on", random_cerDE$deletion))
ggplot(filter(gdf, species == "par"), aes(x = genotype, y = expr)) + 
  geom_point(aes(color = genotype == paste0(random_cerDE$deletion, "delete"),
                 shape = time_point_str)) +
  theme(legend.position = "bottom") +
  ggtitle(paste(random_cerDE$gene_name, "expression in par \nfor all TFdels,\n focus on", random_cerDE$deletion))

# how many samples? Should be at least 6 per species
# cer
filter(gdf, species == "cer" & genotype == paste0(random_cerDE$deletion, "delete")) |> 
  select(sample_name, time_point_str, expr)
# par
filter(gdf, species == "par" & genotype == paste0(random_cerDE$deletion, "delete")) |> 
  select(sample_name, time_point_str, expr)

# in par
random_parDE <- plotdf |> 
  filter(direction == "DE S. paradoxus only" & abs(coef_par) > 2 & pval_par < 1e-5) |> 
  select(gene_name, deletion) |> 
  slice_sample(n = 1)
plotdf |> filter(gene_name == random_parDE$gene_name & deletion == random_parDE$deletion)
gdf <- bind_rows(mutate(bind_cols(infos_TFdel$cer, tibble(expr = counts_TFdel$cer[,random_parDE$gene_name])), species = "cer"),
                 mutate(bind_cols(infos_TFdel$par, tibble(expr = counts_TFdel$par[,random_parDE$gene_name])), species = "par"))
ggplot(filter(gdf, species == "cer"), aes(x = genotype, y = expr)) + 
  geom_point(aes(color = genotype == paste0(random_parDE$deletion, "delete"),
                 shape = time_point_str)) + 
  theme(axis.text.x = element_text(angle = 90), legend.position = "none") +
  ggtitle(paste(random_parDE$gene_name, "expression in cer \nfor all TFdels,\n focus on", random_parDE$deletion))
ggplot(filter(gdf, species == "par"), aes(x = genotype, y = expr)) + 
  geom_point(aes(color = genotype == paste0(random_parDE$deletion, "delete"),
                 shape = time_point_str)) + 
  theme(axis.text.x = element_text(angle = 90), legend.position = "none") +
  ggtitle(paste(random_parDE$gene_name, "expression in par \nfor all TFdels,\n focus on", random_parDE$deletion))
# how many samples?
# cer
filter(gdf, species == "cer" & genotype == paste0(random_parDE$deletion, "delete")) |> 
  select(sample_name, time_point_str, expr) 
# par
filter(gdf, species == "par" & genotype == paste0(random_parDE$deletion, "delete")) |> 
  select(sample_name, time_point_str, expr)

#### repeating test of lfc correlation for hybrid alleles ####

# is the trans-regulatory environment responsible for how different the lfc's 
# of each gene are btwn species?
plotdf <- bind_rows(mutate(TFdeldf_hyc, allele = "cer"),
                    mutate(TFdeldf_hyp, allele = "par")) |> 
  pivot_wider(id_cols = c("gene_name", "deletion"), 
              values_from = c("coef", "pval"),
              names_from = "allele") |> 
  filter((pval_cer < p_thresh & abs(coef_cer) > coef_thresh) | 
           pval_par < p_thresh & abs(coef_par) > coef_thresh) # filtering for significance in at least one species
plotdf$direction <- map2(plotdf$coef_cer, plotdf$coef_par, \(x, y) {
  if ((x > coef_thresh & y > coef_thresh) | (x < -coef_thresh & y < -coef_thresh)) {
    return("DE same direction")
  }
  if ((x > coef_thresh & y < -coef_thresh) | (x < -coef_thresh & y > coef_thresh)) {
    return("DE opposite direction")
  }
  if ((x > coef_thresh) | (x < -coef_thresh)) {
    return("DE S. cerevisiae only")
  }
  if ((y > coef_thresh) | (y < -coef_thresh)) {
    return("DE S. paradoxus only")
  }
}) |> unlist()
ggplot(plotdf, aes(x = coef_cer, y = coef_par)) + 
  geom_point(aes(color = direction), alpha = 0.5) +
  scale_color_discrete(type = c("gold",
                                "mediumseagreen",
                                "orange1",
                                "blue2"), 
                       limits = c("DE same direction",
                                  "DE opposite direction",
                                  "DE S. cerevisiae only",
                                  "DE S. paradoxus only"),
                       labels = c("DE same direction",
                                  "DE opposite direction",
                                  "DE S. cerevisiae only",
                                  "DE S. paradoxus only"))
# those very crazy outliers in one species or the other have very 
# nonsig pvalues in that species---they just happen to be sig in the other species:
plotdf |> filter(abs(coef_cer) > 20) |> select(pval_cer) |> min() # huge, nonsig
plotdf |> filter(abs(coef_cer) > 20) |> select(pval_par) |> min() # tiny, sig
plotdf |> filter(abs(coef_par) > 20) |> select(pval_par) |> min() # huge, nonsig
plotdf |> filter(abs(coef_par) > 20) |> select(pval_cer) |> min() # tiny, sig
# without crazy outliers
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/TFdel_lfc_correlation_hybrid.pdf",
    width = 8, height = 5)
p_point <- ggplot(filter(plotdf, abs(coef_cer) < 20 & abs(coef_par) < 20), 
                  aes(x = coef_cer, y = coef_par)) + 
  geom_point(aes(color = direction), alpha = 0.5) +
  scale_color_discrete(type = c("gold",
                                "mediumseagreen",
                                "orange1",
                                "blue2"), 
                       limits = c("DE same direction",
                                  "DE opposite direction",
                                  "DE S. cerevisiae only",
                                  "DE S. paradoxus only"),
                       labels = c("DE same direction",
                                  "DE opposite direction",
                                  "DE S. cerevisiae only",
                                  "DE S. paradoxus only")) +
  theme_classic() +
  geom_hline(yintercept = coef_thresh, color = "red") +
  geom_hline(yintercept = -coef_thresh, color = "red") +
  geom_vline(xintercept = coef_thresh, color = "red") +
  geom_vline(xintercept = -coef_thresh, color = "red") +
  xlab("S. cerevisiae log2 fold change") +
  ylab("S. paradoxus log2 fold change")
p_count <- ggplot(filter(plotdf, abs(coef_cer) < 20 & abs(coef_par) < 20), 
                  aes(x = direction)) + 
  geom_bar(aes(fill = direction)) +
  scale_fill_discrete(type = c("gold",
                               "mediumseagreen",
                               "orange1",
                               "blue2"), 
                      limits = c("DE same direction",
                                 "DE opposite direction",
                                 "DE S. cerevisiae only",
                                 "DE S. paradoxus only"),
                      labels = c("DE same direction",
                                 "DE opposite direction",
                                 "DE S. cerevisiae only",
                                 "DE S. paradoxus only")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
p <- ggarrange(p_point, p_count, nrow = 1, ncol = 2, common.legend = TRUE,
               legend = "left", widths = c(2, 1))
annotate_figure(p, top = "Log2 fold change upon TF deletion, \nS. paradoxus vs. S. cerevisiae F1 hybrid alleles\nall genes, all TF deletions")
dev.off()
table(plotdf$direction) # TODO: consider adding names to the outliers like my pattern vs level single gene level plots


# TODO: also visualize random pairs again

# TODO: upset plots for genes DE up or down 
# (maybe two different plots, one up one down)
# for cer, par, hyc, hyp as the 4 venn diagram categories
# to see how often the same gene is DE in cer and hyc or par and hyp
# (we already know most DE genes in cer are not the same DE genes in par)


#### Probably Supplement: Identifying candidate trans regulators of CCMs ####
# I'll probably relegate this to the supplement b/c the barplots are so much more intuitive.
# But it is useful for me b/c it lets me see between-species effect sizes
# and estimate the rate at which DE genes tend to arise (by comparing block size and nDE including random modules)

# How often is the set of DE genes enriched in each CCM?
# For each TF deletion, plot the effect size of DE genes by CCM
library(ggbeeswarm)
plotDeletionEffectByModule <- function(.del, .species = c("cer", "par")) {
  if (.species == "cer") {
    coefmat <- TFdel_coef_cer
    pvalmat <- TFdel_pval_cer
  }
  if (.species == "par") {
    coefmat <- TFdel_coef_par
    pvalmat <- TFdel_pval_par
  }
  plotdf <- tibble(CCM = character(0),
                   effect_size = numeric(0),
                   pvalues = numeric(0))
  for (clr in unique(moduledf[moduledf$is_CCM,]$CCM_color)) {
    gene_names <- module_genedf$gene_name[module_genedf$CCM_color == clr]
    effect_sizes <- coefmat |> filter(gene_name %in% gene_names) |> select(all_of(.del)) |> pull()
    pvalues <- pvalmat |> filter(gene_name %in% gene_names) |> select(all_of(.del)) |> pull()
    plotdf <- bind_rows(plotdf, tibble(CCM = rep(clr, length(effect_sizes)), effect_size = effect_sizes, pvalues = pvalues))
  }
  plotdf <- plotdf[plotdf$pvalues < 1e-5,] |> arrange(CCM)
  p <- ggplot(plotdf, aes(x = CCM, y = effect_size)) +
    geom_beeswarm(color = plotdf$CCM) +
    geom_hline(yintercept = 0, color = "orange") +
    theme_classic() +
    xlab("") + 
    ylab("Log Fold Change") + 
    ggtitle(paste(.del, .species))
  return(p)
}
plotDeletionEffectByModule("DAL80", .species = "cer") # seems red and pink genes increase expr in response to its deletion
plotDeletionEffectByModule("DAL80", .species = "par") # red and pink don't have the same relationship in par

# TODO: before splitting into CCMs, which deletions affect the most genes and 
# which have the biggest difference of DE genes between cer and par?
# just plot cer vs par where each point is one deletion
common_deletions <- intersect(colnames(TFdel_coef_cer[-1]), colnames(TFdel_coef_par[-1]))
plotdf <- tibble(deletion = common_deletions)
plotdf$nsig_cer <- map_dbl(plotdf$deletion, \(del) sum(TFdel_pval_cer[, del] < 0.05))
plotdf$nsig_par <- map_dbl(plotdf$deletion, \(del) sum(TFdel_pval_par[, del] < 0.05))
ggplot(plotdf, aes(x = nsig_cer, y = nsig_par, label = deletion)) + 
  geom_label() +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  theme_classic() +
  xlab("Number of DE genes in cerevisiae") +
  ylab("Number of DE genes in paradoxus")

# met28 and zap1 are wildly more influential in paradoxus
plotDeletionEffectByModule("MET28", "par")
plotDeletionEffectByModule("MET28", "cer") # that's wild
plotDeletionEffectByModule("ZAP1", "par")
plotDeletionEffectByModule("ZAP1", "cer") # so is that
# hap3 slightly more in cer
plotDeletionEffectByModule("HAP3", "cer")
plotDeletionEffectByModule("HAP3", "par") # but both are very strong in pink
# SUM1 really only in cer
plotDeletionEffectByModule("SUM1", "cer")
plotDeletionEffectByModule("SUM1", "par") # but CDA1 is still sig up in par, we know from stats tests above (CDA1 isn't in a CCM though I believe)

# getting plots for all TF deletions in both species to ID
# TFs that might regulate whole CCMs
# midnightblue and yellow CCMs especially, as they're the most divergent
# cer
deletions <- colnames(TFdel_coef_cer)[-1]
plotlist_cer <- vector(mode = "list", length = length(deletions))
names(plotlist_cer) <- deletions
for (del in deletions) {
  plotlist_cer[[del]] <- plotDeletionEffectByModule(.del = del, .species = "cer")
}
library(ggpubr)
ggarrange(plotlist = plotlist_cer[1:3], ncol = 1, nrow = 3)
ggarrange(plotlist = plotlist_cer[4:6], ncol = 1, nrow = 3)
ggarrange(plotlist = plotlist_cer[7:9], ncol = 1, nrow = 3)
ggarrange(plotlist = plotlist_cer[10:12], ncol = 1, nrow = 3)
ggarrange(plotlist = plotlist_cer[13:15], ncol = 1, nrow = 3)
ggarrange(plotlist = plotlist_cer[16:18], ncol = 1, nrow = 3)
ggarrange(plotlist = plotlist_cer[19:21], ncol = 1, nrow = 3)
ggarrange(plotlist = plotlist_cer[22:24], ncol = 1, nrow = 3)
ggarrange(plotlist = plotlist_cer[25:27], ncol = 1, nrow = 3) # PHO4 down in midnightblue
ggarrange(plotlist = plotlist_cer[28:30], ncol = 1, nrow = 3) 
ggarrange(plotlist = plotlist_cer[31:33], ncol = 1, nrow = 3)
ggarrange(plotlist = plotlist_cer[34:36], ncol = 1, nrow = 3)
ggarrange(plotlist = plotlist_cer[37:39], ncol = 1, nrow = 3) # SUM1 up in brown
ggarrange(plotlist = plotlist_cer[40:42], ncol = 1, nrow = 3) # URE2 up in blue, green, purple, red

# par
deletions <- colnames(TFdel_coef_par)[-1]
plotlist_par <- vector(mode = "list", length = length(deletions))
names(plotlist_par) <- deletions
for (del in deletions) {
  plotlist_par[[del]] <- plotDeletionEffectByModule(.del = del, .species = "par")
}

ggarrange(plotlist = plotlist_par[1:3], ncol = 1, nrow = 3)
ggarrange(plotlist = plotlist_par[4:6], ncol = 1, nrow = 3) # pink and red maybe up in BAS1?
ggarrange(plotlist = plotlist_par[7:9], ncol = 1, nrow = 3)
ggarrange(plotlist = plotlist_par[10:12], ncol = 1, nrow = 3)
ggarrange(plotlist = plotlist_par[13:15], ncol = 1, nrow = 3) # magenta down GZF3
ggarrange(plotlist = plotlist_par[16:18], ncol = 1, nrow = 3) # pink (and black?) way down in all HAPs
ggarrange(plotlist = plotlist_par[19:21], ncol = 1, nrow = 3)
ggarrange(plotlist = plotlist_par[22:24], ncol = 1, nrow = 3)
ggarrange(plotlist = plotlist_par[25:27], ncol = 1, nrow = 3)
ggarrange(plotlist = plotlist_par[28:30], ncol = 1, nrow = 3) # brown and salmon down in RPN4
ggarrange(plotlist = plotlist_par[31:33], ncol = 1, nrow = 3)
ggarrange(plotlist = plotlist_par[34:36], ncol = 1, nrow = 3)
ggarrange(plotlist = plotlist_par[37:39], ncol = 1, nrow = 3) # blue, magenta, red up in URE2
ggarrange(plotlist = plotlist_par[40:41], ncol = 1, nrow = 2)

# Null model - random modules
# same as plotDeletionEffectByModule but for random modules instead of CCMs
plotDeletionEffectByRandomModule <- function(.del, .species = c("cer", "par")) {
  if (.species == "cer") {
    coefmat <- TFdel_coef_cer
    pvalmat <- TFdel_pval_cer
  }
  if (.species == "par") {
    coefmat <- TFdel_coef_par
    pvalmat <- TFdel_pval_par
  }
  plotdf <- tibble(random_module = character(0),
                   effect_size = numeric(0),
                   pvalues = numeric(0))
  for (mod in moduledf[moduledf$type == "Random",]$module_name) {
    gene_names <- module_genedf$gene_name[module_genedf$random_module == mod]
    effect_sizes <- coefmat |> filter(gene_name %in% gene_names) |> select(.del) |> pull()
    pvalues <- pvalmat |> filter(gene_name %in% gene_names) |> select(.del) |> pull()
    plotdf <- bind_rows(plotdf, tibble(random_module = rep(mod, length(effect_sizes)), effect_size = effect_sizes, pvalues = pvalues, gene_name = gene_names))
  }
  plotdf <- plotdf[plotdf$pvalues < 1e-5,] |> arrange(random_module)
  p <- ggplot(plotdf, aes(x = random_module, y = effect_size)) +
    geom_beeswarm(aes(color = random_module)) +
    geom_hline(yintercept = 0, color = "orange") +
    theme_classic() +
    xlab("") +
    ylab("Log Fold Change") +
    ggtitle(paste(.del, .species))
  return(p)
}
plotDeletionEffectByRandomModule("SUM1", .species = "cer")
plotDeletionEffectByRandomModule("PHO4", .species = "cer")

# creating 2 dataframes (1 cer, 1 par) of effect size distributions for
# each random and CCM module and use the null distributions
# to ID which CCM distributions are asymmetric
# 7 columns: gene name, CCM or random, module name, block size, deletion, effect size, pvalue
# each row is one gene, either as CCM or random (so nrow = 2*nCCMGenes*42deletions)

# Old version that threw out genes without CCMs
# module_genedf_just_CCMs <- module_genedf |> filter(is_CCM)
# block_size_lookup <- tibble(module_name = names(table(module_genedf_just_CCMs$CCM_color)),
#                             block_size = as.numeric(table(module_genedf_just_CCMs$CCM_color))) |> 
#   bind_rows(tibble(module_name = names(table(module_genedf_just_CCMs$random_module)),
#                    block_size = as.numeric(table(module_genedf_just_CCMs$random_module))))
# # cer
# TFdeldf_cer_coef <- module_genedf_just_CCMs |> select(gene_name, CCM_color, random_module) |> 
#   left_join(TFdel_coef_cer, by = "gene_name")
# TFdeldf_cer_coef$name <- "coef"
# TFdeldf_cer_pval <- module_genedf |> filter(is_CCM) |> select(gene_name, CCM_color, random_module) |> 
#   left_join(TFdel_pval_cer, by = "gene_name")
# TFdeldf_cer_pval$name <- "pval"
# TFdeldf_cer <- bind_rows(TFdeldf_cer_coef, TFdeldf_cer_pval) |> 
#   pivot_longer(cols = colnames(TFdel_coef_cer[,-1]), names_to = "deletion") |> 
#   pivot_wider(id_cols = c("gene_name", "CCM_color", "random_module", "deletion")) |> 
#   pivot_longer(cols = c("CCM_color", "random_module"), names_to = "CCM_or_random", values_to = "module_name")
# TFdeldf_cer$CCM_or_random <- map_chr(TFdeldf_cer$CCM_or_random, \(x) strsplit(x, split = "_")[[1]][1])
# TFdeldf_cer <- left_join(TFdeldf_cer, block_size_lookup, by = "module_name")
# # par
# TFdeldf_par_coef <- module_genedf_just_CCMs |> select(gene_name, CCM_color, random_module) |> 
#   left_join(TFdel_coef_par, by = "gene_name")
# TFdeldf_par_coef$name <- "coef"
# TFdeldf_par_pval <- module_genedf |> filter(is_CCM) |> select(gene_name, CCM_color, random_module) |> 
#   left_join(TFdel_pval_par, by = "gene_name")
# TFdeldf_par_pval$name <- "pval"
# TFdeldf_par <- bind_rows(TFdeldf_par_coef, TFdeldf_par_pval) |> 
#   pivot_longer(cols = colnames(TFdel_coef_par[,-1]), names_to = "deletion") |> 
#   pivot_wider(id_cols = c("gene_name", "CCM_color", "random_module", "deletion")) |> 
#   pivot_longer(cols = c("CCM_color", "random_module"), names_to = "CCM_or_random", values_to = "module_name")
# TFdeldf_par$CCM_or_random <- map_chr(TFdeldf_par$CCM_or_random, \(x) strsplit(x, split = "_")[[1]][1])
# TFdeldf_par <- left_join(TFdeldf_par, block_size_lookup, by = "module_name")

# using block size to estimate DE gene discovery rate and identify CCMs with DE enrichment
# cer
plotdf_cer <- TFdeldf_cer |> filter(is_CCM) |> 
  select(c("gene_name", "deletion", "CCM_color", "pval", "coef", "random_module")) |> 
  pivot_longer(cols = c("CCM_color", "random_module"), names_to = "CCM_or_random", values_to = "module_name")
plotdf_cer$CCM_or_random <- map_chr(plotdf_cer$CCM_or_random, \(x) strsplit(x, split = "_")[[1]][1])
# random block sizes include non-CCM genes, so we need to recalculate them
block_size_lookup <- table(plotdf_cer$module_name)/length(unique(plotdf_cer$deletion))
block_size_lookup
plotdf_cer$block_size <- apply(plotdf_cer, 1, \(x) {
  modname <- x["module_name"]
  return(block_size_lookup[which(names(block_size_lookup) == modname)] |> as.numeric())
}) |> unlist()
table(plotdf_cer$gene_name) |> table() # verifying each gene appears 2 * nDeletions (42 for cer) times, once for its random module and once for its CCM
total_nsig <- sum(filter(TFdel_pval_cer, gene_name %in% plotdf_cer$gene_name) < 0.05) # total number of sig genes, verify that it's equal to the nsig sums below
total_nsig
plotdf_cer <- plotdf_cer |> 
  group_by(deletion, module_name, CCM_or_random, block_size) |> 
  summarise(nsig = sum(pval < 0.05)) # we're looking at false discovery rate, so no p-value correction
# verifying the total number of sig genes is the same for CCMs and random, as the same set of genes are counted in both
plotdf_cer |> filter(CCM_or_random == "CCM") |> select(nsig) |> pull() |> sum()
plotdf_cer |> filter(CCM_or_random == "random") |> select(nsig) |> pull() |> sum()
total_nsig # should be equal to two sums above ^
plotdf_cer$point_color <- apply(plotdf_cer, 1, \(x) {
  if (x["CCM_or_random"] == "CCM") {
    return(x["module_name"])
  }
  else {
    return("purple2")
  }
}) |> as.character()
plotlist_cer <- vector(mode = "list", length = length(unique(plotdf_cer$deletion)))
names(plotlist_cer) <- unique(plotdf_cer$deletion)
for (del in unique(plotdf_cer$deletion)) {
  plotdf_del <- plotdf_cer |> filter(deletion == del)
  plotlist_cer[[del]] <- ggplot(plotdf_del, aes(x = block_size, y = nsig)) + 
    geom_point(color = plotdf_del$point_color, aes(shape = CCM_or_random)) +
    theme_classic() + ggtitle(del)
}
ggarrange(plotlist = plotlist_cer[1:4], nrow = 2, ncol = 2, common.legend = TRUE, legend = "none")
ggarrange(plotlist = plotlist_cer[5:8], nrow = 2, ncol = 2, common.legend = TRUE, legend = "none") # red in DAL80
ggarrange(plotlist = plotlist_cer[9:12], nrow = 2, ncol = 2, common.legend = TRUE, legend = "none")
ggarrange(plotlist = plotlist_cer[13:16], nrow = 2, ncol = 2, common.legend = TRUE, legend = "none") # pink in HAP2, maybe tan and greenyellow in GLN3
ggarrange(plotlist = plotlist_cer[17:20], nrow = 2, ncol = 2, common.legend = TRUE, legend = "none") # pink in HAP4 and HAP5
ggarrange(plotlist = plotlist_cer[21:24], nrow = 2, ncol = 2, common.legend = TRUE, legend = "none") # red in MIG1
ggarrange(plotlist = plotlist_cer[25:28], nrow = 2, ncol = 2, common.legend = TRUE, legend = "none") 
ggarrange(plotlist = plotlist_cer[29:32], nrow = 2, ncol = 2, common.legend = TRUE, legend = "none") # pink and tan in RPN4
ggarrange(plotlist = plotlist_cer[33:36], nrow = 2, ncol = 2, common.legend = TRUE, legend = "none") # pink in RTG3 and SKN7, midnightblue in SKN7
ggarrange(plotlist = plotlist_cer[37:40], nrow = 2, ncol = 2, common.legend = TRUE, legend = "none") # purple and red in URE2
ggarrange(plotlist = plotlist_cer[41:42], nrow = 2, ncol = 1, common.legend = TRUE, legend = "none") # pink and tan in ZAP1

# checking interesting examples (with p-value correction)
# red
plotDeletionEffectByModule("MIG1", .species = "cer") # red way down
plotDeletionEffectByModule("NRG1", .species = "cer") # red way down <- MIG1 and NRG1 repress glucose metabolism (Seok Hee Park et al. 1999)
plotDeletionEffectByModule("PHD1", .species = "cer") # not really anything for red
plotDeletionEffectByModule("RGT1", .species = "cer") # red maybe a little up
# greenyellow and tan
plotDeletionEffectByModule("GLN3", .species = "cer") # not that much to go on

# purple in URE2
plotDeletionEffectByModule("URE2", .species = "cer") # purple and red def up

# midnightblue in PHO4
plotDeletionEffectByModule("PHO4", .species = "cer") # hmm yeah

# repeat for par
plotdf_par <- TFdeldf_par |> filter(is_CCM) |> 
  select(c("gene_name", "deletion", "CCM_color", "pval", "coef", "random_module")) |> 
  pivot_longer(cols = c("CCM_color", "random_module"), names_to = "CCM_or_random", values_to = "module_name")
plotdf_par$CCM_or_random <- map_chr(plotdf_par$CCM_or_random, \(x) strsplit(x, split = "_")[[1]][1])
# random block sizes include non-CCM genes, so we need to recalculate them
block_size_lookup <- table(plotdf_par$module_name)/length(unique(plotdf_par$deletion))
block_size_lookup
plotdf_par$block_size <- apply(plotdf_par, 1, \(x) {
  modname <- x["module_name"]
  return(block_size_lookup[which(names(block_size_lookup) == modname)] |> as.numeric())
}) |> unlist()
table(plotdf_par$gene_name) |> table() # verifying each gene appears 2 * nDeletions (41 for par) times, once for its random module and once for its CCM
total_nsig <- sum(filter(TFdel_pval_par, gene_name %in% plotdf_par$gene_name) < 0.05) # total number of sig genes, verify that it's equal to the nsig sums below
total_nsig
plotdf_par <- plotdf_par |> 
  group_by(deletion, module_name, CCM_or_random, block_size) |> 
  summarise(nsig = sum(pval < 0.05)) # we're looking at false discovery rate, so no p-value correction
# verifying the total number of sig genes is the same for CCMs and random, as the same set of genes are counted in both
plotdf_par |> filter(CCM_or_random == "CCM") |> select(nsig) |> pull() |> sum()
plotdf_par |> filter(CCM_or_random == "random") |> select(nsig) |> pull() |> sum()
total_nsig # should be equal to two sums above ^plotdf_par$point_color <- apply(plotdf_par, 1, \(x) {
plotdf_par$point_color <- apply(plotdf_par, 1, \(x) {
  if (x["CCM_or_random"] == "CCM") {
    return(x["module_name"])
  }
  else {
    return("purple2")
  }
}) |> as.character()
plotlist_par <- vector(mode = "list", length = length(unique(plotdf_par$deletion)))
names(plotlist_par) <- unique(plotdf_par$deletion)
for (del in unique(plotdf_par$deletion)) {
  plotdf_del <- plotdf_par |> filter(deletion == del)
  plotlist_par[[del]] <- ggplot(plotdf_del, aes(x = block_size, y = nsig)) + 
    geom_point(color = plotdf_del$point_color, aes(shape = CCM_or_random)) +
    theme_classic() + ggtitle(del)
}
ggarrange(plotlist = plotlist_par[1:4], nrow = 2, ncol = 2, common.legend = TRUE, legend = "none") # green in BAS1 still maybe
ggarrange(plotlist = plotlist_par[5:8], nrow = 2, ncol = 2, common.legend = TRUE, legend = "none") # pink in CHA4, midnightblue in DAL80
ggarrange(plotlist = plotlist_par[9:12], nrow = 2, ncol = 2, common.legend = TRUE, legend = "none") # magenta in GCN4 and GLN3, green in GCR2 and GLN3
ggarrange(plotlist = plotlist_par[13:16], nrow = 2, ncol = 2, common.legend = TRUE, legend = "none") # pink in HAP2 and HAP3
ggarrange(plotlist = plotlist_par[17:20], nrow = 2, ncol = 2, common.legend = TRUE, legend = "none") # pink in HAP4 and HAP5
ggarrange(plotlist = plotlist_par[21:24], nrow = 2, ncol = 2, common.legend = TRUE, legend = "none") # green in MSN2
ggarrange(plotlist = plotlist_par[25:28], nrow = 2, ncol = 2, common.legend = TRUE, legend = "none") # red and midnightblue in RIM101
ggarrange(plotlist = plotlist_par[29:32], nrow = 2, ncol = 2, common.legend = TRUE, legend = "none") # pink in ROX1, green in RTG1, pink tan and black in RPN4 (and it's for real)
ggarrange(plotlist = plotlist_par[33:36], nrow = 2, ncol = 2, common.legend = TRUE, legend = "none") 
ggarrange(plotlist = plotlist_par[37:40], nrow = 2, ncol = 2, common.legend = TRUE, legend = "none") # magenta in URE2
ggarrange(plotlist = plotlist_par[41], nrow = 1, ncol = 1, common.legend = TRUE, legend = "none")
# checking interesting cases (now with multiple testing correction)
# magenta URE2
plotDeletionEffectByModule("URE2", .species = "par") # magenta blue and red all up
plotDeletionEffectByModule("URE2", .species = "cer") # as a reminder, red and blue are also seen in cereivisae but not so much magenta and there's purple
# green in GCR2
plotDeletionEffectByModule("GCR2", .species = "par") # nah, just a lot
# blank and pink definitely have some conserved hap thing
plotDeletionEffectByModule("HAP2", .species = "par")
plotDeletionEffectByModule("HAP3", .species = "par")
plotDeletionEffectByModule("HAP4", .species = "par")
plotDeletionEffectByModule("HAP5", .species = "par")
plotDeletionEffectByModule("HAP2", .species = "cer") # but not quite as pronounced in cer
plotDeletionEffectByModule("HAP3", .species = "cer")
plotDeletionEffectByModule("HAP4", .species = "cer")
plotDeletionEffectByModule("HAP5", .species = "cer")

# TODO: we can also plot block size vs effect size to look
# at magnitude and direction of CCMs versus random effect sizes, but doesn't seem worth it at this
# juncture as the beeswarm plots combined with having a null expectation for discovery rate
# seem to be sufficient for screening for interesting examples

# example of how to visualize effect size directions and magnitudes 
# next to random modules of comparable sizes

# TODO: haven't checked this since adding non-CCMs to the TFdeldfs
# TODO: could make this a function
# test case, URE2 deletion ups magenta CCM in par not cer
# plotdf <- TFdeldf_par |> filter(deletion == "URE2" & pval < 1e-5)
# plotdf$point_color <- apply(plotdf, 1, \(x) {
#   if (x["CCM_or_random"] == "CCM") {
#     return(x["module_name"])
#   }
#   else {
#     return("purple2")
#   }
# }) |> as.character()
# plotdf$block_size_order <- dense_rank(plotdf$block_size)
# ggplot(plotdf, aes(x = block_size_order, y = coef)) + 
#   geom_jitter(color = plotdf$point_color) +
#   geom_hline(yintercept = 0, color = "orange") +
#   theme_classic() + xlab("Module size order (smallest to largest)") +
#   ylab("Log-fold change per gene \nin response to TF deletion") +
#   ggtitle("URE2 paradoxus")
# # URE2 in cerevisiae
# plotdf <- TFdeldf_cer |> filter(deletion == "URE2" & pval < 1e-5)
# plotdf$point_color <- apply(plotdf, 1, \(x) {
#   if (x["CCM_or_random"] == "CCM") {
#     return(x["module_name"])
#   }
#   else {
#     return("purple2")
#   }
# }) |> as.character()
# plotdf$block_size_order <- dense_rank(plotdf$block_size)
# ggplot(plotdf, aes(x = block_size_order, y = coef)) + 
#   geom_jitter(color = plotdf$point_color) +
#   geom_hline(yintercept = 0, color = "orange") +
#   theme_classic() + xlab("Module size (smallest to largest)") +
#   ylab("Log-fold change per gene \nin response to TF deletion") +
#   ggtitle("URE2 cerevisiae") # magenta no longer has many DE. Very hard to see but purple module (colored "purple", not "purple2") has a lot DE in cer

#### Expression profiles of test cases in WT vs deletion ####
### 1) global diverged: MET28 deletion causing massive upregulation across modules in par not cer
# comparing par WT versus MET28 delete in any module
gene_idxs <- module_genedf |> filter(module_name_cer == "d") |> select(gene_name) |> pull()
plotExpressionProfilePair(counts_TFdel$par[infos_TFdel$par$genotype == "WT", gene_idxs],
                          counts_TFdel$par[infos_TFdel$par$genotype == "MET28delete", gene_idxs],
                          infos_TFdel$par[infos_TFdel$par$genotype == "WT",],
                          infos_TFdel$par[infos_TFdel$par$genotype == "MET28delete",],
                          .name1 = "WT", .name2 = "MET28delete",
                          .method = "line", .show_points = TRUE,
                          .normalization = "log2") # try with CCMs a-p, paradoxus tends to be higher (especially e)
# cer for comparison, shouldn't have many DE genes in response to MET28 deletion
plotExpressionProfilePair(counts_TFdel$cer[infos_TFdel$cer$genotype == "WT", gene_idxs],
                          counts_TFdel$cer[infos_TFdel$cer$genotype == "MET28delete", gene_idxs],
                          infos_TFdel$cer[infos_TFdel$cer$genotype == "WT",],
                          infos_TFdel$cer[infos_TFdel$cer$genotype == "MET28delete",],
                          .name1 = "WT", .name2 = "MET28delete",
                          .method = "line", .show_points = TRUE,
                          .normalization = "log2")
# doesn't appear super different

# 2) module-specific conserved: HAP2-5 deletion causing upregulation of pink (and black?) CCM
gene_idxs <- module_genedf |> filter(module_name_cer %in% c("l", "k")) |> select(gene_name) |> pull()
# par
# HAP1
plotExpressionProfilePair(counts_TFdel$par[infos_TFdel$par$genotype == "WT", gene_idxs],
                          counts_TFdel$par[infos_TFdel$par$genotype == "HAP1delete", gene_idxs],
                          infos_TFdel$par[infos_TFdel$par$genotype == "WT",],
                          infos_TFdel$par[infos_TFdel$par$genotype == "HAP1delete",],
                          .name1 = "WT", .name2 = "HAP1delete",
                          .method = "line", .show_points = TRUE,
                          .normalization = "log2")
# HAP2
plotExpressionProfilePair(counts_TFdel$par[infos_TFdel$par$genotype == "WT", gene_idxs],
                          counts_TFdel$par[infos_TFdel$par$genotype == "HAP2delete", gene_idxs],
                          infos_TFdel$par[infos_TFdel$par$genotype == "WT",],
                          infos_TFdel$par[infos_TFdel$par$genotype == "HAP2delete",],
                          .name1 = "WT", .name2 = "HAP2delete",
                          .method = "line", .show_points = TRUE,
                          .normalization = "log2")
# HAP3
plotExpressionProfilePair(counts_TFdel$par[infos_TFdel$par$genotype == "WT", gene_idxs],
                          counts_TFdel$par[infos_TFdel$par$genotype == "HAP3delete", gene_idxs],
                          infos_TFdel$par[infos_TFdel$par$genotype == "WT",],
                          infos_TFdel$par[infos_TFdel$par$genotype == "HAP3delete",],
                          .name1 = "WT", .name2 = "HAP3delete",
                          .method = "line", .show_points = TRUE,
                          .normalization = "log2")
# HAP4
plotExpressionProfilePair(counts_TFdel$par[infos_TFdel$par$genotype == "WT", gene_idxs],
                          counts_TFdel$par[infos_TFdel$par$genotype == "HAP4delete", gene_idxs],
                          infos_TFdel$par[infos_TFdel$par$genotype == "WT",],
                          infos_TFdel$par[infos_TFdel$par$genotype == "HAP4delete",],
                          .name1 = "WT", .name2 = "HAP4delete",
                          .method = "line", .show_points = TRUE,
                          .normalization = "log2")
# HAP5
plotExpressionProfilePair(counts_TFdel$par[infos_TFdel$par$genotype == "WT", gene_idxs],
                          counts_TFdel$par[infos_TFdel$par$genotype == "HAP5delete", gene_idxs],
                          infos_TFdel$par[infos_TFdel$par$genotype == "WT",],
                          infos_TFdel$par[infos_TFdel$par$genotype == "HAP5delete",],
                          .name1 = "WT", .name2 = "HAP5delete",
                          .method = "line", .show_points = TRUE,
                          .normalization = "log2")
# cer
# HAP1
plotExpressionProfilePair(counts_TFdel$cer[infos_TFdel$cer$genotype == "WT", gene_idxs],
                          counts_TFdel$cer[infos_TFdel$cer$genotype == "HAP1delete", gene_idxs],
                          infos_TFdel$cer[infos_TFdel$cer$genotype == "WT",],
                          infos_TFdel$cer[infos_TFdel$cer$genotype == "HAP1delete",],
                          .name1 = "WT", .name2 = "HAP1delete",
                          .method = "line", .show_points = TRUE,
                          .normalization = "log2")
# HAP2
plotExpressionProfilePair(counts_TFdel$cer[infos_TFdel$cer$genotype == "WT", gene_idxs],
                          counts_TFdel$cer[infos_TFdel$cer$genotype == "HAP2delete", gene_idxs],
                          infos_TFdel$cer[infos_TFdel$cer$genotype == "WT",],
                          infos_TFdel$cer[infos_TFdel$cer$genotype == "HAP2delete",],
                          .name1 = "WT", .name2 = "HAP2delete",
                          .method = "line", .show_points = TRUE,
                          .normalization = "log2")
# HAP3
plotExpressionProfilePair(counts_TFdel$cer[infos_TFdel$cer$genotype == "WT", gene_idxs],
                          counts_TFdel$cer[infos_TFdel$cer$genotype == "HAP3delete", gene_idxs],
                          infos_TFdel$cer[infos_TFdel$cer$genotype == "WT",],
                          infos_TFdel$cer[infos_TFdel$cer$genotype == "HAP3delete",],
                          .name1 = "WT", .name2 = "HAP3delete",
                          .method = "line", .show_points = TRUE,
                          .normalization = "log2")
# HAP4
plotExpressionProfilePair(counts_TFdel$cer[infos_TFdel$cer$genotype == "WT", gene_idxs],
                          counts_TFdel$cer[infos_TFdel$cer$genotype == "HAP4delete", gene_idxs],
                          infos_TFdel$cer[infos_TFdel$cer$genotype == "WT",],
                          infos_TFdel$cer[infos_TFdel$cer$genotype == "HAP4delete",],
                          .name1 = "WT", .name2 = "HAP4delete",
                          .method = "line", .show_points = TRUE,
                          .normalization = "log2")
# HAP5
plotExpressionProfilePair(counts_TFdel$cer[infos_TFdel$cer$genotype == "WT", gene_idxs],
                          counts_TFdel$cer[infos_TFdel$cer$genotype == "HAP5delete", gene_idxs],
                          infos_TFdel$cer[infos_TFdel$cer$genotype == "WT",],
                          infos_TFdel$cer[infos_TFdel$cer$genotype == "HAP5delete",],
                          .name1 = "WT", .name2 = "HAP5delete",
                          .method = "line", .show_points = TRUE,
                          .normalization = "log2")
# doesn't seem to be as strong in cerevisiae, especially HAP4, but hard to tell if that's a significant difference

# 3) maybe module-specific diverged: GLN3 deletion upregulates greenyellow and tan in both cer and par but maybe more in cer 
gene_idxs <- module_genedf |> filter(module_name_cer %in% c("e", "d")) |> select(gene_name) |> pull()
# cer
plotExpressionProfilePair(counts_TFdel$cer[infos_TFdel$cer$genotype == "WT", gene_idxs],
                          counts_TFdel$cer[infos_TFdel$cer$genotype == "GLN3delete", gene_idxs],
                          infos_TFdel$cer[infos_TFdel$cer$genotype == "WT",],
                          infos_TFdel$cer[infos_TFdel$cer$genotype == "GLN3delete",],
                          .name1 = "WT", .name2 = "GLN3delete",
                          .method = "line", .show_points = TRUE,
                          .normalization = "log2")
# par
plotExpressionProfilePair(counts_TFdel$par[infos_TFdel$par$genotype == "WT", gene_idxs],
                          counts_TFdel$par[infos_TFdel$par$genotype == "GLN3delete", gene_idxs],
                          infos_TFdel$par[infos_TFdel$par$genotype == "WT",],
                          infos_TFdel$par[infos_TFdel$par$genotype == "GLN3delete",],
                          .name1 = "WT", .name2 = "GLN3delete",
                          .method = "line", .show_points = TRUE,
                          .normalization = "log2") # if anything it's more in paradoxus

# 4) module-specific with no expression divergence: URE2 deletion ups magenta module in par and purple module in cer
# magenta
gene_idxs <- module_genedf |> filter(module_name_cer %in% c("j")) |> select(gene_name) |> pull()
# cer
plotExpressionProfilePair(counts_TFdel$cer[infos_TFdel$cer$genotype == "WT", gene_idxs],
                          counts_TFdel$cer[infos_TFdel$cer$genotype == "URE2delete", gene_idxs],
                          infos_TFdel$cer[infos_TFdel$cer$genotype == "WT",],
                          infos_TFdel$cer[infos_TFdel$cer$genotype == "URE2delete",],
                          .name1 = "WT", .name2 = "URE2delete",
                          .method = "line", .show_points = TRUE,
                          .normalization = "log2")
# par
plotExpressionProfilePair(counts_TFdel$par[infos_TFdel$par$genotype == "WT", gene_idxs],
                          counts_TFdel$par[infos_TFdel$par$genotype == "URE2delete", gene_idxs],
                          infos_TFdel$par[infos_TFdel$par$genotype == "WT",],
                          infos_TFdel$par[infos_TFdel$par$genotype == "URE2delete",],
                          .name1 = "WT", .name2 = "URE2delete",
                          .method = "line", .show_points = TRUE,
                          .normalization = "log2")
# purple
gene_idxs <- module_genedf |> filter(module_name_cer %in% c("c")) |> select(gene_name) |> pull()
# cer
plotExpressionProfilePair(counts_TFdel$cer[infos_TFdel$cer$genotype == "WT", gene_idxs],
                          counts_TFdel$cer[infos_TFdel$cer$genotype == "URE2delete", gene_idxs],
                          infos_TFdel$cer[infos_TFdel$cer$genotype == "WT",],
                          infos_TFdel$cer[infos_TFdel$cer$genotype == "URE2delete",],
                          .name1 = "WT", .name2 = "URE2delete",
                          .method = "line", .show_points = TRUE,
                          .normalization = "log2")
# par
plotExpressionProfilePair(counts_TFdel$par[infos_TFdel$par$genotype == "WT", gene_idxs],
                          counts_TFdel$par[infos_TFdel$par$genotype == "URE2delete", gene_idxs],
                          infos_TFdel$par[infos_TFdel$par$genotype == "WT",],
                          infos_TFdel$par[infos_TFdel$par$genotype == "URE2delete",],
                          .name1 = "WT", .name2 = "URE2delete",
                          .method = "line", .show_points = TRUE,
                          .normalization = "log2") # both species' purple up in URE2 delete

# 5) module-specific and maybe some expression divergence: PHO4 downs midnightblue in cer not par
gene_idxs <- module_genedf |> filter(module_name_cer %in% c("o")) |> select(gene_name) |> pull()
# cer
plotExpressionProfilePair(counts_TFdel$cer[infos_TFdel$cer$genotype == "WT", gene_idxs],
                          counts_TFdel$cer[infos_TFdel$cer$genotype == "PHO4delete", gene_idxs],
                          infos_TFdel$cer[infos_TFdel$cer$genotype == "WT",],
                          infos_TFdel$cer[infos_TFdel$cer$genotype == "PHO4delete",],
                          .name1 = "WT", .name2 = "PHO4delete",
                          .method = "line", .show_points = TRUE,
                          .normalization = "log2")
# par
plotExpressionProfilePair(counts_TFdel$par[infos_TFdel$par$genotype == "WT", gene_idxs],
                          counts_TFdel$par[infos_TFdel$par$genotype == "PHO4delete", gene_idxs],
                          infos_TFdel$par[infos_TFdel$par$genotype == "WT",],
                          infos_TFdel$par[infos_TFdel$par$genotype == "PHO4delete",],
                          .name1 = "WT", .name2 = "PHO4delete",
                          .method = "line", .show_points = TRUE,
                          .normalization = "log2") # maybe lower in par

# important quality control: Do the genes labeled as DE in the glm.nb models actually
# look DE when their expression profiles are compared between WT and deletion?
# for the module specific TFdel effects that we want to check expression
# profiles of, label the genes that are identified as DE in the expression profiles and see
# if they appear to be changing expression profiles in the same direction as their effect sizes
# (Example: DAL80 downregulates some genes in par midnightblue while the midnightblue module appears to go UP in DAL80 delete overall)
inspectTFdelTestsAcrossTimepoints <- function(.gene_idxs, .deletion, .eff_thresh,  .species = c("cer", "par")) {
  .info <- infos_TFdel[[.species]]
  .counts <- counts_TFdel[[.species]][, .gene_idxs]
  gdf <- bind_cols(.info, .counts) |>  
    pivot_longer(cols = colnames(.counts), 
                 names_to = "gene_name", values_to = "expr") |> 
    filter(genotype %in% c("WT", paste0(.deletion, "delete"))) |> 
    group_by(genotype, gene_name, time_point_str) |>
    summarise(mean_expr = mean(expr))
  if (.species == "cer") {
    .TFdel_coef <- TFdel_coef_cer
    .TFdel_pval <- TFdel_pval_cer
  }
  if (.species == "par") {
    .TFdel_coef <- TFdel_coef_par
    .TFdel_pval <- TFdel_pval_par
  }
  gene_tests <- .TFdel_coef[.TFdel_coef$gene_name %in% .gene_idxs, 
                            which(colnames(.TFdel_coef) %in% c("gene_name", .deletion))]
  gene_tests$pval <- .TFdel_pval[.TFdel_pval$gene_name %in% .gene_idxs, 
                                 which(colnames(.TFdel_pval) == .deletion)]
  gdf <- left_join(gdf, gene_tests, by = "gene_name")
  gdf <- rename(gdf, "TFdel"=.deletion)
  if (sign(.eff_thresh) == -1) {
    gdf_sigdown <- filter(gdf, TFdel < .eff_thresh & pval < 1e-5)
    gdf_nonsig <- filter(gdf, TFdel > .eff_thresh & pval > 1e-5)
    p1_title <- c("Significantly downregulated genes")
  }
  if (sign(.eff_thresh) == 1) {
    gdf_sigdown <- filter(gdf, TFdel > .eff_thresh & pval < 1e-5)
    gdf_nonsig <- filter(gdf, TFdel < .eff_thresh & pval > 1e-5)
    p1_title <- c("Significantly upregulated genes")
  }
  # genes detected as significantly downregulated in HAP3 deletion
  p1 <- ggplot(gdf_sigdown, aes(x = time_point_str, y = log2(mean_expr + 1), group = paste(gene_name, genotype))) +
    geom_line(aes(color = genotype)) + theme_classic() +
    ggtitle(p1_title)
  # genes not detected
  p2 <- ggplot(gdf_nonsig, aes(x = time_point_str, y = log2(mean_expr + 1), group = paste(gene_name, genotype))) +
    geom_line(aes(color = genotype)) + theme_classic() +
    ggtitle("Unaltered genes")
  # same gene in WT and deletion, averaged across timepoints
  gdf_sigdown_singlemeasure <- gdf_sigdown |> group_by(genotype, gene_name) |> summarise(mean_expr = mean(mean_expr))
  gdf_sigdown_singlemeasure$sig <- TRUE
  gdf_nonsig_singlemeasure <- gdf_nonsig |> group_by(genotype, gene_name) |> summarise(mean_expr = mean(mean_expr))
  gdf_nonsig_singlemeasure$sig <- FALSE
  gdf_singlemeasure <- bind_rows(gdf_sigdown_singlemeasure, gdf_nonsig_singlemeasure)
  p3 <- ggplot(gdf_singlemeasure, aes(x = genotype, y = log2(mean_expr + 1), group = gene_name)) +
    geom_line(aes(color = sig)) + theme_classic() +
    ggtitle("All genes in module")
  return(list(sig = p1, nonsig = p2, all = p3))
}

# test examples:
# 1) DE test matches expression profile: HAP3 pink either species
# par has slightly more
gene_idxs <- module_genedf |> filter(CCM_color == "pink") |> select(gene_name) |> pull()
plots <- inspectTFdelTestsAcrossTimepoints(.gene_idxs = gene_idxs, 
                                           .deletion = "HAP3",
                                           .eff_thresh = -1,
                                           .species = "par")
ggarrange(plotlist = plots, legend = "bottom", nrow = 1, ncol = 3)
# checks out
# Is this a real module effect? How many genes in a random subset of the same size are this DE?
random_idxs <- sample(module_genedf$gene_name, size = length(gene_idxs), replace = FALSE)
sum(random_idxs %in% gene_idxs)
plots <- inspectTFdelTestsAcrossTimepoints(.gene_idxs = random_idxs, 
                                           .deletion = "HAP3",
                                           .eff_thresh = -1,
                                           .species = "par")
ggarrange(plotlist = plots, legend = "bottom", nrow = 1, ncol = 3) # There can be some, but not very many

# 2A) DE test shows significant effect not replicated in expression profile: midnightblue DE down in DAL80delete in par
gene_idxs <- module_genedf |> filter(CCM_color == "midnightblue") |> select(gene_name) |> pull()
plots <- inspectTFdelTestsAcrossTimepoints(.gene_idxs = gene_idxs, 
                                           .deletion = "DAL80",
                                           .eff_thresh = -.1,
                                           .species = "par")
ggarrange(plotlist = plots, legend = "bottom", nrow = 1, ncol = 3)
# well now there's not even any DE genes detected, probably because we raised the threshold

# 2B) DE test shows significant effect not replicated in expression profile: midnightblue DE down in PHO4delete in cer
gene_idxs <- module_genedf |> filter(CCM_color == "midnightblue") |> select(gene_name) |> pull()
plots <- inspectTFdelTestsAcrossTimepoints(.gene_idxs = gene_idxs, 
                                           .deletion = "PHO4",
                                           .eff_thresh = -.1,
                                           .species = "cer")
ggarrange(plotlist = plots, legend = "bottom", nrow = 1, ncol = 3) |>  annotate_figure(top = "PHO4 cerevisiae, midnightblue CCM")
# there's a strong effect on a few genes, just not enough of them to affect average expression
# check random sample of same size
random_idxs <- sample(module_genedf$gene_name, size = length(gene_idxs), replace = FALSE)
sum(random_idxs %in% gene_idxs)
plots <- inspectTFdelTestsAcrossTimepoints(.gene_idxs = random_idxs, 
                                           .deletion = "PHO4",
                                           .eff_thresh = -.1,
                                           .species = "cer")
ggarrange(plotlist = plots, legend = "bottom", nrow = 1, ncol = 3) |>  annotate_figure(top = "PHO4 cerevisiae, random genes")
# how about selecting all the genes that are down in PHO4 delete regardless of module?
TFdeldf_cer |> filter(deletion == "PHO4" & coef < -coef_thresh & pval < p_thresh) |> select(gene_name, CCM_color, cer_color)

# 3) DE test shows no effect but expression profile appears different: purple module up in URE2delete in paradoxus (only sig DE in cerevisiae)
gene_idxs <- module_genedf |> filter(CCM_color == "purple") |> select(gene_name) |> pull()
plots <- inspectTFdelTestsAcrossTimepoints(.gene_idxs = gene_idxs, 
                                           .deletion = "URE2",
                                           .eff_thresh = .1,
                                           .species = "cer")
# cerevisiae, expected to show an effect
ggarrange(plotlist = plots, legend = "bottom", nrow = 1, ncol = 3) |>  annotate_figure(top = "URE2 cerevisiae, purple CCM")
# moderately visible in cerevisiae
plots <- inspectTFdelTestsAcrossTimepoints(.gene_idxs = gene_idxs, 
                                           .deletion = "URE2",
                                           .eff_thresh = .1,
                                           .species = "par")
# paradoxus, not expected to show an effect
ggarrange(plotlist = plots, legend = "bottom", nrow = 1, ncol = 3) |>  annotate_figure(top = "URE2 paradoxus, purple CCM")
# not visible in paradoxus

# In summary, we can trust the DE tests and if only a few genes out of a module of 50-100 genes
# are DE, it's probably not going to affect average expression of the module
# But the question remains, if X/67 or so genes are strongly DE in a module, is that more than expected by chance?
# Does that count as a module-specific response?

#### Barplots: Which groups of genes affected by a certain TFdel are located in the same module more than expected by chance? ####

# It seems that subsets of modules can be strongly affected by a TFdel without affecting the rest of the module
# Is this simply because by chance 5 or so DE genes happened to be in a module together?
# Or is this a consistent observation: that subsets of modules can be regulated differently than the rest of the module? 
# (I mean yes almost certainly 100+ genes have differences in their regulation),
# but first we need to eliminate null possibility with Exact Tests or some other test

# how about one enrichment test first for whether or not more DE genes are in CCMs than not?
# first just checking proportions:
# fraction of DE hits in CCMs:
inspectTFDEproportion <- function(.deldf_cer, .deldf_par, .m_genedf) {
  DE_cer <- .deldf_cer |> filter(abs(coef) > coef_thresh & pval < p_thresh)
  nDE_in_CCMs_cer <- DE_cer |> left_join(y = select(.m_genedf, gene_name, is_CCM)) |> 
    select(is_CCM) |> pull() |> sum()
  DE_par <- .deldf_par |> filter(abs(coef) > coef_thresh & pval < p_thresh)
  nDE_in_CCMs_par <- DE_par |> left_join(y = select(.m_genedf, gene_name, is_CCM)) |> 
    select(is_CCM) |> pull() |> sum()
  return(list(prop_TFDE_in_CCMs_cer = nDE_in_CCMs_cer/nrow(DE_cer),
              prop_TFDE_in_CCMs_par = nDE_in_CCMs_par/nrow(DE_par)))
}
coef_thresh <- 2
# merge 00
inspectTFDEproportion(TFdeldf_cer, TFdeldf_par, module_genedf00)
sum(module_genedf00$coexpressed == "conserved co-expressed")/nrow(module_genedf00)
# merge 10
inspectTFDEproportion(TFdeldf_cer, TFdeldf_par, module_genedf10)
sum(module_genedf10$coexpressed == "conserved co-expressed")/nrow(module_genedf10)
# merge 25
inspectTFDEproportion(TFdeldf_cer, TFdeldf_par, module_genedf25)
sum(module_genedf25$coexpressed == "conserved co-expressed")/nrow(module_genedf25)
# merge 35
inspectTFDEproportion(TFdeldf_cer, TFdeldf_par, module_genedf35)
sum(module_genedf35$coexpressed == "conserved co-expressed")/nrow(module_genedf35)
# for all merge levels, the proportion of DE genes in CCMs is slightly higher 
# than the proportion of all genes in CCMs, so slight enrichment
# (it's even more strong when there's a stronger coef_thresh cutoff)

# stacked barplots of which module each of the DE genes in each deletion is from so we can see them all on one plot
coef_thresh <- 2
plotTFBarplots <- function(.deldf_cer, .deldf_par, .m_genedf,
                           .useLiteralColors = TRUE) {
  plotdf_cer <-  .deldf_cer |> left_join(y = select(.m_genedf, gene_name, CCM_color),
                                         by = "gene_name")
  plotdf_cer$sig <- apply(plotdf_cer, 1, \(x) {
    coeff <- x["coef"] |> as.numeric()
    pvalue <- x["pval"] |> as.numeric()
    if (coeff > coef_thresh & pvalue < p_thresh) {
      return("up")
    }
    if (coeff < -coef_thresh & pvalue < p_thresh) {
      return("down")
    }
    else {
      return("none")
    }
  })
  # repeat for paradoxus
  plotdf_par <-  .deldf_par |> left_join(y = select(.m_genedf, gene_name, CCM_color),
                                         by = "gene_name")
  plotdf_par$sig <- apply(plotdf_par, 1, \(x) {
    coeff <- x["coef"] |> as.numeric()
    pvalue <- x["pval"] |> as.numeric()
    if (coeff > coef_thresh & pvalue < p_thresh) {
      return("up") 
    }
    if (coeff < -coef_thresh & pvalue < p_thresh) {
      return("down")
    }
    else {
      return("none")
    }
  })
  # separating genes found in both species from genes that are species unique
  plotdf_common <- bind_rows(select(plotdf_cer, gene_name, CCM_color, deletion, sig), 
                             select(plotdf_par, gene_name, CCM_color, deletion, sig)) |> filter(sig != "none")
  plotdf_common <- plotdf_common[duplicated(plotdf_common),] |> 
    mutate(idx =  paste0(gene_name, deletion))
  # filtering out these gene/TF combinations in cer and par
  plotdf_cer <- plotdf_cer |> filter(sig != "none") |> 
    mutate(idx = paste0(gene_name, deletion)) |>
    filter(!(idx %in% plotdf_common$idx))
  plotdf_par <- plotdf_par |> filter(sig != "none") |> 
    mutate(idx = paste0(gene_name, deletion)) |>
    filter(!(idx %in% plotdf_common$idx))
  if (.useLiteralColors) {
    plotdf_cer$CCM_color <- factor(plotdf_cer$CCM_color) |> relevel(ref = "none")
    plotdf_par$CCM_color <- factor(plotdf_par$CCM_color) |> relevel(ref = "none")
    plotdf_common$CCM_color <- factor(plotdf_common$CCM_color) |> relevel(ref = "none")
  }
  # counting
  plotdf_cer <- plotdf_cer |> 
    group_by(CCM_color, deletion) |> 
    summarise(up = sum(sig == "up"), down = -sum(sig == "down"))
  plotdf_par <- plotdf_par |> 
    group_by(CCM_color, deletion) |> 
    summarise(up = sum(sig == "up"), down = -sum(sig == "down"))
  plotdf_common <- plotdf_common |> 
    group_by(CCM_color, deletion) |> 
    summarise(up = sum(sig == "up"), down = -sum(sig == "down"))
  # getting overall max, to ensure equivalent y-axes between plots
  eachtf <- bind_rows(mutate(plotdf_cer, plot = "cer"),
                      mutate(plotdf_par, plot = "par"),
                      mutate(plotdf_common, plot = "common")) |> 
                        group_by(deletion, plot) |> 
    summarise(up = sum(up), down = sum(down)) 
  max_y <- max(abs(c(eachtf$up, eachtf$down)))
  cat(max_y)
  # plotting
  # genes with a TF response only in cer
  p_cer <- plotdf_cer |> 
    ggplot(aes(x = deletion, y = up, fill = CCM_color)) +
    geom_col(position = "stack") +
    geom_col(aes(y = down), position = "stack") + 
    geom_hline(yintercept = 0) +
    ylab("") +
    xlab("Transcription Factor Deletion") +
    scale_x_discrete(limits = common_TFs, labels = common_TFs) +
    #theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none") + 
    ggtitle("genes only regulated in S. cereivisiae") +
    ylim(c(-max_y - 5, max_y + 5))
  # genes with a TF response only in par
  p_par <- plotdf_par |> 
    ggplot(aes(x = deletion, y = up, fill = CCM_color)) +
    geom_col(position = "stack") +
    geom_col(aes(y = down), position = "stack") + 
    geom_hline(yintercept = 0) +
    ylab("") +
    xlab("Transcription Factor Deletion") +
    scale_x_discrete(limits = common_TFs, labels = common_TFs) +
    #theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none") + 
    ggtitle("genes only regulated in S. paradoxus") +
    ylim(c(-max_y - 5, max_y + 5))
  # genes with conserved TF responses in cer and par
  p_common <- plotdf_common |> 
    ggplot(aes(x = deletion, y = up, fill = CCM_color)) +
    geom_col(position = "stack") +
    geom_col(aes(y = down), position = "stack") + 
    geom_hline(yintercept = 0) +
    ylab("") +
    xlab("Transcription Factor Deletion") +
    scale_x_discrete(limits = common_TFs, labels = common_TFs) +
    #theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none") + 
    ggtitle("Genes with conserved regulation between species") +
    ylim(c(-max_y - 5, max_y + 5))
  if (.useLiteralColors) {
    p_cer <- p_cer +
      scale_fill_discrete(limits = plotdf_cer$CCM_color, 
                          type = gsub("none", "grey80", plotdf_cer$CCM_color))
    p_par <- p_par +
      scale_fill_discrete(limits = plotdf_par$CCM_color, 
                          type = gsub("none", "grey80", plotdf_par$CCM_color))
    p_common <- p_common +
      scale_fill_discrete(limits = plotdf_common$CCM_color, 
                          type = gsub("none", "grey80", plotdf_common$CCM_color)) 
  }
  return(list(cer = p_cer, par = p_par, common = p_common))
}

### barplots of DE genes per TF, colored by module: 
# merge 00
p_00 <- plotTFBarplots(TFdeldf_cer, TFdeldf_par, module_genedf00)
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/TFdelBarplots00.pdf",
    width = 6, height = 10)
p <- ggarrange(p_00$cer, 
          p_00$par, 
          p_00$common, nrow = 3, ncol = 1)
annotate_figure(p, left = "number of genes up/down in response to each TF deletion")
dev.off()

# merge 10
p_10 <- plotTFBarplots(TFdeldf_cer, TFdeldf_par, module_genedf10)
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/TFdelBarplots10.pdf",
    width = 6, height = 10)
p <- ggarrange(p_10$cer, 
               p_10$par, 
               p_10$common, nrow = 3, ncol = 1)
annotate_figure(p, left = "number of genes up/down in response to each TF deletion")

dev.off()

# merge 25
p_25 <- plotTFBarplots(TFdeldf_cer, TFdeldf_par, module_genedf25)
pdf("../../aligning_the_molecular_phenotype/paper_figures/TFdel/TFdelBarplots25.pdf",
    width = 6, height = 10)
p <- ggarrange(p_25$cer, 
               p_25$par, 
               p_25$common, nrow = 3, ncol = 1)
annotate_figure(p, left = "number of genes up/down in response to each TF deletion")
dev.off()

# merge 35
p_35 <- plotTFBarplots(TFdeldf_cer, TFdeldf_par, module_genedf35)
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/TFdelBarplots35.pdf",
    width = 6, height = 10)
p <- ggarrange(p_35$cer, 
               p_35$par, 
               p_35$common, nrow = 3, ncol = 1)
annotate_figure(p, left = "number of genes up/down in response to each TF deletion")
dev.off()

# repeating with random modules to see what baseline/null module enrichment looks like
# random 00
p_r00 <- plotTFBarplots(TFdeldf_cer, TFdeldf_par, 
                        mutate(random_module_genedf00,
                        CCM_color = paste(cer_color, par_color)),
                        .useLiteralColors = FALSE)
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/TFdelBarplots_random00.pdf",
    width = 6, height = 10)
p <- ggarrange(p_r00$cer, 
               p_r00$par, 
               p_r00$common, nrow = 3, ncol = 1)
annotate_figure(p, left = "number of genes up/down in response to each TF deletion")
dev.off()

# random 10
p_r10 <- plotTFBarplots(TFdeldf_cer, TFdeldf_par, 
                        mutate(random_module_genedf10,
                               CCM_color = paste(cer_color, par_color)),
                        .useLiteralColors = FALSE)
pdf("../../aligning_the_molecular_phenotype/paper_figures/TFdel/TFdelBarplots_random10.pdf",
    width = 6, height = 10)
p <- ggarrange(p_r10$cer, 
               p_r10$par, 
               p_r10$common, nrow = 3, ncol = 1)
annotate_figure(p, left = "number of genes up/down in response to each TF deletion")
dev.off()

# random 25
p_r25 <- plotTFBarplots(TFdeldf_cer, TFdeldf_par, 
                        mutate(random_module_genedf25,
                               CCM_color = paste(cer_color, par_color)),
                        .useLiteralColors = FALSE)
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/TFdelBarplots_random25.pdf",
    width = 6, height = 10)
p <- ggarrange(p_r25$cer, 
               p_r25$par, 
               p_r25$common, nrow = 3, ncol = 1)
annotate_figure(p, left = "number of genes up/down in response to each TF deletion")
dev.off()

# random 35
p_r35 <- plotTFBarplots(TFdeldf_cer, TFdeldf_par, 
                        mutate(random_module_genedf35,
                               CCM_color = paste(cer_color, par_color)),
                        .useLiteralColors = FALSE)
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/TFdelBarplots_random35.pdf",
    width = 6, height = 10)
p <- ggarrange(p_r35$cer, 
               p_r35$par, 
               p_r35$common, nrow = 3, ncol = 1)
annotate_figure(p, left = "number of genes up/down in response to each TF deletion")
dev.off()

#### getting to know the TF deletions ####
# which module are all the TFs in?
TFdf |> left_join(y = select(module_genedf10, gene_name, CCM_color),
                  by = "gene_name") |> select(common, gene_name, CCM_color) |>
  filter(common %in% common_TFs) |> 
  print(n = nrow(TFdf)) 
TFdf |> left_join(y = select(module_genedf10, gene_name, CCM_color),
                  by = "gene_name") |> 
  filter(common %in% common_TFs) |> 
  select(CCM_color) |> 
  table() # most are divergent, conserved are mostly in red CCM

# which environments are the TFs DE between species in?
# not the world's prettiest plot, but you can mostly see the TF names
pdf("../../aligning_the_molecular_phenotype/paper_figures/TFdel/TFeffectSizes.pdf", width = 5, height = 4)
ggplot(filter(TFbyexpdf, adjusted_effect_size != 0 &
                common %in% common_TFs), aes(x = experiment, y = adjusted_effect_size)) +
  geom_point(aes(group = common, color = common)) +
  geom_line(aes(group = common, color = common)) +
  geom_text(aes(label = common, color = common),
            check_overlap = FALSE, size = 3,
            position = position_nudge(x = 0.1, y = 0.05)) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_discrete(breaks = c("CC", "HAP4", "LowPi", "TFdelxLowN"),
                   labels = c("Urea Shock", "Saturated Growth", "Low Phosphorus", "Low Nitrogen")) +
  ylab("+ = Higher expressed in cerevisiae\n
       - = Higher expressed in paradoxus") +
  geom_hline(yintercept = 0) +
  ggtitle("Expression divergence of transcription factors in this study")
dev.off()
# effect sizes are fairly consistent between environments. Sometimes not enough divergence
# to be detected, but when it's detected it's always in the same direction
# (i.e. a TF that's higher in cer in one environment is also higher in cer in another)

# investigating case examples (it throws warnings for taking the max of single genes but it's fine)
# DAL80, only DE in LowN and LowPi
gene_idx <- TFdf |> filter(common == "DAL80") |> select(gene_name) |> pull()
plotExpressionProfilePair(.cts1 = counts_all2$cer[,gene_idx, drop = FALSE],
                          .cts2 = counts_all2$par[,gene_idx, drop = FALSE],
                          .info1 = info,
                          .info2 = info, .show_points = TRUE, .normalization = "log2", .method = "line")
# PHD1, consistently DE
gene_idx <- TFdf |> filter(common == "PHD1") |> select(gene_name) |> pull()
plotExpressionProfilePair(.cts1 = counts_all2$cer[,gene_idx, drop = FALSE],
                          .cts2 = counts_all2$par[,gene_idx, drop = FALSE],
                          .info1 = info,
                          .info2 = info, .show_points = TRUE, .normalization = "log2", .method = "line")
# MSN2, never DE
gene_idx <- TFdf |> filter(common == "MSN2") |> select(gene_name) |> pull()
plotExpressionProfilePair(.cts1 = counts_all2$cer[,gene_idx, drop = FALSE],
                          .cts2 = counts_all2$par[,gene_idx, drop = FALSE],
                          .info1 = info,
                          .info2 = info, .show_points = TRUE, .normalization = "log2", .method = "line")

# The lack of divergence in certain environments is sometimes due to low expression (high 0 counts lowering pvalue)
# but sometimes it's just not DE in those environments:
# YAP1, DE in CC and LowPi but expressed in all 4
gene_idx <- TFdf |> filter(common == "YAP1") |> select(gene_name) |> pull()
plotExpressionProfilePair(.cts1 = counts_all2$cer[,gene_idx, drop = FALSE],
                          .cts2 = counts_all2$par[,gene_idx, drop = FALSE],
                          .info1 = info,
                          .info2 = info, .show_points = TRUE, .normalization = "log2", .method = "line")
# AFT1, DE in all but LowPi (maybe due to high 0 counts, but it's also just really not visibly DE like the other 3 environments)
gene_idx <- TFdf |> filter(common == "AFT1") |> select(gene_name) |> pull()
plotExpressionProfilePair(.cts1 = counts_all2$cer[,gene_idx, drop = FALSE],
                          .cts2 = counts_all2$par[,gene_idx, drop = FALSE],
                          .info1 = info,
                          .info2 = info, .show_points = TRUE, .normalization = "log2", .method = "line")

# TODO: is hybrid TF expression a nice neat sum of the parents?
# i.e. if a TF is expressed in one species higher than the other, is hybrid TF
# expression in the middle?
# To assess this we need to plot both parents and the hybrid on one plot, 
# which will require the plotExpressionProfileQuartet function


#### Regulator Divergence ####
# This section will explore this question:
# Does divergence in the expression of a regulator TF explain
# divergence in gene expression?

# We know from the previous section that TFs are either up in cereivisae,
# up in paradoxus, or conserved.
# We can add to this knowledge, knowledge of which genes each TF regulates
# and whether those genes are DE between species
# (and what module those genes are in)
# For example, Do TFs that are up in cerevisiae correlate with genes they 
# activate (genes down in their TF deletion) being up in cerevisiae too?

# For each TF, we need to know which genes it activates (up in WT versus deletion),
# represses (down in WT versus deletion), or doesn't interact with (not DE in WT versus deletion),
# and whether any activation/repression relationships have been rewired between species
# (only up/down/conserved in one species and not the other)

# a glorified table basically showing that there's not much of an 
# enrichment/paucity for any category (except there are no genes that have a 
# conserved repressed relationship with one of the TFs that's up in paradoxus?):
ggplot(regdf, aes(x = TF_relationship, y = direction)) + geom_jitter() +
  theme(axis.text.x = element_text(angle = 90))

# merge 10
# summarising for all 4 experiments, confusing for 1 gene to appear 4 times
table(regdf$gene_name) |> table()
plotdf <- regdf |> left_join(y = select(module_genedf10, gene_name, CCM_color),
                             by = "gene_name") |>
  group_by(gene_name, CCM_color, TF_relationship, direction) |> 
  summarise(mean_lfc = mean(effect_size),
            max_lfc = max(effect_size),
            min_lfc = min(effect_size))
plotdf$greatest_lfc <- apply(plotdf, 1, \(x) {
  maxlfc <- x["max_lfc"] |> as.numeric()
  minlfc <- x["min_lfc"] |> as.numeric()
  if (maxlfc >= abs(minlfc)) {
    return(maxlfc)
  }
  if (maxlfc < abs(minlfc)) {
    return(minlfc)
  }
})
# For each TF, identify which genes go up and down in its deletion and how much
# rewiring there is between species, then plot effect sizes of those genes

# TFs that are upregulated in paradoxus
plot_par <- plotdf |> filter(direction == "up_par" & mean_lfc != 0) |> 
  ggplot(aes(x = TF_relationship, y = mean_lfc)) +
  geom_hline(yintercept = 0, color = "gold2") +
  geom_jitter(aes(color = CCM_color)) +
  scale_color_discrete(limits = plotdf$CCM_color, 
                       type = gsub("none", "grey80", plotdf$CCM_color)) +
  geom_boxplot(outlier.shape = NA, alpha = 0) +
  ggtitle("expression of genes regulated by TFs expressed higher in paradoxus") +
  ylab("Log Fold Change\n<- higher in paradoxus, higher in cerevisiae ->") +
  xlab("") +
  theme_classic() +
  scale_x_discrete(labels = c("activated_noRelationship"="activated\ncerevisiae only",
                              "repressed_noRelationship"="repressed\ncerevisiae only",
                              "repressed_repressed"="repressed\nboth species", 
                              "activated_activated"="activated\nboth species",
                              "noRelationship_activated"="activated\nparadoxus only",
                              "noRelationship_repressed"="repressed\nparadoxus only",
                              "activated_repressed"="activated cerevisiae\nrepressed paradoxus",
                              "repressed_activated"="repressed cerevisiae\nactivated paradoxus"),
                   limits = c("activated_activated",
                              "repressed_repressed",
                              "activated_noRelationship",
                              "repressed_noRelationship",
                              "noRelationship_activated",
                              "noRelationship_repressed",
                              "activated_repressed",
                              "repressed_activated")) +
  theme(axis.text.y = element_blank(), legend.position = "none")
# TFs that are upregulated in cerevisiae
plot_cer <- plotdf |> filter(direction == "up_cer" & mean_lfc != 0) |> 
  ggplot(aes(x = TF_relationship, y = mean_lfc)) +
  geom_hline(yintercept = 0, color = "gold2") +
  geom_jitter(aes(color = CCM_color)) +
  scale_color_discrete(limits = plotdf$CCM_color, 
                       type = gsub("none", "grey80", plotdf$CCM_color)) +
  geom_boxplot(outlier.shape = NA, alpha = 0) +
  ggtitle("expression of genes regulated by TFs expressed higher in cerevisiae") +
  ylab("Log Fold Change\n<- higher in paradoxus, higher in cerevisiae ->") +
  xlab("") +
  theme_classic() +
  scale_x_discrete(labels = c("activated_noRelationship"="activated\ncerevisiae only",
                              "repressed_noRelationship"="repressed\ncerevisiae only",
                              "repressed_repressed"="repressed\nboth species", 
                              "activated_activated"="activated\nboth species",
                              "noRelationship_activated"="activated\nparadoxus only",
                              "noRelationship_repressed"="repressed\nparadoxus only",
                              "activated_repressed"="activated cerevisiae\nrepressed paradoxus",
                              "repressed_activated"="repressed cerevisiae\nactivated paradoxus"),
                   limits = c("activated_activated",
                              "repressed_repressed",
                              "activated_noRelationship",
                              "repressed_noRelationship",
                              "noRelationship_activated",
                              "noRelationship_repressed",
                              "activated_repressed",
                              "repressed_activated")) +
  theme(axis.text.y = element_blank(), legend.position = "none")

# TFs that are conserved
plot_cons <- plotdf |> filter(direction == "conserved" & mean_lfc != 0) |> 
  ggplot(aes(x = TF_relationship, y = mean_lfc)) +
  geom_hline(yintercept = 0, color = "gold2") +
  geom_jitter(aes(color = CCM_color)) +
  scale_color_discrete(limits = plotdf$CCM_color, 
                       type = gsub("none", "grey80", plotdf$CCM_color)) +
  geom_boxplot(outlier.shape = NA, alpha = 0) +
  ggtitle("expression of genes regulated by TFs with conserved expression") +
  ylab("Log Fold Change\n<- higher in paradoxus, higher in cerevisiae ->") +
  xlab("") +
  theme_classic() +
  scale_x_discrete(labels = c("activated_noRelationship"="activated\ncerevisiae only",
                              "repressed_noRelationship"="repressed\ncerevisiae only",
                              "repressed_repressed"="repressed\nboth species", 
                              "activated_activated"="activated\nboth species",
                              "noRelationship_activated"="activated\nparadoxus only",
                              "noRelationship_repressed"="repressed\nparadoxus only",
                              "activated_repressed"="activated cerevisiae\nrepressed paradoxus",
                              "repressed_activated"="repressed cerevisiae\nactivated paradoxus"),
                   limits = c("activated_activated",
                              "repressed_repressed",
                              "activated_noRelationship",
                              "repressed_noRelationship",
                              "noRelationship_activated",
                              "noRelationship_repressed",
                              "activated_repressed",
                              "repressed_activated")) +
  theme(axis.text.y = element_blank(), legend.position = "none")

pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/TF_relationship_effect_sizes10.pdf",
    width = 11, height = 15)
ggarrange(plot_cer, plot_par, plot_cons, nrow = 3, ncol = 1)
dev.off()
# This plot is helpful for me to see module enrichment 
# in each TF relationship group and individual gene effect sizes
# at the same time, but too confusing for the main paper probably

# Example: Brown merge10 is up in cerevisiae and is regulated by a TF that is also higher in cereviaise
brown_idxs <- module_genedf10 |> filter(CCM_color == "brown") |> select(gene_name) |> pull()
regdf |> filter(gene_name %in% brown_idxs & 
                  TF_relationship == "activated_activated" &
                  effect_size > 1) |> print(n = 46)
brownDE_idxs <- regdf |> filter(gene_name %in% brown_idxs & 
                                  TF == "GLN3" & effect_size > 1) |>
  select(gene_name) |> pull() |> unique()
plotExpressionProfileQuartet(counts_TFdel$cer[infos_TFdel$cer$genotype == "WT", brownDE_idxs],
                             counts_TFdel$par[infos_TFdel$par$genotype == "WT", brownDE_idxs],
                             counts_TFdel$cer[infos_TFdel$cer$genotype == "GLN3delete", brownDE_idxs],
                             counts_TFdel$par[infos_TFdel$par$genotype == "GLN3delete", brownDE_idxs],
                             infos_TFdel$cer[infos_TFdel$cer$genotype == "WT",],
                             infos_TFdel$par[infos_TFdel$par$genotype == "WT",],
                             infos_TFdel$cer[infos_TFdel$cer$genotype == "GLN3delete",],
                             infos_TFdel$par[infos_TFdel$par$genotype == "GLN3delete",],
                             "brownDE WT cer",
                             "brownDE GLN3delete cer",
                             "brownDE WT par",
                             "brownDE GLN3delete par",
                             "orange1",
                             "blue2",
                             "orange4",
                             "blue4",
                             .method = "line",
                             .show_points = TRUE,
                             .show_confidence_intervals = TRUE,
                             .normalization = "log2")

#### Identifying TF-regulated gene groups ####

# to generate the TFreg.RData dataset, we tested whether each
# module (all pairwise cer-par color combinations) had an enrichment of
# genes DE due to one TF deletion (so if it appeared statistically likely
# that the module was responding to the same TF in concert)

# now we want to know what proportion of those enriched modules are CCMs

# merge 00
sum(module_genedf00$is_CCM)/nrow(module_genedf00) # 44% of genes are in CCMs
n_DE_CCM_cer <- filter(TFenrichdf00, is_CCM) |> select(nDE_cer) |> pull() |> sum()
n_DE_cer <- sum(TFenrichdf00$nDE_cer)
n_DE_CCM_cer/n_DE_cer # 60 of DE cer genes are in CCMs
n_enriched_CCM_cer <- filter(TFenrichdf00, is_CCM & enriched_cer) |> select(nDE_cer) |> pull() |> sum()
n_enriched_cer <- filter(TFenrichdf00, enriched_cer) |> select(nDE_cer) |> pull() |> sum()
n_enriched_CCM_cer/n_enriched_cer # 96% of enriched DE cer gene groups are in CCMs
n_enriched_cer/n_DE_cer # 30% of DE genes are considered significantly enriched

n_DE_CCM_par <- filter(TFenrichdf00, is_CCM) |> select(nDE_par) |> pull() |> sum()
n_DE_par <- sum(TFenrichdf00$nDE_par) 
n_DE_CCM_par/n_DE_par # 52% of DE par genes are in CCMs
n_enriched_CCM_par <- filter(TFenrichdf00, is_CCM & enriched_par) |> select(nDE_par) |> pull() |> sum()
n_enriched_par <- filter(TFenrichdf00, enriched_par) |> select(nDE_par) |> pull() |> sum()
n_enriched_CCM_par/n_enriched_par # 100% of enriched DE par gene groups are in CCMs
n_enriched_par/n_DE_par # 22% of DE genes are considered significantly enriched

# merge 10
sum(module_genedf10$is_CCM)/nrow(module_genedf10) # 49% of genes are in CCMs
n_DE_CCM_cer <- filter(TFenrichdf10, is_CCM) |> select(nDE_cer) |> pull() |> sum()
n_DE_cer <- sum(TFenrichdf10$nDE_cer) 
n_DE_CCM_cer/n_DE_cer # 62% of DE cer genes are in CCMs
n_enriched_CCM_cer <- filter(TFenrichdf10, is_CCM & enriched_cer) |> select(nDE_cer) |> pull() |> sum()
n_enriched_cer <- filter(TFenrichdf10, enriched_cer) |> select(nDE_cer) |> pull() |> sum()
n_enriched_CCM_cer/n_enriched_cer # 92% of enriched DE cer gene groups are in CCMs
n_enriched_cer/n_DE_cer # 32% of DE genes are considered significantly enriched

n_DE_CCM_par <- filter(TFenrichdf10, is_CCM) |> select(nDE_par) |> pull() |> sum()
n_DE_par <- sum(TFenrichdf10$nDE_par) 
n_DE_CCM_par/n_DE_par # 53% of DE par genes are in CCMs
n_enriched_CCM_par <- filter(TFenrichdf10, is_CCM & enriched_par) |> select(nDE_par) |> pull() |> sum()
n_enriched_par <- filter(TFenrichdf10, enriched_par) |> select(nDE_par) |> pull() |> sum()
n_enriched_CCM_par/n_enriched_par # 100% of enriched DE par gene groups are in CCMs
n_enriched_par/n_DE_par # 22% of DE genes are considered significantly enriched

# merge 25
sum(module_genedf25$is_CCM)/nrow(module_genedf25) # 57% of genes are in CCMs
n_DE_CCM_cer <- filter(TFenrichdf25, is_CCM) |> select(nDE_cer) |> pull() |> sum()
n_DE_cer <- sum(TFenrichdf25$nDE_cer) 
n_DE_CCM_cer/n_DE_cer # 70% of DE cer genes are in CCMs
n_enriched_CCM_cer <- filter(TFenrichdf25, is_CCM & enriched_cer) |> select(nDE_cer) |> pull() |> sum()
n_enriched_cer <- filter(TFenrichdf25, enriched_cer) |> select(nDE_cer) |> pull() |> sum()
n_enriched_CCM_cer/n_enriched_cer # 98% of enriched DE cer gene groups are in CCMs

n_DE_CCM_par <- filter(TFenrichdf25, is_CCM) |> select(nDE_par) |> pull() |> sum()
n_DE_par <- sum(TFenrichdf25$nDE_par) 
n_DE_CCM_par/n_DE_par # 65% of DE par genes are in CCMs
n_enriched_CCM_par <- filter(TFenrichdf25, is_CCM & enriched_par) |> select(nDE_par) |> pull() |> sum()
n_enriched_par <- filter(TFenrichdf25, enriched_par) |> select(nDE_par) |> pull() |> sum()
n_enriched_CCM_par/n_enriched_par # 100% of enriched DE par gene groups are in CCMs

# merge 35
sum(module_genedf35$is_CCM)/nrow(module_genedf35) # 54% of genes are in CCMs
n_DE_CCM_cer <- filter(TFenrichdf35, is_CCM) |> select(nDE_cer) |> pull() |> sum()
n_DE_cer <- sum(TFenrichdf35$nDE_cer) 
n_DE_CCM_cer/n_DE_cer # 63% of DE cer genes are in CCMs
n_enriched_CCM_cer <- filter(TFenrichdf35, is_CCM & enriched_cer) |> select(nDE_cer) |> pull() |> sum()
n_enriched_cer <- filter(TFenrichdf35, enriched_cer) |> select(nDE_cer) |> pull() |> sum()
n_enriched_CCM_cer/n_enriched_cer # 89% of enriched DE cer gene groups are in CCMs

n_DE_CCM_par <- filter(TFenrichdf35, is_CCM) |> select(nDE_par) |> pull() |> sum()
n_DE_par <- sum(TFenrichdf35$nDE_par) 
n_DE_CCM_par/n_DE_par # 62% of DE par genes are in CCMs
n_enriched_CCM_par <- filter(TFenrichdf35, is_CCM & enriched_par) |> select(nDE_par) |> pull() |> sum()
n_enriched_par <- filter(TFenrichdf35, enriched_par) |> select(nDE_par) |> pull() |> sum()
n_enriched_CCM_par/n_enriched_par # 100% of enriched DE par gene groups are in CCMs

# random10 for comparison
n_DE_cer <- sum(random_TFenrichdf10$nDE_cer) 
n_enriched_cer <- filter(random_TFenrichdf10, enriched_cer) |> select(nDE_cer) |> pull() |> sum()
n_enriched_cer/n_DE_cer # 4% of DE cer gene groups are sig enriched

n_DE_par <- sum(random_TFenrichdf10$nDE_par) 
n_enriched_par <- filter(random_TFenrichdf10, enriched_par) |> select(nDE_par) |> pull() |> sum()
n_enriched_par/n_DE_par # 1% of DE par gene groups are sig enriched


############################### Archive ######################################## 
#### troubleshooting DAL80 model ####

# Archived b/c poisson model with no interaction term fixes the problem

# Goal: illustrate how it's treating timepoint as a factor, not 
# the interaction term that causes the DAL80 model to go haywire
# then figure out why that is the case

# first looking at random gene's model
gene_idx <- sample(module_genedf25$gene_name, 1)
deletion <- paste0(common_TFs, "delete") |> sample(1)
# in par
test_genedf <- bind_cols(tibble(expr = counts_TFdel$par[,gene_idx]), infos_TFdel$par) |> 
  filter(genotype %in% c("WT", deletion))
test_genedf$genotype <- as.factor(test_genedf$genotype) |> relevel(ref = "WT")
test_genedf$time_point_str <- as.factor(test_genedf$time_point_str) |> relevel(ref = "0 h, YPD")
spmodtp1 <- glm.nb(expr ~ genotype, data = filter(test_genedf, time_point_num == 0),
                   link = log, init.theta = 7)
spmodtp2 <- glm.nb(expr ~ genotype, data = filter(test_genedf, time_point_num == 60),
                   link = log, init.theta = 7)
spmodtp3 <- glm.nb(expr ~ genotype, data = filter(test_genedf, time_point_num == 960),
                   link = log, init.theta = 7)
sprowtp1 <- summary(spmodtp1)$coefficients[paste0("genotype", deletion),, drop = FALSE]
sprowtp2 <- summary(spmodtp2)$coefficients[paste0("genotype", deletion),, drop = FALSE]
sprowtp3 <- summary(spmodtp3)$coefficients[paste0("genotype", deletion),, drop = FALSE]
sprowtp1 # positive effect size is counts biased for cer, negative is counts biased for par
sprowtp2
sprowtp3
slope <- 1/exp(sprowtp3[1]) # change to different timepoint models to see different slopes/intercepts
intercept <- summary(spmodtp3)$coefficients["(Intercept)",, drop = FALSE][1]

test_m <- glm(expr ~ genotype + time_point_str, family = poisson(link = "log"), data = test_genedf)
mrow <- summary(test_m)$coefficients[paste0("genotype", deletion),, drop = FALSE]
slope <- 1/exp(mrow[1])
plotdf <- test_genedf |>
  select(allele, condition, time_point_num, well_flask_ID, expr, genotype) %>%
  mutate(non_unique_sample_name = paste(condition, well_flask_ID)) %>% 
  pivot_wider(id_cols = c("time_point_num", "well_flask_ID"), 
              names_from = genotype, 
              values_from = expr) |> 
  select(all_of(c("time_point_num", "well_flask_ID", deletion, "WT"))) |> 
  drop_na()
colnames(plotdf) <- gsub(deletion, "deletion", colnames(plotdf))
max_expr <- max(select(plotdf, WT, deletion), na.rm = TRUE)
ggplot(plotdf, aes(x = deletion, y = WT)) + 
  geom_point(aes(shape = factor(time_point_num))) + 
  geom_abline(color = "gold", slope = slope, intercept = intercept) + 
  geom_abline(color = "midnightblue", slope = 1, intercept = 0) + 
  xlim(c(0, max_expr)) + ylim(c(0, max_expr)) + 
  theme_classic() + ggtitle("Standard GLM") +
  xlab(paste("Expression in", deletion)) +
  ylab("Expression in WT") +
  theme(legend.title = element_blank()) +
  scale_color_discrete(labels = c("data", "fitted values"))

# Visualizing individual TF deletion models
library(MASS, include.only = "glm.nb")
visualizeTFdelModel <- function(.gene_idx, .genotype, .species) {
  if (.species == "cer") {
    cts <- counts_TFdel$cer
    info_tf <- infos_TFdel$cer
  }
  if (.species == "par") {
    cts <- counts_TFdel$par
    info_tf <- infos_TFdel$par
  }
  gdf <- bind_cols(tibble(expr = cts[,.gene_idx]), info_tf)
  gdf$genotype <- as.factor(gdf$genotype) |> relevel(ref = "WT")
  #gdf <- filter(gdf, genotype %in% c("WT", .genotype))
  gdf$time_point_str <- as.factor(gdf$time_point_str) |> relevel(ref = "0 h, YPD")
  # m <- glm.nb(expr ~ genotype + time_point_str, data = gdf, link = log)
  # m <- glm.nb(expr ~ genotype * time_point_num, data = gdf, link = log)
  # m <- glm(expr ~ genotype * time_point_str, data = gdf, family = poisson(link = "log"))
  m <- glm(expr ~ genotype + time_point_str, data = gdf, family = poisson(link = "log"))
  gdf$p_hat <- m$fitted.values
  coefsum <- summary(m)$coefficients
  coeffs <- coefsum[grepl("^genotype", rownames(coefsum)), "Estimate", drop = FALSE]
  pvals <- coefsum[,4][grepl("^genotype", rownames(coefsum))]
  sigGenotypes <- rownames(coeffs)[(pvals < 1e-5) & (abs(coeffs) > 1)] %>% gsub(pattern = "genotype", replacement = "")
  gdf$sig <- gdf$genotype %in% sigGenotypes
  gdf$sig[gdf$genotype == "WT"] <- NA
  p1 <- ggplot(data = gdf, aes(x = genotype, y = expr)) + 
    geom_point(aes(color = sig)) +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(paste(.gene_idx, "expression in all genotypes"))
  sprow <- summary(m)$coefficients[paste0("genotype", .genotype),, drop = FALSE]
  cat(sprow, "\n")
  slope <- 1/exp(sprow[1])
  intercept <- summary(m)$coefficients["(Intercept)",, drop = FALSE][1]
  plotdf <- gdf %>% mutate(fittedvals = m$fitted.values)%>% select(allele, experiment, time_point_str, well_flask_ID, expr, fittedvals, genotype) %>% 
    pivot_longer(cols = c(expr, fittedvals)) %>% 
    filter(genotype %in% c("WT", .genotype)) |> 
    pivot_wider(id_cols = c("experiment", "time_point_str", "well_flask_ID", "name"),
                names_from = "genotype", values_from = "value")
  colnames(plotdf) <- gsub(.genotype, "deletion", colnames(plotdf))
  max_expr <- max(select(plotdf, WT, deletion), na.rm = TRUE)
  p2 <- ggplot(plotdf, aes(y = WT, x = deletion)) + 
    geom_point(aes(color = name)) + geom_abline(color = "gold", slope = slope, intercept = intercept) + 
    geom_abline(color = "midnightblue", slope = 1, intercept = 0) + 
    xlim(c(0, max_expr)) + ylim(c(0, max_expr)) + 
    theme_classic() + ggtitle("Standard GLM") +
    xlab("Expression in deletion") +
    ylab(paste("Expression in WT")) +
    theme(legend.title = element_blank()) +
    scale_color_discrete(labels = c("data", "fitted values")) +
    annotate("text", x = max_expr/2, y = 10, label = round(slope, 2), size = 2.5)
  return(list("all_TFs" = p1, "TF_of_interest" = p2))
}
# test <- visualizeTFdelModel("YGL071W", "AFT1delete") 
test <- visualizeTFdelModel("YKR034W", "DAL80delete", .species = "par") # "YKR034W" = "DAL80delete"
ggarrange(plotlist = test, nrow = 1, ncol = 2)

test <- visualizeTFdelModel("YKR034W", "STB5delete", .species = "par") # "YKR034W" = "DAL80delete"
ggarrange(plotlist = test, nrow = 1, ncol = 2)
# We only have STB5delete samples from first two timepoints, so they can only be compared to
# WT from first 2 timepoints. That's why STB5 is so low but not considered significant 
# (rightfully so, as DAL80 only appears to be expressed at TP3, so there's no telling 
# if STB5 deletion affects that or not)

test <- visualizeTFdelModel("YIR017C", "MET28delete", .species = "par")
ggarrange(plotlist = test, nrow = 1, ncol = 2)
# ranking TF deletions by pleitropic-ness 
# which pairs of TFdels are most likely to be affected together?
multiple_cer <- which(rowSums(TFdel_pval_cer[,-1] < 1e-5) > 1)
TFdel_groups_cer <- apply(TFdel_pval_cer[multiple_cer, -1], 1, \(x) {
  reduce(.x = names(x)[which(x < 1e-5)], .f = paste)
}) |> table()
TFdel_groups_cer[TFdel_groups_cer > 1] |> sort(decreasing = TRUE)

multiple_par <- which(rowSums(TFdel_pval_par[,-1] < 1e-5) > 1)
TFdel_groups_par <- apply(TFdel_pval_par[multiple_par, -1], 1, \(x) {
  reduce(.x = names(x)[which(x < 1e-5)], .f = paste)
}) |> table()
TFdel_groups_par[TFdel_groups_par > 1] |> sort(decreasing = TRUE)

# for cer: CBF1, BAS1, HAP3, AFT1, GCR2, STB5 combinations in top groups (MET28 is nowhere to be found, even though it complexes with CBF1...)
# interestingly for par: MET28 is in the top, then there are the rest of the ones in cer: CBF1, BAS1, GCR2, STB5, plus MET28

#### Extending subset expression divergence to all 4 environments ####
# Section is vague, just thinking that it might be good to use the F1 hybrid
# data in other experiments, even though we don't have data on the regulatory TF,
# to see if there's any pattern to how the subsets behave across environments.
# Are they more environmentally responsive? Maybe in the species that has regulation
# versus the one that doesn't?

# all of those case examples above just happened to have strongly divergent expression between species
# (probably because we picked examples based on effect size)
# but is it possibly a common pattern? Let's just compare effect size of genes
# that are regulated by ANY TF in EITHER (or both) species versus genes that aren't
# first getting list of genes that are present in any TF's regulated gene group
TFreg_idxs_cer <- TFenrich_genedf |> filter(enriched_cer & DE_cer) |>
  select(gene_name) |> pull()
TFreg_idxs_par <- TFenrich_genedf |> filter(enriched_par & DE_par) |>
  select(gene_name) |> pull()
TFreg_idxs <- c(TFreg_idxs_cer, TFreg_idxs_par) |> unique()
plotdf <- spaldf |> filter(coefficient == "species" & experiment == "all") |> 
  select(gene_name, effect_size, pvalue) |> 
  mutate(TFreg = gene_name %in% TFreg_idxs)
sum(plotdf$TFreg)
sum(!plotdf$TFreg)
plotdf$adj_effect_size <- if_else(plotdf$pvalue < p_thresh,
                                  true = plotdf$effect_size,
                                  false = 0)
ggplot(plotdf, aes(x = TFreg, y = adj_effect_size)) + geom_jitter()
# nothing obvious, and we just care about obvious

# TODO: figure out a systematic way of quantifying the observations
# for the above 4 examples for all TF-regulated gene groups

# possible strategy: for each group of TF-regulated genes in one CCM, 
# compare its expression between species and F1 hybrid alleles 
# (by avg expression correlation, or maybe effect size from single gene models?) to determine
# A) whether its strongly diverged in expression 
# (should be more diverged for gene groups that aren't conserved between species)
# B) what proportion of its divergence is in cis or trans

# first do these steps for case examples above then systematically:

# Step 1: identify TF regulated gene groups (probably limit to groups of 10 or more genes or >=5% of CCM or something)
# Step 1.5: also identify genes never regulated by any TF?
# TODO: we're here, this df enumerates all TF-regulated gene groups:
TFenrich_genedf
# Step 2: get average correlation between species for those TF-regulated gene groups
# between parents and hybrid alleles
# Step 3: Plot those avg correlations in parents vs hybrids colored by CCM
# Step 4: Quantify how many gene groups persist in hybrid (rewired) versus
# don't persist in hybrid (trans-divergent)



# And if that looks confusing, just stick with bar plots and pull out case examples!

# And then go back and remove TF deletion data from network construction
# AND THEN YOURE DONE


#### enrichment exploration ####
# # are most DE genes clustered in modules together or no?
# TFenrichdf25$nDE_cer |> table()
# TFenrichdf25$nDE_par |> table() # looks 50-50
# 
# # simple enrichment test: are more DE genes present in gene groups of <5 or >=5?
# table_counts_cer <- table(TFenrichdf25$nDE_cer) |> as.numeric()
# table_names_cer <- table(TFenrichdf25$nDE_cer) |> names() |> as.numeric()
# # nDE in gene groups <5 in cer:
# sum(table_counts_cer[table_names_cer < 5]*table_names_cer[table_names_cer < 5])
# # nDE in gene groups >= 5 in cer:
# sum(table_counts_cer[table_names_cer >= 5]*table_names_cer[table_names_cer >= 5])
# # many more in small groups in cer
# # repeat for par
# table_counts_par <- table(TFenrichdf25$nDE_par) |> as.numeric()
# table_names_par <- table(TFenrichdf25$nDE_par) |> names() |> as.numeric()
# # nDE in gene groups <5 in cer:
# sum(table_counts_par[table_names_par < 5]*table_names_par[table_names_par < 5])
# # nDE in gene groups >= 5 in cer:
# sum(table_counts_par[table_names_par >= 5]*table_names_par[table_names_par >= 5])
# 
# # total number of enriched groups in each species
# sum(TFenrichdf25$enriched_cer)
# sum(TFenrichdf25$enriched_par)
# 
# # are DE genes more likely to be part of enriched gene groups?
# TFenrichdf25 |> group_by(enriched_cer) |> summarise(nDE = sum(nDE_cer))
# TFenrichdf25 |> group_by(enriched_par) |> summarise(nDE = sum(nDE_par))
# # strongly no, but it obviously depends on our significant enrichment criteria
# # what if you filter for CCMs?
# TFenrichdf25 |> filter(is_CCM) |> group_by(enriched_cer) |> summarise(nDE = sum(nDE_cer))
# TFenrichdf25 |> filter(is_CCM) |> group_by(enriched_par) |> summarise(nDE = sum(nDE_par))
# # oh dang yeah. Vast majority of DEs are in CCMs, at least for enriched groups
# TFenrichdf25 |> filter(enriched_cer) |>  group_by(is_CCM) |> 
#   summarise(nDE_cer_all = sum(nDE_cer),
#             n = n())
# TFenrichdf25 |> filter(enriched_par) |>  group_by(is_CCM) |> 
#   summarise(nDE_par_all = sum(nDE_par),
#             n = n())
# sum(module_genedf25$is_CCM)/nrow(module_genedf25) # 30% of genes are in CCMs
# (filter(TFenrichdf25, is_CCM) |> select(nDE_cer) |> pull() |> sum())/sum(TFenrichdf25$nDE_cer) # 36% of DE cer genes are in CCMs
# (filter(TFenrichdf25, is_CCM & enriched_cer) |> select(nDE_cer) |> pull() |> sum())/
#   (filter(TFenrichdf25, enriched_cer) |> select(nDE_cer) |> pull() |> sum()) # 59% of enriched DE cer gene groups are in CCMs
# # is this just because there are bigger blocks in the CCMs?
# plotdf <- module_genedf25 |> select(cer_color, par_color, is_CCM) |> group_by(cer_color, par_color, is_CCM) |> summarise(n = n())
# ggplot(plotdf, aes(x = is_CCM, y = n)) + geom_jitter(color = plotdf$cer_color) +
#   scale_color_discrete(limits = plotdf$cer_color, type = plotdf$cer_color)
# # CCMs definitely have a larger mean block size, but there are plenty of large blocks that aren't CCMs
# # If we drop small blocks, blocks smaller than 50 say, do we still see CCM enrichment?
# TFenrichdf25 |> filter(block_size > 50 & enriched_cer) |>  group_by(is_CCM) |> 
#   summarise(nDE_cer_all = sum(nDE_cer),
#             n = n()) 
# TFenrichdf25 |> filter(block_size > 50 & enriched_par) |>  group_by(is_CCM) |> 
#   summarise(nDE_par_all = sum(nDE_par),
#             n = n())
# # the enrichment remains. It's not just due to larger block size of CCMs
# # is the opposite true? Are un-enriched DE genes (we'll call them singletons, even though there might be like 3 or so)
# # appearing more in nonCCMs?
# TFenrichdf25 |> filter(!enriched_cer) |> group_by(is_CCM) |> 
#   summarise(nDE_cer_all = sum(nDE_cer),
#             n = n())
# TFenrichdf25 |> filter(!enriched_par) |> group_by(is_CCM) |> 
#   summarise(nDE_par_all = sum(nDE_par),
#             n = n()) # not really considering that more of these groups are not CCMs
# 
# # summary: When a TF affects genes in a certain module (cer/par color combination),
# # it tends to affect slightly more genes when that module is a CCM:
# 
# # plotting nDE for CCMs vs non CCMs where each point is one cer/par color combination
# # cer
# ggplot(filter(TFenrichdf25, block_size > 10 & nDE_cer > 0), aes(x = is_CCM, y = nDE_cer/block_size)) + 
#   geom_jitter(aes(color = deletion))
# ggplot(filter(TFenrichdf25, enriched_cer), aes(x = is_CCM, y = nDE_cer/block_size)) + 
#   #geom_jitter(aes(color = deletion), position = position_jitter(seed = 2)) +
#   geom_text(aes(label = deletion, color = deletion), position = position_jitter(seed = 2)) +
#   theme_classic() +
#   theme(legend.position = "none")
# # par
# ggplot(filter(TFenrichdf25, block_size > 10 & nDE_par > 0), aes(x = is_CCM, y = nDE_par/block_size)) + geom_jitter(aes(color = deletion))
# ggplot(filter(TFenrichdf25, enriched_par), aes(x = is_CCM, y = nDE_par/block_size)) + geom_jitter(aes(color = deletion))
# 

#### case examples for modules with divergent TF behavior - pre WT-only networks ####
### 1) Magenta genes repressed by GLN3 in paradoxus are up in cerevisiae. GLN3 is expressed
#    higher in cerevisiae
# possible interpretation: they lost their ability to be repressed by GLN3 in cerevisiae,
# leading to higher expression

# # first which genes are these and how much of a magenta subset are they?
# magenta_DE_idxs <- regdf |> filter(TF == "GLN3" & CCM_color == "magenta" & effect_size > 0) |> 
#   select(gene_name) |> pull() |> unique()
# magenta_idxs <- module_genedf |> filter(CCM_color == "magenta") |> select(gene_name) |> pull() |> unique()
# sum(magenta_DE_idxs %in% magenta_idxs)/length(magenta_idxs)
# sum(magenta_DE_idxs %in% magenta_idxs)/length(magenta_DE_idxs)
# # What do they look like in WT cerevisiae?
# plotExpressionProfilePair(.cts1 = counts_all2$cer[,setdiff(magenta_idxs, magenta_DE_idxs)],
#                           .cts2 = counts_all2$cer[,magenta_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "S. cerevisiae NonDE",
#                           .name2 = "S. cerevisiae DE",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2")
# # They are an insanely obvious subset of higher expressed magenta genes
# # repeat for paradoxus
# plotExpressionProfilePair(.cts1 = counts_all2$par[,setdiff(magenta_idxs, magenta_DE_idxs)],
#                           .cts2 = counts_all2$par[,magenta_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "S. paradoxus NonDE",
#                           .name2 = "S. paradoxus DE",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2")
# # not really as obvious, the genes are a little higher in LowN and LowPi
# # if that change is due to trans factors, both alleles for the DE_idxs should be more similar in the F1 hybrid
# # F1 hybrid
# plotExpressionProfilePair(.cts1 = counts_all2_allele$cer[, magenta_DE_idxs],
#                           .cts2 = counts_all2_allele$par[, magenta_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "S. cerevisiae F1 hybrid allele",
#                           .name2 = "S. paradoxus F1 hybrid allele",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2")
# # parent
# plotExpressionProfilePair(.cts1 = counts_all2$cer[, magenta_DE_idxs],
#                           .cts2 = counts_all2$par[, magenta_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "S. cerevisiae",
#                           .name2 = "S. paradoxus",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2")
# 
# ### 2) Midnightblue genes activated by PHO4 in cerevisiae are up in cerevisiae. PHO4
# #    is expressed higher in paradoxus, 
# # possible interpretation: they lost their ability to be activated by PHO4 in paradoxus,
# # leading to lower expression
# 
# # first which genes are these and how much of a module subset are they?
# midnightblue_DE_idxs <- regdf |> filter(TF == "PHO4" & CCM_color == "midnightblue" & effect_size > 0) |> 
#   select(gene_name) |> pull() |> unique()
# midnightblue_idxs <- module_genedf |> filter(CCM_color == "midnightblue") |> select(gene_name) |> pull() |> unique()
# sum(midnightblue_DE_idxs %in% midnightblue_idxs)/length(midnightblue_idxs)
# sum(midnightblue_DE_idxs %in% midnightblue_idxs)/length(midnightblue_DE_idxs)
# # What do they look like in WT cerevisiae?
# plotExpressionProfilePair(.cts1 = counts_all2$cer[,setdiff(midnightblue_idxs, midnightblue_DE_idxs)],
#                           .cts2 = counts_all2$cer[,midnightblue_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "S. cerevisiae NonDE",
#                           .name2 = "S. cerevisiae DE",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2") + ggtitle("midnightblue module, WT cerevisiae parent")
# # repeat for paradoxus
# plotExpressionProfilePair(.cts1 = counts_all2$par[,setdiff(midnightblue_idxs, midnightblue_DE_idxs)],
#                           .cts2 = counts_all2$par[,midnightblue_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "S. paradoxus NonDE",
#                           .name2 = "S. paradoxus DE",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2") + ggtitle("midnightblue module, WT paradoxus parent")
# # upregulation in LowPi is still present, but otherwise not DE
# # if that change is due to trans factors, both alleles for the DE_idxs should be more similar in the F1 hybrid
# # F1 hybrid
# plotExpressionProfilePair(.cts1 = counts_all2_allele$cer[, midnightblue_DE_idxs],
#                           .cts2 = counts_all2_allele$par[, midnightblue_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "S. cerevisiae F1 hybrid allele",
#                           .name2 = "S. paradoxus F1 hybrid allele",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2")
# # parent
# plotExpressionProfilePair(.cts1 = counts_all2$cer[, midnightblue_DE_idxs],
#                           .cts2 = counts_all2$par[, midnightblue_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "S. cerevisiae",
#                           .name2 = "S. paradoxus",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2")
# ### 3) Pink genes activated by HAP4 in paradoxus are up in paradoxus. HAP4 is also
# # up in paradoxus
# # possible interpretation: HAP4 activates these genes and its increased expression
# # leads to increased expression of the genes
# 
# pink_DE_idxs <- regdf |> filter(TF == "HAP4" & CCM_color == "pink" & effect_size < 0) |> 
#   select(gene_name) |> pull() |> unique()
# pink_idxs <- module_genedf |> filter(CCM_color == "pink") |> select(gene_name) |> pull() |> unique()
# sum(pink_DE_idxs %in% pink_idxs)/length(pink_idxs)
# sum(pink_DE_idxs %in% pink_idxs)/length(pink_DE_idxs)
# # What do they look like in WT cerevisiae?
# plotExpressionProfilePair(.cts1 = counts_all2$cer[,setdiff(pink_idxs, pink_DE_idxs)],
#                           .cts2 = counts_all2$cer[,pink_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "S. cerevisiae NonDE",
#                           .name2 = "S. cerevisiae DE",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2")
# # repeat for paradoxus
# plotExpressionProfilePair(.cts1 = counts_all2$par[,setdiff(pink_idxs, pink_DE_idxs)],
#                           .cts2 = counts_all2$par[,pink_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "S. paradoxus NonDE",
#                           .name2 = "S. paradoxus DE",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2")
# # in both cases they're up, but more in paradoxus
# # F1 hybrid
# plotExpressionProfilePair(.cts1 = counts_all2_allele$cer[, pink_DE_idxs],
#                           .cts2 = counts_all2_allele$par[, pink_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "S. cerevisiae F1 hybrid allele",
#                           .name2 = "S. paradoxus F1 hybrid allele",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2")
# # parent
# plotExpressionProfilePair(.cts1 = counts_all2$cer[, pink_DE_idxs],
#                           .cts2 = counts_all2$par[, pink_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "S. cerevisiae",
#                           .name2 = "S. paradoxus",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2")
# 
# # 4) 2 Yellow genes (YIL169C aka CSS1 and YKL070W, both unknown function) are strongly
# # up in paradoxus. CSS1 is activated by PHO4 (up in paradoxus) and RPN4 (conserved),
# # YKL070W is repressed by MET28 (conserved), both of these are only relationships in paradoxus
# # possible explanation: These are two very divergent genes, possibly only really expressed
# # in paradoxus
# 
# yellow_DE_idxs <- c("YIL169C", "YKL070W")
# yellow_idxs <- module_genedf |> filter(CCM_color == "yellow") |> select(gene_name) |> pull() |> unique()
# sum(yellow_DE_idxs %in% yellow_idxs)/length(yellow_idxs)
# sum(yellow_DE_idxs %in% yellow_idxs)/length(yellow_DE_idxs)
# # What do they look like in WT cerevisiae?
# plotExpressionProfilePair(.cts1 = counts_all2$cer[,setdiff(yellow_idxs, yellow_DE_idxs)],
#                           .cts2 = counts_all2$cer[,yellow_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "NonDE",
#                           .name2 = "DE",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2") + ggtitle("yellow module w/ 2 DE genes, WT cerevisiae parent")
# # repeat for paradoxus
# plotExpressionProfilePair(.cts1 = counts_all2$par[,setdiff(yellow_idxs, yellow_DE_idxs)],
#                           .cts2 = counts_all2$par[,yellow_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "NonDE",
#                           .name2 = "DE",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2") + ggtitle("yellow module w/ 2 DE genes, WT paradoxus parent")
# # F1 hybrid
# plotExpressionProfilePair(.cts1 = counts_all2_allele$cer[, yellow_DE_idxs],
#                           .cts2 = counts_all2_allele$par[, yellow_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "S. cerevisiae F1 hybrid allele",
#                           .name2 = "S. paradoxus F1 hybrid allele",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2") + ggtitle("2 DE yellow genes, F1 hybrid alleles")
# # parent
# plotExpressionProfilePair(.cts1 = counts_all2$cer[, yellow_DE_idxs],
#                           .cts2 = counts_all2$par[, yellow_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "S. cerevisiae",
#                           .name2 = "S. paradoxus",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2") + ggtitle("2 DE yellow genes, parental expression")
# # Well that's very clearly the result of cis regulation

#### Totally useless Fisher exact threshold calculator ####
# for some reason I thought I would be doing multiple fisher exact tests per TF/module,
# making this a reasonable thing to have coded, but then I realized I'll only do one,
# making this less computationally efficient

# # given marginal totals, returns the minimum number of genes in a given module
# # that need to be DE upon a give TF's deletion to consider that module enriched
# # for genes that are regulated by that TF
# # marginal totals needed:
# # nGenes, total number of genes
# # nDE, the total number of genes that are DE upon a given TF's deletion (DE up or down)
# # modsize, number of genes in the module
# # pthresh, desired pvalue cutoff for significance
# # @output: number of DE genes in that module required to qualify as significant enrichment
# calculateMinDEforEnrichment <- function(.nDE, .modsize, .nGenes, .pthresh) {
#   # helper for calculateMinDEforEnrichment
#   getFisherExactPvalue <- function(.a, .nDE, .modsize, .nGenes) {
#     contingencytab <- rbind(c(.a, .nDE - .a),
#                             c(.modsize - .a, .nGenes - .modsize - .nDE + .a))
#     return(fisher.test(contingencytab, alternative = "greater")$p.value)
#   }
#   findFisherThreshold <- function(.a, .nDE, .modsize, .nGenes, .pthresh) {
#     p <- getFisherExactPvalue(.a = .a, .nDE = .nDE, .modsize = .modsize, .nGenes = .nGenes)
#     if (p > .pthresh) {
#       return(.a)
#     }
#     else {
#       return(findFisherThreshold(.a = (.a - 1), .nDE = .nDE, .modsize = .modsize, .nGenes = .nGenes, .pthresh = .pthresh))
#     }
#   }
#   a <- findFisherThreshold(.a = .nDE, .nDE = .nDE, .modsize = .modsize, .nGenes = .nGenes, .pthresh = .pthresh)
#   return(a + 1) # a is the first value that isn't significant, starting at nDE and working down to 0, so a+1 is the lease extreme value still consider significant, aka the minDE required for enrichment
# }
# # tests for calculateMinDEforEnrichment
# # ACE2 min thresh for enrichment in the 
# testminDE <- calculateMinDEforEnrichment(.nDE = 12, .modsize = 67, .nGenes = 4863, .pthresh = p_thresh)
# fisher.test(rbind(c(testminDE, 12 - testminDE),
#                   c(67 - testminDE, 4863 - 67 - 12 + testminDE)),
#             alternative = "greater")$p.value < p_thresh


#### Probably Archive: Are genes that have lost TF regulation in one species more constitutively expressed in that species versus the other? ####
# # That seems to be the pattern based on these 4 examples (well first 3 really, hard to say with yellow cause there's only 2 of them)
# 
# # TODO: Case example with midnight blue
# # lost reg in cer
# midnightblue_DE_idxs
# 
# # very quick plot of how much effect size changes between the 4 experiments for genes that lost regulation in one species versus the other
# lost_reg_in_cer_idxs <- regdf |> filter(TF_relationship %in% c("noRelationship_activated",
#                                                                "noRelationship_repressed")) |> 
#   select(gene_name) |> pull() |> unique()
# lost_reg_in_par_idxs <- regdf |> filter(TF_relationship %in% c("activated_noRelationship",
#                                                                "repressed_noRelationship")) |> 
#   select(gene_name) |> pull() |> unique()
# conserved_reg_idxs <- regdf |> filter(TF_relationship %in% c("activated_activated",
#                                                              "repressed_repressed")) |> 
#   select(gene_name) |> pull() |> unique()
# 
# # TODO: for each gene in each category, get its expression variation across all 4 environments
# # in each species and see if there's a difference in expression variance between species
# vardf <- bind_rows(tibble(gene_name = lost_reg_in_cer_idxs, type = "lost_in_cer"),
#                    tibble(gene_name = lost_reg_in_par_idxs, type = "lost_in_par"),
#                    tibble(gene_name = conserved_reg_idxs, type = "conserved"))
# vardf$expr_var_cer <- map(vardf$gene_name, \(g) {
#   cts_cer <- counts_all2$cer[,g]
#   return(var(cts_cer))
# }) |> unlist()
# vardf$expr_var_par <- map(vardf$gene_name, \(g) {
#   cts_par <- counts_all2$par[,g]
#   return(var(cts_par))
# }) |> unlist()
# plotdf <- vardf |> pivot_longer(cols = c("expr_var_cer", "expr_var_par"), 
#                                 names_to = "species",
#                                 values_to = "expr_var",
#                                 names_prefix = "expr_var_")
# # lost in cer, should have lower variance in cerevisiae
# ggplot(filter(plotdf, type == "lost_in_cer"), aes(x = species, y = log(expr_var + 1))) + 
#   geom_line(aes(group = gene_name)) + 
#   geom_point(data = filter(plotdf, type == "lost_in_cer") |> group_by(species) |> summarise(median_expr_var = median(expr_var)),
#              aes(x = species, y = log(median_expr_var + 1)), color = "red")
# # lost in par, should have lower variance in paradoxus
# ggplot(filter(plotdf, type == "lost_in_par"), aes(x = species, y = log(expr_var + 1))) + 
#   geom_line(aes(group = gene_name)) + 
#   geom_point(data = filter(plotdf, type == "lost_in_par") |> group_by(species) |> summarise(median_expr_var = median(expr_var)),
#              aes(x = species, y = log(median_expr_var + 1)), color = "red")
# # conserved regulation, should have same variance
# ggplot(filter(plotdf, type == "conserved"), aes(x = species, y = log(expr_var + 1))) + 
#   geom_line(aes(group = gene_name)) +  
#   geom_point(data = filter(plotdf, type == "conserved") |> group_by(species) |> summarise(median_expr_var = median(expr_var)),
#              aes(x = species, y = log(median_expr_var + 1)), color = "red")
# 
# # ummm no not really


#### Bonus examples: any enriched gene group versus randomly selected module subset of the same size ####
# # repeat this section to compare different enriched groups
# 
# # this could be adaptedto work with the fact that we no longer have TFenrich genedfs
# # but we already have a hardcoded version above, and I don't think being able to make it so randomized is worth it
# 
# coinFlip <- sample(c("heads", "tails"), 1)
# if (coinFlip == "heads") { # Heads we do an enriched group in cerevisiae
#   anyEnriched_row <- TFenrichdf10 |> filter(enriched_cer) |> slice_sample(n = 1)
#   anyEnriched_DE_idxs <- TFenrich_genedf10 |> filter(cer_color == as.character(anyEnriched_row$cer_color) &
#                                                        par_color == as.character(anyEnriched_row$par_color) &
#                                                        deletion == as.character(anyEnriched_row$deletion) &
#                                                        DE_cer) |> 
#     select(gene_name) |> pull()
# }
# if (coinFlip == "tails") { # Tails we do an enriched group in paradoxus
#   anyEnriched_row <- TFenrichdf10 |> filter(enriched_par) |> slice_sample(n = 1)
#   anyEnriched_DE_idxs <- TFenrich_genedf10 |> filter(cer_color == as.character(anyEnriched_row$cer_color) &
#                                                        par_color == as.character(anyEnriched_row$par_color) &
#                                                        deletion == as.character(anyEnriched_row$deletion) &
#                                                        DE_par) |> 
#     select(gene_name) |> pull()
# }
# restOfModule_idxs <- module_genedf10 |> filter(cer_color == as.character(anyEnriched_row$cer_color) &
#                                                  par_color == as.character(anyEnriched_row$par_color)) |> 
#   select(gene_name) |> pull() |> setdiff(y = anyEnriched_DE_idxs)
# randomSubset_idxs <- sample(restOfModule_idxs, length(anyEnriched_DE_idxs), replace = FALSE)
# # plotting
# # WT cerevisiae, enriched versus same-size subest of not-enriched from same module
# plotExpressionProfilePair(.cts1 = counts_all2$cer[,randomSubset_idxs],
#                           .cts2 = counts_all2$cer[,anyEnriched_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "not enriched subset of same size",
#                           .name2 = "enriched",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2") + ggtitle(paste(anyEnriched_row$deletion, 
#                                                                    anyEnriched_row$cer_color, 
#                                                                    anyEnriched_row$par_color, 
#                                                                    "module\nWT cerevisiae parent"))
# # repeat for paradoxus
# plotExpressionProfilePair(.cts1 = counts_all2$par[,randomSubset_idxs],
#                           .cts2 = counts_all2$par[,anyEnriched_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "not enriched subset of same size",
#                           .name2 = "enriched",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2") + ggtitle(paste(anyEnriched_row$deletion, 
#                                                                    anyEnriched_row$cer_color, 
#                                                                    anyEnriched_row$par_color, 
#                                                                    "module\nWT paradoxus parent"))
# # comparing enriched genes between species
# # F1 hybrid
# plotExpressionProfilePair(.cts1 = counts_all2_allele$cer[, anyEnriched_DE_idxs],
#                           .cts2 = counts_all2_allele$par[, anyEnriched_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "enriched, S. cerevisiae F1 hybrid allele",
#                           .name2 = "enriched, S. paradoxus F1 hybrid allele",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2") + ggtitle(paste(anyEnriched_row$deletion, 
#                                                                    anyEnriched_row$cer_color, 
#                                                                    anyEnriched_row$par_color, 
#                                                                    "module\nF1 hybrid"))
# # parent
# plotExpressionProfilePair(.cts1 = counts_all2$cer[, anyEnriched_DE_idxs],
#                           .cts2 = counts_all2$par[, anyEnriched_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "enriched, S. cerevisiae parent",
#                           .name2 = "enriched, S. paradoxus parent",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2") + ggtitle(paste(anyEnriched_row$deletion, 
#                                                                    anyEnriched_row$cer_color, 
#                                                                    anyEnriched_row$par_color, 
#                                                                    "module\nparents"))
# # how do the unenriched genes compare?
# # F1 hybrid
# plotExpressionProfilePair(.cts1 = counts_all2_allele$cer[, randomSubset_idxs],
#                           .cts2 = counts_all2_allele$par[, randomSubset_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "not enriched, S. cerevisiae F1 hybrid allele",
#                           .name2 = "not enriched, S. paradoxus F1 hybrid allele",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2") + ggtitle(paste(anyEnriched_row$deletion, 
#                                                                    anyEnriched_row$cer_color, 
#                                                                    anyEnriched_row$par_color, 
#                                                                    "module\nF1 hybrid"))
# # parent
# plotExpressionProfilePair(.cts1 = counts_all2$cer[, randomSubset_idxs],
#                           .cts2 = counts_all2$par[, randomSubset_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "not enriched, S. cerevisiae parent",
#                           .name2 = "not enriched, S. paradoxus parent",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2") + ggtitle(paste(anyEnriched_row$deletion, 
#                                                                    anyEnriched_row$cer_color, 
#                                                                    anyEnriched_row$par_color, 
#                                                                    "module\nparents"))
# # now how about versus the set of genes DE by that TF regardless of module membership. Aka is 
# # limiting to co-regulation groups more informative of what's going on?
# regdf |> filter(TF == as.character(anyEnriched_row$deletion) &
#                   gene_name %in% anyEnriched_DE_idxs) |> select(TF_relationship, gene_name) |> table()
# # genes might vary in whether the TF relationship is conserved or not, but every gene should be DE in the focal species' (and probably in the same direction)
# TF_relationships_in_enrichment <-  filter(regdf, TF == as.character(anyEnriched_row$deletion) &
#                                             gene_name %in% anyEnriched_DE_idxs) |> select(TF_relationship) |> pull() |> unique() 
# regdf |> filter(TF == as.character(anyEnriched_row$deletion) & TF_relationship %in% TF_relationships_in_enrichment) |> select(TF, TF_relationship) |> table()
# all_TFreg_idxs <- regdf |> filter(TF == as.character(anyEnriched_row$deletion) & 
#                                     TF_relationship %in% TF_relationships_in_enrichment) |> 
#   select(gene_name) |> pull() |> unique()
# randomSubset_idxs <- sample(all_TFreg_idxs, length(anyEnriched_DE_idxs), replace = FALSE)
# # same thing as above but now we're using random subset of genes regulated in same direction by same TF, regardless of module
# # WT cerevisiae, enriched versus same-size subest of not-enriched from same module
# plotExpressionProfilePair(.cts1 = counts_all2$cer[,randomSubset_idxs],
#                           .cts2 = counts_all2$cer[,anyEnriched_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "not enriched subset of same size",
#                           .name2 = "enriched",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2") + ggtitle(paste(anyEnriched_row$deletion, 
#                                                                    "versus",  
#                                                                    anyEnriched_row$cer_color, 
#                                                                    anyEnriched_row$par_color, 
#                                                                    "module\nWT cerevisiae parent"))
# # repeat for paradoxus
# plotExpressionProfilePair(.cts1 = counts_all2$par[,randomSubset_idxs],
#                           .cts2 = counts_all2$par[,anyEnriched_DE_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "not enriched subset of same size",
#                           .name2 = "enriched",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2") + ggtitle(paste(anyEnriched_row$deletion, 
#                                                                    "versus",  
#                                                                    anyEnriched_row$cer_color, 
#                                                                    anyEnriched_row$par_color, 
#                                                                    "module\nWT paradoxus parent"))
# # skipping enriched comparison btwn species, as it's the same as above
# # how do the unenriched genes compare?
# # F1 hybrid
# plotExpressionProfilePair(.cts1 = counts_all2_allele$cer[, randomSubset_idxs],
#                           .cts2 = counts_all2_allele$par[, randomSubset_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "not enriched, S. cerevisiae F1 hybrid allele",
#                           .name2 = "not enriched, S. paradoxus F1 hybrid allele",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2") + ggtitle(paste(anyEnriched_row$deletion, 
#                                                                    anyEnriched_row$cer_color, 
#                                                                    anyEnriched_row$par_color, 
#                                                                    "module\nF1 hybrid"))
# # parent
# plotExpressionProfilePair(.cts1 = counts_all2$cer[, randomSubset_idxs],
#                           .cts2 = counts_all2$par[, randomSubset_idxs],
#                           .info1 = info,
#                           .info2 = info,
#                           .name1 = "not enriched, S. cerevisiae parent",
#                           .name2 = "not enriched, S. paradoxus parent",
#                           .show_points = TRUE,
#                           .method = "line",
#                           .normalization = "log2") + ggtitle(paste(anyEnriched_row$deletion, 
#                                                                    anyEnriched_row$cer_color, 
#                                                                    anyEnriched_row$par_color, 
#                                                                    "module\nparents"))
# # Based on these examples, it seems like the enriched subset tends to be a higher-expressed segment of the module. Is this true?
# # just check mean expr of enriched genes versus rest of module, as we don't have glm's within species and do not wish to make them unless absolutely necessary
# plotdf <- module_genedf10 |> select(gene_name, cer_color, par_color) |> 
#   left_join(tibble(gene_name = colnames(counts_all2$cer),
#                    mean_expr_cer = colMeans(counts_all2$cer),
#                    mean_expr_par = colMeans(counts_all2$par)),
#             by = "gene_name")
# # TODO: this left join results in gene-TF combinations getting thrown out if they're not enriched
# # which we don't want because we want to look at a specific TFs
# # data for genes that are both enriched and not enriched
# # so do it better
# enriched_cer_idxs <- filter(TFenrich_genedf, DE_cer & enriched_cer) |> 
#   select(gene_name) |> pull() |> unique()
# plotdf$inTFregGroup_cer <- plotdf$gene_name %in% enriched_cer_idxs
# enriched_par_idxs <- filter(TFenrich_genedf, DE_par & enriched_par) |> 
#   select(gene_name) |> pull() |> unique()
# plotdf$inTFregGroup_par <- plotdf$gene_name %in% enriched_par_idxs
# # in cer
# plotdf |> group_by(cer_color, par_color, inTFregGroup_cer) |> summarise(mean_expr_cer_group = mean(mean_expr_cer)) |> 
#   mutate(group_name = paste0(cer_color, par_color)) |> 
#   ggplot(aes(x = inTFregGroup_cer, y = log(mean_expr_cer_group))) + 
#   geom_jitter(aes(color = group_name)) +
#   theme(legend.position = "none")
# # in par
# plotdf |> group_by(cer_color, par_color, inTFregGroup_par) |> summarise(mean_expr_par_group = mean(mean_expr_par)) |> 
#   mutate(group_name = paste0(cer_color, par_color)) |> 
#   ggplot(aes(x = inTFregGroup_par, y = log(mean_expr_par_group))) + 
#   geom_jitter(aes(color = group_name)) +
#   theme(legend.position = "none")
# # definitely higher in the reg groups, possibly why they were the ones detected as DE
# 
# # TODO: Is my DE criteria too stringent? Is it only possible to achieve
# # as a highly expressed gene? 
# # It's also fine to just acknowledge that
# # my strict criteria leads to only a subset of the DE genes being detected,
# # as I'm not trying to say that this subset is special
# 
#### Mechanisms of co-regulatory divergence ####
# # We now have identified TF-regulated subsets that are strongly enriched
# # in CCMs, even after correcting for CCMs having a larger block size
# 
# # How do these subsets respond to TF deletion in parents and hybrids of both species?
# # They'll obviously be DE in the parent they were detected in, as that's how they were detected
# # But what about the other parent? Seeing as most of the subsets were only
# # detected in one parent, does the other parent's subset also respond
# # to the TF deletion just not as strongly, or is this a divergent response to the TF?
# 
# # Once we establish which one of those scenarios is true,
# # is this a cis change in these subset genes, or a trans change in the TF?
# # to address this, we will compare WT and deletion of subset genes between F1 hybrid alleles
# # If cis, there should a difference between alleles in WT, but not in deletion
# # If trans, there should not be a difference in WT or deletion
# 
# # You know the drill, first we do a randomly-chosen enriched subset,
# # then we quantify it in a way that makes sense,
# # then we quantify all subsets
# 
# # something that seems easier than it is:
# TFenrich_genedf <- TFdeldf_cer |> mutate(DE_cer = abs(coef) > coef_thresh & pval < p_thresh) |> 
#   select(gene_name, deletion, DE_cer)
# TFenrich_genedf <- TFdeldf_par |> mutate(DE_par = abs(coef) > coef_thresh & pval < p_thresh) |> 
#   select(gene_name, deletion, DE_par) |> 
#   full_join(y = TFenrich_genedf, by = c("gene_name", "deletion"))
# TFenrich_genedf <- TFenrich_genedf |> filter(DE_cer | DE_par)
# TFenrich_genedf25 <- left_join(TFenrich_genedf, 
#                                select(module_genedf25, gene_name, cer_color, par_color, CCM_color),
#                                by = "gene_name", relationship = "many-to-one")
# # repeat this section to inspect different enriched groups:
# coinFlip <- sample(c("heads", "tails"), 1)
# if (coinFlip == "heads") { # Heads we do an enriched group in cerevisiae
#   anyEnriched_row <- TFenrichdf10 |> filter(enriched_cer) |> slice_sample(n = 1)
#   anyEnriched_DE_idxs <- TFenrich_genedf25 |> filter(cer_color == as.character(anyEnriched_row$cer_color) &
#                                                        par_color == as.character(anyEnriched_row$par_color) &
#                                                        deletion == as.character(anyEnriched_row$deletion) &
#                                                        DE_cer) |> 
#     select(gene_name) |> pull()
# }
# if (coinFlip == "tails") { # Tails we do an enriched group in paradoxus
#   anyEnriched_row <- TFenrichdf10 |> filter(enriched_par) |> slice_sample(n = 1)
#   anyEnriched_DE_idxs <- TFenrich_genedf25 |> filter(cer_color == as.character(anyEnriched_row$cer_color) &
#                                                        par_color == as.character(anyEnriched_row$par_color) &
#                                                        deletion == as.character(anyEnriched_row$deletion) &
#                                                        DE_par) |> 
#     select(gene_name) |> pull()
# }
# anyEnriched_row |> t() # (note that if the group is enriched in both, we don't 
# # expect to see a difference in how parents respond to the deletion)
# deletion_name <- paste0(anyEnriched_row$deletion, "delete")
# p_cer <- plotExpressionProfilePair(.cts1 = counts_TFdel$cer[infos_TFdel$cer$genotype == "WT", 
#                                                             anyEnriched_DE_idxs],
#                                    .cts2 = counts_TFdel$cer[infos_TFdel$cer$genotype == deletion_name, 
#                                                             anyEnriched_DE_idxs],
#                                    .info1 = infos_TFdel$cer[infos_TFdel$cer$genotype == "WT",],
#                                    .info2 = infos_TFdel$cer[infos_TFdel$cer$genotype == deletion_name,],
#                                    .name1 = "WT",
#                                    .name2 = deletion_name,
#                                    .method = "line", .normalization = "log2", .show_points = TRUE) +
#   ggtitle(paste("enriched subset", anyEnriched_row$cer_color, anyEnriched_row$par_color, 
#                 "\nWT versus", deletion_name, "\nin cerevisiae parent"))
# p_par <- plotExpressionProfilePair(.cts1 = counts_TFdel$par[infos_TFdel$par$genotype == "WT", 
#                                                             anyEnriched_DE_idxs],
#                                    .cts2 = counts_TFdel$par[infos_TFdel$par$genotype == paste0(anyEnriched_row$deletion, "delete"), 
#                                                             anyEnriched_DE_idxs],
#                                    .info1 = infos_TFdel$par[infos_TFdel$par$genotype == "WT",],
#                                    .info2 = infos_TFdel$par[infos_TFdel$par$genotype == paste0(anyEnriched_row$deletion, "delete"),],
#                                    .name1 = "WT",
#                                    .name2 = deletion_name,
#                                    .method = "line", .normalization = "log2", .show_points = TRUE) +
#   ggtitle(paste("enriched subset", anyEnriched_row$cer_color, anyEnriched_row$par_color, 
#                 "\nWT versus", deletion_name, "\nin paradoxus parent"))
# # first plotting WT vs deletion for species where it was detected
# if (coinFlip == "heads") { # if it was detected in cer
#   p_cer
# }
# if (coinFlip == "tails") { # if it was detected in par
#   p_par
# }
# # now plotting WT vs deletion for species where it may or may not have been detected
# if (coinFlip == "heads") { # if it was detected in cer
#   p_par
# }
# if (coinFlip == "tails") { # if it was detected in par
#   p_cer
# }
# # Now is the divergence due to cis or trans-reg mechanisms?
# # plotting hyc vs hyp allele for WT
# p_hyWT <- plotExpressionProfilePair(.cts1 = counts_TFdel_allele$cer[infos_TFdel_allele$cer$genotype == "WT", 
#                                                                     anyEnriched_DE_idxs],
#                                     .cts2 = counts_TFdel_allele$par[infos_TFdel_allele$par$genotype == "WT", 
#                                                                     anyEnriched_DE_idxs],
#                                     .info1 = infos_TFdel_allele$cer[infos_TFdel_allele$cer$genotype == "WT",],
#                                     .info2 = infos_TFdel_allele$par[infos_TFdel_allele$par$genotype == "WT",],
#                                     .name1 = "F1 hybrid, cerevisiae allele",
#                                     .name2 = "F1 hybrid, paradoxus allele",
#                                     .method = "line", .normalization = "log2", .show_points = TRUE) +
#   ggtitle(paste("enriched subset", anyEnriched_row$cer_color, anyEnriched_row$par_color, 
#                 "\nWT in F1 hybrid"))
# p_hyTFdel <- plotExpressionProfilePair(.cts1 = counts_TFdel_allele$cer[infos_TFdel_allele$cer$genotype == deletion_name, 
#                                                                        anyEnriched_DE_idxs],
#                                        .cts2 = counts_TFdel_allele$par[infos_TFdel_allele$par$genotype == deletion_name, 
#                                                                        anyEnriched_DE_idxs],
#                                        .info1 = infos_TFdel_allele$cer[infos_TFdel_allele$cer$genotype == deletion_name,],
#                                        .info2 = infos_TFdel_allele$par[infos_TFdel_allele$par$genotype == deletion_name,],
#                                        .name1 = "F1 hybrid, cerevisiae allele",
#                                        .name2 = "F1 hybrid, paradoxus allele",
#                                        .method = "line", .normalization = "log2", .show_points = TRUE) +
#   ggtitle(paste("enriched subset", anyEnriched_row$cer_color, anyEnriched_row$par_color, 
#                 "\n", deletion_name, "in F1 hybrid"))
# # first plotting WT hyc vs hyp
# p_hyWT # if it's cis, here's where there should be a difference
# p_hyTFdel # neither should have a difference here, as the regulator is gone
