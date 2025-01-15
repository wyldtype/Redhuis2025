sapply(c("tidyverse", "ggpubr", "huxtable"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Barkai_data_analysis/")

# Script to construct Figure 1: individual gene models in two datasets (parent and hybrid) 
# associate divergence with expression variability and possibly expression level

# load single gene model dataset
load(file = "data_files/single_gene_models.RData")

# useful functions
getExprVector <- function(.gene_name, .organism = c("cer", "par", "hyb"), .experiment = c("TFdelxLowN", "CC", "HAP4", "LowPi")) {
  if (.organism %in% c("cer", "par")) {
    info <- spinfo[.experiment] %>% Reduce(f = bind_rows) # allows for multiple experiments
    cts <- spcts[.experiment] %>% Reduce(f = cbind)
    expr <- cts[.gene_name, info$organism == .organism]
  }
  if (.organism == "hyb") {
    info <- alinfo[.experiment] %>% Reduce(f = bind_rows)
    cts <- alcts[.experiment] %>% Reduce(f = cbind)
    expr_hyc <- cts[.gene_name, info$organism == .organism & info$allele == "cer"]
    expr_hyp <- cts[.gene_name, info$organism == .organism & info$allele == "par"]
    expr <- expr_hyc + expr_hyp
  }
  if (.organism == "hyc") {
    info <- alinfo[.experiment] %>% Reduce(f = bind_rows)
    cts <- alcts[.experiment] %>% Reduce(f = cbind)
    expr <- cts[.gene_name, info$organism == "hyb" & info$allele == "cer"]
  }
  if (.organism == "hyp") {
    info <- alinfo[.experiment] %>% Reduce(f = bind_rows)
    cts <- alcts[.experiment] %>% Reduce(f = cbind)
    expr <- cts[.gene_name, info$organism == "hyb" & info$allele == "par"]
  }
  return(expr)
}

getGeneDf <- function(.gene_name, .mode = c("parents", "hybrid"), .experiment = c("TFdelxLowN", "CC", "HAP4", "LowPi")) {
  if (.mode == "parents") {
    info <- spinfo[.experiment] %>% Reduce(f = bind_rows)
    cts <- spcts[.experiment] %>% Reduce(f = cbind)
  }
  if (.mode == "hybrid") {
    info <- alinfo[.experiment] %>% Reduce(f = bind_rows)
    cts <- alcts[.experiment] %>% Reduce(f = cbind)
  }
  genedf <- bind_cols(tibble(expr = cts[.gene_name,]), info) %>% 
    mutate(non_unique_sample_name = paste(condition, well_flask_ID)) %>% 
    select(expr, non_unique_sample_name, allele, experiment, time_point_num, genotype) %>% 
    pivot_wider(id_cols = c(non_unique_sample_name, experiment, time_point_num, genotype), values_from = expr, names_from = allele)
  return(genedf)
}
# tests for getGeneDf
gene_idx <- sample(rownames(spcts[[1]]), 1)
test <- getGeneDf(gene_idx, .mode = "parents")
ggplot(test, aes(x = cer, y = par)) + geom_point(aes(color = experiment)) + ggtitle(paste(gene_idx, "parents")) + geom_abline(color = "gold")
test <- getGeneDf(gene_idx, .mode = "hybrid") # hybrid should be less variable because cer and par measurements were collected in the same cell
ggplot(test, aes(x = cer, y = par)) + geom_point(aes(color = experiment)) + ggtitle(paste(gene_idx, "hybrid")) + geom_abline(color = "gold")
spaldf %>% filter(gene_name == gene_idx)

visualizeGeneExpressionScatterplot <- function(.gdf, .plotname = "individual gene counts", .color_from = "experiment", .log = FALSE) {
  if (!.log) {
    max_expr <- max(c(.gdf$cer, .gdf$par), na.rm = TRUE)
    output <- ggplot(.gdf, aes(x = cer, y = par)) + 
      geom_point(aes(color = pull(.gdf[, as.character(.color_from)]))) + geom_abline(color = "gold") +
      ggtitle(.plotname) + xlab("expression level in cerevisiae") + ylab("expression level in paradoxus") + 
      theme_classic() + xlim(c(0, max_expr)) + ylim(c(0, max_expr)) + theme(legend.title = element_blank())
    return(output)
  }
  if (.log) {
    max_expr <- max(c(log(.gdf$cer), log(.gdf$par)), na.rm = TRUE)
    output <- ggplot(.gdf, aes(x = log(cer), y = log(par))) + 
      geom_point(aes(color = pull(.gdf[, as.character(.color_from)]))) + geom_abline(color = "gold") +
      ggtitle(.plotname) + xlab("expression level in cerevisiae") + ylab("expression level in paradoxus") + 
      theme_classic() + xlim(c(0, max_expr)) + ylim(c(0, max_expr)) + theme(legend.title = element_blank())
    return(output)
  }
}
# TDH3
getGeneDf("YGR192C", .mode = "parents") %>% visualizeGeneExpressionScatterplot(.plotname = "YGR192C - parents")
getGeneDf("YGR192C", .mode = "hybrid") %>% visualizeGeneExpressionScatterplot(.plotname = "YGR192C - hybrid")
getGeneDf("YGR192C", .mode = "parents") %>% visualizeGeneExpressionScatterplot(.plotname = "YGR192C - parents - log scale", .log = TRUE)
# TDH3 WT/YPD only
test <- getGeneDf("YGR192C", .mode = "parents")
test$is_ypd_wt <- test$genotype == "WT" & test$time_point_num <= 0
test %>% visualizeGeneExpressionScatterplot(.plotname = "YGR192C - parents", .color_from = "is_ypd_wt", .log = TRUE)
filter(test, is_ypd_wt) %>% visualizeGeneExpressionScatterplot(.plotname = "YGR192C - parents", .log = TRUE)

# random gene
gene_idx <- sample(rownames(spcts[[1]]), 1)
getGeneDf(gene_idx, .mode = "parents") %>% visualizeGeneExpressionScatterplot(.plotname = paste(gene_idx, "parents", sep = " - "))
getGeneDf(gene_idx, .mode = "hybrid") %>% visualizeGeneExpressionScatterplot(.plotname = paste(gene_idx, "hybrid", sep = " - "))
spaldf %>% filter(gene_name == gene_idx)

# random sig gene
gene_idx <- spaldf %>% filter(sig) %>% select(gene_name) %>% pull() %>% sample(1)
getGeneDf(gene_idx, .mode = "parents") %>% visualizeGeneExpressionScatterplot(.plotname = paste(gene_idx, "parents", sep = " - "), .log = FALSE)
getGeneDf(gene_idx, .mode = "hybrid") %>% visualizeGeneExpressionScatterplot(.plotname = paste(gene_idx, "hybrid", sep = " - "), .log = FALSE)
spaldf %>% filter(gene_name == gene_idx)

# get the set of names of genes that matches your criteria
# @input: .criteria in the form of a 
# @output: vector of gene names
getGeneIdx <- function(.spaldf, .criteria, .mode = c("parents", "hybrid")) {
  # TODO: I think this will be better long term than the unhinged pipeline things I have in Fig 1b
  # But I don't quite know how to define the criteria
}

# re-call sig if you want to change thresholds
p_thresh <- 0.05/(length(unique(spaldf$gene_name))*5)
eff_thresh <- 0.75
spaldf$sig <- spaldf$pvalue < p_thresh & abs(spaldf$effect_size) > eff_thresh & spaldf$isexpressed
sum(spaldf$sig)
spaldf %>% filter(coefficient == "species") %>% group_by(gene_name) %>% summarise(sig_gene = any(sig)) %>% select(sig_gene) %>% pull() %>% sum()
# re-creating genedf, where each row is one gene
genedf <-  spaldf %>% 
  pivot_wider(names_from = c("experiment", "coefficient"), values_from = c("effect_size", "pvalue", "sig"), id_cols = c("gene_name"))
table(genedf$gene_name) %>% table() # checking that every gene is only represented once

###################### Figure 1 ###########################
### Figure 1a: How effect size (log fold change) relates to expression divergence ####
library(MASS, include.only = "glm.nb")
# since it was sooooo fun to set the label positions below, I did YEL060C LowPi for reference
random_sig_gene <- filter(spaldf, gene_name == "YEL060C" & experiment == "LowPi") %>% select(gene_name, experiment, effect_size) %>% slice_sample(n = 1)
# random_sig_gene <- filter(spaldf, sig & experiment != "all") %>% select(gene_name, experiment, effect_size) %>% slice_sample(n = 1)
random_sig_gene # based on its experiment, you'll want to change model
random_sig_gdf <- getGeneDf(random_sig_gene$gene_name, .mode = "parents", .experiment = random_sig_gene$experiment) %>% pivot_longer(cols = c(cer, par), values_to = "expr", names_to = "allele")
random_sig_gdf$allele <- as.factor(random_sig_gdf$allele) %>% relevel(ref = "par")
spmod <- glm.nb(expr ~ time_point_num + allele, data = random_sig_gdf)
sprow <- summary(spmod)$coefficients["allelecer",, drop = FALSE]
slope <- 1/exp(sprow[1])
intercept <- summary(spmod)$coefficients["(Intercept)",, drop = FALSE][1]
plotdf <- random_sig_gdf %>% mutate(fittedvals = spmod$fitted.values) %>% select(allele, expr, fittedvals, genotype, non_unique_sample_name) %>% 
  pivot_longer(cols = c(expr, fittedvals)) %>% 
  pivot_wider(id_cols = c(non_unique_sample_name, name, genotype), names_from = allele, values_from = value)
ggplot(plotdf, aes(x = cer, y = par)) + geom_point(aes(color = name)) + geom_abline(color = "gold", slope = slope, intercept = intercept) + 
  geom_abline(color = "navy", slope = 1, intercept = 0) + 
  xlim(c(0, max(select(plotdf, cer, par)))) + ylim(c(0, max(select(plotdf, cer, par)))) + 
  theme_classic() + ggtitle(paste(random_sig_gene$gene_name, random_sig_gene$experiment, sep = " - ")) + 
  annotate("text", label = paste("Log Fold Change (LFC) = ", round(sprow[1], digits = 2), 
                                 "\nslope = 1/e^(LFC) =", round(slope, digits = 2),
                                 "\npvalue = ", round(sprow[4], digits = 10)), 
           x = 450, y = 100, size = 2) +
  annotate("text", color = "gold", x = 425, y = intercept + 335*slope, label = "Model", size = 3) +
  annotate("text", color = "navy", x = 600, y = 540, label = "y = x", size = 3)

#### Figure 1b: upset plots of nsig per dataset, species and allele ####
library(ComplexHeatmap)
# species
plotdf <- genedf[, paste("sig", names(spcts), "species", sep = "_")] %>% 
  apply(MARGIN = 2, FUN = as.integer) %>% as.matrix()
colnames(plotdf) <- names(spcts)
# replacing NA with 0 because make_comb_mat cannot handle NA for some reason
plotdf[is.na(plotdf)] <- 0
plotdf <- make_comb_mat(plotdf)
p <- UpSet(plotdf, set_order = names(spcts), comb_order = order(comb_size(plotdf)), # this, I will say, makes ggplot look like a freakin' dream. I literally just wanted to rename the histogram y axis and add counts to the top of the bars
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
# accompanying number of samples barplot:
spinfo_full %>% ggplot(aes(y = experiment)) + 
  geom_bar(aes(fill = organism), position = position_dodge()) + 
  scale_fill_manual(values = c("gold", "mediumorchid")) + 
  theme_classic() + xlab("number of samples") + scale_y_discrete(limits = c("LowPi", "HAP4", "CC", "LowN"))
# allele
plotdf <- genedf[, paste("sig", names(spcts), "allele", sep = "_")] %>% 
  apply(MARGIN = 2, FUN = as.integer) %>% as.matrix()
colnames(plotdf) <- names(spcts)
# replacing NA with 0
plotdf[is.na(plotdf)] <- 0
plotdf <- make_comb_mat(plotdf)
p <- UpSet(plotdf, set_order = names(spcts), comb_order = order(comb_size(plotdf)),
           top_annotation = HeatmapAnnotation(
             "number of genes" = anno_barplot(comb_size(plotdf), 
                                              ylim = c(0, max(comb_size(plotdf))*1.1),
                                              border = FALSE, 
                                              gp = gpar(fill = "black"), 
                                              height = unit(4, "cm")
             ), 
             annotation_name_side = "left", 
             annotation_name_rot = 90))
draw(p)
decorate_annotation("number of genes", {
  grid.text(comb_size(plotdf)[column_order(p)], x = seq_along(comb_size(plotdf)), y = unit(comb_size(plotdf)[column_order(p)], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
}) 
# accompanying number of samples barplot:
alinfo_full %>% ggplot(aes(y = experiment)) + 
  geom_bar(aes(fill = allele), position = position_dodge()) + 
  scale_fill_manual(values = c("gold", "mediumorchid")) + 
  theme_classic() + xlab("number of samples") + scale_y_discrete(limits = c("LowPi", "HAP4", "CC", "LowN"))

#### Fig 1c: Examples of genes in different categories of upset plots ####
# divergent (env-specific and constitutively) and non-divergent gene examples and how effect size relates to divergence 
# env-specific divergent YPL019C
esdiv_idx <- spaldf %>% filter(coefficient == "species" & experiment != "all") %>% group_by(gene_name) %>%
  summarise(is_esdiv = sum(sig) < 4 & sum(sig) > 0) %>% 
  filter(is_esdiv) %>% select(gene_name) %>% pull() %>% sample(1)
esdiv_plot_parents <- getGeneDf(esdiv_idx, .mode = "parents") %>% visualizeGeneExpressionScatterplot(.plotname = paste(esdiv_idx, "Diverged,\nEnvironment-Specific"))
esdiv_plot_hybrid <- getGeneDf(esdiv_idx, .mode = "hybrid") %>% visualizeGeneExpressionScatterplot(.plotname = paste(esdiv_idx, "Hybrid alleles"))
esdiv_plot_parents
esdiv_plot_hybrid
spaldf %>% filter(coefficient == "species" & gene_name == esdiv_idx)
spaldf %>% filter(coefficient == "allele" & gene_name == esdiv_idx)
# constitutively diverged
consdiv_idx <- spaldf %>% filter(coefficient == "species" & experiment != "all") %>% group_by(gene_name) %>%
  summarise(is_consdiv = sum(sig) == 4) %>% 
  filter(is_consdiv) %>% select(gene_name) %>% pull() %>% sample(1)
consdiv_plot_parents <- getGeneDf(consdiv_idx, .mode = "parents") %>% visualizeGeneExpressionScatterplot(.plotname = paste(consdiv_idx, "Diverged,\nConstitutively"))
consdiv_plot_hybrid <- getGeneDf(consdiv_idx, .mode = "hybrid") %>% visualizeGeneExpressionScatterplot(.plotname = paste(consdiv_idx, "Hybrid alleles"))
consdiv_plot_parents
consdiv_plot_hybrid
spaldf %>% filter(coefficient == "species" & gene_name == consdiv_idx)
spaldf %>% filter(coefficient == "allele" & gene_name == consdiv_idx)
# not diverged
notdiv_idx <- spaldf %>% filter(coefficient == "species" & experiment != "all") %>% group_by(gene_name) %>%
  summarise(is_notdiv = sum(sig) == 0) %>% 
  filter(is_notdiv) %>% select(gene_name) %>% pull() %>% sample(1)
notdiv_plot_parents <- getGeneDf(notdiv_idx, .mode = "parents") %>% visualizeGeneExpressionScatterplot(.plotname = paste(notdiv_idx, "Not Diverged"))
notdiv_plot_hybrid <- getGeneDf(notdiv_idx, .mode = "hybrid") %>% visualizeGeneExpressionScatterplot(.plotname = paste(notdiv_idx, "Hybrid alleles"))
notdiv_plot_parents
notdiv_plot_hybrid
spaldf %>% filter(coefficient == "species" & gene_name == notdiv_idx)
spaldf %>% filter(coefficient == "allele" & gene_name == notdiv_idx)
# plotting all 3 genes in parents and hybrids
p <- ggarrange(plotlist = list(esdiv_plot_parents, 
                          esdiv_plot_hybrid, 
                          consdiv_plot_parents, 
                          consdiv_plot_hybrid, 
                          notdiv_plot_parents, 
                          notdiv_plot_hybrid), common.legend = TRUE, nrow = 3, ncol = 2)
annotate_figure(p)

#### Fig 1d: some sort of effect size (lfc par to cer) summary for all expressed genes with high enough pvalue in parent and hybrid ####
plotdf <- spaldf %>% filter(pvalue > p_thresh & isexpressed & experiment != "all") %>% pivot_wider(id_cols = c(gene_name, experiment), names_from = coefficient, values_from = effect_size)
ggplot(plotdf, aes(x = sign(species)*log(abs(species) + 1e-9))) + geom_histogram()
plotdf$allele[is.na(plotdf$allele)] <- 0
plotdf$species[is.na(plotdf$species)] <- 0

ggplot(plotdf, aes(x = sign(allele)*log(abs(allele)), y = sign(species)*log(abs(species)))) + 
  geom_point(aes(color = experiment)) + theme_classic()

#### Fig 1e: expression level and variance/max don't explain divergence as measured by effect size ####

################################# Archive ######################################
# parsing glm.nb models into dataframes using avg_comparisons
# Archived b/c it's not as intuitive as log fold change
# converts glm.nb model to one row of an avg comparison df 
# modToMargEffectSpeciesAllele <- function(.mod, .genedf) {
#   avgcomps <- tryCatch({
#     comps <- avg_comparisons(.mod, variables = "allele", newdata = .genedf)
#     output <- tibble(gene_name = .genedf$gene_name[1],
#                      coefficient = .genedf$allele_or_species[1],
#                      estimate = comps$estimate,
#                      pvalue = comps$p.value)
#     return(output)
#   }, error = function(e) {
#     return(NA)
#   }, warning = function(w) {
#     return(NA)
#   })
#   if (is.na(avgcomps)) {
#     output <- tibble(gene_name = .genedf$gene_name[1],
#                      coefficient = .genedf$allele_or_species[1],
#                      estimate = NA,
#                      pvalue = NA)
#     return(output)
#   }
#   return(output)
# }
# parseAvgCompsFromModslist <- function(.modslist, .cts, .info) {
#   good_idxs <- which(!is.na(.modslist))
#   result <- map2(.modslist[good_idxs], names(.modslist)[good_idxs], function(mod, gname) {
#     cat(gname, ":", which(names(.modslist[good_idxs]) == gname), "/", length(good_idxs), "\n")
#     gdf <- bind_cols(tibble(expr = .cts[gname,],
#                             gene_name = gname,
#                             allele_or_species = if_else(all(.info$organism == "hyb"), true = "allele", false = "species")),
#                      .info)
#     return(modToMargEffectSpeciesAllele(mod, gdf))
#   }) %>% Reduce(f = bind_rows)
#   return(result)
# }
# # tests for parseAvgCompsFromModslist
# dataset_idx <- sample(c(1:length(spmods)), 1)
# test_modslist <- spmods[[dataset_idx]]
# # don't run this whole thing:
# test <- parseAvgCompsFromModslist(test_modslist, spcts[[dataset_idx]], spinfo[[dataset_idx]])
# spavgs <- map(c(1:length(spmods)), function(i) {return(parseAvgCompsFromModslist(spmods[[i]], spcts[[i]], spinfo[[i]]))})
# names(spavgs) <- names(spdfs)
# spavgs <- map2(spavgs, names(spavgs), function(x, y) {
#   x$experiment <- y
#   return(x)
# })
# spavg <- Reduce(f = bind_rows, spavgs)

# ###################### computing individual gene metrics ###########################
# Archived b/c it's not particularly interesting to say that expression level/var/etc. doesn't predict gene expression divergence
# could be revived for a supplement that backs up a single sentence in the manuscript or something
# # adding gene metrics for each organism/experiment
for (o in c("cer", "par", "hyb")) {
  for (e in ExperimentNames) {
    metrics <- lapply(GeneNames, \(g) {
      cat(o, e, g, which(GeneNames == g), "/", nGenes, "\n")
      expr <- getExprVector(g, .organism = o, .experiment = e)
      scaled_var <- var(expr/max(expr))
      mean_expr <- mean(expr)
      max_expr <- max(expr)
      return(tibble("scaled_var" = scaled_var, "mean_expr" = mean_expr, "max_expr" = max_expr))
    }) %>% Reduce(f = bind_rows)
    genedf[, paste("scaled_var", e, o, sep = "_")] <- metrics$scaled_var
    genedf[, paste("mean", e, o, sep = "_")] <- metrics$mean_expr
    genedf[, paste("max", e, o, sep = "_")] <- metrics$max_expr
  }
}

# QC: is there a mean-scaled_var relationship?
library(ggpubr)
plots <- lapply(ExperimentNames, \(e) {
  output <- vector(mode = "list", length = 3)
  names(output) <- c("cer", "par", "hyb")
  for (o in c("cer", "par", "hyb")) {
    df <- select(genedf, paste("scaled_var", e, o, sep = "_"), paste("mean", e, o, sep = "_"))
    colnames(df) <- c("scaled_var", "mean")
    p <- ggplot(df, aes(x = log(scaled_var), y = log(mean))) + geom_point() + xlab("Mean") + ylab("Var/Max") + ggtitle(paste(o, e))
    output[[o]] <- p
  }
  return(output)
})
# TFdelxLowN
ggarrange(plotlist = plots[[1]], ncol = 3, nrow = 1, common.legend = TRUE)
# CC
ggarrange(plotlist = plots[[2]], ncol = 3, nrow = 1, common.legend = TRUE)
# HAP4
ggarrange(plotlist = plots[[3]], ncol = 3, nrow = 1, common.legend = TRUE)
# LowPi
ggarrange(plotlist = plots[[4]], ncol = 3, nrow = 1, common.legend = TRUE)
# really no relationship! Also hybrid def has lower variability on average, which is interesting
# TODO: this looks a lot better than when I didn't separate experiments. Repeat with Tau to see if it also looks better and check if scaled var and tau are correlated now

# QC: how often do scaled_vars and mean_expr of each gene correlate between species?
plots <- lapply(ExperimentNames, \(e) {
  df <- select(genedf, paste("scaled_var", e, "cer", sep = "_"), paste("scaled_var", e, "par", sep = "_"))
  colnames(df) <- c("cer", "par")
  p <- ggplot(df, aes(x = cer, y = par)) + geom_point() + xlab("Cer") + ylab("Par") + ggtitle(e)
  return(p)
})
p <- ggarrange(plotlist = plots, ncol = 2, nrow = 2, common.legend = TRUE)
annotate_figure(p, top = text_grob("Expression Variability var(expr/max_expr)"))
# like we've seen before, correlated but far from identical
# How often do means correlate?
plots <- lapply(ExperimentNames, \(e) {
  df <- select(genedf, paste("mean", e, "cer", sep = "_"), paste("mean", e, "par", sep = "_"))
  colnames(df) <- c("cer", "par")
  p <- ggplot(df, aes(x = log(cer), y = log(par))) + geom_point() +
    xlab("Cer") + ylab("Par") + ggtitle(e) + theme_classic()
  return(p)
})
p <- ggarrange(plotlist = plots, ncol = 2, nrow = 2, common.legend = TRUE)
annotate_figure(p, top = text_grob("Mean Expression (log scale)"))
# way stronger correlation

# calculating median expression ratio (cer/par)
med_then_ratio <- map_dbl(GeneNames, \(g) {
  cat(g, which(GeneNames == g), "/", nGenes, "\n")
  output <- median(getExprVector(g, .organism = "cer"))/median(getExprVector(g, .organism = "par"))
  return(output)
})
ratio_then_med <- map_dbl(GeneNames, \(g) {
  cat(g, which(GeneNames == g), "/", nGenes, "\n")
  gdf <- getGeneDf(g, .mode = "parents") %>% mutate(expr_ratio = cer/par)
  ratios <- gdf %>% filter(is.finite(expr_ratio)) %>% select(expr_ratio) %>% pull()
  return(median(ratios))
})
# QC: are the two ways of measuring median ever different?
good_idxs <- which(is.finite(ratio_then_med) & is.finite(med_then_ratio))
plotdf <- tibble(med_then_ratio = med_then_ratio[good_idxs],
                 ratio_then_med = ratio_then_med[good_idxs])
ggplot(plotdf, aes(x = med_then_ratio, y = ratio_then_med)) + geom_point()
sum(is.finite(ratio_then_med))
sum(is.finite(med_then_ratio)) # basically equivalent
# TODO: also check with .mode = "hybrid"

# we'll use ratio_then_med because it feels like it takes advantage of the dataset more (even though we've just shown that it doesn't matter)
genedf$median_expr_ratio <- ratio_then_med

# calculating expression variability and mean metrics
calculateTau <- function(.expr) {
  m <- length(.expr)
  maxexpr <- max(.expr)
  unscaled_tau <- (1 - .expr/maxexpr) %>% sum(na.rm = TRUE)
  return(unscaled_tau/(m - 1))
}
calculateTau(getExprVector(.gene_name = "YGR192C", .organism = "cer"))
calculateTau(getExprVector(.gene_name = "YGR192C", .organism = "par"))
calculateTau(getExprVector(.gene_name = "YGR192C", .organism = "hyb"))

# trying to recreate Gat's no correlation btwn expr level divergence and regulatory dynamic divergence plot
ggplot(genedf, aes(x = median_expr_ratio, y = abs(scaled_var_cer - scaled_var_par))) + geom_point() # alright close enough

# # Making plots based on marginal effects estimates instead of median expression ratios
# # pros: more intuitive axes actual expression counts
# # cons: doesn't control for expression level, so changes in highly expressed genes like TDH3 are heavily emphasized
# # all genes
# ggplot(plotdf, aes(x = species_estimate, y = allele_estimate, label = plot_name)) +
#   geom_text(check_overlap = TRUE, hjust = 0, nudge_x = 300, size = 3, aes(color = pct_cis_bin)) +
#   geom_point(aes(color = pct_cis_bin)) +
#   theme_classic() + ggtitle("average change in expression level\n shifting from cerevisiae to paradoxus") +
#   xlab("between species") + ylab("between hybrid alleles") +
#   theme(legend.position = "none") + geom_rect(xmin = -5000, xmax = 5000, ymin = -5000, ymax = 5000, color = "gold", alpha = 0) +
#   xlim(c(-17000, 17000)) + ylim(c(-17000, 17000))
# # effect +-5000
# ggplot(plotdf, aes(x = species_estimate, y = allele_estimate, label = plot_name)) +
#   geom_text(check_overlap = TRUE, hjust = 0, nudge_x = 300, size = 3, aes(color = pct_cis_bin)) +
#   geom_point(aes(color = pct_cis_bin)) +
#   theme_classic() + ggtitle("average change in expression level\n shifting from cerevisiae to paradoxus") +
#   xlab("between species") + ylab("between hybrid alleles") +
#   theme(legend.position = "none") + geom_rect(xmin = -5000, xmax = 5000, ymin = -5000, ymax = 5000, color = "gold", alpha = 0) +
#   geom_rect(xmin = -1000, xmax = 1000, ymin = -1000, ymax = 1000, color = "mediumorchid", alpha = 0) +
#   xlim(c(-5000, 5000)) + ylim(c(-5000, 5000))
# # effect +-1000
# ggplot(plotdf, aes(x = species_estimate, y = allele_estimate, label = plot_name)) +
#   geom_point(aes(color = pct_cis_bin)) +
#   theme_classic() + ggtitle("average change in expression level\n shifting from cerevisiae to paradoxus") +
#   xlab("between species") + ylab("between hybrid alleles") +
#   theme(legend.text = element_blank(), legend.title = element_blank()) +
#   geom_rect(xmin = -1000, xmax = 1000, ymin = -1000, ymax = 1000, color = "mediumorchid", alpha = 0) +
#   xlim(c(-1000, 1000)) + ylim(c(-1000, 1000))
# # log scale, just looking at magnitude not direction of effect estimates now
# ggplot(plotdf, aes(x = log(abs(species_estimate)), y = log(abs(allele_estimate)), label = gene_name)) +
#   geom_text(check_overlap = TRUE, hjust = 0, nudge_x = 0.1, size = 3, aes(color = pct_cis_bin)) +
#   geom_point(aes(color = pct_cis_bin)) +
#   theme_classic() + ggtitle("average change in expression level shifting\n from cerevisiae to paradoxus (log scale)") +
#   xlab("between species") + ylab("between hybrid alleles") +
#   theme(legend.position = "none")
# 
# # # Predicting Divergence with Expression Level and Variance
# # Archived because it's too close to Gat's Paper and too far from what's interesting about the data---hybrid expression differences illustrate how allelic effects can propagate through the network
# var_cutoff <- quantile(unlist(genedf[,grep("scaled_var", colnames(genedf))]), 0.9, na.rm = TRUE) %>% as.numeric()
# genedf$scaled_var_cer <- colMeans(rbind(genedf$scaled_var_TFdelxLowN_cer, genedf$scaled_var_CC_cer, genedf$scaled_var_HAP4_cer, genedf$scaled_var_LowPi_cer))
# genedf$scaled_var_par <- colMeans(rbind(genedf$scaled_var_TFdelxLowN_par, genedf$scaled_var_CC_par, genedf$scaled_var_HAP4_par, genedf$scaled_var_LowPi_par))
# # between species
# high_var_div <- which(unlist(apply(genedf[, paste("sig", ExperimentNames, "species", sep = "_")], 1, all)) & (genedf$scaled_var_cer > var_cutoff | genedf$scaled_var_par > var_cutoff)) %>% sample(size = 1)
# high_var_notDiv <- which(unlist(apply(!genedf[, paste("sig", ExperimentNames, "species", sep = "_")], 1, any)) & (genedf$scaled_var_cer > var_cutoff | genedf$scaled_var_par > var_cutoff)) %>% sample(size = 1)
# low_var_div <- which(unlist(apply(genedf[, paste("sig", ExperimentNames, "species", sep = "_")], 1, all)) & (genedf$scaled_var_cer <= var_cutoff & genedf$scaled_var_par <= var_cutoff)) %>% sample(size = 1)
# low_var_notDiv <- which(unlist(apply(!genedf[, paste("sig", ExperimentNames, "species", sep = "_")], 1, any)) & (genedf$scaled_var_cer <= var_cutoff & genedf$scaled_var_par <= var_cutoff)) %>% sample(size = 1)
# # parent plots
# max_expr <- lapply(GeneNames[c(high_var_div, high_var_notDiv, low_var_div, low_var_notDiv)], getGeneDf, .mode = "parents") %>% Reduce(f = bind_rows) %>% select(cer, par) %>% pull() %>% max()
# hvd_plot <- visualizeGeneExpressionScatterplot(genedf$gene_name[high_var_div], .mode = "parents") + xlim(c(0, max_expr)) + ylim(c(0, max_expr))
# hvnd_plot <- visualizeGeneExpressionScatterplot(genedf$gene_name[high_var_notDiv], .mode = "parents") + xlim(c(0, max_expr)) + ylim(c(0, max_expr))
# lvd_plot <- visualizeGeneExpressionScatterplot(genedf$gene_name[low_var_div], .mode = "parents") + xlim(c(0, max_expr)) + ylim(c(0, max_expr))
# lvnd_plot <- visualizeGeneExpressionScatterplot(genedf$gene_name[low_var_notDiv], .mode = "parents") + xlim(c(0, max_expr)) + ylim(c(0, max_expr))
# ggarrange(plotlist = list(hvd_plot, hvnd_plot, lvd_plot, lvnd_plot))
# 
# # hybrid plots
# max_expr <- lapply(GeneNames[c(high_var_div, high_var_notDiv, low_var_div, low_var_notDiv)], getGeneDf, .mode = "hybrid") %>% Reduce(f = bind_rows) %>% select(cer, par) %>% pull() %>% max()
# hvd_plot <- visualizeGeneExpressionScatterplot(genedf$gene_name[high_var_div], .mode = "hybrid") + xlim(c(0, max_expr)) + ylim(c(0, max_expr))
# hvnd_plot <- visualizeGeneExpressionScatterplot(genedf$gene_name[high_var_notDiv], .mode = "hybrid") + xlim(c(0, max_expr)) + ylim(c(0, max_expr))
# lvd_plot <- visualizeGeneExpressionScatterplot(genedf$gene_name[low_var_div], .mode = "hybrid") + xlim(c(0, max_expr)) + ylim(c(0, max_expr))
# lvnd_plot <- visualizeGeneExpressionScatterplot(genedf$gene_name[low_var_notDiv], .mode = "hybrid") + xlim(c(0, max_expr)) + ylim(c(0, max_expr))
# ggarrange(plotlist = list(hvd_plot, hvnd_plot, lvd_plot, lvnd_plot))


