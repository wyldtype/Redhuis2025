sapply(c("dplyr", "purrr", "ggplot2", "ggpubr", "tidyr", "data.table", "energy", "WGCNA", "bipartite"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Redhuis2024/")

# load functions and datasets
source(file = "functions_for_figure_scripts.R")

# Figure 1: Most genes have co-expression partners that are conserved between cerevisiae and paradoxus

#### What proportion of genes are in modules/CCMs? ####

# How many modules per species?
table(colors00$cer) # 14, not including grey
table(colors10$cer) # 11
table(colors25$cer) # 6
table(colors35$cer) # 6
table(colors00$par) # 12
table(colors10$par) # 11
table(colors25$par) # 9
table(colors35$par) # 7

# How many CCMs per colorset?
# colors 00:
CCMdf00 %>% filter(p_value < 0.05 &
                     cer_color != "grey" & par_color != "grey") %>% nrow()
# colors 10:
CCMdf10 %>% filter(p_value < 0.05 &
                     cer_color != "grey" & par_color != "grey") %>% nrow()
# colors 25:
CCMdf25 %>% filter(p_value < 0.05 &
                     cer_color != "grey" & par_color != "grey") %>% nrow()
# colors 35:
CCMdf35 %>% filter(p_value < 0.05 &
                     cer_color != "grey" & par_color != "grey") %>% nrow()

# What proportion of genes are in CCMs?
# nConserved/nDivergent for colors00:
sum(module_genedf00$is_CCM)
sum(!module_genedf00$is_CCM)
sum(module_genedf00$is_CCM)/nrow(module_genedf00)
# 2417/4794 = 50% conserved

# nConserved/nDivergent for colors10:
sum(module_genedf10$is_CCM)
sum(!module_genedf10$is_CCM)
sum(module_genedf10$is_CCM)/nrow(module_genedf10)
# 2555/4794 = 53% conserved

# nConserved/nDivergent for colors25:
sum(module_genedf25$is_CCM)
sum(!module_genedf25$is_CCM)
sum(module_genedf25$is_CCM)/nrow(module_genedf25)
# 2762/4794 = 58% conserved

# nConserved/nDivergent for colors35:
sum(module_genedf35$is_CCM)
sum(!module_genedf35$is_CCM)
sum(module_genedf35$is_CCM)/nrow(module_genedf35)
# 2793/4794 = 58% conserved

#### Any plots to use in workflow image, as necessary ####

# expr vs timepoint for random gene pair i and j in each experiment and each species, 
# then showing how correlation between two genes' expression is used to build network and define modules

# Getting a highly correlated pair of genes (highly correlated in cerevisiae)
cor_mat <- cor(counts_all2$cer, use = "pairwise.complete.obs")
diag(cor_mat) <- NA
which(cor_mat == sort(cor_mat, decreasing = TRUE)[26], arr.ind = TRUE) # change index from 1-10ish to see top 10 or so pairs
# correlated_gene_pair <- c("YKR097W", "YJR095W") # these two are super correlated because of a few instances where they're both almost identically highly expressed ("YLR377C"/"YLR174W", "YNL117W"/"YER065C" are like that too)
# correlated_gene_pair <- c("YGR034W", "YBL092W") # these two are both 60S subunit members ("YDR064W"/"YDL083C", "YOR369C"/"YNL178W" are also ribosomal)
# correlated_gene_pair <- c("YJR010W", "YIL074C") # MET3/SER33
# correlated_gene_pair <- c("YER024W", "YCR005C")
correlated_gene_pair <- c("YML063W", "YJR123W") # the most immediately obviously correlated example

# or you could do an actually random pair (but corr is usually low and just looks like a cloud of points, not as good for a figure)
random_gene_pair <- sample(colnames(counts_all2$cer), 2, replace = FALSE)

# getting an obviously low correlated gene pair for figure (but highly expressed to reduce noise for interpretability in the small plot)
which(abs(cor_mat["YBL003C",]) < 0.1 & colMeans(counts_all2$cer) > 1000) # YBL003C random highly expressed gene
lowcor_gene_pair <- c("YBL003C", "YML123C")
colMeans(counts_all2$cer)[which(colnames(counts_all2$cer) %in% lowcor_gene_pair)]

# cer correlated
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/DefiningModules/genePair_correlated_cer.pdf",
    width = 12, height = 2)
plotExpressionProfilePair(.cts1 = counts_all2$cer[,correlated_gene_pair[1], drop = FALSE],
                          .cts2 = counts_all2$cer[,correlated_gene_pair[2], drop = FALSE],
                          .info1 = info,
                          .info2 = info,
                          .name1 = correlated_gene_pair[1],
                          .name2 = correlated_gene_pair[2],
                          .color1 = "black",
                          .color2 = "grey50",
                          .method = "line", .show_points = FALSE, .normalization = "log2")
dev.off()

# par correlated
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/DefiningModules/genePair_correlated_par.pdf",
    width = 12, height = 2)
plotExpressionProfilePair(.cts1 = counts_all2$par[,correlated_gene_pair[1], drop = FALSE],
                          .cts2 = counts_all2$par[,correlated_gene_pair[2], drop = FALSE],
                          .info1 = info,
                          .info2 = info,
                          .name1 = correlated_gene_pair[1],
                          .name2 = correlated_gene_pair[2],
                          .color1 = "black",
                          .color2 = "grey50",
                          .method = "line", .show_points = FALSE, .normalization = "log2")
dev.off()

# cer low cor
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/DefiningModules/genePair_random_cer.pdf",
    width = 12, height = 2)
plotExpressionProfilePair(.cts1 = counts_all2$cer[,lowcor_gene_pair[1], drop = FALSE],
                          .cts2 = counts_all2$cer[,lowcor_gene_pair[2], drop = FALSE],
                          .info1 = info,
                          .info2 = info,
                          .name1 = lowcor_gene_pair[1],
                          .name2 = lowcor_gene_pair[2],
                          .color1 = "black",
                          .color2 = "grey50",
                          .method = "line", .show_points = FALSE, .normalization = "log2")
dev.off()

# par low cor
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/DefiningModules/genePair_random_par.pdf",
    width = 12, height = 2)
plotExpressionProfilePair(.cts1 = counts_all2$par[,lowcor_gene_pair[1], drop = FALSE],
                          .cts2 = counts_all2$par[,lowcor_gene_pair[2], drop = FALSE],
                          .info1 = info,
                          .info2 = info,
                          .name1 = lowcor_gene_pair[1],
                          .name2 = lowcor_gene_pair[2],
                          .color1 = "black",
                          .color2 = "grey50",
                          .method = "line", .show_points = FALSE, .normalization = "log2")
dev.off()

# cer and par gene trees, showing how merge similarity affects module definitions
# figure DefiningModules workflow: comparison of module definitions for cer, par, hyc, and hyp in merge 00, merge 10, merge 25
# cer
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/DefiningModules/geneTrees_cer.pdf",
    width = 5, height = 4)
plotDendroAndColors(dendro = geneTrees$cer, colors = cbind(colors00$cer, colors10$cer, colors25$cer, colors35$cer),
                    groupLabels = c("unmerged", "merged at 0.90", "merged at 0.75", "merged at 0.65"), main = "Saccharomyces cerevisiae",
                    dendroLabels = FALSE)
dev.off()
# par
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/DefiningModules/geneTrees_par.pdf",
    width = 5, height = 4)
plotDendroAndColors(dendro = geneTrees$par, colors = cbind(colors00$par, colors10$par, colors25$par, colors35$par),
                    groupLabels = c("unmerged", "merged at 0.90", "merged at 0.75", "merged at 0.65"), main = "Saccharomyces paradoxus",
                    dendroLabels = FALSE)
dev.off()

#### Parental Bipartite Module Alignments ####

# cer and par, merged at 0.25
is_grey <- (colors25$cer == "grey") | (colors25$par == "grey") 
sum(is_grey)/length(is_grey)
parents_color_counts <- makeBipartiteMatrix(.top_colors = colors25$cer[!is_grey], .bottom_colors = colors25$par[!is_grey])
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/AligningModules/Bipartite_Merge25_notPseudocolored.pdf",
    width = 4, height = 3)
plotweb(parents_color_counts$matrix, col.high = parents_color_counts$colors_high, col.low = parents_color_counts$colors_low, labsize = 1,
        high.lablength = 0, low.lablength = 0)
dev.off()

# corresponding random bipartite, for the plot
is_grey <- (random_colors25$cer == "grey") | (random_colors25$par == "grey") 
sum(is_grey)/length(is_grey)
parents_color_counts <- makeBipartiteMatrix(.top_colors = random_colors25$cer[!is_grey], .bottom_colors = random_colors25$par[!is_grey])
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/AligningModules/Random_Bipartite_Merge25_notPseudocolored.pdf",
    width = 4, height = 3)
plotweb(parents_color_counts$matrix, col.high = parents_color_counts$colors_high, col.low = parents_color_counts$colors_low, labsize = 1,
        high.lablength = 0, low.lablength = 0)
dev.off()

### Supplemental figure: other merge level parental bipartites
# I promise it's not worth trying to make these a single figure within R using recordPlot

# cer and par, unmerged
is_grey <- (colors00$cer == "grey") | (colors00$par == "grey") 
sum(is_grey)/length(is_grey) # ~1/3 of genes are not co-expressed in one species or the other (they also tend to have low variance across environments)
parents_color_counts <- makeBipartiteMatrix(.top_colors = colors00$cer[!is_grey], .bottom_colors = colors00$par[!is_grey])
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/Bipartite_Merge00_notPseudocolored.pdf",
    width = 4, height = 3)
plotweb(parents_color_counts$matrix, col.high = parents_color_counts$colors_high, col.low = parents_color_counts$colors_low, labsize = 1,
        high.lablength = 0, low.lablength = 0)
dev.off()

# corresponding random bipartite, for the plot
is_grey <- (random_colors00$cer == "grey") | (random_colors00$par == "grey") 
sum(is_grey)/length(is_grey) # more now b/c there's less overlap
parents_color_counts <- makeBipartiteMatrix(.top_colors = random_colors00$cer[!is_grey], .bottom_colors = random_colors00$par[!is_grey])
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/Random_Bipartite_Merge00_notPseudocolored.pdf",
    width = 4, height = 3)
plotweb(parents_color_counts$matrix, col.high = parents_color_counts$colors_high, col.low = parents_color_counts$colors_low, labsize = 1,
        high.lablength = 0, low.lablength = 0)
dev.off()

# cer and par, merged at 0.90
is_grey <- (colors10$cer == "grey") | (colors10$par == "grey") 
sum(is_grey)/length(is_grey) # ~1/3 of genes are not co-expressed in one species or the other (they also tend to have low variance across environments)
parents_color_counts <- makeBipartiteMatrix(.top_colors = colors10$cer[!is_grey], .bottom_colors = colors10$par[!is_grey])
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/Bipartite_Merge10_notPseudocolored.pdf",
    width = 4, height = 3)
plotweb(parents_color_counts$matrix, col.high = parents_color_counts$colors_high, col.low = parents_color_counts$colors_low, labsize = 1,
        high.lablength = 0, low.lablength = 0)
dev.off()

# corresponding random bipartite, for the plot
is_grey <- (random_colors10$cer == "grey") | (random_colors10$par == "grey") 
sum(is_grey)/length(is_grey)
parents_color_counts <- makeBipartiteMatrix(.top_colors = random_colors10$cer[!is_grey], .bottom_colors = random_colors10$par[!is_grey])
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/Random_Bipartite_Merge10_notPseudocolored.pdf",
    width = 4, height = 3)
plotweb(parents_color_counts$matrix, col.high = parents_color_counts$colors_high, col.low = parents_color_counts$colors_low, labsize = 1,
        high.lablength = 0, low.lablength = 0)
dev.off()

# cer and par, merged at 0.75
is_grey <- (colors25$cer == "grey") | (colors25$par == "grey") 
sum(is_grey)/length(is_grey) # ~1/3 of genes are not co-expressed in one species or the other (they also tend to have low variance across environments)
parents_color_counts <- makeBipartiteMatrix(.top_colors = colors25$cer[!is_grey], .bottom_colors = colors25$par[!is_grey])
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/Bipartite_Merge25_notPseudocolored.pdf",
    width = 4, height = 3)
plotweb(parents_color_counts$matrix, col.high = parents_color_counts$colors_high, col.low = parents_color_counts$colors_low, labsize = 1,
        high.lablength = 0, low.lablength = 0)
dev.off()

# corresponding random bipartite, for the plot
is_grey <- (random_colors25$cer == "grey") | (random_colors25$par == "grey") 
sum(is_grey)/length(is_grey) # more now b/c there's less overlap
parents_color_counts <- makeBipartiteMatrix(.top_colors = random_colors25$cer[!is_grey], .bottom_colors = random_colors25$par[!is_grey])
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/Random_Bipartite_Merge25_notPseudocolored.pdf",
    width = 4, height = 3)
plotweb(parents_color_counts$matrix, col.high = parents_color_counts$colors_high, col.low = parents_color_counts$colors_low, labsize = 1,
        high.lablength = 0, low.lablength = 0)
dev.off()

# cer and par, merged at 0.65
is_grey <- (colors35$cer == "grey") | (colors35$par == "grey") 
sum(is_grey)/length(is_grey) # ~1/3 of genes are not co-expressed in one species or the other (they also tend to have low variance across environments)
parents_color_counts <- makeBipartiteMatrix(.top_colors = colors35$cer[!is_grey], .bottom_colors = colors35$par[!is_grey])
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/Bipartite_Merge35_notPseudocolored.pdf",
    width = 4, height = 3)
plotweb(parents_color_counts$matrix, col.high = parents_color_counts$colors_high, col.low = parents_color_counts$colors_low, labsize = 1,
        high.lablength = 0, low.lablength = 0)
dev.off()

# corresponding random bipartite, for the plot
is_grey <- (random_colors35$cer == "grey") | (random_colors35$par == "grey") 
sum(is_grey)/length(is_grey) # more now b/c there's less overlap
parents_color_counts <- makeBipartiteMatrix(.top_colors = random_colors35$cer[!is_grey], .bottom_colors = random_colors35$par[!is_grey])
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/Random_Bipartite_Merge35_notPseudocolored.pdf",
    width = 4, height = 3)
plotweb(parents_color_counts$matrix, col.high = parents_color_counts$colors_high, col.low = parents_color_counts$colors_low, labsize = 1,
        high.lablength = 0, low.lablength = 0)
dev.off()

#### hybrid bipartites ####
# hyc and hyp, merged at 90%
is_grey <- (colors10$hyc == "grey") | (colors10$hyp == "grey") 
sum(is_grey)/length(is_grey)
hybrid_color_counts <- makeBipartiteMatrix(.top_colors = colors10$hyc[!is_grey], .bottom_colors = colors10$hyp[!is_grey])
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/Bipartite_HycHyp10_notPseudocolored.pdf",
    width = 4, height = 3)
plotweb(hybrid_color_counts$matrix, col.high = hybrid_color_counts$colors_high, col.low = hybrid_color_counts$colors_low, labsize = 1,
        high.lablength = 0, low.lablength = 0)
dev.off()

#### Are networks robust to leaving one environment out? ####

# how similar are module definitions when you remove any of the environments?

load("data_files/LeaveOneOut.RData")

# exploratory bipartites for each leave-one-out parental network versus
# the full network

# cer, no LowN
noLowN_cer <- makeBipartiteMatrix(.top_colors = networks_noLowN$colors$cer, .bottom_colors = colors00$cer)
plotweb(noLowN_cer$matrix, col.high = noLowN_cer$colors_high, col.low = noLowN_cer$colors_low, labsize = 1,
        high.lablength = 0, low.lablength = 0)
# par, no LowN
noLowN_par <- makeBipartiteMatrix(.top_colors = networks_noLowN$colors$par, .bottom_colors = colors00$par)
plotweb(noLowN_par$matrix, col.high = noLowN_par$colors_high, col.low = noLowN_par$colors_low, labsize = 1,
        high.lablength = 0, low.lablength = 0)

# cer, no HAP4
noHAP4_cer <- makeBipartiteMatrix(.top_colors = networks_noHAP4$colors$cer, .bottom_colors = colors00$cer)
plotweb(noHAP4_cer$matrix, col.high = noHAP4_cer$colors_high, col.low = noHAP4_cer$colors_low, labsize = 1,
        high.lablength = 0, low.lablength = 0)
# par, no HAP4
noHAP4_par <- makeBipartiteMatrix(.top_colors = networks_noHAP4$colors$par, .bottom_colors = colors00$par)
plotweb(noHAP4_par$matrix, col.high = noHAP4_par$colors_high, col.low = noHAP4_par$colors_low, labsize = 1,
        high.lablength = 0, low.lablength = 0)

# cer, no LowPi
noLowPi_cer <- makeBipartiteMatrix(.top_colors = networks_noLowPi$colors$cer, .bottom_colors = colors00$cer)
plotweb(noLowPi_cer$matrix, col.high = noLowPi_cer$colors_high, col.low = noLowPi_cer$colors_low, labsize = 1,
        high.lablength = 0, low.lablength = 0)
# par, no LowPi
noLowPi_par <- makeBipartiteMatrix(.top_colors = networks_noLowPi$colors$par, .bottom_colors = colors00$par)
plotweb(noLowPi_par$matrix, col.high = noLowPi_par$colors_high, col.low = noLowPi_par$colors_low, labsize = 1,
        high.lablength = 0, low.lablength = 0)

# cer, no CC
noCC_cer <- makeBipartiteMatrix(.top_colors = networks_noCC$colors$cer, .bottom_colors = colors00$cer)
plotweb(noCC_cer$matrix, col.high = noCC_cer$colors_high, col.low = noCC_cer$colors_low, labsize = 1,
        high.lablength = 0, low.lablength = 0)
# par, no CC
noCC_par <- makeBipartiteMatrix(.top_colors = networks_noCC$colors$par, .bottom_colors = colors00$par)
plotweb(noCC_par$matrix, col.high = noCC_par$colors_high, col.low = noCC_par$colors_low, labsize = 1,
        high.lablength = 0, low.lablength = 0)

# cer, no Temp
noTemp_cer <- makeBipartiteMatrix(.top_colors = networks_noTemp$colors$cer, .bottom_colors = colors00$cer)
plotweb(noTemp_cer$matrix, col.high = noTemp_cer$colors_high, col.low = noTemp_cer$colors_low, labsize = 1,
        high.lablength = 0, low.lablength = 0)
# par, no Temp
noTemp_par <- makeBipartiteMatrix(.top_colors = networks_noTemp$colors$par, .bottom_colors = colors00$par)
plotweb(noTemp_par$matrix, col.high = noTemp_par$colors_high, col.low = noTemp_par$colors_low, labsize = 1,
        high.lablength = 0, low.lablength = 0)

# modules appear to be larger in the leave-one-outs, but largely the same genes still grouping together

# TODO: simple supplementary figure where x axis is catagorical,
# whichever environment was excluded and y axis is % module conservation
# (# of genes with sig pvalue in leaveoneoutdf)
# and another where y axis is number of modules

#### P-value distributions of module conservation tests ####

### Fig inset: example random module distribution with observed, significant value
callConservedConnectionRaw <- function(.color_name1, .color_name2, .colors1, .colors2, .n_permutations_per_species = 50000) {
  if (length(.colors1) != length(.colors2)) {
    stop("color vectors are not the same length!\n")
  }
  n_connections <- sum(.colors1 == .color_name1 & .colors2 == .color_name2)
  # helper function
  scrambleConnections <- function(.steady_name, .scramble_name, .steady_colors, .scramble_colors) {
    scrambled <- sample(.scramble_colors, size = length(.scramble_colors), replace = FALSE)
    n_connections <- sum(.steady_colors == .steady_name & scrambled == .scramble_name)
    return(n_connections)
  }
  # scrambling first set of colors
  null_connections1 <- sapply(c(1:.n_permutations_per_species), \(dummy) {
    output <- scrambleConnections(.steady_name = .color_name1, 
                                  .scramble_name = .color_name2,
                                  .steady_colors = .colors1,
                                  .scramble_colors = .colors2)
    return(output)
  })
  # scrambling second set of colors
  null_connections2 <- sapply(c(1:.n_permutations_per_species), \(dummy) {
    output <- scrambleConnections(.steady_name = .color_name2, 
                                  .scramble_name = .color_name1,
                                  .steady_colors = .colors2,
                                  .scramble_colors = .colors1)
    return(output)
  })
  null_connections <- c(null_connections1, null_connections2)
  return(list("null" = null_connections,
              "observed" = n_connections))
}
# example ccm
ccm_color_pair <- moduledf25 |> filter(is_CCM) |>
  slice_sample(n = 1) |> 
  select(cer_color, par_color, CCM_color)
nullmodsizes <- callConservedConnectionRaw(.color_name1 = ccm_color_pair$cer_color,
                                           .color_name2 = ccm_color_pair$par_color,
                                           .colors1 = colors25$cer,
                                           .colors2 = colors25$par)

pdf("../../aligning_the_molecular_phenotype/paper_figures/AligningModules/example_permutation_test.pdf",
    height = 4, width = 5)
hist(nullmodsizes$null, xlim = c(0, nullmodsizes$observed + 10),
     xlab = "module size (number of genes)", 
     main = paste("Example module permutation test,\n",
                  ccm_color_pair$CCM_color, "module"))
abline(v = nullmodsizes$observed, col = "red")
dev.off()

# example random
random_color_pair <- random_moduledf25 |>
  slice_sample(n = 1) |> 
  select(cer_color, par_color, CCM_color)
nullmodsizes <- callConservedConnectionRaw(.color_name1 = random_color_pair$cer_color,
                                           .color_name2 = random_color_pair$par_color,
                                           .colors1 = random_colors25$cer,
                                           .colors2 = random_colors25$par)

pdf("../../aligning_the_molecular_phenotype/paper_figures/AligningModules/example_random_permutation_test.pdf",
    height = 4, width = 5)
hist(nullmodsizes$null,
     xlab = "module size (number of genes)", 
     main = paste("Example module permutation test,\nrandom module"))
abline(v = nullmodsizes$observed, col = "red")
dev.off()

# all modules pvalue distribution, merged at 0.75
ccm_p_thresh <- 0.05/nrow(CCMdf25)
plotvec <- CCMdf25 |> filter(cer_color != "grey" & par_color != "grey") |>
  mutate(prop_null_smaller = (1-p_value)) |> select(prop_null_smaller) |> pull() |> 
  cut(c(-Inf, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, (1-ccm_p_thresh), Inf))
max_y <- plotvec |> table() |> max()
pdf("../../aligning_the_molecular_phenotype/paper_figures/AligningModules/CCMpvalues25.pdf",
    width = 4, height = 4)
hist(as.numeric(plotvec), breaks = 0:(nlevels(plotvec)), xaxt = 'n', 
     xlab = "fraction of permuted modules \nsmaller than observed",
     ylab = "count", main = "Observed modules",
     ylim = c(0, max_y))
abline(v = 10, col = "gold")
axis(1, at = c(0.5, 2.5, 4.5, 6.5, 10.5), 
     labels = c(0, 0.2, 0.4, 0.6, paste0(">", round(1 - ccm_p_thresh, 3))))
dev.off()

# Corresponding random modules
ccm_p_thresh <- 0.05/nrow(CCMdf25)
plotvec <- random_CCMdf25 |> filter(cer_color != "grey" & par_color != "grey") |>
  mutate(prop_null_smaller = (1-p_value)) |> select(prop_null_smaller) |> pull() |> 
  cut(c(-Inf, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, (1-ccm_p_thresh), Inf))
pdf("../../aligning_the_molecular_phenotype/paper_figures/AligningModules/Random_CCMpvalues25.pdf",
    width = 4, height = 4)
hist(as.numeric(plotvec), breaks = 0:(nlevels(plotvec)), xaxt = 'n', 
     xlab = "fraction of permuted modules \nsmaller than observed",
     ylab = "count", main = "Random modules",
     ylim = c(0, max_y))
abline(v = 10, col = "gold")
axis(1, at = c(0.5, 2.5, 4.5, 6.5, 10.5), 
     labels = c(0, 0.2, 0.4, 0.6, paste0(">", round(1 - ccm_p_thresh, 3))))
dev.off()

### Supplement: other merge thresholds

# all modules pvalue distribution, unmerged
ccm_p_thresh <- 0.05/nrow(CCMdf00)
plotvec <- CCMdf00 |> filter(cer_color != "grey" & par_color != "grey") |>
  mutate(prop_null_smaller = (1-p_value)) |> select(prop_null_smaller) |> pull() |> 
  cut(c(-Inf, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, (1-ccm_p_thresh), Inf))
max_y <- plotvec |> table() |> max()
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/CCMpvalues00.pdf",
    width = 4, height = 4)
hist(as.numeric(plotvec), breaks = 0:(nlevels(plotvec)), xaxt = 'n', 
     xlab = "fraction of permuted modules \nsmaller than observed",
     ylab = "count", main = "Observed modules",
     ylim = c(0, max_y))
abline(v = 10, col = "gold")
axis(1, at = c(0.5, 2.5, 4.5, 6.5, 10.5), 
     labels = c(0, 0.2, 0.4, 0.6, paste0(">", round(1 - ccm_p_thresh, 4))))
dev.off()

# Corresponding random modules
ccm_p_thresh <- 0.05/nrow(CCMdf00)
plotvec <- random_CCMdf00 |> filter(cer_color != "grey" & par_color != "grey") |>
  mutate(prop_null_smaller = (1-p_value)) |> select(prop_null_smaller) |> pull() |> 
  cut(c(-Inf, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, (1-ccm_p_thresh), Inf))
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/Random_CCMpvalues00.pdf",
    width = 4, height = 4)
hist(as.numeric(plotvec), breaks = 0:(nlevels(plotvec)), xaxt = 'n', 
     xlab = "fraction of permuted modules \nsmaller than observed",
     ylab = "count", main = "Random modules",
     ylim = c(0, max_y))
abline(v = 10, col = "gold")
axis(1, at = c(0.5, 2.5, 4.5, 6.5, 10.5), 
     labels = c(0, 0.2, 0.4, 0.6, paste0(">", round(1 - ccm_p_thresh, 4))))
dev.off()

# merged at 0.90
ccm_p_thresh <- 0.05/nrow(CCMdf10)
plotvec <- CCMdf10 |> filter(cer_color != "grey" & par_color != "grey") |>
  mutate(prop_null_smaller = (1-p_value)) |> select(prop_null_smaller) |> pull() |> 
  cut(c(-Inf, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, (1-ccm_p_thresh), Inf))
max_y <- plotvec |> table() |> max()
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/CCMpvalues10.pdf",
    width = 4, height = 4)
hist(as.numeric(plotvec), breaks = 0:(nlevels(plotvec)), xaxt = 'n', 
     xlab = "fraction of permuted modules \nsmaller than observed",
     ylab = "count", main = "Observed modules",
     ylim = c(0, max_y))
abline(v = 10, col = "gold")
axis(1, at = c(0.5, 2.5, 4.5, 6.5, 10.5), 
     labels = c(0, 0.2, 0.4, 0.6, paste0(">", round(1 - ccm_p_thresh, 4))))
dev.off()

# Corresponding random modules
ccm_p_thresh <- 0.05/nrow(CCMdf10)
plotvec <- random_CCMdf10 |> filter(cer_color != "grey" & par_color != "grey") |>
  mutate(prop_null_smaller = (1-p_value)) |> select(prop_null_smaller) |> pull() |> 
  cut(c(-Inf, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, (1-ccm_p_thresh), Inf))
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/Random_CCMpvalues10.pdf",
    width = 4, height = 4)
hist(as.numeric(plotvec), breaks = 0:(nlevels(plotvec)), xaxt = 'n', 
     xlab = "fraction of permuted modules \nsmaller than observed",
     ylab = "count", main = "Random modules",
     ylim = c(0, max_y))
abline(v = 10, col = "gold")
axis(1, at = c(0.5, 2.5, 4.5, 6.5, 10.5), 
     labels = c(0, 0.2, 0.4, 0.6, paste0(">", round(1 - ccm_p_thresh, 4))))
dev.off()

# all modules pvalue distribution, merged at 0.75
ccm_p_thresh <- 0.05/nrow(CCMdf25)
plotvec <- CCMdf25 |> filter(cer_color != "grey" & par_color != "grey") |>
  mutate(prop_null_smaller = (1-p_value)) |> select(prop_null_smaller) |> pull() |> 
  cut(c(-Inf, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, (1-ccm_p_thresh), Inf))
max_y <- plotvec |> table() |> max()
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/CCMpvalues25.pdf",
    width = 4, height = 4)
hist(as.numeric(plotvec), breaks = 0:(nlevels(plotvec)), xaxt = 'n', 
     xlab = "fraction of permuted modules \nsmaller than observed",
     ylab = "count", main = "Observed modules",
     ylim = c(0, max_y))
abline(v = 10, col = "gold")
axis(1, at = c(0.5, 2.5, 4.5, 6.5, 10.5), 
     labels = c(0, 0.2, 0.4, 0.6, paste0(">", round(1 - ccm_p_thresh, 3))))
dev.off()

# Corresponding random modules
ccm_p_thresh <- 0.05/nrow(CCMdf25)
plotvec <- random_CCMdf25 |> filter(cer_color != "grey" & par_color != "grey") |>
  mutate(prop_null_smaller = (1-p_value)) |> select(prop_null_smaller) |> pull() |> 
  cut(c(-Inf, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, (1-ccm_p_thresh), Inf))
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/Random_CCMpvalues25.pdf",
    width = 4, height = 4)
hist(as.numeric(plotvec), breaks = 0:(nlevels(plotvec)), xaxt = 'n', 
     xlab = "fraction of permuted modules \nsmaller than observed",
     ylab = "count", main = "Random modules",
     ylim = c(0, max_y))
abline(v = 10, col = "gold")
axis(1, at = c(0.5, 2.5, 4.5, 6.5, 10.5), 
     labels = c(0, 0.2, 0.4, 0.6, paste0(">", round(1 - ccm_p_thresh, 3))))
dev.off()

# all modules pvalue distribution, merged at 0.65
ccm_p_thresh <- 0.05/nrow(CCMdf35)
plotvec <- CCMdf35 |> filter(cer_color != "grey" & par_color != "grey") |>
  mutate(prop_null_smaller = (1-p_value)) |> select(prop_null_smaller) |> pull() |> 
  cut(c(-Inf, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, (1-ccm_p_thresh), Inf))
max_y <- plotvec |> table() |> max()
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/CCMpvalues35.pdf",
    width = 4, height = 4)
hist(as.numeric(plotvec), breaks = 0:(nlevels(plotvec)), xaxt = 'n', 
     xlab = "fraction of permuted modules \nsmaller than observed",
     ylab = "count", main = "Observed modules",
     ylim = c(0, max_y))
abline(v = 10, col = "gold")
axis(1, at = c(0.5, 2.5, 4.5, 6.5, 10.5), 
     labels = c(0, 0.2, 0.4, 0.6, paste0(">", round(1 - ccm_p_thresh, 3))))
dev.off()

# Corresponding random modules
ccm_p_thresh <- 0.05/nrow(CCMdf35)
plotvec <- random_CCMdf35 |> filter(cer_color != "grey" & par_color != "grey") |>
  mutate(prop_null_smaller = (1-p_value)) |> select(prop_null_smaller) |> pull() |> 
  cut(c(-Inf, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, (1-ccm_p_thresh), Inf))
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/Random_CCMpvalues35.pdf",
    width = 4, height = 4)
hist(as.numeric(plotvec), breaks = 0:(nlevels(plotvec)), xaxt = 'n', 
     xlab = "fraction of permuted modules \nsmaller than observed",
     ylab = "count", main = "Random modules",
     ylim = c(0, max_y))
abline(v = 10, col = "gold")
axis(1, at = c(0.5, 2.5, 4.5, 6.5, 10.5), 
     labels = c(0, 0.2, 0.4, 0.6, paste0(">", round(1 - ccm_p_thresh, 3))))
dev.off()

#### Supplement: CCM membership cannot predict single gene expression divergence (as measured by abs(LFC)) ####
plotdf <- module_genedf25 |> select(gene_name, cer_color, par_color, is_CCM, max_lfc, coexpressed)
# converting LFC from natural log to log2 fold change (because that's what's typically used in gene expression studies)
plotdf$max_lfc <- plotdf$max_lfc/log(2)
group_means <- plotdf |> group_by(coexpressed) |> summarise(mean_lfc = mean(max_lfc, na.rm = TRUE)) |> select(mean_lfc)
ttest1 <- compare_means(max_lfc ~ coexpressed, data = plotdf, method = "t.test")
ttest1
ttest1 <- ttest1 |> filter(group1 == "conserved co-expressed")
ttest1
ttest1$y.position <- c(5:8)
ttest1$p.format <- round(as.numeric(ttest1$p.format), 3)

# check extreme values
plotdf |> filter(max_lfc > 10)

p <- ggplot(plotdf, aes(x = coexpressed, y = max_lfc)) + 
  geom_jitter(aes(color = coexpressed)) + 
  xlab("Co-expression partners between species") +
  ylab("Expression change between species\n (log fold change)") +
  scale_color_discrete(type = c("gold",
                                "mediumseagreen",
                                "orange1",
                                "blue2",
                                "grey"),
                       limits = c("conserved co-expressed",
                                  "diverged co-expressed",
                                  "S. cerevisiae co-expressed",
                                  "S. paradoxus co-expressed",
                                  "never co-expressed")) +
  scale_x_discrete(limits = c("conserved co-expressed",
                                  "diverged co-expressed",
                                  "S. cerevisiae co-expressed",
                                  "S. paradoxus co-expressed",
                                  "never co-expressed"),
                       labels = c("conserved\n co-expressed",
                                  "diverged \nco-expressed",
                                  "S. cerevisiae\n co-expressed",
                                  "S. paradoxus\n co-expressed",
                                  "never\n co-expressed")) +
  xlab("") +
  ylab("| Log2 fold change |\nmax across 4 environments") +
  theme_classic() +
  theme(legend.position = "none") +
  geom_segment(data = group_means,
               aes(x = seq(from = 0.8, to = 4.8, by = 1), 
                   xend = seq(from = 1.2, to = 5.2, by = 1), 
                   y = as.numeric(mean_lfc), yend = as.numeric(mean_lfc)), color = "red") +
  stat_pvalue_manual(ttest1, label = "t.test p = {p.format}", tip.length = 0.01) +
  scale_y_continuous(position = "right", limits = c(0, 10.5)) + 
  # manually adding vertical arrowheads to indicate the 4 outliers with their names and values next to them
  geom_segment(aes(x = 1, xend = 1, y = 10, yend = 10.5), color = "gold", arrow = arrow(angle = 30, length = unit(0.05, "inches"), type = "closed")) +
  annotate("text", x = 1, y = 9.75, color = "gold", label = "YPL274W, 12.3", size = 2) +
  geom_segment(aes(x = 4, xend = 4, y = 10, yend = 10.5), color = "blue2", arrow = arrow(angle = 30, length = unit(0.05, "inches"), type = "closed")) +
  annotate("text", x = 4, y = 9.75, color = "blue2", label = "YDL037C, 57.4", size = 2) +
  geom_segment(aes(x = 4, xend = 4, y = 9, yend = 9.5), color = "blue2", arrow = arrow(angle = 30, length = unit(0.05, "inches"), type = "closed")) +
  annotate("text", x = 4, y = 8.75, color = "blue2", label = "YFL051C, 10.7", size = 2) +
  geom_segment(aes(x = 3, xend = 3, y = 10, yend = 10.5), color = "orange1", arrow = arrow(angle = 30, length = unit(0.05, "inches"), type = "closed")) +
  annotate("text", x = 3, y = 9.75, color = "orange1", label = "YPR199C, 54.7", size = 2) +
  geom_segment(aes(x = 3, xend = 3, y = 9, yend = 9.5), color = "orange1", arrow = arrow(angle = 30, length = unit(0.05, "inches"), type = "closed")) +
  annotate("text", x = 3, y = 8.75, color = "orange1", label = "YIR041W, 58.6", size = 2)
### Supplement: being part of a conserved, co-expressed module is not sufficient to explain expression divergence
# the easiest way to make sure I'm getting the right means and pvalue put on the ggplot is to do both versions of t-test
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/LFC_vs_module_conservation.pdf",
#     width = 6, height = 3)
p
# dev.off()

# Checking LFC in other module designations
getTtestResult <- function(.df, .mod_gdf) {
  df0 <- select(ungroup(.df), gene_name, max_lfc) |> 
    left_join(select(.mod_gdf, gene_name, cer_color, par_color, is_CCM, coexpressed),
              by = c("gene_name")) |> drop_na()
  ttest_res <- compare_means(max_lfc ~ coexpressed, data = df0, method = "t.test") |> 
    filter(group1 == "conserved co-expressed")
  ttest_res$mean1 <- sapply(ttest_res$group1, \(x) mean(df0$max_lfc[df0$coexpressed == x]))
  ttest_res$mean2 <- sapply(ttest_res$group2, \(x) mean(df0$max_lfc[df0$coexpressed == x]))
  return(ttest_res)
}
# Repeating t-tests with other three module definitions (colors00, colors10, and colors35) to make
# sure lack of connection btwn co-expression and expression divergence is robust to module definitions
# Merge 00
getTtestResult(plotdf, module_genedf00)
# Merge 10
getTtestResult(plotdf, module_genedf10)
# Merge 25
getTtestResult(plotdf, module_genedf25)
# Merge 35
getTtestResult(plotdf, module_genedf35)
# same pattern whichever module designation is used:
# divergent co-expressed genes maybe have slightly higher lfc than conserved co-expressed,
# but otherwise co-expression status cannot predict whether a gene has diverged in overall expression level

# Inset: example of how LFC is calculated
library(MASS, include.only = "glm.nb")
gene_idx <- "YEL060C" # arbitrary CCM gene we looked up the max LFC magnitude for
test_genedf <- bind_cols(tibble(expr = spcts$CC[gene_idx,]), spinfo$CC)
spmod <- glm.nb(expr ~ time_point_num + allele, data = test_genedf, link = log, init.theta = 7)
sprow <- summary(spmod)$coefficients["allelecer",, drop = FALSE]
sprow # positive effect size is counts biased for cer, negative is counts biased for par
slope <- 1/exp(sprow[1])
sprow[1]/log(2) # log2 fold change (model is natural log fold change)
intercept <- summary(spmod)$coefficients["(Intercept)",, drop = FALSE][1]
plotdf <- test_genedf %>% mutate(fittedvals = spmod$fitted.values)%>% select(allele, condition, well_flask_ID, expr, fittedvals, genotype) %>% 
  pivot_longer(cols = c(expr, fittedvals)) %>% 
  mutate(non_unique_sample_name = paste(condition, well_flask_ID)) %>% 
  pivot_wider(id_cols = c(non_unique_sample_name, name, genotype), names_from = allele, values_from = value)
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/LFC_vs_module_conservation_insert.pdf",
    width = 4, height = 3)
ggplot(plotdf, aes(x = cer, y = par)) + geom_point(aes(color = name)) + geom_abline(color = "gold", slope = slope, intercept = intercept) + 
  geom_abline(color = "midnightblue", slope = 1, intercept = 0) + 
  xlim(c(0, max(select(plotdf, cer, par)))) + ylim(c(0, max(select(plotdf, cer, par)))) + 
  theme_classic() + ggtitle("Generalized linear model") +
  annotate("text", x = 400, y = 150, label = "|LFC| = |-2.35| = 2.35\nslope = 1/2^(-2.35) = 5.05\npvalue < 1e-50", size = 2.5) +
  annotate("text", x = 500, y = 450, label = "y = x", size = 2.5, color = "midnightblue") +
  annotate("text", x = 170, y = 630, label = "model", size = 2.5, color = "gold") +
  xlab("Expression in S. cerevisiae (cpm)") +
  ylab("Expression in S. paradoxus (cpm)") +
  theme(legend.title = element_blank()) +
  scale_color_discrete(labels = c("data", "fitted values"))
dev.off()

### Same plot but looking at mean effect size, positive (higher in par) or negative (higher in cer)
plotdf <- left_join(filter(spaldf, coefficient == "species"),
                    select(module_genedf25, gene_name, cer_color, par_color, is_CCM, coexpressed),
                    by = "gene_name", relationship = "many-to-one") |> drop_na()
plotdf <- plotdf |> group_by(gene_name, cer_color, par_color, is_CCM, coexpressed) |> 
  summarise(mean_lfc = mean(effect_size))
plotdf$mean_lfc <- plotdf$mean_lfc/log(2)
group_means <- plotdf |> group_by(coexpressed) |> summarise(mean_lfc = mean(mean_lfc, na.rm = TRUE)) |> select(mean_lfc)

# checking extreme values
plotdf |> filter(abs(mean_lfc) > 10)

# first the plot that includes direction of lfc (higher in cer, or higher in par)
# no statistical tests for this one
p <- ggplot(plotdf, aes(x = coexpressed, y = mean_lfc)) + 
  geom_hline(yintercept = 0, color = "midnightblue") +
  geom_jitter(aes(color = coexpressed)) + 
  geom_text(data = filter(plotdf, abs(mean_lfc) > 10), aes(label = gene_name, color = coexpressed),
            size = 2, nudge_y = 0.75) +
  xlab("Co-expression partners between species") +
  ylab("Expression change between species\n (log fold change)") +
  scale_color_discrete(type = c("gold",
                                "mediumseagreen",
                                "orange1",
                                "blue2",
                                "grey"), 
                       limits = c("conserved co-expressed",
                                  "diverged co-expressed",
                                  "S. cerevisiae co-expressed",
                                  "S. paradoxus co-expressed",
                                  "never co-expressed")) +
  scale_x_discrete(limits = c("conserved co-expressed",
                              "diverged co-expressed",
                              "S. cerevisiae co-expressed",
                              "S. paradoxus co-expressed",
                              "never co-expressed"),
                   labels = c("conserved\n co-expressed",
                              "diverged \nco-expressed",
                              "S. cerevisiae\n co-expressed",
                              "S. paradoxus\n co-expressed",
                              "never\n co-expressed")) +
  xlab("") +
  ylab("Log fold change\nmean across 4 environments\n<-higher in S. cerevisiae, higher in S. paradoxus ->") +
  theme_classic() +
  theme(legend.position = "none") +
  geom_segment(data = group_means,
               aes(x = seq(from = 0.8, to = 4.8, by = 1), 
                   xend = seq(from = 1.2, to = 5.2, by = 1), 
                   y = as.numeric(mean_lfc), yend = as.numeric(mean_lfc)), color = "red") +
  scale_y_continuous(position = "right")

# pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/mean_lfc.pdf",
#     width = 7, height = 4)
p
# dev.off()

# now looking at just magnitude of mean LFC, with statistical tests
plotdf$abs_mean_lfc <- abs(plotdf$mean_lfc)
ttest1 <- compare_means(abs_mean_lfc ~ coexpressed, data = plotdf, method = "t.test")
ttest1
ttest1 <- ttest1 |> filter(group1 == "conserved co-expressed")
ttest1
ttest1$y.position <- seq(from = 7, to = 13, by = 2)
ttest1$p.format <- round(as.numeric(ttest1$p.format), 5)

group_means <- plotdf |> group_by(coexpressed) |> summarise(abs_mean_lfc = mean(abs_mean_lfc, na.rm = TRUE)) |> select(abs_mean_lfc)

p <- ggplot(plotdf, aes(x = coexpressed, y = abs_mean_lfc)) + 
  geom_jitter(aes(color = coexpressed)) + 
  geom_text(data = filter(plotdf, abs(mean_lfc) > 10), aes(label = gene_name, color = coexpressed),
            size = 2, nudge_y = 0.75) +
  xlab("Co-expression partners between species") +
  ylab("Expression change between species\n (log fold change)") +
  scale_color_discrete(type = c("gold",
                                "mediumseagreen",
                                "orange1",
                                "blue2",
                                "grey"), 
                       limits = c("conserved co-expressed",
                                  "diverged co-expressed",
                                  "S. cerevisiae co-expressed",
                                  "S. paradoxus co-expressed",
                                  "never co-expressed")) +
  scale_x_discrete(limits = c("conserved co-expressed",
                              "diverged co-expressed",
                              "S. cerevisiae co-expressed",
                              "S. paradoxus co-expressed",
                              "never co-expressed"),
                   labels = c("conserved\n co-expressed",
                              "diverged \nco-expressed",
                              "S. cerevisiae\n co-expressed",
                              "S. paradoxus\n co-expressed",
                              "never\n co-expressed")) +
  xlab("") +
  ylab("|Log fold change|\n(absolute value of mean \nacross 4 environments)") +
  theme_classic() +
  theme(legend.position = "none") +
  geom_segment(data = group_means,
               aes(x = seq(from = 0.8, to = 4.8, by = 1), 
                   xend = seq(from = 1.2, to = 5.2, by = 1), 
                   y = as.numeric(abs_mean_lfc), yend = as.numeric(abs_mean_lfc)), color = "red") +
  stat_pvalue_manual(ttest1, label = "t.test p = {p.format}", tip.length = 0.01) +
  scale_y_continuous(position = "right")

# pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/mean_lfc_magnitude.pdf",
#     width = 7, height = 4)
# p
# dev.off()

# Checking LFC in other module designations
getTtestResult <- function(.df, .mod_gdf) {
  df0 <- select(ungroup(.df), gene_name, abs_mean_lfc) |> 
    left_join(select(.mod_gdf, gene_name, cer_color, par_color, is_CCM, coexpressed),
              by = c("gene_name")) |> drop_na()
  ttest_res <- compare_means(abs_mean_lfc ~ coexpressed, data = df0, method = "t.test") |> 
    filter(group1 == "conserved co-expressed")
  ttest_res$mean1 <- sapply(ttest_res$group1, \(x) mean(df0$abs_mean_lfc[df0$coexpressed == x]))
  ttest_res$mean2 <- sapply(ttest_res$group2, \(x) mean(df0$abs_mean_lfc[df0$coexpressed == x]))
  return(ttest_res)
}
# Repeating t-tests with other three module definitions (colors00, colors10, and colors35) to make
# sure lack of connection btwn co-expression and expression divergence is robust to module definitions
# Merge 00
getTtestResult(plotdf, module_genedf00)
# Merge 10
getTtestResult(plotdf, module_genedf10)
# Merge 25
getTtestResult(plotdf, module_genedf25)
# Merge 35
getTtestResult(plotdf, module_genedf35)

# when using mean lfc instead of max lfc, diverged is still slightly higher mean than
# conserved co-expressed, but explaining very little variation

################################### Archive ######################################
# ### Bipartites of parental ccms preservation in the hybrid ####
# # @input: .modules = vector of length nGenes giving the module name of each gene
# #         .colors = vector of length nGenes giving the color (single species module) of each gene
# #         .colorseq = order of length unique(.colors)
# makeModColorBipartite <- function(.modules, .colors, .colorseq = c("grey", setdiff(unique(.colors), "grey"))) {
#   color_counts <- makeBipartiteMatrix(.top_colors = .modules,
#                                       .bottom_colors = .colors, .simplify_thresh =  10)
#   color_counts$colors_high
#   # assigning 5 colors to the modules
#   new_colors <- sapply(color_counts$colors_high, \(m) {
#     if (grepl("ccm", m)) {
#       return("gold")
#     }
#     if (grepl("div", m)) {
#       return("mediumseagreen")
#     }
#     if (grepl("cer", m)) {
#       return("orange1")
#     }
#     if (grepl("par", m)) {
#       return("blue2")
#     }
#     if (grepl("nev", m)) {
#       return("grey")
#     }
#   }) |> unlist() |> as.character()
#   mat <- color_counts$matrix
#   modorder <- c(colnames(mat)[grep("nev", colnames(mat))],
#                 colnames(mat)[grep("par", colnames(mat))],
#                 colnames(mat)[grep("cer", colnames(mat))],
#                 colnames(mat)[grep("div", colnames(mat))],
#                 colnames(mat)[grep("ccm", colnames(mat))])
#   return(plotweb(mat, col.high = new_colors,
#                  col.low = rownames(mat), empty = FALSE,
#                  low.lablength = 0,
#                  labsize = 1,
#                  text.rot = 90,
#                  sequence = list(seq.high = modorder,
#                                  seq.low = .colorseq)))
# }
# 
# # parental modules vs hyc
# pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/Bipartite_parentalModulesVsHyc10.pdf", 
#     width = 4, height = 3)
# makeModColorBipartite(module_genedf10$module_name, colors10$hyc)
# dev.off()
# 
# # parental modules vs hyp
# pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/Bipartite_parentalModulesVsHyp10.pdf", 
#     width = 4, height = 3)
# makeModColorBipartite(module_genedf10$module_name, colors10$hyp)
# dev.off()
# 
# # parental modules vs hyb unmerged
# makeModColorBipartite(module_genedf00$module_name, colors00$hyb)
# # parental modules vs hyc unmerged
# makeModColorBipartite(module_genedf00$module_name, colors00$hyc)
# # parental modules vs hyp unmerged
# makeModColorBipartite(module_genedf00$module_name, colors00$hyp)

# #### Are not co-expressed (grey) genes just less variable between environments? ####
# # rationale: if not co-expressed genes have flatter expression vectors,
# # they would tend to have lower correlations with all other genes 
# # and therefore not be placed in modules
# 
# ### Mean vs Var for genes in 4 groups (1 plot per species): 
# # 1) conserved co-expressed, significant configuration model pvalue
# # 2) co-expressed, non significant config model
# # 3) co-expressed only in one species (all of these turned out to have non significant pvalues)
# # 4) never co-expressed
# 
# # Merge 10
# plotdf <- module_genedf10 |> select(gene_name, cer_color, par_color, is_CCM, coexpressed)
# plotdf$mean_cer <- colMeans(counts_all2$cer)
# plotdf$var_cer <- apply(counts_all2$cer, 2, var)
# plotdf$mean_par <- colMeans(counts_all2$par)
# plotdf$var_par <- apply(counts_all2$par, 2, var)
# # plotting
# pcer_meanvar <- ggplot(plotdf, aes(x = log2(mean_cer), y = log2(var_cer))) + 
#   geom_point(aes(color = coexpressed), alpha = 0.5) +
#   theme_classic() +
#   theme(legend.title = element_blank(),
#         legend.position = "left") +
#   xlab("mean expression (log2)") +
#   ylab("variance expression (log2)") +
#   ggtitle("Variance vs mean in S. cerevisiae") +
#   scale_color_discrete(type = c("gold",
#                                 "mediumseagreen",
#                                 "orange1",
#                                 "blue2",
#                                 "grey"),
#                        limits = c("conserved co-expressed",
#                                   "diverged co-expressed",
#                                   "S. cerevisiae co-expressed",
#                                   "S. paradoxus co-expressed",
#                                   "never co-expressed"))
# ppar_meanvar <- ggplot(plotdf, aes(x = log2(mean_par), y = log2(var_par))) + 
#   geom_point(aes(color = coexpressed), alpha = 0.5) +
#   theme_classic() +
#   theme(legend.title = element_blank(),
#         legend.position = "left") +
#   xlab("mean expression (log2)") +
#   ylab("variance expression (log2)") +
#   ggtitle("Variance vs mean in S. paradoxus") +
#   scale_color_discrete(type = c("gold",
#                                 "mediumseagreen",
#                                 "orange1",
#                                 "blue2",
#                                 "grey"),
#                        limits = c("conserved co-expressed",
#                                   "diverged co-expressed",
#                                   "S. cerevisiae co-expressed",
#                                   "S. paradoxus co-expressed",
#                                   "never co-expressed"))
# library(ggExtra)
# pcer_meanvar <- ggMarginal(pcer_meanvar, groupColour = TRUE, groupFill = TRUE, margins = "y")
# ppar_meanvar <- ggMarginal(ppar_meanvar, groupColour = TRUE, groupFill = TRUE, margins = "y")
# 
# pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/mean_var_cer10.pdf",
#     width = 6, height = 4)
# pcer_meanvar
# dev.off()
# pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/mean_var_par10.pdf",
#     width = 6, height = 4)
# ppar_meanvar
# dev.off()
# # no obvious differences in variance, if anything grey cerevisiae genes are more variable,
# # but there's only 120 of them, so a couple might be disporpotionately affecting the distribution
# 
# ### formal test of differences in variance between coexpression groups
# # with all module distinctions
# testVarDifBetweenGroups <- function(.df, .mod_gdf) {
#   df0 <- select(ungroup(.df), gene_name, var_cer, mean_cer, var_par, mean_par) |> 
#     left_join(select(.mod_gdf, gene_name, coexpressed),
#               by = c("gene_name")) |> drop_na() |> 
#     mutate(log2var_cer = log2(var_cer),
#            log2var_par = log2(var_par),
#            log2mean_cer = log2(mean_cer),
#            log2mean_par = log2(mean_par))
#   # cer
#   modcer <- lm(log2var_cer ~ log2mean_cer, data = df0)
#   df0$residuals <- modcer$residuals
#   group_means <- df0 |> group_by(coexpressed) |> 
#     summarise(mean_res = mean(residuals)) |> 
#     select(mean_res, coexpressed)
#   pvals <- compare_means(residuals ~ coexpressed, data = df0) |> 
#     slice(c(1,2,3,4,8,9,10))
#   pvals$y.position <- c(8, 9.5, 11, 12.5, 5, 6.5, 8)
#   pcer <- ggplot(df0, aes(x = coexpressed, y = residuals)) +
#     geom_jitter(aes(color = coexpressed)) + 
#     scale_color_discrete(type = c("gold",
#                                   "mediumseagreen",
#                                   "orange1",
#                                   "blue2",
#                                   "grey"),
#                          limits = c("conserved co-expressed",
#                                     "diverged co-expressed",
#                                     "S. cerevisiae co-expressed",
#                                     "S. paradoxus co-expressed",
#                                     "never co-expressed")) +
#     ylab("residual variance") +
#     guides(x = guide_axis(angle = 90)) +
#     theme_classic() +
#     theme(legend.position = "none") +
#     stat_pvalue_manual(pvals, label = "wilcox.test p = {p.format}", tip.length = 0.01) +
#     geom_segment(data = group_means,
#                  aes(x = seq(from = 0.8, to = 4.8, by = 1),
#                      xend = seq(from = 1.2, to = 5.2, by = 1),
#                      y = as.numeric(mean_res), yend = as.numeric(mean_res)), color = "red") +
#     ggtitle("S. cerevisiae residual variance")
#   # par
#   modpar <- lm(log2var_par ~ log2mean_par, data = df0)
#   df0$residuals <- modpar$residuals
#   group_means <- df0 |> group_by(coexpressed) |> 
#     summarise(mean_res = mean(residuals)) |> 
#     select(mean_res, coexpressed)
#   pvals <- compare_means(residuals ~ coexpressed, data = df0) |> 
#     slice(c(1,2,3,4,8,9,10))
#   pvals$y.position <- c(8, 9.5, 11, 12.5, 5, 6.5, 8)
#   ppar <- ggplot(df0, aes(x = coexpressed, y = residuals)) +
#     geom_jitter(aes(color = coexpressed)) + 
#     scale_color_discrete(type = c("gold",
#                                   "mediumseagreen",
#                                   "orange1",
#                                   "blue2",
#                                   "grey"),
#                          limits = c("conserved co-expressed",
#                                     "diverged co-expressed",
#                                     "S. cerevisiae co-expressed",
#                                     "S. paradoxus co-expressed",
#                                     "never co-expressed")) +
#     ylab("residual variance") +
#     guides(x = guide_axis(angle = 90)) +
#     theme_classic() +
#     theme(legend.position = "none") +
#     stat_pvalue_manual(pvals, label = "wilcox.test p = {p.format}", tip.length = 0.01) +
#     geom_segment(data = group_means,
#                  aes(x = seq(from = 0.8, to = 4.8, by = 1),
#                      xend = seq(from = 1.2, to = 5.2, by = 1),
#                      y = as.numeric(mean_res), yend = as.numeric(mean_res)), color = "red") +
#     ggtitle("S. paradoxus residual variance")
#   
#   return(ggarrange(pcer, ppar, nrow = 1, ncol = 2))
# }
# # merge 00
# testVarDifBetweenGroups(plotdf, module_genedf00)
# # merge 10
# pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/residual_var10.pdf",
#     width = 12, height = 6)
# testVarDifBetweenGroups(plotdf, module_genedf10)
# dev.off()
# # merge 25
# testVarDifBetweenGroups(plotdf, module_genedf25)
# # merge 35
# testVarDifBetweenGroups(plotdf, module_genedf35)
# # The pvalues are quite close whichever module distinction you use
# # (so close in fact, they can look identical in the graph,
# # but I did verify they are different)
# 
# #### Do unsigned network modules contain signed network modules? ####
# 
# # just for exploration's sake
# # we would expect unsigned network modules to contain signed network module pairs,
# # 1 positively correlated module and one negatively correlated module (if the pair exists, it seems rare)
# # is this the case? how similar are signed and unsigned networks?
# 
# # using unmerged networks, as they have the most modules and the best chance of catching this
# signed_cer <- colors00$cer
# signed_par <- colors00$par
# load("data_files/Networks.RData")
# unsigned_cer <- colors00$cer
# unsigned_par <- colors00$par
# 
# # bipartites
# # cer
# cer_signedunsigned_counts <- makeBipartiteMatrix(.top_colors = signed_cer, .bottom_colors = unsigned_cer)
# plotweb(cer_signedunsigned_counts$matrix, col.high = cer_signedunsigned_counts$colors_high, 
#         col.low = cer_signedunsigned_counts$colors_low, labsize = 1,
#         high.lablength = 0, low.lablength = 0)
# # par
# par_signedunsigned_counts <- makeBipartiteMatrix(.top_colors = signed_par, .bottom_colors = unsigned_par)
# plotweb(par_signedunsigned_counts$matrix, col.high = par_signedunsigned_counts$colors_high, 
#         col.low = par_signedunsigned_counts$colors_low, labsize = 1,
#         high.lablength = 0, low.lablength = 0)
# 
# # both species signed for comparison (we'll do this for real in the bipartite section below)
# test <- makeBipartiteMatrix(.top_colors = signed_cer, .bottom_colors = signed_par)
# plotweb(test$matrix, col.high = test$colors_high, 
#         col.low = test$colors_low, labsize = 1,
#         high.lablength = 0, low.lablength = 0)
# 
# # mostly they're just quite similar (note how many fewer connections aka non-zero 
# # entries in matrix there are total versus between-species comparisons)
# sum(cer_signedunsigned_counts$matrix == 0)
# sum(par_signedunsigned_counts$matrix == 0)
# sum(test$matrix == 0)
# 
# # But do unsigned networks contain signed networks? No. Possibly because the recursive
# # module selection won't allow modules large enough to contain pairs of other modules,
# # but it doesn't feel worth it to investigate right now
# 
# # remember to reload the desired networks at the end
# rm(colors00, colors10, colors25, colors35)
# load("data_files/NetworksSigned.RData")
# 

#### troubleshooting bipartite mis-coloring ####
# Archived b/c this problem was solved (plotweb defaults to order groups for the fewest crossovers,
# you can override this behavior by changing method = "cca" to method = "normal")

# # For some very cryptic reason, creating a bipartite with module names
# # instead of colors then providing plot web with the corresponding colors rearranges everything!
# 
# # Example, First changing colors in parents
# parents_color_counts <- makeBipartiteMatrix(.top_colors = colors25$cer, 
#                                             .bottom_colors = colors25$par)
# plotweb(parents_color_counts$matrix, col.high = parents_color_counts$colors_high,
#         col.low = parents_color_counts$colors_low)
# new_colors <- if_else(parents_color_counts$colors_high %in% c("black", "red", "brown", "grey"),
#                       true = "hotpink", false = "mediumseagreen")
# new_colors
# plotweb(parents_color_counts$matrix, col.high = new_colors,
#         col.low = parents_color_counts$colors_low) # this works
# # now parents but using module names for one parent
# parents_color_counts <- makeBipartiteMatrix(.top_colors = module_genedf25$module_name, 
#                                             .bottom_colors = colors25$par)
# parents_color_counts$colors_high
# many_colors <- standardColors(length(parents_color_counts$colors_high))
# plotweb(parents_color_counts$matrix, col.high = many_colors,
#         col.low = parents_color_counts$colors_low)
# which(grepl("nev", parents_color_counts$colors_high))
# many_colors[30] # grey-grey group is steelblue
# # assigning 5 colors to the 60 modules
# new_colors <- sapply(parents_color_counts$colors_high, \(m) {
#   if (grepl("ccm", m)) {
#     return("gold")
#   }
#   if (grepl("div", m)) {
#     return("mediumseagreen")
#   }
#   if (grepl("cer", m)) {
#     return("orange1")
#   }
#   if (grepl("par", m)) {
#     return("blue2")
#   }
#   if (grepl("nev", m)) {
#     return("grey")
#   }
# }) |> unlist() |> as.character()
# parents_color_counts$colors_high[30]
# new_colors[30] # grey-grey is now grey
# plotweb(parents_color_counts$matrix, col.high = new_colors,
#         col.low = parents_color_counts$colors_low) # grey-grey is correctly grey
# 
# # repeat with hybrid
# hybrid_color_counts <- makeBipartiteMatrix(.top_colors = module_genedf25$module_name, 
#                                            .bottom_colors = colors25$hyb)
# hybrid_color_counts$colors_high
# many_colors <- standardColors(length(hybrid_color_counts$colors_high))
# plotweb(hybrid_color_counts$matrix, col.high = many_colors,
#         col.low = hybrid_color_counts$colors_low)
# which(grepl("nev", parents_color_counts$colors_high))
# many_colors[30] # grey-grey group is steelblue
# # assigning 5 colors to the 60 modules
# new_colors <- sapply(hybrid_color_counts$colors_high, \(m) {
#   if (grepl("ccm", m)) {
#     return("gold")
#   }
#   if (grepl("div", m)) {
#     return("mediumseagreen")
#   }
#   if (grepl("cer", m)) {
#     return("orange1")
#   }
#   if (grepl("par", m)) {
#     return("blue2")
#   }
#   if (grepl("nev", m)) {
#     return("grey")
#   }
# }) |> unlist() |> as.character()
# hybrid_color_counts$colors_high[30]
# new_colors[30] # grey-grey is now grey
# plotweb(hybrid_color_counts$matrix, col.high = new_colors,
#         col.low = hybrid_color_counts$colors_low) # grey-grey is correctly grey
# simplemat <- hybrid_color_counts$matrix
# simplemat[simplemat < 10] <- 0
# plotweb(simplemat, col.high = new_colors,
#         col.low = hybrid_color_counts$colors_low) # ah the problem is simplifying
# simplemat[,30]
# # simplifying earlier
# hybrid_color_counts <- makeBipartiteMatrix(.top_colors = module_genedf25$module_name, 
#                                            .bottom_colors = colors25$hyb, .simplify_thresh = 10)
# hybrid_color_counts$colors_high
# many_colors <- standardColors(length(hybrid_color_counts$colors_high))
# plotweb(hybrid_color_counts$matrix, col.high = many_colors,
#         col.low = hybrid_color_counts$colors_low) # grey grey is brown now...
# which(grepl("nev", hybrid_color_counts$colors_high))
# many_colors[30] # grey-grey group should be steelblue
# 
# # Why does eliminating connections change plotweb colors??
# test <- makeBipartiteMatrix(c("red", "red", "green", "red", "grey"), 
#                             c("green", "green", "green", "red", "grey"))
# test$matrix
# plotweb(test$matrix, col.high = test$colors_high, col.low = test$colors_low)
# test$matrix[2, 1] <- 0
# test$matrix
# # eliminating the red-red connection shifts grey to interact with red...
# plotweb(test$matrix, col.high = test$colors_high, col.low = test$colors_low)
# test$matrix # even though in the matrix it still shows it only interacts with grey
# test0 <- makeBipartiteMatrix(c("red", "red", "green", "grey"), 
#                              c("green", "green", "green", "grey"))
# identical(test0$matrix, test$matrix)
# # ah. is the problem groups that have no other connections?
# test1 <- test$matrix[c(1,3),]
# plotweb(test$matrix, col.high = test$colors_high, col.low = test$colors_low, empty = FALSE)
# plotweb(test1, col.high = colnames(test1), col.low = rownames(test1))
# 
# # adding empty = FALSE argument
# hybrid_color_counts <- makeBipartiteMatrix(.top_colors = module_genedf25$module_name, 
#                                            .bottom_colors = colors25$hyb, .simplify_thresh = 10)
# hybrid_color_counts$colors_high
# many_colors <- standardColors(length(hybrid_color_counts$colors_high))
# plotweb(hybrid_color_counts$matrix, col.high = many_colors,
#         col.low = hybrid_color_counts$colors_low) # grey grey is brown now...
# plotweb(hybrid_color_counts$matrix, col.high = many_colors,
#         col.low = hybrid_color_counts$colors_low, empty = FALSE) # colors are correct now