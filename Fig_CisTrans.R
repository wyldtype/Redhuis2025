setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Redhuis2025/")
source("functions_for_figure_scripts.R")
load("data_files/FinalDataframe3Disp.RData")
load("data_files/Cleaned_Counts.RData")
load("data_files/Cleaned_Counts_Allele.RData")

cor_thresh <- 0.8

#### Log2 Fold Change Controls/Sanity Checks ####

### Do level divergers have more correlated l2fc between parents and hybrid alleles?

# R squared for conserved-level genes
mod <- lm(effect_size_allele ~ effect_size_species, data = filter(finaldf, level == "conserved"))
summary(mod)
# R squared for level-diverged genes
mod <- lm(effect_size_allele ~ effect_size_species, data = filter(finaldf, level == "diverged"))
summary(mod) # lower
# Are log2 fold changes less correlated among level diverged genes b/c as they get
# larger they get noisier?

# using LowPi-unique upcer level genes---transest case u no
length(lowpi_upcer_idxs)
idx <- 1
while (idx < length(lowpi_upcer_idxs)) {
  cat(idx, "\n")
  print(annotate_figure(plotEnvironments(lowpi_upcer_idxs[idx], .quartet = TRUE),
                        top = lowpi_upcer_idxs[idx]))
  idx <- idx + 1
}
# conclusion: it's not all 83 genes, but lots do have
# Spar drop off early on and both hybrid alleles climb
# with Scer
# Question: are the ones that have this behavior more
# likely to be diverging in dynamics?
lowpi_levdyn_idxs <- finaldf |> filter(experiment == "LowPi" &
                                         dynamics == "diverged") |> 
  select(gene_name) |> pull() |> 
  intersect(y = lowpi_upcer_idxs)
idx <- 1
while (idx < length(lowpi_levdyn_idxs)) {
  cat(idx, "\n")
  print(annotate_figure(plotEnvironments(lowpi_levdyn_idxs[idx], .quartet = TRUE),
                        top = lowpi_levdyn_idxs[idx]))
  idx <- idx + 1
}
# no, about half of them still

# observation from average expression: level divergers are way
# more likely to have hybrid l2fc in same direction as parents
# can we confirm this?
ggplot(finaldf, aes(x = effect_size_species,
                    y = effect_size_allele)) +
  geom_point(aes(color = group4)) +
  scale_color_discrete(breaks = levdyn_colordf$limits,
                       type = levdyn_colordf$type,
                       limits = levdyn_colordf$limits) +
  facet_wrap(~group4) +
  xlim(c(-5, 5)) +
  ylim(c(-5, 5))
ggMarginal(ggplot(finaldf, aes(x = effect_size_species,
                               y = effect_size_allele)) +
             geom_point(aes(color = group4)) +
             scale_color_discrete(breaks = levdyn_colordf$limits,
                                  type = levdyn_colordf$type,
                                  limits = levdyn_colordf$limits) +
             xlim(c(-5, 5)) +
             ylim(c(-5, 5)),
           groupColour = TRUE, groupFill = TRUE)
# you can kind of see that diverged level genes 
# have stronger effect_size_allele, but can't tell
# that it's in the same direction as parent

# TODO: parent-hybrid lfc correlations for each divergence group


#### Hybrid expression in environment specific vs environment robust level-divergers ####

# TODO: come back to this after finding a better
# way of visualizing hybrid allele l2fc
# (because level divergers should gernally have correlated
# l2fc but it's so noisy it's hard to see this)

# idxs generated in Fig_environmental_patterns for now
length(constitutive_upcer_tagseq_idxs)
length(constitutive_uppar_tagseq_idxs)
length(constitutive_upcer_heatcold_idxs)
length(constitutive_uppar_heatcold_idxs)
length(hap4_upcer_idxs)
length(cc_upcer_idxs)
length(lown_upcer_idxs)
length(lowpi_upcer_idxs)
length(heat_upcer_idxs)
length(cold_upcer_idxs)
length(hap4_uppar_idxs)
length(cc_uppar_idxs)
length(lown_uppar_idxs)
length(lowpi_uppar_idxs)
length(heat_uppar_idxs)
length(cold_uppar_idxs)
env_robust_idxs <- unique(c(constitutive_upcer_tagseq_idxs,
                            constitutive_uppar_tagseq_idxs,
                            constitutive_upcer_heatcold_idxs,
                            constitutive_uppar_heatcold_idxs))
env_specific_idxs <- unique(c(hap4_upcer_idxs, cc_upcer_idxs,
                              lown_upcer_idxs, lowpi_upcer_idxs,
                              heat_upcer_idxs, cold_upcer_idxs,
                              hap4_uppar_idxs, cc_uppar_idxs,
                              lown_uppar_idxs, lowpi_uppar_idxs,
                              heat_uppar_idxs, cold_uppar_idxs))
env_specific_list <- list(HAP4 = c(hap4_upcer_idxs, hap4_uppar_idxs),
                          CC = c(cc_upcer_idxs, cc_uppar_idxs),
                          LowN = c(lown_upcer_idxs, lown_uppar_idxs),
                          LowPi = c(lowpi_upcer_idxs, lowpi_uppar_idxs),
                          Heat = c(heat_upcer_idxs, heat_uppar_idxs),
                          Cold = c(cold_upcer_idxs, cold_uppar_idxs))

### Log2 fold change parents vs hybrids env specific vs constitutive level divergers
plotdf <- filter(finaldf, level == "diverged") |> 
  mutate(env_specific_or_robust = if_else(gene_name %in% env_robust_idxs,
                                          true = "robust",
                                          false = if_else(gene_name %in% env_specific_idxs,
                                                          true = "specific",
                                                          false = "neither")))
plotdf$env_specific_or_robust |> table()
ggplot(plotdf, aes(x = abs(effect_size_species))) +
  geom_density(aes(fill = env_specific_or_robust), alpha = 0.5) +
  facet_wrap(~experiment) +
  theme_classic() # robust have larger effect sizes
for (e in ExperimentNames) {
  e_idxs <- env_specific_list[[e]]
  plotdf <- filter(finaldf, level == "diverged" & experiment == e) |> 
    mutate(env_specific_or_robust = if_else(gene_name %in% env_robust_idxs,
                                            true = "robust",
                                            false = if_else(gene_name %in% e_idxs,
                                                            true = "specific",
                                                            false = "neither")))
  p_scatter <- ggplot(plotdf,
                      aes(x = effect_size_species,
                          y = effect_size_allele)) +
    geom_abline(slope = 1, intercept = 0, color = "black") +
    geom_point(aes(color = env_specific_or_robust)) +
    facet_wrap(~env_specific_or_robust) +
    theme_classic() +
    labs(color = "gene group") +
    xlab("log2 fold change, parents") +
    ylab("log2 fold change, hybrid") +
    ggtitle(e)
  p_dens <- ggplot(plotdf, aes(x = log2(abs(effect_size_allele))*sign(effect_size_allele))) +
    geom_density(aes(fill = env_specific_or_robust), alpha = 0.5) +
    theme_classic() +
    ggtitle(e)
  print(ggarrange(p_scatter, p_dens, nrow = 2, ncol = 1))
}

#### Divergence in level is preserved in the hybrid ####
plotdf <- filter(finaldf, !(level %in% c("low_expr", "biased")))
# scatter plot of LFCs in parent and hybrid with dynamics divergence colored
# R squared for all (non-lowly expressed) genes
mod <- lm(effect_size_allele ~ effect_size_species, data = plotdf)
summary(mod)
# R squared for level-diverged genes
mod <- lm(effect_size_allele ~ effect_size_species, data = filter(plotdf, level == "diverged"))
summary(mod)
# plotting
p_scatter <- ggplot(plotdf,
                    aes(x = effect_size_species,
                        y = effect_size_allele)) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_vline(xintercept = c(-1, 1), color = "black", linetype = "dashed") +
  geom_point(aes(color = group4)) +
  facet_wrap(~experiment) +
  scale_color_discrete(type = levdyn_colordf$type,
                       labels = levdyn_colordf$labels,
                       limits = levdyn_colordf$limits) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(color = "gene group") +
  xlab("log2 fold change, parents") +
  ylab("log2 fold change, hybrid") +
  ggtitle("") +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  stat_cor(cor.coef.name = "R")

# accompanying barplot to show that most of these genes have conserved dynamics
p_bar <- ggplot(plotdf, aes(x = group4)) +
  geom_bar(aes(fill = group4)) +
  geom_text(stat='count', aes(label = after_stat(count)), vjust=-1) +
  scale_fill_discrete(type = levdyn_colordf$type,
                      labels = levdyn_colordf$labels,
                      limits = levdyn_colordf$limits) +
  scale_x_discrete(limits = levdyn_colordf$limits) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        legend.position = "right",
        plot.margin = margin(t = 1, r = 0.1, b = 0.1, l = 0.1, "cm")) +
  ggtitle("genes in each group") +
  facet_wrap(~experiment) +
  ylim(c(0, 4100))

# level and dynamics summary
library(ggExtra)
ggarrange(p_scatter, p_bar, nrow = 1, ncol = 2,
          common.legend = TRUE, widths = c(2, 1), hjust = 1,
          legend = "bottom")

pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/CisTrans/l2fc.pdf"),
    width = 4, height = 4)
p_scatter
dev.off()

pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/CisTrans/bar.pdf"),
    width = 6, height = 4)
p_bar
dev.off()

#### Calculating ortholog correlations ####
getCorelation <- function(.gene_name, .experiment, 
                          .cts1, .cts2, .info1, .info2) {
  common_cols <- intersect(colnames(.cts1[,.info1$experiment == .experiment]), 
                           colnames(.cts2[,.info2$experiment == .experiment]))
  output <- cor(as.numeric(.cts1[.gene_name, common_cols]),
                as.numeric(.cts2[.gene_name, common_cols]))
  return(as.numeric(output))
}
# tests for getCorelation
getCorelation("YGR192C", "LowN",
              .cts1 = collapsed$cer,
              .cts2 = collapsed$par,
              .info1 = info,
              .info2 = info)


# adding corelation columns
finaldf$cor_hybrid <- map2(finaldf$gene_name, finaldf$experiment, getCorelation,
                           .cts1 = collapsed_allele$cer,
                           .cts2 = collapsed_allele$par,
                           .info1 = info_allele,
                           .info2 = info_allele) |> unlist()
finaldf$cor_parents <- map2(finaldf$gene_name, finaldf$experiment, getCorelation,
                            .cts1 = collapsed$cer,
                            .cts2 = collapsed$par,
                            .info1 = info,
                            .info2 = info) |> unlist()
finaldf$cor_scer <- map2(finaldf$gene_name, finaldf$experiment, getCorelation,
                         .cts1 = collapsed_allele$cer,
                         .cts2 = collapsed$cer,
                         .info1 = info_allele,
                         .info2 = info) |> unlist()
finaldf$cor_spar <- map2(finaldf$gene_name, finaldf$experiment, getCorelation,
                         .cts1 = collapsed_allele$par,
                         .cts2 = collapsed$par,
                         .info1 = info_allele,
                         .info2 = info) |> unlist()

#### Trans-varying plasticity shows parental dominance in some environments ####

# TODO: as the two environments with strongest dominance were 
# the longest running experiments,
# limit LowN and HAP4 to 1hr see if the dominance
# still holds up

for (e in unique(finaldf$experiment)) {
  cat(e, "\n")
  p_density <- ggplot(filter(finaldf, experiment == e),
                      aes(x = cor_hybrid)) +
    geom_density(aes(fill = experiment)) +
    geom_vline(xintercept = cor_thresh) +
    theme_classic()
  plotdf <- finaldf |> filter(experiment == e & 
                                cor_hybrid > cor_thresh &
                                dynamics != "conserved")
  p_scat <- ggplot(plotdf, aes(x = cor_scer, y = cor_spar)) +
    geom_point(aes(color = paste(cer, par),
                   shape = paste(cer, par))) +
    theme_classic() +
    geom_vline(xintercept = cor_thresh) +
    geom_hline(yintercept = cor_thresh)
  if (nrow(plotdf) > 0) {
    # collecting top dynamics category for each parental dominance and no dominance
    # Scer dominance
    scer_dom <- plotdf |> filter(cor_scer > cor_thresh &
                                   cor_spar <= cor_thresh) |> 
      group_by(cer, par) |> 
      summarise(n = n()) |> 
      ungroup() |> 
      arrange(desc(n)) |> 
      slice(1)
    gene_idxs <- plotdf |> filter(cer == scer_dom$cer,
                                  par == scer_dom$par &
                                    cor_scer > cor_thresh &
                                    cor_spar <= cor_thresh) |> 
      select(gene_name) |> pull()
    p_scer <- annotate_figure(plotGenes(gene_idxs, .normalization = "scale",
                                        .experiment_name = e, .quartet = TRUE),
                              top = paste0(length(gene_idxs), " genes\nScer dominance"))
    # Spar dominance
    spar_dom <- plotdf |> filter(cor_spar > cor_thresh &
                                   cor_scer <= cor_thresh) |> 
      group_by(cer, par) |> 
      summarise(n = n()) |> 
      ungroup() |> 
      arrange(desc(n)) |> 
      slice(1)
    gene_idxs <- plotdf |> filter(cer == spar_dom$cer,
                                  par == spar_dom$par &
                                    cor_spar > cor_thresh &
                                    cor_scer <= cor_thresh) |> 
      select(gene_name) |> pull()
    p_spar <- annotate_figure(plotGenes(gene_idxs, .normalization = "scale",
                                        .experiment_name = e, .quartet = TRUE),
                              top = paste0(length(gene_idxs), " genes\nSpar dominance"))
    # No parental dominance
    no_dom <- plotdf |> filter(cor_spar <= cor_thresh &
                                 cor_scer <= cor_thresh) |> 
      group_by(cer, par) |> 
      summarise(n = n()) |> 
      ungroup() |> 
      arrange(desc(n)) |> 
      slice(1)
    gene_idxs <- plotdf |> filter(cer == no_dom$cer,
                                  par == no_dom$par &
                                    cor_spar <= cor_thresh &
                                    cor_scer <= cor_thresh) |> 
      select(gene_name) |> pull()
    p_no <- annotate_figure(plotGenes(gene_idxs, .normalization = "scale",
                                      .experiment_name = e, .quartet = TRUE),
                            top = paste0(length(gene_idxs), " genes\nNo parental dominance"))
    p_blank <- ggplot() + theme_void()
    top_row <- ggarrange(p_density, p_scat, nrow = 1, ncol = 2)
    bottom_row <- ggarrange(p_scer, p_spar, p_no, nrow = 1, ncol = 3)
    # assembling plots
    pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/CisTrans/",
                      e, "_dominance.pdf"),
        width = 7, height = 5)
    print(ggarrange(top_row, bottom_row, nrow = 2, ncol = 1))
    dev.off()
  }
  else {
    pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/CisTrans/",
                      e, "_dominance.pdf"),
        width = 3, height = 3)
    print(p_density)
    dev.off()
  }
}

##### Investigation: calculating cor_hybrid as pairs of replicates ####
# test if calculating cor_hybrid using counts that don't collapse replicates
# changes cor_hybrid distribution, especially in LowPi
info_allele_cer <- filter(sample_info_allele, grepl("hyc", sample_name))
info_allele_par <- filter(sample_info_allele, grepl("hyp", sample_name))
counts_allele_cer <- counts_allele[, info_allele_cer$sample_name]
counts_allele_par <- counts_allele[, info_allele_par$sample_name]
colnames(counts_allele_cer) <- gsub(pattern = "hyc", replacement = "hyb",
                                    colnames(counts_allele_cer))
colnames(counts_allele_par) <- gsub(pattern = "hyp", replacement = "hyb",
                                    colnames(counts_allele_par))
finaldf$cor_hybrid_withReps <- map2(finaldf$gene_name, finaldf$experiment, getCorelation,
                           .cts1 = counts_allele_cer,
                           .cts2 = counts_allele_par,
                           .info1 = info_allele_cer,
                           .info2 = info_allele_par) |> unlist()

finaldf |> filter(experiment == "LowPi") |> 
  mutate(mostly_trans = cor_hybrid_withReps > cor_thresh) |> 
  select(mostly_trans) |> table() # with reps
finaldf |> filter(experiment == "LowPi") |> 
  mutate(mostly_trans = cor_hybrid > cor_thresh) |> 
  select(mostly_trans) |> table() # no reps
test_gene <- finaldf |> filter(experiment == "LowPi" &
                    cor_hybrid_noReps > cor_thresh &
                    cor_hybrid < cor_thresh) |> 
  slice_sample(n = 1) |> select(gene_name) |> pull()
# collapsed counts
collapsed_allele$cer[test_gene, info_allele$experiment == "LowPi"]
collapsed_allele$par[test_gene, info_allele$experiment == "LowPi"]
cor(collapsed_allele$cer[test_gene, info_allele$experiment == "LowPi"],
    collapsed_allele$par[test_gene, info_allele$experiment == "LowPi"])
plot(x = collapsed_allele$cer[test_gene, info_allele$experiment == "LowPi"],
     y = collapsed_allele$par[test_gene, info_allele$experiment == "LowPi"])
# replicated counts
test_hyc <- counts_allele[test_gene, sample_info_allele$experiment == "LowPi" &
                sample_info_allele$allele == "cer"]
test_hyp <- counts_allele[test_gene, sample_info_allele$experiment == "LowPi" &
                sample_info_allele$allele == "par"]
cor(as.numeric(test_hyc), as.numeric(test_hyp))
plotdf <- tibble(hyc = as.numeric(test_hyc),
                 hyp = as.numeric(test_hyp),
                 sample_name = names(test_hyc)) |> 
  left_join(sample_info_allele, by = "sample_name")
ggplot(plotdf, aes(x = hyc, y = hyp)) + 
  geom_point(aes(color = time_point_num))
ggplot(filter(finaldf, experiment == "LowPi"), aes(x = cor_hybrid)) +
  geom_density()
ggplot(filter(finaldf, experiment == "LowPi"), aes(x = cor_hybrid_noReps)) +
  geom_density()
mean(test_hyc[grepl(pattern = "TP1", names(test_hyc))])
mean(test_hyc[grepl(pattern = "TP2", names(test_hyc))])
mean(test_hyc[grepl(pattern = "TP3", names(test_hyc))])
mean(test_hyp[grepl(pattern = "TP1", names(test_hyc))])
mean(test_hyp[grepl(pattern = "TP2", names(test_hyc))])
mean(test_hyp[grepl(pattern = "TP3", names(test_hyc))])

#### Investigation: Slowing down Hybrid by spreading out timepoints in Heat/Cold ####
# Based on PCA, hybrid seems to respond to heat stress faster
# what if we remove 60min timepoint from hybrid and make 30min the new final timepoint
# so cor(x = c(parent0, parent15, parent60), y = c(hybrid0, hybrid15, hybrid30))
timedf <- finaldf |> filter(experiment == "Heat")
collapsed_allele_cer_time <- collapsed_allele$cer[, info_allele$time_point_num != 60 &
                                                    info_allele$experiment == "Heat"]
collapsed_allele_par_time <- collapsed_allele$par[, info_allele$time_point_num != 60 &
                                                    info_allele$experiment == "Heat"]
colnames(collapsed_allele_cer_time) <- gsub(pattern = "30", replacement = "60", colnames(collapsed_allele_cer_time))
colnames(collapsed_allele_par_time) <- gsub(pattern = "30", replacement = "60", colnames(collapsed_allele_par_time))
info_allele_time <- filter(info_allele, condition %in% colnames(collapsed_allele_cer_time))

collapsed_cer_time <- collapsed$cer[, info$time_point_num != 30 &
                                      info$experiment == "Heat"]
collapsed_par_time <- collapsed$par[, info$time_point_num != 30 &
                                      info$experiment == "Heat"]
info_time <- filter(info, condition %in% colnames(collapsed_cer_time))

timedf$cor_hybrid <- map2(timedf$gene_name, timedf$experiment, getCorelation,
                           .cts1 = collapsed_allele_cer_time,
                           .cts2 = collapsed_allele_par_time,
                           .info1 = info_allele_time,
                           .info2 = info_allele_time) |> unlist()
timedf$cor_parents <- map2(timedf$gene_name, timedf$experiment, getCorelation,
                            .cts1 = collapsed_cer_time,
                            .cts2 = collapsed_par_time,
                            .info1 = info_time,
                            .info2 = info_time) |> unlist()
timedf$cor_scer <- map2(timedf$gene_name, timedf$experiment, getCorelation,
                         .cts1 = collapsed_allele_cer_time,
                         .cts2 = collapsed_cer_time,
                         .info1 = info_allele_time,
                         .info2 = info_time) |> unlist()
timedf$cor_spar <- map2(timedf$gene_name, timedf$experiment, getCorelation,
                         .cts1 = collapsed_allele_par_time,
                         .cts2 = collapsed_par_time,
                         .info1 = info_allele_time,
                         .info2 = info_time) |> unlist()

plotdf <- timedf |> filter(cor_hybrid > cor_thresh &
                             dynamics != "conserved")
ggplot(plotdf, aes(x = cor_scer, y = cor_spar)) +
  geom_point(aes(color = paste(cer, par),
                 shape = paste(cer, par)))
# now paradoxus parent has strong dominance
# three most common expression shapes amongst high spar corr:
timedf |> filter(cor_hybrid > cor_thresh &
                   dynamics != "conserved" &
                   cor_spar > cor_thresh) |> 
  select(cer, par) |> table()
# 2-1
gene_idxs <- timedf |> filter(cor_hybrid > cor_thresh &
                                dynamics != "conserved" &
                                cor_spar > cor_thresh &
                                cer == 2 & par == 1) |> 
  select(gene_name) |> pull()
plotGenes(gene_idxs, .experiment_name = "Heat",
          .quartet = TRUE)
# 0-1
gene_idxs <- timedf |> filter(cor_hybrid > cor_thresh &
                                dynamics != "conserved" &
                                cor_spar > cor_thresh &
                                cer == 0 & par == 1) |> 
  select(gene_name) |> pull()
plotGenes(gene_idxs, .experiment_name = "Heat",
          .quartet = TRUE)
# 1-0
gene_idxs <- timedf |> filter(cor_hybrid > cor_thresh &
                                dynamics != "conserved" &
                                cor_spar > cor_thresh &
                                cer == 1 & par == 0) |> 
  select(gene_name) |> pull()
plotGenes(gene_idxs, .experiment_name = "Heat",
          .quartet = TRUE)

# does this also work for cold?
# basically no. Kind of but not as distinctly and when we
# plotted avg expression of the parent dominant genes,
# the hybrid just had low variation

# Looking at random genes to see if they really do have less parental dominance
# in Heat/Cold
# Heat/Cold
test_exp <- "Heat"
gene_idxs <- faydf |> filter(cor_scer < cor_thresh &
                               cor_spar < cor_thresh &
                               cor_hybrid > cor_thresh &
                               cer != 0 & par != 0 &
                               experiment == test_exp) |> 
  select(gene_name) |> pull()
plotGenes(gene_idxs, .experiment_name = test_exp,
          .quartet = TRUE, .plotlims = c(4, 10))
test_idx <- sample(gene_idxs, size = 1)
faydf |> filter(gene_name == test_idx & experiment == test_exp) |> 
  select(gene_name, cor_hybrid, cor_parents, cor_scer, cor_spar)
annotate_figure(plotGenes(test_idx, .experiment_name = test_exp,
          .quartet = TRUE, .plotlims = c(4, 10)), top = test_idx)
# LowN
gene_idxs <- finaldf |> filter(cor_scer < cor_thresh &
                               cor_spar < cor_thresh &
                               cor_hybrid > cor_thresh &
                               cer != 0 & par != 0 &
                               experiment == "LowN") |> 
  select(gene_name) |> pull()
test_idx <- sample(gene_idxs, size = 1)
finaldf |> filter(gene_name == test_idx & experiment == "LowN") |> 
  select(gene_name, cor_hybrid, cor_parents, cor_scer, cor_spar)
annotate_figure(plotGenes(test_idx, .experiment_name = "LowN",
                          .quartet = TRUE, .plotlims = c(4, 9)), top = test_idx)

################# Scraps
# "pure cis":
gene_idx <- filter(plotdf, hybrid_cor < 0.8 &
                     dynamics != "conserved" &
                     cer_cor > 0.8 & par_cor > 0.8 &
                     cer == 2 & par == 1) |> select(gene_name) |> pull()
plotGenes(gene_idx[1], .experiment_name = "HAP4",
          .quartet = TRUE, .normalization = "scale")
plotGenes(gene_idx[2], .experiment_name = "HAP4",
          .quartet = TRUE, .normalization = "scale")
plotGenes(gene_idx[3], .experiment_name = "HAP4",
          .quartet = TRUE, .normalization = "scale")
plotGenes(gene_idx[4], .experiment_name = "HAP4",
          .quartet = TRUE, .normalization = "log2")
plotGenes(gene_idx[5], .experiment_name = "HAP4",
          .quartet = TRUE, .normalization = "scale")
plotGenes(gene_idx[6], .experiment_name = "HAP4",
          .quartet = TRUE, .normalization = "scale")
plotGenes(gene_idx[7], .experiment_name = "HAP4",
          .quartet = TRUE, .normalization = "scale")
plotGenes(gene_idx[8], .experiment_name = "HAP4",
          .quartet = TRUE, .normalization = "scale")

# "weird trans"
gene_idx <- filter(plotdf, hybrid_cor > 0.8 &
                     dynamics != "conserved" &
                     cer_cor > 0.8 & par_cor > 0.8 &
                     cer == 2 & par == 1) |> select(gene_name) |> pull()
plotGenes(gene_idx[1], .experiment_name = "HAP4",
          .quartet = TRUE, .normalization = "scale")
plotGenes(gene_idx[2], .experiment_name = "HAP4",
          .quartet = TRUE, .normalization = "scale")
# misclassified as diverged in the parents, which will happen, there's only 2

# Scer dominant trans
gene_idx <- filter(plotdf, hybrid_cor > 0.8 &
                     cer_cor > 0.8 & par_cor < 0.5 &
                     cer == 2 & par == 1) |> select(gene_name) |> pull()
plotGenes(gene_idx, .experiment_name = "HAP4",
          .quartet = TRUE, .normalization = "scale")
# Spar dominant trans
gene_idx <- filter(plotdf, hybrid_cor > 0.8 &
                     cer_cor < 0.5 & par_cor > 0.8 &
                     cer == 2 & par == 1) |> select(gene_name) |> pull()
plotGenes(gene_idx, .experiment_name = "HAP4",
          .quartet = TRUE, .normalization = "scale")

# diverged genes Scer cor vs hybrid cor
# 1-2
ggplot(filter(plotdf, (cer == 1 & par == 2)),
       aes(x = hybrid_cor, y = cer_cor)) +
  geom_point()
ggplot(filter(plotdf, (cer == 1 & par == 2)),
       aes(x = hybrid_cor, y = par_cor)) +
  geom_point()
# 2-1
ggplot(filter(plotdf, (cer == 2 & par == 1)),
       aes(x = cer_cor, y = hybrid_cor)) +
  geom_point()
ggplot(filter(plotdf, (cer == 2 & par == 1)),
       aes(x = par_cor, y = hybrid_cor)) +
  geom_point()
# diverged genes with higher Scer cor than hybrid cor
gene_idx <- plotdf |> filter((cer == 2 & par == 1) &
                               hybrid_cor < cer_cor) |> 
  select(gene_name) |> pull()
plotGenes(.gene_idxs = gene_idx, 
          .quartet = TRUE, .experiment_name = experiment_name) 
# diverged genes with higher Scer cor than hybrid cor
gene_idx <- plotdf |> filter((cer == 2 & par == 1) &
                               hybrid_cor > cer_cor) |> 
  select(gene_name) |> pull()
plotGenes(.gene_idxs = gene_idx, 
          .quartet = TRUE, .experiment_name = experiment_name) 
gene_idx <- plotdf |> filter((cer == 2 & par == 1) &
                               hybrid_cor < cer_cor) |> 
  select(gene_name) |> pull()
plotGenes(.gene_idxs = gene_idx, 
          .quartet = TRUE, .experiment_name = experiment_name) 
# diverged genes with higher Scer cor than hybrid cor
gene_idx <- plotdf |> filter((cer == 2 & par == 1) &
                               hybrid_cor > par_cor) |> 
  select(gene_name) |> pull()
plotGenes(.gene_idxs = gene_idx, 
          .quartet = TRUE, .experiment_name = experiment_name) 
gene_idx <- plotdf |> filter((cer == 2 & par == 1) &
                               hybrid_cor < par_cor) |> 
  select(gene_name) |> pull()
plotGenes(.gene_idxs = gene_idx, 
          .quartet = TRUE, .experiment_name = experiment_name) 

# testing metrics
# random conserved gene, increasing group
gene_idx <- plotdf |> 
  filter(level == "conserved" & cer == 1 & par == 1) |> 
  select(gene_name) |> pull() |> sample(1)
plotdf |> filter(gene_name == gene_idx) |> 
  select(effect_size_species, pvalue_species,
         effect_size_allele, pvalue_allele)
plotGenes(.gene_idxs = gene_idx, .quartet = TRUE, .experiment_name = experiment_name) 

# random level gene
gene_idx <- "YGR192C"
gene_idx <- plotdf |> 
  filter(level == "diverged" & dynamics == "conserved") |> 
  select(gene_name) |> pull() |> sample(1)
gene_idx
plotGenes(.gene_idxs = gene_idx, .quartet = TRUE, .experiment_name = experiment_name) 
# YIL136W YKL216W, YGL028C, YNL195C, pure trans
# YFL030W, YGR192C, YDR111C, YLR179C, YDR007W, YIL107C, YHR071W, YJL217W, cis
# YDR461C-A, YMR107W, cis then trans
# YDR530C, cis, then trans, then cis

# random dynamics gene
gene_idx <- plotdf |> 
  filter(level == "conserved" & dynamics == "diverged") |> 
  select(gene_name) |> pull() |> sample(1)
plotGenes(.gene_idxs = gene_idx, .quartet = TRUE, .experiment_name = experiment_name)

# quantifying all genes
# cor
p <- ggplot(plotdf, aes(x = parent_cor, y = hybrid_cor)) + 
  geom_point(aes(color = group4)) +
  scale_color_discrete(type = levdyn_colordf$type,
                       labels = levdyn_colordf$labels,
                       limits = levdyn_colordf$limits) +
  theme_classic() + 
  xlim(c(-1, 1)) + 
  ylim(c(-1, 1)) +
  theme(legend.position = "none") +
  xlab("parental correlation") +
  ylab("hybrid correlation")
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/cor.pdf"),
    width = 4, height = 4)
ggMarginal(p, groupColour = TRUE, groupFill = FALSE)
dev.off()

#### Investigation: Lightning Bolt Heat/Cold expression ####
### Does Fay et al. 2023 alignment have more clear parental dominance?
load("data_files/Cleaned_Fay_Counts.RData")
load("data_files/Cleaned_Fay_Counts_Allele.RData")

common_conditions_tab_fay <- bind_rows(sample_info_fay, 
                                       sample_info_fay_allele) |> 
  select(condition, organism, allele) |> 
  unique() |> select(condition) |> table()
common_conditions_fay <- names(common_conditions_tab_fay)[common_conditions_tab_fay == 4]
info_fay <- sample_info_fay |> filter(condition %in% common_conditions_fay) |> 
  select(condition, experiment, time_point_str, time_point_num) |> 
  unique()
fay_hyc <- fay_allele[, sample_info_fay_allele$allele == "cer"]
fay_hyp <- fay_allele[, sample_info_fay_allele$allele == "par"]
info_hyb_fay <- sample_info_fay_allele |> 
  mutate(sample_name = gsub(pattern = "S[cp]",
                            replacement = "Hyb",
                            sample_name)) |> 
  select(sample_name, experiment) |> unique()
colnames(fay_hyc) <- gsub(pattern = "Sc",
                          replacement = "Hyb",
                          colnames(fay_hyc))
colnames(fay_hyp) <- gsub(pattern = "Sp",
                          replacement = "Hyb",
                          colnames(fay_hyp))
fay_hyc <- fay_hyc[, info_hyb_fay$sample_name]
fay_hyp <- fay_hyp[, info_hyb_fay$sample_name]
getCorelation("YDR229W", "Heat",
              .cts1 = fay_hyc,
              .cts2 = fay_hyp,
              .info1 = info_hyb_fay, 
              .info2 = info_hyb_fay)

faydf <- expand_grid(gene_name = rownames(fay),
                     experiment = c("Heat", "Cold"))
# hybrid correlations that include replicates
faydf$cor_hybrid <- map2(faydf$gene_name, faydf$experiment, 
                         getCorelation,
                         .cts1 = fay_hyc,
                         .cts2 = fay_hyp,
                         .info1 = info_hyb_fay,
                         .info2 = info_hyb_fay) |> unlist()
# %pct of Fay counts with strong hybrid corr
sum(faydf$cor_hybrid > cor_thresh, na.rm = TRUE)/nrow(faydf)
# finaldf actually has more strong hybrid correlations:
sum(finaldf$cor_hybrid[finaldf$experiment %in% c("Heat", "Cold")] > cor_thresh, na.rm = TRUE)/nrow(finaldf[finaldf$experiment %in% c("Heat", "Cold"),])

# collapsing counts for more correlations w/o replicates
# cer
fay_cer_collapsed <- map(common_conditions_fay, \(cond) {
  rowMeans(fay[,sample_info_fay$condition == cond &
                 sample_info_fay$allele == "cer",
               drop = FALSE])
}) |> purrr::reduce(.f = cbind)
colnames(fay_cer_collapsed) <- common_conditions_fay
rownames(fay_cer_collapsed) <- rownames(fay)
# par
fay_par_collapsed <- map(common_conditions_fay, \(cond) {
  rowMeans(fay[,sample_info_fay$condition == cond &
                 sample_info_fay$allele == "par",
               drop = FALSE])
}) |> purrr::reduce(.f = cbind)
colnames(fay_par_collapsed) <- common_conditions_fay
rownames(fay_par_collapsed) <- rownames(fay)
# hyc
fay_hyc_collapsed <- map(common_conditions_fay, \(cond) {
  rowMeans(fay_allele[,sample_info_fay_allele$condition == cond &
                        sample_info_fay_allele$allele == "cer",
                      drop = FALSE])
}) |> purrr::reduce(.f = cbind)
colnames(fay_hyc_collapsed) <- common_conditions_fay
rownames(fay_hyc_collapsed) <- rownames(fay)
# hyp
fay_hyp_collapsed <- map(common_conditions_fay, \(cond) {
  rowMeans(fay_allele[,sample_info_fay_allele$condition == cond &
                        sample_info_fay_allele$allele == "par",
                      drop = FALSE])
}) |> purrr::reduce(.f = cbind)
colnames(fay_hyp_collapsed) <- common_conditions_fay
rownames(fay_hyp_collapsed) <- rownames(fay)
getCorelation(sample(rownames(fay_hyc_collapsed),
                     size = 1), "Cold",
              .cts1 = fay_hyc_collapsed,
              .cts2 = fay_hyp_collapsed,
              .info1 = info_fay, 
              .info2 = info_fay)
# adding allele and parent correlations
faydf$cor_hybrid_noReps <- map2(faydf$gene_name, faydf$experiment, 
                                getCorelation,
                                .cts1 = fay_hyc_collapsed,
                                .cts2 = fay_hyp_collapsed,
                                .info1 = info_fay,
                                .info2 = info_fay) |> unlist()
# %pct of Fay counts with strong hybrid corr with reps
sum(faydf$cor_hybrid > cor_thresh, na.rm = TRUE)/nrow(faydf)
# %pct of Fay counts with strong hybrid corr without reps
sum(faydf$cor_hybrid_noReps > cor_thresh, na.rm = TRUE)/nrow(faydf) # more
# still not as many as finaldf:
sum(finaldf$cor_hybrid[finaldf$experiment %in% c("Heat", "Cold")] > cor_thresh, na.rm = TRUE)/nrow(finaldf[finaldf$experiment %in% c("Heat", "Cold"),])
# other collapsed cors
faydf$cor_parents <- map2(faydf$gene_name, faydf$experiment, 
                          getCorelation,
                          .cts1 = fay_cer_collapsed,
                          .cts2 = fay_par_collapsed,
                          .info1 = info_fay,
                          .info2 = info_fay) |> unlist()
faydf$cor_scer <- map2(faydf$gene_name, faydf$experiment, 
                       getCorelation,
                       .cts1 = fay_cer_collapsed,
                       .cts2 = fay_hyc_collapsed,
                       .info1 = info_fay,
                       .info2 = info_fay) |> unlist()
faydf$cor_spar <- map2(faydf$gene_name, faydf$experiment, 
                       getCorelation,
                       .cts1 = fay_par_collapsed,
                       .cts2 = fay_hyp_collapsed,
                       .info1 = info_fay,
                       .info2 = info_fay) |> unlist()

faydf <- left_join(faydf, select(filter(finaldf,
                                        experiment == "Heat"),
                                 gene_name, dynamics, cer, par),
                   by = "gene_name")
plotdf <- faydf |> filter(experiment == "Heat" & 
                            cor_hybrid > cor_thresh &
                            cor_parents < 0.25 &
                            dynamics != "conserved")

ggplot(plotdf, aes(x = cor_scer, y = cor_spar)) +
  geom_point(aes(color = paste(cer, par),
                 shape = paste(cer, par)))
# Conclusion: alignment doesn't matter, same pattern

# 3) Lightning bolt
# TODO: identify genes in cer, par, hyc, and hyp with a lightning bolt shape
# defined as 15min is lowest expr and 30min is highest of the 4 timepoints
# also check for reverse lightning bolt---30min is lowest expr and 15min is highest
# Questions:
#   1) Is it just in Heat or also Cold?
#   2) How do non lightning bolt genes adjust their expression? Are there an equal
# number of reverse lightning bolt genes, or are there different sets of genes
# peaking at 30min and dipping at 15min?
#   3) Are hybrid lightning bolts more common among genes with less correlation in parents?
#       or are there plenty of conserved dynamics genes with hybrid lightning bolts too?
#   4) Are hybrid lightning bolts more likely to be seen in both alleles?

# the genes where both hybrid alleles correlate much more with
# each other than either parent by default form this lightning bolt:
test_exp <- "Heat"
lightning_gene_idxs <- finaldf |> filter(cor_scer < cor_thresh &
                                           cor_spar < cor_thresh &
                                           cor_hybrid > cor_thresh &
                                           experiment == test_exp) |> 
  select(gene_name) |> pull()
plotGenes(lightning_gene_idxs, .experiment_name = test_exp,
          .quartet = TRUE)

# What does the parental dominance plot look like in the
# absence of these genes?
plotdf <- finaldf |> filter(experiment == "Heat" & 
                              cor_hybrid > cor_thresh &
                              cor_parents < 0.25 &
                              !(gene_name %in% lightning_gene_idxs) &
                              dynamics != "conserved")

ggplot(plotdf, aes(x = cor_scer, y = cor_spar)) +
  geom_point(aes(color = paste(cer, par),
                 shape = paste(cer, par)))
# Ah yes b/c the lightning idxs were based on having low cor_scer and cor_spar
# this just looks like we removed them
# what if we identify lightning genes differently?

### Identifying lightning genes by timepoint rank
getExprRanks <- function(.gene_name, .experiment, .cts, .info) {
  expr_vec <- .cts[.gene_name, .info$experiment == .experiment]
  output <- rank(-expr_vec) |> list()
  return(output)
}
# # tests for getExprRanks
# getExprRanks(.gene_name = "YGR192C", .experiment = "LowN",
#              .cts = collapsed$cer, .info = info)
# getExprRanks(.gene_name = "YGR192C", .experiment = "LowN",
#              .cts = collapsed$par, .info = info)
# # plot to verify
# plotGenes("YGR192C", .experiment_name = "LowN")
finaldf$cer_ranks <- map2(finaldf$gene_name, finaldf$experiment,
                          .f = getExprRanks, .cts = collapsed$cer,
                          .info = info)
finaldf$par_ranks <- map2(finaldf$gene_name, finaldf$experiment,
                          .f = getExprRanks, .cts = collapsed$par,
                          .info = info)
finaldf$hyc_ranks <- map2(finaldf$gene_name, finaldf$experiment,
                          .f = getExprRanks, .cts = collapsed_allele$cer,
                          .info = info_allele)
finaldf$hyp_ranks <- map2(finaldf$gene_name, finaldf$experiment,
                          .f = getExprRanks, .cts = collapsed_allele$par,
                          .info = info_allele)
# lightning bolts cer
checkLightningBolt <- function(.rank_list) {
  ranks <- unlist(.rank_list)
  if (length(ranks) != 4) {
    stop("Lightning bolts have 4 points!\n")
  }
  if (ranks[2] == 4 & ranks[3] == 1) {
    return("lightning")
  }
  if (ranks[2] == 1 & ranks[3] == 4) {
    return("reverse")
  }
  if (ranks[3] == 1) {
    return("peakTP3")
  }
  else {
    return("none")
  }
}
lightningdf <- finaldf |> filter(experiment %in% c("Heat", "Cold"))
lightningdf$cer_lightning <- map(lightningdf$cer_ranks, 
                                 .f = checkLightningBolt) |> 
  unlist()
lightningdf$par_lightning <- map(lightningdf$par_ranks, 
                                 .f = checkLightningBolt) |> 
  unlist()
lightningdf$hyc_lightning <- map(lightningdf$hyc_ranks, 
                                 .f = checkLightningBolt) |> 
  unlist()
lightningdf$hyp_lightning <- map(lightningdf$hyp_ranks, 
                                 .f = checkLightningBolt) |> 
  unlist()

peak3_gene_idxs_by_rank <- lightningdf |> 
  filter(experiment == "Heat" &
           (hyc_lightning %in% c("lightning", "peakTP3")) &
           (hyp_lightning %in% c("lightning", "peakTP3"))) |> 
  select(gene_name) |> pull()
length(intersect(peak3_gene_idxs_by_rank,
                 lightning_gene_idxs))
length(setdiff(peak3_gene_idxs_by_rank,
               lightning_gene_idxs)) 
# most of rank genes are contained in the lightning genes defined
# by low parent correlation
length(setdiff(lightning_gene_idxs,
               peak3_gene_idxs_by_rank))
# but there are plenty of low parent correlation genes
# that aren't lightning
plotdf <- finaldf |> filter(experiment == "Heat" & 
                              cor_hybrid > cor_thresh &
                              cor_parents < 0.25 &
                              #!(gene_name %in% peak3_gene_idxs_by_rank) &
                              dynamics != "conserved")

ggplot(plotdf, aes(x = cor_scer, y = cor_spar)) +
  geom_point(aes(color = paste(cer, par),
                 shape = paste(cer, par)))

plotGenes(.gene_idxs = peak3_gene_idxs_by_rank,
          .experiment_name = "Heat", .quartet = TRUE)
plotGenes(.gene_idxs = intersect(peak3_gene_idxs_by_rank,
                                 lightning_gene_idxs),
          .experiment_name = "Heat", .quartet = TRUE)
plotGenes(.gene_idxs = setdiff(peak3_gene_idxs_by_rank,
                               lightning_gene_idxs),
          .experiment_name = "Heat", .quartet = TRUE)
# whether or not it has low scer/spar correlation doesn't
# matter too much for the peak3
plotGenes(.gene_idxs = setdiff(lightning_gene_idxs,
                               peak3_gene_idxs_by_rank),
          .experiment_name = "Heat", .quartet = TRUE)
# but without peak3, the expression shapes get much less distinct
plotGenes(.gene_idxs = sample(setdiff(lightning_gene_idxs,
                                      peak3_gene_idxs_by_rank),
                              size = 1),
          .experiment_name = "Heat", .quartet = TRUE)

# The non lightning genes are still genes that peak at TP3

### Checking lightning gene metrics
# parent lightning bolts:
lightningdf |> select(cer_lightning, par_lightning) |> table()
# hybrid lightning bolts:
lightningdf |> select(hyc_lightning, hyp_lightning) |> table()
# heat vs cold:
lightningdf |> filter(cer_lightning == "lightning" |
                        par_lightning == "lightning" |
                        hyc_lightning == "lightning" | 
                        hyp_lightning == "lightning") |> 
  select(experiment) |> table()
# Cold more common in parent lightning:
lightningdf |> filter(cer_lightning == "lightning" |
                        par_lightning == "lightning") |> 
  select(experiment) |> table()
# Heat more common in hybrid lightning:
lightningdf |> filter(hyc_lightning == "lightning" |
                        hyp_lightning == "lightning") |> 
  select(experiment) |> table()
# What clusters are parent lightning genes in?
yesParentsTable <- lightningdf |> filter((cer_lightning == "lightning" |
                                            par_lightning == "lightning") &
                                           experiment == "Cold") |> 
  select(cer, par) |> table()
noParentsTable <- lightningdf |> filter((cer_lightning != "lightning" &
                                           par_lightning != "lightning") &
                                          experiment == "Cold") |> 
  select(cer, par) |> table()
round(yesParentsTable/sum(yesParentsTable), digits = 2)
round(noParentsTable/sum(noParentsTable), digits = 2)
# very few cluster 1 genes it looks like?
lightningdf |> 
  mutate(cluster_1_in_one = (cer == 1 | par == 1),
         lightning_in_one = (cer_lightning == "lightning" |
                               par_lightning == "lightning")) |> 
  select(cluster_1_in_one, lightning_in_one) |> table() # yup
#  What clusters are hybrid lightning genes in?
yesHybridTable <- lightningdf |> filter((hyc_lightning == "lightning" |
                                           hyp_lightning == "lightning") &
                                          experiment == "Heat") |> 
  select(cer, par) |> table()
noHybridTable <- lightningdf |> filter((hyc_lightning != "lightning" &
                                          hyp_lightning != "lightning") &
                                         experiment == "Heat") |> 
  select(cer, par) |> table()
round(yesHybridTable/sum(yesHybridTable), digits = 2)
round(noHybridTable/sum(noHybridTable), digits = 2)
# No obvious differences in cer/par cluster percentages
# hyc/hyp clusters reveal that hybrid lightning genes are much more likely to be
# in hyc/hyp cluster 1 than 2



#### Separating cis and trans diverging level and dynamics genes ####

# Based on correlation and log2 fold change metrics,
# we separate out the genes
# with cis divergence from trans

### trans/level
# level divergers, cluster 1, cer up
# main plot
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/lfc_level_trans.pdf"),
    width = 3, height = 3)
ggplot(filter(plotdf, group4 == "diverged level, conserved dynamics"),
       aes(x = effect_size_species, y = effect_size_allele)) + 
  geom_point(color = levdyn_colordf$type[levdyn_colordf$limits == "diverged level, conserved dynamics"]) +
  theme_classic() + 
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  theme(legend.position = "none") +
  xlab("log2 fold change, parents") +
  ylab("log2 fold change, hybrid") +
  geom_rect(xmin = 0, xmax = 5, ymin = -0.5, ymax = 0.5, color = "black", alpha = 0)
dev.off()
# corresponding avgExpr plot
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/lfc_level_trans_avgExpr.pdf"),
    width = 3, height = 3)
gene_idxs <- plotdf |> filter(cer == 1 & par == 1 &
                                group4 == "diverged level, conserved dynamics" &
                                effect_size_allele > -0.5 &
                                effect_size_allele < 0.5 &
                                effect_size_species > 1) |>
  select(gene_name) |> pull() |> sample(1)
gene_idxs <-  "YDR070C"
annotate_figure(plotGenes(.gene_idxs = gene_idxs, .quartet = TRUE, .experiment_name = experiment_name),
                top = paste(gene_idxs, "example gene,\n level diverging in trans"))
dev.off()

### cis/level
# main plot
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/lfc_level_cis.pdf"),
    width = 3, height = 3)
ggplot(filter(plotdf, group4 == "diverged level, conserved dynamics"),
       aes(x = effect_size_species, y = effect_size_allele)) + 
  geom_point(color = levdyn_colordf$type[levdyn_colordf$limits == "diverged level, conserved dynamics"]) +
  theme_classic() + 
  theme(legend.position = "none") +
  xlab("log2 fold change, parents") +
  ylab("log2 fold change, hybrid") +
  geom_rect(xmin = 0, xmax = 10.5, ymin = 1, ymax = 10.5, color = "black", alpha = 0)
dev.off()
# corresponding avgExpr plot
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/lfc_level_cis_avgExpr.pdf"),
    width = 3, height = 3)
gene_idxs <- plotdf |> filter(cer == 1 & par == 1 &
                                group4 == "diverged level, conserved dynamics" &
                                effect_size_species > 0 &
                                effect_size_allele > 1) |>
  select(gene_name) |> pull()
annotate_figure(plotGenes(.gene_idxs = gene_idxs, .quartet = TRUE, .experiment_name = experiment_name),
                top = paste0(length(gene_idxs), " genes,\nlevel diverging in cis"))
dev.off()

### dynamics/cis
# dynamics divergers, 2-1
# main plot
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/cor_dynamics_cis.pdf"),
    width = 3, height = 3)
ggplot(filter(plotdf, group4 == "conserved level, diverged dynamics"),
       aes(x = parent_cor, y = hybrid_cor)) + 
  geom_point(color = levdyn_colordf$type[levdyn_colordf$limits == "conserved level, diverged dynamics"]) +
  theme_classic() + 
  xlim(c(-1, 1)) + 
  ylim(c(-1, 1)) + 
  theme(legend.position = "none") +
  xlab("correlation, parents") +
  ylab("correlation, hybrid") +
  geom_rect(xmin = -1, xmax = 1, ymin = -1, ymax = -0.75, color = "black", alpha = 0)
dev.off()
# corresponding avgExpr plot
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/cor_dynamics_cis_avgExpr.pdf"),
    width = 3, height = 3)
gene_idxs <- plotdf |> filter(cer == 2 & par == 1 &
                                group4 == "conserved level, diverged dynamics" &
                                hybrid_cor < -0.5) |>
  select(gene_name) |> pull()
length(gene_idxs)
gene_idx <- sample(gene_idxs, 1)
gene_idx <- "YNR061C"
annotate_figure(plotGenes(.gene_idxs = gene_idx, .quartet = TRUE, .experiment_name = experiment_name),
                top = paste(gene_idx, "example gene,\n dynamics diverging in cis"))
dev.off()

### dynamics/trans
# dynamics divergers, 2-1
# main plot
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/cor_dynamics_trans.pdf"),
    width = 3, height = 3)
ggplot(filter(plotdf, group4 == "conserved level, diverged dynamics"),
       aes(x = parent_cor, y = hybrid_cor)) + 
  geom_point(color = levdyn_colordf$type[levdyn_colordf$limits == "conserved level, diverged dynamics"]) +
  theme_classic() + 
  xlim(c(-1, 1)) + 
  ylim(c(-1, 1)) + 
  theme(legend.position = "none") +
  xlab("correlation, parents") +
  ylab("correlation, hybrid") +
  geom_rect(xmin = -1, xmax = 1, ymin = 0, ymax = 1, color = "black", alpha = 0)
dev.off()
# corresponding avgExpr plot
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/cor_dynamics_trans_avgExpr.pdf"),
    width = 3, height = 3)
gene_idxs <- plotdf |> filter(cer == 2 & par == 1 &
                                group4 == "conserved level, diverged dynamics" &
                                hybrid_cor > -0.5) |>
  select(gene_name) |> pull()
annotate_figure(plotGenes(.gene_idxs = gene_idxs, .quartet = TRUE, .experiment_name = experiment_name),
                top = paste0(length(gene_idxs), " genes,\n dynamics diverging in trans"))
dev.off()

### level and dynamics divergers, cis x trans
# main plot, correlation, trans
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/cor_leveldyn_trans.pdf"),
    width = 3, height = 3)
ggplot(filter(plotdf, group4 == "diverged level and dynamics"),
       aes(x = parent_cor, y = hybrid_cor)) + 
  geom_point(color = levdyn_colordf$type[levdyn_colordf$limits == "diverged level and dynamics"]) +
  theme_classic() + 
  theme(legend.position = "none") +
  xlab("correlation, parents") +
  ylab("correlation, hybrid") +
  xlim(c(-1, 1)) +
  ylim(c(-1, 1)) +
  geom_rect(xmin = -1, xmax = 1, ymin = 0, ymax = 1, color = "black", alpha = 0)
dev.off()
# corresponding avgExpr plot
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/cor_leveldyn_cistrans_avgExpr.pdf"),
    width = 3, height = 3)
# these guys are doing too many different things for avgExpr to make sense
gene_idxs <- plotdf |> filter(group4 == "diverged level and dynamics" &
                                hybrid_cor > 0) |>
  select(gene_name) |> pull() |> sample(1)
# YNR050C is a good example gene---goes from cluster 1 - 2 in parents, but it's trans and cer dominant in, cer and hyc consistently higher
gene_idxs <- "YNR050C"
annotate_figure(plotGenes(.gene_idxs = gene_idxs, .quartet = TRUE, .experiment_name = experiment_name),
                top = paste(gene_idxs, "example gene,\nlevel diverging in cis,\ndynamics diverging in trans"))
dev.off()
# also looking for a pure trans lev dyn gene
gene_idxs <- plotdf |> filter(group4 == "diverged level and dynamics" &
                                hybrid_cor > 0) |>
  select(gene_name) |> pull() |> sample(1)
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/cor_leveldyn_trans_avgExpr.pdf"),
    width = 3, height = 3)
gene_idxs <- "YNL007C"
annotate_figure(plotGenes(.gene_idxs = gene_idxs, .quartet = TRUE, .experiment_name = experiment_name),
                top = paste(gene_idxs, "example gene,\nlevel and dynamics \ndiverging in trans"))
dev.off()
# main plot, correlation, outliers
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/cor_leveldyn_outliers.pdf"),
    width = 3, height = 3)
ggplot(filter(plotdf, group4 == "diverged level and dynamics"),
       aes(x = parent_cor, y = hybrid_cor)) + 
  geom_point(color = levdyn_colordf$type[levdyn_colordf$limits == "diverged level and dynamics"]) +
  theme_classic() + 
  theme(legend.position = "none") +
  xlab("correlation, parents") +
  ylab("correlation, hybrid") +
  xlim(c(-1, 1)) +
  ylim(c(-1, 1)) +
  geom_rect(xmin = -1, xmax = 1, ymin = -1, ymax = -0.75, color = "black", alpha = 0)
dev.off()
# corresponding avgExpr plot
gene_idxs <- plotdf |> filter(group4 == "diverged level and dynamics" &
                                hybrid_cor < -0.75) |>
  select(gene_name) |> pull() |> sample(1)
gene_idxs <- "YDL201W"
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/cor_leveldyn_outlier1_avgExpr.pdf"),
    width = 3, height = 3)
annotate_figure(plotGenes(.gene_idxs = gene_idxs, .quartet = TRUE, .experiment_name = experiment_name),
                top = paste(gene_idxs, "example gene\nlevel and dynamics\ndiverging in cis"))
dev.off()

########################## archive #########################
#### Other measures of parent/hybrid allele similarity ####
# calculate Area Between Curves, abc
# takes mean difference between each row in cts1 - each row in cts2
# @input: 2 counts matrices, rows are genes cols are samples, same number of rows, cols are in same order
# @output: vector of length nrow(counts matrix) giving each gene's area between curves
# # if lfc and cor aren't enough, We'll need to make this
# # take "rectangles" between curves to account for timepoints being different distances apart
# # (currently early timepoints have disproportionate weight)
# calculateAreaBetweenCurves <- function(.cts1, .cts2) {
#   abc_mat <- (.cts1 - .cts2)/pmax(.cts1, .cts2)
#   return(rowMeans(abc_mat, na.rm = TRUE))
# } 
# parent_ABC <- calculateAreaBetweenCurves(collapsed$cer[,common_conditions],
#                                          collapsed$par[,common_conditions])
# hybrid_ABC <- calculateAreaBetweenCurves(collapsed_allele$cer[,common_conditions],
#                                          collapsed_allele$par[,common_conditions])
# parent_meanDiff <- map(c(1:nrow(collapsed$cer)), \(i) {
#   abs(mean(collapsed$cer[i, common_conditions]) - 
#         mean(collapsed$par[i, common_conditions]))
# }) |> unlist()
# hybrid_meanDiff <- map(c(1:nrow(collapsed$cer)), \(i) {
#   abs(mean(collapsed_allele$cer[i, common_conditions]) - 
#         mean(collapsed_allele$par[i, common_conditions]))
# }) |> unlist()
# ABCdf <- tibble(gene_name = rownames(collapsed$cer),
#                 parent_abc = parent_ABC,
#                 hybrid_abc = hybrid_ABC,
#                 parent_cor,
#                 hybrid_cor,
#                 parent_meanDiff,
#                 hybrid_meanDiff)
# meanDiff
# p <- ggplot(plotdf, aes(x = log2(parent_meanDiff), y = log2(hybrid_meanDiff))) + 
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point(aes(color = group)) +
#   scale_color_discrete(type = levdyn_colordf$type,
#                        labels = levdyn_colordf$labels,
#                        breaks = levdyn_colordf$breaks) +
#   theme_classic() + 
#   theme(legend.position = "none") +
#   xlim(c(log2(min(c(plotdf$parent_meanDiff, plotdf$hybrid_meanDiff))),
#          log2(max(c(plotdf$parent_meanDiff, plotdf$hybrid_meanDiff))))) +
#   ylim(c(log2(min(c(plotdf$parent_meanDiff, plotdf$hybrid_meanDiff))),
#          log2(max(c(plotdf$parent_meanDiff, plotdf$hybrid_meanDiff))))) +
#   xlab("log2(| parental mean difference |)") +
#   ylab("log2(| hybrid mean difference |)")
# pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/meanDiff.pdf",
#     width = 4, height = 4)
# ggMarginal(p, groupColour = TRUE, groupFill = FALSE)
# dev.off()