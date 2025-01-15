sapply(c("tidyverse", "scales", "ggpubr"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Barkai_data_analysis/")

# Goal: addressing a series of questions of how individual genes are diverging

load(file = "data_files/env_spec_DEdf.RData")
load(file = "data_files/cis_trans_df.RData")
load(file = "data_files/Cleaned_Barkai_Data.RData")
load(file = "data_files/Cleaned_Barkai_Data_AlleleSpecific.RData")

# pre-processing
countsPerMillion <- function(.cts) {
  librarySizes <- colSums(.cts)
  output <- apply(.cts, 1, function(x) {
    normalized <- (x/librarySizes)*1e6
    return(round(normalized))
  })
  return(t(output)) # For some unhinged reason, a vector output of apply across ROWS forms the COLUMNS of a new matrix
}
counts <- countsPerMillion(counts)
counts_allele <- countsPerMillion(counts_allele)

### Question 1: Do cis and trans-diverging genes have similar magnitudes of expression divergence (effect size of allele-associated beta coefficient) and are they expressed at similar levels?
table(spaldf$divergence_mode) # first of all, we can see there are more cis than trans-diverging genes
DivergenceModeColors <- hue_pal()(5) %>% as.character()

# Are their effect sizes the same?
# TODO: run marginal effects to convert effect sizes into actual estimates
# Wittkopp-style plot
ggplot(spaldf, aes(x = species_effect_size, y = allele_effect_size)) + 
  geom_point(aes(color = divergence_mode)) +
  scale_color_manual(values = c("cis"=DivergenceModeColors[1], "compensated"=DivergenceModeColors[2], "conserved"=DivergenceModeColors[3], "trans"=DivergenceModeColors[4], "weirdos"=DivergenceModeColors[5]),
                     breaks = c("cis", "compensated", "conserved", "trans", "weirdos")) + ggtitle("Effect Sizes") + xlab("Species Effect Size") + ylab("Allele Effect Size") + theme_classic()
# mostly about the same effect sizes, but there are a handful of cis-diverging genes with real strong effects!

# cis vs trans species effect size 
plot1 <- ggplot(filter(spaldf, divergence_mode %in% c("cis", "trans")), aes(x = divergence_mode, y = species_effect_size)) + geom_jitter(aes(color = divergence_mode)) +
  scale_color_manual(values = c("cis"=DivergenceModeColors[1], "compensated"=DivergenceModeColors[2], "conserved"=DivergenceModeColors[3], "trans"=DivergenceModeColors[4], "weirdos"=DivergenceModeColors[5]),
                     breaks = c("cis", "compensated", "conserved", "trans", "weirdos")) + 
  ggtitle("Species Effect Sizes") + xlab("Divergence Mode") + ylab("Species Effect Size") + theme_classic()
# allele is obviously smaller in trans
plot2 <- ggplot(filter(spaldf, divergence_mode %in% c("cis", "trans")), aes(x = divergence_mode, y = allele_effect_size)) + geom_jitter(aes(color = divergence_mode)) +
  scale_color_manual(values = c("cis"=DivergenceModeColors[1], "compensated"=DivergenceModeColors[2], "conserved"=DivergenceModeColors[3], "trans"=DivergenceModeColors[4], "weirdos"=DivergenceModeColors[5]),
                     breaks = c("cis", "compensated", "conserved", "trans", "weirdos")) + 
  ggtitle("Allele Effect Sizes") + xlab("Divergence Mode") + ylab("Allele Effect Size") + theme_classic()

ggarrange(plot1, plot2, nrow = 2, ncol = 1, common.legend = TRUE)

# Conclusion: cis has stronger effect. And interestingly there are only a few genes with very strong effect

# Are they expressed at similar levels?
calculateMedianExpression <- function(.gene_name, .organism = c("cer", "par", "hyb")) {
  cts <- counts[.gene_name, sample_info$organism == .organism]
  return(median(cts))
}
# Not sure we'll even need this b/c it's really just a less-good version of our DE models
# more involved because we're condition-matching samples before calculating expression difference
calculateMedianExpressionDifferenceCerMinusPar <- function(.gene_name, .mode = c("parents", "hybrid"), .experiment = c("LowN", "CC", "HAP4", "LowPi")) {
  if (.mode == "parents") {
    cts <- counts
    info <- sample_info
  }
  if (.mode == "hybrid") {
    cts <- counts_allele
    info <- sample_info_allele
  }
  cts_cer <- cts[.gene_name, (info$allele == "cer") & (info$experiment %in% .experiment)]
  si_cer <- info[(info$allele == "cer") & (info$experiment %in% .experiment),]
  cerdf <- bind_cols(tibble(expr = cts_cer),
                     si_cer)
  cts_par <- cts[.gene_name, (info$allele == "par") & (info$experiment %in% .experiment)]
  si_par <- info[(info$allele == "par") & (info$experiment %in% .experiment),]
  pardf <- bind_cols(tibble(expr = cts_par),
                     si_par)
  ctsdf <- bind_rows(cerdf, pardf) %>% mutate(non_unique_sample_name = paste(condition, well_flask_ID)) %>% select(expr, non_unique_sample_name, allele) %>% pivot_wider(id_cols = c(non_unique_sample_name), values_from = expr, names_from = allele)
  return(median(ctsdf$cer - ctsdf$par, na.rm = TRUE))
}
calculateMedianExpressionDifferenceCerMinusPar("YGR192C","parents", "CC")
calculateMedianExpressionDifferenceCerMinusPar(.gene_name = "YGR192C", .mode = "hybrid") # omit experiment to default to all 4

spaldf$median_expr_cer <- map(spaldf$gene_name, calculateMedianExpression, .organism = "cer") %>% unlist()
spaldf$median_expr_par <- map(spaldf$gene_name, calculateMedianExpression, .organism = "par") %>% unlist()
spaldf$median_expr_hyb <- map(spaldf$gene_name, calculateMedianExpression, .organism = "hyb") %>% unlist()

plotdf <- pivot_longer(spaldf, cols = c(median_expr_cer, median_expr_par, median_expr_hyb))

ggplot(filter(plotdf, !is.na(divergence_mode)), aes(x = divergence_mode, y = log(value))) + geom_boxplot(aes(fill = name)) +
  scale_fill_manual(values = c("median_expr_cer"="gold", "median_expr_par"="mediumorchid", "median_expr_hyb"="chartreuse3"),
                    labels = c("cer", "hyb", "par")) + xlab("divergence mode") + ylab("log(median expression)") + theme_classic()

# Conclusion: There really is no bias for highly expressed genes for any divergence mode. Other than genes being called as DE in general (cis/trans/comp/weird vs conserved), but it's very slight

### Question 2: Are environmentally-variable genes more likely to be diverging? (Do any divergence modes have a higher density of genes with highly variable expression? Are genes DE in specific environments overrepresented in genes DE between species/alleles?)
# First calculating Tau, a measure of variability more robust to expression level than variance
calculateTau <- function(.gene_name, .organism = c("cer", "par", "hyb")) {
  expr <- counts[.gene_name, sample_info$organism == .organism]
  m <- length(expr)
  maxexpr <- max(expr)
  unscaled_tau <- (1 - expr/maxexpr) %>% sum(na.rm = TRUE)
  return(unscaled_tau/(m - 1))
}
calculateTau(.gene_name = "YGR192C", .organism = "cer")

spaldf$tau_cer <- map(spaldf$gene_name, calculateTau, .organism = "cer") %>% unlist()
spaldf$tau_par <- map(spaldf$gene_name, calculateTau, .organism = "par") %>% unlist()
spaldf$tau_hyb <- map(spaldf$gene_name, calculateTau, .organism = "hyb") %>% unlist()

plotdf <- pivot_longer(spaldf, cols = c(tau_cer, tau_par, tau_hyb))

ggplot(filter(plotdf, !is.na(divergence_mode)), aes(x = divergence_mode, y = value)) + geom_boxplot(aes(fill = name)) +
  scale_fill_manual(values = c("tau_cer"="gold", "tau_par"="mediumorchid", "tau_hyb"="chartreuse3"),
                    labels = c("cer", "hyb", "par")) + xlab("divergence mode") + ylab("expression variability (tau)") + theme_classic()
# Conclusion: now there isn't even a difference between genes called conserved vs diverged
# But interestingly, Tau is lower across the board fo hybrids (except the weirdos)

# What about for a different measure of variability?
calculateScaledVariance <- function(.gene_name, .organism = c("cer", "par", "hyb")) {
  expr <- counts[.gene_name, sample_info$organism == .organism]
  v <- var(expr)
  m <- mean(expr)
  return(v/m)
}
calculateScaledVariance("YGR192C", "cer")

spaldf$scaledvar_cer <- map(spaldf$gene_name, calculateScaledVariance, .organism = "cer") %>% unlist()
spaldf$scaledvar_par <- map(spaldf$gene_name, calculateScaledVariance, .organism = "par") %>% unlist()
spaldf$scaledvar_hyb <- map(spaldf$gene_name, calculateScaledVariance, .organism = "hyb") %>% unlist()

plotdf <- pivot_longer(spaldf, cols = c(scaledvar_cer, scaledvar_par, scaledvar_hyb))

ggplot(filter(plotdf, !is.na(divergence_mode)), aes(x = divergence_mode, y = log(value))) + geom_boxplot(aes(fill = name)) +
  scale_fill_manual(values = c("scaledvar_cer"="gold", "scaledvar_par"="mediumorchid", "scaledvar_hyb"="chartreuse3"),
                    labels = c("cer", "hyb", "par")) + xlab("divergence mode") + ylab("expression variability (var/mean)") + theme_classic()
# Not as pronounced, but same pattern---hyb is lower than other two

# Are genes that are identified as cis or trans-diverging only diverging in certain environments? I.E. how useful was it to have expr info from 4 different experiments?
# Look at median expr difference in each of the 4 experiments individually to see how many genes have similar expr difference across experiments
# Also look at random individual genes to see if median tells the whole story
cer_minus_par_df <- map(spaldf$gene_name, function(g) {
  cat("currently processing", g, which(spaldf$gene_name == g), "/", nrow(spaldf), "\n")
  output <- tibble(gene_name = g,
                   experiment = rep(c("all", "CC", "LowN", "LowPi", "HAP4"), 2),
                   mode = c(rep("parents", 5), rep("hybrid", 5)))
  output$cer_minus_par <- map(c(1:nrow(output)), function(i) {
    if (output[i, "experiment"] == "all") {
      cer_minus_par <- calculateMedianExpressionDifferenceCerMinusPar(.gene_name = g, .mode = output[i,"mode"])
    }
    if (output[i, "experiment"] != "all") {
      cer_minus_par <- calculateMedianExpressionDifferenceCerMinusPar(.gene_name = g, .mode = output[i,"mode"], .experiment = output[i,"experiment"])
    }
    return(cer_minus_par)
  }) %>% unlist()
  return(output)
}) %>% Reduce(f = bind_rows)

save(cer_minus_par_df, file = "data_files/cer_minus_par_df.RData")

gene_idx <- sample(spaldf$gene_name, 1)
ggplot(filter(plotdf, gene_name == gene_idx), aes(x = experiment, y = cer_minus_par, group = mode)) + geom_line(aes(color = mode)) + theme_classic() + ggtitle("Expression difference between species (cer - par)") + xlab("experiment") + ylab("Median expression\n difference (cer - par)")
test <- plotdf %>% group_by(gene_name) %>% summarise(scaled_range = (max(cer_minus_par) - min(cer_minus_par))/mean(abs(cer_minus_par) + 1),
                                                     unscaled_range = max(cer_minus_par) - min(cer_minus_par))
ggplot(test, aes(x = unscaled_range)) + geom_histogram(bins = 100) + theme_classic()

test <- left_join(test, select(spaldf, c(gene_name, divergence_mode)), by = "gene_name")

ggplot(filter(test, !is.na(divergence_mode)), aes(x = divergence_mode, y = scaled_range)) + geom_boxplot(aes(fill = divergence_mode)) + theme_classic() + xlab("divergence mode") + ylab("cer-par expression\n difference range")

sum(test$unscaled_range > 100)
mean(test$unscaled_range)
median(test$unscaled_range)


# unfortunately the 4601-gene plots crash R
# ggplot(filter(plotdf, mode == "parents"), aes(x = experiment, y = cer_minus_par, group = gene_name)) + geom_line(aes(color = gene_name)) + theme_classic() + ggtitle("Expression difference between species (cer - par)") + xlab("experiment") + ylab("Median expression\n difference (cer - par)")
# ggplot(filter(plotdf, mode == "hybrid"), aes(x = experiment, y = cer_minus_par, group = gene_name)) + geom_line(aes(color = gene_name)) + theme_classic() + ggtitle("Expression difference between hybrid alleles (cer - par)") + xlab("experiment") + ylab("Median expression\n difference (cer - par)")

# TODO: if it seems worth it, run env-specific allele models to see if the same set of genes is ID'd
load(file = "data_files/cis_trans_by_env_DEdf.RData") # TODO: make this
