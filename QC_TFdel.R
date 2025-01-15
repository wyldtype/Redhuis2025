#### QC: comparing power to detect DE genes upon TF deletion in each species ####
sapply(c("dplyr", "purrr", "tidyr", "ggpubr", "readr", "data.table", "ggplot2", "data.table", "msir", "WGCNA", "energy"), require, character.only=TRUE)
setwd("Documents/Wittkopp_Lab/networks/DDivergence/Redhuis2024/")
library(MASS, include.only = "glm.nb")

# set constants
coef_thresh <- 0.25 # minimum magnitude, so if it's < -coef_thresh, we also want that to be detected
p_thresh <- 1e-5

# Goal #1: recapitulate power issue, where paradoxus has many more DE genes but moreso in Redhuis alignment

#### Data cleaning ####
# Redhuis alignment
load(file = "data_files/Cleaned_Redhuis_Data.RData")
load(file = "data_files/Cleaned_Redhuis_Data_AlleleSpecific.RData")

#### Condition-matched counts/info for parents and hybrids ####

### condition-matching samples so that expression similarity can be calculated
counts_cer <- counts[, sample_info$organism == "cer"] |> t()
info_cer <- sample_info %>% filter(sample_name %in% rownames(counts_cer))
counts_par <- counts[, sample_info$organism == "par"] |> t()
info_par <- sample_info %>% filter(sample_name %in% rownames(counts_par))
counts_hyc <- counts_allele[, sample_info_allele$organism == "hyb" & sample_info_allele$allele == "cer"] |> t()
info_hyc <- sample_info_allele %>% filter(sample_name %in% rownames(counts_hyc))
counts_hyp <- counts_allele[, sample_info_allele$organism == "hyb" & sample_info_allele$allele == "par"] |> t()
info_hyp <- sample_info_allele %>% filter(sample_name %in% rownames(counts_hyp))

# Filtering out conditions (combinations of genotype, experiment, and timepoint) that are not represented in cer, par, and hyb
common_conditions <- bind_rows(info_cer, info_par, info_hyc, info_hyp) %>% 
  select(condition, allele, organism) |> unique() %>% 
  pivot_wider(id_cols = c("condition"), names_from = c("allele", "organism"), values_from = "organism") %>% 
  drop_na() %>% 
  select(condition) |> pull()

# filtering down to common conditions
info_cer <- filter(info_cer, condition %in% common_conditions) %>% arrange(experiment, condition, well_flask_ID)
info_par <- filter(info_par, condition %in% common_conditions) %>% arrange(experiment, condition, well_flask_ID)
info_hyc <- filter(info_hyc, condition %in% common_conditions) %>% arrange(experiment, condition, well_flask_ID)
info_hyp <- filter(info_hyp, condition %in% common_conditions) %>% arrange(experiment, condition, well_flask_ID)

# checking conditions are the same (these should both be 1)
sum(unique(info_cer$condition) == unique(info_par$condition))/length(unique(info_cer$condition))
sum(unique(info_hyc$condition) == unique(info_hyp$condition))/length(unique(info_hyc$condition))

# paring down counts to match condition-matched info dfs
counts_cer <- counts_cer[info_cer$sample_name,] # remember genes are COLUMNS in WGCNA
counts_par <- counts_par[info_par$sample_name,]
counts_hyc <- counts_hyc[info_hyc$sample_name,]
counts_hyp <- counts_hyp[info_hyp$sample_name,]

# splitting off TFdel dataset b/c replicates are useful for glm to tell how much unexplained variation there is
infos_TFdel <- list(cer = info_cer[info_cer$experiment == "LowN",],
                    par = info_par[info_par$experiment == "LowN",]) 
counts_TFdel <- list(cer = counts_cer[info_cer$experiment == "LowN",], 
                     par = counts_par[info_par$experiment == "LowN",])
infos_TFdel_allele <- list(cer = info_hyc[info_hyc$experiment == "LowN",],
                           par = info_hyp[info_hyp$experiment == "LowN",]) 
counts_TFdel_allele <- list(cer = counts_hyc[info_hyc$experiment == "LowN",], 
                            par = counts_hyp[info_hyp$experiment == "LowN",])

# TF filtering: TF must have all at least 2 replicates 
# as these lead to inaccurate glm.nb models below
# cer
goodTable_cer <- infos_TFdel$cer |> filter(genotype != "WT") |> 
  select(condition) |> table()
cer_vec <- sapply(goodTable_cer, \(x) {return(x >= 2)})
goodConditions_cer <- rownames(goodTable_cer)[cer_vec]
# par
goodTable_par <- infos_TFdel$par |> filter(genotype != "WT") |> 
  select(condition) |> table()
par_vec <- sapply(goodTable_par, \(x) {return(x >= 2)})
goodConditions_par <- rownames(goodTable_par)[par_vec]
# hyc
goodTable_hyc <- infos_TFdel_allele$cer |> filter(genotype != "WT") |> 
  select(condition) |> table()
hyc_vec <- sapply(goodTable_hyc, \(x) {return(x >= 2)})
goodConditions_hyc <- rownames(goodTable_hyc)[hyc_vec]
# hyp
goodTable_hyp <- infos_TFdel_allele$par |> filter(genotype != "WT") |> 
  select(condition) |> table()
hyp_vec <- sapply(goodTable_hyp, \(x) {return(x >= 2)})
goodConditions_hyp <- rownames(goodTable_hyp)[hyp_vec]
# all
goodConditions <- reduce(list(goodConditions_cer, goodConditions_par, goodConditions_hyc, goodConditions_hyp), 
                         .f = intersect) # also ensuring the same set of TFs between species and hybrid
# verifying these TFs do have all 6 timepoints x replicates in all 4 species/alleles
goodTable_cer[goodConditions] |> min()
goodTable_par[goodConditions] |> min()
goodTable_hyc[goodConditions] |> min()
goodTable_hyp[goodConditions] |> min() # should all have min 2

# filtering
redhuis_cts <- list(cer = counts_TFdel$cer[(infos_TFdel$cer$condition %in% goodConditions) | (infos_TFdel$cer$genotype == "WT"),],
                    par = counts_TFdel$par[(infos_TFdel$par$condition %in% goodConditions) | (infos_TFdel$par$genotype == "WT"),])

redhuis_info <- list(cer = infos_TFdel$cer[(infos_TFdel$cer$condition %in% goodConditions) | (infos_TFdel$cer$genotype == "WT"),],
                     par = infos_TFdel$par[(infos_TFdel$par$condition %in% goodConditions) | (infos_TFdel$par$genotype == "WT"),])

redhuis_cts_allele <- list(cer = counts_TFdel_allele$cer[(infos_TFdel_allele$cer$condition %in% goodConditions) | (infos_TFdel_allele$cer$genotype == "WT"),],
                           par = counts_TFdel_allele$par[(infos_TFdel_allele$par$condition %in% goodConditions) | (infos_TFdel_allele$par$genotype == "WT"),])

redhuis_info_allele <- list(cer = infos_TFdel_allele$cer[(infos_TFdel_allele$cer$condition %in% goodConditions) | (infos_TFdel_allele$cer$genotype == "WT"),],
                            par = infos_TFdel_allele$par[(infos_TFdel_allele$par$condition %in% goodConditions) | (infos_TFdel_allele$par$genotype == "WT"),])

rm(sample_info, sample_info_allele, counts, counts_allele,
   counts_cer, counts_par, counts_hyc, counts_hyp)

# repeat for Barkai
load(file = "data_files/Cleaned_Barkai_Data.RData")
load(file = "data_files/Cleaned_Barkai_Data_AlleleSpecific.RData")

#### Condition-matched counts/info for parents and hybrids ####

### condition-matching samples so that expression similarity can be calculated
counts_cer <- counts[, sample_info$organism == "cer"] |> t()
info_cer <- sample_info %>% filter(sample_name %in% rownames(counts_cer))
counts_par <- counts[, sample_info$organism == "par"] |> t()
info_par <- sample_info %>% filter(sample_name %in% rownames(counts_par))
counts_hyc <- counts_allele[, sample_info_allele$organism == "hyb" & sample_info_allele$allele == "cer"] |> t()
info_hyc <- sample_info_allele %>% filter(sample_name %in% rownames(counts_hyc))
counts_hyp <- counts_allele[, sample_info_allele$organism == "hyb" & sample_info_allele$allele == "par"] |> t()
info_hyp <- sample_info_allele %>% filter(sample_name %in% rownames(counts_hyp))

# Filtering out conditions (combinations of genotype, experiment, and timepoint) that are not represented in cer, par, and hyb
common_conditions <- bind_rows(info_cer, info_par, info_hyc, info_hyp) %>% 
  select(condition, allele, organism) |> unique() %>% 
  pivot_wider(id_cols = c("condition"), names_from = c("allele", "organism"), values_from = "organism") %>% 
  drop_na() %>% 
  select(condition) |> pull()

# filtering down to common conditions
info_cer <- filter(info_cer, condition %in% common_conditions) %>% arrange(experiment, condition, well_flask_ID)
info_par <- filter(info_par, condition %in% common_conditions) %>% arrange(experiment, condition, well_flask_ID)
info_hyc <- filter(info_hyc, condition %in% common_conditions) %>% arrange(experiment, condition, well_flask_ID)
info_hyp <- filter(info_hyp, condition %in% common_conditions) %>% arrange(experiment, condition, well_flask_ID)

# checking conditions are the same (these should both be 1)
sum(unique(info_cer$condition) == unique(info_par$condition))/length(unique(info_cer$condition))
sum(unique(info_hyc$condition) == unique(info_hyp$condition))/length(unique(info_hyc$condition))

# paring down counts to match condition-matched info dfs
counts_cer <- counts_cer[info_cer$sample_name,] # remember genes are COLUMNS in WGCNA
counts_par <- counts_par[info_par$sample_name,]
counts_hyc <- counts_hyc[info_hyc$sample_name,]
counts_hyp <- counts_hyp[info_hyp$sample_name,]

# splitting off TFdel dataset b/c replicates are useful for glm to tell how much unexplained variation there is
infos_TFdel <- list(cer = info_cer[info_cer$experiment == "LowN",],
                    par = info_par[info_par$experiment == "LowN",]) 
counts_TFdel <- list(cer = counts_cer[info_cer$experiment == "LowN",], 
                     par = counts_par[info_par$experiment == "LowN",])
infos_TFdel_allele <- list(cer = info_hyc[info_hyc$experiment == "LowN",],
                           par = info_hyp[info_hyp$experiment == "LowN",]) 
counts_TFdel_allele <- list(cer = counts_hyc[info_hyc$experiment == "LowN",], 
                            par = counts_hyp[info_hyp$experiment == "LowN",])

# TF filtering: TF must have all at least 2 replicates 
# filtering (goodConditions determined in Redhuis, and we verified that barkai included additional TFs without omitting any)
barkai_cts <- list(cer = counts_TFdel$cer[(infos_TFdel$cer$condition %in% goodConditions) | (infos_TFdel$cer$genotype == "WT"),],
                    par = counts_TFdel$par[(infos_TFdel$par$condition %in% goodConditions) | (infos_TFdel$par$genotype == "WT"),])

barkai_info <- list(cer = infos_TFdel$cer[(infos_TFdel$cer$condition %in% goodConditions) | (infos_TFdel$cer$genotype == "WT"),],
                     par = infos_TFdel$par[(infos_TFdel$par$condition %in% goodConditions) | (infos_TFdel$par$genotype == "WT"),])

barkai_cts_allele <- list(cer = counts_TFdel_allele$cer[(infos_TFdel_allele$cer$condition %in% goodConditions) | (infos_TFdel_allele$cer$genotype == "WT"),],
                           par = counts_TFdel_allele$par[(infos_TFdel_allele$par$condition %in% goodConditions) | (infos_TFdel_allele$par$genotype == "WT"),])

barkai_info_allele <- list(cer = infos_TFdel_allele$cer[(infos_TFdel_allele$cer$condition %in% goodConditions) | (infos_TFdel_allele$cer$genotype == "WT"),],
                            par = infos_TFdel_allele$par[(infos_TFdel_allele$par$condition %in% goodConditions) | (infos_TFdel_allele$par$genotype == "WT"),])

rm(list = setdiff(ls(), c("redhuis_cts", "redhuis_cts_allele","redhuis_info", "redhuis_info_allele",
                          "barkai_cts", "barkai_cts_allele","barkai_info", "barkai_info_allele",
                          "coef_thresh", "p_thresh")))

#### model fitting ####
factorizeGenotypeAndTimepoint <- function(.info) {
  if (length(unique(.info$genotype)) > 1) {
    .info$genotype <- as.factor(.info$genotype) %>% relevel(ref = "WT")
  }
  if (length(unique(.info$time_point_str)) > 1) {
    timepoints <- .info$time_point_str %>% unique()
    reftime <- timepoints[which.min(parse_number(timepoints))] # can't just do 0 because LowPi has a -5 fml
    .info$time_point_str <- as.factor(.info$time_point_str) %>% relevel(ref = reftime)
  }
  return(.info)
}
redhuis_info <- lapply(redhuis_info, factorizeGenotypeAndTimepoint)
redhuis_info_allele <- lapply(redhuis_info_allele, factorizeGenotypeAndTimepoint)
barkai_info <- lapply(barkai_info, factorizeGenotypeAndTimepoint)
barkai_info_allele <- lapply(barkai_info_allele, factorizeGenotypeAndTimepoint)

# verify that counts are counts per million
test_cts <- barkai_cts_allele$cer
test_rowIdx <- sample(c(1:nrow(test_cts)), 1)
test_colIdx <- sample(c(1:ncol(test_cts)), 1)
test_count <- test_cts[test_rowIdx, test_colIdx]
((test_count/rowSums(test_cts, na.rm = TRUE)[test_rowIdx])*1e6) %>% round() # what it should be
test_cts[test_rowIdx, test_colIdx] # what it is currently, should be the same, +/- 1 for rounding

# # first exploring detecting which TF deletions a gene is DE in with a positive control
# # Another good example is YLR307W (CDA1), a sporulation gene, has crazy high expression specifically in SUM1delete, a sporulation repressor, so that's cool, but we don't have SUM1data with replicates
# gene_idx <- "YGR192C" # YGR192C (TDH3) in GCR2 delete
# genedf <- bind_cols(tibble(expr = redhuis_cts$cer[,gene_idx]), redhuis_info$cer)
# ggplot(genedf, aes(x = genotype, y = expr)) + geom_point() + scale_x_discrete(breaks = NULL)
# modnb <- glm.nb(expr ~ time_point_str + genotype + time_point_str:genotype, data = genedf, link = log)
# summary(modnb)
# genedf$p_hat <- modnb$fitted.values
# plotdf <- pivot_longer(genedf, cols = c(p_hat, expr))
# ggplot(data = plotdf, aes(x = genotype, y = value)) +
#   geom_point(aes(color = name, shape = as.factor(time_point_str))) +
#   scale_x_discrete(breaks = NULL)
# coeffs <- coef(modnb)[grepl("^genotype", names(coef(modnb)))]
# pvals <- summary(modnb)$coefficients[,4][grepl("^genotype", names(coef(modnb)))]
# sigGenotypes <- names(coeffs[abs(coeffs) > 0.5 & pvals < p_thresh]) %>% gsub(pattern = "genotype", replacement = "")
# genedf$sig <- genedf$genotype %in% sigGenotypes
# genedf$sig[genedf$genotype == "WT"] <- NA
# ggplot(data = genedf, aes(x = genotype, y = expr)) +
#   geom_point(aes(color = sig)) +
#   theme(axis.text.x = element_text(angle = 90))

# now making it a function so we can apply it to all genes
# @input: gene name, counts matrix (where rows are samples, cols are genes),
#         and info dataframe for LowN experiment in one species
# @ouput: list of three elements: 1) pvalues of each coefficient (of length nDeletions)
#                                 2) coefficients (effect sizes) of each TFdel's effect on the focal gene
#                                 3) name of which deletion each coefficient/pvalue corresponds to

getTFdelModResultsByGene <- function(.gene_idx, .cts, .info, .mod_func) {
  gdf <- bind_cols(tibble(expr = .cts[,.gene_idx]), .info)
  gdf$genotype <- factor(gdf$genotype) |> relevel(ref = "WT")
  gdf$time_point_str <- factor(gdf$time_point_str) |> relevel(ref = "0 h, YPD")
  get1genemod <- function(del) {
    gdf2gen <- gdf |> filter(genotype %in% c("WT", del))
    # filter out timepoints not present in deletion
    goodTimes <- gdf2gen |> filter(genotype != "WT") |> select(time_point_str) |> pull() |> unique()
    gdf2gen <- gdf2gen |> filter(time_point_str %in% goodTimes)
    # fitting model
    m <- tryCatch({
      .mod_func(.data = gdf2gen)
    }, error = function(e) {
      return(NA)
    })
    if (!is.na(m)[1]) {
      pval <- summary(m)$coefficients[,4][paste0("genotype", del)]
      coef <- coef(m)[paste0("genotype", del)]
      return(list("pval" = as.numeric(pval),
                  "coef" = as.numeric(coef),
                  "deletion" = gsub("delete", "", as.character(del))))
    }
    else {
      return(list("pval" = NA,
                  "coef" = NA,
                  "deletion" = NA))
    }
  }
  output <- map_dfr(.x = setdiff(unique(gdf$genotype), "WT"),
                    .f = get1genemod)
  output$gene_name <- .gene_idx
  return(output)
}
# tests for getGenesAllTFdelCoeffs
# test_result_cer <- getTFdelModResultsByGene("YGR192C", redhuis_cts$cer, redhuis_info$cer) # TDH3 should be way down when its direct regulator, GCR2, is deleted
# test_result_cer$coef[test_result_cer$deletion == "GCR2"]
# test_result_cer$pval[test_result_cer$deletion == "GCR2"]
# test_result_par <- getTFdelModResultsByGene("YGR192C", redhuis_cts$par, redhuis_info$par)
# test_result_par$coef[test_result_cer$deletion == "GCR2"]
# test_result_par$pval[test_result_cer$deletion == "GCR2"]
# # barkai alignment
# test_result_cer <- getTFdelModResultsByGene("YGR192C", barkai_cts$cer, barkai_info$cer) # TDH3 should be way down when its direct regulator, GCR2, is deleted
# test_result_cer$coef[test_result_cer$deletion == "GCR2"]
# test_result_cer$pval[test_result_cer$deletion == "GCR2"]
# test_result_par <- getTFdelModResultsByGene("YGR192C", barkai_cts$par, barkai_info$par)
# test_result_par$coef[test_result_cer$deletion == "GCR2"]
# test_result_par$pval[test_result_cer$deletion == "GCR2"]
# 
# # YKR034W is DAL80, so this is a negative control
# test_result_cer <- getTFdelModResultsByGene("YKR034W", redhuis_cts$cer, redhuis_info$cer) # TDH3 should be way down when its direct regulator, GCR2, is deleted
# test_result_cer$coef[test_result_cer$deletion == "DAL80"]
# test_result_cer$pval[test_result_cer$deletion == "DAL80"]
# test_result_par <- getTFdelModResultsByGene("YKR034W", redhuis_cts$par, redhuis_info$par)
# test_result_par$coef[test_result_cer$deletion == "DAL80"]
# test_result_par$pval[test_result_cer$deletion == "DAL80"]
# # barkai alignment
# test_result_cer <- getTFdelModResultsByGene("YKR034W", barkai_cts$cer, barkai_info$cer) # TDH3 should be way down when its direct regulator, GCR2, is deleted
# test_result_cer$coef[test_result_cer$deletion == "DAL80"]
# test_result_cer$pval[test_result_cer$deletion == "DAL80"]
# test_result_par <- getTFdelModResultsByGene("YKR034W", barkai_cts$par, barkai_info$par)
# test_result_par$coef[test_result_cer$deletion == "DAL80"]
# test_result_par$pval[test_result_cer$deletion == "DAL80"]

mod_poisglm <- function(.data) {
  m <- glm(expr ~ genotype + time_point_str, data = .data, family = poisson(link = "log"))
  return(m)
}
mod_poisglm_int <- function(.data) {
  m <- glm(expr ~ genotype * time_point_str, data = .data, family = poisson(link = "log"))
  return(m)
}
mod_negbin <- function(.data) {
  m <- glm.nb(expr ~ genotype + time_point_str, data = .data, link = log)
  return(m)
}
mod_negbin_int <- function(.data) {
  m <- glm.nb(expr ~ genotype * time_point_str, data = .data, link = log)
  return(m)
}
mod_notime <- function(.data) {
  m <- glm(expr ~ genotype, data = .data, family = poisson(link = "log"))
  return(m)
} # this is a terrible idea but should have effect sizes exactly the same as standardized mean difference, so it's a good control

# fitting models
TFdeldfs_pois <- map2(.x = list(redhuis_cts$cer, redhuis_cts$par, redhuis_cts_allele$cer, redhuis_cts_allele$par,
                                barkai_cts$cer, barkai_cts$par, barkai_cts_allele$cer, barkai_cts_allele$par),
                      .y = list(redhuis_info$cer, redhuis_info$par, redhuis_info_allele$cer, redhuis_info_allele$par,
                                barkai_info$cer, barkai_info$par, barkai_info_allele$cer, barkai_info_allele$par),
                      \(x, y) {
                        output <- map_dfr(colnames(x), \(g) {
                          cat("working on", g, which(colnames(x) == g), "/", ncol(x), "\n")
                          return(getTFdelModResultsByGene(g, x, y, .mod_func = mod_poisglm))
                        })
                      })
TFdeldfs_poisint <- map2(.x = list(redhuis_cts$cer, redhuis_cts$par, redhuis_cts_allele$cer, redhuis_cts_allele$par,
                                   barkai_cts$cer, barkai_cts$par, barkai_cts_allele$cer, barkai_cts_allele$par),
                         .y = list(redhuis_info$cer, redhuis_info$par, redhuis_info_allele$cer, redhuis_info_allele$par,
                                   barkai_info$cer, barkai_info$par, barkai_info_allele$cer, barkai_info_allele$par),
                         \(x, y) {
                           output <- map_dfr(colnames(x), \(g) {
                             cat("working on", g, which(colnames(x) == g), "/", ncol(x), "\n")
                             return(getTFdelModResultsByGene(g, x, y, .mod_func = mod_poisglm_int))
                           })
                         })
TFdeldfs_negbin <- map2(.x = list(redhuis_cts$cer, redhuis_cts$par, redhuis_cts_allele$cer, redhuis_cts_allele$par,
                                  barkai_cts$cer, barkai_cts$par, barkai_cts_allele$cer, barkai_cts_allele$par),
                        .y = list(redhuis_info$cer, redhuis_info$par, redhuis_info_allele$cer, redhuis_info_allele$par,
                                  barkai_info$cer, barkai_info$par, barkai_info_allele$cer, barkai_info_allele$par),
                        \(x, y) {
                          output <- map_dfr(colnames(x), \(g) {
                            cat("working on", g, which(colnames(x) == g), "/", ncol(x), "\n")
                            return(getTFdelModResultsByGene(g, x, y, .mod_func = mod_negbin))
                          })
                        })
TFdeldfs_negbinint <- map2(.x = list(redhuis_cts$cer, redhuis_cts$par, redhuis_cts_allele$cer, redhuis_cts_allele$par,
                                  barkai_cts$cer, barkai_cts$par, barkai_cts_allele$cer, barkai_cts_allele$par),
                        .y = list(redhuis_info$cer, redhuis_info$par, redhuis_info_allele$cer, redhuis_info_allele$par,
                                  barkai_info$cer, barkai_info$par, barkai_info_allele$cer, barkai_info_allele$par),
                        \(x, y) {
                          output <- map_dfr(colnames(x), \(g) {
                            cat("working on", g, which(colnames(x) == g), "/", ncol(x), "\n")
                            return(getTFdelModResultsByGene(g, x, y, .mod_func = mod_negbin_int))
                          })
                        })
TFdeldfs_notime <- map2(.x = list(redhuis_cts$cer, redhuis_cts$par, redhuis_cts_allele$cer, redhuis_cts_allele$par,
                                  barkai_cts$cer, barkai_cts$par, barkai_cts_allele$cer, barkai_cts_allele$par),
                        .y = list(redhuis_info$cer, redhuis_info$par, redhuis_info_allele$cer, redhuis_info_allele$par,
                                  barkai_info$cer, barkai_info$par, barkai_info_allele$cer, barkai_info_allele$par),
                        \(x, y) {
                          output <- map_dfr(colnames(x), \(g) {
                            cat("working on", g, which(colnames(x) == g), "/", ncol(x), "\n")
                            return(getTFdelModResultsByGene(g, x, y, .mod_func = mod_notime))
                          })
                        })

names(TFdeldfs_pois) <- c("redhuisdf_cer", "redhuisdf_par", "redhuisdf_hyc", "redhuisdf_hyp",
                          "barkaidf_cer", "barkaidf_par", "barkaidf_hyc", "barkaidf_hyp")
names(TFdeldfs_poisint) <- c("redhuisdf_cer", "redhuisdf_par", "redhuisdf_hyc", "redhuisdf_hyp",
                          "barkaidf_cer", "barkaidf_par", "barkaidf_hyc", "barkaidf_hyp")
names(TFdeldfs_negbin) <- c("redhuisdf_cer", "redhuisdf_par", "redhuisdf_hyc", "redhuisdf_hyp",
                            "barkaidf_cer", "barkaidf_par", "barkaidf_hyc", "barkaidf_hyp")
names(TFdeldfs_negbinint) <- c("redhuisdf_cer", "redhuisdf_par", "redhuisdf_hyc", "redhuisdf_hyp",
                               "barkaidf_cer", "barkaidf_par", "barkaidf_hyc", "barkaidf_hyp")
names(TFdeldfs_notime) <- c("redhuisdf_cer", "redhuisdf_par", "redhuisdf_hyc", "redhuisdf_hyp",
                            "barkaidf_cer", "barkaidf_par", "barkaidf_hyc", "barkaidf_hyp")

# converting natural log fold change to log2 fold change
# e^x= 2^y
# log2(e^x) = y
TFdeldfs_pois <- lapply(TFdeldfs_pois, \(x) {
  x$coef <- log2(exp(x$coef))
  return(x)
})
TFdeldfs_poisint <- lapply(TFdeldfs_poisint, \(x) {
  x$coef <- log2(exp(x$coef))
  return(x)
})
TFdeldfs_negbin <- lapply(TFdeldfs_negbin, \(x) {
  x$coef <- log2(exp(x$coef))
  return(x)
})
TFdeldfs_negbinint <- lapply(TFdeldfs_negbinint, \(x) {
  x$coef <- log2(exp(x$coef))
  return(x)
})
TFdeldfs_notime <- lapply(TFdeldfs_notime, \(x) {
  x$coef <- log2(exp(x$coef))
  return(x)
})

# standardized mean difference (https://cran.r-project.org/web/packages/TOSTER/vignettes/SMD_calcs.html)
# with dots colored by p_value significance threshold
# standarized mean difference is calculated as:
# (m2 - m1)/s_av
# s_av = sqrt((sd1^2 + sd2^2)/2) (square root of the mean variance)
# where m1 and sd1 are the mean and sd of counts from WT 
# and m2 and sd2 are the same for the TF deletion
# (m2 - m1 allows the sign to match the effect sizes when WT is reference:
# + means the gene increased expr upon deletion, - means it decreased)
getStandardizedMeanDiffByGene <- function(.gene_idx, .cts, .info) {
  cat("working on", .gene_idx, 
      which(.gene_idx == colnames(redhuis_cts$cer)), 
      "/", ncol(redhuis_cts$cer), "\n")
  gdf <- bind_cols(tibble(expr = .cts[,.gene_idx]), .info)
  gdf$genotype <- factor(gdf$genotype) |> relevel(ref = "WT")
  get1geneSMD <- function(del) {
    gdf2gen <- gdf |> filter(genotype %in% c("WT", del))
    # filter out timepoints not present in deletion
    goodTimes <- gdf2gen |> filter(genotype != "WT") |> select(time_point_str) |> pull() |> unique()
    gdf2gen <- gdf2gen |> filter(time_point_str %in% goodTimes)
    # calculating SMD
    wt_vec <- filter(gdf2gen, genotype == "WT") |> select(expr) |> pull()
    del_vec <- filter(gdf2gen, genotype != "WT") |> select(expr) |> pull()
    m1 <- mean(wt_vec, na.rm = TRUE)
    sd1 <- sd(wt_vec, na.rm = TRUE)
    m2 <- mean(del_vec, na.rm = TRUE)
    sd2 <- sd(del_vec, na.rm = TRUE)
    s_av <- sqrt((sd1^2 + sd2^2)/2)
    smd <- (m2-m1)/s_av
    return(list("smd" = smd,
                "deletion" = gsub("delete", "", del)))
  }
  output <- map_dfr(.x = setdiff(unique(gdf$genotype), "WT"),
                    .f = get1geneSMD)
  output$gene_name <- .gene_idx
  return(output)
}
# # tests for getStandardizedMeanDiffByGene
# # cer
# getStandardizedMeanDiffByGene("YGR192C", redhuis_cts$cer, redhuis_info$cer) |> filter(deletion == "GCR2")
# TFdeldfs_pois$redhuisdf_cer |> filter(deletion == "GCR2" & gene_name == "YGR192C")
# TFdeldfs_poisint$redhuisdf_cer |> filter(deletion == "GCR2" & gene_name == "YGR192C")
# TFdeldfs_negbin$redhuisdf_cer |> filter(deletion == "GCR2" & gene_name == "YGR192C")
# TFdeldfs_negbinint$redhuisdf_cer |> filter(deletion == "GCR2" & gene_name == "YGR192C") # the no interaction term ones are significantly different than with interaction
# # par
# getStandardizedMeanDiffByGene("YGR192C", redhuis_cts$par, redhuis_info$par) |> filter(deletion == "GCR2")
# TFdeldfs_pois$redhuisdf_par |> filter(deletion == "GCR2" & gene_name == "YGR192C")
# TFdeldfs_poisint$redhuisdf_par |> filter(deletion == "GCR2" & gene_name == "YGR192C")
# TFdeldfs_negbin$redhuisdf_par |> filter(deletion == "GCR2" & gene_name == "YGR192C")
# TFdeldfs_negbinint$redhuisdf_par |> filter(deletion == "GCR2" & gene_name == "YGR192C") # pars are more consistent between models. Is this a common pattern?

GeneNames <- colnames(redhuis_cts$cer)
DeletionNames <- setdiff(unique(redhuis_info$cer$genotype), "WT") |> gsub(pattern = "delete", replacement = "")
SMD_cer <- map(GeneNames, getStandardizedMeanDiffByGene, 
               .cts = redhuis_cts$cer,
               .info = redhuis_info$cer) |> reduce(.f = bind_rows)
SMD_par <- map(GeneNames, getStandardizedMeanDiffByGene, 
               .cts = redhuis_cts$par,
               .info = redhuis_info$par) |> reduce(.f = bind_rows)
SMDdf <- tibble(gene_name = rep(GeneNames, length(DeletionNames)),
                deletion = rep(DeletionNames, length(GeneNames))) |> 
  left_join(y = rename(SMD_cer, "cer"="smd"), by = c("gene_name", "deletion")) |> 
  left_join(y = rename(SMD_par, "par"="smd"), by = c("gene_name", "deletion"))


# saving
save(TFdeldfs_pois, TFdeldfs_poisint,
     TFdeldfs_negbin, TFdeldfs_negbinint,
     TFdeldfs_notime,
     SMDdf, file = "data_files/QC_TFdel.RData")
load(file = "data_files/QC_TFdel.RData")


#### effect size distributions by TF ####
# do any TFs have more nsig than others within each dataset?
plotlist <- vector(mode = "list", length = length(TFdeldfs_pois))
for (i in c(1:length(TFdeldfs_pois))) {
  df <- TFdeldfs_pois[[i]]
  plotdf <- df |> filter(abs(coef) > 0.5 & pval < p_thresh) |> 
    group_by(deletion)
  p <- ggplot(plotdf, aes(x = deletion)) +
    geom_bar(aes(fill = deletion)) +
    ggtitle(names(TFdeldfs_pois)[i]) +
    theme(legend.position = "none", 
          axis.text.x = element_text(angle = 90)) +
    ylim(c(0, 1600))
  plotlist[[i]] <- p
}
ggarrange(plotlist = plotlist, ncol = 4, nrow = 2)

# both barkai and redhuis have the discrepancy of par having ~2x as many DE as cer
# hybrids don't, so this isn't an alignment issue
# no obvious single TFs that are driving the pattern
# barkai vs redhuis alingment have slight differences in power, 
# but are mostly similar, especially hybrid 
# (except TEC1, which has more DE in redhuis hyc and hyp)

# but at some point, I didn't have this issue

#### comparing models ####

# overall n DE gene/TF combinations
# cer pois, no interaction
sum(TFdeldfs_pois$redhuisdf_cer$pval < p_thresh & 
      abs(TFdeldfs_pois$redhuisdf_cer$coef) > 0.25)
# par pois, no interaction
sum(TFdeldfs_pois$redhuisdf_par$pval < p_thresh & 
      abs(TFdeldfs_pois$redhuisdf_par$coef) > 0.25) # wait that's not that many more
# cer pois, with interaction
sum(TFdeldfs_poisint$redhuisdf_cer$pval < p_thresh & 
      abs(TFdeldfs_poisint$redhuisdf_cer$coef) > 0.25)
# par pois, with interaction
sum(TFdeldfs_poisint$redhuisdf_par$pval < p_thresh & 
      abs(TFdeldfs_poisint$redhuisdf_par$coef) > 0.25)
# extreme values only
# cer pois, no interaction
sum(TFdeldfs_pois$redhuisdf_cer$pval < p_thresh & 
      abs(TFdeldfs_pois$redhuisdf_cer$coef) > 2)
# par pois, no interaction
sum(TFdeldfs_pois$redhuisdf_par$pval < p_thresh & 
      abs(TFdeldfs_pois$redhuisdf_par$coef) > 2) # ah now it's over twice as many
# cer pois, with interaction
sum(TFdeldfs_poisint$redhuisdf_cer$pval < p_thresh & 
      abs(TFdeldfs_poisint$redhuisdf_cer$coef) > 2)
# par pois, with interaction
sum(TFdeldfs_poisint$redhuisdf_par$pval < p_thresh & 
      abs(TFdeldfs_poisint$redhuisdf_par$coef) > 2)
# cer negbin, no interaction
sum(TFdeldfs_negbin$redhuisdf_cer$pval < p_thresh & 
      abs(TFdeldfs_negbin$redhuisdf_cer$coef) > 2)
# par negbin, no interaction
sum(TFdeldfs_negbin$redhuisdf_par$pval < p_thresh & 
      abs(TFdeldfs_negbin$redhuisdf_par$coef) > 2, na.rm = TRUE)
# cer negbin, with interaction
sum(TFdeldfs_negbinint$redhuisdf_cer$pval < p_thresh & 
      abs(TFdeldfs_negbinint$redhuisdf_cer$coef) > 2)
# par negbin, with interaction
sum(TFdeldfs_negbinint$redhuisdf_par$pval < p_thresh & 
      abs(TFdeldfs_negbinint$redhuisdf_par$coef) > 2, na.rm = TRUE)

# conclusions: no obvious solution, but par seems to have more extreme effect sizes,
# especially in poisson models and models without interaction terms

### density plots of effect sizes in cer versus par
# SMD plot for ground-truth reference
# all SMDs, regardless of size
plotdf <- SMDdf |> pivot_longer(cols = c("cer", "par"),
                                names_to = "species",
                                values_to = "smd") 
ggplot(plotdf, aes(x = smd)) + 
  geom_density(aes(fill = species), alpha = 0.5) # par has stronger effects in both directions, not just lower tail

# now models
# pois, no interaction
plotdf <- TFdeldfs_pois$redhuisdf_cer |> 
  #filter(pval < p_thresh) |> 
  select(gene_name, deletion, coef) |> 
  mutate(species = "cer")
plotdf <- TFdeldfs_pois$redhuisdf_par |> 
  #filter(pval < p_thresh) |> 
  select(gene_name, deletion, coef) |> 
  mutate(species = "par") |> 
  bind_rows(y = plotdf)
ggplot(plotdf, aes(x = (-log(abs(coef)))*sign(coef))) + 
  geom_density(aes(fill = species), alpha = 0.5) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1)
# okayyy there is a bias for paradoxus to have more negative strong effect sizes
# this becomes overwhelming for stronger than effect sizes less than -1:
ggplot(plotdf, aes(x = coef)) + 
  geom_density(aes(fill = species), alpha = 0.5) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) +
  ylim(c(0, 0.1)) + xlim(c(-5, 5))

# pois, with interaction
plotdf <- TFdeldfs_poisint$redhuisdf_cer |> 
  filter(pval < p_thresh) |> 
  select(gene_name, deletion, coef) |> 
  mutate(species = "cer")
plotdf <- TFdeldfs_poisint$redhuisdf_par |> 
  filter(pval < p_thresh) |> 
  select(gene_name, deletion, coef) |> 
  mutate(species = "par") |> 
  bind_rows(y = plotdf)
ggplot(plotdf, aes(x = coef)) + 
  geom_density(aes(fill = species), alpha = 0.5) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1)

# negbin, no interaction
plotdf <- TFdeldfs_negbin$redhuisdf_cer |> 
  filter(pval < p_thresh) |> 
  select(gene_name, deletion, coef) |> 
  mutate(species = "cer")
plotdf <- TFdeldfs_negbin$redhuisdf_par |> 
  filter(pval < p_thresh) |> 
  select(gene_name, deletion, coef) |> 
  mutate(species = "par") |> 
  bind_rows(y = plotdf)
ggplot(plotdf, aes(x = coef)) + 
  geom_density(aes(fill = species), alpha = 0.5) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1)

# negbin, with interaction
plotdf <- TFdeldfs_negbinint$redhuisdf_cer |> 
  filter(pval < p_thresh) |> 
  select(gene_name, deletion, coef) |> 
  mutate(species = "cer")
plotdf <- TFdeldfs_negbinint$redhuisdf_par |> 
  filter(pval < p_thresh) |> 
  select(gene_name, deletion, coef) |> 
  mutate(species = "par") |> 
  bind_rows(y = plotdf)
ggplot(plotdf, aes(x = coef)) + 
  geom_density(aes(fill = species), alpha = 0.5) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1)


#### finding a better coef thresh cutoff ####

# seems like I'm being too stringent with my effect size cutoff
# use SMD of ~1 to figure out what the corresponding coef_thresh would be
# and use that instead
# (or look back at those coef vs smd plots and see where coef starts to "fan out"
# i.e. becomes unrelated to smd)

# control: no accounting for time, should be equivalent to smd
plotdf <- left_join(SMDdf, 
                    rename(TFdeldfs_notime$redhuisdf_cer,
                           c("coef_cer"="coef",
                             "pval_cer"="pval")),
                    by = join_by(gene_name, deletion)) |> 
  left_join(y = rename(TFdeldfs_notime$redhuisdf_par,
                       c("coef_par"="coef",
                         "pval_par"="pval")),
            by = join_by(gene_name, deletion))

ggplot(slice_sample(plotdf, n = 10000), 
       aes(x = cer,
           y = coef_cer)) +
  geom_point(aes(color = pval_cer < p_thresh)) +
  geom_vline(xintercept = c(-1, 1))

### smd versus coef in pois with interaction (most realistic model)

plotdf <- left_join(SMDdf, 
                    rename(TFdeldfs_poisint$redhuisdf_cer,
                           c("coef_cer"="coef",
                             "pval_cer"="pval")),
                    by = join_by(gene_name, deletion)) |> 
  left_join(y = rename(TFdeldfs_poisint$redhuisdf_par,
                       c("coef_par"="coef",
                         "pval_par"="pval")),
            by = join_by(gene_name, deletion))

# cer
ggplot(slice_sample(plotdf, n = 10000), 
       aes(x = cer,
           y = coef_cer)) +
  geom_point(aes(color = pval_cer < p_thresh)) +
  geom_vline(xintercept = c(-1, 1)) + ylim(c(-5, 5))

# par
ggplot(slice_sample(plotdf, n = 10000), 
       aes(x = par,
           y = coef_par)) +
  geom_point(aes(color = pval_par < p_thresh)) +
  geom_vline(xintercept = c(-1, 1)) + ylim(c(-5, 5))


#### barplots ####
coef_thresh <- 0.5
load("data_files/modules.RData")
plotTFBarplots <- function(.deldf_cer, .deldf_par, .m_genedf,
                           .useLiteralColors = TRUE) {
  common_TFs <- intersect(unique(.deldf_cer$deletion), unique(.deldf_par$deletion))
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
### parents
# redhuis alignment, pois with interaction
barplots_redhuis_parents <- plotTFBarplots(TFdeldfs_poisint$redhuisdf_cer, 
                                           TFdeldfs_poisint$redhuisdf_par, module_genedf25)
ggarrange(plotlist = barplots_redhuis_parents, nrow = 3, ncol = 1)
# barkai alignment, pois with interaction
barplots_barkai_parents <- plotTFBarplots(TFdeldfs_poisint$barkaidf_cer, 
                                           TFdeldfs_poisint$barkaidf_par, module_genedf25)
ggarrange(plotlist = barplots_barkai_parents, nrow = 3, ncol = 1)
# redhuis alignment, pois no interaction
barplots_redhuis_parents <- plotTFBarplots(TFdeldfs_pois$redhuisdf_cer, 
                                           TFdeldfs_pois$redhuisdf_par, module_genedf25)
ggarrange(plotlist = barplots_redhuis_parents, nrow = 3, ncol = 1)
# barkai alignment, pois no interaction
barplots_barkai_parents <- plotTFBarplots(TFdeldfs_pois$barkaidf_cer, 
                                          TFdeldfs_pois$barkaidf_par, module_genedf25)
ggarrange(plotlist = barplots_barkai_parents, nrow = 3, ncol = 1)
### hybrids
# redhuis alignment, pois with interaction
barplots_redhuis_hybrid <- plotTFBarplots(TFdeldfs_poisint$redhuisdf_hyc, 
                                           TFdeldfs_poisint$redhuisdf_hyp, module_genedf25)
ggarrange(plotlist = barplots_redhuis_hybrid, nrow = 3, ncol = 1)
# barkai alignment, pois with interaction
barplots_barkai_hybrid <- plotTFBarplots(TFdeldfs_poisint$barkaidf_hyc, 
                                          TFdeldfs_poisint$barkaidf_hyp, module_genedf25)
ggarrange(plotlist = barplots_barkai_hybrid, nrow = 3, ncol = 1)
# redhuis alignment, pois no interaction
barplots_redhuis_hybrid <- plotTFBarplots(TFdeldfs_pois$redhuisdf_hyc, 
                                           TFdeldfs_pois$redhuisdf_hyp, module_genedf25)
ggarrange(plotlist = barplots_redhuis_hybrid, nrow = 3, ncol = 1)
# barkai alignment, pois no interaction
barplots_barkai_hybrid <- plotTFBarplots(TFdeldfs_pois$barkaidf_hyc, 
                                          TFdeldfs_pois$barkaidf_hyp, module_genedf25)
ggarrange(plotlist = barplots_barkai_hybrid, nrow = 3, ncol = 1)

# conclusions: parents show same pattern regardless of which alignment/poisson formula is used
# hybrids have consistent DE genes in all modules across all TFs, which does not seem right

#### Do the hybrids really have DE genes in all modules for all TFs? ####

# parents
# pois with interaction, at least one sig gene x TF combo
# cer sig
random_gene_tf <- TFdeldfs_poisint$redhuisdf_cer |> 
  filter(abs(coef) > 1 & pval < p_thresh) |> 
  slice_sample(n = 1) |> select(gene_name, deletion)
genedf <- bind_cols(tibble(expr = redhuis_cts$cer[,random_gene_tf$gene_name]), redhuis_info$cer)
sigGenotypes <- TFdeldfs_negbin$redhuisdf_cer |> 
  filter(gene_name == random_gene_tf$gene_name &
           abs(coef) > 1 & pval < p_thresh) |> 
  select(deletion) |> pull()
ggplot(data = genedf, aes(x = genotype, y = expr)) +
  geom_point(aes(color = gsub("delete", "", genotype) %in% sigGenotypes)) +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")
# par sig
random_gene_tf <- TFdeldfs_poisint$redhuisdf_par |> 
  filter(abs(coef) > 1 & pval < p_thresh) |> 
  slice_sample(n = 1) |> select(gene_name, deletion)
TFdeldfs_poisint$redhuisdf_par |> 
  filter(gene_name == random_gene_tf$gene_name & 
           deletion == random_gene_tf$deletion)
genedf <- bind_cols(tibble(expr = redhuis_cts$par[,random_gene_tf$gene_name]), redhuis_info$par)
sigGenotypes <- TFdeldfs_poisint$redhuisdf_par |> 
  filter(gene_name == random_gene_tf$gene_name &
           abs(coef) > 1 & pval < p_thresh) |> 
  select(deletion) |> pull()
ggplot(data = genedf, aes(x = genotype, y = expr)) +
  geom_point(aes(color = gsub("delete", "", genotype) %in% sigGenotypes)) +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")

# TODO: why does YHL035C have a coef of 1 in AFT1delete, but there are
# plenty of other genotypes that appear to have the same effect 
# (ARG81/GLN3/PHD1 for example) that have small or negative effect sizes?
# TODO: check if notime model helps with this
genedf <- bind_cols(tibble(expr = redhuis_cts$par[,"YHL035C"]), redhuis_info$par)
ggplot(filter(genedf, genotype %in% c("WT", paste0(c("AFT1", "ARG81", "GLN3", "PHD1"), "delete"))), 
       aes(x = genotype, y = expr)) + geom_point(aes(color = time_point_str))
TFdeldfs_poisint$redhuisdf_par |> 
  filter(gene_name == "YHL035C" & 
           deletion == "AFT1")
TFdeldfs_poisint$redhuisdf_par |> 
  filter(gene_name == "YHL035C" & 
           deletion == "ARG81")
TFdeldfs_poisint$redhuisdf_par |> 
  filter(gene_name == "YHL035C" & 
           deletion == "GLN3")
TFdeldfs_poisint$redhuisdf_par |> 
  filter(gene_name == "YHL035C" & 
           deletion == "PHD1")

# hybrids
# pois with interaction, at least one sig gene x TF combo
# hyc sig
random_gene_tf <- TFdeldfs_poisint$redhuisdf_hyc |> 
  filter(abs(coef) > 1 & pval < p_thresh) |> 
  slice_sample(n = 1) |> select(gene_name, deletion)
# hyc plot
genedf <- bind_cols(tibble(expr = redhuis_cts_allele$cer[,random_gene_tf$gene_name]), redhuis_info_allele$cer)
sigGenotypes <- TFdeldfs_poisint$redhuisdf_hyc |> 
  filter(gene_name == random_gene_tf$gene_name &
           abs(coef) > 1 & pval < p_thresh) |> 
  select(deletion) |> pull()
ggplot(data = genedf, aes(x = genotype, y = expr)) +
  geom_point(aes(color = gsub("delete", "", genotype) %in% sigGenotypes)) +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")
# hyp plot
genedf <- bind_cols(tibble(expr = redhuis_cts_allele$par[,random_gene_tf$gene_name]), redhuis_info_allele$par)
sigGenotypes <- TFdeldfs_poisint$redhuisdf_hyp |> 
  filter(gene_name == random_gene_tf$gene_name &
           abs(coef) > 1 & pval < p_thresh) |> 
  select(deletion) |> pull()
ggplot(data = genedf, aes(x = genotype, y = expr)) +
  geom_point(aes(color = gsub("delete", "", genotype) %in% sigGenotypes)) +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")
# hyp sig
random_gene_tf <- TFdeldfs_poisint$redhuisdf_hyp |> 
  filter(abs(coef) > 1 & pval < p_thresh) |> 
  slice_sample(n = 1) |> select(gene_name, deletion)
# hyp plot
genedf <- bind_cols(tibble(expr = redhuis_cts_allele$par[,random_gene_tf$gene_name]), redhuis_info_allele$par)
sigGenotypes <- TFdeldfs_poisint$redhuisdf_hyp |> 
  filter(gene_name == random_gene_tf$gene_name &
           abs(coef) > 1 & pval < p_thresh) |> 
  select(deletion) |> pull()
ggplot(data = genedf, aes(x = genotype, y = expr)) +
  geom_point(aes(color = gsub("delete", "", genotype) %in% sigGenotypes)) +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")
# hyc plot
genedf <- bind_cols(tibble(expr = redhuis_cts_allele$cer[,random_gene_tf$gene_name]), redhuis_info_allele$cer)
sigGenotypes <- TFdeldfs_poisint$redhuisdf_hyc |> 
  filter(gene_name == random_gene_tf$gene_name &
           abs(coef) > 1 & pval < p_thresh) |> 
  select(deletion) |> pull()
ggplot(data = genedf, aes(x = genotype, y = expr)) +
  geom_point(aes(color = gsub("delete", "", genotype) %in% sigGenotypes)) +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")



#### Ensuring Hybrid LFC measures comparable to parents ####

# main goal: Get a measure of expression change that isn't sensitive to hybrids
#            having less variation between measurements than parents
#            (thus having larger expr change measurements than parents)
# I think this will be achievable using log fold change,
# which *I don't believe* is the same as exp(genotype coefficient) in the glm,
# which in the simplest case (expr ~ genotype) is just the ratio of
# means of expression (not the ratio of means of log2(expr)), 
# which would be comparable between expression levels
# I *believe* LFC is instead this ratio of log2(expr mutant)/log2(expr WT)
# but I need to verify that first

# TDH3 in GCR2delete is a good starting example. If this expr
# change isn't the same for hybrids and parents, nothing is
gene_idx <- "YGR192C"
genedf <- bind_cols(tibble(expr = redhuis_cts$cer[,gene_idx]), redhuis_info$cer) |> 
  filter(time_point_str == "16 h, low N")
genedf$genotype <- as.factor(genedf$genotype) |> relevel(ref = "WT")
sigGenotypes <- TFdeldfs_pois$redhuisdf_cer |> 
  filter(gene_name == gene_idx & abs(coef) > coef_thresh & pval < p_thresh) |> 
  select(deletion) |> pull()
ggplot(genedf, aes(x = genotype, y = expr)) + 
  geom_point(aes(color = gsub("delete", "", genotype) %in% sigGenotypes)) + 
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")
mod <- glm(expr ~ genotype, data = genedf, family = poisson(link = "log"))
mod$coefficients |> sort(decreasing = TRUE) # a cutoff of 1 would unfortunately still get GCR2, but that might be necessary
# now tf vs wt for gln3 specifically
genedf$p_hat <- mod$fitted.values
gcr2_effect_size <- mod$coefficients["genotypeGCR2delete"] |> as.numeric()
slope <- exp(gcr2_effect_size)
intercept <- mod$coefficients["(Intercept)"] |> as.numeric()
plotdf <- pivot_longer(genedf, cols = c(p_hat, expr))
ggplot(data = plotdf, aes(x = genotype, y = value)) +
  geom_point(aes(color = name, shape = as.factor(time_point_str))) +
  scale_x_discrete(breaks = NULL)
plotdf <- plotdf |> filter(genotype %in% c("WT", "GCR2delete")) |> 
  pivot_wider(id_cols = c("time_point_str", "name", "well_flask_ID"), 
              names_from = "genotype",
              values_from = "value")
ggplot(plotdf, aes(x = WT, y = GCR2delete)) + 
  geom_point(aes(color = name, shape = time_point_str)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "gold") +
  xlim(c(0, max(c(plotdf$WT, plotdf$GCR2delete), na.rm = TRUE))) +
  ylim(c(0, max(c(plotdf$WT, plotdf$GCR2delete), na.rm = TRUE)))
# genedf |> filter(genotype %in% c("WT", "GCR2delete")) |> 
#   select(expr, p_hat, genotype, well_flask_ID) |> View()
# mean WT:
genedf |> filter(genotype == "WT") |> 
  select(expr) |> pull() |> mean()
# mean GCR2 delete:
genedf |> filter(genotype == "GCR2delete") |> 
  select(expr) |> pull() |> mean()
# sd WT:
genedf |> filter(genotype == "WT") |> 
  select(expr) |> pull() |> sd()
# sd GCR2 delete:
genedf |> filter(genotype == "GCR2delete") |> 
  select(expr) |> pull() |> sd()

# effect size, log(mean(TFdel)/mean(WT)):
gcr2_effect_size
log(mean(genedf$expr[genedf$genotype == "GCR2delete"])/mean(genedf$expr[genedf$genotype == "WT"]))
# log2 effect size--- a value around 1 means expr doubled
log2(mean(genedf$expr[genedf$genotype == "GCR2delete"])/mean(genedf$expr[genedf$genotype == "WT"]))
log2(exp(gcr2_effect_size))
# "log-scale" fold change (lsfc, name pending): mean(log2(TFdel))/mean(log2(WT)):
mean(log2(genedf$expr[genedf$genotype == "GCR2delete"]))/mean(log2(genedf$expr[genedf$genotype == "WT"]))
# a decrease in expr is anything less than 1 now: lsfc(WT = 2^14, TFdel = 2^13) = 13/14 = 0.9ish

# repeat with par --- lower expressed and slightly smaller magnitude of expr change
genedf <- bind_cols(tibble(expr = redhuis_cts$par[,gene_idx]), redhuis_info$par) |> 
  filter(time_point_str == "16 h, low N")
genedf$genotype <- as.factor(genedf$genotype) |> relevel(ref = "WT")
sigGenotypes <- TFdeldfs_pois$redhuisdf_par |> 
  filter(gene_name == gene_idx & abs(coef) > coef_thresh & pval < p_thresh) |> 
  select(deletion) |> pull()
ggplot(genedf, aes(x = genotype, y = expr)) + 
  geom_point(aes(color = gsub("delete", "", genotype) %in% sigGenotypes)) + 
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")
mod <- glm(expr ~ genotype, data = genedf, family = poisson(link = "log"))
mod$coefficients |> sort(decreasing = TRUE) # a cutoff of 1 would unfortunately still get GCR2, but that might be necessary
# now tf vs wt for gln3 specifically
genedf$p_hat <- mod$fitted.values
gcr2_effect_size <- mod$coefficients["genotypeGCR2delete"] |> as.numeric()
slope <- exp(gcr2_effect_size)
intercept <- mod$coefficients["(Intercept)"] |> as.numeric()
plotdf <- pivot_longer(genedf, cols = c(p_hat, expr))
ggplot(data = plotdf, aes(x = genotype, y = value)) +
  geom_point(aes(color = name, shape = as.factor(time_point_str))) +
  scale_x_discrete(breaks = NULL)
plotdf <- plotdf |> filter(genotype %in% c("WT", "GCR2delete")) |> 
  pivot_wider(id_cols = c("time_point_str", "name", "well_flask_ID"), 
              names_from = "genotype",
              values_from = "value")
ggplot(plotdf, aes(x = WT, y = GCR2delete)) + 
  geom_point(aes(color = name, shape = time_point_str)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "gold") +
  xlim(c(0, max(c(plotdf$WT, plotdf$GCR2delete), na.rm = TRUE))) +
  ylim(c(0, max(c(plotdf$WT, plotdf$GCR2delete), na.rm = TRUE)))
# genedf |> filter(genotype %in% c("WT", "GCR2delete")) |> 
#   select(expr, p_hat, genotype, well_flask_ID) |> View()
# mean WT:
genedf |> filter(genotype == "WT") |> select(expr) |> pull() |> mean()
# mean GCR2delete:
genedf |> filter(genotype == "GCR2delete") |> select(expr) |> pull() |> mean()

# effect size/LFC, log(mean(TFdel)/mean(WT)):
gcr2_effect_size
log(mean(genedf$expr[genedf$genotype == "GCR2delete"])/mean(genedf$expr[genedf$genotype == "WT"]))
# "log-scale" fold change: mean(log2(TFdel))/mean(log2(WT)):
mean(log2(genedf$expr[genedf$genotype == "GCR2delete"]))/mean(log2(genedf$expr[genedf$genotype == "WT"]))
# this is roughly the same as cer

#### repeat with hybrid to see if LFC/effect sizes match parents
# hyc
genedf <- bind_cols(tibble(expr = redhuis_cts_allele$cer[,gene_idx]), redhuis_info_allele$cer) |> 
  filter(time_point_str == "16 h, low N")
genedf$genotype <- as.factor(genedf$genotype) |> relevel(ref = "WT")
sigGenotypes <- TFdeldfs_pois$redhuisdf_hyc |> 
  filter(gene_name == gene_idx & abs(coef) > coef_thresh & pval < p_thresh) |> 
  select(deletion) |> pull()
ggplot(genedf, aes(x = genotype, y = expr)) + 
  geom_point(aes(color = gsub("delete", "", genotype) %in% sigGenotypes)) + 
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")
mod <- glm(expr ~ genotype, data = genedf, family = poisson(link = "log"))
mod$coefficients |> sort(decreasing = TRUE) # a cutoff of 1 would unfortunately still get GCR2, but that might be necessary
# now tf vs wt for gln3 specifically
genedf$p_hat <- mod$fitted.values
gcr2_effect_size <- mod$coefficients["genotypeGCR2delete"] |> as.numeric()
slope <- exp(gcr2_effect_size)
intercept <- mod$coefficients["(Intercept)"] |> as.numeric()
plotdf <- pivot_longer(genedf, cols = c(p_hat, expr))
ggplot(data = plotdf, aes(x = genotype, y = value)) +
  geom_point(aes(color = name, shape = as.factor(time_point_str))) +
  scale_x_discrete(breaks = NULL)
plotdf <- plotdf |> filter(genotype %in% c("WT", "GCR2delete")) |> 
  pivot_wider(id_cols = c("time_point_str", "name", "well_flask_ID"), 
              names_from = "genotype",
              values_from = "value")
ggplot(plotdf, aes(x = WT, y = GCR2delete)) + 
  geom_point(aes(color = name, shape = time_point_str)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "gold") +
  xlim(c(0, max(c(plotdf$WT, plotdf$GCR2delete), na.rm = TRUE))) +
  ylim(c(0, max(c(plotdf$WT, plotdf$GCR2delete), na.rm = TRUE)))
# genedf |> filter(genotype %in% c("WT", "GCR2delete")) |> 
#   select(expr, p_hat, genotype, well_flask_ID) |> View()
# mean WT:
genedf |> filter(genotype == "WT") |> select(expr) |> pull() |> mean()
# mean GCR2delete:
genedf |> filter(genotype == "GCR2delete") |> select(expr) |> pull() |> mean()

# effect size, log(mean(TFdel)/mean(WT)):
gcr2_effect_size
log(mean(genedf$expr[genedf$genotype == "GCR2delete"])/mean(genedf$expr[genedf$genotype == "WT"]))
log2(mean(genedf$expr[genedf$genotype == "GCR2delete"])/mean(genedf$expr[genedf$genotype == "WT"]))
log2(exp(gcr2_effect_size))
# log fold change: mean(log2(TFdel))/mean(log2(WT)):
mean(log2(genedf$expr[genedf$genotype == "GCR2delete"]))/mean(log2(genedf$expr[genedf$genotype == "WT"]))

# repeat with hyp
genedf <- bind_cols(tibble(expr = redhuis_cts_allele$par[,gene_idx]), redhuis_info_allele$par) |> 
  filter(time_point_str == "16 h, low N")
genedf$genotype <- as.factor(genedf$genotype) |> relevel(ref = "WT")
sigGenotypes <- TFdeldfs_pois$redhuisdf_hyp |> 
  filter(gene_name == gene_idx & abs(coef) > coef_thresh & pval < p_thresh) |> 
  select(deletion) |> pull()
ggplot(genedf, aes(x = genotype, y = expr)) + 
  geom_point(aes(color = gsub("delete", "", genotype) %in% sigGenotypes)) + 
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")
mod <- glm(expr ~ genotype, data = genedf, family = poisson(link = "log"))
mod$coefficients |> sort(decreasing = TRUE) # a cutoff of 1 would unfortunately still get GCR2, but that might be necessary
# now tf vs wt for gln3 specifically
genedf$p_hat <- mod$fitted.values
gcr2_effect_size <- mod$coefficients["genotypeGCR2delete"] |> as.numeric()
slope <- exp(gcr2_effect_size)
intercept <- mod$coefficients["(Intercept)"] |> as.numeric()
plotdf <- pivot_longer(genedf, cols = c(p_hat, expr))
ggplot(data = plotdf, aes(x = genotype, y = value)) +
  geom_point(aes(color = name, shape = as.factor(time_point_str))) +
  scale_x_discrete(breaks = NULL)
plotdf <- plotdf |> filter(genotype %in% c("WT", "GCR2delete")) |> 
  pivot_wider(id_cols = c("time_point_str", "name", "well_flask_ID"), 
              names_from = "genotype",
              values_from = "value")
ggplot(plotdf, aes(x = WT, y = GCR2delete)) + 
  geom_point(aes(color = name, shape = time_point_str)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "gold") +
  xlim(c(0, max(c(plotdf$WT, plotdf$GCR2delete), na.rm = TRUE))) +
  ylim(c(0, max(c(plotdf$WT, plotdf$GCR2delete), na.rm = TRUE)))
# genedf |> filter(genotype %in% c("WT", "GCR2delete")) |> 
#   select(expr, p_hat, genotype, well_flask_ID) |> View()
# mean WT:
genedf |> filter(genotype == "WT") |> select(expr) |> pull() |> mean()
# mean GCR2delete:
genedf |> filter(genotype == "GCR2delete") |> select(expr) |> pull() |> mean()

# effect size, log(mean(TFdel)/mean(WT)):
gcr2_effect_size
log(mean(genedf$expr[genedf$genotype == "GCR2delete"])/mean(genedf$expr[genedf$genotype == "WT"]))
log2(mean(genedf$expr[genedf$genotype == "GCR2delete"])/mean(genedf$expr[genedf$genotype == "WT"]))
log2(exp(gcr2_effect_size))
# log fold change: mean(log2(TFdel))/mean(log2(WT)):
mean(log2(genedf$expr[genedf$genotype == "GCR2delete"]))/mean(log2(genedf$expr[genedf$genotype == "WT"]))


# secondary goal: find an appropriate threshold for our chosen
# measure of expr change
# YNL065W is a good example to explore.
# Much stronger effect in GLN3 in hyp (somewhat in hyc)
# than other genotypes, but most genotypes are called sig in glms
gene_idx <- "YNL065W"
# hyp has one main regulator: GLN3
# Then a couple minor ones: AFT1, HAP1, MBP1, ROX1, SOK2
# Then a few dubious ones we wouldn't want to count: GCR2 (maybe real but v noisy), then basically all the other TFs except DAL80, PHD1, and TEC1 are called significant but they shouldn't be
# TODO: figure out what effect size threshold works for this gene in hyp
# then see if it's the same threshold that works for parents (just par for starters)

# Note: YNL065W no longer seems like a good example --- not 
# obviously more DE in GLN3delete (it's really just lowly expressed and noisy)
