sapply(c("tidyr", "dplyr", "readr", "purrr"), require, character.only=TRUE)
library("MASS", include.only = c("glm.nb", "theta.ml"))
library("lme4", include.only = "glmer.nb")
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Barkai_data_analysis/")

# Using single-gene negative binomial GLMs to identify differential expression
# Part I: Identifying genes DE in specific environments (testing each species/experiment separately, 12 datasets total (4 experiments x cer/par/hyb))
# Part II: Identifying genes DE between parental species and hybrid alleles (cis and trans-diverging genes)

# load barkai data as DESeq2 datasets
load("data_files/Cleaned_Barkai_Data_Env_Specific.RData")
spes <- es
spsi <- si
GeneNames <- rownames(spes[[1]])
ExperimentNames <- c("TFdelxLowN", "CC", "HAP4", "LowPi")

# combining two species' counts for each experiment into the same dataset
# filtering down to conditions/replicates (ex. LowN_TP1_rep2) where both species have data
spcts <- vector(mode = "list", length = 4)
names(spcts) <- ExperimentNames
spinfo <- vector(mode = "list", length = 4)
names(spinfo) <- ExperimentNames
for(e in ExperimentNames) {
  cer_cts <- spes[[paste(e, "cer", sep = "_")]]
  par_cts <- spes[[paste(e, "par", sep = "_")]]
  cer_info <- spsi[[paste(e, "cer", sep = "_")]]
  par_info <- spsi[[paste(e, "par", sep = "_")]]
  cer_info$condition_with_rep <- paste(cer_info$condition, cer_info$well_flask_ID, sep = "_") # TODO: this should prob go in cleaning script
  par_info$condition_with_rep <- paste(par_info$condition, par_info$well_flask_ID, sep = "_")
  if(nrow(cer_info) >= nrow(par_info)) {
    infos <- list(larger = cer_info, smaller = par_info)
  }
  if(nrow(cer_info) < nrow(par_info)) {
    infos <- list(larger = par_info, smaller = cer_info)
  }
  larger_samps <- infos$larger$condition_with_rep
  smaller_samps <- infos$smaller$condition_with_rep
  common_samples <- smaller_samps[which(smaller_samps %in% larger_samps)]
  cer_keep <- which(cer_info$condition_with_rep %in% common_samples)
  par_keep <- which(par_info$condition_with_rep %in% common_samples)
  # filtering
  cer_cts <- cer_cts[,cer_keep]
  cer_info <- cer_info[cer_keep,]
  par_cts <- par_cts[,par_keep]
  par_info <- par_info[par_keep,]
  # combining
  spcts[[e]] <- cbind(cer_cts, par_cts)
  spinfo[[e]] <- bind_rows(cer_info, par_info)
  cat("number of times each condition/rep is represented in", e, ":", 
      names(table(table(spinfo[[e]]$condition_with_rep))), "\n") # we want this to be two each time, cer and par
}
# factorize genotype, experiment, and time_point_str (no specific reference level for the last two)
spinfo <- lapply(spinfo, mutate, genotype = as.factor(genotype) %>% relevel(ref = "WT"))
# releveling time_point_str is a tad more tricky b/c the reference is named differently in each experiment (and isn't 0 in LowPi fml)
spinfo <- lapply(spinfo, mutate, time_point_str = as.factor(time_point_str) %>% relevel(ref = time_point_str[which.min(parse_number(time_point_str))]))
spinfo <- lapply(spinfo, mutate, allele = as.factor(allele) %>% relevel(ref = "par")) # we always do cer/par for expression ratios, so we want a big estimate to also mean cer is more strongly expressed

# Single gene divergence models

# okay hear me out... all our samples are condition-matched so why not just do a linear model regressing cer vs par?
# sure the variance is expected to increase with the mean, but that'll only lower our discovery rate not raise it
gene_idx <- sample(GeneNames, 1)
experiment_idx <- sample(ExperimentNames, 1)
gdf <- bind_cols(expr = spcts[[experiment_idx]][gene_idx,],
                 spinfo[[experiment_idx]])
spdf <- gdf %>% select(allele, condition, well_flask_ID, expr, genotype) %>% 
  pivot_longer(cols = c(expr)) %>% 
  mutate(non_unique_sample_name = paste(condition, well_flask_ID)) %>% 
  pivot_wider(id_cols = c(non_unique_sample_name, name, genotype), names_from = allele, values_from = value)
ggplot(spdf, aes(x = cer, y = par)) + geom_point() + geom_abline(color = "gold") + 
  xlim(c(0, max(select(spdf, cer, par)))) + ylim(c(0, max(select(spdf, cer, par)))) + 
  theme_classic() + ggtitle("Standard LM")
mod <- lm(cer ~ par, data = spdf)
summary(mod) # ohhhh right. If cer predicts par, that means they're not divergent. Interpretability is tricky here

# Exploration: what is theta normally estimated as?
thetas <- vector(mode = "list")
for (i in 1:1000) {
  cat(i, "/1000\n")
  gene_idx <- sample(GeneNames, 1)
  experiment_idx <- sample(ExperimentNames, 1)
  gdf <- bind_cols(expr = spcts[[experiment_idx]][gene_idx,],
                   spinfo[[experiment_idx]])
  form <- if_else(table(gdf$genotype) %>% length() > 2, 
                  true = "expr ~ time_point_str + genotype + time_point_str:genotype + allele + allele:genotype", 
                  false = "expr ~ time_point_str + allele")
  pois <- glm(form, family = poisson(link = "log"), data = gdf)
  thetaml <- MASS::theta.ml(pois, weights = pois$weights, trace = TRUE)
  thetas[[paste(gene_idx, experiment_idx, sep = "_")]] <- thetaml[1]
}
hist(unlist(thetas), breaks = 50)
median(unlist(thetas))
max(unlist(thetas))
quantile(unlist(thetas)) # 7 seems like a fine initial theta

# what's the theta of a highly expressed, not very variable gene like TDH3 in hybrid?
gene_idx <- "YGR192C"
experiment_idx <- sample(ExperimentNames, 1)
gdf <- bind_cols(expr = alcts[[experiment_idx]][gene_idx,],
                 alinfo[[experiment_idx]])
form <- if_else(table(gdf$genotype) %>% length() > 2, 
                true = "expr ~ time_point_str + genotype + time_point_str:genotype + allele + allele:genotype", 
                false = "expr ~ time_point_str + allele")
pois <- glm(form, family = poisson(link = "log"), data = gdf)
thetaml <- MASS::theta.ml(pois, weights = pois$weights, trace = TRUE)
thetaml[1] # not that high (parents is even lower)
nb <- glm.nb(form, gdf, link = log)
nb$theta
nb <- glm.nb(form, gdf, link = log, init.theta = 7)
nb$theta

# what kind of gene has a very high theta?
hightheta <- strsplit(names(thetas)[which.max(unlist(thetas))], split = "_") %>% unlist()
gene_idx <- hightheta[1]
experiment_idx <- hightheta[2]
gdf <- bind_cols(expr = spcts[[experiment_idx]][gene_idx,],
                 spinfo[[experiment_idx]]) 
plotdf <- gdf %>% mutate(non_unique_sample_name = paste(condition, well_flask_ID, sep = "_")) %>% 
  pivot_wider(id_cols = c(non_unique_sample_name, genotype, time_point_str), names_from = allele, values_from = expr)
ggplot(plotdf, aes(x = cer, y = par)) + geom_point(aes(color = time_point_str)) + theme(legend.position = "none") + geom_abline(color = "gold")
# will the estimate change if we set init.theta to 7?
form <- if_else(table(gdf$genotype) %>% length() > 2, 
                true = "expr ~ time_point_str + genotype + time_point_str:genotype + allele + allele:genotype", 
                false = "expr ~ time_point_str + allele")
pois <- glm(form, gdf, family = poisson(link = "log"))
MASS::theta.ml(pois, weights = pois$weights, trace = TRUE)
noinitnb <- glm.nb(form, gdf, link = "log")
summary(noinitnb)$coefficients["allelecer",]
noinitnb$theta
yesinitnb <- glm.nb(form, gdf, link = "log", init.theta = 7)
summary(yesinitnb)$coefficients["allelecer",] # nope
yesinitnb$theta # literally makes no difference

# will setting init.theta affect genes that would've failed?
gene_idx <- "YDR243C"
experiment_idx <- "TFdelxLowN"
gdf <- bind_cols(expr = spcts[[experiment_idx]][gene_idx,],
                 spinfo[[experiment_idx]]) 
plotdf <- gdf %>% mutate(non_unique_sample_name = paste(condition, well_flask_ID, sep = "_")) %>% 
  pivot_wider(id_cols = c(non_unique_sample_name, genotype, time_point_str), names_from = allele, values_from = expr)
ggplot(plotdf, aes(x = cer, y = par)) + geom_point(aes(color = time_point_str)) + theme(legend.position = "none") + geom_abline(color = "gold")
form <- if_else(table(gdf$genotype) %>% length() > 2, 
                true = "expr ~ time_point_str + genotype + time_point_str:genotype + allele + allele:genotype", 
                false = "expr ~ time_point_str + allele")
pois <- glm(form, gdf, family = poisson(link = "log"))
MASS::theta.ml(pois, weights = pois$weights, trace = TRUE)
noinitnb <- glm.nb(form, gdf, link = "log")
summary(noinitnb)$coefficients["allelecer",]
noinitnb$theta
yesinitnb <- glm.nb(form, gdf, link = "log", init.theta = 7)
summary(yesinitnb)$coefficients["allelecer",]
yesinitnb$theta # Weird. I could've sworn there was a TFdelxLowN gene that super didn't look divergent but glm.nb gave a high theta (like 1000) that was reduced to like 2.5 when you set init.theta = 5
# bottom line: setting init theta to 7 doesn't have much effect, and it might prevent some abbarant thetas so might as well

# TODO: Issue: LowPi seems to have larger effect sizes for no good reason
# example to check: YDR063W in parents
test_df_lowPi <- bind_cols(expr = spcts[["LowPi"]]["YDR063W",],
                           spinfo[["LowPi"]])
test_df_TFdel <- bind_cols(expr = spcts[["TFdelxLowN"]]["YDR063W",],
                           spinfo[["TFdelxLowN"]])
test_df_HAP4 <- bind_cols(expr = spcts[["HAP4"]]["YDR063W",],
                          spinfo[["HAP4"]])
test_df_lowPi %>% pivot_wider(id_cols = c(condition_with_rep, time_point_num), values_from = expr, names_from = allele) %>% 
  ggplot(aes(x = cer, y = par)) + geom_point(aes(color = time_point_num)) + geom_abline(color = "gold")
test_df_TFdel %>% pivot_wider(id_cols = c(condition_with_rep, time_point_str), values_from = expr, names_from = allele) %>% 
  ggplot(aes(x = cer, y = par)) + geom_point(aes(color = time_point_str)) + geom_abline(color = "gold")
test_df_HAP4 %>% pivot_wider(id_cols = c(condition_with_rep, genotype), values_from = expr, names_from = allele) %>% 
  ggplot(aes(x = cer, y = par)) + geom_point(aes(color = genotype)) + geom_abline(color = "gold")

mod_LowPi <- glm.nb(expr ~ time_point_num + allele, data = test_df_lowPi)
summary(mod_LowPi)$coefficient["allelecer",] # time_point_str leads to larger effect size
mod_TFdel <- glm.nb(expr ~ genotype + time_point_str + allele + genotype:time_point_str + genotype:allele, data = test_df_TFdel)
summary(mod_TFdel)$coefficient["allelecer",]
mod_HAP4 <- glm.nb(expr ~ genotype + time_point_num + genotype:time_point_num + allele, data = test_df_TFdel)
summary(mod_HAP4)$coefficient["allelecer",] # doesn't do the same for HAP4 for some reason
# Solution: we use time_point_num for HAP4, LowPi, and CC

# Issue: low counts obscure divergence of some genes
# ex) YOL052C-A in TFdelxLowN
# possible solution: subtract a basal count number (say 50) from every sample to reduce noise at the low end driving divergence
# Actual solution I used: filtered to make sure the gene in that experiment had expression counts where the 75% percentile was greater than 100
gene_idx <- "YOL052C-A"
experiment_idx <- 1
# before subtracting basal count
gdf <- bind_cols(expr = spcts[[experiment_idx]][gene_idx,],
                 spinfo[[experiment_idx]]) 
spmod <- glm.nb(expr ~ time_point_str * genotype * allele, data = gdf)
sprow <- summary(spmod)$coefficients["allelecer",]
slope <- 1/exp(sprow[1])
intercept <- summary(spmod)$coefficients["(Intercept)",, drop = FALSE][1]
plotdf <- gdf %>% mutate(non_unique_sample_name = paste(condition, well_flask_ID, sep = "_")) %>% 
  pivot_wider(id_cols = c(non_unique_sample_name, genotype, time_point_str), names_from = allele, values_from = expr)
ggplot(plotdf, aes(x = cer, y = par)) + geom_point(aes(color = genotype)) + geom_abline(color = "gold", slope = slope, intercept = intercept) + geom_abline(color = "navy") + 
  xlim(c(0, max(select(plotdf, cer, par)))) + ylim(c(0, max(select(plotdf, cer, par)))) + 
  theme(legend.position = "none") + ggtitle("Standard GLM")
# after subtracting basal count
gdf$expr <- gdf$expr - 50
gdf$expr[gdf$expr < 0] <- 0
condition_matched <- gdf %>% mutate(non_unique_sample_name = paste(condition, well_flask_ID, sep = "_")) %>% 
  pivot_wider(id_cols = c(non_unique_sample_name, genotype, time_point_str), names_from = allele, values_from = expr) %>% 
  filter(cer != 0 & par != 0)
gdf <- condition_matched %>% pivot_longer(cols = c(cer, par), names_to = "allele", values_to = "expr")
gdf$allele <- as.factor(gdf$allele) %>% relevel(ref = "par")
spmod <- glm.nb(expr ~ genotype * allele, data = gdf)
sprow <- summary(spmod)$coefficients["allelecer",]
slope <- 1/exp(sprow[1])
intercept <- summary(spmod)$coefficients["(Intercept)",, drop = FALSE][1]
plotdf <- gdf %>%  
  pivot_wider(id_cols = c(non_unique_sample_name, genotype, time_point_str), names_from = allele, values_from = expr)
ggplot(plotdf, aes(x = cer, y = par)) + geom_point(aes(color = time_point_str)) + geom_abline(color = "gold", slope = slope, intercept = intercept) + geom_abline(color = "navy") + 
  xlim(c(0, max(select(plotdf, cer, par)))) + ylim(c(0, max(select(plotdf, cer, par)))) + 
  theme(legend.position = "none") + ggtitle("Standard GLM")

# glm.nb
gene_idx <- sample(GeneNames, 1)
experiment_idx <- sample(ExperimentNames, 1)
gdf <- bind_cols(expr = spcts[[experiment_idx]][gene_idx,],
                 spinfo[[experiment_idx]])
form <- if_else(table(gdf$genotype) %>% length() > 2, 
                true = "expr ~ time_point_str + genotype + time_point_str:genotype + allele + allele:genotype", 
                false = "expr ~ time_point_str + allele")
spmod <- glm.nb(formula(form), data = gdf)
# if the model fails, it's likely because theta is going to infinity, you can check with this:
# initpois <- glm(form, gdf, family = poisson(link = "log"))
# MASS::theta.ml(initpois, weights = initpois$weights, trace = TRUE)
plotdf <- gdf %>% mutate(fittedvals = spmod$fitted.values)%>% select(allele, condition, well_flask_ID, expr, fittedvals, genotype) %>% 
  pivot_longer(cols = c(expr, fittedvals)) %>% 
  mutate(non_unique_sample_name = paste(condition, well_flask_ID)) %>% 
  pivot_wider(id_cols = c(non_unique_sample_name, name, genotype), names_from = allele, values_from = value)
sprow <- summary(spmod)$coefficients["allelecer",]
slope <- 1/exp(sprow[1])
intercept <- summary(spmod)$coefficients["(Intercept)",, drop = FALSE][1]
ggplot(plotdf, aes(x = cer, y = par)) + geom_point(aes(color = name)) + geom_abline(color = "gold", slope = slope, intercept = intercept) + geom_abline(color = "navy") + 
  xlim(c(0, max(select(plotdf, cer, par)))) + ylim(c(0, max(select(plotdf, cer, par)))) + 
  theme_classic() + ggtitle("Standard GLM")
glm_results <- summary(spmod)$coefficients["allelecer",]
# allele effect (log-fold change of switching from par to cer allele) aka inverse-exponent of the slope of fitted vals line:
allele_effect <- glm_results[1] %>% as.numeric()
allele_effect
# slope:
1/exp(allele_effect)
# calculating rise over run to prove this is the slope:
# rise (par2 - par1) <- two arbitrary par fitted values (par2 > par1), largest and smallest for simplicity
par1 <- sort(spmod$fitted.values[gdf$allele == "par"])[1] %>% as.numeric()
par2 <- sort(spmod$fitted.values[gdf$allele == "par"])[sum(gdf$allele == "par")] %>% as.numeric()
rise <- par2 - par1
# run (corresponding cer2 - cer1) <- the only difference between condition-matched cer and par measurements is the allele effect size
cer1 <- par1*exp(allele_effect)
cer2 <- par2*exp(allele_effect)
run <- cer2 - cer1
# slope:
rise/run

# repeat with model that accounts for pseudoreplication (reps names are only shared by samples that came from the same exact well or flask (ex) rep2_SOK2delete, rep24_WT)
form <- ifelse(experiment_idx %in% c("TFdelxLowN", "HAP4"), 
               yes = "expr ~ time_point_str * genotype * allele + (allele|well_flask_ID)", 
               no = "expr ~ time_point_str * allele + (allele|well_flask_ID)")
gdf$well_flask_ID <- paste(gdf$well_flask_ID, gdf$genotype, gdf$allele, sep = "_")
spmod_mixed <- glmer.nb(formula(form), data = gdf) # This takes HOURS, especially for TFdelxLowN
plotdf <- gdf %>% mutate(fittedvals = fitted(spmod_mixed)) %>% select(allele, condition, well_flask_ID, expr, fittedvals) %>% 
  pivot_longer(cols = c(expr, fittedvals)) %>% 
  mutate(non_unique_sample_name = paste(condition, well_flask_ID)) %>% 
  pivot_wider(id_cols = c(non_unique_sample_name, name), names_from = allele, values_from = value)
ggplot(plotdf, aes(x = cer, y = par)) + geom_point(aes(color = name)) + geom_abline(color = "gold") + 
  xlim(c(0, max(select(plotdf, cer, par)))) + ylim(c(0, max(select(plotdf, cer, par)))) + 
  theme_classic() + ggtitle("Mixed Effects GLM")
glmer_results <- summary(spmod_mixed)$coefficients["allelecer",]
# I expected glmer to tank the p-value but didn't realize it might also lower effect sizes:
rbind(glm_results, glmer_results)

# Where I left mixed models for now: I haven't observed much of a replicate effect that
# would warrant use of a mixed effect model given they're suuuuper computationally expensive,
# and would require the cluster. We're talking overnight to run one glmer.nb versus seconds to run one glm.nb
# So until I know I need one, I won't use one
# TODO: to quantify whether replicates have an effect or not, I could look into doing a 
# QC ANOVA to see if the variance explained by replicate is anywhere as great as
# variance explained by allele (cer or par)

# likelihood ratio test
form <- if_else(experiment_idx %in% c("TFdelxLowN", "HAP4"), 
                true = "expr ~ time_point_str + genotype + time_point_str:genotype", 
                false = "expr ~ time_point_str")
spmod_reduced <- glm.nb(formula(form), data = gdf)
lrt <- lrtest(spmod, spmod_reduced)
lrt$`Pr(>Chisq)`[2]
# we determined that the pvalues obtained from the LRT are less stringent than 
# the standard Wald test, plus a low LRT pvalue was notably less predictive of a strong effect size

# Filtering for specific desired datasets (I think this was to get rid of plain TFdel datasets and replicate datasets, which I ended up doing in env-specific cleaning script)
# es0 <- es
# si0 <- si
# bind_cols(names(es), c(1:length(es))) %>% print(n = length(es))
# keep <- c(3:6, 9:12, 15:18) # hard coding so I don't go insane
# names(es)[keep]
# es <- es[keep]
# si <- si[keep]
# cbind(names(es), names(si))

# setting some values we'll use a lot
GeneNames <- rownames(es[[1]])
DatasetNames <- names(es)
nGenes <- nrow(es[[1]])
nDatasets <- length(es)
pthresh <- 0.05/nGenes

# Skipping b/c we're just treating CC like a stress anyway
# removing CC samples in HU, as HU is causing changes in gene expression but not in a linear way and we're more interested in using CC as a control here
# for (i in c(1:nDatasets)) {
#   if (si[[i]]$experiment[1] == "CC") {
#     keep <- !grepl("HU", si[[i]]$time_point_str)
#     es[[i]] <- es[[i]][,keep]
#     si[[i]] <- si[[i]][keep,]
#   }
# }

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
si <- lapply(si, factorizeGenotypeAndTimepoint)

# normalizing counts
# @input: count matrix (genes are rows, columns are samples like God intended)
# @output: a count matrix normalzied for library size---integer counts in counts-per-million
countsPerMillion <- function(.cts) {
  librarySizes <- colSums(.cts)
  output <- apply(.cts, 1, function(x) {
    normalized <- (x/librarySizes)*1e6
    return(round(normalized))
  })
  return(t(output)) # For some unhinged reason, a vector output of apply across ROWS forms the COLUMNS of a new matrix
}
# tests for countsPerMillion
test_cts <- es[[sample(c(1:length(es)),1)]] # re-run to try different count matrix
test_rowIdx <- sample(c(1:nGenes), 1)
test_colIdx <- sample(c(1:ncol(test_cts)), 1)
test_count <- test_cts[test_rowIdx, test_colIdx]
((test_count/colSums(test_cts)[test_colIdx])*1e6) %>% round() # what it should be
test_cpm <- countsPerMillion(test_cts)
test_cpm[test_rowIdx, test_colIdx] # what it is using our function

# normalizing all count matrices we'll use
es <- lapply(es, countsPerMillion)

##### Exploring GLMs with single gene #####
# finding a candidate gene with enough variation
cts <- es$TFdelxLowN_cer # change these to try different dataset
info <- si$TFdelxLowN_cer
geneVars <- apply(cts, 1, var, na.rm = TRUE)
geneMeans <- apply(cts, 1, mean, na.rm = TRUE)
varying_idxs <- which(geneVars/geneMeans > 100)
gene_idx <- sample(varying_idxs, 1) # YHR138C is a fun one
genedf <- bind_cols(tibble(expr = cts[gene_idx,]), info)

# example of how you check for full rank of design matrix (Needs to be full rank to guarantee each variable value can be measured independently of values in other variables (example: If there's no HAP4delete at LowN TP3, the model can't separate these two variable values)
mm <- stats::model.matrix(expr ~ time_point_num, data = genedf)
q <- qr(mm)
q$rank == ncol(mm) # TRUE = Full Rank

# trying with glm.nb
library(MASS, include.only = c("glm.nb", "theta.ml")) # MASS has a select function that masks dplyr's, so let's just get glm.nb---that's all I use
modnb <- glm.nb(expr ~ time_point_num, data = genedf)
summary(modnb)
genedf$p_hat <- modnb$fitted.values
plotdf <- pivot_longer(genedf, cols = c(expr, p_hat))
ggplot(data = plotdf, aes(x = time_point_num, y = value)) + geom_point(aes(color = name, shape = well_flask_ID))

# linear predictors
beta1 <- modnb$coefficients["time_point_num"] %>% as.numeric()
beta0 <- modnb$coefficients["(Intercept)"] %>% as.numeric()
tibble(maybe = beta1*genedf$time_point_num + beta0,
       definitely = modnb$linear.predictors) # yup that's what linear preditors are

# fitted values
tibble(maybe = exp(modnb$linear.predictors),
       definitely = modnb$fitted.values) # and that's what fitted values are!

# What's a large enough effect size to be important? I.e. what's a large enough beta1?
# Let's see what it is for a gene that's definitely not change expr by timepoint
gene_idx <- sample(c(1:nrow(es[[1]])), 1)
genedf <- bind_cols(tibble(expr = es$LowPi_cer[gene_idx,]),
                    si$LowPi_cer)
ggplot(genedf, aes(x = time_point_num, y = expr)) + geom_point()
modnb <- glm.nb(expr ~ time_point_num, data = genedf)
summary(modnb)
genedf$p_hat <- modnb$fitted.values %>% as.numeric()
plotdf <- pivot_longer(genedf, cols = c(expr, p_hat))
ggplot(data = plotdf, aes(x = time_point_num, y = value)) + geom_point(aes(color = name))
summary(modnb)$coefficients

# seems like anything less than magnitude 0.001 is a good effect size cutoff
MinEffectSizeContinuous <- 0.001

# TFdel example
gene_idx <- sample(varying_idxs, 1) # YLR307W 3199, a sporulation gene, has crazy high expression specifically in SUM1delete, a sporulation repressor, so that's cool
genedf <- bind_cols(tibble(expr = es$TFdelxLowN_cer[gene_idx,]), si$TFdelxLowN_cer)
ggplot(genedf, aes(x = genotype, y = expr)) + geom_point() + scale_x_discrete(breaks = NULL)
modnb <- glm.nb(expr ~ time_point_str + genotype + time_point_str:genotype, data = genedf, link = log)
summary(modnb) 
genedf$p_hat <- modnb$fitted.values
plotdf <- pivot_longer(genedf, cols = c(p_hat, expr))
ggplot(data = plotdf, aes(x = genotype, y = value)) + 
  geom_point(aes(color = name, shape = as.factor(time_point_str))) +
  scale_x_discrete(breaks = NULL)
coeffs <- coef(modnb)[grepl("^genotype", names(coef(modnb)))]
pvals <- summary(modnb)$coefficients[,4][grepl("^genotype", names(coef(modnb)))]
sigGenotypes <- names(coeffs[abs(coeffs) > 1 & pvals < pthresh]) %>% gsub(pattern = "genotype", replacement = "")
genedf$sig <- genedf$genotype %in% sigGenotypes
genedf$sig[genedf$genotype == "WT"] <- NA
ggplot(data = genedf, aes(x = genotype, y = expr)) + 
  geom_point(aes(color = sig)) +
  scale_x_discrete(breaks = NULL)
library(marginaleffects)
# This only works when you find a significant gene:
# avg_comparisons(modnb, variables = list(genotype = c("WT", sigGenotypes)), newdata = datagrid(model = modnb, genotype = genedf$genotype, timepoint = genedf$timepoint))
MinEffectSizeCatagorical <- 1

##### Data Exploration: Estimate Proportion of False Positives in my data #####
# Goal: make sure most genes aren't being called as DE
# Note: estimating FDR this way takes about 1hr per species, 
# and we've already got an estimated FDR
# so only run this block if anything about the model or dataset changes

# A random two-group partition of WT TP1 LowN in each species/hybrid shouldn't have any true DE genes
# So the percent of genome with DE in these null sets will be our FDR
fitSingleGeneNullModel <- function(.df) {
  if (var(.df$expr) == 0) {
    return(1)
  }
  mod <- glm.nb(expr ~ group, data = .df, link = log)
  grouprow <- summary(mod)$coefficients[2,]
  return(as.numeric(grouprow["Pr(>|z|)"]))
}

estimateFDR <- function(.cts, .groups) {
  pvalues <- apply(.cts, 1, function(x) {
    genedf <- tibble(expr = x,
                     group = .groups)
    return(fitSingleGeneNullModel(genedf))
  })
  return(sum(pvalues < 0.05)/length(pvalues))
}

applyFDREstimationToCounts <- function(.cts, niter) {
  fdrs <- map(c(1:niter), function(i) {
    cat(i, "\n")
    groups_null <- sample(c("A","B"), size = ncol(.cts), replace = TRUE, prob = c(0.5, 0.5))
    fdr <- tryCatch({
      estimateFDR(.cts, groups_null)
    }, error = function(e) {return(NA)})
    if (!is.na(fdr)) {
      if (fdr > 0.1) {
        cat("fdr exceeded 10% for sample groupings", groups_null, "\n")
      }
    }
    return(fdr)
  })
  return(fdrs)
}

counts_null_cer <- es$LowN_cer[, si$LowN_cer$genotype == "WT" & si$LowN_cer$time_point_num == 0]
counts_null_par <- es$LowN_par[, si$LowN_par$genotype == "WT" & si$LowN_par$time_point_num == 0]
counts_null_hyb <- es$LowN_hyb[, si$LowN_hyb$genotype == "WT" & si$LowN_hyb$time_point_num == 0]
counts_null <- list(counts_null_cer, counts_null_par, counts_null_hyb)
fdrs <- lapply(counts_null, applyFDREstimationToCounts, niter = 100)

# They're all .06-.07. DESeq defaults to assume that 10% of tests will be false positives, so this is in line with that
median(unlist(fdrs[[1]]))
median(unlist(fdrs[[2]]))
median(unlist(fdrs[[3]]))
median(unlist(fdrs))

mean(unlist(fdrs[[1]]))
mean(unlist(fdrs[[2]]))
mean(unlist(fdrs[[3]]))
mean(unlist(fdrs))

max(unlist(fdrs[[1]]))
max(unlist(fdrs[[2]]))
max(unlist(fdrs[[3]]))
max(unlist(fdrs)) # and none are 76% FDR anymore thankfully, now that I normalized

##### Functions for fitting env-specific models #####


# Helper function for fitTPGenotypeModelAllGenes
# @input: gene data frame including columns named time_point_num, well_flask_ID, and (if .isGxT = TRUE) genotype
# @output: model generated by glm.nb
fitSingleGeneTPGenotypeModel <- function(.df, .gene_name) {
  experiment <- unique(.df$experiment)
  output <- NA
  if (var(.df$expr) == 0) {
    return(output)
  }
  if (experiment == "LowN") { # TFdelxLowN, multiple genotypes and catagorical timepoint (interested in genotype, controlling for timepoint)
    output <- tryCatch({
      mod <- glm.nb(expr ~ time_point_str + genotype + time_point_str:genotype, .df, link = log)
      return(mod)
    }, error = function(e) {
      return(NA)
    }, warning = function(w) {
      return(NA)
    })
  }
  if (experiment == "HAP4") { # HAP4, two genotypes and continuous timepoint (interested in timepoint, controlling for genotype)
    output <- tryCatch({
      mod <- glm.nb(expr ~ genotype + time_point_num + genotype:time_point_num, .df, link = log)
      return(mod)
    }, error = function(e) {
      return(NA)
    }, warning = function(w) {
      return(NA)
    })
  }
  if (experiment %in% c("CC", "LowPi")) { # single genotype and continuous timepoint
    output <- tryCatch({
      mod <- glm.nb(expr ~ time_point_num, .df, link = log)
      return(mod)
    }, error = function(e) {
      return(NA)
    }, warning = function(w) {
      return(NA)
    })
  }
  return(output)
}
# tests for fitSingleGeneTPGenotypeModel
gene_idx <- sample(c(1:nGenes), 1)
dataset_idx <- sample(c(1:length(es)), 1)
test <- bind_cols(tibble(expr = es[[dataset_idx]][gene_idx,]), si[[dataset_idx]])
testmod <- fitSingleGeneTPGenotypeModel(test, GeneNames[gene_idx])
summary(testmod)$coefficients
# this gene failed mapping but seems fine for the individual gene:
gene_idx <- 2048
dataset_idx <- 9
test <- bind_cols(tibble(expr = es[[dataset_idx]][gene_idx,]), si[[dataset_idx]])
testmod <- fitSingleGeneTPGenotypeModel(test, GeneNames[gene_idx])
summary(testmod)$coefficients
# ultra low expressed gene in TFdel, shouldn't have any DE instances:
gene_idx <- 1345
dataset_idx <- 1
test <- bind_cols(tibble(expr = es[[dataset_idx]][gene_idx,]), si[[dataset_idx]])
testmod <- fitSingleGeneTPGenotypeModel(test, GeneNames[gene_idx])
pvals <- summary(testmod)$coefficients[,4] %>% as.numeric() # good nothing's significant here
pvals
sum(pvals < pthresh)

################################# CSCAR issue ###############################
# the warnings/errors from what I can tell are more common in expression sets with a lot of zeros:
test_genedf <- bind_cols(tibble(expr = es$TFdelxLowN_par[495,]), si$TFdelxLowN_par)
test_mod <- glm.nb(expr ~ time_point_num + genotype, data = test_genedf, link = log) # this should run fine
test_genedf$expr[sample(c(1:nrow(test_genedf)), floor(nrow(test_genedf)/3))] <- 0 # setting some of the counts to 0 randomly
test_mod <- glm.nb(expr ~ time_point_num + genotype, data = test_genedf, link = log)


# TODO: it seems like when there are warning messages, the pvalue is super inflated. I wanted converged to indicate whether the warning was thrown or not, but it appears to simply always converge unless there's an error... Instead of converged I want a column that's only TRUE when no warnings about the iteration limit are thrown
# example to replicate this: gene_idx: YMR132C dataset: HAP4_hyb_HAP4del
# bad
test <- bind_cols(tibble(expr = es[[27]][3541,]), si[[27]])
badpois <- glm(expr ~ time_point_num, test, family = poisson(link = "log"))
MASS::theta.ml(badpois, weights = badpois$weights, trace = TRUE)
summary(glm.nb(expr ~ time_point_num, test, init.theta = t0bad))
summary(glm.nb(expr ~ time_point_num, test))
# good
test <- bind_cols(tibble(expr = es[[1]][3541,]), si[[1]])
goodpois <- glm(expr ~ time_point_num, test, family = poisson(link = "log"))
MASS::theta.ml(goodpois, weights = goodpois$weights, trace = TRUE)
summary(glm.nb(expr ~ time_point_num, test))

# seems to have to do with the high number of zero counts driving the delta for each iteration of theta way off course:
# del = score/info for each iter, and score is generally the same for both good and bad
info <- function(n, th, mu, y, w) {
  sum(w * (-trigamma(th + y) + trigamma(th) - 1/th + 2/(mu + th) - (y + th)/(mu + th)^2))
}
info(sum(goodpois$weights), 0.6, goodpois$fitted.values, goodpois$y, goodpois$weights)
info(sum(badpois$weights), 0.6, badpois$fitted.values, badpois$y, badpois$weights)
info_unsummed <- function(n, th, mu, y, w) {
  (-trigamma(th + y) + trigamma(th) - 1/th + 2/(mu + th) - (y + th)/(mu + th)^2)
}
goodkernel <- info_unsummed(sum(goodpois$weights), 0.6, goodpois$fitted.values, goodpois$y, goodpois$weights)
badkernel <- info_unsummed(sum(badpois$weights), 0.6, badpois$fitted.values, badpois$y, badpois$weights)
sum(goodkernel < 0)/length(goodkernel)
sum(badkernel < 0)/length(badkernel)

###############################################################

# Applies fitSingleGeneTPGenotypeModel to entire count matrix
# @input: count matrix and corresponding sample info
# @output: list of same length as nGenes, with the 3 named elements described in fitTPGenotypeModelAllGenes
fitTPGenotypeModelAllGenes <- function(.cts, .info) {
  output <- map(c(1:nrow(.cts)), function(i) {
    x <- .cts[i,]
    gname <- rownames(.cts)[i]
    genedf <- bind_cols(tibble(expr = x), .info)
    result <- fitSingleGeneTPGenotypeModel(genedf, .gene_name = gname)
    return(result)
  })
  names(output) <- rownames(.cts) # ANNA: I haven't tested this line as of 4/18/23, but I want a way to have each nbdf be able to be indexed by gene name (ex nbdfs$TFdelxLowN_cer$YGR192C)
  return(output)
}
# tests for fitTPGenotypeModelAllGenes
test <- fitTPGenotypeModelAllGenes(es$LowPi_cer, si$LowPi_cer) # this dataset has a number of genes where theta can't be estimated (likely due to high zero counts)
sum(is.na(test))
summary(test[[sample(c(1:length(test)), 1)]])$coefficients
# here's one that failed the map implementation initially but now seems (hopefully) fine:
test <- fitTPGenotypeModelAllGenes(es$HAP4_cer, si$HAP4_cer)
sum(is.na(test))
# summary(test[[sample(c(1:length(test)), 1)]])$coefficients # sometimes runs into errors for missing genes

#### running env-specific models on data ####
# Takes about 1 min per dataset (TFdelxLowN ones are longer tho):
modnbs <- map(c(1:length(es)), function(i) {
  cat(i, "/", nDatasets, "\n")
  cat("currently fitting", DatasetNames[i], "\n")
  return(fitTPGenotypeModelAllGenes(es[[i]], si[[i]]))
  })
# add names
names(modnbs) <- DatasetNames

# converting models into dataframes, extracting pvalues and effect sizes for genotype or timepoint coefficients
lookup <- tibble(dataset_name = DatasetNames,
                 organism = c(rep("cer", 4), rep("par", 4), rep("hyb", 4)),
                 experiment = lapply(si, function(x) {return(unique(x$experiment))}) %>% unlist())
parseDFfromModel <- function(.mod, .experiment, .gene_name) {
  if (.experiment == "LowN") {
    generows <- summary(.mod)$coefficients[grepl("^genotype", rownames(summary(.mod)$coefficients)),, drop = FALSE]
    result <- tibble(gene_name = .gene_name,
                     experiment = .experiment,
                     coefficient = gsub("genotype", "", rownames(generows)),
                     effect_size = as.numeric(generows[,"Estimate"]),
                     pvalue = as.numeric(generows[,"Pr(>|z|)"]))
  }
  if (.experiment != "LowN") {
    timerow <- summary(.mod)$coefficients[grepl("^time_point_num$", rownames(summary(.mod)$coefficients)),, drop = FALSE]
    result <- tibble(gene_name = .gene_name,
                     experiment = .experiment,
                     coefficient = "timepoint",
                     effect_size = as.numeric(timerow[,"Estimate"]),
                     pvalue = as.numeric(timerow[,"Pr(>|z|)"]))
  }
  return(result)
}
# tests for parseDFfromModel:
parseDFfromModel(modnbs$TFdelxLowN_cer[[sample(c(1:nGenes), 1)]], .experiment = "LowN", .gene_name = "LowN_test")
parseDFfromModel(modnbs$LowPi_hyb[[sample(c(1:nGenes), 1)]], .experiment = "LowPi", .gene_name = "LowPi_test")
parseDFfromModel(modnbs$HAP4_par[[sample(c(1:nGenes), 1)]], .experiment = "HAP4", .gene_name = "HAP4_test")
parseDFfromModel(modnbs$CC_cer[[sample(c(1:nGenes), 1)]], .experiment = "CC", .gene_name = "CC_test")
# doesn't work for NA (that's ok we'll skip them):
parseDFfromModel(modnbs$LowPi_par[[which(is.na(modnbs$LowPi_par))[1]]], .experiment = "LowPi", .gene_name = "NA_test") 
# applying it to our list of lists, converting to just a list
nbdfs <- lapply(DatasetNames, function(d) {
  cat("currently processing", d, "\n")
  x <- modnbs[[d]]
  organism <- lookup %>% filter(dataset_name == d) %>% select(organism) %>% pull()
  experiment <- lookup %>% filter(dataset_name == d) %>% select(experiment) %>% pull()
  output <- lapply(which(!is.na(x)), function(i) {
    gene_name <- names(x)[i]
    return(parseDFfromModel(x[[i]], .experiment = experiment, .gene_name = gene_name))
  }) %>% Reduce(f = bind_rows)
  output$experiment <- experiment
  output$organism <- organism
  return(output)
})
names(nbdfs) <- DatasetNames

# TODO: also create dfs based on marginal effects
modToMargEffect <- function(.mod, .genedf) {
  experiment <- .genedf$experiment[1]
  gene_name <- .genedf$gene_name[1]
  if (experiment == "LowN") {
    avgcomps <- avg_comparisons(.mod, variables = "genotype", newdata = .genedf)
    output <- tibble(gene_name = gene_name,
                     experiment = experiment,
                     coefficient = gsub(" - WT", "", avgcomps$contrast),
                     estimate = avgcomps$estimate,
                     pvalue = avgcomps$p.value)
  }
  if (.genedf$experiment[1] != "LowN") {
    avgcomps <- avg_comparisons(.mod, variables = "time_point_num", newdata = .genedf)
    output <- tibble(gene_name = gene_name,
                     experiment = experiment,
                     coefficient = "timepoint",
                     estimate = avgcomps$estimate,
                     pvalue = avgcomps$p.value)
  }
  return(output)
}
# tests for modToMargEffect
# LowPi
testmod <- modnbs$LowPi_cer[["YBR162C"]]
genedf <- bind_cols(tibble(expr = es$LowPi_cer["YBR162C",],
                           gene_name = "YBR162C"),
                    si$LowPi_cer)
ggplot(genedf, aes(x = time_point_num, y = expr)) + geom_point() + scale_x_discrete(breaks = NULL)
modToMargEffect(testmod, genedf)
# TFdel
testmod <- modnbs$TFdelxLowN_par[["YBR162C"]]
genedf <- bind_cols(tibble(expr = es$TFdelxLowN_par["YBR162C",],
                           gene_name = "YBR162C"),
                    si$TFdelxLowN_par)
ggplot(genedf, aes(x = genotype, y = expr)) + geom_point() + scale_x_discrete(breaks = NULL) # really looks like nothing should be significant
medf <- modToMargEffect(testmod, genedf) # takes a liiiitle longer
genedf <- left_join(genedf, select(medf, coefficient, estimate, pvalue), by = c("genotype" = "coefficient"))
genedf$sig <- genedf$pvalue < pthresh
genedf$genotype <- as.factor(genedf$genotype) %>% relevel(ref = "WT")
ggplot(genedf, aes(x = genotype, y = expr)) + geom_point(aes(color = sig)) + scale_x_discrete(breaks = NULL)

# DONT RUN: overnight job for sure if you dare, prob not worth it unless the QC below finds it's significantly better than just inspecting the glm coefficient
# medfs <- lapply(DatasetNames, function(d) {
#   cat("currently processing", d, "\n")
#   x <- modnbs[[d]] # list of models, length nGenes with names in GeneNames
#   experiment <- lookup %>% filter(dataset_name == d) %>% select(experiment) %>% pull()
#   organism <- lookup %>% filter(dataset_name == d) %>% select(organism) %>% pull()
#   output <- lapply(which(!is.na(x)), function(i) {
#     gname <- names(x)[i]
#     cat("working on gene:", gname, i, "/", nGenes, "\n")
#     genedf <- bind_cols(tibble(expr = es[[d]][gname,],
#                                gene_name = gname),
#                         si[[d]])
#     return(modToMargEffect(.mod = x[[i]], .genedf = genedf))
#   }) %>% Reduce(f = bind_rows)
#   output$experiment <- experiment
#   output$organism <- organism
#   return(output)
# })
# names(medfs) <- DatasetNames

# creating the big dataset of all the pvalues and effect sizes calculated
# 5000 genes * (3 starvation + 46 TFdel experiments) * 3 species * 1-2 replicates = roughly 710000 rows
resultdf <- Reduce(nbdfs, f = bind_rows) %>% 
  mutate(sig = if_else(experiment == "LowN",
                       true = pvalue < pthresh & abs(effect_size) > MinEffectSizeCatagorical,
                       false = pvalue < pthresh & abs(effect_size) > MinEffectSizeContinuous))

sum(resultdf$sig, na.rm = TRUE)
save(resultdf, file = "data_files/env_spec_DEdf.RData")

#### Quality Control ####
# TODO: Anecdotally it seems like marginal effects is better at calling DE genes when they're resulting in LOWER expression than the reference and the standard wald test is better when they're higher than the reference, so check if there's a bias for negative estimates in medfs and positive estimates in nbdfs
# Env Stress (note these are biased for genes called DE in glm)
qcdf <- tibble(glmeffsize = NULL, glmpval = NULL, meeffsize = NULL, mepval = NULL, experiment = NULL, organism = NULL)
for (i in c(1:1000)) {
  cat (i, "\n")
  dataset_idx <- sample(c("CC_cer", "HAP4_cer", "LowPi_cer",
                          "CC_par", "HAP4_par", "LowPi_par",
                          "CC_hyb", "HAP4_hyb", "LowPi_hyb"), 1) # skipping TFdel for now cause they're soo long
  exprmt <- filter(lookup, dataset_name == dataset_idx) %>% select(experiment) %>% pull() %>% as.character()
  org <- filter(lookup, dataset_name == dataset_idx) %>% select(organism) %>% pull() %>% as.character()
  gene_idx <- resultdf %>% filter(experiment == exprmt & organism == org & sig) %>% 
    select(gene_name) %>% pull() %>% sample(size = 1)
  margeff <- modToMargEffect(modnbs[[dataset_idx]][[gene_idx]], .genedf = bind_cols(tibble(expr = es[[dataset_idx]][gene_idx,],
                                                                                           gene_name = gene_idx),
                                                                                    si[[dataset_idx]]))
  glmcoeff <- nbdfs[[dataset_idx]] %>% filter(gene_name == gene_idx)
  qcdf <- bind_rows(qcdf, tibble(glmeffsize = glmcoeff$effect_size,
                                 glmpval = glmcoeff$pvalue,
                                 meeffsize = margeff$estimate,
                                 mepval = margeff$pvalue,
                                 experiment = exprmt,
                                 organism = org))
}
sum(qcdf$glmpval < pthresh & qcdf$mepval < pthresh)
sum(qcdf$glmpval < pthresh & qcdf$mepval >= pthresh)
ggplot(qcdf, aes(x = log(abs(glmeffsize)), y = log(meeffsize))) + geom_point(aes(col = mepval < pthresh))
# No real bias in sight, but they're def not the same thing (I mean obvs b/c the ME "effect size" is actually an expr estimate)

# TFdel (note these are biased for genes called DE in glm)
qcdf <- tibble(glmeffsize = NULL, glmpval = NULL, meeffsize = NULL, mepval = NULL, experiment = NULL, organism = NULL)
for (i in c(1:10)) {
  cat (i, "\n")
  dataset_idx <- sample(c("TFdelxLowN_cer", "TFdelxLowN_par", "TFdelxLowN_hyb"), 1) # skipping TFdel for now cause they're soo long
  exprmt <- filter(lookup, dataset_name == dataset_idx) %>% select(experiment) %>% pull() %>% as.character()
  org <- filter(lookup, dataset_name == dataset_idx) %>% select(organism) %>% pull() %>% as.character()
  gene_idx <- resultdf %>% filter(experiment == exprmt & organism == org & sig) %>% 
    select(gene_name) %>% pull() %>% sample(size = 1)
  margeff <- modToMargEffect(modnbs[[dataset_idx]][[gene_idx]], .genedf = bind_cols(tibble(expr = es[[dataset_idx]][gene_idx,],
                                                                                           gene_name = gene_idx),
                                                                                    si[[dataset_idx]]))
  glmcoeff <- nbdfs[[dataset_idx]] %>% filter(gene_name == gene_idx)
  qcdf <- bind_rows(qcdf, tibble(glmeffsize = glmcoeff$effect_size,
                                 glmpval = glmcoeff$pvalue,
                                 meeffsize = margeff$estimate,
                                 mepval = margeff$pvalue,
                                 experiment = exprmt,
                                 organism = org))
}
ggplot(qcdf, aes(x = log(abs(glmeffsize)), y = log(meeffsize))) + geom_point(aes(col = mepval < pthresh))
# Ah. Here's the bias. Has to be huge to be detected in margeff

# TODO: Verify that the reason YLR307W was only significant in glm.nb not marginal effects was because of the datagrid excluding conditions and losing power (by running it with avg_comparisons)

# visualizing genes by pvalue
dataset_idx <- sample(c(1:length(es)), 1)
datasetdf <- nbdfs[[dataset_idx]]
datasetdf$mean_expr <- sapply(datasetdf$gene_name, function(g) {
  return(mean(es[[dataset_idx]][g,], na.rm = TRUE))
})
ggplot(datasetdf, aes(x = log(mean_expr), y = -log(pvalue))) + 
  geom_point(alpha = 0.5, aes(color = coefficient, size = abs(effect_size))) + 
  geom_hline(yintercept = -log(pthresh), color = "red") + theme(legend.position = "none")

# now plot nsig genes vs nsamples where each point is a dataset
nsamples <- lapply(es, ncol)
nSig <- map2(nbdfs, names(nbdfs), function(x, y) {
  cat("calculating nsig for", y, "\n")
  if (grepl("LowN", y)) {
    mineffect <- MinEffectSizeCatagorical
    nsig <- x %>% group_by(gene_name) %>% summarise(geneDE = any(pvalue < pthresh & abs(effect_size) > mineffect, na.rm = TRUE)) %>% select(geneDE) %>% pull() %>% sum(na.rm = TRUE)
  } else {
    mineffect <- MinEffectSizeContinuous
    nsig <- x %>% group_by(gene_name) %>% summarise(geneDE = pvalue < pthresh & abs(effect_size) > mineffect) %>% select(geneDE) %>% pull() %>% sum(na.rm = TRUE)
  }
  return(nsig)
})
powerdf <- tibble(n_samples = unlist(nsamples), n_significant = unlist(nSig))
powerdf$experiment <- lapply(si, function(x) {return(unique(x$experiment))}) %>% unlist()
powerdf$organism <- lapply(si, function(x) {return(unique(x$organism))}) %>% unlist()
ggplot(powerdf, aes(x = n_samples, y = n_significant)) + geom_point(aes(color = experiment, shape = organism == "hyb"))
lm(n_significant ~ n_samples, powerdf) %>% summary() # Rsquarred is measily, as is nonsignificant effect size, we prob shouldn't use LowN with WT only

### Visualizing Individual Genes ###

# QC: making sure each TF deletion has its corresponding genes' expression tanked
gene_lookup <- read_delim("data_files/46TFdelLookup.csv", delim = ";", # there are actually 45 lol. WT is the 46th condition
                          col_names = FALSE, col_select = c(1,2))
colnames(gene_lookup) <- c("common_name", "systematic_name")
for (i in c(1:nrow(gene_lookup))) {
  cat("working on", gene_lookup$common_name[i], "\n")
  genedel <- paste0(gene_lookup$common_name[i], "delete")
  gname <- gene_lookup$systematic_name[i]
  gdf <- filter(resultdf, coefficient == genedel & gene_name == gname)
  issig <- gdf %>% select(sig) %>% pull()
  if (!all(issig)) {
    cat(gene_lookup$common_name[i], "not all sig:\n")
    print(gdf[!issig,], n = sum(!issig))
  }
}
# re-run this whole chunk to look at different TFdels
lookup_idx <- which(gene_lookup$common_name == "MBP1") # run with any number 1-45
gene_idx <- gene_lookup$systematic_name[lookup_idx]
genedel <- paste0(gene_lookup$common_name[lookup_idx], "delete")
genedf_cer <- bind_cols(tibble(expr = es$TFdelxLowN_cer[gene_idx,]),
                        si$TFdelxLowN_cer)
sigs <- filter(resultdf, experiment == "LowN" & organism == "cer" & gene_name == gene_idx) %>% select(coefficient, sig)
genedf_cer <- left_join(genedf_cer, sigs, by = c("genotype" = "coefficient"))
genedf_cer$genotype <- as.factor(genedf_cer$genotype) %>% relevel(ref = "WT")
cer_plot <- ggplot(genedf_cer, aes(x = genotype, y = expr)) + geom_point(aes(color = genotype == genedel)) + scale_x_discrete(breaks = NULL) + theme(legend.position = "none")
# par
genedf_par <- bind_cols(tibble(expr = es$TFdelxLowN_par[gene_idx,]),
                        si$TFdelxLowN_par)
sigs <- filter(resultdf, experiment == "LowN" & organism == "par" & gene_name == gene_idx) %>% select(coefficient, sig)
genedf_par <- left_join(genedf_par, sigs, by = c("genotype" = "coefficient"))
genedf_par$genotype <- as.factor(genedf_par$genotype) %>% relevel(ref = "WT")
par_plot <- ggplot(genedf_par, aes(x = genotype, y = expr)) + geom_point(aes(color = genotype == genedel)) + scale_x_discrete(breaks = NULL) + theme(legend.position = "none")
# hyb
genedf_hyb <- bind_cols(tibble(expr = es$TFdelxLowN_hyb[gene_idx,]),
                        si$TFdelxLowN_hyb)
sigs <- filter(resultdf, experiment == "LowN" & organism == "par" & gene_name == gene_idx) %>% select(coefficient, sig)
genedf_hyb <- left_join(genedf_hyb, sigs, by = c("genotype" = "coefficient"))
genedf_hyb$genotype <- as.factor(genedf_hyb$genotype) %>% relevel(ref = "WT")
hyb_plot <- ggplot(genedf_hyb, aes(x = genotype, y = expr)) + geom_point(aes(color = genotype == genedel)) + scale_x_discrete(breaks = NULL) + theme(legend.position = "none")
library(gridExtra)
grid.arrange(cer_plot, par_plot, hyb_plot, nrow = 1)
# List of bad ones: INO4, CBF1, FKH1, MBP1 (sort of), RFX1, TEC1, ZAP1 (sort of, just lowly expressed)
# some of: GCN4 (some samples in paradoxus have horribly high counts, possible contamination), HAP2 (one in cer), ROX1 (in cer), RPN4 (one in cer)
# Good example of one that should def be sig: MET28/YIR017C
# Best ones: MET28, SOK2, URE2, YAP1

# Is MET28 recognized as significantly down-expressed in MET28delete?
resultdf %>% filter(coefficient == "MET28delete" & gene_name == "YIR017C") # not based on coefficient alone...
metmodcer <- modnbs$TFdelxLowN_cer[["YIR017C"]]
genedf <- bind_cols(tibble(expr = es$TFdelxLowN_cer["YIR017C",]),
                    si$TFdelxLowN_cer)
avg_comparisons(metmodcer, variables = list(genotype = c("WT", "MET28delete")), newdata = datagrid(model = metmodcer, genotype = genedf$genotype, timepoint = genedf$timepoint))
# Yes it is!

# my "Manhattan Plot" of number of DE genes in each experimental condition
single_perturbations <- resultdf %>% 
  filter(sig) %>% 
  group_by(organism, experiment, coefficient) %>% 
  summarise(nsig = length(unique(gene_name)), sig_genes = list(unique(gene_name)))
gene_deletion_names <- filter(resultdf, experiment == "LowN" & !(coefficient %in% c("gene_var_zero", "model_fit_failed"))) %>% 
  select(coefficient) %>% pull() %>% unique() %>% sort() %>% as.character()
ggplot(single_perturbations, aes(x = paste(experiment, coefficient), y = nsig, fill = organism)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_x_discrete(limits = c("CC timepoint", "HAP4 timepoint", "LowPi timepoint", paste("LowN", gene_deletion_names)),
                   labels = c("Cell Cycle", "Growth Curve", "Low Pi", gene_deletion_names)) + 
  xlab(element_blank()) + ylab("number of DE genes") + theme_classic() + theme(legend.position = "bottom") + scale_fill_manual(values = c("gold", "chartreuse3", "mediumorchid"), breaks = c("cer", "hyb", "par")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

visualizeGeneExpression <- function(.name, .experiment = c("LowN", "LowPi", "CC", "HAP4"), .organism = c("cer", "par", "hyb")) {
  if (.experiment == "LowN") {
     dataset_name <- paste("TFdelxLowN", .organism, sep = "_")
     genedf <- bind_cols(tibble(expr = es[[dataset_name]][.name,]),
                         si[[dataset_name]])
     sigs <- filter(resultdf, experiment == .experiment & organism == .organism & gene_name == .name) %>% select(coefficient, sig)
     genedf <- left_join(genedf, sigs, by = c("genotype" = "coefficient"))
     genedf$genotype <- as.factor(genedf$genotype) %>% relevel(ref = "WT")
     output <- ggplot(genedf, aes(x = genotype, y = expr)) + 
       geom_point(aes(color = sig)) + 
       ggtitle(.name) + xlab("genotype") + ylab("expression level") + 
       theme_classic() + theme(legend.title = element_blank()) +
       scale_x_discrete(breaks = NULL) +
       scale_color_manual(breaks = c(TRUE, FALSE, NA), values = c("#F8766D", "#619CCF", "#00BA38"), labels = c("DE", "not DE", "WT (ref)"))
  }
  if (.experiment != "LowN") {
    dataset_name <- paste(.experiment, .organism, sep = "_")
    genedf <- bind_cols(tibble(expr = es[[dataset_name]][.name,]),
                        si[[dataset_name]])
    sig <- filter(resultdf, experiment == .experiment & organism == .organism & gene_name == .name) %>% select(sig) %>% pull()
    output <- ggplot(genedf, aes(x = time_point_num, y = expr)) + 
      geom_point(aes(color = sig)) + 
      ggtitle(.name) + xlab("timepoint") + ylab("expression level") + 
      theme_classic() + theme(legend.title = element_blank()) +
      scale_color_manual(breaks = c(TRUE, FALSE, NA), values = c("#F8766D", "#619CCF", "#00BA38"), labels = c("DE", "not DE", "Missing"))
  }
  return(output)
}
visualizeGeneExpression("YGR192C", "HAP4", "hyb")
visualizeGeneExpression("YGR192C", "LowPi", "par")
visualizeGeneExpression("YOL151W", "LowN", "hyb")
visualizeGeneExpression("YKL071W", "CC", "par")

getRandomSigGeneName <- function(.experiment = c("LowN", "LowPi", "CC", "HAP4", "all"), .organism = c("cer", "par", "hyb")) {
  df <- resultdf %>% filter(experiment == .experiment & organism == .organism)
  if (.experiment == "LowN") {
    sig_genes <- df %>% group_by(gene_name) %>% summarise(genesig = any(sig)) %>% ungroup() %>% filter(genesig) %>% select(gene_name) %>% pull() %>% unique()
    return(sample(sig_genes, 1))
  }
  if (.experiment != "LowN") {
    sig_genes <- df %>% filter(sig) %>% select(gene_name) %>% pull() %>% unique()
    return(sample(sig_genes, 1))
  }
}
# tests for getRandomSigGene
random_name <- getRandomSigGeneName("HAP4", "par")
resultdf %>% filter(experiment == "HAP4" & organism == "par" & gene_name == random_name)
random_name <- getRandomSigGeneName("LowPi", "cer")
resultdf %>% filter(experiment == "LowPi" & organism == "cer" & gene_name == random_name)
random_name <- getRandomSigGeneName("LowN", "hyb")
resultdf %>% filter(experiment == "LowN" & organism == "hyb" & gene_name == random_name & pvalue < pthresh & abs(effect_size) > MinEffectSizeCatagorical)
random_name <- getRandomSigGeneName("CC", "par")
resultdf %>% filter(experiment == "CC" & organism == "par" & gene_name == random_name)

getRandomNonsigGeneName <- function(.experiment = c("LowN", "LowPi", "CC", "HAP4"), .organism = c("cer", "par", "hyb")) {
  df <- resultdf %>% filter(experiment == .experiment & organism == .organism)
  if (.experiment == "LowN") {
    nonsig_genes <- df %>% group_by(gene_name) %>% summarise(genenonsig = !any(sig)) %>% ungroup() %>% filter(genenonsig) %>% select(gene_name) %>% pull() %>% unique()
    return(sample(nonsig_genes, 1))
  }
  if (.experiment != "LowN") {
    nonsig_genes <- df %>% filter(!sig) %>% select(gene_name) %>% pull() %>% unique()
    return(sample(nonsig_genes, 1))
  }
}
# tests for getRandomNonsigGene
random_name <- getRandomNonsigGeneName("HAP4", "par")
resultdf %>% filter(experiment == "HAP4" & organism == "par" & gene_name == random_name)
random_name <- getRandomNonsigGeneName("LowPi", "cer")
resultdf %>% filter(experiment == "LowPi" & organism == "cer" & gene_name == random_name)
random_name <- getRandomNonsigGeneName("LowN", "hyb")
resultdf %>% filter(experiment == "LowN" & organism == "hyb" & gene_name == random_name & pvalue < pthresh & abs(effect_size) > MinEffectSizeCatagorical) # filtering for sig b/c it's LowN, shouldn't have anything
random_name <- getRandomNonsigGeneName("CC", "par")
resultdf %>% filter(experiment == "CC" & organism == "par" & gene_name == random_name)

# 1B inset: plotting random sig and non-sig genes, 1 from LowPi and 1 from TFdel (4 plots total)
# LowPi sig
sig_LowPi <- getRandomSigGeneName(.experiment = "LowPi", .organism = "hyb")
visualizeGeneExpression(sig_LowPi, .experiment = "LowPi", .organism = "hyb")
# LowPi nonsig
nonsig_LowPi <- getRandomNonsigGeneName(.experiment = "LowPi", .organism = "hyb")
visualizeGeneExpression(nonsig_LowPi, .experiment = "LowPi", .organism = "hyb")
# HAP4 sig
sig_HAP4 <- getRandomSigGeneName(.experiment = "HAP4", .organism = "hyb")
visualizeGeneExpression(sig_HAP4, .experiment = "HAP4", .organism = "hyb")
# HAP4 nonsig
nonsig_HAP4 <- getRandomNonsigGeneName(.experiment = "HAP4", .organism = "hyb")
visualizeGeneExpression(nonsig_HAP4, .experiment = "HAP4", .organism = "hyb")
# CC sig
sig_CC <- getRandomSigGeneName(.experiment = "CC", .organism = "hyb")
visualizeGeneExpression(sig_CC, .experiment = "CC", .organism = "hyb")
# CC nonsig
nonsig_CC <- getRandomNonsigGeneName(.experiment = "CC", .organism = "hyb")
visualizeGeneExpression(nonsig_HAP4, .experiment = "CC", .organism = "hyb")
# TFdel sig
sig_TFdel <- getRandomSigGeneName(.experiment = "LowN", .organism = "hyb")
visualizeGeneExpression(sig_TFdel, .experiment = "LowN", .organism = "hyb")
# TFdel nonsig
nonsig_TFdel <- getRandomNonsigGeneName(.experiment = "LowN", .organism = "hyb")
visualizeGeneExpression(nonsig_TFdel, .experiment = "LowN", .organism = "hyb")

plot_list <- list(visualizeGeneExpression(sig_LowPi, .experiment = "LowPi", .organism = "hyb"),
                  visualizeGeneExpression(sig_HAP4, .experiment = "HAP4", .organism = "hyb"),
                  visualizeGeneExpression(sig_CC, .experiment = "CC", .organism = "hyb"),
                  visualizeGeneExpression(sig_TFdel, .experiment = "LowN", .organism = "hyb"),
                  visualizeGeneExpression(nonsig_LowPi, .experiment = "LowPi", .organism = "hyb"),
                  visualizeGeneExpression(nonsig_HAP4, .experiment = "HAP4", .organism = "hyb"),
                  visualizeGeneExpression(nonsig_CC, .experiment = "CC", .organism = "hyb"),
                  visualizeGeneExpression(nonsig_TFdel, .experiment = "LowN", .organism = "hyb"))
figure <- ggarrange(plotlist = plot_list, nrow = 2, ncol = 4)
annotate_figure(figure, top = text_grob("Order: LowPi Growth Curve Cell Cycle LowN", color = "black", face = "bold")) + theme(legend.position = "bottom")

# Good luck trying to do this (it's like 10GB compressed)
# save(modnbs, file = "data_files/env_specific_DE.RData")

############ Testing for inter-species and inter-allele DE #############
# species datasets, each sample has a count for each gene that's either cer or par
spcts <- Reduce(f = cbind, list(es$TFdelxLowN_cer, es$CC_cer, es$HAP4_cer, es$LowPi_cer,
                es$TFdelxLowN_par, es$CC_par, es$HAP4_par, es$LowPi_par))
spinfo <- Reduce(f = bind_rows, list(si$TFdelxLowN_cer, si$CC_cer, si$HAP4_cer, si$LowPi_cer,
                                     si$TFdelxLowN_par, si$CC_par, si$HAP4_par, si$LowPi_par))
sum(spinfo$sample_name == colnames(spcts))/nrow(spinfo)

# hybrid allele datasets
load("data_files/Cleaned_Barkai_Data_AlleleSpecific.RData")
alcts <- counts_allele
alinfo <- sample_info_allele
rm(counts_allele, sample_info_allele)
sum(alinfo$sample_name == colnames(alcts))/nrow(alinfo)

# Reduce must have messed up the rownames
rownames(spcts) <- GeneNames

# these are much easier with only two datasets lol
# counts per million
# spcts <- countsPerMillion(spcts) # pretty sure we already did this up above
alcts <- countsPerMillion(alcts)
# factorize genotype, experiment, and time_point_str (no specific reference level for the laste two)
spinfo$genotype <- as.factor(spinfo$genotype) %>% relevel(ref = "WT")
alinfo$genotype <- as.factor(alinfo$genotype) %>% relevel(ref = "WT")
spinfo$experiment <- as.factor(spinfo$experiment) %>% relevel(ref = "LowN")
alinfo$experiment <- as.factor(alinfo$experiment) %>% relevel(ref = "LowN")
spinfo$time_point_str <- as.factor(spinfo$time_point_str)
alinfo$time_point_str <- as.factor(alinfo$time_point_str)
spinfo$allele <- as.factor(spinfo$allele) %>% relevel(ref = "cer")
alinfo$allele <- as.factor(alinfo$allele) %>% relevel(ref = "cer")

modToMargEffectSpeciesAllele <- function(.mod, .genedf) {
  avgcomps <- avg_comparisons(.mod, variables = "allele", newdata = .genedf)
  output <- tibble(gene_name = .genedf$gene_name[1],
                   coefficient = .genedf$allele_or_species[1],
                   estimate = avgcomps$estimate,
                   pvalue = avgcomps$p.value)
  return(output)
}

# test model on a single gene
# allele
random_gene <- sample(GeneNames, 1)
algenedf <- bind_cols(tibble(expr = alcts[random_gene,],
                             gene_name = random_gene,
                             allele_or_species = "allele"),
                      alinfo)
plotdf <- algenedf %>% mutate(non_unique_sample_name = gsub("_hy._", "_hyb_", sample_name)) %>% select(expr, non_unique_sample_name, allele, gene_name) %>% pivot_wider(id_cols = c(non_unique_sample_name, gene_name), values_from = expr, names_from = allele)
ggplot(plotdf, aes(x = cer, y = par)) + geom_point() + geom_abline(color = "gold") + ggtitle(paste(plotdf$gene_name[1], "hybrid alleles"))
almod <- glm.nb(expr ~ genotype + experiment + time_point_str + allele, algenedf, link = log)
summary(almod)
alrow <- summary(almod)$coefficients[grepl("allele", rownames(summary(almod)$coefficients)),, drop = FALSE]
alme <- modToMargEffectSpeciesAllele(.mod = almod, .genedf = algenedf)

# species
spgenedf <- bind_cols(tibble(expr = spcts[random_gene,],
                             gene_name = random_gene,
                             allele_or_species = "species"),
                      spinfo)
# each dot is a condition-matched sample (arbitrary replicates paired, some replicates removed)
plotdf <- spgenedf %>% mutate(non_unique_sample_name = paste(condition, well_flask_ID)) %>% select(expr, non_unique_sample_name, allele, gene_name) %>% pivot_wider(id_cols = c(non_unique_sample_name, gene_name), values_from = expr, names_from = allele)
ggplot(plotdf, aes(x = cer, y = par)) + geom_point() + geom_abline(color = "gold") + ggtitle(paste(plotdf$gene_name[1], "parental species"))
spmod <- glm.nb(expr ~ genotype + experiment + time_point_str + allele, spgenedf, link = log)
summary(spmod)
sprow <- summary(spmod)$coefficients[grepl("allele", rownames(summary(spmod)$coefficients)),, drop = FALSE]
spme <- modToMargEffectSpeciesAllele(.mod = spmod, .genedf = spgenedf)
rbind(alrow, sprow)
bind_rows(alme, spme) # re-run a few times with random genes. Given what this gene looks like in the plots, do you trust either strategy better? YNL124W is a good cis example, YOR087W is a good trans example.

# or visualize counts in aggregate, not condition-matched
genedf <- bind_rows(spgenedf, algenedf)
ggplot(genedf, aes(x = factor(paste(organism, allele), levels = c("cer cer", "par par", "hyb cer", "hyb par")), y = expr)) + 
  geom_boxplot(aes(fill = factor(paste(organism, allele), levels = c("cer cer", "par par", "hyb cer", "hyb par")))) + 
  scale_fill_manual(breaks = c("cer cer", "par par", "hyb cer", "hyb par"),
                      values = c("gold","mediumorchid","lemonchiffon", "plum")) +
  theme(legend.position = "none") + xlab(element_blank())


fitSingleGeneSpeciesAlleleModel <- function(.df, .gene_name) {
  output <- tryCatch({
    mod <- glm.nb(expr ~ genotype + experiment + time_point_str + allele, .df, link = log)
    return(mod)
  }, error = function(e) {
    return(NA)
  }, warning = function(w) {
    return(NA)
  })
  return(output)
}
gene_idx <- sample(GeneNames, 1) 
genedf <- bind_cols(tibble(expr = spcts[gene_idx,]),
                    spinfo)
fitSingleGeneSpeciesAlleleModel(genedf, gene_idx) %>% summary()

# fitting the models and getting the results in two forms - extracting allele coefficient (nbdfs) and computing marginal effect of allele (medfs)
# species
spmods <- map(GeneNames, function(g) {
  cat(which(GeneNames == g), "/", nGenes, "\n")
  cat("currently processing", g, "\n")
  gdf <- bind_cols(tibble(expr = spcts[g,]),
                   spinfo)
  output <- fitSingleGeneSpeciesAlleleModel(gdf, g)
  return(output)
})
names(spmods) <- GeneNames

# should be v fast
spdf <- lapply(which(!is.na(spmods)), function(i) {
  g <- GeneNames[i]
  cat(which(GeneNames == g), "/", nGenes, "\n")
  cat("currently processing", g, "\n")
  mod <- spmods[[i]]
  allelerow <- summary(mod)$coefficients[grepl("allele", rownames(summary(mod)$coefficients)),, drop = FALSE]
  result <- tibble(gene_name = g,
                   coefficient = rownames(allelerow), # precaution to make sure cer is always reference level
                   effect_size = as.numeric(allelerow[,"Estimate"]),
                   pvalue = as.numeric(allelerow[,"Pr(>|z|)"]))
  return(result)
}) %>% Reduce(f = bind_rows)

# should be more like running glm.nb
# spme <- lapply(which(!is.na(spmods)), function(i) {
#   g <- GeneNames[i]
#   cat("working on", g, "\n")
#   cat(which(GeneNames == g), "/", nGenes, "\n")
#   mod <- spmods[[i]]
#   gdf <- bind_cols(tibble(expr = spcts[g,],
#                           gene_name = g,
#                           allele_or_species = "species"),
#                    spinfo)
#   result <- modToMargEffectSpeciesAllele(.mod = mod, .genedf = gdf)
#   return(result)
# }) %>% Reduce(f = bind_rows)

# allele
almods <- map(GeneNames, function(g) {
  cat(which(GeneNames == g), "/", nGenes, "\n")
  cat("currently processing", g, "\n")
  gdf <- bind_cols(tibble(expr = alcts[g,]),
                   alinfo)
  output <- fitSingleGeneSpeciesAlleleModel(gdf, g)
  return(output)
})
names(almods) <- GeneNames

# should be v fast
aldf <- lapply(which(!is.na(almods)), function(i) {
  g <- GeneNames[i]
  cat(which(GeneNames == g), "/", nGenes, "\n")
  cat("currently processing", g, "\n")
  mod <- almods[[i]]
  allelerow <- summary(mod)$coefficients[grepl("allele", rownames(summary(mod)$coefficients)),, drop = FALSE]
  result <- tibble(gene_name = g,
                   coefficient = rownames(allelerow), # precaution to make sure cer is always reference level
                   effect_size = as.numeric(allelerow[,"Estimate"]),
                   pvalue = as.numeric(allelerow[,"Pr(>|z|)"]))
  return(result)
}) %>% Reduce(f = bind_rows)

# should be more like running glm.nb
# alme <- lapply(which(!is.na(almods)), function(i) {
#   g <- GeneNames[i]
#   cat("working on", g, "\n")
#   cat(which(GeneNames == g), "/", nGenes, "\n")
#   mod <- almods[[i]]
#   gdf <- bind_cols(tibble(expr = alcts[g,],
#                           gene_name = g,
#                           allele_or_species = "allele"),
#                    alinfo)
#   result <- modToMargEffectSpeciesAllele(.mod = mod, .genedf = gdf)
#   return(result)
# }) %>% Reduce(f = bind_rows)

spdf$coefficient <- gsub("allele", "species", spdf$coefficient)
spdf <- rename(spdf, species_effect_size = effect_size)
aldf <- rename(aldf, allele_effect_size = effect_size)
# Looking at the most diverged genes: taking the 95th percentile for magnitude of effect size for genes that are over a fairly arbitrary p-value threshold
# Rationale: we're not trying to ID every single DE gene, just get a large and un-biased enough sample to identify patterns
quantile(-log(spdf$pvalue), 0.75, na.rm = TRUE)
quantile(-log(aldf$pvalue), 0.75, na.rm = TRUE)
quantile(spdf$species_effect_size, 0.75, na.rm = TRUE)
quantile(aldf$allele_effect_size, 0.75, na.rm = TRUE)
# TODO: Now that we trust the models more, I don't think we want to have different cutoffs for allele vs species, as they're on the same scale, but we can pick more informed values than these:
spal_pthresh <- 1e-20
spal_effthresh <- 0.1

spdf$species_sig <- spdf$pvalue < spal_pthresh & abs(spdf$species_effect_size) > spal_effthresh
aldf$allele_sig <- aldf$pvalue < spal_pthresh & abs(aldf$allele_effect_size) > spal_effthresh
spaldf <- left_join(select(spdf, gene_name, species_sig, species_effect_size), select(aldf, gene_name, allele_sig, allele_effect_size), by = c("gene_name"))
sum(sign(spaldf$allele_effect_size) == sign(spaldf$species_effect_size), na.rm = TRUE)
sum(sign(spaldf$allele_effect_size) != sign(spaldf$species_effect_size), na.rm = TRUE)

spaldf %>% filter(species_sig | allele_sig) %>% summarise(nSameDirection = sum(sign(allele_effect_size) == sign(species_effect_size), na.rm = TRUE),
                                                          nOppositeDirection = sum(sign(allele_effect_size) != sign(species_effect_size), na.rm = TRUE))

visualizeGeneExpressionSpeciesAlleleBoxplot <- function(.name) {
  algenedf <- bind_cols(tibble(expr = alcts[.name,]),
                        alinfo)
  spgenedf <- bind_cols(tibble(expr = spcts[.name,]),
                        spinfo)
  algenedf$allele_or_species <- "allele"
  spgenedf$allele_or_species <- "species"
  genedf <- bind_rows(algenedf, spgenedf)
  genedf$type <- factor(paste(genedf$allele, genedf$allele_or_species), level = c("cer species", "par species", "cer allele", "par allele"))
  output <- ggplot(genedf, aes(x = type, y = log(expr))) + 
    geom_boxplot(aes(fill = type)) + 
    ggtitle(.name) + xlab(element_blank()) + ylab("expression level") + 
    theme_classic() + scale_fill_manual(values = c("gold","mediumorchid","lemonchiffon", "plum")) +
    theme(legend.position = "none")
  return(output)
}
visualizeGeneExpressionSpeciesAlleleBoxplot(gene_idx)

spaldf$divergence_mode <- map(c(1:nrow(spaldf)), function(i) {
  s <- spaldf[i,, drop = FALSE]
  if (is.na(s$species_sig) | is.na(s$allele_sig)) {
    return(NA)
  }
  if (s$species_sig & !s$allele_sig) {
    return("trans")
  }
  if (s$species_sig & s$allele_sig & (sign(s$species_effect_size) == sign(s$allele_effect_size))) {
    return("cis")
  }
  if (!s$species_sig & s$allele_sig) {
    return("compensated")
  }
  if (s$species_sig & s$allele_sig & (sign(s$species_effect_size) != sign(s$allele_effect_size))) {
    return("weirdos")
  }
  if (!s$species_sig & !s$allele_sig) {
    return("conserved")
  }
}) %>% unlist()
table(spaldf$divergence_mode)

random_trans_gene <- filter(spaldf, divergence_mode == "trans") %>% select(gene_name) %>% pull() %>% sample(size = 1)
visualizeGeneExpressionSpeciesAlleleBoxplot(random_trans_gene)

random_cis_gene <- filter(spaldf, divergence_mode == "cis") %>% select(gene_name) %>% pull() %>% sample(size = 1)
visualizeGeneExpressionSpeciesAlleleBoxplot(random_cis_gene)

random_weirdo_gene <- filter(spaldf, divergence_mode == "weirdos") %>% select(gene_name) %>% pull() %>% sample(size = 1)
visualizeGeneExpressionSpeciesAlleleBoxplot(random_weirdo_gene)

random_compensated_gene <- filter(spaldf, divergence_mode == "compensated") %>% select(gene_name) %>% pull() %>% sample(size = 1)
visualizeGeneExpressionSpeciesAlleleBoxplot(random_compensated_gene)

# Wittkopp plot of median expression of each gene colored by divergence mode (noting that this obscures most of what is so cool about this dataset, namely the environment-specific information)
plotdf <- tibble(gene_name = GeneNames,
                 parental_expr = apply(spcts[, spinfo$allele == "cer"], 1, median)/apply(spcts[, spinfo$allele == "par"], 1, median),
                 hybrid_expr = apply(alcts[, alinfo$allele == "cer"], 1, median)/apply(alcts[, alinfo$allele == "par"], 1, median)) %>% 
  left_join(y = select(spaldf, gene_name, divergence_mode), by = c("gene_name"))
library(scales)
DivergenceModeColors <- hue_pal()(5)
# full
plotfull <- ggplot(filter(plotdf, !is.na(divergence_mode)), aes(x = log(parental_expr), y = log(hybrid_expr))) + geom_point(aes(color = divergence_mode)) + 
  xlim(c(-2.5, 2.5)) + ylim(c(-2.5, 2.5)) + scale_color_manual(values = c("cis"=DivergenceModeColors[1], "compensated"=DivergenceModeColors[2], "conserved"=DivergenceModeColors[3], "trans"=DivergenceModeColors[4], "weirdos"=DivergenceModeColors[5]),
                                                               breaks = c("cis", "compensated", "conserved", "trans", "weirdos")) + ggtitle("Full") + xlab("") + ylab("") + theme(legend.position = "none") + theme_classic()

# trans only
plottrans <- ggplot(filter(plotdf, divergence_mode == "trans"), aes(x = log(parental_expr), y = log(hybrid_expr))) + geom_point(color = hue_pal()(5)[4]) + 
  xlab("") + ylab("") +
  xlim(c(-2.5, 2.5)) + ylim(c(-2.5, 2.5)) + ggtitle("Trans-Diverging") + theme(legend.position = "none") + theme_classic()

# cis only
plotcis <- ggplot(filter(plotdf, divergence_mode == "cis"), aes(x = log(parental_expr), y = log(hybrid_expr))) + geom_point(color = hue_pal()(5)[1]) + 
  xlab("") + ylab("") +
  xlim(c(-2.5, 2.5)) + ylim(c(-2.5, 2.5)) + ggtitle("Cis-Diverging") + theme(legend.position = "none") + theme_classic()

# compensated only
plotcomp <- ggplot(filter(plotdf, divergence_mode == "compensated"), aes(x = log(parental_expr), y = log(hybrid_expr))) + geom_point(color = hue_pal()(5)[2]) +
  xlab("") + ylab("") +
  xlim(c(-2.5, 2.5)) + ylim(c(-2.5, 2.5)) + ggtitle("Compensated") + theme(legend.position = "none") + theme_classic()

# conserved only
plotcons <- ggplot(filter(plotdf, divergence_mode == "conserved"), aes(x = log(parental_expr), y = log(hybrid_expr))) + geom_point(color = hue_pal()(5)[3]) + 
  xlab("") + ylab("") +
  xlim(c(-2.5, 2.5)) + ylim(c(-2.5, 2.5)) + ggtitle("Conserved") + theme(legend.position = "none") + theme_classic()

# weirdos
plotweird <- ggplot(filter(plotdf, divergence_mode == "weirdos"), aes(x = log(parental_expr), y = log(hybrid_expr))) + geom_point(color = hue_pal()(5)[5]) + 
  xlab("") + ylab("") +
  xlim(c(-2.5, 2.5)) + ylim(c(-2.5, 2.5)) + ggtitle("Weirdos") + theme(legend.position = "none") + theme_classic()

library(ggpubr)
figure <- ggarrange(plotlist = list(plotfull, plotcis, plottrans, plotcons, plotcomp, plotweird), nrow = 2, ncol = 3, common.legend = FALSE)
annotate_figure(figure, left = text_grob("Hybrid Median Expression Ratio (cer/par)", color = "black", face = "bold", rot = 90),
                bottom = text_grob("Parent Median Expression Ratio (cer/par)", color = "black", face = "bold")) + theme(legend.position = "none")

save(spaldf, file = "data_files/cis_trans_df.RData")

# TODO if necessary: Seeing if genes are only ID'd as diverging based on one experiment



