########## Options to split datasets by replicates (randomly or evenly between environments) and to output Environmental-specific datasets for WGCNA ########
sapply(c("tidyr", "dplyr", "readr", "purrr"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Redhuis2024/")

# read in count data and sample info
load("data_files/Cleaned_Redhuis_Data.RData")

################## Environment-specific networks ##############################
# partitioning the data for environmental-specific networks, 7 conditions (4 with replicates) for each of the 3 species (30 networks total):
# 1) TF deletions only (the 46 TF del samples in LowN timepoint 0, so YPD)
# 2) Low Nitrogen only (only WT samples at all timepoints)
# 3) TF del x Low N (all LowN samples)
# 4) Cell Cycle
# 5) OD (Hap4)
# 6) Low Phosphorus

# constructing our environment specific counts and sample info lists (this is the best format for WGCNA to loop through each list)
OrganismNames <- c("cer", "par") # no hyb currently --- just need env-specific datasets for cer vs par
ExperimentNames <- c("LowN", "CC", "HAP4", "LowPi") # Edit these to change which experiment datasets are included
lookup_table <- tibble(name = vector(mode = "character", length = 0),
                       experiment = vector(mode = "character", length = 0),
                       organism = vector(mode = "character", length = 0),
                       time_point_str = vector(mode = "list", length = 0),
                       genotype = vector(mode = "list", length = 0),
                       split_along = vector(mode = "character", length = 0))
for (o in OrganismNames) {
  for (e in ExperimentNames) {
    env_name <- paste(e, o, sep = "_") 
    lookup_table_row <- tibble(name = env_name,
                               experiment = e,
                               organism = o,
                               time_point_str = NA,
                               genotype = NA,
                               split_along = NA)
    if (e == "LowN") {
      lookup_table_row$time_point_str <- list(c("16 h, low N", "1 h, low N", "0 h, YPD"))
      lookup_table_row$genotype <- sample_info %>% filter(experiment == "LowN") %>% select(genotype) %>% pull() %>% unique() %>% list()
      lookup_table_row$split_along <- "genotype" # NOT condition, that'll split up samples from different timepoints but from the same well! Splitting by genotype already guarantees an even distribution of timepoints because every genotype/replicate has all 3 timepoints!
    }
    if (e == "CC") {
      lookup_table_row$time_point_str <- sample_info %>% filter(experiment == "CC") %>% select(time_point_str) %>% pull() %>% unique() %>% list()
      lookup_table_row$genotype <- list(c("WT"))
      lookup_table_row$split_along <- "time_point_num"
    }
    if (e == "HAP4") {
      lookup_table_row$time_point_str <- sample_info %>% filter(experiment == "HAP4") %>% select(time_point_str) %>% pull() %>% unique() %>% list()
      lookup_table_row$genotype <- list(c("WT", "HAP4delete"))
    }
    if (e == "LowPi") {
      lookup_table_row$time_point_str <- sample_info %>% filter(experiment == "LowPi") %>% select(time_point_str) %>% pull() %>% unique() %>% list()
      lookup_table_row$genotype <- list(c("WT"))
    }
    lookup_table <- bind_rows(lookup_table, lookup_table_row)
  }
} 
rm(o, e, lookup_table_row)

# functions that will be useful

#### Subset Count Matrix and Sample Info ####
# @input: count matrix and sample info to subset, boolean vector the same length as ncol(count matrix) and nrow(sample info) instructing which samples to keep
# @output: the subsetted count matrix and sample info
subsetCountMatrixAndSampleInfo <- function(cts, infodf, samplesToKeep) {
  return(list(counts = cts[,samplesToKeep], info = infodf[samplesToKeep,]))
}

# tests for subsetCountMatrixAndSampleInfo
CC_cer <- counts[,sample_info$organism == "cer" & # subsetting the old fashioned way
                      sample_info$experiment == "CC"]
si_CC_cer <- sample_info[sample_info$organism == "cer" &
                              sample_info$experiment == "CC",]
keepers <- (sample_info$organism == "cer") & (sample_info$experiment == "CC")
CC_cer_subsetResults <- subsetCountMatrixAndSampleInfo(counts, sample_info, keepers) # subsetting with my new function
sum(CC_cer == CC_cer_subsetResults$counts)/(nrow(CC_cer)*ncol(CC_cer))
sum(si_CC_cer == CC_cer_subsetResults$info)/(ncol(si_CC_cer)*nrow(si_CC_cer)) # both methods are the same
rm(CC_cer, si_CC_cer, keepers, CC_cer_subsetResults)

#### Split Replicates ####
# split replicates of a given dataset evenly into 2 groups per organism/allele to allow an even distribution of samples in each replicate group based on a specified column of sample_info (usually timepoint or genotype)
# @input: 1) counts matrix and sample info dataframe in ONE SPECIES/EXPERIMENT for samples WITH REPLICATES
#         2) name of the sample info column that you want an even distribution of the values in between both replicates
# @output: list of 4 items: rep1 and rep2 count matrices and sample info data frames
splitReplicates <- function(cts, infodf, splitAlong = colnames(sample_info)) {
  rep1_cts <- matrix(nrow = nrow(cts), ncol = 0)
  rep2_cts <- matrix(nrow = nrow(cts), ncol = 0)
  rep1_info <- matrix(nrow = 0, ncol = ncol(infodf))
  rep2_info <- matrix(nrow = 0, ncol = ncol(infodf))
  condition_list <- infodf[, splitAlong] %>% pull()
  conditionsToBalance <- condition_list %>% unique()
  for (sa in conditionsToBalance) {
    idxs <- which(condition_list == sa)
    replist <- infodf[idxs,] %>% select("well_flask_ID") %>% pull() # replist is in same order as idxs
    reps <- unique(replist) # reps is NOT aligned with idxs
    if (length(reps) %% 2 == 0) { # if there's an even number of replicates
      reps1 <- sample(reps, length(reps)/2)
      reps2 <- setdiff(reps, reps1)
    }
    if (length(reps) %% 2 == 1) { # if there's an odd number of replicates, removing a random one
      to_remove_idx <- sample(c(1:length(reps)), 1)
      cat("Odd number: removing", to_remove_idx, "th sample\n")
      reps <- reps[-to_remove_idx]
      reps1 <- sample(reps, length(reps)/2)
      reps2 <- setdiff(reps, reps1)
    }
    idx1 <- idxs[which(replist %in% reps1)]
    idx2 <- idxs[which(replist %in% reps2)]
    if (!setequal(c(idx1, idx2), idxs[which(replist %in% reps)])) { 
      stop("Indices are not one-to-one represented among the two partitions for condition:", sa, "\n", "idx1:", paste0(idx1, sep = " "), "\n", "idx2:", paste0(idx2, sep = " "), "\n", "idxs:", paste0(idxs, sep = " "), "\n")
    }
    rep1_cts <- cbind(rep1_cts, cts[,idx1, drop = FALSE]) # drop = FALSE preserves column names for vectors
    rep2_cts <- cbind(rep2_cts, cts[,idx2, drop = FALSE])
    rep1_info <- rbind(rep1_info, infodf[idx1,, drop = FALSE])
    rep2_info <- rbind(rep2_info, infodf[idx2,, drop = FALSE])
  }
  return(list(rep1_cts = rep1_cts, rep2_cts = rep2_cts, rep1_info = rep1_info, rep2_info = rep2_info))
}

# tests for splitReplicates
LowN_cer <- subsetCountMatrixAndSampleInfo(counts, sample_info,
                                           sample_info$experiment == "LowN" &
                                             sample_info$organism == "cer")

LowN_cer$info %>% select(sample_name, time_point_str, genotype, well_flask_ID)
# before splitting:
test_table <- table(paste(LowN_cer$info$genotype, LowN_cer$info$well_flask_ID))
test_table
cat("percent replicates with all 3 timepoints:", sum((test_table == 3)/length(test_table)), "\n")
# splitting:
test <- splitReplicates(LowN_cer$counts, LowN_cer$info, "genotype")
# after splitting:
dim(test$rep1_info)
dim(test$rep2_info)
test_table1 <- table(paste(test$rep1_info$genotype, test$rep1_info$well_flask_ID))
test_table1
cat("percent replicates with all 3 timepoints:", sum((test_table1 == 3)/length(test_table1)), "\n")
test_table2 <- table(paste(test$rep2_info$genotype, test$rep2_info$well_flask_ID))
test_table2
cat("percent replicates with all 3 timepoints:", sum((test_table2 == 3)/length(test_table2)), "\n")
common_genotypes <- union(test$rep1_info$genotype, test$rep2_info$genotype)
for (gen in common_genotypes) {
  rep1 <- test$rep1_info %>% filter(genotype == gen) %>% select(well_flask_ID) %>% pull() %>% unique()
  rep2 <- test$rep2_info %>% filter(genotype == gen) %>% select(well_flask_ID) %>% pull() %>% unique()
  cat("Genotype:", gen, "1:", rep1, "2:", rep2, "\n")
  if (length(intersect(rep1, rep2)) != 0) {
    stop("one rep present in both datasets!\n")
  }
}
rm(test, test_table, test_table1, test_table2, gen, LowN_cer, common_genotypes)

# partitions correct samples to each es/si dataset (each name column in lookup table)
# @input: dataset name
# @output: corresponding es and si datasets
partitionExperiments <- function(ex) {
  lt <- filter(lookup_table, name == ex)
  keepers <- sample_info$organism == lt$organism &
    sample_info$experiment == lt$experiment &
    sample_info$time_point_str %in% unlist(lt$time_point_str) &
    sample_info$genotype %in% unlist(lt$genotype)
  output <- subsetCountMatrixAndSampleInfo(counts, sample_info, keepers)
  return(output)
}

# partitions an environment-specific es and si into two es's and si's
# @input: name of environment-specific dataset to partition into two replicates
# @output: list of lists with rep1 and rep2 es's and si's
partitionReplicates <- function(env) {
  sa <- lookup_table %>% filter(name == env) %>% select(split_along) %>% pull()
  idx <- which(names(es) == env)
  output <- splitReplicates(es[[idx]], si[[idx]], sa)
  return(list(cts = list(r1 = output$rep1_cts,
                        r2 = output$rep2_cts),
              info = list(r1 = output$rep1_info,
                        r2 = output$rep2_info)))
}

#### Now to actually partitioning the data by organism/experiment ####
partitioned_datasets <- map(lookup_table$name, partitionExperiments)
es <- lapply(partitioned_datasets, function(x) {return(x[[1]])})
si <- lapply(partitioned_datasets, function(x) {return(x[[2]])})
names(es) <- lookup_table$name
names(si) <- lookup_table$name

# Now adding replicate datasets (skipping in 5/11/23 version of the script)
EnvsWithReplicates <- c("LowN", "CC")
EnvsWithReplicatesBySpecies <- lapply(OrganismNames, function(o) {
  return(paste(EnvsWithReplicates, o, sep = "_"))
  }) %>% unlist()
EnvsWithReplicates
EnvsWithReplicatesBySpecies

# partitioning replicates and adding them to es/si as new list elements
split_rep_results <- map(EnvsWithReplicatesBySpecies, partitionReplicates)
names(split_rep_results) <- EnvsWithReplicatesBySpecies
for(env in EnvsWithReplicatesBySpecies) {
  r1str <- paste(env, "rep1", sep = "_")
  r2str <- paste(env, "rep2", sep = "_")
  es[[r1str]] <- split_rep_results[[env]]$cts$r1
  es[[r2str]] <- split_rep_results[[env]]$cts$r2
  si[[r1str]] <- split_rep_results[[env]]$info$r1
  si[[r2str]] <- split_rep_results[[env]]$info$r2
} 
rm(env)

#### Verifying datasets have been succesfully partitioned ####

# basic check: do the counts and info entries all line up and contain the same samples?
checkSameSamplesInCountsAndInfo <- function(cts, info) {
  if (length(colnames(cts)) == 0) {
    stop("colnames not found")
  }
  pct <- sum(colnames(cts) == info$sample_name)/ncol(cts)
  cat("percent sample names the same:", pct, "\n")
  if (pct < 1) {
    stop ("some samples are different!\n")
  }
}
dummy <- map2(es, si, checkSameSamplesInCountsAndInfo) # all good!
rm(dummy)

# For datasets with replicates, making sure every sample is represented
# once and only once between both of the replicate datasets
checkRepsAreCompletePartition <- function(env) {
  # for es
  full <- es[[env]]
  r1 <- es[[paste0(env, "_rep1")]]
  r2 <- es[[paste0(env, "_rep2")]]
  nduplicated <- cbind(full, r1, r2) %>% t() %>% duplicated() %>% sum()
  if (nduplicated != ncol(cbind(r1, r2))) {
    stop (env, " has", nduplicated, " duplicated count samples when it should have", ncol(cbind(r1, r2)), "\n")
  }
  if (sum(duplicated(t(full))) + sum(duplicated(t(r1))) + sum(duplicated(t(r2))) != 0) {
    stop (env, " has duplicated info samples in individual datasets\n")
  }
  
  # for si
  full <- si[[env]]
  r1 <- si[[paste0(env, "_rep1")]]
  r2 <- si[[paste0(env, "_rep2")]]
  nduplicated <- rbind(full, r1, r2) %>% duplicated() %>% sum()
  if (nduplicated != nrow(rbind(r1, r2))) {
    stop (env, " has", nduplicated, " duplicated info samples when it should have", nrow(rbind(r1, r2)), "\n")
  }
  if (sum(duplicated(full)) + sum(duplicated(r1)) + sum(duplicated(r2)) != 0) {
    stop (env, " has duplicated count samples in individual datasets\n")
  }
  return(TRUE)
}
map(EnvsWithReplicatesBySpecies, checkRepsAreCompletePartition) %>% unlist() %>% all()

# most complicated one: make sure that the two replicates have about even
# distribution of samples of each value in the corresponding split_along for that experiment
checkRepsRepresentSplitAlongRelativelyEvenly <- function(env) {
  r1 <- paste0(env, "_rep1")
  r2 <- paste0(env, "_rep2")
  sa <- lookup_table %>% filter(name == env) %>% select(split_along) %>% pull()
  tab1 <- si[[r1]] %>% select(all_of(sa)) %>% table() %>% sort()
  tab2 <- si[[r2]] %>% select(all_of(sa)) %>% table() %>% sort()
  if (length(setdiff(names(tab1), names(tab2))) > 0) {
    cat("categories not found in both replicates:", setdiff(names(tab1), names(tab2)), "\n")
  }
  common_names <- intersect(names(tab1), names(tab2))
  tab1 <- tab1[common_names]
  tab2 <- tab2[common_names]
  tabdiff <- abs(tab1 - tab2)
  if (sum(tabdiff > 1) > 0) {
    cat("reps not balanced:", names(tab1)[which(tabdiff > 3)], "\n")
  }
  else {print("All reps are evenly balanced")}
  return(list(r1 = tab1, r2 = tab2, diff = tabdiff))
}
rep_tables <- map(EnvsWithReplicatesBySpecies, checkRepsRepresentSplitAlongRelativelyEvenly) %>% unlist() %>% all()

# same as above except instead of computationally checking for uneven replicates, 
# it just shows you the table to judge for yourself
inspectSplitAlongDistribution <- function(r1, r2, splitAlong = colnames(sample_info)) {
  split1 <- r1 %>% select(all_of(splitAlong)) %>% pull() %>% paste("rep1", sep = "_")
  split2 <- r2 %>% select(all_of(splitAlong)) %>% pull() %>% paste("rep2", sep = "_")
  rep_table <- table(c(split1, split2))
  return(rep_table[order(names(rep_table))])
}
# examples:
# TFdelxLowN hyb
inspectSplitAlongDistribution(si$LowN_cer_rep1, si$LowN_cer_rep2, splitAlong = "condition")
# HAP2 par
inspectSplitAlongDistribution(si$CC_par_rep1, si$CC_par_rep2, splitAlong = "time_point_str")
# LowN cer
inspectSplitAlongDistribution(si$LowN_par_rep1, si$LowN_par_rep2, splitAlong = "time_point_str")

# making sure each (not checking replicates b/c we already verified they're a good partition)
# dataset only contains variable values that are expected to be contained according to the lookup table
checkDatasetsOnlyHaveCorrectConditions <- function(.nm) {
  gen <- lookup_table %>% filter(name == .nm) %>% select(genotype) %>% pull() %>% unlist()
  ex <- lookup_table %>% filter(name == .nm) %>% select(experiment) %>% pull()
  timept <- lookup_table %>% filter(name == .nm) %>% select(time_point_str) %>% pull() %>% unlist()
  o <- lookup_table %>% filter(name == .nm) %>% select(organism) %>% pull()
  
  isCorrectGen <- si[[.nm]] %>% select(genotype) %>% pull() %in% gen %>% all()
  isCorrectExp <- si[[.nm]] %>% select(experiment) %>% pull() %in% ex %>% all()
  isCorrectTimepoint <- si[[.nm]] %>% select(time_point_str) %>% pull() %in% timept %>% all()
  isCorrectOrganism <- si[[.nm]] %>% select(organism) %>% pull() %in% o %>% all()
  
  if (all(isCorrectGen, isCorrectExp, isCorrectTimepoint, isCorrectOrganism)) {
    return(TRUE)
  }
  else {
    stop ("One or more variables had abberant values present\n", 
          "Genotypes correct:", isCorrectGen, "\n",
          "Experiment correct:", isCorrectExp, "\n",
          "Timepoint correct:", isCorrectTimepoint, "\n",
          "Organism correct:", isCorrectOrganism, "\n")
  }
}
map_lgl(lookup_table$name, checkDatasetsOnlyHaveCorrectConditions) %>% all()

#### Saving Environment-Specific datasets ####
save(es, si, file = "data_files/Cleaned_Redhuis_Data_Env_Specific.RData")

#### Now repeating with hybrid allele-specific datasets hyc and hyp ####
# allele-specific hybrid samples
load("data_files/Cleaned_Redhuis_Data_AlleleSpecific.RData")

# reconstructing lookup table for hyc/hyp
OrganismNames <- c("hyc", "hyp")
ExperimentNames <- c("LowN", "CC", "HAP4", "LowPi") # Edit these to change which experiment datasets are included
lookup_table <- tibble(name = vector(mode = "character", length = 0),
                       experiment = vector(mode = "character", length = 0),
                       organism = vector(mode = "character", length = 0),
                       time_point_str = vector(mode = "list", length = 0),
                       genotype = vector(mode = "list", length = 0),
                       split_along = vector(mode = "character", length = 0))
for (o in OrganismNames) {
  for (e in ExperimentNames) {
    env_name <- paste(e, o, sep = "_") 
    lookup_table_row <- tibble(name = env_name,
                               experiment = e,
                               organism = o,
                               time_point_str = NA,
                               genotype = NA,
                               split_along = NA)
    if (e == "LowN") {
      lookup_table_row$time_point_str <- list(c("16 h, low N", "1 h, low N", "0 h, YPD"))
      lookup_table_row$genotype <- list(c("WT"))
      lookup_table_row$split_along <- "time_point_num"
    }
    if (e == "CC") {
      lookup_table_row$time_point_str <- sample_info_allele %>% filter(experiment == "CC") %>% select(time_point_str) %>% pull() %>% unique() %>% list()
      lookup_table_row$genotype <- list(c("WT"))
      lookup_table_row$split_along <- "time_point_num"
    }
    if (e == "HAP4") {
      lookup_table_row$time_point_str <- sample_info_allele %>% filter(experiment == "HAP4") %>% select(time_point_str) %>% pull() %>% unique() %>% list()
      lookup_table_row$genotype <- list(c("WT", "HAP4delete"))
    }
    if (e == "LowPi") {
      lookup_table_row$time_point_str <- sample_info_allele %>% filter(experiment == "LowPi") %>% select(time_point_str) %>% pull() %>% unique() %>% list()
      lookup_table_row$genotype <- list(c("WT"))
    }
    lookup_table <- bind_rows(lookup_table, lookup_table_row)
  }
} 
rm(o, e, lookup_table_row)

partitionExperimentsAllele <- function(ex) {
  lt <- filter(lookup_table, name == ex)
  hyc_or_hyp <- if_else(sample_info_allele$allele == "cer", true = "hyc", false = "hyp")
  keepers <- hyc_or_hyp == lt$organism &
    sample_info_allele$experiment == lt$experiment &
    sample_info_allele$time_point_str %in% unlist(lt$time_point_str) &
    sample_info_allele$genotype %in% unlist(lt$genotype)
  output <- subsetCountMatrixAndSampleInfo(counts_allele, sample_info_allele, keepers)
  return(output)
}

# partitioning
partitioned_datasets <- map(lookup_table$name, partitionExperimentsAllele)
ales <- lapply(partitioned_datasets, function(x) {return(x[[1]])})
alsi <- lapply(partitioned_datasets, function(x) {return(x[[2]])})
names(ales) <- lookup_table$name
names(alsi) <- lookup_table$name

# QC
dummy <- map2(ales, alsi, checkSameSamplesInCountsAndInfo) # all good!
rm(dummy)

checkDatasetsOnlyHaveCorrectConditionsAllele <- function(.nm) {
  gen <- lookup_table %>% filter(name == .nm) %>% select(genotype) %>% pull() %>% unlist()
  ex <- lookup_table %>% filter(name == .nm) %>% select(experiment) %>% pull()
  timept <- lookup_table %>% filter(name == .nm) %>% select(time_point_str) %>% pull() %>% unlist()
  o <- lookup_table %>% filter(name == .nm) %>% select(organism) %>% pull()
  
  isCorrectGen <- alsi[[.nm]] %>% select(genotype) %>% pull() %in% gen %>% all()
  isCorrectClass <- alsi[[.nm]] %>% select(experiment) %>% pull() %in% ex %>% all()
  isCorrectTimepoint <- alsi[[.nm]] %>% select(time_point_str) %>% pull() %in% timept %>% all()
  org <- alsi[[.nm]] %>% select(allele) %>% pull() 
  isCorrectOrganism <- if_else(org == "cer", true = "hyc", false = "hyp") %in% o %>% all()
  
  if (all(isCorrectGen, isCorrectClass, isCorrectTimepoint, isCorrectOrganism)) {
    return(TRUE)
  }
  else {
    stop ("One or more variables had abberant values present\n", 
          "Genotypes correct:", isCorrectGen, "\n",
          "Experiment correct:", isCorrectClass, "\n",
          "Timepoint correct:", isCorrectTimepoint, "\n",
          "Organism correct:", isCorrectOrganism, "\n")
  }
}
map_lgl(lookup_table$name, checkDatasetsOnlyHaveCorrectConditionsAllele) %>% all()

# saving 
save(ales, alsi, file = "data_files/Cleaned_Redhuis_Data_Env_Specific_Allele.RData")

################### Hybrid Split Networks for WGCNA - Archived for now ###############################
# Archived 3/27/23 because we're just working with summed hybrid expression for now - hyc + hyp
# But this may be useful for WGCNA once we figure out how to rationalize an "allele-specific network" inside a cell with a second, apparently independent "allele-specific network"

# Random hybrid split
# Creating "hybrid split" control data - 50-50 split of hybrid data so that no 
# allele's measurement was in the same environment as the other allele's measurement
# (as a way to decouple common cell environment from common trans regulatory environment)
# (also downsampling parents by 50% to make it even)
# downsampling parents
cer_idx <- which(sample_info$organism == "cer")
cer_sample_idx <- sample(cer_idx, length(cer_idx)/2)

par_idx <- which(sample_info$organism == "par")
par_sample_idx <- sample(par_idx, length(par_idx)/2)

# downsampling one allele (cerevisiae arbitrarily) then choosing the opposite samples
# for the second allele
hyc_idx <- which(sample_info$organism == "hyc")
hyc_sample_idx <- sample(hyc_idx, length(hyc_idx)/2)
hyc_not_sample <- setdiff(hyc_idx, hyc_sample_idx)
# conveniently, the corresponding hyp sample is directly below its hyc sample (of course they are---you literally coded them in the consolidateSpeciesCounts loop)
hyb_names <- gsub("hy.", "hyb", sample_info$sample_name[sample_info$organism == "hyc" | sample_info$organism == "hyp"])
cat("percent hyb names with exactly two entries:", sum(table(hyb_names) == 2)/length(table(hyb_names)))
cat("percent hyp samples directly below their hyc samples:", 
    sum(sample_info$sample_name[hyc_not_sample] == 
          gsub("hyp", "hyc", sample_info$sample_name[hyc_not_sample + 1]))/length(hyc_not_sample))

hyp_sample_idx <- (hyc_not_sample + 1) # long walk for a short drink of water

# actually subsetting data
hybsplit_idx <- c(cer_sample_idx, par_sample_idx, hyc_sample_idx, hyp_sample_idx)
counts_hybsplit <- counts[,hybsplit_idx]
sample_info_hybsplit <- sample_info[hybsplit_idx,]

# saving to .RData file
save(counts_hybsplit, sample_info_hybsplit, file = "data_files/Cleaned_Barkai_Data_HybSplit.RData")

# Hybrid split by replicates
# Filtering for two replicates per condition (condition = genotype + experiment + timepoint) then randomly splitting replicates
# also need the same set of conditions for all 4 species

# getting our condition set
cer_conditions <- sample_info[sample_info$organism == "cer",]$condition
par_conditions <- sample_info[sample_info$organism == "par",]$condition
hyc_conditions <- sample_info[sample_info$organism == "hyc",]$condition
hyp_conditions <- sample_info[sample_info$organism == "hyp",]$condition

all_conditions <- unique(sample_info$condition)
conditions_inAll4 <- reduce(list(cer_conditions, par_conditions, hyc_conditions, hyp_conditions), intersect)

# getting lists of conditions that don't have at least 2 replicates in each species
cer_yes_reps <- names(which(table(cer_conditions) > 1))
par_yes_reps <- names(which(table(par_conditions) > 1))
hyc_yes_reps <- names(which(table(hyc_conditions) > 1))
hyp_yes_reps <- names(which(table(hyp_conditions) > 1))

condition_set <- reduce(list(conditions_inAll4, cer_yes_reps, par_yes_reps, hyp_yes_reps), intersect)

# Creating hybrid split by spliting REPLICATES 50-50 (rep1 to cer, rep2 to par or vice vesrsa for each sample)
# Only 2 replicates are kept for each condition (so samples w/o replicates or w/ excess replicates are thrown out)
cat("percent of sample conditions that will be removed: ", length(condition_set)/length(all_conditions)) # not only that but HAP4 and LowPi experiments are totally gone b/c they have no replicates

table(sample_info[sample_info$condition %in% condition_set,]$condition) # all should have at least 8

# looping through condition set to put one sample from each condition into a new count matrix (if you have a way to do this that doesn't involve 5,000 lines, I'm all ears)
cer_counts_rep1 <- data.frame(matrix(nrow = nrow(counts), ncol = 0))
cer_counts_rep2 <- data.frame(matrix(nrow = nrow(counts), ncol = 0))
par_counts_rep1 <- data.frame(matrix(nrow = nrow(counts), ncol = 0))
par_counts_rep2 <- data.frame(matrix(nrow = nrow(counts), ncol = 0))
hyc_counts_rep1 <- data.frame(matrix(nrow = nrow(counts), ncol = 0))
hyc_counts_rep2 <- data.frame(matrix(nrow = nrow(counts), ncol = 0))
hyp_counts_rep1 <- data.frame(matrix(nrow = nrow(counts), ncol = 0))
hyp_counts_rep2 <- data.frame(matrix(nrow = nrow(counts), ncol = 0))

for (condition in condition_set) {
  idx1 <- sample(c(1,2), 1)
  idx2 <- setdiff(c(1,2), idx1)
  # cer
  cer_idx1 <- which(sample_info$condition == condition & sample_info$organism == "cer")[idx1]
  cer_idx2 <- which(sample_info$condition == condition & sample_info$organism == "cer")[idx2]
  cer_counts_rep1 <- cbind(cer_counts_rep1, counts[, cer_idx1])
  cer_counts_rep2 <- cbind(cer_counts_rep2, counts[, cer_idx2])
  colnames(cer_counts_rep1)[ncol(cer_counts_rep1)] <- sample_info$sample_name[cer_idx1]
  colnames(cer_counts_rep2)[ncol(cer_counts_rep2)] <- sample_info$sample_name[cer_idx2]
  # par
  par_idx1 <- which(sample_info$condition == condition & sample_info$organism == "par")[idx1]
  par_idx2 <- which(sample_info$condition == condition & sample_info$organism == "par")[idx2]
  par_counts_rep1 <- cbind(par_counts_rep1, counts[, par_idx1])
  par_counts_rep2 <- cbind(par_counts_rep2, counts[, par_idx2])
  colnames(par_counts_rep1)[ncol(par_counts_rep1)] <- sample_info$sample_name[par_idx1]
  colnames(par_counts_rep2)[ncol(par_counts_rep2)] <- sample_info$sample_name[par_idx2]
  # hyc
  hyc_idx1 <- which(sample_info$condition == condition & sample_info$organism == "hyc")[idx1] # Note: hyc_counts_rep1 and hyp_counts_rep1 are FROM THE SAME CELL (same for the rep2s)
  hyc_idx2 <- which(sample_info$condition == condition & sample_info$organism == "hyc")[idx2]
  hyc_counts_rep1 <- cbind(hyc_counts_rep1, counts[, hyc_idx1])
  hyc_counts_rep2 <- cbind(hyc_counts_rep2, counts[, hyc_idx2])
  colnames(hyc_counts_rep1)[ncol(hyc_counts_rep1)] <- sample_info$sample_name[hyc_idx1]
  colnames(hyc_counts_rep2)[ncol(hyc_counts_rep2)] <- sample_info$sample_name[hyc_idx2]
  # hyp
  hyp_idx1 <- which(sample_info$condition == condition & sample_info$organism == "hyp")[idx1]
  hyp_idx2 <- which(sample_info$condition == condition & sample_info$organism == "hyp")[idx2]
  hyp_counts_rep1 <- cbind(hyp_counts_rep1, counts[, hyp_idx1])
  hyp_counts_rep2 <- cbind(hyp_counts_rep2, counts[, hyp_idx2])
  colnames(hyp_counts_rep1)[ncol(hyp_counts_rep1)] <- sample_info$sample_name[hyp_idx1]
  colnames(hyp_counts_rep2)[ncol(hyp_counts_rep2)] <- sample_info$sample_name[hyp_idx2]
}

# fianlly we can save this insanity
save(cer_counts_rep1, cer_counts_rep2, 
     par_counts_rep1, par_counts_rep2, 
     hyc_counts_rep1, hyc_counts_rep2, 
     hyp_counts_rep1, hyp_counts_rep2,
     sample_info, file = "data_files/Cleaned_Barkai_Data_RepSplit.RData") # Note: hyc_counts_rep1 and hyp_counts_rep1 are FROM THE SAME CELL (same for the rep2s)

############## Proceed if you dare: Ensuring the same set of sample conditions across species ###############

# NB: This is shaping up to be a stressful process that might not even be necessary
# Make sure you NEED to have the same samples across species (like for plotting them on a PC plot or something) before continuing

# checking if the same number of each timepoint is present between species
tabc <- table(si$cer$TFdel[,c("collection_date", "time_point_str")])
tabp <- table(si$par$TFdel[,c("collection_date", "time_point_str")])
tabh <- table(si$hyb$TFdel[,c("collection_date", "time_point_str")])
# it turns out all but one collection date is shared by all three species, so let's remove that date
multisetIntersect <- function(x, y, ...) {
  sets <- list(x, y, ...)
  output <- Reduce(f = intersect, x = sets)
  return(output)
}
goodDates <- multisetIntersect(si$cer$TFdel$collection_date, 
                               si$par$TFdel$collection_date, 
                               si$hyb$TFdel$collection_date)
# function to streamline removing samples from both es and si based on their values in a column of si
# @input: 
# 1) env-specific counts list of lists es
# 2) corresponding sample information si
# 3) environment of interest
# 4) which column of si to check for good values
# 5) vector of all acceptable values for samples to have in the si column of interest
# @output: new es and si with the offending samples removed
removeSamplesBySIColumn <- function(es_df = es, si_tab = si, env = names(es$cer), 
                                    si_col = colnames(sample_info), goodSamps) {
  es_output <- map2(es_df, si_tab, function(x, y) {
    keep <- y[[env]][,si_col] %in% goodSamps
    x[[env]] <- x[[env]][,keep]
    y[[env]] <- y[[env]][keep,]
    return(x)
  })
  si_output <- lapply(si_tab, function(y) {
    keep <- y[[env]][,si_col] %in% goodSamps
    y[[env]] <- y[[env]][keep,]
    return(y)
  })
  return(list(es_output, si_output))
}
new_es_si <- removeSamplesBySIColumn(es, si, "TFdel", "collection_date", goodDates)
es <- new_es_si[[1]]
si <- new_es_si[[2]]

tabc <- table(si$cer$TFdel[,c("collection_date", "time_point_str")])
tabp <- table(si$par$TFdel[,c("collection_date", "time_point_str")])
tabh <- table(si$hyb$TFdel[,c("collection_date", "time_point_str")])
tab <- matrix(nrow = nrow(tabc), ncol = ncol(tabc))
for (i in 1:length(tab)) {
  tab[i] <- min(c(tabc[i], tabp[i], tabh[i]))
} 
rm(i)

# Where I left off: I think I want to sort the samples based on 1) collection date, 2) timepoint then narrow it down to the same amount of samples per collection date x timepoint in all 3 species
# Sanity check: Each rep# should only have one hit for each timepoint