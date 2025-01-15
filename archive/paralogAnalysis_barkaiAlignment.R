rm(list=ls())
sapply(c("tidyverse", "edgeR", "openxlsx", "Glimma"), require, character.only=TRUE)
setwd("Documents/Wittkopp_Lab/paralogs/barkailab_paper")

# following this Tag-seq sript: https://github.com/ben-laufer/Tag-seq/blob/master/tag-seq.R

#### Table of Contents ####
# 1. Count Data Import + Normalization
# 2. Functions
#    I. Regulatory Correlation
#    II. Paralog Expression Ratios
#    III. dN/dS
# 3. Objects
#    I. Regulatory Correlation
#    II. Paralog Expression Ratios
#    III. dN/dS
# 4. Pretty Plots

##################### 1. Barkai Count Data Import + Normalization #################
counts_strict_cer <- read.xlsx("Supplemental_Table_S11.xlsx", 
                               sheet = 2, na.strings = "")
geneNames <- counts_strict_cer[,1]
counts_strict_cer <- counts_strict_cer[,2:ncol(counts_strict_cer)] %>% 
  magrittr::set_rownames(geneNames)
counts_strict_cer[is.na(counts_strict_cer)] <- 0

counts_strict_par <- read.xlsx("Supplemental_Table_S12.xlsx", 
                               sheet = 2, na.strings = "")
geneNames_par <- counts_strict_par[,1]
sum(geneNames == geneNames_par) #just checking they're in the same order
counts_strict_par <- counts_strict_par[,2:ncol(counts_strict_par)] %>% 
  magrittr::set_rownames(geneNames_par)
counts_strict_par[is.na(counts_strict_par)] <- 0

sample_info <- read.xlsx("../bioSample1to999.xlsx", na.strings="not applicable") %>% 
  bind_rows(read.xlsx("../bioSample1000toEnd.xlsx", na.strings="not applicable"))
desired_sample_order <- sapply(c(1:nrow(sample_info)), function(i) {
  return(which(sample_info$sample_name==colnames(counts_strict_cer)[i]))
  })
sample_info <- sample_info[desired_sample_order,]
sum(sample_info$sample_name == colnames(counts_strict_cer)) # just want to make sure the 1636 conditions are in the same order in the sample info and in the counts table

dgelist0_cer <- DGEList(counts_strict_cer, samples=sample_info) %>% 
  calcNormFactors()
dgelist0_par <- DGEList(counts_strict_par, samples=sample_info) %>% 
  calcNormFactors()

# filter out samples from the wrong species (cer in par reads, par in cer reads)
dgelist1_cer <- dgelist0_cer[, dgelist0_cer$samples$organism != "Saccharomyces paradoxus"]
dgelist1_par <- dgelist0_par[, dgelist0_par$samples$organism != "Saccharomyces cerevisiae"]

# filter out samples with less than a set number of total reads 
# (need to do more research on what an appropriate cutoff is)
libCutoff <- 50000
dgelist2_cer <- dgelist1_cer[,dgelist1_cer$samples$lib.size > libCutoff]
dgelist2_par <- dgelist1_par[,dgelist1_par$samples$lib.size > libCutoff]

par(mfrow = c(2,1), mar = c(1,1,1,1))
hist(dgelist2_cer$samples$lib.size, main = NULL, breaks = 13)
hist(dgelist2_par$samples$lib.size, main = NULL, breaks = 13)
dev.off()

# filter lowly expressed genes (< 1 cpm)
dgelist2_cer$samples$group <- as.factor(dgelist2_cer$samples$experiment)
keep_cer <- filterByExpr(dgelist2_cer, group = dgelist2_cer$samples$group, 
                     lib.size = dgelist2_cer$samples$lib.size)
dgelist2_par$samples$group <- as.factor(dgelist2_par$samples$experiment)
keep_par <- filterByExpr(dgelist2_par, group = dgelist2_par$samples$group, 
                         lib.size = dgelist2_par$samples$lib.size)
keep <- keep_cer | keep_par # keep all genes that at least one species has expr for

dgelist_cer <- dgelist2_cer[keep,,]
dgelist_par <- dgelist2_par[keep,,]

cpm_cer <- cpm(dgelist_cer)
lcpm_cer <- cpm(dgelist_cer, log=TRUE) # used for visualization

cpm_par <- cpm(dgelist_par)
lcpm_par <- cpm(dgelist_par, log=TRUE) # used for visualization

# checking which genes were filtered out
removed <- dgelist2_cer[!keep,,]$counts
large <- removed[order(removed, decreasing=TRUE)][1:20]
large_rows <- numeric(0)
for (i in 1:20) {
       large_rows <- c(large_rows, which(removed == large[i], arr.ind = TRUE)[1])
}
rownames(removed)[large_rows]
hist(dgelist2_cer$counts["YJL038C", 2:ncol(dgelist2_cer$counts)])

# remvoe unnecessary variables to get rid of clutter
rm(list = setdiff(ls(), c("dgelist_cer", "dgelist_par", "cpm_cer", "cpm_par", 
                          "lcpm_cer", "lcpm_par")))

# Barkai Rc and dataset control values
Rc <- read.xlsx("../../../barkailab_paper/Supplemental_Table_S7.xlsx", sheet = 2, na.strings = "")
dataset_control <- read.xlsx("../../../barkailab_paper/Supplemental_Table_S7.xlsx", sheet = 3, na.strings = "")
diverging_genes <- Rc$ORF[Rc$bw_cer_par < 0.2 & dataset_control$cer > 0.5 & dataset_control$par > 0.5]
diverging_genes <- diverging_genes[!is.na(diverging_genes)]
duplicates <- c(marchant$P1[marchant$Duplication != "S"], marchant$P2[marchant$Duplication != "S"])
duplicates <- duplicates[!duplicated(duplicates)]
n_diverging_dupes <- sum(diverging_genes %in% duplicates)
diverging_dupes_enrichment_counts <- matrix(nrow = 2, c(n_diverging_dupes, length(diverging_genes)-n_diverging_dupes, length(duplicates), nrow(Rc)-length(duplicates)))
fisher.test(diverging_dupes_enrichment_counts)
duplicate_enrichment <- tibble("fraction_of_genes" = c(162/(162+119), 3514/(3514+3136)),
                               "group" = c("diverging", "genome"))

# just WGDs
wgds <- c(marchant$P1[marchant$Duplication == "WGD"], marchant$P2[marchant$Duplication == "WGD"])
sum(duplicated(wgds)) # should be 0
n_diverging_wgds <- sum(diverging_genes %in% wgds)
diverging_wgds_enrichment_counts <- matrix(nrow = 2, c(n_diverging_wgds, length(diverging_genes)-n_diverging_wgds, length(wgds), nrow(Rc)-length(wgds)))
fisher.test(diverging_wgds_enrichment_counts)

# just SSDs
ssds <- c(marchant$P1[marchant$Duplication == "SSD"], marchant$P2[marchant$Duplication == "SSD"])
ssds <- ssds[!duplicated(ssds)]
n_diverging_ssds <- sum(diverging_genes %in% ssds)
diverging_ssds_enrichment_counts <- matrix(nrow = 2, c(n_diverging_ssds, length(diverging_genes)-n_diverging_ssds, length(ssds), nrow(Rc)-length(ssds)))
fisher.test(diverging_ssds_enrichment_counts)

###################### 2. Functions ########################

##### I. Regulatory Correlation ####
groupNames <- c("CellCycle", "HAP4andWTonYPD", "SCtoLowPi", "YPD to Low N")

# function that returns the correlation of each gene's expression in an
# experimental group with every other gene in the genome
# input: dgelist d containing single-species counts (cer or par), index of 
# gene of interest, name of experimental group of interest
# output: vector of correlations with all other genes
correlate_with_genome <- function(d, gene_idx, group) {
  counts <- d$counts[-gene_idx, d$samples$group == group]
  gene <- d$counts[gene_idx, d$samples$group == group]
  corrs <- apply(counts, 1, function(x) {
    return(cor(x, gene, method="pearson"))
  })
  return(corrs)
}

# uses correlate_with_genome to extract the median correlation of every 
# gene with every other gene in the genome out of the 4 experiments
correlate_across_experiments <- function(d, gene_idx, groups) {
  corrs_per_expr <- sapply(groups, function(g) {
    return(correlate_with_genome(d, gene_idx, g))
  })
  corrs_across_expr <- apply(corrs_per_expr, 1, median)
}

#### II. Paralog Expression Ratio ####
# Goal: ID pairs whose expression have diverged in tandem between cer + par
# (P1/P2 expr ratio differs btwn cer and par for at least some environments)

# given the SGD names of a pair of paralogs and a gene count matrix (where row names are sgd names),
# returns the expression level ratio of that pair in each column of the count matrix
get_expr_ratio <- function(p1, p2, counts) {
  return(counts[p1,]/counts[p2,])
}

library(INLA)

# function to get p-values for a negative binomial model for paralog pair expression ratios
# a pair is significant if its ratio differs significantly in cer versus par
# input: dgelist and cpm-normalized counts (ROUNDED TO INTEGERS) for each species 
# (containing ONLY samples that you want to include in the model),
# vector of length 2 of SGD IDs for the pair of interest,
# string designating which of the 4 experiments
get_expr_ratio_model <- function(samples_cer, counts_cer, samples_par, counts_par, pairIDs, experimentGroup=
                                   c("CellCycle", "HAP4andWTonYPD", 
                                     "SCtoLowPi", "YPD to Low N" )) {
  
  cat("current pair: ", pairIDs, "(", which(marchant_justparalogs$P1 == pairIDs[1]), "/", nrow(marchant_justparalogs), ")", "\n")
  
  if(!(pairIDs[1] %in% rownames(counts_cer)) | 
     !(pairIDs[2] %in% rownames(counts_cer))) {
    return(tibble("est" = NA, 
                  "p_value" = NA))
  }
  
  #cer
  P1_cer <- counts_cer[pairIDs[1], samples_cer$experiment == experimentGroup]
  P2_cer <- counts_cer[pairIDs[2], samples_cer$experiment == experimentGroup]
  cer_data <- matrix(nrow = length(P1_cer), ncol = 4) %>% as.data.frame()
  cer_data[,1] <- "cer"
  cer_data[,2] <- round(P1_cer)
  cer_data[,3] <- round(P2_cer)
  cer_data[,4] <- samples_cer$time_point[samples_cer$experiment == experimentGroup] %>% parse_number()
  names(cer_data) <- c("species", "P1", "P2", "timepoint")
  
  #par
  P1_par <- counts_par[pairIDs[1], samples_par$experiment == experimentGroup]
  P2_par <- counts_par[pairIDs[2], samples_par$experiment == experimentGroup]
  par_data <- matrix(nrow = length(P1_par), ncol = 4) %>% as.data.frame()
  par_data[,1] <- "par"
  par_data[,2] <- round(P1_par)
  par_data[,3] <- round(P2_par)
  par_data[,4] <- samples_par$time_point[samples_par$experiment == experimentGroup] %>% parse_number()
  names(par_data) <- c("species", "P1", "P2", "timepoint")
  
  d <- rbind(cer_data, par_data)
  d$total_reads <- d$P1 + d$P2
  
  mod <- inla(P1 ~ 1 + species + timepoint, data = d, family = "binomial", Ntrials = total_reads)
  
  coef <- mod$summary.fixed
  fixed_effect_posterior <- mod$marginals.fixed[[2]]
  lower_p <- inla.pmarginal(0, fixed_effect_posterior)
  upper_p <- 1 - inla.pmarginal(0, fixed_effect_posterior)
  post_pred_p <- 2 * (min(lower_p, upper_p))
  
  mod_output <- tibble("est" = coef[2,1], 
                       "p_value" = post_pred_p)
  
  return(mod_output)
}

# version of get_expr_ratio_model that applies to all paralog pairs in marchant
get_expr_ratio_model_allGenes <- function(samples_cer, counts_cer, samples_par, counts_par, experimentGroup=
                                            c("CellCycle", "HAP4andWTonYPD", 
                                              "SCtoLowPi", "YPD to Low N" )) {
  output <- do.call(rbind, mclapply(c(1:nrow(marchant_justparalogs)), function(idx) {
    get_expr_ratio_model(samples_cer, 
                         counts_cer, 
                         samples_par, 
                         counts_par, 
                         c(marchant_justparalogs$P1[idx], marchant_justparalogs$P2[idx]),
                         experimentGroup)
  }, mc.cores = 1))
  return(output)
}

#### III. Calculate K_s and K_a between pairs of paralogs ####

setwd("../../awesome_power/")
library(readxl)

get_cds <- function(sgd, species = c("Scer", "Spar")) {
  filename <- list.files("coding/", pattern = paste0("^.*", sgd, 
                                                     ".*\\Q.\\Ecodon\\Q.\\Emfa"))
  if(length(filename) == 0) {
    cat(sgd, "file not found\n")
    return()
  }
  orthogroup <- read_lines(paste0("coding/", filename))
  cds_start <- grep(paste0(">", species), orthogroup) + 1
  if(length(cds_start) == 0){
    cat(sgd, "doesn't have", species, "for some reason\n")
    return()
  }
  cds_end <- grep(">", orthogroup[cds_start:length(orthogroup)])[1] + cds_start - 2
  if(length(cds_end) == 0 | is.na(cds_end)) {
    cds_end <- length(orthogroup)
  }
  cds <- paste0(orthogroup[cds_start:cds_end], collapse="")
  # ANNA: if you're gonna gsub out the "-" characters anyway, why use condons.mfa (the ortholog alignments)
  # instead of just the .fsa files?!? I'm assuming you just didn't realize that's what it was
  cds <- gsub("-", "", cds)
  return(cds)
}

codons <- list("F" = c("TTT", "TTC"),
               "L" = c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
               "I" = c("ATT", "ATC", "ATA"),
               "M" = "ATG",
               "V" = c("GTT", "GTC", "GTA", "GTG"),
               "S" = c("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
               "P" = c("CCT", "CCC", "CCA", "CCG"),
               "T" = c("ACT", "ACC", "ACA", "ACG"),
               "A" = c("GCT", "GCC", "GCA", "GCG"),
               "Y" = c("TAT", "TAC"),
               "Stop" = c("taa", "tag", "tga"),
               "H" = c("CAT", "CAC"),
               "Q" = c("CAA", "CAG"),
               "N" = c("AAT", "AAC"),
               "K" = c("AAA", "AAG"),
               "D" = c("GAT", "GAC"),
               "E" = c("GAA", "GAG"),
               "C" = c("TGT", "TGC"),
               "W" = "TGG",
               "R" = c("AGA", "AGG", "CGT", "CGC", "CGA", "CGG"),
               "G" = c("GGT", "GGC", "GGA", "GGG"))

# sanity check that I did the codon list correctly - makes sure every 3-letter
# combination of A, C, T, G appears exactly once
for (first in c("A", "T", "C", "G")) {
  for (second in c("A", "T", "C", "G")) {
    for (third in c("A", "T", "C", "G")) {
      query <- paste0(first, second, third)
      times <- length(which((unlist(codons) == query) | 
                              (unlist(codons) == tolower(query))))
      cat(query, "appears", times, "times\n")
    }
  }
}

# Helper function for calculate_ks
# checks if the codons are the same DNA sequence and if not calculates how many
# differences they have (synonomous or nonsynonymous)
# input: two codons to compare
# output: number of differences
check_codon <- function(cdn1, cdn2) {
  if(cdn1 == cdn2) {
    return(0)
  }
  else {
    differences <- 0
    for (i in 1:3) {
      if(substring(cdn1, i, i) != substring(cdn2, i, i)) {
        differences <- differences + 1
      }
    }
    return(differences)
  }
}

# second helper function for calculate_ks
# finds 3-letter name of amino acid given 3-base (DNA) sequence
get_codon_name <- function(codon) {
  for (i in 1:length(codons)) {
    if(toupper(codon) %in% toupper(codons[[i]])) {
      return(names(codons)[i])
    }
  }
}

# got tired of writing this over and over again, so here we go
# returns vector of bits of original string where bits are a given size
break_up_string <- function(string_to_break, size_of_bits) {
  broken <- sapply(seq(from=1, to=nchar(string_to_break), by=size_of_bits), function(i) {
    return(substring(string_to_break, i, i+(size_of_bits-1)))
  })
  return(broken)
}

# helper function for get_aligned_cds
# converts DNA sequence into corresponding amino acid sequence
dna_to_aa <- function(dna_cds) {
  codonized <- break_up_string(dna_cds, 3)
  aas_vec <- character(0)
  for (i in 1:length(codonized)) {
    aa <- get_codon_name(codonized[i])
    aas_vec <- c(aas_vec, aa)
  }
  aas <- paste0(aas_vec, collapse = "") %>% 
    gsub(pattern = "Stop", replacement = "")
  return(aas)
}

# second helper function for get_aligned cds
# given an un-aligned dna sequence of a coding region and its aa alignment, 
# inserts the necessary gaps ("-") in the dna seq to make it aligned
insert_gaps <- function(gene, aa_aligned) {
  broken_aa <- break_up_string(aa_aligned, 1)
  gap_mask <-  broken_aa == "-"
  codonized <- break_up_string(gene, 3)
  gapped_dna <- broken_aa
  codon_idx <- 1
  for (i in 1:length(gap_mask)) {
    if (!gap_mask[i]) {
      gapped_dna[i] <- codonized[codon_idx]
      codon_idx <- codon_idx + 1
    }
    if (gap_mask[i]) {
      gapped_dna[i] <- "---"
    }
  }
  return(paste0(gapped_dna, collapse = ""))
}

library(msa)
# Yet a third helper function for calculate_ks
# This one returns aligned cds strings for pairs that are different lengths
# (which, as it turns out, are most of them)
# input: un-aligned cds of both genes (DNA)
# output: aligned cds of both genes (both are same length, includes "-" characters)
# as well as the number of nonsynonymous differences because it's really easy to calculate
# from the aa alignment
get_aligned_cds <- function(gene1, gene2) {
  aa1 <- dna_to_aa(gene1)
  aa2 <- dna_to_aa(gene2)
  aas_to_align <- AAStringSet(c(aa1, aa2))
  alignment <- msa(aas_to_align, cluster = 100, type = "protein", 
                   method = "ClustalOmega", maxiters = 0, 
                   substitutionMatrix = "BLOSUM40")
  alignment <- capture.output(print(alignment, show=c("alignment", "complete"), showNames=TRUE,
                                    showConsensus=TRUE, halfNrow=9, nameWidth=20))
  
  aligned1 <- alignment[grep("\\Q[1]\\E", alignment)]
  aligned1 <- paste0(gsub(pattern = "\\Q[1] \\E", replacement = "", 
                          x = aligned1), collapse = "")
  aligned1 <- gsub(pattern = " ", replacement = "", x = aligned1)
  aligned2 <- alignment[grep("\\Q[2]\\E", alignment)]
  aligned2 <- paste0(gsub(pattern = "\\Q[2] \\E", replacement = "", 
                          x = aligned2), collapse = "")
  aligned2 <- gsub(pattern = " ", replacement = "", x = aligned2)
  con <- alignment[grep("Con", alignment)]
  con <- paste0(gsub(pattern = "Con ", con, replacement = ""), collapse = "")
  nonsyn <- sum(break_up_string(con, 1) =="?")
  
  dna_aligned1 <- insert_gaps(gene1, aligned1)
  dna_aligned2 <- insert_gaps(gene2, aligned2)
  
  if(nchar(dna_aligned1) != nchar(dna_aligned2)) {
    cat("Aligned sequences are not the same length\n")
    return(NA)
  }
  return(list(dna_aligned1, dna_aligned2, nonsyn))
}

# Given SGD names of a pair of genes, returns K_s, 
# the number of synonomous substitutions per site
# input: CDS of both genes
# output: K_s, K_a, gene_length (for subs/site)
calculate_ks_ka <- function(sgd1, sgd2, species = c("Scer", "Spar")) {
  gene1 <- get_cds(sgd1, species)
  gene2 <- get_cds(sgd2, species)
  gene_length1 <- nchar(gene1)
  gene_length2 <- nchar(gene2)
  if(gene_length1 != gene_length2) {
    aligned_genes <- get_aligned_cds(gene1, gene2)
    aligned1 <- aligned_genes[[1]]
    aligned2 <- aligned_genes[[2]]
    nonsyn <- aligned_genes[[3]]
  }
  codonized1 <- break_up_string(aligned1, 3)
  codonized2 <- break_up_string(aligned2, 3)
  syn <- 0
  syn_sites <- 0
  for (i in 1:length(codonized1)) {
    cdn1 <- codonized1[i]
    cdn2 <- codonized2[i]
    # if there are any "-" gaps in the codon, we don't want it
    # (I could potentially modify this to check if the gap is only in the 3rd 
    # place and if so if it doesn't change the codon but that seems like a lot of work)
    if(grepl("-", cdn1) | grepl("-", cdn2)) {
      next
    }
    # omg there are spaces for some reason?!? I'll have to fix this later
    if(grepl(" ", cdn1) | grepl(" ", cdn2)) {
      next
    }
    #omfg also Ns sometimes fml
    if(grepl("N", cdn1) | grepl("N", cdn2)) {
      next
    }
    # also skip if we have a codon that's not length 3
    if((nchar(cdn1) != 3) | (nchar(cdn2) != 3)) {
      next
    }
    syn_sites <- syn_sites + 1 # only increment if they have codons that can be compared
    ndiffs <- check_codon(cdn1, cdn2)
    if(ndiffs != 0) {
      cdn1_name <- get_codon_name(cdn1)
      cdn2_name <- get_codon_name(cdn2)
      if(cdn1_name == cdn2_name) {
        syn <- syn + ndiffs
      }
    }
  }
  return(tibble("gene_name1" = sgd1,
                "gene_name2" = sgd2,
                "syn" = syn, 
                "nonsyn" = nonsyn, 
                "sites" = nchar(aligned1), 
                "gene1_length" = gene_length1,
                "gene2_length" = gene_length2,
                "syn_sites" = syn_sites))
}

# tests for calculate_ks_ka
tdh3_cer <- get_cds("YGR192C", "Scer")
tdh2_cer <- get_cds("YJR009C", "Scer")
tdh1_cer <- get_cds("YJL052W", "Scer")
calculate_ks_ka(tdh3_cer, tdh1_cer)

gene1 <- get_cds("YGR180C", "Spar") # pair that isn't same length
gene2 <- get_cds("YJL026W", "Spar")
calculate_ks_ka(gene1, gene2)
 
# @input: gene names and correlation method ("pearson", "spearman", "kendall")
# as character vectors
# @output: expression correlation between those two genes across the environments in count matrix
correlate_paralogs <- function(gene1, gene2, corrMethod, species = c("cer", "par")) {
  if (species == "cer") {
    countMatrix <- cpm_cer
  }
  if (species == "par") {
    countMatrix <- cpm_par
  }
  if (!is.na(gene2) & gene1 %in% rownames(countMatrix) & gene2 %in% rownames(countMatrix)) {
    exprVec1 <- countMatrix[rownames(countMatrix) == gene1,] %>% 
      as.numeric()
    exprVec2 <- countMatrix[rownames(countMatrix) == gene2,] %>% 
      as.numeric()
    return(cor(exprVec1, exprVec2, method=corrMethod))
  } 
  else {
    return(NA)
  }
}


#### IV. Identify diverging genes by expression ~ timepoint + species regression ####
library(MASS)

# function that determines if a gene of interest has significantly diverged based 
# on whether the species term coefficient is significant in the linear model:
# Expression level ~ timepoint + species
# input: geneName (sgd standard name), experimentGroup (1 of 4 options), and 
# whether or not to use parental counts (as opposed to hybrid allele counts)
# output: tibble of pvalue, beta coefficients, and goodness of fit for the regression (in the form of % deviance explained)
# (for now, this only includes WT samples, which is probably for the best when looking at divergence of the WT strains)
regress_by_timepoint <- function(geneName, experimentGroup=
                                   c("CellCycle", "HAP4andWTonYPD", 
                                     "SCtoLowPi", "YPD to Low N" )) {
  #format data for regressing
  cer_timepoint_expression_data <- tibble("expression_level" = dgelist_cer$counts[geneName,dgelist_cer$samples$experiment == experimentGroup & dgelist_cer$samples$genotype == "WT"],
                                          "species" = "cer")
  cer_timepoint_expression_data$timepoint <- dgelist_cer$samples$time_point[dgelist_cer$samples$experiment == experimentGroup & dgelist_cer$samples$genotype == "WT"] %>% 
    parse_number()
  
  par_timepoint_expression_data <- tibble("expression_level" = dgelist_par$counts[geneName, dgelist_par$samples$experiment == experimentGroup & dgelist_par$samples$genotype == "WT"],
                                          "species" = "par")
  par_timepoint_expression_data$timepoint <- dgelist_par$samples$time_point[dgelist_par$samples$experiment == experimentGroup & dgelist_par$samples$genotype == "WT"] %>% 
    parse_number()
  timepoint_expression_data <- rbind(cer_timepoint_expression_data, 
                                     par_timepoint_expression_data)
  
  result <- tryCatch({
    mod <- glm.nb(expression_level ~ timepoint + species, data = timepoint_expression_data, na.action = na.exclude)
    output <- tibble("species_beta" = coef(summary(mod))[3,1],
                     "species_pvalue" = coef(summary(mod))[3,4],
                     "percent_deviance_explained" = (1 - mod$deviance/mod$null.deviance),
                     "converged" = mod$converged)
    return(output)
  }, error = function(e) {
    output <- tibble("species_beta" = NA,
                     "species_pvalue" = NA,
                     "percent_deviance_explained" = NA,
                     "converged" = NA)
    return(output)
  })
  return(result)
}

# same thing as regress_by_timepoint, except it uses a poisson distribution instead
# of a negative binomial
regress_by_timepoint_poisson <- function(geneName, experimentGroup=
                                   c("CellCycle", "HAP4andWTonYPD", 
                                     "SCtoLowPi", "YPD to Low N" )) {
  #format data for regressing
  cer_timepoint_expression_data <- tibble("expression_level" = dgelist_cer$counts[geneName,dgelist_cer$samples$experiment == experimentGroup & dgelist_cer$samples$genotype == "WT"],
                                          "species" = "cer")
  cer_timepoint_expression_data$timepoint <- dgelist_cer$samples$time_point[dgelist_cer$samples$experiment == experimentGroup & dgelist_cer$samples$genotype == "WT"] %>% 
    parse_number()
  
  par_timepoint_expression_data <- tibble("expression_level" = dgelist_par$counts[geneName, dgelist_par$samples$experiment == experimentGroup & dgelist_par$samples$genotype == "WT"],
                                          "species" = "par")
  par_timepoint_expression_data$timepoint <- dgelist_par$samples$time_point[dgelist_par$samples$experiment == experimentGroup & dgelist_par$samples$genotype == "WT"] %>% 
    parse_number()
  timepoint_expression_data <- rbind(cer_timepoint_expression_data, 
                                     par_timepoint_expression_data)
  
  result <- tryCatch({
    mod <- glm(expression_level ~ timepoint + species, data = timepoint_expression_data, family = "poisson")
    # output <- tibble("species_beta" = coef(summary(mod))[3,1],
    #                  "species_pvalue" = coef(summary(mod))[3,4],
    #                  "percent_deviance_explained" = (1 - mod$deviance/mod$null.deviance),
    #                  "converged" = mod$converged)
    output <- mod
    return(output)
  }, error = function(e) {
    output <- tibble("species_beta" = NA,
                     "species_pvalue" = NA,
                     "percent_deviance_explained" = NA,
                     "converged" = NA)
    return(output)
  })
  return(result)
}

##################### 3. Objects #############

#### I. Regulatory Correlation ####

# wow these take a while! Like an hour basically. Definitely try to make more
# efficient if possible
# You're calculating the same correlation twice, so maybe fix that first?
corrs_cer <- sapply(c(1:nrow(dgelist_cer)), function(g_idx) {
  cat(g_idx, "/", nrow(dgelist_cer), "done\n")
  return(correlate_across_experiments(
    dgelist_cer[, dgelist_cer$samples$organism == "Saccharomyces cerevisiae"], 
    g_idx, groupNames))
})

corrs_par <- sapply(c(1:nrow(dgelist_par)), function(g_idx) {
  cat(g_idx, "/", nrow(dgelist_par), "done\n")
  return(correlate_across_experiments(
    dgelist_par[, dgelist_par$samples$organism == "Saccharomyces paradoxus"], 
    g_idx, groupNames))
})

rownames(corrs_cer) <- NULL # These are just the rownames for the first gene's partners, so they contain all but YAL001C and are out of order with what's actually in that row
rownames(corrs_par) <- NULL
colnames(corrs_cer) <- rownames(dgelist_cer)
colnames(corrs_par) <- rownames(dgelist_par)

# a metric of regulatory conservation - the between-species correlation of each 
# gene's correlation to all other genes within each species
reg_cons <- sapply(c(1:ncol(corrs_cer)), function(idx) {
  mask <- (is.na(corrs_cer[,idx]) | is.na(corrs_par[,idx]))
  cer <- corrs_cer[!mask, idx]
  par <- corrs_par[!mask, idx]
  corrs <- cor(cer, par, method = "pearson")
  return(corrs)
})

genes <- rownames(dgelist_cer)
names(reg_cons) <- genes
reg_cons <- reg_cons[!is.na(reg_cons)]
redgrave_diverged <- names(reg_cons[reg_cons < 0.2])

reg_div <- sapply(genes, function(gene) {
  reg_div_mask <- genes == gene
  reg_div <- 1 - abs(reg_cons[reg_div_mask])
  return(reg_div)
}) %>% unlist()

# one of two dataframes for downstream analysis. The second will be "pairs_df" - each row is list of paralog pair
genes_df <- tibble("gene" = genes,
                   "reg_div" = reg_div)

#### II. Paralog Expression Ratios ####
# goal: ID diverging pairs 
# (P1/P2 expr ratio is different btwn cer and par in some or all environments)

marchant_justparalogs <- marchant[marchant$Duplication != "S",]
groupNames <- c("CellCycle", "HAP4andWTonYPD", "SCtoLowPi", "YPD to Low N")

# First let's just look in the most normal environment---YPD at 0 h timepoint

WT_YPD_cer <- dgelist_cer[,dgelist_cer$samples$time_point == "0 h, YPD" & 
                            dgelist_cer$samples$genotype == "WT" &
                            dgelist_cer$samples$organism == "Saccharomyces cerevisiae"]
WT_YPD_par <- dgelist_par[,dgelist_par$samples$time_point == "0 h, YPD" & 
                            dgelist_par$samples$genotype == "WT" &
                            dgelist_par$samples$organism == "Saccharomyces paradoxus"]

# sample a random pair of WT samples---1 cer, 1 par
i <- sample(1:ncol(WT_YPD_cer), 1)
j <- sample(1:ncol(WT_YPD_par), 1)
wt_ypd_cer <- WT_YPD_cer$counts[,i]
wt_ypd_par <- WT_YPD_par$counts[,j]

# or, take average from all WT samples
wt_ypd_cer <- data.frame(counts = apply(WT_YPD_cer$counts, 1, mean)) # these have to be data frames instead of vectors so that get_expr_ratio doesn't return an incorrect # of dimensions error
wt_ypd_par <- data.frame(counts = apply(WT_YPD_par$counts, 1, mean))

# whichever count data you use, get tibble of expr ratios for each pair in Marchant list
marchant <- read.table("../marchant_paralogs.txt", header = TRUE)
expr_ratios <- tibble("P1" = marchant$P1,
                      "P2" = marchant$P2)
expr_ratios$cer <- map2_dbl(expr_ratios$P1, expr_ratios$P2, get_expr_ratio, counts = wt_ypd_cer)
expr_ratios$par <- map2_dbl(expr_ratios$P1, expr_ratios$P2, get_expr_ratio, counts = wt_ypd_par)

# get avg of both species' expression and plot
expr_ratios[expr_ratios == Inf] <- NA
expr_ratios$avg <- map2_dbl(expr_ratios$cer, expr_ratios$par, function(x, y) {
  return((x + y)/2)
})
ggplot(data = tibble(expr_ratios)) + geom_point(aes(x = log(cer), y = log(par), color = log(avg))) +
  geom_abline(aes(slope = 1, intercept = 0)) + scale_color_continuous(name = "average expression ratio") + xlab("expression ratio of pair in cerevisiae") +
  ylab("expression ratio of pair in paradoxus")

# orrrr take samples with same collection date
collectionDates <- unique(WT_YPD_cer$samples$collection_date)
collectionDates <- collectionDates[collectionDates %in% unique(WT_YPD_par$samples$collection_date)]
plots <- list()
expr_ratios <- tibble("P1" = marchant$P1,
                      "P2" = marchant$P2)
for (i in 1:length(collectionDates)) {
  cer_idxs <- which(WT_YPD_cer$samples$collection_date == collectionDates[i])
  par_idxs <- which(WT_YPD_par$samples$collection_date == collectionDates[i])
  wt_ypd_cer <- data.frame(counts = apply(data.frame(WT_YPD_cer$counts[,cer_idxs]), 1, mean)) # once again the data.frame() is so that apply doesn't freak out if given an n x 1 data frame
  wt_ypd_par <- data.frame(counts = apply(data.frame(WT_YPD_par$counts[,par_idxs]), 1, mean))
  expr_ratios$cer <- map2_dbl(expr_ratios$P1, expr_ratios$P2, get_expr_ratio, counts = wt_ypd_cer)
  expr_ratios$par <- map2_dbl(expr_ratios$P1, expr_ratios$P2, get_expr_ratio, counts = wt_ypd_par)
  expr_ratios[expr_ratios == Inf] <- NA
  expr_ratios$avg <- map2_dbl(expr_ratios$cer, expr_ratios$par, function(x, y) {
    return((x + y)/2)
  })
  plots[[i]] <- ggplot(data = tibble(expr_ratios)) + geom_point(aes(x = log(cer), y = log(par), color = log(avg))) + # ah this is average of expression RATIO. I think total expression LEVEL would be more informative
    geom_abline(aes(slope = 1, intercept = 0)) + scale_color_continuous()
}

library(gridExtra)
n <- length(plots)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(plots, ncol=nCol))

## significance testing for divergent pairs
library(parallel)
cpm_wt_cer <- cpm(WT_YPD_cer)
cpm_wt_par <- cpm(WT_YPD_par)
pair_pvalues_ypdtolown <- do.call(rbind, mclapply(c(1:nrow(marchant_justparalogs)), function(i) {
  get_expr_ratio_model(dgelist_cer$samples, 
                       cpm_cer, 
                       dgelist_par$samples, 
                       cpm_par, 
                       c(marchant_justparalogs$P1[i], marchant_justparalogs$P2[i]),
                       "YPD to Low N")
  }, mc.cores = 4))

pair_pvalues_cellcycle <- do.call(rbind, mclapply(c(1:nrow(marchant_justparalogs)), function(i) {
  get_expr_ratio_model(dgelist_cer$samples, 
                       cpm_cer, 
                       dgelist_par$samples, 
                       cpm_par, 
                       c(marchant_justparalogs$P1[i], marchant_justparalogs$P2[i]),
                       "CellCycle")
}, mc.cores = 4))

##Apply FDR correction
##First plot distribution of p-values to check for any weirdness
ggplot(pair_pvalues_ypdtolown, aes(x = p_value)) + geom_histogram(bins = 100) + ggtitle("p values")

##If all is well, run FDR correction
pair_pvalues_ypdtolown$q_value <- p.adjust(pair_pvalues_ypdtolown$p_value, method = "BH")
ggplot(pair_pvalues_ypdtolown, aes(x = q_value)) + geom_histogram(bins = 100) + ggtitle("p values")

expr_ratios <- expr_ratios[!is.na(expr_ratios$P2),]
expr_ratios$qvalue_ypdtolown <- pair_pvalues_ypdtolown$q_value
expr_ratios$sig_ypdtolown <- expr_ratios$qvalue_ypdtolown < 0.05/(nrow(expr_ratios))
ggplot(data = tibble(expr_ratios)) + geom_point(aes(x = log(cer), y = log(par), color = sig_ypdtolown)) +
  geom_abline(aes(slope = 1, intercept = 0)) + xlab("expression ratio of pair in cerevisiae") +
  ylab("expression ratio of pair in paradoxus")

# TODO: Henry says HP_ratio is supposed to be the same number for every gene,
# it is the null for two genes being diverged. This model should work for my data as written

# TODO: get true-replicate version of get_expr_ratio_model working for CellCycle
# and SC to Low Pi (b/c they have replicates explicitly stated) then move
# on to other two if needed

CC_dgelist_cer <- dgelist_cer[,dgelist_cer$samples$experiment == "CellCycle"]
CC_dgelist_par <- dgelist_par[,dgelist_par$samples$experiment == "CellCycle"]
CC_good_names_cer <- NULL
# collect all names in cer with exactly 2 replicates
for (i in 1:length(CC_dgelist_cer$samples$sample_name)) {
  query <- gsub("_rep[1-9]", CC_dgelist_cer$samples$sample_name[i], replacement="")
  names <- gsub("_rep[1-9]", CC_dgelist_cer$samples$sample_name, replacement="")
  if (sum(names == query) != 2) {
    cat(sum(names == query), "reps in", CC_dgelist_cer$samples$sample_name[i], "\n")
  }
  if (sum(names == query) == 2) {
    CC_good_names_cer <- c(CC_good_names_cer, query)
  }
}

CC_good_names_cer <- unique(CC_good_names_cer)

# make sure they all have exactly 2 reps in par too
to_remove_cer <- NULL
CC_good_names_par <- NULL
for (i in 1:length(CC_good_names_cer)) {
  query <-gsub("cer", CC_good_names_cer[i], replacement="par")
  names <- gsub("_rep[1-9]", CC_dgelist_par$samples$sample_name, replacement="")
  if (is.na(query)) {
    # TODO: there is a very annoying error where it hits this with "NA not in par" yet I have checked and none of the CC_good_names_cer values are NA!!
    # (might be index is one too high b/c I've removed one value)
    cat(CC_good_names_cer[i], "not in par\n")
  }
  if (sum(names == query) != 2) {
    cat(sum(names == query), "reps in", query, "\n")
    to_remove_cer <- c(to_remove_cer, i)
  }
  if (sum(names == query) == 2) {
    CC_good_names_par <- c(CC_good_names_par, query)
  }
}
CC_good_names_cer <- CC_good_names_cer[-to_remove_cer]

CC_pair_ratios <- sapply(c(1:length(CC_good_names_cer)), function(i) {
  rep1_cer_idx <- which(CC_dgelist_cer$samples$sample_name == paste0(CC_good_names_cer[i], "_rep1"))
  rep2_cer_idx <- which(CC_dgelist_cer$samples$sample_name == paste0(CC_good_names_cer[i], "_rep2"))

  rep1_par_idx <- which(CC_dgelist_par$samples$sample_name == paste0(CC_good_names_par[i], "_rep1"))
  rep2_par_idx <- which(CC_dgelist_par$samples$sample_name == paste0(CC_good_names_par[i], "_rep2"))
  
  CC_samples_cer <- CC_dgelist_cer[,c(rep1_cer_idx, rep2_cer_idx)]
  CC_samples_par <- CC_dgelist_par[,c(rep1_par_idx, rep2_par_idx)]
  
  CC_counts_cer <- cpm(CC_samples_cer)
  CC_counts_par <- cpm(CC_samples_par)
  
  return(get_expr_ratio_model_allGenes(CC_samples_cer, CC_counts_cer, CC_samples_par, CC_counts_par, experimentGroup = "CellCycle"))
  
})

#### III. dS and dN (K_s and K_a) (Anna I hate you for doing this) ####

# incorporating codeml-calculated dS values
marchant <- read.table("../../../marchant_paralogs.txt", header = TRUE)
setwd("~/Documents/Wittkopp_Lab/paralogs/awesome_power/coding/codeml/")
dS <- read.table("codon_aligned/dS.txt", header = FALSE, fill = TRUE)
dS <- dS[!is.na(dS$V3),]
colnames(dS) <- c("P1", "P2", "dS")
dS <- left_join(dS, marchant)
dS$R <- map2_dbl(dS$P1, dS$P2, correlate_paralogs, 
                 corrMethod = "pearson", species = "cer")
dS_unsat <- dS[dS$dS < 1.5,]
dS_kindasat <- dS[dS$dS < 8,]
ggplot(dS, aes(x = dS, y = R)) + geom_point(aes(color = Duplication)) + 
  xlab("Synonomous sequence divergence") + ylab("Regulatory Correlation")


# extract K_s value for each pair of paralogs in marchant
marchant <- read.table("../marchant_paralogs.txt", header = TRUE)
pairs <- marchant[marchant$Duplication == "WGD" | marchant$Duplication == "SSD",]

# these two take a solid 30 min I'd say
ka_ks_cer <- map2_dfr(pairs$P1, pairs$P2, function(x, y) {
  cat("begining on", x, "and", y, "\n")
  return(calculate_ks_ka(x, y, "Scer"))
})

ka_ks_par <- map2_dfr(pairs$P1, pairs$P2, function(x, y) {
  cat("begining on", x, "and", y, "\n")
  return(calculate_ks_ka(x, y, "Spar"))
})

# because calculate_ks_ka is computationally intensive, I initially output ks 
# and ka in a messy list that now needs unpacking
for (i in 1:length(ks_ka_cer)) {
  ks_cer <- ks_ka_cer[[i]][[1]]
  if (!is.na(ks_cer)) {
    ka_cer <- ks_ka_cer[[i]][[2]]
    sites_cer <- ks_ka_cer[[i]][[3]]
    pairs$ks_cer[i] <- ks_cer/sites_cer
    pairs$ka_cer[i] <- ka_cer/sites_cer
    pairs$nsites_cer[i] <- sites_cer
  }
  ks_par <- ks_ka_par[[i]][[1]]
  if (!is.na(ks_par)) {
    ka_par <- ks_ka_par[[i]][[2]]
    sites_par <- ks_ka_par[[i]][[3]]
    pairs$ks_par[i] <- ks_par/sites_par
    pairs$ka_par[i] <- ka_par/sites_par
    pairs$nsites_par[i] <- sites_par
  }
}

#### IV. Identify diverging genes by expression ~ timepoint + species regression ####

# cell cycle
regressions_by_timepoint_CellCycle <- map_dfr(rownames(dgelist_cer), function(x) {
  cat(which(rownames(dgelist_cer)==x), "/", nrow(dgelist_cer), "\n")
  return(regress_by_timepoint(x, "CellCycle"))
})
regressions_by_timepoint_CellCycle$experiment <- "CellCycle"

# testing model assumptions on a totally random gene (*ahem*)
# following the example in chapter 20.2 of elements of applied biostatistics: https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/generalized-linear-models-i-count-data.html
geneNames <- rownames(dgelist_cer)
mod_pois <- regress_by_timepoint_poisson(geneNames[sample(c(1:length(geneNames)), 1)], groupNames[sample(c(1:4),1)])
fake_data <- map_dfr(mod_pois$fitted.values, function(x) {
  return(rpois(1000, x))
}) # I still have NO IDEA why this doesn't make a tibble with 1000 COLUMNS
n <- length(mod_pois$y)
quant_resids <- map_dbl(c(1:n), function(i) {
  frac <- sum(fake_data[,i] < mod_pois$y[i])/1000
  return(frac)
})
sorted_quant_resids <- sort(quant_resids)
percentiles <- c(1:n)/n
theo_quant_resids <- qunif(percentiles)
plot(theo_quant_resids, sorted_quant_resids)

# HAP4andWTonYPD (but just WT, not the Hap4 deletion mutants)
regressions_by_timepoint_HAP4andWTonYPD <- map_dfr(rownames(dgelist_cer), function(x) {
  cat(which(rownames(dgelist_cer)==x), "/", nrow(dgelist_cer), "\n")
  return(regress_by_timepoint(x, "HAP4andWTonYPD"))
})
regressions_by_timepoint_HAP4andWTonYPD$experiment <- "HAP4andWTonYPD"

# SCtoLowPi
regressions_by_timepoint_SCtoLowPi <- map_dfr(rownames(dgelist_cer), function(x) {
  cat(which(rownames(dgelist_cer)==x), "/", nrow(dgelist_cer), "\n")
  return(regress_by_timepoint(x, "SCtoLowPi"))
})
regressions_by_timepoint_SCtoLowPi$experiment <- "SCtoLowPi"

# YPD to Low N
regressions_by_timepoint_YPDtoLowN <- map_dfr(rownames(dgelist_cer), function(x) {
  cat(which(rownames(dgelist_cer)==x), "/", nrow(dgelist_cer), "\n")
  return(regress_by_timepoint(x, "YPD to Low N"))
})
regressions_by_timepoint_YPDtoLowN$experiment <- "YPD to Low N"

# combine so that each gene still has only one row
names(regressions_by_timepoint_CellCycle) <- paste0("CellCycle_", names(regressions_by_timepoint_CellCycle))
names(regressions_by_timepoint_HAP4andWTonYPD) <- paste0("HAP4andWTonYPD_", names(regressions_by_timepoint_HAP4andWTonYPD))
names(regressions_by_timepoint_SCtoLowPi) <- paste0("SCtoLowPi_", names(regressions_by_timepoint_SCtoLowPi))
names(regressions_by_timepoint_YPDtoLowN) <- paste0("YPDtoLowN_", names(regressions_by_timepoint_YPDtoLowN))

regressions_by_timepoint <- cbind(regressions_by_timepoint_CellCycle,
                                  regressions_by_timepoint_HAP4andWTonYPD,
                                  regressions_by_timepoint_SCtoLowPi,
                                  regressions_by_timepoint_YPDtoLowN)

regressions_by_timepoint$gene <- rownames(dgelist_cer)

# now see which genes are significant in which samples (correcting for 5308 * 4 tests)
alpha <- 0.05/(nrow(dgelist_cer) * 4)
diverged_in_all_exps <- regressions_by_timepoint$gene[regressions_by_timepoint$CellCycle_species_pvalue < alpha &
                                                        regressions_by_timepoint$HAP4andWTonYPD_species_pvalue < alpha &
                                                        regressions_by_timepoint$SCtoLowPi_species_pvalue < alpha &
                                                        regressions_by_timepoint$YPDtoLowN_species_pvalue < alpha]

# fisher exact test for duplicate enrichment in genes that diverged in all 4 experiments
n_diverging_dupes <- sum(diverged_in_all_exps %in% duplicates)
diverging_dupes_enrichment_counts <- matrix(nrow = 2, c(n_diverging_dupes, length(diverged_in_all_exps)-n_diverging_dupes, length(duplicates), 6650-length(duplicates)))
fisher.test(diverging_dupes_enrichment_counts)

diverged_in_only_CellCycle <- regressions_by_timepoint$gene[regressions_by_timepoint$CellCycle_species_pvalue < alpha &
                                                              regressions_by_timepoint$HAP4andWTonYPD_species_pvalue >= alpha &
                                                              regressions_by_timepoint$SCtoLowPi_species_pvalue >= alpha &
                                                              regressions_by_timepoint$YPDtoLowN_species_pvalue >= alpha]

diverged_in_only_HAP4andWTonYPD <- regressions_by_timepoint$gene[regressions_by_timepoint$CellCycle_species_pvalue >= alpha &
                                                              regressions_by_timepoint$HAP4andWTonYPD_species_pvalue < alpha &
                                                              regressions_by_timepoint$SCtoLowPi_species_pvalue >= alpha &
                                                              regressions_by_timepoint$YPDtoLowN_species_pvalue >= alpha]

diverged_in_only_SCtoLowPi <- regressions_by_timepoint$gene[regressions_by_timepoint$CellCycle_species_pvalue >= alpha &
                                                                   regressions_by_timepoint$HAP4andWTonYPD_species_pvalue >= alpha &
                                                                   regressions_by_timepoint$SCtoLowPi_species_pvalue < alpha &
                                                                   regressions_by_timepoint$YPDtoLowN_species_pvalue >= alpha]

diverged_in_only_YPDtoLowN <- regressions_by_timepoint$gene[regressions_by_timepoint$CellCycle_species_pvalue >= alpha &
                                                              regressions_by_timepoint$HAP4andWTonYPD_species_pvalue >= alpha &
                                                              regressions_by_timepoint$SCtoLowPi_species_pvalue >= alpha &
                                                              regressions_by_timepoint$YPDtoLowN_species_pvalue < alpha]

diverged_in_only_1_exp <- c(diverged_in_only_CellCycle,
                            diverged_in_only_HAP4andWTonYPD,
                            diverged_in_only_SCtoLowPi,
                            diverged_in_only_YPDtoLowN)

# fisher's exact for just genes that diverged in one experiment
n_diverging_dupes <- sum(diverged_in_only_1_exp %in% duplicates)
diverging_dupes_enrichment_counts <- matrix(nrow = 2, c(n_diverging_dupes, length(diverged_in_only_1_exp)-n_diverging_dupes, length(duplicates), 6650-length(duplicates)))
fisher.test(diverging_dupes_enrichment_counts)
cell_11 <- diverging_dupes_enrichment_counts[1,1]
cell_21 <- diverging_dupes_enrichment_counts[2,1]
cell_22 <- diverging_dupes_enrichment_counts[2,2]
cell_12 <- diverging_dupes_enrichment_counts[1,2]

exact_test_plotting_data <- tibble("fraction_of_genes" = c(cell_11/(cell_11 + cell_21), cell_12/(cell_12 + cell_22)),
                                      "group" = c("...in diverging genes", "...in whole genome"))
ggplot(data = exact_test_plotting_data, aes(x = group, y = fraction_of_genes, fill = group)) + geom_bar(stat = "identity") +
  xlab("") + ylab("Percent of genes that are duplicates")

# TODO: recreate the Mattenberger plot in my SoT talk that compares set of duplicates
# to known singletons (instead of just all the rest of the genes). You can use the marchant list still

################# 4. Pretty Plots ###############
# Example of how reg divergence is calculated
sample_ij_pair <- tibble("i" = cpm_cer[45,],
                         "j" = cpm_cer[3002,])
ggplot(data = sample_ij_pair, aes(x = i, y = j)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0) + xlab("gene i") + ylab("gene j")

sample_cer_par_corr <- tibble("cer" = corrs_cer[,731],
                              "par" = corrs_par[,731])
ggplot(data = sample_cer_par_corr, aes(x = cer, y = par)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0) + xlab("cerevisiae") + ylab("paradoxus")

# Co-expression of candidate pairs with high and low reg divergence
# TODO: cannot currently plot cer and par co-expression profiles because they have different numbers
# of replicates of the same experiment and timepoint. I probably need to separate by experiment and convert
# timepoint into a quantitative variable so that it can be graphed
high_reg_div_pair <- tibble("cer" = cer_high_reg_div,
                            "par" = par_high_reg_div)
cer_high_reg_div <- cpm_cer["YKR091W", (dgelist_cer$samples$organism == "Saccharomyces cerevisiae") &
                              dgelist_cer$samples$experiment == "CellCycle"]
par_high_reg_div <- cpm_par["YKR091W", (dgelist_par$samples$organism == "Saccharomyces paradoxus") &
                              dgelist_par$samples$experiment == "CellCycle"]
par_high_reg_div <- par_high_reg_div[1:83]

ggplot(data=high_reg_div_pair, aes(x=cer, y=par)) + geom_line()
ggplot(data=sample_pair[sample_pair$condition == "SCtoLowPi_135 min, No Pi" |
                          sample_pair$condition == "CellCycle_230 min, YPD",], 
       aes(x=condition, y=value, group=gene_name)) + geom_point(aes(color=gene_name))

# pretty plots ks vs reg_div
m <- lm(ks_vs_reg_div$reg_div ~ ks_vs_reg_div$subs_cer)
summary(m)

# for cerevisiae
ggplot(data = ks_vs_reg_div, aes(x = ks_cer, y = reg_div)) + geom_point() +
  xlab("Synonomous substitutions per site") + 
  ylab("Regulatory divergence")

# for paradoxus
ggplot(data = ks_vs_reg_div, aes(x = ks_par, y = reg_div)) + geom_point() +
  xlab("Synonomous substitutions per site") + 
  ylab("Regulatory divergence")

pairs$expr_corr_cer <- map2_dbl(pairs$P1, pairs$P2, ~ 
                                 correlate_paralogs(.x, .y, "pearson", "cer"))
pairs$expr_corr_par <- map2_dbl(pairs$P1, pairs$P2, ~ 
                                  correlate_paralogs(.x, .y, "pearson", "par"))

# ks vs expr correlation of pair (correlation transformed as log((1+R)/(1-R)) to be compatible with Gu et al.)
ggplot(data = pairs, aes(x = ks_cer, y = log((1+expr_corr_cer)/(1-expr_corr_cer)))) + geom_point() +
  xlab("Synonomous substitutions per site") + 
  ylab("ln((1+R)/(1-R))")

# ks vs expr correlation - just WGDs
ggplot(data = pairs[pairs$Duplication == "WGD",], aes(x = ks_cer, y = log((1+expr_corr_cer)/(1-expr_corr_cer)))) + geom_point() +
  xlab("Synonomous substitutions per site") + 
  ylab("ln((1+R)/(1-R))") +
  ggtitle("WGDs only")

# ka vs expr correlation
ggplot(data = pairs, aes(x = ka_cer, y = log((1+expr_corr_cer)/(1-expr_corr_cer)))) + geom_point() +
  xlab("Non-synonomous substitutions per site") + 
  ylab("ln((1+R)/(1-R))")

m <- lm(log((1+expr_corr_cer)/(1-expr_corr_cer)) ~ ks_cer, data = pairs)
summary(m)

# ks for SSDs and WGDs histogram
ggplot(data = pairs) + 
  geom_histogram(data = pairs[pairs$Duplication == "SSD",], fill = "blue", alpha = 0.5, aes(x = ks_cer)) +
  geom_histogram(data = pairs[pairs$Duplication == "WGD",], fill = "red", alpha = 0.5, aes(x = ks_cer)) +
  xlab("Ks")

# ks vs max(reg div) of pair
ks_ka_cer_df$P1 <- paralog_pairs$P1
ks_ka_cer_df$P2 <- paralog_pairs$P2
ks_ka_cer_df$reg_div_P1 <- sapply(ks_ka_cer_df$P1, function(p) {
  rd <- as.numeric(reg_div_df[reg_div_df$gene == p, "reg_div"])
  if(length(rd) == 0) {
    return(NA)
  }
  else {
    return(rd)
  }
}) %>% unlist()

ks_ka_cer_df$reg_div_P2 <- sapply(ks_ka_cer_df$P2, function(p) {
  rd <- as.numeric(reg_div_df[reg_div_df$gene == p, "reg_div"])
  if(length(rd) == 0) {
    return(NA)
  }
  else {
    return(rd)
  }
}) %>% unlist()

ks_ka_cer_df$max_reg_div <- sapply(c(1:nrow(ks_ka_cer_df)), function(i) {
  return(max(ks_ka_cer_df$reg_div_P1[i], ks_ka_cer_df$reg_div_P2[i]))
})

ks_ka_cer_df$ks_per_site <- ks_ka_cer_df$ks/ks_ka_cer_df$sites
ks_ka_cer_df$ka_per_site <- ks_ka_cer_df$ka/ks_ka_cer_df$sites

ggplot(data = ks_ka_cer_df, aes(x = ks_per_site, y = max_reg_div, color = duplication_type)) + 
  geom_point() +
  xlab("Synonomous substitutions per site") + 
  ylab("Regulatory divergence (maximum of pair)")

# type of duplication vs reg divergence
m <- lm(ks_ka_cer_df$max_reg_div ~ ks_ka_cer_df$ks_per_site)
summary(m)

ks_ka_cer_df$duplication_type <- paralog_pairs$Duplication

p <- ggplot(data=pairs, aes(x=Duplication, y=ks_cer)) + 
  geom_boxplot() +
  geom_dotplot(binaxis = 'y',
               stackdir = 'center', 
               alpha = .75,
               binpositions = 'all',
               binwidth = 0.001,
               fill = "lightseagreen", 
               colour = "lightseagreen") +
  labs(x = "Type of duplication (WGD or SSD)",
       y = "Ks")
p

ggplot(data = ks_ka_cer_df, aes(x = duplication_type, y = max_reg_div)) + geom_boxplot()
