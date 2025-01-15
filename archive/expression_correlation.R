rm(list=ls())
sapply(c("tidyverse", "edgeR", "openxlsx", "Glimma"), require, character.only=TRUE)
setwd("~/Documents/Wittkopp Lab/paralogs/test/output")

# following this Tag-seq sript: https://github.com/ben-laufer/Tag-seq/blob/master/tag-seq.R

####################### Barkai Count Data Import + Normalization
counts_strict_cer <- read.xlsx("../../barkailab_paper/Supplemental_Table_S11.xlsx", 
                               sheet = 2, na.strings = "")
geneNames <- counts_strict_cer[,1]
counts_strict_cer <- counts_strict_cer[,2:ncol(counts_strict_cer)] %>% 
  magrittr::set_rownames(geneNames)
counts_strict_cer[is.na(counts_strict_cer)] <- 0

counts_strict_par <- read.xlsx("../../barkailab_paper/Supplemental_Table_S12.xlsx", 
                               sheet = 2, na.strings = "")
geneNames_par <- counts_strict_par[,1]
sum(geneNames == geneNames_par) #just checking they're in the same order
counts_strict_par <- counts_strict_par[,2:ncol(counts_strict_par)] %>% 
  magrittr::set_rownames(geneNames_par)
counts_strict_par[is.na(counts_strict_par)] <- 0

sample_info <- read.xlsx("../../bioSample1to999.xlsx", na.strings="not applicable") %>% 
  bind_rows(read.xlsx("../../bioSample1000toEnd.xlsx", na.strings="not applicable"))
desired_sample_order <- sapply(c(1:nrow(sample_info)), function(i) {
  return(which(sample_info$sample_name==colnames(counts_strict_cer)[i]))
  })
sample_info <- sample_info[desired_sample_order,]
sum(sample_info$sample_name == colnames(counts_strict_cer))

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

##############################################
################## Barkai Count Data Analysis
##############################################

###################### Functions

# I. calculate regulatory divergence - 1 value per gene
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

# II. Calculate K_s and K_a between pairs of paralogs

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

##################### Objects

# I. Regulatory Divergence

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

# incorporating codeml-calculated dS values
setwd("Documents/Wittkopp_Lab/paralogs/awesome_power/coding/codeml/")
dS <- read.table("codon_aligned/dS.txt", header = FALSE, fill = TRUE)
dS <- dS[!is.na(dS$V3),]
colnames(dS) <- c("P1", "P2", "dS")
dS <- left_join(dS, marchant)
dS$R <- map2_dbl(dS$P1, dS$P2, correlate_paralogs, 
                       corrMethod = "pearson", species = "cer")
dS_unsat <- dS[dS$dS < 1.5,]
dS_kindasat <- dS[dS$dS < 8,]
ggplot(dS_unsat, aes(x = dS, y = R)) + geom_point(aes(color = Duplication)) 

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

reg_div <- sapply(genes, function(gene) {
  reg_div_mask <- genes == gene
  reg_div <- 1 - abs(reg_cons[reg_div_mask])
  return(reg_div)
}) %>% unlist()

# one of two dataframes for downstream analysis. The second will be "pairs_df" - each row is list of paralog pair
genes_df <- tibble("gene" = genes,
                   "reg_div" = reg_div)

# TODO: get Barkai reg divergence numbers from supplement and compare to yours

# II. K_s and K_a

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

# Graphing
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
