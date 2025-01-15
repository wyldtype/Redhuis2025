sapply(c("dplyr", "readr", "tidyr", "purrr", "ggplot2", "openxlsx", "matrixStats"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Redhuis2024/")

fay_cts <- read.csv("data_files/downloaded_from_Fay2023/htseq.ortho.csv")
sgd_lookup <- read_tsv("data_files/downloaded_from_SGD/SGD_features.tab", col_names = FALSE) |> 
  mutate(length = abs(X10-X11) + 1) |> 
  rename("systematic"="X4", "common"="X5") |> 
  select(systematic, common, length) |> 
  filter(!is.na(systematic) & grepl("^Y", systematic))
sgd_lookup$common <- sapply(c(1:nrow(sgd_lookup)), \(i) {
  if (is.na(sgd_lookup$common[i])) {
    return(sgd_lookup$systematic[i])
  } 
  else{
    return(sgd_lookup$common[i])
  }
})
sum(fay_cts$Name %in% sgd_lookup$common)
setdiff(fay_cts$Name, sgd_lookup$common)
# checking if any of the fay genes that aren't in the SGD lookup actually matter (i.e. are actually in the tagseq processed count matrix)
load("data_files/Cleaned_Barkai_Data.RData")
fay_cts <- left_join(fay_cts, select(sgd_lookup, common, systematic), by = c("Name"="common"))
FaynotinSGD <- fay_cts$systematic[which(!(fay_cts$systematic %in% sgd_lookup$systematic))]
setdiff(rownames(counts), sgd_lookup$systematic) # there is no Red not in SGD
FaynotinRed <- setdiff(fay_cts$systematic, rownames(counts))
sum(FaynotinSGD %in% FaynotinRed)/length(FaynotinSGD) # none of the fay genes with no annotation in SGD are in the other counts matrix, so they can safely be ignored
sum(is.na(fay_cts$systematic))
fay_cts <- fay_cts[!(fay_cts$systematic %in% FaynotinSGD),]
sum(is.na(fay_cts$systematic))
sum(duplicated(fay_cts$systematic)) # we'll remove these based on gene length below

fay <- select(fay_cts, c("systematic", "Len.Sc", "Len.Sp",
                         paste0("Sc", 25:32), paste0("Sc", 57:64), # cer rep1 and rep2
                         paste0("Sp", 25:32), paste0("Sp", 57:64), # par rep1 and rep2
                         paste0("Sc", 1:8), paste0("Sc", 33:40), # hyc rep1 and rep2
                         paste0("Sp", 1:8), paste0("Sp", 33:40))) # hyp rep1 and rep2

# create sample info to match format of other data
sample_info_fay <- tibble(sample_name = colnames(fay[,-c(1:3)]),
                          time_point_num = rep(c(0, 15, 30, 60), 16),
                          time_point_str = rep(c(paste("cold", c(0, 15, 30, 60), "min"), 
                                                 paste("heat", c(0, 15, 30, 60), "min")), 8),
                          organism = c(rep("cer", 16),
                                       rep("par", 16),
                                       rep("hyb", 16),
                                       rep("hyb", 16)),
                          allele= c(rep("cer", 16),
                                    rep("par", 16),
                                    rep("cer", 16),
                                    rep("par", 16)),
                          genotype = "WT",
                          well_flask_ID = rep(c(rep("rep1", 8), rep("rep2", 8)), 4))
sample_info_fay$experiment <- if_else(grepl("^heat", sample_info_fay$time_point_str), true = "Heat", false = "Cold")
sample_info_fay$condition <- paste(sample_info_fay$genotype, sample_info_fay$experiment, sample_info_fay$time_point_num, sep = "_")

#### removing small lib sizes ####
sort(colSums(fay[,-c(1:3)], na.rm = TRUE))
fay[colSums(fay[,-c(1:3)], na.rm = TRUE) > 100000)]
fay <- fay[, c(TRUE, TRUE, TRUE, colSums(fay[,-c(1:3)], na.rm = TRUE) > 100000)] # removes Sp26 and Sc64
sample_info_fay <- sample_info_fay |> filter(sample_name %in% colnames(fay))
nrow(sample_info_fay) # should now have 62 samples

# TODO: adapt tag-seq code to compare hybrid allelic libsizes
plotdf <- tibble(sample_name = names(libsizes[grepl("_hy[pc]_", names(libsizes))]),
                 libsize = libsizes[grepl("_hy[pc]_", names(libsizes))]) |> 
  left_join(y = select(filter(sample_info, organism == "hyb"), sample_name, experiment), by = "sample_name")
plotdf$allele <- if_else(grepl(pattern = "_hyc_", plotdf$sample_name),
                         true = "cer", false = "par")
plotdf$sample_name <- gsub("_hy[pc]_", "_hyb_", plotdf$sample_name)
plotdf <- plotdf |> pivot_wider(id_cols = c("sample_name", "experiment"), 
                                names_from = allele,
                                values_from = libsize)
ggplot(plotdf, aes(x = cer, y = par)) + 
  geom_point(aes(color = experiment), alpha = 0.5) +
  geom_text(data = filter(plotdf, cer > par*1.25 | par > cer*1.25),
            aes(label = sample_name), check_overlap = TRUE, color = "green") +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_rect(xmin = 0, xmax = 100000, ymin = 0, ymax = 100000, color = "blue", alpha = 0) +
  geom_vline(xintercept = 100000, color = "blue", alpha = 0.5) +
  geom_hline(yintercept = 100000, color = "blue", alpha = 0.5)

#### QC: comparing gene lengths in Fay annotations vs SGD ####
# some fay annotations for genes have slightly different (or dramatically different) lengths
# than the SGD annotations. Some of this may be differences in how annotations are called 
# (in which case I trust SGD more), but some may also be true genetic differences between
# S288C and T73 cerevisiae strains. Because I'm using S288C for all other datasets, I want
# to limit to the most comparable set of fay genes, the ones with lengths that match SGD

# checking for genes with a fay length outside of a 10% margin of SGD
fay <- left_join(fay, select(sgd_lookup, systematic, length), by = "systematic")
sum(fay$Len.Sc > fay$length*0.9 & fay$Len.Sc < fay$length*1.1, na.rm = TRUE)
sum(fay$Len.Sc < fay$length*0.9 | fay$Len.Sc > fay$length*1.1, na.rm = TRUE)
bad_length_idxs <- fay$Len.Sc < fay$length*0.9 | fay$Len.Sc > fay$length*1.1
fay[bad_length_idxs, c("systematic", "Len.Sc", "Len.Sp", "length")] |> View()

# removing duplicated gene rows---all 25 pairs of duplicate entries 
# had at most one copy that had the same SGD length
duplicate_genes <- fay$systematic[duplicated(fay$systematic)]
length(duplicate_genes)
fay[fay$systematic %in% duplicate_genes, 
    c("systematic", "Len.Sc", "Len.Sp", "length")] |> arrange(systematic)

# remove genes where fay length is outside of a 10% margin of SGD length
fay <- fay[(fay$Len.Sc > fay$length*0.9 & fay$Len.Sc < fay$length*1.1) & !(is.na(fay$Len.Sc)),]
sum(duplicated(fay$systematic)) # only one left
fay[fay$systematic == "YOL103W", 
    c("systematic", "Len.Sc", "Len.Sp", "length")]
fay <- fay[!(fay$systematic == "YOL103W" & fay$Len.Sc == 1701),]
sum(duplicated(fay$systematic))

# swap systematic column to be rownames
rownames(fay) <- fay$systematic
sum(colnames(fay[, setdiff(colnames(fay), c("systematic", "Len.Sc", "Len.Sp", "length"))]) == sample_info_fay$sample_name)
gene_lens <- fay[, c("systematic", "Len.Sc", "Len.Sp")] |> drop_na()
fay <- fay[gene_lens$systematic, setdiff(colnames(fay), c("systematic", "Len.Sc", "Len.Sp", "length"))] |>
  apply(1, as.numeric) |> t()
colnames(fay) <- sample_info_fay$sample_name

#### normalizing ####
# Fay et al. 2023 use DESeq2, so they're using un-normalized counts
# Normalizing to cpm taking account of gene length
# following this paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7373998/
# This is paried-end (one count = one pair), stranded, poly-A selected,
# so we should be safe to compare these counts to our tag-seq cpm counts
# (although if they weren't polyA selected or weren't stranded, direct comparison 
# isn't advisable b/c of biases in which subsets of the transcriptome are counted)

# normalizing with this strategy: 
# for each gene i, first obtain gene_i = 10^6 * (raw counts for gene_i/length(gene_i))
# then divide each gene by the new "library size" (the sum of these ^ across all genes in the sample)

# aka CPM = 10^6 * [(raw counts for gene_i/length(gene_i))/sum_i=1^i=nGenes(raw counts for gene_i/length(gene_i))]

countsPerMillionWithLength <- function(.cts, .lens) {
  rnames <- rownames(.cts)
  cnames <- colnames(.cts)
  genes_over_length <- map(rownames(.cts), \(g) {
    gene_vec <- .cts[g,] |> as.numeric()
    gene_len <- filter(.lens, systematic == g) |> select(length) |> pull()
    return(gene_vec/gene_len)
  }) |> purrr::reduce(.f = rbind)
  lib_sums <- colSums(genes_over_length)
  output <- apply(genes_over_length, 1, \(x) {return(x/lib_sums)}) |> t()
  rownames(output) <- rnames
  colnames(output) <- cnames
  return(round(output*10^6))
}
# tests for countsPerMillionWithLength
test_cts <- fay[,sample_info_fay$allele == "par"]
test_rowIdx <- sample(c(1:nrow(test_cts)), 1)
test_colIdx <- sample(c(1:ncol(test_cts)), 1)
test_count <- test_cts[test_rowIdx, test_colIdx]
test_lens <- gene_lens |> select(systematic, Len.Sp) |> rename("length"="Len.Sp")
test_cpm <- countsPerMillionWithLength(test_cts, test_lens)
# test_cpm_cpm <- countsPerMillionWithLength(test_cpm, test_lens)
round(((test_count/test_lens$length[test_rowIdx])/sum(test_cts[,test_colIdx]/test_lens$length))*10^6) # what it should be
test_cts[test_rowIdx, test_colIdx] # what it is before using our function
test_cpm[test_rowIdx, test_colIdx] # what it is using our function
# test_cpm_cpm[test_rowIdx, test_colIdx] # what it is if you run cpm too many times (unlike the TagSeq normalization, this cannot be run over and over again)

# normalizing (have to separate by species b/c gene lengths are different)
fay_cer <- countsPerMillionWithLength(.cts = fay[,sample_info_fay$allele == "cer"],
                                      .lens = rename(gene_lens[c("systematic", "Len.Sc")], "length"="Len.Sc"))
fay_par <- countsPerMillionWithLength(.cts = fay[,sample_info_fay$allele == "par"],
                                      .lens = rename(gene_lens[c("systematic", "Len.Sp")], "length"="Len.Sp"))

fay <- cbind(fay_cer, fay_par)
sum(colnames(fay) == sample_info_fay$sample_name)
sum(colnames(fay) %in% sample_info_fay$sample_name)/nrow(sample_info_fay)
fay <- fay[,sample_info_fay$sample_name]
sum(colnames(fay) == sample_info_fay$sample_name)

#### PCA ####
load("data_files/Cleaned_Redhuis_Data_AlleleSpecific.RData")
counts <- cbind(counts, counts_allele)
sample_info <- bind_rows(sample_info, sample_info_allele)
# the 0 min samples were grown at 25C in YPD, so they should be comparable to 
# the WT YPD samples from the main count data. So we want to see if the WT 0 min YPD
# samples cluster closer to their own species than they do to their own experiment (fay vs redhuis)
fay_0minYPD <- fay[,sample_info_fay$time_point_num == 0]
sample_info_LowN0minWT <- filter(sample_info, experiment == "LowN" & time_point_num == 0 & genotype == "WT")
LowN_0minYPD <- counts[,sample_info_LowN0minWT$sample_name]
sample_info_HAP40minWT <- filter(sample_info, experiment == "HAP4" & time_point_num == 0 & genotype == "WT")
HAP4_0minYPD <- counts[,sample_info_HAP40minWT$sample_name]
sample_info_LowPi0minWT <- filter(sample_info, experiment == "LowPi" & time_point_num == -5 & genotype == "WT")
LowPi_0minYPD <- counts[,sample_info_LowPi0minWT$sample_name]
sample_info_CC0minWT <- filter(sample_info, experiment == "CC" & time_point_num == 0 & genotype == "WT")
CC_0minYPD <- counts[,sample_info_CC0minWT$sample_name]

pca_info <- tibble(orgallele = c(rep("cercer", 4),
                                  rep("parpar", 4),
                                  rep("hybcer", 4),
                                  rep("hybpar", 4)),
                   experiment = "Temp")
pca_info <- bind_rows(pca_info, tibble(orgallele = paste0(sample_info_LowN0minWT$organism, 
                                                          sample_info_LowN0minWT$allele),
                                       experiment = "LowN"))
pca_info <- bind_rows(pca_info, tibble(orgallele = paste0(sample_info_HAP40minWT$organism, 
                                                          sample_info_HAP40minWT$allele),
                                       experiment = "HAP4"))
pca_info <- bind_rows(pca_info, tibble(orgallele = paste0(sample_info_LowPi0minWT$organism, 
                                                          sample_info_LowPi0minWT$allele),
                                       experiment = "LowPi"))
pca_info <- bind_rows(pca_info, tibble(orgallele = paste0(sample_info_CC0minWT$organism, 
                                                          sample_info_CC0minWT$allele),
                                       experiment = "CC"))

common_genes <- intersect(rownames(fay), rownames(counts))
sum(common_genes %in% rownames(fay)) 
sum(common_genes %in% rownames(counts)) 
# princomps
combcts <- cbind(fay_0minYPD[common_genes,], 
                 LowN_0minYPD[common_genes,],
                 HAP4_0minYPD[common_genes,],
                 LowPi_0minYPD[common_genes,],
                 CC_0minYPD[common_genes,])

covmat <- cov(combcts)
colnames(covmat) <- colnames(combcts)
pca_res <- princomp(covmat)
pcadf <- tibble(pc1 = pca_res$scores[,1], pc2 = pca_res$scores[,2])
pcadf <- bind_cols(pcadf, pca_info)
# summary(pca_res) # <- scroll through components to see that pc1 explains 80% of variance and pc2 explains 11%
ggplot(pcadf, aes(x = pc1, y = pc2)) + geom_point(aes(color = experiment)) +
  xlab("PC1, 80% of variance") + ylab("PC2, 11% of variance")
ggplot(pcadf, aes(x = pc1, y = pc2)) + geom_point(aes(color = orgallele)) +
  xlab("PC1, 80% of variance") + ylab("PC2, 11% of variance")
# pc2 solidly separates the Temp experiment from the others
# pc1 seems to separate the cerevisiae wine strain from s288C and paradoxus/hybrid
# but then again it also somewhat separates the paradoxus and hybrids,
# so it might just be picking up on a different aspect of the experiment
# in conclusion the experiment plays a much bigger role in overall transcriptome
# changes than species differences do

# repeating with genes as points to see if there are any outliers
# clustering genes by 4 tag-seq experiments versus all 6 
# experiments (in each species) and
# compare genes by which ones have the strongest difference in PC1
combcts <- cbind(counts[common_genes,], 
                 fay[common_genes,])
covmat4 <- cov(t(combcts[, c(1:ncol(counts))]))
colnames(covmat4) <- rownames(combcts)
pca_res4 <- princomp(covmat4)
covmat6 <- cov(t(combcts))
colnames(covmat6) <- rownames(combcts)
pca_res6 <- princomp(covmat6)
covmat2 <- cov(t(fay[common_genes,]))
colnames(covmat2) <- rownames(combcts)
pca_res2 <- princomp(covmat2)

pcadf <- tibble(pc1_4exp = pca_res4$scores[,1], pc1_6exp = pca_res6$scores[,1],
                pc2_4exp = pca_res4$scores[,2], pc2_6exp = pca_res6$scores[,2],
                pc1_2exp = pca_res2$scores[,1], pc2_2exp = pca_res2$scores[,2])
pcadf <- bind_cols(pcadf, tibble(gene_name = rownames(combcts)))
ggplot(pcadf, aes(x = pc1_4exp, y = pc1_2exp)) + geom_point()
ggplot(pcadf, aes(x = pc1_6exp, y = pc1_2exp)) + geom_point()
ggplot(pcadf, aes(x = pc1_4exp, y = pc1_6exp)) + geom_point() # Adding or removing the temp experiments barely changes PC1

# Conclusion: not enough variation in the Temp experiments to
# reveal any outlier genes, which also suggests there are no obvious outlier genes
# (main difference between strains are hgt genes in the Scer wine strain, 
# which aren't in our set of genes to compare)

# one gene seems to have higher 4exp PC1 than 6exp
sort((pcadf$pc1_4exp/pcadf$pc1_6exp)[pcadf$pc1_4exp > 1e6 & pcadf$pc1_6exp > 1e6], 
     decreasing = TRUE)[1:10] 
ggplot(pcadf, aes(x = pc1_4exp, y = pc1_6exp)) + 
  geom_point(alpha = 0.25) + 
  geom_point(data = filter(pcadf, gene_name == "YKL096W-A"),
             size = 2, color = "red", alpha = 0.5) # it's another cell wall mannoprotein, CWP2, (but not the same one that's misexpressed in the hybrid, that's YIL011W TIR3)

#### comparing gene mean expression and mean rank ####

# Do most genes still have the same mean expression for fay versus our main
# count data?
fay_cer <- fay[common_genes, sample_info_fay$organism == "cer"]
fay_par <- fay[common_genes, sample_info_fay$organism == "par"]
fay_hyc <- fay[common_genes, sample_info_fay$organism == "hyb" & sample_info_fay$allele == "cer"]
fay_hyp <- fay[common_genes, sample_info_fay$organism == "hyb" & sample_info_fay$allele == "par"]

red_cer <- counts[common_genes, sample_info$organism == "cer"]
red_par <- counts[common_genes, sample_info$organism == "par"]
red_hyc <- counts[common_genes, sample_info$organism == "hyb" & sample_info$allele == "cer"]
red_hyp <- counts[common_genes, sample_info$organism == "hyb" & sample_info$allele == "par"]

# raw means
# cer
plot(log2(rowMeans(fay_cer)), log2(rowMeans(red_cer)))
abline(a = 0, b = 1, col = "red")
# par
plot(log2(rowMeans(fay_par)), log2(rowMeans(red_par)))
abline(a = 0, b = 1, col = "red")
# hyc
plot(log2(rowMeans(fay_hyc)), log2(rowMeans(red_hyc)))
abline(a = 0, b = 1, col = "red")
# hyp
plot(log2(rowMeans(fay_hyp)), log2(rowMeans(red_hyp)))
abline(a = 0, b = 1, col = "red")

# conclusion: noisier than comparing alignments of barkai data,
# but the raw means have some relationship

#### splitting off hybrid and saving ####
fay_allele <- fay[,sample_info_fay$organism == "hyb"]
sample_info_fay_allele <- sample_info_fay |> filter(organism == "hyb")
fay <- fay[,sample_info_fay$organism != "hyb"]
sample_info_fay <- sample_info_fay |> filter(organism != "hyb")

save(fay, sample_info_fay, file = "data_files/Cleaned_Fay_Counts.RData")
save(fay_allele, sample_info_fay_allele, file = "data_files/Cleaned_Fay_Counts_Allele.RData")





