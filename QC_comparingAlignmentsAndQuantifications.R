sapply(c("dplyr", "purrr", "tidyr", "ggpubr", "readr", "data.table", "ggplot2", "DESeq2", "stringr", "openxlsx"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Redhuis2024//")

# goal: compare in-house versus published gene counts for 
# tagseq (Krieger et al. 2020 and Lupo et al. 2021) and rnaseq (Fay et al. 2023) data

#### Reading in Krieger et al. 2020 and Lupo et al. 2021 tagseq data ####
# read in count data (needs to be read in, organized, and normalized)
# (copied from QC_comparingAlignmentsAndQuantifications, which compared this
# alignment to the Krieger et al. 2020 alignment)
tagseq <- list.files("data_files/tagseq_counts/", full.names = TRUE) |> 
  map(read_table, col_names = FALSE, show_col_types = FALSE) |> 
  map(.f = select, X1, X3) |> # X1 are gene names, X3 is the sense strand read count
  purrr::reduce(.f = \(x, y) full_join(x = x, y = y, by = "X1")) |> 
  as.data.frame()

colnames(tagseq) <- c("gene", gsub("_ReadsPerGene.out.tab", "", list.files("data_files/tagseq_counts/", full.names = FALSE)))
QCdf <- tagseq[grepl("N_", tagseq$gene), ]
tagseq <- tagseq[!grepl("N_", tagseq$gene),]
tagseq <- tagseq[!tagseq$gene %in% c("cer_NA", "par_NA"),]
tagseq_cer <- tagseq[grepl("^cer_", tagseq$gene),]
tagseq_par <- tagseq[grepl("^par_", tagseq$gene),]
common_genes <- intersect(gsub("^cer_", "", tagseq_cer$gene),
                          gsub("^par_", "", tagseq_par$gene))
tagseq_cer$gene <- gsub("^cer_", "", tagseq_cer$gene)
tagseq_par$gene <- gsub("^par_", "", tagseq_par$gene)
tagseq_cer <- tagseq_cer[sapply(common_genes, \(x) which(x == tagseq_cer$gene)),]
tagseq_par <- tagseq_par[sapply(common_genes, \(x) which(x == tagseq_par$gene)),]
sum(tagseq_cer$gene == tagseq_par$gene)
length(common_genes)

#### Comparing tagseq counts to previously published count matrix ####
# read in data
krieger_cer <- read.xlsx("data_files/downloaded_from_Krieger2020/Supplemental_Table_S11.xlsx", sheet = 2, na.strings = "")
krieger_par <- read.xlsx("data_files/downloaded_from_Krieger2020/Supplemental_Table_S12.xlsx", sheet = 2, na.strings = "")

cat("percent count data columns with same names in both species: ",
    sum(colnames(krieger_cer)==colnames(krieger_par))/ncol(krieger_cer))

krieger_info <- read.xlsx("data_files/downloaded_from_Krieger2020/bioSample1to999.xlsx", na.strings="not applicable", cols=c(1,4,9,13,14,15,17)) %>%
  bind_rows(read.xlsx("data_files/downloaded_from_Krieger2020/bioSample1000toEnd.xlsx", na.strings="not applicable", cols=c(1,4,9,13,14,15,17)))
colnames(krieger_info) <- c("sample_name", "organism" , "collection_date", "genotype", "experiment","time_point", "well_flask_ID")
cat("percent sample info columns with same names in count data (before matching): ",
    sum(colnames(krieger_cer)[-1]==krieger_info$sample_name)/ncol(krieger_cer[,-1]))
desired_sample_order <- sapply(c(1:nrow(krieger_info)), function(i) {
  return(which(krieger_info$sample_name==colnames(krieger_cer[,-1])[i]))
})
krieger_info <- krieger_info[desired_sample_order,]
cat("percent sample info columns with same names in count data (after matching): ",
    sum(colnames(krieger_cer)[-1]==krieger_info$sample_name)/ncol(krieger_cer[,-1]))

# limiting to common_genes
common_genes_krieger <- intersect(krieger_cer[,1], common_genes)
# cer
sum(is.na(krieger_cer[,1]))
table(krieger_cer[,1] %in% common_genes_krieger, useNA = "always")
krieger_cer <- krieger_cer[krieger_cer[,1] %in% common_genes_krieger,]
sum(is.na(krieger_cer[,1]))
krieger_cer <- krieger_cer[match(krieger_cer[,1], common_genes_krieger),]
sum(krieger_cer[,1] == common_genes_krieger)
# par
krieger_par <- krieger_par[krieger_par[,1] %in% common_genes_krieger,]
krieger_par <- krieger_par[match(krieger_par[,1], common_genes_krieger),]
sum(krieger_par[,1] == common_genes_krieger)
# tagseq
tagseq_cer <- tagseq_cer[tagseq_cer[,1] %in% common_genes_krieger,]
tagseq_par <- tagseq_par[tagseq_par[,1] %in% common_genes_krieger,]
tagseq_cer <- tagseq_cer[match(common_genes_krieger, tagseq_cer[,1]),]
tagseq_par <- tagseq_par[match(common_genes_krieger, tagseq_par[,1]),]
sum(tagseq_cer[,1] == common_genes_krieger)
sum(tagseq_par[,1] == common_genes_krieger)

#### PCA plot ####

# checking that samples from the same experiment
# (CC, LowN, LowPi, HAP4) cluster with each other
# more than each other in all samples

# Scer allele counts
redhuis <- tagseq_cer[,-1]
barkai <- krieger_cer[,-1]

colnames(redhuis) <- paste0("red_", colnames(redhuis))
combcts <- cbind(barkai, redhuis)
covmat <- cov(combcts)
colnames(covmat) <- colnames(combcts)
pca_res <- princomp(covmat)
pcadf_cer <- tibble(pc1 = pca_res$scores[,1], pc2 = pca_res$scores[,2],
                    sample_name = gsub("^red_", "", colnames(combcts)))
pcadf_cer$lab <- if_else(grepl("^red_", colnames(combcts)), true = "redhuis", false = "barkai")
pcadf_cer$experiment <- if_else(grepl("TP[123]", colnames(combcts)), true = "LowN", false = 
                              if_else(grepl("CellCycle", colnames(combcts)), true = "CC", false =
                                        if_else(grepl("HAP4andWTonYPD", colnames(combcts)), true = "HAP4", false = "LowPi")))
# which samples have the farthest pc1/pc2 distance?
pcadf_cer |> select(sample_name, lab, pc1) |>  
  spread(lab, pc1) |> 
  group_by(sample_name) |> 
  summarise(pc1_diff = abs(barkai - redhuis)) |> 
  arrange(desc(pc1_diff)) # pc1
pcadf_cer |> select(sample_name, lab, pc2) |>  
  spread(lab, pc2) |> 
  group_by(sample_name) |> 
  summarise(pc2_diff = abs(barkai - redhuis)) |> 
  arrange(desc(pc2_diff)) # pc2
pcadf_cer$is_WT39rep2 <- grepl("WT_39", pcadf_cer$sample_name) & 
  grepl("rep2", pcadf_cer$sample_name)

ggplot(pcadf_cer, aes(x = pc1, y = pc2)) + geom_point(aes(color = lab))
ggplot(pcadf_cer, aes(x = pc1, y = pc2)) + geom_line(aes(group = sample_name, color = is_WT39rep2))
ggplot(pcadf_cer, aes(x = pc1, y = pc2)) + geom_point(aes(color = experiment))

# Spar allele counts
redhuis <- tagseq_par[,-1]
barkai <- krieger_par[,-1]

colnames(redhuis) <- paste0("red_", colnames(redhuis))
combcts <- cbind(barkai, redhuis)
covmat <- cov(combcts)
colnames(covmat) <- colnames(combcts)
pca_res <- princomp(covmat)
pcadf_par <- tibble(pc1 = pca_res$scores[,1], pc2 = pca_res$scores[,2],
                    sample_name = gsub("^red_", "", colnames(combcts)))
pcadf_par$lab <- if_else(grepl("^red_", colnames(combcts)), true = "redhuis", false = "barkai")
pcadf_par$experiment <- if_else(grepl("TP[123]", colnames(combcts)), true = "LowN", false = 
                              if_else(grepl("CellCycle", colnames(combcts)), true = "CC", false =
                                        if_else(grepl("HAP4andWTonYPD", colnames(combcts)), true = "HAP4", false = "LowPi")))
# which samples have the farthest pc1/pc2 distance?
pcadf_par |> select(sample_name, lab, pc1) |>  
  spread(lab, pc1) |> 
  group_by(sample_name) |> 
  summarise(pc1_diff = abs(barkai - redhuis)) |> 
  arrange(desc(pc1_diff)) # pc1
pcadf_par |> select(sample_name, lab, pc2) |>  
  spread(lab, pc2) |> 
  group_by(sample_name) |> 
  summarise(pc2_diff = abs(barkai - redhuis)) |> 
  arrange(desc(pc2_diff)) # pc2
pcadf_par$is_WT39rep2 <- grepl("WT_39", pcadf_par$sample_name) & 
  grepl("rep2", pcadf_par$sample_name)

ggplot(pcadf_par, aes(x = pc1, y = pc2)) + geom_point(aes(color = lab))
ggplot(pcadf_par, aes(x = pc1, y = pc2)) + geom_line(aes(group = sample_name, color = is_WT39rep2))
ggplot(pcadf_par, aes(x = pc1, y = pc2)) + geom_point(aes(color = experiment))

# Conclusion: all very similar except that WT_39_rep2 timepoint in all 3 organisms

# What's causing those CellCycle WT_39 rep2 samples to look so different?
# rep 1, not weird
plot(x = krieger_cer[,"WT_39_cer_CellCycle_rep1"],
     y = tagseq_cer[,"WT_39_cer_CellCycle_rep1"])
plot(x = krieger_par[,"WT_39_par_CellCycle_rep1"],
     y = tagseq_par[,"WT_39_par_CellCycle_rep1"])
plot(x = krieger_cer[,"WT_39_hyb_CellCycle_rep1"],
     y = tagseq_cer[,"WT_39_hyb_CellCycle_rep1"])
plot(x = krieger_par[,"WT_39_hyb_CellCycle_rep1"],
     y = tagseq_par[,"WT_39_hyb_CellCycle_rep1"])
# rep 2, weird
plot(x = krieger_cer[,"WT_39_cer_CellCycle_rep2"],
     y = tagseq_cer[,"WT_39_cer_CellCycle_rep2"])
plot(x = krieger_par[,"WT_39_par_CellCycle_rep2"],
     y = tagseq_par[,"WT_39_par_CellCycle_rep2"]) # Spar is fine actually
plot(x = krieger_cer[,"WT_39_hyb_CellCycle_rep2"],
     y = tagseq_cer[,"WT_39_hyb_CellCycle_rep2"])
plot(x = krieger_par[,"WT_39_hyb_CellCycle_rep2"],
     y = tagseq_par[,"WT_39_hyb_CellCycle_rep2"])

# Do WT_39_cer/hyb rep2 samples have smaller library sizes?
mask39_krieger <- grepl("WT_39", colnames(krieger_cer)) & grepl("rep2", colnames(krieger_cer))
mask39_tagseq <- grepl("WT_39", colnames(tagseq_cer)) & grepl("rep2", colnames(tagseq_cer))
par(mar = c(5,15,4,1))
barplot(colSums(cbind(krieger_cer[, "WT_39_cer_CellCycle_rep2"],
                      tagseq_cer[, "WT_39_cer_CellCycle_rep2"],
                      krieger_par[, "WT_39_cer_CellCycle_rep2"],
                      tagseq_par[, "WT_39_cer_CellCycle_rep2"],
                      krieger_cer[, "WT_39_hyb_CellCycle_rep2"],
                      tagseq_cer[, "WT_39_hyb_CellCycle_rep2"],
                      krieger_par[, "WT_39_hyb_CellCycle_rep2"],
                      tagseq_par[, "WT_39_hyb_CellCycle_rep2"],
                      krieger_par[, "WT_39_par_CellCycle_rep2"],
                      tagseq_par[, "WT_39_par_CellCycle_rep2"],
                      krieger_cer[, "WT_39_par_CellCycle_rep2"],
                      tagseq_cer[, "WT_39_par_CellCycle_rep2"])),
        names = c("cerSample, krieger, cer allele",
                  "cerSample, tagseq, cer allele",
                  "cerSample, krieger, par allele",
                  "cerSample, tagseq, par allele",
                  "hybSample, krieger, cer allele",
                  "hybSample, tagseq, cer allele",
                  "hybSample, krieger, par allele",
                  "hybSample, tagseq, par allele",
                  "parSample, krieger, par allele",
                  "parSample, tagseq, par allele",
                  "parSample, krieger, cer allele",
                  "parSample, tagseq, cer allele"),
        horiz = T, las = 1)

# It reaaally looks like the cer and hyb WT_39 rep2 samples got mislabled in tagseq
# Cer alleles, switching hyb/cer samples
plot(x = krieger_cer[,"WT_39_cer_CellCycle_rep2"],
     y = tagseq_cer[,"WT_39_hyb_CellCycle_rep2"])
plot(x = krieger_cer[,"WT_39_hyb_CellCycle_rep2"],
     y = tagseq_cer[,"WT_39_cer_CellCycle_rep2"])
# Par alleles, switching hyb/cer samples
plot(x = krieger_par[,"WT_39_cer_CellCycle_rep2"],
     y = tagseq_par[,"WT_39_hyb_CellCycle_rep2"]) # this one is the messiest b/c it's such small libsizes:
krieger_par[,"WT_39_cer_CellCycle_rep2"] |> sum()
tagseq_par[,"WT_39_hyb_CellCycle_rep2"] |> sum()
plot(x = krieger_par[,"WT_39_hyb_CellCycle_rep2"],
     y = tagseq_par[,"WT_39_cer_CellCycle_rep2"])

# Conclusion: WT_39_hyb_CellCycle_rep2 and WT_39_cer_CellCycle_rep2 samples are
# mislabled and should be switched (also re-downloaded from SRA to verify that it was
# also true for the publicly available raw files)

########################### Archive ########################### 
# #### comparing means and ranks of means of each gene's expression barkai vs redhuis ####
# load("data_files/Cleaned_Count_Data.RData")
# redhuis <- redhuis[intersect(rownames(redhuis), rownames(barkai)),
#                    intersect(colnames(redhuis), colnames(barkai))]
# barkai <- barkai[intersect(rownames(redhuis), rownames(barkai)),
#                  intersect(colnames(redhuis), colnames(barkai))]
# redhuis_avg_expr <- apply(redhuis, 1, mean)
# barkai_avg_expr <- apply(barkai, 1, mean)
# plot(rank(redhuis_avg_expr), rank(barkai_avg_expr)) # most genes have quite similar rank between species, but the ones that are different tend to be higher rank in redhuis!
# 
# # splitting into subsets of samples
# # mean
# plotMeansByGroup <- function(.colcat) {
#   red <- apply(redhuis[,grepl(.colcat, colnames(redhuis))], 1, mean)
#   bar <- apply(barkai[,grepl(.colcat, colnames(barkai))], 1, mean)
#   p <- ggplot(tibble(redhuis = log2(red), barkai = log2(bar)), aes(x = redhuis, y = barkai)) +
#     geom_point()
#   return(p)
# }
# # by species/allele
# plotMeansByGroup("_cer_")
# plotMeansByGroup("_par_")
# plotMeansByGroup("_hyc_")
# plotMeansByGroup("_hyp_") # not specific to any species/allele
# # by experiment
# plotMeansByGroup("TP")
# plotMeansByGroup("CellCycle")
# plotMeansByGroup("SCtoLowPi")
# plotMeansByGroup("HAP4andWTonYPD")
# # rank
# plotRanksByGroup <- function(.colcat) {
#   red <- apply(redhuis[,grepl(.colcat, colnames(redhuis))], 1, mean)
#   bar <- apply(barkai[,grepl(.colcat, colnames(barkai))], 1, mean)
#   p <- ggplot(tibble(redhuis = rank(red), barkai = rank(bar)), aes(x = redhuis, y = barkai)) +
#     geom_point()
#   return(p)
# }
# # by species/allele
# plotRanksByGroup("_cer_")
# plotRanksByGroup("_par_")
# plotRanksByGroup("_hyc_")
# plotRanksByGroup("_hyp_") # not specific to any species/allele
# # by experiment
# plotRanksByGroup("TP")
# plotRanksByGroup("CellCycle")
# plotRanksByGroup("SCtoLowPi")
# plotRanksByGroup("HAP4andWTonYPD") # also not specific to any experiment
# 
# # which genes are they? And what do their actual counts look like in individual samples?
# highredgenes <- rownames(redhuis)[which((rank(redhuis_avg_expr)/rank(barkai_avg_expr)) > 3)]
# redhuis[highredgenes,sample(c(1:ncol(redhuis)), 10)]
# barkai[highredgenes,sample(c(1:ncol(barkai)), 10)]
# # nothing crazy, just looks like there's a slight bias to have higher counts of low expressed genes in redhuis
# # making sure there aren't certain redhuis samples driving this. Unlikely b/c it's only seen in a few 100 genes and not all genes
# random_highred <- sample(highredgenes, 1)
# which.max(redhuis[random_highred,]) # rerun to check that all the highest counts aren't from the same sample
# highsampnames <- colnames(redhuis)[apply(redhuis[rownames(redhuis) %in% highredgenes,], 1, which.max)]
# highsampnames |> table() |> sort(decreasing = TRUE) # WT_39_hyp_CellCycle_rep2 is the highest for 
# libsizes <- colSums(counts)
# highsampsdf <- tibble("libsize" = libsizes[names(libsizes) %in% highsampnames], 
#                       "high" = TRUE) |>
#   bind_rows(tibble("libsize" = libsizes[!(names(libsizes) %in% highsampnames)],
#                    "high" = FALSE))
# ggplot(highsampsdf, aes(x = libsize)) + 
#   geom_density(aes(fill = high), alpha = 0.7) +
#   geom_vline(xintercept = 500000)
# # seems to be related to libsize, which makes sense. We haven't filtered for small libsize yet
# 
# #### visualizing cer/par whole genome alignment ####
# dots <- readr::read_tsv("../../test/dots.rdp")
# dots2 <- readr::read_tsv("../../test/dots2.rdp") # to check individual chromosome alignments
# plot(dots, type = "l") # note, the NAs are supposed to be there, they visually separate every start/end target-query pair
# # dataframe describing alignment blocks
# f1df <- readr::read_tsv("../../test/f1.txt")
# # no indels in alignment:
# sum(f1df$length1 == f1df$length2)/nrow(f1df)
# hist(f1df$length1, breaks = 50) # was hoping there would be an obvious small group of tiny sequences I could filter, but no
# sum(parse_number(f1df$name1)==parse_number(f1df$name2))/nrow(f1df) # but there are some alignments (~14%) from one chrom to another, which we're not expecting for these two species
# # do those non-canonical alignments tend to be shorter?
# f1df$same_chrom <- parse_number(f1df$name1)==parse_number(f1df$name2)
# ggplot(f1df, aes(x = length1)) + geom_histogram(aes(group = same_chrom, fill = same_chrom))
# # no! the distributions are identical
# # what's the longest different chromosome alignment?
# f1df |> filter(!same_chrom) |> arrange(desc(length1))
# # this is suspicious. why are they all the same exact length? Actually I can see that in the histogram, the bump above 5000
# f1df |> filter(length1 == 5471)
# f1df |> filter(length1 == 5471 & same_chrom) # none of them are same-chrom alignments
# # IDing r dot plot entries that are 5471 long
# dots5471 <- dots[0,]
# idxs5471 <- vector(mode = "integer", length = 0)
# for (i in which(is.na(dots$seq1))) {
#   start1 <- dots[i-2,"seq1"] |> as.numeric()
#   end1 <- dots[i-1,"seq1"] |> as.numeric()
#   seglength <- end1-start1+1
#   if (seglength == 5471) {
#     cat(i, seglength, "is 5471\n")
#     idxs5471 <- c(idxs5471, i)
#     dots5471 <- bind_rows(dots5471, dots[c(i-2, i-1, i),])
#   }
# }
# plot(dots, type="l")
# plot(dots5471, type="l") # the 5471 segments appear to be specific to 5 chromosomes, and on this scale they are indeed tiny
# # I can retrieve one (or a couple, to see how similar they are)
# # of the sequences and blast them against each yeast's genome
# f1df |> filter(length1 == 5471) # take some of these and awk the .maf file on the command line to print the real sequences
# 
# #### Inspecting annotation file ####
# 
# # TODO: there's now an R package that allows you
# # to do all this in a tidy-style: plyranges
# # that'll probably work better than my series of greps below
# 
# pargff <- read_table("../../test/CBS432.1000bpTES.gff", col_name = FALSE)
# cergff <- read_table("../../test/S288C.1000bpTES.gff", col_name = FALSE)
# 
# # extract unique ID from column 9
# pargff$ID <- map(pargff$X9, str_extract, pattern = "ID=CBS432_.{8}") |> 
#   map(gsub, pattern = "ID=CBS432_", replacement = "") |> 
#   map(gsub, pattern = "G", replacement = "T") |> unlist()
# cergff$ID <- map(cergff$X9, str_extract, pattern = "ID=S288c_.{8}") |> 
#   map(gsub, pattern = "ID=S288c_", replacement = "") |> 
#   map(gsub, pattern = "G", replacement = "T") |> unlist()
# 
# # you can check these manually, but they make the pivot not work below
# intronicIDs_par <- table(pargff$ID)[which(table(pargff$ID) != 4)] |> names()
# intronicIDs_cer <- table(cergff$ID)[which(table(cergff$ID) != 4)] |> names()
# 
# # do gene sequences always end 500bp after corresponding CDS?
# # par
# pardf <- pargff |> filter(X3 %in% c("gene", "CDS") & !(ID %in% intronicIDs_par)) |> 
#   mutate(TES = if_else(X7 == "-", true = X4, false = X5)) |> 
#   select(ID, X3, X7, TES) |> 
#   pivot_wider(id_cols = c("ID", "X7"),
#               names_from = "X3", values_from = c("TES"))
# pardf$TESext <- apply(pardf, 1, \(x) {
#   if (x["X7"] == "-") {
#     return(as.numeric(x["CDS"]) - as.numeric(x["gene"]))
#   }
#   if (x["X7"] == "+") {
#     return(as.numeric(x["gene"]) - as.numeric(x["CDS"]))
#   }
# })
# pardf$TESext |> table() # weee're gonna call that good
# 
# # cer
# cerdf <- cergff |> filter(X3 %in% c("gene", "CDS") & !(ID %in% intronicIDs_cer)) |> 
#   mutate(TES = if_else(X7 == "-", true = X4, false = X5)) |> 
#   select(ID, X3, X7, TES) |> 
#   pivot_wider(id_cols = c("ID", "X7"),
#               names_from = "X3", values_from = c("TES"))
# cerdf$TESext <- apply(cerdf, 1, \(x) {
#   if (x["X7"] == "-") {
#     return(as.numeric(x["CDS"]) - as.numeric(x["gene"]))
#   }
#   if (x["X7"] == "+") {
#     return(as.numeric(x["gene"]) - as.numeric(x["CDS"]))
#   }
# })
# cerdf$TESext |> table()
# 
# # manually checking intronic genes (mRNA rows X5 should be 500bp shorter than + genes and X4 should be 500bp longer than - genes)
# pargff |> filter(ID %in% intronicIDs_par) |> View()
# cergff |> filter(ID %in% intronicIDs_cer) |> View()