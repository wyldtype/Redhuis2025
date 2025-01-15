rm(list=ls())
sapply(c("tidyverse", "edgeR", "openxlsx", "Glimma"), require, character.only=TRUE)
setwd("~/Documents/Wittkopp Lab/paralogs/DDivergence/DataFromPapers/") #change based on your local computer

# read + organize Barkai datasets
cer_csv <- read_csv("counts_strict_cer.csv", na = "NaN")
par_csv <- read_csv("counts_strict_par.csv", na = "NaN")

cer_counts <- cer_csv[,-c(1,3)]
cer_counts[is.na(cer_counts)] <- 0
par_counts <- par_csv[,-c(1,3)]
par_counts[is.na(par_counts)] <- 0

colnames(cer_counts) <- c("gene_name", "P1_bar", "P2_bar", "P3_bar", "C1_bar", "CxP.C.1_bar", "CxP.C.2_bar", "CxP.C.3_bar",
                          "C2_bar", "C3_bar")
colnames(par_counts) <- c("gene_name", "P1_bar", "P2_bar", "P3_bar", "C1_bar", "CxP.P.1_bar", "CxP.P.2_bar", "CxP.P.3_bar",
                          "C2_bar", "C3_bar")
barkai_data <- bind_cols(cer_counts[,c("gene_name", "C1_bar", "C2_bar", "C3_bar", "CxP.C.1_bar", "CxP.C.2_bar", "CxP.C.3_bar")],
                 par_counts[,c("P1_bar", "P2_bar", "P3_bar", "CxP.P.1_bar", "CxP.P.2_bar", "CxP.P.3_bar")])

# read + organize Brem dataset
brem_data <- read_tsv("Metzger2017_supp2.txt", skip = 1, col_names = FALSE)
colnames(brem_data) <- c("gene_name", "C1_brem", "R_brem", "CxR.C_brem", "CxR.R_brem", "C2_brem", "P_brem", "CxP.C_brem",
                         "CxP.P_brem", "C3_brem", "M_brem", "CxM.C_brem", "CxM.M_brem", "C4_brem", "B_brem", "CxB.C_brem", "CxB.B_brem")
brem_data <- brem_data[,c(1,2,6:10,14)]
# join both datasets
gcm <- inner_join(brem_data, barkai_data, by = "gene_name")
gcm_untidy <- magrittr::set_rownames(gcm[,-1], gcm$gene_name) %>% as.matrix()

sample_info <- tibble(lab = c(rep("Brem", 7), rep("Barkai", 12)),
                      species = c("cer", "cer", "par", "hyb_cer", "hyb_par",
                                  "cer", "cer", "cer", "cer", "cer", "hyb_cer",
                                  "hyb_cer", "hyb_cer", "par", "par", "par", "hyb_par",
                                  "hyb_par", "hyb_par"))

# create DGElist object and normalize counts
dgelist0 <- DGEList(counts = gcm_untidy, samples = sample_info) %>% calcNormFactors()
dgelist0$samples
dgelist0$samples$group <- as.factor(dgelist0$samples$lab)
keep <- filterByExpr(dgelist0)
dgelist <- dgelist0[keep,,]
cpm <- cpm(dgelist)
lcpm <- cpm(dgelist, log = TRUE)
ggplot(data = dgelist$samples, aes(x = lib.size)) + geom_histogram(aes(fill = lab))

# create MDS plot (basically a PCA)
glMDSPlot(lcpm, labels=rownames(dgelist$samples), 
          groups=dgelist$samples[,c("lab", "species")], path = "../Barkai_data_analysis/", html = "MDS-downsampled")


