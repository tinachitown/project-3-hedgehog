setwd("/projectnb2/bf528/users/hedgehog/project_3/analyst")
library(dplyr)

# Reading in the data -----------------------------------------------------
## Example data
# DESeq_results <- read.csv("/project/bf528/project_3/results/example_deseq_results.csv")
# limma_results <- read.csv("/project/bf528/project_3/results/example_limma_results.csv")

## Real data
# Beta-naphthoflavone
beta_rna_results <- read.csv("/projectnb2/bf528/users/hedgehog/project_3/programmer/csv_files/resCMC_AHR_deseq_results.csv")
beta_mic_results <- read.csv("/projectnb2/bf528/users/hedgehog/project_3/analyst/limma_results_betaNap.csv")
# Econazole
econ_rna_results <- read.csv("/projectnb2/bf528/users/hedgehog/project_3/programmer/csv_files/resCORN_CARPXR_deseq_results.csv")
econ_mic_results <- read.csv("/projectnb2/bf528/users/hedgehog/project_3/analyst/limma_results_econ.csv")
# Thioacetamide
thio_rna_results <- read.csv("/projectnb2/bf528/users/hedgehog/project_3/programmer/csv_files/resSALINE_Cytotoxic_deseq_results.csv")
thio_mic_results <- read.csv("/projectnb2/bf528/users/hedgehog/project_3/analyst/limma_results_thio.csv")

# affy_map
refseq_mapping <- read.csv("/project/bf528/project_3/refseq_affy_map.csv")
refseq_mapping <- as.data.frame(refseq_mapping) # Mapping probesets / refseq IDs to genes


# Filtering for DEs --------------------------------------------------------
## Filtering criteria:
#   1. absolute FC > 1.5
#   2. p_uncorrected < 0.05
### Beta-naphthoflavone
beta_mic_de <- beta_mic_results %>% 
  dplyr::filter(exp(logFC)>1.5) %>% 
  dplyr::filter(P.Value < 0.05)
beta_rna_de <- beta_rna_results %>% 
  dplyr::filter(2^(log2FoldChange)>1.5) %>% 
  dplyr::filter(pvalue < 0.05)

### Econazole
econ_mic_de <- econ_mic_results %>% 
  dplyr::filter(exp(logFC)>1.5) %>% 
  dplyr::filter(P.Value < 0.05)
econ_rna_de <- econ_rna_results %>% 
  dplyr::filter(2^(log2FoldChange)>1.5) %>% 
  dplyr::filter(pvalue < 0.05)

### Thioacetamide
thio_mic_de <- thio_mic_results %>% 
  dplyr::filter(exp(logFC)>1.5) %>% 
  dplyr::filter(P.Value < 0.05)
thio_rna_de <- thio_rna_results %>% 
  dplyr::filter(2^(log2FoldChange)>1.5) %>% 
  dplyr::filter(pvalue < 0.05)


# Conversion functions ----------------------------------------------------
affy_to_refseq <- function(x){
  refseq <- refseq_mapping %>% 
    dplyr::filter(PROBEID == x) # uses probe ID
  refseq <- refseq[,1]
  if(length(refseq)!= 1){
    refseq <- NA
    return(refseq)
    } 
  else if(length(refseq) == 1){
    return(refseq)
  }
}

refseq_to_gene <- function(x){ 
  gene <- refseq_mapping %>% 
    dplyr::filter(REFSEQ == x) # uses refseq ID
  gene <- unique(gene[,3])
  return(gene)
}

# Getting RefSeqs ---------------------------------------------------------
## Example data
### Converting from affyID -> refseq
# limma_refs <- sapply(limma_results$X, affy_to_refseq, USE.NAMES = F)
# limma_results$REFSEQ <- limma_refs

## Real Data
### Beta-naphthoflavone
beta_mic_refs <- sapply(beta_mic_de$X, affy_to_refseq, USE.NAMES = F)
beta_mic_de$REFSEQ <- beta_mic_refs
### Econazole
econ_mic_refs <- sapply(econ_mic_de$X, affy_to_refseq, USE.NAMES = F)
econ_mic_de$REFSEQ <- econ_mic_refs
### Thioacetamide
thio_mic_refs <- sapply(thio_mic_de$X, affy_to_refseq, USE.NAMES = F)
thio_mic_de$REFSEQ <- thio_mic_refs

# Getting DEGs ------------------------------------------------------------
## Refseq -> Gene
### Beta-naphthoflavone
beta_mic_refs <- beta_mic_de$REFSEQ
beta_mic_refs <- na.omit(beta_mic_refs)
beta_mic_degs <- vapply(beta_mic_refs, refseq_to_gene, character(1), USE.NAMES = F)
beta_mic_degs <- na.omit(beta_mic_degs)
beta_mic_degs <- unique(beta_mic_degs)

beta_rna_degs <- vapply(beta_rna_de$X, refseq_to_gene, character(1), USE.NAMES = F)
beta_rna_degs <- na.omit(beta_rna_degs)
beta_rna_degs <- unique(beta_rna_degs)

### Econazole
#Microarray
econ_mic_refs <- econ_mic_de$REFSEQ
econ_mic_refs <- na.omit(econ_mic_refs)
econ_mic_degs <- vapply(econ_mic_refs, refseq_to_gene, character(1), USE.NAMES = F)
econ_mic_degs <- na.omit(econ_mic_degs)
econ_mic_degs <- unique(econ_mic_degs)

#RNA-Seq
econ_rna_degs <- vapply(econ_rna_de$X, refseq_to_gene, character(1), USE.NAMES = F)
econ_rna_degs <- na.omit(econ_rna_degs)
econ_rna_degs <- unique(econ_rna_degs)

### Thioacetamide
#Microarray
thio_mic_refs <- thio_mic_de$REFSEQ
thio_mic_refs <- na.omit(thio_mic_refs)
thio_mic_degs <- vapply(thio_mic_refs, refseq_to_gene, character(1), USE.NAMES = F)
thio_mic_degs <- na.omit(thio_mic_degs)
thio_mic_degs <- unique(thio_mic_degs)

#RNA-Seq
thio_rna_degs <- vapply(thio_rna_de$X, refseq_to_gene, character(1), USE.NAMES = F)
thio_rna_degs <- na.omit(thio_rna_degs)
thio_rna_degs <- unique(thio_rna_degs)


# Finding agreement in directionality --------------------------
agreement <- function(mic_data, rna_data){
  inter <- intersect(mic_data$REFSEQ, rna_data$X)
  agreed <- c()
  for(element in inter){
    sL <- sign(mic_data$logFC[which(mic_data$REFSEQ == element)])
    sL <- mean(sL)
    sD <- sign(rna_data$log2FoldChange[which(rna_data$X == element)])
    sD <- mean(sD)
    if(sL == sD){
      agreed <- c(element, agreed)
    }
  }
  return(agreed)
}

# Beta-naphthoflavone
beta_agreed <- agreement(beta_mic_de, beta_rna_de)
# Econazole
econ_agreed <- agreement(econ_mic_de, econ_rna_de)
# Thioacetamide
thio_agreed <- agreement(thio_mic_de, thio_rna_de)


# Calculating concordance -------------------------------------------------
## Example data
# # Checking if intersection values are the same
# geneL <- length(intersect(limma_degs, DESeq_degs))
# refseqL <- length(intersect(limma_results$REFSEQ, DESeq_results$X))
# geneL == refseqL
# # because this returns true, I can set n0 to 
# 
# n0 <- length(agreed)
# n1 <- length(DESeq_degs)
# n2 <- length(limma_degs)
# N <- 30000 # number of genes in the genome
# num <- (N*n0) - (n1*n2)
# denom <- N+n0-n1-n2
# corrected_int <- num/denom
# conc <- (2*(corrected_int))/(n1+n2)

## Real data
concordance <- function(agreed, rna_degs, mic_degs){
  n0 <- length(agreed)
  n1 <- length(rna_degs)
  n2 <- length(mic_degs)
  N <- 30000 # number of genes in the genome
  num <- (N*n0) - (n1*n2)
  denom <- N+n0-n1-n2
  corrected_int <- num/denom
  conc <- (2*(corrected_int))/(n1+n2)
  return(conc)
}

# beta
beta_conc <- concordance(beta_agreed, beta_rna_degs, beta_mic_degs)
# econ
econ_conc <- concordance(econ_agreed, econ_rna_degs, econ_mic_degs)
# thio
thio_conc <- concordance(thio_agreed, thio_rna_degs, thio_mic_degs)


# Plots -------------------------------------------------------------------
x <- c(length(beta_rna_degs), length(econ_rna_degs), length(thio_rna_degs))
y <- c(beta_conc, econ_conc, thio_conc)
plot(x,y)


