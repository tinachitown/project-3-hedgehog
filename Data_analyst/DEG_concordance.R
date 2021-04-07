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
  x <- na.omit(x)
  geneL <- c()
  for(refseq in x){
    gene <- refseq_mapping %>%
      dplyr::filter(REFSEQ == refseq) %>% # uses refseq ID
      dplyr::select(SYMBOL)
    gene <- unique(gene[,1])
    if(length(gene)==1){
      geneL <- c(gene, geneL)
    } else{
      gene <- NA
      geneL <- c(gene, geneL)
    }
  }
  geneL <- na.omit(geneL)
  geneL <- unique(geneL)
  return(geneL)
}



# Getting RefSeqs (from microarray data) ----------------------------------
beta_mic_de$REFSEQ <- sapply(beta_mic_de$X, affy_to_refseq, USE.NAMES = F)
econ_mic_de$REFSEQ <- sapply(econ_mic_de$X, affy_to_refseq, USE.NAMES = F)
thio_mic_de$REFSEQ <- sapply(thio_mic_de$X, affy_to_refseq, USE.NAMES = F)

# Getting DEGs ------------------------------------------------------------
# Beta
beta_mic_degs <- refseq_to_gene(beta_mic_de$REFSEQ)
beta_rna_degs <- refseq_to_gene(beta_rna_de$X)
# Econ
econ_mic_degs <- refseq_to_gene(econ_mic_de$REFSEQ)
econ_rna_degs <- refseq_to_gene(econ_rna_de$X)
# Thio
thio_mic_degs <- refseq_to_gene(thio_mic_de$REFSEQ)
thio_rna_degs <- refseq_to_gene(thio_rna_de$X)

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

# Concordance by above/below median ---------------------------------------
## Above Median
### Filtering
# beta-naphthoflavone
abmed_beta_mic <- beta_mic_de %>% 
  filter(AveExpr>median(AveExpr))
abmed_beta_rna <- beta_rna_de %>% 
  filter(baseMean>median(baseMean))
# econazole
abmed_econ_mic <- econ_mic_de %>% 
  filter(AveExpr>median(AveExpr))
abmed_econ_rna <- econ_rna_de %>% 
  filter(baseMean>median(baseMean))
# thioacetamide
abmed_thio_mic <- thio_mic_de %>% 
  filter(AveExpr>median(AveExpr))
abmed_thio_rna <- thio_rna_de %>% 
  filter(baseMean>median(baseMean))
### Finding DEGs
# Beta
abmed_beta_mic_degs <- refseq_to_gene(abmed_beta_mic$REFSEQ)
abmed_beta_rna_degs <- refseq_to_gene(abmed_beta_rna$X)
# Econ
abmed_econ_mic_degs <- refseq_to_gene(abmed_econ_mic$REFSEQ)
abmed_econ_rna_degs <- refseq_to_gene(abmed_econ_rna$X)
# Thio
abmed_thio_mic_degs <- refseq_to_gene(abmed_thio_mic$REFSEQ)
abmed_thio_rna_degs <- refseq_to_gene(abmed_thio_rna$X)
### Calculating concordance
# beta
abmed_beta_agreed <- agreement(abmed_beta_mic, abmed_beta_rna)
abmed_beta_conc <- concordance(abmed_beta_agreed, abmed_beta_rna_degs, abmed_beta_mic_degs)
# econ
abmed_econ_agreed <- agreement(abmed_econ_mic, abmed_econ_rna)
abmed_econ_conc <- concordance(abmed_econ_agreed, abmed_econ_rna_degs, abmed_econ_mic_degs)
# thio
abmed_thio_agreed <- agreement(abmed_thio_mic, abmed_thio_rna)
abmed_thio_conc <- concordance(abmed_thio_agreed, abmed_thio_rna_degs, abmed_thio_mic_degs)

## Below Median
### Filtering
# beta-naphthoflavone
bemed_beta_mic <- beta_mic_de %>% 
  filter(AveExpr<median(AveExpr))
bemed_beta_rna <- beta_rna_de %>% 
  filter(baseMean<median(baseMean))
# econazole
bemed_econ_mic <- econ_mic_de %>% 
  filter(AveExpr<median(AveExpr))
bemed_econ_rna <- econ_rna_de %>% 
  filter(baseMean<median(baseMean))
# thioacetamide
bemed_thio_mic <- thio_mic_de %>% 
  filter(AveExpr<median(AveExpr))
bemed_thio_rna <- thio_rna_de %>% 
  filter(baseMean<median(baseMean))
### Finding DEGs
# Beta
bemed_beta_mic_degs <- refseq_to_gene(bemed_beta_mic$REFSEQ)
bemed_beta_rna_degs <- refseq_to_gene(bemed_beta_rna$X)
# Econ
bemed_econ_mic_degs <- refseq_to_gene(bemed_econ_mic$REFSEQ)
bemed_econ_rna_degs <- refseq_to_gene(bemed_econ_rna$X)
# Thio
bemed_thio_mic_degs <- refseq_to_gene(bemed_thio_mic$REFSEQ)
bemed_thio_rna_degs <- refseq_to_gene(bemed_thio_rna$X)
### Calculating concordance
# beta
bemed_beta_agreed <- agreement(bemed_beta_mic, bemed_beta_rna)
bemed_beta_conc <- concordance(bemed_beta_agreed, bemed_beta_rna_degs, bemed_beta_mic_degs)
# econ
bemed_econ_agreed <- agreement(bemed_econ_mic, bemed_econ_rna)
bemed_econ_conc <- concordance(bemed_econ_agreed, bemed_econ_rna_degs, bemed_econ_mic_degs)
# thio
bemed_thio_agreed <- agreement(bemed_thio_mic, bemed_thio_rna)
bemed_thio_conc <- concordance(bemed_thio_agreed, bemed_thio_rna_degs, bemed_thio_mic_degs)


# Plots -------------------------------------------------------------------
jpeg("concordance_scatter.jpeg")
par(mfrow=c(1,2), oma=c(0,0,2,0))
## RNA-Seq
x1 <- c(length(beta_rna_degs), length(econ_rna_degs), length(thio_rna_degs))
y <- c(beta_conc, econ_conc, thio_conc)
bf_line <- stats::lm(y~x1)
names <- c("NAP", "ECO", "THI")
plot(x1,y,
     xlab = "Number of DEGs",
     ylab = "Concordance",
     main = "RNA-Seq")
lines(x1, fitted(bf_line),lty="dashed")
text(x1, c((y[1:2]+0.01), (y[3]-0.01)), labels = names)
## Microarray
x2 <- c(length(beta_mic_degs), length(econ_mic_degs), length(thio_mic_degs))
y <- c(beta_conc, econ_conc, thio_conc)
bf_line <- stats::lm(y~x2)
names <- c("NAP", "ECO", "THI")
plot(x2,y,
     xlab = "Number of DEGs",
     ylab = "Concordance",
     main = "Microarray")
lines(x2, fitted(bf_line),lty="dashed")
text(x2, c((y[1:2]+0.01), (y[3]-0.01)), labels = names)
mtext("Concordances by Platform", line = 0, outer=TRUE, cex=1.5)
dev.off()

# Barplot
ovr <- c(beta_conc, econ_conc, thio_conc)
abmed <- c(abmed_beta_conc, abmed_econ_conc, abmed_thio_conc)
bemed <- c(bemed_beta_conc, bemed_econ_conc, bemed_thio_conc)
c <- matrix(c(ovr,abmed,bemed), ncol=3, nrow = 3)
colnames(c) <- c("Overall", "Above Median", "Below Median")
jpeg("barplot.jpeg")
par(mfrow=c(1,1), oma=c(0,0,2,0))
barplot(c,beside=T, 
        col=c('blue', 'red', 'green'))
legend("topright", 
       c("NAP", "ECO", "THIO"),
       fill=c("blue", "red", 'green'))
mtext("Concordance by Average Expression Values", line = -1, outer=TRUE, cex=1.5)
dev.off()