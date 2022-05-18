library(tximport)
library(tximeta)
library(SummarizedExperiment)
###---> ThreeDRNAseq R package
library(ThreeDRNAseq)

###---> Denpendency R package
library(tximport)
library(edgeR)
library(limma)
library(RUVSeq)
library(eulerr)
library(gridExtra)
library(grid)
library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)
library(ensembldb)
options(stringsAsFactors=F)

#tutorial
files <- file.path("~/salmon_SRR1005811/quants/quant.sf")
file.exists(files)
tx2gene <- read.delim("~/lncrna_isoforms/tx2gene_grch38_ens94.txt")


dir <- system.file("~/salmon_tutorial/quants", package="tximportData")
file.exists(files)
coldata <- data.frame(files, names="SRR1005811", condition="A", stringsAsFactors=FALSE)
coldata
se <- tximeta(coldata,type = "salmon")
txi <- tximport(files, type="salmon",tx2gene = tx2gene[,c("tx_id", "ensgene")], countsFromAbundance="lengthScaledTPM",ignoreTxVersion = TRUE)

# Write the counts to an object
data <- txi$counts %>% 
  round() %>% 
  data.frame()

#dataset 

txi.seurat <- CreateSeuratObject(counts = txi$counts , min.cells = 3, min.features = 200, project = "SRR1005811")
