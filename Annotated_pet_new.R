source("https://bioconductor.org/biocLite.R")
BiocManager::install("MultiAssayExperiment")
source("https://bioconductor.org/biocLite.R")
BiocManager::install("SingleCellExperiment")
source("https://bioconductor.org/biocLite.R")
BiocManager::install("scater")
source("https://bioconductor.org/biocLite.R")
BiocManager::install("biomaRt")

suppressPackageStartupMessages({
library(MultiAssayExperiment)
library(SingleCellExperiment)
library(dplyr)
library(scater)
library(biomaRt)
library("AnnotationDbi")
library("org.Hs.eg.db")
})

#set working directory
setwd("/mnt/nfs/data/Vincent_Lab/Joel_Chappell/")

###Load MultiAssayExperiment object

#read data
(PetData <- readRDS("EMTAB3929.rds"))

## Extract the gene-level raw counts
Pet_counts <- assays(experiments(PetData)[["gene"]])[["count"]]


## Extract the phenotype data.
phn <- colData(PetData)

#Subset Genes and ERCC
Pet_ERCC<-Pet_counts[65127:65218,]
Pet_counts<-Pet_counts[1:65126,]

#Check for duplicated 
any(duplicated(rownames(Pet_counts)))
any(duplicated(rownames(Pet_ERCC)))

#Annotate gene names 

rownames(Pet_counts) <- gsub('\\..+$', '', rownames(Pet_counts))

rownames(Pet_counts) <- mapIds(org.Hs.eg.db,keys=rownames(Pet_counts),column="SYMBOL",keytype="ENSEMBL",multiVals="first")

#Remove NA's                          
Pet_counts <- subset(Pet_counts, rownames(Pet_counts) != "NA")

#Add ERCC's back in 
Pet_counts<-rbind(Pet_counts, Pet_ERCC)

#Create sce object
stopifnot(all(colnames(Pet_counts) == rownames(phn)))

Pet_sce <- SingleCellExperiment(
  assays = list(
    counts = Pet_counts,
    logcounts = log2((Pet_counts)+1)
  ),
  colData = phn
)

# define feature names in feature_symbol column
rowData(Pet_sce)$feature_symbol <- rownames(Pet_sce)

# remove features with duplicated names
Pet_sce <- Pet_sce[!duplicated(rowData(Pet_sce)$feature_symbol), ] #does remove both?

# define spike-ins
isSpike(Pet_sce, "ERCC") <- grepl("ERCC", rowData(Pet_sce)$feature_symbol)

#Exclude features that are not expressed
keep_features <- rowSums(counts(Pet_sce) > 0) > 0

#table(keep_features)
Pet_sce <- Pet_sce[keep_features, ]

#quick PCA on Embryonic day
plotPCA(Pet_sce, colour_by = "Characteristics.developmental.stage.")

#save object 
save(Pet_sce, file="Pet_sce.rds")







