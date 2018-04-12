#! /usr/local/bin/Rscript
library(TCGAbiolinks)
library(dplyr)
library(DT)#package to create HTML pages from the data
library(SummarizedExperiment)#for colData function
library(DESeq2)
library(genefilter)

cat("----------TCGAbiolinks: START PROCESSING NEW _NORMAL_ QUERY: TCGA-GBM Transcriptome Profiling , data.type = Gene Expression Quantification, workflow = HTSeq - Counts, barcode = c()----------
")

try({query1 <- GDCquery(project= c("TCGA-GBM"), data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow = "HTSeq - Counts", barcode = c(""))

#download the actual data
GDCdownload(query1)
#prepares the data for analysis and putting it into an SummarizedExperiment object
data1 <- GDCprepare(query1)
#create DESeqDataSet for further processing
dds1 <- DESeqDataSet(data1, design = ~ 1)

#filter genes by count threshold (default: 0 counts)
#dds1 <- dds1[ rowSums(counts(dds1)) > 30, ]

#filter out genes
thres <- (ncol(dds1) * 30) / 100
dds1 <- dds1[ rowSums(counts(dds1) > 0) > thres, ]

#filter out samples
thres <- (nrow(dds1) * 30) / 100
dds1 <- dds1[ ,colSums(counts(dds1) > 0) > thres]

#logtransform counts and normalize by library size (=cross-sample comparison)
transdds1 <- vst(dds1)
selected_features1 <- transdds1
metadata1 <- colData(transdds1)
metadata1 <- metadata1[c("project_id", "shortLetterCode")]
transp_metadata1 <-t(as.matrix(metadata1))
write.csv(transp_metadata1, file="./GDCdata/2018-01-22_15-36-26/TCGA-GBM__GeneExpressionQuantification_HTSeq-Counts_metadata.csv")
sel_meta1 <- colData(selected_features1)
sel_meta1 <- sel_meta1[c("project_id", "shortLetterCode")]
transp_sel_meta1 <- t(as.matrix(sel_meta1))
write.csv(transp_sel_meta1, file="./GDCdata/2018-01-22_15-36-26/TCGA-GBM__GeneExpressionQuantification_HTSeq-Counts_metadata_topfeatures.csv")
matrix1 <- t(assay(transdds1)) #creates the gene expression matrix with genes as rows
write.csv(matrix1, file="./GDCdata/2018-01-22_15-36-26/TCGA-GBM__GeneExpressionQuantification_HTSeq-Counts.csv")
final_matrix1 <- t(assay(selected_features1)) #creates the gene expression matrix with genes as rows
write.csv(final_matrix1, file="./GDCdata/2018-01-22_15-36-26/TCGA-GBM__GeneExpressionQuantification_HTSeq-Counts_topfeatures.csv")
})


cat("---------------------TCGAbiolinks: FINISHED QUERY PROCESSING---------------------
")
