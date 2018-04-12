#! /usr/local/bin/Rscript
library(TCGAbiolinks)
library(dplyr)
library(DT)#package to create HTML pages from the data
library(SummarizedExperiment)#for colData function
library(DESeq2)
library(genefilter)

cat("----------TCGAbiolinks: START PROCESSING NEW _NORMAL_ QUERY: TCGA-UCEC", "TCGA-GBM Transcriptome Profiling , data.type = Gene Expression Quantification, workflow = HTSeq - Counts, barcode = c()----------
")

try({query1 <- GDCquery(project= c("TCGA-UCEC", "TCGA-GBM"), data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow = "HTSeq - Counts", barcode = c(""))

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
write.csv(transp_metadata1, file="./GDCdata/2018-03-10_20-43-23/TCGA-UCEC_TCGA-GBM__GeneExpressionQuantification_HTSeq-Counts_metadata.csv")
sel_meta1 <- colData(selected_features1)
sel_meta1 <- sel_meta1[c("project_id", "shortLetterCode")]
transp_sel_meta1 <- t(as.matrix(sel_meta1))
write.csv(transp_sel_meta1, file="./GDCdata/2018-03-10_20-43-23/TCGA-UCEC_TCGA-GBM__GeneExpressionQuantification_HTSeq-Counts_metadata_topfeatures.csv")
matrix1 <- t(assay(transdds1)) #creates the gene expression matrix with genes as rows
write.csv(matrix1, file="./GDCdata/2018-03-10_20-43-23/TCGA-UCEC_TCGA-GBM__GeneExpressionQuantification_HTSeq-Counts.csv")
final_matrix1 <- t(assay(selected_features1)) #creates the gene expression matrix with genes as rows
write.csv(final_matrix1, file="./GDCdata/2018-03-10_20-43-23/TCGA-UCEC_TCGA-GBM__GeneExpressionQuantification_HTSeq-Counts_topfeatures.csv")
})


cat("---------------------TCGAbiolinks: FINISHED QUERY PROCESSING---------------------
")
cat("----------TCGAbiolinks: START PROCESSING NEW _NORMAL_ QUERY: TCGA-UCEC", "TCGA-GBM Transcriptome Profiling , data.type = Clinical, workflow = HTSeq - Counts, barcode = c()----------
")

try({query2 <- GDCquery(project= c("TCGA-UCEC", "TCGA-GBM"), data.category = "Transcriptome Profiling", data.type = "Clinical", workflow = "HTSeq - Counts", barcode = c(""))

#download the actual data
GDCdownload(query2)
#prepares the data for analysis and putting it into an SummarizedExperiment object
data2 <- GDCprepare(query2)
#create DESeqDataSet for further processing
dds2 <- DESeqDataSet(data2, design = ~ 1)

#filter genes by count threshold (default: 0 counts)
#dds2 <- dds2[ rowSums(counts(dds2)) > 30, ]

#filter out genes
thres <- (ncol(dds2) * 30) / 100
dds2 <- dds2[ rowSums(counts(dds2) > 0) > thres, ]

#filter out samples
thres <- (nrow(dds2) * 30) / 100
dds2 <- dds2[ ,colSums(counts(dds2) > 0) > thres]

#logtransform counts and normalize by library size (=cross-sample comparison)
transdds2 <- vst(dds2)
selected_features2 <- transdds2
metadata2 <- colData(transdds2)
metadata2 <- metadata2[c("project_id", "shortLetterCode")]
transp_metadata2 <-t(as.matrix(metadata2))
write.csv(transp_metadata2, file="./GDCdata/2018-03-10_20-43-23/TCGA-UCEC_TCGA-GBM__Clinical_HTSeq-Counts_metadata.csv")
sel_meta2 <- colData(selected_features2)
sel_meta2 <- sel_meta2[c("project_id", "shortLetterCode")]
transp_sel_meta2 <- t(as.matrix(sel_meta2))
write.csv(transp_sel_meta2, file="./GDCdata/2018-03-10_20-43-23/TCGA-UCEC_TCGA-GBM__Clinical_HTSeq-Counts_metadata_topfeatures.csv")
matrix2 <- t(assay(transdds2)) #creates the gene expression matrix with genes as rows
write.csv(matrix2, file="./GDCdata/2018-03-10_20-43-23/TCGA-UCEC_TCGA-GBM__Clinical_HTSeq-Counts.csv")
final_matrix2 <- t(assay(selected_features2)) #creates the gene expression matrix with genes as rows
write.csv(final_matrix2, file="./GDCdata/2018-03-10_20-43-23/TCGA-UCEC_TCGA-GBM__Clinical_HTSeq-Counts_topfeatures.csv")
})


cat("---------------------TCGAbiolinks: FINISHED QUERY PROCESSING---------------------
")
cat("----------TCGAbiolinks: START PROCESSING NEW _CLINICAL_ QUERY: TCGA-UCEC", "TCGA-GBM Clinical , barcode = c()----------
")

try({queryquery3 <- GDCquery(project= c("TCGA-UCEC", "TCGA-GBM"), data.category = "Clinical", barcode = c(""))

#download the actual data
GDCdownload(queryquery3)
clinical <- GDCprepare_clinic(query3, clinical.info = "patient")

for(i in c("admin","radiation","follow_up","drug","new_tumor_event")){
      message(i)
      aux <- GDCprepare_clinic(query3, clinical.info = i)
      if(is.null(aux)) next
      # add suffix manually if it already exists
      replicated <- which(grep("bcr_patient_barcode",colnames(aux), value = T,invert = T) %in% colnames(clinical))
      colnames(aux)[replicated] <- paste0(colnames(aux)[replicated],".",i)
      if(!is.null(aux)) clinical <- merge(clinical,aux,by = "bcr_patient_barcode", all = TRUE)
}

readr::write_csv(clinical,path = "./GDCdata/2018-03-10_20-43-23/TCGA-UCEC_TCGA-GBM__clinical_from_XML.csv") # Save the clinical data into a csv file

})


cat("---------------------TCGAbiolinks: FINISHED QUERY PROCESSING---------------------
")
