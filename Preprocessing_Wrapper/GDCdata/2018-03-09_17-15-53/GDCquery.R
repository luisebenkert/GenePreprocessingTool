#! /usr/local/bin/Rscript
library(TCGAbiolinks)
library(dplyr)
library(DT)#package to create HTML pages from the data
library(SummarizedExperiment)#for colData function
library(DESeq2)
library(genefilter)

cat("----------TCGAbiolinks: START PROCESSING NEW _CLINICAL_ QUERY: TCGA-BRCA Clinical , barcode = c()----------
")

try({queryquery1 <- GDCquery(project= c("TCGA-BRCA"), data.category = "Clinical", barcode = c(""))

#download the actual data
GDCdownload(queryquery1)
clinical <- GDCprepare_clinic(query1, clinical.info = "patient")

for(i in c("admin","radiation","follow_up","drug","new_tumor_event")){
      message(i)
      aux <- GDCprepare_clinic(query1, clinical.info = i)
      if(is.null(aux)) next
      # add suffix manually if it already exists
      replicated <- which(grep("bcr_patient_barcode",colnames(aux), value = T,invert = T) %in% colnames(clinical))
      colnames(aux)[replicated] <- paste0(colnames(aux)[replicated],".",i)
      if(!is.null(aux)) clinical <- merge(clinical,aux,by = "bcr_patient_barcode", all = TRUE)
}

readr::write_csv(clinical,path = "./GDCdata/2018-03-09_17-15-53/TCGA-BRCA__clinical_from_XML.csv") # Save the clinical data into a csv file

})


cat("---------------------TCGAbiolinks: FINISHED QUERY PROCESSING---------------------
")
cat("----------TCGAbiolinks: START PROCESSING NEW _NORMAL_ QUERY: TCGA-BRCA Transcriptome Profiling , data.type = Gene Expression Quantification, workflow = HTSeq - Counts, barcode = c()----------
")

try({query3 <- GDCquery(project= c("TCGA-BRCA"), data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow = "HTSeq - Counts", barcode = c(""))

#download the actual data
GDCdownload(query3)
#prepares the data for analysis and putting it into an SummarizedExperiment object
data3 <- GDCprepare(query3)
#create DESeqDataSet for further processing
dds3 <- DESeqDataSet(data3, design = ~ 1)

#filter genes by count threshold (default: 0 counts)
#dds3 <- dds3[ rowSums(counts(dds3)) > 30, ]

#filter out genes
thres <- (ncol(dds3) * 30) / 100
dds3 <- dds3[ rowSums(counts(dds3) > 0) > thres, ]

#filter out samples
thres <- (nrow(dds3) * 30) / 100
dds3 <- dds3[ ,colSums(counts(dds3) > 0) > thres]

#logtransform counts and normalize by library size (=cross-sample comparison)
transdds3 <- vst(dds3)
selected_features3 <- transdds3
metadata3 <- colData(transdds3)
metadata3 <- metadata3[c("project_id", "shortLetterCode")]
transp_metadata3 <-t(as.matrix(metadata3))
write.csv(transp_metadata3, file="./GDCdata/2018-03-09_17-15-53/TCGA-BRCA__GeneExpressionQuantification_HTSeq-Counts_metadata.csv")
sel_meta3 <- colData(selected_features3)
sel_meta3 <- sel_meta3[c("project_id", "shortLetterCode")]
transp_sel_meta3 <- t(as.matrix(sel_meta3))
write.csv(transp_sel_meta3, file="./GDCdata/2018-03-09_17-15-53/TCGA-BRCA__GeneExpressionQuantification_HTSeq-Counts_metadata_topfeatures.csv")
matrix3 <- t(assay(transdds3)) #creates the gene expression matrix with genes as rows
write.csv(matrix3, file="./GDCdata/2018-03-09_17-15-53/TCGA-BRCA__GeneExpressionQuantification_HTSeq-Counts.csv")
final_matrix3 <- t(assay(selected_features3)) #creates the gene expression matrix with genes as rows
write.csv(final_matrix3, file="./GDCdata/2018-03-09_17-15-53/TCGA-BRCA__GeneExpressionQuantification_HTSeq-Counts_topfeatures.csv")
})


cat("---------------------TCGAbiolinks: FINISHED QUERY PROCESSING---------------------
")
cat("----------TCGAbiolinks: START PROCESSING NEW _NORMAL_ QUERY: TCGA-BRCA Transcriptome Profiling , data.type = Clinical, workflow = HTSeq - Counts, barcode = c()----------
")

try({query4 <- GDCquery(project= c("TCGA-BRCA"), data.category = "Transcriptome Profiling", data.type = "Clinical", workflow = "HTSeq - Counts", barcode = c(""))

#download the actual data
GDCdownload(query4)
#prepares the data for analysis and putting it into an SummarizedExperiment object
data4 <- GDCprepare(query4)
#create DESeqDataSet for further processing
dds4 <- DESeqDataSet(data4, design = ~ 1)

#filter genes by count threshold (default: 0 counts)
#dds4 <- dds4[ rowSums(counts(dds4)) > 30, ]

#filter out genes
thres <- (ncol(dds4) * 30) / 100
dds4 <- dds4[ rowSums(counts(dds4) > 0) > thres, ]

#filter out samples
thres <- (nrow(dds4) * 30) / 100
dds4 <- dds4[ ,colSums(counts(dds4) > 0) > thres]

#logtransform counts and normalize by library size (=cross-sample comparison)
transdds4 <- vst(dds4)
selected_features4 <- transdds4
metadata4 <- colData(transdds4)
metadata4 <- metadata4[c("project_id", "shortLetterCode")]
transp_metadata4 <-t(as.matrix(metadata4))
write.csv(transp_metadata4, file="./GDCdata/2018-03-09_17-15-53/TCGA-BRCA__Clinical_HTSeq-Counts_metadata.csv")
sel_meta4 <- colData(selected_features4)
sel_meta4 <- sel_meta4[c("project_id", "shortLetterCode")]
transp_sel_meta4 <- t(as.matrix(sel_meta4))
write.csv(transp_sel_meta4, file="./GDCdata/2018-03-09_17-15-53/TCGA-BRCA__Clinical_HTSeq-Counts_metadata_topfeatures.csv")
matrix4 <- t(assay(transdds4)) #creates the gene expression matrix with genes as rows
write.csv(matrix4, file="./GDCdata/2018-03-09_17-15-53/TCGA-BRCA__Clinical_HTSeq-Counts.csv")
final_matrix4 <- t(assay(selected_features4)) #creates the gene expression matrix with genes as rows
write.csv(final_matrix4, file="./GDCdata/2018-03-09_17-15-53/TCGA-BRCA__Clinical_HTSeq-Counts_topfeatures.csv")
})


cat("---------------------TCGAbiolinks: FINISHED QUERY PROCESSING---------------------
")
