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

readr::write_csv(clinical,path = "./GDCdata/2018-03-10_16-59-57/TCGA-BRCA__clinical_from_XML.csv") # Save the clinical data into a csv file

})


cat("---------------------TCGAbiolinks: FINISHED QUERY PROCESSING---------------------
")
