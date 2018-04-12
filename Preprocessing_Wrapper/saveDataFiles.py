#!/usr/bin/python

#Script that cares for all R script commands of storing any kind of data into a file.

import sys
import os, stat, subprocess

SCRIPT_WRITEOVERALLMETADATAFILE ="""metadata{0} <- colData(transdds{0})
metadata{0} <- metadata{0}[c("project_id", "shortLetterCode")]
transp_metadata{0} <-t(as.matrix(metadata{0}))
write.csv(transp_metadata{0}, file="{1}_metadata.csv")"""

SCRIPT_WRITEFINALMETADATAFILE ="""sel_meta{0} <- colData(selected_features{0})
sel_meta{0} <- sel_meta{0}[c("project_id", "shortLetterCode")]
transp_sel_meta{0} <- t(as.matrix(sel_meta{0}))
write.csv(transp_sel_meta{0}, file="{1}_metadata_topfeatures.csv")"""

SCRIPT_WRITEOVERALLFILE ="""matrix{0} <- t(assay(transdds{0})) #creates the gene expression matrix with genes as rows
write.csv(matrix{0}, file="{1}.csv")"""

SCRIPT_WRITEFINALFILE ="""final_matrix{0} <- t(assay(selected_features{0})) #creates the gene expression matrix with genes as rows
write.csv(final_matrix{0}, file="{1}_topfeatures.csv")"""

SCRIPT_PREP_LOOP = """for(i in c(\"admin\",\"radiation\",\"follow_up\",\"drug\",\"new_tumor_event\")){{
      message(i)
      aux <- GDCprepare_clinic({0}, clinical.info = i)
      if(is.null(aux)) next
      # add suffix manually if it already exists
      replicated <- which(grep(\"bcr_patient_barcode\",colnames(aux), value = T,invert = T) %in% colnames(clinical))
      colnames(aux)[replicated] <- paste0(colnames(aux)[replicated],\".\",i)
      if(!is.null(aux)) clinical <- merge(clinical,aux,by = \"bcr_patient_barcode\", all = TRUE)
}}"""
SCRIPT_WRITE_XML = "readr::write_csv(clinical,path = \"{0}_clinical_from_XML.csv\") # Save the clinical data into a csv file"


def write_metadatafiles(scriptfile, datasetid, outputfilepath):

    scriptfile.write(SCRIPT_WRITEOVERALLMETADATAFILE.format(datasetid, outputfilepath) + "\n")
    scriptfile.write(SCRIPT_WRITEFINALMETADATAFILE.format(datasetid, outputfilepath) + "\n")

def write_datafiles(scriptfile, datasetid, outputfilepath):

    scriptfile.write(SCRIPT_WRITEOVERALLFILE.format(datasetid, outputfilepath) + "\n")
    scriptfile.write(SCRIPT_WRITEFINALFILE.format(datasetid, outputfilepath) + "\n")

def write_clinicalfile(scriptfile, queryname, outputfilepath):
    scriptfile.write(SCRIPT_PREP_LOOP.format(queryname) + "\n")
    scriptfile.write("\n")
    scriptfile.write(SCRIPT_WRITE_XML.format(outputfilepath) + "\n")
    scriptfile.write("\n")
