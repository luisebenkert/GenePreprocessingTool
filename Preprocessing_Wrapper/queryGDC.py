#!/usr/bin/python

#Script for preparing a GDC query in an R script.

import sys
import os, stat, subprocess
import itertools

SAMPLE_TYPES = {"TP": "Primary solid Tumor", "TR": "Recurrent Solid Tumor",
"TB": "Primary Blood Derived Cancer - Peripheral Blood",
"TRBM": "Recurrent Blood Derived Cancer - Bone Marrow ",
"TAP": "Additional - New Primary", "TM": "Metastatic", "TAM": "Additional Metastatic",
"THOC": "Human Tumor Original Cells", "TBM": "Primary Blood Derived Cancer - Bone Marrow",
"NB": "Blood Derived Normal", "NT": "Solid Tissue Normal", "NBC": "Buccal Cell Normal",
"NEBV": "EBV Immortalized Normal", "NBM": "Bone Marrow Normal", "CELLC": "Control Analyte",
"TRB": "Recurrent Blood Derived Cancer - Peripheral Blood", "CELL": "Cell Lines",
"XP": "Primary Xenograft Tissue", "XCL": "Cell Line Derived Xenograft Tissue"}


SCRIPT_GDCQUERY ="""query{0} <- GDCquery(project= c(\"{1}\"), data.category = \"{2}\"{3})

#download the actual data
GDCdownload(query{0})"""

SCRIPT_PREPARE_NORMAL = """#prepares the data for analysis and putting it into an SummarizedExperiment object
data{0} <- GDCprepare(query{0})"""

SCRIPT_PREPARE_CLINICAL = "clinical <- GDCprepare_clinic({0}, clinical.info = \"patient\")"

def prepare_normalqueryparams(data_type, sample_types, workflow, bc_string):
    dt_string = ""
    wf_string = ""
    output_filename = ""
    st_string = ""

    #map sample type abbreviation into concrete definitions to use for the query
    if (sample_types != [""]):
        st_string = "sample.type = c("
        for s_type in sample_types:
            st_string += "\"" + SAMPLE_TYPES[s_type] + "\","
        st_string = st_string[:-1] + ")" #remove last comma
    else:
        st_string = ""

    if (data_type):
        dt_string = "data.type = \"" + data_type + "\""
        output_filename += "_" + data_type
    if (sample_types != [""]):
        output_filename += "_" + "_".join(sample_types)
    if (workflow):
        wf_string = "workflow = \"" + workflow + "\""
        output_filename += "_" + workflow

    other_params = ", ".join(filter(None, [dt_string, st_string, wf_string, bc_string]))
    if (other_params != ""):
        other_params = ", " + other_params

    return (other_params, output_filename)

def prepare_clinicalqueryparams(bc_string):

    other_params = ", ".join(filter(None, [bc_string]))
    if (other_params != ""):
        other_params = ", " + other_params

    return other_params

def write_GDCquery(scriptfile, datasetid, projects, data_category, other_params):

    scriptfile.write(SCRIPT_GDCQUERY.format(datasetid, projects, data_category, other_params) + "\n")
    scriptfile.write(SCRIPT_PREPARE_NORMAL.format(datasetid))
    scriptfile.write("\n")

def write_GDCclinicalquery(scriptfile, queryname, projects, data_category, other_params):
    scriptfile.write(SCRIPT_GDCQUERY.format(queryname, projects, data_category, other_params) + "\n")
    scriptfile.write(SCRIPT_PREPARE_CLINICAL.format(queryname) + "\n")
    scriptfile.write("\n")

def prepare_overallparameters(projects, barcodes, sample_types):
    #prepare parameter strings
    pj_string = "\", \"".join(projects)

    if (barcodes):
        bc_string = "barcode = c(\"" + "\", \"".join(barcodes) + "\")"
    else:
        bc_string = ""

    return (pj_string, bc_string)
