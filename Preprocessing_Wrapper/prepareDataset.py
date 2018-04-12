#!/usr/bin/python

#This python script creates an R script to query the GDC repository (NOT the legacy archive)
#via TCGAbiolinks for data sets according to the given specifications.
#It creates multiple output files per data category, sample, type etc.
#For example, if you want to download both gene expression and clinical data, you will receive
#two output files, one containing the gene expressions and the other containing
#the clinical data, respectively. The final result files and the generated R script
#will be stored in ./GDCdata/timestamp

#ATTENTION: This script assumes you have the following packages installed in R:
#--> TCGAbiolinks, dplyr, DT, SummarizedExperiment

#ATTENTION! If you specify a data type and have multiple data categories, make sure
#that the provided data types match the categories.
#What will not work: data.category = "Clinical" "Transcriptome Profiling", data.type = "Clinical"
#--> this will yield to queries with data.category = "Transcriptome Profiling", data.type = "Clinical"
#--> as there is no "Clinical" data type for transcriptome profiling, the query will result in an error

#ATTENTION! If you want to get a proper gene expression matrix at the end, make sure to specify a workflow type!
#-->If you do not specify it, data from all workflows will be downloaded
#--> a proper matrix cannot be created, you will only get the downloaded files.

#ATTENTION! Project and data category are mandatory!
#Options to specify:
#project - A list of valid project (it can be more than one) (see table below)
#data_category - A list of valid project (see list with getProjectSummary(project))
#data_type - A list of data types to filter the files to download
#sample_types - A list of abbreviated sample types, e.g. TP, TR (the concrete definitions will be mapped internally for the query). See list below for types.
#workflow - A list of GDC workflow type
#barcode - A list of barcodes to filter the files to download (can be partial barcodes)

#Additional parameters to specify:
#preprocessing - Specifies the preprocessing workflow to be used (currenty: DESeq2 and edgeR). Default is DESeq2.
#filtering_threshold - Specifies minimum percentage of samples of samples with non-zero values for a gene. Default is 30 (thus, genes with missing data in greater than 30% are filtered).
#feature_selection - Specifies the feature selection procedure to use. Default is rowVars from genefilter R package.
#fs_threshold - Specifies the number of genes to select as features. Default is 50.


#EXAMPLE: python prepareDataset.py --project "TCGA-BRCA" --data_category "Clinical" "Transcriptome Profiling" --data_type "Gene Expression Quantification" "Clinical" --workflow "HTSeq - Counts"

import sys
import os, stat, subprocess, datetime
import time
import argparse
import itertools

from saveDataFiles import *
from preprocessing import *
from queryGDC import *
from selectFeatures import *
#from discretization import *


WORKING_DIR = "./GDCdata/"#folder GDCdata is automatically created by TCGAbiolinks. we put all other data into it as well
DATASET_DIR = WORKING_DIR + datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S') +  "/"
SCRIPT_NAME = "GDCquery.R"

SCRIPT_PREAMBLE ="""#! /usr/local/bin/Rscript
library(TCGAbiolinks)
library(dplyr)
library(DT)#package to create HTML pages from the data
library(SummarizedExperiment)#for colData function
library(DESeq2)
library(genefilter)"""

SCRIPT_TRYSTART = "try({"
SCRIPT_TRYEND = "})\n"
SCRIPT_PRINT_START = "cat(\"----------TCGAbiolinks: START PROCESSING NEW {0} QUERY: {1}----------\n\")"
SCRIPT_PRINT_END = "cat(\"---------------------TCGAbiolinks: FINISHED QUERY PROCESSING---------------------\n\")"

def parse_arguments():
    try:
        print "Reading command line arguments..."
        parser = argparse.ArgumentParser(description="Preprocess gene expression data.")
        parser.add_argument("--project",  nargs="+", default=None, required = True,
            help="Specifies project names. At least one must be provided.")
        parser.add_argument("--data_category",  nargs="+", default=[""], required = True,
            help="Specifies data category. At least one must be provided.")
        parser.add_argument('--data_type', nargs="+", default=[""],
            help="Specifies data type. None specified results in multiple output files. Multiple specified results in multiple queries.")
        parser.add_argument('--sample_types', nargs="+", default=[""],
            help="Specifies sample types. Provide the abbreviations, e.g. TP for Tumor solid Primary. The abbreviations will be mapped to the concrete definitions internally.")
        parser.add_argument('--workflow', nargs="+", default=[""],
            help="Specifies workflow. None specified results in multiple output files. Multiple specified results in multiple queries.")
        parser.add_argument('--barcode', nargs="+", default=[""],
            help="Specifies barcodes. None specified results in download of all available barcode files.")
        parser.add_argument('--preprocessing',  nargs='?', default="DESeq2",
            help="Specifies preprocessing strategy. Default will be DESeq2.")
        parser.add_argument('--normalization',  nargs='?', default="vsd",
            help="Specifies normalization and log-transformation approach (vst or rlog). Default will be vsd. Use rlog for data sets with #samples < 100.")
        parser.add_argument('--feature_selection', nargs="?", default="rowVars",
            help="Specifies strategy for feature selection. Default will be selection by rowVars.")
        parser.add_argument('--fs_threshold', nargs="?", default="0",
            help="Specifies number of genes to select. Default will be all genes.")
        parser.add_argument('--filtering_threshold', nargs="?", default="30",
            help="Specifies minimum percentage of samples of samples with non-zero values for a gene. Default is 30 (thus, genes with missing data in greater than 30% are filtered).")
        parser.add_argument('--enable_discretization', nargs="?", default="no",
            help="Specifies whether to use discretization. Default will be no.")
        parser.add_argument('--discretization', nargs="?", default="mean",
            help="Specifies discretization strategy. Default will be to discretize according to the mean value.")

        global projects
        global datacats
        global datatypes
        global sample_types
        global workflows
        global barcodes
        global filtering_threshold
        global normalization
        global fs_threshold
        global preprocessing_pipeline
        global feature_selection
        global discretization
        global enable_discretization

        args = parser.parse_args()
        projects = args.project
        datacats = args.data_category
        datatypes = args.data_type
        sample_types = args.sample_types
        workflows = args.workflow
        barcodes = args.barcode
        preprocessing_pipeline = args.preprocessing
        normalization = args.normalization
        filtering_threshold = args.filtering_threshold
        feature_selection = args.feature_selection
        fs_threshold = args.fs_threshold
        enable_discretization = (args.enable_discretization == "yes")
        discretization = args.discretization

        print "... done."
    except:
        print "Stopping execution. An error occurred: ", sys.exc_info()
        exit()

def check_if_directories_exist():
    try:
        print "Checking if directories " + WORKING_DIR + " and " + DATASET_DIR + " exist..."
        #create working dir - all files produced by this script go here
        if not (os.path.isdir(DATASET_DIR)):
            os.makedirs(DATASET_DIR)
            os.chmod(DATASET_DIR, stat.S_IRWXU)
        print "...done."
    except:
        print "Stopping execution. An error occurred: ", sys.exc_info()
        exit()

def write_processingstart(scriptfile, querytype, projects, data_category, other_params):
    scriptfile.write(SCRIPT_PRINT_START.format(querytype, projects + " " + data_category + " " + other_params.replace("\"", ""))+ "\n")
    scriptfile.write("\n")
    scriptfile.write(SCRIPT_TRYSTART)

def write_preamble(scriptfile):
    scriptfile.write(SCRIPT_PREAMBLE + "\n")
    scriptfile.write("\n")

def write_processingending(scriptfile):
    scriptfile.write(SCRIPT_TRYEND + "\n")
    scriptfile.write("\n")
    scriptfile.write(SCRIPT_PRINT_END + "\n")

def write_normal_processing(scriptfile, index, pj_string, data_category, other_params, outputfile):

    datasetid =  str(index)
    outputfilepath = DATASET_DIR + outputfile

    write_processingstart(scriptfile, "_NORMAL_", pj_string, data_category, other_params)

    write_GDCquery(scriptfile, datasetid, pj_string, data_category, other_params)

    preprocess(scriptfile, datasetid, preprocessing_pipeline, filtering_threshold, normalization)

    apply_featureselection(scriptfile, datasetid, feature_selection, fs_threshold)

    if enable_discretization:
        discretize(scriptfile, datasetid, discretization)
    write_metadatafiles(scriptfile, datasetid, outputfilepath)
    write_datafiles(scriptfile, datasetid, outputfilepath)

    write_processingending(scriptfile)

def write_clinical_processing(scriptfile, index, pj_string, data_category, other_params, outputfile):
    queryname = "query" + str(index)
    outputfilepath = DATASET_DIR + outputfile

    write_processingstart(scriptfile, "_CLINICAL_", pj_string, data_category, other_params)

    write_GDCclinicalquery(scriptfile, queryname, pj_string, data_category, other_params)

    write_clinicalfile(scriptfile, queryname, outputfilepath)

    write_processingending(scriptfile)

def iterate_datasets(scriptfile, pj_string, bc_string, outputfile):

    #create a list with all possible combinations of parameters
    #e.g. data.category = "Clinical" "Transcriptome Profiling", data.type = "Clinical" "Gene Expression Quantification"
    #results in [("Clinical", "Clinical"),
    #   ("Clinical", Gene Expression Quantification"),
    #   ("Transcriptome Profiling", "Clinical"),
    #   ("Transcriptome Profiling", "Gene Expression Quantification")]
    parametercombinations = itertools.product(datacats, datatypes, workflows)

    biospecimendone = False
    clinicaldone = False
    #create a query for each combination
    index = 0
    for combination in parametercombinations:
        index += 1
        data_category = combination[0]
        data_type = combination[1]
        workflow = combination[2]

        clinical_data = (data_category == "Clinical")
        biospecimen_data = (data_category == "Biospecimen")
        #if we have clinical data, we need to apply a different script
        if not (clinical_data or biospecimen_data):
            params_filename = prepare_normalqueryparams(data_type, sample_types, workflow, bc_string)
            other_params = params_filename[0]
            write_normal_processing(scriptfile, index, pj_string, data_category, other_params, (outputfile + params_filename[1]).replace(" ", ""))
        else:
            #we only need to create one query for clinical or biospecimen data, otherwise we
            #get duplicate result files because we have multiple combinations with clinical/biospecimen as data_category
            if (clinical_data and not (clinicaldone)) or (biospecimen_data and not (biospecimendone)):
                other_params = prepare_clinicalqueryparams(bc_string)
                write_clinical_processing(scriptfile, index, pj_string, data_category, other_params, outputfile.replace(" ", ""))
                if (clinical_data):
                    clinicaldone = True
                else:
                    biospecimendone = True

def create_Rscript():
    try:
        print "Create R script for downloading data ..."

        #prepare file name for final data set
        outputfile = "_".join(projects) + "_" + "_".join(barcodes)

        #write the Rscript file TODO: open it at the beginning and write one query for each combination
        with open(DATASET_DIR + SCRIPT_NAME, 'w') as f:

            write_preamble(f)

            overall_params = prepare_overallparameters(projects, barcodes, sample_types)
            print overall_params
            iterate_datasets(f, overall_params[0], overall_params[1], outputfile)

        os.chmod(DATASET_DIR + SCRIPT_NAME, stat.S_IRWXU)

        print "... done."
    except:
        print "Stopping execution. An error occurred: ", sys.exc_info()
        exit()

def execute_Rscript():
    try:
        print "Start downloading data..."
        
        subprocess.call('C:\Users\Luise\OneDrive\Arbeit\Hiwi EPIC\Ressources\Preprocessing_Wrapper\Preprocessing_Wrapper\GDCdata\2018-01-19_12-06-54\GDCquery.R')

        print "... done."
    except:
        print "Stopping execution. An error occurred: ", sys.exc_info()
        exit()


start=time.time()
print "-------------------------------------------"
parse_arguments()
print "-------------------------------------------"
check_if_directories_exist()
print "-------------------------------------------"
create_Rscript()
print "-------------------------------------------"
execute_Rscript()
print "-------------------------------------------"
print "Data download finished. Find your final data set files here: " + DATASET_DIR
end=time.time()
