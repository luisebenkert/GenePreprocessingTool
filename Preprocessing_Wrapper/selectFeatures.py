#!/usr/bin/python

#Script for including feature selection in data processing.

import sys
import os, stat, subprocess

SCRIPT_ROWVARS ="""topVarGenes <- head(order(rowVars(assay(transdds{0})), decreasing = TRUE), {1})
selected_features{0} <- transdds{0}[topVarGenes, ]"""
SCRIPT_NOFS ="""selected_features{0} <- transdds{0}"""

def apply_featureselection(scriptfile, datasetid, approach, fs_threshold):

    if fs_threshold > "0": #only do feature selection if a concrete number of features is specified
        if approach == "rowVars":
            scriptfile.write(SCRIPT_ROWVARS.format(datasetid, fs_threshold) + "\n")
    else:
        scriptfile.write(SCRIPT_NOFS.format(datasetid) + "\n")
