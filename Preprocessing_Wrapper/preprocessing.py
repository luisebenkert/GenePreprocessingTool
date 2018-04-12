#!/usr/bin/python

#Script for including preprocessing, i.e. filtering, normalization and transformation. Uses state-of-the-art tools, e.g. DESeq2.

import sys
import os, stat, subprocess

SCRIPT_DESEQ2 ="""#create DESeqDataSet for further processing
dds{0} <- DESeqDataSet(data{0}, design = ~ 1)

#filter genes by count threshold (default: 0 counts)
#dds{0} <- dds{0}[ rowSums(counts(dds{0})) > {1}, ]

#filter out genes
thres <- (ncol(dds{0}) * {1}) / 100
dds{0} <- dds{0}[ rowSums(counts(dds{0}) > 0) > thres, ]

#filter out samples
thres <- (nrow(dds{0}) * {1}) / 100
dds{0} <- dds{0}[ ,colSums(counts(dds{0}) > 0) > thres]

#logtransform counts and normalize by library size (=cross-sample comparison)
transdds{0} <- vst(dds{0})""" #todo: make threshold for filtering a param and choose between rlog and vst

def preprocess(scriptfile, datasetid, pipeline, filtering_threshold, normalization):


    if pipeline == "DESeq2":
        scriptfile.write(SCRIPT_DESEQ2.format(datasetid, filtering_threshold, normalization) + "\n")
