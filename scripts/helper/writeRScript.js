var projects, dataCategories, dataTypes, sampleTypes, workflows, barcodes, filteringThreshold, normalization, featureSelectionThreshold, preprocessingPipeline, featureSelection, discretization, enableDiscretization, otherParams;

//default for values not included in input settings
barcodes = [];
featureSelection = "rowVars";


function saveArguments() {
    projects = $('#project').val();
    dataCategories = $('#data_category').val();
    dataTypes = $('#data_type').val();
    sampleTypes = $('#sample_types').val();
    workflows = $('#workflow').val();
    normalization = $('#normalization').val();
    filteringThreshold = $('#filtering_threshold').val();
    featureSelectionThreshold = $('#fs_threshold').val();
    preprocessingPipeline = $('#preprocessing').val();
}

function getParameterCombinations(arr) {
    return arr.reduce(function (a, b) {
        return a.map(function (x) {
            return b.map(function (y) {
                return x.concat(y);
            })
        }).reduce(function (a, b) { return a.concat(b) }, [])
    }, [[]])
}

function filterEmptyStringsFromArray(arr) {
    return arr.filter(function (n) { return n != "" });
}

function prepareNormalQueryParams(dataType, sampleTypes, workflow, barcodeString) {
    var dataTypeString = "";
    var sampleTypesString = "";
    var workflowString = "";
    var outputName = "";

    //map sample type abbreviation into concrete definitions to use for the query
    if (sampleTypes.length > 1) {
        sampleTypesString = "sample.type = c(";
        for (var i = 0; i < sampleTypes.length; i++) {
            sampleTypesString.concat("\"", sampleTypes[i], "\",");
        }
        //remove last comma       
        sampleTypesString = sampleTypesString.substring(0, sampleTypesString.length - 1);
        sampleTypesString += ")";
    }

    if (dataType) {
        dataTypeString = "data.type = \"" + dataType + "\"";
        outputName += "_" + dataType;
    }

    if (sampleTypes.length > 1) {
        outputName += "_" + sampleTypes.join("_");
    }

    if (workflow) {
        workflowString = "workflow = \"" + workflow + "\"";
        outputName += "_" + workflow;
    }

    otherParams = filterEmptyStringsFromArray([dataTypeString, sampleTypesString, workflowString, barcodeString]).join(", ");
    if (otherParams != "") {
        otherParams = ", " + otherParams;
    }

    return [otherParams, outputName];
}

function prepareClinicalQueryParams(barcodeString) {
    otherParams = [barcodeString].join(", ");
    if (otherParams != "") {
        otherParams = ", " + otherParams;
    }
    return otherParams;
}

function prepareOverallParams(projects, barcodes, sampleTypes) {
    var projectString = "";
    var barcodeString = "";
    projectString += projects.join("\", \"");

    if (barcodes) {
        barcodeString = "barcode = c(\"";
        barcodeString += barcodes.join("\", \"");
        barcodeString += "\")";
    }
    return [projectString, barcodeString];
}

var SCRIPT_TRYSTART = "try({"
var SCRIPT_TRYEND = "})\n"
var SCRIPT_PRINT_START = `cat(\"----------TCGAbiolinks: START PROCESSING NEW{var0}QUERY: {var1}----------\n\")\n`
var SCRIPT_PRINT_END = "cat(\"---------------------TCGAbiolinks: FINISHED QUERY PROCESSING---------------------\n\")"
var SCRIPT_PREAMBLE = `
#! /usr/local/bin/Rscript
library(TCGAbiolinks)
library(dplyr)
library(DT)#package to create HTML pages from the data
library(SummarizedExperiment)#for colData function
library(DESeq2)
library(genefilter)`;

var SCRIPT_DESEQ2 = `#create DESeqDataSet for further processing
dds{var0} <- DESeqDataSet(data{var0}, design = ~ 1)

#filter genes by count threshold (default: 0 counts)
#dds{var0} <- dds{var0}[ rowSums(counts(dds{var0})) > {var1}, ]

#filter out genes
thres <- (ncol(dds{var0}) * {var1}) / 100
dds{var0} <- dds{var0}[ rowSums(counts(dds{var0}) > 0) > thres, ]

#filter out samples
thres <- (nrow(dds{var0}) * {var1}) / 100
dds{var0} <- dds{var0}[ ,colSums(counts(dds{var0}) > 0) > thres]

#logtransform counts and normalize by library size (=cross-sample comparison)
transdds{var0} <- vst(dds{var0}) `;


var SCRIPT_ROWVARS = `topVarGenes <- head(order(rowVars(assay(transdds{var0})), decreasing = TRUE), {var1})
    selected_features{var0} <- transdds{var0}[topVarGenes, ]`;

var SCRIPT_NOFS = `selected_features{var0} <- transdds{var0}`;

function writeProcessingStart(querytype, projects, dataCategory, otherParams) {
    SCRIPT_PRINT_START = SCRIPT_PRINT_START.replace(/{var0}/g, querytype);
    SCRIPT_PRINT_START = SCRIPT_PRINT_START.replace(/{var1}/g, projects + " " + dataCategory + " " + otherParams.replace(/"/g, ""));
    script += SCRIPT_PRINT_START;
    script += "\n";
    script += SCRIPT_TRYSTART;
}

function writePreamble() {
    script += SCRIPT_PREAMBLE;
    script += "\n";
    script += "\n";
}

function writeProcessingEnding() {
    script += SCRIPT_TRYEND;
    script += "\n\n";
    script += SCRIPT_PRINT_END;
    script += "\n";
}

function preprocess(datasetId, pipeline, filteringThreshold, normalization) {
    if (pipeline == "DESeq2") {
        SCRIPT_DESEQ2 = SCRIPT_DESEQ2.replace(/{var0}/g, datasetId);
        SCRIPT_DESEQ2 = SCRIPT_DESEQ2.replace(/{var1}/g, filteringThreshold);
        script += SCRIPT_DESEQ2;
        script += "\n";
    }
}

function applyFeatureSelection(datasetId, approach, featureSelectionThreshold) {
    if (featureSelectionThreshold > "0") {
        if (approach == "rowVars") {
            SCRIPT_ROWVARS = SCRIPT_ROWVARS.replace(/{var0}/g, datasetId);
            SCRIPT_ROWVARS = SCRIPT_ROWVARS.replace(/{var1}/g, featureSelectionThreshold);
            script += SCRIPT_ROWVARS;
            script += "\n";
        }
    }
    else {
        SCRIPT_NOFS = SCRIPT_NOFS.replace(/{var0}/g, datasetId);
        script += SCRIPT_NOFS;
        script += "\n";
    }
}

var SCRIPT_GDCQUERY = `query{var0} <- GDCquery(project= c(\"{var1}\"), data.category = \"{var2}\"{var3}) \n
#download the actual data
GDCdownload(query{var0}) \n`;

var SCRIPT_PREPARE_NORMAL = `#prepares the data for analysis and putting it into an SummarizedExperiment object
data{var0} <- GDCprepare(query{var0}) \n`;

var SCRIPT_PREPARE_CLINICAL = `clinical <- GDCprepare_clinic({var0}, clinical.info = \"patient\")`

function writeGDCQuery(datasetId, projects, dataCategory, otherParams) {
    SCRIPT_GDCQUERY = SCRIPT_GDCQUERY.replace(/{var0}/g, datasetId);
    SCRIPT_GDCQUERY = SCRIPT_GDCQUERY.replace(/{var1}/g, projects);
    SCRIPT_GDCQUERY = SCRIPT_GDCQUERY.replace(/{var2}/g, dataCategory);
    SCRIPT_GDCQUERY = SCRIPT_GDCQUERY.replace(/{var3}/g, otherParams);
    script += SCRIPT_GDCQUERY;
    script += "\n";

    SCRIPT_PREPARE_NORMAL = SCRIPT_PREPARE_NORMAL.replace(/{var0}/g, datasetId);
    script += SCRIPT_PREPARE_NORMAL;
    script += "\n";
}

function writeGDCClinicalQuery(queryname, projects, dataCategory, otherParams) {
    SCRIPT_GDCQUERY = SCRIPT_GDCQUERY.replace(/{var0}/g, queryname);
    SCRIPT_GDCQUERY = SCRIPT_GDCQUERY.replace(/{var1}/g, projects);
    SCRIPT_GDCQUERY = SCRIPT_GDCQUERY.replace(/{var2}/g, dataCategory);
    SCRIPT_GDCQUERY = SCRIPT_GDCQUERY.replace(/{var3}/g, otherParams);
    script += SCRIPT_GDCQUERY;
    script += "\n";

    SCRIPT_PREPARE_CLINICAL = SCRIPT_PREPARE_CLINICAL.replace(/{var0}/g, queryname);
    script += SCRIPT_PREPARE_CLINICAL;
    script += "\n \n";
}

function discretize(datasetId, approach) {
    //fix this
}

var SCRIPT_WRITEOVERALLMETADATAFILE = `metadata{var0} <- colData(transdds{var0})
metadata{var0} <- metadata{var0}[c("project_id", "shortLetterCode")]
transp_metadata{var0} <- t(as.matrix(metadata{var0}))
write.csv(transp_metadata{var0}, file="{var1}_metadata.csv")`;

var SCRIPT_WRITEFINALMETADATAFILE = `sel_meta{var0} <- colData(selected_features{var0})
sel_meta{var0} <- sel_meta{var0}[c("project_id", "shortLetterCode")]
transp_sel_meta{var0} <- t(as.matrix(sel_meta{var0}))
write.csv(transp_sel_meta{var0}, file="{var1}_metadata_topfeatures.csv")`;

function writeMetadataFiles(datasetId, outputFilePath) {
    SCRIPT_WRITEOVERALLMETADATAFILE = SCRIPT_WRITEOVERALLMETADATAFILE.replace(/{var0}/g, datasetId);
    SCRIPT_WRITEOVERALLMETADATAFILE = SCRIPT_WRITEOVERALLMETADATAFILE.replace(/{var1}/g, outputFilePath);
    script += SCRIPT_WRITEOVERALLMETADATAFILE;
    script += "\n";

    SCRIPT_WRITEFINALMETADATAFILE = SCRIPT_WRITEFINALMETADATAFILE.replace(/{var0}/g, datasetId);
    SCRIPT_WRITEFINALMETADATAFILE = SCRIPT_WRITEFINALMETADATAFILE.replace(/{var1}/g, outputFilePath);
    script += SCRIPT_WRITEFINALMETADATAFILE;
    script += "\n";
}

var SCRIPT_WRITEOVERALLFILE = `matrix{var0} <- t(assay(transdds{var0})) #creates the gene expression matrix with genes as rows
write.csv(matrix{var0}, file="{var1}.csv")`;

var SCRIPT_WRITEFINALFILE = `final_matrix{var0} <- t(assay(selected_features{var0})) #creates the gene expression matrix with genes as rows
write.csv(final_matrix{var0}, file="{var1}_topfeatures.csv")`;


function writeDataFiles(datasetId, outputFilePath) {
    SCRIPT_WRITEOVERALLFILE = SCRIPT_WRITEOVERALLFILE.replace(/{var0}/g, datasetId);
    SCRIPT_WRITEOVERALLFILE = SCRIPT_WRITEOVERALLFILE.replace(/{var1}/g, outputFilePath);
    script += SCRIPT_WRITEOVERALLFILE;
    script += "\n";

    SCRIPT_WRITEFINALFILE = SCRIPT_WRITEFINALFILE.replace(/{var0}/g, datasetId);
    SCRIPT_WRITEFINALFILE = SCRIPT_WRITEFINALFILE.replace(/{var1}/g, outputFilePath);
    script += SCRIPT_WRITEFINALFILE;
    script += "\n";
}

function getTimeStamp() {
    var d = new Date();
    var time = d.getFullYear();
    time += "-"
    time += d.getMonth();
    time += "-"
    time += d.getDay();
    time += "_"
    time += d.getHours();
    time += "-"
    time += d.getMinutes();
    time += "-"
    time += d.getSeconds();
    return time;
}

var WORKING_DIR = "./GDCdata/";
var DATASET_DIR = WORKING_DIR + getTimeStamp() + "/";
var SCRIPT_NAME = "GDCquery.R";

function writeNormalProcessing(index, projectString, dataCategory, otherParams, outputFile) {
    var datasetId = index.toString();
    var outputFilePath = DATASET_DIR + outputFile;

    writeProcessingStart("_NORMAL_", projectString, dataCategory, otherParams);
    writeGDCQuery(datasetId, projectString, dataCategory, otherParams);
    preprocess(datasetId, preprocessingPipeline, filteringThreshold, normalization);
    applyFeatureSelection(datasetId, featureSelection, featureSelectionThreshold);

    if (enableDiscretization) {
        discretize(datasetId, discretization);
    }

    writeMetadataFiles(datasetId, outputFilePath);
    writeDataFiles(datasetId, outputFilePath);

    writeProcessingEnding();
}

var SCRIPT_PREP_LOOP = `for(i in c(\"admin\",\"radiation\",\"follow_up\",\"drug\",\"new_tumor_event\")){{
      message(i)
      aux <- GDCprepare_clinic({var0}, clinical.info = i)
      if(is.null(aux)) next
      # add suffix manually if it already exists
      replicated <- which(grep(\"bcr_patient_barcode\",colnames(aux), value = T,invert = T) %in% colnames(clinical))
      colnames(aux)[replicated]<- paste0(colnames(aux)[replicated], \".\",i)
      if(!is.null(aux)) clinical <- merge(clinical, aux, by = \"bcr_patient_barcode\", all = TRUE)
}}`;

var SCRIPT_WRITE_XML = `readr::write_csv(clinical,path = \"{var0}_clinical_from_XML.csv\") # Save the clinical data into a csv file`;

function writeClinicalFile(queryname, outputFilePath) {
    SCRIPT_PREP_LOOP = SCRIPT_PREP_LOOP.replace(/{var0}/g, queryname);
    script += SCRIPT_PREP_LOOP;
    script += "\n \n";

    SCRIPT_WRITE_XML = SCRIPT_WRITE_XML.replace(/{var0}/g, outputFilePath);
    script += SCRIPT_WRITE_XML;
    script += "\n \n";
}

function writeClinicalProcessing(index, projectString, dataCategory, otherParams, outputFile) {
    var queryname = "query" + index.toString();
    var outputFilePath = DATASET_DIR + outputFile

    writeProcessingStart("_CLINICAL_", projectString, dataCategory, otherParams);
    writeGDCClinicalQuery(queryname, projectString, dataCategory, otherParams);
    writeClinicalFile(queryname, outputFilePath);
    writeProcessingEnding();
}

function iterateDatasets(projectString, barcodeString, outputFile) {
    /*create a list with all possible combinations of parameters
    e.g. data.category = "Clinical" "Transcriptome Profiling", data.type = "Clinical" "Gene Expression Quantification"
    results in:
        [("Clinical", "Clinical"),
        ("Clinical", Gene Expression Quantification"),
        ("Transcriptome Profiling", "Clinical"),
        ("Transcriptome Profiling", "Gene Expression Quantification")]
    */

    var parameterCombinations = [];
    parameterCombinations = getParameterCombinations([dataCategories, dataTypes, workflows]);

    var clinicaldone = false;
    var biospecimendone = false;

    //create a query for each combination   
    for (var i = 0; i < parameterCombinations.length; i++) {
        var dataCategory = parameterCombinations[i][0];
        var dataType = parameterCombinations[i][1];
        var workflow = parameterCombinations[i][2];

        clinicalData = (dataCategory == "Clinical");
        biospecimenData = (dataCategory == "Biospecimen");

        //if we have clinical or biospecimen data, we need to apply a different script
        if (!clinicalData && !biospecimenData) {
            var paramsFilename = prepareNormalQueryParams(dataType, sampleTypes, workflow, barcodeString);
            otherParams = paramsFilename[0];
            writeNormalProcessing(i, projectString, dataCategory, otherParams, (outputFile + paramsFilename[1]).replace(" ", ""));
        }
        else {
            /*we only need to create one query for clinical or biospecimen data, otherwise we get duplicate result files 
            because we have multiple combinations with clinical/biospecimen as data_category */
            if ((clinicalData && !clinicaldone) || (biospecimenData && !biospecimendone)) {
                otherParams = prepareClinicalQueryParams(barcodeString);
                writeClinicalProcessing(i, projectString, dataCategory, otherParams, outputFile.replace(" ", ""));
                if (clinicalData) {
                    clinicaldone = true;
                }
                else {
                    biospecimendone = true;
                }
            }
        }
    }
}

var script = "";
var outputFile = "";

function createRScript(processTasks) {
    outputFile = outputFile.concat(projects.join("_"), "_", barcodes.join("_"));
    writePreamble(script);
    var overallParams = prepareOverallParams(projects, barcodes, sampleTypes);
    iterateDatasets(overallParams[0], overallParams[1], outputFile);
}
