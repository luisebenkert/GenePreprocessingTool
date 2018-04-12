
///////////////WRITE SCRIPT////////////////////

var projects, dataCategories, dataTypes, sampleTypes, workflows, barcodes, filteringThreshold, normalization, featureSelectionThreshold, preprocessingPipeline, featureSelection, discretization, enableDiscretization, otherParams;

function saveArguments() {
    projects = $('#project').val();
    dataCategories = $('#data_category').val();
    dataTypes = $('#data_type').val();
    sampleTypes = $('#sample_types').val();
    workflows = $('#workflow').val();
    prep = $('#project').val();
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

function prepareNormalQueryParams(dataType, sampleTypes, workflow, barcodeString) {
    var dataTypeString = "";
    var sampleTypesString = "";
    var workflowString = "";
    var outputName = "";

    //map sample type abbreviation into concrete definitions to use for the query
    if (sampleTypes != [""]) {
        sampleTypesString = "sample.type = c(";
        for (var i = 0; i < sampleTypes.length; i++) {
            sampleTypesString.concat("\"", sampleTypes[i], "\",");
        }
        //remove last comma
        sampleTypesString = sampleTypesString.substring(0, sampleTypesString.length - 1);
    }

    if (dataType) {
        dataTypeString = "data.type = \"" + dataType + "\"";
        outputName += "_" + dataType;
    }

    if (sampleTypes != [""]) {
        outputName += "_" + sampleTypes.join("_");
    }

    if (workflow) {
        workflowString = "workflow = \"" + workflow + "\"";
        outputName += "_" + workflow;
    }

    otherParams = [dataTypeString, sampleTypesString, workflowString, barcodeString].join(", ");
    if (otherParams != "") {
        otherParams = ", " + otherParams;
    }

    return (otherParams, outputName);
}

function prepareClinicalQueryParams(barcodeString) {
    otherParams = [barcodeString].join(", ");
    if (otherParams != "") {
        otherParams = ", " + otherParams;
    }
    return otherParams;
}

const SCRIPT_TRYSTART = "try({"
const SCRIPT_TRYEND = "})\n"
const SCRIPT_PRINT_START = "cat(\"----------TCGAbiolinks: START PROCESSING NEW ${start0} QUERY: ${start1}----------\n\")"
const SCRIPT_PRINT_END = "cat(\"---------------------TCGAbiolinks: FINISHED QUERY PROCESSING---------------------\n\")"
const SCRIPT_PREAMBLE = `
    #! /usr/local/bin/Rscript
    library(TCGAbiolinks)
    library(dplyr)
    library(DT)#package to create HTML pages from the data
    library(SummarizedExperiment)#for colData function
    library(DESeq2)
    library(genefilter)
    `;

const SCRIPT_DESEQ2 = `
    ""#create DESeqDataSet for further processing
    dds${deseq0} <-DESeqDataSet(data${deseq0}, design = ~1)

    #filter genes by count threshold (default: 0 counts)
    #dds${deseq0} <-dds${deseq0}[rowSums(counts(dds${deseq0})) > ${deseq1}, ]

    #filter out genes
    thres <-(ncol(dds${deseq0}) * ${deseq1}) / 100
    dds${deseq0} <-dds${deseq0}[rowSums(counts(dds${deseq0}) > 0) > thres, ]

    #filter out samples
    thres <-(nrow(dds${deseq0}) * ${deseq1}) / 100
    dds${deseq0} <-dds${deseq0}[, colSums(counts(dds${deseq0}) > 0) > thres]

    #logtransform counts and normalize by library size (=cross-sample comparison)
    transdds${deseq0} <-vst(dds${deseq0}) ""`;

const SCRIPT_ROWVARS = `
    ""topVarGenes <- head(order(rowVars(assay(transdds${rowvars0})), decreasing = TRUE), ${rowvars1})
    selected_features${rowvars0} <- transdds${rowvars0}[topVarGenes, ]""`;

const SCRIPT_NOFS = `
    ""selected_features${nofs0} <- transdds${nofs0}""`;


function writeProcessingStart(script, querytype, projects, dataCategory, otherParams) {
    var start0 = querytype;
    var start1 = projects + " " + dataCategory + " " + otherParams.replace("\"", "");
    script += SCRIPT_PRINT_START;
    script += "\n";
    script += SCRIPT_TRYSTART;
}

function writePreamble(script) {
    script += SCRIPT_PREAMBLE;
    script += "\n";
    script += "\n";
}

function writeProcessingEnding(script) {
    script += SCRIPT_TRYEND;
    script += "\n\n";
    script += SCRIPT_PRINT_END;
    script += "\n";
}

function preprocess(script, datasetId, pipeline, filteringThreshold, normalization) {
    if (pipeline == "DESeq2") {
        var deseq0 = datasetId;
        var deseq1 = filteringThreshold;
        script += SCRIPT_DESEQ2;
        script += "\n";
    }
}

function applyFeatureSelection(script, datasetId, approach, featureSelectionThreshold) {
    if (featureSelectionThreshold > "0") {
        if (approach == "rowVars") {
            var rowvars0 = datasetId;
            var rowvars1 = featureSelectionThreshold;
            script += SCRIPT_ROWVARS;
            script += "\n";
        }
    }
    else {
        var nofs0 = datasetId;
        script += SCRIPT_NOFS;
        script += "\n";
    }
}

const SCRIPT_GDCQUERY = `
    ""query${gdc0} <- GDCquery(project= c(\"${gdc1}\"), data.category = \"${gdc2}\"${gdc3})

    #download the actual data
    GDCdownload(query${gdc0})""`;

const SCRIPT_PREPARE_NORMAL = `
    ""#prepares the data for analysis and putting it into an SummarizedExperiment object
    data${script0} <- GDCprepare(query${script0})""`;

const SCRIPT_PREPARE_CLINICAL = "clinical <- GDCprepare_clinic(${scriptClinical0}, clinical.info = \"patient\")"

function prepareNormalQueryParams(dataType, sampleType, workflow, barcodeString) {

}

function writeGDCQuery() {

}

function writeGDCClinicalQuery() {

}

function discretize(script, datasetId, approach) {
    //fix this
}

function writeNormalProcessing(script, index, projectString, dataCategory, otherParams) {
    var datasetId = index.toString();

    writeProcessingStart(script, "_NORMAL_", projectString, dataCategory, otherParams);
    writeGDCQuery(script, datasetId, projectString, dataCategory, otherParams);
    preprocess(script, datasetId, preprocessingPipeline, filteringThreshold, normalization);
    applyFeatureSelection(script, datasetId, featureSelection, featureSelectionThreshold);

    if (enableDiscretization) {
        discretize(script, datasetId, discretization);
    }

    writeProcessingEnding(script);
}

function writeClinicalProcessing(script, index, projectString, dataCategory, otherParams) {
    var queryname = "query" + index.toString();
    writeProcessingStart(script, "_CLINICAL_", projectString, dataCategory, otherParams);
    writeGDCClinicalQuery(script, queryname, projectString, dataCategory, otherParams);
    writeProcessingEnding(script);
}

function iterateDatasets(script, projectString, barcodeString) {
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
            paramsFilename = prepareNormalQueryParams(dataType, sampleTypes, workflow, barcodeString);
            otherParams = paramsFilename[0];
            writeNormalProcessing(script, i, projectString, dataCategory, otherParams, (output + paramsFilename[1]).replace(" ", ""));
        }
        else {
            /*we only need to create one query for clinical or biospecimen data, otherwise we get duplicate result files 
            because we have multiple combinations with clinical/biospecimen as data_category */
            if ((clinicalData && !clinicaldone) || (biospecimenData && !biospecimendone)) {
                otherParams = prepareClinicalQueryParams(barcodeString);
                writeClinicalProcessing(script, i, projectString, dataCategory, otherParams, output.replace(" ", ""));
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

var finalScript = "";
function createRScript(processTasks) {
    writePreamble(finalScript);
    var overallParams = prepareOverallParams(projects, barcodes, sampleTypes);
    iterateDatasets(finalScript, overallParams[0], overallParams[1]);
    console.log(finalScript);
}
