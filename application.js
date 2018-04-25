////////////MODELER//////////////

var cancelDragDrop = false;

//allow drop on dropzone when drag element is dragged in
function allowDrop(ev) {
    ev.preventDefault();
}

//called when drag is initiated
function drag(ev) {
    if (!cancelDragDrop) {
        ev.dataTransfer.setData("text", ev.target.id);
    }    
}

//called when drag element is dropped on workspace
function onWorkspaceDrop(ev) {
    ev.preventDefault();    

    var dragElementId = ev.dataTransfer.getData("text");
    var dragElement = document.getElementById(dragElementId);    
    var dropzone = getDropzone(ev.target);  

    dropDragElement(dragElement, dropzone);
    removeHoles();
   
}

//called when drag element is dropped back on library
function onLibraryDrop(ev) {
    ev.preventDefault();
    
    var dragElementId = ev.dataTransfer.getData("text");
    var dragElement = document.getElementById(dragElementId);    
    var library = document.getElementById("library");
   
    library.appendChild(dragElement);
    clean(library);
    removeHoles();
}

//cleans useless and whitespace nodes in a nodelist
function clean(node) {
    for (var n = 0; n < node.childElementCount; n++) {
        var child = node.childNodes[n];
        if 
    (
          child.nodeType === 8
          ||
          (child.nodeType === 3 && !/\S/.test(child.nodeValue))
        ) {
            node.removeChild(child);
            n--;
        }
        else if (child.nodeType === 1) {
            clean(child);
        }
    }
}

//remove all holes in the line
function removeHoles() {
    var dropCount = document.getElementById("workspace").childElementCount;

    for (var i = 1; i <= dropCount; i++) {        
        var currentDrop = getElementFromId("drop", i); 
        if (isDropzoneFree(currentDrop)) {            
            var nextDrag = getNextDrag(i + 1, dropCount);            
            if (nextDrag != false) {
                currentDrop.appendChild(nextDrag);
            }            
        }
    }
}

//get the next drag element in line starting from index i
function getNextDrag(i, dropCount) {  
    for (i; i <= dropCount; i++) {        
        var currentDrop = getElementFromId("drop", i);        
        if (!isDropzoneFree(currentDrop)) {           
            var currentDrag = currentDrop.childNodes[0];            
            currentDrop.removeChild(currentDrag);            
            return currentDrag;
        }        
    }  
    return false;
}

//get drag or drop element from id number
function getElementFromId(elementType, idNumber) {
    return document.getElementById(elementType + idNumber);
}

//get the dropzone element when drag element was dropped on occupied dropzone
function getDropzone(dropzone) {
    //check that dragElement was dropped on legit dropzone and not another dragElement
    while (dropzone.id.substring(0, 4) != "drop") {
        dropzone = dropzone.parentNode;
    }
    return dropzone;
}

//returns true if dropzone is free
function isDropzoneFree(dropzone) {
    if (dropzone.childElementCount > 0) {
        return false;
    }
    else {
        return true;
    }
}

//gets id number of drag or drop element name
function getElementIdNumber(e) {   
    return parseInt(e.substring(4, 5));
}

//handles the drop of a drag element
function dropDragElement(draggedElement, activeDropzone) {   
    var firstDrag = false;   

    if (isDropzoneFree(activeDropzone)) {
        activeDropzone.appendChild(draggedElement);       
    }
    else {
        var originDropzone = draggedElement.parentNode;        
        //check if parent is valid dropzone
        if (originDropzone.id.substring(0, 4) != "drop") {
            firstDrag = true;
        }
        else {
            var originDropzoneId = getElementIdNumber(draggedElement.parentNode.id);
        }
        
        var activeDropzoneId = getElementIdNumber(activeDropzone.id);
        var blockingElement = document.getElementById(activeDropzone.childNodes[0].id)
        
        var idDifference = Math.abs(originDropzoneId - activeDropzoneId);        
       

        //move all dragElements down        
        if (idDifference > 1 || firstDrag) {            
            var parentElement = draggedElement.parentNode;
            parentElement.removeChild(draggedElement);
            try {
                moveDragDown(activeDropzone.childNodes[0]);
            }
            catch (error) {
                console.error("Oh no. Something went wrong." + error);
                return false;
            }
            activeDropzone.appendChild(draggedElement);
        }
        //swap both dragElements
        else if (idDifference <= 1) {            
            activeDropzone.removeChild(blockingElement);
            activeDropzone.appendChild(draggedElement);
            originDropzone.appendChild(blockingElement);
        }
    }
}

//moves drag elements down when inserting another inbetween
function moveDragDown(dragElement) {
    var currentDropzone = dragElement.parentNode;
    var nextDropzoneId = parseInt(currentDropzone.id.substring(4, 5)) + 1;
    var nextDropzone = document.getElementById("drop" + nextDropzoneId);
    while (!isDropzoneFree(nextDropzone)) {
        moveDragDown(document.getElementById(nextDropzone.childNodes[0].id));        
    }
    nextDropzone.appendChild(dragElement);
    return true;
}

// clears specified process and returns to start 
function clearProcess() {
    var att = document.getElementById("attributeSetter");
    var attBtn = document.getElementById("attributeButtonBar");
    var cdb = document.getElementById("codeBox");
    var cdbBtn = document.getElementById("codeBoxButtonBar");
    cdb.style.visibility = "hidden";
    cdbBtn.style.visibility = "hidden";
    att.style.visibility = "hidden";
    attBtn.style.visibility = "hidden";

    cancelDragDrop = false;

    var dropCount = document.getElementById("workspace").childElementCount;
    var library = document.getElementById("library");

    // takes all tasks from workspace back to library
    for (var i = 1; i < dropCount; i++) {
        var currentDrop = document.getElementById("drop" + i);
        if (currentDrop.hasChildNodes()) {
            var currentDrag = currentDrop.firstChild;
            library.appendChild(currentDrag);
        }        
    }
}

// checks that given process is correct
function processIsCorrect() {
    var dropCount = document.getElementById("workspace").childElementCount;
    for (var i = 1; i < dropCount; i++) {
        var currentDrop = document.getElementById("drop" + i);
        if (currentDrop.hasChildNodes()) {
            return true;
        }
    }
    return false;
}    

// accepts process and goes onto next step: specification of attributes
function acceptProcess() {
    if (processIsCorrect()) {        
        var att = document.getElementById("attributeSetter");        
        var attBtn = document.getElementById("attributeButtonBar");
        att.style.visibility = "visible";
        attBtn.style.visibility = "visible";
        cancelDragDrop = true;
        setAttributes();
    }    
}

// gets the code and displays it in the codeBox
function generateCode() {
    var cdb = document.getElementById("codeBox");
    var cdbBtn = document.getElementById("codeBoxButtonBar");
    cdb.style.visibility = "visible";
    cdbBtn.style.visibility = "visible";       

    try {
        processTasks = getProcess();
        createRScript(processTasks);
        displayRScript();
        displayButtons();
    }
    catch (err) {
        console.error("Could not parse process" + err);
    }
}

// displays attributes depending on which tasked where specified in the process
function setAttributes() {
    var attHold1 = document.getElementById("attHold1");
    var attHold2 = document.getElementById("attHold2");
    var attHold3 = document.getElementById("attHold3");
    var attHold4 = document.getElementById("attHold4");
    var attHold = document.getElementById("attributeSetter");

    attHold.innerHTML = "";

    var filterAtt = `
                        <div id="filterAtts" class="attribute">
                            <p class="attributeTitle">Filtering Threshold</p>
                            <input type="range" id="filtering_threshold" class="slider" min="0" max="100" step="1" />
                            <input type="text" readonly id="ftValue" class="sliderTextField" maxlength="3">
                        </div>
                    `;
    var featureAtt = `
                        <div id="featureAtts" class="attribute">
                            <p class="attributeTitle">Feature Selection Threshold</p>
                            <input type="range" id="fs_threshold" class="slider" min="0" max="100" step="1" />
                            <input type="text" readonly id="fstValue" class="sliderTextField" maxlength="3">
                        </div>
                    `;
    var normalizeAtt = `
                        <div id="normalizeAtts" class="attribute">
                            <p class="attributeTitle">Normalization</p>
                            <select id="normalization" class="select">
                                <option value="vsd">VSD</option>
                                <option value="rlog">RLog</option>
                            </select>
                        </div>
                    `;
    var preprocessAtt = `
                        <div id="preprocessAtts" class ="attribute">
                            <p class ="attributeTitle">Preprocessing</p>
                            <select id="preprocessing" class ="select">
                                <option value="DESeq2">DESeq2</option>
                                <option value="edgeR">edgeR</option>
                            </select>
                        </div>
                    `;

    var process = getProcess();
    
    for (var i = 0; i < process.length; i++) {
        var task = process[i];
        switch(task) {
            case "Filter":
                attHold.innerHTML += filterAtt;
                break;
            case "Normalize":
                attHold.innerHTML += normalizeAtt;
                break;
            case "Log Transform":
                attHold.innerHTML += preprocessAtt;
                break;
            case "Discretize":
                attHold.innerHTML += featureAtt;
                break;
            default:
                console.error("Weird.");
        } 
    }
}

////////////INPUT SETTINGS//////////////

function setValue(id, value) {
    document.getElementById(id).value = value;
}

// set default values of all input settings
function setDefaultValues() {
    setValue('project', 'TCGA-BRCA');
    setValue('data_category', 'Transcriptome Profiling');
    setValue('data_type');
    setValue('sample_types');
    setValue('workflow');
    setValue('preprocessing', 'DESeq2');
    setValue('normalization', 'vsd');    
    discretizationCheckbox.checked = false;
    ftText.value = ftRange.value = '30';
    fstText.value = fstRange.value = '0';
}

function figurethisout() {
    // Update the value for the Filtering Threshold Range
    var ftRange = document.getElementById('filtering_threshold');
    var ftText = document.getElementById('ftValue');
    var discretizationCheckbox = document.getElementById('enable_discretization');

    ftRange.onchange = function () {
        ftText.value = ftRange.value + " %";
    }

    // Update the value for the Feature Selection Threshold Range
    var fstRange = document.getElementById('fs_threshold');
    var fstText = document.getElementById('fstValue');

    fstRange.onchange = function () {
        fstText.value = fstRange.value + " %";
    }
}


////////////PARSER//////////////

function parseProcess(xml) {
    var parser = new DOMParser();
    var xmlDoc = parser.parseFromString(xml, "text/xml");
    var array = [];

    var process = xmlDoc.getElementsByTagName("bpmn:sequenceFlow");

    if (process.length == 0) {
        alert("Please model a process first.");
        return false;
    }

    var currentSource = "StartEvent";
    var numOfTasks = process.length - 1;
    var i = 0;

    while (array.length <= numOfTasks) {

        // find sequence with currentSource as source    
        if (process[i].getAttribute("sourceRef") == currentSource) {
            currentTarget = process[i].getAttribute("targetRef");

            // Exception handling with EndEvent
            if (currentTarget == "EndEvent") {
                // no tasks
                if (currentSource == "StartEvent") {
                    alert("You must include at least one task or subprocess.");
                    return false;
                }
                    // parsing finished
                else {                    
                    alert("Parsing was successful. See console (F12)");
                    return array;
                }
            }
                // all went well, iterate on
            else {
                array.push(currentTarget);
                currentSource = currentTarget;
                var i = 0;
            }
        }
        else {
            if (i == numOfTasks) {
                alert("Something went wrong");
                return false;
            }
            else {
                i++;
            }
        }
    }

    alert("Something went wrong. Please make sure your process starts with a StartEvent, ends with an EndEvent and includes at least one task or subprocess inbetween.");
    return false;
}

function getProcess() {
    var array = [];
    var dropCount = document.getElementById("workspace").childElementCount;

    for (var i = 1; i < dropCount; i++) {
        var currentDrop = document.getElementById("drop" + i);
        if (currentDrop.hasChildNodes()) {
            var currentDrag = currentDrop.firstChild.textContent;
            array.push(currentDrag);
        }        
    }

    if (array.length < 1) {
        console.error("Please model a process first.");
        return false;
    }

    return array;
}

function download(filename, text) {
    var element = document.createElement('a');
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
    element.setAttribute('download', filename);

    element.style.display = 'none';
    document.body.appendChild(element);

    element.click();

    document.body.removeChild(element);
}

function downloadRScript() {
    download(SCRIPT_NAME, script);
}

function copyScript() {
    var copyText = document.getElementById("textField");
    copyText.select();
    document.execCommand("Copy");
}

function displayButtons() {
    var dwl = document.getElementById("downloadScript-button");
    var cpy = document.getElementById("copyScript-button");
    dwl.style.display = "inline";
    cpy.style.display = "inline";
}

function displayRScript() {
    document.getElementById("textField").innerHTML = script;
}




////////////APPLICATION//////////////

// handle Help popup
var openHelpButton = document.getElementById('openHelp-button');
var helpPopup = document.getElementById('helpPopup');

openHelpButton.onclick = function () {
    helpPopup.style.display = "block";
}

// handle Input Popup
var inputPopup = document.getElementById('inputPopup');

// close popups
window.onclick = function (event) {
    if (event.target == helpPopup) {
        helpPopup.style.display = "none";
    }
    if (event.target == inputPopup) {
        inputPopup.style.display = "none";
    }
}

// dropdown for selecting process
var processSelector = document.getElementById("process-selector");

// Download Script Button
var downloadScriptButton = document.querySelector('#downloadScript-button');
downloadScriptButton.addEventListener('click', function () {
    downloadRScript();
});

// Copy Script Button
var copyScriptButton = document.querySelector('#copyScript-button');
copyScriptButton.addEventListener('click', function () {
    copyScript();
});

// Clear Process Button
var clearProcessButton = document.querySelector('#clearProcess-button');
clearProcessButton.addEventListener('click', function () {
    clearProcess();
});

// Accept Process Button
var acceptProcessButton = document.querySelector('#acceptProcess-button');
acceptProcessButton.addEventListener('click', function () {
    acceptProcess();
});

// Accept Attributes Button
var acceptAttributesButton = document.querySelector('#acceptAttributes-button');
acceptAttributesButton.addEventListener('click', function () {    
    inputPopup.style.display = "block";
    setDefaultValues();   
});

// Process Selector dropdown & button
var selectProcessButton = document.querySelector('#selectProcess-button');
selectProcessButton.addEventListener('click', function () {
    var chosenProcessID = processSelector.options[processSelector.selectedIndex].value;
    if (chosenProcessID > 0) {
        var chosenProcessFileName = dict[chosenProcessID];
        getModel(chosenProcessFileName);
    }
});

// Save Input Settings Button
var saveInputButton = document.getElementById('setInput-button');
saveInputButton.onclick = function () {
    saveArguments();
    inputPopup.style.display = "none";
    generateCode();
}

// Reset Input Settings Button
var resetInputButton = document.getElementById('resetInput-button');
resetInputButton.onclick = function () {
    setDefaultValues();
}







///////////////WRITE SCRIPT////////////////////

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

