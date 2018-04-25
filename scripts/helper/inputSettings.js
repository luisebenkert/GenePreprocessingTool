function setValue(id, value) {
    document.getElementById(id).value = value;
}

// set default values of all input settings
function setDefaultValues() {
    setValue('project', 'TCGA-BRCA');
    setValue('data_category', 'Transcriptome Profiling');
    setValue('data_type', 'Gene Expression Quantification');
    setValue('sample_types', 'TP');
    setValue('workflow', 'HTSeq - Counts');        
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