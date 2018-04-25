function getProcess() {    
    var array = [];
    var dropCount = document.getElementById("workspace").childElementCount;

    for (var i = 1; i <= dropCount; i++) {
        var currentDrop = document.getElementById("drop" + i);
        if (currentDrop.hasChildNodes()) {
            var currentDrag = currentDrop.firstChild.textContent;
            array.push(currentDrag);
        }
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
    showAlert("File is ready for download.");
}

function downloadRScript() {
    download(SCRIPT_NAME, script);
}

function copyScript() {
    var copyText = document.getElementById("textField");
    copyText.select();
    document.execCommand("Copy");
    showAlert("Script was copied to your clipboard");
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
