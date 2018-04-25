function showAlert(alert) {
    var alertBox = document.getElementById("alertBox");
    alertBox.textContent = alert;
    alertBox.style.backgroundColor = "darkseagreen";

    alertBox.style.display = "inline";

    $("#alertBox").stop(true,true).delay(1000).fadeOut(500);    
}

function showError(error) {
    var errorBox = document.getElementById("alertBox");
    errorBox.textContent = error;
    errorBox.style.backgroundColor = "indianred";
    errorBox.style.display = "inline";

    $("#alertBox").delay(1000).fadeOut(500);
}


// handle Help popup
var openHelpButton = document.getElementById("openHelp-button");
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

// Save Input Settings Button
var saveInputButton = document.getElementById('saveInput-button');
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