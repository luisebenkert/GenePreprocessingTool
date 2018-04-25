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
                showError("Something went wrong. Please try again.")
                console.error(error);                
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
    for (var i = 1; i <= dropCount; i++) {
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
    for (var i = 1; i <= dropCount; i++) {
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
    else {
        showError("Please model a process first.");
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
        showAlert("Script was successfully parsed");
    }
    catch (err) {
        showError("Could not parse process. Please try again.");
        console.error(err);
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
        switch (task) {
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
        }
    }
}