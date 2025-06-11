// aflowlib_webapp_entry.js
//
// author: Bob Hanson
// called from aflow_web_interface2.cpp via
// http://aflowlib.duke.edu/users/jmolers/test/entry.php?id=aflow:137d3cb39fa592d3

//BH20160521 12:57:40 PM
//BH20160522 10:03:22 PM finalized
//CO20170623 12:02:22 PM minor changes
//PC20180723 Add Symmetry Operations

//var jsmolDir = "http://aflowlib.duke.edu/users/jmolers/test/jsmol";

document.write("<div id='jmol'></div>");

//Jmol.db._DirectDatabaseCalls["aflowlib.duke.edu"] = "%URL"; // not indicated in JSmolCore.js
//Jmol.db._DirectDatabaseCalls["aflow.duke.edu"] = "%URL"; // not indicated in JSmolCore.js
//Jmol.db._DirectDatabaseCalls["aflowlib.duke.edu"] = "%URL"; // not indicated in JSmolCore.js

//function testfunction() {
//  var all = document.getElementById('symop-container').querySelectorAll('*');
// all.find('input[type=checkbox]:checked').removeAttr('checked');
//console.log('hello worl!');
//}

if (!self.AFLOW) {
  AFLOW = {};
  AFLOW.version = "31028";
  AFLOW.url_WEB =
    "http://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/BCC/Ag3Au1Se2_ICSD_171959";
  AFLOW.label = "Ag3Au1Se2_ICSD_171959";
  AFLOW.spaceGroupNo = 214;
  AFLOW.spaceGroupName = "I4_{1}32";
  AFLOW.cif_sconv = [10.13, 10.13, 10.13, 90.0, 90.0, 90.0];
  AFLOW.cif = [8.773, 8.773, 8.773, 109.471, 109.471, 109.471];
  AFLOW.cif_sprim = [8.773, 8.773, 8.773, 109.471, 109.471, 109.471];
  AFLOW.baderUnitcell = [
    -5.065,
    5.065,
    5.065,
    5.065,
    -5.065,
    5.065,
    5.065,
    5.065,
    -5.065
  ];
  AFLOW.baderVSpecies = ["Ag", "Au", "Se"];
  AFLOW.jsmolDir = "../test/jsmol";
}
AFLOW.url = AFLOW.url_WEB + "/" + AFLOW.label;

// top left message about the space group
AFLOW.cif2oss = function(params) {
  return (
    " font echo 14;color echo white; set echo 3% 97%; echo AFLOW.org consortium (AFLOW v" +
    AFLOW.version +
    ") | entry=" +
    AFLOW.label +
    "  |  | " +
    (AFLOW.spaceGroupNo
      ? "  Spacegroup = " +
        AFLOW.spaceGroupName +
        " (#" +
        AFLOW.spaceGroupNo +
        ")   |"
      : "") +
    " a=" +
    params[0] +
    "\u212B, b=" +
    params[1] +
    "\u212B, c=" +
    params[2] +
    "\u212B    |" +
    " \u03B1=" +
    params[3] +
    "\u00B0, \u03B2=" +
    params[4] +
    "\u00B0, \u03B3=" +
    params[5] +
    "\u00B0 ; "
  );
};

// isosurface creation
AFLOW.colors = [
  "red",
  "green",
  "yellow",
  "blue",
  "orange",
  "white",
  "purple",
  "brown",
  "pink"
];
AFLOW.iso2oss = function(element, cutoff, index) {
  return (
    "  ISOSURFACE  " +
    element +
    " '" +
    AFLOW.url +
    "_Bader_" +
    cutoff +
    "_" +
    element +
    ".jvxl';" +
    "isosurface mesh; color isosurface " +
    AFLOW.colors[index] +
    " translucent;" +
    "VAR charges = load('" +
    AFLOW.url +
    "_abader.out" +
    "').split('=')[2].split('(')[1].split(','); {*}.label = charges; label %[label];"
  );
};

// load command, with optional parameters
//PC -- add a save command to reset the last loaded structure in symmetry part
AFLOW.load = function(root, loadparams) {
  var key = root + loadparams;
  return (
    "if (key != '" +
    key +
    "') {load '" +
    AFLOW.url +
    root +
    ".cif' " +
    (loadparams ? loadparams : "packed") +
    ";" +
    AFLOW.cif2oss(AFLOW["cif" + root]) +
    ";key='" +
    key +
    "'}; save state restorestate; " //CO20190430 - missing semicolon
  );
};

// JSmol info
AFLOW.Info = {
  width: Math.floor(window.innerWidth * 0.6),
  height: Math.floor(window.innerWidth * 0.6),
  debug: false,
  color: "black",
  addSelectionOptions: false,
  serverURL: "https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php",
  use: "HTML5",
  j2sPath: AFLOW.jsmolDir + "/j2s",
  readyFunction: function(applet) {
    Jmol._getElement(applet, "appletdiv").style.border = "1px solid blue";
  },
  script:
    "set zoomlarge false; set platformspeed 3;set antialiasDisplay; frank off; set showUnitCellinfo false; " +
    AFLOW.load("_sconv"),
  //disableJ2SLoadMonitor: true,
  disableInitialConsole: true
};

// thumb wheel widget
AFLOW.input = function(id, val) {
  return (
    "<input class='dim' id='" +
    id +
    "' type='number' value='" +
    (val || 2) +
    "' style='width:60px' >"
  );
};

//symetry operations -- PC20180723
//

function printDigit(x) {
  return Number.parseFloat(x).toFixed(2);
}

function multiplyMinusOne(arr) {
  var newArr = arr.map(function(element) {
    return element * -1;
  });
  return newArr;
}

if (AFLOW.sym2json) {
  var numSymOp = AFLOW.sym2json.length;
  AFLOW.symOp = [];
  var namesSymOp = [];
  var mAxis = [];
  for (var i = 0; i < numSymOp; i++) {
    var axis = AFLOW.sym2json[i]["axis"];
    for (var j = 0; j < axis.length; j++) {
      if (!Number.isInteger(axis[j])) {
        axis[j] = printDigit(axis[j]);
        if (axis[j] == 0.0 || axis[j] == -0.0) {
          axis[j] = 0;
        }
      }
    }
    mAxis[0] = -axis[0];
    mAxis[1] = -axis[1];
    mAxis[2] = -axis[2];
    var angle = AFLOW.sym2json[i]["angle"];
    var checkAngle = "";
    if (angle > 180) {
      checkAngle += ", angle: " + angle;
    }
    var adname = "";
    var uc =
      "" +
      AFLOW.sym2json[i]["Uc"][0] +
      "," +
      AFLOW.sym2json[i]["Uc"][1] +
      "," +
      AFLOW.sym2json[i]["Uc"][2] +
      "";
    var mUc =
      "" +
      multiplyMinusOne(AFLOW.sym2json[i]["Uc"][0], -1) +
      "," +
      multiplyMinusOne(AFLOW.sym2json[i]["Uc"][1], -1) +
      "," +
      multiplyMinusOne(AFLOW.sym2json[i]["Uc"][2], -1) +
      "";
    var script = "select 1.1; rotateSelected [[" + uc + "]] 80;";
    var rotoinversionsScript =
      "select 1.1; rotateSelected [[" +
      mUc +
      "]] 80; invertSelected POINT {0 0 0};";
    var drawAxis =
      "draw axis" +
      i +
      " fixed {" +
      mAxis +
      "/1} {" +
      axis +
      "/1} 200; color $axis" +
      i +
      " yellow;";
    var removeAxis = "hide $axis" + i + ";";
    var drawPlane =
      "isosurface ID plane" +
      i +
      " PLANE {" +
      axis +
      ",0}; color $plane" +
      i +
      " yellow translucent 0.4;";
    var removePlane = "hide $plane" + i + ";";
    var reflectionScript =
      "select 1.1; color atoms translucent 0.8; delay 0.2; select 1.1; color atoms translucent 0;select 1.1; invertSelected PLANE {" +
      axis +
      ",0};";
    var reflectionName = "reflection,  plane: (" + axis + ",0)";
    var inversionName = "Inversion -I";
    var inversionScript = "select 1.1; invertSelected POINT {0 0 0};";
    var ftau = AFLOW.sym2json[i]["ftau"];
    if (ftau[0] > 0.001 || ftau[1] > 0.001 || ftau[2] > 0.001) {
      script += "delay 0.5; translateSelected {" + ftau + "/1};";
      inversionScript += "delay 0.5; translateSelected {" + ftau + "/1};";
      reflectionScript += "delay 0.5; translateSelected {" + ftau + "/1};";
      rotoinversionsScript += "delay 0.5; translateSelected {" + ftau + "/1};";
      adname += " +translation";
      inversionName += " +translation";
      reflectionName += " +translation";
    }
    //numerotate the similar operations
    var name =
      AFLOW.sym2json[i]["type"] +
      " " +
      AFLOW.sym2json[i]["Schoenflies"] +
      adname;
    namesSymOp.push(name);
    var count = 0;
    for (var j = 0; j < namesSymOp.length; j++) {
      if (namesSymOp[j] == name) {
        count += 1;
      }
    }
    name += " " + count;
    AFLOW.symOp[i] = {
      name: name,
      reflectionName: reflectionName,
      inversionName: inversionName,
      script: script,
      reflectionScript: reflectionScript,
      inversionScript: inversionScript,
      rotoinversionsScript: rotoinversionsScript,
      drawAxis: drawAxis,
      removeAxis: removeAxis,
      axis: axis,
      mAxis: mAxis,
      angle: angle,
      drawPlane: drawPlane,
      removePlane: removePlane,
      checkAngle: checkAngle
    };
  }
  var rotations = [];
  var rotoinversions = [];
  for (var i = 0; i < AFLOW.symOp.length; i++) {
    if (AFLOW.symOp[i].name.includes("rotation")) {
      rotations.push(AFLOW.symOp[i].name);
    } else if (AFLOW.symOp[i].name.includes("inversion")) {
      rotoinversions.push(AFLOW.symOp[i].name);
    }
  }
  rotations.sort();
  rotoinversions.sort();
  var rotationsSorted = [];
  var rotoinversionsSorted = [];
  var allOperations = [];
  for (var i = 0; i < rotations.length; i++) {
    for (var j = 0; j < AFLOW.symOp.length; j++) {
      if (AFLOW.symOp[j].name == rotations[i]) {
        AFLOW.symOp[j].name = AFLOW.symOp[j].name.substring(
          0,
          AFLOW.symOp[j].name.length - 1
        );
        rotationsSorted.push(AFLOW.symOp[j]);
        allOperations.push(AFLOW.symOp[j]);
      }
    }
  }
  for (var i = 0; i < rotoinversions.length; i++) {
    for (var j = 0; j < AFLOW.symOp.length; j++) {
      if (AFLOW.symOp[j].name == rotoinversions[i]) {
        AFLOW.symOp[j].name = AFLOW.symOp[j].name.substring(
          0,
          AFLOW.symOp[j].name.length - 1
        );
        rotoinversionsSorted.push(AFLOW.symOp[j]);
        allOperations.push(AFLOW.symOp[j]);
      }
    }
  }
}

$(document).ready(function() {
  // create header, jmol html, options, and save options, then put them all together at the end
  var header =
    "<br /><br /><b>Space Group</b>: " +
    (AFLOW.spaceGroupNo
      ? AFLOW.spaceGroupName + " (#" + AFLOW.spaceGroupNo + ")"
      : "N/A");

  var jmol = Jmol.getAppletHtml("jmolApplet0", AFLOW.Info);

  var options = "";

  Jmol.setButtonCss(null, "style='width:180px'");

  options +=
    "<b>Relaxed Structure:</b><br /><br />" +
    Jmol.jmolButton(jmolApplet0, AFLOW.load(""), "As calculated") +
    "<br />" +
    Jmol.jmolButton(
      jmolApplet0,
      AFLOW.load("_sconv"),
      "Standard conventional"
    ) +
    "<a href=http://www.sciencedirect.com/science/article/pii/S0927025610002697>[info]</a>" +
    "<br />" +
    Jmol.jmolButton(jmolApplet0, AFLOW.load("_sprim"), "Standard primitive") +
    "<a href=http://www.sciencedirect.com/science/article/pii/S0927025610002697>[info]</a>";

  Jmol.setButtonCss(null, "style='width:110px'");

  options +=
    "<br /><br /><b>Supercell:</b><br /><br />" +
    Jmol.jmolButton(jmolApplet0, "key ='';load '' {2 2 2} packed", "2x2x2") +
    Jmol.jmolButton(jmolApplet0, "key ='';load '' fill 20", "20&#8491; box") +
    "<p>" +
    AFLOW.input("dim_1") +
    " X " +
    AFLOW.input("dim_2") +
    " X " +
    AFLOW.input("dim_3") +
    "</p>" +
    "<input type='button' id='build_button' value='Build' style='width:160px'>" +
    Jmol.jmolButton(jmolApplet0, AFLOW.load("_sconv") + ";reset;", "RESET");

  Jmol.setButtonCss(null, "style='width:110px'");

  options +=
    "<br /><br /><b>Visualization:</b><br /><br />" +
    Jmol.jmolButton(
      jmolApplet0,
      "spacefill only;spacefill 23%;wireframe 0.15",
      "Ball & Stick"
    ) +
    Jmol.jmolButton(
      jmolApplet0,
      "spacefill #alt:SETTING van der Waals Spheres",
      "Spacefill"
    ) +
    "<br />" +
    Jmol.jmolCheckbox(jmolApplet0, "spin on", "spin off", "Rotation") +
    Jmol.jmolCheckbox(jmolApplet0, "label  %a ", "labels off", "Labels") +
    Jmol.jmolCheckbox(
      jmolApplet0,
      "background white",
      "background black",
      "Background"
    );
  Jmol.setButtonCss(null, "style='width:30px'");

  options +=
    "<br />axis: " +
    Jmol.jmolButton(jmolApplet0, "moveto axis a", "a") +
    Jmol.jmolButton(jmolApplet0, "moveto axis b", "b") +
    Jmol.jmolButton(jmolApplet0, "moveto axis c", "c");
  Jmol.setButtonCss(null, "style='width:110px'");
  options += Jmol.jmolButton(jmolApplet0, "reset", "RESET");

  options +=
    "<br /><br /><b>Crystallographic Planes:</b><br />" +
    "<p> h:" +
    AFLOW.input("plane_1") +
    " k:" +
    AFLOW.input("plane_2") +
    " l:" +
    AFLOW.input("plane_3") +
    "</p>" +
    "<input type='button' id='plane_button' value='Show plane' style='width:160px'>" +
    Jmol.jmolButton(jmolApplet0, "isosurface off", "RESET");

  //BEGIN BADER ISOSURFACES
  var bader = "";
  if (AFLOW.baderUnitcell) {
    bader += "<br /><br /><b>Bader Isosurfaces:</b><br />";

    Jmol.setButtonCss(null, "style='width:50px'");
    var cutoffs = [20, 30, 40, 50];
    var nSpecies = AFLOW.baderVSpecies.length;
    for (i = 0; i < cutoffs.length; i++) {
      var cutoff = cutoffs[i];
      bader += "<br />Cutoff = 0." + cutoff + "<br />";
      var unitcell = "packed UNITCELL[" + AFLOW.baderUnitcell.join(",") + "]";
      for (var j = 0; j < nSpecies; j++) {
        var spec = AFLOW.baderVSpecies[j];
        bader += Jmol.jmolButton(
          jmolApplet0,
          AFLOW.load("", unitcell) +
            ";isosurface delete;" +
            AFLOW.iso2oss(spec, cutoff, j),
          spec
        );
      }
      var s = "isosurface delete;" + AFLOW.load("", unitcell);
      for (var j = 0; j < nSpecies; j++)
        s += AFLOW.iso2oss(AFLOW.baderVSpecies[j], cutoff, j);
      bader += Jmol.jmolButton(jmolApplet0, s, "All");
    }
    bader += "<br />";
  }
  //END BADER ISOSURFACES

  //BEGIN SYMMETRY -- PC20180723
  //
  var symProperties = "";
  if (AFLOW.sym2json) {
    Jmol.setButtonCss(null, "style='width:200px'");
    symProperties += "<br /><br /><b>Symmetry Operations: </b><br />";
    Jmol.setButtonCss(null, "style='width:60px'");
    symProperties +=
      "<p>" +
      numSymOp +
      " operations calculated with respect to the 'As Calculated' material.</p>";
    symProperties += Jmol.jmolCheckbox(
      jmolApplet0,
      "load append '' {444 666 -1}; select 2.1; color translucent 0.9 [224 224 224];  frame all; center; ",
      "zap 2.1;",
      "Background supercell (for better visualization)"
    );
    symProperties += "<div id='symop-container' class='symop-container'>";
    symProperties += "<span style='font-size:small'> Identity";
    symProperties += "</span><br />";
    for (var k = 0; k < rotationsSorted.length; k++) {
      symProperties +=
        "<span style='font-size:small'>" +
        rotationsSorted[k].name +
        rotationsSorted[k].checkAngle +
        ",   axis: (" +
        rotationsSorted[k].axis +
        "), ";
      symProperties += Jmol.jmolCheckbox(
        jmolApplet0,
        rotationsSorted[k].drawAxis,
        rotationsSorted[k].removeAxis,
        "axis"
      );
      symProperties += Jmol.jmolButton(
        jmolApplet0,
        rotationsSorted[k].script,
        "apply"
      );
      symProperties += "</span><br />";
    }
    for (var k = 0; k < rotoinversionsSorted.length; k++) {
      if (rotoinversionsSorted[k].name.includes("-I")) {
        symProperties += "<span style='font-size:small'>";
        symProperties += rotoinversionsSorted[k].inversionName;
        symProperties += Jmol.jmolCheckbox(
          jmolApplet0,
          "polyhedra ID pol {0 0 0} TO [{0 0 0.25/1}, {-0.25 0 0/1}, {0 0.25 0/1}, {0 0 -0.25/1}, {0.25 0 0/1}, {0 -0.25 0/1}] ; color $pol yellow;",
          "delete $pol",
          "point"
        );
        symProperties += Jmol.jmolButton(
          jmolApplet0,
          rotoinversionsSorted[k].inversionScript,
          "apply"
        );
        symProperties += "</span><br />";
      } else if (rotoinversionsSorted[k].name.includes("S2")) {
        symProperties += "<span style='font-size:small'>";
        symProperties += rotoinversionsSorted[k].reflectionName;
        symProperties += Jmol.jmolCheckbox(
          jmolApplet0,
          rotoinversionsSorted[k].drawPlane,
          rotoinversionsSorted[k].removePlane,
          "plane"
        );
        symProperties += Jmol.jmolButton(
          jmolApplet0,
          rotoinversionsSorted[k].reflectionScript,
          "apply"
        );
        symProperties += "</span><br />";
      } else if (rotoinversionsSorted[k].name.includes("rotoinversion s")) {
        symProperties += "<span style='font-size:small'>";
        symProperties += rotoinversionsSorted[k].reflectionName;
        symProperties += Jmol.jmolCheckbox(
          jmolApplet0,
          rotoinversionsSorted[k].drawPlane,
          rotoinversionsSorted[k].removePlane,
          "plane"
        );
        symProperties += Jmol.jmolButton(
          jmolApplet0,
          rotoinversionsSorted[k].reflectionScript,
          "apply"
        );
        symProperties += "</span><br />";
      }
    }
    for (var k = 0; k < rotoinversionsSorted.length; k++) {
      if (
        rotoinversionsSorted[k].name.includes("S3") ||
        rotoinversionsSorted[k].name.includes("S4") ||
        rotoinversionsSorted[k].name.includes("S6")
      ) {
        symProperties += "<span style='font-size:small'>";
        symProperties += rotoinversionsSorted[k].name;
        symProperties +=
          rotoinversionsSorted[k].checkAngle +
          ",   axis: (" +
          rotoinversionsSorted[k].axis +
          "), ";
        symProperties += Jmol.jmolCheckbox(
          jmolApplet0,
          rotoinversionsSorted[k].drawAxis,
          rotoinversionsSorted[k].removeAxis,
          "axis"
        );
        symProperties += Jmol.jmolCheckbox(
          jmolApplet0,
          "polyhedra ID pol {0 0 0} TO [{0 0 0.25/1}, {-0.25 0 0/1}, {0 0.25 0/1}, {0 0 -0.25/1}, {0.25 0 0/1}, {0 -0.25 0/1}] ; color $pol yellow;",
          "delete $pol",
          "point"
        );
        symProperties += Jmol.jmolButton(
          jmolApplet0,
          rotoinversionsSorted[k].rotoinversionsScript,
          "apply"
        );
        symProperties += "</span><br />";
      }
    }
    symProperties += "</div>";
    Jmol.setButtonCss(null, "style='width:110px'");
    symProperties += "<form id='test'>";
    symProperties += Jmol.jmolButton(
      jmolApplet0,
      "restore state restorestate;",
      "RESET",
      (id = "symop-reset")
    );
  }
  // END SYMMETRY -- PC20180723

  Jmol.setButtonCss(null, "style='width:140px'");

  var saveoptions =
    "<b>Save:</b>" +
    Jmol.jmolButton(jmolApplet0, "write FILE ?", "CIF FILE") +
    Jmol.jmolButton(jmolApplet0, "write STATE ?.spt", "STATE") +
    Jmol.jmolButton(jmolApplet0, "write IMAGE ?.jpg", "JPG") +
    Jmol.jmolButton(jmolApplet0, "write IMAGE ?.png", "PNG") +
    Jmol.jmolButton(jmolApplet0, "write PNGJ ?.png", "PNG+Jmol");

  var html =
    header +
    "<table><tr><td align=center valign=top>" +
    jmol +
    "</td>" +
    "<td><form>" +
    options +
    symProperties +
    bader +
    "</form></td></tr>" +
    "<tr><td align=center>" +
    saveoptions +
    "</td></tr></table>";

  $("#jmol").html(html);

  var resetSym = document.getElementById("symop-reset");
  resetSym.addEventListener("click", function() {
    console.log("enter function");
    var checkboxes = document.getElementsByTagName("input");
    console.log("checkboxes");
    for (var i = 0; i < checkboxes.length; i++) {
      checkboxes[i].checked = false;
    }
  });

  // set non-Jmol button click events
  $("#build_button").click(function() {
    var scriptCommand =
      "key='';load '' {" +
      $("#dim_1").val() +
      " " +
      $("#dim_2").val() +
      " " +
      $("#dim_3").val() +
      "} packed;";
    Jmol.script(jmolApplet0, scriptCommand);
  });
  $("#plane_button").click(function() {
    //BH removed: key='';load '' {444 555 -1} packed;  -- why reload?
    var scriptCommand =
      "isosurface ID 'hklplane' hkl {" +
      $("#plane_1").val() +
      " " +
      $("#plane_2").val() +
      " " +
      $("#plane_3").val() +
      "} colorscheme sets translucent 0.5 green";
    Jmol.script(jmolApplet0, scriptCommand);
  });
});
