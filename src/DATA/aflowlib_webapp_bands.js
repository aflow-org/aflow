// aflowlib_webapp_bands.js
//
// author: Harvey Shi
// edited: Geena Gomez (February 2018)
// edited: Pauline Colinet (May 2018)
//
// JS file for the interactive bands plot and
//    density of states plot on the entry page.

/*
function getMax(arrayOfArray) {
    for (var i = 1; i< arrayOfArray.length; i++) {
      arr = arrayOfArray[i];
      intMax = 0;
      var max = arr.reduce(function(a, b) {
            return Math.max(a, b);
      });
      finalMax = Math.max(max, intMax);
      intMax = max;
      return finalMax
    }
};
*/

$(document).ready(function() {
  var json = d3_bands_data,
    lineLabels = json.kpoint_labels_html,
    lineLocs = json.kpoint_positions,
    title = json.title,
    yExtrema = [-5, 5]; //json.Emin, json.Emax];

  // Reformat line labels
  for (var i = 0; i < lineLabels.length; ++i) {
    var label = lineLabels[i];
    label = label.replace("G", "G");
    label = label.replace("Sigm", "&Sigma;");	//CO20190301
    label = label.replace("_1", "<tspan class='subScript' dy='3'>1</tspan>");
    lineLabels[i] = label;
  }

  var spinP = json.pDOS_data.spin_polarized;
  if (spinP == true) {
    // Spin Polarized Case
    var bandData = [json.bands_data_majority, json.bands_data_minority],
      xValues = function(d) {
        return d[0];
      },
      bandX = [d3.min(bandData[0], xValues), d3.max(bandData[0], xValues)],
      numBands = json.n_bands,
      rows = numBands * 2,
      cols = bandData[0].length,
      maxBandPts = cols,
      bandArray = new Array();
    for (var i = 0; i < rows; ++i) {
      bandArray.push([]);
    }
    for (var i = 0; i < rows; ++i) {
      for (var j = 0; j < cols; ++j) {
        if (i < numBands) {
          bandArray[i].push({
            x: bandData[0][j][0],
            y: bandData[0][j][i + 1]
          });
        } else {
          bandArray[i].push({
            x: bandData[1][j][0],
            y: bandData[1][j][i - numBands + 1]
          });
        }
      }
    }
    var bandXArray = new Array();
    for (var j = 0; j < cols; ++j) {
      bandXArray.push(bandData[0][j][0]);
    }
    var dosData = [
        json.pDOS_data.energy,
        json.pDOS_data.sum_s_majority,
        json.pDOS_data.sum_p_majority,
        json.pDOS_data.sum_d_majority,
        json.pDOS_data.sum_s_minority,
        json.pDOS_data.sum_p_minority,
        json.pDOS_data.sum_d_minority
      ],
      dosNames = [
        "s majority",
        "p majority",
        "d majority",
        "s minority",
        "p minority",
        "d minority"
      ],
      numDos = dosNames.length;
    var dosArray = new Array(),
      rows = dosData.length - 1,
      cols = dosData[0].length;

    //PC20180809 - Calculates x range for DOS plot
    var minDosX = 1000;
    var maxDosX = -1000;
    for (var i = 0; i < dosData[0].length; i++) {
      if ((dosData[0][i] > -5) && (dosData[0][i] < 5)) {
      var arrTemp = [] 
      for (var j = 1; j < dosData.length; j++) {
          arrTemp.push(dosData[j][i]);
       }
      var minArrTemp = Math.min.apply(null, arrTemp);
      minDosX = Math.min(minDosX, minArrTemp);
      var maxArrTemp = Math.max.apply(null, arrTemp);
      maxDosX = Math.max(maxDosX, maxArrTemp);
      }
    }
    maxDosX = maxDosX +1;
    minDosX = minDosX -1;
    finalMaxDosX = Math.max(Math.abs(minDosX),maxDosX);
    for (var i = 0; i < rows; ++i) {
      dosArray.push([]);
    }
    for (var i = 0; i < rows; ++i) {
      for (var j = 0; j < cols; ++j) {
        dosArray[i].push({
          x: dosData[i + 1][j],
          y: dosData[0][j]
        });
      }
    }

    $("#bandLegend").fadeIn(200);
    $("#dosText").css("height", "60px");
    var dosYArray = dosData[0],
      dosColor = [
        "#66ff29",
        "#0662F9",
        "#FA7943",
        "#66ff29",
        "#0662F9",
        "#FA7943"
      ],
      legendNames = ["s", "p", "d"],
      dosData = [
        json.pDOS_data.sum_s_majority,
        json.pDOS_data.sum_p_majority,
        json.pDOS_data.sum_d_majority,
        json.pDOS_data.sum_s_minority,
        json.pDOS_data.sum_p_minority,
        json.pDOS_data.sum_d_minority
      ];
    if (numDos == 4) {
      dosData.push(json.pDOS_data.sum_f);
    }

    var min = d3.min(
      dosData.map(function(array) {
        return d3.min(array);
      })
    );
    var max = d3.max(
      dosData.map(function(array) {
        return d3.max(array);
      })
    );
    //PC20180809
  // dosX = [-10, 10];
    //dosX = [min-1, max+1];
     dosX = [-finalMaxDosX, finalMaxDosX];

    // Transpose matrix
    dosData = dosData[0].map(function(col, i) {
      return dosData.map(function(row) {
        return row[i];
      });
    });
    if (numDos == 4) {
      dosData.push(json.pDOS_data.sum_f);
    }
  } else {
    // Not Spin Polarized Case
    var bandData = json.bands_data,
      numBands = json.n_bands,
      xValues = function(d) {
        return d[0];
      },
      bandX = [d3.min(bandData, xValues), d3.max(bandData, xValues)],
      bandArray = new Array(),
      rows = numBands,
      cols = bandData.length,
      maxBandPts = cols;
    for (var i = 0; i < rows; ++i) {
      bandArray.push([]);
    }
    for (var i = 0; i < rows; ++i) {
      for (var j = 0; j < cols; ++j) {
        bandArray[i].push({
          x: bandData[j][0],
          y: bandData[j][i + 1]
        });
      }
    }
    var bandXArray = new Array();
    for (var j = 0; j < cols; ++j) {
      bandXArray.push(bandData[j][0]);
    }
    var dosData = [
        json.pDOS_data.energy,
        json.pDOS_data.sum_s,
        json.pDOS_data.sum_p,
        json.pDOS_data.sum_d
      ],
      dosNames = json.pDOS_data.orbitals,
      numDos = dosNames.length;
    if (numDos == 4) {
      dosData.push(json.pDOS_data.sum_f);
    }
    var dosArray = new Array(),
      rows = dosData.length - 1,
      cols = dosData[0].length;
    //PC20180809 - Calculates x range for DOS plot
    var maxDosX = -1000;
    for (var i = 0; i < dosData[0].length; i++) {
      if ((dosData[0][i] > -5) && (dosData[0][i] < 5)) {
      var arrTemp = [] 
      for (var j = 1; j < dosData.length; j++) {
          arrTemp.push(dosData[j][i]);
       }
      var maxArrTemp = Math.max.apply(null, arrTemp);
      maxDosX = Math.max(maxDosX, maxArrTemp);
      }
    }
    maxDosX = maxDosX +1;

    for (var i = 0; i < rows; ++i) {
      dosArray.push([]);
    }
    for (var i = 0; i < rows; ++i) {
      for (var j = 0; j < cols; ++j) {
        dosArray[i].push({
          x: dosData[i + 1][j],
          y: dosData[0][j]
        });
      }
    }

    var dosYArray = dosData[0],
      dosColor = ["#66ff29", "#0662F9", "#FA7943", "#19e3e3"],
      dosData = [
        json.pDOS_data.sum_s,
        json.pDOS_data.sum_p,
        json.pDOS_data.sum_d
      ];
    if (numDos == 4) {
      dosData.push(json.pDOS_data.sum_f);
      var dosData3 = dosData[3],
        dosDataFiltered3 = new Array();
    }

    var dosData0 = dosData[0],
      dosData1 = dosData[1],
      dosData2 = dosData[2]; // separating dosData data sets

    var rangeIndexArray = new Array();
    for (i = 0; i < cols; ++i) {
      if (dosYArray[i] <= 10 && dosYArray[i] >= -10) {
        rangeIndexArray.push(i);
      }
    }
    var rangeStart = rangeIndexArray[0],
      rangeEnd = rangeIndexArray.length + rangeStart; // beginning and end of the range of relevant energy value indexes

    var dosDataFiltered0 = new Array(),
      dosDataFiltered1 = new Array(),
      dosDataFiltered2 = new Array();
    for (i = rangeStart; i < rangeEnd; ++i) {
      dosDataFiltered0.push(dosData0[i]);
    }
    for (i = rangeStart; i < rangeEnd; ++i) {
      dosDataFiltered1.push(dosData1[i]);
    }
    for (i = rangeStart; i < rangeEnd; ++i) {
      dosDataFiltered2.push(dosData2[i]);
    }
    if (numDos == 4) {
      for (i = rangeStart; i < rangeEnd; ++i) {
        dosDataFiltered3.push(dosData3[i]);
      }
      var dosDataFiltered = [
        dosDataFiltered0,
        dosDataFiltered1,
        dosDataFiltered2,
        dosDataFiltered3
      ]; // new dosData array with values within relevant indexes
    } else {
      var dosDataFiltered = [
        dosDataFiltered0,
        dosDataFiltered1,
        dosDataFiltered2
      ];
    }

    //for (var i = 0; i< dosData.length; i++) {  
    //}

    var max = d3.max(
      dosDataFiltered.map(function(array) {
        return d3.max(array);
      })
    );
    //PC20180809 
    //dosX = [0, max + 2];
     dosX = [0, maxDosX];

    // [OBSOLETE]    var max = d3.max(dosData.map(function (array) {
    // [OBSOLETE]        return d3.max(array);
    // [OBSOLETE]   }));

    // Transpose matrix
    dosData = dosData[0].map(function(col, i) {
      return dosData.map(function(row) {
        return row[i];
      });
    });
  }
  var bands = d3.select("#bands_wrapper"),
    WIDTH = 900,
    HEIGHT = 550,
    MARGINS = {
      top: 30,
      right: 10,
      bottom: 30,
      left: 60
    },
    d3xR = d3.scale
      .linear()
      .range([MARGINS.left, WIDTH - MARGINS.right])
      .domain([bandX[0], bandX[1]]),
    d3yR = d3.scale
      .linear()
      .range([HEIGHT - MARGINS.top, MARGINS.bottom])
      .domain([yExtrema[0], yExtrema[1]]),
    xAxis = d3.svg
      .axis()
      .scale(d3xR)
      .tickSubdivide(true)
      .tickValues(false),
    yAxis = d3.svg
      .axis()
      .scale(d3yR)
      .tickSize(5)
      .orient("left")
      .tickSubdivide(true);
  var dosOffset = WIDTH;
  bands
    .append("clipPath")
    .attr("id", "clip")
    .append("rect")
    .attr("width", WIDTH - MARGINS.left - MARGINS.right)
    .attr("height", HEIGHT - MARGINS.top - MARGINS.bottom)
    .attr("x", MARGINS.left)
    .attr("y", MARGINS.top);

  bands
    .append("svg:g")
    .attr("class", "x axis")
    .attr("transform", "translate(0," + (HEIGHT - MARGINS.bottom) + ")")
    .call(xAxis);
  // Add k-point labels
  for (var i = 0; i < lineLabels.length; ++i) {
    var xPos;
    if (i == lineLabels.length - 1) {
      xPos = d3xR(d3.max(bandData[0], xValues));
    } else {
      xPos = d3xR(lineLocs[i]);
    }
    var group = bands
      .append("g")
      .attr("class", "lineLabels")
      .attr("data-xPos", xPos)
      .attr("transform", function() {
        return "translate(" + xPos + ",0)";
      });

    group
      .append("text")
      .html(lineLabels[i])
      .attr("x", 0)
      .attr("y", HEIGHT - MARGINS.bottom + 20)
      .attr("text-anchor", "middle");

    group
      .append("svg:line")
      .attr("class", "line")
      .attr("stroke-width", 1)
      .attr("x1", 0)
      .attr("y1", MARGINS.top)
      .attr("x2", 0)
      .attr("y2", HEIGHT - MARGINS.top);
  }

  // Pseudo-clipping mask for k-point lines and labels
  bands
    .append("rect")
    .attr("fill", "white")
    .attr("width", MARGINS.left)
    .attr("height", HEIGHT - MARGINS.bottom)
    .attr("x", 0)
    .attr("y", 0);
  bands
    .append("rect")
    .attr("fill", "white")
    .attr("width", MARGINS.left - 5)
    .attr("height", MARGINS.bottom)
    .attr("x", 0)
    .attr("y", HEIGHT - MARGINS.bottom);
  bands
    .append("rect")
    .attr("fill", "white")
    .attr("width", MARGINS.right)
    .attr("height", HEIGHT - MARGINS.bottom)
    .attr("x", WIDTH - MARGINS.right)
    .attr("y", 0);
  bands
    .append("rect")
    .attr("fill", "white")
    .attr("width", MARGINS.right - 5)
    .attr("height", MARGINS.bottom)
    .attr("x", WIDTH - MARGINS.right + 5)
    .attr("y", HEIGHT - MARGINS.bottom);

  bands
    .append("svg:g")
    .attr("class", "y axis")
    .attr("transform", "translate(" + MARGINS.left + ",0)")
    .call(yAxis);

  bands
    .append("text")
    .attr("class", "titles")
    .attr("x", WIDTH / 2)
    .attr("y", 15)
    .attr("text-anchor", "middle")
    .text(title);

  bands
    .append("text")
    .attr("class", "titles")
    .attr("x", -(HEIGHT / 2 + MARGINS.top + 20))
    .attr("y", 0)
    .attr("dy", ".75em")
    .attr("transform", "rotate(-90)")
    .text("Energy (eV)");

  bands
    .append("svg:svg")
    .attr("id", "bandMove")
    .attr("fill", "transparent")
    .attr("width", WIDTH - MARGINS.left - MARGINS.right)
    .attr("height", HEIGHT - MARGINS.top - MARGINS.bottom)
    .attr("x", MARGINS.left)
    .attr("y", MARGINS.top);

  var line = d3.svg
    .line()
    .interpolate("linear")
    .x(function(d) {
      return d3xR(d.x);
    })
    .y(function(d) {
      return d3yR(d.y);
    });

  if (spinP == true) {
    var minColor = "#ff0d00";
    $(".legendLine.min").css("border-color", "#ff0d00");
    var band_lines = bands
      .selectAll(".band_lines")
      .data(bandArray)
      .enter()
      .append("path")
      .attr("class", "band_lines")
      .attr("stroke", function(d, i) {
        if (i < numBands) {
          return "black";
        } else {
          return minColor;
        }
      })
      .attr("id", function(d, i) {
        if (i < numBands) {
          return "majSpin";
        } else {
          return "minSpin";
        }
      })
      .attr("stroke-width", 1)
      .attr("fill", "none")
      .attr("clip-path", "url(#clip)")
      .attr("d", line);

    var bandTotal = [],
      max = bandData[0].length;
    for (var i = 0; i < max; ++i) {
      bandTotal.push(bandData[0][i].concat(bandData[1][i]));
    }
  } else {
    var band_lines = bands.selectAll(".band_lines").data(bandArray);
    band_lines
      .enter()
      .append("path")
      .attr("class", "band_lines")
      .attr("stroke", "rgba(0,0,0,.7)")
      .attr("stroke-width", 1)
      .attr("fill", "none")
      .attr("clip-path", "url(#clip)")
      .attr("d", line);
    var bandTotal = bandData;
  }

  bands.on("mousemove", mousemoveSpin).on("mouseover", function() {
    $("#dosTip")
      .attr("opacity", 0)
      .attr("cx", 0 + "px")
      .attr("cy", 0 + "px");
    $("#dosText").css("opacity", 0);
  });
  var bisect = d3.bisector(function(d) {
    return d;
  }).right;

  //modif Pauline
  //Creation of two simple arrays (one for majority, one for minority) with the y coordinates of each band.
  var bandMajority = [];
  var bandMinority = [];
  for (i = 0; i < bandTotal.length; i++) {
    var bandM_i = [];
    var bandm_i = [];
    for (j = 1; j <= 0.5 * (bandTotal[0].length - 2); j++) {
      bandM_i.push(bandTotal[i][j]);
      bandm_i.push(bandTotal[i][j + 0.5 * bandTotal[0].length]);
    }
    bandMajority.push(bandM_i);
    bandMinority.push(bandm_i);
  }

  //modif Pauline: new mousemove function
  //mousemove changed so that it takes into account the selected value of the Spin Selection button
  //cf var spinSelectorValue
  //
  //modif Pauline: mousemoveSpin works also for non-spin polarized cases
  //i.e. SpinP != true

  function mousemoveSpin() {
    var mouseX = d3.mouse(this)[0],
      mouseY = d3.mouse(this)[1];
    var xPos = d3xR.invert(mouseX),
      xIndex = bisect(bandXArray, xPos);
    if (xIndex == maxBandPts) {
      --xIndex;
    }

    var bandTot = [];
    if (document.getElementById("spinBandsOptions") == null) {
      var spinSelectorValue = "none";
    } else {
      var spinSelectorValue = document.getElementById("spinBandsOptions").value;
    }

    if (xPos > 0 && xPos < 0.01 + bandXArray[bandXArray.length - 1]) {
      $("#bandTip").attr("visibility", "visible");
      if (spinSelectorValue == "majority") {
        $("#bandTip rect").attr("stroke", "black");
        bandTot = bandMajority;
        bandTotBis = bandMajority;
      } else if (spinSelectorValue == "minority") {
        bandTot = bandMinority;
        bandTotBis = bandMinority;
        $("#bandTip rect").attr("stroke", minColor);
      } else if (spinSelectorValue == "bothS") {
        bandTot = bandMajority;
        bandTotBis = bandMinority;
      } else if (spinSelectorValue == "none") {
        bandTot = bandTotal;
        bandTotBis = bandTotal;
        $("#bandTip rect").attr("stroke", "black");
      }

      var yPos = d3yR.invert(mouseY),
        i = bisect(bandTot[xIndex], yPos),
        d0 = bandTot[xIndex][i - 1],
        d1 = bandTot[xIndex][i],
        index = yPos - d0 > d1 - yPos ? i : i - 1,
        y1 = bandTot[Math.max(xIndex - 1, 0)][index],
        y2 = bandTot[xIndex][index],
        x1 = bandXArray[Math.max(xIndex - 1, 0)],
        x2 = bandXArray[xIndex];
      if (x2 - x1 == 0) {
        var yLocOne = 0;
      } else {
        var yLocOne = y1 + (xPos - x1) * (y2 - y1) / (x2 - x1); // Linear interpolation
      }
      var yPosBis = d3yR.invert(mouseY),
        iBis = bisect(bandTotBis[xIndex], yPosBis),
        d0Bis = bandTotBis[xIndex][iBis - 1],
        d1Bis = bandTotBis[xIndex][iBis],
        indexBis = yPosBis - d0Bis > d1Bis - yPosBis ? iBis : iBis - 1,
        y1Bis = bandTotBis[Math.max(xIndex - 1, 0)][indexBis],
        y2Bis = bandTotBis[xIndex][indexBis],
        x1Bis = bandXArray[Math.max(xIndex - 1, 0)],
        x2Bis = bandXArray[xIndex];
      if (x2Bis - x1Bis == 0) {
        var yLocBis = 0;
      } else {
        var yLocBis =
          y1Bis + (xPos - x1Bis) * (y2Bis - y1Bis) / (x2Bis - x1Bis); // Linear interpolation
      }

      yLoc =
        Math.abs(yPosBis - yLocBis) > Math.abs(yPos - yLocOne)
          ? yLocOne
          : yLocBis;

      if (yLocOne !== yLocBis) {
        if (yLoc === yLocOne) {
          $("#bandTip rect").attr("stroke", "black");
        } else {
          $("#bandTip rect").attr("stroke", minColor);
        }
      }
      if (d3xR(xPos) > 735) {
        $("#bandTip")
          .attr("x", d3xR(xPos) - 108 + "px")
          .attr("y", d3yR(yLoc) - 16 + "px");
        $("#bandTip circle").attr("cx", "104");
        $("#bandTip polygon").attr("points", "96, 16, 89, 10, 89, 22 ");
        $("#bandTip rect").attr("x", 2);
        $("#bandTip text").attr("x", 45);
      } else {
        $("#bandTip")
          .attr("x", d3xR(xPos) - 6 + "px")
          .attr("y", d3yR(yLoc) - 16 + "px");
        $("#bandTip circle").attr("cx", "6");
        $("#bandTip polygon").attr("points", "15, 16, 22, 10, 22, 22 ");
        $("#bandTip rect").attr("x", 21);
        $("#bandTip text").attr("x", 65);
      }
      $("#bandText").text(yLoc.toFixed(4) + " eV");
      checkBounds(mouseX, mouseY, d3yR(yLoc) + 15);
    } else {
      $("#bandTip").attr("visibility", "hidden");
    }
  }

 function mousemove() {
    var mouseX = d3.mouse(this)[0],
      mouseY = d3.mouse(this)[1];
    var xPos = d3xR.invert(mouseX),
      xIndex = bisect(bandXArray, xPos);
    if (xIndex == maxBandPts) {
      --xIndex;
    }

    var yPos = d3yR.invert(mouseY),
      i = bisect(bandTotal[xIndex], yPos),
      d0 = bandTotal[xIndex][i - 1],
      d1 = bandTotal[xIndex][i],
      index = yPos - d0 > d1 - yPos ? i : i - 1,
      y1 = bandTotal[Math.max(xIndex - 1, 0)][index],
      y2 = bandTotal[xIndex][index],
      x1 = bandXArray[Math.max(xIndex - 1, 0)],
      x2 = bandXArray[xIndex];
    if (x2 - x1 == 0) {
      var yLoc = 0;
    } else {
      var yLoc = y1 + (xPos - x1) * (y2 - y1) / (x2 - x1); // Linear interpolation
    }
    if (index > numBands) {
      $("#bandTip rect").attr("stroke", minColor);
    } else {
      $("#bandTip rect").attr("stroke", "black");
    }
    if (d3xR(xPos) > 735) {
      $("#bandTip")
        .attr("x", d3xR(xPos) - 108 + "px")
        .attr("y", d3yR(yLoc) - 16 + "px");
      $("#bandTip circle").attr("cx", "104");
      $("#bandTip polygon").attr("points", "96, 16, 89, 10, 89, 22 ");
      $("#bandTip rect").attr("x", 2);
      $("#bandTip text").attr("x", 45);
    } else {
      $("#bandTip")
        .attr("x", d3xR(xPos) - 6 + "px")
        .attr("y", d3yR(yLoc) - 16 + "px");
      $("#bandTip circle").attr("cx", "6");
      $("#bandTip polygon").attr("points", "15, 16, 22, 10, 22, 22 ");
      $("#bandTip rect").attr("x", 21);
      $("#bandTip text").attr("x", 65);
    }
    $("#bandText").text(yLoc.toFixed(4) + " eV");
    checkBounds(mouseX, mouseY, d3yR(yLoc) + 15);
  }

  var bandXB = [2];
  bandXB[0] = [MARGINS.left, WIDTH - MARGINS.right];
  bandXB[1] = [MARGINS.top, HEIGHT - MARGINS.bottom];
  var heightM = [MARGINS.top + 20, HEIGHT - MARGINS.bottom + 20];
  //modif Pauline, one more conditon on x
  function checkBounds(x, y, toolYpos) {
    if (
      bandXB[0][0] <= x &&
      bandXB[1][0] <= y &&
      x <= bandXB[0][1] &&
      y <= bandXB[1][1] &&
      heightM[0] <= toolYpos &&
      toolYpos <= heightM[1]
    ) {
      $("#bandTip").attr("opacity", 1);
    } else {
      $("#bandTip")
        .attr("opacity", 0)
        .attr("x", "0px")
        .attr("y", "0px");
    }
  }

  bands
    .append("line")
    .attr("class", "centerLine")
    .attr("stroke-width", 1)
    .attr("x1", MARGINS.left)
    .attr("y1", 0)
    .attr("x2", WIDTH - MARGINS.right)
    .attr("y2", 0)
    .attr("transform", function() {
      return "translate(0," + d3yR(0) + ")";
    });

  bands
    .append("svg:line")
    .attr("class", "border")
    .attr("stroke", "black")
    .attr("stroke-width", 1)
    .attr("fill", "none")
    .attr("x1", MARGINS.left)
    .attr("y1", MARGINS.top)
    .attr("x2", WIDTH - MARGINS.right)
    .attr("y2", MARGINS.top);
  bands
    .append("svg:line")
    .attr("class", "border")
    .attr("stroke", "black")
    .attr("stroke-width", 1)
    .attr("fill", "none")
    .attr("x1", WIDTH - MARGINS.right)
    .attr("y1", MARGINS.top)
    .attr("x2", WIDTH - MARGINS.right)
    .attr("y2", HEIGHT - MARGINS.top);
  var bandTip = bands.append("svg").attr("id", "bandTip");
  bandTip
    .append("circle")
    .attr("cx", "6")
    .attr("cy", "16")
    .attr("fill", "rgba(0, 0, 0, .6)")
    .attr("r", "6");
  bandTip
    .append("polygon")
    .attr("points", "15, 16, 22, 10, 22, 22 ")
    .attr("fill", "#333");
  bandTip
    .append("rect")
    .attr("width", 88)
    .attr("height", 28)
    .attr("x", 21)
    .attr("y", 2)
    .attr("rx", 9)
    .attr("ry", 9)
    .style("fill", "rgba(0,0,0,.8)")
    .attr("stroke-width", 2)
    .attr("stroke", "steelblue");
  bandTip
    .append("text")
    .attr("x", 65)
    .attr("y", 21)
    .attr("fill", "white")
    .attr("font-size", "14px")
    .attr("text-anchor", "middle")
    .attr("id", "bandText");
  $("#bandTip").attr("opacity", 0);

  // Dos plot code
  var dos = d3.select("#dos_wrapper"),
    WIDTH = 300,
    HEIGHT = 550,
    MARGINS = {
      top: 30,
      right: 20,
      bottom: 30,
      left: 10
    },
    d3xR2 = d3.scale
      .linear()
      .range([MARGINS.left, WIDTH - MARGINS.right])
      .domain([dosX[0], dosX[1]]),
    xAxis_dos = d3.svg
      .axis()
      .scale(d3xR2)
      .tickSubdivide(true)
      .ticks(5) 
      .tickSize(5),
    yAxis_dos = d3.svg
      .axis()
      .scale(d3yR)
      .orient("left")
      .tickPadding(10)
      .tickSubdivide(true);

  dos.on("mousemove", mousemove).on("mouseover", function() {
    $("#bandTip")
      .attr("opacity", 0)
      .attr("x", 0 + "px")
      .attr("y", 0 + "px");
  });

  var bisect = d3.bisector(function(d) {
    return d;
  }).right;
  var dosXB = [2];

  dosXB[0] = [MARGINS.left, WIDTH - MARGINS.right];
  dosXB[1] = [MARGINS.top, HEIGHT - MARGINS.bottom];
  var widthM = [MARGINS.left + 20, WIDTH - MARGINS.right + 20];

  function checkBounds2(x, y, toolXpos) {
    if (
      dosXB[0][0] <= x &&
      x <= dosXB[0][1] &&
      dosXB[1][0] <= y &&
      y <= dosXB[1][1] &&
      widthM[0] <= toolXpos &&
      toolXpos <= widthM[1]
    ) {
      $("#dosTip").attr("opacity", 1);
      $("#dosText").css("opacity", 1);
    } else {
      $("#dosTip")
        .attr("opacity", 0)
        .attr("cx", 0 + "px")
        .attr("cy", 0 + "px");
      $("#dosText").css("opacity", 0);
    }
  }
  dos.on("mousemove", function() {
    var mouseX = d3.mouse(this)[0],
      mouseY = d3.mouse(this)[1];

    var xPos = d3xR2.invert(mouseX),
      yPos = d3yR.invert(mouseY),
      yIndex = Math.min(dosYArray.length - 1, bisect(dosYArray, yPos)),
      dosUnsorted = dosData[yIndex],
      dosSorted = dosUnsorted.slice(0).sort(function(a, b) {
        return a - b;
      }),
      i = Math.max(0, bisect(dosSorted, xPos) - 1),
      d0 = dosSorted[i],
      d1 = dosSorted[Math.min(i + 1, numDos - 1)],
      xIndex = xPos - d0 > d1 - xPos ? Math.min(i + 1, numDos - 1) : i,
      realIndex = dosUnsorted.indexOf(dosSorted[xIndex]),
      name = dosNames[realIndex],
      x1 = dosData[Math.max(yIndex - 1, 0)][realIndex],
      x2 = dosData[yIndex][realIndex],
      y1 = dosYArray[Math.max(yIndex - 1, 0)],
      y2 = dosYArray[yIndex],
      yLoc = y1 + (xPos - x1) * (y2 - y1) / (x2 - x1); // Linear interpolation

    $("#dosTip")
      .attr("cx", d3xR2(x2) + "px")
      .attr("cy", d3yR(yPos) + "px");
    $("#dosText")
      .html(name + ": " + x2.toFixed(4) + " States/eV")
      .css("border-color", dosColor[realIndex]);
    checkBounds2(mouseX, mouseY, d3xR2(x2) + 21);
  });

  dos
    .append("clipPath")
    .attr("id", "clip2")
    .append("rect")
    .attr("width", WIDTH - MARGINS.left - MARGINS.right)
    .attr("height", HEIGHT - MARGINS.top - MARGINS.bottom)
    .attr("x", MARGINS.left)
    .attr("y", MARGINS.top);

  dos
    .append("svg:g")
    .attr("class", "xAxisDos")
    .attr("transform", "translate(0," + (HEIGHT - MARGINS.bottom) + ")")
    .call(xAxis_dos);

  dos
    .append("text")
    .attr("class", "titles")
    .attr("x", WIDTH / 2)
    .attr("y", 15)
    .attr("text-anchor", "middle")
    .text("Density of States (States/eV)");

  if (spinP == true) {
    // Add Legend Entries
    for (var i = 0; i < legendNames.length; ++i) {
      $("#dosLegend").append(
        '<text  class="legendText">' +
          legendNames[i] +
          '<line class="legendLine" style="border-color:' +
          dosColor[i] +
          ';"/></text>'
      );
    }
  } else {
    // Add Legend Entries
    for (var i = 0; i < dosNames.length; ++i) {
      $("#dosLegend").append(
        '<text  class="legendText">' +
          dosNames[i] +
          '<line class="legendLine" style="border-color:' +
          dosColor[i] +
          ';"/></text>'
      );
    }
  }
  $("#dosLegend").fadeIn(350);

  var line2 = d3.svg
    .line()
    .interpolate("linear")
    .x(function(d) {
      return d3xR2(d.x);
    })
    .y(function(d) {
      return d3yR(d.y);
    });
  dos
    .selectAll(".dos_data")
    .data(dosArray)
    .enter()
    .append("path")
    .attr("class", "dos_data")
    .attr("stroke", function(d, i) {
      return dosColor[i];
    })
    .attr("stroke-width", 1)
    .attr("fill", "none")
    .attr("clip-path", "url(#clip2)")
    .attr("d", line2);

  dos
    .append("line")
    .attr("class", "centerLine")
    .attr("stroke-width", 1)
    .attr("x1", MARGINS.left)
    .attr("y1", 0)
    .attr("x2", WIDTH - MARGINS.right)
    .attr("y2", 0)
    .attr("transform", function() {
      return "translate(0," + d3yR(0) + ")";
    });

  dos
    .append("svg:line")
    .attr("class", "border")
    .attr("stroke", "black")
    .attr("stroke-width", 1)
    .attr("fill", "none")
    .attr("x1", MARGINS.left)
    .attr("y1", MARGINS.top)
    .attr("x2", WIDTH - MARGINS.right)
    .attr("y2", MARGINS.top);
  dos
    .append("svg:line")
    .attr("class", "border")
    .attr("stroke", "black")
    .attr("stroke-width", 1)
    .attr("fill", "none")
    .attr("x1", WIDTH - MARGINS.right)
    .attr("y1", MARGINS.top)
    .attr("x2", WIDTH - MARGINS.right)
    .attr("y2", HEIGHT - MARGINS.top);
  dos
    .append("svg:line")
    .attr("class", "border")
    .attr("stroke", "black")
    .attr("stroke-width", 1)
    .attr("fill", "none")
    .attr("x1", MARGINS.left)
    .attr("y1", MARGINS.top)
    .attr("x2", MARGINS.left)
    .attr("y2", HEIGHT - MARGINS.top);
  dos
    .append("circle")
    .attr("id", "dosTip")
    .attr("fill", "rgba(0, 0, 0, .6)")
    .attr("r", "6");
  $("#dosTip").attr("opacity", "0");

  function zoomFunc() {
    //console.log(d3xR(0), d3xR2(0));
    bands.selectAll(".band_lines").attr("d", line);
    bands.select(".x.axis").call(xAxis);
    bands.select(".y.axis").call(yAxis);
    bands.select(".centerLine").attr("transform", function() {
      return "translate(0," + d3yR(0) + ")";
    });
    bands.selectAll(".lineLabels").attr("transform", function(d, i) {
      return "translate(" + d3xR(lineLocs[i]) + ",0)";
    });
    dos.selectAll("path.dos_data").attr("d", line2);
    dos.select(".xAxisDos").call(xAxis_dos);
    dos.select(".centerLine").attr("transform", function() {
      return "translate(0," + d3yR(0) + ")";
    });
  }
  var width = WIDTH - MARGINS.left - MARGINS.right,
    height = HEIGHT - MARGINS.top - MARGINS.bottom,
    zoom1b = d3.behavior
      .zoom()
      .x(d3xR)
      .y(d3yR)
      .scaleExtent([1, 100])
      .on("zoom", zoomFunc),
    zoom1x = d3.behavior
      .zoom()
      .x(d3xR)
      .scaleExtent([1, 100])
      .on("zoom", zoomFunc),
    zoom1y = d3.behavior
      .zoom()
      .y(d3yR)
      .scaleExtent([1, 100])
      .on("zoom", zoomFunc),
    zoom2b = d3.behavior
      .zoom()
      .x(d3xR2)
      .y(d3yR)
      .scaleExtent([1, 100])
      .on("zoom", zoomFunc),
    zoom2x = d3.behavior
      .zoom()
      .x(d3xR2)
      .scaleExtent([1, 100])
      .on("zoom", zoomFunc),
    zoom2y = d3.behavior
      .zoom()
      .y(d3yR)
      .scaleExtent([1, 100])
      .on("zoom", zoomFunc),
    zoom1 = [zoom1b, zoom1x, zoom1y],
    zoom2 = [zoom2b, zoom2x, zoom2y];

  // Clicking the Reset Zoom button
  $(".reset").click(function() {
    bands.call(zoom1[0].translate([0, 0]).scale(1).event);
    dos.call(zoom2[0].translate([0, 0]).scale(1).event);
  });

  bands.call(zoom1[0]);
  dos.call(zoom2[0]);

  $("#zoomOptions select").change(function() {
    var selectV = $(this).val();
    //console.log("test1");
    switch (selectV) {
      case "both":
        //console.log("test2");
        bands.call(zoom1[0]);
        dos.call(zoom2[0]);
        break;
      case "xOnly":
        bands.call(zoom1[1]);
        dos.call(zoom2[1]);
        break;
      case "yOnly":
        bands.call(zoom1[2]);
        dos.call(zoom2[2]);
        break;
    }
  });
  if (spinP == true) {
    $("#spinBandsOptions select").change(function() {
      var selectVal = $(this).val();
      switch (selectVal) {
        case "majority":
          break;
        case "minority":
          break;
      }
    });
  }
  $("#show").fadeIn(1000);
});
