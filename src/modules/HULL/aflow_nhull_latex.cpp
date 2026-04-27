// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by Nicholas Anderson
// nicholas.anderson@duke.edu

#include "aflow_nhull_latex.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <iosfwd>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <extern/SYMBOLICCPLUSPLUS/symbolicc++.h>

#include "AUROSTD/aurostd_xplotter.h"
#include "AUROSTD/aurostd_xvector.h"

#include "../../flow/aflow_figures.h"
#include "modules/HULL/aflow_nhull.h"
#include "modules/HULL/aflow_nhull_entry.h"
#include "modules/HULL/aflow_nhull_facet.h"
#include "modules/HULL/aflow_nhull_util.h"

using std::cerr;
using std::cout;
using std::deque;
using std::endl;
using std::ifstream;
using std::ios_base;
using std::istream;
using std::istringstream;
using std::map;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::setprecision;
using std::setw;
using std::string;
using std::stringstream;
using std::to_string;
using std::tuple;
using std::vector;

//************************************* Utility Functions:

/// @brief converts full stoichiometric coords (3 species) into cartesian coords (2d)
/// @param full_stoich_enthalpy_coord hull point with first element index included
/// @return cartesion coordinates of this coordinate in ternary plot
/// @authors
/// mod{NHA,20260206,created}
vector<double> stoichToTernaryUnits(const vector<double>& full_stoich_enthalpy_coord) {
  vector<double> new_coords;
  if (full_stoich_enthalpy_coord.size() != 3) {
    cerr << "ERROR IN FUNCTION stoichToTernaryUnits(): input vector is incorrect size!" << endl;
    new_coords = {0, 0};
    return new_coords;
  }
  double a = full_stoich_enthalpy_coord[0];
  double b = full_stoich_enthalpy_coord[1];
  double c = full_stoich_enthalpy_coord[2];
  double x = (0.5) * (2 * b / (a + b + c) + c / (a + b + c));
  double y = (sqrt(3) / 2) * (c / (a + b + c));
  new_coords = {x, y};
  return new_coords;
}

/// @brief Takes a point and returns coordinates for 3D ternary plot (energy coord removed)
/// @authors
/// mod{NHA,20260206,created}
vector<double> stoichToTernaryPoint(vector<double> full_stoich_enthalpy_coord) {
  double energy = full_stoich_enthalpy_coord.back();
  full_stoich_enthalpy_coord.pop_back();
  vector<double> ternary_point = stoichToTernaryUnits(full_stoich_enthalpy_coord);
  ternary_point.push_back(energy);

  return ternary_point;
}

/// @brief Takes an entry and returns coordinates for 3D ternary plot
/// @authors
/// mod{NHA,20260206,created}
aurostd::xvector<double> entryToTernaryPoint(const nhull::Entry& entry, const vector<string>& specie_vector) {
  aurostd::xvector<double> full_stoich = entry.compoundToPoint(specie_vector, true);

  vector<double> ternary_point = stoichToTernaryPoint(xvector2vector(full_stoich));
  return aurostd::vector2xvector(ternary_point);
}

/// @brief takes vector of entries and returns only entries that are hull points
/// @authors
/// mod{NHA,20260206,created}
vector<nhull::Entry> hullEntryParse(vector<nhull::Entry> entries) {
  vector<nhull::Entry> hull_entries;
  for (auto entry : entries) {
    if (entry.is_hull_point) {
      hull_entries.push_back(entry);
    }
  }
  return hull_entries;
}

/// @brief takes prototype from Entry object and formats for latex table
/// @return formatted proto label
/// @authors
/// mod{NHA,20260206,created}
string prototypeFormat(string prototype) {
  string formatted_prototype;
  for (auto char1 : prototype) {
    if (char1 == '_') {
      formatted_prototype += "\\_";
    } else {
      formatted_prototype += char1;
    }
  }
  return formatted_prototype;
}

/// @brief takes double, rounds to nearest int and casts to int
int intRound(double value) {
  return static_cast<int>(std::round(value));
}

/// @brief returns string representation of an entry's compound with no 1's included
/// @authors
/// mod{NHA,20260206,created}
string get_compound_string_reduced(nhull::Entry entry) {
  string final_string;
  for (auto element : entry.compound_map) {
    if (element.second == 1) {
      final_string += element.first;
    } else if (element.second != 0) {
      final_string += element.first + std::to_string(element.second);
    }
  }

  return final_string;
}

/// @brief utility function does latex formatting for compound labels
/// @authors
/// mod{NHA,20260206,created}
string get_compound_string_latex(nhull::Entry entry) {
  string final_string;
  for (auto element : entry.compound_map) {
    if (element.second == 1) {
      final_string += element.first;
    } else {
      final_string += element.first + "$_{" + std::to_string(element.second) + "}$";
    }
  }

  return final_string;
}

/// @brief overload: utility function does formatting for compound labels for latex
/// @authors
/// mod{NHA,20260206,created}
string get_compound_string_latex(string compound) {
  string final_string;
  string number_store;
  for (auto char1 : compound) {
    if (!isdigit(char1) && number_store.empty()) {
      final_string += char1;
    } else if (!isdigit(char1) && !number_store.empty()) {
      if (number_store != "1") {
        final_string += "$_{" + number_store + "}$";
      }
      final_string += char1;
      number_store.clear();
    } else if (isdigit(char1) && !number_store.empty()) {
      number_store += char1;
    } else if (isdigit(char1) && number_store.empty()) {
      number_store += char1;
    }
  }
  if (!number_store.empty()) {
    final_string += "$_{" + number_store + "}$";
  }

  return final_string;
}

/// @brief Utility returns true when all entries fall above 0-energy tie line
/// @authors
/// mod{NHA,20260206,created}
bool allEntriesUnstable(const vector<nhull::Entry>& entries) {
  bool all_unstable = true;

  for (auto entry : entries) {
    if (entry.nspecies != 1 && entry.enthalpy_formation_atom <= 0) {
      all_unstable = false;
      break;
    }
  }

  return all_unstable;
}

/// @brief takes an index and returns its string representation: 2-> binary, 3-> ternary, and so on
/// @authors
/// mod{NHA,20260206,created}
string indexToString(const int i) {
  switch (i) {
    case 2: return "Binaries";
    case 3: return "Ternaries";
    case 4: return "Quaternaries";
    case 5: return "Quinaries";
    case 6: return "Senaries";
    case 7: return "Septenaries";
    case 8: return "Octonaries";
    case 9: return "Nonaries";
    case 10: return "Denaries";
    default: return std::to_string(i);
  }
}

/// @brief takes phase decomp map from an entry and converts to a string for latex
/// @note used specifically for xplotter tables
/// @authors
/// mod{NHA,20260206,created}
string phaseDecompCompact(const nhull::Entry& entry) {
  auto phase_decomp = entry.phaseDecompToMap();
  int itteration = 0;
  string red_cmpd = get_compound_string_reduced(entry);
  string string_decomp;

  string_decomp += +" $\\rightarrow$";
  for (auto compound : phase_decomp) {
    string_decomp += " ";
    if (compound.second != 0) {
      if (itteration != 0) {
        string_decomp += " + ";
      }
      string ratio = std::to_string(round(compound.second * 1000) / 1000).substr(0, 5);
      // take only first 3 decimal places of ratios
      string_decomp += ratio + "$\\cdot$" + get_compound_string_latex(compound.first);
      itteration++;
    }
  }

  return string_decomp;
}

/// @brief converts a vector containing all species of given compound into a string
string specieString(const vector<string>& specie_vector) {
  string specie_string;
  for (auto element : specie_vector) {
    specie_string += element;
  }

  return specie_string;
}

/// @brief converts a vector containing all species of given compound into a string
string vectorString(vector<double> specie_vector) {
  string specie_string;
  for (auto element : specie_vector) {
    specie_string += std::to_string(element) + " ";
  }

  return specie_string;
}

//************************* All functions associated with creating ternary plots

namespace ternary {

  /// @brief adds hull point markings:
  /// @return a stringstream containing latex formating for hull points
  /// @authors
  /// mod{NHA,20260206,created}
  stringstream addHullMarks3D(const vector<string>& specie_vector, const vector<nhull::Entry>& entries) {
    stringstream ss_out;
    string first_specie = specie_vector[0];
    vector<aurostd::xvector<double>> facet_points;
    string delimiter = "                           ";

    // get all facet points from entries
    for (auto entry : entries) {
      if (entry.is_hull_point) {
        aurostd::xvector<double> full_stoich = entryToTernaryPoint(entry, specie_vector);
        full_stoich.shift(0);
        facet_points.push_back(full_stoich);
      }
    }

    ss_out << "\\addplot3+[" << endl;
    ss_out << "only marks," << endl;
    ss_out << "mark=*," << endl;
    ss_out << "mark size=3," << endl;
    ss_out << "point meta=\\thisrow{H_f_meVatom}," << endl;
    ss_out << "nodes near coords*={}," << endl;
    ss_out << R"(visualization depends on={\thisrow{H_f_meVatom} \as \H_f_meVatom},)" << endl;
    ss_out << "] table [z=" << "H_f_meVatom" << "]{" << endl;
    ss_out << "x" << delimiter << "y" << delimiter << "H_f_meVatom" << endl;

    // add facet edges here:
    for (auto point : facet_points) {
      ss_out << point[0] << delimiter << point[1] << delimiter << point[2] << endl;
    }

    ss_out << "};" << endl;

    return ss_out;
  }

  /// @brief add labels to hull points in 3D ternary plot
  /// @return latex formated stringstream
  /// @authors
  /// mod{NHA,20260206,created}
  stringstream addHullPointLabels(const vector<nhull::Entry>& filtered_entries, const vector<string>& specie_vector) {
    stringstream ss_out;
    vector<string> prev_entries;

    ss_out << "\\begin{pgfonlayer}{foreground}" << endl;
    for (auto entry : filtered_entries) {
      if (entry.is_hull_point && entry.nspecies != 1) {
        string curr_string = get_compound_string_latex(entry);
        if (count(prev_entries.begin(), prev_entries.end(), curr_string) == 0) {
          prev_entries.push_back(curr_string);

          aurostd::xvector<double> point = entryToTernaryPoint(entry, specie_vector);
          point.shift(0);

          ss_out << "\\node [anchor=south] at (axis cs:" << point[0] << "," << point[1] << "," << point[2];
          ss_out << ") {\\color{black}\\large{" << curr_string << "}};" << endl;
        }
      }
    }

    ss_out << "\\end{pgfonlayer}{background}" << endl;

    return ss_out;
  }

  /// @brief add connecting lines (facet edges) to ternary plot
  /// @param list of facets
  /// @return stringsteam of the latex formating
  /// @authors
  /// mod{NHA,20260206,created}
  stringstream addFacetEdges3D(const vector<nhull::NhullFacet>& nhull_facets) {
    stringstream ss_out;
    vector<vector<double>> facet_points;
    string delimiter = "                           ";

    for (const auto& facet : nhull_facets) {
      ss_out << "\\addplot3 [" << endl;
      ss_out << "patch," << endl;
      ss_out << "faceted color=black," << endl;
      ss_out << "samples=15," << endl;
      ss_out << "domain=0:1,y domain=-1:1," << endl;
      ss_out << "fill opacity=0.50" << endl;
      ss_out << "] table [z=" << "H_f_meVatom" << "]{" << endl;
      ss_out << "x" << delimiter << "y" << delimiter << "H_f_meVatom" << endl;

      if (facet.vertices.size() != 3) {
        cerr << "ERROR IN FUNCTION addFacets(): " << endl;
      }
      for (auto vertex : facet.vertices) {
        vector<double> full_stoich = nhull::fullCoords(aurostd::xvector2vector(vertex), true);
        vector<double> ternary_point = stoichToTernaryPoint(full_stoich);
        ss_out << ternary_point[0] << delimiter << ternary_point[1] << delimiter << ternary_point[2] << endl;
      }
      ss_out << "};" << endl << endl;
    }

    return ss_out;
  }

  /// @brief specialized plot creation routine for ternaries projected into 3D
  /// @param specie_vector
  /// @param nhull_facets
  /// @param entries
  /// @param number_of_entries size of entries
  /// @return string of latex code to generate
  /// @authors
  /// mod{NHA,20251005,created}
  string createTernary3D(const vector<string> specie_vector, vector<nhull::NhullFacet>& nhull_facets, const vector<nhull::Entry>& entries, const int number_of_entries) {
    stringstream ss_out;
    double min_energy = nhull::Entry::minEnergy(entries);
    double element_label_height = abs(min_energy) * 0.043;
    ss_out << "\\documentclass{article}" << endl;
    ss_out << "\\usepackage{graphicx} % Required for inserting images" << endl;
    ss_out << "\\usepackage[a4paper,margin=0in,landscape]{geometry}" << endl;
    ss_out << "\\usepackage{graphicx} % Required for inserting images" << endl;
    ss_out << "\\usepackage{pgfplots}" << endl;
    ss_out << "\\usepackage{tikz}" << endl;
    ss_out << "\\usetikzlibrary{positioning}" << endl;
    ss_out << "\\usetikzlibrary{pgfplots.ternary}" << endl;
    ss_out << "\\pgfdeclarelayer{background}" << endl;
    ss_out << "\\pgfdeclarelayer{foreground}" << endl;
    ss_out << "\\pgfsetlayers{background,main,foreground}" << endl << endl;

    ss_out << "\\begin{document}" << endl;
    ss_out << "\\pgfplotsset{width=7cm,compat=1.18}" << endl;
    ss_out << "\\begin{tikzpicture}" << endl << endl;

    ss_out << "\\begin{axis}[view={5}{50}," << endl;
    ss_out << "width=.96\\textwidth," << endl;
    ss_out << "height=0.72\\textwidth," << endl;
    ss_out << "title =\\LARGE{~~~~~~~~~~~~~~~" << specieString(specie_vector) << " entries: " << number_of_entries << "},";
    ss_out << "x axis line style={opacity=0}," << endl;
    ss_out << "xtick=\\empty," << endl;
    ss_out << "y axis line style={opacity=0}," << endl;
    ss_out << "ytick=\\empty," << endl;
    ss_out << "z axis line style={opacity=0}," << endl;
    ss_out << "ztick=\\empty," << endl;
    ss_out << "colorbar," << endl;
    ss_out << "colorbar style={" << endl;
    ss_out << "title=\\Large{$H_{\\mathrm{f}}$ (meV/atom)}, " << endl;
    ss_out << "title style={at={(-.2,.89)}, align=right, rotate=90}," << endl;
    ss_out << "ytick distance = 100," << endl;
    ss_out << "yticklabel style={font = \\Large}," << endl;
    ss_out << "at={(.97,1)}," << endl;
    ss_out << "colorbar/width = 0.7cm}," << endl;
    ss_out << "colormap={mymap}{rgb(0pt)=(0,0,1)," << endl;
    ss_out << "rgb(63pt)=(1,0.644,0)}," << endl;
    ss_out << " clip = false]" << endl;

    ss_out << "\\begin{pgfonlayer}{foreground}" << endl;
    ss_out << "\\node [anchor=south] at (axis cs:0,0," << element_label_height << "){\\color{black}\\LARGE{" << specie_vector[0] << "}};" << endl;
    ss_out << "\\node [anchor=south] at (axis cs:1,0," << element_label_height << "){\\color{black}\\LARGE{" << specie_vector[1] << "}};" << endl;
    ss_out << "\\node [anchor=south] at (axis cs:0.5,0.866," << element_label_height << "){\\color{black}\\LARGE{" << specie_vector[2] << "}};" << endl;
    ss_out << "\\addplot3[mesh, draw = gray] coordinates {(0,0,0) (1,0, 0) (0.5,0.866,0) (0,0,0)};" << endl;
    ss_out << "\\end{pgfonlayer}{foreground}" << endl;

    ss_out << addHullPointLabels(entries, specie_vector).str() << endl;

    ss_out << "\\addplot3[mesh, draw = gray] coordinates {(0,0," << min_energy << ") (1,0, " << min_energy << ") (0.5,0.866," << min_energy << ") (0,0," << min_energy << ")};" << endl;
    ss_out << "\\addplot3[mesh, draw = gray] coordinates {(0,0," << min_energy << ") (0,0,0)};" << endl;
    ss_out << "\\addplot3[mesh, draw = gray] coordinates {(1,0," << min_energy << ") (1,0,0)};" << endl;
    ss_out << "\\addplot3[mesh, draw = gray] coordinates {(0.5,0.866," << min_energy << ") (0.5,0.866,0)};" << endl;

    ss_out << addHullMarks3D(specie_vector, entries).str() << endl << endl;
    ss_out << addFacetEdges3D(nhull_facets).str() << endl << endl;

    ss_out << "\\end{axis}" << endl;
    ss_out << "\\end{tikzpicture}" << endl;
    ss_out << "\\end{document}" << endl;

    return ss_out.str();
  }
} // namespace ternary

//********************** Table Functions:

///@brief utility struct for creating hull table pdfs
/// @authors
/// mod{NHA,20260206,created}
struct Compound {
  string hyper_compound; // latex formatted compound string
  bool stable = false; // indicates that one of entries  is stable compound
  vector<nhull::Entry> entries;
  int nspecies;

  Compound() = default;

  Compound(string hyper_compound, nhull::Entry first_entry, int nspecies) : hyper_compound(hyper_compound), nspecies(nspecies) {
    this->entries.push_back(first_entry);
    if (first_entry.is_hull_point) {
      this->stable = true;
    }
  }

  /// @brief takes list of entries and converts to vector<vector<string>> to help make tables
  /// @authors
  /// mod{NHA,20260206,created}
  vector<map<string, vector<string>>> entriesToTableVector() {
    vector<map<string, vector<string>>> data_store;
    //first sort entries from least to greatest dhulls:
    vector<nhull::Entry> sorted_entries = this->entries;
    std::sort(sorted_entries.begin(), sorted_entries.end(), [](const nhull::Entry& a, const nhull::Entry& b) { return a.distance_hull_enthalpy_formation_atom < b.distance_hull_enthalpy_formation_atom; });

    for (auto& entry : sorted_entries) {
      string prototype = "\\verb|" + entry.prototype + "|";
      auto enth = entry.enthalpy_formation_atom;
      auto disthull = entry.distance_hull_enthalpy_formation_atom;
      string enthalpy = std::to_string(round(enth * 1000) / 1000).substr(0, std::to_string(enth).find(".") + 4);
      // take only first 4 decimal places
      string dhull = std::to_string(round(disthull * 1000) / 1000).substr(0, std::to_string(disthull).find(".") + 4);

      vector<string> row = {entry.auid, prototype, enthalpy, dhull};
      map<string, vector<string>> row_color_map;
      if (entry.m_uses_cce) {
        row_color_map["blue"] = row;
      } else {
        row_color_map[""] = row;
      }
      data_store.push_back(row_color_map);
    }
    return data_store;
  }
};

//*************************** implementing plotting function with xplotter.h
namespace nhullxplotter {
  // add hull point markings for xplotter: outputs <x, y, z, m>

  /// @brief Plotting utility returns point coordinates for a vector of entries in 2D or 3D
  /// @param specie_vector contains all elements
  /// @param entries vector of entries
  /// @return list of vectors for each coordinate + energy
  /// authors
  /// mod{NHA,20251005,created}
  tuple<vector<double>, vector<double>, vector<double>, vector<double>> addHullMarks(const vector<string> specie_vector, const vector<nhull::Entry>& entries) {
    int dimension = specie_vector.size();
    vector<double> x;
    vector<double> y;
    vector<double> z;
    vector<double> m;

    if (dimension == 3) {
      // get all facet points from entries
      for (auto entry : entries) {
        if (entry.is_hull_point && entry.lowest_energy_alloy) {
          vector<double> curr_point = aurostd::xvector2vector(entry.compoundToPoint(specie_vector, true));
          x.push_back(curr_point[0]);
          y.push_back(curr_point[1]);
          z.push_back(curr_point[2]);
          m.push_back(curr_point[3]);
        }
      }
    } else if (dimension == 2) {
      // get all facet points from entries
      for (auto entry : entries) {
        vector<double> curr_point = aurostd::xvector2vector(entry.compoundToPoint(specie_vector, false));
        x.push_back(curr_point[0]);
        y.push_back(curr_point[1]);
      }
    } else {
      cerr << "ERROR IN FUNCTION addHullMarks(): dimension exceeded 3; no plot is being generated but this plotting function was called!" << endl;
    }

    return {x, y, z, m};
  }

  /// @brief add connecting lines (facet edges) to xplotter plot
  /// @param nhull_facets
  /// @param specie_vector
  /// @param xplt xplotter object used to create lines
  /// @author
  /// mod{NHA,20251005,created}
  void addFacetEdges(vector<nhull::NhullFacet>& nhull_facets, vector<string> specie_vector, aurostd::xplotter& xplt) {
    int dimension = specie_vector.size();
    if (dimension == 2) {
      //make hull lines thick
      std::map<std::string, std::string> properties = {
          {"style", "ultra thick"}
      };
      for (auto facet : nhull_facets) {
        // add facet points here:
        vector<double> start_point = aurostd::xvector2vector(facet.vertices[0]);
        vector<double> end_point = aurostd::xvector2vector(facet.vertices[1]);

        std::pair<double, double> start = {start_point[0], start_point[1]};
        std::pair<double, double> end = {end_point[0], end_point[1]};

        xplt.line(start, end, properties);
      }
    } else if (dimension == 3) {
      // add facet points here:
      for (auto facet : nhull_facets) {
        vector<vector<double>> plot_vertices;
        // remove energy coordinate
        for (const auto element : facet.vertices) {
          auto tmp = aurostd::xvector2vector(element);
          tmp.pop_back();
          plot_vertices.push_back(tmp);
        }

        vector<double> first_point = nhull::fullCoords(plot_vertices[0], false);
        vector<double> second_point = nhull::fullCoords(plot_vertices[1], false);
        vector<double> third_point = nhull::fullCoords(plot_vertices[2], false);

        vector<double> formatted_points1 = {first_point[0], first_point[1], second_point[0], second_point[1]};
        vector<double> formatted_points2 = {second_point[0], second_point[1], third_point[0], third_point[1]};
        vector<double> formatted_points3 = {first_point[0], first_point[1], third_point[0], third_point[1]};

        xplt.line(formatted_points1);
        xplt.line(formatted_points2);
        xplt.line(formatted_points3);
      }
    } else {
      cerr << "ERROR IN FUNCTION addFacetEdges(): dimension exceeded 3; no plot is being generated but this plotting function was called!" << endl;
    }
  }

  /// @brief add interpolated patches to xplotter plot: ONLY USED IN TERNARY PLOT
  /// @param nhull_facets
  /// @param specie_vector
  /// @param xplt xplotter object used to generate ternary scatter plot
  /// author
  /// mod{NHA,20251005,created}
  void addFacetInterpolate(vector<nhull::NhullFacet>& nhull_facets, vector<string> specie_vector, aurostd::xplotter& xplt) {
    int dimension = specie_vector.size();

    std::map<std::string, std::string> properties = {
        {              "only marks",                             "false"},
        {                   "patch",                              "true"},
        {              "patch type",                          "triangle"},
        {                  "shader",                            "interp"},
        {              "point meta",                      "\\thisrow{m}"},
        {      "nodes near coords*",                                "{}"},
        {"visualization depends on", R"({\thisrow{m} \as \H_f_meVatom})"},
        {                    "mark",                                 "*"}
    };

    // add facet points here:
    for (auto facet : nhull_facets) {
      vector<double> first_point = nhull::fullCoords(aurostd::xvector2vector(facet.vertices[0]), true);
      vector<double> second_point = nhull::fullCoords(aurostd::xvector2vector(facet.vertices[1]), true);
      vector<double> third_point = nhull::fullCoords(aurostd::xvector2vector(facet.vertices[2]), true);

      vector<double> x = {first_point[0], second_point[0], third_point[0]};
      vector<double> y = {first_point[1], second_point[1], third_point[1]};
      vector<double> z = {first_point[2], second_point[2], third_point[2]};
      vector<double> m = {first_point[3], second_point[3], third_point[3]};

      // ensure that added facets are not above zero tie lines
      if (m[0] <= 0 && m[1] <= 0 && m[2] <= 0) {
        xplt.ternary_scatter(x, y, z, m, properties);
      }
    }
  }

  /// @brief add labels and hyperrefs (links) to plot
  /// @param filtered_entries
  /// @param specie_vector
  /// @return vector of strings of latex code that generate the nodes
  /// author
  /// mod{NHA,20251005,created}
  std::vector<string> addNodes3D(const vector<nhull::Entry>& filtered_entries, const vector<string>& specie_vector) {
    std::vector<string> nodes;
    vector<string> prev_entries;
    string aflow_logo = "\\node [opacity=0.3,shift={(-2in,0.5in)},anchor=north] at (axis cs:1,0,0){\\includegraphics[scale=0.2]{aflow4_logo_skinny.png}}";

    for (auto entry : filtered_entries) {
      if (entry.is_hull_point && entry.nspecies != 1) {
        //test_i++;
        string curr_string = get_compound_string_latex(entry);
        if (count(prev_entries.begin(), prev_entries.end(), curr_string) == 0) {
          int index = -1;

          prev_entries.push_back(curr_string);

          aurostd::xvector point = entry.compoundToPoint(specie_vector, true);
          point.shift(0);

          for (int i = point.lrows; i < point.urows; i++) {
            if (point[i] == 0) {
              index = i;
            }
          }

          if (index == 0) {
            stringstream ss_out;
            ss_out << "\\node [rotate=90,anchor=east] at (axis cs:" << point[0] << "," << point[1] << "," << point[2];
            ss_out << ") {\\color{black}\\small{" << curr_string << "~~}};";
            nodes.push_back(ss_out.str());
          } else if (index == 1) {
            stringstream ss_out;
            ss_out << "\\node [rotate=30,anchor=west] at (axis cs:" << point[0] << "," << point[1] << "," << point[2];
            ss_out << ") {\\color{black}\\small{~~" << curr_string << "}};";
            nodes.push_back(ss_out.str());
          } else if (index == 2) {
            stringstream ss_out;
            ss_out << "\\node [rotate=-30,anchor=east] at (axis cs:" << point[0] << "," << point[1] << "," << point[2];
            ss_out << ") {\\color{black}\\small{" << curr_string << "~~}};";
            nodes.push_back(ss_out.str());
          } else if (index == -1) {
            stringstream ss_out;
            ss_out << "\\node at (axis cs:" << point[0] << "," << point[1] << "," << point[2];
            ss_out << ") {\\color{white}\\small{" << curr_string << "~~}};";
            nodes.push_back(ss_out.str());
          }
        }
      }
    }

    return nodes;
  }

  /// @brief add labels and hyperrefs (links) to 2D plot
  /// @param filtered_entries
  /// @param specie_vector
  /// @return vector of strings of latex code that generate the nodes
  /// author
  /// mod{NHA,20251005,created}
  std::vector<string> addNodes2D(const vector<nhull::Entry>& filtered_entries, const vector<string>& specie_vector) {
    std::vector<string> nodes;
    vector<string> prev_entries;
    string aflow_logo = "\\node [opacity=0.3,shift={(-2in,0.5in)},anchor=north] at (axis cs:1,0){\\includegraphics[scale=0.2]{aflow4_logo_skinny.png}}";

    for (auto entry : filtered_entries) {
      if (entry.is_hull_point && entry.nspecies != 1) {
        string curr_string = get_compound_string_latex(entry);
        if (count(prev_entries.begin(), prev_entries.end(), curr_string) == 0) {
          prev_entries.push_back(curr_string);
          aurostd::xvector<double> point = entry.compoundToPoint(specie_vector, true);
          point.shift(0);

          stringstream ss_out;
          ss_out << "\\node [rotate=90,anchor=east] at (axis cs:" << point[1] << "," << point[2];
          ss_out << ") {\\color{black}\\small{" << curr_string << "~~}};";
          nodes.push_back(ss_out.str());
        }
      }
    }
    // nodes.push_back(aflow_logo);

    return nodes;
  }

  /// @brief utility for creating tables. Takes in many entries and returns maps of relevent information.
  /// @param filtered_entries
  /// @return map of strings indexed according to display table
  /// @authors
  /// mod{NHA,20251005,created}
  std::map<string, vector<map<string, vector<string>>>> entriesToCmpdData(const vector<nhull::Entry>& filtered_entries) {
    stringstream ss_out;
    map<string, Compound> compounds;
    map<string, vector<map<string, vector<string>>>> compound_data;

    // sort entries by compound:
    for (const auto& entry : filtered_entries) {
      // std::map<string, int> cmpd_map = entry.compound_map;
      string red_cmpd = get_compound_string_reduced(entry);
      string hyp_cmpd = get_compound_string_latex(entry);
      int nspecies = entry.nspecies;

      if (compounds.find(red_cmpd) != compounds.end()) {
        compounds[red_cmpd].entries.push_back(entry);
        if (entry.is_hull_point) {
          compounds[red_cmpd].stable = true;
        }
      } else {
        compounds[red_cmpd] = Compound(hyp_cmpd, entry, nspecies);
      }
    }

    // convert map of compounds to map of vector<vector<string>>
    for (auto& compound_pair : compounds) {
      string subtitle = "\\textbf{";
      subtitle += compound_pair.second.hyper_compound;
      if (compound_pair.second.stable) {
        subtitle += " (stable)";
      } else {
        subtitle += " (unstable) ";
        subtitle += phaseDecompCompact(compound_pair.second.entries[0]);
      }
      subtitle += "}";
      compound_data[subtitle] = compound_pair.second.entriesToTableVector();
    }

    return compound_data;
  }

  /// @brief filters entries to make plots and tables format nicely
  /// @param entries unfiltered list of entries to plot
  /// @param all_unstable true if all entries are unstable (above zero-energy tie line)
  /// @returns list of filtered entries
  /// @author
  /// mod{NHA,20260324,created}
  vector<nhull::Entry> plotFilter(const vector<nhull::Entry>& entries, bool all_unstable) {
    double lowest_energy = 0;
    double highest_energy = 0;
    vector<nhull::Entry> first_filter_entries = nhull::lightFilter(entries);
    vector<nhull::Entry> second_filter_entries;

    for (const auto& entry : first_filter_entries) {
      bool blacklist_entry = false;
      double curr_energy = entry.enthalpy_formation_atom;

      if (curr_energy > 0 && !all_unstable) {
        blacklist_entry = true;
      }
      if (all_unstable && curr_energy > 1000) {
        blacklist_entry = true;
      }
      if (entry.nspecies == 1) { //ignore all unaries
        blacklist_entry = true;
      }
      if (entry.m_is_artificial) {
        blacklist_entry = true;
      }
      if (!blacklist_entry) {
        second_filter_entries.push_back(entry);
        // update lowest energy:
        if (curr_energy < lowest_energy) {
          lowest_energy = curr_energy;
        }
        if (curr_energy > highest_energy) {
          highest_energy = curr_energy;
        }
      }
    }
    return second_filter_entries;
  }

  /// @brief xplotter implementation of plotting functions
  /// @param entries
  /// @param specie_vector
  /// @param nhull_facets pertain to facets generated by convex hull
  /// @param save_root
  /// @param plot_3d
  /// @param calc_stability_criterion
  /// @author
  /// mod{NHA,20251005,created}
  void toTexXplotter(const vector<nhull::Entry>& entries, vector<string> specie_vector, vector<nhull::NhullFacet> nhull_facets, string save_root, bool plot_3d) {
    vector<nhull::Entry> first_filter_entries = nhull::lightFilter(entries);
    vector<nhull::Entry> second_filter_entries;
    int dimension = specie_vector.size();
    int number_of_entries = first_filter_entries.size();

    // get number of ingested hull points:
    // ignore flagged entries and entries with enthalpies greater than 0
    bool all_unstable = allEntriesUnstable(entries);
    second_filter_entries = plotFilter(entries, all_unstable);

    if (dimension == 2) {
      figures::plotNhullBinary(second_filter_entries, number_of_entries, specie_vector, nhull_facets, save_root);
    } else if (dimension == 3) {
      if (plot_3d) {
        figures::plotNhullTernary3D(second_filter_entries, number_of_entries, specie_vector, nhull_facets, save_root);
      } else {
        figures::plotNhullTernary(second_filter_entries, number_of_entries, specie_vector, nhull_facets, save_root);
      }
    } else {
      // else just create a table:
      figures::plotNhullTable(second_filter_entries, number_of_entries, specie_vector, save_root);
    }
  }
} // namespace nhullxplotter
