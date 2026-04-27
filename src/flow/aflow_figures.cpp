//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2025           *
// *                  Simon Divilov - Duke University 2023                   *
// *                  Hagen Eckert - Duke University 2023                    *
// *                                                                         *
//****************************************************************************
// Written by Simon Divilov and Hagen Eckert, 2023

#include "aflow_figures.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <deque>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <istream>
#include <iterator>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xfile.h"
#include "AUROSTD/aurostd_xoption.h"

#include "../AUROSTD/aurostd_xplotter.h"
#include "../modules/HULL/aflow_nhull_latex.h"
#include "aflow.h"
#include "aflow_defs.h"
#include "aflow_ivasp.h"
#include "modules/HULL/aflow_nhull_entry.h"
#include "modules/HULL/aflow_nhull_facet.h"
#include "structure/aflow_xatom.h"
#include "structure/aflow_xstructure.h"

using namespace std;
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
using std::vector;

namespace figures {
  void plotDOS(const aurostd::xoption& plot_options) {
    static const std::array<string, 16> ORBITALS_LM_LABELS = {
        "s", "p_y", "p_z", "p_x", "d_{xy}", "d_{yz}", "d_{z^2}", "d_{xz}", "d_{x^2-y^2}", "f_{y(3x^2-y^2)}", "f_{xyz}", "f_{yz^2}", "f_{z^3}", "f_{xz^2}", "f_{z(x^2-y^2)}", "f_{x(x^2-3y^2)}"};
    static const std::array<string, 4> ORBITALS_LABELS = {"s", "p", "d", "f"};
    static const std::array<string, 20> ROMAN_NUMERALS = {"I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII", "XVIII", "XIX", "XX"};
    string directory;
    if (!plot_options.getattachedscheme("PLOT::DIRECTORY").empty()) {
      directory = aurostd::CleanFileName(plot_options.getattachedscheme("PLOT::DIRECTORY"));
    } else {
      directory = aurostd::getPWD();
    }
    xDOSCAR xdos;
    if (!xdos.initialize(directory + "/DOSCAR.static")) {
      string message = "Missing DOSCAR.static";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
    xstructure xstr(directory + "/POSCAR.static", IOVASP_POSCAR);
    if (xstr.species.empty()) {
      string message = "Missing POSCAR.static";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
    if (xstr.is_vasp4_poscar_format) {
      vector<string> atom_names = KBIN::ExtractAtomicSpecies(directory);
      deque<_atom> atoms = xstr.atoms;
      for (size_t i = 0; i < atoms.size(); i++) {
        atoms[i].name_is_given = true;
        atoms[i].name = atom_names[i];
        atoms[i].cleanname = atom_names[i];
      }
      xstr.ReplaceAtoms(atoms);
    }
    bool lmResolved = plot_options.flag("PLOT::LMRESOLVED") && xdos.lmResolved;
    long int proj_index = aurostd::string2utype<long int>(plot_options.getattachedscheme("PLOT::PROJ_INDEX"));
    if (proj_index < 0) {
      string message = "Projection index cannot be less than 0";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }
    deque<deque<deque<deque<double>>>> vDOS;
    string projection = aurostd::toupper(plot_options.getattachedscheme("PLOT::PROJECTION"));
    vector<string> proj_labels = {""};
    if (projection == "ATOM") {
      if (lmResolved) {
        vDOS = xdos.vDOS_lm_atom;
      } else {
        vDOS = xdos.vDOS_atom;
      }
      if (proj_index > (long int) vDOS.size() - 1) {
        stringstream message;
        message << "Projection index greater than number of atoms (" << vDOS.size() - 1 << ")";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
      }
      for (size_t iat = 0; iat < xstr.atoms.size(); iat++) {
        proj_labels.push_back(xstr.atoms[iat].cleanname + "(" + aurostd::utype2string<uint>(iat + 1) + ")");
      }
    } else if (projection == "IATOM") {
      xstr.sortAtomsEquivalent();
      xdos.GetVDOSIAtom(xstr);
      if (lmResolved) {
        vDOS = xdos.vDOS_lm_iatom;
      } else {
        vDOS = xdos.vDOS_iatom;
      }
      if (proj_index > (long int) vDOS.size() - 1) {
        stringstream message;
        message << "Projection index greater than number of inequivalent atoms (" << vDOS.size() - 1 << ")";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
      }
      vector<string> atom_names;
      vector<size_t> index;
      for (size_t iat = 0; iat < xstr.iatoms.size(); iat++) {
        atom_names.push_back(xstr.atoms[xstr.iatoms[iat][0]].cleanname);
        aurostd::WithinList(atom_names, atom_names.back(), index);
        proj_labels.push_back(atom_names.back() + "(" + ROMAN_NUMERALS[index.size() - 1] + ")");
      }
    } else if (projection == "SPECIES") {
      xdos.GetVDOSSpecies(xstr);
      if (lmResolved) {
        vDOS = xdos.vDOS_lm_species;
      } else {
        vDOS = xdos.vDOS_species;
      }
      if (proj_index > (long int) vDOS.size() - 1) {
        stringstream message;
        message << "Projection index greater than number of species (" << vDOS.size() - 1 << ")";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
      }
      for (size_t isp = 0; isp < xstr.species.size(); isp++) {
        proj_labels.push_back(xstr.species[isp]);
      }
    } else if (projection == "NONE") {
      vDOS.push_back(xdos.vDOS_atom[0]);
    } else {
      string message = "Invalid projection option";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    uint nskip = aurostd::string2utype<uint>(plot_options.getattachedscheme("PLOT::ENERGY_SKIP"));
    if (nskip < 1) {
      nskip = 1;
    }
    vector<double> energy;
    double efermi = 0.0;
    double energy_max = 0.0;
    double energy_min = 0.0;
    if (plot_options.flag("PLOT::NOSHIFT")) {
      energy = aurostd::getEveryNth(aurostd::deque2vector(xdos.venergy), nskip);
      efermi = xdos.Efermi;
      energy_max = xdos.energy_max;
      energy_min = xdos.energy_min;
    } else {
      energy = aurostd::getEveryNth(aurostd::deque2vector(xdos.venergyEf), nskip);
      energy_max = xdos.energy_max - xdos.Efermi;
      energy_min = xdos.energy_min - xdos.Efermi;
    }
    aurostd::xplotter xplt;
    xplt.axis_properties["xmax"] = aurostd::utype2string<double>(energy_max);
    xplt.axis_properties["xmin"] = aurostd::utype2string<double>(energy_min);
    size_t nspin = vDOS[0][0].size();
    if (nspin == 1) {
      xplt.axis_properties["ymin"] = "0.0";
    }
    vector<double> dos;
    if (projection == "NONE") {
      for (size_t ispin = 0; ispin < nspin; ispin++) {
        dos = aurostd::getEveryNth(aurostd::deque2vector(vDOS[0][0][ispin]), nskip);
        if (ispin) {
          std::transform(dos.begin(), dos.end(), dos.begin(), [](double& x) { return -x; });
        }
        xplt.plot(energy, dos,
                  {
                      {"label", (!ispin) ? "Total" : ""},
                      {"color",               "cyclic0"}
        });
      }
    } else {
      size_t norb = 1;
      vector<string> orb_labels;
      if (plot_options.flag("PLOT::ORBITALS")) {
        norb = vDOS[0].size();
        orb_labels.insert(orb_labels.end(), std::begin(ORBITALS_LABELS), std::end(ORBITALS_LABELS));
      } else if (lmResolved) {
        norb = vDOS[0].size();
        orb_labels.insert(orb_labels.end(), std::begin(ORBITALS_LM_LABELS), std::end(ORBITALS_LM_LABELS));
      }
      for (size_t iproj = proj_index; iproj < vDOS.size(); iproj++) {
        for (size_t iorb = 0; iorb < norb; iorb++) {
          string label;
          if (iorb == 0) {
            label = "Total " + proj_labels[iproj];
          } else {
            label = orb_labels[iorb - 1];
          }
          for (size_t ispin = 0; ispin < nspin; ispin++) {
            dos = aurostd::getEveryNth(aurostd::deque2vector(vDOS[iproj][iorb][ispin]), nskip);
            if (ispin) {
              std::transform(dos.begin(), dos.end(), dos.begin(), [](double& x) { return -x; });
            }
            xplt.plot(energy, dos,
                      {
                          {"label",                  (!ispin) ? label : ""},
                          {"color", "cyclic" + aurostd::utype2string(iorb)}
            });
          }
        }
        if (proj_index) {
          break;
        }
      }
    }
    xplt.vline(efermi, {
                           {"line style", "dashed"},
                           {     "color",   "gray"}
    });
    string title = xdos.title;
    aurostd::StringSubst(title, "_", "\\_");
    xplt.axis_properties["title"] = title;
    xplt.axis_properties["xlabel"] = "Energy";
    xplt.axis_properties["x unit"] = "eV";
    xplt.axis_properties["ylabel"] = "Density of states";
    xplt.axis_properties["y unit"] = "1/eV";
    xplt.axis_properties["legend pos"] = "north east";
    cerr << xplt.toString() << endl;
  }

  /// @brief plotting routine for ternary convex hull systems
  /// @param second_filter_entries list of entries filtered specifically for plotting
  /// @param number_of_entries total number of entries
  /// @param specie_vector the global user-submitted compound
  /// @param nhull_facets hull facets
  /// @authors
  /// mod{NHA,20260206,created}
  void plotNhullTernary(const vector<nhull::Entry>& second_filter_entries, const int number_of_entries, const vector<string>& specie_vector, vector<nhull::NhullFacet>& nhull_facets, const string& save_root) {
    string plot_path = save_root + "_plot.pdf";
    string table_path = save_root + "_table.pdf";
    aurostd::xplotter xplt;
    string title = specieString(specie_vector) + " entries: " + to_string(number_of_entries);

    plotNhullTable(second_filter_entries, number_of_entries, specie_vector, save_root); //create table here

    //added interpolated patches:
    nhullxplotter::addFacetInterpolate(nhull_facets, specie_vector, xplt);

    //add facet edges:
    nhullxplotter::addFacetEdges(nhull_facets, specie_vector, xplt);

    xplt.axis_properties["title"] = title;
    xplt.axis_properties["xlabel"] = specie_vector[0];
    xplt.axis_properties["ylabel"] = specie_vector[1];
    xplt.axis_properties["zlabel"] = specie_vector[2];
    xplt.axis_properties["ylabel style"] = "{font=\\large,at={(axis cs:0,1,0)},anchor=north east,opacity=1}";
    xplt.axis_properties["xlabel style"] = "{font=\\large,at={(axis cs:1,0,0)},anchor=south,opacity=1}";
    xplt.axis_properties["zlabel style"] = "{font=\\large,at={(axis cs:0,0,1)},anchor=north west,opacity=1}";
    xplt.axis_properties["label style"] = "{}";
    xplt.axis_properties["ticks"] = "none";
    xplt.axis_properties["colorbar style"] = "{ylabel={$H_{\\mathrm{f}}$ (meV/atom)}}";
    xplt.axis_properties["colormap name"] = "{mymap}";
    xplt.axis_properties["colormap"] = "{mymap}{rgb(0pt)=(0,0,1); rgb(63pt)=(1,0.644,0)}";
    xplt.axis_properties["clip"] = "false";
    xplt.nodes = nhullxplotter::addNodes3D(second_filter_entries, specie_vector);

    xplt.save(plot_path);
  }

  /// @brief plotting routine for binary convex hull systems
  /// @param second_filter_entries list of entries filtered specifically for plotting
  /// @param number_of_entries total number of entries
  /// @param specie_vector the global user-submitted compound
  /// @param nhull_facets hull facets
  /// @authors
  /// mod{NHA,20260206,created}
  void plotNhullBinary(const vector<nhull::Entry>& second_filter_entries, const int number_of_entries, const vector<string>& specie_vector, vector<nhull::NhullFacet>& nhull_facets, string save_root) {
    string plot_path = save_root + "_plot.pdf";
    string table_path = save_root + "_table.pdf";
    aurostd::xplotter xplt;
    string title = specieString(specie_vector) + " entries: " + to_string(number_of_entries);

    plotNhullTable(second_filter_entries, number_of_entries, specie_vector, save_root); //create table here

    //retrieve plot points:
    auto add_store = nhullxplotter::addHullMarks(specie_vector, second_filter_entries);
    vector<double> x = get<0>(add_store);
    vector<double> y = get<1>(add_store);

    //add facet edges:
    nhullxplotter::addFacetEdges(nhull_facets, specie_vector, xplt);

    //xplt.axis_properties["colorbar"] = "false";
    xplt.axis_properties["ymax"] = "0";
    xplt.axis_properties["xmin"] = "0";
    xplt.axis_properties["xmax"] = "1";
    xplt.axis_properties["title"] = title;
    xplt.axis_properties["xtick"] = "{1,0.8,0.6,0.4,0.2,0}";
    xplt.axis_properties["xticklabels"] = "{{\\LARGE " + specie_vector[1] + "},0.8,0.6,0.4,0.2,{\\LARGE " + specie_vector[0] + "}}";
    xplt.axis_properties["ylabel"] = "{$H_{\\mathrm{f}}$}";
    xplt.axis_properties["y unit"] = "meV/atom";
    xplt.axis_properties["legend pos"] = "north east";
    xplt.axis_properties["colormap name"] = "{mymap}";
    xplt.axis_properties["colormap"] = "{mymap}{rgb(0pt)=(0,0,1); rgb(63pt)=(1,0.644,0)}";
    xplt.nodes = nhullxplotter::addNodes2D(second_filter_entries, specie_vector);

    std::map<std::string, std::string> properties = {
        {"colormap name", "mymap"}
    };

    xplt.scatter(x, y, y, properties);
    xplt.save(plot_path);
  }

  /// @brief plotting routine for ternary convex hull systems
  /// @param second_filter_entries list of entries filtered specifically for plotting
  /// @param number_of_entries total number of entries
  /// @param specie_vector the global user-submitted compound
  /// @param nhull_facets hull facets
  /// @authors
  /// mod{NHA,20260206,created}
  void plotNhullTernary3D(const vector<nhull::Entry>& second_filter_entries, const int number_of_entries, const vector<string>& specie_vector, vector<nhull::NhullFacet>& nhull_facets, string save_root) {
    string plot_path = save_root + "_3Dplot.tex";
    string table_path = save_root + "_table.pdf";

    string store = ternary::createTernary3D(specie_vector, nhull_facets, second_filter_entries, number_of_entries);
    plotNhullTable(second_filter_entries, number_of_entries, specie_vector, save_root); //create table here

    aurostd::string2file(store, plot_path);
  }

  /// @brief creates a table of data for a convex hull run
  /// @param second_filter_entries list of entries filtered specifically for displaying
  /// @param number_of_entries total number of entries
  /// @param specie_vector the global user-submitted compound
  /// @authors
  /// mod{NHA,20260206,created}
  void plotNhullTable(const vector<nhull::Entry>& second_filter_entries, int number_of_entries, const vector<string>& specie_vector, string save_root) {
    string table_path = save_root + "_table.pdf";

    //table data:
    string title = specieString(specie_vector) + " entries: " + to_string(number_of_entries);
    std::map<string, vector<map<string, vector<string>>>> compound_data = nhullxplotter::entriesToCmpdData(second_filter_entries);
    vector<string> column_titles = {"auid", "prototype", "{$H_{\\mathrm{f}}$}", "{$\\Delta H_{\\mathrm{hull}}$}"};
    vector<string> column_units = {"", "", "", ""};
    int font_size = 10;
    auto table = aurostd::xtable(compound_data, column_titles, column_units, font_size, title);

    table.save(table_path);
  }
} // namespace figures
