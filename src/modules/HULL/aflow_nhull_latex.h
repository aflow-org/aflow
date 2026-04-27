// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by Nicholas Anderson
// nicholas.anderson@duke.edu

#ifndef AFLOW_LATEX_H
#define AFLOW_LATEX_H

#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "AUROSTD/aurostd_xhttp.h"
#include "AUROSTD/aurostd_xplotter.h"

#include "aflow_nhull_facet.h"
#include "modules/HULL/aflow_nhull_entry.h"

// main function creates .tex file from Entry inputs and add functions
void toTex(const std::vector<nhull::Entry>& entries, std::vector<std::string> specie_vector, std::vector<nhull::NhullFacet> nhull_facets);

// utility functions
std::string specieString(const std::vector<std::string>& specie_vector);
std::string vectorString(std::vector<double> specie_vector);

namespace nhullxplotter {
  std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> addHullMarks(std::vector<std::string> specie_vector, const std::vector<nhull::Entry>& entries);
  void addFacetEdges(std::vector<nhull::NhullFacet>& nhull_facets, std::vector<std::string> specie_vector, aurostd::xplotter& xplt);
  void addFacetInterpolate(std::vector<nhull::NhullFacet>& nhull_facets, std::vector<std::string> specie_vector, aurostd::xplotter& xplt);
  std::vector<std::string> addNodes3D(const std::vector<nhull::Entry>& filtered_entries, const std::vector<std::string>& specie_vector);
  std::vector<std::string> addNodes2D(const std::vector<nhull::Entry>& filtered_entries, const std::vector<std::string>& specie_vector);
  std::vector<std::vector<std::string>> entriesToVector(const std::vector<nhull::Entry>& second_filter_entries);
  std::map<std::string, std::vector<std::map<std::string, std::vector<std::string>>>> entriesToCmpdData(const std::vector<nhull::Entry>& filtered_entries);
  void toTexXplotter(const std::vector<nhull::Entry>& entries, std::vector<std::string> specie_vector, std::vector<nhull::NhullFacet> nhull_facets, std::string save_root, bool plot_3d);
} // namespace nhullxplotter

namespace ternary {
  std::string createTernary3D(std::vector<std::string> specie_vector, std::vector<nhull::NhullFacet>& nhull_facets, const std::vector<nhull::Entry>& entries, int number_of_entries);
}

#endif // AFLOW_LATEX_H
