//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2025           *
// *                  Simon Divilov - Duke University 2023                   *
// *                  Hagen Eckert - Duke University 2023                    *
// *                                                                         *
//****************************************************************************
// Written by Simon Divilov and Hagen Eckert, 2023

#ifndef _AFLOW_FIGURES_H_
#define _AFLOW_FIGURES_H_

#include <string>
#include <vector>

#include "AUROSTD/aurostd_xoption.h"

#include "modules/HULL/aflow_nhull_entry.h"
#include "modules/HULL/aflow_nhull_facet.h"

namespace figures {
  void plotDOS(const aurostd::xoption&);
  void plotNhullTernary(const std::vector<nhull::Entry>& second_filter_entries, int number_of_entries, const std::vector<std::string>& specie_vector, std::vector<nhull::NhullFacet>& nhull_facets, const std::string& save_root);
  void plotNhullBinary(const std::vector<nhull::Entry>& second_filter_entries, int number_of_entries, const std::vector<std::string>& specie_vector, std::vector<nhull::NhullFacet>& nhull_facets, std::string save_root);
  void plotNhullTernary3D(const std::vector<nhull::Entry>& second_filter_entries, int number_of_entries, const std::vector<std::string>& specie_vector, std::vector<nhull::NhullFacet>& nhull_facets, std::string save_root);
  void plotNhullTable(const std::vector<nhull::Entry>& second_filter_entries, int number_of_entries, const std::vector<std::string>& specie_vector, std::string save_root);
} // namespace figures

#endif
