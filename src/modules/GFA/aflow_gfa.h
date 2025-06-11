// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *           Aflow ERIC PERIM MARTINS - Duke University 2014-2016          *
// *           Aflow DENISE FORD - Duke University 2016-2019                 *
// *                                                                         *
// ***************************************************************************
// Original code written by Eric Perim Martins
// Revised and expanded by Denise Ford

#ifndef AFLOW_GFA_H
#define AFLOW_GFA_H

#include <string>

#include "AUROSTD/aurostd_xoption.h"

namespace pflow {

  // Calculation functions

  void CalculateGFA(aurostd::xoption& vpflow, const std::string& input, const std::string& AE_file_read, double fe_cut);

} // namespace pflow

// ***************************************************************************

#endif
