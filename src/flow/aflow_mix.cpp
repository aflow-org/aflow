// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo - 2009 Duke

#ifndef _AFLOW_MIX_CPP
#define _AFLOW_MIX_CPP

#include <cstdlib>
#include <string>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xerror.h"

#include "aflow.h"
#include "aflow_defs.h"
#include "flow/aflow_ivasp.h"
#include "structure/aflow_xatom.h"

// ***************************************************************************
// HT miscibility wrap up
int MiscibilityCheck(int speciesA, int speciesB) {       // aflow_mix.cpp
  const std::string system = XATOM_AlphabetizationSpecies(GetAtomSymbol(speciesA), GetAtomSymbol(speciesB));
  return MiscibilityCheck(system);
}

int MiscibilityCheck(std::string speciesA, std::string speciesB) {       // aflow_mix.cpp
  const std::string system = XATOM_AlphabetizationSpecies(speciesA, speciesB);
  return MiscibilityCheck(system);
}

// ***************************************************************************
// Pauling miscibility wrap-up
int MiscibilityExperimentsCheck(int speciesA, int speciesB) {// aflow_mix.cpp
  const std::string system = XATOM_AlphabetizationSpecies(GetAtomSymbol(speciesA), GetAtomSymbol(speciesB));
  return MiscibilityExperimentsCheck(system);
}

int MiscibilityExperimentsCheck(std::string speciesA, std::string speciesB) {// aflow_mix.cpp
  const std::string system = XATOM_AlphabetizationSpecies(speciesA, speciesB);
  return MiscibilityExperimentsCheck(system);
}

// ***************************************************************************
// Miedema
int MiscibilityMiedemaCheck(int speciesA, int speciesB) {// aflow_mix.cpp
  double ratio;
  // cerr << std::endl << "DEBUG " << speciesA << " " << speciesB << std::endl;
  // cerr << "DEBUG " << GetAtomSymbol(speciesA) << " " << GetAtomSymbol(speciesB) << std::endl;
  // cerr << "DEBUG " << GetAtomName(speciesA) << " " << GetAtomName(speciesB) << std::endl;
  //  cerr << speciesA << " " << speciesB << std::endl;
  if (speciesA == speciesB) {
    return MISCIBILITY_SYSTEM_MISCIBLE; // A is obviously miscibile with itself
  }
  if (vatom_miedema_phi_star.at(speciesA) == NNN) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "vatom_miedema_phi_star.at(speciesA) undefined", _INPUT_ILLEGAL_);
  } // CO20200624
  if (vatom_miedema_phi_star.at(speciesB) == NNN) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "vatom_miedema_phi_star.at(speciesB) undefined", _INPUT_ILLEGAL_);
  } // CO20200624
  if (vatom_miedema_nws.at(speciesA) == NNN) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "vatom_miedema_nws.at(speciesA) undefined", _INPUT_ILLEGAL_);
  } // CO20200624
  if (vatom_miedema_nws.at(speciesB) == NNN) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "vatom_miedema_nws.at(speciesB) undefined", _INPUT_ILLEGAL_);
  } // CO20200624
  ratio = (vatom_miedema_phi_star.at(speciesA) - vatom_miedema_phi_star.at(speciesB)) / (vatom_miedema_nws.at(speciesA) - vatom_miedema_nws.at(speciesB));
  // cerr << "DEBUG " << vatom_miedema_phi_star.at(speciesA) << " " << vatom_miedema_phi_star.at(speciesB) << std::endl;
  // cerr << "DEBUG " << vatom_miedema_nws.at(speciesA) << " " << vatom_miedema_nws.at(speciesB) << std::endl;
  // cerr << ratio << std::endl;
  if (ratio >= MIEDEMA_MIX_SLOPE) {
    return MISCIBILITY_SYSTEM_MISCIBLE;
  } else {
    return MISCIBILITY_SYSTEM_NOMIX;
  }
  return MISCIBILITY_SYSTEM_UNKNOWN; // impossible because it is one of the other before
}

int MiscibilityMiedemaCheck(std::string speciesA, std::string speciesB) {// aflow_mix.cpp
  return MiscibilityMiedemaCheck(GetAtomNumber(speciesA), GetAtomNumber(speciesB));
}

int MiscibilityMiedemaCheck(std::string system_in) {  // (nomix,unknown,mix)
  std::string speciesA;
  std::string speciesB;
  KBIN::VASP_SplitAlloySpecies(KBIN::VASP_PseudoPotential_CleanName(system_in), speciesA, speciesB);
  return MiscibilityMiedemaCheck(GetAtomNumber(speciesA), GetAtomNumber(speciesB));
}

// ***************************************************************************
// HumeRothery
int MiscibilityHumeRotheryCheck(int speciesA, int speciesB) {// aflow_mix.cpp
  const int delta_valence = std::abs(GetAtomValenceIupac(speciesA) - GetAtomValenceIupac(speciesB));
  const double delta_electronegativity = std::abs(GetAtomElectronegativity(speciesA) - GetAtomElectronegativity(speciesB)) / (GetAtomElectronegativity(speciesA) + GetAtomElectronegativity(speciesB));
  const double delta_radius = std::abs(GetAtomRadius(speciesA) - GetAtomRadius(speciesB)) / (GetAtomRadius(speciesA) + GetAtomRadius(speciesB));
  // cerr << delta_radius << ",";
  // cerr << delta_electronegativity << ",";
  // cerr << delta_valence << ",";
  //  cerr << GetAtomCrystal(speciesA) << "," << GetAtomCrystal(speciesB) << " ";
  if (delta_radius <= 0.15) {                                        // cut-off is 15% on radius mismatch
    if (delta_electronegativity <= 0.15) {                           // cut-off is 15% on electronegativity mismatch
      if (delta_valence <= 2) {                                      // similar valnce, max mismatch is +-2
        if (GetAtomCrystal(speciesA) == GetAtomCrystal(speciesB)) {  // same crystal structure
          return MISCIBILITY_SYSTEM_SOLUTION;
        }
      }
    }
  }

  return MISCIBILITY_SYSTEM_UNKNOWN; // impossible because it is one of the other before
}

int MiscibilityHumeRotheryCheck(std::string speciesA, std::string speciesB) {// aflow_mix.cpp
  return MiscibilityHumeRotheryCheck(GetAtomNumber(speciesA), GetAtomNumber(speciesB));
}

int MiscibilityHumeRotheryCheck(std::string system_in) {  // (nomix,unknown,mix)
  std::string speciesA;
  std::string speciesB;
  KBIN::VASP_SplitAlloySpecies(KBIN::VASP_PseudoPotential_CleanName(system_in), speciesA, speciesB);
  return MiscibilityHumeRotheryCheck(GetAtomNumber(speciesA), GetAtomNumber(speciesB));
}

// ***************************************************************************
#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
