
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <ios>
#include <iostream>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xscalar.h"

#include "aflow.h"
#include "aflow_aflowrc.h"
#include "aflow_defs.h"
#include "aflow_xhost.h"
#include "extern/SYMBOLICCPLUSPLUS/symbolic.h"
#include "flow/aflow_ivasp.h"
#include "flow/aflow_pflow.h"
#include "interfaces/aflow_symbolic.h"  //ME20220124
#include "modules/SYM/aflow_symmetry_spacegroup.h" //DX20180723
#include "modules/SYM/aflow_wyckoff.h"
#include "structure/aflow_xatom.h"
#include "structure/aflow_xstructure.h"

using std::ios_base;
using std::setw;

// **************************************************************************
// PAULING DETECTOR
// **************************************************************************
bool PAULING_WyckoffDetector(const vector<string>& vinput) {
  vector<string> tokens;
  vector<string> elements;
  aurostd::string2tokens(vinput.at(0), tokens);
  for (size_t i = 1; i < tokens.size(); i++) {
    elements.push_back(tokens.at(i));
  }
  for (size_t i = 0; i < elements.size(); i++) {
    //  cout << "scanning " << elements.at(i) << endl;
    for (size_t j = 1; j < vinput.size(); j++) {
      aurostd::string2tokens(vinput[j], tokens);
      if (tokens.size() == 7) {
        if (tokens.at(1) == elements.at(i)) {
          cout << "   ";
          for (uint k = 0; k < 3; k++) {
            cout << aurostd::PaddedPOST(tokens.at(4 + k), 7, " ") << "  ";
          }
          cout << char(i + 65) << "   "; // shift to ascii
          cout << aurostd::PaddedPRE(tokens.at(0), 7, " ") << "  ";
          for (uint k = 1; k < 4; k++) {
            cout << aurostd::PaddedPOST(tokens.at(k), 4, " ") << "  ";
          }
          cout << endl;
        }
      }
    }
  }

  throw aurostd::xerror(__AFLOW_FILE__, XPID + "PAULING_WyckoffDetector():", "Throw for debugging purposes.", _GENERIC_ERROR_);
}

string DeStupidizer(string& strin) {
  aurostd::StringSubstInPlace(strin, "\t", " ");
  string s = " ";
  for (char c = 1; c < 32; c++) {
    s[0] = c;
    aurostd::StringSubstInPlace(strin, s, " ");
  } // destupidization
  aurostd::StringSubstInPlace(strin, "  ", " "); // cleaning
  return strin;
}

// **************************************************************************
// Xstructure operator>>  INPUT_XSTRUCTURE_INPUT
// **************************************************************************
// loads istream into a xstructure: istream >> xstructure
istream& operator>>(istream& cinput, xstructure& a) {
  // this is also a constructor so everything should look well defined
  const bool LDEBUG = (false || XHOST.DEBUG);
  stringstream message;

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
  }
  if (LDEBUG) {
    if (a.iomode == IOAFLOW_AUTO) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOAFLOW_AUTO" << endl;
    }
    if (a.iomode == IOVASP_AUTO) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOVASP_AUTO" << endl;
    }
    if (a.iomode == IOVASP_POSCAR) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOVASP_POSCAR" << endl;
    }
    if (a.iomode == IOVASP_ABCCAR) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOVASP_ABCCAR" << endl;
    }
    if (a.iomode == IOVASP_WYCKCAR) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOVASP_WYCKCAR" << endl;
    }
    if (a.iomode == IOQE_AUTO) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOQE_AUTO" << endl;
    }
    if (a.iomode == IOQE_GEOM) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOQE_GEOM" << endl;
    }
    if (a.iomode == IOAIMS_AUTO) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOAIMS_AUTO" << endl; // CO20171008
    }
    if (a.iomode == IOAIMS_GEOM) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOAIMS_GEOM" << endl; // CO20171008
    }
    if (a.iomode == IOABINIT_GEOM) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOABINIT_GEOM" << endl; // DX20200310
    }
    if (a.iomode == IOELK_GEOM) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOELK_GEOM" << endl; // DX20200310
    }
    if (a.iomode == IOCIF) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOCIF" << endl; // DX20180723
    }
    if (a.iomode == IOATAT_STR) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOATAT_STR" << endl; // SD20220114
    }
    if (a.iomode == IOLMP_DATA) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOLMP_DATA" << endl; // SD20240111
    }
  }

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " definitions" << endl;
  }
  uint iline = 0;
  vector<string> vinput;
  vector<string> tokens;
  aurostd::stream2vectorstring(cinput, vinput);
  // CO20180702 - detect NO input
  string input_no_spaces = aurostd::joinWDelimiter(vinput, "");
  input_no_spaces = aurostd::RemoveWhiteSpaces(input_no_spaces);
  if (input_no_spaces.empty()) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No input", _INPUT_MISSING_); // CO20190629
  } // CO20180702

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " DeStupidizer" << endl;
  }
  // now clean for comments, tabs, double spaces ... etc
  // CO20180409 - fixing for issues with # at the beginning of the line
  string::size_type loc; // CO20180409
  for (size_t i = 1; i < vinput.size(); i++) {
    if (i > 0) {
      // CO20180409
      loc = vinput[i].find("//");
      vinput[i] = vinput[i].substr(0, loc);
      loc = vinput[i].find('#');
      vinput[i] = vinput[i].substr(0, loc);
      loc = vinput[i].find("!");
      vinput[i] = vinput[i].substr(0, loc);
    }
    DeStupidizer(vinput[i]);
  }
  if (vinput.empty()) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No input", _INPUT_MISSING_); // CO20190629
  } // CO20180420

  if (LDEBUG) {
    // CO20230502
    for (size_t i = 0; i < vinput.size(); i++) {
      cerr << __AFLOW_FUNC__ << " vinput[i=" << i << "]=\"" << vinput[i] << "\"" << endl;
    }
  }

  //  for(size_t i=0;i<vinput.size();i++) cerr << "[" << i << "] " <<  vinput[i] << " " << "[X]" << endl;   throw aurostd::xerror(__AFLOW_FILE__,XPID+"sortAtomsTypes():","Throw for debugging purposes.",_GENERIC_ERROR_);
  string stmp;
  bool IOMODE_found = false;
  vector<double> poccaus; // partial occupation local host
  a.partial_occupation_sublattice.clear(); // partial occupation local host

  // CO20190219 - need to clear atoms
  // a.atoms.clear();  //CO20190219 - need to use RemoveAtom() (safe)
  // CO20190219 - really remove atoms
  for (size_t i = a.atoms.size() - 1; i < a.atoms.size(); i--) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " removing atom[" << i << "]" << endl;
    }
    a.RemoveAtom(i);
  }

  // ----------------------------------------------------------------------
  // SEARTH FOR MODES
  //  LDEBUG=true;
  // PAULING PROTO DETECTOR
  if (!IOMODE_found) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PAULING PROTO DETECTOR" << endl;
    }
    aurostd::string2tokens(vinput.at(0), tokens);
    if (!tokens.empty() && (tokens.at(0) == "PAULING" || tokens.at(0) == "pauling" || tokens.at(0) == "Pauling")) {
      // CO20180420 - ask for size() to print core file
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " PAULING PROTO DETECTOR = true" << endl;
      }
      IOMODE_found = true;
      PAULING_WyckoffDetector(vinput);
      tokens.clear();
    }
  }

  // ----------------------------------------------------------------------
  // QUANTUM ESPRESSO FINDER
  if (!IOMODE_found) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " QUANTUM ESPRESSO DETECTOR" << endl;
    }
    uint QE = 0;
    bool QE_ERROR = false;
    if (LDEBUG) {
      for (size_t i = 0; i < vinput.size(); i++) {
        cerr << vinput[i] << endl;
      }
    }
    for (size_t i = 0; i < vinput.size(); i++) {
      QE += aurostd::substring2bool(vinput[i], "&system", true) + aurostd::substring2bool(vinput[i], "&SYSTEM", true); // DX20180123 - added true to clean the spaces in string
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " QUANTUM ESPRESSO DETECTOR QE(&system)=" << QE << endl;
    }
    for (size_t i = 0; i < vinput.size(); i++) {
      QE += aurostd::substring2bool(vinput[i], "ibrav=", true) + aurostd::substring2bool(vinput[i], "IBRAV=", true); // DX20180123 - added true to clean the spaces in string
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " QUANTUM ESPRESSO DETECTOR QE(ibrav)=" << QE << endl;
    }
    for (size_t i = 0; i < vinput.size(); i++) {
      QE += aurostd::substring2bool(vinput[i], "nat=", true) + aurostd::substring2bool(vinput[i], "NAT=", true); // DX20180123 - added true to clean the spaces in string
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " QUANTUM ESPRESSO DETECTOR QE(nat)=" << QE << endl;
    }
    for (size_t i = 0; i < vinput.size(); i++) {
      QE += aurostd::substring2bool(vinput[i], "ntyp=", true) + aurostd::substring2bool(vinput[i], "NTYP=", true); // DX20180123 - added true to clean the spaces in string
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " QUANTUM ESPRESSO DETECTOR QE(ntyp)=" << QE << endl;
    }
    for (size_t i = 0; i < vinput.size(); i++) {
      QE += aurostd::substring2bool(vinput[i], "atomic_positions", true) + aurostd::substring2bool(vinput[i], "ATOMIC_POSITIONS", true);
    }
    // DX20180123 - added true to clean the spaces in string
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " QUANTUM ESPRESSO DETECTOR QE(ATOMIC_POSITIONS)=" << QE << endl;
    }

    for (size_t i = 0; i < vinput.size() && QE == 5; i++) {
      if (aurostd::substring2bool(vinput[i], "ATOMIC_POSITIONS") && !aurostd::substring2bool(vinput[i], "crystal", "CRYSTAL") && !aurostd::substring2bool(vinput[i], "bohr", "BOHR") && // DX added
          !aurostd::substring2bool(vinput[i], "angstrom", "ANGSTROM")) {
        cerr << __AFLOW_FUNC__ << " QE input(1) not supported vinput.at(" << i << ")= \"" << vinput[i] << "\"" << endl;
        QE_ERROR = true;
      }
      if (aurostd::substring2bool(vinput[i], "&system")) {
        //	cerr << __AFLOW_FUNC__ << " QE input(2) not supported vinput.at(" << i << ")= \"" << vinput[i] << "\"" << endl;QE_ERROR=true;
      }
    }
    if (QE == 5 && QE_ERROR) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "QE input errors", _INPUT_MISSING_); // CO20190629
    }
    if (QE == 5 && !QE_ERROR) {
      a.iomode = IOQE_AUTO; // might need further discipline but for now it is ok.. 2013 May SC
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " QUANTUM ESPRESSO DETECTOR = true" << endl;
      }
      IOMODE_found = true;
    }
  }

  // ----------------------------------------------------------------------
  // for CIF input //DX20180723 - add cif reader - START
  if (!IOMODE_found) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CIF DETECTOR" << endl;
    }
    uint CIF = 0;
    if (LDEBUG) {
      for (size_t i = 0; i < vinput.size(); i++) {
        cerr << vinput[i] << endl;
      }
    }
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(vinput[i], "loop_", true)) {
        CIF += 1;
        break;
      }
    }
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(vinput[i], "_atom_site_", true)) {
        CIF += 1;
        break;
      }
    }
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(vinput[i], "_cell_length_", true)) {
        CIF += 1;
        break;
      }
    }
    if (CIF == 3) {
      a.iomode = IOCIF;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " CIF DETECTOR = true" << endl;
      }
      IOMODE_found = true;
    }
  }
  // DX20180723 - add cif reader - END

  // ----------------------------------------------------------------------
  // ABINIT input - START (DX20200310)
  // based on documentation found at https://docs.abinit.org/variables/basic/
  // NOTE: originally had acell and rprim as keywords, but they are not mandatory
  if (!IOMODE_found) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ABINIT DETECTOR" << endl;
    }
    uint ABINIT = 0;
    // find number of atoms (natom)
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(vinput[i], "natom", true)) {
        ABINIT += 1;
        break;
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ABINIT DETECTOR (natom)=" << ABINIT << endl;
    }
    // find types of atoms (typat)
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(vinput[i], "typat", true)) {
        ABINIT += 1;
        break;
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ABINIT DETECTOR (typat)=" << ABINIT << endl;
    }
    // find atom positions (xred, xcart, xangst)
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(vinput[i], "xred", true) || aurostd::substring2bool(vinput[i], "xcart", true) || aurostd::substring2bool(vinput[i], "xangst", true)) {
        ABINIT += 1;
        break;
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ABINIT DETECTOR (xred,xcart,xangst)=" << ABINIT << endl;
    }

    if (ABINIT == 3) {
      a.iomode = IOABINIT_GEOM;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " ABINIT DETECTOR = true" << endl;
      }
      IOMODE_found = true;
    }
  }
  // ABINIT input - END (DX20200310)

  // ----------------------------------------------------------------------
  // ELK input - START (DX20200310)
  if (!IOMODE_found) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ELK DETECTOR" << endl;
    }
    uint ELK = 0;
    // find atoms keyword (atoms)
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(vinput[i], "atoms", true)) {
        ELK += 1;
        break;
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ELK DETECTOR (atoms)=" << ELK << endl;
    }
    // find lattice keyword (avec)
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(vinput[i], "avec", true)) {
        ELK += 1;
        break;
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ELK DETECTOR (avec)=" << ELK << endl;
    }

    if (ELK == 2) {
      a.iomode = IOELK_GEOM;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " ELK DETECTOR = true" << endl;
      }
      IOMODE_found = true;
    }
  }
  // ELK input - END (DX20200310)

  // ----------------------------------------------------------------------
  // ATAT input - START (SD20220117)
  if (!IOMODE_found) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ATAT DETECTOR" << endl;
    }
    uint ATAT = 1;
    // count number of entries for axes and fractional cell vectors, check for correct type
    size_t line = 0;
    for (; line < vinput.size() - 1 && line < 6 && ATAT; line++) {
      aurostd::string2tokens(vinput[line], tokens, " ");
      ATAT = (tokens.size() == 3 && aurostd::isfloat(tokens[0]) && aurostd::isfloat(tokens[1]) && aurostd::isfloat(tokens[2]));
    }
    // count number of entries for atom positions, check for correct type
    for (; line < vinput.size() - 1 && ATAT; line++) {
      aurostd::string2tokens(vinput[line], tokens, " ");
      ATAT = (tokens.size() == 4 && aurostd::isfloat(tokens[0]) && aurostd::isfloat(tokens[1]) && aurostd::isfloat(tokens[2]));
      if (ATAT && !xelement::xelement::isElement(tokens[3])) {
        ATAT = 0;
      }
    }
    if (ATAT == 1) {
      a.iomode = IOATAT_STR;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " ATAT DETECTOR = true" << endl;
      }
      IOMODE_found = true;
    }
  }
  // ATAT input - END (SD20220117)

  // ----------------------------------------------------------------------
  // LMP input - START (SD20240111)
  if (!IOMODE_found) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " LMP DETECTOR" << endl;
    }
    for (size_t i = 0; i < vinput.size() && !IOMODE_found; i++) {
      if (aurostd::substring2bool(vinput[i], "LAMMPS data file via write_data")) {
        a.iomode = IOLMP_DATA;
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " LMP DETECTOR = true" << endl;
        }
        IOMODE_found = true;
      }
    }
  }
  // LMP input - END (SD20240111)

  // ----------------------------------------------------------------------
  // for AIMS input - unfortunately, it's very generic so leave for last
  if (!IOMODE_found) {
    vector<string> tokens_line;
    // since it's so generic, let's be super strict, only look at the first word in the line
    // needs to match atom, atom_frac, or lattice
    for (size_t i = 0; i < vinput.size() && !IOMODE_found; i++) {
      aurostd::string2tokens(vinput[i], tokens_line, " ");
      if (!tokens_line.empty()) {
        IOMODE_found = (IOMODE_found || tokens_line[0] == "atom" || tokens_line[0] == "atom_frac" || tokens_line[0] == "lattice_vector");
      }
    }
    if (IOMODE_found) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " AIMS GEOM DETECTOR = true" << endl;
      }
      a.iomode = IOAIMS_GEOM;
      // if atom_frac, then lattice MUST be provided
      // if atom, no need for lattice
      bool lat_found;
      bool lat_found_anywhere;
      bool lat1_found;
      bool lat2_found;
      bool lat3_found;
      bool atom_found;
      bool atom_found_anywhere;
      bool frac_found_anywhere;
      lat_found = lat_found_anywhere = lat1_found = lat2_found = lat3_found = atom_found = atom_found_anywhere = frac_found_anywhere = false;
      for (size_t i = 0; i < vinput.size(); i++) {
        tokens_line.clear();
        aurostd::string2tokens(vinput[i], tokens_line, " ");
        if (tokens_line.empty()) {
          continue;
        }
        lat_found = (tokens_line[0] == "lattice_vector");
        lat_found_anywhere = (lat_found_anywhere || lat_found);
        frac_found_anywhere = (frac_found_anywhere || tokens_line[0] == "atom_frac");
        atom_found = (aurostd::substring2bool(tokens_line[0], "atom"));
        atom_found_anywhere = (atom_found_anywhere || atom_found || frac_found_anywhere);
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " AIMS line[" << i + 1 << "]: " << vinput[i] << endl;
          cerr << __AFLOW_FUNC__ << " AIMS lat_found=" << lat_found << endl;
          cerr << __AFLOW_FUNC__ << " AIMS lat_found_anywhere=" << lat_found_anywhere << endl;
          cerr << __AFLOW_FUNC__ << " AIMS lat1_found=" << lat1_found << endl;
          cerr << __AFLOW_FUNC__ << " AIMS lat2_found=" << lat2_found << endl;
          cerr << __AFLOW_FUNC__ << " AIMS lat3_found=" << lat3_found << endl;
          cerr << __AFLOW_FUNC__ << " AIMS atom_found=" << atom_found << endl;
          cerr << __AFLOW_FUNC__ << " AIMS atom_found_anywhere=" << atom_found_anywhere << endl;
          cerr << __AFLOW_FUNC__ << " AIMS frac_found_anywhere=" << frac_found_anywhere << endl;
        }
        if (lat_found || atom_found) {
          if (lat_found && tokens_line.size() < 4) {
            // could be more, but not less
            message << " AIMS input error, "; // CO20190629
            message << "lattice_vector "; // CO20190629
            message << "at line[" << i + 1 << "] is ill-defined" << endl; // CO20190629
            message << "line: " << vinput[i] << endl; // CO20190629
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
          }
          if (atom_found && tokens_line.size() < 5) {
            // could be more, but not less (need name/type in last column)
            message << " AIMS input error, "; // CO20190629
            message << "atom position "; // CO20190629
            message << "at line[" << i + 1 << "] "; // CO20190629
            if (tokens_line.size() == 4) {
              message << "is missing the atom name" << endl;
            } // CO20190629
            else {
              message << "is ill-defined" << endl;
            } // CO20190629
            message << "line: " << vinput[i] << endl; // CO20190629
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
          }
        }
      }
      if (!atom_found_anywhere) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "AIMS input error, no atoms found", _INPUT_ERROR_);
        // CO20190629
      }
      a.coord_flag = _COORDS_CARTESIAN_;
      if (lat_found_anywhere || frac_found_anywhere) {
        a.coord_flag = _COORDS_FRACTIONAL_;
        uint lat_count = 0;
        const vector<string> lat;
        for (size_t i = 0; i < vinput.size(); i++) {
          aurostd::string2tokens(vinput[i], tokens_line, " ");
          if (!tokens_line.empty() && tokens_line[0] == "lattice_vector") {
            lat_count++;
            if (LDEBUG) {
              cerr << __AFLOW_FUNC__ << " AIMS lat_count=" << lat_count << " line[" << i << "," << vinput[i] << "]" << endl;
            }
            if (lat_count == 1) {
              lat1_found = true;
            } else if (lat_count == 2) {
              lat2_found = true;
            } else if (lat_count == 3) {
              lat3_found = true;
            } else {
              throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "AIMS input error, too many lattice vectors found",
                                    _INPUT_ERROR_); // CO20190629
            }
          }
        }
        if (!lat1_found || !lat2_found || !lat3_found) {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "AIMS input error, incomplete lattice vector specification (needed if atom_frac found)",
                                _INPUT_ERROR_); // CO20190629
        }
      }
    } else {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " AIMS GEOM DETECTOR = false" << endl;
      }
    }
  }
  // DESPERATE FINDING => VASP
  if (!IOMODE_found) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " VASP DETECTOR" << endl;
    }
    if (a.iomode == IOAFLOW_AUTO) {
      a.iomode = IOVASP_AUTO; // still not found, need something to eat.
      //  if(a.iomode==IOVASP_AUTO) a.iomode=IOVASP_POSCAR;
      IOMODE_found = true;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " VASP DETECTOR = true" << endl;
      }
    }
  }

  // ----------------------------------------------------------------------
  // SOME EXTRA VERBOSE
  if (LDEBUG) {
    if (a.iomode == IOAFLOW_AUTO) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOAFLOW_AUTO" << endl;
    }
  }
  if (LDEBUG) {
    if (a.iomode == IOVASP_AUTO) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOVASP_AUTO" << endl;
    }
  }
  if (LDEBUG) {
    if (a.iomode == IOVASP_POSCAR) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOVASP_POSCAR" << endl;
    }
  }
  if (LDEBUG) {
    if (a.iomode == IOVASP_ABCCAR) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOVASP_ABCCAR" << endl;
    }
  }
  if (LDEBUG) {
    if (a.iomode == IOVASP_WYCKCAR) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOVASP_WYCKCAR" << endl;
    }
  }
  if (LDEBUG) {
    if (a.iomode == IOQE_AUTO) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOQE_AUTO" << endl;
    }
  }
  if (LDEBUG) {
    if (a.iomode == IOQE_GEOM) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOQE_GEOM" << endl;
    }
  }
  if (LDEBUG) {
    if (a.iomode == IOAIMS_AUTO) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOAIMS_AUTO" << endl; // CO20171008
    }
  }
  if (LDEBUG) {
    if (a.iomode == IOAIMS_GEOM) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOAIMS_GEOM" << endl; // CO20171008
    }
  }
  if (LDEBUG) {
    if (a.iomode == IOATAT_STR) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOATAT_STR" << endl; // SD20220114
    }
  }
  if (LDEBUG) {
    if (a.iomode == IOLMP_DATA) {
      cerr << __AFLOW_FUNC__ << " a.iomode = IOLMP_DATA" << endl; // SD20240111
    }
  }
  // ----------------------------------------------------------------------
  // VASP INPUT
  if (a.iomode == IOVASP_AUTO || a.iomode == IOVASP_POSCAR || a.iomode == IOVASP_ABCCAR || a.iomode == IOVASP_WYCKCAR) {
    // VASP POSCAR
    // for variable number of items
    // bool scale_second_flag=false;//,scale_third_flag=false; //CO20180409
    // double scale_second_value=0.0;//,scale_third_value=0.0;; //CO20180409
    //
    a.is_vasp4_poscar_format = false;
    a.is_vasp5_poscar_format = false;
    // -------------- POSCAR EXISTS
    int num_atoms = 0;
    // -------------- TITLE
    // input.getline(stmp,MAX_TITLE_SIZE);title=stmp;
    if (vinput.size() - 1 < iline) {
      message << "missing line[" << iline << "]" << endl; // CO20190629
      for (size_t i = 0; i < vinput.size(); i++) {
        message << vinput[i] << endl; // CO20190629
      }
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
    } // CO20180420 - check for missing lines
    aurostd::RemoveControlCodeCharactersFromString(vinput[iline++], a.title);
    // DX+ME20210525 - need to remove control code characters from input, important for web
    //  -------------- SCALE
    //     input >> a.scale;
    if (vinput.size() - 1 < iline) {
      message << "missing line[" << iline << "]" << endl; // CO20190629
      for (size_t i = 0; i < vinput.size(); i++) {
        message << vinput[i] << endl; // CO20190629
      }
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
    } // CO20180420 - check for missing lines
    stmp = vinput.at(iline++);
    aurostd::StringSubstInPlace(stmp, "\t", " ");
    aurostd::StringSubstInPlace(stmp, "  ", " ");
    aurostd::StringSubstInPlace(stmp, "  ", " ");
    aurostd::string2tokens(stmp, tokens);
    if (tokens.empty()) {
      message << "missing second line in poscar" << endl; // CO20190629
      for (size_t i = 0; i < vinput.size(); i++) {
        message << vinput[i] << endl; // CO20190629
      }
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
    }
    // oss << tokens.size() <<  " = " << tokens.at(0) << endl;throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_);
    a.scale = 0.0;
    if (!tokens.empty()) {
      a.scale = aurostd::string2utype<double>(tokens.at(0));
    }
    if (tokens.size() > 1) {
      /*a.neg_scale_second=true;*/
      a.scale_second = aurostd::string2utype<double>(tokens.at(1));
      a.neg_scale_second = std::signbit(a.scale_second);
    } // CO20180409
    if (tokens.size() > 2) {
      a.scale_third.isentry = true;
      a.scale_third.content_double = aurostd::string2utype<double>(tokens.at(2));
    } // CO20170803 - site tol
    //  oss << a.scale << endl;
    // if(a.neg_scale_second) oss << a.scale_second << endl; //CO20180409
    // if(scale_third_flag) oss << scale_third_value << endl;
    // oss << sstring << endl;
    // -------------- UNIT CELL
    if (vinput.size() - 1 < iline) {
      message << "missing line[" << iline << "]" << endl; // CO20190629
      for (size_t i = 0; i < vinput.size(); i++) {
        message << vinput[i] << endl; // CO20190629
      }
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
    } // CO20180420 - check for missing lines
    stmp = vinput.at(iline++);
    aurostd::string2tokens(stmp, tokens);
    a.iomode = IOVASP_POSCAR;
    if (tokens.size() == 3) {
      a.iomode = IOVASP_POSCAR;
    }
    if (tokens.size() == 6) {
      a.iomode = IOVASP_ABCCAR;
    }
    if (tokens.size() == 7 || tokens.size() == 8) {
      a.iomode = IOVASP_WYCKCAR;
    }
    //    oss << " LDEBUG token.size()=" << tokens.size() << "" << endl;
    // ---------------------------------------------------------------
    if (a.iomode == IOVASP_POSCAR) {
      // oss << __AFLOW_FUNC__ << " Chosen IOVASP_POSCAR" << endl;
      //    input >> a.lattice(1,1) >> a.lattice(1,2) >> a.lattice(1,3);
      a.lattice(1, 1) = aurostd::string2utype<double>(tokens[0]);
      a.lattice(1, 2) = aurostd::string2utype<double>(tokens[1]);
      a.lattice(1, 3) = aurostd::string2utype<double>(tokens[2]);
      stringstream input_tmp;
      if (vinput.size() - 1 < iline) {
        message << "missing line[" << iline << "]" << endl; // CO20190629
        for (size_t i = 0; i < vinput.size(); i++) {
          message << vinput[i] << endl; // CO20190629
        }
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
      } // CO20180420 - check for missing lines
      input_tmp.clear();
      input_tmp.str(vinput.at(iline++));
      input_tmp >> a.lattice(2, 1) >> a.lattice(2, 2) >> a.lattice(2, 3);
      if (vinput.size() - 1 < iline) {
        message << "missing line[" << iline << "]" << endl; // CO20190629
        for (size_t i = 0; i < vinput.size(); i++) {
          message << vinput[i] << endl; // CO20190629
        }
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
      } // CO20180420 - check for missing lines
      input_tmp.clear();
      input_tmp.str(vinput.at(iline++));
      input_tmp >> a.lattice(3, 1) >> a.lattice(3, 2) >> a.lattice(3, 3);
      xvector<double> data(6);
      data = Getabc_angles(a.lattice, DEGREES);
      a.a = data[1];
      a.b = data[2];
      a.c = data[3];
      a.alpha = data[4];
      a.beta = data[5];
      a.gamma = data[6];
    }
    // ---------------------------------------------------------------
    if (a.iomode == IOVASP_ABCCAR) {
      // cerr << __AFLOW_FUNC__ << " Chosen IOVASP_ABCCAR" << endl;
      a.a = aurostd::string2utype<double>(tokens[0]);
      a.b = aurostd::string2utype<double>(tokens[1]);
      a.c = aurostd::string2utype<double>(tokens[2]);
      a.alpha = aurostd::string2utype<double>(tokens[3]);
      a.beta = aurostd::string2utype<double>(tokens[4]);
      a.gamma = aurostd::string2utype<double>(tokens[5]);
      a.lattice = GetClat(a.a, a.b, a.c, a.alpha, a.beta, a.gamma);
    }
    // ---------------------------------------------------------------
    if (a.iomode == IOVASP_WYCKCAR) {
      //     cerr << __AFLOW_FUNC__ << " Chosen IOVASP_WYCKCAR" << endl;
      // DX20170905 - GET SYM_EPS
      vector<string> title_tokens;
      aurostd::string2tokens(a.title, title_tokens);
      if (aurostd::substring2bool(title_tokens, "sym_eps")) {
        a.sym_eps = aurostd::string2utype<double>(title_tokens[title_tokens.size() - 1]);
      }
      // DX20170905 - GET SYM_EPS
      a.a = aurostd::string2utype<double>(tokens[0]);
      a.b = aurostd::string2utype<double>(tokens[1]);
      a.c = aurostd::string2utype<double>(tokens[2]);
      a.alpha = aurostd::string2utype<double>(tokens[3]);
      a.beta = aurostd::string2utype<double>(tokens[4]);
      a.gamma = aurostd::string2utype<double>(tokens[5]);
      a.spacegroupnumber = aurostd::string2utype<int>(tokens[6]);
      a.spacegroupnumberoption = 0; // 1; no option
      if (tokens.size() >= 8) {
        a.spacegroupnumberoption = aurostd::string2utype<int>(tokens[7]);
      }
      a.spacegrouplabel = GetSpaceGroupLabel(a.spacegroupnumber);
      a.spacegroup = GetSpaceGroupName(a.spacegroupnumber, a.directory); // DX20180526 - add directory
      a.lattice = GetClat(a.a, a.b, a.c, a.alpha, a.beta, a.gamma);
    }
    // If scale < 0 then it should be treated as the volume.
    a.neg_scale = false;
    if (a.scale < 0.0) {
      // CO20211130 - the usual scaling factor is the lattice parameter
      // for fcc, GetVol(a.lattice)=0.25
      const double nvol = -1.0 * (a.scale);
      const double ovol = GetVol(a.lattice);
      a.scale = std::pow((nvol / ovol), 1.0 / 3.0);
      a.neg_scale = true;
    }
    a.FixLattices(); // Reciprocal/f2c/c2f
    a.kpoints_k1 = 0;
    a.kpoints_k2 = 0;
    a.kpoints_k3 = 0;
    a.kpoints_kmax = 0;
    a.kpoints_kppra = 0;
    a.kpoints_kscheme = "";
    clear(a.origin);
    // ---------------------------------------------------------------
    // -------------- CHECK VASP4 VASP5 and CHECK DIRECT/CARTESIANS STUFF
    // CO20190629 - shift to line with Direct/Cartesian to find POCC specification
    // don't worry about setting Selective Dynamics yet (but be aware of it)
    if (a.iomode == IOVASP_POSCAR) {
      stmp = vinput[6]; // true for VASP4, change if VASP5 later
    }
    if (a.iomode == IOVASP_ABCCAR) {
      stmp = vinput[4];
    }
    if (a.iomode == IOVASP_WYCKCAR) {
      stmp = vinput[4];
    }
    aurostd::StringSubstInPlace(stmp, "\t", " ");
    aurostd::StringSubstInPlace(stmp, "  ", " ");
    aurostd::StringSubstInPlace(stmp, "  ", " ");
    aurostd::string2tokens(stmp, tokens);
    if (tokens.empty()) {
      message << R"(Missing "Selective Dynamics"/"Direct"/"Cartesian" line)" << endl; // CO20190629
      for (size_t i = 0; i < vinput.size(); i++) {
        message << vinput[i] << endl; // CO20190629
      }
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_MISSING_); // CO20190629
    }

    // VASP 4
    a.is_vasp4_poscar_format = true;
    a.is_vasp5_poscar_format = false;

    if (!(tokens[0][0] == 'S' || tokens[0][0] == 's' || tokens[0][0] == 'D' || tokens[0][0] == 'd' || tokens[0][0] == 'C' || tokens[0][0] == 'c')) {
      // VASP 5
      a.is_vasp4_poscar_format = false;
      a.is_vasp5_poscar_format = true;
      // Kesong Yang identified this bug, we need to shift the line if vasp5, also if selective dynamics
      uint ijline = 0;
      if (a.iomode == IOVASP_POSCAR) {
        ijline = 7;
      }
      if (a.iomode == IOVASP_ABCCAR) {
        ijline = 5;
      }
      if (a.iomode == IOVASP_WYCKCAR) {
        ijline = 5;
      }
      if (tokens[0][0] == 'S' || tokens[0][0] == 's') {
        ijline++;
      } // selective dynamics
      stmp = vinput.at(ijline);
      aurostd::StringSubstInPlace(stmp, "\t", " ");
      aurostd::StringSubstInPlace(stmp, "  ", " ");
      aurostd::StringSubstInPlace(stmp, "  ", " ");
      aurostd::string2tokens(stmp, tokens); // now we should have D/C line
    }

    if (!(tokens[0][0] == 'S' || tokens[0][0] == 's' || tokens[0][0] == 'D' || tokens[0][0] == 'd' || tokens[0][0] == 'C' || tokens[0][0] == 'c')) {
      message << R"(Missing "Selective Dynamics"/"Direct"/"Cartesian" line)" << endl; // CO20190629
      for (size_t i = 0; i < vinput.size(); i++) {
        message << vinput[i] << endl; // CO20190629
      }
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_MISSING_); // CO20190629
    }

    a.partial_occupation_flag = false;

    for (size_t j = 1; j < tokens.size(); j++) {
      if (tokens[j][0] == 'P' || tokens[j][0] == 'p') {
        a.partial_occupation_flag = true;
      }
    }
    // if(a.scale_second==0.0) { //CO20180409
    a.partial_occupation_HNF = 0; // nothing defined //CO20180409
    a.partial_occupation_site_tol = DEFAULT_POCC_SITE_TOL;
    // DEFAULT_PARTIAL_OCCUPATION_TOLERANCE; // nothing defined  //CO20170803 - site tol
    a.partial_occupation_stoich_tol = DEFAULT_POCC_STOICH_TOL;
    // DEFAULT_PARTIAL_OCCUPATION_TOLERANCE; // nothing defined //CO20180409
    // }
    if (a.partial_occupation_flag) //&& a.neg_scale_second) //CO20180409
    {
      // CO20200106 - patching for auto-indenting
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " a.neg_scale_second=" << a.neg_scale_second << ", a.scale_second=" << a.scale_second << endl;
      }
      if (a.neg_scale_second) {
        a.partial_occupation_HNF = (int) (-a.scale_second);
      } // CO20180409
      else {
        a.partial_occupation_site_tol = a.partial_occupation_stoich_tol = a.scale_second; // CO20180409
        if (a.scale_third.isentry) {
          a.partial_occupation_stoich_tol = a.scale_third.content_double;
        }
        // CO20170803 - site tol
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " a.partial_occupation_HNF=" << a.partial_occupation_HNF << endl;
        cerr << __AFLOW_FUNC__ << " a.partial_occupation_site_tol=" << a.partial_occupation_site_tol << endl;
        cerr << __AFLOW_FUNC__ << " a.partial_occupation_stoich_tol=" << a.partial_occupation_stoich_tol << endl;
      }
      // if(a.scale_second==0.0) {//CO20180409
      //   a.partial_occupation_stoich_tol=DEFAULT_PARTIAL_OCCUPATION_TOLERANCE; // nothing defined
      //   a.partial_occupation_site_tol=DEFAULT_PARTIAL_OCCUPATION_TOLERANCE; // nothing defined  //CO20170803 - site tol
      // }
      //      if(abs(a.partial_occupation_stoich_tol)>1e-12) cerr << "a.partial_occupation_stoich_tol=" << a.partial_occupation_stoich_tol << endl; //CO20180409
      //    if(abs(a.partial_occupation_HNF)>0.0) cerr << "a.partial_occupation_HNF=" << a.partial_occupation_HNF << endl;
    }
    // if(a.partial_occupation_flag && a.scale_third.isentry){a.partial_occupation_site_tol=a.scale_third.content_double;} //CO20170803 - site tol

    // last line was last lattice vector/geometry line
    if (a.is_vasp5_poscar_format) {
      if (vinput.size() - 1 < iline) {
        message << "missing line[" << iline << "]" << endl; // CO20190629
        for (size_t i = 0; i < vinput.size(); i++) {
          message << vinput[i] << endl; // CO20190629
        }
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
      } // CO20180420 - check for missing lines
      stmp = vinput.at(iline++); // to skip toward vasp5
    }
    // ---------------------------------------------------------------
    // -------------- ATOMS
    // Number of atoms of each type and total number of atoms.
    if (vinput.size() - 1 < iline) {
      message << "missing line[" << iline << "]" << endl; // CO20190629
      for (size_t i = 0; i < vinput.size(); i++) {
        message << vinput[i] << endl; // CO20190629
      }
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
    } // CO20180420 - check for missing lines
    stmp = vinput.at(iline++);
    // The following is necessary because if the last lattice parameter has
    // tabs/spaces after it then the getline just grabs those spaces.  This
    // second getline will then be called to actually read in the number
    // of each type of atom.
    string tmpns;
    tmpns = aurostd::RemoveWhiteSpaces(stmp);
    if (string(tmpns).empty()) {
      if (vinput.size() - 1 < iline) {
        message << "missing line[" << iline << "]" << endl; // CO20190629
        for (size_t i = 0; i < vinput.size(); i++) {
          message << vinput[i] << endl; // CO20190629
        }
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
      } // CO20180420 - check for missing lines
      stmp = vinput.at(iline++);
    }
    a.num_each_type.clear();
    num_atoms = 0;
    aurostd::string2tokens(stmp, tokens);
    // -------------- check for partial occupation
    vector<string> tokens_i = tokens;
    vector<string> tokens_j;
    vector<string> tokens_k;
    int number;
    double dpocc;
    //    a.partial_occupation_flag=false;
    for (size_t i = 0; i < tokens_i.size() && !a.partial_occupation_flag; i++) {
      if (aurostd::substring2bool(tokens_i[i], "*") || aurostd::substring2bool(tokens_i[i], "+")) {
        a.partial_occupation_flag = true;
      }
    }
    // -------------- no partial occupation
    if (a.partial_occupation_flag == false) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " PARTIAL OCCUPATION = false" << endl;
      }
      for (size_t i = 0; i < tokens_i.size(); i++) {
        number = aurostd::string2utype<int>(tokens_i[i]);
        dpocc = 1.0;
        a.num_each_type.push_back(number);
        num_atoms += number;
        for (uint l = 0; l < (uint) number; l++) {
          poccaus.push_back(dpocc);
        }
        for (uint l = 0; l < (uint) number; l++) {
          a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
        }
      }
    }
    // -------------- yes partial occupation
    if (a.partial_occupation_flag == true) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " PARTIAL OCCUPATION = true" << endl;
      }
      for (size_t i = 0; i < tokens_i.size(); i++) {
        if (!aurostd::substring2bool(tokens_i[i], "*") && !aurostd::substring2bool(tokens_i[i], "+")) {
          // NO POCC KEYWORD
          number = aurostd::string2utype<int>(tokens_i[i]);
          dpocc = 1.0;
          a.num_each_type.push_back(number);
          num_atoms += number;
          for (uint l = 0; l < (uint) number; l++) {
            poccaus.push_back(dpocc);
          }
          for (uint l = 0; l < (uint) number; l++) {
            a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
          }
        } else {
          aurostd::string2tokens(tokens_i[i], tokens_j, "+");
          if (LDEBUG) {
            cerr << "tokens_j.size()=" << tokens_j.size() << endl;
          }
          number = 0;
          int nnumber;
          for (size_t j = 0; j < tokens_j.size(); j++) {
            aurostd::string2tokens(tokens_j[j], tokens_k, "*");
            if (tokens_k.empty()) {
              message << "PARTIAL OCCUPATION error [1] tokens_k.size()==0, no *" << endl; // CO20190629
              for (size_t i = 0; i < vinput.size(); i++) {
                message << vinput[i] << endl; // CO20190629
              }
              throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
            }
            nnumber = aurostd::string2utype<int>(tokens_k.at(0));
            if (tokens_k.size() == 1) {
              dpocc = 1.0;
            }
            if (tokens_k.size() == 2) {
              dpocc = aurostd::string2utype<double>(tokens_k.at(1));
            }
            // SOMETHING FOR ROUNDING DPOCC UP TO FRACTIONS....
            number += nnumber;
            num_atoms += nnumber;
            for (uint l = 0; l < (uint) nnumber; l++) {
              poccaus.push_back(dpocc);
            }
            for (uint l = 0; l < (uint) nnumber; l++) {
              if (tokens_k.size() == 1) {
                a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
              }
              // no specie number
              if (tokens_k.size() == 2) {
                a.partial_occupation_sublattice.push_back(i); // put specie number
              }
            }
            if (tokens_k.size() >= 3) {
              message << "PARTIAL OCCUPATION error [1] tokens_k.size()>=3, too many *" << endl; // CO20190629
              for (size_t i = 0; i < vinput.size(); i++) {
                message << vinput[i] << endl; // CO20190629
              }
              throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
            }
          } // loop on +
          a.num_each_type.push_back(number);
        }
      }
      if (LDEBUG) {
        cerr << "P(" << poccaus.size() << ") = ";
        for (size_t j = 0; j < poccaus.size(); j++) {
          cerr << poccaus.at(j) << " ";
        }
        cerr << endl;
      }
    }
    //  cerr << "num_atoms=" << num_atoms << endl; throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_); // num_atoms is the SUM of the atoms in the numbers
    // -------------- COORDINATE TYPES
    // Type of coordinates (Fractional or Cartesian) - only 1st character matters (D/d or C/c).
    // This line might also be the Selective Dynamics line so you must check for that (S/s).
    a.isd = false;
    string stmp;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " DEBUG [1]" << endl;
    }
    if (vinput.size() - 1 < iline) {
      message << "missing line[" << iline << "]" << endl; // CO20190629
      for (size_t i = 0; i < vinput.size(); i++) {
        message << vinput[i] << endl; // CO20190629
      }
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
    } // CO20180420 - check for missing lines
    stmp = vinput.at(iline++);
    aurostd::StringSubstInPlace(stmp, "\t", " ");
    std::vector<string> stmp_tokens;
    aurostd::string2tokens(stmp, stmp_tokens);
    // Note that if there are spaces at the beginning of the line we have to remove them.
    if (stmp_tokens.empty()) {
      message << "Found blank line on line 7. This line should give coordinate type or selective dynamics." << endl;
      // CO20190629
      for (size_t i = 0; i < vinput.size(); i++) {
        message << vinput[i] << endl; // CO20190629
      }
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
    } else {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " DEBUG [2]" << endl;
      }
      string sstmp = stmp_tokens.at(0);
      a.order_parameter_structure = false;
      if (sstmp[0] == 'S' || sstmp[0] == 's') {
        a.isd = true;
        if (vinput.size() - 1 < iline) {
          message << "missing line[" << iline << "]" << endl; // CO20190629
          for (size_t i = 0; i < vinput.size(); i++) {
            message << vinput[i] << endl; // CO20190629
          }
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
        } // CO20180420 - check for missing lines
        stmp = vinput.at(iline++);
        sstmp = aurostd::RemoveSpaces(stmp);
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " DEBUG [3]" << endl;
      }
      if (sstmp[0] == 'D' || sstmp[0] == 'd') {
        //	cerr << "FRAC" << endl;
        a.coord_type[0] = sstmp[0];
        a.coord_flag = _COORDS_FRACTIONAL_;
      } else {
        if (sstmp[0] == 'C' || sstmp[0] == 'c') {
          //	  cerr << "CART" << endl;
          a.coord_type[0] = sstmp[0];
          a.coord_flag = _COORDS_CARTESIAN_;
          if (a.iomode == IOVASP_WYCKCAR) {
            message << "WYCKOFF mode requires FRACTIONAL coordinates (DIRECT)." << endl; // CO20190629
            for (size_t i = 0; i < vinput.size(); i++) {
              message << vinput[i] << endl; // CO20190629
            }
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
          }
        } else {
          message << "Did not find coordinate type D/d or C/c." << endl; // CO20190629
          for (size_t i = 0; i < vinput.size(); i++) {
            message << vinput[i] << endl; // CO20190629
          }
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
        }
      }
      a.coord_type[1] = '\0';
      if (stmp_tokens.size() > 1) {
        for (size_t ipos = 1; ipos < stmp_tokens.size(); ipos++) {
          if (stmp_tokens[ipos][0] == 'O' || stmp_tokens[ipos][0] == 'o') {
            a.order_parameter_structure = true;
          }
          //	  if(stmp_tokens[ipos][0]=='P' || stmp_tokens[ipos][0]=='p') a.partial_occupation_flag=true;
        }
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " DEBUG a.order_parameter_structure=" << a.order_parameter_structure << endl;
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " DEBUG [4]" << endl;
    }

    // --------------  Basis atom positions and names
    //   int cnt;
    xvector<double> v(3);
    clear(v);
    xvector<double> vzero(3);
    clear(vzero);
    // Set default values: names, sd.
    const string def_name = "H";
    const string def_sd = "TTT";
    // NOTHING TO CLEAR     clear(*cpos);clear(*fpos);clear(*types);
    //  Read in all the atom data.
    // The number of fields on a line can vary but must mean something, which depends on isd.
    // For isd=0:
    //   If nf=3 then there are only atomic positions.
    //   If nf>=4 then there are atomic positions and a name.
    // For isd!=0:
    //   If nf=3 then there are only atomic positions.
    //   If nf=4 then there are atomic positions and a name.
    //   If nf=6 then there are atomic positions and sd.
    //   If nf>=7 then there are atomic positions and sd and a name.

    // clear up a little
    // NOTHING TO CLEAR     for(int i=0;i<=MAX_ATOM_UCELL;i++) (*name)[i].clear();          // vector(0,MAX_ATOM_UCELL)
    // NOTHING TO CLEAR     for(int i=0;i<=MAX_ATOM_UCELL;i++) (*sd)[i].clear();             // vector(0,MAX_ATOM_UCELL)
    // NOTHING TO CLEAR     for(int i=1;i<=MAX_ATOM_UCELL;i++) (*name_is_given)(i)=false; // xvector(1,MAX_ATOM_UCELL)

    int iat = -1;
    uint itype = 0;
    uint j;
    const uint iline_ref = iline;
    for (itype = 0; itype < a.num_each_type.size(); itype++) {
      //  cerr << a.num_each_type[itype] << endl;
      for (j = 0; j < (uint) a.num_each_type[itype]; j++) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " DEBUG [5] itype,j=" << itype << "," << j << endl;
        }
        //	bool plug_atom=true;
        iat++; // it startf from -1 so the first is ZERO
        //  cerr << iat << " type=" << itype << endl;
        //  for(int iat=1;iat<=(num_atoms);iat++) { //[CO20200106 - close bracket for indenting]}
        // cerr << iat << " " << (num_atoms) << endl;
        _atom atom; // create new atom
        string stmp;
        stmp = vinput.at(iline++);
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " " << iline << " " << vinput.size() << "," << iline - iline_ref << "," << num_atoms << "," << stmp << endl;
        }
        if (iline == vinput.size() && (iline - iline_ref < (size_t) num_atoms)) {
          message << "Insufficient number of atom lines." << endl; // CO20190629
          for (size_t i = 0; i < vinput.size(); i++) {
            message << vinput[i] << endl; // CO20190629
          }
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
        }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " DEBUG [6]" << endl;
        }
        stmp = aurostd::RemoveCharacter(stmp, '\t');
        aurostd::StringSubstInPlace(stmp, "\t", " ");
        std::vector<string> stmp_tokens;
        aurostd::string2tokens(stmp, stmp_tokens);
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " DEBUG [6b] stmp_tokens.size()=" << stmp_tokens.size() << endl;
        }
        if (stmp_tokens.size() < 3) {
          message << "Insufficient number of atom entries in atom=" << iline - iline_ref << "" << endl; // CO20190629
          for (size_t i = 0; i < vinput.size(); i++) {
            message << vinput[i] << endl; // CO20190629
          }
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
        }

        int id = 0;
        // Read in the atom positions.
        v(1) = atof(stmp_tokens.at(id++).c_str()); // when not a matter of speed, vectors
        v(2) = atof(stmp_tokens.at(id++).c_str()); // should be indicized with () instead of [],
        v(3) = atof(stmp_tokens.at(id++).c_str()); // b/c () does boundary checks !
        if (a.coord_flag == _COORDS_FRACTIONAL_) {
          //    for(int ii=1;ii<=3;ii++)
          atom.fpos = v;
          atom.cpos = F2C(a.lattice, atom.fpos);
        }
        if (a.coord_flag == _COORDS_CARTESIAN_) {
          atom.cpos = v;
          atom.fpos = C2F(a.lattice, atom.cpos);
        }
        //[CO20200130 - number->basis]atom.number=iat;    // reference position for convasp
        atom.basis = iat; // position in the basis
        atom.ijk(1) = 0;
        atom.ijk(2) = 0;
        atom.ijk(3) = 0; // inside the zero cell...
        atom.corigin(1) = 0.0;
        atom.corigin(2) = 0.0;
        atom.corigin(3) = 0.0; // inside the zero cell
        atom.coord(1) = 0.0;
        atom.coord(2) = 0.0;
        atom.coord(3) = 0.0; // inside the zero cell
        atom.spin = 0.0;
        atom.noncoll_spin.clear(); // DX20171205 - non-collinear spin
        atom.type = itype; // CONVASP_MODE if I want type 0,1,2,3,...
        atom.order_parameter_value = 0;
        atom.order_parameter_atom = false;
        atom.partial_occupation_value = 1.0;
        atom.partial_occupation_flag = false;
        //	plug_atom=true;   // not used
        // NO ORDER PARAMETER
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " DEBUG [7]" << endl;
        }
        if (a.order_parameter_structure == false) {
          // stmp_tokens.size() = 4 (plus possible comments).
          // Read in the names.
          if (stmp_tokens.size() == 4 || (stmp_tokens.size() >= 4 && a.isd == 0)) {
            atom.name = stmp_tokens.at(id++);
            atom.CleanName();
            atom.CleanSpin();
            atom.name_is_given = true;
          }
          // stmp_tokens.size()=6
          if (a.isd != 0 && stmp_tokens.size() == 6) {
            string sdt;
            for (int ic = 1; ic <= 3; ic++) {
              sdt = sdt + stmp_tokens.at(id++);
            }
            atom.sd = sdt;
          }
          // stmp_tokens.size()=7
          if (a.isd == true && stmp_tokens.size() >= 7) {
            string sdt;
            for (int ic = 1; ic <= 3; ic++) {
              sdt = sdt + stmp_tokens.at(id++);
            }
            atom.sd = sdt;
            atom.name = stmp_tokens.at(id++);
            atom.CleanName();
            atom.CleanSpin();
            atom.name_is_given = true;
          }
        }
        // ORDER PARAMETER
        if (a.order_parameter_structure == true) {
          if (stmp_tokens.size() != 5 && stmp_tokens.size() != 6) {
            if (a.order_parameter_structure == true) {
              message << "With the order parameter you must specify " << endl; // CO20190629
              message << "  x y z Name OrderParameter " << endl; // CO20190629
              message << " where x,y,z are the coordinates " << endl; // CO20190629
              message << " Name is the symbol of the atom  " << endl; // CO20190629
              message << " Order parameter is -=none, *=consider, -1,0,1 (integer) values " << endl; // CO20190629
              throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
            }
          }
          if (stmp_tokens.size() == 5 || stmp_tokens.size() == 6) {
            atom.name = stmp_tokens.at(id);
            atom.CleanName();
            atom.CleanSpin();
            atom.name_is_given = true;
            id++;
            if (a.order_parameter_structure == true) {
              if (stmp_tokens.at(id) == "-") {
                atom.order_parameter_value = 0;
                atom.order_parameter_atom = false;
              } else {
                if (stmp_tokens.at(id) == "*") {
                  atom.order_parameter_value = 0;
                  atom.order_parameter_atom = true;
                } else {
                  atom.order_parameter_value = aurostd::string2utype<int>(stmp_tokens.at(id).c_str());
                  atom.order_parameter_atom = true;
                  a.order_parameter_sum += atom.order_parameter_value;
                }
              }
              id++;
            }
          }
        }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " DEBUG [8]" << endl;
        }
        // now plug atom into the atomlist
        //      cerr << atom.name << endl;
        //      cerr << atom.cleanname << endl;
        //      cerr << atom.name_is_given << endl;
        a.atoms.push_back(atom);
      } // iat loop
    }
  } //  end VASP POSCAR MODE

  // ----------------------------------------------------------------------
  // QE INPUT
  if (a.iomode == IOQE_AUTO || a.iomode == IOQE_GEOM) {
    // LDEBUG=true;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " QUANTUM ESPRESSO IOQE_AUTO/GEOM" << endl;
    }
    // START FROM CELL
    a.scale = 1.0; // standard
    a.neg_scale = false; // standard
    int ibrav_value = 0; // DX20180123 - added ibrav
    int nat_value = 0; // DX20180123 - added nat
    bool bohr_lat = false; // DX20180215 - added bohr for lattice
    bool bohr_pos = false; // DX20180215 - added bohr for positions
    bool alat = false; // DX20180215 - added alat
    uint iline = 0;
    uint jline = 0;
    iline = vinput.size();
    // DX20180123 - added nat - START
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "NAT=", true)) {
        vector<string> comma_tokens;
        aurostd::string2tokens(aurostd::RemoveWhiteSpaces(vinput[i]), comma_tokens, ",");
        // DX - it is possible to have multiple fields per line
        for (size_t j = 0; j < comma_tokens.size(); j++) {
          if (aurostd::substring2bool(aurostd::toupper(comma_tokens[j]), "NAT=", true)) {
            vector<string> tokens;
            aurostd::string2tokens(aurostd::RemoveWhiteSpaces(comma_tokens[j]), tokens, "=");
            nat_value = aurostd::string2utype<int>(tokens.at(1));
          }
        }
      }
    }
    // DX20180123 - added nat - END
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "CELL_PARAMETERS")) {
        // DX20170124 - added aurostd::toupper
        if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "ANGSTROM")) {
          iline = i + 1;
          jline = i;
        } else if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "ALAT")) {
          // DX20180215 - added alat bool
          iline = i + 1;
          jline = i; // DX20180215 - added alat bool
          alat = true; // DX20180215 - added alat bool
        }
        // DX20180123 - added bohr - START
        else if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "BOHR")) {
          iline = i + 1;
          jline = i;
          bohr_lat = true;
        }
        // DX20180123 - added bohr - END
      }
      // DX20180123 - added ibrav - START
      else if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "IBRAV=", true)) {
        vector<string> comma_tokens;
        aurostd::string2tokens(aurostd::RemoveWhiteSpaces(vinput[i]), comma_tokens, ",");
        // DX - it is possible to have multiple fields per line
        for (size_t j = 0; j < comma_tokens.size(); j++) {
          if (aurostd::substring2bool(aurostd::toupper(comma_tokens[j]), "IBRAV=", true)) {
            vector<string> tokens;
            aurostd::string2tokens(aurostd::RemoveWhiteSpaces(comma_tokens[j]), tokens, "=");
            ibrav_value = aurostd::string2utype<int>(tokens.at(1));
          }
        }
      }
      // DX20180123 - added ibrav - END
    }
    const xvector<double> parameters(6);
    uint celldm_count = 0;
    bool isabc = false; // distinguish between celldm and a,b,c,cosAB,cosAC,cosBC

    // celldm variant
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "CELLDM")) {
        vector<string> comma_tokens;
        aurostd::string2tokens(aurostd::RemoveWhiteSpaces(vinput[i]), comma_tokens, ",");
        // DX - it is possible to have multiple fields per line
        for (size_t j = 0; j < comma_tokens.size(); j++) {
          if (aurostd::substring2bool(aurostd::toupper(comma_tokens[j]), "CELLDM", true)) {
            // DX20170124 - added aurostd::toupper
            celldm_count++;
            vector<string> tokens;
            aurostd::string2tokens(aurostd::RemoveWhiteSpaces(comma_tokens[j]), tokens, "=");
            // In case the celldms are not in order, or not all are given; need to check number
            vector<string> parentheses_tokens;
            aurostd::string2tokens(aurostd::RemoveWhiteSpaces(tokens.at(0)), parentheses_tokens, "(");
            const uint index = aurostd::string2utype<uint>(aurostd::RemoveCharacter(parentheses_tokens.at(1), ')'));
            parameters[index] = aurostd::string2utype<double>(tokens.at(1));
          }
        }
      }
    }
    // a,b,c,cosAB,cosAC,cosBC variant
    if (celldm_count == 0) {
      isabc = true;
      for (size_t i = 0; i < vinput.size(); i++) {
        if (aurostd::substring2bool(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])), "A=")) {
          vector<string> comma_tokens;
          aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])), comma_tokens, ",");
          // DX - it is possible to have multiple fields per line
          for (size_t j = 0; j < comma_tokens.size(); j++) {
            if (aurostd::substring2bool(aurostd::toupper(comma_tokens[j]), "A=", true)) {
              // DX20170124 - added aurostd::toupper
              vector<string> tokens;
              aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(comma_tokens[j])), tokens, "A=");
              if (tokens.size() == 1) {
                parameters[1] = aurostd::string2utype<double>(tokens.at(0)); // A only
              }
            }
          }
        }
        if (aurostd::substring2bool(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])), "B=")) {
          vector<string> comma_tokens;
          aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])), comma_tokens, ",");
          // DX - it is possible to have multiple fields per line
          for (size_t j = 0; j < comma_tokens.size(); j++) {
            if (aurostd::substring2bool(aurostd::toupper(comma_tokens[j]), "B=", true)) {
              // DX20170124 - added aurostd::toupper
              vector<string> tokens;
              aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(comma_tokens[j])), tokens, "B=");
              if (tokens.size() == 1) {
                parameters[2] = aurostd::string2utype<double>(tokens.at(0)); // B only
              }
            }
          }
        }
        if (aurostd::substring2bool(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])), "C=")) {
          vector<string> comma_tokens;
          aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])), comma_tokens, ",");
          // DX - it is possible to have multiple fields per line
          for (size_t j = 0; j < comma_tokens.size(); j++) {
            if (aurostd::substring2bool(aurostd::toupper(comma_tokens[j]), "C=", true)) {
              // DX20170124 - added aurostd::toupper
              vector<string> tokens;
              aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(comma_tokens[j])), tokens, "C=");
              if (tokens.size() == 1) {
                parameters[3] = aurostd::string2utype<double>(tokens.at(0)); // C only
              }
            }
          }
        }
        if (aurostd::substring2bool(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])), "COSAB=")) {
          vector<string> comma_tokens;
          aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])), comma_tokens, ",");
          // DX - it is possible to have multiple fields per line
          for (size_t j = 0; j < comma_tokens.size(); j++) {
            if (aurostd::substring2bool(aurostd::toupper(comma_tokens[j]), "COSAB=", true)) {
              // DX20170124 - added aurostd::toupper
              vector<string> tokens;
              aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(comma_tokens[j])), tokens, "COSAB=");
              if (tokens.size() == 1) {
                parameters[4] = aurostd::string2utype<double>(tokens.at(0)); // COSAB only
              }
            }
          }
        }
        if (aurostd::substring2bool(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])), "COSAC=")) {
          vector<string> comma_tokens;
          aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])), comma_tokens, ",");
          // DX - it is possible to have multiple fields per line
          for (size_t j = 0; j < comma_tokens.size(); j++) {
            if (aurostd::substring2bool(aurostd::toupper(comma_tokens[j]), "COSAC=", true)) {
              // DX20170124 - added aurostd::toupper
              vector<string> tokens;
              aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(comma_tokens[j])), tokens, "COSAC=");
              if (tokens.size() == 1) {
                parameters[5] = aurostd::string2utype<double>(tokens.at(0)); // COSAC only
              }
            }
          }
        }
        if (aurostd::substring2bool(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])), "COSBC=")) {
          vector<string> comma_tokens;
          aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])), comma_tokens, ",");
          // DX - it is possible to have multiple fields per line
          for (size_t j = 0; j < comma_tokens.size(); j++) {
            if (aurostd::substring2bool(aurostd::toupper(comma_tokens[j]), "COSBC=", true)) {
              // DX20170124 - added aurostd::toupper
              vector<string> tokens;
              aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(comma_tokens[j])), tokens, "COSBC=");
              if (tokens.size() == 1) {
                parameters[6] = aurostd::string2utype<double>(tokens.at(0)); // COSBC only
              }
            }
          }
        }
      }
    }
    // DX20180123 - added ibrav/parameters - START
    if (ibrav_value != 0) {
      a.lattice = pflow::QE_ibrav2lattice(ibrav_value, parameters, isabc);
    }
    // DX20180123 - added ibrav/parameters - END
    else {
      if (iline < vinput.size() - 3) {
        // IN a1 / a2 / a3 rows version
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " QUANTUM ESPRESSO FOUND 3 extra lines after, trying a1/a2/a3 on raws" << endl;
        }
        stringstream input_tmp;
        input_tmp.clear();
        input_tmp.str(vinput.at(iline++));
        input_tmp >> a.lattice(1, 1) >> a.lattice(1, 2) >> a.lattice(1, 3);
        input_tmp.clear();
        input_tmp.str(vinput.at(iline++));
        input_tmp >> a.lattice(2, 1) >> a.lattice(2, 2) >> a.lattice(2, 3);
        input_tmp.clear();
        input_tmp.str(vinput.at(iline++));
        input_tmp >> a.lattice(3, 1) >> a.lattice(3, 2) >> a.lattice(3, 3);
        if (bohr_lat) {
          // DX20180123 - added bohr
          a.lattice = a.lattice * bohr2angstrom; // DX20180123 - added bohr
        } // DX20180123 - added bohr
        if (alat) {
          // DX20180215 - added alat
          if (isabc) {
            a.lattice = a.lattice * parameters[1]; // DX20180123 - added alat
          } else {
            a.lattice = a.lattice * parameters[1] * bohr2angstrom; // DX20180123 - added alat
          }
        } // DX20180215 - added alat
        xvector<double> data(6);
        data = Getabc_angles(a.lattice, DEGREES);
        a.a = data[1];
        a.b = data[2];
        a.c = data[3];
        a.alpha = data[4];
        a.beta = data[5];
        a.gamma = data[6];
      }
    }
    // ---------------------------------------------------------------
    // SOME FOR ABCCAR
    // ...
    // ---------------------------------------------------------------
    // SOME FOR WYCCAR
    // ...
    a.FixLattices(); // Reciprocal/f2c/c2f
    a.kpoints_k1 = 0;
    a.kpoints_k2 = 0;
    a.kpoints_k3 = 0;
    a.kpoints_kmax = 0;
    a.kpoints_kppra = 0;
    a.kpoints_kscheme = "";
    clear(a.origin);

    // cerr << "iline=" << iline << endl;
    // cerr << "vinput.size()=" << vinput.size() << endl;
    // NOW ADD ATOMS
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "ATOMIC_POSITIONS") && aurostd::substring2bool(aurostd::toupper(vinput[i]), "CRYSTAL")) // DX20180124 - added aurostd::toupper
      {
        // CO20200106 - patching for auto-indenting
        iline = i + 1;
        a.coord_flag = _COORDS_FRACTIONAL_;
        jline = iline + nat_value; // DX20180123 - added nat
      }
      if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "ATOMIC_POSITIONS") && aurostd::substring2bool(aurostd::toupper(vinput[i]), "ANGSTROM")) // DX20180124 - added aurostd::toupper
      {
        // CO20200106 - patching for auto-indenting
        iline = i + 1;
        a.coord_flag = _COORDS_CARTESIAN_;
        jline = iline + nat_value; // DX20180123 - added nat
      }
      // DX20180124 -- added Bohr case - START
      if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "ATOMIC_POSITIONS") && aurostd::substring2bool(aurostd::toupper(vinput[i]), "BOHR")) {
        // DX20180124 - added aurostd::toupper
        iline = i + 1;
        a.coord_flag = _COORDS_CARTESIAN_;
        jline = iline + nat_value; // DX20180123 - added nat
        bohr_pos = true;
      }
      // DX20180124 -- added Bohr case - END
    }

    for (size_t i = iline; i < vinput.size() && i < jline; i++) {
      _atom atom; // create new atom
      string stmp;
      xvector<double> v(3);
      clear(v);
      stmp = vinput[i];
      stmp = aurostd::RemoveCharacter(vinput[i], '\t');
      aurostd::StringSubstInPlace(stmp, "\t", " ");
      std::vector<string> stmp_tokens;
      aurostd::string2tokens(stmp, stmp_tokens);
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " DEBUG [6b] stmp_tokens.size()=" << stmp_tokens.size() << endl;
      }
      if (stmp_tokens.size() < 4) {
        message << "Insufficient number of atom entries in atom=" << i << "" << endl; // CO20190629
        for (size_t i = 0; i < vinput.size(); i++) {
          message << vinput[i] << endl; // CO20190629
        }
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
      }
      int id = 0;
      string name;
      // Read in the atom positions.
      name = stmp_tokens.at(id++);
      v(1) = aurostd::string2utype<double>(stmp_tokens.at(id++)); // when not a matter of speed, vectors
      v(2) = aurostd::string2utype<double>(stmp_tokens.at(id++)); // should be indicized with () instead of [],
      v(3) = aurostd::string2utype<double>(stmp_tokens.at(id++)); // b/c () does boundary checks !
      if (a.coord_flag == _COORDS_FRACTIONAL_) {
        atom.fpos = v;
        atom.cpos = F2C(a.lattice, atom.fpos);
      }
      if (a.coord_flag == _COORDS_CARTESIAN_) {
        if (bohr_pos) {
          // DX20180124 - added bohr case
          atom.cpos = bohr2angstrom * v; // DX20180124 - added bohr case
          atom.fpos = C2F(a.lattice, atom.cpos);
        } else {
          // DX20180124 - added bohr case
          atom.cpos = v;
          atom.fpos = C2F(a.lattice, atom.cpos);
        }
      }
      atom.name = name;
      atom.CleanName();
      atom.CleanSpin();
      atom.name_is_given = true;

      //[CO20200130 - number->basis]atom.number=a.atoms.size();    // reference position for convasp
      atom.basis = a.atoms.size(); // position in the basis
      atom.ijk(1) = 0;
      atom.ijk(2) = 0;
      atom.ijk(3) = 0; // inside the zero cell...
      atom.corigin(1) = 0.0;
      atom.corigin(2) = 0.0;
      atom.corigin(3) = 0.0; // inside the zero cell
      atom.coord(1) = 0.0;
      atom.coord(2) = 0.0;
      atom.coord(3) = 0.0; // inside the zero cell
      atom.spin = 0.0;
      atom.noncoll_spin.clear(); // DX20171205 - non-collinear spin
      // FIXED BELOW atom.type=itype;                // CONVASP_MODE if I want type 0,1,2,3,...
      atom.order_parameter_value = 0;
      atom.order_parameter_atom = false;
      atom.partial_occupation_value = 1.0;
      atom.partial_occupation_flag = false;
      // DONE
      a.AddAtom(atom, true); // CO20230319 - add by species
      // NO PARTIAL OCCUPATION
      a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
    }
    // FIX TITLE
    // FIX ITYPE
    // for(size_t i=0;i<a.atoms.size();i++) cerr << "a.atoms.at(i).type=" << a.atoms.at(i).type << endl;
    a.SpeciesPutAlphabetic(); // fight analphabetization
    uint iat = 0;
    string _title; // CO20200731 - NO_TITLE_GIVEN clear
    for (size_t itype = 0; itype < a.num_each_type.size(); itype++) {
      for (uint j = 0; j < (uint) a.num_each_type[itype]; j++) {
        if (j == 0) {
          _title += a.atoms.at(iat).name + aurostd::utype2string(a.num_each_type[itype]);
        }
        a.atoms.at(iat++).type = itype;
      }
    }
    if (!_title.empty()) {
      a.title = _title;
    } // CO20200731 - safer for empty titles
    // for(size_t i=0;i<a.atoms.size();i++) cerr << "a.atoms.at(i).type=" << a.atoms.at(i).type << endl;
    a.partial_occupation_flag = false;
    //    a.xstructure2vasp();
    a.is_vasp4_poscar_format = false;
    a.is_vasp5_poscar_format = false;
    // DONE ?
  } // QE INPUT

  // ----------------------------------------------------------------------
  // ABINIT INPUT (DX20200310)
  if (a.iomode == IOABINIT_AUTO || a.iomode == IOABINIT_GEOM) {
    // ABINIT
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ABINIT READER begin" << endl;
    }

    a.scale = 1.0; // standard
    a.neg_scale = false; // standard

    // ----------------------------------------------------------------------
    // get lattice scaling: acell (optional keyword)
    xvector<double> acell;
    acell(1) = 1.0;
    acell(2) = 1.0;
    acell(3) = 1.0;
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "ACELL", true)) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " ABINIT READER acell line=" << vinput[i] << endl;
        }
        string acell_line = aurostd::toupper(vinput[i]);
        aurostd::StringSubstInPlace(acell_line, "ACELL", "");

        // get units first (then remove from temp string for easy parsing)
        string lattice_vec_unit = "Bohr"; // default is Bohr
        if (aurostd::substring2bool(acell_line, "ANGST", true)) {
          lattice_vec_unit = "Angstrom";
          aurostd::StringSubstInPlace(acell_line, "ANGSTROM", "");
          aurostd::StringSubstInPlace(acell_line, "ANGSTR", "");
          aurostd::StringSubstInPlace(acell_line, "ANGST", "");
        }
        if (aurostd::substring2bool(acell_line, "BOHR", true)) {
          lattice_vec_unit = "Bohr";
          aurostd::StringSubstInPlace(acell_line, "BOHR", "");
        }

        // clean-up
        acell_line = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(acell_line);

        // check if explicitly given or uses multiplication
        // explicit : e.g., "acell 1.0 1.0 1.0"
        // multiplication : e.g., "acell 3*1.0"
        bool multiplication_variant = false;
        if (aurostd::substring2bool(acell_line, "*", true)) {
          multiplication_variant = true;
        }

        if (multiplication_variant) {
          vector<string> tokens;
          aurostd::string2tokens(acell_line, tokens, "*");
          const double factor = aurostd::string2utype<double>(aurostd::RemoveWhiteSpaces(tokens[1]));
          acell = factor * acell;
        } else {
          vector<string> tokens;
          const uint number_tokens = aurostd::string2tokens(acell_line, tokens, " ");
          if (number_tokens == 3) {
            acell(1) = aurostd::string2utype<double>(tokens[0]);
            acell(2) = aurostd::string2utype<double>(tokens[1]);
            acell(3) = aurostd::string2utype<double>(tokens[2]);
          } else {
            message << "Unable to parse the acell line; unexpected format. acell_line = " << acell_line << " tokens: " << aurostd::joinWDelimiter(tokens, ",");
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
          }
        }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " ABINIT READER extracted acell = " << acell << endl;
        }
      }
    }

    // ----------------------------------------------------------------------
    // get lattice (rprim) : given column-major (optional)
    // default: identity (i.e., cube)
    // supported formats:
    //   1) single line 9 fields
    //   2) three lines 3 fields per line
    bool is_lattice_line = false;
    uint lattice_line_count = 0;
    a.lattice = aurostd::identity((double) 0, 3); // abinit default //DX20200521 - new identity format
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "RPRIM", true)) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " ABINIT READER rprim line found = " << vinput[i] << endl;
        }
        is_lattice_line = true;
      }
      if (is_lattice_line) {
        string rprim_line = aurostd::toupper(vinput[i]);
        aurostd::StringSubstInPlace(rprim_line, "RPRIM", "");
        vector<string> tokens;
        const uint number_tokens = aurostd::string2tokens(rprim_line, tokens, " ");
        if (number_tokens == 9) {
          // column-major
          a.lattice(1, 1) = aurostd::frac2dbl(tokens[0]) * acell(1);
          a.lattice(2, 1) = aurostd::frac2dbl(tokens[1]) * acell(1);
          a.lattice(3, 1) = aurostd::frac2dbl(tokens[2]) * acell(1);
          a.lattice(1, 2) = aurostd::frac2dbl(tokens[3]) * acell(2);
          a.lattice(2, 2) = aurostd::frac2dbl(tokens[4]) * acell(2);
          a.lattice(3, 2) = aurostd::frac2dbl(tokens[5]) * acell(2);
          a.lattice(1, 3) = aurostd::frac2dbl(tokens[6]) * acell(3);
          a.lattice(2, 3) = aurostd::frac2dbl(tokens[7]) * acell(3);
          a.lattice(3, 3) = aurostd::frac2dbl(tokens[8]) * acell(3);
          is_lattice_line = false;
          break;
        } else if (number_tokens == 3) {
          lattice_line_count += 1;
          // column-major
          a.lattice(1, lattice_line_count) = aurostd::frac2dbl(tokens[0]) * acell(lattice_line_count);
          a.lattice(2, lattice_line_count) = aurostd::frac2dbl(tokens[1]) * acell(lattice_line_count);
          a.lattice(3, lattice_line_count) = aurostd::frac2dbl(tokens[2]) * acell(lattice_line_count);
          if (lattice_line_count == 3) {
            is_lattice_line = false;
            break;
          }
        } else if (number_tokens == 0) {
          continue;
        } else {
          message << "Unable to parse the rprim line; unexpected format.";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
        }
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ABINIT READER lattice (row-major; AFLOW convention) = " << endl << a.lattice << endl;
    }

    // ----------------------------------------------------------------------
    // get number of atoms
    // NOTE: sometimes multiple keywords can be in the same line
    // natom : number of atoms in unit cell
    // natrd : number of atoms to read in file (useful for symmetry reduced files, e.g., iatom representations)
    uint number_of_atoms = 0;
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "NATOM", true)) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " ABINIT READER natom line = " << vinput[i] << endl;
        }
        const string natom_line = aurostd::toupper(vinput[i]);
        vector<string> tokens;
        aurostd::string2tokens(natom_line, tokens, " ");
        // the loop method protects against multiple keywords per line
        for (size_t t = 0; t < tokens.size(); t++) {
          if (tokens[t] == "NATOM") {
            number_of_atoms = aurostd::string2utype<uint>(tokens[t + 1]);
            break;
          }
        }
        if (number_of_atoms == 0) {
          message << "Unable to parse the natom line; unexpected format. vinput[i] = \"" << natom_line << "\".";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
        }
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ABINIT READER number of atoms; natom = " << number_of_atoms << endl;
    }

    // natrd (optional)
    bool found_natrd = false;
    uint number_of_atoms_to_read = 0;
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "NATRD", true)) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " ABINIT READER natrd line = " << vinput[i] << endl;
        }
        found_natrd = true;
        const string natrd_line = aurostd::toupper(vinput[i]);
        vector<string> tokens;
        aurostd::string2tokens(natrd_line, tokens, " ");
        // the loop method protects against multiple keywords per line
        for (size_t t = 0; t < tokens.size(); t++) {
          if (tokens[t] == "NATRD") {
            number_of_atoms_to_read = aurostd::string2utype<uint>(tokens[t + 1]);
            break;
          }
        }
        if (number_of_atoms_to_read == 0) {
          message << "Unable to parse the natrd line; unexpected format. natrd_line = \"" << natrd_line << "\".";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
        }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " ABINIT READER number of atoms to read; natrd = " << number_of_atoms_to_read << endl;
        }
      }
    }
    if (!found_natrd) {
      number_of_atoms_to_read = number_of_atoms;
    } // if natrd not specified set to natom

    if (number_of_atoms != number_of_atoms_to_read) {
      message << "The natoms != natrd, i.e., unit cell requires atoms (iatoms) to be expanded by symmetry. Functionality not yet supported.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
    }

    // ----------------------------------------------------------------------
    // get number of types (optional)
    bool found_ntypat = false;
    uint number_of_atom_types = 0;
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "NTYPAT", true)) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " ABINIT READER ntypat line = " << vinput[i] << endl;
        }
        found_ntypat = true;
        const string ntypat_line = aurostd::toupper(vinput[i]);
        vector<string> tokens;
        aurostd::string2tokens(ntypat_line, tokens, " ");
        // the loop method protects against multiple keywords per line
        for (size_t t = 0; t < tokens.size(); t++) {
          if (tokens[t] == "NTYPAT") {
            number_of_atom_types = aurostd::string2utype<uint>(tokens[t + 1]);
            break;
          }
        }
        if (number_of_atoms_to_read == 0) {
          message << "Unable to parse the ntypat line; unexpected format. ntypat_line = \"" << ntypat_line << "\".";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
        }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " ABINIT READER number of atom types; ntypat = " << number_of_atom_types << endl;
        }
      }
    }

    // ----------------------------------------------------------------------
    // get atom type index
    // in same order as xred, xcart, or xangst
    vector<uint> atom_types; // follows order of atoms
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "TYPAT", true)) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " ABINIT READER typat line = " << vinput[i] << endl;
        }
        const string typat_line = aurostd::toupper(vinput[i]);
        vector<string> tokens;
        aurostd::string2tokens(typat_line, tokens, " ");
        // the loop method protects against multiple keywords per line
        bool found_typat = false;
        bool is_field_attribute = false;
        for (size_t t = 0; t < tokens.size(); t++) {
          if (tokens[t] == "TYPAT") {
            found_typat = true;
            is_field_attribute = true;
            continue;
          }
          if (is_field_attribute) {
            if (aurostd::substring2bool(tokens[t], "*")) {
              vector<string> sub_tokens;
              aurostd::string2tokens(tokens[t], sub_tokens, "*");
              const uint multiplier = aurostd::string2utype<uint>(sub_tokens[0]);
              const uint type_index = aurostd::string2utype<uint>(sub_tokens[1]);
              for (uint m = 0; m < multiplier; m++) {
                atom_types.push_back(type_index);
              }
            } else if (tokens[t][0] >= '0' && tokens[t][0] <= '9') {
              // is digit
              atom_types.push_back(aurostd::string2utype<uint>(tokens[t]));
            } else {
              // signals a new keyword in the same line
              is_field_attribute = false;
              break;
            }
          }
        }
        if (!found_typat) {
          continue;
        } // found ntypat not typat, try another line
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ABINIT READER atom type sequence (matches order of atom positions); typat = " << aurostd::joinWDelimiter(atom_types, ",") << endl;
    }

    // check the number of unique types matches ntypat
    vector<uint> unique_types = atom_types;
    std::stable_sort(unique_types.begin(), unique_types.end());
    unique_types.erase(std::unique(unique_types.begin(), unique_types.end()), unique_types.end());
    if (found_ntypat && (unique_types.size() != number_of_atom_types)) {
      message << "Number of atom types does not match ntypat variable. typat=" << aurostd::joinWDelimiter(unique_types, ",") << " | ntypat=" << number_of_atom_types;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    if (atom_types.size() != number_of_atoms_to_read) {
      message << "Number of atom types does not match the number of atoms to read. typat=" << aurostd::joinWDelimiter(atom_types, ",") << " | natom=" << number_of_atoms;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    // ----------------------------------------------------------------------
    // get atom type; element (znucl)
    vector<uint> nuclear_charge;
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "ZNUCL", true)) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " ABINIT READER znucl line = " << vinput[i] << endl;
        }
        const string znucl_line = aurostd::toupper(vinput[i]);
        vector<string> tokens;
        aurostd::string2tokens(znucl_line, tokens, " ");
        // the loop method protects against multiple keywords per line
        bool is_field_attribute = false;
        for (size_t t = 0; t < tokens.size(); t++) {
          if (tokens[t] == "ZNUCL") {
            is_field_attribute = true;
            continue;
          }
          if (is_field_attribute) {
            if (tokens[t][0] >= '0' && tokens[t][0] <= '9') {
              // is digit
              nuclear_charge.push_back(aurostd::string2utype<uint>(tokens[t]));
            } else {
              // signals a new keyword in the same line
              is_field_attribute = false;
              break;
            }
          }
        }
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ABINIT READER nuclear charge of each atom type (order corresponds to typat index); znucl = " << aurostd::joinWDelimiter(nuclear_charge, ",") << endl;
    }

    // ----------------------------------------------------------------------
    // normally a mandatory keyword, but the AFLOW-ABINIT writer has not been printing this keyword
    // for backwards compatability, we will throw a warning (for now)
    if (nuclear_charge.empty()) {
      message << "The atom elements/types (znucl) are not specified. Using fictious atoms.";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, std::cerr, _LOGGER_WARNING_);
    }

    // ----------------------------------------------------------------------
    // get atom positions
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ABINIT READER reading atom positions" << endl;
    }
    deque<_atom> atoms_temp;
    bool is_atom_line = false;
    bool is_Bohr_units = false;
    uint atom_line_count = 0;
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "XRED", true)) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " ABINIT READER xred (fractional) line = " << vinput[i] << endl;
        }
        a.coord_flag = _COORDS_FRACTIONAL_;
        is_atom_line = true;
      } else if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "XCART", true)) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " ABINIT READER xcart (Cartesian, unit=Bohr) line = " << vinput[i] << endl;
        }
        a.coord_flag = _COORDS_CARTESIAN_;
        is_Bohr_units = true;
        is_atom_line = true;
      } else if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "XANGST", true)) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " ABINIT READER xangst (Cartesian, unit=Angstrom) line = " << vinput[i] << endl;
        }
        a.coord_flag = _COORDS_CARTESIAN_;
        is_atom_line = true;
      }
      if (is_atom_line) {
        string atom_line = aurostd::toupper(vinput[i]);
        aurostd::StringSubstInPlace(atom_line, "XRED", ""); // for easier parsing
        aurostd::StringSubstInPlace(atom_line, "XCART", ""); // for easier parsing
        aurostd::StringSubstInPlace(atom_line, "XANGST", ""); // for easier parsing
        vector<string> tokens;
        const uint number_tokens = aurostd::string2tokens(atom_line, tokens, " ");
        if (number_tokens == 3) {
          atom_line_count += 1;
          const xvector<double> coordinate;
          coordinate(1) = aurostd::frac2dbl(tokens[0]);
          coordinate(2) = aurostd::frac2dbl(tokens[1]);
          coordinate(3) = aurostd::frac2dbl(tokens[2]);

          _atom atom;
          if (a.coord_flag == _COORDS_FRACTIONAL_) {
            atom.fpos = coordinate;
            atom.cpos = F2C(a.lattice, atom.fpos);
          } else if (a.coord_flag == _COORDS_CARTESIAN_) {
            atom.cpos = coordinate;
            if (is_Bohr_units) {
              atom.cpos *= bohr2angstrom;
            } // AFLOW expects Angstroms
            atom.fpos = C2F(a.lattice, atom.cpos);
          }
          atoms_temp.push_back(atom);
        } else if (number_tokens == (3 * number_of_atoms_to_read)) {
          for (size_t t = 0; t < tokens.size(); t += 3) {
            atom_line_count += 1;
            const xvector<double> coordinate;
            coordinate(1) = aurostd::frac2dbl(tokens[t]);
            coordinate(2) = aurostd::frac2dbl(tokens[t + 1]);
            coordinate(3) = aurostd::frac2dbl(tokens[t + 2]);

            _atom atom;
            if (a.coord_flag == _COORDS_FRACTIONAL_) {
              atom.fpos = coordinate;
              atom.cpos = F2C(a.lattice, atom.fpos);
            } else if (a.coord_flag == _COORDS_CARTESIAN_) {
              atom.cpos = coordinate;
              if (is_Bohr_units) {
                atom.cpos *= bohr2angstrom;
              } // AFLOW expects Angstroms
              atom.fpos = C2F(a.lattice, atom.cpos);
            }
            atoms_temp.push_back(atom);
          }
        } else if (number_tokens == 0) {
          continue;
        } else {
          message << "Unable to parse the atom line; unexpected format. vinput[i] = " << atom_line;
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
        }
        if (atom_line_count == number_of_atoms_to_read) {
          break;
        }
      }
    }

    if (LDEBUG) {
      for (size_t iat = 0; iat < atoms_temp.size(); iat++) {
        cerr << __AFLOW_FUNC__ << " ABINIT READER atom position [" << iat << "] = " << atoms_temp[iat] << endl;
      }
    }

    // ----------------------------------------------------------------------
    // ensure correct number of atoms
    if (atoms_temp.size() != number_of_atoms_to_read) {
      message << "The number of atoms is not commensurate: number of atom positions found = " << atoms_temp.size() << " vs natom/natrd = " << number_of_atoms_to_read;
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, std::cerr, _RUNTIME_ERROR_);
    }

    // ----------------------------------------------------------------------
    // add atoms to xstructure
    // add name/type info as well
    for (size_t i = 0; i < atoms_temp.size(); i++) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " ABINIT READER finding element/type for atom [" << a << "]" << endl;
      }
      xelement::xelement element;
      // use nuclear charge associated atom index
      if (!nuclear_charge.empty()) {
        element = xelement::xelement(nuclear_charge[atom_types[i] - 1]);
      }
      // i-1 since index doesn't start at zero
      //  use arbitrary element names (in order of increasing Z)
      else {
        element = xelement::xelement(atom_types[i]);
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " ABINIT READER element/type for atom [" << a << "] = " << element.symbol << endl;
      }
      atoms_temp[i].name = element.symbol;
      atoms_temp[i].CleanName();
      atoms_temp[i].CleanSpin();
      atoms_temp[i].name_is_given = true;

      //[CO20200130 - number->basis]atoms_temp[i].number=atoms_temp.size();    // reference position for convasp
      atoms_temp[i].basis = atoms_temp.size(); // position in the basis
      atoms_temp[i].ijk(1) = 0;
      atoms_temp[i].ijk(2) = 0;
      atoms_temp[i].ijk(3) = 0; // inside the zero cell...
      atoms_temp[i].corigin(1) = 0.0;
      atoms_temp[i].corigin(2) = 0.0;
      atoms_temp[i].corigin(3) = 0.0; // inside the zero cell
      atoms_temp[i].coord(1) = 0.0;
      atoms_temp[i].coord(2) = 0.0;
      atoms_temp[i].coord(3) = 0.0; // inside the zero cell
      atoms_temp[i].spin = 0.0;
      atoms_temp[i].noncoll_spin.clear(); // DX20171205 - non-collinear spin
      // FIXED BELOW atom.type=itype;                // CONVASP_MODE if I want type 0,1,2,3,...
      atoms_temp[i].order_parameter_value = 0;
      atoms_temp[i].order_parameter_atom = false;
      atoms_temp[i].partial_occupation_value = 1.0;
      atoms_temp[i].partial_occupation_flag = false;
      // DONE
      a.AddAtom(atoms_temp[i], true); // CO20230319 - add by species
      // NO PARTIAL OCCUPATION
      a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ABINIT READER fixing atom information (alphabetize, make basis, set species, etc.)" << endl;
    }
    a.SpeciesPutAlphabetic();
    std::stable_sort(a.atoms.begin(), a.atoms.end(), sortAtomsNames);
    a.MakeBasis();
    a.MakeTypes(); // DX20190508 - otherwise types are not created
    // add system name to title
    a.buildGenericTitle(); // DX20210211
    a.partial_occupation_flag = false;
    a.is_vasp4_poscar_format = false;
    a.is_vasp5_poscar_format = false;

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ABINIT READER - Finished" << endl;
    }
  }

  // ----------------------------------------------------------------------
  // ELK INPUT - START (DX20200310)
  // based on information from http://elk.sourceforge.net/elk.pdf
  if (a.iomode == IOELK_AUTO || a.iomode == IOELK_GEOM) {
    // ELK
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ELK READER begin" << endl;
    }

    a.scale = 1.0; // standard
    a.neg_scale = false; // standard
    a.coord_flag = _COORDS_FRACTIONAL_; // always for ELK

    // ----------------------------------------------------------------------
    // get isotropic lattice scaling; two types
    //   1) scale : isotropic scaling
    //   2) scale1/2/3 : anisotropic scaling
    double isotropic_scaling = 1.0;
    const xvector<double> anisotropic_scaling;
    anisotropic_scaling(1) = 1.0;
    anisotropic_scaling(2) = 1.0;
    anisotropic_scaling(3) = 1.0;

    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "SCALE", true)) {
        const string scaling_title_line = aurostd::toupper(vinput[i]);
        // ----------------------------------------------------------------------
        // anisotropic
        if (aurostd::substring2bool(scaling_title_line, "SCALE1", true) || aurostd::substring2bool(scaling_title_line, "SCALE2", true) || aurostd::substring2bool(scaling_title_line, "SCALE3", true)) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " ELK READER isotropic scaling line found = " << vinput[i] << endl;
          }
          const uint vector_index = aurostd::string2utype<uint>(aurostd::RemoveWhiteSpacesFromTheFrontAndBack(aurostd::StringSubst(scaling_title_line, "SCALE", "")));
          const string scaling_line = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vinput[i + 1]);
          // scaling value is in the next line
          anisotropic_scaling(vector_index) = aurostd::string2utype<double>(scaling_line);

          if (anisotropic_scaling(vector_index) < _ZERO_TOL_) {
            message << "Anisotropic scaling factor is zero or negative for " << vector_index << " component;; unable to parse the scale line; unexpected format. scaling_line = " << scaling_line;
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
          }
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " ELK READER anisotropic scaling for vector " << vector_index << " = " << anisotropic_scaling(vector_index) << endl;
          }
        }
        // ----------------------------------------------------------------------
        // isotropic
        else {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " ELK READER isotropic scaling line found = " << vinput[i] << endl;
          }
          const string scaling_line = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vinput[i + 1]);
          // scaling is the next line
          isotropic_scaling = aurostd::string2utype<double>(scaling_line);
          if (isotropic_scaling < _ZERO_TOL_) {
            message << "Isotropic scaling factor is zero or negative; unable to parse the scale line; unexpected format. scaling_line = " << scaling_line;
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
          }
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " ELK READER isotropic scaling = " << isotropic_scaling << endl;
          }
        }
      }
    }

    // ----------------------------------------------------------------------
    // get lattice (avec) : row-major
    bool is_lattice_line = false;
    uint lattice_line_count = 0;
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "AVEC", true)) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " ELK READER avec line found = " << vinput[i] << endl;
        }
        is_lattice_line = true;
      }
      if (is_lattice_line) {
        string avec_line = aurostd::toupper(vinput[i]);
        aurostd::StringSubstInPlace(avec_line, "AVEC", "");
        vector<string> tokens;
        const uint number_tokens = aurostd::string2tokens(avec_line, tokens, " ");
        if (number_tokens == 9) {
          // row-major
          a.lattice(1, 1) = aurostd::frac2dbl(tokens[0]) * isotropic_scaling * anisotropic_scaling(1);
          a.lattice(1, 2) = aurostd::frac2dbl(tokens[1]) * isotropic_scaling * anisotropic_scaling(1);
          a.lattice(1, 3) = aurostd::frac2dbl(tokens[2]) * isotropic_scaling * anisotropic_scaling(1);
          a.lattice(2, 1) = aurostd::frac2dbl(tokens[3]) * isotropic_scaling * anisotropic_scaling(2);
          a.lattice(2, 2) = aurostd::frac2dbl(tokens[4]) * isotropic_scaling * anisotropic_scaling(2);
          a.lattice(2, 3) = aurostd::frac2dbl(tokens[5]) * isotropic_scaling * anisotropic_scaling(2);
          a.lattice(3, 1) = aurostd::frac2dbl(tokens[6]) * isotropic_scaling * anisotropic_scaling(3);
          a.lattice(3, 2) = aurostd::frac2dbl(tokens[7]) * isotropic_scaling * anisotropic_scaling(3);
          a.lattice(3, 3) = aurostd::frac2dbl(tokens[8]) * isotropic_scaling * anisotropic_scaling(3);
          is_lattice_line = false;
          break;
        } else if (number_tokens == 3) {
          lattice_line_count += 1;
          // row-major
          a.lattice(lattice_line_count, 1) = aurostd::frac2dbl(tokens[0]) * isotropic_scaling * anisotropic_scaling(lattice_line_count);
          a.lattice(lattice_line_count, 2) = aurostd::frac2dbl(tokens[1]) * isotropic_scaling * anisotropic_scaling(lattice_line_count);
          a.lattice(lattice_line_count, 3) = aurostd::frac2dbl(tokens[2]) * isotropic_scaling * anisotropic_scaling(lattice_line_count);
          if (lattice_line_count == 3) {
            is_lattice_line = false;
            break;
          }
        } else if (number_tokens == 0) {
          continue;
        } else {
          message << "Unable to parse the avec line; unexpected format.";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
        }
      }
    }
    a.lattice = a.lattice * bohr2angstrom; // convert to Angstroms
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ELK READER lattice (row-major; AFLOW convention) = " << endl << a.lattice << endl;
    }

    // ----------------------------------------------------------------------
    // get number of atom types (atoms)
    double number_of_atom_types = 0;

    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "ATOMS", true) && !aurostd::substring2bool(aurostd::toupper(vinput[i]), "NATOMS", true)) {
        vector<string> tokens;
        aurostd::string2tokens(vinput[i + 1], tokens, " "); // scaling value is in the next line
        number_of_atom_types = aurostd::string2utype<uint>(aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[0]));

        if (number_of_atom_types == 0) {
          message << "Unable to parse the atoms line; unexpected format. atoms_line = \"" << tokens[0] << "\".";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
        }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " ELK READER number of atom types atoms = " << number_of_atom_types << endl;
        }
        break;
      }
    }

    // ----------------------------------------------------------------------
    // get atoms (elements, positions, etc.)
    uint number_of_element_dot_in_files = 0;
    vector<uint> num_each_type;

    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(aurostd::toupper(vinput[i]), ".IN", true)) {
        number_of_element_dot_in_files++;
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " ELK READER \"ELEMENT.in\" file line found = " << vinput[i] << endl;
        }

        vector<string> tokens;
        // ----------------------------------------------------------------------
        // get element
        aurostd::string2tokens(vinput[i], tokens, ".");
        string element_symbol = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[0]);
        // clean //DX20210409 - remove whitespace before as well
        element_symbol = aurostd::RemoveCharacterFromTheFrontAndBack(element_symbol, '\''); // clean
        element_symbol = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(element_symbol); // clean
        tokens.clear();
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " ELK READER element extracted = " << element_symbol << endl;
        }

        // ----------------------------------------------------------------------
        // get number of that element (+1 line down)
        aurostd::string2tokens(vinput[i + 1], tokens, " "); // note: +1 line down
        const uint num_of_this_type = aurostd::string2utype<uint>(aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[0]));
        // clean
        num_each_type.push_back(num_of_this_type);
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " ELK READER number of " << element_symbol << " atoms = " << num_of_this_type << endl;
        }
        tokens.clear();

        // ----------------------------------------------------------------------
        // get subsequent atom positions
        for (uint iat = 1; iat <= num_of_this_type; iat++) {
          const uint line_index = i + 1 + iat; // keep track of line index
          _atom atom;

          const string atom_line = vinput[line_index];
          const uint number_of_tokens = aurostd::string2tokens(atom_line, tokens, " ");

          // first three fields are atom positions (fractional)
          if (number_of_tokens >= 3) {
            atom.fpos(1) = aurostd::string2utype<double>(tokens[0]);
            atom.fpos(2) = aurostd::string2utype<double>(tokens[1]);
            atom.fpos(3) = aurostd::string2utype<double>(tokens[2]);

            if (number_of_tokens == 6) {
              // subsequent are external magnetic field in Cartesian coords
              // not stored beyond here (for now)
              const xvector<double> magnetic_field;
              magnetic_field(1) = aurostd::string2utype<double>(tokens[3]);
              magnetic_field(2) = aurostd::string2utype<double>(tokens[4]);
              magnetic_field(3) = aurostd::string2utype<double>(tokens[5]);
              atom.noncoll_spin_is_given = true; // DX20210409
              atom.noncoll_spin = magnetic_field; // DX20210409
              if (LDEBUG) {
                cerr << __AFLOW_FUNC__ << " ELK READER magnetic field/spin found : " << atom.noncoll_spin << endl;
              }
            }
          }
          // F2C
          atom.cpos = F2C(a.lattice, atom.fpos);

          atom.name = element_symbol;
          atom.CleanName();
          atom.name_is_given = true;

          //[CO20200130 - number->basis]atom.number=a.atoms.size();    // reference position for convasp
          atom.basis = a.atoms.size(); // position in the basis
          atom.ijk(1) = 0;
          atom.ijk(2) = 0;
          atom.ijk(3) = 0; // inside the zero cell...
          atom.corigin(1) = 0.0;
          atom.corigin(2) = 0.0;
          atom.corigin(3) = 0.0; // inside the zero cell
          atom.coord(1) = 0.0;
          atom.coord(2) = 0.0;
          atom.coord(3) = 0.0; // inside the zero cell
          atom.order_parameter_value = 0;
          atom.order_parameter_atom = false;
          atom.partial_occupation_value = 1.0;
          atom.partial_occupation_flag = false;

          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " ELK READER atom added = " << atom << endl;
          }
          // DONE
          a.AddAtom(atom, true); // CO20230319 - add by species
          // NO PARTIAL OCCUPATION
          a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
        }
        i += num_of_this_type + 1; // no need to read over the same lines
      }
    }

    // ----------------------------------------------------------------------
    // check number of element types is consistent
    if (number_of_atom_types == 0 || number_of_atom_types != number_of_element_dot_in_files) {
      message << "The number of atom types do not match the number of \"ELEMENT.in\" files: atom = " << number_of_atom_types << " vs # .in files = " << number_of_element_dot_in_files;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
    }
    // ----------------------------------------------------------------------
    // check number of total atoms is consistent
    const uint sum_each_type = aurostd::sum(num_each_type);
    if (sum_each_type != a.atoms.size()) {
      message << "The total number of atoms does not match the sum of each atom type: a.atoms.size() = " << a.atoms.size() << " vs sum(num_each_type) = " << sum_each_type;
      ;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ELK READER fixing atom information (alphabetize, make basis, set species, etc.)" << endl;
    }
    a.SpeciesPutAlphabetic();
    std::stable_sort(a.atoms.begin(), a.atoms.end(), sortAtomsNames);
    a.MakeBasis();
    a.MakeTypes();
    // add system name to title
    a.buildGenericTitle(); // DX20210211
    // NO PARTIAL OCCUPATION
    a.partial_occupation_flag = false;
    a.is_vasp4_poscar_format = false;
    a.is_vasp5_poscar_format = false;

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ELK READER - Finished" << endl;
    }
  }

  // ----------------------------------------------------------------------
  // CIF INPUT
  if (a.iomode == IOCIF) {
    // CIF
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CIF [0]" << endl;
    }
    a.scale = 1.0;
    a.neg_scale = false;
    a.lattice = aurostd::eye<double>(3, 3); // CO20190520

    a.spacegroupnumber = 0;
    a.spacegroupnumberoption = 0;
    // get space group
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "_SYMMETRY_INT_TABLES_NUMBER")) {
        // converted to upper to be case insensitive
        vector<string> tokens;
        aurostd::string2tokens(aurostd::RemoveWhiteSpaces(aurostd::toupper(vinput[i])), tokens,
                               "_SYMMETRY_INT_TABLES_NUMBER"); // converted to upper to be case insensitive
        if (tokens.size() == 1) {
          a.spacegroupnumber = aurostd::string2utype<uint>(tokens.at(0));
        }
      } else if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "_SPACE_GROUP_IT_NUMBER")) {
        // converted to upper to be case insensitive
        vector<string> tokens;
        aurostd::string2tokens(aurostd::RemoveWhiteSpaces(aurostd::toupper(vinput[i])), tokens,
                               "_SPACE_GROUP_IT_NUMBER"); // converted to upper to be case insensitive
        if (tokens.size() == 1) {
          a.spacegroupnumber = aurostd::string2utype<uint>(tokens.at(0));
        }
      }
      // DX20190708 - added another space group variant - START
      else if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "_SYMMETRY_SPACE_GROUP_NAME_H-M")) {
        // converted to upper to be case insensitive
        vector<string> tokens;
        aurostd::string2tokens(vinput[i], tokens);
        if (aurostd::toupper(tokens.at(0)) == "_SYMMETRY_SPACE_GROUP_NAME_H-M") {
          tokens.erase(tokens.begin());
          string spacegroupsymbol = aurostd::joinWDelimiter(tokens, "");
          spacegroupsymbol = aurostd::RemoveCharacterFromTheFrontAndBack(spacegroupsymbol, '\''); // clean
          spacegroupsymbol = aurostd::RemoveCharacterFromTheFrontAndBack(spacegroupsymbol, '\"'); // clean
          try {
            a.spacegroupnumber = GetSpaceGroupNumber(spacegroupsymbol);
          } // DX20191029 - added try/catch sequence
          catch (aurostd::xerror& re) {
            if (LDEBUG) {
              message << "Cannot determine space group setting from the Hermann-Mauguin symbol; non-standard setting.";
              pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, std::cerr, _LOGGER_WARNING_);
            }
          } // DX20191029 - added try/catch sequence
        }
      }
      // DX20190708 - added another space group variant - END
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CIF [1]" << endl;
    }
    // get space group setting
    string spacegroup_Hall;
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(aurostd::toupper(vinput[i]), "_SYMMETRY_SPACE_GROUP_NAME_HALL")) {
        // converted to upper to be case insensitive //DX20200521 - missing "_SYMMETRY_" prefix
        vector<string> tokens;
        aurostd::string2tokens(vinput[i], tokens); // DX20200521 - this line was missing
        // DX20190708 - fix Hall reader - START
        if (aurostd::toupper(tokens.at(0)) == "_SYMMETRY_SPACE_GROUP_NAME_HALL") {
          // DX20200521 - changed H-M to HALL
          tokens.erase(tokens.begin());
          spacegroup_Hall = aurostd::joinWDelimiter(tokens, " "); // need a space here for Hall designation
          spacegroup_Hall = aurostd::RemoveCharacterFromTheFrontAndBack(spacegroup_Hall, '\''); // clean
          spacegroup_Hall = aurostd::RemoveCharacterFromTheFrontAndBack(spacegroup_Hall, '\"'); // clean
        }
        // DX20190708 - fix Hall reader - END
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CIF [2]" << endl;
    }
    // DX20191029 - check if space group number is found - START
    if (a.spacegroupnumber == 0) {
      message << "Either space group number was not given or it was given in a non-standard setting.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    // DX20191029 - check if space group number is found - END
    // ME20220124 - Read symmetry operations without consistency checks first
    bool found_setting = false;
    vector<string> spacegroup_symop_xyz;
    bool found_symops = false;
    bool found_symop_id = false; // DX20190708
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(vinput[i], "_space_group_symop_operation_xyz") || aurostd::substring2bool(vinput[i], "_symmetry_equiv_pos_as_xyz")) {
        found_symops = true;
      } else if (aurostd::substring2bool(vinput[i], "_space_group_symop_id") || aurostd::substring2bool(vinput[i], "_symmetry_equiv_pos_site_id")) {
        // DX20190708
        found_symop_id = true;
      } else if (found_symops) {
        if (vinput[i].find("loop_") != string::npos) {
          break; // End of symops
        }
        vector<string> tokens;
        aurostd::string2tokens(vinput[i], tokens, " ");
        // DX20181210 - account for many formats (i.e., x,y,z or 'x, y, z') - START
        if (tokens.empty()) {
          continue;
        }
        if (found_symop_id && !isdigit(tokens[0][0])) {
          break; // End of symops
        }
        if (found_symop_id) {
          tokens.erase(tokens.begin());
        }
        // erase symop index, not needed //DX20190708 - enclose in if-statement
        string symop = aurostd::joinWDelimiter(tokens, "");
        symop = aurostd::RemoveCharacter(symop, '\''); // remove '
        symop = aurostd::RemoveCharacter(symop, '\"'); // remove "
        symop = aurostd::RemoveWhiteSpaces(symop); // remove spaces
        try {
          symop = SYM::reorderWyckoffPosition(symop);
          // DX20190708 - standardize order of equation (variable first, then number)
        } catch (aurostd::xerror& e) {
          break; // End of symops
        }
        spacegroup_symop_xyz.push_back(symop);
        // DX20181210 - account for many formats (i.e., x,y,z or 'x, y, z') - END
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CIF [3]" << endl;
    }
    for (uint setting_number = 1; setting_number <= 2; setting_number++) {
      string setting_string = aurostd::utype2string<uint>(setting_number);
      int general_wyckoff_multiplicity = 0; // general Wyckoff position multiplicity
      vector<string> general_wyckoff_position; // general Wyckoff position equations
      // get general Wyckoff multiplicity and position saved in aflow
      SYM::getGeneralWyckoffMultiplicityAndPosition(a.spacegroupnumber, setting_string, general_wyckoff_multiplicity, general_wyckoff_position);
      // ME20220124 - moved up
      //      SYM::initsgs(setting_string);
      //      using SYM::gl_sgs;
      //      cerr << "find the spacegroupstring" << endl;
      //      string spacegroupstring = gl_sgs[a.spacegroupnumber - 1];
      //      cerr << "spacegroupstring: " << spacegroupstring << endl;
      //      vector<int> wyckoff_multiplicities = SYM::get_multiplicities(spacegroupstring);
      //      cerr << "wyckoff_multiplicity: " << wyckoff_multiplicities[1] << endl;
      //      // get symops from cif
      //      vector<string> spacegroup_symop_xyz;
      //      int multiplicity_count=0;
      //      bool found_symops=false;
      //      cerr << "wyckoff_multiplicities[1]: " << wyckoff_multiplicities[1] << endl;
      //      for(size_t i=0;i<vinput.size();i++) {
      //	if(aurostd::substring2bool(vinput[i],"_space_group_symop_operation_xyz") || aurostd::substring2bool(vinput[i],"_symmetry_equiv_pos_as_xyz")){ // _space_group_symop_operation_xyz supersedes all
      //	  found_symops=true;
      //	}
      //	else if(found_symops && multiplicity_count<wyckoff_multiplicities[1]){
      //	  multiplicity_count+=1;
      //	  vector<string> tokens;
      //	  aurostd::string2tokens(vinput[i],tokens," ");
      //	  if(tokens.size()==2){
      //	    spacegroup_symop_xyz.push_back(tokens[1]);
      //	  }
      //
      //	}
      //      }
      //      // check symops against general wyckoff position to determine setting
      //      //uint setting_number = 1; //default is first setting
      //      vector<string> general_wyckoff_position = SYM::findGeneralWyckoffPosition(spacegroupstring, wyckoff_multiplicities[1]);
      //      cerr << "general position" << endl;
      //      print(general_wyckoff_position);
      //      cerr << "====================================" << endl;
      //      cerr << "cif symop" << endl;
      //      print(spacegroup_symop_xyz);
      //      cerr << "general_wyckoff_position.size(): " << general_wyckoff_position.size() << endl;
      //      cerr << "spacegroup_symop_xyz.size(): " << spacegroup_symop_xyz.size() << endl;
      // check equations in cif
      // vector<string> spacegroup_symop_xyz;
      // int multiplicity_count=0;
      // bool found_symops=false;
      // bool found_symop_id=false; //DX20190708
      // for(size_t i=0;i<vinput.size();i++) {
      //  if(aurostd::substring2bool(vinput[i],"_space_group_symop_operation_xyz") || aurostd::substring2bool(vinput[i],"_symmetry_equiv_pos_as_xyz")){
      //    found_symops=true;
      //  }
      //  else if(aurostd::substring2bool(vinput[i],"_space_group_symop_id") || aurostd::substring2bool(vinput[i],"_symmetry_equiv_pos_site_id")){ //DX20190708
      //    found_symop_id=true;
      //  }
      //  else if(found_symops && multiplicity_count<general_wyckoff_multiplicity){
      //    multiplicity_count+=1;
      //    vector<string> tokens;
      //    aurostd::string2tokens(vinput[i],tokens," ");
      //    //DX20181210 - account for many formats (i.e., x,y,z or 'x, y, z') - START
      //    if(found_symop_id){ tokens.erase(tokens.begin()); } //erase symop index, not needed //DX20190708 - enclose in if-statement
      //    string symop = aurostd::joinWDelimiter(tokens,"");
      //    symop = aurostd::RemoveCharacter(symop,'\''); // remove '
      //    symop = aurostd::RemoveCharacter(symop,'\"'); // remove "
      //    symop = aurostd::RemoveWhiteSpaces(symop); // remove spaces
      //    symop = SYM::reorderWyckoffPosition(symop); //DX20190708 - standardize order of equation (variable first, then number)
      //    spacegroup_symop_xyz.push_back(symop);
      //    //if(tokens.size()==2){
      //    //  spacegroup_symop_xyz.push_back(tokens[1]);
      //    //}
      //    //DX20181210 - account for many formats (i.e., x,y,z or 'x, y, z') - END
      //  }
      //}
      // compare cif and aflow's general position
      uint match_count = 0;
      if (general_wyckoff_position.size() == spacegroup_symop_xyz.size()) {
        for (size_t i = 0; i < general_wyckoff_position.size(); i++) {
          bool match_found = false;
          for (size_t j = 0; j < spacegroup_symop_xyz.size(); j++) {
            if (general_wyckoff_position[i] == spacegroup_symop_xyz[j]) {
              match_found = true;
              match_count += 1;
              break;
            }
          }
          if (!match_found) {
            if (LDEBUG) {
              cerr << "WARNING: Could not match " << i << ": " << general_wyckoff_position[i] << " to any symop in CIF. Trying other setting (if exists)." << endl;
            }
            break;
          }
        }
        if (match_count != general_wyckoff_position.size()) {
          // oss << "ERROR - xstructure::operator>>: Symmetry operations do not match between input operations and space group number/option." << endl;
          // print(general_wyckoff_position);
          // print(spacegroup_symop_xyz);
          // throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_);
          continue; // try a different setting
        } else {
          a.spacegroupnumberoption = setting_number;
          a.spacegroupoption = setting_string; // DX20191029
          found_setting = true;
          break; // found setting
        }
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CIF [4]" << endl;
    }
    if (!found_setting) {
      // ME20220124 - Changed to warning
      message << "Symmetry operations do not match between input operations and space group number/option.";
      // CO20190629
      message << " Building structure using symmetry operations in CIF file with space group P1.";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, std::cerr, _LOGGER_WARNING_);
      a.spacegroupnumber = 1;
      if (LDEBUG) {
        for (size_t i = 0; i < vinput.size(); i++) {
          std::cerr << vinput[i] << endl;
        }
        // CO20190629 // ME20220124 - moved to LDEBUG because outputting the entire CIF makes the error message unreadable
      }
      // throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CIF [5]" << endl;
    }
    // get lattice
    for (size_t i = 0; i < vinput.size(); i++) {
      vector<string> tokens;
      aurostd::string2tokens(vinput[i], tokens, " ");
      if (aurostd::substring2bool(vinput[i], "_cell_length_a") && tokens.size() == 2) {
        a.a = aurostd::string2utype<double>(tokens[1]);
      } else if (aurostd::substring2bool(vinput[i], "_cell_length_b") && tokens.size() == 2) {
        a.b = aurostd::string2utype<double>(tokens[1]);
      } else if (aurostd::substring2bool(vinput[i], "_cell_length_c") && tokens.size() == 2) {
        a.c = aurostd::string2utype<double>(tokens[1]);
      } else if (aurostd::substring2bool(vinput[i], "_cell_angle_alpha") && tokens.size() == 2) {
        a.alpha = aurostd::string2utype<double>(tokens[1]);
      } else if (aurostd::substring2bool(vinput[i], "_cell_angle_beta") && tokens.size() == 2) {
        a.beta = aurostd::string2utype<double>(tokens[1]);
      } else if (aurostd::substring2bool(vinput[i], "_cell_angle_gamma") && tokens.size() == 2) {
        a.gamma = aurostd::string2utype<double>(tokens[1]);
      }
    }
    a.lattice = GetClat(a.a, a.b, a.c, a.alpha, a.beta, a.gamma);
    a.FixLattices();
    a.partial_occupation_flag = false;

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CIF [6]" << endl;
    }
    // ME20220124 - Convert xyz to symbolic representation
    vector<symbolic::Symbolic> spacegroup_symop_symbolic;
    if (!found_setting) {
      const size_t nsym = spacegroup_symop_xyz.size();
      vector<vector<string>> vsymops(nsym);
      for (uint i = 0; i < nsym; ++i) {
        aurostd::string2tokens(spacegroup_symop_xyz[i], vsymops[i], ",");
      }
      spacegroup_symop_symbolic = anrl::equations2SymbolicEquations(vsymops);
    }

    // get atoms
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CIF [7]" << endl;
    }
    vector<string> atom_site_fields;
    bool found_atom_site_labels = false;
    bool already_storing_atoms = false;
    for (size_t i = 0; i < vinput.size(); i++) {
      // DX20181210 - CIF can have multiple loops, reset after each loop and neglect aniso loop - START
      if (aurostd::substring2bool(vinput[i], "loop")) {
        found_atom_site_labels = false;
      }
      if (aurostd::substring2bool(vinput[i], "_atom_site") && !aurostd::substring2bool(vinput[i], "aniso")) {
        if (already_storing_atoms) {
          break;
        }
        // DX20190718 - to handle format of Springer Materials cifs (adds extra fields at the end; do not read them)
        found_atom_site_labels = true;
        atom_site_fields.push_back(vinput[i]);
      } else if (found_atom_site_labels && !aurostd::RemoveWhiteSpaces(vinput[i]).empty() && vinput[i][0] == '_' && !aurostd::substring2bool(vinput[i], "aniso")) {
        // DX20190718 - to account for non-standard fields, e.g., _sm_ (Springer Materials) fields; not standard, but abundant
        if (already_storing_atoms) {
          break;
        }
        // DX20190718 - to handle format of Springer Materials cifs (adds extra fields at the end; do not read them)
        atom_site_fields.push_back(vinput[i]);
      }
      // DX20181210 - CIF can have multiple loops, reset after each loop and neglect aniso_U loop - END
      else if (found_atom_site_labels == true && !aurostd::RemoveWhiteSpaces(vinput[i]).empty()) {
        vector<string> tokens;
        aurostd::string2tokens(vinput[i], tokens, " ");
        // DX20190718 - check tokens first; Springer Materials has extra spaces - START
        for (size_t t = 0; t < tokens.size(); t++) {
          if (aurostd::RemoveWhiteSpaces(tokens[t])[0] == '\'' && aurostd::RemoveWhiteSpaces(tokens[t])[tokens[t].size() - 1] != '\'') {
            if (t + 1 < tokens.size()) {
              if (aurostd::RemoveWhiteSpaces(tokens[t + 1])[0] != '\'' && aurostd::RemoveWhiteSpaces(tokens[t + 1])[tokens[t + 1].size() - 1] == '\'') {
                tokens[t] += "_" + tokens[t + 1];
                tokens.erase(tokens.begin() + t + 1);
              }
            }
          }
        }
        // DX20190718 - check tokens first; Springer Materials has extra spaces - END
        if (tokens.size() == atom_site_fields.size()) {
          _atom atom_tmp;
          wyckoffsite_ITC wyckoff_tmp; // DX20191029
          // ME20220124 - prepare for P1 if setting not found
          if (!found_setting) {
            wyckoff_tmp.letter = "a";
            wyckoff_tmp.multiplicity = 1;
            wyckoff_tmp.site_symmetry = "1";
            const vector<string> eq{"x", "y", "z"};
            wyckoff_tmp.equations.push_back(eq);
          }
          for (size_t t = 0; t < tokens.size(); t++) {
            if (aurostd::substring2bool(atom_site_fields.at(t), "_atom_site_type_symbol")) {
              string name = aurostd::RemoveCharacterFromTheFrontAndBack(tokens[t], '\'');
              // DX20190718 - remove surrounding '' (common in Springer Materials cifs)
              //  ME20220113 - name could have oxidation states (found in newer ICSD CIFs)
              if (name.length() > 2) {
                name = name.substr(0, isdigit(name[1]) ? 1 : 2);
              }
              atom_tmp.name = name;
              atom_tmp.name_is_given = true;
            }
            if (aurostd::substring2bool(atom_site_fields.at(t), "_atom_site_fract_x")) {
              atom_tmp.fpos[1] = aurostd::string2utype<double>(tokens[t]);
            }
            if (aurostd::substring2bool(atom_site_fields.at(t), "_atom_site_fract_y")) {
              atom_tmp.fpos[2] = aurostd::string2utype<double>(tokens[t]);
            }
            if (aurostd::substring2bool(atom_site_fields.at(t), "_atom_site_fract_z")) {
              atom_tmp.fpos[3] = aurostd::string2utype<double>(tokens[t]);
            }
            if (aurostd::substring2bool(atom_site_fields.at(t), "_atom_site_occupancy")) {
              atom_tmp.partial_occupation_value = aurostd::string2utype<double>(tokens[t]);
              wyckoff_tmp.site_occupation = atom_tmp.partial_occupation_value; // DX20191029
              if (aurostd::abs(atom_tmp.partial_occupation_value - 1.0) < 1e-6) {
                atom_tmp.partial_occupation_flag = false;
              } else {
                atom_tmp.partial_occupation_flag = true;
                a.partial_occupation_flag = true;
              }
            }
            if (aurostd::substring2bool(atom_site_fields.at(t), "_atom_site_symmetry_multiplicity")) {
              wyckoff_tmp.multiplicity = aurostd::string2utype<double>(tokens[t]);
            }
            if (aurostd::substring2bool(atom_site_fields.at(t), "_atom_site_Wyckoff_label") || aurostd::substring2bool(atom_site_fields.at(t), "_atom_site_Wyckoff_symbol")) {
              wyckoff_tmp.letter = tokens[t];
            } // ME20220125 - Added _atoms_site_Wyckoff_symbol check
          }
          wyckoff_tmp.type = atom_tmp.name; // DX20191029
          wyckoff_tmp.coord = atom_tmp.fpos; // DX20191029
          if (found_setting && wyckoff_tmp.multiplicity != 0 && !wyckoff_tmp.letter.empty()) {
            SYM::getWyckoffInformation(a.spacegroupnumber, a.spacegroupoption, wyckoff_tmp.letter, wyckoff_tmp.multiplicity, wyckoff_tmp.site_symmetry, wyckoff_tmp.equations);
          }
          a.wyckoff_sites_ITC.push_back(wyckoff_tmp);
          atom_tmp.cpos = a.f2c * atom_tmp.fpos;
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " pre AddAtom() [i=" << i << "]" << endl;
          }
          a.AddAtom(atom_tmp, true); // CO20230319 - add by species
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " post AddAtom() [i=" << i << "]" << endl;
          }
          // ME20220124 - Use symmetry operations when setting unknown
          if (!found_setting) {
            const size_t natoms = a.atoms.size(); // For Wyckoff positions
            _atom at;
            deque<_atom> atoms_symop;
            const symbolic::Symbolic x("x");
            const symbolic::Symbolic y("y");
            const symbolic::Symbolic z("z");
            symbolic::Symbolic result;
            for (size_t i = 0; i < spacegroup_symop_symbolic.size(); i++) {
              const symbolic::Symbolic& eq = spacegroup_symop_symbolic[i];
              const xvector<double>& fpos = atom_tmp.fpos;
              result = eq[x == fpos[1], y == fpos[2], z == fpos[3]];
              at = atom_tmp;
              at.fpos[1] = (double) result(0);
              at.fpos[2] = (double) result(1);
              at.fpos[3] = (double) result(2);
              at.cpos = a.f2c * at.fpos;
              atoms_symop.push_back(at);
            }
            a.AddAtom(atoms_symop, true); // CO20230319 - add by species
            // Add Wyckoff positions
            const size_t natoms_added = a.atoms.size() - natoms;
            for (uint i = 0; i < natoms_added; ++i) {
              a.wyckoff_sites_ITC.push_back(wyckoff_tmp);
            }
          }
          already_storing_atoms = true;
          // DX20190718 - to handle format of Springer Materials cifs (adds extra fields at the end; do not read them)
        } else {
          message << "Unexpected number of input fields based on _atom_site_[] information (tokens=" << tokens.size() << ", atom_sites_[]=" << atom_site_fields.size() << ")." << endl; // CO20190629
          for (size_t i = 0; i < vinput.size(); i++) {
            message << vinput[i] << endl; // CO20190629
          }
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
        }
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CIF [8]" << endl;
    }
    a = WyckoffPOSITIONS(a.spacegroupnumber, a.spacegroupnumberoption, a);
    a.isd = false; // set Selective Dynamics to false
    // DX20191010 - moved loop that used to be here after re-alphabetizing
    a.SpeciesPutAlphabetic(); // DX20190508 - put alphabetic, needed for many AFLOW functions to work properly
    std::stable_sort(a.atoms.begin(), a.atoms.end(), sortAtomsNames); // DX20200312
    std::sort(a.wyckoff_sites_ITC.begin(), a.wyckoff_sites_ITC.end(), sortWyckoffByType);
    // DX20201014 - sort the Wyckoff positions too
    a.MakeBasis(); // DX20200803 - must be after alphabetic sort
    a.MakeTypes(); // DX20190508 - otherwise types are not created //DX20200803 - must be after alphabetic sort
    // DX20191010 - moved this loop - START
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CIF [9]" << endl;
    }
    for (size_t i = 0; i < a.atoms.size(); i++) {
      if (a.atoms[i].partial_occupation_flag == true) {
        poccaus.push_back(a.atoms[i].partial_occupation_value);
        a.partial_occupation_sublattice.push_back(a.atoms[i].type);
      } else {
        poccaus.push_back(1.0);
        a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
      }
    }
    // DX20191010 - moved this loop - END
    a.is_vasp4_poscar_format = false; // DX20190308 - needed or SPECIES section breaks
    a.is_vasp5_poscar_format = false; // DX20190308 - needed or SPECIES section breaks

    // add title, CIFs do not generally have a canonical "title" line, so make one
    a.BringInCell(); // ME20220124
    a.buildGenericTitle(); // DX20210211
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CIF [DONE]" << endl;
    }
  } // CIF INPUT

  // ----------------------------------------------------------------------
  // AIMS INPUT
  if (a.iomode == IOAIMS_AUTO || a.iomode == IOAIMS_GEOM) {
    // AIMS GEOM
    // CO20171008
    // remember, we already did all the debugging above
    // if we get here, we can assume geometry.in is solid!
    // if more debugging needed, test above where we detect iomode
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " AIMS IOAIMS_AUTO/IOAIMS_GEOM" << endl;
    }
    // START FROM CELL
    a.scale = 1.0; // standard
    a.neg_scale = false; // standard
    a.lattice = aurostd::eye<double>(3, 3); // CO20190520

    uint lat_count = 0;
    // get lattice first, if available (c2f,f2c)
    for (size_t i = 0; i < vinput.size(); i++) {
      aurostd::string2tokens(vinput[i], tokens, " ");
      if (tokens[0] == "lattice_vector" && tokens.size() > 3) {
        lat_count++;
        a.lattice(lat_count, 1) = aurostd::string2utype<double>(tokens[1]);
        a.lattice(lat_count, 2) = aurostd::string2utype<double>(tokens[2]);
        a.lattice(lat_count, 3) = aurostd::string2utype<double>(tokens[3]);
      }
      if (lat_count == 3) {
        break;
      }
    }

    a.FixLattices();
    // a.f2c=trasp(a.lattice), a.c2f=inverse(a.f2c);

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " AIMS lattice" << endl;
      cerr << a.lattice << endl;
      cerr << __AFLOW_FUNC__ << " AIMS f2c" << endl;
      cerr << a.f2c << endl;
      cerr << __AFLOW_FUNC__ << " AIMS c2f" << endl;
      cerr << a.c2f << endl;
    }

    deque<_atom> atoms;
    // now get atoms
    bool atom_found;
    _atom atom;
    bool atom_props_search = false; // AIMS stores more atom properties UNDER "atom" and before next "atom"
    for (size_t i = 0; i < vinput.size(); i++) {
      atom_found = false;
      aurostd::string2tokens(vinput[i], tokens, " ");
      if (tokens.empty()) {
        continue;
      }
      if (aurostd::substring2bool(tokens[0], "atom") && tokens.size() > 4) {
        atom_props_search = false;
        if (tokens[0] == "atom") {
          atom_found = true;
          atom.clear();
          atom.cpos[1] = aurostd::string2utype<double>(tokens[1]);
          atom.cpos[2] = aurostd::string2utype<double>(tokens[2]);
          atom.cpos[3] = aurostd::string2utype<double>(tokens[3]);
          atom.fpos = a.c2f * atom.cpos;
          atom.name = atom.cleanname = tokens[4];
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " AIMS atom[" << atom.name << "] found (cartesian):" << endl;
            cerr << "    cpos" << atom.cpos << endl;
            cerr << "    fpos" << atom.fpos << endl;
          }
        } else if (tokens[0] == "atom_frac") {
          atom_found = true;
          atom.clear();
          atom.fpos[1] = aurostd::string2utype<double>(tokens[1]);
          atom.fpos[2] = aurostd::string2utype<double>(tokens[2]);
          atom.fpos[3] = aurostd::string2utype<double>(tokens[3]);
          atom.cpos = a.f2c * atom.fpos;
          atom.name = atom.cleanname = tokens[4];
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " AIMS atom[" << atom.name << "] found (fractional):" << endl;
            cerr << "    fpos" << atom.fpos << endl;
            cerr << "    cpos" << atom.cpos << endl;
          }
        } // ignore else, garbage
        if (atom_found) {
          atom_props_search = true; // so we can look for other props on the next line

          // set some defaults! - START
          atom.name_is_given = true;
          // atom.CleanName();
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " AIMS line=" << vinput[i] << endl;
            cerr << __AFLOW_FUNC__ << " AIMS atom.cleanname=" << atom.cleanname << endl;
          }
          atom.CleanSpin();
          // FIXED BELOW atom.number=atoms.size();    // reference position for convasp
          // FIXED BELOW atom.basis=atoms.size();     // position in the basis //MAKE SURE NOT TO SET BASIS HERE, WILL SCREW UP sortAtomsNames() BELOW
          atom.ijk(1) = 0;
          atom.ijk(2) = 0;
          atom.ijk(3) = 0; // inside the zero cell...
          atom.corigin(1) = 0.0;
          atom.corigin(2) = 0.0;
          atom.corigin(3) = 0.0; // inside the zero cell
          atom.coord(1) = 0.0;
          atom.coord(2) = 0.0;
          atom.coord(3) = 0.0; // inside the zero cell
          atom.spin = 0.0;
          atom.noncoll_spin.clear(); // DX20171205 - non-collinear spin
          // FIXED BELOW atom.type=itype;                // CONVASP_MODE if I want type 0,1,2,3,...
          atom.order_parameter_value = 0;
          atom.order_parameter_atom = false;
          atom.partial_occupation_value = 1.0;
          atom.partial_occupation_flag = false;
          // set some defaults! - STOP

          atoms.push_back(atom);
          // wait till later to add to structure, sort first that we don't screw up species
          // a.AddAtom(atom,true); //CO20230319 - add by species
          // a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);  //default
        }
      } else if (atom_props_search && !atoms.empty()) {
        // look and store additional properties here to atoms.back();
      } // ignore else
    }

    // sort first, then assign types
    // MUST BE STABLE SORT, absolutely critical
    // otherwise the order of equivalent names will be mixed EVERY TIME
    // this will screw up forces in APL
    std::stable_sort(atoms.begin(), atoms.end(), sortAtomsNames); // safe because we do AddAtom() below

    // CO20230319 - DX, I believe this code is obsolete now that I patch AddAtom() with add_species, keep for safety
    // assign types
    uint itype = 0;
    atoms[0].type = itype;
    for (size_t i = 1; i < atoms.size(); i++) {
      if (atoms[i].name != atoms[i - 1].name) {
        itype++;
      }
      atoms[i].type = itype;
    }

    // a.atoms.clear();  //we clear atoms earlier with RemoveAtom()
    // now add atoms in right order (for species, etc.)
    for (size_t i = 0; i < atoms.size(); i++) {
      a.AddAtom(atoms[i], true); // does num_each_type and comp_each_type  //CO20230319 - add by species
      a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
    }

    a.MakeBasis();

    if (LDEBUG) {
      for (size_t i = 0; i < a.num_each_type.size(); i++) {
        cerr << __AFLOW_FUNC__ << " AIMS num_each_type[" << i << "]=" << a.num_each_type[i] << ", ";
        cerr << "comp_each_type[" << i << "]=" << a.comp_each_type[i] << endl;
      }
    }

    // grab title last
    a.title.clear();
    if (!vinput.empty()) {
      const std::size_t pos = vinput[0].find_first_of('#');
      if (pos != std::string::npos) {
        a.title = vinput[0];
        a.title.erase(pos, 1);
        a.title = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(a.title);
      }
    }
    if (a.title.empty()) {
      a.buildGenericTitle();
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " AIMS title=" << a.title << endl;
    }

    a.partial_occupation_flag = false;
    a.is_vasp4_poscar_format = false;
    a.is_vasp5_poscar_format = false;
    // DONE ?
  } // AIMS INPUT
  // ----------------------------------------------------------------------

  // ----------------------------------------------------------------------
  // ATAT INPUT //SD20220114
  // Alloy-Theoretic Automated Toolkit
  // See: https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/manual.pdf
  if (a.iomode == IOATAT_STR) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ATAT IOATAT_STR" << endl;
    }
    a.scale = 1.0;
    a.neg_scale = false;
    const xmatrix<double> axes(3, 3);
    const xmatrix<double> frac_cell(3, 3);

    // read the axes
    size_t line = 0;
    uint vec_count = 1;
    for (; line < vinput.size() && vec_count < 4; line++) {
      aurostd::string2tokens(vinput[line], tokens);
      axes(vec_count, 1) = aurostd::string2utype<double>(tokens[0]);
      axes(vec_count, 2) = aurostd::string2utype<double>(tokens[1]);
      axes(vec_count, 3) = aurostd::string2utype<double>(tokens[2]);
      vec_count++;
    }
    // read the fractional cell vectors
    vec_count = 1;
    for (; line < vinput.size() && vec_count < 4; line++) {
      aurostd::string2tokens(vinput[line], tokens);
      frac_cell(vec_count, 1) = aurostd::string2utype<double>(tokens[0]);
      frac_cell(vec_count, 2) = aurostd::string2utype<double>(tokens[1]);
      frac_cell(vec_count, 3) = aurostd::string2utype<double>(tokens[2]);
      vec_count++;
    }
    a.lattice = axes * frac_cell;
    a.FixLattices();
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ATAT lattice" << endl;
      cerr << a.lattice << endl;
      cerr << __AFLOW_FUNC__ << " ATAT f2c" << endl;
      cerr << a.f2c << endl;
      cerr << __AFLOW_FUNC__ << " ATAT c2f" << endl;
      cerr << a.c2f << endl;
    }

    // read atoms
    deque<_atom> atoms;
    _atom atom;
    const xvector<double> avec(3);
    for (; line < vinput.size() - 1; line++) {
      atom.clear();
      aurostd::string2tokens(vinput[line], tokens);
      avec(1) = aurostd::string2utype<double>(tokens[0]);
      avec(2) = aurostd::string2utype<double>(tokens[1]);
      avec(3) = aurostd::string2utype<double>(tokens[2]);
      atom.name = atom.cleanname = tokens[3];
      atom.cpos = trasp(axes) * avec;
      atom.fpos = a.c2f * atom.cpos;
      atom.name_is_given = (!atom.name.empty());
      atoms.push_back(atom);
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " ATAT atom[" << atom.name << "] found:" << endl;
        cerr << "    fpos" << atom.fpos << endl;
        cerr << "    cpos" << atom.cpos << endl;
      }
    }
    std::stable_sort(atoms.begin(), atoms.end(), sortAtomsNames);

    // CO20230319 - DX, I believe this code is obsolete now that I patch AddAtom() with add_species, keep for safety
    //  assign types
    uint itype = 0;
    atoms[0].type = itype;
    for (size_t i = 1; i < atoms.size(); i++) {
      if (atoms[i].name != atoms[i - 1].name) {
        itype++;
      }
      atoms[i].type = itype;
    }

    // add atoms
    for (size_t i = 0; i < atoms.size(); i++) {
      a.AddAtom(atoms[i], true); // CO20230319 - add by species
      a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
    }
    a.SpeciesPutAlphabetic();

    // add additional attributes
    a.MakeBasis();
    a.MakeTypes();
    a.partial_occupation_flag = false;
    a.is_vasp4_poscar_format = false;
    a.is_vasp5_poscar_format = false;
    a.buildGenericTitle();
  } // ATAT INPUT

  // ----------------------------------------------------------------------
  // LMP INPUT //SD20240111
  if (a.iomode == IOLMP_DATA) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " LMP IOLMP_DATA" << endl;
    }

    // find line indices
    size_t i_units = 0;
    size_t i_natoms = 0;
    size_t i_cell = 0;
    size_t i_atoms = 0;
    for (size_t i = 0; i < vinput.size() && !i_atoms; i++) {
      if (aurostd::substring2bool(vinput[i], "LAMMPS data file via write_data")) {
        i_units = i;
      } else if (aurostd::substring2bool(vinput[i], "atoms")) {
        i_natoms = i;
      } else if (aurostd::substring2bool(vinput[i], "xlo xhi")) {
        i_cell = i;
      } else if (aurostd::substring2bool(vinput[i], "Atoms")) {
        i_atoms = i;
      }
    }
    if (i_natoms == 0 || i_cell == 0 || i_atoms == 0) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Invalid LAMMPS data file format", _FILE_WRONG_FORMAT_);
    }

    // read the units
    a.neg_scale = false;
    const string units = aurostd::RemoveWhiteSpaces(vinput[i_units].substr(vinput[0].find("units = ") + 8));
    if (units == "real" || units == "metal") {
      a.scale = 1.0;
    } else if (units == "electron") {
      a.scale = bohr2angstrom;
    } else if (units == "micro") {
      a.scale = 1E4;
    } else if (units == "nano") {
      a.scale = 1E1;
    } else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Undefined units", _INPUT_UNKNOWN_);
    }

    // read the number of atoms
    aurostd::string2tokens(vinput[i_natoms], tokens);
    const uint natoms = aurostd::string2utype<uint>(tokens[0]);

    // read the cell
    double lx = 0.0;
    double ly = 0.0;
    double lz = 0.0;
    double xy = 0.0;
    double xz = 0.0;
    double yz = 0.0;
    aurostd::string2tokens(vinput[i_cell + 0], tokens);
    lx = aurostd::string2utype<double>(tokens[1]) - aurostd::string2utype<double>(tokens[0]);
    aurostd::string2tokens(vinput[i_cell + 1], tokens);
    ly = aurostd::string2utype<double>(tokens[1]) - aurostd::string2utype<double>(tokens[0]);
    aurostd::string2tokens(vinput[i_cell + 2], tokens);
    lz = aurostd::string2utype<double>(tokens[1]) - aurostd::string2utype<double>(tokens[0]);
    aurostd::string2tokens(vinput[i_cell + 3], tokens);
    if (!tokens.empty()) {
      xy = aurostd::string2utype<double>(tokens[0]);
      xz = aurostd::string2utype<double>(tokens[1]);
      yz = aurostd::string2utype<double>(tokens[2]);
    }
    a.lattice = aurostd::xmatrix<double>{
        {lx, 0.0, 0.0},
        {xy,  ly, 0.0},
        {xz,  yz,  lz}
    };
    a.FixLattices();
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " LMP lattice" << endl;
      cerr << a.lattice << endl;
    }

    // read the atoms
    aurostd::string2tokens(vinput[i_atoms + 2], tokens);
    bool name_is_given = false;
    if (!aurostd::string2utype<uint>(tokens[1])) {
      name_is_given = true;
    }
    deque<_atom> atoms;
    _atom atom;
    const xvector<double> avec(3);
    for (size_t i = i_atoms + 2; i < i_atoms + 2 + natoms; i++) {
      aurostd::string2tokens(vinput[i], tokens);
      atom.clear();
      atom.name_is_given = name_is_given;
      atom.name = atom.cleanname = tokens[1];
      avec(1) = aurostd::string2utype<double>(tokens[2]);
      avec(2) = aurostd::string2utype<double>(tokens[3]);
      avec(3) = aurostd::string2utype<double>(tokens[4]);
      atom.cpos = avec;
      atom.fpos = a.c2f * atom.cpos;
      atoms.push_back(atom);
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " LMP atom[" << atom.name << "] found:" << endl;
        cerr << "    fpos" << atom.fpos << endl;
        cerr << "    cpos" << atom.cpos << endl;
      }
    }
    std::stable_sort(atoms.begin(), atoms.end(), sortAtomsNames);
    for (size_t i = 0; i < atoms.size(); i++) {
      a.AddAtom(atoms[i], true);
      a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
    }
    a.SpeciesPutAlphabetic();
    a.BringInCell();

    // add additional attributes
    a.MakeBasis();
    a.MakeTypes();
    a.partial_occupation_flag = false;
    a.is_vasp4_poscar_format = false;
    a.is_vasp5_poscar_format = false;
    a.buildGenericTitle();
  } // LMP INPUT

  // ----------------------------------------------------------------------

  // COMMON CODE (NEED TO BE PATCHED, THOUGH).
  // FIX NORMAL AND PARTIAL OCCUPAITON
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " COMMON CODE [9]" << endl;
  }
  if (a.partial_occupation_flag == false) {
    // have partial
    a.comp_each_type.clear();
    for (size_t i = 0; i < a.num_each_type.size(); i++) {
      a.comp_each_type.push_back((double) a.num_each_type[i]);
    }
  } else {
    // have partial
    if (poccaus.size() != a.atoms.size()) {
      message << "poccaus.size()=" << poccaus.size() << " a.atoms.size()=" << a.atoms.size() << " " << endl;
      // CO20190629
      for (size_t i = 0; i < vinput.size(); i++) {
        message << vinput[i] << endl; // CO20190629
      }
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
    }
    // create list (empty)
    a.comp_each_type.clear();
    for (size_t i = 0; i < a.num_each_type.size(); i++) {
      a.comp_each_type.push_back(0.0);
    }
    // screen each one
    for (size_t i = 0; i < a.atoms.size(); i++) {
      if (poccaus.at(i) < 1.0 - 1e-5) {
        a.atoms[i].partial_occupation_flag = true;
        a.atoms[i].partial_occupation_value = poccaus.at(i);
        a.comp_each_type.at(a.atoms[i].type) += a.atoms[i].partial_occupation_value;
      } else {
        a.atoms[i].partial_occupation_flag = false;
        a.atoms[i].partial_occupation_value = 1.0;
        a.comp_each_type.at(a.atoms[i].type) += a.atoms[i].partial_occupation_value;
      }
    }
  }
  a.GetStoich(); // CO20170724
  // ---------------------------------------------------------------
  // -------------- CREATE PARTIAL OCCUPATION STUFF
  // cerr << a.atoms.size() << endl;
  // cerr << a.partial_occupation_sublattice.size() << endl;
  // for(size_t i=0;i<a.partial_occupation_sublattice.size();i++) cerr << a.partial_occupation_sublattice.at(i) << endl;
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " PARTIAL OCCUPATION STUFF [10]" << endl;
  }
  for (size_t iatom = 0; iatom < a.atoms.size(); iatom++) {
    vector<uint> partial_occupation_sublattice_iatom;
    if (a.partial_occupation_sublattice.at(iatom) == _pocc_no_sublattice_) {
      partial_occupation_sublattice_iatom.push_back(iatom);
    } else {
      for (size_t jatom = 0; jatom < a.atoms.size(); jatom++) {
        if (a.atoms[jatom].type == a.partial_occupation_sublattice.at(iatom)) {
          partial_occupation_sublattice_iatom.push_back(iatom);
        }
      }
    }
  }
  // ---------------------------------------------------------------
  // REMOVE ORDER PARAMETER
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " REMOVE ORDER PARAMETER [11]" << endl;
  }
  a.order_parameter_atoms.clear();
  for (size_t i = 0; i < a.atoms.size(); i++) {
    if (a.atoms[i].order_parameter_atom == true) {
      a.order_parameter_atoms.push_back(i);
    }
  }
  a.order_parameter_structure = (!a.order_parameter_atoms.empty());
  // ---------------------------------------------------------------
  // -------------- SPECIES
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " SPECIES [12]" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " a.is_vasp4_poscar_format=" << a.is_vasp4_poscar_format << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " a.is_vasp5_poscar_format=" << a.is_vasp5_poscar_format << endl;
  }
  if (a.is_vasp4_poscar_format == true) {
    a.species.clear();
    a.species_pp.clear();
    a.species_pp_type.clear();
    a.species_pp_version.clear();
    a.species_pp_ZVAL.clear();
    a.species_pp_vLDAU.clear();
    a.species_volume.clear();
    a.species_mass.clear();
    // Plug the species
    for (size_t i = 0, j = 0; i < a.num_each_type.size(); j += a.num_each_type.at(i), i++) {
      a.species.push_back(a.atoms.at(j).name);
      a.species_pp.push_back(a.atoms.at(j).name);
      a.species_pp_type.emplace_back("");
      a.species_pp_version.emplace_back("");
      a.species_pp_ZVAL.push_back(0.0);
      a.species_pp_vLDAU.emplace_back();
      a.species_volume.push_back(0.0);
      a.species_mass.push_back(0.0);
    }
  }
  if (a.is_vasp5_poscar_format == true) {
    a.species.clear();
    a.species_pp.clear();
    a.species_pp_type.clear();
    a.species_pp_version.clear();
    a.species_pp_ZVAL.clear();
    a.species_pp_vLDAU.clear();
    a.species_volume.clear();
    a.species_mass.clear();
    if (a.iomode == IOVASP_POSCAR) {
      aurostd::string2tokens((vinput[5]), tokens);
    }
    if (a.iomode == IOVASP_ABCCAR) {
      aurostd::string2tokens((vinput[3]), tokens);
    }
    if (a.iomode == IOVASP_WYCKCAR) {
      aurostd::string2tokens((vinput[3]), tokens);
    }
    if (a.num_each_type.size() != tokens.size()) {
      message << "You need to specify the same number of species and atoms types" << endl; // CO20190629
      message << "      a.num_each_type.size()=" << a.num_each_type.size() << endl; // CO20190629
      message << "      tokens.size()=" << tokens.size() << endl; // CO20190629
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
    }
    for (size_t i = 0; i < tokens.size(); i++) {
      a.species.push_back(tokens[i]);
      a.species_pp.push_back(tokens[i]);
      a.species_pp_type.emplace_back("");
      a.species_pp_version.emplace_back("");
      a.species_pp_ZVAL.push_back(0.0);
      a.species_pp_vLDAU.emplace_back();
      a.species_volume.push_back(0.0);
      a.species_mass.push_back(0.0);
    }
    uint k = 0;
    for (size_t i = 0; i < a.num_each_type.size(); i++) {
      for (uint j = 0; j < (uint) a.num_each_type[i]; j++) {
        if (a.atoms.at(k).name_is_given == false) {
          a.atoms.at(k).name = a.species.at(i);
          a.atoms.at(k).name_is_given = true;
        }
        k++;
      }
    }
  }
  // ---------------------------------------------------------------
  // ALL atoms have been added. Now add the wyckoff ones
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " WYCKCAR [13]" << endl;
  }
  if (a.iomode == IOVASP_WYCKCAR) {
    // cerr << "[" << a.atoms.size() << "]" << endl;
    a = WyckoffPOSITIONS(a.spacegroupnumber, a.spacegroupnumberoption, a);
    a.title = a.title + " (WYCKOFF " + a.spacegroup + " " + a.spacegrouplabel + ")";
    // cerr << "[" << a.atoms.size() << "]" << endl;
  }

  // ---------------------------------------------------------------
  // Make spaces and links inside the qm part
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " QM_SPACE [14]" << endl;
  }
  for (size_t i = 0; i < a.atoms.size(); i++) {
    xvector<double> v(3);
    v.clear();
    _atom atom;
    atom = a.atoms[i];
    atom.cpos.clear();
    atom.fpos.clear();
    a.qm_forces.push_back(v);
    a.qm_positions.push_back(v);
    a.qm_atoms.push_back(atom);
  }

  // ---------------------------------------------------------------
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " WRAPPING [13]" << endl;
  }
  // TOLERANCES ------------------------
  a.equiv_fpos_epsilon = _EQUIV_FPOS_EPS_; // standard but you can change
  // SORT ATOMS (false) -----------------------------
  // AFLOW prefers alphabetic ordering; HOWEVER, sorting by default can cause
  // issues for functions/processes outside of AFLOW (e.g., settings in INCAR).
  // For now, it is safer to sort the atoms inside the particular AFLOW function
  // where it is needed (e.g., symmetry and prototype functions).
  // For readers/writers other than the VASP geometry file, we will always sort
  // alphabetically. //DX+CO20210706
  const bool force_alphabetic_sorting = false;
  if (force_alphabetic_sorting) {
    a.SpeciesPutAlphabetic();
    std::stable_sort(a.atoms.begin(), a.atoms.end(), sortAtomsNames);
  }
  // MAKE BASIS
  a.MakeBasis();
  if (force_alphabetic_sorting) {
    a.MakeTypes();
  } // DX+CO20210706
  // FLAGS -----------------------------
  a.Niggli_calculated = false;
  a.Niggli_avoid = false;
  a.Minkowski_calculated = false;
  a.Minkowski_avoid = false;
  a.LatticeReduction_calculated = false;
  a.LatticeReduction_avoid = false;
  // LATTICE FLAGS -----------------------------
  a.Standard_Lattice_calculated = false;
  a.Standard_Lattice_avoid = false;
  a.Standard_Lattice_primitive = false;
  a.Standard_Lattice_conventional = false;
  a.Standard_Lattice_has_failed = false;
  a.bravais_lattice_type = "";
  a.bravais_lattice_variation_type = "";
  a.bravais_lattice_system = "";
  a.bravais_lattice_lattice_type = "";
  a.bravais_lattice_lattice_variation_type = "";
  a.bravais_lattice_lattice_system = "";
  a.pearson_symbol = "";
  a.reciprocal_lattice_type = "";
  a.reciprocal_lattice_variation_type = "";
  a.bravais_superlattice_type = "";
  a.bravais_superlattice_variation_type = "";
  a.bravais_superlattice_system = "";
  a.pearson_symbol_superlattice = "";
  // QM CALCULATED STUFF
  a.qm_origin.clear();
  a.qm_scale = 0.0;
  a.qm_lattice.clear();
  a.qm_klattice.clear();
  a.qm_f2c.clear();
  a.qm_c2f.clear();
  a.qm_calculated = false;
  a.qm_forces_write = false;
  a.qm_positions_write = false;
  // CALCULATED STUFF
  // DX+CO START
  a.dist_nn_min = AUROSTD_NAN; // CO
  a.sym_eps = AUROSTD_NAN; // DX
  a.sym_eps_calculated = false; // DX
  a.sym_eps_change_count = 0; // DX20180222 - added tolerance count specific to structure
  a.sym_eps_no_scan = false; // DX20210331 - added no scan specific to structure
  // DX+CO END
  a.iatoms_calculated = false;
  a.pgroup_calculated = false;
  a.pgroup_xtal_calculated = false;
  a.pgroupk_Patterson_calculated = false; // DX20200129
  a.pgroupk_calculated = false;
  a.pgroupk_xtal_calculated = false; // DX20171205 - Added pgroupk_xtal
  a.fgroup_calculated = false;
  a.sgroup_calculated = false;
  a.grid_atoms_calculated = false;
  a.lijk_calculated = false;
  // DX20180712 START
  //  ANRL SYMBOLIC MATH
  a.symbolic_math_representation_only = false;
  a.constrained_symmetry_calculation = false;
  // DX20180712 END
  //  OUTPUT STUFF
  a.error_flag = false;
  a.error_string = "";
  a.write_lattice_flag = false;
  a.write_inequivalent_flag = false;
  a.write_klattice_flag = false;
  a.write_DEBUG_flag = false;
#ifdef XSTR_DEBUG
  a.write_lattice_flag = true;
  a.write_inequivalent_flag = true;
  a.write_klattice_flag = true;
  a.write_DEBUG_flag = true;
#endif
  // RF20200310 BEGIN
  for (size_t i = 0; i < a.atoms.size(); i++) {
    // CO20200624
    if (a.atoms[i].cleanname.empty()) {
      a.atoms[i].CleanName();
    }
    //(KBIN::VASP_PseudoPotential_CleanName(a.atoms[i].name)); //CO20200624 - fixed CleanName()
  }
  // RF20200310 END
  //  CHECKS
  if (a.atoms.size() != a.qm_atoms.size()) {
    message << "a.atoms.size()!=a.qm_atoms.size() " << endl; // CO20190629
    for (size_t i = 0; i < vinput.size(); i++) {
      message << vinput[i] << endl; // CO20190629
    }
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
  }
  if (a.atoms.size() != a.qm_forces.size()) {
    message << "a.atoms.size()!=a.qm_forces.size() " << endl; // CO20190629
    for (size_t i = 0; i < vinput.size(); i++) {
      message << vinput[i] << endl; // CO20190629
    }
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
  }
  if (a.atoms.size() != a.qm_positions.size()) {
    message << "a.atoms.size()!=a.qm_positions.size() " << endl; // CO20190629
    for (size_t i = 0; i < vinput.size(); i++) {
      message << vinput[i] << endl; // CO20190629
    }
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // CO20190629
  }
  if (det(a.lattice) < 0.0) {
    // CO20200201
    message << "Found negative determinant for lattice (det()=" << det(a.lattice) << "). Flip your basis."; // CO20200201
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_); // CO20200201
  } // CO20200201

  // ---------------------------------------------------------------
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " DONE [99]" << endl;
  }
  // DONE
  return cinput;
}

// **************************************************************************
// Xstructure operator<< OUTPUT_XSTRUCTURE_OUTPUT
// **************************************************************************
// print an xstructure in a variety of forms
ostream& operator<<(ostream& cout, const xstructure& a) {
  // operator<<
  const bool LDEBUG = (false || XHOST.DEBUG);
  stringstream message;
  int a_iomode = a.iomode;
  //  DEBUG=true;
  // ----------------------------------------------------------------------
  // PUT DEFAULT
  if (a_iomode == IOAFLOW_AUTO) {
    a_iomode = IOVASP_AUTO; // some default
  }
  // ----------------------------------------------------------------------
  // VASP OUTPUT
  if (a_iomode == IOVASP_AUTO || a_iomode == IOVASP_POSCAR || a_iomode == IOVASP_ABCCAR || a_iomode == IOVASP_WYCKCAR) {
    // VASP POSCAR
    cout.setf(std::ios::fixed, std::ios::floatfield);
    const uint _precision_ = _DOUBLE_WRITE_PRECISION_MAX_; // 14; //was 16 SC 10 DM //CO20180515
    cout.precision(_precision_);
    // DX20180618 - Check for symbolic representaion only - START
    if (a.symbolic_math_representation_only) {
      xstructure aa(a);
      cout << aa.PrintSymbolicMathRepresentation();
      return cout;
    }
    // DX20180618 - Check for symbolic representation only - END
    if (a_iomode == IOAFLOW_AUTO) {
      cout << a.title << endl; // << " (AUTO) " << endl;
    }
    if (a_iomode == IOVASP_AUTO) {
      cout << a.title << endl; // << " (AUTO) " << endl;
    }
    if (a_iomode == IOVASP_POSCAR) {
      cout << a.title << endl; // << " (POSCAR) " << endl;
    }
    if (a_iomode == IOVASP_ABCCAR) {
      cout << a.title << endl; // << " (ABCCAR) " << endl;
    }
    // DX20221128 - START
    if (a_iomode == IOVASP_WYCKCAR) {
      cout << a.title << "| SG: ";
      // check if space group has been calculated, otherwise GetSpaceGroupName() could fail //DX20221128
      try {
        cout << GetSpaceGroupName(a.space_group_ITC, a.directory) << " " << a.space_group_ITC << " PG: " << a.point_group_ITC << " BL: " << a.bravais_label_ITC << " | sym_eps: " << a.sym_eps << endl;
        // DX20210526 - extend title
      } catch (aurostd::xerror& re) {
        try {
          cout << GetSpaceGroupName(a.spacegroupnumber, a.directory) << " " << a.spacegroupnumber << " PG: " << a.point_group_ITC << " BL: " << a.bravais_label_ITC << " | sym_eps: " << a.sym_eps << endl;
          // DX20210526 - extend title
        } catch (aurostd::xerror& re) {
          if (LDEBUG) {
            message << "Cannot determine space group name" << endl;
          }
        }
      }
    }
    // DX20221128 - STOP
    if (a.neg_scale == false || a_iomode == IOVASP_WYCKCAR) {
      // DX20210708 - wyccar should always use scale factor
      cout.precision(6); // DM
      cout << a.scale; // << endl; //CO20170630
    } else {
      // oss << a.scale << endl;
      const double s = a.scale;
      const double vol = s * s * s * GetVol(a.lattice);
      cout.precision(6); // DM
      cout << -1 * vol; // << endl; //CO20170630
    }
    cout.precision(_precision_); // SC to cut/paste from matlab in format long
    // CO20170630, add pocc tol to make truly pocc readable by aflow
    if (a.partial_occupation_flag == true) {
      cout.unsetf(ios_base::floatfield);
      cout << " "; //<< std::defaultfloat;
      if (a.neg_scale_second) {
        cout << (-1) * a.partial_occupation_HNF;
      } else {
        cout << a.partial_occupation_site_tol;
      }
      if (true || a.scale_third.isentry) {
        // CO20170803 - stoich tol //always print
        cout << " ";
        cout << a.partial_occupation_stoich_tol;
      }
      cout << std::fixed;
    }
    cout << endl;
    // ----------------------------------------------------------------------
    if (a_iomode == IOVASP_POSCAR || a_iomode == IOVASP_AUTO) {
      for (uint i = 1; i <= 3; i++) {
        for (uint j = 1; j <= 3; j++) {
          cout << " ";
          if (std::abs(a.lattice(i, j)) < 10.0) {
            cout << " ";
          }
          if (!std::signbit(a.lattice(i, j))) {
            cout << " ";
          }
          cout << a.lattice(i, j) << "";
        }
        cout << endl;
      }
    }
    // ----------------------------------------------------------------------
    if (a_iomode == IOVASP_ABCCAR) {
      // DX20210525 - separated a_iomode==IOVASP_WYCKCAR
      cout << " ";
      cout.precision(10); // SC to cut/paste from matlab in format long
      if (std::abs(a.a) < 10.0) {
        cout << " ";
      }
      if (!std::signbit(a.a)) {
        cout << " ";
      }
      cout << a.a << "";
      if (std::abs(a.b) < 10.0) {
        cout << " ";
      }
      if (!std::signbit(a.b)) {
        cout << " ";
      }
      cout << a.b << "";
      if (std::abs(a.c) < 10.0) {
        cout << " ";
      }
      if (!std::signbit(a.c)) {
        cout << " ";
      }
      cout << a.c << "";
      cout.precision(4); // SC to cut/paste from matlab in format long
      if (std::abs(a.alpha) < 10.0) {
        cout << " ";
      }
      if (!std::signbit(a.alpha)) {
        cout << " ";
      }
      cout << a.alpha << "";
      if (std::abs(a.beta) < 10.0) {
        cout << " ";
      }
      if (!std::signbit(a.beta)) {
        cout << " ";
      }
      cout << a.beta << "";
      if (std::abs(a.gamma) < 10.0) {
        cout << " ";
      }
      if (!std::signbit(a.gamma)) {
        cout << " ";
      }
      cout << a.gamma << "";
      cout << endl;
      cout.precision(_precision_); // SC to cut/paste from matlab in format long
    } else if (a_iomode == IOVASP_WYCKCAR) {
      // DX20210525 - note wyccar uses lattice parameters of the conventional cell
      const xvector<double> data = Getabc_angles(a.standard_lattice_ITC, DEGREES);
      cout << " ";
      cout.precision(10); // SC to cut/paste from matlab in format long
      if (std::abs(data(1)) < 10.0) {
        cout << " ";
      }
      if (!std::signbit(data(1))) {
        cout << " ";
      }
      cout << data(1) << "";
      if (std::abs(data(2)) < 10.0) {
        cout << " ";
      }
      if (!std::signbit(data(2))) {
        cout << " ";
      }
      cout << data(2) << "";
      if (std::abs(data(3)) < 10.0) {
        cout << " ";
      }
      if (!std::signbit(data(3))) {
        cout << " ";
      }
      cout << data(3) << "";
      cout.precision(4); // SC to cut/paste from matlab in format long
      if (std::abs(data(4)) < 10.0) {
        cout << " ";
      }
      if (!std::signbit(data(4))) {
        cout << " ";
      }
      cout << data(4) << "";
      if (std::abs(data(5)) < 10.0) {
        cout << " ";
      }
      if (!std::signbit(data(5))) {
        cout << " ";
      }
      cout << data(5) << "";
      if (std::abs(a.gamma) < 10.0) {
        cout << " ";
      }
      if (!std::signbit(data(6))) {
        cout << " ";
      }
      cout << data(6) << "";
      cout << " " << a.space_group_ITC << " " << a.setting_ITC;
      cout << endl;
      cout.precision(_precision_); // SC to cut/paste from matlab in format long
    }
    // ----------------------------------------------------------------------
    if (a.is_vasp4_poscar_format == true) {} // nothing to do

    if (a.is_vasp5_poscar_format == true) {
      for (size_t i = 0; i < a.species_pp.size(); i++) {
        // ME20190308 - species is empty when structure is based on vasp4 POSCAR
        cout << KBIN::VASP_PseudoPotential_CleanName(a.species_pp[i]) << " "; // ME20190308
      }
      cout << endl;
    }
    // DX20210526 - add WYCCAR format - START
    if (a_iomode == IOVASP_WYCKCAR) {
      cout << aurostd::joinWDelimiter(SYM::countWyckoffTypes(a.wyckoff_sites_ITC), " ") << endl;
      // gets the number of Wyckoff positions per type
      cout << "Direct(WYCCAR)" << endl; // wyccar is always in direct/fractional
      const size_t nWyckoff_sites = a.wyckoff_sites_ITC.size();
      double _coord = AUROSTD_MAX_DOUBLE;
      for (uint i = 0; i < nWyckoff_sites; i++) {
        _coord = AUROSTD_MAX_DOUBLE;
        cout << " ";
        for (uint j = 1; j <= 3; j++) {
          _coord = aurostd::roundoff(a.wyckoff_sites_ITC[i].coord(j), pow(10.0, -(double) _precision_));
          if (std::abs(_coord) < 10.0) {
            cout << " ";
          }
          if (!std::signbit(_coord)) {
            cout << " ";
          }
          cout << _coord << " ";
        }
        cout << std::setprecision(_precision_) << std::left << " " << a.wyckoff_sites_ITC[i].type << " " << a.wyckoff_sites_ITC[i].multiplicity << " " << a.wyckoff_sites_ITC[i].letter << " "
             << a.wyckoff_sites_ITC[i].site_symmetry;
        cout << endl;
      }
    }
    // DX20210526 - add WYCCAR format - END
    //  ----------------------------------------------------------------------
    // CO20170630 - fixing for POCC
    //[CO20180705 - we have const str&, so we can't modify atom arrangement, this MUST be done before structure is printed]a.MakeTypes();  //CO20180705 - repetita iuvant
    //[CO20180705 - we have const str&, so we can't modify atom arrangement, this MUST be done before structure is printed]std::stable_sort(a.atoms.begin(),a.atoms.end(),sortAtomsType);  //CO20180705 - this makes it necessary that atoms are properly typed
    if (a_iomode != IOVASP_WYCKCAR) {
      // DX20210611 - do not do for Wyccar, atom count is not the same as number of Wyckoff positiosn
      if (a.partial_occupation_flag == true) {
        // need to figure out the '+'
        uint iatom = 0;
        vector<vector<uint>> vsame_pocc;
        double last_pocc = 0.0;
        if (LDEBUG) {
          for (size_t i = 0; i < a.atoms.size(); i++) {
            cerr << __AFLOW_FUNC__ << " name=" << a.atoms[i].name << " type=" << a.atoms[i].type << " pocc=" << a.atoms[i].partial_occupation_value << endl;
          }
        }
        if (!a.atoms.empty()) {
          for (size_t i = 0; i < a.num_each_type.size(); i++) {
            vsame_pocc.emplace_back(0); // for first atom
            vsame_pocc.back().push_back(0); // first atom
            last_pocc = a.atoms[iatom].partial_occupation_value;
            for (uint j = 0; j < (uint) a.num_each_type[i]; j++) {
              if (aurostd::isequal(a.atoms[iatom].partial_occupation_value, last_pocc, _AFLOW_POCC_ZERO_TOL_)) {
                vsame_pocc.back().back() += 1;
              } // same
              else {
                //'+'
                vsame_pocc.back().push_back(1);
                last_pocc = a.atoms[iatom].partial_occupation_value;
              }
              iatom++;
            }
          }
        }
        if (LDEBUG) {
          for (size_t i = 0; i < vsame_pocc.size(); i++) {
            for (size_t j = 0; j < vsame_pocc[i].size(); j++) {
              cerr << __AFLOW_FUNC__ << " vsame_pocc[" << i << "][" << j << "]=" << vsame_pocc[i][j] << endl;
            }
          }
        }
        iatom = 0;
        for (size_t i = 0; i < vsame_pocc.size(); i++) {
          // if(vsame_pocc[i].size()==1){  //no '+'
          //   oss << vsame_pocc[i][0] << "*";
          //   //oss << std::defaultfloat;
          //   oss.unsetf(ios_base::floatfield);
          //   oss << a.atoms[iatom++].partial_occupation_value << std::fixed << " ";
          // } else {  //need '+'
          for (size_t j = 0; j < vsame_pocc[i].size(); j++) {
            cout << vsame_pocc[i][j] << "*";
            // oss << std::defaultfloat;
            cout.unsetf(ios_base::floatfield);
            cout << a.atoms[iatom].partial_occupation_value << std::fixed << (j != vsame_pocc[i].size() - 1 ? "+" : " ");
            iatom += vsame_pocc[i][j];
          }
          //}
        }
      } else {
        for (size_t i = 0; i < a.num_each_type.size(); i++) {
          cout << a.num_each_type[i] << " ";
        }
      }
      cout << endl;
      if (a.isd) {
        cout << "Selective Dynamics" << endl; // DONE YOYO BUG
      }
      // oss << a.coord_type << endl;

      if (a.coord_flag == _COORDS_FRACTIONAL_) {
        cout << "Direct(" << a.atoms.size() << ") ";
      }
      if (a.coord_flag == _COORDS_CARTESIAN_) {
        cout << "Cartesian(" << a.atoms.size() << ") ";
      }
      //  if(a.partial_occupation_flag==true)  oss << "Pocc ";
      if (a.order_parameter_structure == true) {
        cout << "OrderParameter(" << a.order_parameter_atoms.size() << ") ";
      }
      //      oss << "[";for(size_t i=0;i<a.num_each_type.size();i++) {oss << char('A'+i) << a.num_each_type.at(i);}oss << "] ";
      // CO20170630, the original num_each_type doesn't work here, so we fix
      if (a.partial_occupation_flag == true) {
        cout << "Partial ";
        // oss.precision(_pocc_precision_);  //CO20170630
        // oss << std::defaultfloat;
        const int comp_prec = (int) ceil(log10(1.0 / a.partial_occupation_stoich_tol));
        // ceil ensures we round up above 1 //CO20181226
        cout.precision(comp_prec); // CO20181226
        cout.unsetf(ios_base::floatfield);
        cout << "[";
        for (size_t i = 0; i < a.comp_each_type.size(); i++) {
          cout << char('A' + i) << a.comp_each_type[i];
        }
        cout << "] ";
        cout << std::fixed;
        cout.precision(_precision_); // CO20170630 //CO20181226
      } else {
        cout << "[";
        for (size_t i = 0, k = 0; i < a.num_each_type.size(); k += a.num_each_type.at(i), i++) {
          cout << char(a.atoms.at(k).type + 65) << a.num_each_type.at(i);
        }
        cout << "] ";
      }

      // done
      cout << endl;

      double _coord; // CO20190322 - remove annoying -0.0000000
      for (size_t iat = 0; iat < a.atoms.size(); iat++) {
        cout << " ";
        for (uint j = 1; j <= 3; j++) {
          //	oss << " ";
          if (a.coord_flag == _COORDS_FRACTIONAL_) {
            _coord = a.atoms[iat].fpos(j);
          }
          // CO20190322 - remove annoying -0.0000000
          if (a.coord_flag == _COORDS_CARTESIAN_) {
            _coord = a.atoms[iat].cpos(j);
          }
          // CO20190322 - remove annoying -0.0000000

          _coord = aurostd::roundoff(_coord, pow(10.0, -(double) _precision_)); // CO20190322 - remove annoying -0.0000000
          if (std::abs(_coord) < 10.0) {
            cout << " "; // CO20190322 - remove annoying -0.0000000
          }
          if (!std::signbit(_coord)) {
            cout << " "; // CO20190322 - remove annoying -0.0000000
          }
          cout << _coord << " "; // CO20190322 - remove annoying -0.0000000
        }
        //  cout << aurostd::modulus(a.atoms[iat].cpos) << " ";
        if (a.isd == true) {
          cout << " " << a.atoms[iat].sd[0] << " " << a.atoms[iat].sd[1] << " " << a.atoms[iat].sd[2];
        }
        if (a.atoms[iat].name_is_given == true) {
          cout << " " << a.atoms[iat].name << " ";
          for (uint j = a.atoms[iat].name.length(); j < 5; j++) {
            cout << " ";
          }
        }
        if (a.partial_occupation_flag == true) {
          // oss.precision(_pocc_precision_);  //CO20170630
          // if(a.atoms[iat].partial_occupation_flag==false) oss << "-      ";
          //	if(a.atoms[iat].partial_occupation_flag==true) oss << a.atoms[iat].partial_occupation_value << "  ";// << " (" << iat << "/" << a.partial_occupation_flags.size() << ")";
          // oss << std::defaultfloat;
          cout.unsetf(ios_base::floatfield);
          cout << "pocc=" << a.atoms[iat].partial_occupation_value << "  ";
          cout << std::fixed;
          // oss.precision(_precision_); //CO20170630
        }
        if (a.order_parameter_structure == true) {
          if (a.atoms[iat].order_parameter_atom == false) {
            cout << "- ";
          }
          if (a.atoms[iat].order_parameter_atom == true) {
            cout << a.atoms[iat].order_parameter_value << " ";
          }
          // << " (" << iat << "/" << a.order_parameter_atoms.size() << ")";
        }
        if (a.write_inequivalent_flag == true) {
          cout << " ";
          // ?	if(i<10) oss << "0";
          cout << iat << "[";
          if (a.atoms[iat].equivalent < 10) {
            cout << "0";
          }
          cout << a.atoms[iat].equivalent << "]";
          if (a.atoms[iat].is_inequivalent) {
            cout << "*";
            cout << "_(" << a.atoms[iat].num_equivalents << ") ";
            //<< "  index=" << a.atoms[iat].index_iatoms << "  ";
            //  " v" << a.iatoms.size() << "   burp ";
            // for(size_t jat=0;jat<a.iatoms.size();jat++)  oss << a.iatoms.at(jat).size() << " ";
          }
        }
        if (a.qm_forces_write) {
          if (a.qm_calculated == true) {
            cout << "F *(";
          }
          if (a.qm_calculated == false) {
            cout << "F  (";
          }
          for (uint j = 1; j <= 3; j++) {
            if (std::abs(a.qm_forces.at(iat)(j)) < 10.0) {
              cout << " ";
            }
            if (!std::signbit(a.qm_forces.at(iat)(j))) {
              cout << " ";
            }
            cout << a.qm_forces.at(iat)(j) << " ";
          }
          cout << ")_   ";
        }
        if (a.qm_positions_write) {
          if (a.qm_calculated == true) {
            cout << "P *(";
          }
          if (a.qm_calculated == false) {
            cout << "P  (";
          }
          for (uint j = 1; j <= 3; j++) {
            if (std::abs(a.qm_positions.at(iat)(j)) < 10.0) {
              cout << " ";
            }
            if (!std::signbit(a.qm_positions.at(iat)(j))) {
              cout << " ";
            }
            cout << a.qm_positions.at(iat)(j) << " ";
          }
          cout << ")_   ";
        }
        if (a.write_DEBUG_flag) {
          cout << " s" << a.atoms[iat].spin;
          //[CO20200130 - number->basis]oss << " n"<<a.atoms[iat].number;
          cout << " b" << a.atoms[iat].basis;
          cout << " N(" << a.atoms[iat].cleanname;
          cout << " " << a.atoms[iat].atomic_number << " " << " [" << a.atoms[iat].type << "] ";
          cout << " ijk(" << a.atoms[iat].ijk(1) << "," << a.atoms[iat].ijk(2) << "," << a.atoms[iat].ijk(3) << ")";
        }
        cout << endl;
        cout.flush();
      } // iat
    } // DX20210610 - end Wyccar if-statement
    if (a.write_lattice_flag) {
      cout << "DIRECT LATTICE per raw" << endl;
      for (uint i = 1; i <= 3; i++) {
        for (uint j = 1; j <= 3; j++) {
          cout << " ";
          if (!std::signbit(a.scale * a.lattice(i, j))) {
            cout << " ";
          }
          cout << a.scale * a.lattice(i, j) << " ";
        }
        cout << endl;
      }
    }
    if (a.write_klattice_flag) {
      cout << "RECIPROCAL LATTICE per raw" << endl;
      for (uint i = 1; i <= 3; i++) {
        for (uint j = 1; j <= 3; j++) {
          cout << " ";
          if (!std::signbit(a.klattice(i, j))) {
            cout << " ";
          }
          cout << a.klattice(i, j) << " ";
        }
        cout << endl;
      }
    }
    if (a.write_lattice_flag && a.write_klattice_flag) {
      cout << "ORTOGONALITY (a*b')/2pi=I" << endl;
      xmatrix<double> orto(3, 3);
      orto = ((a.scale * a.lattice)) * trasp(a.klattice) / (2.0 * pi);
      for (uint i = 1; i <= 3; i++) {
        for (uint j = 1; j <= 3; j++) {
          cout << " ";
          if (!std::signbit(orto(i, j))) {
            cout << " ";
          }
          cout << orto(i, j) << " ";
        }
        cout << endl;
      }
    }
    if (a.write_DEBUG_flag) {
      cout << "kpoints_k1,k2,k3 = " << a.kpoints_k1 << "," << a.kpoints_k2 << "," << a.kpoints_k3 << endl;
      cout << "kpoints_kmax     = " << a.kpoints_kmax << endl;
      cout << "kpoints_kppra    = " << a.kpoints_kppra << endl;
      cout << "kpoints_mode     = " << a.kpoints_mode << endl;
      cout << "kpoints_kscheme  = " << a.kpoints_kscheme << endl;
    }
    if (false && a.partial_occupation_flag == true) {
      // KY
      cout << "*******************************" << endl;
      for (int hnf = 1; hnf < 10; hnf++) {
        cout << hnf << " ";
        cout << endl;
      }
    }
    // DX20180618 - Check for symmetry constrained calculation - START
    if (a.constrained_symmetry_calculation) {
      xstructure aa(a);
      cout << aa.PrintSymbolicMathRepresentation();
    }
    // DX20180618 - Check for symmetry constrained calculation - END
    return cout;
  } // END OF VASP
  // ----------------------------------------------------------------------
  //  QUANTUM ESPRESSO OUTPUT
  if (a_iomode == IOQE_AUTO || a_iomode == IOQE_GEOM) {
    // VASP POSCAR
    cout << "! AFLOW::QE BEGIN " << endl;
    const uint depthQE = 27;
    if (a_iomode == IOQE_AUTO) {
      cout << "! " << a.title << endl; //<< " (AUTO)" <<endl;
    }
    if (a_iomode == IOQE_GEOM) {
      cout << "! " << a.title << endl; //<< " (GEOM)" <<endl;
    }
    cout << aurostd::PaddedPOST("&system", depthQE, " ") << " ! // aflow " << endl;
    cout << aurostd::PaddedPOST(" ibrav=0,", depthQE, " ") << " ! // free " << endl;
    cout << aurostd::PaddedPOST(" nat=" + aurostd::utype2string(a.atoms.size()) + ",", depthQE) << " ! // a.atoms.size() " << endl;
    cout << aurostd::PaddedPOST(" ntyp=" + aurostd::utype2string(a.num_each_type.size()), depthQE) << " ! // a.num_each_type.size() " << endl;
    cout << " /" << endl;
    cout.setf(std::ios::fixed, std::ios::floatfield);
    const uint _precision_ = _DOUBLE_WRITE_PRECISION_MAX_; // 14; //was 16 SC 10 DM //CO20180515
    cout.precision(_precision_);
    if (a.coord_flag == _COORDS_FRACTIONAL_) {
      cout << "ATOMIC_POSITIONS (crystal)" << endl;
    }
    if (a.coord_flag == _COORDS_CARTESIAN_) {
      cout << "ATOMIC_POSITIONS (angstrom)" << endl;
    }
    for (size_t iat = 0; iat < a.atoms.size(); iat++) {
      cout << " ";
      if (a.atoms[iat].name_is_given == true) {
        cout << " " << aurostd::PaddedPOST(KBIN::VASP_PseudoPotential_CleanName(a.atoms[iat].name), 5, " ") << " ";
      } else {
        message << "QE needs atoms species names"; // CO20190629
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_MISSING_); // CO20190629
      }
      for (uint j = 1; j <= 3; j++) {
        //  oss << " ";
        if (a.coord_flag == _COORDS_FRACTIONAL_) {
          if (std::abs(a.atoms[iat].fpos(j)) < 10.0) {
            cout << " ";
          }
          if (!std::signbit(a.atoms[iat].fpos(j))) {
            cout << " ";
          }
          cout << a.atoms[iat].fpos(j) << " ";
        }
        if (a.coord_flag == _COORDS_CARTESIAN_) {
          if (std::abs(a.atoms[iat].cpos(j)) < 10.0) {
            cout << " ";
          }
          if (!std::signbit(a.atoms[iat].cpos(j))) {
            cout << " ";
          }
          cout << a.atoms[iat].cpos(j) << " ";
        }
      }
      cout << " ! // " << a.atoms[iat].cleanname << " ";
      if (a.write_inequivalent_flag == true) {
        cout << " ";
        //	if(i<10) oss << "0";
        cout << iat << "[";
        if (a.atoms[iat].equivalent < 10) {
          cout << "0";
        }
        cout << a.atoms[iat].equivalent << "]";
        if (a.atoms[iat].is_inequivalent) {
          cout << "*";
          cout << "_(" << a.atoms[iat].num_equivalents << ") ";
          //<< "  index=" << a.atoms[iat].index_iatoms << "  ";
          //  " v" << a.iatoms.size() << "   burp ";
          // for(size_t jat=0;jat<a.iatoms.size();jat++)  oss << a.iatoms.at(jat).size() << " ";
        }
      }
      cout << endl;
    }
    // ----------------------------------------------------------------------
    cout.precision(_precision_); // SC to cut/paste from matlab in format long
    cout << "CELL_PARAMETERS (angstrom)" << endl;
    {
      for (uint i = 1; i <= 3; i++) {
        for (uint j = 1; j <= 3; j++) {
          cout << " ";
          if (std::abs(a.lattice(i, j)) < 10.0) {
            cout << " ";
          }
          if (!std::signbit(a.lattice(i, j))) {
            cout << " ";
          }
          cout << a.lattice(i, j) * a.scale << ""; // DX20180215 - added scaling factor
        }
        cout << endl;
      }
    }
    cout << "! AFLOW::QE END " << endl; // AZ20240223 changed from # to ! so quantum espresso output works
    return cout;
  }
  // ----------------------------------------------------------------------
  //  CIF OUTPUT
  if (a.iomode == IOCIF) {
    // CIF
    pflow::PrintCIF(cout, a, a.spacegroupnumber, a.setting_ITC);
    // DX20210630 - add setting (otherwise, this will mess up rhl systems by mixing the hex and rhl setting)
    return cout;
  }

  // ----------------------------------------------------------------------
  //  ABINIT OUTPUT
  if (a_iomode == IOABINIT_AUTO || a_iomode == IOABINIT_GEOM) {
    // VASP POSCAR
    xstructure aa(a); // DX20210415 - need to make a copy to rescale
    aa.ReScale(1.0); // DX20210415
    cout << "# AFLOW::ABINIT BEGIN " << endl;
    const uint _precision_ = _DOUBLE_WRITE_PRECISION_MAX_; // 14; //was 16 SC 10 DM //CO20180515
    cout.precision(_precision_);
    cout.setf(std::ios::fixed, std::ios::floatfield);
    if (a_iomode == IOABINIT_AUTO) {
      cout << "# " << aa.title << endl; //<< " (AUTO)" << endl;
    }
    if (a_iomode == IOABINIT_GEOM) {
      cout << "# " << aa.title << endl; //<< " (GEOM)" << endl;
    }
    cout << "acell   " << double(1) << "   " << double(1) << "   " << double(1) << "  ANGSTR" << endl;
    // scaling of the primitive vectors, in Bohr.
    for (uint j = 1; j <= 3; j++) {
      // CO20190908 - manual is misleading, it's row-based// each COLUMN of this array is one primitive translation
      if (j == 1) {
        cout << "rprim";
      }
      if (j == 2) {
        cout << "     ";
      }
      if (j == 3) {
        cout << "     ";
      }
      for (uint i = 1; i <= 3; i++) {
        cout << " ";
        if (std::abs(aa.lattice(j, i)) < 10.0) {
          cout << " "; // CO20190908 - manual is misleading, it's row-based
        }
        if (!std::signbit(aa.lattice(j, i))) {
          cout << " "; // CO20190908 - manual is misleading, it's row-based
        }
        cout << aa.lattice(j, i) << ""; // CO20190908 - manual is misleading, it's row-based
      }
      cout << endl;
    }
    cout << "natom " << aa.atoms.size() << endl;
    // DX20200313 - add atom type info via znucl - START
    cout << "znucl ";
    for (size_t i = 0; i < aa.species.size(); i++) {
      for (size_t e = 0; e < velement.size(); e++) {
        // external variable (see aflow_xelement.h)
        if (velement[e].symbol == KBIN::VASP_PseudoPotential_CleanName(aa.species[i])) {
          cout << e << " "; // index corresponds to Z value
          break;
        }
      }
    }
    cout << endl;
    // DX20200313 - add atom type info via znucl - END
    cout << "typat ";
    //   for(size_t i=0;i<aa.num_each_type.size();i++) oss << a.num_each_type.at(i) << " ";  oss << endl;
    for (size_t i = 0; i < aa.atoms.size(); i++) {
      cout << aa.atoms[i].type + 1 << " ";
    }
    cout << endl;
    if (aa.coord_flag == _COORDS_FRACTIONAL_) {
      cout << "xred " << endl;
    }
    if (aa.coord_flag == _COORDS_CARTESIAN_) {
      cout << "xangst " << endl;
    }
    for (size_t iat = 0; iat < aa.atoms.size(); iat++) {
      cout << "      ";
      for (uint j = 1; j <= 3; j++) {
        if (aa.coord_flag == _COORDS_FRACTIONAL_) {
          if (std::abs(aa.atoms[iat].fpos(j)) < 10.0) {
            cout << " ";
          }
          if (!std::signbit(aa.atoms[iat].fpos(j))) {
            cout << " ";
          }
          cout << aa.atoms[iat].fpos(j) << " ";
        }
        if (aa.coord_flag == _COORDS_CARTESIAN_) {
          if (std::abs(aa.atoms[iat].cpos(j)) < 10.0) {
            cout << " ";
          }
          if (!std::signbit(aa.atoms[iat].cpos(j))) {
            cout << " ";
          }
          cout << aa.atoms[iat].cpos(j) << " ";
        }
      }
      cout << " # " << aa.atoms[iat].cleanname << " ";
      if (aa.write_inequivalent_flag == true) {
        cout << " ";
        //	if(i<10) oss << "0";
        cout << iat << "[";
        if (aa.atoms[iat].equivalent < 10) {
          cout << "0";
        }
        cout << aa.atoms[iat].equivalent << "]";
        if (aa.atoms[iat].is_inequivalent) {
          cout << "*";
          cout << "_(" << aa.atoms[iat].num_equivalents << ") ";
          //<< "  index=" << a.atoms.at(iat).index_iatoms << "  ";
          //  " v" << aa.iatoms.size() << "   burp ";
          // for(size_t jat=0;jat<aa.iatoms.size();jat++)  oss << a.iatoms.at(jat).size() << " ";
        }
      }
      cout << endl;
    }
    cout << "# AFLOW::ABINIT END " << endl;
    return cout;
  }

  // ----------------------------------------------------------------------
  //  ELK OUTPUT //DX20200315
  if (a_iomode == IOELK_AUTO || a_iomode == IOELK_GEOM) {
    // ELK
    xstructure aa(a); // DX20210415 - need to make a copy to rescale
    aa.ReScale(1.0); // DX20210415
    cout << "# AFLOW::ELK BEGIN " << endl;
    const uint _precision_ = _DOUBLE_WRITE_PRECISION_MAX_; // 14; //was 16 SC 10 DM //CO20180515
    cout.precision(_precision_);
    cout.setf(std::ios::fixed, std::ios::floatfield);
    cout << "# " << aa.title << endl;
    cout << endl;
    // scaling factors
    cout << "scale" << endl << " " << aa.scale << endl << endl;
    cout << "scale1" << endl << " 1.0" << endl << endl; // returns unscaled (for now)
    cout << "scale2" << endl << " 1.0" << endl << endl; // returns unscaled (for now)
    cout << "scale3" << endl << " 1.0" << endl << endl; // returns unscaled (for now)

    // lattice, note: convert to atomic units (Bohr)
    cout << "avec" << endl;
    cout << " " << aa.lattice(1) * angstrom2bohr << endl;
    cout << " " << aa.lattice(2) * angstrom2bohr << endl;
    cout << " " << aa.lattice(3) * angstrom2bohr << endl;
    cout << endl;

    // atom info
    cout << "atoms" << endl;
    cout << " " << setw(49) << std::left << aa.species.size();
    cout << ": nspecies" << endl;
    for (size_t i = 0; i < aa.num_each_type.size(); i++) {
      cout << setw(50) << std::left << "\'" + aa.species[i] + ".in\'";
      cout << ": spfname" << endl;
      cout << " " << setw(49) << std::left << aa.num_each_type[i];
      cout << ": natoms; atpos, bfcmt below" << endl;
      for (size_t iat = 0; iat < aa.atoms.size(); iat++) {
        if (aa.atoms[iat].name == aa.species[i]) {
          // atom coordinates
          for (uint j = 1; j <= 3; j++) {
            if (aa.coord_flag == _COORDS_CARTESIAN_) {
              if (std::abs(aa.atoms[iat].fpos(j)) < 10.0) {
                cout << " ";
              }
              if (!std::signbit(aa.atoms[iat].fpos(j))) {
                cout << " ";
              }
              cout << aa.atoms[iat].fpos(j) << " ";
            }
            if (aa.coord_flag == _COORDS_FRACTIONAL_) {
              if (std::abs(aa.atoms[iat].fpos(j)) < 10.0) {
                cout << " ";
              }
              if (!std::signbit(aa.atoms[iat].fpos(j))) {
                cout << " ";
              }
              cout << aa.atoms[iat].fpos(j) << " ";
            }
          }
          // magnetic field //DX20210409 - updated with non-collinear spin
          for (uint j = 1; j <= 3; j++) {
            if (std::abs(aa.atoms[iat].noncoll_spin(j)) < 10.0) {
              cout << " ";
            }
            if (!std::signbit(aa.atoms[iat].noncoll_spin(j))) {
              cout << " ";
            }
            cout << aa.atoms[iat].noncoll_spin(j) << " ";
          }
          cout << endl;
        }
      }
    }
    cout << "# AFLOW::ELK END " << endl;
    return cout;
  }

  // ----------------------------------------------------------------------
  //  AIMS OUTPUT
  if (a_iomode == IOAIMS_AUTO || a_iomode == IOAIMS_GEOM) {
    // VASP POSCAR
    xstructure aa(a);
    aa.ReScale(1.0); // very important because there is NO scale factor in AIMS //CO20180420
    const uint _precision_ = _DOUBLE_WRITE_PRECISION_MAX_; // 14; //was 16 SC 10 DM //CO20180515
    cout.precision(_precision_);
    cout.setf(std::ios::fixed, std::ios::floatfield);
    cout << "# " << aa.title << endl; //<< " (AUTO)" << endl;
    cout << "# AFLOW::AIMS BEGIN " << endl; // come after title
    // DX20180618 - Check for symbolic representaion only - START
    if (aa.symbolic_math_representation_only) {
      cout << aa.PrintSymbolicMathRepresentation();
      cout << "# AFLOW::AIMS END " << endl;
      return cout;
    }
    // DX20180618 - Check for symbolic representation only - END
    for (uint i = 1; i <= 3; i++) {
      // each COLUMN of this array is one primitive translation
      cout << "lattice_vector ";
      for (uint j = 1; j <= 3; j++) {
        if (std::abs(aa.lattice(i, j)) < 100.0) {
          cout << " ";
        }
        if (std::abs(aa.lattice(i, j)) < 10.0) {
          cout << " ";
        }
        if (!std::signbit(aa.lattice(i, j))) {
          cout << " ";
        }
        cout << aa.lattice(i, j) << "";
      }
      cout << endl;
    }
    for (size_t iat = 0; iat < aa.atoms.size(); iat++) {
      if (false || LDEBUG) {
        cerr << "XSTRUCTURE<<: AIMS aa.coord_flag==" << aa.coord_flag << endl;
        cerr << "XSTRUCTURE<<: AIMS _COORDS_FRACTIONAL_==" << _COORDS_FRACTIONAL_ << endl;
        cerr << "XSTRUCTURE<<: AIMS _COORDS_CARTESIAN_==" << _COORDS_CARTESIAN_ << endl;
      }
      cout << (aa.coord_flag == _COORDS_FRACTIONAL_ ? "atom_frac" : "atom") << " ";
      for (uint j = 1; j <= 3; j++) {
        if (aa.coord_flag == _COORDS_FRACTIONAL_) {
          if (std::abs(aa.atoms[iat].fpos(j)) < 100.0) {
            cout << " ";
          }
          if (std::abs(aa.atoms[iat].fpos(j)) < 10.0) {
            cout << " ";
          }
          if (!std::signbit(aa.atoms[iat].fpos(j))) {
            cout << " ";
          }
          cout << aa.atoms[iat].fpos(j) << " ";
        }
        if (aa.coord_flag == _COORDS_CARTESIAN_) {
          if (std::abs(aa.atoms[iat].cpos(j)) < 100.0) {
            cout << " ";
          }
          if (std::abs(aa.atoms[iat].cpos(j)) < 10.0) {
            cout << " ";
          }
          if (!std::signbit(aa.atoms[iat].cpos(j))) {
            cout << " ";
          }
          cout << aa.atoms[iat].cpos(j) << " ";
        }
      }
      cout << " " << aa.atoms[iat].cleanname << " ";
      if (aa.write_inequivalent_flag == true) {
        cout << " # ";
        //	if(i<10) oss << "0";
        cout << iat << "[";
        if (aa.atoms[iat].equivalent < 10) {
          cout << "0";
        }
        cout << aa.atoms[iat].equivalent << "]";
        if (aa.atoms[iat].is_inequivalent) {
          cout << "*";
          cout << "_(" << aa.atoms[iat].num_equivalents << ") ";
          //<< "  index=" << aa.atoms[iat].index_iatoms << "  ";
          //  " v" << aa.iatoms.size() << "   burp ";
          // for(size_t jat=0;jat<aa.iatoms.size();jat++)  oss << aa.iatoms.at(jat).size() << " ";
        }
      }
      cout << endl;
    }
    // DX20180618 - Check for symmetry constrained calculation - START
    if (aa.constrained_symmetry_calculation) {
      cout << aa.PrintSymbolicMathRepresentation();
    }
    // DX20180618 - Check for symmetry constrained calculation - END
    cout << "# AFLOW::AIMS END " << endl;
    return cout;
  }

  // ----------------------------------------------------------------------
  //  ATAT OUTPUT // SD20220123
  //  Alloy-Theoretic Automated Toolkit
  //  See: https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/manual.pdf
  if (a_iomode == IOATAT_STR) {
    // ATAT
    xstructure aa(a);
    const uint _precision_ = _DOUBLE_WRITE_PRECISION_MAX_; // 14; //was 16 SC 10 DM //CO20180515
    cout.precision(_precision_);
    cout.setf(std::ios::fixed, std::ios::floatfield);
    const xmatrix<double> axes = aurostd::eye<double>(3, 3); // set axes to idenity
    // write the axes
    for (uint i = 1; i <= 3; i++) {
      for (uint j = 1; j <= 3; j++) {
        cout << axes(i, j) << " ";
      }
      cout << endl;
    }
    // write the fractional cell vectors (== lattice)
    for (uint i = 1; i <= 3; i++) {
      for (uint j = 1; j <= 3; j++) {
        cout << " ";
        if (std::abs(aa.lattice(i, j)) < 10.0) {
          cout << " ";
        }
        if (!std::signbit(aa.lattice(i, j))) {
          cout << " ";
        }
        cout << aa.lattice(i, j) << "";
      }
      cout << endl;
    }
    // write the atoms
    for (size_t iat = 0; iat < aa.atoms.size(); iat++) {
      cout << " ";
      for (uint i = 1; i <= 3; i++) {
        if (std::abs(aa.atoms[iat].cpos(i)) < 10.0) {
          cout << " ";
        }
        if (!std::signbit(aa.atoms[iat].cpos(i))) {
          cout << " ";
        }
        cout << aa.atoms[iat].cpos(i) << " ";
      }
      if (aa.atoms[iat].name_is_given == true) {
        cout << aa.atoms[iat].cleanname;
      }
      cout << endl;
    }
    return cout;
  }

  // ----------------------------------------------------------------------
  //  LMP OUTPUT // SD20240111
  if (a_iomode == IOLMP_DATA) {
    // LMP
    xstructure aa(a);
    const uint _precision_ = _DOUBLE_WRITE_PRECISION_MAX_; // 14; //was 16 SC 10 DM //CO20180515
    cout.precision(_precision_);
    cout.setf(std::ios::fixed, std::ios::floatfield);
    cout << "LAMMPS data file via write_data, timestep = 0, units = real" << endl << endl;
    cout << aa.atoms.size() << " atoms" << endl;
    cout << aa.species.size() << " atom types" << endl << endl;
    const double lx = aa.a;
    const double xy = aa.b * std::cos(deg2rad * aa.gamma);
    const double xz = aa.c * std::cos(deg2rad * aa.beta);
    const double ly = std::sqrt(std::pow(aa.b, 2.0) - std::pow(xy, 2.0));
    const double yz = (aa.b * aa.c * std::cos(deg2rad * aa.alpha) - xy * xz) / ly;
    const double lz = std::sqrt(std::pow(aa.c, 2.0) - std::pow(xz, 2.0) - std::pow(yz, 2.0));
    // write the cell
    cout << 0.0 << " " << lx << " xlo xhi" << endl;
    cout << 0.0 << " " << ly << " ylo yhi" << endl;
    cout << 0.0 << " " << lz << " zlo zhi" << endl;
    cout << xy << " " << xz << " " << yz << " xy xz yz" << endl << endl;
    // write the masses
    cout << "Masses" << endl << endl;
    for (size_t i = 0; i < aa.species_mass.size(); i++) {
      cout << i + 1 << " " << (aurostd::isequal(aa.species_mass[i], 0.0) ? 1.0 : aa.species_mass[i]) << endl;
    }
    cout << endl;
    // write the atoms
    cout << "Atoms # atomic" << endl << endl;
    for (size_t iat = 0; iat < aa.atoms.size(); iat++) {
      cout << iat + 1 << " ";
      if (aa.atoms[iat].name_is_given == true) {
        cout << aa.atoms[iat].cleanname << " ";
      } else {
        cout << aa.atoms[iat].type + 1 << " ";
      }
      for (uint i = 1; i <= 3; i++) {
        if (std::abs(aa.atoms[iat].cpos(i)) < 10.0) {
          cout << " ";
        }
        if (!std::signbit(aa.atoms[iat].cpos(i))) {
          cout << " ";
        }
        cout << aa.atoms[iat].cpos(i) << " ";
      }
      cout << endl;
    }
    cout << endl;
    // write the velocities
    cout << "Velocities" << endl << endl;
    for (size_t iat = 0; iat < aa.atoms.size(); iat++) {
      cout << iat + 1 << " 0 0 0" << endl;
    }
    return cout;
  }
  // ----------------------------------------------------------------------

  cout << "NOT CODED YET" << endl;
  return cout;
}

// ***************************************************************************
// Function xstructure2qe
// ***************************************************************************
void xstructure::xstructure2qe() {
  ReScale(1.0);
  neg_scale = false;
  coord_flag = _COORDS_FRACTIONAL_;
  iomode = IOQE_GEOM;
  if (title.empty()) {
    buildGenericTitle();
  } // CO20171008 - pushed all of this to a function
}

// ***************************************************************************
// Function xstructure2vasp
// ***************************************************************************
void xstructure::xstructure2vasp() {
  ReScale(1.0); // DX20180425 - needs to be generic and rescale to 1.0; see other xstructure converters
  neg_scale = false;
  coord_flag = _COORDS_FRACTIONAL_;
  iomode = IOVASP_AUTO;
  if (title.empty()) {
    buildGenericTitle();
  } // CO20171008 - pushed all of this to a function
}

// ***************************************************************************
// Function xstructure2itc
// ***************************************************************************
void xstructure::xstructure2itc() {
  // CO20220613
  ReScale(1.0);
  neg_scale = false;
  coord_flag = _COORDS_FRACTIONAL_;
  const char iomode_orig = iomode; // save
  iomode = IOVASP_WYCKCAR;
  if (title.empty()) {
    buildGenericTitle();
  } // CO20171008 - pushed all of this to a function
  if (!partial_occupation_flag) {
    // CO20220715
    (*this).spacegroupnumber = (*this).SpaceGroup_ITC();
    (*this).lattice = (*this).standard_lattice_ITC; // need to update the lattice; may have rotated
  }
  stringstream ss;
  ss << (*this);
  (*this).clear();
  ss >> (*this);
  iomode = iomode_orig;
}

// ***************************************************************************
// Function xstructure2aims
// ***************************************************************************
void xstructure::xstructure2aims() {
  ReScale(1.0);
  neg_scale = false;
  coord_flag = _COORDS_FRACTIONAL_;
  iomode = IOAIMS_GEOM;
  if (title.empty()) {
    buildGenericTitle();
  } // CO20171008
}

// ***************************************************************************
// Function xstructure2abinit
// ***************************************************************************
void xstructure::xstructure2abinit() {
  ReScale(1.0);
  neg_scale = false;
  coord_flag = _COORDS_FRACTIONAL_;
  iomode = IOABINIT_GEOM;
  if (title.empty()) {
    buildGenericTitle();
  } // CO20171008 - pushed all of this to a function
}

// ***************************************************************************
// Function xstructure2cif
// ***************************************************************************
void xstructure::xstructure2cif() {
  // DX20190131
  ReScale(1.0);
  neg_scale = false;
  coord_flag = _COORDS_FRACTIONAL_;
  iomode = IOCIF;
  if (!partial_occupation_flag) {
    if ((*this).spacegroupnumber != 1) { // DX20250319 - only calculate if the CIF value was forced to be P1
      (*this).spacegroupnumber = (*this).SpaceGroup_ITC((*this).sym_eps, -1, (*this).setting_ITC, (*this).sym_eps_no_scan);
      (*this).lattice = (*this).standard_lattice_ITC; // need to update the lattice; may have rotated
    }
  }
  if (title.empty()) {
    buildGenericTitle();
  } // CO20171008 - pushed all of this to a function
  return;
}

// ***************************************************************************
// Function xstructure2abccar
// ***************************************************************************
void xstructure::xstructure2abccar() {
  // DX20190131
  ReScale(1.0);
  neg_scale = false;
  coord_flag = _COORDS_FRACTIONAL_;
  iomode = IOVASP_ABCCAR;
  if (title.empty()) {
    buildGenericTitle();
  } // CO20171008 - pushed all of this to a function
}

// ***************************************************************************
// Function xstructure2elk
// ***************************************************************************
void xstructure::xstructure2elk() {
  // DX20200313
  ReScale(1.0);
  neg_scale = false;
  coord_flag = _COORDS_FRACTIONAL_;
  iomode = IOELK_GEOM;
  if (title.empty()) {
    buildGenericTitle();
  } // CO20171008 - pushed all of this to a function
}

// ***************************************************************************
// Function xstructure2atat
// ***************************************************************************
void xstructure::xstructure2atat() {
  // SD20220123
  ReScale(1.0);
  neg_scale = false;
  coord_flag = _COORDS_FRACTIONAL_;
  iomode = IOATAT_STR;
  if (title.empty()) {
    buildGenericTitle();
  }
}

// ***************************************************************************
// Function xstructure2lmp
// ***************************************************************************
void xstructure::xstructure2lmp() {
  // SD20240111
  ReScale(1.0);
  neg_scale = false;
  coord_flag = _COORDS_FRACTIONAL_;
  iomode = IOLMP_DATA;
  if (title.empty()) {
    buildGenericTitle();
  }
}
