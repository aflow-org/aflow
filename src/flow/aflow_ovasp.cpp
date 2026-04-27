// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_automatic_template.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_hash.h"
#include "AUROSTD/aurostd_time.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xfile.h"
#include "AUROSTD/aurostd_xhttp.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xparser.h"
#include "AUROSTD/aurostd_xparser_json.h"
#include "AUROSTD/aurostd_xscalar.h"
#include "AUROSTD/aurostd_xtensor.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow.h"
#include "aflow_aflowrc.h"
#include "aflow_defs.h"
#include "aflow_xhost.h"
#include "aflowlib/aflowlib.h"
#include "flow/aflow_pflow.h"
#include "flow/aflow_xclasses.h"
#include "modules/SYM/aflow_symmetry.h"
#include "structure/aflow_xatom.h"

using std::cerr;
using std::cout;
using std::deque;
using std::endl;
using std::ifstream;
using std::istream;
using std::istringstream;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::setw;
using std::string;
using std::stringstream;
using std::vector;

using aurostd::xmatrix;
using aurostd::xvector;

//---------------------------------------------------------------------------------
// for bechmark
string time_delay(long double seconds) {
  return aurostd::utype2string<long double>(aurostd::get_seconds(seconds), 6);
}

//---------------------------------------------------------------------------------
// class xOUTCAR
//---------------------------------------------------------------------------------
ostream& operator<<(ostream& oss, const xOUTCAR& xOUT) {  // SC20200330
  const bool LDEBUG = (false || XHOST.DEBUG);
  const long double seconds = aurostd::get_seconds();
  if (LDEBUG) {
    cerr << XPID << "xOUTCAR::operator<<: ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << XPID << "xOUTCAR::operator<<: BEGIN (" << time_delay(seconds) << ")" << endl;
  }
  if (LDEBUG) {
    cerr << XPID << "xOUTCAR::operator<<: filename=[" << xOUT.filename << "]" << endl;
  }
  oss << " filename=" << xOUT.filename << endl;
  oss << " vcontent.size()=" << xOUT.vcontent.size() << endl;
  oss << " SYSTEM=" << xOUT.SYSTEM << endl;
  oss << " NELM=" << xOUT.NELM << endl; // CO20200624
  oss << " NIONS=" << xOUT.NIONS << endl;
  oss << " Efermi=" << xOUT.Efermi << endl;
  oss << " isLSCOUPLING=" << xOUT.isLSCOUPLING << endl;
  oss << " nelectrons=" << xOUT.nelectrons << endl; // AS20200528
  oss << " natoms=" << xOUT.natoms << endl;
  oss << " energy_cell=" << xOUT.energy_cell << endl;
  oss << " energy_atom=" << xOUT.energy_atom << endl;
  oss << " energy_atom=" << xOUT.energy_atom << endl;
  oss << " enthalpy_atom=" << xOUT.enthalpy_atom << endl;
  oss << " eentropy_cell=" << xOUT.eentropy_cell << endl;
  oss << " eentropy_atom=" << xOUT.eentropy_atom << endl;
  oss << " PV_cell=" << xOUT.PV_cell << endl;
  oss << " PV_atom=" << xOUT.PV_atom << endl;
  // xmatrix<double> stress;                                       // for aflowlib_libraries.cpp
  oss << " mag_cell=" << xOUT.mag_cell << endl;
  oss << " mag_atom=" << xOUT.mag_atom << endl;
  //   vector<double> vmag;                                          // for aflowlib_libraries.cpp
  //   vector<xvector<double> > vmag_noncoll;                        //DX20171205 - non-collinear
  oss << " volume_cell=" << xOUT.volume_cell << endl;
  oss << " volume_atom=" << xOUT.volume_atom << endl;
  oss << " pressure=" << xOUT.pressure << endl;
  oss << " pressure_residual=" << xOUT.pressure_residual << endl;
  oss << " Pulay_stress=" << xOUT.Pulay_stress << endl;
  //   vector<aurostd::xvector<double> > vforces;                    // for aflowlib_libraries.cpp
  //   vector<aurostd::xvector<double> > vpositions_cartesian;       // for aflowlib_libraries.cpp
  oss << " ENCUT=" << xOUT.ENCUT << endl;
  oss << " EDIFF=" << xOUT.EDIFF << endl;
  oss << " EDIFFG=" << xOUT.EDIFFG << endl;
  oss << " POTIM=" << xOUT.POTIM << endl;
  oss << " TEIN=" << xOUT.TEIN << endl;
  oss << " TEBEG=" << xOUT.TEBEG << endl;
  oss << " TEEND=" << xOUT.TEEND << endl;
  oss << " SMASS=" << xOUT.SMASS << endl;
  oss << " NPACO=" << xOUT.NPACO << endl;
  oss << " APACO=" << xOUT.APACO << endl;
  oss << " PSTRESS=" << xOUT.PSTRESS << endl;
  oss << " NBANDS=" << xOUT.NBANDS << endl;
  oss << " NKPTS=" << xOUT.NKPTS << endl;
  oss << " NSW=" << xOUT.NSW << endl;
  oss << " NBLOCK=" << xOUT.NBLOCK << endl;
  oss << " KBLOCK=" << xOUT.KBLOCK << endl;
  oss << " IBRION=" << xOUT.IBRION << endl;
  oss << " NFREE=" << xOUT.NFREE << endl;
  oss << " ISIF=" << xOUT.ISIF << endl;
  oss << " IWAVPR=" << xOUT.IWAVPR << endl;
  oss << " ISYM=" << xOUT.ISYM << endl;
  oss << " ISPIN=" << xOUT.ISPIN << endl;
  oss << " EMIN=" << xOUT.EMIN << endl;
  oss << " EMAX=" << xOUT.EMAX << endl;
  oss << " SIGMA=" << xOUT.SIGMA << endl;
  oss << " ISMEAR=" << xOUT.ISMEAR << endl;
  oss << " IALGO=" << xOUT.IALGO << endl;
  oss << " LDIAG=" << xOUT.LDIAG << endl;
  oss << " IMIX=" << xOUT.IMIX << endl;
  oss << " INIMIX=" << xOUT.INIMIX << endl;
  oss << " MIXPRE=" << xOUT.MIXPRE << endl;
  oss << " AMIX=" << xOUT.AMIX << endl;
  oss << " BMIX=" << xOUT.BMIX << endl;
  oss << " AMIX_MAG=" << xOUT.AMIX_MAG << endl;
  oss << " BMIX_MAG=" << xOUT.BMIX_MAG << endl;
  oss << " AMIN=" << xOUT.AMIN << endl;
  oss << " WC=" << xOUT.WC << endl;
  oss << " WEIMIN=" << xOUT.WEIMIN << endl;
  oss << " EBREAK=" << xOUT.EBREAK << endl;
  oss << " DEPER=" << xOUT.DEPER << endl;
  oss << " TIME=" << xOUT.TIME << endl;
  oss << " ENMAX=" << xOUT.ENMAX << "  vENMAX.size()=" << xOUT.vENMAX.size() << ": ";
  for (size_t i = 0; i < xOUT.vENMAX.size(); i++) {
    oss << xOUT.vENMAX[i] << " ";
  }
  oss << endl;
  oss << " ENMIN=" << xOUT.ENMIN << "  vENMIN.size()=" << xOUT.vENMIN.size() << ": ";
  for (size_t i = 0; i < xOUT.vENMIN.size(); i++) {
    oss << xOUT.vENMIN[i] << " ";
  }
  oss << endl;
  oss << " POMASS_sum=" << xOUT.POMASS_sum << " POMASS_min=" << xOUT.POMASS_min << " POMASS_max=" << xOUT.POMASS_max << " vPOMASS.size()=" << xOUT.vPOMASS.size() << ": ";
  for (size_t i = 0; i < xOUT.vPOMASS.size(); i++) {
    oss << xOUT.vPOMASS[i] << " ";
  }
  oss << " " << endl;
  oss << " ZVAL_sum=" << xOUT.ZVAL_sum << " ZVAL_min=" << xOUT.ZVAL_min << " ZVAL_max=" << xOUT.ZVAL_max << " vZVAL.size()=" << xOUT.vZVAL.size() << ": ";
  for (size_t i = 0; i < xOUT.vZVAL.size(); i++) {
    oss << xOUT.vZVAL[i] << " ";
  }
  oss << " " << endl;
  oss << " EATOM_min=" << xOUT.EATOM_min << " EATOM_max=" << xOUT.EATOM_max << " vEATOM.size()=" << xOUT.vEATOM.size() << ": ";
  for (size_t i = 0; i < xOUT.vEATOM.size(); i++) {
    oss << xOUT.vEATOM[i] << " ";
  }
  oss << " " << endl;
  oss << " RCORE_min=" << xOUT.RCORE_min << " RCORE_max=" << xOUT.RCORE_max << " vRCORE.size()=" << xOUT.vRCORE.size() << ": ";
  for (size_t i = 0; i < xOUT.vRCORE.size(); i++) {
    oss << xOUT.vRCORE[i] << " ";
  }
  oss << " " << endl;
  oss << " RWIGS_min=" << xOUT.RWIGS_min << " RWIGS_max=" << xOUT.RWIGS_max << " vRWIGS.size()=" << xOUT.vRWIGS.size() << ": ";
  for (size_t i = 0; i < xOUT.vRWIGS.size(); i++) {
    oss << xOUT.vRWIGS[i] << " ";
  }
  oss << " " << endl;
  oss << " EAUG_min=" << xOUT.EAUG_min << " EAUG_max=" << xOUT.EAUG_max << " vEAUG.size()=" << xOUT.vEAUG.size() << ": ";
  for (size_t i = 0; i < xOUT.vEAUG.size(); i++) {
    oss << xOUT.vEAUG[i] << " ";
  }
  oss << " " << endl;
  oss << " RAUG_min=" << xOUT.RAUG_min << " RAUG_max=" << xOUT.RAUG_max << " vRAUG.size()=" << xOUT.vRAUG.size() << ": ";
  for (size_t i = 0; i < xOUT.vRAUG.size(); i++) {
    oss << xOUT.vRAUG[i] << " ";
  }
  oss << " " << endl;
  oss << " RMAX_min=" << xOUT.RMAX_min << " RMAX_max=" << xOUT.RMAX_max << " vRMAX.size()=" << xOUT.vRMAX.size() << ": ";
  for (size_t i = 0; i < xOUT.vRMAX.size(); i++) {
    oss << xOUT.vRMAX[i] << " ";
  }
  oss << " " << endl;
  oss << " vTITEL.size()=" << xOUT.vTITEL.size() << ": ";
  for (size_t i = 0; i < xOUT.vTITEL.size(); i++) {
    oss << xOUT.vTITEL[i] << " ";
  }
  oss << endl;
  oss << " vLEXCH.size()=" << xOUT.vLEXCH.size() << ": ";
  for (size_t i = 0; i < xOUT.vLEXCH.size(); i++) {
    oss << xOUT.vLEXCH[i] << " ";
  }
  oss << endl;
  oss << " pp_type=" << xOUT.pp_type << endl;
  oss << " species.size()=" << xOUT.species.size() << ": ";
  for (size_t i = 0; i < xOUT.species.size(); i++) {
    oss << xOUT.species[i] << " ";
  }
  oss << endl;
  oss << " species_Z.size()=" << xOUT.species_Z.size() << ": ";
  for (size_t i = 0; i < xOUT.species_Z.size(); i++) {
    oss << xOUT.species_Z[i] << " ";
  }
  oss << endl;
  oss << " species_pp.size()=" << xOUT.species_pp.size() << ": ";
  for (size_t i = 0; i < xOUT.species_pp.size(); i++) {
    oss << xOUT.species_pp[i] << " ";
  }
  oss << endl;
  oss << " species_pp_type.size()=" << xOUT.species_pp_type.size() << ": ";
  for (size_t i = 0; i < xOUT.species_pp_type.size(); i++) {
    oss << xOUT.species_pp_type[i] << " ";
  }
  oss << endl;
  oss << " species_pp_version.size()=" << xOUT.species_pp_version.size() << ": ";
  for (size_t i = 0; i < xOUT.species_pp_version.size(); i++) {
    oss << xOUT.species_pp_version[i] << " ";
  }
  oss << endl;
  oss << " species_pp_AUID.size()=" << xOUT.species_pp_AUID.size() << ": ";
  for (size_t i = 0; i < xOUT.species_pp_AUID.size(); i++) {
    oss << xOUT.species_pp_AUID[i] << " ";
  }
  oss << endl;
  oss << " species_pp_AUID_collisions.size()=" << xOUT.species_pp_AUID_collisions.size() << ": ";
  for (size_t i = 0; i < xOUT.species_pp_AUID_collisions.size(); i++) {
    oss << xOUT.species_pp_AUID_collisions[i] << " ";
  }
  oss << endl;
  oss << " species_pp_groundstate_energy.size()=" << xOUT.species_pp_groundstate_energy.size() << ": ";
  for (size_t i = 0; i < xOUT.species_pp_groundstate_energy.size(); i++) {
    oss << xOUT.species_pp_groundstate_energy[i] << " ";
  }
  oss << endl;
  oss << " species_pp_groundstate_structure.size()=" << xOUT.species_pp_groundstate_structure.size() << ": ";
  for (size_t i = 0; i < xOUT.species_pp_groundstate_structure.size(); i++) {
    oss << xOUT.species_pp_groundstate_structure[i] << " ";
  }
  oss << endl;
  // deque<deque<double> > species_pp_vLDAU;  // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  oss << " isKIN=" << xOUT.isKIN << endl;
  oss << " isMETAGGA=" << xOUT.isMETAGGA << endl;
  oss << " METAGGA=" << xOUT.METAGGA << endl;
  oss << " nweights=" << xOUT.nweights << endl;
  oss << " nkpoints_irreducible=" << xOUT.nkpoints_irreducible << endl;
  // vector<aurostd::xvector<double> > vkpoint_reciprocal;         // kpoints reading
  // vector<aurostd::xvector<double> > vkpoint_cartesian;          // kpoints reading
  oss << " vweights.size()=" << xOUT.vweights.size() << ": ";
  for (size_t i = 0; i < xOUT.vweights.size(); i++) {
    oss << xOUT.vweights[i] << " ";
  }
  oss << endl;
  oss << " calculation_time=" << xOUT.calculation_time << endl;
  oss << " calculation_memory=" << xOUT.calculation_memory << endl;
  oss << " calculation_cores=" << xOUT.calculation_cores << endl;
  // xstructure xstr;                                              // for GetBandGap()
  // EFFECTIVE MASSES
  oss << " band_index.size()=" << xOUT.band_index.size() << ": ";
  for (size_t i = 0; i < xOUT.band_index.size(); i++) {
    oss << xOUT.band_index[i] << " ";
  }
  oss << endl;
  oss << " carrier_spin.size()=" << xOUT.carrier_spin.size() << ": ";
  for (size_t i = 0; i < xOUT.carrier_spin.size(); i++) {
    oss << xOUT.carrier_spin[i] << " ";
  }
  oss << endl;
  oss << " carrier_type.size()=" << xOUT.carrier_type.size() << ": ";
  for (size_t i = 0; i < xOUT.carrier_type.size(); i++) {
    oss << xOUT.carrier_type[i] << " ";
  }
  oss << endl;
  // vector<vector<double> > extrema_cart_coord;
  // vector<vector<double> > effective_mass_axes;
  oss << " equivalent_valley.size()=" << xOUT.equivalent_valley.size() << ": ";
  for (size_t i = 0; i < xOUT.equivalent_valley.size(); i++) {
    oss << xOUT.equivalent_valley[i] << " ";
  }
  oss << endl;
  oss << " effective_mass_DOS.size()=" << xOUT.effective_mass_DOS.size() << ": ";
  for (size_t i = 0; i < xOUT.effective_mass_DOS.size(); i++) {
    oss << xOUT.effective_mass_DOS[i] << " ";
  }
  oss << endl;
  oss << " effective_mass_COND.size()=" << xOUT.effective_mass_COND.size() << ": ";
  for (size_t i = 0; i < xOUT.effective_mass_COND.size(); i++) {
    oss << xOUT.effective_mass_COND[i] << " ";
  }
  oss << endl;
  oss << " mass_elec_dos.size()=" << xOUT.mass_elec_dos.size() << ": ";
  for (size_t i = 0; i < xOUT.mass_elec_dos.size(); i++) {
    oss << xOUT.mass_elec_dos[i] << " ";
  }
  oss << endl;
  oss << " mass_hole_dos.size()=" << xOUT.mass_hole_dos.size() << ": ";
  for (size_t i = 0; i < xOUT.mass_hole_dos.size(); i++) {
    oss << xOUT.mass_hole_dos[i] << " ";
  }
  oss << endl;
  oss << " mass_elec_conduction.size()=" << xOUT.mass_elec_conduction.size() << ": ";
  for (size_t i = 0; i < xOUT.mass_elec_conduction.size(); i++) {
    oss << xOUT.mass_elec_conduction[i] << " ";
  }
  oss << endl;
  oss << " mass_hole_conduction.size()=" << xOUT.mass_hole_conduction.size() << ": ";
  for (size_t i = 0; i < xOUT.mass_hole_conduction.size(); i++) {
    oss << xOUT.mass_hole_conduction[i] << " ";
  }
  oss << endl;
  oss << " conduction_band_min.size()=" << xOUT.conduction_band_min.size() << ": ";
  for (size_t i = 0; i < xOUT.conduction_band_min.size(); i++) {
    oss << xOUT.conduction_band_min[i] << " ";
  }
  oss << endl;
  oss << " conduction_band_min_net=" << xOUT.conduction_band_min_net << endl;
  oss << " valence_band_max.size()=" << xOUT.valence_band_max.size() << ": ";
  for (size_t i = 0; i < xOUT.valence_band_max.size(); i++) {
    oss << xOUT.valence_band_max[i] << " ";
  }
  oss << endl;
  oss << " valence_band_max_net=" << xOUT.valence_band_max_net << endl;
  oss << " Egap.size()=" << xOUT.Egap.size() << ": ";
  for (size_t i = 0; i < xOUT.Egap.size(); i++) {
    oss << xOUT.Egap[i] << " ";
  }
  oss << endl;
  oss << " Egap_net=" << xOUT.Egap_net << endl;
  oss << " Egap_fit.size()=" << xOUT.Egap_fit.size() << ": ";
  for (size_t i = 0; i < xOUT.Egap_fit.size(); i++) {
    oss << xOUT.Egap_fit[i] << " ";
  }
  oss << endl;
  oss << " Egap_fit_net=" << xOUT.Egap_fit_net << endl;
  oss << " Egap_type.size()=" << xOUT.Egap_type.size() << ": ";
  for (size_t i = 0; i < xOUT.Egap_type.size(); i++) {
    oss << xOUT.Egap_type[i] << " ";
  }
  oss << endl;
  oss << " Egap_type_net=" << xOUT.Egap_type_net << endl;
  // ----------------------------------------------------------------------
  if (LDEBUG) {
    cerr << XPID << "xOUTCAR::operator<<: END (" << time_delay(seconds) << ")" << endl;
  }
  return oss;
} // SC20200330

bool xOUTCAR::GetProperties(const string& stringIN, bool QUIET) {
  stringstream sss;
  sss.str(stringIN);
  if (filename.empty()) {
    filename = "string";
  }
  return xOUTCAR::GetProperties(sss, QUIET);
}

bool xOUTCAR::GetPropertiesFile(const string& fileIN, bool QUIET) {
  stringstream sss;
  aurostd::compressfile2stringstream(fileIN, sss);
  filename = fileIN;
  return xOUTCAR::GetProperties(sss, QUIET);
}

bool xOUTCAR::GetPropertiesFile(const string& fileIN, uint natoms_check, bool QUIET) {
  const bool flag = GetPropertiesFile(fileIN, QUIET);
  if (aurostd::abs(natoms_check - natoms) > 0.1) {
    stringstream message;
    message << "natoms_check(" << natoms_check << ")!= (int) natoms(" << natoms << ")";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
  }
  return flag;
}

bool xOUTCAR::GetPropertiesUrlFile(const string& url, const string& file, bool QUIET) {
  const string tmpfile = aurostd::TmpFileCreate("xOUTCAR_GetProperties"); // CO20200502 - threadID
  aurostd::httpGetFileStatus(url + "/" + file, tmpfile);
  const bool out = GetPropertiesFile(tmpfile, QUIET); // CO20200404 - added QUIET
  filename = "url=" + url;  // CO20210315
  aurostd::RemoveFile(tmpfile);
  return out;
}

/// @brief get the correct entries from a line
/// @param line from the OUTCAR
/// @param expected_count expected number of space seperated entries
/// @authors
/// @mod{CO,20170725,created}
/// @mod{SD,20240204,renamed function; added doxy}
/// @return vector of correct entries
/// @note FIRST FIX: NEGATIVE SIGN
/// @note this function should fix the last line
/// @note    volume of cell :      369.80
/// @note      direct lattice vectors                 reciprocal lattice vectors
/// @note    -2.730747137  2.730747137 12.397646334     0.000000000  0.183100073  0.040330236
/// @note     2.730747137 -2.730747137 12.397646334     0.183100073  0.000000000  0.040330236
/// @note     2.730747137  2.730747137-12.397646334     0.183100073  0.183100073  0.000000000
/// @note
/// @note SECOND FIX: LARGE NUMBERS
/// @note this function should fix the last line
/// @note      direct lattice vectors                 reciprocal lattice vectors
/// @note     2.156793936 -3.735676679  0.000000000     0.231825578 -0.133844560  0.000000000
/// @note     2.156793936  3.735676679  0.000000000     0.231825578  0.133844560  0.000000000
/// @note     0.000000000  0.000000000109.286277550     0.000000000  0.000000000  0.009150280
vector<string> xOUTCAR::GetCorrectEntriesFromLine(const string& line, const uint expected_count) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  vector<string> tokens;
  aurostd::string2tokens(line, tokens);
  if (tokens.size() == expected_count) {
    return tokens;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " issuing fix for bad lattice vectors (negative sign) on this line: " << line << endl;
  }
  // try fixing negative sign first
  vector<string> neg_tokens;
  vector<string> _tokens;
  size_t pos;
  bool first_number_negative;
  for (size_t i = 0; i < tokens.size(); i++) {
    first_number_negative = false;
    pos = tokens[i].find('-', 0);
    if (pos != std::string::npos) {
      if (pos == 0) {
        first_number_negative = true;
      }
      aurostd::string2tokens(tokens[i], _tokens, "-");
      neg_tokens.push_back((first_number_negative ? "-" : "") + _tokens[0]);
      for (size_t ii = 1; ii < _tokens.size(); ii++) {
        neg_tokens.push_back("-" + _tokens[ii]);
      }
    } else {
      neg_tokens.push_back(tokens[i]);
    }
  }
  if (neg_tokens.size() == expected_count) {
    return neg_tokens;
  }
  // negative sign fix not enough, now  we need to look at large number problems
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " issuing fix for bad lattice vectors (large numbers) on this line: " << line << endl;
  }
  vector<string> dec_tokens;
  string good_num;
  _tokens.clear();
  // first, let's find the precision, if we can
  for (size_t i = 0; i < neg_tokens.size() && good_num.empty(); i++) {
    aurostd::string2tokens(neg_tokens[i], _tokens, ".");
    if (_tokens.size() == 2) {
      good_num = neg_tokens[i];
    }  // this number is ok
  }
  if (good_num.empty()) {
    dec_tokens.clear();
    return dec_tokens;
  } // we have no hope
  //_tokens contains precision
  const size_t precision = _tokens[1].size();
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " precision of numbers on this line: " << precision << endl;
  }
  // let's build numbers with the right count of digits
  string num;
  bool fidelity = true;
  for (size_t i = 0; i < neg_tokens.size() && fidelity; i++) {
    aurostd::string2tokens(neg_tokens[i], _tokens, ".");
    if (_tokens.size() == 2) {
      dec_tokens.push_back(neg_tokens[i]);
      continue;
    }  // this number is ok
    // we need to parse this number per precision
    num = _tokens[0];
    for (size_t j = 1; j < _tokens.size() - 1 && fidelity; j++) {
      if (_tokens[j].size() < precision + 1) {
        fidelity = false;
        break;
      }  // should be at least as big as precision + 1 (full number after "." of first number + 0 of next number)
      num += ".";
      for (size_t k = 0; k < precision; k++) {
        if (!isdigit(_tokens[j][k])) {
          fidelity = false;
          break;
        }  // check every character is digit
        num += _tokens[j][k];
      }
      dec_tokens.push_back(num);
      num = "";
      for (size_t k = precision; k < _tokens[j].size(); k++) {
        if (!isdigit(_tokens[j][k])) {
          fidelity = false;
          break;
        }  // check every character is digit, we took care of negative signs already
        num += _tokens[j][k];
      }
    }
    num += ".";
    const size_t j = _tokens.size() - 1;
    if (_tokens[j].size() != precision) {
      fidelity = false;
      break;
    }
    for (size_t k = 0; k < _tokens[j].size() && fidelity; k++) {
      if (!isdigit(_tokens[j][k])) {
        fidelity = false;
        break;
      }  // check every character is digit
      num += _tokens[j][k];
    }
    dec_tokens.push_back(num);
  }
  // run through all numbers just to make sure they all have the right precision
  for (size_t i = 0; i < dec_tokens.size() && fidelity; i++) {
    aurostd::string2tokens(dec_tokens[i], _tokens, ".");
    if (_tokens.size() != 2 || _tokens[1].size() != precision) {
      fidelity = false;
      break;
    }
  }
  if (!fidelity || dec_tokens.size() != expected_count) {
    dec_tokens.clear();
    return dec_tokens;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " repaired vector: ";
    for (size_t i = 0; i < dec_tokens.size(); i++) {
      cerr << dec_tokens[i] << " ";
    }
    cerr << endl;
  }
  return dec_tokens;
}

bool xOUTCAR::GetProperties(const stringstream& stringstreamIN, bool QUIET) {
  const bool LDEBUG = (false || XHOST.DEBUG || !QUIET);
  stringstream message;
  const bool force_exit = XHOST.POSTPROCESS; // SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only

  bool ERROR_flag = false;
  clear(); // so it does not mess up vector/deque

  const long double seconds = aurostd::get_seconds();
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BEGIN (" << time_delay(seconds) << ")" << endl;
  }
  stringstream sss;
  sss.str(stringstreamIN.str());
  content = stringstreamIN.str();
  vcontent.clear();
  vector<string> vline;
  vector<string> tokens;
  vector<string> vcontentRED;
  aurostd::string2vectorstring(content, vcontent);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vcontent.size()=" << vcontent.size() << " (" << time_delay(seconds) << ")" << endl;
  }
  for (size_t iline = 0; iline < vcontent.size(); iline++) { // NEW - FROM THE BACK
    const string saus = vcontent[iline];
    if (aurostd::substring2bool(saus, "=")) {
      vcontentRED.push_back(saus);
    }
    if (aurostd::substring2bool(saus, "total")) {
      vcontentRED.push_back(saus); // total
    }
    if (aurostd::substring2bool(saus, "number")) {
      vcontentRED.push_back(saus); // number
    }
    if (aurostd::substring2bool(saus, "bands")) {
      vcontentRED.push_back(saus); // bands
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vcontentRED.size()=" << vcontentRED.size() << " (" << time_delay(seconds) << ")" << endl;
  }

  string line;
  if (filename.empty()) {
    filename = "stringstream";
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " filename=[" << filename << "]" << endl;
  }
  if (LDEBUG) {
    cerr.precision(12);
  }

  // ----------------------------------------------------------------------
  // SC
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD SYSTEM (" << time_delay(seconds) << ")" << endl;
  }
  SYSTEM = "";
  line = "";
  for (size_t iline = 0; iline < vcontent.size(); iline++) { // NEW - FROM THE BACK
    if (aurostd::substring2bool(vcontent[iline], "SYSTEM")) { // VASP
      if (aurostd::substring2bool(vcontent[iline], "=")) { // VASP
        line = vcontent[iline];
        break;
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " line=" << line << endl;
  }
  if (!line.empty()) {
    aurostd::string2tokens(line, tokens, "="); // cerr << tokens.size() << endl;
    SYSTEM = aurostd::RemoveWhiteSpaces(tokens.at(1));
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " SYSTEM=" << SYSTEM << endl;
  }
  // ----------------------------------------------------------------------
  // KY STUFF DONT TOUCH - GET Number of IONS and Fermi
  // ----------------------------------------------------------------------
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD KESONG STUFF (" << time_delay(seconds) << ")" << endl;
  }
  const string anchor_word_NELM = "NELM"; // CO20200624
  const string anchor_word_NIONS = "NIONS";
  const string anchor_word_LSORBIT = "LSORBIT =";
  const string anchor_word_Efermi = "E-fermi";
  const string anchor_word_EFIELD_PEAD = "EFIELD_PEAD"; // CO20210315
  string tmp;
  efield_pead.resize(3);  // CO20210315
  while (getline(sss, line)) {
    if (line.find(anchor_word_NELM) != string::npos) {  // CO20200624
      aurostd::string2tokens(line, tokens, ";");
      if (tokens.empty()) {
        continue;
      }
      tmp = tokens[0];
      aurostd::string2tokens(tmp, tokens, "=");
      if (tokens.size() < 2) {
        continue;
      }
      tmp = tokens[1];
      NELM = aurostd::string2utype<int>(aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tmp));
    }
    if (line.find(anchor_word_NIONS) != string::npos) {
      aurostd::string2tokens(line, tokens, " ");
      NIONS = aurostd::string2utype<int>(tokens.at(tokens.size() - 1)); // last one
    }
    if (line.find(anchor_word_LSORBIT) != string::npos) {
      aurostd::string2tokens(line, tokens, " ");
      const string LSlabel = tokens.at(2);
      if (LSlabel == "T") {
        isLSCOUPLING = true;
      }
      if (LSlabel == "F") {
        isLSCOUPLING = false;
      }
    }
    if (line.find(anchor_word_Efermi) != string::npos) {
      aurostd::string2tokens(line, tokens, " ");
      Efermi = aurostd::string2utype<double>(tokens.at(2));
    }
    if (line.find(anchor_word_EFIELD_PEAD) != string::npos) { // CO20210315 - moved from aflow_ivasp XVASP_INCAR_EFIELD_PEAD()
      aurostd::string2tokens(line, tokens, "=");
      if (tokens.size() < 2) {
        continue;
      }
      tmp = tokens[1];
      aurostd::string2tokens(tmp, tokens);
      if (tokens.size() < 3) {
        continue;
      }
      for (int i = 0; i < 3; i++) {
        efield_pead[i + 1] = aurostd::string2utype<double>(tokens[i]);
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " NELM=" << NELM << endl; // CO20200624
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " NIONS=" << NIONS << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " isLSCOUPLING=" << isLSCOUPLING << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " Efermi=" << Efermi << endl;
  }
  natoms = (double) NIONS;
  // ----------------------------------------------------------------------
  // LOAD EENTROPY DATA
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD EENTROPY DATA (" << time_delay(seconds) << ")" << endl;
  }
  for (int iline = (int) vcontent.size() - 1; iline >= 0; iline--) {  // NEW - FROM THE BACK
    if (aurostd::substring2bool(vcontent[iline], "entropy")) {
      if (aurostd::substring2bool(vcontent[iline], "EENTRO")) {
        if (aurostd::substring2bool(vcontent[iline], "T*S")) {
          line = vcontent[iline];
          break;
        }
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " line=" << line << endl;
  }
  aurostd::string2tokens(line, tokens, "=");
  if (tokens.empty()) {
    message << "Wrong number of entries (entropy) in OUTCAR; line=[ " << line << "]" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  }
  eentropy_cell = 0;
  if (tokens.size() > 1) {
    eentropy_cell = aurostd::string2utype<double>(tokens.at(1));
  }
  eentropy_atom = eentropy_cell / natoms;
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " eentropy_cell=" << eentropy_cell << endl;
  }
  // LOAD ENERGY DATA  // without energy of electrons
  // CO20211109 - notes about which energy to pick
  // we report the energy without entropy (0K energy).
  // the sigma->0 is an extrapolation calculation.
  // the only way to know the energy at sigma=0 is to do the calculation at
  // 0, but then the system would not converge.
  // it does not matter much which energy to pick: raw energy values of VASP are arbitrary.
  // focus on delta energy (just be consistent).
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD ENERGY DATA" << endl;
  }
  line = "";
  for (int iline = (int) vcontent.size() - 1; iline >= 0; iline--) {  // NEW - FROM THE BACK
    if (aurostd::substring2bool(vcontent[iline], "energy")) { // VASP
      if (aurostd::substring2bool(vcontent[iline], "without")) { // VASP
        if (aurostd::substring2bool(vcontent[iline], "entropy")) {
          line = vcontent[iline];
          break;
        }  // VASP
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " line=" << line << endl;
  }
  aurostd::string2tokens(line, tokens, "=");
  if (tokens.empty()) {
    message << "Wrong number of entries (energy_1) in OUTCAR; line=[ " << line << "]" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  }
  if (tokens.size() > 1) {
    line = tokens.at(1);
  }
  aurostd::string2tokens(line, tokens, "e");
  if (tokens.empty()) {
    message << "Wrong number of entries (energy_2) in OUTCAR; line=[ " << line << "]" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  }
  if (!tokens.empty()) {
    energy_cell = aurostd::string2utype<double>(tokens.at(0));
  }
  energy_atom = energy_cell / natoms;
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " energy_cell=" << energy_cell << " energy_atom=" << energy_atom << endl;
  }
  // ----------------------------------------------------------------------
  // LOAD PV DATA (IF PRESENT)
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD PV DATA (" << time_delay(seconds) << ")" << endl;
  }
  line = "";
  for (int iline = (int) vcontent.size() - 1; iline >= 0; iline--) {  // NEW - FROM THE BACK
    if (aurostd::substring2bool(vcontent[iline], "TOTEN")) {  // VASP
      if (aurostd::substring2bool(vcontent[iline], "P")) {  // VASP
        if (aurostd::substring2bool(vcontent[iline], "V")) {  // VASP
          if (!aurostd::substring2bool(vcontent[iline], "VPU")) {    // SOME VPU/CPU stuff ORTHCH  // VASP
            if (!aurostd::substring2bool(vcontent[iline], "CPU")) { // SOME VPU/CPU stuff ORTHCH  // VASP
              line = vcontent[iline];
              break;
            }
          }
        }
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " line=" << line << endl;
  }
  if (line.length() < 10) {
    PV_cell = 0.0;
    PV_atom = PV_cell / natoms;
  } else {
    aurostd::string2tokens(line, tokens, "=");
    if (tokens.size() != 3) {
      message << "Wrong number of entries (PV) in OUTCAR; line=[ " << line << "]" << "   filename=[" << filename << "]";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
      ERROR_flag = true;
    }
    if (tokens.size() > 2) {
      PV_cell = aurostd::string2utype<double>(tokens.at(2));
    }
    PV_atom = PV_cell / natoms;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PV_cell=" << PV_cell << " PV_atom=" << PV_atom << endl;
    }
  }

  // correct if PV

  if (aurostd::abs(PV_atom) < PRESSURE_ZERO_ENTHALPY_ENERGY) {
    enthalpy_cell = energy_cell;  // default without PV
    enthalpy_atom = energy_atom;   // default without PV
  } else {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " pressure=" << pressure << endl;
    }
    // LOAD ENTHALPY DATA IF P!=0
    line = "";
    for (int iline = (int) vcontent.size() - 1; iline >= 0; iline--) {  // NEW
      if (aurostd::substring2bool(vcontent[iline], "TOTEN")) { // VASP
        if (aurostd::substring2bool(vcontent[iline], "enthalpy")) { // VASP
          line = vcontent[iline];
          break;
        }
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " line=" << line << endl;
    }
    aurostd::string2tokens(line, tokens, " ");

    if (!tokens.empty()) {
      for (int i = tokens.size() - 2; i >= 0; i--)
        //    for(size_t i=0;i<tokens.size();i++)
      { // CO20200106 - patching for auto-indenting
        if (aurostd::substring2bool(tokens.at(i), "=")) {
          if (aurostd::substring2bool(tokens.at(i - 1), "TOTEN")) {
            enthalpy_cell = aurostd::string2utype<double>(tokens.at(i + 1));
            enthalpy_atom = enthalpy_cell / natoms;
            if (LDEBUG) {
              cerr << __AFLOW_FUNC__ << " enthalpy_cell=" << enthalpy_cell << " enthalpy_atom=" << enthalpy_atom << endl;
            }
          }
        }
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " enthalpy_cell=" << enthalpy_cell << " enthalpy_atom=" << enthalpy_atom << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " energy_cell=" << energy_cell << " energy_atom=" << energy_atom << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " PV_cell=" << PV_cell << " PV_atom=" << PV_atom << endl;
  }

  // ----------------------------------------------------------------------
  // LOAD PVstress DATA (IF PRESENT)
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD STRESS DATA (" << time_delay(seconds) << ")" << endl;
  }
  line = "";
  for (int iline = (int) vcontent.size() - 1; iline >= 0; iline--) {  // NEW - FROM THE BACK
    if (aurostd::substring2bool(vcontent[iline], "external")) {  // VASP
      if (aurostd::substring2bool(vcontent[iline], "pressure")) {  // VASP
        if (iline > 1) {
          if (aurostd::substring2bool(vcontent[iline - 1], "kB")) { // VASP
            line = vcontent[iline - 1];
            break;
          }
        }
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " line=" << line << endl;
  }
  stress.clear();
  aurostd::string2tokens(line, tokens, " ");
  if (tokens.size() >= 8) {
    stress(1, 1) = aurostd::string2utype<double>(tokens.at(2));
    stress(2, 2) = aurostd::string2utype<double>(tokens.at(3));
    stress(3, 3) = aurostd::string2utype<double>(tokens.at(4));
    stress(1, 2) = aurostd::string2utype<double>(tokens.at(5));
    stress(2, 1) = aurostd::string2utype<double>(tokens.at(5));
    stress(2, 3) = aurostd::string2utype<double>(tokens.at(6));
    stress(3, 2) = aurostd::string2utype<double>(tokens.at(6));
    stress(1, 3) = aurostd::string2utype<double>(tokens.at(7));
    stress(3, 1) = aurostd::string2utype<double>(tokens.at(7));
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " stress=" << endl << stress << endl;
  }
  //   Direction    XX        YY        ZX        XY       YZ       ZX
  //   in kB       -5.93     -5.93    -19.31      0.00      0.00     -0.00
  // ----------------------------------------------------------------------
  // LOAD VOLUME DATA (IF PRESENT)
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD volume DATA (" << time_delay(seconds) << ")" << endl;
  }
  line = "";
  for (int iline = (int) vcontent.size() - 1; iline >= 0; iline--) {  // NEW FROM THE BACK
    if (aurostd::substring2bool(vcontent[iline], "volume")) {
      if (aurostd::substring2bool(vcontent[iline], "of")) {
        if (aurostd::substring2bool(vcontent[iline], "cell")) {
          line = vcontent[iline];
          break;
        }
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " line=" << line << endl;
  }
  volume_cell = 0.0;
  aurostd::string2tokens(line, tokens, " ");
  if (tokens.size() >= 5) {
    if (tokens.at(0) == "volume") {
      volume_cell = aurostd::string2utype<double>(tokens.at(tokens.size() - 1));
    }
  }
  volume_atom = volume_cell / natoms;
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " volume_cell=" << volume_cell << " volume_atom=" << volume_atom << endl;
  }
  // ----------------------------------------------------------------------
  // LOAD PRESSURE_RESIDUAL/PULAY DATA
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD PRESSURE_RESIDUAL/PULAY DATA (" << time_delay(seconds) << ")" << endl;
  }
  line = "";
  for (int iline = (int) vcontent.size() - 1; iline >= 0; iline--) {  // NEW - FROM THE BACK
    if (aurostd::substring2bool(vcontent[iline], "Pullay") || aurostd::substring2bool(vcontent[iline], "pullay") || aurostd::substring2bool(vcontent[iline], "Pulay") ||
        aurostd::substring2bool(vcontent[iline], "pulay")) {  // VASP
      if (aurostd::substring2bool(vcontent[iline], "stress")) { // VASP
        line = vcontent[iline];
        break;
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " line=" << line << endl;
  }
  aurostd::string2tokens(line, tokens, "=");
  if (tokens.empty()) {
    message << "Wrong number of entries (Pulay stress) in OUTCAR; line=[ " << line << "]" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  }

  pressure_residual = 0;
  Pulay_stress = 0;
  if (!tokens.empty()) {
    const string line_Pulay_stress = tokens.at(tokens.size() - 1);
    const string line_pressure_residual = tokens.at(tokens.size() - 2);
    // PRESSURE
    pressure_residual = 0.0;
    aurostd::string2tokens(line_pressure_residual, tokens, " ");
    if (!tokens.empty()) {
      pressure_residual = aurostd::string2utype<double>(tokens.at(0));
    }
    // PULAY
    Pulay_stress = 0.0;
    aurostd::string2tokens(line_Pulay_stress, tokens, " ");
    if (!tokens.empty()) {
      Pulay_stress = aurostd::string2utype<double>(tokens.at(0));
    }
  }

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " pressure_residual=" << pressure_residual << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " Pulay_stress=" << Pulay_stress << endl;
  }

  // ----------------------------------------------------------------------
  // LOAD SPIN
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD SPIN DATA (" << time_delay(seconds) << ")" << endl;
  }
  line = "";
  for (int iline = (int) vcontent.size() - 1; iline >= 0; iline--) {  // NEW FROM THE BACK
    if (aurostd::substring2bool(vcontent[iline], "magnetization")) {
      if (aurostd::substring2bool(vcontent[iline], "number")) {
        line = vcontent[iline];
        break;
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " line=" << line << endl;
  }
  mag_cell = 0.0;
  aurostd::string2tokens(line, tokens, " ");
  if (tokens.size() >= 6) {
    if (tokens.at(0) == "number" && tokens.at(4) == "magnetization" && tokens.size() == 6) {
      mag_cell = aurostd::string2utype<double>(tokens.at(5)); // JX seems to work
    }
  }
  mag_atom = mag_cell / natoms;
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " mag_cell=" << mag_cell << " mag_atom=" << mag_atom << endl;
  }
  // ----------------------------------------------------------------------
  // LOAD SPIN DECOMPOSITION
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD SPIN DECOMPOSITION DATA" << endl;
  }
  line = "";
  uint mline = 0;
  // DX20171205 - Check for magnetization z (non-collinear) - START
  vector<double> vmag_z;
  bool found_magz_line = false;
  for (int iline = (int) vcontent.size() - 1; iline >= 0; iline--) {  // NEW FROM THE BACK
    if (aurostd::substring2bool(vcontent[iline], "magnetization")) {
      if (aurostd::substring2bool(vcontent[iline], "(z)")) {
        mline = iline + 4;
        found_magz_line = true;
        break;
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " magnetization (z) line=" << line << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " magnetization (z) mline=" << mline << endl;
  }
  if (found_magz_line) {
    for (uint iline = mline; iline < mline + natoms; iline++) {
      aurostd::string2tokens(vcontent[iline], tokens, " ");
      if (tokens.size() >= 5) {
        vmag_z.push_back(aurostd::string2utype<double>(tokens.at(tokens.size() - 1)));
      }
    }
  }
  // DX20171205 - Check for magnetization z (non-collinear) - END
  // DX20171205 - Check for magnetization y (non-collinear) - START
  vector<double> vmag_y;
  bool found_magy_line = false;
  for (int iline = (int) vcontent.size() - 1; iline >= 0; iline--) {  // NEW FROM THE BACK
    if (aurostd::substring2bool(vcontent[iline], "magnetization")) {
      if (aurostd::substring2bool(vcontent[iline], "(y)")) {
        mline = iline + 4;
        found_magy_line = true;
        break;
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " magnetization (y) line=" << line << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " magnetization (y) mline=" << mline << endl;
  }
  if (found_magy_line) {
    for (uint iline = mline; iline < mline + natoms; iline++) {
      aurostd::string2tokens(vcontent[iline], tokens, " ");
      if (tokens.size() >= 5) {
        vmag_y.push_back(aurostd::string2utype<double>(tokens.at(tokens.size() - 1)));
      }
    }
  }
  // DX20171205 - Check for magnetization y (non-collinear) - END
  vector<double> vmag_x;  // DX20171205
  bool found_magx_line = false;
  for (int iline = (int) vcontent.size() - 1; iline >= 0; iline--) {  // NEW FROM THE BACK
    if (aurostd::substring2bool(vcontent[iline], "magnetization")) {
      if (aurostd::substring2bool(vcontent[iline], "(x)")) {
        mline = iline + 4;
        found_magx_line = true;
        break;
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " magnetization (x) line=" << line << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " magnetization (x) mline=" << mline << endl;
  }
  if (found_magx_line) {
    for (uint iline = mline; iline < mline + natoms; iline++) {
      aurostd::string2tokens(vcontent[iline], tokens, " ");
      if (tokens.size() >= 5) {
        vmag_x.push_back(aurostd::string2utype<double>(tokens.at(tokens.size() - 1))); // DX20171205 vmag to vmag_x
      }
    }
  }
  // DX20171205 - Non-collinear vs collinear - START
  if (found_magx_line && found_magy_line && found_magz_line) { // non-collinear calculation
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " non-collinear magnetization found." << endl;
    }
    if (vmag_x.size() != vmag_y.size() || vmag_x.size() != vmag_z.size()) {
      message << "Number of magnetization components (x, y, z) are not the same in OUTCAR; filename=[" << filename << "]";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
      ERROR_flag = true;
    }
    for (size_t m = 0; m < vmag_x.size(); m++) {
      xvector<double> mag_xyz;
      mag_xyz(1) = vmag_x[m];
      mag_xyz(2) = vmag_y[m];
      mag_xyz(3) = vmag_z[m];
      vmag_noncoll.push_back(mag_xyz);
    }
  } else if (found_magx_line && !found_magy_line && !found_magz_line) { // collinear calculation (x only)
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " collinear magnetization found." << endl;
    }
    for (size_t m = 0; m < vmag_x.size(); m++) {
      vmag.push_back(vmag_x[m]);
    }
  }
  // ----------------------------------------------------------------------
  // LOAD FORCES/POSITIONS
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD FORCES/POSITIONS DATA (" << time_delay(seconds) << ")" << endl;
  }
  vforces.clear();                    // QM FORCES calculation
  vpositions_cartesian.clear();                 // QM POSITIONS calculation
  for (uint i = 0; i < (uint) natoms; i++) {
    const xvector<double> force(3);
    vforces.push_back(force);
    const xvector<double> position(3);
    vpositions_cartesian.push_back(position);
  }
  for (int iline = (int) vcontent.size() - 1; iline >= 0; iline--) {
    if (aurostd::substring2bool(vcontent[iline], "TOTAL-FORCE (eV/Angst)")) {
      for (size_t iat = 0; iat < (size_t) natoms && iat < vcontent.size(); iat++) {
        aurostd::string2tokens(vcontent[iline + iat + 2], tokens, " ");
        if (tokens.size() < 6) {
          message << "Wrong number of force/positions entries in OUTCAR; line=[ " << vcontent[iline + iat + 2] << "]" << "   filename=[" << filename << "]";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
          ERROR_flag = true;
        }
        vpositions_cartesian.at(iat)[1] = aurostd::string2utype<double>(tokens.at(0));
        vpositions_cartesian.at(iat)[2] = aurostd::string2utype<double>(tokens.at(1));
        vpositions_cartesian.at(iat)[3] = aurostd::string2utype<double>(tokens.at(2));
        vforces.at(iat)[1] = aurostd::string2utype<double>(tokens.at(3));
        vforces.at(iat)[2] = aurostd::string2utype<double>(tokens.at(4));
        vforces.at(iat)[3] = aurostd::string2utype<double>(tokens.at(5));
      }
      iline = -1;
    }
  }

  // ----------------------------------------------------------------------
  // LOAD DIELECTRIC DATA (IF PRESENT)
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD DIELECTRIC DATA (" << time_delay(seconds) << ")" << endl;
  }
  int iline_freq_plasma = 0;
  int iline_df_real = 0;
  int iline_df_imag = 0;
  for (int iline = static_cast<int>(vcontent.size()) - 1; iline >= 0 && (iline_freq_plasma == 0 || iline_df_real == 0 || iline_df_imag == 0); iline--) {
    if (aurostd::substring2bool(vcontent[iline], "plasma frequency squared")) {
      if (aurostd::substring2bool(vcontent[iline], "from intraband transitions")) {
        iline_freq_plasma = iline;
      }
    } else if (aurostd::substring2bool(vcontent[iline], "REAL DIELECTRIC FUNCTION")) {
      if (aurostd::substring2bool(vcontent[iline], "density-density")) {
        iline_df_real = iline;
      }
    } else if (aurostd::substring2bool(vcontent[iline], "IMAGINARY DIELECTRIC FUNCTION")) {
      if (aurostd::substring2bool(vcontent[iline], "density-density")) {
        iline_df_imag = iline;
      }
    }
  }
  aurostd::xmatrix<double> tensor(3, 3);
  if (iline_freq_plasma) {
    for (int ir = 1; ir <= 3; ir++) {
      aurostd::string2tokens(vcontent[iline_freq_plasma + 1 + ir], tokens);
      tensor[ir][1] = std::sqrt(std::abs(aurostd::string2utype<double>(tokens[0])));
      tensor[ir][2] = std::sqrt(std::abs(aurostd::string2utype<double>(tokens[1])));
      tensor[ir][3] = std::sqrt(std::abs(aurostd::string2utype<double>(tokens[2])));
    }
    freq_plasma = tensor;
  }
  if (iline_df_real) {
    int il = 3;
    while (vcontent[iline_df_real + il].length() > 1 && vcontent[iline_df_real + il][2] != 'f') {
      aurostd::string2tokens(vcontent[iline_df_real + il], tokens);
      freq_grid.emplace_back(aurostd::string2utype<double>(tokens[0]));
      tensor[1][1] = aurostd::string2utype<double>(tokens[1]);
      tensor[2][2] = aurostd::string2utype<double>(tokens[2]);
      tensor[3][3] = aurostd::string2utype<double>(tokens[3]);
      tensor[1][2] = aurostd::string2utype<double>(tokens[4]);
      tensor[2][3] = aurostd::string2utype<double>(tokens[5]);
      tensor[1][3] = aurostd::string2utype<double>(tokens[6]);
      tensor[2][1] = tensor[1][2];
      tensor[3][2] = tensor[2][3];
      tensor[3][1] = tensor[1][3];
      dielectric_interband_real.emplace_back(tensor);
      il++;
    }
    dielectric_static = dielectric_interband_real[0];
  }
  if (iline_df_imag) {
    int il = 3;
    while (vcontent[iline_df_imag + il].length() > 1 && vcontent[iline_df_imag + il][2] != 'f') {
      aurostd::string2tokens(vcontent[iline_df_imag + il], tokens);
      tensor[1][1] = aurostd::string2utype<double>(tokens[1]);
      tensor[2][2] = aurostd::string2utype<double>(tokens[2]);
      tensor[3][3] = aurostd::string2utype<double>(tokens[3]);
      tensor[1][2] = aurostd::string2utype<double>(tokens[4]);
      tensor[2][3] = aurostd::string2utype<double>(tokens[5]);
      tensor[1][3] = aurostd::string2utype<double>(tokens[6]);
      tensor[2][1] = tensor[1][2];
      tensor[3][2] = tensor[2][3];
      tensor[3][1] = tensor[1][3];
      dielectric_interband_imag.emplace_back(tensor);
      il++;
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " stress=" << endl << stress << endl;
  }

  // ----------------------------------------------------------------------
  // ENCUT EDIFF EDIFFG POTIM TEIN TEBEG TEEND SMASS NPACO APACO PSTRESS
  // NBANDS NKPTS NSW NBLOCK KBLOCK IBRION NFREE ISIF IWAVPR ISYM TEIN TEBEG ISPIN
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD \"Electronic convergence/Ionic relaxation\" DATA (" << time_delay(seconds) << ")" << endl;
  }
  vline.clear();
  ENCUT = 0.0;
  EDIFF = 0.0;
  EDIFFG = 0.0;
  POTIM = 0.0;
  TEIN = 0.0;
  TEBEG = 0.0;
  TEEND = 0.0;
  SMASS = 0.0;
  NPACO = 0.0;
  APACO = 0.0;
  PSTRESS = 0.0;
  pressure = 0.0;     //
  NBANDS = 0;
  NKPTS = 0;
  NSW = 0;
  NBLOCK = 0;
  KBLOCK = 0;
  IBRION = 0;
  NFREE = 0;
  ISIF = 0;
  IWAVPR = 0;
  ISYM = 0;
  TEIN = 0;
  TEBEG = 0;
  ISPIN = 0; //
  total_energy_change = 0.0;
  for (size_t iline = 0; iline < vcontentRED.size(); iline++) { // CO20200106 - patching for auto-indenting
    if (aurostd::substring2bool(vcontentRED[iline], "ENCUT") && aurostd::substring2bool(vcontentRED[iline], "eV")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "EDIFF") && aurostd::substring2bool(vcontentRED[iline], "stopping-criterion for ELM")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "EDIFFG") && aurostd::substring2bool(vcontentRED[iline], "stopping-criterion for IOM")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "NSW") && aurostd::substring2bool(vcontentRED[iline], "number of steps for IOM")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "NBLOCK") && aurostd::substring2bool(vcontentRED[iline], "inner block")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "KBLOCK") && aurostd::substring2bool(vcontentRED[iline], "outer block")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "IBRION") && aurostd::substring2bool(vcontentRED[iline], "ionic relax")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "NFREE") && aurostd::substring2bool(vcontentRED[iline], "steps in history")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "ISIF") && aurostd::substring2bool(vcontentRED[iline], "stress and relaxation")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "IWAVPR") && aurostd::substring2bool(vcontentRED[iline], "prediction")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "ISYM") && aurostd::substring2bool(vcontentRED[iline], "nonsym")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "POTIM") && aurostd::substring2bool(vcontentRED[iline], "time-step")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "TEIN") && aurostd::substring2bool(vcontentRED[iline], "initial temperature")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "TEBEG") && aurostd::substring2bool(vcontentRED[iline], "temperature during run")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "TEEND") && aurostd::substring2bool(vcontentRED[iline], "temperature during run")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "SMASS") && aurostd::substring2bool(vcontentRED[iline], "Nose mass-parameter")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "NPACO") && aurostd::substring2bool(vcontentRED[iline], "distance and # of slots")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "APACO") && aurostd::substring2bool(vcontentRED[iline], "distance and # of slots")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "PSTRESS") && (aurostd::substring2bool(vcontentRED[iline], "pullay stress") || aurostd::substring2bool(vcontentRED[iline], "pulay stress"))) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "ISPIN") && aurostd::substring2bool(vcontentRED[iline], "spin polarized calculation?")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "NKPTS") && aurostd::substring2bool(vcontentRED[iline], "k-Points")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "NBANDS") && aurostd::substring2bool(vcontentRED[iline], "number of bands")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "total energy-change") && aurostd::substring2bool(vcontentRED[iline], "(2. order)")) {
      vline.push_back(vcontentRED[iline]);
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vline.size()=" << vline.size() << " (" << time_delay(seconds) << ")" << endl;
  }
  if (vline.empty()) {
    message << "Wrong number of \"Ionic relaxation\" in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  }
  for (size_t j = 0; j < vline.size(); j++) {   // to the back
    aurostd::StringSubstInPlace(vline[j], "=", " ");
    aurostd::StringSubstInPlace(vline[j], ";", " ");
    aurostd::StringSubstInPlace(vline[j], ":", " ");
    aurostd::string2tokens(vline[j], tokens, " ");
    for (size_t k = 0; k < tokens.size(); k++) {
      if (tokens[k] == "ENCUT" && k + 1 < tokens.size()) {
        ENCUT = aurostd::string2utype<double>(tokens.at(k + 1));
      }
      if (tokens[k] == "EDIFF" && k + 1 < tokens.size()) {
        EDIFF = aurostd::string2utype<double>(tokens.at(k + 1));
      }
      if (tokens[k] == "EDIFFG" && k + 1 < tokens.size()) {
        EDIFFG = aurostd::string2utype<double>(tokens.at(k + 1));
      }
      if (tokens[k] == "NSW" && k + 1 < tokens.size()) {
        NSW = aurostd::string2utype<int>(tokens.at(k + 1));
      }
      if (tokens[k] == "NBLOCK" && k + 1 < tokens.size()) {
        NBLOCK = aurostd::string2utype<int>(tokens.at(k + 1));
      }
      if (tokens[k] == "KBLOCK" && k + 1 < tokens.size()) {
        KBLOCK = aurostd::string2utype<int>(tokens.at(k + 1));
      }
      if (tokens[k] == "IBRION" && k + 1 < tokens.size()) {
        IBRION = aurostd::string2utype<int>(tokens.at(k + 1));
      }
      if (tokens[k] == "NFREE" && k + 1 < tokens.size()) {
        NFREE = aurostd::string2utype<int>(tokens.at(k + 1));
      }
      if (tokens[k] == "ISIF" && k + 1 < tokens.size()) {
        ISIF = aurostd::string2utype<int>(tokens.at(k + 1));
      }
      if (tokens[k] == "IWAVPR" && k + 1 < tokens.size()) {
        IWAVPR = aurostd::string2utype<int>(tokens.at(k + 1));
      }
      if (tokens[k] == "ISYM" && k + 1 < tokens.size()) {
        ISYM = aurostd::string2utype<int>(tokens.at(k + 1));
      }
      if (tokens[k] == "POTIM" && k + 1 < tokens.size()) {
        POTIM = aurostd::string2utype<double>(tokens.at(k + 1));
      }
      if (tokens[k] == "TEIN" && k + 1 < tokens.size()) {
        TEIN = aurostd::string2utype<double>(tokens.at(k + 1));
      }
      if (tokens[k] == "TEBEG" && k + 1 < tokens.size()) {
        TEBEG = aurostd::string2utype<double>(tokens.at(k + 1));
      }
      if (tokens[k] == "TEEND" && k + 1 < tokens.size()) {
        TEEND = aurostd::string2utype<double>(tokens.at(k + 1));
      }
      if (tokens[k] == "SMASS" && k + 1 < tokens.size()) {
        SMASS = aurostd::string2utype<double>(tokens.at(k + 1));
      }
      if (tokens[k] == "NPACO" && k + 1 < tokens.size()) {
        NPACO = aurostd::string2utype<double>(tokens.at(k + 1));
      }
      if (tokens[k] == "APACO" && k + 1 < tokens.size()) {
        APACO = aurostd::string2utype<double>(tokens.at(k + 1));
      }
      if (tokens[k] == "PSTRESS" && k + 1 < tokens.size()) {
        PSTRESS = aurostd::string2utype<double>(tokens.at(k + 1));
        pressure = PSTRESS;
      }
      if (tokens[k] == "NBANDS" && k + 1 < tokens.size()) {
        NBANDS = aurostd::string2utype<int>(tokens.at(k + 1));
      }
      if (tokens[k] == "NKPTS" && k + 1 < tokens.size()) {
        NKPTS = aurostd::string2utype<int>(tokens.at(k + 1));
      }
      if (tokens[k] == "ISPIN" && k + 1 < tokens.size()) {
        ISPIN = aurostd::string2utype<int>(tokens.at(k + 1));
      }
      if (tokens[k] == "order)" && k + 1 < tokens.size()) {
        total_energy_change = aurostd::string2utype<double>(tokens.at(k + 1));
      }
    }
  }
  if (!vline.empty()) { // CO
    for (size_t j = vline.size() - 1; j > 0; j--) {   // to the front
      aurostd::StringSubstInPlace(vline.at(j), "=", " ");
      aurostd::StringSubstInPlace(vline.at(j), ";", " ");
      aurostd::StringSubstInPlace(vline.at(j), ":", " ");
      aurostd::string2tokens(vline.at(j), tokens, " ");
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ENCUT=" << ENCUT << endl;
    cerr << __AFLOW_FUNC__ << " EDIFF=" << EDIFF << endl;
    cerr << __AFLOW_FUNC__ << " EDIFFG=" << EDIFFG << endl;
    cerr << __AFLOW_FUNC__ << " NSW=" << NSW << endl;
    cerr << __AFLOW_FUNC__ << " NBLOCK=" << NBLOCK << endl;
    cerr << __AFLOW_FUNC__ << " KBLOCK=" << KBLOCK << endl;
    cerr << __AFLOW_FUNC__ << " IBRION=" << IBRION << endl;
    cerr << __AFLOW_FUNC__ << " NFREE=" << NFREE << endl;
    cerr << __AFLOW_FUNC__ << " ISIF=" << ISIF << endl;
    cerr << __AFLOW_FUNC__ << " IWAVPR=" << IWAVPR << endl;
    cerr << __AFLOW_FUNC__ << " ISYM=" << ISYM << endl;
    cerr << __AFLOW_FUNC__ << " POTIM=" << POTIM << endl;
    cerr << __AFLOW_FUNC__ << " TEIN=" << TEIN << endl;
    cerr << __AFLOW_FUNC__ << " TEBEG=" << TEBEG << endl;
    cerr << __AFLOW_FUNC__ << " TEEND=" << TEEND << endl;
    cerr << __AFLOW_FUNC__ << " SMASS=" << SMASS << endl;
    cerr << __AFLOW_FUNC__ << " NPACO=" << NPACO << endl;
    cerr << __AFLOW_FUNC__ << " APACO=" << APACO << endl;
    cerr << __AFLOW_FUNC__ << " PSTRESS=" << PSTRESS << endl;
    cerr << __AFLOW_FUNC__ << " pressure(PSTRESS)=" << pressure << endl;
    cerr << __AFLOW_FUNC__ << " NBANDS=" << NBANDS << endl;
    cerr << __AFLOW_FUNC__ << " NKPTS=" << NKPTS << endl;
    cerr << __AFLOW_FUNC__ << " ISPIN=" << ISPIN << endl;
    cerr << __AFLOW_FUNC__ << " total_energy_change=" << total_energy_change << endl;
  }

  // if pressure correct enthalpy
  // ----------------------------------------------------------------------
  // ENMAX ENMIN POMASS ZVAL EATOM RCORE RWIGS EAUG RMAX (there are many)
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD \"ENMAX ENMIN POMASS ZVAL EATOM RCORE RWIGS EAUG RMAX\" DATA (" << time_delay(seconds) << ")" << endl;
  }
  vline.clear();
  vENMAX.clear();
  vENMIN.clear();
  vPOMASS.clear();
  vZVAL.clear();
  vEATOM.clear();
  vRCORE.clear();
  vRWIGS.clear();
  vEAUG.clear();
  vRAUG.clear();
  vRMAX.clear();
  for (size_t iline = 0; iline < vcontentRED.size(); iline++) {
    aurostd::StringSubstInPlace(vcontentRED[iline], "EMMIN", "ENMIN"); // Kresseeeeeeeeeee
    if (aurostd::substring2bool(vcontentRED[iline], "ENMAX") && aurostd::substring2bool(vcontentRED[iline], "ENMIN")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "POMASS") && aurostd::substring2bool(vcontentRED[iline], "ZVAL")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "EATOM") && aurostd::substring2bool(vcontentRED[iline], "eV")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "RCORE") && aurostd::substring2bool(vcontentRED[iline], "radius")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "RWIGS") && aurostd::substring2bool(vcontentRED[iline], "radius")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "EAUG")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "RAUG") && aurostd::substring2bool(vcontentRED[iline], "sphere")) {
      vline.push_back(vcontentRED[iline]);
    }
    if (aurostd::substring2bool(vcontentRED[iline], "RMAX") && aurostd::substring2bool(vcontentRED[iline], "radius")) {
      vline.push_back(vcontentRED[iline]);
    }
    // AS20200528
    if (aurostd::substring2bool(vcontentRED[iline], "NELECT") && aurostd::substring2bool(vcontentRED[iline], "total number of electrons")) {
      vline.push_back(vcontentRED[iline]);
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vline.size()=" << vline.size() << " (" << time_delay(seconds) << ")" << endl;
  }
  if (vline.empty()) {
    message << "Wrong number of \"ENMAX ENMIN POMASS ZVAL EATOM RCORE RWIGS EAUG RAUG RMAX\" in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  }
  for (size_t j = 0; j < vline.size(); j++) {
    aurostd::StringSubstInPlace(vline[j], "=", " ");
    aurostd::StringSubstInPlace(vline[j], ";", " ");
    aurostd::string2tokens(vline[j], tokens, " ");
    for (size_t k = 0; k < tokens.size(); k++) {
      if (tokens[k] == "ENMAX" && k + 1 < tokens.size()) {
        vENMAX.push_back(aurostd::string2utype<double>(tokens.at(k + 1)));
      }
      if (tokens[k] == "ENMIN" && k + 1 < tokens.size()) {
        vENMIN.push_back(aurostd::string2utype<double>(tokens.at(k + 1)));
      }
      if (tokens[k] == "POMASS" && k + 1 < tokens.size()) {
        vPOMASS.push_back(aurostd::string2utype<double>(tokens.at(k + 1)));
      }
      if (tokens[k] == "ZVAL" && k + 1 < tokens.size()) {
        vZVAL.push_back(aurostd::string2utype<double>(tokens.at(k + 1)));
      }
      if (tokens[k] == "EATOM" && k + 1 < tokens.size()) {
        vEATOM.push_back(aurostd::string2utype<double>(tokens.at(k + 1)));
      }
      if (tokens[k] == "RCORE" && k + 1 < tokens.size()) {
        vRCORE.push_back(aurostd::string2utype<double>(tokens.at(k + 1)));
      }
      if (tokens[k] == "RWIGS" && k == 0 && k + 1 < tokens.size()) {
        vRWIGS.push_back(aurostd::string2utype<double>(tokens.at(k + 1))); // pick the 1st
      }
      if (tokens[k] == "EAUG" && k + 1 < tokens.size()) {
        vEAUG.push_back(aurostd::string2utype<double>(tokens.at(k + 1)));
      }
      if (tokens[k] == "RAUG" && k + 1 < tokens.size()) {
        vRAUG.push_back(aurostd::string2utype<double>(tokens.at(k + 1)));
      }
      if (tokens[k] == "RMAX" && k + 1 < tokens.size()) {
        vRMAX.push_back(aurostd::string2utype<double>(tokens.at(k + 1)));
      }
      if (tokens[k] == "NELECT" && k + 1 < tokens.size()) {
        nelectrons = aurostd::string2utype<int>(tokens.at(k + 1));// AS20200528
      }
    }
  }

  if (vENMIN.size() < vENMAX.size()) {
    vENMIN.push_back(0.0);  // ENMIN MIGHT NOT BE THERE
  }
  if (vRCORE.size() < vENMAX.size()) {
    vRCORE.push_back(0.0);  // RCORE MIGHT NOT BE THERE
  }
  if (vRWIGS.size() < vENMAX.size()) {
    vRWIGS.push_back(0.0);  // RWIGS MIGHT NOT BE THERE
  }
  if (vEAUG.size() < vENMAX.size()) {
    vEAUG.push_back(0.0);    // EAUG MIGHT NOT BE THERE
  }
  if (vRAUG.size() < vENMAX.size()) {
    vRAUG.push_back(0.0);    // RAUG MIGHT NOT BE THERE
  }
  if (vRMAX.size() < vENMAX.size()) {
    vRMAX.push_back(0.0);    // RMAX MIGHT NOT BE THERE
  }

  ENMAX = aurostd::max(vENMAX);
  ENMIN = aurostd::min(vENMIN);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ENMAX=" << ENMAX << " vENMAX.size()=" << vENMAX.size() << ": ";
    for (size_t i = 0; i < vENMAX.size(); i++) {
      cerr << vENMAX[i] << " ";
    }
    cerr << " " << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ENMIN=" << ENMIN << " vENMIN.size()=" << vENMIN.size() << ": ";
    for (size_t i = 0; i < vENMIN.size(); i++) {
      cerr << vENMIN[i] << " ";
    }
    cerr << " " << endl;
  }

  POMASS_sum = aurostd::sum(vPOMASS);
  POMASS_min = aurostd::min(vPOMASS);
  POMASS_max = aurostd::max(vPOMASS);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " POMASS_sum=" << POMASS_sum << " POMASS_min=" << POMASS_min << " POMASS_max=" << POMASS_max << " vPOMASS.size()=" << vPOMASS.size() << ": ";
    for (size_t i = 0; i < vPOMASS.size(); i++) {
      cerr << vPOMASS[i] << " ";
    }
    cerr << " " << endl;
  }

  ZVAL_sum = aurostd::sum(vZVAL);
  ZVAL_min = aurostd::min(vZVAL);
  ZVAL_max = aurostd::max(vZVAL);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ZVAL_sum=" << ZVAL_sum << " ZVAL_min=" << ZVAL_min << " ZVAL_max=" << ZVAL_max << " vZVAL.size()=" << vZVAL.size() << ": ";
    for (size_t i = 0; i < vZVAL.size(); i++) {
      cerr << vZVAL[i] << " ";
    }
    cerr << " " << endl;
  }

  EATOM_min = aurostd::min(vEATOM);
  EATOM_max = aurostd::max(vEATOM);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " EATOM_min=" << EATOM_min << " EATOM_max=" << EATOM_max << " vEATOM.size()=" << vEATOM.size() << ": ";
    for (size_t i = 0; i < vEATOM.size(); i++) {
      cerr << vEATOM[i] << " ";
    }
    cerr << " " << endl;
  }

  RCORE_min = aurostd::min(vRCORE);
  RCORE_max = aurostd::max(vRCORE);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " RCORE_min=" << RCORE_min << " RCORE_max=" << RCORE_max << " vRCORE.size()=" << vRCORE.size() << ": ";
    for (size_t i = 0; i < vRCORE.size(); i++) {
      cerr << vRCORE[i] << " ";
    }
    cerr << " " << endl;
  }

  RWIGS_min = aurostd::min(vRWIGS);
  RWIGS_max = aurostd::max(vRWIGS);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " RWIGS_min=" << RWIGS_min << " RWIGS_max=" << RWIGS_max << " vRWIGS.size()=" << vRWIGS.size() << ": ";
    for (size_t i = 0; i < vRWIGS.size(); i++) {
      cerr << vRWIGS[i] << " ";
    }
    cerr << " " << endl;
  }

  EAUG_min = aurostd::min(vEAUG);
  EAUG_max = aurostd::max(vEAUG);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " EAUG_min=" << EAUG_min << " EAUG_max=" << EAUG_max << " vEAUG.size()=" << vEAUG.size() << ": ";
    for (size_t i = 0; i < vEAUG.size(); i++) {
      cerr << vEAUG[i] << " ";
    }
    cerr << " " << endl;
  }

  RAUG_min = aurostd::min(vRAUG);
  RAUG_max = aurostd::max(vRAUG);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " RAUG_min=" << RAUG_min << " RAUG_max=" << RAUG_max << " vRAUG.size()=" << vRAUG.size() << ": ";
    for (size_t i = 0; i < vRAUG.size(); i++) {
      cerr << vRAUG[i] << " ";
    }
    cerr << " " << endl;
  }

  RMAX_min = aurostd::min(vRMAX);
  RMAX_max = aurostd::max(vRMAX);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " RMAX_min=" << RMAX_min << " RMAX_max=" << RMAX_max << " vRMAX.size()=" << vRMAX.size() << ": ";
    for (size_t i = 0; i < vRMAX.size(); i++) {
      cerr << vRMAX[i] << " ";
    }
    cerr << " " << endl;
  }

  // ----------------------------------------------------------------------
  // KINETIC AND METAGGA DATA
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD KINETIC AND METAGGA DATA (" << time_delay(seconds) << ")" << endl;
  }
  isKIN = false;
  isMETAGGA = false;
  METAGGA = "";
  vline.clear();
  for (size_t iline = 0; iline < vcontent.size(); iline++) {
    if (aurostd::substring2bool(vcontent[iline], "partial kinetic energy density read in")) {
      isKIN = true;
    }
    if (aurostd::substring2bool(vcontent[iline], "METAGGA")) {
      vline.push_back(vcontent[iline]);
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vline.size()=" << vline.size() << endl;
  }

  if (!vline.empty()) {
    aurostd::string2tokens(vline.at(0), tokens, "=");
    if (tokens.size() > 1) {
      aurostd::string2tokens(string(tokens.at(1)), tokens);
      if (!tokens.empty()) {
        if (tokens.at(0) != "F") {
          isMETAGGA = true;
          METAGGA = tokens.at(0);
        }
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " isKIN=" << isKIN << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " isMETAGGA=" << isMETAGGA << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " METAGGA=" << METAGGA << endl;
  }

  // ----------------------------------------------------------------------
  // PSEUDOPOTENTIAL DATA
  // LDEBUG=1;
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD PSEUDOPOTENTIAL DATA (" << time_delay(seconds) << ")" << endl;
  }

  vline.clear();
  vTITEL.clear();
  vLEXCH.clear();
  for (size_t iline = 0; iline < vcontentRED.size(); iline++) {
    if (aurostd::substring2bool(vcontentRED[iline], "=")) {
      if (aurostd::substring2bool(vcontentRED[iline], "TITEL")) {
        aurostd::string2tokens(vcontentRED[iline], tokens, "=");
        if (tokens.size() > 1) {
          vTITEL.push_back(tokens.at(1));
        }
      }
      if (aurostd::substring2bool(vcontentRED[iline], "LEXCH")) {
        if (!aurostd::substring2bool(vcontentRED[iline], "exchange")) {
          aurostd::string2tokens(vcontentRED[iline], tokens, "=");
          if (tokens.size() > 1) {
            vLEXCH.push_back(tokens.at(1));
          }
        }
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vTITEL.size()=" << vTITEL.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vLEXCH.size()=" << vLEXCH.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vEATOM.size()=" << vEATOM.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vRMAX.size()=" << vRMAX.size() << endl;
  }
  if (vTITEL.empty()) {
    message << "Wrong number of pseudopotentials (TITEL) in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  } // CO20200106 - patching for auto-indenting
  if (vLEXCH.empty()) {
    message << "Wrong number of pseudopotentials (LEXCH) in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  } // CO20200106 - patching for auto-indenting
  if (vEATOM.empty()) {
    message << "Wrong number of pseudopotentials (EATOM) in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  } // CO20200106 - patching for auto-indenting
  if (vRMAX.empty()) {
    message << "Wrong number of pseudopotentials (RMAX) in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  }   // CO20200106 - patching for auto-indenting
  if (vTITEL.size() != vLEXCH.size()) {
    message << "Wrong number of pseudopotentials (TITEL/LEXCH) in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  } // CO20200106 - patching for auto-indenting
  if (vLEXCH.size() != vEATOM.size()) {
    message << "Wrong number of pseudopotentials (LEXCH/EATOM) in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  } // CO20200106 - patching for auto-indenting
  if (vEATOM.size() != vRMAX.size()) {
    message << "Wrong number of pseudopotentials (EATOM/RMAX) in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  }   // CO20200106 - patching for auto-indenting

  for (size_t j = 0; j < vTITEL.size(); j++) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " SPECIES(" << j << ") " << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vTITEL[j]=" << vTITEL[j] << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vLEXCH.at(j)=" << vLEXCH.at(j) << endl;
    }
    aurostd::string2tokens(vTITEL[j], tokens, " ");
    pp_type = tokens.at(0);
    if (pp_type == "US" or pp_type == "NC") {
      pp_type = "GGA";
      for (size_t iTITEL = 0; iTITEL < vcontent.size(); iTITEL++) {
        if (aurostd::substring2bool(vcontent[iTITEL], "GGA")) {
          if (aurostd::substring2bool(vcontent[iTITEL], "eV")) {
            pp_type = "LDA";
          }
        }
      }
    }
    if (pp_type == "PAW") {
      pp_type = "PAW_LDA"; // cerr << __AFLOW_FUNC__ << " PAW_LDA" << endl;
    }
    if (isKIN) {
      pp_type += "_KIN";
    }
    if (isMETAGGA) {
      pp_type += ":" + METAGGA;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " pp_type=" << pp_type << endl;
    }
    species.push_back(aurostd::VASP_PseudoPotential_CleanName(tokens.at(1)));
    species_Z.push_back(xelement::symbol2Z(aurostd::VASP_PseudoPotential_CleanName(tokens.at(1))));
    species_pp.push_back(tokens.at(1));
    species_pp_type.push_back(pp_type);
    if (pp_type == "LDA" && tokens.size() < 3) {
      tokens.push_back(DEFAULT_VASP_POTCAR_DATE_POT_LDA);
    }
    if (pp_type == "GGA" && tokens.size() < 3) {
      tokens.push_back(DEFAULT_VASP_POTCAR_DATE_POT_GGA);
    }
    species_pp_version.push_back(tokens.at(1) + ":" + pp_type + ":" + tokens.at(2));

    vTITEL[j] = aurostd::RemoveWhiteSpaces(vTITEL[j]);
    vLEXCH.at(j) = aurostd::RemoveWhiteSpaces(vLEXCH.at(j));
    xPOTCAR xPOT(xPOTCAR_Finder(species_pp_AUID, species_pp_AUID_collisions, vTITEL[j], vLEXCH.at(j), vEATOM.at(j), vRMAX.at(j), LDEBUG)); // FIXES species_pp_AUID,species_pp_AUID_collisions
    species_pp_groundstate_energy.push_back(xPOT.species_pp_groundstate_energy.at(0));
    species_pp_groundstate_structure.push_back(xPOT.species_pp_groundstate_structure.at(0));

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " SPECIES(" << j << ") [pp, type, version, AUID] = " << species.at(species.size() - 1) << " [" << species_Z.at(species_Z.size() - 1) << ", " << species_pp.at(species_pp.size() - 1)
           << ", " << species_pp_type.at(species_pp_type.size() - 1) << ", " << species_pp_version.at(species_pp_version.size() - 1) << ", " << species_pp_AUID.at(species_pp_AUID.size() - 1) << ", "
           << species_pp_groundstate_energy.at(species_pp_groundstate_energy.size() - 1) << ", " << species_pp_groundstate_structure.at(species_pp_groundstate_structure.size() - 1) << "]" << endl;
    }
  }
  // CO20210213 - check types are all the same, if not issue warning/error (mixing is not advisable)
  for (size_t i = 0; i < species_pp_type.size(); i++) {
    if (species_pp_type[i] != pp_type) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Mismatch in species_pp_types (" + species_pp_type[i] + " vs. " + pp_type + ")", *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    }
  }

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " PSEUDOPOTENTIAL type = " << pp_type << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species.size()=" << species.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species_Z.size()=" << species_Z.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species_pp.size()=" << species_pp.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species_pp_type.size()=" << species_pp_type.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species_pp_version.size()=" << species_pp_version.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species_pp_AUID.size()=" << species_pp_AUID.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species_pp_AUID_collisions.size()=" << species_pp_AUID_collisions.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species_pp_groundstate_energy.size()=" << species_pp_groundstate_energy.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species_pp_groundstate_structure.size()=" << species_pp_groundstate_structure.size() << endl;
  }
  if (!species_pp_AUID_collisions.empty()) {
    message << "COLLISION species_pp_AUID_collisions.size()=" << species_pp_AUID_collisions.size();
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  }

  // ----------------------------------------------------------------------
  // LDAU DATA
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD LDAU DATA (" << time_delay(seconds) << ")" << endl;
  }
  for (size_t j = 0; j < species.size(); j++) {
    species_pp_vLDAU.emplace_back();  // make space for LDAU
  }

  vline.clear();
  for (size_t iline = 0; iline < vcontent.size(); iline++) {
    if (!aurostd::substring2bool(vcontent[iline], "POSCAR")) {
      if (!aurostd::substring2bool(vcontent[iline], "SYSTEM")) {
        if (aurostd::substring2bool(vcontent[iline], "LDAU")) {
          if (!aurostd::substring2bool(vcontent[iline], "SYSTEM")) {
            vline.push_back(vcontent[iline]);
          }
        }
      }
    }
  }

  int LDAUT = 0;
  stringstream sdata_ldau;
  if (!vline.empty()) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " LDAU calculation in OUTCAR" << endl;
    }
    vector<int> vLDAUL;
    vector<double> vLDAUU;
    vector<double> vLDAUJ;
    for (size_t j = 0; j < vline.size(); j++) {
      aurostd::string2tokens(vline[j], tokens, "=");
      vline[j] = tokens.at(1);
    }
    for (size_t j = 0; j < vline.size(); j++) {
      aurostd::string2tokens(vline[j], tokens, " ");
      if (j == 0) {
        LDAUT = aurostd::string2utype<int>(vline[j]);
      }
      if (j == 1) {
        for (size_t i = 0; i < tokens.size(); i++) {
          vLDAUL.push_back(aurostd::string2utype<int>(tokens[i]));
        }
      }
      if (j == 2) {
        for (size_t i = 0; i < tokens.size(); i++) {
          vLDAUU.push_back((100 * aurostd::string2utype<double>(tokens[i])) / 100.0);
        }
      }
      if (j == 3) {
        for (size_t i = 0; i < tokens.size(); i++) {
          vLDAUJ.push_back((100 * aurostd::string2utype<double>(tokens[i])) / 100.0);
        }
      }
    }
    if (species_pp_vLDAU.size() != species.size()) {
      message << "species_pp_vLDAU.size()[" << species_pp_vLDAU.size() << "] != species.size()[" << species.size() << "]" << "   filename=[" << filename << "]";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
      ERROR_flag = true;
    }
    if (species_pp_vLDAU.size() != vLDAUL.size()) {
      message << "species_pp_vLDAU.size()[" << species_pp_vLDAU.size() << "] != vLDAUL.size()[" << vLDAUL.size() << "]" << "   filename=[" << filename << "]";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
      ERROR_flag = true;
    }
    if (species_pp_vLDAU.size() != vLDAUU.size()) {
      message << "species_pp_vLDAU.size()[" << species_pp_vLDAU.size() << "] != vLDAUU.size()[" << vLDAUU.size() << "]" << "   filename=[" << filename << "]";
      ERROR_flag = true;
    }
    if (species_pp_vLDAU.size() != vLDAUJ.size()) {
      message << "species_pp_vLDAU.size()[" << species_pp_vLDAU.size() << "] != vLDAUJ.size()[" << vLDAUJ.size() << "]" << "   filename=[" << filename << "]";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
      ERROR_flag = true;
    }
    for (size_t j = 0; j < species.size(); j++) {
      species_pp_vLDAU.at(j).push_back(LDAUT);
      species_pp_vLDAU.at(j).push_back(vLDAUL.at(j));
      species_pp_vLDAU.at(j).push_back(vLDAUU.at(j));
      species_pp_vLDAU.at(j).push_back(vLDAUJ.at(j));
    }

    if (LDEBUG) {
      //
      cerr << __AFLOW_FUNC__ << " LDA_type=" << LDAUT << endl;
      //
      cerr << __AFLOW_FUNC__ << " LDAU_L=";
      for (size_t i = 0; i < vLDAUL.size(); i++) {
        cerr << vLDAUL[i] << ((i < vLDAUL.size() - 1) ? "," : "");
      }
      cerr << endl;
      //
      cerr << __AFLOW_FUNC__ << " LDAU_U=";
      for (size_t i = 0; i < vLDAUU.size(); i++) {
        cerr << vLDAUU[i] << ((i < vLDAUU.size() - 1) ? "," : "");
      }
      cerr << endl;
      //
      cerr << __AFLOW_FUNC__ << " LDAU_J=";
      for (size_t i = 0; i < vLDAUJ.size(); i++) {
        cerr << vLDAUJ[i] << ((i < vLDAUJ.size() - 1) ? "," : "");
      }
      cerr << endl;
    }
    sdata_ldau << aurostd::utype2string(LDAUT) + ";";
    for (size_t i = 0; i < vLDAUL.size(); i++) {
      sdata_ldau << vLDAUL[i] << ((i < vLDAUL.size() - 1) ? "," : "");
    }
    sdata_ldau << ";";
    for (size_t i = 0; i < vLDAUU.size(); i++) {
      sdata_ldau << vLDAUU[i] << ((i < vLDAUU.size() - 1) ? "," : "");
    }
    sdata_ldau << ";";
    for (size_t i = 0; i < vLDAUJ.size(); i++) {
      sdata_ldau << vLDAUJ[i] << ((i < vLDAUJ.size() - 1) ? "," : "");
    }
    // sdata_ldau << ";";
    string_LDAU = sdata_ldau.str();
  } else {
    // CO20210713 - push back 0 for no +U
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " LDAU calculation NOT FOUND in OUTCAR" << endl;
    }
    for (size_t j = 0; j < species.size(); j++) {
      species_pp_vLDAU.at(j).push_back(LDAUT);
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " LDA_type=" << LDAUT << endl;
    }
    //[CO+ME20210713 - keep legacy behavior, only print when non-zero]sdata_ldau << aurostd::utype2string(LDAUT);
    //[CO+ME20210713 - keep legacy behavior, only print when non-zero]string_LDAU=sdata_ldau.str();
  }

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " string_LDAU=" << string_LDAU << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species_pp_vLDAU.size()=" << species_pp_vLDAU.size() << endl;
  }

  // ----------------------------------------------------------------------
  // DOS related values:
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD DOS related values DATA (" << time_delay(seconds) << ")" << endl;
  }
  EMIN = 0.0;
  EMAX = 0.0;
  SIGMA = 0.0;
  ISMEAR = 0;  // for aflowlib_libraries.cpp
  vline.clear();
  for (size_t iline = vcontentRED.size() - 1; iline < vcontentRED.size(); iline--)  // CO20200404
    //   for(size_t iline=0;iline<vcontentRED.size();iline++)
  { // CO20200106 - patching for auto-indenting
    if (aurostd::substring2bool(vcontentRED.at(iline), "EMIN") && aurostd::substring2bool(vcontentRED.at(iline), "energy-range for DOS")) {
      vline.push_back(vcontentRED.at(iline));
    }
    if (aurostd::substring2bool(vcontentRED.at(iline), "EMAX") && aurostd::substring2bool(vcontentRED.at(iline), "energy-range for DOS")) {
      vline.push_back(vcontentRED.at(iline));
    }
    if (aurostd::substring2bool(vcontentRED.at(iline), "ISMEAR") && aurostd::substring2bool(vcontentRED.at(iline), "broadening in eV")) {
      vline.push_back(vcontentRED.at(iline));
    }
    if (aurostd::substring2bool(vcontentRED.at(iline), "SIGMA") && aurostd::substring2bool(vcontentRED.at(iline), "broadening in eV")) {
      vline.push_back(vcontentRED.at(iline));
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vline.size()=" << vline.size() << endl;
  }
  if (vline.empty()) {
    message << "Wrong number of \"DOS related values\" in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  }
  for (size_t j = 0; j < vline.size(); j++) {
    //   if(LDEBUG) cerr << __AFLOW_FUNC__ << " vline.at(" << j << ")=" << vline[j] << endl;
    aurostd::StringSubstInPlace(vline[j], "=", " ");
    aurostd::StringSubstInPlace(vline[j], ";", " ");
    //    if(LDEBUG) cerr << vline[j] << endl;
    aurostd::string2tokens(vline[j], tokens, " ");
    for (size_t k = 0; k < tokens.size(); k++) {
      if (tokens[k] == "EMIN" && k + 1 < tokens.size()) {
        EMIN = aurostd::string2utype<double>(tokens.at(k + 1));
      }
      if (tokens[k] == "EMAX" && k + 1 < tokens.size()) {
        EMAX = aurostd::string2utype<double>(tokens.at(k + 1));
      }
      if (tokens[k] == "ISMEAR" && k + 1 < tokens.size()) {
        ISMEAR = aurostd::string2utype<int>(tokens.at(k + 1));
      }
      if (tokens[k] == "SIGMA" && k + 1 < tokens.size()) {
        SIGMA = aurostd::string2utype<double>(tokens.at(k + 1));
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " EMIN=" << EMIN << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " EMAX=" << EMAX << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ISMEAR=" << ISMEAR << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " SIGMA=" << SIGMA << endl;
  }

  // ----------------------------------------------------------------------
  // Electronic relaxation
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD Electronic relaxation DATA (" << time_delay(seconds) << ")" << endl;
  }
  IALGO = 0;                      // for aflowlib_libraries.cpp
  LDIAG = "";                     // for aflowlib_libraries.cpp
  IMIX = 0;
  INIMIX = 0;
  MIXPRE = 0;     // for aflowlib_libraries.cpp
  AMIX = 0.0;
  BMIX = 0.0;
  AMIX_MAG = 0.0;
  BMIX_MAG = 0.0;
  AMIN = 0.0;
  WC = 0.0;  // for aflowlib_libraries.cpp
  vline.clear();
  for (size_t iline = vcontentRED.size() - 1; iline < vcontentRED.size(); iline--)  // CO20200404  // DOWN
    //  for(size_t iline=0;iline<vcontentRED.size();iline++)  // UP
  { // CO20200106 - patching for auto-indenting
    if (aurostd::substring2bool(vcontentRED.at(iline), "IALGO") && aurostd::substring2bool(vcontentRED.at(iline), "algorithm")) {
      vline.push_back(vcontentRED.at(iline));
    }
    if (aurostd::substring2bool(vcontentRED.at(iline), "LDIAG") && aurostd::substring2bool(vcontentRED.at(iline), "sub-space diagonalisation")) {
      vline.push_back(vcontentRED.at(iline));
    }
    if (aurostd::substring2bool(vcontentRED.at(iline), "IMIX") && aurostd::substring2bool(vcontentRED.at(iline), "mixing-type and parameters")) {
      vline.push_back(vcontentRED.at(iline));
    }
    if (aurostd::substring2bool(vcontentRED.at(iline), "AMIX") && aurostd::substring2bool(vcontentRED.at(iline), "BMIX")) {
      vline.push_back(vcontentRED.at(iline));
    }
    if (aurostd::substring2bool(vcontentRED.at(iline), "AMIX_MAG") && aurostd::substring2bool(vcontentRED.at(iline), "BMIX_MAG")) {
      vline.push_back(vcontentRED.at(iline));
    }
    if (aurostd::substring2bool(vcontentRED.at(iline), "AMIN")) {
      vline.push_back(vcontentRED.at(iline));
    }
    if (aurostd::substring2bool(vcontentRED.at(iline), "WC") && aurostd::substring2bool(vcontentRED.at(iline), "INIMIX") && aurostd::substring2bool(vcontentRED.at(iline), "MIXPRE")) {
      vline.push_back(vcontentRED.at(iline));
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vline.size()=" << vline.size() << endl;
  }
  if (vline.empty()) {
    message << "Wrong number of \"Electronic relaxation\" in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  }
  for (size_t j = 0; j < vline.size(); j++) {
    //  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vline.at(" << j << ")=" << vline[j] << endl;
    aurostd::StringSubstInPlace(vline[j], "=", " ");
    aurostd::StringSubstInPlace(vline[j], ";", " ");
    //    if(LDEBUG) cerr << vline[j] << endl;
    aurostd::string2tokens(vline[j], tokens, " ");
    for (size_t k = 0; k < tokens.size(); k++) {
      if (tokens[k] == "IALGO" && k + 1 < tokens.size()) {
        IALGO = aurostd::string2utype<int>(tokens.at(k + 1));
      }
      if (tokens[k] == "LDIAG" && k + 1 < tokens.size()) {
        LDIAG = (tokens.at(k + 1));
      }
      if (tokens[k] == "IMIX" && k + 1 < tokens.size()) {
        IMIX = aurostd::string2utype<int>(tokens.at(k + 1));
      }
      if (tokens[k] == "AMIX" && k + 1 < tokens.size()) {
        AMIX = aurostd::string2utype<double>(tokens.at(k + 1));
      }
      if (tokens[k] == "BMIX" && k + 1 < tokens.size()) {
        BMIX = aurostd::string2utype<double>(tokens.at(k + 1));
      }
      if (tokens[k] == "AMIX_MAG" && k + 1 < tokens.size()) {
        AMIX_MAG = aurostd::string2utype<double>(tokens.at(k + 1));
      }
      if (tokens[k] == "BMIX_MAG" && k + 1 < tokens.size()) {
        BMIX_MAG = aurostd::string2utype<double>(tokens.at(k + 1));
      }
      if (tokens[k] == "AMIN" && k + 1 < tokens.size()) {
        AMIN = aurostd::string2utype<double>(tokens.at(k + 1));
      }
      if (tokens[k] == "WC" && k + 1 < tokens.size()) {
        WC = aurostd::string2utype<double>(tokens.at(k + 1));
      }
      if (tokens[k] == "INIMIX" && k + 1 < tokens.size()) {
        INIMIX = aurostd::string2utype<int>(tokens.at(k + 1));
      }
      if (tokens[k] == "MIXPRE" && k + 1 < tokens.size()) {
        MIXPRE = aurostd::string2utype<int>(tokens.at(k + 1));
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " IALGO=" << IALGO << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LDIAG=" << LDIAG << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " IMIX=" << IMIX << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " AMIX=" << AMIX << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BMIX=" << BMIX << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " AMIX_MAG=" << AMIX_MAG << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BMIX_MAG=" << BMIX_MAG << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " AMIN=" << AMIN << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " WC=" << WC << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " INIMIX=" << INIMIX << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " MIXPRE=" << MIXPRE << endl;
  }

  // ----------------------------------------------------------------------
  // Intra band minimization
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD Intra band minimization DATA (" << time_delay(seconds) << ")" << endl;
  }
  WEIMIN = 0.0;
  EBREAK = 0.0;
  DEPER = 0.0;
  TIME = 0.0;  // for aflowlib_libraries.cpp
  vline.clear();
  for (size_t iline = vcontentRED.size() - 1; iline < vcontentRED.size(); iline--)  // CO20200404
    //   for(size_t iline=0;iline<vcontentRED.size();iline++)
  { // CO20200106 - patching for auto-indenting
    if (aurostd::substring2bool(vcontentRED.at(iline), "WEIMIN") && aurostd::substring2bool(vcontentRED.at(iline), "energy-eigenvalue tresh-hold")) {
      vline.push_back(vcontentRED.at(iline));
    }
    if (aurostd::substring2bool(vcontentRED.at(iline), "EBREAK") && aurostd::substring2bool(vcontentRED.at(iline), "absolut break condition")) {
      vline.push_back(vcontentRED.at(iline));
    }
    if (aurostd::substring2bool(vcontentRED.at(iline), "DEPER") && aurostd::substring2bool(vcontentRED.at(iline), "relativ break condition")) {
      vline.push_back(vcontentRED.at(iline));
    }
    if (aurostd::substring2bool(vcontentRED.at(iline), "TIME") && aurostd::substring2bool(vcontentRED.at(iline), "timestep for ELM")) {
      vline.push_back(vcontentRED.at(iline));
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vline.size()=" << vline.size() << endl;
  }
  if (vline.empty()) {
    message << "Wrong number of \"Intra band minimization\" in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  }
  for (size_t j = 0; j < vline.size(); j++) {
    //  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vline.at(" << j << ")=" << vline[j] << endl;
    aurostd::StringSubstInPlace(vline[j], "=", " ");
    aurostd::StringSubstInPlace(vline[j], ";", " ");
    //    if(LDEBUG) cerr << vline[j] << endl;
    aurostd::string2tokens(vline[j], tokens, " ");
    for (size_t k = 0; k < tokens.size(); k++) {
      if (tokens[k] == "WEIMIN" && k + 1 < tokens.size()) {
        WEIMIN = aurostd::string2utype<double>(tokens.at(k + 1));
      }
      if (tokens[k] == "EBREAK" && k + 1 < tokens.size()) {
        EBREAK = aurostd::string2utype<double>(tokens.at(k + 1));
      }
      if (tokens[k] == "DEPER" && k + 1 < tokens.size()) {
        DEPER = aurostd::string2utype<double>(tokens.at(k + 1));
      }
      if (tokens[k] == "TIME" && k + 1 < tokens.size()) {
        TIME = aurostd::string2utype<double>(tokens.at(k + 1));
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " WEIMIN=" << WEIMIN << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " EBREAK=" << EBREAK << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " DEPER=" << DEPER << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " TIME=" << TIME;
    cerr << endl;
  }

  // ----------------------------------------------------------------------
  // LOAD NWEIGHTS VKPOINT
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD NWEIGHTS KPOINTLIST DATA (" << time_delay(seconds) << ")" << endl;
  }
  uint nkpoints_line = 0;
  nweights = 0;
  nkpoints_irreducible = 0;
  vkpoint_reciprocal.clear();                    // QM KPOINTLIST calculation
  vkpoint_cartesian.clear();                     // QM KPOINTLIST calculation
  vweights.clear();                              // QM KPOINTLIST calculation
  uint minweight = 1e6;
  vline.clear();
  for (size_t iline = 0; iline < vcontent.size(); iline++) {
    if (aurostd::substring2bool(vcontent[iline], "Found") && aurostd::substring2bool(vcontent[iline], "irreducible k-points")) {
      vline.push_back(vcontent[iline]);
      nkpoints_line = iline + 4;
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vline.size()=" << vline.size() << endl;
  }
  for (size_t j = 0; j < vline.size(); j++) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vline.at(" << j << ")=" << vline[j] << endl;
    }
    aurostd::StringSubstInPlace(vline[j], "=", " ");
    aurostd::StringSubstInPlace(vline[j], ";", " ");
    aurostd::string2tokens(vline[j], tokens, " ");
    for (size_t k = 0; k < tokens.size(); k++) {
      if (tokens[k] == "Found" && k + 1 < tokens.size()) {
        nkpoints_irreducible = aurostd::string2utype<uint>(tokens.at(k + 1));
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " nkpoints_irreducible=" << nkpoints_irreducible << endl;
  }

  if (nkpoints_irreducible) {
    for (size_t iline = nkpoints_line; (iline < nkpoints_line + nkpoints_irreducible && iline < vcontent.size()); iline++) {
      aurostd::string2tokens(vcontent[iline], tokens, " ");
      if (tokens.size() == 4) {
        xvector<double> kpoint(3);
        kpoint[1] = aurostd::string2utype<double>(tokens.at(0));
        kpoint[2] = aurostd::string2utype<double>(tokens.at(1));
        kpoint[3] = aurostd::string2utype<double>(tokens.at(2));
        vkpoint_reciprocal.push_back(kpoint);
        vweights.push_back(aurostd::string2utype<uint>(tokens.at(3)));
        minweight = aurostd::min(minweight, aurostd::string2utype<uint>(tokens.at(3)));
        nweights += aurostd::string2utype<uint>(tokens.at(3));
      }
    }
    nkpoints_line += 3;
    for (size_t iline = nkpoints_line; (iline < nkpoints_line + nkpoints_irreducible && iline < vcontent.size()); iline++) {
      aurostd::string2tokens(vcontent[iline], tokens, " ");
      if (tokens.size() == 4) {
        xvector<double> kpoint(3);
        kpoint[1] = aurostd::string2utype<double>(tokens.at(0));
        kpoint[2] = aurostd::string2utype<double>(tokens.at(1));
        kpoint[3] = aurostd::string2utype<double>(tokens.at(2));
        vkpoint_cartesian.push_back(kpoint); // cerr.precision(20);
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vkpoint_reciprocal.size()=" << vkpoint_reciprocal.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vkpoint_cartesian.size()=" << vkpoint_cartesian.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vweights.size()=" << vweights.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " nweights=" << nweights << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " nkpoints_irreducible=" << nkpoints_irreducible << endl;
  }

  // ----------------------------------------------------------------------
  // LOAD CALCULATION STUFF
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD CALCULATION STUFF DATA (" << time_delay(seconds) << ")" << endl;
  }
  // CALCULATION_CORES
  calculation_cores = 1;
  for (size_t iline = 0; iline < vcontent.size(); iline++) {
    if (aurostd::substring2bool(vcontent[iline], "running")) {
      if (aurostd::substring2bool(vcontent[iline], "nodes")) {
        line = vcontent[iline];
        break;
      }
    }
  }
  if (!line.empty()) {
    aurostd::string2tokens(line, tokens);    // cerr << tokens.at(2) << endl;
    if (tokens.size() > 2) {
      calculation_cores = aurostd::string2utype<uint>(tokens.at(2));
    }
  }
  if (calculation_cores < 1) {
    calculation_cores = 1;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " calculation_cores=" << calculation_cores << endl;
  }
  // CALCULATION_TIME
  calculation_time = 0.0;
  for (int iline = (int) vcontent.size() - 1; iline >= 0; iline--) {  // NEW FROM THE BACK
    if (aurostd::substring2bool(vcontent[iline], "time")) {
      if (aurostd::substring2bool(vcontent[iline], "Total")) {
        line = vcontent[iline];
        break;
      }
    }
  }
  if (line.empty()) {
    message << "In OUTCAR (no calculation_time)" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  }
  aurostd::string2tokens(line, tokens);
  if (tokens.size() > 1) {
    calculation_time = aurostd::string2utype<double>(tokens.at(tokens.size() - 1));
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " calculation_time=" << calculation_time << endl;
  }
  // CALCULATION_MEMORY
  calculation_memory = 0.0;
  for (int iline = (int) vcontent.size() - 1; iline >= 0; iline--) {  // NEW FROM THE BACK
    if (aurostd::substring2bool(vcontent[iline], "MB")) {
      if (aurostd::substring2bool(vcontent[iline], "storing wavefunctions")) {
        line = vcontent[iline];
        break;
      }
    }
  }
  if (line.empty()) {
    message << "In OUTCAR (no calculation_memory)" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  }
  aurostd::string2tokens(line, tokens); //   cerr << tokens.at(3) << endl;
  if (tokens.size() > 3) {
    calculation_memory = aurostd::string2utype<double>(tokens.at(3));
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " calculation_memory=" << calculation_memory << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  }

  // ----------------------------------------------------------------------
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " END (" << time_delay(seconds) << ")" << endl;
  }

  // DONE NOW RETURN
  if (ERROR_flag) {
    message << "ERROR_flag set in xOUTCAR";
    if (force_exit) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    } else {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_, QUIET);
    }
    return false;
  }
  return true;
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
bool xOUTCAR::GetXStructure() {
  // CO20211107 - this is good for OUTCAR.static or .bands which only has one structure inside
  // use GetIonicStepsData() otherwise
  const bool LDEBUG = (false || XHOST.DEBUG);
  xstr.clear(); // DX20191220 - uppercase to lowercase clear

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " Trying to build the xstructure from the OUTCAR" << endl;
  }

  // get lattice
  bool found_lattice;
  bool found_types;
  bool found_positions;
  found_lattice = found_types = found_positions = false;
  string token;
  vector<string> tokens;
  xmatrix<double> lattice(3, 3);
  xmatrix<double> klattice(3, 3);
  deque<int> num_each_type;
  double num = 0.0;
  double natoms = 0.0;
  deque<_atom> atoms;
  xvector<double> cpos(3);
  xvector<double> fpos(3);
  const uint vcontent_size = vcontent.size();
  for (uint iline = vcontent_size - 1; iline < vcontent_size; iline--) {  // NEW FROM THE BACK
    // get lattice
    if (!found_lattice && (iline < vcontent_size - 3) && aurostd::substring2bool(vcontent[iline], "direct")) {
      aurostd::string2tokens(vcontent[iline], tokens);
      if ((tokens.size() > 5) && (tokens[0] == "direct") && (tokens[1] == "lattice") && (tokens[2] == "vectors") && (tokens[3] == "reciprocal") && (tokens[4] == "lattice") && (tokens[5] == "vectors")) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " lattice found!" << endl;
        }
        found_lattice = true;
        //
        tokens = GetCorrectEntriesFromLine(vcontent[iline + 1], 6);
        if (tokens.empty()) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " line with lattice vector is ill-written, see: " << vcontent[iline + 1] << endl;
          }
          return false;
        }
        lattice(1, 1) = aurostd::string2utype<double>(tokens[0]);
        lattice(1, 2) = aurostd::string2utype<double>(tokens[1]);
        lattice(1, 3) = aurostd::string2utype<double>(tokens[2]);
        klattice(1, 1) = aurostd::string2utype<double>(tokens[3]);
        klattice(1, 2) = aurostd::string2utype<double>(tokens[4]);
        klattice(1, 3) = aurostd::string2utype<double>(tokens[5]);
        //
        tokens = GetCorrectEntriesFromLine(vcontent[iline + 2], 6);
        if (tokens.empty()) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " line with lattice vector is ill-written, see: " << vcontent[iline + 2] << endl;
          }
          return false;
        }
        lattice(2, 1) = aurostd::string2utype<double>(tokens[0]);
        lattice(2, 2) = aurostd::string2utype<double>(tokens[1]);
        lattice(2, 3) = aurostd::string2utype<double>(tokens[2]);
        klattice(2, 1) = aurostd::string2utype<double>(tokens[3]);
        klattice(2, 2) = aurostd::string2utype<double>(tokens[4]);
        klattice(2, 3) = aurostd::string2utype<double>(tokens[5]);
        //
        tokens = GetCorrectEntriesFromLine(vcontent[iline + 3], 6);
        if (tokens.empty()) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " line with lattice vector is ill-written, see: " << vcontent[iline + 3] << endl;
          }
          return false;
        }
        lattice(3, 1) = aurostd::string2utype<double>(tokens[0]);
        lattice(3, 2) = aurostd::string2utype<double>(tokens[1]);
        lattice(3, 3) = aurostd::string2utype<double>(tokens[2]);
        klattice(3, 1) = aurostd::string2utype<double>(tokens[3]);
        klattice(3, 2) = aurostd::string2utype<double>(tokens[4]);
        klattice(3, 3) = aurostd::string2utype<double>(tokens[5]);
      }
    }
    // get types
    if (!found_types && aurostd::substring2bool(vcontent[iline], "type")) {
      aurostd::string2tokens(vcontent[iline], tokens);
      if ((tokens.size() > 2) && (tokens[0] == "ions") && (tokens[1] == "per") && (tokens[2] == "type")) {
        aurostd::string2tokens(vcontent[iline], tokens, "=");
        if (tokens.size() != 2) {
          continue;
        }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " types found!" << endl;
        }
        found_types = true;
        token = tokens[1];
        aurostd::string2tokens(token, tokens);
        natoms = 0;
        for (size_t i = 0; i < tokens.size(); i++) {
          num = aurostd::string2utype<double>(tokens[i]);
          natoms += num;
          num_each_type.push_back((int) num);
        }
      }
    }
  }
  if (!found_lattice) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " lattice not found" << endl;
    }
    return false;
  }
  if (!found_types) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " types not found" << endl;
    }
    return false;
  }

  const xmatrix<double> f2c = trasp(lattice);
  const xmatrix<double> c2f = inverse(f2c);

  // need types before positions
  for (int iline = (int) vcontent_size - 1; iline >= 0; iline--) {  // NEW FROM THE BACK
    // get positions
    // from bottom up, cartesian coordinates found first
    // this is fortunate as they carry more precision (fractional needs multiplication)
    if (!found_positions && aurostd::substring2bool(vcontent[iline], "coordinates")) {
      aurostd::string2tokens(vcontent[iline], tokens);
      if (!found_positions && (iline < (int) vcontent_size - natoms) && (tokens.size() > 5) && (tokens[0] == "position") && (tokens[1] == "of") && (tokens[2] == "ions") && (tokens[3] == "in") &&
          (tokens[4] == "cartesian") && (tokens[5] == "coordinates")) {
        found_positions = true;
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " positions (cartesian) found!" << endl;
        }
        atoms.clear();
        for (size_t itype = 0; itype < num_each_type.size(); itype++) {
          for (uint iatom = 0; iatom < (uint) num_each_type[itype]; iatom++) {
            atoms.emplace_back();
            atoms.back().type = itype;
            tokens = GetCorrectEntriesFromLine(vcontent[iline + atoms.size()], 3);
            if (tokens.empty()) {
              if (LDEBUG) {
                cerr << __AFLOW_FUNC__ << " line with atom positions is ill-written, see: " << vcontent[iline + atoms.size()] << endl;
              }
              return false;
            }
            cpos(1) = aurostd::string2utype<double>(tokens[0]);
            cpos(2) = aurostd::string2utype<double>(tokens[1]);
            cpos(3) = aurostd::string2utype<double>(tokens[2]);
            atoms.back().cpos = cpos;
            atoms.back().fpos = c2f * cpos;
          }
        }
      }
      if (!found_positions && (iline < (int) vcontent_size - natoms) && (tokens.size() > 5) && (tokens[0] == "position") && (tokens[1] == "of") && (tokens[2] == "ions") && (tokens[3] == "in") &&
          (tokens[4] == "fractional") && (tokens[5] == "coordinates")) {
        found_positions = true;
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " positions (fractional) found!" << endl;
        }
        atoms.clear();
        for (size_t itype = 0; itype < num_each_type.size(); itype++) {
          for (uint iatom = 0; iatom < (uint) num_each_type[itype]; iatom++) {
            atoms.emplace_back();
            atoms.back().type = itype;
            tokens = GetCorrectEntriesFromLine(vcontent[iline + atoms.size()], 3);
            if (tokens.empty()) {
              if (LDEBUG) {
                cerr << __AFLOW_FUNC__ << " line with atom positions is ill-written, see: " << vcontent[iline + atoms.size()] << endl;
              }
              return false;
            }
            fpos(1) = aurostd::string2utype<double>(tokens[0]);
            fpos(2) = aurostd::string2utype<double>(tokens[1]);
            fpos(3) = aurostd::string2utype<double>(tokens[2]);
            atoms.back().fpos = fpos;
            atoms.back().cpos = f2c * fpos;
          }
        }
      }
    }
  }

  if (!found_positions) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " atom positions not found" << endl;
    }
    return false;
  }

  // default!
  for (size_t i = 0; i < atoms.size(); i++) {
    atoms[i].name = "";
    atoms[i].name_is_given = false;
    atoms[i].CleanName();
    atoms[i].CleanSpin();
    atoms[i].ijk(1) = 0;
    atoms[i].ijk(2) = 0;
    atoms[i].ijk(3) = 0; // inside the zero cell...
    atoms[i].corigin(1) = 0.0;
    atoms[i].corigin(2) = 0.0;
    atoms[i].corigin(3) = 0.0; // inside the zero cell
    atoms[i].coord(1) = 0.0;
    atoms[i].coord(2) = 0.0;
    atoms[i].coord(3) = 0.0; // inside the zero cell
    atoms[i].spin = 0.0;
    atoms[i].order_parameter_value = 0;
    atoms[i].order_parameter_atom = false;
    atoms[i].partial_occupation_value = 1.0;
    atoms[i].partial_occupation_flag = false;
  }

  // occupy xstructure
  xstr.title = SYSTEM;
  xstr.lattice = lattice;
  xstr.klattice = 2.0 * pi * klattice;
  // xstr.num_each_type=num_each_type; //done in AddAtom()
  // xstr.comp_each_type.clear();      //done in AddAtom()
  // for(size_t i=0;i<xstr.num_each_type.size();i++){xstr.comp_each_type.push_back((double)xstr.num_each_type[i]);}  //default!, done in AddAtom()
  // xstr.natoms=natoms; //not an attribute
  xstr.f2c = f2c;
  xstr.c2f = c2f;
  xstr.FixLattices();
  // xstr.atoms=atoms;
  for (size_t i = 0; i < atoms.size(); i++) {
    xstr.AddAtom(atoms[i], false); // CO20230319 - add by type
    xstr.partial_occupation_sublattice.push_back(_pocc_no_sublattice_); // default!
  }
  xstr.MakeBasis();

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " lattice" << endl;
    cerr << lattice << endl;
    cerr << __AFLOW_FUNC__ << " klattice" << endl;
    cerr << 2.0 * pi * klattice << endl;
    cerr << __AFLOW_FUNC__ << " reciprocal of lattice (check)" << endl;
    cerr << ReciprocalLattice(lattice) << endl;
    cerr << __AFLOW_FUNC__ << " natoms = " << natoms << endl;
    for (size_t i = 0; i < num_each_type.size(); i++) {
      cerr << __AFLOW_FUNC__ << " num_each_type[" << i << "] = " << num_each_type[i] << endl;
    }
    for (size_t i = 0; i < atoms.size(); i++) {
      cerr << __AFLOW_FUNC__ << " atom[" << i << "] type=" << atoms[i].type << " cpos=" << atoms[i].cpos << " fpos=" << atoms[i].fpos << endl;
    }
    cerr << __AFLOW_FUNC__ << " full xstructure" << endl;
    cerr << xstr << endl;
  }

  return true;
}

int xOUTCAR::isKPointLine(uint iline) {
  xvector<double> kpoint;
  return isKPointLine(iline, kpoint);
}

int xOUTCAR::isKPointLine(uint iline, xvector<double>& _kpoint) {
  if (iline > vcontent.size() - 2) {
    return 0;
  }  // can't be last line, we check iline+1
  if (!aurostd::substring2bool(vcontent[iline], "k-point")) {
    return 0;
  }
  vector<string> tokens;
  aurostd::string2tokens(vcontent[iline], tokens);
  if (tokens.size() != 6) {
    return 0;
  }
  if (tokens[0] != "k-point") {
    return 0;
  }
  // snag kpoint index
  const int kpt_num = (tokens[1] == "***" ? -1 : aurostd::string2utype<int>(tokens[1]));  // issue with printing kpt_num > 999
  // snag kpoint
  xvector<double> kpoint(3);
  for (size_t i = 3; i < tokens.size(); i++) {
    kpoint[i - 2] = aurostd::string2utype<double>(tokens[i]);
  }
  _kpoint = kpoint;
  // check next line
  aurostd::string2tokens(vcontent[iline + 1], tokens);
  if (tokens.size() != 5) {
    return 0;
  }
  if (tokens[0] != "band") {
    return 0;
  }
  if (tokens[1] != "No.") {
    return 0;
  }
  if (tokens[2] != "band") {
    return 0;
  }
  if (tokens[3] != "energies") {
    return 0;
  }
  if (tokens[4] != "occupation") {
    return 0;
  }
  return kpt_num;
}

bool xOUTCAR::GetStartingKPointLines(vector<uint>& ilines) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  ilines.clear();
  if (!(ISPIN == 1 || ISPIN == 2)) {
    if (!GetProperties(content) || !(ISPIN == 1 || ISPIN == 2)) {
      if (LDEBUG) {
        cerr << XPID << "xOUTCAR::GetStartingKPointLines: GetProperties failed." << endl;
      }
      return false;
    }
  }
  const vector<string> tokens;
  xvector<double> kpoint;
  for (size_t iline = 0; iline < vcontent.size(); iline++) {
    if (aurostd::substring2bool(vcontent[iline], "k-point")) {
      if (1 == isKPointLine(iline, kpoint)) {
        ilines.push_back(iline);
      }
    }
  }
  if (ilines.size() != (uint) ISPIN) {
    if (LDEBUG) {
      cerr << XPID << "xOUTCAR::GetStartingKPointLines: ISPIN does not match starting k-point set counts" << endl;
    }
    return false;
  }
  return true;
}

bool xOUTCAR::GetNextKPointLine(uint& iline) {
  iline++;  // march forward
  xvector<double> kpoint;
  while (iline < vcontent.size() && !isKPointLine(iline, kpoint)) {
    iline++;
  } // will work for k-point 1000 (***)
  if (iline > vcontent.size() - 1) {
    return false;
  }
  return true;
}

bool xOUTCAR::ProcessKPoint(uint iline, double EFERMI, vector<double>& b_energies, vector<double>& b_occs) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  b_energies.clear();
  b_occs.clear();

  if (!isKPointLine(iline)) {
    if (LDEBUG) {
      cerr << XPID << "xOUTCAR::ProcessKPoint: this is NOT a k-point line" << endl;
    }
    return false;
  }
  iline += 2; // march forward to first band

  vector<string> tokens;
  uint iband_found;
  double b_energy;
  double b_occ;
  for (uint iband = 1; iband < (uint) NBANDS + 1; iband++) {
    if (iline > vcontent.size() - 1) {
      if (LDEBUG) {
        cerr << XPID << "xOUTCAR::ProcessKPoint: reached the end of the file" << endl;
      }
      return false;
    }
    aurostd::string2tokens(vcontent[iline], tokens);
    if (tokens.size() != 3) {
      if (LDEBUG) {
        cerr << XPID << "xOUTCAR::ProcessKPoint: odd count of tokens(3) for iband=" << iband << endl;
      }
      return false;
    }
    iband_found = aurostd::string2utype<int>(tokens[0]);
    b_energy = aurostd::string2utype<double>(tokens[1]);
    b_occ = aurostd::string2utype<double>(tokens[2]);
    if (iband_found != iband) {
      if (LDEBUG) {
        cerr << XPID << "xOUTCAR::ProcessKPoint: missing iband=" << iband << endl;
      }
      return false;
    }
    b_energies.push_back(b_energy - EFERMI);  // adjust for E-fermi (STATIC) here!!!!!!!!!!!!!!! (not the OUTCAR Efermi, but one from STATIC hopefully)
    b_occs.push_back(b_occ);
    iline++;
  }
  return true;
}

// if energies are equal, sort by occs next in DESCENDING order
bool xOUTCAR::bandEnergyOccCompare::operator()(const bandEnergyOcc& a, const bandEnergyOcc b) const {
  return (a.energy < b.energy || (aurostd::identical(a.energy, b.energy, energy_tol) && a.occ > b.occ));  // sort by energy ASCENDING, occs DESCENDING
}

bool xOUTCAR::orderBands(vector<double>& b_energies, vector<double>& b_occs, double energy_tol) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  if (b_energies.size() != b_occs.size()) {
    if (LDEBUG) {
      cerr << XPID << "xOUTCAR::orderBands: size of energies != size of occupations" << endl;
    }
    return false;
  }
  if (b_energies.empty()) {
    if (LDEBUG) {
      cerr << XPID << "xOUTCAR::orderBands: no input energies or occupations found" << endl;
    }
    return false;
  }
  const uint len = b_energies.size();
  vector<bandEnergyOcc> beo;
  for (size_t i = 0; i < b_energies.size(); i++) {
    beo.emplace_back();
    beo.back().energy = b_energies[i];
    beo.back().occ = b_occs[i];
  }
  std::sort(beo.begin(), beo.end(), bandEnergyOccCompare(energy_tol));
  b_energies.clear();
  b_occs.clear();
  for (size_t i = 0; i < beo.size(); i++) {
    b_energies.push_back(beo[i].energy);
    b_occs.push_back(beo[i].occ);
  }
  if (b_energies.size() != len || b_occs.size() != len) {
    if (LDEBUG) {
      cerr << XPID << "xOUTCAR::orderBands: size mismatch from input" << endl;
    }
    return false;
  }
  return true;
}

bool xOUTCAR::GetBandEdge(vector<double>& b_energies, vector<double>& b_occs, double EFERMI, uint& iedge, double efermi_tol, double energy_tol, double occ_tol) {
  // MOST IMPORTANT FUNCTION FOR FINDING BANDGAP
  // edge is defined here as valence band maximum
  // since energies run from min to max, conduction band min is iedge+1
  // we go backwards, max to E-fermi (above is really below), to find the first band below E-fermi
  // remember, E-fermi is calculated such that an integral of DOS up to E-fermi == # of electrons
  // but it's only reliable in STATIC and not BANDS, where STATIC resolves energies self-consistently (ICHARG=11)
  // BANDS calc only samples path at specific places in BZ, which may not be representative of the full charge density
  // this discrepancy is important, and quantifies the error in the energies of the BANDS calc
  // E-fermi BANDS may differ from E-fermi STATIC by quite a bit (I've seen as much as 0.5 eV)
  // we already shift energies by the INPUT E-fermi, which should be the one from the STATIC calculation
  // first, we want to see if a sharp edge exists
  // this is where a clear fully occupied (occ=2 for ISPIN==1, occ=1 for ISPIN==2) / unoccupied (occ=0) transition exists
  // return immediately if found
  // otherwise, look for a soft edge
  // soft edges occur in band(s) near E-fermi
  // any occupancy (>0) above E-fermi is GARBAGE (an artifact of smearing/sigma for convergence), so ignore
  // but how do we differentiate garbage from non-garbage?
  // algorithmically, we search for the biggest change in occupancy, and define this to be where the transition occurs
  // WARNING: if NUPDOWN is set manually in INCAR, E-fermi must be calculated manually as well (vasp E-fermi only applies to spin-down)
  // http://cms.mpi.univie.ac.at/vasp-forum/viewtopic.php?f=4&t=7442
  // INPUT: energies + occs, E-fermi from STATIC, efermi_tol (largest shift in energy levels between BANDS and STATIC),
  // energy_tol is resolution of energy level values, and occ_tol is resolution of occupation values
  // OUTPUT: iedge - index of VBT

  const bool LDEBUG = (false || XHOST.DEBUG);

  //////////////////////////////////////////////////////////////////////////////
  // tests of stupidity ROUND 1 - START

  if (b_energies.size() != b_occs.size()) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " size of energies != size of occupations" << endl;
    }
    return false;
  }
  if (b_energies.empty()) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " no input energies or occupations found" << endl;
    }
    return false;
  }

  if (!(ISPIN == 1 || ISPIN == 2) || NKPTS == 0 || NBANDS == 0) {
    if (!GetProperties(content) || !(ISPIN == 1 || ISPIN == 2) || NKPTS == 0 || NBANDS == 0) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " GetProperties failed." << endl;
      }
      return false;
    }
  }

  // error in E-fermi, unless otherwise provided, is the difference between the input E-fermi (STATIC) and the one found in this OUTCAR (BANDS)
  // can be significant!
  // think of this as largest possible shift of energies in BANDS relative to real energy levels
  if (efermi_tol == AUROSTD_NAN) {
    efermi_tol = Efermi - EFERMI;                             // generlly, E-fermi BANDS > E-fermi STATIC
    if (std::signbit(efermi_tol)) {
      efermi_tol = energy_tol;
    }  // just in case
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " tol(E-fermi)=" << efermi_tol << endl;
  }

  const double full_occ = (ISPIN == 1 ? 2.0 : 1.0);

  // sort band energies + occs
  // band index is completely phony, low energy bands can appear above higher energy bands:  see FCC/Cu5Lu1_ICSD_103045 k-point 172 (spin=1)
  // this does NOT seem to be an IO error, just dummy VASP indices
  // two bands with the same energy should also be sorted by occs:  see ORC/Al1Pd1Sm1_ICSD_609058 k-point 21 (spin=1)
  // looks like VASP band indices are completely meaningless!
  if (!orderBands(b_energies, b_occs, energy_tol)) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " cannot sort band energies and occupations, this is VERY odd!" << endl;
    }
    return false;
  }

  // now that everything is sorted...
  // occupations should monotonically INCREASE as energy decreases
  bool found_full_occ = false;
  for (int iband = (int) b_energies.size() - 2; iband >= 0; iband--) { //-2 because we compare iband and iband+1
    if (found_full_occ) {
      if (!(aurostd::abs(full_occ - b_occs[iband]) < occ_tol)) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " drop from full_occ to non-full_occ (iband=" << iband + 1 << "), this is VERY odd!" << endl;
        }
        return false;
      }
    }
    if (aurostd::abs(full_occ - b_occs[iband]) < occ_tol) {
      found_full_occ = true;
    }
    if (std::signbit(b_occs[iband] - b_occs[iband + 1]) && aurostd::isdifferent(b_occs[iband], b_occs[iband + 1], occ_tol)) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " occupation decreased at lower energy (iband=" << iband + 1 << "), this is VERY odd!" << endl;
      }
      return false;
    }
  }

  // tests of stupidity ROUND 1 - STOP
  //////////////////////////////////////////////////////////////////////////////

  // sharp edge - A REAL TRANSITION
  // a sharp edge is a sharp edge is a sharp edge
  // very strictly find full_occ vs. no_occ transition
  // don't worry about whether the edge is above or below E-fermi HERE since E-fermi(STATIC) != E-fermi(BANDS)
  // VASP wouldn't occupy above what it believes to be E-fermi, by definition
  for (int iband = (int) b_energies.size() - 2; iband >= 0; iband--) { //-2 because we compare iband and iband+1
    // if(b_energies[iband]>0.020){continue;}  //adjusted for E-fermi (STATIC) already!!!!!
    if (aurostd::abs(full_occ - b_occs[iband]) < occ_tol && aurostd::abs(b_occs[iband + 1]) < occ_tol) {
      iedge = iband;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " sharp edge found between energies " << b_energies[iedge] << ",";
        cerr << b_energies[iedge + 1] << " (ibands=" << iedge + 1 << "," << iedge + 2 << ")" << endl;
      }
      return true;
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // tests of stupidity ROUND 2 - START
  // we do it in this order because we want sharp edges to be found no matter what
  // see for example:  HEX/H2_ICSD_68271

  // band energy edges should not be near E-fermi
  if (b_energies[0] > -efermi_tol) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " first band energy is near/above E-fermi!" << endl;
    }
    return false;
  }
  if (b_energies[b_energies.size() - 1] < efermi_tol) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " last band energy is near/below E-fermi!" << endl;
    }
    return false;
  }

  // first band should be fully occupied
  if (aurostd::abs(full_occ - b_occs[0]) >= occ_tol) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " first band is not fully occupied!" << endl;
    }
    return false;
  }
  // last band should be fully un-occupied
  if (aurostd::abs(b_occs[b_occs.size() - 1]) >= occ_tol) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " last band is not fully un-occupied!" << endl;
    }
    // return false; //don't exit here, we might still have a legitimate soft edge, see FCC/C6Mn20W3_ICSD_618279
  }

  // tests of stupidity ROUND 2 - STOP
  //////////////////////////////////////////////////////////////////////////////

  // soft edge - a fuzzy transition due to numerical/convergence issues
  // here we care whether we are near E-fermi
  // two issues that must be overcome here - (full_occ vs. non-full_occ) and (E-fermi STATIC vs. E-fermi BANDS)
  // 1. full_occ vs. non-full_occ: find the largest delta in occupation, set that to be the transition
  // 2. E-fermi STATIC vs. E-fermi BANDS: look no higher than the difference of the two (tolerance of E-fermi)
  // so don't look above E-fermi+tol (E-fermi BANDS)
  // VASP wouldn't occupy above what it believes to be E-fermi, by definition
  double occ_delta;
  double occ_delta_max = 0;
  found_full_occ = false;
  for (int iband = (int) b_energies.size() - 2; iband >= 0; iband--) { //-2 because we compare iband and iband+1
    if (b_energies[iband] > efermi_tol) {
      continue;
    }  // adjusted for E-fermi (STATIC) already!!!!!
    if (aurostd::abs(full_occ - b_occs[iband]) < occ_tol) {  // if we find full_occ twice in a row, stop
      if (found_full_occ) {
        break;
      }
      found_full_occ = true;
    }
    occ_delta = b_occs[iband] - b_occs[iband + 1];
    if (occ_delta > occ_delta_max) {
      iedge = iband;
      occ_delta_max = occ_delta;
    }
  }

  if (aurostd::abs(occ_delta_max) < occ_tol) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " no edge found!" << endl;
    }
    return false;
  }

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " soft edge found between energies " << b_energies[iedge] << ",";
    cerr << b_energies[iedge + 1] << " (ibands=" << iedge + 1 << "," << iedge + 2 << ")" << endl;
  }
  return true;
}

bool xOUTCAR::identicalKPoints(vector<xvector<double>>& vkpoints, uint kpt1, uint kpt2, double tol) {
  if (kpt1 == kpt2) {
    return true;
  }
  return identicalKPoints(vkpoints[kpt1], vkpoints[kpt2], tol);
}

bool xOUTCAR::identicalKPoints(xvector<double>& kpoint1, xvector<double>& kpoint2, double tol) {
  if (aurostd::identical(kpoint1, kpoint2, tol)) {
    return true;
  }  // we can be dam sure at this tolerance
  return false;
}

bool xOUTCAR::removeDuplicateKPoints(vector<xvector<double>>& vkpoints, vector<uint>& vikpt) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  vector<uint> vikpt_unique;
  bool found;
  for (size_t i = 0; i < vikpt.size(); i++) {
    if (vikpt[i] > vkpoints.size() - 1) {
      return false;
    }
    found = false;
    for (size_t j = 0; j < vikpt_unique.size() && !found; j++) {
      if (identicalKPoints(vkpoints, vikpt[i], vikpt_unique[j])) {
        if (LDEBUG) {
          cerr << XPID << "xOUTCAR::removeDuplicateKPoints: removing duplicate k-point ";
          cerr << vkpoints[vikpt[i]] << " (kpt=" << vikpt[i] << ")" << endl;
        }
        found = true;
      }
    }
    if (!found) {
      vikpt_unique.push_back(vikpt[i]);
    }
  }
  vikpt = vikpt_unique;
  return true;
}

bool xOUTCAR::removeDuplicateKPoints(vector<vector<xvector<double>>>& vkpoints, vector<uint>& vikpt, vector<uint>& vispin) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  if (vikpt.size() != vispin.size()) {
    return false;
  }  // test of stupidity
  vector<uint> vikpt_unique;
  vector<uint> vispin_unique;
  bool found;
  for (size_t i = 0; i < vikpt.size(); i++) {
    if (vikpt[i] > vkpoints[vispin[i]].size() - 1) {
      return false;
    }
    found = false;
    for (size_t j = 0; j < vikpt_unique.size() && !found; j++) {
      if (identicalKPoints(vkpoints[vispin[i]][vikpt[i]], vkpoints[vispin_unique[j]][vikpt_unique[j]])) {
        if (LDEBUG) {
          cerr << XPID << "xOUTCAR::removeDuplicateKPoints: removing duplicate k-point ";
          cerr << vkpoints[vispin[i]][vikpt[i]] << " (kpt=" << vikpt[i] << ",spin=" << vispin[i] << ")" << endl;
        }
        found = true;
      }
    }
    if (!found) {
      vikpt_unique.push_back(vikpt[i]);
      vispin_unique.push_back(vispin[i]);
    }
  }
  vikpt = vikpt_unique;
  vispin = vispin_unique;
  return true;
}

double xOUTCAR::minimumDistanceKPoints(vector<xvector<double>>& vkpoints, uint ikp1, uint ikp2) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  double dist_min = AUROSTD_MAX_DOUBLE;
  if ((ikp1 > vkpoints.size() - 1) || (ikp2 > vkpoints.size() - 1)) {
    return dist_min;
  }
  if (xstr.atoms.empty()) { // assume it's already loaded
    if (!GetXStructure()) {
      if (LDEBUG) {
        cerr << XPID << "xOUTCAR::minimumDistanceKPoints: GetXStructure failed." << endl;
      }
      return dist_min;
    }
  }
  if (LDEBUG) {
    cerr << XPID << "xOUTCAR::minimumDistanceKPoints: subtracting ";
    cerr << "kpoint[" << ikp1 << "]";
    cerr << " from ";
    cerr << "kpoint[" << ikp2 << "]";
    cerr << endl;
  }
  // special case, identicalKPoints()
  if (ikp1 == ikp2) {
    dist_min = 0.0;
    if (LDEBUG) {
      cerr << XPID << "xOUTCAR::minimumDistanceKPoints: Found identical k-points, distance in cartesian = " << dist_min << endl;
    }
    return dist_min;
  }
  return minimumDistanceKPoints(vkpoints[ikp1], vkpoints[ikp2]);
}

double xOUTCAR::minimumDistanceKPoints(xvector<double>& kpoint1_kl, xvector<double>& kpoint2_kl) {  // these are in units of reciprocal lattice vectors
  const bool LDEBUG = (false || XHOST.DEBUG);
  double dist_min = AUROSTD_MAX_DOUBLE;
  if (LDEBUG) {
    cerr << XPID << "xOUTCAR::minimumDistanceKPoints: kpoint1_kl=" << kpoint1_kl << endl;
    cerr << XPID << "xOUTCAR::minimumDistanceKPoints: kpoint2_kl=" << kpoint2_kl << endl;
  }
  // special case, identicalKPoints()
  if (identicalKPoints(kpoint1_kl, kpoint2_kl)) {
    dist_min = 0.0;
    if (LDEBUG) {
      cerr << XPID << "xOUTCAR::minimumDistanceKPoints: Found identical k-points, distance in cartesian = " << dist_min << endl;
    }
    return dist_min;
  }
  const xmatrix<double> klattice = ReciprocalLattice(xstr);  // repetita iuvant
  const xmatrix<double> metric_tensor = MetricTensor(xstr);
  if (LDEBUG) {
    cerr << XPID << "xOUTCAR::minimumDistanceKPoints: metric tensor" << endl;
    cerr << metric_tensor << endl;
    cerr << XPID << "xOUTCAR::minimumDistanceKPoints: klattice" << endl;
    cerr << klattice << endl;
  }

  const xmatrix<double> kf2c = trasp(klattice);
  const xvector<double> kpoint1 = kf2c * kpoint1_kl;  // units of 2pi/scale
  const xvector<double> kpoint2 = kf2c * kpoint2_kl;  // units of 2pi/scale
  if (LDEBUG) {
    cerr << XPID << "xOUTCAR::minimumDistanceKPoints: subtracting (units of 2*pi/scale) " << endl;
    cerr << "                                 kpoint=" << kpoint1 << endl;
    cerr << "                                 kpoint=" << kpoint2 << endl;
  }

  // make sure to account for skewed cells
  xvector<int> dims = LatticeDimensionSphere(klattice, RadiusSphereLattice(klattice));

  xvector<double> vdist_kcart;
  xvector<double> vdist_kcart_min;
  double dist = 0.0;
  for (int i = -dims[1]; i <= dims[1]; i++) {
    for (int j = -dims[2]; j <= dims[2]; j++) {
      for (int k = -dims[3]; k <= dims[3]; k++) {
        vdist_kcart = (kpoint1 - kpoint2 + double(i) * klattice(1) + double(j) * klattice(2) + double(k) * klattice(3));
        dist = aurostd::modulus(vdist_kcart);
        if (dist < dist_min) {
          vdist_kcart_min = vdist_kcart;
          dist_min = dist;
        }
      }
    }
  }
  if (LDEBUG) {
    cerr << XPID << "xOUTCAR::minimumDistanceKPoints: distance vector in reciprocal space  = " << vdist_kcart_min << endl;
  }

  // convert distance to cartesian in real space
  xvector<double> vdist_dcart;
  double dist_dcart;
  vdist_dcart(1) = sum(metric_tensor(1) * vdist_kcart_min(1));  // remember metric_tensor is symmetric!
  vdist_dcart(2) = sum(metric_tensor(2) * vdist_kcart_min(2));  // remember metric_tensor is symmetric!
  vdist_dcart(3) = sum(metric_tensor(3) * vdist_kcart_min(3));  // remember metric_tensor is symmetric!
  vdist_dcart /= (2 * pi);  // conversion back requires undoing 2pi
  if (LDEBUG) {
    cerr << XPID << "xOUTCAR::minimumDistanceKPoints: distance vector in real space        = " << vdist_dcart << endl;
  }
  dist_dcart = aurostd::modulus(vdist_dcart);
  if (LDEBUG) {
    cerr << XPID << "xOUTCAR::minimumDistanceKPoints: distance between kpoints (angstroms) = " << dist_dcart << endl;
  }
  return dist_dcart;
}

/// @brief gets the band gap information
/// @param EFERMI Fermi energy
/// @param efermi_tol Fermi energy tolerance
/// @param energy_tol eigenenergy tolerance
/// @param occ_tol occupation tolerance
/// @return true if successful
/// @authors
/// @mod{CO,20171004,created function}
/// @mod{SD,20250407,added doxy}
bool xOUTCAR::GetBandGap(double EFERMI, double efermi_tol, double energy_tol, double occ_tol) {
  // some nice examples when debugging - validated by CO20171006
  // aflow --bandgap=/common/ICSD/LIB/RHL/Ag1Ni1O2_ICSD_73974
  // System        :   Ag1Ni1O2_ICSD_73974
  // Spin tag      :   2
  // Fermi level   :  +2.9623e+00
  //                   VBT           CBB           Egap          Egap_fit     Type
  // Majority Spin :  -1.0000e+00   -1.0000e+00   -1.0000e+09   -1.0000e+09   metal
  // Minority Spin :  -8.7090e-01   +1.4751e+00   +2.3460e+00   +4.0754e+00   insulator-indirect
  // Net Result    :  -1.0000e+00   -1.0000e+00   -1.0000e+09   -1.0000e+09   half-metal
  //
  // aflow --bandgap=/common/ICSD/LIB/FCC/Ni1O1_ICSD_92133
  // System        :   Ni1O1_ICSD_92133
  // Spin tag      :   2
  // Fermi level   :  +5.7976e+00
  //                   VBT           CBB           Egap          Egap_fit     Type
  // Majority Spin :  +1.0100e-02   +2.1721e+00   +2.1620e+00   +3.8274e+00   insulator-indirect
  // Minority Spin :  -1.1203e+00   +1.8017e+00   +2.9220e+00   +4.8519e+00   insulator-indirect
  // Net Result    :  +1.0100e-02   +1.8017e+00   +1.7916e+00   +3.3281e+00   insulator-indirect_spin-polarized
  //
  // aflow --bandgap=/common/ICSD/LIB/FCC/Cr1O2_ICSD_186838
  // System        :   Cr1O2_ICSD_186838
  // Spin tag      :   2
  // Fermi level   :  +5.1970e+00
  //                   VBT           CBB           Egap          Egap_fit     Type
  // Majority Spin :  -1.0000e+00   -1.0000e+00   -1.0000e+09   -1.0000e+09   metal
  // Minority Spin :  -4.8020e-01   +2.1972e+00   +2.6774e+00   +4.5221e+00   insulator-indirect
  // Net Result    :  -1.0000e+00   -1.0000e+00   -1.0000e+09   -1.0000e+09   half-metal
  //
  // aflow --bandgap=/common/ICSD/LIB/FCC/Ba2Dy1Nb1O6_ICSD_109156
  // System        :   Ba2Dy1Nb1O6_ICSD_109156
  // Spin tag      :   1
  // Fermi level   :  +2.8134e+00
  //                   VBT           CBB           Egap          Egap_fit     Type
  // Net Result    :  -6.1000e-03   +3.0853e+00   +3.0914e+00   +5.0802e+00   insulator-direct
  //
  // aflow --bandgap=/common/ICSD/LIB/FCC/Si1_ICSD_150530
  // System        :   Si1_ICSD_150530
  // Spin tag      :   1
  // Fermi level   :  +5.6171e+00
  //                   VBT           CBB           Egap          Egap_fit     Type
  // Net Result    :  +7.2000e-03   +6.1720e-01   +6.1000e-01   +1.7353e+00   insulator-indirect
  //
  // aflow --bandgap=/common/ICSD/LIB/BCC/Fe1_ICSD_52258
  // System        :   Fe1_ICSD_52258
  // Spin tag      :   2
  // Fermi level   :  +4.8034e+00
  //                   VBT           CBB           Egap          Egap_fit     Type
  // Majority Spin :  -1.0000e+00   -1.0000e+00   -1.0000e+09   -1.0000e+09   metal
  // Minority Spin :  -1.0000e+00   -1.0000e+00   -1.0000e+09   -1.0000e+09   metal
  // Net Result    :  -1.0000e+00   -1.0000e+00   -1.0000e+09   -1.0000e+09   metal
  //
  // aflow --bandgap=/common/ICSD/LIB/HEX/H2_ICSD_68271
  // WARNING - xOUTCAR::GetBandGap: unable to detect band edge for k-point 1 (spin=2), this system should be rerun with a wider energy range
  //
  //    filename=[stringstream]
  //
  //    pflow::BANDGAP: . failed
  //
  //
  // aflow --bandgap=/common/ICSD/LIB/HEX/Br2Mn1_ICSD_60250
  // System        :   Br2Mn1_ICSD_60250
  // Spin tag      :   2
  // Fermi level   :  +1.8650e-01
  //                   VBT           CBB           Egap          Egap_fit     Type
  // Majority Spin :  -1.3600e-02   +3.3703e+00   +3.3839e+00   +5.4745e+00   insulator-indirect
  // Minority Spin :  -5.8990e-01   +3.4065e+00   +3.9964e+00   +6.3001e+00   insulator-indirect
  // Net Result    :  -1.3600e-02   +3.3703e+00   +3.3839e+00   +5.4745e+00   insulator-indirect

  // repetita iuvant!!!!!!!!
  const bool LDEBUG = (false || XHOST.DEBUG);
  stringstream message;
  const bool force_exit = false; //[CO20210621 - we will now exit gracefully]XHOST.POSTPROCESS; //SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only

  if ((content.empty()) || (vcontent.empty())) {
    message << "xOUTCAR needs to be loaded before." << endl;
    message << "GetProperties(const stringstream&);" << endl;
    message << "GetProperties(const string&);" << endl;
    message << "GetPropertiesFile(const string&);" << endl;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_MISSING_);
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " OUTCAR content found" << endl;
  }

  // quick check if GetProperties() failed
  if (!(ISPIN == 1 || ISPIN == 2) || NKPTS == 0 || NBANDS == 0) {
    if (!GetProperties(content) || !(ISPIN == 1 || ISPIN == 2) || NKPTS == 0 || NBANDS == 0) {
      message << "GetProperties() failed";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " OUTCAR properties retrieved" << endl;
  }

  // GET FERMI LEVEL
  if (EFERMI == AUROSTD_NAN) {
    // VERY IMPORTANT - CO20171009 from discussion with SC
    // we strongly prefer to use Efermi from OUTCAR.static, not OUTCAR.bands
    // bands is not self-consistent (ICHARG=11), only used to determine energy states
    // since it is not self-consistent, it does not fill in states correctly
    // static is the best way to go!

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " Using Efermi from CURRENT OUTCAR (recommended to use Efermi from OUTCAR.static)" << endl;
    }
    EFERMI = Efermi;
  }

  // ----------------------------------------------------------------------
  // GET POSCAR INFO FOR GETBANDGAP()
  // silly to load in POSCAR when all the information is in OUTCAR
  if (xstr.atoms.empty()) { // assume it's already loaded
    if (!GetXStructure()) {
      message << "GetXStructure() failed";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " xstructure built from OUTCAR" << endl;
  }

  double kpt_tol;
  if (xstr.CalculateSymmetry()) {
    kpt_tol = xstr.sym_eps;
  } else {
    kpt_tol = SYM::defaultTolerance(xstr);
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " kpt_tol=" << kpt_tol << endl;
  }
  // ----------------------------------------------------------------------

  vector<uint> starting_lines;
  if (!GetStartingKPointLines(starting_lines)) {
    message << "Unable to grab starting k-point lines";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " reading k-points starting at line";
    cerr << (starting_lines.size() == 2 ? "s: " : ": ");
    cerr << starting_lines[0];
    cerr << (starting_lines.size() == 2 ? " and " + aurostd::utype2string(starting_lines[1]) : "");
    cerr << endl;
  }

  uint iline;
  vector<double> b_energies;
  vector<double> b_occs;
  uint iedge;
  xvector<double> kpoint; // reciprocal

  // resize for spin
  vector<vector<xvector<double>>> vkpoints;
  vkpoints.resize(ISPIN); // reciprocal
  vector<vector<uint>> vedges;
  vedges.resize(ISPIN);
  vector<vector<double>> vVBT;
  vVBT.resize(ISPIN);
  vector<vector<double>> vCBB;
  vCBB.resize(ISPIN);
  vector<uint> empty_channel;
  empty_channel.resize(ISPIN);
  vector<uint> partially_empty_channel;
  partially_empty_channel.resize(ISPIN);

  bool band_edge_found;
  int kpt_found;
  int first_kpt_empty;
  int first_spin_empty;
  first_kpt_empty = first_spin_empty = -1;

  // initialize empty_channel to 0's
  for (uint i = 0; i < (uint) ISPIN; i++) {
    empty_channel[i] = 0;
  }

  for (uint ispin = 0; ispin < (uint) ISPIN; ispin++) {
    iline = starting_lines[ispin];
    for (uint ikpt = 1; ikpt < (uint) NKPTS + 1; ikpt++) {
      // is OUTCAR poorly written?
      kpt_found = isKPointLine(iline, kpoint);
      if ((int) ikpt != kpt_found) {
        if (kpt_found == -1 && ikpt < 1000) {  // ONLY exception, above 999, VASP starting writing ***, so ignore this check starting at 1000
          message << "Missing k-point " << ikpt << " (spin=" << ispin + 1 << ")";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
        }
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " looking at k-point[" << ikpt << "]=" << kpoint << " (spin=" << ispin + 1 << ")" << endl;
      }

      /////////////////////////////////////////////
      // BANDS analysis - START
      vkpoints[ispin].push_back(kpoint);  // push back kpoints even if empty channel

      if (!ProcessKPoint(iline, EFERMI, b_energies, b_occs)) {
        message << "Unable to process k-point " << ikpt << " (spin=" << ispin + 1 << +")";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
      }

      // keep vedges, vVBT, vCBB same size as vkpoints
      vedges[ispin].push_back(AUROSTD_MAX_UINT);
      vVBT[ispin].push_back(AUROSTD_MAX_DOUBLE);
      vCBB[ispin].push_back(AUROSTD_MAX_DOUBLE);

      band_edge_found = GetBandEdge(b_energies, b_occs, EFERMI, iedge, efermi_tol, energy_tol, occ_tol);

      if (band_edge_found) {  // normal case
        if (empty_channel[ispin] == 1) {
          // if we get here, empty channel detected for k-point 1, but not for current k-point
          // therefore, 1 is the outlier
          // inconsistency within channel
          partially_empty_channel[ispin] = 1;
          if (first_kpt_empty == -1) {
            first_kpt_empty = 1;
          }
          if (first_spin_empty == -1) {
            first_spin_empty = ispin;
          }
          // ERROR = __AFLOW_FUNC__ + " unable to detetc band edge for k-point "+aurostd::utype2string(1)+
          //   " (spin="+aurostd::utype2string(ispin+1)+") \n";
          // return false;
        }//[CO20210621 - always save good vedges, even if empty channel] else {
        vedges[ispin].back() = iedge;
        vVBT[ispin].back() = b_energies[iedge];
        vCBB[ispin].back() = b_energies[iedge + 1];
        //[CO20210621 - always save good vedges, even if empty channel]}
      } else {  // special case
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " no band edge found for k-point[" << ikpt << "]=" << kpoint << " (spin=" << ispin + 1 << ")" << endl;
        }
        if (ikpt != 1) {
          // if we get here, kpt==1 presents an edge, but current k-point does not
          // therefore, current k-point is the outlier
          // inconsistency within channel
          partially_empty_channel[ispin] = 1;
          // ERROR = __AFLOW_FUNC__ + " unable to detect band edge for k-point "+aurostd::utype2string(ikpt)+
          //   " (spin="+aurostd::utype2string(ispin+1)+") \n";
          // return false;
        }
        empty_channel[ispin] = 1;
        if (first_kpt_empty == -1) {
          first_kpt_empty = ikpt;
        }
        if (first_spin_empty == -1) {
          first_spin_empty = ispin;
        }
      }

      // BANDS analysis - STOP
      /////////////////////////////////////////////

      // can we find the next kpoint?
      if (ikpt < (uint) NKPTS && !GetNextKPointLine(iline)) {
        message << "Missing k-point " << ikpt + 1 << " (spin=" << ispin + 1 << +")";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
      }
    }
  }

  // first, we differentiate between metals / insulators
  // in metals, bands cross at the edge
  // this manifests as a change in the edge index
  // think of this as a band "appearing"/"disappearing" at the edge
  // or a "conservation" of bands below edge, lack of conservation means crossing, or a metal
  // enum BROAD_TYPES {empty,metal,insulator};
  // enum EMPTY_TYPES {empty_all,empty_partial};
  // enum INSULATOR_TYPES {insulator_direct,insulator_indirect};
  // enum GAP_TYPES {zero_gap,non_zero_gap};

  // resize for spin
  vector<BROAD_TYPES> broad_type;
  broad_type.resize(ISPIN);
  vector<EMPTY_TYPES> empty_type;
  empty_type.resize(ISPIN);
  vector<INSULATOR_TYPES> insulator_type;
  insulator_type.resize(ISPIN);
  vector<GAP_TYPES> gap_type;
  gap_type.resize(ISPIN);
  vector<double> gap;
  gap.resize(ISPIN);
  vector<uint> vimax_VBTs;
  vimax_VBTs.resize(ISPIN); // for each spin
  vector<uint> vimin_CBBs;
  vimin_CBBs.resize(ISPIN); // for each spin

  // broad_type
  for (uint ispin = 0; ispin < (uint) ISPIN; ispin++) {
    //[CO20210621 - check for metal first, even if empty]if(empty_channel[ispin]){broad_type[ispin]=empty;continue;}
    broad_type[ispin] = insulator;  // insulator unless proven otherwise
    for (size_t ikpt = 1; ikpt < vedges[ispin].size(); ikpt++) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " looking at vedges[ispin=" << ispin << "][ikpt-1=" << ikpt - 1 << "]=" << vedges[ispin][ikpt - 1] << " vs. " << vedges[ispin][ikpt] << "=vedges[ispin=" << ispin
             << "][ikpt=" << ikpt << "]" << endl;
      }
      if (vedges[ispin][ikpt - 1] == AUROSTD_MAX_UINT || vedges[ispin][ikpt] == AUROSTD_MAX_UINT) {
        continue;
      }
      if (vedges[ispin][ikpt - 1] != vedges[ispin][ikpt]) {
        broad_type[ispin] = metal; // metal
        break;
      }
    }
    if (empty_channel[ispin] && broad_type[ispin] != metal) {
      broad_type[ispin] = empty;
    }
  }

  for (uint ispin = 0; ispin < (uint) ISPIN; ispin++) {
    if (broad_type[ispin] == empty) {
      message << "Unable to detect band edge for k-point " << first_kpt_empty << " (spin=" << first_spin_empty + 1 << ")," << " this system should be rerun with a wider energy range";
      if (force_exit) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      } else {
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
        return false;
      }
    }
  }

  double max_VBT = 0.0;
  double min_CBB = 0.0;
  double dist = 0.0;
  double dist_min = 0.0;
  uint imax_VBT = 0;
  uint imin_CBB = 0; // absolutes, then our final picks
  vector<uint> vimax_VBTs_equiv;
  vector<uint> vimin_CBBs_equiv; // within tolerance of absolutes
  bool vbt_duplicate_remove;
  bool cbb_duplicate_remove;

  // specific types here
  for (uint ispin = 0; ispin < (uint) ISPIN; ispin++) {
    if (broad_type[ispin] == empty) {
      gap[ispin] = 0.0;
      if (partially_empty_channel[ispin] == 1) {
        empty_type[ispin] = empty_partial;
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " empty (partial) channel found spin=" << ispin + 1 << endl;
        }
      } else {
        empty_type[ispin] = empty_all;
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " empty (all) channel found spin=" << ispin + 1 << endl;
        }
      }
      continue;
    }
    // metal, only one type
    // modify starting here in the future for SEMI-metals (need DOS)
    else if (broad_type[ispin] == metal) {
      gap[ispin] = 0.0;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " metal found in spin=" << ispin + 1 << endl;
      }
      continue;
    } else if (broad_type[ispin] != insulator) {  // test of stupidity
      message << "Unknown material type (!empty && !metal && !insulator) (spin=" << ispin + 1 << ")";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " insulator found in spin=" << ispin + 1 << endl;
    }
    max_VBT = (-1.0) * AUROSTD_MAX_DOUBLE;
    min_CBB = AUROSTD_MAX_DOUBLE;
    // determine if direct/indirect insulator
    // first, get absolute max of VBT/min of CBB
    for (size_t ikpt = 0; ikpt < vkpoints[ispin].size(); ikpt++) {
      if (vVBT[ispin][ikpt] > max_VBT) {
        imax_VBT = ikpt;
        max_VBT = vVBT[ispin][imax_VBT];
      }
      if (vCBB[ispin][ikpt] < min_CBB) {
        imin_CBB = ikpt;
        min_CBB = vCBB[ispin][imin_CBB];
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " absolute VBT=" << max_VBT << " at k-point " << vkpoints[ispin][imax_VBT] << " (kpt=" << imax_VBT << ",spin=" << ispin + 1 << ")" << endl;
      cerr << __AFLOW_FUNC__ << " absolute CBB=" << min_CBB << " at k-point " << vkpoints[ispin][imin_CBB] << " (kpt=" << imin_CBB << ",spin=" << ispin + 1 << ")" << endl;
    }
    vimax_VBTs_equiv.clear();
    vimin_CBBs_equiv.clear();
    // grab any equivalently high/low extrema
    for (size_t ikpt = 0; ikpt < vkpoints[ispin].size(); ikpt++) {
      if (aurostd::abs(max_VBT - vVBT[ispin][ikpt]) < energy_tol) {
        vimax_VBTs_equiv.push_back(ikpt);
        if (LDEBUG && ikpt != imax_VBT) {
          cerr << __AFLOW_FUNC__ << " found equivalent VBT at " << vkpoints[ispin][ikpt] << " (kpt=" << ikpt << ",spin=" << ispin + 1 << ")" << endl;
        }
      }
      if (aurostd::abs(min_CBB - vCBB[ispin][ikpt]) < energy_tol) {
        vimin_CBBs_equiv.push_back(ikpt);
        if (LDEBUG && ikpt != imin_CBB) {
          cerr << __AFLOW_FUNC__ << " found equivalent CBB at " << vkpoints[ispin][ikpt] << " (kpt=" << ikpt << ",spin=" << ispin + 1 << ")" << endl;
        }
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " removing duplicate k-points from VBT search (spin=" << ispin + 1 << ")" << endl;
    }
    vbt_duplicate_remove = removeDuplicateKPoints(vkpoints[ispin], vimax_VBTs_equiv);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " removing duplicate k-points from CBB search (spin=" << ispin + 1 << ")" << endl;
    }
    cbb_duplicate_remove = removeDuplicateKPoints(vkpoints[ispin], vimin_CBBs_equiv);
    if (!vbt_duplicate_remove || !cbb_duplicate_remove) {
      message << "Cannot find equivalent band extrema (spin=" << ispin + 1 << ")";
      if (force_exit) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      } else {
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
        return false;
      }
    }
    // if we found more than one possible extrema, minimize kpoint distance between max/min
    // this simulates the electron trying to reduce momentum needed to conduct
    if (vimax_VBTs_equiv.size() == 1 && vimin_CBBs_equiv.size() == 1) {  // easy case, keep already defined imax/imin
      dist_min = minimumDistanceKPoints(vkpoints[ispin], imax_VBT, imin_CBB);
    } else { // find minimum distance pairs of max/min
      dist_min = AUROSTD_MAX_DOUBLE;
      for (size_t imax = 0; imax < vimax_VBTs_equiv.size(); imax++) {
        for (size_t imin = 0; imin < vimin_CBBs_equiv.size(); imin++) {
          dist = minimumDistanceKPoints(vkpoints[ispin], vimax_VBTs_equiv[imax], vimin_CBBs_equiv[imin]);
          if (dist < dist_min) {
            dist_min = dist;
            imax_VBT = vimax_VBTs_equiv[imax];
            imin_CBB = vimin_CBBs_equiv[imin];
          }
        }
      }
    }
    if (dist_min == AUROSTD_MAX_DOUBLE) {
      message << "Cannot calculate k-point distance between band extrema (spin=" << ispin + 1 << ")";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " dist_min=" << dist_min << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vCBB[ispin=" << ispin << "][imin_CBB=" << imin_CBB << "]=" << vCBB[ispin][imin_CBB] << endl;
      cerr << __AFLOW_FUNC__ << " vVBT[ispin=" << ispin << "][imax_VBT=" << imax_VBT << "]=" << vVBT[ispin][imax_VBT] << endl;
    }
    gap[ispin] = vCBB[ispin][imin_CBB] - vVBT[ispin][imax_VBT];
    vimin_CBBs[ispin] = imin_CBB; // save for later spin loops
    vimax_VBTs[ispin] = imax_VBT; // save for later spin loops
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " gap[ispin=" << ispin << "]=" << gap[ispin] << endl;
    }
    if (gap[ispin] < 0) { // test of stupidity, this should NEVER happen
      message << "Negative band gap found, something broke (spin=" << ispin + 1 << ")";
      if (force_exit) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      } else {
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
        return false;
      }
    }
    gap_type[ispin] = (gap[ispin] < energy_tol ? zero_gap : non_zero_gap);
    insulator_type[ispin] = (aurostd::abs(dist_min) < kpt_tol ? insulator_direct : insulator_indirect);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ";
      cerr << (insulator_type[ispin] == insulator_indirect ? "IN" : "") << "DIRECT insulator ";
      cerr << "gap=" << gap[ispin] << " ";
      cerr << (gap_type[ispin] == zero_gap ? "(zero-gap) " : "");
      cerr << "found in spin=" << ispin + 1 << endl;
    }
  }

  // get spin-channel specific properties
  valence_band_max.clear();
  valence_band_max.resize(ISPIN);
  conduction_band_min.clear();
  conduction_band_min.resize(ISPIN);
  Egap.clear();
  Egap.resize(ISPIN);
  Egap_type.clear();
  Egap_type.resize(ISPIN);
  Egap_fit.clear();
  Egap_fit.resize(ISPIN);

  for (uint ispin = 0; ispin < (uint) ISPIN; ispin++) {
    if (broad_type[ispin] == empty) {
      valence_band_max[ispin] = _METALEDGE_;
      conduction_band_min[ispin] = _METALEDGE_;
      Egap[ispin] = _METALGAP_;
      Egap_type[ispin] = "empty" + string(empty_type[ispin] == empty_partial ? "-partially" : "");
    } else if (broad_type[ispin] == metal) {
      valence_band_max[ispin] = _METALEDGE_;
      conduction_band_min[ispin] = _METALEDGE_;
      Egap[ispin] = _METALGAP_;
      Egap_type[ispin] = "metal";
    } else if (broad_type[ispin] != insulator) { // test of stupidity
      message << "Unknown material type (!empty && !metal && !insulator) (spin=" << ispin + 1 << ")";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    } else {
      valence_band_max[ispin] = vVBT[ispin][vimax_VBTs[ispin]];
      conduction_band_min[ispin] = vCBB[ispin][vimin_CBBs[ispin]];
      Egap[ispin] = (gap_type[ispin] == zero_gap ? 0.0 : gap[ispin]);
      Egap_type[ispin] = "insulator-" + string(insulator_type[ispin] == insulator_indirect ? "in" : "") + "direct" + string(gap_type[ispin] == zero_gap ? "_zero-gap" : "");
    }
    Egap_fit[ispin] = (broad_type[ispin] == empty || broad_type[ispin] == metal ? _METALGAP_ : 1.348 * Egap[ispin] + 0.913);
  }

  // get system-wide properties
  if (broad_type[0] == empty || (ISPIN == 2 && broad_type[1] == empty)) {
    if (ISPIN == 1 || (ISPIN == 2 && broad_type[0] == empty && broad_type[1] == empty)) {
      valence_band_max_net = _METALEDGE_;
      conduction_band_min_net = _METALEDGE_;
      Egap_net = _METALGAP_;
      Egap_type_net = "empty" + string(empty_type[0] == empty_partial || (ISPIN == 2 && empty_type[1] == empty_partial) ? "-partially" : "");
      ;
    } else {
      if (broad_type[0] == empty) { // grab properties of spin-channel 2
        const uint ispin = 1;
        valence_band_max_net = (broad_type[ispin] == metal ? _METALEDGE_ : vVBT[ispin][vimax_VBTs[ispin]]);
        conduction_band_min_net = (broad_type[ispin] == metal ? _METALEDGE_ : vCBB[ispin][vimin_CBBs[ispin]]);
        Egap_net = (broad_type[ispin] == metal ? _METALGAP_ : gap[ispin]);
        Egap_type_net = (broad_type[ispin] == metal ? "metal" : "insulator-" + string(insulator_type[ispin] == insulator_indirect ? "in" : "") + "direct" + string(gap_type[ispin] == zero_gap ? "_zero-gap" : ""));
      } else {  // broad_type[1]==empty, grab properties of spin-channel 1
        const uint ispin = 0;
        valence_band_max_net = (broad_type[ispin] == metal ? _METALEDGE_ : vVBT[ispin][vimax_VBTs[ispin]]);
        conduction_band_min_net = (broad_type[ispin] == metal ? _METALEDGE_ : vCBB[ispin][vimin_CBBs[ispin]]);
        Egap_net = (broad_type[ispin] == metal ? _METALGAP_ : gap[ispin]);
        Egap_type_net = (broad_type[ispin] == metal ? "metal" : "insulator-" + string(insulator_type[ispin] == insulator_indirect ? "in" : "") + "direct" + string(gap_type[ispin] == zero_gap ? "_zero-gap" : ""));
      }
    }
  } else if (!(ISPIN == 2 && !(broad_type[0] == metal && broad_type[0] == broad_type[1]))) { // easy cases
    if (ISPIN == 1) {  // all ISPIN==1 come here
      const uint ispin = 0;
      valence_band_max_net = (broad_type[ispin] == metal ? _METALEDGE_ : vVBT[ispin][vimax_VBTs[ispin]]);
      conduction_band_min_net = (broad_type[ispin] == metal ? _METALEDGE_ : vCBB[ispin][vimin_CBBs[ispin]]);
      Egap_net = (broad_type[ispin] == metal ? _METALGAP_ : gap[ispin]);
      Egap_type_net = (broad_type[ispin] == metal ? "metal" : "insulator-" + string(insulator_type[ispin] == insulator_indirect ? "in" : "") + "direct" + string(gap_type[ispin] == zero_gap ? "_zero-gap" : ""));
    } else { // special case ISPIN==2 where both are metallic
      valence_band_max_net = _METALEDGE_;
      conduction_band_min_net = _METALEDGE_;
      Egap_net = _METALGAP_;
      Egap_type_net = "metal";
    }
  } else if (broad_type[0] == metal || broad_type[1] == metal) {  // special case ISPIN==2, one metallic and one insulating spin-channel means the whole system is half-metallic
    valence_band_max_net = _METALEDGE_;
    conduction_band_min_net = _METALEDGE_;
    Egap_net = _METALGAP_;
    Egap_type_net = "half-metal";
  } else if (!(broad_type[0] == insulator && broad_type[1] == insulator)) {  // test of stupidity
    message << "Unknown material type (!insulator_spin-polarized) (spin-averaged)";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
  } else {  // do more work for 2 spin channels that are both insulating
    // get spin-averaged properties
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " we have a spin-polarized insulator, need to find spin-averaged gap" << endl;
    }
    double max_VBT_net;
    double min_CBB_net;
    max_VBT_net = (-1.0) * AUROSTD_MAX_DOUBLE;
    min_CBB_net = AUROSTD_MAX_DOUBLE;
    uint imax_VBT_net = 0;
    uint imin_CBB_net = 0;
    uint ispin_VBT_net = 0;
    uint ispin_CBB_net = 0;
    vector<uint> vimax_VBTs_equiv_net;
    vector<uint> vimin_CBBs_equiv_net; // within tolerance of absolutes
    vector<uint> vispin_VBTs_net;
    vector<uint> vispin_CBBs_net;
    // first, get absolute max of VBT/min of CBB
    // these will be the max(vVBT[ispin]) and min(vCBB[ispin]) across spins from before
    for (uint ispin = 0; ispin < (uint) ISPIN; ispin++) {
      for (size_t ikpt = 0; ikpt < vkpoints[ispin].size(); ikpt++) {
        if (vVBT[ispin][ikpt] > max_VBT_net) {
          imax_VBT_net = ikpt;
          ispin_VBT_net = ispin;
          max_VBT_net = vVBT[ispin_VBT_net][imax_VBT_net];
        }
        if (vCBB[ispin][ikpt] < min_CBB_net) {
          imin_CBB_net = ikpt;
          ispin_CBB_net = ispin;
          min_CBB_net = vCBB[ispin_CBB_net][imin_CBB_net];
        }
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " absolute VBT_net=" << max_VBT_net << " at k-point " << vkpoints[ispin_VBT_net][imax_VBT_net] << " (kpt=" << imax_VBT_net << ",spin=" << ispin_VBT_net << ")" << endl;
      cerr << __AFLOW_FUNC__ << " absolute CBB_net=" << min_CBB_net << " at k-point " << vkpoints[ispin_VBT_net][imin_CBB_net] << " (kpt=" << imin_CBB_net << ",spin=" << ispin_CBB_net << ")" << endl;
    }
    vimax_VBTs_equiv_net.clear();
    vispin_VBTs_net.clear();
    vimin_CBBs_equiv_net.clear();
    vispin_CBBs_net.clear();
    // grab any equivalently high/low extrema
    for (uint ispin = 0; ispin < (uint) ISPIN; ispin++) {
      for (size_t ikpt = 0; ikpt < vkpoints[ispin].size(); ikpt++) {
        if (aurostd::abs(max_VBT_net - vVBT[ispin][ikpt]) < energy_tol) {
          vimax_VBTs_equiv_net.push_back(ikpt);
          vispin_VBTs_net.push_back(ispin);
          if (LDEBUG && ikpt != imax_VBT_net) {
            cerr << __AFLOW_FUNC__ << " found equivalent VBT_net at " << vkpoints[ispin][ikpt] << " (kpt=" << ikpt << ",spin=" << ispin << ")" << endl;
          }
        }
        if (aurostd::abs(min_CBB_net - vCBB[ispin][ikpt]) < energy_tol) {
          vimin_CBBs_equiv_net.push_back(ikpt);
          vispin_CBBs_net.push_back(ispin);
          if (LDEBUG && ikpt != imin_CBB_net) {
            cerr << __AFLOW_FUNC__ << " found equivalent CBB_net at " << vkpoints[ispin][ikpt] << " (kpt=" << ikpt << ",spin=" << ispin << ")" << endl;
          }
        }
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " removing duplicate k-points from VBT search (spin-averaged)" << endl;
    }
    vbt_duplicate_remove = removeDuplicateKPoints(vkpoints, vimax_VBTs_equiv_net, vispin_VBTs_net);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " removing duplicate k-points from CBB search (spin-averaged)" << endl;
    }
    cbb_duplicate_remove = removeDuplicateKPoints(vkpoints, vimin_CBBs_equiv_net, vispin_CBBs_net);
    if (!vbt_duplicate_remove || !cbb_duplicate_remove) {
      message << "Cannot find equivalent band extrema (spin-averaged)";
      if (force_exit) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      } else {
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
        return false;
      }
    }
    // if we found more than one possible extrema, minimize kpoint distance between max/min
    // this simulates the electron trying to reduce momentum needed to conduct
    if (vimax_VBTs_equiv_net.size() == 1 && vimin_CBBs_equiv_net.size() == 1) {  // easy case, keep already defined imax/imin
      dist_min = minimumDistanceKPoints(vkpoints[vispin_VBTs_net[0]][vimax_VBTs_equiv_net[0]], vkpoints[vispin_CBBs_net[0]][vimin_CBBs_equiv_net[0]]);
    } else { // find minimum distance pairs of max/min
      dist_min = AUROSTD_MAX_DOUBLE;
      for (size_t imax = 0; imax < vimax_VBTs_equiv_net.size(); imax++) {
        for (size_t imin = 0; imin < vimin_CBBs_equiv_net.size(); imin++) {
          dist = minimumDistanceKPoints(vkpoints[vispin_VBTs_net[imax]][vimax_VBTs_equiv_net[imax]], vkpoints[vispin_CBBs_net[imin]][vimin_CBBs_equiv_net[imin]]);
          if (dist < dist_min) {
            dist_min = dist;
            imax_VBT_net = vimax_VBTs_equiv_net[imax];
            ispin_VBT_net = vispin_VBTs_net[imax];
            imin_CBB_net = vimin_CBBs_equiv_net[imin];
            ispin_CBB_net = vispin_CBBs_net[imin];
          }
        }
      }
    }
    if (dist_min == AUROSTD_MAX_DOUBLE) {
      message << "Cannot calculate k-point distance between band extrema (spin-averaged)";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    valence_band_max_net = vVBT[ispin_VBT_net][imax_VBT_net];
    conduction_band_min_net = vCBB[ispin_CBB_net][imin_CBB_net];
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vCBB[ispin_CBB_net=" << ispin_CBB_net << "][imin_CBB_net=" << imin_CBB_net << "]=" << vCBB[ispin_CBB_net][imin_CBB_net] << endl;
      cerr << __AFLOW_FUNC__ << " vVBT[ispin_VBT_net=" << ispin_VBT_net << "][imax_VBT_net=" << imax_VBT_net << "]" << vVBT[ispin_VBT_net][imax_VBT_net] << endl;
    }
    const bool direct_insulator_net = aurostd::abs(dist_min) < kpt_tol;
    const double gap_net = (conduction_band_min_net - valence_band_max_net);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " gap_net=" << gap_net << endl;
    }
    if (gap_net < 0) { // test of stupidity, this should NEVER happen
      message << "Negative band gap found, something broke (spin-averaged)";
      if (force_exit) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      } else {
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
        return false;
      }
    }
    const bool zero_gap_net = gap_net < energy_tol;
    const bool spin_polarized_net = ispin_VBT_net != ispin_CBB_net;
    Egap_net = (zero_gap_net ? 0.0 : gap_net);
    Egap_type_net = "insulator-" + string(direct_insulator_net ? "" : "in") + "direct" + string(zero_gap_net ? "_zero-gap" : "") + string(spin_polarized_net ? "_spin-polarized" : "");
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ";
      cerr << (direct_insulator_net ? "" : "IN") << "DIRECT insulator ";
      cerr << "gap=" << gap_net << " ";
      cerr << (zero_gap_net ? "(zero-gap) " : "");
      cerr << (spin_polarized_net ? "(spin_polarized) " : "");
      cerr << "spin-averaged" << endl;
    }
  }
  Egap_fit_net = (aurostd::substring2bool(Egap_type_net, "empty") || aurostd::substring2bool(Egap_type_net, "metal") ? _METALGAP_ : 1.348 * Egap_net + 0.913);

  return true;
}

/// @brief calculates total dielectric function, EELS, and reflectivity
/// @note See A. M. Fox, Optical Properties of Solids
/// @authors
/// @mod{NHA,20251102,created function}
void xOUTCAR::GetDielectricData() {
  for (size_t i = 0; i < freq_grid.size(); i++) {
    xmatrix<double> dielectric_real;
    xmatrix<double> dielectric_imag;
    double dielectric_iso_real;
    double dielectric_iso_imag;
    double dielectric_norm;
    double refractive_index_iso;
    double extinction_coeff_iso;

    //full:
    dielectric_real = dielectric_interband_real[i] + dielectric_drude_real[i];
    dielectric_imag = dielectric_interband_imag[i] + dielectric_drude_imag[i];
    dielectric_full_real.push_back(dielectric_real);
    dielectric_full_imag.push_back(dielectric_imag);
    //isotropic:
    dielectric_iso_real = dielectric_interband_iso_real[i] + dielectric_drude_iso_real[i];
    dielectric_iso_imag = dielectric_interband_iso_imag[i] + dielectric_drude_iso_imag[i];
    dielectric_full_iso_real.push_back(dielectric_iso_real);
    dielectric_full_iso_imag.push_back(dielectric_iso_imag);
    dielectric_norm = std::sqrt(std::pow(dielectric_iso_real, 2.0) + std::pow(dielectric_iso_imag, 2.0));
    energy_loss_function_iso.emplace_back(dielectric_iso_imag / std::pow(dielectric_norm, 2.0));
    refractive_index_iso = std::sqrt(0.5 * (dielectric_norm + dielectric_iso_real));
    extinction_coeff_iso = std::sqrt(0.5 * (dielectric_norm - dielectric_iso_real));
    reflectivity_iso.emplace_back((std::pow(refractive_index_iso - 1.0, 2.) + std::pow(extinction_coeff_iso, 2.0)) / (std::pow(refractive_index_iso + 1.0, 2.) + std::pow(extinction_coeff_iso, 2.0)));
  }
}

/// @brief gets the optical information
/// @param freq_relax relaxation frequency of the electrons
/// @return true if successful
/// @note See A. M. Fox, Optical Properties of Solids
/// @authors
/// @mod{SD,20250407,created function}
bool xOUTCAR::GetOptical(const double freq_relax) {
  dielectric_drude_real.clear();
  dielectric_drude_imag.clear();
  dielectric_drude_iso_real.clear();
  dielectric_drude_iso_imag.clear();
  if (aurostd::isequal(aurostd::sum(aurostd::abs(freq_plasma)), 0.0)) {
    aurostd::xmatrix<double> dummy(3, 3);
    dielectric_drude_real.resize(freq_grid.size(), dummy);
    dielectric_drude_imag.resize(freq_grid.size(), dummy);
    dielectric_drude_iso_real.resize(freq_grid.size(), 0.0);
    dielectric_drude_iso_imag.resize(freq_grid.size(), 0.0);
  } else {
    freq_plasma_iso = 0.0;
    for (int i = 1; i <= 3; i++) {
      freq_plasma_iso += freq_plasma[i][i] / 3.0;
    }
    xmatrix<double> identity_matrix = aurostd::identity((double) 0, 3);

    const double denergy = freq_grid[1] - freq_grid[0];
    for (const double& freq : freq_grid) {
      xmatrix<double> drude = (freq_plasma * freq_plasma) / (std::pow(freq, 2.0) + std::pow(freq_relax, 2.0));
      double drude_iso = std::pow(freq_plasma_iso, 2.0) / (std::pow(freq, 2.0) + std::pow(freq_relax, 2.0));
      dielectric_drude_real.emplace_back(identity_matrix - drude);
      dielectric_drude_imag.emplace_back(drude * (freq_relax / std::max(freq, denergy)));
      dielectric_drude_iso_real.emplace_back(1.0 - drude_iso);
      dielectric_drude_iso_imag.emplace_back(drude_iso * (freq_relax / std::max(freq, denergy)));
    }
  }

  double dielectric_iso;
  for (const aurostd::xmatrix<double>& dielectric_interband : dielectric_interband_real) {
    dielectric_iso = 0.0;
    for (int i = 1; i <= 3; i++) {
      dielectric_iso += dielectric_interband[i][i] / 3.0;
    }
    dielectric_interband_iso_real.emplace_back(dielectric_iso);
  }
  for (const aurostd::xmatrix<double>& dielectric_interband : dielectric_interband_imag) {
    dielectric_iso = 0.0;
    for (int i = 1; i <= 3; i++) {
      dielectric_iso += dielectric_interband[i][i] / 3.0;
    }
    dielectric_interband_iso_imag.emplace_back(dielectric_iso);
  }

  GetDielectricData();
  return true;
}

/// @brief gets the ionic step data
/// @authors
/// @mod{CO,20211106,created function}
/// @mod{SD,20221208,patched for VASP5 and VASP6}
/// @mod{SD,20241902,fixed stress output}
bool xOUTCAR::GetIonicStepsData() {
  const bool LDEBUG = (false || XHOST.DEBUG);
  vxstr_ionic.clear();
  venergy_ionic.clear();
  vstresses_ionic.clear();

  bool reading_ionic = false;
  bool reading_stresses = false;
  bool reading_lattice = false;
  bool reading_atoms = false;
  bool convert_kBar2eV = false;

  // get species data first
  uint iline = 0;
  uint i = 0;
  vector<string> vtokens;
  string tmp_str;
  deque<string> species_pp;
  deque<int> num_each_type;
  const uint vcontent_size = vcontent.size();
  for (iline = 0; iline < vcontent_size; iline++) {
    if (aurostd::substring2bool(vcontent[iline], "POTCAR:")) {
      tmp_str = vcontent[iline].substr(vcontent[iline].find("POTCAR:") + 6 + 1);  // len(POSCAR:)==6
      aurostd::string2tokens(tmp_str, vtokens, " ");
      // vtokens[1] is the species_pp
      if (!aurostd::WithinList(species_pp, vtokens[1])) {
        species_pp.push_back(vtokens[1]);
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " vcontent[iline=" << iline << "]=\"" << vcontent[iline] << "\"" << endl;
      }
    }
    if (aurostd::substring2bool(vcontent[iline], "ions") && // THIS MUST COME AFTER species_pp IS FILLED
        aurostd::substring2bool(vcontent[iline], "per") && aurostd::substring2bool(vcontent[iline], "type")) {
      aurostd::string2tokens(vcontent[iline], vtokens, "=");
      if (vtokens.size() > 1) {
        tmp_str = vtokens[1];
        aurostd::string2tokens(tmp_str, vtokens, " ");
        if (!vtokens.empty() && vtokens.size() == species_pp.size()) {
          num_each_type = aurostd::vector2deque(aurostd::vectorstring2vectorutype<int>(vtokens));
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " vcontent[iline=" << iline << "]=\"" << vcontent[iline] << "\"" << endl;
          }
          break;
        }
      }
    }
  }
  const uint natoms = aurostd::sum(num_each_type);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species_pp=" << aurostd::joinWDelimiter(species_pp, ",") << endl;
    cerr << __AFLOW_FUNC__ << " num_each_type=" << aurostd::joinWDelimiter(num_each_type, ",") << endl;
    cerr << __AFLOW_FUNC__ << " natoms=" << natoms << endl;
  }

  uint ilattice = 0;
  uint iatom = 0;
  uint itype = 0;
  xstructure xstr;
  double energy = AUROSTD_MAX_DOUBLE;
  xvector<double> stresses(6);
  stresses[stresses.lrows] = AUROSTD_MAX_DOUBLE;
  vector<_atom> vatoms; // must keep like this until we settle types...
  for (iline = 0; iline < vcontent_size; iline++) {
    if (reading_ionic == false) {
      if (aurostd::substring2bool(vcontent[iline], "aborting") && aurostd::substring2bool(vcontent[iline], "loop") && aurostd::substring2bool(vcontent[iline], "because") &&
          aurostd::substring2bool(vcontent[iline], "EDIFF") && aurostd::substring2bool(vcontent[iline], "is") && aurostd::substring2bool(vcontent[iline], "reached")) {
        reading_ionic = true;
        reading_stresses = false;
        reading_lattice = false;
        reading_atoms = false;
        convert_kBar2eV = false;
        // clear everything, set stresses to 3 as an indicator that it's not set
        xstr.clear();
        vatoms.clear();
        ilattice = 0;
        iatom = 0;
        energy = AUROSTD_MAX_DOUBLE;
        stresses[stresses.lrows] = AUROSTD_MAX_DOUBLE; // stresses.resize(3);
        continue;
      }
    } else {
      if (aurostd::substring2bool(vcontent[iline], "Iteration") || (iline == vcontent_size - 1)) {
        reading_ionic = false;
        reading_stresses = false;
        if (xstr.atoms.size() == natoms && ilattice == 3 && iatom == natoms && energy != AUROSTD_MAX_DOUBLE && stresses[stresses.lrows] != AUROSTD_MAX_DOUBLE) { // add additional checks as needed //stresses.rows==6
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " adding a new ionic step" << endl;
          }
          vxstr_ionic.push_back(xstr);
          venergy_ionic.push_back(energy);
          vstresses_ionic.push_back(stresses);
        }
        continue;
      }
      if (aurostd::substring2bool(vcontent[iline], "energy") && aurostd::substring2bool(vcontent[iline], "without") && aurostd::substring2bool(vcontent[iline], "entropy")) {
        aurostd::string2tokens(vcontent[iline], vtokens, "=");
        if (vtokens.size() > 1) {
          tmp_str = vtokens[1];
          aurostd::string2tokens(tmp_str, vtokens, "e");
          if (!vtokens.empty()) {
            energy = aurostd::string2utype<double>(vtokens[0]);
          }
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " vcontent[iline=" << iline << "]=\"" << vcontent[iline] << "\"" << endl;
            cerr << __AFLOW_FUNC__ << " energy=" << energy << endl;
          }
        }
      }
      if (aurostd::substring2bool(vcontent[iline], "FORCE") && aurostd::substring2bool(vcontent[iline], "-STRESS")) {
        reading_stresses = true;
        reading_lattice = false;
        reading_atoms = false;
        convert_kBar2eV = false;
        continue;
      }
      if (reading_stresses) {
        if (aurostd::substring2bool(vcontent[iline], "Total")) {
          aurostd::string2tokens(vcontent[iline], vtokens, " ");
          if (vtokens.size() == 7) {
            // stresses.resize(6);
            for (i = 1; i < 7; i++) {
              stresses[stresses.lrows + i - 1] = -1.0 * aurostd::string2utype<double>(vtokens[i]);
            } // SD20241902 - stress is the negative of the pressure
            if (LDEBUG) {
              cerr << __AFLOW_FUNC__ << " vcontent[iline=" << iline << "]=\"" << vcontent[iline] << "\"" << endl;
              cerr << __AFLOW_FUNC__ << " stresses(eV)=" << stresses << endl; // X Y Z XY YZ ZX
            }
          }
        }
        // per ME's suggestion, let's grab 'in kB' instead of 'Total' (eV) which has more significant figures
        if (aurostd::substring2bool(vcontent[iline], "in") && aurostd::substring2bool(vcontent[iline], "kB")) {
          aurostd::string2tokens(vcontent[iline], vtokens, " ");
          if (vtokens.size() == 8) {
            // stresses.resize(6);
            for (i = 2; i < 8; i++) {
              stresses[stresses.lrows + i - 2] = -1.0 * aurostd::string2utype<double>(vtokens[i]);
            } // SD20241902 - stress is the negative of the pressure
            convert_kBar2eV = true;
            if (LDEBUG) {
              cerr << __AFLOW_FUNC__ << " vcontent[iline=" << iline << "]=\"" << vcontent[iline] << "\"" << endl;
              cerr << __AFLOW_FUNC__ << " stresses(kBar)=" << stresses << endl; // X Y Z XY YZ ZX
            }
          }
        }
      }
      if (aurostd::substring2bool(vcontent[iline], "direct") && aurostd::substring2bool(vcontent[iline], "lattice") && aurostd::substring2bool(vcontent[iline], "vectors") &&
          aurostd::substring2bool(vcontent[iline], "reciprocal")) {
        reading_lattice = true;
        reading_stresses = false;
        reading_atoms = false;
        ilattice = 0;
        xstr.clear();
        continue;
      }
      if (reading_lattice && ilattice < 3) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " vcontent[iline=" << iline << "]=\"" << vcontent[iline] << "\"" << endl;
        }
        vtokens = GetCorrectEntriesFromLine(vcontent[iline], 6); // aurostd::string2tokens(vcontent[iline],vtokens," ");
        if (vtokens.size() != 6) {
          ilattice = 0;
          reading_lattice = false;
          continue;
        } // also contains reciprocal lattice components //ilattice==0 is an indicator that lattice was not read (correctly)
        xstr.lattice[xstr.lattice.lrows + ilattice][xstr.lattice.lcols + 0] = aurostd::string2utype<double>(vtokens[0]);
        xstr.lattice[xstr.lattice.lrows + ilattice][xstr.lattice.lcols + 1] = aurostd::string2utype<double>(vtokens[1]);
        xstr.lattice[xstr.lattice.lrows + ilattice][xstr.lattice.lcols + 2] = aurostd::string2utype<double>(vtokens[2]);
        ilattice++;
        if (ilattice == 3) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " xstr.lattice=" << endl << xstr.lattice << endl;
            cerr << __AFLOW_FUNC__ << " xstr.GetVolume()=" << xstr.GetVolume() << endl;
          }
          xstr.scale = 1.0;
          xstr.FixLattices();
          if (convert_kBar2eV) {
            stresses *= kBar2eV;
          } else {
            stresses /= xstr.GetVolume();
          }
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " stresses(eV/Ang^3)=" << stresses << endl; // X Y Z XY YZ ZX //this is eV per volume of cell
          }
          reading_lattice = false;
        }
      }
      if (aurostd::substring2bool(vcontent[iline], "POSITION") && aurostd::substring2bool(vcontent[iline], "TOTAL-FORCE")) {
        reading_atoms = true;
        reading_stresses = false;
        reading_lattice = false;
        convert_kBar2eV = false;
        iatom = 0;
        vatoms.clear();
        continue;
      }
      if (reading_atoms && iatom < natoms) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " vcontent[iline=" << iline << "]=\"" << vcontent[iline] << "\"" << endl;
        }
        if (aurostd::substring2bool(vcontent[iline], "---------")) {
          continue;
        }
        vtokens = GetCorrectEntriesFromLine(vcontent[iline], 6); // aurostd::string2tokens(vcontent[iline],vtokens," ");
        if (vtokens.size() != 6) {
          iatom = 0;
          reading_atoms = false;
          continue;
        } // also contains force components //iatom==0 is an indicator that atoms were not read (correctly)
        vatoms.emplace_back();
        for (i = 0; i < 3; i++) {
          vatoms.back().cpos[vatoms.back().cpos.lrows + i] = aurostd::string2utype<double>(vtokens[i]);
        }
        for (i = 0; i < 3; i++) {
          vatoms.back().force[vatoms.back().force.lrows + i] = aurostd::string2utype<double>(vtokens[i + 3]);
        }
        iatom++;
        if (iatom == natoms) {
          iatom = 0; // reset
          for (itype = 0; itype < num_each_type.size(); itype++) {
            for (i = 0; i < (uint) num_each_type[itype]; i++) {
              vatoms[iatom].type = itype;
              vatoms[iatom].name = species_pp[itype];
              vatoms[iatom].name_is_given = (!vatoms[iatom].name.empty());
              vatoms[iatom].fpos = xstr.c2f * vatoms[iatom].cpos;
              xstr.AddAtom(vatoms[iatom], false); // CO20230319 - add by type
              xstr.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
              iatom++;
            }
          }
          xstr.title = SYSTEM;
          xstr.MakeBasis();
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " xstr=" << endl << xstr << endl;
          }
          reading_atoms = false;
        }
      }
    }
  }

  return true;
}

/// @brief write ionic step data to JSON
/// @param jo JSON object
/// @param entry aflowlib entry
/// @authors
/// @mod{CO,20211106,created function}
/// @mod{SD+HE,20221208,rewritten using JSON objects}
void xOUTCAR::AddStepsIAPCFG(aurostd::JSON::object& jo, aflowlib::_aflowlib_entry& entry) {
  if (vxstr_ionic.size() != venergy_ionic.size()) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "vxstr_ionic.size()!=venergy_ionic", _FILE_CORRUPT_);
  }
  if (vxstr_ionic.size() != vstresses_ionic.size()) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "vxstr_ionic.size()!=vstresses_ionic", _FILE_CORRUPT_);
  }
  // bool FORMATION_CALC = false;
  vector<string> type;
  vector<xvector<double>> fpos;
  vector<xvector<double>> cpos;
  vector<xvector<double>> force;
  for (size_t istr = 0; istr < vxstr_ionic.size(); istr++) {
    aurostd::JSON::object step(aurostd::JSON::object_types::DICTIONARY);
    step["auid"] = entry.auid;
    step["energy"] = venergy_ionic[istr];
    step["lattice"] = vxstr_ionic[istr].lattice;
    step["stress"] = vstresses_ionic[istr];
    type.clear();
    fpos.clear();
    cpos.clear();
    force.clear();
    for (size_t iatom = 0; iatom < vxstr_ionic[istr].atoms.size(); iatom++) {
      type.push_back(vxstr_ionic[istr].atoms[iatom].cleanname);
      fpos.push_back(vxstr_ionic[istr].atoms[iatom].fpos);
      cpos.push_back(vxstr_ionic[istr].atoms[iatom].cpos);
      force.push_back(vxstr_ionic[istr].atoms[iatom].force);
    }
    aurostd::JSON::object atoms(aurostd::JSON::object_types::DICTIONARY);
    atoms["type"] = type;
    atoms["fpos"] = fpos;
    atoms["cpos"] = cpos;
    atoms["force"] = force;
    step["atoms"] = atoms;
    jo["data"].push_back(step);
  }
}

aurostd::JSON::object xOUTCAR::serialize() const {
  return aurostd::JSON::object({AST_JSON_GETTER(JSON_xOUTCAR_MEMBERS)});
}

xOUTCAR xOUTCAR::deserialize(const aurostd::JSON::object& jo) {
  AST_JSON_SETTER(JSON_xOUTCAR_MEMBERS)
  return *this;
}

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************

// ***************************************************************************
// class xDOSCAR
bool xDOSCAR::GetProperties(const string& stringIN, bool QUIET) {
  stringstream sss;
  sss.str(stringIN);
  if (filename.empty()) {
    filename = "string";
  }
  return xDOSCAR::GetProperties(sss, QUIET);
}

bool xDOSCAR::GetPropertiesFile(const string& fileIN, bool QUIET) {
  stringstream sss;
  filename = fileIN;
  aurostd::compressfile2stringstream(fileIN, sss);
  return xDOSCAR::GetProperties(sss, QUIET);
}

bool xDOSCAR::GetPropertiesUrlFile(const string& url, const string& file, bool QUIET) {
  const string tmpfile = aurostd::TmpFileCreate("xDOSCAR_GetProperties"); // CO20200502 - threadID
  aurostd::httpGetFileStatus(url + "/" + file, tmpfile);
  const bool out = GetPropertiesFile(tmpfile, QUIET);
  filename = "url=" + url; // CO20210315
  aurostd::RemoveFile(tmpfile);
  return out;
}

bool xDOSCAR::GetProperties(const stringstream& stringstreamIN, bool QUIET) {
  const bool LDEBUG = (false || XHOST.DEBUG || !QUIET);
  const bool force_exit = XHOST.POSTPROCESS; // SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only
  stringstream message;

  const bool ERROR_flag = false;
  const long double seconds = aurostd::get_seconds();
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BEGIN (" << time_delay(seconds) << ")" << endl;
  }
  clear(); // so it does not mess up vector/deque
  content = stringstreamIN.str();
  vcontent.clear();
  const vector<string> vline;
  vector<string> tokens;
  aurostd::string2vectorstring(content, vcontent);
  uint ndos = 1; // ME20190614
  if (filename.empty()) {
    filename = "stringstream";
  }
  // ME20190812 - Add checks for broken DOSCARs
  if (vcontent.size() < 7) {
    message << "Broken DOSCAR: no content";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
  }
  for (uint iline = 0; iline < 7; iline++) { // ME20190614 - Read header
    aurostd::string2tokens(vcontent[iline], tokens);
    // cerr << "iline=" << iline << "  " << vcontent[iline] << " tokens.size()=" << tokens.size() << endl;
    // WRONG if(iline==0 && tokens.size()==4)  spin=aurostd::string2utype<int>(tokens.at(tokens.size()-1))-1;
    // ME20190614 START
    if (iline == 0 && tokens.size() >= 4) {
      number_atoms = aurostd::string2utype<uint>(tokens[1]);
      partial = aurostd::string2utype<bool>(tokens[2]);
      if (partial) {
        ndos += number_atoms;
      }
    }
    // ME20190614 END
    if (iline == 1 && tokens.size() >= 5) {
      uint i = 0;
      Vol = aurostd::string2utype<double>(tokens.at(i++));
      lattice(1) = aurostd::string2utype<double>(tokens.at(i++));
      lattice(2) = aurostd::string2utype<double>(tokens.at(i++));
      lattice(3) = aurostd::string2utype<double>(tokens.at(i++));
      POTIM = aurostd::string2utype<double>(tokens.at(i++));
    }
    if (iline == 2) {
      temperature = aurostd::string2utype<double>(vcontent[iline]);
    }
    if (iline == 3) {
      carstring = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vcontent[iline]); // ME20190620 - what kind of DOSCAR?
    }
    if (iline == 4) {
      title = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vcontent[iline]); // ME20190508 - clean title
    }
    if (iline == 5) {
      // cerr << "iline=" << iline << "  " << vcontent[iline] << " tokens.size()=" << tokens.size() << endl;
      energy_max = aurostd::string2utype<double>(tokens.at(0));
      energy_min = aurostd::string2utype<double>(tokens.at(1));
      number_energies = aurostd::string2utype<uint>(tokens.at(2));
      Efermi = aurostd::string2utype<double>(tokens.at(3));
    }
    if (iline == 6 && tokens.size() == 4) {
      spin = 0;
      RWIGS = true;
    }
    if (iline == 6 && tokens.size() == 7) {
      spin = 1;
      RWIGS = true;
    }
    // ME20190614 START
    if (iline == 6) {
      if (tokens.size() == 3) {
        spin = 0;
      }
      if (tokens.size() == 5) {
        spin = 1;
      }
    }
    // ME20190614 END
  }
  // Done reading header
  // ME20190614 START
  // ME20200305 - Check if DOSCAR is broken
  uint number_lines = number_energies + 6;
  if (partial) {
    number_lines += number_atoms * (number_energies + 1);
  }
  if (vcontent.size() < number_lines) {
    message << "Broken DOSCAR - not enough lines.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
  }
  uint norbitals = 0;
  int d = 0;
  int e = 0;
  double dos = 0.0;
  if (partial) {
    aurostd::string2tokens(vcontent.at(number_energies + 7), tokens, " ");
    const int ncols = (int) tokens.size() - 1; // Don't count the energy column
    if (carstring == "PHON") { // APL DOSCAR
      norbitals = ncols;
      isLSCOUPLING = false;
      lmResolved = true;
    } else { // No special case, so VASP DOSCAR
      // Determine whether the DOSCAR is lm-resolved and/or has spin-orbit coupling (SOC)
      if (ncols == 16) {
        // If the DOSCAR has 16 columns, the variables cannot be properly resolved.
        // It could either lm-resolved with f-electrons and without spin polarization,
        // or it is not lm-resolved with f-electrons and SOC. We need another input file
        // to resolve this. The former case is VERY unlikely since calculations with
        // f-electrons typically require spin polarization. So, SOC is the default.
        // First, try and find filename extensions
        vector<string> vstr;
        string ext;
        if (aurostd::substring2bool(filename, "DOSCAR.")) {
          aurostd::string2tokens(filename, vstr, ".");
          for (size_t i = 1; i < vstr.size(); i++) {
            ext += "." + vstr[i];
          }
        } else if (aurostd::substring2bool(filename, "DOSCAR_")) {
          aurostd::string2tokens(filename, vstr, "_");
          for (size_t i = 1; i < vstr.size(); i++) {
            ext += "_" + vstr[i];
          }
        }
        vstr.clear();
        // Try INCAR first because it's the smallest file
        if (aurostd::CompressFileExist("INCAR" + ext) || aurostd::FileExist("INCAR" + ext)) {
          aurostd::compressfile2vectorstring("INCAR" + ext, vstr);
          aurostd::RemoveComments(vstr);
          for (size_t i = 0; i < vstr.size(); i++) {
            if (aurostd::substring2bool(vstr[i], "LSORBIT")) {
              vector<string> tokens;
              aurostd::string2tokens(vstr[i], tokens, " =");
              tokens[1] = aurostd::toupper(tokens[1]);
              if ((tokens[1][0] == 'T') || (aurostd::substring2bool(tokens[1], "true"))) {
                isLSCOUPLING = true;
              } else {
                isLSCOUPLING = false;
              }
              break;
            }
          }
          // No INCAR found. Try vasprun.xml
        } else if (aurostd::CompressFileExist("vasprun.xml" + ext) || aurostd::FileExist("vasprun.xml" + ext)) {
          aurostd::compressfile2vectorstring("vasprun.xml" + ext, vstr);
          for (size_t i = 0; i < vstr.size(); i++) {
            if (aurostd::substring2bool(vstr[i], "LSORBIT")) {
              isLSCOUPLING = !(aurostd::substring2bool(vstr[i], "F")); // Only contains "F" if LSORBIT is false
              break;
            }
          }
          // At last, try OUTCAR
        } else if (aurostd::CompressFileExist("OUTCAR" + ext) || aurostd::FileExist("OUTCAR" + ext)) {
          aurostd::compressfile2vectorstring("vasprun.xml" + ext, vstr);
          for (size_t i = 0; i < vstr.size(); i++) {
            if (aurostd::substring2bool(vstr[i], "LSORBIT")) {
              isLSCOUPLING = !(aurostd::substring2bool(vstr[i], "F")); // Only contains "F" if LSORBIT is false
              break;
            }
          }
          // Nothing found, assume SOC (more likely case)
        } else {
          const string message =
              "Could not determine whether the DOSCAR is lm-resolved"
              " or contains spin-orbit coupling. AFLOW will assume that the"
              " DOSCAR is lm-resolved. If this is not the case, please put an"
              " INCAR" +
              ext + ", a vasprun.xml " + ext +
              ", or an"
              " OUTCAR" +
              ext + " file into the working directory and try again.";
          pflow::logger(__AFLOW_FILE__, "xDOSCAR::GetProperties()", message, std::cerr, _LOGGER_WARNING_, QUIET);
          isLSCOUPLING = true;
        }
        lmResolved = !(isLSCOUPLING); // With 16 columns, it cannot be both
      } else {
        // Otherwise, the number of columns per spin smaller than 10 without
        // SOC, except for lm-resolved DOS with f-electrons (32 columns)
        isLSCOUPLING = ((ncols / (spin + 1) > 10) && (ncols != 32));
      }
      if (isLSCOUPLING) {
        norbitals = ncols / 4;
      } else {
        norbitals = (ncols) / (spin + 1);
      }
      // Since VASP always prints s, p, and d orbitals, lm-resolved
      // DOSCARS contain at least 9 oribtals
      lmResolved = (norbitals > 8);
    }
  }
  if (isLSCOUPLING) { // ME20190620 - LSCOUPLING has four spin channels
    vDOS.assign(ndos, deque<deque<deque<double>>>(norbitals + 1, deque<deque<double>>(4, deque<double>(number_energies, 0.0))));
  } else {
    vDOS.assign(ndos, deque<deque<deque<double>>>(norbitals + 1, deque<deque<double>>(spin + 1, deque<double>(number_energies, 0.0))));
  }
  venergy.resize(number_energies);
  viDOS.assign(spin + 1, deque<double>(number_energies, 0.0));
  // Read data
  for (size_t iline = 6; iline < vcontent.size(); iline++) {
    if (iline == (d + 1) * number_energies + 6 + d) {
      d++;
      // ME20190810 - Safeguard against DOSCARs with additional lines
      if (d == (int) ndos) {
        message << "DOSCAR contains more lines than the header suggests." << endl;
        message << "xDOSCAR object may not be properly populated.";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_oss, _LOGGER_WARNING_, QUIET);
        break;
      }
      e = 0;
    } else {
      aurostd::string2tokens(vcontent[iline], tokens);
      if (d == 0) {
        venergy[e] = aurostd::string2utype<double>(tokens[0]);
        for (uint i = 0; i < (spin + 1); i++) {
          vDOS[0][0][i][e] = aurostd::string2utype<double>(tokens[i + 1]);
          viDOS[i][e] = aurostd::string2utype<double>(tokens[i + spin + 2]);
        }
        e++;
      } else if (isLSCOUPLING) { // ME20190620 - LSCOUPLING has four spin channels
        for (uint i = 0; i < 4 * norbitals; i++) {
          dos = aurostd::string2utype<double>(tokens[i + 1]);
          vDOS[d][i / 4 + 1][i % 4][e] = dos;
          vDOS[d][0][i % 4][e] += dos;
          vDOS[0][i / 4 + 1][i % 4][e] += dos;
        }
        e++;
      } else {
        for (uint i = 0; i < (spin + 1) * norbitals; i++) {
          dos = aurostd::string2utype<double>(tokens[i + 1]);
          vDOS[d][i / (spin + 1) + 1][i % (spin + 1)][e] = dos;
          vDOS[d][0][i % (spin + 1)][e] += dos;
          vDOS[0][i / (spin + 1) + 1][i % (spin + 1)][e] += dos;
        }
        e++;
      }
      // ME20190614 - The procedure below makes too many assumptions about how the
      //  DOSCAR file is structured and fails for f-elements.
    }
    // ME20190614 END
  }
  // ME20190812 - Safeguard against broken DOSCARs
  if ((d + 1 < (int) ndos) || (e < (int) number_energies)) { // ME20191010: needs to be d + 1
    message << "Broken DOSCAR: not enough lines";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
  }
  // fix denergy
  denergy = venergy.at(1) - venergy.at(0);
  venergyEf.clear();
  for (size_t i = 0; i < venergy.size(); i++) {
    venergyEf.push_back(venergy[i] - Efermi);
  }

  // ----------------------------------------------------------------------
  // spin polarization at FERMI level

  double Fup = 0.0;
  double Fdown = 0.0;
  double Minup = 0.0;
  double Mindown = 0.0;
  double Maxup = 0.0;
  double Maxdown = 0.0; // CAMILOFIX
  bool Fermifound = false;
  bool firstenter = true;
  const double zeroTol = 1e-8;
  if (!spin) {
    spinF = 0.0;
  } else {
    for (size_t i = 6; i < vcontent.size(); i++) {
      double energytp;
      aurostd::string2tokens(vcontent[i], tokens, " ");
      if (!tokens.empty() && tokens.size() < 6) {
        energytp = aurostd::string2utype<double>(tokens.at(0));
        if (Efermi - energytp > 0) {
          Minup = aurostd::string2utype<double>(tokens.at(1));
          Mindown = aurostd::string2utype<double>(tokens.at(2));
        } else if (firstenter) {
          Maxup = aurostd::string2utype<double>(tokens.at(1));
          Maxdown = aurostd::string2utype<double>(tokens.at(2));
          firstenter = false;
          Fermifound = true;
        }
        if (Fermifound) {
          Fup = Minup + Maxup;
          Fdown = Mindown + Maxdown;
          break;
        }
      }
    }
    if ((Fup + Fdown) < zeroTol) {
      spinF = 0.0;
    } else {
      spinF = aurostd::abs((Fup - Fdown) / (Fup + Fdown)); // otherwise AFLOW_NAN
    }
  }

  // ----------------------------------------------------------------------
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " title=" << title << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " spin=" << spin << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " Vol=" << Vol << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " lattice=" << lattice << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " POTIM=" << POTIM << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " temperature=" << temperature << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " RWIGS=" << RWIGS << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " Efermi=" << Efermi << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " spinF=" << spinF << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " number_energies=" << number_energies << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " energy_max=" << energy_max << " energy_min=" << energy_min << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " denergy=" << denergy << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " venergy.size()=" << venergy.size() << " venergyEf.size()=" << venergyEf.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vDOS.size()=" << vDOS.size() << ", " << vDOS[0].size() << ", " << vDOS[0][0].size() << ", " << vDOS[0][0][0].size() << std::endl; // ME20190614 - new vDOS format
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " viDOS.size()=" << viDOS.size() << ", " << viDOS[0].size() << ", " << std::endl; // ME20190614 - new viDOS format
  }
  // if(LDEBUG) cerr << __AFLOW_FUNC__ << " vDOSs.at(max).size()=" << vDOSs.at(vDOSs.size()-1).size() << endl;
  // if(LDEBUG) cerr << __AFLOW_FUNC__ << " vDOSp.at(max).size()=" << vDOSp.at(vDOSp.size()-1).size() << endl;
  // if(LDEBUG) cerr << __AFLOW_FUNC__ << " vDOSd.at(max).size()=" << vDOSd.at(vDOSd.size()-1).size() << endl;
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " END (" << time_delay(seconds) << ")" << endl;
  }
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  if (ERROR_flag) {
    message << "ERROR_flag set in xDOSCAR";
    if (force_exit) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    } else {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_, QUIET);
    }
    return false;
  }
  return true;
}

// CO20191217 - copies everything from spin channel 1 to spin channel 2
void xDOSCAR::convertSpinOFF2ON() { // CO20191217
  const bool LDEBUG = (false || XHOST.DEBUG);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
  }

  uint iatom = 0;
  uint iorbital = 0;

  // check that it is truly SPIN-OFF
  if (viDOS.size() != 1) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "viDOS.size()!=1", _INPUT_ERROR_);
  }; // no conversion needed
  for (iatom = 0; iatom < vDOS.size(); iatom++) {
    for (iorbital = 0; iorbital < vDOS[iatom].size(); iorbital++) {
      if (vDOS[iatom][iorbital].size() != 1) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "vDOS[iatom][iorbital].size()!=1", _INPUT_ERROR_);
      }; // no conversion needed
    }
  }
  if (conduction_band_min.size() != 1) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "conduction_band_min.size()!=1", _INPUT_ERROR_);
  }; // no conversion needed
  if (valence_band_max.size() != 1) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "valence_band_max.size()!=1", _INPUT_ERROR_);
  }; // no conversion needed
  if (Egap.size() != 1) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Egap!=1", _INPUT_ERROR_);
  }; // no conversion needed
  if (Egap_fit.size() != 1) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Egap_fit!=1", _INPUT_ERROR_);
  }; // no conversion needed
  if (Egap_type.size() != 1) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Egap_type!=1", _INPUT_ERROR_);
  }; // no conversion needed
  if (spin != 0) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "spin!=0", _INPUT_ERROR_);
  }; // no conversion needed

  // copy everything over
  viDOS.push_back(viDOS.back());
  for (iatom = 0; iatom < vDOS.size(); iatom++) {
    for (iorbital = 0; iorbital < vDOS[iatom].size(); iorbital++) {
      vDOS[iatom][iorbital].push_back(vDOS[iatom][iorbital].back());
    }
  }
  conduction_band_min.push_back(conduction_band_min.back());
  valence_band_max.push_back(valence_band_max.back());
  Egap.push_back(Egap.back());
  Egap_fit.push_back(Egap_fit.back());
  Egap_type.push_back(Egap_type.back());
  spin = 1;
}

void xDOSCAR::addAtomChannel() { // CO20211124
  const bool LDEBUG = (false || XHOST.DEBUG);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
  }

  string ERROR_out;
  if (!checkDOS(ERROR_out)) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, ERROR_out, _INDEX_BOUNDS_);
  } // no conversion needed

  const uint atoms_size = vDOS.size();
  if (atoms_size == 0) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "vDOS.size()==0", _INDEX_BOUNDS_);
  } // no conversion needed
  const uint orbital_size = vDOS.front().size();
  uint spin_size = 0;
  if (orbital_size > 0) {
    spin_size = vDOS.front().front().size();
  }
  uint energy_size = 0;
  if (orbital_size > 0 && spin_size > 0) {
    energy_size = vDOS.front().front().front().size();
  }
  vDOS.emplace_back(0);
  uint iorbital = 0;
  uint ispin = 0;
  if (orbital_size > 0) {
    vDOS.back().resize(orbital_size);
  }
  if (spin_size > 0) {
    for (iorbital = 0; iorbital < orbital_size; iorbital++) {
      vDOS.back()[iorbital].resize(spin_size);
      if (energy_size > 0) {
        for (ispin = 0; ispin < spin_size; ispin++) {
          vDOS.back()[iorbital][ispin].resize(energy_size, 0.0);
        }
      }
    }
  }
  number_atoms += 1;
}

void xDOSCAR::addOrbitalChannel() { // CO20211124
  const bool LDEBUG = (false || XHOST.DEBUG);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
  }

  string ERROR_out;
  if (!checkDOS(ERROR_out)) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, ERROR_out, _INDEX_BOUNDS_);
  } // no conversion needed

  const uint atoms_size = vDOS.size();
  if (atoms_size == 0) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "vDOS.size()==0", _INDEX_BOUNDS_);
  } // no conversion needed
  const uint orbital_size = vDOS.front().size();
  uint spin_size = 0;
  if (orbital_size > 0) {
    spin_size = vDOS.front().front().size();
  }
  uint energy_size = 0;
  if (orbital_size > 0 && spin_size > 0) {
    energy_size = vDOS.front().front().front().size();
  }
  uint iatom = 0;
  uint ispin = 0;
  for (iatom = 0; iatom < atoms_size; iatom++) {
    vDOS[iatom].emplace_back(0);
    if (spin_size > 0) {
      vDOS[iatom].back().resize(spin_size);
    }
    if (energy_size > 0) {
      for (ispin = 0; ispin < spin_size; ispin++) {
        vDOS[iatom].back()[ispin].resize(energy_size, 0.0);
      }
    }
  }
}

void xDOSCAR::resetVDOS() { // CO20211124
  uint iatom = 0;
  uint iorbital = 0;
  uint ispin = 0;
  for (iatom = 0; iatom < vDOS.size(); iatom++) {
    for (iorbital = 0; iorbital < vDOS[iatom].size(); iorbital++) {
      for (ispin = 0; ispin < vDOS[iatom][iorbital].size(); ispin++) {
        std::fill(vDOS[iatom][iorbital][ispin].begin(), vDOS[iatom][iorbital][ispin].end(), 0.0);
      }
    }
  }
}

bool xDOSCAR::checkDOS(string& ERROR_out) const { // CO20191110
  const bool LDEBUG = (false || XHOST.DEBUG);
  stringstream message;

  const uint IENERGY = venergy.size();
  if (IENERGY == 0) {
    message << "GetProperties() failed. IENERGY==0 (venergy.size()==0)";
    ERROR_out = message.str();
    return false;
  }
  const uint IATOM = vDOS.size();
  if (IATOM == 0) {
    message << "GetProperties() failed. IATOM==0 (vDOS.size()==0)";
    ERROR_out = message.str();
    return false;
  }
  const uint IORBITAL = vDOS.front().size();
  if (IORBITAL == 0) {
    message << "GetProperties() failed. IORBITAL==0 (vDOS.front().size()==0)";
    ERROR_out = message.str();
    return false;
  }
  const uint ISPIN = vDOS.front().front().size();
  if (ISPIN == 0) {
    message << "GetProperties() failed. ISPIN==0 (vDOS.front().front().size()==0)";
    ERROR_out = message.str();
    return false;
  }
  for (uint iatom = 0; iatom < IATOM; iatom++) {
    if (vDOS[iatom].size() != IORBITAL) {
      message << "GetProperties() failed. (vDOS[iatom=" << iatom << "].size()!=vDOS.front().size())";
      ERROR_out = message.str();
      return false;
    }
    for (uint iorbital = 0; iorbital < IORBITAL; iorbital++) {
      if (vDOS[iatom][iorbital].size() != ISPIN) {
        message << "GetProperties() failed. (vDOS[iatom=" << iatom << "][iorbital=" << iorbital << "].size()!=vDOS.front().front().size())";
        ERROR_out = message.str();
        return false;
      }
      for (uint ispin = 0; ispin < ISPIN; ispin++) {
        if (vDOS[iatom][iorbital][ispin].size() != IENERGY) {
          message << "GetProperties() failed. (vDOS[iatom=" << iatom << "][iorbital=" << iorbital << "][ispin=" << ispin << "].size()!=venergy.size())";
          ERROR_out = message.str();
          return false;
        }
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " DOSCAR properties retrieved" << endl;
  }
  return true;
}

void xDOSCAR::GetVDOSSpecies(const xstructure& xstr) {
  return GetVDOSSpecies(xstr.num_each_type);
} //CO20191004
void xDOSCAR::GetVDOSSpecies(const deque<int>& num_each_type) { //CO20191004
  bool LDEBUG = (false || XHOST.DEBUG);
  stringstream message;

  if ((content.empty()) || (vcontent.empty())) {
    message << "xDOSCAR needs to be loaded before." << endl;
    message << "       GetProperties(const stringstream&);" << endl;
    message << "       GetProperties(const string&);" << endl;
    message << "       GetPropertiesFile(const string&);";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " DOSCAR content found" << endl;
  }

  string error_string; //keep this function const for plotter
  if (!checkDOS(error_string)) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, error_string, _INPUT_ERROR_);
  }; //quick check if GetProperties() failed

  //this should all work now that we checkDOS()
  uint IENERGY = venergy.size();
  uint IATOM = vDOS_atom.size();
  uint IORBITAL = vDOS_atom.front().size();
  uint IORBITAL_LM = vDOS_lm_atom.front().size();
  uint ISPIN = vDOS_atom.front().front().size();

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " IENERGY=" << IENERGY << endl;
    cerr << __AFLOW_FUNC__ << " IATOM=" << IATOM << endl;
    cerr << __AFLOW_FUNC__ << " IORBITAL=" << IORBITAL << endl;
    cerr << __AFLOW_FUNC__ << " IORBITAL_LM=" << IORBITAL_LM << endl;
    cerr << __AFLOW_FUNC__ << " ISPIN=" << ISPIN << endl;
  }

  //check that num_each_type corresponds to vDOS
  uint atoms_total = 0;
  for (uint iatom = 0; iatom < num_each_type.size(); iatom++) {
    atoms_total += num_each_type[iatom];
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " atoms_total=" << atoms_total << endl;
    cerr << __AFLOW_FUNC__ << " vDOS.size()=" << vDOS.size() << endl;
  }

  if (atoms_total + 1 != IATOM) { //total column
    message << "Input xstructure and DOS mismatch: atoms_total+1!=vDOS.size()" << endl;
    message << "For POCC runs, check that all ARUN.POCC's have the same num_each_type" << endl;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
  }

  //create vDOS_species with appropriate dimensions
  vDOS_species.assign(num_each_type.size() + 1, deque<deque<deque<double>>>(IORBITAL, deque<deque<double>>(ISPIN, deque<double>(number_energies, 0.0))));
  if (lmResolved) {
    vDOS_lm_species.assign(num_each_type.size() + 1, deque<deque<deque<double>>>(IORBITAL_LM, deque<deque<double>>(ISPIN, deque<double>(number_energies, 0.0))));
  }

  uint iatom = 0;
  double dos = 0; //ME trick
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " summing vDOS_species[ispecies] (individuals)" << endl;
  }
  for (uint ispecies = 0; ispecies < num_each_type.size(); ispecies++) {
    for (uint i = 0; i < (uint) num_each_type[ispecies]; i++) {
      for (uint iorbital = 1; iorbital < IORBITAL; iorbital++) { //iorbital==0 is total (VASP), skip that for now
        for (uint ispin = 0; ispin < ISPIN; ispin++) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " ispecies=" << ispecies << ", iorbital=" << iorbital << ", ispin=" << ispin << endl;
          }
          for (uint ienergy = 0; ienergy < IENERGY; ienergy++) {
            dos = vDOS_atom[iatom + 1][iorbital][ispin][ienergy]; //iatom==0 is total (VASP), skip that for now
            vDOS_species[ispecies + 1][iorbital][ispin][ienergy] += dos; //ispecies==0 is total (VASP), skip for now
            vDOS_species[ispecies + 1][0][ispin][ienergy] += dos;
            if (lmResolved && iorbital < IORBITAL_LM) {
              dos = vDOS_lm_atom[iatom + 1][iorbital][ispin][ienergy]; //iatom==0 is total (VASP), skip that for now
              vDOS_lm_species[ispecies + 1][iorbital][ispin][ienergy] += dos; //ispecies==0 is total (VASP), skip for now
              vDOS_lm_species[ispecies + 1][0][ispin][ienergy] += dos;
            }
          }
        }
      }
      iatom++;
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " getting vDOS_species[iorbital] (VASP totals)" << endl;
  }
  for (uint iorbital = 1; iorbital < IORBITAL; iorbital++) { //iorbital==0 is total (VASP), skip for now
    for (uint ispin = 0; ispin < ISPIN; ispin++) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " iorbital=" << iorbital << ", ispin=" << ispin << endl;
      }
      for (uint ienergy = 0; ienergy < IENERGY; ienergy++) {
        dos = vDOS_atom[0][iorbital][ispin][ienergy];
        vDOS_species[0][iorbital][ispin][ienergy] += dos;
        if (lmResolved && iorbital < IORBITAL_LM) {
          dos = vDOS_lm_atom[0][iorbital][ispin][ienergy];
          vDOS_lm_species[0][iorbital][ispin][ienergy] += dos;
        }
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " getting vDOS_species[0] (VASP totals)" << endl;
  }
  for (uint ispin = 0; ispin < ISPIN; ispin++) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ispin=" << ispin << endl;
    }
    for (uint ienergy = 0; ienergy < IENERGY; ienergy++) {
      dos = vDOS_atom[0][0][ispin][ienergy];
      vDOS_species[0][0][ispin][ienergy] += dos;
      if (lmResolved) {
        vDOS_lm_species[0][0][ispin][ienergy] += dos;
      }
    }
  }

  //[SD20230213 - OBSOLETE]return vDOS_species;
}

void xDOSCAR::GetVDOSIAtom(const xstructure& _xstr) { //SD20230227
  if (_xstr.iatoms_calculated) {
    return GetVDOSIAtom(_xstr.iatoms);
  }
  xstructure xstr = _xstr;
  xstr.sortAtomsEquivalent();
  return GetVDOSIAtom(xstr.iatoms);
}
void xDOSCAR::GetVDOSIAtom(const vector<vector<int>>& iatoms_index) { //SD20230227
  bool LDEBUG = (false || XHOST.DEBUG);
  stringstream message;

  if ((content.empty()) || (vcontent.empty())) {
    message << "xDOSCAR needs to be loaded before." << endl;
    message << "       GetProperties(const stringstream&);" << endl;
    message << "       GetProperties(const string&);" << endl;
    message << "       GetPropertiesFile(const string&);";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " DOSCAR content found" << endl;
  }

  string error_string; //keep this function const for plotter
  if (!checkDOS(error_string)) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, error_string, _INPUT_ERROR_);
  }; //quick check if GetProperties() failed

  //this should all work now that we checkDOS()
  uint IENERGY = venergy.size();
  uint IATOM = vDOS_atom.size();
  uint IORBITAL = vDOS_atom.front().size();
  uint IORBITAL_LM = vDOS_lm_atom.front().size();
  uint ISPIN = vDOS_atom.front().front().size();

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " IENERGY=" << IENERGY << endl;
    cerr << __AFLOW_FUNC__ << " IATOM=" << IATOM << endl;
    cerr << __AFLOW_FUNC__ << " IORBITAL=" << IORBITAL << endl;
    cerr << __AFLOW_FUNC__ << " IORBITAL_LM=" << IORBITAL_LM << endl;
    cerr << __AFLOW_FUNC__ << " ISPIN=" << ISPIN << endl;
  }
  //create vDOS_iatom with appropriate dimensions
  vDOS_iatom.assign(iatoms_index.size() + 1, deque<deque<deque<double>>>(IORBITAL, deque<deque<double>>(ISPIN, deque<double>(number_energies, 0.0))));
  if (lmResolved) {
    vDOS_lm_iatom.assign(iatoms_index.size() + 1, deque<deque<deque<double>>>(IORBITAL_LM, deque<deque<double>>(ISPIN, deque<double>(number_energies, 0.0))));
  }

  double dos = 0; //ME trick
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " summing vDOS_iatom[iatom] (individuals)" << endl;
  }
  for (uint iatom = 0; iatom < iatoms_index.size(); iatom++) {
    for (uint iorbital = 1; iorbital < IORBITAL; iorbital++) { //iorbital==0 is total (VASP), skip that for now
      for (uint ispin = 0; ispin < ISPIN; ispin++) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " iatom=" << iatom << ", iorbital=" << iorbital << ", ispin=" << ispin << endl;
        }
        for (uint ienergy = 0; ienergy < IENERGY; ienergy++) {
          dos = vDOS_atom[iatoms_index[iatom][0] + 1][iorbital][ispin][ienergy]; //iatom==0 is total (VASP), skip that for now
          vDOS_iatom[iatom + 1][iorbital][ispin][ienergy] += dos; //ispecies==0 is total (VASP), skip for now
          vDOS_iatom[iatom + 1][0][ispin][ienergy] += dos;
          if (lmResolved && iorbital < IORBITAL_LM) {
            dos = vDOS_lm_atom[iatoms_index[iatom][0] + 1][iorbital][ispin][ienergy]; //iatom==0 is total (VASP), skip that for now
            vDOS_lm_iatom[iatom + 1][iorbital][ispin][ienergy] += dos; //ispecies==0 is total (VASP), skip for now
            vDOS_lm_iatom[iatom + 1][0][ispin][ienergy] += dos;
          }
        }
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " getting vDOS_iatom[iorbital] (VASP totals)" << endl;
  }
  for (uint iorbital = 1; iorbital < IORBITAL; iorbital++) { //iorbital==0 is total (VASP), skip for now
    for (uint ispin = 0; ispin < ISPIN; ispin++) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " iorbital=" << iorbital << ", ispin=" << ispin << endl;
      }
      for (uint ienergy = 0; ienergy < IENERGY; ienergy++) {
        dos = vDOS_atom[0][iorbital][ispin][ienergy];
        vDOS_iatom[0][iorbital][ispin][ienergy] += dos;
        if (lmResolved && iorbital < IORBITAL_LM) {
          dos = vDOS_lm_atom[0][iorbital][ispin][ienergy];
          vDOS_lm_iatom[0][iorbital][ispin][ienergy] += dos;
        }
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " getting vDOS_iatom[0] (VASP totals)" << endl;
  }
  for (uint ispin = 0; ispin < ISPIN; ispin++) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ispin=" << ispin << endl;
    }
    for (uint ienergy = 0; ienergy < IENERGY; ienergy++) {
      dos = vDOS_atom[0][0][ispin][ienergy];
      vDOS_iatom[0][0][ispin][ienergy] += dos;
      if (lmResolved) {
        vDOS_lm_iatom[0][0][ispin][ienergy] += dos;
      }
    }
  }
}

bool xDOSCAR::GetBandGap(double EFERMI, double efermi_tol, double energy_tol, double occ_tol) { // CO20191004
  const bool LDEBUG = (false || XHOST.DEBUG);
  stringstream message;
  const bool force_exit = XHOST.POSTPROCESS; // SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only

  if ((content.empty()) || (vcontent.empty())) {
    message << "xDOSCAR needs to be loaded before." << endl;
    message << "GetProperties(const stringstream&);" << endl;
    message << "GetProperties(const string&);" << endl;
    message << "GetPropertiesFile(const string&);" << endl;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_MISSING_);
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " DOSCAR content found" << endl;
  }

  string error_string;
  if (!checkDOS(error_string)) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, error_string, _INPUT_ERROR_);
  }; // quick check if GetProperties() failed

  // this should all work now that we checkDOS()
  const uint IENERGY = venergy.size();
  const uint IATOM = vDOS.size();
  const uint IORBITAL = vDOS.front().size();
  const uint ISPIN = vDOS.front().front().size();

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " IENERGY=" << IENERGY << endl;
    cerr << __AFLOW_FUNC__ << " IATOM=" << IATOM << endl;
    cerr << __AFLOW_FUNC__ << " IORBITAL=" << IORBITAL << endl;
    cerr << __AFLOW_FUNC__ << " ISPIN=" << ISPIN << endl;
  }

  // GET FERMI LEVEL
  if (EFERMI == AUROSTD_NAN) {
    // VERY IMPORTANT - CO20171009 from discussion with SC
    // we strongly prefer to use Efermi from DOSCAR.static, not DOSCAR.bands
    // bands is not self-consistent (ICHARG=11), only used to determine energy states
    // since it is not self-consistent, it does not fill in states correctly
    // static is the best way to go!

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " Using Efermi from CURRENT DOSCAR (recommended to use Efermi from DOSCAR.static)" << endl;
    }
    EFERMI = Efermi;
  }

  // error in E-fermi, unless otherwise provided, is the difference between the input E-fermi (STATIC) and the one found in this DOSCAR (BANDS)
  // can be significant!
  // think of this as largest possible shift of energies in BANDS relative to real energy levels
  if (efermi_tol == AUROSTD_NAN) {
    efermi_tol = Efermi - EFERMI; // generlly, E-fermi BANDS > E-fermi STATIC
    if (std::signbit(efermi_tol)) {
      efermi_tol = energy_tol;
    } // just in case
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " tol(E-fermi)=" << efermi_tol << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " tol(occ)=" << occ_tol << endl;
  }

  valence_band_max.resize(ISPIN);
  conduction_band_min.resize(ISPIN);
  for (uint ispin = 0; ispin < ISPIN; ispin++) {
    valence_band_max[ispin] = AUROSTD_MAX_DOUBLE;
    conduction_band_min[ispin] = AUROSTD_MAX_DOUBLE;
  }

  bool found_EFERMI = false;
  bool metal_found = false;
  for (uint ispin = 0; ispin < ISPIN; ispin++) {
    found_EFERMI = false;
    for (size_t ienergy = 0; ienergy < vDOS.front().front()[ispin].size(); ienergy++) {
      if (venergy[ienergy] >= EFERMI) {
        if (valence_band_max[ispin] == AUROSTD_MAX_DOUBLE) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " DOS[energy=" << venergy[ienergy] << "]=" << vDOS.front().front()[ispin][ienergy] << " ?< " << occ_tol << endl;
          }
          if (found_EFERMI && (venergy[ienergy] - EFERMI) > efermi_tol) { // found metal for this spin
            if (LDEBUG) {
              cerr << __AFLOW_FUNC__ << " FOUND metal: ispin=" << ispin << endl;
            }
            metal_found = true;
            break;
          }
          if (vDOS.front().front()[ispin][ienergy] < occ_tol) {
            if (LDEBUG) {
              cerr << __AFLOW_FUNC__ << " FOUND valence_band_max[ispin=" << ispin << "]=" << venergy[ienergy] << endl;
            }
            valence_band_max[ispin] = venergy[ienergy];
          }
        } else {
          if (vDOS.front().front()[ispin][ienergy] >= occ_tol) {
            if (LDEBUG) {
              cerr << __AFLOW_FUNC__ << " FOUND conduction_band_min[ispin=" << ispin << "]=" << venergy[ienergy] << endl;
            }
            conduction_band_min[ispin] = venergy[ienergy];
            break;
          }
        }
        found_EFERMI = true;
      }
    }
    if (valence_band_max[ispin] != AUROSTD_MAX_DOUBLE && conduction_band_min[ispin] == AUROSTD_MAX_DOUBLE) {
      message << "Conduction band maximum not found [ispin=" << ispin << "]";
      if (force_exit) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      } else {
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
        return false;
      }
    }
  }

  if (!found_EFERMI) {
    message << " E-fermi not found";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
  }

  if (LDEBUG) {
    for (uint ispin = 0; ispin < ISPIN; ispin++) {
      cerr << __AFLOW_FUNC__ << " valence_band_max[ispin=" << ispin << "]=" << valence_band_max[ispin] << endl;
      cerr << __AFLOW_FUNC__ << " conduction_band_min[ispin=" << ispin << "]=" << conduction_band_min[ispin] << endl;
    }
    cerr << __AFLOW_FUNC__ << " metal_found=" << metal_found << endl;
  }

  //[CO20191004]double _METALGAP = -AUROSTD_NAN, _METALEDGE = -1.0;
  Egap.resize(ISPIN);
  Egap_fit.resize(ISPIN);
  Egap_type.resize(ISPIN);
  for (uint ispin = 0; ispin < ISPIN; ispin++) {
    if (valence_band_max[ispin] != AUROSTD_MAX_DOUBLE && conduction_band_min[ispin] != AUROSTD_MAX_DOUBLE) {
      if (valence_band_max[ispin] > conduction_band_min[ispin]) {
        message << "Negative band gap found in ispin=" << ispin << ", not expected";
        if (force_exit) {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
        } else {
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
          return false;
        }
      }
      Egap[ispin] = conduction_band_min[ispin] - valence_band_max[ispin];
      Egap_fit[ispin] = 1.348 * Egap[ispin] + 0.913;
      Egap_type[ispin] = "insulator";
      if (Egap[ispin] < energy_tol) {
        Egap_type[ispin] += "_zero_gap";
      }
    } else if (valence_band_max[ispin] == AUROSTD_MAX_DOUBLE && conduction_band_min[ispin] == AUROSTD_MAX_DOUBLE) {
      valence_band_max[ispin] = _METALEDGE_;
      conduction_band_min[ispin] = _METALEDGE_;
      Egap[ispin] = _METALGAP_;
      Egap_fit[ispin] = _METALGAP_;
      Egap_type[ispin] = "metal";
    } else {
      message << "Valence_band_max/conduction_band_min not found [ispin=" << ispin << "]";
      if (force_exit) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      } else {
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
        return false;
      }
    }
  }
  if (ISPIN == 1) {
    Egap_net = Egap[0];
    Egap_fit_net = Egap_fit[0];
    Egap_type_net = Egap_type[0];
  } else {
    if (metal_found) {
      Egap_net = _METALGAP_;
      Egap_fit_net = _METALGAP_;
      if (aurostd::substring2bool(Egap_type[0], "metal") && aurostd::substring2bool(Egap_type[1], "metal")) {
        Egap_type_net = "metal";
      } else {
        Egap_type_net = "half-metal";
      }
    } else {
      const double cbm = aurostd::min(conduction_band_min[0], conduction_band_min[1]);
      const double vbm = aurostd::max(valence_band_max[0], valence_band_max[1]);
      if (vbm > cbm) {
        message << "Negative band gap found in net, not expected";
        if (force_exit) {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
        } else {
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
          return false;
        }
      }
      Egap_net = cbm - vbm;
      Egap_fit_net = 1.348 * Egap_net + 0.913;
      Egap_type_net = "insulator";
      if (Egap_net < energy_tol) {
        Egap_type_net += "_zero_gap";
      }
    }
  }

  if (LDEBUG) {
    for (uint ispin = 0; ispin < ISPIN; ispin++) {
      cerr << __AFLOW_FUNC__ << " valence_band_max[ispin=" << ispin << "]=" << valence_band_max[ispin] << endl;
      cerr << __AFLOW_FUNC__ << " conduction_band_min[ispin=" << ispin << "]=" << conduction_band_min[ispin] << endl;
      cerr << __AFLOW_FUNC__ << " Egap[ispin=" << ispin << "]=" << Egap[ispin] << endl;
      cerr << __AFLOW_FUNC__ << " Egap_fit[ispin=" << ispin << "]=" << Egap_fit[ispin] << endl;
      cerr << __AFLOW_FUNC__ << " Egap_type[ispin=" << ispin << "]=" << Egap_type[ispin] << endl;
    }
    cerr << __AFLOW_FUNC__ << " Egap_net=" << Egap_net << endl;
    cerr << __AFLOW_FUNC__ << " Egap_fit_net=" << Egap_fit_net << endl;
    cerr << __AFLOW_FUNC__ << " Egap_type_net=" << Egap_type_net << endl;
  }

  return true;
}

deque<deque<deque<deque<double>>>> xDOSCAR::GetVDOSSpecies(const xstructure& xstr) const {
  return GetVDOSSpecies(xstr.num_each_type);
} // CO20191004
deque<deque<deque<deque<double>>>> xDOSCAR::GetVDOSSpecies(deque<int> num_each_type) const { // CO20191004
  const bool LDEBUG = (false || XHOST.DEBUG);
  stringstream message;

  if ((content.empty()) || (vcontent.empty())) {
    message << "xDOSCAR needs to be loaded before." << endl;
    message << "       GetProperties(const stringstream&);" << endl;
    message << "       GetProperties(const string&);" << endl;
    message << "       GetPropertiesFile(const string&);";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " DOSCAR content found" << endl;
  }

  string error_string; // keep this function const for plotter
  if (!checkDOS(error_string)) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, error_string, _INPUT_ERROR_);
  }; // quick check if GetProperties() failed

  // this should all work now that we checkDOS()
  const uint IENERGY = venergy.size();
  const uint IATOM = vDOS.size();
  const uint IORBITAL = vDOS.front().size();
  const uint ISPIN = vDOS.front().front().size();

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " IENERGY=" << IENERGY << endl;
    cerr << __AFLOW_FUNC__ << " IATOM=" << IATOM << endl;
    cerr << __AFLOW_FUNC__ << " IORBITAL=" << IORBITAL << endl;
    cerr << __AFLOW_FUNC__ << " ISPIN=" << ISPIN << endl;
  }

  // check that num_each_type corresponds to vDOS
  uint atoms_total = 0;
  for (size_t iatom = 0; iatom < num_each_type.size(); iatom++) {
    atoms_total += num_each_type[iatom];
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " atoms_total=" << atoms_total << endl;
    cerr << __AFLOW_FUNC__ << " vDOS.size()=" << vDOS.size() << endl;
  }

  if (atoms_total + 1 != IATOM) { // total column
    message << "Input xstructure and DOS mismatch: atoms_total+1!=vDOS.size()" << endl;
    message << "For POCC runs, check that all ARUN.POCC's have the same num_each_type" << endl;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
  }

  // create vDOS_species with appropriate dimensions
  deque<deque<deque<deque<double>>>> vDOS_species;
  vDOS_species.resize(num_each_type.size() + 1); //+1 for total
  for (size_t ispecies = 0; ispecies < vDOS_species.size(); ispecies++) {
    vDOS_species[ispecies].resize(IORBITAL);
    for (uint iorbital = 0; iorbital < IORBITAL; iorbital++) {
      vDOS_species[ispecies][iorbital].resize(ISPIN);
      for (uint ispin = 0; ispin < ISPIN; ispin++) {
        vDOS_species[ispecies][iorbital][ispin].assign(IENERGY, 0.0);
      }
    }
  }

  uint iatom = 0;
  double dos = 0; // ME trick
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " summing vDOS_species[ispecies] (individuals)" << endl;
  }
  for (size_t ispecies = 0; ispecies < num_each_type.size(); ispecies++) {
    for (uint i = 0; i < (uint) num_each_type[ispecies]; i++) {
      for (uint iorbital = 1; iorbital < IORBITAL; iorbital++) { // iorbital==0 is total (VASP), skip that for now
        for (uint ispin = 0; ispin < ISPIN; ispin++) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " ispecies=" << ispecies << ", iorbital=" << iorbital << ", ispin=" << ispin << endl;
          }
          for (uint ienergy = 0; ienergy < IENERGY; ienergy++) {
            dos = vDOS[iatom + 1][iorbital][ispin][ienergy]; // iatom==0 is total (VASP), skip that for now
            vDOS_species[ispecies + 1][iorbital][ispin][ienergy] += dos; // ispecies==0 is total (VASP), skip for now
            vDOS_species[ispecies + 1][0][ispin][ienergy] += dos;
          }
        }
      }
      iatom++;
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " getting vDOS_species[iorbital] (VASP totals)" << endl;
  }
  for (uint iorbital = 1; iorbital < IORBITAL; iorbital++) { // iorbital==0 is total (VASP), skip for now
    for (uint ispin = 0; ispin < ISPIN; ispin++) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " iorbital=" << iorbital << ", ispin=" << ispin << endl;
      }
      for (uint ienergy = 0; ienergy < IENERGY; ienergy++) {
        dos = vDOS[0][iorbital][ispin][ienergy];
        vDOS_species[0][iorbital][ispin][ienergy] += dos;
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " getting vDOS_species[0] (VASP totals)" << endl;
  }
  for (uint ispin = 0; ispin < ISPIN; ispin++) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ispin=" << ispin << endl;
    }
    for (uint ienergy = 0; ienergy < IENERGY; ienergy++) {
      dos = vDOS[0][0][ispin][ienergy];
      vDOS_species[0][0][ispin][ienergy] += dos;
    }
  }

  return vDOS_species;
}

// ME20190623 BEGIN
ostream& operator<<(ostream& oss, const xDOSCAR& xdos) {
  // Header
  oss << std::setw(4) << xdos.number_atoms << std::setw(4) << xdos.number_atoms << std::setw(4) << xdos.partial << std::setw(4) << 0 << std::endl;
  oss << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  oss << std::setprecision(7) << std::scientific;
  oss << std::setw(15) << xdos.Vol;
  for (int i = 1; i <= 3; i++) {
    oss << std::setw(15) << xdos.lattice[i];
  }
  oss << std::setw(15) << xdos.POTIM << std::endl;
  oss << "  " << std::setprecision(15) << xdos.temperature << std::endl;
  oss << "  " << xdos.carstring << std::endl;
  oss << " " << xdos.title << std::endl;

  stringstream dosline; // Will be reused for projected DOS
  dosline << std::dec << std::fixed << std::setprecision(8) << std::setw(15) << xdos.energy_max << std::fixed << std::setw(15) << xdos.energy_min << std::setprecision(0) << "  " << xdos.number_energies
          << std::setprecision(8) << std::fixed << std::setw(15) << xdos.Efermi << std::setprecision(8) << std::fixed << std::setw(15) << 1.0;
  oss << dosline.str() << std::endl;

  // Data
  oss << std::setprecision(4);
  for (uint e = 0; e < xdos.number_energies; e++) {
    oss << std::setw(12) << xdos.venergy[e];
    // Do not use vDOS.size() because it will mess up DOSCARs with spin-orbit coupling
    for (uint s = 0; s < xdos.spin + 1; s++) {
      oss << std::setw(12) << xdos.vDOS[0][0][s][e];
    }
    for (uint s = 0; s < xdos.spin + 1; s++) {
      oss << std::setw(12) << xdos.viDOS[s][e];
    }
    oss << std::endl;
  }

  for (size_t p = 1; p < xdos.vDOS.size(); p++) {
    oss << dosline.str() << std::endl;
    for (uint e = 0; e < xdos.number_energies; e++) {
      oss << std::setw(12) << xdos.venergy[e];
      for (size_t o = 1; o < xdos.vDOS[p].size(); o++) { // Do not output the total
        for (size_t s = 0; s < xdos.vDOS[p][o].size(); s++) {
          oss << std::setw(12) << xdos.vDOS[p][o][s][e];
        }
      }
      oss << std::endl;
    }
  }
  return oss;
}

// JSON serialization of xDOSCAR
aurostd::JSON::object xDOSCAR::serialize() const {
  return aurostd::JSON::object({AST_JSON_GETTER(JSON_xDOSCAR_MEMBERS)});
}

xDOSCAR xDOSCAR::deserialize(const aurostd::JSON::object& jo) {
  AST_JSON_SETTER(JSON_xDOSCAR_MEMBERS)
  return *this;
}

// ME20190623 END

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************

// ***************************************************************************
// class xEIGENVAL
bool xEIGENVAL::GetProperties(const string& stringIN, bool QUIET) {
  stringstream sss;
  sss.str(stringIN);
  if (filename.empty()) {
    filename = "string";
  }
  return xEIGENVAL::GetProperties(sss, QUIET);
}

bool xEIGENVAL::GetPropertiesFile(const string& fileIN, bool QUIET) {
  stringstream sss;
  if (filename.empty()) {
    filename = fileIN;
  }
  aurostd::compressfile2stringstream(fileIN, sss);
  return xEIGENVAL::GetProperties(sss, QUIET);
}

bool xEIGENVAL::GetPropertiesUrlFile(const string& url, const string& file, bool QUIET) {
  const string tmpfile = aurostd::TmpFileCreate("xEIGENVAL_GetProperties"); // CO20200502 - threadID
  aurostd::httpGetFileStatus(url + "/" + file, tmpfile);
  const bool out = GetPropertiesFile(tmpfile, QUIET);
  filename = "url=" + url; // CO20210315
  aurostd::RemoveFile(tmpfile);
  return out;
}

bool xEIGENVAL::GetProperties(const stringstream& stringstreamIN, bool QUIET) {
  const bool LDEBUG = (false || XHOST.DEBUG || !QUIET);
  stringstream message;
  const bool force_exit = XHOST.POSTPROCESS; // SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only

  const bool ERROR_flag = false;
  const long double seconds = aurostd::get_seconds();
  if (LDEBUG) {
    cerr << "xEIGENVAL::GetProperties: BEGIN (" << time_delay(seconds) << ")" << endl;
  }
  clear(); // so it does not mess up vector/deque
  content = stringstreamIN.str();
  vcontent.clear();
  const vector<string> vline;
  vector<string> tokens;
  // AS20200420 BEGIN]
  //  vasp 5.x.x uses '\n' as a separator between blocks of eigenvalues of
  //  different k-points, while vasp 4.x.x uses " \n" (space+newline).
  //  When a regular string2vectorsting is used, a sequence of "\n\n" is
  //  treated as a single '\n'. As a result separating line disappears
  //  in vcontent and the line indices are inconsistent.
  //  As a result, vweight, vkpoint and venergy vectors are not correctly populated.
  //  Setting consecutive = true in string2vectorstring fixes the issue.
  //
  //  aurostd::string2vectorstring(content,vcontent);
  aurostd::string2vectorstring(content, vcontent, true, false);
  // AS20200420 END]
  if (filename.empty()) {
    filename = "stringstream";
  }
  // crunching to eat the info

  // get parameters
  //  vline.clear();
  for (size_t iline = 0; iline < vcontent.size(); iline++) {
    aurostd::string2tokens(vcontent[iline], tokens);
    //    cerr << "iline=" << iline << "  " << vcontent[iline] << " tokens.size()=" << tokens.size() << endl;
    if (iline == 0 && tokens.size() == 4) {
      number_atoms = aurostd::string2utype<int>(tokens[0]); // ME20190623
      number_loops = aurostd::string2utype<int>(tokens[2]); // ME20190623
      spin = aurostd::string2utype<int>(tokens.at(tokens.size() - 1)) - 1;
    }
    if (iline == 1 && tokens.size() >= 5) {
      uint i = 0;
      Vol = aurostd::string2utype<double>(tokens.at(i++));
      lattice(1) = aurostd::string2utype<double>(tokens.at(i++));
      lattice(2) = aurostd::string2utype<double>(tokens.at(i++));
      lattice(3) = aurostd::string2utype<double>(tokens.at(i++));
      POTIM = aurostd::string2utype<double>(tokens.at(i++));
    }
    if (iline == 2) {
      temperature = aurostd::string2utype<double>(vcontent[iline]);
    }
    if (iline == 3) {
      carstring = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vcontent[iline]); // ME20190620 - what kind of EIGENVAL file?
    }
    if (iline == 4) {
      title = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vcontent[iline]); // ME20190614 - clean title
    }
    if (iline == 5 && tokens.size() >= 3) {
      uint i = 0;
      number_electrons = aurostd::string2utype<uint>(tokens.at(i++));
      number_kpoints = aurostd::string2utype<uint>(tokens.at(i++));
      number_bands = aurostd::string2utype<uint>(tokens.at(i++));
    }
    if (iline >= 7) {
      for (size_t jline = 0; jline < number_kpoints && iline < vcontent.size(); jline++, iline++) { // the iline++ is to get rid of the vacuum
        xvector<double> kpoint(3);
        aurostd::string2tokens(vcontent[iline], tokens);
        if (tokens.size() >= 4) {
          uint i = 0;
          kpoint(1) = aurostd::string2utype<double>(tokens.at(i++));
          kpoint(2) = aurostd::string2utype<double>(tokens.at(i++));
          kpoint(3) = aurostd::string2utype<double>(tokens.at(i++));
          vkpoint.push_back(kpoint);
          vweight.push_back(aurostd::string2utype<double>(tokens.at(i++)));
        }
        iline++; // move to energies
        deque<double> keigenval;
        deque<deque<double>> knenergy;
        double eigen; // ME20190614
        for (size_t kline = 0; kline < number_bands && iline < vcontent.size(); kline++, iline++) {
          aurostd::string2tokens(vcontent[iline], tokens);
          keigenval.clear();
          for (size_t i = 1; i < tokens.size(); i++) {
            // ME20190614 START
            eigen = aurostd::string2utype<double>(tokens[i]);
            if (eigen < energy_min) {
              energy_min = eigen;
            }
            if (eigen > energy_max) {
              energy_max = eigen;
            }
            keigenval.push_back(eigen);
            // ME20190614 END
          }
          knenergy.push_back(keigenval);
        }
        venergy.push_back(knenergy);
      }
    }
  }
  // ----------------------------------------------------------------------
  if (LDEBUG) {
    cerr << "xEIGENVAL::GetProperties: title=" << title << endl;
  }
  if (LDEBUG) {
    cerr << "xEIGENVAL::GetProperties: spin=" << spin << endl;
  }
  if (LDEBUG) {
    cerr << "xEIGENVAL::GetProperties: Vol=" << Vol << endl;
  }
  if (LDEBUG) {
    cerr << "xEIGENVAL::GetProperties: lattice=" << lattice << endl;
  }
  if (LDEBUG) {
    cerr << "xEIGENVAL::GetProperties: POTIM=" << POTIM << endl;
  }
  if (LDEBUG) {
    cerr << "xEIGENVAL::GetProperties: temperature=" << temperature << endl;
  }
  if (LDEBUG) {
    cerr << "xEIGENVAL::GetProperties: number_electrons=" << number_electrons << endl;
  }
  if (LDEBUG) {
    cerr << "xEIGENVAL::GetProperties: number_kpoints=" << number_kpoints << endl;
  }
  if (LDEBUG) {
    cerr << "xEIGENVAL::GetProperties: number_bands=" << number_bands << endl;
  }
  if (LDEBUG) {
    cerr << "xEIGENVAL::GetProperties: vweight.size()=" << vweight.size() << endl;
  }
  if (LDEBUG) {
    cerr << "xEIGENVAL::GetProperties: vkpoint.size()=" << vkpoint.size() << endl;
  }
  if (LDEBUG) {
    cerr << "xEIGENVAL::GetProperties: venergy.size()=" << venergy.size() << endl;
  }
  if (LDEBUG) {
    cerr << "xEIGENVAL::GetProperties: venergy.at(max).size()=" << venergy.at(venergy.size() - 1).size() << endl;
  }
  if (LDEBUG) {
    cerr << "xEIGENVAL::GetProperties: venergy.at(max).at(max).size()=" << venergy.at(venergy.size() - 1).at(venergy.at(venergy.size() - 1).size() - 1).size() << endl;
  }

  // for(size_t i=0;i<venergy.size();i++)
  //   if(LDEBUG) cerr << "xEIGENVAL::GetProperties: venergy.at.(" << i << ").size()=" << venergy.at(i).size() << endl;
  // for(size_t i=0;i<venergy.size();i++)
  //   for(size_t j=0;j<venergy.at(i).size();j++)
  //     if(LDEBUG) cerr << "xEIGENVAL::GetProperties: venergy.at.(" << i << ").at.(" << j << ").size()=" << venergy.at(i).at(j).size() << endl;

  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  if (LDEBUG) {
    cerr << "xEIGENVAL::GetProperties: END (" << time_delay(seconds) << ")" << endl;
  }
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  if (ERROR_flag) {
    message << "ERROR_flag set in xEIGENVAL";
    if (force_exit) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    } else {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_, QUIET);
    }
    return false;
  }
  return true;
}

ostream& operator<<(ostream& oss, const xEIGENVAL& xeigen) {
  // Header
  oss << std::setw(4) << xeigen.number_atoms << std::setw(4) << xeigen.number_atoms << std::setw(4) << xeigen.number_loops << std::setw(4) << (xeigen.spin + 1) << std::endl;
  oss << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  oss << std::setprecision(7) << std::scientific;
  oss << std::setw(15) << xeigen.Vol;
  for (int i = 1; i < 4; i++) {
    oss << std::setw(15) << xeigen.lattice[i];
  }
  oss << std::setw(15) << xeigen.POTIM << std::endl;
  oss << "  " << std::setprecision(15) << xeigen.temperature << std::endl;
  oss << "  " << xeigen.carstring << std::endl;
  oss << " " << xeigen.title << std::endl;
  oss << std::dec << std::setw(4) << xeigen.number_electrons << "  " << std::setw(4) << xeigen.number_kpoints << std::setw(4) << xeigen.number_bands << std::endl;

  // Data
  for (uint k = 0; k < xeigen.number_kpoints; k++) {
    oss << " " << std::endl; // space MUST be there or the EIGENVAL reader will fail
    oss << std::scientific << std::setprecision(8);
    for (int i = 1; i < 4; i++) {
      oss << "  " << xeigen.vkpoint[k][i];
    }
    oss << "  " << xeigen.vweight[k] << std::endl;
    for (uint br = 0; br < xeigen.number_bands; br++) {
      oss << std::dec << std::setprecision(0) << std::setw(4) << (br + 1);
      for (uint s = 0; s < xeigen.spin + 1; s++) {
        oss << std::setprecision(8) << std::fixed << std::setw(15) << xeigen.venergy[k][br][s];
      }
      oss << std::endl;
    }
  }
  return oss;
}

//---------------------------------------------------------------------------------
// GetEffectiveMass: Calculate the eff. mass. via the Harmonic approximation
bool GetEffectiveMass(xOUTCAR& xoutcar, xDOSCAR& xdoscar, xEIGENVAL& xeigenval, xstructure xstr, ostream& oss) {
  ofstream FileMESSAGE;
  return GetEffectiveMass(xoutcar, xdoscar, xeigenval, xstr, FileMESSAGE, oss);
} // CO20200404
bool GetEffectiveMass(xOUTCAR& xoutcar, xDOSCAR& xdoscar, xEIGENVAL& xeigenval, xstructure xstr, ofstream& FileMESSAGE, ostream& oss) { // CO20200404
  // The band gap, number of spins, and some other data are read out of the DOSCAR
  // The extrema (valleys) of the valence and conduction bands are gathered and
  // the curvature of the valley is solved for by fitting the E( kx, ky, kz) data to an ellipse.
  // From the curvature, the effective mass tensor is calculated, and the eigenvalues are found.
  // These eigenvalues are used as m_1, m_2, and m_3 for each valley.
  // The number of equivalent valleys in the Brillouin zone is found from the point group of the
  // reciprocal lattice.
  // The DoS effective masses are computed by averaging over the crystal using these data
  // in the following equation (copy & paste to LaTeX)
  // m^*_\text{carrier} = \left( \frac{\sum_i^\text{Nvalleys} M_i^2 m_{i1} m_{i2} m_{i3}}{ \text{Nvalleys} } \right)^{\frac{1}{3}}
  // and the conductivity effective masses are calculated by (copy & paste to LaTeX)
  // m^*_\text{carrier} = \frac{3 \left( \sum_i^{ \text{Nvalley}} M_i \right)}{\sum_i^{ \text{Nvalley}} M_i \left( \frac{1}{m_{i1}} + \frac{1}{m_{i2}} + \frac{1}{m_{i3}} \right)}
  // The resulting effective masses are written to the vectors
  // xOUTCAR.mass_elec_dos, xOUTCAR.mass_hole_dos, xOUTCAR.mass_elec_conduction, and xOUTCAR.mass_hole_conduction
  // algorithm depends on the energy in doscar.venergy being sorted in ascending order
  ///////////////////////////////////////////////////////////////////////
  stringstream message;

  xmatrix<double> reciprocal_lattice(1, 1, 3, 3);
  vector<vector<vector<kEn_st>>> fit_data_all;
  vector<vector<int>> number_of_valley_list;
  vector<double> band_info_vbt(2, xdoscar.Efermi);
  vector<double> band_info_cbb(2, xdoscar.Efermi);
  vector<int> valley_elec(2);
  vector<int> valley_hole(2);
  const double _ENER_RANGE = 0.026; //[CO20191004]_METALGAP =  -1.0E09;
  bool SPIN_UP = false;
  bool SPIN_DN = false;
  const string compound_name = xdoscar.title;
  const int ispin = xeigenval.spin + 1;

  // Eff. Mass arrays
  xoutcar.band_index.clear();
  xoutcar.carrier_type.clear();
  xoutcar.carrier_spin.clear();
  xoutcar.mass_elec_dos.clear();
  xoutcar.mass_elec_dos.resize(2);
  xoutcar.mass_hole_dos.clear();
  xoutcar.mass_hole_dos.resize(2);
  xoutcar.mass_elec_conduction.clear();
  xoutcar.mass_elec_conduction.resize(2);
  xoutcar.mass_hole_conduction.clear();
  xoutcar.mass_hole_conduction.resize(2);
  if (ispin == 1) {
    if (xoutcar.Egap.at(0) > _METALGAP_) {
      SPIN_UP = true;
      band_info_vbt.at(0) = xoutcar.valence_band_max.at(0) + xdoscar.Efermi;
      band_info_cbb.at(0) = xoutcar.conduction_band_min.at(0) + xdoscar.Efermi;
      band_info_vbt.at(1) = xoutcar.valence_band_max.at(0) + xdoscar.Efermi;
      band_info_cbb.at(1) = xoutcar.conduction_band_min.at(0) + xdoscar.Efermi;
    }
  } else if (ispin == 2) {
    if (xoutcar.Egap.at(0) > _METALGAP_) {
      SPIN_UP = true;
      band_info_vbt.at(0) = xoutcar.valence_band_max.at(0) + xdoscar.Efermi;
      band_info_cbb.at(0) = xoutcar.conduction_band_min.at(0) + xdoscar.Efermi;
    }
    if (xoutcar.Egap.at(1) > _METALGAP_) {
      SPIN_DN = true;
      band_info_vbt.at(1) = xoutcar.valence_band_max.at(1) + xdoscar.Efermi;
      band_info_cbb.at(1) = xoutcar.conduction_band_min.at(1) + xdoscar.Efermi;
    }
  } else {
    message << "EIGENVAL spin value is neither 0 nor 1"; // CO20200404
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_ERROR_); // CO20200404
    return false;
  }
  if (ispin == 1) {
    if (!SPIN_UP) {
      message << "Metallic system encountered"; // CO20200404
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_ERROR_); // CO20200404
      return false;
    }
  } else if (ispin == 2) {
    if (!SPIN_UP and !SPIN_DN) {
      message << "Metallic system encountered"; // CO20200404
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_ERROR_); // CO20200404
      return false;
    }
  }

  if (SPIN_UP or SPIN_DN) { // this disappears
    xstr.FixLattices();
    xstr.CalculateSymmetryPointGroupKLattice();
    reciprocal_lattice = xstr.klattice;
    vector<vector<kEn_st>> allkE_points;
    // allkE_points.at(eigenvalue).at(kpoint).{kEn_st details}
    allkE_points.resize(xeigenval.number_bands);
    xvector<double> temp_recip(1, 3);
    xvector<double> temp_cart(1, 3);
    // loop over KPOINTS (vkpoints)
    for (vector<int>::size_type ix = 0; ix != xeigenval.vkpoint.size(); ++ix) {
      kEn_st kp;
      temp_recip[1] = xeigenval.vkpoint.at(ix)[1];
      temp_recip[2] = xeigenval.vkpoint.at(ix)[2];
      temp_recip[3] = xeigenval.vkpoint.at(ix)[3];
      temp_cart = temp_recip * reciprocal_lattice;
      kp.kpoint(1) = temp_cart[1];
      kp.kpoint(2) = temp_cart[2];
      kp.kpoint(3) = temp_cart[3];
      // loop over EIGENVALUES (venergy)
      for (vector<int>::size_type iy = 0; iy != xeigenval.venergy.at(ix).size(); ++iy) {
        kp.band_index = iy + 1;
        // SPIN UNPOLARIZED
        if (ispin == 1) {
          kp.energy[0] = xeigenval.venergy.at(ix).at(iy).at(0);
          kp.energy[1] = kp.energy[0];
        }
        // SPIN POLARIZED
        else if (ispin == 2) {
          if (SPIN_UP and SPIN_DN) {
            kp.energy[0] = xeigenval.venergy.at(ix).at(iy).at(0);
            kp.energy[1] = xeigenval.venergy.at(ix).at(iy).at(1);
          } else if (SPIN_UP and !SPIN_DN) {
            kp.energy[0] = xeigenval.venergy.at(ix).at(iy).at(0);
            kp.energy[1] = kp.energy[0];
          } else if (!SPIN_UP and SPIN_DN) {
            kp.energy[0] = xeigenval.venergy.at(ix).at(iy).at(1);
            kp.energy[1] = kp.energy[0];
          }
          kp.energy[0] = xeigenval.venergy.at(ix).at(iy).at(0);
          kp.energy[1] = xeigenval.venergy.at(ix).at(iy).at(1);
        }
        allkE_points.at(kp.band_index - 1).push_back(kp);
      }
    }
    for (int spin_idx = 0; spin_idx < ispin; spin_idx++) {
      // CHOOSE BANDS CONTAINING EXTREMES
      vector<int> number_of_valley_list_tmp;
      vector<vector<kEn_st>> fit_data;
      // 1. sort the energy in ascending order
      for (size_t ii = 0; ii < allkE_points.size(); ii++) {
        if (spin_idx == 0) {
          sort(allkE_points[ii].begin(), allkE_points[ii].end(), comparison_kEn_str_up);
        } else if (spin_idx == 1) {
          sort(allkE_points[ii].begin(), allkE_points[ii].end(), comparison_kEn_str_dn);
        }
      }
      // 2. get the points with energy in the range
      for (size_t ii = 0; ii < allkE_points.size(); ii++) {
        // These 'if-then' statements determine the number of pockets in the system
        vector<kEn_st> fit_data_band;
        // CONDUCTION BANDS
        if (aurostd::abs(allkE_points[ii].front().energy[spin_idx] - band_info_cbb.at(spin_idx)) < _ENER_RANGE) {
          const double band_energy_minimum = allkE_points[ii].front().energy[spin_idx];
          for (int jj = 0; jj < _FIT_POINTS_NUMBER; jj++) {
            fit_data_band.push_back(allkE_points[ii].at(jj));
            fit_data_band.back().band_type = 1;
          }
          for (size_t jj = _FIT_POINTS_NUMBER; jj < allkE_points[ii].size(); jj++) {
            if (allkE_points[ii][jj].energy[spin_idx] - band_energy_minimum < _FIT_ENERGY_RANGE) {
              fit_data_band.push_back(allkE_points[ii][jj]);
              fit_data_band.back().band_type = 1;
            } else {
              break;
            }
          }
        }
        // VALENCE BANDS
        else if (aurostd::abs(band_info_vbt.at(spin_idx) - allkE_points[ii].back().energy[spin_idx]) < _ENER_RANGE) {
          const double band_energy_maximum = allkE_points[ii].back().energy[spin_idx];
          for (size_t jj = allkE_points[ii].size() - 1; jj > allkE_points[ii].size() - _FIT_POINTS_NUMBER - 1; jj--) {
            fit_data_band.push_back(allkE_points[ii].at(jj));
            fit_data_band.back().band_type = 0;
          }
          for (int jj = allkE_points[ii].size() - _FIT_POINTS_NUMBER - 1; jj >= 0; jj--) {
            if (band_energy_maximum - allkE_points[ii].at(jj).energy[spin_idx] < _FIT_ENERGY_RANGE) {
              fit_data_band.push_back(allkE_points[ii].at(jj));
              fit_data_band.back().band_type = 0;
            } else {
              break;
            }
          }
        }
        // fit_data contains the spin-specific pocket information
        if (!fit_data_band.empty()) {
          fit_data.push_back(fit_data_band);
        }
      }
      vector<double> max_distance;
      for (int i = 1; i < 4; i++) {
        double max_tmp;
        max_tmp = std::max(aurostd::abs(reciprocal_lattice[1][i]), aurostd::abs(reciprocal_lattice[2][i]));
        max_tmp = std::max(aurostd::abs(reciprocal_lattice[3][i]), max_tmp);
        max_tmp *= _BANDS_PARAMETER_MIN_RATIO;
        max_distance.push_back(max_tmp);
      }
      vector<vector<kEn_st>> fit_data_new;
      for (size_t i = 0; i < fit_data.size(); i++) {
        vector<kEn_st> fit_data_band;
        kEn_st kp;
        xvector<double> pt(1, 3);
        // the first point is closest to the extremes
        kp = fit_data[i].at(0);
        pt[1] = kp.kpoint(1);
        pt[2] = kp.kpoint(2);
        pt[3] = kp.kpoint(3);
        for (size_t ii = 0; ii < fit_data[i].size(); ii++) {
          xvector<double> pt1(1, 3);
          xvector<double> pt_sym(1, 3);
          kEn_st kp1;
          // the first point is closest to the extremes
          kp1 = fit_data[i][ii];
          pt1[1] = kp1.kpoint(1);
          pt1[2] = kp1.kpoint(2);
          pt1[3] = kp1.kpoint(3);
          for (size_t j = 0; j < xstr.pgroupk.size(); j++) {
            pt_sym = pt1 * xstr.pgroupk[j].Uc;
            if (near_to(pt, pt_sym, max_distance)) {
              // compare the distance between the most extreme points and the generated one
              kEn_st k1;
              k1.kpoint(1) = pt_sym[1];
              k1.kpoint(2) = pt_sym[2];
              k1.kpoint(3) = pt_sym[3];
              k1.energy[0] = kp1.energy[0];
              k1.energy[1] = kp1.energy[1];
              k1.band_index = kp1.band_index;
              k1.band_type = kp1.band_type;
              fit_data_band.push_back(k1);
            }
          }
          if (ii == 0) {
            const int number_of_valley = xstr.pgroupk.size() / fit_data_band.size();
            number_of_valley_list_tmp.push_back(number_of_valley);
          }
        }
        vector<kEn_st>::iterator it;
        sort(fit_data_band.begin(), fit_data_band.end(), comparison_kEn_str_position);
        it = unique(fit_data_band.begin(), fit_data_band.end(), is_equal_position_kEn_str);
        fit_data_band.resize(it - fit_data_band.begin());
        if (spin_idx == 1) {
          sort(fit_data_band.begin(), fit_data_band.end(), comparison_kEn_str_band_type_up);
        } else {
          sort(fit_data_band.begin(), fit_data_band.end(), comparison_kEn_str_band_type_dn);
        }
        fit_data_new.push_back(fit_data_band);
      }
      fit_data.clear();
      fit_data = fit_data_new;
      fit_data_all.push_back(fit_data);
      number_of_valley_list.push_back(number_of_valley_list_tmp);
    } // end spin_ndx
    /////////////////////////////////////////////////////////////////////////////
    // Spin Loop
    const vector<int> index_of_extremas_zero;
    int number_of_records = 0;
    for (int spin_idx = 0; spin_idx < ispin; spin_idx++) {
      // fit_data.size(): number of effective masses detected
      vector<vector<kEn_st>> fit_data = fit_data_all.at(spin_idx);
      vector<vector<double>> mass_eff_list;
      for (size_t ii = 0; ii < fit_data.size(); ii++) {
        const kEn_st kp1 = fit_data[ii].at(0);
        //    cout << "kp1 " << kp1.kpoint << endl;
        const int nrow = fit_data[ii].size() - 1;
        const int ncol = 9; // 9 polynomial coefficients, (a) thru (i)
        xvector<double> y_vec(1, nrow); // = En2 - En1
        xvector<double> y_sig(1, nrow);
        xmatrix<double> x_mat(1, 1, nrow, ncol);
        for (int jj = 1; jj < nrow + 1; jj++) {
          y_vec[jj] = fit_data[ii].at(jj).energy[spin_idx] - kp1.energy[spin_idx]; // = En2 - En1
          y_sig[jj] = _SIGMA; // _SIGMA = 1 (default std dev, see aflow.h)
        }
        // least-square fitting to an ellipses equation
        // En = a x^2 + b y^2 + c z^2 + d xy + e xz + f yz + g x + h y + i z + j
        // to get rid of j, we fit the function
        // En2 - En1 = a (x2^2  - x1^2)  +
        //             b (y2^2  - y1^2)  +
        //             c (z2^2  - z1^2)  +
        //             d (x2*y2 - x1*y1) +
        //             e (x2*z2 - x1*z1) +
        //             f (y2*z2 - y1*z1) +
        //             g (x2    - x1)    +
        //             h (y2    - y1)    +
        //             i (z2    - z1)
        for (int jj = 1; jj < nrow + 1; jj++) {
          const kEn_st kp2 = fit_data[ii].at(jj);
          //     cout << "kp2 " << kp2.kpoint << endl;
          x_mat[jj][1] = kp2.kpoint(1) * kp2.kpoint(1) - kp1.kpoint(1) * kp1.kpoint(1); // for(a)
          x_mat[jj][2] = kp2.kpoint(2) * kp2.kpoint(2) - kp1.kpoint(2) * kp1.kpoint(2); // for(b)
          x_mat[jj][3] = kp2.kpoint(3) * kp2.kpoint(3) - kp1.kpoint(3) * kp1.kpoint(3); // for(c)
          x_mat[jj][4] = kp2.kpoint(1) * kp2.kpoint(2) - kp1.kpoint(1) * kp1.kpoint(2); // for(d)
          x_mat[jj][5] = kp2.kpoint(1) * kp2.kpoint(3) - kp1.kpoint(1) * kp1.kpoint(3); // for(e)
          x_mat[jj][6] = kp2.kpoint(2) * kp2.kpoint(3) - kp1.kpoint(2) * kp1.kpoint(3); // for(f)
          x_mat[jj][7] = kp2.kpoint(1) - kp1.kpoint(1); // for(g)
          x_mat[jj][8] = kp2.kpoint(2) - kp1.kpoint(2); // for(h)
          x_mat[jj][9] = kp2.kpoint(3) - kp1.kpoint(3); // for(i)
        }
        // [x_mat] [a,b,c,d,e,f,g,h,i]=[y_vec]
        // check x_mat for columns of full of 0's - SINGULAR MATRIX
        // x_mat[rows][columns] <<-- same as Fortran
        //  cout << std::showpos;
        for (int rows = 1; rows <= nrow; rows++) {
          bool SINGULAR = true;
          for (int cols = 1; cols <= ncol; cols++) {
            //      cout << setw(10) << std::left << std::scientific << x_mat[rows][cols] << "  ";
            if (x_mat[rows][cols] != 0.0) {
              SINGULAR = false;
            }
          }
          //  cout << endl;
          if (SINGULAR) {
            message << "Singular system: ill-defined matrix problem encountered."; // CO20200404
            pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_WARNING_); // CO20200404
            return false;
          }
        }
        //    cout << "=====================================" << endl;
        aurostd::cematrix ECI_matrix(x_mat);
        ECI_matrix.LeastSquare(y_vec, y_sig); // minimization happens here
        // starting @ 1, ending at 3 and 3x3
        xmatrix<double> mass_m(1, 1, 3, 3);
        // looks like upper triangular stuff happens here
        mass_m[1][1] = ECI_matrix.AVec().at(0);
        mass_m[2][2] = ECI_matrix.AVec().at(1);
        mass_m[3][3] = ECI_matrix.AVec().at(2);
        mass_m[1][2] = ECI_matrix.AVec().at(3) * 0.5;
        mass_m[1][3] = ECI_matrix.AVec().at(4) * 0.5;
        mass_m[2][3] = ECI_matrix.AVec().at(5) * 0.5;
        mass_m[2][1] = mass_m[1][2];
        mass_m[3][1] = mass_m[1][3];
        mass_m[3][2] = mass_m[2][3];
        aurostd::cematrix mass_m_ce(mass_m);
        xvector<double> mr = mass_m_ce.EigenValues();
        vector<double> mr_tmp;
        for (int jj = 1; jj <= 3; jj++) {
          mr[jj] = 1.0 * _MASS_FACTOR / mr[jj];
          mr_tmp.push_back(mr[jj]);
        }
        sort(mr_tmp.begin(), mr_tmp.end());
        mass_eff_list.push_back(mr_tmp);
      }

      vector<double> elec_cond_mass_per_valley;
      vector<double> hole_cond_mass_per_valley;
      vector<int> elec_valley_multiplicity;
      vector<int> hole_valley_multiplicity;
      xoutcar.mass_elec_dos.at(spin_idx) = 0.0;
      xoutcar.mass_hole_dos.at(spin_idx) = 0.0;
      valley_elec.at(spin_idx) = 0;
      valley_hole.at(spin_idx) = 0;

      // problem starts here
      for (size_t ii = 0; ii < mass_eff_list.size(); ii++) {
        const vector<double> tempvector(3);
        if (fit_data.at(ii).at(0).band_type == 0) {
          xoutcar.carrier_type.emplace_back("hole");
        } else if (fit_data.at(ii).at(0).band_type == 1) {
          xoutcar.carrier_type.emplace_back("elec");
        }
        xoutcar.band_index.push_back(fit_data.at(ii).at(0).band_index);
        xoutcar.carrier_spin.push_back(spin_idx);
        xoutcar.extrema_cart_coord.push_back(tempvector);
        xoutcar.effective_mass_axes.push_back(tempvector);
        number_of_records++;
        double dos_mass = 1.0;
        double cond_mass = 0.0;
        // diagonalized eff mass tensor
        for (uint jj = 0; jj < 3; jj++) {
          dos_mass *= mass_eff_list[ii].at(jj);
          const double one = 1.0;
          // add the reciprocal of 3 individual masses
          cond_mass += one / mass_eff_list[ii].at(jj);
          // the following causes out of bounds errors
          xoutcar.extrema_cart_coord.back().at(jj) = fit_data_all.at(0).at(ii).at(0).kpoint(jj + 1);
          xoutcar.effective_mass_axes.back().at(jj) = mass_eff_list[ii].at(jj);
        }
        const int number_of_valley = number_of_valley_list.at(spin_idx).at(ii);
        xoutcar.equivalent_valley.push_back(number_of_valley_list.at(spin_idx).at(ii));
        dos_mass *= number_of_valley * number_of_valley;
        if (xoutcar.carrier_type.back() == "elec") {
          xoutcar.mass_elec_dos.at(spin_idx) += dos_mass;
          valley_elec.at(spin_idx)++;
        } else if (xoutcar.carrier_type.back() == "hole") {
          xoutcar.mass_hole_dos.at(spin_idx) += dos_mass;
          valley_hole.at(spin_idx)++;
        }
        dos_mass = std::pow(dos_mass, 1.0 / 3.0);
        double temp = 1.0;
        temp /= cond_mass;
        // flip sum of inverses to denominator
        cond_mass = temp;
        // store cond_mass for the valley and store multiplicity of the valley
        if (xoutcar.carrier_type.back() == "elec") {
          elec_cond_mass_per_valley.push_back(cond_mass);
          elec_valley_multiplicity.push_back(number_of_valley);
        } else if (xoutcar.carrier_type.back() == "hole") {
          hole_cond_mass_per_valley.push_back(cond_mass);
          hole_valley_multiplicity.push_back(number_of_valley);
        }
        // final factor of 3 for single valley total
        cond_mass *= 3.0;
        xoutcar.effective_mass_DOS.push_back(dos_mass);
        xoutcar.effective_mass_COND.push_back(cond_mass);
      } // END: loop over mass_eff_list
      // problem ends here

    }
    // end Spin Loop

  }

  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  return true;
}
/////////////////////////////////////////////////////////////////////////////////////
// spin dependent comparisons
bool comparison_kEn_str_up(const kEn_st& k1, const kEn_st& k2) {
  return static_cast<bool>(k1.energy[0] < k2.energy[0]);
}
bool comparison_kEn_str_dn(const kEn_st& k1, const kEn_st& k2) {
  return static_cast<bool>(k1.energy[1] < k2.energy[1]);
}
bool comparison_kEn_str_band_type_up(const kEn_st& k1, const kEn_st& k2) {
  if (k1.band_type == 1 && k2.band_type == 1) { // both are in conduction bands
    return static_cast<bool>(k1.energy[0] < k2.energy[0]);
  } else {
    if (k1.band_type == 0 && k2.band_type == 0) { // both are in valence bands
      return static_cast<bool>(k1.energy[0] > k2.energy[0]);
    } else {
      return static_cast<bool>(k1.energy[0] < k2.energy[0]);
    }
  }
}
bool comparison_kEn_str_band_type_dn(const kEn_st& k1, const kEn_st& k2) {
  if (k1.band_type == 1 && k2.band_type == 1) {
    return static_cast<bool>(k1.energy[1] < k2.energy[1]);
  } else {
    if (k1.band_type == 0 && k2.band_type == 0) {
      return static_cast<bool>(k1.energy[1] > k2.energy[1]);
    } else {
      return static_cast<bool>(k1.energy[1] < k2.energy[1]);
    }
  }
}
/////////////////////////////////////////////////////////////////////////////////////
bool comparison_kEn_str_position(const kEn_st& k1, const kEn_st& k2) {
  const double _ENER_EPS = 1.0e-4;
  bool flag = false;
  if (aurostd::abs(k1.kpoint(1) - k2.kpoint(1)) < _ENER_EPS) {
    if (aurostd::abs(k1.kpoint(2) - k2.kpoint(2)) < _ENER_EPS) {
      if (k1.kpoint(3) < k2.kpoint(3)) {
        flag = true;
      } else {
        flag = false;
      }
    } else {
      if (k1.kpoint(2) < k2.kpoint(2)) {
        flag = true;
      } else {
        flag = false;
      }
    }
  } else {
    if (k1.kpoint(1) < k2.kpoint(1)) {
      flag = true;
    } else {
      flag = false;
    }
  }
  return flag;
}
bool is_equal_position_kEn_str(const kEn_st& k1, const kEn_st& k2) {
  const double _ENER_EPS = 1.0e-4;
  return isequal(k1.kpoint, k2.kpoint, double(_ENER_EPS));
}
bool near_to(const xvector<double>& k1, const xvector<double>& k2, const vector<double>& max_distance) {
  return (aurostd::abs(k1[1] - k2[1]) < max_distance.at(0) && aurostd::abs(k1[2] - k2[2]) < max_distance.at(1) && aurostd::abs(k1[3] - k2[3]) < max_distance.at(2));
}
//-------------------------------------------------------------------------------------------------
// PrintBandGap: Print the output of xOUTCAR::GetBandGap
// 2014: Camilo E. Calderon
// 2017: Corey Oses
bool PrintBandGap(string& directory, ostream& oss) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  stringstream message; // CO20200404
  stringstream ss_outcar_static("");
  stringstream ss_outcar_bands("");
  string path_outcar_static;
  string path_outcar_bands;
  string path_POSCAR;
  xOUTCAR xoutcar_static;
  xOUTCAR xoutcar_bands;
  const char LastChar = *directory.rbegin();
  if (LastChar == '/') {
    directory.erase(directory.size() - 1);
  }

  const long double seconds = aurostd::get_seconds();
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BEGIN (" << time_delay(seconds) << ")" << endl;
  }

  // OUTCAR_bands
  if (ss_outcar_bands.str().empty() && aurostd::CompressFileExist(directory + "/OUTCAR.bands", path_outcar_bands)) {
    aurostd::compressfile2stringstream(path_outcar_bands, ss_outcar_bands); // CO20200404
  }
  if (ss_outcar_bands.str().empty()) {
    message << "OUTCAR.bands not found"; // CO20200404
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, oss, _LOGGER_WARNING_); // CO20200404
  } // CO20200404
  if (ss_outcar_bands.str().empty() && aurostd::CompressFileExist(directory + "/OUTCAR", path_outcar_bands)) {
    aurostd::compressfile2stringstream(path_outcar_bands, ss_outcar_bands); // CO20200404
  }
  if (ss_outcar_bands.str().empty()) { // CO20200404
    message << "OUTCAR not found"; // CO20200404
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, oss, _LOGGER_ERROR_); // CO20200404
    return false;
  }

  // CO20171002 - using tolerance from symmetry calc - START
  // double tol;
  // if(xstr.CalculateSymmetry()){tol=xstr.sym_eps;}
  // else {tol=SYM::defaultTolerance(xstr);}
  // cerr << tol << endl;
  // CO20171002 - using tolerance from symmetry calc - STOP

  if (!xoutcar_bands.GetProperties(ss_outcar_bands)) {
    return false;
  }
  // try to grab xstr from OUTCAR
  if (!xoutcar_bands.GetXStructure()) {
    if (!aurostd::FileExist(directory + "/POSCAR.bands", path_POSCAR) && !aurostd::CompressFileExist(directory + "/POSCAR.bands", path_POSCAR)) {
      return false;
    }
    const xstructure xstr(path_POSCAR, IOVASP_POSCAR);
    xoutcar_bands.xstr = xstr;
  }

  double EFERMI = xoutcar_bands.Efermi; // hopefully we can grab this from static, otherwise settle on the one in bands

  // OUTCAR_static - try to grab right Efermi
  if (!aurostd::CompressFileExist(directory + "/OUTCAR.static", path_outcar_static)) {
    message << "OUTCAR.static not found, defaulting E-Fermi to that found in OUTCAR.bands"; // CO20200404
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, oss, _LOGGER_WARNING_); // CO20200404
    // return false;
  } else {
    aurostd::compressfile2stringstream(path_outcar_static, ss_outcar_static); // CO20200404
    if (!xoutcar_static.GetProperties(ss_outcar_static)) {
      message << "Defaulting E-Fermi to that found in OUTCAR.bands"; // CO20200404
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, oss, _LOGGER_WARNING_); // CO20200404
      // return false;
    } else {
      EFERMI = xoutcar_static.Efermi;
      if (LDEBUG) {
        cerr << XPID << "xOUTCAR::PrintBandGap: Found E-fermi from OUTCAR.static: " << EFERMI << endl;
      }
    }
  }

  if (!xoutcar_bands.GetBandGap(EFERMI)) {
    return false;
  }

  // could use pointers to save memory, but these objects are so small it doesn't matter
  // make full copies and forget
  const string SYSTEM = xoutcar_bands.SYSTEM;
  vector<double> valence_band_max = xoutcar_bands.valence_band_max;
  vector<double> conduction_band_min = xoutcar_bands.conduction_band_min;
  vector<double> Egap = xoutcar_bands.Egap;
  vector<double> Egap_fit = xoutcar_bands.Egap_fit;
  vector<string> Egap_type = xoutcar_bands.Egap_type;
  const double valence_band_max_net = xoutcar_bands.valence_band_max_net;
  const double conduction_band_min_net = xoutcar_bands.conduction_band_min_net;
  const double Egap_net = xoutcar_bands.Egap_net;
  const double Egap_fit_net = xoutcar_bands.Egap_fit_net;
  const string Egap_type_net = xoutcar_bands.Egap_type_net;

  // SUCCESS!
  if (Egap.size() == 1) {
    oss.precision(4);
    oss << "System        :   " << SYSTEM << endl;
    oss << "Spin tag      :   " << Egap.size() << endl;
    oss << "Fermi level   :  " << std::scientific << std::showpos << EFERMI << endl;
    oss << "                  VBT           CBB           Egap          Egap_fit     Type" << endl;
    oss << "Net Result    :  ";
    oss << setw(14) << std::left << std::scientific << std::showpos << valence_band_max.at(0);
    oss << setw(14) << std::left << std::scientific << std::showpos << conduction_band_min.at(0);
    oss << setw(14) << std::left << std::scientific << std::showpos << Egap.at(0);
    oss << setw(14) << std::left << std::scientific << std::showpos << Egap_fit.at(0);
    oss << setw(20) << std::left << Egap_type.at(0) << endl;
    oss << endl;
  } else if (Egap.size() == 2) {
    oss.precision(4);
    oss << "System        :   " << SYSTEM << endl;
    oss << "Spin tag      :   " << Egap.size() << endl;
    oss << "Fermi level   :  " << std::scientific << std::showpos << EFERMI << endl;
    oss << "                  VBT           CBB           Egap          Egap_fit     Type" << endl;
    oss << "Majority Spin :  ";
    oss << setw(14) << std::left << std::scientific << std::showpos << valence_band_max.at(0);
    oss << setw(14) << std::left << std::scientific << std::showpos << conduction_band_min.at(0);
    oss << setw(14) << std::left << std::scientific << std::showpos << Egap.at(0);
    oss << setw(14) << std::left << std::scientific << std::showpos << Egap_fit.at(0);
    oss << setw(20) << std::left << Egap_type.at(0) << endl;
    oss << "Minority Spin :  ";
    oss << setw(14) << std::left << std::scientific << std::showpos << valence_band_max.at(1);
    oss << setw(14) << std::left << std::scientific << std::showpos << conduction_band_min.at(1);
    oss << setw(14) << std::left << std::scientific << std::showpos << Egap.at(1);
    oss << setw(14) << std::left << std::scientific << std::showpos << Egap_fit.at(1);
    oss << setw(20) << std::left << Egap_type.at(1) << endl;
    oss << "Net Result    :  ";
    oss << setw(14) << std::left << std::scientific << std::showpos << valence_band_max_net;
    oss << setw(14) << std::left << std::scientific << std::showpos << conduction_band_min_net;
    oss << setw(14) << std::left << std::scientific << std::showpos << Egap_net;
    oss << setw(14) << std::left << std::scientific << std::showpos << Egap_fit_net;
    oss << setw(20) << std::left << Egap_type_net << endl;
    oss << endl;
  }
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " END (" << time_delay(seconds) << ")" << endl;
  }
  return true;
}
//-------------------------------------------------------------------------------------------------
// PrintBandGap: Print the output of xDOSCAR::GetBandGap
// 2019: Corey Oses
bool PrintBandGap_DOS(string& directory, ostream& oss) { // CO20191110
  const bool LDEBUG = (false || XHOST.DEBUG);
  stringstream message; // CO20200404
  const char LastChar = *directory.rbegin();
  if (LastChar == '/') {
    directory.erase(directory.size() - 1);
  }

  const long double seconds = aurostd::get_seconds();
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BEGIN (" << time_delay(seconds) << ")" << endl;
  }

  string path_doscar_static;
  stringstream ss_doscar_static;
  if (ss_doscar_static.str().empty() && aurostd::FileExist(directory + "/DOSCAR.static", path_doscar_static)) {
    aurostd::file2stringstream(path_doscar_static, ss_doscar_static); // CO20200404
  }
  if (ss_doscar_static.str().empty() && aurostd::CompressFileExist(directory + "/DOSCAR.static", path_doscar_static)) {
    aurostd::compressfile2stringstream(path_doscar_static, ss_doscar_static); // CO20200404
  }
  if (ss_doscar_static.str().empty()) {
    message << "DOSCAR.static not found"; // CO20200404
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, oss, _LOGGER_WARNING_); // CO20200404
  } // CO20200404
  if (ss_doscar_static.str().empty() && aurostd::FileExist(directory + "/DOSCAR", path_doscar_static)) {
    aurostd::file2stringstream(path_doscar_static, ss_doscar_static); // CO20200404
  }
  if (ss_doscar_static.str().empty() && aurostd::CompressFileExist(directory + "/DOSCAR", path_doscar_static)) {
    aurostd::compressfile2stringstream(path_doscar_static, ss_doscar_static); // CO20200404
  }
  if (ss_doscar_static.str().empty()) { // CO20200404
    message << "DOSCAR not found"; // CO20200404
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, oss, _LOGGER_ERROR_); // CO20200404
    return false;
  }

  xDOSCAR xdoscar;
  if (!xdoscar.GetProperties(ss_doscar_static)) {
    return false;
  }
  if (!xdoscar.GetBandGap()) { // use EFERMI of doscar.static (most accurate)
    return false;
  }

  // could use pointers to save memory, but these objects are so small it doesn't matter
  // make full copies and forget
  const double EFERMI = xdoscar.Efermi;
  const string SYSTEM = xdoscar.title;
  vector<double> valence_band_max = xdoscar.valence_band_max;
  vector<double> conduction_band_min = xdoscar.conduction_band_min;
  vector<double> Egap = xdoscar.Egap;
  vector<double> Egap_fit = xdoscar.Egap_fit;
  vector<string> Egap_type = xdoscar.Egap_type;
  const double valence_band_max_net = xdoscar.valence_band_max_net;
  const double conduction_band_min_net = xdoscar.conduction_band_min_net;
  const double Egap_net = xdoscar.Egap_net;
  const double Egap_fit_net = xdoscar.Egap_fit_net;
  const string Egap_type_net = xdoscar.Egap_type_net;

  // SUCCESS!
  if (Egap.size() == 1) {
    oss.precision(4);
    oss << "System        :   " << SYSTEM << endl;
    oss << "Spin tag      :   " << Egap.size() << endl;
    oss << "Fermi level   :  " << std::scientific << std::showpos << EFERMI << endl;
    oss << "                  VBT           CBB           Egap          Egap_fit     Type" << endl;
    oss << "Net Result    :  ";
    oss << setw(14) << std::left << std::scientific << std::showpos << valence_band_max.at(0);
    oss << setw(14) << std::left << std::scientific << std::showpos << conduction_band_min.at(0);
    oss << setw(14) << std::left << std::scientific << std::showpos << Egap.at(0);
    oss << setw(14) << std::left << std::scientific << std::showpos << Egap_fit.at(0);
    oss << setw(20) << std::left << Egap_type.at(0) << endl;
    oss << endl;
  } else if (Egap.size() == 2) {
    oss.precision(4);
    oss << "System        :   " << SYSTEM << endl;
    oss << "Spin tag      :   " << Egap.size() << endl;
    oss << "Fermi level   :  " << std::scientific << std::showpos << EFERMI << endl;
    oss << "                  VBT           CBB           Egap          Egap_fit     Type" << endl;
    oss << "Majority Spin :  ";
    oss << setw(14) << std::left << std::scientific << std::showpos << valence_band_max.at(0);
    oss << setw(14) << std::left << std::scientific << std::showpos << conduction_band_min.at(0);
    oss << setw(14) << std::left << std::scientific << std::showpos << Egap.at(0);
    oss << setw(14) << std::left << std::scientific << std::showpos << Egap_fit.at(0);
    oss << setw(20) << std::left << Egap_type.at(0) << endl;
    oss << "Minority Spin :  ";
    oss << setw(14) << std::left << std::scientific << std::showpos << valence_band_max.at(1);
    oss << setw(14) << std::left << std::scientific << std::showpos << conduction_band_min.at(1);
    oss << setw(14) << std::left << std::scientific << std::showpos << Egap.at(1);
    oss << setw(14) << std::left << std::scientific << std::showpos << Egap_fit.at(1);
    oss << setw(20) << std::left << Egap_type.at(1) << endl;
    oss << "Net Result    :  ";
    oss << setw(14) << std::left << std::scientific << std::showpos << valence_band_max_net;
    oss << setw(14) << std::left << std::scientific << std::showpos << conduction_band_min_net;
    oss << setw(14) << std::left << std::scientific << std::showpos << Egap_net;
    oss << setw(14) << std::left << std::scientific << std::showpos << Egap_fit_net;
    oss << setw(20) << std::left << Egap_type_net << endl;
    oss << endl;
  }
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " END (" << time_delay(seconds) << ")" << endl;
  }
  return true;
}
//---------------------------------------------------------------------------------
// PrintEffectiveMass: Print the output of xOUTCAR::GetEffectiveMass
// 2014: Camilo E. Calderon
bool PrintEffectiveMass(string& directory, ostream& oss) {
  //[CO20191004]double _METALGAP = -1.0E09;
  bool SPIN_UP = false;
  bool SPIN_DN = false;
  bool EM_TAG;
  stringstream file_OUTCAR;
  stringstream file_DOSCAR;
  stringstream file_EIGENVAL;
  string path_POSCAR;

  const vector<string> vext{".bz2", ".xz", ".gz"};

  xOUTCAR xoutcar;
  xDOSCAR xdoscar;
  xEIGENVAL xeigenval;
  const char LastChar = *directory.rbegin();
  if (LastChar == '/') {
    directory.erase(directory.size() - 1);
  }

  // OUTCAR
  bool found_OUTCAR = false;
  for (size_t iext = 0; iext < vext.size(); iext++) {
    if (!found_OUTCAR && aurostd::FileExist(directory + "/OUTCAR.static" + vext[iext])) {
      found_OUTCAR = true;
      aurostd::compressfile2stringstream(directory + "/OUTCAR.static" + vext[iext], file_OUTCAR); // .EXT
    }
  }
  if (!found_OUTCAR && aurostd::FileExist(directory + "/OUTCAR.static")) {
    found_OUTCAR = true;
    aurostd::file2stringstream(directory + "/OUTCAR.static", file_OUTCAR); // plain text
  }
  if (!found_OUTCAR) {
    return false;
  }

  // DOSCAR
  bool found_DOSCAR = false;
  for (size_t iext = 0; iext < vext.size(); iext++) {
    if (!found_DOSCAR && aurostd::FileExist(directory + "/DOSCAR.static" + vext[iext])) {
      found_DOSCAR = true;
      aurostd::compressfile2stringstream(directory + "/DOSCAR.static" + vext[iext], file_DOSCAR); // .EXT
    }
  }
  if (!found_DOSCAR && aurostd::FileExist(directory + "/DOSCAR.static")) {
    found_DOSCAR = true;
    aurostd::file2stringstream(directory + "/DOSCAR.static", file_DOSCAR); // plain text
  }
  if (!found_DOSCAR) {
    return false;
  }

  // EIGENVAL
  bool found_EIGENVAL = false;
  for (size_t iext = 0; iext < vext.size(); iext++) {
    if (!found_EIGENVAL && aurostd::FileExist(directory + "/EIGENVAL.static" + vext[iext])) {
      found_EIGENVAL = true;
      aurostd::compressfile2stringstream(directory + "/EIGENVAL.static" + vext[iext], file_EIGENVAL); // .EXT
    }
  }
  if (!found_EIGENVAL && aurostd::FileExist(directory + "/EIGENVAL.static")) {
    found_EIGENVAL = true;
    aurostd::file2stringstream(directory + "/EIGENVAL.static", file_EIGENVAL); // plain text
  }
  if (!found_EIGENVAL) {
    return false;
  }

  // POSCAR
  bool found_POSCAR = false;
  for (size_t iext = 0; iext < vext.size(); iext++) {
    if (!found_POSCAR && aurostd::FileExist(directory + "/POSCAR.bands" + vext[iext])) {
      found_POSCAR = true;
      path_POSCAR = directory + "/POSCAR.bands" + vext[iext]; // EXR
    }
  }
  if (!found_POSCAR && aurostd::FileExist(directory + "/POSCAR.bands")) {
    found_POSCAR = true;
    path_POSCAR = directory + "/POSCAR.bands"; // plain text
  }
  if (!found_POSCAR) {
    return false;
  }

  // GET THE DATA
  xoutcar.GetProperties(file_OUTCAR);
  xstructure xstr(path_POSCAR, IOVASP_POSCAR);

  // CO20171002 - using tolerance from symmetry calc - START
  double tol;
  if (xstr.CalculateSymmetry()) {
    tol = xstr.sym_eps;
  } else {
    tol = SYM::defaultTolerance(xstr);
  }
  // CO20171002 - using tolerance from symmetry calc - STOP

  xoutcar.GetBandGap(tol);
  xdoscar.GetProperties(file_DOSCAR);
  xeigenval.GetProperties(file_EIGENVAL);
  // EFFECTIVE MASSES
  EM_TAG = GetEffectiveMass(xoutcar, xdoscar, xeigenval, xstr);
  if (!EM_TAG) {
    oss << "System              : " << xoutcar.SYSTEM << endl;
    oss << string("FAILED") << "   filename[=" << xoutcar.filename << "]" << endl;
    oss << endl;
    return false;
  }
  // METALLIC CHANNELS
  const int ispin = xeigenval.spin + 1;
  if (ispin == 1) {
    if (xoutcar.Egap.at(0) > _METALGAP_) {
      SPIN_UP = true;
    }
  } else if (ispin == 2) {
    if (xoutcar.Egap.at(0) > _METALGAP_) {
      SPIN_UP = true;
    }
    if (xoutcar.Egap.at(1) > _METALGAP_) {
      SPIN_DN = true;
    }
  }
  // GET THE OUTPUT
  if (SPIN_UP or SPIN_DN) {
    oss << "System                 :  " << xoutcar.SYSTEM << endl;
    oss << "Number of records      :  " << xoutcar.carrier_type.size() << endl;
    for (size_t itr0 = 0; itr0 < xoutcar.carrier_type.size(); itr0++) {
      oss << std::noshowpos;
      oss << "** Record  " << itr0 + 1 << endl;
      oss << "   Band index          :  " << xoutcar.band_index.at(itr0) << endl;
      oss << "   Carrier type        :  " << xoutcar.carrier_type[itr0] << endl;
      oss << "   Spin type           :  " << xoutcar.carrier_spin.at(itr0) << endl;
      oss << "   Extrema coord       : ";
      oss << setw(14) << std::left << std::scientific << std::showpos << xoutcar.extrema_cart_coord.at(itr0).at(0);
      oss << setw(14) << std::left << std::scientific << std::showpos << xoutcar.extrema_cart_coord.at(itr0).at(1);
      oss << setw(14) << std::left << std::scientific << std::showpos << xoutcar.extrema_cart_coord.at(itr0).at(2);
      oss << endl;
      oss << "   Principal axes      : ";
      oss << setw(14) << std::left << std::scientific << std::showpos << xoutcar.effective_mass_axes.at(itr0).at(0);
      oss << setw(14) << std::left << std::scientific << std::showpos << xoutcar.effective_mass_axes.at(itr0).at(1);
      oss << setw(14) << std::left << std::scientific << std::showpos << xoutcar.effective_mass_axes.at(itr0).at(2);
      oss << endl;
      oss << std::noshowpos;
      oss << "   Equivalent valleys  :  " << xoutcar.equivalent_valley.at(itr0) << endl;
      oss << std::showpos;
      oss << "   DOS  eff. mass.     : " << xoutcar.effective_mass_DOS.at(itr0) << endl;
      oss << "   COND eff. mass.     : " << xoutcar.effective_mass_COND.at(itr0) << endl;
    }
    oss << "** Carrier Masses " << endl;
    oss << "   DOS  elec eff. mass : " << setw(14) << std::left << std::scientific << std::showpos << xoutcar.mass_elec_dos.at(0) << endl;
    oss << "   DOS  hole eff. mass : " << setw(14) << std::left << std::scientific << std::showpos << xoutcar.mass_hole_dos.at(0) << endl;
    oss << "   COND elec eff. mass : " << setw(14) << std::left << std::scientific << std::showpos << xoutcar.mass_elec_conduction.at(0) << endl;
    oss << "   COND hole eff. mass : " << setw(14) << std::left << std::scientific << std::showpos << xoutcar.mass_hole_conduction.at(0) << endl;
    oss << endl;
  }
  return true;
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// PrintEigCurv
// Subroutine that obtains band gap extrema curvatures of the band structure.
// This is then used to approximate the location of the e.m. ellipsoids in the IBZ
// 2015: Camilo E. Calderon
bool PrintEigCurv(string& directory, ostream& oss) {
  stringstream message;

  vector<vector<vector<vector<vector<double>>>>> branches_bnds;
  vector<vector<vector<xvector<double>>>> branches_kpts;
  vector<vector<xvector<int>>> branches;
  vector<vector<vector<int>>> branches_indx;
  vector<xvector<double>> special_kpts;
  vector<xvector<double>> unique_kpts;
  vector<xvector<double>> unique_kpts_EIG;
  vector<xvector<double>> vkpoint_eig;
  vector<xvector<int>> repeat_kpts_EIG;
  vector<xvector<int>> connect_kpts;
  vector<xvector<int>> connect_kpts_EIG;
  vector<xvector<int>> vrtx_path;
  vector<xvector<int>> ndx_edges;
  vector<int> connect_kpts_num;
  vector<int> repeat_kpts_num;
  stringstream file_EIGENVAL;
  stringstream file_KPOINTS;
  stringstream file_OUTCAR;
  string path_POSCAR;
  const double _ZERO = 0.0;
  int GRIDS;
  xEIGENVAL xeigenval;
  xOUTCAR xoutcar;
  const char LastChar = *directory.rbegin();
  if (LastChar == '/') {
    directory.erase(directory.size() - 1);
  }

  const vector<string> vext{".bz2", ".xz", ".gz"};

  // EIGENVAL
  bool found_EIGENVAL = false;
  for (size_t iext = 0; iext < vext.size(); iext++) {
    if (!found_EIGENVAL && aurostd::FileExist(directory + "/EIGENVAL.bands" + vext[iext])) {
      found_EIGENVAL = true;
      aurostd::compressfile2stringstream(directory + "/EIGENVAL.bands" + vext[iext], file_EIGENVAL); // .EXT
    }
  }
  if (!found_EIGENVAL && aurostd::FileExist(directory + "/EIGENVAL.bands")) {
    found_EIGENVAL = true;
    aurostd::file2stringstream(directory + "/EIGENVAL.bands", file_EIGENVAL); // plain text
  }
  if (!found_EIGENVAL) {
    return false;
  }

  // KPOINTS
  bool found_KPOINTS = false;
  for (size_t iext = 0; iext < vext.size(); iext++) {
    if (!found_KPOINTS && aurostd::FileExist(directory + "/KPOINTS.bands" + vext[iext])) {
      found_KPOINTS = true;
      aurostd::compressfile2stringstream(directory + "/KPOINTS.bands" + vext[iext], file_KPOINTS); // .EXT
    }
  }
  if (!found_KPOINTS && aurostd::FileExist(directory + "/KPOINTS.bands")) {
    found_KPOINTS = true;
    aurostd::file2stringstream(directory + "/KPOINTS.bands", file_KPOINTS); // plain text
  }
  if (!found_KPOINTS) {
    return false;
  }

  // OUTCAR
  bool found_OUTCAR = false;
  for (size_t iext = 0; iext < vext.size(); iext++) {
    if (!found_OUTCAR && aurostd::FileExist(directory + "/OUTCAR.bands" + vext[iext])) {
      found_OUTCAR = true;
      aurostd::compressfile2stringstream(directory + "/OUTCAR.bands" + vext[iext], file_OUTCAR); // .EXT
    }
  }
  if (!found_OUTCAR && aurostd::FileExist(directory + "/OUTCAR.bands")) {
    found_OUTCAR = true;
    aurostd::file2stringstream(directory + "/OUTCAR.bands", file_OUTCAR); // plain text
  }
  if (!found_OUTCAR) {
    return false;
  }

  // POSCAR
  bool found_POSCAR = false;
  for (size_t iext = 0; iext < vext.size(); iext++) {
    if (!found_POSCAR && aurostd::FileExist(directory + "/POSCAR.bands" + vext[iext])) {
      found_POSCAR = true;
      path_POSCAR = directory + "/POSCAR.bands" + vext[iext]; // EXR
    }
  }
  if (!found_POSCAR && aurostd::FileExist(directory + "/POSCAR.bands")) {
    found_POSCAR = true;
    path_POSCAR = directory + "/POSCAR.bands"; // plain text
  }
  if (!found_POSCAR) {
    return false;
  }

  // now have all files

  xeigenval.GetProperties(file_EIGENVAL);
  xoutcar.GetProperties(file_OUTCAR);

  xstructure xstr(path_POSCAR, IOVASP_POSCAR);

  // CO20171002 - using tolerance from symmetry calc - START
  double tol;
  if (xstr.CalculateSymmetry()) {
    tol = xstr.sym_eps;
  } else {
    tol = SYM::defaultTolerance(xstr);
  }
  // CO20171002 - using tolerance from symmetry calc - STOP

  xoutcar.GetBandGap(tol);
  if (xoutcar.Egap_net >= _ZERO) {
    ParseKPOINTS(file_KPOINTS, GRIDS, special_kpts, unique_kpts, repeat_kpts_num);
    AdjacencyList_KPT(special_kpts, unique_kpts, connect_kpts, connect_kpts_num);
    AdjacencyList_EIG(unique_kpts, connect_kpts, connect_kpts_num, xeigenval, unique_kpts_EIG, connect_kpts_EIG, vkpoint_eig);
    RepeatsList(unique_kpts_EIG, repeat_kpts_num, vkpoint_eig, repeat_kpts_EIG);
    VertexPaths(repeat_kpts_EIG, connect_kpts_EIG, repeat_kpts_num, GRIDS, vrtx_path);
    RepeatedEdges(vrtx_path, repeat_kpts_EIG, repeat_kpts_num, ndx_edges);
    VertexBranches(ndx_edges, repeat_kpts_num, repeat_kpts_EIG, branches);
    PathDataStuct(xeigenval, vkpoint_eig, branches, branches_indx, branches_kpts, branches_bnds);
    IBZextrema(xeigenval, vkpoint_eig, branches);
  } else {
    oss << endl;
    message << "Material is metallic"; // CO20200404
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, oss, _LOGGER_WARNING_); // CO20200404
    return false;
  }
  return true;
}
// -----------------------------------------------------------------------------------------------
// ParseKPOINTS
//
// Extract & clean up the KPOINTS.bands information.
// GRIDS            : line density in KPOINTS.bands
// special_kpts     : high symmetry points in KPOINTS.bands
// unique_kpts      : unique special kpts - get rid of repeated instances of each kpt
// repeat_kpts_KPTS : rows contain the indices of where each unique_kpts shows up
// repeat_kpts_num  : size of each row found in repeat_kpts
//
// Camilo E. Calderon, 2015
// -----------------------------------------------------------------------------------------------
bool ParseKPOINTS(stringstream& file_KPOINTS, int& GRIDS, vector<xvector<double>>& special_kpts, vector<xvector<double>>& unique_kpts, vector<int>& repeat_kpts_num) {
  vector<string> StringKpts;
  vector<string> tokens;
  aurostd::string2vectorstring(file_KPOINTS.str(), StringKpts);
  const string line0 = StringKpts.at(2);
  aurostd::string2tokens(line0, tokens);
  if (tokens.at(0).at(0) != 'L') {
    cout << "ParseKPOINTS - KPOINTS file is not in Line-mode: " << endl;
    return false;
  } else {
    aurostd::string2tokens(StringKpts.at(1), tokens);
    GRIDS = aurostd::string2utype<int>(tokens.at(0));
    uint itr1 = 0;
    if (!StringKpts.empty()) {
      for (size_t itr0 = 4; itr0 < StringKpts.size(); itr0++) {
        const xvector<double> tempvec(4);
        const string line1 = StringKpts[itr0];
        aurostd::string2tokens(line1, tokens);
        if (!tokens.empty()) {
          special_kpts.push_back(tempvec);
          special_kpts.at(itr1)[1] = aurostd::string2utype<double>(tokens.at(0));
          special_kpts.at(itr1)[2] = aurostd::string2utype<double>(tokens.at(1));
          special_kpts.at(itr1)[3] = aurostd::string2utype<double>(tokens.at(2));
          special_kpts.at(itr1)[4] = (double) itr1;
          itr1++;
        }
      }
    } else {
      cout << "ParseKPOINTS - No strings found in the KPOINTS file: " << endl;
      return false;
    }
  }
  vector<xvector<double>> tmp_special_kpts;
  for (size_t itr0 = 0; itr0 < special_kpts.size(); itr0++) {
    const xvector<double> tempvec(4);
    tmp_special_kpts.push_back(tempvec);
    tmp_special_kpts.at(itr0)[1] = special_kpts.at(itr0)[1];
    tmp_special_kpts.at(itr0)[2] = special_kpts.at(itr0)[2];
    tmp_special_kpts.at(itr0)[3] = special_kpts.at(itr0)[3];
    tmp_special_kpts.at(itr0)[4] = special_kpts.at(itr0)[4];
  }
  uint itr2 = 0;
  for (size_t itr0 = 0; itr0 < special_kpts.size(); itr0++) {
    if (itr2 >= special_kpts.size()) {
      break;
    }
    uint delcnt = 0;
    unique_kpts.push_back(tmp_special_kpts.at(itr0));
    tmp_special_kpts.erase(tmp_special_kpts.begin() + itr0);
    delcnt += 1;
    for (size_t itr1 = 0; itr1 < tmp_special_kpts.size(); itr1++) {
      int comp = 0;
      for (int itr2 = 1; itr2 <= 3; itr2++) {
        bool MATCH = false;
        CompareDoublesChar(MATCH, unique_kpts.back()[itr2], tmp_special_kpts.at(itr1)[itr2]);
        if (MATCH) {
          comp++;
        } else if (!MATCH) {
          break;
        }
      }
      if (comp == 3) {
        tmp_special_kpts.erase(tmp_special_kpts.begin() + itr1);
        itr1--;
        delcnt++;
      } else {
        continue;
      }
    }
    itr0--;
    itr2 += delcnt;
  }
  for (size_t itr0 = 0; itr0 < unique_kpts.size(); itr0++) {
    xvector<int> tempvec(special_kpts.size());
    const int index1 = unique_kpts[itr0][4];
    int count = 1;
    for (size_t itr1 = 1; itr1 <= unique_kpts.size(); itr1++) {
      tempvec[itr1] = 99999;
    }
    tempvec[count] = (int) unique_kpts[itr0][4];
    count++;
    for (size_t itr1 = 0; itr1 < special_kpts.size(); itr1++) {
      const int index2 = special_kpts[itr1][4];
      int compare = 0;
      for (int itr2 = 1; itr2 <= 3; itr2++) {
        bool MATCH = false;
        CompareDoublesChar(MATCH, unique_kpts[itr0][itr2], special_kpts[itr1][itr2]);
        if (MATCH) {
          compare++;
        } else if (!MATCH) {
          break;
        }
      }
      if (compare == 3 and (index1 != index2)) {
        tempvec[count] = index2;
        count++;
      }
    }
    const xvector<int> kptrepeat(count - 1);
    repeat_kpts_num.push_back(count - 1);
  }
  return true;
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// AdjacencyList_KPT
//
// build the connectivity of each unique kpoint found in KPOINTS.bands file
// connect_kpts     : unique_kpt + index of its nearest neighbors
// connect_kpts_num : number of entries found in each row of connect_kpts
//
// Camilo E. Calderon, 2015
//-------------------------------------------------------------------------------------------------
bool AdjacencyList_KPT(vector<xvector<double>>& special_kpts, vector<xvector<double>>& unique_kpts, vector<xvector<int>>& connect_kpts, vector<int>& connect_kpts_num) {
  for (size_t itr0 = 0; itr0 < unique_kpts.size(); itr0++) {
    int count = 1;
    xvector<int> tempvec0(unique_kpts.size() + 1);
    for (size_t itr1 = 1; itr1 <= unique_kpts.size() + 1; itr1++) {
      tempvec0[itr1] = 99999;
    }
    tempvec0[count] = (int) unique_kpts[itr0][4];
    count++;
    for (size_t itr1 = 0; itr1 < special_kpts.size(); itr1++) {
      int compare0 = 0;
      for (int itr2 = 1; itr2 <= 3; itr2++) {
        bool MATCH = false;
        CompareDoublesChar(MATCH, unique_kpts[itr0][itr2], special_kpts[itr1][itr2]);
        if (MATCH) {
          compare0++;
        } else if (!MATCH) {
          break;
        }
      }
      if (compare0 == 3) {
        if ((int) special_kpts[itr1][4] % 2 == 0) {
          for (size_t itr2 = 0; itr2 < unique_kpts.size(); itr2++) {
            int compare1 = 0;
            for (int itr3 = 1; itr3 <= 3; itr3++) {
              bool MATCH = false;
              CompareDoublesChar(MATCH, unique_kpts[itr2][itr3], special_kpts.at(itr1 + 1)[itr3]);
              if (MATCH) {
                compare1++;
              } else if (!MATCH) {
                break;
              }
            }
            if (compare1 == 3) {
              if (count <= (int) unique_kpts.size() + 1) {
                tempvec0[count] = (int) unique_kpts[itr2][4];
                count++;
              } else {
                cout << "GetKPOINTSAdjacencyList - Wrong symmetry" << endl;
                return false;
              }
            }
          }
        } else if ((int) special_kpts[itr1][4] % 2 == 1) {
          for (size_t itr2 = 0; itr2 < unique_kpts.size(); itr2++) {
            int compare1 = 0;
            for (int itr3 = 1; itr3 <= 3; itr3++) {
              bool MATCH = false;
              CompareDoublesChar(MATCH, unique_kpts[itr2][itr3], special_kpts.at(itr1 - 1)[itr3]);
              if (MATCH) {
                compare1++;
              } else if (!MATCH) {
                break;
              }
            }
            if (compare1 == 3) {
              if (count <= (int) unique_kpts.size() + 1) {
                tempvec0[count] = (int) unique_kpts[itr2][4];
                count++;
              } else {
                cout << "GetKPOINTSAdjacencyList - Wrong symmetry" << endl;
                return false;
              }
            }
          }
        }
      }
    }
    connect_kpts_num.push_back(count - 1); // possibly count-2?
    const xvector<int> kptconnect(count - 1);
    connect_kpts.push_back(kptconnect);
    for (uint itr1 = 1; itr1 <= (uint) count - 1; itr1++) {
      if (tempvec0[itr1] < 99999) {
        connect_kpts.at(itr0)[itr1] = tempvec0[itr1];
      }
    }
  }
  return true;
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// AdjacencyList_EIG
// Camilo E. Calderon, 2015
//-------------------------------------------------------------------------------------------------
bool AdjacencyList_EIG(vector<xvector<double>>& unique_kpts,
                       vector<xvector<int>>& connect_kpts,
                       vector<int>& connect_kpts_num,
                       xEIGENVAL& xeigenval,
                       vector<xvector<double>>& unique_kpts_EIG,
                       vector<xvector<int>>& connect_kpts_EIG,
                       vector<xvector<double>>& vkpoint_eig) {
  for (size_t itr0 = 0; itr0 < xeigenval.vkpoint.size(); itr0++) {
    xvector<double> tempvec(3);
    for (int itr1 = 1; itr1 <= 3; itr1++) {
      if (aurostd::abs(xeigenval.vkpoint[itr0][itr1]) <= 1.0E-15) {
        tempvec[itr1] = 0.0;
      } else {
        tempvec[itr1] = xeigenval.vkpoint[itr0][itr1];
      }
    }
    vkpoint_eig.push_back(tempvec);
  }
  uint start = 0;
  for (size_t itr0 = 0; itr0 < unique_kpts.size(); itr0++) {
    const int templen = connect_kpts_num.at(itr0);
    const xvector<int> tempvec1(templen);
    connect_kpts_EIG.push_back(tempvec1);
    for (uint itr1 = 1; itr1 <= (uint) templen; itr1++) {
      connect_kpts_EIG.back()[itr1] = connect_kpts.at(itr0)[itr1];
    }
    uint endval = 0;
    for (size_t itr1 = start; itr1 < vkpoint_eig.size(); itr1++) {
      int compare = 0;
      for (int itr2 = 1; itr2 <= 3; itr2++) {
        bool MATCH = false;
        CompareDoublesChar(MATCH, unique_kpts[itr0][itr2], vkpoint_eig.at(itr1)[itr2]);
        if (MATCH) {
          compare++;
        } else if (!MATCH) {
          break;
        }
      }
      if (compare == 3) {
        const xvector<double> tempvec2(4);
        unique_kpts_EIG.push_back(tempvec2);
        unique_kpts_EIG.back()[1] = vkpoint_eig.at(itr1)[1];
        unique_kpts_EIG.back()[2] = vkpoint_eig.at(itr1)[2];
        unique_kpts_EIG.back()[3] = vkpoint_eig.at(itr1)[3];
        unique_kpts_EIG.back()[4] = (double) itr1;
        endval = itr1;
        break;
      }
    }
    start = endval + 1;
  }
  for (size_t itr0 = 0; itr0 < connect_kpts.size(); itr0++) {
    for (uint itr1 = 1; itr1 <= (uint) connect_kpts_num.at(itr0); itr1++) {
      for (size_t itr2 = 0; itr2 < connect_kpts.size(); itr2++) {
        if (connect_kpts[itr2][1] == connect_kpts[itr0][itr1]) {
          connect_kpts_EIG.at(itr0)[itr1] = (int) unique_kpts_EIG.at(itr2)[4];
        }
      }
    }
  }
  return true;
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// RepeatsList
//
// repeat_kpts_EIG: list of equivalent unique kpts in the EIGENVAL file. Each unique kpt is listed
// per row. These are equivalent indices, i.e. kpt #0 & kpt #GRID might be identical to each other
// The list is ordered in ascending.
//
// Camilo E. Calderon, 2015
//-------------------------------------------------------------------------------------------------
bool RepeatsList(vector<xvector<double>>& unique_kpts_EIG, vector<int>& repeat_kpts_num, vector<xvector<double>>& vkpoint_eig, vector<xvector<int>>& repeat_kpts_EIG) {
  for (size_t itr0 = 0; itr0 < unique_kpts_EIG.size(); itr0++) {
    int count = 1;
    const xvector<int> tempvec(repeat_kpts_num.at(itr0));
    repeat_kpts_EIG.push_back(tempvec);
    repeat_kpts_EIG.back()[count] = (int) unique_kpts_EIG[itr0][4];
    count++;
    for (size_t itr1 = 0; itr1 < vkpoint_eig.size(); itr1++) {
      int compare = 0;
      for (int itr2 = 1; itr2 <= 3; itr2++) {
        bool MATCH = false;
        CompareDoublesChar(MATCH, unique_kpts_EIG[itr0][itr2], vkpoint_eig[itr1][itr2]);
        if (MATCH) {
          compare++;
        } else if (!MATCH) {
          break;
        }
      }
      if ((compare == 3) and (unique_kpts_EIG[itr0][4] != (int) itr1)) {
        if (count <= repeat_kpts_num.at(itr0)) {
          repeat_kpts_EIG.at(itr0)[count] = (int) itr1;
          count++;
        } else {
          cout << "KPOINTS.bands does not match EIGENVAL.bands" << endl;
          return false;
        }
      }
    }
  }
  return true;
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// VertexPaths
//
// vrtx_path: Array of 4-vectors, where an A-B-C vertex is turned into A-B1-B2-C form. Note
// that Bn are equivalent vertex kpoints, while A & C are the neighboring kpoints.
// These are integers point to the k-point index in the xeigenval.vkpoint array and
// the A-B1 & B2-C pairs are always neighboring special kpoints on the vkpoint array.
//
// Camilo E. Calderon, 2015
//-------------------------------------------------------------------------------------------------
bool VertexPaths(vector<xvector<int>>& repeat_kpts_EIG,
                 vector<xvector<int>>& connect_kpts_EIG,
                 vector<int>& repeat_kpts_num,
                 int& GRIDS,
                 // returns
                 vector<xvector<int>>& vrtx_path) {
  vector<xvector<int>> vrtx_list;
  for (size_t itr0 = 0; itr0 < connect_kpts_EIG.size(); itr0++) {
    const uint templen = repeat_kpts_num.at(itr0) + 1;
    if (templen >= 3) {
      for (uint itr1 = 2; itr1 <= templen - 1; itr1++) {
        for (uint itr2 = 2; itr2 <= templen; itr2++) {
          if (itr1 != itr2 and itr2 > itr1) {
            const xvector<int> tempvec(3);
            vrtx_list.push_back(tempvec);
            vrtx_list.back()[1] = connect_kpts_EIG[itr0][itr1];
            vrtx_list.back()[2] = connect_kpts_EIG[itr0][1];
            vrtx_list.back()[3] = connect_kpts_EIG[itr0][itr2];
          }
        }
      }
    }
  }
  const vector<xvector<int>> tmp_vrtx_path;
  for (size_t itr0 = 0; itr0 < vrtx_list.size(); itr0++) {
    xvector<int> vrtx_segments(4);
    if (aurostd::abs(vrtx_list[itr0][1] - vrtx_list[itr0][2]) == GRIDS - 1) {
      vrtx_segments[1] = vrtx_list[itr0][1];
      vrtx_segments[2] = vrtx_list[itr0][2];
    } else {
      for (size_t itr1 = 0; itr1 < repeat_kpts_EIG.size(); itr1++) {
        if (repeat_kpts_EIG[itr1][1] == vrtx_list[itr0][1]) {
          for (size_t itr2 = 0; itr2 < repeat_kpts_EIG.size(); itr2++) {
            if (repeat_kpts_EIG[itr2][1] == vrtx_list[itr0][2]) {
              for (uint itr3 = 1; itr3 <= (uint) repeat_kpts_num.at(itr1); itr3++) {
                for (uint itr4 = 1; itr4 <= (uint) repeat_kpts_num.at(itr2); itr4++) {
                  if (aurostd::abs(repeat_kpts_EIG[itr1][itr3] - repeat_kpts_EIG[itr2][itr4]) == GRIDS - 1) {
                    vrtx_segments[1] = repeat_kpts_EIG[itr1][itr3];
                    vrtx_segments[2] = repeat_kpts_EIG[itr2][itr4];
                  }
                }
              }
            }
          }
        }
      }
    }
    if (aurostd::abs(vrtx_list[itr0][2] - vrtx_list[itr0][3]) == GRIDS - 1) {
      vrtx_segments[3] = vrtx_list[itr0][2];
      vrtx_segments[4] = vrtx_list[itr0][3];
    } else {
      for (size_t itr1 = 0; itr1 < repeat_kpts_EIG.size(); itr1++) {
        if (repeat_kpts_EIG[itr1][1] == vrtx_list[itr0][2]) {
          for (size_t itr2 = 0; itr2 < repeat_kpts_EIG.size(); itr2++) {
            if (repeat_kpts_EIG[itr2][1] == vrtx_list[itr0][3]) {
              for (uint itr3 = 1; itr3 <= (uint) repeat_kpts_num.at(itr1); itr3++) {
                for (uint itr4 = 1; itr4 <= (uint) repeat_kpts_num.at(itr2); itr4++) {
                  if (aurostd::abs(repeat_kpts_EIG[itr1][itr3] - repeat_kpts_EIG[itr2][itr4]) == GRIDS - 1) {
                    vrtx_segments[3] = repeat_kpts_EIG[itr1][itr3];
                    vrtx_segments[4] = repeat_kpts_EIG[itr2][itr4];
                  }
                }
              }
            }
          }
        }
      }
    }
    vrtx_path.push_back(vrtx_segments);
  }
  return true;
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// RepeatedEdges
//
// ndx_edges: List of unique edges found in the vrtx_paths array. Contains a list of each type of
// A-B connection in the BZ paths, indexed according to their vkpoints position. Used for tracking
// related data in the path building routines.
//
// Camilo E. Calderon, 2015
//-------------------------------------------------------------------------------------------------
bool RepeatedEdges(vector<xvector<int>>& vrtx_path,
                   vector<xvector<int>>& repeat_kpts_EIG,
                   vector<int>& repeat_kpts_num,
                   // returns:
                   vector<xvector<int>>& ndx_edges) {
  vector<vector<xvector<int>>> allpairs;
  for (size_t itr0 = 0; itr0 < repeat_kpts_EIG.size() - 1; itr0++) {
    for (size_t itr1 = itr0 + 1; itr1 < repeat_kpts_EIG.size(); itr1++) {
      vector<xvector<int>> temparray;
      for (int itr2 = 1; itr2 <= repeat_kpts_num.at(itr0); itr2++) {
        for (int itr3 = 1; itr3 <= repeat_kpts_num.at(itr1); itr3++) {
          xvector<int> tempvec(2);
          tempvec[1] = repeat_kpts_EIG[itr0][itr2];
          tempvec[2] = repeat_kpts_EIG[itr1][itr3];
          temparray.push_back(tempvec);
        }
      }
      allpairs.push_back(temparray);
    }
  }
  vector<vector<xvector<int>>> ndx_edges_tmp;
  for (size_t itr0 = 0; itr0 < vrtx_path.size(); itr0++) {
    xvector<int> edge1(2);
    xvector<int> edge2(2);
    xvector<int> edge_type(2);
    int count = 0;
    edge1[1] = vrtx_path[itr0][1];
    edge1[2] = vrtx_path[itr0][2];
    edge2[1] = vrtx_path[itr0][3];
    edge2[2] = vrtx_path[itr0][4];
    for (size_t itr1 = 0; itr1 < allpairs.size(); itr1++) {
      for (size_t itr2 = 0; itr2 < allpairs[itr1].size(); itr2++) {
        xvector<int> pair(2);
        pair[1] = allpairs[itr1][itr2][1];
        pair[2] = allpairs[itr1][itr2][2];
        if ((edge1[1] == pair[1] and edge1[2] == pair[2]) or (edge1[1] == pair[2] and edge1[2] == pair[1])) {
          edge_type[1] = itr1;
          count++;
        }
        if ((edge2[1] == pair[1] and edge2[2] == pair[2]) or (edge2[1] == pair[2] and edge2[2] == pair[1])) {
          edge_type[2] = itr1;
          count++;
        }
      }
    }
    if (count == 2) {
      vector<xvector<int>> edge_cur;
      xvector<int> temp1(3);
      temp1[1] = edge_type[1];
      temp1[2] = edge1[1];
      temp1[3] = edge1[2];
      edge_cur.push_back(temp1);
      xvector<int> temp2(3);
      temp2[1] = edge_type[2];
      temp2[2] = edge2[1];
      temp2[3] = edge2[2];
      edge_cur.push_back(temp2);
      ndx_edges_tmp.push_back(edge_cur);
    }
  }
  xvector<int> ndx_edges0(3);
  ndx_edges0 = ndx_edges_tmp.at(0).at(0);
  ndx_edges.push_back(ndx_edges0);
  for (size_t itr0 = 0; itr0 < ndx_edges_tmp.size(); itr0++) {
    for (size_t itr1 = 0; itr1 < ndx_edges_tmp[itr0].size(); itr1++) {
      bool EDGE = false;
      for (size_t itr2 = 0; itr2 < ndx_edges.size(); itr2++) {
        if (ndx_edges[itr2][1] == ndx_edges_tmp[itr0][itr1][1]) {
          EDGE = true;
          break;
        }
      }
      if (!EDGE) {
        ndx_edges.push_back(ndx_edges_tmp[itr0][itr1]);
      }
    }
  }
  ndx_edges_tmp.clear();
  return true;
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// VertexBranches
//
// Build a data structure that contains every single connection associated with a given vertex
// found in the vrtx_path list
//
// Camilo E. Calderon, 2015
//-------------------------------------------------------------------------------------------------
bool VertexBranches(vector<xvector<int>>& ndx_edges,
                    vector<int>& repeat_kpts_num,
                    vector<xvector<int>>& repeat_kpts_EIG,
                    // returns:
                    vector<vector<xvector<int>>>& branches) {
  vector<int> ndx_edges_row;
  bool NEWEDGE = false; // CAMILOFIX
  for (size_t itr0 = 0; itr0 <= ndx_edges.size(); itr0++) { // NOT A TYPO!!!
    vector<xvector<int>> vertex_edges;
    // bool NEWEDGE ; // CAMILOFIX
    int repeat_kpts_row;
    if (ndx_edges_row.empty()) {
      for (size_t itr1 = 0; itr1 < repeat_kpts_EIG.size(); itr1++) {
        for (uint itr2 = 1; itr2 <= (uint) repeat_kpts_num.at(itr1); itr2++) {
          if (ndx_edges[itr0][2] == repeat_kpts_EIG[itr1][itr2]) {
            xvector<int> edge_crnt(2);
            if (ndx_edges[itr0][2] == repeat_kpts_EIG[itr1][itr2]) {
              edge_crnt[1] = ndx_edges[itr0][2];
              edge_crnt[2] = ndx_edges[itr0][3];
            } else if (ndx_edges[itr0][2] != repeat_kpts_EIG[itr1][itr2]) {
              edge_crnt[1] = ndx_edges[itr0][3];
              edge_crnt[2] = ndx_edges[itr0][2];
            }
            ndx_edges_row.push_back(ndx_edges[itr0][1]);
            vertex_edges.push_back(edge_crnt);
            repeat_kpts_row = itr1;
            NEWEDGE = true;
            goto EDGE2VERTICES;
          }
        }
      }
    } else if (!ndx_edges_row.empty()) {
      for (size_t itr1 = 0; itr1 < ndx_edges.size(); itr1++) {
        xvector<int> edge_crnt(2);
        bool MATCH1 = false;
        for (size_t itr2 = 0; itr2 < ndx_edges_row.size(); itr2++) {
          if (ndx_edges[itr1][1] != ndx_edges_row[itr2]) {
            edge_crnt[1] = ndx_edges[itr1][2];
            edge_crnt[2] = ndx_edges[itr1][3];
            CompareEdges(branches, vertex_edges, edge_crnt, MATCH1);
            if (!MATCH1) {
              for (size_t itr3 = 0; itr3 < repeat_kpts_EIG.size(); itr3++) {
                for (uint itr4 = 1; itr4 <= (uint) repeat_kpts_num.at(itr3); itr4++) {
                  if (edge_crnt[1] == repeat_kpts_EIG[itr3][itr4]) {
                    ndx_edges_row.push_back(ndx_edges[itr1][1]);
                    vertex_edges.push_back(edge_crnt);
                    repeat_kpts_row = itr3;
                    NEWEDGE = true;
                    goto EDGE2VERTICES;
                  }
                }
              }
            } else if (MATCH1) {
              xvector<int> edge_crnt_new(2);
              bool MATCH2 = false;
              edge_crnt_new[1] = ndx_edges[itr1][3];
              edge_crnt_new[2] = ndx_edges[itr1][2];
              CompareEdges(branches, vertex_edges, edge_crnt_new, MATCH2);
              if (!MATCH2) {
                for (size_t itr3 = 0; itr3 < repeat_kpts_EIG.size(); itr3++) {
                  for (uint itr4 = 1; itr4 <= (uint) repeat_kpts_num.at(itr3); itr4++) {
                    if (edge_crnt_new[1] == repeat_kpts_EIG[itr3][itr4]) {
                      ndx_edges_row.push_back(ndx_edges[itr1][1]);
                      vertex_edges.push_back(edge_crnt_new);
                      repeat_kpts_row = itr3;
                      NEWEDGE = true;
                      goto EDGE2VERTICES;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  EDGE2VERTICES:
    if (NEWEDGE) {
      for (size_t itr1 = 0; itr1 < ndx_edges.size(); itr1++) {
        for (size_t itr2 = 0; itr2 < ndx_edges_row.size(); itr2++) {
          if (ndx_edges_row[itr2] != ndx_edges[itr1][1]) {
            for (int itr3 = 2; itr3 <= 3; itr3++) {
              for (int itr4 = 1; itr4 <= repeat_kpts_num.at(repeat_kpts_row); itr4++) {
                if (repeat_kpts_EIG.at(repeat_kpts_row)[itr4] == ndx_edges[itr1][itr3]) {
                  xvector<int> edge_crnt(2);
                  if (repeat_kpts_EIG.at(repeat_kpts_row)[itr4] == ndx_edges[itr1][2]) {
                    edge_crnt[1] = ndx_edges[itr1][2];
                    edge_crnt[2] = ndx_edges[itr1][3];
                  } else if (repeat_kpts_EIG.at(repeat_kpts_row)[itr4] == ndx_edges[itr1][3]) {
                    edge_crnt[1] = ndx_edges[itr1][3];
                    edge_crnt[2] = ndx_edges[itr1][2];
                  }
                  bool MATCH3 = false;
                  CompareEdges(branches, vertex_edges, edge_crnt, MATCH3);
                  if (!MATCH3) {
                    vertex_edges.push_back(edge_crnt);
                  }
                  NEWEDGE = false;
                }
              }
            }
          } else if (ndx_edges_row[itr2] == ndx_edges[itr1][1]) {
            for (int itr3 = 2; itr3 <= 3; itr3++) {
              for (int itr4 = 1; itr4 <= repeat_kpts_num.at(repeat_kpts_row); itr4++) {
                if (repeat_kpts_EIG.at(repeat_kpts_row)[itr4] == ndx_edges[itr1][itr3]) {
                  xvector<int> edge_crnt(2);
                  bool MATCH3 = false;
                  edge_crnt[1] = ndx_edges[itr1][3];
                  edge_crnt[2] = ndx_edges[itr1][2];
                  CompareEdges(branches, vertex_edges, edge_crnt, MATCH3);
                  if (!MATCH3) {
                    for (uint itr5 = 1; itr5 <= (uint) repeat_kpts_num.at(repeat_kpts_row); itr5++) {
                      if (repeat_kpts_EIG.at(repeat_kpts_row)[itr5] == edge_crnt[1]) {
                        vertex_edges.push_back(edge_crnt);
                      }
                    }
                  }
                  NEWEDGE = false;
                }
              }
            }
          }
        }
      }
      branches.push_back(vertex_edges);
    }
  }
  return true;
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// PathDataStuct
//
// Produce the basic data structures needed for band curvature determinations.
// branches_bnds: all energies associated with a given branch-edge associated with a vertex
// branches_kpts: positions of the kpoints
// branches_indx: xeigenval.vkpoint indices of the kpoints associated with the edge
// POSSIBLY USELESS BUT KEEP FOR NOW
//
// Camilo E. Calderon, 2015
//-------------------------------------------------------------------------------------------------
bool PathDataStuct(xEIGENVAL& xeigenval,
                   vector<xvector<double>>& vkpoint_eig,
                   vector<vector<xvector<int>>>& branches,
                   // returns:
                   vector<vector<vector<int>>>& branches_indx,
                   vector<vector<vector<xvector<double>>>>& branches_kpts,
                   vector<vector<vector<vector<vector<double>>>>>& branches_bnds) {
  for (size_t itr0 = 0; itr0 < branches.size(); itr0++) {
    vector<vector<vector<vector<double>>>> ener_edge;
    vector<vector<xvector<double>>> kpts_edge;
    vector<vector<int>> indx_edge;
    for (size_t itr1 = 0; itr1 < branches[itr0].size(); itr1++) {
      if (branches[itr0][itr1][1] > branches[itr0][itr1][2]) {
        vector<vector<vector<double>>> ener_list;
        vector<xvector<double>> kpts_list;
        vector<int> indx_list;
        for (int itr2 = branches[itr0][itr1][1]; itr2 >= branches[itr0][itr1][2]; itr2--) {
          const xvector<double> kpts_cart = vkpoint_eig.at(itr2);
          vector<vector<double>> ener_band;
          kpts_list.push_back(kpts_cart);
          indx_list.push_back(itr2);
          for (uint itr3 = 0; itr3 < xeigenval.number_bands; itr3++) {
            vector<double> ener_spin;
            for (uint itr4 = 0; itr4 < (uint) (xeigenval.spin + 1); itr4++) {
              ener_spin.push_back(xeigenval.venergy.at(itr2).at(itr3).at(itr4));
            }
            ener_band.push_back(ener_spin);
          }
          ener_list.push_back(ener_band);
        }
        ener_edge.push_back(ener_list);
        indx_edge.push_back(indx_list);
        kpts_edge.push_back(kpts_list);
      } else if (branches[itr0][itr1][1] < branches[itr0][itr1][2]) {
        vector<vector<vector<double>>> ener_list;
        vector<xvector<double>> kpts_list;
        vector<int> indx_list;
        for (int itr2 = branches[itr0][itr1][1]; itr2 <= branches[itr0][itr1][2]; itr2++) {
          const xvector<double> kpts_cart = vkpoint_eig.at(itr2);
          vector<vector<double>> ener_band;
          kpts_list.push_back(kpts_cart);
          indx_list.push_back(itr2);
          for (uint itr3 = 0; itr3 < xeigenval.number_bands; itr3++) {
            vector<double> ener_spin;
            for (uint itr4 = 0; itr4 < (uint) (xeigenval.spin + 1); itr4++) {
              ener_spin.push_back(xeigenval.venergy.at(itr2).at(itr3).at(itr4));
            }
            ener_band.push_back(ener_spin);
          }
          ener_list.push_back(ener_band);
        }
        ener_edge.push_back(ener_list);
        indx_edge.push_back(indx_list);
        kpts_edge.push_back(kpts_list);
      }
    }
    branches_bnds.push_back(ener_edge);
    branches_indx.push_back(indx_edge);
    branches_kpts.push_back(kpts_edge);
  }
  return true;
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// IBZextrema
//
// Generate lists of all the IBZ band extrema & their curvatures.
// NAIVE CURVATURES: PROJECT GRIDS TO 1 DIMENSION, AND USE O(4) CENTRAL DIFFERENCE
// REFINE THIS FURTHER
//
// Camilo E. Calderon, 2015
//-------------------------------------------------------------------------------------------------
bool IBZextrema(xEIGENVAL& xeigenval, vector<xvector<double>>& vkpoint_eig, vector<vector<xvector<int>>>& branches) {
  vector<double> all_curves;
  vector<double> curvature;
  for (size_t itr0 = 0; itr0 < branches.size(); itr0++) {
    for (size_t itr1 = 0; itr1 < branches[itr0].size(); itr1++) {
      vector<xvector<int>> beg_edge_ndx;
      vector<xvector<int>> end_edge_ndx;
      for (size_t itr2 = 0; itr2 < branches.size(); itr2++) {
        for (size_t itr3 = 0; itr3 < branches[itr2].size(); itr3++) {
          if (branches[itr2][itr3][1] == branches[itr0][itr1][1]) {
            for (size_t itr4 = 0; itr4 < branches[itr2].size(); itr4++) {
              if (itr3 != itr4) {
                beg_edge_ndx.push_back(branches[itr2][itr4]);
              }
            }
          }
          if (branches[itr2][itr3][1] == branches[itr0][itr1][2]) {
            for (size_t itr4 = 0; itr4 < branches[itr2].size(); itr4++) {
              if (itr3 != itr4) {
                end_edge_ndx.push_back(branches[itr2][itr4]);
              }
            }
          }
        }
      }
      // edge curvatures
      vector<xvector<double>> posvec;
      // for(uint itr3=0; itr3<xeigenval.number_bands; itr3++)
      for (uint itr3 = 0; itr3 < 1; itr3++) { // CO20200106 - patching for auto-indenting
        // for(uint itr4=0; itr4<xeigenval.spin; itr4++)
        for (uint itr4 = 0; itr4 < 1; itr4++) { // CO20200106 - patching for auto-indenting
          xvector<double> eigvec(5);
          xvector<int> ndxvec(5);
          if (branches[itr0][itr1][1] > branches[itr0][itr1][2]) {
            for (int itr2 = branches[itr0][itr1][1] - 2; itr2 >= branches[itr0][itr1][2] + 2; itr2--) {
              ndxvec[1] = itr2 + 2;
              ndxvec[2] = itr2 + 1;
              ndxvec[3] = itr2;
              ndxvec[4] = itr2 - 1;
              ndxvec[5] = itr2 - 2;
              for (int itr2 = 1; itr2 <= 5; itr2++) {
                eigvec[itr2] = xeigenval.venergy.at(ndxvec[itr2]).at(itr3).at(itr4);
                posvec.push_back(vkpoint_eig.at(ndxvec[itr2]));
              }
              NaiveCurvatures(eigvec, posvec, curvature);
              for (size_t itr3 = 0; itr3 < curvature.size(); itr3++) {
                all_curves.push_back(curvature[itr3]);
              }
              curvature.clear();
              posvec.clear();
            }
            // ME20200724 - there used to be an exit here without comment or explanation
            const string message = "branches[itr0][itr1][1] > branches[itr0][itr1][2])";
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
          } else if (branches[itr0][itr1][1] < branches[itr0][itr1][2]) {
            for (int itr2 = branches[itr0][itr1][1] + 2; itr2 <= branches[itr0][itr1][2] - 2; itr2++) {
              ndxvec[1] = itr2 - 2;
              ndxvec[2] = itr2 - 1;
              ndxvec[3] = itr2;
              ndxvec[4] = itr2 + 1;
              ndxvec[5] = itr2 + 2;
              for (int itr2 = 1; itr2 <= 5; itr2++) {
                eigvec[itr2] = xeigenval.venergy.at(ndxvec[itr2]).at(itr3).at(itr4);
                posvec.push_back(vkpoint_eig.at(ndxvec[itr2]));
              }
              NaiveCurvatures(eigvec, posvec, curvature);
              for (size_t itr3 = 0; itr3 < curvature.size(); itr3++) {
                all_curves.push_back(curvature[itr3]);
              }
              curvature.clear();
              posvec.clear();
            }
          }
          // vertex curvatures
          if (!beg_edge_ndx.empty()) {
            for (size_t itr2 = 0; itr2 < beg_edge_ndx.size(); itr2++) {
              xvector<double> eigvec(6);
              xvector<int> ndxvec(6);
              if (beg_edge_ndx[itr2][1] < beg_edge_ndx[itr2][2]) {
                ndxvec[1] = beg_edge_ndx[itr2][1] + 2;
                ndxvec[2] = beg_edge_ndx[itr2][1] + 1;
                ndxvec[3] = branches[itr0][itr1][1];
                if (branches[itr0][itr1][1] > branches[itr0][itr1][2]) {
                  ndxvec[4] = branches[itr0][itr1][1] - 1;
                  ndxvec[5] = branches[itr0][itr1][1] - 2;
                  ndxvec[6] = branches[itr0][itr1][1] - 3;
                } else if (branches[itr0][itr1][1] < branches[itr0][itr1][2]) {
                  ndxvec[4] = branches[itr0][itr1][1] + 1;
                  ndxvec[5] = branches[itr0][itr1][1] + 2;
                  ndxvec[6] = branches[itr0][itr1][1] + 3;
                }
                for (int itr5 = 1; itr5 <= 6; itr5++) {
                  eigvec[itr5] = xeigenval.venergy.at(ndxvec[itr5]).at(itr3).at(itr4);
                  posvec.push_back(vkpoint_eig.at(ndxvec[itr5]));
                }
                NaiveCurvatures(eigvec, posvec, curvature);
                for (size_t itr3 = 0; itr3 < curvature.size(); itr3++) {
                  all_curves.push_back(curvature[itr3]);
                }
                curvature.clear();
                posvec.clear();
              } else if (beg_edge_ndx[itr2][1] > beg_edge_ndx[itr2][2]) {
                ndxvec[1] = beg_edge_ndx[itr2][1] - 2;
                ndxvec[2] = beg_edge_ndx[itr2][1] - 1;
                ndxvec[3] = branches[itr0][itr1][1];
                if (branches[itr0][itr1][1] > branches[itr0][itr1][2]) {
                  ndxvec[4] = branches[itr0][itr1][1] - 1;
                  ndxvec[5] = branches[itr0][itr1][1] - 2;
                  ndxvec[6] = branches[itr0][itr1][1] - 3;
                } else if (branches[itr0][itr1][1] < branches[itr0][itr1][2]) {
                  ndxvec[4] = branches[itr0][itr1][1] + 1;
                  ndxvec[5] = branches[itr0][itr1][1] + 2;
                  ndxvec[6] = branches[itr0][itr1][1] + 3;
                }
                for (int itr5 = 1; itr5 <= 6; itr5++) {
                  eigvec[itr5] = xeigenval.venergy.at(ndxvec[itr5]).at(itr3).at(itr4);
                  posvec.push_back(vkpoint_eig.at(ndxvec[itr5]));
                }
                NaiveCurvatures(eigvec, posvec, curvature);
                for (size_t itr3 = 0; itr3 < curvature.size(); itr3++) {
                  all_curves.push_back(curvature[itr3]);
                }
                curvature.clear();
                posvec.clear();
              }
            }
          } else if (beg_edge_ndx.empty()) { // if no branches @ beg
            xvector<double> eigvec(6);
            xvector<int> ndxvec(6);
            if (branches[itr0][itr1][1] > branches[itr0][itr1][2]) {
              ndxvec[1] = branches[itr0][itr1][1] - 3;
              ndxvec[2] = branches[itr0][itr1][1] - 2;
              ndxvec[3] = branches[itr0][itr1][1] - 1;
              ndxvec[4] = branches[itr0][itr1][1];
              ndxvec[5] = branches[itr0][itr1][1] - 1;
              ndxvec[6] = branches[itr0][itr1][1] - 2;
              for (int itr5 = 1; itr5 <= 6; itr5++) {
                eigvec[itr5] = xeigenval.venergy.at(ndxvec[itr5]).at(itr3).at(itr4);
                posvec.push_back(vkpoint_eig.at(ndxvec[itr5]));
              }
              NaiveCurvatures(eigvec, posvec, curvature);
              for (size_t itr3 = 0; itr3 < curvature.size(); itr3++) {
                all_curves.push_back(curvature[itr3]);
              }
              curvature.clear();
              posvec.clear();
            } else if (branches[itr0][itr1][1] < branches[itr0][itr1][2]) {
              ndxvec[1] = branches[itr0][itr1][1] + 3;
              ndxvec[2] = branches[itr0][itr1][1] + 2;
              ndxvec[3] = branches[itr0][itr1][1] + 1;
              ndxvec[4] = branches[itr0][itr1][1];
              ndxvec[5] = branches[itr0][itr1][1] + 1;
              ndxvec[6] = branches[itr0][itr1][1] + 2;
              for (int itr5 = 1; itr5 <= 6; itr5++) {
                eigvec[itr5] = xeigenval.venergy.at(ndxvec[itr5]).at(itr3).at(itr4);
                posvec.push_back(vkpoint_eig.at(ndxvec[itr5]));
              }
              NaiveCurvatures(eigvec, posvec, curvature);
              for (size_t itr3 = 0; itr3 < curvature.size(); itr3++) {
                all_curves.push_back(curvature[itr3]);
              }
              curvature.clear();
              posvec.clear();
            }
          }
          if (!end_edge_ndx.empty()) {
            for (size_t itr2 = 0; itr2 < end_edge_ndx.size(); itr2++) {
              xvector<double> eigvec(6);
              xvector<int> ndxvec(6);
              if (end_edge_ndx[itr2][1] < end_edge_ndx[itr2][2]) {
                if (branches[itr0][itr1][1] > branches[itr0][itr1][2]) {
                  ndxvec[1] = branches[itr0][itr1][2] + 3;
                  ndxvec[2] = branches[itr0][itr1][2] + 2;
                  ndxvec[3] = branches[itr0][itr1][2] + 1;
                } else if (branches[itr0][itr1][1] < branches[itr0][itr1][2]) {
                  ndxvec[1] = branches[itr0][itr1][2] - 3;
                  ndxvec[2] = branches[itr0][itr1][2] - 2;
                  ndxvec[3] = branches[itr0][itr1][2] - 1;
                }
                ndxvec[4] = branches[itr0][itr1][2];
                ndxvec[5] = end_edge_ndx[itr2][1] + 1;
                ndxvec[6] = end_edge_ndx[itr2][1] + 2;
                for (int itr5 = 1; itr5 <= 6; itr5++) {
                  eigvec[itr5] = xeigenval.venergy.at(ndxvec[itr5]).at(itr3).at(itr4);
                  posvec.push_back(vkpoint_eig.at(ndxvec[itr5]));
                }
                NaiveCurvatures(eigvec, posvec, curvature);
                for (size_t itr3 = 0; itr3 < curvature.size(); itr3++) {
                  all_curves.push_back(curvature[itr3]);
                }
                curvature.clear();
                posvec.clear();
              } else if (end_edge_ndx[itr2][1] > end_edge_ndx[itr2][2]) {
                if (branches[itr0][itr1][1] > branches[itr0][itr1][2]) {
                  ndxvec[1] = branches[itr0][itr1][2] + 3;
                  ndxvec[2] = branches[itr0][itr1][2] + 2;
                  ndxvec[3] = branches[itr0][itr1][2] + 1;
                } else if (branches[itr0][itr1][1] < branches[itr0][itr1][2]) {
                  ndxvec[1] = branches[itr0][itr1][2] - 3;
                  ndxvec[2] = branches[itr0][itr1][2] - 2;
                  ndxvec[3] = branches[itr0][itr1][2] - 1;
                }
                ndxvec[4] = branches[itr0][itr1][2];
                ndxvec[5] = end_edge_ndx[itr2][1] - 1;
                ndxvec[6] = end_edge_ndx[itr2][1] - 2;
                for (int itr5 = 1; itr5 <= 6; itr5++) {
                  eigvec[itr5] = xeigenval.venergy.at(ndxvec[itr5]).at(itr3).at(itr4);
                  posvec.push_back(vkpoint_eig.at(ndxvec[itr5]));
                }
                NaiveCurvatures(eigvec, posvec, curvature);
                for (size_t itr3 = 0; itr3 < curvature.size(); itr3++) {
                  all_curves.push_back(curvature[itr3]);
                }
                curvature.clear();
                posvec.clear();
              }
            }
          } else if (end_edge_ndx.empty()) { // only in disconnected edges
            xvector<double> eigvec(6);
            xvector<int> ndxvec(6);
            if (branches[itr0][itr1][1] > branches[itr0][itr1][2]) {
              ndxvec[1] = branches[itr0][itr1][1] - 3;
              ndxvec[2] = branches[itr0][itr1][1] - 2;
              ndxvec[3] = branches[itr0][itr1][1] - 1;
              ndxvec[4] = branches[itr0][itr1][1];
              ndxvec[5] = branches[itr0][itr1][1] - 1;
              ndxvec[6] = branches[itr0][itr1][1] - 2;
              for (int itr5 = 1; itr5 <= 6; itr5++) {
                eigvec[itr5] = xeigenval.venergy.at(ndxvec[itr5]).at(itr3).at(itr4);
                posvec.push_back(vkpoint_eig.at(ndxvec[itr5]));
              }
              NaiveCurvatures(eigvec, posvec, curvature);
              for (size_t itr3 = 0; itr3 < curvature.size(); itr3++) {
                all_curves.push_back(curvature[itr3]);
              }
              curvature.clear();
              posvec.clear();
            } else if (branches[itr0][itr1][1] < branches[itr0][itr1][2]) {
              ndxvec[1] = branches[itr0][itr1][1] + 3;
              ndxvec[2] = branches[itr0][itr1][1] + 2;
              ndxvec[3] = branches[itr0][itr1][1] + 1;
              ndxvec[4] = branches[itr0][itr1][1];
              ndxvec[5] = branches[itr0][itr1][1] + 1;
              ndxvec[6] = branches[itr0][itr1][1] + 2;
              for (int itr5 = 1; itr5 <= 6; itr5++) {
                eigvec[itr5] = xeigenval.venergy.at(ndxvec[itr5]).at(itr3).at(itr4);
                posvec.push_back(vkpoint_eig.at(ndxvec[itr5]));
              }
              NaiveCurvatures(eigvec, posvec, curvature);
              for (size_t itr3 = 0; itr3 < curvature.size(); itr3++) {
                all_curves.push_back(curvature[itr3]);
              }
              curvature.clear();
              posvec.clear();
            }
          }
        }
      }
    }
  }
  return true;
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// NaiveCurvatures:
//
// Project all distances onto a false x-axis, average the distances between the points to define
// the value for 'h' and then take the 5-point central difference curvature
//
// Camilo E. Calderon, 2015
//-------------------------------------------------------------------------------------------------
// void NaiveCurvatures(xvector<double> eigvec,
void NaiveCurvatures(xvector<double>& eigvec,
                     vector<xvector<double>>& posvec,
                     // returns:
                     vector<double>& curvature) {
  if (posvec.size() == 5) {
    curvature.push_back(StencilLinear1D(posvec, eigvec));
  } else if (posvec.size() == 6) {
    xvector<double> eigenvals(posvec.size() - 1);
    vector<xvector<double>> positions;
    for (size_t itr0 = 0; itr0 < posvec.size() - 1; itr0++) {
      xvector<double> tempvec(3);
      tempvec[1] = posvec[itr0][1];
      tempvec[2] = posvec[itr0][2];
      tempvec[3] = posvec[itr0][3];
      positions.push_back(tempvec);
      eigenvals[itr0 + 1] = eigvec[itr0 + 1];
    }
    curvature.push_back(StencilLinear1D(positions, eigenvals));
    positions.clear();
    for (size_t itr0 = 1; itr0 < posvec.size(); itr0++) {
      xvector<double> tempvec(3);
      tempvec[1] = posvec[itr0][1];
      tempvec[2] = posvec[itr0][2];
      tempvec[3] = posvec[itr0][3];
      positions.push_back(tempvec);
      eigenvals[itr0] = eigvec[itr0 + 1];
    }
    curvature.push_back(StencilLinear1D(positions, eigenvals));
  }
  return;
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// StencilLinear1D:
//
// Five point linear stencil for O(4) central difference curvatures. Returns a double with the
// value for the curvature at the central point
// Uses: f"(@f3) = (-f1+16f2-30f3+16f4-f5)/(12h^2)
//
// Camilo E. Calderon, 2015
//-------------------------------------------------------------------------------------------------
double StencilLinear1D(vector<xvector<double>>& positions, xvector<double>& eigenvals) {
  xvector<double> numer(5);
  xvector<double> posns(5);
  xvector<double> delta(3);
  double denom;
  double dist = 0;
  posns[1] = 0;
  for (size_t itr0 = 1; itr0 < positions.size(); itr0++) {
    delta[1] = pow((positions.at(itr0 - 1)[1] - positions[itr0][1]), 2.0);
    delta[2] = pow((positions.at(itr0 - 1)[2] - positions[itr0][2]), 2.0);
    delta[3] = pow((positions.at(itr0 - 1)[3] - positions[itr0][3]), 2.0);
    dist += pow((delta[1] + delta[2] + delta[3]), 0.5);
    posns[itr0 + 1] = dist;
  }
  numer[1] = -eigenvals[1];
  numer[2] = 16.0 * eigenvals[2];
  numer[3] = -30.0 * eigenvals[3];
  numer[4] = 16.0 * eigenvals[4];
  numer[5] = -eigenvals[5];
  denom = pow((12.0 * pow(((posns[5] - posns[1]) / (4)), 2.0)), -1.0);
  return (numer[1] + numer[2] + numer[3] + numer[4] + numer[5]) / denom;
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// CompareDoublesChar:
//
// Do ASCII character based comparisons between numbers. Useful for parsing-induced troubles.
//
// Camilo E. Calderon, 2015
//-------------------------------------------------------------------------------------------------
void CompareDoublesChar(bool& MATCH, double& number1, double& number2) {
  ostringstream oss1;
  ostringstream oss2;
  oss1.precision(8);
  oss2.precision(8);
  oss1 << std::scientific << number1;
  oss2 << std::scientific << number2;
  if (oss1.str().size() == oss2.str().size()) {
    const int lencomp = oss1.str().size();
    int compare = 0;
    for (size_t itr3 = 0; itr3 < oss1.str().size(); itr3++) {
      if (oss1.str()[itr3] == oss2.str().at(itr3)) {
        compare++;
      }
    }
    if (compare == lencomp) {
      MATCH = true;
      return;
    } else {
      MATCH = false;
      return;
    }
    // following is probably not needed, but keep it.
  } else if (oss1.str().size() != oss2.str().size()) {
    int decloc1 = -1;
    int decloc2 = -1;
    const int len1 = oss1.str().size();
    const int len2 = oss2.str().size();
    for (uint itr0 = 0; itr0 < (uint) aurostd::min(len1, len2); itr0++) {
      if ((int) oss1.str().at(itr0) == 46) {
        decloc1 = itr0;
      }
      if ((int) oss2.str().at(itr0) == 46) {
        decloc2 = itr0;
      }
    }
    if ((decloc1 == -1 and decloc2 > -1) or (decloc1 > -1 and decloc2 == -1) or (decloc1 == -1 and decloc2 == -1)) {
      MATCH = false;
      return;
    } else if (decloc1 != decloc2) {
      MATCH = false;
      return;
    } else if (decloc1 == decloc2) {
      double min_num = 0.0;
      double max_num = 0.0; // CAMILOFIX
      const int min_len = aurostd::min(len1, len2) - 1;
      const int max_len = aurostd::max(len1, len2) - 1;
      if (len1 == min_len + 1) {
        min_num = number1;
        max_num = number2;
      } else if (len2 == min_len + 1) {
        min_num = number2;
        max_num = number1;
      }
      const double hi_lim = min_num + 4.0 * pow(10, -max_len);
      const double lo_lim = min_num - 5.0 * pow(10, -max_len);
      if ((max_num >= lo_lim) and (max_num < hi_lim)) {
        MATCH = true;
        return;
      } else {
        MATCH = false;
        return;
      }
    }
  }
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// CompareEdges
//
// Checks if current edge is already contained in the "branches" array
// input args are current edge & branches array / returns true / false
//
// Camilo E. Calderon, 2015
//-------------------------------------------------------------------------------------------------
void CompareEdges(vector<vector<xvector<int>>>& branches, vector<xvector<int>>& vertex_edges, xvector<int>& test_edge, bool& COMPARE_EDGES) {
  // branches
  if (!branches.empty()) {
    for (size_t itr0 = 0; itr0 < branches.size(); itr0++) {
      for (size_t itr1 = 0; itr1 < branches[itr0].size(); itr1++) {
        if ((branches[itr0][itr1][1] == test_edge[1]) and (branches[itr0][itr1][2] == test_edge[2])) {
          COMPARE_EDGES = true;
          return;
        }
      }
    }
  } else if (branches.empty()) {
    COMPARE_EDGES = false;
    return;
  }
  // vertex_edges
  if (!vertex_edges.empty()) {
    for (size_t itr0 = 0; itr0 < vertex_edges.size(); itr0++) {
      if ((vertex_edges[itr0][1] == test_edge[1]) and (vertex_edges[itr0][2] == test_edge[2])) {
        COMPARE_EDGES = true;
        return;
      }
    }
  } else if (vertex_edges.empty()) {
    COMPARE_EDGES = false;
    return;
  }
}
// JSON serialization of xEIGENVAL
aurostd::JSON::object xEIGENVAL::serialize() const {
  return aurostd::JSON::object({AST_JSON_GETTER(JSON_xEIGENVAL_MEMBERS)});
}

xEIGENVAL xEIGENVAL::deserialize(const aurostd::JSON::object& jo) {
  AST_JSON_SETTER(JSON_xEIGENVAL_MEMBERS)
  return *this;
}

//-------------------------------------------------------------------------------------------------
// ***************************************************************************
// class xPOTCAR
bool xPOTCAR::GetProperties(const string& stringIN, bool QUIET) {
  stringstream sss;
  sss.str(stringIN);
  if (filename.empty()) {
    filename = "string";
  }
  return xPOTCAR::GetProperties(sss, QUIET);
}

bool xPOTCAR::GetPropertiesFile(const string& fileIN, bool QUIET) {
  stringstream sss;
  if (filename.empty()) {
    filename = fileIN;
  }
  aurostd::compressfile2stringstream(fileIN, sss);
  return xPOTCAR::GetProperties(sss, QUIET);
}

bool xPOTCAR::GetPropertiesUrlFile(const string& url, const string& file, bool QUIET) {
  const string tmpfile = aurostd::TmpFileCreate("xPOTCAR_GetProperties"); // CO20200502 - threadID
  aurostd::httpGetFileStatus(url + "/" + file, tmpfile);
  const bool out = GetPropertiesFile(tmpfile, QUIET);
  filename = "url=" + url; // CO20210315
  aurostd::RemoveFile(tmpfile);
  return out;
}

#define PSEUDOPOTENTIAL_GENERATOR_pad 70

bool xPOTCAR::GetProperties(const stringstream& stringstreamIN, bool QUIET) {
  const bool LDEBUG = (false || XHOST.DEBUG || !QUIET);
  stringstream message;
  const bool force_exit = XHOST.POSTPROCESS; // SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only

  bool ERROR_flag = false;
  // XHOST.PSEUDOPOTENTIAL_GENERATOR=true;
  const long double seconds = aurostd::get_seconds();
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BEGIN (" << time_delay(seconds) << ")" << endl;
  }
  clear(); // so it does not mess up vector/deque
  content = stringstreamIN.str();
  if (true && LDEBUG) {
    cerr << __AFLOW_FUNC__ << endl;
    cerr << content;
  }
  vcontent.clear();
  vector<string> vline;
  vector<string> tokens;
  aurostd::string2vectorstring(content, vcontent);
  if (filename.empty()) {
    filename = "stringstream";
  }
  // crunching to eat the info

  if (vcontent.empty()) {
    return false;
  } // CO+ME20200604 - file does not exist, don't spit out warnings for the workshops

  // GET AUID AS SOON AS POSSIBLE
  if (aurostd::FileExist(filename)) {
    AUID = aurostd::file2auid(filename);
  }

  // get parameters
  //  vline.clear();
  for (size_t iline = 0; iline < vcontent.size(); iline++) {
    aurostd::string2tokens(vcontent[iline], tokens);
    if (iline == 4) {
      title = vcontent[iline]; // might be wrong
    }
  }

  // ----------------------------------------------------------------------
  // EATOM RCORE RWIGS EAUG RAUG ENMAX ENMIN POMASS ZVAL RMAX
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD \"EATOM RCORE RWIGS EAUG RAUG ENMAX ENMIN POMASS ZVAL RMAX\" DATA" << endl;
  }
  vline.clear();
  vENMAX.clear();
  vENMIN.clear();
  vPOMASS.clear();
  vZVAL.clear();
  vEATOM.clear();
  vRCORE.clear();
  vRWIGS.clear();
  vEAUG.clear();
  vRAUG.clear();
  vRMAX.clear();
  for (size_t iline = 0; iline < vcontent.size(); iline++) {
    aurostd::StringSubstInPlace(vcontent[iline], "EMMIN", "ENMIN"); // Kresseeeeeeeeeee
    if (aurostd::substring2bool(vcontent[iline], "ENMAX") && aurostd::substring2bool(vcontent[iline], "ENMIN")) {
      vline.push_back(vcontent[iline]);
    }
    if (aurostd::substring2bool(vcontent[iline], "POMASS") && aurostd::substring2bool(vcontent[iline], "ZVAL")) {
      vline.push_back(vcontent[iline]);
    }
    if (aurostd::substring2bool(vcontent[iline], "EATOM") && aurostd::substring2bool(vcontent[iline], "eV")) {
      vline.push_back(vcontent[iline]);
    }
    if (aurostd::substring2bool(vcontent[iline], "RCORE") && aurostd::substring2bool(vcontent[iline], "radius")) {
      vline.push_back(vcontent[iline]);
    }
    if (aurostd::substring2bool(vcontent[iline], "RWIGS") && aurostd::substring2bool(vcontent[iline], "radius")) {
      vline.push_back(vcontent[iline]);
    }
    if (aurostd::substring2bool(vcontent[iline], "EAUG")) {
      vline.push_back(vcontent[iline]);
    }
    if (aurostd::substring2bool(vcontent[iline], "RAUG") && aurostd::substring2bool(vcontent[iline], "sphere")) {
      vline.push_back(vcontent[iline]);
    }
    if (aurostd::substring2bool(vcontent[iline], "RMAX") && aurostd::substring2bool(vcontent[iline], "radius")) {
      vline.push_back(vcontent[iline]);
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vline.size()=" << vline.size() << endl;
  }
  if (vline.empty()) {
    message << "Wrong number of \"EATOM RCORE RWIGS EAUG RAUG ENMAX ENMIN POMASS ZVAL RMAX\" in POTCAR" << "   filename=[" << filename << "]" << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  }
  for (size_t j = 0; j < vline.size(); j++) {
    aurostd::StringSubstInPlace(vline[j], "=", " ");
    aurostd::StringSubstInPlace(vline[j], ";", " ");
    //    if(LDEBUG) cerr << vline[j] << endl;
    aurostd::string2tokens(vline[j], tokens, " ");
    for (size_t k = 0; k < tokens.size(); k++) {
      if (tokens[k] == "ENMAX" && k + 1 < tokens.size()) {
        vENMAX.push_back(aurostd::string2utype<double>(tokens.at(k + 1)));
      }
      if (tokens[k] == "ENMIN" && k + 1 < tokens.size()) {
        vENMIN.push_back(aurostd::string2utype<double>(tokens.at(k + 1)));
      }
      if (tokens[k] == "POMASS" && k + 1 < tokens.size()) {
        vPOMASS.push_back(aurostd::string2utype<double>(tokens.at(k + 1)));
      }
      if (tokens[k] == "ZVAL" && k + 1 < tokens.size()) {
        vZVAL.push_back(aurostd::string2utype<double>(tokens.at(k + 1)));
      }
      if (tokens[k] == "EATOM" && k + 1 < tokens.size()) {
        vEATOM.push_back(aurostd::string2utype<double>(tokens.at(k + 1)));
      }
      if (tokens[k] == "RCORE" && k + 1 < tokens.size()) {
        vRCORE.push_back(aurostd::string2utype<double>(tokens.at(k + 1)));
      }
      if (tokens[k] == "RWIGS" && k == 0 && k + 1 < tokens.size()) {
        vRWIGS.push_back(aurostd::string2utype<double>(tokens.at(k + 1))); // pick the 1st
      }
      if (tokens[k] == "EAUG" && k + 1 < tokens.size()) {
        vEAUG.push_back(aurostd::string2utype<double>(tokens.at(k + 1)));
      }
      if (tokens[k] == "RAUG" && k + 1 < tokens.size()) {
        vRAUG.push_back(aurostd::string2utype<double>(tokens.at(k + 1)));
      }
      if (tokens[k] == "RMAX" && k + 1 < tokens.size()) {
        vRMAX.push_back(aurostd::string2utype<double>(tokens.at(k + 1)));
      }
    }
  }

  if (vENMIN.size() < vENMAX.size()) {
    vENMIN.push_back(0.0); // ENMIN MIGHT NOT BE THERE
  }
  if (vRCORE.size() < vENMAX.size()) {
    vRCORE.push_back(0.0); // RCORE MIGHT NOT BE THERE
  }
  if (vRWIGS.size() < vENMAX.size()) {
    vRWIGS.push_back(0.0); // RWIGS MIGHT NOT BE THERE
  }
  if (vEAUG.size() < vENMAX.size()) {
    vEAUG.push_back(0.0); // EAUG MIGHT NOT BE THERE
  }
  if (vRAUG.size() < vENMAX.size()) {
    vRAUG.push_back(0.0); // RAUG MIGHT NOT BE THERE
  }
  if (vRMAX.size() < vENMAX.size()) {
    vRMAX.push_back(0.0); // RMAX MIGHT NOT BE THERE
  }

  ENMAX = aurostd::max(vENMAX);
  ENMIN = aurostd::min(vENMIN);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ENMAX=" << ENMAX << " vENMAX.size()=" << vENMAX.size() << ": ";
    for (size_t i = 0; i < vENMAX.size(); i++) {
      cerr << vENMAX[i] << " ";
    }
    cerr << " " << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ENMIN=" << ENMIN << " vENMIN.size()=" << vENMIN.size() << ": ";
    for (size_t i = 0; i < vENMIN.size(); i++) {
      cerr << vENMIN[i] << " ";
    }
    cerr << " " << endl;
  }

  POMASS_sum = aurostd::sum(vPOMASS);
  POMASS_min = aurostd::min(vPOMASS);
  POMASS_max = aurostd::max(vPOMASS);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " POMASS_sum=" << POMASS_sum << " POMASS_min=" << POMASS_min << " POMASS_max=" << POMASS_max << " vPOMASS.size()=" << vPOMASS.size() << ": ";
    for (size_t i = 0; i < vPOMASS.size(); i++) {
      cerr << vPOMASS[i] << " ";
    }
    cerr << " " << endl;
  }

  ZVAL_sum = aurostd::sum(vZVAL);
  ZVAL_min = aurostd::min(vZVAL);
  ZVAL_max = aurostd::max(vZVAL);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ZVAL_sum=" << ZVAL_sum << " ZVAL_min=" << ZVAL_min << " ZVAL_max=" << ZVAL_max << " vZVAL.size()=" << vZVAL.size() << ": ";
    for (size_t i = 0; i < vZVAL.size(); i++) {
      cerr << vZVAL[i] << " ";
    }
    cerr << " " << endl;
  }

  EATOM_min = aurostd::min(vEATOM);
  EATOM_max = aurostd::max(vEATOM);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " EATOM_min=" << EATOM_min << " EATOM_max=" << EATOM_max << " vEATOM.size()=" << vEATOM.size() << ": ";
    for (size_t i = 0; i < vEATOM.size(); i++) {
      cerr << vEATOM[i] << " ";
    }
    cerr << " " << endl;
  }

  RCORE_min = aurostd::min(vRCORE);
  RCORE_max = aurostd::max(vRCORE);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " RCORE_min=" << RCORE_min << " RCORE_max=" << RCORE_max << " vRCORE.size()=" << vRCORE.size() << ": ";
    for (size_t i = 0; i < vRCORE.size(); i++) {
      cerr << vRCORE[i] << " ";
    }
    cerr << " " << endl;
  }

  RWIGS_min = aurostd::min(vRWIGS);
  RWIGS_max = aurostd::max(vRWIGS);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " RWIGS_min=" << RWIGS_min << " RWIGS_max=" << RWIGS_max << " vRWIGS.size()=" << vRWIGS.size() << ": ";
    for (size_t i = 0; i < vRWIGS.size(); i++) {
      cerr << vRWIGS[i] << " ";
    }
    cerr << " " << endl;
  }

  EAUG_min = aurostd::min(vEAUG);
  EAUG_max = aurostd::max(vEAUG);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " EAUG_min=" << EAUG_min << " EAUG_max=" << EAUG_max << " vEAUG.size()=" << vEAUG.size() << ": ";
    for (size_t i = 0; i < vEAUG.size(); i++) {
      cerr << vEAUG[i] << " ";
    }
    cerr << " " << endl;
  }

  RAUG_min = aurostd::min(vRAUG);
  RAUG_max = aurostd::max(vRAUG);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " RAUG_min=" << RAUG_min << " RAUG_max=" << RAUG_max << " vRAUG.size()=" << vRAUG.size() << ": ";
    for (size_t i = 0; i < vRAUG.size(); i++) {
      cerr << vRAUG[i] << " ";
    }
    cerr << " " << endl;
  }

  RMAX_min = aurostd::min(vRMAX);
  RMAX_max = aurostd::max(vRMAX);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " RMAX_min=" << RMAX_min << " RMAX_max=" << RMAX_max << " vRMAX.size()=" << vRMAX.size() << ": ";
    for (size_t i = 0; i < vRMAX.size(); i++) {
      cerr << vRMAX[i] << " ";
    }
    cerr << " " << endl;
  }

  // ----------------------------------------------------------------------
  // PSEUDOPOTENTIAL DATA
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD PSEUDOPOTENTIAL DATA" << endl;
  }

  vline.clear();
  vTITEL.clear();
  vLEXCH.clear();
  for (size_t iline = 0; iline < vcontent.size(); iline++) {
    if (aurostd::substring2bool(vcontent[iline], "=")) {
      if (aurostd::substring2bool(vcontent[iline], "TITEL")) {
        aurostd::string2tokens(vcontent[iline], tokens, "=");
        if (tokens.size() > 1) {
          vTITEL.push_back(tokens.at(1));
        }
        //	if(LDEBUG) cerr << __AFLOW_FUNC__ << " TITEL=" << tokens.at(1) << endl;
      }
      if (aurostd::substring2bool(vcontent[iline], "LEXCH")) {
        if (!aurostd::substring2bool(vcontent[iline], "exchange")) {
          aurostd::string2tokens(vcontent[iline], tokens, "=");
          if (tokens.size() > 1) {
            vLEXCH.push_back(tokens.at(1));
          }
          //	if(LDEBUG) cerr << __AFLOW_FUNC__ << " LEXCH=" << tokens.at(1) << endl;
        }
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vTITEL.size()=" << vTITEL.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vLEXCH.size()=" << vLEXCH.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vRMAX.size()=" << vRMAX.size() << endl;
  }

  // ME20200604 - Only exit if force_exit is true
  if (vTITEL.empty()) {
    message << "Wrong number of pseudopotentials (TITEL) in POTCAR" << "   filename=[" << filename << "]" << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  } // CO20200106 - patching for auto-indenting
  if (vLEXCH.empty()) {
    message << "Wrong number of pseudopotentials (LEXCH) in POTCAR" << "   filename=[" << filename << "]" << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  } // CO20200106 - patching for auto-indenting
  if (vEATOM.empty()) {
    message << "Wrong number of pseudopotentials (EATOM) in POTCAR" << "   filename=[" << filename << "]" << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  } // CO20200106 - patching for auto-indenting
  if (vRMAX.empty()) {
    message << "Wrong number of pseudopotentials (RMAX) in POTCAR" << "   filename=[" << filename << "]" << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  } // CO20200106 - patching for auto-indenting
  if (vTITEL.size() != vLEXCH.size()) {
    message << "Wrong number of pseudopotentials (TITEL/LEXCH) in POTCAR" << "   filename=[" << filename << "]" << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  } // CO20200106 - patching for auto-indenting
  if (vLEXCH.size() != vEATOM.size()) {
    message << "Wrong number of pseudopotentials (LEXCH/EATOM) in POTCAR" << "   filename=[" << filename << "]" << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  } // CO20200106 - patching for auto-indenting
  if (vEATOM.size() != vRMAX.size()) {
    message << "Wrong number of pseudopotentials (EATOM/RMAX) in POTCAR" << "   filename=[" << filename << "]" << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
    ERROR_flag = true;
  } // CO20200106 - patching for auto-indenting

  for (size_t j = 0; j < vTITEL.size(); j++) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " SPECIES(" << j << ") " << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vTITEL[j]=" << vTITEL[j] << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vLEXCH.at(j)=" << vLEXCH.at(j) << endl;
    }
    aurostd::string2tokens(vTITEL[j], tokens, " ");
    pp_type = tokens.at(0);
    if (pp_type == "US" or pp_type == "NC") {
      pp_type = "GGA";
      for (size_t iiline = 0; iiline < vcontent.size(); iiline++) {
        if (aurostd::substring2bool(vcontent[iiline], "GGA")) {
          if (aurostd::substring2bool(vcontent[iiline], "eV")) {
            pp_type = "LDA";
          }
        }
      }
    }
    if (tokens.size() > 1) {
      if (pp_type == "PAW") {
        pp_type = "PAW_LDA"; // cerr << __AFLOW_FUNC__ << " PAW_LDA" << endl;
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " pp_type=" << pp_type << " " << "tokens.size()=" << tokens.size() << endl;
      }
      species.push_back(aurostd::VASP_PseudoPotential_CleanName(tokens.at(1)));
      species_Z.push_back(xelement::symbol2Z(aurostd::VASP_PseudoPotential_CleanName(tokens.at(1))));
      species_pp.push_back(tokens.at(1));
      species_pp_type.push_back(pp_type);
      if (pp_type == "LDA" && tokens.size() < 3) {
        tokens.push_back(DEFAULT_VASP_POTCAR_DATE_POT_LDA);
      }
      if (pp_type == "GGA" && tokens.size() < 3) {
        tokens.push_back(DEFAULT_VASP_POTCAR_DATE_POT_GGA);
      }
      species_pp_version.push_back(tokens.at(1) + ":" + pp_type + ":" + tokens.at(2));
    } else { // only  one definition
      pp_type = "AE";
      species.push_back(aurostd::VASP_PseudoPotential_CleanName(tokens.at(0)));
      species_Z.push_back(xelement::symbol2Z(aurostd::VASP_PseudoPotential_CleanName(tokens.at(0))));
      species_pp.push_back(aurostd::VASP_PseudoPotential_CleanName(tokens.at(0))); // same as species
      species_pp_type.push_back(pp_type);
      species_pp_version.push_back(tokens.at(0) + ":" + pp_type + ":" + "BigBang"); // ALL Atomic Potential has been existing since the BigBang
    }

    vTITEL[j] = aurostd::RemoveWhiteSpaces(vTITEL[j]);
    vLEXCH.at(j) = aurostd::RemoveWhiteSpaces(vLEXCH.at(j));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vTITEL[j]=" << vTITEL[j] << endl;
      cerr << __AFLOW_FUNC__ << " vLEXCH.at(j)=" << vLEXCH.at(j) << endl;
      cerr << __AFLOW_FUNC__ << " vEATOM.at(j)=" << vEATOM.at(j) << endl;
      cerr << __AFLOW_FUNC__ << " vRMAX.at(j)=" << vRMAX.at(j) << endl;
    }
    xPOTCAR xPOT(xPOTCAR_Finder(species_pp_AUID, species_pp_AUID_collisions, vTITEL[j], vLEXCH.at(j), vEATOM.at(j), vRMAX.at(j), LDEBUG)); // FIXES species_pp_AUID,species_pp_AUID_collisions
    species_pp_groundstate_energy.push_back(xPOT.species_pp_groundstate_energy.at(0));
    species_pp_groundstate_structure.push_back(xPOT.species_pp_groundstate_structure.at(0));

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " SPECIES(" << j << ") [pp, type, version, AUID] = " << species.at(species.size() - 1) << " [" << species_Z.at(species_Z.size() - 1) << ", " << species_pp.at(species_pp.size() - 1)
           << ", " << species_pp_type.at(species_pp_type.size() - 1) << ", " << species_pp_version.at(species_pp_version.size() - 1) << ", " << species_pp_AUID.at(species_pp_AUID.size() - 1) << ", "
           << species_pp_groundstate_energy.at(species_pp_groundstate_energy.size() - 1) << ", " << species_pp_groundstate_structure.at(species_pp_groundstate_structure.size() - 1) << "]" << endl;
    }
  }

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " PSEUDOPOTENTIAL type = " << pp_type << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species.size()=" << species.size() << ": ";
    for (size_t i = 0; i < species.size(); i++) {
      cerr << species[i] << " ";
    }
    cerr << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species_Z.size()=" << species_Z.size() << ": ";
    for (size_t i = 0; i < species_Z.size(); i++) {
      cerr << species_Z[i] << " ";
    }
    cerr << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species_pp.size()=" << species_pp.size() << ": ";
    for (size_t i = 0; i < species_pp.size(); i++) {
      cerr << species_pp[i] << " ";
    }
    cerr << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species_pp_type.size()=" << species_pp_type.size() << ": ";
    for (size_t i = 0; i < species_pp_type.size(); i++) {
      cerr << species_pp_type[i] << " ";
    }
    cerr << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species_pp_version.size()=" << species_pp_version.size() << ": ";
    for (size_t i = 0; i < species_pp_version.size(); i++) {
      cerr << species_pp_version[i] << " ";
    }
    cerr << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species_pp_AUID.size()=" << species_pp_AUID.size() << ": ";
    for (size_t i = 0; i < species_pp_AUID.size(); i++) {
      cerr << species_pp_AUID[i] << " ";
    }
    cerr << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species_pp_AUID_collisions.size()=" << species_pp_AUID_collisions.size() << ": ";
    for (size_t i = 0; i < species_pp_AUID_collisions.size(); i++) {
      cerr << species_pp_AUID_collisions[i] << " ";
    }
    cerr << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species_pp_groundstate_energy.size()=" << species_pp_groundstate_energy.size() << ": ";
    for (size_t i = 0; i < species_pp_groundstate_energy.size(); i++) {
      cerr << species_pp_groundstate_energy[i] << " ";
    }
    cerr << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species_pp_groundstate_structure.size()=" << species_pp_groundstate_structure.size() << ": ";
    for (size_t i = 0; i < species_pp_groundstate_structure.size(); i++) {
      cerr << species_pp_groundstate_structure[i] << " ";
    }
    cerr << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vTITEL.size()=" << vTITEL.size() << ": ";
    for (size_t i = 0; i < vTITEL.size(); i++) {
      cerr << vTITEL[i] << " ";
    }
    cerr << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vLEXCH.size()=" << vLEXCH.size() << ": ";
    for (size_t i = 0; i < vLEXCH.size(); i++) {
      cerr << vLEXCH[i] << " ";
    }
    cerr << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " filename=" << filename << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " AUID=" << AUID << endl;
  }
  if (!species_pp_AUID_collisions.empty()) {
    cerr << __AFLOW_FUNC__ << " COLLISION species_pp_AUID_collisions.size()=" << species_pp_AUID_collisions.size() << endl;
    ERROR_flag = true;
  }
  if (LDEBUG) {
    cerr.flush();
  }

  // ----------------------------------------------------------------------
  // POTCAR_PAW POTCAR_TYPE POTCAR_KINETIC
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD POTCAR_PAW POTCAR_TYPE" << endl;
  }
  POTCAR_PAW = false;
  POTCAR_KINETIC = false;
  POTCAR_GW = false;
  POTCAR_AE = false;
  if (content.find("PAW") != string::npos) {
    POTCAR_PAW = true;
  }
  if (content.find("GW") != string::npos) {
    POTCAR_GW = true;
  }
  if (content.find("AE potential") != string::npos) {
    POTCAR_AE = true;
  }
  // POTCAR DONE **************************************************
  // CHECK FOR US => LDA/GGA
  bool is_US = false;
  bool is_LDA = false;
  POTCAR_TYPE = "";
  for (size_t i = 0; i < vcontent.size() && POTCAR_TYPE.empty(); i++) { // cerr << vcontent[i] << endl;
    if (aurostd::substring2bool(vcontent[i], "TITEL") && aurostd::substring2bool(vcontent[i], "US")) {
      is_US = true;
    }
    if (aurostd::substring2bool(vcontent[i], "GGA ")) {
      is_LDA = true;
    }
    if (aurostd::substring2bool(vcontent[i], "mkinetic energy-density pseudized")) {
      POTCAR_KINETIC = true;
    }
  }
  if (is_US && is_LDA) {
    POTCAR_TYPE = "LDA";
  }
  if (is_US && !is_LDA) {
    POTCAR_TYPE = "GGA";
  }
  for (size_t i = 0; i < vcontent.size() && POTCAR_TYPE.empty(); i++) {
    if (aurostd::substring2bool(vcontent[i], "TITEL")) {
      aurostd::string2tokens(vcontent[i], tokens, "=");
      vcontent[i] = tokens.at(1);
      aurostd::string2tokens(vcontent[i], tokens, " ");
      if (tokens.at(0) == "PAW") {
        POTCAR_TYPE = "PAW_LDA";
      }
      if (tokens.at(0) == "PAW_GGA") {
        POTCAR_TYPE = "PAW_GGA";
      }
      if (tokens.at(0) == "PAW_PBE") {
        POTCAR_TYPE = "PAW_PBE";
      }
      if (tokens.at(0) == "PAW_RPBE") {
        POTCAR_TYPE = "PAW_RPBE";
      }
      if (tokens.at(0) == "PAW_PBE" && POTCAR_KINETIC) {
        POTCAR_TYPE = "PAW_PBE_KIN";
      }
      if (tokens.at(0) == "PAW" && POTCAR_KINETIC) {
        POTCAR_TYPE = "PAW_LDA_KIN";
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " POTCAR_PAW = " << POTCAR_PAW << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " POTCAR_TYPE = " << POTCAR_TYPE << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " POTCAR_KINETIC = " << POTCAR_KINETIC << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " POTCAR_GW = " << POTCAR_GW << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " POTCAR_AE = " << POTCAR_AE << endl;
  }

  // ----------------------------------------------------------------------
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " WRITE vxpseudopotentials" << endl;
  }

  if (XHOST.PSEUDOPOTENTIAL_GENERATOR && species.size() == 1) {
    xPOTCAR_PURE_Printer((*this), cout);
  }

  // ----------------------------------------------------------------------
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " (BULLSHIT DONT USE) title=" << title << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vcontent.size()=" << vcontent.size() << endl;
  }

  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " END (" << time_delay(seconds) << ")" << endl;
  }
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  if (ERROR_flag) {
    message << "ERROR_flag set in xPOTCAR";
    if (force_exit) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    } else {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_, QUIET);
    }
    return false;
  }
  return true;
}

// JSON serialization of xPOTCAR
aurostd::JSON::object xPOTCAR::serialize() const {
  return aurostd::JSON::object({AST_JSON_GETTER(JSON_xPOTCAR_MEMBERS)});
}

xPOTCAR xPOTCAR::deserialize(const aurostd::JSON::object& jo) {
  AST_JSON_SETTER(JSON_xPOTCAR_MEMBERS)
  return *this;
}

//---------------------------------------------------------------------------------
// class xVASPRUNXML
//---------------------------------------------------------------------------------
bool xVASPRUNXML::GetProperties(const string& stringIN, bool QUIET) {
  stringstream sss;
  sss.str(stringIN);
  if (filename.empty()) {
    filename = "string";
  }
  return xVASPRUNXML::GetProperties(sss, QUIET);
}

bool xVASPRUNXML::GetPropertiesFile(const string& fileIN, bool QUIET) {
  stringstream sss;
  if (filename.empty()) {
    filename = fileIN;
  }
  aurostd::compressfile2stringstream(fileIN, sss);
  return xVASPRUNXML::GetProperties(sss, QUIET);
}

bool xVASPRUNXML::GetPropertiesUrlFile(const string& url, const string& file, bool QUIET) {
  const string tmpfile = aurostd::TmpFileCreate("xVASPRUNXML_GetProperties"); // CO20200502 - threadID
  aurostd::httpGetFileStatus(url + "/" + file, tmpfile);
  const bool out = GetPropertiesFile(tmpfile, QUIET);
  filename = "url=" + url; // CO20210315
  aurostd::RemoveFile(tmpfile);
  return out;
}

bool xVASPRUNXML::GetProperties(const stringstream& stringstreamIN, bool QUIET) {
  const bool LDEBUG = (false || XHOST.DEBUG || !QUIET);
  stringstream message;
  const bool force_exit = XHOST.POSTPROCESS; // SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only
  bool ERROR_flag = false;
  const long double seconds = aurostd::get_seconds();
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BEGIN (" << time_delay(seconds) << ")" << endl;
  }
  clear(); // so it does not mess up vector/deque
  content = stringstreamIN.str();
  vcontent.clear();
  vector<string> tokens;
  aurostd::string2vectorstring(content, vcontent);
  string line;
  if (filename.empty()) {
    filename = "stringstream";
  }
  // crunching to eat the info

  // ----------------------------------------------------------------------
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vcontent.size()=" << vcontent.size() << endl;
  }

  // ----------------------------------------------------------------------
  // LOAD natoms

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD natoms DATA" << endl;
  }
  line = "";
  for (int iline = (int) vcontent.size() - 1; iline >= 0; iline--) { // NEW FROM THE BACK
    if (aurostd::substring2bool(vcontent[iline], "<atoms>")) {
      if (aurostd::substring2bool(vcontent[iline], "</atoms>")) {
        line = vcontent[iline];
        break;
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " line=" << line << endl;
  }
  aurostd::StringSubstInPlace(line, "<atoms>", "");
  aurostd::StringSubstInPlace(line, "</atoms>", "");

  // ----------------------------------------------------------------------
  // LOAD NATOMS
  natoms = 0.0;
  aurostd::string2tokens(line, tokens, ">");
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " tokens.size()=" << tokens.size() << endl;
  }
  natoms = aurostd::string2utype<double>(line);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " natoms=" << natoms << endl;
  }

  // ----------------------------------------------------------------------
  // LOAD FORCES
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD FORCES DATA" << endl;
  }
  vforces.clear(); // QM FORCES calculation
  for (int iline = (int) vcontent.size() - 1; iline >= 0; iline--) {
    if (aurostd::substring2bool(vcontent[iline], "<varray name=\"forces\" >")) {
      for (size_t iat = 0; iat < (size_t) natoms && iat < vcontent.size(); iat++) {
        aurostd::StringSubstInPlace(vcontent[iline + iat + 1], "<v>", "");
        aurostd::StringSubstInPlace(vcontent[iline + iat + 1], "</v>", "");
        // if(LDEBUG) cerr << __AFLOW_FUNC__ << " vcontent[iline+iat+1]=" << vcontent[iline+iat+1] << endl;
        aurostd::string2tokens(vcontent[iline + iat + 1], tokens, " ");
        if (tokens.size() == 3) {
          xvector<double> force(3);
          force[1] = aurostd::string2utype<double>(tokens.at(0));
          force[2] = aurostd::string2utype<double>(tokens.at(1));
          force[3] = aurostd::string2utype<double>(tokens.at(2));
          vforces.push_back(force); // cerr.precision(20);
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " force=" << force << endl;
          }
          // if(LDEBUG) cerr << __AFLOW_FUNC__ << " force=" << force[1]  << " " << force[2] << " " << force[3] << endl;
        } else {
          message << "Error in QM FORCES calculation" << endl;
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
          ERROR_flag = true;
        }
      }
      iline = -1;
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vforces.size()=" << vforces.size() << endl;
  }

  // ----------------------------------------------------------------------
  // LOAD KPOINTLIST
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD KPOINTLIST DATA" << endl;
  }
  vkpoint.clear(); // QM KPOINTLIST calculation
  for (int iline = 0; iline < (int) vcontent.size(); iline++) {
    if (aurostd::substring2bool(vcontent[iline], "<varray name=\"kpointlist\" >")) {
      iline++;
      for (int iat = 0; iline + iat < (int) vcontent.size(); iat++) {
        if (aurostd::substring2bool(vcontent[iline + iat], "</varray>")) {
          iline = (int) vcontent.size();
        } else {
          // cerr << __AFLOW_FUNC__ << " vcontent[iline+iat]=" << vcontent[iline+iat] << endl;
          aurostd::StringSubstInPlace(vcontent[iline + iat], "<v>", "");
          aurostd::StringSubstInPlace(vcontent[iline + iat], "</v>", "");
          aurostd::string2tokens(vcontent[iline + iat], tokens, " ");
          if (tokens.size() == 3) {
            xvector<double> kpoint(3);
            kpoint[1] = aurostd::string2utype<double>(tokens.at(0));
            kpoint[2] = aurostd::string2utype<double>(tokens.at(1));
            kpoint[3] = aurostd::string2utype<double>(tokens.at(2));
            vkpoint.push_back(kpoint); // cerr.precision(20);
            //	    if(LDEBUG) cerr << __AFLOW_FUNC__ << " kpoint=" << kpoint << endl;
          } else {
            message << "Error in QM KPOINTLIST calculation" << endl;
            pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
            ERROR_flag = true;
          }
        }
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vkpoint.size()=" << vkpoint.size() << endl;
  }

  // ----------------------------------------------------------------------
  // LOAD WEIGHTS
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD WEIGHTS DATA" << endl;
  }
  vweights.clear(); // QM WEIGHTS calculation
  for (int iline = 0; iline < static_cast<int>(vcontent.size()); iline++) {
    if (aurostd::substring2bool(vcontent[iline], "<varray name=\"weights\" >")) {
      iline++;
      for (int iat = 0; iline + iat < (int) vcontent.size(); iat++) {
        if (aurostd::substring2bool(vcontent[iline + iat], "</varray>")) {
          iline = (int) vcontent.size();
        } else {
          // cerr << __AFLOW_FUNC__ << " vcontent[iline+iat]=" << vcontent[iline+iat] << endl;
          aurostd::StringSubstInPlace(vcontent[iline + iat], "<v>", "");
          aurostd::StringSubstInPlace(vcontent[iline + iat], "</v>", "");
          aurostd::string2tokens(vcontent[iline + iat], tokens, " ");
          if (tokens.size() == 1) {
            vweights.emplace_back(aurostd::string2utype<double>(tokens.at(0))); // cerr.precision(20);
            //	    if(LDEBUG) cerr << __AFLOW_FUNC__ << " weight=" << weight << endl;
          } else {
            message << "Error in QM WEIGHTS calculation" << endl;
            pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
            ERROR_flag = true;
          }
        }
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vweights.size()=" << vweights.size() << endl;
  }

  // ----------------------------------------------------------------------
  // LOAD STRESS
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD STRESS DATA" << endl;
  }
  stress.clear(); // QM STRESS calculation
  for (int iline = (int) vcontent.size() - 1; iline >= 0; iline--) {
    if (aurostd::substring2bool(vcontent[iline], "<varray name=\"stress\" >")) {
      for (size_t iat = 0; iat < (size_t) 3 && iat < vcontent.size(); iat++) { // only three lines
        aurostd::StringSubstInPlace(vcontent[iline + iat + 1], "<v>", "");
        aurostd::StringSubstInPlace(vcontent[iline + iat + 1], "</v>", "");
        // if(LDEBUG) cerr << __AFLOW_FUNC__ << " vcontent[iline+iat+1]=" << vcontent[iline+iat+1] << endl;
        aurostd::string2tokens(vcontent[iline + iat + 1], tokens, " ");
        if (tokens.size() == 3) {
          stress(iat + 1, 1) = aurostd::string2utype<double>(tokens.at(0));
          stress(iat + 1, 2) = aurostd::string2utype<double>(tokens.at(1));
          stress(iat + 1, 3) = aurostd::string2utype<double>(tokens.at(2));
        } else {
          message << "Error in QM STRESS calculation" << endl;
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
          ERROR_flag = true;
        }
      }
      iline = -1;
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " stress=" << endl << stress << endl;
  }

  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " END (" << time_delay(seconds) << ")" << endl;
  }
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  if (ERROR_flag) {
    message << "ERROR_flag set in xVASPRUNXML";
    if (force_exit) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    } else {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_, QUIET);
    }
    return false;
  }
  return true;
}

// ME20190204 BEGIN
//  GetForces just grabs the forces from the vasprun file with a much faster
//  algorithm than in GetProperties. This is significant for APL and AAPL which
//  process many large vasprun files.
bool xVASPRUNXML::GetForces(const string& stringIN, bool QUIET) {
  stringstream sss;
  sss.str(stringIN);
  if (filename.empty()) {
    filename = "string";
  }
  return xVASPRUNXML::GetForces(sss, QUIET);
}

bool xVASPRUNXML::GetForcesFile(const string& fileIN, bool QUIET) {
  stringstream sss;
  if (filename.empty()) {
    filename = fileIN;
  }
  aurostd::compressfile2stringstream(fileIN, sss);
  return xVASPRUNXML::GetForces(sss, QUIET);
}

// Cannot use const stringstream& because of std::getline
bool xVASPRUNXML::GetForces(stringstream& stringstreamIN, bool QUIET) {
  const bool LDEBUG = (false || XHOST.DEBUG || !QUIET);
  stringstream message;
  const bool force_exit = XHOST.POSTPROCESS; // SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only
  bool ERROR_flag = false;

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
  }

  vforces.clear();
  string line;
  vector<double> tokens;
  xvector<double> force(3);
  while (std::getline(stringstreamIN, line)) {
    if (aurostd::substring2bool(line, "<varray name=\"forces\" >")) {
      while (std::getline(stringstreamIN, line) && !aurostd::substring2bool(line, "</varray>")) {
        line = aurostd::RemoveSubStringFirst(line, "<v>");
        line = aurostd::RemoveSubString(line, "</v>");
        aurostd::string2tokens(line, tokens, " ");
        if (tokens.size() == 3) {
          force[1] = tokens[0];
          force[2] = tokens[1];
          force[3] = tokens[2];
          vforces.push_back(force);
        } else {
          message << "Error in QM FORCES calculation" << endl;
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
          ERROR_flag = true;
        }
      }
      return true;
    }
  }
  ERROR_flag = true; // set true if you get here
  if (ERROR_flag) {
    message << "ERROR_flag set: forces not found";
    if (force_exit) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    } else {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_, QUIET);
    }
    return false;
  }
  return false; // catch all false
}

// JSON serialization of xVASPRUNXML
aurostd::JSON::object xVASPRUNXML::serialize() const {
  return aurostd::JSON::object({AST_JSON_GETTER(JSON_xVASPRUNXML_MEMBERS)});
}

xVASPRUNXML xVASPRUNXML::deserialize(const aurostd::JSON::object& jo) {
  AST_JSON_SETTER(JSON_xVASPRUNXML_MEMBERS)
  return *this;
}

// ME20190204 END

//---------------------------------------------------------------------------------
// class xIBZKPT
//---------------------------------------------------------------------------------
bool xIBZKPT::GetProperties(const string& stringIN, bool QUIET) {
  stringstream sss;
  sss.str(stringIN);
  if (filename.empty()) {
    filename = "string";
  }
  return xIBZKPT::GetProperties(sss, QUIET);
}

bool xIBZKPT::GetPropertiesFile(const string& fileIN, bool QUIET) {
  stringstream sss;
  if (filename.empty()) {
    filename = fileIN;
  }
  aurostd::compressfile2stringstream(fileIN, sss);
  return xIBZKPT::GetProperties(sss, QUIET);
}

bool xIBZKPT::GetPropertiesUrlFile(const string& url, const string& file, bool QUIET) {
  const string tmpfile = aurostd::TmpFileCreate("xIBZKPT_GetProperties"); // CO20200502 - threadID
  aurostd::httpGetFileStatus(url + "/" + file, tmpfile);
  const bool out = GetPropertiesFile(tmpfile, QUIET);
  filename = "url=" + url; // CO20210315
  aurostd::RemoveFile(tmpfile);
  return out;
}

bool xIBZKPT::GetProperties(const stringstream& stringstreamIN, bool QUIET) {
  const bool LDEBUG = (false || XHOST.DEBUG || !QUIET);
  stringstream message;
  const bool force_exit = XHOST.POSTPROCESS; // SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only
  bool ERROR_flag = false;

  const long double seconds = aurostd::get_seconds();
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BEGIN (" << time_delay(seconds) << ")" << endl;
  }
  clear(); // so it does not mess up vector/deque
  content = stringstreamIN.str();
  vcontent.clear();
  vector<string> tokens;
  aurostd::string2vectorstring(content, vcontent);
  if (filename.empty()) {
    filename = "stringstream";
  }
  // crunching to eat the info

  // ----------------------------------------------------------------------
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vcontent.size()=" << vcontent.size() << endl;
  }

  // ----------------------------------------------------------------------
  // LOAD NWEIGHTS VKPOINT
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD VKPOINT DATA" << endl;
  }
  nweights = 0;
  nkpoints_irreducible = 0;
  vkpoint.clear(); // QM VKPOINT calculation
  vweights.clear(); // QM VKPOINT calculation
  for (int iline = 0; iline < (int) vcontent.size(); iline++) {
    if (aurostd::substring2bool(vcontent[iline], "Automatically generated mesh")) {
      iline++;
      nkpoints_irreducible = aurostd::string2utype<double>(vcontent[iline]);
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " nkpoints_irreducible=" << nkpoints_irreducible << endl;
      }
      iline += 2; // skip text
      //  cerr << __AFLOW_FUNC__ << " vcontent[iline]=" << vcontent[iline] << endl;
      for (uint iat = 0; iat < nkpoints_irreducible; iat++) {
        // 	cerr << __AFLOW_FUNC__ << " vcontent[iline+iat]=" << vcontent[iline+iat] << endl;
        aurostd::string2tokens(vcontent[iline + iat], tokens, " ");
        if (tokens.size() == 4) {
          xvector<double> kpoint(3);
          kpoint[1] = aurostd::string2utype<double>(tokens.at(0));
          kpoint[2] = aurostd::string2utype<double>(tokens.at(1));
          kpoint[3] = aurostd::string2utype<double>(tokens.at(2));
          vkpoint.push_back(kpoint); // cerr.precision(20);
          vweights.push_back(aurostd::string2utype<uint>(tokens.at(3)));
          nweights += aurostd::string2utype<uint>(tokens.at(3));
          //	  if(LDEBUG) cerr << __AFLOW_FUNC__ << " kpoint=" << kpoint << " " << "weight=" << aurostd::string2utype<double>(tokens.at(3))<< endl;
        } else {
          message << "Error in QM NWEIGHTS/VKPOINT calculation" << endl;
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_, QUIET);
          ERROR_flag = true;
        }
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vkpoint.size()=" << vkpoint.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vweights.size()=" << vweights.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " nweights=" << nweights << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " nkpoints_irreducible=" << nkpoints_irreducible << endl;
  }

  // ---------------------------------------------------------------------
  // LOAD NTETRAHEDRA WTETRAHEDRA TETRAHEDRA
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD TETRAHEDRA DATA" << endl;
  }
  ntetrahedra = 0;
  wtetrahedra = 0.0;
  vtetrahedra.clear(); // QM TETRAHEDRA calculation
  for (int iline = 0; iline < (int) vcontent.size(); iline++) {
    if (aurostd::substring2bool(vcontent[iline], "Tetrahedra")) {
      iline++;
      aurostd::string2tokens(vcontent[iline], tokens, " ");
      ntetrahedra = aurostd::string2utype<uint>(tokens.at(0));
      wtetrahedra = aurostd::string2utype<double>(tokens.at(1));
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " ntetrahedra=" << ntetrahedra << endl;
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " wtetrahedra=" << wtetrahedra << endl;
      }
      iline += 1; // skip text
      for (uint iat = 0; iat < ntetrahedra; iat++) {
        // 	cerr << __AFLOW_FUNC__ << " vcontent[iline+iat]=" << vcontent[iline+iat] << endl;
        aurostd::string2tokens(vcontent[iline + iat], tokens, " ");
        if (tokens.size() == 5) {
          xvector<int> tetrahedra(5);
          tetrahedra[1] = aurostd::string2utype<double>(tokens.at(0));
          tetrahedra[2] = aurostd::string2utype<double>(tokens.at(1));
          tetrahedra[3] = aurostd::string2utype<double>(tokens.at(2));
          tetrahedra[4] = aurostd::string2utype<double>(tokens.at(3));
          tetrahedra[5] = aurostd::string2utype<double>(tokens.at(4));
          vtetrahedra.push_back(tetrahedra); // cerr.precision(20);
          //	  if(LDEBUG) cerr << __AFLOW_FUNC__ << " tetrahedra=" << tetrahedra << " " << endl;
        } else {
          if (!QUIET) {
            cerr << __AFLOW_FUNC__ << " error in QM NTETRAHEDRA/WTETRAHEDRA/TETRAHEDRA calculation" << endl;
          }
          ERROR_flag = true;
        }
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ntetrahedra=" << ntetrahedra << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " wtetrahedra=" << wtetrahedra << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vtetrahedra.size()=" << vtetrahedra.size() << endl;
  }

  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " END (" << time_delay(seconds) << ")" << endl;
  }
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  if (ERROR_flag) {
    message << "ERROR_flag set in xIBZKPT";
    if (force_exit) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    } else {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_, QUIET);
    }
    return false;
  }
  return true;
}

// JSON serialization of xIBZKPT
aurostd::JSON::object xIBZKPT::serialize() const {
  return aurostd::JSON::object({AST_JSON_GETTER(JSON_xIBZKPT_MEMBERS)});
}

xIBZKPT xIBZKPT::deserialize(const aurostd::JSON::object& jo) {
  AST_JSON_SETTER(JSON_xIBZKPT_MEMBERS)
  return *this;
}

//---------------------------------------------------------------------------------
// class xKPOINTS
//---------------------------------------------------------------------------------
bool xKPOINTS::GetProperties(const string& stringIN, bool QUIET) {
  stringstream sss;
  sss.str(stringIN);
  if (filename.empty()) {
    filename = "string";
  }
  return xKPOINTS::GetProperties(sss, QUIET);
}

bool xKPOINTS::GetPropertiesFile(const string& fileIN, bool QUIET) {
  stringstream sss;
  if (filename.empty()) {
    filename = fileIN;
  }
  aurostd::compressfile2stringstream(fileIN, sss);
  return xKPOINTS::GetProperties(sss, QUIET);
}

bool xKPOINTS::GetPropertiesUrlFile(const string& url, const string& file, bool QUIET) {
  const string tmpfile = aurostd::TmpFileCreate("xKPOINTS_GetProperties"); // CO20200502 - threadID
  aurostd::httpGetFileStatus(url + "/" + file, tmpfile);
  const bool out = GetPropertiesFile(tmpfile, QUIET);
  filename = "url=" + url; // CO20210315
  aurostd::RemoveFile(tmpfile);
  return out;
}

bool xKPOINTS::GetProperties(const stringstream& stringstreamIN, bool QUIET) {
  const bool LDEBUG = (false || XHOST.DEBUG || !QUIET);
  stringstream message;
  const bool force_exit = XHOST.POSTPROCESS; // SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only

  const bool ERROR_flag = false;
  clear(); // so it does not mess up vector/deque

  const long double seconds = aurostd::get_seconds();
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BEGIN (" << time_delay(seconds) << ")" << endl;
  }

  content = stringstreamIN.str();
  vcontent.clear();
  vector<string> tokens;
  aurostd::string2vectorstring(content, vcontent);
  if (filename.empty()) {
    filename = "stringstream";
  }
  // crunching to eat the info
  title = "";
  mode = -1;
  grid_type = "";
  is_KPOINTS_NNN = false;
  is_KPOINTS_PATH = false;
  nnn_kpoints.clear(); // N*N*N
  ooo_kpoints.clear(); // ORIGIN
  nkpoints = 0;
  path_mode = "";
  path = "";
  path_grid = 0;
  vpath.clear();
  // ----------------------------------------------------------------------
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vcontent.size()=" << vcontent.size() << endl;
  }

  // ----------------------------------------------------------------------
  // CHECK IF WITH KPOINTS NUMBERS
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " CHECK IF WITH KPOINTS NUMBERS" << endl;
  }
  if (!is_KPOINTS_NNN && !is_KPOINTS_PATH && vcontent.size() >= 5) {
    if (vcontent.at(2).at(0) == 'M' || vcontent.at(2).at(0) == 'm' || vcontent.at(2).at(0) == 'G' || vcontent.at(2).at(0) == 'g') {
      aurostd::string2tokens(vcontent.at(3), tokens);
      //    if(LDEBUG) cerr << __AFLOW_FUNC__ << " tokens.size()=" << tokens.size() << endl;
      if (tokens.size() == 3) {
        aurostd::string2tokens(vcontent.at(4), tokens);
        //     if(LDEBUG) cerr << __AFLOW_FUNC__ << " tokens.size()=" << tokens.size() << endl;
        if (tokens.size() == 3) {
          is_KPOINTS_NNN = true;
          is_KPOINTS_PATH = false;
          title = vcontent.at(0);
          mode = aurostd::string2utype<int>(vcontent.at(1));
          grid_type = vcontent.at(2);
          aurostd::string2tokens(vcontent.at(3), tokens);
          nnn_kpoints[1] = aurostd::string2utype<int>(tokens.at(0));
          nnn_kpoints[2] = aurostd::string2utype<int>(tokens.at(1));
          nnn_kpoints[3] = aurostd::string2utype<int>(tokens.at(2));
          nkpoints = nnn_kpoints[1] * nnn_kpoints[2] * nnn_kpoints[3];
          aurostd::string2tokens(vcontent.at(4), tokens);
          ooo_kpoints[1] = aurostd::string2utype<double>(tokens.at(0));
          ooo_kpoints[2] = aurostd::string2utype<double>(tokens.at(1));
          ooo_kpoints[3] = aurostd::string2utype<double>(tokens.at(2));
        }
      }
    }
  }

  // ----------------------------------------------------------------------
  // CHECK IF WITH PATH
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " CHECK IF WITH PATH" << endl;
  }
  if (!is_KPOINTS_NNN && !is_KPOINTS_PATH && vcontent.size() >= 5) {
    if (vcontent.at(2).at(0) == 'L' || vcontent.at(2).at(0) == 'l') {
      is_KPOINTS_NNN = false;
      is_KPOINTS_PATH = true;
      title = vcontent.at(0);
      mode = aurostd::string2utype<int>(vcontent.at(1));
      path_grid = aurostd::string2utype<int>(vcontent.at(1));
      grid_type = vcontent.at(2);
      path_mode = vcontent.at(3);
      xvector<double> kpt(3); // ME20190614
      for (size_t iline = 4; iline < vcontent.size(); iline++) {
        aurostd::StringSubstInPlace(vcontent[iline], "!", "@");
        aurostd::StringSubstInPlace(vcontent[iline], "@", "@ "); // CO20210712 - for "!K" vs. "! K"
        if (aurostd::substring2bool(vcontent[iline], "@")) { // avoid removing ! as comment
          aurostd::string2tokens(vcontent[iline], tokens, " ");
          if (tokens.size() >= 5) {
            //	    if(LDEBUG) cerr << __AFLOW_FUNC__ << " tokens.size()=" << tokens.size() << endl;
            for (size_t k = 0; k < tokens.size(); k++) {
              if (tokens[k] == "@" && k + 1 < tokens.size()) {
                vpath.push_back(tokens.at(k + 1));
              } else if (k < 3) { // ME20190614
                kpt[k + 1] = aurostd::string2utype<double>(tokens[k]);
              }
            }
            vkpoints.push_back(kpt); // ME20190614
          }
        }
      }
      // ----------------------------------------------------------------------
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " vpath.size()=" << vpath.size() << endl;
      }
      path = "";
      for (size_t i = 0; i < vpath.size(); i += 2) {
        path += vpath.at(i) + "-" + vpath.at(i + 1) + (i + 2 < vpath.size() ? "," : "");
      }
      // \Gamma-X,X-W,W-K,K-\Gamma,\Gamma-L,L-U,U-W,W-L,L-K,U-X
      // \Gamma-X-W-K-\Gamma-L-U-W-L-K,U-X
      // \Gamma-X-W-K-\Gamma-L-U-W-L-K,U-X
      //	if(LDEBUG) cerr << __AFLOW_FUNC__ << " path=[" << path << "]" << endl;
    }
  }
  // ----------------------------------------------------------------------

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " title=[" << title << "]" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " mode=" << mode << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " grid_type=[" << grid_type << "]" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " nkpoints=" << nkpoints << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " nnn_kpoints=" << nnn_kpoints << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " ooo_kpoints=" << ooo_kpoints << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " is_KPOINTS_NNN=" << is_KPOINTS_NNN << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " is_KPOINTS_PATH=" << is_KPOINTS_PATH << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " path_mode=[" << path_mode << "]" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " path=[" << path << "]" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vpath.size()=" << vpath.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vpath=";
    for (size_t i = 0; i < vpath.size(); i++) {
      cerr << vpath[i] << " ";
    }
    cerr << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " path_grid=" << path_grid << endl;
  }

  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " END (" << time_delay(seconds) << ")" << endl;
  }
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  if (ERROR_flag) {
    message << "ERROR_flag set in xKPOINTS";
    if (force_exit) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    } else {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_, QUIET);
    }
    return false;
  }
  return true;
}

// ME20190623 BEGIN
ostream& operator<<(ostream& oss, const xKPOINTS& xkpts) {
  if (xkpts.is_KPOINTS_NNN) {
    oss << xkpts.title << std::endl;
    oss << xkpts.mode << std::endl;
    oss << xkpts.grid_type << std::endl;
    oss << aurostd::joinWDelimiter(xkpts.nnn_kpoints, " ") << std::endl;
    if (aurostd::iszero(xkpts.ooo_kpoints)) {
      oss << "0 0 0" << std::endl; // No need to print all those decimal places if zero
    } else {
      oss << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
      for (int i = 1; i < 4; i++) {
        oss << std::scientific << std::setprecision(4) << std::setw(9) << xkpts.ooo_kpoints[i];
      }
      oss << std::endl;
    }
  } else if (xkpts.is_KPOINTS_PATH) {
    oss << xkpts.title << std::endl;
    oss << xkpts.path_grid << std::endl;
    oss << xkpts.grid_type << std::endl;
    oss << xkpts.path_mode << std::endl;
    for (size_t i = 0; i < xkpts.vpath.size(); i++) {
      for (int j = 1; j < 4; j++) {
        oss << std::setprecision(4) << std::fixed << std::setw(9) << xkpts.vkpoints[i][j];
      }
      oss << "  ! " << xkpts.vpath[i] << std::endl;
      if ((i % 2 == 1) && (i < xkpts.vpath.size() - 1)) {
        oss << std::endl;
      }
    }
  }
  return oss;
}

string xKPOINTS::createStandardTitlePath(const xstructure& xstr) {
  // Lattice part
  string stdtitle = xstr.bravais_lattice_type;
  if (stdtitle == "CUB") {
    stdtitle += " (simple cubic) ";
  } else if (stdtitle == "FCC") {
    stdtitle += " (face-centered cubic) ";
  } else if (stdtitle == "BCC") {
    stdtitle += " (body-centered cubic) ";
  } else if (stdtitle == "TET") {
    stdtitle += " (tetragonal) ";
  } else if (stdtitle == "BCT") {
    stdtitle += " (body-centered tetragonal) ";
  } else if (stdtitle == "BCT1") {
    stdtitle += " (body-centered tetragonal c < a) ";
  } else if (stdtitle == "BCT2") {
    stdtitle += " (body-centered tetragonal a < c) ";
  } else if (stdtitle == "ORC") {
    stdtitle += " (orthorhombic) ";
  } else if (stdtitle == "ORCF") {
    stdtitle += " (face-centered orthorhombic) ";
  } else if (stdtitle == "ORCF1") {
    stdtitle += " (face-centered orthorhombic 1/a^2 > 1/b^2+1/c^2) ";
  } else if (stdtitle == "ORCF2") {
    stdtitle += " (face-centered orthorhombic 1/a^2 < 1/b^2+1/c^2) ";
  } else if (stdtitle == "ORCF3") {
    stdtitle += " (face-centered orthorhombic 1/a^2 = 1/b^2+1/c^2) ";
  } else if (stdtitle == "ORCI") {
    stdtitle += " (body-centered orthorhombic a < b < c) ";
  } else if (stdtitle == "ORCC") {
    stdtitle += " (C-centered orthorhombic a < b) ";
  } else if (stdtitle == "HEX") {
    stdtitle += " (hexagonal) ";
  } else if (stdtitle == "RHL") {
    stdtitle += " (rhombohedral) ";
  } else if (stdtitle == "RHL1") {
    stdtitle += " (rhombohedral alpha < 90) ";
  } else if (stdtitle == "RHL2") {
    stdtitle += " (rhombohedral alpha > 90) ";
  } else if (stdtitle == "MCL") {
    stdtitle += " (monoclinic) ";
  } else if (stdtitle == "MCLC") {
    stdtitle += " (C-centered monoclinic) ";
  } else if (stdtitle == "MCLC1") {
    stdtitle += " (C-centered monoclinic kgamma > 90) ";
  } else if (stdtitle == "MCLC2") {
    stdtitle += " (C-centered monoclinic kgamma = 90) ";
  } else if (stdtitle == "MCLC3") {
    stdtitle += " (C-centered monoclinic kgamma < 90, bcos(alpha)/c+(bsin(alpha)/a)^2<1) ";
  } else if (stdtitle == "MCLC4") {
    stdtitle += " (C-centered monoclinic kgamma < 90, bcos(alpha)/c+(bsin(alpha)/a)^2=1) ";
  } else if (stdtitle == "MCLC5") {
    stdtitle += " (C-centered monoclinic kgamma < 90, bcos(alpha)/c+(bsin(alpha)/a)^2>1) ";
  } else if ((stdtitle == "TRI") || (aurostd::toupper(stdtitle) == "TRI1A") || (aurostd::toupper(stdtitle) == "TRI1B") || (aurostd::toupper(stdtitle) == "TRI2A") || (aurostd::toupper(stdtitle) == "TRI2B")) {
    stdtitle += " (triclinic) ";
  }

  // Path part
  stdtitle += vpath[0];
  for (size_t i = 2; i < vpath.size(); i += 2) {
    stdtitle += "-" + vpath[i - 1];
    if (vpath[i - 1] != vpath[i]) {
      stdtitle += "|" + vpath[i];
    }
  }
  stdtitle += "-" + vpath.back();
  return stdtitle;
}

// ME20190623 END

// JSON serialization of xKPOINTS
aurostd::JSON::object xKPOINTS::serialize() const {
  return aurostd::JSON::object({AST_JSON_GETTER(JSON_xKPOINTS_MEMBERS)});
}

xKPOINTS xKPOINTS::deserialize(const aurostd::JSON::object& jo) {
  AST_JSON_SETTER(JSON_xKPOINTS_MEMBERS)
  return *this;
}

//---------------------------------------------------------------------------------
// class xCHGCAR
//---------------------------------------------------------------------------------
bool xCHGCAR::GetProperties(const string& stringIN, bool QUIET) {
  stringstream sss;
  sss.str(stringIN);
  if (filename.empty()) {
    filename = "string";
  }
  return xCHGCAR::GetProperties(sss, QUIET);
}

bool xCHGCAR::GetPropertiesFile(const string& fileIN, bool QUIET) {
  stringstream sss;
  if (filename.empty()) {
    filename = fileIN;
  }
  aurostd::compressfile2stringstream(fileIN, sss);
  return xCHGCAR::GetProperties(sss, QUIET);
}

bool xCHGCAR::GetPropertiesUrlFile(const string& url, const string& file, bool QUIET) {
  const string tmpfile = aurostd::TmpFileCreate("xCHGCAR_GetProperties"); // CO20200502 - threadID
  aurostd::httpGetFileStatus(url + "/" + file, tmpfile);
  const bool out = GetPropertiesFile(tmpfile, QUIET);
  filename = "url=" + url; // CO20210315
  aurostd::RemoveFile(tmpfile);
  return out;
}

bool xCHGCAR::GetProperties(const stringstream& stringstreamIN, bool QUIET) {
  const bool LDEBUG = (false || XHOST.DEBUG || !QUIET);
  stringstream message;
  const bool force_exit = XHOST.POSTPROCESS; // SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only

  const bool ERROR_flag = false;
  const long double seconds = aurostd::get_seconds();
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BEGIN (" << time_delay(seconds) << ")" << endl;
  }
  clear(); // so it does not mess up vector/deque
  content = stringstreamIN.str();
  vcontent.clear();
  vector<string> tokens;
  aurostd::string2vectorstring(content, vcontent);
  if (filename.empty()) {
    filename = "stringstream";
  }

  // crunching to eat the info
  uint natoms = 0;
  grid.clear(); // N*N*N
  vstring.clear();
  vvalues.clear();
  tvalues.clear();
  uint index = 5;
  // ----------------------------------------------------------------------
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vcontent.size()=" << vcontent.size() << endl;
  }
  // ----------------------------------------------------------------------
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD GRID " << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vcontent.at(" << index << ")=" << vcontent.at(index) << endl;
  }
  aurostd::string2tokens(vcontent.at(index), tokens);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " tokens.size()=" << tokens.size() << endl;
  }
  for (size_t i = 0; i < tokens.size(); i++) {
    natoms += aurostd::string2utype<uint>(tokens[i]);
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " natoms=" << natoms << endl;
  }
  index += natoms + 1 + 2; // skip direct and atoms space and get grid
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vcontent.at(" << index << ")=" << vcontent.at(index) << endl;
  }
  aurostd::string2tokens(vcontent.at(index), tokens);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " tokens.size()=" << tokens.size() << endl;
  }
  for (size_t i = 0; i < tokens.size(); i++) {
    grid(i + 1) = aurostd::string2utype<uint>(tokens[i]);
  }
  index++;
  const uint size_grid = grid(1) * grid(2) * grid(3);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " grid=[" << grid(1) << "," << grid(2) << "," << grid(3) << "]" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " size_grid=" << size_grid << "]" << endl;
  }
  // ----------------------------------------------------------------------
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD VSTRING " << endl;
  }
  for (size_t i = index; i < vcontent.size(); i++) {
    aurostd::string2tokens(vcontent.at(i), tokens);
    //    if(LDEBUG) cerr << __AFLOW_FUNC__ << " tokens.size()=" << tokens.size() << endl;
    for (size_t j = 0; j < tokens.size(); j++) {
      if (vstring.size() < size_grid) {
        vstring.push_back(tokens[j]);
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vstring.size()=" << vstring.size() << endl;
  }
  // ----------------------------------------------------------------------
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD VVALUES " << endl;
  }
  xvector<double> vvalues_aus(vstring.size());
  for (size_t i = 0; i < vstring.size(); i++) {
    vvalues_aus(i + 1) = aurostd::string2utype<double>(vstring[i]);
  }
  // now copy on the real vvalues which has undefined size
  vvalues = vvalues_aus;
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vvalues.rows=" << vvalues.rows << endl;
  }
  // ----------------------------------------------------------------------
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " LOAD TVALUES " << endl;
  }
  const aurostd::xtensor<double> tvalues_aus(grid); // ME20180705
  int iii = 0;
  for (int i3 = 1; i3 <= grid(3); i3++) { // CO - x is fastest, z is slowest
    for (int i2 = 1; i2 <= grid(2); i2++) {
      for (int i1 = 1; i1 <= grid(1); i1++) { // CO - x is fastest, z is slowest
        tvalues_aus[i1][i2][i3] = vvalues(++iii); // ME20180705  //CO20200404 - using ++iii
      }
    }
  }

  tvalues = tvalues_aus;
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " tvalues.shape[1]=" << tvalues.shape[1] << endl; // ME20180705
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " tvalues.shape[2]=" << tvalues.shape[2] << endl; // ME20180705
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " tvalues.shape[3]=" << tvalues.shape[3] << endl; // ME20180705
  }
  // ----------------------------------------------------------------------
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " END (" << time_delay(seconds) << ")" << endl;
  }
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  if (ERROR_flag) {
    message << "ERROR_flag set in xCHGCAR";
    if (force_exit) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    } else {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_, QUIET);
    }
    return false;
  }
  return true;
}

//-------------------------------------------------------------------------------------------------

// CO20190803 START
//---------------------------------------------------------------------------------
//  class xQMVASP
//---------------------------------------------------------------------------------
bool xQMVASP::GetProperties(const string& stringIN, bool QUIET) { // CO20191110
  stringstream sss;
  sss.str(stringIN);
  if (filename.empty()) {
    filename = DEFAULT_AFLOW_QMVASP_OUT;
  }
  return xQMVASP::GetProperties(sss, QUIET);
}

bool xQMVASP::GetPropertiesFile(const string& fileIN, bool QUIET) { // CO20191110
  stringstream sss;
  if (filename.empty()) {
    filename = fileIN;
  }
  aurostd::compressfile2stringstream(fileIN, sss);
  return xQMVASP::GetProperties(sss, QUIET);
}

bool xQMVASP::GetPropertiesUrlFile(const string& url, const string& file, bool QUIET) { // CO20191110
  const string tmpfile = aurostd::TmpFileCreate("xQMVASP_GetProperties"); // CO20200502 - threadID
  aurostd::httpGetFileStatus(url + "/" + file, tmpfile);
  const bool out = GetPropertiesFile(tmpfile, QUIET);
  filename = "url=" + url; // CO20210315
  aurostd::RemoveFile(tmpfile);
  return out;
}

bool xQMVASP::GetProperties(const stringstream& stringstreamIN, bool QUIET) { // CO20191110
  const bool LDEBUG = (false || XHOST.DEBUG || !QUIET);
  stringstream message;
  const bool force_exit = XHOST.POSTPROCESS; // SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only

  const bool ERROR_flag = false;
  const long double seconds = aurostd::get_seconds();
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BEGIN (" << time_delay(seconds) << ")" << endl;
  }
  clear(); // so it does not mess up vector/deque
  content = stringstreamIN.str();
  vcontent.clear();
  vector<string> tokens;
  vector<string> tokens2;
  stringstream ss_xstr;
  aurostd::string2vectorstring(content, vcontent);
  if (filename.empty()) {
    filename = "stringstream";
  }

  // crunching to eat the info
  H_atom_relax = AUROSTD_NAN; // this will be the LAST relax
  H_atom_static = AUROSTD_NAN;
  bool inside_relax = false;
  bool inside_static = false;
  // ----------------------------------------------------------------------
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vcontent.size()=" << vcontent.size() << endl;
  }
  // ----------------------------------------------------------------------
  for (size_t iline = 0; iline < vcontent.size(); iline++) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vcontent.at(" << iline << ")=" << vcontent[iline] << endl;
    }
    if (aurostd::substring2bool(vcontent[iline], "STOP_relax")) {
      inside_relax = false;
    } else if (aurostd::substring2bool(vcontent[iline], "START_relax")) {
      inside_relax = true;
    } else if (aurostd::substring2bool(vcontent[iline], "STOP_static")) {
      inside_static = false;
    } else if (aurostd::substring2bool(vcontent[iline], "START_static")) {
      inside_static = true;
    } else {
      if (aurostd::substring2bool(vcontent[iline], "H_atom")) {
        // H_atom=-9.859773400000e+00  (eV/at)
        tokens.clear();
        tokens2.clear();
        aurostd::string2tokens(vcontent[iline], tokens, "=");
        if (!tokens.empty()) {
          aurostd::string2tokens(tokens[1], tokens2, " ");
        }
        if (!tokens2.empty()) {
          if (!aurostd::isfloat(tokens2[0])) {
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "H_atom input cannot be parsed (!isfloat)", _INPUT_ERROR_);
          }
          if (inside_relax) {
            H_atom_relax = aurostd::string2utype<double>(tokens2[0]);
          } else if (inside_static) {
            H_atom_static = aurostd::string2utype<double>(tokens2[0]);
          }
        }
        // ME20191218 - substring2bool ignores comments and TOTAL-FORCE is in a comment line
      } else if (vcontent[iline].find("TOTAL-FORCE") != string::npos) { // CO20200106 - patching for auto-indenting
        vforces.clear();
        iline++; // skip first [AFLOW]
        // ME20191219 - ++iline needs to be in the while statement or the loop
        //  will never start
        while (iline < vcontent.size() && !aurostd::substring2bool(vcontent[++iline], "[AFLOW]")) {
          vforces.emplace_back(3);
          aurostd::string2tokens(vcontent[iline], tokens, " ");
          if (tokens.size() == 6) {
            for (int i = 1; i < 4; i++) {
              if (aurostd::isfloat(tokens[i + 2])) {
                vforces.back()[i] = aurostd::string2utype<double>(tokens[i + 2]);
              } else {
                throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Expected force input to be a number", _FILE_CORRUPT_);
              }
            }
          } else {
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Unexpected count of force components", _FILE_CORRUPT_);
          }
        }
      } else if (vcontent[iline].find("VASP_POSCAR_MODE_EXPLICIT") != string::npos) {
        aurostd::StringstreamClean(ss_xstr);
        while (iline < vcontent.size() && !aurostd::substring2bool(vcontent[++iline], "[AFLOW]")) {
          ss_xstr << vcontent[iline] << endl;
        }
        xstr_final.initialize(ss_xstr);
      }
    }
  }

  // ----------------------------------------------------------------------
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " END (" << time_delay(seconds) << ")" << endl;
  }
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  if (ERROR_flag) {
    message << "ERROR_flag set in xQMVASP";
    if (force_exit) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    } else {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_, QUIET);
    }
    return false;
  }
  return true;
}

// JSON serialization of xQMVASP
aurostd::JSON::object xQMVASP::serialize() const {
  return aurostd::JSON::object({AST_JSON_GETTER(JSON_xQMVASP_MEMBERS)});
}

xQMVASP xQMVASP::deserialize(const aurostd::JSON::object& jo) {
  AST_JSON_SETTER(JSON_xQMVASP_MEMBERS)
  return *this;
}

//-------------------------------------------------------------------------------------------------
// CO20190803 STOP

//-------------------------------------------------------------------------------------------------
// CO20190803 STOP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
