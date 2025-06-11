// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo - Kesong Yang - Camilo Calderon - Corey Oses
// fixed for xz - SC 2018

#ifndef _AFLOW_OVASP_CPP_
#define _AFLOW_OVASP_CPP_

#include "aflow.h"

//---------------------------------------------------------------------------------
// for bechmark
string time_delay(long double seconds) { return aurostd::utype2string<long double>(aurostd::get_seconds(seconds),6); }

//---------------------------------------------------------------------------------
// class xOUTCAR
//---------------------------------------------------------------------------------
//ME20200427 - included xStream::initialize
xOUTCAR::xOUTCAR(ostream& oss) : xStream(oss) {free();} //CO20200404 - xStream integration for logging
xOUTCAR::xOUTCAR(ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss) {free();} //CO20200404 - xStream integration for logging
xOUTCAR::xOUTCAR(const string& fileIN,bool QUIET,ostream& oss) : xStream(oss) {m_initialized=initialize(fileIN,QUIET);}  //CO20200404 - xStream integration for logging
xOUTCAR::xOUTCAR(const string& fileIN,ofstream& FileMESSAGE,bool QUIET,ostream& oss) : xStream(FileMESSAGE,oss) {m_initialized=initialize(fileIN,QUIET);}  //CO20200404 - xStream integration for logging

bool xOUTCAR::initialize(const string& fileIN,ostream& oss,bool QUIET) {xStream::initialize(oss);return initialize(fileIN,QUIET);} //CO20200508
bool xOUTCAR::initialize(const string& fileIN,ofstream& FileMESSAGE,ostream& oss,bool QUIET) {xStream::initialize(FileMESSAGE,oss);return initialize(fileIN,QUIET);} //CO20200508
bool xOUTCAR::initialize(const string& fileIN,bool QUIET) {
  free();
  filename = fileIN;
  return GetPropertiesFile(fileIN,QUIET);
}

xOUTCAR::xOUTCAR(const xOUTCAR& b) : xStream(*b.getOFStream(),*b.getOSS()) {free();copy(b);} // copy PUBLIC

xOUTCAR::~xOUTCAR() {xStream::free();free();}  //CO20200404 - xStream integration for logging

void xOUTCAR::free() {
  //------------------------------------------------------------------------------
  m_initialized=false; //CO20200404
  // GetProperties
  content="";                   // for aflowlib_libraries.cpp
  vcontent.clear();             // for aflowlib_libraries.cpp
  filename="";                  // for aflowlib_libraries.cpp
  SYSTEM = "";                  // for aflowlib_libraries.cpp 
  NELM=0;                       // for xwarnings  //CO20200624
  NIONS=0;                      // for aflowlib_libraries.cpp
  Efermi=0.0;                   // for aflowlib_libraries.cpp
  isLSCOUPLING=FALSE;           // for aflowlib_libraries.cpp
  efield_pead.clear();          // for ivasp
  nelectrons = 0;               // AS20200528
  natoms=0.0;                   // for aflowlib_libraries.cpp
  energy_cell=0.0;              // for aflowlib_libraries.cpp
  energy_atom=0.0;              // for aflowlib_libraries.cpp
  enthalpy_cell=0.0;            // for aflowlib_libraries.cpp
  enthalpy_atom=0.0;            // for aflowlib_libraries.cpp
  eentropy_cell=0.0;            // for aflowlib_libraries.cpp
  eentropy_atom=0.0;            // for aflowlib_libraries.cpp
  PV_cell=0.0;                  // for aflowlib_libraries.cpp
  PV_atom=0.0;                  // for aflowlib_libraries.cpp
  stress.clear();               // for aflowlib_libraries.cpp
  mag_cell=0.0;                 // for aflowlib_libraries.cpp
  mag_atom=0.0;                 // for aflowlib_libraries.cpp
  vmag.clear();                 // for aflowlib_libraries.cpp
  vmag_noncoll.clear();         //DX20171205 - non-collinear magnetization 
  volume_cell=0.0;              // for aflowlib_libraries.cpp
  volume_atom=0.0;              // for aflowlib_libraries.cpp
  pressure=0.0;                 // for aflowlib_libraries.cpp // copy of PSTRESS
  pressure_residual=0.0;        // for aflowlib_libraries.cpp
  Pulay_stress=0.0;             // for aflowlib_libraries.cpp
  vforces.clear();              // for aflowlib_libraries.cpp
  vpositions_cartesian.clear(); // for aflowlib_libraries.cpp
  ENCUT=0.0;EDIFF=0.0;EDIFFG=0.0;POTIM=0.0;TEIN=0.0;TEBEG=0.0;TEEND=0.0;SMASS=0.0;NPACO=0.0;APACO=0.0;PSTRESS=0.0;     // 
  NBANDS=0;NKPTS=0;NSW=0;NBLOCK=0;KBLOCK=0;IBRION=0;NFREE=0;ISIF=0;IWAVPR=0;ISYM=0;ISPIN=0; //
  total_energy_change=0.0;
  // DOS related values:
  EMIN=0.0;EMAX=0.0;SIGMA=0.0;ISMEAR=0;  // for aflowlib_libraries.cpp 
  //  Electronic relaxation
  IALGO=0;                      // for aflowlib_libraries.cpp
  LDIAG="";                     // for aflowlib_libraries.cpp
  IMIX=0;INIMIX=0;MIXPRE=0;     // for aflowlib_libraries.cpp
  AMIX=0.0;BMIX=0.0;AMIX_MAG=0.0;BMIX_MAG=0.0;AMIN=0.0;WC=0.0;  // for aflowlib_libraries.cpp
  // Intra band minimization
  WEIMIN=0.0;EBREAK=0.0;DEPER=0.0;TIME=0.0;  // for aflowlib_libraries.cpp
  // begin shared xPOTCAR
  ENMAX=0.0;vENMAX.clear();                                      // for aflowlib_libraries.cpp
  ENMIN=0.0;vENMIN.clear();                                      // for aflowlib_libraries.cpp
  POMASS_sum=0.0;POMASS_min=0.0;POMASS_max=0.0;vPOMASS.clear();  // for aflowlib_libraries.cpp
  ZVAL_sum=0.0;ZVAL_min=0.0;ZVAL_max=0.0;vZVAL.clear();          // for aflowlib_libraries.cpp
  EATOM_min=0.0;EATOM_max=0.0;vEATOM.clear();          // for aflowlib_libraries.cpp
  RCORE_min=0.0;RCORE_max=0.0;vRCORE.clear();          // for aflowlib_libraries.cpp
  RWIGS_min=0.0;RWIGS_max=0.0;vRWIGS.clear();          // for aflowlib_libraries.cpp
  EAUG_min=0.0;EAUG_max=0.0;vEAUG.clear();              // for aflowlib_libraries.cpp
  RAUG_min=0.0;RAUG_max=0.0;vRAUG.clear();              // for aflowlib_libraries.cpp
  RMAX_min=0.0;RMAX_max=0.0;vRMAX.clear();              // unicity
  vTITEL.clear();                                                // unicity
  vLEXCH.clear();                                                // unicity 
  // end shared xPOTCAR
  pp_type="";                   // for aflowlib_libraries.cpp
  species.clear();              // for aflowlib_libraries.cpp
  species_Z.clear();            // for aflowlib_libraries.cpp
  species_pp.clear();           // for aflowlib_libraries.cpp
  species_pp_type.clear();      // for aflowlib_libraries.cpp
  species_pp_version.clear();   // for aflowlib_libraries.cpp
  species_pp_AUID.clear();      // for aflowlib_libraries.cpp
  species_pp_AUID_collisions.clear();   // for aflowlib_libraries.cpp
  species_pp_groundstate_energy.clear();      // for aflowlib_libraries.cpp
  species_pp_groundstate_structure.clear();   // for aflowlib_libraries.cpp
  species_pp_vLDAU.clear();     // for aflowlib_libraries.cpp
  string_LDAU="";               // for aflowlib_libraries.cpp
  isKIN=FALSE;                  // for aflowlib_libraries.cpp
  isMETAGGA=FALSE;              // for aflowlib_libraries.cpp
  METAGGA="";                   // for aflowlib_libraries.cpp
  //------------------------------------------------------------------------------
  nweights=0;nkpoints_irreducible=0;   // for aflowlib_libraries.cpp
  vkpoint_reciprocal.clear();      // for aflowlib_libraries.cpp
  vkpoint_cartesian.clear();       // for aflowlib_libraries.cpp
  vweights.clear();                    // for aflowlib_libraries.cpp
  //------------------------------------------------------------------------------
  calculation_time=0.0;         // for aflowlib_libraries.cpp 
  calculation_memory=0.0;       // for aflowlib_libraries.cpp 
  calculation_cores=1;          // for aflowlib_libraries.cpp 
  //------------------------------------------------------------------------------
  // GetEffectiveMass
  band_index.clear();           // for aflowlib_libraries.cpp 
  carrier_type.clear();         // for aflowlib_libraries.cpp 
  carrier_spin.clear();         // for aflowlib_libraries.cpp 
  extrema_cart_coord.clear();   // for aflowlib_libraries.cpp 
  effective_mass_axes.clear();  // for aflowlib_libraries.cpp 
  equivalent_valley.clear();    // for aflowlib_libraries.cpp 
  effective_mass_DOS.clear();   // for aflowlib_libraries.cpp 
  effective_mass_COND.clear();  // for aflowlib_libraries.cpp 
  mass_elec_dos.clear();        // for aflowlib_libraries.cpp 
  mass_hole_dos.clear();        // for aflowlib_libraries.cpp 
  mass_elec_conduction.clear(); // for aflowlib_libraries.cpp 
  mass_hole_conduction.clear(); // for aflowlib_libraries.cpp 
  //------------------------------------------------------------------------------
  // GetBandGap
  xstr.clear(); //DX20191220 - uppercase to lowercase clear
  conduction_band_min.clear();  // for aflowlib_libraries.cpp 
  valence_band_max.clear();     // for aflowlib_libraries.cpp 
  Egap_type.clear();            // for aflowlib_libraries.cpp 
  Egap_fit.clear();             // for aflowlib_libraries.cpp 
  Egap.clear();                 // for aflowlib_libraries.cpp 
  conduction_band_min_net = 0.0;// for aflowlib_libraries.cpp 
  valence_band_max_net = 0.0;   // for aflowlib_libraries.cpp 
  Egap_type_net="";             // for aflowlib_libraries.cpp 
  //Egap_fit_net = 0.0;           // for aflowlib_libraries.cpp     //CO
  //Egap_net = 0.0;               // for aflowlib_libraries.cpp     //CO
  Egap_fit_net = AUROSTD_NAN;   // for aflowlib_libraries.cpp 
  Egap_net = AUROSTD_NAN;       // for aflowlib_libraries.cpp 

  //[CO20200404 - DUPLICATE]ERROR = "";                   // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]//------------------------------------------------------------------------------
  //[CO20200404 - DUPLICATE]// GetProperties
  //[CO20200404 - DUPLICATE]vcontent.clear();                    // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vmag.clear();                        // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vmag_noncoll.clear();                //DX20171205 - non-collinear magnetization 
  //[CO20200404 - DUPLICATE]vforces.clear();                     // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vpositions_cartesian.clear();        // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]// begin shared xPOTCAR
  //[CO20200404 - DUPLICATE]vENMAX.clear();                      // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vENMIN.clear();                      // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vPOMASS.clear();                     // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vZVAL.clear();                       // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vEATOM.clear();                      // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vRCORE.clear();                      // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vRWIGS.clear();                      // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vEAUG.clear();                       // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vRAUG.clear();                       // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vRMAX.clear();                       // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vTITEL.clear();                      // unicity
  //[CO20200404 - DUPLICATE]vLEXCH.clear();                      // unicity
  //[CO20200404 - DUPLICATE]// end shared xPOTCAR
  //[CO20200404 - DUPLICATE]species.clear();                     // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]species_Z.clear();                   // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]species_pp.clear();                  // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]species_pp_type.clear();             // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]species_pp_version.clear();          // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]species_pp_AUID.clear();             // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]species_pp_AUID_collisions.clear();  // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]species_pp_groundstate_energy.clear();      // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]species_pp_groundstate_structure.clear();   // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]species_pp_vLDAU.clear();            // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vkpoint_reciprocal.clear();          // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vkpoint_cartesian.clear();           // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vweights.clear();                    // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]//------------------------------------------------------------------------------
  //[CO20200404 - DUPLICATE]// GetEffectiveMass
  //[CO20200404 - DUPLICATE]band_index.clear();                  // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]carrier_type.clear();                // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]carrier_spin.clear();                // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]extrema_cart_coord.clear();          // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]effective_mass_axes.clear();         // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]equivalent_valley.clear();           // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]effective_mass_DOS.clear();          // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]effective_mass_COND.clear();         // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]mass_elec_dos.clear();               // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]mass_hole_dos.clear();               // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]mass_elec_conduction.clear();        // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]mass_hole_conduction.clear();        // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]//------------------------------------------------------------------------------
  //[CO20200404 - DUPLICATE]// GetBandGap
  //[CO20200404 - DUPLICATE]xstr.clear(); //DX20191220 - uppercase to lowercase clear
  //[CO20200404 - DUPLICATE]conduction_band_min.clear();         // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]valence_band_max.clear();            // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]Egap_type.clear();                   // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]Egap.clear();                        // for aflowlib_libraries.cpp
}

void xOUTCAR::copy(const xOUTCAR& b) { // copy PRIVATE
  xStream::copy(b);  //CO20200404 - xStream integration for logging
  free();
  m_initialized=b.m_initialized; //CO20200404
  content=b.content;
  vcontent.clear(); for(uint i=0;i<b.vcontent.size();i++) vcontent.push_back(b.vcontent.at(i));  // for aflowlib_libraries.cpp
  filename=b.filename;
  SYSTEM = b.SYSTEM; // camilo
  NELM=b.NELM;  //CO20200624
  NIONS=b.NIONS;
  Efermi=b.Efermi;
  isLSCOUPLING=b.isLSCOUPLING;
  efield_pead=b.efield_pead;  //CO20210315
  nelectrons=b.nelectrons; //AS20200528
  natoms=b.natoms;                              // for aflowlib_libraries.cpp
  energy_cell=b.energy_cell;                    // for aflowlib_libraries.cpp
  energy_atom=b.energy_atom;                    // for aflowlib_libraries.cpp
  enthalpy_cell=b.enthalpy_cell;                // for aflowlib_libraries.cpp
  enthalpy_atom=b.enthalpy_atom;                // for aflowlib_libraries.cpp
  eentropy_cell=b.eentropy_cell;                // for aflowlib_libraries.cpp
  eentropy_atom=b.eentropy_atom;                // for aflowlib_libraries.cpp
  PV_cell=b.PV_cell;                            // for aflowlib_libraries.cpp
  PV_atom=b.PV_atom;                            // for aflowlib_libraries.cpp
  stress=b.stress;                              // for aflowlib_libraries.cpp
  mag_cell=b.mag_cell;                          // for aflowlib_libraries.cpp
  mag_atom=b.mag_atom;                          // for aflowlib_libraries.cpp
  vmag.clear(); for(uint i=0;i<b.vmag.size();i++) vmag.push_back(b.vmag.at(i));  // for aflowlib_libraries.cpp
  vmag_noncoll.clear(); for(uint i=0;i<b.vmag_noncoll.size();i++) vmag_noncoll.push_back(b.vmag_noncoll.at(i));  //DX20171205
  volume_cell=b.volume_cell;                    // for aflowlib_libraries.cpp
  volume_atom=b.volume_atom;                    // for aflowlib_libraries.cpp
  pressure=b.pressure;                          // for aflowlib_libraries.cpp
  pressure_residual=b.pressure_residual;        // for aflowlib_libraries.cpp
  Pulay_stress=b.Pulay_stress;                  // for aflowlib_libraries.cpp
  vforces.clear(); for(uint i=0;i<b.vforces.size();i++) vforces.push_back(b.vforces.at(i));  // for aflowlib_libraries.cpp
  vpositions_cartesian.clear(); for(uint i=0;i<b.vpositions_cartesian.size();i++) vpositions_cartesian.push_back(b.vpositions_cartesian.at(i)); // for aflowlib_libraries.cpp
  ENCUT=b.ENCUT;EDIFF=b.EDIFF;EDIFFG=b.EDIFFG;POTIM=b.POTIM;TEIN=b.TEIN;TEBEG=b.TEBEG;TEEND=b.TEEND;SMASS=b.SMASS;NPACO=b.NPACO;APACO=b.APACO;PSTRESS=b.PSTRESS;     // 
  NBANDS=b.NBANDS;NKPTS=b.NKPTS;NSW=b.NSW;NBLOCK=b.NBLOCK;KBLOCK=b.KBLOCK;IBRION=b.IBRION;NFREE=b.NFREE;ISIF=b.ISIF;IWAVPR=b.IWAVPR;ISYM=b.ISYM;ISPIN=b.ISPIN; //
  total_energy_change=b.total_energy_change;
  // DOS related values:
  EMIN=b.EMIN;  // for aflowlib_libraries.cpp 
  EMAX=b.EMAX;  // for aflowlib_libraries.cpp 
  SIGMA=b.SIGMA;  // for aflowlib_libraries.cpp 
  ISMEAR=b.ISMEAR;  // for aflowlib_libraries.cpp 
  //  Electronic relaxation
  IALGO=b.IALGO;                                // for aflowlib_libraries.cpp
  LDIAG=b.LDIAG;                                // for aflowlib_libraries.cpp
  IMIX=b.IMIX;INIMIX=b.INIMIX;MIXPRE=b.MIXPRE;  // for aflowlib_libraries.cpp
  AMIX=b.AMIX;BMIX=b.BMIX;AMIX_MAG=b.AMIX_MAG;BMIX_MAG=b.BMIX_MAG;AMIN=b.AMIN;WC=b.WC;  // for aflowlib_libraries.cpp
  // Intra band minimization
  WEIMIN=b.WEIMIN;EBREAK=b.EBREAK;DEPER=b.DEPER;TIME=b.TIME;  // for aflowlib_libraries.cpp
  // begin shared xPOTCAR
  ENMAX=b.ENMAX;
  vENMAX.clear(); for(uint i=0;i<b.vENMAX.size();i++) vENMAX.push_back(b.vENMAX.at(i)); // for aflowlib_libraries.cpp
  ENMIN=b.ENMIN;
  vENMIN.clear(); for(uint i=0;i<b.vENMIN.size();i++) vENMIN.push_back(b.vENMIN.at(i)); // for aflowlib_libraries.cpp
  POMASS_sum=b.POMASS_sum;POMASS_min=b.POMASS_min;POMASS_max=b.POMASS_max;
  vPOMASS.clear(); for(uint i=0;i<b.vPOMASS.size();i++) vPOMASS.push_back(b.vPOMASS.at(i)); // for aflowlib_libraries.cpp
  ZVAL_sum=b.ZVAL_sum;ZVAL_min=b.ZVAL_min;ZVAL_max=b.ZVAL_max;
  vZVAL.clear(); for(uint i=0;i<b.vZVAL.size();i++) vZVAL.push_back(b.vZVAL.at(i)); // for aflowlib_libraries.cpp
  EATOM_min=b.EATOM_min;EATOM_max=b.EATOM_max;
  vEATOM.clear(); for(uint i=0;i<b.vEATOM.size();i++) vEATOM.push_back(b.vEATOM.at(i)); // for aflowlib_libraries.cpp
  RCORE_min=b.RCORE_min;RCORE_max=b.RCORE_max;
  vRCORE.clear(); for(uint i=0;i<b.vRCORE.size();i++) vRCORE.push_back(b.vRCORE.at(i)); // for aflowlib_libraries.cpp
  RWIGS_min=b.RWIGS_min;RWIGS_max=b.RWIGS_max;
  vRWIGS.clear(); for(uint i=0;i<b.vRWIGS.size();i++) vRWIGS.push_back(b.vRWIGS.at(i)); // for aflowlib_libraries.cpp
  EAUG_min=b.EAUG_min;EAUG_max=b.EAUG_max;
  vEAUG.clear(); for(uint i=0;i<b.vEAUG.size();i++) vEAUG.push_back(b.vEAUG.at(i)); // for aflowlib_libraries.cpp
  RAUG_min=b.RAUG_min;RAUG_max=b.RAUG_max;
  vRAUG.clear(); for(uint i=0;i<b.vRAUG.size();i++) vRAUG.push_back(b.vRAUG.at(i)); // for aflowlib_libraries.cpp
  RMAX_min=b.RMAX_min;RMAX_max=b.RMAX_max;
  vRMAX.clear(); for(uint i=0;i<b.vRMAX.size();i++) vRMAX.push_back(b.vRMAX.at(i)); // for aflowlib_libraries.cpp
  vTITEL.clear(); for(uint i=0;i<b.vTITEL.size();i++) vTITEL.push_back(b.vTITEL.at(i)); // for aflowlib_libraries.cpp
  vLEXCH.clear(); for(uint i=0;i<b.vLEXCH.size();i++) vLEXCH.push_back(b.vLEXCH.at(i)); // for aflowlib_libraries.cpp
  // end shared xPOTCAR
  pp_type=b.pp_type;                            // for aflowlib_libraries.cpp
  species.clear(); for(uint i=0;i<b.species.size();i++) species.push_back(b.species.at(i)); // for aflowlib_libraries.cpp
  species_Z.clear(); for(uint i=0;i<b.species_Z.size();i++) species_Z.push_back(b.species_Z.at(i)); // for aflowlib_libraries.cpp
  species_pp.clear(); for(uint i=0;i<b.species_pp.size();i++) species_pp.push_back(b.species_pp.at(i)); // for aflowlib_libraries.cpp
  species_pp_type.clear(); for(uint i=0;i<b.species_pp_type.size();i++) species_pp_type.push_back(b.species_pp_type.at(i)); // for aflowlib_libraries.cpp
  species_pp_version.clear(); for(uint i=0;i<b.species_pp_version.size();i++) species_pp_version.push_back(b.species_pp_version.at(i));  // for aflowlib_libraries.cpp
  species_pp_AUID.clear(); for(uint i=0;i<b.species_pp_AUID.size();i++) species_pp_AUID.push_back(b.species_pp_AUID.at(i));  // for aflowlib_libraries.cpp
  species_pp_AUID_collisions.clear(); for(uint i=0;i<b.species_pp_AUID_collisions.size();i++) species_pp_AUID_collisions.push_back(b.species_pp_AUID_collisions.at(i));  // for aflowlib_libraries.cpp
  species_pp_groundstate_energy.clear(); for(uint i=0;i<b.species_pp_groundstate_energy.size();i++) species_pp_groundstate_energy.push_back(b.species_pp_groundstate_energy.at(i));  // for aflowlib_libraries.cpp
  species_pp_groundstate_structure.clear(); for(uint i=0;i<b.species_pp_groundstate_structure.size();i++) species_pp_groundstate_structure.push_back(b.species_pp_groundstate_structure.at(i));  // for aflowlib_libraries.cpp
  species_pp_vLDAU.clear(); for(uint i=0;i<b.species_pp_vLDAU.size();i++) species_pp_vLDAU.push_back(b.species_pp_vLDAU.at(i));
  string_LDAU=b.string_LDAU;                     // for aflowlib_libraries.cpp
  isKIN=b.isKIN;                                 // for aflowlib_libraries.cpp
  isMETAGGA=b.isMETAGGA;                         // for aflowlib_libraries.cpp
  METAGGA=b.METAGGA;                             // for aflowlib_libraries.cpp
  // KPOINTS
  nweights=b.nweights;                           // for aflowlib_libraries.cpp
  nkpoints_irreducible=b.nkpoints_irreducible;   // for aflowlib_libraries.cpp
  vkpoint_reciprocal.clear();for(uint i=0;i<b.vkpoint_reciprocal.size();i++) vkpoint_reciprocal.push_back(b.vkpoint_reciprocal.at(i));      // for aflowlib_libraries.cpp
  vkpoint_cartesian.clear();for(uint i=0;i<b.vkpoint_cartesian.size();i++) vkpoint_cartesian.push_back(b.vkpoint_cartesian.at(i));       // for aflowlib_libraries.cpp
  vweights.clear();for(uint i=0;i<b.vweights.size();i++) vweights.push_back(b.vweights.at(i));                   // for aflowlib_libraries.cpp
  // times
  calculation_time=b.calculation_time;          // for aflowlib_libraries.cpp 
  calculation_memory=b.calculation_memory;      // for aflowlib_libraries.cpp 
  calculation_cores=b.calculation_cores;        // for aflowlib_libraries.cpp 
  //------------------------------------------------------------------------------
  // GetEffectiveMass
  band_index.clear(); for(uint i=0;i<b.band_index.size();i++) band_index.push_back(b.band_index.at(i));
  carrier_type.clear(); for(uint i=0;i<b.carrier_type.size();i++) carrier_type.push_back(b.carrier_type.at(i));
  carrier_spin.clear(); for(uint i=0;i<b.carrier_spin.size();i++) carrier_spin.push_back(b.carrier_spin.at(i));
  extrema_cart_coord.clear(); for(uint i=0;i<b.extrema_cart_coord.size();i++) extrema_cart_coord.push_back(b.extrema_cart_coord.at(i));
  effective_mass_axes.clear(); for(uint i=0;i<b.effective_mass_axes.size();i++) effective_mass_axes.push_back(b.effective_mass_axes.at(i));
  equivalent_valley.clear(); for(uint i=0;i<b.equivalent_valley.size();i++) equivalent_valley.push_back(b.equivalent_valley.at(i));
  effective_mass_DOS.clear(); for(uint i=0;i<b.effective_mass_DOS.size();i++) effective_mass_DOS.push_back(b.effective_mass_DOS.at(i));
  effective_mass_COND.clear(); for(uint i=0;i<b.effective_mass_COND.size();i++) effective_mass_COND.push_back(b.effective_mass_COND.at(i));
  mass_elec_dos.clear(); for(uint i=0;i<b.mass_elec_dos.size();i++) mass_elec_dos.push_back(b.mass_elec_dos.at(i));
  mass_hole_dos.clear(); for(uint i=0;i<b.mass_hole_dos.size();i++) mass_hole_dos.push_back(b.mass_hole_dos.at(i));
  mass_elec_conduction.clear(); for(uint i=0;i<b.mass_elec_conduction.size();i++) mass_elec_conduction.push_back(b.mass_elec_conduction.at(i));
  mass_hole_conduction.clear(); for(uint i=0;i<b.mass_hole_conduction.size();i++) mass_hole_conduction.push_back(b.mass_hole_conduction.at(i));
  //------------------------------------------------------------------------------
  // GetBandGap
  xstr=b.xstr;
  conduction_band_min.clear(); for(uint i=0;i<b.conduction_band_min.size();i++) conduction_band_min.push_back(b.conduction_band_min.at(i));
  valence_band_max.clear(); for(uint i=0;i<b.valence_band_max.size();i++) valence_band_max.push_back(b.valence_band_max.at(i));
  Egap_type.clear(); for(uint i=0;i<b.Egap_type.size();i++) Egap_type.push_back(b.Egap_type.at(i));
  Egap_fit.clear(); for(uint i=0;i<b.Egap_fit.size();i++) Egap_fit.push_back(b.Egap_fit.at(i));
  Egap.clear(); for(uint i=0;i<b.Egap.size();i++) Egap.push_back(b.Egap.at(i));
  conduction_band_min_net = b.conduction_band_min_net;
  valence_band_max_net = b.valence_band_max_net;
  Egap_type_net = b.Egap_type_net;
  Egap_fit_net = b.Egap_net;
  Egap_net = b.Egap_net;
  //[CO20200404 - OBSOLETE]ERROR = b.ERROR ; // camilo
}

const xOUTCAR& xOUTCAR::operator=(const xOUTCAR& b) {  // operator= PUBLIC
  if(this!=&b) {free();copy(b);}
  return *this;
}

void xOUTCAR::clear() {  // clear PRIVATE
  xOUTCAR _temp;
  string filename_aus=filename;
  copy(_temp);
  filename=filename_aus;
}

ostream& operator<<(ostream& oss, const xOUTCAR& xOUT) {  //SC20200330
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  long double seconds=aurostd::get_seconds();
  if(LDEBUG) cerr << XPID << "xOUTCAR::operator<<: ---------------------------------" << endl;
  if(LDEBUG) cerr << XPID << "xOUTCAR::operator<<: BEGIN (" << time_delay(seconds) << ")" << endl;
  if(LDEBUG) cerr << XPID << "xOUTCAR::operator<<: filename=[" << xOUT.filename << "]" << endl;
  oss << " filename=" << xOUT.filename<< endl;
  oss << " vcontent.size()=" << xOUT.vcontent.size() << endl;
  oss << " SYSTEM=" << xOUT.SYSTEM << endl;
  oss << " NELM=" << xOUT.NELM << endl; //CO20200624
  oss << " NIONS=" << xOUT.NIONS << endl;
  oss << " Efermi=" << xOUT.Efermi << endl;
  oss << " isLSCOUPLING=" << xOUT.isLSCOUPLING << endl;
  oss << " nelectrons=" << xOUT.nelectrons << endl; //AS20200528
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
  oss << " ENMAX=" << xOUT.ENMAX << "  vENMAX.size()=" << xOUT.vENMAX.size() << ": "; for(uint i=0;i<xOUT.vENMAX.size();i++) { oss << xOUT.vENMAX.at(i) << " "; } oss << endl;
  oss << " ENMIN=" << xOUT.ENMIN << "  vENMIN.size()=" << xOUT.vENMIN.size() << ": "; for(uint i=0;i<xOUT.vENMIN.size();i++) { oss << xOUT.vENMIN.at(i) << " "; } oss << endl;
  oss << " POMASS_sum=" << xOUT.POMASS_sum << " POMASS_min=" << xOUT.POMASS_min << " POMASS_max=" << xOUT.POMASS_max << " vPOMASS.size()=" << xOUT.vPOMASS.size() << ": "; for(uint i=0;i<xOUT.vPOMASS.size();i++) { oss << xOUT.vPOMASS.at(i) << " "; } oss << " " << endl;
  oss << " ZVAL_sum=" << xOUT.ZVAL_sum << " ZVAL_min=" << xOUT.ZVAL_min << " ZVAL_max=" << xOUT.ZVAL_max << " vZVAL.size()=" << xOUT.vZVAL.size() << ": "; for(uint i=0;i<xOUT.vZVAL.size();i++) { oss << xOUT.vZVAL.at(i) << " "; } oss << " " << endl;
  oss << " EATOM_min=" << xOUT.EATOM_min << " EATOM_max=" << xOUT.EATOM_max << " vEATOM.size()=" << xOUT.vEATOM.size() << ": "; for(uint i=0;i<xOUT.vEATOM.size();i++) { oss << xOUT.vEATOM.at(i) << " "; } oss << " " << endl;
  oss << " RCORE_min=" << xOUT.RCORE_min << " RCORE_max=" << xOUT.RCORE_max << " vRCORE.size()=" << xOUT.vRCORE.size() << ": "; for(uint i=0;i<xOUT.vRCORE.size();i++) { oss << xOUT.vRCORE.at(i) << " "; } oss << " " << endl;
  oss << " RWIGS_min=" << xOUT.RWIGS_min << " RWIGS_max=" << xOUT.RWIGS_max << " vRWIGS.size()=" << xOUT.vRWIGS.size() << ": "; for(uint i=0;i<xOUT.vRWIGS.size();i++) { oss << xOUT.vRWIGS.at(i) << " "; } oss << " " << endl;
  oss << " EAUG_min=" << xOUT.EAUG_min << " EAUG_max=" << xOUT.EAUG_max << " vEAUG.size()=" << xOUT.vEAUG.size() << ": "; for(uint i=0;i<xOUT.vEAUG.size();i++) { oss << xOUT.vEAUG.at(i) << " "; } oss << " " << endl;
  oss << " RAUG_min=" << xOUT.RAUG_min << " RAUG_max=" << xOUT.RAUG_max << " vRAUG.size()=" << xOUT.vRAUG.size() << ": "; for(uint i=0;i<xOUT.vRAUG.size();i++) { oss << xOUT.vRAUG.at(i) << " "; } oss << " " << endl;
  oss << " RMAX_min=" << xOUT.RMAX_min << " RMAX_max=" << xOUT.RMAX_max << " vRMAX.size()=" << xOUT.vRMAX.size() << ": "; for(uint i=0;i<xOUT.vRMAX.size();i++) { oss << xOUT.vRMAX.at(i) << " "; } oss << " " << endl;
  oss << " vTITEL.size()=" << xOUT.vTITEL.size() << ": "; for(uint i=0;i<xOUT.vTITEL.size();i++) { oss << xOUT.vTITEL.at(i) << " "; } oss << endl;
  oss << " vLEXCH.size()=" << xOUT.vLEXCH.size() << ": "; for(uint i=0;i<xOUT.vLEXCH.size();i++) { oss << xOUT.vLEXCH.at(i) << " "; } oss << endl;
  oss << " pp_type=" << xOUT.pp_type << endl;
  oss << " species.size()=" << xOUT.species.size() << ": "; for(uint i=0;i<xOUT.species.size();i++) { oss << xOUT.species.at(i) << " "; } oss << endl;
  oss << " species_Z.size()=" << xOUT.species_Z.size() << ": "; for(uint i=0;i<xOUT.species_Z.size();i++) { oss << xOUT.species_Z.at(i) << " "; } oss << endl;
  oss << " species_pp.size()=" << xOUT.species_pp.size() << ": "; for(uint i=0;i<xOUT.species_pp.size();i++) { oss << xOUT.species_pp.at(i) << " "; } oss << endl;
  oss << " species_pp_type.size()=" << xOUT.species_pp_type.size() << ": "; for(uint i=0;i<xOUT.species_pp_type.size();i++) { oss << xOUT.species_pp_type.at(i) << " "; } oss << endl;
  oss << " species_pp_version.size()=" << xOUT.species_pp_version.size() << ": "; for(uint i=0;i<xOUT.species_pp_version.size();i++) { oss << xOUT.species_pp_version.at(i) << " "; } oss << endl;
  oss << " species_pp_AUID.size()=" << xOUT.species_pp_AUID.size() << ": "; for(uint i=0;i<xOUT.species_pp_AUID.size();i++) { oss << xOUT.species_pp_AUID.at(i) << " "; } oss << endl;
  oss << " species_pp_AUID_collisions.size()=" << xOUT.species_pp_AUID_collisions.size() << ": "; for(uint i=0;i<xOUT.species_pp_AUID_collisions.size();i++) { oss << xOUT.species_pp_AUID_collisions.at(i) << " "; } oss << endl;
  oss << " species_pp_groundstate_energy.size()=" << xOUT.species_pp_groundstate_energy.size() << ": "; for(uint i=0;i<xOUT.species_pp_groundstate_energy.size();i++) { oss << xOUT.species_pp_groundstate_energy.at(i) << " "; } oss << endl;
  oss << " species_pp_groundstate_structure.size()=" << xOUT.species_pp_groundstate_structure.size() << ": "; for(uint i=0;i<xOUT.species_pp_groundstate_structure.size();i++) { oss << xOUT.species_pp_groundstate_structure.at(i) << " "; } oss << endl;
  // deque<deque<double> > species_pp_vLDAU;  // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  oss << " isKIN=" << xOUT.isKIN << endl;
  oss << " isMETAGGA=" << xOUT.isMETAGGA << endl;
  oss << " METAGGA=" << xOUT.METAGGA << endl;
  oss << " nweights=" << xOUT.nweights << endl;
  oss << " nkpoints_irreducible=" << xOUT.nkpoints_irreducible << endl;
  //vector<aurostd::xvector<double> > vkpoint_reciprocal;         // kpoints reading
  //vector<aurostd::xvector<double> > vkpoint_cartesian;          // kpoints reading
  oss << " vweights.size()=" << xOUT.vweights.size() << ": "; for(uint i=0;i<xOUT.vweights.size();i++) { oss << xOUT.vweights.at(i) << " "; } oss << endl;
  oss << " calculation_time=" << xOUT.calculation_time << endl;
  oss << " calculation_memory=" << xOUT.calculation_memory << endl;
  oss << " calculation_cores=" << xOUT.calculation_cores << endl;
  // xstructure xstr;                                              // for GetBandGap()
  // EFFECTIVE MASSES
  oss << " band_index.size()=" << xOUT.band_index.size() << ": "; for(uint i=0;i<xOUT.band_index.size();i++) { oss << xOUT.band_index.at(i) << " "; } oss << endl;
  oss << " carrier_spin.size()=" << xOUT.carrier_spin.size() << ": "; for(uint i=0;i<xOUT.carrier_spin.size();i++) { oss << xOUT.carrier_spin.at(i) << " "; } oss << endl;
  oss << " carrier_type.size()=" << xOUT.carrier_type.size() << ": "; for(uint i=0;i<xOUT.carrier_type.size();i++) { oss << xOUT.carrier_type.at(i) << " "; } oss << endl;
  // vector<vector<double> > extrema_cart_coord;
  // vector<vector<double> > effective_mass_axes;
  oss << " equivalent_valley.size()=" << xOUT.equivalent_valley.size() << ": "; for(uint i=0;i<xOUT.equivalent_valley.size();i++) { oss << xOUT.equivalent_valley.at(i) << " "; } oss << endl;
  oss << " effective_mass_DOS.size()=" << xOUT.effective_mass_DOS.size() << ": "; for(uint i=0;i<xOUT.effective_mass_DOS.size();i++) { oss << xOUT.effective_mass_DOS.at(i) << " "; } oss << endl;
  oss << " effective_mass_COND.size()=" << xOUT.effective_mass_COND.size() << ": "; for(uint i=0;i<xOUT.effective_mass_COND.size();i++) { oss << xOUT.effective_mass_COND.at(i) << " "; } oss << endl;
  oss << " mass_elec_dos.size()=" << xOUT.mass_elec_dos.size() << ": "; for(uint i=0;i<xOUT.mass_elec_dos.size();i++) { oss << xOUT.mass_elec_dos.at(i) << " "; } oss << endl;
  oss << " mass_hole_dos.size()=" << xOUT.mass_hole_dos.size() << ": "; for(uint i=0;i<xOUT.mass_hole_dos.size();i++) { oss << xOUT.mass_hole_dos.at(i) << " "; } oss << endl;
  oss << " mass_elec_conduction.size()=" << xOUT.mass_elec_conduction.size() << ": "; for(uint i=0;i<xOUT.mass_elec_conduction.size();i++) { oss << xOUT.mass_elec_conduction.at(i) << " "; } oss << endl;
  oss << " mass_hole_conduction.size()=" << xOUT.mass_hole_conduction.size() << ": "; for(uint i=0;i<xOUT.mass_hole_conduction.size();i++) { oss << xOUT.mass_hole_conduction.at(i) << " "; } oss << endl;
  oss << " conduction_band_min.size()=" << xOUT.conduction_band_min.size() << ": "; for(uint i=0;i<xOUT.conduction_band_min.size();i++) { oss << xOUT.conduction_band_min.at(i) << " "; } oss << endl;
  oss << " conduction_band_min_net=" << xOUT.conduction_band_min_net << endl;
  oss << " valence_band_max.size()=" << xOUT.valence_band_max.size() << ": "; for(uint i=0;i<xOUT.valence_band_max.size();i++) { oss << xOUT.valence_band_max.at(i) << " "; } oss << endl;
  oss << " valence_band_max_net=" << xOUT.valence_band_max_net << endl;
  oss << " Egap.size()=" << xOUT.Egap.size() << ": "; for(uint i=0;i<xOUT.Egap.size();i++) { oss << xOUT.Egap.at(i) << " "; } oss << endl;
  oss << " Egap_net=" << xOUT.Egap_net << endl;
  oss << " Egap_fit.size()=" << xOUT.Egap_fit.size() << ": "; for(uint i=0;i<xOUT.Egap_fit.size();i++) { oss << xOUT.Egap_fit.at(i) << " "; } oss << endl;
  oss << " Egap_fit_net=" << xOUT.Egap_fit_net << endl;
  oss << " Egap_type.size()=" << xOUT.Egap_type.size() << ": "; for(uint i=0;i<xOUT.Egap_type.size();i++) { oss << xOUT.Egap_type.at(i) << " "; } oss << endl;
  oss << " Egap_type_net=" << xOUT.Egap_type_net << endl;
  //[CO20200404 - OBSOLETE]oss << " ERROR=" << xOUT.ERROR << endl;
  // ----------------------------------------------------------------------
  if(LDEBUG) cerr << XPID << "xOUTCAR::operator<<: END (" << time_delay(seconds) << ")" << endl;
  return oss;   
} //SC20200330

bool xOUTCAR::GetProperties(const string& stringIN,bool QUIET) {
  stringstream sss; sss.str(stringIN);
  if(filename=="") filename="string";
  return xOUTCAR::GetProperties(sss,QUIET);
}

bool xOUTCAR::GetPropertiesFile(const string& fileIN,bool QUIET) {
  stringstream sss;
  aurostd::efile2stringstream(fileIN,sss);
  //[CO20210315 - always set for GetPropertiesFile]if(filename=="") 
  filename=fileIN;
  return xOUTCAR::GetProperties(sss,QUIET);
}

bool xOUTCAR::GetPropertiesFile(const string& fileIN,uint natoms_check,bool QUIET) {
  bool flag=GetPropertiesFile(fileIN,QUIET);
  if(aurostd::abs(natoms_check-(double) natoms)>0.1) { 
    //[CO20200404 - OBSOLETE]cerr << "ERROR xOUTCAR::GetPropertiesFile: natoms_check(" << natoms_check << ")!= (int) natoms(" << natoms << ") ..." << endl;
    stringstream message;
    message << "natoms_check(" << natoms_check << ")!= (int) natoms(" << natoms << ")";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _VALUE_ERROR_);
  }
  return flag;
}

bool xOUTCAR::GetPropertiesUrlFile(const string& url,const string& file,bool QUIET) {
  //[CO20200502 - OBSOLETE]string tmpfile=XHOST.tmpfs+"/_aflow_"+XHOST.user+".pid"+XHOST.ostrPID.str()+".a"+string(AFLOW_VERSION)+".rnd"+aurostd::utype2string(uint((double) std::floor((double)100000*aurostd::ran0())))+".u"+aurostd::utype2string(uint((double) aurostd::get_useconds()))+"_"+file;
  string tmpfile=aurostd::TmpFileCreate("xOUTCAR_GetProperties"); //CO20200502 - threadID
  aurostd::url2file(url+"/"+file,tmpfile,!QUIET);
  bool out=GetPropertiesFile(tmpfile,QUIET); //CO20200404 - added QUIET
  filename="url="+url;  //CO20210315
  aurostd::RemoveFile(tmpfile);
  return out;
}

vector<string> xOUTCAR::GetCorrectPositions(const string& line,uint expected_count) {
  //FIRST FIX : NEGATIVE SIGN
  //this function should fix the last line
  //   volume of cell :      369.80
  //     direct lattice vectors                 reciprocal lattice vectors
  //   -2.730747137  2.730747137 12.397646334     0.000000000  0.183100073  0.040330236
  //    2.730747137 -2.730747137 12.397646334     0.183100073  0.000000000  0.040330236
  //    2.730747137  2.730747137-12.397646334     0.183100073  0.183100073  0.000000000
  //
  //SECOND FIX: LARGE NUMBERS
  //this function should fix the last line
  //     direct lattice vectors                 reciprocal lattice vectors
  //    2.156793936 -3.735676679  0.000000000     0.231825578 -0.133844560  0.000000000
  //    2.156793936  3.735676679  0.000000000     0.231825578  0.133844560  0.000000000
  //    0.000000000  0.000000000109.286277550     0.000000000  0.000000000  0.009150280
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  vector<string> tokens;
  aurostd::string2tokens(line,tokens);
  if(tokens.size()==expected_count){return tokens;}
  if(0||LDEBUG) cerr << __AFLOW_FUNC__ << " issuing fix for bad lattice vectors (negative sign) on this line: " << line << endl;
  //try fixing negative sign first
  vector<string> neg_tokens,_tokens;
  std::size_t pos;
  bool first_number_negative;
  for(uint i=0;i<tokens.size();i++){
    first_number_negative=false;
    pos=tokens[i].find('-',0);
    if(pos!=std::string::npos){
      if(pos==0){first_number_negative=true;}
      aurostd::string2tokens(tokens[i],_tokens,"-");
      neg_tokens.push_back(( first_number_negative ? "-" : "" )+_tokens[0]);
      for(uint ii=1;ii<_tokens.size();ii++){neg_tokens.push_back("-" + _tokens[ii]);}
    } else {neg_tokens.push_back(tokens[i]);}
  }
  if(neg_tokens.size()==expected_count){return neg_tokens;}
  //negative sign fix not enough, now  we need to look at large number problems
  if(0||LDEBUG) cerr << __AFLOW_FUNC__ << " issuing fix for bad lattice vectors (large numbers) on this line: " << line << endl;
  vector<string> dec_tokens;
  string good_num;
  _tokens.clear();
  //first, let's find the precision, if we can
  for(uint i=0;i<neg_tokens.size()&&good_num.empty();i++){
    aurostd::string2tokens(neg_tokens[i],_tokens,".");
    if(_tokens.size()==2){good_num=neg_tokens[i];}  //this number is ok
  }
  if(good_num.empty()){dec_tokens.clear();return dec_tokens;} //we have no hope
  //_tokens contains precision
  uint precision=_tokens[1].size();
  if(0||LDEBUG) cerr << __AFLOW_FUNC__ << " precision of numbers on this line: " << precision << endl;
  //let's build numbers with the right count of digits
  string num;
  bool fidelity=true;
  for(uint i=0;i<neg_tokens.size()&&fidelity;i++){
    aurostd::string2tokens(neg_tokens[i],_tokens,".");
    if(_tokens.size()==2){dec_tokens.push_back(neg_tokens[i]);continue;}  //this number is ok
    //we need to parse this number per precision
    num=_tokens[0];
    for(uint j=1;j<_tokens.size()-1&&fidelity;j++){
      if(_tokens[j].size()<precision+1){fidelity=false;break;}  //should be at least as big as precision + 1 (full number after "." of first number + 0 of next number)
      num+=".";
      for(uint k=0;k<precision;k++){
        if(!isdigit(_tokens[j][k])){fidelity=false;break;}  //check every character is digit
        num+=_tokens[j][k];
      }
      dec_tokens.push_back(num);
      num="";
      for(uint k=precision;k<_tokens[j].size();k++){
        if(!isdigit(_tokens[j][k])){fidelity=false;break;}  //check every character is digit, we took care of negative signs already
        num+=_tokens[j][k];
      }
    }
    num+=".";
    uint j=_tokens.size()-1;
    if(_tokens[j].size()!=precision){fidelity=false;break;}
    for(uint k=0;k<_tokens[j].size()&&fidelity;k++){
      if(!isdigit(_tokens[j][k])){fidelity=false;break;}  //check every character is digit
      num+=_tokens[j][k];
    }
    dec_tokens.push_back(num);
  }
  //run through all numbers just to make sure they all have the right precision
  for(uint i=0;i<dec_tokens.size()&&fidelity;i++){
    aurostd::string2tokens(dec_tokens[i],_tokens,".");
    if(_tokens.size()!=2 || _tokens[1].size()!=precision){fidelity=false;break;}
  }
  if(!fidelity || dec_tokens.size()!=expected_count){dec_tokens.clear();return dec_tokens;}
  if(0||LDEBUG){ 
    cerr << __AFLOW_FUNC__ << " repaired vector: ";
    for(uint i=0;i<dec_tokens.size();i++){cerr << dec_tokens[i] << " ";}
    cerr << endl;
  }
  return dec_tokens;
}


bool xOUTCAR::GetProperties(const stringstream& stringstreamIN,bool QUIET) {
  bool LDEBUG=(FALSE || XHOST.DEBUG || !QUIET);
  stringstream message;
  bool force_exit=XHOST.POSTPROCESS; //SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only

  bool ERROR_flag=FALSE;
  clear(); // so it does not mess up vector/deque

  long double seconds=aurostd::get_seconds();
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN (" << time_delay(seconds) << ")" << endl;
  stringstream sss; sss.str(stringstreamIN.str());
  content=stringstreamIN.str();
  vcontent.clear();
  vector<string> vline,tokens,vcontentRED;
  aurostd::string2vectorstring(content,vcontent);
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vcontent.size()=" << vcontent.size() << " (" << time_delay(seconds) << ")" << endl;
  for(uint iline=0;iline<vcontent.size();iline++) { // NEW - FROM THE BACK
    string saus=vcontent.at(iline);
    if(aurostd::substring2bool(saus,"="))  vcontentRED.push_back(saus);
    if(aurostd::substring2bool(saus,"total")) vcontentRED.push_back(saus); // total
    if(aurostd::substring2bool(saus,"number")) vcontentRED.push_back(saus); // number
    if(aurostd::substring2bool(saus,"bands")) vcontentRED.push_back(saus); // bands
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vcontentRED.size()=" << vcontentRED.size() << " (" << time_delay(seconds) << ")" << endl;

  string line;
  if(filename=="") filename="stringstream";
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " filename=[" << filename << "]" << endl;
  if(LDEBUG) cerr.precision(12);

  // ----------------------------------------------------------------------
  //SC
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;

  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD SYSTEM (" << time_delay(seconds) << ")" << endl;
  SYSTEM="";
  line="";
  for(uint iline=0;iline<vcontent.size();iline++) { // NEW - FROM THE BACK 
    if(aurostd::substring2bool(vcontent.at(iline),"SYSTEM")) // VASP
      if(aurostd::substring2bool(vcontent.at(iline),"=")) { // VASP
        line=vcontent.at(iline);
        break;
      }
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " line=" << line << endl;
  if(line!="") {
    aurostd::string2tokens(line,tokens,"="); // cerr << tokens.size() << endl;
    SYSTEM=aurostd::RemoveWhiteSpaces(tokens.at(1));
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " SYSTEM=" << SYSTEM << endl; 
  // ----------------------------------------------------------------------
  //KY STUFF DONT TOUCH - GET Number of IONS and Fermi
  // ----------------------------------------------------------------------
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD KESONG STUFF (" << time_delay(seconds) << ")" << endl;
  string anchor_word_NELM="NELM"; //CO20200624
  string anchor_word_NIONS="NIONS";
  string anchor_word_LSORBIT="LSORBIT =";
  string anchor_word_Efermi="E-fermi";
  string anchor_word_EFIELD_PEAD="EFIELD_PEAD"; //CO20210315
  string tmp="";
  efield_pead.resize(3);  //CO20210315
  while(getline(sss,line)) {
    if(line.find(anchor_word_NELM) !=string::npos) {  //CO20200624
      aurostd::string2tokens(line,tokens,";");
      if(tokens.size()<1){continue;}
      tmp=tokens[0];
      aurostd::string2tokens(tmp,tokens,"=");
      if(tokens.size()<2){continue;}
      tmp=tokens[1];
      NELM=aurostd::string2utype<int>(aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tmp));
    }
    if(line.find(anchor_word_NIONS) !=string::npos) {
      aurostd::string2tokens(line,tokens," ");
      NIONS=aurostd::string2utype<int>(tokens.at(tokens.size()-1)); //last one
    }
    if(line.find(anchor_word_LSORBIT) !=string::npos) {
      aurostd::string2tokens(line,tokens," ");
      string LSlabel=tokens.at(2);
      if(LSlabel == "T") isLSCOUPLING=TRUE;
      if(LSlabel == "F") isLSCOUPLING=FALSE;
    }
    if(line.find(anchor_word_Efermi) !=string::npos) {
      aurostd::string2tokens(line,tokens," ");
      Efermi=aurostd::string2utype<double>(tokens.at(2));
    }
    if(line.find(anchor_word_EFIELD_PEAD) !=string::npos) { //CO20210315 - moved from aflow_ivasp XVASP_INCAR_EFIELD_PEAD()
      aurostd::string2tokens(line,tokens,"=");
      if(tokens.size()<2){continue;}
      tmp=tokens[1];
      aurostd::string2tokens(tmp,tokens);
      if(tokens.size()<3){continue;}
      for(int i=0;i<3;i++){
        efield_pead[i+1]=aurostd::string2utype<double>(tokens[i]);
      }
    }
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " NELM=" << NELM << endl; //CO20200624
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " NIONS=" << NIONS << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " isLSCOUPLING=" << isLSCOUPLING << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " Efermi=" << Efermi << endl;
  // return TRUE;
  // ----------------------------------------------------------------------
  // DO SC STUFF
  // ----------------------------------------------------------------------
  natoms=(double) NIONS;
  // ----------------------------------------------------------------------
  // LOAD EENTROPY DATA
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD EENTROPY DATA (" << time_delay(seconds) << ")" << endl;
  for(int iline=(int)vcontent.size()-1;iline>=0;iline--)  // NEW - FROM THE BACK 
    if(aurostd::substring2bool(vcontent.at(iline),"entropy"))
      if(aurostd::substring2bool(vcontent.at(iline),"EENTRO"))
        if(aurostd::substring2bool(vcontent.at(iline),"T*S")) {
          line=vcontent.at(iline);
          break;
        }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " line=" << line << endl; 
  aurostd::string2tokens(line,tokens,"=");
  if(tokens.size()==0) {
    //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - "<< __AFLOW_FUNC__ << " wrong number of entries (entropy) in OUTCAR; line=[ " << line << "]" << "   filename=[" << filename << "]" << endl;
    message << "Wrong number of entries (entropy) in OUTCAR; line=[ " << line << "]" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  }
  eentropy_cell=0;
  if(tokens.size()>1) eentropy_cell=aurostd::string2utype<double>(tokens.at(1));
  eentropy_atom=eentropy_cell/natoms;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " eentropy_cell=" << eentropy_cell << endl;
  // LOAD ENERGY DATA  // without energy of electrons
  // CO20211109 - notes about which energy to pick
  // we report the energy without entropy (0K energy).
  // the sigma->0 is an extrapolation calculation.
  // the only way to know the energy at sigma=0 is to do the calculation at
  // 0, but then the system would not converge.
  // it does not matter much which energy to pick: raw energy values of VASP are arbitrary.
  // focus on delta energy (just be consistent).
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD ENERGY DATA" << endl;
  line="";
  for(int iline=(int)vcontent.size()-1;iline>=0;iline--)  // NEW - FROM THE BACK 
    if(aurostd::substring2bool(vcontent.at(iline),"energy")) // VASP
      if(aurostd::substring2bool(vcontent.at(iline),"without")) // VASP
        if(aurostd::substring2bool(vcontent.at(iline),"entropy"))   {
          line=vcontent.at(iline);
          break;
        }  // VASP
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " line=" << line << endl;
  aurostd::string2tokens(line,tokens,"=");
  if(tokens.size()==0) {
    //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - "<< __AFLOW_FUNC__ << " wrong number of entries (energy_1) in OUTCAR; line=[ " << line << "]" << "   filename=[" << filename << "]" << endl;
    message << "Wrong number of entries (energy_1) in OUTCAR; line=[ " << line << "]" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  }
  if(tokens.size()>1) line=tokens.at(1);
  aurostd::string2tokens(line,tokens,"e");
  if(tokens.size()==0) {
    //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - "<< __AFLOW_FUNC__ << " wrong number of entries (energy_2) in OUTCAR; line=[ " << line << "]" << "   filename=[" << filename << "]" << endl;
    message << "Wrong number of entries (energy_2) in OUTCAR; line=[ " << line << "]" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  }
  if(tokens.size()>0) energy_cell=aurostd::string2utype<double>(tokens.at(0));
  energy_atom=energy_cell/natoms;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " energy_cell=" << energy_cell << " energy_atom=" << energy_atom << endl; 
  // ----------------------------------------------------------------------
  // LOAD PV DATA (IF PRESENT)
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD PV DATA (" << time_delay(seconds) << ")" << endl;
  line="";
  for(int iline=(int)vcontent.size()-1;iline>=0;iline--)  // NEW - FROM THE BACK 
    if(aurostd::substring2bool(vcontent.at(iline),"TOTEN"))  // VASP
      if(aurostd::substring2bool(vcontent.at(iline),"P"))  // VASP
        if(aurostd::substring2bool(vcontent.at(iline),"V"))  // VASP
          if(!aurostd::substring2bool(vcontent.at(iline),"VPU"))    // SOME VPU/CPU stuff ORTHCH  // VASP
            if(!aurostd::substring2bool(vcontent.at(iline),"CPU")) { // SOME VPU/CPU stuff ORTHCH  // VASP
              line=vcontent.at(iline);
              break;
            } 
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " line=" << line << endl;
  if(line.length()<10) {
    PV_cell=0.0;
    PV_atom=PV_cell/natoms;
  } else {
    aurostd::string2tokens(line,tokens,"=");
    if(tokens.size()!=3) {
      //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - "<< __AFLOW_FUNC__ << " wrong number of entries (PV) in OUTCAR; line=[ " << line << "]" << "   filename=[" << filename << "]" << endl;
      message << "Wrong number of entries (PV) in OUTCAR; line=[ " << line << "]" << "   filename=[" << filename << "]";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
      ERROR_flag=TRUE;
    }
    if(tokens.size()>2) PV_cell=aurostd::string2utype<double>(tokens.at(2));
    PV_atom=PV_cell/natoms;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " PV_cell=" << PV_cell << " PV_atom=" << PV_atom << endl; 
  }

  // correct if PV

  if(abs(PV_atom)<PRESSURE_ZERO_ENTHALPY_ENERGY) {
    enthalpy_cell=energy_cell;  // default without PV
    enthalpy_atom=energy_atom;   // default without PV
  } else {
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " pressure=" << pressure << endl; 
    // LOAD ENTHALPY DATA IF P!=0
    line="";
    for(int iline=(int)vcontent.size()-1;iline>=0;iline--)  // NEW
      if(aurostd::substring2bool(vcontent.at(iline),"TOTEN")) // VASP
        if(aurostd::substring2bool(vcontent.at(iline),"enthalpy")) { // VASP
          line=vcontent.at(iline);
          break;
        } 
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " line=" << line << endl;
    aurostd::string2tokens(line,tokens," ");

    if(tokens.size()!=0) {
      for(int i=tokens.size()-2;i>=0;i--)
        //    for(uint i=0;i<tokens.size();i++)
      { //CO20200106 - patching for auto-indenting
        if(aurostd::substring2bool(tokens.at(i),"=")) {
          if(aurostd::substring2bool(tokens.at(i-1),"TOTEN")) {
            enthalpy_cell=aurostd::string2utype<double>(tokens.at(i+1));
            enthalpy_atom=enthalpy_cell/natoms;
            if(LDEBUG) cerr << __AFLOW_FUNC__ << " enthalpy_cell=" << enthalpy_cell << " enthalpy_atom=" << enthalpy_atom << endl; 
          }
        }
      }
    }
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " enthalpy_cell=" << enthalpy_cell << " enthalpy_atom=" << enthalpy_atom << endl; 
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " energy_cell=" << energy_cell << " energy_atom=" << energy_atom << endl; 
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " PV_cell=" << PV_cell << " PV_atom=" << PV_atom << endl; 

  // ----------------------------------------------------------------------
  // LOAD PVstress DATA (IF PRESENT)
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD STRESS DATA (" << time_delay(seconds) << ")" << endl;
  line="";
  for(int iline=(int)vcontent.size()-1;iline>=0;iline--)  // NEW - FROM THE BACK 
    if(aurostd::substring2bool(vcontent.at(iline),"external"))  // VASP
      if(aurostd::substring2bool(vcontent.at(iline),"pressure"))  // VASP
        if(iline>1)
          if(aurostd::substring2bool(vcontent.at(iline-1),"kB")) { // VASP
            line=vcontent.at(iline-1);
            break;
          }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " line=" << line << endl;
  stress.clear();
  aurostd::string2tokens(line,tokens," ");
  if(tokens.size()>=8) {
    stress(1,1)=aurostd::string2utype<double>(tokens.at(2));
    stress(2,2)=aurostd::string2utype<double>(tokens.at(3));
    stress(3,3)=aurostd::string2utype<double>(tokens.at(4));
    stress(1,2)=aurostd::string2utype<double>(tokens.at(5));
    stress(2,1)=aurostd::string2utype<double>(tokens.at(5));
    stress(2,3)=aurostd::string2utype<double>(tokens.at(6));
    stress(3,2)=aurostd::string2utype<double>(tokens.at(6));
    stress(1,3)=aurostd::string2utype<double>(tokens.at(7));
    stress(3,1)=aurostd::string2utype<double>(tokens.at(7));
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " stress=" << endl << stress << endl; 
  //   Direction    XX        YY        ZX        XY       YZ       ZX
  //   in kB       -5.93     -5.93    -19.31      0.00      0.00     -0.00
  // ----------------------------------------------------------------------
  // LOAD VOLUME DATA (IF PRESENT)
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD volume DATA (" << time_delay(seconds) << ")" << endl;
  line="";
  for(int iline=(int)vcontent.size()-1;iline>=0;iline--)  // NEW FROM THE BACK
    if(aurostd::substring2bool(vcontent.at(iline),"volume"))
      if(aurostd::substring2bool(vcontent.at(iline),"of"))
        if(aurostd::substring2bool(vcontent.at(iline),"cell")) {
          line=vcontent.at(iline);
          break;
        } 
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " line=" << line << endl; 
  volume_cell=0.0;
  aurostd::string2tokens(line,tokens," ");
  if(tokens.size()>=5) {
    if(tokens.at(0)=="volume") {
      volume_cell=aurostd::string2utype<double>(tokens.at(tokens.size()-1)); 
    }
  }
  volume_atom=volume_cell/natoms;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " volume_cell=" << volume_cell << " volume_atom=" << volume_atom << endl;   
  // ----------------------------------------------------------------------
  // LOAD PRESSURE_RESIDUAL/PULAY DATA
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD PRESSURE_RESIDUAL/PULAY DATA (" << time_delay(seconds) << ")" << endl;
  line="";
  for(int iline=(int)vcontent.size()-1;iline>=0;iline--)  // NEW - FROM THE BACK 
    if(aurostd::substring2bool(vcontent.at(iline),"Pullay") || aurostd::substring2bool(vcontent.at(iline),"pullay") ||
        aurostd::substring2bool(vcontent.at(iline),"Pulay") || aurostd::substring2bool(vcontent.at(iline),"pulay"))  // VASP
      if(aurostd::substring2bool(vcontent.at(iline),"stress")) { // VASP
        line=vcontent.at(iline);
        break;
      } 
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " line=" << line << endl;
  aurostd::string2tokens(line,tokens,"=");
  if(tokens.size()==0) {
    //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - "<< __AFLOW_FUNC__ << " wrong number of entries (Pulay stress) in OUTCAR; line=[ " << line << "]" << "   filename=[" << filename << "]" << endl;
    message << "Wrong number of entries (Pulay stress) in OUTCAR; line=[ " << line << "]" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  }

  pressure_residual=0;
  Pulay_stress=0;
  if(tokens.size()>0) {
    string line_Pulay_stress=tokens.at(tokens.size()-1);
    string line_pressure_residual=tokens.at(tokens.size()-2);
    // PRESSURE
    pressure_residual=0.0;
    aurostd::string2tokens(line_pressure_residual,tokens," ");
    if(tokens.size()>0) pressure_residual=aurostd::string2utype<double>(tokens.at(0));
    //   if(pressure_residual<0.0 && pressure_residual >-10.0) pressure_residual=0.0; // TITANIC PATCH
    // PULAY
    Pulay_stress=0.0;
    aurostd::string2tokens(line_Pulay_stress,tokens," ");
    if(tokens.size()>0) Pulay_stress=aurostd::string2utype<double>(tokens.at(0));
  }

  if(LDEBUG) cerr << __AFLOW_FUNC__ << " pressure_residual=" << pressure_residual << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " Pulay_stress=" << Pulay_stress << endl;

  // ----------------------------------------------------------------------
  // LOAD SPIN
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD SPIN DATA (" << time_delay(seconds) << ")" << endl;
  line="";
  for(int iline=(int)vcontent.size()-1;iline>=0;iline--)  // NEW FROM THE BACK
    if(aurostd::substring2bool(vcontent.at(iline),"magnetization"))
      if(aurostd::substring2bool(vcontent.at(iline),"number")) {
        line=vcontent.at(iline);
        break;
      } 
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " line=" << line << endl; 
  mag_cell=0.0;
  aurostd::string2tokens(line,tokens," ");
  if(tokens.size()>=6) {
    if(tokens.at(0)=="number" && tokens.at(4)=="magnetization" && tokens.size()==6) 
      mag_cell=aurostd::string2utype<double>(tokens.at(5)); //JX seems to work
  }
  mag_atom=mag_cell/natoms;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " mag_cell=" << mag_cell << " mag_atom=" << mag_atom << endl;   
  // ----------------------------------------------------------------------
  // LOAD SPIN DECOMPOSITION
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD SPIN DECOMPOSITION DATA" << endl;
  line="";
  uint mline=0;
  //DX20171205 - Check for magnetization z (non-collinear) - START
  vector<double> vmag_z;
  bool found_magz_line=false;
  for(int iline=(int)vcontent.size()-1;iline>=0;iline--)  // NEW FROM THE BACK
    if(aurostd::substring2bool(vcontent.at(iline),"magnetization"))
      if(aurostd::substring2bool(vcontent.at(iline),"(z)")) {
        mline=iline+4;
        found_magz_line=true;
        break;
      } 
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " magnetization (z) line=" << line << endl; 
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " magnetization (z) mline=" << mline << endl; 
  if(found_magz_line) {
    for(uint iline=mline;iline<mline+natoms;iline++) {
      aurostd::string2tokens(vcontent.at(iline),tokens," ");
      if(tokens.size()>=5)
        vmag_z.push_back(aurostd::string2utype<double>(tokens.at(tokens.size()-1)));
    }
  }
  //DX20171205 - Check for magnetization z (non-collinear) - END
  //DX20171205 - Check for magnetization y (non-collinear) - START
  vector<double> vmag_y;
  bool found_magy_line=false;
  for(int iline=(int)vcontent.size()-1;iline>=0;iline--)  // NEW FROM THE BACK
    if(aurostd::substring2bool(vcontent.at(iline),"magnetization"))
      if(aurostd::substring2bool(vcontent.at(iline),"(y)")) {
        mline=iline+4;
        found_magy_line=true;
        break;
      } 
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " magnetization (y) line=" << line << endl; 
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " magnetization (y) mline=" << mline << endl; 
  if(found_magy_line) {
    for(uint iline=mline;iline<mline+natoms;iline++) {
      aurostd::string2tokens(vcontent.at(iline),tokens," ");
      if(tokens.size()>=5)
        vmag_y.push_back(aurostd::string2utype<double>(tokens.at(tokens.size()-1)));
    }
  }
  //DX20171205 - Check for magnetization y (non-collinear) - END
  vector<double> vmag_x;  //DX20171205
  bool found_magx_line=false;
  for(int iline=(int)vcontent.size()-1;iline>=0;iline--)  // NEW FROM THE BACK
    if(aurostd::substring2bool(vcontent.at(iline),"magnetization"))
      if(aurostd::substring2bool(vcontent.at(iline),"(x)")) {
        mline=iline+4;
        found_magx_line=true;
        break;
      } 
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " magnetization (x) line=" << line << endl; 
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " magnetization (x) mline=" << mline << endl; 
  if(found_magx_line) {
    for(uint iline=mline;iline<mline+natoms;iline++) {
      aurostd::string2tokens(vcontent.at(iline),tokens," ");
      if(tokens.size()>=5)
        vmag_x.push_back(aurostd::string2utype<double>(tokens.at(tokens.size()-1))); //DX20171205 vmag to vmag_x
    }
  }
  //DX20171205 - Non-collinear vs collinear - START
  if(found_magx_line && found_magy_line && found_magz_line){ //non-collinear calculation
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " non-collinear magnetization found." << endl; 
    if(vmag_x.size()!=vmag_y.size() || vmag_x.size()!=vmag_z.size()){
      //[CO20200404 - OBSOLETE]if(!QUIET){cerr << "WARNING - "<< __AFLOW_FUNC__ << " number of magnetization components (x, y, z) are not the same in OUTCAR; filename=[" << filename << "]" << endl;}
      message << "Number of magnetization components (x, y, z) are not the same in OUTCAR; filename=[" << filename << "]";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
      ERROR_flag=TRUE;
    }
    for(uint m=0;m<vmag_x.size();m++){
      xvector<double> mag_xyz; mag_xyz(1)=vmag_x[m]; mag_xyz(2)=vmag_y[m]; mag_xyz(3)=vmag_z[m];
      vmag_noncoll.push_back(mag_xyz);
    }
  }
  else if(found_magx_line && !found_magy_line && !found_magz_line){ //collinear calculation (x only)
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " collinear magnetization found." << endl; 
    for(uint m=0;m<vmag_x.size();m++){
      vmag.push_back(vmag_x[m]);
    }
  }
  //  if(LDEBUG) cerr << __AFLOW_FUNC__ << " mag_cell=" << mag_cell << " mag_atom=" << mag_atom << endl;   
  // ----------------------------------------------------------------------
  // LOAD FORCES/POSITIONS
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD FORCES/POSITIONS DATA (" << time_delay(seconds) << ")" << endl;
  vforces.clear();                    // QM FORCES calculation
  vpositions_cartesian.clear();                 // QM POSITIONS calculation
  for(uint i=0;i<(uint) natoms;i++) {
    xvector<double> force(3);
    vforces.push_back(force);
    xvector<double> position(3);
    vpositions_cartesian.push_back(position);
  }
  for(int iline=(int)vcontent.size()-1;iline>=0;iline--) {
    if(aurostd::substring2bool(vcontent.at(iline),"TOTAL-FORCE (eV/Angst)")) {
      for(uint iat=0;iat<(uint)natoms && iat<vcontent.size();iat++) {
        aurostd::string2tokens(vcontent.at(iline+iat+2),tokens," ");
        //cerr << natoms << " " << tokens.size() << "    " << vcontent.at(iline+iat+2) << endl;
        if(tokens.size()<6) {
          //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - "<< __AFLOW_FUNC__ << " wrong number of force/positions entries in OUTCAR; line=[ " << vcontent.at(iline+iat+2) << "]" << "   filename=[" << filename << "]" << endl;
          message << "Wrong number of force/positions entries in OUTCAR; line=[ " << vcontent.at(iline+iat+2) << "]" << "   filename=[" << filename << "]";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
          ERROR_flag=TRUE;
        }
        vpositions_cartesian.at(iat)[1]=aurostd::string2utype<double>(tokens.at(0));
        vpositions_cartesian.at(iat)[2]=aurostd::string2utype<double>(tokens.at(1));
        vpositions_cartesian.at(iat)[3]=aurostd::string2utype<double>(tokens.at(2));
        vforces.at(iat)[1]=aurostd::string2utype<double>(tokens.at(3));
        vforces.at(iat)[2]=aurostd::string2utype<double>(tokens.at(4));
        vforces.at(iat)[3]=aurostd::string2utype<double>(tokens.at(5));
      }
      iline=-1;
    }
  }
  // ----------------------------------------------------------------------
  // ENCUT EDIFF EDIFFG POTIM TEIN TEBEG TEEND SMASS NPACO APACO PSTRESS  
  // NBANDS NKPTS NSW NBLOCK KBLOCK IBRION NFREE ISIF IWAVPR ISYM TEIN TEBEG ISPIN
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD \"Electronic convergence/Ionic relaxation\" DATA (" << time_delay(seconds) << ")" << endl;
  vline.clear();
  ENCUT=0.0;EDIFF=0.0;EDIFFG=0.0;POTIM=0.0;TEIN=0.0;TEBEG=0.0;TEEND=0.0;SMASS=0.0;NPACO=0.0;APACO=0.0;PSTRESS=0.0;pressure=0.0;     // 
  NBANDS=0;NKPTS=0;NSW=0;NBLOCK=0;KBLOCK=0;IBRION=0;NFREE=0;ISIF=0;IWAVPR=0;ISYM=0;TEIN=0;TEBEG=0;ISPIN=0; //
  total_energy_change=0.0;
  //for(uint iline=vcontent.size()-1;iline>0;iline--)
  for(uint iline=0;iline<vcontentRED.size();iline++)
  { //CO20200106 - patching for auto-indenting
    if(aurostd::substring2bool(vcontentRED.at(iline),"ENCUT") && aurostd::substring2bool(vcontentRED.at(iline),"eV")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"EDIFF") && aurostd::substring2bool(vcontentRED.at(iline),"stopping-criterion for ELM")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"EDIFFG") && aurostd::substring2bool(vcontentRED.at(iline),"stopping-criterion for IOM")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"NSW") && aurostd::substring2bool(vcontentRED.at(iline),"number of steps for IOM")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"NBLOCK") && aurostd::substring2bool(vcontentRED.at(iline),"inner block")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"KBLOCK") && aurostd::substring2bool(vcontentRED.at(iline),"outer block")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"IBRION") && aurostd::substring2bool(vcontentRED.at(iline),"ionic relax")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"NFREE") && aurostd::substring2bool(vcontentRED.at(iline),"steps in history")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"ISIF") && aurostd::substring2bool(vcontentRED.at(iline),"stress and relaxation")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"IWAVPR") && aurostd::substring2bool(vcontentRED.at(iline),"prediction")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"ISYM") && aurostd::substring2bool(vcontentRED.at(iline),"nonsym")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"POTIM") && aurostd::substring2bool(vcontentRED.at(iline),"time-step")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"TEIN") && aurostd::substring2bool(vcontentRED.at(iline),"initial temperature")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"TEBEG") && aurostd::substring2bool(vcontentRED.at(iline),"temperature during run")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"TEEND") && aurostd::substring2bool(vcontentRED.at(iline),"temperature during run")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"SMASS") && aurostd::substring2bool(vcontentRED.at(iline),"Nose mass-parameter")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"NPACO") && aurostd::substring2bool(vcontentRED.at(iline),"distance and # of slots")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"APACO") && aurostd::substring2bool(vcontentRED.at(iline),"distance and # of slots")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"PSTRESS") && (aurostd::substring2bool(vcontentRED.at(iline),"pullay stress") || aurostd::substring2bool(vcontentRED.at(iline),"pulay stress"))) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"ISPIN") && aurostd::substring2bool(vcontentRED.at(iline),"spin polarized calculation?")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"NKPTS") && aurostd::substring2bool(vcontentRED.at(iline),"k-Points")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"NBANDS") && aurostd::substring2bool(vcontentRED.at(iline),"number of bands")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"total energy-change") && aurostd::substring2bool(vcontentRED.at(iline),"(2. order)")) vline.push_back(vcontentRED.at(iline));
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vline.size()=" << vline.size() << " (" << time_delay(seconds) << ")" << endl;
  if(vline.size()==0) {
    //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - "<< __AFLOW_FUNC__ << " wrong number of \"Ionic relaxation\" in OUTCAR" << "   filename=[" << filename << "]" << endl;
    message << "Wrong number of \"Ionic relaxation\" in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  }
  for(uint j=0;j<vline.size();j++) {   // to the back
    //   if(LDEBUG) cerr << __AFLOW_FUNC__ << " vline.at(" << j << ")=" << vline.at(j) << endl;
    aurostd::StringSubst(vline.at(j),"="," ");
    aurostd::StringSubst(vline.at(j),";"," ");
    aurostd::StringSubst(vline.at(j),":"," ");
    //   if(LDEBUG) cerr << vline.at(j) << endl;
    aurostd::string2tokens(vline.at(j),tokens," ");
    for(uint k=0;k<tokens.size();k++) {
      if(tokens.at(k)=="ENCUT" && k+1<tokens.size()) ENCUT=aurostd::string2utype<double>(tokens.at(k+1));
      if(tokens.at(k)=="EDIFF" && k+1<tokens.size()) EDIFF=aurostd::string2utype<double>(tokens.at(k+1));
      if(tokens.at(k)=="EDIFFG" && k+1<tokens.size()) EDIFFG=aurostd::string2utype<double>(tokens.at(k+1));
      if(tokens.at(k)=="NSW" && k+1<tokens.size()) NSW=aurostd::string2utype<int>(tokens.at(k+1));
      if(tokens.at(k)=="NBLOCK" && k+1<tokens.size()) NBLOCK=aurostd::string2utype<int>(tokens.at(k+1));
      if(tokens.at(k)=="KBLOCK" && k+1<tokens.size()) KBLOCK=aurostd::string2utype<int>(tokens.at(k+1));
      if(tokens.at(k)=="IBRION" && k+1<tokens.size()) IBRION=aurostd::string2utype<int>(tokens.at(k+1));
      if(tokens.at(k)=="NFREE" && k+1<tokens.size()) NFREE=aurostd::string2utype<int>(tokens.at(k+1));
      if(tokens.at(k)=="ISIF" && k+1<tokens.size()) ISIF=aurostd::string2utype<int>(tokens.at(k+1));
      if(tokens.at(k)=="IWAVPR" && k+1<tokens.size()) IWAVPR=aurostd::string2utype<int>(tokens.at(k+1));
      if(tokens.at(k)=="ISYM" && k+1<tokens.size()) ISYM=aurostd::string2utype<int>(tokens.at(k+1));
      if(tokens.at(k)=="POTIM" && k+1<tokens.size()) POTIM=aurostd::string2utype<double>(tokens.at(k+1));
      if(tokens.at(k)=="TEIN" && k+1<tokens.size()) TEIN=aurostd::string2utype<double>(tokens.at(k+1));
      if(tokens.at(k)=="TEBEG" && k+1<tokens.size()) TEBEG=aurostd::string2utype<double>(tokens.at(k+1));
      if(tokens.at(k)=="TEEND" && k+1<tokens.size()) TEEND=aurostd::string2utype<double>(tokens.at(k+1));
      if(tokens.at(k)=="SMASS" && k+1<tokens.size()) SMASS=aurostd::string2utype<double>(tokens.at(k+1));
      if(tokens.at(k)=="NPACO" && k+1<tokens.size()) NPACO=aurostd::string2utype<double>(tokens.at(k+1));
      if(tokens.at(k)=="APACO" && k+1<tokens.size()) APACO=aurostd::string2utype<double>(tokens.at(k+1));
      if(tokens.at(k)=="PSTRESS" && k+1<tokens.size()) {
        PSTRESS=aurostd::string2utype<double>(tokens.at(k+1)); 
        pressure=PSTRESS;
      }
      if(tokens.at(k)=="NBANDS" && k+1<tokens.size()) NBANDS=aurostd::string2utype<int>(tokens.at(k+1));
      if(tokens.at(k)=="NKPTS" && k+1<tokens.size()) NKPTS=aurostd::string2utype<int>(tokens.at(k+1));
      if(tokens.at(k)=="ISPIN" && k+1<tokens.size()) ISPIN=aurostd::string2utype<int>(tokens.at(k+1));
      if(tokens.at(k)=="order)" && k+1<tokens.size()) total_energy_change=aurostd::string2utype<double>(tokens.at(k+1));
    }
  }
  if(vline.size()){ //CO
    for(uint j=vline.size()-1;j>0;j--) {   // to the front
      aurostd::StringSubst(vline.at(j),"="," ");
      aurostd::StringSubst(vline.at(j),";"," ");
      aurostd::StringSubst(vline.at(j),":"," ");
      //   if(LDEBUG) cerr << vline.at(j) << endl;
      aurostd::string2tokens(vline.at(j),tokens," ");
      for(uint k=0;k<tokens.size();k++) {
        //    if(tokens.at(k)=="order)" && k+1<tokens.size()) total_energy_change=aurostd::string2utype<double>(tokens.at(k+1));
      }
    }
  }
  if(LDEBUG) {
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
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD \"ENMAX ENMIN POMASS ZVAL EATOM RCORE RWIGS EAUG RMAX\" DATA (" << time_delay(seconds) << ")" << endl;
  vline.clear();
  vENMAX.clear();vENMIN.clear();vPOMASS.clear();vZVAL.clear();vEATOM.clear();vRCORE.clear();vRWIGS.clear();vEAUG.clear();vRAUG.clear();vRMAX.clear();
  for(uint iline=0;iline<vcontentRED.size();iline++) {
    aurostd::StringSubst(vcontentRED.at(iline),"EMMIN","ENMIN"); // Kresseeeeeeeeeee
    if(aurostd::substring2bool(vcontentRED.at(iline),"ENMAX") && aurostd::substring2bool(vcontentRED.at(iline),"ENMIN")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"POMASS") && aurostd::substring2bool(vcontentRED.at(iline),"ZVAL")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"EATOM") && aurostd::substring2bool(vcontentRED.at(iline),"eV")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"RCORE") && aurostd::substring2bool(vcontentRED.at(iline),"radius")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"RWIGS") && aurostd::substring2bool(vcontentRED.at(iline),"radius")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"EAUG")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"RAUG") && aurostd::substring2bool(vcontentRED.at(iline),"sphere")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"RMAX") && aurostd::substring2bool(vcontentRED.at(iline),"radius")) vline.push_back(vcontentRED.at(iline));
    //AS20200528
    if(aurostd::substring2bool(vcontentRED.at(iline),"NELECT") && aurostd::substring2bool(vcontentRED.at(iline),"total number of electrons")) vline.push_back(vcontentRED.at(iline));
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vline.size()=" << vline.size() << " (" << time_delay(seconds) << ")" << endl;
  if(vline.size()==0) {
    //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - "<< __AFLOW_FUNC__ << " wrong number of \"ENMAX ENMIN POMASS ZVAL EATOM RCORE RWIGS EAUG RAUG RMAX\" in OUTCAR" << "   filename=[" << filename << "]" << endl;
    message << "Wrong number of \"ENMAX ENMIN POMASS ZVAL EATOM RCORE RWIGS EAUG RAUG RMAX\" in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  }
  for(uint j=0;j<vline.size();j++) {
    // if(LDEBUG) cerr << __AFLOW_FUNC__ << " vline.at(" << j << ")=" << vline.at(j) << endl;
    aurostd::StringSubst(vline.at(j),"="," ");
    aurostd::StringSubst(vline.at(j),";"," ");
    //    if(LDEBUG) cerr << vline.at(j) << endl;
    aurostd::string2tokens(vline.at(j),tokens," ");
    for(uint k=0;k<tokens.size();k++) {
      if(tokens.at(k)=="ENMAX" && k+1<tokens.size()) vENMAX.push_back(aurostd::string2utype<double>(tokens.at(k+1)));
      if(tokens.at(k)=="ENMIN" && k+1<tokens.size()) vENMIN.push_back(aurostd::string2utype<double>(tokens.at(k+1)));
      if(tokens.at(k)=="POMASS" && k+1<tokens.size()) vPOMASS.push_back(aurostd::string2utype<double>(tokens.at(k+1)));
      if(tokens.at(k)=="ZVAL" && k+1<tokens.size()) vZVAL.push_back(aurostd::string2utype<double>(tokens.at(k+1)));
      if(tokens.at(k)=="EATOM" && k+1<tokens.size()) vEATOM.push_back(aurostd::string2utype<double>(tokens.at(k+1)));
      if(tokens.at(k)=="RCORE" && k+1<tokens.size()) vRCORE.push_back(aurostd::string2utype<double>(tokens.at(k+1)));
      if(tokens.at(k)=="RWIGS" && k==0 && k+1<tokens.size()) vRWIGS.push_back(aurostd::string2utype<double>(tokens.at(k+1))); // pick the 1st
      if(tokens.at(k)=="EAUG" && k+1<tokens.size()) vEAUG.push_back(aurostd::string2utype<double>(tokens.at(k+1)));
      if(tokens.at(k)=="RAUG" && k+1<tokens.size()) vRAUG.push_back(aurostd::string2utype<double>(tokens.at(k+1)));
      if(tokens.at(k)=="RMAX" && k+1<tokens.size()) vRMAX.push_back(aurostd::string2utype<double>(tokens.at(k+1)));
      if(tokens.at(k)=="NELECT" && k+1<tokens.size()) nelectrons=aurostd::string2utype<int>(tokens.at(k+1));//AS20200528
    }
  }

  if(vENMIN.size()<vENMAX.size()) vENMIN.push_back(0.0);  // ENMIN MIGHT NOT BE THERE
  if(vRCORE.size()<vENMAX.size()) vRCORE.push_back(0.0);  // RCORE MIGHT NOT BE THERE
  if(vRWIGS.size()<vENMAX.size()) vRWIGS.push_back(0.0);  // RWIGS MIGHT NOT BE THERE
  if(vEAUG.size()<vENMAX.size()) vEAUG.push_back(0.0);    // EAUG MIGHT NOT BE THERE
  if(vRAUG.size()<vENMAX.size()) vRAUG.push_back(0.0);    // RAUG MIGHT NOT BE THERE
  if(vRMAX.size()<vENMAX.size()) vRMAX.push_back(0.0);    // RMAX MIGHT NOT BE THERE

  ENMAX=aurostd::max(vENMAX);
  ENMIN=aurostd::min(vENMIN);
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " ENMAX=" << ENMAX << " vENMAX.size()=" << vENMAX.size() << ": "; for(uint i=0;i<vENMAX.size();i++) { cerr << vENMAX.at(i) << " "; } cerr << " " << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " ENMIN=" << ENMIN << " vENMIN.size()=" << vENMIN.size() << ": "; for(uint i=0;i<vENMIN.size();i++) { cerr << vENMIN.at(i) << " "; } cerr << " " << endl;}

  POMASS_sum=aurostd::sum(vPOMASS);POMASS_min=aurostd::min(vPOMASS);POMASS_max=aurostd::max(vPOMASS);
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " POMASS_sum=" << POMASS_sum << " POMASS_min=" << POMASS_min << " POMASS_max=" << POMASS_max << " vPOMASS.size()=" << vPOMASS.size() << ": "; for(uint i=0;i<vPOMASS.size();i++) { cerr << vPOMASS.at(i) << " "; } cerr << " " << endl;}

  ZVAL_sum=aurostd::sum(vZVAL);ZVAL_min=aurostd::min(vZVAL);ZVAL_max=aurostd::max(vZVAL);
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " ZVAL_sum=" << ZVAL_sum << " ZVAL_min=" << ZVAL_min << " ZVAL_max=" << ZVAL_max << " vZVAL.size()=" << vZVAL.size() << ": "; for(uint i=0;i<vZVAL.size();i++) { cerr << vZVAL.at(i) << " "; } cerr << " " << endl;}

  EATOM_min=aurostd::min(vEATOM);EATOM_max=aurostd::max(vEATOM);
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " EATOM_min=" << EATOM_min << " EATOM_max=" << EATOM_max << " vEATOM.size()=" << vEATOM.size() << ": "; for(uint i=0;i<vEATOM.size();i++) { cerr << vEATOM.at(i) << " "; } cerr << " " << endl;}

  RCORE_min=aurostd::min(vRCORE);RCORE_max=aurostd::max(vRCORE);
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " RCORE_min=" << RCORE_min << " RCORE_max=" << RCORE_max << " vRCORE.size()=" << vRCORE.size() << ": "; for(uint i=0;i<vRCORE.size();i++) { cerr << vRCORE.at(i) << " "; } cerr << " " << endl;}

  RWIGS_min=aurostd::min(vRWIGS);RWIGS_max=aurostd::max(vRWIGS);
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " RWIGS_min=" << RWIGS_min << " RWIGS_max=" << RWIGS_max << " vRWIGS.size()=" << vRWIGS.size() << ": "; for(uint i=0;i<vRWIGS.size();i++) { cerr << vRWIGS.at(i) << " "; } cerr << " " << endl;}

  EAUG_min=aurostd::min(vEAUG);EAUG_max=aurostd::max(vEAUG);
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " EAUG_min=" << EAUG_min << " EAUG_max=" << EAUG_max << " vEAUG.size()=" << vEAUG.size() << ": "; for(uint i=0;i<vEAUG.size();i++) { cerr << vEAUG.at(i) << " "; } cerr << " " << endl;}

  RAUG_min=aurostd::min(vRAUG);RAUG_max=aurostd::max(vRAUG);
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " RAUG_min=" << RAUG_min << " RAUG_max=" << RAUG_max << " vRAUG.size()=" << vRAUG.size() << ": "; for(uint i=0;i<vRAUG.size();i++) { cerr << vRAUG.at(i) << " "; } cerr << " " << endl;}

  RMAX_min=aurostd::min(vRMAX);RMAX_max=aurostd::max(vRMAX);
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " RMAX_min=" << RMAX_min << " RMAX_max=" << RMAX_max << " vRMAX.size()=" << vRMAX.size() << ": "; for(uint i=0;i<vRMAX.size();i++) { cerr << vRMAX.at(i) << " "; } cerr  << " " << endl;}

  // ----------------------------------------------------------------------
  // KINETIC AND METAGGA DATA
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD KINETIC AND METAGGA DATA (" << time_delay(seconds) << ")" << endl;
  isKIN=FALSE;
  isMETAGGA=FALSE;
  METAGGA="";
  vline.clear();
  for(uint iline=0;iline<vcontent.size();iline++) {
    if(aurostd::substring2bool(vcontent.at(iline),"partial kinetic energy density read in")) { isKIN=TRUE; }
    if(aurostd::substring2bool(vcontent.at(iline),"METAGGA")) { vline.push_back(vcontent.at(iline)); }
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vline.size()=" << vline.size() << endl;

  if(vline.size()>0) {
    aurostd::string2tokens(vline.at(0),tokens,"=");
    if(tokens.size()>1) {
      aurostd::string2tokens(string(tokens.at(1)),tokens);
      if(tokens.size()>0) {
        if(tokens.at(0)!="F") {
          isMETAGGA=TRUE;
          METAGGA=tokens.at(0);
        }	
      }
    }
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " isKIN=" << isKIN << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " isMETAGGA=" << isMETAGGA << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " METAGGA=" << METAGGA << endl;

  // ----------------------------------------------------------------------
  // PSEUDOPOTENTIAL DATA
  // LDEBUG=1;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD PSEUDOPOTENTIAL DATA (" << time_delay(seconds) << ")" << endl;

  vline.clear();
  vTITEL.clear();
  vLEXCH.clear();
  for(uint iline=0;iline<vcontentRED.size();iline++)  {
    if(aurostd::substring2bool(vcontentRED.at(iline),"=")) {
      if(aurostd::substring2bool(vcontentRED.at(iline),"TITEL")) {
        aurostd::string2tokens(vcontentRED.at(iline),tokens,"=");
        if(tokens.size()>1) vTITEL.push_back(tokens.at(1));
        //	if(LDEBUG) cerr << __AFLOW_FUNC__ << " TITEL=" << tokens.at(1) << endl;
      }
      if(aurostd::substring2bool(vcontentRED.at(iline),"LEXCH")) {
        if(!aurostd::substring2bool(vcontentRED.at(iline),"exchange")) {
          aurostd::string2tokens(vcontentRED.at(iline),tokens,"=");
          if(tokens.size()>1) vLEXCH.push_back(tokens.at(1));
          //	if(LDEBUG) cerr << __AFLOW_FUNC__ << " LEXCH=" << tokens.at(1) << endl;
        }
      }
    }
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vTITEL.size()=" << vTITEL.size() << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vLEXCH.size()=" << vLEXCH.size() << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vEATOM.size()=" << vEATOM.size() << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vRMAX.size()=" << vRMAX.size() << endl;
  if(vTITEL.size()==0) {
    //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - "<< __AFLOW_FUNC__ << " wrong number of pseudopotentials (TITEL) in OUTCAR" << "   filename=[" << filename << "]" << endl;
    message << "Wrong number of pseudopotentials (TITEL) in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;  
  } //CO20200106 - patching for auto-indenting
  if(vLEXCH.size()==0) {
    //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - "<< __AFLOW_FUNC__ << " wrong number of pseudopotentials (LEXCH) in OUTCAR" << "   filename=[" << filename << "]" << endl;
    message << "Wrong number of pseudopotentials (LEXCH) in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;  
  } //CO20200106 - patching for auto-indenting
  if(vEATOM.size()==0) {
    message << "Wrong number of pseudopotentials (EATOM) in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  } //CO20200106 - patching for auto-indenting
  if(vRMAX.size()==0) {
    //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - "<< __AFLOW_FUNC__ << " wrong number of pseudopotentials (RMAX) in OUTCAR" << "   filename=[" << filename << "]" << endl;
    message << "Wrong number of pseudopotentials (RMAX) in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  }   //CO20200106 - patching for auto-indenting
  if(vTITEL.size()!=vLEXCH.size()) {
    //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - "<< __AFLOW_FUNC__ << " wrong number of pseudopotentials (TITEL/LEXCH) in OUTCAR" << "   filename=[" << filename << "]" << endl;
    message << "Wrong number of pseudopotentials (TITEL/LEXCH) in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  } //CO20200106 - patching for auto-indenting
  if(vLEXCH.size()!=vEATOM.size()) {
    //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - "<< __AFLOW_FUNC__ << " wrong number of pseudopotentials (LEXCH/EATOM) in OUTCAR" << "   filename=[" << filename << "]" << endl;
    message << "Wrong number of pseudopotentials (LEXCH/EATOM) in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  } //CO20200106 - patching for auto-indenting
  if(vEATOM.size()!=vRMAX.size()) {
    //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - "<< __AFLOW_FUNC__ << " wrong number of pseudopotentials (EATOM/RMAX) in OUTCAR" << "   filename=[" << filename << "]" << endl;
    message << "Wrong number of pseudopotentials (EATOM/RMAX) in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  }   //CO20200106 - patching for auto-indenting

  for(uint j=0;j<vTITEL.size();j++) {
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " SPECIES(" << j << ") " << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " vTITEL.at(j)=" << vTITEL.at(j) << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " vLEXCH.at(j)=" << vLEXCH.at(j) << endl;
    aurostd::string2tokens(vTITEL.at(j),tokens," ");
    pp_type=tokens.at(0);
    if(pp_type=="US" or pp_type=="NC") {
      pp_type="GGA";
      for(uint iTITEL=0;iTITEL<vcontent.size();iTITEL++) {
        if(aurostd::substring2bool(vcontent.at(iTITEL),"GGA")) {
          if(aurostd::substring2bool(vcontent.at(iTITEL),"eV")) {
            pp_type="LDA";
          }
        }
      }
    }
    if(pp_type=="PAW") pp_type="PAW_LDA"; // cerr << __AFLOW_FUNC__ << " PAW_LDA" << endl;
    if(isKIN) pp_type+="_KIN";
    if(isMETAGGA) pp_type+=":"+METAGGA;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " pp_type=" << pp_type << endl;
    species.push_back(KBIN::VASP_PseudoPotential_CleanName(tokens.at(1)));
    species_Z.push_back(xelement::symbol2Z(KBIN::VASP_PseudoPotential_CleanName(tokens.at(1))));
    species_pp.push_back(tokens.at(1));
    species_pp_type.push_back(pp_type);
    if(pp_type=="LDA" && tokens.size()<3) tokens.push_back(DEFAULT_VASP_POTCAR_DATE_POT_LDA);
    if(pp_type=="GGA" && tokens.size()<3) tokens.push_back(DEFAULT_VASP_POTCAR_DATE_POT_GGA);
    species_pp_version.push_back(tokens.at(1)+":"+pp_type+":"+tokens.at(2));


    vTITEL.at(j)=aurostd::RemoveWhiteSpaces(vTITEL.at(j));
    vLEXCH.at(j)=aurostd::RemoveWhiteSpaces(vLEXCH.at(j));   
    xPOTCAR xPOT(xPOTCAR_Finder(species_pp_AUID,species_pp_AUID_collisions,vTITEL.at(j),vLEXCH.at(j),vEATOM.at(j),vRMAX.at(j),LDEBUG)); // FIXES species_pp_AUID,species_pp_AUID_collisions
    species_pp_groundstate_energy.push_back(xPOT.species_pp_groundstate_energy.at(0));
    species_pp_groundstate_structure.push_back(xPOT.species_pp_groundstate_structure.at(0));

    if(LDEBUG) cerr << __AFLOW_FUNC__ << " SPECIES(" << j << ") [pp, type, version, AUID] = "
      << species.at(species.size()-1) << " ["
        << species_Z.at(species_Z.size()-1) << ", "
        << species_pp.at(species_pp.size()-1) << ", "
        << species_pp_type.at(species_pp_type.size()-1) << ", "
        << species_pp_version.at(species_pp_version.size()-1) << ", "
        << species_pp_AUID.at(species_pp_AUID.size()-1) <<  ", "
        << species_pp_groundstate_energy.at(species_pp_groundstate_energy.size()-1) <<  ", "
        << species_pp_groundstate_structure.at(species_pp_groundstate_structure.size()-1) << "]" << endl;
  }
  //CO20210213 - check types are all the same, if not issue warning/error (mixing is not advisable)
  for(uint i=0;i<species_pp_type.size();i++){
    if(species_pp_type[i]!=pp_type){
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,"Mismatch in species_pp_types ("+species_pp_type[i]+" vs. "+pp_type+")",*p_FileMESSAGE,*p_oss,_LOGGER_WARNING_,QUIET);
    }
  }


  if(LDEBUG) cerr << __AFLOW_FUNC__ << " PSEUDOPOTENTIAL type = " << pp_type << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " species.size()=" << species.size() << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " species_Z.size()=" << species_Z.size() << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " species_pp.size()=" << species_pp.size() << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " species_pp_type.size()=" << species_pp_type.size() << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " species_pp_version.size()=" << species_pp_version.size() << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " species_pp_AUID.size()=" << species_pp_AUID.size() << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " species_pp_AUID_collisions.size()=" << species_pp_AUID_collisions.size() << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " species_pp_groundstate_energy.size()=" << species_pp_groundstate_energy.size() << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " species_pp_groundstate_structure.size()=" << species_pp_groundstate_structure.size() << endl;
  if(species_pp_AUID_collisions.size()) {
    //[CO20200404 - OBSOLETE]cerr << __AFLOW_FUNC__ << " COLLISION species_pp_AUID_collisions.size()=" << species_pp_AUID_collisions.size() << endl;
    message << "COLLISION species_pp_AUID_collisions.size()=" << species_pp_AUID_collisions.size();
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  }

  // ----------------------------------------------------------------------
  // LDAU DATA
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD LDAU DATA (" << time_delay(seconds) << ")" << endl;
  for(uint j=0;j<species.size();j++) 
    species_pp_vLDAU.push_back(deque<double>());  // make space for LDAU

  vline.clear();
  for(uint iline=0;iline<vcontent.size();iline++) 
    if(!aurostd::substring2bool(vcontent.at(iline),"POSCAR"))
      if(!aurostd::substring2bool(vcontent.at(iline),"SYSTEM"))   
        if(aurostd::substring2bool(vcontent.at(iline),"LDAU"))
          if(!aurostd::substring2bool(vcontent.at(iline),"SYSTEM"))
            vline.push_back(vcontent.at(iline));

  int LDAUT=0;
  stringstream sdata_ldau;
  if(vline.size()!=0) {
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " LDAU calculation in OUTCAR" << endl;
    vector<int> vLDAUL;vector<double> vLDAUU,vLDAUJ;
    for(uint j=0;j<vline.size();j++) {
      aurostd::string2tokens(vline.at(j),tokens,"=");
      vline.at(j)=tokens.at(1);
    }
    for(uint j=0;j<vline.size();j++) {
      aurostd::string2tokens(vline.at(j),tokens," ");
      if(j==0) LDAUT=aurostd::string2utype<int>(vline.at(j));
      if(j==1) for(uint i=0;i<tokens.size();i++) vLDAUL.push_back(aurostd::string2utype<int>(tokens.at(i)));
      if(j==2) for(uint i=0;i<tokens.size();i++) vLDAUU.push_back(((int) 100*aurostd::string2utype<double>(tokens.at(i)))/100.0);  
      if(j==3) for(uint i=0;i<tokens.size();i++) vLDAUJ.push_back(((int) 100*aurostd::string2utype<double>(tokens.at(i)))/100.0);
    }
    if(species_pp_vLDAU.size()!=species.size()) {
      //[CO20200404 - OBSOLETE]if(!QUIET) cerr << __AFLOW_FUNC__ << " ERROR - species_pp_vLDAU.size()[" << species_pp_vLDAU.size() << "] != species.size()[" << species.size() << "]" << "   filename=[" << filename << "]" << endl;
      message << "species_pp_vLDAU.size()[" << species_pp_vLDAU.size() << "] != species.size()[" << species.size() << "]" << "   filename=[" << filename << "]";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
      ERROR_flag=TRUE;
    }
    if(species_pp_vLDAU.size()!=vLDAUL.size()) {
      //[CO20200404 - OBSOLETE]if(!QUIET) cerr << __AFLOW_FUNC__ << " ERROR - species_pp_vLDAU.size()[" << species_pp_vLDAU.size() << "] != vLDAUL.size()[" << vLDAUL.size() << "]" << "   filename=[" << filename << "]" << endl;
      message << "species_pp_vLDAU.size()[" << species_pp_vLDAU.size() << "] != vLDAUL.size()[" << vLDAUL.size() << "]" << "   filename=[" << filename << "]";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
      ERROR_flag=TRUE;
    }
    if(species_pp_vLDAU.size()!=vLDAUU.size()) {
      //[CO20200404 - OBSOLETE]if(!QUIET) cerr << __AFLOW_FUNC__ << " ERROR - species_pp_vLDAU.size()[" << species_pp_vLDAU.size() << "] != vLDAUU.size()[" << vLDAUU.size() << "]" << "   filename=[" << filename << "]" << endl;
      message << "species_pp_vLDAU.size()[" << species_pp_vLDAU.size() << "] != vLDAUU.size()[" << vLDAUU.size() << "]" << "   filename=[" << filename << "]";
      ERROR_flag=TRUE;
    }
    if(species_pp_vLDAU.size()!=vLDAUJ.size()) {
      //[CO20200404 - OBSOLETE]if(!QUIET) cerr << __AFLOW_FUNC__ << " ERROR - species_pp_vLDAU.size()[" << species_pp_vLDAU.size() << "] != vLDAUJ.size()[" << vLDAUJ.size() << "]" << "   filename=[" << filename << "]" << endl;
      message << "species_pp_vLDAU.size()[" << species_pp_vLDAU.size() << "] != vLDAUJ.size()[" << vLDAUJ.size() << "]" << "   filename=[" << filename << "]";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
      ERROR_flag=TRUE;
    }
    for(uint j=0;j<species.size();j++) {
      species_pp_vLDAU.at(j).push_back(LDAUT);
      species_pp_vLDAU.at(j).push_back(vLDAUL.at(j));
      species_pp_vLDAU.at(j).push_back(vLDAUU.at(j));
      species_pp_vLDAU.at(j).push_back(vLDAUJ.at(j));
    }

    if(LDEBUG) { 
      //
      cout << __AFLOW_FUNC__ << " LDA_type=" << LDAUT << endl;
      //
      cout << __AFLOW_FUNC__ << " LDAU_L=";
      for(uint i=0;i<vLDAUL.size();i++) cout << vLDAUL.at(i) << ((i<vLDAUL.size()-1)?",":""); 
      cout << endl;
      //
      cout << __AFLOW_FUNC__ << " LDAU_U=";
      for(uint i=0;i<vLDAUU.size();i++) cout << vLDAUU.at(i) << ((i<vLDAUU.size()-1)?",":""); 
      cout << endl;
      //
      cout << __AFLOW_FUNC__ << " LDAU_J=";
      for(uint i=0;i<vLDAUJ.size();i++) cout << vLDAUJ.at(i) << ((i<vLDAUJ.size()-1)?",":"");
      cout << endl;
    }
    sdata_ldau << aurostd::utype2string(LDAUT)+";";
    for(uint i=0;i<vLDAUL.size();i++) sdata_ldau << vLDAUL.at(i) << ((i<vLDAUL.size()-1)?",":"");
    sdata_ldau << ";";
    for(uint i=0;i<vLDAUU.size();i++) sdata_ldau << vLDAUU.at(i) << ((i<vLDAUU.size()-1)?",":"");
    sdata_ldau << ";";
    for(uint i=0;i<vLDAUJ.size();i++) sdata_ldau << vLDAUJ.at(i) << ((i<vLDAUJ.size()-1)?",":"");
    //sdata_ldau << ";";
    string_LDAU=sdata_ldau.str();
  }else{
    //CO20210713 - push back 0 for no +U
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " LDAU calculation NOT FOUND in OUTCAR" << endl;
    for(uint j=0;j<species.size();j++) {
      species_pp_vLDAU.at(j).push_back(LDAUT);
    }
    if(LDEBUG) { cout << __AFLOW_FUNC__ << " LDA_type=" << LDAUT << endl;}
    //[CO+ME20210713 - keep legacy behavior, only print when non-zero]sdata_ldau << aurostd::utype2string(LDAUT);
    //[CO+ME20210713 - keep legacy behavior, only print when non-zero]string_LDAU=sdata_ldau.str();
  }

  if(LDEBUG) cerr << __AFLOW_FUNC__ << " string_LDAU=" << string_LDAU << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " species_pp_vLDAU.size()=" << species_pp_vLDAU.size() << endl;

  // ----------------------------------------------------------------------
  // DOS related values:
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD DOS related values DATA (" << time_delay(seconds) << ")" << endl;
  EMIN=0.0;EMAX=0.0;SIGMA=0.0;ISMEAR=0;  // for aflowlib_libraries.cpp 
  vline.clear();
  for(uint iline=vcontentRED.size()-1;iline<vcontentRED.size();iline--)  //CO20200404
    //  for(uint iline=0;iline<vcontentRED.size();iline++)
  { //CO20200106 - patching for auto-indenting
    if(aurostd::substring2bool(vcontentRED.at(iline),"EMIN") && aurostd::substring2bool(vcontentRED.at(iline),"energy-range for DOS")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"EMAX") && aurostd::substring2bool(vcontentRED.at(iline),"energy-range for DOS")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"ISMEAR") && aurostd::substring2bool(vcontentRED.at(iline),"broadening in eV")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"SIGMA") && aurostd::substring2bool(vcontentRED.at(iline),"broadening in eV")) vline.push_back(vcontentRED.at(iline));
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vline.size()=" << vline.size() << endl;
  if(vline.size()==0) {
    //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - "<< __AFLOW_FUNC__ << " wrong number of \"DOS related values\" in OUTCAR" << "   filename=[" << filename << "]" << endl;
    message << "Wrong number of \"DOS related values\" in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  }
  for(uint j=0;j<vline.size();j++) {
    //   if(LDEBUG) cerr << __AFLOW_FUNC__ << " vline.at(" << j << ")=" << vline.at(j) << endl;
    aurostd::StringSubst(vline.at(j),"="," ");
    aurostd::StringSubst(vline.at(j),";"," ");
    //    if(LDEBUG) cerr << vline.at(j) << endl;
    aurostd::string2tokens(vline.at(j),tokens," ");
    for(uint k=0;k<tokens.size();k++) {
      if(tokens.at(k)=="EMIN" && k+1<tokens.size()) EMIN=aurostd::string2utype<double>(tokens.at(k+1));
      if(tokens.at(k)=="EMAX" && k+1<tokens.size()) EMAX=aurostd::string2utype<double>(tokens.at(k+1));
      if(tokens.at(k)=="ISMEAR" && k+1<tokens.size()) ISMEAR=aurostd::string2utype<int>(tokens.at(k+1));
      if(tokens.at(k)=="SIGMA" && k+1<tokens.size()) SIGMA=aurostd::string2utype<double>(tokens.at(k+1));
    }
  }
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " EMIN=" << EMIN << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " EMAX=" << EMAX << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " ISMEAR=" << ISMEAR << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " SIGMA=" << SIGMA << endl;}

  // ----------------------------------------------------------------------
  // Electronic relaxation
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD Electronic relaxation DATA (" << time_delay(seconds) << ")" << endl;
  IALGO=0;                      // for aflowlib_libraries.cpp
  LDIAG="";                     // for aflowlib_libraries.cpp
  IMIX=0;INIMIX=0;MIXPRE=0;     // for aflowlib_libraries.cpp
  AMIX=0.0;BMIX=0.0;AMIX_MAG=0.0;BMIX_MAG=0.0;AMIN=0.0;WC=0.0;  // for aflowlib_libraries.cpp
  vline.clear();
  for(uint iline=vcontentRED.size()-1;iline<vcontentRED.size();iline--)  //CO20200404  // DOWN
    // for(uint iline=0;iline<vcontentRED.size();iline++)  // UP
  { //CO20200106 - patching for auto-indenting
    if(aurostd::substring2bool(vcontentRED.at(iline),"IALGO") && aurostd::substring2bool(vcontentRED.at(iline),"algorithm")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"LDIAG") && aurostd::substring2bool(vcontentRED.at(iline),"sub-space diagonalisation")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"IMIX") && aurostd::substring2bool(vcontentRED.at(iline),"mixing-type and parameters")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"AMIX") && aurostd::substring2bool(vcontentRED.at(iline),"BMIX")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"AMIX_MAG") && aurostd::substring2bool(vcontentRED.at(iline),"BMIX_MAG")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"AMIN")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"WC") && aurostd::substring2bool(vcontentRED.at(iline),"INIMIX") && aurostd::substring2bool(vcontentRED.at(iline),"MIXPRE")) vline.push_back(vcontentRED.at(iline));
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vline.size()=" << vline.size() << endl;
  if(vline.size()==0) {
    //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - "<< __AFLOW_FUNC__ << " wrong number of \"Electronic relaxation\" in OUTCAR" << "   filename=[" << filename << "]" << endl;
    message << "Wrong number of \"Electronic relaxation\" in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  }
  for(uint j=0;j<vline.size();j++) {
    //  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vline.at(" << j << ")=" << vline.at(j) << endl;
    aurostd::StringSubst(vline.at(j),"="," ");
    aurostd::StringSubst(vline.at(j),";"," ");
    //    if(LDEBUG) cerr << vline.at(j) << endl;
    aurostd::string2tokens(vline.at(j),tokens," ");
    for(uint k=0;k<tokens.size();k++) {
      if(tokens.at(k)=="IALGO" && k+1<tokens.size()) IALGO=aurostd::string2utype<int>(tokens.at(k+1));
      if(tokens.at(k)=="LDIAG" && k+1<tokens.size()) LDIAG=(tokens.at(k+1));
      if(tokens.at(k)=="IMIX" && k+1<tokens.size()) IMIX=aurostd::string2utype<int>(tokens.at(k+1));
      if(tokens.at(k)=="AMIX" && k+1<tokens.size()) AMIX=aurostd::string2utype<double>(tokens.at(k+1));
      if(tokens.at(k)=="BMIX" && k+1<tokens.size()) BMIX=aurostd::string2utype<double>(tokens.at(k+1));
      if(tokens.at(k)=="AMIX_MAG" && k+1<tokens.size()) AMIX_MAG=aurostd::string2utype<double>(tokens.at(k+1));
      if(tokens.at(k)=="BMIX_MAG" && k+1<tokens.size()) BMIX_MAG=aurostd::string2utype<double>(tokens.at(k+1));
      if(tokens.at(k)=="AMIN" && k+1<tokens.size()) AMIN=aurostd::string2utype<double>(tokens.at(k+1));
      if(tokens.at(k)=="WC" && k+1<tokens.size()) WC=aurostd::string2utype<double>(tokens.at(k+1));
      if(tokens.at(k)=="INIMIX" && k+1<tokens.size()) INIMIX=aurostd::string2utype<int>(tokens.at(k+1));
      if(tokens.at(k)=="MIXPRE" && k+1<tokens.size()) MIXPRE=aurostd::string2utype<int>(tokens.at(k+1));
    }
  }
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " IALGO=" << IALGO << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " LDIAG=" << LDIAG << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " IMIX=" << IMIX << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " AMIX=" << AMIX << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " BMIX=" << BMIX << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " AMIX_MAG=" << AMIX_MAG << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " BMIX_MAG=" << BMIX_MAG << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " AMIN=" << AMIN << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " WC=" << WC << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " INIMIX=" << INIMIX << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " MIXPRE=" << MIXPRE << endl;}

  // ----------------------------------------------------------------------
  // Intra band minimization
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD Intra band minimization DATA (" << time_delay(seconds) << ")" << endl;
  WEIMIN=0.0;EBREAK=0.0;DEPER=0.0;TIME=0.0;  // for aflowlib_libraries.cpp
  vline.clear();
  for(uint iline=vcontentRED.size()-1;iline<vcontentRED.size();iline--)  //CO20200404
    //  for(uint iline=0;iline<vcontentRED.size();iline++)
  { //CO20200106 - patching for auto-indenting
    if(aurostd::substring2bool(vcontentRED.at(iline),"WEIMIN") && aurostd::substring2bool(vcontentRED.at(iline),"energy-eigenvalue tresh-hold")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"EBREAK") && aurostd::substring2bool(vcontentRED.at(iline),"absolut break condition")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"DEPER") && aurostd::substring2bool(vcontentRED.at(iline),"relativ break condition")) vline.push_back(vcontentRED.at(iline));
    if(aurostd::substring2bool(vcontentRED.at(iline),"TIME") && aurostd::substring2bool(vcontentRED.at(iline),"timestep for ELM")) vline.push_back(vcontentRED.at(iline));
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vline.size()=" << vline.size() << endl;
  if(vline.size()==0) {
    //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - "<< __AFLOW_FUNC__ << " wrong number of \"Intra band minimization\" in OUTCAR" << "   filename=[" << filename << "]" << endl;
    message << "Wrong number of \"Intra band minimization\" in OUTCAR" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  }
  for(uint j=0;j<vline.size();j++) {
    //  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vline.at(" << j << ")=" << vline.at(j) << endl;
    aurostd::StringSubst(vline.at(j),"="," ");
    aurostd::StringSubst(vline.at(j),";"," ");
    //    if(LDEBUG) cerr << vline.at(j) << endl;
    aurostd::string2tokens(vline.at(j),tokens," ");
    for(uint k=0;k<tokens.size();k++) {
      if(tokens.at(k)=="WEIMIN" && k+1<tokens.size()) WEIMIN=aurostd::string2utype<double>(tokens.at(k+1));
      if(tokens.at(k)=="EBREAK" && k+1<tokens.size()) EBREAK=aurostd::string2utype<double>(tokens.at(k+1));
      if(tokens.at(k)=="DEPER" && k+1<tokens.size()) DEPER=aurostd::string2utype<double>(tokens.at(k+1));
      if(tokens.at(k)=="TIME" && k+1<tokens.size()) TIME=aurostd::string2utype<double>(tokens.at(k+1));
    }
  }
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " WEIMIN=" << WEIMIN << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " EBREAK=" << EBREAK << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " DEPER=" << DEPER << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " TIME=" << TIME; cerr << endl;}

  // ----------------------------------------------------------------------
  // LOAD NWEIGHTS VKPOINT
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD NWEIGHTS KPOINTLIST DATA (" << time_delay(seconds) << ")" << endl;
  uint nkpoints_line=0;
  nweights=0;
  nkpoints_irreducible=0;
  vkpoint_reciprocal.clear();                    // QM KPOINTLIST calculation
  vkpoint_cartesian.clear();                     // QM KPOINTLIST calculation
  vweights.clear();                              // QM KPOINTLIST calculation
  uint minweight=1e6;
  vline.clear();
  for(uint iline=0;iline<vcontent.size();iline++)
    //  for(uint iline=0;iline<vcontent.size();iline++)
  { //CO20200106 - patching for auto-indenting
    if(aurostd::substring2bool(vcontent.at(iline),"Found") && aurostd::substring2bool(vcontent.at(iline),"irreducible k-points")) {vline.push_back(vcontent.at(iline));nkpoints_line=iline+4;}
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vline.size()=" << vline.size() << endl;
  //[CO20200404 - breaks for all OUTCAR.bands]if(vline.size()==0) {
  //[CO20200404 - breaks for all OUTCAR.bands]  //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - "<< __AFLOW_FUNC__ << " wrong number of \" LOAD NWEIGHTS VKPOINT\" in OUTCAR" << "   filename=[" << filename << "]" << endl;
  //[CO20200404 - breaks for all OUTCAR.bands]  message << "Wrong number of \" LOAD NWEIGHTS VKPOINT\" in OUTCAR" << "   filename=[" << filename << "]";
  //[CO20200404 - breaks for all OUTCAR.bands]  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
  //[CO20200404 - breaks for all OUTCAR.bands]  ERROR_flag=TRUE;
  for(uint j=0;j<vline.size();j++) {
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " vline.at(" << j << ")=" << vline.at(j) << endl;
    aurostd::StringSubst(vline.at(j),"="," ");
    aurostd::StringSubst(vline.at(j),";"," ");
    //    if(LDEBUG) cerr << vline.at(j) << endl;
    aurostd::string2tokens(vline.at(j),tokens," ");
    for(uint k=0;k<tokens.size();k++) {
      if(tokens.at(k)=="Found" && k+1<tokens.size()) nkpoints_irreducible=aurostd::string2utype<uint>(tokens.at(k+1));
    }
  }
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " nkpoints_irreducible=" << nkpoints_irreducible << endl;}

  if(nkpoints_irreducible) {
    for(uint iline=nkpoints_line;(iline<nkpoints_line+nkpoints_irreducible && iline<vcontent.size());iline++) {
      aurostd::string2tokens(vcontent.at(iline),tokens," ");
      //  cerr << vcontent.at(iline) << endl;
      if(tokens.size()==4) {
        xvector<double> kpoint(3);
        kpoint[1]=aurostd::string2utype<double>(tokens.at(0));
        kpoint[2]=aurostd::string2utype<double>(tokens.at(1));
        kpoint[3]=aurostd::string2utype<double>(tokens.at(2));
        vkpoint_reciprocal.push_back(kpoint); // cerr.precision(20);
        vweights.push_back(aurostd::string2utype<uint>(tokens.at(3)));
        minweight=aurostd::min(minweight,aurostd::string2utype<uint>(tokens.at(3)));
        nweights+=aurostd::string2utype<uint>(tokens.at(3));
      } else {
        //	if(!QUIET) cerr << __AFLOW_FUNC__ << " error in NWEIGHTS/VKPOINT calculation tokens.size()=" << tokens.size() << endl;
      }
    }
    nkpoints_line+=3;
    for(uint iline=nkpoints_line;(iline<nkpoints_line+nkpoints_irreducible && iline<vcontent.size());iline++) {
      aurostd::string2tokens(vcontent.at(iline),tokens," ");
      // cerr << vcontent.at(iline) << endl;
      if(tokens.size()==4) {
        xvector<double> kpoint(3);
        kpoint[1]=aurostd::string2utype<double>(tokens.at(0));
        kpoint[2]=aurostd::string2utype<double>(tokens.at(1));
        kpoint[3]=aurostd::string2utype<double>(tokens.at(2));
        vkpoint_cartesian.push_back(kpoint); // cerr.precision(20);
      } else {
        //	if(!QUIET) cerr << __AFLOW_FUNC__ << " error in NWEIGHTS/VKPOINT calculation tokens.size()=" << tokens.size() << endl;
      }
    }
  }
  //  nweights=nweights/minweight; 
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vkpoint_reciprocal.size()=" << vkpoint_reciprocal.size() << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vkpoint_cartesian.size()=" << vkpoint_cartesian.size() << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vweights.size()=" << vweights.size() << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " nweights=" << nweights << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " nkpoints_irreducible=" << nkpoints_irreducible << endl;

  // ----------------------------------------------------------------------
  // LOAD CALCULATION STUFF
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD CALCULATION STUFF DATA (" << time_delay(seconds) << ")" << endl;
  // CALCULATION_CORES
  calculation_cores=1;
  for(uint iline=0;iline<vcontent.size();iline++) 
    if(aurostd::substring2bool(vcontent.at(iline),"running"))
      if(aurostd::substring2bool(vcontent.at(iline),"nodes")) {
        line=vcontent.at(iline);
        break;
      }
  if(!line.empty()) {
    aurostd::string2tokens(line,tokens);    // cerr << tokens.at(2) << endl;
    if(tokens.size()>2) calculation_cores=aurostd::string2utype<uint>(tokens.at(2));
  }
  if(calculation_cores<1) calculation_cores=1; 
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " calculation_cores=" << calculation_cores << endl;
  // CALCULATION_TIME
  calculation_time=0.0;
  for(int iline=(int)vcontent.size()-1;iline>=0;iline--){  // NEW FROM THE BACK
    if(aurostd::substring2bool(vcontent.at(iline),"time")){
      if(aurostd::substring2bool(vcontent.at(iline),"Total")) {
        line=vcontent.at(iline);
        break;
      } 
    }
  }
  if(line.empty()) {
    //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - "<< __AFLOW_FUNC__ << " in OUTCAR (no calculation_time)" << "   filename=[" << filename << "]" << endl;
    message << "In OUTCAR (no calculation_time)" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  }
  aurostd::string2tokens(line,tokens);
  if(tokens.size()>1) calculation_time=aurostd::string2utype<double>(tokens.at(tokens.size()-1));
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " calculation_time=" << calculation_time << endl;
  // CALCULATION_MEMORY 
  calculation_memory=0.0;
  for(int iline=(int)vcontent.size()-1;iline>=0;iline--)  // NEW FROM THE BACK
    if(aurostd::substring2bool(vcontent.at(iline),"MB"))
      if(aurostd::substring2bool(vcontent.at(iline),"storing wavefunctions")) {
        line=vcontent.at(iline);
        break;
      } 
  if(line.empty()) {
    //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - "<< __AFLOW_FUNC__ << " in OUTCAR (no calculation_memory)" << "   filename=[" << filename << "]" << endl;
    message << "In OUTCAR (no calculation_memory)" << "   filename=[" << filename << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  }
  aurostd::string2tokens(line,tokens); //   cerr << tokens.at(3) << endl;
  if(tokens.size()>3) calculation_memory=aurostd::string2utype<double>(tokens.at(3));
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " calculation_memory=" << calculation_memory << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ---------------------------------" << endl;

  // ----------------------------------------------------------------------
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " END (" << time_delay(seconds) << ")" << endl;

  // DONE NOW RETURN
  //[CO20200404 - OBSOLETE]if(ERROR_flag && !QUIET) cerr << "WARNING - "<< __AFLOW_FUNC__ << " ERROR_flag set in xOUTCAR" << endl;
  if(ERROR_flag){
    message << "ERROR_flag set in xOUTCAR";
    if(force_exit){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);}
    else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_,QUIET);}
    return FALSE;
  }
  return TRUE;
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
bool xOUTCAR::GetXStructure() {
  //CO20211107 - this is good for OUTCAR.static or .bands which only has one structure inside
  //use GetIonicStepsData() otherwise
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  xstr.clear(); //DX20191220 - uppercase to lowercase clear

  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " Trying to build the xstructure from the OUTCAR" << endl;}

  //get lattice
  bool found_lattice,found_types,found_positions;
  found_lattice=found_types=found_positions=false;
  string line = "", token = "";
  vector<string> tokens;
  xmatrix<double> lattice(3,3),klattice(3,3);
  deque<int> num_each_type;
  double num = 0.0, natoms = 0.0;
  deque<_atom> atoms;
  xvector<double> cpos(3),fpos(3);
  uint vcontent_size=vcontent.size();
  for(uint iline=vcontent_size-1;iline<vcontent_size;iline--){  // NEW FROM THE BACK
    //get lattice
    if(!found_lattice && (iline<vcontent_size-3) && aurostd::substring2bool(vcontent[iline],"direct")) {
      aurostd::string2tokens(vcontent[iline],tokens);
      if((tokens.size()>5) && (tokens[0] == "direct") && (tokens[1] == "lattice") && (tokens[2] == "vectors") &&
          (tokens[3] == "reciprocal") && (tokens[4] == "lattice") && (tokens[5] == "vectors")){
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " lattice found!" << endl;}
        found_lattice=true;
        //
        tokens=GetCorrectPositions(vcontent[iline+1],6);
        if(!tokens.size()){
          if(LDEBUG) {cerr << __AFLOW_FUNC__ << " line with lattice vector is ill-written, see: " << vcontent[iline+1] << endl;}
          return false;
        }
        lattice(1,1)=aurostd::string2utype<double>(tokens[0]);
        lattice(1,2)=aurostd::string2utype<double>(tokens[1]);
        lattice(1,3)=aurostd::string2utype<double>(tokens[2]);
        klattice(1,1)=aurostd::string2utype<double>(tokens[3]);
        klattice(1,2)=aurostd::string2utype<double>(tokens[4]);
        klattice(1,3)=aurostd::string2utype<double>(tokens[5]);
        //
        tokens=GetCorrectPositions(vcontent[iline+2],6);
        if(!tokens.size()){
          if(LDEBUG) {cerr << __AFLOW_FUNC__ << " line with lattice vector is ill-written, see: " << vcontent[iline+2] << endl;}
          return false;
        }
        lattice(2,1)=aurostd::string2utype<double>(tokens[0]);
        lattice(2,2)=aurostd::string2utype<double>(tokens[1]);
        lattice(2,3)=aurostd::string2utype<double>(tokens[2]);
        klattice(2,1)=aurostd::string2utype<double>(tokens[3]);
        klattice(2,2)=aurostd::string2utype<double>(tokens[4]);
        klattice(2,3)=aurostd::string2utype<double>(tokens[5]);
        //
        tokens=GetCorrectPositions(vcontent[iline+3],6);
        if(!tokens.size()){
          if(LDEBUG) {cerr << __AFLOW_FUNC__ << " line with lattice vector is ill-written, see: " << vcontent[iline+3] << endl;}
          return false;
        }
        lattice(3,1)=aurostd::string2utype<double>(tokens[0]);
        lattice(3,2)=aurostd::string2utype<double>(tokens[1]);
        lattice(3,3)=aurostd::string2utype<double>(tokens[2]);
        klattice(3,1)=aurostd::string2utype<double>(tokens[3]);
        klattice(3,2)=aurostd::string2utype<double>(tokens[4]);
        klattice(3,3)=aurostd::string2utype<double>(tokens[5]);
      }
    }
    //get types
    if(!found_types && aurostd::substring2bool(vcontent[iline],"type")) {
      aurostd::string2tokens(vcontent[iline],tokens);
      if((tokens.size()>2) && (tokens[0]=="ions") && (tokens[1]=="per") && (tokens[2]=="type")){
        aurostd::string2tokens(vcontent[iline],tokens,"=");
        if(tokens.size()!=2){continue;}
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " types found!" << endl;}
        found_types=true;
        token=tokens[1];
        aurostd::string2tokens(token,tokens);
        natoms=0;
        for(uint i=0;i<tokens.size();i++){
          num=aurostd::string2utype<double>(tokens[i]);
          natoms+=num;
          num_each_type.push_back((int)num);
        }
      }
    }
  }
  if(!found_lattice){
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " lattice not found" << endl;}
    return false;
  }
  if(!found_types){
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " types not found" << endl;}
    return false;
  }

  xmatrix<double> f2c=trasp(lattice), c2f=inverse(f2c);

  //need types before positions
  for(int iline=(int)vcontent_size-1;iline>=0;iline--){  // NEW FROM THE BACK
    //get positions
    //from bottom up, cartesian coordinates found first
    //this is fortunate as they carry more precision (fractional needs multiplication)
    if(!found_positions && aurostd::substring2bool(vcontent[iline],"coordinates")){
      aurostd::string2tokens(vcontent[iline],tokens);
      if(!found_positions && (iline<(int)vcontent_size-natoms) && (tokens.size()>5) && (tokens[0] == "position") && (tokens[1] == "of") && (tokens[2] == "ions") &&
          (tokens[3] == "in") && (tokens[4] == "cartesian") && (tokens[5] == "coordinates")){
        found_positions=true;
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " positions (cartesian) found!" << endl;}
        atoms.clear();
        for(uint itype=0;itype<num_each_type.size();itype++){
          for(uint iatom=0;iatom<(uint)num_each_type[itype];iatom++){
            atoms.push_back(_atom());
            atoms.back().type=itype;
            tokens=GetCorrectPositions(vcontent[iline+atoms.size()],3);
            if(!tokens.size()){
              if(LDEBUG) {cerr << __AFLOW_FUNC__ << " line with atom positions is ill-written, see: " << vcontent[iline+atoms.size()] << endl;}
              return false;
            }
            cpos(1)=aurostd::string2utype<double>(tokens[0]);
            cpos(2)=aurostd::string2utype<double>(tokens[1]);
            cpos(3)=aurostd::string2utype<double>(tokens[2]);
            atoms.back().cpos=cpos;
            atoms.back().fpos=c2f*cpos;
          }
        }
      }
      if(!found_positions && (iline<(int)vcontent_size-natoms) && (tokens.size()>5) && (tokens[0] == "position") && (tokens[1] == "of") && (tokens[2] == "ions") &&
          (tokens[3] == "in") && (tokens[4] == "fractional") && (tokens[5] == "coordinates")){
        found_positions=true;
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " positions (fractional) found!" << endl;}
        atoms.clear();
        for(uint itype=0;itype<num_each_type.size();itype++){
          for(uint iatom=0;iatom<(uint)num_each_type[itype];iatom++){
            atoms.push_back(_atom());
            atoms.back().type=itype;
            tokens=GetCorrectPositions(vcontent[iline+atoms.size()],3);
            if(!tokens.size()){
              if(LDEBUG) {cerr << __AFLOW_FUNC__ << " line with atom positions is ill-written, see: " << vcontent[iline+atoms.size()] << endl;}
              return false;
            }
            fpos(1)=aurostd::string2utype<double>(tokens[0]);
            fpos(2)=aurostd::string2utype<double>(tokens[1]);
            fpos(3)=aurostd::string2utype<double>(tokens[2]);
            atoms.back().fpos=fpos;
            atoms.back().cpos=f2c*fpos;
          }
        }
      }
    }
  }

  if(!found_positions){
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " atom positions not found" << endl;}
    return false;
  }

  //default!
  for(uint i=0;i<atoms.size();i++){
    atoms[i].name="";
    atoms[i].name_is_given=false;
    atoms[i].CleanName();
    atoms[i].CleanSpin();
    atoms[i].ijk(1)=0;atoms[i].ijk(2)=0;atoms[i].ijk(3)=0; // inside the zero cell...
    atoms[i].corigin(1)=0.0;atoms[i].corigin(2)=0.0;atoms[i].corigin(3)=0.0; // inside the zero cell
    atoms[i].coord(1)=0.0;atoms[i].coord(2)=0.0;atoms[i].coord(3)=0.0; // inside the zero cell
    atoms[i].spin=0.0;
    atoms[i].order_parameter_value=0;
    atoms[i].order_parameter_atom=FALSE;
    atoms[i].partial_occupation_value=1.0;
    atoms[i].partial_occupation_flag=FALSE;
  }

  //occupy xstructure
  xstr.title=SYSTEM;
  xstr.lattice=lattice;
  xstr.klattice=2.0*pi*klattice;
  //xstr.num_each_type=num_each_type; //done in AddAtom()
  //xstr.comp_each_type.clear();      //done in AddAtom()
  //for(uint i=0;i<xstr.num_each_type.size();i++){xstr.comp_each_type.push_back((double)xstr.num_each_type[i]);}  //default!, done in AddAtom()
  //xstr.natoms=natoms; //not an attribute
  xstr.f2c=f2c;
  xstr.c2f=c2f;
  xstr.FixLattices();
  //xstr.atoms=atoms;
  for(uint i=0;i<atoms.size();i++){
    xstr.AddAtom(atoms[i],false); //CO20230319 - add by type
    xstr.partial_occupation_sublattice.push_back(_pocc_no_sublattice_); //default!
  }
  xstr.MakeBasis();

  if(LDEBUG) {
    cerr << __AFLOW_FUNC__ << " lattice" << endl;
    cerr << lattice << endl;
    cerr << __AFLOW_FUNC__ << " klattice" << endl;
    cerr << 2.0*pi*klattice << endl;
    cerr << __AFLOW_FUNC__ << " reciprocal of lattice (check)" << endl;
    cerr << ReciprocalLattice(lattice) << endl;
    cerr << __AFLOW_FUNC__ << " natoms = " << natoms << endl;
    for(uint i=0;i<num_each_type.size();i++){
      cerr << __AFLOW_FUNC__ << " num_each_type[" << i << "] = " << num_each_type[i] << endl;
    }
    for(uint i=0;i<atoms.size();i++){
      cerr << __AFLOW_FUNC__ << " atom[" << i << "] type=" << atoms[i].type << " cpos=" << atoms[i].cpos << " fpos=" << atoms[i].fpos << endl;
    }
    cerr << __AFLOW_FUNC__ << " full xstructure" << endl;
    cerr << xstr << endl;
  }

  return true;
}

int xOUTCAR::isKPointLine(uint iline){
  xvector<double> kpoint;
  return isKPointLine(iline,kpoint);
}

int xOUTCAR::isKPointLine(uint iline,xvector<double>& _kpoint){
  if(iline>vcontent.size()-2){return 0;}  //can't be last line, we check iline+1
  if(!aurostd::substring2bool(vcontent[iline],"k-point")){return 0;}
  vector<string> tokens;
  aurostd::string2tokens(vcontent[iline],tokens);
  if(tokens.size()!=6){return 0;}
  if(tokens[0]!="k-point"){return 0;}
  //snag kpoint index
  int kpt_num=(tokens[1]=="***" ? -1 : aurostd::string2utype<int>(tokens[1]));  //issue with printing kpt_num > 999
  //snag kpoint
  xvector<double> kpoint(3);
  for(uint i=3;i<tokens.size();i++){kpoint[i-2]=aurostd::string2utype<double>(tokens[i]);} 
  _kpoint=kpoint;
  //check next line
  aurostd::string2tokens(vcontent[iline+1],tokens);
  if(tokens.size()!=5){return 0;}
  if(tokens[0]!="band"){return 0;}
  if(tokens[1]!="No."){return 0;}
  if(tokens[2]!="band"){return 0;}
  if(tokens[3]!="energies"){return 0;}
  if(tokens[4]!="occupation"){return 0;}
  return kpt_num;
}

bool xOUTCAR::GetStartingKPointLines(vector<uint>& ilines) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  ilines.clear();
  if(!(ISPIN==1 || ISPIN==2)){
    if(!GetProperties(content) || !(ISPIN==1 || ISPIN==2)){
      if(LDEBUG) {cerr << XPID << "xOUTCAR::GetStartingKPointLines: GetProperties failed." << endl;}
      return false;
    }
  }
  vector<string> tokens;
  xvector<double> kpoint;
  for(uint iline=0;iline<vcontent.size();iline++){
    if(aurostd::substring2bool(vcontent[iline],"k-point")){
      if(1==isKPointLine(iline,kpoint)){ilines.push_back(iline);}
    }
  }
  if(ilines.size()!=(uint)ISPIN){
    if(LDEBUG) {cerr << XPID << "xOUTCAR::GetStartingKPointLines: ISPIN does not match starting k-point set counts" << endl;}
    return false;
  }
  return true;
}

bool xOUTCAR::GetNextKPointLine(uint& iline){
  iline++;  //march forward
  xvector<double> kpoint;
  while(iline<vcontent.size() && !isKPointLine(iline,kpoint)){iline++;} //will work for k-point 1000 (***)
  if(iline>vcontent.size()-1){return false;}
  return true;
}

bool xOUTCAR::ProcessKPoint(uint iline,double EFERMI,vector<double>& b_energies,vector<double>& b_occs){
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  b_energies.clear(); b_occs.clear();

  if(!isKPointLine(iline)){
    if(LDEBUG) {cerr << XPID << "xOUTCAR::ProcessKPoint: this is NOT a k-point line" << endl;}
    return false;
  }
  iline+=2; //march forward to first band

  vector<string> tokens;
  uint iband_found;
  double b_energy, b_occ;
  for(uint iband=1;iband<(uint)NBANDS+1;iband++){
    if(iline>vcontent.size()-1){
      if(LDEBUG) {cerr << XPID << "xOUTCAR::ProcessKPoint: reached the end of the file" << endl;}
      return false;
    }
    aurostd::string2tokens(vcontent[iline],tokens);
    if(tokens.size()!=3){
      if(LDEBUG) {cerr << XPID << "xOUTCAR::ProcessKPoint: odd count of tokens(3) for iband=" << iband << endl;}
      return false;
    }
    iband_found=aurostd::string2utype<int>(tokens[0]);
    b_energy=aurostd::string2utype<double>(tokens[1]);
    b_occ=aurostd::string2utype<double>(tokens[2]);
    if(iband_found!=iband){
      if(LDEBUG) {cerr << XPID << "xOUTCAR::ProcessKPoint: missing iband=" << iband << endl;}
      return false;
    }
    b_energies.push_back(b_energy-EFERMI);  //adjust for E-fermi (STATIC) here!!!!!!!!!!!!!!! (not the OUTCAR Efermi, but one from STATIC hopefully)
    b_occs.push_back(b_occ);
    iline++;
  }
  return true;
}

//if energies are equal, sort by occs next in DESCENDING order
bool xOUTCAR::bandEnergyOccCompare::operator()(const bandEnergyOcc& a,const bandEnergyOcc b) const {
  return (a.energy<b.energy || (aurostd::identical(a.energy,b.energy,energy_tol) && a.occ>b.occ));  //sort by energy ASCENDING, occs DESCENDING
}

bool xOUTCAR::orderBands(vector<double>& b_energies,vector<double>& b_occs,double energy_tol){
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  if(b_energies.size() != b_occs.size()){
    if(LDEBUG) {cerr << XPID << "xOUTCAR::orderBands: size of energies != size of occupations" << endl;} 
    return false;
  }
  if(b_energies.size()==0){
    if(LDEBUG) {cerr << XPID << "xOUTCAR::orderBands: no input energies or occupations found" << endl;} 
    return false;
  }
  uint len=b_energies.size();
  vector<bandEnergyOcc> beo;
  for(uint i=0;i<b_energies.size();i++){
    beo.push_back(bandEnergyOcc());
    beo.back().energy=b_energies[i];
    beo.back().occ=b_occs[i];
  }
  std::sort(beo.begin(),beo.end(),bandEnergyOccCompare(energy_tol));
  b_energies.clear();
  b_occs.clear();
  for(uint i=0;i<beo.size();i++){
    b_energies.push_back(beo[i].energy);
    b_occs.push_back(beo[i].occ);
  }
  if(b_energies.size()!=len || b_occs.size()!=len){
    if(LDEBUG) {cerr << XPID << "xOUTCAR::orderBands: size mismatch from input" << endl;} 
    return false;
  }
  return true;
}

bool xOUTCAR::GetBandEdge(vector<double>& b_energies,vector<double>& b_occs,double EFERMI,uint& iedge,double efermi_tol,double energy_tol,double occ_tol){
  //MOST IMPORTANT FUNCTION FOR FINDING BANDGAP
  //edge is defined here as valence band maximum
  //since energies run from min to max, conduction band min is iedge+1
  //we go backwards, max to E-fermi (above is really below), to find the first band below E-fermi
  //remember, E-fermi is calculated such that an integral of DOS up to E-fermi == # of electrons
  //but it's only reliable in STATIC and not BANDS, where STATIC resolves energies self-consistently (ICHARG=11)
  //BANDS calc only samples path at specific places in BZ, which may not be representative of the full charge density
  //this discrepancy is important, and quantifies the error in the energies of the BANDS calc
  //E-fermi BANDS may differ from E-fermi STATIC by quite a bit (I've seen as much as 0.5 eV)
  //we already shift energies by the INPUT E-fermi, which should be the one from the STATIC calculation
  //first, we want to see if a sharp edge exists
  //this is where a clear fully occupied (occ=2 for ISPIN==1, occ=1 for ISPIN==2) / unoccupied (occ=0) transition exists
  //return immediately if found
  //otherwise, look for a soft edge
  //soft edges occur in band(s) near E-fermi
  //any occupancy (>0) above E-fermi is GARBAGE (an artifact of smearing/sigma for convergence), so ignore
  //but how do we differentiate garbage from non-garbage?
  //algorithmically, we search for the biggest change in occupancy, and define this to be where the transition occurs
  //WARNING: if NUPDOWN is set manually in INCAR, E-fermi must be calculated manually as well (vasp E-fermi only applies to spin-down)
  //http://cms.mpi.univie.ac.at/vasp-forum/viewtopic.php?f=4&t=7442
  //INPUT: energies + occs, E-fermi from STATIC, efermi_tol (largest shift in energy levels between BANDS and STATIC),
  //energy_tol is resolution of energy level values, and occ_tol is resolution of occupation values
  //OUTPUT: iedge - index of VBT

  bool LDEBUG=(FALSE || XHOST.DEBUG);

  //////////////////////////////////////////////////////////////////////////////
  //tests of stupidity ROUND 1 - START

  if(b_energies.size() != b_occs.size()){
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " size of energies != size of occupations" << endl;} 
    return false;
  }
  if(b_energies.size()==0){
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " no input energies or occupations found" << endl;} 
    return false;
  }

  if(!(ISPIN==1 || ISPIN==2) || NKPTS==0 || NBANDS==0){
    if(!GetProperties(content) || !(ISPIN==1 || ISPIN==2) || NKPTS==0 || NBANDS==0){
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " GetProperties failed." << endl;}
      return false;
    }
  }

  //error in E-fermi, unless otherwise provided, is the difference between the input E-fermi (STATIC) and the one found in this OUTCAR (BANDS)
  //can be significant!
  //think of this as largest possible shift of energies in BANDS relative to real energy levels
  if(efermi_tol==AUROSTD_NAN){
    efermi_tol=Efermi-EFERMI;                             //generlly, E-fermi BANDS > E-fermi STATIC
    if(std::signbit(efermi_tol)){efermi_tol=energy_tol;}  //just in case
  }
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " tol(E-fermi)=" << efermi_tol << endl;}

  double full_occ = (ISPIN==1 ? 2.0 : 1.0);

  //sort band energies + occs
  //band index is completely phony, low energy bands can appear above higher energy bands:  see FCC/Cu5Lu1_ICSD_103045 k-point 172 (spin=1)
  //this does NOT seem to be an IO error, just dummy VASP indices
  //two bands with the same energy should also be sorted by occs:  see ORC/Al1Pd1Sm1_ICSD_609058 k-point 21 (spin=1)
  //looks like VASP band indices are completely meaningless!
  if(!orderBands(b_energies,b_occs,energy_tol)){
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " cannot sort band energies and occupations, this is VERY odd!" << endl;}
    return false;
  }

  //now that everything is sorted...
  //occupations should monotonically INCREASE as energy decreases
  bool found_full_occ=false;
  for(int iband=(int)b_energies.size()-2;iband>=0;iband--){ //-2 because we compare iband and iband+1
    if(found_full_occ){
      if(!(abs(full_occ-b_occs[iband])<occ_tol)){
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " drop from full_occ to non-full_occ (iband=" << iband+1 << "), this is VERY odd!" << endl;}
        return false;
      }
    }
    if(abs(full_occ-b_occs[iband])<occ_tol){found_full_occ=true;}
    if(std::signbit(b_occs[iband]-b_occs[iband+1]) && aurostd::isdifferent(b_occs[iband],b_occs[iband+1],occ_tol)){
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " occupation decreased at lower energy (iband=" << iband+1 << "), this is VERY odd!" << endl;}
      return false;
    }
  }

  //tests of stupidity ROUND 1 - STOP
  //////////////////////////////////////////////////////////////////////////////

  //sharp edge - A REAL TRANSITION
  //a sharp edge is a sharp edge is a sharp edge
  //very strictly find full_occ vs. no_occ transition
  //don't worry about whether the edge is above or below E-fermi HERE since E-fermi(STATIC) != E-fermi(BANDS)
  //VASP wouldn't occupy above what it believes to be E-fermi, by definition
  for(int iband=(int)b_energies.size()-2;iband>=0;iband--){ //-2 because we compare iband and iband+1
    //if(b_energies[iband]>0.020){continue;}  //adjusted for E-fermi (STATIC) already!!!!!
    if(abs(full_occ-b_occs[iband])<occ_tol && abs(b_occs[iband+1])<occ_tol){
      iedge=iband;
      if(LDEBUG) {
        cerr << __AFLOW_FUNC__ << " sharp edge found between energies " << b_energies[iedge] << ",";
        cerr << b_energies[iedge+1] << " (ibands=" << iedge+1 << "," << iedge+2 << ")" << endl;
      }
      return true;
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  //tests of stupidity ROUND 2 - START
  //we do it in this order because we want sharp edges to be found no matter what
  //see for example:  HEX/H2_ICSD_68271

  //band energy edges should not be near E-fermi
  if(b_energies[0]>-efermi_tol){
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " first band energy is near/above E-fermi!" << endl;}
    return false;
  }
  if(b_energies[b_energies.size()-1]<efermi_tol){
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " last band energy is near/below E-fermi!" << endl;}
    return false;
  }

  //first band should be fully occupied
  if(abs(full_occ-b_occs[0])>=occ_tol){
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " first band is not fully occupied!" << endl;}
    return false;
  }
  //last band should be fully un-occupied
  if(abs(b_occs[b_occs.size()-1])>=occ_tol){
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " last band is not fully un-occupied!" << endl;}
    //return false; //don't exit here, we might still have a legitimate soft edge, see FCC/C6Mn20W3_ICSD_618279
  }

  //tests of stupidity ROUND 2 - STOP
  //////////////////////////////////////////////////////////////////////////////

  //soft edge - a fuzzy transition due to numerical/convergence issues
  //here we care whether we are near E-fermi
  //two issues that must be overcome here - (full_occ vs. non-full_occ) and (E-fermi STATIC vs. E-fermi BANDS)
  //1. full_occ vs. non-full_occ: find the largest delta in occupation, set that to be the transition
  //2. E-fermi STATIC vs. E-fermi BANDS: look no higher than the difference of the two (tolerance of E-fermi)
  //so don't look above E-fermi+tol (E-fermi BANDS)
  //VASP wouldn't occupy above what it believes to be E-fermi, by definition
  double occ_delta,occ_delta_max=0;
  found_full_occ=false;
  for(int iband=(int)b_energies.size()-2;iband>=0;iband--){ //-2 because we compare iband and iband+1
    if(b_energies[iband]>efermi_tol){continue;}  //adjusted for E-fermi (STATIC) already!!!!!
    if(abs(full_occ-b_occs[iband])<occ_tol){  //if we find full_occ twice in a row, stop
      if(found_full_occ){break;}
      found_full_occ=true;
    }
    occ_delta=b_occs[iband]-b_occs[iband+1];
    if(occ_delta>occ_delta_max){
      iedge=iband;
      occ_delta_max=occ_delta;
    }
  }

  if(abs(occ_delta_max)<occ_tol){
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " no edge found!" << endl;}
    return false;
  }

  if(LDEBUG) {
    cerr << __AFLOW_FUNC__ << " soft edge found between energies " << b_energies[iedge] << ",";
    cerr << b_energies[iedge+1] << " (ibands=" << iedge+1 << "," << iedge+2 << ")" << endl;
  }
  return true;
}

bool xOUTCAR::identicalKPoints(vector<xvector<double> >& vkpoints,uint kpt1,uint kpt2,double tol){
  if(kpt1==kpt2){return true;}
  return identicalKPoints(vkpoints[kpt1],vkpoints[kpt2],tol);
}

bool xOUTCAR::identicalKPoints(xvector<double>& kpoint1,xvector<double>& kpoint2,double tol){
  if(aurostd::identical(kpoint1,kpoint2,tol)){return true;}  //we can be dam sure at this tolerance
  return false;
}

bool xOUTCAR::removeDuplicateKPoints(vector<xvector<double> >& vkpoints,vector<uint>& vikpt){
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  vector<uint> vikpt_unique;
  bool found;
  for(uint i=0;i<vikpt.size();i++){
    if(vikpt[i]>vkpoints.size()-1){return false;}
    found=false;
    for(uint j=0;j<vikpt_unique.size()&&!found;j++){
      if(identicalKPoints(vkpoints,vikpt[i],vikpt_unique[j])){
        if(LDEBUG) {
          cerr << XPID << "xOUTCAR::removeDuplicateKPoints: removing duplicate k-point ";
          cerr << vkpoints[vikpt[i]] << " (kpt=" << vikpt[i] << ")" << endl;
        }
        found=true;
      }
    }
    if(!found){vikpt_unique.push_back(vikpt[i]);}
  }
  vikpt=vikpt_unique;
  return true;
}

bool xOUTCAR::removeDuplicateKPoints(vector<vector<xvector<double> > >& vkpoints,vector<uint>& vikpt,vector<uint>& vispin){
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  if(vikpt.size()!=vispin.size()){return false;}  //test of stupidity
  vector<uint> vikpt_unique,vispin_unique;
  bool found;
  for(uint i=0;i<vikpt.size();i++){
    if(vikpt[i]>vkpoints[vispin[i]].size()-1){return false;}
    found=false;
    for(uint j=0;j<vikpt_unique.size()&&!found;j++){
      if(identicalKPoints(vkpoints[vispin[i]][vikpt[i]],vkpoints[vispin_unique[j]][vikpt_unique[j]])){
        if(LDEBUG) {
          cerr << XPID << "xOUTCAR::removeDuplicateKPoints: removing duplicate k-point ";
          cerr << vkpoints[vispin[i]][vikpt[i]] << " (kpt=" << vikpt[i] << ",spin=" << vispin[i] << ")" << endl;
        }
        found=true;
      }
    }
    if(!found){vikpt_unique.push_back(vikpt[i]);vispin_unique.push_back(vispin[i]);}
  }
  vikpt=vikpt_unique;
  vispin=vispin_unique;
  return true;
}

double xOUTCAR::minimumDistanceKPoints(vector<xvector<double> >& vkpoints,uint ikp1,uint ikp2){
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  double dist_min=AUROSTD_MAX_DOUBLE;
  if((ikp1>vkpoints.size()-1) || (ikp2>vkpoints.size()-1)){return dist_min;}
  if(!xstr.atoms.size()){ //assume it's already loaded
    if(!GetXStructure()){
      if(LDEBUG) {cerr << XPID << "xOUTCAR::minimumDistanceKPoints: GetXStructure failed." << endl;}
      return dist_min;
    }
  }
  if(LDEBUG) {
    cerr << XPID << "xOUTCAR::minimumDistanceKPoints: subtracting ";
    cerr << "kpoint[" << ikp1 << "]";
    cerr << " from ";
    cerr << "kpoint[" << ikp2 << "]";
    cerr << endl;
  }
  //special case, identicalKPoints()
  if(ikp1==ikp2){
    dist_min=0.0;
    if(LDEBUG) {cerr << XPID << "xOUTCAR::minimumDistanceKPoints: Found identical k-points, distance in cartesian = " << dist_min << endl;}
    return dist_min;
  }
  return minimumDistanceKPoints(vkpoints[ikp1],vkpoints[ikp2]);
}

double xOUTCAR::minimumDistanceKPoints(xvector<double>& kpoint1_kl,xvector<double>& kpoint2_kl){  //these are in units of reciprocal lattice vectors
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  double dist_min=AUROSTD_MAX_DOUBLE;
  if(LDEBUG){
    cerr << XPID << "xOUTCAR::minimumDistanceKPoints: kpoint1_kl=" << kpoint1_kl << endl;
    cerr << XPID << "xOUTCAR::minimumDistanceKPoints: kpoint2_kl=" << kpoint2_kl << endl;
  }
  //special case, identicalKPoints()
  if(identicalKPoints(kpoint1_kl,kpoint2_kl)){
    dist_min=0.0;
    if(LDEBUG) {cerr << XPID << "xOUTCAR::minimumDistanceKPoints: Found identical k-points, distance in cartesian = " << dist_min << endl;}
    return dist_min;
  }
  xmatrix<double> klattice=ReciprocalLattice(xstr);  //repetita iuvant
  xmatrix<double> metric_tensor=MetricTensor(xstr);
  if(LDEBUG) {
    cerr << XPID << "xOUTCAR::minimumDistanceKPoints: metric tensor" << endl;
    cerr << metric_tensor << endl;
    cerr << XPID << "xOUTCAR::minimumDistanceKPoints: klattice" << endl;
    cerr << klattice << endl;
  }

  xmatrix<double> kf2c=trasp(klattice);
  xvector<double> kpoint1=kf2c*kpoint1_kl;  //units of 2pi/scale
  xvector<double> kpoint2=kf2c*kpoint2_kl;  //units of 2pi/scale
  if(LDEBUG) {
    cerr << XPID << "xOUTCAR::minimumDistanceKPoints: subtracting (units of 2*pi/scale) " << endl;
    cerr << "                                 kpoint=" << kpoint1 << endl;
    cerr << "                                 kpoint=" << kpoint2 << endl;
  }

  //make sure to account for skewed cells
  xvector<int> dims=LatticeDimensionSphere(klattice,RadiusSphereLattice(klattice));

  xvector<double> vdist_kcart,vdist_kcart_min;
  double dist=0.0;
  for(int i=-dims[1];i<=dims[1];i++){
    for(int j=-dims[2];j<=dims[2];j++){
      for(int k=-dims[3];k<=dims[3];k++){
        vdist_kcart=(kpoint1-kpoint2+
            double(i)*klattice(1)+double(j)*klattice(2)+double(k)*klattice(3));
        dist=aurostd::modulus(vdist_kcart);
        if(dist<dist_min){
          vdist_kcart_min=vdist_kcart;
          dist_min=dist;
        }
      }
    }
  }
  if(LDEBUG) {cerr << XPID << "xOUTCAR::minimumDistanceKPoints: distance vector in reciprocal space  = " << vdist_kcart_min << endl;}

  //convert distance to cartesian in real space
  xvector<double> vdist_dcart;
  double dist_dcart;
  vdist_dcart(1)=sum(metric_tensor(1)*vdist_kcart_min(1));  //remember metric_tensor is symmetric!
  vdist_dcart(2)=sum(metric_tensor(2)*vdist_kcart_min(2));  //remember metric_tensor is symmetric!
  vdist_dcart(3)=sum(metric_tensor(3)*vdist_kcart_min(3));  //remember metric_tensor is symmetric!
  vdist_dcart/=(2*pi);  //conversion back requires undoing 2pi
  if(LDEBUG) {cerr << XPID << "xOUTCAR::minimumDistanceKPoints: distance vector in real space        = " << vdist_dcart << endl;}
  dist_dcart=aurostd::modulus(vdist_dcart);
  if(LDEBUG) {cerr << XPID << "xOUTCAR::minimumDistanceKPoints: distance between kpoints (angstroms) = " << dist_dcart << endl;}
  return dist_dcart;
}

bool xOUTCAR::GetBandGap(double EFERMI,double efermi_tol,double energy_tol,double occ_tol) { //CO20171004 - redoing camilo's function more robustly

  //some nice examples when debugging - validated by CO20171006
  //aflow --bandgap=/common/ICSD/LIB/RHL/Ag1Ni1O2_ICSD_73974
  //System        :   Ag1Ni1O2_ICSD_73974
  //Spin tag      :   2
  //Fermi level   :  +2.9623e+00
  //                  VBT           CBB           Egap          Egap_fit     Type
  //Majority Spin :  -1.0000e+00   -1.0000e+00   -1.0000e+09   -1.0000e+09   metal
  //Minority Spin :  -8.7090e-01   +1.4751e+00   +2.3460e+00   +4.0754e+00   insulator-indirect
  //Net Result    :  -1.0000e+00   -1.0000e+00   -1.0000e+09   -1.0000e+09   half-metal
  //
  //aflow --bandgap=/common/ICSD/LIB/FCC/Ni1O1_ICSD_92133
  //System        :   Ni1O1_ICSD_92133
  //Spin tag      :   2
  //Fermi level   :  +5.7976e+00
  //                  VBT           CBB           Egap          Egap_fit     Type
  //Majority Spin :  +1.0100e-02   +2.1721e+00   +2.1620e+00   +3.8274e+00   insulator-indirect
  //Minority Spin :  -1.1203e+00   +1.8017e+00   +2.9220e+00   +4.8519e+00   insulator-indirect
  //Net Result    :  +1.0100e-02   +1.8017e+00   +1.7916e+00   +3.3281e+00   insulator-indirect_spin-polarized
  //
  //aflow --bandgap=/common/ICSD/LIB/FCC/Cr1O2_ICSD_186838
  //System        :   Cr1O2_ICSD_186838
  //Spin tag      :   2
  //Fermi level   :  +5.1970e+00
  //                  VBT           CBB           Egap          Egap_fit     Type
  //Majority Spin :  -1.0000e+00   -1.0000e+00   -1.0000e+09   -1.0000e+09   metal
  //Minority Spin :  -4.8020e-01   +2.1972e+00   +2.6774e+00   +4.5221e+00   insulator-indirect
  //Net Result    :  -1.0000e+00   -1.0000e+00   -1.0000e+09   -1.0000e+09   half-metal
  //
  //aflow --bandgap=/common/ICSD/LIB/FCC/Ba2Dy1Nb1O6_ICSD_109156
  //System        :   Ba2Dy1Nb1O6_ICSD_109156
  //Spin tag      :   1
  //Fermi level   :  +2.8134e+00
  //                  VBT           CBB           Egap          Egap_fit     Type
  //Net Result    :  -6.1000e-03   +3.0853e+00   +3.0914e+00   +5.0802e+00   insulator-direct
  //
  //aflow --bandgap=/common/ICSD/LIB/FCC/Si1_ICSD_150530
  //System        :   Si1_ICSD_150530
  //Spin tag      :   1
  //Fermi level   :  +5.6171e+00
  //                  VBT           CBB           Egap          Egap_fit     Type
  //Net Result    :  +7.2000e-03   +6.1720e-01   +6.1000e-01   +1.7353e+00   insulator-indirect
  //
  //aflow --bandgap=/common/ICSD/LIB/BCC/Fe1_ICSD_52258
  //System        :   Fe1_ICSD_52258
  //Spin tag      :   2
  //Fermi level   :  +4.8034e+00
  //                  VBT           CBB           Egap          Egap_fit     Type
  //Majority Spin :  -1.0000e+00   -1.0000e+00   -1.0000e+09   -1.0000e+09   metal
  //Minority Spin :  -1.0000e+00   -1.0000e+00   -1.0000e+09   -1.0000e+09   metal
  //Net Result    :  -1.0000e+00   -1.0000e+00   -1.0000e+09   -1.0000e+09   metal
  //
  //aflow --bandgap=/common/ICSD/LIB/HEX/H2_ICSD_68271
  //WARNING - xOUTCAR::GetBandGap: unable to detect band edge for k-point 1 (spin=2), this system should be rerun with a wider energy range
  //
  //   filename=[stringstream]
  //
  //   pflow::BANDGAP: . failed
  //
  //[OBSOLETE]aflow --bandgap=/common/ICSD/LIB/HEX/H2_ICSD_68271 (SHOULD GIVE WARNING NOW)
  //[OBSOLETE]System        :   H2_ICSD_68271
  //[OBSOLETE]Spin tag      :   2
  //[OBSOLETE]Fermi level   :  -5.6146e+00
  //[OBSOLETE]                  VBT           CBB           Egap          Egap_fit     Type
  //[OBSOLETE]Majority Spin :  -2.4930e-01   +7.1586e+00   +7.4079e+00   +1.0899e+01   insulator-indirect
  //[OBSOLETE]Minority Spin :  -1.0000e+00   -1.0000e+00   -1.0000e+09   -1.0000e+09   empty
  //[OBSOLETE]Net Result    :  -2.4930e-01   +7.1586e+00   +7.4079e+00   +1.0899e+01   insulator-indirect
  //
  //aflow --bandgap=/common/ICSD/LIB/HEX/Br2Mn1_ICSD_60250
  //System        :   Br2Mn1_ICSD_60250
  //Spin tag      :   2
  //Fermi level   :  +1.8650e-01
  //                  VBT           CBB           Egap          Egap_fit     Type
  //Majority Spin :  -1.3600e-02   +3.3703e+00   +3.3839e+00   +5.4745e+00   insulator-indirect
  //Minority Spin :  -5.8990e-01   +3.4065e+00   +3.9964e+00   +6.3001e+00   insulator-indirect
  //Net Result    :  -1.3600e-02   +3.3703e+00   +3.3839e+00   +5.4745e+00   insulator-indirect

  //repetita iuvant!!!!!!!!
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream message;
  bool force_exit=false; //[CO20210621 - we will now exit gracefully]XHOST.POSTPROCESS; //SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only

  if((content == "") || (vcontent.size() == 0)) {
    //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " xOUTCAR needs to be loaded before. \n"
    //[CO20200404 - OBSOLETE]  "        GetProperties(const stringstream&); \n"
    //[CO20200404 - OBSOLETE]  "        GetProperties(const string&); \n"
    //[CO20200404 - OBSOLETE]  "        GetPropertiesFile(const string&); \n";
    //[CO20200404 - OBSOLETE]return FALSE;
    message << "xOUTCAR needs to be loaded before." << endl;
    message << "GetProperties(const stringstream&);" << endl;
    message << "GetProperties(const string&);" << endl;
    message << "GetPropertiesFile(const string&);" << endl;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _INPUT_MISSING_);
  }
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " OUTCAR content found" << endl;}

  //quick check if GetProperties() failed
  if(!(ISPIN==1 || ISPIN==2) || NKPTS==0 || NBANDS==0){
    if(!GetProperties(content) || !(ISPIN==1 || ISPIN==2) || NKPTS==0 || NBANDS==0){
      //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " GetProperties failed. \n";
      //[CO20200404 - OBSOLETE]return false;
      message << "GetProperties() failed";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
  }
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " OUTCAR properties retrieved" << endl;}

  // GET FERMI LEVEL
  if(EFERMI==AUROSTD_NAN){
    //VERY IMPORTANT - CO20171009 from discussion with SC
    //we strongly prefer to use Efermi from OUTCAR.static, not OUTCAR.bands
    //bands is not self-consistent (ICHARG=11), only used to determine energy states
    //since it is not self-consistent, it does not fill in states correctly
    //static is the best way to go!

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " Using Efermi from CURRENT OUTCAR (recommended to use Efermi from OUTCAR.static)" << endl;}
    EFERMI=Efermi;
  }

  // ----------------------------------------------------------------------
  // GET POSCAR INFO FOR GETBANDGAP()
  // silly to load in POSCAR when all the information is in OUTCAR
  if(!xstr.atoms.size()){ //assume it's already loaded
    if(!GetXStructure()){
      //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " GetXStructure failed. \n";
      //[CO20200404 - OBSOLETE]return false;
      message << "GetXStructure() failed";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
  }
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " xstructure built from OUTCAR" << endl;}

  double kpt_tol;
  if(xstr.CalculateSymmetry()){kpt_tol=xstr.sym_eps;}
  else {kpt_tol=SYM::defaultTolerance(xstr);}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " kpt_tol=" << kpt_tol << endl;}
  // ----------------------------------------------------------------------

  vector<uint> starting_lines;
  if(!GetStartingKPointLines(starting_lines)){
    //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " Unable to grab starting k-point lines. \n";
    //[CO20200404 - OBSOLETE]return false;
    message << "Unable to grab starting k-point lines";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _FILE_CORRUPT_);
  }
  if(LDEBUG) {
    cerr << __AFLOW_FUNC__ << " reading k-points starting at line";
    cerr << (starting_lines.size()==2 ? "s: " : ": ");
    cerr << starting_lines[0];
    cerr << (starting_lines.size()==2 ? " and "+aurostd::utype2string(starting_lines[1]) : "");
    cerr << endl;
  }

  uint iline;
  vector<double> b_energies, b_occs;
  uint iedge;
  xvector<double> kpoint; //reciprocal

  //resize for spin
  vector<vector<xvector<double> > > vkpoints; vkpoints.resize(ISPIN); //reciprocal
  vector<vector<uint> > vedges; vedges.resize(ISPIN);
  vector<vector<double> > vVBT; vVBT.resize(ISPIN);
  vector<vector<double> > vCBB; vCBB.resize(ISPIN);
  vector<uint> empty_channel; empty_channel.resize(ISPIN);
  vector<uint> partially_empty_channel; partially_empty_channel.resize(ISPIN);

  bool band_edge_found;
  int kpt_found;
  int first_kpt_empty,first_spin_empty;
  first_kpt_empty=first_spin_empty=-1;

  //initialize empty_channel to 0's
  for(uint i=0;i<(uint)ISPIN;i++){empty_channel[i]=0;}

  for(uint ispin=0;ispin<(uint)ISPIN;ispin++){
    iline=starting_lines[ispin];
    for(uint ikpt=1;ikpt<(uint)NKPTS+1;ikpt++){
      //is OUTCAR poorly written?
      kpt_found=isKPointLine(iline,kpoint);
      if((int)ikpt!=kpt_found){
        if(kpt_found==-1 && ikpt<1000){  //ONLY exception, above 999, VASP starting writing ***, so ignore this check starting at 1000
          //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " missing k-point "+aurostd::utype2string(ikpt)+
          //[CO20200404 - OBSOLETE]  " (spin="+aurostd::utype2string(ispin+1)+") \n";
          //[CO20200404 - OBSOLETE]return false;
          message << "Missing k-point " << ikpt << " (spin=" << ispin+1 << ")";
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _FILE_CORRUPT_);
        }
      }
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " looking at k-point[" << ikpt << "]=" << kpoint << " (spin=" << ispin+1 << ")" << endl;}

      /////////////////////////////////////////////
      //BANDS analysis - START
      vkpoints[ispin].push_back(kpoint);  //push back kpoints even if empty channel

      if(!ProcessKPoint(iline,EFERMI,b_energies,b_occs)){
        //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " unable to process k-point "+aurostd::utype2string(ikpt)+
        //[CO20200404 - OBSOLETE]  " (spin="+aurostd::utype2string(ispin+1)+") \n";
        //[CO20200404 - OBSOLETE]return false;
        message << "Unable to process k-point " << ikpt << " (spin=" << ispin+1 << +")";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _FILE_CORRUPT_);
      }

      //keep vedges, vVBT, vCBB same size as vkpoints
      vedges[ispin].push_back(AUROSTD_MAX_UINT);
      vVBT[ispin].push_back(AUROSTD_MAX_DOUBLE);
      vCBB[ispin].push_back(AUROSTD_MAX_DOUBLE);

      band_edge_found=GetBandEdge(b_energies,b_occs,EFERMI,iedge,efermi_tol,energy_tol,occ_tol);

      if(band_edge_found){  //normal case
        if(empty_channel[ispin]==1){  
          //if we get here, empty channel detected for k-point 1, but not for current k-point
          //therefore, 1 is the outlier
          //inconsistency within channel
          partially_empty_channel[ispin]=1;
          if(first_kpt_empty==-1){first_kpt_empty=1;}
          if(first_spin_empty==-1){first_spin_empty=ispin;}
          //ERROR = __AFLOW_FUNC__ + " unable to detetc band edge for k-point "+aurostd::utype2string(1)+
          //  " (spin="+aurostd::utype2string(ispin+1)+") \n";
          //return false;
        }//[CO20210621 - always save good vedges, even if empty channel] else {
        vedges[ispin].back()=iedge;
        vVBT[ispin].back()=b_energies[iedge];
        vCBB[ispin].back()=b_energies[iedge+1];
        //[CO20210621 - always save good vedges, even if empty channel]}
      }else{  //special case
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " no band edge found for k-point[" << ikpt << "]=" << kpoint << " (spin=" << ispin+1 << ")" << endl;}
        if(ikpt!=1){  
          //if we get here, kpt==1 presents an edge, but current k-point does not
          //therefore, current k-point is the outlier
          //inconsistency within channel
          partially_empty_channel[ispin]=1;
          //ERROR = __AFLOW_FUNC__ + " unable to detect band edge for k-point "+aurostd::utype2string(ikpt)+
          //  " (spin="+aurostd::utype2string(ispin+1)+") \n";
          //return false;
        }
        empty_channel[ispin]=1;
        if(first_kpt_empty==-1){first_kpt_empty=ikpt;}
        if(first_spin_empty==-1){first_spin_empty=ispin;}
      }

      //BANDS analysis - STOP
      /////////////////////////////////////////////

      //can we find the next kpoint?
      if(ikpt<(uint)NKPTS && !GetNextKPointLine(iline)){
        //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " missing k-point "+aurostd::utype2string(ikpt+1)+
        //[CO20200404 - OBSOLETE]  " (spin="+aurostd::utype2string(ispin+1)+") \n";
        //[CO20200404 - OBSOLETE]return false;
        message << "Missing k-point " << ikpt+1 << " (spin=" << ispin+1 << +")";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _FILE_CORRUPT_);
      }
    }
  }

  //CO20171007 - this code can handle empty types (returning empty(-partial))
  //but after talking with SC, I decided to simply return false!
  //these runs should be rerun - they are garbage!
  if(0){  //CO20210621 - moved below
    if(empty_channel[0]==1 || (ISPIN==2 && empty_channel[1]==1)){
      //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " unable to detect band edge for k-point "+aurostd::utype2string(first_kpt_empty)+
      //[CO20200404 - OBSOLETE]  " (spin="+aurostd::utype2string(first_spin_empty+1)+"),"+
      //[CO20200404 - OBSOLETE]  " this system should be rerun with a wider energy range \n";
      message << "Unable to detect band edge for k-point " << first_kpt_empty << " (spin=" << first_spin_empty+1 << ")," << " this system should be rerun with a wider energy range";
      if(force_exit){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);}
      else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);return false;}
    }
  }

  //first, we differentiate between metals / insulators
  //in metals, bands cross at the edge
  //this manifests as a change in the edge index
  //think of this as a band "appearing"/"disappearing" at the edge
  //or a "conservation" of bands below edge, lack of conservation means crossing, or a metal
  //enum BROAD_TYPES {empty,metal,insulator};
  //enum EMPTY_TYPES {empty_all,empty_partial};
  //enum INSULATOR_TYPES {insulator_direct,insulator_indirect};
  //enum GAP_TYPES {zero_gap,non_zero_gap};

  //resize for spin
  vector<BROAD_TYPES> broad_type; broad_type.resize(ISPIN);
  vector<EMPTY_TYPES> empty_type; empty_type.resize(ISPIN);
  vector<INSULATOR_TYPES> insulator_type; insulator_type.resize(ISPIN);
  vector<GAP_TYPES> gap_type; gap_type.resize(ISPIN);
  vector<double> gap; gap.resize(ISPIN);
  vector<uint> vimax_VBTs;vimax_VBTs.resize(ISPIN); //for each spin
  vector<uint> vimin_CBBs;vimin_CBBs.resize(ISPIN); //for each spin

  //broad_type
  for(uint ispin=0;ispin<(uint)ISPIN;ispin++){
    //[CO20210621 - check for metal first, even if empty]if(empty_channel[ispin]){broad_type[ispin]=empty;continue;}
    broad_type[ispin]=insulator;  //insulator unless proven otherwise
    for(uint ikpt=1;ikpt<vedges[ispin].size();ikpt++){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking at vedges[ispin=" << ispin << "][ikpt-1=" << ikpt-1 << "]=" << vedges[ispin][ikpt-1] << " vs. " << vedges[ispin][ikpt] << "=vedges[ispin=" << ispin << "][ikpt=" << ikpt << "]" << endl;}
      if(vedges[ispin][ikpt-1]==AUROSTD_MAX_UINT||vedges[ispin][ikpt]==AUROSTD_MAX_UINT){continue;}
      if(vedges[ispin][ikpt-1]!=vedges[ispin][ikpt]){
        broad_type[ispin]=metal; //metal
        break;
      }
    }
    if(empty_channel[ispin] && broad_type[ispin]!=metal){broad_type[ispin]=empty;}
  }

  for(uint ispin=0;ispin<(uint)ISPIN;ispin++){
    if(broad_type[ispin]==empty){
      message << "Unable to detect band edge for k-point " << first_kpt_empty << " (spin=" << first_spin_empty+1 << ")," << " this system should be rerun with a wider energy range";
      if(force_exit){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);}
      else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);return false;}
    }
  }

  double max_VBT=0.0,min_CBB=0.0;
  double dist=0.0,dist_min=0.0;
  uint imax_VBT=0,imin_CBB=0; //absolutes, then our final picks
  vector<uint> vimax_VBTs_equiv,vimin_CBBs_equiv; //within tolerance of absolutes
  bool vbt_duplicate_remove,cbb_duplicate_remove;

  //specific types here
  for(uint ispin=0;ispin<(uint)ISPIN;ispin++){
    if(broad_type[ispin]==empty){
      gap[ispin]=0.0;
      if(partially_empty_channel[ispin]==1){
        empty_type[ispin]=empty_partial;
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " empty (partial) channel found spin=" << ispin+1 << endl;}
      } else {
        empty_type[ispin]=empty_all;
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " empty (all) channel found spin=" << ispin+1 << endl;}
      }
      continue;
    }
    //metal, only one type
    //modify starting here in the future for SEMI-metals (need DOS)
    else if(broad_type[ispin]==metal){
      gap[ispin]=0.0;
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " metal found in spin=" << ispin+1 << endl;}
      continue;
    }
    else if(broad_type[ispin]!=insulator){  //test of stupidity
      //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " unknown material type (!empty && !metal && !insulator) (spin="+
      //[CO20200404 - OBSOLETE]  aurostd::utype2string(ispin+1)+") \n";
      //[CO20200404 - OBSOLETE]return false;
      message << "Unknown material type (!empty && !metal && !insulator) (spin=" << ispin+1 << ")";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " insulator found in spin=" << ispin+1 << endl;}
    max_VBT = (-1.0) * AUROSTD_MAX_DOUBLE;
    min_CBB =          AUROSTD_MAX_DOUBLE;
    //determine if direct/indirect insulator
    //first, get absolute max of VBT/min of CBB
    for(uint ikpt=0;ikpt<vkpoints[ispin].size();ikpt++){
      if(vVBT[ispin][ikpt]>max_VBT){
        imax_VBT=ikpt;
        max_VBT=vVBT[ispin][imax_VBT];
      }
      if(vCBB[ispin][ikpt]<min_CBB){
        imin_CBB=ikpt;
        min_CBB=vCBB[ispin][imin_CBB];
      }
    }
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " absolute VBT=" << max_VBT << " at k-point " << vkpoints[ispin][imax_VBT] << " (kpt=" << imax_VBT << ",spin=" << ispin+1 << ")" << endl;
      cerr << __AFLOW_FUNC__ << " absolute CBB=" << min_CBB << " at k-point " << vkpoints[ispin][imin_CBB] << " (kpt=" << imin_CBB << ",spin=" << ispin+1 << ")" << endl;
    }
    vimax_VBTs_equiv.clear(); vimin_CBBs_equiv.clear();
    //grab any equivalently high/low extrema
    for(uint ikpt=0;ikpt<vkpoints[ispin].size();ikpt++){
      if(abs(max_VBT-vVBT[ispin][ikpt])<energy_tol){
        vimax_VBTs_equiv.push_back(ikpt);
        if(LDEBUG && ikpt!=imax_VBT){cerr << __AFLOW_FUNC__ << " found equivalent VBT at " << vkpoints[ispin][ikpt] << " (kpt=" << ikpt << ",spin=" << ispin+1 << ")" << endl;}
      }
      if(abs(min_CBB-vCBB[ispin][ikpt])<energy_tol){
        vimin_CBBs_equiv.push_back(ikpt);
        if(LDEBUG && ikpt!=imin_CBB){cerr << __AFLOW_FUNC__ << " found equivalent CBB at " << vkpoints[ispin][ikpt] << " (kpt=" << ikpt << ",spin=" << ispin+1 << ")" << endl;}
      }
    }
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " removing duplicate k-points from VBT search (spin=" << ispin+1 << ")" << endl;}
    vbt_duplicate_remove=removeDuplicateKPoints(vkpoints[ispin],vimax_VBTs_equiv);
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " removing duplicate k-points from CBB search (spin=" << ispin+1 << ")" << endl;}
    cbb_duplicate_remove=removeDuplicateKPoints(vkpoints[ispin],vimin_CBBs_equiv);
    if(!vbt_duplicate_remove || !cbb_duplicate_remove){
      //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " cannot find equivalent band extrema (spin="+
      //[CO20200404 - OBSOLETE]  aurostd::utype2string(ispin+1)+") \n";
      message << "Cannot find equivalent band extrema (spin=" << ispin+1 << ")";
      if(force_exit){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);}
      else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);return false;}
    }
    //if we found more than one possible extrema, minimize kpoint distance between max/min
    //this simulates the electron trying to reduce momentum needed to conduct
    if(vimax_VBTs_equiv.size()==1 && vimin_CBBs_equiv.size()==1){  //easy case, keep already defined imax/imin
      dist_min=minimumDistanceKPoints(vkpoints[ispin],imax_VBT,imin_CBB);
    }
    else { //find minimum distance pairs of max/min
      dist_min=AUROSTD_MAX_DOUBLE;
      for(uint imax=0;imax<vimax_VBTs_equiv.size();imax++){
        for(uint imin=0;imin<vimin_CBBs_equiv.size();imin++){
          dist=minimumDistanceKPoints(vkpoints[ispin],vimax_VBTs_equiv[imax],vimin_CBBs_equiv[imin]);
          if(dist<dist_min){
            dist_min=dist;
            imax_VBT=vimax_VBTs_equiv[imax];
            imin_CBB=vimin_CBBs_equiv[imin];
          }
        }
      }
    }
    if(dist_min==AUROSTD_MAX_DOUBLE){
      //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " cannot calculate k-point distance between band extrema (spin="+
      //[CO20200404 - OBSOLETE]  aurostd::utype2string(ispin+1)+") \n";
      //[CO20200404 - OBSOLETE]return false;
      message << "Cannot calculate k-point distance between band extrema (spin=" << ispin+1 << ")";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " dist_min=" << dist_min << endl;}
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " vCBB[ispin=" << ispin << "][imin_CBB=" << imin_CBB << "]=" << vCBB[ispin][imin_CBB] << endl;
      cerr << __AFLOW_FUNC__ << " vVBT[ispin=" << ispin << "][imax_VBT=" << imax_VBT << "]=" << vVBT[ispin][imax_VBT] << endl;
    }
    gap[ispin]=vCBB[ispin][imin_CBB]-vVBT[ispin][imax_VBT];
    vimin_CBBs[ispin]=imin_CBB; //save for later spin loops
    vimax_VBTs[ispin]=imax_VBT; //save for later spin loops
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " gap[ispin=" << ispin << "]=" << gap[ispin] << endl;}
    if(gap[ispin]<0){ //test of stupidity, this should NEVER happen
      //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " negative band gap found, something broke (spin="+
      //[CO20200404 - OBSOLETE]  aurostd::utype2string(ispin+1)+") \n";
      //[CO20200404 - OBSOLETE]return false;
      message << "Negative band gap found, something broke (spin=" << ispin+1 << ")";
      if(force_exit){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);}
      else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);return false;}
    }
    gap_type[ispin]=(gap[ispin] < energy_tol ? zero_gap : non_zero_gap);
    insulator_type[ispin]=(abs(dist_min)<kpt_tol ? insulator_direct : insulator_indirect);
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ";
      cerr << (insulator_type[ispin]==insulator_indirect ? "IN" : "" ) << "DIRECT insulator ";
      cerr << "gap=" << gap[ispin] << " ";
      cerr << (gap_type[ispin]==zero_gap ? "(zero-gap) " : "");
      cerr << "found in spin=" << ispin+1 << endl;
    }
  }

  //get spin-channel specific properties
  valence_band_max.clear();     valence_band_max.resize(ISPIN);
  conduction_band_min.clear();  conduction_band_min.resize(ISPIN);
  Egap.clear();                 Egap.resize(ISPIN);
  Egap_type.clear();            Egap_type.resize(ISPIN);
  Egap_fit.clear();             Egap_fit.resize(ISPIN);

  //[CO20191004]double _METALGAP = -AUROSTD_NAN, _METALEDGE = -1.0;
  for(uint ispin=0;ispin<(uint)ISPIN;ispin++){
    if(broad_type[ispin]==empty){
      valence_band_max[ispin]=_METALEDGE_;
      conduction_band_min[ispin]=_METALEDGE_;
      Egap[ispin]=_METALGAP_;
      Egap_type[ispin]="empty"+string(empty_type[ispin]==empty_partial ? "-partially" : "");
    } else if(broad_type[ispin]==metal){
      valence_band_max[ispin]=_METALEDGE_;
      conduction_band_min[ispin]=_METALEDGE_;
      Egap[ispin]=_METALGAP_;
      Egap_type[ispin]="metal";
    } else if(broad_type[ispin]!=insulator){ //test of stupidity
      //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " unknown material type (!empty && !metal && !insulator) (spin="+
      //[CO20200404 - OBSOLETE]  aurostd::utype2string(ispin+1)+") \n";
      //[CO20200404 - OBSOLETE]return false;
      message << "Unknown material type (!empty && !metal && !insulator) (spin=" << ispin+1 << ")";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    } else {
      valence_band_max[ispin]=vVBT[ispin][vimax_VBTs[ispin]];
      conduction_band_min[ispin]=vCBB[ispin][vimin_CBBs[ispin]];
      Egap[ispin]=(gap_type[ispin]==zero_gap ? 0.0 : gap[ispin]);
      Egap_type[ispin]="insulator-"+
        string(insulator_type[ispin]==insulator_indirect ? "in" : "")+"direct"+
        string(gap_type[ispin]==zero_gap ? "_zero-gap" : "");
    }
    Egap_fit[ispin]=(broad_type[ispin]==empty || broad_type[ispin]==metal ? _METALGAP_ : 1.348 * Egap[ispin] + 0.913);
  }

  //get system-wide properties
  if(broad_type[0]==empty || (ISPIN==2 && broad_type[1]==empty)){
    if(ISPIN==1 || (ISPIN==2 && broad_type[0]==empty && broad_type[1]==empty)){
      valence_band_max_net=_METALEDGE_;
      conduction_band_min_net=_METALEDGE_;
      Egap_net=_METALGAP_;
      Egap_type_net="empty"+string(empty_type[0]==empty_partial || (ISPIN==2 && empty_type[1]==empty_partial) ? "-partially" : "");;
    } else {
      if(broad_type[0]==empty){ //grab properties of spin-channel 2
        uint ispin=1;
        valence_band_max_net=(broad_type[ispin]==metal ? _METALEDGE_ : vVBT[ispin][vimax_VBTs[ispin]]);
        conduction_band_min_net=(broad_type[ispin]==metal ? _METALEDGE_ : vCBB[ispin][vimin_CBBs[ispin]]);
        Egap_net=(broad_type[ispin]==metal ? _METALGAP_ : gap[ispin]);
        Egap_type_net=(broad_type[ispin]==metal ? "metal" : "insulator-"+
            string(insulator_type[ispin]==insulator_indirect ? "in" : "")+"direct"+
            string(gap_type[ispin]==zero_gap ? "_zero-gap" : ""));
      } else {  //broad_type[1]==empty, grab properties of spin-channel 1
        uint ispin=0;
        valence_band_max_net=(broad_type[ispin]==metal ? _METALEDGE_ : vVBT[ispin][vimax_VBTs[ispin]]);
        conduction_band_min_net=(broad_type[ispin]==metal ? _METALEDGE_ : vCBB[ispin][vimin_CBBs[ispin]]);
        Egap_net=(broad_type[ispin]==metal ? _METALGAP_ : gap[ispin]);
        Egap_type_net=(broad_type[ispin]==metal ? "metal" : "insulator-"+
            string(insulator_type[ispin]==insulator_indirect ? "in" : "")+"direct"+
            string(gap_type[ispin]==zero_gap ? "_zero-gap" : ""));
      }
    }
  } else if(!(ISPIN==2 && !(broad_type[0]==metal && broad_type[0]==broad_type[1]))){ //easy cases
    if(ISPIN==1){  //all ISPIN==1 come here
      uint ispin=0;
      valence_band_max_net=(broad_type[ispin]==metal ? _METALEDGE_ : vVBT[ispin][vimax_VBTs[ispin]]);
      conduction_band_min_net=(broad_type[ispin]==metal ? _METALEDGE_ : vCBB[ispin][vimin_CBBs[ispin]]);
      Egap_net=(broad_type[ispin]==metal ? _METALGAP_ : gap[ispin]);
      Egap_type_net=(broad_type[ispin]==metal ? "metal" : "insulator-"+
          string(insulator_type[ispin]==insulator_indirect ? "in" : "")+"direct"+
          string(gap_type[ispin]==zero_gap ? "_zero-gap" : ""));
    } else { //special case ISPIN==2 where both are metallic
      valence_band_max_net=_METALEDGE_;
      conduction_band_min_net=_METALEDGE_;
      Egap_net=_METALGAP_;
      Egap_type_net="metal";
    }
  } else if(broad_type[0]==metal || broad_type[1]==metal) {  //special case ISPIN==2, one metallic and one insulating spin-channel means the whole system is half-metallic
    valence_band_max_net=_METALEDGE_;
    conduction_band_min_net=_METALEDGE_;
    Egap_net=_METALGAP_;
    Egap_type_net="half-metal";
  } else if(!(broad_type[0]==insulator && broad_type[1]==insulator)){  //test of stupidity
    //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " unknown material type (!insulator_spin-polarized) (spin-averaged) \n";
    //[CO20200404 - OBSOLETE]return false;
    message << "Unknown material type (!insulator_spin-polarized) (spin-averaged)";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
  } else {  //do more work for 2 spin channels that are both insulating
    //get spin-averaged properties
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " we have a spin-polarized insulator, need to find spin-averaged gap" << endl;}
    double max_VBT_net,min_CBB_net;
    max_VBT_net = (-1.0) * AUROSTD_MAX_DOUBLE;
    min_CBB_net =          AUROSTD_MAX_DOUBLE;
    uint imax_VBT_net = 0, imin_CBB_net = 0, ispin_VBT_net = 0, ispin_CBB_net = 0;
    vector<uint> vimax_VBTs_equiv_net,vimin_CBBs_equiv_net; //within tolerance of absolutes
    vector<uint> vispin_VBTs_net,vispin_CBBs_net;
    //first, get absolute max of VBT/min of CBB
    //these will be the max(vVBT[ispin]) and min(vCBB[ispin]) across spins from before
    for(uint ispin=0;ispin<(uint)ISPIN;ispin++){
      for(uint ikpt=0;ikpt<vkpoints[ispin].size();ikpt++){
        if(vVBT[ispin][ikpt]>max_VBT_net){
          imax_VBT_net=ikpt;ispin_VBT_net=ispin;
          max_VBT_net=vVBT[ispin_VBT_net][imax_VBT_net];
        }
        if(vCBB[ispin][ikpt]<min_CBB_net){
          imin_CBB_net=ikpt;ispin_CBB_net=ispin;
          min_CBB_net=vCBB[ispin_CBB_net][imin_CBB_net];
        }
      }
    }
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " absolute VBT_net=" << max_VBT_net << " at k-point " << vkpoints[ispin_VBT_net][imax_VBT_net] << " (kpt=" << imax_VBT_net << ",spin=" << ispin_VBT_net << ")" << endl;
      cerr << __AFLOW_FUNC__ << " absolute CBB_net=" << min_CBB_net << " at k-point " << vkpoints[ispin_VBT_net][imin_CBB_net] << " (kpt=" << imin_CBB_net << ",spin=" << ispin_CBB_net << ")" << endl;
    }
    vimax_VBTs_equiv_net.clear(); vispin_VBTs_net.clear(); vimin_CBBs_equiv_net.clear(); vispin_CBBs_net.clear();
    //grab any equivalently high/low extrema
    for(uint ispin=0;ispin<(uint)ISPIN;ispin++){
      for(uint ikpt=0;ikpt<vkpoints[ispin].size();ikpt++){
        if(abs(max_VBT_net-vVBT[ispin][ikpt])<energy_tol){
          vimax_VBTs_equiv_net.push_back(ikpt);
          vispin_VBTs_net.push_back(ispin);
          if(LDEBUG && ikpt!=imax_VBT_net){cerr << __AFLOW_FUNC__ << " found equivalent VBT_net at " << vkpoints[ispin][ikpt] << " (kpt=" << ikpt << ",spin=" << ispin << ")" << endl;}
        }
        if(abs(min_CBB_net-vCBB[ispin][ikpt])<energy_tol){
          vimin_CBBs_equiv_net.push_back(ikpt);
          vispin_CBBs_net.push_back(ispin);
          if(LDEBUG && ikpt!=imin_CBB_net){cerr << __AFLOW_FUNC__ << " found equivalent CBB_net at " << vkpoints[ispin][ikpt] << " (kpt=" << ikpt << ",spin=" << ispin << ")" << endl;}
        }
      }
    }
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " removing duplicate k-points from VBT search (spin-averaged)" << endl;}
    vbt_duplicate_remove=removeDuplicateKPoints(vkpoints,vimax_VBTs_equiv_net,vispin_VBTs_net);
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " removing duplicate k-points from CBB search (spin-averaged)" << endl;}
    cbb_duplicate_remove=removeDuplicateKPoints(vkpoints,vimin_CBBs_equiv_net,vispin_CBBs_net);
    if(!vbt_duplicate_remove || !cbb_duplicate_remove){
      //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " cannot find equivalent band extrema (spin-averaged) \n";
      //[CO20200404 - OBSOLETE]return false;
      message << "Cannot find equivalent band extrema (spin-averaged)";
      if(force_exit){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);}
      else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);return false;}
    }
    //if we found more than one possible extrema, minimize kpoint distance between max/min
    //this simulates the electron trying to reduce momentum needed to conduct
    if(vimax_VBTs_equiv_net.size()==1 && vimin_CBBs_equiv_net.size()==1){  //easy case, keep already defined imax/imin
      dist_min=minimumDistanceKPoints(vkpoints[vispin_VBTs_net[0]][vimax_VBTs_equiv_net[0]],vkpoints[vispin_CBBs_net[0]][vimin_CBBs_equiv_net[0]]);
    }
    else { //find minimum distance pairs of max/min
      dist_min=AUROSTD_MAX_DOUBLE;
      for(uint imax=0;imax<vimax_VBTs_equiv_net.size();imax++){
        for(uint imin=0;imin<vimin_CBBs_equiv_net.size();imin++){
          dist=minimumDistanceKPoints(vkpoints[vispin_VBTs_net[imax]][vimax_VBTs_equiv_net[imax]],vkpoints[vispin_CBBs_net[imin]][vimin_CBBs_equiv_net[imin]]);
          if(dist<dist_min){
            dist_min=dist;
            imax_VBT_net=vimax_VBTs_equiv_net[imax];ispin_VBT_net=vispin_VBTs_net[imax];
            imin_CBB_net=vimin_CBBs_equiv_net[imin];ispin_CBB_net=vispin_CBBs_net[imin];
          }
        }
      }
    }
    if(dist_min==AUROSTD_MAX_DOUBLE){
      //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " cannot calculate k-point distance between band extrema (spin-averaged) \n";
      //[CO20200404 - OBSOLETE]return false;
      message << "Cannot calculate k-point distance between band extrema (spin-averaged)";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    valence_band_max_net=vVBT[ispin_VBT_net][imax_VBT_net];
    conduction_band_min_net=vCBB[ispin_CBB_net][imin_CBB_net];
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " vCBB[ispin_CBB_net=" << ispin_CBB_net << "][imin_CBB_net=" << imin_CBB_net << "]=" << vCBB[ispin_CBB_net][imin_CBB_net] << endl;
      cerr << __AFLOW_FUNC__ << " vVBT[ispin_VBT_net=" << ispin_VBT_net << "][imax_VBT_net=" << imax_VBT_net << "]" << vVBT[ispin_VBT_net][imax_VBT_net] << endl;
    }
    bool direct_insulator_net=abs(dist_min)<kpt_tol;
    double gap_net=(conduction_band_min_net-valence_band_max_net);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " gap_net=" << gap_net << endl;}
    if(gap_net<0){ //test of stupidity, this should NEVER happen
      //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " negative band gap found, something broke (spin-averaged) \n";
      //[CO20200404 - OBSOLETE]return false;
      message << "Negative band gap found, something broke (spin-averaged)";
      if(force_exit){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);}
      else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);return false;}
    }
    bool zero_gap_net=gap_net < energy_tol;
    bool spin_polarized_net=ispin_VBT_net!=ispin_CBB_net;
    Egap_net=(zero_gap_net ? 0.0 : gap_net);
    Egap_type_net="insulator-"+
      string(direct_insulator_net ? "" : "in" )+"direct"+
      string(zero_gap_net ? "_zero-gap" : "")+
      string(spin_polarized_net ? "_spin-polarized" : "");
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ";
      cerr << (direct_insulator_net ? "" : "IN" ) << "DIRECT insulator ";
      cerr << "gap=" << gap_net << " ";
      cerr << (zero_gap_net ? "(zero-gap) " : "");
      cerr << (spin_polarized_net ? "(spin_polarized) " : "");
      cerr << "spin-averaged" << endl;
    }
  }
  Egap_fit_net=(aurostd::substring2bool(Egap_type_net,"empty") || aurostd::substring2bool(Egap_type_net,"metal") ? _METALGAP_ : 1.348 * Egap_net + 0.913);

  return true;
}

//[CO20200404 - OBSOLETE]bool xOUTCAR::GetBandGap_Camilo(double kpt_tol) {
//[CO20200404 - OBSOLETE]  string line0, line1, line2, line3, line4, line5;
//[CO20200404 - OBSOLETE]  vector<string> tokens1, tokens2, tokens3, tokens4, tokens5;
//[CO20200404 - OBSOLETE]  vector<double> CBB, VBT;
//[CO20200404 - OBSOLETE]  xmatrix<double> KlatticeTmp(3,3), direct_lattice(3,3), metric_tensor(3,3);  //CO20171002
//[CO20200404 - OBSOLETE]  uint NBANDS, ISPIN, NKPT, NSW;
//[CO20200404 - OBSOLETE]  //bool c_nkpt=TRUE, c_nbands=TRUE, c_system=TRUE, c_ispin=TRUE, c_efermi=TRUE;
//[CO20200404 - OBSOLETE]  bool c_nkpt=TRUE, c_nbands=TRUE, c_ispin=TRUE, c_efermi=TRUE;
//[CO20200404 - OBSOLETE]  bool c_nsw=TRUE, c_cellvol=TRUE, c_dlattice=TRUE;
//[CO20200404 - OBSOLETE]  bool SPIN_UP = FALSE, SPIN_DN = FALSE;
//[CO20200404 - OBSOLETE]  //double  CELLVOL=0; CELLVOLCUTOFF=0,   //CO20171002 - GARBAGE A != A^-1
//[CO20200404 - OBSOLETE]  double EMIN = -100.0, EMAX = 100.0;
//[CO20200404 - OBSOLETE]  //[CO20191004]double _METALGAP = -1.0E+09, _METALEDGE = -1.0; 
//[CO20200404 - OBSOLETE]  double _ENER_EPS = 1.0E-08, _ELEC_EPS = 1.0E-02, _ZERO = 0.0;
//[CO20200404 - OBSOLETE]
//[CO20200404 - OBSOLETE]  // change this into an error message
//[CO20200404 - OBSOLETE]  if((content == "") or (vcontent.size() == 0)) {
//[CO20200404 - OBSOLETE]    ERROR = "xOUTCAR::GetBandGap(void): xOUTCAR needs to be loaded before. \n"
//[CO20200404 - OBSOLETE]      "        GetProperties(const stringstream&); \n"
//[CO20200404 - OBSOLETE]      "        GetProperties(const string&); \n"
//[CO20200404 - OBSOLETE]      "        GetPropertiesFile(const string&);";
//[CO20200404 - OBSOLETE]    return FALSE;
//[CO20200404 - OBSOLETE]    message << "xOUTCAR needs to be loaded before." << endl;
//[CO20200404 - OBSOLETE]    message << "GetProperties(const stringstream&);" << endl;
//[CO20200404 - OBSOLETE]    message << "GetProperties(const string&);" << endl;
//[CO20200404 - OBSOLETE]    message << "GetPropertiesFile(const string&);" << endl;
//[CO20200404 - OBSOLETE]    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _INPUT_MISSING_);
//[CO20200404 - OBSOLETE]  }
//[CO20200404 - OBSOLETE]
//[CO20200404 - OBSOLETE]  uint tagcount = 0;
//[CO20200404 - OBSOLETE]  // get only what is necessary for the Egap information
//[CO20200404 - OBSOLETE]  for(uint ii=0;  ii<vcontent.size(); ii++) {
//[CO20200404 - OBSOLETE]    if(c_nkpt) {
//[CO20200404 - OBSOLETE]      if(aurostd::substring2bool(vcontent.at(ii),"NKPTS")) {
//[CO20200404 - OBSOLETE]        line0=vcontent.at(ii);
//[CO20200404 - OBSOLETE]        aurostd::string2tokens(line0,tokens1);
//[CO20200404 - OBSOLETE]        NKPT = aurostd::string2utype<uint>(tokens1.at(3));
//[CO20200404 - OBSOLETE]        c_nkpt = FALSE;
//[CO20200404 - OBSOLETE]        tagcount++;
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]    }
//[CO20200404 - OBSOLETE]    if(c_nbands) {
//[CO20200404 - OBSOLETE]      if(aurostd::substring2bool(vcontent.at(ii),"NBANDS=")) {
//[CO20200404 - OBSOLETE]        line0=vcontent.at(ii);
//[CO20200404 - OBSOLETE]        aurostd::string2tokens(line0,tokens1);
//[CO20200404 - OBSOLETE]        NBANDS = aurostd::string2utype<uint>(tokens1.at(8));
//[CO20200404 - OBSOLETE]        c_nbands = FALSE;
//[CO20200404 - OBSOLETE]        tagcount++;
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]    }
//[CO20200404 - OBSOLETE]    // [OBSOLETE ] if(c_system) {
//[CO20200404 - OBSOLETE]    // [OBSOLETE ]   if(aurostd::substring2bool(vcontent.at(ii),"SYSTEM")) {
//[CO20200404 - OBSOLETE]    // [OBSOLETE ] 	line0=vcontent.at(ii);
//[CO20200404 - OBSOLETE]    // [OBSOLETE ] 	aurostd::string2tokens(line0,tokens1);
//[CO20200404 - OBSOLETE]    // [OBSOLETE ] 	if(tokens1.at(1) == "=") {
//[CO20200404 - OBSOLETE]    // [OBSOLETE ] 	  SYSTEM = aurostd::string2utype<string>(tokens1.at(2));
//[CO20200404 - OBSOLETE]    // [OBSOLETE ] 	  c_system = FALSE;
//[CO20200404 - OBSOLETE]    // [OBSOLETE ] 	  tagcount++;
//[CO20200404 - OBSOLETE]    // [OBSOLETE ] 	}
//[CO20200404 - OBSOLETE]    // [OBSOLETE ]   }
//[CO20200404 - OBSOLETE]    // [OBSOLETE ] }
//[CO20200404 - OBSOLETE]    if(c_ispin) {
//[CO20200404 - OBSOLETE]      if(aurostd::substring2bool(vcontent.at(ii),"ISPIN")) {
//[CO20200404 - OBSOLETE]        line0=vcontent.at(ii);
//[CO20200404 - OBSOLETE]        aurostd::string2tokens(line0,tokens1);
//[CO20200404 - OBSOLETE]        ISPIN = aurostd::string2utype<uint>(tokens1.at(2));
//[CO20200404 - OBSOLETE]        c_ispin = FALSE;
//[CO20200404 - OBSOLETE]        tagcount++;
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]    }
//[CO20200404 - OBSOLETE]    if(c_nsw) {
//[CO20200404 - OBSOLETE]      if(aurostd::substring2bool(vcontent.at(ii),"NSW")) {
//[CO20200404 - OBSOLETE]        line0=vcontent.at(ii);
//[CO20200404 - OBSOLETE]        aurostd::string2tokens(line0,tokens1);
//[CO20200404 - OBSOLETE]        NSW = aurostd::string2utype<uint>(tokens1.at(2));
//[CO20200404 - OBSOLETE]        if(NSW > 0) {
//[CO20200404 - OBSOLETE]          ERROR = "Number of ionic steps (NSW) > 0, exiting bandgap determination";
//[CO20200404 - OBSOLETE]          return FALSE;
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]        c_nsw = FALSE;
//[CO20200404 - OBSOLETE]        tagcount++;
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]    }
//[CO20200404 - OBSOLETE]    if(c_efermi) {
//[CO20200404 - OBSOLETE]      if(aurostd::substring2bool(vcontent.at(ii),"E-fermi")) {
//[CO20200404 - OBSOLETE]        line0=vcontent.at(ii);
//[CO20200404 - OBSOLETE]        aurostd::string2tokens(line0,tokens1);
//[CO20200404 - OBSOLETE]        Efermi = aurostd::string2utype<double>(tokens1.at(2));
//[CO20200404 - OBSOLETE]        c_efermi = FALSE;
//[CO20200404 - OBSOLETE]        tagcount++;
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]    }
//[CO20200404 - OBSOLETE]    if(c_cellvol) {
//[CO20200404 - OBSOLETE]      if(aurostd::substring2bool(vcontent.at(ii),"volume")) {
//[CO20200404 - OBSOLETE]        line0=vcontent.at(ii);
//[CO20200404 - OBSOLETE]        aurostd::string2tokens(line0,tokens1);
//[CO20200404 - OBSOLETE]        if((tokens1.at(1) == "of") and (tokens1.at(2) == "cell")) {
//[CO20200404 - OBSOLETE]          //CELLVOL = aurostd::string2utype<double>(tokens1.at(4));
//[CO20200404 - OBSOLETE]          //CELLVOLCUTOFF = 0.05*pow(CELLVOL,1.0/3.0);  //CO20171002 - GARBAGE A != A^-1
//[CO20200404 - OBSOLETE]          c_cellvol = FALSE;
//[CO20200404 - OBSOLETE]          tagcount++;
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]    }
//[CO20200404 - OBSOLETE]    //if(c_dlattice) //CO get most up to date (closer to bottom of file) lattice vectors
//[CO20200404 - OBSOLETE]    //if(aurostd::substring2bool(vcontent.at(ii),"Lattice")) //CO
//[CO20200404 - OBSOLETE]    if(aurostd::substring2bool(vcontent.at(ii),"direct"))
//[CO20200404 - OBSOLETE]    { //CO20200106 - patching for auto-indenting
//[CO20200404 - OBSOLETE]      line0=vcontent.at(ii);
//[CO20200404 - OBSOLETE]      aurostd::string2tokens(line0,tokens1);
//[CO20200404 - OBSOLETE]      //if((tokens1.at(0) == "Lattice") and (tokens1.at(1) == "vectors:")) //CO this is SYM-tag dependent
//[CO20200404 - OBSOLETE]      if((tokens1.at(0) == "direct") and (tokens1.at(1) == "lattice") and (tokens1.at(2) == "vectors") and
//[CO20200404 - OBSOLETE]          (tokens1.at(3) == "reciprocal") and (tokens1.at(4) == "lattice") and (tokens1.at(5) == "vectors")) //direct lattice vectors reciprocal lattice vectors
//[CO20200404 - OBSOLETE]      { //CO20200106 - patching for auto-indenting
//[CO20200404 - OBSOLETE]        //line1=vcontent.at(ii+2);    //CO
//[CO20200404 - OBSOLETE]        line1=vcontent.at(ii+1);
//[CO20200404 - OBSOLETE]        //aurostd::string2tokens(line1,tokens1);  //CO20170725
//[CO20200404 - OBSOLETE]        tokens1=GetCorrectPositions(line1,6); //CO20170725
//[CO20200404 - OBSOLETE]        if(!tokens1.size()){
//[CO20200404 - OBSOLETE]          ERROR = "line with lattice vector is ill-written, see: "+line1;
//[CO20200404 - OBSOLETE]          return false;
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]        //direct_lattice(1,1) = aurostd::string2utype<double>(tokens1.at(3).substr(0,tokens1.at(3).size()-1));    //CO
//[CO20200404 - OBSOLETE]        //direct_lattice(1,2) = aurostd::string2utype<double>(tokens1.at(4).substr(0,tokens1.at(4).size()-1));    //CO
//[CO20200404 - OBSOLETE]        //direct_lattice(1,3) = aurostd::string2utype<double>(tokens1.at(5).substr(0,tokens1.at(5).size()-1));    //CO
//[CO20200404 - OBSOLETE]        direct_lattice(1,1) = aurostd::string2utype<double>(tokens1.at(0));
//[CO20200404 - OBSOLETE]        direct_lattice(1,2) = aurostd::string2utype<double>(tokens1.at(1));
//[CO20200404 - OBSOLETE]        direct_lattice(1,3) = aurostd::string2utype<double>(tokens1.at(2));
//[CO20200404 - OBSOLETE]        //CO just grab reciprocal lattice vectors as calculated by vasp, faster than calculating yourself
//[CO20200404 - OBSOLETE]        KlatticeTmp(1,1) = aurostd::string2utype<double>(tokens1.at(3));
//[CO20200404 - OBSOLETE]        KlatticeTmp(1,2) = aurostd::string2utype<double>(tokens1.at(4));
//[CO20200404 - OBSOLETE]        KlatticeTmp(1,3) = aurostd::string2utype<double>(tokens1.at(5));
//[CO20200404 - OBSOLETE]        //line1=vcontent.at(ii+3);    //CO
//[CO20200404 - OBSOLETE]        line1=vcontent.at(ii+2);
//[CO20200404 - OBSOLETE]        //aurostd::string2tokens(line1,tokens1);  //CO20170725
//[CO20200404 - OBSOLETE]        tokens1=GetCorrectPositions(line1,6); //CO20170725
//[CO20200404 - OBSOLETE]        if(!tokens1.size()){
//[CO20200404 - OBSOLETE]          ERROR = "line with lattice vector is ill-written, see: "+line1;
//[CO20200404 - OBSOLETE]          return false;
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]        //direct_lattice(2,1) = aurostd::string2utype<double>(tokens1.at(3).substr(0,tokens1.at(3).size()-1));    //CO
//[CO20200404 - OBSOLETE]        //direct_lattice(2,2) = aurostd::string2utype<double>(tokens1.at(4).substr(0,tokens1.at(4).size()-1));    //CO
//[CO20200404 - OBSOLETE]        //direct_lattice(2,3) = aurostd::string2utype<double>(tokens1.at(5).substr(0,tokens1.at(5).size()-1));    //CO
//[CO20200404 - OBSOLETE]        direct_lattice(2,1) = aurostd::string2utype<double>(tokens1.at(0));
//[CO20200404 - OBSOLETE]        direct_lattice(2,2) = aurostd::string2utype<double>(tokens1.at(1));
//[CO20200404 - OBSOLETE]        direct_lattice(2,3) = aurostd::string2utype<double>(tokens1.at(2));
//[CO20200404 - OBSOLETE]        //CO just grab reciprocal lattice vectors as calculated by vasp, faster than calculating yourself
//[CO20200404 - OBSOLETE]        KlatticeTmp(2,1) = aurostd::string2utype<double>(tokens1.at(3));
//[CO20200404 - OBSOLETE]        KlatticeTmp(2,2) = aurostd::string2utype<double>(tokens1.at(4));
//[CO20200404 - OBSOLETE]        KlatticeTmp(2,3) = aurostd::string2utype<double>(tokens1.at(5));
//[CO20200404 - OBSOLETE]        //line1=vcontent.at(ii+4);    //CO
//[CO20200404 - OBSOLETE]        line1=vcontent.at(ii+3);
//[CO20200404 - OBSOLETE]        //aurostd::string2tokens(line1,tokens1);  //CO20170725
//[CO20200404 - OBSOLETE]        tokens1=GetCorrectPositions(line1,6); //CO20170725
//[CO20200404 - OBSOLETE]        if(!tokens1.size()){
//[CO20200404 - OBSOLETE]          ERROR = "line with lattice vector is ill-written, see: "+line1;
//[CO20200404 - OBSOLETE]          return false;
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]        //direct_lattice(3,1) = aurostd::string2utype<double>(tokens1.at(3).substr(0,tokens1.at(3).size()-1));    //CO
//[CO20200404 - OBSOLETE]        //direct_lattice(3,2) = aurostd::string2utype<double>(tokens1.at(4).substr(0,tokens1.at(4).size()-1));    //CO
//[CO20200404 - OBSOLETE]        //direct_lattice(3,3) = aurostd::string2utype<double>(tokens1.at(5).substr(0,tokens1.at(5).size()-1));    //CO
//[CO20200404 - OBSOLETE]        direct_lattice(3,1) = aurostd::string2utype<double>(tokens1.at(0));
//[CO20200404 - OBSOLETE]        direct_lattice(3,2) = aurostd::string2utype<double>(tokens1.at(1));
//[CO20200404 - OBSOLETE]        direct_lattice(3,3) = aurostd::string2utype<double>(tokens1.at(2));
//[CO20200404 - OBSOLETE]        //CO just grab reciprocal lattice vectors as calculated by vasp, faster than calculating yourself
//[CO20200404 - OBSOLETE]        KlatticeTmp(3,1) = aurostd::string2utype<double>(tokens1.at(3));
//[CO20200404 - OBSOLETE]        KlatticeTmp(3,2) = aurostd::string2utype<double>(tokens1.at(4));
//[CO20200404 - OBSOLETE]        KlatticeTmp(3,3) = aurostd::string2utype<double>(tokens1.at(5));
//[CO20200404 - OBSOLETE]        //KlatticeTmp = ReciprocalLattice(direct_lattice); // scale is 1  //CO already fetched
//[CO20200404 - OBSOLETE]
//[CO20200404 - OBSOLETE]        //CO20171002
//[CO20200404 - OBSOLETE]        metric_tensor=MetricTensor(direct_lattice,1.0);
//[CO20200404 - OBSOLETE]
//[CO20200404 - OBSOLETE]        if(c_dlattice) { //CO increment tagcount once
//[CO20200404 - OBSOLETE]          c_dlattice = FALSE;
//[CO20200404 - OBSOLETE]          tagcount++;
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]    }
//[CO20200404 - OBSOLETE]    if(tagcount > 7) break;
//[CO20200404 - OBSOLETE]  }
//[CO20200404 - OBSOLETE]  if(tagcount < 7) {
//[CO20200404 - OBSOLETE]    ERROR = "OUTCAR file lacks information needed for Egap determination";
//[CO20200404 - OBSOLETE]    return FALSE;
//[CO20200404 - OBSOLETE]  }
//[CO20200404 - OBSOLETE]  //-------------------------------------------------------------------------------------------------
//[CO20200404 - OBSOLETE]  //   BEGIN Egap DETERMINATIONS
//[CO20200404 - OBSOLETE]  //-------------------------------------------------------------------------------------------------
//[CO20200404 - OBSOLETE]  xmatrix<double> cbb_kpoints(ISPIN,3);
//[CO20200404 - OBSOLETE]  xmatrix<double> vbt_kpoints(ISPIN,3);
//[CO20200404 - OBSOLETE]  conduction_band_min.resize(ISPIN);
//[CO20200404 - OBSOLETE]  valence_band_max.resize(ISPIN);
//[CO20200404 - OBSOLETE]  Egap_type.resize(ISPIN);
//[CO20200404 - OBSOLETE]  Egap_fit.resize(ISPIN);
//[CO20200404 - OBSOLETE]  Egap.resize(ISPIN);
//[CO20200404 - OBSOLETE]  VBT.resize(ISPIN);
//[CO20200404 - OBSOLETE]  CBB.resize(ISPIN);
//[CO20200404 - OBSOLETE]  for(int ii=0; ii<(int)ISPIN; ii++) {
//[CO20200404 - OBSOLETE]    VBT.at(ii) = EMIN;
//[CO20200404 - OBSOLETE]    CBB.at(ii) = EMAX;
//[CO20200404 - OBSOLETE]  }
//[CO20200404 - OBSOLETE]  //-------------------------------------------------------------------------------------------------
//[CO20200404 - OBSOLETE]  // SPIN UNPOLARIZED
//[CO20200404 - OBSOLETE]  if(ISPIN == 1) {
//[CO20200404 - OBSOLETE]    xvector<double> kptdist(ISPIN); //CO20171002 - changing this to be in ANGSTROMS!
//[CO20200404 - OBSOLETE]    vector<uint> strtkpt, band_ndx;
//[CO20200404 - OBSOLETE]    for(uint ii=0; ii<vcontent.size(); ii++) {
//[CO20200404 - OBSOLETE]      if(aurostd::substring2bool(vcontent.at(ii),"k-point   1")) {
//[CO20200404 - OBSOLETE]        strtkpt.push_back(ii);
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]    }
//[CO20200404 - OBSOLETE]    if(strtkpt.size() == 0) {
//[CO20200404 - OBSOLETE]      ERROR = "Corrupt OUTCAR: 'k-point   1' not found.";
//[CO20200404 - OBSOLETE]      return FALSE;
//[CO20200404 - OBSOLETE]    }
//[CO20200404 - OBSOLETE]    for(int ii=1; ii<=(int)ISPIN; ii++) kptdist(ii) = 1.0E09;
//[CO20200404 - OBSOLETE]    for(uint ii=strtkpt.at(strtkpt.size()-1)+3; ii<vcontent.size(); ii++) {
//[CO20200404 - OBSOLETE]      line0 = vcontent.at(ii);
//[CO20200404 - OBSOLETE]      for(uint jj=0; jj<NKPT; jj++) {
//[CO20200404 - OBSOLETE]        line3 = vcontent.at(ii-3);
//[CO20200404 - OBSOLETE]        for(uint kk=1; kk<NBANDS; kk++) {
//[CO20200404 - OBSOLETE]          line1 = vcontent.at(ii-1);
//[CO20200404 - OBSOLETE]          line2 = vcontent.at(ii);
//[CO20200404 - OBSOLETE]          aurostd::string2tokens(line1,tokens1);
//[CO20200404 - OBSOLETE]          aurostd::string2tokens(line2,tokens2);
//[CO20200404 - OBSOLETE]          double ene_1 = aurostd::string2utype<double>(tokens1.at(1));
//[CO20200404 - OBSOLETE]          double ene_2 = aurostd::string2utype<double>(tokens2.at(1));
//[CO20200404 - OBSOLETE]          double occ_1 = aurostd::string2utype<double>(tokens1.at(2));
//[CO20200404 - OBSOLETE]          double occ_2 = aurostd::string2utype<double>(tokens2.at(2));
//[CO20200404 - OBSOLETE]          // SHARP EDGE
//[CO20200404 - OBSOLETE]          if((abs(occ_1-2.00)<_ELEC_EPS) and (abs(occ_2)<_ELEC_EPS) and (Efermi>ene_1) and (Efermi<ene_2)) {
//[CO20200404 - OBSOLETE]            band_ndx.push_back(kk);
//[CO20200404 - OBSOLETE]            SPIN_UP = TRUE;
//[CO20200404 - OBSOLETE]            if(ene_1 > VBT.at(0)) {
//[CO20200404 - OBSOLETE]              VBT.at(0) = ene_1;
//[CO20200404 - OBSOLETE]              aurostd::string2tokens(line3,tokens3);
//[CO20200404 - OBSOLETE]              vbt_kpoints(1,1) = aurostd::string2utype<double>(tokens3.at(3));
//[CO20200404 - OBSOLETE]              vbt_kpoints(1,2) = aurostd::string2utype<double>(tokens3.at(4));
//[CO20200404 - OBSOLETE]              vbt_kpoints(1,3) = aurostd::string2utype<double>(tokens3.at(5));
//[CO20200404 - OBSOLETE]            }
//[CO20200404 - OBSOLETE]            if(ene_2 < CBB.at(0)) {
//[CO20200404 - OBSOLETE]              CBB.at(0) = ene_2;
//[CO20200404 - OBSOLETE]              aurostd::string2tokens(line3,tokens3);
//[CO20200404 - OBSOLETE]              cbb_kpoints(1,1) = aurostd::string2utype<double>(tokens3.at(3));
//[CO20200404 - OBSOLETE]              cbb_kpoints(1,2) = aurostd::string2utype<double>(tokens3.at(4));
//[CO20200404 - OBSOLETE]              cbb_kpoints(1,3) = aurostd::string2utype<double>(tokens3.at(5));
//[CO20200404 - OBSOLETE]            }
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          // SOFT  EDGE
//[CO20200404 - OBSOLETE]          else if((occ_1>_ELEC_EPS) and (occ_2<_ELEC_EPS) and (Efermi>ene_1) and (Efermi<ene_2)) {
//[CO20200404 - OBSOLETE]            band_ndx.push_back(kk);
//[CO20200404 - OBSOLETE]            SPIN_UP = TRUE;
//[CO20200404 - OBSOLETE]            if(ene_1 > VBT.at(0)) {
//[CO20200404 - OBSOLETE]              VBT.at(0) = ene_1;
//[CO20200404 - OBSOLETE]              aurostd::string2tokens(line3,tokens3);
//[CO20200404 - OBSOLETE]              vbt_kpoints(1,1) = aurostd::string2utype<double>(tokens3.at(3));
//[CO20200404 - OBSOLETE]              vbt_kpoints(1,2) = aurostd::string2utype<double>(tokens3.at(4));
//[CO20200404 - OBSOLETE]              vbt_kpoints(1,3) = aurostd::string2utype<double>(tokens3.at(5));
//[CO20200404 - OBSOLETE]            }
//[CO20200404 - OBSOLETE]            if(ene_2 < CBB.at(0)) {
//[CO20200404 - OBSOLETE]              CBB.at(0) = ene_2;
//[CO20200404 - OBSOLETE]              aurostd::string2tokens(line3,tokens3);
//[CO20200404 - OBSOLETE]              cbb_kpoints(1,1) = aurostd::string2utype<double>(tokens3.at(3));
//[CO20200404 - OBSOLETE]              cbb_kpoints(1,2) = aurostd::string2utype<double>(tokens3.at(4));
//[CO20200404 - OBSOLETE]              cbb_kpoints(1,3) = aurostd::string2utype<double>(tokens3.at(5));
//[CO20200404 - OBSOLETE]            }
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          ii++;
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]        if(!SPIN_UP) break;
//[CO20200404 - OBSOLETE]        ii+=3;
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]      break;
//[CO20200404 - OBSOLETE]    }
//[CO20200404 - OBSOLETE]    // Metallic state via band index check
//[CO20200404 - OBSOLETE]    for(uint ii=1; ii<band_ndx.size(); ii++) {
//[CO20200404 - OBSOLETE]      if(band_ndx.at(ii-1) != band_ndx.at(ii)) {
//[CO20200404 - OBSOLETE]        SPIN_UP = FALSE;
//[CO20200404 - OBSOLETE]        break;
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]    }
//[CO20200404 - OBSOLETE]    // reciprocal space distance between 2 kpoints - units are (A^-1)
//[CO20200404 - OBSOLETE]    xvector<double> vdist_kcart,vdist_kcart_min;
//[CO20200404 - OBSOLETE]    double dist,dist_min=AUROSTD_MAX_DOUBLE;
//[CO20200404 - OBSOLETE]    for(int ii=-2; ii<=2; ii++) {
//[CO20200404 - OBSOLETE]      for(int jj=-2; jj<=2; jj++) {
//[CO20200404 - OBSOLETE]        for(int kk=-2; kk<=2; kk++) {
//[CO20200404 - OBSOLETE]          vdist_kcart=(vbt_kpoints(1)-cbb_kpoints(1) + (double)ii*KlatticeTmp(1)+(double)jj*KlatticeTmp(2)+(double)kk*KlatticeTmp(3));
//[CO20200404 - OBSOLETE]          dist=aurostd::modulus(vdist_kcart);
//[CO20200404 - OBSOLETE]          if(dist<dist_min){
//[CO20200404 - OBSOLETE]            vdist_kcart_min=vdist_kcart;
//[CO20200404 - OBSOLETE]            dist_min=dist;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          //CO20171002
//[CO20200404 - OBSOLETE]          // kptdist is basically an energy 
//[CO20200404 - OBSOLETE]          //kptdist(1) = min(kptdist(1),modulus(vbt_kpoints(1)-cbb_kpoints(1) + 
//[CO20200404 - OBSOLETE]          //			      (double)ii*KlatticeTmp(1)+(double)jj*KlatticeTmp(2)+(double)kk*KlatticeTmp(3)));
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]    }
//[CO20200404 - OBSOLETE]    //CO20171002
//[CO20200404 - OBSOLETE]    xvector<double> vdist_dcart;
//[CO20200404 - OBSOLETE]    vdist_dcart(1)=sum(metric_tensor(1)*vdist_kcart_min(1));  //remember metric_tensor is symmetric!
//[CO20200404 - OBSOLETE]    vdist_dcart(2)=sum(metric_tensor(2)*vdist_kcart_min(2));  //remember metric_tensor is symmetric!
//[CO20200404 - OBSOLETE]    vdist_dcart(3)=sum(metric_tensor(3)*vdist_kcart_min(3));  //remember metric_tensor is symmetric!
//[CO20200404 - OBSOLETE]    kptdist(1)=aurostd::modulus(vdist_dcart);
//[CO20200404 - OBSOLETE]    // OUTPUT
//[CO20200404 - OBSOLETE]    if(SPIN_UP) {
//[CO20200404 - OBSOLETE]      conduction_band_min.at(0) = CBB.at(0) - Efermi;
//[CO20200404 - OBSOLETE]      valence_band_max.at(0) = VBT.at(0) - Efermi;
//[CO20200404 - OBSOLETE]      conduction_band_min_net = conduction_band_min.at(0);
//[CO20200404 - OBSOLETE]      valence_band_max_net = valence_band_max.at(0);
//[CO20200404 - OBSOLETE]      if(CBB.at(0)-VBT.at(0) > _ENER_EPS) { // u1
//[CO20200404 - OBSOLETE]        //[CO20171002 - GARBAGE A != A^-1]if      (kptdist(1) <= CELLVOLCUTOFF) Egap_type.at(0) = "insulator_direct";
//[CO20200404 - OBSOLETE]        //[CO20171002 - GARBAGE A != A^-1]else if(kptdist(1) >  CELLVOLCUTOFF) Egap_type.at(0) = "insulator_indirect";
//[CO20200404 - OBSOLETE]        if      (kptdist(1) <= kpt_tol) Egap_type.at(0) = "insulator_direct";
//[CO20200404 - OBSOLETE]        else if(kptdist(1) >  kpt_tol) Egap_type.at(0) = "insulator_indirect";
//[CO20200404 - OBSOLETE]        Egap_type_net = Egap_type.at(0);
//[CO20200404 - OBSOLETE]        Egap.at(0)    = CBB.at(0) - VBT.at(0);
//[CO20200404 - OBSOLETE]        Egap_net      = Egap.at(0);
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]      else if((abs(CBB.at(0)-VBT.at(0)) <= _ENER_EPS) and (abs(CBB.at(0)-VBT.at(0)) >= 0)) { // u2
//[CO20200404 - OBSOLETE]        //[CO20171002 - GARBAGE A != A^-1]if      (kptdist(1) <= CELLVOLCUTOFF) Egap_type.at(0) = "insulator_direct";
//[CO20200404 - OBSOLETE]        //[CO20171002 - GARBAGE A != A^-1]else if(kptdist(1) >  CELLVOLCUTOFF) Egap_type.at(0) = "insulator_indirect";
//[CO20200404 - OBSOLETE]        if      (kptdist(1) <= kpt_tol) Egap_type.at(0) = "insulator_direct";
//[CO20200404 - OBSOLETE]        else if(kptdist(1) >  kpt_tol) Egap_type.at(0) = "insulator_indirect";
//[CO20200404 - OBSOLETE]        Egap_type_net = Egap_type.at(0);
//[CO20200404 - OBSOLETE]        Egap.at(0)    = _ZERO;
//[CO20200404 - OBSOLETE]        Egap_net      = Egap.at(0);
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]      //WSETYAWAN
//[CO20200404 - OBSOLETE]      Egap_fit.at(0) = 1.348 * Egap.at(0) + 0.913;
//[CO20200404 - OBSOLETE]      Egap_fit_net   = Egap_fit.at(0);
//[CO20200404 - OBSOLETE]    }
//[CO20200404 - OBSOLETE]    else if(!SPIN_UP) { // u3
//[CO20200404 - OBSOLETE]      conduction_band_min.at(0) = _METALEDGE_;
//[CO20200404 - OBSOLETE]      valence_band_max.at(0)    = _METALEDGE_;
//[CO20200404 - OBSOLETE]      conduction_band_min_net = conduction_band_min.at(0);
//[CO20200404 - OBSOLETE]      valence_band_max_net    = valence_band_max.at(0);
//[CO20200404 - OBSOLETE]      Egap_type.at(0) = "metal";
//[CO20200404 - OBSOLETE]      Egap_type_net   = Egap_type.at(0);
//[CO20200404 - OBSOLETE]      Egap.at(0)      = _METALGAP_;
//[CO20200404 - OBSOLETE]      Egap_net        = Egap.at(0);
//[CO20200404 - OBSOLETE]      //WSETYAWAN
//[CO20200404 - OBSOLETE]      Egap_fit.at(0) = Egap.at(0);
//[CO20200404 - OBSOLETE]      Egap_fit_net   = Egap_fit.at(0);
//[CO20200404 - OBSOLETE]    }
//[CO20200404 - OBSOLETE]  }
//[CO20200404 - OBSOLETE]  //-------------------------------------------------------------------------------------------------
//[CO20200404 - OBSOLETE]  // SPIN POLARIZED
//[CO20200404 - OBSOLETE]  else if(ISPIN == 2) {
//[CO20200404 - OBSOLETE]    bool SPIN_UP_ALRT = FALSE;
//[CO20200404 - OBSOLETE]    bool SPIN_DN_ALRT = FALSE;
//[CO20200404 - OBSOLETE]    vector<uint> band_ndx_up, band_ndx_dn;
//[CO20200404 - OBSOLETE]    xvector<double> kptdist(2*ISPIN); //CO20171002 - changing this to be in ANGSTROMS!
//[CO20200404 - OBSOLETE]    for(int ii=1; ii<=2*(int)ISPIN; ii++) { kptdist(ii) = 1.0E09; }
//[CO20200404 - OBSOLETE]    // SPIN UP
//[CO20200404 - OBSOLETE]    for(uint ii=0; ii<vcontent.size(); ii++) {
//[CO20200404 - OBSOLETE]      if(aurostd::substring2bool(vcontent.at(ii),"spin component 1")) {
//[CO20200404 - OBSOLETE]        if(aurostd::substring2bool(vcontent.at(ii+1),"k-point   1 :")) {
//[CO20200404 - OBSOLETE]          ii+=4;
//[CO20200404 - OBSOLETE]          for(uint jj=1; jj<=NKPT; jj++) {
//[CO20200404 - OBSOLETE]            line3 = vcontent.at(ii-3);
//[CO20200404 - OBSOLETE]            for(uint kk=1; kk<NBANDS; kk++) {
//[CO20200404 - OBSOLETE]              line1 = vcontent.at(ii-1); // band #1
//[CO20200404 - OBSOLETE]              line2 = vcontent.at(ii);   // band #2
//[CO20200404 - OBSOLETE]              aurostd::string2tokens(line1,tokens1);
//[CO20200404 - OBSOLETE]              aurostd::string2tokens(line2,tokens2);
//[CO20200404 - OBSOLETE]              double ene_1 = aurostd::string2utype<double>(tokens1.at(1));
//[CO20200404 - OBSOLETE]              double ene_2 = aurostd::string2utype<double>(tokens2.at(1));
//[CO20200404 - OBSOLETE]              double occ_1 = aurostd::string2utype<double>(tokens1.at(2));
//[CO20200404 - OBSOLETE]              double occ_2 = aurostd::string2utype<double>(tokens2.at(2));
//[CO20200404 - OBSOLETE]              // SHARP EDGE
//[CO20200404 - OBSOLETE]              if((abs(occ_1-1.00)<_ELEC_EPS) and (abs(occ_2)<_ELEC_EPS) and (Efermi>ene_1) and (Efermi<ene_2)) {
//[CO20200404 - OBSOLETE]                band_ndx_up.push_back(kk);
//[CO20200404 - OBSOLETE]                SPIN_UP = TRUE;
//[CO20200404 - OBSOLETE]                if(ene_1 > VBT.at(0)) {
//[CO20200404 - OBSOLETE]                  VBT.at(0) = ene_1;
//[CO20200404 - OBSOLETE]                  aurostd::string2tokens(line3,tokens3);
//[CO20200404 - OBSOLETE]                  vbt_kpoints(1,1) = aurostd::string2utype<double>(tokens3.at(3));
//[CO20200404 - OBSOLETE]                  vbt_kpoints(1,2) = aurostd::string2utype<double>(tokens3.at(4));
//[CO20200404 - OBSOLETE]                  vbt_kpoints(1,3) = aurostd::string2utype<double>(tokens3.at(5));
//[CO20200404 - OBSOLETE]                }
//[CO20200404 - OBSOLETE]                if(ene_2 < CBB.at(0)) {
//[CO20200404 - OBSOLETE]                  CBB.at(0) = ene_2;
//[CO20200404 - OBSOLETE]                  aurostd::string2tokens(line3,tokens3);
//[CO20200404 - OBSOLETE]                  cbb_kpoints(1,1) = aurostd::string2utype<double>(tokens3.at(3));
//[CO20200404 - OBSOLETE]                  cbb_kpoints(1,2) = aurostd::string2utype<double>(tokens3.at(4));
//[CO20200404 - OBSOLETE]                  cbb_kpoints(1,3) = aurostd::string2utype<double>(tokens3.at(5));
//[CO20200404 - OBSOLETE]                }
//[CO20200404 - OBSOLETE]              }
//[CO20200404 - OBSOLETE]              // SOFT  EDGE
//[CO20200404 - OBSOLETE]              else if((occ_1>_ELEC_EPS) and (occ_2<_ELEC_EPS) and (Efermi>ene_1) and (Efermi<ene_2)) {
//[CO20200404 - OBSOLETE]                band_ndx_up.push_back(kk);
//[CO20200404 - OBSOLETE]                SPIN_UP = TRUE;
//[CO20200404 - OBSOLETE]                if(ene_1 > VBT.at(0)) {
//[CO20200404 - OBSOLETE]                  VBT.at(0) = ene_1;
//[CO20200404 - OBSOLETE]                  aurostd::string2tokens(line3,tokens3);
//[CO20200404 - OBSOLETE]                  vbt_kpoints(1,1) = aurostd::string2utype<double>(tokens3.at(3));
//[CO20200404 - OBSOLETE]                  vbt_kpoints(1,2) = aurostd::string2utype<double>(tokens3.at(4));
//[CO20200404 - OBSOLETE]                  vbt_kpoints(1,3) = aurostd::string2utype<double>(tokens3.at(5));
//[CO20200404 - OBSOLETE]                }
//[CO20200404 - OBSOLETE]                if(ene_2 < CBB.at(0)) {
//[CO20200404 - OBSOLETE]                  CBB.at(0) = ene_2;
//[CO20200404 - OBSOLETE]                  aurostd::string2tokens(line3,tokens3);
//[CO20200404 - OBSOLETE]                  cbb_kpoints(1,1) = aurostd::string2utype<double>(tokens3.at(3));
//[CO20200404 - OBSOLETE]                  cbb_kpoints(1,2) = aurostd::string2utype<double>(tokens3.at(4));
//[CO20200404 - OBSOLETE]                  cbb_kpoints(1,3) = aurostd::string2utype<double>(tokens3.at(5));
//[CO20200404 - OBSOLETE]                }
//[CO20200404 - OBSOLETE]              }
//[CO20200404 - OBSOLETE]              // special case: empty spin channel
//[CO20200404 - OBSOLETE]              else if((kk == 1) and (ene_1-Efermi > _ENER_EPS)) {
//[CO20200404 - OBSOLETE]                SPIN_UP_ALRT = TRUE;
//[CO20200404 - OBSOLETE]                SPIN_UP      = FALSE;
//[CO20200404 - OBSOLETE]                break;
//[CO20200404 - OBSOLETE]              }
//[CO20200404 - OBSOLETE]              ii++;
//[CO20200404 - OBSOLETE]            }
//[CO20200404 - OBSOLETE]            if(!SPIN_UP) break;
//[CO20200404 - OBSOLETE]            ii+=3;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          break;
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]    }
//[CO20200404 - OBSOLETE]    // SPIN DN
//[CO20200404 - OBSOLETE]    for(uint ii=0; ii<vcontent.size(); ii++) {
//[CO20200404 - OBSOLETE]      if(aurostd::substring2bool(vcontent.at(ii),"spin component 2")) {
//[CO20200404 - OBSOLETE]        if(aurostd::substring2bool(vcontent.at(ii+1),"k-point   1 :")) {
//[CO20200404 - OBSOLETE]          ii+=4; // starts at band #2
//[CO20200404 - OBSOLETE]          for(uint jj=1; jj<=NKPT; jj++) {
//[CO20200404 - OBSOLETE]            line3 = vcontent.at(ii-3);
//[CO20200404 - OBSOLETE]            for(uint kk=1; kk<NBANDS; kk++) {
//[CO20200404 - OBSOLETE]              line1 = vcontent.at(ii-1); // band #1
//[CO20200404 - OBSOLETE]              line2 = vcontent.at(ii);   // band #2
//[CO20200404 - OBSOLETE]              aurostd::string2tokens(line1,tokens1);
//[CO20200404 - OBSOLETE]              aurostd::string2tokens(line2,tokens2);
//[CO20200404 - OBSOLETE]              double ene_1 = aurostd::string2utype<double>(tokens1.at(1));
//[CO20200404 - OBSOLETE]              double ene_2 = aurostd::string2utype<double>(tokens2.at(1));
//[CO20200404 - OBSOLETE]              double occ_1 = aurostd::string2utype<double>(tokens1.at(2));
//[CO20200404 - OBSOLETE]              double occ_2 = aurostd::string2utype<double>(tokens2.at(2));
//[CO20200404 - OBSOLETE]              // SHARP EDGE
//[CO20200404 - OBSOLETE]              if((abs(occ_1-1.00)<_ELEC_EPS) and (abs(occ_2)<_ELEC_EPS) and (Efermi>ene_1) and (Efermi<ene_2)) {
//[CO20200404 - OBSOLETE]                band_ndx_dn.push_back(kk);
//[CO20200404 - OBSOLETE]                SPIN_DN = TRUE;
//[CO20200404 - OBSOLETE]                if(ene_1 > VBT.at(1)) {
//[CO20200404 - OBSOLETE]                  VBT.at(1) = ene_1;
//[CO20200404 - OBSOLETE]                  aurostd::string2tokens(line3,tokens3);
//[CO20200404 - OBSOLETE]                  vbt_kpoints(1,1) = aurostd::string2utype<double>(tokens3.at(3));
//[CO20200404 - OBSOLETE]                  vbt_kpoints(1,2) = aurostd::string2utype<double>(tokens3.at(4));
//[CO20200404 - OBSOLETE]                  vbt_kpoints(1,3) = aurostd::string2utype<double>(tokens3.at(5));
//[CO20200404 - OBSOLETE]                }
//[CO20200404 - OBSOLETE]                if(ene_2 < CBB.at(1)) {
//[CO20200404 - OBSOLETE]                  CBB.at(1) = ene_2;
//[CO20200404 - OBSOLETE]                  aurostd::string2tokens(line3,tokens3);
//[CO20200404 - OBSOLETE]                  cbb_kpoints(1,1) = aurostd::string2utype<double>(tokens3.at(3));
//[CO20200404 - OBSOLETE]                  cbb_kpoints(1,2) = aurostd::string2utype<double>(tokens3.at(4));
//[CO20200404 - OBSOLETE]                  cbb_kpoints(1,3) = aurostd::string2utype<double>(tokens3.at(5));
//[CO20200404 - OBSOLETE]                }
//[CO20200404 - OBSOLETE]              }
//[CO20200404 - OBSOLETE]              // SOFT  EDGE
//[CO20200404 - OBSOLETE]              else if((occ_1>_ELEC_EPS) and (occ_2<_ELEC_EPS) and (Efermi>ene_1) and (Efermi<ene_2)) {
//[CO20200404 - OBSOLETE]                band_ndx_dn.push_back(kk);
//[CO20200404 - OBSOLETE]                SPIN_DN = TRUE;
//[CO20200404 - OBSOLETE]                if(ene_1 > VBT.at(1)) {
//[CO20200404 - OBSOLETE]                  VBT.at(1) = ene_1;
//[CO20200404 - OBSOLETE]                  aurostd::string2tokens(line3,tokens3);
//[CO20200404 - OBSOLETE]                  vbt_kpoints(1,1) = aurostd::string2utype<double>(tokens3.at(3));
//[CO20200404 - OBSOLETE]                  vbt_kpoints(1,2) = aurostd::string2utype<double>(tokens3.at(4));
//[CO20200404 - OBSOLETE]                  vbt_kpoints(1,3) = aurostd::string2utype<double>(tokens3.at(5));
//[CO20200404 - OBSOLETE]                }
//[CO20200404 - OBSOLETE]                if(ene_2 < CBB.at(1)) {
//[CO20200404 - OBSOLETE]                  CBB.at(1) = ene_2;
//[CO20200404 - OBSOLETE]                  aurostd::string2tokens(line3,tokens3);
//[CO20200404 - OBSOLETE]                  cbb_kpoints(1,1) = aurostd::string2utype<double>(tokens3.at(3));
//[CO20200404 - OBSOLETE]                  cbb_kpoints(1,2) = aurostd::string2utype<double>(tokens3.at(4));
//[CO20200404 - OBSOLETE]                  cbb_kpoints(1,3) = aurostd::string2utype<double>(tokens3.at(5));
//[CO20200404 - OBSOLETE]                }
//[CO20200404 - OBSOLETE]              }
//[CO20200404 - OBSOLETE]              // special case: empty spin channel
//[CO20200404 - OBSOLETE]              else if((kk == 1) and (ene_1-Efermi > _ENER_EPS)) {
//[CO20200404 - OBSOLETE]                SPIN_DN_ALRT = TRUE;
//[CO20200404 - OBSOLETE]                SPIN_DN = FALSE;
//[CO20200404 - OBSOLETE]                break;
//[CO20200404 - OBSOLETE]              }
//[CO20200404 - OBSOLETE]              ii++;
//[CO20200404 - OBSOLETE]            }
//[CO20200404 - OBSOLETE]            if(!SPIN_DN) break;
//[CO20200404 - OBSOLETE]            ii+=3;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          break;
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]    }
//[CO20200404 - OBSOLETE]    // Metallic state via band index check
//[CO20200404 - OBSOLETE]    for(uint ii=1; ii<band_ndx_up.size(); ii++) {
//[CO20200404 - OBSOLETE]      if(band_ndx_up.at(ii-1) != band_ndx_up.at(ii)) {
//[CO20200404 - OBSOLETE]        SPIN_UP = FALSE;
//[CO20200404 - OBSOLETE]        break;
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]    }
//[CO20200404 - OBSOLETE]    for(uint ii=1; ii<band_ndx_dn.size(); ii++) {
//[CO20200404 - OBSOLETE]      if(band_ndx_dn.at(ii-1) != band_ndx_dn.at(ii)) {
//[CO20200404 - OBSOLETE]        SPIN_DN = FALSE;
//[CO20200404 - OBSOLETE]        break;
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]    }
//[CO20200404 - OBSOLETE]    //CO20171002
//[CO20200404 - OBSOLETE]    xvector<double> vdist1_kcart,vdist1_kcart_min,vdist2_kcart,vdist2_kcart_min,vdist3_kcart,vdist3_kcart_min,vdist4_kcart,vdist4_kcart_min;
//[CO20200404 - OBSOLETE]    double dist1,dist2,dist3,dist4,dist1_min,dist2_min,dist3_min,dist4_min;
//[CO20200404 - OBSOLETE]    dist1_min=dist2_min=dist3_min=dist4_min=AUROSTD_MAX_DOUBLE;
//[CO20200404 - OBSOLETE]    // reciprocal space distance between 2 kpoints - units are (A^-1)
//[CO20200404 - OBSOLETE]    for(int ii=-2; ii<=2; ii++) {
//[CO20200404 - OBSOLETE]      for(int jj=-2; jj<=2; jj++) {
//[CO20200404 - OBSOLETE]        for(int kk=-2; kk<=2; kk++) {
//[CO20200404 - OBSOLETE]          // within spin channels
//[CO20200404 - OBSOLETE]          vdist1_kcart=(vbt_kpoints(1)-cbb_kpoints(1) + (double)ii*KlatticeTmp(1)+(double)jj*KlatticeTmp(2)+(double)kk*KlatticeTmp(3));
//[CO20200404 - OBSOLETE]          dist1=aurostd::modulus(vdist1_kcart);
//[CO20200404 - OBSOLETE]          if(dist1<dist1_min){
//[CO20200404 - OBSOLETE]            vdist1_kcart_min=vdist1_kcart;
//[CO20200404 - OBSOLETE]            dist1_min=dist1;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          vdist2_kcart=(vbt_kpoints(1)-cbb_kpoints(1) + (double)ii*KlatticeTmp(1)+(double)jj*KlatticeTmp(2)+(double)kk*KlatticeTmp(3));
//[CO20200404 - OBSOLETE]          dist2=aurostd::modulus(vdist2_kcart);
//[CO20200404 - OBSOLETE]          if(dist2<dist2_min){
//[CO20200404 - OBSOLETE]            vdist2_kcart_min=vdist2_kcart;
//[CO20200404 - OBSOLETE]            dist2_min=dist2;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          // across spin channels
//[CO20200404 - OBSOLETE]          vdist3_kcart=(vbt_kpoints(1)-cbb_kpoints(2) + (double)ii*KlatticeTmp(1)+(double)jj*KlatticeTmp(2)+(double)kk*KlatticeTmp(3));
//[CO20200404 - OBSOLETE]          dist3=aurostd::modulus(vdist3_kcart);
//[CO20200404 - OBSOLETE]          if(dist3<dist3_min){
//[CO20200404 - OBSOLETE]            vdist3_kcart_min=vdist3_kcart;
//[CO20200404 - OBSOLETE]            dist3_min=dist3;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          vdist4_kcart=(vbt_kpoints(2)-cbb_kpoints(1) + (double)ii*KlatticeTmp(1)+(double)jj*KlatticeTmp(2)+(double)kk*KlatticeTmp(3));
//[CO20200404 - OBSOLETE]          dist4=aurostd::modulus(vdist4_kcart);
//[CO20200404 - OBSOLETE]          if(dist4<dist4_min){
//[CO20200404 - OBSOLETE]            vdist4_kcart_min=vdist4_kcart;
//[CO20200404 - OBSOLETE]            dist4_min=dist4;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          //CO20171002
//[CO20200404 - OBSOLETE]          //// within spin channels
//[CO20200404 - OBSOLETE]          //kptdist(1) = min(kptdist(1),modulus(vbt_kpoints(1)-cbb_kpoints(1) + 
//[CO20200404 - OBSOLETE]          //			      (double)ii*KlatticeTmp(1)+(double)jj*KlatticeTmp(2)+(double)kk*KlatticeTmp(3)));
//[CO20200404 - OBSOLETE]          //kptdist(2) = min(kptdist(2),modulus(vbt_kpoints(2)-cbb_kpoints(2) + 
//[CO20200404 - OBSOLETE]          //			      (double)ii*KlatticeTmp(1)+(double)jj*KlatticeTmp(2)+(double)kk*KlatticeTmp(3)));
//[CO20200404 - OBSOLETE]          //// across spin channels
//[CO20200404 - OBSOLETE]          //kptdist(3) = min(kptdist(3),modulus(vbt_kpoints(1)-cbb_kpoints(2) + 
//[CO20200404 - OBSOLETE]          //			      (double)ii*KlatticeTmp(1)+(double)jj*KlatticeTmp(2)+(double)kk*KlatticeTmp(3)));
//[CO20200404 - OBSOLETE]          //kptdist(4) = min(kptdist(4),modulus(vbt_kpoints(2)-cbb_kpoints(1) + 
//[CO20200404 - OBSOLETE]          //			      (double)ii*KlatticeTmp(1)+(double)jj*KlatticeTmp(2)+(double)kk*KlatticeTmp(3)));
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]    }
//[CO20200404 - OBSOLETE]    //CO20171002
//[CO20200404 - OBSOLETE]    xvector<double> vdist1_dcart,vdist2_dcart,vdist3_dcart,vdist4_dcart;
//[CO20200404 - OBSOLETE]    vdist1_dcart(1)=sum(metric_tensor(1)*vdist1_kcart_min(1));  //remember metric_tensor is symmetric!
//[CO20200404 - OBSOLETE]    vdist1_dcart(2)=sum(metric_tensor(2)*vdist1_kcart_min(2));  //remember metric_tensor is symmetric!
//[CO20200404 - OBSOLETE]    vdist1_dcart(3)=sum(metric_tensor(3)*vdist1_kcart_min(3));  //remember metric_tensor is symmetric!
//[CO20200404 - OBSOLETE]    kptdist(1)=modulus(vdist1_dcart);
//[CO20200404 - OBSOLETE]
//[CO20200404 - OBSOLETE]    vdist2_dcart(1)=sum(metric_tensor(1)*vdist2_kcart_min(1));  //remember metric_tensor is symmetric!
//[CO20200404 - OBSOLETE]    vdist2_dcart(2)=sum(metric_tensor(2)*vdist2_kcart_min(2));  //remember metric_tensor is symmetric!
//[CO20200404 - OBSOLETE]    vdist2_dcart(3)=sum(metric_tensor(3)*vdist2_kcart_min(3));  //remember metric_tensor is symmetric!
//[CO20200404 - OBSOLETE]    kptdist(2)=modulus(vdist2_dcart);
//[CO20200404 - OBSOLETE]
//[CO20200404 - OBSOLETE]    vdist3_dcart(1)=sum(metric_tensor(1)*vdist3_kcart_min(1));  //remember metric_tensor is symmetric!
//[CO20200404 - OBSOLETE]    vdist3_dcart(2)=sum(metric_tensor(2)*vdist3_kcart_min(2));  //remember metric_tensor is symmetric!
//[CO20200404 - OBSOLETE]    vdist3_dcart(3)=sum(metric_tensor(3)*vdist3_kcart_min(3));  //remember metric_tensor is symmetric!
//[CO20200404 - OBSOLETE]    kptdist(3)=modulus(vdist3_dcart);
//[CO20200404 - OBSOLETE]
//[CO20200404 - OBSOLETE]    vdist4_dcart(1)=sum(metric_tensor(1)*vdist4_kcart_min(1));  //remember metric_tensor is symmetric!
//[CO20200404 - OBSOLETE]    vdist4_dcart(2)=sum(metric_tensor(2)*vdist4_kcart_min(2));  //remember metric_tensor is symmetric!
//[CO20200404 - OBSOLETE]    vdist4_dcart(3)=sum(metric_tensor(3)*vdist4_kcart_min(3));  //remember metric_tensor is symmetric!
//[CO20200404 - OBSOLETE]    kptdist(4)=modulus(vdist4_dcart);
//[CO20200404 - OBSOLETE]    // ---------------------------------------------------------------------------------------------------
//[CO20200404 - OBSOLETE]    // These are systems that are actually not polarized or only have 1 usable spin channel
//[CO20200404 - OBSOLETE]    if(SPIN_UP_ALRT or SPIN_DN_ALRT) {
//[CO20200404 - OBSOLETE]      conduction_band_min.resize(1);
//[CO20200404 - OBSOLETE]      valence_band_max.resize(1);
//[CO20200404 - OBSOLETE]      Egap_type.resize(1);
//[CO20200404 - OBSOLETE]      Egap.resize(1);
//[CO20200404 - OBSOLETE]      if(SPIN_DN_ALRT) {
//[CO20200404 - OBSOLETE]        if      (kptdist(1) <= kpt_tol) Egap_type.at(0) = "insulator_direct";
//[CO20200404 - OBSOLETE]        else if(kptdist(1) >  kpt_tol) Egap_type.at(0) = "insulator_indirect";
//[CO20200404 - OBSOLETE]        //[CO20171002 - GARBAGE A != A^-1]if      (kptdist(1) <= CELLVOLCUTOFF) Egap_type.at(0) = "insulator_direct";
//[CO20200404 - OBSOLETE]        //[CO20171002 - GARBAGE A != A^-1]else if(kptdist(1) >  CELLVOLCUTOFF) Egap_type.at(0) = "insulator_indirect";
//[CO20200404 - OBSOLETE]        if(CBB.at(0)-VBT.at(0) > _ENER_EPS) {
//[CO20200404 - OBSOLETE]          conduction_band_min.at(0) = CBB.at(0) - Efermi;
//[CO20200404 - OBSOLETE]          valence_band_max.at(0)    = VBT.at(0) - Efermi;
//[CO20200404 - OBSOLETE]          conduction_band_min_net   = conduction_band_min.at(0);
//[CO20200404 - OBSOLETE]          valence_band_max_net      = valence_band_max.at(0);
//[CO20200404 - OBSOLETE]          Egap_type_net = Egap_type.at(0);
//[CO20200404 - OBSOLETE]          Egap.at(0)    = CBB.at(0)-VBT.at(0);
//[CO20200404 - OBSOLETE]          Egap_net      = Egap.at(0);
//[CO20200404 - OBSOLETE]          //WSETYAWAN fit
//[CO20200404 - OBSOLETE]          Egap_fit.at(0) = 1.348 * Egap.at(0) + 0.913;
//[CO20200404 - OBSOLETE]          Egap_fit_net   = Egap_fit.at(0);
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]        else if((abs(CBB.at(0)-VBT.at(0)) <= _ENER_EPS) and (abs(CBB.at(0)-VBT.at(0)) >= 0)) {
//[CO20200404 - OBSOLETE]          conduction_band_min.at(0) = CBB.at(0) - Efermi;
//[CO20200404 - OBSOLETE]          valence_band_max.at(0)    = VBT.at(0) - Efermi;
//[CO20200404 - OBSOLETE]          conduction_band_min_net   = conduction_band_min.at(0);
//[CO20200404 - OBSOLETE]          valence_band_max_net      = valence_band_max.at(0);
//[CO20200404 - OBSOLETE]          Egap_type_net = Egap_type.at(0);
//[CO20200404 - OBSOLETE]          Egap.at(0)    = _ZERO;
//[CO20200404 - OBSOLETE]          Egap_net      = Egap.at(0);
//[CO20200404 - OBSOLETE]          //WSETYAWAN fit
//[CO20200404 - OBSOLETE]          Egap_fit.at(0) = 1.348 * Egap.at(0) + 0.913;
//[CO20200404 - OBSOLETE]          Egap_fit_net   = Egap_fit.at(0);
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]        else {
//[CO20200404 - OBSOLETE]          conduction_band_min.at(0) = _METALEDGE_;
//[CO20200404 - OBSOLETE]          valence_band_max.at(0)    = _METALEDGE_;
//[CO20200404 - OBSOLETE]          conduction_band_min_net   = conduction_band_min.at(0);
//[CO20200404 - OBSOLETE]          valence_band_max_net      = valence_band_max.at(0);
//[CO20200404 - OBSOLETE]          Egap_type.at(0) = "metal";
//[CO20200404 - OBSOLETE]          Egap_type_net = Egap_type.at(0);
//[CO20200404 - OBSOLETE]          Egap.at(0)    = _METALGAP_;
//[CO20200404 - OBSOLETE]          Egap_net      = Egap.at(0);
//[CO20200404 - OBSOLETE]          //WSETYAWAN fit
//[CO20200404 - OBSOLETE]          Egap_fit.at(0) = _METALGAP_;
//[CO20200404 - OBSOLETE]          Egap_fit_net   = Egap_fit.at(0);
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]      else if(SPIN_UP_ALRT) {
//[CO20200404 - OBSOLETE]        if      (kptdist(2) <= kpt_tol) Egap_type.at(0) = "insulator_direct";
//[CO20200404 - OBSOLETE]        else if(kptdist(2) >  kpt_tol) Egap_type.at(0) = "insulator_indirect";
//[CO20200404 - OBSOLETE]        //[CO20171002 - GARBAGE A != A^-1]if      (kptdist(2) <= CELLVOLCUTOFF) Egap_type.at(0) = "insulator_direct";
//[CO20200404 - OBSOLETE]        //[CO20171002 - GARBAGE A != A^-1]else if(kptdist(2) >  CELLVOLCUTOFF) Egap_type.at(0) = "insulator_indirect";
//[CO20200404 - OBSOLETE]        if(CBB.at(1)-VBT.at(1) > _ENER_EPS) {
//[CO20200404 - OBSOLETE]          conduction_band_min.at(0) = CBB.at(1) - Efermi;
//[CO20200404 - OBSOLETE]          valence_band_max.at(0)    = VBT.at(1) - Efermi;
//[CO20200404 - OBSOLETE]          conduction_band_min_net   = conduction_band_min.at(0);
//[CO20200404 - OBSOLETE]          valence_band_max_net      = valence_band_max.at(0);
//[CO20200404 - OBSOLETE]          Egap_type_net = Egap_type.at(0);
//[CO20200404 - OBSOLETE]          Egap.at(0)    = CBB.at(1)-VBT.at(1);
//[CO20200404 - OBSOLETE]          Egap_net      = Egap.at(0);
//[CO20200404 - OBSOLETE]          //WSETYAWAN fit
//[CO20200404 - OBSOLETE]          Egap_fit.at(0) = 1.348 * Egap.at(0) + 0.913;
//[CO20200404 - OBSOLETE]          Egap_fit_net   = Egap_fit.at(0);
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]        else if((abs(CBB.at(1)-VBT.at(1)) <= _ENER_EPS) and (abs(CBB.at(1)-VBT.at(1)) >= 0)) {
//[CO20200404 - OBSOLETE]          conduction_band_min.at(0) = CBB.at(1) - Efermi;
//[CO20200404 - OBSOLETE]          valence_band_max.at(0)    = VBT.at(1) - Efermi;
//[CO20200404 - OBSOLETE]          conduction_band_min_net   = conduction_band_min.at(0);
//[CO20200404 - OBSOLETE]          valence_band_max_net      = valence_band_max.at(0);
//[CO20200404 - OBSOLETE]          Egap_type_net = Egap_type.at(0);
//[CO20200404 - OBSOLETE]          Egap.at(0)    = _ZERO;
//[CO20200404 - OBSOLETE]          Egap_net      = Egap.at(0);
//[CO20200404 - OBSOLETE]          //WSETYAWAN fit
//[CO20200404 - OBSOLETE]          Egap_fit.at(0) = 1.348 * Egap.at(0) + 0.913;
//[CO20200404 - OBSOLETE]          Egap_fit_net   = Egap_fit.at(0);
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]        else {
//[CO20200404 - OBSOLETE]          conduction_band_min.at(0) = _METALEDGE_;
//[CO20200404 - OBSOLETE]          valence_band_max.at(0)    = _METALEDGE_;
//[CO20200404 - OBSOLETE]          conduction_band_min_net   = conduction_band_min.at(0);
//[CO20200404 - OBSOLETE]          valence_band_max_net      = valence_band_max.at(0);
//[CO20200404 - OBSOLETE]          Egap_type.at(0) = "metal";
//[CO20200404 - OBSOLETE]          Egap_type_net = Egap_type.at(0);
//[CO20200404 - OBSOLETE]          Egap.at(0)    = _METALGAP_;
//[CO20200404 - OBSOLETE]          Egap_net      = Egap.at(0);
//[CO20200404 - OBSOLETE]          //WSETYAWAN fit
//[CO20200404 - OBSOLETE]          Egap_fit.at(0) = _METALGAP_;
//[CO20200404 - OBSOLETE]          Egap_fit_net   = Egap_fit.at(0);
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]    }
//[CO20200404 - OBSOLETE]    // ---------------------------------------------------------------------------------------------------
//[CO20200404 - OBSOLETE]    else if(!SPIN_UP_ALRT and !SPIN_DN_ALRT) {
//[CO20200404 - OBSOLETE]      conduction_band_min.at(0) = CBB.at(0) - Efermi;
//[CO20200404 - OBSOLETE]      conduction_band_min.at(1) = CBB.at(1) - Efermi;
//[CO20200404 - OBSOLETE]      valence_band_max.at(0)    = VBT.at(0) - Efermi;
//[CO20200404 - OBSOLETE]      valence_band_max.at(1)    = VBT.at(1) - Efermi;
//[CO20200404 - OBSOLETE]      conduction_band_min_net   = min(conduction_band_min);
//[CO20200404 - OBSOLETE]      valence_band_max_net      = max(valence_band_max);
//[CO20200404 - OBSOLETE]      if(SPIN_UP and SPIN_DN) {
//[CO20200404 - OBSOLETE]        // Egap value
//[CO20200404 - OBSOLETE]        for(uint ii=0; ii<2; ii++) {
//[CO20200404 - OBSOLETE]          Egap.at(ii)     = CBB.at(ii) - VBT.at(ii);
//[CO20200404 - OBSOLETE]          Egap_fit.at(ii) = 1.348 * Egap.at(ii) + 0.913;
//[CO20200404 - OBSOLETE]          if(Egap.at(ii)     < _ENER_EPS) Egap.at(ii)     = _ZERO;
//[CO20200404 - OBSOLETE]          if(Egap_fit.at(ii) < _ENER_EPS) Egap_fit.at(ii) = _ZERO;
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]        // Gap types
//[CO20200404 - OBSOLETE]        if      (kptdist(1) <= kpt_tol) Egap_type.at(0) = "insulator_direct";
//[CO20200404 - OBSOLETE]        else if(kptdist(1) > kpt_tol) Egap_type.at(0) = "insulator_indirect";
//[CO20200404 - OBSOLETE]        if      (kptdist(2) <= kpt_tol) Egap_type.at(1) = "insulator_direct";
//[CO20200404 - OBSOLETE]        else if(kptdist(2) >  kpt_tol) Egap_type.at(1) = "insulator_indirect";
//[CO20200404 - OBSOLETE]        //[CO20171002 - GARBAGE A != A^-1]if      (kptdist(1) <= CELLVOLCUTOFF) Egap_type.at(0) = "insulator_direct";
//[CO20200404 - OBSOLETE]        //[CO20171002 - GARBAGE A != A^-1]else if(kptdist(1) >  CELLVOLCUTOFF) Egap_type.at(0) = "insulator_indirect";
//[CO20200404 - OBSOLETE]        //[CO20171002 - GARBAGE A != A^-1]if      (kptdist(2) <= CELLVOLCUTOFF) Egap_type.at(1) = "insulator_direct";
//[CO20200404 - OBSOLETE]        //[CO20171002 - GARBAGE A != A^-1]else if(kptdist(2) >  CELLVOLCUTOFF) Egap_type.at(1) = "insulator_indirect";
//[CO20200404 - OBSOLETE]        ///////////////////////////////////////////////////
//[CO20200404 - OBSOLETE]        // SET 1: VBT1<CBB1 && VBT2=CBB2
//[CO20200404 - OBSOLETE]        if((CBB.at(0)-VBT.at(0)>_ENER_EPS) and (abs(CBB.at(1)-VBT.at(1))<_ENER_EPS)) {
//[CO20200404 - OBSOLETE]          if((VBT.at(0)-VBT.at(1)>_ENER_EPS) and (VBT.at(0)-CBB.at(1)>_ENER_EPS)) { // 1a
//[CO20200404 - OBSOLETE]            Egap_type_net = "metal";
//[CO20200404 - OBSOLETE]            Egap_net      = _METALGAP_ ;
//[CO20200404 - OBSOLETE]            Egap_fit_net  = _METALGAP_ ;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          else if((abs(VBT.at(0)-VBT.at(1))<_ENER_EPS) and (abs(VBT.at(0)-CBB.at(1))<_ENER_EPS)) { // 1b
//[CO20200404 - OBSOLETE]            if  (abs(kptdist(1)-kptdist(3)) < kpt_tol) Egap_type_net = Egap_type.at(1);
//[CO20200404 - OBSOLETE]            else if(kptdist(1)-kptdist(3)  < kpt_tol) Egap_type_net = Egap_type.at(1);
//[CO20200404 - OBSOLETE]            else if(kptdist(1)-kptdist(3)  > kpt_tol) {
//[CO20200404 - OBSOLETE]              //[CO20171002 - GARBAGE A != A^-1]if  (abs(kptdist(1)-kptdist(3)) < CELLVOLCUTOFF) Egap_type_net = Egap_type.at(1);
//[CO20200404 - OBSOLETE]              //[CO20171002 - GARBAGE A != A^-1]else if(kptdist(1)-kptdist(3)  < CELLVOLCUTOFF) Egap_type_net = Egap_type.at(1);
//[CO20200404 - OBSOLETE]              //[CO20171002 - GARBAGE A != A^-1]else if(kptdist(1)-kptdist(3)  > CELLVOLCUTOFF) {  //[CO20200106 - close bracket for indenting]}
//[CO20200404 - OBSOLETE]              if      (Egap_type.at(0) == "insulator_indirect") Egap_type_net = "insulator_direct"  ;
//[CO20200404 - OBSOLETE]              else if(Egap_type.at(0) == "insulator_direct")   Egap_type_net = "insulator_indirect";
//[CO20200404 - OBSOLETE]            }
//[CO20200404 - OBSOLETE]            Egap_net     = _ZERO;
//[CO20200404 - OBSOLETE]            Egap_fit_net = _ZERO;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          else if((VBT.at(1)-VBT.at(0)>_ENER_EPS) and (CBB.at(0)-CBB.at(1)>_ENER_EPS)) { // 1c
//[CO20200404 - OBSOLETE]            Egap_type_net = Egap_type.at(1);
//[CO20200404 - OBSOLETE]            Egap_net     = _ZERO;
//[CO20200404 - OBSOLETE]            Egap_fit_net = _ZERO;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          else if((VBT.at(1)-VBT.at(0)>_ENER_EPS) and (abs(VBT.at(1)-CBB.at(0))<_ENER_EPS)) { // 1d
//[CO20200404 - OBSOLETE]            if  (abs(kptdist(2)-kptdist(4)) < kpt_tol) Egap_type_net = Egap_type.at(1);
//[CO20200404 - OBSOLETE]            else if(kptdist(2)-kptdist(4)  < kpt_tol) Egap_type_net = Egap_type.at(1);
//[CO20200404 - OBSOLETE]            else if(kptdist(2)-kptdist(4)  > kpt_tol) {
//[CO20200404 - OBSOLETE]              //[CO20171002 - GARBAGE A != A^-1]if  (abs(kptdist(2)-kptdist(4)) < CELLVOLCUTOFF) Egap_type_net = Egap_type.at(1);
//[CO20200404 - OBSOLETE]              //[CO20171002 - GARBAGE A != A^-1]else if(kptdist(2)-kptdist(4)  < CELLVOLCUTOFF) Egap_type_net = Egap_type.at(1);
//[CO20200404 - OBSOLETE]              //[CO20171002 - GARBAGE A != A^-1]else if(kptdist(2)-kptdist(4)  > CELLVOLCUTOFF) {  //[CO20200106 - close bracket for indenting]}
//[CO20200404 - OBSOLETE]              if      (Egap_type.at(1) == "insulator_indirect") Egap_type_net = "insulator_direct"  ;
//[CO20200404 - OBSOLETE]              else if(Egap_type.at(1) == "insulator_direct")   Egap_type_net = "insulator_indirect";
//[CO20200404 - OBSOLETE]            }
//[CO20200404 - OBSOLETE]            Egap_net     = _ZERO;
//[CO20200404 - OBSOLETE]            Egap_fit_net = _ZERO;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          if((VBT.at(1)-CBB.at(0)>_ENER_EPS) and (CBB.at(1)-CBB.at(0)>_ENER_EPS)) { // 1e
//[CO20200404 - OBSOLETE]            Egap_type_net = "metal";
//[CO20200404 - OBSOLETE]            Egap_net      = _METALGAP_ ;
//[CO20200404 - OBSOLETE]            Egap_fit_net  = _METALGAP_ ;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]        ///////////////////////////////////////////////////
//[CO20200404 - OBSOLETE]        // SET 2: VBT1=CBB1 && VBT2<CBB2
//[CO20200404 - OBSOLETE]        else if((abs(CBB.at(0)-VBT.at(0))<_ENER_EPS) and (CBB.at(1)-VBT.at(1)>_ENER_EPS)) {
//[CO20200404 - OBSOLETE]          if((VBT.at(1)-VBT.at(0)>_ENER_EPS) and (VBT.at(1)-CBB.at(0)>_ENER_EPS)) { // 2a
//[CO20200404 - OBSOLETE]            Egap_type_net = "metal";
//[CO20200404 - OBSOLETE]            Egap_net      = _METALGAP_ ;
//[CO20200404 - OBSOLETE]            Egap_fit_net  = _METALGAP_ ;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          else if((abs(VBT.at(1)-VBT.at(0))<_ENER_EPS) and (abs(VBT.at(1)-CBB.at(0))<_ENER_EPS)) { // 2b
//[CO20200404 - OBSOLETE]            if  (abs(kptdist(1)-kptdist(4)) < kpt_tol) Egap_type_net = Egap_type.at(0);
//[CO20200404 - OBSOLETE]            else if(kptdist(1)-kptdist(4)  < kpt_tol) Egap_type_net = Egap_type.at(0);
//[CO20200404 - OBSOLETE]            else if(kptdist(1)-kptdist(4)  > kpt_tol) {
//[CO20200404 - OBSOLETE]              //[CO20171002 - GARBAGE A != A^-1]if  (abs(kptdist(1)-kptdist(4)) < CELLVOLCUTOFF) Egap_type_net = Egap_type.at(0);
//[CO20200404 - OBSOLETE]              //[CO20171002 - GARBAGE A != A^-1]else if(kptdist(1)-kptdist(4)  < CELLVOLCUTOFF) Egap_type_net = Egap_type.at(0);
//[CO20200404 - OBSOLETE]              //[CO20171002 - GARBAGE A != A^-1]else if(kptdist(1)-kptdist(4)  > CELLVOLCUTOFF) {  //[CO20200106 - close bracket for indenting]}
//[CO20200404 - OBSOLETE]              if      (Egap_type.at(0) == "insulator_indirect") Egap_type_net = "insulator_direct"  ;
//[CO20200404 - OBSOLETE]              else if(Egap_type.at(0) == "insulator_direct")   Egap_type_net = "insulator_indirect";
//[CO20200404 - OBSOLETE]            }
//[CO20200404 - OBSOLETE]            Egap_net     = _ZERO;
//[CO20200404 - OBSOLETE]            Egap_fit_net = _ZERO;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          else if((VBT.at(0)-VBT.at(1)>_ENER_EPS) and (CBB.at(1)-CBB.at(0)>_ENER_EPS)) { // 2c
//[CO20200404 - OBSOLETE]            Egap_type_net = Egap_type.at(0);
//[CO20200404 - OBSOLETE]            Egap_net     = _ZERO;
//[CO20200404 - OBSOLETE]            Egap_fit_net = _ZERO;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          else if((VBT.at(0)-VBT.at(1)>_ENER_EPS) and (abs(VBT.at(0)-CBB.at(1))<_ENER_EPS)) { // 2d
//[CO20200404 - OBSOLETE]            if  (abs(kptdist(1)-kptdist(3)) < kpt_tol) Egap_type_net = Egap_type.at(0);
//[CO20200404 - OBSOLETE]            else if(kptdist(1)-kptdist(3)  < kpt_tol) Egap_type_net = Egap_type.at(0);
//[CO20200404 - OBSOLETE]            else if(kptdist(1)-kptdist(3)  > kpt_tol) {
//[CO20200404 - OBSOLETE]              //[CO20171002 - GARBAGE A != A^-1]if  (abs(kptdist(1)-kptdist(3)) < CELLVOLCUTOFF) Egap_type_net = Egap_type.at(0);
//[CO20200404 - OBSOLETE]              //[CO20171002 - GARBAGE A != A^-1]else if(kptdist(1)-kptdist(3)  < CELLVOLCUTOFF) Egap_type_net = Egap_type.at(0);
//[CO20200404 - OBSOLETE]              //[CO20171002 - GARBAGE A != A^-1]else if(kptdist(1)-kptdist(3)  > CELLVOLCUTOFF) {  //[CO20200106 - close bracket for indenting]}
//[CO20200404 - OBSOLETE]              if      (Egap_type.at(0) == "insulator_indirect") Egap_type_net = "insulator_direct"  ;
//[CO20200404 - OBSOLETE]              else if(Egap_type.at(0) == "insulator_direct"  ) Egap_type_net = "insulator_indirect";
//[CO20200404 - OBSOLETE]            }
//[CO20200404 - OBSOLETE]            Egap_net     = _ZERO;
//[CO20200404 - OBSOLETE]            Egap_fit_net = _ZERO;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          if((VBT.at(0)-CBB.at(1)>_ENER_EPS) and (CBB.at(0)-CBB.at(1)>_ENER_EPS)) { // 2e
//[CO20200404 - OBSOLETE]            Egap_type_net = "metal";
//[CO20200404 - OBSOLETE]            Egap_net      = _METALGAP_ ;
//[CO20200404 - OBSOLETE]            Egap_fit_net  = _METALGAP_ ;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]        ///////////////////////////////////////////////////
//[CO20200404 - OBSOLETE]        // SET 3: VBT1=CBB1 && VBT2=CBB2
//[CO20200404 - OBSOLETE]        else if( (abs(CBB.at(0)-VBT.at(0))<_ENER_EPS) and (abs(CBB.at(1)-VBT.at(1))<_ENER_EPS) ) { // CAMILOFIX
//[CO20200404 - OBSOLETE]          if((VBT.at(0)-VBT.at(1)>_ENER_EPS) and ((CBB.at(0)-CBB.at(1))>_ENER_EPS)) { // 3a
//[CO20200404 - OBSOLETE]            Egap_type_net = "metal";
//[CO20200404 - OBSOLETE]            Egap_net      = _METALGAP_ ;
//[CO20200404 - OBSOLETE]            Egap_fit_net  = _METALGAP_ ;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          else if( (abs(VBT.at(0)-VBT.at(1))<_ENER_EPS) and (abs(CBB.at(0)-CBB.at(1))<_ENER_EPS) ) { // 3b // CAMILOFIX
//[CO20200404 - OBSOLETE]            if((kptdist(1)<=kpt_tol) or (kptdist(2)<=kpt_tol) or
//[CO20200404 - OBSOLETE]                (kptdist(3)<=kpt_tol) or (kptdist(4)<=kpt_tol)) {
//[CO20200404 - OBSOLETE]              //[CO20171002 - GARBAGE A != A^-1]if((kptdist(1)<=CELLVOLCUTOFF) or (kptdist(2)<=CELLVOLCUTOFF) or
//[CO20200404 - OBSOLETE]              //[CO20171002 - GARBAGE A != A^-1]   (kptdist(3)<=CELLVOLCUTOFF) or (kptdist(4)<=CELLVOLCUTOFF)) { //[CO20200106 - close bracket for indenting]}
//[CO20200404 - OBSOLETE]              Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            }
//[CO20200404 - OBSOLETE]            else Egap_type_net = "insulator_indirect";
//[CO20200404 - OBSOLETE]            Egap_net     = _ZERO;
//[CO20200404 - OBSOLETE]            Egap_fit_net = _ZERO;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          if((VBT.at(1)-VBT.at(0)>_ENER_EPS) and ((CBB.at(1)-CBB.at(0))>_ENER_EPS)) { // 3c
//[CO20200404 - OBSOLETE]            Egap_type_net = "metal";
//[CO20200404 - OBSOLETE]            Egap_net      = _METALGAP_ ;
//[CO20200404 - OBSOLETE]            Egap_fit_net  = _METALGAP_ ;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]        ///////////////////////////////////////////////////
//[CO20200404 - OBSOLETE]        // SET 4: VBT1<CBB1 && VBT2<CBB2
//[CO20200404 - OBSOLETE]        else if((CBB.at(0)-VBT.at(0)>_ENER_EPS) and (CBB.at(1)-VBT.at(1)>_ENER_EPS)) {
//[CO20200404 - OBSOLETE]          if(VBT.at(0)-CBB.at(1)>_ENER_EPS) { // 4a
//[CO20200404 - OBSOLETE]            Egap_type_net = "metal";
//[CO20200404 - OBSOLETE]            Egap_net      = _METALGAP_ ;
//[CO20200404 - OBSOLETE]            Egap_fit_net  = _METALGAP_ ;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          else if(VBT.at(1)-CBB.at(0)>_ENER_EPS) { // 4b
//[CO20200404 - OBSOLETE]            Egap_type_net = "metal";
//[CO20200404 - OBSOLETE]            Egap_net      = _METALGAP_ ;
//[CO20200404 - OBSOLETE]            Egap_fit_net  = _METALGAP_ ;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          else if(abs(VBT.at(0)-CBB.at(1))<_ENER_EPS) { // 4c
//[CO20200404 - OBSOLETE]            if(kptdist(3)<=kpt_tol) Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            //[CO20171002 - GARBAGE A != A^-1]if(kptdist(3)<=CELLVOLCUTOFF) Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            else Egap_type_net = "insulator_indirect";
//[CO20200404 - OBSOLETE]            Egap_net     = _ZERO;
//[CO20200404 - OBSOLETE]            Egap_fit_net = 1.348 * Egap_net + 0.913;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          else if(abs(VBT.at(1)-CBB.at(0))<_ENER_EPS) { // 4d
//[CO20200404 - OBSOLETE]            if(kptdist(4)<=kpt_tol) Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            //[CO20171002 - GARBAGE A != A^-1]if(kptdist(4)<=CELLVOLCUTOFF) Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            else Egap_type_net = "insulator_indirect";
//[CO20200404 - OBSOLETE]            Egap_net     = _ZERO;
//[CO20200404 - OBSOLETE]            Egap_fit_net = 1.348 * Egap_net + 0.913;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          else if((CBB.at(0)-CBB.at(1)>_ENER_EPS) and (VBT.at(0)-VBT.at(1)>_ENER_EPS)) { // 4e
//[CO20200404 - OBSOLETE]            if(kptdist(3)<=kpt_tol) Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            //[CO20171002 - GARBAGE A != A^-1]if(kptdist(3)<=CELLVOLCUTOFF) Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            else Egap_type_net   = "insulator_indirect";
//[CO20200404 - OBSOLETE]            Egap_net     = CBB.at(1) - VBT.at(0);
//[CO20200404 - OBSOLETE]            Egap_fit_net = 1.348 * Egap_net + 0.913;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          else if((CBB.at(1)-CBB.at(0)>_ENER_EPS) and (VBT.at(1)-VBT.at(0)>_ENER_EPS)) { // 4f
//[CO20200404 - OBSOLETE]            if(kptdist(4)<=kpt_tol) Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            //[CO20171002 - GARBAGE A != A^-1]if(kptdist(4)<=CELLVOLCUTOFF) Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            else Egap_type_net = "insulator_indirect";
//[CO20200404 - OBSOLETE]            Egap_net     = CBB.at(0) - VBT.at(1);
//[CO20200404 - OBSOLETE]            Egap_fit_net = 1.348 * Egap_net + 0.913;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          else if((VBT.at(0)-VBT.at(1)>_ENER_EPS) and (abs(CBB.at(0)-CBB.at(1))<_ENER_EPS)) { // 4g
//[CO20200404 - OBSOLETE]            if(kptdist(1)<=kpt_tol) Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            //[CO20171002 - GARBAGE A != A^-1]if(kptdist(1)<=CELLVOLCUTOFF) Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            else Egap_type_net = "insulator_indirect";
//[CO20200404 - OBSOLETE]            Egap_net     = CBB.at(0) - VBT.at(0);
//[CO20200404 - OBSOLETE]            Egap_fit_net = 1.348 * Egap_net + 0.913;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          else if((VBT.at(1)-VBT.at(0)>_ENER_EPS) and (abs(CBB.at(0)-CBB.at(1))<_ENER_EPS)) { // 4h
//[CO20200404 - OBSOLETE]            if(kptdist(2)<=kpt_tol) Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            //[CO20171002 - GARBAGE A != A^-1]if(kptdist(2)<=CELLVOLCUTOFF) Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            else Egap_type_net = "insulator_indirect";
//[CO20200404 - OBSOLETE]            Egap_net     = CBB.at(1) - VBT.at(1);
//[CO20200404 - OBSOLETE]            Egap_fit_net = 1.348 * Egap_net + 0.913;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          else if((VBT.at(0)-VBT.at(1)>_ENER_EPS) and (CBB.at(1)-CBB.at(0)>_ENER_EPS)) { // 4i
//[CO20200404 - OBSOLETE]            if(kptdist(1)<=kpt_tol) Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            //[CO20171002 - GARBAGE A != A^-1]if(kptdist(1)<=CELLVOLCUTOFF) Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            else Egap_type_net = "insulator_indirect";
//[CO20200404 - OBSOLETE]            Egap_net     = CBB.at(0) - VBT.at(0);
//[CO20200404 - OBSOLETE]            Egap_fit_net = 1.348 * Egap_net + 0.913;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          else if((VBT.at(1)-VBT.at(0)>_ENER_EPS) and (CBB.at(0)-CBB.at(1)>_ENER_EPS)) { // 4j
//[CO20200404 - OBSOLETE]            if(kptdist(2)<=kpt_tol) Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            //[CO20171002 - GARBAGE A != A^-1]if(kptdist(2)<=CELLVOLCUTOFF) Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            else Egap_type_net = "insulator_indirect";
//[CO20200404 - OBSOLETE]            Egap_net     = CBB.at(1) - VBT.at(1);
//[CO20200404 - OBSOLETE]            Egap_fit_net = 1.348 * Egap_net + 0.913;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          else if( (abs(VBT.at(0)-VBT.at(1))<_ENER_EPS) and (abs(CBB.at(0)-CBB.at(1))<_ENER_EPS) ) { // 4k // CAMILOFIX
//[CO20200404 - OBSOLETE]            if((kptdist(1)<=kpt_tol) or (kptdist(2)<=kpt_tol) or
//[CO20200404 - OBSOLETE]                (kptdist(3)<=kpt_tol) or (kptdist(4)<=kpt_tol)) {
//[CO20200404 - OBSOLETE]              //[CO20171002 - GARBAGE A != A^-1]if((kptdist(1)<=CELLVOLCUTOFF) or (kptdist(2)<=CELLVOLCUTOFF) or
//[CO20200404 - OBSOLETE]              //[CO20171002 - GARBAGE A != A^-1]   (kptdist(3)<=CELLVOLCUTOFF) or (kptdist(4)<=CELLVOLCUTOFF)) { //[CO20200106 - close bracket for indenting]}
//[CO20200404 - OBSOLETE]              Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            }
//[CO20200404 - OBSOLETE]            else Egap_type_net = "insulator_indirect";
//[CO20200404 - OBSOLETE]            Egap_net     = CBB.at(0) - VBT.at(0);
//[CO20200404 - OBSOLETE]            Egap_fit_net = 1.348 * Egap_net + 0.913;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          else if((abs(VBT.at(0)-VBT.at(1))<_ENER_EPS) and (CBB.at(1)-CBB.at(0)>_ENER_EPS)) { // 4l
//[CO20200404 - OBSOLETE]            if(kptdist(1)<=kpt_tol) Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            //[CO20171002 - GARBAGE A != A^-1]if(kptdist(1)<=CELLVOLCUTOFF) Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            else Egap_type_net = "insulator_indirect";
//[CO20200404 - OBSOLETE]            Egap_net     = CBB.at(0) - VBT.at(0);
//[CO20200404 - OBSOLETE]            Egap_fit_net = 1.348 * Egap_net + 0.913;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]          else if((abs(VBT.at(0)-VBT.at(1))<_ENER_EPS) and (CBB.at(0)-CBB.at(1)>_ENER_EPS)) { // 4m
//[CO20200404 - OBSOLETE]            if(kptdist(2)<=kpt_tol) Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            //[CO20171002 - GARBAGE A != A^-1]if(kptdist(2)<=CELLVOLCUTOFF) Egap_type_net = "insulator_direct";
//[CO20200404 - OBSOLETE]            else Egap_type_net = "insulator_indirect";
//[CO20200404 - OBSOLETE]            Egap_net     = CBB.at(1) - VBT.at(1);
//[CO20200404 - OBSOLETE]            Egap_fit_net = 1.348 * Egap_net + 0.913;
//[CO20200404 - OBSOLETE]          }
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]      else if( SPIN_UP and !SPIN_DN) { // m1
//[CO20200404 - OBSOLETE]        if      (kptdist(1) <= kpt_tol) Egap_type.at(0) = "insulator_direct";
//[CO20200404 - OBSOLETE]        else if(kptdist(1) >  kpt_tol) Egap_type.at(0) = "insulator_indirect";
//[CO20200404 - OBSOLETE]        //[CO20171002 - GARBAGE A != A^-1]if      (kptdist(1) <= CELLVOLCUTOFF) Egap_type.at(0) = "insulator_direct";
//[CO20200404 - OBSOLETE]        //[CO20171002 - GARBAGE A != A^-1]else if(kptdist(1) >  CELLVOLCUTOFF) Egap_type.at(0) = "insulator_indirect";
//[CO20200404 - OBSOLETE]        conduction_band_min.at(1) = _METALEDGE_;
//[CO20200404 - OBSOLETE]        valence_band_max.at(1)    = _METALEDGE_;
//[CO20200404 - OBSOLETE]        conduction_band_min_net = conduction_band_min.at(1);
//[CO20200404 - OBSOLETE]        valence_band_max_net    = valence_band_max.at(1);
//[CO20200404 - OBSOLETE]        Egap_type.at(1) = "metal";
//[CO20200404 - OBSOLETE]        Egap_type_net   = Egap_type.at(1);
//[CO20200404 - OBSOLETE]        Egap.at(0)      = CBB.at(0) - VBT.at(0);
//[CO20200404 - OBSOLETE]        Egap.at(1)      = _METALGAP_;
//[CO20200404 - OBSOLETE]        Egap_net        = Egap.at(1);
//[CO20200404 - OBSOLETE]        //WSETYAWAN fit
//[CO20200404 - OBSOLETE]        Egap_fit.at(0)  = 1.348 * Egap.at(0) + 0.913;
//[CO20200404 - OBSOLETE]        Egap_fit.at(1)  = _METALGAP_;
//[CO20200404 - OBSOLETE]        Egap_fit_net    = Egap_fit.at(1);
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]      else if(!SPIN_UP and  SPIN_DN) { // m2
//[CO20200404 - OBSOLETE]        if      (kptdist(2) <= kpt_tol) Egap_type.at(1) = "insulator_direct";
//[CO20200404 - OBSOLETE]        else if(kptdist(2) >  kpt_tol) Egap_type.at(1) = "insulator_indirect";
//[CO20200404 - OBSOLETE]        //[CO20171002 - GARBAGE A != A^-1]if      (kptdist(2) <= CELLVOLCUTOFF) Egap_type.at(1) = "insulator_direct";
//[CO20200404 - OBSOLETE]        //[CO20171002 - GARBAGE A != A^-1]else if(kptdist(2) >  CELLVOLCUTOFF) Egap_type.at(1) = "insulator_indirect";
//[CO20200404 - OBSOLETE]        conduction_band_min.at(0) = _METALEDGE_;
//[CO20200404 - OBSOLETE]        valence_band_max.at(0)    = _METALEDGE_;
//[CO20200404 - OBSOLETE]        conduction_band_min_net = conduction_band_min.at(0);
//[CO20200404 - OBSOLETE]        valence_band_max_net    = valence_band_max.at(0);
//[CO20200404 - OBSOLETE]        Egap_type.at(0) = "metal";
//[CO20200404 - OBSOLETE]        Egap_type_net   = Egap_type.at(0);
//[CO20200404 - OBSOLETE]        Egap.at(0)      = _METALGAP_;
//[CO20200404 - OBSOLETE]        Egap.at(1)      = CBB.at(1) - VBT.at(1);
//[CO20200404 - OBSOLETE]        Egap_net        = Egap.at(0);
//[CO20200404 - OBSOLETE]        //WSETYAWAN fit
//[CO20200404 - OBSOLETE]        Egap_fit.at(0)  = _METALGAP_;
//[CO20200404 - OBSOLETE]        Egap_fit.at(1)  = 1.348 * Egap.at(1) + 0.913;
//[CO20200404 - OBSOLETE]        Egap_fit_net    = Egap_fit.at(0);
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]      // Full metal
//[CO20200404 - OBSOLETE]      else if(!SPIN_UP and !SPIN_DN) { // m3
//[CO20200404 - OBSOLETE]        for(int ii=0; ii<(int)ISPIN; ii++) {
//[CO20200404 - OBSOLETE]          conduction_band_min.at(ii) = _METALEDGE_;
//[CO20200404 - OBSOLETE]          valence_band_max.at(ii) = _METALEDGE_;
//[CO20200404 - OBSOLETE]        }
//[CO20200404 - OBSOLETE]        conduction_band_min_net = _METALEDGE_;
//[CO20200404 - OBSOLETE]        valence_band_max_net    = _METALEDGE_;
//[CO20200404 - OBSOLETE]        Egap.at(0) = _METALGAP_;
//[CO20200404 - OBSOLETE]        Egap.at(1) = _METALGAP_;
//[CO20200404 - OBSOLETE]        Egap_type.at(0) = "metal";
//[CO20200404 - OBSOLETE]        Egap_type.at(1) = "metal";
//[CO20200404 - OBSOLETE]        Egap_type_net   = "metal";
//[CO20200404 - OBSOLETE]        Egap_net = _METALGAP_;
//[CO20200404 - OBSOLETE]        //WSETYAWAN fit
//[CO20200404 - OBSOLETE]        Egap_fit.at(0) = _METALGAP_;
//[CO20200404 - OBSOLETE]        Egap_fit.at(1) = _METALGAP_;
//[CO20200404 - OBSOLETE]        Egap_fit_net   = _METALGAP_;
//[CO20200404 - OBSOLETE]      }
//[CO20200404 - OBSOLETE]    }
//[CO20200404 - OBSOLETE]  }
//[CO20200404 - OBSOLETE]  // ----------------------------------------------------------------------
//[CO20200404 - OBSOLETE]  // DONE NOW RETURN
//[CO20200404 - OBSOLETE]  //  if(ERROR_flag && !QUIET) return FALSE;
//[CO20200404 - OBSOLETE]  return TRUE;
//[CO20200404 - OBSOLETE]}

/// @brief gets the ionic step data
/// @authors
/// @mod{CO,20211106,created function}
/// @mod{SD,20221208,patched for VASP5 and VASP6}
bool xOUTCAR::GetIonicStepsData(){
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  vxstr_ionic.clear();
  venergy_ionic.clear();
  vstresses_ionic.clear();

  bool reading_ionic=false;
  bool reading_stresses=false;
  bool reading_lattice=false;
  bool reading_atoms=false;
  bool convert_kBar2eV=false;

  //get species data first
  uint iline=0,i=0;
  vector<string> vtokens;
  string tmp_str="";
  deque<string> species_pp;
  deque<int> num_each_type;
  uint vcontent_size=vcontent.size();
  for(iline=0;iline<vcontent_size;iline++){
    if(aurostd::substring2bool(vcontent[iline],"POTCAR:")){
      tmp_str=vcontent[iline].substr(vcontent[iline].find("POTCAR:")+6+1);  //len(POSCAR:)==6
      aurostd::string2tokens(tmp_str,vtokens," ");
      //vtokens[1] is the species_pp
      if(!aurostd::WithinList(species_pp,vtokens[1])){species_pp.push_back(vtokens[1]);}
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " vcontent[iline=" << iline << "]=\"" << vcontent[iline] << "\"" << endl;
      }
    }
    if(aurostd::substring2bool(vcontent[iline],"ions") && //THIS MUST COME AFTER species_pp IS FILLED
        aurostd::substring2bool(vcontent[iline],"per") &&
        aurostd::substring2bool(vcontent[iline],"type")){
      aurostd::string2tokens(vcontent[iline],vtokens,"=");
      if(vtokens.size()>1){
        tmp_str=vtokens[1];
        aurostd::string2tokens(tmp_str,vtokens," ");
        if(vtokens.size()>0 && vtokens.size()==species_pp.size()){
          num_each_type=aurostd::vector2deque(aurostd::vectorstring2vectorutype<int>(vtokens));
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << " vcontent[iline=" << iline << "]=\"" << vcontent[iline] << "\"" << endl;
          }
          break;
        }
      }
    }
  }
  uint natoms=aurostd::sum(num_each_type);
  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " species_pp=" << aurostd::joinWDelimiter(species_pp,",") << endl;
    cerr << __AFLOW_FUNC__ << " num_each_type=" << aurostd::joinWDelimiter(num_each_type,",") << endl;
    cerr << __AFLOW_FUNC__ << " natoms=" << natoms << endl;
  }

  uint ilattice=0,iatom=0,itype=0;
  xstructure xstr;
  double energy=AUROSTD_MAX_DOUBLE;
  xvector<double> stresses(6);stresses[stresses.lrows]=AUROSTD_MAX_DOUBLE;
  vector<_atom> vatoms; //must keep like this until we settle types...
  for(iline=0;iline<vcontent_size;iline++){
    if(reading_ionic==false){
      if(aurostd::substring2bool(vcontent[iline],"aborting") &&
           aurostd::substring2bool(vcontent[iline],"loop") &&
           aurostd::substring2bool(vcontent[iline],"because") &&
           aurostd::substring2bool(vcontent[iline],"EDIFF") &&
           aurostd::substring2bool(vcontent[iline],"is") &&
           aurostd::substring2bool(vcontent[iline],"reached")){
        reading_ionic=true;
        reading_stresses=false;reading_lattice=false;reading_atoms=false;convert_kBar2eV=false;
        //clear everything, set stresses to 3 as an indicator that it's not set
        xstr.clear();vatoms.clear();
        ilattice=0;iatom=0;
        energy=AUROSTD_MAX_DOUBLE;
        stresses[stresses.lrows]=AUROSTD_MAX_DOUBLE; //stresses.resize(3);
        continue;
      }
    }else{
      if(aurostd::substring2bool(vcontent[iline],"Iteration")||(iline==vcontent_size-1)){
        reading_ionic=false;
        reading_stresses=false;
        if(xstr.atoms.size()==natoms && ilattice==3 && iatom==natoms && energy!=AUROSTD_MAX_DOUBLE && stresses[stresses.lrows]!=AUROSTD_MAX_DOUBLE){  //add additional checks as needed //stresses.rows==6
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " adding a new ionic step" << endl;}
          vxstr_ionic.push_back(xstr);
          venergy_ionic.push_back(energy);
          vstresses_ionic.push_back(stresses);
        }
        continue;
      }
      if(aurostd::substring2bool(vcontent[iline],"energy") && 
          aurostd::substring2bool(vcontent[iline],"without") &&
          aurostd::substring2bool(vcontent[iline],"entropy")){
        aurostd::string2tokens(vcontent[iline],vtokens,"=");
        if(vtokens.size()>1){
          tmp_str=vtokens[1];
          aurostd::string2tokens(tmp_str,vtokens,"e");
          if(vtokens.size()>0){energy=aurostd::string2utype<double>(vtokens[0]);}
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << " vcontent[iline=" << iline << "]=\"" << vcontent[iline] << "\"" << endl;
            cerr << __AFLOW_FUNC__ << " energy=" << energy << endl;
          }
        }
      }
      if(aurostd::substring2bool(vcontent[iline],"FORCE") && aurostd::substring2bool(vcontent[iline],"-STRESS")){reading_stresses=true;reading_lattice=false;reading_atoms=false;convert_kBar2eV=false;continue;}
      if(reading_stresses){
        if(aurostd::substring2bool(vcontent[iline],"Total")){
          aurostd::string2tokens(vcontent[iline],vtokens," ");
          if(vtokens.size()==7){
            //stresses.resize(6);
            for(i=1;i<7;i++){stresses[stresses.lrows+i-1]=aurostd::string2utype<double>(vtokens[i]);}
            if(LDEBUG){
              cerr << __AFLOW_FUNC__ << " vcontent[iline=" << iline << "]=\"" << vcontent[iline] << "\"" << endl;
              cerr << __AFLOW_FUNC__ << " stresses(eV)=" << stresses << endl; //X Y Z XY YZ ZX //this is eV*volume of cell
            }
          }
        }
        //per ME's suggestion, let's grab 'in kB' instead of 'Total' (eV) which has more significant figures
        if(aurostd::substring2bool(vcontent[iline],"in") && aurostd::substring2bool(vcontent[iline],"kB")){
          aurostd::string2tokens(vcontent[iline],vtokens," ");
          if(vtokens.size()==8){
            //stresses.resize(6);
            for(i=2;i<8;i++){stresses[stresses.lrows+i-2]=aurostd::string2utype<double>(vtokens[i]);}
            convert_kBar2eV=true;
            if(LDEBUG){
              cerr << __AFLOW_FUNC__ << " vcontent[iline=" << iline << "]=\"" << vcontent[iline] << "\"" << endl;
              cerr << __AFLOW_FUNC__ << " stresses(kBar)=" << stresses << endl; //X Y Z XY YZ ZX
            }
          }
        }
      }
      if(aurostd::substring2bool(vcontent[iline],"direct") && 
          aurostd::substring2bool(vcontent[iline],"lattice") &&
          aurostd::substring2bool(vcontent[iline],"vectors") &&
          aurostd::substring2bool(vcontent[iline],"reciprocal")){reading_lattice=true;reading_stresses=false;reading_atoms=false;ilattice=0;xstr.clear();continue;}
      if(reading_lattice && ilattice<3){
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " vcontent[iline=" << iline << "]=\"" << vcontent[iline] << "\"" << endl;}
        vtokens=GetCorrectPositions(vcontent[iline],6); //aurostd::string2tokens(vcontent[iline],vtokens," ");
        if(vtokens.size()!=6){ilattice=0;reading_lattice=false;}  //also contains reciprocal lattice components //ilattice==0 is an indicator that lattice was not read (correctly)
        xstr.lattice[xstr.lattice.lrows+ilattice][xstr.lattice.lcols+0]=aurostd::string2utype<double>(vtokens[0]);
        xstr.lattice[xstr.lattice.lrows+ilattice][xstr.lattice.lcols+1]=aurostd::string2utype<double>(vtokens[1]);
        xstr.lattice[xstr.lattice.lrows+ilattice][xstr.lattice.lcols+2]=aurostd::string2utype<double>(vtokens[2]);
        ilattice++;
        if(ilattice==3){
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " xstr.lattice=" << endl << xstr.lattice << endl;}
          xstr.scale=1.0;
          xstr.FixLattices();
          if(convert_kBar2eV){
            if(LDEBUG){
              cerr << __AFLOW_FUNC__ << " xstr.GetVolume()=" << xstr.GetVolume() << endl;
              cerr << __AFLOW_FUNC__ << " stresses(kBar)=" << stresses << endl; //X Y Z XY YZ ZX
              cerr << __AFLOW_FUNC__ << " kBar2eV=" << kBar2eV << endl;
            }
            stresses*=(kBar2eV*xstr.GetVolume());
            if(LDEBUG){
              cerr << __AFLOW_FUNC__ << " stresses(eV)=" << stresses << endl; //X Y Z XY YZ ZX //this is eV*volume of cell
            }
          }
          reading_lattice=false;
        }
      }
      if(aurostd::substring2bool(vcontent[iline],"POSITION") && aurostd::substring2bool(vcontent[iline],"TOTAL-FORCE")){reading_atoms=true;reading_stresses=false;reading_lattice=false;convert_kBar2eV=false;iatom=0;vatoms.clear();continue;}
      if(reading_atoms && iatom<natoms){
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " vcontent[iline=" << iline << "]=\"" << vcontent[iline] << "\"" << endl;}
        if(aurostd::substring2bool(vcontent[iline],"---------")){continue;}
        vtokens=GetCorrectPositions(vcontent[iline],6); //aurostd::string2tokens(vcontent[iline],vtokens," ");
        if(vtokens.size()!=6){iatom=0;reading_atoms=false;}  //also contains force components //iatom==0 is an indicator that atoms were not read (correctly)
        vatoms.push_back(_atom());
        for(i=0;i<3;i++){vatoms.back().cpos[vatoms.back().cpos.lrows+i]=aurostd::string2utype<double>(vtokens[i]);}
        for(i=0;i<3;i++){vatoms.back().force[vatoms.back().force.lrows+i]=aurostd::string2utype<double>(vtokens[i+3]);}
        iatom++;
        if(iatom==natoms){
          iatom=0; //reset
          for(itype=0;itype<num_each_type.size();itype++){
            for(i=0;i<(uint)num_each_type[itype];i++){
              vatoms[iatom].type=itype;
              vatoms[iatom].name=species_pp[itype];
              vatoms[iatom].name_is_given=(!vatoms[iatom].name.empty());
              vatoms[iatom].fpos=xstr.c2f*vatoms[iatom].cpos;
              xstr.AddAtom(vatoms[iatom],false);  //CO20230319 - add by type
              xstr.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
              iatom++;
            }
          }
          xstr.title=SYSTEM;
          xstr.MakeBasis();
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " xstr=" << endl << xstr << endl;}
          reading_atoms=false;
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
  if (vxstr_ionic.size() != venergy_ionic.size()) {throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "vxstr_ionic.size()!=venergy_ionic", _FILE_CORRUPT_);}
  if (vxstr_ionic.size() != vstresses_ionic.size()) {throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "vxstr_ionic.size()!=vstresses_ionic", _FILE_CORRUPT_);}
  //bool FORMATION_CALC = false;
  vector<string> type;
  vector<xvector<double>> fpos, force;
  for(size_t istr = 0; istr < vxstr_ionic.size(); istr++) {
    aurostd::JSON::object step(aurostd::JSON::object_types::DICTIONARY);
    step["auid"] = entry.auid;
    step["energy"] = venergy_ionic[istr];
    step["lattice"] = vxstr_ionic[istr].lattice;
    step["stress"] = vstresses_ionic[istr];
    type.clear();
    fpos.clear();
    force.clear();
    for(size_t iatom = 0; iatom < vxstr_ionic[istr].atoms.size(); iatom++) {
      type.push_back(vxstr_ionic[istr].atoms[iatom].cleanname);
      fpos.push_back(vxstr_ionic[istr].atoms[iatom].fpos);
      force.push_back(vxstr_ionic[istr].atoms[iatom].force);
    }
    aurostd::JSON::object atoms(aurostd::JSON::object_types::DICTIONARY);
    atoms["type"] = type;
    atoms["position"] = fpos;
    atoms["force"] = force;
    step["atoms"] = atoms;
  //  entry.enthalpy_cell = entry.enthalpy_atom = venergy_ionic[istr];
  //  entry.enthalpy_atom /= vxstr_ionic[istr].atoms.size();
  //  FORMATION_CALC = aflowlib::LIB2RAW_Calculate_FormationEnthalpy(entry, vxstr_ionic[istr], __AFLOW_FUNC__);
  //  if (FORMATION_CALC) {
  //    step["enthalpy_formation_atom"] = entry.enthalpy_formation_atom;
  //  }
    jo["data"].push_back(step);
  }

}

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************

// ***************************************************************************
// class xDOSCAR
//ME20200427 - included xStream::initialize
xDOSCAR::xDOSCAR(ostream& oss) : xStream(oss) {free();} //CO20200404 - xStream integration for logging
xDOSCAR::xDOSCAR(ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss) {free();} //CO20200404 - xStream integration for logging
xDOSCAR::xDOSCAR(const string& fileIN,bool QUIET,ostream& oss) : xStream(oss) {m_initialized=initialize(fileIN,QUIET);}  //CO20200404 - xStream integration for logging
xDOSCAR::xDOSCAR(const string& fileIN,ofstream& FileMESSAGE,bool QUIET,ostream& oss) : xStream(FileMESSAGE,oss) {m_initialized=initialize(fileIN,QUIET);}  //CO20200404 - xStream integration for logging

bool xDOSCAR::initialize(const string& fileIN,ostream& oss,bool QUIET) {xStream::initialize(oss);return initialize(fileIN,QUIET);} //CO20200508
bool xDOSCAR::initialize(const string& fileIN,ofstream& FileMESSAGE,ostream& oss,bool QUIET) {xStream::initialize(FileMESSAGE,oss);return initialize(fileIN,QUIET);} //CO20200508
bool xDOSCAR::initialize(const string& fileIN, bool QUIET) {
  free();
  filename = fileIN;
  return GetPropertiesFile(fileIN, QUIET);
}

xDOSCAR::xDOSCAR(const xDOSCAR& b) : xStream(*b.getOFStream(),*b.getOSS()) {free();copy(b);} // copy PUBLIC

xDOSCAR::~xDOSCAR() {xStream::free();free();} //CO20191110 //CO20200404 - xStream integration for logging

void xDOSCAR::free() {  //CO20191110
  m_initialized=false; //CO20200404
  content="";
  vcontent.clear();
  filename="";
  title="";
  spin=0;
  Vol=0.0;POTIM=0.0;
  lattice.clear();
  temperature=0.0;
  RWIGS=FALSE;
  Efermi=0.0;
  //spinF=0.0;  //CO
  spinF=AUROSTD_NAN;
  energy_max=0.0;
  energy_min=0.0;
  number_energies=0;
  number_atoms=0; //CO20191010
  partial=false;  //ME20190614
  denergy=0.0;
  venergy.clear();
  venergy.clear();
  vDOS.clear();viDOS.clear();
  // vDOSs.clear();vDOSp.clear();vDOSd.clear(); OBSOLETE ME20190614 - not used
  isLSCOUPLING = false;  //ME20190620
  lmResolved = false;  //ME20190620
  carstring=""; //CO20191010
  //[CO20200404 - OBSOLETE]ERROR="";
  conduction_band_min.clear();
  conduction_band_min_net=AUROSTD_NAN;
  valence_band_max.clear();
  valence_band_max_net=AUROSTD_NAN;
  Egap.clear();
  Egap_net=AUROSTD_NAN;
  Egap_fit.clear();
  Egap_fit_net=AUROSTD_NAN;
  Egap_type.clear();
  Egap_type_net="";
}

void xDOSCAR::copy(const xDOSCAR& b) { // copy PRIVATE
  xStream::copy(b);  //CO20200404 - xStream integration for logging
  free();
  m_initialized=b.m_initialized; //CO20200404 - xStream integration for logging
  content=b.content;
  vcontent.clear(); for(uint i=0;i<b.vcontent.size();i++) vcontent.push_back(b.vcontent.at(i));  // for aflowlib_libraries.cpp
  filename=b.filename;
  title=b.title;
  spin=b.spin;
  Vol=b.Vol;
  POTIM=b.POTIM;
  lattice=b.lattice;
  temperature=b.temperature;
  RWIGS=b.RWIGS;
  Efermi=b.Efermi;
  spinF=b.spinF;
  energy_max=b.energy_max;
  energy_min=b.energy_min;
  number_energies=b.number_energies;
  number_atoms=b.number_atoms;  //CO20191010
  partial = b.partial;  //ME20190614
  denergy=b.denergy;  
  venergy.clear(); for(uint ienergy=0;ienergy<b.venergy.size();ienergy++) venergy.push_back(b.venergy[ienergy]);  //CO20191004
  venergyEf.clear(); for(uint ienergy=0;ienergy<b.venergyEf.size();ienergy++) venergyEf.push_back(b.venergyEf[ienergy]); //CO20191004
  viDOS.clear();
  for(uint ispin=0;ispin<b.viDOS.size();ispin++){ //CO20191004
    viDOS.push_back(deque<double>(0));
    for(uint ienergy=0;ienergy<b.viDOS[ispin].size();ienergy++) viDOS[ispin].push_back(b.viDOS[ispin][ienergy]);
  }
  vDOS.clear();   //ME20190614
  for(uint iatom=0;iatom<b.vDOS.size();iatom++){ //CO20191004
    vDOS.push_back(deque<deque<deque<double> > >(0));
    for(uint iorbital=0;iorbital<b.vDOS[iatom].size();iorbital++){
      vDOS[iatom].push_back(deque<deque<double> >(0));
      for(uint ispin=0;ispin<b.vDOS[iatom][iorbital].size();ispin++){
        vDOS[iatom][iorbital].push_back(deque<double>(0));
        for(uint ienergy=0;ienergy<b.vDOS[iatom][iorbital][ispin].size();ienergy++) vDOS[iatom][iorbital][ispin].push_back(b.vDOS[iatom][iorbital][ispin][ienergy]);
      }
    }
  }
  //vDOS.clear(); for(uint i=0;i<b.vDOS.size();i++) vDOS.push_back(b.vDOS.at(i)); OBSOLETE ME20190614
  //vDOSs.clear(); for(uint i=0;i<b.vDOSs.size();i++) vDOSs.push_back(b.vDOSs.at(i));  OBSOLETE ME20190614 - not used
  //vDOSs.clear(); for(uint i=0;i<b.vDOSs.size();i++) vDOSs.push_back(b.vDOSs.at(i));  OBSOLETE ME20190614 - not used
  //vDOSp.clear(); for(uint i=0;i<b.vDOSp.size();i++) vDOSp.push_back(b.vDOSp.at(i));  OBSOLETE ME20190614 - not used
  //vDOSd.clear(); for(uint i=0;i<b.vDOSd.size();i++) vDOSd.push_back(b.vDOSd.at(i));  OBSOLETE ME20190614 - not used
  isLSCOUPLING = b.isLSCOUPLING;  //ME20190620
  lmResolved = b.lmResolved;  //ME20190620
  carstring=b.carstring;  //CO20191010
  //[CO20200404 - OBSOLETE]ERROR=b.ERROR;
  conduction_band_min.clear(); for(uint ispin=0;ispin<b.conduction_band_min.size();ispin++){conduction_band_min.push_back(b.conduction_band_min[ispin]);}
  conduction_band_min_net=b.conduction_band_min_net;
  valence_band_max.clear(); for(uint ispin=0;ispin<b.valence_band_max.size();ispin++){valence_band_max.push_back(b.valence_band_max[ispin]);}
  valence_band_max_net=b.valence_band_max_net;
  Egap.clear(); for(uint ispin=0;ispin<b.Egap.size();ispin++){Egap.push_back(b.Egap[ispin]);}
  Egap_net=b.Egap_net;
  Egap_fit.clear(); for(uint ispin=0;ispin<b.Egap_fit.size();ispin++){Egap_fit.push_back(b.Egap_fit[ispin]);}
  Egap_fit_net=b.Egap_fit_net;
  Egap_type.clear(); for(uint ispin=0;ispin<b.Egap_type.size();ispin++){Egap_type.push_back(b.Egap_type[ispin]);}
  Egap_type_net=b.Egap_type_net;
}

const xDOSCAR& xDOSCAR::operator=(const xDOSCAR& b) {  // operator= PUBLIC
  if(this!=&b) {free();copy(b);}
  return *this;
}

void xDOSCAR::clear() {  // clear PRIVATE
  xDOSCAR _temp;
  string filename_aus=filename;
  copy(_temp);
  filename=filename_aus;
}

bool xDOSCAR::GetProperties(const string& stringIN,bool QUIET) {
  stringstream sss; sss.str(stringIN);
  if(filename=="") filename="string";
  return xDOSCAR::GetProperties(sss,QUIET);
}

bool xDOSCAR::GetPropertiesFile(const string& fileIN,bool QUIET) {
  stringstream sss;
  //[CO20210315 - always set for GetPropertiesFile]if(filename=="") 
  filename=fileIN;
  aurostd::efile2stringstream(fileIN,sss);
  return xDOSCAR::GetProperties(sss,QUIET);
}

bool xDOSCAR::GetPropertiesUrlFile(const string& url,const string& file,bool QUIET) {
  //[CO20200502 - OBSOLETE]string tmpfile=XHOST.tmpfs+"/_aflow_"+XHOST.user+".pid"+XHOST.ostrPID.str()+".a"+string(AFLOW_VERSION)+".rnd"+aurostd::utype2string(uint((double) std::floor((double)100000*aurostd::ran0())))+".u"+aurostd::utype2string(uint((double) aurostd::get_useconds()))+"_"+file;
  string tmpfile=aurostd::TmpFileCreate("xDOSCAR_GetProperties"); //CO20200502 - threadID
  aurostd::url2file(url+"/"+file,tmpfile,!QUIET);
  bool out=GetPropertiesFile(tmpfile,QUIET);
  filename="url="+url;  //CO20210315
  aurostd::RemoveFile(tmpfile);
  return out;
}

bool xDOSCAR::GetProperties(const stringstream& stringstreamIN,bool QUIET) {
  bool LDEBUG=(FALSE || XHOST.DEBUG || !QUIET);
  bool force_exit=XHOST.POSTPROCESS; //SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only
  stringstream message;

  bool ERROR_flag=FALSE;
  long double seconds=aurostd::get_seconds();
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN (" << time_delay(seconds) << ")" << endl;
  clear(); // so it does not mess up vector/deque
  content=stringstreamIN.str();
  vcontent.clear();
  vector<string> vline,tokens;
  aurostd::string2vectorstring(content,vcontent);
  string line;
  uint ndos = 1;  //ME20190614
  if(filename=="") filename="stringstream";
  //ME20190812 - Add checks for broken DOSCARs
  if (vcontent.size() < 7) {
    //[CO20200404 - OBSOLETE]string function = "xDOSCAR::GetProperties()";
    message << "Broken DOSCAR: no content";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _FILE_ERROR_);
  }
  for(uint iline = 0; iline < 7;iline++) { //ME20190614 - Read header
    aurostd::string2tokens(vcontent.at(iline),tokens);
    // cerr << "iline=" << iline << "  " << vcontent.at(iline) << " tokens.size()=" << tokens.size() << endl;
    // WRONG if(iline==0 && tokens.size()==4)  spin=aurostd::string2utype<int>(tokens.at(tokens.size()-1))-1;
    //ME20190614 START
    if (iline == 0 && tokens.size() >= 4) {
      number_atoms = aurostd::string2utype<uint>(tokens[1]);
      partial = aurostd::string2utype<bool>(tokens[2]);
      if (partial) ndos += number_atoms;
    }
    //ME20190614 END
    if(iline==1 && tokens.size()>=5) {
      uint i=0;
      Vol=aurostd::string2utype<double>(tokens.at(i++));
      lattice(1)=aurostd::string2utype<double>(tokens.at(i++));
      lattice(2)=aurostd::string2utype<double>(tokens.at(i++));
      lattice(3)=aurostd::string2utype<double>(tokens.at(i++));
      POTIM=aurostd::string2utype<double>(tokens.at(i++));
    }
    if(iline==2) temperature=aurostd::string2utype<double>(vcontent.at(iline));
    if(iline==3) carstring = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vcontent.at(iline));  //ME20190620 - what kind of DOSCAR?
    if(iline==4) title=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vcontent.at(iline));  //ME20190508 - clean title
    if(iline==5) {
      // cerr << "iline=" << iline << "  " << vcontent.at(iline) << " tokens.size()=" << tokens.size() << endl;
      energy_max=aurostd::string2utype<double>(tokens.at(0));
      energy_min=aurostd::string2utype<double>(tokens.at(1));
      number_energies=aurostd::string2utype<uint>(tokens.at(2));
      Efermi=aurostd::string2utype<double>(tokens.at(3));
    }
    if(iline==6 && tokens.size()==4) {spin=0;RWIGS=TRUE;}
    if(iline==6 && tokens.size()==7) {spin=1;RWIGS=TRUE;}
    //ME20190614 START
    if (iline == 6) {
      if (tokens.size() == 3) spin = 0;
      if (tokens.size() == 5) spin = 1;
    }
    //ME20190614 END
  }
  // Done reading header
  //ME20190614 START
  //ME20200305 - Check if DOSCAR is broken
  uint number_lines = number_energies + 6;
  if (partial) number_lines += number_atoms * (number_energies + 1);
  if (vcontent.size() < number_lines) {
    message << "Broken DOSCAR - not enough lines.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
  }
  uint norbitals = 0;
  int d = 0, e = 0;
  double dos = 0.0;
  if (partial) {
    aurostd::string2tokens(vcontent.at(number_energies + 7), tokens, " ");
    int ncols = (int) tokens.size() - 1;  // Don't count the energy column
    if (carstring == "PHON") {  // APL DOSCAR
      norbitals = ncols;
      isLSCOUPLING = false;
      lmResolved = true;
    } else {  // No special case, so VASP DOSCAR
      // Determine whether the DOSCAR is lm-resolved and/or has spin-orbit coupling (SOC)
      if (ncols == 16) {
        // If the DOSCAR has 16 columns, the variables cannot be properly resolved.
        // It could either lm-resolved with f-electrons and without spin polarization,
        // or it is not lm-resolved with f-electrons and SOC. We need another input file
        // to resolve this. The former case is VERY unlikely since calculations with
        // f-electrons typically require spin polarization. So, SOC is the default.
        // First, try and find filename extensions
        vector<string> vstr;
        string ext = "";
        if (aurostd::substring2bool(filename, "DOSCAR.")) { 
          aurostd::string2tokens(filename, vstr, ".");
          for (uint i = 1; i < vstr.size(); i++) ext += "." + vstr[i];
        } else if (aurostd::substring2bool(filename, "DOSCAR_")) {
          aurostd::string2tokens(filename, vstr, "_");
          for (uint i = 1; i < vstr.size(); i++) ext += "_" + vstr[i];
        }
        vstr.clear();
        // Try INCAR first because it's the smallest file
        if (aurostd::EFileExist("INCAR" + ext) || aurostd::FileExist("INCAR" + ext)) {
          aurostd::efile2vectorstring("INCAR" + ext, vstr);
          aurostd::RemoveComments(vstr);
          for (uint i = 0; i < vstr.size(); i++) {
            if (aurostd::substring2bool(vstr[i], "LSORBIT")) {
              vector<string> tokens;
              aurostd::string2tokens(vstr[i], tokens, " =");
              tokens[1] = aurostd::toupper(tokens[1]);
              if ((tokens[1][0] == 'T') || (aurostd::substring2bool(tokens[1], "TRUE"))) {
                isLSCOUPLING = true;
              } else {
                isLSCOUPLING = false;
              }
              break;
            }
          }
          // No INCAR found. Try vasprun.xml
        } else if (aurostd::EFileExist("vasprun.xml" + ext) || aurostd::FileExist("vasprun.xml" + ext)) {
          aurostd::efile2vectorstring("vasprun.xml" + ext, vstr);
          for (uint i = 0; i < vstr.size(); i++) {
            if (aurostd::substring2bool(vstr[i], "LSORBIT")) {
              isLSCOUPLING = !(aurostd::substring2bool(vstr[i], "F"));  // Only contains "F" if LSORBIT is false
              break;
            }
          }
          // At last, try OUTCAR
        } else if (aurostd::EFileExist("OUTCAR" + ext) || aurostd::FileExist("OUTCAR" + ext)) {
          aurostd::efile2vectorstring("vasprun.xml" + ext, vstr);
          for (uint i = 0; i < vstr.size(); i++) {
            if (aurostd::substring2bool(vstr[i], "LSORBIT")) {
              isLSCOUPLING = !(aurostd::substring2bool(vstr[i], "F"));  // Only contains "F" if LSORBIT is false
              break;
            }
          }
          // Nothing found, assume SOC (more likely case)
        } else {
          string message = "Could not determine whether the DOSCAR is lm-resolved"
            " or contains spin-orbit coupling. AFLOW will assume that the"
            " DOSCAR is lm-resolved. If this is not the case, please put an"
            " INCAR" + ext + ", a vasprun.xml " + ext + ", or an"
            " OUTCAR" +  ext + " file into the working directory and try again.";
          pflow::logger(__AFLOW_FILE__, "xDOSCAR::GetProperties()", message, std::cerr, _LOGGER_WARNING_,QUIET);
          isLSCOUPLING = true;
        }
        lmResolved = !(isLSCOUPLING);  // With 16 columns, it cannot be both
      } else {
        // Otherwise, the number of columns per spin smaller than 10 without
        // SOC, except for lm-resolved DOS with f-electrons (32 columns)
        isLSCOUPLING = ((ncols/(spin + 1) > 10) && (ncols != 32));
      }
      if (isLSCOUPLING) norbitals = ncols/4;
      else norbitals = (ncols)/(spin + 1);
      // Since VASP always prints s, p, and d orbitals, lm-resolved
      // DOSCARS contain at least 9 oribtals
      lmResolved = (norbitals > 8);
    }
  }
  if (isLSCOUPLING) {  //ME20190620 - LSCOUPLING has four spin channels
    vDOS.assign(ndos, deque<deque<deque<double> > >(norbitals + 1, deque<deque<double> >(4, deque<double>(number_energies, 0.0))));
  } else {
    vDOS.assign(ndos, deque<deque<deque<double> > >(norbitals + 1, deque<deque<double> >(spin + 1, deque<double>(number_energies, 0.0))));
  }
  venergy.resize(number_energies);
  viDOS.assign(spin + 1, deque<double>(number_energies, 0.0));
  // Read data
  for (uint iline = 6; iline < vcontent.size(); iline++) {
    if (iline == (d + 1) * number_energies + 6 + d) {
      d++;
      //ME20190810 - Safeguard against DOSCARs with additional lines
      if (d == (int) ndos) {
        message << "DOSCAR contains more lines than the header suggests." << endl;
        message << "xDOSCAR object may not be properly populated.";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_oss, _LOGGER_WARNING_,QUIET);
        break;
      }
      e = 0;
    } else {
      aurostd::string2tokens(vcontent.at(iline),tokens);
      if (d == 0) {
        venergy[e] = aurostd::string2utype<double>(tokens[0]);
        for (uint i = 0; i < (spin + 1); i++) {
          vDOS[0][0][i][e] = aurostd::string2utype<double>(tokens[i + 1]);
          viDOS[i][e] = aurostd::string2utype<double>(tokens[i + spin + 2]);
        }
        e++;
      } else if (isLSCOUPLING) {  //ME20190620 - LSCOUPLING has four spin channels
        for (uint i = 0; i < 4 * norbitals; i++) {
          dos = aurostd::string2utype<double>(tokens[i + 1]);
          vDOS[d][i/4 + 1][i % 4][e] = dos;
          vDOS[d][0][i % 4][e] += dos;
          vDOS[0][i/4 + 1][i % 4][e] += dos;
        }
        e++;
      } else {
        for (uint i = 0; i < (spin + 1) * norbitals; i++) {
          dos = aurostd::string2utype<double>(tokens[i + 1]);
          vDOS[d][i/(spin + 1) + 1][i % (spin + 1)][e] = dos;
          vDOS[d][0][i % (spin + 1)][e] += dos;
          vDOS[0][i/(spin + 1) + 1][i % (spin + 1)][e] += dos;
        }
        e++;
      }
      //ME20190614 - The procedure below makes too many assumptions about how the
      // DOSCAR file is structured and fails for f-elements.
      //[OBSOLETE]    if(iline>=6 && iline < number_energies + 6 ) {
      //[OBSOLETE]      aurostd::string2tokens(vcontent.at(iline),tokens);
      //[OBSOLETE]      uint i=0;
      //[OBSOLETE]      if(tokens.size()==3 || tokens.size()==5) {
      //[OBSOLETE]	RWIGS=FALSE;
      //[OBSOLETE]	if(tokens.size()==3) spin=0;
      //[OBSOLETE]	if(tokens.size()==5) spin=1;
      //[OBSOLETE]	if (venergy.size() <= number_energies) venergy.push_back(aurostd::string2utype<double>(tokens.at(i++)));  //ME20190507
      //[OBSOLETE]	//vDOS
      //[OBSOLETE]	deque<double> entryDOS;
      //[OBSOLETE]	entryDOS.push_back(aurostd::string2utype<double>(tokens.at(i++)));
      //[OBSOLETE]	if(tokens.size()==5) entryDOS.push_back(aurostd::string2utype<double>(tokens.at(i++)));
      //[OBSOLETE]	vDOS.push_back(entryDOS);
      //[OBSOLETE]	//viDOS
      //[OBSOLETE]	deque<double> entryiDOS;
      //[OBSOLETE]	entryiDOS.push_back(aurostd::string2utype<double>(tokens.at(i++)));
      //[OBSOLETE]	if(tokens.size()==5) entryiDOS.push_back(aurostd::string2utype<double>(tokens.at(i++)));
      //[OBSOLETE]	viDOS.push_back(entryiDOS);
      //[OBSOLETE]      }
      //[OBSOLETE]      if(tokens.size()==4 || tokens.size()==7) {
      //[OBSOLETE]	RWIGS=TRUE;
      //[OBSOLETE]	if(tokens.size()==4) spin=0;
      //[OBSOLETE]	if(tokens.size()==7) spin=1;
      //[OBSOLETE]	venergy.push_back(aurostd::string2utype<double>(tokens.at(i++)));
      //[OBSOLETE]	//vDOS s p d
      //[OBSOLETE]	deque<double> entryDOSs,entryDOSp,entryDOSd;
      //[OBSOLETE]	entryDOSs.push_back(aurostd::string2utype<double>(tokens.at(i++)));
      //[OBSOLETE]	entryDOSp.push_back(aurostd::string2utype<double>(tokens.at(i++)));
      //[OBSOLETE]	entryDOSd.push_back(aurostd::string2utype<double>(tokens.at(i++)));
      //[OBSOLETE]	if(tokens.size()==7) entryDOSs.push_back(aurostd::string2utype<double>(tokens.at(i++)));
      //[OBSOLETE]	if(tokens.size()==7) entryDOSp.push_back(aurostd::string2utype<double>(tokens.at(i++)));
      //[OBSOLETE]	if(tokens.size()==7) entryDOSd.push_back(aurostd::string2utype<double>(tokens.at(i++)));
      //[OBSOLETE]	vDOSs.push_back(entryDOSs);
      //[OBSOLETE]	vDOSp.push_back(entryDOSp);
      //[OBSOLETE]	vDOSd.push_back(entryDOSd);
      //[OBSOLETE]      }
      //[OBSOLETE]    }
    }
    //ME20190614 END
  }
  //ME20190812 - Safeguard against broken DOSCARs
  if ((d + 1 < (int) ndos) || (e < (int) number_energies)) {  //ME20191010: needs to be d + 1
    //[CO20200404 - OBSOLETE]string function = "xDOSCAR::GetProperties()";
    message << "Broken DOSCAR: not enough lines";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _FILE_ERROR_);
  }
  // fix denergy
  denergy=venergy.at(1)-venergy.at(0);
  venergyEf.clear();
  for(uint i=0;i<venergy.size();i++)
    venergyEf.push_back(venergy.at(i)-Efermi);

  // ----------------------------------------------------------------------
  // spin polarization at FERMI level

  double Fup=0.0,Fdown=0.0,Minup=0.0,Mindown=0.0,Maxup=0.0,Maxdown=0.0; // CAMILOFIX
  bool Fermifound=FALSE,firstenter=TRUE;
  double zeroTol=1e-8;
  if(!spin) {
    spinF=0.0;
  } else {
    for(uint i=6;i<vcontent.size();i++) {
      double energytp;
      aurostd::string2tokens(vcontent.at(i),tokens," ");
      if(tokens.size()>0 && tokens.size()<6) {
        energytp=aurostd::string2utype<double>(tokens.at(0));
        if(Efermi-energytp>0) {
          Minup=aurostd::string2utype<double>(tokens.at(1));Mindown=aurostd::string2utype<double>(tokens.at(2));
        } else if(firstenter) {
          Maxup=aurostd::string2utype<double>(tokens.at(1));Maxdown=aurostd::string2utype<double>(tokens.at(2));
          firstenter=FALSE;
          Fermifound=TRUE;
        }
        if(Fermifound) {
          Fup=Minup+Maxup;Fdown=Mindown+Maxdown;
          break;
        }
      }
    }
    if((Fup+Fdown)<zeroTol) {
      spinF=0.0;
    } else {
      spinF=aurostd::abs((Fup-Fdown)/(Fup+Fdown)); // otherwise AFLOW_NAN
    }
  }

  if(0) cout << " spinF = " << ((spinF!=AUROSTD_NAN)?aurostd::utype2string(spinF,5):"unavailable") << endl;

  // ----------------------------------------------------------------------
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " title=" << title << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " spin=" << spin << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " Vol=" << Vol << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " lattice=" << lattice << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " POTIM=" << POTIM << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " temperature=" << temperature << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " RWIGS=" << RWIGS << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " Efermi=" << Efermi << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " spinF=" << spinF << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " number_energies=" << number_energies << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " energy_max=" << energy_max << " energy_min=" << energy_min << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " denergy=" << denergy << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " venergy.size()=" << venergy.size() << " venergyEf.size()=" << venergyEf.size() << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vDOS.size()=" << vDOS.size() << ", " << vDOS[0].size() << ", " << vDOS[0][0].size() << ", " << vDOS[0][0][0].size() << std::endl;  //ME20190614 - new vDOS format
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " viDOS.size()=" << viDOS.size() << ", " << viDOS[0].size() << ", " << std::endl;  //ME20190614 - new viDOS format
  //if(LDEBUG) cerr << __AFLOW_FUNC__ << " vDOS.size()=" << vDOS.size() << " vDOS.at(max).size()=" << vDOS.at(vDOS.size()-1).size() << endl;  OBSOLETE ME20190614
  //if(LDEBUG) cerr << __AFLOW_FUNC__ << " viDOS.size()=" << viDOS.size() << " viDOS.at(max).size()=" << viDOS.at(viDOS.size()-1).size() << endl;  OBSOLETE ME20190614
  //if(LDEBUG) cerr << __AFLOW_FUNC__ << " vDOSs.size()=" << vDOSs.size() << endl;  OBSOLETE ME20190614
  // if(LDEBUG) cerr << __AFLOW_FUNC__ << " vDOSs.at(max).size()=" << vDOSs.at(vDOSs.size()-1).size() << endl;
  //if(LDEBUG) cerr << __AFLOW_FUNC__ << " vDOSp.size()=" << vDOSp.size() << endl;  OBSOLETE ME20190614
  // if(LDEBUG) cerr << __AFLOW_FUNC__ << " vDOSp.at(max).size()=" << vDOSp.at(vDOSp.size()-1).size() << endl;
  //if(LDEBUG) cerr << __AFLOW_FUNC__ << " vDOSd.size()=" << vDOSd.size() << endl;  OBSOLETE ME20190614
  // if(LDEBUG) cerr << __AFLOW_FUNC__ << " vDOSd.at(max).size()=" << vDOSd.at(vDOSd.size()-1).size() << endl;
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " END (" << time_delay(seconds) << ")" << endl;  
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  //[CO20200404 - OBSOLETE]if(ERROR_flag && !QUIET) cerr << "WARNING - xDOSCAR::GetProperties: ERROR_flag set in xDOSCAR" << endl;
  if(ERROR_flag){
    message << "ERROR_flag set in xDOSCAR";
    if(force_exit){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);}
    else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_,QUIET);}
    return FALSE;
  }
  return TRUE;
}

//CO20191217 - copies everything from spin channel 1 to spin channel 2
void xDOSCAR::convertSpinOFF2ON() { //CO20191217
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

  uint iatom=0,iorbital=0;

  //check that it is truly SPIN-OFF
  if(viDOS.size()!=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"viDOS.size()!=1",_INPUT_ERROR_);}; //no conversion needed
  for(iatom=0;iatom<vDOS.size();iatom++){
    for(iorbital=0;iorbital<vDOS[iatom].size();iorbital++){
      if(vDOS[iatom][iorbital].size()!=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"vDOS[iatom][iorbital].size()!=1",_INPUT_ERROR_);}; //no conversion needed
    }
  }
  if(conduction_band_min.size()!=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"conduction_band_min.size()!=1",_INPUT_ERROR_);}; //no conversion needed
  if(valence_band_max.size()!=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"valence_band_max.size()!=1",_INPUT_ERROR_);}; //no conversion needed
  if(Egap.size()!=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Egap!=1",_INPUT_ERROR_);}; //no conversion needed
  if(Egap_fit.size()!=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Egap_fit!=1",_INPUT_ERROR_);}; //no conversion needed
  if(Egap_type.size()!=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Egap_type!=1",_INPUT_ERROR_);}; //no conversion needed
  if(spin!=0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"spin!=0",_INPUT_ERROR_);}; //no conversion needed

  //copy everything over
  viDOS.push_back(viDOS.back());
  for(iatom=0;iatom<vDOS.size();iatom++){
    for(iorbital=0;iorbital<vDOS[iatom].size();iorbital++){
      vDOS[iatom][iorbital].push_back(vDOS[iatom][iorbital].back());
    }
  }
  conduction_band_min.push_back(conduction_band_min.back());
  valence_band_max.push_back(valence_band_max.back());
  Egap.push_back(Egap.back());
  Egap_fit.push_back(Egap_fit.back());
  Egap_type.push_back(Egap_type.back());
  spin=1;
}

void xDOSCAR::addAtomChannel(){  //CO20211124
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

  string ERROR_out="";
  if(!checkDOS(ERROR_out)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,ERROR_out,_INDEX_BOUNDS_);} //no conversion needed

  uint atoms_size=vDOS.size();if(atoms_size==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"vDOS.size()==0",_INDEX_BOUNDS_);} //no conversion needed
  uint orbital_size=vDOS.front().size();
  uint spin_size=0;if(orbital_size>0){spin_size=vDOS.front().front().size();}
  uint energy_size=0;if(orbital_size>0 && spin_size>0){energy_size=vDOS.front().front().front().size();}
  vDOS.push_back(deque<deque<deque<double> > >(0));
  uint iorbital=0,ispin=0;
  if(orbital_size>0){vDOS.back().resize(orbital_size);}
  if(spin_size>0){
    for(iorbital=0;iorbital<orbital_size;iorbital++){
      vDOS.back()[iorbital].resize(spin_size);
      if(energy_size>0){
        for(ispin=0;ispin<spin_size;ispin++){vDOS.back()[iorbital][ispin].resize(energy_size,0.0);}
      }
    }
  }
  number_atoms+=1;
}

void xDOSCAR::addOrbitalChannel(){  //CO20211124
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

  string ERROR_out="";
  if(!checkDOS(ERROR_out)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,ERROR_out,_INDEX_BOUNDS_);} //no conversion needed

  uint atoms_size=vDOS.size();if(atoms_size==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"vDOS.size()==0",_INDEX_BOUNDS_);} //no conversion needed
  uint orbital_size=vDOS.front().size();
  uint spin_size=0;if(orbital_size>0){spin_size=vDOS.front().front().size();}
  uint energy_size=0;if(orbital_size>0 && spin_size>0){energy_size=vDOS.front().front().front().size();}
  uint iatom=0,ispin=0;
  for(iatom=0;iatom<atoms_size;iatom++){
    vDOS[iatom].push_back(deque<deque<double> >(0));
    if(spin_size>0){vDOS[iatom].back().resize(spin_size);}
    if(energy_size>0){
      for(ispin=0;ispin<spin_size;ispin++){vDOS[iatom].back()[ispin].resize(energy_size,0.0);}
    }
  }
}

void xDOSCAR::resetVDOS(){  //CO20211124
  uint iatom=0,iorbital=0,ispin=0;
  for(iatom=0;iatom<vDOS.size();iatom++){
    for(iorbital=0;iorbital<vDOS[iatom].size();iorbital++){
      for(ispin=0;ispin<vDOS[iatom][iorbital].size();ispin++){
        std::fill(vDOS[iatom][iorbital][ispin].begin(),vDOS[iatom][iorbital][ispin].end(),0.0);
      }
    }
  }
}

bool xDOSCAR::checkDOS(string& ERROR_out) const { //CO20191110
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream message;

  uint IENERGY=venergy.size();
  if(IENERGY==0){
    message << "GetProperties() failed. IENERGY==0 (venergy.size()==0)";
    ERROR_out=message.str();return false;
  }
  uint IATOM=vDOS.size();
  if(IATOM==0){
    message << "GetProperties() failed. IATOM==0 (vDOS.size()==0)";
    ERROR_out=message.str();return false;
  }
  uint IORBITAL=vDOS.front().size();
  if(IORBITAL==0){
    message << "GetProperties() failed. IORBITAL==0 (vDOS.front().size()==0)";
    ERROR_out=message.str();return false;
  }
  uint ISPIN=vDOS.front().front().size();
  if(ISPIN==0){
    message << "GetProperties() failed. ISPIN==0 (vDOS.front().front().size()==0)";
    ERROR_out=message.str();return false;
  }
  for(uint iatom=0;iatom<IATOM;iatom++){
    if(vDOS[iatom].size()!=IORBITAL){
      message << "GetProperties() failed. (vDOS[iatom=" << iatom << "].size()!=vDOS.front().size())";
      ERROR_out=message.str();return false;
    }
    for(uint iorbital=0;iorbital<IORBITAL;iorbital++){
      if(vDOS[iatom][iorbital].size()!=ISPIN){
        message << "GetProperties() failed. (vDOS[iatom=" << iatom << "][iorbital=" << iorbital << "].size()!=vDOS.front().front().size())";
        ERROR_out=message.str();return false;
      }
      for(uint ispin=0;ispin<ISPIN;ispin++){
        if(vDOS[iatom][iorbital][ispin].size()!=IENERGY){
          message << "GetProperties() failed. (vDOS[iatom=" << iatom << "][iorbital=" << iorbital << "][ispin=" << ispin << "].size()!=venergy.size())";
          ERROR_out=message.str();return false;
        }
      }
    }
  }
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " DOSCAR properties retrieved" << endl;}
  return true;
}

bool xDOSCAR::GetBandGap(double EFERMI,double efermi_tol,double energy_tol,double occ_tol) { //CO20191004
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream message;
  bool force_exit=XHOST.POSTPROCESS; //SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only

  if((content == "") || (vcontent.size() == 0)) {
    //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " xDOSCAR needs to be loaded before. \n"
    //[CO20200404 - OBSOLETE]  "        GetProperties(const stringstream&); \n"
    //[CO20200404 - OBSOLETE]  "        GetProperties(const string&); \n"
    //[CO20200404 - OBSOLETE]  "        GetPropertiesFile(const string&); \n";
    //[CO20200404 - OBSOLETE]return FALSE;
    message << "xDOSCAR needs to be loaded before." << endl;
    message << "GetProperties(const stringstream&);" << endl;
    message << "GetProperties(const string&);" << endl;
    message << "GetPropertiesFile(const string&);" << endl;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _INPUT_MISSING_);
  }
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " DOSCAR content found" << endl;}

  string error_string="";
  if(!checkDOS(error_string)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,error_string,_INPUT_ERROR_);};  //quick check if GetProperties() failed

  //this should all work now that we checkDOS()
  uint IENERGY=venergy.size();
  uint IATOM=vDOS.size();
  uint IORBITAL=vDOS.front().size();
  uint ISPIN=vDOS.front().front().size();

  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " IENERGY=" << IENERGY << endl;
    cerr << __AFLOW_FUNC__ << " IATOM=" << IATOM << endl;
    cerr << __AFLOW_FUNC__ << " IORBITAL=" << IORBITAL << endl;
    cerr << __AFLOW_FUNC__ << " ISPIN=" << ISPIN << endl;
  }

  // GET FERMI LEVEL
  if(EFERMI==AUROSTD_NAN){
    //VERY IMPORTANT - CO20171009 from discussion with SC
    //we strongly prefer to use Efermi from DOSCAR.static, not DOSCAR.bands
    //bands is not self-consistent (ICHARG=11), only used to determine energy states
    //since it is not self-consistent, it does not fill in states correctly
    //static is the best way to go!

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " Using Efermi from CURRENT DOSCAR (recommended to use Efermi from DOSCAR.static)" << endl;}
    EFERMI=Efermi;
  }

  //error in E-fermi, unless otherwise provided, is the difference between the input E-fermi (STATIC) and the one found in this DOSCAR (BANDS)
  //can be significant!
  //think of this as largest possible shift of energies in BANDS relative to real energy levels
  if(efermi_tol==AUROSTD_NAN){
    efermi_tol=Efermi-EFERMI;                             //generlly, E-fermi BANDS > E-fermi STATIC
    if(std::signbit(efermi_tol)){efermi_tol=energy_tol;}  //just in case
  }
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " tol(E-fermi)=" << efermi_tol << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " tol(occ)=" << occ_tol << endl;}

  valence_band_max.resize(ISPIN);
  conduction_band_min.resize(ISPIN);
  for(uint ispin=0;ispin<ISPIN;ispin++){
    valence_band_max[ispin]=AUROSTD_MAX_DOUBLE;
    conduction_band_min[ispin]=AUROSTD_MAX_DOUBLE;
  }

  bool found_EFERMI=false;
  bool metal_found=false;
  for(uint ispin=0;ispin<ISPIN;ispin++){
    found_EFERMI=false;
    for(uint ienergy=0;ienergy<vDOS.front().front()[ispin].size();ienergy++){
      if(venergy[ienergy]>=EFERMI){
        if(valence_band_max[ispin]==AUROSTD_MAX_DOUBLE){
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " DOS[energy=" << venergy[ienergy] << "]=" << vDOS.front().front()[ispin][ienergy] << " ?< " << occ_tol << endl;}
          if(found_EFERMI && (venergy[ienergy]-EFERMI)>efermi_tol){  //found metal for this spin
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " FOUND metal: ispin=" << ispin << endl;}
            metal_found=true;
            break;
          }
          if(vDOS.front().front()[ispin][ienergy]<occ_tol){
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " FOUND valence_band_max[ispin=" << ispin << "]=" << venergy[ienergy] << endl;}
            valence_band_max[ispin]=venergy[ienergy];
          }
        }else{
          if(vDOS.front().front()[ispin][ienergy]>=occ_tol){
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " FOUND conduction_band_min[ispin=" << ispin << "]=" << venergy[ienergy] << endl;}
            conduction_band_min[ispin]=venergy[ienergy];
            break;
          }
        }
        found_EFERMI=true;
      }
    }
    if(valence_band_max[ispin]!=AUROSTD_MAX_DOUBLE && conduction_band_min[ispin]==AUROSTD_MAX_DOUBLE){
      //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " conduction band maximum not found [ispin="+aurostd::utype2string(ispin)+"] \n";
      message << "Conduction band maximum not found [ispin=" << ispin << "]";
      if(force_exit){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);}
      else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);return false;}
    }
  }

  if(!found_EFERMI){
    //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " E-fermi not found \n";
    //[CO20200404 - OBSOLETE]return false;
    message << " E-fermi not found";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
  }

  if(LDEBUG){
    for(uint ispin=0;ispin<ISPIN;ispin++){
      cerr << __AFLOW_FUNC__ << " valence_band_max[ispin=" << ispin << "]=" << valence_band_max[ispin] << endl;
      cerr << __AFLOW_FUNC__ << " conduction_band_min[ispin=" << ispin << "]=" << conduction_band_min[ispin] << endl;
    }
    cerr << __AFLOW_FUNC__ << " metal_found=" << metal_found << endl;
  }

  //[CO20191004]double _METALGAP = -AUROSTD_NAN, _METALEDGE = -1.0;
  Egap.resize(ISPIN);
  Egap_fit.resize(ISPIN);
  Egap_type.resize(ISPIN);
  for(uint ispin=0;ispin<ISPIN;ispin++){
    if(valence_band_max[ispin]!=AUROSTD_MAX_DOUBLE && conduction_band_min[ispin]!=AUROSTD_MAX_DOUBLE){
      if(valence_band_max[ispin]>conduction_band_min[ispin]){
        //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " Negative band gap found in ispin="+aurostd::utype2string(ispin)+", not expected \n";
        message << "Negative band gap found in ispin=" << ispin << ", not expected";
        if(force_exit){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);}
        else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);return false;}
      }
      Egap[ispin]=conduction_band_min[ispin]-valence_band_max[ispin];
      Egap_fit[ispin]=1.348 * Egap[ispin] + 0.913;
      Egap_type[ispin]="insulator";
      if(Egap[ispin]<energy_tol){Egap_type[ispin]+="_zero_gap";}
    }else if(valence_band_max[ispin]==AUROSTD_MAX_DOUBLE && conduction_band_min[ispin]==AUROSTD_MAX_DOUBLE){
      valence_band_max[ispin]=_METALEDGE_;
      conduction_band_min[ispin]=_METALEDGE_;
      Egap[ispin]=_METALGAP_;
      Egap_fit[ispin]=_METALGAP_;
      Egap_type[ispin]="metal";
    }else{
      //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " valence_band_max/conduction_band_min not found [ispin="+aurostd::utype2string(ispin)+"] \n";
      message << "Valence_band_max/conduction_band_min not found [ispin=" << ispin << "]";
      if(force_exit){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);}
      else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);return false;}
    }
  }
  if(ISPIN==1){
    Egap_net=Egap[0];
    Egap_fit_net=Egap_fit[0];
    Egap_type_net=Egap_type[0];
  }else{
    if(metal_found){
      Egap_net=_METALGAP_;
      Egap_fit_net=_METALGAP_;
      if(aurostd::substring2bool(Egap_type[0],"metal") && aurostd::substring2bool(Egap_type[1],"metal")){Egap_type_net="metal";}
      else{Egap_type_net="half-metal";}
    }else{
      double cbm=aurostd::min(conduction_band_min[0],conduction_band_min[1]);
      double vbm=aurostd::max(valence_band_max[0],valence_band_max[1]);
      if(vbm>cbm){
        //[CO20200404 - OBSOLETE]ERROR = __AFLOW_FUNC__ + " Negative band gap found in net, not expected \n";
        message << "Negative band gap found in net, not expected";
        if(force_exit){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);}
        else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);return false;}
      }
      Egap_net=cbm-vbm;
      Egap_fit_net=1.348 * Egap_net + 0.913;
      Egap_type_net="insulator";
      if(Egap_net<energy_tol){Egap_type_net+="_zero_gap";}
    }
  }

  if(LDEBUG){
    for(uint ispin=0;ispin<ISPIN;ispin++){
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

deque<deque<deque<deque<double> > > > xDOSCAR::GetVDOSSpecies(const xstructure& xstr) const {return GetVDOSSpecies(xstr.num_each_type);} //CO20191004
deque<deque<deque<deque<double> > > > xDOSCAR::GetVDOSSpecies(deque<int> num_each_type) const { //CO20191004
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream message;

  if((content == "") || (vcontent.size() == 0)) {
    message << "xDOSCAR needs to be loaded before." << endl;
    message << "       GetProperties(const stringstream&);" << endl;
    message << "       GetProperties(const string&);" << endl;
    message << "       GetPropertiesFile(const string&);";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_);
  }
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " DOSCAR content found" << endl;}

  string error_string=""; //keep this function const for plotter
  if(!checkDOS(error_string)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,error_string,_INPUT_ERROR_);};  //quick check if GetProperties() failed

  //this should all work now that we checkDOS()
  uint IENERGY=venergy.size();
  uint IATOM=vDOS.size();
  uint IORBITAL=vDOS.front().size();
  uint ISPIN=vDOS.front().front().size();

  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " IENERGY=" << IENERGY << endl;
    cerr << __AFLOW_FUNC__ << " IATOM=" << IATOM << endl;
    cerr << __AFLOW_FUNC__ << " IORBITAL=" << IORBITAL << endl;
    cerr << __AFLOW_FUNC__ << " ISPIN=" << ISPIN << endl;
  }

  //check that num_each_type corresponds to vDOS
  uint atoms_total=0;
  for(uint iatom=0;iatom<num_each_type.size();iatom++){atoms_total+=num_each_type[iatom];}
  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " atoms_total=" << atoms_total << endl;
    cerr << __AFLOW_FUNC__ << " vDOS.size()=" << vDOS.size() << endl;
  }

  if(atoms_total+1!=IATOM){ //total column
    message << "Input xstructure and DOS mismatch: atoms_total+1!=vDOS.size()" << endl;
    message << "For POCC runs, check that all ARUN.POCC's have the same num_each_type" << endl;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INDEX_MISMATCH_);
  }

  //create vDOS_species with appropriate dimensions
  deque<deque<deque<deque<double> > > > vDOS_species;
  vDOS_species.resize(num_each_type.size()+1);  //+1 for total
  for(uint ispecies=0;ispecies<vDOS_species.size();ispecies++){
    vDOS_species[ispecies].resize(IORBITAL);
    for(uint iorbital=0;iorbital<IORBITAL;iorbital++){
      vDOS_species[ispecies][iorbital].resize(ISPIN);
      for(uint ispin=0;ispin<ISPIN;ispin++){
        vDOS_species[ispecies][iorbital][ispin].assign(IENERGY,0.0);
      }
    }
  }

  uint iatom=0;
  double dos=0; //ME trick
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " summing vDOS_species[ispecies] (individuals)" << endl;}
  for(uint ispecies=0;ispecies<num_each_type.size();ispecies++){
    for(uint i=0;i<(uint)num_each_type[ispecies];i++){
      for(uint iorbital=1;iorbital<IORBITAL;iorbital++){  //iorbital==0 is total (VASP), skip that for now
        for(uint ispin=0;ispin<ISPIN;ispin++){
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " ispecies=" << ispecies << ", iorbital=" << iorbital << ", ispin=" << ispin << endl;}
          for(uint ienergy=0;ienergy<IENERGY;ienergy++){
            dos=vDOS[iatom+1][iorbital][ispin][ienergy];  //iatom==0 is total (VASP), skip that for now
            vDOS_species[ispecies+1][iorbital][ispin][ienergy]+=dos;  //ispecies==0 is total (VASP), skip for now
            vDOS_species[ispecies+1][0][ispin][ienergy]+=dos;
          }
        }
      }
      iatom++;
    }
  }
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " getting vDOS_species[iorbital] (VASP totals)" << endl;}
  for(uint iorbital=1;iorbital<IORBITAL;iorbital++){  //iorbital==0 is total (VASP), skip for now
    for(uint ispin=0;ispin<ISPIN;ispin++){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " iorbital=" << iorbital << ", ispin=" << ispin << endl;}
      for(uint ienergy=0;ienergy<IENERGY;ienergy++){
        dos=vDOS[0][iorbital][ispin][ienergy];
        vDOS_species[0][iorbital][ispin][ienergy]+=dos;
      }
    }
  }
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " getting vDOS_species[0] (VASP totals)" << endl;}
  for(uint ispin=0;ispin<ISPIN;ispin++){
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " ispin=" << ispin << endl;}
    for(uint ienergy=0;ienergy<IENERGY;ienergy++){
      dos=vDOS[0][0][ispin][ienergy];
      vDOS_species[0][0][ispin][ienergy]+=dos;
    }
  }

  //  for(uint ispecies=0;ispecies<num_each_type.size();ispecies++){
  //    for(uint i=0;i<(uint)num_each_type[ispecies];i++){
  //      for(uint iorbital=1;iorbital<vDOS[iatom+1].size();iorbital++){  //iatom/iorbital total is 0
  //        for(uint ispin=0;ispin<vDOS[iatom+1][iorbital].size();ispin++){
  //          if(LDEBUG){cerr << __AFLOW_FUNC__ << " ispecies=" << ispecies << ", iorbital=" << iorbital << ", ispin=" << ispin << endl;}
  //          for(uint ienergy=0;ienergy<vDOS[iatom+1][iorbital][ispin].size();ienergy++){
  //            dos=vDOS[iatom+1][iorbital][ispin][ienergy];
  //            vDOS_species[ispecies+1][ispin][ienergy]+=dos;
  //          }
  //        }
  //      }
  //      iatom++;
  //    }
  //  }

  return vDOS_species;

}

//ME20190623 BEGIN
ostream& operator<<(ostream& oss, const xDOSCAR& xdos) {
  // Header
  oss << std::setw(4) << xdos.number_atoms
    << std::setw(4) << xdos.number_atoms
    << std::setw(4) << xdos.partial
    << std::setw(4) << 0 << std::endl;
  oss << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  oss << std::setprecision(7) << std::scientific;
  oss << std::setw(15) << xdos.Vol;
  for (int i = 1; i <= 3; i++) oss << std::setw(15) << xdos.lattice[i];
  oss << std::setw(15) << xdos.POTIM << std::endl;
  oss << "  " << std::setprecision(15) << xdos.temperature << std::endl;
  oss << "  " << xdos.carstring << std::endl;
  oss << " " << xdos.title << std::endl;

  stringstream dosline;  // Will be reused for projected DOS
  dosline << std::dec << std::fixed << std::setprecision(8) << std::setw(15) << xdos.energy_max
    << std::fixed << std::setw(15) << xdos.energy_min
    << std::setprecision(0) << "  " << xdos.number_energies
    << std::setprecision(8) << std::fixed << std::setw(15) << xdos.Efermi
    << std::setprecision(8) << std::fixed << std::setw(15) << 1.0;
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

  for (uint p = 1; p < xdos.vDOS.size(); p++) {
    oss << dosline.str() << std::endl;
    for (uint e = 0; e < xdos.number_energies; e++) {
      oss << std::setw(12) << xdos.venergy[e];
      for (uint o = 1; o < xdos.vDOS[p].size(); o++) {  // Do not output the total
        for (uint s = 0; s < xdos.vDOS[p][o].size(); s++) {
          oss << std::setw(12) << xdos.vDOS[p][o][s][e];
        }
      }
      oss << std::endl;
    }
  }
  return oss;
}
//ME20190623 END

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************

// ***************************************************************************
// class xEIGENVAL
//ME20200427 - included xStream::initialize
xEIGENVAL::xEIGENVAL(ostream& oss) : xStream(oss) {free();} //CO20200404 - xStream integration for logging
xEIGENVAL::xEIGENVAL(ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss) {free();} //CO20200404 - xStream integration for logging
xEIGENVAL::xEIGENVAL(const string& fileIN,bool QUIET,ostream& oss) : xStream(oss) {m_initialized=initialize(fileIN,QUIET);} //CO20200404 - xStream integration for logging
xEIGENVAL::xEIGENVAL(const string& fileIN,ofstream& FileMESSAGE,bool QUIET,ostream& oss) : xStream(FileMESSAGE,oss) {m_initialized=initialize(fileIN,QUIET);} //CO20200404 - xStream integration for logging

bool xEIGENVAL::initialize(const string& fileIN,ostream& oss,bool QUIET) {xStream::initialize(oss);return initialize(fileIN,QUIET);} //CO20200508
bool xEIGENVAL::initialize(const string& fileIN,ofstream& FileMESSAGE,ostream& oss,bool QUIET) {xStream::initialize(FileMESSAGE,oss);return initialize(fileIN,QUIET);} //CO20200508
bool xEIGENVAL::initialize(const string& fileIN, bool QUIET) {
  free();
  filename = fileIN;
  return GetPropertiesFile(fileIN,QUIET);
}

xEIGENVAL::xEIGENVAL(const xEIGENVAL& b) : xStream(*b.getOFStream(),*b.getOSS()) {free();copy(b);} // copy PUBLIC

xEIGENVAL::~xEIGENVAL() {xStream::free();free();} //CO20191110 //CO20200404 - xStream integration for logging

void xEIGENVAL::free() {
  m_initialized=false; //CO20200404
  content="";
  vcontent.clear(); 
  filename="";
  title="";
  spin=0;;
  Vol=0.0;POTIM=0.0;
  lattice.clear();
  temperature=0.0;
  number_electrons=0;
  number_kpoints=0;
  number_bands=0;
  energy_max = -1e30;  //ME20190614
  energy_min = 1e30;  //ME20190614
  vweight.clear();
  vkpoint.clear();
  venergy.clear();
  carstring = "";  //ME20190620
  number_atoms = 0;  //ME20190623
  number_loops = 0;  //ME20190623

  //[CO20200404 - DUPLICATE]vcontent.clear(); 
  //[CO20200404 - DUPLICATE]lattice.clear();
  //[CO20200404 - DUPLICATE]vweight.clear();
  //[CO20200404 - DUPLICATE]for(uint i=0;i<vkpoint.size();i++) vkpoint.at(i).clear(); 
  //[CO20200404 - DUPLICATE]vkpoint.clear();
  //[CO20200404 - DUPLICATE]for(uint i=0;i<venergy.size();i++) 
  //[CO20200404 - DUPLICATE]  for(uint j=0;j<venergy.at(i).size();j++) 
  //[CO20200404 - DUPLICATE]    venergy.at(i).at(j).clear(); 
  //[CO20200404 - DUPLICATE]for(uint i=0;i<venergy.size();i++) 
  //[CO20200404 - DUPLICATE]  venergy.at(i).clear(); 
  //[CO20200404 - DUPLICATE]venergy.clear();
}

void xEIGENVAL::copy(const xEIGENVAL& b) { // copy PRIVATE
  xStream::copy(b);  //CO20200404 - xStream integration for logging
  free();
  m_initialized=b.m_initialized; //CO20200404 - xStream integration for logging
  content=b.content;
  vcontent.clear(); 
  for(uint i=0;i<b.vcontent.size();i++) vcontent.push_back(b.vcontent.at(i));  // for aflowlib_libraries.cpp
  filename=b.filename;
  title=b.title;
  spin=b.spin;
  Vol=b.Vol;
  lattice=b.lattice;
  POTIM=b.POTIM;
  temperature=b.temperature;
  number_electrons=b.number_electrons;
  number_kpoints=b.number_kpoints;
  number_bands=b.number_bands;
  energy_min = b.energy_min;  //ME20190614
  energy_max = b.energy_max;  //ME20190614
  vweight.clear(); 
  for(uint i=0;i<b.vweight.size();i++) vweight.push_back(b.vweight.at(i));
  vkpoint.clear(); 
  for(uint i=0;i<b.vkpoint.size();i++) vkpoint.push_back(b.vkpoint.at(i));
  venergy.clear(); 
  for(uint i=0;i<b.venergy.size();i++) {
    deque<deque<double> > temp;
    for(uint j=0;j<b.venergy.at(i).size();j++) 
      temp.push_back(b.venergy.at(i).at(j));
    venergy.push_back(temp);
  }  
  carstring = b.carstring;  //ME20190620
  number_atoms = b.number_atoms;  //ME20190623
  number_loops = b.number_loops;  //ME20190623
}

const xEIGENVAL& xEIGENVAL::operator=(const xEIGENVAL& b) {  // operator= PUBLIC
  if(this!=&b) {free();copy(b);}
  return *this;
}

void xEIGENVAL::clear() {  // clear PRIVATE
  xEIGENVAL _temp;
  string filename_aus=filename;
  copy(_temp);
  filename=filename_aus;
}

bool xEIGENVAL::GetProperties(const string& stringIN,bool QUIET) {
  stringstream sss; sss.str(stringIN);
  if(filename=="") filename="string";
  return xEIGENVAL::GetProperties(sss,QUIET);
}

bool xEIGENVAL::GetPropertiesFile(const string& fileIN,bool QUIET) {
  stringstream sss;
  if(filename=="") filename=fileIN;
  aurostd::efile2stringstream(fileIN,sss);
  return xEIGENVAL::GetProperties(sss,QUIET);
}

bool xEIGENVAL::GetPropertiesUrlFile(const string& url,const string& file,bool QUIET) {
  //[CO20200502 - OBSOLETE]string tmpfile=XHOST.tmpfs+"/_aflow_"+XHOST.user+".pid"+XHOST.ostrPID.str()+".a"+string(AFLOW_VERSION)+".rnd"+aurostd::utype2string(uint((double) std::floor((double)100000*aurostd::ran0())))+".u"+aurostd::utype2string(uint((double) aurostd::get_useconds()))+"_"+file;
  string tmpfile=aurostd::TmpFileCreate("xEIGENVAL_GetProperties"); //CO20200502 - threadID
  aurostd::url2file(url+"/"+file,tmpfile,!QUIET);
  bool out=GetPropertiesFile(tmpfile,QUIET);
  filename="url="+url;  //CO20210315
  aurostd::RemoveFile(tmpfile);
  return out;
}

bool xEIGENVAL::GetProperties(const stringstream& stringstreamIN,bool QUIET) {
  bool LDEBUG=(FALSE || XHOST.DEBUG || !QUIET);
  stringstream message;
  bool force_exit=XHOST.POSTPROCESS; //SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only

  bool ERROR_flag=FALSE;
  long double seconds=aurostd::get_seconds();
  if(LDEBUG) cerr << "xEIGENVAL::GetProperties: BEGIN (" << time_delay(seconds) << ")" << endl;
  clear(); // so it does not mess up vector/deque
  content=stringstreamIN.str();
  vcontent.clear();
  vector<string> vline,tokens;
  //AS20200420 BEGIN]
  // vasp 5.x.x uses '\n' as a separator between blocks of eigenvalues of
  // different k-points, while vasp 4.x.x uses " \n" (space+newline).
  // When a regular string2vectorsting is used, a sequence of "\n\n" is
  // treated as a single '\n'. As a result separating line disappears
  // in vcontent and the line indices are inconsistent.
  // As a result, vweight, vkpoint and venergy vectors are not correctly populated.
  // Setting consecutive = true in string2vectorstring fixes the issue.
  //
  // aurostd::string2vectorstring(content,vcontent);
  aurostd::string2vectorstring(content,vcontent,true,false);
  //AS20200420 END]
  string line;
  if(filename=="") filename="stringstream";
  // crunching to eat the info

  // get parameters
  //  vline.clear();
  for(uint iline=0;iline<vcontent.size();iline++) {
    aurostd::string2tokens(vcontent.at(iline),tokens);
    //    cerr << "iline=" << iline << "  " << vcontent.at(iline) << " tokens.size()=" << tokens.size() << endl;
    if(iline==0 && tokens.size()==4) {
      number_atoms = aurostd::string2utype<int>(tokens[0]);  //ME20190623
      number_loops = aurostd::string2utype<int>(tokens[2]);  //ME20190623
      spin=aurostd::string2utype<int>(tokens.at(tokens.size()-1))-1;
    }
    if(iline==1 && tokens.size()>=5) {
      uint i=0;
      Vol=aurostd::string2utype<double>(tokens.at(i++));
      lattice(1)=aurostd::string2utype<double>(tokens.at(i++));
      lattice(2)=aurostd::string2utype<double>(tokens.at(i++));
      lattice(3)=aurostd::string2utype<double>(tokens.at(i++));
      POTIM=aurostd::string2utype<double>(tokens.at(i++));
    }
    if(iline==2) temperature=aurostd::string2utype<double>(vcontent.at(iline));
    if(iline==3) carstring = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vcontent.at(iline)); //ME20190620 - what kind of EIGENVAL file?
    if(iline==4) title=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vcontent.at(iline));  //ME20190614 - clean title
    if(iline==5 && tokens.size()>=3) {
      uint i=0;
      number_electrons=aurostd::string2utype<uint>(tokens.at(i++));
      number_kpoints=aurostd::string2utype<uint>(tokens.at(i++));
      number_bands=aurostd::string2utype<uint>(tokens.at(i++));
    }
    if(iline>=7) {
      for(uint jline=0;jline<number_kpoints&&iline<vcontent.size();jline++,iline++) { // the iline++ is to get rid of the vacuum
        xvector<double> kpoint(3);
        aurostd::string2tokens(vcontent.at(iline),tokens);
        if(tokens.size()>=4) {
          uint i=0;
          kpoint(1)=aurostd::string2utype<double>(tokens.at(i++));
          kpoint(2)=aurostd::string2utype<double>(tokens.at(i++));
          kpoint(3)=aurostd::string2utype<double>(tokens.at(i++));
          vkpoint.push_back(kpoint);
          vweight.push_back(aurostd::string2utype<double>(tokens.at(i++)));
        }
        iline++; // move to energies
        deque<double> keigenval;
        deque<deque<double> > knenergy;
        double eigen;  //ME20190614
        for(uint kline=0;kline<number_bands&&iline<vcontent.size();kline++,iline++) {
          aurostd::string2tokens(vcontent.at(iline),tokens);
          keigenval.clear();
          for(uint i=1;i<tokens.size();i++) {
            //ME20190614 START
            eigen = aurostd::string2utype<double>(tokens[i]);
            if (eigen < energy_min) energy_min = eigen;
            if (eigen > energy_max) energy_max = eigen;
            keigenval.push_back(eigen);
            //[OBSOLETE] keigenval.push_back(aurostd::string2utype<double>(tokens.at(i)));  // build one eigenvalue
            //ME20190614 END
          }
          knenergy.push_back(keigenval);
        }
        venergy.push_back(knenergy);
      }
    }
  }
  // ----------------------------------------------------------------------
  if(LDEBUG) cerr << "xEIGENVAL::GetProperties: title=" << title << endl;
  if(LDEBUG) cerr << "xEIGENVAL::GetProperties: spin=" << spin << endl;
  if(LDEBUG) cerr << "xEIGENVAL::GetProperties: Vol=" << Vol << endl;
  if(LDEBUG) cerr << "xEIGENVAL::GetProperties: lattice=" << lattice << endl;
  if(LDEBUG) cerr << "xEIGENVAL::GetProperties: POTIM=" << POTIM << endl;
  if(LDEBUG) cerr << "xEIGENVAL::GetProperties: temperature=" << temperature << endl;
  if(LDEBUG) cerr << "xEIGENVAL::GetProperties: number_electrons=" << number_electrons << endl;
  if(LDEBUG) cerr << "xEIGENVAL::GetProperties: number_kpoints=" << number_kpoints << endl;
  if(LDEBUG) cerr << "xEIGENVAL::GetProperties: number_bands=" << number_bands << endl;
  if(LDEBUG) cerr << "xEIGENVAL::GetProperties: vweight.size()=" << vweight.size() << endl;
  if(LDEBUG) cerr << "xEIGENVAL::GetProperties: vkpoint.size()=" << vkpoint.size() << endl;
  if(LDEBUG) cerr << "xEIGENVAL::GetProperties: venergy.size()=" << venergy.size() << endl;
  if(LDEBUG) cerr << "xEIGENVAL::GetProperties: venergy.at(max).size()=" << venergy.at(venergy.size()-1).size() << endl;
  if(LDEBUG) cerr << "xEIGENVAL::GetProperties: venergy.at(max).at(max).size()=" 
    << venergy.at(venergy.size()-1).at(venergy.at(venergy.size()-1).size()-1).size() << endl;

  // for(uint i=0;i<venergy.size();i++)
  //   if(LDEBUG) cerr << "xEIGENVAL::GetProperties: venergy.at.(" << i << ").size()=" << venergy.at(i).size() << endl;
  // for(uint i=0;i<venergy.size();i++)
  //   for(uint j=0;j<venergy.at(i).size();j++)
  //     if(LDEBUG) cerr << "xEIGENVAL::GetProperties: venergy.at.(" << i << ").at.(" << j << ").size()=" << venergy.at(i).at(j).size() << endl;

  // ----------------------------------------------------------------------
  // DONE NOW RETURN  
  if(LDEBUG) cerr << "xEIGENVAL::GetProperties: END (" << time_delay(seconds) << ")" << endl;
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  //[CO20200404 - OBSOLETE]if(ERROR_flag && !QUIET) cerr << "WARNING - xEIGENVAL::GetProperties: ERROR_flag set in xEIGENVAL" << endl;
  if(ERROR_flag){
    message << "ERROR_flag set in xEIGENVAL";
    if(force_exit){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);}
    else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_,QUIET);}
    return FALSE;
  }
  return TRUE;
}

ostream& operator<<(ostream& oss, const xEIGENVAL& xeigen) {
  // Header
  oss << std::setw(4) << xeigen.number_atoms
    << std::setw(4) << xeigen.number_atoms
    << std::setw(4) << xeigen.number_loops
    << std::setw(4) << (xeigen.spin + 1) << std::endl;
  oss << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  oss << std::setprecision(7) << std::scientific;
  oss << std::setw(15) << xeigen.Vol;
  for (int i = 1; i < 4; i++) oss << std::setw(15) << xeigen.lattice[i];
  oss << std::setw(15) << xeigen.POTIM << std::endl;
  oss << "  " << std::setprecision(15) << xeigen.temperature << std::endl;
  oss << "  " << xeigen.carstring << std::endl;
  oss << " " << xeigen.title << std::endl;
  oss << std::dec << std::setw(4) << xeigen.number_electrons
    << "  " << std::setw(4) << xeigen.number_kpoints
    << std::setw(4) << xeigen.number_bands << std::endl;

  // Data
  for (uint k = 0; k < xeigen.number_kpoints; k++) {
    oss << " " << std::endl;  // space MUST be there or the EIGENVAL reader will fail
    oss << std::scientific << std::setprecision(8);
    for (int i = 1; i < 4; i++) oss << "  " << xeigen.vkpoint[k][i];
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
bool GetEffectiveMass(xOUTCAR& xoutcar, xDOSCAR& xdoscar, xEIGENVAL& xeigenval, xstructure xstr,ostream& oss) {ofstream FileMESSAGE;return GetEffectiveMass(xoutcar,xdoscar,xeigenval,xstr,FileMESSAGE,oss);} //CO20200404
bool GetEffectiveMass(xOUTCAR& xoutcar, xDOSCAR& xdoscar, xEIGENVAL& xeigenval, xstructure xstr,ofstream& FileMESSAGE,ostream& oss) { //CO20200404
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

  xmatrix<double> reciprocal_lattice(1,1,3,3);
  vector<vector<vector<kEn_st> > > fit_data_all;
  vector<vector<int> > number_of_valley_list;
  vector<double> band_info_vbt(2, xdoscar.Efermi);
  vector<double> band_info_cbb(2, xdoscar.Efermi);
  vector<int> valley_elec(2), valley_hole(2);
  double _ENER_RANGE = 0.026; //[CO20191004]_METALGAP =  -1.0E09;
  bool SPIN_UP = FALSE, SPIN_DN = FALSE;
  string compound_name = xdoscar.title;
  int ispin = xeigenval.spin+1;

  // Eff. Mass arrays
  xoutcar.band_index.clear();
  xoutcar.carrier_type.clear();
  xoutcar.carrier_spin.clear();
  xoutcar.mass_elec_dos.clear()       ; xoutcar.mass_elec_dos.resize(2);
  xoutcar.mass_hole_dos.clear()       ; xoutcar.mass_hole_dos.resize(2);
  xoutcar.mass_elec_conduction.clear(); xoutcar.mass_elec_conduction.resize(2);
  xoutcar.mass_hole_conduction.clear(); xoutcar.mass_hole_conduction.resize(2);
  if(ispin == 1) {
    if(xoutcar.Egap.at(0) > _METALGAP_) {
      SPIN_UP = TRUE;
      band_info_vbt.at(0) = xoutcar.valence_band_max.at(0)    + xdoscar.Efermi;
      band_info_cbb.at(0) = xoutcar.conduction_band_min.at(0) + xdoscar.Efermi;
      band_info_vbt.at(1) = xoutcar.valence_band_max.at(0)    + xdoscar.Efermi;
      band_info_cbb.at(1) = xoutcar.conduction_band_min.at(0) + xdoscar.Efermi;
    }
  }
  else if(ispin == 2) {
    if(xoutcar.Egap.at(0) > _METALGAP_) {
      SPIN_UP = TRUE;
      band_info_vbt.at(0) = xoutcar.valence_band_max.at(0)    + xdoscar.Efermi;
      band_info_cbb.at(0) = xoutcar.conduction_band_min.at(0) + xdoscar.Efermi;
    }
    if(xoutcar.Egap.at(1) > _METALGAP_) {
      SPIN_DN = TRUE;
      band_info_vbt.at(1) = xoutcar.valence_band_max.at(1)    + xdoscar.Efermi;
      band_info_cbb.at(1) = xoutcar.conduction_band_min.at(1) + xdoscar.Efermi;
    }
  }
  else {
    //[CO20200404 - OBSOLETE]xoutcar.ERROR = "EIGENVAL spin value is neither 0 nor 1."; // probably unnecessary
    message << "EIGENVAL spin value is neither 0 nor 1"; //CO20200404
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_ERROR_); //CO20200404
    return FALSE;
  }
  if(ispin == 1) {
    if(!SPIN_UP) {
      //[CO20200404 - OBSOLETE]xoutcar.ERROR = "Metallic system encountered";
      message << "Metallic system encountered"; //CO20200404
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_ERROR_); //CO20200404
      return FALSE;
    }
  }
  else if(ispin == 2) {
    if(!SPIN_UP and !SPIN_DN) {
      //[CO20200404 - OBSOLETE]xoutcar.ERROR = "Metallic system encountered";
      message << "Metallic system encountered"; //CO20200404
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_ERROR_); //CO20200404
      return FALSE;
    }
  }

  if(SPIN_UP or SPIN_DN) { // this disappears
    xstr.FixLattices();
    xstr.CalculateSymmetryPointGroupKLattice();
    reciprocal_lattice = xstr.klattice;
    vector<vector<kEn_st> > allkE_points;
    // allkE_points.at(eigenvalue).at(kpoint).{kEn_st details}
    allkE_points.resize(xeigenval.number_bands);
    xvector<double> temp_recip(1,3), temp_cart(1,3);
    // loop over KPOINTS (vkpoints)
    for(vector<int>::size_type ix=0; ix != xeigenval.vkpoint.size(); ++ix) {
      kEn_st kp;
      temp_recip[1] = xeigenval.vkpoint.at(ix)[1];
      temp_recip[2] = xeigenval.vkpoint.at(ix)[2];
      temp_recip[3] = xeigenval.vkpoint.at(ix)[3];
      temp_cart     = temp_recip * reciprocal_lattice;
      kp.kpoint(1)  = temp_cart[1];
      kp.kpoint(2)  = temp_cart[2];
      kp.kpoint(3)  = temp_cart[3];
      // loop over EIGENVALUES (venergy)
      for(vector<int>::size_type iy=0; iy != xeigenval.venergy.at(ix).size(); ++iy) {
        kp.band_index = iy+1;
        // SPIN UNPOLARIZED
        if(ispin == 1) {
          kp.energy[0]  = xeigenval.venergy.at(ix).at(iy).at(0);
          kp.energy[1]  = kp.energy[0];
        }
        // SPIN POLARIZED
        else if(ispin == 2) {
          if( SPIN_UP and  SPIN_DN) {
            kp.energy[0]  = xeigenval.venergy.at(ix).at(iy).at(0);
            kp.energy[1]  = xeigenval.venergy.at(ix).at(iy).at(1);
          }
          else if( SPIN_UP and !SPIN_DN) {
            kp.energy[0]  = xeigenval.venergy.at(ix).at(iy).at(0);
            kp.energy[1]  = kp.energy[0];
          }
          else if(!SPIN_UP and  SPIN_DN) {
            kp.energy[0]  = xeigenval.venergy.at(ix).at(iy).at(1);
            kp.energy[1]  = kp.energy[0];
          }
          kp.energy[0]  = xeigenval.venergy.at(ix).at(iy).at(0);
          kp.energy[1]  = xeigenval.venergy.at(ix).at(iy).at(1);
        }
        allkE_points.at(kp.band_index-1).push_back(kp);
      }
    }
    for(int spin_idx=0; spin_idx<ispin; spin_idx++) {
      // CHOOSE BANDS CONTAINING EXTREMES
      vector<int> number_of_valley_list_tmp;
      vector< vector<kEn_st> > fit_data;
      // 1. sort the energy in ascending order
      for(uint ii=0; ii<allkE_points.size(); ii++) {
        if      (spin_idx == 0) sort(allkE_points.at(ii).begin(), allkE_points.at(ii).end(), comparison_kEn_str_up);
        else if(spin_idx == 1) sort(allkE_points.at(ii).begin(), allkE_points.at(ii).end(), comparison_kEn_str_dn);
      }
      // 2. get the points with energy in the range
      for(uint ii=0; ii<allkE_points.size(); ii++) {
        // These 'if-then' statements determine the number of pockets in the system
        vector<kEn_st> fit_data_band;
        // CONDUCTION BANDS
        if(abs(allkE_points.at(ii).front().energy[spin_idx] - band_info_cbb.at(spin_idx)) < _ENER_RANGE) {
          double band_energy_minimum = allkE_points.at(ii).front().energy[spin_idx];
          for(int jj=0; jj<_FIT_POINTS_NUMBER; jj++) {
            fit_data_band.push_back(allkE_points.at(ii).at(jj));
            fit_data_band.back().band_type = 1;
          }
          for(uint jj=_FIT_POINTS_NUMBER; jj<allkE_points.at(ii).size(); jj++) {
            if(allkE_points.at(ii).at(jj).energy[spin_idx] - band_energy_minimum < _FIT_ENERGY_RANGE) {
              fit_data_band.push_back(allkE_points.at(ii).at(jj));
              fit_data_band.back().band_type = 1;
            }
            else {
              break;
            }
          }
        }
        // VALENCE BANDS
        else if(abs(band_info_vbt.at(spin_idx) - allkE_points.at(ii).back().energy[spin_idx]) < _ENER_RANGE) {
          double band_energy_maximum = allkE_points.at(ii).back().energy[spin_idx];
          for(uint jj=allkE_points.at(ii).size()-1; jj>allkE_points.at(ii).size() - _FIT_POINTS_NUMBER-1; jj--) {
            fit_data_band.push_back(allkE_points.at(ii).at(jj));
            fit_data_band.back().band_type = 0;
          }
          for(int jj=allkE_points.at(ii).size() - _FIT_POINTS_NUMBER-1; jj>=0; jj--) {
            if(band_energy_maximum-allkE_points.at(ii).at(jj).energy[spin_idx] < _FIT_ENERGY_RANGE) {
              fit_data_band.push_back(allkE_points.at(ii).at(jj));
              fit_data_band.back().band_type = 0;
            }
            else {
              break;
            }
          }
        }
        // fit_data contains the spin-specific pocket information
        if(fit_data_band.size() > 0) {
          fit_data.push_back(fit_data_band);
        }
      }
      vector<double> max_distance;
      for(int i=1; i<4; i++) {
        double max_tmp;
        max_tmp  = std::max(abs(reciprocal_lattice[1][i]), abs(reciprocal_lattice[2][i]));
        max_tmp  = std::max(abs(reciprocal_lattice[3][i]), max_tmp);
        max_tmp *= _BANDS_PARAMETER_MIN_RATIO;
        max_distance.push_back(max_tmp);
      }
      vector<vector<kEn_st> > fit_data_new;
      for(uint i=0; i< fit_data.size(); i++) {
        vector<kEn_st> fit_data_band;
        kEn_st kp;
        xvector<double> pt(1,3);
        // the first point is closest to the extremes
        kp    = fit_data.at(i).at(0);
        pt[1] = kp.kpoint(1);
        pt[2] = kp.kpoint(2);
        pt[3] = kp.kpoint(3);
        for(uint ii=0; ii<fit_data.at(i).size(); ii++) {
          xvector<double> pt1(1,3), pt_sym(1,3);
          kEn_st kp1;
          // the first point is closest to the extremes
          kp1    = fit_data.at(i).at(ii);
          pt1[1] = kp1.kpoint(1);
          pt1[2] = kp1.kpoint(2);
          pt1[3] = kp1.kpoint(3);
          for(uint j=0; j<xstr.pgroupk.size(); j++) {
            pt_sym = pt1 * xstr.pgroupk.at(j).Uc;
            if(near_to(pt, pt_sym, max_distance)) {
              // compare the distance between the most extreme points and the generated one
              kEn_st k1;
              k1.kpoint(1)  = pt_sym[1];
              k1.kpoint(2)  = pt_sym[2];
              k1.kpoint(3)  = pt_sym[3];
              k1.energy[0]  = kp1.energy[0];
              k1.energy[1]  = kp1.energy[1];
              k1.band_index = kp1.band_index;
              k1.band_type  = kp1.band_type;
              fit_data_band.push_back(k1);
            }
          }
          if(ii == 0) {
            int number_of_valley=xstr.pgroupk.size()/fit_data_band.size();
            number_of_valley_list_tmp.push_back(number_of_valley);
          }
        }
        vector<kEn_st>::iterator it;
        sort(fit_data_band.begin(), fit_data_band.end(), comparison_kEn_str_position);
        it = unique(fit_data_band.begin(), fit_data_band.end(), is_equal_position_kEn_str);
        fit_data_band.resize(it-fit_data_band.begin());
        if(spin_idx == 1) {
          sort(fit_data_band.begin(), fit_data_band.end(), comparison_kEn_str_band_type_up);
        }
        else {
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
    vector<int> index_of_extremas_zero;  
    int number_of_records = 0;
    for(int spin_idx=0; spin_idx<ispin; spin_idx++) {
      // fit_data.size(): number of effective masses detected
      vector<vector<kEn_st> > fit_data = fit_data_all.at(spin_idx);
      vector<vector<double> > mass_eff_list;
      for(uint ii=0; ii<fit_data.size(); ii++) {
        kEn_st kp1    = fit_data.at(ii).at(0);
        //    cout << "kp1 " << kp1.kpoint << endl;
        int    nrow   = fit_data.at(ii).size()-1;
        int    ncol   = 9; // 9 polynomial coefficients, (a) thru (i)
        xvector<double> y_vec(1, nrow); // = En2 - En1
        xvector<double> y_sig(1, nrow);
        xmatrix<double> x_mat(1, 1, nrow, ncol);
        for(int jj=1; jj<nrow+1; jj++) {
          y_vec[jj] = fit_data.at(ii).at(jj).energy[spin_idx] - kp1.energy[spin_idx]; // = En2 - En1
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
        for(int jj=1; jj<nrow+1; jj++) {
          kEn_st kp2   = fit_data.at(ii).at(jj);
          //     cout << "kp2 " << kp2.kpoint << endl;
          x_mat[jj][1] = kp2.kpoint(1)*kp2.kpoint(1) - kp1.kpoint(1)*kp1.kpoint(1); // for(a)
          x_mat[jj][2] = kp2.kpoint(2)*kp2.kpoint(2) - kp1.kpoint(2)*kp1.kpoint(2); // for(b)
          x_mat[jj][3] = kp2.kpoint(3)*kp2.kpoint(3) - kp1.kpoint(3)*kp1.kpoint(3); // for(c)
          x_mat[jj][4] = kp2.kpoint(1)*kp2.kpoint(2) - kp1.kpoint(1)*kp1.kpoint(2); // for(d)
          x_mat[jj][5] = kp2.kpoint(1)*kp2.kpoint(3) - kp1.kpoint(1)*kp1.kpoint(3); // for(e)
          x_mat[jj][6] = kp2.kpoint(2)*kp2.kpoint(3) - kp1.kpoint(2)*kp1.kpoint(3); // for(f)
          x_mat[jj][7] = kp2.kpoint(1)               - kp1.kpoint(1);               // for(g)
          x_mat[jj][8] = kp2.kpoint(2)               - kp1.kpoint(2);               // for(h)
          x_mat[jj][9] = kp2.kpoint(3)               - kp1.kpoint(3);               // for(i)
        }
        // [x_mat] [a,b,c,d,e,f,g,h,i]=[y_vec]
        // check x_mat for columns of full of 0's - SINGULAR MATRIX
        // x_mat[rows][columns] <<-- same as Fortran
        //  cout << std::showpos;
        for(int rows=1; rows <= nrow; rows++) {
          bool SINGULAR = TRUE;
          for(int cols=1; cols <= ncol; cols++) {
            //      cout << setw(10) << std::left << std::scientific << x_mat[rows][cols] << "  ";
            if(x_mat[rows][cols] != 0.0) {
              SINGULAR = FALSE;
            }
          }
          //  cout << endl;
          if(SINGULAR) {
            //[CO20200404 - OBSOLETE]xoutcar.ERROR = "Singular system: ill-defined matrix problem encountered.";
            message << "Singular system: ill-defined matrix problem encountered."; //CO20200404
            pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_WARNING_); //CO20200404
            return FALSE;
          }
        }
        //    cout << "=====================================" << endl;
        aurostd::cematrix ECI_matrix(x_mat);
        ECI_matrix.LeastSquare(y_vec, y_sig); // minimization happens here
        // starting @ 1, ending at 3 and 3x3
        xmatrix<double> mass_m(1,1,3,3);
        // looks like upper triangular stuff happens here
        mass_m[1][1] = ECI_matrix.AVec().at(0);
        mass_m[2][2] = ECI_matrix.AVec().at(1);
        mass_m[3][3] = ECI_matrix.AVec().at(2);
        mass_m[1][2] = ECI_matrix.AVec().at(3)*0.5;
        mass_m[1][3] = ECI_matrix.AVec().at(4)*0.5;
        mass_m[2][3] = ECI_matrix.AVec().at(5)*0.5;
        mass_m[2][1] = mass_m[1][2];
        mass_m[3][1] = mass_m[1][3];
        mass_m[3][2] = mass_m[2][3];
        aurostd::cematrix    mass_m_ce(mass_m);
        xvector<double> mr = mass_m_ce.EigenValues();
        vector<double>  mr_tmp;
        for(int jj=1; jj<=3; jj++) {
          mr[jj] = 1.0*_MASS_FACTOR/mr[jj];
          mr_tmp.push_back(mr[jj]);
        }
        sort(mr_tmp.begin(),mr_tmp.end());
        mass_eff_list.push_back(mr_tmp);
      }

      vector<double> elec_cond_mass_per_valley;
      vector<double> hole_cond_mass_per_valley;
      vector<int>    elec_valley_multiplicity;
      vector<int>    hole_valley_multiplicity;
      xoutcar.mass_elec_dos.at(spin_idx) = 0.0; 
      xoutcar.mass_hole_dos.at(spin_idx) = 0.0;
      valley_elec.at(spin_idx) = 0; 
      valley_hole.at(spin_idx) = 0;

      // problem starts here
      for(uint ii=0; ii<mass_eff_list.size(); ii++) {
        vector<double> tempvector(3);
        if      (fit_data.at(ii).at(0).band_type == 0) xoutcar.carrier_type.push_back("hole");
        else if(fit_data.at(ii).at(0).band_type == 1) xoutcar.carrier_type.push_back("elec");
        xoutcar.band_index.push_back(fit_data.at(ii).at(0).band_index);
        xoutcar.carrier_spin.push_back(spin_idx);
        xoutcar.extrema_cart_coord.push_back(tempvector);
        xoutcar.effective_mass_axes.push_back(tempvector);
        number_of_records++;
        double dos_mass  = 1.0;
        double cond_mass = 0.0;
        // diagonalized eff mass tensor
        for(uint jj=0; jj<3; jj++) {
          dos_mass   *= mass_eff_list.at(ii).at(jj);
          double one  = 1.0;
          // add the reciprocal of 3 individual masses
          cond_mass                                 += one / mass_eff_list.at(ii).at(jj);
          // the following causes out of bounds errors 
          xoutcar.extrema_cart_coord.back().at(jj)   = fit_data_all.at(0).at(ii).at(0).kpoint(jj+1);
          xoutcar.effective_mass_axes.back().at(jj)  = mass_eff_list.at(ii).at(jj);
        }
        int number_of_valley=number_of_valley_list.at(spin_idx).at(ii);
        xoutcar.equivalent_valley.push_back(number_of_valley_list.at(spin_idx).at(ii));
        dos_mass *= number_of_valley*number_of_valley;
        if(xoutcar.carrier_type.back() == "elec") {
          xoutcar.mass_elec_dos.at(spin_idx) += dos_mass;
          valley_elec.at(spin_idx)++;
        }
        else if(xoutcar.carrier_type.back() == "hole") {
          xoutcar.mass_hole_dos.at(spin_idx) += dos_mass;
          valley_hole.at(spin_idx)++;
        }
        dos_mass = std::pow((double) dos_mass,(double) 1.0/3.0);
        double temp = 1.0;
        temp /= cond_mass;
        // flip sum of inverses to denominator
        cond_mass = temp;
        // store cond_mass for the valley and store multiplicity of the valley
        if      (xoutcar.carrier_type.back() == "elec") {
          elec_cond_mass_per_valley.push_back(cond_mass);
          elec_valley_multiplicity.push_back(number_of_valley);
        }
        else if(xoutcar.carrier_type.back() == "hole") {
          hole_cond_mass_per_valley.push_back(cond_mass);
          hole_valley_multiplicity.push_back(number_of_valley);
        }
        // final factor of 3 for single valley total
        cond_mass *= 3.0;
        xoutcar.effective_mass_DOS.push_back(dos_mass);
        xoutcar.effective_mass_COND.push_back(cond_mass);
      } // END: loop over mass_eff_list
      // problem ends here

      // Calculate the total electron conductivity effective mass
      //       double numerator   = 0.0; double denominator = 0.0;
      //       for(vector<int>::size_type ix=0; ix != elec_cond_mass_per_valley.size(); ++ix) {
      //          numerator   += elec_valley_multiplicity.at(ix);
      //          denominator += elec_valley_multiplicity.at(ix)/elec_cond_mass_per_valley.at(ix);
      //       }
      //       numerator *= 3;
      //       xoutcar.mass_elec_conduction.at(spin_idx) = numerator/denominator;
      //       // Calculate the total hole conductivity effective mass
      //       numerator   = 0.0;
      //       denominator = 0.0;
      //       for(vector<int>::size_type ix=0; ix != hole_cond_mass_per_valley.size(); ++ix) {
      //          numerator   += hole_valley_multiplicity.at(ix);
      //          denominator += hole_valley_multiplicity.at(ix)/hole_cond_mass_per_valley.at(ix);
      //       }
      //       numerator *= 3;
      //       xoutcar.mass_hole_conduction.at(spin_idx) = numerator / denominator;
    }
    // end Spin Loop

    // DOS electron effective mass for all spins
    //    for(int spin_idx=0; spin_idx<ispin; spin_idx++) {
    //       xoutcar.mass_elec_dos.at(spin_idx) /= valley_elec.at(spin_idx);
    //       xoutcar.mass_elec_dos.at(spin_idx)  = std::pow((double)xoutcar.mass_elec_dos.at(spin_idx),(double)1.0/3.0);
    //    }
    //    // DOS hole effective mass for all spins
    //    for(int spin_idx=0; spin_idx<ispin; spin_idx++) {
    //       xoutcar.mass_hole_dos.at(spin_idx) /= valley_hole.at(spin_idx);
    //       xoutcar.mass_hole_dos.at(spin_idx)  = std::pow((double)xoutcar.mass_hole_dos.at(spin_idx),(double)  1.0/3.0);
    //    }
  }
  // if(ispin==1) {
  //    xoutcar.mass_elec_dos.pop_back();
  //    xoutcar.mass_hole_dos.pop_back();
  //    xoutcar.mass_elec_conduction.pop_back();
  //    xoutcar.mass_hole_conduction.pop_back();
  // } // END: EffMass 

  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  return TRUE;
}
/////////////////////////////////////////////////////////////////////////////////////
// spin dependent comparisons
bool comparison_kEn_str_up(const kEn_st& k1, const kEn_st& k2) {
  return static_cast<bool>(k1.energy[0] < k2.energy[0]);
}
bool comparison_kEn_str_dn(const kEn_st& k1, const kEn_st& k2) {
  return static_cast<bool>(k1.energy[1] < k2.energy[1]);
}
bool comparison_kEn_str_band_type_up(const kEn_st & k1, const kEn_st & k2) {
  if(k1.band_type == 1 && k2.band_type == 1) {  // both are in conduction bands
    return static_cast<bool>(k1.energy[0] < k2.energy[0]);
  }
  else {
    if(k1.band_type == 0 && k2.band_type == 0) { // both are in valence bands
      return static_cast<bool>(k1.energy[0] > k2.energy[0]);
    }
    else {
      return static_cast<bool>(k1.energy[0] < k2.energy[0]);
    }
  }
}
bool comparison_kEn_str_band_type_dn(const kEn_st & k1, const kEn_st & k2) {
  if(k1.band_type == 1 && k2.band_type == 1) {
    return static_cast<bool>(k1.energy[1] < k2.energy[1]);
  }
  else {
    if(k1.band_type == 0 && k2.band_type == 0) {
      return static_cast<bool>(k1.energy[1] > k2.energy[1]);
    }
    else {
      return static_cast<bool>(k1.energy[1] < k2.energy[1]);
    }
  }
}
/////////////////////////////////////////////////////////////////////////////////////
bool comparison_kEn_str_position(const kEn_st & k1, const kEn_st & k2) {
  const double _ENER_EPS=1.0e-4;
  bool flag=FALSE;
  if(abs(k1.kpoint(1) - k2.kpoint(1)) < _ENER_EPS) {
    if(abs(k1.kpoint(2) - k2.kpoint(2)) < _ENER_EPS) {
      if(k1.kpoint(3) < k2.kpoint(3)) {
        flag = TRUE;
      }
      else {
        flag = FALSE;
      }
    }
    else {
      if(k1.kpoint(2) < k2.kpoint(2)) {
        flag = TRUE;
      }
      else {
        flag = FALSE;
      }
    }
  }
  else {
    if(k1.kpoint(1) < k2.kpoint(1)) {
      flag = TRUE;
    }
    else {
      flag = FALSE;
    }
  }  
  return flag;
}
bool is_equal_position_kEn_str(const kEn_st & k1, const kEn_st & k2) {
  const double _ENER_EPS=1.0e-4;
  return isequal(k1.kpoint, k2.kpoint,double(_ENER_EPS));
}
bool near_to(const xvector<double> & k1, const xvector<double> & k2, const vector<double> & max_distance) {
  return static_cast<bool>(abs(k1[1]-k2[1])<max_distance.at(0) && 
      abs(k1[2]-k2[2])<max_distance.at(1) && 
      abs(k1[3]-k2[3])<max_distance.at(2));
}
//-------------------------------------------------------------------------------------------------
// PrintBandGap: Print the output of xOUTCAR::GetBandGap
// 2014: Camilo E. Calderon
// 2017: Corey Oses
bool PrintBandGap(string& directory, ostream &oss) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream message;  //CO20200404
  stringstream ss_outcar_static(""),ss_outcar_bands("");
  string path_outcar_static="",path_outcar_bands="",path_POSCAR="";
  xOUTCAR xoutcar_static,xoutcar_bands;
  char LastChar = *directory.rbegin();
  if(LastChar == '/') directory.erase(directory.size()-1);

  long double seconds=aurostd::get_seconds();
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN (" << time_delay(seconds) << ")" << endl;}

  // OUTCAR_bands
  //[CO20200404 - OBSOLETE]if(!ss_outcar_bands.str().length() && aurostd::FileExist(directory+"/OUTCAR.bands",path_outcar_bands)) aurostd::file2stringstream(path_outcar_bands, ss_outcar_bands); //CO20200404
  if(!ss_outcar_bands.str().length() && aurostd::EFileExist(directory+"/OUTCAR.bands",path_outcar_bands)) aurostd::efile2stringstream(path_outcar_bands, ss_outcar_bands); //CO20200404
  if(!ss_outcar_bands.str().length()) {
    message << "OUTCAR.bands not found"; //CO20200404
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, oss, _LOGGER_WARNING_); //CO20200404
    //[CO20200404 - OBSOLETE]oss << "WARNING " << __AFLOW_FUNC__ << " OUTCAR.bands not found here: " << directory << endl;
    //[CO20200404 - OBSOLETE]oss << endl;
    //[CO20200404 - OBSOLETE]return FALSE;
  }  //CO20200404
  //[CO20200404 - OBSOLETE]if(!ss_outcar_bands.str().length() && aurostd::FileExist(directory+"/OUTCAR",path_outcar_bands)) aurostd::file2stringstream(path_outcar_bands, ss_outcar_bands); //CO20200404
  if(!ss_outcar_bands.str().length() && aurostd::EFileExist(directory+"/OUTCAR",path_outcar_bands)) aurostd::efile2stringstream(path_outcar_bands, ss_outcar_bands); //CO20200404
  if(!ss_outcar_bands.str().length()) {  //CO20200404
    message << "OUTCAR not found"; //CO20200404
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, oss, _LOGGER_ERROR_); //CO20200404
    return FALSE;
  }

  //CO20171002 - using tolerance from symmetry calc - START
  //double tol;
  //if(xstr.CalculateSymmetry()){tol=xstr.sym_eps;}
  //else {tol=SYM::defaultTolerance(xstr);}
  //cerr << tol << endl;
  //CO20171002 - using tolerance from symmetry calc - STOP

  if(!xoutcar_bands.GetProperties(ss_outcar_bands)){
    //[CO20200404 - OBSOLETE]pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, xoutcar_bands.ERROR, oss, _LOGGER_ERROR_); //CO20200404
    //[CO20200404 - OBSOLETE]oss << "WARNING " << __AFLOW_FUNC__ << " " << xoutcar_bands.ERROR << endl << "   filename=[" << path_outcar_bands << "]" << endl;
    //[CO20200404 - OBSOLETE]oss << endl;
    return FALSE;
  }
  //try to grab xstr from OUTCAR
  if(!xoutcar_bands.GetXStructure()){
    if(!aurostd::FileExist(directory+"/POSCAR.bands",path_POSCAR) && !aurostd::EFileExist(directory+"/POSCAR.bands",path_POSCAR)) {
      //[CO20200404 - OBSOLETE]pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, xoutcar_bands.ERROR, oss, _LOGGER_ERROR_); //CO20200404
      //[CO20200404 - OBSOLETE]oss << "WARNING " << __AFLOW_FUNC__ << " " << xoutcar_bands.ERROR << endl << "   filename=[" << path_outcar_bands << "]" << endl;
      //[CO20200404 - OBSOLETE]oss << endl;
      return FALSE;
    }
    xstructure xstr(path_POSCAR,IOVASP_POSCAR);
    xoutcar_bands.xstr=xstr;
  }

  double EFERMI=xoutcar_bands.Efermi;  //hopefully we can grab this from static, otherwise settle on the one in bands

  // OUTCAR_static - try to grab right Efermi
  //[CO20200404 - OBSOLETE]if(!aurostd::FileExist(directory+"/OUTCAR.static",path_outcar_static) && !aurostd::EFileExist(directory+"/OUTCAR.static",path_outcar_static))
  if(!aurostd::EFileExist(directory+"/OUTCAR.static",path_outcar_static)){
    message << "OUTCAR.static not found, defaulting E-Fermi to that found in OUTCAR.bands";  //CO20200404
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, oss, _LOGGER_WARNING_); //CO20200404
    //[CO20200404 - OBSOLETE]oss << "WARNING " << __AFLOW_FUNC__ << " OUTCAR.static not found here: " << directory << endl;
    //[CO20200404 - OBSOLETE]oss << "WARNING " << __AFLOW_FUNC__ << " Defaulting E-Fermi to that found in OUTCAR.bands" << endl;
    //[CO20200404 - OBSOLETE]oss << endl;
    //return FALSE;
  } else {
    //[CO20200404 - OBSOLETE]if(!xoutcar_static.GetProperties(ss_outcar_static) && aurostd::FileExist(directory+"/OUTCAR.static",path_outcar_static)) aurostd::file2stringstream(path_outcar_static, ss_outcar_static); //CO20200404
    aurostd::efile2stringstream(path_outcar_static, ss_outcar_static); //CO20200404
    if(!xoutcar_static.GetProperties(ss_outcar_static)){
      //[CO20200404 - OBSOLETE]message << xoutcar_bands.ERROR << endl;  //CO20200404
      message << "Defaulting E-Fermi to that found in OUTCAR.bands"; //CO20200404
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, oss, _LOGGER_WARNING_); //CO20200404
      //[CO20200404 - OBSOLETE]oss << "WARNING " << __AFLOW_FUNC__ << " " << xoutcar_static.ERROR << endl << "   filename=[" << path_outcar_static << "]" << endl;
      //[CO20200404 - OBSOLETE]oss << "WARNING " << __AFLOW_FUNC__ << " Defaulting E-Fermi to that found in OUTCAR.bands" << endl;
      //[CO20200404 - OBSOLETE]oss << endl;
      //return FALSE;
    } else {
      EFERMI=xoutcar_static.Efermi;
      if(LDEBUG) {cerr << XPID << "xOUTCAR::PrintBandGap: Found E-fermi from OUTCAR.static: " << EFERMI << endl;}
    }
  }

  if(!xoutcar_bands.GetBandGap(EFERMI)){
    //[CO20200404 - OBSOLETE]pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, xoutcar_bands.ERROR, oss, _LOGGER_ERROR_); //CO20200404
    //[CO20200404 - OBSOLETE]oss << "WARNING " << __AFLOW_FUNC__ << " " << xoutcar_bands.ERROR << endl << "   filename=[" << path_outcar_bands << "]" << endl;
    //[CO20200404 - OBSOLETE]oss << endl;
    return FALSE;
  }

  //could use pointers to save memory, but these objects are so small it doesn't matter
  //make full copies and forget
  string SYSTEM=xoutcar_bands.SYSTEM;
  vector<double> valence_band_max=xoutcar_bands.valence_band_max;
  vector<double> conduction_band_min=xoutcar_bands.conduction_band_min;
  vector<double> Egap=xoutcar_bands.Egap;
  vector<double> Egap_fit=xoutcar_bands.Egap_fit;
  vector<string> Egap_type=xoutcar_bands.Egap_type;
  double valence_band_max_net=xoutcar_bands.valence_band_max_net;
  double conduction_band_min_net=xoutcar_bands.conduction_band_min_net;
  double Egap_net=xoutcar_bands.Egap_net;
  double Egap_fit_net=xoutcar_bands.Egap_fit_net;
  string Egap_type_net=xoutcar_bands.Egap_type_net;

  //SUCCESS!
  if(Egap.size() == 1) { 
    oss.precision(4);
    oss << "System        :   " << SYSTEM << endl; 
    oss << "Spin tag      :   " << Egap.size() << endl;
    oss << "Fermi level   :  "  << std::scientific << std::showpos << EFERMI << endl;
    oss << "                  VBT           CBB           Egap          Egap_fit     Type" << endl;
    oss << "Net Result    :  ";
    oss << setw(14) << std::left << std::scientific << std::showpos << valence_band_max.at(0);
    oss << setw(14) << std::left << std::scientific << std::showpos << conduction_band_min.at(0);
    oss << setw(14) << std::left << std::scientific << std::showpos << Egap.at(0);
    oss << setw(14) << std::left << std::scientific << std::showpos << Egap_fit.at(0);
    oss << setw(20) << std::left << Egap_type.at(0) << endl;
    oss << endl;
  }
  else if(Egap.size() == 2) { 
    oss.precision(4);
    oss << "System        :   " << SYSTEM << endl; 
    oss << "Spin tag      :   " << Egap.size() << endl;
    oss << "Fermi level   :  "  << std::scientific << std::showpos << EFERMI << endl;
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
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " END (" << time_delay(seconds) << ")" << endl;}
  return TRUE;
}
//-------------------------------------------------------------------------------------------------
// PrintBandGap: Print the output of xDOSCAR::GetBandGap
// 2019: Corey Oses
bool PrintBandGap_DOS(string& directory, ostream &oss) { //CO20191110
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream message;  //CO20200404
  char LastChar = *directory.rbegin();
  if(LastChar == '/') directory.erase(directory.size()-1);

  long double seconds=aurostd::get_seconds();
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN (" << time_delay(seconds) << ")" << endl;}

  string path_doscar_static="";
  stringstream ss_doscar_static;
  if(!ss_doscar_static.str().length() && aurostd::FileExist(directory+"/DOSCAR.static",path_doscar_static)) aurostd::file2stringstream(path_doscar_static, ss_doscar_static); //CO20200404
  if(!ss_doscar_static.str().length() && aurostd::EFileExist(directory+"/DOSCAR.static",path_doscar_static)) aurostd::efile2stringstream(path_doscar_static, ss_doscar_static); //CO20200404
  if(!ss_doscar_static.str().length()) {
    message << "DOSCAR.static not found";  //CO20200404
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, oss, _LOGGER_WARNING_); //CO20200404
    //[CO20200404 - OBSOLETE]oss << "ERROR " << __AFLOW_FUNC__ << " DOSCAR.static not found here: " << directory << endl;
    //[CO20200404 - OBSOLETE]oss << endl;
    //[CO20200404 - OBSOLETE]return FALSE;
  }  //CO20200404
  if(!ss_doscar_static.str().length() && aurostd::FileExist(directory+"/DOSCAR",path_doscar_static)) aurostd::file2stringstream(path_doscar_static, ss_doscar_static);  //CO20200404
  if(!ss_doscar_static.str().length() && aurostd::EFileExist(directory+"/DOSCAR",path_doscar_static)) aurostd::efile2stringstream(path_doscar_static, ss_doscar_static);  //CO20200404
  if(!ss_doscar_static.str().length()) { //CO20200404
    message << "DOSCAR not found"; //CO20200404
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, oss, _LOGGER_ERROR_); //CO20200404
    return FALSE;
  }

  xDOSCAR xdoscar;
  if(!xdoscar.GetProperties(ss_doscar_static)){
    //[CO20200404 - OBSOLETE]pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, xdoscar.ERROR, oss, _LOGGER_ERROR_); //CO20200404
    //[CO20200404 - OBSOLETE]oss << "ERROR " << __AFLOW_FUNC__ << " " << xdoscar.ERROR << endl << "   filename=[" << path_doscar_static << "]" << endl;
    //[CO20200404 - OBSOLETE]oss << endl;
    return FALSE;
  }
  if(!xdoscar.GetBandGap()){ //use EFERMI of doscar.static (most accurate)
    //[CO20200404 - OBSOLETE]pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, xdoscar.ERROR, oss, _LOGGER_ERROR_); //CO20200404
    //[CO20200404 - OBSOLETE]oss << "ERROR " << __AFLOW_FUNC__ << " " << xdoscar.ERROR << endl << "   filename=[" << path_doscar_static << "]" << endl;
    //[CO20200404 - OBSOLETE]oss << endl;
    return FALSE;
  }

  //could use pointers to save memory, but these objects are so small it doesn't matter
  //make full copies and forget
  double EFERMI=xdoscar.Efermi;
  string SYSTEM=xdoscar.title;
  vector<double> valence_band_max=xdoscar.valence_band_max;
  vector<double> conduction_band_min=xdoscar.conduction_band_min;
  vector<double> Egap=xdoscar.Egap;
  vector<double> Egap_fit=xdoscar.Egap_fit;
  vector<string> Egap_type=xdoscar.Egap_type;
  double valence_band_max_net=xdoscar.valence_band_max_net;
  double conduction_band_min_net=xdoscar.conduction_band_min_net;
  double Egap_net=xdoscar.Egap_net;
  double Egap_fit_net=xdoscar.Egap_fit_net;
  string Egap_type_net=xdoscar.Egap_type_net;

  //SUCCESS!
  if(Egap.size() == 1) { 
    oss.precision(4);
    oss << "System        :   " << SYSTEM << endl; 
    oss << "Spin tag      :   " << Egap.size() << endl;
    oss << "Fermi level   :  "  << std::scientific << std::showpos << EFERMI << endl;
    oss << "                  VBT           CBB           Egap          Egap_fit     Type" << endl;
    oss << "Net Result    :  ";
    oss << setw(14) << std::left << std::scientific << std::showpos << valence_band_max.at(0);
    oss << setw(14) << std::left << std::scientific << std::showpos << conduction_band_min.at(0);
    oss << setw(14) << std::left << std::scientific << std::showpos << Egap.at(0);
    oss << setw(14) << std::left << std::scientific << std::showpos << Egap_fit.at(0);
    oss << setw(20) << std::left << Egap_type.at(0) << endl;
    oss << endl;
  }
  else if(Egap.size() == 2) { 
    oss.precision(4);
    oss << "System        :   " << SYSTEM << endl; 
    oss << "Spin tag      :   " << Egap.size() << endl;
    oss << "Fermi level   :  "  << std::scientific << std::showpos << EFERMI << endl;
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
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " END (" << time_delay(seconds) << ")" << endl;}
  return TRUE;
}
//---------------------------------------------------------------------------------
// PrintEffectiveMass: Print the output of xOUTCAR::GetEffectiveMass 
// 2014: Camilo E. Calderon
bool PrintEffectiveMass(string& directory, ostream &oss) {
  //[CO20191004]double _METALGAP = -1.0E09;
  bool SPIN_UP = FALSE, SPIN_DN = FALSE;
  bool EM_TAG;
  stringstream file_OUTCAR, file_DOSCAR, file_EIGENVAL;
  string path_POSCAR;

  deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");

  xOUTCAR xoutcar;
  xDOSCAR xdoscar;
  xEIGENVAL xeigenval;
  char LastChar = *directory.rbegin();
  if(LastChar == '/') directory.erase(directory.size()-1);

  // OUTCAR
  bool found_OUTCAR=FALSE;
  for(uint iext=0;iext<vext.size();iext++) { 
    if(!found_OUTCAR && aurostd::FileExist(directory+"/OUTCAR.static"+vext.at(iext))) {
      found_OUTCAR=TRUE;
      aurostd::efile2stringstream(directory+"/OUTCAR.static"+vext.at(iext), file_OUTCAR); // .EXT
    }
  }
  if(!found_OUTCAR && aurostd::FileExist(directory+"/OUTCAR.static")) {
    found_OUTCAR=TRUE;
    aurostd::file2stringstream(directory+"/OUTCAR.static", file_OUTCAR); // plain text
  }
  if(!found_OUTCAR) return FALSE;

  // DOSCAR
  bool found_DOSCAR=FALSE;
  for(uint iext=0;iext<vext.size();iext++) { 
    if(!found_DOSCAR && aurostd::FileExist(directory+"/DOSCAR.static"+vext.at(iext))) {
      found_DOSCAR=TRUE;
      aurostd::efile2stringstream(directory+"/DOSCAR.static"+vext.at(iext), file_DOSCAR); // .EXT
    }
  }
  if(!found_DOSCAR && aurostd::FileExist(directory+"/DOSCAR.static")) {
    found_DOSCAR=TRUE;
    aurostd::file2stringstream(directory+"/DOSCAR.static", file_DOSCAR); // plain text
  }
  if(!found_DOSCAR) return FALSE;

  // EIGENVAL
  bool found_EIGENVAL=FALSE;
  for(uint iext=0;iext<vext.size();iext++) { 
    if(!found_EIGENVAL && aurostd::FileExist(directory+"/EIGENVAL.static"+vext.at(iext))) {
      found_EIGENVAL=TRUE;
      aurostd::efile2stringstream(directory+"/EIGENVAL.static"+vext.at(iext), file_EIGENVAL);  // .EXT
    }
  }
  if(!found_EIGENVAL && aurostd::FileExist(directory+"/EIGENVAL.static")) {
    found_EIGENVAL=TRUE;
    aurostd::file2stringstream(directory+"/EIGENVAL.static", file_EIGENVAL); // plain text
  }
  if(!found_EIGENVAL) return FALSE;

  // POSCAR
  bool found_POSCAR=FALSE;
  for(uint iext=0;iext<vext.size();iext++) { 
    if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR.bands"+vext.at(iext))) {
      found_POSCAR=TRUE;
      path_POSCAR=directory+"/POSCAR.bands"+vext.at(iext); // EXR
    }
  }
  if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR.bands")) {
    found_POSCAR=TRUE;
    path_POSCAR=directory+"/POSCAR.bands"; // plain text

  }
  if(!found_POSCAR) return FALSE;

  // GET THE DATA
  xoutcar.GetProperties(file_OUTCAR);
  xstructure xstr(path_POSCAR,IOVASP_POSCAR);

  //CO20171002 - using tolerance from symmetry calc - START
  double tol;
  if(xstr.CalculateSymmetry()){tol=xstr.sym_eps;}
  else {tol=SYM::defaultTolerance(xstr);}
  //CO20171002 - using tolerance from symmetry calc - STOP

  xoutcar.GetBandGap(tol);
  xdoscar.GetProperties(file_DOSCAR);
  xeigenval.GetProperties(file_EIGENVAL);
  // EFFECTIVE MASSES
  EM_TAG = GetEffectiveMass(xoutcar,xdoscar,xeigenval, xstr);
  if(!EM_TAG) {
    oss << "System              : " << xoutcar.SYSTEM << endl; 
    //[CO20200404 - OBSOLETE]oss << xoutcar.ERROR << "   filename[=" << xoutcar.filename << "]" << endl;
    oss << string("FAILED") << "   filename[=" << xoutcar.filename << "]" << endl;
    oss << endl;
    return FALSE;
  }
  // METALLIC CHANNELS
  int ispin = xeigenval.spin+1;
  if(ispin == 1) {
    if(xoutcar.Egap.at(0) > _METALGAP_) SPIN_UP = TRUE;
  }
  else if(ispin == 2) {
    if(xoutcar.Egap.at(0) > _METALGAP_) SPIN_UP = TRUE;
    if(xoutcar.Egap.at(1) > _METALGAP_) SPIN_DN = TRUE;
  }
  // GET THE OUTPUT
  if(SPIN_UP or SPIN_DN) {
    oss << "System                 :  " << xoutcar.SYSTEM << endl; 
    oss << "Number of records      :  " << xoutcar.carrier_type.size() << endl;
    for(uint itr0=0; itr0<xoutcar.carrier_type.size(); itr0++) {
      oss << std::noshowpos;
      oss << "** Record  " << itr0+1 << endl;
      oss << "   Band index          :  " << xoutcar.band_index.at(itr0)   << endl;
      oss << "   Carrier type        :  " << xoutcar.carrier_type.at(itr0) << endl;
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
      oss << "   Equivalent valleys  :  " << xoutcar.equivalent_valley.at(itr0)   << endl;
      oss << std::showpos;
      oss << "   DOS  eff. mass.     : "  << xoutcar.effective_mass_DOS.at(itr0)  << endl; 
      oss << "   COND eff. mass.     : "  << xoutcar.effective_mass_COND.at(itr0) << endl;
    }
    oss << "** Carrier Masses " << endl;
    oss << "   DOS  elec eff. mass : " 
      << setw(14) << std::left << std::scientific << std::showpos << xoutcar.mass_elec_dos.at(0) << endl;
    oss << "   DOS  hole eff. mass : " 
      << setw(14) << std::left << std::scientific << std::showpos << xoutcar.mass_hole_dos.at(0) << endl;
    oss << "   COND elec eff. mass : " 
      << setw(14) << std::left << std::scientific << std::showpos << xoutcar.mass_elec_conduction.at(0) << endl;
    oss << "   COND hole eff. mass : " 
      << setw(14) << std::left << std::scientific << std::showpos << xoutcar.mass_hole_conduction.at(0) << endl;
    oss << endl;
  }
  return TRUE;
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// PrintEigCurv
// Subroutine that obtains band gap extrema curvatures of the band structure.
// This is then used to approximate the location of the e.m. ellipsoids in the IBZ
// 2015: Camilo E. Calderon
bool PrintEigCurv(string& directory, ostream &oss) {
  stringstream message;

  vector<vector<vector<vector<vector<double> > > > > branches_bnds ;
  vector<vector<vector<xvector<double> > > > branches_kpts ;
  vector<vector<xvector<int> > > branches ;
  vector<vector<vector<int> > > branches_indx ;
  vector<xvector<double> > special_kpts, unique_kpts, unique_kpts_EIG, vkpoint_eig ;
  vector<xvector<int> > repeated_edges, repeat_kpts_EIG ;
  vector<xvector<int> > connect_kpts, connect_kpts_EIG, vrtx_path, ndx_edges ;
  vector<int> connect_kpts_num, repeat_kpts_num ; 
  stringstream file_EIGENVAL, file_KPOINTS, file_OUTCAR ;
  string path_POSCAR;
  double _ZERO=0.0 ;
  int GRIDS ;
  xEIGENVAL xeigenval ;
  xOUTCAR   xoutcar ;
  char LastChar = *directory.rbegin();
  if(LastChar == '/') directory.erase(directory.size()-1);

  deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");

  // EIGENVAL
  bool found_EIGENVAL=FALSE;
  for(uint iext=0;iext<vext.size();iext++) { 
    if(!found_EIGENVAL && aurostd::FileExist(directory+"/EIGENVAL.bands"+vext.at(iext))) {
      found_EIGENVAL=TRUE;
      aurostd::efile2stringstream(directory+"/EIGENVAL.bands"+vext.at(iext), file_EIGENVAL);  // .EXT
    }
  }
  if(!found_EIGENVAL && aurostd::FileExist(directory+"/EIGENVAL.bands")) {
    found_EIGENVAL=TRUE;
    aurostd::file2stringstream(directory+"/EIGENVAL.bands", file_EIGENVAL); // plain text
  }
  if(!found_EIGENVAL) return FALSE;

  // KPOINTS
  bool found_KPOINTS=FALSE;
  for(uint iext=0;iext<vext.size();iext++) { 
    if(!found_KPOINTS && aurostd::FileExist(directory+"/KPOINTS.bands"+vext.at(iext))) {
      found_KPOINTS=TRUE;
      aurostd::efile2stringstream(directory+"/KPOINTS.bands"+vext.at(iext), file_KPOINTS);  // .EXT
    }
  }
  if(!found_KPOINTS && aurostd::FileExist(directory+"/KPOINTS.bands")) {
    found_KPOINTS=TRUE;
    aurostd::file2stringstream(directory+"/KPOINTS.bands", file_KPOINTS); // plain text
  }
  if(!found_KPOINTS) return FALSE;

  // OUTCAR
  bool found_OUTCAR=FALSE;
  for(uint iext=0;iext<vext.size();iext++) { 
    if(!found_OUTCAR && aurostd::FileExist(directory+"/OUTCAR.bands"+vext.at(iext))) {
      found_OUTCAR=TRUE;
      aurostd::efile2stringstream(directory+"/OUTCAR.bands"+vext.at(iext), file_OUTCAR);  // .EXT
    }
  }
  if(!found_OUTCAR && aurostd::FileExist(directory+"/OUTCAR.bands")) {
    found_OUTCAR=TRUE;
    aurostd::file2stringstream(directory+"/OUTCAR.bands", file_OUTCAR); // plain text
  }
  if(!found_OUTCAR) return FALSE;

  // POSCAR
  bool found_POSCAR=FALSE;
  for(uint iext=0;iext<vext.size();iext++) { 
    if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR.bands"+vext.at(iext))) {
      found_POSCAR=TRUE;
      path_POSCAR=directory+"/POSCAR.bands"+vext.at(iext); // EXR
    }
  }
  if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR.bands")) {
    found_POSCAR=TRUE;
    path_POSCAR=directory+"/POSCAR.bands"; // plain text

  }
  if(!found_POSCAR) return FALSE;

  // now have all files

  xeigenval.GetProperties(file_EIGENVAL) ;
  xoutcar.GetProperties(file_OUTCAR) ;

  xstructure xstr(path_POSCAR,IOVASP_POSCAR);

  //CO20171002 - using tolerance from symmetry calc - START
  double tol;
  if(xstr.CalculateSymmetry()){tol=xstr.sym_eps;}
  else {tol=SYM::defaultTolerance(xstr);}
  //CO20171002 - using tolerance from symmetry calc - STOP

  xoutcar.GetBandGap(tol) ;
  if(xoutcar.Egap_net >= _ZERO) {
    ParseKPOINTS      (file_KPOINTS,GRIDS,special_kpts,unique_kpts,repeat_kpts_num) ;
    AdjacencyList_KPT (special_kpts,unique_kpts,connect_kpts,connect_kpts_num) ;
    AdjacencyList_EIG (unique_kpts,connect_kpts,connect_kpts_num,xeigenval,unique_kpts_EIG,connect_kpts_EIG,vkpoint_eig) ;
    RepeatsList       (unique_kpts_EIG,repeat_kpts_num,vkpoint_eig,repeat_kpts_EIG) ;
    VertexPaths       (repeat_kpts_EIG,connect_kpts_EIG,repeat_kpts_num,GRIDS,vrtx_path) ;
    RepeatedEdges     (vrtx_path,repeat_kpts_EIG,repeat_kpts_num,ndx_edges) ;
    VertexBranches    (ndx_edges,repeat_kpts_num,repeat_kpts_EIG,branches) ;
    PathDataStuct     (xeigenval,vkpoint_eig,branches,branches_indx,branches_kpts,branches_bnds) ;
    IBZextrema        (xeigenval,vkpoint_eig,branches) ;
  } else {
    oss << endl ;
    //[CO20200404 - OBSOLETE]xoutcar.ERROR = "Material is metallic." ;
    message << "Material is metallic";  //CO20200404
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, oss, _LOGGER_WARNING_); //CO20200404
    return FALSE ;
  }
  return TRUE ;
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
bool ParseKPOINTS(stringstream& file_KPOINTS,
    int& GRIDS,
    vector<xvector<double> >& special_kpts,
    vector<xvector<double> >& unique_kpts,
    vector<int>& repeat_kpts_num) {

  vector<string> StringKpts, tokens;
  aurostd::string2vectorstring(file_KPOINTS.str(),StringKpts) ;
  string line0 = StringKpts.at(2) ;
  aurostd::string2tokens(line0,tokens) ;
  if(tokens.at(0).at(0) != 'L') {
    cout << "ParseKPOINTS - KPOINTS file is not in Line-mode: " << endl ;
    return FALSE ;
  } else {
    aurostd::string2tokens(StringKpts.at(1),tokens) ;
    GRIDS = aurostd::string2utype<int>(tokens.at(0)) ;
    uint itr1 = 0 ;
    if(StringKpts.size() > 0) {
      for(uint itr0=4; itr0<StringKpts.size(); itr0++) {
        xvector<double> tempvec(4) ;
        string line1 = StringKpts.at(itr0) ;
        aurostd::string2tokens(line1,tokens) ;
        if(tokens.size()>0) {
          special_kpts.push_back(tempvec) ;
          special_kpts.at(itr1)[1] = aurostd::string2utype<double>(tokens.at(0)) ;
          special_kpts.at(itr1)[2] = aurostd::string2utype<double>(tokens.at(1)) ;
          special_kpts.at(itr1)[3] = aurostd::string2utype<double>(tokens.at(2)) ;
          special_kpts.at(itr1)[4] = (double)itr1 ;
          itr1++ ;
        }
      }
    } else {
      cout << "ParseKPOINTS - No strings found in the KPOINTS file: " << endl ;
      return FALSE ;
    }
  }
  vector<xvector<double> > tmp_special_kpts ;
  for(uint itr0=0;itr0<special_kpts.size(); itr0++) {
    xvector<double> tempvec(4) ;
    tmp_special_kpts.push_back(tempvec) ;
    tmp_special_kpts.at(itr0)[1] = special_kpts.at(itr0)[1] ;
    tmp_special_kpts.at(itr0)[2] = special_kpts.at(itr0)[2] ;
    tmp_special_kpts.at(itr0)[3] = special_kpts.at(itr0)[3] ;
    tmp_special_kpts.at(itr0)[4] = special_kpts.at(itr0)[4] ;
  }
  uint itr2=0 ;
  for(uint itr0=0; itr0<special_kpts.size(); itr0++) {
    if(itr2 >= special_kpts.size()) {
      break ;
    }
    uint delcnt=0 ;
    unique_kpts.push_back(tmp_special_kpts.at(itr0)) ;
    tmp_special_kpts.erase(tmp_special_kpts.begin() + itr0) ;
    delcnt+=1 ;
    for(uint itr1=0; itr1<tmp_special_kpts.size(); itr1++) {
      int comp = 0 ;
      for(int itr2=1; itr2<=3; itr2++) {
        bool MATCH = FALSE ;
        CompareDoublesChar(MATCH, unique_kpts.back()[itr2], tmp_special_kpts.at(itr1)[itr2]) ;
        if(MATCH) {
          comp++ ;
        } else if(!MATCH) {
          break ;
        }
      }
      if(comp == 3) {
        tmp_special_kpts.erase(tmp_special_kpts.begin() + itr1) ;
        itr1-- ;
        delcnt++ ;
      } else {
        continue ;
      }
    }
    itr0-- ;
    itr2+=delcnt ;
  }
  for(uint itr0=0; itr0<unique_kpts.size(); itr0++) {
    xvector<int> tempvec(special_kpts.size()) ;
    int index1  = unique_kpts.at(itr0)[4] ;
    int count = 1 ;
    for(uint itr1=1; itr1<=unique_kpts.size(); itr1++) {
      tempvec[itr1] = 99999 ;
    }
    tempvec[count] = (int)unique_kpts.at(itr0)[4] ;
    count++ ;
    for(uint itr1=0; itr1<special_kpts.size(); itr1++) {
      int index2  = special_kpts.at(itr1)[4] ;
      int compare = 0 ;
      for(int itr2=1; itr2<=3; itr2++) {
        bool MATCH = FALSE ;
        CompareDoublesChar(MATCH, unique_kpts.at(itr0)[itr2], special_kpts.at(itr1)[itr2]) ;
        if(MATCH) {
          compare++ ;
        } else if(!MATCH) {
          break ;
        }
      }
      if(compare == 3 and (index1 != index2)) {
        tempvec[count] = index2 ;
        count++ ;
      }
    }
    xvector<int> kptrepeat(count-1) ;
    repeat_kpts_num.push_back(count-1) ;
  }
  return TRUE ;
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
bool AdjacencyList_KPT(vector<xvector<double> >& special_kpts,
    vector<xvector<double> >& unique_kpts,
    vector<xvector<int> >& connect_kpts,
    vector<int>& connect_kpts_num) {

  for(uint itr0=0; itr0<unique_kpts.size(); itr0++) {
    int count = 1 ;
    xvector<int> tempvec0(unique_kpts.size()+1) ;
    for(uint itr1=1; itr1<=unique_kpts.size()+1; itr1++) {
      tempvec0[itr1] = 99999 ;
    }
    tempvec0[count] = (int)unique_kpts.at(itr0)[4] ;
    count++ ;
    for(uint itr1=0; itr1<special_kpts.size(); itr1++) {
      int compare0 = 0 ;
      for(int itr2=1; itr2<=3; itr2++) {
        bool MATCH = FALSE ;
        CompareDoublesChar(MATCH, unique_kpts.at(itr0)[itr2], special_kpts.at(itr1)[itr2]) ;
        if(MATCH) {
          compare0++ ;
        } else if(!MATCH) {
          break ;
        }
      }
      if(compare0 == 3) {
        if((int)special_kpts.at(itr1)[4] % 2 == 0) {
          for(uint itr2=0; itr2<unique_kpts.size(); itr2++) {
            int compare1 = 0 ;
            for(int itr3=1; itr3<=3; itr3++) {
              bool MATCH = FALSE ;
              CompareDoublesChar(MATCH, unique_kpts.at(itr2)[itr3], special_kpts.at(itr1+1)[itr3]) ;
              if(MATCH) {
                compare1++ ;
              } else if(!MATCH) {
                break ;
              }
            }
            if(compare1 == 3) {
              if(count <= (int)unique_kpts.size()+1) {
                tempvec0[count] = (int)unique_kpts.at(itr2)[4] ;
                count++ ;
              } else {
                cout << "GetKPOINTSAdjacencyList - Wrong symmetry" << endl ;
                return FALSE ;
              }
            }
          }
        } else if((int)special_kpts.at(itr1)[4] % 2 == 1) {
          for(uint itr2=0; itr2<unique_kpts.size(); itr2++) {
            int compare1 = 0 ;
            for(int itr3=1; itr3<=3; itr3++) {
              bool MATCH = FALSE ;
              CompareDoublesChar(MATCH, unique_kpts.at(itr2)[itr3], special_kpts.at(itr1-1)[itr3]) ;
              if(MATCH) {
                compare1++ ;
              } else if(!MATCH) {
                break ;
              }
            }
            if(compare1 == 3) {
              if(count <= (int)unique_kpts.size()+1) {
                tempvec0[count] = (int)unique_kpts.at(itr2)[4] ;
                count++ ;
              } else {
                cout << "GetKPOINTSAdjacencyList - Wrong symmetry" << endl ;
                return FALSE ;
              }
            }
          }
        }
      }
    }
    connect_kpts_num.push_back(count-1) ; // possibly count-2?
    xvector<int> kptconnect(count-1) ;
    connect_kpts.push_back(kptconnect) ;
    for(uint itr1=1; itr1<=(uint)count-1; itr1++) {
      if(tempvec0[itr1] < 99999) {
        connect_kpts.at(itr0)[itr1] = tempvec0[itr1] ;
      }
    }
  }
  return TRUE ;
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// AdjacencyList_EIG
// Camilo E. Calderon, 2015
//-------------------------------------------------------------------------------------------------
bool AdjacencyList_EIG(vector<xvector<double> >& unique_kpts,
    vector<xvector<int> >& connect_kpts,
    vector<int>& connect_kpts_num,
    xEIGENVAL& xeigenval,
    vector<xvector<double> >& unique_kpts_EIG,
    vector<xvector<int> >& connect_kpts_EIG,
    vector<xvector<double> >& vkpoint_eig) {

  for(uint itr0=0; itr0<xeigenval.vkpoint.size(); itr0++) {
    xvector<double> tempvec(3) ;
    for(int itr1=1; itr1<=3; itr1++) {
      if(abs(xeigenval.vkpoint.at(itr0)[itr1]) <= 1.0E-15) {
        tempvec[itr1] = 0.0 ;
      } else {
        tempvec[itr1] = xeigenval.vkpoint.at(itr0)[itr1] ;
      }
    }
    vkpoint_eig.push_back(tempvec) ;
  }
  uint start = 0 ;
  for(uint itr0=0; itr0<unique_kpts.size(); itr0++) {
    int templen = connect_kpts_num.at(itr0) ;
    xvector<int> tempvec1(templen) ;
    connect_kpts_EIG.push_back(tempvec1) ;
    for(uint itr1=1; itr1<=(uint)templen; itr1++) {
      connect_kpts_EIG.back()[itr1] = connect_kpts.at(itr0)[itr1] ;
    }
    uint endval = 0 ;
    for(uint itr1=start; itr1<vkpoint_eig.size(); itr1++) {
      int compare = 0 ;
      for(int itr2=1; itr2<=3; itr2++) {
        bool MATCH = FALSE ;
        CompareDoublesChar(MATCH, unique_kpts.at(itr0)[itr2], vkpoint_eig.at(itr1)[itr2]) ;
        if(MATCH) {
          compare++ ;
        } else if(!MATCH) {
          break ;
        }
      }
      if(compare == 3) {
        xvector<double> tempvec2(4) ;
        unique_kpts_EIG.push_back(tempvec2) ;
        unique_kpts_EIG.back()[1] = vkpoint_eig.at(itr1)[1] ;
        unique_kpts_EIG.back()[2] = vkpoint_eig.at(itr1)[2] ;
        unique_kpts_EIG.back()[3] = vkpoint_eig.at(itr1)[3] ;
        unique_kpts_EIG.back()[4] = (double)itr1 ;
        endval = itr1 ;
        break ;
      }
    }
    start = endval+1 ;
  }
  for(uint itr0=0; itr0<connect_kpts.size(); itr0++) {
    for(uint itr1=1; itr1<=(uint)connect_kpts_num.at(itr0); itr1++) {
      for(uint itr2=0; itr2<connect_kpts.size(); itr2++) {
        if(connect_kpts.at(itr2)[1] == connect_kpts.at(itr0)[itr1]) {
          connect_kpts_EIG.at(itr0)[itr1] = (int)unique_kpts_EIG.at(itr2)[4] ;
        }
      }
    }
  }
  return TRUE ;
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
bool RepeatsList(vector<xvector<double> >& unique_kpts_EIG,
    vector<int>& repeat_kpts_num,
    vector<xvector<double> >& vkpoint_eig,
    vector<xvector<int> >& repeat_kpts_EIG) {

  for(uint itr0=0; itr0<unique_kpts_EIG.size(); itr0++) {
    int count =1 ;
    xvector<int> tempvec(repeat_kpts_num.at(itr0)) ;
    repeat_kpts_EIG.push_back(tempvec) ;
    repeat_kpts_EIG.back()[count] = (int)unique_kpts_EIG.at(itr0)[4] ;
    count++ ;
    for(uint itr1=0; itr1<vkpoint_eig.size(); itr1++) {
      int compare = 0 ;
      for(int itr2=1; itr2<=3; itr2++) {
        bool MATCH = FALSE ;
        CompareDoublesChar(MATCH, unique_kpts_EIG.at(itr0)[itr2], vkpoint_eig.at(itr1)[itr2]) ;
        if(MATCH) {
          compare++ ;
        } else if(!MATCH) {
          break ;
        }
      }
      if((compare == 3) and (unique_kpts_EIG.at(itr0)[4] != (int)itr1)) {
        if(count <= repeat_kpts_num.at(itr0)) {
          repeat_kpts_EIG.at(itr0)[count] = (int)itr1 ;
          count++ ;
        } else {
          cout << "KPOINTS.bands does not match EIGENVAL.bands" << endl ;
          return FALSE ;
        }
      }
    }
  }
  return TRUE ;
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
bool VertexPaths(vector<xvector<int> >& repeat_kpts_EIG,
    vector<xvector<int> >& connect_kpts_EIG,
    vector<int>& repeat_kpts_num,
    int& GRIDS,
    // returns
    vector<xvector<int> >& vrtx_path) {

  vector<xvector<int> > vrtx_list ;
  for(uint itr0=0; itr0<connect_kpts_EIG.size(); itr0++) {
    uint templen = repeat_kpts_num.at(itr0)+1 ;
    if(templen >= 3) {
      for(uint itr1=2; itr1<=templen-1; itr1++) {
        for(uint itr2=2; itr2<=templen; itr2++) {
          if(itr1 != itr2 and itr2 > itr1) {
            xvector<int> tempvec(3) ;
            vrtx_list.push_back(tempvec) ;
            vrtx_list.back()[1] = connect_kpts_EIG.at(itr0)[itr1] ;
            vrtx_list.back()[2] = connect_kpts_EIG.at(itr0)[1] ;
            vrtx_list.back()[3] = connect_kpts_EIG.at(itr0)[itr2] ;
          }
        }
      }
    }
  }
  vector<xvector<int> > tmp_vrtx_path ;
  for(uint itr0=0; itr0<vrtx_list.size(); itr0++) {
    xvector<int> vrtx_segments(4) ;
    if(abs(vrtx_list.at(itr0)[1]-vrtx_list.at(itr0)[2]) == GRIDS-1) {
      vrtx_segments[1] = vrtx_list.at(itr0)[1] ;
      vrtx_segments[2] = vrtx_list.at(itr0)[2] ;
    } else {
      for(uint itr1=0; itr1<repeat_kpts_EIG.size(); itr1++) {
        if(repeat_kpts_EIG.at(itr1)[1] == vrtx_list.at(itr0)[1]) {
          for(uint itr2=0; itr2<repeat_kpts_EIG.size(); itr2++) {
            if(repeat_kpts_EIG.at(itr2)[1] == vrtx_list.at(itr0)[2]) {
              for(uint itr3=1; itr3<=(uint)repeat_kpts_num.at(itr1); itr3++) {
                for(uint itr4=1; itr4<=(uint)repeat_kpts_num.at(itr2); itr4++) {
                  if(abs(repeat_kpts_EIG.at(itr1)[itr3]-repeat_kpts_EIG.at(itr2)[itr4]) == GRIDS-1) {
                    vrtx_segments[1] = repeat_kpts_EIG.at(itr1)[itr3] ;
                    vrtx_segments[2] = repeat_kpts_EIG.at(itr2)[itr4] ;
                  }
                }
              }
            }
          }
        }
      }
    }
    if(abs(vrtx_list.at(itr0)[2]-vrtx_list.at(itr0)[3]) == GRIDS-1) {
      vrtx_segments[3] = vrtx_list.at(itr0)[2] ;
      vrtx_segments[4] = vrtx_list.at(itr0)[3] ;
    } else {
      for(uint itr1=0; itr1<repeat_kpts_EIG.size(); itr1++) {
        if(repeat_kpts_EIG.at(itr1)[1] == vrtx_list.at(itr0)[2]) {
          for(uint itr2=0; itr2<repeat_kpts_EIG.size(); itr2++) {
            if(repeat_kpts_EIG.at(itr2)[1] == vrtx_list.at(itr0)[3]) {
              for(uint itr3=1; itr3<=(uint)repeat_kpts_num.at(itr1); itr3++) {
                for(uint itr4=1; itr4<=(uint)repeat_kpts_num.at(itr2); itr4++) {
                  if(abs(repeat_kpts_EIG.at(itr1)[itr3]-repeat_kpts_EIG.at(itr2)[itr4]) == GRIDS-1) {
                    vrtx_segments[3] = repeat_kpts_EIG.at(itr1)[itr3] ;
                    vrtx_segments[4] = repeat_kpts_EIG.at(itr2)[itr4] ;
                  }
                }
              }
            }
          }
        }
      }
    }
    vrtx_path.push_back(vrtx_segments) ;
  }
  return TRUE ;
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
bool RepeatedEdges(vector<xvector<int> >& vrtx_path,
    vector<xvector<int> >& repeat_kpts_EIG,
    vector<int>& repeat_kpts_num,
    // returns:
    vector<xvector<int> >& ndx_edges) {

  vector<vector<xvector<int> > > allpairs ;
  for(uint itr0=0; itr0<repeat_kpts_EIG.size()-1; itr0++) {
    for(uint itr1=itr0+1; itr1<repeat_kpts_EIG.size(); itr1++) {
      vector<xvector<int> > temparray ;
      for(int itr2=1; itr2<=repeat_kpts_num.at(itr0); itr2++) {
        for(int itr3=1; itr3<=repeat_kpts_num.at(itr1); itr3++) {
          xvector<int> tempvec(2) ;
          tempvec[1] = repeat_kpts_EIG.at(itr0)[itr2] ;
          tempvec[2] = repeat_kpts_EIG.at(itr1)[itr3] ;
          temparray.push_back(tempvec) ;
        }
      }
      allpairs.push_back(temparray) ;
    }
  }
  vector<vector<xvector<int> > > ndx_edges_tmp ;
  for(uint itr0=0; itr0<vrtx_path.size(); itr0++) {
    xvector<int> edge1(2), edge2(2), edge_type(2) ;
    int count=0 ;
    edge1[1] = vrtx_path.at(itr0)[1] ;
    edge1[2] = vrtx_path.at(itr0)[2] ;
    edge2[1] = vrtx_path.at(itr0)[3] ;
    edge2[2] = vrtx_path.at(itr0)[4] ;
    for(uint itr1=0; itr1<allpairs.size(); itr1++) {
      for(uint itr2=0; itr2<allpairs.at(itr1).size(); itr2++) {
        xvector<int> pair(2) ;
        pair[1] = allpairs.at(itr1).at(itr2)[1] ;
        pair[2] = allpairs.at(itr1).at(itr2)[2] ;
        if((edge1[1] == pair[1] and edge1[2] == pair[2]) or
            (edge1[1] == pair[2] and edge1[2] == pair[1])) {
          edge_type[1] = itr1 ;
          count++ ;
        }
        if((edge2[1] == pair[1] and edge2[2] == pair[2]) or
            (edge2[1] == pair[2] and edge2[2] == pair[1])) {
          edge_type[2] = itr1 ;
          count++ ;
        }
      }
    }
    if(count == 2) {
      vector<xvector<int> > edge_cur ;
      xvector<int> temp1(3) ;
      temp1[1] = edge_type[1] ;
      temp1[2] = edge1[1] ;
      temp1[3] = edge1[2] ;
      edge_cur.push_back(temp1) ;
      xvector<int> temp2(3) ;
      temp2[1] = edge_type[2] ;
      temp2[2] = edge2[1] ;
      temp2[3] = edge2[2] ;
      edge_cur.push_back(temp2) ;
      ndx_edges_tmp.push_back(edge_cur) ;
    }
  }
  xvector<int> ndx_edges0(3) ;
  ndx_edges0 = ndx_edges_tmp.at(0).at(0) ;
  ndx_edges.push_back(ndx_edges0) ;
  for(uint itr0=0; itr0<ndx_edges_tmp.size(); itr0++) {
    for(uint itr1=0; itr1<ndx_edges_tmp.at(itr0).size(); itr1++) {
      bool EDGE = FALSE ;
      for(uint itr2=0; itr2<ndx_edges.size(); itr2++) {
        if(ndx_edges.at(itr2)[1] == ndx_edges_tmp.at(itr0).at(itr1)[1]) {
          EDGE = TRUE ;
          break ;
        }
      }
      if(!EDGE) {
        ndx_edges.push_back(ndx_edges_tmp.at(itr0).at(itr1)) ;
      }
    }
  }
  ndx_edges_tmp.clear() ;
  return TRUE ;
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
bool VertexBranches(vector<xvector<int> >& ndx_edges,
    vector<int>& repeat_kpts_num,
    vector<xvector<int> >& repeat_kpts_EIG,
    // returns:
    vector<vector<xvector<int> > >& branches) {

  vector<int> ndx_edges_row ;
  bool NEWEDGE=FALSE ; // CAMILOFIX
  for(uint itr0=0; itr0<=ndx_edges.size(); itr0++) { // NOT A TYPO!!!
    vector<xvector<int> > vertex_edges ;
    //bool NEWEDGE ; // CAMILOFIX
    int repeat_kpts_row ;
    if(ndx_edges_row.size() == 0) {
      for(uint itr1=0; itr1<repeat_kpts_EIG.size(); itr1++) {
        for(uint itr2=1; itr2<=(uint)repeat_kpts_num.at(itr1); itr2++) {
          if(ndx_edges.at(itr0)[2] == repeat_kpts_EIG.at(itr1)[itr2]) {
            xvector<int> edge_crnt(2) ;
            if(ndx_edges.at(itr0)[2] == repeat_kpts_EIG.at(itr1)[itr2]) {
              edge_crnt[1] = ndx_edges.at(itr0)[2] ;
              edge_crnt[2] = ndx_edges.at(itr0)[3] ;
            } else if(ndx_edges.at(itr0)[2] != repeat_kpts_EIG.at(itr1)[itr2]) {
              edge_crnt[1] = ndx_edges.at(itr0)[3] ;
              edge_crnt[2] = ndx_edges.at(itr0)[2] ;
            }
            ndx_edges_row.push_back(ndx_edges.at(itr0)[1]) ;
            vertex_edges.push_back(edge_crnt) ;
            repeat_kpts_row = itr1 ;
            NEWEDGE = TRUE ;
            goto EDGE2VERTICES ;
          }
        }
      }
    }
    else if(ndx_edges_row.size() > 0) {
      for(uint itr1=0; itr1<ndx_edges.size(); itr1++) {
        xvector<int> edge_crnt(2) ;
        bool MATCH1 = FALSE ;
        for(uint itr2=0; itr2<ndx_edges_row.size(); itr2++) {
          if(ndx_edges.at(itr1)[1] != ndx_edges_row.at(itr2)) {
            edge_crnt[1] = ndx_edges.at(itr1)[2] ;
            edge_crnt[2] = ndx_edges.at(itr1)[3] ;
            CompareEdges(branches, vertex_edges, edge_crnt, MATCH1) ;
            if(!MATCH1) {
              for(uint itr3=0; itr3<repeat_kpts_EIG.size(); itr3++) {
                for(uint itr4=1; itr4<=(uint)repeat_kpts_num.at(itr3); itr4++) {
                  if(edge_crnt[1] == repeat_kpts_EIG.at(itr3)[itr4]) {
                    ndx_edges_row.push_back(ndx_edges.at(itr1)[1]) ;
                    vertex_edges.push_back(edge_crnt) ;
                    repeat_kpts_row = itr3 ;
                    NEWEDGE = TRUE ;
                    goto EDGE2VERTICES ;
                  }
                }
              }
            } else if(MATCH1) {
              xvector<int> edge_crnt_new(2) ;
              bool MATCH2 = FALSE ;
              edge_crnt_new[1] = ndx_edges.at(itr1)[3] ;
              edge_crnt_new[2] = ndx_edges.at(itr1)[2] ;
              CompareEdges(branches, vertex_edges, edge_crnt_new, MATCH2) ;
              if(!MATCH2) {
                for(uint itr3=0; itr3<repeat_kpts_EIG.size(); itr3++) {
                  for(uint itr4=1; itr4<=(uint)repeat_kpts_num.at(itr3); itr4++) {
                    if(edge_crnt_new[1] == repeat_kpts_EIG.at(itr3)[itr4]) {
                      ndx_edges_row.push_back(ndx_edges.at(itr1)[1]) ;
                      vertex_edges.push_back(edge_crnt_new) ;
                      repeat_kpts_row = itr3 ;
                      NEWEDGE = TRUE ;
                      goto EDGE2VERTICES ;
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
    if(NEWEDGE) {
      for(uint itr1=0; itr1<ndx_edges.size(); itr1++) {
        for(uint itr2=0; itr2<ndx_edges_row.size(); itr2++) {
          if( ndx_edges_row.at(itr2) != ndx_edges.at(itr1)[1] ) {
            for(int itr3=2; itr3<=3; itr3++) {
              for(int itr4=1; itr4<=repeat_kpts_num.at(repeat_kpts_row); itr4++) {
                if(repeat_kpts_EIG.at(repeat_kpts_row)[itr4] == ndx_edges.at(itr1)[itr3]) {
                  xvector<int> edge_crnt(2) ;
                  if(repeat_kpts_EIG.at(repeat_kpts_row)[itr4] == ndx_edges.at(itr1)[2]) {
                    edge_crnt[1] = ndx_edges.at(itr1)[2] ;
                    edge_crnt[2] = ndx_edges.at(itr1)[3] ;
                  } else if(repeat_kpts_EIG.at(repeat_kpts_row)[itr4] == ndx_edges.at(itr1)[3]) {
                    edge_crnt[1] = ndx_edges.at(itr1)[3] ;
                    edge_crnt[2] = ndx_edges.at(itr1)[2] ;
                  }
                  bool MATCH3 = FALSE ;
                  CompareEdges(branches, vertex_edges, edge_crnt, MATCH3) ;
                  if(!MATCH3) {
                    vertex_edges.push_back(edge_crnt) ;
                  }
                  NEWEDGE = FALSE ;
                }
              }
            }
          } else if( ndx_edges_row.at(itr2) == ndx_edges.at(itr1)[1] ) {
            for(int itr3=2; itr3<=3; itr3++) {
              for(int itr4=1; itr4<=repeat_kpts_num.at(repeat_kpts_row); itr4++) {
                if(repeat_kpts_EIG.at(repeat_kpts_row)[itr4] == ndx_edges.at(itr1)[itr3]) {
                  xvector<int> edge_crnt(2) ;
                  bool MATCH3 = FALSE ;
                  edge_crnt[1] = ndx_edges.at(itr1)[3] ;
                  edge_crnt[2] = ndx_edges.at(itr1)[2] ;
                  CompareEdges(branches, vertex_edges, edge_crnt, MATCH3) ;
                  if(!MATCH3) {
                    for(uint itr5=1; itr5<=(uint)repeat_kpts_num.at(repeat_kpts_row); itr5++) {
                      if(repeat_kpts_EIG.at(repeat_kpts_row)[itr5] == edge_crnt[1]) {
                        vertex_edges.push_back(edge_crnt) ;
                      }
                    }
                  }
                  NEWEDGE = FALSE ;
                }
              }
            }
          }
        }
      }
      branches.push_back(vertex_edges) ;
    }
  }
  return TRUE ;
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
    vector<xvector<double> >& vkpoint_eig,
    vector<vector<xvector<int> > >& branches,
    // returns:
    vector<vector<vector<int> > >& branches_indx,
    vector<vector<vector<xvector<double> > > >& branches_kpts,
    vector<vector<vector<vector<vector<double> > > > >& branches_bnds) {

  for(uint itr0=0; itr0<branches.size(); itr0++) {
    vector<vector<vector<vector<double> > > > ener_edge ;
    vector<vector<xvector<double> > > kpts_edge ;
    vector<vector<int> > indx_edge ;
    for(uint itr1=0; itr1<branches.at(itr0).size(); itr1++) {
      if(branches.at(itr0).at(itr1)[1] > branches.at(itr0).at(itr1)[2]) {
        vector<vector<vector<double> > > ener_list ;
        vector<xvector<double> > kpts_list ;
        vector<int> indx_list ;
        for(int itr2=branches.at(itr0).at(itr1)[1]; itr2>=branches.at(itr0).at(itr1)[2]; itr2--) {
          xvector<double> kpts_cart = vkpoint_eig.at(itr2) ;
          vector<vector<double> > ener_band ;
          kpts_list.push_back(kpts_cart) ;
          indx_list.push_back(itr2) ;
          for(uint itr3=0; itr3<xeigenval.number_bands; itr3++) {
            vector<double> ener_spin ;
            for(uint itr4=0; itr4<(uint)(xeigenval.spin+1); itr4++) {
              ener_spin.push_back(xeigenval.venergy.at(itr2).at(itr3).at(itr4)) ;
            }
            ener_band.push_back(ener_spin) ;
          }
          ener_list.push_back(ener_band) ;
        }
        ener_edge.push_back(ener_list) ;
        indx_edge.push_back(indx_list) ;
        kpts_edge.push_back(kpts_list) ;
      } else if(branches.at(itr0).at(itr1)[1] < branches.at(itr0).at(itr1)[2]) {
        vector<vector<vector<double> > > ener_list ;
        vector<xvector<double> > kpts_list ;
        vector<int> indx_list ;
        for(int itr2=branches.at(itr0).at(itr1)[1]; itr2<=branches.at(itr0).at(itr1)[2]; itr2++) {
          xvector<double> kpts_cart = vkpoint_eig.at(itr2) ;
          vector<vector<double> > ener_band ;
          kpts_list.push_back(kpts_cart) ;
          indx_list.push_back(itr2) ;
          for(uint itr3=0; itr3<xeigenval.number_bands; itr3++) {
            vector<double> ener_spin ;
            for(uint itr4=0; itr4<(uint)(xeigenval.spin+1); itr4++) {
              ener_spin.push_back(xeigenval.venergy.at(itr2).at(itr3).at(itr4)) ;
            }
            ener_band.push_back(ener_spin) ;
          }
          ener_list.push_back(ener_band) ;
        }
        ener_edge.push_back(ener_list) ;
        indx_edge.push_back(indx_list) ;
        kpts_edge.push_back(kpts_list) ;
      }
    }
    branches_bnds.push_back(ener_edge) ;
    branches_indx.push_back(indx_edge) ;
    branches_kpts.push_back(kpts_edge) ;
  }
  return TRUE ;
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
bool IBZextrema(xEIGENVAL& xeigenval,
    vector<xvector<double> >& vkpoint_eig,
    vector<vector<xvector<int> > >& branches) {

  vector<double> all_curves ;
  vector<double> curvature ;
  for(uint itr0=0; itr0<branches.size(); itr0++) {
    for(uint itr1=0; itr1<branches.at(itr0).size(); itr1++) {
      vector<xvector<int> > beg_edge_ndx, end_edge_ndx ;
      for(uint itr2=0; itr2<branches.size(); itr2++) {
        for(uint itr3=0; itr3<branches.at(itr2).size(); itr3++) {
          if(branches.at(itr2).at(itr3)[1] == branches.at(itr0).at(itr1)[1]) {
            for(uint itr4=0; itr4<branches.at(itr2).size(); itr4++) {
              if(itr3 != itr4) {
                beg_edge_ndx.push_back(branches.at(itr2).at(itr4)) ;
              }
            }
          }
          if(branches.at(itr2).at(itr3)[1] == branches.at(itr0).at(itr1)[2]) {
            for(uint itr4=0; itr4<branches.at(itr2).size(); itr4++) {
              if(itr3 != itr4) {
                end_edge_ndx.push_back(branches.at(itr2).at(itr4)) ;
              }
            }
          }
        }
      }
      // edge curvatures
      vector<xvector<double> > posvec ;
      // for(uint itr3=0; itr3<xeigenval.number_bands; itr3++)
      for(uint itr3=0; itr3<1; itr3++)
      {   //CO20200106 - patching for auto-indenting
        // for(uint itr4=0; itr4<xeigenval.spin; itr4++)
        for(uint itr4=0; itr4<1; itr4++)
        { //CO20200106 - patching for auto-indenting
          xvector<double> eigvec(5) ;
          xvector<int> ndxvec(5) ;
          if(branches.at(itr0).at(itr1)[1] > branches.at(itr0).at(itr1)[2]) {
            for(int itr2=branches.at(itr0).at(itr1)[1]-2; itr2>=branches.at(itr0).at(itr1)[2]+2; itr2--) {
              ndxvec[1] = itr2+2 ;
              ndxvec[2] = itr2+1 ;
              ndxvec[3] = itr2   ;
              ndxvec[4] = itr2-1 ;
              ndxvec[5] = itr2-2 ;
              for(int itr2=1; itr2<=5; itr2++) {
                eigvec[itr2] = xeigenval.venergy.at(ndxvec[itr2]).at(itr3).at(itr4) ;
                posvec.push_back(vkpoint_eig.at(ndxvec[itr2])) ;
              }
              NaiveCurvatures(eigvec,posvec,curvature) ;
              for(uint itr3=0; itr3<curvature.size(); itr3++) all_curves.push_back(curvature.at(itr3)) ;
              curvature.clear() ; posvec.clear() ;
            }
            //ME20200724 - there used to be an exit here without comment or explanation
            string message = "branches.at(itr0).at(itr1)[1] > branches.at(itr0).at(itr1)[2])";
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
          } else if(branches.at(itr0).at(itr1)[1] < branches.at(itr0).at(itr1)[2]) {
            for(int itr2=branches.at(itr0).at(itr1)[1]+2; itr2<=branches.at(itr0).at(itr1)[2]-2; itr2++) {
              ndxvec[1] = itr2-2 ;
              ndxvec[2] = itr2-1 ;
              ndxvec[3] = itr2   ;
              ndxvec[4] = itr2+1 ;
              ndxvec[5] = itr2+2 ;
              for(int itr2=1; itr2<=5; itr2++) {
                eigvec[itr2] = xeigenval.venergy.at(ndxvec[itr2]).at(itr3).at(itr4) ;
                posvec.push_back(vkpoint_eig.at(ndxvec[itr2])) ;
              }
              NaiveCurvatures(eigvec,posvec,curvature) ;
              for(uint itr3=0; itr3<curvature.size(); itr3++) all_curves.push_back(curvature.at(itr3)) ;
              curvature.clear() ; posvec.clear() ;
            }
          }
          // vertex curvatures
          if(beg_edge_ndx.size() != 0) {
            for(uint itr2=0; itr2<beg_edge_ndx.size(); itr2++) {
              xvector<double> eigvec(6) ;
              xvector<int> ndxvec(6) ;
              if(beg_edge_ndx.at(itr2)[1] < beg_edge_ndx.at(itr2)[2]) {
                ndxvec[1] = beg_edge_ndx.at(itr2)[1]+2 ;
                ndxvec[2] = beg_edge_ndx.at(itr2)[1]+1 ;
                ndxvec[3] = branches.at(itr0).at(itr1)[1] ;
                if(branches.at(itr0).at(itr1)[1] > branches.at(itr0).at(itr1)[2]) {
                  ndxvec[4] = branches.at(itr0).at(itr1)[1]-1 ;
                  ndxvec[5] = branches.at(itr0).at(itr1)[1]-2 ;
                  ndxvec[6] = branches.at(itr0).at(itr1)[1]-3 ;
                } else if(branches.at(itr0).at(itr1)[1] < branches.at(itr0).at(itr1)[2]) {
                  ndxvec[4] = branches.at(itr0).at(itr1)[1]+1 ;
                  ndxvec[5] = branches.at(itr0).at(itr1)[1]+2 ;
                  ndxvec[6] = branches.at(itr0).at(itr1)[1]+3 ;
                }
                for(int itr5=1; itr5<=6; itr5++) {
                  eigvec[itr5] = xeigenval.venergy.at(ndxvec[itr5]).at(itr3).at(itr4) ;
                  posvec.push_back(vkpoint_eig.at(ndxvec[itr5])) ;
                }
                NaiveCurvatures(eigvec,posvec,curvature) ;
                for(uint itr3=0; itr3<curvature.size(); itr3++) all_curves.push_back(curvature.at(itr3)) ;
                curvature.clear() ; posvec.clear() ;
              } else if(beg_edge_ndx.at(itr2)[1] > beg_edge_ndx.at(itr2)[2]) {
                ndxvec[1] = beg_edge_ndx.at(itr2)[1]-2 ;
                ndxvec[2] = beg_edge_ndx.at(itr2)[1]-1 ;
                ndxvec[3] = branches.at(itr0).at(itr1)[1] ;
                if(branches.at(itr0).at(itr1)[1] > branches.at(itr0).at(itr1)[2]) {
                  ndxvec[4] = branches.at(itr0).at(itr1)[1]-1 ;
                  ndxvec[5] = branches.at(itr0).at(itr1)[1]-2 ;
                  ndxvec[6] = branches.at(itr0).at(itr1)[1]-3 ;
                } else if(branches.at(itr0).at(itr1)[1] < branches.at(itr0).at(itr1)[2]) {
                  ndxvec[4] = branches.at(itr0).at(itr1)[1]+1 ;
                  ndxvec[5] = branches.at(itr0).at(itr1)[1]+2 ;
                  ndxvec[6] = branches.at(itr0).at(itr1)[1]+3 ;
                }
                for(int itr5=1; itr5<=6; itr5++) {
                  eigvec[itr5] = xeigenval.venergy.at(ndxvec[itr5]).at(itr3).at(itr4) ;
                  posvec.push_back(vkpoint_eig.at(ndxvec[itr5])) ;
                }
                NaiveCurvatures(eigvec,posvec,curvature) ;
                for(uint itr3=0; itr3<curvature.size(); itr3++) all_curves.push_back(curvature.at(itr3)) ;
                curvature.clear() ; posvec.clear() ;
              }
            }
          } else if(beg_edge_ndx.size() == 0) { // if no branches @ beg
            xvector<double> eigvec(6) ;
            xvector<int> ndxvec(6) ;
            if(branches.at(itr0).at(itr1)[1] > branches.at(itr0).at(itr1)[2]) {
              ndxvec[1] = branches.at(itr0).at(itr1)[1]-3 ;
              ndxvec[2] = branches.at(itr0).at(itr1)[1]-2 ;
              ndxvec[3] = branches.at(itr0).at(itr1)[1]-1 ;
              ndxvec[4] = branches.at(itr0).at(itr1)[1] ;
              ndxvec[5] = branches.at(itr0).at(itr1)[1]-1 ;
              ndxvec[6] = branches.at(itr0).at(itr1)[1]-2 ;
              for(int itr5=1; itr5<=6; itr5++) {
                eigvec[itr5] = xeigenval.venergy.at(ndxvec[itr5]).at(itr3).at(itr4) ;
                posvec.push_back(vkpoint_eig.at(ndxvec[itr5])) ;
              }
              NaiveCurvatures(eigvec,posvec,curvature) ;
              for(uint itr3=0; itr3<curvature.size(); itr3++) all_curves.push_back(curvature.at(itr3)) ;
              curvature.clear() ; posvec.clear() ;
            } else if(branches.at(itr0).at(itr1)[1] < branches.at(itr0).at(itr1)[2]) {
              ndxvec[1] = branches.at(itr0).at(itr1)[1]+3 ;
              ndxvec[2] = branches.at(itr0).at(itr1)[1]+2 ;
              ndxvec[3] = branches.at(itr0).at(itr1)[1]+1 ;
              ndxvec[4] = branches.at(itr0).at(itr1)[1] ;
              ndxvec[5] = branches.at(itr0).at(itr1)[1]+1 ;
              ndxvec[6] = branches.at(itr0).at(itr1)[1]+2 ;
              for(int itr5=1; itr5<=6; itr5++) {
                eigvec[itr5] = xeigenval.venergy.at(ndxvec[itr5]).at(itr3).at(itr4) ;
                posvec.push_back(vkpoint_eig.at(ndxvec[itr5])) ;
              }
              NaiveCurvatures(eigvec,posvec,curvature) ;
              for(uint itr3=0; itr3<curvature.size(); itr3++) all_curves.push_back(curvature.at(itr3)) ;
              curvature.clear() ; posvec.clear() ;
            }
          }
          if(end_edge_ndx.size() != 0) {
            for(uint itr2=0; itr2<end_edge_ndx.size(); itr2++) {
              xvector<double> eigvec(6) ;
              xvector<int> ndxvec(6) ;
              if(end_edge_ndx.at(itr2)[1] < end_edge_ndx.at(itr2)[2]) {
                if(branches.at(itr0).at(itr1)[1] > branches.at(itr0).at(itr1)[2]) {
                  ndxvec[1] = branches.at(itr0).at(itr1)[2]+3 ;
                  ndxvec[2] = branches.at(itr0).at(itr1)[2]+2 ;
                  ndxvec[3] = branches.at(itr0).at(itr1)[2]+1 ;
                } else if(branches.at(itr0).at(itr1)[1] < branches.at(itr0).at(itr1)[2]) {
                  ndxvec[1] = branches.at(itr0).at(itr1)[2]-3 ;
                  ndxvec[2] = branches.at(itr0).at(itr1)[2]-2 ;
                  ndxvec[3] = branches.at(itr0).at(itr1)[2]-1 ;
                }
                ndxvec[4] = branches.at(itr0).at(itr1)[2] ;
                ndxvec[5] = end_edge_ndx.at(itr2)[1]+1 ;
                ndxvec[6] = end_edge_ndx.at(itr2)[1]+2 ;
                for(int itr5=1; itr5<=6; itr5++) {
                  eigvec[itr5] = xeigenval.venergy.at(ndxvec[itr5]).at(itr3).at(itr4) ;
                  posvec.push_back(vkpoint_eig.at(ndxvec[itr5])) ;
                }
                NaiveCurvatures(eigvec,posvec,curvature) ;
                for(uint itr3=0; itr3<curvature.size(); itr3++) all_curves.push_back(curvature.at(itr3)) ;
                curvature.clear() ; posvec.clear() ;
              } else if(end_edge_ndx.at(itr2)[1] > end_edge_ndx.at(itr2)[2]) {
                if(branches.at(itr0).at(itr1)[1] > branches.at(itr0).at(itr1)[2]) {
                  ndxvec[1] = branches.at(itr0).at(itr1)[2]+3 ;
                  ndxvec[2] = branches.at(itr0).at(itr1)[2]+2 ;
                  ndxvec[3] = branches.at(itr0).at(itr1)[2]+1 ;
                } else if(branches.at(itr0).at(itr1)[1] < branches.at(itr0).at(itr1)[2]) {
                  ndxvec[1] = branches.at(itr0).at(itr1)[2]-3 ;
                  ndxvec[2] = branches.at(itr0).at(itr1)[2]-2 ;
                  ndxvec[3] = branches.at(itr0).at(itr1)[2]-1 ;
                }
                ndxvec[4] = branches.at(itr0).at(itr1)[2] ;
                ndxvec[5] = end_edge_ndx.at(itr2)[1]-1 ;
                ndxvec[6] = end_edge_ndx.at(itr2)[1]-2 ;
                for(int itr5=1; itr5<=6; itr5++) {
                  eigvec[itr5] = xeigenval.venergy.at(ndxvec[itr5]).at(itr3).at(itr4) ;
                  posvec.push_back(vkpoint_eig.at(ndxvec[itr5])) ;
                }
                NaiveCurvatures(eigvec,posvec,curvature) ;
                for(uint itr3=0; itr3<curvature.size(); itr3++) all_curves.push_back(curvature.at(itr3)) ;
                curvature.clear() ; posvec.clear() ;
              }
            }
          } else if(end_edge_ndx.size() == 0) { // only in disconnected edges
            xvector<double> eigvec(6) ;
            xvector<int> ndxvec(6) ;
            if(branches.at(itr0).at(itr1)[1] > branches.at(itr0).at(itr1)[2]) {
              ndxvec[1] = branches.at(itr0).at(itr1)[1]-3 ;
              ndxvec[2] = branches.at(itr0).at(itr1)[1]-2 ;
              ndxvec[3] = branches.at(itr0).at(itr1)[1]-1 ;
              ndxvec[4] = branches.at(itr0).at(itr1)[1] ;
              ndxvec[5] = branches.at(itr0).at(itr1)[1]-1 ;
              ndxvec[6] = branches.at(itr0).at(itr1)[1]-2 ;
              for(int itr5=1; itr5<=6; itr5++) {
                eigvec[itr5] = xeigenval.venergy.at(ndxvec[itr5]).at(itr3).at(itr4) ;
                posvec.push_back(vkpoint_eig.at(ndxvec[itr5])) ;
              }
              NaiveCurvatures(eigvec,posvec,curvature) ;
              for(uint itr3=0; itr3<curvature.size(); itr3++) all_curves.push_back(curvature.at(itr3)) ;
              curvature.clear() ; posvec.clear() ;
            } else if(branches.at(itr0).at(itr1)[1] < branches.at(itr0).at(itr1)[2]) {
              ndxvec[1] = branches.at(itr0).at(itr1)[1]+3 ;
              ndxvec[2] = branches.at(itr0).at(itr1)[1]+2 ;
              ndxvec[3] = branches.at(itr0).at(itr1)[1]+1 ;
              ndxvec[4] = branches.at(itr0).at(itr1)[1] ;
              ndxvec[5] = branches.at(itr0).at(itr1)[1]+1 ;
              ndxvec[6] = branches.at(itr0).at(itr1)[1]+2 ;
              for(int itr5=1; itr5<=6; itr5++) {
                eigvec[itr5] = xeigenval.venergy.at(ndxvec[itr5]).at(itr3).at(itr4) ;
                posvec.push_back(vkpoint_eig.at(ndxvec[itr5])) ;
              }
              NaiveCurvatures(eigvec,posvec,curvature) ;
              for(uint itr3=0; itr3<curvature.size(); itr3++) all_curves.push_back(curvature.at(itr3)) ;
              curvature.clear() ; posvec.clear() ;
            }
          }
        }
      }
    }
  }
  return TRUE ;
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
//void NaiveCurvatures(xvector<double> eigvec,
void NaiveCurvatures(xvector<double>& eigvec,
    vector<xvector<double> >& posvec,
    // returns:
    vector<double>& curvature) {

  if(posvec.size() == 5) {
    curvature.push_back(StencilLinear1D(posvec,eigvec) ) ;
  } else if(posvec.size() == 6) {
    xvector<double> eigenvals(posvec.size()-1) ;
    vector<xvector<double> > positions ;
    for(uint itr0=0; itr0<posvec.size()-1; itr0++) {
      xvector<double> tempvec(3) ;
      tempvec[1] = posvec.at(itr0)[1] ;
      tempvec[2] = posvec.at(itr0)[2] ;
      tempvec[3] = posvec.at(itr0)[3] ;
      positions.push_back(tempvec) ;
      eigenvals[itr0+1] = eigvec[itr0+1] ;
    }
    curvature.push_back(StencilLinear1D(positions,eigenvals) ) ;
    positions.clear() ;
    for(uint itr0=1; itr0<posvec.size(); itr0++) {
      xvector<double> tempvec(3) ;
      tempvec[1] = posvec.at(itr0)[1] ;
      tempvec[2] = posvec.at(itr0)[2] ;
      tempvec[3] = posvec.at(itr0)[3] ;
      positions.push_back(tempvec) ;
      eigenvals[itr0] = eigvec[itr0+1] ;
    }
    curvature.push_back(StencilLinear1D(positions,eigenvals) ) ;
  }
  return ;
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
double StencilLinear1D(vector<xvector<double> >& positions, xvector<double>& eigenvals) {
  xvector<double> numer(5) ;
  xvector<double> posns(5) ;
  xvector<double> delta(3) ;
  double denom, dist = 0 ;
  posns[1] = 0 ;
  for(uint itr0=1; itr0<positions.size(); itr0++) {
    delta[1] = pow((positions.at(itr0-1)[1]-positions.at(itr0)[1]),2.0) ;
    delta[2] = pow((positions.at(itr0-1)[2]-positions.at(itr0)[2]),2.0) ;
    delta[3] = pow((positions.at(itr0-1)[3]-positions.at(itr0)[3]),2.0) ;
    dist += pow((delta[1]+delta[2]+delta[3]),0.5) ;
    posns[itr0+1] = dist ;
  }
  numer[1] =       -eigenvals[1] ;
  numer[2] =  16.0* eigenvals[2] ;
  numer[3] = -30.0* eigenvals[3] ;
  numer[4] =  16.0* eigenvals[4] ;
  numer[5] =       -eigenvals[5] ;
  denom    = pow((12.0*pow(((posns[5]-posns[1])/(4)),2.0)),-1.0) ;
  return (numer[1]+numer[2]+numer[3]+numer[4]+numer[5])/denom ;
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

  ostringstream oss1 ;
  ostringstream oss2 ;
  oss1.precision(8) ;
  oss2.precision(8) ;
  oss1 << std::scientific << number1 ;
  oss2 << std::scientific << number2 ;
  if(oss1.str().size() == oss2.str().size()) {
    int lencomp = oss1.str().size() ;
    int compare = 0 ;
    for(uint itr3=0; itr3<oss1.str().size(); itr3++) {
      if(oss1.str().at(itr3) == oss2.str().at(itr3)) {
        compare++ ;
      }
    }
    if(compare == lencomp) {
      MATCH = TRUE ;
      return ;
    } else {
      MATCH = FALSE ;
      return ;
    }
    // following is probably not needed, but keep it.
  } else if(oss1.str().size() != oss2.str().size()) {
    int decloc1 = -1, decloc2 = -1 ;
    int len1 = oss1.str().size() ;
    int len2 = oss2.str().size() ;
    for(uint itr0=0; itr0<(uint)min(len1,len2); itr0++) {
      if((int)oss1.str().at(itr0) == 46) {
        decloc1 = itr0 ;
      }
      if((int)oss2.str().at(itr0) == 46) {
        decloc2 = itr0 ;
      }
    }
    if((decloc1 == -1 and decloc2  > -1) or
        (decloc1  > -1 and decloc2 == -1) or
        (decloc1 == -1 and decloc2 == -1)) {
      MATCH = FALSE ;
      return ;
    } else if(decloc1 != decloc2) {
      MATCH = FALSE ;
      return ;
    } else if(decloc1 == decloc2) {
      double min_num=0.0, max_num=0.0 ; // CAMILOFIX
      int min_len = min(len1,len2) - 1 ;
      int max_len = max(len1,len2) - 1 ;
      if(len1 == min_len+1) {
        min_num = number1 ;
        max_num = number2 ;
      } else if(len2 == min_len+1) {
        min_num = number2 ;
        max_num = number1 ;
      }
      double hi_lim = min_num + 4.0 * pow(10,-max_len) ;
      double lo_lim = min_num - 5.0 * pow(10,-max_len) ;
      if((max_num >= lo_lim) and (max_num < hi_lim)) {
        MATCH = TRUE ;
        return ;
      } else {
        MATCH = FALSE ;
        return ;
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
void CompareEdges (vector<vector<xvector<int> > >& branches,
    vector<xvector<int> >& vertex_edges,
    xvector<int>& test_edge,
    bool& COMPARE_EDGES) {

  // branches
  if(branches.size() > 0) {
    for(uint itr0=0; itr0<branches.size(); itr0++) {
      for(uint itr1=0; itr1<branches.at(itr0).size(); itr1++) {
        if((branches.at(itr0).at(itr1)[1] == test_edge[1]) and
            (branches.at(itr0).at(itr1)[2] == test_edge[2])) {
          COMPARE_EDGES = TRUE ;
          return ;
        }
      }
    }
  } else if(branches.size() == 0) {
    COMPARE_EDGES = FALSE ;
    return ;
  }
  // vertex_edges
  if(vertex_edges.size() > 0) {
    for(uint itr0=0; itr0<vertex_edges.size(); itr0++) {
      if((vertex_edges.at(itr0)[1] == test_edge[1]) and
          (vertex_edges.at(itr0)[2] == test_edge[2])) {
        COMPARE_EDGES = TRUE ;
        return ;
      }
    }
  } else if(vertex_edges.size() == 0) {
    COMPARE_EDGES = FALSE ;
    return ;
  }
}
//-------------------------------------------------------------------------------------------------
// ***************************************************************************
// class xPOTCAR
//ME20200427 - included xStream::initialize
xPOTCAR::xPOTCAR(ostream& oss) : xStream(oss) {free();} //CO20200404 - xStream integration for logging
xPOTCAR::xPOTCAR(ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss) {free();} //CO20200404 - xStream integration for logging
xPOTCAR::xPOTCAR(const string& fileIN,bool QUIET,ostream& oss) : xStream(oss) {m_initialized=initialize(fileIN,QUIET);} //CO20200404 - xStream integration for logging
xPOTCAR::xPOTCAR(const string& fileIN,ofstream& FileMESSAGE,bool QUIET,ostream& oss) : xStream(FileMESSAGE,oss) {m_initialized=initialize(fileIN,QUIET);} //CO20200404 - xStream integration for logging

bool xPOTCAR::initialize(const string& fileIN,ostream& oss,bool QUIET) {xStream::initialize(oss);return initialize(fileIN,QUIET);} //CO20200508
bool xPOTCAR::initialize(const string& fileIN,ofstream& FileMESSAGE,ostream& oss,bool QUIET) {xStream::initialize(FileMESSAGE,oss);return initialize(fileIN,QUIET);} //CO20200508
bool xPOTCAR::initialize(const string& fileIN, bool QUIET) {
  free();
  filename = fileIN;
  return GetPropertiesFile(fileIN, QUIET);
}

xPOTCAR::xPOTCAR(const xPOTCAR& b) : xStream(*b.getOFStream(),*b.getOSS()) {free();copy(b);} // copy PUBLIC

xPOTCAR::~xPOTCAR() {xStream::free();free();} //CO20191110 //CO20200404 - xStream integration for logging

void xPOTCAR::free() {
  m_initialized=false; //CO20200404
  content="";                   // for aflowlib_libraries.cpp
  vcontent.clear();             // for aflowlib_libraries.cpp
  filename="";                  // for aflowlib_libraries.cpp
  title="";                     // for aflowlib_libraries.cpp
  POTCAR_PAW=FALSE;             // for aflowlib_libraries.cpp
  POTCAR_TYPE="";               // for aflowlib_libraries.cpp
  POTCAR_KINETIC=FALSE;         // for aflowlib_libraries.cpp
  POTCAR_GW=FALSE;              // for aflowlib_libraries.cpp
  POTCAR_AE=FALSE;              // for aflowlib_libraries.cpp
  // begin shared xPOTCAR
  ENMAX=0.0;vENMAX.clear();                             // for aflowlib_libraries.cpp
  ENMIN=0.0;vENMIN.clear();                             // for aflowlib_libraries.cpp
  POMASS_sum=0.0;POMASS_min=0.0;POMASS_max=0.0;vPOMASS.clear();  // for aflowlib_libraries.cpp
  ZVAL_sum=0.0;ZVAL_min=0.0;ZVAL_max=0.0;vZVAL.clear(); // for aflowlib_libraries.cpp
  EATOM_min=0.0;EATOM_max=0.0;vEATOM.clear();           // for aflowlib_libraries.cpp
  RCORE_min=0.0;RCORE_max=0.0;vRCORE.clear();           // for aflowlib_libraries.cpp
  RWIGS_min=0.0;RWIGS_max=0.0;vRWIGS.clear();           // for aflowlib_libraries.cpp
  EAUG_min=0.0;EAUG_max=0.0;vEAUG.clear();              // for aflowlib_libraries.cpp
  RAUG_min=0.0;RAUG_max=0.0;vRAUG.clear();              // for aflowlib_libraries.cpp
  RMAX_min=0.0;RMAX_max=0.0;vRMAX.clear();              // unicity
  vTITEL.clear();                                       // unicity
  vLEXCH.clear();                                       // unicity 
  // end shared xPOTCAR
  pp_type="";                   // for aflowlib_libraries.cpp
  species.clear();              // for aflowlib_libraries.cpp
  species_Z.clear();            // for aflowlib_libraries.cpp
  species_pp.clear();           // for aflowlib_libraries.cpp
  species_pp_type.clear();      // for aflowlib_libraries.cpp
  species_pp_version.clear();   // for aflowlib_libraries.cpp
  species_pp_AUID.clear();      // for aflowlib_libraries.cpp
  species_pp_AUID_collisions.clear();   // for aflowlib_libraries.cpp
  species_pp_groundstate_energy.clear();      // for aflowlib_libraries.cpp
  species_pp_groundstate_structure.clear();      // for aflowlib_libraries.cpp
  // objects/functions for references energies
  // [OBSOLETE] vsymbol.clear();                              // for pseudopotential references
  // [OBSOLETE] vname.clear();                                // for pseudopotential references
  // [OBSOLETE] vdate.clear();                                // for pseudopotential references
  AUID="";                                      // for pseudopotential references

  //[CO20200404 - DUPLICATE]vcontent.clear(); 
  //[CO20200404 - DUPLICATE]species.clear();                    // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]species_Z.clear();                  // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]species_pp.clear();                 // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]species_pp_type.clear();            // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]species_pp_version.clear();         // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]species_pp_AUID.clear();   // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]species_pp_AUID_collisions.clear();   // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]species_pp_groundstate_energy.clear();   // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]species_pp_groundstate_structure.clear();      // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]// begin shared xPOTCAR
  //[CO20200404 - DUPLICATE]vENMAX.clear();                      // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vENMIN.clear();                      // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vPOMASS.clear();                     // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vZVAL.clear();                       // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vEATOM.clear();                      // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vRCORE.clear();                      // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vRWIGS.clear();                      // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vEAUG.clear();                       // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vRAUG.clear();                       // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vRMAX.clear();                       // for aflowlib_libraries.cpp
  //[CO20200404 - DUPLICATE]vTITEL.clear();                      // unicity
  //[CO20200404 - DUPLICATE]vLEXCH.clear();                      // unicity
  //[CO20200404 - DUPLICATE]// [OBSOLETE] vsymbol.clear();                     // for pseudopotential references
  //[CO20200404 - DUPLICATE]// [OBSOLETE] vname.clear();                       // for pseudopotential references
  //[CO20200404 - DUPLICATE]// [OBSOLETE] vdate.clear();                       // for pseudopotential references
  //[CO20200404 - DUPLICATE]// end shared xPOTCAR
}

void xPOTCAR::copy(const xPOTCAR& b) { // copy PRIVATE
  xStream::copy(b);  //CO20200404 - xStream integration for logging
  free();
  m_initialized=b.m_initialized; //CO20200404 - xStream integration for logging
  content=b.content;
  vcontent.clear(); for(uint i=0;i<b.vcontent.size();i++) vcontent.push_back(b.vcontent.at(i));  // for aflowlib_libraries.cpp
  filename=b.filename;
  title=b.title;
  POTCAR_PAW=b.POTCAR_PAW;
  POTCAR_TYPE=b.POTCAR_TYPE;
  POTCAR_KINETIC=b.POTCAR_KINETIC;
  POTCAR_GW=b.POTCAR_GW;
  POTCAR_AE=b.POTCAR_AE;
  // begin shared xPOTCAR
  ENMAX=b.ENMAX;
  vENMAX.clear(); for(uint i=0;i<b.vENMAX.size();i++) vENMAX.push_back(b.vENMAX.at(i)); // for aflowlib_libraries.cpp
  ENMIN=b.ENMIN;
  vENMIN.clear(); for(uint i=0;i<b.vENMIN.size();i++) vENMIN.push_back(b.vENMIN.at(i)); // for aflowlib_libraries.cpp
  POMASS_sum=b.POMASS_sum;POMASS_min=b.POMASS_min;POMASS_max=b.POMASS_max;
  vPOMASS.clear(); for(uint i=0;i<b.vPOMASS.size();i++) vPOMASS.push_back(b.vPOMASS.at(i)); // for aflowlib_libraries.cpp
  ZVAL_sum=b.ZVAL_sum;ZVAL_min=b.ZVAL_min;ZVAL_max=b.ZVAL_max;
  vZVAL.clear(); for(uint i=0;i<b.vZVAL.size();i++) vZVAL.push_back(b.vZVAL.at(i)); // for aflowlib_libraries.cpp
  EATOM_min=b.EATOM_min;EATOM_max=b.EATOM_max;
  vEATOM.clear(); for(uint i=0;i<b.vEATOM.size();i++) vEATOM.push_back(b.vEATOM.at(i)); // for aflowlib_libraries.cpp
  RCORE_min=b.RCORE_min;RCORE_max=b.RCORE_max;
  vRCORE.clear(); for(uint i=0;i<b.vRCORE.size();i++) vRCORE.push_back(b.vRCORE.at(i)); // for aflowlib_libraries.cpp
  RWIGS_min=b.RWIGS_min;RWIGS_max=b.RWIGS_max;
  vRWIGS.clear(); for(uint i=0;i<b.vRWIGS.size();i++) vRWIGS.push_back(b.vRWIGS.at(i)); // for aflowlib_libraries.cpp
  EAUG_min=b.EAUG_min;EAUG_max=b.EAUG_max;
  vEAUG.clear(); for(uint i=0;i<b.vEAUG.size();i++) vEAUG.push_back(b.vEAUG.at(i)); // for aflowlib_libraries.cpp
  RAUG_min=b.RAUG_min;RAUG_max=b.RAUG_max;
  vRAUG.clear(); for(uint i=0;i<b.vRAUG.size();i++) vRAUG.push_back(b.vRAUG.at(i)); // for aflowlib_libraries.cpp
  RMAX_min=b.RMAX_min;RMAX_max=b.RMAX_max;
  vRMAX.clear(); for(uint i=0;i<b.vRMAX.size();i++) vRMAX.push_back(b.vRMAX.at(i)); // for aflowlib_libraries.cpp
  vTITEL.clear(); for(uint i=0;i<b.vTITEL.size();i++) vTITEL.push_back(b.vTITEL.at(i)); // for aflowlib_libraries.cpp
  vLEXCH.clear(); for(uint i=0;i<b.vLEXCH.size();i++) vLEXCH.push_back(b.vLEXCH.at(i)); // for aflowlib_libraries.cpp
  // end shared xPOTCAR
  pp_type=b.pp_type;                   // for aflowlib_libraries.cpp
  species.clear(); for(uint i=0;i<b.species.size();i++) species.push_back(b.species.at(i)); // for aflowlib_libraries.cpp
  species_Z.clear(); for(uint i=0;i<b.species_Z.size();i++) species_Z.push_back(b.species_Z.at(i)); // for aflowlib_libraries.cpp
  species_pp.clear(); for(uint i=0;i<b.species_pp.size();i++) species_pp.push_back(b.species_pp.at(i)); // for aflowlib_libraries.cpp
  species_pp_type.clear(); for(uint i=0;i<b.species_pp_type.size();i++) species_pp_type.push_back(b.species_pp_type.at(i)); // for aflowlib_libraries.cpp
  species_pp_version.clear(); for(uint i=0;i<b.species_pp_version.size();i++) species_pp_version.push_back(b.species_pp_version.at(i));  // for aflowlib_libraries.cpp
  species_pp_AUID.clear(); for(uint i=0;i<b.species_pp_AUID.size();i++) species_pp_AUID.push_back(b.species_pp_AUID.at(i));  // for aflowlib_libraries.cpp
  species_pp_AUID_collisions.clear(); for(uint i=0;i<b.species_pp_AUID_collisions.size();i++) species_pp_AUID_collisions.push_back(b.species_pp_AUID_collisions.at(i));  // for aflowlib_libraries.cpp
  species_pp_groundstate_energy.clear(); for(uint i=0;i<b.species_pp_groundstate_energy.size();i++) species_pp_groundstate_energy.push_back(b.species_pp_groundstate_energy.at(i)); // for pseudopotential references
  species_pp_groundstate_structure.clear(); for(uint i=0;i<b.species_pp_groundstate_structure.size();i++) species_pp_groundstate_structure.push_back(b.species_pp_groundstate_structure.at(i)); // for pseudopotential references
  // objects/functions for references energies
  // [OBSOLETE] vsymbol.clear(); for(uint i=0;i<b.vsymbol.size();i++) vsymbol.push_back(b.vsymbol.at(i)); // for pseudopotential references
  // [OBSOLETE] vname.clear(); for(uint i=0;i<b.vname.size();i++) vname.push_back(b.vname.at(i)); // for pseudopotential references
  // [OBSOLETE] vdate.clear(); for(uint i=0;i<b.vdate.size();i++) vdate.push_back(b.vdate.at(i)); // for pseudopotential references
  AUID=b.AUID; // for pseudopotential references
}

const xPOTCAR& xPOTCAR::operator=(const xPOTCAR& b) {  // operator= PUBLIC
  if(this!=&b) {free();copy(b);}
  return *this;
}

void xPOTCAR::clear() {  // clear PRIVATE
  xPOTCAR _temp;
  string filename_aus=filename;
  copy(_temp);
  filename=filename_aus;
}

bool xPOTCAR::GetProperties(const string& stringIN,bool QUIET) {
  stringstream sss; sss.str(stringIN);
  if(filename=="") filename="string";
  return xPOTCAR::GetProperties(sss,QUIET);
}

bool xPOTCAR::GetPropertiesFile(const string& fileIN,bool QUIET) {
  stringstream sss;
  if(filename=="") filename=fileIN;
  aurostd::efile2stringstream(fileIN,sss);
  return xPOTCAR::GetProperties(sss,QUIET);
}

bool xPOTCAR::GetPropertiesUrlFile(const string& url,const string& file,bool QUIET) {
  //[CO20200502 - OBSOLETE]string tmpfile=XHOST.tmpfs+"/_aflow_"+XHOST.user+".pid"+XHOST.ostrPID.str()+".a"+string(AFLOW_VERSION)+".rnd"+aurostd::utype2string(uint((double) std::floor((double)100000*aurostd::ran0())))+".u"+aurostd::utype2string(uint((double) aurostd::get_useconds()))+"_"+file;
  string tmpfile=aurostd::TmpFileCreate("xPOTCAR_GetProperties"); //CO20200502 - threadID
  aurostd::url2file(url+"/"+file,tmpfile,!QUIET);
  bool out=GetPropertiesFile(tmpfile,QUIET);
  filename="url="+url;  //CO20210315
  aurostd::RemoveFile(tmpfile);
  return out;
}

#define PSEUDOPOTENTIAL_GENERATOR_pad 70

bool xPOTCAR::GetProperties(const stringstream& stringstreamIN,bool QUIET) {
  bool LDEBUG=(FALSE || XHOST.DEBUG || !QUIET);
  stringstream message;
  bool force_exit=XHOST.POSTPROCESS; //SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only

  bool ERROR_flag=FALSE;
  // XHOST.PSEUDOPOTENTIAL_GENERATOR=TRUE;
  long double seconds=aurostd::get_seconds();
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN (" << time_delay(seconds) << ")" << endl;
  clear(); // so it does not mess up vector/deque
  content=stringstreamIN.str();
  if(1&&LDEBUG){cerr << __AFLOW_FUNC__ << endl;cerr << content;}
  vcontent.clear();
  vector<string> vline,tokens;
  aurostd::string2vectorstring(content,vcontent);
  string line;
  if(filename=="") filename="stringstream";
  // crunching to eat the info

  if(vcontent.size()==0){return false;} //CO+ME20200604 - file does not exist, don't spit out warnings for the workshops

  // GET AUID AS SOON AS POSSIBLE
  if(aurostd::FileExist(filename)) {
    AUID=aurostd::file2auid(filename);
  }

  // get parameters
  //  vline.clear();
  for(uint iline=0;iline<vcontent.size();iline++) {
    aurostd::string2tokens(vcontent.at(iline),tokens);
    //    cerr << "iline=" << iline << "  " << vcontent.at(iline) << " tokens.size()=" << tokens.size() << endl;
    if(iline==4) title=vcontent.at(iline);  // might be wrong
  } 

  // ----------------------------------------------------------------------
  // EATOM RCORE RWIGS EAUG RAUG ENMAX ENMIN POMASS ZVAL RMAX
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD \"EATOM RCORE RWIGS EAUG RAUG ENMAX ENMIN POMASS ZVAL RMAX\" DATA" << endl;
  vline.clear();
  vENMAX.clear();vENMIN.clear();vPOMASS.clear();vZVAL.clear();vEATOM.clear();vRCORE.clear();vRWIGS.clear();vEAUG.clear();vRAUG.clear();vRMAX.clear();
  for(uint iline=0;iline<vcontent.size();iline++) {
    aurostd::StringSubst(vcontent.at(iline),"EMMIN","ENMIN"); // Kresseeeeeeeeeee
    if(aurostd::substring2bool(vcontent.at(iline),"ENMAX") && aurostd::substring2bool(vcontent.at(iline),"ENMIN")) vline.push_back(vcontent.at(iline));
    if(aurostd::substring2bool(vcontent.at(iline),"POMASS") && aurostd::substring2bool(vcontent.at(iline),"ZVAL")) vline.push_back(vcontent.at(iline));
    if(aurostd::substring2bool(vcontent.at(iline),"EATOM") && aurostd::substring2bool(vcontent.at(iline),"eV")) vline.push_back(vcontent.at(iline));
    if(aurostd::substring2bool(vcontent.at(iline),"RCORE") && aurostd::substring2bool(vcontent.at(iline),"radius")) vline.push_back(vcontent.at(iline));
    if(aurostd::substring2bool(vcontent.at(iline),"RWIGS") && aurostd::substring2bool(vcontent.at(iline),"radius")) vline.push_back(vcontent.at(iline));
    if(aurostd::substring2bool(vcontent.at(iline),"EAUG")) vline.push_back(vcontent.at(iline));
    if(aurostd::substring2bool(vcontent.at(iline),"RAUG") && aurostd::substring2bool(vcontent.at(iline),"sphere")) vline.push_back(vcontent.at(iline));
    if(aurostd::substring2bool(vcontent.at(iline),"RMAX") && aurostd::substring2bool(vcontent.at(iline),"radius")) vline.push_back(vcontent.at(iline));
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vline.size()=" << vline.size() << endl;
  if(vline.size()==0) {
    //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - " << __AFLOW_FUNC__ << " wrong number of \"EATOM RCORE RWIGS EAUG RAUG ENMAX ENMIN POMASS ZVAL RMAX\" in POTCAR" << "   filename=[" << filename << "]" << endl;ERROR_flag=TRUE;
    message << "Wrong number of \"EATOM RCORE RWIGS EAUG RAUG ENMAX ENMIN POMASS ZVAL RMAX\" in POTCAR" << "   filename=[" << filename << "]" << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  }
  for(uint j=0;j<vline.size();j++) {
    aurostd::StringSubst(vline.at(j),"="," ");
    aurostd::StringSubst(vline.at(j),";"," ");
    //    if(LDEBUG) cerr << vline.at(j) << endl;
    aurostd::string2tokens(vline.at(j),tokens," ");
    for(uint k=0;k<tokens.size();k++) {
      if(tokens.at(k)=="ENMAX" && k+1<tokens.size()) vENMAX.push_back(aurostd::string2utype<double>(tokens.at(k+1)));
      if(tokens.at(k)=="ENMIN" && k+1<tokens.size()) vENMIN.push_back(aurostd::string2utype<double>(tokens.at(k+1)));
      if(tokens.at(k)=="POMASS" && k+1<tokens.size()) vPOMASS.push_back(aurostd::string2utype<double>(tokens.at(k+1)));
      if(tokens.at(k)=="ZVAL" && k+1<tokens.size()) vZVAL.push_back(aurostd::string2utype<double>(tokens.at(k+1)));
      if(tokens.at(k)=="EATOM" && k+1<tokens.size()) vEATOM.push_back(aurostd::string2utype<double>(tokens.at(k+1)));
      if(tokens.at(k)=="RCORE" && k+1<tokens.size()) vRCORE.push_back(aurostd::string2utype<double>(tokens.at(k+1)));
      if(tokens.at(k)=="RWIGS" && k==0 && k+1<tokens.size()) vRWIGS.push_back(aurostd::string2utype<double>(tokens.at(k+1))); // pick the 1st
      if(tokens.at(k)=="EAUG" && k+1<tokens.size()) vEAUG.push_back(aurostd::string2utype<double>(tokens.at(k+1)));
      if(tokens.at(k)=="RAUG" && k+1<tokens.size()) vRAUG.push_back(aurostd::string2utype<double>(tokens.at(k+1)));
      if(tokens.at(k)=="RMAX" && k+1<tokens.size()) vRMAX.push_back(aurostd::string2utype<double>(tokens.at(k+1)));
    }
  }

  if(vENMIN.size()<vENMAX.size()) vENMIN.push_back(0.0);  // ENMIN MIGHT NOT BE THERE
  if(vRCORE.size()<vENMAX.size()) vRCORE.push_back(0.0);  // RCORE MIGHT NOT BE THERE
  if(vRWIGS.size()<vENMAX.size()) vRWIGS.push_back(0.0);  // RWIGS MIGHT NOT BE THERE
  if(vEAUG.size()<vENMAX.size()) vEAUG.push_back(0.0);    // EAUG MIGHT NOT BE THERE
  if(vRAUG.size()<vENMAX.size()) vRAUG.push_back(0.0);    // RAUG MIGHT NOT BE THERE
  if(vRMAX.size()<vENMAX.size()) vRMAX.push_back(0.0);    // RMAX MIGHT NOT BE THERE

  ENMAX=aurostd::max(vENMAX);
  ENMIN=aurostd::min(vENMIN);
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " ENMAX=" << ENMAX << " vENMAX.size()=" << vENMAX.size() << ": "; for(uint i=0;i<vENMAX.size();i++) { cerr << vENMAX.at(i) << " "; } cerr << " " << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " ENMIN=" << ENMIN << " vENMIN.size()=" << vENMIN.size() << ": "; for(uint i=0;i<vENMIN.size();i++) { cerr << vENMIN.at(i) << " "; } cerr << " " << endl;}

  POMASS_sum=aurostd::sum(vPOMASS);POMASS_min=aurostd::min(vPOMASS);POMASS_max=aurostd::max(vPOMASS);
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " POMASS_sum=" << POMASS_sum << " POMASS_min=" << POMASS_min << " POMASS_max=" << POMASS_max << " vPOMASS.size()=" << vPOMASS.size() << ": "; for(uint i=0;i<vPOMASS.size();i++) { cerr << vPOMASS.at(i) << " "; } cerr << " " << endl;}

  ZVAL_sum=aurostd::sum(vZVAL);ZVAL_min=aurostd::min(vZVAL);ZVAL_max=aurostd::max(vZVAL);
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " ZVAL_sum=" << ZVAL_sum << " ZVAL_min=" << ZVAL_min << " ZVAL_max=" << ZVAL_max << " vZVAL.size()=" << vZVAL.size() << ": "; for(uint i=0;i<vZVAL.size();i++) { cerr << vZVAL.at(i) << " "; } cerr << " " << endl;}

  EATOM_min=aurostd::min(vEATOM);EATOM_max=aurostd::max(vEATOM);
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " EATOM_min=" << EATOM_min << " EATOM_max=" << EATOM_max << " vEATOM.size()=" << vEATOM.size() << ": "; for(uint i=0;i<vEATOM.size();i++) { cerr << vEATOM.at(i) << " "; } cerr << " " << endl;}

  RCORE_min=aurostd::min(vRCORE);RCORE_max=aurostd::max(vRCORE);
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " RCORE_min=" << RCORE_min << " RCORE_max=" << RCORE_max << " vRCORE.size()=" << vRCORE.size() << ": "; for(uint i=0;i<vRCORE.size();i++) { cerr << vRCORE.at(i) << " "; } cerr << " " << endl;}

  RWIGS_min=aurostd::min(vRWIGS);RWIGS_max=aurostd::max(vRWIGS);
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " RWIGS_min=" << RWIGS_min << " RWIGS_max=" << RWIGS_max << " vRWIGS.size()=" << vRWIGS.size() << ": "; for(uint i=0;i<vRWIGS.size();i++) { cerr << vRWIGS.at(i) << " "; } cerr << " " << endl;}

  EAUG_min=aurostd::min(vEAUG);EAUG_max=aurostd::max(vEAUG);
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " EAUG_min=" << EAUG_min << " EAUG_max=" << EAUG_max << " vEAUG.size()=" << vEAUG.size() << ": "; for(uint i=0;i<vEAUG.size();i++) { cerr << vEAUG.at(i) << " "; } cerr << " " << endl;}

  RAUG_min=aurostd::min(vRAUG);RAUG_max=aurostd::max(vRAUG);
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " RAUG_min=" << RAUG_min << " RAUG_max=" << RAUG_max << " vRAUG.size()=" << vRAUG.size() << ": "; for(uint i=0;i<vRAUG.size();i++) { cerr << vRAUG.at(i) << " "; } cerr << " " << endl;}

  RMAX_min=aurostd::min(vRMAX);RMAX_max=aurostd::max(vRMAX);
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " RMAX_min=" << RMAX_min << " RMAX_max=" << RMAX_max << " vRMAX.size()=" << vRMAX.size() << ": "; for(uint i=0;i<vRMAX.size();i++) { cerr << vRMAX.at(i) << " "; } cerr << " " << endl;}


  // ----------------------------------------------------------------------
  // PSEUDOPOTENTIAL DATA
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD PSEUDOPOTENTIAL DATA" << endl;

  vline.clear();
  vTITEL.clear();
  vLEXCH.clear();
  for(uint iline=0;iline<vcontent.size();iline++)  {
    if(aurostd::substring2bool(vcontent.at(iline),"=")) {
      if(aurostd::substring2bool(vcontent.at(iline),"TITEL")) {
        aurostd::string2tokens(vcontent.at(iline),tokens,"=");
        if(tokens.size()>1) vTITEL.push_back(tokens.at(1));
        //	if(LDEBUG) cerr << __AFLOW_FUNC__ << " TITEL=" << tokens.at(1) << endl;
      }
      if(aurostd::substring2bool(vcontent.at(iline),"LEXCH")) {
        if(!aurostd::substring2bool(vcontent.at(iline),"exchange")) {
          aurostd::string2tokens(vcontent.at(iline),tokens,"=");
          if(tokens.size()>1) vLEXCH.push_back(tokens.at(1));
          //	if(LDEBUG) cerr << __AFLOW_FUNC__ << " LEXCH=" << tokens.at(1) << endl;
        }
      }
    }
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vTITEL.size()=" << vTITEL.size() << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vLEXCH.size()=" << vLEXCH.size() << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vRMAX.size()=" << vRMAX.size() << endl;

  //ME20200604 - Only exit if force_exit is true
  if(vTITEL.size()==0) {
    //[CO20200404 - OBSOLETE]if(!QUIET) cerr << "WARNING - " << __AFLOW_FUNC__ << " wrong number of pseudopotentials (TITEL) in POTCAR" << "   filename=[" << filename << "]" << endl;ERROR_flag=TRUE;
    message << "Wrong number of pseudopotentials (TITEL) in POTCAR" << "   filename=[" << filename << "]" << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  } //CO20200106 - patching for auto-indenting
  if(vLEXCH.size()==0) {
    //[CO20200604 - OBSOLETE]if(!QUIET) cerr << "WARNING - " << __AFLOW_FUNC__ << " wrong number of pseudopotentials (LEXCH) in POTCAR" << "   filename=[" << filename << "]" << endl;ERROR_flag=TRUE;
    message << "Wrong number of pseudopotentials (LEXCH) in POTCAR" << "   filename=[" << filename << "]" << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  } //CO20200106 - patching for auto-indenting
  if(vEATOM.size()==0) {
    //[CO20200604 - OBSOLETE]if(!QUIET) cerr << "WARNING - " << __AFLOW_FUNC__ << " wrong number of pseudopotentials (EATOM) in POTCAR" << "   filename=[" << filename << "]" << endl;ERROR_flag=TRUE;
    message << "Wrong number of pseudopotentials (EATOM) in POTCAR" << "   filename=[" << filename << "]" << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  } //CO20200106 - patching for auto-indenting
  if(vRMAX.size()==0) {
    //[CO20200604 - OBSOLETE]if(!QUIET) cerr << "WARNING - " << __AFLOW_FUNC__ << " wrong number of pseudopotentials (RMAX) in POTCAR" << "   filename=[" << filename << "]" << endl;ERROR_flag=TRUE;
    message << "Wrong number of pseudopotentials (RMAX) in POTCAR" << "   filename=[" << filename << "]" << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  }   //CO20200106 - patching for auto-indenting
  if(vTITEL.size()!=vLEXCH.size()) {
    //[CO20200604 - OBSOLETE]if(!QUIET) cerr << "WARNING - " << __AFLOW_FUNC__ << " wrong number of pseudopotentials (TITEL/LEXCH) in POTCAR" << "   filename=[" << filename << "]" << endl;ERROR_flag=TRUE;
    message << "Wrong number of pseudopotentials (TITEL/LEXCH) in POTCAR" << "   filename=[" << filename << "]" << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  } //CO20200106 - patching for auto-indenting
  if(vLEXCH.size()!=vEATOM.size()) {
    //[CO20200604 - OBSOLETE]if(!QUIET) cerr << "WARNING - " << __AFLOW_FUNC__ << " wrong number of pseudopotentials (LEXCH/EATOM) in POTCAR" << "   filename=[" << filename << "]" << endl;ERROR_flag=TRUE;
    message << "Wrong number of pseudopotentials (LEXCH/EATOM) in POTCAR" << "   filename=[" << filename << "]" << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  } //CO20200106 - patching for auto-indenting
  if(vEATOM.size()!=vRMAX.size()) {
    //[CO20200604 - OBSOLETE]if(!QUIET) cerr << "WARNING - " << __AFLOW_FUNC__ << " wrong number of pseudopotentials (EATOM/RMAX) in POTCAR" << "   filename=[" << filename << "]" << endl;ERROR_flag=TRUE;
    message << "Wrong number of pseudopotentials (EATOM/RMAX) in POTCAR" << "   filename=[" << filename << "]" << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
    ERROR_flag=TRUE;
  }   //CO20200106 - patching for auto-indenting

  for(uint j=0;j<vTITEL.size();j++) {
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " SPECIES(" << j << ") " << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " vTITEL.at(j)=" << vTITEL.at(j) << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " vLEXCH.at(j)=" << vLEXCH.at(j) << endl;
    aurostd::string2tokens(vTITEL.at(j),tokens," ");
    pp_type=tokens.at(0);
    if(pp_type=="US" or pp_type=="NC") {
      pp_type="GGA";
      for(uint iiline=0;iiline<vcontent.size();iiline++) 
        if(aurostd::substring2bool(vcontent.at(iiline),"GGA"))
          if(aurostd::substring2bool(vcontent.at(iiline),"eV"))
            pp_type="LDA";
    }
    if(tokens.size()>1) {
      if(pp_type=="PAW") pp_type="PAW_LDA"; // cerr << __AFLOW_FUNC__ << " PAW_LDA" << endl;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " pp_type=" << pp_type << " " << "tokens.size()=" << tokens.size() << endl;
      species.push_back(KBIN::VASP_PseudoPotential_CleanName(tokens.at(1)));
      species_Z.push_back(xelement::symbol2Z(KBIN::VASP_PseudoPotential_CleanName(tokens.at(1))));
      species_pp.push_back(tokens.at(1));
      species_pp_type.push_back(pp_type);
      if(pp_type=="LDA" && tokens.size()<3) tokens.push_back(DEFAULT_VASP_POTCAR_DATE_POT_LDA);
      if(pp_type=="GGA" && tokens.size()<3) tokens.push_back(DEFAULT_VASP_POTCAR_DATE_POT_GGA);
      species_pp_version.push_back(tokens.at(1)+":"+pp_type+":"+tokens.at(2));
    } else { // only  one definition
      pp_type="AE";
      species.push_back(KBIN::VASP_PseudoPotential_CleanName(tokens.at(0))); 
      species_Z.push_back(xelement::symbol2Z(KBIN::VASP_PseudoPotential_CleanName(tokens.at(0)))); 
      species_pp.push_back(KBIN::VASP_PseudoPotential_CleanName(tokens.at(0))); // same as species
      species_pp_type.push_back(pp_type);
      species_pp_version.push_back(tokens.at(0)+":"+pp_type+":"+"BigBang"); // ALL Atomic Potential has been existing since the BigBang
    }

    vTITEL.at(j)=aurostd::RemoveWhiteSpaces(vTITEL.at(j));
    vLEXCH.at(j)=aurostd::RemoveWhiteSpaces(vLEXCH.at(j));   
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " vTITEL.at(j)=" << vTITEL.at(j) << endl;
      cerr << __AFLOW_FUNC__ << " vLEXCH.at(j)=" << vLEXCH.at(j) << endl;
      cerr << __AFLOW_FUNC__ << " vEATOM.at(j)=" << vEATOM.at(j) << endl;
      cerr << __AFLOW_FUNC__ << " vRMAX.at(j)=" << vRMAX.at(j) << endl;
    }
    xPOTCAR xPOT(xPOTCAR_Finder(species_pp_AUID,species_pp_AUID_collisions,vTITEL.at(j),vLEXCH.at(j),vEATOM.at(j),vRMAX.at(j),LDEBUG)); // FIXES species_pp_AUID,species_pp_AUID_collisions
    species_pp_groundstate_energy.push_back(xPOT.species_pp_groundstate_energy.at(0));
    species_pp_groundstate_structure.push_back(xPOT.species_pp_groundstate_structure.at(0));

    if(LDEBUG) cerr << __AFLOW_FUNC__ << " SPECIES(" << j << ") [pp, type, version, AUID] = "
      << species.at(species.size()-1) << " ["
        << species_Z.at(species_Z.size()-1) << ", "
        << species_pp.at(species_pp.size()-1) << ", "
        << species_pp_type.at(species_pp_type.size()-1) << ", "
        << species_pp_version.at(species_pp_version.size()-1) << ", "
        << species_pp_AUID.at(species_pp_AUID.size()-1) <<  ", "
        << species_pp_groundstate_energy.at(species_pp_groundstate_energy.size()-1) <<  ", "
        << species_pp_groundstate_structure.at(species_pp_groundstate_structure.size()-1) << "]" << endl;
  }

  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " PSEUDOPOTENTIAL type = " << pp_type << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " species.size()=" << species.size() << ": "; for(uint i=0;i<species.size();i++) { cerr << species.at(i) << " "; } cerr << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " species_Z.size()=" << species_Z.size() << ": "; for(uint i=0;i<species_Z.size();i++) { cerr << species_Z.at(i) << " "; } cerr << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " species_pp.size()=" << species_pp.size() << ": "; for(uint i=0;i<species_pp.size();i++) { cerr << species_pp.at(i) << " "; } cerr << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " species_pp_type.size()=" << species_pp_type.size() << ": "; for(uint i=0;i<species_pp_type.size();i++) { cerr << species_pp_type.at(i) << " "; } cerr << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " species_pp_version.size()=" << species_pp_version.size() << ": "; for(uint i=0;i<species_pp_version.size();i++) { cerr << species_pp_version.at(i) << " "; } cerr << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " species_pp_AUID.size()=" << species_pp_AUID.size() << ": "; for(uint i=0;i<species_pp_AUID.size();i++) { cerr << species_pp_AUID.at(i) << " "; } cerr << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " species_pp_AUID_collisions.size()=" << species_pp_AUID_collisions.size() << ": "; for(uint i=0;i<species_pp_AUID_collisions.size();i++) { cerr << species_pp_AUID_collisions.at(i) << " "; } cerr << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " species_pp_groundstate_energy.size()=" << species_pp_groundstate_energy.size() << ": "; for(uint i=0;i<species_pp_groundstate_energy.size();i++) { cerr << species_pp_groundstate_energy.at(i) << " "; } cerr << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " species_pp_groundstate_structure.size()=" << species_pp_groundstate_structure.size() << ": "; for(uint i=0;i<species_pp_groundstate_structure.size();i++) { cerr << species_pp_groundstate_structure.at(i) << " "; } cerr << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " vTITEL.size()=" << vTITEL.size() << ": "; for(uint i=0;i<vTITEL.size();i++) { cerr << vTITEL.at(i) << " "; } cerr << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " vLEXCH.size()=" << vLEXCH.size() << ": "; for(uint i=0;i<vLEXCH.size();i++) { cerr << vLEXCH.at(i) << " "; } cerr << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " filename=" << filename << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " AUID=" << AUID << endl;}
  if(species_pp_AUID_collisions.size()) {
    cerr << __AFLOW_FUNC__ << " COLLISION species_pp_AUID_collisions.size()=" << species_pp_AUID_collisions.size() << endl;
    ERROR_flag=TRUE;
  }
  if(LDEBUG) cerr.flush();

  // ----------------------------------------------------------------------
  // POTCAR_PAW POTCAR_TYPE POTCAR_KINETIC
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD POTCAR_PAW POTCAR_TYPE" << endl;
  POTCAR_PAW=FALSE;
  POTCAR_KINETIC=FALSE;
  POTCAR_GW=FALSE;
  POTCAR_AE=FALSE;
  if(content.find("PAW")!=string::npos) POTCAR_PAW=TRUE;
  if(content.find("GW")!=string::npos) POTCAR_GW=TRUE;
  if(content.find("AE potential")!=string::npos) POTCAR_AE=TRUE;
  // POTCAR DONE **************************************************
  // CHECK FOR US => LDA/GGA
  bool is_US=FALSE,is_LDA=FALSE;
  POTCAR_TYPE="";
  for(uint i=0;i<vcontent.size()&&POTCAR_TYPE.empty();i++) { // cerr << vcontent.at(i) << endl;
    if(aurostd::substring2bool(vcontent.at(i),"TITEL") && aurostd::substring2bool(vcontent.at(i),"US")) is_US=TRUE;
    if(aurostd::substring2bool(vcontent.at(i),"GGA ")) is_LDA=TRUE;
    if(aurostd::substring2bool(vcontent.at(i),"mkinetic energy-density pseudized")) POTCAR_KINETIC=TRUE;
  }
  if(is_US && is_LDA) POTCAR_TYPE="LDA";
  if(is_US && !is_LDA) POTCAR_TYPE="GGA";
  for(uint i=0;i<vcontent.size()&&POTCAR_TYPE.empty();i++) {
    if(aurostd::substring2bool(vcontent.at(i),"TITEL")) {
      aurostd::string2tokens(vcontent.at(i),tokens,"=");vcontent.at(i)=tokens.at(1);aurostd::string2tokens(vcontent.at(i),tokens," ");
      if(tokens.at(0)=="PAW") POTCAR_TYPE="PAW_LDA";
      if(tokens.at(0)=="PAW_GGA") POTCAR_TYPE="PAW_GGA";
      if(tokens.at(0)=="PAW_PBE") POTCAR_TYPE="PAW_PBE";
      if(tokens.at(0)=="PAW_RPBE") POTCAR_TYPE="PAW_RPBE";
      if(tokens.at(0)=="PAW_PBE" && POTCAR_KINETIC) POTCAR_TYPE="PAW_PBE_KIN";
      if(tokens.at(0)=="PAW" && POTCAR_KINETIC) POTCAR_TYPE="PAW_LDA_KIN";
    }
  }    
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " POTCAR_PAW = " << POTCAR_PAW << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " POTCAR_TYPE = " << POTCAR_TYPE << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " POTCAR_KINETIC = " << POTCAR_KINETIC << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " POTCAR_GW = " << POTCAR_GW << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " POTCAR_AE = " << POTCAR_AE << endl;

  // ----------------------------------------------------------------------
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " WRITE vxpseudopotentials" << endl;

  if(XHOST.PSEUDOPOTENTIAL_GENERATOR && species.size()==1) xPOTCAR_PURE_Printer((*this),cout);

  // ----------------------------------------------------------------------
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " (BULLSHIT DONT USE) title=" << title << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vcontent.size()=" << vcontent.size() << endl;

  // ----------------------------------------------------------------------
  // DONE NOW RETURN  
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " END (" << time_delay(seconds) << ")" << endl;
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  //[CO20200404 - OBSOLETE]if(ERROR_flag && !QUIET) cerr << "WARNING - " << __AFLOW_FUNC__ << " ERROR_flag set in xPOTCAR" << endl;
  if(ERROR_flag){
    message << "ERROR_flag set in xPOTCAR";
    if(force_exit){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);}
    else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_,QUIET);}
    return FALSE;
  }
  return TRUE;
}

//---------------------------------------------------------------------------------
// class xVASPRUNXML
//---------------------------------------------------------------------------------
xVASPRUNXML::xVASPRUNXML(ostream& oss) : xStream(oss) {free();} //CO20200404 - xStream integration for logging
xVASPRUNXML::xVASPRUNXML(ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss) {free();} //CO20200404 - xStream integration for logging
xVASPRUNXML::xVASPRUNXML(const string& fileIN,bool QUIET,ostream& oss) : xStream(oss) {m_initialized=initialize(fileIN,QUIET);}  //CO20200404 - xStream integration for logging
xVASPRUNXML::xVASPRUNXML(const string& fileIN,ofstream& FileMESSAGE,bool QUIET,ostream& oss) : xStream(FileMESSAGE,oss) {m_initialized=initialize(fileIN,QUIET);}  //CO20200404 - xStream integration for logging

bool xVASPRUNXML::initialize(const string& fileIN,ostream& oss,bool QUIET) {xStream::initialize(oss);return initialize(fileIN,QUIET);} //CO20200508
bool xVASPRUNXML::initialize(const string& fileIN,ofstream& FileMESSAGE,ostream& oss,bool QUIET) {xStream::initialize(FileMESSAGE,oss);return initialize(fileIN,QUIET);} //CO20200508
bool xVASPRUNXML::initialize(const string& fileIN,bool QUIET) {
  free();
  filename = fileIN;
  return GetPropertiesFile(fileIN,QUIET);
}

xVASPRUNXML::xVASPRUNXML(const xVASPRUNXML& b) : xStream(*b.getOFStream(),*b.getOSS()) {free();copy(b);} // copy PUBLIC

xVASPRUNXML::~xVASPRUNXML() {xStream::free();free();}  //CO20200404 - xStream integration for logging

void xVASPRUNXML::free() {
  m_initialized=false; //CO20200404
  //------------------------------------------------------------------------------
  // GetProperties
  content="";                   // for aflowlib_libraries.cpp
  vcontent.clear();             // for aflowlib_libraries.cpp
  filename="";                  // for aflowlib_libraries.cpp
  natoms=0.0;                   // for aflowlib_libraries.cpp
  stress.clear();               // for aflowlib_libraries.cpp
  vkpoint.clear();          // for aflowlib_libraries.cpp
  vweights.clear();             // for aflowlib_libraries.cpp
  vforces.clear();              // for aflowlib_libraries.cpp
}        

void xVASPRUNXML::copy(const xVASPRUNXML& b) { // copy PRIVATE
  content=b.content;
  vcontent.clear(); for(uint i=0;i<b.vcontent.size();i++) vcontent.push_back(b.vcontent.at(i));  // for aflowlib_libraries.cpp
  filename=b.filename;
  natoms=b.natoms;                   // for aflowlib_libraries.cpp
  stress=b.stress;                   // for aflowlib_libraries.cpp
  vkpoint.clear(); for(uint i=0;i<b.vkpoint.size();i++) vkpoint.push_back(b.vkpoint.at(i));  // for aflowlib_libraries.cpp
  vweights.clear(); for(uint i=0;i<b.vweights.size();i++) vweights.push_back(b.vweights.at(i));  // for aflowlib_libraries.cpp
  vforces.clear(); for(uint i=0;i<b.vforces.size();i++) vforces.push_back(b.vforces.at(i));  // for aflowlib_libraries.cpp
}

const xVASPRUNXML& xVASPRUNXML::operator=(const xVASPRUNXML& b) {  // operator= PUBLIC
  if(this!=&b) {free();copy(b);}
  return *this;
}

void xVASPRUNXML::clear() {  // clear PRIVATE
  xVASPRUNXML _temp;
  string filename_aus=filename;
  copy(_temp);
  filename=filename_aus;
}

bool xVASPRUNXML::GetProperties(const string& stringIN,bool QUIET) {
  stringstream sss; sss.str(stringIN);
  if(filename=="") filename="string";
  return xVASPRUNXML::GetProperties(sss,QUIET);
}

bool xVASPRUNXML::GetPropertiesFile(const string& fileIN,bool QUIET) {
  stringstream sss;
  if(filename=="") filename=fileIN;
  aurostd::efile2stringstream(fileIN,sss);
  return xVASPRUNXML::GetProperties(sss,QUIET);
}

bool xVASPRUNXML::GetPropertiesUrlFile(const string& url,const string& file,bool QUIET) {
  //[CO20200502 - OBSOLETE]string tmpfile=XHOST.tmpfs+"/_aflow_"+XHOST.user+".pid"+XHOST.ostrPID.str()+".a"+string(AFLOW_VERSION)+".rnd"+aurostd::utype2string(uint((double) std::floor((double)100000*aurostd::ran0())))+".u"+aurostd::utype2string(uint((double) aurostd::get_useconds()))+"_"+file;
  string tmpfile=aurostd::TmpFileCreate("xVASPRUNXML_GetProperties"); //CO20200502 - threadID
  aurostd::url2file(url+"/"+file,tmpfile,!QUIET);
  bool out=GetPropertiesFile(tmpfile,QUIET);
  filename="url="+url;  //CO20210315
  aurostd::RemoveFile(tmpfile);
  return out;
}

bool xVASPRUNXML::GetProperties(const stringstream& stringstreamIN,bool QUIET) {
  bool LDEBUG=(FALSE || XHOST.DEBUG || !QUIET);
  stringstream message;
  bool force_exit=XHOST.POSTPROCESS; //SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only
  bool ERROR_flag=FALSE;
  long double seconds=aurostd::get_seconds();
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN (" << time_delay(seconds) << ")" << endl;
  clear(); // so it does not mess up vector/deque
  content=stringstreamIN.str();
  vcontent.clear();
  vector<string> tokens;
  aurostd::string2vectorstring(content,vcontent);
  string line;
  if(filename=="") filename="stringstream";
  // crunching to eat the info

  // ----------------------------------------------------------------------
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vcontent.size()=" << vcontent.size() << endl;

  // ----------------------------------------------------------------------
  // LOAD natoms

  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD natoms DATA" << endl;
  line="";
  for(int iline=(int)vcontent.size()-1;iline>=0;iline--)  // NEW FROM THE BACK
    if(aurostd::substring2bool(vcontent.at(iline),"<atoms>"))
      if(aurostd::substring2bool(vcontent.at(iline),"</atoms>"))
      {line=vcontent.at(iline);break;} 
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " line=" << line << endl;
  aurostd::StringSubst(line,"<atoms>","");
  aurostd::StringSubst(line,"</atoms>","");

  // ----------------------------------------------------------------------
  // LOAD NATOMS
  natoms=0.0;
  aurostd::string2tokens(line,tokens,">");
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " tokens.size()=" << tokens.size() << endl; 
  natoms=aurostd::string2utype<double>(line);
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " natoms=" << natoms << endl;

  // ----------------------------------------------------------------------
  // LOAD FORCES
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD FORCES DATA" << endl;
  vforces.clear();                    // QM FORCES calculation
  for(int iline=(int)vcontent.size()-1;iline>=0;iline--) {
    if(aurostd::substring2bool(vcontent.at(iline),"<varray name=\"forces\" >")) {
      for(uint iat=0;iat<(uint)natoms && iat<vcontent.size();iat++) {
        aurostd::StringSubst(vcontent.at(iline+iat+1),"<v>","");
        aurostd::StringSubst(vcontent.at(iline+iat+1),"</v>","");
        // if(LDEBUG) cerr << __AFLOW_FUNC__ << " vcontent.at(iline+iat+1)=" << vcontent.at(iline+iat+1) << endl;
        aurostd::string2tokens(vcontent.at(iline+iat+1),tokens," ");
        if(tokens.size()==3) {
          xvector<double> force(3);
          force[1]=aurostd::string2utype<double>(tokens.at(0));
          force[2]=aurostd::string2utype<double>(tokens.at(1));
          force[3]=aurostd::string2utype<double>(tokens.at(2));
          vforces.push_back(force); // cerr.precision(20);
          if(LDEBUG) cerr << __AFLOW_FUNC__ << " force=" << force << endl;
          // if(LDEBUG) cerr << __AFLOW_FUNC__ << " force=" << force[1]  << " " << force[2] << " " << force[3] << endl;
        } else {
          //[CO20200604 - OBSOLETE]if(!QUIET) cerr << __AFLOW_FUNC__ << " error in QM FORCES calculation" << endl;
          message << "Error in QM FORCES calculation" << endl;
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
          ERROR_flag=TRUE;
        }
      }
      iline=-1;
    }
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vforces.size()=" << vforces.size() << endl;

  // ----------------------------------------------------------------------
  // LOAD KPOINTLIST
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD KPOINTLIST DATA" << endl;
  vkpoint.clear();                    // QM KPOINTLIST calculation
  for(int iline=0;iline<(int)vcontent.size();iline++) {
    if(aurostd::substring2bool(vcontent.at(iline),"<varray name=\"kpointlist\" >")) {
      iline++;
      for(int iat=0;iline+iat<(int)vcontent.size();iat++) {
        if(aurostd::substring2bool(vcontent.at(iline+iat),"</varray>")) {
          iline=(int)vcontent.size();
        } else {
          // cerr << __AFLOW_FUNC__ << " vcontent.at(iline+iat)=" << vcontent.at(iline+iat) << endl;
          aurostd::StringSubst(vcontent.at(iline+iat),"<v>","");
          aurostd::StringSubst(vcontent.at(iline+iat),"</v>","");
          aurostd::string2tokens(vcontent.at(iline+iat),tokens," ");
          if(tokens.size()==3) {
            xvector<double> kpoint(3);
            kpoint[1]=aurostd::string2utype<double>(tokens.at(0));
            kpoint[2]=aurostd::string2utype<double>(tokens.at(1));
            kpoint[3]=aurostd::string2utype<double>(tokens.at(2));
            vkpoint.push_back(kpoint); // cerr.precision(20);
            //	    if(LDEBUG) cerr << __AFLOW_FUNC__ << " kpoint=" << kpoint << endl;
          } else {
            //[CO20200604 - OBSOLETE]if(!QUIET) cerr << __AFLOW_FUNC__ << " error in QM KPOINTLIST calculation" << endl;
            message << "Error in QM KPOINTLIST calculation" << endl;
            pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
            ERROR_flag=TRUE;
          }
        }
      }
    }
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vkpoint.size()=" << vkpoint.size() << endl;

  // ----------------------------------------------------------------------
  // LOAD WEIGHTS
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD WEIGHTS DATA" << endl;
  vweights.clear();                    // QM WEIGHTS calculation
  for(int iline=0;iline<(int)vcontent.size();iline++) {
    if(aurostd::substring2bool(vcontent.at(iline),"<varray name=\"weights\" >")) {
      iline++;
      for(int iat=0;iline+iat<(int)vcontent.size();iat++) {
        if(aurostd::substring2bool(vcontent.at(iline+iat),"</varray>")) {
          iline=(int)vcontent.size();
        } else {
          // cerr << __AFLOW_FUNC__ << " vcontent.at(iline+iat)=" << vcontent.at(iline+iat) << endl;
          aurostd::StringSubst(vcontent.at(iline+iat),"<v>","");
          aurostd::StringSubst(vcontent.at(iline+iat),"</v>","");
          aurostd::string2tokens(vcontent.at(iline+iat),tokens," ");
          if(tokens.size()==1) {
            vweights.push_back(aurostd::string2utype<double>(tokens.at(0))); // cerr.precision(20);
            //	    if(LDEBUG) cerr << __AFLOW_FUNC__ << " weight=" << weight << endl;
          } else {
            //[CO20200604 - OBSOLETE]if(!QUIET) cerr << __AFLOW_FUNC__ << " error in QM WEIGHTS calculation" << endl;
            message << "Error in QM WEIGHTS calculation" << endl;
            pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
            ERROR_flag=TRUE;
          }
        }
      }
    }
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vweights.size()=" << vweights.size() << endl;

  // ----------------------------------------------------------------------
  // LOAD STRESS
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD STRESS DATA" << endl;
  stress.clear();                    // QM STRESS calculation
  for(int iline=(int)vcontent.size()-1;iline>=0;iline--) {
    if(aurostd::substring2bool(vcontent.at(iline),"<varray name=\"stress\" >")) {
      for(uint iat=0;iat<(uint)3 && iat<vcontent.size();iat++) { // only three lines
        aurostd::StringSubst(vcontent.at(iline+iat+1),"<v>","");
        aurostd::StringSubst(vcontent.at(iline+iat+1),"</v>","");
        // if(LDEBUG) cerr << __AFLOW_FUNC__ << " vcontent.at(iline+iat+1)=" << vcontent.at(iline+iat+1) << endl;
        aurostd::string2tokens(vcontent.at(iline+iat+1),tokens," ");
        if(tokens.size()==3) {
          stress(iat+1,1)=aurostd::string2utype<double>(tokens.at(0));
          stress(iat+1,2)=aurostd::string2utype<double>(tokens.at(1));
          stress(iat+1,3)=aurostd::string2utype<double>(tokens.at(2));
        } else {
          //[CO20200604 - OBSOLETE]if(!QUIET) cerr << __AFLOW_FUNC__ << " error in QM STRESS calculation" << endl;
          message << "Error in QM STRESS calculation" << endl;
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
          ERROR_flag=TRUE;
        }
      }
      iline=-1;
    }
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " stress=" << endl << stress << endl;

  // ----------------------------------------------------------------------
  // DONE NOW RETURN  
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " END (" << time_delay(seconds) << ")" << endl;
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  //[CO20200604 - OBSOLETE]if(ERROR_flag && !QUIET) cerr << XPID << "WARNING - " << __AFLOW_FUNC__ << " ERROR_flag set in xVASPRUNXML" << endl;
  if(ERROR_flag){
    message << "ERROR_flag set in xVASPRUNXML";
    if(force_exit){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);}
    else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_,QUIET);}
    return FALSE;
  }
  return TRUE;
}

//ME20190204 BEGIN
// GetForces just grabs the forces from the vasprun file with a much faster
// algorithm than in GetProperties. This is significant for APL and AAPL which
// process many large vasprun files.
bool xVASPRUNXML::GetForces(const string& stringIN, bool QUIET) {
  stringstream sss; sss.str(stringIN);
  if (filename == "") filename="string";
  return xVASPRUNXML::GetForces(sss, QUIET);
}

bool xVASPRUNXML::GetForcesFile(const string& fileIN, bool QUIET) {
  stringstream sss;
  if (filename == "") filename = fileIN;
  aurostd::efile2stringstream(fileIN, sss);
  return xVASPRUNXML::GetForces(sss, QUIET);
}

// Cannot use const stringstream& because of std::getline
bool xVASPRUNXML::GetForces(stringstream& stringstreamIN, bool QUIET) {
  bool LDEBUG=(FALSE || XHOST.DEBUG || !QUIET);
  stringstream message;
  bool force_exit=XHOST.POSTPROCESS; //SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only
  bool ERROR_flag=FALSE;

  if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

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
          //[CO20200404 - OBSOLETE]if (!QUIET) cerr << XPID << "xVASPRUNXML::GetForces: error in QM FORCES calculation" << endl;
          //[CO20200404 - OBSOLETE]return false;
          message << "Error in QM FORCES calculation" << endl;
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
          ERROR_flag=TRUE;
        }
      }
      return true;
    }
  }
  ERROR_flag=TRUE;  //set true if you get here
  //[CO20200404 - OBSOLETE]if (!QUIET) cerr << XPID << "xVASPRUNXML::GetForces: forces not found" << endl;
  //[CO20200404 - OBSOLETE]return false;
  if(ERROR_flag){
    message << "ERROR_flag set: forces not found";
    if(force_exit){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);}
    else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_,QUIET);}
    return FALSE;
  }
  return false; //catch all false
}
//ME20190204 END

//---------------------------------------------------------------------------------
// class xIBZKPT
//---------------------------------------------------------------------------------
xIBZKPT::xIBZKPT(ostream& oss) : xStream(oss) {free();} //CO20200404 - xStream integration for logging
xIBZKPT::xIBZKPT(ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss) {free();} //CO20200404 - xStream integration for logging
xIBZKPT::xIBZKPT(const string& fileIN,bool QUIET,ostream& oss) : xStream(oss) {m_initialized=initialize(fileIN,QUIET);}  //CO20200404 - xStream integration for logging
xIBZKPT::xIBZKPT(const string& fileIN,ofstream& FileMESSAGE,bool QUIET,ostream& oss) : xStream(FileMESSAGE,oss) {m_initialized=initialize(fileIN,QUIET);}  //CO20200404 - xStream integration for logging

bool xIBZKPT::initialize(const string& fileIN,ostream& oss,bool QUIET) {xStream::initialize(oss);return initialize(fileIN,QUIET);} //CO20200508
bool xIBZKPT::initialize(const string& fileIN,ofstream& FileMESSAGE,ostream& oss,bool QUIET) {xStream::initialize(FileMESSAGE,oss);return initialize(fileIN,QUIET);} //CO20200508
bool xIBZKPT::initialize(const string& fileIN,bool QUIET) {
  free();
  filename = fileIN;
  return GetPropertiesFile(fileIN,QUIET);
}

xIBZKPT::xIBZKPT(const xIBZKPT& b) : xStream(*b.getOFStream(),*b.getOSS()) {free();copy(b);} // copy PUBLIC

xIBZKPT::~xIBZKPT() {xStream::free();free();}  //CO20200404 - xStream integration for logging

void xIBZKPT::free() {
  m_initialized=false; //CO20200404
  //------------------------------------------------------------------------------
  // GetProperties
  content="";                   // for aflowlib_libraries.cpp
  vcontent.clear();             // for aflowlib_libraries.cpp
  filename="";                  // for aflowlib_libraries.cpp
  nweights=0;                   // for aflowlib_libraries.cpp
  nkpoints_irreducible=0;       // for aflowlib_libraries.cpp
  vkpoint.clear();              // for aflowlib_libraries.cpp
  vweights.clear();             // for aflowlib_libraries.cpp
  ntetrahedra=0;                // for aflowlib_libraries.cpp
  wtetrahedra=0.0;              // for aflowlib_libraries.cpp
  vtetrahedra.clear();          // for aflowlib_libraries.cpp
}        

void xIBZKPT::copy(const xIBZKPT& b) { // copy PRIVATE
  content=b.content;
  vcontent.clear(); 
  for(uint i=0;i<b.vcontent.size();i++) vcontent.push_back(b.vcontent.at(i));  // for aflowlib_libraries.cpp
  filename=b.filename;
  nweights=b.nweights;                   // for aflowlib_libraries.cpp
  nkpoints_irreducible=b.nkpoints_irreducible;                   // for aflowlib_libraries.cpp
  vkpoint.clear(); 
  for(uint i=0;i<b.vkpoint.size();i++) vkpoint.push_back(b.vkpoint.at(i));  // for aflowlib_libraries.cpp
  vweights.clear(); for(uint i=0;i<b.vweights.size();i++) vweights.push_back(b.vweights.at(i));  // for aflowlib_libraries.cpp
  ntetrahedra=b.ntetrahedra;                   // for aflowlib_libraries.cpp
  wtetrahedra=b.wtetrahedra;                   // for aflowlib_libraries.cpp
  vtetrahedra.clear(); 
  for(uint i=0;i<b.vtetrahedra.size();i++) vtetrahedra.push_back(b.vtetrahedra.at(i));  // for aflowlib_libraries.cpp
}

const xIBZKPT& xIBZKPT::operator=(const xIBZKPT& b) {  // operator= PUBLIC
  if(this!=&b) {free();copy(b);}
  return *this;
}

void xIBZKPT::clear() {  // clear PRIVATE
  xIBZKPT _temp;
  string filename_aus=filename;
  copy(_temp);
  filename=filename_aus;
}

bool xIBZKPT::GetProperties(const string& stringIN,bool QUIET) {
  stringstream sss; sss.str(stringIN);
  if(filename=="") filename="string";
  return xIBZKPT::GetProperties(sss,QUIET);
}

bool xIBZKPT::GetPropertiesFile(const string& fileIN,bool QUIET) {
  stringstream sss;
  if(filename=="") filename=fileIN;
  aurostd::efile2stringstream(fileIN,sss);
  return xIBZKPT::GetProperties(sss,QUIET);
}

bool xIBZKPT::GetPropertiesUrlFile(const string& url,const string& file,bool QUIET) {
  //[CO20200502 - OBSOLETE]string tmpfile=XHOST.tmpfs+"/_aflow_"+XHOST.user+".pid"+XHOST.ostrPID.str()+".a"+string(AFLOW_VERSION)+".rnd"+aurostd::utype2string(uint((double) std::floor((double)100000*aurostd::ran0())))+".u"+aurostd::utype2string(uint((double) aurostd::get_useconds()))+"_"+file;
  string tmpfile=aurostd::TmpFileCreate("xIBZKPT_GetProperties"); //CO20200502 - threadID
  aurostd::url2file(url+"/"+file,tmpfile,!QUIET);
  bool out=GetPropertiesFile(tmpfile,QUIET);
  filename="url="+url;  //CO20210315
  aurostd::RemoveFile(tmpfile);
  return out;
}

bool xIBZKPT::GetProperties(const stringstream& stringstreamIN,bool QUIET) {
  bool LDEBUG=(FALSE || XHOST.DEBUG || !QUIET);
  stringstream message;
  bool force_exit=XHOST.POSTPROCESS; //SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only
  bool ERROR_flag=FALSE;

  long double seconds=aurostd::get_seconds();
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN (" << time_delay(seconds) << ")" << endl;
  clear(); // so it does not mess up vector/deque
  content=stringstreamIN.str();
  vcontent.clear();
  vector<string> tokens;
  aurostd::string2vectorstring(content,vcontent);
  string line;
  if(filename=="") filename="stringstream";
  // crunching to eat the info

  // ----------------------------------------------------------------------
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vcontent.size()=" << vcontent.size() << endl;

  // ----------------------------------------------------------------------
  // LOAD NWEIGHTS VKPOINT
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD VKPOINT DATA" << endl;
  nweights=0;
  nkpoints_irreducible=0;
  vkpoint.clear();                    // QM VKPOINT calculation
  vweights.clear();                       // QM VKPOINT calculation
  for(int iline=0;iline<(int)vcontent.size();iline++) {
    if(aurostd::substring2bool(vcontent.at(iline),"Automatically generated mesh")) {
      iline++;
      nkpoints_irreducible=aurostd::string2utype<double>(vcontent.at(iline));
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " nkpoints_irreducible=" << nkpoints_irreducible << endl;
      iline+=2; // skip text
      //  cerr << __AFLOW_FUNC__ << " vcontent.at(iline)=" << vcontent.at(iline) << endl;
      for(uint iat=0;iat<nkpoints_irreducible;iat++) {
        // 	cerr << __AFLOW_FUNC__ << " vcontent.at(iline+iat)=" << vcontent.at(iline+iat) << endl;
        aurostd::string2tokens(vcontent.at(iline+iat),tokens," ");
        if(tokens.size()==4) {
          xvector<double> kpoint(3);
          kpoint[1]=aurostd::string2utype<double>(tokens.at(0));
          kpoint[2]=aurostd::string2utype<double>(tokens.at(1));
          kpoint[3]=aurostd::string2utype<double>(tokens.at(2));
          vkpoint.push_back(kpoint); // cerr.precision(20);
          vweights.push_back(aurostd::string2utype<uint>(tokens.at(3)));
          nweights+=aurostd::string2utype<uint>(tokens.at(3));
          //	  if(LDEBUG) cerr << __AFLOW_FUNC__ << " kpoint=" << kpoint << " " << "weight=" << aurostd::string2utype<double>(tokens.at(3))<< endl;
        } else {
          //[CO20200604 - OBSOLETE]if(!QUIET) cerr << __AFLOW_FUNC__ << " error in QM NWEIGHTS/VKPOINT calculation" << endl;
          message << "Error in QM NWEIGHTS/VKPOINT calculation" << endl;
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_,QUIET);
          ERROR_flag=TRUE;
        }
      }
    }
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vkpoint.size()=" << vkpoint.size() << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vweights.size()=" << vweights.size() << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " nweights=" << nweights << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " nkpoints_irreducible=" << nkpoints_irreducible << endl;

  // ---------------------------------------------------------------------
  // LOAD NTETRAHEDRA WTETRAHEDRA TETRAHEDRA
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD TETRAHEDRA DATA" << endl;
  ntetrahedra=0;
  wtetrahedra=0.0;
  vtetrahedra.clear();                    // QM TETRAHEDRA calculation
  for(int iline=0;iline<(int)vcontent.size();iline++) {
    if(aurostd::substring2bool(vcontent.at(iline),"Tetrahedra")) {
      iline++;
      aurostd::string2tokens(vcontent.at(iline),tokens," ");
      ntetrahedra=aurostd::string2utype<uint>(tokens.at(0));
      wtetrahedra=aurostd::string2utype<double>(tokens.at(1));
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " ntetrahedra=" << ntetrahedra << endl;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " wtetrahedra=" << wtetrahedra << endl;
      iline+=1; // skip text
      for(uint iat=0;iat<ntetrahedra;iat++) {
        // 	cerr << __AFLOW_FUNC__ << " vcontent.at(iline+iat)=" << vcontent.at(iline+iat) << endl;
        aurostd::string2tokens(vcontent.at(iline+iat),tokens," ");
        if(tokens.size()==5) {
          xvector<int> tetrahedra(5);
          tetrahedra[1]=aurostd::string2utype<double>(tokens.at(0));
          tetrahedra[2]=aurostd::string2utype<double>(tokens.at(1));
          tetrahedra[3]=aurostd::string2utype<double>(tokens.at(2));
          tetrahedra[4]=aurostd::string2utype<double>(tokens.at(3));
          tetrahedra[5]=aurostd::string2utype<double>(tokens.at(4));
          vtetrahedra.push_back(tetrahedra); // cerr.precision(20);
          //	  if(LDEBUG) cerr << __AFLOW_FUNC__ << " tetrahedra=" << tetrahedra << " " << endl;
        } else {
          if(!QUIET) cerr << __AFLOW_FUNC__ << " error in QM NTETRAHEDRA/WTETRAHEDRA/TETRAHEDRA calculation" << endl;
          ERROR_flag=TRUE;
        }
      }
    }
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ntetrahedra=" << ntetrahedra << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " wtetrahedra=" << wtetrahedra << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vtetrahedra.size()=" << vtetrahedra.size() << endl;

  // ----------------------------------------------------------------------
  // DONE NOW RETURN  
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " END (" << time_delay(seconds) << ")" << endl;
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  //[CO20200604 - OBSOLETE]if(ERROR_flag && !QUIET) cerr << XPID << "WARNING - " << __AFLOW_FUNC__ << " ERROR_flag set in xIBZKPT" << endl;
  //[CO20200604 - OBSOLETE]if(ERROR_flag) return FALSE;
  if(ERROR_flag){
    message << "ERROR_flag set in xIBZKPT";
    if(force_exit){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);}
    else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_,QUIET);}
    return FALSE;
  }
  return TRUE;
}

//---------------------------------------------------------------------------------
// class xKPOINTS
//---------------------------------------------------------------------------------
xKPOINTS::xKPOINTS(ostream& oss) : xStream(oss) {free();} //CO20200404 - xStream integration for logging
xKPOINTS::xKPOINTS(ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss) {free();} //CO20200404 - xStream integration for logging
xKPOINTS::xKPOINTS(const string& fileIN,bool QUIET,ostream& oss) : xStream(oss) {m_initialized=initialize(fileIN,QUIET);}  //CO20200404 - xStream integration for logging
xKPOINTS::xKPOINTS(const string& fileIN,ofstream& FileMESSAGE,bool QUIET,ostream& oss) : xStream(FileMESSAGE,oss) {m_initialized=initialize(fileIN,QUIET);}  //CO20200404 - xStream integration for logging

bool xKPOINTS::initialize(const string& fileIN,ostream& oss,bool QUIET) {xStream::initialize(oss);return initialize(fileIN,QUIET);} //CO20200508
bool xKPOINTS::initialize(const string& fileIN,ofstream& FileMESSAGE,ostream& oss,bool QUIET) {xStream::initialize(FileMESSAGE,oss);return initialize(fileIN,QUIET);} //CO20200508
bool xKPOINTS::initialize(const string& fileIN,bool QUIET) {
  free();
  filename = fileIN;
  return GetPropertiesFile(fileIN,QUIET);
}

xKPOINTS::xKPOINTS(const xKPOINTS& b) : xStream(*b.getOFStream(),*b.getOSS()) {free();copy(b);} // copy PUBLIC

xKPOINTS::~xKPOINTS() {xStream::free();free();}  //CO20200404 - xStream integration for logging

void xKPOINTS::free() {
  m_initialized=false; //CO20200404
  //------------------------------------------------------------------------------
  // constructor
  content="";                   // for aflowlib_libraries.cpp
  vcontent.clear();             // for aflowlib_libraries.cpp
  filename="";                  // for aflowlib_libraries.cpp
  title="";                     // for aflowlib_libraries.cpp
  mode=-1;                      // for aflowlib_libraries.cpp
  grid_type="";                 // for aflowlib_libraries.cpp
  is_KPOINTS_NNN=FALSE;         // for aflowlib_libraries.cpp
  is_KPOINTS_PATH=FALSE;        // for aflowlib_libraries.cpp
  nnn_kpoints.clear();          // for aflowlib_libraries.cpp
  ooo_kpoints.clear();          // for aflowlib_libraries.cpp
  nkpoints=0;                   // for aflowlib_libraries.cpp
  path_mode="";                 // for aflowlib_libraries.cpp
  path="";                      // for aflowlib_libraries.cpp
  vpath.clear();                // for aflowlib_libraries.cpp
  vkpoints.clear();             //ME20190614
  path_grid=0;                  // for aflowlib_libraries.cpp
}        

void xKPOINTS::copy(const xKPOINTS& b) { // copy PRIVATE
  content=b.content;
  vcontent.clear(); 
  for(uint i=0;i<b.vcontent.size();i++) vcontent.push_back(b.vcontent.at(i));  // for aflowlib_libraries.cpp
  filename=b.filename;
  title=b.title;
  mode=b.mode;
  grid_type=b.grid_type;
  is_KPOINTS_NNN=b.is_KPOINTS_NNN;
  is_KPOINTS_PATH=b.is_KPOINTS_PATH;
  nnn_kpoints=b.nnn_kpoints;
  ooo_kpoints=b.ooo_kpoints;
  nkpoints=b.nkpoints;
  path_mode=b.path_mode;
  path=b.path;
  vpath=b.vpath;
  vkpoints = b.vkpoints;  //ME20190614
  path_grid=b.path_grid;
}

const xKPOINTS& xKPOINTS::operator=(const xKPOINTS& b) {  // operator= PUBLIC
  if(this!=&b) {free();copy(b);}
  return *this;
}

void xKPOINTS::clear() {  // clear PRIVATE
  xKPOINTS _temp;
  string filename_aus=filename;
  copy(_temp);
  filename=filename_aus;
}

bool xKPOINTS::GetProperties(const string& stringIN,bool QUIET) {
  stringstream sss; sss.str(stringIN);
  if(filename=="") filename="string";
  return xKPOINTS::GetProperties(sss,QUIET);
}

bool xKPOINTS::GetPropertiesFile(const string& fileIN,bool QUIET) {
  stringstream sss;
  if(filename=="") filename=fileIN;
  aurostd::efile2stringstream(fileIN,sss);
  return xKPOINTS::GetProperties(sss,QUIET);
}

bool xKPOINTS::GetPropertiesUrlFile(const string& url,const string& file,bool QUIET) {
  //[CO20200502 - OBSOLETE]string tmpfile=XHOST.tmpfs+"/_aflow_"+XHOST.user+".pid"+XHOST.ostrPID.str()+".a"+string(AFLOW_VERSION)+".rnd"+aurostd::utype2string(uint((double) std::floor((double)100000*aurostd::ran0())))+".u"+aurostd::utype2string(uint((double) aurostd::get_useconds()))+"_"+file;
  string tmpfile=aurostd::TmpFileCreate("xKPOINTS_GetProperties"); //CO20200502 - threadID
  aurostd::url2file(url+"/"+file,tmpfile,!QUIET);
  bool out=GetPropertiesFile(tmpfile,QUIET);
  filename="url="+url;  //CO20210315
  aurostd::RemoveFile(tmpfile);
  return out;
}

bool xKPOINTS::GetProperties(const stringstream& stringstreamIN,bool QUIET) {
  bool LDEBUG=(FALSE || XHOST.DEBUG || !QUIET);
  stringstream message;
  bool force_exit=XHOST.POSTPROCESS; //SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only

  bool ERROR_flag=FALSE;
  clear(); // so it does not mess up vector/deque

  long double seconds=aurostd::get_seconds();
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN (" << time_delay(seconds) << ")" << endl;

  content=stringstreamIN.str();
  vcontent.clear();
  vector<string> tokens;
  aurostd::string2vectorstring(content,vcontent);
  string line;
  if(filename=="") filename="stringstream";
  // crunching to eat the info
  title="";
  mode=-1;
  grid_type="";
  is_KPOINTS_NNN=FALSE;
  is_KPOINTS_PATH=FALSE;
  nnn_kpoints.clear(); // N*N*N
  ooo_kpoints.clear(); // ORIGIN
  nkpoints=0;
  path_mode="";
  path="";
  path_grid=0;
  vpath.clear();
  // ----------------------------------------------------------------------
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vcontent.size()=" << vcontent.size() << endl;

  // ----------------------------------------------------------------------
  // CHECK IF WITH KPOINTS NUMBERS
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK IF WITH KPOINTS NUMBERS" << endl;
  if(!is_KPOINTS_NNN && !is_KPOINTS_PATH && vcontent.size()>=5) {
    if(vcontent.at(2).at(0)=='M'||vcontent.at(2).at(0)=='m' || vcontent.at(2).at(0)=='G'||vcontent.at(2).at(0)=='g') {
      aurostd::string2tokens(vcontent.at(3),tokens);
      //    if(LDEBUG) cerr << __AFLOW_FUNC__ << " tokens.size()=" << tokens.size() << endl;
      if(tokens.size()==3) {
        aurostd::string2tokens(vcontent.at(4),tokens);
        //     if(LDEBUG) cerr << __AFLOW_FUNC__ << " tokens.size()=" << tokens.size() << endl;
        if(tokens.size()==3) {
          is_KPOINTS_NNN=TRUE;
          is_KPOINTS_PATH=FALSE;
          title=vcontent.at(0);
          mode=aurostd::string2utype<int>(vcontent.at(1));
          grid_type=vcontent.at(2);
          aurostd::string2tokens(vcontent.at(3),tokens);
          nnn_kpoints[1]=aurostd::string2utype<int>(tokens.at(0));
          nnn_kpoints[2]=aurostd::string2utype<int>(tokens.at(1));
          nnn_kpoints[3]=aurostd::string2utype<int>(tokens.at(2));
          nkpoints=nnn_kpoints[1]*nnn_kpoints[2]*nnn_kpoints[3];
          aurostd::string2tokens(vcontent.at(4),tokens);
          ooo_kpoints[1]=aurostd::string2utype<double>(tokens.at(0));
          ooo_kpoints[2]=aurostd::string2utype<double>(tokens.at(1));
          ooo_kpoints[3]=aurostd::string2utype<double>(tokens.at(2));
        }
      }
    }
  }

  // ----------------------------------------------------------------------
  // CHECK IF WITH PATH
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK IF WITH PATH" << endl;
  if(!is_KPOINTS_NNN && !is_KPOINTS_PATH && vcontent.size()>=5) {
    if(vcontent.at(2).at(0)=='L'||vcontent.at(2).at(0)=='l') {
      is_KPOINTS_NNN=FALSE;
      is_KPOINTS_PATH=TRUE;
      title=vcontent.at(0);
      mode=aurostd::string2utype<int>(vcontent.at(1));
      path_grid=aurostd::string2utype<int>(vcontent.at(1));
      grid_type=vcontent.at(2);
      path_mode=vcontent.at(3);
      xvector<double> kpt(3);  //ME20190614
      for(uint iline=4;iline<vcontent.size();iline++) {
        aurostd::StringSubst(vcontent.at(iline),"!","@");
        aurostd::StringSubst(vcontent.at(iline),"@","@ ");  //CO20210712 - for "!K" vs. "! K"
        if(aurostd::substring2bool(vcontent.at(iline),"@")) { // avoid removing ! as comment
          aurostd::string2tokens(vcontent.at(iline),tokens," ");
          if(tokens.size()>=5) {
            //	    if(LDEBUG) cerr << __AFLOW_FUNC__ << " tokens.size()=" << tokens.size() << endl;
            for(uint k=0;k<tokens.size();k++) {
              if(tokens.at(k)=="@" && k+1<tokens.size()) {
                vpath.push_back(tokens.at(k+1));
              } else if (k < 3) {  //ME20190614
                kpt[k+1] = aurostd::string2utype<double>(tokens[k]);
              }
            }
            vkpoints.push_back(kpt);  //ME20190614
          }
        }
      }
      // ----------------------------------------------------------------------
      //     if(LDEBUG) cerr << __AFLOW_FUNC__ << " vpath.size()=" << vpath.size() << endl;
      if(0) { // old
        path=vpath.at(0);
        for(uint i=1;i<vpath.size();i++) {
          if(i+1<vpath.size()) { 
            if(vpath.at(i)==vpath.at(i+1)) {
              path+="-"+vpath.at(i);
              i++;
            } else {
              path+="-"+vpath.at(i)+","+vpath.at(i+1);
              i++;
            }
          } else {
            path+="-"+vpath.at(i);
          }
        }
        // \Gamma-X,X-W,W-K,K-\Gamma,\Gamma-L,L-U,U-W,W-L,L-K,U-X
        // \Gamma-X-W-K-\Gamma-L-U-W-L-K,U-X
        // \Gamma-X-W-K-\Gamma-L-U-W-L-K,U-X
        //	if(LDEBUG) cerr << __AFLOW_FUNC__ << " path=[" << path << "]" << endl;
      }
      // ----------------------------------------------------------------------
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " vpath.size()=" << vpath.size() << endl;
      if(1) { // new
        path="";
        for(uint i=0;i<vpath.size();i+=2) 
          path+=vpath.at(i)+"-"+vpath.at(i+1)+(i+2<vpath.size()?",":"");
        // \Gamma-X,X-W,W-K,K-\Gamma,\Gamma-L,L-U,U-W,W-L,L-K,U-X
        // \Gamma-X-W-K-\Gamma-L-U-W-L-K,U-X
        // \Gamma-X-W-K-\Gamma-L-U-W-L-K,U-X
        //	if(LDEBUG) cerr << __AFLOW_FUNC__ << " path=[" << path << "]" << endl;
      }
    }
  }
  // ----------------------------------------------------------------------

  if(LDEBUG) cerr << __AFLOW_FUNC__ << " title=[" << title << "]" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " mode=" << mode << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " grid_type=[" << grid_type << "]" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " nkpoints=" << nkpoints << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " nnn_kpoints=" << nnn_kpoints << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " ooo_kpoints=" << ooo_kpoints << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " is_KPOINTS_NNN=" << is_KPOINTS_NNN << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " is_KPOINTS_PATH=" << is_KPOINTS_PATH << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " path_mode=[" << path_mode << "]" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " path=[" << path << "]" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vpath.size()=" << vpath.size() << endl;
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " vpath="; for(uint i=0;i<vpath.size();i++) cerr << vpath.at(i) << " "; cerr << endl;}
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " path_grid=" << path_grid << endl;

  // ----------------------------------------------------------------------
  // DONE NOW RETURN  
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " END (" << time_delay(seconds) << ")" << endl;
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  //[CO20200404 - OBSOLETE]if(ERROR_flag && !QUIET) cerr << XPID << "WARNING - " << __AFLOW_FUNC__ << " ERROR_flag set in xKPOINTS" << endl;
  //[CO20200404 - OBSOLETE]if(ERROR_flag) return FALSE;
  if(ERROR_flag){
    message << "ERROR_flag set in xKPOINTS";
    if(force_exit){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);}
    else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_,QUIET);}
    return FALSE;
  }
  return TRUE;
}

//ME20190623 BEGIN
ostream& operator<<(ostream& oss, const xKPOINTS& xkpts) {
  if (xkpts.is_KPOINTS_NNN) {
    oss << xkpts.title << std::endl;
    oss << xkpts.mode << std::endl;
    oss << xkpts.grid_type << std::endl;
    oss << aurostd::joinWDelimiter(xkpts.nnn_kpoints, " ") << std::endl;
    if (aurostd::iszero(xkpts.ooo_kpoints)) {
      oss << "0 0 0" << std::endl;  // No need to print all those decimal places if zero
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
    for (uint i = 0; i < xkpts.vpath.size(); i++) {
      for (int j = 1; j < 4; j++) {
        oss << std::setprecision(4) << std::fixed << std::setw(9) << xkpts.vkpoints[i][j];
      }
      oss << "  ! " << xkpts.vpath[i] << std::endl;
      if ((i % 2 == 1) && (i < xkpts.vpath.size() - 1)) oss << std::endl;
    }
  }
  return oss;
}

string xKPOINTS::createStandardTitlePath(const xstructure& xstr) {
  // Lattice part
  string stdtitle = xstr.bravais_lattice_type;
  if (stdtitle == "CUB") stdtitle += " (simple cubic) ";
  else if (stdtitle == "FCC") stdtitle += " (face-centered cubic) ";
  else if (stdtitle == "BCC") stdtitle += " (body-centered cubic) ";
  else if (stdtitle == "TET") stdtitle += " (tetragonal) ";
  else if (stdtitle == "BCT") stdtitle += " (body-centered tetragonal) ";
  else if (stdtitle == "BCT1") stdtitle += " (body-centered tetragonal c < a) ";
  else if (stdtitle == "BCT2") stdtitle += " (body-centered tetragonal a < c) ";
  else if (stdtitle == "ORC") stdtitle += " (orthorhombic) ";
  else if (stdtitle == "ORCF") stdtitle += " (face-centered orthorhombic) ";
  else if (stdtitle == "ORCF1") stdtitle += " (face-centered orthorhombic 1/a^2 > 1/b^2+1/c^2) ";
  else if (stdtitle == "ORCF2") stdtitle += " (face-centered orthorhombic 1/a^2 < 1/b^2+1/c^2) ";
  else if (stdtitle == "ORCF3") stdtitle += " (face-centered orthorhombic 1/a^2 = 1/b^2+1/c^2) ";
  else if (stdtitle == "ORCI") stdtitle += " (body-centered orthorhombic a < b < c) ";
  else if (stdtitle == "ORCC") stdtitle += " (C-centered orthorhombic a < b) ";
  else if (stdtitle == "HEX") stdtitle += " (hexagonal) ";
  else if (stdtitle == "RHL") stdtitle += " (rhombohedral) ";
  else if (stdtitle == "RHL1") stdtitle += " (rhombohedral alpha < 90) ";
  else if (stdtitle == "RHL2") stdtitle += " (rhombohedral alpha > 90) ";
  else if (stdtitle == "MCL") stdtitle += " (monoclinic) ";
  else if (stdtitle == "MCLC") stdtitle += " (C-centered monoclinic) ";
  else if (stdtitle == "MCLC1") stdtitle += " (C-centered monoclinic kgamma > 90) ";
  else if (stdtitle == "MCLC2") stdtitle += " (C-centered monoclinic kgamma = 90) ";
  else if (stdtitle == "MCLC3") stdtitle += " (C-centered monoclinic kgamma < 90, bcos(alpha)/c+(bsin(alpha)/a)^2<1) ";
  else if (stdtitle == "MCLC4") stdtitle += " (C-centered monoclinic kgamma < 90, bcos(alpha)/c+(bsin(alpha)/a)^2=1) ";
  else if (stdtitle == "MCLC5") stdtitle += " (C-centered monoclinic kgamma < 90, bcos(alpha)/c+(bsin(alpha)/a)^2>1) ";
  else if ((stdtitle == "TRI") || (aurostd::toupper(stdtitle) == "TRI1A") || (aurostd::toupper(stdtitle) == "TRI1B") || (aurostd::toupper(stdtitle) == "TRI2A") || (aurostd::toupper(stdtitle) == "TRI2B")) stdtitle += " (triclinic) ";

  // Path part
  stdtitle += vpath[0];
  for (uint i = 2; i < vpath.size(); i+= 2) {
    stdtitle += "-" + vpath[i-1];
    if (vpath[i-1] != vpath[i]) stdtitle += "|" + vpath[i];
  }
  stdtitle += "-" + vpath.back();
  return stdtitle;
}

//ME20190623 END

//---------------------------------------------------------------------------------
// class xCHGCAR
//---------------------------------------------------------------------------------
xCHGCAR::xCHGCAR(ostream& oss) : xStream(oss) {free();} //CO20200404 - xStream integration for logging
xCHGCAR::xCHGCAR(ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss) {free();} //CO20200404 - xStream integration for logging
xCHGCAR::xCHGCAR(const string& fileIN,bool QUIET,ostream& oss) : xStream(oss) {m_initialized=initialize(fileIN,QUIET);}  //CO20200404 - xStream integration for logging
xCHGCAR::xCHGCAR(const string& fileIN,ofstream& FileMESSAGE,bool QUIET,ostream& oss) : xStream(FileMESSAGE,oss) {m_initialized=initialize(fileIN,QUIET);}  //CO20200404 - xStream integration for logging

bool xCHGCAR::initialize(const string& fileIN,ostream& oss,bool QUIET) {xStream::initialize(oss);return initialize(fileIN,QUIET);} //CO20200508
bool xCHGCAR::initialize(const string& fileIN,ofstream& FileMESSAGE,ostream& oss,bool QUIET) {xStream::initialize(FileMESSAGE,oss);return initialize(fileIN,QUIET);} //CO20200508
bool xCHGCAR::initialize(const string& fileIN,bool QUIET) {
  free();
  filename = fileIN;
  return GetPropertiesFile(fileIN,QUIET);
}

xCHGCAR::xCHGCAR(const xCHGCAR& b) : xStream(*b.getOFStream(),*b.getOSS()) {free();copy(b);} // copy PUBLIC

xCHGCAR::~xCHGCAR() {xStream::free();free();}  //CO20200404 - xStream integration for logging

void xCHGCAR::free() {
  m_initialized=false; //CO20200404
  //------------------------------------------------------------------------------
  // constructor
  content="";                   // for aflowlib_libraries.cpp
  vcontent.clear();             // for aflowlib_libraries.cpp
  filename="";                  // for aflowlib_libraries.cpp
  grid.clear();                 // N*N*N triplet of grid
  vstring.clear();              // ORIGIN xvector of values
  vvalues.clear();              // ORIGIN xvector of values
  tvalues.clear();              // ORIGIN xtensor of values
}        

void xCHGCAR::copy(const xCHGCAR& b) { // copy PRIVATE
  content=b.content;
  vcontent.clear(); 
  for(uint i=0;i<b.vcontent.size();i++) vcontent.push_back(b.vcontent.at(i));  // for aflowlib_libraries.cpp
  filename=b.filename;
  grid.clear();grid=b.grid;
  vstring.clear();vstring=b.vstring;
  vvalues.clear();vvalues=b.vvalues;
  tvalues.clear();tvalues=b.tvalues;
}

const xCHGCAR& xCHGCAR::operator=(const xCHGCAR& b) {  // operator= PUBLIC
  if(this!=&b) {free();copy(b);}
  return *this;
}

void xCHGCAR::clear() {  // clear PRIVATE
  xCHGCAR _temp;
  string filename_aus=filename;
  copy(_temp);
  filename=filename_aus;
}

bool xCHGCAR::GetProperties(const string& stringIN,bool QUIET) {
  stringstream sss; sss.str(stringIN);
  if(filename=="") filename="string";
  return xCHGCAR::GetProperties(sss,QUIET);
}

bool xCHGCAR::GetPropertiesFile(const string& fileIN,bool QUIET) {
  stringstream sss;
  if(filename=="") filename=fileIN;
  aurostd::efile2stringstream(fileIN,sss);
  return xCHGCAR::GetProperties(sss,QUIET);
}

bool xCHGCAR::GetPropertiesUrlFile(const string& url,const string& file,bool QUIET) {
  //[CO20200502 - OBSOLETE]string tmpfile=XHOST.tmpfs+"/_aflow_"+XHOST.user+".pid"+XHOST.ostrPID.str()+".a"+string(AFLOW_VERSION)+".rnd"+aurostd::utype2string(uint((double) std::floor((double)100000*aurostd::ran0())))+".u"+aurostd::utype2string(uint((double) aurostd::get_useconds()))+"_"+file;
  string tmpfile=aurostd::TmpFileCreate("xCHGCAR_GetProperties"); //CO20200502 - threadID
  aurostd::url2file(url+"/"+file,tmpfile,!QUIET);
  bool out=GetPropertiesFile(tmpfile,QUIET);
  filename="url="+url;  //CO20210315
  aurostd::RemoveFile(tmpfile);
  return out;
}

bool xCHGCAR::GetProperties(const stringstream& stringstreamIN,bool QUIET) {
  bool LDEBUG=(FALSE || XHOST.DEBUG || !QUIET);
  stringstream message;
  bool force_exit=XHOST.POSTPROCESS; //SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only

  bool ERROR_flag=FALSE;
  long double seconds=aurostd::get_seconds();
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN (" << time_delay(seconds) << ")" << endl;
  clear(); // so it does not mess up vector/deque
  content=stringstreamIN.str();
  vcontent.clear();
  vector<string> tokens;
  aurostd::string2vectorstring(content,vcontent);
  string line;
  if(filename=="") filename="stringstream";

  // crunching to eat the info
  uint natoms=0;
  grid.clear(); // N*N*N
  vstring.clear(); 
  vvalues.clear(); 
  tvalues.clear();
  uint index=5;
  // ----------------------------------------------------------------------
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vcontent.size()=" << vcontent.size() << endl;
  // ----------------------------------------------------------------------
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD GRID " << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vcontent.at(" << index << ")=" << vcontent.at(index) << endl;
  aurostd::string2tokens(vcontent.at(index),tokens);
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " tokens.size()=" << tokens.size() << endl;
  for(uint i=0;i<tokens.size();i++) {
    natoms+=aurostd::string2utype<uint>(tokens.at(i));
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " natoms=" << natoms << endl;
  index+=natoms+1+2; // skip direct and atoms space and get grid
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vcontent.at(" << index << ")=" << vcontent.at(index) << endl;
  aurostd::string2tokens(vcontent.at(index),tokens);
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " tokens.size()=" << tokens.size() << endl;
  for(uint i=0;i<tokens.size();i++) {
    grid(i+1)=aurostd::string2utype<uint>(tokens.at(i));
  }
  index++;
  uint size_grid=grid(1)*grid(2)*grid(3);
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " grid=[" << grid(1) << "," << grid(2) << "," << grid(3) << "]" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " size_grid=" << size_grid << "]" << endl;
  // ----------------------------------------------------------------------
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD VSTRING " << endl;
  for(uint i=index;i<vcontent.size();i++) {
    aurostd::string2tokens(vcontent.at(i),tokens);
    //    if(LDEBUG) cerr << __AFLOW_FUNC__ << " tokens.size()=" << tokens.size() << endl;
    for(uint j=0;j<tokens.size();j++) {
      if(vstring.size()<size_grid) vstring.push_back(tokens.at(j));
    }
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vstring.size()=" << vstring.size() << endl;
  // ----------------------------------------------------------------------
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD VVALUES " << endl;
  xvector<double> vvalues_aus(vstring.size());
  for(uint i=0;i<vstring.size();i++) {
    vvalues_aus(i+1)=aurostd::string2utype<double>(vstring.at(i));
  }
  // now copy on the real vvalues which has undefined size
  vvalues=vvalues_aus;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vvalues.rows=" << vvalues.rows << endl;
  // ----------------------------------------------------------------------
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOAD TVALUES " << endl;
  //[OBSOLETE ME20180705]xtensor3<double> tvalues_aus(grid(1),grid(2),grid(3));
  xtensor<double> tvalues_aus(grid); //ME20180705
  int iii=0;
  for(int i3=1;i3<=grid(3);i3++) {  //CO - x is fastest, z is slowest
    for(int i2=1;i2<=grid(2);i2++) {
      for(int i1=1;i1<=grid(1);i1++) {  //CO - x is fastest, z is slowest
        tvalues_aus[i1][i2][i3]=vvalues(++iii); //ME20180705  //CO20200404 - using ++iii
      }
    }
  }

  tvalues=tvalues_aus;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " tvalues.shape[1]=" << tvalues.shape[1] << endl; //ME20180705
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " tvalues.shape[2]=" << tvalues.shape[2] << endl; //ME20180705
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " tvalues.shape[3]=" << tvalues.shape[3] << endl; //ME20180705
  //[OBSOLETE ME20180705]if(LDEBUG) cerr << __AFLOW_FUNC__ << " tvalues.index[1]=" << tvalues.index[1] << endl;
  //[OBSOLETE ME20180705]if(LDEBUG) cerr << __AFLOW_FUNC__ << " tvalues.index[2]=" << tvalues.index[2] << endl;
  //[OBSOLETE ME20180705]if(LDEBUG) cerr << __AFLOW_FUNC__ << " tvalues.index[3]=" << tvalues.index[3] << endl;
  // ----------------------------------------------------------------------   
  // ----------------------------------------------------------------------
  // DONE NOW RETURN  
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " END (" << time_delay(seconds) << ")" << endl;
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  //[CO20200404 - OBSOLETE]if(ERROR_flag && !QUIET) cerr << XPID << "WARNING - " << __AFLOW_FUNC__ << " ERROR_flag set in xCHGCAR" << endl;
  //[CO20200404 - OBSOLETE]if(ERROR_flag) return FALSE;
  if(ERROR_flag){
    message << "ERROR_flag set in xCHGCAR";
    if(force_exit){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);}
    else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_,QUIET);}
    return FALSE;
  }
  return TRUE;
}

//-------------------------------------------------------------------------------------------------

//CO20190803 START
//---------------------------------------------------------------------------------
// class xQMVASP
//---------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// constructor
//ME20200427 - included xStream::initialize
xQMVASP::xQMVASP(ostream& oss) : xStream(oss) {free();} //CO20200404 - xStream integration for logging
xQMVASP::xQMVASP(ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss) {free();} //CO20200404 - xStream integration for logging
xQMVASP::xQMVASP(const string& fileIN,bool QUIET,ostream& oss) : xStream(oss) {m_initialized=initialize(fileIN,QUIET);} //CO20200404 - xStream integration for logging
xQMVASP::xQMVASP(const string& fileIN,ofstream& FileMESSAGE,bool QUIET,ostream& oss) : xStream(FileMESSAGE,oss) {m_initialized = initialize(fileIN, QUIET);} //CO20200404 - xStream integration for logging

bool xQMVASP::initialize(const string& fileIN,ostream& oss,bool QUIET) {xStream::initialize(oss);return initialize(fileIN,QUIET);} //CO20200508
bool xQMVASP::initialize(const string& fileIN,ofstream& FileMESSAGE,ostream& oss,bool QUIET) {xStream::initialize(FileMESSAGE,oss);return initialize(fileIN,QUIET);} //CO20200508
bool xQMVASP::initialize(const string& fileIN, bool QUIET) {
  free();
  filename = fileIN;
  return GetPropertiesFile(fileIN, QUIET);
}

xQMVASP::xQMVASP(const xQMVASP& b) : xStream(*b.getOFStream(),*b.getOSS()) {free();copy(b);} // copy PUBLIC

xQMVASP::~xQMVASP() {xStream::free();free();} //CO20191110 //CO20200404 - xStream integration for logging

void xQMVASP::free() { //CO20191110
  m_initialized=false; //CO20200404
  //------------------------------------------------------------------------------
  content="";                   // for aflowlib_libraries.cpp
  vcontent.clear();             // for aflowlib_libraries.cpp
  filename="";                  // for aflowlib_libraries.cpp
  H_atom_relax=AUROSTD_NAN;
  H_atom_static=AUROSTD_NAN;
  vforces.clear();              //CO20191112
}

void xQMVASP::copy(const xQMVASP& b) { // copy PRIVATE //CO20191110
  xStream::copy(b);  //CO20200404 - xStream integration for logging
  free();
  m_initialized=b.m_initialized; //CO20200404 - xStream integration for logging
  content=b.content;
  vcontent.clear();for(uint i=0;i<b.vcontent.size();i++){vcontent.push_back(b.vcontent[i]);}  // for aflowlib_libraries.cpp
  filename=b.filename;
  H_atom_relax=b.H_atom_relax;
  H_atom_static=b.H_atom_static;
  vforces.clear();for(uint i=0;i<b.vforces.size();i++){vforces.push_back(b.vforces[i]);}  //CO20191112
}

const xQMVASP& xQMVASP::operator=(const xQMVASP& b) {  // operator= PUBLIC //CO20191110
  if(this!=&b) {free();copy(b);}
  return *this;
}

void xQMVASP::clear() {  // clear PRIVATE  //CO20191110
  xQMVASP _temp;
  string filename_aus=filename;
  copy(_temp);
  filename=filename_aus;
}

bool xQMVASP::GetProperties(const string& stringIN,bool QUIET) { //CO20191110
  stringstream sss; sss.str(stringIN);
  if(filename=="") filename=DEFAULT_AFLOW_QMVASP_OUT;
  return xQMVASP::GetProperties(sss,QUIET);
}

bool xQMVASP::GetPropertiesFile(const string& fileIN,bool QUIET) { //CO20191110
  stringstream sss;
  if(filename=="") filename=fileIN;
  aurostd::efile2stringstream(fileIN,sss);
  return xQMVASP::GetProperties(sss,QUIET);
}

bool xQMVASP::GetPropertiesUrlFile(const string& url,const string& file,bool QUIET) {  //CO20191110
  //[CO20200502 - OBSOLETE]string tmpfile=XHOST.tmpfs+"/_aflow_"+XHOST.user+".pid"+XHOST.ostrPID.str()+".a"+string(AFLOW_VERSION)+".rnd"+aurostd::utype2string(uint((double) std::floor((double)100000*aurostd::ran0())))+".u"+aurostd::utype2string(uint((double) aurostd::get_useconds()))+"_"+file;
  string tmpfile=aurostd::TmpFileCreate("xQMVASP_GetProperties"); //CO20200502 - threadID
  aurostd::url2file(url+"/"+file,tmpfile,!QUIET);
  bool out=GetPropertiesFile(tmpfile,QUIET);
  filename="url="+url;  //CO20210315
  aurostd::RemoveFile(tmpfile);
  return out;
}

bool xQMVASP::GetProperties(const stringstream& stringstreamIN,bool QUIET) { //CO20191110
  bool LDEBUG=(FALSE || XHOST.DEBUG || !QUIET);
  stringstream message;
  bool force_exit=XHOST.POSTPROCESS; //SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only

  bool ERROR_flag=FALSE;
  long double seconds=aurostd::get_seconds();
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN (" << time_delay(seconds) << ")" << endl;
  clear(); // so it does not mess up vector/deque
  content=stringstreamIN.str();
  vcontent.clear();
  vector<string> tokens,tokens2;
  aurostd::string2vectorstring(content,vcontent);
  if(filename=="") filename="stringstream";

  // crunching to eat the info
  H_atom_relax=AUROSTD_NAN; //this will be the LAST relax
  H_atom_static=AUROSTD_NAN;
  bool inside_relax=false,inside_static=false;
  // ----------------------------------------------------------------------
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vcontent.size()=" << vcontent.size() << endl;
  // ----------------------------------------------------------------------
  for(uint iline=0;iline<vcontent.size();iline++){
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " vcontent.at(" << iline << ")=" << vcontent.at(iline) << endl;
    if(aurostd::substring2bool(vcontent[iline],"STOP_relax")){inside_relax=false;}
    else if(aurostd::substring2bool(vcontent[iline],"START_relax")){inside_relax=true;}
    else if(aurostd::substring2bool(vcontent[iline],"STOP_static")){inside_static=false;}
    else if(aurostd::substring2bool(vcontent[iline],"START_static")){inside_static=true;}
    else {
      if(aurostd::substring2bool(vcontent[iline],"H_atom")){
        //H_atom=-9.859773400000e+00  (eV/at)
        tokens.clear();tokens2.clear();
        aurostd::string2tokens(vcontent[iline],tokens,"=");
        if(tokens.size()){aurostd::string2tokens(tokens[1],tokens2," ");}
        if(tokens2.size()){
          if(!aurostd::isfloat(tokens2[0])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"H_atom input cannot be parsed (!isfloat)",_INPUT_ERROR_);}
          if(inside_relax){H_atom_relax=aurostd::string2utype<double>(tokens2[0]);}
          else if(inside_static){H_atom_static=aurostd::string2utype<double>(tokens2[0]);}
        }
        //ME20191218 - substring2bool ignores comments and TOTAL-FORCE is in a comment line
      }
      // [OBSOLETE - ME20191219] else if(aurostd::substring2bool(vcontent[iline],"TOTAL-FORCE")) //CO20191112 - forces for APL
      else if(vcontent[iline].find("TOTAL-FORCE") != string::npos)
      { //CO20200106 - patching for auto-indenting
        vforces.clear();
        iline++;  //skip first [AFLOW]
        //ME20191219 - ++iline needs to be in the while statement or the loop
        // will never start
        while(iline<vcontent.size() && aurostd::substring2bool(vcontent[++iline],"[AFLOW]")==FALSE){
          vforces.push_back(xvector<double>(3));
          aurostd::string2tokens(vcontent[iline],tokens," ");
          if(tokens.size()==6){
            for(int i=1;i<4;i++){
              if(aurostd::isfloat(tokens[i+2])){
                vforces.back()[i]=aurostd::string2utype<double>(tokens[i+2]);
              } else {
                throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Expected force input to be a number",_FILE_CORRUPT_);
              }
            }
          } else {
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unexpected count of force components",_FILE_CORRUPT_);
          }
        }
      }
    }
  }

  // ----------------------------------------------------------------------   
  // ----------------------------------------------------------------------
  // DONE NOW RETURN  
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " END (" << time_delay(seconds) << ")" << endl;
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  //[CO20200404 - OBSOLETE]if(ERROR_flag && !QUIET) cerr << "WARNING - " << __AFLOW_FUNC__ << " ERROR_flag set in xQMVASP" << endl;
  if(ERROR_flag){
    message << "ERROR_flag set in xQMVASP";
    if(force_exit){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);}
    else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_,QUIET);}
    return FALSE;
  }
  return TRUE;
}

//-------------------------------------------------------------------------------------------------
//CO20190803 STOP

//CO20190803 START
//---------------------------------------------------------------------------------
// class xPLASMONICS
//---------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// constructor
xPLASMONICS::xPLASMONICS(ostream& oss) : xStream(oss) {free();} //CO20200404 - xStream integration for logging
xPLASMONICS::xPLASMONICS(ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss) {free();} //CO20200404 - xStream integration for logging
xPLASMONICS::xPLASMONICS(const string& fileIN,bool QUIET,ostream& oss) : xStream(oss) {m_initialized=initialize(fileIN,QUIET);} //CO20200404 - xStream integration for logging
xPLASMONICS::xPLASMONICS(const string& fileIN,ofstream& FileMESSAGE,bool QUIET,ostream& oss) : xStream(FileMESSAGE,oss) {m_initialized = initialize(fileIN, QUIET);} //CO20200404 - xStream integration for logging

bool xPLASMONICS::initialize(const string& fileIN,ostream& oss,bool QUIET) {xStream::initialize(oss);return initialize(fileIN,QUIET);} //CO20200508
bool xPLASMONICS::initialize(const string& fileIN,ofstream& FileMESSAGE,ostream& oss,bool QUIET) {xStream::initialize(FileMESSAGE,oss);return initialize(fileIN,QUIET);} //CO20200508
bool xPLASMONICS::initialize(const string& fileIN, bool QUIET) {
  free();
  filename = fileIN;
  return GetPropertiesFile(fileIN, QUIET);
}

xPLASMONICS::xPLASMONICS(const xPLASMONICS& b) : xStream(*b.getOFStream(),*b.getOSS()) {free();copy(b);} // copy PUBLIC

xPLASMONICS::~xPLASMONICS() {xStream::free();free();} //CO20191110 //CO20200404 - xStream integration for logging

void xPLASMONICS::free() { //CO20191110
  m_initialized=false; //CO20200404
  //------------------------------------------------------------------------------
  content="";                   // for aflowlib_libraries.cpp
  vcontent.clear();             // for aflowlib_libraries.cpp
  filename="";                  // for aflowlib_libraries.cpp
  eps.clear();
  venergy.clear();
  veels.clear();
  vdielectric.clear();
}

void xPLASMONICS::copy(const xPLASMONICS& b) { // copy PRIVATE //CO20191110
  xStream::copy(b);  //CO20200404 - xStream integration for logging
  free();
  m_initialized=b.m_initialized; //CO20200404 - xStream integration for logging
  content=b.content;
  vcontent.clear();for(uint i=0;i<b.vcontent.size();i++){vcontent.push_back(b.vcontent[i]);}  // for aflowlib_libraries.cpp
  filename=b.filename;
  eps=b.eps;
  venergy.clear();for(uint i=0;i<b.venergy.size();i++){venergy.push_back(b.venergy[i]);}
  veels.clear();for(uint i=0;i<b.veels.size();i++){veels.push_back(b.veels[i]);}
  vdielectric.clear();for(uint i=0;i<b.vdielectric.size();i++){vdielectric.push_back(b.vdielectric[i]);}
}

const xPLASMONICS& xPLASMONICS::operator=(const xPLASMONICS& b) {  // operator= PUBLIC //CO20191110
  if(this!=&b) {free();copy(b);}
  return *this;
}

void xPLASMONICS::clear() {  // clear PRIVATE  //CO20191110
  string filename_aus=filename;
  free();
  filename=filename_aus;
}

void xPLASMONICS::getEPS() {  //CO20211120
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " filename=" << filename << endl;}
  if(filename.find(DEFAULT_AFLOW_PLASMONICS_FILE)!=string::npos){
    //get eps and store it
    string _eps=filename;
    _eps=aurostd::basename(_eps);
    aurostd::StringSubst(_eps,DEFAULT_AFLOW_PLASMONICS_FILE+"_",""); //remove DEFAULT_AFLOW_PLASMONICS_FILE+'_'
    aurostd::StringSubst(_eps,DEFAULT_AFLOW_PLASMONICS_FILE,""); //remove DEFAULT_AFLOW_PLASMONICS_FILE
    _eps=_eps.substr(0,_eps.find("_")); //remove everything before '_'
    _eps=_eps.substr(0,_eps.find("."+POCC_OUT_FILE)); //remove everything before '.out'
    _eps=_eps.substr(0,_eps.find(POCC_OUT_FILE)); //remove everything before 'out'
    if(aurostd::isfloat(_eps)){eps=_eps;}
  }
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " eps=" << eps << endl;}
}

bool xPLASMONICS::GetProperties(const string& stringIN,bool QUIET) { //CO20191110
  stringstream sss; sss.str(stringIN);
  if(filename=="") filename=DEFAULT_AFLOW_PLASMONICS_FILE+"_0.01."+POCC_OUT_FILE;
  return xPLASMONICS::GetProperties(sss,QUIET);
}

bool xPLASMONICS::GetPropertiesFile(const string& fileIN,bool QUIET) { //CO20191110
  stringstream sss;
  //[CO20211121 - getEPS() needs the right filename]if(filename=="") 
  filename=fileIN;
  aurostd::efile2stringstream(fileIN,sss);
  return xPLASMONICS::GetProperties(sss,QUIET);
}

bool xPLASMONICS::GetPropertiesUrlFile(const string& url,const string& file,bool QUIET) {  //CO20191110
  //[CO20200502 - OBSOLETE]string tmpfile=XHOST.tmpfs+"/_aflow_"+XHOST.user+".pid"+XHOST.ostrPID.str()+".a"+string(AFLOW_VERSION)+".rnd"+aurostd::utype2string(uint((double) std::floor((double)100000*aurostd::ran0())))+".u"+aurostd::utype2string(uint((double) aurostd::get_useconds()))+"_"+file;
  string tmpfile=aurostd::TmpFileCreate("xPLASMONICS_GetProperties"); //CO20200502 - threadID
  aurostd::url2file(url+"/"+file,tmpfile,!QUIET);
  bool out=GetPropertiesFile(tmpfile,QUIET);
  filename="url="+url;  //CO20210315
  aurostd::RemoveFile(tmpfile);
  return out;
}

bool xPLASMONICS::GetProperties(const stringstream& stringstreamIN,bool QUIET) { //CO20191110
  bool LDEBUG=(FALSE || XHOST.DEBUG || !QUIET);
  stringstream message;
  bool force_exit=XHOST.POSTPROCESS; //SC wants to exit here so we can fix the problem  // ME20200604 - do not exit with generate_aflowin_only

  bool ERROR_flag=FALSE;
  long double seconds=aurostd::get_seconds();
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN (" << time_delay(seconds) << ")" << endl;
  clear(); // so it does not mess up vector/deque
  content=stringstreamIN.str();
  vcontent.clear();
  vector<string> vlines,tokens,tokens2;
  aurostd::string2vectorstring(content,vcontent);
  vlines=aurostd::RemoveComments(vcontent);
  if(filename=="") filename="stringstream";

  // ----------------------------------------------------------------------
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vcontent.size()=" << vcontent.size() << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " vlines.size()=" << vlines.size() << endl;
  // ----------------------------------------------------------------------

  getEPS(); //extract eps from filename if possible

  // crunching to eat the info
  vector<string> vtokens;
  xcomplex<double> xcomp_tmp;
  venergy.clear();veels.clear();vdielectric.clear();  //clear
  for(uint iline=0;iline<vlines.size();iline++){
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " vlines[iline=" << iline << "]=\"" << vlines[iline] << "\"" << endl;}
    aurostd::string2tokens(vlines[iline],vtokens," ");
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " vtokens=" << aurostd::joinWDelimiter(vtokens,"|") << endl;}
    //there should be 5 columns
    if(vtokens.size()!=5){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"vcols!=5 (="+aurostd::utype2string(vtokens.size())+")",_FILE_NOT_FOUND_);}
    //The file has 5 columns:
    //column1 = energy (eV)
    //column2 = do not consider this
    //column3 = EELS spectrum
    //column4 = real part of dielectric function
    //column5 = imaginary part of dielectric function
    venergy.push_back(aurostd::string2utype<double>(vtokens[0]));
    veels.push_back(aurostd::string2utype<double>(vtokens[2]));
    xcomp_tmp.re=aurostd::string2utype<double>(vtokens[3]);
    xcomp_tmp.im=aurostd::string2utype<double>(vtokens[4]);
    vdielectric.push_back(xcomp_tmp);
  }

  // ----------------------------------------------------------------------   
  // ----------------------------------------------------------------------
  // DONE NOW RETURN  
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " END (" << time_delay(seconds) << ")" << endl;
  // ----------------------------------------------------------------------
  // DONE NOW RETURN
  //[CO20200404 - OBSOLETE]if(ERROR_flag && !QUIET) cerr << "WARNING - " << __AFLOW_FUNC__ << " ERROR_flag set in xPLASMONICS" << endl;
  if(ERROR_flag){
    message << "ERROR_flag set in xPLASMONICS";
    if(force_exit){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);}
    else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_,QUIET);}
    return FALSE;
  }
  return TRUE;
}

//-------------------------------------------------------------------------------------------------
//CO20190803 STOP


#endif //  _AFLOW_OVASP_CPP_
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
