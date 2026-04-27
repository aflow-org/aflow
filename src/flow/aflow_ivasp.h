
#ifndef AFLOW_IVASP_H
#define AFLOW_IVASP_H

#include <deque>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xvector.h"

#include "flow/aflow_xclasses.h"

namespace KBIN {

  bool VASP_Produce_INPUT(_xvasp& xvasp, const std::string& AflowIn, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags, _vflags& vflags, bool load_POSCAR_from_xvasp = false);
  bool VASP_Modify_INPUT(_xvasp& xvasp, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags, _vflags& vflags);
  bool VASP_Produce_and_Modify_INPUT(_xvasp& xvasp, const std::string& AflowIn, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags, _vflags& vflags, bool load_POSCAR_from_xvasp = false); // CO20180418
  bool VASP_Write_ppAUID_FILE(const std::string& directory, const std::vector<std::string>& vppAUIDs, const std::vector<std::string>& species);
  bool VASP_Write_ppAUID_FILE(const std::string& directory, const std::deque<std::string>& vppAUIDs, const std::deque<std::string>& species);
  bool VASP_Write_ppAUID_AFLOWIN(const std::string& directory, const std::vector<std::string>& vppAUIDs, const std::vector<std::string>& species);
  bool VASP_Write_ppAUID_AFLOWIN(const std::string& directory, const std::deque<std::string>& vppAUIDs, const std::deque<std::string>& species);
  bool VASP_Write_INPUT(_xvasp& xvasp, _vflags& vflags, const std::string& ext_module = ""); // AS20210302
  bool VASP_Produce_INCAR(_xvasp& xvasp, const std::string& AflowIn, std::ofstream& FileERROR, _aflags& aflags, _kflags& kflags, _vflags& vflags);
  bool VASP_Modify_INCAR(_xvasp& xvasp, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags, _vflags& vflags);
  void VASP_CleanUp_INCAR(_xvasp& xvasp);
  bool VASP_Reread_INCAR(_xvasp& xvasp); // CO20210315
  bool VASP_Reread_INCAR(_xvasp& xvasp, std::ofstream& FileMESSAGE, _aflags& aflags);
  bool VASP_Produce_POSCAR(_xvasp& xvasp, const std::string& AflowIn, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags, _vflags& vflags);
  bool VASP_Produce_POSCAR(_xvasp& xvasp);
  bool VASP_Modify_POSCAR(_xvasp& xvasp, const std::string& AflowIn, std::ofstream& FileMESSAGE, _aflags& aflags, _vflags& vflags);
  void convertPOSCARFormat(_xvasp&, const _kflags&); // ME20190220 //CO20210713 - aflags //SD20220923 - updated for vasp6.3
  bool VASP_Convert_Unit_Cell(_xvasp&, _vflags&, _aflags&, std::ofstream&, std::ostringstream&); // ME20181216
  bool VASP_Reread_POSCAR(_xvasp& xvasp); // CO20210315
  bool VASP_Reread_POSCAR(_xvasp& xvasp, std::ofstream& FileMESSAGE, _aflags& aflags);
  bool VASP_Produce_KPOINTS(_xvasp& xvasp, const std::string& AflowIn, std::ofstream& FileERROR, _aflags& aflags, _kflags& kflags, _vflags& vflags);
  bool VASP_Modify_KPOINTS(_xvasp& xvasp, std::ofstream& FileERROR, _aflags& aflags, _vflags& vflags);
  bool VASP_Reread_KPOINTS(_xvasp& xvasp); // CO20210315
  bool VASP_Reread_KPOINTS(_xvasp& xvasp, std::ofstream& FileMESSAGE, _aflags& aflags);
  bool VASP_Find_DATA_POTCAR(const std::string& species_pp, std::string& FilePotcar, std::string& DataPotcar, std::string& AUIDPotcar);
  bool VASP_Find_FILE_POTCAR(const std::string& species_pp, std::string& FilePotcar, std::string& DataPotcar, std::string& AUIDPotcar);
  bool VASP_Produce_POTCAR(_xvasp& xvasp, const std::string& AflowIn, std::ofstream& FileERROR, _aflags& aflags, _kflags& kflags, _vflags& vflags);
  bool VASP_Modify_POTCAR(_xvasp& xvasp, std::ofstream& FileERROR, _aflags& aflags, _vflags& vflags);
  bool VASP_Reread_POTCAR(_xvasp& xvasp, std::ofstream& FileMESSAGE, _aflags& aflags);
  bool VASP_PseudoPotential_CleanName_TEST(); // CO20190712
  uint VASP_SplitAlloySpecies(std::string alloy_in, std::vector<std::string>& speciesX);
  uint VASP_SplitAlloySpecies(std::string alloy_in, std::vector<std::string>& speciesX, std::vector<double>& natomsX);
  bool VASP_SplitAlloySpecies(std::string alloy_in, std::string& specieA, std::string& specieB);
  bool VASP_SplitAlloySpecies(std::string alloy_in, std::string& specieA, std::string& specieB, std::string& specieC);
  bool VASP_SplitAlloySpecies(std::vector<std::string> alloy, std::vector<std::string>& speciesA, std::vector<std::string>& speciesB);
  bool VASP_SplitAlloySpecies(std::vector<std::string> alloy, std::vector<std::string>& speciesA, std::vector<std::string>& speciesB, std::vector<std::string>& speciesC);
  uint VASP_SplitAlloyPseudoPotentials(std::string alloy_in, std::vector<std::string>& species_ppX);
  uint VASP_SplitAlloyPseudoPotentials(std::string alloy_in, std::vector<std::string>& species_ppX, std::vector<double>& natomsX);
  bool VASP_SplitAlloyPseudoPotentials(std::string alloy, std::string& species_ppA, std::string& species_ppB);
  bool VASP_SplitAlloyPseudoPotentials(std::string alloy, std::string& species_ppA, std::string& species_ppB, std::string& species_ppC);
  bool VASP_SplitAlloyPseudoPotentials(std::vector<std::string> alloy, std::vector<std::string>& species_ppsA, std::vector<std::string>& species_ppsB);
  bool VASP_SplitAlloyPseudoPotentials(std::vector<std::string> alloy, std::vector<std::string>& species_ppsA, std::vector<std::string>& species_ppsB, std::vector<std::string>& species_ppsC);
  void XVASP_Get_NPAR_NCORE(const _xvasp& xvasp, const _aflags& aflags, int& NPAR, int& NCORE); // CO20210315
  void XVASP_MPI_Autotune(_xvasp& xvasp, _aflags& aflags, bool VERBOSE);
  void XVASP_INCAR_System_Auto(_xvasp& xvasp, bool VERBOSE);
  void XVASP_INCAR_Relax_ON(_xvasp& xvasp, bool VERBOSE);
  void XVASP_INCAR_Relax_ON(_xvasp& xvasp, _vflags& vflags, int number); // for steps
  void XVASP_INCAR_Static_ON(_xvasp& xvasp, _vflags& vflags);
  void XVASP_INCAR_Relax_Static_ON(_xvasp& xvasp, _vflags& vflags);
  void XVASP_INCAR_Bands_ON(_xvasp& xvasp, _vflags& vflags);
  void XVASP_INCAR_Dielectric_ON(_xvasp& xvasp, _vflags& vflags, _aflags& aflags);
  void XVASP_INCAR_RWIGS_Static(_xvasp& xvasp, _vflags& vflags, std::ofstream& FileMESSAGE, bool OPERATION);
  void XVASP_INCAR_Precision(_xvasp& xvasp, _vflags& vflags);
  void XVASP_INCAR_Metagga(_xvasp& xvasp, _vflags& vflags);
  void XVASP_INCAR_Ivdw(_xvasp& xvasp, _vflags& vflags);
  void XVASP_INCAR_ABMIX(_xvasp& xvasp, _vflags& vflags);
  int XVASP_INCAR_GetNBANDS_AFLOW3(const _xvasp& xvasp, const _aflags& aflags, bool ispin = true); // CO20210315 - spin==true is safer
  int XVASP_INCAR_GetNBANDS(const _xvasp& xvasp, const _aflags& aflags, bool ispin = true);
  std::string INCAR_IALGO2ALGO(int ialgo); // CO20210315
  bool XVASP_INCAR_Read_MAGMOM(_xvasp& xvasp); // CO20210315
  bool XVASP_INCAR_PREPARE_GENERIC(const std::string& command, _xvasp& xvasp, const _vflags& vflags, const std::string& svalue, int ivalue, double dvalue, bool bvalue);
  void XVASP_INCAR_ADJUST_ICHARG(_xvasp&, _vflags&, _aflags&, int, bool write_incar, std::ofstream&); // ME20191028 //CO20210315 - write_incar
  void XVASP_INCAR_SPIN_REMOVE_RELAX(_xvasp& xvasp, _aflags& aflags, _vflags& vflags, int step, bool write_incar, std::ofstream& FileMESSAGE); // CO20210315 - write_incar
  void XVASP_KPOINTS_IBZKPT_UPDATE(_xvasp& xvasp, _aflags& aflags, _vflags& vflags, int step, bool write_incar, std::ofstream& FileMESSAGE); // CO20210315 - write_incar
  void XVASP_INCAR_LDAU_OFF(_xvasp& xvasp, bool VERBOSE);
  void XVASP_INCAR_LDAU_ON(_xvasp& xvasp, _vflags& vflags, uint type);
  void XVASP_INCAR_LDAU_ADIABATIC(_xvasp& xvasp, int step);
  void XVASP_INCAR_LDAU_CUTOFF(_xvasp& xvasp, bool VERBOSE);
  void XVASP_INCAR_REMOVE_ENTRY(_xvasp& xvasp, const std::string& entry, const std::string& COMMENT, bool VERBOSE); // CO20200624 - using aurostd::kvpair2bool() now
  void XVASP_INCAR_REMOVE_ENTRY(_xvasp& xvasp, const std::vector<std::string>& entries, const std::string& COMMENT, bool VERBOSE); // CO20200624 - using aurostd::kvpair2bool() now

  bool AFLOWIN_REMOVE(const std::string& aflowin_file, const std::string& keyword, const std::string& comment); // CO20210314
  bool AFLOWIN_REMOVE(const std::string& aflowin_file, const std::vector<std::string>& vkeywords, const std::string& comment); // CO20210314
  bool AFLOWIN_REMOVE(const std::string& aflowin_file, const std::string& keyword, const std::string& keyword2avoid, const std::string& comment); // CO20210314
  bool AFLOWIN_REMOVE(const std::string& aflowin_file, const std::vector<std::string>& vkeywords, const std::vector<std::string>& vkeywords2ignore, const std::string& comment); // CO20210314
  void AFLOWIN_ADD(const std::string& aflowin_file, const std::stringstream& streamin, const std::string& comment); // CO20210315
  void AFLOWIN_ADD(const std::string& aflowin_file, const std::ostringstream& streamin, const std::string& comment); // CO20210315
  void AFLOWIN_ADD(const std::string& aflowin_file, const std::string& line, const std::string& comment); // CO20210315
  void AFLOWIN_ADD(const std::string& aflowin_file, const std::vector<std::string>& vlines2add, const std::string& comment); // CO20210315

  void XVASP_KPOINTS_KPOINTS(_xvasp& xvasp, std::ofstream& FileMESSAGE, bool VERBOSE);
  void XVASP_KPOINTS_KPOINTS(_xvasp& xvasp);
  // bool XVASP_KPOINTS_EVEN(_xvasp& xvasp); TO REMOVE
  // bool XVASP_KPOINTS_ODD(_xvasp& xvasp); TO REMOVE
  bool XVASP_KPOINTS_IncludesGamma(const _xvasp& xvasp); // CO20210315
  bool XVASP_KPOINTS_OPERATION(_xvasp& xvasp, const std::string& operation);
  // bool XVASP_KPOINTS_Kshift_Gamma_EVEN(_xvasp& xvasp); TO REMOVE
  // bool XVASP_KPOINTS_Kshift_Gamma_ODD(_xvasp& xvasp); TO REMOVE
  // bool XVASP_KPOINTS_Kscheme(_xvasp& xvasp,string kscheme);
  bool XVASP_KPOINTS_isAutoMesh(const _xvasp& xvasp); // CO20210315
  bool XVASP_KPOINTS_string2numbers(_xvasp& xvasp); // CO20210315 - cleaned up
  bool XVASP_KPOINTS_numbers2string(_xvasp& xvasp); // CO20210315 - cleaned up
  void XVASP_Afix_Clean(const _xvasp& xvasp, const std::string& preserve_name);
  // the following functions are all associated with XVASP_Afix()
  bool XVASP_Afix_NBANDS(_xvasp& xvasp, const _aflags& aflags, const _vflags& vflags, bool increase = true); // CO20200624 - this is not a general Afix, this can only be used inside Afix_GENERIC
  bool XVASP_Afix_NBANDS(_xvasp& xvasp, const _aflags& aflags, const _vflags& vflags, int& nbands, bool increase = true); // CO20200624 - this is not a general Afix, this can only be used inside Afix_GENERIC
  bool XVASP_Afix_POTIM(_xvasp& xvasp, const _vflags& vflags); // CO20200624 - this is not a general Afix, this can only be used inside Afix_GENERIC
  bool XVASP_Afix_POTIM(_xvasp& xvasp, const _vflags& vflags, double& potim); // CO20200624 - this is not a general Afix, this can only be used inside Afix_GENERIC
  bool XVASP_Afix_NELM(_xvasp& xvasp, const _vflags& vflags); // CO20200624 - this is not a general Afix, this can only be used inside Afix_GENERIC
  bool XVASP_Afix_NELM(_xvasp& xvasp, const _vflags& vflags, int& nelm); // CO20200624 - this is not a general Afix, this can only be used inside Afix_GENERIC
  bool XVASP_Afix_EFIELD_PEAD(_xvasp& xvasp, const _vflags& vflags); // CO20200624 - this is not a general Afix, this can only be used inside Afix_GENERIC
  bool XVASP_Afix_EFIELD_PEAD(_xvasp& xvasp, const _vflags& vflags, aurostd::xvector<double>& E); // CO20200624 - this is not a general Afix, this can only be used inside Afix_GENERIC
  bool XVASP_Afix_ApplyFix(const std::string& fix, aurostd::xoption& xfixed, _xvasp& xvasp, _kflags& kflags, _vflags& vflags, _aflags& aflags, std::ofstream& FileMESSAGE); // CO20200624 - adding submode so we don't need to make a bunch of spin-off functions
  bool XVASP_Afix_IgnoreFix(const std::string& _fix, const _vflags& vflags); // CO20210315
  bool XVASP_Afix(const std::string& mode, int& submode, bool try_last_ditch_efforts, aurostd::xoption& xfixed, _xvasp& xvasp, _kflags& kflags, _vflags& vflags, _aflags& aflags, std::ofstream& FileMESSAGE); // CO20200624 - adding submode so we don't need to make a bunch of spin-off functions

  std::string ExtractSystemName(const std::string& directory); // ME20200217
  std::string ExtractSystemNameFromAFLOWIN(const std::string& directory); // ME20200217
  std::string ExtractSystemNameFromVASP(const std::string& directory); // ME20200217
  xstructure GetPOSCARXStructureFromDirectory(const std::string& directory, const int iomode, const int index); // SD20220228
  bool ExtractPOSCARStringStreamFromDirectory(const std::string& directory, std::stringstream& poscar, const int index); // SD20220228
  xstructure GetPOSCARXStructureFromAFLOWIN(const std::string& AflowIn, const int iomode, const int index); // SD20220228
  bool ExtractPOSCARStringStreamFromAFLOWIN(const std::string& AflowIn, std::stringstream& poscar, const int index); // SD20220228
  double ExtractEfermiOUTCAR(std::string directory);
  xstructure GetMostRelaxedStructure(std::string directory); // CO20180627
  std::vector<std::string> ExtractAtomicSpecies(const std::string& directory);
  std::vector<std::string> ExtractAtomicSpecies(const std::string& directory, std::ofstream& FileMESSAGE);

} // namespace KBIN

#endif // AFLOW_IVASP_H
