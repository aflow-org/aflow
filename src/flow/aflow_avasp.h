
#ifndef AFLOW_AVASP_H
#define AFLOW_AVASP_H

#include <deque>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xoption.h"

#include "flow/aflow_xclasses.h"

struct _AVASP_PROTO {
  std::vector<std::string> ucell;
  std::deque<int> vkppra;
  std::vector<double> vpressure;
  aurostd::xoption vparams;
};

void PARAMS2xvasp(_AVASP_PROTO *PARAMS, _xvasp &xvasp);  // CO20210624 - avoid duplicate code: AVASP_MakePrototype_AFLOWIN() and AVASP_MakePrototypeICSD_AFLOWIN()
bool AVASP_MakePrototype_AFLOWIN(_AVASP_PROTO *PARAMS);
bool AVASP_MakePrototypeICSD_AFLOWIN(_AVASP_PROTO *PARAMS, bool flag_AFLOW_IN_ONLY_IF_MISSING);
void AVASP_Get_LDAU_Parameters(std::string species, bool &LDAU, std::vector<std::string> &vLDAUspecies, std::vector<uint> &vLDAUtype, std::vector<int> &vLDAUL, std::vector<double> &vLDAUU, std::vector<double> &vLDAUJ);
std::string AVASP_Get_PseudoPotential_PAW_PBE_KIN(std::string species);
std::string AVASP_Get_PseudoPotential_PAW_PBE(std::string species);
std::string AVASP_Get_PseudoPotential_PAW_GGA(std::string species);
std::string AVASP_Get_PseudoPotential_PAW_LDA_KIN(std::string species);
std::string AVASP_Get_PseudoPotential_PAW_LDA(std::string species);
std::string AVASP_Get_PseudoPotential_PBE(std::string species);
std::string AVASP_Get_PseudoPotential_GGA(std::string species);
std::string AVASP_Get_PseudoPotential_LDA(std::string species);
bool AVASP_populateXVASP(const _aflags &aflags, const _kflags &kflags, const _vflags &vflags, _xvasp &xvasp);
void AVASP_populateXVASP_ARUN(const _aflags &, const _kflags &, const _vflags &, _xvasp &);  // ME20181030
void setStatic(_xvasp &);  // ME20181102
void setPreserveUnitCell(_xvasp &);   // ME20181102
void AVASP_fix_volumes_masses_XVASP(_xvasp &, bool skip_volume = false); // ME20181103 //CO20181226
bool AVASP_MakeSingleAFLOWIN(_xvasp &xvasp_in, std::stringstream &_aflowin, bool flag_WRITE, int = -1, bool flag_PRINT = true);   // last is pthread number, if <0 then serial
bool AVASP_MakeSingleAFLOWIN(_xvasp &xvasp_in, bool flag_WRITE, int = -1, bool flag_PRINT = true);  // last is pthread number, if <0 then serial
bool AVASP_MakeSingleAFLOWIN(_xvasp &xvasp_in, int = -1, bool flag_PRINT = true);  // last is pthread number, if <0 then serial
bool AVASP_DefaultValuesBinary_AFLOWIN(_xvasp &xvasp);
bool AVASP_MakeSinglePOSCAR(_xvasp &xvaspin);
bool Alloys_LibraryU(std::vector<std::string> &alloy, std::vector<std::string> &pseudosA, std::vector<std::string> &pseudosB);
bool Alloys_LibraryG(std::vector<std::string> &alloy, std::vector<std::string> &pseudosA, std::vector<std::string> &pseudosB);
bool Alloys_LibraryX(std::vector<std::string> &alloy, std::vector<std::string> &pseudosA, std::vector<std::string> &pseudosB);

#endif // AFLOW_AVASP_H
