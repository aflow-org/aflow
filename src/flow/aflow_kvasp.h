
#ifndef AFLOW_KVASP_H
#define AFLOW_KVASP_H

#include <cstddef>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
#include <unordered_map>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xvector.h"

#include "flow/aflow_xclasses.h"

namespace KBIN {

  _kflags VASP_Get_Kflags_from_AflowIN(const std::string &AflowIn, _aflags &aflags, std::ostream &oss = std::cout);
  _kflags VASP_Get_Kflags_from_AflowIN(const std::string &AflowIn, std::ofstream &FileMESSAGE, _aflags &aflags, std::ostream &oss = std::cout);
  _vflags VASP_Get_Vflags_from_AflowIN(const std::string &AflowIn, _aflags &aflags, _kflags &kflags, std::ostream &oss = std::cout);
  _vflags VASP_Get_Vflags_from_AflowIN(const std::string &AflowIn, std::ofstream &FileMESSAGE, _aflags &aflags, _kflags &kflags, std::ostream &oss = std::cout);
  bool VASP_Fix_Machine_Kflags_from_AflowIN(std::ofstream &FileMESSAGE, _aflags &aflags, _kflags &kflags, _vflags &vflags);
  bool VASP_Directory(std::ofstream &FileERROR, _aflags &aflags, _kflags &kflags);
  void VASP_BackupOriginal(_aflags aflags);
  void VASP_InitWarnings(std::unordered_map<std::string, bool> &vasp_warnings);
  void VASP_ProcessWarnings(_xvasp &xvasp, _aflags &aflags, _kflags &kflags, size_t &offset, std::unordered_map<std::string, bool> &vasp_warnings, aurostd::xoption &xmessage, aurostd::xoption &xwarning, std::ofstream &FileMESSAGE); // CO20210315
  void VASP_ProcessWarnings(
      _xvasp &xvasp, _aflags &aflags, _kflags &kflags, size_t &offset, std::unordered_map<std::string, bool> &vasp_warnings, aurostd::xoption &xmessage, aurostd::xoption &xwarning, aurostd::xoption &xmonitor, std::ofstream &FileMESSAGE); // CO20210315
  bool VASP_Error2Fix(const std::string &error, _xvasp &xvasp, aurostd::xoption &xwarning, aurostd::xoption &xfixed, _aflags &aflags, _kflags &kflags, _vflags &vflags, std::ofstream &FileMESSAGE); // CO20210315
  bool VASP_Error2Fix(const std::string &error, const std::string &mode, _xvasp &xvasp, aurostd::xoption &xwarning, aurostd::xoption &xfixed, _aflags &aflags, _kflags &kflags, _vflags &vflags, std::ofstream &FileMESSAGE); // CO20210315
  bool VASP_Error2Fix(const std::string &error, int &submode, _xvasp &xvasp, aurostd::xoption &xwarning, aurostd::xoption &xfixed, _aflags &aflags, _kflags &kflags, _vflags &vflags, std::ofstream &FileMESSAGE); // CO20210315
  bool VASP_Error2Fix(const std::string &error, const std::string &mode, int &submode, _xvasp &xvasp, aurostd::xoption &xwarning, aurostd::xoption &xfixed, _aflags &aflags, _kflags &kflags, _vflags &vflags, std::ofstream &FileMESSAGE); // CO20210315
  bool VASP_Error2Fix(const std::string &error, bool try_last_ditch_efforts, _xvasp &xvasp, aurostd::xoption &xwarning, aurostd::xoption &xfixed, _aflags &aflags, _kflags &kflags, _vflags &vflags, std::ofstream &FileMESSAGE); // CO20210315
  bool VASP_Error2Fix(const std::string &error, const std::string &mode, bool try_last_ditch_efforts, _xvasp &xvasp, aurostd::xoption &xwarning, aurostd::xoption &xfixed, _aflags &aflags, _kflags &kflags, _vflags &vflags, std::ofstream &FileMESSAGE); // CO20210315
  bool VASP_Error2Fix(const std::string &error, int &submode, bool try_last_ditch_efforts, _xvasp &xvasp, aurostd::xoption &xwarning, aurostd::xoption &xfixed, _aflags &aflags, _kflags &kflags, _vflags &vflags, std::ofstream &FileMESSAGE); // CO20210315
  bool VASP_Error2Fix(
      const std::string &error, const std::string &mode, int &submode, bool try_last_ditch_efforts, _xvasp &xvasp, aurostd::xoption &xwarning, aurostd::xoption &xfixed, _aflags &aflags, _kflags &kflags, _vflags &vflags, std::ofstream &FileMESSAGE); // CO20210315
  bool VASP_FixErrors(_xvasp &xvasp, aurostd::xoption &xmessage, aurostd::xoption &xwarning, aurostd::xoption &xfixed, _aflags &aflags, _kflags &kflags, _vflags &vflags, std::ofstream &FileMESSAGE); // CO20210315
  bool VASP_Run(_xvasp &xvasp, _aflags &aflags, _kflags &kflags, _vflags &vflags, std::ofstream &FileMESSAGE);
  bool VASP_Run(_xvasp &xvasp, _aflags &aflags, _kflags &kflags, _vflags &vflags, std::string relaxA, std::string relaxB, bool qmwrite, std::ofstream &FileMESSAGE);
  bool VASP_Run(_xvasp &xvasp, _aflags &aflags, _kflags &kflags, _vflags &vflags, std::string relaxA, bool qmwrite, std::ofstream &FileMESSAGE);
  bool VASP_RunFinished(_xvasp &xvasp, _aflags &aflags, std::ofstream &FileMESSAGE, bool = false);
  void WaitFinished(_xvasp &xvasp, _aflags &aflags, std::ofstream &FileMESSAGE, uint max_count = AUROSTD_MAX_UINT, bool = false); // CO20201220 - added max_count
  void VASP_Error(const _xvasp &xvasp, const std::string &message1 = "", const std::string &message2 = "", const std::string &message3 = "");
  void VASP_Error(const _xvasp &xvasp, std::ofstream &FileMESSAGE, const std::string &message1 = "", const std::string &message2 = "", const std::string &message3 = "");
  std::string VASP_Analyze(_xvasp &xvasp, bool qmwrite);
  void VASP_Backup(_xvasp &xvasp, bool qmwrite, const std::string &relax); // CO20210315
  void VASP_CONTCAR_Save(const _xvasp &xvasp, const std::string &relax = "breakpoint"); // CO20210315
  void VASP_CONTCAR_Save(const std::string &directory, const std::string &relax = "breakpoint"); // CO20210315
  void VASP_Recycle(const _xvasp &xvasp, const std::string &relax); // CO20210315
  void VASP_Recycle(const _xvasp &xvasp, int relax_number); // CO20210315
  void VASP_RecycleExtraFile(const _xvasp &xvasp, const std::string &xfile, const std::string &relax); // CO20210315
  void VASP_RecycleExtraFile(const _xvasp &xvasp, const std::string &xfile, int relax_number); // CO20210315
  uint VASP_getNELM(const std::string &outcar); // CO20200624
  uint VASP_getNSTEPS(const std::string &oszicar); // CO20200624
  bool VASP_OSZICARUnconverging(const std::string &dir, uint cutoff = 3);
  bool VASP_OSZICARUnconverged(const std::string &oszicar, const std::string &outcar);
  void GetStatDiel(std::string &outcar, aurostd::xvector<double> &eigr, aurostd::xvector<double> &eigi); // CAMILO
  void GetDynaDiel(std::string &outcar, aurostd::xvector<double> &eigr, aurostd::xvector<double> &eigi); // CAMILO
  std::string BIN2VASPVersion(const std::string &binfile);
  std::string OUTCAR2VASPVersion(const std::string &outcar);

} // namespace KBIN

#endif // AFLOW_KVASP_H
