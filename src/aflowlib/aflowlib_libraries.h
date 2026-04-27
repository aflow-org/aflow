
#ifndef AFLOWLIB_LIBRARIES_H
#define AFLOWLIB_LIBRARIES_H

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflowlib/aflowlib_web_interface.h"
#include "flow/aflow_xclasses.h"

namespace aflowlib {
// ***************************************************************************
// XPLUG/XUPDATE STUFF
  class ARun : public xStream {  // CO20210302
  private:
    ARun(const std::string& dir, std::ostream& oss = std::cout) : xStream(oss), m_initialized(!dir.empty()), m_dir(dir) {}
    ARun(const std::string& dir, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout) : xStream(FileMESSAGE, oss), m_initialized(!dir.empty()), m_dir(dir) {}
    void clear();

      // attributes
    bool m_initialized = false;
    std::string m_dir;
    std::string m_aflowin;
    std::string m_LOCK;
    bool m_completed = false;

  public:
    void free();
  };

// ***************************************************************************
// to create/analyze the tokens and load the stuff up for qhull
  bool TokenPresentAFLOWLIB(const std::string& line, const std::string& query);

// ***************************************************************************
// LIBS to RAWS of each entry
  uint GetSpeciesDirectory(const std::string& directory, std::vector<std::string>& vspecies);

  void XFIX_LIBRARY_ALL(const std::string& LIBRARY_IN, std::vector<std::string>);
  void XPLUG(const aurostd::xoption& xplug_opts);
  void XPLUG_Check(const std::vector<std::filesystem::path>& vdir_all, std::vector<std::filesystem::path>& vdir_zip, const bool clean, const double stale_threshold_days);
  void XPLUG_Zip(const std::vector<std::filesystem::path>& vdir_zip, const int zip_size_gb, const std::filesystem::path& relative_path, const std::string& prefix);
  std::string LIB2RAW_CheckProjectFromDirectory(const std::string& directory);
  bool LIB2RAW_ALL(const std::string& options, bool overwrite);

  void LIB2RAW_FileNeeded(const std::string& directory_LIB, const std::string& fileLIB, const std::string& directory_RAW, const std::string& fileRAW, std::vector<std::string>& vfiles, const std::string& MESSAGE);
  void CleanDirectoryLIB(std::string& directory);  // CO20200624
  void setAURL(aflowlib::_aflowlib_entry& aflowlib_data, const std::string& directory_LIB, bool LOCAL = false);

  void LIB2RAW(const std::string& options, bool flag_FORCE, bool LOCAL = false);
  bool AddFileNameBeforeExtension(const std::string& _file, const std::string& addendum, std::string& out_file); // CO20171025
  bool LIB2RAW_Calculate_FormationEnthalpy(aflowlib::_aflowlib_entry& data, const xstructure& xstr, const std::string& MESSAGE); // CO20200731
  void LIB2RAW_Loop_Thermodynamics(const std::string& directory_LIB, const std::string& directory_RAW, std::vector<std::string>& vfiles, aflowlib::_aflowlib_entry& data, const std::string& MESSAGE, bool LOCAL = false);
  void LIB2RAW_Loop_Static(const std::string& directory_LIB, const std::string& directory_RAW, std::vector<std::string>& vfiles, aflowlib::_aflowlib_entry& data, const std::string& MESSAGE);
  void LIB2RAW_Loop_Bands(const std::string& directory_LIB, const std::string& directory_RAW, std::vector<std::string>& vfiles, aflowlib::_aflowlib_entry& data, const std::string& MESSAGE);
  void LIB2RAW_Loop_Magnetic(const std::string& directory_LIB, const std::string& directory_RAW, std::vector<std::string>& vfiles, aflowlib::_aflowlib_entry& data, const std::string& MESSAGE);
  void LIB2RAW_Loop_Dielectric(const std::string& directory_LIB, const std::string& directory_RAW, std::vector<std::string>& vfiles, aflowlib::_aflowlib_entry& data, const std::string& MESSAGE);
  void LIB2RAW_Loop_Bader(const std::string& directory_LIB, const std::string& directory_RAW, std::vector<std::string>& vfiles, aflowlib::_aflowlib_entry& data, const std::string& MESSAGE);
  void LIB2RAW_Loop_AGL(const std::string& directory_LIB, const std::string& directory_RAW, std::vector<std::string>& vfiles, aflowlib::_aflowlib_entry& data, const std::string& MESSAGE);
  void LIB2RAW_Loop_AEL(const std::string& directory_LIB, const std::string& directory_RAW, std::vector<std::string>& vfiles, aflowlib::_aflowlib_entry& data, const std::string& MESSAGE);
  void LIB2RAW_Loop_APL(const std::string& directory_LIB, const std::string& directory_RAW, std::vector<std::string>& vfiles, aflowlib::_aflowlib_entry& data, const std::string& MESSAGE);
  void LIB2RAW_Loop_QHA(const std::string& directory_LIB, const std::string& directory_RAW, std::vector<std::string>& vfiles, aflowlib::_aflowlib_entry& data, const std::string& MESSAGE);
  void LIB2RAW_Loop_LOCK(const std::string& directory_LIB, const std::string& directory_RAW, std::vector<std::string>& vfiles, aflowlib::_aflowlib_entry& data, const std::string& MESSAGE);
  void LIB2RAW_Loop_POCC(const std::string& directory_LIB, const std::string& directory_RAW, std::vector<std::string>& vfiles, aflowlib::_aflowlib_entry& data, const std::string& MESSAGE);
  void LIB2RAW_Loop_PATCH(const std::string& directory_LIB, const std::string& directory_RAW, std::vector<std::string>& vfiles, aflowlib::_aflowlib_entry& data, const std::string& MESSAGE);
  bool LIB2LIB(const std::string& options, bool overwrite, bool LOCAL = false); // CT20181212

  bool VaspFileExist(const std::string& str_dir, const std::string& FILE);
  std::string vaspfile2stringstream(const std::string& str_dir, const std::string& FILE, std::stringstream& streamFILE);
  std::string vaspfile2stringstream(const std::string& str_dir, const std::string& FILE);

  void insertStoichStats(const std::vector<std::string> vstats, const aurostd::xvector<double>& nspecies_xv, const aurostd::xvector<double>& stoich_xv, std::vector<double>& vfeatures);  // CO20201111

} // namespace aflowlib

#endif // AFLOWLIB_LIBRARIES_H
