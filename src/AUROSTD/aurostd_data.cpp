// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by hagen.eckert@duke.edu

#ifndef _AUROSTD_DATA_CPP_
#define _AUROSTD_DATA_CPP_

#ifndef _AUROSTD_INCBIN_H_
#define INCBIN_PREFIX afdata
#define INCBIN_SILENCE_BITCODE_WARNING
#include "aurostd_incbin.h"

#endif // _AUROSTD_INCBIN_H_

// Include files either as text (INCTXT) or binary (INCBIN) into the aflow binary
// see AUROSTD/aurostd_incbin.h
extern "C" {
  // strings used for unittests
  INCTXT(UnitTest, "DATA/unit_test.json");
  INCBIN(FileNomadLogo, "DATA/NOMAD_LOGO.png");
  INCBIN(FileAflowLogoFull, "DATA/AFLOW_LOGO_FULL.pdf");
  INCBIN(FileAflowLogoSkinny, "DATA/AFLOW_LOGO_SKINNY.pdf");
  INCTXT(CCEPython, "DATA/aflow_cce_python.py");
  INCTXT(ChullJupyterPlotter, "DATA/aflow_chull_jupyter_plotter.py");
  INCTXT(ChullJupyterRequirements, "DATA/aflow_chull_jupyter_requirements.txt");
  INCTXT(ChullJupyterJson, "DATA/aflow_chull_jupyter.json");
  INCTXT(ChullPython, "DATA/aflow_chull_python.py");
  INCTXT(SymPython, "DATA/aflow_sym_python.py");
  INCTXT(XtalfinderPython, "DATA/aflow_xtalfinder_python.py");
  INCTXT(WebappBandsJS, "DATA/aflowlib_webapp_bands.js"); // seems unused //HE20230413
  INCTXT(WebappEntryJS, "DATA/aflowlib_webapp_entry.js"); // seems unused //HE20230413
}

static const std::map<std::string, std::pair<char *, unsigned int>> aflow_data_file_list_binary = {
    {"NOMAD_LOGO.png", {(char *) &afdataFileNomadLogoData[0], afdataFileNomadLogoSize}},
    {"AFLOW_LOGO_FULL.pdf", {(char *) &afdataFileAflowLogoFullData[0], afdataFileAflowLogoFullSize}},
    {"AFLOW_LOGO_SKINNY.pdf", {(char *) &afdataFileAflowLogoSkinnyData[0], afdataFileAflowLogoSkinnySize}}
};

static const std::map<std::string, std::string> aflow_data_file_list_text = {
    {"unit_test.json", afdataUnitTestData},
    {"aflow_cce_python.py", afdataCCEPythonData},
    {"aflow_chull_jupyter_plotter.py", afdataChullJupyterPlotterData},
    {"aflow_chull_jupyter_requirements.txt", afdataChullJupyterRequirementsData},
    {"aflow_chull_jupyter.json", afdataChullJupyterJsonData},
    {"aflow_chull_python.py", afdataChullPythonData},
    {"aflow_sym_python.py", afdataSymPythonData},
    {"aflow_xtalfinder_python.py", afdataXtalfinderPythonData},
    {"aflowlib_webapp_bands.js", afdataWebappBandsJSData},
    {"aflowlib_webapp_entry.js", afdataWebappEntryJSData}
    };

namespace aurostd {

  // Helper functions to retrieve embedded data //HE20230413
  namespace EmbData {

    /// @brief return data for unittests
    /// @return content of `DATA/unit_test.json`
    JSON::object get_unit_test() {
      return JSON::loadString(afdataUnitTestData);
    }

    /// @brief get the content of an embedded file
    /// @param filename name of the file (as in the DATA folder)
    /// @return content of the embedded file
    std::string get_content(const std::string &filename) {
      if (aflow_data_file_list_text.count(filename)) {
        return aflow_data_file_list_text.at(filename);
      } else {
        string message = "File \"" + filename + " \" was not embedded";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
      }
    }

    /// @brief write a embedded file into the filesystem
    /// @param filename name of the file (as in the DATA folder)
    /// @param target_path target location
    /// @authors
    /// @mod{HE,20230413,create function}
    void save_to_file(const std::string &filename, const std::string &target_path) {
      if (aflow_data_file_list_binary.count(filename)) {
        ofstream output;
        output.open(target_path.c_str(), std::ios::out | std::ios::binary);
        output.write(aflow_data_file_list_binary.at(filename).first, aflow_data_file_list_binary.at(filename).second);
        output.close();
      } else if (aflow_data_file_list_text.count(filename)) {
        aurostd::string2file(aflow_data_file_list_text.at(filename), target_path);
      } else {
        string message = "File \"" + filename + " \" was not embedded";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
      }
    }

  } // end namespace EmbData
} // end namespace aurostd

#endif  //_AUROSTD_DATA_CPP_