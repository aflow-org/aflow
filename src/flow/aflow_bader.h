// ***************************************************************************
// *                                                                         *
// *              AFlow KEVIN RASCH - Duke University 2013                   *
// *                                                                         *
// ***************************************************************************
// aflow_bader.h
// functions written by KEVIN RASCH
// 2013: kevin.rasch@duke.edu

#ifndef _AFLOW_BADER_H_
#define _AFLOW_BADER_H_

// ***************************************************************************
#include <deque>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd_xoption.h"

// perform Bader analysis using the program from Henkelman Group at UT, Austin
namespace bader_functions {
  std::string BaderCalc(aurostd::xoption vpflow);
  bool BaderCalc(aurostd::xoption& vpflow, const std::string& bader_options, const std::string& directory, std::ostream& oss);
  bool BaderCalc(aurostd::xoption& vpflow,
                 const std::string& bader_options,
                 const std::string& prototype,
                 const std::vector<std::string>& vspecies,
                 const std::deque<int>& num_each_species,
                 const std::vector<double>& vZVAL,
                 const std::vector<double>& cutoffs,
                 const std::vector<int>& downsample_ratios,
                 const std::string& directory,
                 std::ostream& oss);
  bool BaderCalc(const std::string& bader_options,
                 const std::string& prototype,
                 const std::vector<std::string>& vspecies,
                 const std::deque<int>& num_each_species,
                 const std::vector<double>& vZVAL,
                 const std::vector<double>& cutoffs,
                 const std::vector<int>& downsample_ratios,
                 const std::string& directory,
                 std::ostream& oss);
  void FixDirectory(std::string& directory);
  bool Flags2BaderCommands(aurostd::xoption& vpflow, std::string& bader_options, std::ostream& oss);
  bool getPushCommand(const std::string& misc_option, std::string& push_command, std::ostream& oss);
  bool listORrange2vec(const std::string& misc_option, std::vector<int>& vout, std::ostream& oss);
  bool BaderExtensionFound(const std::string& FileNameIn, std::string& FileNameOut, const std::string& directory);
  bool BaderExtensionFound(const std::string& FileNameIn, const std::string& directory);
  void adjust_header(std::string& new_header, std::stringstream& FileIN_ss);
  std::string prepare_CHGCAR_4_Jmol(aurostd::xoption vpflow);
  bool get_species_string(std::string& outcar_file, std::string& species_string, const std::string& dir_to_look, const std::string& file, std::ostream& oss);
  bool prepare_CHGCAR_4_Jmol(const std::string& _chgcar_file, std::string& species_header, bool zip_file, std::ostream& oss);
  bool prepare_CHGCAR_4_Jmol(std::string& _chgcar_file, std::string& species_header, std::ostream& oss);
  bool prepare_CHGCAR_4_Jmol(std::string& _chgcar_file, std::string& species_header);
  //[CO20180220 - moved to aurostd]bool efile2tempfile(std::string _FileNameIn, std::string& FileNameOut);
} // namespace bader_functions

// CHGCAR2JVXL function
namespace pflow {
  std::string CHGCAR2JVXL(aurostd::xoption& vpflow);
  bool CHGCAR2JVXL(std::vector<std::string>& chgcar_files, const std::vector<double>& cutoffs, const std::vector<int>& downsample_ratios, const bool& cyclic, std::ostream& oss);
  bool CHGCAR2JVXL(std::vector<std::string>& chgcar_files, const std::vector<double>& cutoffs, const std::vector<int>& downsample_ratios, std::vector<std::string>& output_files, const bool& cyclic, std::ostream& oss);
  bool CHGCAR2JVXL(std::vector<std::string>& chgcar_files, const std::vector<double>& cutoffs, const std::vector<int>& downsample_ratios, std::ostream& oss);
  bool CHGCAR2JVXL(std::vector<std::string>& chgcar_files, const std::vector<double>& cutoffs, const bool& cyclic, std::ostream& oss);
  bool CHGCAR2JVXL(std::vector<std::string>& chgcar_files, const std::vector<double>& cutoffs, std::ostream& oss);
  bool CHGCAR2JVXL(std::vector<std::string>& chgcar_files, const std::vector<double>& cutoffs, const std::vector<int>& downsample_ratios, std::vector<std::string>& output_files, std::ostringstream& oss);
  bool CHGCAR2JVXL(std::vector<std::string>& chgcar_files, const std::vector<double>& cutoffs, std::vector<std::string>& output_files, std::ostringstream& oss);
  bool CHGCAR2JVXL(std::string chgcar_file, const double cutoff, const int downsample_ratio, const std::string output_file, std::ostream& oss);
  bool CHGCAR2JVXL(std::string chgcar_file, const double& cutoff, const int& downsample_ratio, std::ostream& oss);
  bool CHGCAR2JVXL(std::string chgcar_file, const double& cutoff, std::string& output_file, std::ostream& oss);
  bool CHGCAR2JVXL(std::string chgcar_file, const double& cutoff, std::ostream& oss);
  std::string CHGCAR2JVXL_get_output_filename(std::string chgcar_file, const double& cutoff, const int& downsample_ratio);
  std::string CHGCAR2JVXL_get_output_filename(std::string chgcar_file, const double& cutoff);
} // namespace pflow

#endif

// ***************************************************************************
// *                                                                         *
// *              AFlow KEVIN RASCH - Duke University 2013                   *
// *                                                                         *
// ***************************************************************************
