//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2025           *
// *                  Simon Divilov - Duke University 2025                   *
// *                                                                         *
//****************************************************************************
// Written by Simon Divilov, 2025
// simon.divilov@duke.edu
//
#include "modules/CUMULANTS/aflow_cumulants.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <deque>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xcombos.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xfile.h"
#include "AUROSTD/aurostd_xfit.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xparser_json.h"
#include "AUROSTD/aurostd_xscalar.h"
#include "AUROSTD/aurostd_xtensor.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow.h"
#include "aflow_aflowrc.h"
#include "aflow_defs.h"
#include "aflow_init.h"
#include "flow/aflow_pflow.h"
#include "flow/aflow_xclasses.h"
#include "modules/POCC/aflow_pocc.h"
#include "structure/aflow_xstructure.h"

using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::vector;

using aurostd::xmatrix;
using aurostd::xtensor;
using aurostd::xvector;

namespace cumulants {
  /// @class CumulantsCalculator
  ///
  /// @brief calculator to evaluate the cumulants and then calculate the chemical or coherent spinodal
  ///
  /// The CumulantsCalculator class is made to be called from the command line (`aflow --cumulants --usage`).
  /// Primarily, the user passes the directory containing POCC directories of the same parent lattice but at different species concentrations.
  ///
  /// Workflow of the calculator:
  ///   (1) Using umbral calculus, evaluate the cumulants from the POCC energies up to a cut-off specified by the user.
  ///   (2) Using least squares, fit a multivariate polynomial, of order specified by the user, for each cumulant, on the concentration grid specified by the user.
  ///   (3a) Analytically evaluate the second derivative of the cumulants, and construct the bordered Hessian of the free energy D^2F/Dx^2.
  ///   (3b) If POCC elastic tensors are available, also constructs the bordered Hessian of the elastic energy in the continuum limit.
  ///   (4) At each concentration, using Brent's method, find the temperature for which det(D^2F/Dx^2) = 0.
  ///   (5) Optionally, to recover the correct local concavity of the bulk free energy, apply the disorder viscosity correction (DVC)
  ///
  /// @authors
  /// @mod{SD,20250701,created class}
  ///
  /// @see
  /// @doi{TBA}
  CumulantsCalculator CumulantsCalculator::fromOptions(const aurostd::xoption& cumulants_opts) {
    CumulantsCalculator cc;
    cc.readOptions(cumulants_opts);
    return cc;
  }
} // namespace cumulants

namespace cumulants {
  /// @brief read the cumulants flags
  ///
  /// @param cumulants_opts xoptions object
  ///
  /// @authors
  /// @mod{SD,20250701,created function}
  void CumulantsCalculator::readOptions(const aurostd::xoption& cumulants_opts) {
    if (cumulants_opts.getattachedscheme("CUMULANTS::DIRECTORY").empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "CUMULANTS::DIRECTORY flag must be specified", _INPUT_MISSING_);
    }
    directory_ = std::filesystem::path(aurostd::CleanFileName(cumulants_opts.getattachedscheme("CUMULANTS::DIRECTORY")));
    if (!std::filesystem::exists(directory_)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Directory does not exist", _INPUT_ILLEGAL_);
    }
    num_cumulants_ = aurostd::string2utype<int>(cumulants_opts.getattachedscheme("CUMULANTS::NUM_CUMULANTS"));
    if (num_cumulants_ < 1) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "CUMULANTS::NUM_CUMULANTS must be a positive integer", _INPUT_ILLEGAL_);
    }
    vector<string> tokens;
    aurostd::string2tokens(cumulants_opts.getattachedscheme("CUMULANTS::POLY_ORDERS"), tokens, ",");
    poly_orders_ = aurostd::vectorstring2vectorutype<unsigned int>(tokens);
    if (poly_orders_.size() != num_cumulants_) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Mismatch between number of cumulants and amount of polynomial orders given", _INPUT_ILLEGAL_);
    }
    aurostd::string2tokens(cumulants_opts.getattachedscheme("CUMULANTS::CONC_START"), tokens, ",");
    conc_start_ = aurostd::vectorstring2vectorutype<double>(tokens);
    aurostd::string2tokens(cumulants_opts.getattachedscheme("CUMULANTS::CONC_STOP"), tokens, ",");
    conc_stop_ = aurostd::vectorstring2vectorutype<double>(tokens);
    if (conc_start_.size() != conc_stop_.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Starting and stoping concentration vectors must have the same length", _INPUT_ILLEGAL_);
    }
    conc_step_ = aurostd::string2utype<double>(cumulants_opts.getattachedscheme("CUMULANTS::CONC_STEP"));
    pocc_temp_ = aurostd::string2utype<double>(cumulants_opts.getattachedscheme("CUMULANTS::POCC_TEMP"));
    min_temp_ = aurostd::string2utype<double>(cumulants_opts.getattachedscheme("CUMULANTS::MIN_TEMP"));
    temp_range_ = aurostd::logspace(std::log10(min_temp_), std::log10(MAX_TEMP), 200);
    aurostd::string2tokens(cumulants_opts.getattachedscheme("CUMULANTS::POLY_ORDERS_ELAS"), tokens, ",");
    poly_orders_elas_ = aurostd::vectorstring2vectorutype<unsigned int>(tokens);
    if (poly_orders_elas_.size() != 0) {
      if (poly_orders_elas_.size() != 2) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Elastic polynomial orders need to be length 2", _INPUT_ILLEGAL_);
      }
      aurostd::string2tokens(cumulants_opts.getattachedscheme("CUMULANTS::KVEC"), tokens, ",");
      kvec_ = aurostd::vector2xvector(aurostd::vectorstring2vectorutype<double>(tokens));
      if (kvec_.rows != 3) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Wavevector needs to be length 3", _INPUT_ILLEGAL_);
      }
    }
    if (cumulants_opts.flag("CUMULANTS::FILTER")) {
      filter_outliers_ = true;
    }
    if (cumulants_opts.flag("CUMULANTS::DVC")) {
      dvc_ = true;
    }
    initialized_ = true;
  }

  /// @brief read the POCC data from all the directories
  ///
  /// @authors
  /// @mod{SD,20250701,created function}
  void CumulantsCalculator::readPOccData() {
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Reading POCC directories from: " + directory_.string(), aflags_, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    double energy = 0.0;
    const double volume = 0.0;
    std::deque<string> species;
    const vector<std::filesystem::directory_entry> sorted_directory_entries = [this] {
      vector<std::filesystem::directory_entry> v{std::filesystem::directory_iterator(directory_), std::filesystem::directory_iterator()};
      sort(v.begin(), v.end(), [](const auto& a, const auto& b) {return a.path().filename() < b.path().filename();});
      return v;
    }();
    for (const auto& dir_entry : sorted_directory_entries) {
      if (dir_entry.is_directory() && std::filesystem::exists(dir_entry.path() / "aflow.pocc.structures_unique.out.xz")) {
        xstructure xstr_pocc = pocc::extractPARTCAR(aurostd::file2string(dir_entry.path() / "aflow.in"));
        if (species.empty()) {
          species = xstr_pocc.species;
          alloy_name_ = aurostd::joinWDelimiter(xstr_pocc.GetElements(true), "");
        } else {
          if (species != xstr_pocc.species) {
            const string message = "POCC species are not equal across runs [expected=" + aurostd::joinWDelimiter(species, ",") + ", got=" + aurostd::joinWDelimiter(xstr_pocc.species, ",") + "]";
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
          }
        }
        const string conc = aurostd::joinWDelimiter((vector(xstr_pocc.stoich_each_type.begin(), xstr_pocc.stoich_each_type.end())), ",");
        if (xstr_pocc.stoich_each_type.size() != conc_start_.size()) {
          const string message = "POCC concentration does not have the expected length [expected=" + aurostd::utype2string(conc_start_.size()) +
                                 ", got=" + aurostd::utype2string(xstr_pocc.stoich_each_type.size()) + "]";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
        }
        concs_str_.emplace_back(conc);
        pocc_data_[conc] = aurostd::JSON::object(aurostd::JSON::object_types::DICTIONARY);
        pocc_data_[conc]["conc"] = xstr_pocc.stoich_each_type;
        if (!poly_orders_elas_.empty()) {
          vector<xmatrix<double>> pocc_rot = {aurostd::eye<double>(3)};
          const vector<pocc::POccUnit> pocc_sites = pocc::getPOccSites(xstr_pocc);
          xstructure xstr_nonpocc = pocc::createNonPOccStructure(xstr_pocc, pocc_sites);
          xstr_nonpocc.CalculateSymmetryPointGroupCrystal();
          for (size_t isym = 0; isym < xstr_nonpocc.pgroup_xtal.size(); isym++) {
            if (aurostd::RemoveWhiteSpaces(xstr_nonpocc.pgroup_xtal[isym].str_type) == "rotation") {
              pocc_rot.emplace_back(xstr_nonpocc.pgroup_xtal[isym].Uc);
            }
          }
          pocc_data_[conc]["rot"] = pocc_rot;
        }
        vector<int> pocc_idx;
        vector<unsigned int> pocc_dg;
        vector<double> pocc_energy;
        vector<double> pocc_volume;
        vector<xmatrix<double>> pocc_elas_tensor;
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Reading POCC concentration [" + conc + "] from: " + dir_entry.path().string(), aflags_, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
        for (const auto& dir_pocc_entry : std::filesystem::directory_iterator(dir_entry.path())) {
          if (dir_pocc_entry.is_directory()) {
            pocc_idx.emplace_back(aurostd::string2utype<int>(aurostd::substring2string(dir_pocc_entry.path().filename().string(), "_", "_")));
            const xQMVASP qmvasp(dir_pocc_entry.path() / DEFAULT_AFLOW_QMVASP_OUT);
            if (qmvasp.H_atom_static != AUROSTD_NAN) {
              energy = qmvasp.H_atom_static;
            } else if (qmvasp.H_atom_relax != AUROSTD_NAN) {
              energy = qmvasp.H_atom_relax;
            } else {
              energy = 0.0;
              pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "POCC directory " + dir_pocc_entry.path().string() + " is invalid", aflags_, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
            }
            pocc_energy.emplace_back(energy);
            if (energy == 0.0) {
              pocc_dg.emplace_back(0);
            } else {
                pocc_dg.emplace_back(static_cast<unsigned int>(pocc::getDGFromXStructureTitle(qmvasp.xstr_final.title)));
            }
            if (!poly_orders_elas_.empty()) {
              pocc_volume.emplace_back(qmvasp.xstr_final.GetVolume() / static_cast<double>(qmvasp.xstr_final.atoms.size()));
              xmatrix<double> elas_tensor(6, 6);
              const aurostd::JSON::object elas_data = aurostd::JSON::loadFile(dir_pocc_entry.path() / "AEL_elastic_tensor.json.xz")["elastic_stiffness_tensor"];
              for (const auto& elas_entry : elas_data) {
                const std::pair<int, int> idx = std::make_pair(elas_entry.first[2] - '0', elas_entry.first[3] - '0');
                elas_tensor(idx.first, idx.second) = static_cast<double>(elas_entry.second);
              }
              pocc_elas_tensor.emplace_back(elas_tensor);
            }
          }
        }
        pocc_data_[conc]["dg"] = vector<unsigned int>(pocc_dg.size());
        pocc_data_[conc]["energy"] = vector<double>(pocc_energy.size());
        for (size_t i = 0; i < pocc_idx.size(); i++) {
          pocc_data_[conc]["dg"][pocc_idx[i] - 1] = pocc_dg[i];
          pocc_data_[conc]["energy"][pocc_idx[i] - 1] = pocc_energy[i];
        }
        if (!poly_orders_elas_.empty()) {
          pocc_data_[conc]["volume"] = vector<double>(pocc_volume.size());
          pocc_data_[conc]["elas_tensor"] = vector<xmatrix<double>>(pocc_elas_tensor.size());
          for (size_t i = 0; i < pocc_idx.size(); i++) {
            pocc_data_[conc]["volume"][pocc_idx[i] - 1] = pocc_volume[i];
            pocc_data_[conc]["elas_tensor"][pocc_idx[i] - 1] = pocc_elas_tensor[i];
          }
        }
      }
    }
  }

  /// @brief filter the POCC data to remove outliers
  ///
  /// @authors
  /// @mod{SD,20250701,created function}
  void CumulantsCalculator::filterPOccData() {
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Filtering POCC data outliers", aflags_, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    for (const string& conc : concs_str_) {
      vector<unsigned int> pocc_dg = pocc_data_[conc]["dg"];
      vector<double> pocc_energy = pocc_data_[conc]["energy"];
      vector<double> pocc_volume;
      vector<xmatrix<double>> pocc_elas_tensor;
      if (!poly_orders_elas_.empty()) {
        pocc_volume = pocc_data_[conc]["volume"];
        pocc_elas_tensor = pocc_data_[conc]["elas_tensor"];
      }
      vector<double> pocc_energy_flatten;
      pocc_energy_flatten.reserve(std::accumulate(pocc_dg.begin(), pocc_dg.end(), 0U));
      for (size_t i = 0; i < pocc_dg.size(); i++) {
        pocc_energy_flatten.insert(pocc_energy_flatten.end(), pocc_dg[i], pocc_energy[i]);
      }
      const xvector<double> xenergy = aurostd::vector2xvector(pocc_energy_flatten);
      double q1;
      double med;
      double q3;
      aurostd::getQuartiles(xenergy, q1, med, q3);
      const double mad = aurostd::getMAD(xenergy, med);
      double zscore = 0.0;
      vector<unsigned int> pocc_dg_filter;
      pocc_dg_filter.reserve(pocc_energy.size());
      vector<double> pocc_energy_filter;
      pocc_energy_filter.reserve(pocc_energy.size());
      vector<double> pocc_volume_filter;
      pocc_volume_filter.reserve(pocc_energy.size());
      vector<xmatrix<double>> pocc_elas_tensor_filter;
      pocc_elas_tensor_filter.reserve(pocc_energy.size());
      for (size_t i = 0; i < pocc_energy.size(); i++) {
        zscore = 0.6745 * std::abs(pocc_energy[i] - med) / mad;
        if (pocc_dg[i] == 0 || zscore >= 3.5) {
          continue;
        }
        pocc_dg_filter.emplace_back(pocc_dg[i]);
        pocc_energy_filter.emplace_back(pocc_energy[i]);
        if (!poly_orders_elas_.empty()) {
          pocc_volume_filter.emplace_back(pocc_volume[i]);
          pocc_elas_tensor_filter.emplace_back(pocc_elas_tensor[i]);
        }
      }
      pocc_data_[conc]["dg"] = pocc_dg_filter;
      pocc_data_[conc]["energy"] = pocc_energy_filter;
      pocc_data_[conc]["volume"] = pocc_volume_filter;
      pocc_data_[conc]["elas_tensor"] = pocc_elas_tensor_filter;
    }
  }

  /// @brief calculate the moments and cumulants
  ///
  /// @authors
  /// @mod{SD,20250701,created function}
  void CumulantsCalculator::calculateMomentsCumulants() {
    const std::function<double(const int, const int)> fact_func = [](const int i, const int j) { return (j != 0) ? aurostd::factorial(i + 1 - j) : aurostd::factorial(i); };
    for (const string& conc : concs_str_) {
      const xvector<double> pocc_dg = pocc_data_[conc]["dg"];
      const xvector<double> pocc_energy = pocc_data_[conc]["energy"];
      vector<double> pocc_moment(num_cumulants_ + 1);
      vector<double> pocc_cumulant(num_cumulants_ + 1);
      const double pocc_dg_total = aurostd::sum(pocc_dg);
      for (size_t ic = 0; ic < num_cumulants_ + 1; ic++) {
        pocc_moment[ic] = aurostd::sum(aurostd::elementwise_product(pocc_dg, aurostd::pow(pocc_energy, static_cast<double>(ic)))) / pocc_dg_total;
      }
      for (size_t ic = 1; ic < num_cumulants_ + 1; ic++) {
        xmatrix<double> cumulant_mat(ic, ic);
        for (int i = cumulant_mat.lrows; i <= cumulant_mat.urows; i++) {
          for (int j = cumulant_mat.lcols; j <= cumulant_mat.ucols; j++) {
            if (i - cumulant_mat.lrows + 1 - j + cumulant_mat.lcols < 0) {
              continue;
            }
            cumulant_mat(i, j) = pocc_moment[i - cumulant_mat.lrows + 1 - j + cumulant_mat.lcols] / fact_func(i - cumulant_mat.lrows, j - cumulant_mat.lcols);
          }
        }
        pocc_cumulant[ic] = std::pow(-1.0, ic - 1) * aurostd::factorial(ic - 1) * aurostd::det(cumulant_mat);
      }
      pocc_data_[conc]["moment"] = pocc_moment;
      pocc_data_[conc]["cumulant"] = pocc_cumulant;
    }
  }

  /// @brief calculate the concentration grid
  ///
  /// @authors
  /// @mod{SD,20250701,created function}
  void CumulantsCalculator::calculateConcentrationGrid() {
    int num_pnts = 0;
    int tmp = 0;
    for (size_t i = 0; i < conc_start_.size(); i++) {
      tmp = static_cast<int>(std::abs((conc_start_[i] - conc_stop_[i]) / conc_step_)) + 1;
      num_pnts = std::max(tmp, num_pnts);
    }
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Calculating concentration grid [npts=" + aurostd::utype2string(num_pnts) + "]", aflags_, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    vector<vector<double>> conc_grid(conc_start_.size());
    for (size_t i = 0; i < conc_start_.size(); i++) {
      conc_grid[i] = aurostd::xvector2vector(aurostd::linspace(conc_start_[i], conc_stop_[i], num_pnts));
    }
    conc_grid_.reserve(conc_start_.size());
    for (size_t i = 0; i < conc_grid[0].size(); i++) {
      aurostd::xvector<double> xv(conc_grid.size());
      for (size_t j = 0; j < conc_grid.size(); j++) {
        xv[j + xv.lrows] = conc_grid[j][i];
      }
      conc_grid_.emplace_back(xv);
    }
  }

  /// @brief calculate the POCC probabilities
  ///
  /// @authors
  /// @mod{SD,20250701,created function}
  void CumulantsCalculator::calculateProbabilities() {
    for (const string& conc : concs_str_) {
      const xvector<long double> pocc_dg = pocc_data_[conc]["dg"];
      const xvector<long double> pocc_energy = pocc_data_[conc]["energy"];
      xvector<long double> pocc_prob = aurostd::elementwise_product(pocc_dg, aurostd::exp(-1.0 * pocc_energy / (KBOLTZEV * pocc_temp_)));
      pocc_prob = pocc_prob / aurostd::sum(pocc_prob);
      pocc_data_[conc]["prob"] = pocc_prob;
    }
  }

  /// @brief calculate the POCC averaged quantities
  ///
  /// @authors
  /// @mod{SD,20250701,created function}
  void CumulantsCalculator::calculatePOccAveraged() {
    for (const string& conc : concs_str_) {
      xvector<double> pocc_prob = pocc_data_[conc]["prob"];
      const xvector<double> pocc_energy = pocc_data_[conc]["energy"];
      pocc_data_[conc]["energy_avg"] = aurostd::sum(aurostd::elementwise_product(pocc_prob, pocc_energy));
      if (!poly_orders_elas_.empty()) {
        const xvector<double> pocc_volume = pocc_data_[conc]["volume"];
        pocc_data_[conc]["volume_avg"] = aurostd::sum(aurostd::elementwise_product(pocc_prob, pocc_volume));
        const vector<xmatrix<double>> pocc_rot = pocc_data_[conc]["rot"];
        vector<xmatrix<double>> pocc_elas_tensor = pocc_data_[conc]["elas_tensor"];
        xmatrix<double> pocc_elas_tensor_avg(6, 6);
        for (size_t ipocc = 0; ipocc < pocc_elas_tensor.size(); ipocc++) {
          const xtensor<double> elas_tensor = aurostd::convert_voigt_to_tensor(pocc_elas_tensor[ipocc]);
          xtensor<double> elas_tensor_rot(vector<int>{3, 3, 3, 3});
          for (const xmatrix<double>& rot : pocc_rot) {
            for (int i = 1; i <= 3; i++)
            for (int j = 1; j <= 3; j++)
            for (int k = 1; k <= 3; k++)
            for (int l = 1; l <= 3; l++) {
              for (int a = 1; a <= 3; a++)
              for (int b = 1; b <= 3; b++)
              for (int c = 1; c <= 3; c++)
              for (int d = 1; d <= 3; d++) {
                elas_tensor_rot[i][j][k][l] += rot[i][a] * rot[j][b] * rot[k][c] * rot[l][d] * elas_tensor[a][b][c][d];
              }
            }
          }
          elas_tensor_rot /= static_cast<double>(pocc_rot.size());
          pocc_elas_tensor_avg += pocc_prob[ipocc + pocc_prob.lrows] * aurostd::convert_tensor_to_voigt(elas_tensor_rot);
        }
        pocc_data_[conc]["elas_tensor_avg"] = pocc_elas_tensor_avg;
      }
    }
  }

  /// @brief generate data for the multivariante polynomial fitting
  ///
  /// @authors
  /// @mod{SD,20250701,created function}
  void CumulantsCalculator::generatePolyFitInputs() {
    pf_inputs_concs_.reserve(concs_str_.size());
    for (const string& conc : concs_str_) {
      const xvector<double> pocc_conc = pocc_data_[conc]["conc"];
      pf_inputs_concs_.emplace_back(pocc_conc);
    }
    xvector<double> energy_avg(concs_str_.size());
    pf_outputs_energy_avg_.reserve(1);
    for (size_t i = 0; i < concs_str_.size(); i++) {
      energy_avg[i + energy_avg.lrows] = aurostd::string2utype<double>(pocc_data_[concs_str_[i]]["energy_avg"].toString());
    }
    pf_outputs_energy_avg_.emplace_back(energy_avg);
    pf_outputs_cumulants_.reserve(num_cumulants_);
    for (size_t ic = 1; ic < num_cumulants_ + 1; ic++) {
      xvector<double> xv(concs_str_.size());
      for (size_t i = 0; i < concs_str_.size(); i++) {
        xv[i + xv.lrows] = aurostd::string2utype<double>(pocc_data_[concs_str_[i]]["cumulant"][ic].toString());
      }
      pf_outputs_cumulants_.emplace_back(xv);
    }
    if (!poly_orders_elas_.empty()) {
      xvector<double> volume_avg(concs_str_.size());
      pf_outputs_volume_avg_.reserve(1);
      for (size_t i = 0; i < concs_str_.size(); i++) {
        volume_avg[i + volume_avg.lrows] = aurostd::string2utype<double>(pocc_data_[concs_str_[i]]["volume_avg"].toString());
      }
      pf_outputs_volume_avg_.emplace_back(volume_avg);
      pf_outputs_elas_tensor_avg_.reserve(36);
      for (int irow = 1; irow <= 36; irow++) {
        xvector<double> xv(concs_str_.size());
        for (size_t i = 0; i < concs_str_.size(); i++) {
          xmatrix<double> flat = pocc_data_[concs_str_[i]]["elas_tensor_avg"];
          flat.reshape(36);
          xv[i + xv.lrows] = flat[irow][1];
        }
        pf_outputs_elas_tensor_avg_.emplace_back(xv);
      }
    }
  }

  /// @brief calculate the multivariante polynomial fitting of the data
  ///
  /// @param inputs to the polynomial fit
  /// @param outputs of the polynomial fit
  /// @param poly_orders the order of the polynomial for each input vector
  /// @param poly_coeffs coefficients of the polynomial fit for each input vector
  /// @param poly_exps exponents of the polynomial fit for each input vector
  ///
  /// @authors
  /// @mod{SD,20250701,created function}
  void CumulantsCalculator::calculatePolyFit(const vector<xvector<double>>& inputs,
                                             const vector<xvector<double>>& outputs,
                                             const vector<unsigned int>& poly_orders,
                                             vector<xvector<double>>& poly_coeffs,
                                             vector<vector<xvector<double>>>& poly_exps) {
    poly_coeffs.resize(poly_orders.size());
    poly_exps.resize(poly_orders.size());
    aurostd::xcombos xc;
    for (size_t ipo = 0; ipo < poly_orders.size(); ipo++) {
      xc = aurostd::xcombos(poly_orders[ipo] + 1, inputs[0].rows, 'E', true);
      xmatrix<double> A(inputs.size(), static_cast<int>(std::pow(poly_orders[ipo] + 1, inputs[0].rows)));
      int icol = 0;
      vector<xvector<double>> exps_all;
      while (xc.increment()) {
        icol++;
        const vector<int>& combo = xc.getCombo();
        xvector<double> exps(combo.size());
        for (size_t i = 0; i < combo.size(); i++) {
          exps[i + exps.lrows] = static_cast<double>(combo[i]);
        }
        exps_all.emplace_back(exps);
        xvector<double> col(inputs.size());
        for (size_t i = 0; i < inputs.size(); i++) {
          col[i + col.lrows] = aurostd::elements_product(aurostd::pow(inputs[i], exps));
        }
        A.setcol(col, icol);
      }
      poly_coeffs[ipo] = aurostd::LinearLeastSquares(A, outputs[ipo], true);
      poly_exps[ipo] = exps_all;
    }
  }

  /// @brief calculate the values of a multivariante polynomial fit
  ///
  /// @param poly_coeffs coefficients of the polynomial fit for each input vector
  /// @param poly_exps exponents of the polynomial fit for each input vector
  /// @param inputs to the polynomial fit
  ///
  /// @return values of the polynomial fit for each input vector
  ///
  /// @authors
  /// @mod{SD,20250704,created function}
  vector<xvector<double>> CumulantsCalculator::calculatePolyFitValues(const vector<xvector<double>>& poly_coeffs, const vector<vector<xvector<double>>>& poly_exps, const vector<xvector<double>>& inputs) {
    vector<xvector<double>> poly_values(poly_coeffs.size());
    for (size_t ipo = 0; ipo < poly_coeffs.size(); ipo++) {
      xvector<double> vals(inputs.size());
      for (size_t i = 0; i < inputs.size(); i++) {
        xvector<double> inp_exps(poly_exps[ipo].size());
        for (size_t j = 0; j < poly_exps[ipo].size(); j++) {
          inp_exps[j + inp_exps.lrows] = aurostd::elements_product(aurostd::pow(inputs[i], poly_exps[ipo][j]));
        }
        vals[i + vals.lrows] = aurostd::scalar_product(poly_coeffs[ipo], inp_exps);
      }
      poly_values[ipo] = vals;
    }
    return poly_values;
  }

  /// @brief calculate the entropy stability matrix at a particular concentration
  ///
  /// @param conc concentration to evaluate the stability matrix
  /// @param stab_mat stability matrix
  ///
  /// @authors
  /// @mod{SD,20250704,created function}
  ///
  /// @note defined with the factor of kB
  void CumulantsCalculator::calculateStabilityMatrixEntropy(const xvector<double>& conc, xmatrix<double>& stab_mat) {
    stab_mat = xmatrix<double>(conc.rows - 1, conc.rows - 1);
    for (int i = stab_mat.lrows; i <= stab_mat.urows; i++) {
      for (int j = stab_mat.lcols; j <= stab_mat.ucols; j++) {
        if (i == j) {
          stab_mat[i][j] = 1.0 / conc(i) + 1.0 / conc(conc.urows);
        } else {
          stab_mat[i][j] = 1.0 / conc(conc.urows);
        }
      }
    }
    stab_mat = -1.0 * KBOLTZEV * stab_mat;
  }

  /// @brief calculate the energy stability matrix at a particular concentration
  ///
  /// @param poly_coeffs coefficients of the polynomial fit for each input vector
  /// @param poly_exps exponents of the polynomial fit for each input vector
  /// @param conc concentration to evaluate the stability matrix
  /// @param stab_mats stability matrix for each input vector
  ///
  /// @note we are computing: D^2F/(Dx_iDx_j) \equiv D_ij = F_ij - F_in - F_nj + F_nn; i != n, j != n
  /// @note where F is the free energy and the concentration vector has size n, so D^2F/Dx^2 is a (n-1)-by-(n-1) matrix
  ///
  /// @authors
  /// @mod{SD,20250704,created function}
  void CumulantsCalculator::calculateStabilityMatrixEnergy(const vector<xvector<double>>& poly_coeffs, const vector<vector<xvector<double>>>& poly_exps, const xvector<double>& conc, vector<xmatrix<double>>& stab_mats) {
    stab_mats.clear();
    stab_mats.reserve(poly_coeffs.size());
    // Compute D_nn, two derivatives with respect to x_n
    for (size_t ipo = 0; ipo < poly_coeffs.size(); ipo++) {
      xmatrix<double> stab_mat(conc.rows - 1, conc.rows - 1);
      xvector<double> pc_nn = poly_coeffs[ipo];
      vector<xvector<double>> pe_nn = poly_exps[ipo];
      // For any term that has x_n, subtract 2 from the n-th exponent, if the exponent is less than 2, set the coefficient to zero
      for (size_t k = 0; k < pe_nn.size(); k++) {
        if (pe_nn[k][pe_nn[k].urows] < 2) {
          pc_nn[k + pc_nn.lrows] = 0.0;
          continue;
        }
        pc_nn[k + pc_nn.lrows] *= pe_nn[k][pe_nn[k].urows] * (pe_nn[k][pe_nn[k].urows] - 1.0);
        pe_nn[k][pe_nn[k].urows] -= 2.0;
      }
      xvector<double> inp_exps_nn(pe_nn.size());
      for (size_t k = 0; k < pe_nn.size(); k++) {
        inp_exps_nn[k + inp_exps_nn.lrows] = aurostd::elements_product(aurostd::pow(conc, pe_nn[k]));
      }
      // We only need to compute D_nn once
      const double D_nn = aurostd::scalar_product(pc_nn, inp_exps_nn);

      // Compute D_ij, derivative with respect to x_i and x_j
      // Loop over i and j, which go {1...n-1}, we need to compute D_ij, D_in, D_jn
      for (int i = stab_mat.lrows; i <= stab_mat.urows; i++) {
        for (int j = stab_mat.lcols; j <= stab_mat.ucols; j++) {
          stab_mat(i, j) = D_nn;
          xvector<double> pc_ij = poly_coeffs[ipo];
          vector<xvector<double>> pe_ij = poly_exps[ipo];
          if (i == j) {
            // Compute D_ii
            // For any term that has x_i, subtract 2 from the i-th exponent, if the exponent is less than 2, set the coefficient to zero
            for (size_t k = 0; k < pe_ij.size(); k++) {
              if (pe_ij[k](i) < 2) {
                pc_ij[k + pc_ij.lrows] = 0.0;
                continue;
              }
              pc_ij[k + pc_ij.lrows] *= pe_ij[k][i] * (pe_ij[k][i] - 1.0);
              pe_ij[k][i] -= 2.0;
            }
          } else {
            // Compute D_ij
            // For any term that has x_i and x_j, subtract 1 from the i-th exponent and subtract 1 from the j-th exponent, if the exponents are less than 1, set the coefficient to zero
            for (size_t k = 0; k < pe_ij.size(); k++) {
              if (pe_ij[k][i] == 0 || pe_ij[k][j] == 0) {
                pc_ij[k + pc_ij.lrows] = 0.0;
                continue;
              }
              pc_ij[k + pc_ij.lrows] *= pe_ij[k][i] * pe_ij[k][j];
              pe_ij[k][i] -= 1.0;
              pe_ij[k][j] -= 1.0;
            }
          }
          xvector<double> inp_exps_ij(pe_ij.size());
          for (size_t k = 0; k < pe_ij.size(); k++) {
            inp_exps_ij[k + inp_exps_ij.lrows] = aurostd::elements_product(aurostd::pow(conc, pe_ij[k]));
          }
          stab_mat[i][j] += aurostd::scalar_product(pc_ij, inp_exps_ij);
          // Compute D_in
          // For any term that has x_i and x_n, subtract 1 from the i-th exponent and subtract 1 from the n-th exponent, if the exponents are less than 1, set the coefficient to zero
          xvector<double> pc_in = poly_coeffs[ipo];
          vector<xvector<double>> pe_in = poly_exps[ipo];
          for (size_t k = 0; k < pe_in.size(); k++) {
            if (pe_in[k][i] == 0 || pe_in[k](pe_in[k].urows) == 0) {
              pc_in[k + pc_in.lrows] = 0.0;
              continue;
            }
            pc_in[k + pc_in.lrows] *= pe_in[k][i] * pe_in[k][pe_in[k].urows];
            pe_in[k][i] -= 1.0;
            pe_in[k][pe_in[k].urows] -= 1.0;
          }
          xvector<double> inp_exps_in(pe_in.size());
          for (size_t k = 0; k < pe_in.size(); k++) {
            inp_exps_in[k + inp_exps_in.lrows] = aurostd::elements_product(aurostd::pow(conc, pe_in[k]));
          }
          stab_mat[i][j] -= aurostd::scalar_product(pc_in, inp_exps_in);
          // Compute D_nj
          // For any term that has x_n and x_j, subtract 1 from the n-th exponent and subtract 1 from the j-th exponent, if the exponents are less than 1, set the coefficient to zero
          xvector<double> pc_nj = poly_coeffs[ipo];
          vector<xvector<double>> pe_nj = poly_exps[ipo];
          for (size_t k = 0; k < pe_nj.size(); k++) {
            if (pe_nj[k][j] == 0 || pe_nj[k][pe_nj[k].urows] == 0) {
              pc_nj[k + pc_nj.lrows] = 0.0;
              continue;
            }
            pc_nj[k + pc_nj.lrows] *= pe_nj[k][j] * pe_nj[k][pe_nj[k].urows];
            pe_nj[k][j] -= 1.0;
            pe_nj[k][pe_nj[k].urows] -= 1.0;
          }
          xvector<double> inp_exps_nj(pe_nj.size());
          for (size_t k = 0; k < pe_nj.size(); k++) {
            inp_exps_nj[k + inp_exps_nj.lrows] = aurostd::elements_product(aurostd::pow(conc, pe_nj[k]));
          }
          stab_mat[i][j] -= aurostd::scalar_product(pc_nj, inp_exps_nj);
        }
      }
      stab_mats.emplace_back(stab_mat);
    }
  }

  /// @brief calculate the elastic energy stability matrix at a particular concentration
  ///
  /// @param poly_coeffs_elas_tensor coefficients of the polynomial fit for each element of the elastic tensor (Voigt notation)
  /// @param poly_exps_elas_tensor exponents of the polynomial fit for each of the elastic tensor (Voigt notation)
  /// @param poly_coeffs_volume coefficients of the polynomial fit for the volume
  /// @param poly_exps_volume exponents of the polynomial fit for the volume
  /// @param conc concentration to evaluate the stability matrix
  /// @param stab_mats stability matrix for each input vector
  ///
  /// @note we are computing: E_ij = c_abcd * \eta_iab * \eta_jcd - (c_abcd * k_b * \eta_icd) * (c_acbd * k_c * k_d)^(-1) * (c_abcd * k_a * \eta_jcd)
  /// @note a,b,c,d \in [1,3]; i,j \in [1,n-1]
  /// @note see Eq. 12.9 in "Configurational Thermodynamics of Solid Solutions"
  ///
  /// @doi{10.1016/S0081-1947(08)60360-4}
  ///
  /// @authors
  /// @mod{SD,20250828,created function}
  void CumulantsCalculator::calculateStabilityMatrixElastic(const vector<xvector<double>>& poly_coeffs_elas_tensor, const vector<vector<xvector<double>>>& poly_exps_elas_tensor,
                                                     const vector<xvector<double>>& poly_coeffs_volume, const vector<vector<xvector<double>>>& poly_exps_volume,
                                                     const xvector<double>& conc, xmatrix<double>& stab_mat) {
    vector<xvector<double>> elas_tensor_vec = calculatePolyFitValues(poly_coeffs_elas_tensor, poly_exps_elas_tensor, vector<xvector<double>>{conc});
    xmatrix<double> elas_tensor_flat(36, 1);
    for (size_t i = 0; i < elas_tensor_vec.size(); i++) {
        elas_tensor_flat[i + elas_tensor_flat.lrows][1] = elas_tensor_vec[i][elas_tensor_vec[i].lrows];
    }
    elas_tensor_flat.reshape(6, 6);
    const xtensor<double> elas_tensor = GPa2eV * convert_voigt_to_tensor(elas_tensor_flat); // eV/Ang^3
    const double volume = calculatePolyFitValues(poly_coeffs_volume, poly_exps_volume, vector<xvector<double>>{conc})[0][1];
    xvector<double> cce_tensor_flat(stab_mat.rows); // coefficient of chemical expansion
    // Compute D_n, one derivative with respect to x_n
    xvector<double> pc_n = poly_coeffs_volume[0];
    vector<xvector<double>> pe_n = poly_exps_volume[0];
    // For any term that has x_n, subtract 1 from the n-th exponent, if the exponent is less than 1, set the coefficient to zero
    for (size_t k = 0; k < pe_n.size(); k++) {
      if (pe_n[k][pe_n[k].urows] < 1) {
        pc_n[k + pc_n.lrows] = 0.0;
        continue;
      }
      pc_n[k + pc_n.lrows] *= pe_n[k][pe_n[k].urows];
      pe_n[k][pe_n[k].urows] -= 1.0;
    }
    xvector<double> inp_exps_n(pe_n.size());
    for (size_t k = 0; k < pe_n.size(); k++) {
      inp_exps_n[k + inp_exps_n.lrows] = aurostd::elements_product(aurostd::pow(conc, pe_n[k]));
    }
    // We only need to compute D_n once
    const double D_n = aurostd::scalar_product(pc_n, inp_exps_n);
    // Loop over i, which goes {1...n-1}
    for (int i = cce_tensor_flat.lrows; i <= cce_tensor_flat.urows; i++) {
      // Compute D_i
      // For any term that has x_i, subtract 1 from the i-th exponent, if the exponents are less than 1, set the coefficient to zero
      xvector<double> pc_i = poly_coeffs_volume[0];
      vector<xvector<double>> pe_i = poly_exps_volume[0];
      for (size_t k = 0; k < pe_i.size(); k++) {
        if (pe_i[k][i] == 0) {
          pc_i[k + pc_i.lrows] = 0.0;
          continue;
        }
        pc_i[k + pc_i.lrows] *= pe_i[k][i];
        pe_i[k][i] -= 1.0;
      }
      xvector<double> inp_exps_i(pe_i.size());
      for (size_t k = 0; k < pe_i.size(); k++) {
        inp_exps_i[k + inp_exps_i.lrows] = aurostd::elements_product(aurostd::pow(conc, pe_i[k]));
      }
      cce_tensor_flat[i] = aurostd::scalar_product(pc_i, inp_exps_i) - D_n;
    }
    // Calculate CCE tensor in the isotropic limit
    const xtensor<double> cce_tensor(vector<int>{cce_tensor_flat.rows, 3, 3});
    for (int i = cce_tensor_flat.lrows; i <= cce_tensor_flat.urows; i++) {
      for (int a = 1; a <= 3; a++)
      for (int b = 1; b <= 3; b++) {
        if (a == b) {
          cce_tensor[i][a][b] = cce_tensor_flat[i] / (3.0 * volume);
        }
      }
    }
    // Calculate k-independent tensor
    xmatrix<double> k_indep_tensor(stab_mat.rows, stab_mat.rows);
    for (int i = k_indep_tensor.lrows; i <= k_indep_tensor.urows; i++)
    for (int j = k_indep_tensor.lcols; j <= k_indep_tensor.ucols; j++) {
      for (int a = 1; a <= 3; a++)
      for (int b = 1; b <= 3; b++)
      for (int c = 1; c <= 3; c++)
      for (int d = 1; d <= 3; d++) {
        k_indep_tensor[i][j] += elas_tensor[a][b][c][d] * cce_tensor[i][a][b] * cce_tensor[j][c][d];
      }
    }
    // Calculate Green's tensor
    xmatrix<double> green_tensor(3, 3);
    for (int a = 1; a <= 3; a++)
    for (int b = 1; b <= 3; b++)
    for (int c = 1; c <= 3; c++)
    for (int d = 1; d <= 3; d++) {
      green_tensor[a][b] += elas_tensor[a][c][b][d] * kvec_[c] * kvec_[d]; // note the order
    }
    green_tensor = aurostd::inverse(green_tensor);
    // Calculate the left and right k-dependent tensors
    xmatrix<double> k_dep_l_tensor(stab_mat.rows, 3);
    xmatrix<double> k_dep_r_tensor(stab_mat.rows, 3);
    for (int i = k_dep_l_tensor.lrows; i <= k_dep_l_tensor.urows; i++) {
      for (int a = 1; a <= 3; a++)
      for (int b = 1; b <= 3; b++)
      for (int c = 1; c <= 3; c++)
      for (int d = 1; d <= 3; d++) {
        k_dep_l_tensor[i][a] += elas_tensor[a][b][c][d] * kvec_[b] * cce_tensor[i][c][d];
        k_dep_r_tensor[i][b] += elas_tensor[a][b][c][d] * kvec_[a] * cce_tensor[i][c][d];
      }
    }
    // Combine the terms
    for (int i = stab_mat.lrows; i <= stab_mat.urows; i++)
    for (int j = stab_mat.lcols; j <= stab_mat.ucols; j++) {
      stab_mat[i][j] = volume * k_indep_tensor[i][j];
      for (int a = 1; a <= 3; a++)
      for (int b = 1; b <= 3; b++) {
        stab_mat[i][j] -= volume * (k_dep_l_tensor[i][a] * green_tensor[a][b] * k_dep_r_tensor[j][b]);
      }
    }
  }

  /// @brief calculate the critical temperature on the concentration grid
  ///
  /// @authors
  /// @mod{SD,20250704,created function}
  void CumulantsCalculator::calculateCritTemp() {
    crit_temps_.clear();
    crit_temps_.reserve(conc_grid_.size());
    for (const xvector<double>& conc : conc_grid_) {
      xmatrix<double> stab_mat_entropy;
      vector<xmatrix<double>> stab_mats_cumulants;
      vector<xmatrix<double>> stab_mats_energy_avg;
      calculateStabilityMatrixEntropy(conc, stab_mat_entropy);
      calculateStabilityMatrixEnergy(poly_coeffs_cumulants_, poly_exps_cumulants_, conc, stab_mats_cumulants);
      const std::function<xmatrix<double>(double)> stab_mat_energy = [stab_mats_cumulants](double temp) {
        xmatrix<double> stab_mat(stab_mats_cumulants[0].rows, stab_mats_cumulants[0].cols);
        for (size_t ipo = 0; ipo < stab_mats_cumulants.size(); ipo++) {
          stab_mat += stab_mats_cumulants[ipo] * std::pow(-1.0 / (KBOLTZEV * temp), static_cast<double>(ipo)) * std::pow(aurostd::factorial(ipo + 1), -1.0);
        }
        return stab_mat;
      };
      xmatrix<double> stab_mat_elastic(stab_mat_entropy.rows, stab_mat_entropy.cols);
      if (!poly_orders_elas_.empty()) {
        calculateStabilityMatrixElastic(poly_coeffs_elas_tensor_avg_, poly_exps_elas_tensor_avg_, poly_coeffs_volume_avg_, poly_exps_volume_avg_, conc, stab_mat_elastic);
      }
      calculateStabilityMatrixEnergy(poly_coeffs_energy_avg_, poly_exps_energy_avg_, conc, stab_mats_energy_avg);
      const xmatrix<double> stab_mat_energy_corr = dvc_factor_ * stab_mats_energy_avg[0];
      const std::function<double(double)> stab_mat_det = [stab_mat_energy, stab_mat_elastic, stab_mat_energy_corr, stab_mat_entropy](double temp) {
        const double det = aurostd::det(stab_mat_energy(temp) + stab_mat_elastic - stab_mat_energy_corr - temp * stab_mat_entropy);
        return det;
      };

      double crit_temp;
      try {
        aurostd::findZeroBrent(min_temp_, MAX_TEMP, stab_mat_det, crit_temp);
      } catch (aurostd::xerror& excpt) {
        crit_temp = 0.0;
      }
      crit_temps_.emplace_back(crit_temp);
    }
  }

  /// @brief returns the bolean if the critical temperature as a function of the concentration is all zeros or there are zeros in the interior
  ///
  /// @return true if the critical temperature is zero everywhere or in the interior
  ///
  /// @authors
  /// @mod{SD,20250704,created function}
  /// @mod{SD,20251218,updated to check interior points}
  bool CumulantsCalculator::isZeros() {
    const size_t zeros = std::count(crit_temps_.begin(), crit_temps_.end(), 0.0); // no need to worry about floating-point
    if (zeros == crit_temps_.size()) {
      return true;
    }
    if (crit_temps_.size() < 3) {
      return false;
    }

    size_t first_non = crit_temps_.size();
    for (size_t i = 0; i < crit_temps_.size(); i++) {
      if (crit_temps_[i] != 0.0) {
          first_non = i; 
          break; 
      }
    }
    size_t last_non = crit_temps_.size();
    for (size_t i = crit_temps_.size(); i-- > 0; ) {
        if (crit_temps_[i] != 0.0) {
          last_non = i;
          break;
        }
    }

    // interior zeros cannot exist
    if (last_non == crit_temps_.size() || first_non >= last_non) {
      return false;
    }

    // check for interior zeros
    for (size_t i = first_non + 1; i < last_non; i++) {
        if (crit_temps_[i] == 0.0) {
          return true;
        }
    }
    return false;
  }

  /// @brief calculate, self-consistently, the critical temperature on the concentration grid
  ///
  /// @authors
  /// @mod{SD,20250704,created function}
  void CumulantsCalculator::calculateCritTempSC() {
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, string("Calculating critical temperature on the concentration grid ") + (dvc_ ? "WITH" : "WITHOUT") + " DVC", aflags_, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    calculateCritTemp();
    if (dvc_) {
      while (dvc_factor_ < 1.0 && !isZeros()) {
        dvc_factor_ += 1e-2;
        calculateCritTemp();
      }
      if (dvc_factor_ > 1) {
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "DVC factor at maximum value", aflags_, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        dvc_factor_ = 1;
      }
      dvc_factor_ /= 2.0;
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "DVC factor converged to " + aurostd::utype2string(dvc_factor_, 3), aflags_, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      calculateCritTemp();
    }
  }

  /// @brief calculate the spinodal using the cumulant expansion
  ///
  /// @authors
  /// @mod{SD,20250701,created function}
  void CumulantsCalculator::calculateSpinodal() {
    readPOccData();
    if (filter_outliers_) {
      filterPOccData();
    }
    calculateMomentsCumulants();
    calculateProbabilities();
    calculatePOccAveraged();
    calculateConcentrationGrid();
    generatePolyFitInputs();
    calculatePolyFit(pf_inputs_concs_, pf_outputs_cumulants_, poly_orders_, poly_coeffs_cumulants_, poly_exps_cumulants_);
    calculatePolyFit(pf_inputs_concs_, pf_outputs_energy_avg_, vector{poly_orders_[0]}, poly_coeffs_energy_avg_, poly_exps_energy_avg_);
    if (!poly_orders_elas_.empty()) {
      calculatePolyFit(pf_inputs_concs_, pf_outputs_volume_avg_, vector{poly_orders_elas_[0]}, poly_coeffs_volume_avg_, poly_exps_volume_avg_);
      calculatePolyFit(pf_inputs_concs_, pf_outputs_elas_tensor_avg_, vector<unsigned int>(36, poly_orders_elas_[1]), poly_coeffs_elas_tensor_avg_, poly_exps_elas_tensor_avg_);
    }
    calculateCritTempSC();
  }

  /// @brief print the output of the calculations
  ///
  /// @authors
  /// @mod{SD,20250706,created function}
  void CumulantsCalculator::printOutput() const {
    stringstream output;
    output << AFLOWIN_SEPARATION_LINE << endl;
    for (size_t i = 0; i < conc_grid_.size(); i++) {
      output << std::setw(3) << i + 1 << " [" << conc_grid_[i] << "] " << std::fixed << std::setprecision(3) << std::setw(8) << crit_temps_[i] << endl;
    }
    output << AFLOWIN_SEPARATION_LINE << endl;
    cout << output.str() << endl;
  }

  /// @brief write the output of the calculations
  ///
  /// @authors
  /// @mod{SD,20250706,created function}
  void CumulantsCalculator::writeOutput() const {
    const string filename = "aflow_cumulants_" + alloy_name_ + ".json";
    const std::filesystem::path filepath = std::filesystem::absolute(std::filesystem::current_path()) / filename;
    const aurostd::JSON::object jo = getOutput();
    jo.saveFile(filepath);
  }

  /// @brief get the output of the calculations
  ///
  /// @authors
  /// @mod{SD,20251010,created function}
  aurostd::JSON::object CumulantsCalculator::getOutput() const {
    aurostd::JSON::object jo(aurostd::JSON::object_types::DICTIONARY);
    jo["alloy_name"] = alloy_name_;
    jo["pocc_data"] = pocc_data_;
    jo["filter_outliers"] = filter_outliers_;
    jo["dvc"] = dvc_;
    jo["num_cumulants"] = num_cumulants_;
    jo["poly_orders"] = poly_orders_;
    jo["conc_grid"] = conc_grid_;
    jo["pocc_temp"] = pocc_temp_;
    jo["dvc_factor"] = dvc_factor_;
    jo["crit_temps"] = crit_temps_;
    jo["poly_orders_elas"] = poly_orders_elas_;
    jo["kvec"] = kvec_;
    return jo;
  }
} // namespace cumulants

namespace cumulants {
  /// @brief Given a set of POCC calculations at different concentrations computes the spinodal
  ///
  /// @param cumulants_opts command line options
  ///
  /// @authors
  /// @mod{SD,20250701,created function}
  void cumulantsMain(const aurostd::xoption& cumulants_opts) {
    if (cumulants_opts.flag("CUMULANTS::USAGE")) {
      displayUsage();
      return;
    }
    cumulants::CumulantsCalculator cc = cumulants::CumulantsCalculator::fromOptions(cumulants_opts);
    cc.calculateSpinodal();
    if (cumulants_opts.flag("CUMULANTS::PRINT")) {
      cc.printOutput();
    }
    if (cumulants_opts.flag("CUMULANTS::WRITE")) {
      cc.writeOutput();
    }
  }

  /// @brief displays the usage commands and options
  ///
  /// @authors
  /// @mod{SD,20250701,created function}
  void displayUsage() {
    const vector<string> usage_options = {
        "aflow --cumulants --directory=dir [cumulants_options]",
        " ",
        "The directory structure is assumed to be the following:",
        "dir/",
        "├── POCC_at_x1/",
        "│   ├── ARUN_1/",
        "│   ├── ARUN_2/",
        "│   ├── ARUN_3/",
        "│   :",
        "│   └── ARUN_N1",
        "├── POCC_at_x2/",
        "│   ├── ARUN_1/",
        "│   ├── ARUN_2/",
        "│   ├── ARUN_3/",
        "│   :",
        "│   └── ARUN_N2/",
        ":",
        ":",
        "├── POCC_at_xG/",
        "│   ├── ARUN_1/",
        "│   ├── ARUN_2/",
        "│   ├── ARUN_3/",
        "│   :",
        "│   └── ARUN_NG/",
        " ",
        "cumulants_options:",
        " ",
        "GENERAL OPTIONS:",
        "--usage",
        "--print",
        "--write",
        "--filter_outliers|--filter",
        "--dvc",
        " ",
        "CHEMICAL SPINODAL OPTIONS:",
        "--num_cumulants=|--nc=2",
        "--poly_orders=|--orders=|--po=1,1",
        "--conc_start=0.9,0.1",
        "--conc_stop=0.1,0.9",
        "--conc_step=0.025",
        "--pocc_temp=2000",
        "--min_temp=200",
        " ",
        "COHERENT SPINODAL OPTIONS:",
        "--kvec=1.0,0.0,0.0",
        "--poly_orders_el=|--orders_el=|--po_el=2,4",
    };
    init::MessageOption("--usage", "CUMULANTS()", usage_options);
  }
} // namespace cumulants
