//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2025           *
// *                  Simon Divilov - Duke University 2025                   *
// *                                                                         *
//****************************************************************************
// Written by Simon Divilov, 2025
// simon.divilov@duke.edu
//
#ifndef _AFLOW_CUMULANTS_H_
#define _AFLOW_CUMULANTS_H_

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xparser_json.h"
#include "AUROSTD/aurostd_xvector.h"

#include "flow/aflow_xclasses.h"

namespace cumulants {
  class CumulantsCalculator : public xStream {
  public:
    CumulantsCalculator(std::ostream& oss = std::cout) : xStream(oss) {}
    CumulantsCalculator(std::ofstream& FileMESSAGE, std::ostream& oss = std::cout) : xStream(FileMESSAGE, oss) {}

    void initialize();
    static CumulantsCalculator fromOptions(const aurostd::xoption& cumulants_opts);

    void readOptions(const aurostd::xoption& cumulants_opts);
    void calculateSpinodal();
    void printOutput() const;
    void writeOutput() const;
    aurostd::JSON::object getOutput() const;

  private:
    static constexpr double MAX_TEMP = 10000; ///< upper bound temperature to search solutions (units: K)
    bool initialized_ = false; ///< CumulantsCalculator is initialized
    _aflags aflags_; ///< AFLOW flags
    std::filesystem::path directory_; ///< directory of POCC directories
    std::string alloy_name_; ///< alloy name of the POCC data
    std::vector<std::string> concs_str_; ///< concentration strings of the POCC data
    aurostd::JSON::object pocc_data_ = aurostd::JSON::object(aurostd::JSON::object_types::DICTIONARY); ///< JSON object of POCC data
    bool filter_outliers_ = false; ///< filter outlier POCC data
    bool dvc_ = false; ///< apply the disorder viscosity correction
    bool print_ = false; ///< print the critical temperatures
    int num_cumulants_ = 0; ///< number of cumulants to compute
    std::vector<unsigned int> poly_orders_; ///< order of the polynomial for each cumulant
    std::vector<double> conc_start_; ///< concentration starting vector
    std::vector<double> conc_stop_; ///< concentration stoping vector
    double conc_step_ = 0.0; ///< concentration step size
    std::vector<aurostd::xvector<double>> conc_grid_; ///< concentration grid
    double pocc_temp_ = 0.0; ///< POCC temperature (units: K)
    double min_temp_ = 0.0; ///< minimum temperature to consider when searching for solution (units: K)
    aurostd::xvector<double> temp_range_; ///< range at which to evaluate the cumulant expansion (units: K)
    std::vector<aurostd::xvector<double>> pf_inputs_concs_; ///< concentrations for the polynomial regression
    std::vector<aurostd::xvector<double>> pf_outputs_cumulants_; ///< cumulants for the polynomial regression
    std::vector<aurostd::xvector<double>> pf_outputs_energy_avg_; ///< averaged energy for the polynomial regression
    std::vector<aurostd::xvector<double>> pf_outputs_volume_avg_; ///< averaged vol for the polynomial regression
    std::vector<aurostd::xvector<double>> poly_coeffs_cumulants_; ///< polynomial coefficients for cumulants
    std::vector<aurostd::xvector<double>> poly_coeffs_energy_avg_; ///< polynomial coefficients for averaged energy per atom
    std::vector<aurostd::xvector<double>> poly_coeffs_volume_avg_; ///< polynomial coefficients for averaged volume per atom
    std::vector<std::vector<aurostd::xvector<double>>> poly_exps_cumulants_; ///< polynomial exponents for cumulants
    std::vector<std::vector<aurostd::xvector<double>>> poly_exps_energy_avg_; ///< polynomial exponents for averaged energy per atom
    std::vector<std::vector<aurostd::xvector<double>>> poly_exps_volume_avg_; ///< polynomial exponents for averaged volume per atom
    double dvc_factor_ = 0.0; ///< value of the disorder viscosity correction factor
    std::vector<double> crit_temps_; ///< critical temperatures on the concentration grid
    aurostd::xvector<double> kvec_; ///< directional cosine wavevector
    std::vector<unsigned int> poly_orders_elas_; ///< order of the polynomial for the volume and elastic tensor, respectively
    std::vector<aurostd::xvector<double>> pf_outputs_elas_tensor_avg_; ///< averaged elastic tensor (Voigt notation) for the polynomial regression
    std::vector<aurostd::xvector<double>> poly_coeffs_elas_tensor_avg_; ///< polynomial coefficients for averaged elastic tensor
    std::vector<std::vector<aurostd::xvector<double>>> poly_exps_elas_tensor_avg_; ///< polynomial exponents for averaged elastic tensor

    void readPOccData();
    void filterPOccData();
    void calculateMomentsCumulants();
    void calculateConcentrationGrid();
    void calculateProbabilities();
    void calculatePOccAveraged();
    void generatePolyFitInputs();
    void calculatePolyFit(const std::vector<aurostd::xvector<double>>& inputs,
                          const std::vector<aurostd::xvector<double>>& outputs,
                          const std::vector<unsigned int>& poly_orders,
                          std::vector<aurostd::xvector<double>>& poly_coeffs,
                          std::vector<std::vector<aurostd::xvector<double>>>& poly_exps);
    std::vector<aurostd::xvector<double>> calculatePolyFitValues(const std::vector<aurostd::xvector<double>>& poly_coeffs,
                                                                 const std::vector<std::vector<aurostd::xvector<double>>>& poly_exps,
                                                                 const std::vector<aurostd::xvector<double>>& inputs);
    void calculateStabilityMatrixEntropy(const aurostd::xvector<double>& conc, aurostd::xmatrix<double>& stab_mat);
    void calculateStabilityMatrixEnergy(const std::vector<aurostd::xvector<double>>& poly_coeffs,
                                        const std::vector<std::vector<aurostd::xvector<double>>>& poly_exps,
                                        const aurostd::xvector<double>& conc,
                                        std::vector<aurostd::xmatrix<double>>& stab_mats);
    void calculateStabilityMatrixElastic(const std::vector<aurostd::xvector<double>>& poly_coeffs_elas_tensor,
                                         const std::vector<std::vector<aurostd::xvector<double>>>& poly_exps_elas_tensor,
                                         const std::vector<aurostd::xvector<double>>& poly_coeffs_volume,
                                         const std::vector<std::vector<aurostd::xvector<double>>>& poly_exps_volume,
                                         const aurostd::xvector<double>& conc,
                                         aurostd::xmatrix<double>& stab_mat);
    void calculateCritTemp();
    bool isZeros();
    void calculateCritTempSC();
  };
} // namespace cumulants

namespace cumulants {
  void cumulantsMain(const aurostd::xoption& cumulants_opts);
  void displayUsage();
} // namespace cumulants

#endif
