//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                  Simon Divilov - Duke University 2022                   *
// *                                                                         *
//****************************************************************************
// Written by Simon Divilov, 2022
// simon.divilov@duke.edu
// 
#ifndef _AFLOW_QCA_H_
#define _AFLOW_QCA_H_

#include "aflow.h"
#include "aflow_pocc.h"

#define CONC_SHIFT 0.01 // concentration shift away from 0
#define QCA_FILE_PREFIX string("aflow_qca_")
#define QCA_AFLOW_TAG string("[AFLOW_QCA]")

namespace qca {
  class QuasiChemApproxCalculator : public xStream {
    public:
      // constructors - START
      QuasiChemApproxCalculator(ostream& oss=cout);
      QuasiChemApproxCalculator(ofstream& FileMESSAGE, ostream& oss=cout);
      QuasiChemApproxCalculator(const aurostd::xoption& qca_flags, ostream& oss=cout);
      QuasiChemApproxCalculator(const aurostd::xoption& qca_flags, ofstream& FileMESSAGE, ostream& oss=cout);
      QuasiChemApproxCalculator(const QuasiChemApproxCalculator& b);
      // constructors - END
      ~QuasiChemApproxCalculator();
      const QuasiChemApproxCalculator& operator=(const QuasiChemApproxCalculator& b);
      void clear();

      // initializers - START
      bool initialize(ostream& oss);
      bool initialize(ofstream& FileMESSAGE, ostream& oss);
      bool initialize();
      bool initialize(const aurostd::xoption& qca_flags, ostream& oss);
      bool initialize(const aurostd::xoption& qca_flags, ofstream& FileMESSAGE, ostream& oss);
      bool initialize(const aurostd::xoption& qca_flags);
      // initializers - END

      void printParams();
      void errorFix();
      void calculateBinodal();
      void writeData();
      void readData();
      void plotData();

    private:
      void free();
      void copy(const QuasiChemApproxCalculator& b);

      bool m_initialized; ///< QuasiChemApproxCalculator is initialized
      _aflags m_aflags; ///< AFLOW flags
      uint min_sleep; ///< minimum number of seconds to sleep in between checks of ATAT
      string print; ///< output format
      bool screen_only; ///< output the data to terminal only
      bool image_only; ///< read the written data and creates and image
      bool use_sg; ///< compare initial and final xstructures only using their space groups
      bool calc_binodal; ///< calculate the binodal curve
      string rootdirpath; ///< path to the root directory where all calculations will be made
      string aflowlibpath; ///< path to the parent directory where the aflowlib output files are stored
      string plattice; ///< parent lattice of the alloy
      vector<string> elements; ///< elements in the alloy 
      int aflow_max_num_atoms; ///< maximum number of atoms in the cluster expansion calculated by AFLOW
      int max_num_atoms; ///< maximum number of atoms in the cluster expansion
      int conc_npts; ///< number of points used to evaluate the macroscopic concentration
      bool conc_curve; ///< the macroscopic concentration curve is defined
      vector<double> conc_curve_range; ///< concentration endpoints used to evaluate the macroscopic concentration, DIM: 2*Ne
      int temp_npts; ///< number of points used to evaluate the temperature range
      vector<double> temp_range; ///< temperature endpoints used to evaluate the temperature range, UNIT: K | DIM: 2
      double cv_cut; ///< coefficient of variation cut-off, UNIT: eV
      string alloyname; ///< elemental name of the alloy
      string rundirpath; ///< path to the directory where AFLOW is running
      vector<xstructure> vstr_aflow; ///< xstructures from AFLOW runs
      string lat_atat; ///< ATAT lattice
      vector<xstructure> vstr_ce; ///< xstructures from cluster expansion ATAT runs
      vector<int> mapstr; ///< xstructure map between AFLOW and ATAT
      xvector<int> skipstr; ///< skip xstructure match if the value is non-zero
      unsigned long int nelem; ///< number of unique elements
      unsigned long int ncluster; ///< number of clusters, DIM: Nc
      unsigned long int nconc; ///< number of concentration points, DIM: Nx
      xvector<double> beta; ///< inverse temperature, UNIT: unitless | DIM: Nt
      xmatrix<double> soln0; ///< initial guess for the probability distribution
      double cv_cluster; ///< coefficient of variation of the cluster expansion, UNIT: eV
      xvector<int> num_atom_cluster; ///< number of atoms of the clusters, DIM: Nc
      xmatrix<double> conc_cluster; ///< concentration of the clusters, UNIT: unitless | DIM: Nc, Ne
      xmatrix<double> num_elem_cluster; ///< number of each element of the clusters, UNIT: unitless | DIM: Nc, Ne
      xvector<double> excess_energy_cluster; ///< excess energy of the clusters, UNIT: eV | DIM: Nc
      xvector<unsigned long long int> degeneracy_cluster; ///< degeneracy of the clusters, DIM: Nc
      xmatrix<double> conc_macro; ///< macroscopic concentration of the alloy, UNIT: unitless | DIM: Nx, Ne
      xvector<double> temp; ///< temperature range, UNIT: K | DIM: Nt
      xmatrix<double> prob_ideal_cluster; ///< ideal (high-T) probability of the clusters as a function of concentration and temperature, DIM: Nx, Nc
      vector<xmatrix<double>> prob_cluster; ///< equilibrium probability of the clusters as a function of concentration and temperature, DIM: Nx, Nc, Nt
      std::pair<double, double> param_ec; ///< relative entropy at the equi-concentration and the transition temperature, UNIT: unitless, K
      xmatrix<double> rel_s; ///< relative entropy as a function of concentration and temperature, UNIT: unitless | DIM: Nx, Nt
      xvector<double> binodal_curve; ///< binodal curve as a function of the concentration, UNIT: K | DIM: Nx


      void readQCAFlags(const aurostd::xoption& qca_flags);
      string getLatForATAT(bool scale=false) const;
      void readAFLOWXstructures();
      vector<xstructure> getATATXstructures(const int max_num_atoms=0, bool fromfile=false) const;
      void calculateMapForXstructures(const vector<xstructure>& vstr1, const vector<xstructure>& vstr2);
      void generateFilesForATAT();
      void runATAT();
      void readCVCluster();
      void calculateNumAtomCluster();
      void calculateConcentrationCluster();
      void readExcessEnergyCluster();
      void setCongruentClusters();
      void calculateDegeneracyCluster();
      void calculateConcentrationMacro();
      void calculateTemperatureRange();
      void calculateProbabilityIdealCluster();
      void checkProbabilityIdeal() const;
      void calculateProbabilityCluster1D(int iix, const int it);
      void calculateProbabilityClusterND(int iix, const int it);
      void calculateProbabilityCluster();
      double getProbabilityConstraint(const int it, const int ix, const int ie, const int ideq, const xvector<double>& xvar) const;
      void checkProbabilityEquilibrium() const;
      void calculateRelativeEntropyEC();
      void calculateRelativeEntropy();
      void calculateBinodalCurve();
  };
}

namespace qca {
  void quasiChemicalApprox(const aurostd::xoption& vpflow);
  void displayUsage();
}

#endif
