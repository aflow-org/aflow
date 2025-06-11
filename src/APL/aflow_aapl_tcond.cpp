// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *            Aflow MARCO ESTERS - Duke University 2019-2021               *
// *                                                                         *
// ***************************************************************************
// Written by Marco Esters, 2019.
//
// This class calculates the thermal conductivity of a material using the
// Boltzmann Transport Equation (BTE).

#include "aflow_apl.h"

#define _DEBUG_AAPL_TCOND_ false

static const int max_iter = 250;  // Maximum number of iterations for the iterative BTE solution
static const aurostd::xcomplex<double> iONE(0.0, 1.0);  // imaginary number
static const double TCOND_ITER_THRESHOLD = 1e-4;  // Convergence criterion for thermal conductivity

using aurostd::xcomplex;
using aurostd::xmatrix;
using aurostd::xvector;
using aurostd::xerror;
using std::vector;
using namespace std::placeholders;

/************************** CONSTRUCTOR/DESTRUCTOR **************************/

namespace apl {

  //Constructor///////////////////////////////////////////////////////////////
  TCONDCalculator::TCONDCalculator() {
    free();
  }

  TCONDCalculator::TCONDCalculator(PhononCalculator& pc, const aurostd::xoption& opts) {
    _pc = &pc;
    _qm = &_pc->getQMesh();  // This pointer is only defined to make the code more readable
    _pc_set = true;
    initialize(opts);
  }

  //Copy Constructor//////////////////////////////////////////////////////////
  TCONDCalculator::TCONDCalculator(const TCONDCalculator& that) {
    if (this != &that) free();
    copy(that);
  }

  TCONDCalculator& TCONDCalculator::operator=(const TCONDCalculator& that) {
    if (this != &that) free();
    copy(that);
    return *this;
  }

  void TCONDCalculator::copy(const TCONDCalculator& that) {
    if (this == &that) return;
    boundary_grain_size = that.boundary_grain_size;
    calc_boundary = that.boundary_grain_size;
    calc_cumulative = that.calc_cumulative;
    calc_isotope = that.calc_isotope;
    calc_rta_only = that.calc_rta_only;
    eigenvectors = that.eigenvectors;
    freq = that.freq;
    grueneisen_avg = that.grueneisen_avg;
    grueneisen_mode = that.grueneisen_mode;
    gvel = that.gvel;
    _initialized = that._initialized;
    intr_trans_probs = that.intr_trans_probs;
    intr_trans_probs_iso = that.intr_trans_probs_iso;
    nBranches = that.nBranches;
    nIQPs = that.nIQPs;
    nQPs = that.nQPs;
    _pc = that._pc;
    _pc_set = that._pc_set;
    phase_space = that.phase_space;
    processes = that.processes;
    processes_iso = that.processes_iso;
    _qm = that._qm;
    scattering_rates_anharm = that.scattering_rates_anharm;
    scattering_rates_boundary = that.scattering_rates_boundary;
    scattering_rates_isotope = that.scattering_rates_isotope;
    scattering_rates_total = that.scattering_rates_total;
    temperatures = that.temperatures;
    thermal_conductivity = that.thermal_conductivity;
  }

  //Destructor////////////////////////////////////////////////////////////////
  TCONDCalculator::~TCONDCalculator() {
    free();
  }

  //free//////////////////////////////////////////////////////////////////////
  void TCONDCalculator::free() {
    boundary_grain_size = 0.0;
    calc_boundary = false;
    calc_cumulative = false;
    calc_isotope = false;
    calc_rta_only = false;
    eigenvectors.clear();
    freq.clear();
    grueneisen_avg.clear();
    grueneisen_mode.clear();
    gvel.clear();
    _initialized = false;
    intr_trans_probs.clear();
    intr_trans_probs_iso.clear();
    nBranches = 0;
    nIQPs = 0;
    nQPs = 0;
    _pc = NULL;
    _pc_set = false;
    phase_space.clear();
    processes.clear();
    processes_iso.clear();
    _qm = NULL;
    scattering_rates_anharm.clear();
    scattering_rates_boundary.clear();
    scattering_rates_isotope.clear();
    scattering_rates_total.clear();
    temperatures.clear();
    thermal_conductivity.clear();
  }

  //clear/////////////////////////////////////////////////////////////////////
  void TCONDCalculator::clear(PhononCalculator& pc) {
    free();
    _pc = &pc;
    _qm = &_pc->getQMesh();
    _pc_set = true;
  }

  void TCONDCalculator::initialize(const aurostd::xoption& opts) {
    string message = "";
    if (!_pc_set) {
      message = "PhononCalculator pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    // Set up phonon parameters
    nBranches = _pc->getNumberOfBranches();
    vector<int> grid;
    aurostd::string2tokens(opts.getattachedscheme("THERMALGRID"), grid, " xX");
    if (grid.size() != 3) {
      message = "Wrong format for the q-point grid for thermal conductivity calculations.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    _pc->initialize_qmesh(grid, true, true);
    _qm->makeIrreducible();
    nQPs = _qm->getnQPs();
    nIQPs = _qm->getnIQPs();

    // Set up temperatures
    double tstart = aurostd::string2utype<double>(opts.getattachedscheme("TSTART"));
    double tend = aurostd::string2utype<double>(opts.getattachedscheme("TEND"));
    double tstep = aurostd::string2utype<double>(opts.getattachedscheme("TSTEP"));

    if (tstart > tend) {
      message = "Tstart cannot be higher than Tend.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }
    if (tstep < _ZERO_TOL_) {
      message = "Tstep cannot be zero or negative.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }
    for (double t = tstart; t <= tend; t += tstep) temperatures.push_back(t);

    // Set up calculation options
    calc_boundary = opts.flag("BOUNDARY");
    if (calc_boundary) boundary_grain_size = aurostd::string2utype<double>(opts.getattachedscheme("NANO_SIZE"));
    calc_cumulative = opts.flag("CUMULATIVEK");
    calc_isotope = opts.flag("ISOTOPE");
    calc_rta_only = (opts.getattachedscheme("BTE") == "RTA");

    // Frequencies and group velocities
    message = "Calculating frequencies and group velocities.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());
    gvel = _pc->calculateGroupVelocitiesOnMesh(freq, eigenvectors);
    _initialized = true;
  }

}  // namespace apl

/******************************* MAIN FUNCTION ******************************/

namespace apl {

  //calculateThermalConductvity///////////////////////////////////////////////
  // The main function that calculates the thermal conductivity tensor, the
  // Grueneisen parameters, and the scattering phase space.
  void TCONDCalculator::calculateThermalConductivity() {
    if (!_initialized) {
      string message = "Not initialized.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }

    // Transition probabilities
    calculateTransitionProbabilities();

    // Thermal conductivity tensor and scattering rates
    uint ntemps = temperatures.size();
    thermal_conductivity.clear();
    thermal_conductivity.assign(ntemps, xmatrix<double>(3, 3));
    scattering_rates_anharm.clear();
    scattering_rates_anharm.resize(ntemps, vector<vector<double> >(nIQPs, vector<double> (nBranches)));
    scattering_rates_total.clear();
    scattering_rates_total.resize(ntemps, vector<vector<double> >(nIQPs, vector<double> (nBranches)));
    // Only need little groups for full BTE
    if (!calc_rta_only && !_qm->littleGroupsCalculated()) _qm->calculateLittleGroups();

    for (uint t = 0; t < temperatures.size(); t++) {
      thermal_conductivity[t] = calculateThermalConductivityTensor(temperatures[t], scattering_rates_total[t], scattering_rates_anharm[t]);
    }

  }

  void TCONDCalculator::writeOutputFiles(const string& directory) {
    // Frequencies and group velocities
    string filename = aurostd::CleanFileName(directory + "/" + DEFAULT_AAPL_FILE_PREFIX + DEFAULT_AAPL_GVEL_FILE);
    _pc->writeGroupVelocitiesToFile(filename, gvel, freq, "THz");

    // Irreducible q-points
    filename = directory + "/" + DEFAULT_AAPL_FILE_PREFIX + DEFAULT_AAPL_IRRQPTS_FILE;
    _qm->writeIrredQpoints(filename);

    // Grueneisen parameters
    filename = aurostd::CleanFileName(directory + "/" + DEFAULT_AAPL_FILE_PREFIX + DEFAULT_AAPL_GRUENEISEN_FILE);
    writeGrueneisen(filename);

    // Phase space
    filename = aurostd::CleanFileName(directory + "/" + DEFAULT_AAPL_FILE_PREFIX + DEFAULT_AAPL_PS_FILE);
    writePhaseSpace(filename);

    // Scattering rates
    filename = aurostd::CleanFileName(directory + "/" + DEFAULT_AAPL_FILE_PREFIX + DEFAULT_AAPL_RATES_FILE);
    writeTempDepOutput(filename, "SCATTERING_RATES", "1/ps", temperatures, scattering_rates_total);
    filename = aurostd::CleanFileName(directory + "/" + DEFAULT_AAPL_FILE_PREFIX + DEFAULT_AAPL_RATES_3RD_FILE);
    writeTempDepOutput(filename, "SCATTERING_RATES_ANHARMONIC", "1/ps", temperatures, scattering_rates_anharm);
    if (calc_isotope) {
      filename = aurostd::CleanFileName(directory + "/" + DEFAULT_AAPL_FILE_PREFIX + DEFAULT_AAPL_ISOTOPE_FILE);
      writeTempIndepOutput(filename, "SCATTERING_RATES_ISOTOPE", "1/ps", scattering_rates_isotope);
    }
    if (calc_boundary) {
      string filename = aurostd::CleanFileName(directory + "/" + DEFAULT_AAPL_FILE_PREFIX + DEFAULT_AAPL_BOUNDARY_FILE);
      writeTempIndepOutput(filename, "SCATTERING_RATES_ISOTOPE", "1/ps", scattering_rates_boundary);
    }

    // Thermal conductivity
    filename = aurostd::CleanFileName(directory + "/" + DEFAULT_AAPL_FILE_PREFIX + DEFAULT_AAPL_TCOND_FILE);
    writeThermalConductivity(filename);
  }


} // namespace apl

/******************************** GRUENEISEN ********************************/

namespace apl {

  //calculateGrueneisenParameters/////////////////////////////////////////////
  // Calculates the mode and average Grueneisen parameters.
  void TCONDCalculator::calculateGrueneisenParameters() {
    if (!_initialized) {
      string message = "Not initialized.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    // Grueneisen parameters
    string message = "Calculating Grueneisen parameters.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());

    vector<vector<vector<xcomplex<double> > > > phases = calculatePhases(false);  // false: not conjugate
    grueneisen_mode = calculateModeGrueneisen(phases);
    phases.clear();
    grueneisen_avg.clear();
    grueneisen_avg.resize(temperatures.size());
    for (uint t = 0; t < temperatures.size(); t++) {
      grueneisen_avg[t] = calculateAverageGrueneisen(temperatures[t], grueneisen_mode);
    }
  }

  //calculateModeGrueneisen///////////////////////////////////////////////////
  // Calculates the Grueneisen parameters for each mode.
  vector<vector<double> > TCONDCalculator::calculateModeGrueneisen(const vector<vector<vector<xcomplex<double> > > >& phases) {
    // Prepare and precompute
    vector<vector<double> > grueneisen(nIQPs, vector<double>(nBranches));

    const vector<vector<double> >& ifcs = _pc->getAnharmonicForceConstants(3);
    const Supercell& scell = _pc->getSupercell();

    // Inverse masses
    const vector<vector<int> >& clusters = _pc->getClusters(3);
    uint nclusters = clusters.size();
    vector<double> invmasses(nclusters);
    for (uint c = 0; c < nclusters; c++) {
      double mass = 1.0;
      for (int i = 0; i < 2; i++) mass *= scell.getAtomMass(clusters[c][i]);
      invmasses[c] = 1/sqrt(mass);
    }

    // Cartesian indices to avoid running xcombos multiple times
    vector<vector<int> > cart_indices;
    aurostd::xcombos cart(3, 3, 'E', true);
    while (cart.increment()) cart_indices.push_back(cart.getCombo());
    uint ncart = cart_indices.size();

    // Prepare precomputation of eigenvalue products
    int natoms = (int) _pc->getInputCellStructure().atoms.size();
    vector<int> atpowers(2, 1);
    vector<vector<int> > at_eigen;
    aurostd::xcombos at_combos(natoms, 2, 'E' , true);
    while (at_combos.increment()) at_eigen.push_back(at_combos.getCombo());
    uint nateigen = at_eigen.size();
    vector<vector<xcomplex<double> > > eigenprods(nateigen, vector<xcomplex<double> >(ncart));

    // Get distances to avoid running minimizeDistance multiple times
    const xstructure& scell_xstr = scell.getSupercellStructure();
    uint natoms_sc = scell_xstr.atoms.size();
    vector<vector<xvector<double> > > min_dist(natoms, vector<xvector<double> >(natoms_sc));
    for (int i = 0; i < natoms; i++) {
      const xvector<double>& iat_cpos = scell_xstr.atoms[scell.pc2scMap(i)].cpos;
      for (uint j = 0; j < natoms_sc; j++) {
        min_dist[i][j] = SYM::minimizeDistanceCartesianMethod(scell_xstr.atoms[j].cpos, iat_cpos, scell_xstr.lattice);
      }
    }

    // Initialize variables
    int at1_pc = 0, at2_sc = 0, at2_pc = 0, at3_sc = 0, e = 0, q = 0;
    double ifc_prod = 0.0;
    uint c = 0, crt = 0;
    xcomplex<double> prefactor, eigen, g_mode;

    // Start calculation
    for (int iq = 0; iq < nIQPs; iq++) {
      q = _qm->getIbzqpts()[iq];
      for (int br = 0; br < nBranches; br++) {
        if (freq[q][br] > _FLOAT_TOL_) {
          g_mode.re = 0.0;
          g_mode.im = 0.0;

          // Precompute eigenvalue products
          for (c = 0; c < nateigen; c++) {
            for (crt = 0; crt < ncart; crt++) {
              e = at_eigen[c][0] * 3 + cart_indices[crt][0] + 1;
              eigen = conj(eigenvectors[q][e][br + 1]);
              e = at_eigen[c][1] * 3 + cart_indices[crt][1] + 1;
              eigen *= eigenvectors[q][e][br + 1];
              eigenprods[c][crt] = eigen;
            }
          }

          for (c = 0; c < nclusters; c++) {
            at1_pc = scell.sc2pcMap(clusters[c][0]);
            at2_sc = clusters[c][1];
            at2_pc = scell.sc2pcMap(at2_sc);
            at3_sc = clusters[c][2];
            prefactor = invmasses[c] * phases[at1_pc][at2_sc][q];
            e = at1_pc * natoms + at2_pc;
            for (crt = 0; crt < ncart; crt++) {
              ifc_prod = ifcs[c][crt] * min_dist[at1_pc][at3_sc][cart_indices[crt][2] + 1];
              // Perform multiplication explicitly in place instead of using xcomplex.
              // This is three times faster because constructors and destructors are not called.
              g_mode.re += ifc_prod * (prefactor.re * eigenprods[e][crt].re - prefactor.im * eigenprods[e][crt].im);
              g_mode.im += ifc_prod * (prefactor.re * eigenprods[e][crt].im + prefactor.im * eigenprods[e][crt].re);
            }
          }
          g_mode *= -10.0 * au2nmTHz/(6.0 * std::pow(freq[q][br], 2));
          if (g_mode.im > _FLOAT_TOL_) {  // _ZERO_TOL_ is too tight
            stringstream message;
            message << "Grueneisen parameter at mode " << iq << ", " << br << " is not real (" << g_mode.re << ", " << g_mode.im << ").";
            pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS(), _LOGGER_WARNING_);
          }
          grueneisen[iq][br] = g_mode.re;
        } else {
          grueneisen[iq][br] = 0.0;
        }
      }
    }

    return grueneisen;
  }

  //calculateAverageGrueneisen////////////////////////////////////////////////
  // Calculates the average Grueneisen parameter for a specific temperature.
  double TCONDCalculator::calculateAverageGrueneisen(double T,
      const vector<vector<double> >& grueneisen_mode) {
    vector<vector<double> > occ = getOccupationNumbers(T);
    int iq = 0;
    double c, c_tot = 0, g_tot = 0;
    double prefactor = 1E24 * std::pow(PLANCKSCONSTANT_hbar_THz, 2)/(KBOLTZ * std::pow(T, 2));
    for (int q = 0; q < nQPs; q++) {
      iq = _qm->getIbzqpt(q);
      for (int br = 0; br < nBranches; br++) {
        if (freq[q][br] > _FLOAT_TOL_) {
          c = prefactor * occ[q][br] * (1.0 + occ[q][br]) * std::pow(freq[q][br], 2);
          c_tot += c;
          g_tot += grueneisen_mode[iq][br] * c;
        }
      }
    }
    return (g_tot/c_tot);
  }

}  // namespace apl


/************************* TRANSITION PROBABILITIES *************************/

namespace apl {

  //calculateTransitionProbabilities//////////////////////////////////////////
  // Calculates the transition probabilities for three-phonon, isotope, and
  // boundary scattering. Also calculates the scattering phase space.
  void TCONDCalculator::calculateTransitionProbabilities() {
    string message = "Calculating transition probabilities.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());
#ifdef AFLOW_MULTITHREADS_ENABLE
    xthread::xThread xt(_pc->getNCPUs());
#endif
    _qm->generateTetrahedra();
    // The conjugate is necessary because the three-phonon scattering processes
    // will be calculated for g - q' - q" = G
    vector<vector<vector<xcomplex<double> > > > phases = calculatePhases(true);  // true: conjugate

    // Three-phonon transition probabilities
    message = "Calculating transition probabilities for 3-phonon scattering processes.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());
    message = "Transition Probabilities";
    processes.clear();
    processes.resize(nIQPs);
    intr_trans_probs.clear();
    intr_trans_probs.resize(nIQPs);
    // Phase space for each (1) q-point, (2) branch, (3) type (AAA, AAO, etc.), and (4) (normal, umklapp)
    phase_space.clear();
    phase_space.resize(nIQPs, vector<vector<vector<double> > >(nBranches, vector<vector<double> >(4, vector<double>(2, 0.0))));
#ifdef AFLOW_MULTITHREADS_ENABLE
    xt.setProgressBar(*_pc->getOSS());
    std::function<void(int, vector<vector<vector<vector<double> > > >&,
        const vector<vector<vector<xcomplex<double> > > >&)> fn = std::bind(&TCONDCalculator::calculateTransitionProbabilitiesPhonon, this, _1, _2, _3);
    xt.run(nIQPs, fn, phase_space, phases);
    xt.unsetProgressBar();
#else
    for (int i = 0; i < nIQPs; i++) calculateTransitionProbabilitiesPhonon(i, phase_space, phases);
#endif

    if (calc_isotope) {
      message = "Calculating isotope transition probabilities.";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());

      // Test if isotope scattering is possible
      const xstructure& pcell = _pc->getInputCellStructure();
      uint natoms = pcell.atoms.size();
      uint at = 0;
      for (at = 0; at < natoms; at++) {
        if (GetPearsonCoefficient(pcell.atoms[at].atomic_number) > 0) break;
      }

      if (at == natoms) {
        calc_isotope = false;
        message = "There are no atoms with isotopes of different masses. Isotope scattering will be turned off.";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());
      } else {
        processes_iso.resize(nIQPs);
        intr_trans_probs_iso.resize(nIQPs);
        scattering_rates_isotope.clear();
        scattering_rates_isotope.resize(nIQPs, vector<double>(nBranches));

#ifdef AFLOW_MULTITHREADS_ENABLE
        std::function<void(int)> fn = std::bind(&TCONDCalculator::calculateTransitionProbabilitiesIsotope, this, _1);
        xt.run(nIQPs, fn);
#else
        for (int i = 0; i < nIQPs; i++) calculateTransitionProbabilitiesIsotope(i);
#endif
      }
    }

    if (calc_boundary) {
      scattering_rates_boundary.clear();
      scattering_rates_boundary = calculateTransitionProbabilitiesBoundary();
    }
  }

  //calculatePhases///////////////////////////////////////////////////////////
  // Calculates the phase factors for each atom in the supercell. The
  // conjugate is used for the scattering matrices whereas the non-conjugate
  // is used for the Grueneisen parameters. Calculating the phases ahead of
  // time decreases runtimes considerably.
  vector<vector<vector<xcomplex<double> > > > TCONDCalculator::calculatePhases(bool conjugate) {
    vector<xvector<double> > qpts = _pc->getQMesh().getQPointsCPOS();
    const xstructure& scell = _pc->getSuperCellStructure();
    const xstructure& pcell = _pc->getInputCellStructure();
    const vector<int>& sc2pcMap = _pc->getSupercell()._sc2pcMap;
    const vector<int>& pc2scMap = _pc->getSupercell()._pc2scMap;
    uint niatoms = pcell.atoms.size();
    uint natoms = scell.atoms.size();
    vector<vector<vector<xcomplex<double> > > > phases(niatoms, vector<vector<xcomplex<double> > >(natoms, vector<xcomplex<double> >(nQPs)));

    int at_eq = 0, at_eq_sc = 0, iat_sc = 0;
    xvector<double> min_vec(3);
    for (uint iat = 0; iat < niatoms; iat++) {
      iat_sc = pc2scMap[iat];
      const xvector<double>& iat_cpos = scell.atoms[iat_sc].cpos;
      for (uint at = 0; at < natoms; at++) {
        min_vec = SYM::minimizeDistanceCartesianMethod(scell.atoms[at].cpos, iat_cpos, scell.lattice);
        at_eq = sc2pcMap[at];
        at_eq_sc = pc2scMap[at_eq];
        min_vec += iat_cpos;
        min_vec -= scell.atoms[at_eq_sc].cpos;
        for (int q = 0; q < nQPs; q++) {
          if (conjugate) phases[iat][at][q] = exp(-iONE * scalar_product(qpts[q], min_vec));
          else phases[iat][at][q] = exp(iONE * scalar_product(qpts[q], min_vec));
        }
      }
    }
    return phases;
  }

  //calculateTransitionProbabilitiesPhonon////////////////////////////////////
  // Calculates the intrinsic transition probabilities and the scattering
  // phase space for three-phonon scattering processes. It uses the inversion
  // symmetry of the q-point grid and the transposition symmetry of the
  // scattering matrix elements to reduce the computational cost.
  void TCONDCalculator::calculateTransitionProbabilitiesPhonon(
      int i,
      vector<vector<vector<vector<double> > > >& phase_space,
      const vector<vector<vector<xcomplex<double> > > >& phases) {
    // Prepare and precompute
    const Supercell& scell = _pc->getSupercell();
    const vector<vector<int> >& clusters = _pc->getClusters(3);
    uint nclusters = clusters.size();

    // Inverse masses
    vector<double> invmasses(nclusters);
    for (uint c = 0; c < nclusters; c++) {
      double mass = 1.0;
      for (int o = 0; o < 3; o++) mass *= scell.getAtomMass(clusters[c][o]);
      invmasses[c] = 1/sqrt(mass);
    }

    // Cartesian indices to avoid running xcombos multiple times
    const vector<vector<double> >& ifcs = _pc->getAnharmonicForceConstants(3);
    vector<vector<int> > cart_indices;
    aurostd::xcombos cart(3, 3, 'E', true);
    while (cart.increment()) cart_indices.push_back(cart.getCombo());
    uint ncart = cart_indices.size();

    // Branch indices to avoid running xcombos multiple times
    vector<vector<int> > branches;
    aurostd::xcombos branch_combos(nBranches, 3, 'E', true);
    while (branch_combos.increment()) branches.push_back(branch_combos.getCombo());
    uint nbr = branches.size();

    // Units are chosen so that probabilities are in THz (1/ps)
    const double probability_prefactor = std::pow(au2nmTHz * 10.0, 2) * PLANCKSCONSTANTAMU_hbar_THz * PI/4.0;
    const double ps_prefactor = 2.0/(3.0 * std::pow(nBranches, 3) * nQPs);

    // Prepare precomputation of eigenvalue products
    int natoms = (int) _pc->getInputCellStructure().atoms.size();
    vector<int> atpowers(3, 1);
    vector<vector<int> > at_eigen;
    aurostd::xcombos at_combos(natoms, 3, 'E' , true);
    while (at_combos.increment()) at_eigen.push_back(at_combos.getCombo());
    uint nateigen = at_eigen.size();
    vector<vector<xcomplex<double> > > eigenprods(nateigen, vector<xcomplex<double> >(ncart));
    for (int j = 0; j < 2; j++) atpowers[j] = aurostd::powint(natoms, 2 - j);

    // Precompute the indices of -q to use for inversion symmetry
    vector<int> q_minus(nQPs);
    for (int q = 0; q < nQPs; q++) q_minus[q] = _qm->getQPointIndex(-_qm->getQPoint(q).fpos);

    // Initialize variables
    xcomplex<double> matrix, prefactor, eigen;
    vector<vector<double> > weights(3, vector<double>(nQPs)), frequencies(3, vector<double>(nQPs));
    vector<int> qpts(3), proc(2), lastq(nQPs);
    xvector<double> fpos_diff(3);
    vector<bool> is_umklapp(nQPs);
    int iat = 0, j = 0, e = 0, q = 0, p = 0, w = 0, lq = 0, b = 0;
    uint c = 0, crt = 0, br = 0;
    double transprob = 0.0, freq_ref = 0.0, prod = 0.0;
    bool calc = true;

    // Start calculation
    qpts[0] = _qm->getIbzqpts()[i];
    // Get the q-point q" that fulfills q - q' - q" = G. Due to the inversion
    // symmetry of the q-point grid, q + q' - q" = G does not need to be
    // evaluated since for each q' there is also a (-q') on the grid.
    for (q = 0; q < nQPs; q++) {
      fpos_diff = _qm->getQPoint(qpts[0]).fpos - _qm->getQPoint(q).fpos;
      is_umklapp[q] = !inCell(fpos_diff, _ZERO_TOL_, 0.5, -0.5);
      lastq[q] = _qm->getQPointIndex(fpos_diff);
    }

    for (br = 0; br < nbr; br++) {
      freq_ref = freq[qpts[0]][branches[br][0]];
      // Prepare weight calculation for d(w +/- w' - w"). The first two terms
      // describe the + process whereas the last term describes the - process.
      // The + process requires two calculations to exploit transposition
      // symmetry as d(w + w'- w") and d(w + w" - w') in general do not have
      // the same weight.
      for (q = 0; q < nQPs; q++) {
        lq = lastq[q];
        frequencies[0][q] = -freq[q][branches[br][1]] + freq[lq][branches[br][2]];
        frequencies[1][q] =  freq[q][branches[br][1]] - freq[lq][branches[br][2]];
        frequencies[2][q] =  freq[q][branches[br][1]] + freq[lq][branches[br][2]];
      }

      for (j = 0; j < 3; j++) getWeightsLT(freq_ref, frequencies[j], weights[j]);

      for (q = 0; q < nQPs; q++) {
        lq = lastq[q];
        // Transposition symmetry: only use the processes that are unique permuations
        calc = (q < lq);
        // If the integration weights of all scattering processes are zero,
        // there is no need to calculate the scattering matrix.
        if (calc) {
          for (j = 0; j < 3; j++) {
            if (weights[j][q] > _ZERO_TOL_) break;
          }
          calc = (j < 3);
        }
        // Calculate contributions to the scattering phase space
        if (calc) {
          qpts[1] = q;
          qpts[2] = lq;
          p = 0;
          for (j = 0; j < 3; j++) {
            if (branches[br][j] > 2) p++;
          }
          j = (is_umklapp[q_minus[q]]?1:0);
          phase_space[i][branches[br][0]][p][j] += weights[0][q];
          j = (is_umklapp[q_minus[lq]]?1:0);
          phase_space[i][branches[br][0]][p][j] += weights[1][q];
          // No need for the factor 1/2 since permutations are eliminated.
          // This intrinsically prevents double-counting.
          j = (is_umklapp[q]?1:0);
          phase_space[i][branches[br][0]][p][j] += weights[2][q];
          // If any frequency in the process is zero or not real, the process
          // does not contribute to the thermal conductivity tensor (it does
          // contribute to the scattering phase space though).
          for (j = 0; j < 3; j++) {
            if (freq[qpts[j]][branches[br][j]] < _FLOAT_TOL_) break;
          }
          calc = (j == 3);
        }

        if (calc) {
          // Precompute eigenvalue products
          for (c = 0; c < nateigen; c++) {
            for (crt = 0; crt < ncart; crt++) {
              e = at_eigen[c][0] * 3 + cart_indices[crt][0] + 1;
              eigen = eigenvectors[qpts[0]][e][branches[br][0] + 1];
              for (j = 1; j < 3; j++) {
                e = at_eigen[c][j] * 3 + cart_indices[crt][j] + 1;
                eigen *= conj(eigenvectors[qpts[j]][e][branches[br][j] + 1]);
              }
              eigenprods[c][crt] = eigen;
            }
          }

          // Calculate scattering matrix for the process
          matrix.re = 0.0;
          matrix.im = 0.0;
          for (c = 0; c < nclusters; c++) {
            iat = scell.sc2pcMap(clusters[c][0]);
            prefactor.re = invmasses[c];
            prefactor.im = 0.0;
            for (j = 1; j < 3; j++) prefactor *= phases[iat][clusters[c][j]][qpts[j]];
            e = 0;
            for (j = 0; j < 3; j++) e += scell.sc2pcMap(clusters[c][j]) * atpowers[j];
            for (crt = 0; crt < ncart; crt++) {
              // Perform multiplication explicitly in place instead of using xcomplex.
              // This is three times faster because constructors and destructors are not called.
              matrix.re += ifcs[c][crt] * (prefactor.re * eigenprods[e][crt].re - prefactor.im * eigenprods[e][crt].im);
              matrix.im += ifcs[c][crt] * (prefactor.re * eigenprods[e][crt].im + prefactor.im * eigenprods[e][crt].re);
            }
          }
          prod = magsqr(matrix);

          // Only go on if the scattering matrix is not zero.
          if (prod > _ZERO_TOL_) {
            for (j = 0; j < 3; j++) prod /= freq[qpts[j]][branches[br][j]];
            for (j = 0; j < 3; j++) {
              transprob = prod * probability_prefactor * weights[j][q];
              if (transprob > _ZERO_TOL_) {
                // The process information needs to be stored. For q-points, only
                // one index is necessary since the last index follows from
                // momentum conservation.
                // The branches will be stored in a combined index to save memory.
                // Instead of storing the sign as an integer, the q-point index
                // will be signed to save memory. Adding 1 to the index is done
                // to have a clear sign indication for q = 0. This will be reversed
                // in getProcess().
                if (j == 2) {  // -
                  proc[0] = -(q + 1);
                  proc[1] = br;
                  processes[i].push_back(proc);
                  intr_trans_probs[i].push_back(transprob);
                  // Transposition symmetry
                  proc[0] = -(lq + 1);
                  proc[1] = branches[br][0] * aurostd::powint(nBranches, 2) + branches[br][2] * nBranches + branches[br][1];
                  processes[i].push_back(proc);
                  intr_trans_probs[i].push_back(transprob);
                } else {  // +
                  // Inversion symmetry
                  if (j == 0) {
                    proc[0] = q_minus[q] + 1;
                    proc[1] = br;
                  } else{
                    proc[0] = q_minus[lq] + 1;
                    proc[1] = branches[br][0] * aurostd::powint(nBranches, 2) + branches[br][2] * nBranches + branches[br][1];
                  }
                  processes[i].push_back(proc);
                  intr_trans_probs[i].push_back(transprob);
                }
              }
            }
          }
        }
      }
    }
    // Finish up phase space calculation
    w = _qm->getWeights()[i];
    for (b = 0; b < nBranches; b++) {
      for (p = 0; p < 4; p++) {
        for (j = 0; j < 2; j++) {
          phase_space[i][b][p][j] *= ps_prefactor * w;
        }
      }
    }
  }

  //calculateTransitionProbabilitiesIsotope///////////////////////////////////
  // Calculates the intrinsic transition probabilities/scattering rates of
  // the isotope scattering processes.
  void TCONDCalculator::calculateTransitionProbabilitiesIsotope(int i) {
    // Prepare
    const xstructure& pcell = _pc->getInputCellStructure();
    uint natoms = pcell.atoms.size();
    vector<double> pearson(natoms);
    uint at = 0;
    for (at = 0; at < natoms; at++) pearson[at] = GetPearsonCoefficient(pcell.atoms[at].atomic_number);

    vector<double> frequencies(nQPs), weights(nQPs);
    vector<int> proc(2);
    double prefactor = 0.0, eigsqr = 0.0, rate = 0.0;
    int q1 = 0, q2 = 0, br1 = 0, br2 = 0, b = 0, e = 0;
    xcomplex<double> eig;

    q1 = _qm->getIbzqpts()[i];
    for (br1 = 0; br1 < nBranches; br1++) {
      prefactor = freq[q1][br1] * freq[q1][br1] * PI/2.0;
      for (br2 = 0; br2 < nBranches; br2++) {
        for (q2 = 0; q2 < nQPs; q2++) frequencies[q2] = freq[q2][br2];
        getWeightsLT(freq[q1][br1], frequencies, weights);
        for (q2 = 0; q2 < nQPs; q2++) {
          // Only processes with non-zero weights need to be considered.
          if (weights[q2] > _ZERO_TOL_) {
            eigsqr = 0.0;
            for (at = 0; at < natoms; at++) {
              if (pearson[at] != 0) {
                eig.re = 0.0;
                eig.im = 0.0;
                e = 3 * at;
                for (int i = 1; i < 4; i++) {
                  // Perform multiplication explicitly in place instead of using xcomplex.
                  // This is three times faster because constructors and destructors are not called.
                  eig.re += eigenvectors[q1][e + i][br1 + 1].re * eigenvectors[q2][e + i][br2 + 1].re;
                  eig.re += eigenvectors[q1][e + i][br1 + 1].im * eigenvectors[q2][e + i][br2 + 1].im;
                  eig.im += eigenvectors[q1][e + i][br1 + 1].re * eigenvectors[q2][e + i][br2 + 1].im;
                  eig.im -= eigenvectors[q1][e + i][br1 + 1].im * eigenvectors[q2][e + i][br2 + 1].re;
                }
                eigsqr += pearson[at] * magsqr(eig);
              }
            }
            rate = prefactor * weights[q2] * eigsqr;
            // Store branches into combined index to save memory
            b = nBranches * br1 + br2;
            proc[0] = q2;
            proc[1] = b;
            intr_trans_probs_iso[i].push_back(rate);
            processes_iso[i].push_back(proc);
            scattering_rates_isotope[i][br1] += rate;
          }
        }
      }
    }
  }

  // calculateTransitionProbabilitiesBoundary/////////////////////////////////
  // Calculates the intrinsic transition probabilities/scattering rates of
  // the grain boundary scattering processes.
  vector<vector<double> > TCONDCalculator::calculateTransitionProbabilitiesBoundary() {
    int br = 0, iq = 0, q = 0;
    vector<vector<double> > rates(nIQPs, vector<double>(nBranches));

    stringstream message;
    message << "Calculating grain boundary transition probabilities with a grain size of " << boundary_grain_size << " nm.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());
    for (iq = 0; iq < nIQPs; iq++) {
      for (br = 0; br < nBranches; br++) {
        q = _qm->getIbzqpts()[iq];
        rates[iq][br] = aurostd::modulus(gvel[q][br])/boundary_grain_size;
      }
    }
    return rates;
  }

  //getWeightsLT//////////////////////////////////////////////////////////////
  // Calculate the integration weights using the linear tetrahedron method.
  // fij is a helper function that speeds up the calculation. It is not part
  // of the class to allow for more efficient inlining.
  // Moving this function into the LTMethod class caused a significant runtime
  // increase, so it should stay here until the speed issues can be resolved.
  double fij(double fi, double fj, double f) {
    return (f - fj)/(fi - fj);
  }

  void TCONDCalculator::getWeightsLT(double freq_ref, const vector<double>& frequencies, vector<double>& weights) {
    for (int q = 0; q < nQPs; q++) weights[q] = 0;
    const vector<vector<int> >& corners = _qm->getTetrahedra();
    double vol = _qm->getVolumePerTetrahedron();
    int ntet = _qm->getnTetrahedra();

    double g = 0.0, tmp = 0.0;
    int i = 0, j = 0, ii = 0, jj = 0;

    double f[4];
    int index_sort[4], index[4];

    for (i = 0; i < ntet; i++) {
      for (j = 0; j < 4; j++) {
        f[j] = frequencies[corners[i][j]];
        index[j] = corners[i][j];
      }

      // Sort
      for (ii = 0; ii < 4; ii++) index_sort[ii] = ii;

      for (ii = 0; ii < 4; ii++) {
        tmp = f[ii];
        jj = ii;
        while (jj > 0 && tmp < f[jj - 1]) {
          f[jj] = f[jj - 1];
          index_sort[jj] = index_sort[jj - 1];
          jj--;
        }
        f[jj] = tmp;
        index_sort[jj] = ii;
      }

      double f1 = f[0];
      double f2 = f[1];
      double f3 = f[2];
      double f4 = f[3];

      // This only happens when all corners have the same energy within machine
      // epsilon. In that case, freq_ref lies outside the tetrahedron and does
      // not contribute to the integration. However, due to epsilon, freq_ref
      // may be just between two of the frequencies in the tetrahedron, causing
      // numerical instabilities if this tetrahedron is not skipped.
      if (aurostd::isequal(f1, f4)) continue;

      int q1 = index[index_sort[0]];
      int q2 = index[index_sort[1]];
      int q3 = index[index_sort[2]];
      int q4 = index[index_sort[3]];

      double I1 = 0.0;
      double I2 = 0.0;
      double I3 = 0.0;
      double I4 = 0.0;

      if (f3 <= freq_ref && freq_ref < f4) {
        g = std::pow(f4 - freq_ref, 2) / ((f4 - f1) * (f4 - f2) * (f4 - f3));

        I1 = g * fij(f1, f4, freq_ref);
        I2 = g * fij(f2, f4, freq_ref);
        I3 = g * fij(f3, f4, freq_ref);
        I4 = g * (fij(f4, f1, freq_ref) + fij(f4, f2, freq_ref) + fij(f4, f3, freq_ref));

      } else if (f2 <= freq_ref && freq_ref < f3) {
        g = (f2 - f1 + 2.0 * (freq_ref - f2) - (f4 + f3 - f2 - f1)
            * std::pow(freq_ref - f2, 2) / ((f3 - f2) * (f4 - f2))) / ((f3 - f1) * (f4 - f1));

        I1 = g * fij(f1, f4, freq_ref) + fij(f1, f3, freq_ref) * fij(f3, f1, freq_ref) * fij(f2, f3, freq_ref) / (f4 - f1);
        I2 = g * fij(f2, f3, freq_ref) + std::pow(fij(f2, f4, freq_ref), 2) * fij(f3, f2, freq_ref) / (f4 - f1);
        I3 = g * fij(f3, f2, freq_ref) + std::pow(fij(f3, f1, freq_ref), 2) * fij(f2, f3, freq_ref) / (f4 - f1);
        I4 = g * fij(f4, f1, freq_ref) + fij(f4, f2, freq_ref) * fij(f2, f4, freq_ref) * fij(f3, f2, freq_ref) / (f4 - f1);

      } else if (f1 <= freq_ref && freq_ref < f2) {
        g = std::pow(freq_ref - f1, 2) / ((f2 - f1) * (f3 - f1) * (f4 - f1));

        I1 = g * (fij(f1, f2, freq_ref) + fij(f1, f3, freq_ref) + fij(f1, f4, freq_ref));
        I2 = g * fij(f2, f1, freq_ref);
        I3 = g * fij(f3, f1, freq_ref);
        I4 = g * fij(f4, f1, freq_ref);

      }
      weights[q1] += vol * I1;
      weights[q2] += vol * I2;
      weights[q3] += vol * I3;
      weights[q4] += vol * I4;
    }
  }

}  // namespace apl


/************************************ BTE ***********************************/

namespace apl {

  //getProcess////////////////////////////////////////////////////////////////
  // Information on scattering processes is stored in a sparse vector with
  // only one index for sign and q-points, and one for branches. This function
  // restores the full information.
  void TCONDCalculator::getProcess(const vector<int>& process, vector<int>& qpts,
      vector<int>& branches, int& sign) {
    // Signs and q-points
    // Signs are stored in the q-index instead of using an extra integer.
    xvector<double> qsum = _qm->getQPoint(qpts[0]).fpos;
    if (process[0] < 0) {
      sign = 1;
      qpts[1] = -process[0] - 1;
      qpts[2] = _qm->getQPointIndex(_qm->getQPoint(qpts[0]).fpos - _qm->getQPoint(qpts[1]).fpos);
    } else {
      sign = 0;
      qpts[1] = process[0] - 1;
      qpts[2] = _qm->getQPointIndex(qsum + _qm->getQPoint(qpts[1]).fpos);
    }

    // Branches
    int br = process[1];
    int pw = 0;
    for (int i = 0; i < 2; i++) {
      pw = aurostd::powint(nBranches, 2 - i);
      branches[i] = br/pw;
      br = br % pw;
    }
    branches[2] = br;
  }

  //calculateThermalConductivityTensor////////////////////////////////////////
  // Calculates the thermal conductivity tensor, and the total and anharmonic
  // scattering rates for a specific temperature. The rates are passed by
  // reference so that they can be written into output files later.
  xmatrix<double> TCONDCalculator::calculateThermalConductivityTensor(double T,
      vector<vector<double> >& rates_total,
      vector<vector<double> >& rates_anharm) {
    stringstream message;
    message << "Calculating thermal conductivity for " << T << " K.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());
    // Bose-Einstein distribution
    vector<vector<double> > occ = getOccupationNumbers(T);

    message << "Calculating scattering rates.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());
    rates_total = calculateTotalRates(occ, rates_anharm);

    message << "Calculating RTA";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());
    vector<vector<xvector<double> > > mfd = getMeanFreeDispRTA(rates_total);
    xmatrix<double> tcond = calcTCOND(T, occ, mfd); // RTA solution

    // Iteration for the full BTE.
    if (!calc_rta_only) {
      xmatrix<double> tcond_prev(3, 3), diff(3, 3);
      int num_iter = 1;
      double norm = 0.0;
      message << "Begin SCF for the Boltzmann transport equation.";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());
      ostream& oss = *_pc->getOSS();
      oss << std::setiosflags(std::ios::fixed | std::ios::right);
      oss << std::setw(15) << "Iteration";
      oss << std::setiosflags(std::ios::fixed | std::ios::right);
      oss << std::setw(25) << "Rel. Change in Norm" << std::endl;
      do {
        tcond_prev = tcond;
        getMeanFreeDispFull(rates_total, occ, mfd);

        tcond = calcTCOND(T, occ, mfd);
        // Calculate relative changes to the Frobenius norm instead of just
        // the norm. That way, less iterations for high thermal conductivity
        // materials are required.
        norm = frobenius_norm(tcond_prev - tcond)/frobenius_norm(tcond_prev);
        oss << std::setiosflags(std::ios::fixed | std::ios::right);
        oss << std::setw(15) << num_iter;
        oss << std::setiosflags(std::ios::fixed | std::ios::right);
        oss << std::setw(25) << std::dec << (norm) << std::endl;
        num_iter++;
      } while ((std::abs(norm) > TCOND_ITER_THRESHOLD) && (num_iter <= max_iter));
      if (num_iter > max_iter) {
        message << "Thermal conductivity did not converge within " << max_iter << " iterations ";
        message << "at " << T << " K.";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
    }
    return tcond;
  }

  //getOccupationNumbers//////////////////////////////////////////////////////
  // Calculates the Bose-Einstein distribution for all phonons at a specific
  // temperature.
  vector<vector<double> > TCONDCalculator::getOccupationNumbers(double T) {
    vector<vector<double> > occ(nQPs, vector<double>(nBranches, 0.0));
    for (int q = 0; q < nQPs; q++) {
      for (int br = 0; br < nBranches; br++) {
        occ[q][br] = 1.0/(exp(BEfactor_hbar_THz * freq[q][br]/T) - 1.0);
      }
    }
    return occ;
  }

  //getOccupationTerm/////////////////////////////////////////////////////////
  // Calculates the temperature-dependent prefactor for each scattering
  // process.
  double TCONDCalculator::getOccupationTerm(const vector<vector<double> >& occ, int sign,
      const vector<int>& qpts, const vector<int>& branches) {
    if (sign) {  // -
      return (1.0 + occ[qpts[1]][branches[1]] + occ[qpts[2]][branches[2]])/2.0;
    } else {  // +
      return occ[qpts[1]][branches[1]] - occ[qpts[2]][branches[2]];
    }
  }


  //calculateTotalRates///////////////////////////////////////////////////////
  // Calculates the total scattering rates based on three-phonon, isotope,
  // and boundary scattering.
  vector<vector<double> > TCONDCalculator::calculateTotalRates(const vector<vector<double> >& occ,
      vector<vector<double> >& rates_anharm) {
    vector<vector<double> > rates = calculateAnharmonicRates(occ);
    rates_anharm = rates;

    if (calc_isotope) {
      for (int iq = 0; iq < nIQPs; iq++) {
        for (int br = 0; br < nBranches; br++) {
          rates[iq][br] += scattering_rates_isotope[iq][br];
        }
      }
    }

    if (calc_boundary) {
      for (int iq = 0; iq < nIQPs; iq++) {
        for (int br = 0; br < nBranches; br++) {
          rates[iq][br] += scattering_rates_boundary[iq][br];
        }
      }
    }

    return rates;
  }

  //calculateAnharmonicRates//////////////////////////////////////////////////
  // Calculates the three-phonon scattering rates for a specific temperature.
  // Since there are a lot of processes, threading makes the calculation
  // significantly faster.
  vector<vector<double> > TCONDCalculator::calculateAnharmonicRates(const vector<vector<double> >& occ) {
    vector<vector<double> > rates(nIQPs, vector<double>(nBranches, 0.0));
#ifdef AFLOW_MULTITHREADS_ENABLE
    xthread::xThread xt(_pc->getNCPUs());
    std::function<void(int, const vector<vector<double> >&,
        vector<vector<double> >&)> fn = std::bind(&TCONDCalculator::calcAnharmRates, this, _1, _2, _3);
    xt.run(nIQPs, fn, occ, rates);
#else
    for (int i = 0; i < nIQPs; i++) calcAnharmRates(i, occ, rates);
#endif
    return rates;
  }

  void TCONDCalculator::calcAnharmRates(int i, const vector<vector<double> >& occ, vector<vector<double> >& rates) {
    vector<int> qpts(3), branches(3);
    int sign = -1;
    qpts[0] = _qm->getIbzqpts()[i];
    for (uint p = 0, nprocs = processes[i].size(); p < nprocs; p++) {
      getProcess(processes[i][p], qpts, branches, sign);
      rates[i][branches[0]] += intr_trans_probs[i][p] * getOccupationTerm(occ, sign, qpts, branches);
    }
  }

  //getMeanFreeDispRTA////////////////////////////////////////////////////////
  // Calculates the mean free displacement vector for the RTA.
  vector<vector<xvector<double> > > TCONDCalculator::getMeanFreeDispRTA(const vector<vector<double> >& rates) {
    xvector<double> xvec(3);
    vector<vector<xvector<double> > > mfd(nQPs, vector<xvector<double> >(nBranches, xvec));
    int iq = -1;
    for (int q = 0; q < nQPs; q++) {
      iq = _qm->getIbzqpt(q);
      for (int br = 0; br < nBranches; br++) {
        if (rates[iq][br] > 0.0) {
          mfd[q][br] = gvel[q][br] * freq[q][br]/rates[iq][br];
        }
      }
    }
    return mfd;
  }

  //calcTCOND/////////////////////////////////////////////////////////////////
  // Calculates the thermal conductivity tensor.
  xmatrix<double> TCONDCalculator::calcTCOND(double T, const vector<vector<double> >& occ,
      const vector<vector<xvector<double> > >& mfd) {
    xmatrix<double> tcond(3, 3);
    bool cumulative = calc_cumulative;
    double prefactor = 1E24 * std::pow(PLANCKSCONSTANT_hbar_THz, 2)/(KBOLTZ * std::pow(T, 2) * nQPs * Volume(_pc->getInputCellStructure()));
    for (int q = 0; q < nQPs; q++) {
      for (int br = 0; br < nBranches; br++) {
        bool include = true;
        if (freq[q][br] < _FLOAT_TOL_) {
          include = false;
        } else if (cumulative) {  // Only include processes below a certain mean free path if cumulative
          double mfpath = scalar_product(mfd[q][br], gvel[q][br])/aurostd::modulus(gvel[q][br]);
          include = !(mfpath > boundary_grain_size);
        }
        if (include) {
          double x = occ[q][br] * (occ[q][br] + 1) * freq[q][br];
          for (int i = 1; i < 4; i++) {
            for (int j = 1; j < 4; j++) {
              tcond[i][j] += x * gvel[q][br][i] * mfd[q][br][j];
            }
          }
        }
      }
    }
    tcond = prefactor * tcond;
    return tcond;
  }

  //getMeanFreeDispFull///////////////////////////////////////////////////////
  // Calculates/corrects the mean free displacement vector for the iterative
  // solution of the BTE. Since there are a lot of processes, threading
  // speeds up the calculations considerably.
  void TCONDCalculator::getMeanFreeDispFull(const vector<vector<double> >& rates,
      const vector<vector<double> >& occ,
      vector<vector<xvector<double> > >& mfd) {

    vector<vector<xvector<double> > > delta(nIQPs, vector<xvector<double> >(nBranches, xvector<double>(3)));
#ifdef AFLOW_MULTITHREADS_ENABLE
    xthread::xThread xt(_pc->getNCPUs());
    std::function<void(int, const vector<vector<double> >&,
        const vector<vector<xvector<double> > >&,
        vector<vector<xvector<double> > >&)> fn = std::bind(&TCONDCalculator::calculateDelta, this, _1, _2, _3, _4);
    xt.run(nIQPs, fn, occ, mfd, delta);
#else
    for (int i = 0; i < nIQPs; i++) calculateDelta(i, occ, mfd, delta);
#endif

    correctMFD(rates, delta, mfd);
  }

  //calculateDelta////////////////////////////////////////////////////////////
  // Calculates the correction vector (delta) to the mean free displacement.
  // Only irreducible q-points need to be calculated since deltas of
  // equivalent q-points are related by symmetry.
  void TCONDCalculator::calculateDelta(int i, const vector<vector<double> >& occ,
      const vector<vector<xvector<double> > >& mfd,
      vector<vector<xvector<double> > >& delta) {
    xvector<double> correction(3);
    vector<int> qpts(3), branches(3);
    int sign = -1;

    qpts[0] = _qm->getIbzqpts()[i];
    for (uint p = 0, nprocs = processes[i].size(); p < nprocs; p++) {
      getProcess(processes[i][p], qpts, branches, sign);
      correction = mfd[qpts[2]][branches[2]];
      if (processes[i][p][0] < 0) correction += mfd[qpts[1]][branches[1]];
      else correction -= mfd[qpts[1]][branches[1]];
      correction *= getOccupationTerm(occ, sign, qpts, branches);
      delta[i][branches[0]] += intr_trans_probs[i][p] * correction;
    }

    if (calc_isotope) {
      int q2 = 0, br1 = 0, br2 = 0;
      for (uint p = 0, nprocs = processes_iso[i].size(); p < nprocs; p++) {
        q2 = processes_iso[i][p][0];
        br1 = processes_iso[i][p][1]/nBranches;
        br2 = processes_iso[i][p][1] % nBranches;
        delta[i][br1] += intr_trans_probs_iso[i][p] * mfd[q2][br2];
      }
    }

    // Symmetrize
    int symop = 0;
    const vector<_sym_op>& pgroup = _qm->getReciprocalCell().pgroup;
    xmatrix<double> Uc(3, 3);
    const vector<int>& little_group = _qm->getLittleGroup(i);
    uint nsym = little_group.size();
    for (uint isym = 0; isym < nsym; isym++) {
      symop = little_group[isym];
      Uc += pgroup[symop].Uc;
    }
    Uc = 1.0/nsym * Uc;
    for (int br = 0; br < nBranches; br++) {
      delta[i][br] = Uc * delta[i][br];
    }
  }

  //correctMFD////////////////////////////////////////////////////////////////
  // Corrects the mean free displacement vectors of all phonon modes.
  void TCONDCalculator::correctMFD(const vector<vector<double> >& rates,
      const vector<vector<xvector<double> > >& delta,
      vector<vector<xvector<double> > >& mfd) {
    const vector<_sym_op>& pgroup = _qm->getReciprocalCell().pgroup;
    xvector<double> xvec(3);
    int iq = -1;
    for (int q = 0; q < nQPs; q++) {
      const xmatrix<double>& Uc = pgroup[_qm->getQPoint(q).symop].Uc;
      iq = _qm->getIbzqpt(q);
      for (int br = 0; br < nBranches; br++) {
        if (rates[iq][br] > 0.0) mfd[q][br] = (freq[q][br] * gvel[q][br] + Uc * delta[iq][br])/rates[iq][br];
        else mfd[q][br] = xvec;
      }
    }
  }

}  // namespace apl

/********************************* FILE I/O *********************************/

namespace apl {

  //writeTempIndepOutput//////////////////////////////////////////////////////
  // Writes temperature-independent output files.
  void TCONDCalculator::writeTempIndepOutput(const string& filename, string keyword,
      const string& unit, const vector<vector<double> >& data) {
    if (data.size() == 0) return;  // Nothing to write
    stringstream output;
    output << AFLOWIN_SEPARATION_LINE << std::endl;
    string key = "[AAPL_" + aurostd::toupper(aurostd::StringSubst(keyword, " ", "_")) + "]";
    if (!_pc->_system.empty()) {
      output << key << "SYSTEM=" << _pc->_system << std::endl;
    }
    output << key << "START" << std::endl;
    output << std::setiosflags(std::ios::fixed | std::ios::right);
    output << std::setw(10) << "# Q-point"
      << std::setw(20) << " "
      << keyword << " (" << unit << ")" << std::endl;
    writeDataBlock(output, data);
    output << key << "STOP" << std::endl;
    output << AFLOWIN_SEPARATION_LINE << std::endl;
    aurostd::stringstream2file(output, filename);
    if (!aurostd::FileExist(filename)) {
      stringstream message;
      message << "Could not write file " << filename << ".";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
  }

  //writeTempDepOutput////////////////////////////////////////////////////////
  // Writes temperature-dependent output files.
  void TCONDCalculator::writeTempDepOutput(const string& filename, string keyword, const string& unit,
      const vector<double>& temps, const vector<vector<vector<double> > >& data) {
    if (data.size() == 0) return; // Nothing to write
    stringstream output;
    output << AFLOWIN_SEPARATION_LINE << std::endl;
    string key = "[AAPL_" + aurostd::toupper(aurostd::StringSubst(keyword, " ", "_")) + "]";
    if (!_pc->_system.empty()) {
      output << key << "SYSTEM=" << _pc->_system << std::endl;
    }
    output << key << "START" << std::endl;
    for (uint t = 0; t < temps.size(); t++) {
      output << key << "T = " << std::fixed << std::setprecision(2) << temps[t] << " K" << std::endl;
      output << std::setw(10) << "# Q-point"
        << std::setw(20) << " "
        << keyword << " (" << unit << ")" << std::endl;
      writeDataBlock(output, data[t]);
    }
    output << key << "STOP" << std::endl;
    output << AFLOWIN_SEPARATION_LINE << std::endl;
    aurostd::stringstream2file(output, filename);
    if (!aurostd::FileExist(filename)) {
      stringstream message;
      message << "Could not write file " << filename << ".";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
  }

  //writeDataBlock////////////////////////////////////////////////////////////
  // Writes a block of data for each phonon mode into a stream. Used by
  // output file writers.
  void TCONDCalculator::writeDataBlock(stringstream& output,
      const vector<vector<double> >& data) {
    for (uint q = 0; q < data.size(); q++) {
      output << std::setiosflags(std::ios::fixed | std::ios::right);
      output << std::setw(10) << q;
      for (uint br = 0; br < data[q].size(); br++) {
        output << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
        output << std::setw(20) << std::setprecision(10) << std::scientific << data[q][br];
      }
      output << std::endl;
    }
  }

  //writePhaseSpace///////////////////////////////////////////////////////////
  // Writes the scattering phase space into an output file. Numbers are
  // converted into fs to get "nicer" numbers. OOO should be zero, but it is
  // output regardless in case there is a problem.
  void TCONDCalculator::writePhaseSpace(const string& filename) {
    if (phase_space.size() == 0) return;  // Nothing to write
    // Calculate totals
    vector<double> ps_procs(4), ps_nu(2);
    vector<vector<double> > ps_modes(nIQPs, vector<double>(nBranches));
    double ps_total = 0.0, ps = 0.0;
    for (int i = 0; i < nIQPs; i++) {
      for (int b = 0; b < nBranches; b++) {
        for (int p = 0; p < 4; p++) {
          for (int s = 0; s < 2; s++) {
            ps = 1000 * phase_space[i][b][p][s];  // Convert to fs
            ps_total += ps;
            ps_nu[s] += ps;
            ps_procs[p] += ps;
            ps_modes[i][b] += ps;
          }
        }
      }
    }

    stringstream output;
    output << AFLOWIN_SEPARATION_LINE << std::endl;
    output << "# 3-phonon scattering phase space (in fs)" << std::endl;
    output << "[AAPL_SCATTERING_PHASE_SPACE]SYSTEM=" << _pc->_system << std::endl;
    output << "[AAPL_TOTAL_SCATTERING_PHASE_SPACE]START" << std::endl;
    output << std::setiosflags(std::ios::left | std::ios::fixed | std::ios::showpoint);
    output << std::setw(15) << "total" << std::setw(10) << std::setprecision(8) << std::dec << ps_total << std::endl;
    output << std::setw(15) << "total normal" << std::setw(10) << std::setprecision(8) << std::dec << ps_nu[0] << std::endl;
    output << std::setw(15) << "total umklapp" << std::setw(10) << std::setprecision(8) << std::dec << ps_nu[1] << std::endl;
    output << std::setw(15) << "total AAA" << std::setw(10) << std::setprecision(8) << std::dec << ps_procs[0] << std::endl;
    output << std::setw(15) << "total AAO" << std::setw(10) << std::setprecision(8) << std::dec << ps_procs[1] << std::endl;
    output << std::setw(15) << "total AOO" << std::setw(10) << std::setprecision(8) << std::dec << ps_procs[2] << std::endl;
    output << std::setw(15) << "total OOO" << std::setw(10) << std::setprecision(8) << std::dec << ps_procs[3] << std::endl;
    output << "[AAPL_TOTAL_SCATTERING_PHASE_SPACE]STOP" << std::endl;
    output << "[AAPL_SCATTERING_PHASE_SPACE]START" << std::endl;
    for (int i = 0; i < nIQPs; i++) {
      for (int b = 0; b < nBranches; b++) {
        output << std::setw(17) << std::setprecision(10) << std::dec << ps_modes[i][b];
      }
      output << std::endl;
    }
    output << "[AAPL_SCATTERING_PHASE_SPACE]STOP" << std::endl;
    output << AFLOWIN_SEPARATION_LINE << std::endl;
    aurostd::stringstream2file(output, filename);
  }

  //writeGrueneisen///////////////////////////////////////////////////////////
  // Outputs the temperature-dependent average Grueneisen parameters and the
  // mode Grueneisen parameters into a file.
  void TCONDCalculator::writeGrueneisen(const string& filename) {
    if (grueneisen_mode.size() == 0) return;  // Nothing to write
    stringstream output;
    output << AFLOWIN_SEPARATION_LINE << std::endl;
    output << "[AAPL_GRUENEISEN]SYSTEM=" << _pc->_system << std::endl;
    output << "[AAPL_GRUENEISEN_AVERAGE]START" << std::endl;
    output << std::setiosflags(std::ios::right | std::ios::fixed | std::ios::showpoint);
    output << std::setw(8) << "# T (K)"
      << std::setw(23) << "Grueneisen parameter" << std::endl;
    for (uint t = 0; t < grueneisen_avg.size(); t++) {
      output << std::setw(8) << std::fixed << std::setprecision(2) << temperatures[t];
      output << std::setw(23) << std::setprecision(10) << std::dec << grueneisen_avg[t] << std::endl;
    }
    output << "[AAPL_GRUENEISEN_AVERAGE]STOP" << std::endl;
    output << "[AAPL_GRUENEISEN_MODE]START" << std::endl;
    for (int i = 0; i < nIQPs; i++) {
      for (int b = 0; b < nBranches; b++) {
        output << std::setw(17) << std::setprecision(10) << std::dec << grueneisen_mode[i][b];
      }
      output << std::endl;
    }
    output << "[AAPL_GRUENEISEN_MODE]STOP" << std::endl;
    output << AFLOWIN_SEPARATION_LINE << std::endl;
    aurostd::stringstream2file(output, filename);
  }

  //writeThermalConductivity////////////////////////////////////////////////////
  // Outputs the thermal conductivity tensor into a file.
  void TCONDCalculator::writeThermalConductivity(const string& filename) {
    if (thermal_conductivity.size() == 0) return;  // Nothing to write
    stringstream output;

    // Header
    output << AFLOWIN_SEPARATION_LINE << std::endl;
    if (!_pc->_system.empty()) {
      output << "[AAPL_THERMAL_CONDUCTIVITY]SYSTEM=" << _pc->_system << std::endl;
    }
    output << "[AAPL_THERMAL_CONDUCTIVITY]START" << std::endl;
    output << std::setw(8) << "# T (K)"
      << std::setw(75) << " "
      << "Thermal Conductivity (W/m*K)" << std::endl;
    string xyz[3] = {"x", "y", "z"};
    output << "#" << std::setw(8) << " ";
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        output << std::setw(7) << " ";
        // Write columns first to make compatible with gnuplot
        string k = "k(" + xyz[j] + "," + xyz[i] + ")";
        output << k;
        output << std::setw(7) << " ";
      }
    }
    output << std::endl;

    // Body
    for (uint t = 0; t < temperatures.size(); t++) {
      output << std::setw(8) << std::fixed << std::setprecision(2) << temperatures[t];
      for (int i = 1; i < 4; i++) {
        for (int j = 1; j < 4; j++) {
          output << std::setw(2) << " ";
          // Write columns first to make compatible with gnuplot
          output << std::setprecision(10) << std::scientific << thermal_conductivity[t][j][i];
          output << std::setw(2) << " ";
        }
      }
      output << std::endl;
    }
    output << "[AAPL_THERMAL_CONDUCTIVITY]STOP" << std::endl;
    output << AFLOWIN_SEPARATION_LINE << std::endl;

    aurostd::stringstream2file(output, filename);
    if (!aurostd::FileExist(filename)) {
      string message = "Could not write thermal conductivities to file.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
  }

}  // namespace apl

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *            Aflow MARCO ESTERS - Duke University 2019-2021               *
// *                                                                         *
// ***************************************************************************
