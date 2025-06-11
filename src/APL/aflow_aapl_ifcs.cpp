// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *            Aflow MARCO ESTERS - Duke University 2017-2021               *
// *                                                                         *
// ***************************************************************************
// Written by Marco Esters, 2018. Based on work by Jose J. Plata (AFLOW AAPL,
// DOI: 10.1038/s41524-017-0046-7) and Jesus Carrete (ShengBTE, 
// DOI: 10.1016/j.cpc.2014.02.015).
//
// This class calculates the anharmonic interatomic force constants (IFCs). It
// reads the forces of the inequivalent distortions first, and then calculates
// the forces of the equivalent distortions. These forces are used to
// calculate the IFCs of the inequivalent clusters that are then symmetrized
// in an iterative procedure using sum rules.
//
// See aflow_apl.h for descriptions of the classes and their members, and for
// the struct _linearCombinations.

#include "aflow_apl.h"

#define _DEBUG_AAPL_IFCS_ false

using std::vector;
using aurostd::xcombos;
using aurostd::xerror;

static const string _CLUSTER_SET_FILE_[2] = {"clusterSet_3rd.xml", "clusterSet_4th.xml"};

// tform represents a tensor transformation containing the index and the
// coefficients. vector<vector<int> > holds the indices, vector<double>
// the coefficients.
typedef vector<std::pair<vector<int>, vector<double> > > tform;
// v4int defined for brevity
typedef vector<vector<vector<vector<int> > > > v4int;

// Coefficients for the finite difference method
static const double C12[2][3] = {{0.5,  0.0, -0.5},
  {1.0, -2.0,  1.0}};
static const double C3[5] = {-0.5, 1.0, 0.0, -1.0, 0.5};

/************************ CONSTRUCTORS/DESTRUCTOR ***************************/

namespace apl {

  //Constructors//////////////////////////////////////////////////////////////
  // Default constructor
  AnharmonicIFCs::AnharmonicIFCs(ostream& oss) : xStream(oss) {
    free();
    clst = ClusterSet(oss);
    directory = "./";
  }

  AnharmonicIFCs::AnharmonicIFCs(ofstream& mf, ostream& oss) : xStream(mf,oss) {
    free();
    clst = ClusterSet(mf, oss);
    directory = "./";
  }

  //Copy constructors
  AnharmonicIFCs::AnharmonicIFCs(const AnharmonicIFCs& that) : xStream(*that.getOFStream(),*that.getOSS()) {
    if (this != &that) free();
    copy(that);
  }

  const AnharmonicIFCs& AnharmonicIFCs::operator=(const AnharmonicIFCs& that) {
    if (this != &that) free();
    copy(that);
    return *this;
  }

  //copy//////////////////////////////////////////////////////////////////////
  void AnharmonicIFCs::copy(const AnharmonicIFCs& that) {
    if (this == &that) return;
    xStream::copy(that);
    cart_indices = that.cart_indices;
    clst = that.clst;
    directory = that.directory;
    distortion_magnitude = that.distortion_magnitude;
    force_constants = that.force_constants;
    initialized = that.initialized;
    max_iter = that.max_iter;
    mixing_coefficient = that.mixing_coefficient;
    order = that.order;
    sumrule_threshold = that.sumrule_threshold;
    _useZeroStateForces = that._useZeroStateForces;
    xInputs = that.xInputs;
  }

  //Destructor////////////////////////////////////////////////////////////////
  AnharmonicIFCs::~AnharmonicIFCs() {
    xStream::free();
    free();
  }

  //free//////////////////////////////////////////////////////////////////////
  // Clears all vectors and resets all values.
  void AnharmonicIFCs::free() {
    cart_indices.clear();
    clst.clear();
    directory = "";
    distortion_magnitude = 0.0;
    force_constants.clear();
    initialized = false;
    max_iter = 0;
    mixing_coefficient = 0.0;
    order = 0;
    sumrule_threshold = 0.0;
    _useZeroStateForces = false;
    xInputs.clear();
  }

  //clear/////////////////////////////////////////////////////////////////////
  void AnharmonicIFCs::clear() {
    free();
  }


  //setOptions////////////////////////////////////////////////////////////////
  void AnharmonicIFCs::setOptions(double dmag, int iter, double mix, double threshold, bool zero) {
    distortion_magnitude = dmag;
    max_iter = iter;
    mixing_coefficient = mix;
    sumrule_threshold = threshold;
    _useZeroStateForces = zero;
  }

  //directory/////////////////////////////////////////////////////////////////
  const string& AnharmonicIFCs::getDirectory() const {
    return directory;
  }

  void AnharmonicIFCs::setDirectory(const string& dir) {
    directory = dir;
    clst.directory = dir;
  }

  //initialize////////////////////////////////////////////////////////////////
  // Initializes the anharmonic IFC calculator by building the ClusterSet.
  void AnharmonicIFCs::initialize(const Supercell& scell, int _order, const aurostd::xoption& opts, ofstream& mf, ostream& oss) {
    xStream::initialize(mf, oss);
    initialize(scell, _order, opts);
  }

  void AnharmonicIFCs::initialize(const Supercell& scell, int _order, const aurostd::xoption& opts) {
    string message = "";
    // Initialize IFC parameters
    order = _order;
    distortion_magnitude = aurostd::string2utype<double>(opts.getattachedscheme("DMAG"));
    max_iter = aurostd::string2utype<int>(opts.getattachedscheme("SUMRULE_MAX_ITER"));
    mixing_coefficient = aurostd::string2utype<double>(opts.getattachedscheme("MIXING_COEFFICIENT"));
    sumrule_threshold = aurostd::string2utype<double>(opts.getattachedscheme("SUMRULE"));
    _useZeroStateForces = opts.flag("ZEROSTATE");

    // Initialize cluster
    vector<string> tokens;
    aurostd::string2tokens(opts.getattachedscheme("CUT_RAD"), tokens, ",");
    if (tokens.size() < (uint) order - 2) {
      message = "Not enough parameters for CUT_RAD";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    double cut_rad = aurostd::string2utype<double>(tokens[order - 3]);
    aurostd::string2tokens(opts.getattachedscheme("CUT_SHELL"), tokens, ",");
    if (tokens.size() < (uint) order - 2) {
      message = "Not enough parameters for CUT_SHELL";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    int cut_shell = aurostd::string2utype<int>(tokens[order - 3]);
    clst.initialize(scell, _order, cut_shell, cut_rad);
    string clust_hib_file = directory + "/" + DEFAULT_AAPL_FILE_PREFIX + _CLUSTER_SET_FILE_[_order-3];
    bool awakeClusterSet = aurostd::EFileExist(clust_hib_file);
    if (awakeClusterSet) {
      try {
        clst.readClusterSetFromFile(clust_hib_file);
      } catch (aurostd::xerror& excpt) {
        awakeClusterSet = false;
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, excpt.buildMessageString(), directory, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
      }
    }
    if (!awakeClusterSet) {
      clst.build();
      clst.buildDistortions();
      clst.writeClusterSetToFile(clust_hib_file);
    }
    initialized = true;
  }

  //xStream initializers
  void AnharmonicIFCs::initialize(ostream& oss) {
    xStream::initialize(oss);
  }

  void AnharmonicIFCs::initialize(ofstream& mf, ostream& oss) {
    xStream::initialize(mf, oss);
  }

  //getOrder//////////////////////////////////////////////////////////////////
  int AnharmonicIFCs::getOrder() const {
    return order;
  }

}  // namespace apl

/***************************** DFT CALCULATIONS *****************************/

namespace apl {

  bool AnharmonicIFCs::runVASPCalculations(_xinput& xinput, _aflags& aflags, _kflags& kflags, _xflags& xflags) {
    string message = "";
    if (order > 4) {
      message = "Not implemented for order > 4.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }
    if (!initialized) {
      message = "Not initialized.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }

    bool stagebreak = false;
    xinput.xvasp.AVASP_arun_mode = __AFLOW_FUNC__;
    stringstream _logger;
    _logger << "Managing directories for ";
    if (order == 3) {
      _logger << "3rd";
    } else if (order == 4) {
      _logger << "4th";
    }
    _logger << " order IFCs.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, _logger, directory, *p_FileMESSAGE, *p_oss);

    // Determine the number of runs so the run ID in the folder name can be
    // padded with the appropriate number of zeros.
    int nruns = 0;
    for (uint ineq = 0; ineq < clst.ineq_distortions.size(); ineq++) {
      nruns += clst.ineq_distortions[ineq].distortions.size();
    }
    if (order == 4) {
      for (uint ineq = 0; ineq < clst.higher_order_ineq_distortions.size(); ineq++) {
        nruns += clst.higher_order_ineq_distortions[ineq].distortions.size();
      }
    }

    xInputs.assign(nruns, xinput);

    int idxRun = 0;
    for (uint ineq = 0; ineq < clst.ineq_distortions.size(); ineq++) {
      const vector<int>& atoms = clst.ineq_distortions[ineq].atoms;
      const _ineq_distortions& idist = clst.ineq_distortions[ineq];
      for (uint dist = 0; dist < idist.distortions.size(); dist++) {
        const vector<int>& distortions = idist.distortions[dist][0];

        //ME20190109 - add title
        xstructure& xstr = xInputs[idxRun].getXStr();
        LightCopy(clst.scell, xstr);
        xstr.title = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xstr.title);
        if (xstr.title.empty()) {
          xstr.buildGenericTitle(true, false);
        }
        xstr.title += " AAPL supercell=" + aurostd::joinWDelimiter(clst.sc_dim, 'x');

        // Set up runname and generate distorted structure
        xInputs[idxRun].xvasp.AVASP_arun_runname = buildRunName(distortions, atoms, idxRun, nruns);
        applyDistortions(xInputs[idxRun], clst.distortion_vectors, distortions, atoms);

        // Create aflow.in files if they don't exist. Stagebreak is true as soon
        // as one aflow.in file was created.
        stagebreak = (createAflowInPhonons(aflags, kflags, xflags, xInputs[idxRun]) || stagebreak);
        idxRun++;
      }
    }
    if (order == 4) {
      for (uint ineq = 0; ineq < clst.higher_order_ineq_distortions.size(); ineq++) {
        const _ineq_distortions& idist = clst.higher_order_ineq_distortions[ineq];
        const vector<int>& atoms = idist.atoms;
        for (uint dist = 0; dist < idist.distortions.size(); dist++) {
          xInputs[idxRun] = xinput;
          const vector<int>& distortions = idist.distortions[dist][0];
          xstructure& xstr = xInputs[idxRun].getXStr();
          xstr.title = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xstr.title);

          if (xstr.title.empty()) {
            xstr.buildGenericTitle(true, false);
          }
          xstr.title += " AAPL supercell=" + aurostd::joinWDelimiter(clst.sc_dim, 'x');

          //ME20190113 - make sure that POSCAR has the correct format
          if ((kflags.KBIN_MPI && (kflags.KBIN_BIN.find("46") != string::npos)) ||
              (kflags.KBIN_MPI && (kflags.KBIN_MPI_BIN.find("46") != string::npos))) {
            xstr.is_vasp5_poscar_format = false;
          }

          // Set up runname and generate distorted structure
          xInputs[idxRun].xvasp.AVASP_arun_runname = buildRunName(distortions, atoms, idxRun, nruns);
          xInputs[idxRun].xvasp.AVASP_arun_runname += "_4th";
          applyDistortions(xInputs[idxRun], clst.distortion_vectors, distortions, atoms, 2.0);

          // Create aflow.in files if they don't exist. Stagebreak is true as soon
          // as one aflow.in file was created.
          stagebreak = (createAflowInPhonons(aflags, kflags, xflags, xInputs[idxRun]) || stagebreak);
          idxRun++;
        }
      }
    }

    return stagebreak;
  }

  string AnharmonicIFCs::buildRunName(const vector<int>& distortions,
      const vector<int>& atoms, int run, int nruns) {
    stringstream runname;
    runname << order << "_";

    // Run ID with padding
    runname << aurostd::PaddedNumString(run + 1, aurostd::getZeroPadding(nruns)) << "_";

    // Atom and distortion IDs
    for (int at = 0; at < order - 1; at++) {
      if (at > 0) runname << "-";
      runname << "A" << atoms[at] << "D" << distortions[at];
    }

    return runname.str();
  }

  void AnharmonicIFCs::applyDistortions(_xinput& xinp,
      const vector<xvector<double> >& distortion_vectors,
      const vector<int>& distortions,
      const vector<int>& atoms, double scale) {
    xstructure& xstr = xinp.getXStr();  //ME20190109
    for (uint at = 0; at < atoms.size(); at++) {
      int atsc = atoms[at];
      int dist_index = distortions[at];
      xvector<double> dist_cart = distortion_vectors[dist_index];
      while (((at + 1) < atoms.size()) && (atoms[at] == atoms[at+1])) {
        at++;
        dist_index = distortions[at];
        dist_cart += distortion_vectors[dist_index];
      }
      // Normalize dist_cart coordinates to 1 so that distortions have the same magnitude
      for (int i = 1; i < 4; i++) {
        if (abs(dist_cart(i)) < _ZERO_TOL_) {
          dist_cart(i) = 0.0;
        } else {
          dist_cart(i) *= (scale/std::abs(dist_cart(i)));
        }
      }
      dist_cart *= distortion_magnitude;
      //ME20190109 - Add to title
      xstr.title += " atom=" + aurostd::utype2string<int>(atoms[at]);
      std::stringstream distortion; // ME20190112 - need stringstream for nicer formatting
      distortion << " distortion=["
        << std::setprecision(3) << dist_cart[1] << ","
        << std::setprecision(3) << dist_cart[2] << ","
        << std::setprecision(3) << dist_cart[3] << "]"; //ME20190112
      xstr.title += distortion.str();
      xstr.atoms[atsc].cpos += dist_cart;
      xstr.atoms[atsc].fpos = xstr.c2f * xstr.atoms[atsc].cpos;
    }
  }


}  // namespace apl

/****************************** CALCULATE IFCs ******************************/

namespace apl {

  bool AnharmonicIFCs::calculateForceConstants() {
    cart_indices = getCartesianIndices();

    // Read forces
    // outfileFoundAnywherePhonons detects whether no calculations have been
    // run (no error message) whereas outfileFounEverywherePhonons tries to
    // read the force files. The latter outputs messages, which is not desired
    // when directories have just been created.
    if (!outfileFoundAnywherePhonons(xInputs)) return false;
    if (!outfileFoundEverywherePhonons(xInputs, directory, *p_FileMESSAGE, *p_oss)) return false;
    if (_useZeroStateForces) {
      vector<string> dirs;
      aurostd::DirectoryLS(directory, dirs);
      uint ndir = dirs.size();
      uint d = 0;
      for (d = 0; d < ndir; d++) {
        if (aurostd::IsDirectory(directory + "/" + dirs[d]) && aurostd::substring2bool(dirs[d], "ZEROSTATE")) {
          break;
        }
      }
      if (d == ndir) {
        string message = "Could not find ZEROSTATE directory.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
      } else {
        _xinput zerostate(xInputs[0]);
        xstructure& xstr = zerostate.getXStr();
        LightCopy(clst.scell, xstr);
        zerostate.setDirectory(aurostd::CleanFileName(directory + "/" + dirs[d]));
        subtractZeroStateForces(xInputs, zerostate);
      }
    }
    vector<vector<vector<xvector<double> > > > force_tensors = storeForces(xInputs);

    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Calculating anharmonic IFCs.", directory, *p_FileMESSAGE, *p_oss);
    vector<vector<double> > ifcs_unsym = calculateUnsymmetrizedIFCs(clst.ineq_distortions, force_tensors);
    force_tensors.clear();

    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Symmetrizing IFCs.", directory, *p_FileMESSAGE, *p_oss);
    force_constants = symmetrizeIFCs(ifcs_unsym);
    return true;
  }

  const vector<vector<double> >& AnharmonicIFCs::getForceConstants() const {
    return force_constants;
  }

  vector<vector<int> > AnharmonicIFCs::getClusters() const {
    uint nclusters = clst.clusters.size();
    vector<vector<int> > clusters(nclusters, vector<int>(order));
    for (uint i = 0; i < nclusters; i++) {
      for (int j = 0; j < order; j++) {
        clusters[i][j] = clst.clusters[i].atoms[j];
      }
    }
    return clusters;
  }

}  // namespace apl

/*************************** INITIAL CALCULATIONS ***************************/

namespace apl {

  //getCartesianIndices///////////////////////////////////////////////////////
  // Returns a list of Cartesian indices. Since they are looped over
  // frequently, it is quicker to calculate them once at the beginning.
  vector<vector<int> > AnharmonicIFCs::getCartesianIndices() {
    vector<vector<int> > indices;
    xcombos crt_ind(3, order, 'E', true);
    while (crt_ind.increment()) {
      indices.push_back(crt_ind.getCombo());
    }
    return indices;
  }

  // BEGIN Forces
  //storeForces///////////////////////////////////////////////////////////////
  // Stores the forces from the VASP calculations. Each item in the vector
  // holds the force tensor for a set of distorted atoms.
  vector<vector<vector<xvector<double> > > > AnharmonicIFCs::storeForces(vector<_xinput>& xInp) {
    vector<vector<vector<xvector<double> > > > force_tensors(clst.ineq_distortions.size());
    int idxRun = 0;
    for (uint id = 0; id < clst.ineq_distortions.size(); id++) {
      force_tensors[id] = getForces(id, idxRun, xInp);
    }
    if (order == 4) addHigherOrderForces(force_tensors, idxRun, xInp);
    return force_tensors;
  }

  //getForces/////////////////////////////////////////////////////////////////
  // Retrieves all forces from the calculations. Also transforms the forces
  // into the forces of the equivalent distortions.
  vector<vector<xvector<double> > > AnharmonicIFCs::getForces(int id, int& idxRun,
      vector<_xinput>& xInp) {
    const _ineq_distortions& ineq_dists = clst.ineq_distortions[id];
    int natoms = (int) clst.scell.atoms.size();
    vector<int> powers(order - 1, 1);
    int ndist = 1;
    for (int i = 0; i < order - 1; i++) {
      ndist *= 6;
      for (int j = 0; j < order - 2 - i; j++) {
        powers[i] *= 6;
      }
    }
    vector<vector<xvector<double> > > force_tensor(ndist, vector<xvector<double> >(natoms, xvector<double>(3)));

    int attrans = 0, fg = 0, index = 0;
    for (uint ineq = 0; ineq < ineq_dists.distortions.size(); ineq++) {
      // For the inequivalent distortion, just read the forces from VASP 
      const vector<xvector<double> >& qmforces = xInp[idxRun].getXStr().qm_forces;
      index = 0;
      for (int ind = 0; ind < order - 1; ind++) {
        index += powers[ind] * ineq_dists.distortions[ineq][0][ind];
      }
      for (int at = 0; at < natoms; at++) { 
        force_tensor[index][at] = qmforces[at];
      }
      for (uint i = 1; i < ineq_dists.distortions[ineq].size(); i++) {
        // For each equivalent distortion, transform the forces using symmetry
        vector<xvector<double> > qmforces_trans(natoms, xvector<double>(3));
        index = 0;
        for (int ind = 0; ind < order - 1; ind++) {
          index += powers[ind] * ineq_dists.distortions[ineq][i][ind];
        }
        fg = ineq_dists.rotations[ineq][i];
        const vector<int>& transformation_map = ineq_dists.transformation_maps[ineq][i];
        for (int at = 0; at < natoms; at++) {
          attrans = getTransformedAtom(transformation_map, at);
          force_tensor[index][at] = clst.pcell.fgroup[fg].Uc * qmforces[attrans];
        }
      }
      idxRun++;
    }
    return force_tensor;
  }

  void AnharmonicIFCs::addHigherOrderForces(vector<vector<vector<xvector<double> > > >& force_tensor,
      int& idxRun, vector<_xinput>& xInp) {
    const vector<_ineq_distortions>& ineq_dists = clst.higher_order_ineq_distortions;
    uint ndist = force_tensor[0].size();
    uint natoms = clst.scell.atoms.size();
    vector<xvector<double> > forces(natoms, xvector<double>(3));
    for (uint ineq = 0; ineq < ineq_dists.size(); ineq++) {
      uint idist;
      int at = ineq_dists[ineq].atoms[0];
      for (idist = 0; idist < clst.ineq_distortions.size(); idist++) {
        int a;
        for (a = 0; a < order - 1; a++) {
          if (clst.ineq_distortions[idist].atoms[a] != at) break;
        }
        if (a == order - 1) break;
      }
      for (int i = 0; i < 6; i++) force_tensor[idist].push_back(forces);

      for (uint dist = 0; dist < ineq_dists[ineq].distortions.size(); dist++) {
        const vector<xvector<double> >& qmforces = xInp[idxRun].getXStr().qm_forces;
        int d = ineq_dists[ineq].distortions[dist][0][0];
        force_tensor[idist][ndist + d][at] = qmforces[at];
        int fg = 0;
        for (uint i = 1; i < ineq_dists[ineq].distortions[dist].size(); i++) {
          d = ineq_dists[ineq].distortions[dist][i][0];
          fg = ineq_dists[ineq].rotations[dist][i];
          force_tensor[idist][ndist + d][at] = clst.pcell.fgroup[fg].Uc * qmforces[at];
        }
        xInp[idxRun].clear();
        idxRun++;
      }
    }
  }

  //getTransformedAtom////////////////////////////////////////////////////////
  // Retrieves the atom that is obtained by the transformation in the given
  // symmetry map.
  int AnharmonicIFCs::getTransformedAtom(const vector<int>& symmap, const int& at) {
    for (uint i = 0; i < symmap.size(); i++) {
      if (symmap[i] == at) {
        return i;
      }
    }
    // If the for-loop runs until the end, the atom was not found
    stringstream message;
    message << "Could not transform atom " << at;
    throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
  }
  //END Forces

  // BEGIN Calculate unsymmetrized force constants
  //calculateUnsymmetrizedIFCs////////////////////////////////////////////////
  // Calculates the IFCs of the inequivalent clusters from the forces.
  vector<vector<double> >
    AnharmonicIFCs::calculateUnsymmetrizedIFCs(const vector<_ineq_distortions>& idist,
        const vector<vector<vector<xvector<double> > > >& forces) {
      vector<vector<double> > ifcs(clst.clusters.size(), vector<double>(cart_indices.size()));
      int at = 0, cl = 0, ic = 0;
      double denom = std::pow(distortion_magnitude, order - 1);
      for (uint f = 0; f < forces.size(); f++) {
        for (uint c = 0; c < idist[f].clusters.size(); c++) {
          ic = idist[f].clusters[c];
          cl = clst.ineq_clusters[ic][0];
          at = clst.getCluster(cl).atoms[order - 1];
          for (int cart = 0; cart < clst.nifcs; cart++) {
            ifcs[cl][cart] = finiteDifference(forces[f], at, cart_indices[cart], idist[f].atoms)/denom;
          }
        }
      }
      return ifcs;
    }

  //finiteDifference//////////////////////////////////////////////////////////
  // Calculates the numerator for a specific IFC from the forces using the
  // central difference method. The denominator is calculated separately.
  //
  // See https://www.geometrictools.com/Documentation/FiniteDifferences.pdf
  // for formulas of higher order derivatives.
  double AnharmonicIFCs::finiteDifference(const vector<vector<xvector<double> > >& forces, int at,
      const vector<int>& cart_ind, const vector<int>& atoms) {
    vector<int> powers(order - 1, 1);
    for (int i = 0; i < order - 2; i++) {
      for (int j = 0; j < order - 2 - i; j++) {
        powers[i] *= 6;
      }
    }

    vector<int> derivatives;
    int count = 0;
    for (uint a = 0; a < atoms.size(); a++) {
      count = 1;
      while (((int) a + count < order - 1) &&
          (atoms[a + count] == atoms[a]) &&
          (cart_ind[a + count] == cart_ind[a])) {
        count++;
      }
      derivatives.push_back(count);
      a += count - 1;
    }

    double diff = 0.0;
    vector<int> atm(order - 1);
    if (derivatives[0] == 3) {
      vector<int> dists(5);
      int pwr = 0;
      for (int i = 0; i < order - 1; i++) pwr += powers[i];
      dists[0] = clst.nifcs + cart_ind[0];
      dists[1] = cart_ind[0] * pwr;
      // No need to occupy dists[2] because C3[2] is zero
      dists[3] = (cart_ind[0] + 3) * pwr;
      dists[4] = clst.nifcs + cart_ind[0] + 3;
      for (int i = 0; i < 5; i++) {
        if (C3[i] != 0.0) diff -= C3[i] * forces[dists[i]][at][cart_ind[order - 1] + 1];
      }
    } else {
      uint nder = derivatives.size();
      vector<int> index(nder);
      xcombos ind(3, nder, 'E', true);
      double coeff = 0.0;
      int d = 0, dist = 0;
      while (ind.increment()) {
        coeff = 1.0;
        index = ind.getCombo();
        for (uint i = 0; i < nder; i++) {
          coeff *= C12[derivatives[i] - 1][index[i]];
        }
        if (coeff != 0.0) {
          dist = 0;
          d = 0;
          for (uint i = 0; i < nder; i++) {
            for (int j = 0; j < derivatives[i]; j++) {
              if (index[i] == 1) dist += powers[d] * (cart_ind[d] + 3 * j);
              else dist += powers[d] * (cart_ind[d] + 3 * index[i]/2);
              d++;
            }
          }
          diff -= coeff * forces[dist][at][cart_ind[order - 1] + 1];
        }
      }
    }
    return diff;
  }
  // END Calculate unsymmetrized force constants

}  // namespace apl

/****************************** SYMMETRIZATION ******************************/

namespace apl {

  //symmetrizeIFCs////////////////////////////////////////////////////////////
  // Symmetrizes the IFCs using an iterative procedure to ensure 
  // that the force constants sum to zero.
  // 1. The IFCs of the inequivalent clusters will be symmetrized according to
  //    the linear combinations found while determining ClusterSets.
  // 2. The force constants will be transformed to the other clusters using
  //    the symmetry of the crystal.
  // 3. Determine the deviations of the IFC sums from zero.
  // 4. If at least one deviation is above the chosen threshold, correct the
  //    linearly dependent IFCs. If none are above the threshold or too many
  //    iterations have been performed, stop the iteration procedure.
  //
  // Check the typedefs at the beginning of the file for tform and v4int
  vector<vector<double> > AnharmonicIFCs::symmetrizeIFCs(vector<vector<double> > ifcs) {
    // Initialize tensors
    vector<vector<int> > reduced_clusters = getReducedClusters();
    vector<vector<double> > dev_from_zero(reduced_clusters.size(), vector<double>(cart_indices.size()));
    vector<vector<double> > abssum = dev_from_zero;

    // Tensor transformations
    vector<vector<tform> > transformations(clst.ineq_clusters.size());
    v4int eq_ifcs(clst.ineq_clusters.size());
    getTensorTransformations(eq_ifcs, transformations);

    // Do iterations
    int num_iter = 0;
    double max_err = 0.0;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Begin SCF for anharmonic force constants.", directory, *p_FileMESSAGE, *p_oss);
    *p_oss << std::setiosflags(std::ios::fixed | std::ios::right);
    *p_oss << std::setw(15) << "Iteration";
    *p_oss << std::setiosflags(std::ios::fixed | std::ios::right);
    *p_oss << std::setw(20) << "Abs. max. error" << std::endl;
    do {
      // 1. Symmetrize using linear combinations
      applyLinCombs(ifcs);  

      // 2. Transform IFCs using symmetry and permutations
      transformIFCs(transformations, ifcs);

      // 3. Determine deviations from zero
      calcSums(reduced_clusters, ifcs, dev_from_zero, abssum);

      max_err = -1E30;
      for (uint i = 0; i < dev_from_zero.size(); i++) {
        for (uint j = 0; j < dev_from_zero[i].size(); j++) {
          if (std::abs(dev_from_zero[i][j]) > max_err) max_err = std::abs(dev_from_zero[i][j]);
        }
      }

      *p_oss << std::setiosflags(std::ios::fixed | std::ios::right);
      *p_oss << std::setw(15) << num_iter;
      *p_oss << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
      *p_oss << std::setw(20) << max_err << std::endl;

      // 4. Correct IFCs
      if (max_err > sumrule_threshold) {
        correctIFCs(ifcs, dev_from_zero, abssum, reduced_clusters, eq_ifcs);
      }
      num_iter++;
    } while ((num_iter <= max_iter) && (max_err > sumrule_threshold));
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "End SCF for anharmonic force constants.", directory, *p_FileMESSAGE, *p_oss);
    if (num_iter > max_iter) {
      stringstream message;
      message << "Anharmonic force constants did not converge within " << max_iter << " iterations.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    } else {
      return ifcs;
    }
  }

  // BEGIN Initializers
  //getReducedClusters////////////////////////////////////////////////////////
  // Determines, for each set of inequivalent clusters, a uinque set of
  // clusters that do not contain the last atom of the clusters. This set is
  // important for the sum rules as they frequently require summations over
  // the last atom within a set of inequivalent clusters.
  vector<vector<int> > AnharmonicIFCs::getReducedClusters() {
    vector<vector<int> > reduced_clusters;
    vector<int> cluster(order - 1);
    for (uint c = 0; c < clst.clusters.size(); c++) {
      const vector<int>& atoms = clst.clusters[c].atoms;
      bool append = true;
      for (uint r = 0; r < reduced_clusters.size(); r++) {
        append = false;
        for (int i = 0; i < order - 1; i++) {
          if (atoms[i] != reduced_clusters[r][i]) {
            append = true;
            i = order;
          }
        }
        if (!append) {  // If append stays false, the reduced cluster is not new
          r = reduced_clusters.size();
        }
      }
      if (append) {
        for (int i = 0; i < order - 1; i++) cluster[i] = atoms[i];
        reduced_clusters.push_back(cluster);
      }
    }
    return reduced_clusters;
  }

  //getTensorTransformations//////////////////////////////////////////////////
  // This algorithm does two things: it generates the tensor transformations
  // for each cluster to transform the IFCs of the inequivalent clusters; and
  // it generates a list of equivalent IFCs for each inequivalent cluster to
  // calculate the corrections.
  //
  // Check the typedefs at the beginning of the file for tform and v4int
  void AnharmonicIFCs::getTensorTransformations(v4int& eq_ifcs,
      vector<vector<tform> >& transformations) {
    for (uint ineq = 0; ineq < clst.ineq_clusters.size(); ineq++) {
      vector<tform> transform(clst.ineq_clusters[ineq].size() - 1);
      vector<vector<vector<int> > > eq(clst.nifcs, vector<vector<int> >(clst.ineq_clusters[ineq].size()));
      int ind = 0;
      for (int crt = 0; crt < clst.nifcs; crt++) {
        eq[ind][0].push_back(crt);
        ind++;
      }
      for (uint c = 1; c < clst.ineq_clusters[ineq].size(); c++) {
        tform trf;
        _cluster cluster_trans = clst.getCluster(clst.ineq_clusters[ineq][c]);
        int fg = 0, perm = 0, rw = 0, cl = 0, p = 0;
        vector<int> atoms_trans = cluster_trans.atoms;
        atoms_trans[0] = clst.sc2pcMap[atoms_trans[0]];  // transfer to pcell
        fg = cluster_trans.fgroup;
        perm = cluster_trans.permutation;
        for (int itrans = 0; itrans < clst.nifcs; itrans++) {
          std::pair<vector<int>, vector<double> > t;
          int ind_orig = 0;
          for (int iorig = 0; iorig < clst.nifcs; iorig++) {
            double coeff = 1.0;
            for (int o = 0; o < order; o++) {
              rw = cart_indices[itrans][o] + 1;
              p = clst.permutations[perm][o];
              cl = cart_indices[iorig][p] + 1;
              coeff *= clst.pcell.fgroup[fg].Uc[rw][cl];
              if (abs(coeff) < _ZERO_TOL_) {
                coeff = 0.0;
                o = order;
              }
            }
            if (abs(coeff) > _ZERO_TOL_) {
              t.first.push_back(iorig);
              t.second.push_back(coeff);
              eq[ind_orig][c].push_back(itrans);
            }
            ind_orig++;
          }
          trf.push_back(t);
        }
        transform[c-1] = trf;
      }
      eq_ifcs[ineq] = eq;
      transformations[ineq] =  transform;
    }
  }

  // END Initializers

  // BEGIN Iterations
  //applyLinCombs/////////////////////////////////////////////////////////////
  // Sets the lineraly dependent IFCs according to the obtained linear
  // combinations.
  void AnharmonicIFCs::applyLinCombs(vector<vector<double> >& ifcs) {
    int c = 0, cart_ind = 0, cart_ind_indep = 0;
    for (uint ineq = 0; ineq < clst.ineq_clusters.size(); ineq++) {
      c = clst.ineq_clusters[ineq][0];
      _linearCombinations lcomb = clst.linear_combinations[ineq];
      for (uint d = 0; d < lcomb.dependent.size(); d++) {
        cart_ind = lcomb.dependent[d];
        ifcs[c][cart_ind] = 0.0;  // reset
        for (uint lc = 0; lc < lcomb.indices[d].size(); lc++) {
          cart_ind_indep = lcomb.indices[d][lc];
          ifcs[c][cart_ind] += lcomb.coefficients[d][lc] * ifcs[c][cart_ind_indep];
        }
      }
    }
  }

  //transformIFCs/////////////////////////////////////////////////////////////
  // Transforms all IFCs using symmetry. The first two loops go over all
  // equivalent clusters and transform the inequivalent clusters using tensor
  // transformations.
  //
  // See the top of this file for the typedef of tform.
  void AnharmonicIFCs::transformIFCs(const vector<vector<tform> >& transformations,
      vector<vector<double> >& ifcs) {
    int clst_orig = 0, clst_trans = 0, cart_indices_orig = 0;
    double coeff = 0.0;
    for (uint ineq = 0; ineq < clst.ineq_clusters.size(); ineq++) {
      clst_orig = clst.ineq_clusters[ineq][0];
      for (uint c = 1; c < clst.ineq_clusters[ineq].size(); c++) {
        clst_trans = clst.ineq_clusters[ineq][c];
        for (int itrans = 0; itrans < clst.nifcs; itrans++) {
          ifcs[clst_trans][itrans] = 0.0;  // reset
          const std::pair<vector<int>, vector<double> >& transf = transformations[ineq][c-1][itrans];
          for (uint t = 0; t < transf.first.size(); t++) {
            cart_indices_orig = transf.first[t];
            coeff = transf.second[t];
            ifcs[clst_trans][itrans] += coeff * ifcs[clst_orig][cart_indices_orig];
          }
        }
      }
    }
  }

  //calcSums//////////////////////////////////////////////////////////////////
  // Determines the deviations from zero and the sum of the absolute values
  // for each reduced cluster. Both are used for the correction of the IFCs
  // while the deviations from zero are also used as errors for the sum rules.
  void AnharmonicIFCs::calcSums(const vector<vector<int> >& reduced_clusters,
      const vector<vector<double> >& ifcs,
      vector<vector<double> >& dev_from_zero,
      vector<vector<double> >& abssum) {
    dev_from_zero.assign(reduced_clusters.size(), vector<double>(clst.nifcs, 0));
    abssum.assign(reduced_clusters.size(), vector<double>(clst.nifcs, 0));
    for (uint i = 0; i < clst.ineq_clusters.size(); i++) {
      for (uint j = 0; j < clst.ineq_clusters[i].size(); j++) {
        uint r = findReducedCluster(reduced_clusters, clst.getCluster(clst.ineq_clusters[i][j]).atoms);
        for (int c = 0; c < clst.nifcs; c++) {
          dev_from_zero[r][c] += ifcs[clst.ineq_clusters[i][j]][c];
          abssum[r][c] += abs(ifcs[clst.ineq_clusters[i][j]][c]);
        }
      }
    }
  }

  //correctIFCs///////////////////////////////////////////////////////////////
  // Corrects the IFCs the comply with the sum rule using weighted averages.
  //
  // Check the typedef at the beginning of the file for v5int.
  void AnharmonicIFCs::correctIFCs(vector<vector<double> >& ifcs,
      const vector<vector<double> >& dev_from_zero, 
      const vector<vector<double> >& abssum,
      const vector<vector<int> >& reduced_clusters,
      const v4int& eq_ifcs) {
    vector<int> eq;
    for (uint ineq = 0; ineq < clst.ineq_clusters.size(); ineq++) {
      uint nclusters = clst.ineq_clusters[ineq].size();
      int ic = clst.ineq_clusters[ineq][0];
      // Calculate correction terms
      vector<vector<double> > correction_terms(nclusters);
      for (uint c = 0; c < nclusters; c++) {
        correction_terms[c] = getCorrectionTerms(clst.ineq_clusters[ineq][c],
            reduced_clusters, ifcs,
            dev_from_zero, abssum);
      }

      // Correct the linearly independent IFCs
      const _linearCombinations& lcomb = clst.linear_combinations[ineq];
      const vector<int>& indep = lcomb.independent;
      for (uint i = 0; i < indep.size(); i++) {
        int neq = 0, cl = 0;
        double corrected_ifc = 0.0;
        for (uint c = 0; c < nclusters; c++) {
          cl = clst.ineq_clusters[ineq][c];
          eq = eq_ifcs[ineq][indep[i]][c];
          uint eqsize = eq.size();
          for (uint e = 0; e < eqsize; e++) {
            if (ifcs[cl][eq[e]] != 0.0) {
              neq++;
              corrected_ifc += correction_terms[c][eq[e]] * ifcs[ic][indep[i]] / ifcs[cl][eq[e]];
            }
          }
          for (uint dep = 0; dep < lcomb.indep2depMap[i].size(); dep++) {
            eq = eq_ifcs[ineq][lcomb.indep2depMap[i][dep]][c];
            for (uint e = 0; e < eqsize; e++) {
              if (ifcs[cl][eq[e]] != 0.0) {
                neq++;
                corrected_ifc += correction_terms[c][eq[e]] * ifcs[ic][indep[i]] / ifcs[cl][eq[e]];
              }
            }
          }
        }
        if (neq > 0) {
          ifcs[ic][indep[i]] = mixing_coefficient * ifcs[ic][indep[i]] + (1 - mixing_coefficient) * corrected_ifc/((double) neq);
        }
      }
    }
  }

  //getCorrectionTerms////////////////////////////////////////////////////////
  // Calculates the correction term for each IFC.
  vector<double>
    AnharmonicIFCs::getCorrectionTerms(const int& c,
        const vector<vector<int> >& reduced_clusters,
        const vector<vector<double> >& ifcs,
        const vector<vector<double> >& dev_from_zero,
        const vector<vector<double> >& abssum) {
      vector<double> correction_terms = ifcs[c];
      vector<double> correction(clst.nifcs);
      for (int i = 0; i < clst.nifcs; i++) correction[i] = std::abs(correction_terms[i]);
      uint r = findReducedCluster(reduced_clusters, clst.getCluster(c).atoms);
      for (int crt = 0; crt < clst.nifcs; crt++) {
        if (abssum[r][crt] != 0.0) correction[crt] *= dev_from_zero[r][crt]/abssum[r][crt];
        else correction[crt] = 0.0;
        correction_terms[crt] -= correction[crt];
      }
      return correction_terms;
    }

  uint AnharmonicIFCs::findReducedCluster(const vector<vector<int> >& reduced_clusters,
      const vector<int>& atoms) {
    for (uint r = 0; r < reduced_clusters.size(); r++) {
      int at = 0;
      for (at = 0; at < order - 1; at++) {
        if (reduced_clusters[r][at] != atoms[at]) break;
      }
      if (at == order - 1) return r;
    }
    return -1;
  }

  // END Iterations

}  // namespace apl

/********************************* FILE I/O *********************************/

namespace apl {

  // BEGIN Write files
  //writeIFCsToFile///////////////////////////////////////////////////////////
  // Writes the AnharmonicIFCs object and minimal structure information to an
  // XML file.
  void AnharmonicIFCs::writeIFCsToFile(const string& filename) {
    stringstream output;
    // Header
    output << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << std::endl;
    output << "<anharmonicifcs>" << std::endl;

    output << writeParameters();
    output << writeIFCs();
    output << "</anharmonicifcs>" << std::endl;
    aurostd::stringstream2file(output, filename);
    if (!aurostd::FileExist(filename)) {
      string message = "Could not write tensor to file.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
  }

  //writeParameters///////////////////////////////////////////////////////////
  // Writes the calculation parameters and minimal structure information to
  // the XML file.
  string AnharmonicIFCs::writeParameters() {
    stringstream parameters;
    string tab = " ";

    // Info about calculation run
    parameters << tab << "<generator>" << std::endl;
    parameters << tab << tab << "<i name=\"aflow_version\" type=\"string\">" << AFLOW_VERSION << "</i>" << std::endl;
    string time = aflow_get_time_string();
    if (time[time.size() - 1] == '\n') time.erase(time.size() - 1);
    parameters << tab << tab << "<i name=\"date\" type=\"string\">" << time << "</i>" << std::endl;
    //ME20200428 - We do not compare checksums anymore
    //parameters << tab << tab << "<i name=\"checksum\" file=\"" << _AFLOWIN_;
    //parameters << "\" type=\"" << APL_CHECKSUM_ALGO << "\">" << std::hex << aurostd::getFileCheckSum(directory + "/" + _AFLOWIN_ + "", APL_CHECKSUM_ALGO);  //ME20190219
    //parameters.unsetf(std::ios::hex);  //ME20190125 - Remove hexadecimal formatting
    //parameters  << "</i>" << std::endl;
    parameters << tab << "</generator>" << std::endl;

    // Distortion magnitude
    parameters << tab << "<distortion_magnitude units=\"Angstrom\" cs=\"cartesian\">";
    parameters << distortion_magnitude << "</distortion_magnitude>" << std::endl;

    //Order
    parameters << tab << "<order>" << order << "</order>" << std::endl;

    // Iteration parameters
    parameters << tab << "<iteration>" << std::endl;
    // std::dec prevents hexadecimal output
    parameters << tab << tab << "<max_iter>" << std::dec << max_iter << "</max_iter>" << std::endl;
    parameters << tab << tab << "<mixing_coefficient>";
    parameters << std::setprecision(8) << mixing_coefficient;
    parameters << "</mixing_coefficient>" << std::endl;
    parameters << tab << tab << "<sumrule_threshold>";
    parameters << std::setprecision(15) << sumrule_threshold;
    parameters << "</sumrule_threshold>" << std::endl;
    parameters << tab << "</iteration>" << std::endl;

    // Structure
    parameters << tab << "<structure units=\"Angstrom\" cs=\"fractional\">" << std::endl;
    parameters << tab << tab << "<varray name=\"pcell lattice\">" << std::endl;
    for (int i = 1; i < 4; i++) {
      parameters << tab << tab << tab << "<v>";
      for (int j = 1; j < 4; j++) {
        parameters << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
        parameters << std::setprecision(8);
        parameters << std::setw(15) << clst.pcell.lattice[i][j];
      }
      parameters << "</v>" << std::endl;
    }
    parameters << tab << tab << "</varray>" << std::endl;
    parameters << tab << tab << "<varray name=\"positions\">" << std::endl;
    for (uint i = 0; i < clst.pcell.atoms.size(); i++) {
      int t = clst.pcell.atoms[i].type;
      parameters << tab << tab << tab << "<v species=\"" << clst.pcell.species[t] << "\">";
      for (int j = 1; j < 4; j++) {
        parameters << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
        parameters << std::setprecision(8);
        parameters << std::setw(15) << clst.pcell.atoms[i].fpos[j];
      }
      parameters << "</v>" << std::endl;
    }
    parameters << tab << tab << "</varray>" << std::endl;
    parameters << tab << tab << "<varray name=\"supercell\">" << std::endl;
    parameters << tab << tab << tab << "<v>";
    for (int i = 1; i < 4; i++) {
      parameters << tab << clst.sc_dim[i];
    }
    parameters << "</v>" << std::endl;
    parameters << tab << tab << "</varray>" << std::endl;
    parameters << tab << "</structure>" << std::endl;
    return parameters.str();
  }

  //writeIFCs/////////////////////////////////////////////////////////////////
  // Writes the force constants part of the XML file.
  string AnharmonicIFCs::writeIFCs() {
    stringstream ifcs;
    string tab = " ";
    int precision = 15;
    int extra = 5; // first digit, decimal point, minus sign + 2 spaces
    double max_ifc = -1E30;
    for (uint i = 0; i < force_constants.size(); i++) {
      for (uint j = 0; j < force_constants[i].size(); j++) {
        if (std::abs(force_constants[i][j]) > max_ifc) max_ifc = std::abs(force_constants[i][j]);
      }
    }
    int width = (int) log10(max_ifc);
    if (width < 0) {
      width = precision + extra;
    } else {
      width += precision + extra;
    }

    ifcs << tab << "<force_constants>" << std::endl;
    for (uint c = 0; c < clst.clusters.size(); c++) {
      ifcs << tab << tab << "<varray atoms=\""
        << aurostd::joinWDelimiter(clst.clusters[c].atoms, " ") << "\">" << std::endl;
      int crt = 0;
      for (int i = 0; i < clst.nifcs/9; i++) {
        ifcs << tab << tab << tab << "<varray slice=\"" << cart_indices[crt][0];
        for (int j = 1; j < order - 2; j++) {
          ifcs << " " << cart_indices[crt][j];
        }
        ifcs << "\">" << std::endl;
        for (int j = 0; j < 3; j++) {
          ifcs << tab << tab << tab << tab << "<v>";
          for (int k = 0; k < 3; k++) {
            ifcs << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
            ifcs << std::setprecision(precision);
            ifcs << std::setw(width) << force_constants[c][crt];
            crt++;
          }
          ifcs << "</v>" << std::endl;
        }
        ifcs << tab << tab << tab << "</varray>" << std::endl;
      }
      ifcs << tab << tab << "</varray>" << std::endl;
    }
    ifcs << tab << "</force_constants>" << std::endl;
    return ifcs.str();
  }

  // END Write files

} // namespace apl

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow MARCO ESTERS - Duke University 2018-2021           *
// *                                                                         *
// ***************************************************************************
