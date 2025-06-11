// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *            Aflow MARCO ESTERS - Duke University 2018-2021               *
// *                                                                         *
// ***************************************************************************
// Written by Marco Esters, 2018. Based on work by Jose J. Plata (AFLOW AAPL,
// DOI: 10.1038/s41524-017-0046-7) and Jesus Carrete (ShengBTE, 
// DOI: 10.1016/j.cpc.2014.02.015).
//
// This class is a generalized cluster class for N-order force constants.
// It calculates inequivalent atom clusters and distortions, and establishes
// transformation relationships (rotation matrices, linear combinations)
// between the equivalent and inequivalent clusters.
//
// See aflow_apl.h for descriptions of the class members, and for the structs
// structs _cluster and _ineq_distortions.

#include "aflow_apl.h"
using std::vector;
using aurostd::xcombos;
using aurostd::xerror;

#define _DEBUG_AAPL_CLUSTERS_ false

/************************** CONSTRUCTORS/DESTRUCTOR *************************/

namespace apl {

  //Constructors//////////////////////////////////////////////////////////////
  // Default constructor
  ClusterSet::ClusterSet(ostream& oss) : xStream(oss) {
    free();
    directory = "./";
  }

  ClusterSet::ClusterSet(ofstream& mf, ostream& oss) : xStream(mf,oss) {
    free();
    directory = "./";
  }

  ClusterSet::ClusterSet(const Supercell& supercell, int _order, int cut_shell,
      double cut_rad, const string& dir, ofstream& mf, ostream& oss) : xStream(mf,oss) {
    free();
    directory = dir;
    initialize(supercell, _order, cut_shell, cut_rad);
  }

  //From file
  ClusterSet::ClusterSet(const string& filename, const Supercell& supercell, int _order,
      int cut_shell, double cut_rad, const string& dir, ofstream& mf, ostream& oss) : xStream(mf,oss) {
    free();  // Clear old vectors
    directory = dir;
    initialize(supercell, _order, cut_shell, cut_rad);
    readClusterSetFromFile(filename);
  }

  //Copy Constructors/////////////////////////////////////////////////////////
  ClusterSet::ClusterSet(const ClusterSet& that) : xStream(*that.getOFStream(),*that.getOSS()) {
    if (this != &that) free();
    copy(that);
  }

  const ClusterSet& ClusterSet::operator=(const ClusterSet& that) {
    if (this != &that) free();
    copy(that);
    return *this;
  }

  void ClusterSet::copy(const ClusterSet& that) {
    if (this == &that) return;
    xStream::copy(that);
    clusters = that.clusters;
    coordination_shells = that.coordination_shells;
    cutoff = that.cutoff;
    directory = that.directory;
    distortion_vectors = that.distortion_vectors;
    higher_order_ineq_distortions = that.higher_order_ineq_distortions;
    ineq_clusters = that.ineq_clusters;
    ineq_distortions = that.ineq_distortions;
    linear_combinations = that.linear_combinations;
    nifcs = that.nifcs;
    order = that.order;
    pcell = that.pcell;
    pc2scMap = that.pc2scMap;
    permutations = that.permutations;
    scell = that.scell;
    sc2pcMap = that.sc2pcMap;
    sc_dim = that.sc_dim;
    symmetry_map = that.symmetry_map;
  }

  //Destructor////////////////////////////////////////////////////////////////
  ClusterSet::~ClusterSet() {
    xStream::free();
    free();
  }

  //free//////////////////////////////////////////////////////////////////////
  // Clears all vectors and resets all values.
  void ClusterSet::free() {
    clusters.clear();
    coordination_shells.clear();
    cutoff = 0.0;
    directory = "";
    distortion_vectors.clear();
    higher_order_ineq_distortions.clear();
    ineq_clusters.clear();
    ineq_distortions.clear();
    linear_combinations.clear();
    nifcs = 0;
    order = 0;
    pc2scMap.clear();
    pcell.clear();
    permutations.clear();
    scell.clear();
    sc_dim.clear();
    sc2pcMap.clear();
    symmetry_map.clear();
  }


  //clear/////////////////////////////////////////////////////////////////////
  // Creates an empty ClusterSet object.
  void ClusterSet::clear() {
    free();
  }

  //initialize////////////////////////////////////////////////////////////////
  // Initialize basic parameters.
  void ClusterSet::initialize(const Supercell& supercell, int _order, int cut_shell, double cut_rad, ofstream& mf, ostream& oss) {
    xStream::initialize(mf, oss);
    initialize(supercell, _order, cut_shell, cut_rad);
  }

  void ClusterSet::initialize(const Supercell& supercell, int _order, int cut_shell, double cut_rad) {
    stringstream message;
    if (_order > 1) {
      order = _order;
    } else {
      message << "Cluster order must be larger than 1 (is " << _order << ").";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_RANGE_);
    }
    nifcs = aurostd::powint(3, order);
    scell = supercell.getSupercellStructure();
    pcell = supercell.getInputStructure();
    sc_dim = supercell.scell_dim;
    if (sc_dim.rows != 3) {
      message << "Cannot use non-diagonal supercells for AAPL.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_ILLEGAL_);
    }
    pc2scMap = supercell._pc2scMap;
    sc2pcMap = supercell._sc2pcMap;
    symmetry_map = getSymmetryMap();
    if (cut_shell > 0) {
      double max_rad = getMaxRad(scell, cut_shell);
      if (max_rad > cut_rad) {
        message << "Cutoff has been modified to " << max_rad << " Angstroms.";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, directory, *p_FileMESSAGE, *p_oss);
        cut_rad = max_rad;
      }
    }
    cutoff = cut_rad;
    permutations = getPermutations(order);

    // Clear cluster data
    clusters.clear();
    coordination_shells.clear();
    higher_order_ineq_distortions.clear();
    ineq_clusters.clear();
    ineq_distortions.clear();
    linear_combinations.clear();
  }

  //xStream initializers
  void ClusterSet::initialize(ostream& oss) {
    xStream::initialize(oss);
  }

  void ClusterSet::initialize(ofstream& mf, ostream& oss) {
    xStream::initialize(mf, oss);
  }

}  // namespace apl

/******************************* INTERFACE **********************************/

namespace apl {

  //getCluster////////////////////////////////////////////////////////////////
  // Returns a cluster from the list of clusters
  const _cluster& ClusterSet::getCluster(const int& i) const {
    return clusters[i];
  }

}

/************************** INITIAL CALCULATIONS ****************************/

namespace apl {

  //getSymmetryMap////////////////////////////////////////////////////////////
  // Creates a symmetry map for a list of atoms in the supercell using scaled
  // symmetry operations. Note that the supercell and the primtive cell must
  // have the same symmetry.
  vector<vector<int> > ClusterSet::getSymmetryMap() {
    bool LDEBUG = (FALSE || XHOST.DEBUG || _DEBUG_AAPL_CLUSTERS_);
    string message = "";
    // Check if the symmetry of the supercell and the primitive cell are the
    // same by comparing the crystal point groups.
    if (!scell.pgroup_xtal_calculated) {
      scell.CalculateSymmetryPointGroupCrystal(false);
      if (!scell.pgroup_xtal_calculated) {
        message = "Could not calculate the point group of the supercell.";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
    }
    if (!pcell.pgroup_xtal_calculated) {
      pcell.CalculateSymmetryPointGroupCrystal(false);
      if (!pcell.pgroup_xtal_calculated) {
        message = "Could not calculate the point group of the primitive cell.";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
    }

    if (SYM::PointGroupsIdentical(scell.pgroup_xtal, pcell.pgroup_xtal, scell.sym_eps, false)) {
      // Make sure that the factor group is calculated
      if (!pcell.fgroup_calculated) {
        pcell.CalculateSymmetryFactorGroup(false);
        if (!pcell.fgroup_calculated) {
          message = "Could not calculate the factor group of the primitive cell.";
          throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
        }
      }
      bool mapped = false;
      // Minimum distance for pcell is the same for scell and cheaper to compute
      pcell.dist_nn_min = SYM::minimumDistance(pcell.atoms, pcell.lattice);
      double tol = _ZERO_TOL_;
      bool skewed = SYM::isLatticeSkewed(scell.lattice, pcell.dist_nn_min, tol);
      uint fgsize = pcell.fgroup.size();
      uint natoms = scell.atoms.size();
      vector<vector<int> > sym_map(fgsize, vector<int>(natoms,  -1));

      for (uint fg = 0; fg < fgsize; fg++) {
        // Scale the translation vector to supercell dimensions
        xvector<double> ftau_scaled(3), fpos_scaled(3);
        for (int i = 1; i < 4; i++) {
          ftau_scaled[i] = pcell.fgroup[fg].ftau[i]/sc_dim[i];
        }
        // Map atoms
        for (uint atsc = 0; atsc < natoms; atsc++) {
          mapped = false;
          fpos_scaled = (pcell.fgroup[fg].Uf * scell.atoms[atsc].fpos) + ftau_scaled;
          for (uint at_map = 0; at_map < natoms && !mapped; at_map++) {
            if (SYM::FPOSMatch(fpos_scaled, scell.atoms[at_map].fpos,
                  scell.lattice, scell.f2c, skewed, scell.sym_eps)) { //DX20190619 - lattice and f2c as input
              sym_map[fg][atsc] = at_map;
              mapped = true;
            }
          }
          if (!mapped) {
            message = "At least one atom of the supercell could not be mapped.";
            if (LDEBUG) {
            }
            throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
          }
        }
      }
      return sym_map;
    } else {
      message = "The point group of supercell is different than the point group of the input structure.";
      message += " This feature is not implemented yet.";
      if (LDEBUG) {
        std::cerr << "ClusterSet::getSymmetryMap: Point group mismatch. Primitive cell: ";
        std::cerr << pcell.point_group_Hermann_Mauguin << "; supercell: " << scell.point_group_Hermann_Mauguin << std::endl;
      }
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
  }

  ///getMaxRad////////////////////////////////////////////////////////////////
  // Returns the largest coordination shell for all inequivalent atoms using
  // cut_shell coordination shells.
  double ClusterSet::getMaxRad(const xstructure& cell, int cut_shell){
    bool LDEBUG = (FALSE || XHOST.DEBUG || _DEBUG_AAPL_CLUSTERS_) || _DEBUG_AAPL_CLUSTERS_;
    if (cut_shell > 0) {
      int countshell = 0;
      double shell_rad = 0.0;
      double max_rad = 0.0;
      uint natoms = cell.atoms.size();
      xvector<double> distances(natoms - 1, 0);
      for (uint i = 0; i < cell.iatoms.size(); i++) {
        countshell = 0;
        int at1 = cell.iatoms[i][0];
        for (uint at2 = 0; at2 < natoms; at2++) {
          //DX20190613 [OBSOLETE - changed function name] distances[at2] = SYM::minimumCartesianDistance(cell.atoms[at1].cpos,
          //DX20190613 [OBSOLETE - changed function name]                                                cell.atoms[at2].cpos,
          //DX20190613 [OBSOLETE - changed function name]                                                cell.lattice);
          distances[at2] = aurostd::modulus(SYM::minimizeDistanceCartesianMethod(cell.atoms[at1].cpos,
                cell.atoms[at2].cpos,
                cell.lattice)); //DX20190613
        }
        aurostd::sort(distances);
        for (uint j = 1; j < natoms; j++) {
          if (distances[j] > distances[j - 1] + _APL_SHELL_TOL_ ){
            countshell++;
          }
          if (countshell == cut_shell) {
            shell_rad = distances[j] + _APL_SHELL_TOL_;
            j = natoms;
            if (shell_rad > max_rad){
              max_rad = shell_rad;
            }
          }
        }
      }
      if (LDEBUG) {
        std::cerr << "ClusterSet::getMaxRad: Modified cutoff: " << max_rad << " Angstrom." << std::endl;
      }
      return max_rad;
    } else {
      return 0.0;
    }
  }

  //buildShells///////////////////////////////////////////////////////////////
  // Creates coordination shells around the atoms of the primitive cell.
  void ClusterSet::buildShells() {
    bool LDEBUG = (FALSE || XHOST.DEBUG || _DEBUG_AAPL_CLUSTERS_);
    string message = "Building coordination shells.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, directory, *p_FileMESSAGE, *p_oss);
    int at1sc = 0;
    vector<int> shell;
    coordination_shells.clear();
    uint natoms_pcell = pcell.atoms.size();
    uint natoms_scell = scell.atoms.size();
    for (uint at1 = 0; at1 < natoms_pcell; at1++) {
      shell.clear();
      at1sc = pc2scMap[at1];
      shell.push_back(at1sc);
      for (uint at2 = 0; at2 < natoms_scell; at2++) {
        //DX20190613 [OBSOLETE - changed function name] xvector<double> vec; xvector<int> ijk;
        //DX20190613 [OBSOLETE - changed function name] double cart_dist = SYM::minimumCartesianDistance(scell.atoms[at1sc].cpos,
        //DX20190613 [OBSOLETE - changed function name]                                                  scell.atoms[at2].cpos,
        //DX20190613 [OBSOLETE - changed function name]                                                  scell.lattice, vec, ijk);
        double cart_dist = aurostd::modulus(SYM::minimizeDistanceCartesianMethod(scell.atoms[at1sc].cpos,
              scell.atoms[at2].cpos,
              scell.lattice)); //DX20190613
        if (cart_dist < cutoff) {
          shell.push_back(at2);
        }
      }
      coordination_shells.push_back(shell);
      if (LDEBUG) {
        uint shellsize = shell.size();
        std::cerr << "ClusterSet::buildShells: Found " << shellsize - 1;
        std::cerr << " atoms in shell around atom " << shell[0] << "." << std::endl;
        std::cerr << "ClusterSet::buildShells: Atoms in shell:";
        for (uint i = 1; i < shellsize; i++) {
          std::cerr << " " << shell[i];
        }
        std::cerr << std::endl;
      }
    }
  }

}  // namespace apl

/********************************** BUILD ***********************************/

namespace apl {

  //build/////////////////////////////////////////////////////////////////////
  // Builds the inequivalent clusters.
  void ClusterSet::build() {
    bool LDEBUG = (FALSE || XHOST.DEBUG || _DEBUG_AAPL_CLUSTERS_);
    stringstream message;
    message << "Building clusters of order " << order << ".";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, directory, *p_FileMESSAGE, *p_oss);
    buildShells();
    clusters = buildClusters();
    getInequivalentClusters(clusters, ineq_clusters);

    if (LDEBUG) {
      std::cerr << "ClusterSet::build: Inequivalent clusters built:" << std::endl;
      for (uint i = 0; i < ineq_clusters.size(); i++) {
        std::cerr << "Inequivalent cluster " << i << std::endl;
        for (uint j = 0; j < ineq_clusters[i].size(); j++) {
          for (uint k = 0; k < getCluster(ineq_clusters[i][j]).atoms.size(); k++) {
            std::cerr << getCluster(ineq_clusters[i][j]).atoms[k] << " ";
          }
          std::cerr << std::endl;
        }
        std::cerr << std::endl;
      }
    }
  }

  //getPermutations///////////////////////////////////////////////////////////
  // Initialize permutations.
  vector<vector<int> > ClusterSet::getPermutations(int order) {
    vector<int> permut_base(order);
    vector<vector<int> > permuts;
    for (int i = 0; i < order; i++) {
      permut_base[i] = i;
    }
    xcombos combos(permut_base);
    while (combos.increment()) {
      permuts.push_back(combos.getCombo());
    }
    return permuts;
  }

  //buildClusters/////////////////////////////////////////////////////////////
  // Builds a list of clusters around the central atoms. No symmetry reduction
  // is performed here.
  vector<_cluster> ClusterSet::buildClusters() {
    bool LDEBUG = (FALSE || XHOST.DEBUG || _DEBUG_AAPL_CLUSTERS_);
    vector<_cluster> cluster_list;
    vector<int> atoms_in_cluster(order), clst;
    bool repeat = true;
    for (uint s = 0; s < coordination_shells.size(); s++) {
      atoms_in_cluster[0] = coordination_shells[s][0];
      xcombos clst_combos(coordination_shells[s].size() - 1, order - 1, 'E', repeat);
      while (clst_combos.increment()) {
        bool append = true;
        clst = clst_combos.getCombo();
        for (int o = 1; o < order; o++) {
          int at = coordination_shells[s][clst[o-1] + 1];
          for (int d = 1; d < o; d++) {
            int atc = atoms_in_cluster[d];
            //DX20190613 [OBSOLETE - changed function name] double mdist = SYM::minimumCartesianDistance(scell.atoms[at].cpos,
            //DX20190613 [OBSOLETE - changed function name]                                              scell.atoms[atc].cpos,
            //DX20190613 [OBSOLETE - changed function name]                                              scell.lattice);
            double mdist = aurostd::modulus(SYM::minimizeDistanceCartesianMethod(scell.atoms[at].cpos,
                  scell.atoms[atc].cpos,
                  scell.lattice)); //DX20190613
            if (mdist > cutoff) {
              append = false;
              d = o;
              o = order;
            }
          }
          atoms_in_cluster[o] = at;
        }
        if (append) {
          _cluster cluster;
          cluster.atoms = atoms_in_cluster;
          cluster_list.push_back(cluster);
          if (LDEBUG) {
            std::cerr << "ClusterSet::buildClusters: Built cluster";
            for (uint c = 0; c < cluster.atoms.size(); c++) {
              std::cerr << " " << cluster.atoms[c];
            }
            std::cerr << std::endl;
          }
        }
      }
    }
    return cluster_list;
  }

  // BEGIN getInequivalentClusters functions
  //getInequivalentClusters///////////////////////////////////////////////////
  // Builds a list of inequivalent clusters and their transformations. This is
  // done in two steps: First, equivalenceCluster determines if the cluster in
  // question can be transformed into an inequivalent cluster. However, for
  // the force constants, it is important to know how the inequivalent cluster
  // transforms into the equivalent cluster, which is determined by getSymOp.
  // While it is possible to determine this in one step, the two-step
  // procedure is significantly faster, especially with many coordination
  // shells.
  void ClusterSet::getInequivalentClusters(vector<_cluster>& clusters,
      vector<vector<int> > & ineq_clst) {
    string message = "Determining inequivalent clusters.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, directory, *p_FileMESSAGE, *p_oss);
    int equivalent_clst = -1;
    vector<int> unique_atoms;
    vector<vector<int> > compositions;
    uint nspecies = scell.species.size();
    ineq_clst.clear();
    for (uint c = 0; c < clusters.size(); c++) {
      // Determine if cluster is equivalent to an inequivalent cluster
      equivalent_clst = equivalenceCluster(clusters[c].atoms, unique_atoms,
          compositions, ineq_clst);
      if (equivalent_clst > -1 ) {
        // Equivalent - determine transformation
        getSymOp(clusters[c], getCluster(ineq_clusters[equivalent_clst][0]).atoms);
        ineq_clst[equivalent_clst].push_back(c);
      } else {
        // Not equivalent - create new inequivalent cluster
        clusters[c].fgroup = -1;
        clusters[c].permutation = -1;
        vector<int> ineq_init;
        ineq_init.push_back(c);
        ineq_clst.push_back(ineq_init);
        int unique = getNumUniqueAtoms(clusters[c].atoms);
        unique_atoms.push_back(unique);
        if (nspecies > 1) {
          vector<int> comp = getComposition(clusters[c].atoms);
          compositions.push_back(comp);
        }
      }
    }
  }

  //getNumUniqueAtoms/////////////////////////////////////////////////////////
  // Gets the number of unique atoms in a cluster. This greatly speeds up the
  // determination of inequivalent clusters because clusters with different
  // numbers of unique atoms cannot be equivalent and thus do not need to be
  // compared with each other.
  int ClusterSet::getNumUniqueAtoms(const vector<int>& atoms) {
    int unique = 1;
    for (int at = 1; at < order; at++) {
      bool found = false;
      for (int i = 0; i < at; i++) {
        if (atoms[at] == atoms[at - i]) {
          found = true;
          i = at;
        }
      }
      if (!found) {
        unique++;
      }
    }
    return unique;
  }

  //getComposition////////////////////////////////////////////////////////////
  // Determines the composition of a cluster. This greatly speeds up the
  // determination of inequivalent clusters because clusters with different
  // compoistion cannot be equivalent and thus do not need to be compared
  // with each other.
  vector<int> ClusterSet::getComposition(const vector<int>& atoms) {
    vector<int> composition(scell.species.size());
    int type = -1;
    for (int at = 0; at < order; at++) {
      type = scell.atoms[atoms[at]].type;
      composition[type]++;
    }
    return composition;
  }

  //sameCompositions//////////////////////////////////////////////////////////
  // Compares the compositions of two clusters. If they are not the same, the
  // symmetry check can be skipped.
  bool ClusterSet::sameComposition(const vector<int>& comp1,
      const vector<int>& comp2) {
    for (uint s = 0; s < scell.species.size(); s++) {
      if (comp1[s] != comp2[s]) {
        return false;
      }
    }
    return true;
  }

  //equivalenceCluster////////////////////////////////////////////////////////
  // Determines whether a cluster is equivalent to another inequivalent
  // cluster. The output is the index of the inequivalent cluster the cluster
  // it is equivalent to (-1 if none).
  int ClusterSet::equivalenceCluster(const vector<int>& atoms,
      const vector<int>& unique_atoms,
      const vector<vector<int> >& compositions,
      const vector<vector<int> >& ineq_clst) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || _DEBUG_AAPL_CLUSTERS_);
    if (LDEBUG) {
      std::cerr << "ClusterSet::equivalenceCluster: Determining equivalent cluster for:";
      for (int a = 0; a < order; a++) {
        std::cerr << " " << atoms[a];
      }
      std::cerr << std::endl;
    }
    uint nspecies = scell.species.size();
    int unique = getNumUniqueAtoms(atoms);
    vector<int> composition;
    if (nspecies > 1) {
      composition = getComposition(atoms);
    }
    vector<int> cluster_transl(order), cluster_transformed(order);
    if (ineq_clst.size() > 0) {
      int perm_index = -1;
      for (uint fg = 0; fg < pcell.fgroup.size(); fg++) {
        for (int at = 0; at < order; at++) {
          cluster_transformed[at] = symmetry_map[fg][atoms[at]];
        }
        for (int at = 0; at < order; at++) {
          cluster_transl = translateToPcell(cluster_transformed, at);
          for (uint ineq = 0; ineq < ineq_clst.size(); ineq++) {
            bool same_composition = true;
            if (nspecies > 1) {
              same_composition = sameComposition(composition, compositions[ineq]);
            } else {
              same_composition = true;
            }
            if ((unique == unique_atoms[ineq]) && (same_composition)) {  
              perm_index = comparePermutations(getCluster(ineq_clst[ineq][0]).atoms, cluster_transl);
              if (perm_index != -1) {
                return ineq;
              }
            }
          }
        }
      }
    }
    return -1;
  }

  //translateToPcell//////////////////////////////////////////////////////////
  // Translates the cluster so that the atom with index "at" is inside the
  // primitive cell.
  vector<int> ClusterSet::translateToPcell(const vector<int>& atoms, int at) {
    xvector<double> translation_vector, fpos_transl;
    vector<int> atoms_translated(order);
    bool update = false;
    int atompc = sc2pcMap[atoms[at]];
    int atomsc = pc2scMap[atompc];
    translation_vector = scell.atoms[atomsc].fpos - scell.atoms[atoms[at]].fpos;
    for (uint i = 0; i < atoms.size(); i++) {
      if ((int)i == at) {
        atoms_translated[i] = atomsc;
      } else {
        fpos_transl = translation_vector + scell.atoms[atoms[i]].fpos;
        for (uint j = 0; j < scell.atoms.size(); j++) {
          // No need to compare when atoms are of different type
          if (scell.atoms[j].type == scell.atoms[atoms[i]].type) {
            if (SYM::FPOSMatch(fpos_transl, scell.atoms[j].fpos, scell.lattice,
                  scell.f2c, update, scell.sym_eps)) { //DX20190619 - lattice and f2c as input
              atoms_translated[i] = j;
              j = scell.atoms.size();
            }
          }
        }
      }
    }
    return atoms_translated;
  }

  //comparePermutations///////////////////////////////////////////////////////
  // Compares if clusters are equivalent by permutation and returns the index
  // of the permutation (-1 if none). The algorithm requires the permutations
  // to be in lexicographical order, which is done by xcombos.
  int ClusterSet::comparePermutations(const vector<int>& cluster_ineq,
      const vector<int>& cluster_compare) {
    int index = 0, iper = 0;
    uint depth = 0, jump = 0;
    uint n = cluster_ineq.size();
    uint jump_size = permutations.size()/n;
    while ((jump < n - depth) && (depth < n)) {
      if (atomsMatch(permutations[iper], cluster_compare, cluster_ineq, depth)) {
        index += jump_size * jump;
        depth++;
        if (depth != n) {jump_size /= (n - depth);}
        jump = 0;
      } else {
        jump++;
        iper += jump_size;
      }
    }
    if (depth == n) {
      return index;
    } else {
      return -1;
    }
  }

  //atomsMatch////////////////////////////////////////////////////////////////
  // Auxilliary function for comparePermutations to make the code more
  // legible. Returns true if the indices after permutation are the same.
  bool ClusterSet::atomsMatch(const vector<int>& permutation, const vector<int>& cluster_perm,
      const vector<int>& cluster_orig, const int& depth) {
    return (cluster_perm[permutation[depth]] == cluster_orig[depth]);
  }

  //getSymOp//////////////////////////////////////////////////////////////////
  // Determines the factor group and the permutation index that transforms
  // the inequivalent cluster into the equivalent cluster.
  void ClusterSet::getSymOp(_cluster& cluster, const vector<int>& atoms_orig) {
    vector<int> cluster_transformed(order), cluster_transl(order);
    int perm_index = -1;
    for (uint fg = 0; fg < pcell.fgroup.size(); fg++) {
      for (int at = 0; at < order; at++) {
        cluster_transformed[at] = symmetry_map[fg][atoms_orig[at]];
      }
      for (int at = 0; at < order; at++) {
        cluster_transl = translateToPcell(cluster_transformed, at);
        perm_index = comparePermutations(cluster.atoms, cluster_transl);
        if (perm_index != -1) {
          cluster.fgroup = fg;
          cluster.permutation = perm_index;
          at = order;
          fg = pcell.fgroup.size();
        }
      }
    }
  }
  // END getInequivalentClusters functions

}  // namespace apl

/******************************* DISTORTIONS ********************************/

namespace apl {

  //buildDistortions//////////////////////////////////////////////////////////
  // Builds the inequivalent distortions for the AAPL calculations.
  void ClusterSet::buildDistortions() {
    string message = "Getting inequivalent distortions.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, directory, *p_FileMESSAGE, *p_oss);
    distortion_vectors = getCartesianDistortionVectors();
    linear_combinations = getLinearCombinations();
    ineq_distortions = initializeIneqDists();
    vector<vector<int> > test_distortions;
    for (uint i = 0; i < ineq_distortions.size(); i++) {
      test_distortions = getTestDistortions(ineq_distortions[i].clusters);
      getInequivalentDistortions(test_distortions, ineq_distortions[i]);
    }
    if (order > 3) {
      higher_order_ineq_distortions = getHigherOrderDistortions();
    }
  }

  //getCartesianDistortionVectors/////////////////////////////////////////////
  // Builds vectors along the plus and minus direction of the Cartesian axes.
  vector<xvector<double> > ClusterSet::getCartesianDistortionVectors() {
    bool LDEBUG = (FALSE || XHOST.DEBUG || _DEBUG_AAPL_CLUSTERS_);
    vector<xvector<double> > dist_vecs;
    xvector<double> vec(3);
    for (int i = 1; i < 4; i++) {
      for (int j = 1; j < 4; j++) {
        if (i == j) {
          vec(j) = 1.0;
        } else {
          vec(j) = 0.0;
        }
      }
      dist_vecs.push_back(vec);
    }
    for (int i = 0; i < 3; i++) {
      vec = -dist_vecs[i];
      dist_vecs.push_back(vec);
    }
    if (LDEBUG) {
      std::cerr << "ClusterSet::getCartesianDistortionVectors: Built distortion vectors" << std::endl;
      for (uint i = 0; i < dist_vecs.size(); i++) {
        for (int j = 1; j < 4; j++) {
          std::cerr << dist_vecs[i][j] << " ";
        }
        std::cerr << std::endl;
      }
    }
    return dist_vecs;
  }

  //initializeIneqDists///////////////////////////////////////////////////////
  // Initializes the inequivalent distortion objects and groups all
  // inequivalent clusters that require the same distortions.
  vector<_ineq_distortions> ClusterSet::initializeIneqDists() {
    vector<_ineq_distortions> ineq_dists;
    for (uint ineq = 0; ineq < ineq_clusters.size(); ineq++) {
      int same = sameDistortions(getCluster(ineq_clusters[ineq][0]), ineq_dists);
      if (same == -1) {
        _ineq_distortions dists;
        dists.clusters.clear();
        dists.clusters.push_back(ineq);
        dists.atoms.resize(order - 1);
        for (int at = 0; at < order - 1; at++) {
          dists.atoms[at] = getCluster(ineq_clusters[ineq][0]).atoms[at];
        }
        ineq_dists.push_back(dists);
      } else {
        ineq_dists[same].clusters.push_back(ineq);
      }
    }
    return ineq_dists;
  }

  //sameDistortions///////////////////////////////////////////////////////////
  // Checks if the atoms that would be distorted for the inequivalent cluster
  // are already distorted by another cluster. Returns -1 if none are found.
  int ClusterSet::sameDistortions(const _cluster& clst,
      const vector<_ineq_distortions>& ineq_dists) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || _DEBUG_AAPL_CLUSTERS_);
    int same = -1;
    for (uint ineq = 0; ineq < ineq_dists.size(); ineq++) {
      same = (int) ineq;
      for (int at = 0; at < order - 1; at++) {
        if (clst.atoms[at] != ineq_dists[ineq].atoms[at]) {
          same = -1;
          at = order;
        }
      }
      if (same > -1) {
        ineq = ineq_dists.size(); 
      }
    }
    if (LDEBUG) {
      std::cerr << "ClusterSet::sameDistortions: same = " << same << std::endl;
    }
    return same;
  }

  //getTestDistortions////////////////////////////////////////////////////////
  // Gets the test distortions from the linear dependences of the interatomic
  // force constants. Only the distortions that belong to linearly independent
  // force constants are kept.
  vector<vector<int> > ClusterSet::getTestDistortions(const vector<int>& clusters) {
    vector<int> dists;
    uint max_dists = (uint) nifcs/3;

    // Determine for which IFCs distortions need to be calculated
    for (uint clst = 0; clst < clusters.size(); clst++) {
      int c = clusters[clst];
      vector<int> independent = linear_combinations[c].independent;
      for (uint i = 0; i < independent.size(); i++) {
        int dist = independent[i]/3;
        bool append = true;
        for (uint d = 0; d < dists.size(); d++) {
          if (dist == dists[d]) {
            append = false;
            d = max_dists;
          }
        }
        if (append) {
          dists.push_back(dist);
        }
        if (dists.size() == max_dists) {
          i = independent.size();
          clst = clusters.size();
        }
      }
    }

    // Determine the actual distortions that are necessary
    vector<vector<int> > test_distortions;
    vector<int> distortions(order - 1), distortions_base(order -1), signs;
    for (uint d = 0; d < dists.size(); d++) {
      int dist = dists[d];
      for (int i = order - 2; i >= 0; i--) {
        distortions_base[i] = dist % 3;
        dist /= 3;
      }
      // Bit enumerator indicating the signs of the distortions
      aurostd::xcombos bitenum(2, order - 1, 'E', true);
      while (bitenum.increment()) {
        signs = bitenum.getCombo();
        for (int i = 0; i < order - 1; i++) {
          distortions[i] = 3 * signs[i] + distortions_base[i];
        }
        test_distortions.push_back(distortions);
      }
    }
    return test_distortions;
  }

  // BEGIN Inequivalent distortions functions
  //getInequivalentDistortions////////////////////////////////////////////////
  // Determines a set of inequivalent distortions and the distortions that
  // are equivalent by symmetry.
  void ClusterSet::getInequivalentDistortions(const vector<vector<int> >& dists,
      _ineq_distortions& ineq_dists) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || _DEBUG_AAPL_CLUSTERS_);
    for (uint td = 0; td < dists.size(); td++) {
      if (ineq_dists.distortions.size() == 0) {
        appendDistortion(ineq_dists, dists[td]);
        continue;
      }
      if (!allZeroDistortions(ineq_dists.atoms, dists[td])) {
        int eq = -1;
        for (uint fg = 0; fg < pcell.fgroup.size(); fg++) {
          if (allAtomsMatch(fg, ineq_dists.atoms)) {
            eq = equivalenceDistortions(pcell.fgroup[fg].Uc, dists[td],
                ineq_dists.distortions, ineq_dists.atoms);
            if (eq > -1) {
              appendDistortion(ineq_dists, dists[td], eq, fg);
              if (LDEBUG) {
                std::cerr << "ClusterSet::getInequivalentDistortions: Distortion ";
                for (uint d = 0; d < dists[td].size(); d++) {
                  std::cerr << " " << dists[td][d];
                }
                std::cerr << " equivalent to";
                for (uint d = 0; d < dists[td].size(); d++) {
                  std::cerr << " " << ineq_dists.distortions[eq][0][d];
                }
                std::cerr << " using factor group " << fg << "." << std::endl;
              }
              fg = pcell.fgroup.size();
            }
          }
        }
        if (eq == -1) {
          appendDistortion(ineq_dists, dists[td]);
          if (LDEBUG) {
            std::cerr << "ClusterSet::getInequivalentDistortions: No equivalent distortion found for";
            for (uint d = 0; d < dists[td].size(); d++) {
              std::cerr << " " << dists[td][d];
            }
            std::cerr << "." << std::endl;
          }
        }
      }
    }
  }

  //appendDistortion//////////////////////////////////////////////////////////
  // Append the distortion either to its equivalent distortion or creates a
  // new inequivalent distortion.
  void ClusterSet::appendDistortion(_ineq_distortions& ineq_dists, vector<int> dist,
      const int& eq, const int& fg) {
    vector<int> transformation_map;
    if (eq == -1) {
      // No equivalent distortion found - create a new inequivalent distortion
      vector<vector<int> > init_dist, init_map;
      vector<int> init_trans;
      init_dist.push_back(dist);
      init_map.push_back(transformation_map);
      init_trans.push_back(-1);
      ineq_dists.distortions.push_back(init_dist);
      ineq_dists.rotations.push_back(init_trans);
      ineq_dists.transformation_maps.push_back(init_map);
    } else {
      // Equivalent distortion found - append to inequivalent distortion
      ineq_dists.distortions[eq].push_back(dist);
      ineq_dists.rotations[eq].push_back(fg);
      transformation_map = getTransformationMap(fg, ineq_dists.atoms[0]);
      ineq_dists.transformation_maps[eq].push_back(transformation_map);
    }
  }

  //allZeroDistortions////////////////////////////////////////////////////////
  // Returns whether the set of distortions cancel each other out (i.e. the
  // total distortions are all zero vectors) since the force constant is then
  // zero and does not need to be considered. This algorithm assumes that
  // atoms are in lexicographical order, which is true when the clusters were
  // created using xcombos.
  bool ClusterSet::allZeroDistortions(const vector<int>& atoms, const vector<int>& dist) {
    int atom = 0, at_count = 0;
    uint natoms = atoms.size();
    for (uint at = 0; at < natoms; at++) {
      atom = atoms[at];
      at_count = 1;
      while ((at + at_count < natoms) && (atoms[at+at_count] == atom)) {
        at_count++;
      }
      if (at_count % 2 == 1) { // Cannot cancel out with odd number of instances
        return false;
      } else {
        xvector<double> total_distortion(3);
        for (uint i = at; i < at + at_count; i++) {
          int d = dist[i];
          total_distortion += distortion_vectors[d];
        }
        if (!aurostd::iszero(total_distortion, _ZERO_TOL_)) {
          return false;
        }
        at += (at_count - 1);
      }
    }
    return true;
  }

  //allAtomsMatch/////////////////////////////////////////////////////////////
  // Tests whether all atoms transform into themselves under factor group fg.
  bool ClusterSet::allAtomsMatch(const int& fg, const vector<int>& atoms) {
    vector<int> atoms_transformed(order -1);
    for (int at = 0; at < order - 1; at++) {
      atoms_transformed[at] = symmetry_map[fg][atoms[at]];
    }
    atoms_transformed = translateToPcell(atoms_transformed, 0);
    for (int at = 0; at < order - 1; at++) {
      if (atoms_transformed[at] != atoms[at]) {
        return false;
      }
    }
    return true;
  }

  //equivalenceDistortions////////////////////////////////////////////////////
  // Determines whether a distortion is equivalent to another inequivalent
  // distortion and returns the index of that distortion (-1 if not equivalent
  // to any).
  int ClusterSet::equivalenceDistortions(const xmatrix<double>& Uc,
      const vector<int>& dist,
      const vector<vector<vector<int> > >& ineq_dists,
      const vector<int>& atoms) {
    uint distsize = dist.size();
    uint natoms = atoms.size();
    vector<xvector<double> > trans_dist(distsize);
    int eq = -1;
    for (uint d = 0; d < distsize; d++) {
      int v = dist[d];
      trans_dist[d] = distortion_vectors[v];
    }
    for (uint ineq = 0; ineq < ineq_dists.size(); ineq++) {
      for (uint at = 0; at < natoms; at++) {
        int at_count = 1;
        xvector<double> vec_dist, vec_comp;
        vec_dist = trans_dist[at];
        int v = ineq_dists[ineq][0][at];
        vec_comp = distortion_vectors[v];
        while ((at + at_count < natoms) && (atoms[at] == atoms[at+at_count])) {
          vec_dist += trans_dist[at+at_count];
          int v = ineq_dists[ineq][0][at+at_count];
          vec_comp += distortion_vectors[v];
          at_count++;
        }
        vec_comp = Uc * vec_comp;

        // Normalize so that each component is -1, 0, or 1.
        for (int i = 1; i < 4; i++) {
          if (abs(vec_dist(i)) > _ZERO_TOL_) vec_dist(i) /= abs(vec_dist(i));
          if (abs(vec_comp(i)) > _ZERO_TOL_) vec_comp(i) /= abs(vec_comp(i));
        }

        if (aurostd::iszero(vec_dist - vec_comp, _ZERO_TOL_)) {
          eq = ineq;
          at += (at_count - 1);
        } else {
          eq = -1;
          at = natoms;
        }
      }
      if (eq > -1) {
        ineq = ineq_dists.size();
      }
    }
    return eq;
  }

  //getTransformationMap//////////////////////////////////////////////////////
  // Obtains the transformation map for the equivalent distortion.
  vector<int> ClusterSet::getTransformationMap(const int& fg, const int& atom) {
    uint natoms = scell.atoms.size();
    vector<int> transformation_map(natoms, -1);
    int atom_trans = symmetry_map[fg][atom];
    xvector<double> transl = scell.atoms[atom_trans].cpos - scell.atoms[atom].cpos;
    transl -= pcell.fgroup[fg].ctau;
    // If the symmetry operation does not include a translation, the already
    // determined symmetry map can be used. Otherwise, the transformation map
    // needs to be determined using the supercell factor group
    if (aurostd::iszero(transl, _ZERO_TOL_)) {
      transformation_map = symmetry_map[fg];
    } else {
      double tol = _ZERO_TOL_;
      bool skewed = SYM::isLatticeSkewed(scell.lattice, pcell.dist_nn_min, tol);
      for (uint at = 0; at < natoms; at++) {
        atom_trans = symmetry_map[fg][at];
        xvector<double> fpos = scell.c2f * (scell.atoms[atom_trans].cpos - transl);
        for (uint at_map = 0; at_map < natoms; at_map++) {
          if (SYM::FPOSMatch(fpos, scell.atoms[at_map].fpos,
                scell.lattice, scell.f2c, skewed, scell.sym_eps)) { //DX20190619 - lattice and f2c as input
            transformation_map[at] = at_map;
            at_map = natoms;
          }
        }
      }
    }
    return transformation_map;
  }
  // END Inequivalent distortions functions

  //getHigherOrderDistortions/////////////////////////////////////////////////
  vector<_ineq_distortions> ClusterSet::getHigherOrderDistortions() {
    vector<_ineq_distortions> ineq_dists;
    if (order != 4) return ineq_dists;  // Not implemented for higher order and not necessary for lower
    ineq_dists.resize(pcell.iatoms.size());
    vector<vector<int> > dists(6, vector<int>(1));
    for (uint i = 0; i < 6; i++) dists[i][0] = i;
    for (uint iat = 0; iat < pcell.iatoms.size(); iat++) {
      int at = pc2scMap[pcell.iatoms[iat][0]];
      vector<int> atoms(1, at);
      ineq_dists[iat].atoms = atoms; 
      int cl = 0;
      for (cl = 0; cl < (int) clusters.size(); cl++) {
        int a = 0;
        for (a = 0; a < order; a++) {
          if (clusters[cl].atoms[a] != at) break;
        }
        if (a == order) break;
      }
      vector<int> clusters(1, cl);
      ineq_dists[iat].clusters = clusters;
      for (uint td = 0; td < dists.size(); td++) {
        int eq = -1;
        for (uint fg = 0; fg < pcell.fgroup.size(); fg++) {
          eq = equivalenceDistortions(pcell.fgroup[fg].Uc, dists[td], ineq_dists[iat].distortions, atoms);
          if (eq > -1) {
            appendDistortion(ineq_dists[iat], dists[td], eq, fg);
            fg = pcell.fgroup.size();
          }
        }
        if (eq == -1) appendDistortion(ineq_dists[iat], dists[td]);
      }
    }
    return ineq_dists;
  }

}  // namespace apl


/*************************** LINEAR COMBINATIONS ****************************/

namespace apl {

  //getLinearCombinations/////////////////////////////////////////////////////
  // Obtains the linear combinations of the IFCs. See aflow.apl.h for a
  // description of the _linearCombinations struct.
  vector<_linearCombinations> ClusterSet::getLinearCombinations() {
    vector<_linearCombinations> lincombs(ineq_clusters.size());
    for (uint ineq = 0; ineq < ineq_clusters.size(); ineq++) {
      _linearCombinations lcomb;
      // Build coefficient matrix
      vector<vector<int> > symops = getInvariantSymOps(getCluster(ineq_clusters[ineq][0]));
      vector<vector<double> > coeff_mat = buildCoefficientMatrix(symops);
      if (coeff_mat.size() > 0) {
        // Determine the reduced row echelon form of the coefficiebt
        // matrix and extract the linearly dependent and independnet
        // coefficients.
        vector<vector<double> > rref = getRREF(coeff_mat);
        int shift = 0;
        for (uint r = 0; r < rref.size(); r++) {
          while ((abs(rref[r][r+shift]) < _ZERO_TOL_) &&
              (r + shift < rref[r].size())) {
            lcomb.independent.push_back(r + shift);
            shift++;
          }
          bool allzero = true;
          vector<double> coefficients;
          vector<int> indices;
          for (uint c = r + shift; c < rref[r].size(); c++) {
            if (abs(rref[r][c]) > _ZERO_TOL_) {
              allzero = false;
              if (c != r + shift) {
                // Write each linearly dependent IFC as a linear
                // combination of the independent ones
                coefficients.push_back(-rref[r][c]);
                indices.push_back(c);
              }
            }
          }
          if (allzero) {
            r = rref.size();
          } else {
            lcomb.dependent.push_back(r + shift);
            lcomb.indices.push_back(indices);
            lcomb.coefficients.push_back(coefficients);
          }
        }
        // Create a map that points the linearly independent IFCs to the IFCs
        // that depend on them. This is used for the symmetrization of the
        // anharmonic IFCs.
        lcomb.indep2depMap.resize(lcomb.independent.size());
        for (uint i = 0; i < lcomb.indices.size(); i++) {
          for (uint j = 0; j < lcomb.indices[i].size(); j++) {
            for (uint ind = 0; ind < lcomb.independent.size(); ind++) {
              if (lcomb.indices[i][j] == lcomb.independent[ind]) {
                lcomb.indep2depMap[ind].push_back(lcomb.dependent[i]);
              }
            }
          }
        }
      } else {
        for (int crt = 0; crt < nifcs; crt++) {
          lcomb.independent.push_back(crt);
        }
        lcomb.indep2depMap.resize(lcomb.independent.size());
      }
      lincombs[ineq] = lcomb;
    }
    return lincombs;
  }

  //getSymOps/////////////////////////////////////////////////////////////////
  // Gets all the symmetry operations permutation that transform a cluster of
  // atoms into itself or a permutation of itself. This algorithm cannot use
  // ClusterSet::comparePermutations because it needs to find all permutations
  // and not just the first.
  vector<vector<int> > ClusterSet::getInvariantSymOps(const _cluster& ineq_clst) {
    vector<vector<int> > sym;
    vector<int> cluster_transformed(order), cluster_transl(order);
    for (uint fg = 0; fg < pcell.fgroup.size(); fg++) {
      for (int at = 0; at < order; at++) {
        cluster_transformed[at] = symmetry_map[fg][ineq_clst.atoms[at]];
      }
      vector<int> previous_permutation;
      int at = 0;
      for (uint iperm = 0; iperm < permutations.size(); iperm++) {
        bool permutation_found = true;
        vector<int> permutation = permutations[iperm];
        if ((iperm ==  0) || (permutation[0] != previous_permutation[0])) {
          cluster_transl = translateToPcell(cluster_transformed, at);
          at++;
        }
        // Skip identity because it has only zeros in the coefficient matrix 
        if ((fg > 0) || (iperm > 0)) {
          vector<int> cluster_permut(permutation.size());
          for (uint i = 0; i < cluster_transl.size(); i++) {
            cluster_permut[i] = cluster_transl[permutation[i]];
            if (cluster_permut[i] != ineq_clst.atoms[i]) {
              permutation_found = false;
              i = cluster_transl.size();
            }
          }
          if (permutation_found) {
            vector<int> fgperm(2);
            fgperm[0] = fg;
            fgperm[1] = iperm;
            sym.push_back(fgperm);
          }
        }
        previous_permutation = permutation;
      }
    }
    return sym;
  }

  //buildCoefficientMatrix////////////////////////////////////////////////////
  // Builds the coefficient matrix that will be used to obtain the RREF.
  vector<vector<double> > ClusterSet::buildCoefficientMatrix(const vector<vector<int> >& symops) {
    vector<vector<double> > coeff_mat;
    vector<double> row;
    bool append = true;
    int rw = 0, cl = 0, perm = 0, fg = 0, p = 0;
    // Build Cartesian indices - faster and easier to read than running
    // an xcombo within an xcombo
    vector<vector<int> > cart_indices;
    aurostd::xcombos crt_ind(3, order, 'E', true);
    while (crt_ind.increment()) {
      cart_indices.push_back(crt_ind.getCombo());
    }

    for (uint sop = 0; sop < symops.size(); sop++) {
      fg = symops[sop][0];
      perm = symops[sop][1];
      for (int r = 0; r < nifcs; r++) {
        row.resize(nifcs, 0);
        append = false;
        for (int c = 0; c < nifcs; c++) {
          row[c] = 1.0;
          for (int o = 0; o < order; o++) {
            rw = cart_indices[r][o] + 1;
            p = permutations[perm][o];
            cl = cart_indices[c][p] + 1;
            row[c] *= pcell.fgroup[fg].Uc[rw][cl];
            if (abs(row[c]) < _ZERO_TOL_) {
              row[c] = 0.0;
              o = order;
            }
          }
          if (r == c) {
            row[c] -= 1.0;
          }
          if (abs(row[c]) < _ZERO_TOL_) {
            row[c] = 0.0;
          } else {
            append = true;
          }
        }
        if (append) {
          coeff_mat.push_back(row);
        }
      }
    }
    return coeff_mat;
  }

  //getRREF///////////////////////////////////////////////////////////////////
  // Obtains the reduced row echelon form (RREF) of the coefficient matrix
  // using Gaussian elimination.
  vector<vector<double> > ClusterSet::getRREF(vector<vector<double> > mat) {
    uint currentrow = 0, nrows = 0, ncols = 0;
    double swap = 0.0;
    nrows = mat.size();
    ncols = mat[0].size();
    currentrow = 0;
    for (uint col = 0; col < ncols; col++) { 
      for (uint row = 0; row < nrows; row++) {
        if (abs(mat[row][col]) < _ZERO_TOL_) {
          mat[row][col] = 0.0;
        }
      }
      for (uint row = currentrow + 1; row < nrows; row++) {
        if (abs(mat[row][col]) - abs(mat[currentrow][col]) > _ZERO_TOL_) {
          for (uint c = col; c < ncols; c++) {
            swap = mat[currentrow][c];
            mat[currentrow][c] = mat[row][c];
            mat[row][c] = swap;
          }
        }
      }
      if (abs(mat[currentrow][col]) > _ZERO_TOL_) { 
        for (uint i = ncols - 1; i > col; i--) {
          mat[currentrow][i] /= mat[currentrow][col];
        }
        mat[currentrow][col] = 1.0;
        for (uint row = 0; row < nrows; row++) {
          if (row != currentrow) {
            for (uint i = ncols - 1; i > col; i--) {
              mat[row][i] -= mat[row][col] * mat[currentrow][i]/mat[currentrow][col];
            }
            mat[row][col] = 0.0;
          }
        }
        if (currentrow < nrows - 1) {
          currentrow++;
        }
      }
    }
    return mat;
  }

}  // namespace apl

/********************************* FILE I/O *********************************/

namespace apl {

  // BEGIN Write files
  //writeClusterSetToFile/////////////////////////////////////////////////////
  // Writes the ClusterSet object to an XML file.
  void ClusterSet::writeClusterSetToFile(const string& filename) {
    string message = "Writing ClusterSet to file " + filename + ".";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, directory, *p_FileMESSAGE, *p_oss);
    std::stringstream output;
    // Header
    output << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << std::endl;
    output << "<cluster_set>" << std::endl;

    output << writeParameters();
    output << writeClusters(clusters);
    output << writeInequivalentClusters();
    output << writeInequivalentDistortions();
    if (order == 4) output << writeHigherOrderDistortions();
    output << "</cluster_set>" << std::endl;
    aurostd::stringstream2file(output, filename);
    if (!aurostd::FileExist(filename)) {
      message = "Could not write ClusterSet to file.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
  }

  //writeParameters///////////////////////////////////////////////////////////
  // Writes the calculation parameters and minimal structure information of
  // the XML file.
  string ClusterSet::writeParameters() {
    std::stringstream parameters;
    string tab = " ";

    // Info about calculation run
    parameters << tab << "<generator>" << std::endl;
    parameters << tab << tab << "<i name=\"aflow_version\" type=\"string\">" << AFLOW_VERSION << "</i>" << std::endl;
    string time = aflow_get_time_string();
    if (time[time.size() - 1] == '\n') time.erase(time.size() - 1);
    parameters << tab << tab << "<i name=\"date\" type=\"string\">" << time << "</i>" << std::endl;
    parameters << tab << tab << "<i name=\"checksum\" file=\"" << _AFLOWIN_;
    parameters << "\" type=\"" << APL_CHECKSUM_ALGO << "\">" << std::hex << aurostd::getFileCheckSum(directory + "/" + _AFLOWIN_, APL_CHECKSUM_ALGO);
    parameters.unsetf(std::ios::hex);  //ME20190125 - Remove hexadecimal formatting
    parameters  << "</i>" << std::endl;
    parameters << tab << "</generator>" << std::endl;

    // Parameters
    parameters << tab << "<order>" << order << "</order>" << std::endl;
    parameters << tab << "<cutoff units=\"Angstrom\">" << std::setprecision(12) << cutoff << "</cutoff>" << std::endl;

    // Structure
    parameters << tab << "<structure units=\"Angstrom\" cs=\"fractional\">" << std::endl;
    parameters << tab << tab << "<varray name=\"pcell lattice\">" << std::endl;
    for (int i = 1; i < 4; i++) {
      parameters << tab << tab << tab << "<v>";
      for (int j = 1; j < 4; j++) {
        parameters << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
        parameters << std::setprecision(8);
        parameters << std::setw(15) << pcell.lattice[i][j];
      }
      parameters << "</v>" << std::endl;
    }
    parameters << tab << tab << "</varray>" << std::endl;
    parameters << tab << tab << "<varray name=\"positions\">" << std::endl;
    for (uint i = 0; i < pcell.atoms.size(); i++) {
      int t = pcell.atoms[i].type;
      parameters << tab << tab << tab << "<v species=\"" << pcell.species[t] << "\">";
      for (int j = 1; j < 4; j++) {
        parameters << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
        parameters << std::setprecision(8);
        parameters << std::setw(15) << pcell.atoms[i].fpos[j];
      }
      parameters << "</v>" << std::endl;
    }
    parameters << tab << tab << "</varray>" << std::endl;
    parameters << tab << tab << "<varray name=\"supercell\">" << std::endl;
    parameters << tab << tab << tab << "<v>";
    for (int i = 1; i < 4; i++) {
      parameters << tab << sc_dim[i];
    }
    parameters << "</v>" << std::endl;
    parameters << tab << tab << "</varray>" << std::endl;
    parameters << tab << "</structure>" << std::endl;
    return parameters.str();
  }

  //writeClusters/////////////////////////////////////////////////////////////
  // Writes the _cluster objects inside a set of clusters.
  string ClusterSet::writeClusters(const vector<_cluster>& clusters) {
    string tab = " ";
    stringstream clst;
    clst << "<clusters>" << std::endl;
    for (uint c = 0; c < clusters.size(); c++) {
      clst << tab << "<cluster>" << std::endl;
      clst << tab << tab << "<atoms>";
      for (int at = 0; at < order; at++) {
        clst << " " << clusters[c].atoms[at];
      }
      clst << " </atoms>" << std::endl;
      clst << tab << tab << "<fgroup>";
      clst << clusters[c].fgroup << "</fgroup>" << std::endl;
      clst << tab << tab << "<permutation>";
      clst << clusters[c].permutation << "</permutation>" << std::endl;
      clst << tab << "</cluster>" << std::endl;
    }
    clst << "</clusters>" << std::endl;
    return clst.str();
  }

  //writeInequivalentClusters/////////////////////////////////////////////////
  // Writes the inequivalent clusters and the linear combinations.
  string ClusterSet::writeInequivalentClusters() {
    string tab = " ";
    stringstream clusters;
    clusters << tab << "<inequivalent_clusters>" << std::endl;
    for (uint i = 0; i < ineq_clusters.size(); i++) {
      clusters << tab << tab << "<ineq_cluster ID=\"" << i << "\">" << std::endl;
      clusters << tab << tab << tab << "<varray name=\"clusters\">"
        << aurostd::joinWDelimiter(ineq_clusters[i], " ")
        << "</varray>" << std::endl;
      clusters << writeLinearCombinations(linear_combinations[i]);
      clusters << tab << tab << "</ineq_cluster>" << std::endl;
    }
    clusters << tab << "</inequivalent_clusters>" << std::endl;
    return clusters.str();
  }

  //writeLinearCombinations///////////////////////////////////////////////////
  // Writes the _linearCombinations object.
  string ClusterSet::writeLinearCombinations(const _linearCombinations& lcombs) {
    string tab = " ";
    string basetab = tab + tab + tab;
    stringstream lc;
    lc << basetab << "<linear_combinations>" << std::endl;

    lc << basetab << tab << "<independent>";
    for (uint i = 0; i < lcombs.independent.size(); i++) {
      lc << " " << lcombs.independent[i];
    }
    lc << " </independent>" << std::endl;

    lc << basetab << tab << "<dependent>";
    for (uint i = 0; i < lcombs.dependent.size(); i++) {
      lc << " " << lcombs.dependent[i];
    }
    lc << " </dependent>" << std::endl;

    lc << basetab << tab << "<varray name=\"indices\">" << std::endl;
    for (uint i = 0; i < lcombs.indices.size(); i++) {
      lc << basetab << tab << tab << "<v>";
      for (uint j = 0; j < lcombs.indices[i].size(); j++) {
        lc << " " << lcombs.indices[i][j];
      }
      lc << " </v>" << std::endl;
    }
    lc << basetab << tab << "</varray>" << std::endl;

    lc << basetab << tab << "<varray name=\"coefficients\">" << std::endl;
    for (uint i = 0; i < lcombs.coefficients.size(); i++) {
      lc << basetab << tab << tab << "<v>";
      for (uint j = 0; j < lcombs.coefficients[i].size(); j++) {
        lc << " " << lcombs.coefficients[i][j];
      }
      lc << " </v>" << std::endl;
    }
    lc << basetab << tab << "</varray>" << std::endl;

    lc << basetab << tab << "<varray name=\"indep2depMap\">" << std::endl;
    for (uint i = 0; i < lcombs.indep2depMap.size(); i++) {
      lc << basetab << tab << tab << "<v>";
      for (uint j = 0; j < lcombs.indep2depMap[i].size(); j++) {
        lc << " " << lcombs.indep2depMap[i][j];
      }
      lc << " </v>" << std::endl;
    }
    lc << basetab << tab << "</varray>" << std::endl;

    lc << basetab << "</linear_combinations>" << std::endl;
    return lc.str();
  }

  //writeInequivalentDistortions//////////////////////////////////////////////
  // Writes the inequivalent distortions.
  string ClusterSet::writeInequivalentDistortions() {
    string tab = " ";
    stringstream distortions;
    distortions << tab << "<inequivalent_distortions>" << std::endl;
    for (uint i = 0; i < ineq_distortions.size(); i++) {
      distortions << tab << tab << "<ineq_dist ID=\"" << i << "\">" << std::endl;
      distortions << writeIneqDist(ineq_distortions[i]);
      distortions << tab << tab << "</ineq_dist>" << std::endl;
    } 
    distortions << tab << "</inequivalent_distortions>" << std::endl;
    return distortions.str();
  }

  //writeIneqDist/////////////////////////////////////////////////////////////
  // Writes a single _ineq_distortions object.
  string ClusterSet::writeIneqDist(const _ineq_distortions& ineq_dist) {
    string tab = " ";
    string basetab = tab + tab + tab;
    stringstream idist;

    idist << basetab << "<atoms>";
    for (uint at = 0; at < ineq_dist.atoms.size(); at++) {
      idist << " " << ineq_dist.atoms[at];
    }
    idist << " </atoms>" << std::endl;

    idist << basetab << "<clusters>";
    for (uint c = 0; c < ineq_dist.clusters.size(); c++) {
      idist << " " << ineq_dist.clusters[c];
    }
    idist << " </clusters>" << std::endl;

    idist << basetab << "<distortions>" << std::endl;
    for (uint d = 0; d < ineq_dist.distortions.size(); d++) {
      idist << basetab << tab << "<varray ID=\"" << d << "\">" << std::endl;
      for (uint i = 0; i < ineq_dist.distortions[d].size(); i++) {
        idist << basetab << tab << tab << "<v>";
        for (uint j = 0; j < ineq_dist.distortions[d][i].size(); j++) {
          idist << " " << ineq_dist.distortions[d][i][j];
        }
        idist << " </v>" << std::endl;
      }
      idist << basetab << tab << "</varray>" << std::endl;
    }
    idist << basetab << "</distortions>" << std::endl;

    idist << basetab << "<rotations>" << std::endl;
    idist << basetab << tab << "<varray>" << std::endl;
    for (uint r = 0; r < ineq_dist.rotations.size(); r++) {
      idist << basetab << tab << tab << "<v>";
      for (uint i = 0; i < ineq_dist.rotations[r].size(); i++) {
        idist << " " << ineq_dist.rotations[r][i];
      }
      idist << " </v>" << std::endl;
    }
    idist << basetab << tab << "</varray>" << std::endl;
    idist << basetab << "</rotations>" << std::endl;

    idist << basetab << "<transformation_maps>" << std::endl;
    for (uint t = 0; t < ineq_dist.transformation_maps.size(); t++) {
      idist << basetab << tab << "<varray ID=\"" << t << "\">" << std::endl;
      for (uint i = 0; i < ineq_dist.transformation_maps[t].size(); i++) {
        idist << basetab << tab << tab << "<v>";
        for (uint j = 0; j < ineq_dist.transformation_maps[t][i].size(); j++) {
          idist << " " << ineq_dist.transformation_maps[t][i][j];
        }
        idist << " </v>" << std::endl;
      }
      idist << basetab << tab << "</varray>" << std::endl;
    }
    idist << basetab << "</transformation_maps>" << std::endl;

    return idist.str();
  }

  string ClusterSet::writeHigherOrderDistortions() {
    string tab = " ";
    stringstream distortions;
    distortions << tab << "<higher_order_distortions>" << std::endl;
    for (uint i = 0; i < higher_order_ineq_distortions.size(); i++) {
      distortions << tab << tab << "<ineq_dist ID=\"" << i << "\">" << std::endl;
      distortions << writeIneqDist(higher_order_ineq_distortions[i]);
      distortions << tab << tab << "</ineq_dist>" << std::endl;
    }
    distortions << tab << "</higher_order_distortions>" << std::endl;
    return distortions.str();
  }

  // END Write files

  // BEGIN Read files

  //readClusterSetFromFile////////////////////////////////////////////////////
  // Reads a ClusterSet from an XML file.
  void ClusterSet::readClusterSetFromFile(const string& filename) {
    // Open file and handle exceptions
    stringstream message;

    if (!aurostd::EFileExist(filename) && !aurostd::FileExist(filename)) {
      message << "Could not open file " << filename << ". File not found.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
    vector<string> vlines;
    aurostd::efile2vectorstring(filename, vlines);
    if (vlines.size() == 0) {
      message << "Cannot open file " << filename << ". File empty or corrupt.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    }

    // Start reading
    uint line_count = 0;
    string line = vlines[line_count++];

    // Check that this is a valid xml file
    if (line.find("xml") == string::npos) {
      message << "File is not a valid xml file.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
    }

    //Check if xml file can be used to read the ClusterSet object
    if (checkCompatibility(line_count, vlines)) {
      clusters = readClusters(line_count, vlines);
      readInequivalentClusters(line_count, vlines);
      readInequivalentDistortions(line_count, vlines);
      if (order == 4) readHigherOrderDistortions(line_count, vlines);
    } else {
      message << "The settings in the hibernate file and the aflow.in file are incompatible.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    // Necessary for AAPL calculations; since these calculations are
    // quick, they are not stored in the output file
    distortion_vectors = getCartesianDistortionVectors();
    permutations = getPermutations(order);
    nifcs = aurostd::powint(3, order);
  }

  //checkCompatibility////////////////////////////////////////////////////////
  // Checks if the hibernate XML file is compatible with the aflow.in file.
  // If the checksum in the XML file is the same as the checksum of the
  // aflow.in file, then the parameters are the same. If not, the function
  // checks if the parameters relevant for the ClusterSet (supercell, order,
  // cutoff) are the same. This prevents the ClusterSet from being
  // recalculated when only post-processing parameters are changed.
  bool ClusterSet::checkCompatibility(uint& line_count, const vector<string>& vlines) {
    string line = "";
    std::stringstream message;
    bool compatible = true;
    int t = 0;
    vector<string> tokens;
    uint vsize = vlines.size();

    // Compare checksum
    while (true) {
      if (line_count == vsize) {
        message << "Checksum not found in hibernate file.";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
      }
      line = vlines[line_count++];
      if (line.find("checksum") != string::npos) {
        break;
      }
    }

    t = line.find_first_of(">") + 1;
    aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
    if (strtoul(tokens[0].c_str(), NULL, 16) != aurostd::getFileCheckSum(directory + "/" + _AFLOWIN_, APL_CHECKSUM_ALGO)) {  //ME20190219
      message << "The " << _AFLOWIN_ << " file has been changed from the hibernated state. ";

      tokens.clear();

      // Compare order
      while (compatible) {
        if (line_count == vsize) {
          message << "Could not find order tag. ";
          compatible = false;
        }
        line = vlines[line_count++];
        if (line.find("order") != string::npos) {
          t = line.find_first_of(">") + 1;
          aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
          int ord = aurostd::string2utype<int>(tokens[0]);
          tokens.clear();
          if (ord != order) {
            message << "Hibernate file and aflow.in have different order. ";
            compatible = false;
          }
          break;
        }
      }

      // Compare cutoff
      while (compatible) {
        if (line_count == vsize) {
          message << "Could not find cutoff tag. ";
          compatible = false;
        }
        line = vlines[line_count++];
        if (line.find("cutoff") != string::npos) {
          t = line.find_first_of(">") + 1;
          aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
          double cut = aurostd::string2utype<double>(tokens[0]);
          tokens.clear();
          if (abs(cut - cutoff) > _ZERO_TOL_) {
            message << "Hibernate file and aflow.in have different cutoff. ";
            compatible = false;
          }
          break;
        }
      }

      // Compare supercells
      //// Lattice
      while (compatible) {
        if (line_count == vsize) {
          message << "Could not find structure tag. ";
          compatible = false;
        }
        line = vlines[line_count++];
        if (line.find("structure") != string::npos) {
          break;
        }
      }

      while (compatible) {
        if (line_count == vsize) {
          message << "Could not find primitive lattice vectors. ";
          compatible = false;
        }
        line = vlines[line_count++];
        if (line.find("varray name=\"pcell lattice\"") != string::npos) {
          xmatrix<double> latt(3, 3);
          for (int i = 1; i < 4; i++) {
            line = vlines[line_count++];
            if (line_count == vsize) {
              message << "pcell lattice tag is corrupt. ";
            }
            t = line.find_first_of(">") + 1;
            aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
            for (int j = 1; j < 3; j++) {
              latt(i, j) = aurostd::string2utype<double>(tokens[j - 1]);
            }
            tokens.clear();
          }
          line = vlines[line_count++];
          if (line.find("</varray>") == string::npos) {
            message << "pcell lattice tag is corrupt. ";
            compatible = false;
          } else {
            break;
          }
          if (latt != pcell.lattice) {
            message << "Hibernate file and aflow.in do not have the same lattice. ";
            compatible = false;
          }
        }
      }

      //// Atomic positions
      while (compatible) {
        if (line_count == vsize) {
          message << "Could not find atomic positions. ";
          compatible = false;
        }
        line = vlines[line_count++];
        if (line.find("varray name=\"positions\"") != string::npos) {
          // First check for tag corruption and extract everything
          vector<string> species;
          vector<xvector<double> > positions;
          while (compatible) {
            if (line_count == vsize) {
              message << "positions tag is corrupt. ";
              compatible = false;
            }
            line = vlines[line_count++];
            if (line.find("species=\"") != string::npos) {
              // Extract species
              t = line.find_first_of("\"") + 1;
              aurostd::string2tokens(line.substr(t, line.find_last_of("\"") - t), tokens, string(" "));
              species.push_back(tokens[0]);
              tokens.clear();
              // Extract positions
              t = line.find_first_of(">") + 1;
              aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
              xvector<double> fpos(3);
              for (int i = 1; i < 4; i++) {
                fpos(i) = aurostd::string2utype<double>(tokens[i-1]);
              }
              positions.push_back(fpos);
              tokens.clear();
            } else if (line.find("/varray") != string::npos) {
              break;
            } else {
              message << "positions tag is corrupt. ";
              compatible = false;
            }
          }
          // Now compare with the primitive cell from the aflow.in file
          uint nspecies = species.size();
          if (nspecies == pcell.atoms.size()) {
            for (uint sp = 0; sp < nspecies; sp++) {
              int type = pcell.atoms[sp].type;
              string spec = pcell.species[type];
              xvector<double> pos = pcell.atoms[sp].fpos;
              if (species[sp] != spec) {
                message << "The structures in the hibernate file and aflow.in ";
                message << "have different species. ";
                compatible = false;
                sp = nspecies;
              } else if (positions[sp] != pos) {
                message << "The structures in the hibernate file and aflow.in ";
                message << "have different atomic positions.";
                compatible = false;
                sp = nspecies;
              }
            }
          } else {
            message << "The structures in the hibernate file and in aflow.in ";
            message << "do not have the same number of atoms. ";
            compatible = false;
          }
          break;
        }
      }

      //// Supercell dimensions
      while (compatible) {
        if (line_count == vsize) {
          message << "Could not find supercell dimensions. ";
          compatible = false;
        }
        line = vlines[line_count++];
        if (line.find("varray name=\"supercell\"") != string::npos) {
          line = vlines[line_count++];
          t = line.find_first_of(">") + 1;
          aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
          xvector<int> sc(3);
          for (int i = 1; i < 4; i++) {
            sc(i) = aurostd::string2utype<double>(tokens[i-1]);
          }
          if (sc == sc_dim) {
            break;
          } else {
            message << "The supercells in the hibernate file and in aflow.in ";
            message << "have different dimensions. ";
            compatible = false;
          }
        }
      }

      if (compatible) {
        message << "The relevant settings appear to be the same, ";
        message << "so the ClusterSet will not be determined again. ";
        message << "Make sure that the changes do not impact the IFCs.";
      } else {
        message << "The ClusterSet needs to be determined again.";
      }
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, directory, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
    }
    return compatible;
  }


  //readInequivalentClusters//////////////////////////////////////////////////
  // Reads the inequivalent clusters.
  void ClusterSet::readInequivalentClusters(uint& line_count,
      const vector<string>& vlines) {
    string line = "", message = "";
    vector<string> tokens;
    int t = 0;
    uint vsize = vlines.size();

    // Find inequivalent_clusters tag
    while (true) {
      if (line_count == vsize) {
        message = "inequivalent_clusters tag not found";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
      }
      line = vlines[line_count++];
      if (line.find("inequivalent_clusters") != string::npos) {
        break;
      }
    }

    // Read inequivalent clusters
    while (line.find("/inequivalent_clusters") == string::npos) {
      if (line_count == vsize) {
        message = "inequivalent_clusters tag incomplete";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
      }
      line = vlines[line_count++];
      if (line.find("ineq_cluster") != string::npos) {
        while (line.find("/ineq_cluster") == string::npos) {
          if (line_count == vsize) {
            message = "ineq_cluster tag incomplete";
            throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
          }
          line = vlines[line_count++];
          if (line.find("varray name=\"clusters\"") != string::npos) {
            vector<int> clst;
            t = line.find_first_of(">") + 1;
            aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), clst, string(" "));
            ineq_clusters.push_back(clst);
          } else if (line.find("linear_combinations") != string::npos) {
            _linearCombinations lcomb = readLinearCombinations(line_count, vlines);
            linear_combinations.push_back(lcomb);
          }
        }
      }
    }
  }

  //readClusters//////////////////////////////////////////////////////////////
  // Reads a set of _cluster objects for the inequivalent clusters.
  vector<_cluster> ClusterSet::readClusters(uint& line_count,
      const vector<string>& vlines) {
    string message = "", line = "";
    vector<_cluster> clusters;
    vector<string> tokens;
    int t = 0;
    uint vsize = vlines.size();

    while (line.find("/clusters") == string::npos) {
      if (line_count == vsize) {
        message = "clusters tag incomplete";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
      } 
      line = vlines[line_count++];
      if (line.find("<cluster>") != string::npos) {
        _cluster clst;
        while(line.find("</cluster>") == string::npos) {
          if (line_count == vsize) {
            message = "cluster tag incomplete";
            throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
          }
          line = vlines[line_count++];
          if (line.find("atoms") != string::npos) {
            tokens.clear();
            t = line.find_first_of(">") + 1;
            aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
            for (uint at = 0; at < tokens.size(); at++) {
              clst.atoms.push_back(aurostd::string2utype<int>(tokens[at]));
            }
          } else if (line.find("fgroup") != string::npos) {
            tokens.clear();
            t = line.find_first_of(">") + 1;
            aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
            clst.fgroup = aurostd::string2utype<int>(tokens[0]);
          } else if (line.find("permutation") != string::npos) {
            tokens.clear();
            t = line.find_first_of(">") + 1;
            aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
            clst.permutation = aurostd::string2utype<int>(tokens[0]);
          }
        }
        clusters.push_back(clst);
      }
    }

    return clusters;
  }

  //readLinearCombinations////////////////////////////////////////////////////
  // Reads a single _linearCombinations object.
  _linearCombinations ClusterSet::readLinearCombinations(uint& line_count,
      const vector<string>& vlines) {
    string message = "", line = "";
    _linearCombinations lcombs;
    vector<string> tokens;
    int t = 0;
    uint vsize = vlines.size();

    while (line.find("/linear_combinations") == string::npos) {
      if (line_count == vsize) {
        message = "linear_combinations tag incomplete";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
      }
      line = vlines[line_count++];
      if (line.find("independent") != string::npos) {
        tokens.clear();
        t = line.find_first_of(">") + 1;
        aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        for (uint i = 0; i < tokens.size(); i++) {
          lcombs.independent.push_back(aurostd::string2utype<int>(tokens[i]));
        }
      } else if (line.find("dependent") != string::npos) {
        tokens.clear();
        t = line.find_first_of(">") + 1;
        aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        for (uint i = 0; i < tokens.size(); i++) {
          lcombs.dependent.push_back(aurostd::string2utype<int>(tokens[i]));
        }
      } else if (line.find("indices") != string::npos) {
        while(true) {
          if (line_count == vsize) {
            message = "indices varray incomplete";
            throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
          }
          line = vlines[line_count++];
          if (line.find("/varray") != string::npos) {
            break;
          } else {
            vector<int> indices;
            tokens.clear();
            t = line.find_first_of(">") + 1;
            aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
            for (uint i = 0; i < tokens.size(); i++) {
              indices.push_back(aurostd::string2utype<int>(tokens[i]));
            }
            lcombs.indices.push_back(indices);
          }
        }
      } else if (line.find("coefficients") != string::npos) {
        while(true) {
          if (line_count == vsize) {
            message = "coefficients varray incomplete";
            throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
          }
          line = vlines[line_count++];
          if (line.find("/varray") != string::npos) {
            break;
          } else {
            vector<double> coefficients;
            tokens.clear();
            t = line.find_first_of(">") + 1;
            aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
            for (uint i = 0; i < tokens.size(); i++) {
              coefficients.push_back(aurostd::string2utype<double>(tokens[i]));
            }
            lcombs.coefficients.push_back(coefficients);
          }
        }
      } else if (line.find("indep2depMap") != string::npos) {
        while(true) {
          if (line_count == vsize) {
            message = "indep2depMap varray incomplete";
            throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
          }
          line = vlines[line_count++];
          if (line.find("/varray") != string::npos) {
            break;
          } else {
            vector<int> indep2depMap;
            tokens.clear();
            t = line.find_first_of(">") + 1;
            aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
            for (uint i = 0; i < tokens.size(); i++) {
              indep2depMap.push_back(aurostd::string2utype<int>(tokens[i]));
            }
            lcombs.indep2depMap.push_back(indep2depMap);
          }
        }
      }
    }

    return lcombs;
  }

  //readInequivalentDistortions///////////////////////////////////////////////
  // Reads the inequivalent distortions.
  void ClusterSet::readInequivalentDistortions(uint& line_count,
      const vector<string>& vlines) {
    string line = "", message = "";
    uint vsize = vlines.size();

    while (line.find("/inequivalent_distortions") == string::npos) {
      if (line_count == vsize) {
        message = "inequivalent_distortions tag incomplete";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
      }
      line = vlines[line_count++];
      if (line.find("ineq_dist") != string::npos) {
        _ineq_distortions idist = readIneqDist(line_count, vlines);
        ineq_distortions.push_back(idist);
      }
    }
  }

  _ineq_distortions ClusterSet::readIneqDist(uint& line_count, const vector<string>& vlines) {
    string line = "", message = "";
    uint vsize = vlines.size();
    vector<string> tokens;
    int t = 0;

    _ineq_distortions idist;
    while (line.find("/ineq_dist") == string::npos) {
      if (line_count == vsize) {
        message = "ineq_dist tag incomplete";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
      }
      line = vlines[line_count++];
      if (line.find("atoms") != string::npos) {
        tokens.clear();
        t = line.find_first_of(">") + 1;
        aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        for (uint i = 0; i < tokens.size(); i++) {
          idist.atoms.push_back(aurostd::string2utype<int>(tokens[i]));
        }
      } else if (line.find("clusters") != string::npos) {
        tokens.clear();
        t = line.find_first_of(">") + 1;
        aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        for (uint i = 0; i < tokens.size(); i++) {
          idist.clusters.push_back(aurostd::string2utype<int>(tokens[i]));
        }
      } else if (line.find("distortions") != string::npos) {
        while (line.find("/distortions") == string::npos) {
          if (line_count == vsize) {
            message = "ineq_dist tag incomplete";
            throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
          }
          line = vlines[line_count++];
          if (line.find("varray") != string::npos) {
            vector<vector<int> > distortions;
            while (true) {
              if (line_count == vsize) {
                message = "incomplete distortions varray";
                throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
              }
              line = vlines[line_count++];
              if (line.find("/varray") != string::npos) {
                break;
              } else {
                vector<int> dist;
                tokens.clear();
                t = line.find_first_of(">") + 1;
                aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
                for (uint i = 0; i < tokens.size(); i++) {
                  dist.push_back(aurostd::string2utype<int>(tokens[i]));
                }
                distortions.push_back(dist);
              }
            }
            idist.distortions.push_back(distortions);
          }
        }
      } else if (line.find("rotations") != string::npos) {
        while (line.find("/rotations") == string::npos) {
          if (line_count == vsize) {
            message = "rotations tag incomplete";
            throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
          }
          line = vlines[line_count++];
          if (line.find("varray") != string::npos) {
            while (true) {
              if (line_count == vsize) {
                message = "incomplete rotations varray";
                throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
              }
              line = vlines[line_count++];
              if (line.find("/varray") != string::npos) {
                break;
              } else {
                vector<int> rot;
                tokens.clear();
                t = line.find_first_of(">") + 1;
                aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
                for (uint i = 0; i < tokens.size(); i++) {
                  rot.push_back(aurostd::string2utype<int>(tokens[i]));
                }
                idist.rotations.push_back(rot);
              }
            }
          }
        }
      } else if (line.find("transformation_maps") != string::npos) {
        while (line.find("/transformation_maps") == string::npos) {
          if (line_count == vsize) {
            message = "transformation_maps tag incomplete";
            throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
          }
          line = vlines[line_count++];
          if (line.find("varray") != string::npos) {
            vector<vector<int> > maps;
            while (true) {
              if (line_count == vsize) {
                message = "incomplete distortions varray";
                throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
              }
              line = vlines[line_count++];
              if (line.find("/varray") != string::npos) {
                break;
              } else {
                vector<int> map;
                tokens.clear();
                t = line.find_first_of(">") + 1;
                aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
                for (uint i = 0; i < tokens.size(); i++) {
                  map.push_back(aurostd::string2utype<int>(tokens[i]));
                }
                maps.push_back(map);
              }
            }
            idist.transformation_maps.push_back(maps);
          }
        }
      }
    }
    return idist;
  }

  void ClusterSet::readHigherOrderDistortions(uint& line_count,
      const vector<string>& vlines) {
    string line = "", message = "";
    uint vsize = vlines.size();

    while (line.find("/higher_order_distortions") == string::npos) {
      if (line_count == vsize) {
        message = "higher_order_distortions tag incomplete";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
      }
      line = vlines[line_count++];
      if (line.find("ineq_dist") != string::npos) {
        _ineq_distortions idist = readIneqDist(line_count, vlines);
        higher_order_ineq_distortions.push_back(idist);
      }
    }
  }

  // END Read files

}  // namespace apl

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow MARCO ESTERS - Duke University 2018-2021           *
// *                                                                         *
// ***************************************************************************
