// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

#include "aflow_apl.h"

#define _SYM_ZERO_TOL_LOOSE_ 0.05

//CO START
#define ERROR_VERBOSE false
#define CENTER_PRIM false
#define MAP_VERBOSE false
#define COPY_XSTRUCTURE_FULL false
//CO END

using namespace std;

static const xcomplex<double> iONE(0.0, 1.0);  //ME20200116

namespace apl {

  // ///////////////////////////////////////////////////////////////////////////

  Supercell::Supercell(ostream& oss) : xStream(oss) {
    free();
  }

  Supercell::Supercell(ofstream& mf, ostream& oss) : xStream(mf,oss) {
    free();
  }

  //ME20200102 - Refactored
  Supercell::Supercell(const xstructure& _xstr, ofstream& mf, const string& directory, ostream& oss) : xStream(mf,oss) {  //CO20181226
    free();
    _directory = directory;
    initialize(_xstr);
  }

  //ME20200212 - read from a state file
  Supercell::Supercell(const string& filename, ofstream& mf, const string& directory, ostream& oss) : xStream(mf,oss) {
    free();
    _directory = directory;
    initialize(filename);
  }

  Supercell::Supercell(const Supercell& that) : xStream(*that.getOFStream(),*that.getOSS()) {
    if (this != &that) free();
    copy(that);
  }

  Supercell& Supercell::operator=(const Supercell& that) {
    if (this != &that) free();
    copy(that);
    return *this;
  }

  void Supercell::copy(const Supercell& that) {
    if (this == &that) return;
    xStream::copy(that);
    _directory = that._directory;
    _inStructure = that._inStructure;
    _inStructure_original = that._inStructure_original;  //CO
    _inStructure_light = that._inStructure_light;        //CO
    _pcStructure = that._pcStructure;                    //CO
    _scStructure = that._scStructure;                    //CO
    _scStructure_original = that._scStructure_original;  //CO
    _scStructure_light = that._scStructure_light;        //CO
    _pc2scMap.clear();
    _pc2scMap = that._pc2scMap;
    _sc2pcMap.clear();
    _sc2pcMap = that._sc2pcMap;
    //CO START
    _skew = that._skew;
    _derivative_structure = that._derivative_structure;
    _sym_eps = that._sym_eps;
    //CO END
    _isShellRestricted = that._isShellRestricted;
    _maxShellRadius.clear();
    _maxShellRadius = that._maxShellRadius;
    _isConstructed = that._isConstructed;  //CO
    _initialized = that._initialized;
    phase_vectors = that.phase_vectors;  //ME20200116
  }

  // Destructor
  Supercell::~Supercell() {
    xStream::free();
    free();
  }

  void Supercell::free() {
    _directory = "";
    _inStructure.clear();
    _inStructure_original.clear();
    _inStructure_light.clear();
    _pcStructure.clear();
    _scStructure.clear();
    _scStructure_original.clear();
    _scStructure_light.clear();
    _pc2scMap.clear();
    _sc2pcMap.clear();
    _skew = false;
    _derivative_structure = false;
    _sym_eps = 0.0;
    _isShellRestricted = false;
    _maxShellRadius.clear();
    _isConstructed = false;
    _initialized = false;
    phase_vectors.clear();
    _scStructure.clear();
  }

  void Supercell::clear() {
    free();
  }

  // ///////////////////////////////////////////////////////////////////////////

  void Supercell::initialize(const string& filename, ofstream& mf, ostream& oss) {
    xStream::initialize(mf, oss);
    initialize(filename);
  }

  void Supercell::initialize(const string& filename) {
    string message = "";
    if (!aurostd::EFileExist(filename)) {
      message = "Could not find file " + filename + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }

    // Read data
    vector<string> vlines;
    aurostd::efile2vectorstring(filename, vlines);
    uint nlines = vlines.size();
    uint iline = 0;
    stringstream poscar;

    // Supercell dimensions
    xvector<int> dims;
    for (iline = 0; iline < nlines; iline++) {
      if (aurostd::substring2bool(vlines[iline], "SUPERCELL=")) {
        vector<string> tokens;
        aurostd::string2tokens(vlines[iline], tokens, "=");
        vector<int> vdims;
        aurostd::string2tokens(tokens[1], vdims);
        if (vdims.size() == 3 || vdims.size() == 9) {
          dims = aurostd::vector2xvector(vdims);
          break;
        }
      }
    }
    if (iline == nlines) {
      message = "SUPERCELL tag not found or incomplete.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    }

    // Structure
    xstructure xstr;
    for (iline = 0; iline < nlines; iline++) {
      if (aurostd::substring2bool(vlines[iline], "INPUT_STRUCTURE=START")) {
        while ((++iline < nlines) && !aurostd::substring2bool(vlines[iline], "INPUT_STRUCTURE=STOP")) {
          poscar << vlines[iline] << std::endl;
        }
        if (iline < nlines) {
          xstr = xstructure(poscar);
          break;
        }
      }
    }
    if (iline == nlines) {
      message = "INPUT_STRUCTURE tag not found or incomplete.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    }

    // Build
    initialize(xstr, false);
    build(dims, false);
  }

  //ME20200315 - Added VERBOSE to prevent excessive file output when
  // reading from state file
  void Supercell::initialize(const xstructure& _xstr, ofstream& mf, bool VERBOSE, ostream& oss) {
    xStream::initialize(mf, oss);
    initialize(_xstr, VERBOSE);
  }

  void Supercell::initialize(const xstructure& _xstr, bool VERBOSE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string tmp_dir = _directory;  // Do not delete directory
    clear();
    _directory = tmp_dir;

    //CO20190121 - need to sort by equivalent atoms
    //Discovered with help from Xiaoyu Wang of Eva Zurek's group (UBuffalo)
    xstructure xstr(_xstr);

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " input structure" << std::endl;
      cerr << _inStructure << std::endl;
    }

    xstr.sortAtomsEquivalent(); //CO20190218

    // Copy
    _inStructure = xstr;
    _inStructure.CleanStructure();

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " this is the structure to be analyzed (after sorting via iatoms)" << std::endl;
      cerr << _inStructure << std::endl;
    }

    //CO, DO NOT MODIFY THE STRUCTURE BELOW HERE, THIS INCLUDES RESCALE(), BRINGINCELL(), SHIFORIGINATOM(), etc.
    if (VERBOSE) {
      stringstream message;
      message << "Estimating the symmetry of structure and calculating the input structure. Please be patient."; //primitive cell." << apl::endl; //CO20180216 - we do NOT primitivize unless requested via [VASP_FORCE_OPTION]CONVERT_UNIT_CELL
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
    }
    calculateWholeSymmetry(_inStructure, VERBOSE);
    if(LDEBUG){ //CO20190218
      bool write_inequivalent_flag=_inStructure.write_inequivalent_flag;
      _inStructure.write_inequivalent_flag=true;
      cerr << __AFLOW_FUNC__ << " checking iatoms" << std::endl;
      cerr << _inStructure << std::endl;
      _inStructure.write_inequivalent_flag=write_inequivalent_flag;
    }
    _pcStructure = calculatePrimitiveStructure(); //calculate and store primitive cell

    //store some mutual properties
    //need this for getInputStructureLight()
    //#if CENTER_PRIM
    _inStructure_original = _inStructure;
    LightCopy(_inStructure, _inStructure_light);
    //_inStructure_atoms_original = _inStructure.atoms;
    //#endif
    _skew = SYM::isLatticeSkewed(_inStructure.lattice, _inStructure.dist_nn_min, _inStructure.sym_eps);
    _sym_eps = _inStructure.sym_eps;

    clearSupercell();
    _initialized = true;
  }

  //xStream initializers
  void Supercell::initialize(ostream& oss) {
    xStream::initialize(oss);
  }

  void Supercell::initialize(ofstream& mf, ostream& oss) {
    xStream::initialize(mf, oss);
  }

  // ///////////////////////////////////////////////////////////////////////////

  void Supercell::clearSupercell() {
    _scStructure.info = "not constructed";
    _isConstructed = FALSE;
    _isShellRestricted = FALSE;
    _pc2scMap.clear();
    _sc2pcMap.clear();
    //_skew = FALSE;  //CO, same for _sc and _pc (DO NOT RESET)
    //_derivative_structure = FALSE;  //CO, same for _sc and _pc (DO NOT RESET)
    //_sym_eps = AUROSTD_NAN; //CO, same for _sc and _pc (DO NOT RESET)
    _maxShellRadius.clear();
    _maxShellID = -1;
    phase_vectors.clear();  //ME20200116
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20200315 - Added VERBOSE to prevent excessive file output when reading
  // from state file
  void Supercell::calculateWholeSymmetry(xstructure& xstr, bool VERBOSE) {
    //CO20170804 - we want to append symmetry stuff to ofstream
    _kflags kflags;
    //ME20200330 - Only write for verbose output
    kflags.KBIN_SYMMETRY_PGROUP_WRITE=VERBOSE;
    kflags.KBIN_SYMMETRY_FGROUP_WRITE=VERBOSE;
    kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE=VERBOSE;
    kflags.KBIN_SYMMETRY_PGROUPK_WRITE=VERBOSE;
    kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE=VERBOSE;
    kflags.KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE=VERBOSE;
    kflags.KBIN_SYMMETRY_SGROUP_WRITE=VERBOSE;
    kflags.KBIN_SYMMETRY_IATOMS_WRITE=VERBOSE;
    kflags.KBIN_SYMMETRY_AGROUP_WRITE=VERBOSE;

    //DX CAN REMOVE string options = "";

    xstr.LatticeReduction_avoid = TRUE;
    xstr.sgroup_radius = 8.0;

    //CO20170804 - we want to append symmetry stuff to ofstream
    //if (!pflow::CalculateFullSymmetry(af, xstr))
    if (!pflow::PerformFullSymmetry(xstr,*p_FileMESSAGE,_directory,kflags,VERBOSE,*p_oss)) //CO20181226
    { //CO20200106 - patching for auto-indenting
      string message = "Symmetry routine failed.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  void Supercell::reset() {
    _scStructure = _inStructure;

    _scStructure.atoms.clear();

    for (uint i = 0; i < _scStructure.num_each_type.size(); i++)
      _scStructure.num_each_type[i] = 0;

    for (uint i = 0; i < _scStructure.iatoms.size(); i++)
      _scStructure.iatoms[i].clear();

    _scStructure.agroup.clear();
    _scStructure.fgroup.clear();

    _pc2scMap.clear();
    _sc2pcMap.clear();

    _isConstructed = FALSE;
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20191225
  // Determine the supercell dimensions of the supercell for different methods
  xvector<int> Supercell::determineSupercellDimensions(const aurostd::xoption& opts) {
    stringstream message;
    if (!_initialized) {
      message << "Not initialized.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }

    xvector<int> dims(3);
    string method = opts.getattachedscheme("SUPERCELL::METHOD");
    if (method.empty()) {
      message << "Supercell method empty.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }
    string value = opts.getattachedscheme("SUPERCELL::VALUE");
    if (value.empty()) {
      message << "Supercell value empty.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }
    if (method == "SUPERCELL") {
      vector<int> tokens;
      aurostd::string2tokens(value, tokens, " xX,;");
      dims = aurostd::vector2xvector(tokens);
    } else if (method == "MINATOMS") {
      int minatoms = aurostd::string2utype<int>(value);
      int natoms = (int) _inStructure.atoms.size();
      // ME20200516 - Use the shortest lattice vector as the starting point.
      // The initial sphere needs to be inside the unit cell to catch 1x1x1
      // cells.
      double radius = std::min(std::min(_inStructure.a, _inStructure.b), _inStructure.c)/2.0 - 0.1;
      for ( ; natoms < minatoms; radius += 0.01) {
        dims = LatticeDimensionSphere(_inStructure.lattice, radius);
        natoms = dims[1] * dims[2] * dims[3] * (int) _inStructure.atoms.size();
      }
      if (opts.flag("SUPERCELL::VERBOSE")) {
        message << "Radius=" << aurostd::PaddedPOST(aurostd::utype2string<double>(radius, 3), 4)
          << " supercell=" << dims[1] << "x" << dims[2] << "x" << dims[3]
          << " natoms=" << natoms;
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
      }
    } else if (method == "MINATOMS_UNIFORM") {
      int minatoms = aurostd::string2utype<int>(value);
      int natoms = (int) _inStructure.atoms.size();
      int Ni = 0;
      for (Ni = 1; natoms < minatoms; Ni++) {
        natoms = (int) std::pow(Ni, 3) * _inStructure.atoms.size();
      }
      if (Ni > 1) Ni--; //ME20200521
      dims[1] = Ni; dims[2] = Ni; dims[3] = Ni;
      if (opts.flag("SUPERCELL::VERBOSE")) {
        message << "Ni=" << Ni
          << " supercell=" << Ni << "x" << Ni << "x" << Ni
          << " natoms=" << natoms;
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
      }
    } else if (method == "SHELLS") {
      int shells = aurostd::string2utype<int>(value);
      if (opts.flag("SUPERCELL::VERBOSE")) {
        message << "Searching for suitable cell to handle " << shells << " shells...";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
      }
      bool get_full_shells = true;
      dims = getSupercellDimensionsShell(shells, get_full_shells);
    } else {
      message << "Unknown supercell method " + method + ".";
      aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }
    return dims;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void Supercell::build(aurostd::xoption& opts, bool VERBOSE) {
    opts.flag("SUPERCELL::VERBOSE", VERBOSE);
    xvector<int> dims = determineSupercellDimensions(opts);
    build(dims);
  }

  void Supercell::build(int nx, int ny, int nz, bool VERBOSE) {
    xvector<int> dims(3);
    dims[1] = nx; dims[2] = ny; dims[3] = nz;
    build(dims, VERBOSE);
  }

  void Supercell::build(const xvector<int>& dims, bool VERBOSE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;
    if (!_initialized) {
      message << "Not initialized.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    scell_dim = dims;
    _derivative_structure = ((dims.rows == 9) || !aurostd::identical(dims));
    bool get_full_sym = false; //GETFULLSYMBASIS;  //( GETFULLSYMBASIS || _derivative_structure ); //CO

    // Print info
    if (dims.rows == 3) {
      string dimstring = aurostd::joinWDelimiter(dims, "x");
      _scStructure.info = "Supercell " + dimstring;
      if (VERBOSE) {
        message << "The supercell is going to build as " << dimstring;
        message << " (" << (dims[1] * dims[2] * dims[3] * _inStructure.atoms.size()) << " atoms).";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
      }
    } else {
      xmatrix<int> scmat = aurostd::reshape(dims, 3, 3);
      string dimstring = aurostd::xmatDouble2String(aurostd::xmatrixutype2double(scmat), 0);
      _scStructure.info = "Supercell " + dimstring;
      if (VERBOSE) {
        message << "The supercell is going to build as [" << dimstring << "]";
        message << " (" << (aurostd::det(scmat) * _inStructure.atoms.size()) << " atoms).";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
      }
    }

    if (VERBOSE && _derivative_structure) {
      message << "Derivative structure detected, be patient as we calculate the symmetry of the supercell.";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
    }

    // Get supercell
    //_scStructure = GetSuperCell(_inStructure, dims, _sc2pcMap, _pc2scMap, TRUE, _derivative_structure);  //now gets symmetries too! no need for full_basis (just a check)
    _scStructure = GetSuperCell(_inStructure, dims, _sc2pcMap, _pc2scMap, TRUE, get_full_sym, false, true);  //now gets symmetries too! no need for full_basis (just a check) //CO20190409 - force_supercell_matrix==false as we might have a derivative structure, force_strict_pc2scMap==true because we want to map to true primitive cell, no equivalent atoms

    // Setup output flags
    _scStructure.write_inequivalent_flag = TRUE;

    // OK.
    if (VERBOSE) {
      message << "Supercell successfully created.";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
    }
    _isConstructed = TRUE;

    _scStructure_original = _scStructure;  //COPY EVERYTHING ONCE, will be very slow
    LightCopy(_scStructure, _scStructure_light);
    //_scStructure_atoms_original = _scStructure.atoms;

    //ME20200116 - calculate phase vectors; significantly speeds up post-processing
    calculatePhaseVectors();

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " this is the supercell to be analyzed" << std::endl; //CO20190218
      cerr << _scStructure << std::endl;
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  void Supercell::trimStructure(int n, const xvector<double>& a,
      const xvector<double>& b, const xvector<double>& c, bool constructSymmetry) {
    // Clear some arrays we will rebuild...
    reset();

    // Create lattice of supercell
    _scStructure.lattice(1, 1) = a(1);
    _scStructure.lattice(2, 1) = b(1);
    _scStructure.lattice(3, 1) = c(1);
    _scStructure.lattice(1, 2) = a(2);
    _scStructure.lattice(2, 2) = b(2);
    _scStructure.lattice(3, 2) = c(2);
    _scStructure.lattice(1, 3) = a(3);
    _scStructure.lattice(2, 3) = b(3);
    _scStructure.lattice(3, 3) = c(3);
    _scStructure.FixLattices();

    // Trim from grid p x p x p, where p = 2 * n
    int p = 2 * n;
    _atom atom;
    xvector<double> cshift(3);

    for (uint ia = 0; ia < _inStructure.iatoms.size(); ia++) {
      for (uint iia = 0; iia < _inStructure.iatoms[ia].size(); iia++) {
        // Replicate this atom by given mesh...
        for (int i = 0; i <= p; i++)
          for (int j = 0; j <= p; j++)
            for (int k = 0; k <= p; k++) {
              // Create position of new atoms...
              atom = _inStructure.atoms[_inStructure.iatoms[ia][iia]];
              cshift = (((double)i) * _inStructure.lattice(1) +
                  ((double)j) * _inStructure.lattice(2) +
                  ((double)k) * _inStructure.lattice(3));
              atom.cpos = atom.cpos + cshift;
              atom.fpos = C2F(_scStructure.lattice, atom.cpos);

              // Add only atoms inside the cell
              if (atom.fpos(1) > 1.0 - _FLOAT_TOL_ ||
                  atom.fpos(2) > 1.0 - _FLOAT_TOL_ ||
                  atom.fpos(3) > 1.0 - _FLOAT_TOL_ ||
                  atom.fpos(1) < 0.0 - _FLOAT_TOL_ ||
                  atom.fpos(2) < 0.0 - _FLOAT_TOL_ ||
                  atom.fpos(3) < 0.0 - _FLOAT_TOL_) continue;

              // Increase the number of atoms of this type...
              _scStructure.num_each_type[atom.type]++;

              // Mark this atom as equivalent or not....
              if (_scStructure.iatoms[ia].empty()) {
                atom.equivalent = _scStructure.atoms.size();
                atom.is_inequivalent = TRUE;
              } else {
                atom.equivalent = _scStructure.iatoms[ia][0];
                atom.is_inequivalent = FALSE;
              }

              // Add it to the list of all atoms...
              _scStructure.atoms.push_back(atom);

              // Add its ID number to the list of equivalent atoms of this type...
              _scStructure.iatoms[ia].push_back(_scStructure.atoms.size() - 1);

              // Add its site point group...
              if (constructSymmetry) {
                _scStructure.agroup.push_back(_inStructure.agroup[_inStructure.iatoms[ia][iia]]);
                // We have to correct the Uf for each symop since we have changed the lattice...
                //SC formula - great help!
                vector<_sym_op>::iterator soi = _scStructure.agroup.back().begin();
                for (; soi != _scStructure.agroup.back().end(); soi++)
                  soi->Uf = _scStructure.c2f * soi->Uc * _scStructure.f2c;
              }

              // Update our mapping arrays...
              _sc2pcMap.push_back(_inStructure.iatoms[ia][iia]);
              if (i == 0 && j == 0 && k == 0) _pc2scMap.push_back(_scStructure.atoms.size() - 1);
            }
      }
    }

    // Feed the factor group list (not efficient in this order, but we have all
    // similar operations in order just shifted...
    if (constructSymmetry) {
      for (uint l = 0; l < _inStructure.fgroup.size(); l++) {
        for (int i = 0; i < p; i++)
          for (int j = 0; j < p; j++)
            for (int k = 0; k < p; k++) {
              // Create position of new atoms...
              cshift = (((double)i) * _inStructure.lattice(1) +
                  ((double)j) * _inStructure.lattice(2) +
                  ((double)k) * _inStructure.lattice(3));

              _sym_op symOp = _inStructure.fgroup[l];
              symOp.ctau = symOp.ctau + cshift;
              symOp.ftau = C2F(_scStructure.lattice, symOp.ctau);

              if (symOp.ftau(1) > 1.0 - _FLOAT_TOL_ ||
                  symOp.ftau(2) > 1.0 - _FLOAT_TOL_ ||
                  symOp.ftau(3) > 1.0 - _FLOAT_TOL_ ||
                  symOp.ftau(1) < 0.0 - _FLOAT_TOL_ ||
                  symOp.ftau(2) < 0.0 - _FLOAT_TOL_ ||
                  symOp.ftau(3) < 0.0 - _FLOAT_TOL_) continue;

              // We have to correct the Uf for each symop since we have changed the lattice...
              //SC formula - great help!
              symOp.Uf = _scStructure.c2f * symOp.Uc * _scStructure.f2c;

              _scStructure.fgroup.push_back(symOp);
            }
      }
    }

    // Order atoms in VASP like style
    int start = 0;
    for (int i = 0; i < (int)_scStructure.num_each_type.size(); i++) {
      int end = start + _scStructure.num_each_type[i];
      for (int j = start; j < end - 1; j++) {
        if (aurostd::abs(_scStructure.atoms[j].fpos(1)) < _FLOAT_TOL_) _scStructure.atoms[j].fpos(1) = 0.0;
        if (aurostd::abs(_scStructure.atoms[j].fpos(2)) < _FLOAT_TOL_) _scStructure.atoms[j].fpos(2) = 0.0;
        if (aurostd::abs(_scStructure.atoms[j].fpos(3)) < _FLOAT_TOL_) _scStructure.atoms[j].fpos(3) = 0.0;

        for (int k = j + 1; k < end; k++) {
          if (aurostd::abs(_scStructure.atoms[k].fpos(1)) < _FLOAT_TOL_) _scStructure.atoms[k].fpos(1) = 0.0;
          if (aurostd::abs(_scStructure.atoms[k].fpos(2)) < _FLOAT_TOL_) _scStructure.atoms[k].fpos(2) = 0.0;
          if (aurostd::abs(_scStructure.atoms[k].fpos(3)) < _FLOAT_TOL_) _scStructure.atoms[k].fpos(3) = 0.0;

          if (_scStructure.atoms[k].fpos(1) < _scStructure.atoms[j].fpos(1)) {
            _atom ta = _scStructure.atoms[j];
            _scStructure.atoms[j] = _scStructure.atoms[k];
            _scStructure.atoms[k] = ta;

            for (int l = 0; l < (int)_inStructure.atoms.size(); l++) {
              if (_pc2scMap[l] == j)
                _pc2scMap[l] = k;
              else if (_pc2scMap[l] == k)
                _pc2scMap[l] = j;
            }

            int ti = _sc2pcMap[j];
            _sc2pcMap[j] = _sc2pcMap[k];
            _sc2pcMap[k] = ti;
          } else if (_scStructure.atoms[k].fpos(1) == _scStructure.atoms[j].fpos(1)) {
            if (_scStructure.atoms[k].fpos(2) < _scStructure.atoms[j].fpos(2)) {
              _atom ta = _scStructure.atoms[j];
              _scStructure.atoms[j] = _scStructure.atoms[k];
              _scStructure.atoms[k] = ta;

              for (int l = 0; l < (int)_inStructure.atoms.size(); l++) {
                if (_pc2scMap[l] == j)
                  _pc2scMap[l] = k;
                else if (_pc2scMap[l] == k)
                  _pc2scMap[l] = j;
              }

              int ti = _sc2pcMap[j];
              _sc2pcMap[j] = _sc2pcMap[k];
              _sc2pcMap[k] = ti;
            } else if (_scStructure.atoms[k].fpos(2) == _scStructure.atoms[j].fpos(2)) {
              if (_scStructure.atoms[k].fpos(3) < _scStructure.atoms[j].fpos(3)) {
                _atom ta = _scStructure.atoms[j];
                _scStructure.atoms[j] = _scStructure.atoms[k];
                _scStructure.atoms[k] = ta;

                for (int l = 0; l < (int)_inStructure.atoms.size(); l++) {
                  if (_pc2scMap[l] == j)
                    _pc2scMap[l] = k;
                  else if (_pc2scMap[l] == k)
                    _pc2scMap[l] = j;
                }

                int ti = _sc2pcMap[j];
                _sc2pcMap[j] = _sc2pcMap[k];
                _sc2pcMap[k] = ti;
              }
            }
          }
        }
      }
      start += _scStructure.num_each_type[i];
    }

    // Setup symmetry flags
    if (constructSymmetry) {
      _scStructure.pgroup_xtal_calculated = FALSE;
      _scStructure.pgroupk_calculated = FALSE;
      _scStructure.pgroupk_xtal_calculated = FALSE;
      _scStructure.pgroup_calculated = TRUE;
      _scStructure.fgroup_calculated = TRUE;
      _scStructure.sgroup_calculated = FALSE;
      _scStructure.agroup_calculated = TRUE;
    } else {
      _scStructure.pgroup_xtal_calculated = FALSE;
      _scStructure.pgroupk_calculated = FALSE;
      _scStructure.pgroupk_xtal_calculated = FALSE;
      _scStructure.pgroup_calculated = FALSE;
      _scStructure.fgroup_calculated = FALSE;
      _scStructure.sgroup_calculated = FALSE;
      _scStructure.agroup_calculated = TRUE;
    }

    // Setup output flags
    _scStructure.write_inequivalent_flag = TRUE;

    // Set the information about this construction
    _scStructure.info = "Trimmed supercell";

    // OK.
    //_logger << "Done." << apl::endl;
    _isConstructed = TRUE;
    //cout << _scStructure << std::endl;
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20200116 - rebase to primitive
  // Does not capture rotated conventional cells yet, but does work for AFLOW's
  // standard conventional unit cells.
  bool Supercell::projectToPrimitive() {
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    vector<int> pc2sc, sc2pc;
    xstructure pcell;
    if (!_pcStructure.iatoms_calculated) SYM::CalculateInequivalentAtoms(_pcStructure);
    LightCopy(_pcStructure, pcell);  // No need for symmetry

    // The original structure may be a rotated primitive cell. Transform the
    // primitive cell so that they overlap or else the mapping will not work
    if (_inStructure_original.atoms.size() == pcell.atoms.size()) {
      xmatrix<double> U = trasp(_inStructure_original.scale*_inStructure_original.lattice) * inverse(trasp(pcell.scale*_pcStructure.lattice));
      // For a rotation matrix, trasp(U) * U = I
      if (aurostd::isidentity(trasp(U) * U)) {
        pcell = Rotate(pcell, U);
      }
    }

    bool mapped = getMaps(pcell, _inStructure_original, _scStructure, pc2sc, sc2pc);
    if (!mapped) {
      // When calculating the primitive cell, the positions of the atoms
      // may alternate between x and 1 - x. Test if this is the case here.
      if (LDEBUG) std::cerr << __AFLOW_FUNC__ << " Not mapped succesfully. Try shifting atoms." << std::endl;
      xvector<double> ones(3); ones.set(1.0);
      for (uint at = 0; at < pcell.atoms.size(); at++) {
        pcell.atoms[at].fpos = ones - pcell.atoms[at].fpos;
        pcell.atoms[at].cpos = pcell.f2c * pcell.atoms[at].fpos;
      }
      pcell.BringInCell();
      mapped = getMaps(pcell, _inStructure_original, _scStructure, pc2sc, sc2pc);
    }
    if (mapped) {
      // Sort the iatoms correctly or the non-analytical correction will not work
      // This needs to be done from scratch because the sequence of the atoms may
      // have shifted between the original and the primitive cell.
      xvector<double> fpos(3);
      const vector<vector<int> >& iatoms_pc = _pcStructure.iatoms;
      const vector<vector<int> >& iatoms_oc = _inStructure_original.iatoms;
      uint niatoms_pc = iatoms_pc.size();
      uint niatoms_oc = iatoms_oc.size();
      pcell.iatoms.clear();
      pcell.iatoms.resize(niatoms_pc);
      for (uint iatpc = 0; ((iatpc < niatoms_pc) && mapped); iatpc++) {
        fpos = _inStructure_original.c2f * pcell.atoms[iatoms_pc[iatpc][0]].cpos;
        uint iatoc = 0;
        for (iatoc = 0; iatoc < niatoms_oc; iatoc++) {
          uint i = 0;
          for (i = 0; i < iatoms_oc[iatoc].size(); i++) {
            if (SYM::FPOSMatch(fpos, _inStructure_original.atoms[iatoms_oc[iatoc][i]].fpos,
                  _inStructure_original.lattice, _inStructure_original.f2c, _skew, _sym_eps)) {
              pcell.iatoms[iatoc] = _pcStructure.iatoms[iatpc];
              for (uint j = 0; j < iatoms_pc[iatpc].size(); j++) pcell.atoms[iatoms_pc[iatpc][j]].index_iatoms = iatoc;
              break;
            }
          }
          if (i < iatoms_oc[iatoc].size()) break;
        }
        if (iatoc == niatoms_oc) {
          if (LDEBUG) std::cerr << __AFLOW_FUNC__ << " Did not map iatoms of primitive cell (failed for " << iatpc << ")." << std::endl;
          mapped = false;
        }
      }
    }

    if (mapped) {
      _inStructure = pcell;
      _sc2pcMap = sc2pc;
      _pc2scMap = pc2sc;
      calculatePhaseVectors();
      return true;
    } else {
      return false;
    }
  }

  //ME20200116 - rebase to original
  void Supercell::projectToOriginal() {
    vector<int> pc2sc, sc2pc;
    if (getMaps(_inStructure_original, _inStructure_original, _scStructure, pc2sc, sc2pc)) {
      _inStructure = _inStructure_original;
      _pc2scMap = pc2sc;
      _sc2pcMap = sc2pc;
      calculatePhaseVectors();
    } else {
      // If the mapping fails for the original structure,
      // something went seriously wrong
      string message = "Mapping between original structure and supercell failed.";
      message += " This is likely a bug in the code.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
  }

  //ME20200116
  // Recalculates _sc2pcMap and _pc2scMap based on the projected cell (pcell)
  // using the original cell (ocell).
  // Returns true when mapping is successful.
  bool Supercell::getMaps(const xstructure& pcell, const xstructure& ocell, const xstructure& scell,
      vector<int>& pc2sc, vector<int>& sc2pc) {
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    uint natoms_pc = pcell.atoms.size();
    uint natoms_oc = ocell.atoms.size();
    uint natoms_sc = scell.atoms.size();

    // Check that the cell sizes make sense
    if ((natoms_pc > natoms_oc) || (natoms_oc > natoms_sc)) {
      if (LDEBUG) {
        std::cerr << __AFLOW_FUNC__ << " Cannot map a smaller cell to a larger cell." << std::endl;
      }
      return false;
    }

    pc2sc.clear(); pc2sc.resize(natoms_pc);
    sc2pc.clear(); sc2pc.resize(natoms_sc);
    xvector<double> fpos(3);

    // sc2pcMap - convert scell positions into primitive coordinates
    // and do straightforward fpos matching
    for (uint i = 0; i < natoms_sc; i++) {
      fpos = pcell.c2f * scell.atoms[i].cpos;
      uint j = 0;
      for (j = 0; j < natoms_pc; j++) {
        if (SYM::FPOSMatch(fpos, pcell.atoms[j].fpos, pcell.lattice, pcell.f2c, _skew, _sym_eps)) {
          sc2pc[i] = j;
          break;
        }
      }
      if (j == natoms_pc) {
        if (LDEBUG) {
          std::cerr << __AFLOW_FUNC__ << " sc2pcMap failed for atom " << i << "." << std::endl;
        }
        return false;
      }
    }

    // pc2scMap - more complex. Since the projected cell may not be fully
    // inside original cell, the primitive atoms need to be mapped to the
    // original cell first. Then, pc2scMap can be built using simple cpos
    // matching. Since _pc2scMap may not represent the original structure
    // anymore, it cannot be used to speed up the process.
    for (uint i = 0; i < natoms_pc; i++) {
      uint j = 0;
      fpos = ocell.c2f * pcell.atoms[i].cpos;
      for (j = 0; j < natoms_oc; j++) {
        if (SYM::FPOSMatch(fpos, ocell.atoms[j].fpos, ocell.lattice, ocell.f2c, _skew, _sym_eps)) {
          uint k = 0;
          for (k = 0; k < natoms_sc; k++) {
            if (aurostd::isequal(ocell.atoms[j].cpos, scell.atoms[k].cpos)) break;
          }
          // If scell was created from ocell, this should never happen
          if (k == natoms_sc) {
            if (LDEBUG) {
              std::cerr << __AFLOW_FUNC__ << " pc2scMap failed for atom " << i << "."
                << " Could not map original atom " << j << " to supercell." << std::endl;
            }
            return false;
          }
          pc2sc[i] = k;
          break;
        }
      }
      if (j == natoms_oc) {
        if (LDEBUG) {
          std::cerr << __AFLOW_FUNC__ << " pc2scMap failed for atom " << i << "." << std::endl;
        }
        return false;
      }
    }

    // Mapping successful
    return true;
  }

  // ///////////////////////////////////////////////////////////////////////////

  // ME20200516
  // Calculates the supercell dimensions to fit cut_shell number of coordination
  // shells into the supercell.
  xvector<int> Supercell::getSupercellDimensionsShell(uint cut_shell, bool get_full_shells) {
    _inStructure.grid_atoms.clear();
    uint natoms_iat = _inStructure.iatoms.size();
    vector<bool> shells_in_scell(natoms_iat, false);
    bool all_shells_in_scell = false;

    // Initialize variables that are needed during the loop
    vector<double> distances;
    vector<int> gridatoms_indices;
    xvector<int> dims(3), dims_prev(3);
    xvector<double> cpos_image, a_component, ab_component;
    xmatrix<double> scell_lattice = _inStructure.lattice;
    xmatrix<double> scell_matrix = aurostd::identity((double) 0, 3);
    vector<xvector<double> > l1, l2, l3;
    uint iat = 0, gat = 0, at = 0, ngridatoms = 0, countshell = 0;

    // Set the starting radius to half the length of the smallest lattice vector
    double radius = std::min(std::min(_inStructure.a, _inStructure.b), _inStructure.c)/2.0 - 0.1;
    // Expand the supercell until the coordination shells fit completely into the unit cell
    for ( ; !all_shells_in_scell; radius += 0.01) {
      dims = LatticeDimensionSphere(_inStructure.lattice, radius);
      if (dims == dims_prev) continue;  // No need to check if the supercell dimensions have not changed
      dims_prev = dims;
      scell_matrix[1][1] = dims[1];
      scell_matrix[2][2] = dims[2];
      scell_matrix[3][3] = dims[3];
      scell_lattice = scell_matrix * _inStructure.lattice;
      l1.clear(); l2.clear(); l3.clear();
      for (int i = -1; i <= 1; i++) {
        l1.push_back(i * scell_lattice(1));
        l2.push_back(i * scell_lattice(2));
        l3.push_back(i * scell_lattice(3));
      }
      _inStructure.GenerateGridAtoms(0, dims[1] - 1, 0, dims[2] - 1, 0, dims[3] - 1);
      ngridatoms = _inStructure.grid_atoms.size();
      distances.clear();
      distances.resize(ngridatoms);
      gridatoms_indices.clear();
      gridatoms_indices.resize(ngridatoms);
      for (iat = 0; iat < natoms_iat; iat++) {
        if (!shells_in_scell[iat]) {  // Only check for atoms for which the cell was not big enough yet
          const xvector<double>& cpos_iat = _inStructure.atoms[_inStructure.iatoms[iat][0]].cpos;
          // Determine how many partial coordination shells are inside the cell
          // by finding the number of "unique" interatomic distances.
          for (gat = 0; gat < ngridatoms; gat++) {
            gridatoms_indices[gat] = gat;
            distances[gat] = aurostd::modulus(SYM::minimizeDistanceCartesianMethod(cpos_iat, _inStructure.grid_atoms[gat].cpos, scell_lattice));
          }
          aurostd::sort(distances, gridatoms_indices);
          countshell = 0;
          for (gat = 1; (gat < ngridatoms) && (countshell < cut_shell); gat++) {
            if (distances[gat] > distances[gat - 1] + _APL_SHELL_TOL_) countshell++;
          }
          if (get_full_shells && countshell == cut_shell) {
            // For full shells, every atom of the last coordination shell
            // needs to be inside the supercell. If there is more than one
            // periodic image that has the same distance to the central atom,
            // this condition is not fulfilled.
            bool full_shell = true;
            for (uint i = gat; full_shell && (i < ngridatoms) && (distances[i] < distances[gat] + _APL_SHELL_TOL_); i++) {
              at = gridatoms_indices[i];
              uint image_count = 0;
              for (double nx = -1; full_shell && (nx <= 1); nx++) {
                a_component = _inStructure.grid_atoms[at].cpos + l1[nx + 1];
                for (double ny = -1; full_shell && (ny <= 1); ny++) {
                  ab_component = a_component + l2[ny + 1];
                  for (double nz = -1; full_shell && (nz <= 1); nz++) {
                    cpos_image = ab_component + l3[nz + 1];
                    if (aurostd::modulus(cpos_image - cpos_iat) < distances[gat] + _APL_SHELL_TOL_) {
                      image_count++;
                      full_shell = (image_count < 2);
                    }
                  }
                }
              }
            }
            shells_in_scell[iat] = full_shell;
          } else {
            shells_in_scell[iat] = (countshell == cut_shell);
          }
        }
        if (!shells_in_scell[iat]) break;  // No need to check other atoms - it has to fit for all
      }
      all_shells_in_scell = (iat == natoms_iat);
    }
    _inStructure.grid_atoms.clear();
    return dims;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void Supercell::setupShellRestrictions(int MAX_NN_SHELLS) {
    uint niatoms = _scStructure.iatoms.size();
    uint natoms = _scStructure.atoms.size();
    uint iat = 0, at = 0;
    vector<double> shell_radii(niatoms);
    int countshell = 0;
    vector<double> distances(natoms);
    for (uint iat = 0; iat < niatoms; iat++) {
      countshell = 0;
      const xvector<double>& cpos_iat = _scStructure.atoms[_scStructure.iatoms[iat][0]].cpos;
      for (at = 0; at < natoms; at++) {
        distances[at] = aurostd::modulus(SYM::minimizeDistanceCartesianMethod(_scStructure.atoms[at].cpos, cpos_iat, _scStructure.lattice));
      }
      aurostd::sort(distances);
      for (at = 1; (at < natoms) && (countshell < MAX_NN_SHELLS); at++) {
        if (distances[at] > distances[at - 1] + _APL_SHELL_TOL_) countshell++;
      }
      if (countshell < MAX_NN_SHELLS) {
        stringstream message;
        message << "The supercell is too small to set up shell restrictions for " << MAX_NN_SHELLS << " shells.";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        break;
      }
      shell_radii[iat] = distances[at];
    }
    _isShellRestricted = ((iat == niatoms) && (countshell == MAX_NN_SHELLS));
    if (_isShellRestricted) {
      _maxShellRadius = shell_radii;
      _maxShellID = MAX_NN_SHELLS;
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20190715 - added const to use function with const Supercell &
  bool Supercell::isShellRestricted() const {
    return _isShellRestricted;
  }

  //ME20190715 - added const to use function with const Supercell &
  int Supercell::getMaxShellID() const {
    return _maxShellID;
  }

  // ///////////////////////////////////////////////////////////////////////////

  bool Supercell::isConstructed() {
    return _isConstructed;
  }

  // ///////////////////////////////////////////////////////////////////////////

  const xstructure& Supercell::getSupercellStructure() const {
    return _scStructure;
  }

  // ///////////////////////////////////////////////////////////////////////////

  //CO, only what's necessary (no heavy symmetry copies)
  const xstructure& Supercell::getSupercellStructureLight() const {
    return _scStructure_light;
  }
  // ///////////////////////////////////////////////////////////////////////////

  xstructure Supercell::calculatePrimitiveStructure() const { //CO20180409
    xstructure pcStructure;
    LightCopy(_inStructure, pcStructure);  //ME20211019 - do not copy symmetry
    //ME20200324 - Setting LatticeReduction_avoid to true can results in
    // primitive cells with slightly different lattice parameters, especially
    // for monoclinic cells. This error can propagate and break mappings from
    // the conventional to the primitive cell.
    pcStructure.LatticeReduction_avoid = false;
    pcStructure.Standard_Primitive_UnitCellForm();
    pcStructure.CleanStructure();
    return pcStructure;
  }

  const xstructure& Supercell::getPrimitiveStructure() const {
    return _pcStructure;
  }

  const xstructure& Supercell::getInputStructure() const {
    return _inStructure;
  }

  //ME20200117
  const xstructure& Supercell::getOriginalStructure() const {
    return _inStructure_original;
  }

  // ///////////////////////////////////////////////////////////////////////////

  const xstructure& Supercell::getInputStructureLight() const {
    return _inStructure_light;
    //xstructure a;
    //stringstream POSCAR;
    //POSCAR.str("");
    //POSCAR << _inStructure;
    //POSCAR >> a;
    //enable inequivalent flag to work
    //for(uint i=0;i<a.atoms.size();i++){
    //  a.atoms[i].equivalent=_inStructure.atoms[i].equivalent;
    //  a.atoms[i].is_inequivalent=_inStructure.atoms[i].is_inequivalent;
    //  a.atoms[i].num_equivalents=_inStructure.atoms[i].num_equivalents;
    //}
    //enable inequivalent flag to work
    //a.write_inequivalent_flag = _inStructure.write_inequivalent_flag;
    //a.info = _inStructure.info;
    //return a;
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20190715 - added const to use function with const Supercell &
  uint Supercell::getNumberOfAtoms() const {
    return _scStructure.atoms.size();
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20190715 - added const to use function with const Supercell &
  uint Supercell::getNumberOfUniqueAtoms() const {
    return _scStructure.iatoms.size();
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20190715 - added const to use function with const Supercell &
  uint Supercell::getNumberOfEquivalentAtomsOfType(int i) const { //CO20190218
#ifndef __OPTIMIZE
    if (i >= (int)_scStructure.iatoms.size()) {
      string message = "Wrong index " + aurostd::utype2string<int>(i) + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }
#endif
    return _scStructure.iatoms[i].size();
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20190715 - added const to use function with const Supercell &
  int Supercell::getUniqueAtomID(int i) const {
#ifndef __OPTIMIZE
    if (i >= (int)_scStructure.iatoms.size()) {
      string message = "Wrong index " + aurostd::utype2string<int>(i) + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }
#endif
    return _scStructure.iatoms[i][0];
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20190715 - added const to use function with const Supercell &
  int Supercell::getUniqueAtomID(int i, int j) const {
#ifndef __OPTIMIZE
    if (i >= (int)_scStructure.iatoms.size()) {
      string message = "Wrong index[1] " + aurostd::utype2string<int>(i) + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }

    if (j >= (int)_scStructure.iatoms[i].size()) {
      string message = "Wrong index[2] " + aurostd::utype2string<int>(i) + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }
#endif
    return _scStructure.iatoms[i][j];
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20190715 - added const to use function with const Supercell &
  const _atom& Supercell::getUniqueAtom(int i) const {
    return _scStructure.atoms[getUniqueAtomID(i)];
  }

  // ///////////////////////////////////////////////////////////////////////////

  int Supercell::atomGoesTo(const _sym_op& symOp, int atomID, int centerID, bool translate) { //CO20190218
    //CO
    //change so that if we can retrieve from fullsymbasis,  we do so
    //this functions looks at symop, and asks by applying it, which atom does atomID become?
    //in fgroup, look at basis_atoms_map, return atom at index atomID
#ifndef __OPTIMIZE
    if (atomID >= (int)_scStructure.atoms.size()) {
      string message = "Wrong atomID index " + aurostd::utype2string<int>(atomID) + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }

    if (centerID >= (int)_scStructure.atoms.size()) {
      string message = "Wrong centerID index " + aurostd::utype2string<int>(centerID) + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }
#endif

    //ME20191219 - use basis_atoms_map
    if (symOp.basis_map_calculated) return symOp.basis_atoms_map[atomID];

    // Get the center atom center...
    if (translate && symOp.is_agroup) center(centerID);

    // Transform atom...
    //DX _atom rotatedAtom = SYM::ApplyAtom(_scStructure.atoms[atomID], symOp,
    //DX                                   _scStructure, TRUE, FALSE, _derivative_structure);  //CO no roff
    _atom rotatedAtom;
    if (!SYM::ApplyAtomValidate(_scStructure.atoms[atomID], rotatedAtom, symOp, _scStructure, _skew, TRUE, FALSE, _sym_eps)) {
      string message = "Illegitimate mapping.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    // Find its id...
    int l = 0;
    for (; l < (int)_scStructure.atoms.size(); l++) {
      if (SYM::FPOSMatch(rotatedAtom.fpos, _scStructure.atoms[l].fpos, _scStructure.lattice, _scStructure.f2c, _skew, _sym_eps)) {  //CO NEW, default to symmetry tolerance
        break;
      }
    }

    if (l == (int)_scStructure.atoms.size()) {
#if ERROR_VERBOSE
      int l = 0;
      for (; l < (int)_scStructure.atoms.size(); l++) {
        printXVector(_scStructure.atoms[atomID].fpos, false);
        cout << " -> ";
        printXVector(rotatedAtom.fpos, false);
        cout << " | ";
        printXVector(_scStructure.atoms[l].fpos, false);
        cout << " | ";
        cout << aurostd::modulus(rotatedAtom.fpos - _scStructure.atoms[l].fpos) << std::endl;
      }
#endif
      string message = "Mapping failed.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    // Move center back to zero atom...
    //if( translate && symOp.is_agroup ) center(0);
    if (translate && symOp.is_agroup) center_original();

#if MAP_VERBOSE
    bool will_translate = symOp.is_agroup && translate;
    cerr << "where: " << atomID << " " << centerID << " " << will_translate << " " << l << std::endl;
    cerr << "atomID : " << _scStructure.atoms[atomID] << std::endl;
    cerr << symOp << std::endl;
#endif

    return l;
  }

  // ///////////////////////////////////////////////////////////////////////////

  int Supercell::atomComesFrom(const _sym_op& symOp, int atomID, int centerID, bool translate) { //CO20190218
    //CO
    //this function does the opposite (to above)
    //in basis_atoms_map, return the index of the atom that is atomID
#ifndef __OPTIMIZE
    if (atomID >= (int)_scStructure.atoms.size()) {
      string message = "Wrong atomID index " + aurostd::utype2string<int>(atomID) + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }

    if (centerID >= (int)_scStructure.atoms.size()) {
      string message = "Wrong centerID index " + aurostd::utype2string<int>(centerID) + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }
#endif

    //ME20191219 - use basis_atoms_map
    if (symOp.basis_map_calculated) {
      int natoms = (int) _scStructure.atoms.size();
      for (int at = 0; at < natoms; at++) {
        if (symOp.basis_atoms_map[at] == atomID) return at;
      }
      // If the code makes it past the for-loop, the mapping failed
      string message = "Mapping failed.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    // Get the center atom center...
    if (translate && symOp.is_agroup) center(centerID);

    // Find it
    int l = 0;
    _atom rotatedAtom;  //CO
    for (; l < (int)_scStructure.atoms.size(); l++) {
      //DX _atom rotatedAtom = SYM::ApplyAtom(_scStructure.atoms[l], symOp, _scStructure, TRUE, FALSE, _derivative_structure);  //CO no roff
      if (!SYM::ApplyAtomValidate(_scStructure.atoms[l], rotatedAtom, symOp, _scStructure, _skew, TRUE, FALSE, _sym_eps)) {
        string message = "Illegitimate mapping.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      if (SYM::FPOSMatch(rotatedAtom.fpos, _scStructure.atoms[atomID].fpos, _scStructure.lattice, _scStructure.f2c, _skew, _sym_eps)) {  //CO NEW, default to symmetry tolerance
        break;
      }
    }

    if (l == (int)_scStructure.atoms.size()) {
#if ERROR_VERBOSE
      int l = 0;
      for (; l < (int)_scStructure.atoms.size(); l++) {
        //DX _atom rotatedAtom = SYM::ApplyAtom(_scStructure.atoms[l], symOp, _scStructure, TRUE, FALSE, _derivative_structure);  //CO no roff
        _atom rotatedAtom = SYM::ApplyAtom(_scStructure.atoms[l], symOp, _scStructure, TRUE, FALSE);  //CO no roff
        printXVector(rotatedAtom.fpos, false);
        cout << " | ";
        printXVector(_scStructure.atoms[atomID].fpos, false);
        cout << " | ";
        cout << aurostd::modulus(rotatedAtom.fpos - _scStructure.atoms[atomID].fpos) << std::endl;
      }
#endif
      string message = "Mapping failed.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    // Move center back to zero atom...
    //if( translate && symOp.is_agroup ) center(0);
    if (translate && symOp.is_agroup) center_original();

#if MAP_VERBOSE
    bool will_translate = symOp.is_agroup && translate;
    cerr << "wherefrom: " << atomID << " " << centerID << " " << will_translate << " " << l << std::endl;
    cerr << "atomID : " << _scStructure.atoms[atomID] << std::endl;
    cerr << symOp << std::endl;
#endif

    return l;
  }

  // ///////////////////////////////////////////////////////////////////////////

  const _sym_op& Supercell::getSymOpWhichMatchAtoms(int whichAtomID, int toAtomID, int GROUP) {
    //go through all fgroups, look at basis_atoms_map at index whichatomID, find toAtomID
#ifndef __OPTIMIZE
    if (whichAtomID >= (int)_scStructure.atoms.size()) {
      string message = "Wrong atom1ID index " + aurostd::utype2string<int>(whichAtomID) + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }

    if (toAtomID >= (int)_scStructure.atoms.size()) {
      string message = "Wrong atom2ID index " + aurostd::utype2string<int>(toAtomID) + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }
#endif

    vector<_sym_op>* symPool = NULL;
    if (GROUP == _PGROUP_)
      symPool = &_scStructure.pgroup;
    else if (GROUP == _FGROUP_)
      symPool = &_scStructure.fgroup;
    else {
      string message = "Unknown group type.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    // Apply all symmetry operations on atom1 and find which one produce atom2
    uint iSymOp = 0;
    _atom newAtom;  //CO
    for (; iSymOp < symPool->size(); iSymOp++) {
      //DX _atom newAtom = SYM::ApplyAtom(_scStructure.atoms[whichAtomID], (*symPool)[iSymOp],
      //DX                               _scStructure, TRUE, FALSE, _derivative_structure);  //CO no roff
      if (!SYM::ApplyAtomValidate(_scStructure.atoms[whichAtomID], newAtom, (*symPool)[iSymOp], _scStructure, _skew, TRUE, FALSE, _sym_eps)) {
        string message = "Illegitimate mapping.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      if (SYM::FPOSMatch(newAtom.fpos, _scStructure.atoms[toAtomID].fpos, _scStructure.lattice, _scStructure.f2c, _skew, _sym_eps)) {  //CO NEW, default to symmetry tolerance
        break;
      }
    }

    if (iSymOp == symPool->size()) {
#if ERROR_VERBOSE
      //CO some helpful output
      cout << "_scStructure.fgroup.size()=" << _scStructure.fgroup.size() << std::endl;
      cout << "symPool->size()=" << symPool->size() << std::endl;
      cout << "ATOM " << whichAtomID << ": " << _scStructure.atoms[whichAtomID].fpos << std::endl;
      cout << "ATOM " << toAtomID << ": " << _scStructure.atoms[toAtomID].fpos << std::endl;

      uint l = 0;
      for (; l < symPool->size(); l++) {
        cout << "i=" << l << std::endl;
        cout << (*symPool)[l] << std::endl;
        printXVector(_scStructure.atoms[whichAtomID].fpos, false);
        cout << " -> ";
        //DX _atom newAtom = SYM::ApplyAtom(_scStructure.atoms[whichAtomID], (*symPool)[l],
        //DX                               _scStructure, TRUE, FALSE, _derivative_structure);  //CO no roff
        _atom newAtom = SYM::ApplyAtom(_scStructure.atoms[whichAtomID], (*symPool)[l],
            _scStructure, TRUE, FALSE);  //CO no roff
        printXVector(newAtom.fpos, false);
        cout << " | ";
        printXVector(_scStructure.atoms[toAtomID].fpos, false);
        cout << " | ";
        cout << aurostd::modulus(newAtom.fpos - _scStructure.atoms[toAtomID].fpos) << std::endl;
      }
#endif
      string message = "Mapping failed.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

#if MAP_VERBOSE
    cerr << "getSymOp: " << whichAtomID << " " << toAtomID << " " << std::endl;
    cerr << "whichAtomID: " << _scStructure.atoms[whichAtomID] << std::endl;
    cerr << "toAtomID: " << _scStructure.atoms[toAtomID] << std::endl;
    cerr << (*symPool)[iSymOp] << std::endl;
#endif

    return (*symPool)[iSymOp];
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20190715 - added const to use function with const Supercell &
  int Supercell::pc2scMap(int i) const {
    return _pc2scMap[i];
  }

  //ME20190715 - added const to use function with const Supercell &
  int Supercell::sc2pcMap(int i) const {
    return _sc2pcMap[i];
  }

  // ///////////////////////////////////////////////////////////////////////////

  //CO START
  void Supercell::center(int i) {
    xvector<double> origin(3), frigin(3);
#if CENTER_PRIM
    //_inStructure.ShiftOriginToAtom(_sc2pcMap[i]);
    //_inStructure.BringInCell();
    origin = _inStructure_original.atoms[_sc2pcMap[i]].cpos;
    frigin = _inStructure_original.atoms[_sc2pcMap[i]].fpos;
    _inStructure.origin = origin;
    for (uint i = 0; i < _inStructure.atoms.size(); i++) {
      _inStructure.atoms[i].cpos = _inStructure_original.atoms[i].cpos - origin;
      _inStructure.atoms[i].fpos = _inStructure_original.atoms[i].fpos - frigin;
      _inStructure.atoms[i] = BringInCell(_inStructure.atoms[i], _inStructure.lattice);
    }
#endif
    //_scStructure.ShiftOriginToAtom(i);
    //_scStructure.BringInCell();
    origin = _scStructure_original.atoms[i].cpos;
    frigin = _scStructure_original.atoms[i].fpos;
    _scStructure.origin = origin;
    for (uint i = 0; i < _scStructure.atoms.size(); i++) {
      _scStructure.atoms[i].cpos = _scStructure_original.atoms[i].cpos - origin;
      _scStructure.atoms[i].fpos = _scStructure_original.atoms[i].fpos - frigin;
      _scStructure.atoms[i] = BringInCell(_scStructure.atoms[i], _scStructure.lattice);
    }
  }

  void Supercell::center_original(void) {
    //just a (fast) undo for center(atom);
    //refer to ShiftOriginToAtom in case more properties need to be updated
#if CENTER_PRIM
    _inStructure.origin = _inStructure_original.origin;
    for (uint i = 0; i < _inStructure.atoms.size(); i++) {
      _inStructure.atoms[i].fpos = _inStructure_original.atoms[i].fpos;
      _inStructure.atoms[i].cpos = _inStructure_original.atoms[i].cpos;
    }
#endif
    _scStructure.origin = _scStructure_original.origin;
    for (uint i = 0; i < _scStructure.atoms.size(); i++) {
      _scStructure.atoms[i].fpos = _scStructure_original.atoms[i].fpos;
      _scStructure.atoms[i].cpos = _scStructure_original.atoms[i].cpos;
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20190715 - added const to use function with const Supercell &
  const vector<_sym_op>& Supercell::getFGROUP(void) const {
    return _scStructure.fgroup;
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20191219
  void Supercell::getFullBasisAGROUP() {
    string message = "Calculating the full basis for the site point groups of the supercell.";
    message += " This may take a few minutes for high-symmetry structures.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
    if (!SYM::CalculateSitePointGroup_EquivalentSites(_scStructure, _sym_eps)) {
      message = "Could not calculate the bases of the site point groups.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
  }

  //ME20191219
  bool Supercell::fullBasisCalculatedAGROUP() {
    uint natoms = _scStructure.atoms.size();
    for (uint at = 0; at < natoms; at++) {
      const vector<_sym_op>& agroup = getAGROUP(at);
      for (uint symop = 0; symop < agroup.size(); symop++) {
        if (!agroup[symop].basis_map_calculated) return false;
      }
    }
    return true;
  }

  //ME20190715 - added const to use function with const Supercell &
  const vector<vector<_sym_op> >& Supercell::getAGROUP(void) const {
    return _scStructure.agroup;
  }
  //ME20190715 - added const to use function with const Supercell &
  const vector<_sym_op>& Supercell::getAGROUP(int i) const {
    return _scStructure.agroup[i];
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20190715 - added const to use function with const Supercell &
  double Supercell::getEPS(void) const {
    return _sym_eps;
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20190715 - added const to use function with const Supercell &
  bool Supercell::isDerivativeStructure(void) const {
    return _derivative_structure;
  }
  //CO END

  // ///////////////////////////////////////////////////////////////////////////

  //ME20190715 - added const to use function with const Supercell &
  string Supercell::getUniqueAtomSymbol(int i) const {
#ifndef __OPTIMIZE
    if (i >= (int)_scStructure.iatoms.size()) {
      string message = "Wrong index " + aurostd::utype2string<int> (i) + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }
#endif
    return _scStructure.atoms[_scStructure.iatoms[i][0]].cleanname;
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20190715 - added const to use function with const Supercell &
  double Supercell::getUniqueAtomMass(int i) const {
#ifndef __OPTIMIZE
    if (i >= (int)_scStructure.iatoms.size()) {
      string message = "Wrong index " + aurostd::utype2string<int> (i) + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }
#endif
    //return( GetAtomMass(_scStructure.atoms[_scStructure.iatoms[i][0]].cleanname) * KILOGRAM2AMU ); //JAHNATEK ORIGINAL
    //CO START
    //double mass = GetAtomMass(_scStructure.atoms[_scStructure.iatoms[i][0]].cleanname); ME20190111 - too slow since version 3.216
    double mass = GetAtomMass(_scStructure.atoms[_scStructure.iatoms[i][0]].atomic_number);  //ME20190111
    if (mass == NNN) {
      string message = "Unknown atom types.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }
    return mass * KILOGRAM2AMU;
    //CO END
  }

  // ///////////////////////////////////////////////////////////////////////////
  //ME20190715 - added const to use function with const Supercell &
  double Supercell::getAtomMass(int i) const {
#ifndef __OPTIMIZE
    if (i >= (int)_scStructure.atoms.size()) {
      string message = "Wrong index " + aurostd::utype2string<int> (i) + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }
#endif
    //return (GetAtomMass(_scStructure.atoms[i].cleanname) * KILOGRAM2AMU); ME20190111 - too slow since version 3.216
    return (GetAtomMass(_scStructure.atoms[i].atomic_number) * KILOGRAM2AMU);  //ME20190111
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20190715 - added const to use function with const Supercell &
  int Supercell::getAtomNumber(int i) const {
#ifndef __OPTIMIZE
    if (i >= (int)_scStructure.atoms.size()) {
      string message = "Wrong index " + aurostd::utype2string<int> (i) + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }
#endif
    return ((int)GetAtomNumber(_scStructure.atoms[i].cleanname));
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20200116 - Calculate the real space vectors for the phases once.
  // Calculating them once for all atoms is very quick and significantly speeds
  // up dynamical matrix calculations.
  void Supercell::calculatePhaseVectors() {
    uint natoms_in = _inStructure.atoms.size();
    uint natoms_sc = _scStructure.atoms.size();
    phase_vectors.clear();
    phase_vectors.resize(natoms_in, vector<vector<xvector<double> > >(natoms_sc));
    xvector<double> rf(3), rc(3), delta(3), pf(3), pc(3);
    double rshell = 0.0;
    int at1pc = -1, at2pc = -1, at1sc = -1, at2sc = -1, centerIDsc = -1;

    for (uint centerID = 0; centerID < natoms_in; centerID++) {
      centerIDsc = pc2scMap(centerID);
      at2pc = sc2pcMap(centerIDsc);
      at2sc = pc2scMap(at2pc);
      for (uint atomID = 0; atomID < natoms_sc; atomID++) {
        at1pc = sc2pcMap(atomID);
        at1sc = pc2scMap(at1pc);
        // Get the nearest image of atomID and determine the shell radius
        rc = SYM::minimizeDistanceCartesianMethod(_scStructure.atoms[atomID].cpos, _scStructure.atoms[centerIDsc].cpos, _scStructure.lattice);
        rf = _scStructure.c2f * rc;
        rshell = aurostd::modulus(rc);
        delta = _scStructure.atoms[at1sc].cpos - _scStructure.atoms[at2sc].cpos;
        // Get the phase vectors for all atoms that sit on the shell
        if (!_isShellRestricted || (rshell <= _maxShellRadius[centerIDsc] + _FLOAT_TOL_)) {
          for (int ii = -1; ii <= 1; ii++) {
            for (int jj = -1; jj <= 1; jj++) {
              for (int kk = -1; kk <= 1; kk++) {
                pf[1] = rf[1] + ii;
                pf[2] = rf[2] + jj;
                pf[3] = rf[3] + kk;
                pc = _scStructure.f2c * pf;
                if (aurostd::isequal(aurostd::modulus(pc), rshell, _FLOAT_TOL_)) {
                  pc -= delta;
                  phase_vectors[centerID][atomID].push_back(pc);
                }
              }
            }
          }
        }
      }
    }
  }

  // ///////////////////////////////////////////////////////////////////////////
  //ME20180827 -- overloaded to calculate derivatives
  bool Supercell::calcShellPhaseFactor(int atomID, int centerID, const xvector<double>& qpoint,
      xcomplex<double>& phase) {
    xvector<xcomplex<double> > placeholder;
    int i;
    return calcShellPhaseFactor(atomID, centerID, qpoint, phase, i, placeholder, false);
  }

  //ME20200116 - rewritten to accommodate new phase vectors
  bool Supercell::calcShellPhaseFactor(int atomID, int centerID, const xvector<double>& qpoint,
      xcomplex<double>& phase, int& multiplicity,
      xvector<xcomplex<double> >& derivative, bool calc_derivative) {
    centerID = sc2pcMap(centerID);
    const vector<xvector<double> >& vec = phase_vectors[centerID][atomID];
    multiplicity = (int) vec.size();
    phase.re = 0.0;
    phase.im = 0.0;
    if (calc_derivative) {
      for (int i = 1; i < 4; i++) {
        derivative[i].re = 0.0;
        derivative[i].im = 0.0;
      }
    }

    xcomplex<double> p;
    for (int i = 0; i < multiplicity; i++) {
      p = exp(iONE * scalar_product(qpoint, vec[i]))/((double) multiplicity);
      phase += p;
      if (calc_derivative) {
        for (int j = 1; j < 4; j++) derivative[j] += vec[i][j] * p * iONE;
      }
    }

    return (multiplicity > 0);
  }

}  // namespace apl

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
