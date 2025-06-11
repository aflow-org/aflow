// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// This file contains the ForceConstantCalculator class, which calculates
// harmonic force constants using the direct method or gamma-point density
// functional perturbation theory.

#include "aflow_apl.h"
#define _DEBUG_APL_HARM_IFCS_ false

using std::vector;
using std::string;

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                         CONSTRUCTORS/DESTRUCTORS                         //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  ForceConstantCalculator::ForceConstantCalculator(ostream& oss) : xStream(oss) {
    free();
    _directory = "./";
  }

  ForceConstantCalculator::ForceConstantCalculator(Supercell& sc, ofstream& mf, ostream& oss) : xStream(mf,oss) {
    free();
    _supercell = &sc;
    _sc_set = true;
    xStream::initialize(mf, oss);
    _directory = _supercell->_directory;
  }

  ForceConstantCalculator::ForceConstantCalculator(Supercell& sc, const aurostd::xoption& opts, ofstream& mf, ostream& oss) : xStream(mf,oss) {
    free();
    _supercell = &sc;
    _sc_set = true;
    xStream::initialize(mf, oss);
    _directory = _supercell->_directory;
    initialize(opts);
  }

  ForceConstantCalculator::ForceConstantCalculator(const ForceConstantCalculator& that) : xStream(*that.getOFStream(),*that.getOSS()) {
    if (this != &that) free();
    copy(that);
  }

  ForceConstantCalculator& ForceConstantCalculator::operator=(const ForceConstantCalculator& that) {
    if (this != &that) free();
    copy(that);
    return *this;
  }

  ForceConstantCalculator::~ForceConstantCalculator() {
    free();
  }

  void ForceConstantCalculator::clear(Supercell& sc) {
    free();
    _supercell = &sc;
    _directory = _supercell->_directory;
  }

  void ForceConstantCalculator::copy(const ForceConstantCalculator& that) {
    if (this == &that) return;
    xStream::copy(that);
    _bornEffectiveChargeTensor = that._bornEffectiveChargeTensor;
    _dielectricTensor = that._dielectricTensor;
    _directory = that._directory;
    _forceConstantMatrices = that._forceConstantMatrices;
    _initialized = that._initialized;
    _isPolarMaterial = that._isPolarMaterial;
    _method = that._method;
    _sc_set = that._sc_set;
    _supercell = that._supercell;
    xInputs = that.xInputs;
    _calculateZeroStateForces = that._calculateZeroStateForces;
    AUTO_GENERATE_PLUS_MINUS = that.AUTO_GENERATE_PLUS_MINUS;
    DISTORTION_MAGNITUDE = that.DISTORTION_MAGNITUDE;
    DISTORTION_INEQUIVONLY = that.DISTORTION_INEQUIVONLY;
    DISTORTION_SYMMETRIZE = that.DISTORTION_SYMMETRIZE;
    GENERATE_ONLY_XYZ = that.GENERATE_ONLY_XYZ;
    USER_GENERATE_PLUS_MINUS = that.USER_GENERATE_PLUS_MINUS;
  }

  void ForceConstantCalculator::free() {
    xInputs.clear();
    _bornEffectiveChargeTensor.clear();
    _dielectricTensor.clear();
    _directory = "";
    _forceConstantMatrices.clear();
    _initialized = false;
    _isPolarMaterial = false;
    _method = "";
    _sc_set = false;
    _supercell = NULL;
    _calculateZeroStateForces = false;
    AUTO_GENERATE_PLUS_MINUS = true;   //CO
    DISTORTION_MAGNITUDE = 0.0;
    DISTORTION_INEQUIVONLY = true;   //CO20190116
    DISTORTION_SYMMETRIZE = true;   //CO20190116
    GENERATE_ONLY_XYZ = false;
    USER_GENERATE_PLUS_MINUS = false;  //CO
  }

  void ForceConstantCalculator::initialize(const aurostd::xoption& opts, ofstream& mf, ostream& oss) {
    xStream::initialize(mf, oss);
    initialize(opts);
  }

  void ForceConstantCalculator::initialize(const aurostd::xoption& opts) {
    string message = "";
    if (!_sc_set) {
      message = "Supercell pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    xInputs.clear();
    _method = opts.getattachedscheme("ENGINE");
    _isPolarMaterial = opts.flag("POLAR");
    if (_method == "DM") {  // Direct method - see README_AFLOW_APL.TXT
      _calculateZeroStateForces = opts.flag("ZEROSTATE");
      string autopm = opts.getattachedscheme("DPM");
      AUTO_GENERATE_PLUS_MINUS = (!autopm.empty() && ((autopm[0] == 'A') || (autopm[0] == 'a')));
      DISTORTION_MAGNITUDE = aurostd::string2utype<double>(opts.getattachedscheme("DMAG"));
      DISTORTION_INEQUIVONLY = opts.flag("DINEQUIV_ONLY");
      DISTORTION_SYMMETRIZE = opts.flag("DSYMMETRIZE");
      GENERATE_ONLY_XYZ = opts.flag("DXYZONLY");
      USER_GENERATE_PLUS_MINUS = opts.flag("DPM");
    } else if (_method != "LR") {  // Linear response - see README_AFLOW_APL.TXT
      message = "Unknown method " + _method + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }
    _initialized = true;
  }

  //xStream initializers
  void ForceConstantCalculator::initialize(ostream& oss) {
    xStream::initialize(oss);
  }

  void ForceConstantCalculator::initialize(ofstream& mf, ostream& oss) {
    xStream::initialize(mf, oss);
  }

}  // namespace apl

namespace apl {

  const vector<vector<xmatrix<double> > >& ForceConstantCalculator::getForceConstants() const {
    return _forceConstantMatrices;
  }

  const vector<xmatrix<double> >& ForceConstantCalculator::getBornEffectiveChargeTensor() const {
    return _bornEffectiveChargeTensor;
  }

  const xmatrix<double>& ForceConstantCalculator::getDielectricTensor() const {
    return _dielectricTensor;
  }

  bool ForceConstantCalculator::isPolarMaterial() const {
    return _isPolarMaterial;
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                             FORCE CONSTANTS                              //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  bool ForceConstantCalculator::runVASPCalculations(_xinput& xInput, _aflags& _aflowFlags,
      _kflags& _kbinFlags, _xflags& _xFlags, string& AflowIn) {
    string message = "";
    if (!_initialized) {
      message = "Not initialized";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    if (!_sc_set) {
      message = "Supercell pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }

    if (!_supercell->isConstructed()) {
      message = "The supercell structure has not been initialized yet.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }

    if (_method.empty()) {
      message = "Calculation method not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }

    xInputs.clear();
    xInput.xvasp.AVASP_arun_mode = "APL";

    bool stagebreak = false;
    if (_method == "DM") { // Direct method - see README_AFLOW_APL.TXT
      stagebreak = runVASPCalculationsDM(xInput, _aflowFlags, _kbinFlags, _xFlags, AflowIn);
    } else if (_method == "LR") {  // Linear response - see README_AFLOW_APL.TXT
      xInputs.push_back(xInput);
      stagebreak = runVASPCalculationsLR(xInputs[0], _aflowFlags, _kbinFlags, _xFlags, AflowIn);
    } else {
      return false;
    }

    if (_isPolarMaterial) {
      xInputs.push_back(xInput);
      stagebreak = (runVASPCalculationsBE(xInputs.back(), _aflowFlags, _kbinFlags, _xFlags, AflowIn, xInputs.size()) || stagebreak);
    }
    return stagebreak;
  }

  // Runs the force constant calculator (main post-processing engine)
  bool ForceConstantCalculator::run() {
    string message = "";
    if (!_initialized) {
      message = "Not initialized";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    if (!_sc_set) {
      message = "Supercell pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    // Check if supercell is already built
    if (!_supercell->isConstructed()) {
      message = "The supercell structure has not been initialized yet.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }

    if (_method.empty()) {
      message = "Calculation method not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }

    if (xInputs.size() == 0) {
      message = "No DFT calculations found.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    // First pass - check if any of the calculations ran (gives no error message)
    if (!outfileFoundAnywherePhonons(xInputs)) return false;

    // Read Born effective charges and dielectric tensor
    if (_isPolarMaterial && !calculateBornChargesDielectricTensor(xInputs.back())) return false;

    if (_method == "DM") {  // Direct method - see README_AFLOW_APL.TXT
      if (!calculateForceConstantsDM()) return false;
    } else if (_method == "LR") {  // Linear response - see README_AFLOW_APL.TXT
      if (!readForceConstantsFromVasprun(xInputs[0])) return false;
    } else {
      return false;
    }

    //ME20191219 - atomGoesTo and atomComesFrom can now use basis_atoms_map.
    // Calculating the full basis ahead of time is much faster than calculating all
    // symmetry operations on-the-fly.
    if (!_supercell->fullBasisCalculatedAGROUP()) _supercell->getFullBasisAGROUP();

    // Symmetrization of the force-constant matrices
    symmetrizeForceConstantMatrices();

    // Force the force-constant matrices to obey the sum-rule conditions
    correctSumRules();

    return true;
  }

  // Symmetrizes the force constant matrices using site point group symmetry
  void ForceConstantCalculator::symmetrizeForceConstantMatrices() {
    bool LDEBUG=(FALSE || _DEBUG_APL_HARM_IFCS_ || XHOST.DEBUG);
    string message = "";
    // Test of stupidity...
    if (!_supercell->getSupercellStructure().agroup_calculated) {
      message = "The site groups have not been calculated yet.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    //CO START
    if (_supercell->getEPS() == AUROSTD_NAN) {
      message = "Need to define symmetry tolerance.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    //CO END

    message = "Symmetrizing the force constant matrices.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);

    vector<xmatrix<double> > row;
    for (uint i = 0; i < _supercell->getNumberOfAtoms(); i++) {
      const vector<_sym_op>& agroup = _supercell->getAGROUP(i);  //CO //CO20190218
      if (agroup.size() == 0) {
        message = "Site point group operations are missing.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
      }

      for (uint j = 0; j < _supercell->getNumberOfAtoms(); j++) {
        if(LDEBUG){ //CO20190218
          cerr << __AFLOW_FUNC__ << " compare original m=" << std::endl;
          cerr << _forceConstantMatrices[i][j] << std::endl;
        }
        xmatrix<double> m(3, 3); //CO20190218
        for (uint symOpID = 0; symOpID < agroup.size(); symOpID++) {
          const _sym_op& symOp = agroup[symOpID];

          try {
            // int l = _supercell.atomComesFrom(symOp, j, i, FALSE);  //CO NEW //CO20190218
            //ME20191219 - atomGoesTo now uses basis_atoms_map; keep translation option in case
            // the basis has not been calculated for some reason
            int l = _supercell->atomGoesTo(symOp, j, i, true); //JAHNATEK ORIGINAL //CO20190218
            m = m + (inverse(symOp.Uc) * _forceConstantMatrices[i][l] * symOp.Uc);  //JAHNATEK ORIGINAL //CO20190218
            //m = m + (symOp.Uc * _forceConstantMatrices[i][l] * inverse(symOp.Uc));  //CO NEW //CO20190218
            if(LDEBUG){ //CO20190218
              std::cerr << __AFLOW_FUNC__ << " atom[" << l << "].cpos=" << _supercell->getSupercellStructure().atoms[l].cpos << std::endl;
              std::cerr << __AFLOW_FUNC__ << " atom[" << j << "].cpos=" << _supercell->getSupercellStructure().atoms[j].cpos << std::endl;
              std::cerr << __AFLOW_FUNC__ << " agroup(" << l << " -> " << j << ")=" << std::endl;
              std::cerr << symOp.Uc << std::endl;
              std::cerr << __AFLOW_FUNC__ << " forceConstantMatrices[i=" << i << "][l=" << l << "]=" << std::endl;
              std::cerr << _forceConstantMatrices[i][l] << std::endl;
              std::cerr << __AFLOW_FUNC__ << " with new m=" << std::endl;
              std::cerr << (symOp.Uc * _forceConstantMatrices[i][l] * inverse(symOp.Uc)) << std::endl;
            }
          } catch (aurostd::xerror& e) {
            message = "Mapping problem " + aurostd::utype2string<int>(j);
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
          }
        }
        m = ( 1.0 / agroup.size() ) * m; //CO20190218
        row.push_back(m);
      }
      _forceConstantMatrices[i] = row;
      row.clear();
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20200504 - this function needs to be rewritten to be more clear
  // Enfoces acoustic sum rules
  void ForceConstantCalculator::correctSumRules() {
    //ME20200504
    // sum appears to be the self-interaction term (diagonal terms) of the
    // force constant matrices. See http://cmt.dur.ac.uk/sjc/thesis_prt/node83.html
    xmatrix<double> sum(3, 3);//, sum2(3, 3); OBSOLETE ME20200504 - not used

    for (uint i = 0; i < _supercell->getNumberOfAtoms(); i++) {
      //ME20200504 - sums are not used or cleared before they are used
      //[OBSOLETE] // Get SUMs
      //[OBSOLETE] for (int j = 0; j < _supercell->getNumberOfAtoms(); j++) {
      //[OBSOLETE]   if (i != j) {
      //[OBSOLETE]     sum = sum + _forceConstantMatrices[i][j];
      //[OBSOLETE]     //sum2 = sum2 + trasp(_forceConstantMatrices[j][i]);
      //[OBSOLETE]   }
      //[OBSOLETE] }

      // Correct SUM2
      //ME20200504 - This appears to enforce the invariance of the force constants
      // upon permutations
      for (uint j = 0; j < _supercell->getNumberOfAtoms(); j++) {
        if (i == j) continue;
        _forceConstantMatrices[i][j] = 0.5 * (_forceConstantMatrices[i][j] + trasp(_forceConstantMatrices[j][i]));
        _forceConstantMatrices[j][i] = trasp(_forceConstantMatrices[i][j]);
      }

      // Get SUMs again
      sum.clear();
      //sum2.clear(); OBSOLETE ME20200504 - not used
      for (uint j = 0; j < _supercell->getNumberOfAtoms(); j++) {
        if (i != j) {
          sum = sum + _forceConstantMatrices[i][j];
          //sum2 = sum2 + trasp(_forceConstantMatrices[j][i]);  // OBSOLETE ME20200504 - not used
        }
      }

      // Correct SUM1 to satisfied
      //ME20200504 - Self-interaction term
      _forceConstantMatrices[i][i] = -sum;
    }
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                          BORN/DIELECTRIC TENSOR                          //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  // Sets up the calculations that determine the Born effective charges and
  // the dielectric tensor
  bool ForceConstantCalculator::runVASPCalculationsBE(_xinput& xInput, _aflags& _aflowFlags,
      _kflags& _kbinFlags, _xflags& _xFlags, string& _AflowIn, uint ncalcs) {
    bool stagebreak = false;

    xInput.setXStr(_supercell->getInputStructure());
    xInput.getXStr().title = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xInput.getXStr().title);
    if(xInput.getXStr().title.empty()){xInput.getXStr().buildGenericTitle(true,false);}
    xInput.getXStr().title+=" Born effective charges/dielectric tensor";

    // For VASP, use the standardized aflow.in creator
    if(xInput.AFLOW_MODE_VASP) {
      xInput.xvasp.AVASP_arun_runname = aurostd::utype2string<uint>(ncalcs) + "_" + _AFLOW_APL_BORN_EPSILON_RUNNAME_;
      xInput.xvasp.aopts.flag("APL_FLAG::AVASP_LR", false);
      xInput.xvasp.aopts.flag("APL_FLAG::AVASP_BORN", true);
      // Set POSCAR to VASP5 format
      xInput.getXStr().is_vasp4_poscar_format = false;
      xInput.getXStr().is_vasp5_poscar_format = true;
      stagebreak = (createAflowInPhonons(_aflowFlags, _kbinFlags, _xFlags, xInput) || stagebreak);
    } else if (xInput.AFLOW_MODE_AIMS) {
      string runname = _AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_;
      xInput.setDirectory(_directory + "/" + runname );
      if (!filesExistPhonons(xInput)) {
        string message = "Creating " + xInput.getDirectory();
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
        createAflowInPhononsAIMS(_aflowFlags, _kbinFlags, _xFlags, _AflowIn, xInput, *p_FileMESSAGE);
        stagebreak = true;
      }
    }
    return stagebreak;
  }

  //////////////////////////////////////////////////////////////////////////////

  // Calculates the dielectric tensor and Born effective charges
  bool ForceConstantCalculator::calculateBornChargesDielectricTensor(const _xinput& xinpBE) {
    stringstream message;
    // Parse effective charges from OUTCAR
    if (xinpBE.AFLOW_MODE_VASP) {
      string directory = xinpBE.xvasp.Directory;
      string infilename = directory + "/OUTCAR.static";

      if (!aurostd::EFileExist(infilename, infilename)) {
        infilename = directory + string("/OUTCAR");
        if (!aurostd::EFileExist(infilename, infilename)) {
          message << "The OUTCAR file in " << directory << " is missing.";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
          return false;
        }
      }
    } else if (xinpBE.AFLOW_MODE_AIMS) {
      message << "This functionality has not been implemented yet.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    } else {
      return false;
    }

    if (xinpBE.AFLOW_MODE_VASP) readBornEffectiveChargesFromOUTCAR(xinpBE);
    else if (xinpBE.AFLOW_MODE_AIMS) readBornEffectiveChargesFromAIMSOUT();

    // Enforce ASR (Acoustic sum rules)
    symmetrizeBornEffectiveChargeTensors();

    // Parse epsilon from OUTCAR
    if(xinpBE.AFLOW_MODE_VASP) readDielectricTensorFromOUTCAR(xinpBE);
    if(xinpBE.AFLOW_MODE_AIMS) readDielectricTensorFromAIMSOUT();

    message << "Dielectric tensor: ";
    for (int a = 1; a <= 3; a++)
      for (int b = 1; b <= 3; b++)
        message << std::fixed << std::setw(5) << std::setprecision(3) << _dielectricTensor(a, b) << " ";

    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
    return true;
  }

  //////////////////////////////////////////////////////////////////////////////
  void ForceConstantCalculator::readBornEffectiveChargesFromAIMSOUT(void) {
    string message = "This functionality has not been implemented yet.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
  }

  void ForceConstantCalculator::readBornEffectiveChargesFromOUTCAR(const _xinput& xinp) {  //ME20190113
    string directory = xinp.xvasp.Directory;  //ME20190113
    string message = "";

    //CO START
    string infilename = directory + string("/OUTCAR.static");

    if (!aurostd::EFileExist(infilename, infilename)) {
      // We already know that one of the files exists, so
      // if OUTCAR.static was not found, it must be OUTCAR
      infilename = directory + string("/OUTCAR");
    }

    // Open our file
    //CO START
    vector<string> vlines;
    aurostd::efile2vectorstring(infilename, vlines);
    if (!vlines.size()) {
      //CO END
      message = "Cannot open input file OUTCAR.static.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }

    string line = "";
    uint line_count = 0;  //CO
    string KEY = ""; //ME20181226
    if (DEFAULT_APL_USE_LEPSILON) { //ME20181226
      KEY = string("BORN EFFECTIVE CHARGES (in e, cummulative output)");//ME20181226
    } else { //ME20181226
      KEY = string("BORN EFFECTIVE CHARGES (including local field effects)"); //ME20181226
    } //ME20181226

    while (true) {
      // Get line
      //CO START
      if (line_count == vlines.size()) {
        message = "No information on Born effective charges in OUTCAR file.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
      }
      line = vlines[line_count++];  //CO

      // Check for our key line
      if (line.size() < KEY.size()) continue;
      if (line.find(KEY) != string::npos) break;
    }
    // Read in all ...
    xmatrix<double> m(3, 3);
    vector<string> tokens;
    //CO START
    line = vlines[line_count++]; // Skip line "----------------...."
    for (uint i = 0; i < _supercell->getInputStructure().atoms.size(); i++) {
      // Get atom ID but not use it...
      line = vlines[line_count++];

      // Get its charge tensor
      for (int j = 1; j <= 3; j++) {
        line = vlines[line_count++];
        aurostd::string2tokens(line, tokens, string(" "));
        m(j, 1) = aurostd::string2utype<double>(tokens.at(1));
        m(j, 2) = aurostd::string2utype<double>(tokens.at(2));
        m(j, 3) = aurostd::string2utype<double>(tokens.at(3));
        tokens.clear();
      }

      // Store it
      _bornEffectiveChargeTensor.push_back(m);
    }
    //CO END
  }

  //////////////////////////////////////////////////////////////////////////////
  void ForceConstantCalculator::symmetrizeBornEffectiveChargeTensors(void) {
    //CO START
    // Test of stupidity...
    stringstream message;
    if (_supercell->getEPS() == AUROSTD_NAN) {
      message << "Symmetry tolerance not defined.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    //CO END
    // Show charges
    message << "Input born effective charge tensors (for primitive cell):";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
    for (uint i = 0; i < _bornEffectiveChargeTensor.size(); i++) {
      int id = i;
      message << "Atom [" << aurostd::PaddedNumString(id, 3) << "] ("
        << std::setw(2) << _supercell->getInputStructure().atoms[id].cleanname
        << ") Born effective charge = ";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
      for (int a = 1; a <= 3; a++)
        for (int b = 1; b <= 3; b++)
          message << std::fixed << std::setw(5) << std::setprecision(3) << _bornEffectiveChargeTensor[i](a, b) << " ";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
    }

    // Step1
    for (uint i = 0; i < _supercell->getNumberOfUniqueAtoms(); i++) {
      int basedUniqueAtomID = _supercell->getUniqueAtomID(i);

      xmatrix<double> sum(3, 3);
      for (uint j = 0; j < _supercell->getNumberOfEquivalentAtomsOfType(i); j++) { //CO20190218
        try {  //CO
          const _sym_op& symOp = _supercell->getSymOpWhichMatchAtoms(_supercell->getUniqueAtomID(i, j), basedUniqueAtomID, _FGROUP_);
          sum += inverse(symOp.Uc) * _bornEffectiveChargeTensor[_supercell->sc2pcMap(_supercell->getUniqueAtomID(i, j))] * symOp.Uc;
        }
        //CO START
        catch (aurostd::xerror& e) {
          stringstream message;
          message << "Mapping problem " << _supercell->getUniqueAtomID(i, j) << " <-> " << basedUniqueAtomID << "?";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
        }
        //CO END
      }

      sum = (1.0 / _supercell->getNumberOfEquivalentAtomsOfType(i)) * sum; //CO20190218

      for (uint j = 0; j < _supercell->getNumberOfEquivalentAtomsOfType(i); j++) //CO20190218
        _bornEffectiveChargeTensor[_supercell->sc2pcMap(_supercell->getUniqueAtomID(i, j))] = sum;
    }

    // Step2
    vector<xmatrix<double> > newbe = _bornEffectiveChargeTensor;
    const vector<vector<_sym_op> >& agroup = _supercell->getAGROUP();  //CO
    for (uint i = 0; i < _supercell->getNumberOfAtoms(); i++) {
      // Translate the center to this atom
      _supercell->center(i);

      xmatrix<double> sum(3, 3);
      for (uint symOpID = 0; symOpID < agroup[i].size(); symOpID++) {
        const _sym_op& symOp = agroup[i][symOpID];
        sum = sum + (inverse(symOp.Uc) * _bornEffectiveChargeTensor[_supercell->sc2pcMap(i)] * symOp.Uc);
      }
      newbe[_supercell->sc2pcMap(i)] = (1.0 / agroup[i].size()) * sum;
    }
    // Translate the center back
    _supercell->center_original();  //CO

    _bornEffectiveChargeTensor.clear();
    _bornEffectiveChargeTensor = newbe;
    newbe.clear();

    // Step 3
    message << "Forcing the acoustic sum rule (ASR). Resulting born effective charges (for the supercell):";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);

    xmatrix<double> sum(3, 3);
    for (uint i = 0; i < _bornEffectiveChargeTensor.size(); i++)
      sum += _bornEffectiveChargeTensor[i];
    sum = (1.0 / _bornEffectiveChargeTensor.size()) * sum;
    for (uint i = 0; i < _bornEffectiveChargeTensor.size(); i++)
      _bornEffectiveChargeTensor[i] -= sum;

    // Make list only for unique atoms
    for (uint i = 0; i < _supercell->getNumberOfUniqueAtoms(); i++)
      newbe.push_back(_bornEffectiveChargeTensor[_supercell->sc2pcMap(_supercell->getUniqueAtomID(i))]);
    _bornEffectiveChargeTensor.clear();
    _bornEffectiveChargeTensor = newbe;
    newbe.clear();

    // Show charges
    for (uint i = 0; i < _supercell->getNumberOfUniqueAtoms(); i++) {
      int id = _supercell->getUniqueAtomID(i);
      message << "Atom [" << aurostd::PaddedNumString(id, 3) << "] ("
        << std::setw(2) << _supercell->getSupercellStructure().atoms[id].cleanname
        << ") Born effective charge = ";
      for (int a = 1; a <= 3; a++)
        for (int b = 1; b <= 3; b++)
          message << std::fixed << std::setw(5) << std::setprecision(3) << _bornEffectiveChargeTensor[i](a, b) << " ";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  void ForceConstantCalculator::readDielectricTensorFromAIMSOUT(void) {
    string message = "This functionality has not been implemented yet.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
  }

  //////////////////////////////////////////////////////////////////////////////
  void ForceConstantCalculator::readDielectricTensorFromOUTCAR(const _xinput& xinp) {  //ME20190113
    string message = "";
    string directory = xinp.xvasp.Directory;  //ME20190113

    //CO START
    string infilename = directory + string("/OUTCAR.static");
    if (!aurostd::EFileExist(infilename, infilename)) {
      // We already checked outside if one of the files exists, so if
      // it is not OUTCAR.static, it must be OUTCAR
      infilename = directory + string("/OUTCAR");
    }

    // Open our file
    vector<string> vlines;
    uint line_count = 0;
    string line;
    aurostd::efile2vectorstring(infilename, vlines);
    if (!vlines.size()) {
      message = "Cannot open input file OUTCAR.";
      aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
    //CO END

    // Find
    string KEY = ""; //ME20181226
    if (DEFAULT_APL_USE_LEPSILON) { //ME20181226
      KEY = string("MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)"); //ME20181226
    } else { //ME20181226
      KEY = string("MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects)"); //ME20181226
    }//ME20181226

    while (true) {
      // Get line
      //CO START
      if (line_count == vlines.size()) {
        message = "No information on dielectric tensor in OUTCAR.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
      }
      line = vlines[line_count++];
      //CO END

      // Check for our key line
      if (line.size() < KEY.size()) continue;
      if (line.find(KEY) != string::npos) break;
    }

    // Read in all ...
    //CO START
    line = vlines[line_count++]; // Skip line "----------------...."

    // Get it
    vector<string> tokens;
    for (int j = 1; j <= 3; j++) {
      line = vlines[line_count++];
      aurostd::string2tokens(line, tokens, string(" "));
      _dielectricTensor(j, 1) = aurostd::string2utype<double>(tokens.at(0));
      _dielectricTensor(j, 2) = aurostd::string2utype<double>(tokens.at(1));
      _dielectricTensor(j, 3) = aurostd::string2utype<double>(tokens.at(2));
      tokens.clear();
    }
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                               FILE OUTPUT                                //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  // Writes the results into xml files
  void ForceConstantCalculator::hibernate() {
    string message = "";
    if (!_initialized) {
      message = "Not initialized";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    if (!_sc_set) {
      message = "Supercell pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    string base = _directory + "/" + DEFAULT_APL_FILE_PREFIX;
    string filename = aurostd::CleanFileName(base + DEFAULT_APL_HARMIFC_FILE);
    message = "Writing harmonic IFCs into file " + filename + ".";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
    writeHarmonicIFCs(filename);
    if (_isPolarMaterial) {
      filename = aurostd::CleanFileName(base + DEFAULT_APL_POLAR_FILE);
      message = "Writing harmonic IFCs into file " + filename + ".";
      writeBornChargesDielectricTensor(filename);
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
    }
  }

  // Writes the harmonic force constants into an xml file
  void ForceConstantCalculator::writeHarmonicIFCs(const string& filename) {
    stringstream outfile;
    string tab = " ";

    // Header
    outfile << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << std::endl;
    outfile << "<apl>" << std::endl;
    outfile << tab << "<generator>" << std::endl;
    outfile << tab << tab << "<i name=\"aflow_version\" type=\"string\">" << AFLOW_VERSION << "</i>" << std::endl;
    string time = aflow_get_time_string();
    if (time[time.size() - 1] == '\n') time.erase(time.size() - 1);
    outfile << tab << tab << "<i name=\"date\" type=\"string\">" << time << "</i>" << std::endl;
    // OBSOLETE ME20200427 - we do not compare checksums anymore
    //outfile << tab << tab << "<i name=\"checksum\" file=\"" << _AFLOWIN_ << "\" type=\"" << APL_CHECKSUM_ALGO << "\">"
    //  << std::hex << aurostd::getFileCheckSum(_directory + "/" + _AFLOWIN_ + "", APL_CHECKSUM_ALGO) << "</i>" << std::endl;  //ME20190219
    //outfile.unsetf(std::ios::hex); //CO20190116 - undo hex immediately
    outfile << tab << "</generator>" << std::endl;

    // Force constants
    outfile << tab << "<fcms units=\"eV/Angstrom^2\" cs=\"cartesian\" rows=\""
      << _forceConstantMatrices.size() << "\" cols=\""
      << _forceConstantMatrices[0].size() << "\">" << std::endl;

    outfile << tab << tab << "<varray>" << std::endl;
    for (uint i = 0; i < _forceConstantMatrices.size(); i++) {
      outfile << tab << tab << tab << "<varray row=\"" << i << "\">" << std::endl;
      for (uint j = 0; j < _forceConstantMatrices[i].size(); j++) {
        outfile << tab << tab << tab << tab << "<matrix row=\"" << i
          << "\" col=\"" << j << "\">" << std::endl;
        for (int k = 1; k <= 3; k++) {
          outfile << tab << tab << tab << tab << tab << "<v>";
          for (int l = 1; l <= 3; l++) {
            outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
            outfile << setprecision(15);
            outfile << setw(24) << std::scientific << _forceConstantMatrices[i][j](k, l) << " ";
          }
          outfile << "</v>" << std::endl;
        }
        outfile << tab << tab << tab << tab << "</matrix>" << std::endl;
      }
      outfile << tab << tab << tab << "</varray>" << std::endl;
    }
    outfile << tab << tab << "</varray>" << std::endl;
    outfile << tab << "</fcms>" << std::endl;
    outfile << "</apl>" << std::endl;

    aurostd::stringstream2file(outfile, filename);
    if (!aurostd::FileExist(filename)) {
      string message = "Cannot open output file " + filename + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
  }

  // Writes the Born effective charges and the dielectric tensor into an xml file
  void ForceConstantCalculator::writeBornChargesDielectricTensor(const string& filename) {
    stringstream outfile;
    string tab = " ";

    // Header
    outfile << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << std::endl;
    outfile << "<apl>" << std::endl;
    outfile << tab << "<generator>" << std::endl;
    string time = aflow_get_time_string();
    if (time[time.size() - 1] == '\n') time.erase(time.size() - 1);
    outfile << tab << tab << "<i name=\"date\" type=\"string\">" << time << "</i>" << std::endl;
    // OBSOLETE ME20200428 - Checksums are not used anymore
    //outfile << tab << tab << "<i name=\"checksum\" file=\"" << _AFLOWIN_ << "\" type=\"" << APL_CHECKSUM_ALGO << "\">"
    //  << std::hex << aurostd::getFileCheckSum(_directory + "/" + _AFLOWIN_ + "", APL_CHECKSUM_ALGO) << "</i>" << std::endl;  //ME20190219
    //outfile.unsetf(std::ios::hex); //CO20190116 - undo hex immediately
    outfile << tab << "</generator>" << std::endl;

    // Born effective charge tensors
    outfile << tab << "<born units=\"a.u.\" cs=\"cartesian\">" << std::endl;
    outfile << tab << tab << "<varray>" << std::endl;
    for (uint i = 0; i < _bornEffectiveChargeTensor.size(); i++) {
      int id = _supercell->getUniqueAtomID(i);
      outfile << tab << tab << tab << "<matrix type=\"" << _supercell->getSupercellStructure().atoms[id].cleanname << "\">" << std::endl;
      for (int k = 1; k <= 3; k++) {
        outfile << tab << tab << tab << tab << "<v>";
        for (int l = 1; l <= 3; l++) {
          outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
          outfile << setprecision(8);
          //ME20181030 - fixed prevents hexadecimal output
          outfile << setw(15) << std::fixed << _bornEffectiveChargeTensor[i](k, l) << " ";
        }
        outfile << "</v>" << std::endl;
      }
      outfile << tab << tab << tab << "</matrix>" << std::endl;
    }
    outfile << tab << tab << "</varray>" << std::endl;
    outfile << tab << "</born>" << std::endl;

    // Dielectric tensor
    outfile << tab << "<epsilon units=\"a.u.\" cs=\"cartesian\">" << std::endl;
    outfile << tab << tab << "<matrix>" << std::endl;
    for (int k = 1; k <= 3; k++) {
      outfile << tab << tab << tab << "<v>";
      for (int l = 1; l <= 3; l++) {
        outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
        outfile << setprecision(8);
        //ME20181030 - fixed prevents hexadecimal output
        outfile << setw(15) << std::fixed << _dielectricTensor(k, l) << " ";
      }
      outfile << "</v>" << std::endl;
    }
    outfile << tab << tab << "</matrix>" << std::endl;
    outfile << tab << "</epsilon>" << std::endl;
    outfile << "</apl>" << std::endl;

    aurostd::stringstream2file(outfile, filename);
    if (!aurostd::FileExist(filename)) {
      string message = "Cannot open output file " + filename + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
  }

  void ForceConstantCalculator::saveState(const string& filename) {
    string message = "";
    if (!_sc_set) {
      message = "Supercell pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    if (xInputs.size() == 0) return;  // Nothing to write
    message = "Saving state of the force constant calculator into " + aurostd::CleanFileName(filename) + ".";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
    stringstream out;
    string tag = "[APL_FC_CALCULATOR]";
    out << AFLOWIN_SEPARATION_LINE << std::endl;
    out << tag << "ENGINE=" << _method << std::endl;
    if (xInputs[0].AFLOW_MODE_VASP) out << tag << "AFLOW_MODE=VASP" << std::endl;
    else if (xInputs[0].AFLOW_MODE_AIMS) out << tag << "AFLOW_MODE=AIMS" << std::endl;
    out << AFLOWIN_SEPARATION_LINE << std::endl;
    out << tag << "SUPERCELL=" << _supercell->scell_dim << std::endl;
    out << tag << "INPUT_STRUCTURE=START" << std::endl;
    out << _supercell->getInputStructure();  // No endl necessary
    out << tag << "INPUT_STRUCTURE=STOP" << std::endl;
    out << AFLOWIN_SEPARATION_LINE << std::endl;

    // Distortion parameters for the direct method
    if (_method == "DM") {  // Direct method - see README_AFLOW_APL.TXT
      out << tag << "DISTORTION_MAGNITUDE=" << DISTORTION_MAGNITUDE << std::endl;
      out << tag << "DISTORTION_INEQUIVONLY=" << DISTORTION_INEQUIVONLY << std::endl;
      out << tag << "DISTORTIONS=START" << std::endl;
      int idxRun = 0;
      for (uint i = 0; i < _uniqueDistortions.size(); i++) {
        for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
          out << i << " " << _uniqueDistortions[i][j] << " " << xInputs[idxRun++].xvasp.AVASP_arun_runname;
          if (vvgenerate_plus_minus[i][j]) out << " " << xInputs[idxRun++].xvasp.AVASP_arun_runname;
          out << std::endl;
        }
      }
      out << tag << "DISTORTIONS=STOP" << std::endl;
      out << AFLOWIN_SEPARATION_LINE << std::endl;
      out << tag << "ZEROSTATE=" << _calculateZeroStateForces << std::endl;
      if (_calculateZeroStateForces) out << tag << "ZEROSTATE_RUNNAME=" << xInputs[idxRun++].xvasp.AVASP_arun_runname << std::endl;
      out << AFLOWIN_SEPARATION_LINE << std::endl;
      out << tag << "POLAR=" << _isPolarMaterial << std::endl;
      if (_isPolarMaterial) out << tag << "POLAR_RUNNAME=" << xInputs[idxRun].xvasp.AVASP_arun_runname << std::endl;
      out << AFLOWIN_SEPARATION_LINE << std::endl;
    } else if (_method == "LR") {  // Linear response - see README_AFLOW_APL.TXT
      out << tag << "POLAR=" << _isPolarMaterial << std::endl;
      out << AFLOWIN_SEPARATION_LINE << std::endl;
    }

    aurostd::stringstream2file(out, filename);
    if (!aurostd::FileExist(filename)) {
      message = "Could not save state into file " + filename + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }

  }

  // ME20200501
  // The state reader is designed to read the directory structure and supercell
  // structures from a prior run for post-processing. It cannot create APL
  // aflow.in files and should only be used to read forces for force constant
  // calculations. It is still in development and has only been tested with VASP.
  void ForceConstantCalculator::readFromStateFile(const string& filename) {
    string message = "";
    if (!_sc_set) {
      message = "Supercell pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    message = "Reading state of the phonon calculator from " + filename + ".";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
    if (!aurostd::EFileExist(filename)) {
      message = "Could not find file " + filename + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }

    // Find the ENGINE tag and whether the calculations are VASP or AIMS calculations
    vector<string> vlines, tokens;
    aurostd::efile2vectorstring(filename, vlines);
    uint nlines = vlines.size();
    uint iline = 0;
    _xinput xInput;
    while (++iline < nlines) {
      if (aurostd::substring2bool(vlines[iline], "AFLOW_MODE")) {;
        aurostd::string2tokens(vlines[iline], tokens, "=");
        if (tokens.size() != 2) {
          message = "Tag for AFLOW_MODE is broken.";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
        }
        if (tokens[1] == "VASP") {
          xInput.AFLOW_MODE_VASP = true;
          break;
        } else if (tokens[1] == "AIMS") {
          xInput.AFLOW_MODE_AIMS = true;
          break;
        } else {
          message = "Unknown AFLOW_MODE " + tokens[1] + ".";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
        }
      } else if (aurostd::substring2bool(vlines[iline], "ENGINE")) {
        aurostd::string2tokens(vlines[iline], tokens, "=");
        if (tokens.size() != 2) {
          message = "Tag for ENGINE is broken.";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
        }
        if ((tokens[1] != "LR") && (tokens[1] != "DM")) {
          message = "Unknown value for ENGINE " + tokens[1] + ".";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
        } else {
          _method = tokens[1];
        }
      }
    }
    if (iline >= nlines) {
      if (_method.empty()) message = "ENGINE tag is missing.";
      else message = "AFLOW_MODE tag is missing.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    }

    // Defaults
    xInput.xvasp.AVASP_arun_mode = "APL";
    _isPolarMaterial = DEFAULT_APL_POLAR;
    if (_method == "DM") {
      _calculateZeroStateForces = DEFAULT_APL_ZEROSTATE;
      DISTORTION_MAGNITUDE = DEFAULT_APL_DMAG;
      DISTORTION_INEQUIVONLY = DEFAULT_APL_DINEQUIV_ONLY;
      _uniqueDistortions.clear();
      vvgenerate_plus_minus.clear();
    }

    // Read
    xInputs.clear();
    iline = 0;
    if (_method == "DM") {
      while (++iline < nlines) {
        if (aurostd::substring2bool(vlines[iline], "DISTORTION_MAGNITUDE=")) {
          tokens.clear();
          aurostd::string2tokens(vlines[iline], tokens, "=");
          if (tokens.size() != 2) {
            message = "Tag for DISTORTION_MAGNITUDE is broken.";
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
          }
          DISTORTION_MAGNITUDE = aurostd::string2utype<double>(tokens[1]);
        } else if (aurostd::substring2bool(vlines[iline], "DISTORTION_INEQUIVONLY=")) {
          tokens.clear();
          aurostd::string2tokens(vlines[iline], tokens, "=");
          if (tokens.size() != 2) {
            message = "Tag for DISTORTION_INEQUIVONLY correction is broken.";
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
          }
          DISTORTION_INEQUIVONLY = aurostd::string2utype<bool>(tokens[1]);
        } else if (aurostd::substring2bool(vlines[iline], "DISTORTIONS=START")) {
          xvector<double> distortion(3);
          uint idist = 0;
          while ((iline++ < nlines) && !aurostd::substring2bool(vlines[iline], "DISTORTIONS=STOP")) {
            tokens.clear();
            aurostd::string2tokens(vlines[iline], tokens);
            if ((tokens.size() < 5) || (tokens.size() > 7)) {
              message = "Broken line in DISTORTIONS.";
              throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
            }
            // Distortions
            idist = aurostd::string2utype<uint>(tokens[0]);
            if (idist + 1 > _uniqueDistortions.size()) {
              _uniqueDistortions.push_back(vector<xvector<double> >(0));
              vvgenerate_plus_minus.push_back(vector<bool>(0));
            }
            for (int i = 1; i < 4; i++) distortion[i] = aurostd::string2utype<double>(tokens[i]);
            _uniqueDistortions[idist].push_back(distortion);
            xInputs.push_back(xInput);
            xInputs.back().xvasp.AVASP_arun_runname = tokens[4];
            if (tokens.size() == 5) {
              vvgenerate_plus_minus[idist].push_back(false);
            } else {
              vvgenerate_plus_minus[idist].push_back(true);
              xInputs.push_back(xInput);
              xInputs.back().xvasp.AVASP_arun_runname = tokens[5];
            }
          }
        } else if (aurostd::substring2bool(vlines[iline], "ZEROSTATE=")) {
          tokens.clear();
          aurostd::string2tokens(vlines[iline], tokens, "=");
          if (tokens.size() != 2) {
            message = "Tag for ZEROSTATE calculation is broken.";
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
          }
          _calculateZeroStateForces = aurostd::string2utype<bool>(tokens[1]);
          if (_calculateZeroStateForces) {
            iline++;
            if (iline == nlines) {
              message = "Runname for ZEROSTATE calculation is missing.";
              throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
            }
            tokens.clear();
            aurostd::string2tokens(vlines[iline], tokens, "=");
            if (tokens.size() != 2) {
              message = "Runname tag for ZEROSTATE calculation is broken.";
              throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
            }
            xInputs.push_back(xInput);
            xInputs.back().setXStr(_supercell->getSupercellStructureLight());
            xInputs.back().xvasp.AVASP_arun_runname = tokens[1];
          }
        } else if (aurostd::substring2bool(vlines[iline], "POLAR=")) {
          tokens.clear();
          aurostd::string2tokens(vlines[iline], tokens, "=");
          if (tokens.size() != 2) {
            message = "Tag for POLAR is broken.";
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
          }
          _isPolarMaterial = aurostd::string2utype<bool>(tokens[1]);
          if (_isPolarMaterial) {
            iline++;
            if (iline == nlines) {
              message = "Runname for polar correction is missing.";
              throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
            }
            tokens.clear();
            aurostd::string2tokens(vlines[iline], tokens, "=");
            if (tokens.size() != 2) {
              message = "Runname tag for polar correction is broken.";
              throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
            }
            xInputs.push_back(xInput);
            xInputs.back().setXStr(_supercell->getInputStructureLight());
            xInputs.back().xvasp.AVASP_arun_runname = tokens[1];
          }
        }
      }

      // Done reading - apply distortions to structures
      int idxRun = 0;
      for (uint i = 0; i < _uniqueDistortions.size(); i++) {
        int idAtom = (DISTORTION_INEQUIVONLY ? _supercell->getUniqueAtomID(i) : i );
        for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
          for (uint k = 0; k < (vvgenerate_plus_minus[i][j] ? 2 : 1); k++) {
            xInputs[idxRun].setXStr(_supercell->getSupercellStructureLight());
            xstructure& xstr = xInputs[idxRun].getXStr();
            xstr.atoms[idAtom].cpos += ((k == 0) ? 1.0 : -1.0 ) * DISTORTION_MAGNITUDE * _uniqueDistortions[i][j];
            xstr.atoms[idAtom].fpos = xstr.c2f * xstr.atoms[idAtom].cpos;
            idxRun++;
          }
        }
      }
    } else if (_method == "LR") {
      // Set xInput for the linear response calculation
      xInputs.push_back(xInput);
      xInputs[0].setXStr(_supercell->getSupercellStructureLight());
      xInputs[0].xvasp.AVASP_arun_runname = "1_" + _AFLOW_APL_DFPT_RUNNAME_;

      while (++iline < nlines) {
        if (aurostd::substring2bool(vlines[iline], "POLAR=")) {
          aurostd::string2tokens(vlines[iline], tokens, "=");
          if (tokens.size() != 2) {
            string message = "Tag for POLAR is broken.";
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
          }
          _isPolarMaterial = aurostd::string2utype<bool>(tokens[1]);
          if (_isPolarMaterial) {
            xInputs.push_back(xInput);
            xInputs[1].setXStr(_supercell->getInputStructureLight());
            xInputs[1].xvasp.AVASP_arun_runname = "2_" + _AFLOW_APL_BORN_EPSILON_RUNNAME_;
          }
        }
      }
    }

    // Set directories
    string dir = "";
    for (uint i = 0; i < xInputs.size(); i++) {
      const _xvasp& xvasp = xInputs[i].xvasp;
      dir = _directory + "/ARUN." + xvasp.AVASP_arun_mode + "_" + xvasp.AVASP_arun_runname + "/";
      xInputs[i].setDirectory(aurostd::CleanFileName(dir));
    }
  }

  // ///////////////////////////////////////////////////////////////////////////
  // OBSOLETE ME20200504  - not used
  //[OBSOLETE] void ForceConstantCalculator::printForceConstantMatrices(ostream& os) {
  //[OBSOLETE]   int units = 1;
  //[OBSOLETE]   double conversionFactor = 1.0;

  //[OBSOLETE]   switch (units) {
  //[OBSOLETE]     case (1):
  //[OBSOLETE]       os << "FORCE CONSTANT MATRICES in eV/A^2:" << std::endl;
  //[OBSOLETE]       conversionFactor = 1.0;
  //[OBSOLETE]       break;
  //[OBSOLETE]     case (2):
  //[OBSOLETE]       os << "FORCE CONSTANT MATRICES in 10 Dyn/cm:" << std::endl;
  //[OBSOLETE]       conversionFactor = 1602.17733;
  //[OBSOLETE]       break;
  //[OBSOLETE]     case (3):
  //[OBSOLETE]       os << "FORCE CONSTANT MATRICES in N/m:" << std::endl;
  //[OBSOLETE]       conversionFactor = 16.0217733;
  //[OBSOLETE]       break;
  //[OBSOLETE]   }
  //[OBSOLETE]   os << std::endl;

  //[OBSOLETE]   for (int i = 0; i < _supercell->getNumberOfAtoms(); i++) {
  //[OBSOLETE]     for (int k = 0; k < _supercell->getNumberOfAtoms(); k++) {
  //[OBSOLETE]       os << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  //[OBSOLETE]       os << setprecision(4);
  //[OBSOLETE]       os << "- MATRIX: " << i + 1 << "/" << k + 1 << " " << k + 1 << "/" << i + 1 << std::endl;
  //[OBSOLETE]       for (int m = 1; m <= 3; m++) {
  //[OBSOLETE]         for (int n = 1; n <= 3; n++)
  //[OBSOLETE]           os << setw(10) << (conversionFactor * _forceConstantMatrices[k][i](m, n)) << " ";
  //[OBSOLETE]         os << " ";
  //[OBSOLETE]         for (int n = 1; n <= 3; n++)
  //[OBSOLETE]           os << setw(10) << (conversionFactor * _forceConstantMatrices[i][k](n, m)) << " ";
  //[OBSOLETE]         os << std::endl;
  //[OBSOLETE]       }
  //[OBSOLETE]       os << std::endl;
  //[OBSOLETE]     }
  //[OBSOLETE]   }
  //[OBSOLETE] }

  //[OBSOLETE] // ///////////////////////////////////////////////////////////////////////////

  //[OBSOLETE] void ForceConstantCalculator::printFCShellInfo(ostream& os) {
  //[OBSOLETE]   int units = 4;
  //[OBSOLETE]   double conversionFactor = 1.0;

  //[OBSOLETE]   switch (units) {
  //[OBSOLETE]     case (1):
  //[OBSOLETE]       os << "FORCE CONSTANT MATRICES in eV/A^2:" << std::endl;
  //[OBSOLETE]       conversionFactor = 1.0;
  //[OBSOLETE]       break;
  //[OBSOLETE]     case (2):
  //[OBSOLETE]       os << "FORCE CONSTANT MATRICES in 10 Dyn/cm:" << std::endl;
  //[OBSOLETE]       conversionFactor = 1602.17733;
  //[OBSOLETE]       break;
  //[OBSOLETE]     case (3):
  //[OBSOLETE]       os << "FORCE CONSTANT MATRICES in N/m:" << std::endl;
  //[OBSOLETE]       conversionFactor = 16.0217733;
  //[OBSOLETE]       break;
  //[OBSOLETE]     case (4):
  //[OBSOLETE]       os << "FORCE CONSTANT MATRICES in 10^3 Dyn/cm:" << std::endl;
  //[OBSOLETE]       conversionFactor = 16.0217733;
  //[OBSOLETE]       break;
  //[OBSOLETE]   }
  //[OBSOLETE]   os << std::endl;

  //[OBSOLETE]   int maxshell = _supercell->getMaxShellID();
  //[OBSOLETE]   if (maxshell == -1) maxshell = 25;
  //[OBSOLETE]   std::vector<ShellHandle> sh;
  //[OBSOLETE]   for (int i = 0; i < _supercell->getNumberOfUniqueAtoms(); i++) {
  //[OBSOLETE]     ShellHandle s;
  //[OBSOLETE]     sh.push_back(s);
  //[OBSOLETE]     sh.back().init(_supercell->getInputStructure(),
  //[OBSOLETE]         _supercell->getInputStructure().iatoms[i][0],
  //[OBSOLETE]         maxshell);
  //[OBSOLETE]     sh[i].splitBySymmetry();
  //[OBSOLETE]     sh[i].mapStructure(_supercell->getSupercellStructure(), _supercell->getUniqueAtomID(i));
  //[OBSOLETE]   }

  //[OBSOLETE]   for (int i = 0; i < _supercell->getNumberOfUniqueAtoms(); i++) {
  //[OBSOLETE]     sh[i].printReport(cout);
  //[OBSOLETE]     for (int ishell = 0; ishell <= sh[i].getLastOccupiedShell(); ishell++) {
  //[OBSOLETE]       for (int isubshell = 0; isubshell < sh[i].getNumberOfSubshells(ishell); isubshell++) {
  //[OBSOLETE]         const deque<_atom>& atomsAtSameShell = sh[i].getAtomsAtSameShell(ishell, isubshell);
  //[OBSOLETE]         cout << "SHELL " << ishell << " " << isubshell << std::endl;

  //[OBSOLETE]         for (uint ai = 0; ai < atomsAtSameShell.size(); ai++) {
  //[OBSOLETE]           int nb = atomsAtSameShell[ai].number;
  //[OBSOLETE]           cout << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  //[OBSOLETE]           cout << setprecision(4);
  //[OBSOLETE]           cout << "- MATRIX: " << i << "/" << nb << " " << nb << "/" << i << std::endl;
  //[OBSOLETE]           for (int m = 1; m <= 3; m++) {
  //[OBSOLETE]             for (int n = 1; n <= 3; n++)
  //[OBSOLETE]               cout << setw(10) << (conversionFactor * _forceConstantMatrices[nb][i](m, n)) << " ";
  //[OBSOLETE]             cout << " ";
  //[OBSOLETE]             for (int n = 1; n <= 3; n++)
  //[OBSOLETE]               cout << setw(10) << (conversionFactor * _forceConstantMatrices[i][nb](n, m)) << " ";
  //[OBSOLETE]             cout << std::endl;
  //[OBSOLETE]           }
  //[OBSOLETE]           cout << std::endl;
  //[OBSOLETE]         }
  //[OBSOLETE]       }
  //[OBSOLETE]     }
  //[OBSOLETE]   }

  //[OBSOLETE]   // Clear
  //[OBSOLETE]   for (uint i = 0; i < sh.size(); i++)
  //[OBSOLETE]     sh[i].clear();
  //[OBSOLETE]   sh.clear();
  //[OBSOLETE] }

}  // namespace apl


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                              Direct Method                               //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                            VASP CALCULATIONS                             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  bool ForceConstantCalculator::runVASPCalculationsDM(_xinput& xInput, _aflags& _aflowFlags,
      _kflags& _kbinFlags, _xflags& _xFlags, string& _AflowIn) {
    stringstream message;
    bool stagebreak = false;

    // Determine the distortion vectors
    estimateUniqueDistortions(_supercell->getSupercellStructure(), _uniqueDistortions);

    //CO START
    vvgenerate_plus_minus.clear();  //CO //CO20181226  //ME20191029
    bool generate_plus_minus;           //CO
    //bool         check_minus_needed = ( AUTO_GENERATE_PLUS_MINUS && !USER_GENERATE_PLUS_MINUS && !_supercell.isDerivativeStructure() );
    //bool check_minus_needed = (AUTO_GENERATE_PLUS_MINUS && !USER_GENERATE_PLUS_MINUS);  OBSOLETE ME20181028 - this overrides DPM=OFF
    //[CO20181212]bool check_minus_needed = AUTO_GENERATE_PLUS_MINUS;  //ME20181028
    int ncalcs = 0;  //ME20190107 - total number of calculations for padding
    if (_calculateZeroStateForces) ncalcs++;  //ME20190112
    if (_isPolarMaterial) ncalcs++;  //ME20190112

    for (uint i = 0; i < _uniqueDistortions.size(); i++) {
      vvgenerate_plus_minus.push_back(vector<bool>(0)); //CO20181226
      for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
        //      vvgenerate_plus_minus.push_back(true);  //assume we need plus/minus OBSOLETE ME20181028 - this overrides DPM=OFF
        //ME20190107 - Calculate "need minus" here
        if (AUTO_GENERATE_PLUS_MINUS) {
          vvgenerate_plus_minus.back().push_back(needMinus(i, j, DISTORTION_INEQUIVONLY)); //CO20190218
        } else {
          vvgenerate_plus_minus.back().push_back(USER_GENERATE_PLUS_MINUS);  //ME20181028
        }
        if (vvgenerate_plus_minus[i][j]) {
          ncalcs += 2;
        } else {
          ncalcs += 1;
        }
      }
    }
    //CO END
    //ME20181022 START
    // Generate calculation directories
    for (uint i = 0; i < _uniqueDistortions.size(); i++) {
      for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
        //CO START
        //[CO20181212]if (check_minus_needed)
        //[ME20190107] if(AUTO_GENERATE_PLUS_MINUS)  //CHECK -
        //[ME20190107] {  //CO20200106 - patching for auto-indenting
        //[ME20190107]  vvgenerate_plus_minus[i][j] = needMinus(i, j);
        //[ME20190107]  if (!vvgenerate_plus_minus[i][j]) {_logger << "No negative distortion needed for distortion [atom=" << i << ",direction=" << j << "]." << apl::endl;}
        //[ME20190107]}
        generate_plus_minus = vvgenerate_plus_minus[i][j];
        if (AUTO_GENERATE_PLUS_MINUS && !generate_plus_minus) {
          message << "No negative distortion needed for distortion [atom=" << i << ",direction=" << j << "].";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
        }
        for (uint k = 0; k < (generate_plus_minus ? 2 : 1); k++) {
          //CO END
          // Copy settings from common case
          xInputs.push_back(xInput);
          uint idxRun = xInputs.size() - 1;
          uint idAtom = (DISTORTION_INEQUIVONLY ? _supercell->getUniqueAtomID(i) : i); //CO20190218

          // Create run ID
          //ME20190107 - added padding
          string runname = aurostd::PaddedNumString(idxRun + 1, aurostd::getZeroPadding(ncalcs)) + "_";  //ME20190112
          runname += "A" + aurostd::utype2string<int>(idAtom) + "D" + aurostd::utype2string<int>(j); //CO20190218

          if (generate_plus_minus) {  //CO
            runname = runname + ((k == 0) ? "P" : "M");
          }
          xInputs[idxRun].xvasp.AVASP_arun_runname = runname;
          // Apply the unique distortion to one inequvalent atom
          // This distortion vector is stored in Cartesian form, hence use C2F before applying
          xInputs[idxRun].setXStr(_supercell->getSupercellStructureLight()); //CO faster, only what's necessary here
          xstructure& xstr = xInputs[idxRun].getXStr(); //ME20190109 - Declare to make code more legible
          //CO20190114 - it is very silly to try to add in fpos
          //add to cpos, then convert to fpos
          xstr.atoms[idAtom].cpos += ((k == 0) ? 1.0 : -1.0) * DISTORTION_MAGNITUDE * _uniqueDistortions[i][j];
          xstr.atoms[idAtom].fpos = xstr.c2f * xstr.atoms[idAtom].cpos;
          //[CO20190114 - OBSOLETE]xInputs[idxRun].getXStr().atoms[idAtom].fpos = xInputs[idxRun].getXStr().atoms[idAtom].fpos + C2F(xInputs[idxRun].getXStr().lattice, ((k == 0) ? 1.0 : -1.0) * DISTORTION_MAGNITUDE * _uniqueDistortions[i][j]);
          //[CO20190114 - OBSOLETE]xInputs[idxRun].getXStr().atoms[idAtom].cpos = F2C(xInputs[idxRun].getXStr().lattice,
          //[CO20190114 - OBSOLETE]                                                 xInputs[idxRun].getXStr().atoms[idAtom].fpos);

          //clean title //CO20181226
          //[CO20190131 - moved up]xstructure& xstr = xInputs[idxRun].getXStr(); //ME20190109 - Declare to make code more legible
          xstr.title = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xstr.title); //CO20181226, ME20190109
          if(xstr.title.empty()){xstr.buildGenericTitle(true,false);} //CO20181226, ME20190109
          xstr.title += " APL supercell=";
          if (_supercell->scell_dim.rows == 3) {
            xstr.title += aurostd::joinWDelimiter(_supercell->scell_dim, "x"); //ME20190109
          } else {
            xmatrix<double> scmat = aurostd::xmatrixutype2double(aurostd::reshape(_supercell->scell_dim, 3, 3));
            xstr.title += "[" + xmatDouble2String(scmat, 0) + "]";
          }
          xstr.title += " atom=" + aurostd::utype2string<int>(idAtom); //ME20190109
          xvector<double> dist_cart = DISTORTION_MAGNITUDE * _uniqueDistortions[i][j];  //ME20190112
          xstr.title += " distortion=[" + aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(dist_cart, 3), ',') + "]";

          // For VASP, use the standardized aflow.in creator
          if (_kbinFlags.AFLOW_MODE_VASP){
            // Change format of POSCAR
            //ME20190228 - OBSOLETE for two reasons:
            // 1. This method is not robust
            // 2. This will be taken care of when the actual POSCAR is generated
            // [OBSOLETE - 20190228] if ((!_kbinFlags.KBIN_MPI && (_kbinFlags.KBIN_BIN.find("46") != string::npos)) ||
            // [OBSOLETE - 20190228]    (_kbinFlags.KBIN_MPI && (_kbinFlags.KBIN_MPI_BIN.find("46") != string::npos))) {
            // [OBSOLETE - 20190228]  xInputs[idxRun].getXStr().is_vasp5_poscar_format = false;
            // [OBSOLETE - 20190228] }
            stagebreak = (createAflowInPhonons(_aflowFlags, _kbinFlags, _xFlags, xInputs[idxRun]) || stagebreak);
          }
          // For AIMS, use the old method until we have AVASP_populateXAIMS
          if (_kbinFlags.AFLOW_MODE_AIMS) {
            string runname = ARUN_DIRECTORY_PREFIX + "APL_" + aurostd::utype2string<int>(idxRun) + "A" + aurostd::utype2string<uint>(_supercell->getUniqueAtomID(i)) + "D" + aurostd::utype2string<uint>(j);
            xInputs[idxRun].setDirectory(_directory + "/" + runname);
            if (!filesExistPhonons(xInputs[idxRun])) {
              message << "Creating " << xInputs[idxRun].getDirectory();
              pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
              createAflowInPhononsAIMS(_aflowFlags, _kbinFlags, _xFlags, _AflowIn, xInputs[idxRun], *p_FileMESSAGE);
              stagebreak = true;
            }
          }
        }
      }
    }

    // Add zero state if requested
    if (_calculateZeroStateForces) {
      // Copy settings from common case
      xInputs.push_back(xInput);
      uint idxRun = xInputs.size() - 1;
      // Create run ID //ME20181226
      xInputs[idxRun].xvasp.AVASP_arun_runname = aurostd::PaddedNumString(idxRun+1, aurostd::getZeroPadding(ncalcs)) + "_ZEROSTATE"; //ME20181226, ME20190112

      // Get structure
      xInputs[idxRun].setXStr(_supercell->getSupercellStructureLight()); //CO

      //ME20190108 - Set title
      xInputs[idxRun].getXStr().title=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xInputs[idxRun].getXStr().title); //CO20181226
      if(xInputs[idxRun].getXStr().title.empty()){xInputs[idxRun].getXStr().buildGenericTitle(true,false);} //CO20181226
      xInputs[idxRun].getXStr().title += " APL supercell=" ;
      if (_supercell->scell_dim.rows == 3) {
        xInputs[idxRun].getXStr().title += aurostd::joinWDelimiter(_supercell->scell_dim, "x"); //ME20190109
      } else {
        xmatrix<double> scmat = aurostd::xmatrixutype2double(aurostd::reshape(_supercell->scell_dim, 3, 3));
        xInputs[idxRun].getXStr().title += "[" + xmatDouble2String(scmat, 0) + "]";
      }
      xInputs[idxRun].getXStr().title += " undistorted";
      xInputs[idxRun].xvasp.aopts.flag("APL_FLAG::IS_ZEROSTATE", true);  //ME20191029
      // For VASP, use the standardized aflow.in creator
      if(_kbinFlags.AFLOW_MODE_VASP){
        stagebreak = (createAflowInPhonons(_aflowFlags, _kbinFlags, _xFlags, xInputs[idxRun]) || stagebreak); //ME20181226
      }
      // For AIMS, use the old method until we have AVASP_populateXAIMS //ME20181226
      if(_kbinFlags.AFLOW_MODE_AIMS){
        string runname = ARUN_DIRECTORY_PREFIX + "APL_" + aurostd::utype2string<uint>(idxRun) + "ZEROSTATE";
        xInputs[idxRun].setDirectory(_directory + "/" + runname);
        if (!filesExistPhonons(xInputs[idxRun])) {
          message << "Creating " << xInputs[idxRun].getDirectory();
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
          createAflowInPhononsAIMS(_aflowFlags, _kbinFlags, _xFlags, _AflowIn, xInputs[idxRun], *p_FileMESSAGE);
          stagebreak = true;
        }
      }
    }

    return stagebreak;
  }

  void ForceConstantCalculator::estimateUniqueDistortions(const xstructure& xstr,
      vector<vector<xvector<double> > >& uniqueDistortions) {
    //CO NOTES ON THIS FUNCTION
    // - this function creates symmetrically unique distortion vectors for each iatom
    // - you can have at most 3 unique (orthogonal) distortions per atom, but probably fewer considering symmetry
    // - distortions are relative to the lattice vectors
    // - we first generate different types of basic distortions (along lattice vectors, diagonal, body diagonal) in fractional coordinates
    // - convert fractional to cartesian, now distortions truly represent those along lattice vectors, etc.
    // - despite that the vector order changes with every loop, these basic distortions are not changed
    // - for each iatom,
    //   - we build three lists:
    //     - allDistortionsOfAtom (horrible name) is simply a list of symmetrically equivalent (orthogonal) distortions
    //       to each basic distortion (by rotations of the crystal) - this is true because we clear allDistortionsOfAtom with every loop,
    //       so original distortionVector is unchanged
    //     - uniqueDistortions does not play a role in this loop right now, so ignore
    //     - testVectorDim is a count of symmetrically equivalent distortions for each basic distortion
    //   - we sort the basic distortions by testVectorDim (largest first), i.e. the basic distortions that generate the highest count of
    //     equivalent distortions go first (you could consider it like the most "natural" choices for distortions based on crystal symmetry)

    stringstream message;
    // Is there a list of inequivalent atoms?
    if (DISTORTION_INEQUIVONLY && !xstr.iatoms_calculated) { //CO20190218
      message << "The list of the inequivalent atoms is missing.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    // Clear old stuff
    for (uint i = 0; i < uniqueDistortions.size(); i++) {
      uniqueDistortions[i].clear();
    }
    uniqueDistortions.clear();

    // Determine the distortion vectors for this system
    if (xstr.agroup_calculated) {
      // Test directions for distortion (in cartesian coordinates)
      vector<xvector<double> > testDistortions;
      xvector<double> testVec(3);

      if (GENERATE_ONLY_XYZ) {
        // Test distortion (in cartesian coordinates) along x, y, and z
        // - they have to be orthonormal
        testVec(1) = 1.0;
        testVec(2) = 0.0;
        testVec(3) = 0.0;
        testDistortions.push_back(testVec);
        testVec(1) = 0.0;
        testVec(2) = 1.0;
        testVec(3) = 0.0;
        testDistortions.push_back(testVec);
        testVec(1) = 0.0;
        testVec(2) = 0.0;
        testVec(3) = 1.0;
        testDistortions.push_back(testVec);
      } else {
        // Add lattice vectors
        testVec(1) = 1.0;
        testVec(2) = 0.0;
        testVec(3) = 0.0;
        testDistortions.push_back(testVec);
        testVec(1) = 0.0;
        testVec(2) = 1.0;
        testVec(3) = 0.0;
        testDistortions.push_back(testVec);
        testVec(1) = 0.0;
        testVec(2) = 0.0;
        testVec(3) = 1.0;
        testDistortions.push_back(testVec);

        // Add diagonal vectors
        testVec(1) = 1.0;
        testVec(2) = 1.0;
        testVec(3) = 0.0;
        testDistortions.push_back(testVec);
        testVec(1) = 1.0;
        testVec(2) = 0.0;
        testVec(3) = 1.0;
        testDistortions.push_back(testVec);
        testVec(1) = 0.0;
        testVec(2) = 1.0;
        testVec(3) = 1.0;
        testDistortions.push_back(testVec);

        // Add body diagonal vectors
        testVec(1) = 1.0;
        testVec(2) = 1.0;
        testVec(3) = 1.0;
        testDistortions.push_back(testVec);

        // Convert to cartesian representation
        for (uint idTestVector = 0; idTestVector < testDistortions.size(); idTestVector++) {
          testDistortions[idTestVector] = F2C(xstr.lattice, testDistortions[idTestVector]);
        }

        // Norm to one (so it will be scaled by user magnitude of distortion)
        for (uint idTestVector = 0; idTestVector < testDistortions.size(); idTestVector++) {
          double len = aurostd::modulus(testDistortions[idTestVector]);
          testDistortions[idTestVector](1) /= len;
          testDistortions[idTestVector](2) /= len;
          testDistortions[idTestVector](3) /= len;
        }
      }

      // Loop over inequivalent atoms
      vector<xvector<double> > uniqueDistortionsOfAtom;
      vector<xvector<double> > allDistortionsOfAtom;
      uint natoms = DISTORTION_INEQUIVONLY ? xstr.iatoms.size() : xstr.atoms.size();
      for (uint i = 0; i < natoms; i++) {
        int atomID = (DISTORTION_INEQUIVONLY ? xstr.iatoms[i][0] : i); //CO20190218
        //cout << "atomID = " << atomID << std::endl; //CO20190218
        if (xstr.agroup[atomID].size() == 0) { //CO20190218
          message << "Site point group operations are missing.";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
        }

        // Loop over test directions vectors - we count the number of unique
        // which can be obtain from the each test vector
        vector<int> testVectorDim;
        for (uint idTestVector = 0; idTestVector < testDistortions.size(); idTestVector++) {
          // Get the list of all equivalent distortions
          testDistortion(testDistortions[idTestVector],
              xstr.agroup[atomID], allDistortionsOfAtom,
              uniqueDistortionsOfAtom,true); //CO20190114 - maximize distortion impact first
          testVectorDim.push_back(allDistortionsOfAtom.size());
          allDistortionsOfAtom.clear();
          uniqueDistortionsOfAtom.clear();
        }

        // Now, we order all test vectors according to their dimensionality
        // Those with the highest score go first
        for (uint j = 0; j < testDistortions.size() - 1; j++) {
          for (uint k = j + 1; k < testDistortions.size(); k++) {
            if (testVectorDim[k] > testVectorDim[j]) {
              xvector<double> temp = testDistortions[k];
              testDistortions[k] = testDistortions[j];
              testDistortions[j] = temp;
              int itemp = testVectorDim[k];
              testVectorDim[k] = testVectorDim[j];
              testVectorDim[j] = itemp;
            }
          }
        }

        // Now we are going again, but slightly different, we count all
        // generated directions together until the total count is lower than three
        for (uint idTestVector = 0; idTestVector < testDistortions.size(); idTestVector++) {
          testDistortion(testDistortions[idTestVector],
              xstr.agroup[atomID], allDistortionsOfAtom,
              uniqueDistortionsOfAtom,DISTORTION_SYMMETRIZE);  //CO20190114 - if DISTORTION_SYMMETRIZE==false, we want all 3 directions (don't integrate equivalent distortions into count)
          if (allDistortionsOfAtom.size() >= 3) break;
        }

        //cout << "XXXXX  Number of unique distortion vectors for atom ["
        //     << atomID << "] = " << uniqueDistortionsOfAtom.size() << std::endl; //CO20190218
        uniqueDistortions.push_back(uniqueDistortionsOfAtom);
        // Free useless stuff
        allDistortionsOfAtom.clear();
        uniqueDistortionsOfAtom.clear();
      }
      // Free useless stuff
      testDistortions.clear();
    } else {
      message << "The list of the site point group operations is missing.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    // Print some information
    int dof = 0;
    for (uint i = 0; i < _uniqueDistortions.size(); i++){ //CO20200106 - wrapping with guard
      dof += _uniqueDistortions[i].size();
    }
    message << "Found " << dof << " degree(s) of freedom.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
    uint natoms = DISTORTION_INEQUIVONLY ? _supercell->getNumberOfUniqueAtoms() : _supercell->getNumberOfAtoms();
    for (uint i = 0; i < natoms; i++) {  //CO20200212 - int->uint
      uint id = (DISTORTION_INEQUIVONLY ? _supercell->getUniqueAtomID(i) : i); //CO20190218
      for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
        message << "Atom [" << aurostd::PaddedNumString(id, 3) << "] ("
          << std::setw(2) << _supercell->getSupercellStructure().atoms[id].cleanname
          << ") will be displaced in direction ["
          << std::fixed << std::setw(5) << std::setprecision(3) << _uniqueDistortions[i][j](1) << ","
          << std::fixed << std::setw(5) << std::setprecision(3) << _uniqueDistortions[i][j](2) << ","
          << std::fixed << std::setw(5) << std::setprecision(3) << _uniqueDistortions[i][j](3) << "].";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////

  void ForceConstantCalculator::testDistortion(const xvector<double>& distortionVector,
      const vector<_sym_op>& symPool,
      vector<xvector<double> >& allDistortionsOfAtom,
      vector<xvector<double> >& uniqueDistortionsOfAtom,
      bool integrate_equivalent_distortions) {  //CO20190114
    // Test if it is unique direction
    // Use the Gram–Schmidt method for orthogonalizing, if the final vectors
    // is non zero length -> is unique
    xvector<double> testDistortion = distortionVector;
    for (uint k = 0; k < allDistortionsOfAtom.size(); k++) {
      testDistortion = testDistortion - aurostd::getVectorProjection(testDistortion, allDistortionsOfAtom[k]);
    }
    if (aurostd::modulus(testDistortion) > _FLOAT_TOL_) {
      // Normalize vector
      testDistortion = testDistortion / aurostd::modulus(testDistortion);
      // Store this vector (in Cartesian form!)
      allDistortionsOfAtom.push_back(testDistortion);
      uniqueDistortionsOfAtom.push_back(testDistortion);
      //cout << "new unique distortion vector: " << testDistortion << std::endl;
    }

    if(integrate_equivalent_distortions) { //CO20181226
      // Apply all symmetry operations on test vector and generate all
      // unique directions and store them for future comparison
      for (uint iSymOp = 0; iSymOp < symPool.size(); iSymOp++) {
        testDistortion = symPool[iSymOp].Uc * distortionVector;

        for (uint k = 0; k < allDistortionsOfAtom.size(); k++) {
          testDistortion = testDistortion - aurostd::getVectorProjection(testDistortion, allDistortionsOfAtom[k]);
        }
        if (aurostd::modulus(testDistortion) > _FLOAT_TOL_) {
          testDistortion = testDistortion / aurostd::modulus(testDistortion);
          allDistortionsOfAtom.push_back(testDistortion);
          //cout << "new symmetry generated distortion vector: " << testDistortion << std::endl;
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////

  //CO START
  bool ForceConstantCalculator::needMinus(uint atom_index, uint distortion_index, bool inequiv_only) { //CO20190218
    //bool need_minus = true;
    const vector<_sym_op>& agroup = _supercell->getAGROUP( inequiv_only ? _supercell->getUniqueAtomID(atom_index) : atom_index);  //CO20190116
    //[CO20190131 OBSOLETE]uint atom_index = _supercell->getUniqueAtomID(ineq_atom_indx);
    //[CO20190131 OBSOLETE]vector<_sym_op> agroup = _supercell->getAGROUP(atom_index);

    xvector<double> rotated_distortion, distortion_sum;
    //cerr << agroup.size() << std::endl;
    //cerr << "distortion  : " << _uniqueDistortions[atom_index][distortion_index] << std::endl;
    for (uint i = 0; i < agroup.size(); i++) {
      rotated_distortion = agroup[i].Uc * _uniqueDistortions[atom_index][distortion_index];
      //cerr << "rdistortion : " << rotated_distortion << std::endl;
      if (identical(_uniqueDistortions[atom_index][distortion_index], -rotated_distortion, _FLOAT_TOL_)) {  //just mimicking Jahnatek tolerance here
        return FALSE;
        //need_minus = false;
        //break;
      } //else {
      //  continue;
      //}
    }
    return TRUE;
    //return need_minus;
  }
  //CO END

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                             FORCE CONSTANTS                              //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  bool ForceConstantCalculator::calculateForceConstantsDM() {
    // Get all forces required for the construction of force-constant matrices
    if (!calculateForceFields()) return false;

    //ME20191219 - atomGoesTo and atomComesFrom can now use basis_atoms_map.
    // Calculating the full basis ahead of time is much faster than calculating all
    // symmetry operations on-the-fly.
    if (!_supercell->fullBasisCalculatedAGROUP()) _supercell->getFullBasisAGROUP();

    // For construction of the force-constant matrices we need three
    // independent distortions. Hence, calculate the remaining distortions and
    // forces by the symmetry (if necessary).
    completeForceFields();

    // Ensure that all distortion vectors are along the cartesian directions
    projectToCartesianDirections();

    // Construct the matrix of force-constant matrices for all atoms based
    // on the force fields for the inequivalent atoms
    buildForceConstantMatrices();

    // Store data into DYNMAT file format - vasp like
    string filename = aurostd::CleanFileName(_directory + "/" + DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_DYNMAT_FILE);
    writeDYNMAT(filename);
    return true;
  }

  bool ForceConstantCalculator::calculateForceFields() {
    bool LDEBUG=(FALSE || _DEBUG_APL_HARM_IFCS_ || XHOST.DEBUG);
    // Extract all forces ////////////////////////////////////////////////////

    //first pass, just find if outfile is found ANYWHERE
    if(!outfileFoundAnywherePhonons(xInputs)) return false;

    //second pass, make sure it's everywhere!
    if (!outfileFoundEverywherePhonons(xInputs, _directory, *p_FileMESSAGE, *p_oss, _isPolarMaterial)) return false;

    // Remove zero state forces if necessary
    if (_calculateZeroStateForces) {
      subtractZeroStateForces(xInputs, _isPolarMaterial);
    }

    // Store forces //////////////////////////////////////////////////////////

    bool generate_plus_minus = false;  //ME20190129

    int idxRun = 0;
    uint natoms = DISTORTION_INEQUIVONLY ? _supercell->getNumberOfUniqueAtoms() : _supercell->getNumberOfAtoms();
    for (uint i = 0; i < natoms; i++) {  //CO20200212 - int->uint
      vector<vector<xvector<double> > > forcesForOneAtomAndAllDistortions;
      for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
        generate_plus_minus = vvgenerate_plus_minus[i][j];  //CO20181226
        vector<xvector<double> > forcefield;
        xvector<double> drift(3);
        for (uint k = 0; k < _supercell->getNumberOfAtoms(); k++) {
          xvector<double> force(3);
          force(1) = xInputs[idxRun].getXStr().qm_forces[k](1);
          force(2) = xInputs[idxRun].getXStr().qm_forces[k](2);
          force(3) = xInputs[idxRun].getXStr().qm_forces[k](3);
          if(LDEBUG) { cerr << __AFLOW_FUNC__ << " force[idistortion=" << i << ",atom=" << k << ",+]=" << xInputs[idxRun].getXStr().qm_forces[k] << std::endl;} //CO20190218
          if (generate_plus_minus) {  //CO
            force(1) = 0.5 * (force(1) - xInputs[idxRun + 1].getXStr().qm_forces[k](1));
            force(2) = 0.5 * (force(2) - xInputs[idxRun + 1].getXStr().qm_forces[k](2));
            force(3) = 0.5 * (force(3) - xInputs[idxRun + 1].getXStr().qm_forces[k](3));
            if(LDEBUG) { cerr << __AFLOW_FUNC__ << " force[idistortion=" << i << ",atom=" << k << ",-]=" << xInputs[idxRun + 1].getXStr().qm_forces[k] << std::endl;} //CO20190218
          }
          if(LDEBUG) { cerr << __AFLOW_FUNC__ << " force[idistortion=" << i << ",atom=" << k << ",AVG]=" << force << std::endl;} //CO20190218
          forcefield.push_back(force);
          drift = drift + force;
        }

        if (generate_plus_minus) {  //CO
          idxRun += 2;
        } else {
          idxRun++;
        }

        // Remove drift
        if(LDEBUG) { cerr << __AFLOW_FUNC__ << " drift[idistortion=" << i << ",total]=" << drift << std::endl;} //CO20190218
        drift(1) = drift(1) / forcefield.size();
        drift(2) = drift(2) / forcefield.size();
        drift(3) = drift(3) / forcefield.size();
        if(LDEBUG) { cerr << __AFLOW_FUNC__ << " drift[idistortion=" << i << ",AVG]=" << drift << std::endl;} //CO20190218
        for (uint k = 0; k < forcefield.size(); k++) {
          forcefield[k] = forcefield[k] - drift;
          if(LDEBUG) { cerr << __AFLOW_FUNC__ << " force[idistortion=" << i << ",atom=" << k << ",-drift]=" << forcefield[k] << std::endl;} //CO20190218
        }

        // Store
        forcesForOneAtomAndAllDistortions.push_back(forcefield);
        forcefield.clear();
      }
      _uniqueForces.push_back(forcesForOneAtomAndAllDistortions);
      forcesForOneAtomAndAllDistortions.clear();
    }

    return true;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void ForceConstantCalculator::completeForceFields() {
    stringstream message;
    //CO START
    // Test of stupidity...
    if (_supercell->getEPS() == AUROSTD_NAN) {
      message << "Need to define symmetry tolerance.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    //CO END
    // Show info
    message << "Calculating the missing force fields by symmetry.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);

    // Let's go
    for (uint i = 0; i < (DISTORTION_INEQUIVONLY ? _supercell->getNumberOfUniqueAtoms() : _supercell->getNumberOfAtoms()); i++) { //CO20190218
      // We need to have 3 linearly independent distortions
      if (_uniqueDistortions[i].size() != 3) {
        vector<xvector<double> > allDistortionsOfAtom;
        vector<xvector<double> > testForce;
        vector<vector<xvector<double> > > forcePool;
        xvector<double> testVec(3), testVec0(3);

        int atomID = (DISTORTION_INEQUIVONLY ? _supercell->getUniqueAtomID(i) : i); //CO20190218
        const vector<_sym_op>& agroup = _supercell->getAGROUP(atomID); //CO20190218

        // Generate next independent distortion by symmetry operations...
        uint currentSizeDistortions = _uniqueDistortions[i].size(); //CO20190218
        for (uint idistor = 0; idistor < currentSizeDistortions; idistor++) { //CO20190218
          // Apply all symmetry operations and check if it is independent
          for (uint symOpID = 0; symOpID < agroup.size(); symOpID++) {
            const _sym_op& symOp = agroup[symOpID];

            testForce.clear();
            for (uint k = 0; k < _supercell->getNumberOfAtoms(); k++) {
              try {
                //ME20191219 - atomGoesTo now uses basis_atoms_map; keep translation option in case
                // the basis has not been calculated for some reason
                int l = _supercell->atomComesFrom(symOp, k, atomID, true); //CO20190218
                testForce.push_back(symOp.Uc * _uniqueForces[i][idistor][l]);
              } catch (aurostd::xerror& e) {
                message << "Mapping problem ? <-> " << k << ".";
                throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
              }
            }

            // Get distortion vector (it is in the cartesian form) and apply symmetry rotation
            testVec = symOp.Uc * _uniqueDistortions[i][idistor];

            // Orthogonalize new rotated distortion vector on all accepted distortion vectors
            for (uint k = 0; k < allDistortionsOfAtom.size(); k++) {
              for (uint l = 0; l < _supercell->getNumberOfAtoms(); l++) {
                testForce[l] = testForce[l] - aurostd::getModeratedVectorProjection(forcePool[k][l], testVec, allDistortionsOfAtom[k]);
              }
              testVec = testVec - aurostd::getVectorProjection(testVec, allDistortionsOfAtom[k]);
            }

            // If the remaining vector is non-zero length, it is new independent direction, hence store it...
            if (aurostd::modulus(testVec) > _FLOAT_TOL_) {
              // Normalize to unit length
              double testVectorLength = aurostd::modulus(testVec);
              for (uint l = 0; l < _supercell->getNumberOfAtoms(); l++) {
                testForce[l] = testForce[l] / testVectorLength;
              }
              testVec = testVec / testVectorLength;
              allDistortionsOfAtom.push_back(testVec);

              // We suppose the symOpID == 0 is E (Identity) operation, hence the first
              // independent vector is already the calculated vector, hence new forces are not need to store
              if (allDistortionsOfAtom.size() > 1) {
                // Store new distortion
                _uniqueDistortions[i].push_back(testVec);
                // Store new force field
                _uniqueForces[i].push_back(testForce);
              }

              // Store for next orthogonalization procedure
              forcePool.push_back(testForce);
            }
            if (_uniqueDistortions[i].size() == 3) break;
          }
          if (_uniqueDistortions[i].size() == 3) break;
        }
        allDistortionsOfAtom.clear();
        for (uint ii = 0; ii < forcePool.size(); ii++) forcePool[ii].clear();
        forcePool.clear();
      }

      // I hope this will never happen...
      if (_uniqueDistortions[i].size() != 3) {
        string message = "Cannot complete force fields by symmetry.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  void ForceConstantCalculator::projectToCartesianDirections() {
    bool LDEBUG=(FALSE || _DEBUG_APL_HARM_IFCS_ || XHOST.DEBUG);
    for (uint i = 0; i < (DISTORTION_INEQUIVONLY ? _supercell->getNumberOfUniqueAtoms() : _supercell->getNumberOfAtoms()); i++) { //CO20190218
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " looking at displaced atom[idistortion=" << i << "]" << std::endl;} //CO20190218
      // Construct transformation matrix A
      xmatrix<double> A(3, 3), U(3, 3);
      for (uint j = 0; j < 3; j++) {
        // Ensure it is unit length
        _uniqueDistortions[i][j] = _uniqueDistortions[i][j] / aurostd::modulus(_uniqueDistortions[i][j]);
        if(LDEBUG){ //CO20190218
          cerr << __AFLOW_FUNC__ << " checking if uniqueDistortion[" << i << "][" << j << "] is unit length: ";
          cerr << "modulus(" << _uniqueDistortions[i][j] << ")=" << aurostd::modulus(_uniqueDistortions[i][j]) << std::endl;
        }

        // Copy to rows of U^T matrix
        for (uint k = 1; k <= 3; k++) {
          U(j + 1, k) = _uniqueDistortions[i][j](k);
        }
      }
      A = inverse(U);
      //CO20190116 - I believe U is an orthonormal matrix, as it defines a 3D axis
      //hence A = trasp(U) as well (faster)
      //keep for now

      if(LDEBUG){ //CO20190218
        cerr << __AFLOW_FUNC__ << " distortion matrix U(distortion,direction):" << std::endl;
        cerr << U << std::endl;
        cerr << __AFLOW_FUNC__ << " inverse matrix A:" << std::endl;
        cerr << A << std::endl;
      }

      // Update unique distortion vectors
      //CO20190116 - using trasp(A) instead of A because _uniqueDistortions[i][0] is a vector, not a matrix (as m is below)
      // we are really applying A * U == I,
      // so use A below (not trasp(A))
      _uniqueDistortions[i][0] = trasp(A) * _uniqueDistortions[i][0];
      _uniqueDistortions[i][1] = trasp(A) * _uniqueDistortions[i][1];
      _uniqueDistortions[i][2] = trasp(A) * _uniqueDistortions[i][2];

      if(LDEBUG){ //CO20190218
        cerr << __AFLOW_FUNC__ << " new cartesian-direction-projected uniqueDistortions[" << i << "][0]=" << _uniqueDistortions[i][0] << std::endl;
        cerr << __AFLOW_FUNC__ << " new cartesian-direction-projected uniqueDistortions[" << i << "][1]=" << _uniqueDistortions[i][1] << std::endl;
        cerr << __AFLOW_FUNC__ << " new cartesian-direction-projected uniqueDistortions[" << i << "][2]=" << _uniqueDistortions[i][2] << std::endl;
        //CO20190116 - cerr << __AFLOW_FUNC__ << " testing: trasp(A) * U should give same as above: trasp(A) * U = " << std::endl;  //U ~ m below
        //CO20190116 - cerr << trasp(A) * U << std::endl;
        cerr << __AFLOW_FUNC__ << " testing: A * U should give same as above: A * U = " << std::endl;  //U ~ m below //DUH A = inverse(U), so A*U = I
        cerr << A * U << std::endl;
      }

      // Update forces
      xmatrix<double> m(3, 3);
      for (uint j = 0; j < _supercell->getNumberOfAtoms(); j++) {
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " looking at supercell atom[" << j << "]" << std::endl;} //CO20190218
        for (int k = 0; k < 3; k++)
          for (int l = 1; l <= 3; l++)
            m(k + 1, l) = _uniqueForces[i][k][j](l);
        if(LDEBUG){ //CO20190218
          cerr << __AFLOW_FUNC__ << " BEFORE m = " << std::endl;
          cerr << m << std::endl;
        }
        // m = A * m * U; ??? I am not sure...
        m = A * m;
        // m = trasp(A) * m;  //CO NEW, treat forces exactly as distortion //CO20190116 - wrong, see above, trasp(A) is only for vectors
        if(LDEBUG){ //CO20190218
          cerr << __AFLOW_FUNC__ << " AFTER m = " << std::endl;
          cerr << m << std::endl;
        }
        for (int k = 0; k < 3; k++)
          for (int l = 1; l <= 3; l++)
            _uniqueForces[i][k][j](l) = m(k + 1, l);
      }
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  void ForceConstantCalculator::buildForceConstantMatrices() {
    bool LDEBUG=(FALSE || _DEBUG_APL_HARM_IFCS_ || XHOST.DEBUG);
    stringstream message;
    // Test of stupidity...
    if (DISTORTION_INEQUIVONLY && !_supercell->getSupercellStructure().fgroup_calculated) { //CO20190218
      message << "The factor group has not been calculated yet.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    //CO START
    if (DISTORTION_INEQUIVONLY && _supercell->getEPS() == AUROSTD_NAN) { //CO20190218
      message << "Need to define symmetry tolerance.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    //CO END

    // Clear old matrices
    for (uint i = 0; i < _forceConstantMatrices.size(); i++)
      _forceConstantMatrices[i].clear();
    _forceConstantMatrices.clear();

    //CO20190116 - BIG BUG
    //do NOT push_back() with forceConstantMatrices
    //Jahnatek assumed that iatoms structure was [0 1 2 3] [4 5 6] (in order)
    //therefore, pushing back meant keeping forceConstantMatrices in order of supercell atoms
    //this is not necessarily true, as the mappings could be out of order
    //therefore, we create the vector of the necessary dimensions, and put the row in the right place
    //CO20190131 UPDATE - this is NOT the only part of the code for which this dependency (iatoms sorted) exists

    for (uint i = 0; i < _supercell->getNumberOfAtoms(); i++) {
      _forceConstantMatrices.push_back(vector<xmatrix<double> >(0));
      for (uint j = 0; j < _supercell->getNumberOfAtoms(); j++) {
        _forceConstantMatrices.back().push_back(xmatrix<double>(3,3,1,1));
      }
    }


    //
    message << "Calculating the force constant matrices.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);

    // We have a party. Let's fun with us...
    //vector<xmatrix<double> > row; //JAHNATEK ORIGINAL //CO20190218
    for (uint i = 0; i < (DISTORTION_INEQUIVONLY ? _supercell->getNumberOfUniqueAtoms() : _supercell->getNumberOfAtoms()); i++) { //CO20190218
      // Get the number of this atom in the whole list
      int basedAtomID = (DISTORTION_INEQUIVONLY ? _supercell->getUniqueAtomID(i) : i); //CO20190218

      // This is easy. We know everything. Just construct a set of matrices.
      xmatrix<double> m(3, 3, 1, 1);
      for (uint j = 0; j < _supercell->getNumberOfAtoms(); j++) {
        for (uint k = 0; k < _uniqueDistortions[i].size(); k++) {
          double distortionLength = aurostd::modulus(DISTORTION_MAGNITUDE * _uniqueDistortions[i][k]);
          // FCM element = -F/d, but we will omit minus, because next force transformations are better
          // done without it, and in construction of dyn. matrix we will add it to the sum
          //ME20200212 - we store _forceConstantMatrices in a human-readable file, so they should
          // represent the actual force constants, not an AFLOW-customized construct
          m(k + 1, 1) = _uniqueForces[i][k][j](1) / distortionLength;
          m(k + 1, 2) = _uniqueForces[i][k][j](2) / distortionLength;
          m(k + 1, 3) = _uniqueForces[i][k][j](3) / distortionLength;
          //cout << i << " " << k << " " << j << " "; printXVector(_uniqueForces[i][k][j]);
        }
        //printXMatrix(m);
        //row.push_back(m); //JAHNATEK ORIGINAL //CO20190218
        _forceConstantMatrices[basedAtomID][j] = m;  //CO NEW //CO20190218
        if(LDEBUG){ //CO20190218
          cerr << __AFLOW_FUNC__ << " adding m to forceConstantMatrices[" << basedAtomID << "][" << j << "]=" << std::endl;
          cerr << m << std::endl;
        }
      }
      //_forceConstantMatrices.push_back(row);  //JAHNATEK ORIGINAL //CO20190218
      //row.clear();  //JAHNATEK ORIGINAL //CO20190218

      if(DISTORTION_INEQUIVONLY){ //CO20190218
        _sym_op symOp;  //CO
        // Calculate rows for next equivalent atoms starting 1 (structure of iatoms)... //CO20190218
        for (uint j = 1; j < _supercell->getNumberOfEquivalentAtomsOfType(i); j++) { //CO20190218
          try {
            //CO20190116 - we want to map the forces of the inequivalent atoms (for which we ran vasp) onto the equivalent ones
            //hence, we need the FGROUP that takes us from the inequivalent atom to the equivalent
            //then, we need to find the atom which, upon application of that symop, becomes k (below)
            //symOp = _supercell->getSymOpWhichMatchAtoms(basedAtomID, _supercell->getUniqueAtomID(i, j), _FGROUP_);  //CO NEW //CO20190218
            symOp = _supercell->getSymOpWhichMatchAtoms(_supercell->getUniqueAtomID(i, j), basedAtomID, _FGROUP_);  //JAHNATEK ORIGINAL //CO20190218
            //const _sym_op& symOp = _supercell->getSymOpWhichMatchAtoms(_supercell->getUniqueAtomID(i,j),basedAtomID,_FGROUP_); //JAHNATEK ORIGINAL
            //cout << basedAtomID << " -> " << _supercell->getUniqueAtomID(i,j) << " " << symOp.str_type << " shift:"; printXVector(symOp.ftau);
            //printXVector(_supercell->getSupercellStructure().atoms[basedAtomID].fpos);
            //printXVector(_supercell->getSupercellStructure().atoms[_supercell->getUniqueAtomID(i,j)].fpos);
          } catch (aurostd::xerror& e) {
            message << "Mapping problem " << _supercell->getUniqueAtomID(i, j) << " <-> " << basedAtomID << "?"; //CO20190218
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
          }

          for (uint k = 0; k < _supercell->getNumberOfAtoms(); k++) {
            try {
              //CO20190116 - read atomComesFrom() as: applying symOp to l makes k
              // int l = _supercell.atomComesFrom(symOp, k, _supercell->getUniqueAtomID(i, j));  //CO NEW //CO20190218
              int l = _supercell->atomGoesTo(symOp, k, _supercell->getUniqueAtomID(i, j)); //JAHNATEK ORIGINAL //CO20190218
              //cout << "MAP " << k << " <-> " << l << std::endl;
              //row.push_back(inverse(symOp.Uc) * _forceConstantMatrices[basedAtomID][l] * symOp.Uc); //JAHNATEK ORIGINAL //CO20190218
              //row.push_back(symOp.Uc * _forceConstantMatrices[basedAtomID][l] * inverse(symOp.Uc)); //CO NEW  //JAHNATEK ORIGINAL //CO20190218
              //m = symOp.Uc * _forceConstantMatrices[basedAtomID][l] * inverse(symOp.Uc);  //CO NEW //CO20190218
              m = inverse(symOp.Uc) * _forceConstantMatrices[basedAtomID][l] * symOp.Uc;  //JAHNATEK ORIGINAL //CO20190218
              _forceConstantMatrices[_supercell->getUniqueAtomID(i, j)][k] = m;  //CO NEW //CO20190218
              if(LDEBUG){ //CO20190218
                cerr << __AFLOW_FUNC__ << " adding m to forceConstantMatrices[" << _supercell->getUniqueAtomID(i, j) << "][" << k << "]=" << std::endl;
                cerr << m << std::endl;
              }
            } catch (aurostd::xerror& e) {
              message << "Mapping problem " << k << " <-> ?.";
              throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Mapping failed.");
            }
          }
          //_forceConstantMatrices.push_back(row);  //JAHNATEK ORIGINAL //CO20190218
          //row.clear();  //JAHNATEK ORIGINAL //CO20190218
        }
      }
      //row.clear();  //JAHNATEK ORIGINAL //CO20190218
    }

    // Test of correctness
    if (_forceConstantMatrices.size() != _supercell->getNumberOfAtoms()) {
      message << "Some problem with the application of factor group operations.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    //ME20200211 - force constants are -F/d, not F/d
    for (uint i = 0; i < _forceConstantMatrices.size(); i++) {
      for (uint j = 0; j < _forceConstantMatrices.size(); j++) {
        _forceConstantMatrices[i][j] = -_forceConstantMatrices[i][j];
      }
    }
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                               FILE OUTPUT                                //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  // Writes the forces into a VASP DYNMAT format
  void ForceConstantCalculator::writeDYNMAT(const string& filename) {
    string message = "Writing forces into file " + filename + ".";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);

    stringstream outfile;

    // 1st line
    outfile << _supercell->getNumberOfUniqueAtoms() << " ";
    outfile << _supercell->getNumberOfAtoms() << " ";
    int dof = 0;
    for (uint i = 0; i < _uniqueDistortions.size(); i++)
      dof += _uniqueDistortions[i].size();
    outfile << dof << std::endl;

    // 2nd line
    outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    outfile << setprecision(3);
    for (uint i = 0; i < _supercell->getNumberOfUniqueAtoms(); i++) {
      if (i != 0) outfile << " ";
      outfile << _supercell->getUniqueAtomMass(i);
    }
    outfile << std::endl;

    // forces + 1 line info about distortion
    for (uint i = 0; i < _supercell->getNumberOfUniqueAtoms(); i++) {
      for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
        // line info
        outfile << (_supercell->getUniqueAtomID(i) + 1) << " ";
        outfile << (j + 1) << " ";
        xvector<double> shift(3);
        shift = DISTORTION_MAGNITUDE * _uniqueDistortions[i][j];
        outfile << setprecision(3);
        outfile << shift(1) << " " << shift(2) << " " << shift(3) << std::endl;
        // forces
        outfile << setprecision(6);
        for (uint k = 0; k < _supercell->getNumberOfAtoms(); k++)
          outfile << setw(15) << _uniqueForces[i][j][k](1)
            << setw(15) << _uniqueForces[i][j][k](2)
            << setw(15) << _uniqueForces[i][j][k](3) << std::endl;
      }
    }

    aurostd::stringstream2file(outfile, filename);
    if (!aurostd::FileExist(filename)) {
      message = "Cannot open output file " + filename + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
  }

  // [OBSOLETE]  //////////////////////////////////////////////////////////////////////////////

  // OBSOLETE ME20200504 - Not used
  // [OBSOLETE] // This is the interface to phonopy code

  // [OBSOLETE] void ForceConstantCalculator::writeFORCES() {
  // [OBSOLETE]   string message = "Writing forces into file FORCES.";
  // [OBSOLETE]   pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);

  // [OBSOLETE]   xstructure ix;
  // [OBSOLETE]   string filename = "SPOSCAR";
  // [OBSOLETE]   if (!aurostd::FileEmpty(filename)) {
  // [OBSOLETE]     message = "Reading " + filename;
  // [OBSOLETE]     pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
  // [OBSOLETE]     stringstream SPOSCAR;
  // [OBSOLETE]     aurostd::efile2stringstream(filename, SPOSCAR);
  // [OBSOLETE]     SPOSCAR >> ix;
  // [OBSOLETE]   } else {
  // [OBSOLETE]     ix = _supercell->getSupercellStructure();
  // [OBSOLETE]   }

  // [OBSOLETE]   stringstream outfile;

  // [OBSOLETE]   // 1st line
  // [OBSOLETE]   int dof = 0;
  // [OBSOLETE]   for (uint i = 0; i < _uniqueDistortions.size(); i++)
  // [OBSOLETE]     dof += _uniqueDistortions[i].size();
  // [OBSOLETE]   outfile << dof << std::endl;

  // [OBSOLETE]   // forces + 1 line info about distortion
  // [OBSOLETE]   outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  // [OBSOLETE]   for (int i = 0; i < _supercell->getNumberOfUniqueAtoms(); i++) {
  // [OBSOLETE]     for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
  // [OBSOLETE]       // line info
  // [OBSOLETE]       outfile << (_supercell->getUniqueAtomID(i) + 1) << " ";
  // [OBSOLETE]       xvector<double> shift(3);
  // [OBSOLETE]       shift = C2F(_supercell->getSupercellStructure().lattice, DISTORTION_MAGNITUDE * _uniqueDistortions[i][j]);
  // [OBSOLETE]       outfile << setprecision(6);
  // [OBSOLETE]       outfile << shift(1) << " " << shift(2) << " " << shift(3) << std::endl;
  // [OBSOLETE]       // forces
  // [OBSOLETE]       outfile << setprecision(6);
  // [OBSOLETE]       for (int k = 0; k < _supercell->getNumberOfAtoms(); k++) {
  // [OBSOLETE]         int l = 0;
  // [OBSOLETE]         for (; l < _supercell->getNumberOfAtoms(); l++)
  // [OBSOLETE]           if ((aurostd::abs(ix.atoms[k].cpos(1) - _supercell->getSupercellStructure().atoms[l].cpos(1)) < _FLOAT_TOL_) &&
  // [OBSOLETE]               (aurostd::abs(ix.atoms[k].cpos(2) - _supercell->getSupercellStructure().atoms[l].cpos(2)) < _FLOAT_TOL_) &&
  // [OBSOLETE]               (aurostd::abs(ix.atoms[k].cpos(3) - _supercell->getSupercellStructure().atoms[l].cpos(3)) < _FLOAT_TOL_))
  // [OBSOLETE]             break;
  // [OBSOLETE]         //CO, not really mapping error, just mismatch between structure read in (ix) and current supercell structure (should be exact)
  // [OBSOLETE]         if (l == _supercell->getNumberOfAtoms()) {
  // [OBSOLETE]           cout << k << std::endl;
  // [OBSOLETE]           cout << ix.atoms[k].fpos << std::endl;
  // [OBSOLETE]           cout << ix.atoms[k].cpos << std::endl;
  // [OBSOLETE]           message = "Mapping error.";
  // [OBSOLETE]           throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
  // [OBSOLETE]         }

  // [OBSOLETE]         outfile << setw(15) << _uniqueForces[i][j][l](1) << " "
  // [OBSOLETE]           << setw(15) << _uniqueForces[i][j][l](2) << " "
  // [OBSOLETE]           << setw(15) << _uniqueForces[i][j][l](3) << std::endl;
  // [OBSOLETE]       }
  // [OBSOLETE]     }
  // [OBSOLETE]   }

  // [OBSOLETE]   filename = "FORCES";
  // [OBSOLETE]   aurostd::stringstream2file(outfile, filename);
  // [OBSOLETE]   if (!aurostd::FileExist(filename)) {
  // [OBSOLETE]     message = "Cannot open output file.";
  // [OBSOLETE]     throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
  // [OBSOLETE]   }

  // [OBSOLETE] }

  // [OBSOLETE] // ///////////////////////////////////////////////////////////////////////////

  // [OBSOLETE] void ForceConstantCalculator::writeXCrysDenForces() {
  // [OBSOLETE]   string message = "Writing forces into file XCrysDenForces.";
  // [OBSOLETE]   pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
  // [OBSOLETE]   _supercell->center_original();  //CO

  // [OBSOLETE]   stringstream outfile;  //CO
  // [OBSOLETE]   // forces + 1 line info about distortion
  // [OBSOLETE]   for (int i = 0; i < _supercell->getNumberOfUniqueAtoms(); i++) {
  // [OBSOLETE]     for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
  // [OBSOLETE]       //string s = "FORCES_A" + aurostd::utype2string<uint>(_supercell->getUniqueAtomID(i)) + "D" + aurostd::utype2string<uint>(j) + ".xsf"; //CO
  // [OBSOLETE]       outfile.str("");  //CO
  // [OBSOLETE]       //ofstream outfile(s.c_str(), ios_base::out); //CO

  // [OBSOLETE]       outfile << "CRYSTAL" << std::endl;
  // [OBSOLETE]       outfile << "PRIMVEC 1" << std::endl;
  // [OBSOLETE]       outfile << _supercell->getSupercellStructure().lattice << std::endl;
  // [OBSOLETE]       outfile << "CONVEC 1" << std::endl;
  // [OBSOLETE]       outfile << _supercell->getSupercellStructure().lattice << std::endl;
  // [OBSOLETE]       outfile << "PRIMCOORD 1" << std::endl;
  // [OBSOLETE]       outfile << _supercell->getNumberOfAtoms() << " 1" << std::endl;

  // [OBSOLETE]       xvector<double> shift(3);
  // [OBSOLETE]       shift = C2F(_supercell->getSupercellStructure().lattice, DISTORTION_MAGNITUDE * _uniqueDistortions[i][j]);

  // [OBSOLETE]       outfile << setprecision(6);
  // [OBSOLETE]       for (int k = 0; k < _supercell->getNumberOfAtoms(); k++) {
  // [OBSOLETE]         outfile << _supercell->getAtomNumber(k) << " ";
  // [OBSOLETE]         outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  // [OBSOLETE]         outfile << setprecision(8);
  // [OBSOLETE]         xvector<double> r = F2C(_supercell->getSupercellStructure().lattice,
  // [OBSOLETE]             _supercell->getSupercellStructure().atoms[k].fpos);
  // [OBSOLETE]         outfile << setw(15) << r(1) << setw(15) << r(2) << setw(15) << r(3) << " ";
  // [OBSOLETE]         // this is strange...
  // [OBSOLETE]         //outfile << setw(15) << _superCellStructure.atoms[k].cpos << " ";

  // [OBSOLETE]         // Scale force, it is expected in Hartree/Angs.
  // [OBSOLETE]         xvector<double> f = hartree2eV * _uniqueForces[i][j][k];

  // [OBSOLETE]         outfile << setw(15) << f(1)
  // [OBSOLETE]           << setw(15) << f(2)
  // [OBSOLETE]           << setw(15) << f(3) << std::endl;
  // [OBSOLETE]       }

  // [OBSOLETE]       string filename = "FORCES_A" + aurostd::utype2string<uint>(_supercell->getUniqueAtomID(i)) + "D" + aurostd::utype2string<uint>(j) + ".xsf";
  // [OBSOLETE]       aurostd::stringstream2file(outfile, filename);
  // [OBSOLETE]       if (!aurostd::FileExist(filename)) {
  // [OBSOLETE]         string message = "Cannot create " + filename + " file.";
  // [OBSOLETE]         throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
  // [OBSOLETE]       }

  // [OBSOLETE]     }
  // [OBSOLETE]   }
  // [OBSOLETE] }

}  // namespace apl


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                            Linear Response                               //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  //////////////////////////////////////////////////////////////////////////////
  // We will use VASP5.2+ to calculate Born effective charge tensors and the
  // dielectric constant matrix in the primitive cell with very high precision
  // Both values are needed by the non-analytical term of the dynamical matrix
  // to capture the TO-LO splitting of optical phonon branches of polar systems.
  bool ForceConstantCalculator::runVASPCalculationsLR(_xinput& xInput, _aflags& _aflowFlags,
      _kflags& _kbinFlags, _xflags& _xFlags, string& _AflowIn) {
    bool stagebreak = false;

    xInput.setXStr(_supercell->getSupercellStructureLight());
    xInput.getXStr().title=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xInput.getXStr().title);
    if(xInput.getXStr().title.empty()){xInput.getXStr().buildGenericTitle(true,false);}
    xInput.getXStr().title+=" linear response";

    // For VASP, use the standardized aflow.in creator
    if(xInput.AFLOW_MODE_VASP) {
      xInput.xvasp.AVASP_arun_runname = "1_" + _AFLOW_APL_DFPT_RUNNAME_;  //ME20200213
      xInput.xvasp.aopts.flag("APL_FLAG::AVASP_BORN", false);
      xInput.xvasp.aopts.flag("APL_FLAG::AVASP_LR", true);

      // Set POSCAR to VASP5 format
      xInput.getXStr().is_vasp4_poscar_format = false;
      xInput.getXStr().is_vasp5_poscar_format = true;
      stagebreak = (createAflowInPhonons(_aflowFlags, _kbinFlags, _xFlags, xInput) || stagebreak);
    }
    // For AIMS, use the old method until we have AVASP_populateXAIMS
    if(xInput.AFLOW_MODE_AIMS) {
      string runname = _AFLOW_APL_DFPT_DIRECTORY_NAME_;  //ME20200213
      xInput.setDirectory(_directory + "/" + runname);
      if (!filesExistPhonons(xInput)) {
        string message = "Creating " + xInput.getDirectory();
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
        createAflowInPhononsAIMS(_aflowFlags, _kbinFlags, _xFlags, _AflowIn, xInput, *p_FileMESSAGE);
        stagebreak = true;
      }
    }
    return stagebreak;
  }

  //////////////////////////////////////////////////////////////////////////////
  //ME20200211
  bool ForceConstantCalculator::readForceConstantsFromVasprun(_xinput& xinp) {
    stringstream message;
    message << "Reading force constants from vasprun.xml";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);

    // Read vasprun.xml
    string filename = aurostd::CleanFileName(xinp.getDirectory() + "/vasprun.xml.static");
    if (!aurostd::EFileExist(filename)) {
      filename = aurostd::CleanFileName(xinp.getDirectory() + "/vasprun.xml");
      if (aurostd::EFileExist(filename)) {
        message << "Could not find vasprun.xml file for linear response calculations.";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
        return false;
      }
    }
    vector<string> vlines;
    aurostd::efile2vectorstring(filename, vlines);
    uint nlines = vlines.size();

    // Read Hessian
    vector<vector<double> > hessian;
    vector<double> row;
    vector<string> line;
    uint iline = 0;
    for (iline = 0; iline < nlines; iline++) {
      if (aurostd::substring2bool(vlines[iline], "hessian")) {
        while ((++iline != nlines) && !aurostd::substring2bool(vlines[iline], "</varray>")) {
          aurostd::string2tokens(vlines[iline], line);
          row.clear();
          for (uint i = 1; i < (line.size() - 1); i++) {
            row.push_back(aurostd::string2utype<double>(line[i]));
          }
          hessian.push_back(row);
        }
        if (iline < nlines) break;
      }
    }

    // Check that the file was read successfully.
    if (iline == nlines) {
      message << "Hessian tag not found or incomplete.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    }
    uint natoms = _supercell->getSupercellStructure().atoms.size();
    uint nhessian = hessian.size();
    if (nhessian != 3 * natoms) {
      message << "Hessian matrix does not have the correct number of rows (has "
        << nhessian << ", should have " << (3 * natoms) << ")." << std::endl;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    }
    uint i = 0;
    for (i = 0; i < nhessian; i++) {
      if (hessian[i].size() != 3 * natoms) break;
    }
    if (i != nhessian) {
      message << "Row " << (i + 1) << " of the Hessian matrix does not have the correct number of columns"
        << " (has " << hessian[i].size() << ", should have " << (3 * natoms) << ")." << std::endl;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    }

    _forceConstantMatrices.clear();
    _forceConstantMatrices.resize(natoms, vector<xmatrix<double> >(natoms, xmatrix<double>(3, 3)));
    double mass = 0.0;
    // Convert Hessian matrix into force constants
    for (uint i = 0; i < natoms; i++) {
      for (uint j = 0; j < natoms; j++) {
        // Hessian matrix is normalized by masses, so multiply to get IFCs
        mass = std::sqrt(_supercell->getAtomMass(i) * _supercell->getAtomMass(j));
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            _forceConstantMatrices[i][j][k+1][l+1] = -mass * hessian[3 * i + k][3 * j + l];
          }
        }
      }
    }
    return true;
  }

}  // namespace apl

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
