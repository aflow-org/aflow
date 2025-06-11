// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

#include "aflow_apl.h"

#define _DEBUG_APL_PHONCALC_ false  //CO20190116

using std::string;
using std::vector;
using aurostd::xvector;
using namespace std::placeholders;

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                         CONSTRUCTORS/DESTRUCTORS                         //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  PhononCalculator::PhononCalculator(ostream& oss) : xStream(oss) {
    free();
    _qm = QMesh(oss);
    _supercell = Supercell(oss);
    setDirectory("./");
    _ncpus = 1;
  }

  PhononCalculator::PhononCalculator(ofstream& mf, ostream& oss) : xStream(mf,oss) {
    free();
    _qm = QMesh(mf, oss);
    _supercell = Supercell(mf, oss);
    setDirectory("./");
    _ncpus = 1;
  }

  PhononCalculator::PhononCalculator(const PhononCalculator& that) : xStream(*that.getOFStream(),*that.getOSS()) {
    if (this != &that) free();
    copy(that);
  }

  // Copy constructors
  PhononCalculator& PhononCalculator::operator=(const PhononCalculator& that) {
    if (this != &that) free();
    copy(that);
    return *this;
  }

  void PhononCalculator::copy(const PhononCalculator& that) {
    if (this == &that) return;
    xStream::copy(that);
    _bornEffectiveChargeTensor = that._bornEffectiveChargeTensor;
    _dielectricTensor = that._dielectricTensor;
    _directory = that._directory;
    _forceConstantMatrices = that._forceConstantMatrices;
    _gammaEwaldCorr = that._gammaEwaldCorr;
    _inverseDielectricTensor = that._inverseDielectricTensor;
    _isGammaEwaldPrecomputed = that._isGammaEwaldPrecomputed;
    _ncpus = that._ncpus;
    _qm = that._qm;
    _recsqrtDielectricTensorDeterminant = that._recsqrtDielectricTensorDeterminant;
    _supercell = that._supercell;
    _system = that._system;
  }

  PhononCalculator::~PhononCalculator() {
    xStream::free();
    free();
  }

  void PhononCalculator::free() {
    _bornEffectiveChargeTensor.clear();
    _dielectricTensor.clear();
    _directory = "";
    _forceConstantMatrices.clear();
    _gammaEwaldCorr.clear();
    _inverseDielectricTensor.clear();
    _isGammaEwaldPrecomputed = false;
    _isPolarMaterial = false;
    _ncpus = 0;
    _recsqrtDielectricTensorDeterminant = 0.0;
    _qm.clear();
    _system = "";
    _supercell.clear();
  }


  void PhononCalculator::clear() {
    free();
  }

  //xStream initializers
  void PhononCalculator::initialize(ostream& oss) {
    _qm.initialize(oss);
    _supercell.initialize(oss);
  }

  void PhononCalculator::initialize(ofstream& mf, ostream& oss) {
    _qm.initialize(mf, oss);
    _supercell.initialize(mf, oss);
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                               INTERFACE                                  //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  QMesh& PhononCalculator::getQMesh() {
    return _qm;
  }

  Supercell& PhononCalculator::getSupercell() {
    return _supercell;
  }

  const xstructure& PhononCalculator::getInputCellStructure() const {
    return _supercell.getInputStructure();
  }

  const xstructure& PhononCalculator::getSuperCellStructure() const {
    return _supercell.getSupercellStructure();
  }

  uint PhononCalculator::getNumberOfBranches() const {
    return 3 * _supercell.getInputStructure().atoms.size();
  }

  //ME20190614
  string PhononCalculator::getDirectory() const {
    return _directory;
  }

  int PhononCalculator::getNCPUs() const {
    return _ncpus;
  }

  //ME20200206
  bool PhononCalculator::isPolarMaterial() const {
    return _isPolarMaterial;
  }

  const vector<vector<xmatrix<double> > >& PhononCalculator::getHarmonicForceConstants() const {
    return _forceConstantMatrices;
  }

  const vector<vector<double> >& PhononCalculator::getAnharmonicForceConstants(int order) const {
    return anharmonicIFCs[order - 3];
  }

  const vector<vector<int> >& PhononCalculator::getClusters(int order) const {
    return clusters[order - 3];
  }

}  // namespace apl

// ///////////////////////////////////////////////////////////////////////////

namespace apl {

  void PhononCalculator::setDirectory(const string& dir) {
    _directory = dir;
    _qm._directory = dir;
    _supercell._directory = dir;
  }

  void PhononCalculator::setNCPUs(const _kflags& kfl) {
    _ncpus = KBIN::get_NCPUS(kfl);
  }

  void PhononCalculator::setPolarMaterial(bool polar) {
    _isPolarMaterial = polar;
  }

}  // namespace apl

// ///////////////////////////////////////////////////////////////////////////

namespace apl {

  void PhononCalculator::initialize_qmesh(const vector<int>& grid, bool include_inversions, bool gamma_centered) {
    initialize_qmesh(aurostd::vector2xvector(grid), include_inversions, gamma_centered);
  }

  void PhononCalculator::initialize_qmesh(const xvector<int>& grid, bool include_inversions, bool gamma_centered) {
    _qm.initialize(grid, _supercell.getInputStructure(), include_inversions, gamma_centered);
  }

  void PhononCalculator::initialize_supercell(const xstructure& xstr, bool verbose) {
    _supercell.initialize(xstr, verbose);//AS20200908
  }
  void PhononCalculator::initialize_supercell(const string& filename) {
    _supercell.initialize(filename);
  }
}

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                            FORCE CONSTANTS                               //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  void PhononCalculator::setHarmonicForceConstants(const ForceConstantCalculator& fc) {
    setHarmonicForceConstants(fc.getForceConstants(), fc.getBornEffectiveChargeTensor(),
        fc.getDielectricTensor(), fc.isPolarMaterial());//AS20201208
  }

  // AS20201204 BEGIN
  void PhononCalculator::setHarmonicForceConstants(
      const vector<vector<xmatrix<double> > > &IFC,
      const vector<xmatrix<double> > &bornEffectiveChargeTensor,
      const xmatrix<double> &dielectricTensor,
      bool isPolar)
  {
    string message="";
    // check if the input IFCs have correct size
    uint natoms = _supercell.getNumberOfAtoms();
    if (IFC.size() != natoms){
      message = "The supplied IFC has the wrong size: ";
      message += aurostd::utype2string(IFC.size()) + " instead of ";
      message += aurostd::utype2string(natoms);
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    else {
      for (uint i=0; i<natoms; i++){
        if (IFC[i].size() != natoms){
          message = "The supplied IFC["+aurostd::utype2string(i);
          message += "] has the wrong size: ";
          message += aurostd::utype2string(IFC.size()) + " instead of ";
          message += aurostd::utype2string(natoms);
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
        }

        for (uint j=0; j<natoms; j++){
          if ((IFC[i][j].rows != 3) || (IFC[i][j].cols != 3)){
            message = "The supplied IFC["+aurostd::utype2string(i)+"][";
            message += aurostd::utype2string(j)+"] has the wrong size: ";
            message += aurostd::utype2string(IFC[i][j].rows)+"x"+aurostd::utype2string(IFC[i][j].cols);
            message += " instead of 3x3";
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
          }
        }
      }
    }

    _forceConstantMatrices = IFC;
    _isPolarMaterial = isPolar;
    uint niatoms = _supercell.getNumberOfUniqueAtoms();
    if (isPolar){
      if (bornEffectiveChargeTensor.size() != niatoms){
        message = "The supplied bornEffectiveChargeTensor has the wrong size: ";
        message += aurostd::utype2string(bornEffectiveChargeTensor.size()) + " instead of ";
        message += aurostd::utype2string(niatoms);
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
      else {
        for (uint i=0; i<niatoms; i++){
          if ((bornEffectiveChargeTensor[i].rows != 3) ||
              (bornEffectiveChargeTensor[i].cols != 3)){
            message = "The supplied bornEffectiveChargeTensor["+aurostd::utype2string(i);
            message += "] has the wrong size: ";
            message += aurostd::utype2string(bornEffectiveChargeTensor[i].rows);
            message += "x"+aurostd::utype2string(bornEffectiveChargeTensor[i].cols);
            message += " instead of 3x3";
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
          }
        }
      }

      if ((dielectricTensor.rows != 3) || (dielectricTensor.cols != 3)){
        message = "The supplied dielectricTensor has the wrong size: ";
        message += aurostd::utype2string(dielectricTensor.rows)+"x"+aurostd::utype2string(dielectricTensor.cols);
        message += " instead of 3x3";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }

      _bornEffectiveChargeTensor = bornEffectiveChargeTensor;
      _dielectricTensor = dielectricTensor;
      _inverseDielectricTensor = inverse(dielectricTensor);
      _recsqrtDielectricTensorDeterminant = 1.0/sqrt(determinant(_dielectricTensor));
    }
  }

  void PhononCalculator::setHarmonicForceConstants(const vector<vector<xmatrix<double> > > &IFC)
  {
    vector<xmatrix<double> > bornEffectiveChargeTensor;//dummy object, material is not polar
    xmatrix<double> dielectricTensor;//dummy object, material is not polar

    // set only IFC, material is not polar
    setHarmonicForceConstants(IFC, bornEffectiveChargeTensor, dielectricTensor, false);
  }
  // AS20201204 END

  void PhononCalculator::awake() {
    string base = _directory + "/" + DEFAULT_APL_FILE_PREFIX;
    readHarmonicIFCs(aurostd::CleanFileName(base + DEFAULT_APL_HARMIFC_FILE));
    if (_isPolarMaterial) readBornChargesDielectricTensor(aurostd::CleanFileName(base + DEFAULT_APL_POLAR_FILE));
  }

  void PhononCalculator::readHarmonicIFCs(const string& hibfile) {
    //CO, we already checked that it exists before, just open
    vector<string> vlines;                           //CO
    aurostd::efile2vectorstring(hibfile, vlines);  //CO //ME20181226
    string message = "";

    //CO START
    if (!vlines.size()) {
      message = "Cannot open output file " + hibfile + "."; //ME20181226
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }

    string line;
    uint line_count = 0;
    vector<string> tokens;

    // Test of xml...
    line = vlines[line_count++];
    //getline(infile, line);
    if (line.find("xml") == string::npos) {
      string message = "Not an xml file.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
    }
    //CO END

    // Check if a version string is in the xml file. If not, the force constants
    // follow an older, incompatible format and need to be recalculated
    while (true) {
      if (line_count == vlines.size()) {
        message = "The format for harmonic force constants has changed and is incomptable with the format found in this file.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
      }
      line = vlines[line_count++];  //CO
      if (line.find("aflow_version") != string::npos) break;
    }

    // Get force constant matrices
    while (true) {
      if (line_count == vlines.size()) { //CO
        message = "Cannot find <fcms> tag.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
      }
      line = vlines[line_count++];  //CO
      if (line.find("fcms") != string::npos)
        break;
    }
    //CO START
    line = vlines[line_count++];
    line = vlines[line_count++];
    //CO END
    vector<xmatrix<double> > row;
    xmatrix<double> m(3, 3);
    while (true) {
      if (line_count == vlines.size()) { //CO
        message = "Incomplete <fcms> tag.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
      }
      line = vlines[line_count++];  //CO
      if (line.find("</varray>") != string::npos) {
        _forceConstantMatrices.push_back(row);
        row.clear();
        line = vlines[line_count++];  //CO
        if (line.find("</varray>") != string::npos)
          break;
        line = vlines[line_count++];  //CO
      }

      for (int k = 1; k <= 3; k++) {
        line = vlines[line_count++];  //CO
        int t = line.find_first_of(">") + 1;
        aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        m(k, 1) = aurostd::string2utype<double>(tokens.at(0));
        m(k, 2) = aurostd::string2utype<double>(tokens.at(1));
        m(k, 3) = aurostd::string2utype<double>(tokens.at(2));
        tokens.clear();
      }
      row.push_back(m);
      line = vlines[line_count++];  //CO
    }
  }

  void PhononCalculator::readBornChargesDielectricTensor(const string& hibfile) {
    string message = "";
    if (!aurostd::EFileExist(hibfile)) {
      message = "Cannot find file " + hibfile + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }

    vector<string> vlines;
    aurostd::efile2vectorstring(hibfile, vlines);

    //CO START
    if (!vlines.size()) {
      message = "Cannot open output file " + hibfile + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }

    string line = "";
    uint line_count = 0;
    vector<string> tokens;

    // Test of xml...
    line = vlines[line_count++];
    if (line.find("xml") == string::npos) {
      message = "Not an xml file.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
    }

    // Get Born effective charge tensors
    xmatrix<double> m(3, 3);
    while (true) {
      if (line_count == vlines.size()) {
        message = "Cannot find <born> tag.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
      }
      line = vlines[line_count++];  //CO
      if (line.find("born") != string::npos) break;
    }

    line = vlines[line_count++];  //CO
    while (true) {
      if (line_count == vlines.size()) { //CO
        message = "Incomplete <born> tag.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
      }
      line = vlines[line_count++];  //CO
      if (line.find("</varray>") != string::npos) break;
      for (int k = 1; k <= 3; k++) {
        line = vlines[line_count++];  //CO
        int t = line.find_first_of(">") + 1;
        aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        m(k, 1) = aurostd::string2utype<double>(tokens.at(0));
        m(k, 2) = aurostd::string2utype<double>(tokens.at(1));
        m(k, 3) = aurostd::string2utype<double>(tokens.at(2));
        tokens.clear();
      }
      _bornEffectiveChargeTensor.push_back(m);
      line = vlines[line_count++];  //CO
    }

    // Get dielectric constant tensor
    while (true) {
      if (line_count == vlines.size()) {
        message = "Cannot find <epsilon> tag.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
      }
      line = vlines[line_count++];  //CO
      if (line.find("epsilon") != string::npos) break;
    }

    line = vlines[line_count++];  //CO
    for (int k = 1; k <= 3; k++) {
      line = vlines[line_count++];  //CO
      int t = line.find_first_of(">") + 1;
      aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
      _dielectricTensor(k, 1) = aurostd::string2utype<double>(tokens.at(0));
      _dielectricTensor(k, 2) = aurostd::string2utype<double>(tokens.at(1));
      _dielectricTensor(k, 3) = aurostd::string2utype<double>(tokens.at(2));
      tokens.clear();
    }
    _inverseDielectricTensor = inverse(_dielectricTensor);
    _recsqrtDielectricTensorDeterminant = 1.0 / sqrt(determinant(_dielectricTensor));
  }

  void PhononCalculator::setAnharmonicForceConstants(const AnharmonicIFCs& fc) {
    int order = fc.getOrder();
    int nifcs = (int) anharmonicIFCs.size();
    vector<vector<double> > dummy_ifc;
    vector<vector<int> > dummy_clst;
    for (int i = nifcs; i < order - 2; i++) {
      anharmonicIFCs.push_back(dummy_ifc);
      clusters.push_back(dummy_clst);
    }
    anharmonicIFCs[order - 3] = fc.getForceConstants();
    clusters[order - 3] = fc.getClusters();
  }

  void PhononCalculator::readAnharmonicIFCs(string filename) {
    filename = aurostd::CleanFileName(filename);
    string message = "";
    if (!aurostd::EFileExist(filename)) {
      message = "Could not open file " + filename + ". File not found.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }

    vector<string> vlines;
    aurostd::efile2vectorstring(filename, vlines);
    uint nlines = vlines.size();
    if (nlines == 0) {
      message = "Cannot open file " + filename + ". File empty or corrupt.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    }

    uint line_count = 0;
    string line = vlines[line_count++];

    // Check that this is a valid xml file
    if (line.find("xml") == string::npos) {
      message = "File is not a valid xml file.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
    }

    int t = 0;
    vector<string> tokens;

    // Read order
    uint order = 0;
    while (true) {
      if (line_count == nlines) {
        message = "order tag not found.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
      }
      line = vlines[line_count++];
      if (line.find("order") != string::npos) {
        t = line.find_first_of(">") + 1;
        order = aurostd::string2utype<uint>(line.substr(t, line.find_last_of("<") - t));
        break;
      }
    }

    // Read IFCs
    while (true) {
      if (line_count == nlines) {
        message = "force_constants tag not found.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
      }
      line = vlines[line_count++];
      if (line.find("force_constants") != string::npos) break;
    }

    vector<vector<double> > ifcs;
    vector<vector<int> > clst;
    uint nanharm = anharmonicIFCs.size();
    for (uint i = nanharm; i < order - 2; i++) {
      anharmonicIFCs.push_back(ifcs);
      clusters.push_back(clst);
    }
    vector<double> fc;
    vector<int> cl;
    while (line.find("/force_constants") == string::npos) {
      if (line_count == nlines) {
        message = "force_constants tag incomplete.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
      }
      line = vlines[line_count++];
      if (line.find("atoms") != string::npos) {
        // New tensor and cluster
        fc.clear();
        cl.clear();
        t = line.find_first_of("\"") + 1;
        aurostd::string2tokens(line.substr(t, line.find_last_of("\"") - t), cl, " ");
        clst.push_back(cl);
      } else if (line.find("slice") != string::npos) {
        // Populate tensor
        for (int i = 0; i < 3; i++) {
          tokens.clear();
          line = vlines[line_count++];
          t = line.find_first_of(">") + 1;
          aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, " ");
          for (int j = 0; j < 3; j++) fc.push_back(aurostd::string2utype<double>(tokens[j]));
        }
        line_count++;  // Skip /varray
      } else if (line.find("</varray>") != string::npos) {
        ifcs.push_back(fc);
      }
    }

    if (ifcs.size() != clst.size()) {
      message = "Number of IFCs is different from the number of clusters.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    }

    // Populate
    anharmonicIFCs[order - 3] = ifcs;
    clusters[order - 3] = clst;
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                           DYNAMICAL MATRIX                               //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  //ME20180827 - Overloaded to calculate derivative and eigenvectors for AAPL
  //ME20200206 - Added variants for the case near the Gamma point where the
  // non-analytical correction also needs a direction.
  xvector<double> PhononCalculator::getFrequency(const xvector<double>& kpoint, const IPCFreqFlags& flags) {
    return getFrequency(kpoint, kpoint, flags);
  }

  xvector<double> PhononCalculator::getFrequency(const xvector<double>& kpoint, const xvector<double>& kpoint_nac, const IPCFreqFlags& flags) {
    uint _nBranches = getNumberOfBranches();
    xmatrix<xcomplex<double> > placeholder_eigen(_nBranches, _nBranches, 1, 1);
    return getFrequency(kpoint, kpoint_nac, flags, placeholder_eigen);
  }

  //ME20190624 - get eigenvectors and frequencies
  xvector<double> PhononCalculator::getFrequency(const xvector<double>& kpoint, const IPCFreqFlags& flags,
      xmatrix<xcomplex<double> >& eigenvectors) {
    return getFrequency(kpoint, kpoint, flags, eigenvectors);
  }

  xvector<double> PhononCalculator::getFrequency(const xvector<double>& kpoint, const xvector<double>& kpoint_nac,
      const IPCFreqFlags& flags, xmatrix<xcomplex<double> >& eigenvectors) {
    vector<xvector<double> > placeholder_gvel;
    return getFrequency(kpoint, kpoint_nac, flags, eigenvectors, placeholder_gvel, false);
  }

  xvector<double> PhononCalculator::getFrequency(const xvector<double>& kpoint,
      const IPCFreqFlags& flags, xmatrix<xcomplex<double> >& eigenvectors,
      vector<xvector<double> >& gvel, bool calc_derivative) {
    return getFrequency(kpoint, kpoint, flags, eigenvectors, gvel, calc_derivative);
  }

  xvector<double> PhononCalculator::getFrequency(const xvector<double>& kpoint, const xvector<double>& kpoint_nac,
      const IPCFreqFlags& flags, xmatrix<xcomplex<double> >& eigenvectors,
      vector<xvector<double> >& gvel, bool calc_gvel) {
    // Compute frequency(omega) from eigenvalues [in eV/A/A/atomic_mass_unit]
    vector<xmatrix<xcomplex<double> > > dDynMat;
    xvector<double> omega = getEigenvalues(kpoint, kpoint_nac, eigenvectors, dDynMat, calc_gvel);

    // Get value of conversion factor
    double conversionFactor = getFrequencyConversionFactor(apl::RAW | apl::OMEGA, flags);

    // Transform values to desired format
    for (int i = omega.lrows; i <= omega.urows; i++) {
      if (omega(i) < 0) {
        if (flags & ALLOW_NEGATIVE)
          omega(i) = -sqrt(-omega(i));
        else
          omega(i) = 0.0;
      } else {
        omega(i) = sqrt(omega(i));
      }

      // Convert to desired units
      omega(i) *= conversionFactor;
    }

    if (calc_gvel) {
      double conversionFactorGvel = getFrequencyConversionFactor(flags, apl::THZ | apl::OMEGA);
      uint nbranches = getNumberOfBranches();
      gvel.clear();
      gvel.resize(nbranches, xvector<double>(3));
      xvector<xcomplex<double> > eigenvec(3), eigenvec_conj(3);
      xcomplex<double> prod;
      for (uint br = 0; br < nbranches; br++) {
        if (omega[br + 1] > _FLOAT_TOL_) {
          eigenvec = eigenvectors.getcol(br + 1);
          eigenvec_conj = conj(eigenvec);
          for (uint i = 1; i < 4; i++) {
            prod = eigenvec_conj * (dDynMat[i - 1] * eigenvec);
            gvel[br][i] = au2nmTHz * prod.re/(2.0 * omega[br + 1] * conversionFactorGvel);
          }
        }
      }
    }

    // Return
    return (omega);
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20200108 - replaced with constants in xscalar
  double PhononCalculator::getFrequencyConversionFactor(IPCFreqFlags inFlags, IPCFreqFlags outFlags) {
    double conversionFactor = 1.0;

    // Conversion from eV/A/A/atomic_mass_unit -> something
    if (inFlags & apl::RAW) {
      if (outFlags & apl::RAW) {
        // Transform to eV/A/A/atomic_mass_unit
        conversionFactor = 1.0;
      } else if (outFlags & apl::HERTZ) {
        // Transform to s-1; sqrt(EV_TO_JOULE / (ANGSTROM_TO_METER*ANGSTROM_TO_METER) / AMU_TO_KG);
        conversionFactor = au2Hz;
      } else if (outFlags & apl::THZ) {
        // Transform to THz; (in Hertz) / 1E12;
        conversionFactor = au2Hz * Hz2THz;
      } else if (outFlags & apl::RECIPROCAL_CM) {
        // Transform to cm-1; 1/lambda(m) = freq.(s-1) / light_speed(m/s);
        conversionFactor = au2rcm;
      } else if (outFlags & apl::MEV) {
        // Transform to meV; E(eV) = h(eV.s) * freq(s-1); h[(from J.s) -> (eV.s)] = 4.1356673310E-15
        conversionFactor = 1000 * au2eV;
      } else {
        string message = "Not implemented conversion.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
      }
    }

    // Conversion from THz -> something
    else if (inFlags & apl::THZ) {
      if (outFlags & apl::RAW) {
        // Transform to eV/A/A/atomic_mass_unit
        conversionFactor = 1.0 / (au2Hz * Hz2THz);
      } else if (outFlags & apl::THZ) {
        conversionFactor = 1.0;
      } else if (outFlags & apl::MEV) {
        conversionFactor = 1000 * PLANCKSCONSTANTEV_h * THz2Hz;
      } else {
        string message = "Not implemented conversion.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
      }
    }

    // Nothing suits?
    else {
      string message = "Not implemented conversion.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }

    //
    if ((outFlags & OMEGA) && !(inFlags & OMEGA))
      conversionFactor *= 2.0 * M_PI;
    if (!(outFlags & OMEGA) && (inFlags & OMEGA))
      conversionFactor /= 2.0 * M_PI;

    //
    return (conversionFactor);
  }

  // ///////////////////////////////////////////////////////////////////////////

  xvector<double> PhononCalculator::getEigenvalues(const xvector<double>& kpoint,
      const xvector<double>& kpoint_nac,
      xmatrix<xcomplex<double> >& eigenvectors,
      vector<xmatrix<xcomplex<double> > >& dDynMat,
      bool calc_derivative) {
    // Get dynamical matrix
    xmatrix<xcomplex<double> > dynamicalMatrix = getDynamicalMatrix(kpoint, kpoint_nac, dDynMat, calc_derivative);

    // Diagonalize
    xvector<double> eigenvalues(dynamicalMatrix.rows, 1);

    //ME20180828; OBSOLETE ME20190815 - use Jacobi algorithm in aurostd::xmatrix, which
    // is much, much faster than aplEigensystems for large systems
    //    apl::aplEigensystems e;
    //    e.eigen_calculation(dynamicalMatrix, eigenvalues, eigenvectors, APL_MV_EIGEN_SORT_VAL_ASC);

    eigenvalues = jacobiHermitian(dynamicalMatrix, eigenvectors);  //ME20190815

    return eigenvalues;
  }

  //  // ///////////////////////////////////////////////////////////////////////////
  //ME20180827 - Overloaded to calculate derivative for AAPL
  //ME20200206 - Added variants for the case near the Gamma point where the
  // non-analytical correction also needs a direction. While dynamical matrices
  // are not used directly, these functions are helpful debugging tools.
  xmatrix<xcomplex<double> > PhononCalculator::getDynamicalMatrix(const xvector<double>& kpoint) {
    return getDynamicalMatrix(kpoint, kpoint);
  }

  xmatrix<xcomplex<double> > PhononCalculator::getDynamicalMatrix(const xvector<double>& kpoint, const xvector<double>& kpoint_nac) {
    vector<xmatrix<xcomplex<double> > > placeholder;
    return getDynamicalMatrix(kpoint, kpoint_nac, placeholder, false);
  }

  xmatrix<xcomplex<double> > PhononCalculator::getDynamicalMatrix(const xvector<double>& kpoint,
      const xvector<double>& kpoint_nac,
      vector<xmatrix<xcomplex<double> > >& dDynMat,
      bool calc_derivative) {
    uint scAtomsSize = _supercell.getSupercellStructure().atoms.size();
    uint pcAtomsSize = _supercell.getInputStructure().atoms.size();

    uint _nBranches = getNumberOfBranches();
    xmatrix<xcomplex<double> > dynamicalMatrix(_nBranches, _nBranches, 1, 1);
    xmatrix<xcomplex<double> > dynamicalMatrix0(_nBranches, _nBranches, 1, 1);

    xcomplex<double> phase;
    double value = 0.0;
    //ME20180828 - Prepare derivative calculation
    xvector<xcomplex<double> > derivative(3);
    vector<xmatrix<xcomplex<double> > > dDynMat_NAC;
    if (calc_derivative) {  // reset dDynMat
      dDynMat.clear();
      xmatrix<xcomplex<double> > mat(_nBranches, _nBranches, 1, 1);
      dDynMat.assign(3, mat);
    }

    // Calculate nonanalytical contribution
    xmatrix<xcomplex<double> > dynamicalMatrixNA(_nBranches, _nBranches, 1, 1);
    if (_isPolarMaterial)
      dynamicalMatrixNA = getNonanalyticalTermWang(kpoint_nac, dDynMat_NAC, calc_derivative);

    // Loop over primitive cell
    xcomplex<double> nac;
    for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
      uint isc1 = _supercell.pc2scMap(ipc1);

      for (uint isc2 = 0; isc2 < scAtomsSize; isc2++) {
        uint ipc2 = _supercell.sc2pcMap(isc2);
        int neq = 0;  // Important for NAC derivative
        if (_supercell.calcShellPhaseFactor(isc2, isc1, kpoint, phase, neq, derivative, calc_derivative)) {  //ME20180827
          for (int ix = 1; ix <= 3; ix++) {
            for (int iy = 1; iy <= 3; iy++) {
              value = 0.5 * (_forceConstantMatrices[isc1][isc2](ix, iy) + _forceConstantMatrices[isc2][isc1](iy, ix));
              if (!aurostd::iszero(value)) {
                dynamicalMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) += value * phase;
                dynamicalMatrix0(3 * ipc1 + ix, 3 * ipc2 + iy) += value;
                if (calc_derivative) {
                  for (int d = 0; d < 3; d++) {
                    dDynMat[d](3 * ipc1 + ix, 3 * ipc2 + iy) += value * derivative[d+1];
                  }
                }
              }
              if (_isPolarMaterial) {
                dynamicalMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) += dynamicalMatrixNA(3 * ipc1 + ix, 3 * ipc2 + iy) * phase;
                dynamicalMatrix0(3 * ipc1 + ix, 3 * ipc2 + iy) += dynamicalMatrixNA(3 * ipc1 + ix, 3 * ipc2 + iy);
                if (calc_derivative && !aurostd::iszero(kpoint)) {
                  for (int d = 0; d < 3; d++) {
                    nac = ((double) neq) * phase * dDynMat_NAC[d](3 * ipc1 + ix, 3 * ipc2 + iy);
                    dDynMat[d](3 * ipc1 + ix, 3 * ipc2 + iy) += nac;
                  }
                }
              }
            }
          }
        }
      }
    }

    // Subtract the sum of all "forces" from the central atom, this is like an automatic sum rule...
    for (uint i = 0; i < pcAtomsSize; i++) {
      for (uint j = 0; j < pcAtomsSize; j++) {
        for (int ix = 1; ix <= 3; ix++) {
          for (int iy = 1; iy <= 3; iy++) {
            dynamicalMatrix(3 * i + ix, 3 * i + iy) = dynamicalMatrix(3 * i + ix, 3 * i + iy) - dynamicalMatrix0(3 * i + ix, 3 * j + iy);
          }
        }
      }
    }

    // Get correction for polar materials
    //if( _isPolarMaterial )
    // dynamicMatrix += getNonanalyticalTermGonze(kpoint);

    // Make it hermitian
    for (uint i = 0; i <= pcAtomsSize - 1; i++) {
      for (uint j = 0; j <= i; j++) {
        for (int ix = 1; ix <= 3; ix++) {
          for (int iy = 1; iy <= 3; iy++) {
            dynamicalMatrix(3 * i + ix, 3 * j + iy) += conj(dynamicalMatrix(3 * j + iy, 3 * i + ix));
            dynamicalMatrix(3 * i + ix, 3 * j + iy) *= 0.5;
            dynamicalMatrix(3 * j + iy, 3 * i + ix) = conj(dynamicalMatrix(3 * i + ix, 3 * j + iy));
            if (calc_derivative) {
              for (int d = 0; d < 3; d++) {
                dDynMat[d](3 * i + ix, 3 * j + iy) += conj(dDynMat[d](3 * j + iy, 3 * i + ix));
                dDynMat[d](3 * i + ix, 3 * j + iy) *= 0.5;
                dDynMat[d](3 * j + iy, 3 * i + ix) = conj(dDynMat[d](3 * i + ix, 3 * j + iy));
              }
            }
          }
        }
      }
    }

    // Divide by masses
    for (uint i = 0; i < pcAtomsSize; i++) {
      double mass_i = _supercell.getAtomMass(_supercell.pc2scMap(i));
      for (uint j = 0; j < pcAtomsSize; j++) {
        double mass_j = _supercell.getAtomMass(_supercell.pc2scMap(j));
        for (int ix = 1; ix <= 3; ix++) {
          for (int iy = 1; iy <= 3; iy++) {
            dynamicalMatrix(3 * i + ix, 3 * j + iy) *= 1.0 / sqrt(mass_i * mass_j);
            if (calc_derivative) {
              for (int d = 0; d < 3; d++) {
                dDynMat[d](3 * i + ix, 3 * j + iy) *= 1.0/sqrt(mass_i * mass_j);
              }
            }
          }
        }
      }
    }

    return dynamicalMatrix;
  }

  ///////////////////////////////////////////////////////////////////////////

  // Y. Wang et.al, J. Phys.:Condens. Matter 22, 202201 (2010)
  // DOI: 10.1088/0953-8984/22/20/202201

  //ME20180827 - Overloaded to calculate derivative for AAPL
  //ME20200207 - This function assummed that Born charges were stored for each type,
  // but it is actually stored for each iatom.
  xmatrix<xcomplex<double> > PhononCalculator::getNonanalyticalTermWang(const xvector<double>& _q) {
    vector<xmatrix<xcomplex<double> > > placeholder;
    return getNonanalyticalTermWang(_q, placeholder, false);
  }

  xmatrix<xcomplex<double> > PhononCalculator::getNonanalyticalTermWang(const xvector<double>& _q,
      vector<xmatrix<xcomplex<double> > >& derivative,
      bool calc_derivative) {
    const xstructure& sc = _supercell.getSupercellStructureLight();           //CO
    const xstructure& pc = _supercell.getInputStructure();  //CO  //ME20200207 - grab input structure (need iatoms)

    // to correct the q=\Gamma as a limit
    xvector<double> q(_q);
    if (aurostd::modulus(q) < _FLOAT_TOL_) {
      q(1) = _FLOAT_TOL_ * 1.001;
    }

    uint pcAtomsSize = pc.atoms.size();
    uint pcIAtomsSize = pc.iatoms.size();

    uint _nBranches = getNumberOfBranches();
    xmatrix<xcomplex<double> > dynamicalMatrix(_nBranches, _nBranches);

    if (aurostd::modulus(q) > _FLOAT_TOL_) {
      if (calc_derivative) {  // reset derivative
        derivative.clear();
        xmatrix<xcomplex<double> > mat(_nBranches, _nBranches, 1, 1);
        derivative.assign(3, mat);
      }

      // Calculation
      double fac0 = hartree2eV * bohr2angstrom;  // from a.u. to eV/A  //ME20200206 - replaced with xscalar constants
      double volume = det(pc.lattice);
      double fac1 = 4.0 * PI / volume;
      double nbCells = det(sc.lattice) / volume;

      // Precompute product of q-point with charge tensor
      vector<xvector<double> > qZ(pcIAtomsSize);
      for (uint at = 0; at < pcIAtomsSize; at++) qZ[at] = q * _bornEffectiveChargeTensor[at];

      double dotprod = scalar_product(q, _dielectricTensor * q);
      double prefactor = fac0 * fac1/(dotprod * nbCells);
      for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
        int iat1 = pc.atoms[ipc1].index_iatoms;
        for (uint ipc2 = 0; ipc2 < pcAtomsSize; ipc2++) {
          int iat2 = pc.atoms[ipc2].index_iatoms;
          for (int ix = 1; ix <= 3; ix++) {
            for (int iy = 1; iy <= 3; iy++) {
              //int typei = pc.atoms[ipc1].type;
              //int typej = pc.atoms[ipc2].type;
              //double borni = (q * _bornEffectiveChargeTensor[typei])(ix);
              //double bornj = (q * _bornEffectiveChargeTensor[typej])(iy);
              double borni = qZ[iat1][ix];
              double bornj = qZ[iat2][iy];
              dynamicalMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) = prefactor * borni * bornj;
              if (calc_derivative) {
                for (int d = 0; d < 3; d++) {
                  xcomplex<double> coeff(0, 0);
                  //coeff += borni * _bornEffectiveChargeTensor[ipc1](iy, d + 1);
                  //coeff += bornj * _bornEffectiveChargeTensor[ipc2](ix, d + 1);
                  coeff += borni * _bornEffectiveChargeTensor[iat2](iy, d + 1);
                  coeff += bornj * _bornEffectiveChargeTensor[iat1](ix, d + 1);
                  coeff -= 2 * borni * bornj * scalar_product(_dielectricTensor(d + 1), q)/dotprod;
                  derivative[d](3 * ipc1 + ix, 3 * ipc2 + iy) = prefactor * coeff;
                }
              }
            }
          }
        }
      }
    }

    return dynamicalMatrix;
  }

  // ///////////////////////////////////////////////////////////////////////////

  // X. Gonze et al., Phys. Rev. B 50, 13035 (1994)
  // X. Gonze and Ch. Lee, Phys. Rev. B 55, 10355 (1997)
  //ME20200504 - This function does not appear to be working!

  xmatrix<xcomplex<double> > PhononCalculator::getNonanalyticalTermGonze(const xvector<double> kpoint) {
    uint pcAtomsSize = _supercell.getInputStructure().atoms.size();

    if (!_isGammaEwaldPrecomputed) {
      xvector<double> zero(3);
      xmatrix<xcomplex<double> > dynamicalMatrix0(getEwaldSumDipoleDipoleContribution(zero, false));

      _gammaEwaldCorr.clear();
      for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
        xmatrix<xcomplex<double> > sum(3, 3);
        for (uint ipc2 = 0; ipc2 < pcAtomsSize; ipc2++) {
          for (int ix = 1; ix <= 3; ix++)
            for (int iy = 1; iy <= 3; iy++)
              sum(ix, iy) += dynamicalMatrix0(3 * ipc1 + ix, 3 * ipc2 + iy);
        }
        _gammaEwaldCorr.push_back(sum);
      }

      _isGammaEwaldPrecomputed = true;
    }

    //
    xmatrix<xcomplex<double> > dynamicalMatrix(getEwaldSumDipoleDipoleContribution(kpoint));

    for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
      for (int ix = 1; ix <= 3; ix++)
        for (int iy = 1; iy <= 3; iy++)
          dynamicalMatrix(3 * ipc1 + ix, 3 * ipc1 + iy) -= _gammaEwaldCorr[ipc1](ix, iy);
    }

    //
    return dynamicalMatrix;
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20200207 - This function assummed that Born charges were stored for each type,
  // but it is actually stored for each iatom.
  //ME20200504 - This function does not appear to be working!
  xmatrix<xcomplex<double> > PhononCalculator::getEwaldSumDipoleDipoleContribution(const xvector<double> qpoint, bool includeTerm1) {
    // Definitions
    const xstructure& sc = _supercell.getSupercellStructureLight();           //CO
    const xstructure& pc = _supercell.getInputStructure();  //CO  //ME20200207 - grab input structure (need iatoms)

    uint pcAtomsSize = pc.atoms.size();

    uint _nBranches = getNumberOfBranches();
    xmatrix<xcomplex<double> > dynamicalMatrix(_nBranches, _nBranches);

    double gmax = 14.0;
    double lambda = 1.0;
    double lambda2 = lambda * lambda;
    double lambda3 = lambda2 * lambda;
    double geg = gmax * lambda2 * 4.0;

    // Reciprocal Space
    xmatrix<double> klattice = trasp(ReciprocalLattice(pc.lattice));

    // Grid
    int n1 = (int)(sqrt(geg) / aurostd::modulus(klattice(1))) + 1;
    int n2 = (int)(sqrt(geg) / aurostd::modulus(klattice(2))) + 1;
    int n3 = (int)(sqrt(geg) / aurostd::modulus(klattice(3))) + 1;

    // Calculation
    double fac0 = hartree2eV * bohr2angstrom;  // from a.u. to eV/A  //ME20200207 - replaced with xscalar constants
    double SQRTPI = sqrt(PI);
    double volume = det(pc.lattice);
    double fac = 4.0 * PI / volume;
    xcomplex<double> iONE(0.0, 1.0);

    // Term 1 - Reciprocal space sum

    if (includeTerm1) {
      for (int m1 = -n1; m1 <= n1; m1++) {
        for (int m2 = -n2; m2 <= n2; m2++) {
          for (int m3 = -n3; m3 <= n3; m3++) {
            xvector<double> g = m1 * klattice(1) + m2 * klattice(2) + m3 * klattice(3) + qpoint;

            geg = scalar_product(g, _dielectricTensor * g);

            if (aurostd::abs(geg) > _FLOAT_TOL_ && geg / lambda2 / 4.0 < gmax) {
              double fac2 = fac * exp(-geg / lambda2 / 4.0) / geg;

              for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
                //xvector<double> zag = g * _bornEffectiveChargeTensor[pc.atoms[ipc1].type];
                int iat1 = pc.atoms[ipc1].index_iatoms;
                xvector<double> zag = g * _bornEffectiveChargeTensor[iat1];

                for (uint ipc2 = 0; ipc2 < pcAtomsSize; ipc2++) {
                  //xvector<double> zbg = g * _bornEffectiveChargeTensor[pc.atoms[ipc2].type];
                  int iat2 = pc.atoms[ipc2].index_iatoms;
                  xvector<double> zbg = g * _bornEffectiveChargeTensor[iat2];

                  //xcomplex<double> e;
                  //(void)_supercell.calcShellPhaseFactor(ipc2,ipc1,g,e);
                  //xcomplex<double> facg = fac2 * e;
                  xcomplex<double> facg = fac2 * exp(iONE * scalar_product(g, sc.atoms[ipc2].cpos - sc.atoms[ipc1].cpos));

                  for (int ix = 1; ix <= 3; ix++) {
                    for (int iy = 1; iy <= 3; iy++) {
                      dynamicalMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) += fac0 * facg * zag(ix) * zbg(iy);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    // Term 2 - Real space sum
    //for(int m1 = -n1; m1 <= n1; m1++)
    //  for(int m2 = -n2; m2 <= n2; m2++)
    //    for(int m3 = -n2; m3 <= n3; m3++) {
    //      xvector<double> rc = m1 * pc.lattice(1) + m2 * pc.lattice(2)
    //        + m3 * pc.lattice(3);

    //      //xvector<double> zero(3);
    //      //xvector<double> rf = _supercell.getFPositionItsNearestImage(rc,zero,pc.lattice);
    //      //rc = F2C(pc.lattice,rf);

    //      if( aurostd::modulus(rc) < _FLOAT_TOL_ ) continue;

    //      //
    //      xvector<double> delta = _inverseDielectricTensor * rc;
    //      double D = sqrt( scalar_product(delta,rc) );

    //      //
    //      xmatrix<double> H(3,3);
    //      xvector<double> x = lambda * delta;
    //      double y = lambda * D;
    //      double y2 = y * y;
    //      double ym2 = 1.0 / y2;
    //      double emy2dpi = 2.0 * exp( -y2 ) / SQRTPI;
    //      double erfcdy = erfc(y) / y;
    //      double c1 = ym2 * ( 3.0 * erfcdy * ym2 + ( emy2dpi * ( 3.0 * ym2 + 2.0 ) ) );
    //      double c2 = ym2 * ( erfcdy + emy2dpi );
    //      for(int a = 1; a <= 3; a++)
    //        for(int b = 1; b <= 3; b++) {
    //          H(a,b) = x(a) * x(b) * c1 - _inverseDielectricTensor(a,b) * c2;
    //        }

    //      //
    //      xcomplex<double> e = exp( iONE * scalar_product(qpoint,rc) );
    //      xcomplex<double> fac = fac0 * lambda3 * _recsqrtDielectricTensorDeterminant * e;

    //      //
    //      for(uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
    //        xmatrix<double> zh = _bornEffectiveChargeTensor[pc.atoms[ipc1].type] * H;

    //        for(uint ipc2 = 0; ipc2 < pcAtomsSize; ipc2++) {
    //          xmatrix<double> zhz = zh * _bornEffectiveChargeTensor[pc.atoms[ipc2].type];

    //          for(int ix = 1; ix <= 3; ix++)
    //            for(int iy = 1; iy <= 3; iy++)
    //              dynamicalMatrix(3*ipc1+ix,3*ipc2+iy) -= fac * zhz(ix,iy);
    //        }
    //      }
    //    }

    // Term 2
    uint scAtomsSize = sc.atoms.size();
    for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
      uint isc1 = _supercell.pc2scMap(ipc1);

      for (uint isc2 = 0; isc2 < scAtomsSize; isc2++) {
        uint ipc2 = _supercell.sc2pcMap(isc2);

        xvector<double> rc = SYM::minimizeDistanceCartesianMethod(sc.atoms[isc2].cpos, sc.atoms[isc1].cpos, sc.lattice);
        xvector<double> rf = F2C(sc.lattice, rc);

        if (aurostd::modulus(rc) < _FLOAT_TOL_) continue;

        //
        xvector<double> delta = _inverseDielectricTensor * rc;
        double D = sqrt(scalar_product(delta, rc));

        //
        xmatrix<double> H(3, 3);
        xvector<double> x = lambda * delta;
        double y = lambda * D;
        double y2 = y * y;
        double ym2 = 1.0 / y2;
        double emy2dpi = 2.0 * exp(-y2) / SQRTPI;
        double erfcdy = erfc(y) / y;
        double c1 = ym2 * (3.0 * erfcdy * ym2 + (emy2dpi * (3.0 * ym2 + 2.0)));
        double c2 = ym2 * (erfcdy + emy2dpi);
        for (int a = 1; a <= 3; a++) {
          for (int b = 1; b <= 3; b++) {
            H(a, b) = x(a) * x(b) * c1 - _inverseDielectricTensor(a, b) * c2;
          }
        }

        //xmatrix<double> za = _bornEffectiveChargeTensor[pc.atoms[ipc1].type];
        //xmatrix<double> zb = _bornEffectiveChargeTensor[pc.atoms[ipc2].type];
        int iat1 = pc.atoms[ipc1].index_iatoms;
        int iat2 = pc.atoms[ipc2].index_iatoms;
        xmatrix<double> za = _bornEffectiveChargeTensor[iat1];
        xmatrix<double> zb = _bornEffectiveChargeTensor[iat2];
        xmatrix<double> zhz = za * H * zb;

        //
        xcomplex<double> e;  // = exp( iONE * scalar_product(qpoint,rc) );
        (void)_supercell.calcShellPhaseFactor(isc2, isc1, qpoint, e);

        //
        xcomplex<double> fac = fac0 * lambda3 * _recsqrtDielectricTensorDeterminant * e;
        for (int ix = 1; ix <= 3; ix++)
          for (int iy = 1; iy <= 3; iy++)
            dynamicalMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) -= fac * zhz(ix, iy);
      }
    }

    // Term 3 - Limiting contribution

    double facterm3 = fac0 * 4.0 * lambda3 * _recsqrtDielectricTensorDeterminant / (3.0 * SQRTPI);
    for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
      //xmatrix<double> z = _bornEffectiveChargeTensor[pc.atoms[ipc1].type];
      int iat1 = pc.atoms[ipc1].index_iatoms;
      xmatrix<double> z = _bornEffectiveChargeTensor[iat1];
      xmatrix<double> zez = z * _inverseDielectricTensor * z;

      for (int ix = 1; ix <= 3; ix++)
        for (int iy = 1; iy <= 3; iy++)
          dynamicalMatrix(3 * ipc1 + ix, 3 * ipc1 + iy) -= facterm3 * zez(ix, iy);
    }

    //
    return dynamicalMatrix;
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                           GROUP VELOCITIES                               //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  // Calculate the group velocities on a q-point mesh using the Hellmann-Feynman
  // theorem. Group velocities will be in km/s (nm THz).
  vector<vector<xvector<double> > > PhononCalculator::calculateGroupVelocitiesOnMesh() {
    vector<vector<double> > freqs_placeholder;
    vector<xmatrix<xcomplex<double> > > eigenvectors_placeholder;
    return calculateGroupVelocitiesOnMesh(freqs_placeholder, eigenvectors_placeholder);
  }

  vector<vector<xvector<double> > > PhononCalculator::calculateGroupVelocitiesOnMesh(vector<vector<double> >& freqs) {
    vector<xmatrix<xcomplex<double> > > eigenvectors_placeholder;
    return calculateGroupVelocitiesOnMesh(freqs, eigenvectors_placeholder);
  }

  vector<vector<xvector<double> > > PhononCalculator::calculateGroupVelocitiesOnMesh(vector<vector<double> >& freqs, vector<xmatrix<xcomplex<double> > >& eigenvectors) {
    string message = "";
    if (!_supercell.isConstructed()) {
      message = "Supercell not constructed yet.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    int nQPs = (int) _qm.getnQPs();
    if (nQPs == 0) {
      message = "Mesh has no q-points.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }

    freqs.clear();
    eigenvectors.clear();

    uint nbranches = getNumberOfBranches();
    freqs.resize(nQPs);
    eigenvectors.resize(nQPs, xmatrix<xcomplex<double> >(nbranches, nbranches));
    vector<vector<xvector<double> > > gvel(nQPs);
#ifdef AFLOW_MULTITHREADS_ENABLE
    xthread::xThread xt(_ncpus);
    std::function<void(int, vector<vector<double> >&,
        vector<xmatrix<xcomplex<double> > >&,
        vector<vector<xvector<double> > >&)> fn = std::bind(&PhononCalculator::calculateGroupVelocitiesThread, this, _1, _2, _3, _4);
    xt.run(nQPs, fn, freqs, eigenvectors, gvel);
#else
    for (int i = 0; i < nQPs; i++) calculateGroupVelocitiesThread(i, freqs, eigenvectors, gvel);
#endif
    return gvel;
  }

  void PhononCalculator::calculateGroupVelocitiesThread(int q,
      vector<vector<double> >& freqs, vector<xmatrix<xcomplex<double> > >& eigenvectors, vector<vector<xvector<double> > >& gvel) {
    xvector<double> f;
    f = getFrequency(_qm.getQPoint(q).cpos, apl::THZ | apl::OMEGA, eigenvectors[q], gvel[q]);
    freqs[q] = aurostd::xvector2vector(f);
  }

  //writeGroupVelocitiesToFile////////////////////////////////////////////////
  // Writes the group velocities into a file. Each row belongs to a q-point,
  // and each column triplet belongs to a phonon branch.
  void PhononCalculator::writeGroupVelocitiesToFile(const string& filename,
      const vector<vector<xvector<double> > >& gvel) {
    vector<vector<double> > freqs;
    writeGroupVelocitiesToFile(filename, gvel, freqs);
  }
  void PhononCalculator::writeGroupVelocitiesToFile(const string& filename,
      const vector<vector<xvector<double> > >& gvel, const vector<vector<double> >& freqs, const string& unit) {
    if (gvel.size() == 0) return;  // Nothing to write
    string message = "";
    stringstream output;

    uint nBranches = getNumberOfBranches();
    uint nQPs = _qm.getnQPs();

    // Consistency check
    if (gvel.size() != nQPs) {
      message = "Number of group velocities is not equal to the number of q-points.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    for (uint q = 0; q < nQPs; q++) {
      if (gvel[q].size() != nBranches) {
        message = "Number of group velocities for q-point " + aurostd::utype2string<uint>(q) + " is not equal to the number branches.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
    }

    if (freqs.size() > 0) {
      if (freqs.size() != nQPs) {
        message = "Number of frequencies is not equal to the number of q-points.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
      for (uint q = 0; q < nQPs; q++) {
        if (freqs[q].size() != nBranches) {
          message = "Number of frequencies for q-point " + aurostd::utype2string<uint>(q) + " is not equal to the number branches.";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
        }
      }
    }

    // Header
    output << AFLOWIN_SEPARATION_LINE << std::endl;
    if (!_system.empty()) {
      output << "[APL_GROUP_VELOCITY]SYSTEM=" << _system << std::endl;
      output << AFLOWIN_SEPARATION_LINE << std::endl;
    }

    output << "[APL_GROUP_VELOCITY]START" << std::endl;
    output << std::setiosflags(std::ios::fixed | std::ios::right);
    output << std::setw(10) << "# Q-point";
    output << std::setw(20) << " ";
    output << "Group Velocity (km/s)" << std::endl;

    // Body
    for (uint q = 0; q < nQPs; q++) {
      output << std::setiosflags(std::ios::fixed | std::ios::right);
      output << std::setw(10) << q;
      for (uint br = 0; br < nBranches; br++) {
        for (uint i = 1; i < 4; i++) {
          output << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
          output << std::setw(20) << std::setprecision(10) << std::scientific << gvel[q][br][i];
        }
        output << std::setw(5) << " ";
      }
      output << std::endl;
    }

    output << "[APL_GROUP_VELOCITY]STOP" << std::endl;
    output << AFLOWIN_SEPARATION_LINE << std::endl;

    // Write frequencies if provided
    if (freqs.size() > 0) {
      output << "[APL_FREQUENCY]STOP" << std::endl;
      output << std::setiosflags(std::ios::fixed | std::ios::right);
      output << std::setw(10) << "# Q-point"
        << std::setw(20) << " "
        << "Frequency " << (unit.empty()?"":("(" + unit + ")")) << std::endl;
      for (uint q = 0; q < nQPs; q++) {
        output << std::setiosflags(std::ios::fixed | std::ios::right);
        output << std::setw(10) << q;
        for (uint br = 0; br < nBranches; br++) {
          output << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
          output << std::setw(20) << std::setprecision(10) << std::scientific << freqs[q][br];
        }
        output << std::endl;
      }
      output << "[APL_FREQUENCY]STOP" << std::endl;
      output << AFLOWIN_SEPARATION_LINE << std::endl;
    }

    // Write q-points
    output << "[APL_QPOINTS]START" << std::endl;
    output << std::setiosflags(std::ios::fixed | std::ios::right);
    output << std::setw(10) << "# Index";
    output << std::setw(20) << " ";
    output << "Q-points (fractional)" << std::endl;
    // Body
    for (uint q = 0; q < nQPs; q++) {
      output << std::setiosflags(std::ios::fixed | std::ios::right);
      output << std::setw(10) << q;
      for (uint i = 1; i < 4; i++) {
        output << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
        output << std::setw(20) << std::setprecision(10) << std::scientific << _qm.getQPoint(q).fpos[i];
      }
      output << std::endl;
    }
    output << "[APL_QPOINTS]STOP" << std::endl;
    output << AFLOWIN_SEPARATION_LINE << std::endl;

    // Write to file
    aurostd::stringstream2file(output, filename);
    if (!aurostd::FileExist(filename)) {
      message = "Could not write group velocities to file.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
  }

}  // namespace apl

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
