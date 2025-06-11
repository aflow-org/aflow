// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// MAKEFILE FOR AFLOW_APL
// Written by Michal Jahnatek

#ifndef _AFLOW_APL_H_
#define _AFLOW_APL_H_

// ME20200516 - Tolerance for what is consindered to be a different coordination shell
#define _APL_SHELL_TOL_ 0.1

// Aflow core libraries
#include "../aflow.h"
#include "../AUROSTD/aurostd.h"

//ME20190219 - Define the checksum algorithm used for APL hibernate files
#define APL_CHECKSUM_ALGO  string("Fletcher32")

// ***************************************************************************

// Interface
namespace apl {

  void validateParametersAPL(xoption&, const _aflags&, ofstream&, ostream& oss=std::cout);
  void validateParametersSupercellAPL(xoption&);
  void validateParametersDispersionsAPL(xoption&);
  void validateParametersDosAPL(xoption&, const _aflags&, ofstream&, ostream& oss=std::cout);
  void validateParametersAAPL(xoption&, const _aflags&, ofstream&, ostream& oss=std::cout);
  void validateParametersQHA(xoption& aaplopts, const _aflags& aflags, ofstream& FileMESSAGE, ostream& oss); //AS20200709
  bool createAflowInPhonons(const _aflags&, const _kflags&, const _xflags&, _xinput&); // ME20190108
  void createAflowInPhononsAIMS(_aflags&, _kflags&, _xflags&, string&, _xinput&, ofstream&);
  bool filesExistPhonons(_xinput&);
  bool outfileFoundAnywherePhonons(vector<_xinput>&);
  bool outfileFoundEverywherePhonons(vector<_xinput>&, const string&, ofstream&, ostream&, bool=false);  // ME20191029
  bool readForcesFromDirectory(_xinput&);  // ME20200219
  void subtractZeroStateForces(vector<_xinput>&, bool);
  void subtractZeroStateForces(vector<_xinput>&, _xinput&);  // ME20190114

  bool APL_Get_AflowInName(string& AflowInName, const string& directory_LIB); //ME20210927
}

// ***************************************************************************
// BEGIN: Supplemental classes for APL, AAPL, and QHA
// ***************************************************************************

namespace apl {

  class Supercell : public xStream {
    private:
      xstructure _inStructure;
      xstructure _inStructure_original;  //CO
      xstructure _inStructure_light;     //CO, does not include HEAVY symmetry stuff
      //deque<_atom> _inStructure_atoms_original; //CO
      xstructure _pcStructure;           //CO20180406 - for the path
      xstructure _scStructure;
      xstructure _scStructure_original;  //CO
      xstructure _scStructure_light;     //CO, does not include HEAVY symmetry stuff
      //deque<_atom> _scStructure_atoms_original; //CO
      //CO START
      bool _skew;                  //SYM::isLatticeSkewed(), same for pc and sc
      bool _derivative_structure;  //vs. simple expanded_lattice, derivative structure's lattice has LESS symmetry, so be careful ApplyAtom()'ing
      double _sym_eps;             //same for pc and sc
      //CO END
      bool _isShellRestricted;
      int _maxShellID;
      vector<double> _maxShellRadius;
      bool _isConstructed;
      bool _initialized;
      vector<vector<vector<xvector<double> > > > phase_vectors;  // ME20200116

      void calculateWholeSymmetry(xstructure&, bool=true);
      xstructure calculatePrimitiveStructure() const;
      bool getMaps(const xstructure&, const xstructure&, const xstructure&, vector<int>&, vector<int>&);  // ME20200117
      void free();
      void copy(const Supercell&);

    public:
      Supercell(ostream& oss=std::cout);
      Supercell(ofstream&, ostream& os=std::cout);
      Supercell(const xstructure&, ofstream&, const string& directory="./", ostream& oss=std::cout); //CO20181226
      Supercell(const string&, ofstream&, const string& directory="./", ostream& os=std::cout);  // ME20200112
      Supercell(const Supercell&);
      Supercell& operator=(const Supercell&);
      ~Supercell();
      void clear();
      void initialize(const string&, ofstream&, ostream& oss=std::cout);
      void initialize(const string&);  // ME20200212
      void initialize(const xstructure&, ofstream&, bool=true, ostream& oss=std::cout);
      void initialize(const xstructure&, bool=true);  // ME20191225
      void initialize(ostream& oss);
      void initialize(ofstream& mf, ostream& oss);
      void clearSupercell();
      //void LightCopy(const xstructure& a, xstructure& b);  // OBSOLETE ME20200220
      bool isConstructed();
      void reset();
      xvector<int> determineSupercellDimensions(const aurostd::xoption&);  // ME20191225
      void build(aurostd::xoption&, bool = true);  // ME20191225
      void build(const xvector<int>&, bool = true);  // ME20191225
      void build(int, int, int, bool = TRUE);
      void trimStructure(int, const xvector<double>&,
          const xvector<double>&, const xvector<double>&,
          bool = true);
      bool projectToPrimitive();  // ME20200117
      void projectToOriginal();  // ME20200117
      xvector<int> getSupercellDimensionsShell(uint, bool);
      void setupShellRestrictions(int);
      //ME20190715 BEGIN - added const to getter functions so they can be used with const Supercell &
      bool isShellRestricted() const;
      int getMaxShellID() const;
      uint getNumberOfAtoms() const;
      uint getNumberOfUniqueAtoms() const;
      uint getNumberOfEquivalentAtomsOfType(int) const; //CO20190218
      int getUniqueAtomID(int) const;
      int getUniqueAtomID(int, int) const;
      const _atom& getUniqueAtom(int) const;
      string getUniqueAtomSymbol(int) const;
      double getUniqueAtomMass(int) const;
      double getAtomMass(int) const;
      int getAtomNumber(int) const;
      //ME20190715 END
      const xstructure& getSupercellStructure() const;
      const xstructure& getSupercellStructureLight() const;
      const xstructure& getPrimitiveStructure() const;
      const xstructure& getInputStructure() const;
      const xstructure& getInputStructureLight() const;
      const xstructure& getOriginalStructure() const;  // ME20200117
      int atomGoesTo(const _sym_op&, int, int, bool = TRUE); //CO20190218
      int atomComesFrom(const _sym_op&, int, int, bool = TRUE); //CO20190218
      const _sym_op& getSymOpWhichMatchAtoms(int, int, int);
      void calculatePhaseVectors();  // ME20200117
      bool calcShellPhaseFactor(int, int, const xvector<double>&, xcomplex<double>&);
      bool calcShellPhaseFactor(int, int, const xvector<double>&, xcomplex<double>&,
          int&, xvector<xcomplex<double> >&, bool);  //ME20180828
      int pc2scMap(int) const;
      int sc2pcMap(int) const;
      void center(int);
      //CO START
      void center_original(void);
      //corey
      void getFullBasisAGROUP();  //ME20191218
      bool fullBasisCalculatedAGROUP();  //ME20191218
      const vector<vector<_sym_op> >& getAGROUP(void) const;
      const vector<_sym_op>& getFGROUP(void) const;
      const vector<_sym_op>& getAGROUP(int) const;
      bool isDerivativeStructure() const;
      double getEPS(void) const;
      //ME20190715 END
      //CO END
      // **** BEGIN JJPR *****
      xvector<int> scell_dim;
      vector<int> _pc2scMap;
      vector<int> _sc2pcMap;
      // **** END  JJPR *****
      string _directory;  // for the logger
  };

}  // namespace apl

// ***************************************************************************

namespace apl {

  struct _qpoint {
    xvector<double> cpos;  // Cartesian position of the q-point
    xvector<double> fpos;  // Fractional coordinates of the q-point
    int irredQpt;  // The irreducible q-point this q-point belongs to
    int ibzqpt;  // The index of the irreducible q-point in the _ibzqpt vector
    int symop;  // Symmetry operation to transform the irreducible q-point into this q-point
  };

  struct _kcell {
    xmatrix<double> lattice;  // The reciprocal lattice vectors
    xmatrix<double> rlattice;  // The real space lattice
    xmatrix<double> c2f;  // Conversion matrix from Cartesian to fractional
    xmatrix<double> f2c;  // Conversion matrix from fractional to Cartesian
    bool skewed;  // Is the lattice skewed?
    vector<_sym_op> pgroup;  // The point group operations of the reciprocal cell
  };

  class QMesh : public xStream {
    public:
      QMesh(ostream& oss=std::cout);
      QMesh(ofstream&, ostream& os=std::cout);
      QMesh(const xvector<int>&, const xstructure&, ofstream&, bool include_inversions=true, bool gamma_centered=true, const string& directory="./", ostream& oss=std::cout);
      QMesh(const vector<int>&, const xstructure&, ofstream&, bool include_inversions=true, bool gamma_centered=true, const string& directory="./", ostream& oss=std::cout);
      QMesh(const QMesh&);
      QMesh& operator=(const QMesh&);
      ~QMesh();

      void clear();
      void clear_tetrahedra();

      void initialize(const vector<int>&, const xstructure& xs, ofstream&, bool=true, bool=true, ostream& oss=std::cout);
      void initialize(const vector<int>&, const xstructure& xs, bool=true, bool=true);
      void initialize(const xvector<int>&, const xstructure& xs, ofstream&, bool=true, bool=true, ostream& oss=std::cout);
      void initialize(const xvector<int>&, const xstructure& xs, bool=true, bool=true);
      void initialize(ostream& oss);
      void initialize(ofstream& mf, ostream& oss);

      void makeIrreducible();
      void calculateLittleGroups();  // ME20200109
      void writeQpoints(string, bool=true);
      void writeIrredQpoints(string, bool=true);

      string _directory;

      int getnIQPs() const;
      int getnQPs() const;
      int getGrid(int) const;
      const xvector<int>& getGrid() const;
      const _qpoint& getIrredQPoint(int) const;
      const _qpoint& getIrredQPoint(int, int, int) const;
      vector<xvector<double> > getIrredQPointsCPOS() const;
      vector<xvector<double> > getIrredQPointsFPOS() const;
      int getIrredQPointIndex(int) const;
      int getIrredQPointIndex(int, int, int) const;
      const _qpoint& getQPoint(int) const;
      const _qpoint& getQPoint(int, int, int) const;
      const _qpoint& getQPoint(const xvector<double>&) const;  //ME20190813
      int getQPointIndex(xvector<double>) const;  //ME20190813
      int getQPointIndex(int, int, int) const;
      vector<xvector<double> > getQPointsCPOS() const;
      vector<xvector<double> > getQPointsFPOS() const;
      int getIbzqpt(int) const;
      int getIbzqpt(int, int, int) const;
      const vector<int>& getIbzqpts() const;
      const vector<_qpoint>& getPoints() const;
      const _kcell& getReciprocalCell() const;
      bool isShifted() const;  //ME20190813
      const xvector<double>& getShift() const;
      const vector<int>& getWeights() const;
      bool initialized() const;
      bool isReduced() const;
      bool isGammaCentered() const;
      bool littleGroupsCalculated() const;  // ME20200109
      const vector<int>& getLittleGroup(int) const;  // ME20200109

      // Tetrahedron method
      void generateTetrahedra();
      void makeIrreducibleTetrahedra();

      const vector<vector<int> >& getTetrahedra() const;
      const vector<int>& getTetrahedron(int) const;
      const vector<int>& getIrredTetrahedron(int) const;
      int getTetrahedronCorner(int, int) const;
      vector<vector<int> > getIrreducibleTetrahedra() const;
      vector<vector<int> > getIrreducibleTetrahedraIbzqpt() const;
      int getnTetrahedra() const;
      int getnIrredTetrahedra() const;
      double getVolumePerTetrahedron() const;
      const vector<int>& getWeightsTetrahedra() const;
      int getWeightTetrahedron(int) const;
      bool isReducedTetrahedra() const;

    private:
      void free();
      void copy(const QMesh&);

      vector<int> _ibzqpts;  // The indices of the irreducible q-points
      bool _initialized;  // Indicates whether the QMesh object has been intialized
      bool _isGammaCentered;  // Indicates whether the includes the Gamma point
      vector<vector<int> > _littleGroups;  // The little groups of the irreducible q-points
      bool _littleGroupsCalculated;  // Indicates whether the little groups have been calculated
      int _nIQPs;  // The number of irreducible q-points
      int _nQPs;  // The number of q-points
      xvector<int> _qptGrid;  // The dimensions of the q-point mesh
      vector<vector<vector<int> > > _qptMap;  // Maps a q-point triplet to a q-point index
      vector<_qpoint> _qpoints;  // The q-points of the mesh
      _kcell _recCell;  // The reciprocal cell
      bool _reduced;  // Indicates whether the q-point mesh has been reduced
      bool _shifted;  // Indicates whether the q-point mesh has been shifted
      xvector<double> _shift;  // The shift vector of the mesh
      vector<int> _weights;  // The weights of each irreducible q-point

      void setGrid(const xvector<int>&);
      void setupReciprocalCell(xstructure, bool);
      void generateGridPoints(bool);
      void shiftMesh(const xvector<double>&);
      void moveToBZ(xvector<double>&) const;

      // Tetrahedron method
      vector<vector<int> > _tetrahedra;  // The corners of the tetrahedra
      vector<int> _irredTetrahedra;  // List of irreducible tetrahedra
      bool _reducedTetrahedra; // Indicates whether the tetrahedra are reduced
      int _nTetra;  // The number of tetrahedra
      int _nIrredTetra;  // The number of irreducible tetrahedra - ME20190625
      double _volumePerTetrahedron;  // The relative volume of each tetrahedron
      vector<int> _weightsTetrahedra;  // The weights of each irreducible tetrahedron

      vector<vector<xvector<int> > > initializeTetrahedra();
      void findMostCompactTetrahedra(vector<vector<xvector<int> > >&);
      void generateAllTetrahedra(const vector<vector<xvector<int> > >&);
  };

}  // namespace apl

// ***************************************************************************
// END: Supplemental classes for APL, AAPL, and QHA
// ***************************************************************************

// ***************************************************************************
// BEGIN: Automatic Phonon Library (APL)
// ***************************************************************************

// ***************************************************************************

#define _AFLOW_APL_BORN_EPSILON_RUNNAME_ string("LRBE")  // ME20190108
#define _AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_ string(ARUN_DIRECTORY_PREFIX + "APL_" + _AFLOW_APL_BORN_EPSILON_RUNNAME_) // ME20190108
//#define _AFLOW_APL_FORCEFIELDS_RUNNAME_ string("LRFF")  // ME20190108  // OBSOLETE ME20200213 - the calculation does not use force fields
//#define _AFLOW_APL_FORCEFIELDS_DIRECTORY_NAME_ string(ARUN_DIRECTORY_PREFIX + "APL_" + _AFLOW_APL_FORCEFIELDS_RUNNAME_) // ME20190108  // OBSOLETE ME20200213
#define _AFLOW_APL_DFPT_RUNNAME_ string("DFPT")  // ME20200213
#define _AFLOW_APL_DFPT_DIRECTORY_NAME_ string(ARUN_DIRECTORY_PREFIX + "APL_" + _AFLOW_APL_DFPT_RUNNAME_) // ME20200213

namespace apl {

  class ForceConstantCalculator : public xStream {
    protected:
      Supercell* _supercell;
    private:
      bool _sc_set;
      bool _initialized;

      vector<_xinput> xInputs;
      string _method;

      // Calculate forces at no distortion - since for some structure
      // (not well relaxed, or with other problems) these forces have to be
      // known and substracted from the calculated forces with distortion
      bool _calculateZeroStateForces;

      // For each atom of supercell, there is a full force field
      vector<vector<xmatrix<double> > > _forceConstantMatrices;
      // Stuff for polar materials
      bool _isPolarMaterial;
      // For each atom there is a matrix 3x3 of Born effective charge
      vector<xmatrix<double> > _bornEffectiveChargeTensor;
      // Dielectric tensor
      xmatrix<double> _dielectricTensor;

      void free();
      void copy(const ForceConstantCalculator&);

      // Force constants
      bool runVASPCalculationsDM(_xinput&, _aflags&, _kflags&, _xflags&, string&);
      bool runVASPCalculationsLR(_xinput&, _aflags&, _kflags&, _xflags&, string&);
      bool calculateForceConstants(); // ME20200211
      bool calculateForceConstantsDM();
      bool readForceConstantsFromVasprun(_xinput&);
      void symmetrizeForceConstantMatrices();
      void correctSumRules();

      //void printForceConstantMatrices(ostream&);  // OBSOLETE ME20200504 - not used
      //void printFCShellInfo(ostream&);  // OBSOLETE ME20200504 - not used

      // Direct method
      bool AUTO_GENERATE_PLUS_MINUS;
      bool USER_GENERATE_PLUS_MINUS;
      bool GENERATE_ONLY_XYZ;
      bool DISTORTION_SYMMETRIZE; //CO20190108
      double DISTORTION_MAGNITUDE;
      bool DISTORTION_INEQUIVONLY; //CO20190108
      // For each inequivalent atom, there is a set of unique distortions
      vector<vector<xvector<double> > > _uniqueDistortions;
      // For each inequivalent atom and unique distortion, there is a field
      // of forces (for each atom of the supercell)
      vector<vector<vector<xvector<double> > > > _uniqueForces;
      vector<vector<bool> > vvgenerate_plus_minus;  //ME20191029

      void estimateUniqueDistortions(const xstructure&,
          vector<vector<xvector<double> > >&);
      void testDistortion(const xvector<double>&, const vector<_sym_op>&,
          vector<xvector<double> >&,
          vector<xvector<double> >&,
          bool integrate_equivalent_distortions=true);  //CO20190114
      bool needMinus(uint atom_index, uint distortion_index, bool inequiv_only=true);  //CO //CO20190218
      bool calculateForceFields();  // ME20190412  //ME20191029
      void completeForceFields();
      void projectToCartesianDirections();
      void buildForceConstantMatrices();

      // Born charges + dielectric tensor
      bool runVASPCalculationsBE(_xinput&, _aflags&, _kflags&, _xflags&, string&, uint);
      bool calculateBornChargesDielectricTensor(const _xinput&);  // ME20191029
      void readBornEffectiveChargesFromAIMSOUT(void);
      void readBornEffectiveChargesFromOUTCAR(const _xinput&);  //ME20190113
      void symmetrizeBornEffectiveChargeTensors(void);
      void readDielectricTensorFromAIMSOUT(void);
      void readDielectricTensorFromOUTCAR(const _xinput&);  // ME20190113

    public:
      ForceConstantCalculator(ostream& oss=std::cout);
      ForceConstantCalculator(Supercell&, ofstream&, ostream& os=std::cout);
      ForceConstantCalculator(Supercell&, const aurostd::xoption&, ofstream&, ostream& os=std::cout);
      ForceConstantCalculator(const ForceConstantCalculator&);
      ForceConstantCalculator& operator=(const ForceConstantCalculator&);
      ~ForceConstantCalculator();
      void clear(Supercell&);
      void initialize(const aurostd::xoption&, ofstream&, ostream& oss=std::cout);
      void initialize(const aurostd::xoption&);
      void initialize(const xvector<int>&, const xstructure& xs, bool=true, bool=true);
      void initialize(ostream& oss);
      void initialize(ofstream& mf, ostream& oss);

      bool runVASPCalculations(_xinput&, _aflags&, _kflags&, _xflags&, string&);

      bool run();  // ME20191029
      void hibernate();

      const vector<vector<xmatrix<double> > >& getForceConstants() const;
      const vector<xmatrix<double> >& getBornEffectiveChargeTensor() const;
      const xmatrix<double>& getDielectricTensor() const;
      bool isPolarMaterial() const;

      string _directory;

      void writeHarmonicIFCs(const string&);
      void writeBornChargesDielectricTensor(const string&);
      void writeDYNMAT(const string&);
      void saveState(const string&);  // ME20200112
      void readFromStateFile(const string&);  // ME20200112
  };

}  // namespace apl

// ***************************************************************************

namespace apl {

  enum IPCFreqFlags {
    NONE = 0L,
    ALLOW_NEGATIVE = 1L << 1,
    OMEGA = 1L << 2,
    RAW = 1L << 3,  // eV/A/A/atomic_mass_unit
    HERTZ = 1L << 4,
    THZ = 1L << 5,
    RECIPROCAL_CM = 1L << 6,
    MEV = 1L << 7
  };
  inline IPCFreqFlags operator&(const IPCFreqFlags& __a, const IPCFreqFlags& __b) {
    return IPCFreqFlags(static_cast<int>(__a) & static_cast<int>(__b));
  }
  inline IPCFreqFlags operator|(const IPCFreqFlags& __a, const IPCFreqFlags& __b) {
    return IPCFreqFlags(static_cast<int>(__a) | static_cast<int>(__b));
  }
  inline IPCFreqFlags operator|=(IPCFreqFlags& __a, const IPCFreqFlags& __b) {
    return (__a = (__a | __b));
  }

}  // namespace apl

// ***************************************************************************

namespace apl {

  class AnharmonicIFCs;  // Forward declaration
  class PhononCalculator : public xStream {
    private:
      // USER PARAMETERS
      string _directory;  // for loggers
      int _ncpus;

      QMesh _qm;
      Supercell _supercell;

      // harmonic IFCs
      vector<vector<xmatrix<double> > > _forceConstantMatrices;

      // Stuff for polar materials
      bool _isPolarMaterial;
      // For each atom there is a matrix 3x3 of Born effective charge
      vector<xmatrix<double> > _bornEffectiveChargeTensor;
      // Dielectric tensor
      xmatrix<double> _dielectricTensor;
      // Precomputed values used in non-analytical term (Gonze)
      xmatrix<double> _inverseDielectricTensor;
      double _recsqrtDielectricTensorDeterminant;
      // Precomputed Ewald sum at Gamma point
      bool _isGammaEwaldPrecomputed;
      vector<xmatrix<xcomplex<double> > > _gammaEwaldCorr;
      // Anharmonic IFCs
      vector<vector<vector<double> > > anharmonicIFCs;
      vector<vector<vector<int> > > clusters;

      void copy(const PhononCalculator&);  // ME20191228
      void free();

      void readHarmonicIFCs(const string&);
      void readBornChargesDielectricTensor(const string&);

      xmatrix<xcomplex<double> > getNonanalyticalTermWang(const xvector<double>&);
      xmatrix<xcomplex<double> > getNonanalyticalTermWang(const xvector<double>&,
          vector<xmatrix<xcomplex<double> > >&, bool=true);  // ME20180829
      xmatrix<xcomplex<double> > getNonanalyticalTermGonze(const xvector<double>);
      xmatrix<xcomplex<double> > getEwaldSumDipoleDipoleContribution(const xvector<double>, bool = true);

      void calculateGroupVelocitiesThread(int, vector<vector<double> >&, vector<xmatrix<xcomplex<double> > >&, vector<vector<xvector<double> > >&);

    public:
      PhononCalculator(ostream& oss=std::cout);
      PhononCalculator(ofstream&, ostream& oss=std::cout);
      PhononCalculator(const PhononCalculator&);
      PhononCalculator& operator=(const PhononCalculator&);
      ~PhononCalculator();
      void clear();
      void initialize(ostream& oss);
      void initialize(ofstream& mf, ostream& oss);

      string _system;  // ME20190614 - for VASP-style output files

      // Getter functions
      Supercell& getSupercell();
      QMesh& getQMesh();
      const xstructure& getInputCellStructure() const;
      const xstructure& getSuperCellStructure() const;
      uint getNumberOfBranches() const;
      string getDirectory() const;
      int getNCPUs() const;
      bool isPolarMaterial() const;  // ME20200206
      const vector<vector<xmatrix<double> > >& getHarmonicForceConstants() const;
      const vector<vector<double> >& getAnharmonicForceConstants(int) const;
      const vector<vector<int> >& getClusters(int) const;

      // Set functions
      void setDirectory(const string&);
      void setNCPUs(const _kflags&);
      void setPolarMaterial(bool);

      // Initializers
      void initialize_qmesh(const vector<int>&, bool=true, bool=true);
      void initialize_qmesh(const xvector<int>&, bool=true, bool=true);
      void initialize_supercell(const xstructure&, bool verbose=true);//AS20200908
      void initialize_supercell(const string&);

      // IFCs
      void setHarmonicForceConstants(const ForceConstantCalculator&);
      void setHarmonicForceConstants(const vector<vector<xmatrix<double> > > &);//AS20201204
      void setHarmonicForceConstants(const vector<vector<xmatrix<double> > > &,
          const vector<xmatrix<double> >&, const xmatrix<double> &, bool isPolar=true);//AS20201204
      void awake();
      void setAnharmonicForceConstants(const AnharmonicIFCs&);
      void readAnharmonicIFCs(string);

      // Dynamical Matrix/Frequencies
      xvector<double> getEigenvalues(const xvector<double>&, const xvector<double>&,
          xmatrix<xcomplex<double> >&, vector<xmatrix<xcomplex<double> > >&, bool=true);  // ME20180827
      xmatrix<xcomplex<double> > getDynamicalMatrix(const xvector<double>&);
      xmatrix<xcomplex<double> > getDynamicalMatrix(const xvector<double>&, const xvector<double>&);  // ME20200206
      xmatrix<xcomplex<double> > getDynamicalMatrix(const xvector<double>&, const xvector<double>&,
          vector<xmatrix<xcomplex<double> > >&, bool=true);  // ME20180827
      xvector<double> getFrequency(const xvector<double>&, const IPCFreqFlags&);  // ME20180827
      xvector<double> getFrequency(const xvector<double>&, const xvector<double>&, const IPCFreqFlags&);  // ME20200206
      xvector<double> getFrequency(const xvector<double>&, const IPCFreqFlags&, xmatrix<xcomplex<double> >&);  // ME20190624
      xvector<double> getFrequency(const xvector<double>&, const xvector<double>&, const IPCFreqFlags&, xmatrix<xcomplex<double> >&);  // ME20200206
      xvector<double> getFrequency(const xvector<double>&, const IPCFreqFlags&, xmatrix<xcomplex<double> >&,
          vector<xvector<double> >&, bool=true);  // ME20180827
      xvector<double> getFrequency(const xvector<double>&, const xvector<double>&, const IPCFreqFlags&, xmatrix<xcomplex<double> >&,
          vector<xvector<double> >&, bool=true);  // ME20200206
      double getFrequencyConversionFactor(IPCFreqFlags, IPCFreqFlags);

      // Group velocities
      vector<vector<xvector<double> > > calculateGroupVelocitiesOnMesh();
      vector<vector<xvector<double> > > calculateGroupVelocitiesOnMesh(vector<vector<double> >&);
      vector<vector<xvector<double> > > calculateGroupVelocitiesOnMesh(vector<vector<double> >&, vector<xmatrix<xcomplex<double> > >&);
      void writeGroupVelocitiesToFile(const string&, const vector<vector<xvector<double> > >&);
      void writeGroupVelocitiesToFile(const string&, const vector<vector<xvector<double> > >&, const vector<vector<double> >&, const string& unit="THz");
  };

}  // namespace apl

// ***************************************************************************

namespace apl {

  class PathBuilder {
    public:
      enum StoreEnumType { RECIPROCAL_LATTICE,
        CARTESIAN_LATTICE };
      enum ModeEnumType { SINGLE_POINT_MODE,
        COUPLE_POINT_MODE };

    private:
      std::vector<aurostd::xvector<double> > _path;
      std::vector<aurostd::xvector<double> > _points;
      std::vector<std::string> _labels;
      aurostd::xmatrix<double> reciprocalLattice;
      aurostd::xmatrix<double> cartesianLattice;
      int _pointsVectorDimension;
      int _pointsVectorStartingIndex;
      uint _nPointsPerSubPath;
      ModeEnumType _mode;
      StoreEnumType _store;

    private:
      void buildPath();
      void free();
      void copy(const PathBuilder&);

    public:
      PathBuilder();
      PathBuilder(ModeEnumType);
      PathBuilder(const PathBuilder&);
      PathBuilder& operator=(const PathBuilder&);
      ~PathBuilder();
      void clear();
      void addPoint(const std::string& l, int dim, ...);
      void addPoint(const std::string&, const aurostd::xvector<double>&);
      void transform(const aurostd::xmatrix<double>&);
      void pointsAreDefinedFor(const xstructure&, StoreEnumType);
      void transformPointsFor(const xstructure&, StoreEnumType);
      void defineCustomPoints(const string&,const string&,const Supercell&,bool CARESTIAN_COORDS=false);
      void takeAflowElectronicPath(const string&,const Supercell&);//, const xstructure&, const xstructure&);
      void setMode(ModeEnumType);
      void setStore(StoreEnumType);
      const StoreEnumType& getStore() const;  //ME20190614
      void setDensity(int);
      int getDensity();
      uint getPathSize();
      uint getPointSize();
      aurostd::xvector<double> getPoint(uint);
      uint getPointIndexOnPath(uint);
      std::string getPointLabel(uint);
      std::vector<aurostd::xvector<double> > getPath();
      std::vector<aurostd::xvector<double> > getPath(ModeEnumType, const string&);
      double getPathLength();
      double getPathLength(uint);
      xKPOINTS createKPOINTS(const Supercell&);  //ME20190614
  };
}  // namespace apl

// ***************************************************************************

namespace apl {

  class PhononDispersionCalculator {
    private:
      PhononCalculator* _pc;
      bool _pc_set;
      PathBuilder _pb;
      void copy(const PhononDispersionCalculator&);
      void free();
      std::vector<xvector<double> > _qpoints;
      std::vector<xvector<double> > _freqs;
      IPCFreqFlags _frequencyFormat;
      double _temperature;  // ME20190614
      //[OBSOLETE PN20180705]vector<double> path;       //[PINKU]
      //[OBSOLETE PN20180705]vector<int> path_segment;  //[PINKU]
      void calculateInOneThread(int);
      bool isExactQPoint(const xvector<double>&, const xmatrix<double>&);
      string _system;

    public:
      PhononDispersionCalculator();
      PhononDispersionCalculator(PhononCalculator&);
      PhononDispersionCalculator(const PhononDispersionCalculator&);
      PhononDispersionCalculator& operator=(const PhononDispersionCalculator&);
      ~PhononDispersionCalculator();
      void clear(PhononCalculator&);
      void initPathCoords(const string&,const string&,int,bool=false);  //CO20180406
      void initPathLattice(const string&, int);
      void setPath(const string&);
      void calc(const IPCFreqFlags);
      void writePDIS(const string&);
      std::vector<xvector<double> > get_qpoints() { return _qpoints; }  //[PN]
      //ME20190614 START
      xEIGENVAL createEIGENVAL();
      void writePHEIGENVAL(const string&);
      void writePHKPOINTS(const string&);
      xKPOINTS getPHKPOINTS();//AS20201110
      //ME20190614 STOP
      //[OBSOLETE PN20180705]std::vector<double> get_path() { return path; }                   //[PN]
      //[OBSOLETE PN20180705]std::vector<int> get_path_segment() { return path_segment; }      //[PN]
  };

}  // namespace apl

// ***************************************************************************

namespace apl {

  class DOSCalculator {
    protected:
      PhononCalculator* _pc;
      bool _pc_set;
      string _bzmethod;  //ME20190423
      std::vector<aurostd::xvector<double> > _qpoints;
      //std::vector<int> _qweights;  OBSOLETE ME20190423
      std::vector<aurostd::xvector<double> > _freqs;
      double _minFreq;
      double _maxFreq;
      double _stepDOS;
      double _halfStepDOS;
      std::vector<double> _bins;
      std::vector<double> _dos;
      std::vector<double> _idos;  //ME20190614
      std::vector<xmatrix<xcomplex<double> > > _eigen;  //ME20190624 - eigenvectors for projected DOS
      std::vector<vector<vector<double> > > _projectedDOS; //ME20190614 - projectedDOS.at(atom).at(direction).at(frequency)
      std::vector<xvector<double> > _projections;  //ME20190626 - the projection directions for the DOS in Cartesian coordinates
      double _temperature;  //ME20190614
      //CO START
      //private:
      void copy(const DOSCalculator&);
      void calculateInOneThread(int);
      //CO END
      void calculateFrequencies();
      void smearWithGaussian(vector<double>&, vector<double>&, double, double);  //ME20190614
      void calcDosRS();
      void calcDosLT();

    public:
      DOSCalculator();
      DOSCalculator(PhononCalculator&, const xoption&);//AS20201203 input parameters are passed via xoption
      DOSCalculator(const DOSCalculator&);
      DOSCalculator& operator=(const DOSCalculator&);
      ~DOSCalculator();
      void clear(PhononCalculator&);
      void initialize(const xoption&);//AS20201203 input parameters are passed via xoption
      void calc(int, bool VERBOSE=true);
      void calc(int, double, bool VERBOSE=true);
      void calc(int, double, double, double, bool VERBOSE=true);  //ME20200203
      void writePDOS(const string&);
      void writePDOS(string, string);  //[PN]
      xDOSCAR createDOSCAR() const;  //ME20190614
      void writePHDOSCAR(const string&);  //ME20190614
      // Interface IDOSCalculator
      const std::vector<double>& getBins() const;  //ME20200108 - added const
      const std::vector<double>& getDOS() const;   //ME20200108 - added const
      const std::vector<double>& getIDOS() const;  //ME20200210
      const std::vector<xvector<double> >& getFreqs() const;  //AS20200312
      bool hasImaginaryFrequencies() const;  //ME20200108 - added const
      double getMinFreq() const;  // ME20210927
      double getMaxFreq() const;  // ME20210927
      const xstructure& getInputStructure() const;  //ME20210927
      uint getNumberOfBranches() const;  //ME20210927
      string _system;
    private:
      void free();
  };

}  // namespace apl

// ***************************************************************************

namespace apl {

  enum ThermalPropertiesUnits { eV,
    meV,
    ueV,
    eVK,
    meVK,
    ueVK,
    kB };

  class ThermalPropertiesCalculator : public xStream {
    private:
      std::vector<double> _freqs_0K;
      std::vector<double> _dos_0K;
      string system;

      void free();
      void copy(const ThermalPropertiesCalculator&);

      double getStepDOS(const vector<double>&);
      double getScalingFactor(const ThermalPropertiesUnits&);

    public:
      ThermalPropertiesCalculator(ostream& oss=std::cout);
      ThermalPropertiesCalculator(ofstream&, ostream& os=std::cout);
      ThermalPropertiesCalculator(const DOSCalculator&, ofstream&, const string& directory="./", ostream& os=std::cout);
      ThermalPropertiesCalculator(const xDOSCAR&, ofstream&, const string& directory="./", ostream& os=std::cout);
      ThermalPropertiesCalculator(const ThermalPropertiesCalculator&);
      ThermalPropertiesCalculator& operator=(const ThermalPropertiesCalculator&);
      ~ThermalPropertiesCalculator();
      void clear();
      void initialize(ostream& oss);
      void initialize(ofstream& mf, ostream& oss);

      uint natoms;
      vector<double> temperatures;
      vector<double> Cv;
      vector<double> Fvib;
      vector<double> Svib;
      vector<double> U;
      double U0;
      string _directory;

      void initialize(const xDOSCAR&, ofstream&, ostream& oss=std::cout);
      void initialize(const vector<double>&, const vector<double>&, ofstream&, const string& system="", ostream& oss=std::cout);
      void initialize(const vector<double>&, const vector<double>&, const string& system="");
      void calculateThermalProperties(double, double, double);
      void addPoint(double, const xDOSCAR&);
      void addPoint(double, const vector<double>&, const vector<double>&);

      double getZeroPointEnergy();
      double getInternalEnergy(double, ThermalPropertiesUnits=apl::meV);
      double getInternalEnergy(double, const vector<double>&, const vector<double>&, ThermalPropertiesUnits=apl::meV);
      double getVibrationalFreeEnergy(double, ThermalPropertiesUnits=apl::meV);
      double getVibrationalFreeEnergy(double, const vector<double>&, const vector<double>&, ThermalPropertiesUnits=apl::meV);
      double getVibrationalEntropy(double, ThermalPropertiesUnits=apl::meV);
      double getVibrationalEntropy(double, const vector<double>&, const vector<double>&, ThermalPropertiesUnits=apl::kB);
      double getVibrationalEntropy(double, double, double, ThermalPropertiesUnits=apl::kB);
      double getIsochoricSpecificHeat(double, ThermalPropertiesUnits=apl::kB);
      double getIsochoricSpecificHeat(double, const vector<double>&, const vector<double>&, ThermalPropertiesUnits=apl::kB);

      void writePropertiesToFile(string, filetype ft=txt_ft);
      void addToAPLOut(stringstream&);
      string getPropertiesFileString(filetype ft=txt_ft);
  };
}  // namespace apl

// ***************************************************************************

namespace apl {

  class AtomicDisplacements {
    protected:
      PhononCalculator* _pc;
      bool _pc_set;

    private:
      void free();
      void copy(const AtomicDisplacements&);

      vector<vector<vector<xvector<xcomplex<double> > > > > _eigenvectors;
      vector<vector<double> > _frequencies;
      vector<vector<xmatrix<xcomplex<double> > > > _displacement_matrices;
      vector<vector<vector<xvector<xcomplex<double> > > > > _displacement_modes;
      vector<_qpoint> _qpoints;
      vector<double> _temperatures;

      void calculateEigenvectors();
      void calculateEigenvectorsInThread(int);
      void calculateMeanSquareDisplacementMatrices();
      void calculateModeDisplacements();
      double getOccupationNumber(double, double);

    public:
      AtomicDisplacements();
      AtomicDisplacements(PhononCalculator&);
      AtomicDisplacements(const AtomicDisplacements&);
      AtomicDisplacements& operator=(const AtomicDisplacements&);
      ~AtomicDisplacements();
      void clear(PhononCalculator&);

      void calculateMeanSquareDisplacements(double, double, double);
      void calculateModeDisplacements(const vector<xvector<double> >& qpts, bool=true);

      const vector<double>& getTemperatures() const;
      const vector<vector<xmatrix<xcomplex<double> > > >& getDisplacementMatrices() const;
      vector<vector<xvector<double> > > getDisplacementVectors() const;
      const vector<vector<vector<xvector<xcomplex<double> > > > >& getModeDisplacements() const;

      vector<vector<vector<double> > > createDisplacementsXcrysden(const Supercell&, double, int, int, int);
      void getOrientedDisplacementsVsim(xstructure&, vector<vector<vector<xvector<xcomplex<double> > > > >&, double);

      void writeMeanSquareDisplacementsToFile(string);
      void writeSceneFileXcrysden(string, const xstructure&, const vector<vector<vector<double> > >&, int);
      void writeSceneFileVsim(string, const xstructure&, const vector<vector<vector<xvector<xcomplex<double> > > > >&);
  };

  void createAtomicDisplacementSceneFile(const aurostd::xoption& vpflow, ostream& oss=std::cout);
  void createAtomicDisplacementSceneFile(const aurostd::xoption& vpflow, ofstream&, ostream& oss=std::cout);
}

// ***************************************************************************
// END: Automatic Phonon Library (APL)
// ***************************************************************************

// ***************************************************************************
// BEGIN ME: Automatic Anharmonic Phonon Library (AAPL)
// ***************************************************************************

namespace apl {

  // _cluster holds a single cluster
  struct _cluster {
    vector<int> atoms;  // List of atoms inside the cluster
    int fgroup;  // Index pointing to the factor group that transforms the cluster into another cluster
    int permutation;  // Index pointing to the permutation that transforms the cluster into another cluster
  };

  // _ineq_distortions contains a list of inequivalent distortion and its equivalent
  // distortions for a given set of atoms
  struct _ineq_distortions {
    vector<int> atoms;  // A list of atoms involved in the distortions
    vector<int> clusters;  // A list of cluster sets that use these distortions for their force constant calculations
    vector<vector<vector<int> > > distortions; // Map of distortions. The distortions vectors need to be defined elsewhere. 
    vector<vector<int> > rotations;  // The factor group that holds the rotation to transform the distortions
    vector<vector<vector<int> > > transformation_maps;  // A map containing the transformation of the atoms for each distortion
  };

  // _linearCombinations is a structure to store linear combinations.
  struct _linearCombinations {
    vector<vector<int> > indices;  // Cartesian indices of each linear combination
    vector<vector<double> > coefficients;  // Coefficients of each linear combination
    vector<int> independent;  // The linearly independent values
    vector<int> dependent;  // The linearly dependent values
    // indep2depMap maps the independent coefficients to the coefficients that
    // depend on them. This is used for the IFC correction method.
    vector<vector<int> > indep2depMap;
  };

  class ClusterSet : public xStream {
    // See aflow_aapl_cluster.cpp for detailed descriptions of the functions
    public:
      ClusterSet(ostream& oss=std::cout);
      ClusterSet(ofstream&, ostream& oss=std::cout);
      ClusterSet(const Supercell&, int, int, double, const string&, ofstream&, ostream& oss=std::cout);  // Constructor
      ClusterSet(const string&, const Supercell&, int, int, double, const string&, ofstream&, ostream& oss=std::cout);  // From file
      ClusterSet(const ClusterSet&);  // Constructor from another ClusterSet instance
      ~ClusterSet();  // Destructor
      const ClusterSet& operator=(const ClusterSet&);  // Copy constructor
      void clear();
      void initialize(const Supercell&, int, int, double, ofstream&, ostream& oss=std::cout);
      void initialize(const Supercell&, int, int, double);
      void initialize(ostream& oss);
      void initialize(ofstream& mf, ostream& oss);
      void readClusterSetFromFile(const string&);

      vector<_cluster> clusters;
      vector<vector<int> > coordination_shells;  // Contains all coordinate shells. Central atoms is index 0.
      double cutoff;  // Cutoff radius in Angstroms
      string directory;  // Directory for logging
      vector<xvector<double> > distortion_vectors;  // List of distortion vectors
      vector<_ineq_distortions> higher_order_ineq_distortions;  //ME20190531 - for 3rd derivatives of higher order processes
      vector<vector<int> > ineq_clusters;  // Clusters rearranged into sets of equivalent clusters.  //ME20190520
      vector<_ineq_distortions> ineq_distortions; // List of inequivalent distortions
      vector<_linearCombinations> linear_combinations;  // List of linear combinations of the IFCs
      int nifcs;  // Number of force constants for each set of atoms.
      int order;  // Order of the cluster, i.e. the order of the force constant to be calculated.
      xstructure pcell;  // Structure of the primitive cell.
      vector<int> pc2scMap;  // Atom map from the primitive cell to the supercell.
      vector<vector<int> > permutations;  // List of possible permutations for the cluster
      xstructure scell;  // Structure of the supercell.
      vector<int> sc2pcMap;  // Atom map from the supercell to the primitive cell.
      xvector<int> sc_dim;  // Dimensions of the supercell.
      vector<vector<int> > symmetry_map;  // Symmetry atom map for the atoms in the clusters

      const _cluster& getCluster(const int& i) const;  //ME20190520
      void build();
      void buildDistortions();
      void writeClusterSetToFile(const string&);

    private:
      void free();
      void copy(const ClusterSet&);

      double getMaxRad(const xstructure&, int);
      void buildShells();
      vector<_cluster> buildClusters();
      vector<vector<int> > getSymmetryMap();
      vector<vector<int> > getPermutations(int);

      // Clusters
      void getInequivalentClusters(vector<_cluster>&, vector<vector<int> >&);
      int getNumUniqueAtoms(const vector<int>&);
      vector<int> getComposition(const vector<int>&);
      bool sameComposition(const vector<int>&, const vector<int>&);
      int equivalenceCluster(const vector<int>&, const vector<int>&,
          const vector<vector<int> >&, const vector<vector<int> >&);
      vector<int> translateToPcell(const vector<int>&, int);
      int comparePermutations(const vector<int>&, const vector<int>&);
      bool atomsMatch(const vector<int>&, const vector<int>&, const vector<int>&, const int&);
      void getSymOp(_cluster&, const vector<int>&);

      // Distortions
      vector<xvector<double> > getCartesianDistortionVectors();
      vector<_ineq_distortions> initializeIneqDists();
      int sameDistortions(const _cluster&, const vector<_ineq_distortions>&);
      vector<vector<int> > getTestDistortions(const vector<int>&);
      void getInequivalentDistortions(const vector<vector<int> >&, _ineq_distortions&);
      void appendDistortion(_ineq_distortions&, vector<int>,
          const int& eq=-1, const int& fg=-1);
      bool allZeroDistortions(const vector<int>&, const vector<int>&);
      bool allAtomsMatch(const int&, const vector<int>&);
      int equivalenceDistortions(const xmatrix<double>&, const vector<int>&,
          const vector<vector<vector<int> > >&, const vector<int>&);
      vector<int> getTransformationMap(const int&, const int&);
      vector<_ineq_distortions> getHigherOrderDistortions();

      // Linear Combinations
      vector<_linearCombinations> getLinearCombinations();
      vector<vector<int> > getInvariantSymOps(const _cluster&);
      vector<vector<double> > buildCoefficientMatrix(const vector<vector<int> >&);
      vector<vector<double> > getRREF(vector<vector<double> >);

      // File I/O
      string writeParameters();
      string writeInequivalentClusters();
      string writeClusters(const vector<_cluster>&);
      string writeLinearCombinations(const _linearCombinations&);
      string writeInequivalentDistortions();
      string writeIneqDist(const _ineq_distortions&);
      string writeHigherOrderDistortions();

      bool checkCompatibility(uint&, const vector<string>&);
      void readInequivalentClusters(uint&, const vector<string>&);
      vector<_cluster> readClusters(uint&, const vector<string>&);
      _linearCombinations readLinearCombinations(uint&, const vector<string>&);
      void readInequivalentDistortions(uint&, const vector<string>&);
      _ineq_distortions readIneqDist(uint&, const vector<string>&);
      void readHigherOrderDistortions(uint&, const vector<string>&);
  };
}  // namespace apl

// ***************************************************************************

namespace apl {

  class AnharmonicIFCs : public xStream {
    // See aflow_aapl_ifcs.cpp for detailed descriptions of the functions
    public:
      AnharmonicIFCs(ostream& oss=std::cout);
      AnharmonicIFCs(ofstream&, ostream& oss=std::cout);
      AnharmonicIFCs(const AnharmonicIFCs&);
      const AnharmonicIFCs& operator=(const AnharmonicIFCs&);
      ~AnharmonicIFCs();
      void clear();
      void initialize(const Supercell&, int, const aurostd::xoption&, ofstream&, ostream& oss=std::cout);
      void initialize(const Supercell&, int, const aurostd::xoption&);
      void initialize(ostream& oss);
      void initialize(ofstream& mf, ostream& oss);

      void setOptions(double, int, double, double, bool);
      const string& getDirectory() const;
      void setDirectory(const string&);
      int getOrder() const;

      bool runVASPCalculations(_xinput&, _aflags&, _kflags&, _xflags&);
      bool calculateForceConstants();
      const vector<vector<double> >& getForceConstants() const;
      vector<vector<int> > getClusters() const;
      void writeIFCsToFile(const string&);

    private:
      ClusterSet clst;

      vector<_xinput> xInputs;
      bool _useZeroStateForces;
      bool initialized;
      string directory;
      vector<vector<int> > cart_indices;  // A list of all Cartesian indices
      double distortion_magnitude;  // The magnitude of the distortions in Angstroms
      vector<vector<double> > force_constants;  // Symmetrized IFCs - ME20190520
      int max_iter;  // Number of iterations for the sum rules
      double mixing_coefficient;  // The mixing coefficient for the SCF procedure
      int order;  // The order of the IFCs
      double sumrule_threshold;  // Convergence threshold for the sum rules

      void free();
      void copy(const AnharmonicIFCs&);

      string buildRunName(const vector<int>&, const vector<int>&, int, int);
      void applyDistortions(_xinput&, const vector<xvector<double> >&, const vector<int>&, const vector<int>&, double=1.0);

      vector<vector<int> > getCartesianIndices();

      vector<vector<vector<xvector<double> > > > storeForces(vector<_xinput>&);
      vector<vector<xvector<double> > > getForces(int, int&, vector<_xinput>&);
      int getTransformedAtom(const vector<int>&, const int&);
      void addHigherOrderForces(vector<vector<vector<xvector<double> > > >&, int&, vector<_xinput>&);
      vector<vector<double> > calculateUnsymmetrizedIFCs(const vector<_ineq_distortions>&,
          const vector<vector<vector<xvector<double> > > >&);
      double finiteDifference(const vector<vector<xvector<double> > >&, int,
          const vector<int>&, const vector<int>&);

      // Symmetrization Functions
      vector<vector<double> > symmetrizeIFCs(vector<vector<double> >);
      typedef vector<std::pair<vector<int>, vector<double> > > tform;
      typedef vector<vector<vector<vector<int> > > > v4int;
      void getTensorTransformations(v4int&, vector<vector<tform> >&);
      vector<vector<int> > getReducedClusters();
      void applyLinCombs(vector<vector<double> >&);
      void transformIFCs(const vector<vector<tform> >&, vector<vector<double> >&);
      void applyPermutations(vector<int>, vector<vector<double> >&);
      void calcSums(const vector<vector<int> >&, const vector<vector<double> >&,
          vector<vector<double> >&, vector<vector<double> >&);
      void correctIFCs(vector<vector<double> >&, const vector<vector<double> >&,
          const vector<vector<double> >&,
          const vector<vector<int> >&, const v4int&);
      vector<double> getCorrectionTerms(const int&,
          const vector<vector<int> >&,
          const vector<vector<double> >&,
          const vector<vector<double> >&,
          const vector<vector<double> >&);
      uint findReducedCluster(const vector<vector<int> >&, const vector<int>&);

      // File I/O
      string writeParameters();
      string writeIFCs();
      bool checkCompatibility(uint&, const vector<string>&);
      vector<vector<double> > readIFCs(uint&, const vector<string>&);
  };

}  //namespace apl

// ***************************************************************************

namespace apl {

  class TCONDCalculator {
    // See aflow_aapl_tcond.cpp for detailed descriptions of the functions
    public:
      TCONDCalculator();
      TCONDCalculator(PhononCalculator&, const aurostd::xoption&);
      TCONDCalculator(const TCONDCalculator&);
      TCONDCalculator& operator=(const TCONDCalculator&);
      ~TCONDCalculator();
      void clear(PhononCalculator&);
      void initialize(const aurostd::xoption&);

      double boundary_grain_size;
      bool calc_boundary;
      bool calc_cumulative;
      bool calc_isotope;
      bool calc_rta_only;
      vector<xmatrix<xcomplex<double> > > eigenvectors;  // The eigenvectors at each q-point
      vector<vector<double> > freq;  // The frequencies at each q-point
      vector<vector<xvector<double> > > gvel;  // The group velocities
      int nBranches;  // The number of branches in the phonon spectrum
      int nIQPs;  // The total number of irreducible q-points in the grid
      int nQPs;  // The total number of q-points in the grid
      vector<vector<vector<int> > > processes;  // The sign, q-point and branch indices of the scattering processes
      vector<vector<double> > intr_trans_probs;  // The intrinsic transition probabilities
      vector<vector<vector<int> > > processes_iso;  // The q-point and branch indices of the isotope scattering processes
      vector<vector<double> > intr_trans_probs_iso;  // The intrinsic transition probabilities for isotope processes
      // Scattering rates
      vector<vector<vector<double> > > scattering_rates_total;  // total
      vector<vector<vector<double> > > scattering_rates_anharm;  // anharmonic scattering
      vector<vector<double> > scattering_rates_boundary;  // boundary scattering
      vector<vector<double> > scattering_rates_isotope;  // isotope scattering

      vector<vector<vector<vector<double> > > > phase_space;  // Scattering phase space

      vector<vector<double> > grueneisen_mode;  // Mode Grueneisen parameter
      vector<double>  grueneisen_avg;  // Averaged Grueneisen parameter for each temperature

      vector<double> temperatures;  // The temperatures for the thermal conductivity calculations
      vector<xmatrix<double> > thermal_conductivity;  // The thermal conductivity values

      void calculateThermalConductivity();
      void calculateGrueneisenParameters();
      void writeOutputFiles(const string&);

    private:
      PhononCalculator* _pc;  // Reference to the phonon calculator
      QMesh* _qm;
      bool _pc_set;
      bool _initialized;

      void free();
      void copy(const TCONDCalculator&);

      vector<vector<double> > calculateModeGrueneisen(const vector<vector<vector<xcomplex<double> > > >& phases);
      double calculateAverageGrueneisen(double T, const vector<vector<double> >&);

      void getWeightsLT(double, const vector<double>&, vector<double>&);
      void calculateTransitionProbabilities();
      vector<vector<vector<xcomplex<double> > > > calculatePhases(bool=false);
      void calculateTransitionProbabilitiesPhonon(int, vector<vector<vector<vector<double> > > >&,
          const vector<vector<vector<xcomplex<double> > > >&);
      void calculateTransitionProbabilitiesIsotope(int);
      vector<vector<double> > calculateTransitionProbabilitiesBoundary();
      void getProcess(const vector<int>&, vector<int>&, vector<int>&, int&);
      xmatrix<double> calculateThermalConductivityTensor(double, vector<vector<double> >&, vector<vector<double> >&);
      vector<vector<double> > getOccupationNumbers(double);
      vector<vector<double> > calculateAnharmonicRates(const vector<vector<double> >&);
      vector<vector<double> > calculateTotalRates(const vector<vector<double> >&, vector<vector<double> >&);
      double getOccupationTerm(const vector<vector<double> >&, int, const vector<int>&, const vector<int>&);
      void calcAnharmRates(int, const vector<vector<double> >&, vector<vector<double> >&);
      vector<vector<xvector<double> > > getMeanFreeDispRTA(const vector<vector<double> >&);
      xmatrix<double> calcTCOND(double, const vector<vector<double> >&,
          const vector<vector<xvector<double> > >&);
      void getMeanFreeDispFull(const vector<vector<double> >&,
          const vector<vector<double> >&, vector<vector<xvector<double> > >&);
      void calculateDelta(int, const vector<vector<double> >&, const vector<vector<xvector<double> > >&, vector<vector<xvector<double> > >&);
      void correctMFD(const vector<vector<double> >&, const vector<vector<xvector<double> > >&, vector<vector<xvector<double> > >&);

      void writeTempIndepOutput(const string&, string, const string&, const vector<vector<double> >&);
      void writeTempDepOutput(const string&, string, const string&, const vector<double>&, const vector<vector<vector<double> > >&);
      void writeDataBlock(stringstream&, const vector<vector<double> >&);
      void writePhaseSpace(const string&);
      void writeGrueneisen(const string&);
      void writeThermalConductivity(const string&);
  };

}  // namespace apl

// ***************************************************************************
// END ME: Automatic Anharmonic Phonon Library (AAPL)
// ***************************************************************************

// ***************************************************************************
// BEGIN AS: Quasi-Harmonic Approximation (QHA)
// ***************************************************************************

//AS20200513 BEGIN
#define QHA_ARUN_MODE "QHA" // used in filename
namespace apl
{
  enum EOSmethod {EOS_MURNAGHAN, EOS_SJ, EOS_BIRCH_MURNAGHAN2, EOS_BIRCH_MURNAGHAN3,
    EOS_BIRCH_MURNAGHAN4};
  enum QHAmethod {QHA_CALC, QHA3P_CALC, SCQHA_CALC, QHANP_CALC};
  enum QHAtype   {QHA_FD, QHA_EOS, QHA_TE};

  bool QHA_Get_AflowInName(string &AflowInName, const string &directory_LIB);
  string EOSmethod2label(EOSmethod eos_method);//AS20210518
  string QHAmethod2label(QHAmethod qha_method);//AS20210518
  bool hasImaginary(const string& filename, const string &QHA_method);//AS20210813
  //void linkAPLtoQHA();//AS20201216 OBSOLETE

  /// Calculates QHA-related properties
  class QHA : public xStream {
    public:
      QHA(ostream& oss=std::cout);
      QHA(const QHA& qha);
      QHA(const xstructure &struc, _xinput &xinput, xoption &qhaopts, xoption &aplopts,
          ofstream &FileMESSAGE, ostream &oss=std::cout);
      void initialize(ostream& oss);
      void initialize(ofstream& mf, ostream& oss);
      void initialize(const xstructure &struc, _xinput &xinput, xoption &qhaopts,
          xoption &aplopts, ofstream &FileMESSAGE, ostream &oss);
      ~QHA();
      const QHA& operator=(const QHA &qha);
      void run(_xflags &xflags, _aflags &aflags, _kflags &kflags);
      void clear();
      double calcFrequencyFit(double V, xvector<double> &xomega);
      double calcGrueneisen(double V, xvector<double> &xomega);
      double calcGrueneisen(double V, xvector<double> &xomega, double &w);
      double calcGrueneisenFD(const xvector<double> &xomega);
      void   calcCVandGP(double T, double &CV, double &GP);
      void   calcCVandGPfit(double T, double V, double &CV, double &GP);
      double calcGPinfFit(double V);
      double calcVibFreeEnergy(double T, int id);
      double calcFreeEnergyFit(double T, double V, EOSmethod eos_method,
          QHAmethod method, uint contrib);
      double calcElectronicFreeEnergy(double T, int id);
      double calcChemicalPotential(double T, int Vid);
      double calcIDOS(double e, double T, xEIGENVAL &eig);
      xvector<double> calcElectronicFreeEnergySommerfeld(double T);
      xvector<double> calcElectronicSpecificHeatSommerfeld(double T);
      xvector<double> calcFreeEnergy(double T, QHAmethod qha_method, uint contrib);
      xvector<double> calcDOSatEf();
      double calcInternalEnergyFit(double T, double V, EOSmethod method);
      xvector<double> fitToEOSmodel(xvector<double> &V, xvector<double> &E, EOSmethod method);
      xvector<double> fitToEOSmodel(xvector<double> &E, EOSmethod method);
      double evalEOSmodel(double V, const xvector<double> &p, EOSmethod eos_method);
      double calcBulkModulus(double V, const xvector<double> &parameters, EOSmethod method);
      double calcBprime(double V, const xvector<double> &parameters, EOSmethod method);
      double calcEOS2Pressure(double V, const xvector<double> &parameters, EOSmethod method);
      double calcEquilibriumVolume(const xvector<double> &parameters, EOSmethod method);
      double calcEntropy(double T, double V, EOSmethod eos_method,
          QHAmethod method, uint contrib);
      double getEqVolumeT(double T, EOSmethod eos_method, QHAmethod method, uint contrib);
      double calcThermalExpansion(double T, EOSmethod eos_method, QHAmethod method, uint contrib);
      double calcIsochoricSpecificHeat(double T, double V, EOSmethod eos_method, 
          QHAmethod qha_method, uint contrib);
      // QHA3P and SCQHA and QHANP
      double extrapolateFrequency(double V, const xvector<double> &xomega, QHAmethod qha_method);
      double extrapolateGrueneisen(double V, const xvector<double> &xomega, QHAmethod qha_method);
      // QHA3P
      double calcVibFreeEnergyTaylorExpansion(double T, int Vid, QHAmethod qha_method);
      double calcInternalEnergyTaylorExpansion(double T, double V, QHAmethod qha_method);
      // SCQHA
      double calcVPgamma(double T, double V);
      double calcSCQHAequilibriumVolume(double T, EOSmethod method);
      double calcSCQHAequilibriumVolume(double T, double Vguess, xvector<double> &fit_params, EOSmethod method);
      void   runSCQHA(EOSmethod method, bool all_iterations_self_consistent=true,
          const string &directory=".");
      // output
      void   writeThermalProperties(EOSmethod eos_method, QHAmethod qha_method, 
          const string &directory=".");
      void   writeFVT(const string &directory=".");
      void   writeGPpath(double V, const string &directory=".");
      void   writeAverageGPfiniteDifferences(const string &directory=".");
      void   writeGPmeshFD(const string &directory=".");
      void   writeFrequencies(const string &directory=".");
      void   writeTphononDispersions(EOSmethod eos_method, QHAmethod qha_method,
          const string &directory=".");
      void   writeQHAresults(const string &directory=".");
      // members
      xoption apl_options;
      xoption qha_options;
      string system_title;
      double EOS_volume_at_equilibrium;
      double EOS_energy_at_equilibrium;
      double EOS_bulk_modulus_at_equilibrium;
      double EOS_Bprime_at_equilibrium;
    private:
      bool isEOS;
      bool isGP_FD;
      bool ignore_imaginary;
      bool isQHA, isQHA3P, isSCQHA, isQHANP;
      bool isInitialized;
      bool includeElectronicContribution;
      bool doSommerfeldExpansion;
      int Ntemperatures;
      int N_GPvolumes;   ///< number of volumes/calculations for finite difference calc
      int N_EOSvolumes;  ///< number of volumes/calculations for EOS calc
      int N_QHANPvolumes; ///< number of volumes/calculations for QHANP calc
      int Nbranches;       ///< number of phonon dispersion branches
      int NatomsOrigCell;  ///< number of atoms in original cell
      int Nelectrons;
      int TaylorExpansionOrder;
      double gp_distortion;
      //int NatomsSupercell; ///< number of atoms in supercell
      xstructure origStructure;
      vector<double> Temperatures;
      vector<int> ph_disp_temperatures;///< temperatures for T-dependent phonon dispersions
      vector<double> GPvolumes; ///< a set of volumes for FD Grueneisen calculation
      vector<double> EOSvolumes; ///< a set of volumes for EOS calculation
      vector<double> QHANPvolumes; ///< a set of volumes for QHANP calculation
      vector<double> coefGPVolumes; ///< multiplication coefficient w.r.t initial volume
      vector<double> coefEOSVolumes;
      vector<double> coefQHANPVolumes;
      xvector<double> DOS_Ef;
      // data necessary to calculate thermodynamic properties
      vector<double> Efermi_V; ///< Fermi energy vs V
      vector<double> E0_V;     ///< total energy vs V
      vector<xEIGENVAL> static_eigvals;
      vector<xIBZKPT>   static_ibzkpts;
      vector<vector<double> > pdos_V; ///< phonon DOS
      vector<int> qpWeights;
      vector<xvector<double> > qPoints;
      // data needed for Grueneisen parameter calculation
      xmatrix<double> gp_fit_matrix;
      vector<vector<vector<double> > > omegaV_mesh;
      vector<vector<vector<double> > > omegaV_mesh_EOS;
      vector<vector<vector<double> > > omegaV_mesh_QHANP;
      vector<xEIGENVAL> gp_ph_dispersions;
      vector<xEIGENVAL> eos_ph_dispersions;
      vector<ThermalPropertiesCalculator> eos_vib_thermal_properties;
      //
      vector<string> subdirectories_apl_eos;
      vector<string> subdirectories_apl_gp;
      vector<string> subdirectories_apl_qhanp;
      vector<string> subdirectories_static;
      vector<string> arun_runnames_apl_eos;
      vector<string> arun_runnames_apl_gp;
      vector<string> arun_runnames_apl_qhanp;
      vector<string> arun_runnames_static;
      _xinput xinput;
      string currentDirectory;
      // methods
      // related to static DFT calculations
      void createSubdirectoriesStaticRun(const _xflags &xflags, const _aflags &aflags,
          const _kflags &kflags, const vector<vector<bool> > &list);
      int  checkStaticCalculations(vector<vector<bool> > &file_is_present);
      void printMissingStaticFiles(const vector<vector<bool> > & list,
          const vector<string> &subdirectories);
      bool readStaticCalculationsData();
      // related to APL calculations
      void createSubdirectoriesAPLRun(const _xflags &xflags, const _aflags &aflags,
          const _kflags &kflags, const vector<vector<bool> > &list, QHAtype qhatype);
      int  checkAPLCalculations(vector<vector<bool> > &file_is_present, QHAtype qhatype);
      void printMissingAPLFiles(const vector<vector<bool> > & list, QHAtype qhatype);
      bool readAPLCalculationData(const vector<string> &subdirectories, _kflags &kflags,
          QHAtype type);
      bool runAPL(_xflags &xflags, _aflags &aflags, _kflags &kflags, QHAtype qhatype);
      // mandatory
      void free();
      void copy(const QHA &qha);
  };
}
//AS20200513 END

// ***************************************************************************
// END AS: Quasi-Harmonic Approximation (QHA)
// ***************************************************************************

#endif  // _AFLOW_APL_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
