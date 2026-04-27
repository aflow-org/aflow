// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// MAKEFILE FOR AFLOW_APL
// Written by Michal Jahnatek

#ifndef _AFLOW_APL_H_
#define _AFLOW_APL_H_

#include <fstream>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xcomplex.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xparser_json.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow.h"
#include "aflow_aflowrc.h"
#include "flow/aflow_support_types.h"
#include "flow/aflow_xclasses.h"
#include "structure/aflow_xatom.h"
#include "structure/aflow_xstructure.h"

#define _APL_SHELL_TOL_ 0.1

// ME20190219 - Define the checksum algorithm used for APL hibernate files
#define APL_CHECKSUM_ALGO std::string("Fletcher32")

// ***************************************************************************

// Interface
namespace apl {
  void validateParametersAPL(aurostd::xoption&, const _aflags&, std::ofstream&, std::ostream& oss = std::cout);
  void validateParametersSupercellAPL(aurostd::xoption&);
  void validateParametersDispersionsAPL(aurostd::xoption&);
  void validateParametersDosAPL(aurostd::xoption&, const _aflags&, std::ofstream&, std::ostream& oss = std::cout);
  void validateParametersAAPL(aurostd::xoption&, const _aflags&, std::ofstream&, std::ostream& oss = std::cout);
  void validateParametersQHA(aurostd::xoption& aaplopts, const _aflags& aflags, std::ofstream& FileMESSAGE, std::ostream& oss); // AS20200709
  bool createAflowInPhonons(const _aflags&, const _kflags&, const _xflags&, _xinput&); // ME20190108
  void createAflowInPhononsAIMS(_aflags&, _kflags&, _xflags&, std::string&, _xinput&, std::ofstream&);
  bool filesExistPhonons(_xinput&);
  bool outfileFoundAnywherePhonons(std::vector<_xinput>&);
  bool outfileFoundEverywherePhonons(std::vector<_xinput>&, const std::string&, std::ofstream&, std::ostream&, bool = false);  // ME20191029
  bool readForcesFromDirectory(_xinput&);  // ME20200219
  void subtractZeroStateForces(std::vector<_xinput>&, bool);
  void subtractZeroStateForces(std::vector<_xinput>&, _xinput&);  // ME20190114

  bool APL_Get_AflowInName(std::string& AflowInName, const std::string& directory_LIB); // ME20210927

// ***************************************************************************
// BEGIN: Supplemental classes for APL, AAPL, and QHA
// ***************************************************************************

  class Supercell : public xStream {
  private:
    xstructure _inStructure;
    xstructure _inStructure_original;  // CO
    xstructure _inStructure_light;     // CO, does not include HEAVY symmetry stuff
      // deque<_atom> _inStructure_atoms_original; //CO
    xstructure _pcStructure;           // CO20180406 - for the path
    xstructure _scStructure;
    xstructure _scStructure_original;  // CO
    xstructure _scStructure_light;     // CO, does not include HEAVY symmetry stuff
      // deque<_atom> _scStructure_atoms_original; //CO
      // CO START
    bool _skew = false;                  // SYM::isLatticeSkewed(), same for pc and sc
    bool _derivative_structure = false;  // vs. simple expanded_lattice, derivative structure's lattice has LESS symmetry, so be careful ApplyAtom()'ing
    double _sym_eps = 0.0;             // same for pc and sc
      // CO END
    bool _isShellRestricted = false;
    int _maxShellID = -1;
    std::vector<double> _maxShellRadius;
    bool _isConstructed = false;
    bool _initialized = false;
    std::vector<std::vector<std::vector<aurostd::xvector<double>>>> phase_vectors;  // ME20200116

    void calculateWholeSymmetry(xstructure&, bool = true);
    [[nodiscard]] xstructure calculatePrimitiveStructure() const;
    bool getMaps(const xstructure&, const xstructure&, const xstructure&, std::vector<int>&, std::vector<int>&) const;  // ME20200117

  public:
    Supercell(std::ostream& oss = std::cout) : xStream(oss) {}
    Supercell(std::ofstream& mf, std::ostream& oss = std::cout) : xStream(mf, oss) {}
    Supercell(const xstructure&, std::ofstream&, const std::string& directory = "./", std::ostream& oss = std::cout); // CO20181226
    Supercell(const std::string&, std::ofstream&, const std::string& directory = "./", std::ostream& os = std::cout);  // ME20200112
    void clear();
    void initialize(const std::string&, std::ofstream&, std::ostream& oss = std::cout);
    void initialize(const std::string&);  // ME20200212
    void initialize(const xstructure&, std::ofstream&, bool = true, std::ostream& oss = std::cout);
    void initialize(const xstructure&, bool = true);  // ME20191225
    void clearSupercell();
    [[nodiscard]] bool isConstructed() const;
    void reset();
    aurostd::xvector<int> determineSupercellDimensions(const aurostd::xoption&);  // ME20191225
    void build(aurostd::xoption&, bool = true);  // ME20191225
    void build(const aurostd::xvector<int>&, bool = true);  // ME20191225
    void build(int, int, int, bool = true);
    void trimStructure(int, const aurostd::xvector<double>&, const aurostd::xvector<double>&, const aurostd::xvector<double>&, bool = true);
    bool projectToPrimitive();  // ME20200117
    void projectToOriginal();  // ME20200117
    aurostd::xvector<int> getSupercellDimensionsShell(uint, bool);
    void setupShellRestrictions(int);
      // ME20190715 BEGIN - added const to getter functions so they can be used with const Supercell &
    [[nodiscard]] bool isShellRestricted() const;
    [[nodiscard]] int getMaxShellID() const;
    [[nodiscard]] uint getNumberOfAtoms() const;
    [[nodiscard]] uint getNumberOfUniqueAtoms() const;
    [[nodiscard]] uint getNumberOfEquivalentAtomsOfType(int) const; // CO20190218
    [[nodiscard]] int getUniqueAtomID(int) const;
    [[nodiscard]] int getUniqueAtomID(int, int) const;
    [[nodiscard]] const _atom& getUniqueAtom(int) const;
    [[nodiscard]] std::string getUniqueAtomSymbol(int) const;
    [[nodiscard]] double getUniqueAtomMass(int) const;
    [[nodiscard]] double getAtomMass(int) const;
    [[nodiscard]] int getAtomNumber(int) const;
      // ME20190715 END
    [[nodiscard]] const xstructure& getSupercellStructure() const;
    [[nodiscard]] const xstructure& getSupercellStructureLight() const;
    [[nodiscard]] const xstructure& getPrimitiveStructure() const;
    [[nodiscard]] const xstructure& getInputStructure() const;
    [[nodiscard]] const xstructure& getInputStructureLight() const;
    [[nodiscard]] const xstructure& getOriginalStructure() const;  // ME20200117
    int atomGoesTo(const _sym_op&, int, int, bool = true); // CO20190218
    int atomComesFrom(const _sym_op&, int, int, bool = true); // CO20190218
    const _sym_op& getSymOpWhichMatchAtoms(int, int, int);
    void calculatePhaseVectors();  // ME20200117
    bool calcShellPhaseFactor(int, int, const aurostd::xvector<double>&, aurostd::xcomplex<double>&);
    bool calcShellPhaseFactor(int, int, const aurostd::xvector<double>&, aurostd::xcomplex<double>&, int&, aurostd::xvector<aurostd::xcomplex<double>>&, bool);  // ME20180828
    [[nodiscard]] int pc2scMap(int) const;
    [[nodiscard]] int sc2pcMap(int) const;
    void center(int);
      // CO START
    void center_original();
      // corey
    void getFullBasisAGROUP();  // ME20191218
    [[nodiscard]] bool fullBasisCalculatedAGROUP() const;  // ME20191218
    [[nodiscard]] const std::vector<std::vector<_sym_op>>& getAGROUP() const;
    [[nodiscard]] const std::vector<_sym_op>& getFGROUP() const;
    [[nodiscard]] const std::vector<_sym_op>& getAGROUP(int) const;
    [[nodiscard]] bool isDerivativeStructure() const;
    [[nodiscard]] double getEPS() const;
      // ME20190715 END
      // CO END
      //  **** BEGIN JJPR *****
    aurostd::xvector<int> scell_dim;
    std::vector<int> _pc2scMap;
    std::vector<int> _sc2pcMap;
      // **** END  JJPR *****
    std::string _directory;  // for the logger
  };

// ***************************************************************************

  struct _qpoint {
    aurostd::xvector<double> cpos;  // Cartesian position of the q-point
    aurostd::xvector<double> fpos;  // Fractional coordinates of the q-point
    int irredQpt;  // The irreducible q-point this q-point belongs to
    int ibzqpt;  // The index of the irreducible q-point in the _ibzqpt vector
    int symop;  // Symmetry operation to transform the irreducible q-point into this q-point
  };

  struct _kcell {
    aurostd::xmatrix<double> lattice;  // The reciprocal lattice vectors
    aurostd::xmatrix<double> rlattice;  // The real space lattice
    aurostd::xmatrix<double> c2f;  // Conversion matrix from Cartesian to fractional
    aurostd::xmatrix<double> f2c;  // Conversion matrix from fractional to Cartesian
    bool skewed;  // Is the lattice skewed?
    std::vector<_sym_op> pgroup;  // The point group operations of the reciprocal cell
  };

  class QMesh : public xStream {
  public:
    QMesh(std::ostream& oss = std::cout) : xStream(oss) {}
    QMesh(std::ofstream& mf, std::ostream& oss = std::cout) : xStream(mf, oss) {}
    QMesh(const aurostd::xvector<int>&, const xstructure&, std::ofstream&, bool include_inversions = true, bool gamma_centered = true, const std::string& directory = "./", std::ostream& oss = std::cout);
    QMesh(const std::vector<int>&, const xstructure&, std::ofstream&, bool include_inversions = true, bool gamma_centered = true, const std::string& directory = "./", std::ostream& oss = std::cout);

    void clear();
    void clear_tetrahedra();

    void initialize(const std::vector<int>&, const xstructure& xs, std::ofstream&, bool = true, bool = true, std::ostream& oss = std::cout);
    void initialize(const std::vector<int>&, const xstructure& xs, bool = true, bool = true);
    void initialize(const aurostd::xvector<int>&, const xstructure& xs, std::ofstream&, bool = true, bool = true, std::ostream& oss = std::cout);
    void initialize(const aurostd::xvector<int>&, const xstructure& xs, bool = true, bool = true);
    void initialize(std::ostream& oss);
    void initialize(std::ofstream& mf, std::ostream& oss);

    void makeIrreducible();
    void calculateLittleGroups();  // ME20200109
    void writeQpoints(std::string, bool = true);
    void writeIrredQpoints(std::string, bool = true);

    std::string _directory;

    [[nodiscard]] int getnIQPs() const;
    [[nodiscard]] int getnQPs() const;
    [[nodiscard]] int getGrid(int) const;
    [[nodiscard]] const aurostd::xvector<int>& getGrid() const;
    [[nodiscard]] const _qpoint& getIrredQPoint(int) const;
    [[nodiscard]] const _qpoint& getIrredQPoint(int, int, int) const;
    [[nodiscard]] std::vector<aurostd::xvector<double>> getIrredQPointsCPOS() const;
    [[nodiscard]] std::vector<aurostd::xvector<double>> getIrredQPointsFPOS() const;
    [[nodiscard]] int getIrredQPointIndex(int) const;
    [[nodiscard]] int getIrredQPointIndex(int, int, int) const;
    [[nodiscard]] const _qpoint& getQPoint(int) const;
    [[nodiscard]] const _qpoint& getQPoint(int, int, int) const;
    [[nodiscard]] const _qpoint& getQPoint(const aurostd::xvector<double>&) const;  // ME20190813
    [[nodiscard]] int getQPointIndex(aurostd::xvector<double>) const;  // ME20190813
    [[nodiscard]] int getQPointIndex(int, int, int) const;
    [[nodiscard]] std::vector<aurostd::xvector<double>> getQPointsCPOS() const;
    [[nodiscard]] std::vector<aurostd::xvector<double>> getQPointsFPOS() const;
    [[nodiscard]] int getIbzqpt(int) const;
    [[nodiscard]] int getIbzqpt(int, int, int) const;
    [[nodiscard]] const std::vector<int>& getIbzqpts() const;
    [[nodiscard]] const std::vector<_qpoint>& getPoints() const;
    [[nodiscard]] const _kcell& getReciprocalCell() const;
    [[nodiscard]] bool isShifted() const;  // ME20190813
    [[nodiscard]] const aurostd::xvector<double>& getShift() const;
    [[nodiscard]] const std::vector<int>& getWeights() const;
    [[nodiscard]] bool initialized() const;
    [[nodiscard]] bool isReduced() const;
    [[nodiscard]] bool isGammaCentered() const;
    [[nodiscard]] bool littleGroupsCalculated() const;  // ME20200109
    [[nodiscard]] const std::vector<int>& getLittleGroup(int) const;  // ME20200109

      // Tetrahedron method
    void generateTetrahedra();
    void makeIrreducibleTetrahedra();

    [[nodiscard]] const std::vector<std::vector<int>>& getTetrahedra() const;
    [[nodiscard]] const std::vector<int>& getTetrahedron(int) const;
    [[nodiscard]] const std::vector<int>& getIrredTetrahedron(int) const;
    [[nodiscard]] int getTetrahedronCorner(int, int) const;
    [[nodiscard]] std::vector<std::vector<int>> getIrreducibleTetrahedra() const;
    [[nodiscard]] std::vector<std::vector<int>> getIrreducibleTetrahedraIbzqpt() const;
    [[nodiscard]] int getnTetrahedra() const;
    [[nodiscard]] int getnIrredTetrahedra() const;
    [[nodiscard]] double getVolumePerTetrahedron() const;
    [[nodiscard]] const std::vector<int>& getWeightsTetrahedra() const;
    [[nodiscard]] int getWeightTetrahedron(int) const;
    [[nodiscard]] bool isReducedTetrahedra() const;

  private:
    std::vector<int> _ibzqpts;  // The indices of the irreducible q-points
    bool _initialized = false;  // Indicates whether the QMesh object has been intialized
    bool _isGammaCentered = false;  // Indicates whether the includes the Gamma point
    std::vector<std::vector<int>> _littleGroups;  // The little groups of the irreducible q-points
    bool _littleGroupsCalculated = false;  // Indicates whether the little groups have been calculated
    int _nIQPs = 0;  // The number of irreducible q-points
    int _nQPs = 0;  // The number of q-points
    aurostd::xvector<int> _qptGrid;  // The dimensions of the q-point mesh
    std::vector<std::vector<std::vector<int>>> _qptMap;  // Maps a q-point triplet to a q-point index
    std::vector<_qpoint> _qpoints;  // The q-points of the mesh
    _kcell _recCell;  // The reciprocal cell
    bool _reduced = false;  // Indicates whether the q-point mesh has been reduced
    bool _shifted = false;  // Indicates whether the q-point mesh has been shifted
    aurostd::xvector<double> _shift;  // The shift vector of the mesh
    std::vector<int> _weights;  // The weights of each irreducible q-point

    void setGrid(const aurostd::xvector<int>&);
    void setupReciprocalCell(xstructure, bool);
    void generateGridPoints(bool);
    void shiftMesh(const aurostd::xvector<double>&);
    void moveToBZ(aurostd::xvector<double>&) const;

      // Tetrahedron method
    std::vector<std::vector<int>> _tetrahedra;  // The corners of the tetrahedra
    std::vector<int> _irredTetrahedra;  // List of irreducible tetrahedra
    bool _reducedTetrahedra = false; // Indicates whether the tetrahedra are reduced
    int _nTetra = 0;  // The number of tetrahedra
    int _nIrredTetra = 0;  // The number of irreducible tetrahedra - ME20190625
    double _volumePerTetrahedron = 0.0;  // The relative volume of each tetrahedron
    std::vector<int> _weightsTetrahedra;  // The weights of each irreducible tetrahedron

    std::vector<std::vector<aurostd::xvector<int>>> initializeTetrahedra();
    void findMostCompactTetrahedra(std::vector<std::vector<aurostd::xvector<int>>>&) const;
    void generateAllTetrahedra(const std::vector<std::vector<aurostd::xvector<int>>>&);
  };

// ***************************************************************************
// END: Supplemental classes for APL, AAPL, and QHA
// ***************************************************************************

// ***************************************************************************
// BEGIN: Automatic Phonon Library (APL)
// ***************************************************************************

// ***************************************************************************

#define _AFLOW_APL_BORN_EPSILON_RUNNAME_ std::string("LRBE")  // ME20190108
#define _AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_ std::string(ARUN_DIRECTORY_PREFIX + "APL_" + _AFLOW_APL_BORN_EPSILON_RUNNAME_) // ME20190108
#define _AFLOW_APL_DFPT_RUNNAME_ std::string("DFPT")  // ME20200213
#define _AFLOW_APL_DFPT_DIRECTORY_NAME_ std::string(ARUN_DIRECTORY_PREFIX + "APL_" + _AFLOW_APL_DFPT_RUNNAME_) // ME20200213

  class ForceConstantCalculator : public xStream {
  protected:
    Supercell* _supercell = nullptr;

  private:
    bool _sc_set = false;
    bool _initialized = false;

    std::vector<_xinput> xInputs;
    std::string _method;

      // Calculate forces at no distortion - since for some structure
      // (not well relaxed, or with other problems) these forces have to be
      // known and substracted from the calculated forces with distortion
    bool _calculateZeroStateForces = false;

      // For each atom of supercell, there is a full force field
    std::vector<std::vector<aurostd::xmatrix<double>>> _forceConstantMatrices;
      // Stuff for polar materials
    bool _isPolarMaterial = false;
      // For each atom there is a matrix 3x3 of Born effective charge
    std::vector<aurostd::xmatrix<double>> _bornEffectiveChargeTensor;
      // Dielectric tensor
    aurostd::xmatrix<double> _dielectricTensor;

      // Force constants
    bool runVASPCalculationsDM(_xinput&, _aflags&, _kflags&, _xflags&, std::string&);
    bool runVASPCalculationsLR(_xinput&, _aflags&, _kflags&, _xflags&, std::string&);
    bool calculateForceConstants(); // ME20200211
    bool calculateForceConstantsDM();
    bool readForceConstantsFromVasprun(_xinput&);
    void symmetrizeForceConstantMatrices();
    void correctSumRules();

      // Direct method
    bool AUTO_GENERATE_PLUS_MINUS = true;
    bool USER_GENERATE_PLUS_MINUS = false;
    bool GENERATE_ONLY_XYZ = false;
    bool DISTORTION_SYMMETRIZE = true; // CO20190108
    double DISTORTION_MAGNITUDE = 0.0;
    bool DISTORTION_INEQUIVONLY = true; // CO20190108
      // For each inequivalent atom, there is a set of unique distortions
    std::vector<std::vector<aurostd::xvector<double>>> _uniqueDistortions;
      // For each inequivalent atom and unique distortion, there is a field
      // of forces (for each atom of the supercell)
    std::vector<std::vector<std::vector<aurostd::xvector<double>>>> _uniqueForces;
    std::vector<std::vector<bool>> vvgenerate_plus_minus;  // ME20191029

    void estimateUniqueDistortions(const xstructure&, std::vector<std::vector<aurostd::xvector<double>>>&);
    void testDistortion(const aurostd::xvector<double>&,
                        const std::vector<_sym_op>&,
                        std::vector<aurostd::xvector<double>>&,
                        std::vector<aurostd::xvector<double>>&,
                        bool integrate_equivalent_distortions = true);  // CO20190114
    bool needMinus(uint atom_index, uint distortion_index, bool inequiv_only = true);  // CO //CO20190218
    bool calculateForceFields();  // ME20190412  //ME20191029
    void completeForceFields();
    void projectToCartesianDirections();
    void buildForceConstantMatrices();

      // Born charges + dielectric tensor
    bool runVASPCalculationsBE(_xinput&, _aflags&, _kflags&, _xflags&, std::string&, uint);
    bool calculateBornChargesDielectricTensor(const _xinput&);  // ME20191029
    void readBornEffectiveChargesFromAIMSOUT();
    void readBornEffectiveChargesFromOUTCAR(const _xinput&);  // ME20190113
    void symmetrizeBornEffectiveChargeTensors();
    void readDielectricTensorFromAIMSOUT();
    void readDielectricTensorFromOUTCAR(const _xinput&);  // ME20190113

  public:
    ForceConstantCalculator(std::ostream& oss = std::cout) : xStream(oss), _directory("./") {}
    ForceConstantCalculator(Supercell& sc, std::ostream& oss = std::cout) : xStream(oss), _supercell(&sc), _sc_set(true), _directory(sc._directory) {}
    ForceConstantCalculator(Supercell& sc, std::ofstream& mf, std::ostream& oss = std::cout) : xStream(mf, oss), _supercell(&sc), _sc_set(true), _directory(sc._directory) {}
    ForceConstantCalculator(Supercell&, const aurostd::xoption&, std::ofstream&, std::ostream& os = std::cout);
    void clear();
    void clear(Supercell&);
    void initialize(const aurostd::xoption&, std::ofstream&, std::ostream& oss = std::cout);
    void initialize(const aurostd::xoption&);
    void initialize(const aurostd::xvector<int>&, const xstructure& xs, bool = true, bool = true);

    bool runVASPCalculations(_xinput&, _aflags&, _kflags&, _xflags&, std::string&);

    bool run();  // ME20191029
    void hibernate();

    [[nodiscard]] const std::vector<std::vector<aurostd::xmatrix<double>>>& getForceConstants() const;
    [[nodiscard]] const std::vector<aurostd::xmatrix<double>>& getBornEffectiveChargeTensor() const;
    [[nodiscard]] const aurostd::xmatrix<double>& getDielectricTensor() const;
    [[nodiscard]] bool isPolarMaterial() const;

    std::string _directory;

    void writeHarmonicIFCs(const std::string&);
    void writeBornChargesDielectricTensor(const std::string&);
    void writeDYNMAT(const std::string&);
    void saveState(const std::string&);  // ME20200112
    void readFromStateFile(const std::string&);  // ME20200112
  };

// ***************************************************************************

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

// ***************************************************************************

  class AnharmonicIFCs;  // Forward declaration
  class PhononCalculator : public xStream {
  private:
      // USER PARAMETERS
    std::string _directory;  // for loggers
    int _ncpus = 0;

    QMesh _qm;
    Supercell _supercell;

      // harmonic IFCs
    std::vector<std::vector<aurostd::xmatrix<double>>> _forceConstantMatrices;

      // Stuff for polar materials
    bool _isPolarMaterial = false;
      // For each atom there is a matrix 3x3 of Born effective charge
    std::vector<aurostd::xmatrix<double>> _bornEffectiveChargeTensor;
      // Dielectric tensor
    aurostd::xmatrix<double> _dielectricTensor;
      // Precomputed values used in non-analytical term (Gonze)
    aurostd::xmatrix<double> _inverseDielectricTensor;
    double _recsqrtDielectricTensorDeterminant = 0.0;
      // Precomputed Ewald sum at Gamma point
    bool _isGammaEwaldPrecomputed = false;
    std::vector<aurostd::xmatrix<aurostd::xcomplex<double>>> _gammaEwaldCorr;
      // Anharmonic IFCs
    std::vector<std::vector<std::vector<double>>> anharmonicIFCs;
    std::vector<std::vector<std::vector<int>>> clusters;

    void readHarmonicIFCs(const std::string&);
    void readBornChargesDielectricTensor(const std::string&);

    aurostd::xmatrix<aurostd::xcomplex<double>> getNonanalyticalTermWang(const aurostd::xvector<double>&);
    aurostd::xmatrix<aurostd::xcomplex<double>> getNonanalyticalTermWang(const aurostd::xvector<double>&, std::vector<aurostd::xmatrix<aurostd::xcomplex<double>>>&, bool = true);  // ME20180829
    aurostd::xmatrix<aurostd::xcomplex<double>> getNonanalyticalTermGonze(const aurostd::xvector<double>);
    aurostd::xmatrix<aurostd::xcomplex<double>> getEwaldSumDipoleDipoleContribution(const aurostd::xvector<double>, bool = true);

    void calculateGroupVelocitiesThread(int, std::vector<std::vector<double>>&, std::vector<aurostd::xmatrix<aurostd::xcomplex<double>>>&, std::vector<std::vector<aurostd::xvector<double>>>&);

  public:
    PhononCalculator(std::ostream& oss = std::cout) : xStream(oss), _qm(oss), _supercell(oss), _ncpus(1), _directory("./") {}
    PhononCalculator(std::ofstream& mf, std::ostream& oss = std::cout) : xStream(mf, oss), _qm(mf, oss), _supercell(mf, oss), _ncpus(1), _directory("./") {}
    void clear();

    std::string _system;  // ME20190614 - for VASP-style output files

      // Getter functions
    Supercell& getSupercell();
    QMesh& getQMesh();
    [[nodiscard]] const xstructure& getInputCellStructure() const;
    [[nodiscard]] const xstructure& getSuperCellStructure() const;
    [[nodiscard]] uint getNumberOfBranches() const;
    [[nodiscard]] std::string getDirectory() const;
    [[nodiscard]] int getNCPUs() const;
    [[nodiscard]] bool isPolarMaterial() const;  // ME20200206
    [[nodiscard]] const std::vector<std::vector<aurostd::xmatrix<double>>>& getHarmonicForceConstants() const;
    [[nodiscard]] const std::vector<std::vector<double>>& getAnharmonicForceConstants(int) const;
    [[nodiscard]] const std::vector<std::vector<int>>& getClusters(int) const;

      // Set functions
    void setDirectory(const std::string&);
    void setNCPUs(const _kflags&);
    void setPolarMaterial(bool);

      // Initializers
    void initialize_qmesh(const std::vector<int>&, bool = true, bool = true);
    void initialize_qmesh(const aurostd::xvector<int>&, bool = true, bool = true);
    void initialize_supercell(const xstructure&, bool verbose = true);// AS20200908
    void initialize_supercell(const std::string&);

      // IFCs
    void setHarmonicForceConstants(const ForceConstantCalculator&);
    void setHarmonicForceConstants(const std::vector<std::vector<aurostd::xmatrix<double>>>&);// AS20201204
    void setHarmonicForceConstants(const std::vector<std::vector<aurostd::xmatrix<double>>>&, const std::vector<aurostd::xmatrix<double>>&, const aurostd::xmatrix<double>&, bool isPolar = true);// AS20201204
    void awake();
    void setAnharmonicForceConstants(const AnharmonicIFCs&);
    void readAnharmonicIFCs(std::string);

      // Dynamical Matrix/Frequencies
    aurostd::xvector<double> getEigenvalues(const aurostd::xvector<double>&, const aurostd::xvector<double>&, aurostd::xmatrix<aurostd::xcomplex<double>>&, std::vector<aurostd::xmatrix<aurostd::xcomplex<double>>>&, bool = true); // ME20180827
    aurostd::xmatrix<aurostd::xcomplex<double>> getDynamicalMatrix(const aurostd::xvector<double>&);
    aurostd::xmatrix<aurostd::xcomplex<double>> getDynamicalMatrix(const aurostd::xvector<double>&, const aurostd::xvector<double>&); // ME20200206
    aurostd::xmatrix<aurostd::xcomplex<double>> getDynamicalMatrix(const aurostd::xvector<double>&, const aurostd::xvector<double>&, std::vector<aurostd::xmatrix<aurostd::xcomplex<double>>>&, bool = true); // ME20180827
    aurostd::xvector<double> getFrequency(const aurostd::xvector<double>&, const IPCFreqFlags&); // ME20180827
    aurostd::xvector<double> getFrequency(const aurostd::xvector<double>&, const aurostd::xvector<double>&, const IPCFreqFlags&); // ME20200206
    aurostd::xvector<double> getFrequency(const aurostd::xvector<double>&, const IPCFreqFlags&, aurostd::xmatrix<aurostd::xcomplex<double>>&); // ME20190624
    aurostd::xvector<double> getFrequency(const aurostd::xvector<double>&, const aurostd::xvector<double>&, const IPCFreqFlags&, aurostd::xmatrix<aurostd::xcomplex<double>>&); // ME20200206
    aurostd::xvector<double> getFrequency(const aurostd::xvector<double>&, const IPCFreqFlags&, aurostd::xmatrix<aurostd::xcomplex<double>>&, std::vector<aurostd::xvector<double>>&, bool = true); // ME20180827
    aurostd::xvector<double> getFrequency(const aurostd::xvector<double>&, const aurostd::xvector<double>&, const IPCFreqFlags&, aurostd::xmatrix<aurostd::xcomplex<double>>&, std::vector<aurostd::xvector<double>>&, bool = true); // ME20200206
    double getFrequencyConversionFactor(IPCFreqFlags, IPCFreqFlags);

    // Group velocities
    std::vector<std::vector<aurostd::xvector<double>>> calculateGroupVelocitiesOnMesh();
    std::vector<std::vector<aurostd::xvector<double>>> calculateGroupVelocitiesOnMesh(std::vector<std::vector<double>>&);
    std::vector<std::vector<aurostd::xvector<double>>> calculateGroupVelocitiesOnMesh(std::vector<std::vector<double>>&, std::vector<aurostd::xmatrix<aurostd::xcomplex<double>>>&);
    void writeGroupVelocitiesToFile(const std::string&, const std::vector<std::vector<aurostd::xvector<double>>>&);
    void writeGroupVelocitiesToFile(const std::string&, const std::vector<std::vector<aurostd::xvector<double>>>&, const std::vector<std::vector<double>>&, const std::string& unit = "THz");
  };

  // ***************************************************************************

  class PathBuilder {
  public:
    enum StoreEnumType { RECIPROCAL_LATTICE, CARTESIAN_LATTICE };
    enum ModeEnumType { SINGLE_POINT_MODE, COUPLE_POINT_MODE };

  private:
    std::vector<aurostd::xvector<double>> _path;
    std::vector<aurostd::xvector<double>> _points;
    std::vector<std::string> _labels;
    aurostd::xmatrix<double> reciprocalLattice;
    aurostd::xmatrix<double> cartesianLattice;
    int _pointsVectorDimension = 0;
    int _pointsVectorStartingIndex = 0;
    uint _nPointsPerSubPath = 0;
    ModeEnumType _mode = SINGLE_POINT_MODE;
    StoreEnumType _store = CARTESIAN_LATTICE;

    void buildPath();

  public:
    PathBuilder() = default;
    PathBuilder(ModeEnumType mode) : _mode(mode) {}
    void clear();
    void addPoint(const std::string& l, int dim, ...);
    void addPoint(const std::string&, const aurostd::xvector<double>&);
    void transform(const aurostd::xmatrix<double>&);
    void pointsAreDefinedFor(const xstructure&, StoreEnumType);
    void transformPointsFor(const xstructure&, StoreEnumType);
    void defineCustomPoints(const std::string&, const std::string&, const Supercell&, bool CARESTIAN_COORDS = false);
    void takeAflowElectronicPath(const std::string&, const Supercell&);//, const xstructure&, const xstructure&);
    void setMode(ModeEnumType);
    void setStore(StoreEnumType);
    [[nodiscard]] const StoreEnumType& getStore() const;  // ME20190614
    void setDensity(int);
    [[nodiscard]] int getDensity() const;
    uint getPathSize();
    uint getPointSize();
    aurostd::xvector<double> getPoint(uint);
    uint getPointIndexOnPath(uint);
    std::string getPointLabel(uint);
    std::vector<aurostd::xvector<double>> getPath();
    std::vector<aurostd::xvector<double>> getPath(ModeEnumType, const std::string&);
    double getPathLength();
    double getPathLength(uint);
    xKPOINTS createKPOINTS(const Supercell&);  // ME20190614
  };

// ***************************************************************************

  class PhononDispersionCalculator {
  private:
    PhononCalculator* _pc = nullptr;
    bool _pc_set = false;
    PathBuilder _pb;
    std::vector<aurostd::xvector<double>> _qpoints;
    std::vector<aurostd::xvector<double>> _freqs;
    IPCFreqFlags _frequencyFormat = apl::NONE;
    double _temperature = 0.0; // ME20190614
    void calculateInOneThread(int);
    bool isExactQPoint(const aurostd::xvector<double>&, const aurostd::xmatrix<double>&);
    std::string _system;

  public:
    PhononDispersionCalculator() = default;
    PhononDispersionCalculator(PhononCalculator& pc) : _pc(&pc), _pc_set(true), _system(pc._system) {}
    void clear();
    void clear(PhononCalculator&);
    void initPathCoords(const std::string&, const std::string&, int, bool = false);  // CO20180406
    void initPathLattice(const std::string&, int);
    void setPath(const std::string&);
    void calc(const IPCFreqFlags);
    void writePDIS(const std::string&);
    std::vector<aurostd::xvector<double>> get_qpoints() { return _qpoints; } //[PN]
    // ME20190614 START
    xEIGENVAL createEIGENVAL();
    void writePHEIGENVAL(const std::string&);
    void writePHKPOINTS(const std::string&);
    xKPOINTS getPHKPOINTS(); // AS20201110
    // ME20190614 STOP
  };

  // ***************************************************************************

  class DOSCalculator {
  protected:
    PhononCalculator* _pc = nullptr;
    bool _pc_set = false;
    std::string _bzmethod;  // ME20190423
    std::vector<aurostd::xvector<double>> _qpoints;
    std::vector<aurostd::xvector<double>> _freqs;
    double _minFreq = AUROSTD_NAN;
    double _maxFreq = AUROSTD_NAN;
    double _stepDOS = 0.0;
    double _halfStepDOS = 0.0;
    std::vector<double> _bins;
    std::vector<double> _dos;
    std::vector<double> _idos; // ME20190614
    std::vector<aurostd::xmatrix<aurostd::xcomplex<double>>> _eigen; // ME20190624 - eigenvectors for projected DOS
    std::vector<std::vector<std::vector<double>>> _projectedDOS; // ME20190614 - projectedDOS.at(atom).at(direction).at(frequency)
    std::vector<aurostd::xvector<double>> _projections; // ME20190626 - the projection directions for the DOS in Cartesian coordinates
    double _temperature = 0.0; // ME20190614
    void calculateInOneThread(int);
    void calculateFrequencies();
    void smearWithGaussian(std::vector<double>&, std::vector<double>&, double, double); // ME20190614
    void calcDosRS();
    void calcDosLT();

  public:
    DOSCalculator() = default;
    DOSCalculator(PhononCalculator&, const aurostd::xoption&); // AS20201203 input parameters are passed via xoption
    void clear();
    void clear(PhononCalculator&);
    void initialize(const aurostd::xoption&); // AS20201203 input parameters are passed via xoption
    void calc(int, bool VERBOSE = true);
    void calc(int, double, bool VERBOSE = true);
    void calc(int, double, double, double, bool VERBOSE = true); // ME20200203
    void writePDOS(const std::string&);
    void writePDOS(std::string, std::string); //[PN]
    [[nodiscard]] xDOSCAR createDOSCAR() const; // ME20190614
    void writePHDOSCAR(const std::string&); // ME20190614
    // Interface IDOSCalculator
    [[nodiscard]] const std::vector<double>& getBins() const; // ME20200108 - added const
    [[nodiscard]] const std::vector<double>& getDOS() const; // ME20200108 - added const
    [[nodiscard]] const std::vector<double>& getIDOS() const; // ME20200210
    [[nodiscard]] const std::vector<aurostd::xvector<double>>& getFreqs() const; // AS20200312
    [[nodiscard]] bool hasImaginaryFrequencies() const; // ME20200108 - added const
    [[nodiscard]] double getMinFreq() const; // ME20210927
    [[nodiscard]] double getMaxFreq() const; // ME20210927
    [[nodiscard]] const xstructure& getInputStructure() const; // ME20210927
    [[nodiscard]] uint getNumberOfBranches() const; // ME20210927
    std::string _system;
  };

  // ***************************************************************************

  enum ThermalPropertiesUnits { eV, meV, ueV, eVK, meVK, ueVK, kB };

  class ThermalPropertiesCalculator : public xStream {
  private:
    std::vector<double> _freqs_0K;
    std::vector<double> _dos_0K;
    std::string system;

    double getStepDOS(const std::vector<double>&);
    double getScalingFactor(const ThermalPropertiesUnits&);

  public:
    ThermalPropertiesCalculator(std::ostream& oss = std::cout) : xStream(oss) {}
    ThermalPropertiesCalculator(std::ofstream& mf, std::ostream& oss = std::cout) : xStream(mf, oss) {}
    ThermalPropertiesCalculator(const DOSCalculator&, std::ofstream&, const std::string& directory = "./", std::ostream& os = std::cout);
    ThermalPropertiesCalculator(const xDOSCAR&, std::ofstream&, const std::string& directory = "./", std::ostream& os = std::cout);
    void clear();

    uint natoms = 0;
    std::vector<double> temperatures;
    std::vector<double> Cv;
    std::vector<double> Fvib;
    std::vector<double> Svib;
    std::vector<double> U;
    double U0 = 0.0;
    std::string _directory;

    void initialize(const xDOSCAR&, std::ofstream&, std::ostream& oss = std::cout);
    void initialize(const std::vector<double>&, const std::vector<double>&, std::ofstream&, const std::string& system = "", std::ostream& oss = std::cout);
    void initialize(const std::vector<double>&, const std::vector<double>&, const std::string& system = "");
    void calculateThermalProperties(double, double, double);
    void addPoint(double, const xDOSCAR&);
    void addPoint(double, const std::vector<double>&, const std::vector<double>&);

    double getZeroPointEnergy();
    double getInternalEnergy(double, ThermalPropertiesUnits = apl::meV);
    double getInternalEnergy(double, const std::vector<double>&, const std::vector<double>&, ThermalPropertiesUnits = apl::meV);
    double getVibrationalFreeEnergy(double, ThermalPropertiesUnits = apl::meV);
    double getVibrationalFreeEnergy(double, const std::vector<double>&, const std::vector<double>&, ThermalPropertiesUnits = apl::meV);
    double getVibrationalEntropy(double, ThermalPropertiesUnits = apl::meV);
    double getVibrationalEntropy(double, const std::vector<double>&, const std::vector<double>&, ThermalPropertiesUnits = apl::kB);
    double getVibrationalEntropy(double, double, double, ThermalPropertiesUnits = apl::kB);
    double getIsochoricSpecificHeat(double, ThermalPropertiesUnits = apl::kB);
    double getIsochoricSpecificHeat(double, const std::vector<double>&, const std::vector<double>&, ThermalPropertiesUnits = apl::kB);

    void writePropertiesToFile(std::string, filetype ft = txt_ft);
    void addToAPLOut(std::stringstream&);
    std::string getPropertiesFileString(filetype ft = txt_ft);
    aurostd::JSON::object getPropertiesJSON();
  };

// ***************************************************************************

  class AtomicDisplacements {
  protected:
    PhononCalculator* _pc = nullptr;
    bool _pc_set = false;

  private:
    std::vector<std::vector<std::vector<aurostd::xvector<aurostd::xcomplex<double>>>>> _eigenvectors;
    std::vector<std::vector<double>> _frequencies;
    std::vector<std::vector<aurostd::xmatrix<aurostd::xcomplex<double>>>> _displacement_matrices;
    std::vector<std::vector<std::vector<aurostd::xvector<aurostd::xcomplex<double>>>>> _displacement_modes;
    std::vector<_qpoint> _qpoints;
    std::vector<double> _temperatures;

    void calculateEigenvectors();
    void calculateEigenvectorsInThread(int);
    void calculateMeanSquareDisplacementMatrices();
    void calculateModeDisplacements();
    double getOccupationNumber(double, double);

  public:
    AtomicDisplacements() = default;
    AtomicDisplacements(PhononCalculator& pc) : _pc(&pc), _pc_set(true) {}
    void clear();
    void clear(PhononCalculator&);

    void calculateMeanSquareDisplacements(double, double, double);
    void calculateModeDisplacements(const std::vector<aurostd::xvector<double>>& qpts, bool = true);

    [[nodiscard]] const std::vector<double>& getTemperatures() const;
    [[nodiscard]] const std::vector<std::vector<aurostd::xmatrix<aurostd::xcomplex<double>>>>& getDisplacementMatrices() const;
    [[nodiscard]] std::vector<std::vector<aurostd::xvector<double>>> getDisplacementVectors() const;
    [[nodiscard]] const std::vector<std::vector<std::vector<aurostd::xvector<aurostd::xcomplex<double>>>>>& getModeDisplacements() const;

    std::vector<std::vector<std::vector<double>>> createDisplacementsXcrysden(const Supercell&, double, int, int, int);
    void getOrientedDisplacementsVsim(xstructure&, std::vector<std::vector<std::vector<aurostd::xvector<aurostd::xcomplex<double>>>>>&, double);

    void writeMeanSquareDisplacementsToFile(std::string);
    void writeSceneFileXcrysden(std::string, const xstructure&, const std::vector<std::vector<std::vector<double>>>&, int);
    void writeSceneFileVsim(std::string, const xstructure&, const std::vector<std::vector<std::vector<aurostd::xvector<aurostd::xcomplex<double>>>>>&);
  };

  void createAtomicDisplacementSceneFile(const aurostd::xoption& vpflow, std::ostream& oss = std::cout);
  void createAtomicDisplacementSceneFile(const aurostd::xoption& vpflow, std::ofstream&, std::ostream& oss = std::cout);

  // ***************************************************************************
  // END: Automatic Phonon Library (APL)
  // ***************************************************************************

  // ***************************************************************************
  // BEGIN ME: Automatic Anharmonic Phonon Library (AAPL)
  // ***************************************************************************

  // _cluster holds a single cluster
  struct _cluster {
    std::vector<int> atoms; // List of atoms inside the cluster
    int fgroup; // Index pointing to the factor group that transforms the cluster into another cluster
    int permutation; // Index pointing to the permutation that transforms the cluster into another cluster
  };

  // _ineq_distortions contains a list of inequivalent distortion and its equivalent
  // distortions for a given set of atoms
  struct _ineq_distortions {
    std::vector<int> atoms;  // A list of atoms involved in the distortions
    std::vector<int> clusters;  // A list of cluster sets that use these distortions for their force constant calculations
    std::vector<std::vector<std::vector<int>>> distortions; // Map of distortions. The distortions vectors need to be defined elsewhere.
    std::vector<std::vector<int>> rotations;  // The factor group that holds the rotation to transform the distortions
    std::vector<std::vector<std::vector<int>>> transformation_maps;  // A map containing the transformation of the atoms for each distortion
  };

  // _linearCombinations is a structure to store linear combinations.
  struct _linearCombinations {
    std::vector<std::vector<int>> indices;  // Cartesian indices of each linear combination
    std::vector<std::vector<double>> coefficients;  // Coefficients of each linear combination
    std::vector<int> independent;  // The linearly independent values
    std::vector<int> dependent;  // The linearly dependent values
    // indep2depMap maps the independent coefficients to the coefficients that
    // depend on them. This is used for the IFC correction method.
    std::vector<std::vector<int>> indep2depMap;
  };

  class ClusterSet : public xStream {
    // See aflow_aapl_cluster.cpp for detailed descriptions of the functions
  public:
    ClusterSet(std::ostream& oss = std::cout) : xStream(oss), directory("./") {}
    ClusterSet(std::ofstream& mf, std::ostream& oss = std::cout) : xStream(mf, oss), directory("./") {}
    ClusterSet(const Supercell&, int, int, double, const std::string&, std::ofstream&, std::ostream& oss = std::cout);  // Constructor
    ClusterSet(const std::string&, const Supercell&, int, int, double, const std::string&, std::ofstream&, std::ostream& oss = std::cout);  // From file
    void clear();
    void initialize(const Supercell&, int, int, double, std::ofstream&, std::ostream& oss = std::cout);
    void initialize(const Supercell&, int, int, double);
    void initialize(std::ostream& oss);
    void initialize(std::ofstream& mf, std::ostream& oss);
    void readClusterSetFromFile(const std::string&);

    std::vector<_cluster> clusters;
    std::vector<std::vector<int>> coordination_shells; // Contains all coordinate shells. Central atoms is index 0.
    double cutoff = 0.0; // Cutoff radius in Angstroms
    std::string directory; // Directory for logging
    std::vector<aurostd::xvector<double>> distortion_vectors; // List of distortion vectors
    std::vector<_ineq_distortions> higher_order_ineq_distortions; // ME20190531 - for 3rd derivatives of higher order processes
    std::vector<std::vector<int>> ineq_clusters; // Clusters rearranged into sets of equivalent clusters.  //ME20190520
    std::vector<_ineq_distortions> ineq_distortions; // List of inequivalent distortions
    std::vector<_linearCombinations> linear_combinations; // List of linear combinations of the IFCs
    int nifcs = 0; // Number of force constants for each set of atoms.
    int order = 0; // Order of the cluster, i.e. the order of the force constant to be calculated.
    xstructure pcell; // Structure of the primitive cell.
    std::vector<int> pc2scMap; // Atom map from the primitive cell to the supercell.
    std::vector<std::vector<int>> permutations; // List of possible permutations for the cluster
    xstructure scell; // Structure of the supercell.
    std::vector<int> sc2pcMap; // Atom map from the supercell to the primitive cell.
    aurostd::xvector<int> sc_dim; // Dimensions of the supercell.
    std::vector<std::vector<int>> symmetry_map; // Symmetry atom map for the atoms in the clusters

    [[nodiscard]] const _cluster& getCluster(const int& i) const; // ME20190520
    void build();
    void buildDistortions();
    void writeClusterSetToFile(const std::string&);

  private:
    double getMaxRad(const xstructure&, int);
    void buildShells();
    std::vector<_cluster> buildClusters();
    std::vector<std::vector<int>> getSymmetryMap();
    std::vector<std::vector<int>> getPermutations(int);

      // Clusters
    void getInequivalentClusters(std::vector<_cluster>&, std::vector<std::vector<int>>&);
    [[nodiscard]] int getNumUniqueAtoms(const std::vector<int>&) const;
    std::vector<int> getComposition(const std::vector<int>&);
    [[nodiscard]] bool sameComposition(const std::vector<int>&, const std::vector<int>&) const;
    int equivalenceCluster(const std::vector<int>&, const std::vector<int>&, const std::vector<std::vector<int>>&, const std::vector<std::vector<int>>&);
    std::vector<int> translateToPcell(const std::vector<int>&, int);
    int comparePermutations(const std::vector<int>&, const std::vector<int>&);
    bool atomsMatch(const std::vector<int>&, const std::vector<int>&, const std::vector<int>&, const int&);
    void getSymOp(_cluster&, const std::vector<int>&);

    // Distortions
    std::vector<aurostd::xvector<double>> getCartesianDistortionVectors();
    std::vector<_ineq_distortions> initializeIneqDists();
    [[nodiscard]] int sameDistortions(const _cluster&, const std::vector<_ineq_distortions>&) const;
    std::vector<std::vector<int>> getTestDistortions(const std::vector<int>&);
    void getInequivalentDistortions(const std::vector<std::vector<int>>&, _ineq_distortions&);
    void appendDistortion(_ineq_distortions&, std::vector<int>, const int& eq = -1, const int& fg = -1);
    bool allZeroDistortions(const std::vector<int>&, const std::vector<int>&);
    bool allAtomsMatch(const int&, const std::vector<int>&);
    int equivalenceDistortions(const aurostd::xmatrix<double>&, const std::vector<int>&, const std::vector<std::vector<std::vector<int>>>&, const std::vector<int>&);
    std::vector<int> getTransformationMap(const int&, const int&);
    std::vector<_ineq_distortions> getHigherOrderDistortions();

    // Linear Combinations
    std::vector<_linearCombinations> getLinearCombinations();
    std::vector<std::vector<int>> getInvariantSymOps(const _cluster&);
    std::vector<std::vector<double>> buildCoefficientMatrix(const std::vector<std::vector<int>>&);
    std::vector<std::vector<double>> getRREF(std::vector<std::vector<double>>);

    // File I/O
    std::string writeParameters();
    std::string writeInequivalentClusters();
    [[nodiscard]] std::string writeClusters(const std::vector<_cluster>&) const;
    std::string writeLinearCombinations(const _linearCombinations&);
    std::string writeInequivalentDistortions();
    std::string writeIneqDist(const _ineq_distortions&);
    std::string writeHigherOrderDistortions();

    bool checkCompatibility(uint&, const std::vector<std::string>&);
    void readInequivalentClusters(uint&, const std::vector<std::string>&);
    std::vector<_cluster> readClusters(uint&, const std::vector<std::string>&);
    _linearCombinations readLinearCombinations(uint&, const std::vector<std::string>&);
    void readInequivalentDistortions(uint&, const std::vector<std::string>&);
    _ineq_distortions readIneqDist(uint&, const std::vector<std::string>&);
    void readHigherOrderDistortions(uint&, const std::vector<std::string>&);
  };

// ***************************************************************************

  class AnharmonicIFCs : public xStream {
    // See aflow_aapl_ifcs.cpp for detailed descriptions of the functions
  public:
    AnharmonicIFCs(std::ostream& oss = std::cout) : xStream(oss), clst(oss), directory("./") {}
    AnharmonicIFCs(std::ofstream& mf, std::ostream& oss = std::cout) : xStream(mf, oss), clst(oss), directory("./") {}
    void clear();
    void initialize(const Supercell&, int, const aurostd::xoption&, std::ofstream&, std::ostream& oss = std::cout);
    void initialize(const Supercell&, int, const aurostd::xoption&);

    void setOptions(double, int, double, double, bool);
    [[nodiscard]] const std::string& getDirectory() const;
    void setDirectory(const std::string&);
    [[nodiscard]] int getOrder() const;

    bool runVASPCalculations(_xinput&, _aflags&, _kflags&, _xflags&);
    bool calculateForceConstants();
    [[nodiscard]] const std::vector<std::vector<double>>& getForceConstants() const;
    [[nodiscard]] std::vector<std::vector<int>> getClusters() const;
    void writeIFCsToFile(const std::string&);

  private:
    ClusterSet clst;

    std::vector<_xinput> xInputs;
    bool _useZeroStateForces = false;
    bool initialized = false;
    std::string directory;
    std::vector<std::vector<int>> cart_indices; // A list of all Cartesian indices
    double distortion_magnitude = 0.0; // The magnitude of the distortions in Angstroms
    std::vector<std::vector<double>> force_constants; // Symmetrized IFCs - ME20190520
    int max_iter = 0; // Number of iterations for the sum rules
    double mixing_coefficient = 0.0; // The mixing coefficient for the SCF procedure
    int order = 0; // The order of the IFCs
    double sumrule_threshold = 0.0; // Convergence threshold for the sum rules

    [[nodiscard]] std::string buildRunName(const std::vector<int>&, const std::vector<int>&, int, int) const;
    void applyDistortions(_xinput&, const std::vector<aurostd::xvector<double>>&, const std::vector<int>&, const std::vector<int>&, double = 1.0) const;

    [[nodiscard]] std::vector<std::vector<int>> getCartesianIndices() const;

    std::vector<std::vector<std::vector<aurostd::xvector<double>>>> storeForces(std::vector<_xinput>&);
    std::vector<std::vector<aurostd::xvector<double>>> getForces(int, int&, std::vector<_xinput>&);
    int getTransformedAtom(const std::vector<int>&, const int&);
    void addHigherOrderForces(std::vector<std::vector<std::vector<aurostd::xvector<double>>>>&, int&, std::vector<_xinput>&);
    std::vector<std::vector<double>> calculateUnsymmetrizedIFCs(const std::vector<_ineq_distortions>&, const std::vector<std::vector<std::vector<aurostd::xvector<double>>>>&);
    [[nodiscard]] double finiteDifference(const std::vector<std::vector<aurostd::xvector<double>>>&, int, const std::vector<int>&, const std::vector<int>&) const;

    // Symmetrization Functions
    std::vector<std::vector<double>> symmetrizeIFCs(std::vector<std::vector<double>>);
    typedef std::vector<std::pair<std::vector<int>, std::vector<double>>> tform;
    typedef std::vector<std::vector<std::vector<std::vector<int>>>> v4int;
    void getTensorTransformations(v4int&, std::vector<std::vector<tform>>&);
    std::vector<std::vector<int>> getReducedClusters();
    void applyLinCombs(std::vector<std::vector<double>>&);
    void transformIFCs(const std::vector<std::vector<tform>>&, std::vector<std::vector<double>>&);
    void applyPermutations(std::vector<int>, std::vector<std::vector<double>>&);
    void calcSums(const std::vector<std::vector<int>>&, const std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);
    void correctIFCs(std::vector<std::vector<double>>&, const std::vector<std::vector<double>>&, const std::vector<std::vector<double>>&, const std::vector<std::vector<int>>&, const v4int&);
    std::vector<double> getCorrectionTerms(const int&, const std::vector<std::vector<int>>&, const std::vector<std::vector<double>>&, const std::vector<std::vector<double>>&, const std::vector<std::vector<double>>&);
    [[nodiscard]] uint findReducedCluster(const std::vector<std::vector<int>>&, const std::vector<int>&) const;

      // File I/O
    std::string writeParameters();
    std::string writeIFCs();
    bool checkCompatibility(uint&, const std::vector<std::string>&);
    std::vector<std::vector<double>> readIFCs(uint&, const std::vector<std::string>&);
  };

// ***************************************************************************

  class TCONDCalculator {
    // See aflow_aapl_tcond.cpp for detailed descriptions of the functions
  public:
    TCONDCalculator() = default;
    TCONDCalculator(PhononCalculator&, const aurostd::xoption&);
    void clear();
    void clear(PhononCalculator&);
    void initialize(const aurostd::xoption&);

    double boundary_grain_size = 0.0;
    bool calc_boundary = false;
    bool calc_cumulative = false;
    bool calc_isotope = false;
    bool calc_rta_only = false;
    std::vector<aurostd::xmatrix<aurostd::xcomplex<double>>> eigenvectors; // The eigenvectors at each q-point
    std::vector<std::vector<double>> freq; // The frequencies at each q-point
    std::vector<std::vector<aurostd::xvector<double>>> gvel; // The group velocities
    int nBranches = 0; // The number of branches in the phonon spectrum
    int nIQPs = 0; // The total number of irreducible q-points in the grid
    int nQPs = 0; // The total number of q-points in the grid
    std::vector<std::vector<std::vector<int>>> processes; // The sign, q-point and branch indices of the scattering processes
    std::vector<std::vector<double>> intr_trans_probs; // The intrinsic transition probabilities
    std::vector<std::vector<std::vector<int>>> processes_iso; // The q-point and branch indices of the isotope scattering processes
    std::vector<std::vector<double>> intr_trans_probs_iso; // The intrinsic transition probabilities for isotope processes
    // Scattering rates
    std::vector<std::vector<std::vector<double>>> scattering_rates_total; // total
    std::vector<std::vector<std::vector<double>>> scattering_rates_anharm; // anharmonic scattering
    std::vector<std::vector<double>> scattering_rates_boundary; // boundary scattering
    std::vector<std::vector<double>> scattering_rates_isotope; // isotope scattering

    std::vector<std::vector<std::vector<std::vector<double>>>> phase_space; // Scattering phase space

    std::vector<std::vector<double>> grueneisen_mode; // Mode Grueneisen parameter
    std::vector<double> grueneisen_avg; // Averaged Grueneisen parameter for each temperature

    std::vector<double> temperatures; // The temperatures for the thermal conductivity calculations
    std::vector<aurostd::xmatrix<double>> thermal_conductivity; // The thermal conductivity values

    void calculateThermalConductivity();
    void calculateGrueneisenParameters();
    void writeOutputFiles(const std::string&);

  private:
    PhononCalculator* _pc = nullptr;  // Reference to the phonon calculator
    QMesh* _qm = nullptr;
    bool _pc_set = false;
    bool _initialized = false;

    std::vector<std::vector<double>> calculateModeGrueneisen(const std::vector<std::vector<std::vector<aurostd::xcomplex<double>>>>& phases);
    double calculateAverageGrueneisen(double T, const std::vector<std::vector<double>>&);

    void getWeightsLT(double, const std::vector<double>&, std::vector<double>&);
    void calculateTransitionProbabilities();
    std::vector<std::vector<std::vector<aurostd::xcomplex<double>>>> calculatePhases(bool = false);
    void calculateTransitionProbabilitiesPhonon(int, std::vector<std::vector<std::vector<std::vector<double>>>>&, const std::vector<std::vector<std::vector<aurostd::xcomplex<double>>>>&);
    void calculateTransitionProbabilitiesIsotope(int);
    std::vector<std::vector<double>> calculateTransitionProbabilitiesBoundary();
    void getProcess(const std::vector<int>&, std::vector<int>&, std::vector<int>&, int&);
    aurostd::xmatrix<double> calculateThermalConductivityTensor(double, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);
    std::vector<std::vector<double>> getOccupationNumbers(double);
    std::vector<std::vector<double>> calculateAnharmonicRates(const std::vector<std::vector<double>>&);
    std::vector<std::vector<double>> calculateTotalRates(const std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);
    double getOccupationTerm(const std::vector<std::vector<double>>&, int, const std::vector<int>&, const std::vector<int>&);
    void calcAnharmRates(int, const std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);
    std::vector<std::vector<aurostd::xvector<double>>> getMeanFreeDispRTA(const std::vector<std::vector<double>>&);
    aurostd::xmatrix<double> calcTCOND(double, const std::vector<std::vector<double>>&, const std::vector<std::vector<aurostd::xvector<double>>>&);
    void getMeanFreeDispFull(const std::vector<std::vector<double>>&, const std::vector<std::vector<double>>&, std::vector<std::vector<aurostd::xvector<double>>>&);
    void calculateDelta(int, const std::vector<std::vector<double>>&, const std::vector<std::vector<aurostd::xvector<double>>>&, std::vector<std::vector<aurostd::xvector<double>>>&);
    void correctMFD(const std::vector<std::vector<double>>&, const std::vector<std::vector<aurostd::xvector<double>>>&, std::vector<std::vector<aurostd::xvector<double>>>&);

    void writeTempIndepOutput(const std::string&, std::string, const std::string&, const std::vector<std::vector<double>>&);
    void writeTempDepOutput(const std::string&, std::string, const std::string&, const std::vector<double>&, const std::vector<std::vector<std::vector<double>>>&);
    void writeDataBlock(std::stringstream&, const std::vector<std::vector<double>>&);
    void writePhaseSpace(const std::string&);
    void writeGrueneisen(const std::string&);
    void writeThermalConductivity(const std::string&);
  };

// ***************************************************************************
// END ME: Automatic Anharmonic Phonon Library (AAPL)
// ***************************************************************************

// ***************************************************************************
// BEGIN AS: Quasi-Harmonic Approximation (QHA)
// ***************************************************************************

// AS20200513 BEGIN
#define QHA_ARUN_MODE "QHA" // used in filename
  enum EOSmethod { EOS_MURNAGHAN, EOS_SJ, EOS_BIRCH_MURNAGHAN2, EOS_BIRCH_MURNAGHAN3, EOS_BIRCH_MURNAGHAN4 };
  enum QHAmethod { QHA_CALC, QHA3P_CALC, SCQHA_CALC, QHANP_CALC };
  enum QHAtype { QHA_FD, QHA_EOS, QHA_TE };

  bool QHA_Get_AflowInName(std::string& AflowInName, const std::string& directory_LIB);
  std::string EOSmethod2label(EOSmethod eos_method);// AS20210518
  std::string QHAmethod2label(QHAmethod qha_method);// AS20210518
  bool hasImaginary(const std::string& filename, const std::string& QHA_method);// AS20210813

  /// Calculates QHA-related properties
  class QHA : public xStream {
  public:
    QHA(std::ostream& oss = std::cout) : xStream(oss) {}
    QHA(const xstructure& struc, _xinput& xinput, aurostd::xoption& qhaopts, aurostd::xoption& aplopts, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
    void initialize(const xstructure& struc, _xinput& xinput, aurostd::xoption& qhaopts, aurostd::xoption& aplopts, std::ofstream& FileMESSAGE, std::ostream& oss);
    void run(_xflags& xflags, _aflags& aflags, _kflags& kflags);
    void clear();
    double calcFrequencyFit(double V, aurostd::xvector<double>& xomega);
    double calcGrueneisen(double V, aurostd::xvector<double>& xomega);
    double calcGrueneisen(double V, aurostd::xvector<double>& xomega, double& w);
    double calcGrueneisenFD(const aurostd::xvector<double>& xomega);
    void calcCVandGP(double T, double& CV, double& GP);
    void calcCVandGPfit(double T, double V, double& CV, double& GP);
    double calcGPinfFit(double V);
    double calcVibFreeEnergy(double T, int id);
    double calcFreeEnergyFit(double T, double V, EOSmethod eos_method, QHAmethod method, uint contrib);
    double calcElectronicFreeEnergy(double T, int id);
    double calcChemicalPotential(double T, int Vid);
    double calcIDOS(double e, double T, xEIGENVAL& eig);
    aurostd::xvector<double> calcElectronicFreeEnergySommerfeld(double T);
    aurostd::xvector<double> calcElectronicSpecificHeatSommerfeld(double T);
    aurostd::xvector<double> calcFreeEnergy(double T, QHAmethod qha_method, uint contrib);
    aurostd::xvector<double> calcDOSatEf();
    double calcInternalEnergyFit(double T, double V, EOSmethod method);
    aurostd::xvector<double> fitToEOSmodel(aurostd::xvector<double>& V, aurostd::xvector<double>& E, EOSmethod method);
    aurostd::xvector<double> fitToEOSmodel(aurostd::xvector<double>& E, EOSmethod method);
    double evalEOSmodel(double V, const aurostd::xvector<double>& p, EOSmethod eos_method);
    double calcBulkModulus(double V, const aurostd::xvector<double>& parameters, EOSmethod method);
    double calcBprime(double V, const aurostd::xvector<double>& parameters, EOSmethod method);
    double calcEOS2Pressure(double V, const aurostd::xvector<double>& parameters, EOSmethod method);
    double calcEquilibriumVolume(const aurostd::xvector<double>& parameters, EOSmethod method);
    double calcEntropy(double T, double V, EOSmethod eos_method, QHAmethod method, uint contrib);
    double getEqVolumeT(double T, EOSmethod eos_method, QHAmethod method, uint contrib);
    double calcThermalExpansion(double T, EOSmethod eos_method, QHAmethod method, uint contrib);
    double calcIsochoricSpecificHeat(double T, double V, EOSmethod eos_method, QHAmethod qha_method, uint contrib);
    // QHA3P and SCQHA and QHANP
    double extrapolateFrequency(double V, const aurostd::xvector<double>& xomega, QHAmethod qha_method);
    double extrapolateGrueneisen(double V, const aurostd::xvector<double>& xomega, QHAmethod qha_method);
    // QHA3P
    double calcVibFreeEnergyTaylorExpansion(double T, int Vid, QHAmethod qha_method);
    double calcInternalEnergyTaylorExpansion(double T, double V, QHAmethod qha_method);
    // SCQHA
    double calcVPgamma(double T, double V);
    double calcSCQHAequilibriumVolume(double T, EOSmethod method);
    double calcSCQHAequilibriumVolume(double T, double Vguess, aurostd::xvector<double>& fit_params, EOSmethod method);
    void runSCQHA(EOSmethod method, bool all_iterations_self_consistent = true, const std::string& directory = ".");
    // output
    void writeThermalProperties(EOSmethod eos_method, QHAmethod qha_method, const std::string& directory = ".");
    void writeFVT(const std::string& directory = ".");
    void writeGPpath(double V, const std::string& directory = ".");
    void writeAverageGPfiniteDifferences(const std::string& directory = ".");
    void writeGPmeshFD(const std::string& directory = ".");
    void writeFrequencies(const std::string& directory = ".");
    void writeTphononDispersions(EOSmethod eos_method, QHAmethod qha_method, const std::string& directory = ".");
    void writeQHAresults(const std::string& directory = ".");
      // members
    aurostd::xoption apl_options;
    aurostd::xoption qha_options;
    std::string system_title;
    double EOS_volume_at_equilibrium = 0.0;
    double EOS_energy_at_equilibrium = AUROSTD_NAN;
    double EOS_bulk_modulus_at_equilibrium = 0;
    double EOS_Bprime_at_equilibrium = 0;

  private:
    bool isEOS = false;
    bool isGP_FD = false;
    bool ignore_imaginary = false;
    bool isQHA = false;
    bool isQHA3P = false;
    bool isSCQHA = false;
    bool isQHANP = false;
    bool isInitialized = false;
    bool includeElectronicContribution = false;
    bool doSommerfeldExpansion = false;
    int Ntemperatures = 0;
    int N_GPvolumes = 3;   ///< number of volumes/calculations for finite difference calc
    int N_EOSvolumes = 0;  ///< number of volumes/calculations for EOS calc
    int N_QHANPvolumes = 0; ///< number of volumes/calculations for QHANP calc
    int Nbranches = 0;       ///< number of phonon dispersion branches
    int NatomsOrigCell = 0;  ///< number of atoms in original cell
    int Nelectrons = 0;
    int TaylorExpansionOrder = 0;
    double gp_distortion = 0.0;
    // int NatomsSupercell; ///< number of atoms in supercell
    xstructure origStructure;
    std::vector<double> Temperatures;
    std::vector<int> ph_disp_temperatures;///< temperatures for T-dependent phonon dispersions
    std::vector<double> GPvolumes; ///< a set of volumes for FD Grueneisen calculation
    std::vector<double> EOSvolumes; ///< a set of volumes for EOS calculation
    std::vector<double> QHANPvolumes; ///< a set of volumes for QHANP calculation
    std::vector<double> coefGPVolumes; ///< multiplication coefficient w.r.t initial volume
    std::vector<double> coefEOSVolumes;
    std::vector<double> coefQHANPVolumes;
    aurostd::xvector<double> DOS_Ef;
    // data necessary to calculate thermodynamic properties
    std::vector<double> Efermi_V; ///< Fermi energy vs V
    std::vector<double> E0_V; ///< total energy vs V
    std::vector<xEIGENVAL> static_eigvals;
    std::vector<xIBZKPT> static_ibzkpts;
    std::vector<std::vector<double>> pdos_V; ///< phonon DOS
    std::vector<int> qpWeights;
    std::vector<aurostd::xvector<double>> qPoints;
    // data needed for Grueneisen parameter calculation
    aurostd::xmatrix<double> gp_fit_matrix;
    std::vector<std::vector<std::vector<double>>> omegaV_mesh;
    std::vector<std::vector<std::vector<double>>> omegaV_mesh_EOS;
    std::vector<std::vector<std::vector<double>>> omegaV_mesh_QHANP;
    std::vector<xEIGENVAL> gp_ph_dispersions;
    std::vector<xEIGENVAL> eos_ph_dispersions;
    std::vector<ThermalPropertiesCalculator> eos_vib_thermal_properties;
      //
    std::vector<std::string> subdirectories_apl_eos;
    std::vector<std::string> subdirectories_apl_gp;
    std::vector<std::string> subdirectories_apl_qhanp;
    std::vector<std::string> subdirectories_static;
    std::vector<std::string> arun_runnames_apl_eos;
    std::vector<std::string> arun_runnames_apl_gp;
    std::vector<std::string> arun_runnames_apl_qhanp;
    std::vector<std::string> arun_runnames_static;
    _xinput xinput;
    std::string currentDirectory;
      // methods
      // related to static DFT calculations
    void createSubdirectoriesStaticRun(const _xflags& xflags, const _aflags& aflags, const _kflags& kflags, const std::vector<std::vector<bool>>& list);
    int checkStaticCalculations(std::vector<std::vector<bool>>& file_is_present);
    void printMissingStaticFiles(const std::vector<std::vector<bool>>& list, const std::vector<std::string>& subdirectories);
    bool readStaticCalculationsData();
      // related to APL calculations
    void createSubdirectoriesAPLRun(const _xflags& xflags, const _aflags& aflags, const _kflags& kflags, const std::vector<std::vector<bool>>& list, QHAtype qhatype);
    int checkAPLCalculations(std::vector<std::vector<bool>>& file_is_present, QHAtype qhatype);
    void printMissingAPLFiles(const std::vector<std::vector<bool>>& list, QHAtype qhatype);
    bool readAPLCalculationData(const std::vector<std::string>& subdirectories, _kflags& kflags, QHAtype type);
    bool runAPL(_xflags& xflags, _aflags& aflags, _kflags& kflags, QHAtype qhatype);
  };
// AS20200513 END

// ***************************************************************************
// END AS: Quasi-Harmonic Approximation (QHA)
// ***************************************************************************
} // namespace apl

#endif  // _AFLOW_APL_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
