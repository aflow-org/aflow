// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

#ifndef _AFLOW_APENNSY_H_
#define _AFLOW_APENNSY_H_
// ---------------------------------------------------------------------------
#include "aflow.h"
#include "aflowlib.h"

// DEFINITIONS !!!

#define XGNDSTATE_HOLES 1    // 1 alloys matrix with holes
#define MAX_NUM_SPECIES 4

#define AFLOWLIB_PROJECT_GNDSTATE        string("materials.duke.edu/AFLOWDATA/LIB2_RAW")

#define APENNSY_INF 1e9
#define APENNSY_LATEX_DOC_TYPE        "[aps,amssymb]{article}" // "revtex4-1"
#define APENNSY_LATEX_DOC_TYPE_TWOCOLUMS        "[aps,amssymb]{article}" // "revtex4-1"   // two columns
#define APENNSY_LATEX_FONT string("\\footnotesize")
#define APENNSY_LATEX_FONT_SMALL string("\\small")
#define APENNSY_LATEX_FONT_FOOTNOTESIZE string("\\footnotesize")
#define APENNSY_LATEX_FONT_TINY string("\\tiny")

#define C0000000 (0.0)             // 0
#define C1000000 (1.0)             // 1

#define CEPSILON 0.0001

#ifdef ZUNGER_PROTO
#define _STR5_BETA1   string("\\beta1")
#define _STR6_BETA2   string("\\beta2")
#define _STR9_ALPHA1  string("\\alpha1")
#define _STR10_ALPHA2 string("\\alpha2")
#define _STR13_Z1     string("Z1")
#define _STR14_Z2     string("Z2")
#define _STR15_Z3     string("Z3")
#define _STR16_W3     string("W3")
#define _STR18_W3     string("W3")
#define _STR17_W2     string("W2")
#define _STR19_Y1     string("Y1")
#define _STR20_Y2     string("Y2")
#define _STR21_Y3     string("Y3")
#define _STR23_CH     string("CH")
#define _STR27_V1     string("V1")
#define _STR28_V2     string("V2")
#define _STR29_V3     string("V3")
#else
#define _STR5_BETA1   string("fcc_{AB2}^{[001]}")
#define _STR6_BETA2   string("fcc_{AB2}^{[001]}")
#define _STR9_ALPHA1  string("fcc_{AB2}^{[111]}")
#define _STR10_ALPH2  string("fcc_{AB2}^{[111]}")
#define _STR13_Z1     string("fcc_{AB3}^{[001]}")
#define _STR14_Z2     string("fcc_{A2B2}^{[001]}")
#define _STR15_Z3     string("fcc_{AB3}^{[001]}")
#define _STR16_W3     string("fcc_{AB3}^{[113]}")
#define _STR18_W3     string("fcc_{AB3}^{[113]}")
#define _STR17_W2     string("fcc_{A2B2}^{[113]}")
#define _STR19_Y1     string("fcc_{AB3}^{[011]}")
#define _STR20_Y2     string("fcc_{A2B2}^{[011]}")
#define _STR21_Y3     string("fcc_{AB3}^{[011]}")
//#define _STR23_CH   string("NbAs")
#define _STR23_CH     string("fcc_{A2B2}^{[201]}")
#define _STR27_V1     string("fcc_{AB3}^{[111]}")
#define _STR28_V2     string("fcc_{A2B2}^{[111]}")
#define _STR29_V3     string("fcc_{AB3}^{[111]}")

#define _STR60     string("bcc_{AB}^{[101]}")
#define _STR61     string("bcc_{AB}^{[001]}")
#define _STR62     string("bcc_{AB2}^{[211]}")
#define _STR63     string("bcc_{AB2}^{[211]}")
#define _STR64     string("bcc_{AB2}^{[011]}")
#define _STR65     string("bcc_{AB2}^{[011]}") 
#define _STR66     string("bcc_{AB2}^{[001]}") // C11_b
#define _STR67     string("bcc_{AB2}^{[001]}") // C11_b
#define _STR76     string("bcc_{AB3}^{[001]}") 
#define _STR78     string("bcc_{AB3}^{[001]}")  
#define _STR77     string("bcc_{A2B2}^{[001]}") // gammaCuTi

#endif

#define dENERGY_CUTOFF 0.015
#define dVOLUME_CUTOFF 0.010 // 0.070
#define dBONDSD_CUTOFF 0.010 //0.100

//#define LIB_MODE1 1
#define LIB_MODE2 2
#define LIB_MODEN 4
#define LIB_MODEX 5

#define APENNSY_INT_DEPTH 3
#define APENNSY_STR_DEPTH 14
#define APENNSY_DBL_DEPTH 8
#define APENNSY_LSTR_DEPTH 12

#define apennsy_CUT_OFF 0.00001
#define apennsy_epsilon 0.00001
#define apennsy_TSTART_IN_DEFAULT 20.0
#define _MATLAB_OUTPUT_MODE_      1
#define _GNUPLOT_OUTPUT_MODE_     2
#define _MATPLOTLIB_OUTPUT_MODE_  3

class APENNSY_Parameters;
bool strsub(string strzero, string subA);
bool strsub(string strzero, string subA, string subB);
bool strsub(string strzero, string subA, string subB, string subC);
string ElementSymbolToString(string stringin);
int Apennsymain(vector<string> &argv,vector<string> &cmds);
string README_AFLOW_APENNSY_TXT(void);
string CleanNameStructure(string stringin);
bool _is_alphabetic(string alloy);
uint SGtoNSG(string sgroup);
bool xphaseseparation(const xvector<int>& converhull, const int& nalloy,const APENNSY_Parameters &params);
xvector<int> xconvexhull(int nalloy,const APENNSY_Parameters &params);
xvector<int> xminenergy(int nalloy,const APENNSY_Parameters &params);
double GNDdistance(const int& i, const int& j,const APENNSY_Parameters &params);
xvector<double> GNDdistance(const int& j,const APENNSY_Parameters &params);
double TIELINEdistance(const int& i, const int& j,const APENNSY_Parameters &params);
xvector<double> TIELINEdistance(const int& j,const APENNSY_Parameters &params);
double MINEdistance(const int& i, const int& j,const APENNSY_Parameters &params);
xvector<double> MINEdistance(const int& j,const APENNSY_Parameters &params);
string APENNSY_MiscibilityString(int flag);

// ----------------------------------
// CLASSES

class APENNSY_Parameters {
 public:
  friend class aflowlib::_aflowlib_entry;
  // trivial constructurs/destuctors/operators
  APENNSY_Parameters();                              // default, just allocate
  ~APENNSY_Parameters();                             // kill everything
  APENNSY_Parameters(const APENNSY_Parameters& b);    // constructor copy
  const APENNSY_Parameters& operator=(const APENNSY_Parameters &b); // copy
  void clear(void);                                 // clear
  // FUNCTIONS
  bool LoadDynamicMemory(bool _verbose);             // Load dynamic memory
  // bool LoadStructuresHQTC(_aflags& aflags);    // HTQC structures
  // bool LoadStructuresGUS(_aflags& aflags);     // GUS structures
  bool LoadLibrary(_aflags &aflags);           // Library structures
  bool LibLoadAlloysLIB2(_aflags &aflags);     // LOAD LibraryX
  bool LibLoadAlloysLIB2U(_aflags &aflags);     // LOAD LibraryU
  bool LibLoadAlloysLIB2PGM(_aflags &aflags);     // LOAD LibraryPGM
  bool LibLoadAlloysALLOY(string,_aflags &aflags); // LOAD one system
  bool SplitAlloySpecies(void);                     // Fix Species names
  bool SplitAlloyPseudoPotentials(void);            // Fix PseudoPotentials Names
  bool FixGndStatesNamesConcentrations(bool verb);  // xphases.h
  void CheckAlloyRelaxationStructures(void);        // check identical structures
  string MatlabGndStatesNamesConcentrations(bool);  // xphases.h
  string APENNSY_Data(bool _verbose,_aflags &aflags);                 // apennsy.cpp
  string APENNSY_UNCLE(bool _verbose,_aflags &aflags);                // apennsy.cpp
  string APENNSY_Web(bool _verbose,_aflags &aflags);                  // apennsy.cpp
  string APENNSY_Reference(bool _verbose,_aflags &aflags);            // apennsy.cpp
  string APENNSY_EnergyList(bool _verbose,_aflags &aflags); // apennsy.cpp with optional label
  string APENNSY_ConvexHull(bool _verbose,_aflags &aflags,uint mode); // apennsy.cpp
  string APENNSY_PS_EnergyList(_aflags &aflags);                      // apennsy.cpp
  string APENNSY_HistogramList(_aflags &aflags);                      // apennsy.cpp
  string APENNSY_Rules(_aflags &aflags);                              // apennsy.cpp
  string APENNSY_Structures(_aflags &aflags);                         // apennsy.cpp
  string APENNSY_Order(_aflags &aflags);                              // apennsy.cpp
  string APENNSY_SimulationInformation(_aflags &aflags,string mode);  // apennsy.cpp
  string APENNSY_Miscibility(_aflags &aflags);                        // apennsy.cpp
  string APENNSY_Miscibility_Experiments(_aflags &aflags);            // apennsy.cpp
  string APENNSY_Miscibility_Miedema(_aflags &aflags);                // apennsy.cpp
  string APENNSY_Miscibility_HumeRothery(_aflags &aflags);            // apennsy.cpp
  string APENNSY_Miscibility_Table(_aflags &aflags);                  // apennsy.cpp
  string APENNSY_Miscibility_Statistics(_aflags &aflags);             // apennsy.cpp
  string APENNSY_StructureVolumes(bool _verbose,_aflags &aflags);     // apennsy.cpp
// [OBSOLETE]  string APENNSY_ProtoCheck(bool _verbose,_aflags &aflags);           // apennsy.cpp
  void MakeRankLib(void);
 // CONTENT
  bool PseudopotentialNoclean;                       // to overcome some issues.
  //  uint Nconcentrations[MAX_NUM_SPECIES+1];       // Nconcentrations[number of species] became ZConcentrations.at(k).size()-1
  uint LIBstructures_calcs;                          // LIBstructures_calcs
  uint LIBstructures_holes;                          // LIBstructures_holes
  int LIB_MODE;                                      // LIB_MODE
  bool _flag_ce_cswap_;                              // ce_cswap for morons
  bool _flag_load_xstructures;                       // will load structures
  bool Load_ALL,Load_FCC,Load_BCC,Load_HCP;          // subset of load
  stringstream alloysmesg;                           // Messages
  vector<string> alloys,alloysRAW;                   // (alloy.size())
  vector<string> speciesA,speciesB,speciesAB;        // (alloy.size())
  vector<string> pseudosA,pseudosB;                  // (alloy.size())
  vector<int> Alloy2MiscibilityHT;                   // (alloy.size())
  vector<int> Alloy2MiscibilityExperiments;          // (alloy.size())
  vector<int> Alloy2MiscibilityMiedema;              // (alloy.size())
  vector<int> Alloy2MiscibilityHumeRothery;          // (alloy.size())
  vector<string> PhaseDiagramNameLeft;               // (alloy.size())
  vector<string> PhaseDiagramNameRight;              // (alloy.size())
  vector<double> AlloyOrder;                         // (alloy.size()) min of order_per_atom for each alloy
  //  vector<double> C[MAX_NUM_SPECIES+1];               // (Nconcentrations)
  vector<vector< int> > ZGNDlibrary;    // (alloy.size(),ZConcentrations.size()), (alloy.size(),*)  // contains the index of the GNDs to be used only if >=0 (-1, no gndstate).
  //  vector<vector<uint> > DGNDlibrary;     // (alloy.size(),ZConcentrations.size()), (alloy.size(),*)
  vector<vector<uint> > ZMINElibrary;                // (alloy.size(),ZConcentrations.size())
  vector<vector<double> > ZConcentrations;           // (alloy.size(),dynamic) replace
  vector<vector<vector<bool> > > AlloyStructureIdentity;  // (alloy.size(),NstructuresHTQC,NstructuresHTQC)
  vector<vector<vector<int> > > RankLib;                  // (alloy.size(),ZConcentrations.size(),*)
  vector<vector<vector<int> > > RankLibInt;               // (alloy.size(),ZConcentrations.size(),*)
  vector<vector<aflowlib::_aflowlib_entry> > ZLibrary;  // (alloy.size(),Nstructures_dynamic)
  // ALibrary stuff (properties for each Nalloys)    // (alloy.size())
  vector<double> enthalpy_formation_atom_max;        // (alloy.size()) max of enthalpy_formation_atom (negative)
  vector<double> entropic_temperature_max;           // (alloy.size()) max of entropic_temperature (positive)
  vector<double> alloy_order;                        // (alloy.size()) min of order_per_atom for each alloy
  vector<uint> alloy_holes;                          // (alloy.size())
  vector<double> alloy_deltamuBArich;                // (alloy.size()) delta_muB^A+
  vector<double> alloy_deltamuABrich;                // (alloy.size()) delta_muA^B+
  // REFERENCES
  // vector<uint> A_pureA;
  // vector<uint> A_pureB;
 private:                                           //
  void free();                                      // free space
  void copy(const APENNSY_Parameters& b);            // copy
};

// ----------------------------------------------------------------------------
// aflow_apennsy_vaspin.cpp

namespace apennsy_std {
  void LatexSideMarginThesisON(void);
  void LatexSideMarginThesisOFF(void);
  void ZungerConversionChart(void);
  void read_VASPIN(ostringstream* l,int number,
		   const char* _labels,const char* _strukturbericht,
		   const char* _superstr,const char* _superlatt,
		   const char* _proto,bool VERB);
  void MakeVaspIn(bool VERB);
  void ZungerConversionChart(void);
}

// ----------------------------------------------------------------------------

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
