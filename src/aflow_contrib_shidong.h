// ***************************************************************************
// *                                                                         *
// *               AFlow SHIDONG WANG - Duke University 2010-2011            *
// *                                                                         *
// ***************************************************************************
// aflow_contrib_shidong.cpp
// functions written by
// 2010-2011: shidong.wang@duke.edu

#ifndef _AFLOW_CONTRIB_SHIDONG_H_
#define _AFLOW_CONTRIB_SHIDONG_H_

// put your prototypes here

#include "aflow.h"
#include "aflow_contrib_shidong_cluster_expansion.h"
#include "aflow_contrib_shidong_auxiliary.h"

#define _COMMENT_CHAR '#'

const string _FILENAME="clusters.dat";
const string _CORFILENAME="cor.dat";
const string _ECIFILENAME="eci.dat";

const int _SLname_length = 8;
const int _SLname_pos = 5;

// cross-validation
//const string _LOCVFILE="LOCV.dat";
//const string _SCOREFILE="ga_score.dat";
//const int _TRIAL_NUM = 50; // number of configuration used in cross-validation function
const int _TRIAL_NUM = 150;
const double _TRIAL_RATIO = 0.3;

const int _LEFTOUT_NUM = 5;
//const int _LEFTOUT_NUM = 10;

// stepwise algorithm
//const int _INITIAL_NUM = 4; // number of different trials
//const int _INITIAL_NUM = 6; // number of different trials
//const int _INITIAL_NUM = 20; // number of different trials
const int _INITIAL_NUM = 10; // number of different trials
const double _SCORE_ZERO = 1.0e-8; // zero score, exactly fit

// genetic algorithm
const int _POPULATION = 30;
const double _MUTATION_PROPABILITY = 0.2;
const int _GENE_NUM = 20;
const double _REPRODUCTION_RATE = 0.6;
const double _SP = 2.0;
const int _EVOLUTION_NUM = 100;
const int _NO_REPLACEMENT_NUM = 2;
const int _SAME_BEST_CHROMOSOME_NUM = 5;

bool ACEGetStructureType(const string& structure_type_in);
vector<cestructure>  ReadInFitStructure(istream & os, string & structure_type);
void WriteOutFitStructure(ostream & os, string & alloy_name,
        vector<cestructure> & str_list);

//*************************************
// Convex Hull functions
//*************************************
vector<int> CEConvexHull(vector< vector<double> > & points,
        const string & option="min");
bool comparison_points(const vector<double>  & point1,
        const vector<double> & point2);
vector<int> GroundStateCandidate(vector< vector<double> > & points);

void ACEGroundStateConvexHull(ifstream & os, vector<cestructure> str_list);
struct cepair
{
    string name;
    double num;
};
vector<string> GetNewHullState(vector<cestructure> & str_list);
vector<string> GetNewGroundState(vector<cestructure> & str_list);
bool comparison_pair(cepair pari1, cepair pair2);

//*************************************
// Optimize ECI
//*************************************
struct fittness
{
    double score;
    double rank;
    int index;
    int size;
};
double ACECrossValidation(vector<cestructure> & str_list, ceECIcluster & ECIcluster,
        double energy_A, double energy_B);

double ACECrossValidation(int left_out_num,
        vector<cestructure> & str_list, ceECIcluster & ECIcluster);

double ACECrossValidation( vector< vector<int> > leftout_list_set,
        vector<cestructure> & str_list, ceECIcluster & ECIcluster);

double ACECrossValidationRandom(int left_out_num,
        vector<cestructure> & str_list, ceECIcluster & ECIcluster);

vector< vector<int> >  ACESetUpCrossValidationSetsRandom(int left_out_num,
        int total_str);

void ACEOptimalECIClusters(vector<cestructure> & str_list,
        ceECIcluster & ECIcluster);
vector<int> StepWise( vector<cestructure> & str_list,
        ceECIcluster & ECIcluster);
vector<int> GeneticAlgorithm( vector<cestructure> & str_list,
        ceECIcluster & ECIcluster);
void PrintChromosome(ostream & os, bool *chromosome[], int col, int row);
vector<int> GetTrialSet(bool *chromosome, int col);
bool comparison_fitscore(const fittness & score1, const fittness & score2);

vector<int> LeftOutSet(int total_size, vector<int> kept_list);

void PrintECI(ostream & os, ceECIcluster & ECIcluster);

//*************************************
// Get ECI
//*************************************
void ACEGetECI(vector<cestructure> & str_list, ceECIcluster & cluster1);
void ACEGetECI(vector<cestructure> & str_list, vector<int> fit_list, ceECIcluster & cluster1);
vector<double> ACEGetECIAverage(vector< vector<int> > str_leftout_set_list,
        vector<cestructure> & str_list, ceECIcluster & ECIcluster);

void ACEGenerateSL(xmatrix<int> N_mat, string structure_type, ostream & os, ostream & os1);

void ACESLProperties(istream & os, string & structure_type,
        const ceallclusters & allcluster1,
        ceECIcluster & ecicluster1);
void ACESLProperties(ostream & os, string & structure_type,
        int cell_nr_min, int cell_nr_max, const ceallclusters & allcluster1,
        ceECIcluster & ecicluster1);

void ACESLProperties_Readin_Corfile(istream & os, string & structure_type,
        const ceallclusters & allcluster1,
        ceECIcluster & ecicluster1);

void ACESLCorrelations(istream & os, string & structure_type,
        const ceallclusters & allcluster1,
        ceECIcluster & ecicluster1);

double ACECrossValidation( vector< vector<int> > leftout_list_set,
        vector<cestructure> & str_list, ceECIcluster & ECIcluster);

vector< vector<int> >  ACESetUpCrossValidationSetsRandom(int left_out_num,
        int total_str);

//*************************************
// Generate Superlattices
//*************************************
void ACEGenerateSL(xmatrix<int> N_mat, string structure_type, ostream & os, ostream & os1);

void ACESLProperties(istream & os, string & structure_type,
        const ceallclusters & allcluster1,
        ceECIcluster & ecicluster1, double temperature,
        double energy_A, double energy_B);
void ACESLProperties(ostream & os, string & structure_type,
        int cell_nr_min, int cell_nr_max, const ceallclusters & allcluster1,
        ceECIcluster & ecicluster1, double temperature,
        double energy_A, double energy_B);

void ACESLProperties_Readin_Corfile(istream & os, string & structure_type,
        const ceallclusters & allcluster1,
        ceECIcluster & ecicluster1, double temperature,
        double energy_A, double energy_B);

// [OBSOLETE] const string _AFLOW_INPUT_TEMPLATE = "aflow.in_template";

void GenerateAflowInputFile(string structure_type,
        string AlloyName, string SL_name,
        vector<_ceatom> atom_species, bool mpi_flag);

//*************************************
// Get Special Quasirandom Structure (SQS)
//*************************************
struct SLSQS
{
    double weight;
    ceSL SQS;
    double stoich_b;
    int cell_nr;
};

void ACESLPropertiesSQS_Readin_Corfile(istream & os, string & structure_type,
        const ceallclusters & allcluster1,
        ceECIcluster & ecicluster1,
        int site_num, int NNNum,
        int min_SLcell_nr, int max_SLcell_nr,
        vector<_ceatom> atom_species);
bool comparison_SLSQS(const vector<SLSQS> sqs1_list, const vector<SLSQS> sqs2_list);
bool is_equal_structure(ceSL sl1, ceSL sl2);

//*************************************
// Get Carrier Effective Mass
//*************************************

// [OBSOLETE] struct kEn_st{
// [OBSOLETE]  double kx,ky,kz;
// [OBSOLETE]  double energy[2];
// [OBSOLETE]  int band_index;
// [OBSOLETE]  int band_type; // 0 -- valence band; 1 -- conduction band
// [OBSOLETE]};

// [OBSOLETE]struct band_st
// [OBSOLETE]{
// [OBSOLETE]    double Efermi;
// [OBSOLETE]    double vbt[2];
// [OBSOLETE]    double cbb[2];
// [OBSOLETE]    double band_gap;
// [OBSOLETE]};

// range of energy point to fit the ellipse curve
// [OBSOLETE] const double _FIT_ENERGY_RANGE = 0.026; // eV range of band
// [OBSOLETE] const int _FIT_POINTS_NUMBER = 8; // minimum fit points in Irreducible BZ

//range of band extremes to determine the number of bands for effective mass calculations
// [OBSOLETE] const double _BANDS_ENERGY_RANGE = 0.026; // eV

// used to determine cluster of points
// can be changed to other values
// [OBSOLETE]const double _MIN_RATIO = 0.2;

// factor unit
// mass is in unit of electron mass
// [OBSOLETE] const double _MASS_FACTOR = 3.80998; // hbar^2*10^{20}/(2.0*me*eV)

// [OBSOLETE]bool getEfermiBandGap(string & file_dos, int ispin, string & compound_name, band_st &band_information);

// [OBSOLETE]bool getFitDataForMass(const string & file_kE, 
// [OBSOLETE]		       const string & file_lattice,
// [OBSOLETE]		       const band_st & band_information,
// [OBSOLETE]		       vector<vector<int> > & number_of_valley_list, 
// [OBSOLETE]		       double energy_range, 
// [OBSOLETE]		       vector<vector<vector<kEn_st> > > &fit_data_all
// [OBSOLETE]		       );
// [OBSOLETE]
// [OBSOLETE]bool fitToEllipsisEquation(vector<vector<kEn_st> > & fit_data, int spin_idx, vector<vector<double> > &mass_eff);

// [OBSOLETE]bool comparison_kEn_str_up(const kEn_st & k1, const kEn_st & k2); // up_spin
// [OBSOLETE]bool comparison_kEn_str_down(const kEn_st & k1, const kEn_st & k2); // down_spin
// [OBSOLETE]bool comparison_kEn_str_position(const kEn_st & k1, const kEn_st & k2);
// [OBSOLETE]bool comparison_kEn_str_band_type_up(const kEn_st & k1, const kEn_st & k2);
// [OBSOLETE]bool comparison_kEn_str_band_type_down(const kEn_st & k1, const kEn_st & k2);
// [OBSOLETE]bool near_to(const xvector<double> & k1, const xvector<double> & k2, const vector<double> & max_distance);
// [OBSOLETE]bool is_equal_position_kEn_str(const kEn_st & k1, const kEn_st & k2);

#endif

// ***************************************************************************
// *                                                                         *
// *               AFlow SHIDONG WANG - Duke University 2010-2011            *
// *                                                                         *
// ***************************************************************************
