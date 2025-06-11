// ***************************************************************************
// *                                                                         *
// *               AFlow SHIDONG WANG - Duke University 2010-2011            *
// *                                                                         *
// ***************************************************************************
// aflow_contrib_shidong.cpp
// functions written by
// 2010-2011: shidong.wang@duke.edu
// Some constants

#ifndef _AFLOW_CONTRIB_SHIDONG_CLUSTER_EXPANSION_H_
#define _AFLOW_CONTRIB_SHIDONG_CLUSTER_EXPANSION_H_

#include "aflow.h"
//#include "aflow_pflow.h"
#include <cstdlib>
#include <ctime>
#include <algorithm>

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// constant.h

#define _SKIPPED_LINES 2 // skipped lines in input files
#define _EQUAL_DOUBLE 1.0e-9 // two doulbes are equal if difference is smaller than it
#define _EQUAL_ATOM 1.0e-4 // two doulbes are equal if difference is smaller than it
#define _DIM 3 // dimensionality of the system
#define _RADIUS 2 //range to get atoms related by translational symmetry operator
#define _SIGMA 1.0 // default standard deviation of input data
#define _SITE_NUM_1 3
#define _SITE_NUM_2 8
#define _NNNUM_1 4
#define _NNNUM_2 2
#define _DIST_NUM_1 4
#define _DIST_NUM_2 2
#define _NSTEP 50 // step in calculating alloy properties

// this depends on the definition of structure in the database
#define _FCC_BEGIN 1 // range of FCC structures label with number
#define _FCC_END 29 // range of FCC structures label with number
#define _BCC_BEGIN 58 // range of FCC structures label with number
#define _BCC_END 86 // range of FCC structures label with number
#define _HCP_BEGIN 115 // range of FCC structures label with number
#define _HCP_END 596 // range of FCC structures label with number

// Exit status
#define _EXIT_WRONGTYPE 50
#define _EXIT_NOSTRUCTURE 51
#define _EXIT_RANGE_ERROR 52
#define _EXIT_NO_ATOM_IN_CLUSTER 53
#define _EXIT_NO_SYMMETRIC_OPERATOR 54
#define _EXIT_NO_CONVERGENCE 55
#define _EXIT_RANK_NOT_MATCH 56
#define _EXIT_FAIL 57
#define _EXIT_CORRELATION_NOT_FOUND 58
#define _EXIT_CLUSTER_NOT_FOUND 58
#define _EXIT_NO_INPUTFILE 59
#define _EXIT_SVD_FAIL 101
#define _EXIT_FIT_NUM_TOO_LARGE 102

enum {_A_ATOM = -1, _B_ATOM = 1};

const int _DEFAULT_PRECISION=6; // default precison of iostream object

struct _ceatom {
  string name;
  double volume;
};

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// cexstructure.h

// xstructre is too big and only a small fraction information is needed
class cecstructure {
  private:
    ;
  public:
    //std::deque<int> num_each_type;
    //std::deque<_atom> atoms;    

    deque<_atom> atoms;    
    vector<int> indices; // index of atom in atom_list in ceallcluster class


    cecstructure();
    cecstructure(const cecstructure & structure_in);
    ~cecstructure();

    cecstructure & operator=(const cecstructure & structure_in);

};

//////////////////////////////////////////////////////////////
// cexstructure class -- crystal structure
//////////////////////////////////////////////////////////////
// store the informations in xstructre relevent to CE
// an object of cexstructure class should be always initialized
// by an object of xstructure or similiar classes

class cexstructure : public cecstructure {
  private:
    ;

  public:
    xmatrix<double> lattice;
    xmatrix<double> klattice;
    xmatrix<double> c2f;
    xmatrix<double> f2c;
    string prototype;
    std::deque<int> num_each_type;

    cexstructure();
    cexstructure(xstructure & xstr);
    ~cexstructure();

    // For dimensions other than 3D, change xstructure
    // to the correspoding one
    void SetStructure(xstructure & xstr);
    xstructure SetXstructure();

};

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// cecluster.h 

class cecluster {
  // cluster data type

  private:
    ;
  public:
    string pair_name;
    int equivalent_num;
    int site_num;
    int NNNum;
    int dist;
    int index;
    cecstructure structure;

    cecluster();
    cecluster(const cecluster & cecluster_in);
    ~cecluster();

    // reload operators
    cecluster & operator=(const cecluster & cecluster_in);

    // get the largest neareast neighbor shell index in a cluster
    int GetNNNum(vector<double> NN_distance);
    double ClusterDistance(); // get the sum of distance of any two atoms in a cluster
};



cecluster Tetrahedron();
cecluster Octahedron();
cecluster DoubleTetrahedron();
cecluster Triplet();
cecluster Pair();


vector<int> AllCombination41(int num, int total_num, int index);
vector<int> AllCombination42(int num, int total_num, vector<int> & str);
unsigned long int CombinationNr(int num, int total_num);

bool is_equal(const double f1, const double f2);
bool is_equal(const _atom& atom1, const _atom& atom2);
bool is_equal_fpos(const _atom& atom1, const _atom& atom2);
bool is_equal(const cecluster & cluster1, const cecluster & cluster2);
bool is_equal_crystal(const cecluster & cluster1, const cecluster & cluster2);


// comparison of two clusters
bool comparison_cluster(const cecluster& cluster1, const cecluster& cluster2);


// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ceallcluster.h

//const int _NNNUM_MAX = 9;
//const int _NNUP_LIST[10] = {0, _NNNUM_MAX, 5, 4, 4, 3, 2, 2, 2, 2};
// point pair triple ...
//const int _NNNUM_MAX = 2;
//const int _NNUP_LIST[6] = {0, 2, 2, 2, 2, 2};

// parameters for reading in clusters
// only read in those cluster with dist smaller than values given below
const int _DIST_LIST[10] = {4, 4, 4, 4, 3, 3, 3, 3, 3, 3};

bool comparison_cluster_list(const vector<cecluster> & cluster_list_1,
    const vector<cecluster> & cluster_list_2);

// classes

class ceallclusters {
  // clusters in a structure

  private:
    string name; // base structure name fcc/bcc/hcp
    xstructure structure; // base structure
    int num_rep_cluster;
    int num_total_cluster;

    vector<double> NN_distance;
    vector<int> NN_shell_num;
    deque<_atom> atom_list; // all atoms upto NNNum NN shell

    // return the number of all cluster
    int GenerateAllCluster(int SiteNum_low, int SiteNum_up,
        int NNNum_low, int NNNum_up);
    // return the number of rep_cluster
    int GetRepresentCluster();

    int GetNearestNeighbor(int NNNum); // return the number of atoms generated
    int GetNearestNeighbor();

    void GetClusterDistance(); // add the third index in pair_name
    int EquivalentCluster(const vector<_sym_op> pgroup);

  public:
    vector<cecluster> rep_cluster; // representive clusters
    vector< vector<cecluster> > all_cluster; // all clusters grouped by symmetry

    ceallclusters();
    ceallclusters(string & str_name);
    ceallclusters(const ceallclusters & cluster1);
    ~ceallclusters();


    // overloaded operator methods
    friend ostream & operator<<(ostream & os, const ceallclusters & cestr);
    friend ostream & operator>>(ofstream & os, const ceallclusters & cestr);
    ceallclusters & operator=(const ceallclusters & cestr);

    // functions

    bool ReadIn(string & filename); // read data from a file
    bool WriteFile(string & filename, string & stat); // write data to a file


    string Name() const {return name;};
    xstructure  Structure() const {return structure;};
    vector<double>  NNDistance() {return NN_distance;};
    vector<int>  NNShellNum() {return NN_shell_num;};
    deque<_atom> AtomList() const {return atom_list;};

    void PrintRepCluster() const;
    void PrintAllCluster() const;

    // generate all clusters and get the representive ones
    // return the number of the representive clusters
    int SetCluster( const int& SiteNum_low_in, const int& SiteNum_up_in,
        const int& NNNum_low_in, const int& NNNum_up_in);
    int SetCluster( string & filename );

    string GetNamebyCluster(cecluster & cluster_in) const;
    int GetIndexbyCluster(cecluster & cluster_in) const;

};

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ceECIcluster.h

const double _ECI_ZERO = 1.0e-8; // zero value ECI

class ceECIcluster
{
  private:
    string name; // structure name;

    vector<int> allECI_cluster; // clusters used in ECI given by index in rep_cluster list
    vector<int> ECI_cluster; // clusters used in ECI given by index in rep_cluster list
    vector<double> ECI; // Effecitve Cluster Interactions
    double chisq;

    const ceallclusters * str_cluster; // all clusters of a structure
    bool str_cluster_calculated;


  public:
    ceECIcluster();
    ~ceECIcluster();

    // copy constructor
    ceECIcluster(ceECIcluster & ECIcluster);

    // overloaded operator methods
    ceECIcluster & operator=(ceECIcluster & ECIcluster);

    // functions
    void SetUp(ceallclusters & cluster1);

    void SetStrCluster(const ceallclusters & ceallcluster_in);

    ceallclusters GetStrCluster() {cerr << "calculated? " << str_cluster_calculated << endl; return *str_cluster;};

    void GetECICluster();
    void PrintOutECICluster();
    void PrintECI(ostream & os);

    void GetECICluster(vector<int> & ECI_cluster_in);

    vector<int> AllECICluster(){return allECI_cluster;};
    vector<int> ECICluster(){return ECI_cluster;};
    vector<double> ECIValue(){return ECI;};
    double ChiSQ(){return chisq;};

    void SetECI(vector<double> & ECI_in){ECI = ECI_in;};
    void SetChiSQ(double chisq_in){chisq=chisq_in;};
    void DeleteZeroECI();

    string GetClusterNameByIndex(int index);
    void PrintOutProbability();

    void WriteFile(ostream & os);
    void ReadIn(istream & os);
};

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// cestructure.h

const double _pi = 3.1415926536;
const int _NOCORRELATON_FOUND = -1;

class cestructure
{
  // structure used in cluster expansion calculation
  protected:
    string name;
    //xstructure structure;
    cexstructure structure;
    vector<double> ECI_correlation; // correlations of ECI clusters
    vector<int> ECI_equivalent_num;
    vector<int> ECI_cluster; // clusters used in ECI given by index in rep_cluster list
    vector<double> allECI_correlation;
    vector<int> allECI_equivalent_num;
    vector<int> allECI_cluster; // clusters used in ECI given by index in rep_cluster list

    vector<double> ECI; // Effecitve Cluster Interactions
    double chisq; // error in least squre fit

    const ceallclusters * str_cluster; // all clusters of a structure
    bool str_cluster_calculated;

    // energy is the value to be fit by ECI
    // it can be either internal energy, formation energy or any other
    // quantities which can be fit by ECI, e.g., thermal conductance

    double energy_in;  // this is the fit quantity

    double stoich_b;
    double energy;

    // functions
    int NameToValue(string name);

  public:
    cestructure();
    cestructure(string & str_name, double stoich_b,
        double fit_quantity);
    cestructure(xstructure & xstr);
    virtual ~cestructure();

    // overloaded operator methods
    friend ostream & operator<<(ostream & os, const cestructure & cestr);

    // functions
    vector<double> GetCorrelation(vector<int> & index_list); // get the correlations

    void SetUp(const ceallclusters & cluster1, ceECIcluster & cluster2);
    void SetUp(const ceallclusters & cluster1, ceECIcluster & cluster2,
        istream & corfilein);

    void SetStrCluster(const ceallclusters & ceallcluster_in);
    void SetStoichB(double stoich_b_in){stoich_b = stoich_b_in;};

    ceallclusters GetStrCluster() {cerr << "calculated? " << str_cluster_calculated << endl; return *str_cluster;};

    virtual void GetAllECICorrelation();
    //void GetAllECICorrelation();
    void GetECICorrelation();
    void SetAllECICluster(ceECIcluster & subcluster);
    void SetECICluster(ceECIcluster & subcluster);
    void GetAllECICorrelation(istream & filein);
    int GetAllECICorrelation(istream & filein, int pos);

    void PrintOutCorrelation();
    void WriteFile(ostream & os);
    bool ReadIn(istream & fin);
    int ReadIn(istream & fin, int pos);

    vector<int> AllECICluster(){return allECI_cluster;};
    vector<int> ECICluster(){return ECI_cluster;};
    vector<int> ECIEquivalentNum(){return ECI_equivalent_num;};

    void SetECI(vector<double> ECI_in){ECI=ECI_in;};
    void SetChiSQ(double chisq_in){chisq=chisq_in;};

    vector<double> ECICorrelation() {return ECI_correlation;};
    vector<double> ECIValue(){return ECI;};

    string Name() {return name;}
    double StoichB(){ return stoich_b;};
    double EnergyIn(){ return energy_in;};
    double Energy(){ return energy;};
    void GetEnergy();
    void PrintOutComparison(ostream & os); // print the cacluated and input (formation) eneragies

    cexstructure Structure(){return structure;};
};


// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ceralloy.h

// randow alloy class derived from structure class

class ceralloy : public cestructure
{
  public:
    ceralloy();
    //ceralloy(string & str_name, double stoich_b, double energy, double formation_energy,
    //        double entropy, double volumn, double temperature=0.0);
    //ceralloy(string & str_name, double stoich_b, double temperature=0.0);
    ceralloy(string & str_name, double stoich_b);
    ceralloy(xstructure & xstr);
    ~ceralloy();

    //virtual void GetECICorrelation();
    //virtual void GetAllECICorrelation();
    void GetAllECICorrelation();

    void SetUp(ceallclusters & ceallcluster_in, ceECIcluster & ceECIcluster);

    friend ostream & operator<<(ostream & os, const ceralloy & cestr);
};


// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// cesubcluster.h

struct Baker_str
{
  int index; // cluster index BEGINING FROM 1!!
  int n_num;
  int m_num;

};

class cesubcluster
{
  // subclusters of a base cluster
  private:

    vector<string> name;
    vector<int> base_cluster; // index of base cluster
    vector<Baker_str> overlap_cluster; // overlapped subclusters
    vector< vector<Baker_str> > sub_cluster; // subclusters of all base cluster

    const ceallclusters * str_cluster; // all clusters of a structure
    bool str_cluster_calculated;

    bool overlap_cluster_calculated;

    vector<Baker_str> GetSubCluster(int base_index);

  public:
    cesubcluster();
    //cesubcluster(string & cluster_name);
    ~cesubcluster();

    void SetStrCluster(const ceallclusters & ceallcluster_in);
    void SetStrCluster(const ceallclusters * ceallcluster_in);
    ceallclusters GetStrCluster() {cerr << "calculated? " << str_cluster_calculated << endl; return *str_cluster;};

    void SetBaseCluster(int index);
    void SetBaseCluster(vector<int> index);

    void AddBaseCluster(int index);
    void AddBaseCluster(vector<int> index);

    void GetAllSubCluster();

    void GetOverlapCluster();

    vector<Baker_str> GetOverlapCluster(int index_in);

    vector<Baker_str> GetOverlapCluster(
        const vector<Baker_str> & subcluster_list1,
        const vector<Baker_str> & subcluster_list2);

    vector<Baker_str> MergeTwoSubCluster(
        const vector<Baker_str> & subcluster_list1,
        const vector<Baker_str> & subcluster_list2);

    vector<int> BaseCluster() {return base_cluster;};
    vector<Baker_str> OverlapCluster() {return overlap_cluster;};
    vector< vector<Baker_str> > SubCluster() {return sub_cluster;};

    void SetUp(vector<int> base_cluster, ceallclusters & cluster1);
};


bool is_equal(const Baker_str & cluster1, const Baker_str & cluster2);
bool comparison_subcluster(const Baker_str& subcluster1, const Baker_str& subcluster2);


// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ceSL.h

// superlattice class derived from cestructure class

struct SQScluster
{
  string pair_name;
  int site_num;
  int index;
  double correlation;
  double correlation_random_structure;
  double weight; // weight in the figure of meritof SQScluster
};

class ceSL : public cestructure
{
  // superlattice lattice vectors are defined as
  // r1 = N1 * r01, r2 = N2 * r02, r3 = N3 * r03
  // name is defined as bSL[N1]-[N2]-[N3]_0x[xxxx][hexadecimal base]
  // b stands for binary. the number[xxxx] is in the hexadecimal base
  //
  private:
    vector< vector<int> > N;
    int cell_nr;

    bool sqs_calculated;
    bool is_sqs; // sps or not
    string base;
    vector<int> config;

    deque<_atom> base_atoms; // base atoms in base structures

    // SQS
    // the smaller the is figure of merit, the better is the chance that
    // the structure is an SQS
    double SQS_figure_of_merit;
    vector<SQScluster> SQScompare_cluster_list;

  public:
    ceSL();
    ceSL(int N1, int N2, int N3, string & config);
    ceSL(string & str_name);
    ceSL(string & str_name, double stoich_b,
        double fit_quantity);
    ceSL(string & str_name, xstructure & xstr);
    ~ceSL();

    friend ostream & operator<<(ostream & os, const ceSL & cestr);

    bool NameToStructure(string & SLname);
    string BaseStructure(){return base;};

    string Name() {return name;};
    vector< vector<int> > LatticeVectors() {return N;};

    // set the lattice vector and base atoms
    //void GetStructure(int N1, int N2, int N3);
    void GetStructure(xmatrix<int> Nmat_in, int atom_nr);

    // atom_config specifies the positions of B type atom
    void SetConfig(vector<int> & atom_config);

    void SetName(string _name) { name = _name;};

    void StructureToName();


    vector<int> Configure(){return config;};

    double IsSQS(int site_num, int NNNum);
    void OutputSQS(ostream &oss);
    vector<SQScluster> SQScompareClusterList() {
      return SQScompare_cluster_list;
    }

    //void SetUp(int N1_in, int N2_in, int N3_in, vector<int> & atom_config);
    void SetUp(xmatrix<int> N_in, int cell_nr, vector<int> & atom_config);
    void SetUp(string & SLname);

    void PrintStructure(ostream & os);
    void PrintStructure(ostream & os, vector<_ceatom> atom_species );

    int AtomNr(){return structure.atoms.size();};
    int CellNr(){return cell_nr;};

    //xstructure Structure(){return structure;};
    cexstructure Structure(){return structure;};

    bool IsPrimitiveCell();

    //virtual void GetECICorrelation();
};

string IntToHexString(int from);
vector<int> GetDivisors(int num);
bool comparison_atomtype(_atom atom1, _atom atom2);

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// 
#endif

// ***************************************************************************
// *                                                                         *
// *               AFlow SHIDONG WANG - Duke University 2010-2011            *
// *                                                                         *
// ***************************************************************************
