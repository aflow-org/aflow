// ***************************************************************************
// *                                                                         *
// *               AFlow SHIDONG WANG - Duke University 2010-2011            *
// *                                                                         *
// ***************************************************************************
// aflow_contrib_shidong.cpp
// functions written by
// 2010-2011: shidong.wang@duke.edu

#ifndef _AFLOW_CONTRIB_SHIDONG_CLUSTER_EXPANSION_CPP_
#define _AFLOW_CONTRIB_SHIDONG_CLUSTER_EXPANSION_CPP_

#include "aflow_contrib_shidong_cluster_expansion.h"

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// cexstructure.cpp

//////////////////////////////////////////////////////////////
// cecstructure class -- cluster structure
//////////////////////////////////////////////////////////////

cecstructure::cecstructure() {
  atoms.clear();
  indices.clear();
}
cecstructure::~cecstructure() {
  atoms.clear();
  indices.clear();
}

cecstructure::cecstructure(const cecstructure & structure_in) {
  atoms = structure_in.atoms;
  indices  = structure_in.indices;
}

cecstructure & cecstructure::operator=(const cecstructure & structure_in) {
  atoms = structure_in.atoms;
  indices  = structure_in.indices;
  return *this;
}

//////////////////////////////////////////////////////////////
// cexstructure class -- crystal structure
//////////////////////////////////////////////////////////////

cexstructure::cexstructure() : cecstructure() {
}

cexstructure::cexstructure(xstructure & xstr) {
  SetStructure(xstr);
}

cexstructure::~cexstructure() {
}

void cexstructure::SetStructure(xstructure & xstr) {
  xstr.scale = 1.0;
  xstr.FixLattices();
  atoms = xstr.atoms;
  lattice = xstr.lattice;
  klattice = xstr.klattice;
  c2f = xstr.c2f;
  f2c = xstr.f2c;
  prototype = xstr.prototype;
  num_each_type = xstr.num_each_type;

}

xstructure cexstructure::SetXstructure() {
  xstructure xstr;
  xstr.lattice = lattice;
  xstr.atoms = atoms;
  xstr.prototype = prototype;
  xstr.num_each_type = num_each_type;
  return xstr;
}


// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// cecluster.cpp

//////////////////////////////////////////////////////////////
// cecluster class
//////////////////////////////////////////////////////////////

cecluster::cecluster() : structure() {
  // default constructor
  pair_name ="";
  equivalent_num = 0;
  site_num = 0;
  NNNum = 0;
  index = 0;
}

cecluster::cecluster(const cecluster & cecluster_in) {
  // copy constructor
  pair_name = cecluster_in.pair_name;
  equivalent_num = cecluster_in.equivalent_num;
  site_num = cecluster_in.site_num;
  NNNum = cecluster_in.NNNum;
  dist = cecluster_in.dist;
  index = cecluster_in.index;
  structure = cecluster_in.structure;
}

cecluster::~cecluster() {
  // default destructor
}

cecluster & cecluster::operator=(const cecluster & cecluster_in) {
  // copy constructor
  pair_name = cecluster_in.pair_name;
  equivalent_num = cecluster_in.equivalent_num;
  site_num = cecluster_in.site_num;
  NNNum = cecluster_in.NNNum;
  dist = cecluster_in.dist;
  index = cecluster_in.index;
  structure = cecluster_in.structure;
  return *this;
}

int cecluster::GetNNNum(vector<double> NN_distance) {
  // get the largest neareast neighbor shell index in a cluster

  vector< vector<int> > str;
  const int atom_num = 2;
  double dist, dist_max = -1.0;

  if(structure.atoms.size() == 1) return 1; // only one atom in the cluster


  int total_num = structure.atoms.size();
  int str_size = CombinationNr(atom_num, total_num);
  vector<int> atom_config;
  for(int i=0; i<str_size; i++) {
    if(i == 0) {
      atom_config = AllCombination41(atom_num, total_num, i+1);
    } else {
      atom_config = AllCombination42(atom_num, total_num, atom_config);
    }
    dist = modulus(structure.atoms.at(atom_config.at(0)).cpos
        - structure.atoms.at(atom_config.at(1)).cpos);
    dist_max = max(dist, dist_max); // get the largest distance in the cluster
  }

  for(uint i=0; i<NN_distance.size(); i++) {
    if(is_equal(dist_max, NN_distance.at(i))) {
      return i;
    }
  }

  return -1; // not found in the NN_distance

}

double cecluster::ClusterDistance() {
  // get the total distance between any two atoms in a cluster

  double dist;
  vector< vector<int> > str; // all possible combinations of two atoms
  int total_num, num;


  num = 2; // two atoms involved
  total_num = structure.atoms.size();


  if(total_num == 1) {
    dist = modulus(structure.atoms.at(0).cpos);
  } else {

    int str_size = CombinationNr(num, total_num);
    dist = 0;
    vector<int> atom_config;
    for(int i=0; i<str_size; i++) {
      if(i==0) {
        atom_config = AllCombination41(num, total_num, i+1);
      } else {
        atom_config = AllCombination42(num, total_num, atom_config);
      }
      dist = dist + modulus(structure.atoms.at(atom_config.at(0)).cpos
          - structure.atoms.at(atom_config.at(1)).cpos);
    }



  }

  return dist;
}

//////////////////////////////////////////////////////////////

vector<int> AllCombination41(int num, int total_num, int index) {
  // Knuth's algorithm R

  int c[num+1];

  for(int i=0; i<num; i++) {
    c[i] = i;
  }
  c[num] = total_num;

  vector<int> str;

  long int count = 0;
  int j;


R2:
  ++count;
  if(count == index) {
    str.clear();
    for(int i=0; i< num; i++) {
      //str.push_back(c[i]);
      str.push_back(c[num-i-1]);
    }
    return str;
  }

  if(num%2) {
    if(c[0] + 1 < c[1]) {
      ++c[0];
      goto R2;
    } else {
      j = 2;
      goto R4;
    }
  } else {
    if(c[0] > 0) {
      --c[0];
      goto R2;
    } else {
      j = 2;
      goto R5;
    }
  }

R4:
  if(c[j-1] >= j) {
    c[j-1] = c[j-2];
    c[j-2] = j-2;
    goto R2;
  } else {
    ++j;
  }

R5:
  if(c[j-1] + 1 < c[j]) {
    c[j-2] = c[j-1];
    ++c[j-1];
    goto R2;
  } else {
    ++j;
    if(j <= num) {
      goto R4;
    } else {
      return str;
    }
  }

}

vector<int> AllCombination42(int num, int total_num, vector<int> & str_in) {
  // Knuth's algorithm R

  int c[num+1];

  //for(int i=0; i<num; i++) {
  //    c[i] = i;
  //}

  if(int(str_in.size()) != num) {
    throw aurostd::xerror(__AFLOW_FILE__,"AllCombination42():","input str must have a size equal to "+aurostd::utype2string(num),_INPUT_ILLEGAL_); //CO20200624
  }

  for(int i=0; i<num; i++) {
    c[i] = str_in.at(num-i-1);
  }
  c[num] = total_num;

  vector<int> str;

  long int count = 0;
  int j;
  int next = 2; // only find the next combination


R2:
  ++count;
  if(count == next) {
    str.clear();
    for(int i=0; i< num; i++) {
      //str.push_back(c[i]);
      str.push_back(c[num-i-1]);
    }
    return str;
  }

  if(num%2) {
    if(c[0] + 1 < c[1]) {
      ++c[0];
      goto R2;
    } else {
      j = 2;
      goto R4;
    }
  } else {
    if(c[0] > 0) {
      --c[0];
      goto R2;
    } else {
      j = 2;
      goto R5;
    }
  }

R4:
  if(c[j-1] >= j) {
    c[j-1] = c[j-2];
    c[j-2] = j-2;
    goto R2;
  } else {
    ++j;
  }

R5:
  if(c[j-1] + 1 < c[j]) {
    c[j-2] = c[j-1];
    ++c[j-1];
    goto R2;
  } else {
    ++j;
    if(j <= num) {
      goto R4;
    } else {
      return str;
    }
  }

}

unsigned long int CombinationNr(int num, int total_num) {
  // number of combinations of "num" out of "total_num"
  // needed to add some warning when the combination number is too large
  // to be held in int format

  unsigned long int comNr=1;

  for(int i=0; i<num; i++) {
    comNr *= (total_num-i);
  }
  for(int i=0; i<num; i++) {
    comNr /= (num-i);
  }

  return comNr;
}

// some special clusters used in CVM
cecluster Octahedron() {
  _atom atom_tmp;
  cecluster cluster1;

  atom_tmp.cpos(1) = 0.0;
  atom_tmp.cpos(2) = 0.0;
  atom_tmp.cpos(3) = 0.0;
  cluster1.structure.atoms.push_back(atom_tmp);

  atom_tmp.cpos(1) = 0.5;
  atom_tmp.cpos(2) = 0.0;
  atom_tmp.cpos(3) = 0.5;
  cluster1.structure.atoms.push_back(atom_tmp);

  atom_tmp.cpos(1) = -0.5;
  atom_tmp.cpos(2) = 0.0;
  atom_tmp.cpos(3) = 0.5;
  cluster1.structure.atoms.push_back(atom_tmp);

  atom_tmp.cpos(1) = 0.0;
  atom_tmp.cpos(2) = 0.5;
  atom_tmp.cpos(3) = 0.5;
  cluster1.structure.atoms.push_back(atom_tmp);

  atom_tmp.cpos(1) = 0.0;
  atom_tmp.cpos(2) = -0.5;
  atom_tmp.cpos(3) = 0.5;
  cluster1.structure.atoms.push_back(atom_tmp);

  atom_tmp.cpos(1) = 0.0;
  atom_tmp.cpos(2) = 0.0;
  atom_tmp.cpos(3) = 1.0;
  cluster1.structure.atoms.push_back(atom_tmp);

  return cluster1;
}

cecluster DoubleTetrahedron() {
  _atom atom_tmp;
  cecluster cluster1;

  atom_tmp.cpos(1) = 0.0;
  atom_tmp.cpos(2) = 0.0;
  atom_tmp.cpos(3) = 0.0;
  cluster1.structure.atoms.push_back(atom_tmp);

  atom_tmp.cpos(1) = -0.5;
  atom_tmp.cpos(2) = -0.5;
  atom_tmp.cpos(3) = 0.0;
  cluster1.structure.atoms.push_back(atom_tmp);

  atom_tmp.cpos(1) = -0.5;
  atom_tmp.cpos(2) = 0.5;
  atom_tmp.cpos(3) = 0.0;
  cluster1.structure.atoms.push_back(atom_tmp);

  atom_tmp.cpos(1) = 0.0;
  atom_tmp.cpos(2) = -0.5;
  atom_tmp.cpos(3) = 0.5;
  cluster1.structure.atoms.push_back(atom_tmp);

  atom_tmp.cpos(1) = 0.0;
  atom_tmp.cpos(2) = 0.5;
  atom_tmp.cpos(3) = 0.5;
  cluster1.structure.atoms.push_back(atom_tmp);

  atom_tmp.cpos(1) = -0.5;
  atom_tmp.cpos(2) = 0.0;
  atom_tmp.cpos(3) = 0.5;
  cluster1.structure.atoms.push_back(atom_tmp);

  return cluster1;
}

cecluster Tetrahedron() {
  _atom atom_tmp;
  cecluster cluster1;

  atom_tmp.cpos(1) = 0.0;
  atom_tmp.cpos(2) = 0.0;
  atom_tmp.cpos(3) = 0.0;
  cluster1.structure.atoms.push_back(atom_tmp);

  atom_tmp.cpos(1) = -0.5;
  atom_tmp.cpos(2) = 0.5;
  atom_tmp.cpos(3) = 0.0;
  cluster1.structure.atoms.push_back(atom_tmp);

  atom_tmp.cpos(1) = -0.5;
  atom_tmp.cpos(2) = 0.0;
  atom_tmp.cpos(3) = 0.5;
  cluster1.structure.atoms.push_back(atom_tmp);

  atom_tmp.cpos(1) = 0.0;
  atom_tmp.cpos(2) = 0.5;
  atom_tmp.cpos(3) = 0.5;
  cluster1.structure.atoms.push_back(atom_tmp);

  return cluster1;
}

cecluster Triplet() {
  _atom atom_tmp;
  cecluster cluster1;

  atom_tmp.cpos(1) = 0.0;
  atom_tmp.cpos(2) = 0.0;
  atom_tmp.cpos(3) = 0.0;
  cluster1.structure.atoms.push_back(atom_tmp);

  atom_tmp.cpos(1) = -0.5;
  atom_tmp.cpos(2) = 0.0;
  atom_tmp.cpos(3) = 0.5;
  cluster1.structure.atoms.push_back(atom_tmp);

  atom_tmp.cpos(1) = 0.0;
  atom_tmp.cpos(2) = 0.5;
  atom_tmp.cpos(3) = 0.5;
  cluster1.structure.atoms.push_back(atom_tmp);

  return cluster1;
}

cecluster Pair() {
  _atom atom_tmp;
  cecluster cluster1;

  atom_tmp.cpos(1) = 0.0;
  atom_tmp.cpos(2) = 0.0;
  atom_tmp.cpos(3) = 0.0;
  cluster1.structure.atoms.push_back(atom_tmp);

  atom_tmp.cpos(1) = -0.5;
  atom_tmp.cpos(2) = 0.0;
  atom_tmp.cpos(3) = 0.5;
  cluster1.structure.atoms.push_back(atom_tmp);

  return cluster1;
}

bool comparison_cluster(const cecluster& cluster1, const cecluster& cluster2) {
  // comparison of two clusters
  // "smaller" cluster has "small" pair_name
  // sorted by descending order
  // sorted then by atom position

  int size1 = cluster1.structure.atoms.size();
  int size2 = cluster2.structure.atoms.size();

  if(size1 != size2) {
    return size1 < size2;
  } else {
    if(cluster1.pair_name != cluster2.pair_name) {
      return (cluster1.pair_name < cluster2.pair_name);
    } else {
      double dist1=0.0;
      double dist2=0.0;
      int shell1=-1;
      int shell2=-1;
      _atom atom1, atom2;

      for(int i=0; i<size1; i++) {
        atom1 = cluster1.structure.atoms.at(i);
        atom2 = cluster2.structure.atoms.at(i);

        dist1 += modulus(atom1.cpos);
        dist2 += modulus(atom2.cpos);
        shell1 = max(shell1, atom1.shell);
        shell2 = max(shell2, atom2.shell);
      }

      return (shell1 < shell2) && (dist1 < dist2);
    }
  }

  //return (cluster1.pair_name < cluster2.pair_name);
}

//////////////////////////////////////////////////////////////
// namespace functions
//////////////////////////////////////////////////////////////

bool is_equal(const double f1, const double f2) {
  return (abs(f1 - f2) < _EQUAL_DOUBLE);
}
bool is_equal(const _atom& atom1, const _atom& atom2) {
  return (abs(modulus(atom1.cpos-atom2.cpos)) < _EQUAL_ATOM);
}
bool is_equal_fpos(const _atom& atom1, const _atom& atom2) {
  return (abs(modulus(atom1.fpos-atom2.fpos)) < _EQUAL_ATOM);
}
bool is_equal(const cecluster & cluster1, const cecluster & cluster2) {
  bool flag_cluster = true, flag_atom;
  _atom atom_tmp;

  if(cluster1.structure.atoms.size() != cluster2.structure.atoms.size()) return false;

  for(uint i = 0; i<cluster1.structure.atoms.size(); i++) {
    atom_tmp = cluster1.structure.atoms.at(i);
    flag_atom = false;
    for(uint j = 0; j<cluster2.structure.atoms.size(); j++) {
      flag_atom = flag_atom || is_equal(atom_tmp, cluster2.structure.atoms.at(j));
    }
    flag_cluster = flag_cluster && flag_atom;
  }

  return flag_cluster;
}

bool is_equal_crystal(const cecluster & cluster1,
    const cecluster & cluster2) {
  // two clusters are equivalent if they are related by a translational operator
  _atom atom_tmp, atom_ref;
  cecluster cluster_tmp;

  if(cluster1.structure.atoms.size() != cluster2.structure.atoms.size()) return false;

  for(uint i=0; i<cluster1.structure.atoms.size(); i++) {
    atom_ref = cluster1.structure.atoms.at(i);
    cluster_tmp.structure.atoms.clear();
    for(uint j=0; j<cluster1.structure.atoms.size(); j++) {
      atom_tmp  = cluster1.structure.atoms.at(j);
      atom_tmp.cpos = atom_tmp.cpos - atom_ref.cpos;
      cluster_tmp.structure.atoms.push_back(atom_tmp);
    }

    if(is_equal(cluster_tmp, cluster2)) return true;
  }

  return false;
}



// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ceallcluster.cpp

//////////////////////////////////////////////////////////////
// ceallclusters class
//////////////////////////////////////////////////////////////

ceallclusters::ceallclusters::ceallclusters() : structure() {
  // default constructor
  name ="";
  rep_cluster.clear();
  all_cluster.clear();
  num_rep_cluster = 0;
  num_total_cluster = 0;
  NN_distance.clear();
  NN_shell_num.clear();
}

ceallclusters::ceallclusters(string & name_in) : structure() {
  name = name_in;
  rep_cluster.clear();
  all_cluster.clear();
  num_rep_cluster = 0;
  num_total_cluster = 0;
  NN_distance.clear();
  NN_shell_num.clear();

  // get the strucutre
  structure=aflowlib::PrototypeLibraries(cerr, name,"",LIBRARY_MODE_HTQC); // no parameters
  structure.FixLattices();
  structure.CalculateSymmetryPointGroupCrystal(); // only the crystal symmetry is needed

  cerr << "title " << structure.title << endl;
  cerr << "proto type " << structure.prototype << endl;
  cerr << "info " << structure.info << endl;
  cerr << "coord_type " << structure.coord_type << endl;
  cerr << "atom number " << structure.atoms.size() << endl;

  cerr << "pgroup \n";
  cerr << "pgroup size " << structure.pgroup.size() << endl;
  cerr << "pgroup calculated " << structure.pgroup_calculated << endl;
  cerr << "crystal family " << structure.crystal_family << endl;
  cerr << "crystal system " << structure.crystal_system << endl;
  cerr << "point group class " << structure.point_group_crystal_class << endl;
  cerr << "point group structure " << structure.point_group_structure << endl;
  cerr << "space group structure " << structure.spacegroup << endl;
  cerr << endl;

}

ceallclusters::ceallclusters(const ceallclusters & allcluster_in) : structure(allcluster_in.structure) {
  name = allcluster_in.name;
  rep_cluster = allcluster_in.rep_cluster;
  all_cluster = allcluster_in.all_cluster;
  num_rep_cluster = allcluster_in.num_rep_cluster;
  num_total_cluster = allcluster_in.num_total_cluster;
  NN_distance = allcluster_in.NN_distance;
  NN_shell_num = allcluster_in.NN_shell_num;

  // get the strucutre
  //structure=aflowlib::PrototypeLibraries(cerr, name,"",LIBRARY_MODE_HTQC); // no parameters;
  //structure.FixLattices();
  //structure.CalculateSymmetryPointGroupCrystal(); // only the crystal symmetry is needed

  //cerr << "From copy constructor\n";
  //cerr << "title " << structure.title << endl;
  //cerr << "proto type " << structure.prototype << endl;
  //cerr << "info " << structure.info << endl;
  //cerr << "coord_type " << structure.coord_type << endl;
  //cerr << "atom number " << structure.atoms.size() << endl;

  //cerr << "pgroup \n";
  //cerr << "pgroup size " << structure.pgroup.size() << endl;
  //cerr << "pgroup calculated " << structure.pgroup_calculated << endl;
  //cerr << "crystal family " << structure.crystal_family << endl;
  //cerr << "crystal system " << structure.crystal_system << endl;
  //cerr << "point group class " << structure.point_group_crystal_class << endl;
  //cerr << "point group structure " << structure.point_group_structure << endl;
  //cerr << "space group structure " << structure.spacegroup << endl;
  //cerr << endl;

}

ceallclusters & ceallclusters::operator=(const ceallclusters & allcluster_in) {
  structure=allcluster_in.structure;
  name = allcluster_in.name;
  rep_cluster = allcluster_in.rep_cluster;
  all_cluster = allcluster_in.all_cluster;
  num_rep_cluster = allcluster_in.num_rep_cluster;
  num_total_cluster = allcluster_in.num_total_cluster;
  NN_distance = allcluster_in.NN_distance;
  NN_shell_num = allcluster_in.NN_shell_num;

  // get the strucutre
  //structure=aflowlib::PrototypeLibraries(cerr, name,"",LIBRARY_MODE_HTQC); // no parameters;
  //structure.FixLattices();
  //structure.CalculateSymmetryPointGroupCrystal(); // only the crystal symmetry is needed

  //cerr << "From operator = \n";
  //cerr << "title " << structure.title << endl;
  //cerr << "proto type " << structure.prototype << endl;
  //cerr << "info " << structure.info << endl;
  //cerr << "coord_type " << structure.coord_type << endl;
  //cerr << "atom number " << structure.atoms.size() << endl;

  //cerr << "pgroup \n";
  //cerr << "pgroup size " << structure.pgroup.size() << endl;
  //cerr << "pgroup calculated " << structure.pgroup_calculated << endl;
  //cerr << "crystal family " << structure.crystal_family << endl;
  //cerr << "crystal system " << structure.crystal_system << endl;
  //cerr << "point group class " << structure.point_group_crystal_class << endl;
  //cerr << "point group structure " << structure.point_group_structure << endl;
  //cerr << "space group structure " << structure.spacegroup << endl;
  //cerr << endl;

  return *this;

}

ceallclusters::~ceallclusters() {
  // default destructor
  rep_cluster.clear();
  all_cluster.clear();
  NN_distance.clear();
  NN_shell_num.clear();
}

ostream & operator<<(ostream & os, const ceallclusters & cestr) {

  os << "structure name " << cestr.name << endl;
  os << cestr.structure << endl;

  os.precision(6);
  os.unsetf(ios_base::floatfield);

  os << "    shell      distance       num\n";
  for(uint i = 0; i < cestr.NN_distance.size(); i++) {
    //cerr.setf(ios_base::fixed);
    //cerr.precision(4);
    os << setw(6)
      << i
      << setw(16)
      << cestr.NN_distance.at(i)
      << setw(10)
      << cestr.NN_shell_num.at(i)
      << endl;
  }

  return os;
}

int ceallclusters::SetCluster( const int& SiteNum_low, const int& SiteNum_up,
    const int& NNNum_low, const int& NNNum_up) {

  clock_t start, end;
  double diff;

  //////////////////////////////////////////////////////////////
  start = clock();
  //////////////////////////////////////////////////////////////

  GetNearestNeighbor(NNNum_up);
  cerr << "Number of atom " << atom_list.size() << endl;

  //////////////////////////////////////////////////////////////
  end = clock();
  diff = (end - start)/ CLOCKS_PER_SEC;
  cerr << "GetNearestNeighbor Time: " << diff << " s\n";
  //////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////
  start = clock();
  //////////////////////////////////////////////////////////////

  GenerateAllCluster(SiteNum_low, SiteNum_up, NNNum_low, NNNum_up);

  //////////////////////////////////////////////////////////////
  end = clock();
  diff = (end - start)/ CLOCKS_PER_SEC;
  cerr << "GetAllSubcluster Time: " << diff << " s\n";
  //////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////
  start = clock();
  //////////////////////////////////////////////////////////////

  GetRepresentCluster();

  //////////////////////////////////////////////////////////////
  end = clock();
  diff = (end - start)/ CLOCKS_PER_SEC;
  cerr << "GetRepresentCluster Time: " << diff << " s\n";
  //////////////////////////////////////////////////////////////


  PrintRepCluster();
  PrintAllCluster();

  return num_rep_cluster;
}

int ceallclusters::SetCluster( string & filename) {

  clock_t start, end;
  double diff;

  //////////////////////////////////////////////////////////////
  start = clock();
  //////////////////////////////////////////////////////////////

  ReadIn(filename);

  cerr << "Number of atom " << atom_list.size() << endl;
  //////////////////////////////////////////////////////////////
  end = clock();
  diff = (end - start)/ CLOCKS_PER_SEC;
  cerr << "GetRepresentCluster Time: " << diff << " s\n";
  //////////////////////////////////////////////////////////////


  //PrintRepCluster();
  //PrintAllCluster();

  return num_rep_cluster;
}

int ceallclusters::GetNearestNeighbor(int NNNum) {
  // get the neareast neighbor distances up to NNNum-th order
  //
  // The nearest neighbor distance is calculated by populating atoms
  // next to the one in origin
  // The maximum index number needed to get all nearest neighbor distance lesser
  // lesser than NNNum is sqrt(NNNum*2)
  // see http://www.chem.lsu.edu/htdocs/people/sfwatkins/MERLOT/cubic_neighbors/cubic_near_neighbors.html
  // for details

  vector<double> NN_distance_all;
  _atom atom_tmp, atom_orig = structure.atoms.at(0);
  int grid_num = sqrt(2.0*double(NNNum));


  // more than one atoms in the unit cells
  NN_distance_all.clear();
  NN_distance.clear();
  for(uint i = 0; i < structure.atoms.size(); i++) {
    NN_distance_all.push_back(modulus(structure.atoms.at(i).cpos-atom_orig.cpos));
  }

  structure.GenerateGridAtoms(grid_num);

  //    cout.setf(ios_base::fixed);
  //    cout.precision(10);
  for(uint i = 0; i < structure.grid_atoms.size(); i++) {
    NN_distance_all.push_back(modulus(structure.grid_atoms.at(i).cpos-atom_orig.cpos));
  }

  // sort and delete the duplicate values
  sort(NN_distance_all.begin(), NN_distance_all.end());

  int j=0;
  NN_distance.push_back(NN_distance_all.at(0));
  for(uint i = 0; i< NN_distance_all.size(); i++) { // delete duplicate items
    if(!is_equal(NN_distance_all.at(i), NN_distance_all.at(j))) {
      NN_distance.push_back(NN_distance_all.at(i));
      j = i;
    }

  }

  NN_distance.resize(NNNum+1);

  // get all atoms with distance no larger than the largest value in NN_distance
  // by search all possible atoms. may not be efficient.

  atom_list.clear();
  for(uint i = 0; i < NN_distance.size(); i++) {
    NN_shell_num.push_back(0);
    for(uint j = 0; j< structure.grid_atoms.size(); j++) {
      if(is_equal(modulus(structure.grid_atoms.at(j).cpos-atom_orig.cpos),
            NN_distance.at(i))) {
        atom_tmp = structure.grid_atoms.at(j);
        atom_tmp.shell = i;
        atom_list.push_back(atom_tmp);
        NN_shell_num.at(i)++;
        continue;
      }
    }
  }

  return atom_list.size();
}

int ceallclusters::GenerateAllCluster(int SiteNum_low, int SiteNum_up,
    int NNNum_low, int NNNum_up) {
  // generate all cluster in a structure

  int site_index, NN_index;
  cecluster cluster_tmp;
  int NN_up;

  vector< vector<int> > str;


  if(SiteNum_up == 1) {
    cluster_tmp.pair_name = "1 0";
    cluster_tmp.equivalent_num = 1;
    cluster_tmp.NNNum = 0;
    cluster_tmp.site_num = 1;
    cluster_tmp.structure.atoms.push_back(atom_list.at(0));
    cluster_tmp.structure.indices.push_back(0);
    rep_cluster.push_back(cluster_tmp);
  } else {

    for(site_index = SiteNum_low; site_index < SiteNum_up + 1; site_index++) {

      // clean up the cluster_tmp;
      cluster_tmp.pair_name.clear();
      cluster_tmp.structure.atoms.clear();
      cluster_tmp.structure.indices.clear();
      cluster_tmp.equivalent_num = 1;
      cluster_tmp.site_num = site_index;

      if(site_index == 1) {
        // one atom has no nearest neighbor
        cluster_tmp.pair_name = "1 0";
        cluster_tmp.structure.atoms.push_back(atom_list.at(0));
        cluster_tmp.structure.indices.push_back(0);
        rep_cluster.push_back(cluster_tmp);
        continue;
      }

      cluster_tmp.pair_name.append(aurostd::utype2string(site_index));
      cluster_tmp.pair_name.push_back(' ');

      int num_config = site_index-1;
      int total_num = atom_list.size()-1;
      int str_size = CombinationNr(num_config, total_num);

      NN_up = NNNum_up;

      cluster_tmp.structure.atoms.clear();
      cluster_tmp.structure.indices.clear();
      cluster_tmp.equivalent_num = 1;

      vector<int> atom_config;
      for(int i=0; i< str_size; i++) {
        if(i==0) {
          atom_config=AllCombination41(num_config, total_num, i+1);
        } else {
          atom_config=AllCombination42(num_config, total_num, atom_config);
        }

        cluster_tmp.structure.atoms.clear();
        cluster_tmp.structure.indices.clear();
        for(uint j=0; j<atom_config.size(); j++) {
          //cout << atom_config.at(j) << " ";

          cluster_tmp.structure.atoms.push_back(atom_list.at(atom_config.at(j)+1));
          cluster_tmp.structure.indices.push_back(atom_config.at(j)+1);
        }


        cluster_tmp.structure.atoms.push_back(atom_list.at(0)); // include the atom at origin
        cluster_tmp.structure.indices.push_back(0);
        cluster_tmp.NNNum = cluster_tmp.GetNNNum(NN_distance);


        for(NN_index = NNNum_low; NN_index < NN_up + 1; NN_index++) {
          if(cluster_tmp.NNNum == NN_index) {

            cluster_tmp.pair_name.append(aurostd::utype2string(cluster_tmp.NNNum)); // neareset neighbor shell index

            rep_cluster.push_back(cluster_tmp);

            // erase the old NNN_num in pair_name to store the new one in next step
            cluster_tmp.pair_name.erase(cluster_tmp.pair_name.size()-1, cluster_tmp.pair_name.size()-1);
            //break;
          }
        }
      }
    }

  }

  // sort rep_cluster by pair_name values
  sort(rep_cluster.begin(), rep_cluster.end(), &comparison_cluster);

  GetClusterDistance(); // add the third index in pair_name

  sort(rep_cluster.begin(), rep_cluster.end(), &comparison_cluster);

  num_total_cluster = rep_cluster.size();

  return num_total_cluster;
}

int ceallclusters::GetRepresentCluster() {
  // get the representive clusters

  // keep only one of those symmetric atoms
  // divide clusters into equivlant symmetry groups


  vector<cecluster> rep_cluster_tmp;
  _atom atom_tmp, atom_cmp;
  vector< vector<cecluster> > all_cluster_tmp;
  vector<cecluster> rep_cluster_tmp1;


  // compare every two clusters
  // may not be most efficient way to do this
  rep_cluster_tmp.push_back(rep_cluster.at(0));
  rep_cluster_tmp.at(0).equivalent_num--; // avoid double count of (1 0 1) cluster
  rep_cluster_tmp1.clear();
  all_cluster_tmp.push_back(rep_cluster_tmp1);

  int num;

  cerr << "total number of clusters " << rep_cluster.size() << endl;
  num = EquivalentCluster(structure.pgroup_xtal); // crystal group

  cerr << "representive cluster num " << rep_cluster.size() << " num=" << num << endl;

  for(uint i=0; i<all_cluster.size(); i++) {

    for(uint j=0; j<all_cluster.at(i).size(); j++) {
      all_cluster.at(i).at(j).equivalent_num = rep_cluster.at(i).equivalent_num;
      all_cluster.at(i).at(j).dist = rep_cluster.at(i).dist;
      all_cluster.at(i).at(j).pair_name.push_back(' ');
      all_cluster.at(i).at(j).pair_name.append(aurostd::utype2string(i+1));
      all_cluster.at(i).at(j).index = i+1;
    }

    // same pair name forms in both rep_cluster and all_cluster
    rep_cluster.at(i).pair_name.push_back(' ');
    rep_cluster.at(i).pair_name.append(aurostd::utype2string(i+1));
    rep_cluster.at(i).index = i+1;
  }


  num_rep_cluster = rep_cluster.size();
  return num_rep_cluster;

}

void ceallclusters::PrintRepCluster() const
{
  cerr << "Representive clusters \n";
  cerr << "Cluster list size : " << rep_cluster.size() << endl;


  for(uint i = 0; i < rep_cluster.size(); i++) {
    cerr << "Cluster Name (" << rep_cluster.at(i).pair_name << ")\n";
    cerr << "Equivalent clusters " << rep_cluster.at(i).equivalent_num << "\n";
    cerr << "Site_num " << rep_cluster.at(i).site_num
      << " " << "NNNnum " << rep_cluster.at(i).NNNum
      << " " << "dist " << rep_cluster.at(i).dist
      << " " << "index " << rep_cluster.at(i).index << endl;
    for(uint j = 0; j < rep_cluster.at(i).structure.atoms.size();j++) {
      cerr << setw(6)
        << rep_cluster.at(i).structure.atoms.at(j).cpos(1) << " "
        << setw(6)
        << rep_cluster.at(i).structure.atoms.at(j).cpos(2) << " "
        << setw(6)
        << rep_cluster.at(i).structure.atoms.at(j).cpos(3) << "\n";
    }

    // output the indices of atom in atom_list
    cerr << "Indices of these atoms in atom_lists are: \n";
    for(uint j = 0; j < rep_cluster.at(i).structure.indices.size();j++) {
      cerr << setw(6)
        << rep_cluster.at(i).structure.indices.at(j) ;
    }
    cerr << "\n";
  }

}

void ceallclusters::PrintAllCluster() const
{
  cerr << "All clusters \n";
  cerr << "Cluster type number : " << all_cluster.size() << endl;

  int count;

  count = 0;
  for(uint k = 0; k < all_cluster.size(); k++) {
    for(uint i=0; i< all_cluster.at(k).size(); i++) {
      cerr << "Cluster Name (" << all_cluster.at(k).at(i).pair_name << ")\n";
      cerr << "Equivalent clusters " << all_cluster.at(k).at(i).equivalent_num << "\n";
      cerr << "Site_num " << all_cluster.at(k).at(i).site_num
        << " " << "NNNnum " << all_cluster.at(k).at(i).NNNum
        << " " << "dist " << all_cluster.at(k).at(i).dist
        << " " << "index " << all_cluster.at(k).at(i).index << endl;
      for(uint j = 0; j < all_cluster.at(k).at(i).structure.atoms.size();j++) {
        cerr << setw(6)
          << all_cluster.at(k).at(i).structure.atoms.at(j).cpos(1) << " "
          << setw(6)
          << all_cluster.at(k).at(i).structure.atoms.at(j).cpos(2) << " "
          << setw(6)
          << all_cluster.at(k).at(i).structure.atoms.at(j).cpos(3) << "\n";
      }

      // output the indices of atom in atom_list
      cerr << "Indices of these atoms in atom_lists are: \n";
      for(uint j = 0; j < all_cluster.at(k).at(i).structure.indices.size();j++) {
        cerr << setw(6)
          << all_cluster.at(k).at(i).structure.indices.at(j);
      }
      cerr << "\n";

      count++;
    }
  }

  cerr << "Number of clusters outputed is " << count << endl;

}

int ceallclusters::EquivalentCluster(const vector<_sym_op> pgroup) {
  // separate the clusters into sets of equivalent cluster
  // original cluster list is stored in rep_cluster defined in the class
  // output: rep_cluster contains the representive cluster list
  // all_cluster contains equivalent cluster sets

  cecluster cluster1, cluster2;
  bool flag_eq;
  vector<cecluster> rep_cluster_tmp, cluster_eq_list_tmp;
  vector< vector<cecluster> > all_cluster_tmp;

  int lcol = 1;
  //int ucol = cluster2.structure.atoms.size();
  int ucol;

  if(pgroup.size() == 0) {
    throw aurostd::xerror(__AFLOW_FILE__,"ceallclusters::EquivalentCluster():","no symmetric operator",_RUNTIME_ERROR_); //CO20200624
  }

  int lrow = 1;
  int urow = 3; // number of coordinate elements

  vector<int> index_list; // stored cluster index list

  while (rep_cluster.size() != 0) {
    cluster1 = rep_cluster.at(0);
    ucol = cluster1.structure.atoms.size();
    xmatrix<double> coord_atoms(lrow, lcol, urow, ucol);

    _atom atom_tmp;

    // set the cartesian coordinate matrix for all atom
    for(uint i=0; i<cluster1.structure.atoms.size(); i++) {
      atom_tmp = cluster1.structure.atoms.at(i);
      for(uint j=1; j<_DIM+1; j++) {
        coord_atoms[j][i+1] = atom_tmp.cpos[j];
      }
    }


    cluster_eq_list_tmp.clear();

    cluster_eq_list_tmp.push_back(cluster1);

    rep_cluster_tmp.push_back(cluster1); // store the representive cluster

    xmatrix<double> symmetry_op;
    xmatrix<double> coord_sym;
    cecluster cluster_tmp;
    bool flag;

    for(uint i=0; i<pgroup.size(); i++) {
      // set the symmetric operator matrix
      symmetry_op = pgroup.at(i).Uc;

      coord_sym = symmetry_op*coord_atoms;

      cluster_tmp.structure.atoms.clear();

      index_list.clear();
      flag = false;

      // store the trasferred atoms in a temporary cluster object
      for(uint j=0; j<cluster1.structure.atoms.size(); j++) {
        for(uint k=0; k<3; k++) {
          atom_tmp.cpos[k+1] = coord_sym[k+1][j+1];
        }
        cluster_tmp.structure.atoms.push_back(atom_tmp);
      }

      for(uint j=1; j<rep_cluster.size(); j++) {
        flag_eq = false;
        cluster2 = rep_cluster.at(j);

        if(cluster1.structure.atoms.size()
            < cluster2.structure.atoms.size()
            || (cluster1.pair_name != cluster2.pair_name)) {
          if(j== 1) {
            flag = true;
          }
          break;
        }

        if(cluster1.structure.atoms.size()
            == cluster2.structure.atoms.size()) {
          flag_eq = flag_eq || is_equal_crystal(cluster2, cluster_tmp);
        } else {
          flag_eq = false;
        }

        //if(flag_eq) return j; //return true
        if(flag_eq) {
          cluster_eq_list_tmp.push_back(cluster2); // equivalent clusters
          index_list.push_back(j);
        }

      }

      //cerr << "index_list size " << index_list.size() << endl;
      //cerr << "flag " << flag << endl;

      // delete found cluster from rep_cluster
      for(int j=index_list.size()-1; j>=0; j--) {
        rep_cluster.erase(rep_cluster.begin()+index_list.at(j));
      }

      if(flag) {
        break;
      }

    }

    rep_cluster_tmp.back().equivalent_num =
      cluster_eq_list_tmp.size();
    all_cluster_tmp.push_back(cluster_eq_list_tmp); // store

    rep_cluster.erase(rep_cluster.begin());

    //cerr << "rep_cluster size " << rep_cluster.size() << endl;

  }

  rep_cluster.clear();
  all_cluster.clear();
  rep_cluster = rep_cluster_tmp;
  all_cluster = all_cluster_tmp;

  return rep_cluster.size();

}

void ceallclusters::GetClusterDistance() {
  // catalog cluster further by the distance between two atoms in a cluster
  // modified the rep_cluster
  // rep_cluster must be in ascending order with respect to pair_name


  vector<double> distance;
  double dist;
  int counter_same_old_label=0; // cluster with same old label
  string label_current; // current pair_name
  bool flag;
  int pos; // position of the last cluster with same old pair_name

  pos = 0;
  label_current = rep_cluster.at(0).pair_name;
  for(uint i = 0; i < rep_cluster.size(); i++) {

    dist = 0.0;
    if( rep_cluster.at(i).pair_name == label_current && (i != rep_cluster.size()-1)) {
      counter_same_old_label++;

      dist = rep_cluster.at(i).ClusterDistance();
      flag = false;
      for(uint k=0; k<distance.size(); k++) {
        flag = flag || is_equal(dist, distance.at(k));
      }
      if(!flag) distance.push_back(dist);

    } else {
      // assign the third label to pair_name

      if(i == rep_cluster.size() -1)counter_same_old_label++;
      for(int j = pos; j < pos + counter_same_old_label; j++) {
        dist = rep_cluster.at(j).ClusterDistance();
        for(uint k = 0; k < distance.size(); k++) {
          if(is_equal(dist, distance.at(k))) {
            rep_cluster.at(j).pair_name.push_back(' ');
            rep_cluster.at(j).pair_name.append(aurostd::utype2string(k+1));
            rep_cluster.at(j).dist = k+1;
            break;
          }
        }
      }


      label_current = rep_cluster.at(i).pair_name;
      pos = pos + counter_same_old_label;
      distance.clear();

      counter_same_old_label = 1;

      dist = 0;
      dist = rep_cluster.at(i).ClusterDistance();

      distance.push_back(dist);
    }

  }

}

string ceallclusters::GetNamebyCluster(cecluster & cluster_in) const
{
  // get the name of provided structure

  cecluster cluster_cmp;
  bool flag;
  string name;

  flag = false;

  for(uint i=0; i< all_cluster.size(); i++) {
    for(uint j=0; j<all_cluster.at(i).size(); j++) {
      cluster_cmp = all_cluster.at(i).at(j);
      if(cluster_cmp.structure.atoms.size()
          != cluster_in.structure.atoms.size()) {
        flag = false;
        break;
      } else {
        if(is_equal_crystal(cluster_in, cluster_cmp)) {
          name = cluster_cmp.pair_name;
          flag = true;
          break;
        }
      }
    }
    if(flag) break;
  }

  return name;
}

int ceallclusters::GetIndexbyCluster(cecluster & cluster_in) const
{
  // get the name of provided structure

  cecluster cluster_cmp;
  bool flag;
  int index=-1;

  flag = false;

  for(uint i=0; i< all_cluster.size(); i++) {
    for(uint j=0; j<all_cluster.at(i).size(); j++) {
      cluster_cmp = all_cluster.at(i).at(j);
      if(cluster_cmp.structure.atoms.size()
          != cluster_in.structure.atoms.size()) {
        flag = false;
        break;
      } else {
        if(is_equal_crystal(cluster_in, cluster_cmp)) {
          index = cluster_cmp.index;
          flag = true;
          break;
        }
      }
    }
    if(flag) break;
  }

  return index;
}

bool ceallclusters::ReadIn(string & filename) {
  // read data (rep_cluster and all_cluster from a file

  ifstream fin;

  fin.open(filename.c_str());

  if(!fin) { // cannot open file
    throw aurostd::xerror(__AFLOW_FILE__,"ceallclusters::ReadIn():","file "+filename+" cannot be opened",_FILE_CORRUPT_); //CO20200624
  }

  int rep_cluster_num;
  fin >> rep_cluster_num;


  for(int i = 0; i< rep_cluster_num; i++) {
    cecluster cluster1;
    vector<cecluster> cluster_list;

    fin.ignore(1,'\n');
    getline(fin, cluster1.pair_name, '\n');

    fin >> cluster1.equivalent_num;
    fin >> cluster1.site_num
      >> cluster1.NNNum
      >> cluster1.dist
      >> cluster1.index;

    _atom atom1;
    int atom_index;

    for(int i1=0; i1<cluster1.site_num; i1++) {
      fin >> atom1.cpos(1)
        >> atom1.cpos(2)
        >> atom1.cpos(3)
        >> atom_index;
      cluster1.structure.atoms.push_back(atom1);
      cluster1.structure.indices.push_back(atom_index);
    }

    int equivalent_num;

    if(cluster1.dist <= _DIST_LIST[cluster1.site_num-1]) {
      rep_cluster.push_back(cluster1);
    }

    equivalent_num = cluster1.equivalent_num;

    cluster1.structure.atoms.clear();
    cluster1.structure.indices.clear();
    for(int j=0; j< equivalent_num; j++) {
      // read in all_cluster

      fin.ignore(1,'\n');
      getline(fin, cluster1.pair_name);
      fin >> cluster1.equivalent_num;
      fin >> cluster1.site_num
        >> cluster1.NNNum
        >> cluster1.dist
        >> cluster1.index;

      cluster1.structure.atoms.clear();
      cluster1.structure.indices.clear();
      for(int i1=0; i1<cluster1.site_num; i1++) {
        fin >> atom1.cpos(1)
          >> atom1.cpos(2)
          >> atom1.cpos(3)
          >> atom_index;
        cluster1.structure.atoms.push_back(atom1);
        cluster1.structure.indices.push_back(atom_index);
      }

      if(cluster1.dist <= _DIST_LIST[cluster1.site_num-1]) {
        cluster_list.push_back(cluster1);
      }
    }

    if(cluster1.dist <= _DIST_LIST[cluster1.site_num-1]) {
      all_cluster.push_back(cluster_list);
    }

  }

  // read NN_distance, NN_shell_num, and atom_list

  int NN_distance_num;
  NN_distance.clear();
  NN_shell_num.clear();

  fin >> NN_distance_num;
  for(int i=0; i< NN_distance_num; i++) {
    double dist;
    int  shell;
    fin >> dist >> shell;
    NN_distance.push_back(dist);
    NN_shell_num.push_back(shell);
  }

  int atom_list_num;
  atom_list.clear();
  fin >> atom_list_num;
  for(int i=0; i<atom_list_num; i++) {
    _atom atom1;
    fin >> atom1.cpos[1]
      >> atom1.cpos[2]
      >> atom1.cpos[3];
    atom_list.push_back(atom1);
  }

  fin.close();

  sort(rep_cluster.begin(), rep_cluster.end(), &comparison_cluster);
  sort(all_cluster.begin(), all_cluster.end(), &comparison_cluster_list);

  // reassign the last index
  for(uint i=0; i<rep_cluster.size(); i++) {
    string pair_name;
    pair_name = rep_cluster.at(i).pair_name;
    while (pair_name.at(pair_name.size()-1) != ' ') {
      pair_name.erase(pair_name.end()-1);
    }
    pair_name.append(aurostd::utype2string(i+1));
    rep_cluster.at(i).pair_name = pair_name;
    rep_cluster.at(i).index = i+1;

    for(uint j=0; j<all_cluster.at(i).size(); j++) {
      all_cluster.at(i).at(j).pair_name = pair_name;
      all_cluster.at(i).at(j).index = i+1;
    }


  }


  return true;
}

bool ceallclusters::WriteFile(string & filename, string & stat) {
  // write rep_cluster and all_cluster to a file

  ofstream fout;

  cerr << "stat " << stat << endl;
  if(stat == "new") {
    // overwrite the file
    fout.open(filename.c_str()); // default mode
  } else {
    // append to the file
    fout.open(filename.c_str(), ios_base::out | ios_base::app);
  }

  fout << rep_cluster.size() << endl;
  fout.precision(6);
  for(uint i=0; i<rep_cluster.size(); i++) {
    fout << rep_cluster.at(i).pair_name << endl;
    fout << rep_cluster.at(i).equivalent_num << endl;
    fout << setw(8)
      << rep_cluster.at(i).site_num
      << setw(8)
      << rep_cluster.at(i).NNNum
      << setw(8)
      << rep_cluster.at(i).dist
      << setw(8)
      << rep_cluster.at(i).index
      << endl;
    for(int i1=0; i1<rep_cluster.at(i).site_num; i1++) {
      fout << setw(16)
        << rep_cluster.at(i).structure.atoms.at(i1).cpos[1]
        << setw(16)
        << rep_cluster.at(i).structure.atoms.at(i1).cpos[2]
        << setw(16)
        << rep_cluster.at(i).structure.atoms.at(i1).cpos[3]
        << setw(16)
        << rep_cluster.at(i).structure.indices.at(i1)
        << endl;
    }

    for(int j=0; j<rep_cluster.at(i).equivalent_num; j++) {

      fout << all_cluster.at(i).at(j).pair_name << endl;
      fout << all_cluster.at(i).at(j).equivalent_num << endl;
      fout << setw(8)
        << all_cluster.at(i).at(j).site_num
        << setw(8)
        << all_cluster.at(i).at(j).NNNum
        << setw(8)
        << all_cluster.at(i).at(j).dist
        << setw(8)
        << all_cluster.at(i).at(j).index
        << endl;
      for(int i1=0; i1<all_cluster.at(i).at(j).site_num; i1++) {
        fout << setw(16)
          << all_cluster.at(i).at(j).structure.atoms.at(i1).cpos(1)
          << setw(16)
          << all_cluster.at(i).at(j).structure.atoms.at(i1).cpos(2)
          << setw(16)
          << all_cluster.at(i).at(j).structure.atoms.at(i1).cpos(3)
          << setw(16)
          << all_cluster.at(i).at(j).structure.indices.at(i1)
          << endl;
      }
    }

  }


  // write NN_distance, NN_shell_num, and atom_list

  fout << NN_distance.size() << endl;
  for(uint i=0; i<NN_distance.size(); i++) {
    fout << setw(8)
      << NN_distance.at(i)
      << setw(8)
      << NN_shell_num.at(i)
      << endl;
  }

  fout << atom_list.size() << endl;
  for(uint i=0; i<atom_list.size(); i++) {
    fout << setw(16)
      << atom_list.at(i).cpos(1)
      << setw(16)
      << atom_list.at(i).cpos(2)
      << setw(16)
      << atom_list.at(i).cpos(3)
      << endl;
  }

  fout.close();

  string commend;

  return true;

}

bool comparison_cluster_list(const vector<cecluster> & cluster_list_1,
    const vector<cecluster> & cluster_list_2) {
  // comparison of two clusters lists grouped by the same pair_name
  // "smaller" cluster list has "small" pair_name
  // sorted by descending order
  // sorted then by atom position

  int size1 = cluster_list_1.size();
  int size2 = cluster_list_2.size();

  if(size1 == 0) {
    if(size2 !=0) {
      return false;
    } else {
      return true;
    }
  } else {
    if(size2 == 0) {
      return true;
    } else {
      return (cluster_list_1.at(0).pair_name < cluster_list_2.at(0).pair_name);
    }
  }
}



// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ceECIcluster.cpp

//////////////////////////////////////////////////////////////
// class ceECIcluster
//////////////////////////////////////////////////////////////

ceECIcluster::ceECIcluster() {
  // default constructor
  allECI_cluster.clear();
  ECI_cluster.clear();
  ECI.clear();
  str_cluster_calculated = false;
  str_cluster = 0; // not assign any memory
}

ceECIcluster::ceECIcluster(ceECIcluster & ECIcluster) {
  // copy constructor
  //ECI_cluster = ECIcluster.ECI_cluster;
  allECI_cluster = ECIcluster.ECI_cluster;
  ECI_cluster = ECIcluster.ECI_cluster;
  ECI = ECIcluster.ECI;
  str_cluster_calculated = ECIcluster.str_cluster_calculated;
  str_cluster = ECIcluster.str_cluster; // not assign any memory
}

ceECIcluster::~ceECIcluster() {
  // default destructor;
  allECI_cluster.clear();
  ECI_cluster.clear();
  ECI.clear();
  str_cluster = 0;
}


ceECIcluster & ceECIcluster::operator=(ceECIcluster & ECIcluster) {
  //ECI_cluster = ECIcluster.ECI_cluster;
  allECI_cluster = ECIcluster.allECI_cluster;
  ECI_cluster = ECIcluster.ECI_cluster;
  ECI = ECIcluster.ECI;
  str_cluster_calculated = ECIcluster.str_cluster_calculated;
  str_cluster = ECIcluster.str_cluster; // not assign any memory
  return *this;
}

void ceECIcluster::SetUp(ceallclusters & cluster1) {
  // set up everything
  SetStrCluster(cluster1);
  GetECICluster();
  //GetECICluster();
}

void ceECIcluster::SetStrCluster(const ceallclusters & ceallcluster_in) {

  if(!str_cluster_calculated) {
    str_cluster = &ceallcluster_in;
    str_cluster_calculated = true;
  }
}


//void ceECIcluster::GetECICluster()
//{
//    // keep only clusters used in fit in all_cluster
//    // In this function, the criterion for choosing clusters are set by hand
//    // using the predefine cutoff
//    // In future version, the clusters will be generated and adjusted
//    // automatically
//    // have to include all subclusters in the CVM in order to the
//    // the probabilities of all clusters in CVM
//
//    int index;
//    vector<int> fit_cluster_list;
//    cecluster cluster1;
//
//    fit_cluster_list = ECI_cluster;
//
//    bool flag;
//
//
//    // output the name of fit cluster
//    cerr << "The ECI & CVM clusters are" << endl;
//    for(uint i=0; i<ECI_cluster.size(); i++) {
//        index = ECI_cluster.at(i);
//
//        cerr << "i " << i
//            << " | "
//            << " index "
//            << index
//            << " | "
//            << " " << str_cluster->rep_cluster.at(index).pair_name << endl;
//    }
//
//    ECI_cluster = fit_cluster_list;
//
//}

void ceECIcluster::GetECICluster() {
  // initially set the ECI_cluster
  // use all cluster in the pool

  int index;

  ECI_cluster.clear();

  for(uint i=0; i<str_cluster->rep_cluster.size(); i++) {
    ECI_cluster.push_back(i);
  }

  // output the name of fit cluster
  cerr << "The ECI fit clusters are" << endl;
  for(uint i=0; i<ECI_cluster.size(); i++) {
    index = ECI_cluster.at(i);

    cerr << "i " << i
      << " | "
      << " index "
      << index
      << " | "
      << " " << str_cluster->rep_cluster.at(index).pair_name << endl;
  }

  allECI_cluster = ECI_cluster;

}


void ceECIcluster::GetECICluster(vector<int> & ECI_cluster_in) {
  ECI_cluster = ECI_cluster_in;

  if(allECI_cluster.size() < ECI_cluster_in.size()) {
    allECI_cluster = ECI_cluster;
  }

  //// output the name of fit cluster
  //cerr << "From assignment, the ECI fit clusters are" << endl;
  //for(uint i=0; i<ECI_cluster.size(); i++) {
  //    int index;
  //    index = ECI_cluster.at(i);

  //    cerr << "i " << i
  //        << " | "
  //        << " index "
  //        << index
  //        << " | "
  //        << " " << str_cluster->rep_cluster.at(index).pair_name << endl;
  //}

}

void ceECIcluster::PrintOutECICluster() {
  // initially set the ECI_cluster
  // use all cluster in the pool

  int index;

  //// output the name of fit cluster
  //cerr << "The ECI fit clusters are" << endl;
  //for(uint i=0; i<allECI_cluster.size(); i++) {
  //    index = allECI_cluster.at(i);

  //    cerr << "i " << i
  //        << " | "
  //        << " index "
  //        << index
  //        << " | "
  //        << " " << str_cluster->rep_cluster.at(index).pair_name << endl;
  //}

  // output the name of fit cluster
  cerr << "The ECI fit clusters are" << endl;
  for(uint i=0; i<ECI_cluster.size(); i++) {
    index = ECI_cluster.at(i);

    cerr << "i " << i
      << " | "
      << " index "
      << index
      << " | "
      << " " << str_cluster->rep_cluster.at(index).pair_name << endl;
  }


}


string ceECIcluster::GetClusterNameByIndex(int index) {
  return str_cluster->rep_cluster.at(index).pair_name;
}


void ceECIcluster::PrintECI(ostream & os) {
  ios_base::fmtflags old_stat = os.setf(ios_base::fixed, ios_base::floatfield);
  os.precision(6);

  os << "ECIs " << ECI_cluster.size() << endl;

  os << setw(12)
    << "0 0 0 0"
    << setw(12)
    << ECIValue().at(0)
    << endl;
  for(uint i=0; i<ECI_cluster.size(); i++) {
    int index = ECI_cluster.at(i);
    os << setw(12)
      << str_cluster->rep_cluster.at(index).pair_name
      << setw(12)
      << ECI.at(i+1)
      << endl;
  }

  os.setf(old_stat, ios_base::floatfield);
}

void ceECIcluster::WriteFile(ostream & os) {
  // only write ECI, ECICluster
  // as CVM method will be dropped off

  ios_base::fmtflags old_stat = os.setf(ios_base::fixed, ios_base::floatfield);
  os.precision(12);
  int width=12+6;

  os << ECI.size() << endl;

  // ECI values
  for(uint i=0; i<ECI.size(); i++) {
    os << setw(width)
      << ECI.at(i)
      << endl;
  }

  // ECI clusters
  os << ECI_cluster.size() << endl;
  for(uint i=0; i<ECI_cluster.size(); i++) {
    int index = ECI_cluster.at(i);

    int atom_num = str_cluster->rep_cluster.at(index).structure.atoms.size();

    os << atom_num << endl;

    for(int i1=0; i1 < atom_num; i1++) {
      os << setw(width)
        << str_cluster->rep_cluster.at(index).structure.atoms.at(i1).cpos[1]
        << setw(width)
        << str_cluster->rep_cluster.at(index).structure.atoms.at(i1).cpos[2]
        << setw(width)
        << str_cluster->rep_cluster.at(index).structure.atoms.at(i1).cpos[3]
        << endl;
    }
  }

  os.setf(old_stat, ios_base::floatfield);

}

void ceECIcluster::ReadIn(istream & os) {

  int ECI_size;
  os >> ECI_size;

  if(str_cluster == 0) {
    throw aurostd::xerror(__AFLOW_FILE__,"ceECIcluster::ReadIn():","clusters must be supplied",_FILE_CORRUPT_); //CO20200624
  }

  ECI.clear();
  for(int i=0; i<ECI_size; i++) {
    double eci_tmp;
    os >> eci_tmp;
    ECI.push_back(eci_tmp);
  }

  int ECI_cluster_num;
  os >> ECI_cluster_num;

  int  atom_num;
  ECI_cluster.clear();
  for(int i=0; i<ECI_cluster_num; i++) {
    os >> atom_num;

    cecluster cluster_tmp;
    _atom atom_tmp;
    for(int i1=0; i1<atom_num; i1++) {
      os >> atom_tmp.cpos[1]
        >> atom_tmp.cpos[2]
        >> atom_tmp.cpos[3];
      cluster_tmp.structure.atoms.push_back(atom_tmp);
    }

    int index;

    for(uint j1=0; j1<str_cluster->rep_cluster.size(); j1++) {
      if(is_equal(cluster_tmp, str_cluster->rep_cluster.at(j1))) {
        index = j1;
        break;
      }
    }

    ECI_cluster.push_back(index);
  }

  if(allECI_cluster.size() <= ECI_cluster.size()) {
    allECI_cluster = ECI_cluster;
  }
}

void ceECIcluster::DeleteZeroECI() {
  // delete clusters with zero ECI's
  // but keep (0 0 0 0)

  vector<double> ECI_tmp;
  vector<int> ECIcluster_tmp;

  ECI_tmp.push_back(ECI.at(0));
  for(uint i=1; i<ECI.size(); i++) {
    if(abs(ECI.at(i)) > _ECI_ZERO) {
      ECI_tmp.push_back(ECI.at(i));
      ECIcluster_tmp.push_back(ECI_cluster.at(i-1));
    }
  }

  ECI_cluster = ECIcluster_tmp;
  ECI = ECI_tmp;
}



// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// cestructure.cpp

//////////////////////////////////////////////////////////////
// cestructure class
//////////////////////////////////////////////////////////////

cestructure::cestructure() : structure() {
  // default constructor
  name = "";
  energy_in = 0.0;
  stoich_b = 0.0;
  energy = 0.0;
  ECI_correlation.clear();
  ECI_cluster.clear();
  ECI_equivalent_num.clear();
  allECI_correlation.clear();
  allECI_cluster.clear();
  allECI_equivalent_num.clear();

  str_cluster = 0; // void the pointer
}

cestructure::cestructure(string & str_name, double stoich_b_in,
    double fit_quantity) : structure() {
  name = str_name;
  stoich_b = stoich_b_in;
  energy_in = fit_quantity;

  energy = 0.0;
  ECI_correlation.clear();
  ECI_cluster.clear();
  ECI_equivalent_num.clear();
  allECI_correlation.clear();
  allECI_cluster.clear();
  allECI_equivalent_num.clear();

  // set the volumn of the structure
  //structure=aflowlib::PrototypeLibraries(cerr, name,"",LIBRARY_MODE_HTQC); // no parameters;
  //int volumn_in=1.0;
  //structure.SetVolume(volumn_in);
  //structure.FixLattices();

  // set the volumn of the structure
  xstructure xstr;
  xstr=aflowlib::PrototypeLibraries(cerr, name,"",LIBRARY_MODE_HTQC); // no parameters;
  structure.SetStructure(xstr);

  str_cluster = 0; // void the pointer
}

cestructure::cestructure(xstructure & xstr) : structure(xstr) {
  // default constructor
  name = structure.prototype;
  energy_in = 0.0;
  stoich_b = 0.0;
  energy = 0.0;
  ECI_correlation.clear();
  ECI_cluster.clear();

  str_cluster = 0; // void the pointer


  //structure.FixLattices();
}

cestructure::~cestructure() {
  // default destructor
  ECI_correlation.clear();
  ECI_cluster.clear();

  str_cluster = 0; // void the pointer
}

ostream & operator<<(ostream & os, const cestructure & cestr) {
  // entropy and free energy are not calculated by CVM
  // but by MC
  // therefore, they are not output here

  os.setf(ios_base::fixed, ios_base::floatfield);
  os.precision(8);

  os << setw(6)
    << cestr.name
    << setw(12)
    << cestr.stoich_b
    << setw(12)
    << cestr.energy_in;

  return os;
}

void cestructure::SetUp(const ceallclusters & cluster1, ceECIcluster & cluster2) {
  // set up everthing
  SetStrCluster(cluster1);
  SetAllECICluster(cluster2);
  GetAllECICorrelation();
  PrintOutCorrelation();
}

void cestructure::SetUp(const ceallclusters & cluster1, ceECIcluster & cluster2,
    istream & corfilein) {
  // set up everthing
  SetStrCluster(cluster1);
  SetAllECICluster(cluster2);
  GetAllECICorrelation(corfilein);
  PrintOutCorrelation();
}

void cestructure::SetStrCluster(const ceallclusters & ceallcluster_in) {

  //if(!str_cluster_calculated) {
  str_cluster = &ceallcluster_in;
  str_cluster_calculated = true;
  //}
}

void cestructure::SetAllECICluster(ceECIcluster & subcluster) {
  // set ECI_cluster and CVM_cluster

  allECI_cluster = subcluster.AllECICluster();
  ECI_cluster = allECI_cluster;
}

void cestructure::SetECICluster(ceECIcluster & subcluster) {
  // set ECI_cluster and CVM_cluster

  ECI_cluster = subcluster.ECICluster();
}

vector<double> cestructure::GetCorrelation(vector<int> & index_list) {
  // get the correlation function
  // for those clusters with index in index_list

  vector<double> correlation_list;

  string s;

  double cors[index_list.size()]; // correlations of clusters in index_list

  // get the configuration and the equivalent configurations by
  // shifting of origin in primitive unit cell

  // first get the atom positions, lattice vectors and so on of each fit structure
  //xstructure str_tmp;
  cexstructure str_tmp;

  // scale of aflow lib structure is set somewhere to 3.12
  // need to set it back to 1.0
  //structure.scale = 1.0;

  str_tmp = structure;
  //str_tmp.FixLattices();

  deque<_atom> atom_list = str_cluster->AtomList();

  vector<int> configuration;

  for(uint i=0; i<index_list.size(); i++) {
    cors[i] = 0.0;
  }

  // first group atoms in atom_list wrt base atoms
  vector< vector<int> > equiv_atoms;
  deque<_atom> atoms_org;

  // To be consistent, set all fpos to non-negative value
  for(uint i1=0; i1<structure.atoms.size(); i1++) {
    for(uint l1=1; l1<4; l1++) {
      if(str_tmp.atoms.at(i1).fpos[l1] < 0) {
        str_tmp.atoms.at(i1).fpos[l1] +=
          abs(int(str_tmp.atoms.at(i1).fpos[l1])) + 1.0;
      }
    }

    str_tmp.atoms.at(i1).cpos = F2C(str_tmp.lattice, str_tmp.atoms.at(i1).fpos);
  }
  atoms_org = structure.atoms;

  for(uint i=0; i< str_tmp.atoms.size(); i++) {
    vector<int> vec_tmp;
    equiv_atoms.push_back(vec_tmp);
    equiv_atoms.at(i).clear();
  }

  double g11, g12, g13;
  for(uint i1=0; i1 < atom_list.size(); i1++) {

    _atom atom_tmp;

    //str_tmp.FixLattices();

    for(uint l2=0;l2 < str_tmp.atoms.size(); l2++) {
      atom_tmp.cpos = atom_list.at(i1).cpos -str_tmp.atoms.at(l2).cpos;

      g11 = abs(str_tmp.klattice[1][1]*atom_tmp.cpos[1]
          +str_tmp.klattice[1][2]*atom_tmp.cpos[2]
          +str_tmp.klattice[1][3]*atom_tmp.cpos[3])/(2.0*_pi)+_EQUAL_DOUBLE*0.1;
      g12 = abs(str_tmp.klattice[2][1]*atom_tmp.cpos[1]
          +str_tmp.klattice[2][2]*atom_tmp.cpos[2]
          +str_tmp.klattice[2][3]*atom_tmp.cpos[3])/(2.0*_pi)+_EQUAL_DOUBLE*0.1;
      g13 = abs(str_tmp.klattice[3][1]*atom_tmp.cpos[1]
          +str_tmp.klattice[3][2]*atom_tmp.cpos[2]
          +str_tmp.klattice[3][3]*atom_tmp.cpos[3])/(2.0*_pi)+_EQUAL_DOUBLE*0.1;


      if(( abs(g11 - int(g11)) < _EQUAL_DOUBLE )
          &&(abs(g12 - int(g12)) < _EQUAL_DOUBLE )
          &&(abs(g13 - int(g13)) < _EQUAL_DOUBLE )
        ) {
        equiv_atoms.at(l2).push_back(i1);
        break;
      }
    }

  }

  //cerr << "atom config\n";
  //for(uint l1 = 0; l1<equiv_atoms.size(); l1++) {
  //    cerr << structure.atoms.at(l1).name << " : " ;
  //    for(uint l2 = 0; l2<equiv_atoms.at(l1).size(); l2++) {
  //        cerr << equiv_atoms.at(l1).at(l2) << " ";
  //    }
  //    cerr << endl;
  //}

  for(uint i=0; i<structure.atoms.size(); i++) {
    // number of equivalent configurations is the same
    // as the number of atom in the primitive unit cell

    str_tmp.atoms = atoms_org;
    deque<_atom> atoms_tmp;
    atoms_tmp = atoms_org;
    for(uint i1=0; i1<structure.atoms.size(); i1++) {
      for(uint l1=1; l1<4; l1++) {
        str_tmp.atoms.at(i1).cpos[l1] = structure.atoms.at(i1).cpos[l1]
          - structure.atoms.at(i).cpos[l1];
      }

      str_tmp.atoms.at(i1).fpos = C2F(str_tmp.lattice, str_tmp.atoms.at(i1).cpos);

    }

    //cerr << "base atoms\n";

    // shift the atoms back to the primitive unit cell
    for(uint i1=0; i1<structure.atoms.size(); i1++) {
      for(uint l1=1; l1<4; l1++) {
        // fpos of atom in the primitive unit cell
        // is always positive
        // only shift those atom with negative fpos elements
        if(str_tmp.atoms.at(i1).fpos[l1] < -1.0*_EQUAL_DOUBLE) {
          str_tmp.atoms.at(i1).fpos[l1] +=
            abs(int(str_tmp.atoms.at(i1).fpos[l1])) + 1.0;
        }
      }

      str_tmp.atoms.at(i1).cpos = F2C(str_tmp.lattice, str_tmp.atoms.at(i1).fpos);

      //cerr << str_tmp.atoms.at(i1) << " " << str_tmp.atoms.at(i1).name << endl << endl;
    }

    // restore the order of atoms in unit cell with new configuration
    for(uint i1=0; i1 < atoms_tmp.size(); i1++) {
      for(uint i2=0; i2 < str_tmp.atoms.size(); i2++) {
        if(is_equal(atoms_tmp.at(i1), str_tmp.atoms.at(i2))) {
          atoms_tmp.at(i1).name = str_tmp.atoms.at(i2).name;
          break;
        }
      }
    }

    str_tmp.atoms = atoms_tmp;

    //for(uint i1=0; i1 < atoms_org.size(); i1++) {
    //    str_tmp.atoms.at(i1).cpos = F2C(str_tmp.lattice, str_tmp.atoms.at(i1).fpos);
    //    cerr << str_tmp.atoms.at(i1) << " " << str_tmp.atoms.at(i1).name << endl;
    //}

    for(uint i1=0; i1 < equiv_atoms.size(); i1++) {
      for(uint i2=0; i2 < equiv_atoms.at(i1).size(); i2++) {
        atom_list.at(equiv_atoms.at(i1).at(i2)).name = str_tmp.atoms.at(i1).name;
      }
    }

    //cerr << "\n configuration\n";
    //for(uint i1=0; i1 < atom_list.size(); i1++) {
    //    cerr << atom_list.at(i1).name << " ";
    //}
    //cerr << endl;


    configuration.clear();
    for(uint i1=0; i1<atom_list.size(); i1++) {
      char ch = atom_list.at(i1).name.at(0);
      switch(ch) {
        case 'A':
        case 'a':
          configuration.push_back(_A_ATOM);
          break;
        case 'B':
        case 'b':
          configuration.push_back(_B_ATOM);
          break;
      }
    }

    //// output configurations
    //for(uint i1=0; i1<atom_list.size(); i1++) {
    //    cerr << setw(4)
    //        << atom_list.at(i1).name << " ";
    //}
    //cerr << endl;

    //for(uint i1=0; i1<configuration.size(); i1++) {
    //    cerr << setw(4)
    //        << configuration.at(i1) << " ";
    //}
    //cerr << endl;


    // calcuate the correlations in the current configuration
    for(uint j=0; j<index_list.size(); j++) {

      for(uint j1=0; j1<str_cluster->all_cluster.at(j).size(); j1++) {

        cecluster cluster1 = str_cluster->all_cluster.at(j).at(j1);
        int cor_tmp = 1;
        for(uint j2=0; j2<cluster1.structure.indices.size(); j2++) {
          int ind = cluster1.structure.indices.at(j2);
          cor_tmp *= configuration.at(ind);
        }
        cors[j] += double(cor_tmp);

      }

    }


  }

  for(uint i=0; i<index_list.size(); i++) {
    cors[i] /= double(str_cluster->rep_cluster.at(i).equivalent_num);
    cors[i] /= double(structure.atoms.size());
    correlation_list.push_back(cors[i]);
  }

  // add the uniform term
  correlation_list.insert(correlation_list.begin(), 1.0);

  return correlation_list;

}

void cestructure::GetAllECICorrelation() {
  // get correlations for ECI clusters

  allECI_correlation = GetCorrelation(allECI_cluster);

  // set the equivalent_num

  allECI_equivalent_num.clear();

  allECI_equivalent_num.push_back(1); // "0 0 0 0" cluster
  for(uint i=0; i<ECI_cluster.size(); i++) {
    int index = ECI_cluster.at(i);
    allECI_equivalent_num.push_back(str_cluster->rep_cluster.at(index).equivalent_num);
  }

  ECI_cluster = allECI_cluster;
  ECI_correlation = allECI_correlation;
  ECI_equivalent_num = allECI_equivalent_num;

}

void cestructure::GetECICorrelation() {
  // get correlations for ECI clusters

  //ECI_correlation = GetCorrelation(ECI_cluster);

  if(allECI_correlation.size() == 0) {
    GetAllECICorrelation();
  }

  ECI_correlation.clear();

  ECI_correlation.push_back(1.0); // "0 0 0 0" correlation
  for(uint i=0; i<ECI_cluster.size(); i++) {
    for(uint j=0; j<allECI_cluster.size(); j++) {
      if(ECI_cluster.at(i) == allECI_cluster.at(j)) {

        // ******************************************************
        // first item in ECI_corrleation is the correlation of "0 0 0 0"
        // this structure is not stored in EC_cluster and ECI_cluster!
        // that's why j+1 is used here
        // ******************************************************

        ECI_correlation.push_back(allECI_correlation.at(j+1));
        break;
      }
    }
  }

  ECI_equivalent_num.clear();

  // set the equivalent_num
  ECI_equivalent_num.push_back(1); // "0 0 0 0" cluster
  for(uint i=0; i<ECI_cluster.size(); i++) {
    int index = ECI_cluster.at(i);
    ECI_equivalent_num.push_back(str_cluster->rep_cluster.at(index).equivalent_num);
  }


}

void cestructure::GetAllECICorrelation(istream & filein) {
  // get correlations for ECI clusters

  ReadIn(filein);

  //// set the equivalent_num

  allECI_equivalent_num.clear();

  allECI_equivalent_num.push_back(1); // "0 0 0 0" cluster
  for(uint i=0; i<allECI_cluster.size(); i++) {
    int index = allECI_cluster.at(i);
    allECI_equivalent_num.push_back(str_cluster->rep_cluster.at(index).equivalent_num);
  }

  ECI_cluster = allECI_cluster;
  ECI_correlation = allECI_correlation;
  ECI_equivalent_num = allECI_equivalent_num;


}

int cestructure::GetAllECICorrelation(istream & filein, int pos) {
  // get correlations for ECI clusters

  int pos_new;
  pos_new = ReadIn(filein, pos);

  if(pos_new != _NOCORRELATON_FOUND) {

    allECI_equivalent_num.clear();

    allECI_equivalent_num.push_back(1); // "0 0 0 0" cluster
    for(uint i=0; i<allECI_cluster.size(); i++) {
      int index = allECI_cluster.at(i);
      allECI_equivalent_num.push_back(str_cluster->rep_cluster.at(index).equivalent_num);
    }

  }

  return pos_new;

}


void cestructure::PrintOutCorrelation() {
  cerr << endl;
  cerr << "ECI correlation list size " << ECI_correlation.size() << endl;

  cerr << "structure name " << name << endl;
  for(uint i=0; i<ECI_correlation.size(); i++) {
    if(i==0) {
      cerr << setw(12)
        <<"0 0 0 0";
    } else {
      int index = ECI_cluster.at(i-1);
      cerr << setw(12)
        << str_cluster->rep_cluster.at(index).pair_name;
    }
    cerr.setf(ios_base::fixed, ios_base::floatfield);
    cerr.precision(6);
    cerr << setw(12)
      << "eq cluster "
      << setw(4)
      << ECI_equivalent_num.at(i)
      << setw(6)
      <<" cor "
      << setw(10)
      << ECI_correlation.at(i) << endl;
  }


  cerr << endl;

}

void cestructure::WriteFile(ostream & os) {

  os << setw(10)
    << name
    << setw(10)
    //<< ECI_correlation.size()
    << ECI_correlation.size()
    << endl;


  for(uint i=0; i<ECI_correlation.size(); i++) {
    if(i==0) {
      os << setw(12)
        <<"0 0 0 0" ;
    } else {
      int index = ECI_cluster.at(i-1);
      os << setw(12)
        << str_cluster->rep_cluster.at(index).pair_name ;
    }
    os.setf(ios_base::fixed, ios_base::floatfield);
    os.precision(8);
    os << setw(12)
      << ECI_correlation.at(i) << endl;
  }

  os.precision(6);

}


bool cestructure::ReadIn(istream & fin) {
  // read data (correlations of ECI_cluster from a file)


  string str_name;
  int cor_num;
  bool flag;

  fin.seekg(0, ios_base::beg);

  while(fin >> str_name >> cor_num) {


    if((str_name != name) || (cor_num != int(allECI_cluster.size()) + 1)) {
      //getline(fin, str_name);
      flag = false;
      for(int i = 0; i< cor_num; i++) {

        string pair_name;
        int site_num, NNNum, dist, index;
        double cor;

        fin >> site_num
          >> NNNum
          >> dist
          >> index
          >> cor;

      }
      continue;
    }

    allECI_correlation.clear();
    for(uint i=0; i<allECI_cluster.size()+1; i++) {
      allECI_correlation.push_back(0.0);
    }

    for(int i = 0; i< cor_num; i++) {
      string pair_name;
      int site_num, NNNum, dist, index;
      double cor;

      fin >> site_num
        >> NNNum
        >> dist
        >> index
        >> cor;

      allECI_correlation.at(index) = cor;
    }

    flag = true;
    break;
  }


  return flag;
}

int cestructure::ReadIn(istream & fin, int pos) {
  // read data (correlations of ECI_cluster from a file)

  string str_name;
  int cor_num;

  int length_new = _NOCORRELATON_FOUND;

  // if pos < 0, read the file in sequence
  if(fin.eof()) {
    // clear the eof bit in order to read the file
    // after failing to find the correlations of an SL
    fin.clear();
  }
  if(pos >= 0) {
    fin.seekg(pos, ios_base::beg);
  }

  while (!fin.eof()) {

    fin >> str_name >> cor_num;

    if((str_name != name) || (cor_num != int(allECI_cluster.size()) + 1)) {
      getline(fin, str_name);
      for(int i = 0; i< cor_num; i++) {

        getline(fin, str_name);
      }
      continue;
    }


    allECI_correlation.clear();
    for(uint i=0; i<allECI_cluster.size()+1; i++) {
      allECI_correlation.push_back(0.0);
    }

    for(int i = 0; i< cor_num; i++) {
      string pair_name;
      int site_num, NNNum, dist, index;
      double cor;

      fin >> site_num
        >> NNNum
        >> dist
        >> index
        >> cor;

      allECI_correlation.at(index) = cor;
    }

    length_new = fin.tellg();

    break;
  }

  return length_new;
}

int cestructure::NameToValue(string name) {
  // assign 'spin' value to different atoms
  // only valid for two atoms
  int val;
  if(tolower(name[0]) == 'a')val = _A_ATOM;
  if(tolower(name[0]) == 'b')val = _B_ATOM;
  return val;
}

void cestructure::GetEnergy() {

  energy = 0.0;

  //cerr << "ECI correlation size " << ECI_correlation.size() << endl;
  //cerr << "ECI size " << ECI.size() << endl;
  //cerr << "ECI equivalent size " << ECI_equivalent_num.size() << endl;

  for(uint j=0; j<ECI_correlation.size(); j++) {
    //cerr << "j = " << j << " " << ECI.at(j)
    //    << " " << ECI_correlation.at(j)
    //    << " " << ECI_equivalent_num.at(j) << endl;
    energy += ECI.at(j)*ECI_correlation.at(j)
      *ECI_equivalent_num.at(j);
  }

}


void cestructure::PrintOutComparison(ostream & os) {
  ios_base::fmtflags stat_old;
  stat_old=os.setf(ios_base::fixed, ios_base::floatfield);
  os.precision(6);

  double de;
  if(energy_in != 0) {
    de = (energy-energy_in)/energy_in;
  } else {
    de = (energy-energy_in)/energy;
  }


  os << setw(6)
    << name
    << setw(12)
    << stoich_b
    << setw(12)
    << energy
    << setw(12)
    << energy_in
    << setw(12)
    << de
    << endl;

  // restore to the old setting
  os.setf(stat_old, ios_base::floatfield);
}



// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ceralloy.cpp

//////////////////////////////////////////////////////////////
// class ceralloy
// derived from cestructure
// for random alloy
//////////////////////////////////////////////////////////////

ceralloy::ceralloy() : cestructure() {
}

//ceralloy::ceralloy(string & str_name, double stoich_b,
//        double energy, double formation_energy,
//        double entropy, double volumn,
//        double temperature) : cestructure(str_name, stoich_b,
//            energy, formation_energy,
//            entropy, volumn,
//            temperature)
//{
//}

//ceralloy::ceralloy(string & str_name, double stoich_b,
//        double temperature) : cestructure(str_name, stoich_b,
//            0.0, 0.0, 0.0, 0.0, 1.0, temperature)
ceralloy::ceralloy(string & str_name, double stoich_b) : cestructure(str_name, stoich_b, 0.0) {
}

ceralloy::ceralloy(xstructure & xstr) : cestructure(xstr) {
}

ceralloy::~ceralloy() {
}

ostream & operator<<(ostream & os, const ceralloy & cestr) {
  // entropy and free energy are not calculated by CVM
  // but by MC
  // therefore, they are not output here

  os.setf(ios_base::fixed, ios_base::floatfield);
  os.precision(8);

  os << setw(6)
    << cestr.name
    << setw(12)
    << cestr.stoich_b
    << setw(12)
    << cestr.energy_in
    << setw(12)
    << cestr.energy;

  return os;
}

void ceralloy::SetUp(ceallclusters & allcluster, ceECIcluster & ECIcluster) {

  name = name+"-r";
  SetStrCluster(allcluster);
  //SetECICluster(ECIcluster);
  //GetECICorrelation();

  SetAllECICluster(ECIcluster);
  SetECICluster(ECIcluster);
  GetAllECICorrelation();

  PrintOutCorrelation();

  //SetaNum(ECIcluster.aNum());
  SetECI(ECIcluster.ECIValue());
}

void ceralloy::GetAllECICorrelation() {
  // get correlations for ECI & CVM clusters

  allECI_correlation.clear();
  allECI_correlation.push_back(1.0); // "0 0 0 0" correlation
  for(uint i=0; i<allECI_cluster.size(); i++) {
    int index = allECI_cluster.at(i);
    int site_num = str_cluster->rep_cluster.at(index).site_num;
    double tmp = (1.0 - stoich_b)*_A_ATOM + stoich_b*_B_ATOM;
    double cor = pow(tmp, site_num);
    allECI_correlation.push_back(cor);
  }

  // get correlations for ECI clusters
  ECI_correlation.clear();

  ECI_correlation.push_back(1.0); // "0 0 0 0" correlation
  for(uint i=0; i<ECI_cluster.size(); i++) {
    for(uint j=0; j<allECI_cluster.size(); j++) {
      if(ECI_cluster.at(i) == allECI_cluster.at(j)) {

        // ******************************************************
        // first item in ECI_corrleation is the correlation of "0 0 0 0"
        // this structure is not stored in EC_cluster and CVM_cluster!
        // that's why j+1 is used here
        // ******************************************************

        ECI_correlation.push_back(allECI_correlation.at(j+1));
        break;
      }
    }
  }

  // set the equivalent_num

  ECI_equivalent_num.clear();

  ECI_equivalent_num.push_back(1); // "0 0 0 0" cluster
  for(uint i=0; i<ECI_cluster.size(); i++) {
    int index = ECI_cluster.at(i);
    ECI_equivalent_num.push_back(str_cluster->rep_cluster.at(index).equivalent_num);
  }

}


// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// cesubcluster.cpp


//////////////////////////////////////////////////////////////
// class cesubcluster
//////////////////////////////////////////////////////////////


cesubcluster::cesubcluster() {
  name.clear();
  base_cluster.clear();
  overlap_cluster.clear();
  sub_cluster.clear();
  str_cluster_calculated = false;
  str_cluster = 0; // not assign any memory
  overlap_cluster_calculated = false;
}

cesubcluster::~cesubcluster() {
  name.clear();
  base_cluster.clear();
  overlap_cluster.clear();
  sub_cluster.clear();
  str_cluster = 0; // not assign any memory
}

void cesubcluster::SetStrCluster(const ceallclusters & ceallcluster_in) {

  if(!str_cluster_calculated) {
    str_cluster = &ceallcluster_in;
    str_cluster_calculated = true;
  }
}

void cesubcluster::SetStrCluster(const ceallclusters * ceallcluster_in) {

  if(!str_cluster_calculated) {
    str_cluster = ceallcluster_in;
    str_cluster_calculated = true;
  }
}

void cesubcluster::SetUp(vector<int> base_cluster, ceallclusters & cluster1) {
  SetBaseCluster(base_cluster.at(0));
  SetStrCluster(cluster1);

  for(uint i=1; i<base_cluster.size(); i++) {
    AddBaseCluster(base_cluster.at(i));
  }
  GetOverlapCluster();


}

void cesubcluster::SetBaseCluster(int base_cluster_nr) {

  if(base_cluster.size() !=0) {
    throw aurostd::xerror(__AFLOW_FILE__,"cesubcluster::SetBaseCluster():","already has a base",_INPUT_ILLEGAL_); //CO20200624
  } else {
    base_cluster.push_back(base_cluster_nr);
  }
}

void cesubcluster::SetBaseCluster(vector<int> base_cluster_nr) {
  //if(str_cluster_calculated) {
  base_cluster = base_cluster_nr;
  //} else {
  //    throw aurostd::xerror(__AFLOW_FILE__,"cesubcluster::SetBaseCluster():","must get all clusters of a structure first",_INPUT_ILLEGAL_); //CO20200624
  //}
}

void cesubcluster::AddBaseCluster(int base_cluster_nr) {

  base_cluster.push_back(base_cluster_nr);
}

void cesubcluster::AddBaseCluster(vector<int> base_cluster_nr) {

  for(uint i=0; i<base_cluster_nr.size(); i++) {
    AddBaseCluster(base_cluster_nr.at(i));
  }
}

vector<Baker_str> cesubcluster::GetSubCluster(int base_cluster) {
  // get all subclusters of a base cluster

  // get base cluster structure
  cecluster cluster1 = str_cluster->rep_cluster.at(base_cluster-1);
  vector<Baker_str> sub_cluster;
  Baker_str subcluster_tmp;

  //cesubcluster subcluster_tmp;

  int total_atom_num = cluster1.structure.atoms.size();

  vector< vector<int> > str;
  int atom_index, atom_num;
  string cluster_name;

  cecluster cluster2, cluster_rep;

  bool flag=false;
  bool flag_found=false; if(flag_found) {};
  _atom atom_tmp, atom_ref;

  int dist = cluster1.dist;
  int site_num = cluster1.site_num; if(site_num) {};
  int NNNum = cluster1.NNNum; if(NNNum) {};
  int dist_tmp;
  int index_tmp = cluster1.index; if(index_tmp) {};

  // site_num = cluster1.site_num;
  // NNNum = cluster1.NNNum;
  // dist = cluster1.dist;
  // index_tmp = cluster1.index;

  for(int i=total_atom_num-1; i>= 0; i--) {

    int str_size = CombinationNr(i+1, total_atom_num);

    vector<int> atom_config;
    for(int j=0; j<str_size; j++) {

      if(j==0) {
        atom_config = AllCombination41(i+1, total_atom_num, j+1);
      } else {
        atom_config = AllCombination42(i+1, total_atom_num, atom_config);
      }

      atom_num = atom_config.size();
      cluster2.structure.atoms.clear();


      for(int k=0; k<atom_num; k++) {
        atom_index = atom_config.at(k);
        atom_tmp = cluster1.structure.atoms.at(atom_index);
        if(k == 0) {
          atom_ref = atom_tmp;
        }
        for(int k1=1; k1<=3; k1++) {
          // move the cluster with one atom at origin due to the structure of clusters stored in str_cluster.rep_cluster and str_cluster.all_cluster
          atom_tmp.cpos[k1] -= atom_ref.cpos[k1];
        }
        cluster2.structure.atoms.push_back(atom_tmp);
      }

      flag_found = false;
      bool flag_1, flag_2;
      for(uint k=0; k<str_cluster->all_cluster.size(); k++) {

        for(uint k1=0; k1<str_cluster->all_cluster.at(k).size(); k1++) {

          // [OBSOLETE] set but not used 	  int index_tmp;
          site_num = str_cluster->all_cluster.at(k).at(k1).site_num;
          NNNum = str_cluster->all_cluster.at(k).at(k1).NNNum;
          dist_tmp = str_cluster->all_cluster.at(k).at(k1).dist;
          index_tmp = str_cluster->all_cluster.at(k).at(k1).index;

          // clusters with same atom number keep only those smaller ones, i.e. short distance
          // cluster with smaller atom number keep all

          flag_1 = (str_cluster->all_cluster.at(k).at(k1).site_num == total_atom_num)
            &&(str_cluster->all_cluster.at(k).at(k1).NNNum == cluster1.NNNum)
            && (dist_tmp == dist);
          flag_2 = str_cluster->all_cluster.at(k).at(k1).site_num < total_atom_num;


          if(flag_1 || flag_2) {

            cluster_name = str_cluster->all_cluster.at(k).at(k1).pair_name;

            // check if the representive cluster is stored or not
            flag = false;

            for(uint k2=0; k2<sub_cluster.size(); k2++) {
              if(str_cluster->rep_cluster.at(sub_cluster.at(k2).index-1).pair_name
                  == cluster_name) {
                // if found add the number by 1
                if(is_equal(str_cluster->all_cluster.at(k).at(k1), cluster2)) {
                  sub_cluster.at(k2).n_num++;
                }
                flag = true;
                break; //stop k2 loop
              }
            } // k2

            if(flag == false) {
              // if not found add the subcluster

              int n_num;

              if(is_equal (str_cluster->all_cluster.at(k).at(k1), cluster2)) {
                n_num= 1;
              } else {
                n_num = 0;
              }

              int m_tmp = str_cluster->rep_cluster.at(k).equivalent_num
                /str_cluster->rep_cluster.at(k).structure.atoms.size();

              subcluster_tmp.index = k+1;
              subcluster_tmp.m_num = m_tmp;
              subcluster_tmp.n_num = n_num;
              sub_cluster.push_back(subcluster_tmp);

            }

          } // k1
        }

      } //k

    } //j

  } //i


  sort(sub_cluster.begin(), sub_cluster.end(), &comparison_subcluster);


  cerr << "basic cluster" << base_cluster << endl;
  cerr << "sub_cluster size " << sub_cluster.size() << endl;
  for(uint i=0; i<sub_cluster.size(); i++) {
    cerr << sub_cluster.at(i).index
      << " "
      //<< str_cluster->rep_cluster.at(sub_cluster.at(i)).pair_name
      //<< " "
      << sub_cluster.at(i).m_num
      << " "
      << sub_cluster.at(i).n_num
      << " "
      << endl;
  }

  return sub_cluster;

}


void cesubcluster::GetAllSubCluster() {
  for(uint i=0; i<base_cluster.size(); i++) {
    sub_cluster.push_back(GetSubCluster(base_cluster.at(i)));
  }
}


vector<Baker_str> cesubcluster::GetOverlapCluster(int index_in) {
  // get overlapped subclusters of equivalent base structure

  _atom atom1, atom2;
  deque<_atom> atom_list, atom_list_old;
  cecluster cluster1, cluster2, cluster_overlap;
  string subcluster_name;

  vector<Baker_str> subcluster_list_out; // overlapped subclusters

  string pair_name;

  //int site_num, NNNum, dist;
  int index;

  // get the overlapped subclusters

  cerr << "index in " << index_in << endl;
  index = base_cluster.at(index_in-1);
  cerr << "index  " << index << endl;

  cluster1 = str_cluster->all_cluster.at(index-1).at(0); // base_cluster


  for(uint i=0; i< str_cluster->all_cluster.at(index-1).size(); i++) {
    cluster2 = str_cluster->all_cluster.at(index-1).at(i);

    atom_list.clear();
    for(uint j=0; j<cluster1.structure.atoms.size(); j++) {
      atom1 = cluster1.structure.atoms.at(j);
      for(uint k=0; k<cluster2.structure.atoms.size(); k++) {
        atom2 = cluster2.structure.atoms.at(k);
        if(is_equal(atom1, atom2)) {
          atom_list.push_back(atom1);
          break;
        }
      }
    }

    cluster_overlap.structure.atoms = atom_list;

    int index_tmp= str_cluster->GetIndexbyCluster(cluster_overlap);

    bool flag=false;
    for(uint j=0; j<subcluster_list_out.size(); j++) {
      if(index_tmp == subcluster_list_out.at(j).index) {
        // already stored
        flag = true;
        break;
      }
    }
    if(!flag) { // not stored
      for(uint j=0; j<sub_cluster.at(index_in-1).size(); j++) {
        if(index_tmp == sub_cluster.at(index_in-1).at(j).index) {
          subcluster_list_out.push_back(sub_cluster.at(index_in-1).at(j));
          break;
        }
      }
    }
  }

  //overlap_cluster = subcluster_list_out;

  sort(subcluster_list_out.begin(), subcluster_list_out.end(), &comparison_subcluster);

  cerr << "overlap_cluster size " << subcluster_list_out.size() << endl;
  for(uint i=0; i<subcluster_list_out.size(); i++) {
    cerr << subcluster_list_out.at(i).index
      << " "
      << subcluster_list_out.at(i).m_num
      << " "
      << subcluster_list_out.at(i).n_num
      << endl;
  }



  return subcluster_list_out;

}


vector<Baker_str> cesubcluster::GetOverlapCluster(
    const vector<Baker_str> & subcluster_list1,
    const vector<Baker_str> & subcluster_list2) {
  // get the overlap clusters of two base clusters
  // only those clusters contribute the to configurational entropy




  vector<Baker_str> subcluster_list_new1;
  vector<Baker_str> subcluster_list_new2;
  Baker_str subcluster1;
  Baker_str subcluster2;

  int size1=subcluster_list1.size();
  int size2=subcluster_list2.size();

  for(int i=0; i<size1; i++) {
    subcluster1 = subcluster_list1.at(i);
    if(subcluster1.n_num != 0) {
      subcluster_list_new1.push_back(subcluster1);
    }
  }

  for(int i=0; i<size2; i++) {
    subcluster2 = subcluster_list2.at(i);
    for(uint j=0; j<subcluster_list_new1.size(); j++) {
      subcluster1 = subcluster_list_new1.at(j);
      if((subcluster2.n_num != 0)&&
          (subcluster1.index == subcluster2.index)) {
        subcluster_list_new2.push_back(subcluster2);
        break;
      }
    }
  }


  // store two basic clusters

  bool flag;
  flag = false;
  for(uint i=1; i<subcluster_list_new2.size(); i++) {
    if(subcluster_list_new2.at(i).index == subcluster_list1.at(size1-1).index) {
      flag = true;
    }
  }
  if(!flag) {
    subcluster_list_new2.push_back(subcluster_list1.at(size1-1));
  }

  flag = false;
  for(uint i=1; i<subcluster_list_new2.size(); i++) {
    if(subcluster_list_new2.at(i).index == subcluster_list2.at(size2-1).index) {
      flag = true;
    }
  }
  if(!flag) {
    subcluster_list_new2.push_back(subcluster_list2.at(size2-1));
  }

  sort(subcluster_list_new2.begin(), subcluster_list_new2.end(), &comparison_subcluster);

  //cerr << "Overlap\n";
  //for(uint i=0; i<subcluster_list_new2.size(); i++) {
  //    cerr << "n_num " << subcluster_list_new2.at(i).n_num << endl;
  //    cerr << "m_num " << subcluster_list_new2.at(i).m_num << endl;
  //    cerr << "index " << subcluster_list_new2.at(i).index << endl;
  //    cerr << endl;
  //}

  return subcluster_list_new2;


}

void cesubcluster::GetOverlapCluster() {
  // get the overlapped clusters of basic clusters
  // names of which are stored in base_name_list

  int base_index=0; if(base_index) {};
  vector< vector<Baker_str> > subcluster_list;
  vector<Baker_str> overlap_list, overlap_list_out;
  vector< vector<Baker_str> > subcluster_list_all;

  GetAllSubCluster();

  for(uint i=0; i<base_cluster.size(); i++) {
    // first get all subclusters and overlapped subcluster of equivalent basic cluster

    base_index = base_cluster.at(i);

    overlap_list = GetOverlapCluster(i+1);

    subcluster_list.push_back(overlap_list);

    overlap_list_out = MergeTwoSubCluster(overlap_list, overlap_list_out);
  }

  // second get overlapped subcluters of different basic cluster
  int num=2, total_num, size;
  total_num = base_cluster.size();
  vector<int> config;
  size = CombinationNr(num, total_num);
  for(int i=0; i<size; i++) {
    if(i==0) {
      config = AllCombination41(num, total_num, i+1);
    } else {
      config = AllCombination42(num, total_num, config);
    }

    overlap_list = GetOverlapCluster( sub_cluster.at(config.at(0)),
        sub_cluster.at(config.at(1)) );
    overlap_list_out = MergeTwoSubCluster(overlap_list, overlap_list_out);
  }

  sort(overlap_list_out.begin(), overlap_list_out.end(), &comparison_subcluster);

  overlap_cluster = overlap_list_out;

  // output subcluster_list
  cerr << "all overlap cluster size " << overlap_cluster.size()<<endl;
  for(uint i=0; i<overlap_cluster.size(); i++) {
    cerr << "n_num " << overlap_cluster.at(i).n_num << endl;
    cerr << "m_num " << overlap_cluster.at(i).m_num << endl;
    cerr << "index " << overlap_cluster.at(i).index << endl;
    cerr << endl;
  }
  cerr << endl;

}

vector<Baker_str> cesubcluster::MergeTwoSubCluster(
    const vector<Baker_str> & subcluster_list1,
    const vector<Baker_str> & subcluster_list2) {
  // Merge subcluster lists of two basic clusters

  vector<Baker_str> subcluster_list_new;
  subcluster_list_new = subcluster_list1;
  Baker_str subcluster_tmp;
  bool flag;

  for(uint i=0; i<subcluster_list2.size(); i++) {
    subcluster_tmp = subcluster_list2.at(i);
    flag = false;
    for(uint j=0; j<subcluster_list_new.size(); j++) {
      if(is_equal(subcluster_tmp, subcluster_list_new.at(j))) {
        flag = true;
        break;
      }
    }
    if(!flag) {
      subcluster_list_new.push_back(subcluster_tmp);
    }
  }

  sort(subcluster_list_new.begin(), subcluster_list_new.end(), &comparison_subcluster);

  //cerr << "After Merge\n";
  //for(uint i=0; i<subcluster_list_new.size(); i++) {
  //    cerr << "n_num " << subcluster_list_new.at(i).n_num << endl;
  //    cerr << "m_num " << subcluster_list_new.at(i).m_num << endl;
  //    cerr << "base cluster " << subcluster_list_new.at(i).base_cluster.pair_name << endl;
  //    cerr << "sub cluster " << subcluster_list_new.at(i).sub_cluster.pair_name << endl;
  //    cerr << endl;
  //}

  return subcluster_list_new;


}



bool comparison_subcluster(const Baker_str& subcluster1,
    const Baker_str& subcluster2) {
  // comparison of two clusters
  // "smaller" cluster has "small" pair_name
  // sorted by descending order

  return (subcluster1.index < subcluster2.index);
}



bool is_equal(const Baker_str & cluster1, const Baker_str & cluster2) {
  return (cluster1.index == cluster2.index);
}

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ceSL.cpp


ceSL::ceSL() : cestructure() {
  N.clear();
  config.clear();
  base.clear();
  name.clear();
  base_atoms.clear();
  SQS_figure_of_merit = -1.0;
  SQScompare_cluster_list.clear();
}

ceSL::ceSL(string & str_name) : cestructure(str_name, 0.0, 0.0) {

  for(uint i=0; i<name.size(); i++) {
    base.push_back(tolower(name.at(i)));
  }

  base = str_name;
  name.clear();
  base_atoms.clear();

  if(!(base == "fcc" || base == "bcc" || base =="hcp" )) {
    throw aurostd::xerror(__AFLOW_FILE__,"ceSL::ceSL():","structure must be one of fcc/bcc/hcp",_INPUT_ILLEGAL_); //CO20200624
  }
}

ceSL::ceSL(string & str_name, double stoich_b, 
    double fit_quantity) : cestructure( str_name, stoich_b, fit_quantity) 
{
  for(uint i=0; i<name.size(); i++) {
    base.push_back(tolower(name.at(i)));
  }

  base = str_name;
  name.clear();
  base_atoms.clear();

  if(!(base == "fcc" || base == "bcc" || base =="hcp" )) {
    throw aurostd::xerror(__AFLOW_FILE__,"ceSL::ceSL():","structure must be one of fcc/bcc/hcp",_INPUT_ILLEGAL_); //CO20200624
  }

}

ceSL::ceSL(string & str_name, xstructure & xstr) : cestructure(xstr) {
  base = str_name;
  name.clear();
  base_atoms.clear();

  xstructure str_base; // get the structure of structure_type (fcc/bcc/hcp)
  str_base.coord_flag = _COORDS_CARTESIAN_;


  if(base == "fcc") {
    str_base=aflowlib::PrototypeLibraries(cerr, aurostd::utype2string(_FCC_BEGIN),"",LIBRARY_MODE_HTQC); // no parameters
  } else if(base == "bcc") {
    str_base=aflowlib::PrototypeLibraries(cerr, aurostd::utype2string(_BCC_BEGIN),"",LIBRARY_MODE_HTQC); // no parameters
  } else if(base == "hcp") {
    str_base=aflowlib::PrototypeLibraries(cerr, aurostd::utype2string(_HCP_BEGIN),"",LIBRARY_MODE_HTQC); // no parameters
  }

  str_base.FixLattices();
  xmatrix<double> N_mat(1, 1, 3, 3);
  xmatrix<int> N_in(1, 1, 3, 3);
  int cell_nr_in = xstr.atoms.size();

  N_mat = structure.lattice * inverse(str_base.lattice);

  for(int i=1; i<=3; i++) {
    for(int j=1; j<=3; j++) {
      N_in[i][j] = (int(ceil(N_mat[i][j])));
    }
  }

  GetStructure(N_in, cell_nr_in);

  xstr.CalculateSymmetryFactorGroup();

  vector<int> atom_config;
  _atom atom1, atom2;
  bool flag = true;

  for(uint i=0; i<structure.atoms.size(); i++) {
    atom1 = structure.atoms.at(i);
    flag = true;

    for(uint j=0; j<xstr.fgroup.size() && flag; j++) {
      atom2.fpos = xstr.fgroup.at(j).Uf*atom1.fpos + xstr.fgroup.at(j).ftau;
      for(uint k=0; k<xstr.atoms.size() && flag; k++) {
        if(is_equal_fpos(atom2, xstr.atoms.at(k))) {
          if(xstr.atoms.at(k).name == xstr.species.at(0)) {
            structure.atoms.at(i).name = "A";
          } else if(xstr.atoms.at(k).name == xstr.species.at(1)) {
            structure.atoms.at(i).name = "B";
            atom_config.push_back(i);
          }
          flag = false;
        }
      }
    }

  }

  SetConfig(atom_config);
  StructureToName();
}

ceSL::~ceSL() {
  //config.clear();
}

ostream & operator<<(ostream & os, const ceSL & cestr) {
  // entropy and free energy are not calculated by CVM
  // but by MC
  // therefore, they are not output here

  ios_base::fmtflags stat_old;
  stat_old=os.setf(ios_base::fixed, ios_base::floatfield);
  os.precision(8);

  os << setw(18)
    << cestr.name
    << setw(12)
    << cestr.stoich_b
    << setw(12)
    << cestr.energy;

  // restore to the old setting
  os.setf(stat_old, ios_base::floatfield);

  return os;
}

void ceSL::SetUp(xmatrix<int> N_in, int cell_nr_in, vector<int> & atom_config) {
  GetStructure(N_in, cell_nr_in);
  SetConfig(atom_config);
  StructureToName();
}

void ceSL::SetUp(string & SLname) {
  if(NameToStructure(SLname)) {
    xmatrix<int> N_mat(1, 1, 3, 3);
    for(uint i=0; i<N.size(); i++) {
      for(uint j=0; j<N.size(); j++) {
        N_mat[i+1][j+1] = N.at(i).at(j);
      }
    }
    int cell_nr_in = cell_nr;
    GetStructure(N_mat, cell_nr_in);
    SetConfig(config);
  } else {
    throw aurostd::xerror(__AFLOW_FILE__,"ceSL::SetUp():","cannot find superlattice with name "+SLname,_INPUT_MISSING_); //CO20200624
  }
}

void ceSL::GetStructure(xmatrix<int> Nmat_in, int cell_nr_in) {
  // get the basic structure of a set of superlattice

  xstructure str_base, str_SL; // get the structure of structure_type (fcc/bcc/hcp)
  str_base.coord_flag = _COORDS_CARTESIAN_;

  if(base == "fcc") {
    str_base=aflowlib::PrototypeLibraries(cerr, aurostd::utype2string(_FCC_BEGIN),"",LIBRARY_MODE_HTQC); // no parameters);
  } else if(base == "bcc") {
    str_base=aflowlib::PrototypeLibraries(cerr, aurostd::utype2string(_BCC_BEGIN),"",LIBRARY_MODE_HTQC); // no parameters);
  } else if(base == "hcp") {
    str_base=aflowlib::PrototypeLibraries(cerr, aurostd::utype2string(_HCP_BEGIN),"",LIBRARY_MODE_HTQC); // no parameters);
  }

  //cerr << "Generate superlattice from base structure type " << base << endl;
  //cerr << str_base.lattice << endl;

  str_base.FixLattices();

  base_atoms = str_base.atoms;

  int atom_nr = cell_nr_in*str_base.atoms.size();
  cell_nr = cell_nr_in;

  _atom atom_tmp;

  str_SL.lattice = aurostd::xdouble(Nmat_in)*str_base.lattice;

  N.clear();
  for(int i=1; i<4; i++) {
    vector<int> N_tmp;
    for(int j=1; j<4; j++) {
      N_tmp.push_back(Nmat_in[i][j]);
    }
    N.push_back(N_tmp);
  }

  //cerr << endl;
  //cerr << str_base.lattice << endl;
  //cerr << endl;
  //cerr << str_SL.lattice << endl;
  //cerr << endl;
  //cerr << inverse(xdouble(Nmat_in)) << endl;
  //cerr << endl;

  xmatrix<double> Ninv(1, 1, 3, 3);
  Ninv = inverse(xdouble(Nmat_in));

  //cerr << "atom_nr " << atom_nr << endl;
  int Nmax = atom_nr-1;
  int atom_count = 0;

  // atom at origin is always in the base
  for(int l1 = 1; l1 < 4; l1 ++) {
    atom_tmp.fpos[l1] = 0.0;
    atom_tmp.cpos[l1] = 0.0;
  }

  str_SL.atoms.push_back(atom_tmp);

  //cerr << "Nmax " << Nmax << endl;
  //cerr << Ninv << endl;
  for(int i1=-Nmax; i1 < Nmax + 1; i1++) {
    for(int i2=-Nmax; i2 < Nmax + 1; i2++) {
      for(int i3=-Nmax; i3 < Nmax + 1; i3++) {
        bool flag = true;

        // in fractional coordinates
        xmatrix<double> x(1, 1, 3, 3);
        for(int l1=1; l1 < 4; l1++) {
          x[1][l1] = double(i1)*Ninv[1][l1];
          x[2][l1] = double(i2)*Ninv[2][l1];
          x[3][l1] = double(i3)*Ninv[3][l1];
        }

        double x1=0.0, x2=0.0, x3=0.0;
        for(int l1=1; l1 < 4; l1++) {
          x1 += x[l1][1];
          x2 += x[l1][2];
          x3 += x[l1][3];
        }


        flag = flag && (1.0 - abs(x1) > _EQUAL_DOUBLE)
          && (1.0 - abs(x2) > _EQUAL_DOUBLE)
          && (1.0 - abs(x3) > _EQUAL_DOUBLE) ;

        if(flag) {

          for(uint k=0; k<str_base.atoms.size(); k++) {
            for(int l1=1; l1 < 4; l1++) {
              atom_tmp.cpos[l1] = x1*str_SL.lattice[1][l1]
                + x2*str_SL.lattice[2][l1]
                + x3*str_SL.lattice[3][l1];
              atom_tmp.name = "A";
            }

            atom_tmp.fpos = C2F(str_SL.lattice, atom_tmp.cpos);

            for(uint l1=1; l1<4; l1++) {
              if(atom_tmp.fpos[l1] < 0) {
                atom_tmp.fpos[l1] +=
                  abs(int(atom_tmp.fpos[l1])) + 1.0;
              }
            }

            atom_tmp.cpos = F2C(str_SL.lattice, atom_tmp.fpos);

            bool flag_nodup = true;
            for(uint l1=0; l1<str_SL.atoms.size(); l1++) {
              // no duplicated atoms in base atoms

              if(is_equal(atom_tmp, str_SL.atoms.at(l1))) {
                flag_nodup = false;
                break;
              }
            }

            if(flag_nodup) {
              str_SL.atoms.push_back(atom_tmp);
              ++atom_count;
            }
          }

          if(atom_count == atom_nr) {
            goto out_of_loop;
          }

          //cerr << "atom " << atom_tmp << endl;
          //cerr << "cpos " << atom_tmp.cpos[1];
          //cerr << " " << atom_tmp.cpos[2];
          //cerr << " " << atom_tmp.cpos[3] << endl;

        }

      }
    }
  }

  //cerr << "atom count " << atom_count << endl;

out_of_loop: ;

             structure.SetStructure(str_SL);
}

void ceSL::SetConfig(vector<int> & atom_config) {
  // set the configurate in the superlattice
  // for binary

  config = atom_config;

  // reset all atom to type A
  for(uint k=0; k<structure.atoms.size(); k++)
  {
    structure.atoms.at(k).name = "A";
  }

  for(uint k=0; k<atom_config.size(); k++) {
    structure.atoms.at(atom_config.at(k)).name = "B";
  }

  int total_atom = structure.atoms.size();

  structure.num_each_type.clear();
  structure.num_each_type.push_back(total_atom-atom_config.size());
  structure.num_each_type.push_back(atom_config.size());
  // structure.comp_each_type.clear();
  // structure.comp_each_type.push_back((double) total_atom-atom_config.size());
  // structure.comp_each_type.push_back((double) atom_config.size());

  stoich_b = double(atom_config.size())/double(total_atom);
}

void ceSL::StructureToName() {
  // get the name
  // For binary, 0 stands for atom "A" and 1 for "B"
  // and the name is the number in decimal base
  //
  // vector<int> config stores the position of B type atoms
  //

  int sum=0;
  for(uint i=0; i< config.size(); i++) {
    sum += int(pow(2.0, double(config.at(i))));
  }

  string prefix = base+"_bSL";
  name.clear();
  name = prefix + name;
  for(uint i=0; i<N.size(); i++) {
    for(uint j=0; j<N.at(0).size(); j++) {
      name = name+aurostd::utype2string(N.at(i).at(j));
      name.push_back('_');
    }
  }
  name = name+aurostd::utype2string(cell_nr);
  name.push_back('_');
  name.push_back('0');
  name.push_back('x');
  name = name+IntToHexString(sum);

}

bool ceSL::NameToStructure(string & SLname) {
  string prefix = "SL";
  bool flag_fmt = true;
  uint name_size = SLname.size();
  size_t SL_pos;
  SL_pos = SLname.find("SL");
  if(( uint(SL_pos) > name_size )
      && (SLname.find("0x") > name_size )
      && (SLname.find("_") > name_size )
    ) {
    flag_fmt = false;
  }

  if(!flag_fmt) {
    cerr << "ceSL::NameToStructure: Superlattice name must have the format "
      << "xSL(n11)_(n12)_(n13)_(n21)_(n22)_(n23)_(n32)_(n32)_(n33)_(cell_num)_0x(xxxx) \n";
    return false;
  }

  string name_tmp = SLname;

  name_tmp.erase(0, int(SL_pos)+2);

  string num;
  N.clear();
  vector<int> Ns;
  for(int j=0; j<3; j++) {
    Ns.clear();
    for(int j2=0; j2<3; j2++) {

      num.clear();
      for(uint i=0; i< name_tmp.size(); i++) {
        if(name_tmp.at(i) != '_') {
          num.push_back(name_tmp.at(i));
        } else {
          name_tmp.erase(0, i+1);
          break;
        }
      }

      Ns.push_back(aurostd::string2utype<int>(num));
    }
    N.push_back(Ns);

  }

  num.clear();
  for(uint i=0; i< name_tmp.size(); i++) {
    if(name_tmp.at(i) != '_') {
      num.push_back(name_tmp.at(i));
    } else {
      name_tmp.erase(0, i+1);
      break;
    }
  }
  cell_nr = aurostd::string2utype<int>(num);

  name_tmp.erase(0,2); // delete 0x

  config.clear();
  int pos = -1;

  for(int i=name_tmp.size()-1; i>=0; i--) {

    int n;
    char ch = name_tmp.at(i);
    num.clear();
    num.push_back(ch);
    switch (ch) {
      case 'A' :
      case 'a' :
        n = 10;
        break;
      case 'B':
      case 'b':
        n = 11;
        break;
      case 'C':
      case 'c':
        n = 12;
        break;
      case 'D':
      case 'd':
        n = 13;
        break;
      case 'E':
      case 'e':
        n = 14;
        break;
      case 'F':
      case 'f':
        n = 15;
        break;
      default:
        n = aurostd::string2utype<int>(num);
    }

    for(int j=0; j<4; j++) { // which position has value 1
      int val = n%2;
      ++pos;

      if(val == 1) {
        config.push_back(pos);
      }

      n = n/2;
    }
  }

  name = SLname;

  stoich_b = double(config.size())/double(cell_nr);

  //cerr << "SL name " << name << endl;

  //cerr << "config \n";
  //for(uint i=0; i<config.size(); i++) {
  //    cerr << config.at(i) << " ";
  //}
  //cerr << endl;


  //cerr << "N1 " << N1 << " ";
  //cerr << "N2 " << N2 << " ";
  //cerr << "N3 " << N3 << " ";
  //cerr << endl;


  return flag_fmt;
}


string IntToHexString(int from) {
  std::stringstream temp;
  temp.precision(20);
  temp << std::hex << from;
  string to;
  temp >> to;
  return to;
}

void ceSL::PrintStructure(ostream & os) {

  os << name << endl;
  os << "trans_matrix = ";
  for(uint i1=0; i1<N.size(); i1++) {
    for(uint i2=0; i2<N.at(0).size(); i2++) {
      os << setw(4)
        << N.at(i1).at(i2);
    }
    os << " | ";
  }
  os << endl;
  os << "cell num = " << cell_nr << endl;
  os << "Stoich_b = " << stoich_b << endl;
  //os << structure.GetVolume() << endl;
  //os << "Superlattice lattice vector \n";
  ios_base::fmtflags old_stat=os.setf(ios_base::fixed, ios_base::floatfield);
  os.precision(14);
  for(int i=1; i<=3; i++) {
    os
      << setw(19)
      << structure.lattice[i][1]
      << setw(19)
      << structure.lattice[i][2]
      << setw(19)
      << structure.lattice[i][3]
      << endl;
  }
  os << setw(4)
    << structure.num_each_type.at(0)
    << setw(4)
    << structure.num_each_type.at(1)
    << endl;
  //os << "Catesian(" << structure.atoms.size() <<")" << endl;
  //for(uint i=0; i<structure.atoms.size(); i++) {
  //    os << setw(19)
  //        << structure.atoms.at(i).cpos[1]
  //        << setw(19)
  //        << structure.atoms.at(i).cpos[2]
  //        << setw(19)
  //        << structure.atoms.at(i).cpos[3]
  //        << setw(4)
  //        << structure.atoms.at(i).name
  //        << endl;
  //}

  sort(structure.atoms.begin(), structure.atoms.end(), &comparison_atomtype);
  os << "Direct(" << structure.atoms.size() <<")" << endl;
  for(uint i=0; i<structure.atoms.size(); i++) {
    os << setw(19)
      << structure.atoms.at(i).fpos[1]
      << setw(19)
      << structure.atoms.at(i).fpos[2]
      << setw(19)
      << structure.atoms.at(i).fpos[3]
      << setw(4)
      << structure.atoms.at(i).name
      << endl;
  }
  os.precision(_DEFAULT_PRECISION);
  os.setf(old_stat, ios_base::floatfield);

}

void ceSL::PrintStructure(ostream & os, vector<_ceatom> atom_species) {

  // only valid for binary alloy
  // because only two atom species are specified here

  int A_atom_nr, B_atom_nr;
  B_atom_nr = int(stoich_b * cell_nr);
  A_atom_nr = cell_nr - B_atom_nr;

  double vol;
  vol = A_atom_nr*atom_species.at(0).volume + B_atom_nr*atom_species.at(1).volume;
  vol = (vol > 0) ? -vol : vol;

  os.precision(14);
  ios_base::fmtflags old_stat=os.setf(ios_base::fixed, ios_base::floatfield);

  os << name << endl;
  os << vol << endl;
  for(int i=1; i<=3; i++) {
    os
      << setw(19)
      << structure.lattice[i][1]
      << setw(19)
      << structure.lattice[i][2]
      << setw(19)
      << structure.lattice[i][3]
      << endl;
  }
  os << setw(4)
    << structure.num_each_type.at(0)
    << setw(4)
    << structure.num_each_type.at(1)
    << endl;
  //os << "Catesian(" << structure.atoms.size() <<")" << endl;
  //for(uint i=0; i<structure.atoms.size(); i++) {
  //    os << setw(19)
  //        << structure.atoms.at(i).cpos[1]
  //        << setw(19)
  //        << structure.atoms.at(i).cpos[2]
  //        << setw(19)
  //        << structure.atoms.at(i).cpos[3]
  //        << setw(4)
  //        << structure.atoms.at(i).name
  //        << endl;
  //}

  sort(structure.atoms.begin(), structure.atoms.end(), &comparison_atomtype);

  os << "Direct(" << structure.atoms.size() <<")" << endl;
  for(uint i=0; i<structure.atoms.size(); i++) {
    string name;
    if(structure.atoms.at(i).name == "A") {
      name = atom_species.at(0).name;
    } else if(structure.atoms.at(i).name == "B") {
      name = atom_species.at(1).name;
    }
    os << setw(19)
      << structure.atoms.at(i).fpos[1]
      << setw(19)
      << structure.atoms.at(i).fpos[2]
      << setw(19)
      << structure.atoms.at(i).fpos[3]
      << setw(4)
      << name
      << endl;
  }
  os.precision(_DEFAULT_PRECISION);
  os.setf(old_stat, ios_base::floatfield);

}

double ceSL::IsSQS(int site_num, int NNNum) {

  // Special Quasirandom Structure (SQS)
  // Zunger, Wei, Ferreira, and Bernard, PRL, 65, 353(1990)
  // Output a figure of merit defined as \Sum weight*(cor_i - random_cor_i)^2
  // the smaller the is figure of merit, the better is the chance that
  // the structure is an SQS
  // check upto the correlation of "site_num" cluster with largest Nearest
  // neighbor number of value "NNNum"

  SQScluster cmp_cluster;
  SQS_figure_of_merit = 0.0;
  SQScompare_cluster_list.clear();

  // get the list of clusters which correlations will be compared to
  // the ones in random structure
  // The clusters are those with number of site no larger than "site_num"
  // and the largest NNNum equal or smaller than "NNNum"
  // (0 0 0 0) is excluded

  cecluster cluster_iter;

  for(uint i=0; i< str_cluster->rep_cluster.size(); i++) {
    cluster_iter = str_cluster->rep_cluster.at(i);
    if(cluster_iter.site_num <= site_num &&
        cluster_iter.site_num > 1) {
      // exclude cluster (0 0 0 0) and (1 0 1 1)
      // as their correlations are always the same in both
      // SQS and random structure
      if(cluster_iter.NNNum <= NNNum) {
        cmp_cluster.pair_name = cluster_iter.pair_name;
        cmp_cluster.index = cluster_iter.index;
        cmp_cluster.site_num = cluster_iter.site_num;
        SQScompare_cluster_list.push_back(cmp_cluster);
      }
    }
  }

  // get the correlations and set weights
  vector<int>::iterator iter;
  vector<int> cal_cor_list; // cluster whose correlations needed to be calculated
  vector<int> index_list;
  int size_tmp = SQScompare_cluster_list.size();
  for(uint i=0; i<SQScompare_cluster_list.size(); i++) {
    // whether correlation is needed to be calculated
    int index = SQScompare_cluster_list.at(i).index;
    iter = find(allECI_cluster.begin(), allECI_cluster.end(), index);
    if(iter < allECI_cluster.end()) {
      // correlation has already been calculated
      int pos = iter - allECI_cluster.begin();
      if(pos < int(allECI_correlation.size())) {
        SQScompare_cluster_list.at(i).correlation = allECI_correlation.at(pos);
      } else {
        cal_cor_list.push_back(i);
        index_list.push_back(index);
      }
    } else {
      // not be calculated
      cal_cor_list.push_back(i);
      index_list.push_back(index);
    }

    // get the correlations in random strucutre
    double tmp = (1.0 - stoich_b)*_A_ATOM + stoich_b*_B_ATOM;
    double cor = pow(tmp, SQScompare_cluster_list.at(i).site_num);
    SQScompare_cluster_list.at(i).correlation_random_structure = cor;

    // set the weight
    // For N's cluster, the cluster with index n has a weight N-n
    // better algorithm needed
    // in the futhur version, weighs may be able to be given
    // by user
    // Let (2 1 1 2) to be the most important cluster
    //SQScompare_cluster_list.at(i).weight =
    //    (i==0) ? pow(size_tmp, 3): pow((size_tmp - i -1), 2);
    SQScompare_cluster_list.at(i).weight = pow(10.0, (size_tmp - i ));
  }

  vector<double> cor_list;
  if(cal_cor_list.size() > 0) {
    cor_list = GetCorrelation(index_list);
    for(uint i=0; i<cal_cor_list.size(); i++) {
      int cluster_index = cal_cor_list.at(i);
      SQScompare_cluster_list.at(cluster_index).correlation =
        cor_list.at(i);
    }
  }

  // get the figure of merit

  for(uint i=0; i<SQScompare_cluster_list.size();i++) {
    double cor = SQScompare_cluster_list.at(i).correlation;
    double cor_random = SQScompare_cluster_list.at(i).correlation_random_structure;
    double weight = SQScompare_cluster_list.at(i).weight;

    SQS_figure_of_merit += (cor-cor_random)*(cor-cor_random)*weight;
  }


  return SQS_figure_of_merit;

}

void ceSL::OutputSQS(ostream & oss) {

  oss << setw(20)
    <<"Strucutre "
    << setw(12)
    << name
    << endl;
  oss << setw(20)
    << "concentration"
    << setw(20)
    << stoich_b << endl;
  oss << setw(20)
    << "atom number"
    << setw(20)
    << cell_nr << endl;
  oss << "-------------------------------------------------------------\n";
  oss << setw(10)
    << "cluster"
    << setw(30)
    << "          correlation"
    << endl;

  oss << "         "
    << setw(20)
    << "CE::SQS"
    << setw(20)
    << "random"
    << endl;
  oss << "-------------------------------------------------------------\n";
  for(uint i=0; i<SQScompare_cluster_list.size(); i++) {
    oss << setw(10)
      << SQScompare_cluster_list.at(i).pair_name
      << setw(20)
      << SQScompare_cluster_list.at(i).correlation
      << setw(20)
      << SQScompare_cluster_list.at(i).correlation_random_structure
      << endl;
  }

}

bool ceSL::IsPrimitiveCell() {

  // use the criterion proposed by Ferreria, Wei, and Zunger

  bool flag = true;
  _atom atom1, atom2, atom3, atom4;

  //cerr << "atom list\n";
  //for(uint i=0; i < structure.atoms.size(); i++) {
  //    atom1 = structure.atoms.at(i);
  //        cerr << atom1.fpos[1]
  //            << " "
  //            << atom1.fpos[2]
  //            << " "
  //            << atom1.fpos[3];
  //    cerr << " " << atom1.name;
  //cerr << endl;
  //}
  //cerr << endl;

  for(uint i=0; i < structure.atoms.size(); i++) {

    atom1 = structure.atoms.at(i); // check whether atoms.at(j) is a translation
    if(modulus(atom1.fpos) ==0) {
      break;
    }

    //cerr << "atom 1 ";
    //    cerr << atom1.fpos[1]
    //        << " "
    //        << atom1.fpos[2]
    //        << " "
    //        << atom1.fpos[3];
    //cerr << " " << atom1.name << endl;

    bool flag1 = true;
    for(uint j=0; j < structure.atoms.size(); j++) {

      // check for all base atoms

      atom2 = structure.atoms.at(j);
      atom4 = atom2;

      //cerr << "old atom2 " ;
      //cerr << atom2.fpos[1]
      //    << " "
      //    << atom2.fpos[2]
      //    << " "
      //    << atom2.fpos[3];
      //    cerr << atom2.name << endl;

      atom4.fpos = atom1.fpos + atom4.fpos;


      for(int l1=1; l1<4; l1++) {
        atom4.fpos[l1] = atom4.fpos[l1] - int(atom4.fpos[l1]+_EQUAL_DOUBLE*0.1);
      }

      //cerr << atom4.fpos[1]
      //    << " "
      //    << atom4.fpos[2]
      //    << " "
      //    << atom4.fpos[3]
      //    << endl;

      for(uint l1=1; l1<4; l1++) {
        if(atom4.fpos[l1] < 0) {
          atom4.fpos[l1] += abs(int(atom4.fpos[l1])) + 1.0;
        }
      }

      //cerr << atom4.fpos[1]
      //    << " "
      //    << atom4.fpos[2]
      //    << " "
      //    << atom4.fpos[3]
      //    << endl;

      // get the new atom's type

      for(uint l1=0; l1<structure.atoms.size();l1++) {
        atom3 = structure.atoms.at(l1);
        if(( modulus(atom4.fpos - atom3.fpos)) < _EQUAL_DOUBLE) {
          atom4.name = atom3.name;
          break;
        }
      }

      // cerr << "atom4 " << atom4.name << "  atom2 " << atom2.name << endl;

      // flag1 = true : not primitive unit cell

      flag1 = flag1 &&(atom4.name == atom2.name );
      cerr << "flag1 " << flag1 << endl;
    }

    if(flag1) {
      flag = false; // not be a primitive unit cell
      break;
    }
  }

  return flag; // default value flag = true, be a primitive unit cell
}

vector<int> GetDivisors(int num) {
  vector<int> divisors;

  // calculate all divisor of a given integer num
  for(int i = num; i>=1; i--) {
    if(num % i == 0) {
      divisors.push_back(i);
    }
  }

  return divisors;

}

bool comparison_atomtype(_atom atom1, _atom atom2) {
  return(atom1.name < atom2.name );
}


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
