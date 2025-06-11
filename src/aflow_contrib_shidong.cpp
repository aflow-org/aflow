// ***************************************************************************
// *                                                                         *
// *               AFlow SHIDONG WANG - Duke University 2010-2011            *
// *                                                                         *
// ***************************************************************************
// aflow_contrib_shidong_funs.cpp
// functions written by
// 2010-2011: shidong.wang@duke.edu

#ifndef _AFLOW_CONTRIB_SHIDONG_FUNCS_CPP_
#define _AFLOW_CONTRIB_SHIDONG_FUNCS_CPP_

#include "aflow_contrib_shidong.h"
#include "aflow.h"

// ***************************************************************************
// check the validity of the structure type, which should be one of
// fcc, bcc, hcp
// ***************************************************************************
bool ACEGetStructureType(const string& structure_type_in) {
  bool ACEFlag;

  string structure_type;
  structure_type = structure_type_in;

  for (uint i = 0; i < structure_type.size(); i++) {
    structure_type[i] = tolower(structure_type[i]);
  }

  if( structure_type == "fcc" || structure_type == "bcc" || structure_type == "hcp") {
    ACEFlag = true;
  } else {
    throw aurostd::xerror(__AFLOW_FILE__,"ACEGetStructureType():","The structure given must be one of fcc[FCC]/bcc[BCC]/hcp[HCP]",_INPUT_ILLEGAL_); //CO20200624
  }

  return ACEFlag;
}

// ***************************************************************************
// Get the data of initial training set
// File format
// alloy_name structure_name x_b property_to_be_fit
//
// ***************************************************************************
vector<cestructure>  ReadInFitStructure(istream & ins, string & structure_type) {
  vector<cestructure> str_out_list;
  uint LineCounter = 1;
  string line_content;
  string alloy_name, structure_name;
  double stoich_b, fit_quantity;
  cestructure str_tmp;

  stringstream ss;

  //while ( getline(cin, line_content)  )
  while ( getline(ins, line_content)  ) {
    if( line_content.at(0) != _COMMENT_CHAR ) {
      ss.clear();
      ss << line_content;
      ss >> alloy_name >> structure_name >> stoich_b
        >> fit_quantity;

      if( structure_name.at(0) == structure_type[0] ) {
        // type begin with a letter

        if( int(structure_name.size()) < _SLname_length
            || structure_name.substr(_SLname_pos, 2) != "SL" ) {
          str_tmp = cestructure(structure_name, stoich_b,
              fit_quantity);
          str_out_list.push_back(str_tmp);
        } else {
          ceSL str_SL_tmp;
          string str_type = structure_type;
          str_SL_tmp = ceSL(str_type, stoich_b,
              fit_quantity);
          str_SL_tmp.SetUp(structure_name);
          str_out_list.push_back(str_SL_tmp);
        }

      } else  if(0 <= structure_name.at(0) &&  structure_name.at(0) <= '9') {
        // type present as number
        int typenum = atoi(structure_name.c_str());
        switch (structure_type[0]) {
          case 'f':
            if( _FCC_BEGIN <= typenum && typenum <= _FCC_END ) {
              str_tmp = cestructure(structure_name, stoich_b,
                  fit_quantity);
              str_out_list.push_back(str_tmp);
            }
            break;
          case 'b':
            if( _BCC_BEGIN <= typenum && typenum <= _BCC_END ) {
              str_tmp = cestructure(structure_name, stoich_b,
                  fit_quantity);
              str_out_list.push_back(str_tmp);
            }
            break;
          case 'h':
            if( _HCP_BEGIN <= typenum && typenum <= _HCP_END ) {
              str_tmp = cestructure(structure_name, stoich_b,
                  fit_quantity);
              str_out_list.push_back(str_tmp);
            }
            break;
        }
      }

    }
    LineCounter++;
    //cout << LineCounter << " " << structure_name << endl;
  }


  if( str_out_list.size() == 0) {
    throw aurostd::xerror(__AFLOW_FILE__,"ReadInFitStructure():","No "+structure_type+" structure is found. Check your input file",_FILE_CORRUPT_); //CO20200624
  }

  return str_out_list;
}


void WriteOutFitStructure(ostream & os, string & alloy_name,
    vector<cestructure> & str_list) {

  for (uint i=0; i<str_list.size(); i++) {
    os << "(" << alloy_name << ") " << str_list.at(i) << endl;
  }
}

//*************************************
// Convex Hull functions
//*************************************
vector<int> GroundStateCandidate(vector< vector<double> > & points) {

  vector<int> ground_states;

  // first sort the input list
  // add the index of the old list as the third values
  vector< vector<double> > points_sorted;

  points_sorted = points;
  for (uint i=0; i<points.size(); i++) {
    points_sorted.at(i).push_back(i);
  }

  sort(points_sorted.begin(), points_sorted.end(), &comparison_points);

  //cerr << endl;
  //for (uint i=0; i< points_sorted.size(); i++) {
  //    cerr.precision(8);
  //    cerr << setw(18)
  //        <<points_sorted.at(i).at(0)
  //        << setw(18)
  //        <<points_sorted.at(i).at(1)
  //        << setw(18)
  //        <<points_sorted.at(i).at(2)
  //        << endl;
  //}
  //cerr << endl;

  int pt_min;
  pt_min = 0;

  //const double _EQUAL_DB_HERE = 1.0e-6;
  double _EQUAL_NDB_HERE = -1.0e-6;
  double _EQUAL_DB_HERE = 1.0e-6;
  for (uint i=0; i<points_sorted.size(); i++) {
    if( abs(points_sorted.at(i).at(0) - points_sorted.at(pt_min).at(0)) < _EQUAL_DB_HERE ) {
      //if( points_sorted.at(i).at(1) < points_sorted.at(pt_min).at(1) )
      if( points_sorted.at(i).at(1) - points_sorted.at(pt_min).at(1)
          < _EQUAL_NDB_HERE ) {
        pt_min = i;
      }
    } else {

      ground_states.push_back(int(points_sorted.at(pt_min).at(2)));
      pt_min = i;
    }
  }

  if(pt_min == int(points_sorted.size())-1 ) {
    ground_states.push_back(int(points_sorted.at(pt_min).at(2)));
  }

  //cerr.precision(10);
  //cerr << endl;
  //for (uint i=0; i<ground_states.size(); i++) {
  //    int index = ground_states.at(i);
  //    cerr << setw(18)
  //        <<points_sorted.at(index).at(0)
  //        << setw(18)
  //        <<points_sorted.at(index).at(1)
  //        << setw(18)
  //        <<points_sorted.at(index).at(2)
  //        << endl;
  //}
  //cerr << endl;

  return ground_states;

}

vector<int> CEConvexHull(vector< vector<double> > & points, const string & option) {
  // calculate the convex hull of input point list with size [N, 2]
  // output the index of those points in convex hull
  // option = {min, max, all}
  // min: get hull of minima
  // max: get hull of maxima
  // all: get hull of both minima and maxima


  vector<int> hull_indices;

  // first sort the input list
  // add the index of the old list as the third values
  vector< vector<double> > points_sorted;

  points_sorted = points;
  for (uint i=0; i<points.size(); i++) {
    points_sorted.at(i).push_back(i);
  }

  sort(points_sorted.begin(), points_sorted.end(), &comparison_points);

  //cerr.precision(8);
  //for (uint i=0; i< points_sorted.size(); i++) {
  //    cerr << setw(12)
  //        <<points_sorted.at(i).at(0)
  //        << setw(12)
  //        <<points_sorted.at(i).at(1)
  //        << setw(12)
  //        <<points_sorted.at(i).at(2)
  //        << endl;
  //}

  // get all minima and maxima
  int pt_min, pt_max;
  pt_min = 0;
  pt_max = 0;
  vector<int> minima, maxima;

  double _EQUAL_NDB_HERE = -1.0e-6;
  double _EQUAL_DB_HERE = 1.0e-6;
  for (uint i=1; i<points_sorted.size(); i++) {
    if( abs(points_sorted.at(i).at(0) - points_sorted.at(pt_min).at(0))
        < _EQUAL_DB_HERE ) {
      if( points_sorted.at(i).at(1) - points_sorted.at(pt_min).at(1)
          < _EQUAL_NDB_HERE ) {
        pt_min = i;
      } else if(
          points_sorted.at(pt_min).at(1) -
          points_sorted.at(i).at(1) < _EQUAL_NDB_HERE
          ) {
        pt_max = i;
      }
    } else {
      minima.push_back(pt_min);
      maxima.push_back(pt_max);
      pt_min = i;
      pt_max = i;
    }
  }

  if(pt_min == int(points_sorted.size())-1 ) {
    minima.push_back(pt_min);
  }
  if(pt_max == int(points_sorted.size())-1 ) {
    maxima.push_back(pt_max);
  }

  //cerr << "minima\n";
  //for (uint i=0; i<minima.size(); i++) {
  //    int index = minima.at(i);
  //    cerr << setw(12)
  //        <<points_sorted.at(index).at(0)
  //        << setw(12)
  //        <<points_sorted.at(index).at(1)
  //        << setw(12)
  //        <<points_sorted.at(index).at(2)
  //        << endl;
  //}
  //cerr << "maxima\n";
  //for (uint i=0; i<maxima.size(); i++) {
  //    int index = maxima.at(i);
  //    cerr << setw(12)
  //        <<points_sorted.at(index).at(0)
  //        << setw(12)
  //        <<points_sorted.at(index).at(1)
  //        << setw(12)
  //        <<points_sorted.at(index).at(2)
  //        << endl;
  //}

  // calculate the hull
  // from definition and it is not efficient

  vector<int> hull_min_list;
  vector<int> hull_max_list;
  int hull = 0;
  hull_min_list.push_back(hull);


  for (uint j = 1; j<minima.size()-1; j++) {
    bool flag_hull = true;
    for (uint i=j+1; i<minima.size(); i++) {
      int pt1, pt2, pt3;
      pt1 = minima.at(hull);
      pt2 = minima.at(j);
      pt3 = minima.at(i);
      double aa, bb, aa1, bb1;
      aa = ( points_sorted.at(pt3).at(0) - points_sorted.at(pt1).at(0) );
      bb = ( points_sorted.at(pt3).at(1) - points_sorted.at(pt1).at(1) );
      aa1 = ( points_sorted.at(pt2).at(0) - points_sorted.at(pt1).at(0) );
      bb1 = points_sorted.at(pt1).at(1) + bb/aa*aa1;

      flag_hull = flag_hull && (bb1 > points_sorted.at(pt2).at(1) );
    }

    if(flag_hull) {
      hull = j;
      hull_min_list.push_back(hull);
    }
  }


  hull_min_list.push_back(minima.size()-1);

  //cerr << "hull points\n";
  //for (uint i=0; i<hull_min_list.size(); i++) {
  //    int index1, index2;
  //    index1 = hull_min_list.at(i);
  //    index2 = minima.at(index1);
  //    cerr << setw(12)
  //        <<points_sorted.at(index2).at(0)
  //        << setw(12)
  //        <<points_sorted.at(index2).at(1)
  //        << setw(12)
  //        <<points_sorted.at(index2).at(2)
  //        << endl;
  //}

  // maxima

  hull = 0;
  hull_max_list.push_back(hull);
  for (uint j = 1; j<maxima.size()-1; j++) {
    bool flag_hull = true;
    for (uint i=j+1; i<maxima.size(); i++) {
      int pt1, pt2, pt3;
      pt1 = maxima.at(hull);
      pt2 = maxima.at(j);
      pt3 = maxima.at(i);
      double aa, bb, aa1, bb1;
      aa = ( points_sorted.at(pt3).at(0) - points_sorted.at(pt1).at(0) );
      bb = ( points_sorted.at(pt3).at(1) - points_sorted.at(pt1).at(1) );
      aa1 = ( points_sorted.at(pt2).at(0) - points_sorted.at(pt1).at(0) );
      bb1 = points_sorted.at(pt1).at(1) + bb/aa*aa1;

      flag_hull = flag_hull && (bb1 < points_sorted.at(pt2).at(1) );
    }

    if(flag_hull) {
      hull = j;
      hull_max_list.push_back(hull);
    }
  }


  hull_max_list.push_back(maxima.size()-1);

  //cerr << "hull points\n";
  //for (uint i=0; i<hull_max_list.size(); i++) {
  //    int index1, index2;
  //    index1 = hull_max_list.at(i);
  //    index2 = maxima.at(index1);
  //    cerr << setw(12)
  //        <<points_sorted.at(index2).at(0)
  //        << setw(12)
  //        <<points_sorted.at(index2).at(1)
  //        << setw(12)
  //        <<points_sorted.at(index2).at(2)
  //        << endl;
  //}

  if( option == "min" || option == "all" ) {
    for (uint i=0; i<hull_min_list.size(); i++) {
      int index1, index2;
      index1 = hull_min_list.at(i);
      index2 = minima.at(index1);
      hull_indices.push_back(int(points_sorted.at(index2).at(2)));
    }

  }

  if( option == "max" || option == "all" ) {
    for (uint i=0; i<hull_max_list.size(); i++) {
      int index1, index2;
      index1 = hull_max_list.at(i);
      index2 = maxima.at(index1);
      hull_indices.push_back(int(points_sorted.at(index2).at(2)));
    }
  }


  //cerr << "hull points \n";
  //for (uint i=0; i<hull_indices.size(); i++) {
  //    int index = hull_indices.at(i);
  //    cerr << setw(12)
  //        <<points.at(index).at(0)
  //        << setw(12)
  //        <<points.at(index).at(1)
  //        << endl;
  //}

  return hull_indices;

}

bool comparison_points(const vector<double> & point1, const vector<double> & point2) {
  //return ( point1.at(0) - point2.at(0) < _EQUAL_DOUBLE );
  //return ( point1.at(0) < point2.at(0) );

  double _EQUAL_DB_HERE = 1.0e-6;
  if( abs(point1.at(0) - point2.at(0)) > _EQUAL_DB_HERE ) {
    return ( point1.at(0) < point2.at(0) );
  } else {
    return ( point1.at(2) < point2.at(2) );
  }
}

void ACEGroundStateConvexHull(ifstream & os, vector<cestructure> str_list) {

  // get ground states and convex hull
  // return the names of the ground state in hull file

  vector<vector<double> > alist;
  vector<double> blist;
  vector<string> name_list;
  vector<string> hull_name_list;
  for (uint i=0; i<str_list.size(); i++) {
    blist.clear();
    blist.push_back(str_list.at(i).StoichB());
    blist.push_back(str_list.at(i).Energy());
    alist.push_back(blist);
    name_list.push_back(str_list.at(i).Name());

  }

  vector< vector<double> > SLres_list;
  vector<double> res_tmp;
  vector<string> SLname_list;
  string SL_name;

  SLres_list = alist;
  SLname_list = name_list;

  double stoich_b, energy;
  while ( os >> SL_name >> stoich_b >> energy ) {

    res_tmp.clear();

    res_tmp.push_back(stoich_b);
    res_tmp.push_back(energy);

    SLres_list.push_back(res_tmp);
    SLname_list.push_back(SL_name);

  }

  // get points at hull
  vector<int> hull_pts;
  string option="min";
  hull_pts = CEConvexHull(SLres_list, option);

  // get the ground states
  vector<int> ground_state_candidates;
  ground_state_candidates = GroundStateCandidate(SLres_list);

  ofstream total_pt_file, hull_pt_file, groundstate_file, gnuplot_file;
  total_pt_file.open(_TOTALPTFILE.c_str());
  hull_pt_file.open(_HULLPTFILE.c_str());
  groundstate_file.open(_GNDPTFILE.c_str());

  //gnuplot_file.open(_GNUPLOTPTFILE.c_str());

  total_pt_file.setf(ios_base::fixed, ios_base::floatfield);
  hull_pt_file.setf(ios_base::fixed, ios_base::floatfield);
  total_pt_file.precision(6);
  hull_pt_file.precision(6);
  groundstate_file.precision(6);

  for (uint i=0; i<SLres_list.size(); i++) {
    total_pt_file
      << setw(18)
      << SLname_list.at(i)
      << setw(12)
      << SLres_list.at(i).at(0)
      << setw(12)
      << SLres_list.at(i).at(1)
      << endl;
  }
  for (uint i=0; i<hull_pts.size(); i++) {
    int index = hull_pts.at(i);
    hull_pt_file
      << setw(18)
      << SLname_list.at(index)
      << setw(12)
      << SLres_list.at(index).at(0)
      << setw(12)
      << SLres_list.at(index).at(1)
      << endl;
    hull_name_list.push_back(SLname_list.at(index));
  }
  for (uint i=0; i<ground_state_candidates.size(); i++) {
    int index = ground_state_candidates.at(i);
    groundstate_file << setw(18)
      << SLname_list.at(index)
      << setw(12)
      << SLres_list.at(index).at(0)
      << setw(12)
      << SLres_list.at(index).at(1)
      << endl;
  }

  total_pt_file.close();
  hull_pt_file.close();
  groundstate_file.close();

}

vector<string> GetNewHullState(vector<cestructure> & str_list) {

  vector<string> name_list;
  vector<cepair> state_list;
  cepair state_tmp;

  ifstream fin;
  fin.open(_HULLPTFILE.c_str());
  string name;
  double stoich_b, fit_quantity;

  //while ( !fin.eof() )
  while ( fin >> name >> stoich_b >> fit_quantity ) {

    bool flag = false;
    for (uint i=0; i < str_list.size(); i++) {
      if  (name == str_list.at(i).Name()) {
        flag = true;
        break;
      }
    }

    if( !flag ) {
      state_tmp.name = name;
      state_tmp.num = fit_quantity;

      state_list.push_back(state_tmp);
    }

    //getline(fin, name);
  }

  //state_list.pop_back();

  sort(state_list.begin(), state_list.end(), &comparison_pair);

  for (uint i=0; i<state_list.size(); i++) {
    name_list.push_back(state_list.at(i).name);
  }

  return name_list;
}

vector<string> GetNewGroundState(vector<cestructure> & str_list) {

  vector<string> name_list;
  vector<cepair> state_list;
  cepair state_tmp;

  ifstream fin;
  fin.open(_GNDPTFILE.c_str());
  string name;
  double stoich_b, fit_quantity;

  while ( fin >> name >> stoich_b >> fit_quantity) {

    bool flag = false;
    for (uint i=0; i < str_list.size(); i++) {
      if  (name == str_list.at(i).Name()) {
        flag = true;
        break;
      }
    }

    if( !flag ) {
      state_tmp.name = name;
      state_tmp.num = fit_quantity;

      state_list.push_back(state_tmp);
    }

  }

  sort(state_list.begin(), state_list.end(), &comparison_pair);

  for (uint i=0; i<state_list.size(); i++) {
    name_list.push_back(state_list.at(i).name);
  }

  return name_list;
}

bool comparison_pair(cepair pair1, cepair pair2) {
  return ( pair1.num < pair2.num );
}

// ***************************************************************************
// Get Optimal ECI's
// ***************************************************************************

void ACEOptimalECIClusters(vector<cestructure> & str_list,
    ceECIcluster & ECIcluster) {
  // obtain the optimal set of clusters for ECI and ECI's
  // store it in ECIcluster.ECI_cluster
  // and ECIcluster.ECI

  vector<int> ECIcluster_opt;

  // before optimizing
  cerr << "before optimazation: ECI and CVM cluster\n";
  ECIcluster.PrintOutECICluster();

  //ECIcluster_opt = StepWise(str_list, ECIcluster);
  ECIcluster_opt = GeneticAlgorithm(str_list, ECIcluster);

  ECIcluster.DeleteZeroECI();

  // after optimizing
  cerr << "after optimazation: ECI and CVM cluster\n";
  ECIcluster.PrintOutECICluster();

  // reset ECI, ECI cluster, correlations
  for (uint i=0; i<str_list.size(); i++) {
    str_list.at(i).SetECICluster(ECIcluster);
    str_list.at(i).GetECICorrelation();
  }
}

vector<int> StepWise( vector<cestructure> & str_list,
    ceECIcluster & ECIcluster) {
  // forward/backward greedy search algorithm

  double score_opt_global;
  vector<int> trial_set_opt_global;

  // set up the cross validating sets
  vector< vector<int> > str_leftout_set_list;
  int total_str = str_list.size();
  //int left_out_num = _LEFTOUT_NUM;
  //left_out_num = min(int(total_str*_TRIAL_RATIO)+1,left_out_num);
  int left_out_num = aurostd::min(int(total_str*_TRIAL_RATIO)+1,_LEFTOUT_NUM);

  str_leftout_set_list = ACESetUpCrossValidationSetsRandom(
      left_out_num, total_str);

  score_opt_global = 1.0e4;

  srand(time(0));

  for (int k=0; k<_INITIAL_NUM; k++) {
    // avoid local minimum
    // cluster pool is stored in ECIcluster.ECI_cluster
    ceECIcluster ECIcluster_tmp;
    ECIcluster_tmp = ECIcluster;

    // generate the initial trial set by randomly selection
    // some clusters

    int total_cluster = ECIcluster.ECICluster().size();
    int trial_set_size, leftout_set_size;
    vector<int> trial_set;

    //trial_set_size = rand()%total_cluster + 1 ;
    trial_set_size = aurostd::max(5, rand()%(total_cluster/2) + 1);
    cerr << "SetWise: initial trial set size " << trial_set_size << endl;

    double score_opt;
    vector<int> trial_set_opt;

    score_opt = 1.0e4;

    while (int(trial_set.size()) < trial_set_size) {
      int index = rand()%total_cluster;
      bool flag = true;
      for (uint i=0; i<trial_set.size(); i++) {
        if( index == trial_set.at(i) ) {
          flag = false;
          break;
        }
      }
      if(flag) {
        trial_set.push_back(index);
      }

    }


    sort(trial_set.begin(), trial_set.end());

    cerr << "Inital set \n";
    for (uint i=0; i<trial_set.size(); i++) {
      cerr << trial_set.at(i) << " ";
    }
    cerr << endl;

    ECIcluster_tmp.GetECICluster(trial_set);

    // get correlation functions of ECI cluster
    for (uint i=0; i<str_list.size(); i++) {
      str_list.at(i).SetECICluster(ECIcluster_tmp);
      str_list.at(i).GetECICorrelation();
      //str_list.at(i).PrintOutCorrelation();
    }


    double  score;
    vector<int> trial_set_orig, leftout_set;
    bool flag;

    score = ACECrossValidation(str_leftout_set_list,
        str_list, ECIcluster_tmp);

    //cerr << "str left out list size " << str_leftout_set_list.size() << endl;
    //for (uint l=0; l<str_leftout_set_list.size(); l++) {
    //    for (uint l1=0; l1<str_leftout_set_list.at(l).size(); l1++) {
    //        cerr << str_leftout_set_list.at(l).at(l1) << " ";
    //    }
    //    cerr << endl;
    //}
    //cerr << endl;


    if( score < score_opt ) {
      score_opt = score;
      trial_set_opt = trial_set;
    }


    flag = false; // initial value to go into the loop

    //bool forward_flag, backward_flag;
    //forward_flag = true;
    //backward_flag = true;

    while ( !flag ) {

      leftout_set_size = total_cluster - trial_set_size;


      trial_set_orig = trial_set;

      leftout_set=LeftOutSet(total_cluster, trial_set);

      bool flag1 = true;
      if( leftout_set_size != 0 ) {
        // forward search

        //cerr << "SetWise: left out clusters size " << leftout_set_size << endl;
        cerr << "forward search !\n";

        for (int i=0; i< leftout_set_size; i++) {
          // add one cluster each time

          trial_set = trial_set_orig;
          trial_set.push_back(leftout_set.at(i));

          sort(trial_set.begin(), trial_set.end());

          ECIcluster_tmp.GetECICluster(trial_set);

          // get correlation functions of ECI cluster
          for (uint i=0; i<str_list.size(); i++) {
            str_list.at(i).SetECICluster(ECIcluster_tmp);
            str_list.at(i).GetECICorrelation();
            //str_list.at(i).PrintOutCorrelation();
          }

          score = ACECrossValidation(str_leftout_set_list,
              str_list, ECIcluster_tmp);

          //cerr << "score " << score << endl;
          //cerr << "score opt " << score_opt << endl;

          // keep only the shortest optimal ECI cluster with
          // the same score
          if(score < score_opt && score_opt > _SCORE_ZERO) {
            score_opt = score;
            trial_set_opt = trial_set;
            flag1 = false;
          } else if(score_opt < _SCORE_ZERO) {
            break;
          }
        }
      }

      if( trial_set_size != 1 ) {
        // backward search

        cerr << "backward search !\n";

        for (int i=0; i< trial_set_size; i++) {
          // add one cluster each time

          trial_set = trial_set_orig;
          trial_set.erase(trial_set.begin()+i);

          //cerr << "backward \n";
          //for ( uint i = 0; i<trial_set.size(); i++) {
          //    cerr << "i " << i << " " << trial_set.at(i) << endl;
          //}

          //sort(trial_set.begin(), trial_set.end());

          ECIcluster_tmp.GetECICluster(trial_set);

          // get correlation functions of ECI cluster
          for (uint i=0; i<str_list.size(); i++) {
            str_list.at(i).SetECICluster(ECIcluster_tmp);
            str_list.at(i).GetECICorrelation();
            //str_list.at(i).PrintOutCorrelation();
          }

          score = ACECrossValidation(str_leftout_set_list,
              str_list, ECIcluster_tmp);

          //cerr << "score " << score << endl;
          //cerr << "score opt " << score_opt << endl;

          if(score <= score_opt) {
            score_opt = score;
            trial_set_opt = trial_set;
            flag1 = false;
          }
        }
      }

      flag = flag1;
      trial_set = trial_set_opt;
      trial_set_size = trial_set.size();

      cerr << "SetWise: trial set size " << trial_set_size << endl;
      cerr << "current best cluster size " << trial_set_opt.size() << endl;
      cerr << "score_opt " << score_opt << endl;
      cerr << "current best cluster \n";
      for (uint i=0; i<trial_set_opt.size(); i++) {
        cerr << trial_set_opt.at(i) << " ";
      }
      cerr << endl;
    }

    cerr << "current optimal cluster: "
      << trial_set_opt.size() << endl;
    cerr << "score opt " << score_opt << endl;
    for (uint i=0; i<trial_set_opt.size(); i++) {
      cerr << trial_set_opt.at(i) << " ";
    }
    cerr << endl;

    if(score_opt < score_opt_global) {
      // global minimum
      score_opt_global = score_opt;
      trial_set_opt_global = trial_set_opt;
    }

    // if the score is exactly zero, stop the loop
    if( score_opt_global < _SCORE_ZERO) {
      break;
    }
  }

  cerr << "global score opt " << score_opt_global << endl;
  cerr << "global optimal cluster size " << trial_set_opt_global.size() << endl;
  for (uint i=0; i<trial_set_opt_global.size(); i++) {
    cerr << trial_set_opt_global.at(i) << " ";
  }
  cerr << endl;

  // store ECI clusters and averaged ECI
  ECIcluster.GetECICluster(trial_set_opt_global);

  vector<double> ECI_opt;
  ECI_opt = ACEGetECIAverage( str_leftout_set_list,
      str_list,  ECIcluster);

  ECIcluster.SetECI(ECI_opt);
  ECIcluster.SetChiSQ(score_opt_global);

  return trial_set_opt_global;
}


vector<int> GeneticAlgorithm( vector<cestructure> & str_list,
    ceECIcluster & ECIcluster) {
  // genetic algorithm using rank selection and elitism
  // randomly choose the crossover position

  // set up the cross validating sets
  vector< vector<int> > str_leftout_set_list;
  int total_str = str_list.size();
  int left_out_num = _LEFTOUT_NUM;

  str_leftout_set_list = ACESetUpCrossValidationSetsRandom(
      left_out_num, total_str);


  // set up paramters
  int total_cluster = ECIcluster.ECICluster().size(); // size of gene pool
  int population = _POPULATION;
  int crossover_position;
  vector<int> trial_set_opt_global;
  int total_gene_num = aurostd::min(_GENE_NUM, total_cluster);
  int best=0;
  vector<int> trial_set;

  bool chromosome[population][total_cluster]; // use array representive
  vector<fittness> fit_score;

  ofstream myfile;
  myfile.open(_SCOREFILE.c_str());

  double score;

  srand(time(0));

  for (int i=0; i<population; i++) {
    //for (int i1=0; i1<2; i1++) {
    //    // always include (1 0 1 1), (2 1 1 2)
    //    chromosome[i][i1] = true;
    //}
    //for (int i1=2; i1<total_cluster; i1++) {
    //    chromosome[i][i1] = false;
    //}
    for (int i1=0; i1<total_cluster; i1++) {
      chromosome[i][i1] = false;
    }
  }

  cerr << "Initial total gene_num "
    << total_gene_num << endl;

  int evolution_count = 0;

  // first generation
  for (int i=0; i<population; i++) {
    bool flag1 = false;
    while (!flag1) {
      int gene_num = aurostd::max(6, rand()%total_gene_num);


      for (int i1=0; i1<gene_num; i1++) {
        int gene_index = rand()%total_cluster;
        bool flag;
        flag = true;
        for (int i2=0; i2<i1; i2++) {
          if( chromosome[i][i1] == gene_index ) {
            //duplicate
            flag = false;
            break;
          }
        }
        if(flag) {
          chromosome[i][gene_index] = true;
        }
      }
      if( i == 0 ) {
        flag1 = true;
      } else {
        bool flag2 = true;

        // flag2 = ture : no duplication
        //       = false: duplication

        for (int i1=0; i1<i; i1++) { // no duplication of chromosomes
          flag2 = true;
          for (int i2=0; i2<total_cluster; i2++) {
            flag2 = flag2 && (!(chromosome[i][i2]^chromosome[i1][i2]));
          }
          flag2 = !flag2;

          if(!flag2) {
            break;
          }
        }

        if(flag2) {
          flag1 = true;
        } else {
          flag1 = false;
        }

      }
    }
  }

  //cerr << "Initial population\n";
  //for (int i=0; i<population; i++) {
  //    for (int i1=0; i1<total_cluster; i1++) {
  //        cerr << chromosome[i][i1];
  //    }
  //    cerr << endl;
  //}


  // assign fittness score

  fit_score.clear();
  for (int i=0; i<population; i++) {
    trial_set = GetTrialSet(chromosome[i], total_cluster);

    //for (int l=0; l<total_cluster; l++) {
    //    cerr << chromosome[i][l];
    //}
    //cerr << endl;
    //for (uint l=0; l< trial_set.size(); l++) {
    //    cerr << trial_set.at(l) << " ";
    //}

    ceECIcluster ECIcluster_tmp;
    ECIcluster_tmp = ECIcluster;

    ECIcluster_tmp.GetECICluster(trial_set);

    // get correlation functions of ECI cluster
    for (uint i1=0; i1<str_list.size(); i1++) {
      str_list.at(i1).SetECICluster(ECIcluster_tmp);
      str_list.at(i1).GetECICorrelation();
      //str_list.at(i1).PrintOutCorrelation();
    }


    score = ACECrossValidation(str_leftout_set_list, str_list, ECIcluster_tmp);

    fittness fit_score_tmp;
    fit_score_tmp.score = score;
    fit_score_tmp.rank = 0.;
    fit_score_tmp.index = i;
    fit_score_tmp.size = trial_set.size();

    fit_score.push_back(fit_score_tmp);
  }


  sort(fit_score.begin(), fit_score.end(), &comparison_fitscore);

  // output scores to a file
  cerr.precision(6);
  for (uint i=0; i<fit_score.size(); i++) {
    myfile << setw(8)
      << evolution_count
      << setw(12)
      <<fit_score.at(i).score
      << endl ;
  }

  int no_replacement_num = 0;
  int same_best_chromosome_num = 0;

  while ( no_replacement_num < _NO_REPLACEMENT_NUM
      && evolution_count < _EVOLUTION_NUM
      && same_best_chromosome_num < _SAME_BEST_CHROMOSOME_NUM) {
    // reproduction

    //////////////////////////////////////////////////////////////
    // selection
    // linear-ranking
    // http://www.geatbx.com/docu/algindex-02.html#P240_15119
    // assign rank-based fittness
    // fit_new = 2-SP+2*(SP-1)*(Pos-1)/(N-1)
    // SP: selective pressure in [1.0, 2.0]
    //
    //////////////////////////////////////////////////////////////


    for (uint i=0; i<fit_score.size(); i++) {
      fit_score.at(i).rank = 2.0-_SP + 2.0*(_SP-1.0)*double(i)/(population-1);
    }

    // convert fit_score to selective probability
    double sum;
    sum = 0.0;
    for (uint i=0; i<fit_score.size(); i++) {
      sum += fit_score.at(i).rank;
    }
    for (uint i=0; i<fit_score.size(); i++) {
      fit_score.at(i).rank /= sum;
    }

    //cerr << "score list \n";
    //cerr.precision(6);
    //for (uint i=0; i<fit_score.size(); i++) {
    //    cerr << fit_score.at(i).score << " "
    //        << fit_score.at(i).rank << " "
    //        << fit_score.at(i).index << endl;
    //}

    // build the reproduction pool
    int reproduction_num = int(_REPRODUCTION_RATE*population) +
      int(_REPRODUCTION_RATE*population) % 2; // always be even

    vector<fittness> fit_score_offspring;
    bool parent[reproduction_num][total_cluster];
    vector<int> reproduction_pool;

    for (int i=0; i<reproduction_num; i++) {

      bool flag = false;
      int count=0;
      while (!flag) {
        double r = double(rand()) / double(RAND_MAX);
        double sum_tmp=0.0;

        for (int i1=fit_score.size()-1; i1>=0; i1--) {
          sum_tmp += fit_score.at(i1).rank;
          if( r < sum_tmp ) {
            count = i1;
            break;
          }
        }
        flag = true;

        for (uint i=0; i<reproduction_pool.size(); i++) {
          if( count == reproduction_pool.at(i) ) {
            flag = false;
            break;
          }
        }

      }

      reproduction_pool.push_back(count);

      for (int i1=0; i1<total_cluster; i1++) {
        int index = fit_score.at(count).index;
        parent[i][i1] = chromosome[index][i1];
      }

    }


    //// print out the parents
    //cerr << "parents\n";
    //for (int i1=0; i1<reproduction_num; i1++) {
    //    for (int i2=0; i2<total_cluster; i2++) {
    //        cerr << parent[i1][i2];
    //    }
    //    cerr << endl;
    //}

    //// print out reproduction pool
    //for (uint i1=0; i1<reproduction_pool.size(); i1++) {
    //    cerr << reproduction_pool.at(i1) << " ";
    //}
    //cerr << endl;

    // reproduction of offsprings
    // crossover
    // randomly choose the crossover position

    //crossover_position = max(2, rand()%total_cluster);
    crossover_position = (total_cluster-2)/2;
    for (int i1=2; i1<reproduction_num/2; i1++) {
      int r;
      r = rand() % reproduction_pool.size();
      int index1, index2;
      index1 = r;

      reproduction_pool.erase(reproduction_pool.begin()+r);

      r = rand() % reproduction_pool.size();
      index2 = r;


      for (int i2=crossover_position-1; i2<total_cluster; i2++) {
        bool flag = parent[index1][i2];
        parent[index1][i2] = parent[index2][i2];
        parent[index2][i2] = flag;
      }


      reproduction_pool.erase(reproduction_pool.begin()+r);

    }


    //// after mating
    //cerr << "cross over position " <<crossover_position<<endl;
    //cerr << "after mating\n";
    //for (int i1=0; i1<reproduction_num; i1++) {
    //    for (int i2=0; i2<total_cluster; i2++) {
    //        cerr << parent[i1][i2];
    //    }
    //    cerr << endl;
    //}


    // mutation
    for (int j=0; j<reproduction_num; j++) {

      double rd;

      for (int j1=0; j1<total_cluster; j1++) {

        rd = double(rand()) / double(RAND_MAX);
        if( rd < _MUTATION_PROPABILITY ) {
          parent[j][j1] = !(parent[j][j1]);
        }

      }

    }


    //// after mutation
    //cerr << "after mutation \n";
    //for (int i1=0; i1<reproduction_num; i1++) {
    //    for (int i2=0; i2<total_cluster; i2++) {
    //        cerr << parent[i1][i2];
    //    }
    //    cerr << endl;
    //}


    // check score

    for (int j1=0; j1<reproduction_num; j1++) {

      trial_set = GetTrialSet(parent[j1], total_cluster);

      ceECIcluster ECIcluster_tmp;
      ECIcluster_tmp = ECIcluster;

      ECIcluster_tmp.GetECICluster(trial_set);

      // get correlation functions of ECI cluster
      for (uint i1=0; i1<str_list.size(); i1++) {
        str_list.at(i1).SetECICluster(ECIcluster_tmp);
        str_list.at(i1).GetECICorrelation();
        //str_list.at(i1).PrintOutCorrelation();
      }

      score = ACECrossValidation(str_leftout_set_list, str_list, ECIcluster_tmp);

      fittness fit_score_tmp;
      fit_score_tmp.score = score;
      fit_score_tmp.rank = 0.0;
      fit_score_tmp.index = j1;
      fit_score_tmp.size = trial_set.size();

      fit_score_offspring.push_back(fit_score_tmp);
    }


    sort(fit_score_offspring.begin(), fit_score_offspring.end(), &comparison_fitscore);

    //cerr << "offspring score\n";
    //for (uint i=0; i<fit_score_offspring.size(); i++) {
    //    cerr << fit_score_offspring.at(i).score << " "
    //        << fit_score_offspring.at(i).index << endl;
    //}

    //cerr << "offspring population\n";
    //for (int i=0; i<reproduction_num; i++) {
    //    for (int i1=0; i1<total_cluster; i1++) {
    //        cerr << parent[i][i1];
    //    }
    //    cerr << endl;
    //}


    //cerr << "before replacement population score\n";
    //for (uint i=0; i<fit_score.size(); i++) {
    //    cerr << fit_score.at(i).score << " "
    //        << fit_score.at(i).index << endl;
    //}
    //cerr << "Initial population\n";
    //for (int i=0; i<population; i++) {
    //    for (int i1=0; i1<total_cluster; i1++) {
    //        cerr << chromosome[i][i1];
    //    }
    //    cerr << endl;
    //}

    // replacement
    for (int i=fit_score_offspring.size()-1; i>=0; i--) {
      double score_offspring = fit_score_offspring.at(i).score;
      int index_offspring = fit_score_offspring.at(i).index;
      int size_offspring = fit_score_offspring.at(i).size;

      //////////////////////////////////////////////////////////////
      // make sure that the offspring is not in the population
      // flag = true: not found in the population
      //      = false: found in the population
      //////////////////////////////////////////////////////////////
      bool flag=true;
      for (int i1=0; i1<population; i1++) {
        flag = true;
        for (int i2=0; i2<total_cluster; i2++) {
          flag = flag && (!(parent[index_offspring][i2]^chromosome[i1][i2]));
        }
        flag = !flag;
        if(!flag) {
          //cerr << "offspring " << index_offspring
          //    << " is the same as parent " << i1 << endl;
          break;
        }
      }

      if(flag) {

        for (uint i1=0; i1<fit_score.size(); i1++) {

          double score_parent = fit_score.at(i1).score;
          int index_parent = fit_score.at(i1).index;

          if( score_offspring < score_parent ) {
            //cerr << "offspring " << index_offspring
            //    << " replace parent "  << index_parent << endl;
            for (int j=0; j<total_cluster; j++) {
              chromosome[index_parent][j] = parent[index_offspring][j];
            }
            fit_score.at(i1).score = score_offspring;
            fit_score.at(i1).size = size_offspring;
            break;
          }
        }

        no_replacement_num = 0;
      } else {
        ++no_replacement_num;
      }

    }

    sort(fit_score.begin(), fit_score.end(), &comparison_fitscore);

    //cerr << "population score\n";
    //for (uint i=0; i<fit_score.size(); i++) {
    //    cerr << fit_score.at(i).score << " "
    //        << fit_score.at(i).index << endl;
    //}

    //cerr << "next generation\n";
    //for (int i=0; i<population; i++) {
    //    for (int i1=0; i1<total_cluster; i1++) {
    //        cerr << chromosome[i][i1];
    //    }
    //    cerr << endl;
    //}

    int best_old = best;
    best = fit_score.back().index;
    if(best_old == best) {
      ++same_best_chromosome_num;
    } else {
      same_best_chromosome_num = 0;
    }
    cerr << "evolution " << evolution_count
      << " current best " << best
      << " score " << fit_score.back().score
      << " size " << fit_score.back().size
      << endl;
    cerr << "current best chromosome \n";
    int count_pa = 0;
    for (int i1=0; i1<total_cluster; i1++) {
      cerr << chromosome[best][i1];
      if( chromosome[best][i1] ) {
        ++count_pa;
      }
    }
    cerr << " size " << count_pa << endl;
    cerr << "no replacement number " << no_replacement_num << endl;
    cerr << "same best chromosome number " << same_best_chromosome_num << endl;

    // output scores to a file
    ++evolution_count;
    cerr.precision(6);
    for (uint i=0; i<fit_score.size(); i++) {
      myfile << setw(8)
        << evolution_count
        << setw(12)
        <<fit_score.at(i).score
        << endl ;
    }

    if( fit_score.back().score < _SCORE_ZERO ) {
      break;
    }
  }

  best = fit_score.back().index;
  cerr << "best " << best << " score " << fit_score.back().score << endl;

  trial_set_opt_global = GetTrialSet(chromosome[best], total_cluster);

  cerr << "best chromosome \n";
  for (int i1=0; i1<total_cluster; i1++) {
    cerr << chromosome[best][i1];
  }
  cerr << endl;
  cerr << "best ECI cluster: size " << trial_set_opt_global.size() << endl;
  for (uint i1=0; i1<trial_set_opt_global.size(); i1++) {
    cerr << trial_set_opt_global.at(i1) << " ";
  }
  cerr << endl;

  myfile.close();

  // store ECI clusters and averaged ECI
  ECIcluster.GetECICluster(trial_set_opt_global);
  ECIcluster.PrintOutECICluster();
  vector<double> ECI_opt;

  cerr << "Average ECI cluster size " << ECIcluster.ECICluster().size() << endl;

  ECI_opt = ACEGetECIAverage( str_leftout_set_list,
      str_list,  ECIcluster);
  ECIcluster.SetECI(ECI_opt);
  ECIcluster.SetChiSQ(fit_score.back().score);

  return trial_set_opt_global;
}

vector<int> GetTrialSet(bool *chromosome, int col) {
  vector<int> trial_set;
  for (int i=0; i<col; i++) {
    if( chromosome[i] ) {
      trial_set.push_back(i);
    }
  }

  return trial_set;
}

void PrintChromosome(ostream & os, bool *chromosome[], int col, int row) {

  for (int i=0; i<row; i++) {
    for (int i1=0; i1<col; i1++) {
      os << chromosome[i][i1];
    }
    os << endl;
  }
}

bool comparison_fitscore(const fittness & score1, const fittness & score2) {
  // true: case 1. score1.score > score2.score
  //       case 2. if score1.score == score2.score or
  //               score1.score and score2.score are effectively zero
  //               score1.size > score2.size
  // false: otherwise

  bool flag;
  flag = (score1.score > _SCORE_ZERO) || (score2.score > _SCORE_ZERO);

  if( (score1.score != score2.score) && flag) {
    return (score1.score > score2.score);
  } else {
    if( score1.size == score2.size ) {
      return (score1.score > score2.score);
    } else {
      return (score1.size >= score2.size);
    }
  }

}

// ***************************************************************************
// Calculate ECI's
// ***************************************************************************

void ACEGetECI(vector<cestructure> & str_list, ceECIcluster & cluster1) {
  // calculate ECI from the input structure list
  // first get matrice X, Y, and Y_sigma

  int nrow = str_list.size();
  int ncol = str_list.at(0).ECICorrelation().size();
  xvector<double> y_vec(1, nrow);
  xvector<double> y_sig(1, nrow);
  xmatrix<double> x_mat(1, 1, nrow, ncol);

  for(int i=0; i<nrow; i++) { // ground state energy
    y_vec[i+1] = str_list.at(i).EnergyIn();
  }


  // here set standard deviations of y's to 1
  // can be supplied as input
  for(int i=1; i<= nrow; i++) {
    y_sig[i] = _SIGMA;
  }

  //cerr << "ncol " << ncol << endl;
  //cerr << "nrow " << nrow << endl;


  for(int i=0; i< nrow;i++) { // different structures
    for (int j=0; j< ncol; j++) { // different clusters

      x_mat[i+1][j+1] = str_list.at(i).ECICorrelation().at(j)
        * str_list.at(i).ECIEquivalentNum().at(j);
    }
  }


  aurostd::cematrix ECI_matrix(x_mat);
  ECI_matrix.LeastSquare(y_vec, y_sig);

  vector<double> ECI = ECI_matrix.AVec();
  double chisq = ECI_matrix.ChiSQ();

  cluster1.SetECI(ECI);
  cluster1.SetChiSQ(chisq);

  //double chisq=cluster1.GetECI(x_mat, y_vec, y_sig);

  for (uint i=0; i<str_list.size(); i++) {
    str_list.at(i).SetECI(ECI);
    str_list.at(i).SetChiSQ(chisq);
  }

  //cout << "#ECI \n";
  //for (uint i=0; i<ECI.size(); i++) {
  //    string name;

  //    if( i==0 ) {
  //        cout << setw(10)
  //            << "0 0 0 0";
  //    } else {
  //        int index = cluster1.ECICluster().at(i-1);
  //        name = cluster1.GetClusterNameByIndex(index);
  //        cout << setw(10)
  //            << name;
  //    }
  //    cout.setf(ios_base::fixed, ios_base::floatfield);
  //    cout.precision(8);
  //    cout << setw(16)
  //        << ECI.at(i)
  //        << endl;
  //}

  //cout << "# chisq " << cluster1.ChiSQ() << endl;
}


void ACEGetECI(vector<cestructure> & str_list, vector<int> fit_list, ceECIcluster & cluster1) {
  // calculate ECI from the input structure list
  // only those with index in fit_list are used as fit set
  // those left out are used to asset the predictive power
  // first get matrice X, Y, and Y_sigma
  // ECI, chisq has been stored in cluster1

  int nrow = fit_list.size();
  int ncol = str_list.at(0).ECICorrelation().size();
  xvector<double> y_vec(1, nrow);
  xvector<double> y_sig(1, nrow);
  xmatrix<double> x_mat(1, 1, nrow, ncol);

  for(int i=0; i<nrow; i++) { // ground state energy
    int index = fit_list.at(i);
    y_vec[i+1] = str_list.at(index).EnergyIn();
  }


  // here set standard deviations of y's to 1
  // can be supplied as input
  for(int i=1; i<= nrow; i++) {
    y_sig[i] = _SIGMA;
  }

  //cerr << "ncol " << ncol << endl;
  //cerr << "nrow " << nrow << endl;


  for(int i=0; i< nrow;i++) { // different structures
    int index = fit_list.at(i);
    for (int j=0; j< ncol; j++) { // different clusters

      x_mat[i+1][j+1] = str_list.at(index).ECICorrelation().at(j)
        * str_list.at(index).ECIEquivalentNum().at(j);
    }
  }


  aurostd::cematrix ECI_matrix(x_mat);
  ECI_matrix.LeastSquare(y_vec, y_sig);

  vector<double> ECI = ECI_matrix.AVec();
  double chisq = ECI_matrix.ChiSQ();

  //for (uint i=0; i<ECI.size(); i++) {
  //    ECI.at(i) = ECI.at(i)/1000.0;
  //}

  cluster1.SetECI(ECI);
  cluster1.SetChiSQ(chisq);

  //double chisq=cluster1.GetECI(x_mat, y_vec, y_sig);


  for (uint i=0; i<str_list.size(); i++) {
    str_list.at(i).SetECI(ECI);
    str_list.at(i).SetChiSQ(chisq);
  }

  //cout << "#ECI \n";
  //for (uint i=0; i<ECI.size(); i++) {
  //    string name;

  //    if( i==0 ) {
  //        cout << setw(10)
  //            << "0 0 0 0";
  //    } else {
  //        int index = cluster1.ECICluster().at(i-1);
  //        name = cluster1.GetClusterNameByIndex(index);
  //        cout << setw(10)
  //            << name;
  //    }
  //    cout.setf(ios_base::fixed, ios_base::floatfield);
  //    cout.precision(8);
  //    cout << setw(16)
  //        << ECI.at(i)
  //        << endl;
  //}

  //cout << "# chisq " << cluster1.ChiSQ() << endl;
}

vector<int> LeftOutSet(int total_size, vector<int> kept_set) {
  // get a complementary set from a base set {0 ... total_size-1}

  vector<int> comp_set;

  for (int i=0; i < total_size; i++) {
    comp_set.push_back(i);
  }

  for (uint i=0; i<kept_set.size(); i++) {
    for (uint j=0; j<comp_set.size(); j++) {
      if( kept_set.at(i) == comp_set.at(j) ) {
        comp_set.erase(comp_set.begin()+j);
        break;
      }
    }
  }

  return comp_set;

}

vector<double> ACEGetECIAverage(vector< vector<int> > str_leftout_set_list,
    vector<cestructure> & str_list, ceECIcluster & ECIcluster) {

  // calculate the averaged ECI

  vector<int> fit_list, leftout_list;
  int total_size = str_list.size();
  vector<double> ECI_ave, ECI_new;

  //int total_str = str_list.size();
  int total_config;
  total_config = str_leftout_set_list.size();

  // get correlation functions of ECI cluster
  for (uint i=0; i<str_list.size(); i++) {
    str_list.at(i).SetECICluster(ECIcluster);
    str_list.at(i).GetECICorrelation();
  }

  for (int i=0; i<total_config; i++) {

    leftout_list = str_leftout_set_list.at(i);

    fit_list=LeftOutSet(total_size, leftout_list);

    ACEGetECI(str_list, fit_list, ECIcluster);

    ECI_new = ECIcluster.ECIValue();
    if(i == 0) {
      ECI_ave = ECI_new;
    } else {
      for (uint i=0; i<ECI_ave.size(); i++) {
        ECI_ave.at(i) += ECI_new.at(i);
      }
    }

  }

  vector<double> ECI_tmp, ECIcluster_tmp;

  for (uint i=0; i<ECI_ave.size(); i++) {
    ECI_ave.at(i) /= total_config;
  }

  return ECI_ave;
}

// ****************************************************
// Cross Validation
// ****************************************************

double ACECrossValidation(vector< vector<int> > str_leftout_set_list,
    vector<cestructure> & str_list, ceECIcluster & ECIcluster) {
  // Select the optimal ECIs by cross validation scores
  // leftout_list sets are given by left_out_set_list

  vector<int> fit_list, leftout_list;
  int total_size = str_list.size();
  double total_score;

  ofstream myfile;
  myfile.open(_LOCVFILE.c_str());
  myfile << "Left some out cross-validation\n";


  //int total_str = str_list.size();
  double pre_error;
  int total_config;
  total_config = str_leftout_set_list.size();

  //cerr << "In Validation function: str left out list size "
  //    << str_leftout_set_list.size() << endl;
  //for (uint l=0; l<str_leftout_set_list.size(); l++) {
  //    for (uint l1=0; l1<str_leftout_set_list.at(l).size(); l1++) {
  //        cerr << str_leftout_set_list.at(l).at(l1) << " ";
  //    }
  //    cerr << endl;
  //}
  //cerr << endl;

  total_score = 0.0;
  for (int i=0; i<total_config; i++) {

    //leftout_list.clear();
    //for(int i=0; i<left_out_num; i++) {
    //    leftout_list.push_back(rand()%total_size);
    //}
    leftout_list = str_leftout_set_list.at(i);

    fit_list=LeftOutSet(total_size, leftout_list);

    myfile << "Left out ";
    //cerr << "Left out ";
    for (uint k=0; k<leftout_list.size(); k++) {
      int nr = leftout_list.at(k);
      myfile <<  str_list.at(nr).Name() << " ";

      //cerr << str_list.at(nr).Name() << " ";
    }
    myfile << endl;
    //cerr << endl;

    ACEGetECI(str_list, fit_list, ECIcluster);

    myfile << "fit chisq = " << ECIcluster.ChiSQ() << endl;

    for (uint i=0; i<str_list.size(); i++) {
      str_list.at(i).GetEnergy();
    }

    ECIcluster.PrintECI(myfile);

    for (uint k=0; k<leftout_list.size(); k++) {
      int nr = leftout_list.at(k);
      pre_error = (str_list.at(nr).Energy() - str_list.at(nr).EnergyIn())
        / str_list.at(nr).EnergyIn();

      myfile << "structure "
        << str_list.at(nr).Name()
        << " predict error  = "
        <<  pre_error
        << endl;
    }

    double score = 0.0;
    for (uint k=0; k<str_list.size(); k++) {
      double de = (str_list.at(k).Energy() - str_list.at(k).EnergyIn());
      score += de*de;
    }

    score = sqrt(score)/double(str_list.size());
    total_score += score;

    myfile << "CV score " << score << endl;
    //cerr << "CV score " << score << endl;

    for (uint i=0; i<str_list.size(); i++) {
      str_list.at(i).PrintOutComparison(myfile);
    }

    myfile << endl;
  }

  total_score /= double(total_config);
  myfile << "total score " << total_score << endl;

  myfile.close();

  return total_score;
}

vector< vector<int> >  ACESetUpCrossValidationSetsRandom(int left_out_num,
    int total_str) {
  // Select the optimal ECIs by cross validation scores
  // leftout_list is chosen by random number

  vector< vector<int> > str_leftout_set_list;

  vector<int>  leftout_list;

  //int total_config = _TRIAL_RATIO*CombinationNr(left_out_num, total_str);
  int combination_nr = CombinationNr(left_out_num, total_str);

  int total_config;
  if( total_str > 10 && left_out_num > 2) {
    if( combination_nr > _TRIAL_NUM ) {
      total_config = _TRIAL_NUM;
    } else {
      total_config = combination_nr;
    }
  } else {
    total_config =  aurostd::min(_TRIAL_NUM, total_str);
  }

  srand(time(0));

  //double total_score;
  //total_score = 0.0;
  for (int i=0; i<total_config; i++) {

    leftout_list.clear();
    for(int i=0; i<left_out_num; i++) {
      leftout_list.push_back(rand()%total_str);
    }

    str_leftout_set_list.push_back(leftout_list);

  }

  cerr << "SetUp: str left out list size " << str_leftout_set_list.size() << endl;
  for (uint l=0; l<str_leftout_set_list.size(); l++) {
    for (uint l1=0; l1<str_leftout_set_list.at(l).size(); l1++) {
      cerr << str_leftout_set_list.at(l).at(l1) << " ";
    }
    cerr << endl;
  }
  cerr << endl;

  return str_leftout_set_list;
}


// ********************************************************
// Superlattices
// ********************************************************

//////////////////////////////////////////////////////////////
// Get SL with unit cell atom number between cell_nr_min
// and cell_nr_max
//////////////////////////////////////////////////////////////

void ACESLProperties(ostream & os, string & structure_type,
    int cell_nr_min, int cell_nr_max, const ceallclusters & allcluster1,
    ceECIcluster & ecicluster1) {
  ceSL sl1;
  sl1 = ceSL(structure_type);
  vector<int> atom_config;
  atom_config.push_back(0); // the first atom is always type A

  ofstream mySLfile;
  ofstream mySLfilename;
  ofstream mySLfileresult;
  mySLfile.open(_SLFILE.c_str());

  mySLfilename.open(_SLFILENAME.c_str());
  mySLfilename.precision(8);
  mySLfilename.setf(ios_base::fixed, ios_base::floatfield);

  ofstream mySLfilecor;
  mySLfilecor.open(_SLFILECOR.c_str());

  os << endl;

  int vol;
  xmatrix<int> trans_SL(1, 1, 3, 3); // to get the lattice vector of SL

  xmatrix<double> lattice_old(1, 1, 3, 3);
  vector< vector<double> > lattice_old_list;

  xstructure str_base;
  str_base=aflowlib::PrototypeLibraries(cerr, aurostd::utype2string(_FCC_BEGIN),"",LIBRARY_MODE_HTQC); // no parameters
  str_base.FixLattices();
  str_base.CalculateSymmetryPointGroup();

  xvector<double> b1(1, 3), b2(1, 3), b3(1, 3);

  for (int nr=cell_nr_min; nr<cell_nr_max+1; nr++) {

    int nr1, nr2, nr3;
    double r1, r2;
    r1 = std::pow((double) 4.0*double(nr)*double(nr), (double) 1.0/3.0);
    r2 = r1/3.0 + sqrt( r1*r1/9 + 4.0*double(nr)*double(nr) );

    nr1 = int(sqrt(r1));
    nr2 = int(sqrt(r2));
    nr3 = nr;

    for (int l11=-nr1; l11<nr1+1; l11++) {
      for (int l12=-nr1; l12 < nr1+1; l12++) {
        for (int l13=-nr1; l13 < nr1+1; l13++) {

          b1[1] = l11*str_base.lattice[1][1]
            +l12*str_base.lattice[2][1]
            +l13*str_base.lattice[3][1];
          b1[2] = l11*str_base.lattice[1][2]
            +l12*str_base.lattice[2][2]
            +l13*str_base.lattice[3][2];
          b1[3] = l11*str_base.lattice[1][3]
            +l12*str_base.lattice[2][3]
            +l13*str_base.lattice[3][3];


          for (int l21=-nr2; l21 < nr2+1; l21++) {
            for (int l22=-nr2; l22 < nr2+1; l22++) {
              for (int l23=-nr2; l23 < nr2+1; l23++) {

                b2[1] = l21*str_base.lattice[1][1]
                  +l22*str_base.lattice[2][1]
                  +l23*str_base.lattice[3][1];
                b2[2] = l21*str_base.lattice[1][2]
                  +l22*str_base.lattice[2][2]
                  +l23*str_base.lattice[3][2];
                b2[3] = l21*str_base.lattice[1][3]
                  +l22*str_base.lattice[2][3]
                  +l23*str_base.lattice[3][3];


                if( l11 == l21 && l12 == l22 && l13 == l23 ) {
                  continue;
                }

                //if( scalar_product(b1,b2) < 0 ) {
                //    continue;
                //}

                if( aurostd::modulus(b1) > aurostd::modulus(b2) ) {
                  continue;
                }

                double leng_max = aurostd::max(aurostd::modulus(b1), aurostd::modulus(b2));
                if( (aurostd::modulus(b2-b1) < leng_max) || (aurostd::modulus(b2+b1) < leng_max ) ) {
                  continue;
                }


                for (int l31=-nr3; l31 < nr3+1; l31++) {
                  for (int l32=-nr3; l32 < nr3+1; l32++) {
                    for (int l33=-nr3; l33 < nr3+1; l33++) {


                      b3[1] = l31*str_base.lattice[1][1]
                        +l32*str_base.lattice[2][1]
                        +l33*str_base.lattice[3][1];
                      b3[2] = l31*str_base.lattice[1][2]
                        +l32*str_base.lattice[2][2]
                        +l33*str_base.lattice[3][2];
                      b3[3] = l31*str_base.lattice[1][3]
                        +l32*str_base.lattice[2][3]
                        +l33*str_base.lattice[3][3];

                      //if( scalar_product(b1, b3) < 0
                      //        || scalar_product(b2, b3) < 0) {
                      //    continue;
                      //}

                      if( aurostd::modulus(b2) > aurostd::modulus(b3) ) {
                        continue;
                      }

                      double leng_max;
                      leng_max = aurostd::max(aurostd::modulus(b1), aurostd::modulus(b3));
                      if( (aurostd::modulus(b3-b1) < leng_max) || (aurostd::modulus(b3+b1) < leng_max ) ) {
                        continue;
                      }
                      leng_max = aurostd::max(aurostd::modulus(b2), aurostd::modulus(b3));
                      if( (aurostd::modulus(b3-b2) < leng_max) || (aurostd::modulus(b3+b2) < leng_max ) ) {
                        continue;
                      }


                      vol = (l11*l22 - l12*l21)*l33 + (l12*l23-l13*l22)*l31
                        + (l13*l21 - l11*l23)*l32;
                      //if( abs(vol) != nr ) {
                      //    continue;
                      //}

                      if( vol != nr ) {
                        // keep only vol > 0
                        continue;
                      }

                      trans_SL[1][1] = l11;
                      trans_SL[1][2] = l12;
                      trans_SL[1][3] = l13;
                      trans_SL[2][1] = l21;
                      trans_SL[2][2] = l22;
                      trans_SL[2][3] = l23;
                      trans_SL[3][1] = l31;
                      trans_SL[3][2] = l32;
                      trans_SL[3][3] = l33;

                      //cerr << "trans_SL \n";
                      //cerr << trans_SL << endl;

                      sl1.GetStructure(trans_SL, nr);

                      if( sl1.AtomNr() != nr ) {
                        continue;
                      }

                      // check whether sl1 has been found or not
                      // follow the criterion given by Ferreira, Wei, and Zunger
                      // Int. J. High Performance Comp. Appl. 5, 34 (1991)

                      //sl1.Structure().FixLattices();

                      //cerr << sl1.Name() << endl;
                      //cerr << sl1.Structure().lattice << endl << endl;

                      xmatrix<double> grot(1, 1, 3, 3); // point-group*k-lattice
                      bool flag_eq = false;

                      for (uint k1=0; k1 < lattice_old_list.size(); k1++) {
                        lattice_old[1][1] = lattice_old_list.at(k1).at(0);
                        lattice_old[1][2] = lattice_old_list.at(k1).at(1);
                        lattice_old[1][3] = lattice_old_list.at(k1).at(2);
                        lattice_old[2][1] = lattice_old_list.at(k1).at(3);
                        lattice_old[2][2] = lattice_old_list.at(k1).at(4);
                        lattice_old[2][3] = lattice_old_list.at(k1).at(5);
                        lattice_old[3][1] = lattice_old_list.at(k1).at(6);
                        lattice_old[3][2] = lattice_old_list.at(k1).at(7);
                        lattice_old[3][3] = lattice_old_list.at(k1).at(8);

                        for (uint l1=0; l1< str_base.pgroup.size(); l1++) {
                          grot = sl1.Structure().klattice*str_base.pgroup.at(l1).Uc/(2.0*_pi);
                          double a1, a2, a3; // lattice *dot* reciprocal lattice

                          bool flag1=true;
                          for (int l2=1; l2 < 4; l2++) {
                            a1 = abs(grot[l2][1]*lattice_old[1][1]
                                + grot[l2][2]*lattice_old[1][2]
                                + grot[l2][3]*lattice_old[1][3]) + _EQUAL_DOUBLE*0.1;
                            a2 = abs(grot[l2][1]*lattice_old[2][1]
                                + grot[l2][2]*lattice_old[2][2]
                                + grot[l2][3]*lattice_old[2][3]) + _EQUAL_DOUBLE*0.1;
                            a3 = abs(grot[l2][1]*lattice_old[3][1]
                                + grot[l2][2]*lattice_old[3][2]
                                + grot[l2][3]*lattice_old[3][3]) + _EQUAL_DOUBLE*0.1;

                            // check if a1, a2, a2 are integer

                            //cerr << endl << grot << endl << endl;
                            //cerr << lattice_old << endl << endl;
                            //cerr << a1 << " " << a2 << " " << a3 << endl;

                            flag1 = flag1 && ((abs(a1 - int(a1)) < _EQUAL_DOUBLE)
                                && (abs(a2 - int(a2)) < _EQUAL_DOUBLE)
                                && (abs(a3 - int(a3)) < _EQUAL_DOUBLE));

                          }

                          flag_eq = flag1;

                          if(flag_eq) {
                            break;
                          }

                        }

                        if(flag_eq) {
                          break;
                        }

                      }

                      if( flag_eq ) {
                        continue;
                      } else {
                        // store the lattice
                        vector<double> a_vec;
                        for (int l1=1; l1<4; l1++) {
                          for (int l2=1; l2<4; l2++) {
                            a_vec.push_back(sl1.Structure().lattice[l1][l2]);
                          }
                        }
                        lattice_old_list.push_back(a_vec);
                      }

                      for (int i=1; i<sl1.AtomNr(); i++) {
                        int total_size = CombinationNr(i,sl1.AtomNr()-1);

                        for (int j=0; j<total_size; j++) {

                          int index = j;
                          if( index ==0 ) {
                            atom_config = AllCombination41(i, total_size, index+1);
                          } else {
                            atom_config = AllCombination42(i, total_size, atom_config);
                          }

                          //for (uint k=0; k<atom_config.size(); k++) {
                          //    cerr << atom_config.at(k) << " ";
                          //}
                          //cerr << endl;

                          for(uint i1=0; i1<atom_config.size(); i1++) {
                            atom_config.at(i1) += 1;
                          }

                          sl1.SetConfig(atom_config);
                          sl1.StructureToName();

                          // restore to the original atom_config for getting correct combinations
                          for(uint i1=0; i1<atom_config.size(); i1++) {
                            atom_config.at(i1) -= 1;
                          }

                          if( !sl1.IsPrimitiveCell() ) {
                            continue;
                          }

                          sl1.PrintStructure(mySLfile);

                          sl1.SetStrCluster(allcluster1);
                          sl1.SetAllECICluster(ecicluster1);

                          sl1.GetAllECICorrelation();
                          sl1.WriteFile(mySLfilecor);

                          mySLfilename << setw(12)
                            << sl1.Name()
                            << setw(12)
                            << sl1.StoichB()
                            << endl;

                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }


  mySLfile.close();
  mySLfilename.close();
  mySLfilecor.close();

}

void ACESLProperties_Readin_Corfile(istream & os, string & structure_type,
    const ceallclusters & allcluster1,
    ceECIcluster & ecicluster1) {

  // read correlation from the file
  // not as fast as calculating correlation when
  // the number of superllatice is very large slow
  // due to the slow I/O

  ceSL sl1;
  sl1 = ceSL(structure_type);
  string SL_name;

  ofstream mySLfileresult;
  mySLfileresult.open(_SLFILERESULT.c_str());

  ifstream corfilein;
  string corfilename = _SLFILECOR;
  corfilein.open(corfilename.c_str());

  clock_t start, end;
  double diff;
  //////////////////////////////////////////////////////////////
  start = clock();
  //////////////////////////////////////////////////////////////

  sl1.SetStrCluster(allcluster1);
  sl1.SetAllECICluster(ecicluster1);
  sl1.SetECICluster(ecicluster1);

  int pos = 0;
  double stoich_b;
  while ( os >> SL_name >> stoich_b) {

    //cerr << "SL_name " << SL_name << endl;

    // To save time, sl1 does not contain all
    // necessary information like lattice vectors
    // this form of sl1 can be only used here

    sl1.SetName(SL_name);
    sl1.SetStoichB(stoich_b);

    int pos_new;
    pos_new = sl1.GetAllECICorrelation(corfilein, pos);
    if( pos_new == _NOCORRELATON_FOUND ) {
      cerr << SL_name
        << " : correlations are not found. Ignore this superlattice.\n";
      continue;
    }
    pos = pos_new;

    vector<double> ECI = ecicluster1.ECIValue();
    sl1.SetECI(ECI);

    sl1.GetECICorrelation();
    //sl1.PrintOutCorrelation();


    sl1.GetEnergy();

    mySLfileresult << sl1 << endl;

  }
  //////////////////////////////////////////////////////////////
  end = clock();
  diff = (end - start)/ CLOCKS_PER_SEC;
  cerr << "Calculate SL Time: " << diff << " s\n";
  //////////////////////////////////////////////////////////////


  corfilein.close();
  mySLfileresult.close();

}

void ACESLPropertiesSQS_Readin_Corfile(istream & os, string & structure_type,
    const ceallclusters & allcluster1,
    ceECIcluster & ecicluster1, int site_num, int NNNum,
    int min_SLcell_nr, int max_SLcell_nr,
    vector<_ceatom> atom_species) {

  // get SQS's with different concentration and number of atom inside a unit cell

  ceSL sl1;
  sl1 = ceSL(structure_type);
  string SL_name;

  ifstream corfilein;
  string corfilename = _SLFILECOR;
  corfilein.open(corfilename.c_str());

  ofstream sqsfileout;
  string sqsfilename = _SLFILESQS;
  sqsfileout.open(sqsfilename.c_str());

  vector< vector<SLSQS> > SQS_list;
  vector<SLSQS> SQS_list_tmp;
  SQS_list_tmp.resize(2);
  SLSQS sqs_tmp;

  clock_t start, end;
  double diff;
  //////////////////////////////////////////////////////////////
  start = clock();
  //////////////////////////////////////////////////////////////

  sl1.SetStrCluster(allcluster1);
  sl1.SetAllECICluster(ecicluster1);
  sl1.SetECICluster(ecicluster1);

  const double _MAX_SQS_WEIGHT = 1.0e8;
  int pos = 0;
  double stoich_b;
  bool flag;
  while ( os >> SL_name >> stoich_b) {

    sl1.SetUp(SL_name);
    sl1.SetStoichB(stoich_b);

    if( sl1.CellNr() < min_SLcell_nr ) {
      continue;
    } else if( sl1.CellNr() > max_SLcell_nr) {
      break;
    }

    int pos_new;
    pos_new = sl1.GetAllECICorrelation(corfilein, pos);
    pos = pos_new;

    SQS_list_tmp.clear();
    SQS_list_tmp.resize(2);

    double SQSweight;
    SQSweight = sl1.IsSQS(site_num, NNNum);

    flag = true;

    for (uint i=0; i<SQS_list.size() && flag; i++) {

      if( abs(sl1.StoichB() - SQS_list.at(i).at(0).stoich_b)
          < _EQUAL_DOUBLE
          && sl1.CellNr() == SQS_list.at(i).at(0).cell_nr ) {
        // same concentration and cell number

        for (uint j=0; j<SQS_list.at(i).size(); j++) {
          // store new SL if its weigh is smaller

          if(  SQSweight <= SQS_list.at(i).at(j).weight ) {
            // smaller weight

            if( is_equal_structure(sl1, SQS_list.at(i).at(j).SQS) ) {
              // store only structure with different
              // correlations, that is, not equivalent
              // symmetrically
              flag = false;
              break;
            }

            for (uint k=SQS_list.at(i).size()-1; k>j; k--) {
              SQS_list.at(i).at(k) = SQS_list.at(i).at(k-1);
            }

            sqs_tmp.weight = SQSweight;
            sqs_tmp.SQS = sl1;
            sqs_tmp.stoich_b = sl1.StoichB();
            sqs_tmp.cell_nr = sl1.CellNr();
            SQS_list.at(i).at(j) = sqs_tmp;
            flag = false;
            break;
          }
        }
        flag = false;
      }

    }

    if( flag ) {
      // add new SQS item with different concentration and/or cell number
      sqs_tmp.weight = SQSweight;
      sqs_tmp.SQS = sl1;
      sqs_tmp.stoich_b = sl1.StoichB();
      sqs_tmp.cell_nr = sl1.CellNr();
      SQS_list_tmp.at(0) = sqs_tmp;
      SQS_list_tmp.at(1) = sqs_tmp;
      SQS_list_tmp.at(1).weight = _MAX_SQS_WEIGHT;
      SQS_list.push_back(SQS_list_tmp);
    }
  }

  //////////////////////////////////////////////////////////////
  end = clock();
  diff = (end - start)/ CLOCKS_PER_SEC;
  cerr << "Calculate SL Time: " << diff << " s\n";
  //////////////////////////////////////////////////////////////

  corfilein.close();

  if( SQS_list.size() == 0 ) {
    cerr << "No superlattice with number of atoms in a unit cell \n"
      << "between "
      << min_SLcell_nr
      << " and "
      << max_SLcell_nr
      << ". Please run the command aflow --sl "
      << min_SLcell_nr
      << " "
      << max_SLcell_nr
      << " to get those superlattices."
      << endl;
  } else {
    sort(SQS_list.begin(), SQS_list.end(), &comparison_SLSQS);

    vector<int> cell_nr_list;
    for (uint i = 0; i<SQS_list.size(); i++) {
      cell_nr_list.push_back(SQS_list.at(i).at(0).cell_nr);
    }
    sort(cell_nr_list.begin(), cell_nr_list.end());
    int cell_nr_min = cell_nr_list.at(0);
    int cell_nr_max = cell_nr_list.back();

    if( cell_nr_min != min_SLcell_nr || cell_nr_max != max_SLcell_nr ) {
      sqsfileout  << "Obtain the two best condidates of SQS's of superlattice \n"
        << "with number of atoms in a unit cell between "
        << cell_nr_min
        << " and "
        << cell_nr_max
        << endl;

      if( cell_nr_min > min_SLcell_nr ) {
        sqsfileout << "Please run the command aflow --sl "
          << min_SLcell_nr
          << " "
          << cell_nr_min
          << " to get those superlattices."
          << endl;
      }
      if( cell_nr_max < max_SLcell_nr ) {
        sqsfileout << "Please run the command aflow --sl "
          << cell_nr_max
          << " "
          << max_SLcell_nr
          << " to get those superlattices."
          << endl;
      }
    }

    // output SQS_list
    for (uint i=0; i<SQS_list.size(); i++) {
      for (uint j=0; j<SQS_list.at(i).size(); j++) {
        sqsfileout << "------------------------------------------------------------\n";
        sqsfileout << setw(20)
          << "weight "
          << setw(20)
          << SQS_list.at(i).at(j).weight << endl;
        SQS_list.at(i).at(j).SQS.OutputSQS(sqsfileout);
        SQS_list.at(i).at(j).SQS.PrintStructure(sqsfileout, atom_species);
        sqsfileout << "------------------------------------------------------------\n";
        sqsfileout << endl;
      }
    }
    sqsfileout.close();
  }

}

bool comparison_SLSQS(const vector<SLSQS> sqs1_list, const vector<SLSQS> sqs2_list) {
  bool flag = false;
  if( sqs1_list.at(0).stoich_b < sqs2_list.at(0).stoich_b ) {
    flag = true;
  } else if( sqs1_list.at(0).stoich_b == sqs2_list.at(0).stoich_b ) {
    if( sqs1_list.at(0).cell_nr < sqs2_list.at(0).cell_nr) {
      flag = true;
    } else {
      flag = false;
    }
  } else {
    flag = false;
  }

  return flag;
}

bool is_equal_structure(ceSL sl1, ceSL sl2) {
  // two SL are equivalent symmetrically or not

  if( sl1.SQScompareClusterList().size() !=
      sl2.SQScompareClusterList().size() ) {
    return false;
  } else {
    for (uint i=0; i<sl1.SQScompareClusterList().size(); i++) {
      if( sl1.SQScompareClusterList().at(i).pair_name !=
          sl2.SQScompareClusterList().at(i).pair_name ||
          sl1.SQScompareClusterList().at(i).correlation !=
          sl2.SQScompareClusterList().at(i).correlation ) {
        return false;
      }
    }
    return true;
  }
}

void ACESLProperties(istream & os, string & structure_type,
    const ceallclusters & allcluster1,
    ceECIcluster & ecicluster1) {

  // calculate correlation

  ceSL sl1;
  sl1 = ceSL(structure_type);
  string SL_name;

  ofstream mySLfileresult;
  mySLfileresult.open(_SLFILERESULT.c_str());

  //ofstream mySLfile;
  ////mySLfile.open("SL_tmp.dat");
  //mySLfile.open(_SLFILE.c_str());

  //ifstream corfilein;
  //string corfilename = _SLFILECOR;
  //corfilein.open(corfilename.c_str());

  vector<double> ECI;

  double stoich_b;
  while ( os >> SL_name >> stoich_b) {

    //cerr << "SL_name " << SL_name << endl;

    sl1.SetUp(SL_name);
    //sl1.PrintStructure(mySLfile);

    sl1.SetStrCluster(allcluster1);
    sl1.SetAllECICluster(ecicluster1);

    sl1.GetAllECICorrelation();
    //sl1.PrintOutCorrelation();

    ECI = ecicluster1.ECIValue();
    sl1.SetECI(ECI);

    sl1.GetECICorrelation();
    sl1.GetEnergy();

    mySLfileresult << sl1 << endl;

  }

  //mySLfile.close();

  //corfilein.close();
  mySLfileresult.close();

}

void ACESLCorrelations(istream & os, string & structure_type,
    const ceallclusters & allcluster1,
    ceECIcluster & ecicluster1) {
  // calculate the correlations of all candidate clusters
  // store in slcor.dat to accelerate the computation
  // ecicluster1 must be a ceECIcluster object without
  // assigning eci clusters
  ceSL sl1;
  sl1 = ceSL(structure_type);
  string SL_name;
  double stoich_b;

  ofstream mySLfilecor;
  mySLfilecor.open(_SLFILECOR.c_str());

  //ofstream mySLfileresult;
  ////mySLfileresult.open(_SLFILERESULT.c_str());

  vector<double> ECI;
  while ( os >> SL_name >> stoich_b) {

    sl1.SetUp(SL_name);

    sl1.SetStrCluster(allcluster1);
    sl1.SetAllECICluster(ecicluster1);

    sl1.GetAllECICorrelation();

    // writefile must be used with ECI_cluster with
    // all eci cluster set, not with the fitted ECI clusters
    sl1.WriteFile(mySLfilecor);

    //sl1.PrintOutCorrelation();

    sl1.SetECICluster(ecicluster1);

    ECI = ecicluster1.ECIValue();
    cerr << ECI.size() << endl;
    sl1.SetECI(ECI);

    //sl1.GetECICorrelation();
    //sl1.GetEnergy();

    //mySLfileresult << sl1 << endl;

  }

  mySLfilecor.close();
  //mySLfileresult.close();

}


void GenerateAflowInputFile(string structure_type,
    string AlloyName,
    string SL_name, vector<_ceatom> atom_species, bool mpi_flag) {

  ofstream input_file;
  string aflowin=_AFLOWIN_;
  input_file.open(aflowin.c_str());

  input_file << AFLOWIN_SEPARATION_LINE << endl;
  input_file << "[AFLOW]                                                                               " << endl;
  input_file << "[AFLOW]                     .o.        .o88o. oooo                                    " << endl;
  input_file << "[AFLOW]                    .888.       888 `` `888                                    " << endl;
  input_file << "[AFLOW]                   .8'888.     o888oo   888   .ooooo.  oooo oooo    ooo        " << endl;
  input_file << "[AFLOW]                  .8' `888.     888     888  d88' `88b  `88. `88.  .8'         " << endl;
  input_file << "[AFLOW]                 .88ooo8888.    888     888  888   888   `88..]88..8'          " << endl;
  input_file << "[AFLOW]                .8'     `888.   888     888  888   888    `888'`888'           " << endl;
  input_file << "[AFLOW]               o88o     o8888o o888o   o888o `Y8bod8P'     `8'  `8'  .in       " << endl;
  input_file << "[AFLOW]                                                                               " << endl;
  input_file << AFLOWIN_SEPARATION_LINE << endl;
  input_file << "[AFLOW] * Stefano Curtarolo - (AFLOW V" << string(AFLOW_VERSION) << ") " << endl;
  input_file << "[AFLOW] * Dane Morgan - Wahyu Setyawan - Gus Hart - Michal Jahnatek - Shidong Wang - Ohad Levy " << endl;
  input_file << AFLOWIN_SEPARATION_LINE << endl;
  input_file << "[AFLOW] Aflow automatically generated (aflow_sfunc.cpp) " << endl;
  input_file << AFLOWIN_SEPARATION_LINE << endl;
  input_file << AFLOWIN_SEPARATION_LINE << endl;
  input_file << "[AFLOW]SYSTEM=" << AlloyName << endl;
  input_file << AFLOWIN_SEPARATION_LINE << endl;
  input_file << "[AFLOW] input file for aflow " << endl;
  input_file << "[AFLOW_MODE=VASP] " << endl;
  input_file << AFLOWIN_SEPARATION_LINE << endl;
  input_file << "[AFLOW_MODE_ZIP=bzip2] " << endl;
  if( mpi_flag ) {
    input_file << "#[AFLOW_MODE_BINARY=vasp46s] " << endl;
  } else {
    input_file << "[AFLOW_MODE_BINARY=vasp46s] " << endl;
  }
  input_file << AFLOWIN_SEPARATION_LINE << endl;
  input_file << AFLOWIN_SEPARATION_LINE << endl;
  if( mpi_flag ) {
    input_file << "[AFLOW_MODE_MPI]" << endl;
  } else {
    input_file << "#[AFLOW_MODE_MPI]" << endl;
  }
  input_file << "#[AFLOW_MODE_MPI_MODE]NCPUS=MAX " << endl;
  input_file << "[AFLOW_MODE_MPI_MODE]NCPUS=4" << endl;
#ifdef MPI_LAM
  input_file << "[AFLOW_MODE_MPI_MODE]START=\"lamboot\" " << endl;
  input_file << "[AFLOW_MODE_MPI_MODE]STOP=\"lamhalt\" " << endl;
#endif
  input_file << "[AFLOW_MODE_MPI_MODE]COMMAND =\"mpirun -np\" " << endl;
  input_file << "[AFLOW_MODE_MPI_MODE]AUTOTUNE " << endl;
  input_file << "[AFLOW_MODE_MPI_MODE]BINARY=\"mpivasp46s\" " << endl;
  input_file << AFLOWIN_SEPARATION_LINE << endl;
  input_file << "[AFLOW_SYMMETRY]CALC " << endl;
  input_file << "#[AFLOW_SYMMETRY]SGROUP_WRITE " << endl;
  input_file << "#[AFLOW_SYMMETRY]SGROUP_RADIUS=7.77 " << endl;
  input_file << AFLOWIN_SEPARATION_LINE << endl;
  //DX20210122 [OBSOLETE] input_file << "#[AFLOW_NEIGHBOURS]CALC " << endl;
  //DX20210122 [OBSOLETE] input_file << "[AFLOW_NEIGHBOURS]RADIUS=7.7 " << endl;
  //DX20210122 [OBSOLETE] input_file << "[AFLOW_NEIGHBOURS]DRADIUS=0.1 " << endl;
  //DX20210122 [OBSOLETE] input_file << AFLOWIN_SEPARATION_LINE << endl;
  input_file << "[AFLOW] # Phonons calculations are not implemented yet. They will be included in aflow4." << endl;
  input_file << "#[AFLOW_PHONONS]CALC " << endl;
  input_file << AFLOWIN_SEPARATION_LINE << endl;
  input_file << "#[VASP_RUN]GENERATE                                   // GENERATE STATIC RELAX=N RELAX_STATIC=N STATIC_BANDS RELAX_STATIC_BANDS=N " << endl;
  input_file << "#[VASP_RUN_RELAX=2] " << endl;
  input_file << "[VASP_RUN_RELAX_STATIC=2] " << endl;
  input_file << "#[VASP_RUN_RELAX_STATIC_BANDS=2] " << endl;
  input_file << "#[VASP_FORCE_OPTION]NEGLECT_NOMIX " << endl;
  input_file << "#[VASP_FORCE_OPTION]CHGCAR=OFF                         // ON | OFF (default ON)" << endl;
  input_file << "[VASP_FORCE_OPTION]CHGCAR=ON                         // ON | OFF (default ON)" << endl;
  input_file << "[VASP_FORCE_OPTION]KPOINTS=KEEPK " << endl;
  input_file << "#[VASP_FORCE_OPTION]KPOINTS=EVEN                      // EVEN | ODD (default none)" << endl;
  input_file << "#[VASP_FORCE_OPTION]KPOINTS=KSHIFT_GAMMA_EVEN         // EVEN | ODD (default none)" << endl;
  input_file << "#[VASP_FORCE_OPTION]KPOINTS=KSCHEME_MONKHORST_PACK    // MONKHORST_PACK | GAMMA (manual)" << endl;
  input_file << "#[VASP_FORCE_OPTION]KPOINTS=GAMMA " << endl;
  input_file << "#[VASP_FORCE_OPTION]KPOINTS=IBZKPT " << endl;
  input_file << "[VASP_FORCE_OPTION]SYM=ON                             // ON | OFF  (default ON)" << endl;
  input_file << "[VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=potpaw_PBE   // pot_LDA | pot_GGA | potpaw_LDA | potpaw_GGA | potpaw_PBE | potpaw_LDA_KIN | potpaw_PBE_KIN  " << endl;
  input_file << "[VASP_FORCE_OPTION]NBANDS                             // Estimate Bands (better than VASP)" << endl;
  input_file << "[VASP_FORCE_OPTION]SPIN=ON,REMOVE_RELAX_1             // (ON | OFF  (default ON)), REMOVE_RELAX_1 | _2" << endl;
  input_file << "#[VASP_FORCE_OPTION]AUTO_MAGMOM=ON                    // ON | OFF (default OFF)" << endl;
  input_file << "[VASP_FORCE_OPTION]RELAX_MODE=ENERGY                  // ENERGY | FORCES | ENERGY_FORCES | FORCES_ENERGY (default ENERGY) " << endl;
  input_file << "[VASP_FORCE_OPTION]PREC=MEDIUM                        // (LOW | MEDIUM | NORMAL | HIGH | ACCURATE), PRESERVED (default=MEDIUM)" << endl;
  input_file << "[VASP_FORCE_OPTION]ALGO=FAST                          // (NORMAL | VERYFAST | FAST | ALL | DAMPED), PRESERVED (default=NORMAL)" << endl;
  input_file << "[VASP_FORCE_OPTION]RELAX " << endl;
  input_file << "#[VASP_FORCE_OPTION]NOTUNE " << endl;
  input_file << "[VASP_FORCE_OPTION]TYPE=METAL                     // METAL | INSULATOR | SEMICONDUCTOR | DEFAULT (default DEFAULT) " << endl;
  input_file << "#[VASP_FORCE_OPTION]CONVERT_UNIT_CELL=something       // SPRIM, SCONV, NIGGLI, MINK, INCELL, INCOMPACT, WS, CART, FRAC" << endl;
  input_file << "#[VASP_FORCE_OPTION]VOLUME+=10.0 " << endl;
  input_file << "#[VASP_FORCE_OPTION]VOLUME*=1.05 " << endl;
  input_file << AFLOWIN_SEPARATION_LINE << endl;
  input_file << AFLOWIN_SEPARATION_LINE << endl;
  input_file << "[VASP_INCAR_MODE_EXPLICIT]START " << endl;
  input_file << "SYSTEM=" << AlloyName << endl;
  input_file << "NELM = 120" << endl;
  input_file << "NELMIN=2" << endl;
  input_file << "LPLANE=.TRUE." << endl;
  input_file << "LREAL=.FALSE." << endl;
  input_file << "LSCALU=.FALSE." << endl;
  input_file << "PSTRESS=000                                           # for hand modification" << endl;
  input_file << "#NBANDS=XX                                            # for hand modification" << endl;
  input_file << "#IALGO=48                                             # for hand modification" << endl;
  input_file << "[VASP_INCAR_MODE_EXPLICIT]STOP " << endl;
  input_file << AFLOWIN_SEPARATION_LINE << endl;
  input_file << "[VASP_KPOINTS_MODE_IMPLICIT] " << endl;
  input_file << "[VASP_KPOINTS_FILE]KSCHEME=M " << endl;
  input_file << "[VASP_KPOINTS_FILE]KPPRA=2000" << endl;
  input_file << "[VASP_KPOINTS_FILE]STATIC_KSCHEME=M " << endl;
  input_file << "[VASP_KPOINTS_FILE]STATIC_KPPRA=2000" << endl;
  input_file << "[VASP_KPOINTS_FILE]BANDS_LATTICE=BCC" << endl;
  input_file << "[VASP_KPOINTS_FILE]BANDS_GRID=16" << endl;
  input_file << AFLOWIN_SEPARATION_LINE << endl;
  input_file << AFLOWIN_SEPARATION_LINE << endl;
  input_file << "[VASP_POSCAR_MODE_EXPLICIT]START " << endl;

  // output SL structure
  ceSL sl1;
  sl1 = ceSL(structure_type);

  sl1.SetUp(SL_name);
  sl1.PrintStructure(input_file, atom_species);

  // output the left part

  input_file << "[VASP_POSCAR_MODE_EXPLICIT]STOP " << endl;
  input_file << AFLOWIN_SEPARATION_LINE << endl;
  input_file << "[VASP_POTCAR_MODE_IMPLICIT] " << endl;
  input_file << "[VASP_POTCAR_FILE]" << atom_species.at(0).name << endl;
  input_file << "[VASP_POTCAR_FILE]" << atom_species.at(1).name << endl;
  input_file << "[AFLOW] potpaw_PBE: "
    << atom_species.at(0).name
    << " " << atom_species.at(1).name << endl;
  input_file << "[VASP_FORCE_OPTION]LDAU_PARAMETERS="
    << atom_species.at(0).name
    <<","
    << atom_species.at(1).name
    << ";1,0;4.0,4.0;0.0,0.0"
    << endl;
  input_file.close();

}

// ***************************************************************************
//  Miscellaneous functions
// ***************************************************************************

void CheckAllInputFileExistence(string structure_type) {
  // if no cluster data exists, ask user to
  // obtain it first
  string filename;
  filename = structure_type+_FILENAME;

  ifstream myfile;
  myfile.open(filename.c_str());

  bool flag_cluster_cal;
  if( !myfile.is_open()) {

    cerr << "The cluster data file is missing\n"
      << "Please use the command\n"
      << "    aflow --cluster structure_type minimun_site_num maximum_site_num minimum_nearest_neighbor maximum_nearest_neighbor \n"
      << "to generate it\n";

    myfile.close();
    flag_cluster_cal = true;

  } else {
    flag_cluster_cal = false;
    myfile.close();
  }

  // if no superlattice input files exist, ask user to
  // obtain them first
  ifstream mySLfilein;
  bool flag_SL_cal = true;
  mySLfilein.open(_SLFILENAME.c_str());
  if(! mySLfilein.is_open() ) {
    flag_SL_cal = true;
    mySLfilein.close();
  } else {
    mySLfilein.close();

    ifstream mySLcorfilein;
    mySLcorfilein.open(_SLFILECOR.c_str());

    if( ! mySLcorfilein.is_open() ) {

      flag_SL_cal = true;
      mySLcorfilein.close();

    } else {
      flag_SL_cal = false;
      mySLcorfilein.close();

    }

  }
  if( flag_SL_cal ) {
    cerr << "One or more input files for superlattice calculations are missing.\n"
      << "Please run the command \n"
      << "    aflow --sl structure_type mininum_atom_in_unit_cell maximum_atom_in_unit_cell \n"
      << "to obtained them\n";
  }

  if( flag_cluster_cal || flag_SL_cal ) {
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"no input",_INPUT_MISSING_); //CO20200624
  } else {
    cerr << "All input files are found. Calculation begins.\n"
      << endl;
  }

}


// [OBSOLETE]// carrier effective mass calculations
// [OBSOLETE]bool _OLD_pflow::EffectiveMass(string directory_name,
// [OBSOLETE]			   string suffix,
// [OBSOLETE]			   double& band_gap,
// [OBSOLETE]			   vector<double>& mass_electron_dos,
// [OBSOLETE]			   vector<double>& mass_hole_dos,
// [OBSOLETE]			   vector<double>& mass_electron_conduction,
// [OBSOLETE]			   vector<double>& mass_hole_conduction,
// [OBSOLETE]			   const bool& osswrite,
// [OBSOLETE]			   ostream& oss) {
// [OBSOLETE]  // aflow --effective_mass directory_name
// [OBSOLETE]  // OR
// [OBSOLETE]  // aflow --em directory_name
// [OBSOLETE]  // input: directory containing vasp EIGENVAL, DOSCAR, OUTCAR
// [OBSOLETE]  // *suffix.bz2 generated by aflow will be checked by default
// [OBSOLETE]  // The band gap, number of spins, and some other data are read out of the DOSCAR by getEfermiBandGap()
// [OBSOLETE]  // The extrema (valleys) of the valence and conduction bands are gathered by getFitDataForMass()
// [OBSOLETE]  // The curvature of the valley is solved for by fitting an ellipse to the data by fitToEllipsisEquation()
// [OBSOLETE]  // From the curvature, the effective mass tensor is calculated, and the eigenvalues are found.  
// [OBSOLETE]  // These eigenvalues are used as m_1, m_2, and m_3 for each valley.
// [OBSOLETE]  // The number of equivalent valleys in the Brillouin zone is found from the point group of the 
// [OBSOLETE]  // reciprocal lattice.
// [OBSOLETE]  // The DoS effective masses are computed by averaging over the crystal using these data 
// [OBSOLETE]  // in the following equation (copy & paste to LaTeX)
// [OBSOLETE]  // m^*_\text{carrier} = \left( \frac{\sum_i^\text{Nvalleys} M_i^2 m_{i1} m_{i2} m_{i3}}{ \text{Nvalleys} } \right)^{\frac{1}{3}}
// [OBSOLETE]  // and the conductivity effective masses are calculated by (copy & paste to LaTeX)
// [OBSOLETE]  // m^*_\text{carrier} = \frac{3 \left( \sum_i^{ \text{Nvalley}} M_i \right)}{\sum_i^{ \text{Nvalley}} M_i \left( \frac{1}{m_{i1}} + \frac{1}{m_{i2}} + \frac{1}{m_{i3}} \right)}

// [OBSOLETE]  // The resulting effective masses are written to the vectors mass_electron_dos, mass_hole_dos, mass_electron_conduction, and mass_hole_conduction passed in the function arguments.
// [OBSOLETE]  // Finally, osswrite=TRUE will produce a summary of the individual valley statistics to be written to oss

// [OBSOLETE]  if(osswrite) oss.setf(std::ios::fixed);
// [OBSOLETE]  if(osswrite) oss.precision(2);          // limit the number of decimals in the output
// [OBSOLETE]
// [OBSOLETE]  bool needed_files_flag = true;
// [OBSOLETE]
// [OBSOLETE]  // empty and initialize the containers that were passed into pflow::EffectiveMass()
// [OBSOLETE]  mass_electron_dos.clear();
// [OBSOLETE]  mass_electron_dos.push_back(0);
// [OBSOLETE]  mass_electron_dos.push_back(0); // to make space for UP and DOWN
// [OBSOLETE]  mass_hole_dos.clear();
// [OBSOLETE]  mass_hole_dos.push_back(0);
// [OBSOLETE]  mass_hole_dos.push_back(0); // to make space for UP and DOWN
// [OBSOLETE]  mass_electron_conduction.clear();
// [OBSOLETE]  mass_electron_conduction.push_back(0);
// [OBSOLETE]  mass_electron_conduction.push_back(0); // to make space for UP and DOWN
// [OBSOLETE]  mass_hole_conduction.clear();
// [OBSOLETE]  mass_hole_conduction.push_back(0);
// [OBSOLETE]  mass_hole_conduction.push_back(0); // to make space for UP and DOWN

// [OBSOLETE]  vector<int> valley_electron(2), valley_hole(2); // used for the DoS effective mass denominator
// [OBSOLETE]  string eigenval_file, outcar_file, doscar_file, poscar_file;

// [OBSOLETE]  // prepare a list of the input files
// [OBSOLETE]  vector<string> files;
// [OBSOLETE]  files.push_back("EIGENVAL");
// [OBSOLETE]  files.push_back("DOSCAR");
// [OBSOLETE]  files.push_back("POSCAR");
// [OBSOLETE]  vector<string> destination;
// [OBSOLETE]  destination.push_back(aurostd::TmpFileCreate("EIGENVAL"));
// [OBSOLETE]  destination.push_back(aurostd::TmpFileCreate("DOSCAR"));
// [OBSOLETE]  destination.push_back(aurostd::TmpFileCreate("POSCAR"));

// [OBSOLETE]  // find and decompress (if necessary) the input files; return FALSE if they aren't found
// [OBSOLETE]  for (uint i=0; i<files.size(); i++) {
// [OBSOLETE]    string source_suffix = directory_name+"/"+files.at(i)+suffix;
// [OBSOLETE]    string source_suffix_gz = directory_name+"/"+files.at(i)+suffix+".gz";
// [OBSOLETE]    string source_suffix_bz2 = directory_name+"/"+files.at(i)+suffix+".bz2";
// [OBSOLETE]    string source = directory_name+"/"+files.at(i);
// [OBSOLETE]    if(aurostd::FileExist(source_suffix)) aurostd::execute("cat "+source_suffix+" > "+destination.at(i));
// [OBSOLETE]    if(aurostd::FileExist(source_suffix_gz)) aurostd::execute("zcat "+source_suffix_gz+" > "+destination.at(i));
// [OBSOLETE]    if(aurostd::FileExist(source_suffix_bz2)) aurostd::execute("bzcat "+source_suffix_bz2+" > "+destination.at(i));
// [OBSOLETE]    if(aurostd::FileExist(source)) aurostd::execute("cat "+source+" > "+destination.at(i));
// [OBSOLETE]    if(!aurostd::FileExist(source_suffix) && !aurostd::FileExist(source_suffix_gz) && !aurostd::FileExist(source_suffix_bz2) && !aurostd::FileExist(source)) {
// [OBSOLETE]      needed_files_flag = false;
// [OBSOLETE]      break;
// [OBSOLETE]    }
// [OBSOLETE]  }
// [OBSOLETE]  if( needed_files_flag ) {
// [OBSOLETE]    eigenval_file = destination.at(0);
// [OBSOLETE]    doscar_file = destination.at(1);
// [OBSOLETE]    poscar_file = destination.at(2);
// [OBSOLETE]  } else {
// [OBSOLETE]    cerr << "one or more of EIGENVAL, DOSCAR, POSCAR, ";
// [OBSOLETE]    cerr << "EIGENVAL, DOSCAR, POSCAR is missing!\nAborting\n";
// [OBSOLETE]    return FALSE;
// [OBSOLETE]  }
// [OBSOLETE]  
// [OBSOLETE]  // prepare to read input files
// [OBSOLETE]  ifstream fin_eigenval, fin_outcar, fin_doscar;
// [OBSOLETE]  
// [OBSOLETE]  fin_eigenval.open(eigenval_file.c_str()); 
// [OBSOLETE]  int Nions, Ntmp, ispin;
// [OBSOLETE]  fin_eigenval >> Nions >> Nions >> Ntmp >> ispin;
// [OBSOLETE] 
// [OBSOLETE]  // **********************
// [OBSOLETE]  // get Fermi energy/valence band top/conduction band bottom from DOSCAR
// [OBSOLETE]  // **********************
// [OBSOLETE]  
// [OBSOLETE]  // prepare storage for the results of getEfermiBandGap()
// [OBSOLETE]  string compound_name;
// [OBSOLETE]  band_st band_information;
// [OBSOLETE]  
// [OBSOLETE]  if( !getEfermiBandGap(doscar_file, ispin, compound_name, band_information) ) {
// [OBSOLETE]    cerr << "getEfermiBandGap failed." << endl;
// [OBSOLETE]    cerr << XPID << "pflow::EffectiveMass exiting." << endl;
// [OBSOLETE]    return FALSE;
// [OBSOLETE]  }
// [OBSOLETE]
// [OBSOLETE]  if(band_information.band_gap > 0.0) {    // if system is a semiconductor
// [OBSOLETE]    // **********************
// [OBSOLETE]    // storage for the points around extremum to be fitted
// [OBSOLETE]    vector<vector<int> > number_of_valley_list;
// [OBSOLETE]    vector<vector<vector<kEn_st> > > fit_data_all;
// [OBSOLETE]
// [OBSOLETE]    // include additional bands within this range in eV
// [OBSOLETE]    double energy_range = 0.026;
// [OBSOLETE]
// [OBSOLETE]    // get the data for fitting; if there is a problem inside getFitDataForMass, then return false to say
// [OBSOLETE]    // we can't compute anything meaningful
// [OBSOLETE]    if( !getFitDataForMass(eigenval_file, poscar_file, band_information, number_of_valley_list, energy_range, fit_data_all) ) { 
// [OBSOLETE]      cerr << "getFitDataForMass failed." << endl;
// [OBSOLETE]      cerr << XPID << "pflow::EffectiveMass exiting." << endl;
// [OBSOLETE]      return FALSE;
// [OBSOLETE]    };
// [OBSOLETE]
// [OBSOLETE]    vector<int> index_of_extremas_zero;  
// [OBSOLETE]    int number_of_records = 0;
// [OBSOLETE]    vector<string> spin_name;
// [OBSOLETE]    spin_name.push_back("+");
// [OBSOLETE]    spin_name.push_back("-");
// [OBSOLETE]
// [OBSOLETE]    if(osswrite) oss << "Compound Name: " << compound_name << endl;
// [OBSOLETE]
// [OBSOLETE]    for (int spin_idx=0; spin_idx<ispin; spin_idx++) {
// [OBSOLETE]      vector<vector<kEn_st> > fit_data = fit_data_all.at(spin_idx);
// [OBSOLETE]      vector<vector<double> > mass_eff_list;
// [OBSOLETE]
// [OBSOLETE]      // fit the data for each valley of 1 spin
// [OBSOLETE]      if( !fitToEllipsisEquation(fit_data, spin_idx, mass_eff_list) ) {
// [OBSOLETE]	cerr << "fitToEllipsisEquation failed." << endl;
// [OBSOLETE]	cerr << XPID << "pflow::EffectiveMass exiting." << endl;
// [OBSOLETE]	return FALSE;
// [OBSOLETE]      }
// [OBSOLETE]      
// [OBSOLETE]      vector<double> electron_cond_mass_per_valley; // storage for the averaged effective mass of 1 valley 
// [OBSOLETE]      // for use in computing the total effective mass
// [OBSOLETE]      vector<int> electron_valley_multiplicity;     // store the rank of the constellation for a valley so
// [OBSOLETE]      // the index matches that of electron_cond_mass_per_valley
// [OBSOLETE]      vector<double> hole_cond_mass_per_valley;
// [OBSOLETE]
// [OBSOLETE]      vector<int> hole_valley_multiplicity;
// [OBSOLETE]
// [OBSOLETE]      mass_electron_dos.at(spin_idx) = 0.0; 
// [OBSOLETE]      mass_hole_dos.at(spin_idx) = 0.0;
// [OBSOLETE]      valley_electron.at(spin_idx) = 0; 
// [OBSOLETE]      valley_hole.at(spin_idx) = 0;
// [OBSOLETE]
// [OBSOLETE]      // loop over the valleys of 1 spin
// [OBSOLETE]      for (uint i=0; i<mass_eff_list.size(); i++) {
// [OBSOLETE]	// assign a carrier type to the valley
// [OBSOLETE]	string carrier_type;
// [OBSOLETE]	if(fit_data.at(i).at(0).band_type == 0) {
// [OBSOLETE]	  // valence band
// [OBSOLETE]	  carrier_type = "hole";
// [OBSOLETE]	} else if(fit_data.at(i).at(0).band_type == 1) {
// [OBSOLETE]	  // conduction band
// [OBSOLETE]	  carrier_type = "electron";
// [OBSOLETE]	}
// [OBSOLETE]
// [OBSOLETE]	number_of_records++;
// [OBSOLETE]
// [OBSOLETE]	//if(osswrite) oss << endl;
// [OBSOLETE]	if(osswrite) oss << number_of_records << ". Band index: " << fit_data.at(i).at(0).band_index << endl;
// [OBSOLETE]	if(osswrite) oss << "   Carrier type: " << carrier_type << endl;
// [OBSOLETE]	if(osswrite) oss << "   Carrier Spin: " << spin_name.at(spin_idx) << endl;
// [OBSOLETE]	if(osswrite) oss << "   Extrema Reciprocal Space Cartesian Coordinates: ";
// [OBSOLETE]	if(osswrite) oss << setprecision(3) << fit_data_all.at(0).at(i).at( 0  ).kpoint(1)  << ",  ";
// [OBSOLETE]	if(osswrite) oss << setprecision(3) << fit_data_all.at(0).at(i).at( 0  ).kpoint(2)  << ",  ";
// [OBSOLETE]	if(osswrite) oss << setprecision(3) << fit_data_all.at(0).at(i).at( 0  ).kpoint(3)  << ",  ";
// [OBSOLETE]	if(osswrite) oss << endl;
// [OBSOLETE]	if(osswrite) oss << "   Effective masses along principle axes (m0): ";
// [OBSOLETE]	
// [OBSOLETE]	double dos_mass = 1.0;
// [OBSOLETE]	double cond_mass = 0.0; // stores the valley's conduction effective mass
// [OBSOLETE]	if(osswrite) oss << "(" ;
// [OBSOLETE]	for (uint j=0; j<3; j++) {
// [OBSOLETE]	  dos_mass *= mass_eff_list.at(i).at(j);
// [OBSOLETE]	  double one = 1.0;
// [OBSOLETE]	  cond_mass += one / mass_eff_list.at(i).at(j); // add the reciprocal of 3 individual masses
// [OBSOLETE]	  if(osswrite) oss << setprecision(2) << mass_eff_list.at(i).at(j);
// [OBSOLETE]	  if(j<2) {
// [OBSOLETE]	    if(osswrite) oss << ",";
// [OBSOLETE]	  }
// [OBSOLETE]	}
// [OBSOLETE]	if(osswrite) oss << ")" << endl;
// [OBSOLETE]
// [OBSOLETE]	int number_of_valley = number_of_valley_list.at(spin_idx).at(i);
// [OBSOLETE]
// [OBSOLETE]	if(osswrite) oss << "   Number of equivalent valleys: " << number_of_valley << endl;
// [OBSOLETE]	
// [OBSOLETE]	dos_mass *= number_of_valley*number_of_valley;
// [OBSOLETE]
// [OBSOLETE]	if( carrier_type == "electron") {
// [OBSOLETE]	  mass_electron_dos.at(spin_idx) += dos_mass;
// [OBSOLETE]	  valley_electron.at(spin_idx)++;
// [OBSOLETE]	} else if( carrier_type == "hole") {
// [OBSOLETE]	  mass_hole_dos.at(spin_idx) += dos_mass;
// [OBSOLETE]	  valley_hole.at(spin_idx)++;
// [OBSOLETE]	}
// [OBSOLETE]
// [OBSOLETE]	dos_mass = std::pow((double) dos_mass, (double) 1.0/3.0);
// [OBSOLETE]	double temp = 1.0;
// [OBSOLETE]	temp /= cond_mass;
// [OBSOLETE]	cond_mass = temp; // flip sum of inverses to denominator
// [OBSOLETE]
// [OBSOLETE]	// store cond_mass for the valley and store multiplicity of the valley
// [OBSOLETE]	if( carrier_type == "electron") {
// [OBSOLETE]	  electron_cond_mass_per_valley.push_back(cond_mass);
// [OBSOLETE]	  electron_valley_multiplicity.push_back(number_of_valley);
// [OBSOLETE]	} else if( carrier_type == "hole" ) {
// [OBSOLETE]	  hole_cond_mass_per_valley.push_back(cond_mass);
// [OBSOLETE]	  hole_valley_multiplicity.push_back(number_of_valley);
// [OBSOLETE]	}
// [OBSOLETE]
// [OBSOLETE]	cond_mass *= 3.0; // final factor of 3 for single valley total
// [OBSOLETE]	if(osswrite) oss << "   DOS effective mass (m0): " << setprecision(2) << dos_mass << endl;
// [OBSOLETE]	if(osswrite) oss << "   Conduction mass (m0): " << setprecision(2) << cond_mass << endl;
// [OBSOLETE]      } // END: loop over mass_eff_list
// [OBSOLETE]      
// [OBSOLETE]      // DEBUG:  list the pieces needed to compute the carrier average conduction effective mass
// [OBSOLETE]      
// [OBSOLETE]	// for (vector<int>::size_type ix = 0; ix != electron_cond_mass_per_valley.size(); ++ix) {
// [OBSOLETE]	// cout << "m_e = "<< electron_cond_mass_per_valley.at(ix) << endl;
// [OBSOLETE]	// cout << "      multiplicity = " << electron_valley_multiplicity.at(ix) << endl;
// [OBSOLETE]	// }
// [OBSOLETE]      
// [OBSOLETE]      // Calculate the total electron conductivity effective mass
// [OBSOLETE]      double numerator = 0.0;
// [OBSOLETE]      double denominator = 0.0;
// [OBSOLETE]      for (vector<int>::size_type ix = 0; ix != electron_cond_mass_per_valley.size(); ++ix) {
// [OBSOLETE]	numerator += electron_valley_multiplicity.at(ix);
// [OBSOLETE]	denominator += electron_valley_multiplicity.at(ix)/electron_cond_mass_per_valley.at(ix);
// [OBSOLETE]      } 
// [OBSOLETE]      numerator *= 3;
// [OBSOLETE]      mass_electron_conduction.at(spin_idx) = numerator / denominator;
// [OBSOLETE]
// [OBSOLETE]      // Calculate the total hole conductivity effective mass
// [OBSOLETE]      numerator = 0.0;
// [OBSOLETE]      denominator = 0.0;
// [OBSOLETE]      for (vector<int>::size_type ix = 0; ix != hole_cond_mass_per_valley.size(); ++ix) {
// [OBSOLETE]	numerator += hole_valley_multiplicity.at(ix);
// [OBSOLETE]	denominator += hole_valley_multiplicity.at(ix)/hole_cond_mass_per_valley.at(ix);
// [OBSOLETE]      }
// [OBSOLETE]      numerator *= 3;
// [OBSOLETE]      mass_hole_conduction.at(spin_idx) = numerator / denominator;
// [OBSOLETE]    } // END:  loop over spins
// [OBSOLETE]    
// [OBSOLETE]    if(osswrite) oss << endl;
// [OBSOLETE]    if(osswrite) oss << "Egap: " << band_information.band_gap << " eV" << endl;
// [OBSOLETE]
// [OBSOLETE]    // write the DOS electron effective mass for all spins
// [OBSOLETE]    for (int spin_idx=0; spin_idx<ispin; spin_idx++) {
// [OBSOLETE]      if( spin_idx == 0 ) {
// [OBSOLETE]	if(osswrite) oss << "DOS electron effective mass (m0): {" ;
// [OBSOLETE]      } else if(spin_idx == 1) {
// [OBSOLETE]	if(osswrite) oss << ", ";
// [OBSOLETE]      }
// [OBSOLETE]      mass_electron_dos.at(spin_idx) /= valley_electron.at(spin_idx);
// [OBSOLETE]      mass_electron_dos.at(spin_idx) = std::pow((double) mass_electron_dos.at(spin_idx),(double)  1.0/3.0);
// [OBSOLETE]      if(osswrite) oss << setprecision(2) << mass_electron_dos.at(spin_idx);
// [OBSOLETE]    }
// [OBSOLETE]    if(osswrite) oss << "}" << endl;
// [OBSOLETE]    // write the DOS hole effective mass for all spins
// [OBSOLETE]    for (int spin_idx=0; spin_idx<ispin; spin_idx++) {
// [OBSOLETE]      if( spin_idx == 0 ) {
// [OBSOLETE]	if(osswrite) oss << "DOS hole effective mass (m0): {" ;
// [OBSOLETE]      } else if(spin_idx == 1) {
// [OBSOLETE]	if(osswrite) oss << ", ";
// [OBSOLETE]      }
// [OBSOLETE]      mass_hole_dos.at(spin_idx) /= valley_hole.at(spin_idx);
// [OBSOLETE]      mass_hole_dos.at(spin_idx) = std::pow((double) mass_hole_dos.at(spin_idx),(double)  1.0/3.0);
// [OBSOLETE]      if(osswrite) oss << setprecision(2) << mass_hole_dos.at(spin_idx) ;
// [OBSOLETE]    }
// [OBSOLETE]    if(osswrite) oss << "}" << endl;
// [OBSOLETE]    // write the conduction electron effective mass for all spins
// [OBSOLETE]    for (int spin_idx=0; spin_idx<ispin; spin_idx++) {
// [OBSOLETE]      if( spin_idx == 0 ) {
// [OBSOLETE]	if(osswrite) oss << "Cond. electron effective mass (m0): {" ;
// [OBSOLETE]      } else if(spin_idx == 1) {
// [OBSOLETE]	if(osswrite) oss << ", ";
// [OBSOLETE]      }
// [OBSOLETE]      if(osswrite) oss << setprecision(2) << mass_electron_conduction.at(spin_idx);
// [OBSOLETE]    }
// [OBSOLETE]    if(osswrite) oss << "}" << endl;
// [OBSOLETE]    // write the conduction hole effective mass for all spins
// [OBSOLETE]    for (int spin_idx=0; spin_idx<ispin; spin_idx++) {
// [OBSOLETE]      if( spin_idx == 0 ) {
// [OBSOLETE]	if(osswrite) oss << "Cond. hole effective mass (m0): {" ;
// [OBSOLETE]      } else if(spin_idx == 1) {
// [OBSOLETE]	if(osswrite) oss << ", ";
// [OBSOLETE]      }
// [OBSOLETE]      if(osswrite) oss << setprecision(2) << mass_hole_conduction.at(spin_idx) ;
// [OBSOLETE]    }
// [OBSOLETE]    if(osswrite) oss << "}" << endl;
// [OBSOLETE]  } else { // metal
// [OBSOLETE]    if(osswrite) oss << compound_name << endl;
// [OBSOLETE]    if(osswrite) oss << "Metal" << endl;
// [OBSOLETE]  }
// [OBSOLETE]
// [OBSOLETE]  // clean up
// [OBSOLETE]  for (uint i=0; i<destination.size(); i++) aurostd::RemoveFile(destination.at(i));
// [OBSOLETE]  
// [OBSOLETE]  fin_eigenval.close();
// [OBSOLETE] 
// [OBSOLETE]  if(ispin==1) { // if the calculation isn't spin polarized, remove the second entry
// [OBSOLETE]    mass_electron_dos.pop_back();
// [OBSOLETE]    mass_hole_dos.pop_back();
// [OBSOLETE]    mass_electron_conduction.pop_back();
// [OBSOLETE]    mass_hole_conduction.pop_back();
// [OBSOLETE]  }
// [OBSOLETE]  band_gap=band_information.band_gap;
// [OBSOLETE]
// [OBSOLETE]  // DEBUG: test the contents of the electron/hole up/down effective mass containers
// [OBSOLETE]  
// [OBSOLETE]    // cout << endl << endl << "Testing vector contents" << endl;
// [OBSOLETE]    // for ( vector<int>::size_type ix = 0; ix != mass_electron_dos.size(); ++ix) {
// [OBSOLETE]    // cout << "electron DoS m_{eff} ";
// [OBSOLETE]    // if( ix == 0 && mass_electron_dos.size() > 1 ) cout << "spin up   ";
// [OBSOLETE]    // if( ix == 1 && mass_electron_dos.size() > 1 ) cout << "spin down ";
// [OBSOLETE]    // cout << " = " << setprecision(10) << mass_electron_dos.at(ix) << endl;
// [OBSOLETE]    // }
// [OBSOLETE]    // for ( vector<int>::size_type ix = 0; ix != mass_hole_dos.size(); ++ix) {
// [OBSOLETE]    // cout << "hole DoS m_{eff} ";
// [OBSOLETE]    // if( ix == 0 && mass_electron_dos.size() > 1 ) cout << "spin up   "; 
// [OBSOLETE]    // if( ix == 1 && mass_electron_dos.size() > 1 ) cout << "spin down ";
// [OBSOLETE]    // cout << " = " << setprecision(10) << mass_hole_dos.at(ix) << endl;
// [OBSOLETE]    // }
// [OBSOLETE]   // for ( vector<int>::size_type ix = 0; ix != mass_electron_conduction.size(); ++ix) {
// [OBSOLETE]    // cout << "electron conductivity m_{eff} ";
// [OBSOLETE]    // if( ix == 0 && mass_electron_dos.size() > 1 ) cout << "spin up   "; 
// [OBSOLETE]    // if( ix == 1 && mass_electron_dos.size() > 1 ) cout << "spin down ";    
// [OBSOLETE]    // cout << " = " << setprecision(10) << mass_electron_conduction.at(ix) << endl;
// [OBSOLETE]    // }
// [OBSOLETE]    // for ( vector<int>::size_type ix = 0; ix != mass_hole_conduction.size(); ++ix) {
// [OBSOLETE]    // cout << "hole conductivity m_{eff} ";
// [OBSOLETE]    // if( ix == 0 && mass_electron_dos.size() > 1 ) cout << "spin up   "; 
// [OBSOLETE]    // if( ix == 1 && mass_electron_dos.size() > 1 ) cout << "spin down ";
// [OBSOLETE]    // cout << " = " << setprecision(10) << mass_hole_conduction.at(ix) << endl;
// [OBSOLETE]    // }
// [OBSOLETE]  
// [OBSOLETE]  
// [OBSOLETE]  return TRUE; // success!
// [OBSOLETE]}


// [OBSOLETE]bool _OLD_pflow::EffectiveMass(vector<string> &argv,
// [OBSOLETE]			   string directory_name,
// [OBSOLETE]			   ostream& oss) {
// [OBSOLETE]  // aflow --effective_mass directory_name
// [OBSOLETE]  // OR
// [OBSOLETE]  // aflow --em directory_name
// [OBSOLETE]  if(argv.size()) {;} // phony just to keep argv busy no complaining about unused
// [OBSOLETE]  double band_gap;
// [OBSOLETE]  vector<double> mass_electron_dos; // if size==1 nospin, if size==2 spin UP/DOWN
// [OBSOLETE]  vector<double> mass_hole_dos; // if size==1 nospin, if size==2 spin UP/DOWN
// [OBSOLETE]  vector<double> mass_electron_conduction; // if size==1 nospin, if size==2 spin UP/DOWN
// [OBSOLETE]  vector<double> mass_hole_conduction; // if size==1 nospin, if size==2 spin UP/DOWN
// [OBSOLETE]
// [OBSOLETE]  return _OLD_pflow::EffectiveMass(directory_name, ".static", band_gap, mass_electron_dos, mass_hole_dos, 
// [OBSOLETE]			       mass_electron_conduction, mass_hole_conduction, TRUE, oss);
// [OBSOLETE]}

// [OBSOLETE]bool fitToEllipsisEquation(vector<vector<kEn_st> > & fit_data, 
// [OBSOLETE]			   int spin_idx, 
// [OBSOLETE]			   vector<vector<double> > &mass_eff) {
// [OBSOLETE]  // **********************
// [OBSOLETE]  // least-square fitting to an ellipses equation
// [OBSOLETE]  // En = a x^2 + b y^2 + c z^2 + d xy + e xz + f yz + g x + h y + i z + j
// [OBSOLETE]  // to get rid of j, we fit the function
// [OBSOLETE]  // En1 - En2 = a (x1^2 - x2^2) + b (y1^2 - y2^2) + c (z1^2 - z2^2)
// [OBSOLETE]  //             + d (x1*y1 - x2*y2) + e (x1*z1 - x2*z2) + f (y1*z1 - y2*z2)
// [OBSOLETE]  //             + g (x1 - x2) + h (y1 - y2) + i (z1 - z2)
// [OBSOLETE]  for (uint i=0; i<fit_data.size(); i++) {
// [OBSOLETE]    kEn_st kp1 = fit_data.at(i).at(0);
// [OBSOLETE]    double kp1xx, kp1yy, kp1zz, kp1xy, kp1xz, kp1yz;
// [OBSOLETE]    kp1xx = kp1.kpoint(1) * kp1.kpoint(1); kp1yy = kp1.kpoint(2) * kp1.kpoint(2); kp1zz = kp1.kpoint(3) * kp1.kpoint(3);
// [OBSOLETE]   kp1xy = kp1.kpoint(1) * kp1.kpoint(2); kp1xz = kp1.kpoint(1) * kp1.kpoint(3); kp1yz = kp1.kpoint(2) * kp1.kpoint(3);
// [OBSOLETE]   
// [OBSOLETE]    //cout << "size " << fit_data.at(i).size() << endl;
// [OBSOLETE]    int nrow = fit_data.at(i).size()-1; // one less data
// [OBSOLETE]    int ncol = 9 ; // number of parameters in ellipses equation
// [OBSOLETE]    xvector<double> y_vec(1, nrow);
// [OBSOLETE]    xvector<double> y_sig(1, nrow);
// [OBSOLETE]    xmatrix<double> x_mat(1, 1, nrow, ncol);
// [OBSOLETE]    
// [OBSOLETE]    for (int j=1; j<nrow+1; j++) {
// [OBSOLETE]      y_vec[j] = fit_data.at(i).at(j).energy[spin_idx] - kp1.energy[spin_idx];
// [OBSOLETE]      y_sig[j] = _SIGMA;
// [OBSOLETE]    }
// [OBSOLETE]    
// [OBSOLETE]    for (int j=1; j<nrow+1; j++) {
// [OBSOLETE]      kEn_st kp = fit_data.at(i).at(j);
// [OBSOLETE]      x_mat[j][1] = kp.kpoint(1) * kp.kpoint(1) - kp1xx; // for a, see the description above
// [OBSOLETE]      x_mat[j][2] = kp.kpoint(2) * kp.kpoint(2) - kp1yy; // for b
// [OBSOLETE]      x_mat[j][3] = kp.kpoint(3) * kp.kpoint(3) - kp1zz; // for c
// [OBSOLETE]      x_mat[j][4] = kp.kpoint(1) * kp.kpoint(2) - kp1xy; // for d
// [OBSOLETE]      x_mat[j][5] = kp.kpoint(1) * kp.kpoint(3) - kp1xz; // for e
// [OBSOLETE]      x_mat[j][6] = kp.kpoint(2) * kp.kpoint(3) - kp1yz; // for f
// [OBSOLETE]      x_mat[j][7] = kp.kpoint(1) - kp1.kpoint(1);        // for g
// [OBSOLETE]      x_mat[j][8] = kp.kpoint(2) - kp1.kpoint(2);        // for h
// [OBSOLETE]      x_mat[j][9] = kp.kpoint(3) - kp1.kpoint(3);        // for i
// [OBSOLETE]    }
// [OBSOLETE]
// [OBSOLETE]    // check x_mat for columns of full of 0's; 
// [OBSOLETE]    // if there is a column of 0's then the matrix problem
// [OBSOLETE]    // [x_mat] [a,b,c,d,e,f,g,h,i] = [ y_vec ]
// [OBSOLETE]    // is ill-defined
// [OBSOLETE]    for ( int row = 1; row < 10; row++) {
// [OBSOLETE]      bool is_singular = true; // assume each entry in the colum will be 0
// [OBSOLETE]      for ( int column = 1; column < nrow+1; column++) {
// [OBSOLETE]	if( x_mat[column][row] != 0.0 ) { // test each entry
// [OBSOLETE]	  is_singular = false;             // if one non-zero entry is encountered, the column isn't singular
// [OBSOLETE]	}
// [OBSOLETE]      }
// [OBSOLETE]      if( is_singular ) { // report and exit because there will be errors
// [OBSOLETE]	cerr << "Column " << row << " in elliptical fitting is singular."  << endl;
// [OBSOLETE]	cerr << x_mat << endl;
// [OBSOLETE]	return FALSE;
// [OBSOLETE]      }
// [OBSOLETE]    }
// [OBSOLETE]
// [OBSOLETE]    aurostd::cematrix ECI_matrix(x_mat);
// [OBSOLETE]    ECI_matrix.LeastSquare(y_vec, y_sig);
// [OBSOLETE]    
// [OBSOLETE]    xmatrix<double> mass_m(1,1,3,3);
// [OBSOLETE]    mass_m[1][1] = ECI_matrix.AVec().at(0);
// [OBSOLETE]    mass_m[2][2] = ECI_matrix.AVec().at(1);
// [OBSOLETE]    mass_m[3][3] = ECI_matrix.AVec().at(2);
// [OBSOLETE]    mass_m[1][2] = ECI_matrix.AVec().at(3)*0.5;
// [OBSOLETE]    mass_m[1][3] = ECI_matrix.AVec().at(4)*0.5;
// [OBSOLETE]    mass_m[2][3] = ECI_matrix.AVec().at(5)*0.5;
// [OBSOLETE]    mass_m[2][1] = mass_m[1][2];
// [OBSOLETE]    mass_m[3][1] = mass_m[1][3];
// [OBSOLETE]    mass_m[3][2] = mass_m[2][3];
// [OBSOLETE]    
// [OBSOLETE]    // DEBUG:  print out the unique components of the effective mass tensor
// [OBSOLETE]    
// [OBSOLETE]      // cout << "mass_m[1][1]= " << ECI_matrix.AVec().at(0) << endl;
// [OBSOLETE]      // cout << "mass_m[2][2]= " << ECI_matrix.AVec().at(1) << endl;
// [OBSOLETE]      // cout << "mass_m[3][3]= " << ECI_matrix.AVec().at(2) << endl;
// [OBSOLETE]      // cout << "mass_m[1][2]= " << ECI_matrix.AVec().at(3)*0.5 << endl;
// [OBSOLETE]      // cout << "mass_m[1][3]= " << ECI_matrix.AVec().at(4)*0.5 << endl;
// [OBSOLETE]      // cout << "mass_m[2][3]= " << ECI_matrix.AVec().at(5)*0.5 << endl;
// [OBSOLETE]    
// [OBSOLETE]
// [OBSOLETE]    aurostd::cematrix mass_m_ce(mass_m);
// [OBSOLETE]    xvector<double> mr = mass_m_ce.EigenValues();
// [OBSOLETE]    vector<double> mr_tmp;
// [OBSOLETE]
// [OBSOLETE]    for (int i=1; i<=3; i++) {
// [OBSOLETE]      //cerr << "mr[" << i << "]=" << mr[i] << endl;
// [OBSOLETE]      mr[i] =  1.0*_MASS_FACTOR/mr[i];
// [OBSOLETE]      mr_tmp.push_back(mr[i]);
// [OBSOLETE]    }
// [OBSOLETE]    sort(mr_tmp.begin(), mr_tmp.end());
// [OBSOLETE]    mass_eff.push_back(mr_tmp);
// [OBSOLETE]  }
// [OBSOLETE]
// [OBSOLETE]  // DEBUG:  print out all the effective masses
// [OBSOLETE]  
// [OBSOLETE]     // for ( vector<int>::size_type ix = 0; ix != mass_eff.size(); ++ix) {
// [OBSOLETE]     // int i = ix;
// [OBSOLETE]     // for ( vector<int>::size_type iy = 0; iy != mass_eff.at(ix).size(); ++iy) {
// [OBSOLETE]     // int j = iy;
// [OBSOLETE]     // cerr << "mass_eff[" << i << "][" << j << "]=" << mass_eff.at(ix).at(iy) << endl;
// [OBSOLETE]     // }
// [OBSOLETE]     // }
// [OBSOLETE] 
// [OBSOLETE]  return TRUE;
// [OBSOLETE]}
// [OBSOLETE]
// [OBSOLETE]bool getEfermiBandGap(string & file_dos, int ispin, string & compound_name, band_st &band_information ) {
// [OBSOLETE]
// [OBSOLETE]  double energy_min, energy_max, Epoint_number, Efermi, atmp;
// [OBSOLETE]  ifstream fin_doscar;
// [OBSOLETE]
// [OBSOLETE]  double dos_up_t, dos_down_t;   // dos
// [OBSOLETE]  double idos_up_t, idos_down_t; // integrated dos
// [OBSOLETE]  double _EQUAL_DB_HERE = 1.0e-6;
// [OBSOLETE]
// [OBSOLETE]  bool vbt_up_flag=false, cbb_up_flag=false; //valence band top and conduction band bottom
// [OBSOLETE]  bool vbt_down_flag=false, cbb_down_flag=false; //valence band top and conduction band bottom
// [OBSOLETE]  double vbt_up, cbb_up; //valence band top and conduction band bottom
// [OBSOLETE]  double vbt_down, cbb_down; //valence band top and conduction band bottom
// [OBSOLETE]
// [OBSOLETE]  string line_content;
// [OBSOLETE]  double energy;
// [OBSOLETE]
// [OBSOLETE]  stringstream ss;
// [OBSOLETE]    
// [OBSOLETE]  fin_doscar.open(file_dos.c_str());
// [OBSOLETE]  // skip some lines of header
// [OBSOLETE]  for (uint i=0; i<4; i++) { 
// [OBSOLETE]    getline(fin_doscar, line_content);
// [OBSOLETE]  }
// [OBSOLETE]
// [OBSOLETE]  fin_doscar >> compound_name; // get compound name
// [OBSOLETE]
// [OBSOLETE]  fin_doscar >> energy_max >> energy_min >> Epoint_number >> Efermi >> atmp;
// [OBSOLETE]  // DEBUG:
// [OBSOLETE]  //cerr << "ispin=" << ispin << endl;
// [OBSOLETE]  //cerr << "getEfermiBandGap Efermi=" << Efermi << endl;
// [OBSOLETE]  //cerr << "Epoint_number=" << Epoint_number << endl;
// [OBSOLETE]
// [OBSOLETE]  // loop over the points of energy in DOSCAR
// [OBSOLETE]  for(uint i=0; i<Epoint_number; i++) {
// [OBSOLETE]    getline(fin_doscar, line_content);
// [OBSOLETE]    //cerr << i << line_content << endl;
// [OBSOLETE]    ss.clear();
// [OBSOLETE]    ss.str(line_content);
// [OBSOLETE]    if(ispin == 1) { // on spin polarization
// [OBSOLETE]      ss >> energy >> dos_up_t >> idos_up_t;
// [OBSOLETE]      dos_down_t = dos_up_t;
// [OBSOLETE]    } else if(ispin == 2) {
// [OBSOLETE]      ss >> energy >> dos_up_t >> dos_down_t >> idos_up_t >> idos_down_t;
// [OBSOLETE]    }
// [OBSOLETE]    
// [OBSOLETE]    // only work when energy is ordered abscendently
// [OBSOLETE]    // vbt_*_flag, cbb_*_flag are true when the band gap is larger than zero
// [OBSOLETE]    if( energy <= Efermi) { 
// [OBSOLETE]      if( abs(dos_up_t) < _EQUAL_DB_HERE ) {
// [OBSOLETE]	vbt_up_flag = true;
// [OBSOLETE]      } else {
// [OBSOLETE]	vbt_up = energy;
// [OBSOLETE]	vbt_up_flag = false;
// [OBSOLETE]	cbb_up_flag = false;
// [OBSOLETE]      }
// [OBSOLETE]      if(abs(dos_down_t) < _EQUAL_DB_HERE ) {
// [OBSOLETE]	vbt_down_flag = true;
// [OBSOLETE]      } else {
// [OBSOLETE]	vbt_down = energy;
// [OBSOLETE]	vbt_down_flag = false;
// [OBSOLETE]	cbb_down_flag = false;
// [OBSOLETE]      }
// [OBSOLETE]    } else {
// [OBSOLETE]      if(abs(dos_up_t) < _EQUAL_DB_HERE) {
// [OBSOLETE]	vbt_up_flag = true;
// [OBSOLETE]      } else {
// [OBSOLETE]	if( ! cbb_up_flag ) cbb_up = energy;
// [OBSOLETE]	cbb_up_flag = true;
// [OBSOLETE]      }	
// [OBSOLETE]      if(abs(dos_down_t) < _EQUAL_DB_HERE ) {
// [OBSOLETE]	vbt_down_flag = true;
// [OBSOLETE]      } else {
// [OBSOLETE]	if( ! cbb_down_flag ) cbb_down = energy;
// [OBSOLETE]	cbb_down_flag = true;
// [OBSOLETE]      }
// [OBSOLETE]    }
// [OBSOLETE]  }
// [OBSOLETE]  band_information.Efermi = Efermi;
// [OBSOLETE]  for (int spin_idx=0; spin_idx<2; spin_idx++) {
// [OBSOLETE]    band_information.vbt[spin_idx] = Efermi;
// [OBSOLETE]    band_information.cbb[spin_idx] = Efermi;
// [OBSOLETE]  }
// [OBSOLETE]    
// [OBSOLETE]  if( ispin == 1 ) {
// [OBSOLETE]    if( vbt_up_flag && cbb_up_flag ) {
// [OBSOLETE]      band_information.vbt[0] = vbt_up;
// [OBSOLETE]      band_information.cbb[0] = cbb_up;
// [OBSOLETE]      band_information.vbt[1] = vbt_up;
// [OBSOLETE]      band_information.cbb[1] = cbb_up;
// [OBSOLETE]    }
// [OBSOLETE]  } else if( ispin == 2 ) {
// [OBSOLETE]    if( vbt_up_flag && cbb_up_flag ) {
// [OBSOLETE]      band_information.vbt[0] = vbt_up;
// [OBSOLETE]      band_information.cbb[0] = cbb_up;
// [OBSOLETE]    }
// [OBSOLETE]    if( vbt_down_flag && cbb_down_flag ) {
// [OBSOLETE]      band_information.vbt[1] = vbt_down;
// [OBSOLETE]      band_information.cbb[1] = cbb_down;
// [OBSOLETE]    }
// [OBSOLETE]  }
// [OBSOLETE]  // system band gap is the smaller of the spins' band gap
// [OBSOLETE]  band_information.band_gap = std::min( band_information.cbb[0] - band_information.vbt[0], band_information.cbb[1] - band_information.vbt[1]);
// [OBSOLETE]  // if negative, the valence band energy was larger than conduction band and it's a metal
// [OBSOLETE]  band_information.band_gap = std::max(0.0, band_information.band_gap);
// [OBSOLETE]    
// [OBSOLETE]  return TRUE;
// [OBSOLETE]}


// [OBSOLETE]bool getFitDataForMass(const string & file_kE, 
// [OBSOLETE]		       const string & file_lattice,
// [OBSOLETE]		       const band_st & band_information,
// [OBSOLETE]		       vector<vector<int> > & number_of_valley_list, 
// [OBSOLETE]		       double energy_range,
// [OBSOLETE]		       vector<vector<vector<kEn_st> > > &fit_data_all) {  // get the fit data for effective mass calculations
// [OBSOLETE]  // the data are those kpoints with energy within energy_range of the extremum.
// [OBSOLETE]
// [OBSOLETE]  ifstream fin_poscar, fin_eigenval;
// [OBSOLETE]  fin_poscar.open(file_lattice.c_str());
// [OBSOLETE]  fin_eigenval.open(file_kE.c_str());
// [OBSOLETE]    
// [OBSOLETE]  string line_content;
// [OBSOLETE]
// [OBSOLETE]  // ****************************
// [OBSOLETE]  // Get lattices from POSCAR
// [OBSOLETE]  xmatrix<double> lattice(1,1,3,3), reciprocal_lattice(1,1,3,3);
// [OBSOLETE]  xstructure xstr(fin_poscar, IOVASP_POSCAR) ;
// [OBSOLETE]  fin_poscar.close();
// [OBSOLETE]
// [OBSOLETE]  xstr.FixLattices();
// [OBSOLETE]  xstr.CalculateSymmetryPointGroupKlattice();
// [OBSOLETE]  reciprocal_lattice = xstr.klattice;
// [OBSOLETE]
// [OBSOLETE]  // ****************************
// [OBSOLETE]  // Get energy dispersion relations
// [OBSOLETE]  // get the spin polariation
// [OBSOLETE]  // consume 1st line in EIGENVAL
// [OBSOLETE]  int Nions, Ntmp, ispin;
// [OBSOLETE]  fin_eigenval >> Nions >> Nions >> Ntmp >> ispin;
// [OBSOLETE]
// [OBSOLETE]  vector<double> Efermi_bandGap;
// [OBSOLETE]
// [OBSOLETE]
// [OBSOLETE]  stringstream ss;
// [OBSOLETE]  for (int i=0; i<5; i++) { // ignore lines 2-5 of the header
// [OBSOLETE]    getline(fin_eigenval, line_content);
// [OBSOLETE]  }
// [OBSOLETE]
// [OBSOLETE]  double kpoint_number, band_number, a_tmp;
// [OBSOLETE]  fin_eigenval >> a_tmp >> kpoint_number >> band_number;
// [OBSOLETE]  
// [OBSOLETE]  vector<vector<kEn_st> > allkE_points;
// [OBSOLETE]  allkE_points.resize(band_number);
// [OBSOLETE]  
// [OBSOLETE]  xvector<double> kp_x(1,3), kp_car(1,3);
// [OBSOLETE]  
// [OBSOLETE]  for (int i=0; i<kpoint_number; i++) {
// [OBSOLETE]    getline(fin_eigenval, line_content);
// [OBSOLETE]    kEn_st kp;
// [OBSOLETE]    double tmp_noneed;
// [OBSOLETE]    fin_eigenval >> kp.kpoint(1) >> kp.kpoint(2) >> kp.kpoint(3) >> tmp_noneed;
// [OBSOLETE]    kp_x[1] = kp.kpoint(1); kp_x[2] = kp.kpoint(2); kp_x[3] = kp.kpoint(3);
// [OBSOLETE]    kp_car = kp_x *reciprocal_lattice; // to Cartesian coordinates
// [OBSOLETE]    kp.kpoint(1) = kp_car[1]; kp.kpoint(2) = kp_car[2]; kp.kpoint(3) = kp_car[3];
// [OBSOLETE]    for (int j=0; j<band_number; j++) {
// [OBSOLETE]      if( ispin == 1 ) {
// [OBSOLETE]	fin_eigenval >> kp.band_index >> kp.energy[0];
// [OBSOLETE]	kp.energy[1] = kp.energy[0];
// [OBSOLETE]      } else {
// [OBSOLETE]	fin_eigenval >> kp.band_index >> kp.energy[0] >> kp.energy[1];
// [OBSOLETE]      }
// [OBSOLETE]      allkE_points.at(kp.band_index-1).push_back(kp); // in EIGENVAL band index begins with 1
// [OBSOLETE]    }
// [OBSOLETE]  }
// [OBSOLETE]  
// [OBSOLETE]  // allkE_points contains the data for 1 band in each entry
// [OBSOLETE]  //    size is no of bands
// [OBSOLETE]
// [OBSOLETE]  // allkE_points.at(i) is the vector of all kpoints of the ith band   
// [OBSOLETE]  //    size is no of kpoints
// [OBSOLETE]
// [OBSOLETE]  // allkE_points.at(i).at(j) is the data for the jth kpoint of the ith band
// [OBSOLETE]  //    data is special data structure kEn_st defined in aflow_contrib_shidong.h
// [OBSOLETE]  
// [OBSOLETE]  for (int spin_idx=0; spin_idx<ispin; spin_idx++) {
// [OBSOLETE]    
// [OBSOLETE]    vector<int> number_of_valley_list_tmp;
// [OBSOLETE]    
// [OBSOLETE]    //**********************
// [OBSOLETE]    // choose bands containing extremes
// [OBSOLETE]    // first sort the energy in ascending order
// [OBSOLETE]    for (uint i=0; i<allkE_points.size(); i++) {
// [OBSOLETE]      if(spin_idx == 0) {
// [OBSOLETE]	// spin up
// [OBSOLETE]	sort (allkE_points.at(i).begin(), allkE_points.at(i).end(), comparison_kEn_str_up);
// [OBSOLETE]      } else {
// [OBSOLETE]	// spin down
// [OBSOLETE]	sort (allkE_points.at(i).begin(), allkE_points.at(i).end(), comparison_kEn_str_down);
// [OBSOLETE]      }
// [OBSOLETE]    }
// [OBSOLETE]      
// [OBSOLETE]    vector< vector<kEn_st> > fit_data;
// [OBSOLETE]    
// [OBSOLETE]    // then get the points with energy in the range
// [OBSOLETE]    for (uint i=0; i<allkE_points.size(); i++) { // loop over no. of bands
// [OBSOLETE]      vector<kEn_st> fit_data_band;
// [OBSOLETE]      
// [OBSOLETE]      if( abs(allkE_points.at(i).front().energy[spin_idx] - band_information.cbb[spin_idx]) < energy_range ) {
// [OBSOLETE]	// if the ith band's lowest energy is within energy_range of the conduction band bottom
// [OBSOLETE]	// allkE_points.at(i) is sorted by energy, lowest first.
// [OBSOLETE]	double band_energy_minimum = allkE_points.at(i).front().energy[spin_idx];
// [OBSOLETE]	for (int j=0; j<_FIT_POINTS_NUMBER; j++) { //  collect at least _FIT_POINTS_NUMBER of data
// [OBSOLETE]	  //cerr << "storing point with kz=" << allkE_points.at(i).at(j).kpoint(3) << endl;
// [OBSOLETE]	  fit_data_band.push_back(allkE_points.at(i).at(j));
// [OBSOLETE]	  fit_data_band.back().band_type = 1; // conduction band
// [OBSOLETE]	}
// [OBSOLETE]	  
// [OBSOLETE]	for (uint j=_FIT_POINTS_NUMBER; j<allkE_points.at(i).size(); j++) { // 
// [OBSOLETE]	  if( allkE_points.at(i).at(j).energy[spin_idx] - band_energy_minimum < _FIT_ENERGY_RANGE ) {
// [OBSOLETE]	    
// [OBSOLETE]	    fit_data_band.push_back(allkE_points.at(i).at(j));
// [OBSOLETE]	    fit_data_band.back().band_type = 1; // conduction band
// [OBSOLETE]	  } else {
// [OBSOLETE]	    break;
// [OBSOLETE]	  }
// [OBSOLETE]	}
// [OBSOLETE]      } else if( abs(band_information.vbt[spin_idx] - allkE_points.at(i).back().energy[spin_idx]) < energy_range) {
// [OBSOLETE]	// valence bands
// [OBSOLETE]	double band_energy_maximum = allkE_points.at(i).back().energy[spin_idx];
// [OBSOLETE]	for (uint j=allkE_points.at(i).size()-1; j>allkE_points.at(i).size() - _FIT_POINTS_NUMBER-1; j--) {
// [OBSOLETE]	  fit_data_band.push_back(allkE_points.at(i).at(j));
// [OBSOLETE]	  fit_data_band.back().band_type = 0; // valence band
// [OBSOLETE]	}
// [OBSOLETE]	
// [OBSOLETE]	for (int j=allkE_points.at(i).size() - _FIT_POINTS_NUMBER-1; j>=0; j--) {
// [OBSOLETE]	  if( band_energy_maximum - allkE_points.at(i).at(j).energy[spin_idx]  < _FIT_ENERGY_RANGE ) {
// [OBSOLETE]	    fit_data_band.push_back(allkE_points.at(i).at(j));
// [OBSOLETE]	    fit_data_band.back().band_type = 0; // valence band
// [OBSOLETE]	  } else {
// [OBSOLETE]	    break;
// [OBSOLETE]	  }
// [OBSOLETE]	}
// [OBSOLETE]      }
// [OBSOLETE]      
// [OBSOLETE]      //if(fit_data_band.size() > 0 ) fit_data.push_back(fit_data_band);
// [OBSOLETE]      if(fit_data_band.size() > 0 ) {
// [OBSOLETE]	fit_data.push_back(fit_data_band);
// [OBSOLETE]      }
// [OBSOLETE]    }
// [OBSOLETE]      
// [OBSOLETE]    //**********************
// [OBSOLETE]    // find all points related by symmetric operators
// [OBSOLETE]    vector<double> max_distance;
// [OBSOLETE]    for (int i=1; i<4; i++) {
// [OBSOLETE]      double max_tmp;
// [OBSOLETE]      max_tmp = std::max(abs(reciprocal_lattice[1][i]), abs(reciprocal_lattice[2][i]));
// [OBSOLETE]      max_tmp = std::max(abs(reciprocal_lattice[3][i]), max_tmp);
// [OBSOLETE]      max_tmp *= _MIN_RATIO;
// [OBSOLETE]      max_distance.push_back(max_tmp);
// [OBSOLETE]    }
// [OBSOLETE]      
// [OBSOLETE]    vector<vector<kEn_st> > fit_data_new;
// [OBSOLETE]    for (uint i=0; i< fit_data.size(); i++) {
// [OBSOLETE]      vector<kEn_st> fit_data_band;
// [OBSOLETE]      kEn_st kp;
// [OBSOLETE]      xvector<double> pt(1,3);
// [OBSOLETE]      kp = fit_data.at(i).at(0); // the first point is closest to the extremes
// [OBSOLETE]      pt[1] = kp.kpoint(1);
// [OBSOLETE]      pt[2] = kp.kpoint(2);
// [OBSOLETE]      pt[3] = kp.kpoint(3);
// [OBSOLETE]      for (uint ii=0; ii<fit_data.at(i).size(); ii++) {
// [OBSOLETE]	xvector<double> pt1(1,3), pt_sym(1,3);
// [OBSOLETE]	kEn_st kp1;
// [OBSOLETE]	kp1 = fit_data.at(i).at(ii); // the first point is closest to the extremes
// [OBSOLETE]	pt1[1] = kp1.kpoint(1); pt1[2] = kp1.kpoint(2); pt1[3] = kp1.kpoint(3);
// [OBSOLETE]	for (uint j=0; j<xstr.pgroupk.size(); j++) {
// [OBSOLETE]	  pt_sym = pt1 * xstr.pgroupk.at(j).Uc;
// [OBSOLETE]	  if(near_to(pt, pt_sym, max_distance)) {
// [OBSOLETE]	    // compare the distance between the most extreme points and the generated one
// [OBSOLETE]	    kEn_st k1;
// [OBSOLETE]	    k1.kpoint(1) = pt_sym[1];
// [OBSOLETE]	    k1.kpoint(2) = pt_sym[2];
// [OBSOLETE]	    k1.kpoint(3) = pt_sym[3];
// [OBSOLETE]	    k1.energy[0] = kp1.energy[0];
// [OBSOLETE]	    k1.energy[1] = kp1.energy[1];
// [OBSOLETE]	    k1.band_index = kp1.band_index;
// [OBSOLETE]	    k1.band_type = kp1.band_type;
// [OBSOLETE]	    fit_data_band.push_back(k1);
// [OBSOLETE]	  }
// [OBSOLETE]	}
// [OBSOLETE]	if( ii == 0 ) {
// [OBSOLETE]	  int number_of_valley = xstr.pgroupk.size()/fit_data_band.size();
// [OBSOLETE]	  number_of_valley_list_tmp.push_back(number_of_valley);
// [OBSOLETE]	}
// [OBSOLETE]      }
// [OBSOLETE]      
// [OBSOLETE]      vector<kEn_st>::iterator it;
// [OBSOLETE]      
// [OBSOLETE]      sort (fit_data_band.begin(), fit_data_band.end(), comparison_kEn_str_position);
// [OBSOLETE]      
// [OBSOLETE]      it = unique(fit_data_band.begin(), fit_data_band.end(), is_equal_position_kEn_str);
// [OBSOLETE]      fit_data_band.resize(it - fit_data_band.begin());
// [OBSOLETE]      
// [OBSOLETE]      if(spin_idx == 1) {
// [OBSOLETE]	// spin up
// [OBSOLETE]	sort (fit_data_band.begin(), fit_data_band.end(), comparison_kEn_str_band_type_up);
// [OBSOLETE]      } else {
// [OBSOLETE]	// spin down
// [OBSOLETE]	sort (fit_data_band.begin(), fit_data_band.end(), comparison_kEn_str_band_type_down);
// [OBSOLETE]      }
// [OBSOLETE]      
// [OBSOLETE]      fit_data_new.push_back(fit_data_band);
// [OBSOLETE]    }
// [OBSOLETE]    
// [OBSOLETE]    fit_data.clear();
// [OBSOLETE]    fit_data = fit_data_new;
// [OBSOLETE]    
// [OBSOLETE]    fit_data_all.push_back(fit_data);
// [OBSOLETE]    number_of_valley_list.push_back(number_of_valley_list_tmp);
// [OBSOLETE]  }
// [OBSOLETE]  
// [OBSOLETE]  fin_eigenval.close();
// [OBSOLETE]  
// [OBSOLETE]  return TRUE;
// [OBSOLETE]}

#endif
