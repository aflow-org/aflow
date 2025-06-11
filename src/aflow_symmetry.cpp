// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
//
// Acknowledgements to Dane Morgan and Anton Van der Ven for discussions and
// suggestions.

#include "aflow.h"

//#include "xlibs.h"

#define cdebug cerr
//#define _EPS_ 0.005   // messes up JX monoclinics
#define _EPS_ 0.02 //JX
// #define _EPS_ 0.001 // shidong
//#define _EPS_ 0.02  // seems to work in the mean time of a self correcting one
#define _EPS_roundoff_ _DOUBLE_TOL_ //CO20200731 //1.0e-8
#define COMPILE_SLIM
#define DEBUG_MINIMUM_DISTANCE false //DX+CO20210114

using aurostd::isdiagonal;
//DX+CO START
using aurostd::isidentity;
//DX+CO END
using aurostd::isequal;

//DX+CO START
#include "aflow_symmetry_spacegroup.h"
//DX+CO END

#define  _NO_SCALE_LATTICE_PGROUP_

bool DEBUG_SYMMETRY=FALSE;

// ============================================================================
//DX+CO START
// ----------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------- GENERAL SYMMETRY FUNCTIONS
// These functions are used elsewhere in AFLOW, namely 
// the space group files. DO NOT DELETE. 
// If these functions become obsolete in this file, please move 
// over to aflow_symmetry_spacegroup_functions.cpp and 
// aflow_symmetry_spacegroup.h.  - David Hicks (d.hicks@duke.edu)

// **********************************************************************************************************************
// minimumDistance
// **********************************************************************************************************************
// Find the minimum interatomic distance within a given crystal structure 
namespace SYM {
  double minimumDistance(const xstructure& xstr){
    //xstructure a(xstr); a.ReScale(1.0);
    //xmatrix<double> lattice=xstr.scale*xstr.lattice;
    // -------------------------------------------------------------------
    // need to copy atoms and bring in cell, otherwise minimum distance
    // calculation may not work
    deque<_atom> atoms = xstr.atoms; //DX20210502
    uint natoms=atoms.size(); //DX20210502
    for(uint i=0;i<natoms;i++){ BringInCellInPlace(atoms[i], xstr.lattice); } //DX20210502
    return minimumDistance(atoms,xstr.lattice,xstr.scale);
  }
} // namespace SYM

namespace SYM {
  double minimumDistance(const deque<_atom>& atoms){  //CO20190808
    //for NON periodic systems, use the default with lattice otherwise: minimumDistance(const deque<_atom>& atoms,const xmatrix<double>& lattice,double scale)
    double dist=0,dist_min=AUROSTD_MAX_DOUBLE;
    for(uint i=0;i<atoms.size()-1;i++){
      for(uint j=i+1;j<atoms.size();j++){
        dist=aurostd::modulus(atoms[i].cpos-atoms[j].cpos);
        if(dist<dist_min){dist_min=dist;}
      }
    }
    return dist_min;
  }
  double minimumDistance(const deque<_atom>& atoms,const xmatrix<double>& lattice,double scale){

    //DX20201130 - commented nested LDEBUG if-statements, notable speed increase (a few seconds for systems with >50 atoms)
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    double min_dist=AUROSTD_MAX_DOUBLE;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " INITIAL [0] " << min_dist << endl;

    double radius=RadiusSphereLattice(lattice);
    xvector<int> dims(3);
    dims=LatticeDimensionSphere(lattice,radius);
    int dim=max(dims);

    //length of lattice vectors
    //for(int a=1;a<=dim;a++){
    for(int i=1;i<=3;i++){
      min_dist = aurostd::min(min_dist,aurostd::modulus(lattice(i))); // if loop is brought back, do a*lattice(i)
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " lattice_vector[" << i << "]: min_dist=" << min_dist << endl;}
    }
    //}
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " LATTICE VECTORS [1] " << min_dist << endl;

    //DX20180508 - FASTER MIN CART DISTANCE CALCULATOR - START
    //DX20180508 - only calculate multiplication once (time-saver)
    vector<vector<xvector<double> > > lattice_lengths;
    vector<vector<int> > lattice_indices;
    vector<xvector<double> > latt;
    vector<int> index;
    for(int a=1;a<=dim;a++){latt.push_back(a*lattice(1));index.push_back(a);} //DX calc once and store
    lattice_lengths.push_back(latt); latt.clear();
    lattice_indices.push_back(index); index.clear();
    for(int b=1;b<=dim;b++){latt.push_back(b*lattice(2));index.push_back(b);} //DX calc once and store
    lattice_lengths.push_back(latt); latt.clear();
    lattice_indices.push_back(index); index.clear();
    for(int c=1;c<=dim;c++){latt.push_back(c*lattice(3));index.push_back(c);} //DX calc once and store
    lattice_lengths.push_back(latt); latt.clear();
    lattice_indices.push_back(index); index.clear();

    //combos of two lattice vectors
    for(uint i=0;i<lattice_lengths.size()-1;i++){ //CO20190808 - .size()-1
      for(uint j=i+1;j<lattice_lengths.size();j++){    
        for(uint a=0;a<lattice_lengths[i].size();a++){
          for(uint b=0;b<lattice_lengths[j].size();b++){
            min_dist = aurostd::min(min_dist,aurostd::modulus(lattice_lengths[i][a]+lattice_lengths[j][b]));
            min_dist = aurostd::min(min_dist,aurostd::modulus(lattice_lengths[i][a]-lattice_lengths[j][b]));
#if DEBUG_MINIMUM_DISTANCE
            if(LDEBUG) {cerr << __AFLOW_FUNC__ << " lattice_vectors: a=" << lattice_indices[i][a] << ",b=" << lattice_indices[j][b] << "; min_dist=" << min_dist << endl;}
#endif
          }
        }
      }
    }
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " COMBOS OF 2 LATTICE VECTORS [2] " << min_dist << endl;

    ////combos of two lattice vectors
    //for(int a=1;a<=dim;a++){
    //  for(int b=1;b<=dim;b++){
    //    for(int i=1;i<=3;i++){
    //      for(int j=i+1;j<=3;j++){
    //        min_dist = aurostd::min(min_dist,aurostd::modulus(a*lattice(i)+b*lattice(j)));
    //        min_dist = aurostd::min(min_dist,aurostd::modulus(a*lattice(i)-b*lattice(j)));
    //        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " lattice_vectors: a=" << a << ",b=" << b << "; min_dist=" << min_dist << endl;}
    //      }
    //    }
    //  }
    //}
    //if(LDEBUG) cerr << __AFLOW_FUNC__ << " COMBOS OF 2 LATTICE VECTORS [2] " << min_dist << endl;

    //combos of three lattice vectors 
    lattice_lengths.clear();
    lattice_indices.clear();
    latt.clear();
    index.clear();
    for(int a=1;a<=dims[1];a++){latt.push_back(a*lattice(1));index.push_back(a);} //DX calc once and store
    lattice_lengths.push_back(latt); latt.clear();
    lattice_indices.push_back(index); index.clear();
    for(int b=1;b<=dims[2];b++){latt.push_back(b*lattice(2));index.push_back(b);} //DX calc once and store
    lattice_lengths.push_back(latt); latt.clear();
    lattice_indices.push_back(index); index.clear();
    for(int c=1;c<=dims[3];c++){latt.push_back(c*lattice(3));index.push_back(c);} //DX calc once and store
    lattice_lengths.push_back(latt); latt.clear();
    lattice_indices.push_back(index); index.clear();

    //combos of three lattice vectors 
    for(uint i=0;i<lattice_lengths.size()-1;i++){ //CO20190808 - .size()-1
      for(uint j=i+1;j<lattice_lengths.size();j++){    
        for(uint a=0;a<lattice_lengths[i].size();a++){
          for(uint b=0;b<lattice_lengths[j].size();b++){
            xvector<double> added = lattice_lengths[i][a]+lattice_lengths[j][b];
            xvector<double> subtracted = lattice_lengths[i][a]+lattice_lengths[j][b];
            for(uint k=j+1;k<lattice_lengths.size();k++){    
              for(uint c=0;c<lattice_lengths[k].size();c++){
                min_dist = aurostd::min(min_dist,aurostd::modulus(added+lattice_lengths[k][c]));
                min_dist = aurostd::min(min_dist,aurostd::modulus(subtracted-lattice_lengths[k][c]));
                min_dist = aurostd::min(min_dist,aurostd::modulus(added-lattice_lengths[k][c]));
                min_dist = aurostd::min(min_dist,aurostd::modulus(subtracted+lattice_lengths[k][c]));
#if DEBUG_MINIMUM_DISTANCE
                if(LDEBUG) {cerr << __AFLOW_FUNC__ << " lattice_vectors: a=" << lattice_indices[i][a] << ",b=" << lattice_indices[j][b] << ",c=" << lattice_indices[k][c] << "; min_dist=" << min_dist << endl;}
#endif
              }
            }
          }  
        }
      }
    }


    //for(int a=1;a<=dims[1];a++){
    //  for(int b=1;b<=dims[2];b++){
    //    for(int c=1;c<=dims[3];c++){
    //      for(int i=1;i<=3;i++){
    //        for(int j=i+1;j<=3;j++){
    //          for(int k=j+1;k<=3;k++){
    //            min_dist = aurostd::min(min_dist,aurostd::modulus(a*lattice(i)+b*lattice(j)+c*lattice(k)));
    //            min_dist = aurostd::min(min_dist,aurostd::modulus(a*lattice(i)-b*lattice(j)-c*lattice(k)));
    //            min_dist = aurostd::min(min_dist,aurostd::modulus(a*lattice(i)+b*lattice(j)-c*lattice(k)));
    //            min_dist = aurostd::min(min_dist,aurostd::modulus(a*lattice(i)-b*lattice(j)+c*lattice(k)));
    //            if(LDEBUG) {cerr << __AFLOW_FUNC__ << " lattice_vectors: a=" << a << ",b=" << b << ",c=" << c << "; min_dist=" << min_dist << endl;}
    //            //min_dist = aurostd::min(min_dist,aurostd::modulus(-lattice(i)-lattice(j)-lattice(k)));  //same as -(+++)
    //            //min_dist = aurostd::min(min_dist,aurostd::modulus(-lattice(i)+lattice(j)+lattice(k)));  //same as -(+--)
    //            //min_dist = aurostd::min(min_dist,aurostd::modulus(-lattice(i)-lattice(j)+lattice(k)));  //same as -(++-)
    //            //min_dist = aurostd::min(min_dist,aurostd::modulus(-lattice(i)+lattice(j)-lattice(k)));  //same as -(+-+)
    //          }
    //        }
    //      }  
    //    }
    //  }
    //}
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " COMBOS OF 3 LATTICE VECTORS [3] " << min_dist << endl;

    //distance between each atom
    //DX20171023
    // Since we are finding the minimum interatomic distance, we want to find a distance smaller than the length of the minimum lattice vector 
    // distance.  Centering on each atom, we find the lattice dimensions sphere of the minimum distance.  If any atom is closer 
    // than the minimum distance, it becomes the radius for the new lattice dimensions sphere for the next iteration. 


    //DX20180508 - FASTER MIN CART DISTANCE CALCULATOR - START
    //DX20180508 - only calculate multiplication once (time-saver)
    vector<xvector<double> > l1, l2, l3;
    vector<int> a_index, b_index, c_index;
    for(int a=-dims[1];a<=dims[1];a++){l1.push_back(a*lattice(1));a_index.push_back(a);} //DX calc once and store
    for(int b=-dims[2];b<=dims[2];b++){l2.push_back(b*lattice(2));b_index.push_back(b);} //DX calc once and store
    for(int c=-dims[3];c<=dims[3];c++){l3.push_back(c*lattice(3));c_index.push_back(c);} //DX calc once and store

    xvector<double> tmp;

    for(uint i=0; i<atoms.size()-1; i++){ //CO20190808 - .size()-1
      // Cannot reduce more than (1,1,1), so don't recalculate
      if(!(dims[1]==1 && dims[2]==1 && dims[3]==1)){
        dims=LatticeDimensionSphere(lattice,min_dist);
        l1.clear(); l2.clear(); l3.clear();
        a_index.clear(); b_index.clear(); c_index.clear();
        for(int a=-dims[1];a<=dims[1];a++){l1.push_back(a*lattice(1));a_index.push_back(a);} //DX calc once and store
        for(int b=-dims[2];b<=dims[2];b++){l2.push_back(b*lattice(2));b_index.push_back(b);} //DX calc once and store
        for(int c=-dims[3];c<=dims[3];c++){l3.push_back(c*lattice(3));c_index.push_back(c);} //DX calc once and store
        //cerr << "dims: " << dims << endl;
      }
      for(uint k=i+1; k<atoms.size(); k++){
        xvector<double> incell_dist = atoms[k].cpos-atoms[i].cpos;
        //DX20180423 - running vector in each loop saves computations; fewer duplicate operations
        for(uint m=0;m<l1.size();m++){
          xvector<double> a_component = incell_dist + l1[m];    //DX : coord1-coord2+a*lattice(1)
          for(uint n=0;n<l2.size();n++){
            xvector<double> ab_component = a_component + l2[n]; //DX : coord1-coord2+a*lattice(1) + (b*lattice(2))
            for(uint p=0;p<l3.size();p++){
              tmp = ab_component + l3[p];                       //DX : coord1-coord2+a*lattice(1) + (b*lattice(2)) + (c*lattice(3))
              min_dist=aurostd::min(min_dist,aurostd::modulus(tmp));
#if DEBUG_MINIMUM_DISTANCE
              if(LDEBUG) {cerr << __AFLOW_FUNC__ << " atoms[" << i << "," << k << "]: a=" << a_index[m] << ",b=" << b_index[n] << ",c=" << c_index[p] << "; min_dist=" << min_dist << "; this_dist=" << aurostd::modulus(tmp) << endl;}
#endif
            }
          }
        }
      }
    }


    //for(uint i=0; i<atoms.size(); i++){
    //  // Cannot reduce more than (1,1,1), so don't recalculate
    //  if(!(dims[1]==1 && dims[2]==1 && dims[3]==1)){
    //    dims=LatticeDimensionSphere(lattice,min_dist);
    //    //cerr << "dims: " << dims << endl;
    //  }
    //  for(uint k=i+1; k<atoms.size(); k++){
    //    for(int a=-dims[1];a<=dims[1];a++){
    //      for(int b=-dims[2];b<=dims[2];b++){
    //        for(int c=-dims[3];c<=dims[3];c++){
    //          min_dist=aurostd::min(min_dist,aurostd::modulus(atoms.at(k).cpos-atoms.at(i).cpos+a*lattice(1)+b*lattice(2)+c*lattice(3)));
    //          if(LDEBUG) {cerr << __AFLOW_FUNC__ << " atoms[" << i << "," << k << "]: a=" << a << ",b=" << b << ",c=" << c << "; min_dist=" << min_dist << endl;}
    //        }
    //      }
    //    }
    //  }
    //}
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " DIST BETWEEN ATOMS [4] " << min_dist << endl;

    //rescale
    min_dist*=scale;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " RESCALED [5] " << min_dist << endl;

    //if(min_dist<_XPROTO_TOO_CLOSE_ERROR_){
    //  stringstream message;
    //  message << "Atoms appear to be overlapping (min_dist=" << min_dist << "<" << _XPROTO_TOO_CLOSE_ERROR_ << ")";
    //  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_);
    //}

    return min_dist;
  }
} // namespace SYM

// **********************************************************************************************************************
// defaultTolerance
// **********************************************************************************************************************
// Set the default tolerance; based on minimum interatomic distance
namespace SYM { 
  double defaultTolerance(const xstructure& xstr){
    double min_dist = xstr.dist_nn_min; //CO20180409
    if(min_dist == AUROSTD_NAN || min_dist == AUROSTD_MAX_DOUBLE){min_dist=SYM::minimumDistance(xstr);} //CO20180409 //DX20210615 - add AUROSTD_MAX_DOUBLE
    //if(xstr.dist_nn_min == AUROSTD_NAN){xstr.MinDist();}  //CO20180409
    //min_dist = xstr.dist_nn_min; //CO20180409
    double tolerance = min_dist/100.0;
    //double tolerance = min_dist/10.0;
    return tolerance;
  }
} //namespace SYM

// **********************************************************************************************************************
// checkAngle
// **********************************************************************************************************************
// Tolerance for angles; converted to a distance
namespace SYM {
  bool checkAngle(xvector<double>& v1, xvector<double>& v2, double input_angle, double tolerance){ //CO20190520 - removed pointers for bools and doubles, added const where possible
    bool is_deg = false;
    return checkAngle(v1, v2, input_angle, is_deg, tolerance);
  }
} // namespace SYM

namespace SYM {
  bool checkAngle(xvector<double>& v1, xvector<double>& v2, double input_angle, bool& is_deg, double tolerance){ //CO20190520 - removed pointers for bools and doubles, added const where possible
    double mod_v1 = aurostd::modulus(v1);
    double mod_v2 = aurostd::modulus(v2);
    double avg_vec_mod = (mod_v1+mod_v2)/(2.0);
    double angle_diff = aurostd::angle(v1,v2)-input_angle;
    if(is_deg){
      angle_diff *= Pi_r/180.0;
    }
    return (aurostd::abs(avg_vec_mod*sin(angle_diff)) < tolerance);
  }
} // namespace SYM

namespace SYM {
  bool checkAngle(double& mod_v1, double& mod_v2, double angle1, double angle2, double tolerance){ //CO20190520 - removed pointers for bools and doubles, added const where possible
    bool is_deg = false;
    return checkAngle(mod_v1, mod_v2, angle1, angle2, is_deg, tolerance);
  }
} // namespace SYM

namespace SYM {
  bool checkAngle(double& mod_v1, double& mod_v2, double angle1, double angle2, bool& is_deg, double tolerance){ //CO20190520 - removed pointers for bools and doubles, added const where possible
    double avg_vec_mod = (mod_v1+mod_v2)/(2.0);
    double angle_diff = angle1-angle2;
    if(is_deg){
      angle_diff *= Pi_r/180.0;
    }
    return (aurostd::abs(avg_vec_mod*sin(angle_diff)) < tolerance);
  }
} // namespace SYM

// **********************************************************************************************************************
// change_tolerance
// **********************************************************************************************************************
// Tolerance scan for symmetry analysis
namespace SYM {
  //DX20170905 [OBSOLETE] bool change_tolerance(xstructure& xstr, double tolerance, double& orig_tolerance_old, int& count , double& min_dist, bool& no_scan) //CO20190520 - removed pointers for bools and doubles, added const where possible
  bool change_tolerance(xstructure& xstr, double& tolerance, double& min_dist, bool& no_scan) //CO20190520 - removed pointers for bools and doubles, added const where possible,  //DX20190524 - need pointer for change tolerance, otherwise it will not update
  { //CO20200106 - patching for auto-indenting
    // Scans between 5*orig_tol/10.0 on upper end and 5*orig_tol/10 on the right
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    //DX20180526 [OBSOLETE] string directory=aurostd::execute2string("pwd"); //DX20180426 - added current working directory
    string directory=xstr.directory; //DX20180426 - added current working directory
    //DX20180222 [OBSOLETE] static unsigned int count = 0;
    uint cycle_count = 0; // keeps track of number of full cycle iteration
    double step = 0.05; // logarithmic step size
    double max_range = 1.0; // scan an order of magnitude in each direction
    double sign = 1.0;
    double orig_range = 0.0;
    double range = 0.0;
    double change_count_max = 41; //DX20210406

    double incomming_tolerance = tolerance; //store current tolerance
    double orig_tolerance = 0; 

    // Calculate original tolerance based on the number of times this function has been called
    cycle_count = (xstr.sym_eps_change_count+1)/2; //DX20180222 - is now system specific (changed count to xstr.sym_eps_change_count)
    if(cycle_count<=20){
      orig_range = step*cycle_count;
      if(xstr.sym_eps_change_count%2==0){  //if even, scanned down in previous iteration //DX20180222 - is now system specific (changed count to xstr.sym_eps_change_count)
        sign = -1.0;
      }
      else {  //if odd, scanned up in previous iteration
        sign = 1.0;
      }   
      orig_tolerance = tolerance / (std::pow(10.0,(sign*orig_range)));
    }
    else {
      orig_tolerance = tolerance;
    }
    if(orig_range<=max_range && cycle_count<=20){ 
      // Calculate next tolerance value
      xstr.sym_eps_change_count += 1; //DX20180222 - is now system specific (changed count to xstr.sym_eps_change_count)
      cycle_count = (xstr.sym_eps_change_count+1)/2; //DX20180222 - is now system specific (changed count to xstr.sym_eps_change_count)
      range = step*((double)cycle_count);
      if(xstr.sym_eps_change_count%2==0){  //if even, scan down //DX20180222 - is now system specific (changed count to xstr.sym_eps_change_count)
        sign = -1.0;
      }   
      else {  //if odd, scan up
        sign = 1.0;
      }    
      if(!no_scan && xstr.sym_eps_change_count<=change_count_max){ //DX20210330 - check cycle count
        if(range<=max_range){
          tolerance = std::pow(10.0,(std::log10(orig_tolerance)+(sign*range)));
          if(tolerance >= min_dist){ //if larger than min distance, force lower scan
            sign = -1.0;
            tolerance = std::pow(10.0,(std::log10(orig_tolerance)+(sign*range)));
            xstr.sym_eps_change_count += 1; // can no longer scan up, so we increase the count //DX20180222 - is now system specific (changed count to xstr.sym_eps_change_count)
          }
          if(tolerance < _ZERO_TOL_){
            sign = 1.0;
            xstr.sym_eps_change_count += 1; // can no longer scan down, so we increase the count //DX20180222 - is now system specific (changed count to xstr.sym_eps_change_count)
            cycle_count = (xstr.sym_eps_change_count+1)/2; //DX20180222 - is now system specific (changed count to xstr.sym_eps_change_count)
            range = step*((double)cycle_count);
            tolerance = std::pow(10.0,(std::log10(orig_tolerance)+(sign*range)));
          }
          if(LDEBUG) {
            cerr << "SYM::change_tolerance: Inconsistent symmetry at this tolerance " << incomming_tolerance << ", changing tolerance to: " << tolerance << " [str.title=" << xstr.title << "] [dir=" << directory << "]." << endl; //DX20220926
          } 
          return TRUE; 
        }
        else {
          cerr << "SYM::change_tolerance WARNING: Inconsistent symmetry, tolerance range (" << std::pow(10.0,(std::log10(orig_tolerance)-(max_range))) << " to " << std::pow(10.0,(std::log10(orig_tolerance)+(max_range))) << ") tested [str.title=" << xstr.title << "] [dir=" << directory << "]." << endl; //DX20180426 - changed xstr.directory to directory (pwd)
          //count -=1;
          xstr.sym_eps_change_count += 1; //DX20210330 - need to increase the count here too
          no_scan = true;
          tolerance = orig_tolerance; //DX20170906
          return FALSE;
        }
      }
      //count -=1;
      no_scan = true;
      cerr << "SYM::change_tolerance WARNING: Inconsistent symmetry, but the tolerance scan is suppressed. [str.title=" << xstr.title << "] [dir=" << directory << "]" << endl; //DX20180426 - changed xstr.directory to directory (pwd)
      tolerance = orig_tolerance; //DX20170906
      return FALSE;
    }
    //return TRUE;
    no_scan = true;
    tolerance = orig_tolerance; //DX20170906
    //cerr << "SYM::change_tolerance WARNING: Inconsistent symmetry, but the tolerance scan is suppressed. [dir=" << xstr.directory << "]" << endl;
    return FALSE;
  }
} // namespace SYM

// **********************************************************************************************************************
// break_up_by_type
// **********************************************************************************************************************
// Break up a vector of atoms into types
namespace SYM{
  vector<vector<_atom> > break_up_by_type(vector<_atom> expanded_crystal){
    vector<int> used_types;
    vector<_atom> temp_atom_group;
    vector<vector<_atom > > out;
    uint count=0;
    bool NEW = true;
    int type;
    while(count < expanded_crystal.size()){
      for(uint i=0;i<expanded_crystal.size();i++){
        NEW = true;
        for(uint j=0;j<used_types.size();j++){
          if(expanded_crystal[i].type == used_types[j]){
            NEW = false;
          }
        }
        if(NEW == true){
          type = expanded_crystal[i].type;
          used_types.push_back(type);
          for(uint ii=0;ii<expanded_crystal.size();ii++){
            if(expanded_crystal[ii].type == type){
              temp_atom_group.push_back(expanded_crystal[ii]);
              count++;
            }
          }
          out.push_back(temp_atom_group);
          temp_atom_group.clear();
        }
      }
    }
    return out;
  }
}

namespace SYM {
  deque<deque<_atom> > break_up_by_type(deque<_atom>& expanded_crystal){
    vector<int> used_types;
    deque<_atom> temp_atom_group;
    deque<deque<_atom > > out;
    uint count=0;
    bool NEW = true;
    int type;
    while(count < expanded_crystal.size()){
      for(uint i=0;i<expanded_crystal.size();i++){
        NEW = true;
        for(uint j=0;j<used_types.size();j++){
          if(expanded_crystal[i].type == used_types[j]){
            NEW = false;
          }
        }
        if(NEW == true){
          type = expanded_crystal[i].type;
          used_types.push_back(type);
          for(uint ii=0;ii<expanded_crystal.size();ii++){
            if(expanded_crystal[ii].type == type){
              temp_atom_group.push_back(expanded_crystal[ii]);
              count++;
            }
          }
          out.push_back(temp_atom_group);
          temp_atom_group.clear();
        }
      }
    }
    return out;
  }
} // namespace SYM

// ******************************************************************************
// AtomsMapped
// ******************************************************************************
// Determine if atoms fractional coordinates are the same (also considers cases 
// where atoms are close to the edge of the unit cell wall). Requires atoms to 
// be the same type
namespace SYM {
  //overload //DX20190620
  bool AtomsMapped(const _atom& a, const _atom& b, const xmatrix<double>& lattice, bool skew, double tol){
    xmatrix<double> f2c=trasp(lattice);
    return AtomsMapped(a, b, lattice, f2c, skew, tol);
  }

  bool AtomsMapped(const _atom& a, const _atom& b, const xmatrix<double>& lattice, const xmatrix<double>& f2c, bool skew, double tol){ //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
    if(a.spin_is_given){ //DX20170921 - magnetic symmetry
      bool spins_match = (aurostd::abs(a.spin-b.spin)<=tol);
      //return a.type==b.type&&spins_match&&FPOSMatch(a,b,lattice,f2c,skew,tol); //DX20190619 - lattice and f2c as input, remove "Atom" prefix from name
      return a.type==b.type&&(aurostd::isequal(a.partial_occupation_value,b.partial_occupation_value,_AFLOW_POCC_ZERO_TOL_))&&spins_match&&FPOSMatch(a,b,lattice,f2c,skew,tol); //DX20190619 - lattice and f2c as input, remove "Atom" prefix from name
    }
    else if(a.noncoll_spin_is_given){ //DX20171205 - magnetic sym (non-collinear)
      bool spins_match = (aurostd::abs(a.noncoll_spin(1)-b.noncoll_spin(1))<=tol && aurostd::abs(a.noncoll_spin(2)-b.noncoll_spin(2))<=tol &&
          aurostd::abs(a.noncoll_spin(3)-b.noncoll_spin(3))<=tol);
      //return a.type==b.type&&spins_match&&FPOSMatch(a,b,lattice,f2c,skew,tol); //DX20190619 - lattice and f2c as input, remove "Atom" prefix from name
      return a.type==b.type&&(aurostd::isequal(a.partial_occupation_value,b.partial_occupation_value,_AFLOW_POCC_ZERO_TOL_))&&spins_match&&FPOSMatch(a,b,lattice,f2c,skew,tol); //DX20190619 - lattice and f2c as input, remove "Atom" prefix from name
    }
    else {
      //return FPOSMatch(a,b,lattice,f2c,skew,tol)&&a.type==b.type; //DX20190619 - lattice and f2c as input, remove "Atom" prefix from name 
      return FPOSMatch(a,b,lattice,f2c,skew,tol)&&a.type==b.type&&(aurostd::isequal(a.partial_occupation_value,b.partial_occupation_value,_AFLOW_POCC_ZERO_TOL_)); //DX20190619 - lattice and f2c as input, remove "Atom" prefix from name 
    }
  }
} // namespace SYM

// **********************************************************************************************************************
// Is lattice Skewed?
// **********************************************************************************************************************
// If the lattice too skewed to used the PBC function? 
namespace SYM {
  bool isLatticeSkewed(const xmatrix<double>& lattice, double& min_dist, double tol){ //CO20190520 - removed pointers for bools and doubles, added const where possible
    double max_skew = aurostd::max(aurostd::abs(aurostd::scalar_product(lattice(1),lattice(2))/(aurostd::modulus(lattice(1))*aurostd::modulus(lattice(2)))),
        aurostd::max(aurostd::abs(aurostd::scalar_product(lattice(2),lattice(3))/(aurostd::modulus(lattice(2))*aurostd::modulus(lattice(3)))),
          aurostd::abs(aurostd::scalar_product(lattice(3),lattice(1)/(aurostd::modulus(lattice(3))*aurostd::modulus(lattice(1))))))); //DX20200724 - SYM::DotPro to aurostd::scalar_product
    double skew_tol = (1.0-aurostd::abs(max_skew))*min_dist;
    XHOST.SKEW_TOL=skew_tol; //DX20171019
    return (skew_tol<tol);
  }
} // namespace SYM

//DX20190612 START

// ******************************************************************************
// minimize distance - Cartesian method (ROBUST)
// ******************************************************************************
// Calculates the minimum distance between a pair of Cartesian coordinates. 
// Here, we consider all images of the second coordinate and find the minimum distance
// with respect to the first coordinate. 
// This is the alternative to the "fractional method" (below). 
// If the tolerance is larger than the skewed tolerance, then fractional method 
// may give incorrect minimum distances and this function must be used instead.
// Note: This function takes a bit longer than fractional method, but it is robust.

namespace SYM {
  xvector<double> minimizeDistanceCartesianMethod(const xvector<double>& cpos1, const xvector<double>& cpos2, const xmatrix<double>& lattice){
    xvector<int> ijk;
    return minimizeDistanceCartesianMethod(cpos1,cpos2,lattice,ijk);
  }
}

namespace SYM {
  xvector<double> minimizeDistanceCartesianMethod(const xvector<double>& cpos1, const xvector<double>& cpos2, const xmatrix<double>& lattice, xvector<int>& ijk){

    // ---------------------------------------------------------------------------
    // calculate original distance between points and make temporary variables //DX20210315 - moved to top
    xvector<double> incell_dist = cpos1-cpos2;
    xvector<double> min_vec = incell_dist;
    double min_mod = aurostd::modulus(min_vec);

    if(min_mod<_ZERO_TOL_){ return min_vec; } //DX20210316 - return right away

    // ---------------------------------------------------------------------------
    // determine dimensions to produce uniform sphere; account for skewness 
    double radius=min_mod; //DX20210315 - use distance between atoms, significant speed up for skewed cells (used to be RadiusSphereLattice(lattice); slow)
    xvector<int> dims(3);

    // ---------------------------------------------------------------------------
    // multiply lattice vectors by index once and store - faster
    // DX20210315 - used resetLatticeDimensions, its faster, and it starts with
    // the zeroth (center) cell
    vector<xvector<double> > l1, l2, l3;
    vector<int> a_index, b_index, c_index;
    resetLatticeDimensions(lattice,radius,dims,l1,l2,l3,a_index,b_index,c_index);

    xvector<double> a_component, ab_component, tmp; //DX20210315 - define outside loop
    double mod_tmp = AUROSTD_MAX_DOUBLE; //DX20210315 - initialize

    // ---------------------------------------------------------------------------
    // consider other periodic images within the grid-sphere 
    // vector stored in each loop to save computation; fewer duplicate operations
    for(uint i=0;i<l1.size();i++){
      a_component = incell_dist + l1[i];    //DX : cpos1-cpos2+a*lattice(1)
      for(uint j=0;j<l2.size();j++){
        ab_component = a_component + l2[j]; //DX : cpos1-cpos2+a*lattice(1) + (b+lattice(2))
        for(uint k=0;k<l3.size();k++){
          tmp = ab_component + l3[k];       //DX : cpos1-cpos2+a*lattice(1) + (b+lattice(2)) + (c+lattice(3))
          mod_tmp = aurostd::modulus(tmp);
          if(mod_tmp<min_mod){
            min_mod = mod_tmp;
            min_vec = tmp;
            ijk(1) = a_index[i];
            ijk(2) = b_index[j];
            ijk(3) = c_index[k];
            if(min_mod<_ZERO_TOL_){ return min_vec; } //DX20210315 - return right away
            // ---------------------------------------------------------------------------
            // diminishing dims: if minimum distance changed, then we may not need to
            // search as far; reset loop and search again based on new minimum distance //DX20210315 (significant speed up)
            if(!(dims[1]==1 && dims[2]==1 && dims[3]==1)){
              resetLatticeDimensions(lattice,min_mod,dims,l1,l2,l3,a_index,b_index,c_index);
              i=j=k=0;
            }
          } 
        }
      }
    }

    //DX20180423 - FASTER MIN CART DISTANCE CALCULATOR - END
    //DX20180423 [OBSOLETE] min_vec = coord1-coord2;
    //DX20180423 [OBSOLETE] double min_mod = aurostd::modulus(min_vec);
    //DX20180423 [OBSOLETE] for(int a=-dims[1];a<=dims[1];a++){
    //DX20180423 [OBSOLETE]   for(int b=-dims[2];b<=dims[2];b++){
    //DX20180423 [OBSOLETE]     for(int c=-dims[3];c<=dims[3];c++){
    //DX20180423 [OBSOLETE]       tmp = coord1-coord2+a*lattice(1)+b*lattice(2)+c*lattice(3);
    //DX20180423 [OBSOLETE]       mod_tmp = aurostd::modulus(tmp);
    //DX20180423 [OBSOLETE]       if(mod_tmp<min_mod){
    //DX20180423 [OBSOLETE]         min_mod = mod_tmp;
    //DX20180423 [OBSOLETE]         min_vec = tmp;
    //DX20180423 [OBSOLETE]         ijk(1) = a;
    //DX20180423 [OBSOLETE]         ijk(2) = b;
    //DX20180423 [OBSOLETE]         ijk(3) = c;
    //DX20180423 [OBSOLETE]       } 
    //DX20180423 [OBSOLETE]     }
    //DX20180423 [OBSOLETE]   }
    //DX20180423 [OBSOLETE] }
    //DX20180903 needs to return min_vec : return min_mod;

    return min_vec; //DX20180903 needs to return min_vec
  }
} //namespace SYM

// ******************************************************************************
// minimize distance - fractional method (FAST)
// ******************************************************************************
// Minimizes fractional space distance vector component by component.  
// This may not give the the true minimum distance when transforming back to a 
// E and take the modulus of minimum cpos distanceuclidean space (Cartesian). 
// HOWEVER, this method is safe to use if the mapping tolerance is small.  
// A heuristic maximum tolerance is the suprenum of 
// the normalized off-diagonal elements of the metric tensor (this is a well-tested quanitity).
// For more information, contact David Hicks (d.hicks@duke.edu).

// ---------------------------------------------------------------------------
// overloads 
namespace SYM {
  xvector<double> minimizeDistanceFractionalMethod(const xvector<double>& fpos1, const xvector<double>& fpos2){
    xvector<int> ijk;
    xvector<double> fdiff = fpos1-fpos2;
    return minimizeDistanceFractionalMethod(fdiff, ijk);
  }
} /// namespace SYM

namespace SYM {
  xvector<double> minimizeDistanceFractionalMethod(const xvector<double>& fdiff){
    //DX20190904 [OBSOLETE] xvector<int> ijk;
    //DX20190904 [OBSOLETE] return minimizeDistanceFractionalMethod(fdiff, ijk);

    // ---------------------------------------------------------------------------
    // loop over each component and bring close to origin 
    // (i.e., between -0.5 and 0.5 in fractional coordinates)

    //DX20190904 START
    double tolerance=0.0;
    double upper_bound=0.5;
    double lower_bound=-0.5;

    xvector<double> min_diff = fdiff;
    BringInCellInPlace(min_diff, tolerance, upper_bound, lower_bound);

    return min_diff;
    //DX20190904 END
  }
} /// namespace SYM

// ---------------------------------------------------------------------------
namespace SYM {
  xvector<double> minimizeDistanceFractionalMethod(const xvector<double>& fdiff, xvector<int>& ijk){

    // ---------------------------------------------------------------------------
    // loop over each component and bring close to origin 
    // (i.e., between -0.5 and 0.5 in fractional coordinates)

    xvector<double> min_diff = minimizeDistanceFractionalMethod(fdiff); //DX20190904

    for(int i=min_diff.lrows; i<=min_diff.urows; i++){ //DX20190904
      ijk(i) = (int)(min_diff[i]-fdiff[i]);
    }

    //DX20190416 [OBSOLETE] xvector<double> tmp = v_in;
    //DX20190904 [OBSOLETE] for(uint i=1; i<4; i++){
    //DX20190904 [OBSOLETE]   while(min_diff(i)>0.5){
    //DX20190904 [OBSOLETE]     min_diff(i) = min_diff(i)-1;
    //DX20190904 [OBSOLETE]     ijk(i)-=1;
    //DX20190904 [OBSOLETE]   }
    //DX20190904 [OBSOLETE]   while(min_diff(i)<=-0.5){
    //DX20190904 [OBSOLETE]     min_diff(i) = min_diff(i)+1;
    //DX20190904 [OBSOLETE]     ijk(i)+=1;
    //DX20190904 [OBSOLETE]   }
    //DX20190904 [OBSOLETE] }
    return min_diff;
  }
}

//DX20190612 END
//DX20190613 [OBSOLETE] // ******************************************************************************
//DX20190613 [OBSOLETE] // Minimize Cartesian Distance 
//DX20190613 [OBSOLETE] // ******************************************************************************
//DX20190613 [OBSOLETE] // Calculates the minimum distance between a pair of Cartesian coordinates. 
//DX20190613 [OBSOLETE] // Here, we consider all images of the second coordinate and find the minimum distance
//DX20190613 [OBSOLETE] // with respect to the first coordinate. 
//DX20190613 [OBSOLETE] // This is the alternative to the PBC function (below). If the tolerance is larger than 
//DX20190613 [OBSOLETE] // the skewed tolerance, then PBC is not robust enough and this function is used instead.
//DX20190613 [OBSOLETE] // Note: This function takes a bit longer than PBC, but we do not loose robustness.
//DX20190613 [OBSOLETE] namespace SYM {
//DX20190613 [OBSOLETE]   bool minimizeCartesianDistance(const xvector<double>& coord1, const xvector<double>& coord2, xvector<double>& out, 
//DX20190613 [OBSOLETE]                                  const xmatrix<double>& c2f, const xmatrix<double>& f2c, double tol){ //CO20190520 - removed pointers for bools and doubles, added const where possible
//DX20190613 [OBSOLETE]     xvector<int> ijk;
//DX20190613 [OBSOLETE]     bool restriction = false;
//DX20190613 [OBSOLETE]     return minimizeCartesianDistance(coord1, coord2, out, c2f, f2c, ijk, restriction, tol);
//DX20190613 [OBSOLETE]   }
//DX20190613 [OBSOLETE] } /// namespace SYM
//DX20190613 [OBSOLETE] 
//DX20190613 [OBSOLETE] namespace SYM {
//DX20190613 [OBSOLETE]   bool minimizeCartesianDistance(const xvector<double>& coord1, const xvector<double>& coord2, xvector<double>& out, 
//DX20190613 [OBSOLETE]                                  const xmatrix<double>& c2f, const xmatrix<double>& f2c, xvector<int>& ijk, bool& restriction, double tol){ //CO20190520 - removed pointers for bools and doubles, added const where possible
//DX20190613 [OBSOLETE]     xmatrix<double> lattice = trasp(f2c);
//DX20190613 [OBSOLETE]     xvector<double> min_vec;
//DX20190613 [OBSOLETE]     double min_mod=minimumCartesianDistance(coord1,coord2,lattice,min_vec,ijk);
//DX20190613 [OBSOLETE]     if(restriction){
//DX20190613 [OBSOLETE]       if(aurostd::abs(ijk(1))>3 || aurostd::abs(ijk(2))>3 || aurostd::abs(ijk(3))>3){
//DX20190613 [OBSOLETE]         return false;
//DX20190613 [OBSOLETE]       }
//DX20190613 [OBSOLETE]     }
//DX20190613 [OBSOLETE]     out = c2f*min_vec;
//DX20190613 [OBSOLETE]     return (min_mod<tol);
//DX20190613 [OBSOLETE]   }
//DX20190613 [OBSOLETE] } //namespace SYM
//DX20190613 [OBSOLETE] 
//DX20190613 [OBSOLETE] namespace SYM {
//DX20190613 [OBSOLETE] // ******************************************************************************
//DX20190613 [OBSOLETE] // minimumCartesianDistance()
//DX20190613 [OBSOLETE] // ******************************************************************************
//DX20190613 [OBSOLETE] // Finds minimum distance between two cartesian coordinates in periodic cell
//DX20190613 [OBSOLETE] // min_vec is full cartesian vector distances
//DX20190613 [OBSOLETE] // ijk are the indices of the cell of the minimum distance image
//DX20190613 [OBSOLETE] //ME20180730: minimumCartesianVector returns the smallest vector between two atoms.
//DX20190613 [OBSOLETE] double minimumCartesianDistance(const xvector<double>& coord1, const xvector<double>& coord2, const xmatrix<double>& lattice){
//DX20190613 [OBSOLETE]     xvector<double> min_vec;
//DX20190613 [OBSOLETE]     xvector<int> ijk;
//DX20190613 [OBSOLETE]     return minimumCartesianDistance(coord1,coord2,lattice,min_vec,ijk);
//DX20190613 [OBSOLETE] }
//DX20190613 [OBSOLETE] 
//DX20190613 [OBSOLETE] double minimumCartesianDistance(const xvector<double>& coord1, const xvector<double>& coord2, const xmatrix<double>& lattice, 
//DX20190613 [OBSOLETE]     xvector<double>& min_vec,xvector<int>& ijk){
//DX20190613 [OBSOLETE]    min_vec = minimumCartesianVector(coord1, coord2, lattice, ijk);
//DX20190613 [OBSOLETE]    return aurostd::modulus(min_vec);
//DX20190613 [OBSOLETE] }
//DX20190613 [OBSOLETE]     
//DX20190613 [OBSOLETE] xvector<double> minimumCartesianVector(const xvector<double>& coord1, const xvector<double>& coord2,
//DX20190613 [OBSOLETE]                                        const xmatrix<double>& lattice) {
//DX20190613 [OBSOLETE]   xvector<int> ijk;
//DX20190613 [OBSOLETE]   xvector<double> min_vec = minimumCartesianVector(coord1, coord2, lattice, ijk);
//DX20190613 [OBSOLETE]   return min_vec;
//DX20190613 [OBSOLETE] }
//DX20190613 [OBSOLETE] 
//DX20190613 [OBSOLETE] xvector<double> minimumCartesianVector(const xvector<double>& coord1, const xvector<double>& coord2,
//DX20190613 [OBSOLETE]                                        const xmatrix<double>& lattice, xvector<int>& ijk) {
//DX20190613 [OBSOLETE]     xvector<double> tmp, min_vec;
//DX20190613 [OBSOLETE]     double mod_tmp;
//DX20190613 [OBSOLETE] 
//DX20190613 [OBSOLETE]     double radius=RadiusSphereLattice(lattice);
//DX20190613 [OBSOLETE]     xvector<int> dims(3);
//DX20190613 [OBSOLETE]     dims=LatticeDimensionSphere(lattice,radius);
//DX20190613 [OBSOLETE] 
//DX20190613 [OBSOLETE]     //DX20180423 - FASTER MIN CART DISTANCE CALCULATOR - START
//DX20190613 [OBSOLETE]     //DX20180423 - only calculate multiplication once (time-saver)
//DX20190613 [OBSOLETE]     vector<xvector<double> > l1, l2, l3;
//DX20190613 [OBSOLETE]     vector<int> a_index, b_index, c_index;
//DX20190613 [OBSOLETE]     for(int a=-dims[1];a<=dims[1];a++){l1.push_back(a*lattice(1));a_index.push_back(a);} //DX calc once and store
//DX20190613 [OBSOLETE]     for(int b=-dims[2];b<=dims[2];b++){l2.push_back(b*lattice(2));b_index.push_back(b);} //DX calc once and store
//DX20190613 [OBSOLETE]     for(int c=-dims[3];c<=dims[3];c++){l3.push_back(c*lattice(3));c_index.push_back(c);} //DX calc once and store
//DX20190613 [OBSOLETE] 
//DX20190613 [OBSOLETE]     xvector<double> incell_dist = coord1-coord2;
//DX20190613 [OBSOLETE] 
//DX20190613 [OBSOLETE]     min_vec = coord1-coord2;
//DX20190613 [OBSOLETE]     double min_mod = aurostd::modulus(min_vec);
//DX20190613 [OBSOLETE] 
//DX20190613 [OBSOLETE]     //DX20180423 - running vector in each loop saves computations; fewer duplicate operations
//DX20190613 [OBSOLETE]     for(uint i=0;i<l1.size();i++){
//DX20190613 [OBSOLETE]       xvector<double> a_component = incell_dist + l1[i];    //DX : coord1-coord2+a*lattice(1)
//DX20190613 [OBSOLETE]       for(uint j=0;j<l2.size();j++){
//DX20190613 [OBSOLETE]         xvector<double> ab_component = a_component + l2[j]; //DX : coord1-coord2+a*lattice(1) + (b+lattice(2))
//DX20190613 [OBSOLETE]         for(uint k=0;k<l3.size();k++){
//DX20190613 [OBSOLETE]           tmp = ab_component + l3[k];                       //DX : coord1-coord2+a*lattice(1) + (b+lattice(2)) + (c+lattice(3))
//DX20190613 [OBSOLETE]           mod_tmp = aurostd::modulus(tmp);
//DX20190613 [OBSOLETE]           if(mod_tmp<min_mod){
//DX20190613 [OBSOLETE]             min_mod = mod_tmp;
//DX20190613 [OBSOLETE]             min_vec = tmp;
//DX20190613 [OBSOLETE]             ijk(1) = a_index[i];
//DX20190613 [OBSOLETE]             ijk(2) = b_index[j];
//DX20190613 [OBSOLETE]             ijk(3) = c_index[k];
//DX20190613 [OBSOLETE]           } 
//DX20190613 [OBSOLETE]         }
//DX20190613 [OBSOLETE]       }
//DX20190613 [OBSOLETE]     }
//DX20190613 [OBSOLETE]     //DX20180423 - FASTER MIN CART DISTANCE CALCULATOR - END
//DX20190613 [OBSOLETE] 
//DX20190613 [OBSOLETE]     //DX20180423 [OBSOLETE] min_vec = coord1-coord2;
//DX20190613 [OBSOLETE]     //DX20180423 [OBSOLETE] double min_mod = aurostd::modulus(min_vec);
//DX20190613 [OBSOLETE]     //DX20180423 [OBSOLETE] for(int a=-dims[1];a<=dims[1];a++){
//DX20190613 [OBSOLETE]     //DX20180423 [OBSOLETE]   for(int b=-dims[2];b<=dims[2];b++){
//DX20190613 [OBSOLETE]     //DX20180423 [OBSOLETE]     for(int c=-dims[3];c<=dims[3];c++){
//DX20190613 [OBSOLETE]     //DX20180423 [OBSOLETE]       tmp = coord1-coord2+a*lattice(1)+b*lattice(2)+c*lattice(3);
//DX20190613 [OBSOLETE]     //DX20180423 [OBSOLETE]       mod_tmp = aurostd::modulus(tmp);
//DX20190613 [OBSOLETE]     //DX20180423 [OBSOLETE]       if(mod_tmp<min_mod){
//DX20190613 [OBSOLETE]     //DX20180423 [OBSOLETE]         min_mod = mod_tmp;
//DX20190613 [OBSOLETE]     //DX20180423 [OBSOLETE]         min_vec = tmp;
//DX20190613 [OBSOLETE]     //DX20180423 [OBSOLETE]         ijk(1) = a;
//DX20190613 [OBSOLETE]     //DX20180423 [OBSOLETE]         ijk(2) = b;
//DX20190613 [OBSOLETE]     //DX20180423 [OBSOLETE]         ijk(3) = c;
//DX20190613 [OBSOLETE]     //DX20180423 [OBSOLETE]       } 
//DX20190613 [OBSOLETE]     //DX20180423 [OBSOLETE]     }
//DX20190613 [OBSOLETE]     //DX20180423 [OBSOLETE]   }
//DX20190613 [OBSOLETE]     //DX20180423 [OBSOLETE] }
//DX20190613 [OBSOLETE]     //DX20180903 needs to return min_vec : return min_mod;
//DX20190613 [OBSOLETE]     return min_vec; //DX20180903 needs to return min_vec
//DX20190613 [OBSOLETE]     }
//DX20190613 [OBSOLETE] } //namespace SYM
//DX20190613 [OBSOLETE] 
//DX20190613 [OBSOLETE] // ******************************************************************************
//DX20190613 [OBSOLETE] // Periodic Boundary Conditions (PBC)
//DX20190613 [OBSOLETE] // ******************************************************************************
//DX20190613 [OBSOLETE] // Minimizes fractional vector component by component.  This may not give the 
//DX20190613 [OBSOLETE] // the true minimum distance in Cartesian. HOWEVER, it is safe to assume if the
//DX20190613 [OBSOLETE] // mapping tolerance is small.  A heuristic maximum tolerance is the suprenum of 
//DX20190613 [OBSOLETE] // the off-diagonal elements of the metric tensor (this is a well-tested quanitity).
//DX20190613 [OBSOLETE] // For more information, contact David Hicks (d.hicks@duke.edu).
//DX20190613 [OBSOLETE] namespace SYM {
//DX20190613 [OBSOLETE]   bool PBC(xvector<double>& v_in){
//DX20190613 [OBSOLETE]     xvector<int> ijk;
//DX20190613 [OBSOLETE]     bool restriction = false;
//DX20190613 [OBSOLETE]     return PBC(v_in, ijk, restriction);
//DX20190613 [OBSOLETE]   }
//DX20190613 [OBSOLETE] } /// namespace SYM
//DX20190613 [OBSOLETE] 
//DX20190613 [OBSOLETE] namespace SYM {
//DX20190613 [OBSOLETE]   bool PBC(xvector<double>& v_in, xvector<int>& ijk, bool& restriction){
//DX20190613 [OBSOLETE]     //DX20190416 [OBSOLETE] xvector<double> tmp = v_in;
//DX20190613 [OBSOLETE]     for(uint i=1; i<4; i++){
//DX20190613 [OBSOLETE]       while(v_in(i)>0.5){
//DX20190613 [OBSOLETE]         v_in(i) = v_in(i)-1;
//DX20190613 [OBSOLETE]         ijk(i)-=1;
//DX20190613 [OBSOLETE]       }
//DX20190613 [OBSOLETE]       while(v_in(i)<=-0.5){
//DX20190613 [OBSOLETE]         v_in(i) = v_in(i)+1;
//DX20190613 [OBSOLETE]         ijk(i)+=1;
//DX20190613 [OBSOLETE]       }
//DX20190613 [OBSOLETE]       if(aurostd::abs(ijk(i))>3 && restriction){
//DX20190613 [OBSOLETE]         //DEBUGGER
//DX20190613 [OBSOLETE]         //cerr << "ijk(i): " << ijk(i) << endl;
//DX20190613 [OBSOLETE]         //cerr << "ERROR (PBC Function): More than three cells away!" << endl;
//DX20190613 [OBSOLETE]         return false;
//DX20190613 [OBSOLETE]       }
//DX20190613 [OBSOLETE]     }
//DX20190613 [OBSOLETE]     return true;
//DX20190613 [OBSOLETE]   }
//DX20190613 [OBSOLETE] } /// namespace SYM

// ******************************************************************************
// FPOSMatch
// ******************************************************************************
// Determine if atoms in fractional coordinates are the same (also considers cases
// where atoms are close to the edge of the unit cell wall)

// --------------------------------------------------------
// FPOSMatch: match element of deque<_atom> to _atom 

// overload //DX20190620 
namespace SYM {
  bool FPOSMatch(const deque<_atom>& atom_set, const _atom& atom2, uint& match_type, const xmatrix<double>& lattice, bool skew, double tol){
    xmatrix<double> f2c=trasp(lattice);
    return FPOSMatch(atom_set,atom2,match_type,lattice,f2c,skew,tol);
  }

  bool FPOSMatch(const deque<_atom>& atom_set, const _atom& atom2, uint& match_type, const xmatrix<double>& lattice, 
      const xmatrix<double>& f2c, bool skew, double tol){ //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190620 - lattice and f2c as input and remove "Atom" prefix in name
    for(uint i=0;i<atom_set.size();i++){
      if(FPOSMatch(atom_set[i],atom2,lattice,f2c,skew,tol)){ //DX20190620 - lattice and f2c as input and remove "Atom" prefix in name
        match_type = atom_set[i].type;
        return TRUE;
      }
    }
    return FALSE;
  }
} //namespace SYM

// --------------------------------------------------------
// FPOSMatch: match _atom to _atom 

// overload //DX20190620 
namespace SYM {
  bool FPOSMatch(const _atom& atom1, const _atom& atom2, const xmatrix<double>& lattice, bool skew, double tol){ 
    xmatrix<double> f2c=trasp(lattice);
    return FPOSMatch(atom1,atom2,lattice,f2c,skew,tol);
  }

  bool FPOSMatch(const _atom& atom1, const _atom& atom2, const xmatrix<double>& lattice, 
      const xmatrix<double>& f2c, bool skew, double tol){ //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190620 - lattice and f2c as input and remove "Atom" prefix in name
    return FPOSMatch(atom1.fpos,atom2.fpos,lattice,f2c,skew,tol); //DX20190619 - added lattice to input and remove "Atom" prefix in name
    //xvector<double> fdiff = atom1.fpos - atom2.fpos;
    //if(skew){
    //return minimizeCartesianDistance(atom1.cpos,atom2.cpos,fdiff,c2f,f2c,tol);
    //}
    //else {
    //PBC(fdiff);
    //}
    //return (aurostd::modulus(f2c*fdiff)<tol);
  }
} // namespace SYM

// --------------------------------------------------------
// FPOSMatch: match fpos to fpos
// overload //DX20190620
namespace SYM {
  bool FPOSMatch(const xvector<double>& fpos1, const xvector<double>& fpos2, const xmatrix<double>& lattice, bool skew, double tol){
    xmatrix<double> f2c=trasp(lattice);
    return FPOSMatch(fpos1,fpos2,lattice,f2c,skew,tol);
  }

  bool FPOSMatch(const xvector<double>& fpos1, const xvector<double>& fpos2, const xmatrix<double>& lattice,
      const xmatrix<double>& f2c, bool skew, double tol){ //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190620 - lattice and f2c as input and remove "Atom" prefix in name

    //DX20190613 [OBSOLETE] xvector<double> fdiff = fpos1 - fpos2;
    if(XHOST.SKEW_TEST){
      //DX20190613 [OBOSLETE] xvector<double> pbc_fdiff = fdiff;
      //DX20190613 [OBOSLETE] PBC(pbc_fdiff);
      xvector<double> min_fdiff = minimizeDistanceFractionalMethod(fpos1,fpos2); //DX20190613
      double min_fdiff_dist = aurostd::modulus(f2c*min_fdiff);
      xvector<double> cpos1 = f2c*fpos1;
      xvector<double> cpos2 = f2c*fpos2;
      //DX20190613 [OBOSLETE] minimizeCartesianDistance(cpos1,cpos2,fdiff,c2f,f2c,tol);
      xvector<double> min_cdiff = minimizeDistanceCartesianMethod(cpos1, cpos2, lattice); //DX20190613
      double min_cdiff_dist = aurostd::modulus(min_cdiff); //DX20190613 - changed variable names
      stringstream message; 
      if((min_cdiff_dist<=tol)==(min_fdiff_dist<=tol) && aurostd::abs(min_cdiff_dist-min_fdiff_dist)<_ZERO_TOL_){
        message << __AFLOW_FUNC__ << " minimum distances equal, and mappings same -- globally optimized: " << min_cdiff_dist << " | bring-in-cell: " << min_fdiff_dist << " || tol: " << tol << " || skew_tol: " << XHOST.SKEW_TOL;
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, cerr, _LOGGER_MESSAGE_);
      }
      if((min_cdiff_dist<=tol)==(min_fdiff_dist<=tol) && min_cdiff_dist<=tol && aurostd::abs(min_cdiff_dist-min_fdiff_dist)>_ZERO_TOL_){
        message << __AFLOW_FUNC__ << " WARNING-MAP: minimum distances unequal, but mapping outcome same -- globally optimized: " << min_cdiff_dist << " | bring-in-cell: " << min_fdiff_dist << " || tol: " << tol << " || skew_tol: " << XHOST.SKEW_TOL;
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, cerr, _LOGGER_WARNING_);
      }
      if((min_cdiff_dist<=tol)==(min_fdiff_dist<=tol) && min_cdiff_dist>tol && aurostd::abs(min_cdiff_dist-min_fdiff_dist)>_ZERO_TOL_){
        message << __AFLOW_FUNC__ << " WARNING-NOMAP: minimum distances unequal, but mapping outcome same -- globally optimized: " << min_cdiff_dist << " | bring-in-cell: " << min_fdiff_dist << " || tol: " << tol << " || skew_tol: " << XHOST.SKEW_TOL;
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, cerr, _LOGGER_WARNING_);
      }
      if((min_cdiff_dist<=tol)!=(min_fdiff_dist<=tol) && aurostd::abs(min_cdiff_dist-min_fdiff_dist)>_ZERO_TOL_){
        message << __AFLOW_FUNC__ << "ERROR: minimum distances unequal, and mappings unequal -- globally optimized: " << min_cdiff_dist << " | bring-in-cell: " << min_fdiff_dist << " || tol: " << tol << " || skew_tol: " << XHOST.SKEW_TOL;
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_GENERIC_ERROR_);
      }
      if((min_cdiff_dist-XHOST.SKEW_TOL<_ZERO_TOL_)!=(min_fdiff_dist-XHOST.SKEW_TOL<_ZERO_TOL_)){
        message << __AFLOW_FUNC__ << "THRESHOLD ERROR: minimum distances unequal, and mappings unequal -- globally optimized: " << min_cdiff_dist << " | bring-in-cell: " << min_fdiff_dist << " || tol: " << tol << " || skew_tol: " << XHOST.SKEW_TOL;
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_GENERIC_ERROR_);
      }
    }
    double min_dist = aurostd::modulus(CPOSDistFromFPOS(fpos1,fpos2,lattice,f2c,skew)); //DX20190620
    //DX20190620 [moved into CPOSDistFromFPOS()] if(skew){
    //DX20190620 [moved into CPOSDistFromFPOS()]   xvector<double> cpos1 = f2c*fpos1;
    //DX20190620 [moved into CPOSDistFromFPOS()]   xvector<double> cpos2 = f2c*fpos2;
    //DX20190620 [moved into CPOSDistFromFPOS()]   //DX20190613 [OBSOLETE] return minimizeCartesianDistance(cpos1,cpos2,fdiff,c2f,f2c,tol);
    //DX20190620 [moved into CPOSDistFromFPOS()]   xvector<double> min_cdiff = minimizeDistanceCartesianMethod(cpos1,cpos2,lattice); //DX20190613
    //DX20190620 [moved into CPOSDistFromFPOS()]   min_dist = aurostd::modulus(min_cdiff); //DX20190613
    //DX20190620 [moved into CPOSDistFromFPOS()] }
    //DX20190620 [moved into CPOSDistFromFPOS()] else {
    //DX20190620 [moved into CPOSDistFromFPOS()]   if(VERBOSE){ cerr << __AFLOW_FUNC__ << " fpos1-fpos2=" << (fpos1-fpos2) << endl;}
    //DX20190620 [moved into CPOSDistFromFPOS()]   //DX20190613 [OBSOLETE] PBC(fdiff);
    //DX20190620 [moved into CPOSDistFromFPOS()]   xvector<double> min_fdiff = minimizeDistanceFractionalMethod(fpos1,fpos2); //DX20190613
    //DX20190620 [moved into CPOSDistFromFPOS()]   if(VERBOSE){ cerr << __AFLOW_FUNC__ << " min_fdiff=" << min_fdiff << endl;}
    //DX20190620 [moved into CPOSDistFromFPOS()]   min_dist = aurostd::modulus(f2c*min_fdiff); //DX20190613
    //DX20190620 [moved into CPOSDistFromFPOS()]   if(VERBOSE){ cerr << __AFLOW_FUNC__ << " min_dist=" << min_dist << endl;}
    //DX20190620 [moved into CPOSDistFromFPOS()] }
    return (min_dist<tol); //DX20190613
  }
} // namespace SYM

// ******************************************************************************
// FPOSDistance
// ******************************************************************************
namespace SYM {
  xvector<double> FPOSDistFromFPOS(const xvector<double>& fpos1,const xvector<double>& fpos2,
      const xmatrix<double>& lattice,bool skew){
    xmatrix<double> f2c=trasp(lattice);xmatrix<double> c2f=inverse(trasp(lattice));
    return FPOSDistFromFPOS(fpos1,fpos2,lattice,c2f,f2c,skew);
  }
  xvector<double> FPOSDistFromFPOS(const xvector<double>& fpos1,const xvector<double>& fpos2,
      const xmatrix<double>& lattice,const xmatrix<double>& c2f,const xmatrix<double>& f2c,bool skew){  //CO20190525
    bool VERBOSE=FALSE; //DX20201210
    xvector<double> min_fdiff;
    if(skew){
      xvector<double> cpos1 = f2c*fpos1;
      xvector<double> cpos2 = f2c*fpos2;
      xvector<double> min_cdiff = minimizeDistanceCartesianMethod(cpos1,cpos2,lattice); //DX20190613
      min_fdiff = c2f*min_cdiff;
      if(VERBOSE){
        cerr << __AFLOW_FUNC__ << " fpos1=" << fpos1 << endl;
        cerr << __AFLOW_FUNC__ << " fpos2=" << fpos2 << endl;
      }
    }
    else {
      xvector<double> min_fdiff = minimizeDistanceFractionalMethod(fpos1,fpos2); //DX20190613
      if(VERBOSE){
        cerr << __AFLOW_FUNC__ << " fpos1=" << fpos1 << endl;
        cerr << __AFLOW_FUNC__ << " fpos2=" << fpos2 << endl;
        cerr << __AFLOW_FUNC__ << " fpos1-fpos2=" << (fpos1-fpos2) << endl;
        cerr << __AFLOW_FUNC__ << " min_fdiff=" << min_fdiff << endl;
      }
    }
    return min_fdiff;
  }
}

// ******************************************************************************
// CPOSDistFromFPOS
// ******************************************************************************
// minimize the difference between two FPOS positions and returns the minimum
// distance in Cartesian coordinates
// (why fpos input?: minimizing fpos is faster (minimizeDistanceFractionalMethod()); 
//  we want to promote speed-up by avoiding c2f*fpos operation twice in else-statement)
namespace SYM {
  xvector<double> CPOSDistFromFPOS(const xvector<double>& fpos1,const xvector<double>& fpos2,
      const xmatrix<double>& lattice,bool skew){
    xmatrix<double> f2c=trasp(lattice);
    return CPOSDistFromFPOS(fpos1,fpos2,lattice,f2c,skew);
  }
  xvector<double> CPOSDistFromFPOS(const xvector<double>& fpos1,const xvector<double>& fpos2,
      const xmatrix<double>& lattice,const xmatrix<double>& f2c,bool skew){  //CO20190525
    bool VERBOSE=FALSE; //using LDEBUG would pollute output
    xvector<double> min_cdiff;
    if(skew){
      xvector<double> cpos1 = f2c*fpos1;
      xvector<double> cpos2 = f2c*fpos2;
      min_cdiff = minimizeDistanceCartesianMethod(cpos1,cpos2,lattice); //DX20190613
      if(VERBOSE){ //DX20201210 - only one if-statement, otherwise expensive
        cerr << __AFLOW_FUNC__ << " fpos1=" << fpos1 << endl;
        cerr << __AFLOW_FUNC__ << " fpos2=" << fpos2 << endl;
        cerr << __AFLOW_FUNC__ << " cpos1-cpos2=" << (cpos1-cpos2) << endl;
        cerr << __AFLOW_FUNC__ << " min_cdiff=" << min_cdiff << endl;
      }
    }
    else {
      xvector<double> min_fdiff = minimizeDistanceFractionalMethod(fpos1,fpos2); //DX20190613
      min_cdiff = f2c*min_fdiff;
      if(VERBOSE){ //DX20201210 - only one if-statement, otherwise expensive
        cerr << __AFLOW_FUNC__ << " fpos1=" << fpos1 << endl;
        cerr << __AFLOW_FUNC__ << " fpos2=" << fpos2 << endl;
        cerr << __AFLOW_FUNC__ << " fpos1-fpos2=" << (fpos1-fpos2) << endl;
        cerr << __AFLOW_FUNC__ << " min_fdiff=" << min_fdiff << endl;
        cerr << __AFLOW_FUNC__ << " min_cdiff=" << min_cdiff << endl;
      }
    }
    return min_cdiff;
  }
}

// ******************************************************************************
// validateAtomPosition
// ******************************************************************************
// applies c2f and f2c to atom.fpos and atom.cpos to see if there's a correspondence
namespace SYM {
  bool validateAtomPosition(const _atom& atom,const xmatrix<double>& c2f,const xmatrix<double>& f2c,bool skew,double& _eps_){ //CO20190520 - removed pointers for bools and doubles, added const where possible

    return validateAtomPosition(atom.cpos,atom.fpos,c2f,f2c,skew,_eps_);
  }
  bool validateAtomPosition(const xvector<double>& cpos,const xvector<double>& fpos,const xmatrix<double>& c2f,const xmatrix<double>& f2c,bool skew,double& _eps_){ //CO20190520 - removed pointers for bools and doubles, added const where possible
    //DX20190613 [OBSOLETE] xvector<double> fdiff = (c2f*cpos)-fpos;
    double min_dist = 1e9;
    if(skew){
      //DX20190613 [OBSOLETE] minimizeCartesianDistance(cpos,tmp,fdiff,c2f,f2c,_eps_);
      xvector<double> tmp_cpos = f2c*fpos;
      xmatrix<double> lattice = trasp(f2c); //DX20190613 
      xvector<double> min_cdiff = minimizeDistanceCartesianMethod(cpos,tmp_cpos,lattice); //DX20190613
      min_dist = aurostd::modulus(min_cdiff); //DX20190613
    }
    else {
      //DX20190613 [OBSOLETE] PBC(fdiff);
      xvector<double> tmp_fpos = c2f*cpos; //DX20190613
      xvector<double> min_fdiff = minimizeDistanceFractionalMethod(tmp_fpos,fpos); //DX20190613
      min_dist = aurostd::modulus(f2c*min_fdiff); //DX20190613
    }
    return min_dist<=_eps_; //DX20190613
  }
} // namespace SYM

// ******************************************************************************
// MapAtomWithBasis (vector)
// ******************************************************************************
// Map atoms and get basis map information
// fast==FALSE ensures we map one-to-one, otherwise you need to check yourself later
// MapAtomWithBasis (xvector<_atom>)
namespace SYM {
  //[CO20190515 - not needed and is ambiguous with overload]bool MapAtomWithBasis(vector<_atom>& vec, _atom& a, bool map_types, deque<uint>& index_to_check, xmatrix<double>& c2f,
  //[CO20190515 - not needed and is ambiguous with overload]                     xmatrix<double>& f2c, bool skew, double tol,bool fast){
  //[CO20190515 - not needed and is ambiguous with overload]  uint mapped_index=0;
  //[CO20190515 - not needed and is ambiguous with overload]  return MapAtomWithBasis(vec, a, map_types, index_to_check, c2f, f2c, skew, tol, mapped_index,fast);
  //[CO20190515 - not needed and is ambiguous with overload]}
} // namespace SYM

// MapAtomWithBasis (xvector<_atom>) with mapping index
namespace SYM {
  //overload //DX20190620
  bool MapAtomWithBasis(const vector<_atom>& vec, const _atom& a, bool map_types, deque<uint>& index_to_check, const xmatrix<double>& lattice, 
      bool skew, double tol, uint& mapped_index,bool fast){ 
    xmatrix<double> f2c=trasp(lattice);
    return MapAtomWithBasis(vec, a, map_types, index_to_check, lattice, f2c, skew, tol, mapped_index, fast);
  }

  bool MapAtomWithBasis(const vector<_atom>& vec, const _atom& a, bool map_types, deque<uint>& index_to_check, const xmatrix<double>& lattice, 
      const xmatrix<double>& f2c, bool skew, double tol, uint& mapped_index,bool fast){ //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input, remove "Atom" prefix from name
    int count=0;

    if(map_types){
      for(uint i=0;i<index_to_check.size();i++){
        //_atom b = vec[index_to_check[i]];
        if(a.type==vec[index_to_check[i]].type){ //DX20170731 - Speed increase
          if(AtomsMapped(a,vec[index_to_check[i]],lattice,f2c,skew,tol)){ //type specific //DX20190619 - lattice and f2c as input 
            count++;
            mapped_index=i;
            if(fast){return TRUE;} //DX20170731 - Speed increase, check one-to-one after
          }
        }
      }
    } else {
      for(uint i=0;i<index_to_check.size();i++){
        //_atom b = vec[index_to_check[i]];
        if(FPOSMatch(a,vec[index_to_check[i]],lattice,skew,tol)){ //type specific //DX20190619 - lattice and f2c as input, remove "Atom" prefix from name
          count++;
          mapped_index=i;
          if(fast){return TRUE;} //DX20170731 - Speed increase, check one-to-one after
        }
      }
    }
    if(count == 1){
      return TRUE;
    }
    return FALSE;
  }
} // namespace SYM

// MapAtomWithBasis (deque<_atom>)
namespace SYM {
  //[CO20190515 - not needed and is ambiguous with overload]bool MapAtomWithBasis(deque<_atom>& vec, _atom& a, bool map_types, deque<uint>& index_to_check, xmatrix<double>& c2f, xmatrix<double>& f2c, 
  //[CO20190515 - not needed and is ambiguous with overload]                     bool skew, double tol,bool fast){
  //[CO20190515 - not needed and is ambiguous with overload]  uint mapped_index=0;
  //[CO20190515 - not needed and is ambiguous with overload]  return MapAtomWithBasis(vec, a, map_types, index_to_check, c2f, f2c, skew, tol, mapped_index,fast);
  //[CO20190515 - not needed and is ambiguous with overload]}
} // namespace SYM

// MapAtomWithBasis (deque<_atom>) with mappping index
namespace SYM {
  //overload //DX20190620
  bool MapAtomWithBasis(const deque<_atom>& deq, const _atom& a, bool map_types, deque<uint>& index_to_check, const xmatrix<double>& lattice, 
      bool skew, double tol, uint& mapped_index,bool fast){ 
    xmatrix<double> f2c=trasp(lattice);
    return MapAtomWithBasis(deq, a, map_types, index_to_check, lattice, f2c, skew, tol, mapped_index, fast);
  }

  bool MapAtomWithBasis(const deque<_atom>& deq, const _atom& a, bool map_types, deque<uint>& index_to_check, const xmatrix<double>& lattice,
      const xmatrix<double>& f2c, bool skew, double tol, uint& mapped_index,bool fast){ //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
    int count=0;

    if(map_types){
      for(uint i=0;i<index_to_check.size();i++){
        //_atom b = deq[index_to_check[i]];
        if(a.type==deq[index_to_check[i]].type){ //DX20170731 - Speed increase
          if(AtomsMapped(a,deq[index_to_check[i]],lattice,f2c,skew,tol)){ //DX20190619 - lattice and f2c as input
            count++;
            mapped_index=i;
            if(fast){return TRUE;} //DX20170731 - Speed increase, check one-to-one after
          }
        }
      }
    } else {
      for(uint i=0;i<index_to_check.size();i++){
        //_atom b = deq[index_to_check[i]];
        if(FPOSMatch(a,deq[index_to_check[i]],lattice,f2c,skew,tol)){ //DX20190619 - lattice and f2c as input, remove "Atom" prefix from name
          count++;
          mapped_index=i;
          if(fast){return TRUE;} //DX20170731 - Speed increase, check one-to-one after
        }
      }
    }
    if(count == 1){
      return TRUE;
    }
    return FALSE;
  }
} // namespace SYM

// ******************************************************************************
// MapAtom (Set and single)
// ******************************************************************************
// Map an atom to a set of atoms
// MapAtom (deque<_atom> and _atom)
namespace SYM {
  //DX20190620
  bool MapAtom(const deque<_atom>& a_deq, const _atom& b, bool map_types, const xmatrix<double>& lattice, bool skew, double tol){ //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - replace f2c and c2f with lattice
    xmatrix<double> f2c=trasp(lattice);
    return MapAtom(a_deq, b, map_types, lattice, f2c, skew, tol);
  } 

  bool MapAtom(const deque<_atom>& a_deq, const _atom& b, bool map_types, const xmatrix<double>& lattice, const xmatrix<double>& f2c, bool skew, double tol){ //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
    for(uint i=0;i<a_deq.size();i++){
      //_atom a = a_deq[i];
      if(MapAtom(a_deq[i], b, map_types, lattice, f2c, skew, tol)){ //DX20190619 - replace f2c and c2f with lattice
        return TRUE;
      }
    }
    return FALSE;
  }
} // namespace SYM

// MapAtom (vector<_atom> and _atom)
namespace SYM {
  //overload //DX20190620
  bool MapAtom(const vector<_atom>& a_vec, const _atom& b, bool map_types, const xmatrix<double>& lattice, bool skew, double tol){
    xmatrix<double> f2c=trasp(lattice);
    return MapAtom(a_vec, b, map_types, lattice, f2c, skew, tol);
  }

  bool MapAtom(const vector<_atom> a_vec, const _atom b, bool map_types, const xmatrix<double>& lattice, const xmatrix<double>& f2c, bool skew, double tol){ //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
    for(uint i=0;i<a_vec.size();i++){
      //_atom a = a_vec[i];
      if(MapAtom(a_vec[i], b, map_types, lattice, f2c, skew, tol)){ //DX20190619 - replace f2c and c2f with lattice
        return TRUE;
      }
    }
    return FALSE;
  }
} // namespace SYM

// ******************************************************************************
// MapAtom (single and single)
// ******************************************************************************
// Map one atom to another. Check to see if there is a SINGLE instance where two 
// atoms are equivalent (both type and coord are the same)
// MapAtom (xvector and xvector)
// NOT TYPE SPECIFIC BY IT'S VERY NATURE
namespace SYM {
  //overload //DX20190620
  bool MapAtom(const xvector<double>& a, const xvector<double>& b, const xmatrix<double>& lattice, bool skew, double tol){
    xmatrix<double> f2c=trasp(lattice);
    return MapAtom(a, b, lattice, f2c, skew, tol);
  }

  bool MapAtom(const xvector<double>& a, const xvector<double>& b, const xmatrix<double>& lattice, const xmatrix<double>& f2c, bool skew, double tol){ //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
    //_atom a_atom; a_atom.fpos=a; a_atom.cpos=f2c*a;
    //_atom b_atom; b_atom.fpos=b; b_atom.cpos=f2c*b;
    if(FPOSMatch(a, b, lattice, f2c, skew, tol)){ //DX20190619 - lattice and f2c as input, remove "Atom" prefix from name
      return TRUE;
    }
    return FALSE;
  }
} // namespace SYM

// MapAtom (_atom and _atom)
namespace SYM { 
  bool MapAtom(const _atom& a, const _atom& b, bool map_types, const xmatrix<double>& lattice, bool skew, double tol){ //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - replace f2c and c2f with lattice
    xmatrix<double> f2c=trasp(lattice);
    return MapAtom(a, b, map_types, lattice, f2c, skew, tol);
  }

  bool MapAtom(const _atom& a, const _atom& b, bool map_types, const xmatrix<double>& lattice, const xmatrix<double>& f2c, bool skew, double tol){ //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
    if(map_types){
      return AtomsMapped(a,b,lattice,f2c,skew,tol); //type specific //DX20190619 - replace f2c and c2f with lattice
    } else {
      return FPOSMatch(a.fpos,b.fpos,lattice,f2c,skew,tol); //DX20190619 - lattice and f2c as input, remove "Atom" prefix from name
    }
  }
} // namespace SYM

//DX20190905 [OBSOLETE] // ******************************************************************************
//DX20190905 [OBSOLETE] // mod_one (Modify double by 1; to keep in unit cell)
//DX20190905 [OBSOLETE] // ******************************************************************************
//DX20190905 [OBSOLETE] // Bring a component of coordinate in the cell, based on a tolerance
//DX20190905 [OBSOLETE] namespace SYM {
//DX20190905 [OBSOLETE]   double mod_one(double d){
//DX20190905 [OBSOLETE]     if(d==INFINITY || d!=d || d==-INFINITY ){
//DX20190905 [OBSOLETE]       cerr << "SYM::mod_one: ERROR: (+-)INF or NAN value" << endl;
//DX20190905 [OBSOLETE]       return d;
//DX20190905 [OBSOLETE]     }
//DX20190905 [OBSOLETE]     while(d>=1.0-_ZERO_TOL_){ //1e-10
//DX20190905 [OBSOLETE]       d=d-1.0;
//DX20190905 [OBSOLETE]     }
//DX20190905 [OBSOLETE]     while(d<-_ZERO_TOL_){ //1e-10
//DX20190905 [OBSOLETE]       d=d+1.0;
//DX20190905 [OBSOLETE]     }
//DX20190905 [OBSOLETE]     return d;
//DX20190905 [OBSOLETE]   }
//DX20190905 [OBSOLETE] }

//DX20190905 [OBSOLETE] // ******************************************************************************
//DX20190905 [OBSOLETE] // mod_one_atom (Modify fpos and ijk by 1; to keep in unit cell)
//DX20190905 [OBSOLETE] // ******************************************************************************
//DX20190905 [OBSOLETE] // Bring atom.fpos and atom.ijk in the cell, based on a tolerance
//DX20190905 [OBSOLETE] namespace SYM {
//DX20190905 [OBSOLETE]   _atom mod_one_atom(const _atom& atom_in){
//DX20190905 [OBSOLETE]     _atom atom;
//DX20190905 [OBSOLETE]     atom=atom_in;
//DX20190905 [OBSOLETE]     atom.ijk=atom_in.ijk;
//DX20190905 [OBSOLETE]     for(uint i=1;i<(uint)atom.fpos.rows+1;i++){
//DX20190905 [OBSOLETE]       if(atom.fpos[i]==INFINITY || atom.fpos[i]!=atom.fpos[i] || atom.fpos[i]==-INFINITY ){
//DX20190905 [OBSOLETE]         cerr << "SYM::mod_one_atom:ERROR: (+-)INF or NAN value" << endl;
//DX20190905 [OBSOLETE]         return atom;
//DX20190905 [OBSOLETE]       }
//DX20190905 [OBSOLETE]       while(atom.fpos[i]>=1-_ZERO_TOL_){
//DX20190905 [OBSOLETE]         atom.fpos[i]-=1.0;
//DX20190905 [OBSOLETE]         atom.ijk[i]++;
//DX20190905 [OBSOLETE]       }
//DX20190905 [OBSOLETE]       while(atom.fpos[i]<-_ZERO_TOL_){
//DX20190905 [OBSOLETE]         atom.fpos[i]+=1.0;
//DX20190905 [OBSOLETE]         atom.ijk[i]--;
//DX20190905 [OBSOLETE]       }
//DX20190905 [OBSOLETE]     }
//DX20190905 [OBSOLETE]     return atom;
//DX20190905 [OBSOLETE]   }
//DX20190905 [OBSOLETE] }

//DX20190905 [OBSOLETE] // ******************************************************************************
//DX20190905 [OBSOLETE] // mod_one_xvec (Modify xvector compoenents by 1; to keep in unit cell)
//DX20190905 [OBSOLETE] // ******************************************************************************
//DX20190905 [OBSOLETE] // Bring a coordinate in the cell, based on vector tolerance
//DX20190905 [OBSOLETE] namespace SYM {
//DX20190905 [OBSOLETE]   xvector<double> mod_one_xvec(xvector<double> a){
//DX20190905 [OBSOLETE]     xvector<double> b;
//DX20190905 [OBSOLETE]     b(1) = mod_one(a(1));
//DX20190905 [OBSOLETE]     b(2) = mod_one(a(2));
//DX20190905 [OBSOLETE]     b(3) = mod_one(a(3));
//DX20190905 [OBSOLETE]     return b;
//DX20190905 [OBSOLETE]   }
//DX20190905 [OBSOLETE] }

// ******************************************************************************
// Check for identity
// ******************************************************************************
// A set of symmetry operations should at least contain the idenity operation
namespace SYM {
  bool CheckForIdentity(const xstructure& xstr){
    //we can check just for the string because if they are stored in xstructure then they were typed
    for(uint i=0; i<xstr.pgroup.size(); i++){
      if(xstr.pgroup[i].str_Hermann_Mauguin == "1"){
        return TRUE;
      }
    } 
    return FALSE;
  }
}

// ******************************************************************************
// checkSuperCellLatticePoints 
// ******************************************************************************
// Check if the supercell has lattice points in the correct positions given the symmetry 
// of the primitive cell
namespace SYM {
  bool checkSuperCellLatticePoints(xstructure& xstr, int& num_lattice_points, char& centering, uint& expand_size){
    bool skew = isLatticeSkewed(xstr.lattice,xstr.dist_nn_min,xstr.sym_eps);
    bool reverse_tested = false;

    vector<xvector<double> > conv_lattice_points;
    xvector<double> tmp;
    tmp(1)=0.0; tmp(2)=0.0; tmp(3)=0.0; // (0,0,0)
    conv_lattice_points.push_back(tmp);

    if(num_lattice_points == 4){ //Face centered
      tmp(1)=0.5; tmp(2)=0.0; tmp(3)=0.0; // (1/2,0,0)
      conv_lattice_points.push_back(tmp);
      tmp(1)=0.0; tmp(2)=0.5; tmp(3)=0.0; // (0,1/2,0)
      conv_lattice_points.push_back(tmp);
      tmp(1)=0.0; tmp(2)=0.0; tmp(3)=0.5; // (0,0,1/2)
      conv_lattice_points.push_back(tmp);
    }
    else if(num_lattice_points == 3){ //Rhombohedral
      // obverse setting
      tmp(1)=2.0/3.0; tmp(2)=1.0/3.0; tmp(3)=1.0/3.0; // (2/3,1/3,1/3)
      conv_lattice_points.push_back(tmp);
      tmp(1)=1.0/3.0; tmp(2)=2.0/3.0; tmp(3)=2.0/3.0; // (1/3,2/3,2/3)
      conv_lattice_points.push_back(tmp);
    }
    else if(num_lattice_points == 2 && centering == 'I'){ //Body-centered
      tmp(1)=0.5; tmp(2)=0.5; tmp(3)=0.5; // (1/2,1/2,1/2)
      conv_lattice_points.push_back(tmp);
    }
    else if(num_lattice_points == 2 && (centering == 'A' || centering == 'B' || centering == 'C')){ //Base-centered
      tmp(1)=0.5; tmp(2)=0.0; tmp(3)=0.0; // (1/2,0,0)
      conv_lattice_points.push_back(tmp);
    }  

    // Compress lattice fractional coordinates by the size of the expansion 
    vector<xvector<double> > compressed_lattice_points;
    for(uint i=0;i<conv_lattice_points.size();i++){
      compressed_lattice_points.push_back(conv_lattice_points[i]/(double)expand_size);
    }  

    // Check if the compressed lattice points are in the set of supercell translations
    for(uint j=0;j<compressed_lattice_points.size();j++){
      bool matched = false;
      for(uint ii=0;ii<xstr.fgroup.size();ii++){
        if(aurostd::identical(xstr.fgroup[ii].Uf,xstr.fgroup[0].Uf)){ //DX20171207 use default xmatrix identical tol
          //DX20190613 [OBSOLETE] xvector<double> fdiff = xstr.fgroup[ii].ftau - compressed_lattice_points[j];
          double min_dist = 1e9;
          if(skew){
            xvector<double> tmp1 = xstr.f2c*xstr.fgroup[ii].ftau;
            xvector<double> tmp2 = xstr.f2c*compressed_lattice_points[j];
            xvector<double> min_cdiff = minimizeDistanceCartesianMethod(tmp1,tmp2,xstr.lattice); //DX20190613
            min_dist = aurostd::modulus(min_cdiff); //DX20190613
            //DX20190613 [OBSOLETE] minimizeCartesianDistance(tmp1,tmp2,fdiff,xstr.c2f,xstr.f2c,xstr.sym_eps);
          }
          else {
            xvector<double> fdiff = xstr.fgroup[ii].ftau - compressed_lattice_points[j]; //DX20190613
            xvector<double> min_fdiff = minimizeDistanceFractionalMethod(fdiff); //DX20190613
            min_dist = aurostd::modulus(xstr.f2c*min_fdiff); //DX20190613
            //DX20190613 [OBSOLETE] PBC(fdiff);
          }
          //DX20190613 [OBSOLETE] if(aurostd::modulus(xstr.f2c*fdiff)<=xstr.sym_eps)
          if(min_dist<=xstr.sym_eps) //DX20190613
          { //CO20200106 - patching for auto-indenting
            matched = true;
            break;
          }
        }
      }
      if(matched == false){
        if(num_lattice_points == 3 && !reverse_tested){ // Need to test reverse setting before returning
          // reverse setting
          tmp(1)=1.0/3.0; tmp(2)=2.0/3.0; tmp(3)=1.0/3.0; // (1/3,2/3,1/3)
          conv_lattice_points.push_back(tmp);
          tmp(1)=2.0/3.0; tmp(2)=1.0/3.0; tmp(3)=2.0/3.0; // (2/3,1/3,2/3)
          conv_lattice_points.push_back(tmp);
          reverse_tested = true;
        }
        else {
          // Not a supercell of a conventional lattice
          return false;
        }
      }
    }
    return true;
  }
}

// ******************************************************************************
// ComparePointGroupAndSpaceGroupString
// ******************************************************************************
// Check that the point group and space group are commensurate. (Take away 
// translation component and centering information of space group --> point group)
namespace SYM {
  bool ComparePointGroupAndSpaceGroupString(xstructure& xstr, int& multiplicity_of_primitive, bool& derivative_structure){
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(multiplicity_of_primitive<0){
      cerr << "SYM::ComparePointGroupAndSpaceGroupString: ERROR: Multiplicity of original to primitive cell is negative [dir=" << xstr.directory << "]. " << endl;
    }
    string space_group_string = "";
    if(xstr.space_group_ITC != 0){ 
      space_group_string = GetSpaceGroupName(xstr.space_group_ITC);   
    }
    else {
      derivative_structure = false;
      if(LDEBUG) {
        cerr << "SYM::ComparePointGroupAndSpaceGroupString: WARNING: The space group was not found; cannot compare space group string and the point group of the crystal [dir=" << xstr.directory << "]." << endl;
      }
      return false;
    }

    // Determine lattice centering
    int num_lattice_pts = 0;
    char centering = space_group_string[0];
    if(centering == 'P'){ //simple
      num_lattice_pts = 1;
    }
    else if(centering == 'A' || centering == 'B' || centering == 'C' || centering == 'I'){ //base-/body-centered
      num_lattice_pts = 2;
    }
    else if(centering == 'R'){ //rhombohedral
      num_lattice_pts = 3;
    }
    else if(centering == 'F'){ //face-centered
      num_lattice_pts = 4;
    }

    // Determine if multiplicity of lattice points is consistent with a uniform expansion
    // Check for primitive
    uint n_prim=0;
    uint uniform_prim_lattice_points=0;
    while(uniform_prim_lattice_points<(uint)multiplicity_of_primitive){
      n_prim++;
      uniform_prim_lattice_points = std::pow(n_prim,3.0);
    }
    if(uniform_prim_lattice_points == (uint)multiplicity_of_primitive){
      derivative_structure = false;
    }    

    // Check for conventional
    if(derivative_structure && num_lattice_pts != 1){
      uint n_conv=0;
      uint uniform_conv_lattice_points=0;
      while(uniform_conv_lattice_points<(uint)multiplicity_of_primitive){
        n_conv++;
        uniform_conv_lattice_points = num_lattice_pts*std::pow(n_conv,3.0);
      }
      if(uniform_conv_lattice_points == (uint)multiplicity_of_primitive){
        if(checkSuperCellLatticePoints(xstr,num_lattice_pts, centering, n_conv)){
          derivative_structure = false;
        }
      }    
    }

    string tmp_point_group = "";
    bool screw_operator = false;
    for(uint i=0;i<space_group_string.size();i++){
      if(i > 0){ // Do not need lattice type
        // Reduce glide to mirror
        if(space_group_string[i] == 'a' || space_group_string[i] == 'b' || space_group_string[i] == 'c' ||
            space_group_string[i] == 'd' || space_group_string[i] == 'e' || space_group_string[i] == 'n' ){
          tmp_point_group += "m";
        }
        // Reduce screw to rotation
        else if(space_group_string[i] == '_'){
          screw_operator = true;
        }
        else if(space_group_string[i] == '}'){
          screw_operator = false;
        }
        else if(screw_operator == false){
          if(space_group_string[i] == '1' && space_group_string.size()>3){
            continue;
          } 
          else {
            tmp_point_group += space_group_string[i]; 
          } 
        }
      }
    }
    if(tmp_point_group == xstr.point_group_Hermann_Mauguin){
      return true;
    }
    else { // 2 and m are sometimes interchanged (similar "ranking" in Hermann Mauguin convention)
      string tmp_point_group_switch = "";
      for(uint i=0;i<tmp_point_group.size();i++){
        if(tmp_point_group[i] == '2'){
          tmp_point_group_switch += "m";
        }
        else if(tmp_point_group[i] == 'm'){
          tmp_point_group_switch += "2";
        }
        else {
          tmp_point_group_switch += tmp_point_group[i];
        }
      }
      if(tmp_point_group_switch == xstr.point_group_Hermann_Mauguin){
        return true;
      }
    }
    return false;
  }
} // namespace SYM
//DX+CO NEW - END
// ============================================================================

// ----------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------- SYMMETRY APPLICATION OPERATIONS
// Function SYM::ApplyAtom
//
// This routine takes one atom (type _atom) and applies the
// point group symmetry.. Depending on the flag, it brings back
// the atom to the unit cell (hence the calculation is done in
// fractional coordinates). SC Aug2007
//CO Apr2017 (wow 10 years later...)
// added additional flags for roff and ignoreFractionalOperation
// the roff flag toggles whether or not to roundoff the positions of the atom (if near 0, make 0)
// usually, we do NOT want to round off the actual number, just how it prints (0 vs. 1e-15)
// ignoreFractionalOperation toggles whether to transform by Uf of the _sym_op
// we want this ON if we're propagating the symmetry of a primitive cell to a (non-uniform) supercell (i.e., 4x4x2)
// with such expansions (derivative structure), the lattice loses symmetry, but the (unit cell) cystal does not
// we cannot simply apply Uf to the positions of the atom
// instead, we apply Uc (which always works), and do c2f
namespace SYM {
  //DX+CO START
  bool ApplyAtomValidate(const _atom &atom_in,_atom& atom_out,
      const _sym_op &symop,const xstructure& a){
    bool _incell_=false;
    bool roff=true; //CO, hate this, but it was the default, NEVER fudge positions
    return ApplyAtomValidate(atom_in,atom_out,symop,a,_incell_,roff);
  }
  bool ApplyAtomValidate(const _atom &atom_in,_atom& atom_out,
      const _sym_op &symop,const xstructure& a,bool _incell_,bool roff){
    //xmatrix<double> lattice = a.lattice;
    //xmatrix<double> f2c(3,3);f2c=a.f2c;//trasp(a.lattice); //  cpos=f2c*fpos;
    //xmatrix<double> c2f(3,3);c2f=a.c2f;//inverse(trasp(a.lattice));

    double _eps_ = _EPS_;
    if(a.sym_eps!=AUROSTD_NAN){ //Tolerance came from user or was calculated
      _eps_=a.sym_eps;
    }
    else {
      _eps_=_EPS_; 
    }
    double min_dist = a.dist_nn_min;
    if(min_dist==AUROSTD_NAN){min_dist=SYM::minimumDistance(a);}
    //if(a.dist_nn_min == AUROSTD_NAN){a.MinDist();}  //CO20180409
    //  min_dist = a.dist_nn_min; //CO20180409
    bool skew = isLatticeSkewed(a.lattice,min_dist,_eps_);
    return ApplyAtomValidate(atom_in,atom_out,symop,a.lattice,a.c2f,a.f2c,skew,_incell_,roff,_eps_);
  }
  bool ApplyAtomValidate(const _atom &atom_in,_atom& atom_out,
      const _sym_op &symop,const xstructure& a,bool skew, bool _incell_,bool roff, double _eps_){
    return ApplyAtomValidate(atom_in,atom_out,symop,a.lattice,a.c2f,a.f2c,skew,_incell_,roff,_eps_);
  }
  bool ApplyAtomValidate(const _atom &atom_in,_atom& atom_out,
      const _sym_op &symop,const xmatrix<double>& lattice,const xmatrix<double>& c2f, const xmatrix<double>& f2c,
      bool skew, bool _incell_,bool roff, double _eps_){
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    bool in_valid = validateAtomPosition(atom_in,c2f,f2c,skew,_eps_);
    if(!in_valid){
      if(LDEBUG) {
        cerr << "a.fpos " << atom_in.fpos << endl;
        cerr << "a.cpos " << atom_in.cpos << endl;
        cerr << "c2f " << c2f << endl;
        cerr << "f2c " << f2c << endl;
        xvector<double>fdiff=(c2f*atom_in.cpos)-atom_in.fpos;
        cerr << "(c2f*cpos)-fpos = " << (c2f*atom_in.cpos)-atom_in.fpos << endl;
        xvector<double> min_fdiff = minimizeDistanceFractionalMethod(fdiff); //DX20190613
        cerr << "minimizeDistanceFractionalMethod(fdiff) = " << min_fdiff << endl; //DX20190613
        cerr << "modulus(f2c*min_fdiff) = " << modulus(f2c*min_fdiff) << " <? " << _eps_ <<endl; //DX20190613
        //DX20190613 [OBSOLETE] PBC(fdiff);
        //DX20190613 [OBSOLETE] cerr << "PBC(fdiff) = " << fdiff << endl;
        //DX20160613 [OBSOLETE] cerr << "modulus(f2c*fdiff) = " << modulus(f2c*fdiff) << " <? " << _eps_ <<endl;
      }
      return false;
    }
    atom_out = ApplyAtom(atom_in,symop,lattice,c2f,f2c,skew,_incell_,roff,false,_eps_);
    bool out_valid = validateAtomPosition(atom_out,c2f,f2c,skew,_eps_);
    if(!out_valid){
      return false;
    }
    return true; 
  }
}
//DX+CO END
namespace SYM {
  _atom ApplyAtom(const _atom &atom_in,const _sym_op &symop,const xstructure& str) {
    bool _incell_=false;
    return SYM::ApplyAtom(atom_in,symop,str,_incell_);
  }

  _atom ApplyAtom(const _atom &atom_in,const _sym_op &symop,const xstructure& str,bool _incell_) {
    bool roff=true;
    return ApplyAtom(atom_in,symop,str,_incell_,roff);
  }

  //DX+CO START
  _atom ApplyAtom(const _atom &atom_in,const _sym_op &symop,const xstructure& str,bool _incell_,bool roff) {
    bool validatePosition=true;
    return ApplyAtom(atom_in,symop,str,_incell_,roff,validatePosition);
  }

  _atom ApplyAtom(const _atom &atom_in,const _sym_op &symop,const xstructure& str,bool _incell_,bool roff,bool validatePosition) {
    //CO - WE SHOULD PULL FROM XSTR INSTEAD OF RECALCUALTING EVERYTIME ?
    xmatrix<double> lattice = str.lattice;
    xmatrix<double> f2c(3,3);f2c=str.f2c;//trasp(str.lattice); //  cpos=f2c*fpos;
    xmatrix<double> c2f(3,3);c2f=str.c2f;//inverse(trasp(str.lattice));

    double _eps_ = _EPS_;
    if(str.sym_eps!=AUROSTD_NAN){ //Tolerance came from user or was calculated
      _eps_=str.sym_eps;
    }
    else {
      _eps_=_EPS_; 
    }
    double min_dist = str.dist_nn_min;
    if(min_dist == AUROSTD_NAN){min_dist=SYM::minimumDistance(str);}
    bool skew = isLatticeSkewed(lattice,min_dist,_eps_);

    return ApplyAtom(atom_in,symop,lattice,c2f,f2c,skew,_incell_,roff,validatePosition,_eps_);
  }

  _atom ApplyAtom(const _atom &atom_in,const _sym_op &symop,const xmatrix<double>& lattice,const xmatrix<double>& c2f, const xmatrix<double>& f2c,bool skew, bool _incell_,bool roff, bool validatePosition, double tolerance) {
    return ApplyAtom_20161115(atom_in,symop,lattice,c2f,f2c,skew,_incell_,roff,validatePosition,tolerance);
  }

  _atom ApplyAtom_20161115(const _atom &atom_in,const _sym_op &symop,const xmatrix<double>& lattice,const xmatrix<double>& c2f, const xmatrix<double>& f2c,bool skew, bool _incell_,bool roff,bool validatePosition, double _eps_) {

    stringstream message;
    _atom atom;
    atom=atom_in;  // copies all the info !
    //[OBSOLETE] check check=symop.is_pgroup+symop.is_fgroup+symop.is_sgroup+symop.is_agroup;
    int check = int(symop.is_pgroup)+int(symop.is_pgroupk)+int(symop.is_fgroup)+int(symop.is_pgroup_xtal)+int(symop.is_pgroupk_xtal)+int(symop.is_sgroup)+int(symop.is_agroup); //DX+CO //DX20191207 - added pgroupk_xtal

    bool identity=false;
    //bool roff=TRUE;

    if(check!=1) {
      message << "Error [1]" << endl;
      message << "  symop.is_pgroup=" << symop.is_pgroup << endl;
      message << "  symop.is_pgroupk=" << symop.is_pgroupk << endl;
      message << "  symop.is_pgroup_xtal=" << symop.is_pgroup_xtal << endl;
      message << "  symop.is_pgroupk_xtal=" << symop.is_pgroupk_xtal << endl; //DX20191207 - added pgroupk_xtal
      message << "  symop.is_fgroup=" << symop.is_fgroup << endl;
      message << "  symop.is_sgroup=" << symop.is_sgroup << endl;
      message << "  symop.is_agroup=" << symop.is_agroup << endl;
      message << "Error in symop.is_pgroup symop.is_fgroup symop.is_sgroup symop.is_agroup " << endl;
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
    }
    //check original atom first
    /*xvector<double> fdiff = (c2f*atom.cpos)-atom.fpos;
      xvector<double> tmp;
      if(skew){
      tmp = f2c*atom.fpos;
      minimizeCartesianDistance(atom.cpos,tmp,fdiff,c2f,f2c,_eps_);
      }
      else {
      PBC(fdiff);
      }
      if(modulus(f2c*fdiff)>_eps_) {*/
    if(validatePosition){
      if(!validateAtomPosition(atom,c2f,f2c,skew,_eps_)){
        message << "Error [2]" << endl;
        message << "  atom_in.cpos=" << atom_in.cpos << endl;
        message << "  atom_in.fpos=" << atom_in.fpos << endl;
        message << "  f2c=" << endl << f2c << endl;
        message << "  f2c*atom_in.fpos=" << f2c*atom.fpos << endl;
        message << "  c2f=" << endl << c2f << endl;
        message << "  c2f*atom_in.cpos=" << c2f*atom.cpos << endl;
        message << "  f2c*c2f=" << endl << f2c*c2f << endl;
        message << "  modulus = " << modulus(f2c*((c2f*atom.cpos)-atom.fpos)) << endl;
        message << "  EPS=" << _eps_ << endl;
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
      }
    }

    atom.ijk=atom_in.ijk;
    if(symop.is_pgroup==TRUE || symop.is_pgroup_xtal==TRUE || symop.is_pgroupk==TRUE || symop.is_pgroupk_xtal==TRUE || symop.is_agroup==TRUE) { //DX20191207 - added pgroupk_xtal
      // if(symop.is_pgroup==TRUE) cerr << "point group symmetry" << endl;
      // if(symop.is_agroup==TRUE) cerr << "site point group symmetry" << endl
      //if(symop.str_Hermann_Mauguin != "1"){
      if(!isidentity(symop.Uf)){
        atom.cpos=symop.Uc*atom_in.cpos;  //if(roff) roundoff(atom.cpos,_EPS_roundoff_);
        atom.fpos=symop.Uf*atom_in.fpos;  //if(roff) roundoff(atom.fpos,_EPS_roundoff_);
      }
      //  cerr << "ratom: " << atom.fpos << endl;

      //if(ignoreFractionalOperation){
      //  atom.fpos=c2f*atom.cpos;
      //} else {  //[CO20200106 - close bracket for indenting]}
      // atom.cpos=f2c*atom.fpos;
      //fdiff = (c2f*atom.cpos)-atom.fpos;
      //if(skew){
      //tmp = f2c*atom.fpos;
      //minimizeCartesianDistance(atom.cpos,tmp,fdiff,c2f,f2c,_eps_);
      //}
      //else {
      //PBC(fdiff);
      //}
      //if(modulus(f2c*fdiff)>_eps_) {

      if(validatePosition){
        if(!validateAtomPosition(atom,c2f,f2c,skew,_eps_)){
          message << "Error [3]" << endl;
          message << "  atom_in.cpos=" << atom_in.cpos << endl;
          message << "  atom_in.fpos=" << atom_in.fpos << endl;
          message << "  _sym_op=" << endl;
          message << symop << endl;
          message << "  atom.cpos=" << atom.cpos << endl;
          message << "  atom.fpos=" << atom.fpos << endl;
          message << "  f2c=" << endl << f2c << endl;
          message << "  f2c*atom.fpos=" << f2c*atom.fpos << endl;
          message << "  c2f=" << endl << c2f << endl;
          message << "  c2f*atom.cpos=" << c2f*atom.cpos << endl;
          message << "  f2c*c2f=" << endl << f2c*c2f << endl;
          message << "  modulus = " << modulus(f2c*((c2f*atom.cpos)-atom.fpos)) << endl;
          message << "  EPS=" << _eps_ << endl;
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
        }
      }
      //}
      if(_incell_) atom.ijk=xint(floor(atom.fpos));
      atom.basis=atom_in.basis;
    }
      else if(symop.is_fgroup==TRUE) {
        // cerr << "factor group symmetry" << endl;
        //if(symop.str_Hermann_Mauguin != "1")
        if(!isidentity(symop.Uf))
        { //[CO20200106 - close bracket for indenting]}
        atom.cpos=symop.Uc*atom_in.cpos; //symop.ctau+symop.Uc*atom_in.cpos;if(roff) roundoff(atom.cpos,_EPS_roundoff_);
        atom.fpos=symop.Uf*atom_in.fpos; //symop.ftau+symop.Uf*atom_in.fpos;if(roff) roundoff(atom.fpos,_EPS_roundoff_);
      }
      if(!(symop.ftau[1]<=_ZERO_TOL_ && symop.ftau[2]<=_ZERO_TOL_ && symop.ftau[3]<=_ZERO_TOL_)) {
        atom.cpos=atom.cpos+symop.ctau;
        atom.fpos=atom.fpos+symop.ftau;
      }
      //if(symop.str_Hermann_Mauguin == "1" && symop.ftau[1]<=_ZERO_TOL_ && symop.ftau[2]<=_ZERO_TOL_ && symop.ftau[3]<=_ZERO_TOL_)
      if(isidentity(symop.Uf) && symop.ftau[1]<=_ZERO_TOL_ && symop.ftau[2]<=_ZERO_TOL_ && symop.ftau[3]<=_ZERO_TOL_)
      { //[CO20200106 - close bracket for indenting]}
      identity=true;
    }
    //if(roff) roundoff(atom.fpos,_EPS_roundoff_);
    //if(roff) roundoff(atom.cpos,_EPS_roundoff_);
    //atom.fpos=symop.ftau+symop.Uf*atom_in.fpos;if(roff) roundoff(atom.fpos,_EPS_roundoff_);
    //atom.cpos=symop.ctau+symop.Uc*atom_in.cpos;if(roff) roundoff(atom.cpos,_EPS_roundoff_);
    //  atom.cpos=f2c*atom.fpos;

    //if(ignoreFractionalOperation){
    //  atom.fpos=c2f*atom.cpos;
    //} else {  //[CO20200106 - close bracket for indenting]}
    //fdiff = (c2f*atom.cpos)-atom.fpos;
    //if(skew){
    //tmp = f2c*atom.fpos;
    //minimizeCartesianDistance(atom.cpos,tmp,fdiff,c2f,f2c,_eps_);
    //}
    //else {
    //PBC(fdiff);
    //}
    //if(modulus(atom.cpos-f2c*atom.fpos)>_eps_)
    //if(modulus(f2c*fdiff)>_eps_) 
    //{

    if(validatePosition){
      if(!validateAtomPosition(atom,c2f,f2c,skew,_eps_)){
        message << "Error [4]" << endl;
        message << "  atom_in.cpos=" << atom_in.cpos << endl;
        message << "  atom_in.fpos=" << atom_in.fpos << endl;
        message << "  _sym_op=" << endl;
        message << symop << endl;
        message << "  atom.cpos=" << atom.cpos << endl;
        message << "  atom.fpos=" << atom.fpos << endl;
        message << "  f2c=" << endl << f2c << endl;
        message << "  f2c*atom.fpos=" << f2c*atom.fpos << endl;
        message << "  c2f=" << endl << c2f << endl;
        message << "  c2f*atom.cpos=" << c2f*atom.cpos << endl;
        message << "  f2c*c2f=" << endl << f2c*c2f << endl;
        message << "  modulus = " << modulus(f2c*((c2f*atom.cpos)-atom.fpos)) << endl;
        message << "  EPS=" << _eps_ << endl;
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
      }
    }
    //}
    if(_incell_) atom.ijk=xint(floor(atom.fpos));
    if(symop.basis_map_calculated){ //FIXING FOR ME CO20180420
      if((uint) atom_in.basis>=symop.basis_atoms_map.size()) {
        message << "Error [5]" << endl;
        message << "  atom_in.basis=" << atom_in.basis << endl;
        message << "  symop.basis_atoms_map.size()=" << symop.basis_atoms_map.size() << endl;
        message << "  symop.basis_types_map.size()=" << symop.basis_types_map.size() << endl;
        message << "  fgroup is mapping out of range" << endl;
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
      }
      atom.basis=symop.basis_atoms_map.at(atom_in.basis);
    }
  }
    else if(symop.is_sgroup==TRUE) {
      // cerr << "space group symmetry" << endl;
      //if(symop.str_Hermann_Mauguin != "1")
      if(!isidentity(symop.Uf))
      { //CO20200106 - patching for auto-indenting
        atom.cpos=symop.Uc*atom_in.cpos; //symop.ctau+symop.Uc*atom_in.cpos;if(roff) roundoff(atom.cpos,_EPS_roundoff_);
        atom.fpos=symop.Uf*atom_in.fpos; //symop.ftau+symop.Uf*atom_in.fpos;if(roff) roundoff(atom.fpos,_EPS_roundoff_);
      }
      if(!(symop.ftau[1]<=_ZERO_TOL_ && symop.ftau[2]<=_ZERO_TOL_ && symop.ftau[3]<=_ZERO_TOL_)) {
        atom.cpos=atom.cpos+symop.ctau;
        atom.fpos=atom.fpos+symop.ftau;
      }
      if(!(symop.ftrasl[1]<=_ZERO_TOL_ && symop.ftrasl[2]<=_ZERO_TOL_ && symop.ftrasl[3]<=_ZERO_TOL_)) {
        atom.cpos=atom.cpos+symop.ctrasl;
        atom.fpos=atom.fpos+symop.ftrasl;
      }
      //if(symop.str_Hermann_Mauguin == "1" && symop.ftau[1]<=_ZERO_TOL_ && symop.ftau[2]<=_ZERO_TOL_ && symop.ftau[3]<=_ZERO_TOL_ &&
      if(isidentity(symop.Uf) && symop.ftau[1]<=_ZERO_TOL_ && symop.ftau[2]<=_ZERO_TOL_ && symop.ftau[3]<=_ZERO_TOL_ &&
          symop.ftrasl[1]<=_ZERO_TOL_ && symop.ftrasl[2]<=_ZERO_TOL_ && symop.ftrasl[3]<=_ZERO_TOL_){
        identity=true;
      }
      //if(roff) roundoff(atom.fpos,_EPS_roundoff_);
      //if(roff) roundoff(atom.cpos,_EPS_roundoff_);

      //atom.fpos=symop.ftau+symop.ftrasl+symop.Uf*atom_in.fpos;if(roff) roundoff(atom.fpos,_EPS_roundoff_);
      //atom.cpos=symop.ctau+symop.ctrasl+symop.Uc*atom_in.cpos;if(roff) roundoff(atom.cpos,_EPS_roundoff_);
      // atom.cpos=f2c*atom.fpos;

      //if(ignoreFractionalOperation){
      //  atom.fpos=c2f*atom.cpos;
      //} else {  //[CO20200106 - close bracket for indenting]}
      //fdiff = (c2f*atom.cpos)-atom.fpos;
      //if(skew){
      //tmp = f2c*atom.fpos;
      //minimizeCartesianDistance(atom.cpos,tmp,fdiff,c2f,f2c,_eps_);
      //}
      //else {
      //PBC(fdiff);
      //}

      //if(modulus(atom.cpos-f2c*atom.fpos)>_eps_)
      //if(modulus(f2c*fdiff)>_eps_) 
      //{

      //CO+DX20170810
      //do not validate for sgroups, requires large search radius in the worse case (minimizeCartesianDistance)
      //if you did want to validate, need to change validateAtomPosition and subtract ctrasl/ftrasl first
      //Uc*(cpos-ctrasl) && Uf*(fpos-ftrasl)
      //this is SAFE because fgroups already validated, and sgroup is simply a shift to another cell
      if(0&&validatePosition){  
        if(!validateAtomPosition(atom,c2f,f2c,skew,_eps_)){
          message << "Error [6]" << endl;
          message << "  atom_in.cpos=" << atom_in.cpos << endl;
          message << "  atom_in.fpos=" << atom_in.fpos << endl;
          message << "  _sym_op=" << endl;
          message << symop << endl;
          message << "  atom.cpos=" << atom.cpos << endl;
          message << "  atom.fpos=" << atom.fpos << endl;
          message << "  f2c=" << endl << f2c << endl;
          message << "  f2c*atom.fpos=" << f2c*atom.fpos << endl;
          message << "  c2f=" << endl << c2f << endl;
          message << "  c2f*atom.cpos=" << c2f*atom.cpos << endl;
          message << "  f2c*c2f=" << endl << f2c*c2f << endl;
          message << "  modulus = " << modulus(f2c*((c2f*atom.cpos)-atom.fpos)) << endl;
          message << "  EPS=" << _eps_ << endl;
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
        }
      }
      //}
      if(_incell_) atom.ijk=xint(floor(atom.fpos));
      if((uint) atom_in.basis>=symop.basis_atoms_map.size()) {
        message << "Error [7]" << endl;
        message << "  atom_in.basis=" << atom_in.basis << endl;
        message << "  symop.basis_atoms_map.size()=" << symop.basis_atoms_map.size() << endl;
        message << "  symop.basis_types_map.size()=" << symop.basis_types_map.size() << endl;
        message << "sgroup is mapping out of range" << endl;
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
      }
      atom.basis=symop.basis_atoms_map.at(atom_in.basis);
    }
    if(roff) roundoff(atom.cpos,_EPS_roundoff_);
    if(roff) roundoff(atom.fpos,_EPS_roundoff_);
    //CO, fix this for sym_eps
    if(_incell_){
      if(!identity){
        atom=BringInCell(atom,lattice);
      }
    }
    //[OBSOLETE] if(_incell_) atom=BringInCell(atom,str.lattice);
    return atom;
  }
  //DX+CO END

#ifndef COMPILE_SLIM
  _atom ApplyAtom_20160101(const _atom &atom_in,const _sym_op &symop,const xstructure& str,bool _incell_) {
    stringstream message;
    _atom atom;
    atom=atom_in;  // copies all the info !
    char check=symop.is_pgroup+symop.is_fgroup+symop.is_sgroup+symop.is_agroup;
    xmatrix<double> f2c(3,3);f2c=trasp(str.lattice); //  cpos=f2c*fpos;
    xmatrix<double> c2f(3,3);c2f=inverse(trasp(str.lattice));

    bool roff=TRUE;

    if(check!=1) {
      message << "symop.is_pgroup=" << symop.is_pgroup << endl;
      message << "symop.is_pgroup_xtal=" << symop.is_pgroup_xtal << endl;
      message << "symop.is_fgroup=" << symop.is_fgroup << endl;
      message << "symop.is_sgroup=" << symop.is_sgroup << endl;
      message << "symop.is_agroup=" << symop.is_agroup << endl;
      message << "Error in symop.is_pgroup symop.is_fgroup symop.is_sgroup symop.is_agroup " << endl;
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
    }
    atom.ijk=atom_in.ijk;
    if(symop.is_pgroup==TRUE || symop.is_agroup==TRUE) {
      // if(symop.is_pgroup==TRUE) cerr << "point group symmetry" << endl;
      // if(symop.is_agroup==TRUE) cerr << "site point group symmetry" << endl;
      atom.fpos=symop.Uf*atom_in.fpos;if(roff) roundoff(atom.fpos,_EPS_roundoff_);
      atom.cpos=symop.Uc*atom_in.cpos;if(roff) roundoff(atom.cpos,_EPS_roundoff_);
      // atom.cpos=f2c*atom.fpos;
      if(modulus(atom.cpos-f2c*atom.fpos)>_EPS_) {
        message << "Error [1]" << endl;
        message << "  atom.cpos=" << atom.cpos << endl;
        message << "  atom.fpos=" << atom.fpos << endl;
        message << "  f2c=" << endl << f2c << endl;
        message << "  f2c*atom.fpos=" << f2c*atom.fpos << endl;
        message << "  c2f=" << endl << c2f << endl;
        message << "  c2f*atom.cpos=" << c2f*atom.cpos << endl;
        message << "  f2c*c2f=" << endl << f2c*c2f << endl;
        message << "  modulus = " << modulus(atom.cpos-f2c*atom.fpos) << endl;
        message << "  EPS=" << _EPS_ << endl;
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
      }
      if(_incell_) atom.ijk=xint(floor(atom.fpos));
      atom.basis=atom_in.basis;
    }
    if(symop.is_fgroup==TRUE) {
      // cerr << "factor group symmetry" << endl;
      atom.fpos=symop.ftau+symop.Uf*atom_in.fpos;if(roff) roundoff(atom.fpos,_EPS_roundoff_);
      atom.cpos=symop.ctau+symop.Uc*atom_in.cpos;if(roff) roundoff(atom.cpos,_EPS_roundoff_);
      //  atom.cpos=f2c*atom.fpos;
      if(modulus(atom.cpos-f2c*atom.fpos)>_EPS_) {
        message << atom.fpos << " | " << symop.ftau << " | " << symop.ftrasl << endl;
        message << atom.cpos << " | " << symop.ctau << " | " << symop.ctrasl << endl;
        message << f2c*atom.fpos << endl;
        message << modulus(atom.cpos-f2c*atom.fpos) << " " << _EPS_ << endl;
        message << "Error [2]" << endl;
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
      }
      if(_incell_) atom.ijk=xint(floor(atom.fpos));
      if((uint) atom_in.basis>=symop.basis_atoms_map.size()) {
        message << "atom_in.basis=" << atom_in.basis << endl;
        message << "symop.basis_atoms_map.size()=" << symop.basis_atoms_map.size() << endl;
        message << "symop.basis_types_map.size()=" << symop.basis_types_map.size() << endl;
        message << "fgroup is mapping out of range" << endl;
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
      }
      atom.basis=symop.basis_atoms_map.at(atom_in.basis);
    }
    if(symop.is_sgroup==TRUE) {
      // cerr << "space group symmetry" << endl;
      atom.fpos=symop.ftau+symop.ftrasl+symop.Uf*atom_in.fpos;if(roff) roundoff(atom.fpos,_EPS_roundoff_);
      atom.cpos=symop.ctau+symop.ctrasl+symop.Uc*atom_in.cpos;if(roff) roundoff(atom.cpos,_EPS_roundoff_);
      // atom.cpos=f2c*atom.fpos;
      if(modulus(atom.cpos-f2c*atom.fpos)>_EPS_) {
        message << atom.fpos << " | " << symop.ftau << " | " << symop.ftrasl << endl;
        message << atom.cpos << " | " << symop.ctau << " | " << symop.ctrasl << endl;
        message << f2c*atom.fpos << endl;
        message << modulus(atom.cpos-f2c*atom.fpos) << " " << _EPS_ << endl;
        message << "Error [3]" << endl;
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
      }
      if(_incell_) atom.ijk=xint(floor(atom.fpos));
      if((uint) atom_in.basis>=symop.basis_atoms_map.size()) {
        message << "atom_in.basis=" << atom_in.basis << endl;
        message << "symop.basis_atoms_map.size()=" << symop.basis_atoms_map.size() << endl;
        message << "symop.basis_types_map.size()=" << symop.basis_types_map.size() << endl;
        message << "sgroup is mapping out of range" << endl;
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
      }
      atom.basis=symop.basis_atoms_map.at(atom_in.basis);
    }
    if(roff) roundoff(atom.cpos,_EPS_roundoff_);
    if(roff) roundoff(atom.fpos,_EPS_roundoff_);
    if(_incell_) atom=BringInCell(atom,str.lattice);
    return atom;
  }
#endif

} // namespace SYM

namespace SYM {
  xvector<double> ApplyCpos(const xvector<double> &cpos_in,const _sym_op &symop,const xstructure& str,bool _incell_) {
    _atom atom;
    atom.cpos=cpos_in; atom.fpos=C2F(str.lattice,atom.cpos);
    atom=SYM::ApplyAtom(atom,symop,str,_incell_);
    return atom.cpos;
  }
} // namespace SYM

namespace SYM {
  xvector<double> ApplyCpos(const xvector<double> &cpos_in,const _sym_op &symop,const xstructure& str) {
    return SYM::ApplyCpos(cpos_in,symop,str,FALSE);
  }
} // namespace SYM

namespace SYM {
  xvector<double> ApplyFpos(const xvector<double> &fpos_in,const _sym_op &symop,const xstructure& str,bool _incell_) {
    _atom atom;
    atom.fpos=fpos_in; atom.cpos=F2C(str.lattice,atom.fpos);
    atom=SYM::ApplyAtom(atom,symop,str,_incell_);
    return atom.fpos;
  }
} // namespace SYM

namespace SYM {
  xvector<double> ApplyFpos(const xvector<double> &fpos_in,const _sym_op &symop,const xstructure& str) {
    return SYM::ApplyFpos(fpos_in,symop,str,FALSE);
  }
} // namespace SYM

namespace SYM {
  xvector<int> ApplyIJK(const xvector<int> &ijk_in,const _sym_op &symop,const xstructure& str) {

    stringstream message;
    xvector<double> fijk(3),cijk(3);
    xvector<int>    ijk(3);
    bool roff=TRUE;

    for(int i=1;i<=3;i++) {
      fijk[i]=ijk_in[i];
      cijk[i]=ijk_in[1]*str.lattice[1][i]+ijk_in[2]*str.lattice[2][i]+ijk_in[3]*str.lattice[3][i];
    }

    if(symop.is_pgroup==TRUE || symop.is_agroup==TRUE) {
      fijk=symop.Uf*fijk;if(roff) roundoff(fijk,_EPS_roundoff_);
      for(int i=1;i<=3;i++) ijk[i]=(int) fijk[i];

      cijk=symop.Uc*cijk;if(roff) roundoff(cijk,_EPS_roundoff_);
      if(modulus(fijk-cijk)>0.01) {
        message << "Mismatch error modulus(fijk-cijk)" << endl;
        message << fijk << endl;
        message << C2F(str.lattice,cijk) << endl;
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
      }
      for(int i=1;i<=3;i++) {
        if(abs(fijk[i]-ijk[i])>0.01) {
          message << "Mismatch error" << endl;
          message << ijk_in << endl;
          message << fijk << endl;
          message << ijk << endl;
          message << symop.Uf << endl;
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
        }
      }
    }
    if(symop.is_fgroup==TRUE) {message << "Not defined for fgroup"; throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);}
    if(symop.is_sgroup==TRUE) {message << "Not defined for sgroup"; throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);}
    return ijk;
  }
} // namespace SYM

namespace SYM {
  int ApplyL(const int &l,const _sym_op &symop,const xstructure& str) {
    if(str.lijk_calculated==FALSE) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"str.lijk_calculated must be calculated",_RUNTIME_ERROR_);
    }
    return ijk2l(str,SYM::ApplyIJK(l2ijk(str,l),symop,str));
  }
} // namespace SYM

namespace SYM {
  xstructure ApplyXstructure(const _sym_op &symop,const xstructure& str,bool _incell_) {
    xstructure str_out(str);
    for(uint iat=0;iat<str.atoms.size();iat++)
      str_out.atoms.at(iat)=SYM::ApplyAtom(str.atoms.at(iat),symop,str,_incell_);
    return str_out;
  }
} // namespace SYM

namespace SYM {
  xstructure ApplyXstructure(const _sym_op &symop,const xstructure& str) {
    return SYM::ApplyXstructure(symop,str,FALSE);
  }
} // namespace SYM

// ----------------------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------------------- EQUIVALENCE
// Function SYM::AtomsEquivalent
//
// This routine seels through all the space group and tells if two atoms
// (or positions) are equivalent
//  TYPE SPECIFIC
//DX+CO START
namespace SYM {
  bool AtomsEquivalent(xstructure& str,_atom& atom1,_atom& atom2) {
    return AtomsEquivalent_20161115(str,atom1,atom2,str.sym_eps);
  }

  bool AtomsEquivalent(xstructure& str, _atom& atom1, _atom& atom2, double& eps) {
    return AtomsEquivalent_20161115(str,atom1,atom2,eps);
  }

  bool AtomsEquivalent_20161115(xstructure& str,_atom& atom1,_atom& atom2) {
    double _eps_=AUROSTD_NAN; //DX =_EPS_;
    if(str.sym_eps!=AUROSTD_NAN){ //Tolerance came from user or was calculated
      _eps_=str.sym_eps;
    }
    else {
      _eps_=defaultTolerance(str);
      //_eps_=_EPS_;
    }
    return SYM::AtomsEquivalent(str,atom1,atom2,_eps_);
  }

#ifndef COMPILE_SLIM
  bool AtomsEquivalent_20160101(xstructure& str,_atom& atom1,_atom& atom2) {
    double eps = _EPS_;
    return SYM::AtomsEquivalent(str,atom1,atom2,eps);
    //return SYM::AtomsEquivalent(str,atom1,atom2,str.sym_eps);
  }
#endif

  bool AtomsEquivalent_20161115(xstructure& str, _atom& atom1, _atom& atom2, double& eps) {
    //we scale eps when possible (given lattice)
    //xmatrix<double> f2c=trasp(str.lattice); //co necessary?
    //xmatrix<double> c2f=inverse(trasp(str.lattice));  //co necessary?
    double min_dist = AUROSTD_NAN;
    if(str.dist_nn_min == AUROSTD_NAN){str.MinDist();}
    min_dist = str.dist_nn_min;
    bool skew = isLatticeSkewed(str.lattice,min_dist,eps);
    return AtomsEquivalent(str,atom1,atom2,skew,eps);
  }

#ifndef COMPILE_SLIM
  bool AtomsEquivalent_20160101(xstructure& str, _atom& atom1, _atom& atom2, double& eps) {
    // if(modulus(atom1.fpos-atom2.fpos)<ep s&& atom1.type==atom2.type) return TRUE;
    //DX20210111 [OBSOLETE - sgroup is not needed] if(!str.pgroup_calculated || !str.fgroup_calculated || !str.sgroup_calculated)
    if(!str.pgroup_calculated || !str.fgroup_calculated) { //DX20210111 - remove sgroup
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Symmetry must have been calculated",_RUNTIME_ERROR_);
    }
    _atom tatom;
    for(uint fg=0;fg<str.fgroup.size();fg++) {
      tatom=SYM::ApplyAtom(atom2,str.fgroup[fg],str);
      if(modulus(BringInCell(atom1.fpos-tatom.fpos))<eps) return TRUE;
    }
    //   for(uint sg=0;sg<str.sgroup.size();sg++) {
    //     tatom=SYM::ApplyAtom(atom2,str.sgroup[sg],str);
    //     if(modulus(atom1.cpos-tatom.cpos)<eps) return TRUE;
    //   }
    return FALSE;
  }
#endif

  bool AtomsEquivalent(xstructure& str, _atom& a, _atom& b, bool skew, double tol){ //CO20190520 - removed pointers for bools and doubles, added const where possible
    //DX20210111 [OBSOLETE - sgroup is not needed] if(!str.pgroup_calculated || !str.fgroup_calculated || !str.sgroup_calculated)
    if(!str.pgroup_calculated || !str.fgroup_calculated) { //DX20210111 - removed sgroup
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Symmetry must have been calculated",_RUNTIME_ERROR_);
    }
    _atom tatom;
    for(uint fg=0;fg<str.fgroup.size();fg++) {
      //tatom=SYM::ApplyAtom(b,str.fgroup[fg],str);
      if(!ApplyAtomValidate(b,tatom,str.fgroup[fg],str.lattice,str.c2f,str.f2c,skew,TRUE,false,tol)){
        return FALSE;
      }
      //tatom.fpos=BringInCell(tatom.fpos);
      if(AtomsMapped(a,tatom,str.lattice,str.f2c,skew,tol)){  //TYPE SPECIFIC BY DEFAULT //DX20190619 - lattice and f2c as input
        return TRUE;
      }
    }
    return FALSE;    
  }
}

namespace SYM {
  bool AtomsEquivalent_Basis(xstructure& str, int atom1_indx,int atom2_indx){
    //DX20210111 [OBSOLETE - sgroup is not needed] if(!str.pgroup_calculated || !str.fgroup_calculated || !str.sgroup_calculated)
    if(!str.pgroup_calculated || !str.fgroup_calculated) { //DX20210111 - remove sgroup
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Symmetry must have been calculated",_RUNTIME_ERROR_);
    }
    for(uint fg=0;fg<str.fgroup.size();fg++) {
      if(str.fgroup[fg].basis_atoms_map[atom1_indx]==atom2_indx || str.fgroup[fg].basis_atoms_map[atom2_indx]==atom1_indx){
        return TRUE;
      }
    }
    return FALSE;
  }
}
//DX+CO END

namespace SYM {
  bool CposEquivalent(const xstructure& str,const xvector<double> &cpos1,const xvector<double>& cpos2,const double& eps) {
    if(!str.pgroup_calculated || !str.fgroup_calculated || !str.sgroup_calculated) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Symmetry must have been calculated",_RUNTIME_ERROR_);
    }
    _atom iatom1,iatom2,tatom;
    iatom1.cpos=cpos1;iatom1.fpos=C2F(str.lattice,iatom1.cpos);
    iatom2.cpos=cpos2;iatom2.fpos=C2F(str.lattice,iatom2.cpos);
    for(uint sg=0;sg<str.sgroup.size();sg++) {
      tatom=SYM::ApplyAtom(iatom2,str.sgroup[sg],str);
      if(modulus(iatom1.cpos-tatom.cpos)<eps) return TRUE;
    }
    return FALSE;
  }
}

namespace SYM {
  bool FposEquivalent(const xstructure& str,const xvector<double>& fpos1,const xvector<double>& fpos2,const double& eps) {
    if(!str.pgroup_calculated || !str.fgroup_calculated || !str.sgroup_calculated) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Symmetry must have been calculated",_RUNTIME_ERROR_);
    }
    _atom iatom1,iatom2,tatom;
    iatom1.fpos=fpos1;iatom1.cpos=F2C(str.lattice,iatom1.fpos);
    iatom2.fpos=fpos2;iatom2.cpos=F2C(str.lattice,iatom2.fpos);
    for(uint sg=0;sg<str.sgroup.size();sg++) {
      tatom=SYM::ApplyAtom(iatom2,str.sgroup[sg],str);
      if(modulus(iatom1.fpos-tatom.fpos)<eps) return TRUE;
    }
    return FALSE;
  }
}

namespace SYM {
  bool CposEquivalent(const xstructure& str,const xvector<double>& cpos1,const xvector<double>& cpos2) {
    return SYM::CposEquivalent(str,cpos1,cpos2,_EPS_);
  }
}

namespace SYM {
  bool FposEquivalent(const xstructure& str,const xvector<double>& fpos1,const xvector<double>& fpos2) {
    return SYM::FposEquivalent(str,fpos1,fpos2,_EPS_);
  }
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------
// Function SYM::TypePointGroupOperation
// ----------------------------------------------------------------------------------------------------------------------------------------------------
//
// this function returns information about a point group
// symmetry operation - SC aug2007

namespace SYM {
  //bool TypePointGroupOperation(const xmatrix<double> &Uc,string &_string, bool &_inversion,double &_angle,xvector<double> &_axis,xmatrix<double> &_generator)
  bool TypePointGroupOperation(const xmatrix<double> &Uc, const xmatrix<double> &Uf, string &_string, bool &_inversion,double &_angle,xvector<double> &_axis,xmatrix<double> &_generator, xvector<double>& _generator_coefficients, xmatrix<xcomplex<double> >& _SU2_matrix, xvector<xcomplex<double> >& _su2_coefficients, double _eps_) //DX20180117 - added SU(2) and su(2) coefficients
  { //CO20200106 - patching for auto-indenting
    // calculate the type of symmetry.
    // you need to put your right hand on the origin,
    // and rotate as screwing a scredriver along direction r(1),r(2),r(3),
    // for theta angles.
    // theta along r is equivalent to 2pi-theta along -r... try to plot.
    ostringstream aus;
    //double _eps_=_EPS_,x,y,z; //DX20171024
    double x,y,z;
    // getting angle and vector of rotation
    bool inversion=FALSE;
    bool _roundoff_U=TRUE;
    bool _roundoff_rrr=TRUE;
    xmatrix<double>  A(3,3),U(3,3),I(3,3),II(3,3),Uexp(3,3);
    // U=trasp(Uc);  // U=inverse(U);
    U=Uc;
    // A rotation matrix has eigenvalues 1, cos_theta+i*sin_theta, cos_theta-i*sin_theta
    // and it is orthonormal U^-1=U^T
    // The trace of a matrix, which is the sum of the diagonal elements, is independent of the
    // coordinate system used, so long as that is orthonormal. This implies that the trace is
    // the sum of the eigenvalues. Here that is 1+cos_heta+i*sin_theta+cos_theta-i*sin_theta
    // trace=1 + 2*cos_theta .  Therefore, theta=acos((trace(U)-1)/2)
    // Note that this gives 2 possibilities, because acos gives angles between 0 to pi (and not pi to 2pi)
    xvector<double> r(3),ro(3);
    double theta;
    I=aurostd::identity(I);II=-I;
    if(_roundoff_U) roundoff(U,_EPS_roundoff_);
    //  cerr << " U= " << endl << U << endl;
    //DX20171207 Use Uf, elements are -1,0,1 : if(identical(U,I,_eps_))                                           // check unity  // fix the problem at theta=0.0
    if(identical(Uf,I))
    { //CO20200106 - patching for auto-indenting
      theta=0.0;r(1)=0.0;r(2)=0.0;r(3)=0.0;clear(A);  // ANGLE IS ZERO and axis is not defined
      inversion=FALSE;
      aus << "unity I       ";
      _inversion=inversion;_angle=theta*rad2deg;_axis(1)=r(1);_axis(2)=r(2);_axis(3)=r(3);_string=aus.str();_generator=A;
      _generator_coefficients(1)=r(1);_generator_coefficients(2)=r(2);_generator_coefficients(3)=r(3); //DX20171206 - added generator coefficents
      ComplexSU2Rotations(_SU2_matrix,_su2_coefficients,theta,_axis); //DX20180117 - added SU2 and su2 coefficients
      return TRUE;
    }
    //DX20171207 Use Uf, elements are -1,0,1 : if(identical(U,II,_eps_))                                           // check unity  // fix the problem at theta=0.0
    if(identical(Uf,II))
    { //CO20200106 - patching for auto-indenting
      theta=0.0;r(1)=0.0;r(2)=0.0;r(3)=0.0;clear(A);  // ANGLE IS ZERO and axis is not defined
      inversion=TRUE;
      aus << "inversion -I  ";
      _inversion=inversion;_angle=theta*rad2deg;_axis(1)=r(1);_axis(2)=r(2);_axis(3)=r(3);_string=aus.str();_generator=A;
      _generator_coefficients(1)=r(1);_generator_coefficients(2)=r(2);_generator_coefficients(3)=r(3); //DX20171206 - added generator coefficents
      ComplexSU2Rotations(_SU2_matrix,_su2_coefficients,theta,_axis); //DX20180117 - added SU2 and su2 coefficients
      return TRUE;
    }
    // not unity and not inversion, check for inversions
    if(det(U)<0) {                                                     // contains rotoinversion
      inversion=TRUE;
      U=-U; // flip the matrix
    }
    // check if theta = 0 or pi (0 already done)
    // if angle is pi, then U*U=I, hence inv(U)=U, but inv(U)=Ut, hence U=Ut.... symmetric
    //DX20171207 Use Uf, elements are -1,0,1 : if(identical(U,trasp(U),_eps_))    // if symmetric, theta = pi, because two operations U*U brings 2pi=I
    if(identical(Uf*Uf,I))
    { //CO20200106 - patching for auto-indenting
      double x,y,z;
      x=(U(1,1)+1.0)/2.0;             // definition of x
      if(x>0.0) x=sqrt(x);           // can only take square root of positive number, always TRUE for orthogonal matrix
      else x=0.0;                     // in case matrix has become de-orthogonalised
      y=(U(2,2)+1.0)/2.0;             // definition of y
      if(y>0.0) y=sqrt(y);           // can only take square root of positive number, always TRUE for orthogonal matrix
      else y=0.0;                     // in case matrix has become de-orthogonalised
      z=(U(3,3)+1.0)/2.0;             // definition of z
      if(z>0.0) z=sqrt(z);           // can only take square root of positive number, always TRUE for orthogonal matrix
      else z=0.0;                     // in case matrix has become de-orthogonalised
      if((abs(x)<_eps_) && !(abs(y)<_eps_) && !(abs(z)<_eps_)) { // implements  last 6 rows of above table
        if(!(U(2,3)>0.0)) y=-y;
      } else if((abs(y)<_eps_) && !(abs(z)<_eps_)) {
        if(!(U(1,3)>0.0)) z=-z;
      } else if((abs(z)<_eps_)) {
        if(!(U(1,2)>0.0)) x=-x;
      }
      theta=pi;
      r(1)=x;r(2)=y;r(3)=z;clear(A);
      if(r(3)<0.0) {r=-r;theta=2*pi-theta;}   //theta along r is equivalent to 2pi-theta along -r... try to plot.  CHOICE
      if(_roundoff_rrr) roundoff(r,_EPS_roundoff_);
      if(inversion==FALSE) aus << "rotation      ";//theta = " << (theta*rad2deg<100?" ":"") << theta*rad2deg << "  around r = [" << r << "]";
      if(inversion==TRUE ) aus << "rotoinversion ";//theta = " << (theta*rad2deg<100?" ":"") << theta*rad2deg << "  around r = [" << r << "]";
      A(1,1)=0;A(1,2)=-z;A(1,3)=y;A(2,1)=z;A(2,2)=0;A(2,3)=-x;A(3,1)=-y;A(3,2)=x;A(3,3)=0;A=A*theta;Uexp=aurostd::exp(A);
      //  cerr << Uexp << endl << max(aurostd::abs(Uexp-U)) << endl;
      _inversion=inversion;_angle=theta*rad2deg;_axis(1)=r(1);_axis(2)=r(2);_axis(3)=r(3);_string=aus.str();_generator=A;
      _generator_coefficients(1)=r(1);_generator_coefficients(2)=r(2);_generator_coefficients(3)=r(3); //DX20171206 - added generator coefficents
      ComplexSU2Rotations(_SU2_matrix,_su2_coefficients,theta,_axis); //DX20180117 - added SU2 and su2 coefficients
      return TRUE;
    } else {                          // not unity, not inversion, not theta=pi, generic case
      theta=((aurostd::trace(U)-1.0)/2.0);
      if(theta> 1.0) theta=1.0;
      if(theta<-1.0) theta=-1.0;
      theta=acos(theta);   // but could be somewhere else (mirror in the horrizontal plane).
      // theta=acos((aurostd::trace(U)-1.0)/2.0);  // but could be somewhere else (mirror in the horrizontal plane).
      // Here is another way to find the axis, if q is not 0deg  or 180deg. If M=mij is the matrix and a is the normalized axis, then
      // 2 a sin_theta= (m32-m23,m13-m31,m21-m12). Therefore, (m32-m23,m13-m31,m21-m12)
      // is an unnormalized axis, and its length is 2 sin_theta.
      // try to get the vector ... unnormalized
      // http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToAngle/index.htm
      // x=(m21-m12)/sqrt((m21-m12)^2+(m02-m20)^2+(m10-m01)^2)
      // y=(m02-m20)/sqrt((m21-m12)^2+(m02-m20)^2+(m10-m01)^2)
      // z=(m10-m01)/sqrt((m21-m12)^2+(m02-m20)^2+(m10-m01)^2)
      r(1)=U(3,2)-U(2,3);
      r(2)=U(1,3)-U(3,1);
      r(3)=U(2,1)-U(1,2);
      r=r/modulus(r);
      // make a choice. I choose to be atheist and to take only the rotation vectors with z positive ! (SC)
      if(r(3)<0.0) {r=-r;theta=2*pi-theta;}   //theta along r is equivalent to 2pi-theta along -r... try to plot.  CHOICE
      if(_roundoff_rrr) roundoff(r,_EPS_roundoff_);
      x=r(1);y=r(2);z=r(3);
      if(inversion==FALSE) aus << "rotation      ";//theta = " << (theta*rad2deg<100?" ":"") << theta*rad2deg << "  around r = [" << r << "]";
      if(inversion==TRUE ) aus << "rotoinversion ";//theta = " << (theta*rad2deg<100?" ":"") << theta*rad2deg << "  around r = [" << r << "]";
      // cerr << endl << U << endl;
      A(1,1)=0;A(1,2)=-z;A(1,3)=y;A(2,1)=z;A(2,2)=0;A(2,3)=-x;A(3,1)=-y;A(3,2)=x;A(3,3)=0;A=A*theta;Uexp=aurostd::exp(A);
      // cerr << Uexp << endl << max(aurostd::abs(Uexp-U)) << endl;
      _inversion=inversion;_angle=theta*rad2deg;_axis(1)=r(1);_axis(2)=r(2);_axis(3)=r(3);_string=aus.str();_generator=A;
      _generator_coefficients(1)=r(1);_generator_coefficients(2)=r(2);_generator_coefficients(3)=r(3); //DX20171206 - added generator coefficents
      ComplexSU2Rotations(_SU2_matrix,_su2_coefficients,theta,_axis); //DX20180117 - added SU2 and su2 coefficients
      // cerr << _generator << endl;
      return TRUE;
    }
    cerr << " U= " << endl << U << endl << " theta=" << theta << endl << " r=" << r << endl;
    cerr << endl << "ERROR SYM::TypePointGroupOperation, you should never be here... [0] !" << endl;
    return FALSE;
  }
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------
// Function SYM::TypePointGroupOperationInternational
// ----------------------------------------------------------------------------------------------------------------------------------------------------
// this function returns the International symbol = Hermann-Mauguin notation & Schönflies notation SC0210

namespace SYM {
  //bool TypePointGroupOperationInternational(const xmatrix<double>& Uc,string& _stringHM,string& _stringSC,const bool& _inversion,const double& _angle,
  //					    const xvector<double>& _axis,const xmatrix<double>& _generator)
  bool TypePointGroupOperationInternational(const xmatrix<double>& Uc,string& _stringHM,string& _stringSC,const bool& _inversion,const double& _angle,
      const xvector<double>& _axis,const xmatrix<double>& _generator, xvector<double>& _generator_coefficients, 
      xmatrix<xcomplex<double> >& _SU2_matrix, xvector<xcomplex<double> >& _su2_coefficients, double _eps_) //DX20171206 - added generator coefficents //DX20180117 - added SU(2) and su(2) coefficients
  { //CO20200106 - patching for auto-indenting
    //double _eps_=_EPS_; //DX20171024
    _stringHM="XXXXX";
    _stringSC="XXXXX";

    // get diagonal stuff    1 and -1
    if(isdiagonal(Uc,_eps_)) {
      if(isequal(Uc[1][1], 1.0,_eps_) && isequal(Uc[2][2], 1.0,_eps_) && isequal(Uc[3][3], 1.0,_eps_)) { _stringHM= "1"; _stringSC="1"; return TRUE;}
      if(isequal(Uc[1][1],-1.0,_eps_) && isequal(Uc[2][2], 1.0,_eps_) && isequal(Uc[3][3], 1.0,_eps_)) { _stringHM= "m"; _stringSC="s"; return TRUE;}
      if(isequal(Uc[1][1], 1.0,_eps_) && isequal(Uc[2][2],-1.0,_eps_) && isequal(Uc[3][3], 1.0,_eps_)) { _stringHM= "m"; _stringSC="s"; return TRUE;}
      if(isequal(Uc[1][1], 1.0,_eps_) && isequal(Uc[2][2], 1.0,_eps_) && isequal(Uc[3][3],-1.0,_eps_)) { _stringHM= "m"; _stringSC="s"; return TRUE;}  
      if(isequal(Uc[1][1],-1.0,_eps_) && isequal(Uc[2][2],-1.0,_eps_) && isequal(Uc[3][3],-1.0,_eps_)) { _stringHM="-1"; _stringSC="i"; return TRUE;}
    }
    if(det(Uc)>0.0) {
      if(isequal(_angle, 60.0,_eps_)) { _stringHM= "6"; _stringSC="C6"; return TRUE;}
      if(isequal(_angle, 90.0,_eps_)) { _stringHM= "4"; _stringSC="C4"; return TRUE;}
      if(isequal(_angle,120.0,_eps_)) { _stringHM= "3"; _stringSC="C3"; return TRUE;}
      if(isequal(_angle,180.0,_eps_)) { _stringHM= "2"; _stringSC="C2"; return TRUE;}
      if(isequal(_angle,240.0,_eps_)) { _stringHM= "3"; _stringSC="C3"; return TRUE;}
      if(isequal(_angle,270.0,_eps_)) { _stringHM= "4"; _stringSC="C4"; return TRUE;}
      if(isequal(_angle,300.0,_eps_)) { _stringHM= "6"; _stringSC="C6"; return TRUE;}
    }
    if(det(Uc)<0.0) {
      if(isequal(_angle, 60.0,_eps_)) { _stringHM= "-6"; _stringSC="S6"; return TRUE;}
      if(isequal(_angle, 90.0,_eps_)) { _stringHM= "-4"; _stringSC="S4"; return TRUE;}
      if(isequal(_angle,120.0,_eps_)) { _stringHM= "-3"; _stringSC="S3"; return TRUE;}
      if(isequal(_angle,180.0,_eps_)) { _stringHM= "-2"; _stringSC="S2"; return TRUE;}
      if(isequal(_angle,240.0,_eps_)) { _stringHM= "-3"; _stringSC="S3"; return TRUE;}
      if(isequal(_angle,270.0,_eps_)) { _stringHM= "-4"; _stringSC="S4"; return TRUE;}
      if(isequal(_angle,300.0,_eps_)) { _stringHM= "-6"; _stringSC="S6"; return TRUE;}
    }
    return FALSE;
    cerr << "_inversion=" << _inversion << endl;
    cerr << "_angle=" << _angle << endl;
    cerr << "_axis=" << _axis << endl;
    cerr << "_generator=" << _generator << endl;
    cerr << "_generator_coefficients=" << _generator_coefficients << endl;
    cerr << "_SU2_matrix=" << _SU2_matrix << endl; //DX20180117 - add SU(2)
    cerr << "_su2_coeffiients=" << _su2_coefficients << endl; //DX20180117 - add su(2) coefficients
    cerr << endl << "ERROR SYM::TypePointGroupOperationInternational, you should never be here... [0] !" << endl;
    return FALSE;
  }
}

uint HowMany(vector<string> vstr,string scheck) {
  uint out=0;
  for(uint i=0;i<vstr.size();i++)
    if(scheck==vstr.at(i)) out++;
  return out;
}

bool MapOperations(string s1,string s2) {
  vector<string> tokens1,tokens2;
  aurostd::string2tokens(s1,tokens1);
  aurostd::string2tokens(s2,tokens2);
  if(s1.size()!=s2.size()) return FALSE;
  uint maps=0;
  for(uint i1=0;i1<tokens1.size();i1++)
    for(uint i2=0;i2<tokens2.size();i2++)
      if(tokens1.at(i1)!="X" && tokens2.at(i2)!="X")
        if(tokens1.at(i1)==tokens2.at(i2)) { // nail
          tokens1.at(i1)="X";tokens2.at(i2)="X";
          maps++;
        }
  //  for(uint i1=0;i1<tokens1.size();i1++) cerr << tokens1.at(i1) << endl;
  //   for(uint i1=0;i1<tokens1.size();i1++) cerr << tokens2.at(i1) << endl;
  //   cerr << maps << endl;
  if(maps==tokens1.size()) return TRUE;
  return FALSE;
}

//DX+CO START
// ----------------------------------------------------------------------------------------------------------------------------------------------------
// Function SYM::PointGroupsIdentical
// ----------------------------------------------------------------------------------------------------------------------------------------------------
// this functions comapares whether the point groups are identical, ignores Uf
namespace SYM {
  bool PointGroupsIdentical(const vector<_sym_op>& vpg1,const vector<_sym_op>& vpg2, double eps, bool is_same_lattice){ //DX20171207 - added is_same_lattice
    if(vpg1.size()!=vpg2.size()){
      return FALSE;
    }
    bool found;
    //DX20171207 - compare Uf vs Uc : START
    //DX20171207 - If the pgroups come from the same lattice, then the Ufs would be exactly the same 
    // (Uf comprised of integers; no tolerance needed).  If lattices are not the same, use Uc and eps.
    if(is_same_lattice){
      for(uint i=0;i<vpg1.size();i++){
        found=false;
        for(uint j=0;j<vpg2.size() && !found;j++){
          if(aurostd::identical(vpg1[i].Uf,vpg2[j].Uf)){ //DX20171207 - Use Uf (only integers) not Uc and use xmatrix identical eps
            found=true;
          }
        }
        if(!found){
          cerr << "SYM::PointGroupsIdentical: DID NOT FIND ";
          cerr << vpg1[i] << endl;
          return FALSE;
        }
      }
    }
    else {
      for(uint i=0;i<vpg1.size();i++){
        found=false;
        for(uint j=0;j<vpg2.size() && !found;j++){
          if(aurostd::identical(vpg1[i].Uc,vpg2[j].Uc,eps)){
            found=true;
          }
        }
        if(!found){
          cerr << "SYM::PointGroupsIdentical: DID NOT FIND ";
          cerr << vpg1[i] << endl;
          return FALSE;
        }
      }
    }
    //DX20171207 - compare Uf vs Uc : START
    return TRUE;
  }
}
//DX+CO END

//GG START
// ----------------------------------------------------------------------------------------------------------------------------------------------------
// Function SYM::CalculateQuaternion
// ----------------------------------------------------------------------------------------------------------------------------------------------------
// this function fills in the quaternion attribute(s) to _sym_op
namespace SYM {
  bool CalculateQuaternion(_sym_op& a){
    //GG START WORKING HERE

    a.quaternion_vector.clear(); a.quaternion_matrix.clear();  //CO20170706 - clear first
    xvector<double> quaternion_vector(4); // creating the vector for 4x1
    xmatrix<double> quaternion_matrix(4, 4); // creating the matrix for 4x4
    a.quaternion_vector=quaternion_vector;
    a.quaternion_matrix=quaternion_matrix;

    double angleD, firstcomponent, secondcomponent, thirdcomponent, fourthcomponent; // declaring the variables of the components
    angleD = a.angle*((pi)/(180)); // converting angle from degrees to radians

    // Calculating components of quaternion 
    firstcomponent  = 1*cos(angleD/2);
    secondcomponent = a.axis[1]*sin(angleD/2);
    thirdcomponent  = a.axis[2]*sin(angleD/2);
    fourthcomponent = a.axis[3]*sin(angleD/2);

    // 4x1 Quaternion Vector
    a.quaternion_vector(1)=firstcomponent;
    a.quaternion_vector(2)=secondcomponent;
    a.quaternion_vector(3)=thirdcomponent;
    a.quaternion_vector(4)=fourthcomponent;

    // 4x4 Quaternion Matrix
    // Components of the first row
    a.quaternion_matrix(1, 1)= firstcomponent;
    a.quaternion_matrix(1, 2)= secondcomponent;
    a.quaternion_matrix(1, 3)= thirdcomponent;
    a.quaternion_matrix(1, 4)= fourthcomponent;

    // Components of the second row
    a.quaternion_matrix(2, 1)= -1*secondcomponent;
    a.quaternion_matrix(2, 2)= firstcomponent;
    a.quaternion_matrix(2, 3)= -1*fourthcomponent;
    a.quaternion_matrix(2, 4)= thirdcomponent;

    // Components of the third row
    a.quaternion_matrix(3, 1)= -1*thirdcomponent;
    a.quaternion_matrix(3, 2)= fourthcomponent;
    a.quaternion_matrix(3, 3)= firstcomponent;
    a.quaternion_matrix(3, 4)= -1*secondcomponent;

    // Components of the fourth row
    a.quaternion_matrix(4, 1)= -1*fourthcomponent;
    a.quaternion_matrix(4, 2)= -1*thirdcomponent;
    a.quaternion_matrix(4, 3)= secondcomponent;
    a.quaternion_matrix(4, 4)= firstcomponent;

    return TRUE;
  }
}
//GG STOP

//DX20180117 START
// ----------------------------------------------------------------------------------------------------------------------------------------------------
// Function SYM::ComplexSU2Rotations
// ----------------------------------------------------------------------------------------------------------------------------------------------------
// this function fills in the SU(2) and su(2) coefficient attribute to _sym_op
namespace SYM {
  bool ComplexSU2Rotations(xmatrix<xcomplex<double> > & _SU2_matrix, xvector<xcomplex<double> >& _su2_coefficients, double& theta, xvector<double>& _axis){

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    xcomplex<double> q11(1.0*cos(theta/2.0),_axis(3)*sin(theta/2.0));                    //DX20180117 - added SU(2) matrix
    xcomplex<double> q12(_axis(2)*sin(theta/2.0),_axis(1)*sin(theta/2.0));               //DX20180117 - added SU(2) matrix
    xcomplex<double> q21(-_axis(2)*sin(theta/2.0),_axis(1)*sin(theta/2.0));              //DX20180117 - added SU(2) matrix
    xcomplex<double> q22(1.0*cos(theta/2.0),-_axis(3)*sin(theta/2.0));                   //DX20180117 - added SU(2) matrix
    _SU2_matrix(1,1)=q11;_SU2_matrix(1,2)=q12;_SU2_matrix(2,1)=q21;_SU2_matrix(2,2)=q22; //DX20180117 - added SU(2) matrix
    xcomplex<double> coefficient1(0.0,_axis(1)/2.0);_su2_coefficients(1)=coefficient1;   //DX20180117 - added su(2) coefficients
    xcomplex<double> coefficient2(0.0,_axis(2)/2.0);_su2_coefficients(2)=coefficient2;   //DX20180117 - added su(2) coefficients
    xcomplex<double> coefficient3(0.0,_axis(3)/2.0);_su2_coefficients(3)=coefficient3;   //DX20180117 - added su(2) coefficients

    if(LDEBUG) {
      xmatrix<xcomplex<double> > complex_matrix(2,2);
      xcomplex<double> tmp(0.0,_axis(3)/2.0);
      complex_matrix(1,1) = theta*tmp;
      xcomplex<double> tmp1(_axis(2)/2.0,_axis(1)/2.0);
      complex_matrix(1,2) = theta*tmp1;
      xcomplex<double> tmp2(-_axis(2)/2.0,_axis(1)/2.0);
      complex_matrix(2,1) = theta*tmp2;
      xcomplex<double> tmp3(0.0,-_axis(3)/2.0);
      complex_matrix(2,2) = theta*tmp3;
      xmatrix<xcomplex<double> > exp_matrix = aurostd::exp(complex_matrix);
      //[CO20230319 - not sure what this is?]if(!(abs(_SU2_matrix(1,1).re-exp_matrix(1,1).re)<1e-3 && abs(_SU2_matrix(1,1).im-exp_matrix(1,1).im)<1e-3 &&
      //[CO20230319 - not sure what this is?]      abs(_SU2_matrix(1,2).re-exp_matrix(1,2).re)<1e-3 && abs(_SU2_matrix(1,2).im-exp_matrix(1,2).im)<1e-3 &&
      //[CO20230319 - not sure what this is?]      abs(_SU2_matrix(2,1).re-exp_matrix(2,1).re)<1e-3 && abs(_SU2_matrix(2,1).im-exp_matrix(2,1).im)<1e-3 &&
      //[CO20230319 - not sure what this is?]      abs(_SU2_matrix(2,2).re-exp_matrix(2,2).re)<1e-3 && abs(_SU2_matrix(2,2).im-exp_matrix(2,2).im)<1e-3)){
      //[CO20230319 - not sure what this is?]  //DX20200217 - warning to error
      //[CO20230319 - not sure what this is?]  stringstream message;
      //[CO20230319 - not sure what this is?]  message << "Lie algebra does not connect back to Lie group, i.e., exp(theta*su(2)) != SU(2):" << endl;
      //[CO20230319 - not sure what this is?]  message << "SU2:" << _SU2_matrix(1,1) << " " << _SU2_matrix(1,2) << " " << _SU2_matrix(2,1) << " " << _SU2_matrix(2,2) << endl;
      //[CO20230319 - not sure what this is?]  message << "exp_matrix:" << exp_matrix(1,1) << " " << exp_matrix(1,2) << " " << exp_matrix(2,1) << " " << exp_matrix(2,2) << endl;
      //[CO20230319 - not sure what this is?]  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
      //[CO20230319 - not sure what this is?]}
    }
    return true;
  }
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------
// Function SYM::CalculatePointGroupCrystal
// ----------------------------------------------------------------------------------------------------------------------------------------------------
// this function returns the PointGroup of the Crystal by searching through the Hermann-Mauguin combinations
//DX+CO START
namespace SYM {
  bool CalculatePointGroupCrystal(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss, string format) {
    double _eps_=AUROSTD_NAN; //DX =_EPS_;
    if(a.sym_eps!=AUROSTD_NAN){ //Tolerance came from user or was calculated
      _eps_=a.sym_eps;
    }
    else {
      _eps_=defaultTolerance(a);
      // Calculate point group/space group to find correct tolerance for system
      //SYM::CalculatePointGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss);       
      //_eps_=a.sym_eps;
    }
    return SYM::CalculatePointGroupCrystal(FileMESSAGE,a,aflags,_write_,osswrite,oss,_eps_,format);}
} // namespace SYM
//DX+CO END

//DX20170814 - Split up function so we can call look-up table separately - START
namespace SYM {
  bool CalculatePointGroupCrystal(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format) {
    //return CalculatePointGroupCrystal_20160801(FileMESSAGE,a,aflags,_write_,osswrite,oss,_eps_,format);
    return CalculatePointGroupCrystal_20170814(FileMESSAGE,a,aflags,_write_,osswrite,oss,_eps_,format);
  }

  bool CalculatePointGroupCrystal_20170814(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_, string format) {
    // Obtain the structure tolerance
    a.sym_eps=_eps_; //DX

    ostringstream aus;
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    DEBUG_SYMMETRY=DEBUG_SYMMETRY || LDEBUG;    
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUP_CRYSTAL Symmetry: BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);

    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal [0]" << endl;
    // this function returns the PointGroup of the Crystal by searching through the Hermann-Mauguin combinations

    bool Krun=TRUE;
    xvector<double> ctau(3),ftau(3),ctrasl(3),ftrasl(3);  // translation
    xmatrix<double> Uf(3,3),Uc(3,3);                      // matrices
    // ---------------------------------------------------------------------------
    // create placeholder for atoms_map for AddSymmetryToStructure() (no atom mappings for pgroups)
    std::vector<int> basis_atoms_map(a.atoms.size());           // will map each on each
    for(uint i=0;i<a.atoms.size();i++) basis_atoms_map[i]=i;    // identically map each over each
    std::vector<int> basis_types_map(a.atoms.size());           // will map each on each
    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal [0d]" << endl;
    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal a.atoms.size()=" << a.atoms.size() << endl;
    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal a.species.size()=" << a.species.size() << endl;
    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal basis_types_map.size()=" << basis_types_map.size() << endl;
    for(uint i=0;i<a.atoms.size();i++) basis_types_map[i]=a.atoms[i].type;    // identically map each over each
    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal [0e]" << endl;
    string message="PGROUP_XTAL";

    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal [1]" << endl;

    if(a.pgroup_calculated==FALSE) Krun=Krun && SYM::CalculatePointGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss);        // NEED POINT GROUP
    if(a.fgroup_calculated==FALSE) Krun=Krun && SYM::CalculateFactorGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss);       // NEED FACTOR GROUP

    // this one removes the double matrices so it takes care of non primitive systems
    for(uint ip=0;ip<a.fgroup.size();ip++) {
      clear(ctau);clear(ftau);clear(ctrasl);clear(ftrasl);
      Uf=a.fgroup[ip].Uf;Uc=a.fgroup[ip].Uc;ftau=a.fgroup[ip].ftau;ctau=a.fgroup[ip].ctau;
      bool sym_found=FALSE;
      for(uint ii=0;ii<a.pgroup_xtal.size()&&!sym_found;ii++)
        sym_found=(identical(Uf,a.pgroup_xtal[ii].Uf));     // look in all the list of operations  //DX20171207 - Use xmatrix identical eps
      if(sym_found==FALSE) {                                    // new operation, generate and save it
        SYM::AddSymmetryToStructure(a,Uc,Uf,ctau,ftau,ctrasl,ftrasl,basis_atoms_map,basis_types_map,true,_PGROUP_XTAL_,FALSE); 
      }
    }
    a.pgroup_xtal_calculated=TRUE;

    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal [2]" << endl;

    // PGROUP_XTAL
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " " << SEPARATION_LINE_DASH_SHORT << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    //cerr << "********************************************************" << endl;
    //cerr << "a.pgroup_xtal.size(): " << a.pgroup_xtal.size() << endl;
    for(uint kk=1;kk<=a.pgroup_xtal.size();kk++) {  // shift by 1 to use [kk-1]...
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " " ;
      aus << a.pgroup_xtal[kk-1].str_type << " theta=";
      if(a.pgroup_xtal[kk-1].angle<100) aus << " ";
      if(a.pgroup_xtal[kk-1].angle<10)  aus << " ";
      aus << a.pgroup_xtal[kk-1].angle << " " << " r=(" << a.pgroup_xtal[kk-1].axis << ")";
      aus << "    HM=" <<  aurostd::PaddedPRE(a.pgroup_xtal[kk-1].str_Hermann_Mauguin,2," ");     // remember vectors start from 0
      aus << "    S=" <<  aurostd::PaddedPRE(a.pgroup_xtal[kk-1].str_Schoenflies,2," ");          // remember vectors start from 0
      aus    << endl;  // remember vectors start from 0
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    }
    return PointGroupLookUpTable(FileMESSAGE,a,aflags,_write_,osswrite,oss,format);
  }
} 

// ----------------------------------------------------------------------------------------------------------------------------------------------------
// Function SYM::CalculatePointGroupKPatterson()
// ----------------------------------------------------------------------------------------------------------------------------------------------------
// Returns the point group (in reciprocal space) of the Patterson function (i.e., Patterson symmetry)
// Patterson symmetry represents the vector set symmetry of a crystal (see ITC-A pg. 19)
// The Patterson symmetry is 
//   1) symmorphic : only point group symmetries, no translations and
//   2) centrosymmetric : contains inversion center
// The Patterson symmetry is the same as the Laue symmetry (see GetLaueLabel() in aflow_xatom.cpp), 
// but this function explicity calculates the symmetry operators (matrices, axis-angle, etc.) for the particular representation.
// In the AFLOW symmetry workflow, this analysis is calculated by 
//   1) removing translations from the factor group and adding inversion symmetry or
//   2) adding inversion symmetry to pgroupk_xtal
// NOTE: This analysis is for reciprocal space only, so Uf is with respect to the klattice (reciprocal lattice)
namespace SYM {  
  bool CalculatePointGroupKPatterson(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss, string format) {
    double _eps_=AUROSTD_NAN; //DX =_EPS_;
    if(a.sym_eps!=AUROSTD_NAN){ //Tolerance came from user or was calculated
      _eps_=a.sym_eps;
    }
    else {
      _eps_=defaultTolerance(a);
    }
    return SYM::CalculatePointGroupKPatterson(FileMESSAGE,a,aflags,_write_,osswrite,oss,_eps_,format);}
} // namespace SYM

namespace SYM {  
  bool CalculatePointGroupKPatterson(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_, string format) {

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    DEBUG_SYMMETRY=(DEBUG_SYMMETRY || LDEBUG);    

    ostringstream aus;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUPK_PATTERSON Symmetry: BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);

    // ---------------------------------------------------------------------------
    // obtain the structure tolerance
    a.sym_eps=_eps_; 

    if(DEBUG_SYMMETRY) cerr << "DEBUG: " << __AFLOW_FUNC__ << " [0]" << endl;

    bool Krun=TRUE;
    xvector<double> ctau(3),ftau(3),ctrasl(3),ftrasl(3);        // translation
    clear(ctau);clear(ftau);clear(ctrasl);clear(ftrasl);
    xmatrix<double> Uf(3,3),Uc(3,3);                            // matrices
    // ---------------------------------------------------------------------------
    // create placeholder for atoms_map for AddSymmetryToStructure() (no atom mappings for pgroups)
    std::vector<int> basis_atoms_map(a.atoms.size());           // will map each on each
    for(uint i=0;i<a.atoms.size();i++) basis_atoms_map[i]=i;    // identically map each over each
    std::vector<int> basis_types_map(a.atoms.size());           // will map each on each
    if(DEBUG_SYMMETRY) cerr << "DEBUG: " << __AFLOW_FUNC__ << " [0d]" << endl;
    if(DEBUG_SYMMETRY) cerr << "DEBUG: " << __AFLOW_FUNC__ << " a.atoms.size()=" << a.atoms.size() << endl;
    if(DEBUG_SYMMETRY) cerr << "DEBUG: " << __AFLOW_FUNC__ << " a.species.size()=" << a.species.size() << endl;
    if(DEBUG_SYMMETRY) cerr << "DEBUG: " << __AFLOW_FUNC__ << " basis_types_map.size()=" << basis_types_map.size() << endl;
    for(uint i=0;i<a.atoms.size();i++) basis_types_map[i]=a.atoms[i].type;    // identically map each over each
    if(DEBUG_SYMMETRY) cerr << "DEBUG: " << __AFLOW_FUNC__ << " [0e]" << endl;
    string message="PGROUPK_PATTERSON";

    if(DEBUG_SYMMETRY) cerr << "DEBUG: " << __AFLOW_FUNC__ << " [1]" << endl;

    if(a.pgroup_calculated==FALSE) Krun=Krun && SYM::CalculatePointGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss);        // NEED POINT GROUP
    if(a.fgroup_calculated==FALSE) Krun=Krun && SYM::CalculateFactorGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss);       // NEED FACTOR GROUP
    if(a.pgroup_xtal_calculated==FALSE) Krun=Krun && SYM::CalculatePointGroupCrystal(FileMESSAGE,a,aflags,_write_,osswrite,oss); // NEED POINT GROUP CRYSTAL
    if(a.pgroupk_xtal_calculated==FALSE) Krun=Krun && SYM::CalculatePointGroupKCrystal(FileMESSAGE,a,aflags,_write_,osswrite,oss); // NEED POINT GROUP KCRYSTAL

    xmatrix<double> inversion_symmetry_matrix = -aurostd::identity((double)0,3);

    // ---------------------------------------------------------------------------
    // check if pgroup_xtal contains inversion symmetry already
    if(DEBUG_SYMMETRY) cerr << "DEBUG: " << __AFLOW_FUNC__ << " [2a] Check if pgroup_xtal contains inverison symmetry (i.e., centrosymmetric) " << endl;
    bool contains_inversion = false;
    for(uint ip=0;ip<a.pgroupk_xtal.size()&&!contains_inversion;ip++) {
      if(aurostd::identical(inversion_symmetry_matrix,a.pgroupk_xtal[ip].Uc)){
        contains_inversion = true;
      }
    }

    // ---------------------------------------------------------------------------
    // copy pgroup_xtal to pgroupk_Patterson; equivalent 
    if(contains_inversion){
      if(DEBUG_SYMMETRY) cerr << "DEBUG: " << __AFLOW_FUNC__ << " [3a] point group is already centrosymmetric, same as POINT GROUP CRYSTAL " << endl;
      _sym_op symop;
      for(uint k=0;k<a.pgroupk_xtal.size();k++) {
        SYM::AddSymmetryToStructure(
            a,
            a.pgroupk_xtal[k].Uc,
            a.pgroupk_xtal[k].Uf,
            a.pgroupk_xtal[k].ctau,
            a.pgroupk_xtal[k].ftau,
            a.pgroupk_xtal[k].ctrasl,
            a.pgroupk_xtal[k].ftrasl,
            a.pgroupk_xtal[k].basis_atoms_map,
            a.pgroupk_xtal[k].basis_types_map,
            false,
            _PGROUPK_PATTERSON_);
      }
    }
    else{
      // ---------------------------------------------------------------------------
      // calculation modes (default: 1)
      //   1) add inversion symmetries to pgroup_xtal (fast)
      //   2) explicitly calculate symmetry of vector set (good for verification/debugging) [TO-DO]
      uint calculation_mode = 1;

      // ---------------------------------------------------------------------------
      // add inversion symmetry to pgroup_xtal
      if(calculation_mode==1){
        if(DEBUG_SYMMETRY) cerr << "DEBUG: " << __AFLOW_FUNC__ << " [3b] calculation_mode: add inversion symmetry to pgroup_xtal (fast, default) " << endl;
        std::vector<_sym_op> pgroupk_Patterson;                             // rotations/inversions operations

        // ---------------------------------------------------------------------------
        // create Uc and Uf versions of inverse
        xmatrix<double> Uc_inv=inversion_symmetry_matrix;
        xmatrix<double> Uf_inv=inverse(trasp(a.klattice))*Uc_inv*trasp(a.klattice); // i.e., c2f*Uc_inv*f2c for klattice

        for(uint ip=0;ip<a.pgroupk_xtal.size();ip++) {
          Uf=a.pgroupk_xtal[ip].Uf;Uc=a.pgroupk_xtal[ip].Uc;
          bool sym_found=FALSE;
          for(uint ii=0;ii<a.pgroupk_Patterson.size()&&!sym_found;ii++){
            sym_found=(identical(Uf,a.pgroupk_Patterson[ii].Uf));     // look in all the list of operations  //DX20171207 - Use xmatrix identical eps
          }
          if(sym_found==FALSE) {                                    // new operation, generate and save it
            SYM::AddSymmetryToStructure(a,Uc,Uf,ctau,ftau,ctrasl,ftrasl,basis_atoms_map,basis_types_map,true,_PGROUPK_PATTERSON_,FALSE); 
          }
          Uf=Uf_inv*a.pgroupk_xtal[ip].Uf;Uc=Uc_inv*a.pgroupk_xtal[ip].Uc; //multiply by inversion matrix
          sym_found=FALSE;
          for(uint ii=0;ii<a.pgroupk_Patterson.size()&&!sym_found;ii++){
            sym_found=(identical(Uf,a.pgroupk_Patterson[ii].Uf));     // look in all the list of operations  //DX20171207 - Use xmatrix identical eps
          }
          if(sym_found==FALSE) {                                    // new operation, generate and save it
            SYM::AddSymmetryToStructure(a,Uc,Uf,ctau,ftau,ctrasl,ftrasl,basis_atoms_map,basis_types_map,true,_PGROUPK_PATTERSON_,FALSE); 
          }
        }
        a.pgroupk_Patterson_calculated=TRUE;
      }
      // ---------------------------------------------------------------------------
      // calculate symmetry of vector set [TO-DO]
      else if(calculation_mode==2){
        if(DEBUG_SYMMETRY) cerr << "DEBUG: " << __AFLOW_FUNC__ << " [3c] calculation_mode: explicitly calculate symmetry of vector set " << endl;
        // TO-DO
      }
    }
    if(DEBUG_SYMMETRY) cerr << "DEBUG: " << __AFLOW_FUNC__ << " [4]" << endl;

    // ---------------------------------------------------------------------------
    // PGROUPK_PATTERSON
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " " << SEPARATION_LINE_DASH_SHORT << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    for(uint kk=1;kk<=a.pgroupk_Patterson.size();kk++) {  // shift by 1 to use [kk-1]...
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " " ;
      aus << a.pgroupk_Patterson[kk-1].str_type << " theta=";
      if(a.pgroupk_Patterson[kk-1].angle<100) aus << " ";
      if(a.pgroupk_Patterson[kk-1].angle<10)  aus << " ";
      aus << a.pgroupk_Patterson[kk-1].angle << " " << " r=(" << a.pgroupk_Patterson[kk-1].axis << ")";
      aus << "    HM=" <<  aurostd::PaddedPRE(a.pgroupk_Patterson[kk-1].str_Hermann_Mauguin,2," ");     // remember vectors start from 0
      aus << "    S=" <<  aurostd::PaddedPRE(a.pgroupk_Patterson[kk-1].str_Schoenflies,2," ");          // remember vectors start from 0
      aus    << endl;  // remember vectors start from 0
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    }
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: unique point group operations " << a.pgroupk_Patterson.size() << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: END " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    if(_write_) Krun=Krun && KBIN_SymmetryWrite(FileMESSAGE,a,aflags,_PGROUPK_PATTERSON_,osswrite,oss,format);
    string pgname = "";  
    string operations = "";
    bool point_group_valid = PointGroupMap(a, pgname, operations, _PGROUPK_PATTERSON_);

    // ---------------------------------------------------------------------------
    // check Patterson point group label against Laue label (should be the same)
    string Laue_symbol = GetLaueLabel(a.point_group_Hermann_Mauguin);
    if(pgname != Laue_symbol){
      cerr << "WARNING:: Patterson symmetry label does not match the Laue symbol:" << endl;
      cerr << "Patterson point group: " << pgname << endl;
      cerr << "Laue point group: " << Laue_symbol << endl;
    }
    return (point_group_valid && (pgname==Laue_symbol));

  }
}


namespace SYM {
  bool PointGroupMap(xstructure& a, string& pgname, string& operations, char group){
    // LOOK up table for the point group operations of the crystal
    vector<string> vops;vops.clear();
    if(group==_PGROUP_){
      for(uint kk=1;kk<=a.pgroup.size();kk++) {  // shift by 1 to use [kk-1]...
        vops.push_back(a.pgroup[kk-1].str_Hermann_Mauguin);
        operations=operations+a.pgroup[kk-1].str_Hermann_Mauguin;
        if(kk<a.pgroup.size()) operations=operations+" ";
      }
    }
    if(group==_PGROUP_XTAL_){
      for(uint kk=1;kk<=a.pgroup_xtal.size();kk++) {  // shift by 1 to use [kk-1]...
        vops.push_back(a.pgroup_xtal[kk-1].str_Hermann_Mauguin);
        operations=operations+a.pgroup_xtal[kk-1].str_Hermann_Mauguin;
        if(kk<a.pgroup_xtal.size()) operations=operations+" ";
      }
    }
    if(group==_PGROUPK_PATTERSON_){
      for(uint kk=1;kk<=a.pgroupk_Patterson.size();kk++) {  // shift by 1 to use [kk-1]...
        vops.push_back(a.pgroupk_Patterson[kk-1].str_Hermann_Mauguin);
        operations=operations+a.pgroupk_Patterson[kk-1].str_Hermann_Mauguin;
        if(kk<a.pgroupk_Patterson.size()) operations=operations+" ";
      }
    }
    bool pg_found=FALSE;

    vector<string> pgf,pgn;
    // TRICLINIC
    pgf.push_back("1");pgn.push_back("1"); // 1
    pgf.push_back("1 -1");pgn.push_back("-1"); // 2
    // MONOCLINIC
    pgf.push_back("1 2");pgn.push_back("2"); // 3-5
    pgf.push_back("1 m");pgn.push_back("m"); // 6-9
    pgf.push_back("1 2 m -1");pgn.push_back("2/m"); // 10-15
    // ORTHORHOMBIC
    pgf.push_back("1 2 2 2");pgn.push_back("222"); // 16-24
    pgf.push_back("1 m m 2");pgn.push_back("mm2"); // 25-46
    pgf.push_back("1 m m 2 m 2 2 -1");pgn.push_back("mmm"); // 47-74
    // TETRAGONAL
    pgf.push_back("1 4 2 4");pgn.push_back("4"); // 75-80
    pgf.push_back("1 2 -4 -4");pgn.push_back("-4"); // 81-82
    pgf.push_back("1 4 2 4 m -4 -1 -4");pgn.push_back("4/m"); // 83-88
    pgf.push_back("1 4 2 4 2 2 2 2");pgn.push_back("422"); // 89-98
    pgf.push_back("1 m -2 4 m 2 4 -2");pgn.push_back("4mm"); // 99-110
    pgf.push_back("1 m 2 -4 m 2 -4 2");pgn.push_back("-42m"); // 111-122
    pgf.push_back("1 m -2 4 m 2 4 -2 m 2 2 -4 2 -1 -4 2");pgn.push_back("4/mmm"); // 123-142
    // HEXAGONAL-TRIGONAL
    pgf.push_back("1 3 3");pgn.push_back("3"); // 143-146
    pgf.push_back("1 3 3 -1 -3 -3");pgn.push_back("-3"); // 147-148
    pgf.push_back("1 2 2 3 3 2");pgn.push_back("32"); // 149-155
    pgf.push_back("1 3 3 m -2 -2");pgn.push_back("3m"); // 156-161
    pgf.push_back("1 m -2 3 3 -2 -1 2 2 -3 -3 2");pgn.push_back("-3m"); // 162-167

    // HEXAGONAL-HEXAGONAL
    pgf.push_back("1 3 3 2 6 6");pgn.push_back("6"); // 168-173
    pgf.push_back("1 3 3 m -6 -6");pgn.push_back("-6"); // 174
    pgf.push_back("1 3 3 2 6 6 m -6 -6 -1 -3 -3");pgn.push_back("6/m"); // 175-176
    pgf.push_back("1 3 3 2 6 6 2 2 2 2 2 2");pgn.push_back("622"); // 177-182
    pgf.push_back("1 -2 m 3 3 -2 2 -2 m 6 6 -2");pgn.push_back("6mm"); // 183-186
    pgf.push_back("1 m -2 2 m 2 3 -6 3 -6 -2 2");pgn.push_back("-6m2"); // 187-190
    pgf.push_back("1 m -2 2 m 2 3 -6 3 -6 -2 2 2 -1 -2 2 m 2 6 -3 6 -3 -2 2");pgn.push_back("6/mmm"); // 	191-194
    // CUBIC
    pgf.push_back("1 2 3 3 3 3 2 2 3 3 3 3");pgn.push_back("23"); //  195-199
    pgf.push_back("1 3 3 2 3 3 3 3 2 3 3 2 -1 -3 -3 m -3 -3 -3 -3 m -3 -3 m");pgn.push_back("m-3"); // 200-206
    pgf.push_back("1 4 2 4 2 3 4 3 3 4 3 2 2 2 2 2 4 3 2 3 3 4 3 2");pgn.push_back("432"); // 207-214
    pgf.push_back("1 -2 -2 3 3 -2 2 -2 -4 3 3 -4 -4 3 2 -4 3 -4 -4 2 -2 3 3 -2");pgn.push_back("-43m"); // 215-220
    pgf.push_back("1 -2 2 m -2 3 -3 4 3 -2 4 -3 2 -2 -1 2 -3 -3 -2 3 4 4 3 -2 -3 -4 3 2 2 3 -4 -3 2 -4 -3 3 4 m 2 -4 -3 2 3 -4 m 4 -4 2");pgn.push_back("m-3m"); // 221-230
    //  pgf.push_back("1 -2 -2 3 3 -2 -2 2 3 -4 -4 3 3 -4 -2 3 2 -4 -4 3 3 -2 -4 2 -1 2 2 -3 -3 2 2 m -3 4 4 -3 -3 4 2 -3 -2 4 4 -3 -3 2 4 -2");pgn.push_back("m-3m"); // 221-230
    // done
    for(uint i=0;i<pgf.size();i++) aurostd::StringSubst(pgf.at(i),"-2","m"); // -2 are mirrors orthogonal to the axis....

    // -------------------------------------------------------------------- scanning  ORDER=48
    // find point group
    aurostd::StringSubst(operations,"-2","m");

    for(uint i=0;i<pgf.size()&&!pg_found;i++) {
      if(MapOperations(operations,pgf.at(i))) {pg_found=TRUE;pgname=pgn.at(i);}
    }

    return pg_found;
  }
}

namespace SYM {
  bool PointGroupLookUpTable(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format){
    ostringstream aus;
    bool LDEBUG=(FALSE || XHOST.DEBUG || DEBUG_SYMMETRY);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    bool Krun=TRUE;
    string message="PGROUP_XTAL";

    // LOOK up table for the point group operations of the crystal
    //DX20170916 [OBSOLETE]    vector<string> vops;vops.clear();
    //DX20170916 [OBSOLETE]   string operations="";
    //DX20170916 [OBSOLETE]    for(uint kk=1;kk<=a.pgroup_xtal.size();kk++) {  // shift by 1 to use [kk-1]...
    //DX20170916 [OBSOLETE]      vops.push_back(a.pgroup_xtal[kk-1].str_Hermann_Mauguin);
    //DX20170916 [OBSOLETE]      operations=operations+a.pgroup_xtal[kk-1].str_Hermann_Mauguin;
    //DX20170916 [OBSOLETE]      if(kk<a.pgroup_xtal.size()) operations=operations+" ";
    //DX20170916 [OBSOLETE]    }

    a.crystal_family="";a.crystal_system="";a.point_group_crystal_class="";
    a.point_group_Shoenflies="";a.point_group_Hermann_Mauguin="";a.point_group_orbifold="";
    a.point_group_type="";a.point_group_order="";a.point_group_structure="";
    //DX20170906 [OBSOLETE] bool pg_found=FALSE;
    string pgname = ""; //DX20170906
    string operations = ""; //DX20170906

    //DX20170916 [OBSOLETE]    vector<string> pgf,pgn;
    //DX20170916 [OBSOLETE]    // TRICLINIC
    //DX20170916 [OBSOLETE]    pgf.push_back("1");pgn.push_back("1"); // 1
    //DX20170916 [OBSOLETE]    pgf.push_back("1 -1");pgn.push_back("-1"); // 2
    //DX20170916 [OBSOLETE]    // MONOCLINIC
    //DX20170916 [OBSOLETE]    pgf.push_back("1 2");pgn.push_back("2"); // 3-5
    //DX20170916 [OBSOLETE]    pgf.push_back("1 m");pgn.push_back("m"); // 6-9
    //DX20170916 [OBSOLETE]    pgf.push_back("1 2 m -1");pgn.push_back("2/m"); // 10-15
    //DX20170916 [OBSOLETE]    // ORTHORHOMBIC
    //DX20170916 [OBSOLETE]    pgf.push_back("1 2 2 2");pgn.push_back("222"); // 16-24
    //DX20170916 [OBSOLETE]    pgf.push_back("1 m m 2");pgn.push_back("mm2"); // 25-46
    //DX20170916 [OBSOLETE]    pgf.push_back("1 m m 2 m 2 2 -1");pgn.push_back("mmm"); // 47-74
    //DX20170916 [OBSOLETE]    // TETRAGONAL
    //DX20170916 [OBSOLETE]    pgf.push_back("1 4 2 4");pgn.push_back("4"); // 75-80
    //DX20170916 [OBSOLETE]    pgf.push_back("1 2 -4 -4");pgn.push_back("-4"); // 81-82
    //DX20170916 [OBSOLETE]    pgf.push_back("1 4 2 4 m -4 -1 -4");pgn.push_back("4/m"); // 83-88
    //DX20170916 [OBSOLETE]    pgf.push_back("1 4 2 4 2 2 2 2");pgn.push_back("422"); // 89-98
    //DX20170916 [OBSOLETE]    pgf.push_back("1 m -2 4 m 2 4 -2");pgn.push_back("4mm"); // 99-110
    //DX20170916 [OBSOLETE]    pgf.push_back("1 m 2 -4 m 2 -4 2");pgn.push_back("-42m"); // 111-122
    //DX20170916 [OBSOLETE]    pgf.push_back("1 m -2 4 m 2 4 -2 m 2 2 -4 2 -1 -4 2");pgn.push_back("4/mmm"); // 123-142
    //DX20170916 [OBSOLETE]    // HEXAGONAL-TRIGONAL
    //DX20170916 [OBSOLETE]    pgf.push_back("1 3 3");pgn.push_back("3"); // 143-146
    //DX20170916 [OBSOLETE]    pgf.push_back("1 3 3 -1 -3 -3");pgn.push_back("-3"); // 147-148
    //DX20170916 [OBSOLETE]    pgf.push_back("1 2 2 3 3 2");pgn.push_back("32"); // 149-155
    //DX20170916 [OBSOLETE]    pgf.push_back("1 3 3 m -2 -2");pgn.push_back("3m"); // 156-161
    //DX20170916 [OBSOLETE]    pgf.push_back("1 m -2 3 3 -2 -1 2 2 -3 -3 2");pgn.push_back("-3m"); // 162-167
    //DX20170916 [OBSOLETE]
    //DX20170916 [OBSOLETE]    // HEXAGONAL-HEXAGONAL
    //DX20170916 [OBSOLETE]    pgf.push_back("1 3 3 2 6 6");pgn.push_back("6"); // 168-173
    //DX20170916 [OBSOLETE]    pgf.push_back("1 3 3 m -6 -6");pgn.push_back("-6"); // 174
    //DX20170916 [OBSOLETE]    pgf.push_back("1 3 3 2 6 6 m -6 -6 -1 -3 -3");pgn.push_back("6/m"); // 175-176
    //DX20170916 [OBSOLETE]    pgf.push_back("1 3 3 2 6 6 2 2 2 2 2 2");pgn.push_back("622"); // 177-182
    //DX20170916 [OBSOLETE]    pgf.push_back("1 -2 m 3 3 -2 2 -2 m 6 6 -2");pgn.push_back("6mm"); // 183-186
    //DX20170916 [OBSOLETE]    pgf.push_back("1 m -2 2 m 2 3 -6 3 -6 -2 2");pgn.push_back("-6m2"); // 187-190
    //DX20170916 [OBSOLETE]    pgf.push_back("1 m -2 2 m 2 3 -6 3 -6 -2 2 2 -1 -2 2 m 2 6 -3 6 -3 -2 2");pgn.push_back("6/mmm"); // 	191-194
    //DX20170916 [OBSOLETE]    // CUBIC
    //DX20170916 [OBSOLETE]    pgf.push_back("1 2 3 3 3 3 2 2 3 3 3 3");pgn.push_back("23"); //  195-199
    //DX20170916 [OBSOLETE]    pgf.push_back("1 3 3 2 3 3 3 3 2 3 3 2 -1 -3 -3 m -3 -3 -3 -3 m -3 -3 m");pgn.push_back("m-3"); // 200-206
    //DX20170916 [OBSOLETE]    pgf.push_back("1 4 2 4 2 3 4 3 3 4 3 2 2 2 2 2 4 3 2 3 3 4 3 2");pgn.push_back("432"); // 207-214
    //DX20170916 [OBSOLETE]    pgf.push_back("1 -2 -2 3 3 -2 2 -2 -4 3 3 -4 -4 3 2 -4 3 -4 -4 2 -2 3 3 -2");pgn.push_back("-43m"); // 215-220
    //DX20170916 [OBSOLETE]    pgf.push_back("1 -2 2 m -2 3 -3 4 3 -2 4 -3 2 -2 -1 2 -3 -3 -2 3 4 4 3 -2 -3 -4 3 2 2 3 -4 -3 2 -4 -3 3 4 m 2 -4 -3 2 3 -4 m 4 -4 2");pgn.push_back("m-3m"); // 221-230
    //DX20170916 [OBSOLETE]    //  pgf.push_back("1 -2 -2 3 3 -2 -2 2 3 -4 -4 3 3 -4 -2 3 2 -4 -4 3 3 -2 -4 2 -1 2 2 -3 -3 2 2 m -3 4 4 -3 -3 4 2 -3 -2 4 4 -3 -3 2 4 -2");pgn.push_back("m-3m"); // 221-230
    //DX20170916 [OBSOLETE]    // done
    //DX20170916 [OBSOLETE]    for(uint i=0;i<pgf.size();i++) aurostd::StringSubst(pgf.at(i),"-2","m"); // -2 are mirrors orthogonal to the axis....
    //DX20170916 [OBSOLETE]
    //DX20170916 [OBSOLETE]    // -------------------------------------------------------------------- scanning  ORDER=48
    //DX20170916 [OBSOLETE]    // find point group
    //DX20170916 [OBSOLETE]    string pgname;
    //DX20170916 [OBSOLETE]    aurostd::StringSubst(operations,"-2","m");
    //DX20170916 [OBSOLETE]
    //DX20170916 [OBSOLETE]    for(uint i=0;i<pgf.size()&&!pg_found;i++) {
    //DX20170916 [OBSOLETE]      if(MapOperations(operations,pgf.at(i))) {pg_found=TRUE;pgname=pgn.at(i);}
    //DX20170916 [OBSOLETE]    }
    //DX20170916 [OBSOLETE]    //  cerr << pg_found << endl;
    //DX20170916 [OBSOLETE]    //   cerr << pgname << endl;

    bool pg_found = PointGroupMap(a,pgname,operations,_PGROUP_XTAL_); //DX20170906
    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " point group symbol: " << pgname << endl; } //DX20210327
    // -------------------------------------------------------------------- scanning  ORDER=48
    if(pgname=="m-3m") {
      pg_found=TRUE;
      a.crystal_family="cubic";a.crystal_system="cubic";a.point_group_crystal_class="hexoctahedral";
      a.point_group_Shoenflies="O_h";a.point_group_Hermann_Mauguin="m-3m";a.point_group_orbifold="*432";
      a.point_group_type="centrosymmetric";a.point_group_order="48";a.point_group_structure="2 x symmetric";
    }
    // -------------------------------------------------------------------- scanning  ORDER=24
    if(pgname=="-43m") {
      pg_found=TRUE;
      a.crystal_family="cubic";a.crystal_system="cubic";a.point_group_crystal_class="tetrahedral";
      a.point_group_Shoenflies="T_d";a.point_group_Hermann_Mauguin="-43m";a.point_group_orbifold="*332";
      a.point_group_type="none";a.point_group_order="24";a.point_group_structure="symmetric"; //DX20180827 - point group type changed from - to none
    }
    if(pgname=="432") {
      pg_found=TRUE;
      a.crystal_family="cubic";a.crystal_system="cubic";a.point_group_crystal_class="gyroidal";
      a.point_group_Shoenflies="O";a.point_group_Hermann_Mauguin="432";a.point_group_orbifold="432";
      a.point_group_type="enantiomorphic";a.point_group_order="24";a.point_group_structure="symmetric";
    }
    if(pgname=="m-3") {
      pg_found=TRUE;
      a.crystal_family="cubic";a.crystal_system="cubic";a.point_group_crystal_class="diploidal";
      a.point_group_Shoenflies="T_h";a.point_group_Hermann_Mauguin="m-3";a.point_group_orbifold="3*2";
      a.point_group_type="centrosymmetric";a.point_group_order="24";a.point_group_structure="2 x alternating";
    }
    if(pgname=="6/mmm") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="hexagonal";a.point_group_crystal_class="dihexagonal-dipyramidal";
      a.point_group_Shoenflies="D_6h";a.point_group_Hermann_Mauguin="6/mmm";a.point_group_orbifold="*622";
      a.point_group_type="centrosymmetric";a.point_group_order="24";a.point_group_structure="2 x dihedral";
    }
    // -------------------------------------------------------------------- scanning  ORDER=16
    if(pgname=="4/mmm") {
      pg_found=TRUE;
      a.crystal_family="tetragonal";a.crystal_system="tetragonal";a.point_group_crystal_class="ditetragonal-dipyramidal";
      a.point_group_Shoenflies="D_4h";a.point_group_Hermann_Mauguin="4/mmm";a.point_group_orbifold="*422";
      a.point_group_type="centrosymmetric";a.point_group_order="16";a.point_group_structure="2 x dihedral";
    }
    // -------------------------------------------------------------------- scanning  ORDER=12
    if(pgname=="23") {
      pg_found=TRUE;
      a.crystal_family="cubic";a.crystal_system="cubic";a.point_group_crystal_class="tetrahedral";
      a.point_group_Shoenflies="T";a.point_group_Hermann_Mauguin="23";a.point_group_orbifold="332";
      a.point_group_type="enantiomorphic";a.point_group_order="12";a.point_group_structure="alternating";
    }
    if(pgname=="-6m2" || pgname=="-62m") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="hexagonal";a.point_group_crystal_class="ditrigonal-dipyramidal";
      a.point_group_Shoenflies="D_3h";a.point_group_Hermann_Mauguin="-6m2";a.point_group_orbifold="*322";
      a.point_group_type="none";a.point_group_order="12";a.point_group_structure="dihedral"; //DX20180827 - point group type changed from - to none
    }
    if(pgname=="6mm") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="hexagonal";a.point_group_crystal_class="dihexagonal-pyramidal";
      a.point_group_Shoenflies="C_6v";a.point_group_Hermann_Mauguin="6mm";a.point_group_orbifold="*66";
      a.point_group_type="polar";a.point_group_order="12";a.point_group_structure="dihedral";
    }
    if(pgname=="622") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="hexagonal";a.point_group_crystal_class="hexagonal-trapezoidal";
      a.point_group_Shoenflies="D_6";a.point_group_Hermann_Mauguin="622";a.point_group_orbifold="622";
      a.point_group_type="enantiomorphic";a.point_group_order="12";a.point_group_structure="dihedral";
    }
    if(pgname=="6/m") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="hexagonal";a.point_group_crystal_class="hexagonal-dipyramidal";
      a.point_group_Shoenflies="C_6h";a.point_group_Hermann_Mauguin="6/m";a.point_group_orbifold="6*";
      a.point_group_type="centrosymmetric";a.point_group_order="12";a.point_group_structure="2 × cyclic";
    }
    if(pgname=="-3m" || pgname=="-3m1" || pgname=="-31m") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="trigonal";a.point_group_crystal_class="ditrigonal-scalahedral";
      a.point_group_Shoenflies="D_3d";a.point_group_Hermann_Mauguin="-3m";a.point_group_orbifold="2*3";
      a.point_group_type="centrosymmetric";a.point_group_order="12";a.point_group_structure="dihedral";
    }
    // -------------------------------------------------------------------- scanning  ORDER=8
    if(pgname=="-42m" || pgname=="-4m2") {
      pg_found=TRUE;
      a.crystal_family="tetragonal";a.crystal_system="tetragonal";a.point_group_crystal_class="tetragonal-scalenoidal";
      a.point_group_Shoenflies="D_2d";a.point_group_Hermann_Mauguin="-42m";a.point_group_orbifold="2*2";
      a.point_group_type="none";a.point_group_order="8";a.point_group_structure="dihedral"; //DX20180827 - point group type changed from - to none
    }
    if(pgname=="4mm") {
      pg_found=TRUE;
      a.crystal_family="tetragonal";a.crystal_system="tetragonal";a.point_group_crystal_class="ditetragonal-pyramidal";
      a.point_group_Shoenflies="C_4v";a.point_group_Hermann_Mauguin="4mm";a.point_group_orbifold="*44";
      a.point_group_type="polar";a.point_group_order="8";a.point_group_structure="dihedral";
    }
    if(pgname=="422") {
      pg_found=TRUE;
      a.crystal_family="tetragonal";a.crystal_system="tetragonal";a.point_group_crystal_class="tetragonal-trapezoidal";
      a.point_group_Shoenflies="D_4";a.point_group_Hermann_Mauguin="422";a.point_group_orbifold="422";
      a.point_group_type="enantiomorphic";a.point_group_order="8";a.point_group_structure="dihedral";
    }
    if(pgname=="4/m") {
      pg_found=TRUE;
      a.crystal_family="tetragonal";a.crystal_system="tetragonal";a.point_group_crystal_class="tetragonal-dipyramidal";
      a.point_group_Shoenflies="C_4h";a.point_group_Hermann_Mauguin="4/m";a.point_group_orbifold="4*";
      a.point_group_type="centrosymmetric";a.point_group_order="8";a.point_group_structure="2 x cyclic";
    }
    if(pgname=="mmm") { //
      pg_found=TRUE;
      a.crystal_family="orthorhombic";a.crystal_system="orthorhombic";a.point_group_crystal_class="orthorhombic-bipyramidal";
      a.point_group_Shoenflies="D_2h";a.point_group_Hermann_Mauguin="mmm";a.point_group_orbifold="*222";
      a.point_group_type="centrosymmetric";a.point_group_order="8";a.point_group_structure="2 x dihedral";
    }
    // -------------------------------------------------------------------- scanning  ORDER=6
    if(pgname=="-6") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="hexagonal";a.point_group_crystal_class="trigonal-dipyramidal";
      a.point_group_Shoenflies="C_3h";a.point_group_Hermann_Mauguin="-6";a.point_group_orbifold="3*";
      a.point_group_type="none";a.point_group_order="6";a.point_group_structure="cyclic"; //DX20180827 - point group type changed from - to none
    }
    if(pgname=="6") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="hexagonal";a.point_group_crystal_class="hexagonal-pyramidal";
      a.point_group_Shoenflies="C_6";a.point_group_Hermann_Mauguin="6";a.point_group_orbifold="66";
      a.point_group_type="enantiomorphic polar";a.point_group_order="6";a.point_group_structure="cyclic";
    }
    if(pgname=="3m" || pgname=="3m1" || pgname=="31m") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="trigonal";a.point_group_crystal_class="ditrigonal-pyramidal";
      a.point_group_Shoenflies="C_3v";a.point_group_Hermann_Mauguin="3m";a.point_group_orbifold="*33";
      a.point_group_type="polar";a.point_group_order="6";a.point_group_structure="dihedral";
    }
    if(pgname=="32" || pgname=="321" || pgname=="312") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="trigonal";a.point_group_crystal_class="trigonal-trapezoidal";
      a.point_group_Shoenflies="D_3";a.point_group_Hermann_Mauguin="32";a.point_group_orbifold="322";
      a.point_group_type="enantiomorphic";a.point_group_order="6";a.point_group_structure="dihedral";
    }
    if(pgname=="-3") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="trigonal";a.point_group_crystal_class="rhombohedral";
      a.point_group_Shoenflies="S_6";a.point_group_Hermann_Mauguin="-3";a.point_group_orbifold="3x";
      a.point_group_type="centrosymmetric";a.point_group_order="6";a.point_group_structure="cyclic";
    }
    // -------------------------------------------------------------------- scanning  ORDER=4
    if(pgname=="-4") {
      pg_found=TRUE;
      a.crystal_family="tetragonal";a.crystal_system="tetragonal";a.point_group_crystal_class="tetragonal-disphenoidal";
      a.point_group_Shoenflies="S_4";a.point_group_Hermann_Mauguin="-4";a.point_group_orbifold="2x";
      a.point_group_type="none";a.point_group_order="4";a.point_group_structure="cyclic"; //DX20180827 - point group type changed from - to none
    }
    if(pgname=="4") {
      pg_found=TRUE;
      a.crystal_family="tetragonal";a.crystal_system="tetragonal";a.point_group_crystal_class="tetragonal-pyramidal";
      a.point_group_Shoenflies="C_4";a.point_group_Hermann_Mauguin="4";a.point_group_orbifold="44";
      a.point_group_type="enantiomorphic polar";a.point_group_order="4";a.point_group_structure="2 x cyclic";
    }
    if(pgname=="mm2") {
      pg_found=TRUE;
      a.crystal_family="orthorhombic";a.crystal_system="orthorhombic";a.point_group_crystal_class="orthorhombic-pyramidal";
      a.point_group_Shoenflies="C_2v";a.point_group_Hermann_Mauguin="mm2";a.point_group_orbifold="*22";
      a.point_group_type="polar";a.point_group_order="4";a.point_group_structure="dihedral";
    }
    if(pgname=="222") {
      pg_found=TRUE;
      a.crystal_family="orthorhombic";a.crystal_system="orthorhombic";a.point_group_crystal_class="orthorhombic-sphenoidal";
      a.point_group_Shoenflies="D_2";a.point_group_Hermann_Mauguin="222";a.point_group_orbifold="222";
      a.point_group_type="enantiomorphic";a.point_group_order="4";a.point_group_structure="dihedral";
    }
    if(pgname=="2/m") {
      pg_found=TRUE;
      a.crystal_family="monoclinic";a.crystal_system="monoclinic";a.point_group_crystal_class="monoclinic-prismatic";
      a.point_group_Shoenflies="C_2h";a.point_group_Hermann_Mauguin="2/m";a.point_group_orbifold="2*";
      a.point_group_type="centrosymmetric";a.point_group_order="4";a.point_group_structure="2 x cyclic";
    }
    // -------------------------------------------------------------------- scanning  ORDER=3
    if(pgname=="3") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="trigonal";a.point_group_crystal_class="trigonal-pyramidal";
      a.point_group_Shoenflies="C_3";a.point_group_Hermann_Mauguin="3";a.point_group_orbifold="33";
      a.point_group_type="enantiomorphic polar";a.point_group_order="3";a.point_group_structure="cyclic";
    }
    // -------------------------------------------------------------------- scanning  ORDER=2
    if(pgname=="m") {
      pg_found=TRUE;
      a.crystal_family="monoclinic";a.crystal_system="monoclinic";a.point_group_crystal_class="monoclinic-domatic";
      a.point_group_Shoenflies="C_S";a.point_group_Hermann_Mauguin="m";a.point_group_orbifold="1*";
      a.point_group_type="polar";a.point_group_order="2";a.point_group_structure="cyclic";
    }
    if(pgname=="2") {
      pg_found=TRUE;
      a.crystal_family="monoclinic";a.crystal_system="monoclinic";a.point_group_crystal_class="monoclinic-sphenoidal";
      a.point_group_Shoenflies="C_2";a.point_group_Hermann_Mauguin="2";a.point_group_orbifold="22";
      a.point_group_type="enantiomorphic polar";a.point_group_order="2";a.point_group_structure="cyclic";
    }
    if(pgname=="-1") {
      pg_found=TRUE;
      a.crystal_family="triclinic";a.crystal_system="triclinic";a.point_group_crystal_class="triclinic-pinacoidal";
      a.point_group_Shoenflies="C_i";a.point_group_Hermann_Mauguin="-1";a.point_group_orbifold="x";
      a.point_group_type="centrosymmetric";a.point_group_order="2";a.point_group_structure="cyclic";
    }
    // -------------------------------------------------------------------- scanning  ORDER=1
    if(pgname=="1") {
      pg_found=TRUE;
      a.crystal_family="triclinic";a.crystal_system="triclinic";a.point_group_crystal_class="triclinic-pedial";
      a.point_group_Shoenflies="C_1";a.point_group_Hermann_Mauguin="1";a.point_group_orbifold="11";
      a.point_group_type="enantiomorphic polar";a.point_group_order="1";a.point_group_structure="trivial";
    }


    // ------------------------------------------------------------------------------------------- PRINTING AND LEAVING
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message
      << " " << SEPARATION_LINE_DASH_SHORT << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message
      << " Symmetry: unique point group operations " << a.pgroup_xtal.size() << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: " << operations << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Crystal Family = " << a.crystal_family << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Crystal System = " << a.crystal_system << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Crystal Class = " << a.point_group_crystal_class << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Schoenflies notation = " << a.point_group_Shoenflies << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Hermann Mauguin notation = " << a.point_group_Hermann_Mauguin << endl; //"   - " << a.title << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Orbifold = " << a.point_group_orbifold << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Point Group Type = " << a.point_group_type << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Point Group Order = " << a.point_group_order << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Point Group Structure = " << a.point_group_structure << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " " << SEPARATION_LINE_DASH_SHORT << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: END " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message
      << " " << SEPARATION_LINE_DASH_SHORT << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    if(_write_) Krun=Krun && KBIN_SymmetryWrite(FileMESSAGE,a,aflags,_PGROUP_XTAL_,osswrite,oss,format);
    //DX20170906 [OBSOLETE] return Krun;
    return pg_found;
  }
} // namespace SYM
//DX20170814 - Split up function so we can call look-up table separately - END

namespace SYM {
  bool CalculatePointGroupCrystal_20160801(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_, string format) {
    // Obtain the structure tolerance
    a.sym_eps=_eps_; //DX

    ostringstream aus;
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    DEBUG_SYMMETRY=DEBUG_SYMMETRY || LDEBUG;    
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUP_CRYSTAL Symmetry: BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);

    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal [0]" << endl;
    // this function returns the PointGroup of the Crystal by searching through the Hermann-Mauguin combinations

    bool Krun=TRUE;
    xvector<double> ctau(3),ftau(3),ctrasl(3),ftrasl(3);  // translation
    xmatrix<double> Uf(3,3),Uc(3,3);                      // matrices
    std::vector<int> basis_atoms_map(a.atoms.size());           // will map each on each
    for(uint i=0;i<a.atoms.size();i++) basis_atoms_map[i]=i;    // identically map each over each
    std::vector<int> basis_types_map(a.atoms.size());           // will map each on each
    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal [0d]" << endl;
    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal a.atoms.size()=" << a.atoms.size() << endl;
    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal a.species.size()=" << a.species.size() << endl;
    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal basis_types_map.size()=" << basis_types_map.size() << endl;
    for(uint i=0;i<a.atoms.size();i++) basis_types_map[i]=a.atoms[i].type;    // identically map each over each
    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal [0e]" << endl;
    string message="PGROUP_XTAL";

    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal [1]" << endl;

    if(a.pgroup_calculated==FALSE) Krun=Krun && SYM::CalculatePointGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss);        // NEED POINT GROUP
    if(a.fgroup_calculated==FALSE) Krun=Krun && SYM::CalculateFactorGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss);       // NEED FACTOR GROUP

    // this one removes the double matrices so it takes care of non primitive systems
    for(uint ip=0;ip<a.fgroup.size();ip++) {
      clear(ctau);clear(ftau);clear(ctrasl);clear(ftrasl);
      Uf=a.fgroup[ip].Uf;Uc=a.fgroup[ip].Uc;ftau=a.fgroup[ip].ftau;ctau=a.fgroup[ip].ctau;
      bool sym_found=FALSE;
      for(uint ii=0;ii<a.pgroup_xtal.size()&&!sym_found;ii++)
        sym_found=(identical(Uf,a.pgroup_xtal[ii].Uf));     // look in all the list of operations  //DX20171207 - Use xmatrix identical eps
      if(sym_found==FALSE) {                                    // new operation, generate and save it
        SYM::AddSymmetryToStructure(a,Uc,Uf,ctau,ftau,ctrasl,ftrasl,basis_atoms_map,basis_types_map,true,_PGROUP_XTAL_,FALSE); 
      }
    }
    a.pgroup_xtal_calculated=TRUE;

    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal [2]" << endl;

    // PGROUP_XTAL
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " " << SEPARATION_LINE_DASH_SHORT << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    //cerr << "********************************************************" << endl;
    //cerr << "a.pgroup_xtal.size(): " << a.pgroup_xtal.size() << endl;
    for(uint kk=1;kk<=a.pgroup_xtal.size();kk++) {  // shift by 1 to use [kk-1]...
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " " ;
      aus << a.pgroup_xtal[kk-1].str_type << " theta=";
      if(a.pgroup_xtal[kk-1].angle<100) aus << " ";
      if(a.pgroup_xtal[kk-1].angle<10)  aus << " ";
      aus << a.pgroup_xtal[kk-1].angle << " " << " r=(" << a.pgroup_xtal[kk-1].axis << ")";
      aus << "    HM=" <<  aurostd::PaddedPRE(a.pgroup_xtal[kk-1].str_Hermann_Mauguin,2," ");     // remember vectors start from 0
      aus << "    S=" <<  aurostd::PaddedPRE(a.pgroup_xtal[kk-1].str_Schoenflies,2," ");          // remember vectors start from 0
      aus    << endl;  // remember vectors start from 0
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    }

    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal [3]" << endl;

    // LOOK up table for the point group operations of the crystal
    vector<string> vops;vops.clear();
    string operations="";
    for(uint kk=1;kk<=a.pgroup_xtal.size();kk++) {  // shift by 1 to use [kk-1]...
      vops.push_back(a.pgroup_xtal[kk-1].str_Hermann_Mauguin);
      operations=operations+a.pgroup_xtal[kk-1].str_Hermann_Mauguin;
      if(kk<a.pgroup_xtal.size()) operations=operations+" ";
    }
    a.crystal_family="";a.crystal_system="";a.point_group_crystal_class="";
    a.point_group_Shoenflies="";a.point_group_Hermann_Mauguin="";a.point_group_orbifold="";
    a.point_group_type="";a.point_group_order="";a.point_group_structure="";
    bool pg_found=FALSE;

    vector<string> pgf,pgn;
    // TRICLINIC
    pgf.push_back("1");pgn.push_back("1"); // 1
    pgf.push_back("1 -1");pgn.push_back("-1"); // 2
    // MONOCLINIC
    pgf.push_back("1 2");pgn.push_back("2"); // 3-5
    pgf.push_back("1 m");pgn.push_back("m"); // 6-9
    pgf.push_back("1 2 m -1");pgn.push_back("2/m"); // 10-15
    // ORTHORHOMBIC
    pgf.push_back("1 2 2 2");pgn.push_back("222"); // 16-24
    pgf.push_back("1 m m 2");pgn.push_back("mm2"); // 25-46
    pgf.push_back("1 m m 2 m 2 2 -1");pgn.push_back("mmm"); // 47-74
    // TETRAGONAL
    pgf.push_back("1 4 2 4");pgn.push_back("4"); // 75-80
    pgf.push_back("1 2 -4 -4");pgn.push_back("-4"); // 81-82
    pgf.push_back("1 4 2 4 m -4 -1 -4");pgn.push_back("4/m"); // 83-88
    pgf.push_back("1 4 2 4 2 2 2 2");pgn.push_back("422"); // 89-98
    pgf.push_back("1 m -2 4 m 2 4 -2");pgn.push_back("4mm"); // 99-110
    pgf.push_back("1 m 2 -4 m 2 -4 2");pgn.push_back("-42m"); // 111-122
    pgf.push_back("1 m -2 4 m 2 4 -2 m 2 2 -4 2 -1 -4 2");pgn.push_back("4/mmm"); // 123-142
    // HEXAGONAL-TRIGONAL
    pgf.push_back("1 3 3");pgn.push_back("3"); // 143-146
    pgf.push_back("1 3 3 -1 -3 -3");pgn.push_back("-3"); // 147-148
    pgf.push_back("1 2 2 3 3 2");pgn.push_back("32"); // 149-155
    pgf.push_back("1 3 3 m -2 -2");pgn.push_back("3m"); // 156-161
    pgf.push_back("1 m -2 3 3 -2 -1 2 2 -3 -3 2");pgn.push_back("-3m"); // 162-167

    // HEXAGONAL-HEXAGONAL
    pgf.push_back("1 3 3 2 6 6");pgn.push_back("6"); // 168-173
    pgf.push_back("1 3 3 m -6 -6");pgn.push_back("-6"); // 174
    pgf.push_back("1 3 3 2 6 6 m -6 -6 -1 -3 -3");pgn.push_back("6/m"); // 175-176
    pgf.push_back("1 3 3 2 6 6 2 2 2 2 2 2");pgn.push_back("622"); // 177-182
    pgf.push_back("1 -2 m 3 3 -2 2 -2 m 6 6 -2");pgn.push_back("6mm"); // 183-186
    pgf.push_back("1 m -2 2 m 2 3 -6 3 -6 -2 2");pgn.push_back("-6m2"); // 187-190
    pgf.push_back("1 m -2 2 m 2 3 -6 3 -6 -2 2 2 -1 -2 2 m 2 6 -3 6 -3 -2 2");pgn.push_back("6/mmm"); // 	191-194
    // CUBIC
    pgf.push_back("1 2 3 3 3 3 2 2 3 3 3 3");pgn.push_back("23"); //  195-199
    pgf.push_back("1 3 3 2 3 3 3 3 2 3 3 2 -1 -3 -3 m -3 -3 -3 -3 m -3 -3 m");pgn.push_back("m-3"); // 200-206
    pgf.push_back("1 4 2 4 2 3 4 3 3 4 3 2 2 2 2 2 4 3 2 3 3 4 3 2");pgn.push_back("432"); // 207-214
    pgf.push_back("1 -2 -2 3 3 -2 2 -2 -4 3 3 -4 -4 3 2 -4 3 -4 -4 2 -2 3 3 -2");pgn.push_back("-43m"); // 215-220
    pgf.push_back("1 -2 2 m -2 3 -3 4 3 -2 4 -3 2 -2 -1 2 -3 -3 -2 3 4 4 3 -2 -3 -4 3 2 2 3 -4 -3 2 -4 -3 3 4 m 2 -4 -3 2 3 -4 m 4 -4 2");pgn.push_back("m-3m"); // 221-230
    //  pgf.push_back("1 -2 -2 3 3 -2 -2 2 3 -4 -4 3 3 -4 -2 3 2 -4 -4 3 3 -2 -4 2 -1 2 2 -3 -3 2 2 m -3 4 4 -3 -3 4 2 -3 -2 4 4 -3 -3 2 4 -2");pgn.push_back("m-3m"); // 221-230
    // done
    for(uint i=0;i<pgf.size();i++) aurostd::StringSubst(pgf.at(i),"-2","m"); // -2 are mirrors orthogonal to the axis....

    // -------------------------------------------------------------------- scanning  ORDER=48
    // find point group
    string pgname;
    aurostd::StringSubst(operations,"-2","m");

    for(uint i=0;i<pgf.size()&&!pg_found;i++) {
      if(MapOperations(operations,pgf.at(i))) {pg_found=TRUE;pgname=pgn.at(i);}
    }
    //  cerr << pg_found << endl;
    //   cerr << pgname << endl;

    // -------------------------------------------------------------------- scanning  ORDER=48
    if(pgname=="m-3m") {
      pg_found=TRUE;
      a.crystal_family="cubic";a.crystal_system="cubic";a.point_group_crystal_class="hexoctahedral";
      a.point_group_Shoenflies="O_h";a.point_group_Hermann_Mauguin="m-3m";a.point_group_orbifold="*432";
      a.point_group_type="centrosymmetric";a.point_group_order="48";a.point_group_structure="2 x symmetric";
    }
    // -------------------------------------------------------------------- scanning  ORDER=24
    if(pgname=="-43m") {
      pg_found=TRUE;
      a.crystal_family="cubic";a.crystal_system="cubic";a.point_group_crystal_class="tetrahedral";
      a.point_group_Shoenflies="T_d";a.point_group_Hermann_Mauguin="-43m";a.point_group_orbifold="*332";
      a.point_group_type="none";a.point_group_order="24";a.point_group_structure="symmetric"; //DX20180827 - point group type changed from - to none
    }
    if(pgname=="432") {
      pg_found=TRUE;
      a.crystal_family="cubic";a.crystal_system="cubic";a.point_group_crystal_class="gyroidal";
      a.point_group_Shoenflies="O";a.point_group_Hermann_Mauguin="432";a.point_group_orbifold="432";
      a.point_group_type="enantiomorphic";a.point_group_order="24";a.point_group_structure="symmetric";
    }
    if(pgname=="m-3") {
      pg_found=TRUE;
      a.crystal_family="cubic";a.crystal_system="cubic";a.point_group_crystal_class="diploidal";
      a.point_group_Shoenflies="T_h";a.point_group_Hermann_Mauguin="m-3";a.point_group_orbifold="3*2";
      a.point_group_type="centrosymmetric";a.point_group_order="24";a.point_group_structure="2 x alternating";
    }
    if(pgname=="6/mmm") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="hexagonal";a.point_group_crystal_class="dihexagonal-dipyramidal";
      a.point_group_Shoenflies="D_6h";a.point_group_Hermann_Mauguin="6/mmm";a.point_group_orbifold="*622";
      a.point_group_type="centrosymmetric";a.point_group_order="24";a.point_group_structure="2 x dihedral";
    }
    // -------------------------------------------------------------------- scanning  ORDER=16
    if(pgname=="4/mmm") {
      pg_found=TRUE;
      a.crystal_family="tetragonal";a.crystal_system="tetragonal";a.point_group_crystal_class="ditetragonal-dipyramidal";
      a.point_group_Shoenflies="D_4h";a.point_group_Hermann_Mauguin="4/mmm";a.point_group_orbifold="*422";
      a.point_group_type="centrosymmetric";a.point_group_order="16";a.point_group_structure="2 x dihedral";
    }
    // -------------------------------------------------------------------- scanning  ORDER=12
    if(pgname=="23") {
      pg_found=TRUE;
      a.crystal_family="cubic";a.crystal_system="cubic";a.point_group_crystal_class="tetrahedral";
      a.point_group_Shoenflies="T";a.point_group_Hermann_Mauguin="23";a.point_group_orbifold="332";
      a.point_group_type="enantiomorphic";a.point_group_order="12";a.point_group_structure="alternating";
    }
    if(pgname=="-6m2" || pgname=="-62m") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="hexagonal";a.point_group_crystal_class="ditrigonal-dipyramidal";
      a.point_group_Shoenflies="D_3h";a.point_group_Hermann_Mauguin="-6m2";a.point_group_orbifold="*322";
      a.point_group_type="none";a.point_group_order="12";a.point_group_structure="dihedral"; //DX20180827 - point group type changed from - to none
    }
    if(pgname=="6mm") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="hexagonal";a.point_group_crystal_class="dihexagonal-pyramidal";
      a.point_group_Shoenflies="C_6v";a.point_group_Hermann_Mauguin="6mm";a.point_group_orbifold="*66";
      a.point_group_type="polar";a.point_group_order="12";a.point_group_structure="dihedral";
    }
    if(pgname=="622") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="hexagonal";a.point_group_crystal_class="hexagonal-trapezoidal";
      a.point_group_Shoenflies="D_6";a.point_group_Hermann_Mauguin="622";a.point_group_orbifold="622";
      a.point_group_type="enantiomorphic";a.point_group_order="12";a.point_group_structure="dihedral";
    }
    if(pgname=="6/m") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="hexagonal";a.point_group_crystal_class="hexagonal-dipyramidal";
      a.point_group_Shoenflies="C_6h";a.point_group_Hermann_Mauguin="6/m";a.point_group_orbifold="6*";
      a.point_group_type="centrosymmetric";a.point_group_order="12";a.point_group_structure="2 × cyclic";
    }
    if(pgname=="-3m" || pgname=="-3m1" || pgname=="-31m") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="trigonal";a.point_group_crystal_class="ditrigonal-scalahedral";
      a.point_group_Shoenflies="D_3d";a.point_group_Hermann_Mauguin="-3m";a.point_group_orbifold="2*3";
      a.point_group_type="centrosymmetric";a.point_group_order="12";a.point_group_structure="dihedral";
    }
    // -------------------------------------------------------------------- scanning  ORDER=8
    if(pgname=="-42m" || pgname=="-4m2") {
      pg_found=TRUE;
      a.crystal_family="tetragonal";a.crystal_system="tetragonal";a.point_group_crystal_class="tetragonal-scalenoidal";
      a.point_group_Shoenflies="D_2d";a.point_group_Hermann_Mauguin="-42m";a.point_group_orbifold="2*2";
      a.point_group_type="none";a.point_group_order="8";a.point_group_structure="dihedral"; //DX20180827 - point group type changed from - to none
    }
    if(pgname=="4mm") {
      pg_found=TRUE;
      a.crystal_family="tetragonal";a.crystal_system="tetragonal";a.point_group_crystal_class="ditetragonal-pyramidal";
      a.point_group_Shoenflies="C_4v";a.point_group_Hermann_Mauguin="4mm";a.point_group_orbifold="*44";
      a.point_group_type="polar";a.point_group_order="8";a.point_group_structure="dihedral";
    }
    if(pgname=="422") {
      pg_found=TRUE;
      a.crystal_family="tetragonal";a.crystal_system="tetragonal";a.point_group_crystal_class="tetragonal-trapezoidal";
      a.point_group_Shoenflies="D_4";a.point_group_Hermann_Mauguin="422";a.point_group_orbifold="422";
      a.point_group_type="enantiomorphic";a.point_group_order="8";a.point_group_structure="dihedral";
    }
    if(pgname=="4/m") {
      pg_found=TRUE;
      a.crystal_family="tetragonal";a.crystal_system="tetragonal";a.point_group_crystal_class="tetragonal-dipyramidal";
      a.point_group_Shoenflies="C_4h";a.point_group_Hermann_Mauguin="4/m";a.point_group_orbifold="4*";
      a.point_group_type="centrosymmetric";a.point_group_order="8";a.point_group_structure="2 x cyclic";
    }
    if(pgname=="mmm") { //
      pg_found=TRUE;
      a.crystal_family="orthorhombic";a.crystal_system="orthorhombic";a.point_group_crystal_class="orthorhombic-bipyramidal";
      a.point_group_Shoenflies="D_2h";a.point_group_Hermann_Mauguin="mmm";a.point_group_orbifold="*222";
      a.point_group_type="centrosymmetric";a.point_group_order="8";a.point_group_structure="2 x dihedral";
    }
    // -------------------------------------------------------------------- scanning  ORDER=6
    if(pgname=="-6") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="hexagonal";a.point_group_crystal_class="trigonal-dipyramidal";
      a.point_group_Shoenflies="C_3h";a.point_group_Hermann_Mauguin="-6";a.point_group_orbifold="3*";
      a.point_group_type="none";a.point_group_order="6";a.point_group_structure="cyclic"; //DX20180827 - point group type changed from - to none
    }
    if(pgname=="6") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="hexagonal";a.point_group_crystal_class="hexagonal-pyramidal";
      a.point_group_Shoenflies="C_6";a.point_group_Hermann_Mauguin="6";a.point_group_orbifold="66";
      a.point_group_type="enantiomorphic polar";a.point_group_order="6";a.point_group_structure="cyclic";
    }
    if(pgname=="3m" || pgname=="3m1" || pgname=="31m") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="trigonal";a.point_group_crystal_class="ditrigonal-pyramidal";
      a.point_group_Shoenflies="C_3v";a.point_group_Hermann_Mauguin="3m";a.point_group_orbifold="*33";
      a.point_group_type="polar";a.point_group_order="6";a.point_group_structure="dihedral";
    }
    if(pgname=="32" || pgname=="321" || pgname=="312") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="trigonal";a.point_group_crystal_class="trigonal-trapezoidal";
      a.point_group_Shoenflies="D_3";a.point_group_Hermann_Mauguin="32";a.point_group_orbifold="322";
      a.point_group_type="enantiomorphic";a.point_group_order="6";a.point_group_structure="dihedral";
    }
    if(pgname=="-3") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="trigonal";a.point_group_crystal_class="rhombohedral";
      a.point_group_Shoenflies="S_6";a.point_group_Hermann_Mauguin="-3";a.point_group_orbifold="3x";
      a.point_group_type="centrosymmetric";a.point_group_order="6";a.point_group_structure="cyclic";
    }
    // -------------------------------------------------------------------- scanning  ORDER=4
    if(pgname=="-4") {
      pg_found=TRUE;
      a.crystal_family="tetragonal";a.crystal_system="tetragonal";a.point_group_crystal_class="tetragonal-disphenoidal";
      a.point_group_Shoenflies="S_4";a.point_group_Hermann_Mauguin="-4";a.point_group_orbifold="2x";
      a.point_group_type="none";a.point_group_order="4";a.point_group_structure="cyclic"; //DX20180827 - point group type changed from - to none
    }
    if(pgname=="4") {
      pg_found=TRUE;
      a.crystal_family="tetragonal";a.crystal_system="tetragonal";a.point_group_crystal_class="tetragonal-pyramidal";
      a.point_group_Shoenflies="C_4";a.point_group_Hermann_Mauguin="4";a.point_group_orbifold="44";
      a.point_group_type="enantiomorphic polar";a.point_group_order="4";a.point_group_structure="2 x cyclic";
    }
    if(pgname=="mm2") {
      pg_found=TRUE;
      a.crystal_family="orthorhombic";a.crystal_system="orthorhombic";a.point_group_crystal_class="orthorhombic-pyramidal";
      a.point_group_Shoenflies="C_2v";a.point_group_Hermann_Mauguin="mm2";a.point_group_orbifold="*22";
      a.point_group_type="polar";a.point_group_order="4";a.point_group_structure="dihedral";
    }
    if(pgname=="222") {
      pg_found=TRUE;
      a.crystal_family="orthorhombic";a.crystal_system="orthorhombic";a.point_group_crystal_class="orthorhombic-sphenoidal";
      a.point_group_Shoenflies="D_2";a.point_group_Hermann_Mauguin="222";a.point_group_orbifold="222";
      a.point_group_type="enantiomorphic";a.point_group_order="4";a.point_group_structure="dihedral";
    }
    if(pgname=="2/m") {
      pg_found=TRUE;
      a.crystal_family="monoclinic";a.crystal_system="monoclinic";a.point_group_crystal_class="monoclinic-prismatic";
      a.point_group_Shoenflies="C_2h";a.point_group_Hermann_Mauguin="2/m";a.point_group_orbifold="2*";
      a.point_group_type="centrosymmetric";a.point_group_order="4";a.point_group_structure="2 x cyclic";
    }
    // -------------------------------------------------------------------- scanning  ORDER=3
    if(pgname=="3") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="trigonal";a.point_group_crystal_class="trigonal-pyramidal";
      a.point_group_Shoenflies="C_3";a.point_group_Hermann_Mauguin="3";a.point_group_orbifold="33";
      a.point_group_type="enantiomorphic polar";a.point_group_order="3";a.point_group_structure="cyclic";
    }
    // -------------------------------------------------------------------- scanning  ORDER=2
    if(pgname=="m") {
      pg_found=TRUE;
      a.crystal_family="monoclinic";a.crystal_system="monoclinic";a.point_group_crystal_class="monoclinic-domatic";
      a.point_group_Shoenflies="C_S";a.point_group_Hermann_Mauguin="m";a.point_group_orbifold="1*";
      a.point_group_type="polar";a.point_group_order="2";a.point_group_structure="cyclic";
    }
    if(pgname=="2") {
      pg_found=TRUE;
      a.crystal_family="monoclinic";a.crystal_system="monoclinic";a.point_group_crystal_class="monoclinic-sphenoidal";
      a.point_group_Shoenflies="C_2";a.point_group_Hermann_Mauguin="2";a.point_group_orbifold="22";
      a.point_group_type="enantiomorphic polar";a.point_group_order="2";a.point_group_structure="cyclic";
    }
    if(pgname=="-1") {
      pg_found=TRUE;
      a.crystal_family="triclinic";a.crystal_system="triclinic";a.point_group_crystal_class="triclinic-pinacoidal";
      a.point_group_Shoenflies="C_i";a.point_group_Hermann_Mauguin="-1";a.point_group_orbifold="x";
      a.point_group_type="centrosymmetric";a.point_group_order="2";a.point_group_structure="cyclic";
    }
    // -------------------------------------------------------------------- scanning  ORDER=1
    if(pgname=="1") {
      pg_found=TRUE;
      a.crystal_family="triclinic";a.crystal_system="triclinic";a.point_group_crystal_class="triclinic-pedial";
      a.point_group_Shoenflies="C_1";a.point_group_Hermann_Mauguin="1";a.point_group_orbifold="11";
      a.point_group_type="enantiomorphic polar";a.point_group_order="1";a.point_group_structure="trivial";
    }


    // ------------------------------------------------------------------------------------------- PRINTING AND LEAVING
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message
      << " " << SEPARATION_LINE_DASH_SHORT << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message
      << " Symmetry: unique point group operations " << a.pgroup_xtal.size() << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: " << operations << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Crystal Family = " << a.crystal_family << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Crystal System = " << a.crystal_system << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Crystal Class = " << a.point_group_crystal_class << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Schoenflies notation = " << a.point_group_Shoenflies << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Hermann Mauguin notation = " << a.point_group_Hermann_Mauguin << endl; //"   - " << a.title << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Orbifold = " << a.point_group_orbifold << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Point Group Type = " << a.point_group_type << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Point Group Order = " << a.point_group_order << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Point Group Structure = " << a.point_group_structure << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " " << SEPARATION_LINE_DASH_SHORT << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: END " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message
      << " " << SEPARATION_LINE_DASH_SHORT << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    if(_write_) Krun=Krun && KBIN_SymmetryWrite(FileMESSAGE,a,aflags,_PGROUP_XTAL_,osswrite,oss,format);
    return Krun;
  }
} // namespace SYM
//DX+CO END

#ifndef COMPILE_SLIM
namespace SYM {
  bool CalculatePointGroupCrystal_20160101(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_) {
    ostringstream aus;
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    DEBUG_SYMMETRY=DEBUG_SYMMETRY || LDEBUG;    
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUP_CRYSTAL Symmetry: BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);

    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal [0]" << endl;
    // this function returns the PointGroup of the Crystal by searching through the Hermann-Mauguin combinations

    bool Krun=TRUE;
    xvector<double> ctau(3),ftau(3),ctrasl(3),ftrasl(3);  // translation
    xmatrix<double> Uf(3,3),Uc(3,3);                      // matrices
    std::vector<int> basis_atoms_map(a.atoms.size());           // will map each on each
    for(uint i=0;i<a.atoms.size();i++) basis_atoms_map[i]=i;    // identically map each over each
    std::vector<int> basis_types_map(a.atoms.size());           // will map each on each
    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal [0d]" << endl;
    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal a.atoms.size()=" << a.atoms.size() << endl;
    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal a.species.size()=" << a.species.size() << endl;
    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal basis_types_map.size()=" << basis_types_map.size() << endl;
    for(uint i=0;i<a.atoms.size();i++) basis_types_map[i]=a.atoms[i].type;    // identically map each over each
    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal [0e]" << endl;
    string message="PGROUP_XTAL";
    //  double _eps_=_EPS_;

    xstructure b(a),t(a);   // for backup... will use at the end to plug back the postions
    // a.GetPrimitive(); nooooo I should not change the structure as it willl bust the factor group later
    a.BringInCell();
    // a.BringInWignerSeitz(); //I could but I would loose some.

    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal [1]" << endl;

    // aflags.QUIET=TRUE;
    if(a.pgroup_calculated==FALSE) Krun=Krun && SYM::CalculatePointGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss);        // NEED POINT GROUP
    if(a.fgroup_calculated==FALSE) Krun=Krun && SYM::CalculateFactorGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss);       // NEED FACTOR GROUP
    // if(a.sgroup_calculated==FALSE) Krun=Krun && SYM::CalculateSpaceGroup(FileMESSAGE,a,aflags,FALSE,osswrite,oss);          // NEED SPACE GROUP
    // if(a.iatoms_calculated==FALSE) Krun=Krun && SYM::CalculateInequivalentAtoms(FileMESSAGE,a,aflags,_write_,osswrite,oss); // NEED EQUIV ATOMS
    // if(a.iatoms_calculated==FALSE) Krun=Krun && SYM::CalculateSitePointGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss);    // NEED EQUIV ATOMS

    // this one removes the double matrices so it takes care of non primitive systems
    for(uint ip=0;ip<a.fgroup.size();ip++) {
      clear(ctau);clear(ftau);clear(ctrasl);clear(ftrasl);
      // if(modulus(a.fgroup[kk].ftau)<0.1)
      Uf=a.fgroup[ip].Uf;Uc=a.fgroup[ip].Uc;ftau=a.fgroup[ip].ftau;ctau=a.fgroup[ip].ctau;
      bool sym_found=FALSE;
      for(uint ii=0;ii<a.pgroup_xtal.size()&&!sym_found;ii++)
        //   sym_found=(identical(Uf,a.pgroup_xtal[ii].Uf,_eps_) && identical(ftau,a.pgroup_xtal[ii].ftau,_eps_));
        sym_found=(identical(Uf,a.pgroup_xtal[ii].Uf,_eps_));     // look in all the list of operations
      if(sym_found==FALSE) {                                    // new operation, generate and save it
        // [UNUSED] uint kk; // warning: variable ‘kk’ set but not used 
        // [UNUSED] kk=
        //DX+CO START
        SYM::AddSymmetryToStructure(a,Uc,Uf,ctau,ftau,ctrasl,ftrasl,basis_atoms_map,basis_types_map,true,_PGROUP_XTAL_);   // in kk the number of pgroup_xtal
        //DX+CO END
        // aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUP_XTL " << a.pgroup_xtal[kk-1].str_type
        // << " theta="; if(a.pgroup_xtal[kk-1].angle<100) aus << " "; if(a.pgroup_xtal[kk-1].angle<10) aus << " ";
        // aus << a.pgroup_xtal[kk-1].angle << " " << " r=(" << a.pgroup_xtal[kk-1].axis << ")"
        // << endl;  // remember vectors start from 0
        // aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
      }
    }
    a.pgroup_xtal_calculated=TRUE;

    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal [2]" << endl;

    // PGROUP_XTAL
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " " << SEPARATION_LINE_DASH_SHORT << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    for(uint kk=1;kk<=a.pgroup_xtal.size();kk++) {  // shift by 1 to use [kk-1]...
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " " ;
      aus << a.pgroup_xtal[kk-1].str_type << " theta=";
      if(a.pgroup_xtal[kk-1].angle<100) aus << " ";
      if(a.pgroup_xtal[kk-1].angle<10)  aus << " ";
      aus << a.pgroup_xtal[kk-1].angle << " " << " r=(" << a.pgroup_xtal[kk-1].axis << ")";
      aus << "    HM=" <<  aurostd::PaddedPRE(a.pgroup_xtal[kk-1].str_Hermann_Mauguin,2," ");     // remember vectors start from 0
      aus << "    S=" <<  aurostd::PaddedPRE(a.pgroup_xtal[kk-1].str_Schoenflies,2," ");          // remember vectors start from 0
      aus    << endl;  // remember vectors start from 0
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    }

    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculatePointGroupCrystal [3]" << endl;

    // LOOK up tabe for the point group operations of the crystal
    vector<string> vops;vops.clear();
    string operations="";
    for(uint kk=1;kk<=a.pgroup_xtal.size();kk++) {  // shift by 1 to use [kk-1]...
      vops.push_back(a.pgroup_xtal[kk-1].str_Hermann_Mauguin);
      operations=operations+a.pgroup_xtal[kk-1].str_Hermann_Mauguin;
      if(kk<a.pgroup_xtal.size()) operations=operations+" ";
    }
    a.crystal_family="";a.crystal_system="";a.point_group_crystal_class="";
    a.point_group_Shoenflies="";a.point_group_Hermann_Mauguin="";a.point_group_orbifold="";
    a.point_group_type="";a.point_group_order="";a.point_group_structure="";
    bool pg_found=FALSE;

    vector<string> pgf,pgn;
    // TRICLINIC
    pgf.push_back("1");pgn.push_back("1"); // 1
    pgf.push_back("1 -1");pgn.push_back("-1"); // 2
    // MONOCLINIC
    pgf.push_back("1 2");pgn.push_back("2"); // 3-5
    pgf.push_back("1 m");pgn.push_back("m"); // 6-9
    pgf.push_back("1 2 m -1");pgn.push_back("2/m"); // 10-15
    // ORTHORHOMBIC
    pgf.push_back("1 2 2 2");pgn.push_back("222"); // 16-24
    pgf.push_back("1 m m 2");pgn.push_back("mm2"); // 25-46
    pgf.push_back("1 m m 2 m 2 2 -1");pgn.push_back("mmm"); // 47-74
    // TETRAGONAL
    pgf.push_back("1 4 2 4");pgn.push_back("4"); // 75-80
    pgf.push_back("1 2 -4 -4");pgn.push_back("-4"); // 81-82
    pgf.push_back("1 4 2 4 m -4 -1 -4");pgn.push_back("4/m"); // 83-88
    pgf.push_back("1 4 2 4 2 2 2 2");pgn.push_back("422"); // 89-98
    pgf.push_back("1 m -2 4 m 2 4 -2");pgn.push_back("4mm"); // 99-110
    pgf.push_back("1 m 2 -4 m 2 -4 2");pgn.push_back("-42m"); // 111-122
    pgf.push_back("1 m -2 4 m 2 4 -2 m 2 2 -4 2 -1 -4 2");pgn.push_back("4/mmm"); // 123-142
    // HEXAGONAL-TRIGONAL
    pgf.push_back("1 3 3");pgn.push_back("3"); // 143-146
    pgf.push_back("1 3 3 -1 -3 -3");pgn.push_back("-3"); // 147-148
    pgf.push_back("1 2 2 3 3 2");pgn.push_back("32"); // 149-155
    pgf.push_back("1 3 3 m -2 -2");pgn.push_back("3m"); // 156-161
    pgf.push_back("1 m -2 3 3 -2 -1 2 2 -3 -3 2");pgn.push_back("-3m"); // 162-167

    // HEXAGONAL-HEXAGONAL
    pgf.push_back("1 3 3 2 6 6");pgn.push_back("6"); // 168-173
    pgf.push_back("1 3 3 m -6 -6");pgn.push_back("-6"); // 174
    pgf.push_back("1 3 3 2 6 6 m -6 -6 -1 -3 -3");pgn.push_back("6/m"); // 175-176
    pgf.push_back("1 3 3 2 6 6 2 2 2 2 2 2");pgn.push_back("622"); // 177-182
    pgf.push_back("1 -2 m 3 3 -2 2 -2 m 6 6 -2");pgn.push_back("6mm"); // 183-186
    pgf.push_back("1 m -2 2 m 2 3 -6 3 -6 -2 2");pgn.push_back("-6m2"); // 187-190
    pgf.push_back("1 m -2 2 m 2 3 -6 3 -6 -2 2 2 -1 -2 2 m 2 6 -3 6 -3 -2 2");pgn.push_back("6/mmm"); // 	191-194
    // CUBIC
    pgf.push_back("1 2 3 3 3 3 2 2 3 3 3 3");pgn.push_back("23"); //  195-199
    pgf.push_back("1 3 3 2 3 3 3 3 2 3 3 2 -1 -3 -3 m -3 -3 -3 -3 m -3 -3 m");pgn.push_back("m-3"); // 200-206
    pgf.push_back("1 4 2 4 2 3 4 3 3 4 3 2 2 2 2 2 4 3 2 3 3 4 3 2");pgn.push_back("432"); // 207-214
    pgf.push_back("1 -2 -2 3 3 -2 2 -2 -4 3 3 -4 -4 3 2 -4 3 -4 -4 2 -2 3 3 -2");pgn.push_back("-43m"); // 215-220
    pgf.push_back("1 -2 2 m -2 3 -3 4 3 -2 4 -3 2 -2 -1 2 -3 -3 -2 3 4 4 3 -2 -3 -4 3 2 2 3 -4 -3 2 -4 -3 3 4 m 2 -4 -3 2 3 -4 m 4 -4 2");pgn.push_back("m-3m"); // 221-230
    //  pgf.push_back("1 -2 -2 3 3 -2 -2 2 3 -4 -4 3 3 -4 -2 3 2 -4 -4 3 3 -2 -4 2 -1 2 2 -3 -3 2 2 m -3 4 4 -3 -3 4 2 -3 -2 4 4 -3 -3 2 4 -2");pgn.push_back("m-3m"); // 221-230
    // done
    for(uint i=0;i<pgf.size();i++) aurostd::StringSubst(pgf.at(i),"-2","m"); // -2 are mirrors orthogonal to the axis....

    // -------------------------------------------------------------------- scanning  ORDER=48
    // find point group
    string pgname;
    aurostd::StringSubst(operations,"-2","m");

    for(uint i=0;i<pgf.size()&&!pg_found;i++) {
      if(MapOperations(operations,pgf.at(i))) {pg_found=TRUE;pgname=pgn.at(i);}
    }
    //  cerr << pg_found << endl;
    //   cerr << pgname << endl;

    // -------------------------------------------------------------------- scanning  ORDER=48
    if(pgname=="m-3m") {
      pg_found=TRUE;
      a.crystal_family="cubic";a.crystal_system="cubic";a.point_group_crystal_class="hexoctahedral";
      a.point_group_Shoenflies="O_h";a.point_group_Hermann_Mauguin="m-3m";a.point_group_orbifold="*432";
      a.point_group_type="centrosymmetric";a.point_group_order="48";a.point_group_structure="2 x symmetric";
    }
    // -------------------------------------------------------------------- scanning  ORDER=24
    if(pgname=="-43m") {
      pg_found=TRUE;
      a.crystal_family="cubic";a.crystal_system="cubic";a.point_group_crystal_class="tetrahedral";
      a.point_group_Shoenflies="T_d";a.point_group_Hermann_Mauguin="-43m";a.point_group_orbifold="*332";
      a.point_group_type="none";a.point_group_order="24";a.point_group_structure="symmetric"; //DX20180827 - point group type changed from - to none
    }
    if(pgname=="432") {
      pg_found=TRUE;
      a.crystal_family="cubic";a.crystal_system="cubic";a.point_group_crystal_class="gyroidal";
      a.point_group_Shoenflies="O";a.point_group_Hermann_Mauguin="432";a.point_group_orbifold="432";
      a.point_group_type="enantiomorphic";a.point_group_order="24";a.point_group_structure="symmetric";
    }
    if(pgname=="m-3") {
      pg_found=TRUE;
      a.crystal_family="cubic";a.crystal_system="cubic";a.point_group_crystal_class="diploidal";
      a.point_group_Shoenflies="T_h";a.point_group_Hermann_Mauguin="m-3";a.point_group_orbifold="3*2";
      a.point_group_type="centrosymmetric";a.point_group_order="24";a.point_group_structure="2 x alternating";
    }
    if(pgname=="6/mmm") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="hexagonal";a.point_group_crystal_class="dihexagonal-dipyramidal";
      a.point_group_Shoenflies="D_6h";a.point_group_Hermann_Mauguin="6/mmm";a.point_group_orbifold="*622";
      a.point_group_type="centrosymmetric";a.point_group_order="24";a.point_group_structure="2 x dihedral";
    }
    // -------------------------------------------------------------------- scanning  ORDER=16
    if(pgname=="4/mmm") {
      pg_found=TRUE;
      a.crystal_family="tetragonal";a.crystal_system="tetragonal";a.point_group_crystal_class="ditetragonal-dipyramidal";
      a.point_group_Shoenflies="D_4h";a.point_group_Hermann_Mauguin="4/mmm";a.point_group_orbifold="*422";
      a.point_group_type="centrosymmetric";a.point_group_order="16";a.point_group_structure="2 x dihedral";
    }
    // -------------------------------------------------------------------- scanning  ORDER=12
    if(pgname=="23") {
      pg_found=TRUE;
      a.crystal_family="cubic";a.crystal_system="cubic";a.point_group_crystal_class="tetrahedral";
      a.point_group_Shoenflies="T";a.point_group_Hermann_Mauguin="23";a.point_group_orbifold="332";
      a.point_group_type="enantiomorphic";a.point_group_order="12";a.point_group_structure="alternating";
    }
    if(pgname=="-6m2" || pgname=="-62m") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="hexagonal";a.point_group_crystal_class="ditrigonal-dipyramidal";
      a.point_group_Shoenflies="D_3h";a.point_group_Hermann_Mauguin="-6m2";a.point_group_orbifold="*322";
      a.point_group_type="none";a.point_group_order="12";a.point_group_structure="dihedral"; //DX20180827 - point group type changed from - to none
    }
    if(pgname=="6mm") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="hexagonal";a.point_group_crystal_class="dihexagonal-pyramidal";
      a.point_group_Shoenflies="C_6v";a.point_group_Hermann_Mauguin="6mm";a.point_group_orbifold="*66";
      a.point_group_type="polar";a.point_group_order="12";a.point_group_structure="dihedral";
    }
    if(pgname=="622") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="hexagonal";a.point_group_crystal_class="hexagonal-trapezoidal";
      a.point_group_Shoenflies="D_6";a.point_group_Hermann_Mauguin="622";a.point_group_orbifold="622";
      a.point_group_type="enantiomorphic";a.point_group_order="12";a.point_group_structure="dihedral";
    }
    if(pgname=="6/m") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="hexagonal";a.point_group_crystal_class="hexagonal-dipyramidal";
      a.point_group_Shoenflies="C_6h";a.point_group_Hermann_Mauguin="6/m";a.point_group_orbifold="6*";
      a.point_group_type="centrosymmetric";a.point_group_order="12";a.point_group_structure="2 × cyclic";
    }
    if(pgname=="-3m" || pgname=="-3m1" || pgname=="-31m") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="trigonal";a.point_group_crystal_class="ditrigonal-scalahedral";
      a.point_group_Shoenflies="D_3d";a.point_group_Hermann_Mauguin="-3m";a.point_group_orbifold="2*3";
      a.point_group_type="centrosymmetric";a.point_group_order="12";a.point_group_structure="dihedral";
    }
    // -------------------------------------------------------------------- scanning  ORDER=8
    if(pgname=="-42m" || pgname=="-4m2") {
      pg_found=TRUE;
      a.crystal_family="tetragonal";a.crystal_system="tetragonal";a.point_group_crystal_class="tetragonal-scalenoidal";
      a.point_group_Shoenflies="D_2d";a.point_group_Hermann_Mauguin="-42m";a.point_group_orbifold="2*2";
      a.point_group_type="none";a.point_group_order="8";a.point_group_structure="dihedral"; //DX20180827 - point group type changed from - to none
    }
    if(pgname=="4mm") {
      pg_found=TRUE;
      a.crystal_family="tetragonal";a.crystal_system="tetragonal";a.point_group_crystal_class="ditetragonal-pyramidal";
      a.point_group_Shoenflies="C_4v";a.point_group_Hermann_Mauguin="4mm";a.point_group_orbifold="*44";
      a.point_group_type="polar";a.point_group_order="8";a.point_group_structure="dihedral";
    }
    if(pgname=="422") {
      pg_found=TRUE;
      a.crystal_family="tetragonal";a.crystal_system="tetragonal";a.point_group_crystal_class="tetragonal-trapezoidal";
      a.point_group_Shoenflies="D_4";a.point_group_Hermann_Mauguin="422";a.point_group_orbifold="422";
      a.point_group_type="enantiomorphic";a.point_group_order="8";a.point_group_structure="dihedral";
    }
    if(pgname=="4/m") {
      pg_found=TRUE;
      a.crystal_family="tetragonal";a.crystal_system="tetragonal";a.point_group_crystal_class="tetragonal-dipyramidal";
      a.point_group_Shoenflies="C_4h";a.point_group_Hermann_Mauguin="4/m";a.point_group_orbifold="4*";
      a.point_group_type="centrosymmetric";a.point_group_order="8";a.point_group_structure="2 x cyclic";
    }
    if(pgname=="mmm") { //
      pg_found=TRUE;
      a.crystal_family="orthorhombic";a.crystal_system="orthorhombic";a.point_group_crystal_class="orthorhombic-bipyramidal";
      a.point_group_Shoenflies="D_2h";a.point_group_Hermann_Mauguin="mmm";a.point_group_orbifold="*222";
      a.point_group_type="centrosymmetric";a.point_group_order="8";a.point_group_structure="2 x dihedral";
    }
    // -------------------------------------------------------------------- scanning  ORDER=6
    if(pgname=="-6") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="hexagonal";a.point_group_crystal_class="trigonal-dipyramidal";
      a.point_group_Shoenflies="C_3h";a.point_group_Hermann_Mauguin="-6";a.point_group_orbifold="3*";
      a.point_group_type="none";a.point_group_order="6";a.point_group_structure="cyclic"; //DX20180827 - point group type changed from - to none
    }
    if(pgname=="6") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="hexagonal";a.point_group_crystal_class="hexagonal-pyramidal";
      a.point_group_Shoenflies="C_6";a.point_group_Hermann_Mauguin="6";a.point_group_orbifold="66";
      a.point_group_type="enantiomorphic polar";a.point_group_order="6";a.point_group_structure="cyclic";
    }
    if(pgname=="3m" || pgname=="3m1" || pgname=="31m") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="trigonal";a.point_group_crystal_class="ditrigonal-pyramidal";
      a.point_group_Shoenflies="C_3v";a.point_group_Hermann_Mauguin="3m";a.point_group_orbifold="*33";
      a.point_group_type="polar";a.point_group_order="6";a.point_group_structure="dihedral";
    }
    if(pgname=="32" || pgname=="321" || pgname=="312") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="trigonal";a.point_group_crystal_class="trigonal-trapezoidal";
      a.point_group_Shoenflies="D_3";a.point_group_Hermann_Mauguin="32";a.point_group_orbifold="322";
      a.point_group_type="enantiomorphic";a.point_group_order="6";a.point_group_structure="dihedral";
    }
    if(pgname=="-3") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="trigonal";a.point_group_crystal_class="rhombohedral";
      a.point_group_Shoenflies="S_6";a.point_group_Hermann_Mauguin="-3";a.point_group_orbifold="3x";
      a.point_group_type="centrosymmetric";a.point_group_order="6";a.point_group_structure="cyclic";
    }
    // -------------------------------------------------------------------- scanning  ORDER=4
    if(pgname=="-4") {
      pg_found=TRUE;
      a.crystal_family="tetragonal";a.crystal_system="tetragonal";a.point_group_crystal_class="tetragonal-disphenoidal";
      a.point_group_Shoenflies="S_4";a.point_group_Hermann_Mauguin="-4";a.point_group_orbifold="2x";
      a.point_group_type="none";a.point_group_order="4";a.point_group_structure="cyclic"; //DX20180827 - point group type changed from - to none
    }
    if(pgname=="4") {
      pg_found=TRUE;
      a.crystal_family="tetragonal";a.crystal_system="tetragonal";a.point_group_crystal_class="tetragonal-pyramidal";
      a.point_group_Shoenflies="C_4";a.point_group_Hermann_Mauguin="4";a.point_group_orbifold="44";
      a.point_group_type="enantiomorphic polar";a.point_group_order="4";a.point_group_structure="2 x cyclic";
    }
    if(pgname=="mm2") {
      pg_found=TRUE;
      a.crystal_family="orthorhombic";a.crystal_system="orthorhombic";a.point_group_crystal_class="orthorhombic-pyramidal";
      a.point_group_Shoenflies="C_2v";a.point_group_Hermann_Mauguin="mm2";a.point_group_orbifold="*22";
      a.point_group_type="polar";a.point_group_order="4";a.point_group_structure="dihedral";
    }
    if(pgname=="222") {
      pg_found=TRUE;
      a.crystal_family="orthorhombic";a.crystal_system="orthorhombic";a.point_group_crystal_class="orthorhombic-sphenoidal";
      a.point_group_Shoenflies="D_2";a.point_group_Hermann_Mauguin="222";a.point_group_orbifold="222";
      a.point_group_type="enantiomorphic";a.point_group_order="4";a.point_group_structure="dihedral";
    }
    if(pgname=="2/m") {
      pg_found=TRUE;
      a.crystal_family="monoclinic";a.crystal_system="monoclinic";a.point_group_crystal_class="monoclinic-prismatic";
      a.point_group_Shoenflies="C_2h";a.point_group_Hermann_Mauguin="2/m";a.point_group_orbifold="2*";
      a.point_group_type="centrosymmetric";a.point_group_order="4";a.point_group_structure="2 x cyclic";
    }
    // -------------------------------------------------------------------- scanning  ORDER=3
    if(pgname=="3") {
      pg_found=TRUE;
      a.crystal_family="hexagonal";a.crystal_system="trigonal";a.point_group_crystal_class="trigonal-pyramidal";
      a.point_group_Shoenflies="C_3";a.point_group_Hermann_Mauguin="3";a.point_group_orbifold="33";
      a.point_group_type="enantiomorphic polar";a.point_group_order="3";a.point_group_structure="cyclic";
    }
    // -------------------------------------------------------------------- scanning  ORDER=2
    if(pgname=="m") {
      pg_found=TRUE;
      a.crystal_family="monoclinic";a.crystal_system="monoclinic";a.point_group_crystal_class="monoclinic-domatic";
      a.point_group_Shoenflies="C_S";a.point_group_Hermann_Mauguin="m";a.point_group_orbifold="1*";
      a.point_group_type="polar";a.point_group_order="2";a.point_group_structure="cyclic";
    }
    if(pgname=="2") {
      pg_found=TRUE;
      a.crystal_family="monoclinic";a.crystal_system="monoclinic";a.point_group_crystal_class="monoclinic-sphenoidal";
      a.point_group_Shoenflies="C_2";a.point_group_Hermann_Mauguin="2";a.point_group_orbifold="22";
      a.point_group_type="enantiomorphic polar";a.point_group_order="2";a.point_group_structure="cyclic";
    }
    if(pgname=="-1") {
      pg_found=TRUE;
      a.crystal_family="triclinic";a.crystal_system="triclinic";a.point_group_crystal_class="triclinic-pinacoidal";
      a.point_group_Shoenflies="C_i";a.point_group_Hermann_Mauguin="-1";a.point_group_orbifold="x";
      a.point_group_type="centrosymmetric";a.point_group_order="2";a.point_group_structure="cyclic";
    }
    // -------------------------------------------------------------------- scanning  ORDER=1
    if(pgname=="1") {
      pg_found=TRUE;
      a.crystal_family="triclinic";a.crystal_system="triclinic";a.point_group_crystal_class="triclinic-pedial";
      a.point_group_Shoenflies="C_1";a.point_group_Hermann_Mauguin="1";a.point_group_orbifold="11";
      a.point_group_type="enantiomorphic polar";a.point_group_order="1";a.point_group_structure="trivial";
    }


    // ------------------------------------------------------------------------------------------- PRINTING AND LEAVING
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message
      << " " << SEPARATION_LINE_DASH_SHORT << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message
      << " Symmetry: unique point group operations " << a.pgroup_xtal.size() << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: " << operations << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Crystal Family = " << a.crystal_family << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Crystal System = " << a.crystal_system << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Crystal Class = " << a.point_group_crystal_class << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Schoenflies notation = " << a.point_group_Shoenflies << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Hermann Mauguin notation = " << a.point_group_Hermann_Mauguin << endl; //"   - " << a.title << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Orbifold = " << a.point_group_orbifold << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Point Group Type = " << a.point_group_type << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Point Group Order = " << a.point_group_order << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: Point Group Structure = " << a.point_group_structure << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " " << SEPARATION_LINE_DASH_SHORT << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: END " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message
      << " " << SEPARATION_LINE_DASH_SHORT << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    if(_write_) Krun=Krun && KBIN_SymmetryWrite(FileMESSAGE,a,aflags,_PGROUP_XTAL_,osswrite,oss);
    return Krun;
  }
} // namespace SYM
//DX+CO END
#endif

// ----------------------------------------------------------------------------------------------------------------------------------------------------
// Function SYM::AddSymmetryToStructure
//
// this function adds the operation to the list of
// pgroup operations in the xstructure structure -  SC aug2007
namespace SYM {
  uint AddSymmetryToStructure(xstructure &a,const uint& iat,
      const xmatrix<double> &Uc,const xmatrix<double> &Uf,const xvector<double> &ctau,const xvector<double> &ftau,
      const xvector<double> &ctrasl,const xvector<double> &ftrasl,
      const std::vector<int> &basis_atoms_map,
      const std::vector<int> &basis_types_map,
      bool basis_map_calculated,
      char group) {
    //DX+CO START
    bool roff=TRUE;
    return AddSymmetryToStructure(a,iat,Uc,Uf,ctau,ftau,ctrasl,ftrasl,basis_atoms_map,basis_types_map,basis_map_calculated,group,roff);
  }
  uint AddSymmetryToStructure(xstructure &a,const uint& iat,
      const xmatrix<double> &Uc,const xmatrix<double> &Uf,const xvector<double> &ctau,const xvector<double> &ftau,
      const xvector<double> &ctrasl,const xvector<double> &ftrasl,
      const std::vector<int> &basis_atoms_map,
      const std::vector<int> &basis_types_map,
      bool basis_map_calculated,
      char group,
      bool roff) {
    //DX+CO END
    xmatrix<double> _Uc(3,3);_Uc=Uc;             // here we make copies so we do not mess up input
    xmatrix<double> _Uf(3,3);_Uf=Uf;             // here we make copies so we do not mess up input
    string _string,_stringHM,_stringSC; bool _inversion;
    xvector<double> _axis(3);
    xvector<double> _generator_coefficients(3); //DX20171206 - added generator coefficients
    xmatrix<double> _generator(3,3);
    xmatrix<xcomplex<double> > _SU2_matrix(2,2); //DX20180117 - added SU(2) matrix
    xvector<xcomplex<double> > _su2_coefficients(3); //DX20180117 - added su(2) coefficients
    double _angle;
    _sym_op symop;
    symop.basis_atoms_map.clear(); // empty
    for(uint i=0;i<basis_atoms_map.size();i++) symop.basis_atoms_map.push_back(basis_atoms_map.at(i));
    symop.basis_types_map.clear(); // empty
    for(uint i=0;i<basis_types_map.size();i++) symop.basis_types_map.push_back(basis_types_map.at(i));
    //DX+CO START
    symop.basis_map_calculated=basis_map_calculated;
    //DX+CO END
    if(group==_PGROUP_ || group==_PGROUP_XTAL_ || group==_FGROUP_ || group==_SGROUP_ || group==_AGROUP_ || group==_PGROUPK_ || group==_PGROUPK_XTAL_ || group==_PGROUPK_PATTERSON_) { //DX20171205 - Added pgroupk_xtal //DX20200129 - added Patterson symmetry
      //DX+CO START
      if(roff) {roundoff(_Uc,_EPS_roundoff_);roundoff(_Uf,_EPS_roundoff_);}                       // Uc cleanup from roundoff errors
      if(roff) {roundoff(ctau,_EPS_roundoff_);roundoff(ftau,_EPS_roundoff_);}
      if(roff) {roundoff(ctrasl,_EPS_roundoff_);roundoff(ftrasl,_EPS_roundoff_);}
      //DX+CO END
      symop.Uc=_Uc;symop.Uf=_Uf;symop.ctau=ctau;symop.ftau=ftau;symop.ctrasl=ctrasl;symop.ftrasl=ftrasl;
      SYM::TypePointGroupOperation(_Uc,_Uf,_string,_inversion,_angle,_axis,_generator,_generator_coefficients,_SU2_matrix,_su2_coefficients,a.sym_eps); // extract information in //DX added eps //DX20171206 - generator coefficients //DX20171207 - added Uf //DX20180117 - add SU(2) and su(2) coefficients
      SYM::TypePointGroupOperationInternational(_Uc,_stringHM,_stringSC,_inversion,_angle,_axis,_generator,_generator_coefficients,_SU2_matrix,_su2_coefficients,a.sym_eps); // extract information //DX added eps //DX20171206 - generator coefficients //DX20180117 - add SU(2) and su(2) coefficients
      symop.generator=_generator;symop.generator_coefficients=_generator_coefficients;symop.angle=_angle;symop.axis=_axis; //DX20171206 - generator coefficients
      symop.SU2_matrix=_SU2_matrix;symop.su2_coefficients=_su2_coefficients; //DX20180117 - add SU(2) and su(2) coefficients
      symop.str_type=_string;symop.flag_inversion=_inversion;
      symop.str_Hermann_Mauguin=_stringHM;
      symop.str_Schoenflies=_stringSC;
      //GG
      SYM::CalculateQuaternion(symop);
      //GG
    }
    if(group==_PGROUP_) {
      a.pgroup_calculated=TRUE;
      clear(symop.ctau);clear(symop.ftau);clear(symop.ctrasl);clear(symop.ftrasl);              // no translation on point group
      symop.is_pgroup=TRUE;symop.is_fgroup=FALSE;symop.is_sgroup=FALSE;symop.is_agroup=FALSE;symop.is_pgroupk=FALSE;symop.is_pgroup_xtal=FALSE;symop.is_pgroupk_xtal=FALSE;symop.is_pgroupk_Patterson=FALSE;
      a.pgroup.push_back(symop);
      return a.pgroup.size();  // it returns the number of operations saved
    }
    if(group==_PGROUP_XTAL_) {
      a.pgroup_xtal_calculated=TRUE;
      clear(symop.ctau);clear(symop.ftau);clear(symop.ctrasl);clear(symop.ftrasl);              // no translation on point group
      symop.is_pgroup=FALSE;symop.is_fgroup=FALSE;symop.is_sgroup=FALSE;symop.is_agroup=FALSE;symop.is_pgroupk=FALSE;symop.is_pgroup_xtal=TRUE;symop.is_pgroupk_xtal=FALSE;symop.is_pgroupk_Patterson=FALSE;
      a.pgroup_xtal.push_back(symop);
      return a.pgroup_xtal.size();  // it returns the number of operations saved
    }
    if(group==_PGROUPK_PATTERSON_) { //DX20200129
      a.pgroupk_Patterson_calculated=TRUE;
      clear(symop.ctau);clear(symop.ftau);clear(symop.ctrasl);clear(symop.ftrasl);              // no translation on point group
      symop.is_pgroup=FALSE;symop.is_fgroup=FALSE;symop.is_sgroup=FALSE;symop.is_agroup=FALSE;symop.is_pgroupk=FALSE;symop.is_pgroup_xtal=FALSE;symop.is_pgroupk_xtal=FALSE;symop.is_pgroupk_Patterson=TRUE;
      a.pgroupk_Patterson.push_back(symop);
      return a.pgroupk_Patterson.size();  // it returns the number of operations saved
    }
    if(group==_FGROUP_) {
      a.fgroup_calculated=TRUE;
      symop.is_pgroup=FALSE;symop.is_fgroup=TRUE;symop.is_sgroup=FALSE;symop.is_agroup=FALSE;symop.is_pgroupk=FALSE;symop.is_pgroup_xtal=FALSE;symop.is_pgroupk_xtal=FALSE;symop.is_pgroupk_Patterson=FALSE;
      a.fgroup.push_back(symop);
      return a.fgroup.size();  // it returns the number of operations saved
    }
    if(group==_SGROUP_) {
      a.sgroup_calculated=TRUE;
      symop.is_pgroup=FALSE;symop.is_fgroup=FALSE;symop.is_sgroup=TRUE;symop.is_agroup=FALSE;symop.is_pgroupk=FALSE;symop.is_pgroup_xtal=FALSE;symop.is_pgroupk_xtal=FALSE;symop.is_pgroupk_Patterson=FALSE;
      a.sgroup.push_back(symop);
      return a.sgroup.size();  // it returns the number of operations saved
    }
    if(group==_AGROUP_) {
      a.agroup_calculated=TRUE;
      clear(symop.ctau);clear(symop.ftau);clear(symop.ctrasl);clear(symop.ftrasl);              // no translation on site point group
      symop.is_pgroup=FALSE;symop.is_fgroup=FALSE;symop.is_sgroup=FALSE;symop.is_agroup=TRUE;symop.is_pgroupk=FALSE;symop.is_pgroup_xtal=FALSE;symop.is_pgroupk_xtal=FALSE;symop.is_pgroupk_Patterson=FALSE;
      symop.site=iat; //DX20170803
      a.agroup.at(iat).push_back(symop);
      return a.agroup.at(iat).size();   // it returns the number of operations saved
    }
    if(group==_PGROUPK_) {
      a.pgroupk_calculated=TRUE;
      clear(symop.ctau);clear(symop.ftau);clear(symop.ctrasl);clear(symop.ftrasl);              // no translation on point group
      symop.is_pgroup=FALSE;symop.is_fgroup=FALSE;symop.is_sgroup=FALSE;symop.is_agroup=FALSE;symop.is_pgroupk=TRUE;symop.is_pgroup_xtal=FALSE;symop.is_pgroupk_xtal=FALSE;symop.is_pgroupk_Patterson=FALSE;
      a.pgroupk.push_back(symop);
      return a.pgroupk.size();  // it returns the number of operations saved
    }
    //DX20171205 - Added pgroupk_xtal - START
    if(group==_PGROUPK_XTAL_) {
      a.pgroupk_xtal_calculated=TRUE;
      clear(symop.ctau);clear(symop.ftau);clear(symop.ctrasl);clear(symop.ftrasl);              // no translation on point group
      symop.is_pgroup=FALSE;symop.is_fgroup=FALSE;symop.is_sgroup=FALSE;symop.is_agroup=FALSE;symop.is_pgroupk=FALSE;symop.is_pgroup_xtal=FALSE;symop.is_pgroupk_xtal=TRUE;symop.is_pgroupk_Patterson=FALSE;
      a.pgroupk_xtal.push_back(symop);
      return a.pgroupk_xtal.size();  // it returns the number of operations saved
    }
    //DX20171205 - Added pgroupk_xtal - END
    return (uint) 0;
  }
} // namespace SYM

// obvious variations
namespace SYM {
  //DX+CO START
  uint AddSymmetryToStructure(xstructure &a,const xmatrix<double> &Uc,const xmatrix<double> &Uf,const xvector<double> &ctau,const xvector<double> &ftau,
      const xvector<double> &ctrasl,const xvector<double> &ftrasl,const std::vector<int> &basis_atoms_map,const std::vector<int> &basis_types_map,bool basis_map_calculated,char group) {

    bool roff=TRUE;
    return AddSymmetryToStructure(a,Uc,Uf,ctau,ftau,ctrasl,ftrasl,basis_atoms_map,basis_types_map,basis_map_calculated,group,roff);
  }
  uint AddSymmetryToStructure(xstructure &a,const xmatrix<double> &Uc,const xmatrix<double> &Uf,const xvector<double> &ctau,const xvector<double> &ftau,
      const xvector<double> &ctrasl,const xvector<double> &ftrasl,const std::vector<int> &basis_atoms_map,const std::vector<int> &basis_types_map,bool basis_map_calculated,char group,bool roff) {
    return SYM::AddSymmetryToStructure(a,0,Uc,Uf,ctau,ftau,ctrasl,ftrasl,basis_atoms_map,basis_types_map,basis_map_calculated,group,roff);
    //DX+CO END
  }
} // namespace SYM

namespace SYM {
  //DX+CO START
  uint AddSymmetryToStructure(xstructure &a,const uint& iat,const xmatrix<double> &Uc,const xmatrix<double> &Uf,const std::vector<int> &basis_atoms_map,const std::vector<int> &basis_types_map,bool basis_map_calculated,char group) {
    bool roff=TRUE;
    return AddSymmetryToStructure(a,iat,Uc,Uf,basis_atoms_map,basis_types_map,basis_map_calculated,group,roff);
  }
  uint AddSymmetryToStructure(xstructure &a,const uint& iat,const xmatrix<double> &Uc,const xmatrix<double> &Uf,const std::vector<int> &basis_atoms_map,const std::vector<int> &basis_types_map,bool basis_map_calculated,char group,bool roff) {
    xvector<double> xxx(3);xxx.clear();
    return SYM::AddSymmetryToStructure(a,iat,Uc,Uf,xxx,xxx,xxx,xxx,basis_atoms_map,basis_types_map,basis_map_calculated,group,roff);
  }
  //DX+CO END
} // namespace SYM

namespace SYM {
  //DX+CO START
  uint AddSymmetryToStructure(xstructure &a,const xmatrix<double> &Uc,const xmatrix<double> &Uf,const std::vector<int> &basis_atoms_map,const std::vector<int> &basis_types_map,bool basis_map_calculated,char group) {
    bool roff=TRUE;
    return AddSymmetryToStructure(a,Uc,Uf,basis_atoms_map,basis_types_map,basis_map_calculated,group,roff);
  }
  uint AddSymmetryToStructure(xstructure &a,const xmatrix<double> &Uc,const xmatrix<double> &Uf,const std::vector<int> &basis_atoms_map,const std::vector<int> &basis_types_map,bool basis_map_calculated,char group,bool roff) {
    xvector<double> xxx(3);xxx.clear();
    return SYM::AddSymmetryToStructure(a,0,Uc,Uf,xxx,xxx,xxx,xxx,basis_atoms_map,basis_types_map,basis_map_calculated,group,roff);
    //DX+CO END
  }
} // namespace SYM

// ----------------------------------------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------------------- POINT GROUP OPERATIONS
// SYM::CalculatePointGroup
//
// This function calculate the whole point group and saves in the aflow.pgroup file
// Written by SC, aug 2007

//DX+CO START
namespace SYM {
  bool CalculatePointGroup(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format) {
    double _eps_=AUROSTD_NAN;//=_EPS_;
    if(a.sym_eps!=AUROSTD_NAN){ //Tolerance came from user or was calculated
      _eps_=a.sym_eps;
    }
    else {
      // Calculate the space group and store final tolerance
      // modify this stuff at your own risk - CO
      _eps_ = defaultTolerance(a);
      //DX [OBSOLETE] a.SpaceGroup_ITC(_eps_,-1);
      a.sym_eps=_eps_;
      a.sym_eps_calculated=true;
    }
    return SYM::CalculatePointGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss,_eps_,format);
  }
} // namespace SYM

namespace SYM {
  bool CalculatePointGroup(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss, double _eps_, string format) {
    return CalculatePointGroup_20160801(FileMESSAGE,a,aflags,_write_,osswrite,oss,_eps_,format);
  }
} // namespace SYM

namespace SYM {
  bool CalculatePointGroup_20160801(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format) { 
    bool LDEBUG=(FALSE || XHOST.DEBUG); //CO20190520
    stringstream message;
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}
    // Obtain the structure tolerance
    //DX20180526 [OBSOLETE] string directory=aurostd::execute2string("pwd"); //DX20180426 - added current working directory
    a.sym_eps=_eps_; //DX

    // AFLOW_FUNCTION_IMPLEMENTATION
    if(DEBUG_SYMMETRY) cerr << __AFLOW_FUNC__ << " DEBUG" << endl;
    // ------------------------------------------------------------------------------
    // Some of this routine is inspired by AVDV structure::calc_point_group() code.
    //SC made modificatios for speed and consistency with aflow architecture
    //a.MinkowskiBasisReduction();
    ostringstream aus;
    bool Krun=TRUE;
    string pgroup_type="PGROUP";
    if(a.title=="KLATTICE") pgroup_type="PGROUP_KLATTICE";
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << pgroup_type << " Symmetry: BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    // generic variables
    //  double _eps_=_EPS_; optional
    //SC ORIG double _eps_angle_=10.0*_eps_,_eps_abc_,_eps_vol_;
    bool is_deg=true;                                        // Angles found are in degrees //DX
    double _eps_vol_;
    bool sym_found;
    xvector<double> rrr(3),ddd(3);                    	     // for printing
    // ---------------------------------------------------------------------------
    // create placeholder for atoms_map for AddSymmetryToStructure() (no atom mappings for pgroups)
    std::vector<int> basis_atoms_map(a.atoms.size());        // will map each on each
    for(uint i=0;i<a.atoms.size();i++) basis_atoms_map[i]=i; // identically map each over each
    std::vector<int> basis_types_map(a.atoms.size());        // will map each on each
    for(uint i=0;i<a.atoms.size();i++) basis_types_map[i]=a.atoms[i].type; // identically map each over each
    // clean up pgroup                                       // just initialize
    a.pgroup.clear();
    a.pgroup_calculated=FALSE;                               // just initialize

    xmatrix<double> lattice(3,3),temp_clattice(3,3),temp_flattice(3,3);
    //  a.MinkowskiBasisReduction(); //must
    // (JJPR BUG) if(a.LatticeReduction_avoid==FALSE) a.LatticeReduction(); // only once
    //a.FixLattices();a.BringInCell(); //JJPR BUG             // bring atoms in new basis  // useless but useful JJPR BUG
    xmatrix<double> f2c(3,3);f2c=trasp(a.lattice);
    xmatrix<double> c2f(3,3);c2f=inverse(trasp(a.lattice));
#ifdef _NO_SCALE_LATTICE_PGROUP_
    lattice=(a.lattice);                                  // copies so it does not messup
#else
    lattice=a.scale*(a.lattice);                          // copies so it does not messup
#endif
    xmatrix<double> Act(3,3),Bct(3,3),Uc(3,3);               // for generation of pgroup
    xmatrix<double> Adt(3,3),Bdt(3,3),Uf(3,3);               // for generation of pgroup_ijk
    xvector<double> a1(3),a2(3),a3(3);                       // lattice vectors
    xvector<double> ctau(3),ftau(3),ctrasl(3),ftrasl(3);     // translation
    a1=lattice(1);a2=lattice(2);a3=lattice(3);               // a1,a2,a3 are the rows of the lattice matrix
    Act=lattice;
    Adt=aurostd::identity((double) 0,3);  //ME20200123 - new identity format
    xvector<double> clatticedata(6),temp_clatticedata(6);    // triplet position
    clatticedata=_Getabc_angles(lattice,DEGREES);
    vector<xvector<double>*> grid_clattice; xvector<double> *grid_clattice_ptr;  // grid for rrr points
    vector<xvector<double>*> grid_flattice; xvector<double> *grid_flattice_ptr;  // grid for ijk points
    double na1,na2,na3,amin,latticevol,temp_clatticevol;
    // [UNUSED] double amax; // warning: variable ‘amax’ set but not used
    na1=modulus(a1);na2=modulus(a2);na3=modulus(a3);
    // [UNUSED] amax=max(na1,na2,na3);  // warning: variable ‘amax’ set but not used
    amin=min(na1,na2,na3);
    //SC ORIG _eps_abc_=_eps_*amin;           // relative
    _eps_vol_=_eps_*amin*amin*amin; // relative
    latticevol=abs(det(lattice));
    // ------------------------------------------------------------------------------
    // make a lattice parallelepiped that encompasses a sphere with
    // radius = 3.0*largest latparam dimensions for a sphere
    xvector<int> dims(3);
    double radius=0;
    // radius=3.0*max(na1,na2,na3);
    //  radius=na1+na2+na3;
    radius=RadiusSphereLattice(lattice);  // smarter !
    dims=LatticeDimensionSphere(lattice,radius);
    // for(i=1;i<=3;i++) if(dims(i)>3)dims(i)=3; // check for safety
    //DX ===
    //for(uint i=1;i<=3;i++){
    //  dims(i)=3;
    //}
    //DX ===
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << pgroup_type << " Symmetry: inside sphere, dimensions = [" << dims[1] << "," << dims[2] << "," << dims[3] << "]   radius=" << radius << "  " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    // ------------------------------------------------------------------------------
    // seek for lattice points within the sphere with radius
    grid_clattice.clear();grid_flattice.clear();
    // insert r1,r2,r3
    for(uint i=1;i<=3;i++) {
      grid_clattice_ptr = new xvector<double>(3);
      *grid_clattice_ptr=Act(i);
      grid_clattice.push_back(grid_clattice_ptr);  // TRICK so I get the identity first
      grid_flattice_ptr = new xvector<double>(3);
      *grid_flattice_ptr=Adt(i);
      grid_flattice.push_back(grid_flattice_ptr);  // TRICK so I get the identity first
    }
    //  cerr << I3 << endl;
    for(int i=-dims[1];i<=dims[1];i++){
      for(int j=-dims[2];j<=dims[2];j++){
        for(int k=-dims[3];k<=dims[3];k++){
          if(!(i==0 && j==0 && k==0)){
            if((!(i==1 && j==0 && k==0)) && (!(i==0 && j==1 && k==0)) && (!(i==0 && j==0 && k==1))) { // to avoid plugging a1,a2,a2 back
              rrr=((double)i)*a1+((double)j)*a2+((double)k)*a3;
              ddd(1)=(double) i;ddd(2)=(double) j;ddd(3)=(double) k;
              //keep only the lattice points within the sphere with radius
              if(modulus(rrr)<radius && 
                  //DX START : Speed up vectors to check
                  //((modulus(rrr)>na1_min && modulus(rrr)<na1_max) || 
                  // (modulus(rrr)>na2_min && modulus(rrr)<na2_max) || 
                  // (modulus(rrr)>na3_min && modulus(rrr)<na3_max)))
                (aurostd::abs(modulus(rrr)-na1) < _eps_ ||
                 aurostd::abs(modulus(rrr)-na2) < _eps_ ||
                 aurostd::abs(modulus(rrr)-na3) < _eps_))
                 {  //CO20200106 - patching for auto-indenting
                   //DX END : Speed up vectors to check
                   grid_clattice_ptr = new xvector<double>(3);        // SAVE THEM ALL
                   *grid_clattice_ptr=rrr;                            // SAVE THEM ALL
                   grid_clattice.push_back(grid_clattice_ptr);        // SAVE THEM ALL
                   grid_flattice_ptr = new xvector<double>(3);        // SAVE THEM ALL
                   *grid_flattice_ptr=ddd;                            // SAVE THEM ALL
                   grid_flattice.push_back(grid_flattice_ptr);        // SAVE THEM ALL
                 }
            }
          }
        }
      }
    }
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << pgroup_type << " Symmetry: found " << grid_clattice.size() << " vectors inside radius" << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    // ------------------------------------------------------------------------------
    // for each set of three lattice points within the sphere see which one has the
    // same sets of lengths and angles as the original lattice unit cell vectors.
    //[CO20190520 - testing setUf and setUc]_sym_op test; //delete me
    //[CO20190520 - testing setUf and setUc]bool isequal; //delete me
    uint ii=0,jj=0,kk=0;
    for(uint i=0;i<grid_clattice.size();i++){
      for(ii=1;ii<=3;ii++) temp_clattice(1,ii)=(*grid_clattice[i])(ii);
      for(ii=1;ii<=3;ii++) temp_flattice(1,ii)=(*grid_flattice[i])(ii);
      for(uint j=0;j<grid_clattice.size();j++){
        for(ii=1;ii<=3;ii++) temp_clattice(2,ii)=(*grid_clattice[j])(ii);
        for(ii=1;ii<=3;ii++) temp_flattice(2,ii)=(*grid_flattice[j])(ii);
        for(uint k=0;k<grid_clattice.size();k++){	
          if(i!=j && i!=k && j!=k){
            for(ii=1;ii<=3;ii++) temp_clattice(3,ii)=(*grid_clattice[k])(ii);
            for(ii=1;ii<=3;ii++) temp_flattice(3,ii)=(*grid_flattice[k])(ii);
            temp_clatticevol=abs(det(temp_clattice));
            if(abs(temp_clatticevol-latticevol)<=_eps_vol_ && temp_clatticevol > 0.0) {  // check volume !!! FAST !
              temp_clatticedata=_Getabc_angles(temp_clattice,DEGREES);
              jj++;//if(!mod(jj,100000)) cerr << temp_clatticedata << endl;
              // compare the temp_clattice... and lattice... to see if they are the same lattice
              // that is do the lattice vectors have the same lengths and the same angles
              if(abs(clatticedata[1]-temp_clatticedata[1]) < _eps_)                // check lattices
                if(abs(clatticedata[2]-temp_clatticedata[2]) < _eps_)              // check lattices
                  if(abs(clatticedata[3]-temp_clatticedata[3]) < _eps_)            // check lattices
                    //DX =======
                    if(checkAngle(clatticedata[2],clatticedata[3],clatticedata[4],temp_clatticedata[4],is_deg,_eps_))        // check angles //DX
                      if(checkAngle(clatticedata[3],clatticedata[1],clatticedata[5],temp_clatticedata[5],is_deg,_eps_))      // check angles //DX
                        if(checkAngle(clatticedata[1],clatticedata[2],clatticedata[6],temp_clatticedata[6],is_deg,_eps_)) {  // check angles //DX
                          //DX =======
                          // SYMMETRY OPERATION POINT_GROUP and FACTOR_GROUP
                          // get the matrix that relates the two lattice vectors
                          // A=U*B but in A and B we plug vectors as columns
                          // watch lattice is per row Uc=A*inv(B)
                          // A is the lattice (vectors per colum), B is the test lattice (epr column)
                          // Uc is the point_group operation which operates AFTER the vector (row)
                          // as: new_vector_row=old_vector_row*Uc
                          // point_group is the list of all the Uc !!!
                          // since A and B are saved as trasposed, I get
                          // U=trasp(A)*inv(trasp(U))
                          Bct=temp_clattice;Uc=trasp(Act)*inverse(trasp(Bct)); // Act=lattice;
                          Bdt=temp_flattice;Uf=trasp(Adt)*inverse(trasp(Bdt)); // Adt=I3... could remove !

                          //[CO20190520 - testing setUf and setUc]if(LDEBUG) { //delete me
                          //[CO20190520 - testing setUf and setUc]  test.setUf(Uf,lattice);
                          //[CO20190520 - testing setUf and setUc]  isequal=aurostd::isequal(Uc,test.Uc);
                          //[CO20190520 - testing setUf and setUc]  if(!isequal){cerr << "BAD2!" << endl;}

                          //[CO20190520 - testing setUf and setUc]  cerr << "Uc=" << endl;cerr << Uc << endl;
                          //[CO20190520 - testing setUf and setUc]  //cerr << "c2f*Uc*f2c=" << endl;cerr << c2f*Uc*f2c << endl;
                          //[CO20190520 - testing setUf and setUc]  //cerr << "f2c*Uc*c2f=" << endl;cerr << f2c*Uc*c2f << endl;
                          //[CO20190520 - testing setUf and setUc]  cerr << "test.Uc=" << endl;cerr << test.Uc << endl;

                          //[CO20190520 - testing setUf and setUc]  test.setUc(Uc,lattice);
                          //[CO20190520 - testing setUf and setUc]  isequal=aurostd::isequal(Uf,test.Uf);
                          //[CO20190520 - testing setUf and setUc]  if(!isequal){cerr << "BAD1!" << endl;}

                          //[CO20190520 - testing setUf and setUc]  cerr << "Uf=" << endl;cerr << Uf << endl;
                          //[CO20190520 - testing setUf and setUc]  //cerr << "c2f*Uf*f2c=" << endl;cerr << c2f*Uf*f2c << endl;
                          //[CO20190520 - testing setUf and setUc]  //cerr << "f2c*Uf*c2f=" << endl;cerr << f2c*Uf*c2f << endl;
                          //[CO20190520 - testing setUf and setUc]  cerr << "test.Uf=" << endl;cerr << test.Uf << endl;
                          //[CO20190520 - testing setUf and setUc]}

                          // check whether this symmetry operation is new or not
                          sym_found=FALSE;
                          for(ii=0;ii<a.pgroup.size()&&!sym_found;ii++){
                            sym_found=identical(Uf,a.pgroup[ii].Uf);       // look in all the list of operations  //DX20171207 - Use Uf (only integers) not Uc and use xmatrix identical eps
                          }
                          // if the symmetry operation is new, add it to the pointgroup array
                          // and update all info about the sym_op object
                          if(sym_found==FALSE) {                                 // new operation, generate and save it
                            clear(ctau);clear(ftau);clear(ctrasl);clear(ftrasl);
                            kk=SYM::AddSymmetryToStructure(a,Uc,Uf,ctau,ftau,ctrasl,ftrasl,basis_atoms_map,basis_types_map,false,_PGROUP_,FALSE);  // kk number pgroups
                            aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << pgroup_type << " " << a.pgroup[kk-1].str_type
                              << " theta=";
                            if(a.pgroup[kk-1].angle<100) aus << " ";
                            if(a.pgroup[kk-1].angle<10) aus << " ";
                            aus << a.pgroup[kk-1].angle << " " << " r=(" << a.pgroup[kk-1].axis << ")";
                            aus << "    HM=" <<  aurostd::PaddedPRE(a.pgroup[kk-1].str_Hermann_Mauguin,2," ");     // remember vectors start from 0
                            aus << "    S=" <<  aurostd::PaddedPRE(a.pgroup[kk-1].str_Schoenflies,2," ");          // remember vectors start from 0
                            aus << endl;  // remember vectors start from 0
                            //			    aus << "JJPR =" << Uf << endl;
                            if(0) { // to fix the JJPR BUG
                              if(a.pgroup.size()==28) {
                                cerr << clatticedata[1] << " " << clatticedata[2] << " " << clatticedata[3] << " "
                                  << clatticedata[4] << " " << clatticedata[5] << " " << clatticedata[6] << " " << endl;
                                cerr << "Adt" << endl << Adt << endl;
                                cerr << "Act" << endl << Act << endl;
                                cerr << "Bdt" << endl << Bdt << endl;
                                cerr << "Bct" << endl << Bct << endl;
                                cerr << "Uf" << endl << Uf << endl;
                                cerr << "Uc" << endl << Uc << endl;
                              }
                            }
                            //DX _eps_ to _ZERO_TOL_
                            if(aurostd::sum(aurostd::abs(f2c*Uf*inverse(f2c)-Uc))>_ZERO_TOL_) //DX used to be _eps_
                            {message << "Uf error[1] [dir=" << a.directory << "]"; throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);}
                            if(aurostd::sum(aurostd::abs(c2f*Uc*inverse(c2f)-Uf))>_ZERO_TOL_) //DX used to be _eps_
                            {message << "Uc error[2] [dir=" << a.directory << "]"; throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);}
                            if(aurostd::sum(aurostd::abs(FF2CC(a.lattice,Uf)-Uc))>_ZERO_TOL_) //DX used to be _eps_
                            {message << "Uc error[3] [dir=" << a.directory << "]"; throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);}
                            if(aurostd::sum(aurostd::abs(CC2FF(a.lattice,Uc)-Uf))>_ZERO_TOL_) //DX used to be _eps_
                            {message << "Uc error[4] [dir=" << a.directory << "]"; throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);}
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
                            //DX _eps_ to _ZERO_TOL_
                          }
                        }
            }
          }
        }
      }
    }
    // ------------------------------------------------------------------------------
    a.pgroup_calculated=TRUE;
    // ------------------------------------------------------------------------------
    //for(uint k=0;k<a.pgroup.size();k++) {
    //  Uc=a.pgroup[k].Uc;
    //  Uf=a.pgroup[k].Uf;
    //roundoff(a.pgroup[k].Uc,_EPS_roundoff_); //CO+DX COME BACK
    //for(i=1;i<=3;i++)
    //  for(j=1;j<=3;j++)
    //	if(abs((*a.pgroup[k])(i,j))<_EPS_roundoff_) (*a.pgroup[k])(i,j)=0.0;
    //for(i=0;i<9;i++) if(abs((*a.pgroup[k])(int(i/3)+1,mod(i,3)+1))<1.0e-10) (*a.pgroup[k])(int(i/3)+1,mod(i,3)+1)=0.0;
    //}
    // ------------------------------------------------------------------------------
    // going out, destroy useless things
    for(uint i=0;i<grid_clattice.size();i++){ 
      delete grid_clattice[i]; 
    }
    grid_clattice.clear();
    for(uint i=0;i<grid_flattice.size();i++){ 
      delete grid_flattice[i]; 
    }
    grid_flattice.clear();
    // ------------------------------------------------------------------------------
    // Printing and leaving
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << pgroup_type << " Symmetry: unique point group operations " << a.pgroup.size() << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << pgroup_type << " Symmetry: END " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    if(_write_) Krun=Krun && KBIN_SymmetryWrite(FileMESSAGE,a,aflags,_PGROUP_,osswrite,oss,format);
    string pgname = ""; //DX20170906
    string operations = ""; //DX20170906
    return PointGroupMap(a, pgname, operations, _PGROUP_); //DX20170906
  }
} // namespace SYM

//DX+CO END

#ifndef COMPILE_SLIM
namespace SYM {
  //DX bool CalculatePointGroup(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_)        // AFLOW_FUNCTION_IMPLEMENTATION
  bool CalculatePointGroup_20160101(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_)        // AFLOW_FUNCTION_IMPLEMENTATION
  { //CO20200106 - patching for auto-indenting
    stringstream message;
    if(DEBUG_SYMMETRY) cerr << __AFLOW_FUNC__ << " DEBUG" << endl;
    // ------------------------------------------------------------------------------
    // Some of this routine is inspired by AVDV structure::calc_point_group() code.
    //SC made modificatios for speed and consistency with aflow architecture
    //a.MinkowskiBasisReduction();
    ostringstream aus;
    bool Krun=TRUE;
    string pgroup_type="PGROUP";
    if(a.title=="KLATTICE") pgroup_type="PGROUP_KLATTICE";
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << pgroup_type << " Symmetry: BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    // generic variables
    //  double _eps_=_EPS_; optional
    double _eps_angle_=10.0*_eps_,_eps_abc_,_eps_vol_;
    bool sym_found;
    xvector<double> rrr(3),ddd(3);                    // for printing
    std::vector<int> basis_atoms_map(a.atoms.size());       // will map each on each
    for(uint i=0;i<a.atoms.size();i++) basis_atoms_map[i]=i; // identically map each over each
    std::vector<int> basis_types_map(a.atoms.size());       // will map each on each
    for(uint i=0;i<a.atoms.size();i++) basis_types_map[i]=a.atoms[i].type; // identically map each over each
    // clean up pgroup                                // just initialize
    a.pgroup.clear();
    a.pgroup_calculated=FALSE;                        // just initialize

    xmatrix<double> lattice(3,3),temp_clattice(3,3),temp_flattice(3,3);
    //  a.MinkowskiBasisReduction(); //must
    // (JJPR BUG) if(a.LatticeReduction_avoid==FALSE) a.LatticeReduction(); // only once
    a.FixLattices();a.BringInCell(); //JJPR BUG                       // bring atoms in new basis  // useless but useful JJPR BUG
    xmatrix<double> f2c(3,3);f2c=trasp(a.lattice);
    xmatrix<double> c2f(3,3);c2f=inverse(trasp(a.lattice));
#ifdef _NO_SCALE_LATTICE_PGROUP_
    lattice=(a.lattice);                                  // copies so it does not messup
#else
    lattice=a.scale*(a.lattice);                          // copies so it does not messup
#endif
    xmatrix<double> Act(3,3),Bct(3,3),Uc(3,3);            // for generation of pgroup
    xmatrix<double> Adt(3,3),Bdt(3,3),Uf(3,3);            // for generation of pgroup_ijk
    xvector<double> a1(3),a2(3),a3(3);                    // lattice vectors
    xvector<double> ctau(3),ftau(3),ctrasl(3),ftrasl(3);  // translation
    a1=lattice(1);a2=lattice(2);a3=lattice(3);            // a1,a2,a3 are the rows of the lattice matrix
    Act=lattice;
    Adt=aurostd::identity((double) 0,3);  //ME20200123 - new identity format
    xvector<double> clatticedata(6),temp_clatticedata(6); // triplet position
    clatticedata=_Getabc_angles(lattice,DEGREES);
    vector<xvector<double>*> grid_clattice; xvector<double> *grid_clattice_ptr;  // grid for rrr points
    vector<xvector<double>*> grid_flattice; xvector<double> *grid_flattice_ptr;  // grid for ijk points
    double na1,na2,na3,amin,latticevol,temp_clatticevol;
    // [UNUSED] double amax; // warning: variable ‘amax’ set but not used
    na1=modulus(a1);na2=modulus(a2);na3=modulus(a3);
    // [UNUSED] amax=max(na1,na2,na3);  // warning: variable ‘amax’ set but not used
    amin=min(na1,na2,na3);
    _eps_abc_=_eps_*amin;           // relative
    _eps_vol_=_eps_*amin*amin*amin; // relative
    latticevol=abs(det(lattice));
    // ------------------------------------------------------------------------------
    // make a lattice parallelepiped that encompasses a sphere with
    // radius = 3.0*largest latparam dimensions for a sphere
    xvector<int> dims(3);
    double radius=0;
    // radius=3.0*max(na1,na2,na3);
    //  radius=na1+na2+na3;
    radius=RadiusSphereLattice(lattice);  // smarter !
    dims=LatticeDimensionSphere(lattice,radius);
    // for(i=1;i<=3;i++) if(dims(i)>3)dims(i)=3; // check for safety
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << pgroup_type << " Symmetry: inside sphere, dimensions = [" << dims[1] << "," << dims[2] << "," << dims[3] << "]   radius=" << radius << "  " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    // ------------------------------------------------------------------------------
    // seek for lattice points within the sphere with radius
    grid_clattice.clear();grid_flattice.clear();
    // insert r1,r2,r3
    for(uint i=1;i<=3;i++) {
      grid_clattice_ptr = new xvector<double>(3);
      *grid_clattice_ptr=Act(i);
      grid_clattice.push_back(grid_clattice_ptr);  // TRICK so I get the identity first
      grid_flattice_ptr = new xvector<double>(3);
      *grid_flattice_ptr=Adt(i);
      grid_flattice.push_back(grid_flattice_ptr);  // TRICK so I get the identity first
    }
    //  cerr << I3 << endl;
    for(int i=-dims[1];i<=dims[1];i++) {
      for(int j=-dims[2];j<=dims[2];j++) {
        for(int k=-dims[3];k<=dims[3];k++) {
          if(!(i==0 && j==0 && k==0)) {
            if((!(i==1 && j==0 && k==0)) && (!(i==0 && j==1 && k==0)) && (!(i==0 && j==0 && k==1))) { // to avoid plugging a1,a2,a2 back
              rrr=((double)i)*a1+((double)j)*a2+((double)k)*a3;
              ddd(1)=(double) i;ddd(2)=(double) j;ddd(3)=(double) k;
              //keep only the lattice points within the sphere with radius
              if(modulus(rrr)<radius) {
                grid_clattice_ptr = new xvector<double>(3);        // SAVE THEM ALL
                *grid_clattice_ptr=rrr;                            // SAVE THEM ALL
                grid_clattice.push_back(grid_clattice_ptr);        // SAVE THEM ALL
                grid_flattice_ptr = new xvector<double>(3);        // SAVE THEM ALL
                *grid_flattice_ptr=ddd;                            // SAVE THEM ALL
                grid_flattice.push_back(grid_flattice_ptr);        // SAVE THEM ALL
              }
            }
          }
        }
      }
    }
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << pgroup_type << " Symmetry: found " << grid_clattice.size() << " vectors inside radius" << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    // ------------------------------------------------------------------------------
    // for each set of three lattice points within the sphere see which one has the
    // same sets of lengths and angles as the original lattice unit cell vectors.
    uint ii=0,jj=0,kk=0;
    for(uint i=0;i<grid_clattice.size();i++) {
      for(ii=1;ii<=3;ii++) temp_clattice(1,ii)=(*grid_clattice[i])(ii);
      for(ii=1;ii<=3;ii++) temp_flattice(1,ii)=(*grid_flattice[i])(ii);
      for(uint j=0;j<grid_clattice.size();j++) {
        for(ii=1;ii<=3;ii++) temp_clattice(2,ii)=(*grid_clattice[j])(ii);
        for(ii=1;ii<=3;ii++) temp_flattice(2,ii)=(*grid_flattice[j])(ii);
        for(uint k=0;k<grid_clattice.size();k++) {	
          if(i!=j && i!=k && j!=k) {
            for(ii=1;ii<=3;ii++) temp_clattice(3,ii)=(*grid_clattice[k])(ii);
            for(ii=1;ii<=3;ii++) temp_flattice(3,ii)=(*grid_flattice[k])(ii);
            temp_clatticevol=abs(det(temp_clattice));
            if(abs(temp_clatticevol-latticevol)<=_eps_vol_) {  // check volume !!! FAST !
              temp_clatticedata=_Getabc_angles(temp_clattice,DEGREES);
              jj++;//if(!mod(jj,100000)) cerr << temp_clatticedata << endl;
              // compare the temp_clattice... and lattice... to see if they are the same lattice
              // that is do the lattice vectors have the same lengths and the same angles
              if(abs(clatticedata[1]-temp_clatticedata[1]) < _eps_abc_)                // check lattices
                if(abs(clatticedata[2]-temp_clatticedata[2]) < _eps_abc_)              // check lattices
                  if(abs(clatticedata[3]-temp_clatticedata[3]) < _eps_abc_)            // check lattices
                    if(abs(clatticedata[4]-temp_clatticedata[4]) < _eps_angle_)        // check angles
                      if(abs(clatticedata[5]-temp_clatticedata[5]) < _eps_angle_)      // check angles
                        if(abs(clatticedata[6]-temp_clatticedata[6]) < _eps_angle_) {  // check angles
                          // SYMMETRY OPERATION POINT_GROUP and FACTOR_GROUP
                          // get the matrix that relates the two lattice vectors
                          // A=U*B but in A and B we plug vectors as columns
                          // watch lattice is per row Uc=A*inv(B)
                          // A is the lattice (vectors per colum), B is the test lattice (epr column)
                          // Uc is the point_group operation which operates AFTER the vector (row)
                          // as: new_vector_row=old_vector_row*Uc
                          // point_group is the list of all the Uc !!!
                          // since A and B are saved as trasposed, I get
                          // U=trasp(A)*inv(trasp(U))
                          Bct=temp_clattice;Uc=trasp(Act)*inverse(trasp(Bct)); // Act=lattice;
                          Bdt=temp_flattice;Uf=trasp(Adt)*inverse(trasp(Bdt)); // Adt=I3... could remove !
                          // check whether this symmetry operation is new or not
                          sym_found=FALSE;
                          for(ii=0;ii<a.pgroup.size()&&!sym_found;ii++)
                            sym_found=identical(Uc,a.pgroup[ii].Uc,_eps_);       // look in all the list of operations
                          // if the symmetry operation is new, add it to the pointgroup array
                          // and update all info about the sym_op object
                          if(sym_found==FALSE) {                                 // new operation, generate and save it
                            clear(ctau);clear(ftau);clear(ctrasl);clear(ftrasl);
                            kk=SYM::AddSymmetryToStructure(a,Uc,Uf,ctau,ftau,ctrasl,ftrasl,basis_atoms_map,basis_types_map,false,_PGROUP_);  // kk number pgroups
                            aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << pgroup_type << " " << a.pgroup[kk-1].str_type
                              << " theta=";
                            if(a.pgroup[kk-1].angle<100) aus << " ";
                            if(a.pgroup[kk-1].angle<10) aus << " ";
                            aus << a.pgroup[kk-1].angle << " " << " r=(" << a.pgroup[kk-1].axis << ")";
                            aus << "    HM=" <<  aurostd::PaddedPRE(a.pgroup[kk-1].str_Hermann_Mauguin,2," ");     // remember vectors start from 0
                            aus << "    S=" <<  aurostd::PaddedPRE(a.pgroup[kk-1].str_Schoenflies,2," ");          // remember vectors start from 0
                            aus << endl;  // remember vectors start from 0
                            //			    aus << "JJPR =" << Uf << endl;
                            if(0) { // to fix the JJPR BUG
                              if(a.pgroup.size()==28) {
                                cerr << clatticedata[1] << " " << clatticedata[2] << " " << clatticedata[3] << " "
                                  << clatticedata[4] << " " << clatticedata[5] << " " << clatticedata[6] << " " << endl;
                                cerr << "Adt" << endl << Adt << endl;
                                cerr << "Act" << endl << Act << endl;
                                cerr << "Bdt" << endl << Bdt << endl;
                                cerr << "Bct" << endl << Bct << endl;
                                cerr << "Uf" << endl << Uf << endl;
                                cerr << "Uc" << endl << Uc << endl;
                              }
                            }
                            if(aurostd::sum(aurostd::abs(f2c*Uf*inverse(f2c)-Uc))>_eps_)
                            {message << "Uf error[1] [dir=" << a.directory << "]"; throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);}
                            if(aurostd::sum(aurostd::abs(c2f*Uc*inverse(c2f)-Uf))>_eps_)
                            {message << "Uf error[2] [dir=" << a.directory << "]"; throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);}
                            if(aurostd::sum(aurostd::abs(FF2CC(a.lattice,Uf)-Uc))>_eps_)
                            {message << "Uf error[3] [dir=" << a.directory << "]"; throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);}
                            if(aurostd::sum(aurostd::abs(CC2FF(a.lattice,Uc)-Uf))>_eps_)
                            {message << "Uf error[4] [dir=" << a.directory << "]"; throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);}
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
                          }
                        }
            }
          }
        }
      }
    }
    // ------------------------------------------------------------------------------
    a.pgroup_calculated=TRUE;
    // ------------------------------------------------------------------------------
    for(uint k=0;k<a.pgroup.size();k++) {
      Uc=a.pgroup[k].Uc;
      Uf=a.pgroup[k].Uf;
      roundoff(a.pgroup[k].Uc,_EPS_roundoff_);
      //for(i=1;i<=3;i++)
      //  for(j=1;j<=3;j++)
      //	if(abs((*a.pgroup[k])(i,j))<_EPS_roundoff_) (*a.pgroup[k])(i,j)=0.0;
      //for(i=0;i<9;i++) if(abs((*a.pgroup[k])(int(i/3)+1,mod(i,3)+1))<1.0e-10) (*a.pgroup[k])(int(i/3)+1,mod(i,3)+1)=0.0;
    }
    // ------------------------------------------------------------------------------
    // going out, destroy useless things
    for(uint i=0;i<grid_clattice.size();i++)
      delete grid_clattice[i];
    grid_clattice.clear();
    for(uint i=0;i<grid_flattice.size();i++)
      delete grid_flattice[i];
    grid_flattice.clear();
    // ------------------------------------------------------------------------------
    // Printing and leaving
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << pgroup_type << " Symmetry: unique point group operations " << a.pgroup.size() << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << pgroup_type << " Symmetry: END " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    if(_write_) Krun=Krun && KBIN_SymmetryWrite(FileMESSAGE,a,aflags,_PGROUP_,osswrite,oss);
    return Krun;
  }
} // namespace SYM
#endif

int PointGroup_HITS(xstructure &a,double _eps_) {
  int count=0;
  string message="PGROUP";
  if(a.title=="KLATTICE") message="PGROUP_KLATTICE";
  double _eps_angle_=10.0*_eps_,_eps_abc_,_eps_vol_;
  xvector<double> rrr(3),ddd(3);
  std::vector<int> basis_atoms_map(a.atoms.size());
  for(uint i=0;i<a.atoms.size();i++) basis_atoms_map[i]=i;
  std::vector<int> basis_types_map(a.atoms.size());
  for(uint i=0;i<a.atoms.size();i++) basis_types_map[i]=a.atoms[i].type;
  a.pgroup.clear();
  a.pgroup_calculated=FALSE;
  xmatrix<double> lattice(3,3),temp_clattice(3,3),temp_flattice(3,3);
  // (JJPR BUG) if(a.LatticeReduction_avoid==FALSE) a.LatticeReduction();
  xmatrix<double> f2c(3,3);f2c=trasp(a.lattice);
  xmatrix<double> c2f(3,3);c2f=inverse(trasp(a.lattice));
#ifdef _NO_SCALE_LATTICE_PGROUP_
  lattice=(a.lattice);
#else
  lattice=a.scale*(a.lattice);
#endif
  xmatrix<double> Act(3,3),Bct(3,3),Uc(3,3);
  xmatrix<double> Adt(3,3),Bdt(3,3),Uf(3,3);
  xvector<double> a1(3),a2(3),a3(3);
  xvector<double> ctau(3),ftau(3),ctrasl(3),ftrasl(3);
  a1=lattice(1);a2=lattice(2);a3=lattice(3);
  Act=lattice;
  Adt=aurostd::identity((double) 0,3);  //ME20200123 - new identity format
  xvector<double> clatticedata(6),temp_clatticedata(6);
  clatticedata=_Getabc_angles(lattice,DEGREES);
  vector<xvector<double>*> grid_clattice; xvector<double> *grid_clattice_ptr;
  vector<xvector<double>*> grid_flattice; xvector<double> *grid_flattice_ptr;
  double na1,na2,na3,amin,latticevol,temp_clatticevol;
  // [UNUSED] double amax; // warning: variable ‘amax’ set but not used
  na1=modulus(a1);na2=modulus(a2);na3=modulus(a3);
  // [UNUSED] amax=max(na1,na2,na3);  // warning: variable ‘amax’ set but not used
  amin=min(na1,na2,na3);
  _eps_abc_=_eps_*amin;
  _eps_vol_=_eps_*amin*amin*amin;
  latticevol=abs(det(lattice));
  xvector<int> dims(3);
  double radius=0;
  radius=RadiusSphereLattice(lattice);
  dims=LatticeDimensionSphere(lattice,radius);
  // --------------------------------------------
  grid_clattice.clear();grid_flattice.clear();
  for(uint i=1;i<=3;i++) {
    grid_clattice_ptr = new xvector<double>(3);
    *grid_clattice_ptr=Act(i);
    grid_clattice.push_back(grid_clattice_ptr);
    grid_flattice_ptr = new xvector<double>(3);
    *grid_flattice_ptr=Adt(i);
    grid_flattice.push_back(grid_flattice_ptr);
  }
  for(int i=-dims[1];i<=dims[1];i++) {
    for(int j=-dims[2];j<=dims[2];j++) {
      for(int k=-dims[3];k<=dims[3];k++) {
        if(!(i==0 && j==0 && k==0)) {
          if((!(i==1 && j==0 && k==0)) && (!(i==0 && j==1 && k==0)) && (!(i==0 && j==0 && k==1))) {
            rrr=((double)i)*a1+((double)j)*a2+((double)k)*a3;
            ddd(1)=(double) i;ddd(2)=(double) j;ddd(3)=(double) k;
            if(modulus(rrr)<radius) {
              grid_clattice_ptr = new xvector<double>(3);
              *grid_clattice_ptr=rrr;
              grid_clattice.push_back(grid_clattice_ptr);
              grid_flattice_ptr = new xvector<double>(3);
              *grid_flattice_ptr=ddd;
              grid_flattice.push_back(grid_flattice_ptr);
            }
          }
        }
      }
    }
  }
  uint ii=0,jj=0;
  for(uint i=0;i<grid_clattice.size();i++) {
    for(ii=1;ii<=3;ii++) temp_clattice(1,ii)=(*grid_clattice[i])(ii);
    for(ii=1;ii<=3;ii++) temp_flattice(1,ii)=(*grid_flattice[i])(ii);
    for(uint j=0;j<grid_clattice.size();j++) {
      for(ii=1;ii<=3;ii++) temp_clattice(2,ii)=(*grid_clattice[j])(ii);
      for(ii=1;ii<=3;ii++) temp_flattice(2,ii)=(*grid_flattice[j])(ii);
      for(uint k=0;k<grid_clattice.size();k++) {
        if(i!=j && i!=k && j!=k) {
          for(ii=1;ii<=3;ii++) temp_clattice(3,ii)=(*grid_clattice[k])(ii);
          for(ii=1;ii<=3;ii++) temp_flattice(3,ii)=(*grid_flattice[k])(ii);
          temp_clatticevol=abs(det(temp_clattice));
          if(abs(temp_clatticevol-latticevol)<=_eps_vol_) {
            temp_clatticedata=_Getabc_angles(temp_clattice,DEGREES);
            jj++;
            if(abs(clatticedata[1]-temp_clatticedata[1]) < _eps_abc_)
              if(abs(clatticedata[2]-temp_clatticedata[2]) < _eps_abc_)
                if(abs(clatticedata[3]-temp_clatticedata[3]) < _eps_abc_)
                  if(abs(clatticedata[4]-temp_clatticedata[4]) < _eps_angle_)
                    if(abs(clatticedata[5]-temp_clatticedata[5]) < _eps_angle_)
                      if(abs(clatticedata[6]-temp_clatticedata[6]) < _eps_angle_) {
                        count++;
                      }
          }
        }
      }
    }
  }
  return count;
}

void  PointGroupHistogramCheck() {
  double b=0;
  xstructure a(cin,IOAFLOW_AUTO);
  while(b<0.101) {
    cout<<b<<" ";
    cout<<PointGroup_HITS(a,b)<<endl;
    b+=0.001;
  }
}
vector<double> PointGroupHistogramCheck(xstructure& a) {
  int record=0;
  double tolerance=0;
  vector<double> result;
  result.push_back(0.001);
  while(tolerance<0.101) {
    if(record&&record!=PointGroup_HITS(a,tolerance)) {
      result.push_back(tolerance);
      record=PointGroup_HITS(a,tolerance);
    }
    tolerance+=0.001;
    if(!record) record=PointGroup_HITS(a,tolerance);
  }
  if(0) {
    cerr<<"SIZE: "<<result.size()<<endl;
    cerr<<"EACH: "<<endl;
    for(uint i=0;i<result.size();i++) cerr<<result.at(i)<<endl;
  }
  if(result.size()>1)  return result;
  result.clear();result.push_back(0.001);result.push_back(0.0025);result.push_back(0.005);result.push_back(0.01);result.push_back(0.025); return result;
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------------------- POINT GROUP OPERATIONS
// SYM::CalculatePointGroupKLattice
//
// This function calculate the whole point group and saves in the aflow.pgroup file
// Written by SC, dec 09
//
//ME20200114 - made function name capitalization more consistent with other functions
namespace SYM {
  bool CalculatePointGroupKLattice(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss, string format) {  // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    bool Krun=TRUE;
    //DX+CO START
    string pgroup_type="pgroupk";
    a.pgroupk.clear();  //CO20190204
    xstructure aa;aa.ReScale(1.0);         //contains FixLattices()
    aa.lattice=a.klattice;aa.ReScale(1.0); //contains FixLattices()
    aa.sym_eps = a.sym_eps; //DX NEED TO PROPOGATE sym_eps
    //a.FixLattices();  //don't change a
    //aa.FixLattices(); //done in ReScale()
    //aa.scale=1.0; //not sure why it's just set and not ReScaled
    //DX+CO END
    aa.title="KLATTICE";
    _atom atom;
    aa.AddAtom(atom,false); // just something in the origin;  //CO20230319 - add by type
    //  cerr << aa << endl;
    //DX20170808 - New klattice routine [OBSOLETE] Krun=Krun && SYM::CalculatePointGroup(FileMESSAGE,aa,aflags,FALSE,osswrite,oss);
    Krun=Krun && SYM::TransformSymmetryFromRealToReciprocal(FileMESSAGE,a,aa,aflags,osswrite,oss,pgroup_type); //DX20170808 - New klattice routine
    // ------------------------------------------------------------------------------
    a.pgroupk_calculated=TRUE;
    // ------------------------------------------------------------------------------
    _sym_op symop;
    for(uint k=0;k<aa.pgroup.size();k++) {
      symop=aa.pgroup.at(k);
      // cerr << symop << endl;
      //      SYM::AddSymmetryToStructure(a,symop.Uc,symop.Uf,symop.ctau,symop.ftau,symop.ctrasl,symop.ftrasl,symop.basis_atoms_map,symop.basis_types_map,_PGROUPK_);
      SYM::AddSymmetryToStructure(a,symop.Uc,symop.Uf,symop.ctau,symop.ftau,symop.ctrasl,symop.ftrasl,symop.basis_atoms_map,symop.basis_types_map,false,_PGROUPK_);
    }
    if(_write_) Krun=Krun && KBIN_SymmetryWrite(FileMESSAGE,a,aflags,_PGROUPK_,osswrite,oss,format);
    if(LDEBUG) cerr << "DEBUG: SYM::CalculatePointGroupKLattice: a.pgroupk.size()=" << a.pgroupk.size() << endl;
    return Krun;
  }
} // namespace SYM

// SYM::CalculatePointGroupKCrystal
namespace SYM {
  bool CalculatePointGroupKCrystal(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss, string format) {  // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    bool Krun=TRUE;
    //DX+CO START
    string pgroup_type="pgroupk_xtal";
    a.pgroupk_xtal.clear(); //CO20190204
    xstructure aa;aa.ReScale(1.0);         //contains FixLattices()
    aa.lattice=a.klattice;aa.ReScale(1.0); //contains FixLattices()
    aa.sym_eps = a.sym_eps; //DX NEED TO PROPOGATE sym_eps
    //a.FixLattices();  //don't change a
    //aa.FixLattices(); //done in ReScale()
    //aa.scale=1.0; //not sure why it's just set and not ReScaled
    //DX+CO END
    aa.title="KCRYSTAL";
    _atom atom;
    aa.AddAtom(atom,false); // just something in the origin;  //CO20230319 - add by species
    //  cerr << aa << endl;
    //DX20170808 - New klattice routine [OBSOLETE] Krun=Krun && SYM::CalculatePointGroup(FileMESSAGE,aa,aflags,FALSE,osswrite,oss);
    if(a.pgroup_xtal_calculated==FALSE) Krun=Krun && SYM::CalculatePointGroupCrystal(FileMESSAGE,a,aflags,_write_,osswrite,oss); // NEED POINT GROUP CRYSTAL //HE20220912
    Krun=Krun && SYM::TransformSymmetryFromRealToReciprocal(FileMESSAGE,a,aa,aflags,osswrite,oss,pgroup_type); //DX20170808 - New klattice routine
    // ------------------------------------------------------------------------------
    a.pgroupk_xtal_calculated=TRUE;
    // ------------------------------------------------------------------------------
    _sym_op symop;
    for(uint k=0;k<aa.pgroup_xtal.size();k++) {
      symop=aa.pgroup_xtal.at(k);
      // cerr << symop << endl;
      //      SYM::AddSymmetryToStructure(a,symop.Uc,symop.Uf,symop.ctau,symop.ftau,symop.ctrasl,symop.ftrasl,symop.basis_atoms_map,symop.basis_types_map,_PGROUPK_);
      SYM::AddSymmetryToStructure(a,symop.Uc,symop.Uf,symop.ctau,symop.ftau,symop.ctrasl,symop.ftrasl,symop.basis_atoms_map,symop.basis_types_map,false,_PGROUPK_XTAL_);
    }
    if(_write_) Krun=Krun && KBIN_SymmetryWrite(FileMESSAGE,a,aflags,_PGROUPK_XTAL_,osswrite,oss,format);
    if(LDEBUG) cerr << "DEBUG: SYM::CalculatePointGroupKCrystal: a.pgroupk_xtal.size()=" << a.pgroupk_xtal.size() << endl;
    return Krun;
  }
} // namespace SYM

// ----------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------- KLATTICE POINT GROUP OPERATIONS
// SYM::TransformSymmetryFromRealToReciprocal
//
// Uses the fact the reciprocal space is a change of basis from the real space.
// Transforms the real space symmetry elements into the reciprocal space symmetry elements.
// Otherwise, need to calculate explicitly, but what would the tolerance be in reciprocal space?
// The tolerance would be an ellipsoid, not a uniform sphere (i.e. direction dependent).
// This procedure circumvents this problem and ensures the number of real space point group operations 
// is always equal to the reciprocal space point group operations.
//DX20170808
//
namespace SYM {
  bool TransformSymmetryFromRealToReciprocal(ofstream &FileMESSAGE, xstructure& real_xstr, xstructure& reciprocal_xstr,
      _aflags& aflags, const bool& osswrite, ostream& oss, string& pgroup_type){
    ostringstream aus;
    xvector<double> ctau(3),ftau(3),ctrasl(3),ftrasl(3);                    // translation
    clear(ctau);clear(ftau);clear(ctrasl);clear(ftrasl);
    // ---------------------------------------------------------------------------
    // create placeholder for atoms_map for AddSymmetryToStructure() (no atom mappings for pgroups)
    std::vector<int> basis_atoms_map(reciprocal_xstr.atoms.size());        // will map each on each
    for(uint i=0;i<reciprocal_xstr.atoms.size();i++) basis_atoms_map[i]=i; // identically map each over each
    std::vector<int> basis_types_map(reciprocal_xstr.atoms.size());        // will map each on each
    for(uint i=0;i<reciprocal_xstr.atoms.size();i++) basis_types_map[i]=reciprocal_xstr.atoms[i].type; // identically map each over each
    // clean up pgroupk                                                     // just initialize
    reciprocal_xstr.pgroup.clear();
    reciprocal_xstr.pgroup_calculated=FALSE;                               // just initialize

    //pgroupk version
    if(pgroup_type == "pgroupk" && real_xstr.pgroup_calculated){
      string message="PGROUP_KLATTICE";
      // clean up pgroupk                                                     // just initialize
      reciprocal_xstr.pgroup.clear();
      reciprocal_xstr.pgroup_calculated=FALSE;                               // just initialize
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
      uint kk=0;
      xmatrix<double> Uf, Uc;

      xmatrix<double> trasp_real_f2c = trasp(real_xstr.f2c);
      xmatrix<double> trasp_real_c2f = trasp(real_xstr.c2f);
      for(uint i=0;i<real_xstr.pgroup.size();i++){
        //Uf = trasp_real_f2c*real_xstr.pgroup[i].Uc*trasp_real_c2f;
        //Uc = reciprocal_xstr.f2c*Uf*reciprocal_xstr.c2f;
        Uf = trasp(inverse(real_xstr.pgroup[i].Uf)); //DX20170905 - Faster/more direct calculation of reciprocal Uf (see D.E. Sands Vectors and Tensors in Cryst.)
        Uc = real_xstr.pgroup[i].Uc; //DX20170905 - Faster/more direct calculation of reciprocal Uf (see D.E. Sands Vectors and Tensors in Cryst.)
        kk=SYM::AddSymmetryToStructure(reciprocal_xstr,Uc,Uf,ctau,ftau,ctrasl,ftrasl,basis_atoms_map,basis_types_map,false,_PGROUP_,FALSE);  // kk number pgroups
        aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " " << reciprocal_xstr.pgroup[kk-1].str_type
          << " theta=";
        if(reciprocal_xstr.pgroup[kk-1].angle<100) aus << " ";
        if(reciprocal_xstr.pgroup[kk-1].angle<10) aus << " ";
        aus << reciprocal_xstr.pgroup[kk-1].angle << " " << " r=(" << reciprocal_xstr.pgroup[kk-1].axis << ")";
        aus << "    HM=" <<  aurostd::PaddedPRE(reciprocal_xstr.pgroup[kk-1].str_Hermann_Mauguin,2," ");     // remember vectors start from 0
        aus << "    S=" <<  aurostd::PaddedPRE(reciprocal_xstr.pgroup[kk-1].str_Schoenflies,2," ");          // remember vectors start from 0
        aus << endl;  // remember vectors start from 0
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
      }
      for(uint k=0;k<reciprocal_xstr.pgroup.size();k++) {
        Uc=reciprocal_xstr.pgroup[k].Uc;
        Uf=reciprocal_xstr.pgroup[k].Uf;
        roundoff(reciprocal_xstr.pgroup[k].Uc,_EPS_roundoff_);
      }
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: unique point group operations " << reciprocal_xstr.pgroup.size() << endl;// Message(__AFLOW_FILE__,aflags) << endl;
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: END " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
      return true;
    }
    //pgroupk_xtal version
    else if(pgroup_type == "pgroupk_xtal" && real_xstr.pgroup_xtal_calculated){

      string message="PGROUPK_XTAL";

      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
      uint kk=0;
      xmatrix<double> Uf, Uc;

      xmatrix<double> trasp_real_f2c = trasp(real_xstr.f2c);
      xmatrix<double> trasp_real_c2f = trasp(real_xstr.c2f);
      for(uint i=0;i<real_xstr.pgroup_xtal.size();i++){
        //Uf = trasp_real_f2c*real_xstr.pgroup_xtal[i].Uc*trasp_real_c2f;
        //Uc = reciprocal_xstr.f2c*Uf*reciprocal_xstr.c2f;
        Uf = trasp(inverse(real_xstr.pgroup_xtal[i].Uf)); //DX20170905 - Faster/more direct calculation of reciprocal Uf (see D.E. Sands Vectors and Tensors in Cryst.)
        Uc = real_xstr.pgroup_xtal[i].Uc; //DX20170905 - Faster/more direct calculation of reciprocal Uf (see D.E. Sands Vectors and Tensors in Cryst.)
        kk=SYM::AddSymmetryToStructure(reciprocal_xstr,Uc,Uf,ctau,ftau,ctrasl,ftrasl,basis_atoms_map,basis_types_map,false,_PGROUP_XTAL_,FALSE);  // kk number pgroups
        aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " " << reciprocal_xstr.pgroup_xtal[kk-1].str_type
          << " theta=";
        if(reciprocal_xstr.pgroup_xtal[kk-1].angle<100) aus << " ";
        if(reciprocal_xstr.pgroup_xtal[kk-1].angle<10) aus << " ";
        aus << reciprocal_xstr.pgroup_xtal[kk-1].angle << " " << " r=(" << reciprocal_xstr.pgroup_xtal[kk-1].axis << ")";
        aus << "    HM=" <<  aurostd::PaddedPRE(reciprocal_xstr.pgroup_xtal[kk-1].str_Hermann_Mauguin,2," ");     // remember vectors start from 0
        aus << "    S=" <<  aurostd::PaddedPRE(reciprocal_xstr.pgroup_xtal[kk-1].str_Schoenflies,2," ");          // remember vectors start from 0
        aus << endl;  // remember vectors start from 0
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
      }
      for(uint k=0;k<reciprocal_xstr.pgroup_xtal.size();k++) {
        Uc=reciprocal_xstr.pgroup_xtal[k].Uc;
        Uf=reciprocal_xstr.pgroup_xtal[k].Uf;
        roundoff(reciprocal_xstr.pgroup_xtal[k].Uc,_EPS_roundoff_);
      }
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: unique point group operations " << reciprocal_xstr.pgroup_xtal.size() << endl;// Message(__AFLOW_FILE__,aflags) << endl;
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << message << " Symmetry: END " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
      return true;
    }
    return false;
  }
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------- FACTOR GROUP OPERATIONS
// SYM::CalculateFactorGroup
//
// This function calculate the whole factor group and saves in the aflow.fgroup file
// Written by SC, feb 2010
//

//DX+CO START

//DX20190905 [OBSOLETE] double BringInCell(const double& x) {
//DX20190905 [OBSOLETE]   return BringInCell_20161115(x);
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] double BringInCell_20161115(const double& x) {
//DX20190905 [OBSOLETE]   return SYM::mod_one(x);
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] #ifndef COMPILE_SLIM
//DX20190905 [OBSOLETE] double BringInCell_20160101(const double& x) {
//DX20190905 [OBSOLETE]   //  if(x>0.0) { y=x-(double)floor(x); return y; }
//DX20190905 [OBSOLETE]   // else { if(x<0.0) y=x+floor(-x)+1.0;}
//DX20190905 [OBSOLETE]   if(x> _EPS_) return x-floor(x);
//DX20190905 [OBSOLETE]   if(x<-_EPS_) return x+floor(-x)+1.0;
//DX20190905 [OBSOLETE]   //  if(abs(x)<0.0001) y=0.0;
//DX20190905 [OBSOLETE]   // if(abs(x-1.0)<0.0001) y=0.0;
//DX20190905 [OBSOLETE]   return 0.0;
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] #endif
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] #ifndef COMPILE_SLIM
//DX20190905 [OBSOLETE] double BringInCell_20160101(const double& x,double tolerance) {
//DX20190905 [OBSOLETE]   //  if(x>0.0) { y=x-(double)floor(x); return y; }
//DX20190905 [OBSOLETE]   // else { if(x<0.0) y=x+floor(-x)+1.0;}
//DX20190905 [OBSOLETE]   if(x>tolerance) {return x-floor(x);}
//DX20190905 [OBSOLETE]   if(x<tolerance) {return x+floor(-x)+1.0;}
//DX20190905 [OBSOLETE]   //  if(abs(x)<0.0001) y=0.0;
//DX20190905 [OBSOLETE]   // if(abs(x-1.0)<0.0001) y=0.0;
//DX20190905 [OBSOLETE]   return 0.0;
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] #endif
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] xvector<double> BringInCell2_20161115(const xvector<double>& v_in) {
//DX20190905 [OBSOLETE]   return BringInCell2(v_in);
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] #ifndef COMPILE_SLIM
//DX20190905 [OBSOLETE] xvector<double> BringInCell2_20160101(const xvector<double>& v_in,double tolerance) {
//DX20190905 [OBSOLETE]   xvector<double> v_out(v_in.urows,v_in.lrows);
//DX20190905 [OBSOLETE]   for(int i=v_out.lrows;i<=v_out.urows;i++) {
//DX20190905 [OBSOLETE]     //v_out(i)=BringInCell(v_in(i),tolerance);
//DX20190905 [OBSOLETE]     v_out(i)=BringInCell_20160101(v_in(i),tolerance);
//DX20190905 [OBSOLETE]     if(abs(v_out(i))<tolerance) {v_out(i)=0.0;}
//DX20190905 [OBSOLETE]     if(abs(v_out(i)-1.0)<tolerance) {v_out(i)=0.0;}
//DX20190905 [OBSOLETE]   }
//DX20190905 [OBSOLETE]   return v_out;
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] #endif
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] xvector<double> BringInCell2(const xvector<double>& v_in) {
//DX20190905 [OBSOLETE]   return SYM::mod_one_xvec(v_in);
//DX20190905 [OBSOLETE] }

namespace SYM {
  //xstructure and _sym_op
  bool getFullSymBasis(const xstructure& a, _sym_op& symOp,bool map_types,vector<int>& basis_atoms_map,vector<int>& basis_types_map){
    double tolerance;
    if(a.sym_eps!=AUROSTD_NAN){
      tolerance=a.sym_eps;
    } else {
      tolerance=defaultTolerance(a);
    }
    return getFullSymBasis(a,symOp,map_types,tolerance,basis_atoms_map,basis_types_map);
  }
  bool getFullSymBasis(const xstructure& a, _sym_op& symOp,bool map_types,double tolerance,vector<int>& basis_atoms_map,vector<int>& basis_types_map){ //CO20190520 - removed pointers for bools and doubles, added const where possible
    double min_dist=a.dist_nn_min; //CO20180409
    if(min_dist == AUROSTD_NAN){min_dist=minimumDistance(a);} //CO20180409
    bool skew = isLatticeSkewed(a.lattice,min_dist,tolerance); //CO20180409
    return getFullSymBasis(a,symOp,map_types,skew,tolerance,basis_atoms_map,basis_types_map);
  }
  bool getFullSymBasis(const xstructure& a, _sym_op& symOp,bool map_types,bool skew,double tolerance,vector<int>& basis_atoms_map,vector<int>& basis_types_map){ //CO20190520 - removed pointers for bools and doubles, added const where possible
    return getFullSymBasis(a.atoms,a.lattice,a.c2f,a.f2c,symOp,map_types,skew,tolerance,basis_atoms_map,basis_types_map);
  }
  //atoms, c2f, f2c and _sym_op
  bool getFullSymBasis(const deque<_atom>& atoms,const xmatrix<double>& lattice,const xmatrix<double>& c2f, const xmatrix<double>& f2c, _sym_op& symOp,bool map_types,bool skew,double tolerance,vector<int>& basis_atoms_map,vector<int>& basis_types_map){  //MAIN FUNCTION //CO20190520 - removed pointers for bools and doubles, added const where possible
    bool LDEBUG=(FALSE || XHOST.DEBUG);

    bool fast=true;  //DX COME BACK HERE, MAKE THIS AN INPUT //CO, no need, if fast doesn't work, we have bigger problems

    //CO20180420 - fixing for ME
    bool basis_map_calculated_orig=symOp.basis_map_calculated;
    symOp.basis_map_calculated=false;
    symOp.basis_map_calculated=basis_map_calculated_orig;  //CO20140420 - make sure to include me before any return/throw!!!!!

    //bool cont = true;
    uint count=0;
    deque<_atom> transformedcrystal;
    deque<uint> index_to_check;
    _atom tmp;

    for(uint k=0;k<atoms.size();k++){index_to_check.push_back(k);}
    basis_atoms_map.clear();
    basis_types_map.clear();

    // ===== Check if applying the symmetry element along with internal translation maps to another atom ===== //
    for(uint k=0;k<atoms.size();k++){
      //CO IF PROBLEMS, turn off incell, and modify fpos manually (commented out line)
      //DX tmp=ApplyAtom(atoms[k],symOp,lattice,c2f,f2c,skew,true,false,ignoreFractionalOperation,tolerance); //incell, but no roff
      if(!ApplyAtomValidate(atoms[k],tmp,symOp,lattice,c2f,f2c,skew,true,false,tolerance)){ //incell, but no roff
        // Applying c2f and f2c to the Cartesian and fractional positions give different results (a tolerance issue) 
        symOp.basis_map_calculated=basis_map_calculated_orig;  //CO20140420 - make sure to include me before any return/throw!!!!!
        return FALSE;
      }
      //DX TEST
      //tmp.fpos=SYM::mod_one_xvec(tmp.fpos);
      // === If identity operator and T=transp(0,0,0), point should not change; don't use mod_one_xvec (tolerance can mess this up) === //
      //if(symOp.str_Hermann_Mauguin == "1" && abs(symOp.ftau[1])<=_ZERO_TOL_ && abs(symOp.ftau[2])<=_ZERO_TOL_ && abs(symOp.ftau[3])<=_ZERO_TOL_ ){
      //  tmp.fpos=atoms[k].fpos; //(Uf*atoms[k].fpos + ftau); //CO
      //  }
      //// === Else bring back in the cell === //
      //else {
      //tmp.fpos = (SYM::mod_one_xvec((symOp.Uf*atoms[k].fpos + symOp.ftau))); 
      //}
      //tmp.cpos = f2c*tmp.fpos; 
      //tmp.name = atoms[k].name;
      //tmp.type = atoms[k].type;
      //DX TEST
      uint mapped_index=0;
      if(MapAtomWithBasis(atoms, tmp, map_types, index_to_check, lattice, f2c, skew, tolerance, mapped_index,true)){ //DX20190619 - lattice and f2c as input
        basis_atoms_map.push_back(index_to_check[mapped_index]);
        basis_types_map.push_back(atoms[index_to_check[mapped_index]].type);
        if(fast){ //If you use fast, you should have verified that no atoms overlap (which is most of the time, hopefully)
          //DX20170731 - fast also inherently checks one-to-one
          index_to_check.erase(index_to_check.begin()+mapped_index);
        }
        count++;
      }
      else {
        /*if(LDEBUG) {
          cout << "FAILED 1 " << endl;
          cout << "HM = " << str_Hermann_Mauguin << endl;
          cout << "rotation " << Uf << endl;
          cout << "shift " << ftau << endl;
          cout << "tol = " << tolerance << endl;
          cout << "atom in consideration " << atoms[k].fpos << " -> (transformed) " << tmp.fpos << endl;
          for(uint i=0;i<atoms.size();i++){
          xvector<double> fdiff = atoms.at(i).fpos - tmp.fpos;
          if(skew){
          minimizeCartesianDistance(atoms[i].cpos,tmp.cpos,fdiff,c2f,f2c,tolerance); 
          }
          else {
          PBC(fdiff);
          }
          cerr << "normal atom " << atoms.at(i).fpos << " | " << tmp.fpos << ", diff = " << fdiff << " (mod=" << aurostd::modulus(f2c*fdiff) << ")" << endl;
          }
          }*/
        symOp.basis_map_calculated=basis_map_calculated_orig;  //CO20140420 - make sure to include me before any return/throw!!!!!
        return FALSE;
      }
    }
    if(count && (count == atoms.size())){
      /*cerr << symOp << endl;
        for(uint k=0;k<atoms.size();k++){
        cerr << "atom :" << atoms[k].fpos  << endl;
        cerr << "ratom:" << transformedcrystal[k].fpos  << endl;
        }*/
      //DX20170731 - Check one-to-one speed increase : START
      if(!fast){ //DX20170731 - If not fast, need to check if one-to-one
        vector<int> mappings = basis_atoms_map;
        for(uint i=0;i<mappings.size();i++){
          for(uint j=i+1;j<mappings.size();j++){
            if(mappings[i] == mappings[j]){
              symOp.basis_map_calculated=basis_map_calculated_orig;  //CO20140420 - make sure to include me before any return/throw!!!!!
              return FALSE;
            }
          }
          mappings.erase(mappings.begin()+i);
          i--;
        } 
      }
      //DX20170731 - Check one-to-one speed increase : END
      symOp.basis_map_calculated=basis_map_calculated_orig;  //CO20140420 - make sure to include me before any return/throw!!!!!
      return TRUE;
    }
    if(LDEBUG) {
      cerr << "SYM::getFullSymBasis: FAILED 2 - total mappings = " << count << ", should be " << atoms.size() << endl;
    }
    symOp.basis_map_calculated=basis_map_calculated_orig;  //CO20140420 - make sure to include me before any return/throw!!!!!
    return FALSE;
  }
  //DX+CO END

  //DX+CO START
  // warning: this variant gets wrong basis_atoms_map, it flips indices between transformed crystal and original (but everything else should be the same)
  bool getFullSymBasis_20170729(const deque<_atom>& atoms,const xmatrix<double>& lattice,const xmatrix<double>& c2f, const xmatrix<double>& f2c, _sym_op& symOp,bool map_types,bool skew,double tolerance,vector<int>& basis_atoms_map,vector<int>& basis_types_map){  //MAIN FUNCTION //CO20190520 - removed pointers for bools and doubles, added const where possible
    bool LDEBUG=(FALSE || XHOST.DEBUG);

    bool fast=true;  //DX COME BACK HERE, MAKE THIS AN INPUT //CO, no need, if fast doesn't work, we have bigger problems

    //bool cont = true;
    uint count=0;
    deque<_atom> transformedcrystal;
    deque<uint> index_to_check;
    _atom tmp;

    // ===== Check if applying the symmetry element along with internal translation maps to another atom ===== //
    for(uint k=0;k<atoms.size();k++){
      //CO IF PROBLEMS, turn off incell, and modify fpos manually (commented out line)
      //DX tmp=ApplyAtom(atoms[k],symOp,lattice,c2f,f2c,skew,true,false,ignoreFractionalOperation,tolerance); //incell, but no roff
      if(!ApplyAtomValidate(atoms[k],tmp,symOp,lattice,c2f,f2c,skew,true,false,tolerance)){ //incell, but no roff
        // Applying c2f and f2c to the Cartesian and fractional positions give different results (a tolerance issue) 
        return FALSE;
      }
      //DX TEST
      //tmp.fpos=SYM::mod_one_xvec(tmp.fpos);
      // === If identity operator and T=transp(0,0,0), point should not change; don't use mod_one_xvec (tolerance can mess this up) === //
      /*if(symOp.str_Hermann_Mauguin == "1" && abs(symOp.ftau[1])<=_ZERO_TOL_ && abs(symOp.ftau[2])<=_ZERO_TOL_ && abs(symOp.ftau[3])<=_ZERO_TOL_ ){
        tmp.fpos=atoms[k].fpos; //(Uf*atoms[k].fpos + ftau); //CO
        }
      // === Else bring back in the cell === //
      else {
      tmp.fpos = (SYM::mod_one_xvec((symOp.Uf*atoms[k].fpos + symOp.ftau))); 
      }
      tmp.cpos = f2c*tmp.fpos; 
      tmp.name = atoms[k].name;
      tmp.type = atoms[k].type;*/
      //DX TEST
      if(MapAtom(atoms,tmp,map_types,lattice,f2c,skew,tolerance)){ //DX20190619 - lattice and f2c as input
        transformedcrystal.push_back(tmp);
        index_to_check.push_back(k);
      }
      else {
        /*if(LDEBUG) {
          cout << "FAILED 1 " << endl;
          cout << "HM = " << str_Hermann_Mauguin << endl;
          cout << "rotation " << Uf << endl;
          cout << "shift " << ftau << endl;
          cout << "tol = " << tolerance << endl;
          cout << "atom in consideration " << atoms[k].fpos << " -> (transformed) " << tmp.fpos << endl;
          for(uint i=0;i<atoms.size();i++){
          xvector<double> fdiff = atoms.at(i).fpos - tmp.fpos;
          if(skew){
          minimizeCartesianDistance(atoms[i].cpos,tmp.cpos,fdiff,c2f,f2c,tolerance); 
          }
          else {
          PBC(fdiff);
          }
          cerr << "normal atom " << atoms.at(i).fpos << " | " << tmp.fpos << ", diff = " << fdiff << " (mod=" << aurostd::modulus(f2c*fdiff) << ")" << endl;
          }
          }*/
        return FALSE;
      }
    }
    basis_atoms_map.clear();
    basis_types_map.clear();
    // ===== If each atom is mapped to another atom (i.e. onto), we need to check if each mapping is unique (i.e. one-to-one) ===== //
    for(uint ix=0;ix<atoms.size();ix++){
      uint mapped_index=0;
      if(MapAtomWithBasis(transformedcrystal, atoms[ix], map_types,index_to_check, lattice, f2c, skew, tolerance, mapped_index,true)){ //DX20190619 - lattice and f2c as input, remove "Atom" prefix from name
        basis_atoms_map.push_back(index_to_check[mapped_index]);
        basis_types_map.push_back(atoms[index_to_check[mapped_index]].type);
        if(fast){ //If you use fast, you should have verified that no atoms overlap (which is most of the time, hopefully)
          index_to_check.erase(index_to_check.begin()+mapped_index);
        }
        count++;
      }
    }
    if(count && (count == atoms.size())){
      /*cerr << symOp << endl;
        for(uint k=0;k<atoms.size();k++){
        cerr << "atom :" << atoms[k].fpos  << endl;
        cerr << "ratom:" << transformedcrystal[k].fpos  << endl;
        }*/
      return TRUE;
    }
    if(LDEBUG) {
      cerr << "SYM::getFullSymBasis: FAILED 2 - total mappings = " << count << ", should be " << atoms.size() << endl;
    }
    return FALSE;
  }
}
//DX+CO END

//DX+CO START
namespace SYM {
  bool CalculateFactorGroup(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format) {
    double _eps_;//=_EPS_;
    if(a.sym_eps!=AUROSTD_NAN){ //Tolerance came from user or was calculated
      _eps_=a.sym_eps;
    }
    else {
      _eps_=defaultTolerance(a);
    }
    return SYM::CalculateFactorGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss,_eps_,format);
  }
} // namespace SYM

namespace SYM {
  bool CalculateFactorGroup(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss, double _eps_,string format) {
    return SYM::CalculateFactorGroup_20160801(FileMESSAGE,a,aflags,_write_,osswrite,oss,_eps_,format);
  }
} // namespace SYM

namespace SYM {
  bool CalculateFactorGroup_20160801(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;
    DEBUG_SYMMETRY=DEBUG_SYMMETRY || LDEBUG;    
    // Obtain the structure tolerance
    a.sym_eps=_eps_;
    bool skew = isLatticeSkewed(a.lattice,a.dist_nn_min,_eps_);

    if(DEBUG_SYMMETRY) cerr << __AFLOW_FUNC__ << " DEBUG" << endl;
    // ------------------------------------------------------------------------------
    ostringstream aus;
    bool Krun=TRUE,sym_found;
    _atom hatom,tatom;
    //xvector<double> fshift(3),cshift(3),temp(3),ctau(3),ftau(3),ctrasl(3),ftrasl(3);
    //xmatrix<double> Uf(3,3),Uc(3,3),identity(3,3);identity[1][1]=1.0;identity[2][2]=1.0;identity[3][3]=1.0;
    xmatrix<double> identity(3,3);identity[1][1]=1.0;identity[2][2]=1.0;identity[3][3]=1.0;

    std::vector<int> basis_atoms_map(a.atoms.size());
    std::vector<int> basis_types_map(a.atoms.size());

    xmatrix<double> f2c(3,3);f2c=trasp(a.lattice);
    xmatrix<double> c2f(3,3);c2f=inverse(trasp(a.lattice));

    if(a.pgroup_calculated==FALSE) Krun=Krun && SYM::CalculatePointGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss,format); // NEED POINT GROUP

    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    a.fgroup.clear();
    a.fgroup_calculated=FALSE;

    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: generating translations " << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " DEBUG breaking up by types" << endl;}
    deque<deque<_atom> > atoms_by_type = break_up_by_type(a.atoms);
    if(atoms_by_type.size()==0){
      cerr << __AFLOW_FUNC__ << " ERROR: Structure could not be broken up by atoms [dir=" << a.directory << "]" << endl;
      return false;
    }
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " DEBUG breaking up by types DONE" << endl;}
    uint smallest_group = atoms_by_type[0].size();
    uint index_for_smallest_group = 0;
    for(uint i=1;i<atoms_by_type.size();i++){
      if(atoms_by_type[i].size() < smallest_group){
        smallest_group = atoms_by_type[i].size();
        index_for_smallest_group = i;
      }
    }
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " DEBUG grabbed index for smallest group" << endl;}

    // ===== Loop over symmetry elements ===== //
    //xmatrix<double> R;
    //xvector<double> T;
    _atom tmp;
    //CO START
    _sym_op symOp; //CO
    symOp.is_fgroup=1;  //CO
    for(uint i=0;i<a.atoms.size();i++){
      symOp.basis_atoms_map.push_back(0);
      symOp.basis_types_map.push_back(0);
    }
    //CO END
    for(uint pg=0;pg<a.pgroup.size();pg++){
      //R = a.pgroup[pg].Uf;
      // ===== Use the smallest group of an atom type to find the possible translations ===== //
      //CO START
      symOp.Uf=a.pgroup[pg].Uf;
      symOp.Uc=a.pgroup[pg].Uc;
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " DEBUG starting pg=" << pg << endl;}
      //symOp.str_Hermann_Mauguin=a.pgroup[pg].str_Hermann_Mauguin; //no longer necessary
      //CO END
      for(uint j=0;j<atoms_by_type[index_for_smallest_group].size();j++){
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " DEBUG j=" << j << endl;}
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " DEBUG index_for_smallest_group=" << index_for_smallest_group << endl;}
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " DEBUG atoms_by_type.size()=" << atoms_by_type.size() << endl;}
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " DEBUG atoms_by_type[index_for_smallest_group].size()=" << atoms_by_type[index_for_smallest_group].size() << endl;}
        //DX20190905 [OBSOLETE-no more mod_one_xvec] symOp.ftau = mod_one_xvec(atoms_by_type[index_for_smallest_group][0].fpos - a.pgroup[pg].Uf*atoms_by_type[index_for_smallest_group][j].fpos);
        symOp.ftau = atoms_by_type[index_for_smallest_group][0].fpos - a.pgroup[pg].Uf*atoms_by_type[index_for_smallest_group][j].fpos; //DX20190905 - uses new bring in cell function
        BringInCellInPlace(symOp.ftau); //DX20190905 - uses new bring in cell function
        //CO START
        symOp.ctau=a.f2c*symOp.ftau;
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " DEBUG about to test symop" << endl;}
        //CO END
        //if(getFullSymBasis(a.atoms,a.pgroup[pg].Uf,a.c2f,a.f2c,a.pgroup[pg].str_Hermann_Mauguin,symOp.ftau,skew,_eps_,basis_atoms_map,basis_types_map))
        if(getFullSymBasis(a.atoms,a.lattice,a.c2f,a.f2c,symOp,TRUE,skew,_eps_,basis_atoms_map,basis_types_map))
        { //CO20200106 - patching for auto-indenting
          //Uc=a.pgroup[pg].Uc;Uf=a.pgroup[pg].Uf; ftau=symOp.ftau; // for safety try to check out inverse
          //DX20190905 [OBSOLETE-no more mod_one_xvec] symOp.ftau=mod_one_xvec(symOp.ftau);  // moves fractional coordinates to [0.0 to 1.0[
          BringInCellInPlace(symOp.ftau);  // moves fractional coordinates to [0.0 to 1.0[ //DX20190905 - uses new BringInCell function
          for(uint i=1;i<=3;i++) {
            if(symOp.ftau[i]>1.0-_ZERO_TOL_){ 
              symOp.ftau[i]=0.0; //DX roundoff
            }
          }
          symOp.ctau=f2c*symOp.ftau;clear(symOp.ctrasl);clear(symOp.ftrasl);
          sym_found=FALSE;
          for(uint ii=0;ii<a.fgroup.size()&&!sym_found;ii++){
            sym_found=(aurostd::identical(symOp.Uf,a.fgroup[ii].Uf) &&   //DX20171207 - Use xmatrix identical eps
                aurostd::identical(symOp.ctau,a.fgroup[ii].ctau,_eps_) && //DX20171207 - Use ctau not ftau
                aurostd::identical(basis_atoms_map,a.fgroup[ii].basis_atoms_map,0)); // look in all the list of operations
          }
          if(sym_found==FALSE) {                                          // new operation, generate and save it
            uint kk;
            kk=SYM::AddSymmetryToStructure(a,symOp.Uc,symOp.Uf,symOp.ctau,symOp.ftau,symOp.ctrasl,symOp.ftrasl,basis_atoms_map,basis_types_map,true,_FGROUP_,FALSE);        // in kk there is the number _0 point group
            aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP " << a.fgroup[kk-1].str_type;
            aus << " theta="; if(a.fgroup[kk-1].angle<100) aus << " "; if(a.fgroup[kk-1].angle<10) aus << " ";
            aus << a.fgroup[kk-1].angle << " " << " r=(" << a.fgroup[kk-1].axis << ")";
            // aus      << "  fshift=(" << fshift << ")  ";
            aus << "  ftau=(" << a.fgroup[kk-1].ftau << ")   ";
            aus << " basis_atoms_map = [ "; for(uint n=0;n<a.fgroup[kk-1].basis_atoms_map.size();n++) aus << a.fgroup[kk-1].basis_atoms_map[n] << " "; aus << "] ";
            aus << endl;  // remember vectors start from 0
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss); 
            basis_atoms_map.clear();
            basis_types_map.clear();
          }
        }
        else {
          basis_atoms_map.clear();
          basis_types_map.clear();
        }
      }
    }

    // ------------------------------------------------------------------------------
    a.fgroup_calculated=TRUE;
    // ------------------------------------------------------------------------------
    // check for non primitive
    if(1) {
      uint unities=0;
      for(uint ii=0;ii<a.fgroup.size();ii++)
        if(aurostd::identical(a.fgroup[ii].Uf,a.fgroup[0].Uf)){unities++;}    // check unity 0 is always unity //DX used to be _eps_/10.0  //DX20171207 - Use xmatrix identical eps
      if(unities>1) {
        aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: Cell not primitive: as big as " << unities << endl;// Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
      }
      for(uint ii=0+1;ii<a.fgroup.size();ii++) // avoid origin
        if(aurostd::identical(a.fgroup[ii].Uf,a.fgroup[0].Uf)) { //DX used to be _eps_/10.0  //DX20171207 - Use xmatrix identical eps
          aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: non-primitive internal translation ftau=(" << a.fgroup[ii].ftau << ") " << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
        }
    }
    // Printing and leaving
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: unique factor group operations " << a.fgroup.size() << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: END " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);

    if(_write_) Krun=Krun && KBIN_SymmetryWrite(FileMESSAGE,a,aflags,_FGROUP_,osswrite,oss,format);
    return Krun;
  }
} // namespace SYM
//DX+CO END

#ifndef COMPILE_SLIM
namespace SYM {
  bool CalculateFactorGroup_20160101(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_) {
    // AFLOW_FUNCTION_IMPLEMENTATION
    if(DEBUG_SYMMETRY) cerr << "DEBUG: SYM::CalculateFactorGroup" << endl;
    // ------------------------------------------------------------------------------
    ostringstream aus;
    bool Krun=TRUE,sym_found;
    uint num_suc_maps;
    _atom hatom,tatom;
    xvector<double> fshift(3),cshift(3),temp(3),ctau(3),ftau(3),ctrasl(3),ftrasl(3);
    xmatrix<double> Uf(3,3),Uc(3,3),identity(3,3);identity[1][1]=1.0;identity[2][2]=1.0;identity[3][3]=1.0;
    // double _eps_=_EPS_; optional

    std::vector<int> basis_atoms_map(a.atoms.size());
    std::vector<int> basis_types_map(a.atoms.size());
    // lattice fix to reduce the number of tested operations
    //JJPR BUG if(a.LatticeReduction_avoid==FALSE) a.LatticeReduction(); // better not repeat. but it is checked in
    //JJPR BUG if(a.Minkowski_calculated==FALSE) a.MinkowskiBasisReduction();
    a.FixLattices();a.BringInCell();
    // a.BringInWignerSeitz();

    xmatrix<double> f2c(3,3);f2c=trasp(a.lattice);
    xmatrix<double> c2f(3,3);c2f=inverse(trasp(a.lattice));

    if(a.pgroup_calculated==FALSE) Krun=Krun && SYM::CalculatePointGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss,_eps_); // NEED POINT GROUP
    _eps_=5*_eps_; // factor group builds up errors

    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    // clean up fgroup    // just initialize
    a.fgroup.clear();
    a.fgroup_calculated=FALSE;

    vector<xvector<double> > vfshifts;vfshifts.clear();

    // prepare all the possible shifts

    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: generating translations " << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);

    // 143  Cu2S3Si1_ICSD_24132 ICSD_24132    HEX
    // 144  Te1Zn1_ICSD_80076                 HEX    OK *************
    // 147  451  hP9                                 OK *************
    // 148  463  hR26                                OK *************
    // Bi8Te9_ICSD_106330
    //  int ijk=2;  // Y. La Page 1982 J. Appl. Cryst. 15, 255-259  need to cycle -2,-1,0,1,2

    // cerr << Getabc_angles(a.lattice,DEGREES) << endl;
    // int ijk=1;  // Y. La Page 1982 J. Appl. Cryst. 15, 255-259  need to cycle -2,-1,0,1,2
    // cerr << a << endl;

    xvector<double> corigin,forigin;
    //  corigin=a.atoms[0].cpos;
    // forigin=a.atoms[0].fpos;

    //  double radius=RadiusSphereLattice(a.lattice);  // smarter !
    fshift.clear();  vfshifts.push_back(fshift); // no translation is always good translation    //SC20120710...

    double ops[]={0.0,1./4,1./2,3./4,1./3,2./3,1./6,5./6,1.0,-1./4,-1./2,-3./4,-1./3,-2./3,-1./6,-5./6,-1.0};  // they are 15 0--14
    if(0)
      for(int ic=0;ic<17;ic++)
        for(int jc=0;jc<17;jc++)
          for(int kc=0;kc<17;kc++) {
            fshift[1]=ops[ic];fshift[2]=ops[jc];fshift[3]=ops[kc];
            if(0) vfshifts.push_back(fshift);
            if(1) {
              // fshift=BringInCell(fshift);
              for(uint ii=0;ii<vfshifts.size()&&!sym_found;ii++)
                sym_found=identical(fshift,vfshifts.at(ii),_eps_);  // look in all the list of operations
              if(sym_found==FALSE) vfshifts.push_back(fshift);
            }
          }

    if(0) {
      int ijk=2;
      for(int ic=-ijk;ic<=ijk;ic++)
        for(int jc=-ijk;jc<=ijk;jc++)
          for(int kc=-ijk;kc<=ijk;kc++)
            for(uint i=0;i<a.atoms.size();i++) {
              fshift[1]=(double)ic;fshift[2]=(double)jc;fshift[3]=(double)kc;
              vfshifts.push_back(fshift);
            }
    }

    if(1) {
      for(uint pg=0;pg<a.pgroup.size();pg++) {
        xstructure b(a);
        b=SYM::ApplyXstructure(a.pgroup[pg],a,FALSE);
        b.BringInCell(_eps_);
        // try to map a in b
        for(uint iat=0;iat<a.atoms.size();iat++)
          for(uint jat=0;jat<b.atoms.size();jat++) {
            if(a.atoms.size()<1000) {
              //  if(a.atoms[iat].type==b.atoms[jat].type  && iat==0) // shifts must work for all atoms, so it must work for atom0 too. TRY ALL OF FEW ATOMS Jan2012
              if(a.atoms[iat].type==b.atoms[jat].type) // shifts must work for all atoms, but picking them all helps to find more..
              { //CO20200106 - patching for auto-indenting
                fshift=a.atoms[iat].fpos-b.atoms[jat].fpos; sym_found=FALSE;
                fshift=BringInCell(fshift,_eps_/20.0);
                // fshift=BringInCell2(fshift);
                // fshift=BringInCell2(fshift,0.01);
                for(uint ii=0;ii<vfshifts.size()&&!sym_found;ii++) sym_found=identical(fshift,vfshifts.at(ii),_eps_/10.0);  // look in all the list of operations
                if(sym_found==FALSE) { vfshifts.push_back(fshift);}
              }
            } else {
              if(a.atoms[iat].type==b.atoms[jat].type  && iat==0) // shifts must work for all atoms, so it must work for atom0 too. TRY ONLY ONE Jan2012
                // if(a.atoms[iat].type==b.atoms[jat].type) // shifts must work for all atoms, but picking them all helps to find more..
              { //CO20200106 - patching for auto-indenting
                fshift=a.atoms[iat].fpos-b.atoms[jat].fpos; sym_found=FALSE;
                fshift=BringInCell(fshift,_eps_/20.0);
                // fshift=BringInCell2(fshift);
                for(uint ii=0;ii<vfshifts.size()&&!sym_found;ii++) sym_found=identical(fshift,vfshifts.at(ii),_eps_/10.0);  // look in all the list of operations
                if(sym_found==FALSE) { vfshifts.push_back(fshift);}
              }
            }
          }
      }
    }

    //  cerr << vfshifts.size() << endl;

    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: found test translations = " << vfshifts.size() << endl;
    for(uint ifshift=0;ifshift<vfshifts.size();ifshift++) { //SC
      //   aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: translations " << vfshifts.at(ifshift) << endl;
    }
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);

    // all symmetry operations are done within the fractional coordinate system
    // since translations back into the unit cell are straightforward
    // for each point group operation, apply it to the crystal and store the transformed
    // coordinates in atoms_rotated
    for(uint pg=0;pg<a.pgroup.size();pg++) {
      //   cerr << a.pgroup[pg].Uf << endl;
      std::deque<_atom> atoms_rotated;
      // make this so I do not have to reproduce all the time
      atoms_rotated.clear();
      for(uint i=0;i<a.atoms.size();i++) {
        hatom=SYM::ApplyAtom(a.atoms[i],a.pgroup[pg],a,FALSE); // store transformed atoms... for speed
        atoms_rotated.push_back(hatom);
      }
      for(uint ifshift=0;ifshift<vfshifts.size();ifshift++) { //SC
        //     cerr << vfshifts.at(ifshift) << endl;
        // consider all internal shifts that move an atom of the original structure to an atom of the transformed structure.
        // Then see if that translation maps the transformed crystal onto the original crystal.
        // create a basis_atoms_map and basis_types_map vector graph
        // operation is a.pgroup[pg]+vfshifts.at(ifshift)
        fshift=vfshifts.at(ifshift);     //      fshift=a.atoms[i].fpos-atoms_rotated[j].fpos; // try a translation
        // clean up basis mapping
        for(uint n=0;n<a.atoms.size();n++) { basis_atoms_map[n]=0; basis_types_map[n]=0; }
        //	if(LDEBUG) cerr << "basis_atoms_map.size()=" << basis_atoms_map.size() << endl;
        //	if(LDEBUG) cerr << "basis_types_map.size()=" << basis_types_map.size() << endl;

        num_suc_maps=0; // we store the number of successful maps with the same fshift
        for(uint n=0;n<a.atoms.size();n++) {
          for(uint m=0;m<atoms_rotated.size();m++) {
            if(a.atoms[n].type==atoms_rotated[m].type) { // we fshift them all
              temp=a.atoms[n].fpos-atoms_rotated[m].fpos-fshift;
              //	temp=BringInCell(temp,_eps_); // moves fractional coordinates to [0.0 to 1.0[
              temp=BringInCell2(temp); // moves fractional coordinates to [0.0 to 1.0[
              if(abs(temp[1])<_eps_ && abs(temp[2])<_eps_ && abs(temp[3])<_eps_) {
                num_suc_maps++;
                basis_atoms_map[atoms_rotated[m].basis]=a.atoms[n].basis;
                basis_types_map[atoms_rotated[m].basis]=a.atoms[n].type;
              }
            }
          }
        }    
        if(abs(fshift[1])<_eps_ && abs(fshift[2])<_eps_ && abs(fshift[3])<_eps_ && num_suc_maps==2*a.atoms.size()) num_suc_maps=a.atoms.size();  //SC20120710...
        // cerr << vfshifts.at(ifshift) << " " << num_suc_maps/2 << endl;
        // 	    if(0) if(modulus(BringInCell(temp,_eps_)-BringInCell2(temp))>0.1) {
        // 		cerr << "ERROR" << endl;
        // 		cerr << "temp =" << temp << endl;
        // 		cerr << "BringInCell(temp)=" << BringInCell(temp) << endl;
        // 		cerr << "BringInCell2(temp)=" << BringInCell2(temp) << endl;
        // 	    }

        if(num_suc_maps==a.atoms.size()) { // here if the shift maps all the atoms in the same atom,
          // check whether the symmetry operation already exists in the factorgroup array
          Uc=a.pgroup[pg].Uc;Uf=a.pgroup[pg].Uf; ftau=fshift; // for safety try to check out inverse                                  
          //	  if(iorder==1) {Uc=inverse(a.pgroup[pg].Uc);Uf=inverse(a.pgroup[pg].Uf);ftau=-inverse(a.pgroup[pg].Uf)*fshift;}
          ftau=BringInCell2(ftau);  // moves fractional coordinates to [0.0 to 1.0[
          for(uint i=1;i<=3;i++) if(ftau[i]>0.995) ftau[i]=0.0; // roundoff
          ctau=f2c*ftau;clear(ctrasl);clear(ftrasl);
          sym_found=FALSE;
          //	  cerr << "FOUND CANDIDATE" << endl;
          for(uint ii=0;ii<a.fgroup.size()&&!sym_found;ii++)
            // 2010-07-11	  //
            sym_found=(aurostd::identical(Uf,a.fgroup[ii].Uf,_eps_) && aurostd::identical(ftau,a.fgroup[ii].ftau,_eps_) && aurostd::identical(basis_atoms_map,a.fgroup[ii].basis_atoms_map,0)); // look in all the list of operations
          //	  sym_found=((aurostd::identical(Uf,a.fgroup[ii].Uf,_eps_) && aurostd::identical(ftau,a.fgroup[ii].ftau,_eps_)) || aurostd::identical(basis_atoms_map,a.fgroup[ii].basis_atoms_map,0)); // if basis identical then simple translation of bravais
          // if the symmetry operation is new, add it to the factorgroup array
          // and update all info about the sym_op object
          if(sym_found==FALSE) {                                          // new operation, generate and save it
            uint kk;
            kk=SYM::AddSymmetryToStructure(a,Uc,Uf,ctau,ftau,ctrasl,ftrasl,basis_atoms_map,basis_types_map,true,_FGROUP_);        // in kk there is the number _0 point group
            aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP " << a.fgroup[kk-1].str_type;
            aus << " theta="; if(a.fgroup[kk-1].angle<100) aus << " "; if(a.fgroup[kk-1].angle<10) aus << " ";
            aus << a.fgroup[kk-1].angle << " " << " r=(" << a.fgroup[kk-1].axis << ")";
            // aus	<< "  fshift=(" << fshift << ")  ";
            aus	<< "  ftau=(" << a.fgroup[kk-1].ftau << ")   ";
            aus << " basis_atoms_map = [ "; for(uint n=0;n<a.fgroup[kk-1].basis_atoms_map.size();n++) aus << a.fgroup[kk-1].basis_atoms_map[n] << " "; aus << "] ";
            //  aus << " basis_types_map = [ "; for(uint n=0;n<a.fgroup[kk-1].basis_types_map.size();n++) aus << a.fgroup[kk-1].basis_types_map[n] << " "; aus << "] ";
            //	  if(iorder==0) aus << " (found) ";
            // if(iorder==1) aus << " (inverted) ";
            aus << endl;  // remember vectors start from 0
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
            // for(n=0;n<a.atoms.size();n++) cout << basis_atoms_map[n] << " "; cerr << endl;
            // for(n=0;n<a.atoms.size();n++) cout << basis_types_map[n] << " "; cerr << endl;
          }
        }
      }
      atoms_rotated.clear();
    }

    // ------------------------------------------------------------------------------
    a.fgroup_calculated=TRUE;
    // ------------------------------------------------------------------------------
    // check for non primitive
    if(1) {
      uint unities=0;
      for(uint ii=0;ii<a.fgroup.size();ii++)
        if(aurostd::identical(a.fgroup[ii].Uf,a.fgroup[0].Uf,_eps_/10.0)) {unities++;}    // check unity 0 is always unity
      if(unities>1) {
        aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: Cell not primitive: as big as " << unities << endl;// Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
      }
      for(uint ii=0+1;ii<a.fgroup.size();ii++) // avoid origin
        if(aurostd::identical(a.fgroup[ii].Uf,a.fgroup[0].Uf,_eps_/10.0)) {
          aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: non-primitive internal translation ftau=(" << a.fgroup[ii].ftau << ") " << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
        }
    }
    // Printing and leaving
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: unique factor group operations " << a.fgroup.size() << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: END " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);

    if(_write_) Krun=Krun && KBIN_SymmetryWrite(FileMESSAGE,a,aflags,_FGROUP_,osswrite,oss);
    return Krun;
  }
} // namespace SYM

namespace SYM {
  bool CalculateFactorGroup_20100203(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss) {
    // AFLOW_FUNCTION_IMPLEMENTATION
    // ------------------------------------------------------------------------------
    ostringstream aus;
    bool Krun=TRUE,sym_found;
    uint num_suc_maps;
    _atom hatom,tatom;
    xvector<double> fshift(3),cshift(3),temp(3),ctau(3),ftau(3),ctrasl(3),ftrasl(3);
    xmatrix<double> Uf(3,3),Uc(3,3);
    double _eps_=_EPS_;
    std::vector<int> basis_atoms_map(a.atoms.size());
    std::vector<int> basis_types_map(a.atoms.size());
    // lattice fix to reduce the number of tested operations
    //JJPR BUG if(a.LatticeReduction_avoid==FALSE) a.LatticeReduction(); // better not repeat. but it is checked in
    //JJPR BUG a.MinkowskiBasisReduction();
    //JJPR BUG a.BringInWignerSeitz();

    a.FixLattices();a.BringInCell();

    xmatrix<double> f2c(3,3);f2c=trasp(a.lattice);
    xmatrix<double> c2f(3,3);c2f=inverse(trasp(a.lattice));

    if(a.pgroup_calculated==FALSE) Krun=Krun && SYM::CalculatePointGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss); // NEED POINT GROUP

    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    // clean up fgroup    // just initialize
    a.fgroup.clear();
    a.fgroup_calculated=FALSE;

    vector<xvector<double> > vfshifts;
    vfshifts.clear();
    // prepare all the possible shifts

    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: generating translations " << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);

    // 143  Cu2S3Si1_ICSD_24132 ICSD_24132    HEX
    // 144  Te1Zn1_ICSD_80076                 HEX    OK *************
    // 147  451  hP9                                 OK *************
    // 148  463  hR26                                OK *************
    // Bi8Te9_ICSD_106330
    //  int ijk=2;  // Y. La Page 1982 J. Appl. Cryst. 15, 255-259  need to cycle -2,-1,0,1,2

    cerr << a.lattice << endl;
    //  int ijk=1;  // Y. La Page 1982 J. Appl. Cryst. 15, 255-259  need to cycle -2,-1,0,1,2
    // cerr << a << endl;

    xvector<double> corigin,forigin;
    //  corigin=a.atoms[0].cpos;
    // forigin=a.atoms[0].fpos;

    double radius=RadiusSphereLattice(a.lattice);  // smarter !

    if(2) {
      //    for(int ic=-ijk;ic<=ijk;ic++)
      //       for(int jc=-ijk;jc<=ijk;jc++)
      // 	for(int kc=-ijk;kc<=ijk;kc++) {
      // 	  //	  fshift[1]=ic*a.lattice(1);fshift[2]=jc*a.lattice(2);fshift[3]=kc*a.lattice(3);
      // 	  vfshifts.push_back(fshift);
      // 	}

      for(uint pg=0;pg<a.pgroup.size();pg++) {
        for(uint j=0;j<a.atoms.size();j++) { //SC
          //	hatom=SYM::ApplyAtom(a.atoms[j],a.pgroup[pg],a,FALSE);
          hatom=a.atoms[j];
          hatom.cpos=a.pgroup[pg].Uc*(a.atoms[j].cpos-corigin)+corigin;
          hatom.fpos=a.pgroup[pg].Uf*(a.atoms[j].fpos-forigin)+forigin;
          for(uint i=0;i<a.atoms.size();i++)  //SC
            if(a.atoms[i].type==a.atoms[j].type) // && a.atoms[i].type==a.atoms[0].type) // I can choose the one with fewer atoms..
            { //CO20200106 - patching for auto-indenting
              cshift=a.atoms[i].cpos-hatom.cpos;

              if(modulus(cshift)<=radius+0.0001) {
                fshift=c2f*cshift;
                //   fshift=BringInCell(fshift);
                //	    if(abs(fshift[1])<0.001 || fshift[1]>0.999) fshift[1]=0.0;
                //  if(abs(fshift[2])<0.001 || fshift[2]>0.999) fshift[2]=0.0;
                //  if(abs(fshift[3])<0.001 || fshift[3]>0.999) fshift[3]=0.0;
                sym_found=FALSE;
                //  fshift=BringInCell(fshift/3.0);
                for(uint ii=0;ii<vfshifts.size()&&!sym_found;ii++)
                  sym_found=identical(fshift,vfshifts.at(ii),_eps_);  // look in all the list of operations
                if(sym_found==FALSE) vfshifts.push_back(fshift);
              }
            }
        }
      }
    }
    cerr << vfshifts.size() << endl;

    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: found test translations = " << vfshifts.size() << endl;
    for(uint ifshift=0;ifshift<vfshifts.size();ifshift++) { //SC
      // aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: translations " << vfshifts.at(ifshift) << endl;
    }
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);

    // all symmetry operations are done within the fractional coordinate system
    // since translations back into the unit cell are straightforward
    // for each point group operation, apply it to the crystal and store the transformed
    // coordinates in atoms_rotated
    for(uint pg=0;pg<a.pgroup.size();pg++) {
      std::deque<_atom> atoms_rotated;
      // make this so I do not have to reproduce all the time
      atoms_rotated.clear();
      for(uint i=0;i<a.atoms.size();i++) {
        //     hatom=SYM::ApplyAtom(a.atoms[i],a.pgroup[pg],a,FALSE); // store transformed atoms... for speed
        hatom=a.atoms[i];
        hatom.cpos=a.pgroup[pg].Uc*(a.atoms[i].cpos-corigin)+corigin;
        hatom.fpos=a.pgroup[pg].Uf*(a.atoms[i].fpos-forigin)+forigin;
        atoms_rotated.push_back(hatom);
      }
      for(uint ifshift=0;ifshift<vfshifts.size();ifshift++) { //SC
        // consider all internal shifts that move an atom of the original structure to an atom of the transformed structure.
        // Then see if that translation maps the transformed crystal onto the original crystal.
        // create a basis_atoms_map basis_types_map vector graph
        // operation is a.pgroup[pg]+vfshifts.at(ifshift)
        fshift=vfshifts.at(ifshift);     //      fshift=a.atoms[i].fpos-atoms_rotated[j].fpos; // try a translation
        // clean up basis mapping
        for(uint n=0;n<a.atoms.size();n++) basis_atoms_map[n]=0;
        for(uint n=0;n<a.atoms.size();n++) { basis_atoms_map[n]=0; basis_types_map[n]=0; }

        num_suc_maps=0; // we store the number of successful maps with the same fshift
        for(uint n=0;n<a.atoms.size();n++) {
          for(uint m=0;m<atoms_rotated.size();m++) {
            if(a.atoms[n].type==atoms_rotated[m].type) { // we fshift them all
              temp=a.atoms[n].fpos-atoms_rotated[m].fpos-fshift;
              temp=BringInCell(temp); // moves fractional coordinates to [0.0 to 1.0[
              if(modulus(temp)<1.001*sqrt(3)*_eps_) {
                num_suc_maps++;
                basis_atoms_map[atoms_rotated[m].basis]=a.atoms[n].basis;
                basis_types_map[atoms_rotated[m].basis]=a.atoms[n].type;
              }
            }
          }
        }    

        if(num_suc_maps==a.atoms.size()) { // here if the shift maps all the atoms in the same atom,
          uint iorder=0;
          for(iorder=0;iorder<1;iorder++)
          {
            // check whether the symmetry operation already exists in the factorgroup array
            if(iorder==0) {Uc=a.pgroup[pg].Uc;Uf=a.pgroup[pg].Uf; ftau=fshift;} // for safety try to check out inverse                                  
            //	  if(iorder==1) {Uc=inverse(a.pgroup[pg].Uc);Uf=inverse(a.pgroup[pg].Uf);ftau=-inverse(a.pgroup[pg].Uf)*fshift;}
            ftau=BringInCell(ftau);  // moves fractional coordinates to [0.0 to 1.0[
            for(uint i=1;i<=3;i++) if(ftau[i]>0.995) ftau[i]=0.0; // roundoff
            ctau=f2c*ftau;clear(ctrasl);clear(ftrasl);
            sym_found=FALSE;
            for(uint ii=0;ii<a.fgroup.size()&&!sym_found;ii++)
              sym_found=(identical(Uf,a.fgroup[ii].Uf,_eps_) && identical(ftau,a.fgroup[ii].ftau,_eps_));  // look in all the list of operations
            // if the symmetry operation is new, add it to the factorgroup array
            // and update all info about the sym_op object
            if(sym_found==FALSE) {                                          // new operation, generate and save it
              uint kk;
              kk=SYM::AddSymmetryToStructure(a,Uc,Uf,ctau,ftau,ctrasl,ftrasl,basis_atoms_map,basis_types_map,true,_FGROUP_);        // in kk there is the number _0 point group
              aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP " << a.fgroup[kk-1].str_type;
              aus << " theta="; if(a.fgroup[kk-1].angle<100) aus << " "; if(a.fgroup[kk-1].angle<10) aus << " ";
              aus << a.fgroup[kk-1].angle << " " << " r=(" << a.fgroup[kk-1].axis << ")";
              aus	<< "  fshift=(" << a.fgroup[kk-1].ftau << ")  ";
              for(uint n=0;n<a.fgroup[kk-1].basis_atoms_map.size();n++) aus << a.fgroup[kk-1].basis_atoms_map[n] << " ";
              aus << "  ";
              //  for(uint n=0;n<a.fgroup[kk-1].basis_types_map.size();n++) aus << a.fgroup[kk-1].basis_types_map[n] << " ";
              if(iorder==0) aus << " (found) ";
              if(iorder==1) aus << " (inverted) ";
              aus << endl;  // remember vectors start from 0
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
              // for(n=0;n<a.atoms.size();n++) cout << basis_atoms_map[n] << " "; cerr << endl;
              // for(n=0;n<a.atoms.size();n++) cout << basis_types_map[n] << " "; cerr << endl;
            }
          }
        }
      }
      atoms_rotated.clear();
    }

    //   if(0) {
    //     for(uint m=0;m<a.atoms.size();m++) { basis_atoms_map[m]=0;
    //       for(uint n=0;n<a.atoms.size();n++)
    // 	if(a.atoms[n].type==a.atoms[m].type) {
    // 	  temp=a.atoms[n].fpos-Uf*a.atoms[m].fpos-ftau;temp=BringInCell(temp);
    // 	  if(modulus(temp)<2*_eps_) basis_atoms_map[a.atoms[m].basis]=a.atoms[n].basis;
    // 	}
    //     }
    //   }

    //   // add next ones
    //   aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: composition " << endl;
    //   aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    //   for(uint prr=0;prr<4;prr++)
    //   for(uint fg1=1;fg1<a.fgroup.size();fg1++) {
    //     for(uint fg2=1;fg2<a.fgroup.size();fg2++) {
    //       // create y2=U2*(U1x+t1)+t2
    //       Uc=a.fgroup[fg2].Uc*a.fgroup[fg1].Uc;Uf=a.fgroup[fg2].Uf*a.fgroup[fg1].Uf;fshift=a.fgroup[fg2].Uf*a.fgroup[fg1].ftau+a.fgroup[fg2].ftau;

    //       std::deque<_atom> atoms_rotated;
    //     atoms_rotated.clear();
    //       for(uint i=0;i<a.atoms.size();i++) {
    // 	hatom=a.atoms[i];
    // 	hatom.cpos=Uc*a.atoms[i].cpos;hatom.fpos=Uf*a.atoms[i].fpos;
    // 	atoms_rotated.push_back(hatom);
    //       }

    //       {
    //  //    for(uint ifshift=0;ifshift<vfshifts.size();ifshift++) { //SC  //[CO20200106 - close bracket for indenting]}
    // //       fshift=vfshifts.at(ifshift);     //      fshift=a.atoms[i].fpos-atoms_rotated[j].fpos; // try a translation



    //      // check it
    //       for(uint n=0;n<a.atoms.size();n++) basis_atoms_map[n]=0;
    //       num_suc_maps=0; // we store the number of successful maps with the same fshift
    //       for(uint n=0;n<a.atoms.size();n++) {
    // 	for(uint m=0;m<atoms_rotated.size();m++) {
    // 	  if(a.atoms[n].type==atoms_rotated[m].type) { // we fshift them all
    // 	    temp=a.atoms[n].fpos-atoms_rotated[m].fpos-fshift;
    // 	    temp=BringInCell(temp); // moves fractional coordinates to [0.0 to 1.0[
    // 	    if(modulus(temp)<sqrt(3)*_eps_) {
    // 	      num_suc_maps++; basis_atoms_map[atoms_rotated[m].basis]=a.atoms[n].basis;
    // 	    }
    // 	  }
    // 	}
    //       }    

    //       if(num_suc_maps==a.atoms.size()) { // here if the shift maps all the atoms in the same atom,
    // 	ftau=BringInCell(ftau);  // moves fractional coordinates to [0.0 to 1.0[
    // 	for(uint i=1;i<=3;i++) if(ftau[i]>0.995) ftau[i]=0.0; // roundoff
    // 	ctau=f2c*ftau;clear(ctrasl);clear(ftrasl);
    // 	sym_found=FALSE;
    // 	for(uint ii=0;ii<a.fgroup.size()&&!sym_found;ii++)
    // 	  sym_found=(identical(Uf,a.fgroup[ii].Uf,_eps_) && identical(ftau,a.fgroup[ii].ftau,_eps_));  // look in all the list of operations
    // 	// if the symmetry operation is new, add it to the factorgroup array
    // 	// and update all info about the sym_op object
    // 	if(sym_found==FALSE) {                                          // new operation, generate and save it
    // 	  uint kk;
    // 	  kk=SYM::AddSymmetryToStructure(a,Uc,Uf,ctau,ftau,ctrasl,ftrasl,basis_atoms_map,basis_types_map,true,_FGROUP_);        // in kk there is the number _0 point group
    // 	  aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP " << a.fgroup[kk-1].str_type;
    // 	  aus << " theta="; if(a.fgroup[kk-1].angle<100) aus << " "; if(a.fgroup[kk-1].angle<10) aus << " ";
    // 	  aus << a.fgroup[kk-1].angle << " " << " r=(" << a.fgroup[kk-1].axis << ")";
    // 	  aus	<< "  fshift=(" << a.fgroup[kk-1].ftau << ")  ";
    // 	  for(uint n=0;n<a.fgroup[kk-1].basis_atoms_map.size();n++) aus << a.fgroup[kk-1].basis_atoms_map[n] << " ";
    // 	  aus << "  (" <<  fg1 << "," << fg2 << ")  ";
    // 	  aus << endl;  // remember vectors start from 0
    // 	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    // 	  // for(n=0;n<a.atoms.size();n++) cout << basis_atoms_map[n] << " "; cerr << endl;
    // 	}
    // 	}
    //       }
    //    atoms_rotated.clear();

    //     }
    //   }

    // ------------------------------------------------------------------------------
    a.fgroup_calculated=TRUE;
    // ------------------------------------------------------------------------------
    // Printing and leaving
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: unique factor group operations " << a.fgroup.size() << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: END " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    if(_write_) Krun=Krun && KBIN_SymmetryWrite(FileMESSAGE,a,aflags,_FGROUP_,osswrite,oss);
    return Krun;
  }
} // namespace SYM


namespace SYM {
  //bool CalculateFactorGroup_20100202(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss)
  bool _CalculateFactorGroup(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss)
  { //CO20200106 - patching for auto-indenting
    // ------------------------------------------------------------------------------
    // Some of this routine is taken from AVDV structure::calc_factor_group() code.
    //SC made modificatios for speed and consistency with aflow architecture
    ostringstream aus;
    bool Krun=TRUE,sym_found;
    uint num_suc_maps;
    _atom hatom;
    xvector<double> fshift(3),temp(3),ctau(3),ftau(3),ctrasl(3),ftrasl(3);
    xmatrix<double> Uf(3,3),Uc(3,3);
    double _eps_=_EPS_;
    std::vector<int> basis_atoms_map(a.atoms.size());
    std::vector<int> basis_types_map(a.atoms.size());
    // lattice fix to reduce the number of tested operations
    // (JJPR BUG) if(a.LatticeReduction_avoid==FALSE) a.LatticeReduction(); // better not repeat. but it is checked in
    a.FixLattices();a.BringInCell();

    xmatrix<double> f2c(3,3);f2c=trasp(a.lattice);
    xmatrix<double> c2f(3,3);c2f=inverse(trasp(a.lattice));
    if(a.pgroup_calculated==FALSE) Krun=Krun && SYM::CalculatePointGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss); // NEED POINT GROUP

    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    // clean up fgroup                                // just initialize
    a.fgroup.clear();
    a.fgroup_calculated=FALSE;

    int ijk=2;  // Y. La Page 1982 J. Appl. Cryst. 15, 255-259  need to cycle -2,-1,0,1,2
    vector<xvector<double> > vfshifts;
    for(int ic=-ijk;ic<=ijk;ic++)
      for(int jc=-ijk;jc<=ijk;jc++)
        for(int kc=-ijk;kc<=ijk;kc++)
          for(uint i=0;i<a.atoms.size();i++) {
            if(a.atoms[i].type==a.atoms[0].type) {  // a reference type (fgroup must work for all the atoms)
              //fshift=a.atoms[i].fpos-atomlist[j].fpos; // try a translation
              fshift=a.atoms[i].fpos;
              fshift[1]+=(double)ic;fshift[2]+=(double)jc;fshift[3]+=(double)kc;
              //	vfshifts.push_back(BringInCell(fshift));
              vfshifts.push_back(fshift);
            }
          }
    fshift.clear();

    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: found test translations = " << vfshifts.size() << endl;
    for(uint ifshift=0;ifshift<vfshifts.size();ifshift++) { //SC
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: translations " << vfshifts.at(ifshift) << endl;
    }
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);

    // all symmetry operations are done within the fractional coordinate system
    // since translations back into the unit cell are straightforward
    // for each point group operation, apply it to the crystal and store the transformed
    // coordinates in atomlist
    for(uint pg=0;pg<a.pgroup.size();pg++) {
      std::deque<_atom> atomlist;
      for(uint i=0;i<a.atoms.size();i++) {
        hatom=SYM::ApplyAtom(a.atoms[i],a.pgroup[pg],a,FALSE);
        atomlist.push_back(hatom);                  // store transformed str
      }
      // consider all internal shifts that move an atom of the original structure (e.g. the first
      // atom) to an atom of the transformed structure. Then see if that translation maps the
      // transformed crystal onto the original crystal.
      // create a basis_atoms_map basis_types_map  vector graph
      for(uint i=0;i<a.atoms.size();i++) { // try all the possible operations
        if(a.atoms[i].type==a.atoms[0].type) {  // fix type
          for(uint ifshift=0;ifshift<vfshifts.size();ifshift++) { //SC
            // clean up basis mapping
            for(uint n=0;n<a.atoms.size();n++) basis_atoms_map[n]=0;
            for(uint n=0;n<a.atoms.size();n++) basis_types_map[n]=0;
            fshift=vfshifts.at(ifshift)-atomlist[i].fpos; // try a translation
            // cerr << fshift << " type1=" << a.atoms[iref].type << " type2=" << atomlist[i].type << endl;
            num_suc_maps=0; // we store the number of successful maps with the same fshift
            for(uint n=0;n<a.atoms.size();n++) {
              for(uint m=0;m<atomlist.size();m++) {
                if(a.atoms[n].type==atomlist[m].type) {
                  // we fshift them all
                  temp=a.atoms[n].fpos-atomlist[m].fpos-fshift;
                  temp=BringInCell(temp); // moves fractional coordinates to [0.0 to 1.0[
                  if(abs(temp[1])<0.005 && abs(temp[2])<0.005 && abs(temp[3])<0.005) {
                    num_suc_maps++;
                    basis_atoms_map[atomlist[m].basis]=a.atoms[n].basis;
                    basis_types_map[atomlist[m].basis]=a.atoms[n].type;
                  }
                }
              }
            }    
            if(num_suc_maps==a.atoms.size()) { // here if the shift maps all the atoms in the same atom,
              cerr << "FOUND" << endl;
              // for(uint i=1;i<=2;i++)
              {
                // check whether the symmetry operation already exists in the factorgroup array
                //		if(i==1)
                { Uc=a.pgroup[pg].Uc;Uf=a.pgroup[pg].Uf; ftau=fshift;}
                //	if(i==2) { Uc=inverse(a.pgroup[pg].Uc);Uf=inverse(a.pgroup[pg].Uf); ftau=-inverse(a.pgroup[pg].Uf)*fshift;}
                ftau=BringInCell(ftau);  // moves fractional coordinates to [0.0 to 1.0[
                for(uint i=1;i<=3;i++) if(ftau[i]>0.995) ftau[i]=0.0; // roundoff
                ctau=f2c*ftau;
                clear(ctrasl);clear(ftrasl);
                sym_found=FALSE;
                for(uint ii=0;ii<a.fgroup.size()&&!sym_found;ii++)
                  sym_found=(identical(Uf,a.fgroup[ii].Uf,_eps_) && identical(ftau,a.fgroup[ii].ftau,_eps_));  // look in all the list of operations
                // if the symmetry operation is new, add it to the factorgroup array
                // and update all info about the sym_op object
                if(sym_found==FALSE) {                                          // new operation, generate and save it
                  uint kk;
                  kk=SYM::AddSymmetryToStructure(a,Uc,Uf,ctau,ftau,ctrasl,ftrasl,basis_atoms_map,basis_types_map,true,_FGROUP_);        // in kk there is the number _0 point group
                  aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP " << a.fgroup[kk-1].str_type
                    << " theta="; if(a.fgroup[kk-1].angle<100) aus << " "; if(a.fgroup[kk-1].angle<10) aus << " ";
                  aus << a.fgroup[kk-1].angle << " " << " r=(" << a.fgroup[kk-1].axis << ")"
                    << "  fshift=(" << fshift << ")"
                    << "  ftau=(" << a.fgroup[kk-1].ftau << ")"
                    << endl;  // remember vectors start from 0
                  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
                  // for(n=0;n<a.atoms.size();n++) cout << basis_atoms_map[n] << " "; cerr << endl;
                  // for(n=0;n<a.atoms.size();n++) cout << basis_types_map[n] << " "; cerr << endl;
                }
              }
            }
          }
        }
      }
      atomlist.clear();
    }
    // ------------------------------------------------------------------------------
    a.fgroup_calculated=TRUE;
    // ------------------------------------------------------------------------------
    // Printing and leaving
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: unique factor group operations " << a.fgroup.size() << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: END " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    if(_write_) Krun=Krun && KBIN_SymmetryWrite(FileMESSAGE,a,aflags,_FGROUP_,osswrite,oss);
    return Krun;
  }
} // namespace SYM


// ----------------------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------- FACTOR GROUP OPERATIONS
// SYM::CalculateFactorGroup_20201001
namespace SYM {
  //bool CalculateFactorGroup_20201001(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss) // AFLOW_FUNCTION_IM
  bool __CalculateFactorGroup(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss)
  { //CO20200106 - patching for auto-indenting
    // ------------------------------------------------------------------------------
    // Some of this routine is taken from AVDV structure::calc_factor_group() code.
    //SC made modificatios for speed and consistency with aflow architecture
    ostringstream aus;
    bool Krun=TRUE,sym_found;
    uint iref,num_suc_maps;
    _atom hatom;
    xvector<double> fshift(3),temp(3),ctau(3),ftau(3),ctrasl(3),ftrasl(3);
    xmatrix<double> Uf(3,3),Uc(3,3);
    double _eps_=_EPS_;
    std::vector<int> basis_atoms_map(a.atoms.size());
    std::vector<int> basis_types_map(a.atoms.size());
    // lattice fix to reduce the number of tested operations
    // (JJPR BUG) if(a.LatticeReduction_avoid==FALSE) a.LatticeReduction(); // better not repeat. but it is checked in
    a.FixLattices();a.BringInCell();

    xmatrix<double> f2c(3,3);f2c=trasp(a.lattice);
    xmatrix<double> c2f(3,3);c2f=inverse(trasp(a.lattice));
    if(a.pgroup_calculated==FALSE) Krun=Krun && SYM::CalculatePointGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss); // NEED POINT GROUP

    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    // clean up fgroup                                // just initialize
    a.fgroup.clear();
    a.fgroup_calculated=FALSE;

    int ijk=2;
    vector<xvector<double> > vfshifts;
    for(int ic=-ijk;ic<=ijk;ic++)
      for(int jc=-ijk;jc<=ijk;jc++)
        for(int kc=-ijk;kc<=ijk;kc++) {
          fshift[1]=(double)ic;fshift[2]=(double)jc;fshift[3]=(double)kc;
          // fshift=BringInCell(fshift);  // no shift here
          sym_found=FALSE;
          for(uint ii=0;ii<vfshifts.size()&&!sym_found;ii++)
            sym_found=identical(fshift,vfshifts.at(ii),_eps_); // look in all the list of operations
          if(sym_found==FALSE) vfshifts.push_back(fshift);
        }
    fshift.clear();

    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: found test translations = " << vfshifts.size() << endl;
    for(uint ifshift=0;ifshift<vfshifts.size();ifshift++) { //SC
      // aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: translations " << vfshifts.at(ifshift) << endl;
    }
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);

    // all symmetry operations are done within the fractional coordinate system
    // since translations back into the unit cell are straightforward
    // for each point group operation, apply it to the crystal and store the transformed
    // coordinates in atomlist
    for(uint pg=0;pg<a.pgroup.size();pg++) {
      std::deque<_atom> atomlist;
      for(uint i=0;i<a.atoms.size();i++) {
        hatom=SYM::ApplyAtom(a.atoms[i],a.pgroup[pg],a,FALSE);
        atomlist.push_back(hatom);                  // store transformed str
      }
      // consider all internal shifts that move an atom of the original structure (e.g. the first
      // atom) to an atom of the transformed structure. Then see if that translation maps the
      // transformed crystal onto the original crystal.
      // create a basis_atoms_map basis_types_map vector graph
      //   iref=0;  // anton
      for(uint ifshift=0;ifshift<vfshifts.size();ifshift++) { //SC
        for(iref=0;iref<a.atoms.size();iref++) { //SC
          // cout << iref << endl;cout.flush();
          for(uint i=0;i<a.atoms.size();i++) {
            if(a.atoms[iref].type==a.atoms[i].type) {
              // clean up basis mapping
              for(uint n=0;n<a.atoms.size();n++) basis_atoms_map[n]=0;
              for(uint n=0;n<a.atoms.size();n++) basis_types_map[n]=0;
              // creating translation shift
              fshift=a.atoms[iref].fpos-atomlist[i].fpos; // try a translation
              fshift=fshift+vfshifts.at(ifshift);
              // cerr << fshift << " type1=" << a.atoms[iref].type << " type2=" << atomlist[i].type << endl;
              num_suc_maps=0; // we store the number of successful maps with the same fshift
              for(uint n=0;n<a.atoms.size();n++) {
                for(uint m=0;m<atomlist.size();m++) {
                  if(a.atoms[n].type==atomlist[m].type) {
                    // we fshift them all
                    temp=a.atoms[n].fpos-atomlist[m].fpos-fshift;
                    temp=BringInCell(temp); // moves fractional coordinates to [0.0 to 1.0[
                    if(abs(temp[1])<0.005 && abs(temp[2])<0.005 && abs(temp[3])<0.005) {
                      num_suc_maps++;
                      basis_atoms_map[atomlist[m].basis]=a.atoms[n].basis;
                      basis_types_map[atomlist[m].basis]=a.atoms[n].type;
                    }
                  }
                }
              }    
              if(num_suc_maps==a.atoms.size()) { // here if the shift maps all the atoms in the same atom,
                // for(uint i=1;i<=2;i++)
                {
                  // check whether the symmetry operation already exists in the factorgroup array
                  //		if(i==1)
                  { Uc=a.pgroup[pg].Uc;Uf=a.pgroup[pg].Uf; ftau=fshift;}
                  //	if(i==2) { Uc=inverse(a.pgroup[pg].Uc);Uf=inverse(a.pgroup[pg].Uf); ftau=-inverse(a.pgroup[pg].Uf)*fshift;}
                  ftau=BringInCell(ftau);  // moves fractional coordinates to [0.0 to 1.0[
                  for(uint i=1;i<=3;i++) if(ftau[i]>0.995) ftau[i]=0.0; // roundoff
                  ctau=f2c*ftau;
                  clear(ctrasl);clear(ftrasl);
                  sym_found=FALSE;
                  for(uint ii=0;ii<a.fgroup.size()&&!sym_found;ii++)
                    sym_found=(identical(Uf,a.fgroup[ii].Uf,_eps_) && identical(ftau,a.fgroup[ii].ftau,_eps_));  // look in all the list of operations
                  // if the symmetry operation is new, add it to the factorgroup array
                  // and update all info about the sym_op object
                  if(sym_found==FALSE) {                                          // new operation, generate and save it
                    uint kk;
                    kk=SYM::AddSymmetryToStructure(a,Uc,Uf,ctau,ftau,ctrasl,ftrasl,basis_atoms_map,basis_types_map,true,_FGROUP_);        // in kk there is the number _0 point group
                    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP " << a.fgroup[kk-1].str_type
                      << " theta="; if(a.fgroup[kk-1].angle<100) aus << " "; if(a.fgroup[kk-1].angle<10) aus << " ";
                    aus << a.fgroup[kk-1].angle << " " << " r=(" << a.fgroup[kk-1].axis << ")"
                      << "  fshift=(" << a.fgroup[kk-1].ftau << ")"
                      << endl;  // remember vectors start from 0
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
                    // for(n=0;n<a.atoms.size();n++) cout << basis_atoms_map[n] << " "; cerr << endl;
                    // for(n=0;n<a.atoms.size();n++) cout << basis_types_map[n] << " "; cerr << endl;
                  }
                }
              }
            }
          }
        }
      }
      atomlist.clear();
    }
    // ------------------------------------------------------------------------------
    a.fgroup_calculated=TRUE;
    // ------------------------------------------------------------------------------
    // Printing and leaving
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: unique factor group operations " << a.fgroup.size() << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: END " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    if(_write_) Krun=Krun && KBIN_SymmetryWrite(FileMESSAGE,a,aflags,_FGROUP_,osswrite,oss);
    return Krun;
  }
} // namespace SYM


// ----------------------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------- FACTOR GROUP OPERATIONS
// SYM::CalculateFactorGroup_20200708
//
// This function calculate the whole factor group and saves in the aflow.fgroup file
// Written by SC, aug 2007
//
namespace SYM {
  bool CalculateFactorGroup_20200708(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss) {        // AFLOW_FUNCTION_IMPLEMENTATION
    // ------------------------------------------------------------------------------
    // Some of this routine is taken from AVDV structure::calc_factor_group() code.
    //SC made modificatios for speed and consistency with aflow architecture
    ostringstream aus;
    bool Krun=TRUE,sym_found;
    uint iref,num_suc_maps;
    _atom hatom;
    xvector<double> fshift(3),temp(3),ctau(3),ftau(3),ctrasl(3),ftrasl(3);
    xmatrix<double> Uf(3,3),Uc(3,3);
    double _eps_=_EPS_;
    std::vector<int> basis_atoms_map(a.atoms.size());
    std::vector<int> basis_types_map(a.atoms.size());
    // lattice fix to reduce the number of tested operations
    // (JJPR BUG) if(a.LatticeReduction_avoid==FALSE) a.LatticeReduction(); // better not repeat. but it is checked in
    a.FixLattices();a.BringInCell();

    xmatrix<double> f2c(3,3);f2c=trasp(a.lattice);
    xmatrix<double> c2f(3,3);c2f=inverse(trasp(a.lattice));
    if(a.pgroup_calculated==FALSE) Krun=Krun && SYM::CalculatePointGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss); // NEED POINT GROUP

    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    // clean up fgroup                                // just initialize
    a.fgroup.clear();
    a.fgroup_calculated=FALSE;

    // all symmetry operations are done within the fractional coordinate system
    // since translations back into the unit cell are straightforward
    // for each point group operation, apply it to the crystal and store the transformed
    // coordinates in atomlist
    for(uint pg=0; pg<a.pgroup.size(); pg++) {
      std::deque<_atom> atomlist;
      for(uint i=0;i<a.atoms.size();i++) {
        hatom=SYM::ApplyAtom(a.atoms[i],a.pgroup[pg],a,FALSE);
        atomlist.push_back(hatom);                  // store transformed str
      }
      // consider all internal shifts that move an atom of the original structure (e.g. the first
      // atom) to an atom of the transformed structure. Then see if that translation maps the
      // transformed crystal onto the original crystal.
      // create a basis_atoms_map basis_types_map vector graph
      //   for(iref=0;iref<a.atoms.size();iref++)
      {
        // cout << iref << endl;cout.flush();
        iref=0;
        for(uint i=0;i<a.atoms.size();i++) {
          if(a.atoms[iref].type==a.atoms[i].type) {
            // clean up basis mapping
            for(uint n=0;n<a.atoms.size();n++) basis_atoms_map[n]=0;
            for(uint n=0;n<a.atoms.size();n++) basis_types_map[n]=0;
            // creating translation shift
            fshift=a.atoms[iref].fpos-atomlist[i].fpos; // try a translation
            // cerr << fshift << " type1=" << a.atoms[iref].type << " type2=" << atomlist[i].type << endl;
            num_suc_maps=0; // we store the number of successful maps with the same fshift
            for(uint n=0;n<a.atoms.size();n++) {
              for(uint m=0;m<atomlist.size();m++) {
                if(a.atoms[n].type==atomlist[m].type) {
                  // we fshift them all
                  temp=a.atoms[n].fpos-atomlist[m].fpos-fshift;
                  temp=BringInCell(temp); // moves fractional coordinates to [0.0 to 1.0[
                  if(abs(temp[1])<0.0005 && abs(temp[2])<0.0005 && abs(temp[3])<0.0005) {
                    num_suc_maps++;
                    basis_atoms_map[atomlist[m].basis]=a.atoms[n].basis;
                    basis_types_map[atomlist[m].basis]=a.atoms[n].type;
                  }
                }
              }
            }	
            if(num_suc_maps==a.atoms.size()) {
              // here if the shift maps all the atoms in the same atom,
              fshift=BringInCell(fshift);  // moves fractional coordinates to [0.0 to 1.0[
              // check whether the symmetry operation already exists in the factorgroup array
              Uc=a.pgroup[pg].Uc;
              Uf=a.pgroup[pg].Uf;
              ftau=fshift;
              ctau=f2c*fshift;
              clear(ctrasl);
              clear(ftrasl);
              sym_found=FALSE;
              for(uint ii=0;ii<a.fgroup.size()&&!sym_found;ii++)
                sym_found=(identical(Uf,a.fgroup[ii].Uf,_eps_) &&
                    identical(ftau,a.fgroup[ii].ftau,_eps_));    // look in all the list of operations
              // if the symmetry operation is new, add it to the factorgroup array
              // and update all info about the sym_op object
              if(sym_found==FALSE)                                          // new operation, generate and save it
                //	  if(sym_found==FALSE || a.fgroup.size() == 0)
              { //CO20200106 - patching for auto-indenting
                uint kk;
                kk=SYM::AddSymmetryToStructure(a,Uc,Uf,ctau,ftau,ctrasl,ftrasl,basis_atoms_map,basis_types_map,true,_FGROUP_);        // in kk there is the number _0 point group
                aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP " << a.fgroup[kk-1].str_type
                  << " theta="; if(a.fgroup[kk-1].angle<100) aus << " "; if(a.fgroup[kk-1].angle<10) aus << " ";
                aus << a.fgroup[kk-1].angle << " " << " r=(" << a.fgroup[kk-1].axis << ")"
                  << "  ftau=(" << a.fgroup[kk-1].ftau << ")"
                  << endl;  // remember vectors start from 0
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
                // for(n=0;n<a.atoms.size();n++) cout << basis_atoms_map[n] << " "; cerr << endl;
                // for(n=0;n<a.atoms.size();n++) cout << basis_types_map[n] << " "; cerr << endl;
              }
            }

          }
        }
      }
      atomlist.clear();
    }
    // ------------------------------------------------------------------------------
    a.fgroup_calculated=TRUE;
    // ------------------------------------------------------------------------------
    // Printing and leaving
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: unique factor group operations " << a.fgroup.size() << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: END " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    if(_write_) Krun=Krun && KBIN_SymmetryWrite(FileMESSAGE,a,aflags,_FGROUP_,osswrite,oss);
    return Krun;
  }
} // namespace SYM
#endif

// ----------------------------------------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------------------- SPACE GROUP OPERATIONS
// SYM::CalculateSpaceGroup
//
// This function calculate the whole space group and saves in the aflow.fgroup file
// Written by SC, aug 2007
//
namespace SYM {
  //DX+CO START
  bool CalculateSpaceGroup(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format) {
    return CalculateSpaceGroup_20160801(FileMESSAGE,a,aflags,_write_,osswrite,oss,format);
  }
}
namespace SYM {
  bool CalculateSpaceGroup_20160801(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format) {        // AFLOW_FUNCTION_IMPLEMENTATION
    // ------------------------------------------------------------------------------
    // Some of this routine is taken from AVDV structure::calc_factor_group() code.
    //SC made modificatios for speed and consistency with aflow architecture
    ostringstream aus;
    bool Krun=TRUE,sym_found;
    _atom hatom;
    // lattice fix to reduce the number of tested operations
    //  a.LatticeReduction(); better not repeat
    // lattice fixed, now start
    //DX+CO START
    //DX a.FixLattices();a.BringInCell();
    //DX+CO END


    a.FixLattices();
    a.BringInCell();
    xvector<double> a1(3),a2(3),a3(3);                    // lattice vectors
    a1=a.lattice(1);a2=a.lattice(2);a3=a.lattice(3);      // a1,a2,a3 are the rows of the lattice matrix
    xvector<double> ctau(3),ftau(3),ctrasl(3),ftrasl(3);  // translations
    xmatrix<double> Uf(3,3),Uc(3,3);
    xvector<int> dims(3);
    std::vector<int> basis_atoms_map(a.atoms.size());       // will map like the fgroup
    std::vector<int> basis_types_map(a.atoms.size());       // will map like the fgroup
    xmatrix<double> f2c(3,3);f2c=trasp(a.lattice);
    xmatrix<double> c2f(3,3);c2f=inverse(trasp(a.lattice));
    if(a.pgroup_calculated==FALSE) Krun=Krun && SYM::CalculatePointGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss,format);  // NEED POINT GROUP
    if(a.fgroup_calculated==FALSE) Krun=Krun && SYM::CalculateFactorGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss,format); // NEED FACTOR GROUP

    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "SGROUP Symmetry: BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    // clean up sgroup                                // just initialize
    a.sgroup.clear();
    a.sgroup_calculated=FALSE;
    if(a.sgroup_radius<=0.1) {
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "SGROUP Symmetry: picking default normalized radius = " << _calculate_symmetry_default_sgroup_radius_<< "  " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
      a.sgroup_radius=_calculate_symmetry_default_sgroup_radius_*MaxStructureLattice(a);
    }
    dims=LatticeDimensionSphere(a.lattice,a.sgroup_radius);
    // for(i=1;i<=3;i++) if(dims(i)>3)dims(i)=3;
    a.sgroup_radius_dims=dims;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "SGROUP Symmetry: inside sphere (radius=" << a.sgroup_radius << ", dims = [" << dims[1] << "," << dims[2] << "," << dims[3] << "]) " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    // all symmetry operations are done within the fractional coordinate system
    // since translations back into the unit cell are straightforward
    // for each point group operation, apply it to the crystal and store the transformed
    // coordinates in atomlist
    //  dims.clear();
    for(int i=-dims[1];i<=dims[1];i++){
      for(int j=-dims[2];j<=dims[2];j++){
        for(int k=-dims[3];k<=dims[3];k++){
          ctrasl=((double)i)*a1+((double)j)*a2+((double)k)*a3;
          ftrasl(1)=(double)i;ftrasl(2)=(double) j;ftrasl(3)=(double) k;
          for(uint fg=0; fg<a.fgroup.size(); fg++){
            Uc=a.fgroup[fg].Uc;
            Uf=a.fgroup[fg].Uf;
            ctau=a.fgroup[fg].ctau;
            ftau=a.fgroup[fg].ftau;
            basis_atoms_map=a.fgroup[fg].basis_atoms_map;
            basis_types_map=a.fgroup[fg].basis_types_map;
            sym_found=FALSE;
            // no check required, as they are different by construction (not really)
            // 	  for(uint ii=0;ii<a.sgroup.size()&&!sym_found;ii++)
            // 	    //	    sym_found=(identical(Uf,a.sgroup[ii].Uf,_EPS_) && identical(ftau,a.sgroup[ii].ftau,_EPS_) && identical(ftrasl,a.sgroup[ii].ftrasl,_EPS_));    // look in all the list of operations
            // 	    sym_found=(identical(Uf,a.sgroup[ii].Uf,_EPS_) && identical(ftau+ftrasl,a.sgroup[ii].ftau+a.sgroup[ii].ftrasl,_EPS_));    // look in all the list of operations
            // if the symmetry operation is new, add it to the factorgroup array
            // and update all info about the sym_op object
            if(sym_found==FALSE)                                          // new operation, generate and save it
              //	  if(sym_found==FALSE || a.sgroup.size() == 0)
            { //CO20200106 - patching for auto-indenting
              uint kk;
              kk=SYM::AddSymmetryToStructure(a,Uc,Uf,ctau,ftau,ctrasl,ftrasl,basis_atoms_map,basis_types_map,true,_SGROUP_,FALSE);        // in kk there is the number _0 point group
              if(kk>KBIN_SYMMETRY_SGROUP_MAX_NUMBER) {
                aus << "EEEEE  ERROR: SYM::CalculateSpaceGroup cannot deal with more than " << KBIN_SYMMETRY_SGROUP_MAX_NUMBER << " space group operations." << endl;
                aus << "EEEEE  ERROR: You need to reduce the radius of the space group calculation (SGROUP_RADIUS=XXX) or" << endl;
                aus << "EEEEE  ERROR: modify KBIN_SYMMETRY_SGROUP_MAX_NUMBER in aflow.h and recompile." << endl;
                aus << "EEEEE  ERROR: WARNING: by increasing KBIN_SYMMETRY_SGROUP_MAX_NUMBER you can use a lot of memory" << endl;
                aus << "EEEEE  ERROR: and disk space if you decide to save the space group file (SGROUP_WRITE)." << endl;
                aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET,osswrite);
                return FALSE;
              }
              if(0) { // some printing
                aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "Sgroup " << a.sgroup[kk-1].str_type
                  << " theta="; if(a.sgroup[kk-1].angle<100) aus << " "; if(a.sgroup[kk-1].angle<10) aus << " ";
                aus << a.sgroup[kk-1].angle << " " << " r=(" << a.sgroup[kk-1].axis << ")"
                  << "  ftau=(" << a.sgroup[kk-1].ftau << ")"
                  << "  ftrasl=(" << a.sgroup[kk-1].ftrasl << ")"
                  << endl;  // remember vectors start from 0
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
              }
            }
          }
        }
      }
    }
    // ------------------------------------------------------------------------------
    a.sgroup_calculated=TRUE;
    // ------------------------------------------------------------------------------
    // Printing and leaving
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "SGROUP Symmetry: unique space group operations " << a.sgroup.size() << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "SGROUP Symmetry: END " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    // ------------------------------------------------------------------------------
    if(_write_) Krun=Krun && KBIN_SymmetryWrite(FileMESSAGE,a,aflags,_SGROUP_,osswrite,oss,format);
    return Krun;
  }
} // namespace SYM

#ifndef COMPILE_SLIM
namespace SYM {
  bool CalculateSpaceGroup_20160101(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss) {        // AFLOW_FUNCTION_IMPLEMENTATION
    // ------------------------------------------------------------------------------
    // Some of this routine is taken from AVDV structure::calc_factor_group() code.
    //SC made modificatios for speed and consistency with aflow architecture
    ostringstream aus;
    bool Krun=TRUE,sym_found;
    _atom hatom;
    // lattice fix to reduce the number of tested operations
    //  a.LatticeReduction(); better not repeat
    // lattice fixed, now start
    a.FixLattices();a.BringInCell();
    xvector<double> a1(3),a2(3),a3(3);                    // lattice vectors
    a1=a.lattice(1);a2=a.lattice(2);a3=a.lattice(3);      // a1,a2,a3 are the rows of the lattice matrix
    xvector<double> ctau(3),ftau(3),ctrasl(3),ftrasl(3);  // translations
    xmatrix<double> Uf(3,3),Uc(3,3);
    xvector<int> dims(3);
    std::vector<int> basis_atoms_map(a.atoms.size());       // will map like the fgroup
    std::vector<int> basis_types_map(a.atoms.size());       // will map like the fgroup
    xmatrix<double> f2c(3,3);f2c=trasp(a.lattice);
    xmatrix<double> c2f(3,3);c2f=inverse(trasp(a.lattice));
    if(a.pgroup_calculated==FALSE) Krun=Krun && SYM::CalculatePointGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss);  // NEED POINT GROUP
    if(a.fgroup_calculated==FALSE) Krun=Krun && SYM::CalculateFactorGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss); // NEED FACTOR GROUP

    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "SGROUP Symmetry: BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    // clean up sgroup                                // just initialize
    a.sgroup.clear();
    a.sgroup_calculated=FALSE;
    if(a.sgroup_radius<=0.1) {
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "SGROUP Symmetry: picking default normalized radius = " << _calculate_symmetry_default_sgroup_radius_<< "  " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
      a.sgroup_radius=_calculate_symmetry_default_sgroup_radius_*MaxStructureLattice(a);
    }
    dims=LatticeDimensionSphere(a.lattice,a.sgroup_radius);
    // for(i=1;i<=3;i++) if(dims(i)>3)dims(i)=3;
    a.sgroup_radius_dims=dims;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "SGROUP Symmetry: inside sphere (radius=" << a.sgroup_radius << ", dims = [" << dims[1] << "," << dims[2] << "," << dims[3] << "]) " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    // all symmetry operations are done within the fractional coordinate system
    // since translations back into the unit cell are straightforward
    // for each point group operation, apply it to the crystal and store the transformed
    // coordinates in atomlist
    //  dims.clear();
    for(int i=-dims[1];i<=dims[1];i++) {
      for(int j=-dims[2];j<=dims[2];j++) {
        for(int k=-dims[3];k<=dims[3];k++) {
          ctrasl=((double)i)*a1+((double)j)*a2+((double)k)*a3;
          ftrasl(1)=(double)i;ftrasl(2)=(double) j;ftrasl(3)=(double) k;
          for(uint fg=0; fg<a.fgroup.size(); fg++) {
            Uc=a.fgroup[fg].Uc;
            Uf=a.fgroup[fg].Uf;
            ctau=a.fgroup[fg].ctau;
            ftau=a.fgroup[fg].ftau;
            basis_atoms_map=a.fgroup[fg].basis_atoms_map;
            basis_types_map=a.fgroup[fg].basis_types_map;
            sym_found=FALSE;
            // no check required, as they are different by construction (not really)
            // 	  for(uint ii=0;ii<a.sgroup.size()&&!sym_found;ii++)
            // 	    //	    sym_found=(identical(Uf,a.sgroup[ii].Uf,_EPS_) && identical(ftau,a.sgroup[ii].ftau,_EPS_) && identical(ftrasl,a.sgroup[ii].ftrasl,_EPS_));    // look in all the list of operations
            // 	    sym_found=(identical(Uf,a.sgroup[ii].Uf,_EPS_) && identical(ftau+ftrasl,a.sgroup[ii].ftau+a.sgroup[ii].ftrasl,_EPS_));    // look in all the list of operations
            // if the symmetry operation is new, add it to the factorgroup array
            // and update all info about the sym_op object
            if(sym_found==FALSE)                                          // new operation, generate and save it
              //	  if(sym_found==FALSE || a.sgroup.size() == 0)
            { //CO20200106 - patching for auto-indenting
              uint kk;
              kk=SYM::AddSymmetryToStructure(a,Uc,Uf,ctau,ftau,ctrasl,ftrasl,basis_atoms_map,basis_types_map,true,_SGROUP_);        // in kk there is the number _0 point group
              if(kk>KBIN_SYMMETRY_SGROUP_MAX_NUMBER) {
                aus << "EEEEE  ERROR: SYM::CalculateSpaceGroup cannot deal with more than " << KBIN_SYMMETRY_SGROUP_MAX_NUMBER << " space group operations." << endl;
                aus << "EEEEE  ERROR: You need to reduce the radius of the space group calculation (SGROUP_RADIUS=XXX) or" << endl;
                aus << "EEEEE  ERROR: modify KBIN_SYMMETRY_SGROUP_MAX_NUMBER in aflow.h and recompile." << endl;
                aus << "EEEEE  ERROR: WARNING: by increasing KBIN_SYMMETRY_SGROUP_MAX_NUMBER you can use a lot of memory" << endl;
                aus << "EEEEE  ERROR: and disk space if you decide to save the space group file (SGROUP_WRITE)." << endl;
                aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET,osswrite);
                return FALSE;
              }
              if(0) { // some printing
                aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "Sgroup " << a.sgroup[kk-1].str_type
                  << " theta="; if(a.sgroup[kk-1].angle<100) aus << " "; if(a.sgroup[kk-1].angle<10) aus << " ";
                aus << a.sgroup[kk-1].angle << " " << " r=(" << a.sgroup[kk-1].axis << ")"
                  << "  ftau=(" << a.sgroup[kk-1].ftau << ")"
                  << "  ftrasl=(" << a.sgroup[kk-1].ftrasl << ")"
                  << endl;  // remember vectors start from 0
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
              }
            }
          }
        }
      }
    }
    // ------------------------------------------------------------------------------
    a.sgroup_calculated=TRUE;
    // ------------------------------------------------------------------------------
    // Printing and leaving
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "SGROUP Symmetry: unique space group operations " << a.sgroup.size() << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "SGROUP Symmetry: END " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    // ------------------------------------------------------------------------------
    if(_write_) Krun=Krun && KBIN_SymmetryWrite(FileMESSAGE,a,aflags,_SGROUP_,osswrite,oss);
    return Krun;
  }
} // namespace SYM
#endif
//DX+CO END

// ----------------------------------------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------------------------- EQUIVALENT ATOMS
//DX+CO START
// return the number of inequivalent atoms
namespace SYM {
  //ME20200207
  bool CalculateInequivalentAtoms(xstructure& a) {
    ofstream FileMessage;
    _aflags aflags;
    return CalculateInequivalentAtoms(FileMessage, a, aflags, false, false, std::cout, "");
  }

  bool CalculateInequivalentAtoms(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss, string format) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool rely_on_basis=TRUE; //HUGE SPEED UP
    return CalculateInequivalentAtoms(FileMESSAGE,a,rely_on_basis,aflags,_write_,osswrite,oss,format);
  }
  bool CalculateInequivalentAtoms(ofstream &FileMESSAGE,xstructure &a,bool rely_on_basis,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss, string format) {        // AFLOW_FUNCTION_IMPLEMENTATION
    double _eps_=AUROSTD_NAN;
    if(a.sym_eps!=AUROSTD_NAN){ //Tolerance came from user or was calculated
      _eps_=a.sym_eps;
    }
    else {
      _eps_=defaultTolerance(a);
    }
    return CalculateInequivalentAtoms(FileMESSAGE,a,rely_on_basis,aflags,_write_,osswrite,oss,_eps_,format);
  }
  bool CalculateInequivalentAtoms(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool rely_on_basis=TRUE; //HUGE SPEED UP
    return CalculateInequivalentAtoms(FileMESSAGE,a,rely_on_basis,aflags,_write_,osswrite,oss,_eps_,format);
  }
}

namespace SYM {
  bool CalculateInequivalentAtoms(ofstream &FileMESSAGE,xstructure &a,bool rely_on_basis,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format) {        // AFLOW_FUNCTION_IMPLEMENTATION
    return CalculateInequivalentAtoms_20160801(FileMESSAGE,a,rely_on_basis,aflags,_write_,osswrite,oss,_eps_,format);
  }
  bool CalculateInequivalentAtoms_20160801(ofstream &FileMESSAGE,xstructure &a,bool rely_on_basis,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format) {        // AFLOW_FUNCTION_IMPLEMENTATION
    // Obtain the structure tolerance
    a.sym_eps=_eps_;
    bool skew = isLatticeSkewed(a.lattice,a.dist_nn_min,_eps_);

    ostringstream aus;
    bool Krun=TRUE;

    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "IATOMS ATOMS:  BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);

    a.LatticeReduction_avoid=TRUE; // so it does not mess up the min angles SC20100115
    if(a.pgroup_calculated==FALSE) Krun=Krun && SYM::CalculatePointGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss,format);  // NEED POINT GROUP
    if(a.fgroup_calculated==FALSE) Krun=Krun && SYM::CalculateFactorGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss,format); // NEED FACTOR GROUP
    //DX20210112 [OBSOLETE - sgroup not used in this function] if(a.sgroup_calculated==FALSE) Krun=Krun && SYM::CalculateSpaceGroup(FileMESSAGE,a,aflags,FALSE,osswrite,oss,format);    // NEED SPACE GROUP
    //DX NO MORE a.sgroup_calculated=true; //DX NEED TO DO THIS BETTER

    // ------------------------------------------------------------------------------
    // create table
    for(uint iat=0;iat<a.atoms.size();iat++) {
      a.atoms.at(iat).equivalent=iat;
      a.atoms.at(iat).is_inequivalent=TRUE;
    }
    // ------------------------------------------------------------------------------
    // clear things
    for(uint i=0;i<a.iatoms.size();i++){
      a.iatoms.at(i).clear();
    }
    a.iatoms.clear();
    // ------------------------------------------------------------------------------
    // renormalize table
    //MUCH FASTER, recycle data from fgroup basis
    if(rely_on_basis){
      for(uint iat1=0;iat1<a.atoms.size();iat1++) {
        for(uint iat2=iat1+1;iat2<a.atoms.size();iat2++) {
          //if(iat2!=iat1) {  //ALREADY TRUE  //[CO20200106 - close bracket for indenting]}
          if(a.atoms[iat1].type==a.atoms[iat2].type){
            if(AtomsEquivalent_Basis(a,iat1,iat2)){
              a.atoms.at(iat2).equivalent=a.atoms.at(iat1).equivalent;
              a.atoms.at(iat2).is_inequivalent=FALSE;
              iat2=a.atoms.size();
            }
          }
        }
      }
    } else {
      for(uint iat1=0;iat1<a.atoms.size();iat1++) {
        for(uint iat2=iat1+1;iat2<a.atoms.size();iat2++) {
          //if(iat2!=iat1) {  //ALREADY TRUE  //[CO20200106 - close bracket for indenting]}
          if(a.atoms[iat1].type==a.atoms[iat2].type){
            if(AtomsEquivalent(a,a.atoms.at(iat1),a.atoms.at(iat2),skew,_eps_)){
              a.atoms.at(iat2).equivalent=a.atoms.at(iat1).equivalent;
              a.atoms.at(iat2).is_inequivalent=FALSE;
              iat2=a.atoms.size();
            }
          }
        }
      }
    }
    // ------------------------------------------------------------------------------
    // create the list of inequivalent (rows)
    for(uint iat=0;iat<a.atoms.size();iat++)
      if(a.atoms.at(iat).is_inequivalent)
        a.iatoms.push_back(std::vector<int>(1,iat));
    // create the list of inequivalent (cols)
    for(uint iat1=0;iat1<a.atoms.size();iat1++) {
      if(!a.atoms.at(iat1).is_inequivalent) {
        for(uint iat2=0;iat2<a.iatoms.size();iat2++)
          if(a.atoms.at(iat1).equivalent==a.iatoms.at(iat2).at(0)) {
            a.iatoms.at(iat2).push_back(iat1);
            //ME20200207 - all atoms should be mapped with index_iatoms
            a.atoms[iat1].index_iatoms = iat2;
          }
      }
    }
    // ------------------------------------------------------------------------------
    // save the number of equivalents
    uint iequivalent=0;
    for(uint iat=0;iat<a.atoms.size();iat++) {
      if(a.atoms.at(iat).is_inequivalent) {
        a.atoms.at(iat).num_equivalents=a.iatoms.at(iequivalent).size();
        a.atoms.at(iat).index_iatoms=iequivalent;
        iequivalent++;
      }
    }
    // ------------------------------------------------------------------------------
    //  // DEBUG so far it works
    if(0) { // DEBUG so far it works
      for(uint iat1=0;iat1<a.iatoms.size();iat1++) {
        for(uint iat2=0;iat2<a.iatoms.at(iat1).size();iat2++)
          cout << a.iatoms.at(iat1).at(iat2) << " ";
        cout << endl;
      }
    }
    // ------------------------------------------------------------------------------
    a.iatoms_calculated=TRUE;
    // ------------------------------------------------------------------------------
    // Printing and leaving
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "IATOMS ATOMS: number of inequivalent atoms = " << a.iatoms.size() << "   " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    for(uint iat1=0;iat1<a.iatoms.size();iat1++) {
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "IATOMS ATOMS: [";
      if(a.iatoms.at(iat1).at(0)<10){ 
        aus << "0"; 
      }
      aus << a.iatoms.at(iat1).at(0) << "] ";
      for(uint iat2=0;iat2<a.iatoms.at(iat1).size();iat2++){
        aus << a.iatoms.at(iat1).at(iat2) << " ";
      }
      aus << endl;
    }
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "IATOMS ATOMS: END " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    // ------------------------------------------------------------------------------
    if(_write_) Krun=Krun && KBIN_SymmetryWrite(FileMESSAGE,a,aflags,_IATOMS_,osswrite,oss,format);
    return Krun;
  }
} // namespace SYM
//DX+CO END

#ifndef COMPILE_SLIM
namespace SYM {
  bool CalculateInequivalentAtoms_20160101(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss) {        // AFLOW_FUNCTION_IMPLEMENTATION
    // return the number of inequivalent atoms
    //  str=ReScale(BringInCell(_str),1.0);
    ostringstream aus;
    bool Krun=TRUE;
    a=BringInCell(a);
    //  a=BringInWignerSeitz(a);

    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "IATOMS ATOMS:  BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    double eps=MinStructureLattice(a)/100.0;

    a.LatticeReduction_avoid=TRUE; // so it does not mess up the min angles SC20100115
    if(a.pgroup_calculated==FALSE) Krun=Krun && SYM::CalculatePointGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss);  // NEED POINT GROUP
    if(a.fgroup_calculated==FALSE) Krun=Krun && SYM::CalculateFactorGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss); // NEED FACTOR GROUP
    if(a.sgroup_calculated==FALSE) Krun=Krun && SYM::CalculateSpaceGroup(FileMESSAGE,a,aflags,FALSE,osswrite,oss);    // NEED SPACE GROUP

    // ------------------------------------------------------------------------------
    // create table
    for(uint iat=0;iat<a.atoms.size();iat++) {
      a.atoms.at(iat).equivalent=iat;
      a.atoms.at(iat).is_inequivalent=TRUE;
    }
    // ------------------------------------------------------------------------------
    // clear things
    for(uint i=0;i<a.iatoms.size();i++)
      a. iatoms.at(i).clear();
    a.iatoms.clear();
    // ------------------------------------------------------------------------------
    // renormalize table
    for(uint iat1=0;iat1<a.atoms.size();iat1++) {
      //  a.atoms.at(iat1).num_equivalents=0;
      for(uint iat2=iat1+1;iat2<a.atoms.size();iat2++) {
        if(iat2!=iat1) {
          if(SYM::AtomsEquivalent(a,a.atoms.at(iat1),a.atoms.at(iat2),eps)) {
            a.atoms.at(iat2).equivalent=a.atoms.at(iat1).equivalent;
            //	  if(a.atoms.at(iat1).num_equivalents==0) a.atoms.at(iat1).num_equivalents++;
            //	  a.atoms.at(iat1).num_equivalents++;
            a.atoms.at(iat2).is_inequivalent=FALSE;
            iat2=a.atoms.size();
          }
        }
      }
    }
    // ------------------------------------------------------------------------------
    // create the list of inequivalent (rows)
    for(uint iat=0;iat<a.atoms.size();iat++)
      if(a.atoms.at(iat).is_inequivalent)
        a.iatoms.push_back(std::vector<int>(1,iat));
    // create the list of inequivalent (cols)
    for(uint iat1=0;iat1<a.atoms.size();iat1++) {
      if(!a.atoms.at(iat1).is_inequivalent) {
        for(uint iat2=0;iat2<a.iatoms.size();iat2++)
          if(a.atoms.at(iat1).equivalent==a.iatoms.at(iat2).at(0)) {
            a.iatoms.at(iat2).push_back(iat1);
          }
      }
    }
    // ------------------------------------------------------------------------------
    // save the number of equivalents
    uint iequivalent=0;
    for(uint iat=0;iat<a.atoms.size();iat++) {
      if(a.atoms.at(iat).is_inequivalent) {
        a.atoms.at(iat).num_equivalents=a.iatoms.at(iequivalent).size();
        a.atoms.at(iat).index_iatoms=iequivalent;
        iequivalent++;
      }
    }
    // ------------------------------------------------------------------------------
    //  // DEBUG so far it works
    if(0) { // DEBUG so far it works
      for(uint iat1=0;iat1<a.iatoms.size();iat1++) {
        for(uint iat2=0;iat2<a.iatoms.at(iat1).size();iat2++)
          cout << a.iatoms.at(iat1).at(iat2) << " ";
        cout << endl;
      }
    }
    // ------------------------------------------------------------------------------
    a.iatoms_calculated=TRUE;
    // ------------------------------------------------------------------------------
    // Printing and leaving
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "IATOMS ATOMS: number of inequivalent atoms = " << a.iatoms.size() << "   " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    for(uint iat1=0;iat1<a.iatoms.size();iat1++) {
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "IATOMS ATOMS: [";
      if(a.iatoms.at(iat1).at(0)<10){ 
        aus << "0"; 
      }
      aus << a.iatoms.at(iat1).at(0) << "] ";
      for(uint iat2=0;iat2<a.iatoms.at(iat1).size();iat2++){
        aus << a.iatoms.at(iat1).at(iat2) << " ";
      }
      aus << endl;
    }
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "IATOMS ATOMS: END " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    // ------------------------------------------------------------------------------
    if(_write_) Krun=Krun && KBIN_SymmetryWrite(FileMESSAGE,a,aflags,_IATOMS_,osswrite,oss);
    return Krun;
  }
} // namespace SYM
#endif
//DX+CO END

// ----------------------------------------------------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------------------ SITE POINT GROUP OPERATIONS
// SYM::CalculateSitePointGroup
//
// This function calculate the site point group for every atom in the structure
//
//DX+CO START
namespace SYM {
  bool CalculateSitePointGroup(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss, string format) {
    //bool iatoms_only=TRUE;
    int CALCULATION_MODE=0; //default
    return CalculateSitePointGroup(FileMESSAGE,a,CALCULATION_MODE,aflags,_write_,osswrite,oss,format);
  }
  bool CalculateSitePointGroup(ofstream &FileMESSAGE,xstructure &a,int CALCULATION_MODE,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format) {
    double _eps_=AUROSTD_NAN;//=_EPS_;
    if(a.sym_eps!=AUROSTD_NAN){ //Tolerance came from user or was calculated
      _eps_=a.sym_eps;
    }
    else {
      _eps_=defaultTolerance(a);
    }
    return SYM::CalculateSitePointGroup_20160801(FileMESSAGE,a,CALCULATION_MODE,aflags,_write_,osswrite,oss,_eps_,format);
  }
} // namespace SYM

namespace SYM {
  bool CalculateSitePointGroup_20160801(ofstream &FileMESSAGE,xstructure &a,int CALCULATION_MODE,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;
    DEBUG_SYMMETRY=DEBUG_SYMMETRY || LDEBUG;
    string directory=aurostd::getPWD(); //DX20180426 - added current working directory  //[CO20191112 - OBSOLETE]aurostd::execute2string("pwd")
    //CALCULATION_MODE == 0 (default) - calculate iatoms first, then propagate to equivalent atoms using fgroups, get full basis for all
    //CALCULATION_MODE == 1 - calculate iatoms first, then propagate to equivalent atoms using fgroups, get full basis for iatoms ONLY
    //CALCULATION_MODE == 2 - calculate all atoms standard routine (go through pgroups), get full basis for all
    if(!(CALCULATION_MODE == 0 || CALCULATION_MODE == 1 || CALCULATION_MODE == 2)) {
      message << "Invalid calculation mode.  Must be 0, 1, or 2. [dir=" << a.directory << "]";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_);
    }

    // Obtain the structure tolerance
    a.sym_eps=_eps_; 
    //xstructure b(a);  //NEVER COPY AN XSTRUCTURE UNLESS IT'S ABSOLUTELY NECESSARY
    deque<_atom> b_atoms=a.atoms;
    bool skew = isLatticeSkewed(a.lattice,a.dist_nn_min,_eps_);

    ostringstream aus;
    bool Krun=TRUE;
    std::vector<int> basis_atoms_map;       // will map like the fgroup
    std::vector<int> basis_types_map;       // will map like the fgroup

    if(a.pgroup_calculated==FALSE) Krun=Krun && SYM::CalculatePointGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss,format);        // NEED POINT GROUP
    if(a.fgroup_calculated==FALSE) Krun=Krun && SYM::CalculateFactorGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss,format);       // NEED FACTOR GROUP
    if(a.pgroup_xtal_calculated==FALSE) Krun=Krun && SYM::CalculatePointGroupCrystal(FileMESSAGE,a,aflags,_write_,osswrite,oss,format);    // NEED POINT GROUP CRYSTAL
    //if(a.sgroup_calculated==FALSE) Krun=Krun && SYM::CalculateSpaceGroup(FileMESSAGE,a,aflags,FALSE,osswrite,oss);          // NEED SPACE GROUP
    if(a.iatoms_calculated==FALSE) Krun=Krun && SYM::CalculateInequivalentAtoms(FileMESSAGE,a,aflags,_write_,osswrite,oss,format); // NEED EQUIV ATOMS
    //a.write_inequivalent_flag=TRUE;

    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP Symmetry: BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);

    // ------------------------------------------------------------------------------
    // clear a little and create the framework
    for(uint i=0;i<a.agroup.size();i++) {a.agroup.at(i).clear();}
    a.agroup.clear();
    for(uint i=0;i<a.atoms.size();i++){
      a.agroup.push_back(vector<_sym_op>(0));
    }
    // ------------------------------------------------------------------------------
    // seeking for the
    uint iat;
    uint kk=0;
    bool sym_found;
    _sym_op aSymOp,fSymOp;
    _atom tatom;

    bool iatoms_only=(CALCULATION_MODE==0 || CALCULATION_MODE==1);
    uint atoms_range=(iatoms_only ? a.iatoms.size() : a.atoms.size());  //CO, if iatoms_only, only go over iatoms, all atoms otherwise 
    xvector<double> origin(3),frigin(3);
    string ATOM_INDICATOR=( iatoms_only ? "ineq atom #=" : "atom #=" );
    for(uint iiat=0;iiat<atoms_range;iiat++) {
      iat=( iatoms_only ? a.iatoms.at(iiat).at(0) : iiat );
      //if(iatoms_only){
      //  iat=a.iatoms.at(iiat).at(0);
      //} else {
      //  iat=iiat;
      //}
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP ---------------------------------------------------------------------------" << endl;
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP " << "[site="; if(iat<10) aus << "0"; aus << iat << "] ";
      aus << ATOM_INDICATOR << iiat << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
      if(0) {if(iat<10) oss << "0"; oss << iat << "["; if(a.atoms.at(iat).equivalent<10) oss << "0";  oss << a.atoms.at(iat).equivalent << "]";}
      //reset from original before flipping everything around
      origin=a.atoms[iat].cpos;
      frigin=a.atoms[iat].fpos;
      for(uint ii=0;ii<b_atoms.size();ii++){
        //go back to original first, then subtract new origin
        b_atoms[ii].cpos=a.atoms[ii].cpos-origin;
        b_atoms[ii].fpos=a.atoms[ii].fpos-frigin;
        //now bring in cell
        b_atoms[ii]=BringInCell(b_atoms[ii],a.lattice);
      }
      //b.ShiftOriginToAtom(iat);
      //b.BringInCell();

      //DX+CO, we should use PGROUP_XTAL instead of PGROUP to reduce the number of checks
      //changing pgroup to pgroup_xtal should NOT affect results, as generateFullSymBasis matches atoms in this case
      //just reduce number of checks overall
      //site symmetry is crystal symmetry dependent (match atoms), not just lattice symmetry (match points)
      //crystal symmetry <= lattice symmetry
      for(uint pg=0;pg<a.pgroup_xtal.size();pg++){  //DX+CO
        if(getFullSymBasis(b_atoms,a.lattice,a.c2f,a.f2c,a.pgroup_xtal[pg],TRUE,skew,_eps_,basis_atoms_map,basis_types_map)){ //DX+CO
          // found a site point group
          sym_found=FALSE;
          for(uint ii=0;ii<a.agroup.at(iat).size()&&!sym_found;ii++){
            sym_found=identical(a.pgroup_xtal[pg].Uf,a.agroup.at(iat)[ii].Uf);       // look in all the list of operations //DX+CO  //DX20171207 - Use Uf (only integers) not Uc and use xmatrix identical eps
          }
          // if the symmetry operation is new, add it to the pointgroup array
          // and update all info about the sym_op object
          if(sym_found==FALSE) {                                 // new operation, generate and save it
            kk=SYM::AddSymmetryToStructure(a,iat,a.pgroup_xtal[pg].Uc,a.pgroup_xtal[pg].Uf,basis_atoms_map,basis_types_map,true,_AGROUP_,FALSE);  // kk number pgroup_xtal //DX uncommented //DX+CO
            // copy over the equivalent

            // some verbose  
            //	  cerr << "F " << iiat << " " <<  pg << "/" << a.pgroup_xtal.size() << " " << kk << endl; //DX+CO
            aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP " << "[site="; if(iat<10) aus << "0"; aus << iat << "] ";
            aus << a.agroup[iat][kk-1].str_type << " theta=";                                        // remember vectors start from 0
            if(a.agroup[iat][kk-1].angle<100) aus << " ";                                            // remember vectors start from 0
            if(a.agroup[iat][kk-1].angle<10) aus << " ";                                             // remember vectors start from 0
            aus << a.agroup[iat][kk-1].angle << " " << " r=(" << a.agroup[iat][kk-1].axis << ")";    // remember vectors start from 0
            aus << "    HM=" <<  aurostd::PaddedPRE(a.agroup[iat][kk-1].str_Hermann_Mauguin,2," ");  // remember vectors start from 0
            aus << "    S=" <<  aurostd::PaddedPRE(a.agroup[iat][kk-1].str_Schoenflies,2," ");       // remember vectors start from 0
            aus << endl;                                                                             // remember vectors start from 0
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
          }
        }
      }
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP " << "[site="; if(iat<10) aus << "0"; aus << iat << "] ";
      //aus << "unique site point group operations " << a.agroup.at(iat).size() << " (#pg=" << pgroup.size() << ")  " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
      aus << "unique site point group operations " << a.agroup.at(iat).size() << " (#pg_xtal=" << a.pgroup_xtal.size() << ")  " << endl;// Message(__AFLOW_FILE__,aflags) << endl;  //DX+CO
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);    
    }
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP " << SEPARATION_LINE_DASH_SHORT << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);    
    // ------------------------------------------------------------------------------
    if(iatoms_only){
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP SYMMETRY: propagating to equivalent sites" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);    
      if(!CalculateSitePointGroup_EquivalentSites(a,CALCULATION_MODE!=1,_eps_)){
        return FALSE; 
      }
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP " << SEPARATION_LINE_DASH_SHORT << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);    
    }
    //since iatoms calculated, for (CALCULATION_MODE == 2) we can do a quick validatity check
    if(CALCULATION_MODE == 2){
      uint count_iat;
      for(uint i=0;i<a.iatoms.size();i++){
        count_iat=a.agroup[a.iatoms[i][0]].size();
        for(uint j=1;j<a.iatoms[i].size();j++){
          if(count_iat!=a.agroup[a.iatoms[i][j]].size()){
            if(DEBUG_SYMMETRY){                                
              cerr << "SYM::CalculateSitePointGroup warning - Mismatch in agroup count between atom " << a.iatoms[i][0] << " and atom " << a.iatoms[i][j] << " [dir=" << a.directory << "]" << endl;
            }
            return FALSE;
          }
        }
      }
    }
    // ------------------------------------------------------------------------------
    a.agroup_calculated=TRUE;
    // ------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------
    // summary
    // ------------------------------------------------------------------------------
    // Printing and leaving
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP number of inequivalent atoms = " << a.iatoms.size() << "   " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    for(uint iat1=0;iat1<a.iatoms.size();iat1++) {
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP [";
      if(a.iatoms.at(iat1).at(0)<10) aus << "0"; 
      aus << a.iatoms.at(iat1).at(0) << "] ";
      for(uint iat2=0;iat2<a.iatoms.at(iat1).size();iat2++)
        aus << a.iatoms.at(iat1).at(iat2) << "(" << a.agroup.at(a.iatoms.at(iat1).at(iat2)).size() << ")" << " ";
      aus << endl;
    }
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP Symmetry: END " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    if(_write_) Krun=Krun && KBIN_SymmetryWrite(FileMESSAGE,a,aflags,_AGROUP_,osswrite,oss,format);
    return Krun;
  }
} // namespace SYM
//DX+CO END

#ifndef COMPILE_SLIM
namespace SYM {
  bool CalculateSitePointGroup_20160101(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_) {        // AFLOW_FUNCTION_IMPLEMENTATION
    // Written by SC, dec 2007
    ostringstream aus;
    bool Krun=TRUE;
    std::vector<int> basis_atoms_map(a.atoms.size());       // will map like the fgroup
    std::vector<int> basis_types_map(a.atoms.size());       // will map like the fgroup

    xstructure b(a),t(a);   // for backup... will use at the end to plug back the postions
    a.BringInCell();
    a.BringInWignerSeitz();

    if(a.pgroup_calculated==FALSE) Krun=Krun && SYM::CalculatePointGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss);        // NEED POINT GROUP
    if(a.fgroup_calculated==FALSE) Krun=Krun && SYM::CalculateFactorGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss);       // NEED FACTOR GROUP
    if(a.sgroup_calculated==FALSE) Krun=Krun && SYM::CalculateSpaceGroup(FileMESSAGE,a,aflags,FALSE,osswrite,oss);          // NEED SPACE GROUP
    if(a.iatoms_calculated==FALSE) Krun=Krun && SYM::CalculateInequivalentAtoms(FileMESSAGE,a,aflags,_write_,osswrite,oss); // NEED EQUIV ATOMS
    a.write_inequivalent_flag=TRUE;

    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP Symmetry: BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);

    //  double _eps_=MaxStructureLattice(a)/1000.0;
    //  double _eps_=_EPS_; optional

    // ------------------------------------------------------------------------------
    // clear a little and create the framework
    for(uint i=0;i<a.agroup.size();i++) a.agroup.at(i).clear();
    a.agroup.clear();
    for(uint i=0;i<a.atoms.size();i++)
      a.agroup.push_back(vector<_sym_op>(0));
    // ------------------------------------------------------------------------------
    // seeking for the
    uint iat;
    uint isequal,kk=0;
    bool sym_found;

    for(uint iiat=0;iiat<a.iatoms.size();iiat++) {
      iat=a.iatoms.at(iiat).at(0);
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP ---------------------------------------------------------------------------" << endl;
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP " << "[site="; if(iat<10) aus << "0"; aus << iat << "] ";
      aus << "ineq atom #=" << iiat << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
      if(0) {if(iat<10) oss << "0"; oss << iat << "["; if(a.atoms.at(iat).equivalent<10) oss << "0";  oss << a.atoms.at(iat).equivalent << "]";}
      a.ShiftOriginToAtom(iat);
      //    a.BringInWignerSeitz();
      a.BringInCell();

      for(uint pg=0;pg<a.pgroup.size();pg++) {
        std::deque<_atom> atomlist;
        for(uint i=0;i<a.atoms.size();i++)
          t.atoms[i]=SYM::ApplyAtom(a.atoms[i],a.pgroup[pg],a,FALSE);
        t.BringInCell();
        //  t.BringInWignerSeitz();
        isequal=0;
        for(uint i1=0;i1<basis_atoms_map.size();i1++) 
          basis_atoms_map[i1]=0;
        for(uint i1=0;i1<a.atoms.size();i1++)
          for(uint i2=0;i2<a.atoms.size();i2++)
            if(modulus(a.atoms[i1].fpos-t.atoms[i2].fpos)<_eps_) {
              basis_atoms_map[i1]=i2;
              isequal++;
            }
        if(isequal==a.atoms.size()) {
          // found a site point group
          sym_found=FALSE;
          for(uint ii=0;ii<a.agroup.at(iat).size()&&!sym_found;ii++)
            sym_found=identical(a.pgroup[pg].Uc,a.agroup.at(iat)[ii].Uc,_eps_);       // look in all the list of operations
          // if the symmetry operation is new, add it to the pointgroup array
          // and update all info about the sym_op object
          if(sym_found==FALSE) {                                 // new operation, generate and save it
            // kk=SYM::AddSymmetryToStructure(a,iat,a.pgroup[pg].Uc,a.pgroup[pg].Uf,basis_atoms_map,basis_types_map,_AGROUP_);  // kk number pgroups
            // copy over the equivalent
            for(uint iat2=0;iat2<a.iatoms.at(iiat).size();iat2++)
              kk=SYM::AddSymmetryToStructure(a,a.iatoms.at(iiat).at(iat2),a.pgroup[pg].Uc,a.pgroup[pg].Uf,basis_atoms_map,basis_types_map,true,_AGROUP_);  // kk number pgroups
            //aus << a.iatoms.at(iiat).at(iat2) << " ";

            // some verbose  
            //	  cerr << "F " << iiat << " " <<  pg << "/" << a.pgroup.size() << " " << kk << endl;
            aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP " << "[site="; if(iat<10) aus << "0"; aus << iat << "] ";
            aus << a.agroup[iat][kk-1].str_type << " theta=";                                        // remember vectors start from 0
            if(a.agroup[iat][kk-1].angle<100) aus << " ";                                            // remember vectors start from 0
            if(a.agroup[iat][kk-1].angle<10) aus << " ";                                             // remember vectors start from 0
            aus << a.agroup[iat][kk-1].angle << " " << " r=(" << a.agroup[iat][kk-1].axis << ")";    // remember vectors start from 0
            aus << "    HM=" <<  aurostd::PaddedPRE(a.agroup[iat][kk-1].str_Hermann_Mauguin,2," ");  // remember vectors start from 0
            aus << "    S=" <<  aurostd::PaddedPRE(a.agroup[iat][kk-1].str_Schoenflies,2," ");       // remember vectors start from 0
            aus << endl;                                                                             // remember vectors start from 0
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
          }
        }
      }
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP " << "[site="; if(iat<10) aus << "0"; aus << iat << "] ";
      aus << "unique site point group operations " << a.agroup.at(iat).size() << " (#pg=" << a.pgroup.size() << ")  " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);    
    }
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP ---------------------------------------------------------------------------" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);    
    // ------------------------------------------------------------------------------
    for(uint iat=0;iat<a.atoms.size();iat++) {
      a.atoms.at(iat).cpos=b.atoms.at(iat).cpos;
      a.atoms.at(iat).fpos=b.atoms.at(iat).fpos;
    }
    // ------------------------------------------------------------------------------
    a.agroup_calculated=TRUE;
    // ------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------
    // summary
    // ------------------------------------------------------------------------------
    // Printing and leaving
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP number of inequivalent atoms = " << a.iatoms.size() << "   " << endl;// Message(__AFLOW_FILE__,aflags) << endl;
    for(uint iat1=0;iat1<a.iatoms.size();iat1++) {
      aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP [";
      if(a.iatoms.at(iat1).at(0)<10) aus << "0";
      aus << a.iatoms.at(iat1).at(0) << "] ";
      for(uint iat2=0;iat2<a.iatoms.at(iat1).size();iat2++)
        aus << a.iatoms.at(iat1).at(iat2) << "(" << a.agroup.at(a.iatoms.at(iat1).at(iat2)).size() << ")" << " ";
      aus << endl;
    }
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP Symmetry: END " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    if(_write_) Krun=Krun && KBIN_SymmetryWrite(FileMESSAGE,a,aflags,_AGROUP_,osswrite,oss);
    return Krun;
  }
} // namespace SYM
#endif

//DX+CO START 
namespace SYM {
  //CO
  //ADD FLAG TO IGNORE GETFULLSYMBASISERRORS, KEEP MOVING
  //CO
  bool CalculateSitePointGroup_EquivalentSites(xstructure &a,double _eps_){
    bool get_full_basis=true;
    return CalculateSitePointGroup_EquivalentSites(a,get_full_basis,_eps_);
  }
  bool CalculateSitePointGroup_EquivalentSites(xstructure &a,bool get_full_basis,double _eps_){
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    DEBUG_SYMMETRY=DEBUG_SYMMETRY || LDEBUG;
    string directory=aurostd::getPWD(); //DX20180426 - added current working directory  //[CO20191112 - OBSOLETE]aurostd::execute2string("pwd")

    // Obtain the structure tolerance
    a.sym_eps=_eps_;
    //xstructure b(a);  //NEVER COPY AN XSTRUCTURE UNLESS IT'S ABSOLUTELY NECESSARY
    deque<_atom> b_atoms=a.atoms;
    bool skew = isLatticeSkewed(a.lattice,a.dist_nn_min,_eps_);

    //assumes that agroup already found for agroup.at(i).at(0), for all inequivalent sites
    if(a.pgroup_calculated==FALSE) return FALSE;       // NEED FACTOR GROUP
    if(a.fgroup_calculated==FALSE) return FALSE;       // NEED FACTOR GROUP
    if(a.pgroup_xtal_calculated==FALSE) return FALSE;       // NEED FACTOR GROUP
    if(a.iatoms_calculated==FALSE) return FALSE;       // NEED IATOMS
    if((!a.fgroup.size())||(!a.iatoms.size())||(!a.agroup.size())) return FALSE;
    for(uint iiat=0;iiat<a.iatoms.size();iiat++) {
      if(!a.agroup.at(a.iatoms.at(iiat).at(0)).size()) {return FALSE;}   // NEED AGROUP FOR INEQUIVALENT ATOMS
      for(uint ii=1;ii<a.iatoms.at(iiat).size();ii++){
        a.agroup.at(a.iatoms.at(iiat).at(ii)).clear();  //clear otherwise
      }
    }

    uint iat=0;//,at_indx;
    uint eat=0;
    _sym_op aSymOp,fSymOp;

    bool found_fgroup=false;
    bool fgroup_map;
    _atom tatom;
    xvector<double> origin(3),frigin(3);
    for(uint iiat=0;iiat<a.iatoms.size();iiat++) {
      iat=a.iatoms[iiat][0];
      for(uint ieat=1;ieat<a.iatoms[iiat].size();ieat++){
        found_fgroup=false;
        eat=a.iatoms[iiat][ieat];
        for(uint fg=0;fg<a.fgroup.size() && !found_fgroup;fg++){
          if(a.fgroup[fg].basis_map_calculated){
            //DX [OBSOLETE] fgroup_map=(uint)a.fgroup[fg].basis_atoms_map[eat]==iat; //DX20170731 - Incorrect, want centered at iatom not eatom
            fgroup_map=(uint)a.fgroup[fg].basis_atoms_map[iat]==eat; //DX20170731 find fgroup centered on iatom
          } else {
            if(!SYM::ApplyAtomValidate(a.atoms[iat],tatom,a.fgroup[fg],a.lattice,a.c2f,a.f2c,skew,FALSE,FALSE,_eps_)){
              return FALSE;
            }
            //DX20190905 [OBSOLETE-no more mod_one_xvec] tatom.fpos=mod_one_xvec(tatom.fpos);
            BringInCellInPlace(tatom.fpos); //DX20190905 - uses new BringInCell function
            fgroup_map=MapAtom(tatom, a.atoms[eat], FALSE, a.lattice, a.f2c, skew, _eps_); //DX20190619 - lattice and f2c as input
          }
          if(fgroup_map){
            fSymOp=a.fgroup[fg];
            found_fgroup=true;
            //break;
          }
        }
        if(found_fgroup){
          for(uint ag=0;ag<a.agroup[iat].size();ag++){
            aSymOp=a.agroup[iat][ag];
            aSymOp.Uc= fSymOp.Uc * aSymOp.Uc * inverse(fSymOp.Uc);    //convert basis to that of fgroup //similarity transformation
            aSymOp.Uf= fSymOp.Uf * aSymOp.Uf * inverse(fSymOp.Uf);    //convert basis to that of fgroup //similarity transformation
            if(get_full_basis){
              //update basis_atoms_map and basis_types_map
              //reset from original before flipping everything around
              origin=a.atoms[eat].cpos;
              frigin=a.atoms[eat].fpos;
              for(uint ii=0;ii<b_atoms.size();ii++){
                //go back to original first, then subtract new origin
                b_atoms[ii].cpos=a.atoms[ii].cpos-origin;
                b_atoms[ii].fpos=a.atoms[ii].fpos-frigin;
                //now bring in cell
                b_atoms[ii]=BringInCell(b_atoms[ii],a.lattice);
              }
              //b.ShiftOriginToAtom(eat);
              //b.BringInCell();
              if(!getFullSymBasis(b_atoms,a.lattice,a.c2f,a.f2c,aSymOp,TRUE,skew,_eps_,aSymOp.basis_atoms_map,aSymOp.basis_types_map)){
                //if(derivative_structure){
                //  continue;
                //} else {
                if(DEBUG_SYMMETRY){ 
                  cerr << __AFLOW_FUNC__ << " warning[1] - Cannot find full atom/types basis [dir=" << a.directory << "]" << endl;
                }
                return FALSE; //CO REMOVE
                //}
              }
            } else {
              for(uint iii=0;iii<a.atoms.size();iii++){
                aSymOp.basis_atoms_map.push_back(0);
                aSymOp.basis_types_map.push_back(0);
              }
              aSymOp.basis_map_calculated=false;
            }
            //a.agroup[eat].push_back(aSymOp);
            SYM::AddSymmetryToStructure(a,eat,aSymOp.Uc,aSymOp.Uf,aSymOp.ctau,aSymOp.ftau,aSymOp.ctrasl,aSymOp.ftrasl,aSymOp.basis_atoms_map,aSymOp.basis_types_map,aSymOp.basis_map_calculated,_AGROUP_,FALSE);  //CO20170706 - make sure quaternion is updated
          }
        } else {
          //if(!found_fgroup){  //[CO20200106 - close bracket for indenting]}
          if(DEBUG_SYMMETRY){
            cerr << __AFLOW_FUNC__ << " warning[1a] - Cannot find fgroup mapping between atom " << eat << " and atom " << iat << " [dir=" << directory << "]" << endl; //DX20180426 - changed a.directory to directory (pwd)
          }
          return FALSE;
          //cerr << "Not throwing though, just applying agroups of iatom " << iat << " to atom " << eat << endl;
          //a.agroup[eat]=a.agroup.at(iat);
          //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_);
        }
      } 
    }
    //else {  //[CO20200106 - close bracket for indenting]}
    //if basis_atoms_map was not found/stored, we need to find map from iatom to eatom manually
    //_atom tatom;
    //for(uint iiat=0;iiat<a.iatoms.size();iiat++) {
    //iat=a.iatoms.at(iiat).at(0);
    ////cerr << "IATOM ORIG " << a.atoms.at(iat).fpos << endl;
    ////perform transformation for equivalent atoms and apply
    //for(uint ii=1;ii<a.iatoms.at(iiat).size();ii++){
    //found_fgroup=false;
    //eat=a.iatoms.at(iiat).at(ii);
    ////cerr << "EATOM ORIG " << a.atoms.at(eat).fpos << endl;
    ////find fgroup that maps iatom onto the equivalent atom
    //for(uint fg=0;fg<a.fgroup.size() && !found_fgroup;fg++) {
    ////DX tatom=SYM::ApplyAtom(a.atoms.at(iat),a.fgroup[fg],a,FALSE,FALSE,derivative_structure);
    //if(!SYM::ApplyAtomValidate(a.atoms.at(iat),tatom,a.fgroup[fg],a.lattice,a.c2f,a.f2c,skew,FALSE,FALSE,_eps_)){
    //// Applying c2f and f2c to the Cartesian and fractional positions give different results (a tolerance issue) 
    //return FALSE;
    //}
    //tatom.fpos=mod_one_xvec(tatom.fpos);
    ////cerr << "IATOM ROTATED " << tatom.fpos << endl;
    //if(MapAtom(tatom, a.atoms.at(eat), FALSE,a.c2f, a.f2c, skew, _eps_)){ 
    //fSymOp=a.fgroup[fg];
    //found_fgroup=true;
    ////break;
    //}
    //}
    //if(found_fgroup){
    //for(uint ag=0;ag<a.agroup.at(iat).size();ag++){
    //aSymOp=a.agroup.at(iat).at(ag);
    //aSymOp.Uc= fSymOp.Uc * aSymOp.Uc * inverse(fSymOp.Uc);    //convert basis to that of fgroup
    //aSymOp.Uf= fSymOp.Uf * aSymOp.Uf * inverse(fSymOp.Uf);    //convert basis to that of fgroup
    //if(get_full_basis){
    ////update basis_atoms_map and basis_types_map
    ////reset from original before flipping everything around
    //for(uint ii=0;ii<b.atoms.size();ii++){
    //b.atoms.at(ii).cpos=a.atoms.at(ii).cpos;
    //b.atoms.at(ii).fpos=a.atoms.at(ii).fpos;
    //}
    //b.ShiftOriginToAtom(eat);
    //b.BringInCell();
    //if(!getFullSymBasis(b.atoms,a.lattice,a.c2f,a.f2c,aSymOp,TRUE,skew,_eps_,aSymOp.basis_atoms_map,aSymOp.basis_types_map)){
    ////if(derivative_structure){
    ////  continue;
    ////} else {
    //cerr << __AFLOW_FUNC__ << " error[2] - Cannot find full atom/types basis" << endl;
    //return FALSE;
    ////}
    //}
    //} else {
    //for(uint iii=0;iii<a.atoms.size();iii++){
    //aSymOp.basis_atoms_map.push_back(0);
    //aSymOp.basis_types_map.push_back(0);
    //}
    //aSymOp.basis_map_calculated=false;
    //}
    //a.agroup.at(eat).push_back(aSymOp);
    //}
    //} else {
    ////if(!found_fgroup){
    //cerr << __AFLOW_FUNC__ << " error[1b] - Cannot find fgroup mapping between atom " << eat << " and atom " << iat << endl;
    //return FALSE;
    ////cerr << "Not throwing though, just applying agroups of iatom " << iat << " to atom " << eat << endl;
    ////a.agroup[eat]=a.agroup.at(iat);
    ////throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_);
    //}
    //} 
    //}
    //}
    return TRUE;
  }
}
//DX+CO END


// --------------------------------------------------------------------------
namespace SYM {
  //void CalculateSitePointGroup2(xstructure &a,vector<vector<uint> > &Bsitesym, bool ComMidss)
  void CalculateSitePointGroup2(xstructure &a, bool ComMidss)
  { //CO20200106 - patching for auto-indenting
    //This function calculates site symmetry of all basis in POSCAR.
    //First, it calculates the point group of the lattice.
    //Then, 4-vectors: a1 a2 a3 of the original lattice vectors + Gamma are used
    //to verify that the rotated versions of these 4-vectors map to the lattice points
    //after the symmetry operation is done to the 4vectors.

    //The application of symmetry operation is done after the basis is used as origin.
    //Then, the map checking is done after the original (0,0,0) is used as origin.
    //After the symmetry operation, a translation is done to bring the 4-vectors within
    //cutoff sphere of the original lattice points. The translation integers are determined
    //by calculating the integers that would bring the new Gamma point inside the wigner-seitz
    //cell of the original unit cell.

    //written by wahyu@alumni.duke.edu (2009) with the block to calculate the point group taken from
    //SYM::CalculatePointGroup function.

    stringstream message;

    uint i,j,ib,jp;
    ostringstream aus;
    // generic variables
    double _eps_=_EPS_,_eps_angle_=10.0*_EPS_,_eps_abc_,_eps_vol_;
    bool sym_found;
    xvector<double> rrr(3),ddd(3);                    // for printing
    // ---------------------------------------------------------------------------
    // create placeholder for atoms_map for AddSymmetryToStructure() (no atom mappings for pgroups)
    std::vector<int> basis_atoms_map(a.atoms.size());       // will map each on each
    for(uint i=0;i<a.atoms.size();i++) basis_atoms_map[i]=i; // identically map each over each
    std::vector<int> basis_types_map(a.atoms.size());       // will map each on each
    for(uint i=0;i<a.atoms.size();i++) basis_types_map[i]=a.atoms[i].type; // identically map each over each

    // clean up pgroup                                // just initialize
    a.pgroup.clear();
    a.pgroup_calculated=FALSE;                        // just initialize
    xmatrix<double> lattice(3,3),irlatt(3,3),temp_clattice(3,3),temp_flattice(3,3);
    xmatrix<double> f2c(3,3);f2c=trasp(a.lattice);
    xmatrix<double> c2f(3,3);c2f=inverse(trasp(a.lattice));
#ifdef _NO_SCALE_LATTICE_PGROUP_
    lattice=(a.lattice);                                  // copies so it does not messup
#else
    lattice=a.scale*(a.lattice);                          // copies so it does not messup
#endif
    xmatrix<double> Act(3,3),Bct(3,3),Uc(3,3);            // for generation of pgroup
    xmatrix<double> Adt(3,3),Bdt(3,3),Uf(3,3);            // for generation of pgroup_ijk
    xvector<double> a1(3),a2(3),a3(3);                    // lattice vectors
    xvector<double> ctau(3),ftau(3),ctrasl(3),ftrasl(3);  // translation
    xvector<double> clatticedata(6),temp_clatticedata(6);
    a1=lattice(1);a2=lattice(2);a3=lattice(3);            // a1,a2,a3 are the rows of the lattice matrix
    Act=lattice;
    Adt=aurostd::identity((double) 0,3);  //ME20200123 - new identity format
    clatticedata=_Getabc_angles(lattice,DEGREES);
    vector<xvector<double>*> grid_clattice; xvector<double> *grid_clattice_ptr;  // grid for rrr points
    vector<xvector<double>*> grid_flattice; xvector<double> *grid_flattice_ptr;  // grid for ijk points
    double na1,na2,na3,amin,latticevol,temp_clatticevol;
    // [UNUSED] double amax; // warning: variable ‘amax’ set but not used
    na1=modulus(a1);na2=modulus(a2);na3=modulus(a3);
    // [UNUSED] amax=max(na1,na2,na3); // warning: variable ‘amax’ set but not used
    amin=min(na1,na2,na3);
    _eps_abc_=_eps_*amin;
    _eps_vol_=_eps_*amin*amin*amin;
    latticevol=abs(det(lattice));
    // ------------------------------------------------------------------------------
    // make a lattice parallelepiped
    xvector<int> dims(3);
    double radius=0;
    std::ofstream flog;
    flog.open("agroup2.log");
    radius=1.5*RadiusSphereLattice(lattice);
    dims=LatticeDimensionSphere(lattice,radius);
    flog<<"radius = "<<radius<<endl
      <<"dims: "<<dims<<endl
      <<"tolerance eps, epsabc: "<<_eps_<<", "<<_eps_abc_<<endl;
    grid_clattice.clear();grid_flattice.clear();
    for(uint i=1;i<=3;i++) {
      grid_clattice_ptr = new xvector<double>(3);
      *grid_clattice_ptr=Act(i);
      grid_clattice.push_back(grid_clattice_ptr);  // TRICK so I get the identity first
      grid_flattice_ptr = new xvector<double>(3);
      *grid_flattice_ptr=Adt(i);
      grid_flattice.push_back(grid_flattice_ptr);  // TRICK so I get the identity first
    }
    uint ii,jj,kk;
    for(int i=-dims[1];i<=dims[1];i++) {
      for(int j=-dims[2];j<=dims[2];j++) {
        for(int k=-dims[3];k<=dims[3];k++) {
          if(!(i==0 && j==0 && k==0)) {
            if((!(i==1 && j==0 && k==0)) && (!(i==0 && j==1 && k==0)) && (!(i==0 && j==0 && k==1))) { // to avoid plugging a1,a2,a2 back
              rrr=((double)i)*a1+((double)j)*a2+((double)k)*a3;
              ddd(1)=(double) i;ddd(2)=(double) j;ddd(3)=(double) k;
              if(modulus(rrr)<(radius+3*_eps_abc_)) {
                grid_clattice_ptr = new xvector<double>(3);        // SAVE THEM ALL
                *grid_clattice_ptr=rrr;                            // SAVE THEM ALL
                grid_clattice.push_back(grid_clattice_ptr);        // SAVE THEM ALL
                grid_flattice_ptr = new xvector<double>(3);        // SAVE THEM ALL
                *grid_flattice_ptr=ddd;                            // SAVE THEM ALL
                grid_flattice.push_back(grid_flattice_ptr);        // SAVE THEM ALL
              }
            }
          }
        }
      }
    }
    // for each set of three lattice points within the sphere see which one has the
    // same sets of lengths and angles as the original lattice unit cell vectors.
    ii=0;jj=0;kk=0;
    for(uint i=0;i<grid_clattice.size();i++) {
      for(ii=1;ii<=3;ii++) temp_clattice(1,ii)=(*grid_clattice[i])(ii);
      for(ii=1;ii<=3;ii++) temp_flattice(1,ii)=(*grid_flattice[i])(ii);
      for(uint j=0;j<grid_clattice.size();j++) {
        for(ii=1;ii<=3;ii++) temp_clattice(2,ii)=(*grid_clattice[j])(ii);
        for(ii=1;ii<=3;ii++) temp_flattice(2,ii)=(*grid_flattice[j])(ii);
        for(uint k=0;k<grid_clattice.size();k++) {
          if(i!=j && i!=k && j!=k) {
            for(ii=1;ii<=3;ii++) temp_clattice(3,ii)=(*grid_clattice[k])(ii);
            for(ii=1;ii<=3;ii++) temp_flattice(3,ii)=(*grid_flattice[k])(ii);
            temp_clatticevol=abs(det(temp_clattice));
            if(abs(temp_clatticevol-latticevol)<=_eps_vol_) {  // check volume !!! FAST !
              temp_clatticedata=_Getabc_angles(temp_clattice,DEGREES);
              jj++;//if(!mod(jj,100000)) cerr << temp_clatticedata << endl;
              // compare the temp_clattice... and lattice... to see if they are the same lattice
              // that is do the lattice vectors have the same lengths and the same angles
              if(abs(clatticedata[1]-temp_clatticedata[1]) < _eps_abc_)                // check lattices
                if(abs(clatticedata[2]-temp_clatticedata[2]) < _eps_abc_)              // check lattices
                  if(abs(clatticedata[3]-temp_clatticedata[3]) < _eps_abc_)            // check lattices
                    if(abs(clatticedata[4]-temp_clatticedata[4]) < _eps_angle_)        // check angles
                      if(abs(clatticedata[5]-temp_clatticedata[5]) < _eps_angle_)      // check angles
                        if(abs(clatticedata[6]-temp_clatticedata[6]) < _eps_angle_) {  // check angles
                          // SYMMETRY OPERATION POINT_GROUP and FACTOR_GROUP
                          // get the matrix that relates the two lattice vectors
                          // A=U*B but in A and B we plug vectors as columns
                          // watch lattice is per row Uc=A*inv(B)
                          // A is the lattice (vectors per colum), B is the test lattice (epr column)
                          // Uc is the point_group operation which operates AFTER the vector (row)
                          // as: new_vector_row=old_vector_row*Uc
                          // point_group is the list of all the Uc !!!
                          // since A and B are saved as trasposed, I get
                          // U=trasp(A)*inv(trasp(U))
                          Bct=temp_clattice;Uc=trasp(Act)*inverse(trasp(Bct)); // Act=lattice;
                          Bdt=temp_flattice;Uf=trasp(Adt)*inverse(trasp(Bdt)); // Adt=I3... could remove !
                          // check whether this symmetry operation is new or not
                          sym_found=FALSE;
                          for(ii=0;ii<a.pgroup.size()&&!sym_found;ii++)
                            sym_found=identical(Uf,a.pgroup[ii].Uf);       // look in all the list of operations  //DX20171207 - Use Uf (only integers) not Uc and use xmatrix identical eps
                          // if the symmetry operation is new, add it to the pointgroup array
                          // and update all info about the sym_op object
                          if(sym_found==FALSE) {                                 // new operation, generate and save it
                            clear(ctau);clear(ftau);clear(ctrasl);clear(ftrasl);
                            kk=SYM::AddSymmetryToStructure(a,Uc,Uf,ctau,ftau,ctrasl,ftrasl,basis_atoms_map,basis_types_map,false,_PGROUP_);  // kk number pgroups
                            aus << a.pgroup[kk-1].str_type
                              << " theta=";
                            if(a.pgroup[kk-1].angle<100) aus << " ";
                            if(a.pgroup[kk-1].angle<10) aus << " ";
                            aus << a.pgroup[kk-1].angle << " " << " r=(" << a.pgroup[kk-1].axis << ")"
                              << endl;  // remember vectors start from 0
                            if(aurostd::sum(aurostd::abs(f2c*Uf*inverse(f2c)-Uc))>_eps_)
                            {message << "Uf error[1]"; throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);}
                            if(aurostd::sum(aurostd::abs(c2f*Uc*inverse(c2f)-Uf))>_eps_)
                            {message << "Uf error[2]"; throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);}
                            if(aurostd::sum(aurostd::abs(FF2CC(a.lattice,Uf)-Uc))>_eps_)
                            {message << "Uf error[3]"; throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);}
                            if(aurostd::sum(aurostd::abs(CC2FF(a.lattice,Uc)-Uf))>_eps_)
                            {message << "Uf error[4]"; throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);}
                          }
                        }
            }
          }
        }
      }
    }
    // ------------------------------------------------------------------------------
    a.pgroup_calculated=TRUE;
    // ------------------------------------------------------------------------------
    for(uint k=0;k<a.pgroup.size();k++) {
      Uc=a.pgroup[k].Uc;
      Uf=a.pgroup[k].Uf;
      roundoff(a.pgroup[k].Uc,_EPS_roundoff_);
      //for(i=1;i<=3;i++)
      //  for(j=1;j<=3;j++)
      //  if(abs((*a.pgroup[k])(i,j))<_EPS_roundoff_) (*a.pgroup[k])(i,j)=0.0;
      //for(i=0;i<9;i++) if(abs((*a.pgroup[k])(int(i/3)+1,mod(i,3)+1))<1.0e-10) (*a.pgroup[k])(int(i/3)+1,mod(i,3)+1)=0.0;
    }

    cout<<"Point symmtery operations:"<<endl;
    for(i=0;i<a.pgroup.size();i++) {
      cout<<"pointsym_op "<<i<<"/"<<a.pgroup.size()<<": "
        <<a.pgroup[i].str_type<<", "<<a.pgroup[i].angle<<" about "
        <<a.pgroup[i].axis<<endl
        <<a.pgroup[i].Uc<<endl;
      flog<<"pointsym_op "<<i<<"/"<<a.pgroup.size()<<": "
        <<a.pgroup[i].str_type<<", "<<a.pgroup[i].angle<<" about "
        <<a.pgroup[i].axis<<endl
        <<a.pgroup[i].Uc<<endl;
    }

    bool foundG,founda1,founda2,founda3;
    uint Nbasis,it;
    double n1,n2,n3;
    Nbasis=a.atoms.size();
    xmatrix<double> rlatt(3,3),rlattb4Uc(3,3),rlatttmp(3,3),mG(3,3),mGb4Uc(3,3);
    vector<xvector<double> > basis(Nbasis);
    vector<vector<uint > > Bsitesym(Nbasis);//site symmetry of basis
    for(ib=0;ib<Nbasis;ib++) {
      basis[ib]=a.atoms[ib].cpos;
      flog<<"basis["<<ib<<"]: "<<basis[ib]<<endl;
    }

    mGb4Uc(1,1)=0.0; mGb4Uc(1,2)=0.0; mGb4Uc(1,3)=0.0;
    mGb4Uc(3)=mGb4Uc(2)=mGb4Uc(1);
    a1=lattice(1); a2=lattice(2); a3=lattice(3);
    for(ib=0;ib<Nbasis;ib++) {
      Bsitesym[ib].clear();
      //r=r-A
      flog<<"--- with respect to basis["<<ib<<"]: "<<basis[ib]<<endl;
      flog<<"rlatt row-lattice vectors orig:"<<endl<<lattice<<endl;
      for(it=1;it<4;it++) rlattb4Uc(1,it)=lattice(1,it)-basis[ib][it];
      for(it=1;it<4;it++) rlattb4Uc(2,it)=lattice(2,it)-basis[ib][it];
      for(it=1;it<4;it++) rlattb4Uc(3,it)=lattice(3,it)-basis[ib][it];
      for(it=1;it<4;it++) mGb4Uc(1,it)=-basis[ib][it];
      mGb4Uc(3)=mGb4Uc(2)=mGb4Uc(1);
      flog<<"rlatt r=r-basis["<<ib<<"]:"<<endl<<rlattb4Uc<<endl
        <<"gamma g=g-basis["<<ib<<"]: "<<mGb4Uc(1)<<endl;
      for(jp=0;jp<a.pgroup.size();jp++) {
        //r=Ur
        rlatttmp=a.pgroup[jp].Uc*(trasp(rlattb4Uc));
        rlatt=trasp(rlatttmp);
        rlatttmp=a.pgroup[jp].Uc*(trasp(mGb4Uc));
        mG=trasp(rlatttmp);
        flog<<"try sym_op "<<jp<<"/"<<a.pgroup.size()<<", rotation matrix U:"<<endl<<a.pgroup[jp].Uc<<endl;
        flog<<"rlatt r=Ur:"<<endl<<rlatt<<endl
          <<"gamma  g=U*g: "<<mG(1)<<endl;
        //r=r+A
        for(it=1;it<4;it++) rlatt(1,it)=rlatt(1,it)+basis[ib][it];
        for(it=1;it<4;it++) rlatt(2,it)=rlatt(2,it)+basis[ib][it];
        for(it=1;it<4;it++) rlatt(3,it)=rlatt(3,it)+basis[ib][it];
        for(it=1;it<4;it++) mG(1,it)=mG(1,it)+basis[ib][it];
        flog<<"rlatt r=r+basis["<<ib<<"]:"<<endl<<rlatt<<endl
          <<"gamma  g=g+basis["<<ib<<"]: "<<mG(1)<<endl;
        //bring inside radius
        irlatt=inverse(lattice);
        mG(2,1)=0.0; mG(2,2)=0.0; mG(2,3)=0.0;
        for(it=1;it<4;it++) mG(2,1)=mG(2,1) + mG(1,it)*irlatt(it,1);
        for(it=1;it<4;it++) mG(2,2)=mG(2,2) + mG(1,it)*irlatt(it,2);
        for(it=1;it<4;it++) mG(2,3)=mG(2,3) + mG(1,it)*irlatt(it,3);
        n1=(int)(round(mG(2,1))); n2=(int)(round(mG(2,2))); n3=(int)(round(mG(2,3)));
        flog<<"Gamma cartesian  coordinate: "<<mG(1)<<endl
          <<"Gamma fractional coordinate: "<<mG(2)<<endl
          <<"shifting along -(n1*a1+n2*a2+n3*a3) with n1 n2 n3: "<<n1<<" "<<n2<<" "<<n3<<endl;
        for(it=1;it<4;it++) rlatt(1,it)=rlatt(1,it) - (n1*a1(it)+n2*a2(it)+n3*a3(it));
        for(it=1;it<4;it++) rlatt(2,it)=rlatt(2,it) - (n1*a1(it)+n2*a2(it)+n3*a3(it));
        for(it=1;it<4;it++) rlatt(3,it)=rlatt(3,it) - (n1*a1(it)+n2*a2(it)+n3*a3(it));
        for(it=1;it<4;it++) mG(1,it)=mG(1,it) - (n1*a1(it)+n2*a2(it)+n3*a3(it));
        if(modulus(rlatt(1))>(radius+3*_eps_abc_) ||
            modulus(rlatt(2))>(radius+3*_eps_abc_) ||
            modulus(rlatt(3))>(radius+3*_eps_abc_)
          ) {
          flog<<"ERROR: one of the rlatt points is outside radius, algorithm 1.5*Rsphere failed."<<endl;
          flog.close();
          //deallocate pointer
          for(uint i=0;i<grid_clattice.size();i++)
            delete grid_clattice[i];
          grid_clattice.clear();
          for(uint i=0;i<grid_flattice.size();i++)
            delete grid_flattice[i];
          grid_flattice.clear();
          message << "One of the rlatt points is outside radius, algorithm 1.5*Rsphere failed.";
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_RANGE_);
        }

        flog<<"rlatt row-lattice vectors new coord: "<<endl<<rlatt<<endl
          <<"Gamma-point new cart coord: "<<mG(1)<<endl;
        //there is a chance that one of the rlatt rows is zero
        //or rlatt rows are coplanar
        //check mapping
        sym_found=false;
        foundG=false;founda1=false;founda2=false;founda3=false;
        for(i=0;i<grid_clattice.size();i++) {
          if( !foundG && (modulus(mG(1))<_eps_abc_) ) {
            foundG=true;flog<<"*foundG*"<<endl;}
          if( !foundG && (abs(mG(1,1)-(*grid_clattice[i])(1)) < _eps_abc_) &&
              (abs(mG(1,2)-(*grid_clattice[i])(2)) < _eps_abc_) &&
              (abs(mG(1,3)-(*grid_clattice[i])(3)) < _eps_abc_) ) {
            foundG=true;flog<<"**foundG**"<<endl;}
          if( !founda1 && (modulus(rlatt(1))<_eps_abc_) ) {
            founda1=true;flog<<"*founda1*"<<endl;}
          if( !founda1 && (abs(rlatt(1,1)-(*grid_clattice[i])(1)) < _eps_abc_) &&
              (abs(rlatt(1,2)-(*grid_clattice[i])(2)) < _eps_abc_) &&
              (abs(rlatt(1,3)-(*grid_clattice[i])(3)) < _eps_abc_) ) {
            founda1=true;flog<<"**founda1**"<<endl;}
          if( !founda2 && (modulus(rlatt(2))<_eps_abc_) ) {
            founda2=true;flog<<"*founda2*"<<endl;}
          if( !founda2 && (abs(rlatt(2,1)-(*grid_clattice[i])(1)) < _eps_abc_) &&
              (abs(rlatt(2,2)-(*grid_clattice[i])(2)) < _eps_abc_) &&
              (abs(rlatt(2,3)-(*grid_clattice[i])(3)) < _eps_abc_) ) {
            founda2=true;flog<<"**founda2**"<<endl;}
          if( !founda3 && (modulus(rlatt(3))<_eps_abc_) ) {
            founda3=true;flog<<"*founda3*"<<endl;}
          if( !founda3 && (abs(rlatt(3,1)-(*grid_clattice[i])(1)) < _eps_abc_) &&
              (abs(rlatt(3,2)-(*grid_clattice[i])(2)) < _eps_abc_) &&
              (abs(rlatt(3,3)-(*grid_clattice[i])(3)) < _eps_abc_) ) {
            founda3=true;flog<<"**founda3**"<<endl;}
          if(foundG && founda1 && founda2 && founda3) {
            sym_found=true;
            flog<<"***sym_found***"<<endl;
            Bsitesym[ib].push_back(jp); break;
          }
        }

      }//loop for pgroup
    }//loop for basis

    for(i=0;i<Nbasis;i++) {
      cout<<"Site symmetry operations for basis["<<i<<"]:"<<endl;
      for(j=0;j<Bsitesym[i].size();j++) {
        cout<<"sitesym_op "<<j<<"/"<<Bsitesym[i].size()<<": "
          <<a.pgroup[Bsitesym[i][j]].str_type<<", "
          <<a.pgroup[Bsitesym[i][j]].angle<<" about "<<a.pgroup[Bsitesym[i][j]].axis<<endl
          <<a.pgroup[Bsitesym[i][j]].Uc<<endl;
      }
    }

    uint iA,iB,jpA,jpB;
    vector<vector<vector<uint> > > AMBsitesym(Nbasis-1);
    xvector<double> AMB(3);//cartesian coord on middle point of A&B
    for(iA=0;iA<Nbasis-1;iA++) {
      AMBsitesym[iA].resize(Nbasis);
      for(iB=0;iB<Nbasis;iB++)  AMBsitesym[iA][iB].clear();
    }
    if(ComMidss) {
      for(iA=0;iA<Nbasis-1;iA++) {
        for(jpA=0;jpA<Bsitesym[iA].size();jpA++) {
          for(iB=iA+1;iB<Nbasis;iB++) {
            for(jpB=0;jpB<Bsitesym[iB].size();jpB++) {
              if(Bsitesym[iA][jpA]==Bsitesym[iB][jpB]) {
                jp=Bsitesym[iB][jpB];//index of symop in pgroup that is common for site A and B.
                cerr<<"common site symop found: "<<jp<<endl;
                for(it=1;it<4;it++) AMB[it]=(basis[iA][it]+basis[iB][it])/2.0;
                //r=r-A
                for(it=1;it<4;it++) rlattb4Uc(1,it)=lattice(1,it)-AMB[it];
                for(it=1;it<4;it++) rlattb4Uc(2,it)=lattice(2,it)-AMB[it];
                for(it=1;it<4;it++) rlattb4Uc(3,it)=lattice(3,it)-AMB[it];
                for(it=1;it<4;it++) mGb4Uc(1,it)=-AMB[it];
                mGb4Uc(3)=mGb4Uc(2)=mGb4Uc(1);
                //r=Ur
                rlatttmp=a.pgroup[jp].Uc*(trasp(rlattb4Uc));
                rlatt=trasp(rlatttmp);
                rlatttmp=a.pgroup[jp].Uc*(trasp(mGb4Uc));
                mG=trasp(rlatttmp);
                //r=r+A
                for(it=1;it<4;it++) rlatt(1,it)=rlatt(1,it)+AMB[it];
                for(it=1;it<4;it++) rlatt(2,it)=rlatt(2,it)+AMB[it];
                for(it=1;it<4;it++) rlatt(3,it)=rlatt(3,it)+AMB[it];
                for(it=1;it<4;it++) mG(1,it)=mG(1,it)+AMB[it];
                //bring inside radius
                irlatt=inverse(lattice);
                mG(2,1)=0.0; mG(2,2)=0.0; mG(2,3)=0.0;
                for(it=1;it<4;it++) mG(2,1)=mG(2,1) + mG(1,it)*irlatt(it,1);
                for(it=1;it<4;it++) mG(2,2)=mG(2,2) + mG(1,it)*irlatt(it,2);
                for(it=1;it<4;it++) mG(2,3)=mG(2,3) + mG(1,it)*irlatt(it,3);
                n1=(int)(round(mG(2,1))); n2=(int)(round(mG(2,2))); n3=(int)(round(mG(2,3)));
                //flog<<"Gamma cartesian  coordinate: "<<mG(1)<<endl
                //<<"Gamma fractional coordinate: "<<mG(2)<<endl
                //<<"shifting along -(n1*a1+n2*a2+n3*a3) with n1 n2 n3: "<<n1<<" "<<n2<<" "<<n3<<endl;
                for(it=1;it<4;it++) rlatt(1,it)=rlatt(1,it) - (n1*a1(it)+n2*a2(it)+n3*a3(it));
                for(it=1;it<4;it++) rlatt(2,it)=rlatt(2,it) - (n1*a1(it)+n2*a2(it)+n3*a3(it));
                for(it=1;it<4;it++) rlatt(3,it)=rlatt(3,it) - (n1*a1(it)+n2*a2(it)+n3*a3(it));
                for(it=1;it<4;it++) mG(1,it)=mG(1,it) - (n1*a1(it)+n2*a2(it)+n3*a3(it));
                if(modulus(rlatt(1))>(radius+3*_eps_abc_) ||
                    modulus(rlatt(2))>(radius+3*_eps_abc_) ||
                    modulus(rlatt(3))>(radius+3*_eps_abc_)
                  ) {
                  flog<<"ERROR mid: one of the rlatt points is outside radius, algorithm 1.5*Rsphere failed."<<endl;
                  flog.close();
                  //deallocate pointer
                  for(uint i=0;i<grid_clattice.size();i++)
                    delete grid_clattice[i];
                  grid_clattice.clear();
                  for(uint i=0;i<grid_flattice.size();i++)
                    delete grid_flattice[i];
                  grid_flattice.clear();
                  message << "Mid: One of the rlatt points is outside radius, algorithm 1.5*Rsphere failed.";
                  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_RANGE_);
                }
                //check mapping
                sym_found=false;
                foundG=false;founda1=false;founda2=false;founda3=false;
                for(i=0;i<grid_clattice.size();i++) {
                  if( !foundG && (modulus(mG(1))<_eps_abc_) ) {
                    foundG=true;}
                  if( !foundG && (abs(mG(1,1)-(*grid_clattice[i])(1)) < _eps_abc_) &&
                      (abs(mG(1,2)-(*grid_clattice[i])(2)) < _eps_abc_) &&
                      (abs(mG(1,3)-(*grid_clattice[i])(3)) < _eps_abc_) ) {
                    foundG=true;}
                  if( !founda1 && (modulus(rlatt(1))<_eps_abc_) ) {
                    founda1=true;}
                  if( !founda1 && (abs(rlatt(1,1)-(*grid_clattice[i])(1)) < _eps_abc_) &&
                      (abs(rlatt(1,2)-(*grid_clattice[i])(2)) < _eps_abc_) &&
                      (abs(rlatt(1,3)-(*grid_clattice[i])(3)) < _eps_abc_) ) {
                    founda1=true;}
                  if( !founda2 && (modulus(rlatt(2))<_eps_abc_) ) {
                    founda2=true;}
                  if( !founda2 && (abs(rlatt(2,1)-(*grid_clattice[i])(1)) < _eps_abc_) &&
                      (abs(rlatt(2,2)-(*grid_clattice[i])(2)) < _eps_abc_) &&
                      (abs(rlatt(2,3)-(*grid_clattice[i])(3)) < _eps_abc_) ) {
                    founda2=true;}
                  if( !founda3 && (modulus(rlatt(3))<_eps_abc_) ) {
                    founda3=true;}
                  if( !founda3 && (abs(rlatt(3,1)-(*grid_clattice[i])(1)) < _eps_abc_) &&
                      (abs(rlatt(3,2)-(*grid_clattice[i])(2)) < _eps_abc_) &&
                      (abs(rlatt(3,3)-(*grid_clattice[i])(3)) < _eps_abc_) ) {
                    founda3=true;}
                  if(foundG && founda1 && founda2 && founda3) {
                    sym_found=true;
                    AMBsitesym[iA][iB].push_back(jp); break;
                  }
                }

              }//find common A&B
            }//jpB
          }//iB
        }//jpA
      }//iA
    } //ComMidss

    for(i=0;i<Nbasis;i++) {
      cout<<"Index of site symop for basis["<<i<<"]:(total "<<Bsitesym[i].size()<<")";
      for(j=0;j<Bsitesym[i].size();j++) {
        cout<<" "<<Bsitesym[i][j];
      }
      cout<<endl;
    }

    for(iA=0;iA<Nbasis-1;iA++) {
      for(iB=iA+1;iB<Nbasis;iB++) {
        cout<<"index common site sym_op "<<iA<<"-mid-"<<iB<<":(total "<<AMBsitesym[iA][iB].size()<<")";
        for(jp=0;jp<AMBsitesym[iA][iB].size();jp++) cout<<" "<<AMBsitesym[iA][iB][jp];
        cout<<endl;
      }
    }

    flog.close();

    // ------------------------------------------------------------------------------
    // going out, destroy useless things
    for(uint i=0;i<grid_clattice.size();i++)
      delete grid_clattice[i];
    grid_clattice.clear();
    for(uint i=0;i<grid_flattice.size();i++)
      delete grid_flattice[i];
    grid_flattice.clear();

  }
} // namespace SYM

//DX20170802 START: Xgroups to JSON
// --------------------------------------------------------------------------
// ------------------------------------------------------------- WRITE GROUPS TO JSON
// SymmetryToJson
string SymmetryToJson(vector<_sym_op>& group, char& mode){
  string eendl="";
  bool roff=true; //round off
  stringstream sss;
  stringstream sscontent_json;
  vector<string> vcontent_json;

  string group_str = "";
  if(mode==_PGROUP_){ group_str = "pgroup"; }
  if(mode==_PGROUPK_){ group_str = "pgroupk"; }
  if(mode==_FGROUP_){ group_str = "fgroup"; }
  if(mode==_SGROUP_){ group_str = "sgroup"; }
  if(mode==_AGROUP_){ group_str = "agroup"; }
  if(mode==_PGROUP_XTAL_){ group_str = "pgroup_xtal"; }
  if(mode==_PGROUPK_XTAL_){ group_str = "pgroupk_xtal"; }
  if(mode==_PGROUPK_PATTERSON_){ group_str = "pgroupk_Patterson"; }

  sss << "[" << eendl;

  for(uint i=0;i<group.size();i++){
    // group
    if(group_str.size()){
      sscontent_json << "\"group\":\"" << group_str << "\"" << eendl;
    } else {
      if(PRINT_NULL_JSON){ sscontent_json << "\"group\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // site (agroup only)
    if(group[i].ctau.rows && group_str == "agroup"){
      sscontent_json << "\"site\":" << group[i].site << eendl;
    } else if (group_str == "agroup"){
      if(PRINT_NULL_JSON){ sscontent_json << "\"site\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // type
    if(group[i].str_type.size()){
      sscontent_json << "\"type\":\"" << aurostd::RemoveWhiteSpacesFromTheBack(group[i].str_type) << "\"" << eendl;
    } else {
      if(PRINT_NULL_JSON){ sscontent_json << "\"type\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // Hermann-Mauguin
    if(group[i].str_Hermann_Mauguin.size()){
      sscontent_json << "\"Hermann_Mauguin\":\"" << group[i].str_Hermann_Mauguin << "\"" << eendl;
    } else {
      if(PRINT_NULL_JSON){ sscontent_json << "\"Hermann_Mauguin\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // Schoenflies
    if(group[i].str_Schoenflies.size()){
      sscontent_json << "\"Schoenflies\":\"" << group[i].str_Schoenflies << "\"" << eendl;
    } else {
      if(PRINT_NULL_JSON){ sscontent_json << "\"Schoenflies\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // Uc
    if(group[i].Uc.lrows){
      sscontent_json << "\"Uc\":[" << aurostd::xmatDouble2String(group[i].Uc,5,roff) << "]" << eendl;
    } else {
      if(PRINT_NULL_JSON){ sscontent_json << "\"Uc\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // Uf
    if(group[i].Uf.lrows){
      sscontent_json << "\"Uf\":[" << aurostd::xmatDouble2String(group[i].Uf,1,roff) << "]" << eendl;
    } else {
      if(PRINT_NULL_JSON){ sscontent_json << "\"Uf\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // generator 
    if(group[i].generator.lrows){
      sscontent_json << "\"generator\":[" << aurostd::xmatDouble2String(group[i].generator,5,roff) << "]" << eendl;
    } else {
      if(PRINT_NULL_JSON){ sscontent_json << "\"generator\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    //DX20171207 - added generator_coefficients - START
    // generator coefficients
    if(group[i].generator.lrows){
      sscontent_json << "\"generator_coefficients\":[" << aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(group[i].generator_coefficients,5,roff),",") << "]" << eendl; //DX20180726 - added roff
    } else {
      if(PRINT_NULL_JSON){ sscontent_json << "\"generator_coefficients\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    //DX20171207 - added generator_coefficients - END

    //DX20180117 - added SU2_matrix - START
    // SU(2) matrix
    if(group[i].SU2_matrix.lrows){
      sscontent_json << "\"SU2_matrix\":[" << "[" << aurostd::xcomplex2json(group[i].SU2_matrix(1,1)) << "," << aurostd::xcomplex2json(group[i].SU2_matrix(1,2)) << "]" << "," << eendl; 
      sscontent_json << "[" << aurostd::xcomplex2json(group[i].SU2_matrix(2,1)) << "," << aurostd::xcomplex2json(group[i].SU2_matrix(2,2)) << "]" << "]" << eendl; 
    } else {
      if(PRINT_NULL_JSON){ sscontent_json << "\"SU2_matrix\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    //DX20180117 - added SU2_matrix - END

    //DX20180117 - added su2_coefficients - START
    // su(2) coefficients
    if(group[i].su2_coefficients.lrows){
      sscontent_json << "\"su2_coefficients\":[" <<  aurostd::xcomplex2json(group[i].su2_coefficients(1)) << "," << aurostd::xcomplex2json(group[i].su2_coefficients(2)) << "," << aurostd::xcomplex2json(group[i].su2_coefficients(3)) << "]" << eendl; 
    } else {
      if(PRINT_NULL_JSON){ sscontent_json << "\"su2_coefficients\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    //DX20180117 - added su2_coefficients - END

    // angle
    if(group[i].angle!=AUROSTD_NAN){
      sscontent_json << "\"angle\":" << group[i].angle << eendl;
    } else {
      if(PRINT_NULL_JSON){ sscontent_json << "\"angle\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // axis
    if(group[i].axis.lrows){
      sscontent_json << "\"axis\":[" << aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(group[i].axis,5,roff),",") << "]" << eendl; //DX20180726 - added roff
    } else {
      if(PRINT_NULL_JSON){ sscontent_json << "\"axis\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // quaternion vector
    if(group[i].quaternion_vector.lrows){
      sscontent_json << "\"quaternion_vector\":[" << aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(group[i].quaternion_vector,5,roff),",") << "]" << eendl;
    } else {
      if(PRINT_NULL_JSON){ sscontent_json << "\"quaternion_vector\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // quaternion matrix
    if(group[i].quaternion_matrix.lrows){
      sscontent_json << "\"quaternion_matrix\":[" << aurostd::xmatDouble2String(group[i].quaternion_matrix,5,roff) << "]" << eendl;
    } else {
      if(PRINT_NULL_JSON){ sscontent_json << "\"quaternion_matrix\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // inversion
    if(group[i].flag_inversion == true || group[i].flag_inversion == false){
      sscontent_json << "\"inversion\":";
      if(group[i].flag_inversion){
        sscontent_json << "true" << eendl;
      }
      else if(!group[i].flag_inversion){
        sscontent_json << "false" << eendl;
      }
    } else {
      if(PRINT_NULL_JSON){ sscontent_json << "\"inversion\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // ctau
    if(group[i].ctau.rows && (group_str == "fgroup" || group_str == "sgroup")){
      sscontent_json << "\"ctau\":[" << aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(group[i].ctau,5,roff),",") << "]" << eendl;
    } else if (group_str == "fgroup" || group_str == "sgroup"){
      if(PRINT_NULL_JSON){ sscontent_json << "\"ctau\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // ftau
    if(group[i].ftau.rows && (group_str == "fgroup" || group_str == "sgroup")){
      sscontent_json << "\"ftau\":[" << aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(group[i].ftau,5,roff),",") << "]" << eendl;
    } else if (group_str == "fgroup" || group_str == "sgroup"){
      if(PRINT_NULL_JSON){ sscontent_json << "\"ftau\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // basis_atoms_map
    if(group[i].basis_atoms_map.size() && (group_str == "fgroup" || group_str == "sgroup")){
      sscontent_json << "\"basis_atoms_map\":[" << aurostd::joinWDelimiter(group[i].basis_atoms_map,",") << "]" << eendl;
    } else if (group_str == "fgroup" || group_str == "sgroup"){
      if(PRINT_NULL_JSON){ sscontent_json << "\"basis_atoms_map\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // basis_types_map
    if(group[i].basis_types_map.size() && (group_str == "fgroup" || group_str == "sgroup")){
      sscontent_json << "\"basis_types_map\":[" << aurostd::joinWDelimiter(group[i].basis_types_map,",") << "]" << eendl;
    } else if (group_str == "fgroup" || group_str == "sgroup"){
      if(PRINT_NULL_JSON){ sscontent_json << "\"basis_types_map\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // ctrasl
    if(group[i].ctrasl.rows && group_str == "sgroup"){
      sscontent_json << "\"ctrasl\":[" << aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(group[i].ctrasl,5,roff),",") << "]" << eendl;
    } else if (group_str == "sgroup"){
      if(PRINT_NULL_JSON){ sscontent_json << "\"ctrasl\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // ftrasl
    if(group[i].ftrasl.rows && group_str == "sgroup"){
      sscontent_json << "\"ftrasl\":[" << aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(group[i].ftrasl,5,roff),",") << "]" << eendl;
    } else if (group_str == "sgroup"){
      if(PRINT_NULL_JSON){ sscontent_json << "\"ftrasl\":null" << eendl;}
    }
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // Put into json symmetry operation object
    sss << "{" << aurostd::joinWDelimiter(vcontent_json,",")  << "}";
    vcontent_json.clear();
    if(i != group.size()-1){
      sss << ",";
    }
  }
  sss << "]" << eendl;
  return sss.str();
}
//DX20170802 END: Xgroups to JSON

//DX20170802 START: agroups to JSON
// --------------------------------------------------------------------------
// ------------------------------------------------------------- WRITE AROUPS TO JSON
string AgroupSymmetryToJson(vector<vector<_sym_op> >& group, char& mode){
  string eendl="";
  stringstream sss;
  stringstream sscontent_json;
  vector<string> vcontent_json;

  sss << "[" << eendl;
  for(uint i=0;i<group.size();i++){
    sscontent_json << SymmetryToJson(group[i], mode) << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
  }
  sss << aurostd::joinWDelimiter(vcontent_json,",") << eendl;
  vcontent_json.clear();
  sss << "]" << eendl;
  return sss.str();
}
//DX20170802 END: agroups to JSON

//DX20170802 START: Equivalent atoms to JSON
// --------------------------------------------------------------------------
// ------------------------------------------------------------- WRITE EQUIVALENT ATOMS TO JSON
string EquivalentAtomsToJson(vector<vector<int> >& iatoms){
  string eendl="";
  stringstream sss;
  stringstream sscontent_json;
  vector<string> vcontent_json;

  // iatoms
  if(iatoms.size()){
    sscontent_json << "\"inequivalent_atoms\":[" << eendl;
    for(uint iat1=0;iat1<iatoms.size();iat1++) {
      sscontent_json << iatoms[iat1][0] << eendl;
      if(iat1 != iatoms.size()-1){
        sscontent_json << "," << eendl;
      }
    }
    sscontent_json << "]" << eendl;
  } else {
    if(PRINT_NULL_JSON){ sscontent_json << "\"inequivalent_atoms\":null" << eendl;}
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // equivalent atom sets
  if(iatoms.size()){
    sscontent_json << "\"equivalent_sets\":[" << eendl;
    for(uint iat1=0;iat1<iatoms.size();iat1++) {
      sscontent_json << "[" << eendl;
      for(uint iat2=0;iat2<iatoms[iat1].size();iat2++){
        sscontent_json << iatoms[iat1][iat2] << eendl;
        if(iat2 != iatoms[iat1].size()-1){
          sscontent_json << "," << eendl;
        }
      }
      sscontent_json << "]" << eendl;
      if(iat1 != iatoms.size()-1){
        sscontent_json << "," << eendl;
      }
    }
    sscontent_json << "]" << eendl;
  } else {
    if(PRINT_NULL_JSON){ sscontent_json << "\"equivalent_sets\":null" << eendl;}
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  sss << "{" << aurostd::joinWDelimiter(vcontent_json,",")  << "}";
  vcontent_json.clear();
  return sss.str();
}
//DX20170802 END: Equivalent atoms to JSON

// --------------------------------------------------------------------------
// ------------------------------------------------------------- WRITE GROUPS
// Function KBIN_SymmetryWrite
//
// This function writes aflow.pgroup.out, aflow.fgroup.out and
// aflow.sgroup.out files on the directory - SC aug 2007

//CO20171025 - redundant
////DX20170802 START: Adding symmetry output formatting option
//bool KBIN_SymmetryWrite(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,char mode,const bool& osswrite,ostream& oss) {
//  string format = "txt"; //Default output is aflow.xgroup.out (a text output)
//  return KBIN_SymmetryWrite(FileMESSAGE,a,aflags,mode,osswrite,oss,format);
//}
////DX20170802 END: Adding symmetry output formatting option

bool KBIN_SymmetryWrite(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,char mode,const bool& osswrite,ostream& oss,const string& format) { //DX20170802
  stringstream message;
  ostringstream aus;
  xvector<double> aux_rrr(9),aux_ijk(9);
  bool Krun=TRUE;
  ofstream FileOUTPUT;string FileNameOUTPUT;
  // ------------------------------------------------------------------------------
  // writing aflow.pgroup
  if(mode==_PGROUP_) aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUP Symmetry: writing BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
  if(mode==_PGROUPK_) aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUPK Symmetry: writing BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
  if(mode==_FGROUP_) aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: writing BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
  if(mode==_SGROUP_) aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "SGROUP Symmetry: writing BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
  if(mode==_AGROUP_) aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP Symmetry: writing BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
  if(mode==_IATOMS_) aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "IATOMS Symmetry: writing BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
  if(mode==_PGROUP_XTAL_) aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUP_XTAL Symmetry: writing BEGIN " << Message(__AFLOW_FILE__,aflags) << endl;
  if(mode==_PGROUPK_XTAL_) aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUPK_XTAL Symmetry: writing BEGIN " << Message(__AFLOW_FILE__,aflags) << endl; //DX20171205 - Added pgroupk_xtal
  if(mode==_PGROUPK_PATTERSON_) aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUPK_PATTERSON Symmetry: writing BEGIN " << Message(__AFLOW_FILE__,aflags) << endl; //DX20200129
  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
  if(aurostd::toupper(format)=="TXT" || aurostd::toupper(format)=="TEXT"){ //DX20200206
    if(mode==_PGROUP_) FileNameOUTPUT=aflags.Directory+"/"+DEFAULT_AFLOW_PGROUP_OUT;
    if(mode==_PGROUPK_) FileNameOUTPUT=aflags.Directory+"/"+DEFAULT_AFLOW_PGROUPK_OUT;
    if(mode==_FGROUP_) FileNameOUTPUT=aflags.Directory+"/"+DEFAULT_AFLOW_FGROUP_OUT;
    if(mode==_SGROUP_) FileNameOUTPUT=aflags.Directory+"/"+DEFAULT_AFLOW_SGROUP_OUT;
    if(mode==_AGROUP_) FileNameOUTPUT=aflags.Directory+"/"+DEFAULT_AFLOW_AGROUP_OUT;
    if(mode==_IATOMS_) FileNameOUTPUT=aflags.Directory+"/"+DEFAULT_AFLOW_IATOMS_OUT;
    if(mode==_PGROUP_XTAL_) FileNameOUTPUT=aflags.Directory+"/"+DEFAULT_AFLOW_PGROUP_XTAL_OUT;
    if(mode==_PGROUPK_XTAL_) FileNameOUTPUT=aflags.Directory+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_OUT; //DX20171205 - Added pgroupk_xtal
    if(mode==_PGROUPK_PATTERSON_) FileNameOUTPUT=aflags.Directory+"/"+DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT; //DX20200129
  }
  if(aurostd::toupper(format)=="JSON"){ //DX20200206
    if(mode==_PGROUP_) FileNameOUTPUT=aflags.Directory+"/"+DEFAULT_AFLOW_PGROUP_JSON;
    if(mode==_PGROUPK_) FileNameOUTPUT=aflags.Directory+"/"+DEFAULT_AFLOW_PGROUPK_JSON;
    if(mode==_FGROUP_) FileNameOUTPUT=aflags.Directory+"/"+DEFAULT_AFLOW_FGROUP_JSON;
    if(mode==_SGROUP_) FileNameOUTPUT=aflags.Directory+"/"+DEFAULT_AFLOW_SGROUP_JSON;
    if(mode==_AGROUP_) FileNameOUTPUT=aflags.Directory+"/"+DEFAULT_AFLOW_AGROUP_JSON;
    if(mode==_IATOMS_) FileNameOUTPUT=aflags.Directory+"/"+DEFAULT_AFLOW_IATOMS_JSON;
    if(mode==_PGROUP_XTAL_) FileNameOUTPUT=aflags.Directory+"/"+DEFAULT_AFLOW_PGROUP_XTAL_JSON;
    if(mode==_PGROUPK_XTAL_) FileNameOUTPUT=aflags.Directory+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_JSON; //DX20171205 - Added pgroupk_xtal
    if(mode==_PGROUPK_PATTERSON_) FileNameOUTPUT=aflags.Directory+"/"+DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON; //DX20200129
  }
  FileOUTPUT.open(FileNameOUTPUT.c_str(),std::ios::out);
  if(aurostd::toupper(format)=="TXT" || aurostd::toupper(format)=="TEXT"){ //DX20200206
    FileOUTPUT << SEPARATION_LINE_DASH << endl;
  } //DX20170802
  if(mode==_PGROUP_) {
    if(a.pgroup.empty()){ //DX20210327 - check if empty
      message << "No PGROUP (lattice point group) operations! Symmetry calculation failed (bug).";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
    }
    if(aurostd::toupper(format)=="TXT" || aurostd::toupper(format)=="TEXT"){ //DX20200206
      FileOUTPUT << "AFLOW point group file, operations are as a=U*b (cols vectors), (Uc/Uf for cartesian/fractional) " << endl;
      FileOUTPUT << a.pgroup.size() << "    point group operations " << endl;
      for(uint k=0;k<a.pgroup.size();k++) {
        for(int i=0;i<9;i++) aux_rrr(i+1)=a.pgroup[k].Uc(int(i/3)+1,mod(i,3)+1);                 // put in rows
        for(int i=0;i<9;i++) aux_ijk(i+1)=a.pgroup[k].Uf(int(i/3)+1,mod(i,3)+1);                  // put in rows
        if(k==0) FileOUTPUT << SEPARATION_LINE_DASH << endl;
        FileOUTPUT << " Operation number = " << k+1 << " / " << a.pgroup.size() << endl;
        FileOUTPUT << a.pgroup[k]; // << endl;
        FileOUTPUT << SEPARATION_LINE_DASH << endl;
        FileOUTPUT.flush();
      }
    }                                                       //DX20170802
    if(aurostd::toupper(format)=="JSON"){ //DX20200206
      FileOUTPUT << SymmetryToJson(a.pgroup,mode);          //DX20170802
      FileOUTPUT.flush();
    }                                                       //DX20170802
    FileOUTPUT.clear();FileOUTPUT.close();
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUP Symmetry: writing END " << Message(__AFLOW_FILE__,aflags) << endl;
  }
  if(mode==_PGROUP_XTAL_) {
    if(a.pgroup_xtal.empty()){ //DX20210327 - check if empty
      message << "No PGROUP_XTAL (crystal point group) operations! Symmetry calculation failed (bug).";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
    }
    if(aurostd::toupper(format)=="TXT" || aurostd::toupper(format)=="TEXT"){ //DX20200206
      FileOUTPUT << "AFLOW crystal point group file, operations are as a=U*b (cols vectors), (Uc/Uf for cartesian/fractional) " << endl;
      FileOUTPUT << a.pgroup_xtal.size() << "    crystal point group operations " << endl;
      for(uint k=0;k<a.pgroup_xtal.size();k++) {
        for(int i=0;i<9;i++) aux_rrr(i+1)=a.pgroup_xtal[k].Uc(int(i/3)+1,mod(i,3)+1);                 // put in rows
        for(int i=0;i<9;i++) aux_ijk(i+1)=a.pgroup_xtal[k].Uf(int(i/3)+1,mod(i,3)+1);                  // put in rows
        if(k==0) FileOUTPUT << SEPARATION_LINE_DASH << endl;
        FileOUTPUT << " Operation number = " << k+1 << " / " << a.pgroup_xtal.size() << endl;
        FileOUTPUT << a.pgroup_xtal[k]; // << endl;
        FileOUTPUT << SEPARATION_LINE_DASH << endl;
        FileOUTPUT.flush();
      }
    }                                                       //DX20170802
    if(aurostd::toupper(format)=="JSON"){ //DX20200206
      FileOUTPUT << SymmetryToJson(a.pgroup_xtal,mode);     //DX20170802
      FileOUTPUT.flush();                                   //DX20170802
    }                                                       //DX20170802
    FileOUTPUT.clear();FileOUTPUT.close();
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUP_XTAL Symmetry: writing END " << Message(__AFLOW_FILE__,aflags) << endl;
  }
  //DX20171205 - Added pgroupk_xtal - START
  if(mode==_PGROUPK_XTAL_) {
    if(a.pgroupk_xtal.empty()){ //DX20210327 - check if empty
      message << "No PGROUPK_XTAL (dual of crystal point group) operations! Symmetry calculation failed (bug).";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
    }
    if(aurostd::toupper(format)=="TXT" || aurostd::toupper(format)=="TEXT"){ //DX20200206
      FileOUTPUT << "AFLOW dual of crystal point group file, operations are as a=U*b (cols vectors), (Uc/Uf for cartesian/fractional) " << endl;
      FileOUTPUT << a.pgroupk_xtal.size() << "    dual of crystal point group operations " << endl;
      for(uint k=0;k<a.pgroupk_xtal.size();k++) {
        for(int i=0;i<9;i++) aux_rrr(i+1)=a.pgroupk_xtal[k].Uc(int(i/3)+1,mod(i,3)+1);                 // put in rows
        for(int i=0;i<9;i++) aux_ijk(i+1)=a.pgroupk_xtal[k].Uf(int(i/3)+1,mod(i,3)+1);                  // put in rows
        if(k==0) FileOUTPUT << SEPARATION_LINE_DASH << endl;
        FileOUTPUT << " Operation number = " << k+1 << " / " << a.pgroupk_xtal.size() << endl;
        FileOUTPUT << a.pgroupk_xtal[k]; // << endl;
        FileOUTPUT << SEPARATION_LINE_DASH << endl;
        FileOUTPUT.flush();
      }
    }                                                       //DX20170802
    if(aurostd::toupper(format)=="JSON"){ //DX20200206
      FileOUTPUT << SymmetryToJson(a.pgroupk_xtal,mode);    //DX20170802
      FileOUTPUT.flush();                                   //DX20170802
    }                                                       //DX20170802
    FileOUTPUT.clear();FileOUTPUT.close();
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUPK_XTAL Symmetry: writing END " << Message(__AFLOW_FILE__,aflags) << endl;
  }
  //DX20171205 - Added pgroupk_xtal - END
  if(mode==_PGROUPK_) {
    if(a.pgroupk.empty()){ //DX20210327 - check if empty
      message << "No PGROUPK (dual of lattice point group) operations! Symmetry calculation failed (bug).";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
    }
    if(aurostd::toupper(format)=="TXT" || aurostd::toupper(format)=="TEXT"){ //DX20200206
      FileOUTPUT << "AFLOW point group klattice file, operations are as a=U*b (cols vectors), (Uc/Uf for cartesian/fractional) " << endl;
      FileOUTPUT << a.pgroupk.size() << "    point group operations " << endl;
      for(uint k=0;k<a.pgroupk.size();k++) {
        for(int i=0;i<9;i++) aux_rrr(i+1)=a.pgroupk[k].Uc(int(i/3)+1,mod(i,3)+1);                 // put in rows
        for(int i=0;i<9;i++) aux_ijk(i+1)=a.pgroupk[k].Uf(int(i/3)+1,mod(i,3)+1);                  // put in rows
        if(k==0) FileOUTPUT << SEPARATION_LINE_DASH << endl;
        FileOUTPUT << " Operation number = " << k+1 << " / " << a.pgroupk.size() << endl;
        FileOUTPUT << a.pgroupk[k]; // << endl;
        FileOUTPUT << SEPARATION_LINE_DASH << endl;
        FileOUTPUT.flush();
      }
    }                                                       //DX20170802
    if(aurostd::toupper(format)=="JSON"){ //DX20200206
      FileOUTPUT << SymmetryToJson(a.pgroupk,mode);         //DX20170802
      FileOUTPUT.flush();                                   //DX20170802
    }                                                       //DX20170802
    FileOUTPUT.clear();FileOUTPUT.close();
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUPK Symmetry: writing END " << Message(__AFLOW_FILE__,aflags) << endl;
  }
  if(mode==_PGROUPK_PATTERSON_) { //DX20200129
    if(a.pgroupk_Patterson.empty()){ //DX20210327 - check if empty
      message << "No PGROUPK_PATTERSON (Patterson point group) operations! Symmetry calculation failed (bug).";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
    }
    if(aurostd::toupper(format)=="TXT" || aurostd::toupper(format)=="TEXT"){ //DX20200206
      FileOUTPUT << "AFLOW Patterson point group file, operations are as a=U*b (cols vectors), (Uc/Uf for cartesian/fractional) " << endl;
      FileOUTPUT << a.pgroupk_Patterson.size() << "    Patterson point group operations " << endl;
      for(uint k=0;k<a.pgroupk_Patterson.size();k++) {
        for(int i=0;i<9;i++) aux_rrr(i+1)=a.pgroupk_Patterson[k].Uc(int(i/3)+1,mod(i,3)+1);                 // put in rows
        for(int i=0;i<9;i++) aux_ijk(i+1)=a.pgroupk_Patterson[k].Uf(int(i/3)+1,mod(i,3)+1);                  // put in rows
        if(k==0) FileOUTPUT << SEPARATION_LINE_DASH << endl;
        FileOUTPUT << " Operation number = " << k+1 << " / " << a.pgroupk_Patterson.size() << endl;
        FileOUTPUT << a.pgroupk_Patterson[k]; // << endl;
        FileOUTPUT << SEPARATION_LINE_DASH << endl;
        FileOUTPUT.flush();
      }
    }                                                       //DX20170802
    if(aurostd::toupper(format)=="JSON"){ //DX20200206
      FileOUTPUT << SymmetryToJson(a.pgroupk_Patterson,mode);     //DX20170802
      FileOUTPUT.flush();                                   //DX20170802
    }                                                       //DX20170802
    FileOUTPUT.clear();FileOUTPUT.close();
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUPK_PATTERSON Symmetry: writing END " << Message(__AFLOW_FILE__,aflags) << endl;
  }
  if(mode==_FGROUP_) {
    if(a.fgroup.empty()){ //DX20210327 - check if empty
      message << "No FGROUP (factor group representative, unit cell symmetry) operations! Symmetry calculation failed (bug).";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
    }
    if(aurostd::toupper(format)=="TXT" || aurostd::toupper(format)=="TEXT"){ //DX20200206
      FileOUTPUT << "AFLOW factor group file, operations are as a=U*b+tau (cols vectors), (Uc/Uf and ctau/ftau for cartesian/fractional)" << endl;
      FileOUTPUT << a.fgroup.size() << "    factor group operations " << endl;
      for(uint k=0;k<a.fgroup.size();k++) {
        for(int i=0;i<9;i++) aux_rrr(i+1)=a.fgroup[k].Uc(int(i/3)+1,mod(i,3)+1);                       // put in rows
        for(int i=0;i<9;i++) aux_ijk(i+1)=a.fgroup[k].Uf(int(i/3)+1,mod(i,3)+1);                    // put in rows
        if(k==0) FileOUTPUT << SEPARATION_LINE_DASH << endl;
        FileOUTPUT << " Operation number = " << k+1 << " / " << a.fgroup.size() << endl;
        FileOUTPUT << a.fgroup[k]; // << endl;
        FileOUTPUT << SEPARATION_LINE_DASH << endl;
        FileOUTPUT.flush();
      }
    }                                                       //DX20170802
    if(aurostd::toupper(format)=="JSON"){ //DX20200206
      FileOUTPUT << SymmetryToJson(a.fgroup,mode);          //DX20170802
      FileOUTPUT.flush();                                   //DX20170802
    }                                                       //DX20170802
    FileOUTPUT.clear();FileOUTPUT.close();
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: writing END " << Message(__AFLOW_FILE__,aflags) << endl;
  }
  if(mode==_SGROUP_) {
    if(a.sgroup.empty()){ //DX20210327 - check if empty
      message << "No SGROUP (space group) operations! Symmetry calculation failed (bug).";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
    }
    if(aurostd::toupper(format)=="TXT" || aurostd::toupper(format)=="TEXT"){ //DX20200206
      FileOUTPUT << "AFLOW space group file, operations are as a=U*b+tau+trasl (cols vectors), (Uc/Uf,ctau/ftau,ctrasl/ftrasl for cartesian/fractional)" << endl;
      FileOUTPUT << a.sgroup.size() << "    space group operations " << endl;
      FileOUTPUT << a.sgroup_radius << "    radius of space group " << endl;
      FileOUTPUT << a.sgroup_radius_dims << "    dimension of radius of space group " << endl;
      for(uint k=0;k<a.sgroup.size();k++) {
        for(int i=0;i<9;i++) aux_rrr(i+1)=a.sgroup[k].Uc(int(i/3)+1,mod(i,3)+1);                       // put in rows
        for(int i=0;i<9;i++) aux_ijk(i+1)=a.sgroup[k].Uf(int(i/3)+1,mod(i,3)+1);                    // put in rows
        if(k==0) FileOUTPUT << SEPARATION_LINE_DASH << endl;
        FileOUTPUT << " Operation number = " << k+1 << " / " << a.sgroup.size() << endl;
        FileOUTPUT << a.sgroup[k]; // << endl;
        FileOUTPUT << SEPARATION_LINE_DASH << endl;
        FileOUTPUT.flush();
      }
    }                                                       //DX20170802
    if(aurostd::toupper(format)=="JSON"){ //DX20200206
      FileOUTPUT << SymmetryToJson(a.sgroup,mode);          //DX20170802
      FileOUTPUT.flush();                                   //DX20170802
    }                                                       //DX20170802
    FileOUTPUT.clear();FileOUTPUT.close();
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "SGROUP Symmetry: writing END " << Message(__AFLOW_FILE__,aflags) << endl;
  }
  if(mode==_AGROUP_) {
    if(a.agroup.empty()){ //DX20210327 - check if empty
      message << "No AGROUP (site point group) operations! Symmetry calculation failed (bug).";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
    }
    if(aurostd::toupper(format)=="TXT" || aurostd::toupper(format)=="TEXT"){ //DX20200206
      FileOUTPUT << "AFLOW site point group file, operations (centered on the site) are as a=U*b (cols vectors), (Uc/Uf for cartesian/fractional) " << endl;
      for(uint iat=0;iat<a.atoms.size();iat++) {
        FileOUTPUT << " Site="  << iat << endl;
        FileOUTPUT << a.agroup.at(iat).size() << "   site point group operations " << endl;
        for(uint k=0;k<a.agroup.at(iat).size();k++) {
          for(int i=0;i<9;i++) aux_rrr(i+1)=a.agroup.at(iat)[k].Uc(int(i/3)+1,mod(i,3)+1);                       // put in rows
          for(int i=0;i<9;i++) aux_ijk(i+1)=a.agroup.at(iat)[k].Uf(int(i/3)+1,mod(i,3)+1);                    // put in rows
          if(k==0) FileOUTPUT << SEPARATION_LINE_DASH << endl;
          FileOUTPUT << " Site = " << iat << endl;
          FileOUTPUT << " Operation number = " << k+1 << " / " << a.agroup.at(iat).size() << endl;
          FileOUTPUT << a.agroup.at(iat)[k]; // << endl;
          FileOUTPUT << SEPARATION_LINE_DASH << endl;
          FileOUTPUT.flush();
        }
      }
    }                                                       //DX20170802
    if(aurostd::toupper(format)=="JSON"){ //DX20200206
      FileOUTPUT << AgroupSymmetryToJson(a.agroup,mode);    //DX20170802
      FileOUTPUT.flush();                                   //DX20170802
    }                                                       //DX20170802
    FileOUTPUT.clear();FileOUTPUT.close();
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP Symmetry: writing END " << Message(__AFLOW_FILE__,aflags) << endl;
  }
  if(mode==_IATOMS_) {
    if(a.iatoms.empty()){ //DX20210327 - check if empty
      message << "No IATOMS (inequivalent atoms)! Symmetry calculation failed (bug).";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
    }
    if(aurostd::toupper(format)=="TXT" || aurostd::toupper(format)=="TEXT"){ //DX20200206
      FileOUTPUT << "Equivalent atoms file " << endl;
      FileOUTPUT << SEPARATION_LINE_DASH << endl;
      for(uint iat1=0;iat1<a.iatoms.size();iat1++) {
        FileOUTPUT << " [" << a.iatoms.at(iat1).at(0) << "]  ";
        for(uint iat2=0;iat2<a.iatoms.at(iat1).size();iat2++)
          FileOUTPUT << a.iatoms.at(iat1).at(iat2) << " ";
        FileOUTPUT << endl;
      }
      FileOUTPUT << SEPARATION_LINE_DASH << endl;
      bool temp=a.write_inequivalent_flag;
      a.write_inequivalent_flag=TRUE;
      FileOUTPUT << a;
      a.write_inequivalent_flag=temp;
      FileOUTPUT << SEPARATION_LINE_DASH << endl;
      FileOUTPUT.flush();
    }                                                       //DX20170802
    if(aurostd::toupper(format)=="JSON"){ //DX20200206
      FileOUTPUT << EquivalentAtomsToJson(a.iatoms);        //DX20170802
      FileOUTPUT.flush();                                   //DX20170802
    }                                                       //DX20170802
    FileOUTPUT.clear();FileOUTPUT.close();
    aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "IATOMS Symmetry: writing END " << Message(__AFLOW_FILE__,aflags) << endl;

  }
  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
  return Krun;
}

//DX20170803 - Print symmetry to screen - START
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
// ------------------------------------------------------------ PRINT GROUPS TO SCREEN
// Function KBIN_SymmetryToScreen
//
// This function prints to screen all the symmetry elements for a given structure 
bool KBIN_SymmetryToScreen(xstructure& a, const string& format, ostream& oss, char mode){
  stringstream message;

  // OUT format
  if(aurostd::toupper(format)=="TXT" || aurostd::toupper(format)=="TEXT"){ //DX20200206
    //xvector<double> aux_rrr(9),aux_ijk(9);  //OBSOLETE ME20210402 - not used
    if(mode == '\0' || mode == _PGROUP_){
      if(a.pgroup.empty()){ //DX20210327 - check if empty
        message << "No PGROUP (lattice point group) operations! Symmetry calculation failed (bug).";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
      }
      oss << "AFLOW point group, operations are as a=U*b (cols vectors), (Uc/Uf for cartesian/fractional) " << endl;
      oss << a.pgroup.size() << "    point group operations " << endl;
      oss << SEPARATION_LINE_DASH << endl;
      for(uint k=0;k<a.pgroup.size();k++) {
        //for(int i=0;i<9;i++) aux_rrr(i+1)=a.pgroup[k].Uc(int(i/3)+1,mod(i,3)+1);                 // put in rows
        //for(int i=0;i<9;i++) aux_ijk(i+1)=a.pgroup[k].Uf(int(i/3)+1,mod(i,3)+1);                  // put in rows
        oss << " Operation number = " << k+1 << " / " << a.pgroup.size() << endl;
        oss << a.pgroup[k]; // << endl;
        oss << SEPARATION_LINE_DASH << endl;
        oss.flush();
      }
    }
    if(mode == '\0' || mode == _PGROUPK_){
      if(a.pgroupk.empty()){ //DX20210327 - check if empty
        message << "No PGROUPK (dual of lattice point group) operations! Symmetry calculation failed (bug).";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
      }
      oss << "AFLOW klattice point group, operations are as a=U*b (cols vectors), (Uc/Uf for cartesian/fractional) " << endl;
      oss << a.pgroupk.size() << "    point group operations " << endl;
      oss << SEPARATION_LINE_DASH << endl;
      for(uint k=0;k<a.pgroupk.size();k++) {
        //for(int i=0;i<9;i++) aux_rrr(i+1)=a.pgroupk[k].Uc(int(i/3)+1,mod(i,3)+1);                 // put in rows
        //for(int i=0;i<9;i++) aux_ijk(i+1)=a.pgroupk[k].Uf(int(i/3)+1,mod(i,3)+1);                  // put in rows
        oss << " Operation number = " << k+1 << " / " << a.pgroupk.size() << endl;
        oss << a.pgroupk[k]; // << endl;
        oss << SEPARATION_LINE_DASH << endl;
        oss.flush();
      }
    }
    if(mode == '\0' || mode == _FGROUP_){
      if(a.fgroup.empty()){ //DX20210327 - check if empty
        message << "No FGROUP (factor group representative, unit cell symmetry) operations! Symmetry calculation failed (bug).";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
      }
      oss << "AFLOW factor group, operations are as a=U*b+tau (cols vectors), (Uc/Uf and ctau/ftau for cartesian/fractional)" << endl;
      oss << a.fgroup.size() << "    factor group operations " << endl;
      oss << SEPARATION_LINE_DASH << endl;
      for(uint k=0;k<a.fgroup.size();k++) {
        //for(int i=0;i<9;i++) aux_rrr(i+1)=a.fgroup[k].Uc(int(i/3)+1,mod(i,3)+1);                       // put in rows
        //for(int i=0;i<9;i++) aux_ijk(i+1)=a.fgroup[k].Uf(int(i/3)+1,mod(i,3)+1);                    // put in rows
        oss << " Operation number = " << k+1 << " / " << a.fgroup.size() << endl;
        oss << a.fgroup[k]; // << endl;
        oss << SEPARATION_LINE_DASH << endl;
      }
    }
    if(mode == '\0' || mode == _PGROUP_XTAL_){
      if(a.pgroup_xtal.empty()){ //DX20210327 - check if empty
        message << "No PGROUP_XTAL (crystal point group) operations! Symmetry calculation failed (bug).";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
      }
      oss << "AFLOW crystal point group, operations are as a=U*b (cols vectors), (Uc/Uf for cartesian/fractional) " << endl;
      oss << a.pgroup_xtal.size() << "    point group operations " << endl;
      oss << SEPARATION_LINE_DASH << endl;
      for(uint k=0;k<a.pgroup_xtal.size();k++) {
        //for(int i=0;i<9;i++) aux_rrr(i+1)=a.pgroup_xtal[k].Uc(int(i/3)+1,mod(i,3)+1);                 // put in rows
        //for(int i=0;i<9;i++) aux_ijk(i+1)=a.pgroup_xtal[k].Uf(int(i/3)+1,mod(i,3)+1);                  // put in rows
        oss << " Operation number = " << k+1 << " / " << a.pgroup_xtal.size() << endl;
        oss << a.pgroup_xtal[k]; // << endl;
        oss << SEPARATION_LINE_DASH << endl;
        oss.flush();
      }
    }
    //DX20171205 - Added pgroupk_xtal - START
    if(mode == '\0' || mode == _PGROUPK_XTAL_){
      if(a.pgroupk_xtal.empty()){ //DX20210327 - check if empty
        message << "No PGROUPK_XTAL (dual of crystal point group) operations! Symmetry calculation failed (bug).";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
      }
      oss << "AFLOW dual of crystal point group, operations are as a=U*b (cols vectors), (Uc/Uf for cartesian/fractional) " << endl;
      oss << a.pgroupk_xtal.size() << "    dual of crystal point group operations " << endl;
      oss << SEPARATION_LINE_DASH << endl;
      for(uint k=0;k<a.pgroupk_xtal.size();k++) {
        //for(int i=0;i<9;i++) aux_rrr(i+1)=a.pgroupk_xtal[k].Uc(int(i/3)+1,mod(i,3)+1);                 // put in rows
        //for(int i=0;i<9;i++) aux_ijk(i+1)=a.pgroupk_xtal[k].Uf(int(i/3)+1,mod(i,3)+1);                  // put in rows
        oss << " Operation number = " << k+1 << " / " << a.pgroupk_xtal.size() << endl;
        oss << a.pgroupk_xtal[k]; // << endl;
        oss << SEPARATION_LINE_DASH << endl;
        oss.flush();
      }
    }
    //DX20171205 - Added pgroupk_xtal - END
    //DX20200129 - Patterson symmetry - START
    if(mode == '\0' || mode == _PGROUPK_PATTERSON_){
      if(a.pgroupk_Patterson.empty()){ //DX20210327 - check if empty
        message << "No PGROUPK_PATTERSON (Patterson point group) operations! Symmetry calculation failed (bug).";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
      }
      oss << "AFLOW Patterson point group, operation are as a=U*b (cols vectors), (Uc/Uf for cartesian/fractional) " << endl;
      oss << a.pgroupk_Patterson.size() << "    point group operations " << endl;
      oss << SEPARATION_LINE_DASH << endl;
      for(uint k=0;k<a.pgroupk_Patterson.size();k++) {
        //for(int i=0;i<9;i++) aux_rrr(i+1)=a.pgroupk_Patterson[k].Uc(int(i/3)+1,mod(i,3)+1);                 // put in rows
        //for(int i=0;i<9;i++) aux_ijk(i+1)=a.pgroupk_Patterson[k].Uf(int(i/3)+1,mod(i,3)+1);                  // put in rows
        oss << " Operation number = " << k+1 << " / " << a.pgroupk_Patterson.size() << endl;
        oss << a.pgroupk_Patterson[k]; // << endl;
        oss << SEPARATION_LINE_DASH << endl;
        oss.flush();
      }
    }
    //DX20200129 - Patterson symmetry - END
    if(mode == '\0' || mode == _SGROUP_){
      if(a.sgroup.empty()){ //DX20210327 - check if empty
        message << "No SGROUP (space group) operations! Symmetry calculation failed (bug).";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
      }
      oss << "AFLOW space group, operations are as a=U*b+tau+trasl (cols vectors), (Uc/Uf,ctau/ftau,ctrasl/ftrasl for cartesian/fractional)" << endl;
      oss << a.sgroup.size() << "    space group operations " << endl;
      oss << a.sgroup_radius << "    radius of space group " << endl;
      oss << a.sgroup_radius_dims << "    dimension of radius of space group " << endl;
      oss << SEPARATION_LINE_DASH << endl;
      for(uint k=0;k<a.sgroup.size();k++) {
        //for(int i=0;i<9;i++) aux_rrr(i+1)=a.sgroup[k].Uc(int(i/3)+1,mod(i,3)+1);                       // put in rows
        //for(int i=0;i<9;i++) aux_ijk(i+1)=a.sgroup[k].Uf(int(i/3)+1,mod(i,3)+1);                    // put in rows
        oss << " Operation number = " << k+1 << " / " << a.sgroup.size() << endl;
        oss << a.sgroup[k]; // << endl;
        oss << SEPARATION_LINE_DASH << endl;
        oss.flush();
      }
    }
    if(mode == '\0' || mode == _AGROUP_){
      if(a.agroup.empty()){ //DX20210327 - check if empty
        message << "No AGROUP (site point group) operations! Symmetry calculation failed (bug).";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
      }
      oss << "AFLOW site point group, operations (centered on the site) are as a=U*b (cols vectors), (Uc/Uf for cartesian/fractional) " << endl;
      for(uint iat=0;iat<a.atoms.size();iat++) {
        oss << " Site="  << iat << endl;
        oss << a.agroup.at(iat).size() << "   site point group operations " << endl;
        oss << SEPARATION_LINE_DASH << endl;
        for(uint k=0;k<a.agroup.at(iat).size();k++) {
          //for(int i=0;i<9;i++) aux_rrr(i+1)=a.agroup.at(iat)[k].Uc(int(i/3)+1,mod(i,3)+1);                       // put in rows
          //for(int i=0;i<9;i++) aux_ijk(i+1)=a.agroup.at(iat)[k].Uf(int(i/3)+1,mod(i,3)+1);                    // put in rows
          oss << " Site = " << iat << endl;
          oss << " Operation number = " << k+1 << " / " << a.agroup.at(iat).size() << endl;
          oss << a.agroup.at(iat)[k]; // << endl;
          oss << SEPARATION_LINE_DASH << endl;
        }
      }   
    }   
    if(mode == '\0' || mode == _IATOMS_){
      if(a.iatoms.empty()){ //DX20210327 - check if empty
        message << "No IATOMS (inequivalent atoms)! Symmetry calculation failed (bug).";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
      }
      oss << "Equivalent atoms" << endl;
      oss << SEPARATION_LINE_DASH << endl;
      for(uint iat1=0;iat1<a.iatoms.size();iat1++) {
        oss << " [" << a.iatoms.at(iat1).at(0) << "]  ";
        for(uint iat2=0;iat2<a.iatoms.at(iat1).size();iat2++)
          oss << a.iatoms.at(iat1).at(iat2) << " ";
        oss << endl;
      }
      oss << SEPARATION_LINE_DASH << endl;
      bool temp=a.write_inequivalent_flag;
      a.write_inequivalent_flag=TRUE;
      oss << a;
      a.write_inequivalent_flag=temp;
    }
    return TRUE;
  }

  // JSON format
  if(aurostd::toupper(format)=="JSON"){ //DX20200206
    stringstream sscontent_json;
    vector<string> vcontent_json;
    if(mode == '\0' || mode == _PGROUP_){
      char tmp_mode = _PGROUP_;
      sscontent_json << "\"pgroup\":" << SymmetryToJson(a.pgroup,tmp_mode);
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    }
    if(mode == '\0' || mode == _PGROUPK_){
      char tmp_mode = _PGROUPK_;
      sscontent_json << "\"pgroupk\":" << SymmetryToJson(a.pgroupk,tmp_mode);
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    }
    if(mode == '\0' || mode == _FGROUP_){
      char tmp_mode = _FGROUP_;
      sscontent_json << "\"fgroup\":" << SymmetryToJson(a.fgroup,tmp_mode);
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    }
    if(mode == '\0' || mode == _PGROUP_XTAL_){
      char tmp_mode = _PGROUP_XTAL_;
      sscontent_json << "\"pgroup_xtal\":" << SymmetryToJson(a.pgroup_xtal,tmp_mode);
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    }
    if(mode == '\0' || mode == _PGROUPK_XTAL_){                                   //DX20171205 - Added pgroupk_xtal
      char tmp_mode = _PGROUPK_XTAL_;
      sscontent_json << "\"pgroupk_xtal\":" << SymmetryToJson(a.pgroupk_xtal,tmp_mode); //DX20171205 - Added pgroupk_xtal
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");            //DX20171205 - Added pgroupk_xtal
    }                                                                                   //DX20171205 - Added pgroupk_xtal
    //DX20200129 - Patterson - BEGIN
    if(mode == '\0' || mode == _PGROUPK_PATTERSON_){
      char tmp_mode = _PGROUPK_PATTERSON_;
      sscontent_json << "\"pgroupk_Patterson\":" << SymmetryToJson(a.pgroupk_Patterson,tmp_mode);
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    }
    //DX20200129 - Patterson - END
    if(mode == '\0' || mode == _SGROUP_){
      char tmp_mode = _SGROUP_;
      sscontent_json << "\"sgroup\":" << SymmetryToJson(a.sgroup,tmp_mode);
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    }
    if(mode == '\0' || mode == _AGROUP_){
      char tmp_mode = _AGROUP_;
      sscontent_json << "\"agroup\":" << AgroupSymmetryToJson(a.agroup,tmp_mode);
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    }
    if(mode == '\0' || mode == _IATOMS_){
      sscontent_json << "\"iatoms\":" << EquivalentAtomsToJson(a.iatoms);
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    }
    oss << "{" << aurostd::joinWDelimiter(vcontent_json,",")  << "}" << endl;
    return TRUE;
  }
  return FALSE;
}
//DX20170803 - Print symmetry to screen - END
//ME20210402 - Special web output for symmetry
bool KBIN_SymmetryToScreenWeb(xstructure& a, ostream& oss, char mode) {
  stringstream sscontent_txt, sscontent_json;
  vector<string> vcontent_txt;
  switch(mode) {
    case _PGROUP_:
      sscontent_txt << "AFLOW point group, operations are as a=U*b (cols vectors), (Uc/Uf for cartesian/fractional) " << endl;
      sscontent_txt << a.pgroup.size() << "    point group operations " << endl;
      sscontent_txt << SEPARATION_LINE_DASH << endl;
      for(uint k=0;k<a.pgroup.size();k++) {
        sscontent_txt << " Operation number = " << k+1 << " / " << a.pgroup.size() << endl;
        sscontent_txt << a.pgroup[k];
        sscontent_txt << SEPARATION_LINE_DASH << endl;
      }
      sscontent_json << SymmetryToJson(a.pgroup, mode);
      break;
    case _PGROUPK_:
      sscontent_txt << "AFLOW klattice point group, operations are as a=U*b (cols vectors), (Uc/Uf for cartesian/fractional) " << endl;
      sscontent_txt << a.pgroupk.size() << "    point group operations " << endl;
      sscontent_txt << SEPARATION_LINE_DASH << endl;
      for(uint k=0;k<a.pgroupk.size();k++) {
        sscontent_txt << " Operation number = " << k+1 << " / " << a.pgroupk.size() << endl;
        sscontent_txt << a.pgroupk[k];
        sscontent_txt << SEPARATION_LINE_DASH << endl;
      }
      sscontent_json << SymmetryToJson(a.pgroupk, mode);
      break;
    case _FGROUP_:
      sscontent_txt << "AFLOW factor group, operations are as a=U*b+tau (cols vectors), (Uc/Uf and ctau/ftau for cartesian/fractional)" << endl;
      sscontent_txt << a.fgroup.size() << "    factor group operations " << endl;
      sscontent_txt << SEPARATION_LINE_DASH << endl;
      for(uint k=0;k<a.fgroup.size();k++) {
        sscontent_txt << " Operation number = " << k+1 << " / " << a.fgroup.size() << endl;
        sscontent_txt << a.fgroup[k];
        sscontent_txt << SEPARATION_LINE_DASH << endl;
      }
      sscontent_json << SymmetryToJson(a.fgroup, mode);
      break;
    case _PGROUP_XTAL_:
      sscontent_txt << "AFLOW crystal point group, operations are as a=U*b (cols vectors), (Uc/Uf for cartesian/fractional) " << endl;
      sscontent_txt << a.pgroup_xtal.size() << "    point group operations " << endl;
      sscontent_txt << SEPARATION_LINE_DASH << endl;
      for(uint k=0;k<a.pgroup_xtal.size();k++) {
        sscontent_txt << " Operation number = " << k+1 << " / " << a.pgroup_xtal.size() << endl;
        sscontent_txt << a.pgroup_xtal[k];
        sscontent_txt << SEPARATION_LINE_DASH << endl;
      }
      sscontent_json << SymmetryToJson(a.pgroup_xtal, mode);
      break;
    case _PGROUPK_XTAL_:
      sscontent_txt << "AFLOW dual of crystal point group, operations are as a=U*b (cols vectors), (Uc/Uf for cartesian/fractional) " << endl;
      sscontent_txt << a.pgroupk_xtal.size() << "    dual of crystal point group operations " << endl;
      sscontent_txt << SEPARATION_LINE_DASH << endl;
      for(uint k=0;k<a.pgroupk_xtal.size();k++) {
        sscontent_txt << " Operation number = " << k+1 << " / " << a.pgroupk_xtal.size() << endl;
        sscontent_txt << a.pgroupk_xtal[k];
        sscontent_txt << SEPARATION_LINE_DASH << endl;
      }
      sscontent_json << SymmetryToJson(a.pgroupk_xtal, mode);
      break;
    case _PGROUPK_PATTERSON_:
      sscontent_txt << "AFLOW Patterson point group, operations are as a=U*b (cols vectors), (Uc/Uf for cartesian/fractional) " << endl;
      sscontent_txt << a.pgroupk_Patterson.size() << "    point group operations " << endl;
      sscontent_txt << SEPARATION_LINE_DASH << endl;
      for(uint k=0;k<a.pgroupk_Patterson.size();k++) {
        sscontent_txt << " Operation number = " << k+1 << " / " << a.pgroupk_Patterson.size() << endl;
        sscontent_txt << a.pgroupk_Patterson[k];
        sscontent_txt << SEPARATION_LINE_DASH << endl;
      }
      sscontent_json << SymmetryToJson(a.pgroupk_Patterson, mode);
      break;
    case _SGROUP_:
      sscontent_txt << "AFLOW space group, operations are as a=U*b+tau+trasl (cols vectors), (Uc/Uf,ctau/ftau,ctrasl/ftrasl for cartesian/fractional)" << endl;
      sscontent_txt << a.sgroup.size() << "    space group operations " << endl;
      sscontent_txt << a.sgroup_radius << "    radius of space group " << endl;
      sscontent_txt << a.sgroup_radius_dims << "    dimension of radius of space group " << endl;
      for(uint k=0;k<a.sgroup.size();k++) {
        sscontent_txt << " Operation number = " << k+1 << " / " << a.sgroup.size() << endl;
        sscontent_txt << a.sgroup[k];
        sscontent_txt << SEPARATION_LINE_DASH << endl;
      }
      sscontent_json << SymmetryToJson(a.sgroup, mode);
      break;
    case _AGROUP_:
      sscontent_txt << "AFLOW site point group, operations (centered on the site) are as a=U*b (cols vectors), (Uc/Uf for cartesian/fractional) " << endl;
      for(uint iat=0;iat<a.atoms.size();iat++) {
        sscontent_txt << " Site="  << iat << endl;
        sscontent_txt << a.agroup.at(iat).size() << "   site point group operations " << endl;
        sscontent_txt << SEPARATION_LINE_DASH << endl;
        for(uint k=0;k<a.agroup[iat].size();k++) {
          sscontent_txt << " Operation number = " << k+1 << " / " << a.agroup[iat].size() << endl;
          sscontent_txt << a.agroup[iat][k];
          sscontent_txt << SEPARATION_LINE_DASH << endl;
        }
      }
      sscontent_json << AgroupSymmetryToJson(a.agroup, mode);
      break;
    default:
      return false;
  }
  aurostd::stream2vectorstring(sscontent_txt, vcontent_txt);
  oss << "{\"txt\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vcontent_txt, "\"", "\""), ",") << "],"
    << "\"json\":" << sscontent_json.str() << "}" << std::endl;
  return true;
}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
// ------------------------------------------------------------ STEP SYMMETRY
// Function StepSymmetryPerform
//
// This function makes all the symmetry step starting from a structure and
// a bunch of flags
//DX+CO START
bool KBIN_StepSymmetryPerform(xstructure& a,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,const bool& osswrite,ostream& oss) {
  return KBIN_StepSymmetryPerform_20161205(a,AflowIn,FileMESSAGE,aflags,kflags,osswrite,oss);
}
bool KBIN_StepSymmetryPerform_20161205(xstructure& a,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,const bool& osswrite,ostream& oss) {
  stringstream message; //DX20210703
  ostringstream aus;
  bool Krun=TRUE;
  if(aurostd::substring2bool(AflowIn,"[AFLOW_SYMMETRY]CALC",TRUE)) {
    kflags.KBIN_SYMMETRY_CALCULATION=TRUE;
    kflags.KBIN_SYMMETRY_SGROUP_WRITE=aurostd::substring2bool(AflowIn,"[AFLOW_SYMMETRY]SGROUP_WRITE",TRUE);
    kflags.KBIN_SYMMETRY_SGROUP_RADIUS=aurostd::substring2utype<double>(AflowIn,"[AFLOW_SYMMETRY]SGROUP_RADIUS=",TRUE);
  }
  if(aurostd::substring2bool(AflowIn,"[VASP_SYMMETRY]CALC",TRUE)) {
    kflags.KBIN_SYMMETRY_CALCULATION=TRUE;
    kflags.KBIN_SYMMETRY_SGROUP_WRITE=aurostd::substring2bool(AflowIn,"[VASP_SYMMETRY]SGROUP_WRITE",TRUE);
    kflags.KBIN_SYMMETRY_SGROUP_RADIUS=aurostd::substring2utype<double>(AflowIn,"[VASP_SYMMETRY]SGROUP_RADIUS=",TRUE);
  }
  if(kflags.KBIN_SYMMETRY_CALCULATION==TRUE) {
    kflags.KBIN_SYMMETRY_PGROUP_WRITE=TRUE;
    kflags.KBIN_SYMMETRY_PGROUPK_WRITE=FALSE; //default
    kflags.KBIN_SYMMETRY_FGROUP_WRITE=TRUE;
    kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE=TRUE;
    kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE=TRUE; //DX20171205 - Added pgroupk_xtal
    kflags.KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE=TRUE; //DX20200129
    kflags.KBIN_SYMMETRY_IATOMS_WRITE=TRUE;
    kflags.KBIN_SYMMETRY_AGROUP_WRITE=TRUE;
    // calculate the symmetry operations
    message << "Calculating the full set of symmetry operations."; //DX20210703
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_); //DX20210703
    Krun=(Krun && pflow::PerformFullSymmetry(a,kflags.KBIN_SYMMETRY_EPS,kflags.KBIN_SYMMETRY_NO_SCAN,true,FileMESSAGE,aflags,kflags,osswrite,oss));
    message << "Finished calculating the full set of symmetry operations."; //DX20210703
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_); //DX20210703
  }
  if(a.bravais_lattice_type.empty()){
    // calculate the lattice type/variation
    message << "Calculating the lattice information (type, variation, etc.). This may take some time, please be patient."; //DX20210703
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_); //DX20210703
    a.GetRealLatticeType(); //CO+DX20210616 - needed for AEL/AGL
    message << "Finished calculating the lattice information."; //DX20210703
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_); //DX20210703
  }
  return Krun;
}
//DX+CO END

#ifndef COMPILE_SLIM
bool KBIN_StepSymmetryPerform_20160101(xstructure& a,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,const bool& osswrite,ostream& oss) {
  ostringstream aus;
  bool Krun=TRUE;
  if(aurostd::substring2bool(AflowIn,"[AFLOW_SYMMETRY]CALC",TRUE)) {
    kflags.KBIN_SYMMETRY_CALCULATION=TRUE;
    kflags.KBIN_SYMMETRY_SGROUP_WRITE=aurostd::substring2bool(AflowIn,"[AFLOW_SYMMETRY]SGROUP_WRITE",TRUE);
    kflags.KBIN_SYMMETRY_SGROUP_RADIUS=aurostd::substring2utype<double>(AflowIn,"[AFLOW_SYMMETRY]SGROUP_RADIUS=",TRUE);
  }
  if(aurostd::substring2bool(AflowIn,"[VASP_SYMMETRY]CALC",TRUE)) {
    kflags.KBIN_SYMMETRY_CALCULATION=TRUE;
    kflags.KBIN_SYMMETRY_SGROUP_WRITE=aurostd::substring2bool(AflowIn,"[VASP_SYMMETRY]SGROUP_WRITE",TRUE);
    kflags.KBIN_SYMMETRY_SGROUP_RADIUS=aurostd::substring2utype<double>(AflowIn,"[VASP_SYMMETRY]SGROUP_RADIUS=",TRUE);
  }
  if(kflags.KBIN_SYMMETRY_CALCULATION==TRUE) {
    // point group
    kflags.KBIN_SYMMETRY_PGROUP_WRITE=kflags.KBIN_SYMMETRY_CALCULATION;
    Krun=(Krun && SYM::CalculatePointGroup(FileMESSAGE,a,aflags,kflags.KBIN_SYMMETRY_PGROUP_WRITE,osswrite,oss));
    // factor group
    kflags.KBIN_SYMMETRY_FGROUP_WRITE=kflags.KBIN_SYMMETRY_CALCULATION;
    Krun=(Krun && SYM::CalculateFactorGroup(FileMESSAGE,a,aflags,kflags.KBIN_SYMMETRY_FGROUP_WRITE,osswrite,oss));
    // point group
    kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE=kflags.KBIN_SYMMETRY_CALCULATION;
    Krun=(Krun && SYM::CalculatePointGroupCrystal(FileMESSAGE,a,aflags,kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE,osswrite,oss));
    // space group
    if(kflags.KBIN_SYMMETRY_SGROUP_RADIUS>0.0) {
      if(!aflags.QUIET) aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "POSCAR SGROUP: found RADIUS="<<kflags.KBIN_SYMMETRY_SGROUP_RADIUS<<" " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    } else {
      kflags.KBIN_SYMMETRY_SGROUP_RADIUS=KBIN_SYMMETRY_SGROUP_RADIUS_DEFAULT;
      if(!aflags.QUIET) aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "POSCAR SGROUP: Default RADIUS="<<kflags.KBIN_SYMMETRY_SGROUP_RADIUS<<" " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    }
    a.sgroup_radius=kflags.KBIN_SYMMETRY_SGROUP_RADIUS;
    Krun=(Krun && SYM::CalculateSpaceGroup(FileMESSAGE,a,aflags,kflags.KBIN_SYMMETRY_SGROUP_WRITE,osswrite,oss));
    // inequivalente atoms
    kflags.KBIN_SYMMETRY_IATOMS_WRITE=kflags.KBIN_SYMMETRY_CALCULATION;
    Krun=(Krun && SYM::CalculateInequivalentAtoms(FileMESSAGE,a,aflags,kflags.KBIN_SYMMETRY_IATOMS_WRITE,osswrite,oss));
    // site point group
    kflags.KBIN_SYMMETRY_AGROUP_WRITE=kflags.KBIN_SYMMETRY_CALCULATION;
    Krun=(Krun && SYM::CalculateSitePointGroup(FileMESSAGE,a,aflags,kflags.KBIN_SYMMETRY_AGROUP_WRITE,osswrite,oss));
  }

  return Krun;
}
#endif

// ******************************************************************************
// SYM::writePythonScript() //DX20201228
// ******************************************************************************
namespace SYM {
  void writePythonScript(ostream& oss){

    // Writes AFLOW-SYM Python script in a subdirectory


    string directory = aurostd::getPWD();
    string sym_python_subdir = "AFLOW_SYM_PYTHON";
    string python_dir = directory + "/" + sym_python_subdir;

    aurostd::DirectoryMake(python_dir);

    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Writing out python script to: "+python_dir, oss, _LOGGER_NOTICE_);

    aurostd::EmbData::save_to_file("aflow_sym_python.py", python_dir+"/"+"aflow_sym_python.py");
  }
}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
