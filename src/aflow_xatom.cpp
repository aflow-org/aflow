// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - 2007-2019
#ifndef _AFLOW_XATOM_CPP
#define _AFLOW_XATOM_CPP
#include "aflow.h"
#include "aflow_pflow.h"
#include "aflow_symmetry_spacegroup.h" //DX20180723
#include "AUROSTD/aurostd_xscalar.h"
#include "aflow_compare_structure.h" //CO20180409
#include "aflow_chull.h" //HE20210408
#include "aflow_symbolic.h"  //ME20220124

#define _calculate_symmetry_default_sgroup_radius_   2.0
#define PLATON_MIN_VOLUME_PER_ATOM   6.0   // for symmetry calculation
#define _PRIM_MULTITHREAD_MIN_ATOMS_THRESHOLD_ 50

#define PLATON_TOLERANCE_ANGLE 1.0
#define PLATON_TOLERANCE_D1 0.25
#define PLATON_TOLERANCE_D2 0.25
#define PLATON_TOLERANCE_D3 0.25

#define _EPS_ 0.02

#define DEBUG_LATTICE_DIMENSIONS false //DX20210401

//#define _pocc_precision_ 5  //CO20170630 - we should really be fixing default vs. fixed
// tested right
// ./aflow --zval=/common/LIB3/LIB/AgCdZn/TFCC013.ABC
// aflow.3.1.184 --zval=/common/LIB3/LIB/AgCdZn/TFCC013.ABC
// ./aflow --zval_cell=/common/LIB3/LIB/AgCdZn/TFCC013.ABC
// aflow.3.1.184 --zval_cell=/common/LIB3/LIB/AgCdZn/TFCC013.ABC
// ./aflow --zval_atom=/common/LIB3/LIB/AgCdZn/TFCC013.ABC
// aflow.3.1.184 --zval_atom=/common/LIB3/LIB/AgCdZn/TFCC013.ABC
// ./aflow --pomass=/common/LIB3/LIB/AgCdZn/TFCC013.ABC
// aflow.3.1.184 --pomass=/common/LIB3/LIB/AgCdZn/TFCC013.ABC
// ./aflow --pomass_cell=/common/LIB3/LIB/AgCdZn/TFCC013.ABC
// aflow.3.1.184 --pomass_cell=/common/LIB3/LIB/AgCdZn/TFCC013.ABC
// ./aflow --pomass_atom=/common/LIB3/LIB/AgCdZn/TFCC013.ABC
// aflow.3.1.184 --pomass_atom=/common/LIB3/LIB/AgCdZn/TFCC013.ABC

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// _XATOM_PROTOTYPES_

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// _ATOM
// look into aflow.h for the definitions

// constructors
_atom::_atom() {free();} 

// destructor
_atom::~_atom() {free();}


void _atom::free() { // PRIVATE //HE20220826 changed constructor to use free() 
  fpos.clear();
  cpos.clear();
  corigin.clear();
  coord.clear();
  fpos_equation.clear(); //DX20180607 - symbolic math for atom positions
  cpos_equation.clear(); //DX20180607 - symbolic math for atom positions
  spin=0.0;
  spin_is_given=FALSE; //DX20170921 - magnetic sym
  noncoll_spin.clear(); //DX20171205 - magnetic sym (non-collinear)
  noncoll_spin_is_given=FALSE; //DX20171205 - magnetic sym (non-collinear)
  mass=0.0;
  type=0;
  name="";
  name_is_given=FALSE;
  cleanname="";
  info=0;  // (RHT)
  atomic_number=0;
  //[CO20200130 - number->basis]number=0;
  sd="";
  ijk.clear();
  isincell=FALSE;
  basis=0;
  reference=0.0;
  ireference=0;
  equivalent=-1;
  is_inequivalent=TRUE;
  num_equivalents=0;
  index_iatoms=0;
  order_parameter_value=0;
  order_parameter_atom=FALSE;
  partial_occupation_value=1.0;
  partial_occupation_flag=FALSE;
  shell=0;
  force.clear(); //CO20211107
  verbose=FALSE;
  print_RHT=false;  //CO20190405 //true; //CHANGE THIS BACK TO FALSE WHEN DONE DEBUGGING  (RHT)
  print_cartesian=FALSE;
}

void _atom::copy(const _atom& b) { // copy PRIVATE
  fpos=b.fpos;
  cpos=b.cpos;
  corigin=b.corigin;
  coord=b.coord;
  fpos_equation=b.fpos_equation; //DX20180607 - symbolic math for atom positions
  cpos_equation=b.cpos_equation; //DX20180607 - symbolic math for atom positions
  spin=b.spin;
  spin_is_given=b.spin_is_given; //DX20170921 - magnetic sym
  noncoll_spin=b.noncoll_spin; //DX20171205 - magnetic sym (non-collinear)
  noncoll_spin_is_given=b.noncoll_spin_is_given; //DX20171205 - magnetic sym (non-collinear)
  mass=b.mass;
  type=b.type;
  name=b.name;
  name_is_given=b.name_is_given;
  cleanname=b.cleanname;
  info=b.info;  // (RHT)
  atomic_number=b.atomic_number;
  //[CO20200130 - number->basis]number=b.number;
  sd=b.sd;
  ijk=b.ijk;
  isincell=b.isincell;
  basis=b.basis;
  reference=b.reference;
  ireference=b.ireference;
  equivalent=b.equivalent;
  is_inequivalent=b.is_inequivalent;
  num_equivalents=b.num_equivalents;
  index_iatoms=b.index_iatoms;
  order_parameter_value=b.order_parameter_value;
  order_parameter_atom=b.order_parameter_atom;
  partial_occupation_value=b.partial_occupation_value;
  partial_occupation_flag=b.partial_occupation_flag;
  shell=b.shell;
  force=b.force;  //CO20211107
  verbose=b.verbose;
  print_RHT=b.print_RHT;  // (RHT)
  print_cartesian=b.print_cartesian;
}

const _atom& _atom::operator=(const _atom& b) {  // operator= PUBLIC
  if(this!=&b) {copy(b);}
  return *this;
}

_atom::_atom(const _atom& b) { // copy PUBLIC
  //  free(); *this=b;
  copy(b);
}

void _atom::clear(){free();}

ostream& operator<<(ostream& oss,const _atom& atom) {
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(10);
  if(atom.print_RHT==true) {
    //oss << "ATOM COUT-RHT" << endl;
    oss << atom.coord << " " << atom.name;
  }
  else {
    if(atom.verbose==TRUE) {
      oss << " " << endl;
      oss << "ATOM COUT" << endl;
      oss << "type=" << atom.type << endl;
      oss << "spin=" << atom.spin << endl;
      oss << "spin_is_given=" << atom.spin_is_given << endl; //DX20170921 - magnetic sym
      oss << "noncoll_spin=" << atom.noncoll_spin << endl;                   //DX20171205 - magnetic sym (non-collinear)
      oss << "noncoll_spin_is_given=" << atom.noncoll_spin_is_given << endl; //DX20171205 - magnetic sym (non-collinear)
      oss << "mass=" << atom.mass << endl;
      oss << "name=" << atom.name << endl;
      oss << "info=" << atom.info << endl;
      oss << "cleanname=" << atom.cleanname << endl;
      oss << "atomic_number=" << atom.atomic_number << endl;
      oss << "basis=" << atom.basis << endl;  //[CO20200130 - number->basis]oss << "number=" << atom.number << endl;
      oss << "name_is_given=" << atom.name_is_given << endl;
      oss << "sd=" <<  atom.sd << endl;
      oss << "print_cartesian" << atom.print_cartesian << endl;
      oss << "print_RHT" << atom.print_RHT << endl;
      oss << "corigin" << atom.corigin(1) << " " << atom.corigin(2) << " " << atom.corigin(3) << endl;
      oss << "coord" << atom.coord(1) << " " << atom.coord(2) << " " << atom.coord(3) << endl;

      oss << "fpos_equation" << aurostd::joinWDelimiter(atom.fpos_equation," ") << endl; //DX20180607 - symbolic math for atom positions //DX20191218 - join with delimiter in case empty
      oss << "cpos_equation" << aurostd::joinWDelimiter(atom.cpos_equation," ") << endl; //DX20180607 - symbolic math for atom positions //DX20191218 - join with delimiter in case empty
      oss << "isincell=" << atom.isincell << endl;
      oss << "reference=" << atom.reference << endl;
      oss << "ireference=" << atom.ireference << endl;
      oss << "equivalent=" << atom.equivalent << endl;
      oss << "is_inequivalent=" << atom.is_inequivalent << endl;
      oss << "num_equivalents=" << atom.num_equivalents << endl;
      oss << "index_iatoms=" << atom.index_iatoms << endl;
      oss << "verbose=" << atom.verbose << endl;
      oss << "order_parameter_value=" << atom.order_parameter_value << endl;
      oss << "order_parameter_atom=" << atom.order_parameter_atom << endl;
      oss << "partial_occupation_value=" << atom.partial_occupation_value << endl;
      oss << "partial_occupation_flag=" << atom.partial_occupation_flag << endl;
      oss << "nearest_neighbor_shell_num= " << atom.shell << endl;
      oss << "force= " << atom.force << endl;
    }
    if(atom.print_cartesian==TRUE) {
      if(atom.verbose) oss << "cartesian" << endl;
      oss << "C " << atom.cpos(1) << " " << atom.cpos(2) << " " << atom.cpos(3);// << endl;
    } else {
      if(atom.verbose) oss << "fractional" << endl;
      oss << "F " << atom.fpos(1) << " " << atom.fpos(2) << " " << atom.fpos(3);//  << endl;
    }
    oss << " T=" << atom.type;
    oss << " B=" << atom.basis;
    //[CO20200130 - number->basis]oss << " N=" << atom.number;
    //  oss << setw(1);
    oss << " ijk=[" << atom.ijk(1) << "," << atom.ijk(2) << "," << atom.ijk(3) << "]";
    if(atom.verbose) oss << " " << endl;
  }
  return oss;
}

//DX20190214 [OBSOLETE]bool isequalRHT(const _atom& a, const _atom& b,double tol) {
//DX20190214 [OBSOLETE]  bool out = false;
//DX20190214 [OBSOLETE]  //DX+CO START
//DX20190214 [OBSOLETE]  if(abs(a.fpos(1)-b.fpos(1)) < tol && 
//DX20190214 [OBSOLETE]     abs(a.fpos(2)-b.fpos(2)) < tol && 
//DX20190214 [OBSOLETE]     abs(a.fpos(3)-b.fpos(3)) < tol &&
//DX20190214 [OBSOLETE]     //DX+CO END
//DX20190214 [OBSOLETE]     a.name == b.name) {
//DX20190214 [OBSOLETE]    out = true;
//DX20190214 [OBSOLETE]  }
//DX20190214 [OBSOLETE]  return out;
//DX20190214 [OBSOLETE]}


void _atom::CleanName(void) {
  //CO20200624 - cleaning up function
  //the old function had a sign problem and was inefficient
  cleanname=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(name);
  aurostd::VASP_PseudoPotential_CleanName_InPlace(cleanname); //DX20200907 - changed from KBIN to aurostd (since this is now in xparser)
  if(cleanname.size()>3){cleanname=cleanname.substr(0,3);}  //cannot be longer than 3 characters
  for(uint i=cleanname.size()-1;i<cleanname.size();i--){  //go backwards and clean anything that isn't between A-Z and a-z 
    //65-90 is A-Z
    //97-122 is a-z
    if(!((cleanname[i]>=65 && cleanname[i]<=90)||(cleanname[i]>=97 && cleanname[i]<=122))){
      cleanname.erase(cleanname.begin()+i);
    }
  }
  //the old function insists on correcting if we have aL instead of Al...
  if(cleanname.size()>0 && cleanname[0]>=97 && cleanname[0]<=122){cleanname[0]-=('a'-'A');}
  for(uint i=1;i<cleanname.size();i++){
    if(cleanname[i]>=65 && cleanname[i]<=90){cleanname[i]+=('a'-'A');}
  }

  //[CO20200624 - OBSOLETE]// This function cleanup the name from VASP stuff
  //[CO20200624 - OBSOLETE]string name1,name2;
  //[CO20200624 - OBSOLETE]name1=name+" ";
  //[CO20200624 - OBSOLETE]name2=name+" ";
  //[CO20200624 - OBSOLETE]name1=name1.substr(0,1);  //grab first 2 characters
  //[CO20200624 - OBSOLETE]//65-90 is A-Z
  //[CO20200624 - OBSOLETE]//97-122 is a-z
  //[CO20200624 - OBSOLETE]if(!((name1[0]>=65 && name1[0]<=90)||(name1[0]>=97 && name1[0]<=122))) name1="";
  //[CO20200624 - OBSOLETE]if((name1[0]>=97 && name1[0]<=122)) name1[0]-=-97+65; //THE SIGN LOOKS WRONG!!! - if we have 97, then it would subtract -32: giving us 129
  //[CO20200624 - OBSOLETE]name2=name2.substr(1,1);  //grab second and third character
  //[CO20200624 - OBSOLETE]if(!((name2[0]>=65 && name2[0]<=90)||(name2[0]>=97 && name2[0]<=122))) name2="";
  //[CO20200624 - OBSOLETE]if((name2[0]>=65 && name2[0]<=90)) name2[0]+=-97+65;  //THE SIGN LOOKS WRONG!!! - if we have 65, then it would add -32: giving us 33 
  //[CO20200624 - OBSOLETE]cleanname=name1+name2;
  //[CO20200624 - OBSOLETE]atomic_number=0;

  for(uint j=0;j< vatom_symbol.size();j++) if(cleanname==vatom_symbol.at(j)) atomic_number=j;
}

void _atom::CleanSpin(void) {
  spin=0.0;
  spin_is_given=FALSE; //DX20170921 - magnetic sym
  noncoll_spin.clear();            //DX20171205 - magnetic sym (non-collinear)
  noncoll_spin_is_given=FALSE; //DX20171205 - magnetic sym (non-collinear)
  if(!XHOST.READ_SPIN_FROM_ATOMLABEL) {return;} //SD20220316
  if(name.find("+")!=string::npos) {spin=atof(name.substr(name.find("+")).c_str()); spin_is_given=TRUE;} //DX20170921 - magnetic sym
  if(name.find("-")!=string::npos) {spin=atof(name.substr(name.find("-")).c_str()); spin_is_given=TRUE;} //DX20170921 - magnetic sym
}

void _atom::ClearSymmetry(void) { //CO20190219
  (*this).equivalent=-1;
  (*this).is_inequivalent=TRUE;
  (*this).num_equivalents=0;
}

std::vector<string> vatom_symbol(NUM_ELEMENTS);        // store starting from ONE
std::vector<string> vatom_name(NUM_ELEMENTS);          // store starting from ONE
std::vector<double> vatom_mass(NUM_ELEMENTS);          // store starting from ONE
std::vector<double> vatom_volume(NUM_ELEMENTS);        // store starting from ONE
std::vector<int> vatom_valence_iupac(NUM_ELEMENTS);    // store starting from ONE http://en.wikipedia.org/wiki/Valence_(chemistry)
std::vector<int> vatom_valence_std(NUM_ELEMENTS);      // store starting from ONE http://en.wikipedia.org/wiki/Valence_(chemistry)
std::vector<double> vatom_miedema_phi_star(NUM_ELEMENTS);  // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28  
std::vector<double> vatom_miedema_nws(NUM_ELEMENTS);       // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
std::vector<double> vatom_miedema_Vm(NUM_ELEMENTS);        // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
std::vector<double> vatom_miedema_gamma_s(NUM_ELEMENTS);   // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
std::vector<double> vatom_miedema_BVm(NUM_ELEMENTS);       // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
// for lanthines from J.A. Alonso and N.H. March. Electrons in Metals and Alloys, Academic Press, London (1989) (except La)
std::vector<double> vatom_radius(NUM_ELEMENTS);        // store starting from ONE
std::vector<double> vatom_radius_covalent(NUM_ELEMENTS);// store starting from ONE //DX+CO20170904 
std::vector<double> vatom_electronegativity(NUM_ELEMENTS);        // store starting from ONE
std::vector<string> vatom_crystal(NUM_ELEMENTS);        // store starting from ONE
std::vector<double> vatom_xray_scatt(NUM_ELEMENTS);         // store starting from ONE
std::vector<double> vatom_pettifor_scale(NUM_ELEMENTS);         // store starting from ONE Chemical Scale Pettifor Solid State Communications 51 31-34 1984
std::vector<double> vatom_pearson_coefficient(NUM_ELEMENTS);    //ME20181020 Pearson mass deviation coefficient

void atoms_initialize(void) {
  for(int i=0;i<NUM_ELEMENTS;i++) {       // clear
    vatom_xray_scatt.at(i)=(double)i-1;      // shift+1
    vatom_mass.at(i)=AMU2KILOGRAM;      // masses in kilos
    vatom_volume.at(i)=NNN;             // atomic volume in A^3 from the FCC vasp table and/or successive calculations
    vatom_valence_iupac.at(i)=NNN;            // IUPAC Maximum number of univalent atoms that may combine with an atom of the element under consideration, or with a fragment, or for which an atom of this element can be substituted.
    vatom_valence_std.at(i)=NNN;            // stanmdard: number electrons minus closed shell at leff (noble gas)
    vatom_miedema_phi_star.at(i)=NNN;                 // Miedema Rule Table 1a Physica 100B (1980) 1-28   (phi^\star in (V))
    vatom_miedema_nws.at(i)=NNN;                      // Miedema Rule Table 1a Physica 100B (1980) 1-28   n_{ws}^{1/3} in (d.u.)^1/3
    vatom_miedema_Vm.at(i)=NNN;                       // Miedema Rule Table 1a Physica 100B (1980) 1-28   V_m^{2/3} in (cm^2)
    vatom_miedema_gamma_s.at(i)=NNN;                  // Miedema Rule Table 1a Physica 100B (1980) 1-28   \gamma_s^0 in (mJ/m^2)
    vatom_miedema_BVm.at(i)=NNN;                      // Miedema Rule Table 1a Physica 100B (1980) 1-28   BV_m (kJ/mole)
    vatom_radius.at(i)=NNN;             // Saxena (nm)
    vatom_radius_covalent.at(i)=NNN;    // Codero (Angstroms) //DX+CO20170904
    vatom_electronegativity.at(i)=NNN;  // Saxena
    vatom_crystal.at(i)="nnn";          // Ashcroft-Mermin
  }

  // Xray_scatt_vector All data collected from the NIST online tables
  // http://physics.nist.gov/PhysRefData/FFast/html/form.html
  // All data are ideally for f1 values for Cu-alpha (wavelength=1.5418A, E=8.0416keV).
  // These are for E=7.9026keV (Cu-alpha is wavelength=1.5418A, E=8.0416keV).

  // All data collected from the online tables:
  // http://www-cxro.lbl.gov/optical_constants/pert_form.html
  // All data are f1 values for Cu-alpha (wavelength=1.5418A, E=8.0416keV].

  int i;
  // ROW 0
  i=0; vatom_symbol[i]="XX"; vatom_name[i]="UNDEFINED";    vatom_mass[i]*=0.00;      vatom_volume[i]=NNN;      vatom_valence_std[i]=0; vatom_valence_iupac[i]=0;   vatom_miedema_phi_star[i]=NNN;   vatom_miedema_nws[i]=NNN;   vatom_miedema_Vm[i]=NNN;  vatom_miedema_gamma_s[i]= NNN;   vatom_miedema_BVm[i]=NNN;   vatom_radius[i]=NNN;    vatom_radius_covalent[i]=NNN;   vatom_electronegativity[i]=NNN;  vatom_crystal[i]="nnn"; vatom_pettifor_scale[i]=0; vatom_xray_scatt[i]=0; vatom_pearson_coefficient[i]=0.0;

  // ROW 1
  // s-electron systems
  i++; vatom_symbol[i]="H";  vatom_name[i]="Hydrogen";     vatom_mass[i]*=1.0079;    vatom_volume[i]= 0.75110; vatom_valence_std[i]=1; vatom_valence_iupac[i]=1;  vatom_miedema_phi_star[i]=5.2;   vatom_miedema_nws[i]=1.5;   vatom_miedema_Vm[i]=NNN;  vatom_miedema_gamma_s[i]= NNN;   vatom_miedema_BVm[i]=NNN;   vatom_radius[i]=0.046;  vatom_radius_covalent[i]=0.31;   vatom_electronegativity[i]=2.10;  vatom_crystal[i]="hex";  vatom_pettifor_scale[i]=0; vatom_xray_scatt[i]=1.000;vatom_pearson_coefficient[i]=0.00011460743;  // H volume wrong *dimer*   MIEDEMA = PAUL VAN DER PUT book
  i++; vatom_symbol[i]="He"; vatom_name[i]="Helium";       vatom_mass[i]*=4.0026;    vatom_volume[i]= -1.000;  vatom_valence_std[i]=0; vatom_valence_iupac[i]=0;  vatom_miedema_phi_star[i]=NNN;   vatom_miedema_nws[i]=NNN;   vatom_miedema_Vm[i]=NNN;  vatom_miedema_gamma_s[i]= NNN;   vatom_miedema_BVm[i]=NNN;   vatom_radius[i]=NNN;    vatom_radius_covalent[i]=0.28;   vatom_electronegativity[i]=NNN;   vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=0; vatom_xray_scatt[i]=2.000; vatom_pearson_coefficient[i]=8.32328E-8;  // He

  // ROW2
  // s-electron systems
  i++; vatom_symbol[i]="Li"; vatom_name[i]="Lithium";      vatom_mass[i]*=6.941;     vatom_volume[i]=20.24110; vatom_valence_std[i]=1; vatom_valence_iupac[i]=1;  vatom_miedema_phi_star[i]=2.85;  vatom_miedema_nws[i]=0.98;  vatom_miedema_Vm[i]=5.5;  vatom_miedema_gamma_s[i]= 530;   vatom_miedema_BVm[i]=1.5;   vatom_radius[i]=0.152;  vatom_radius_covalent[i]=1.28;   vatom_electronegativity[i]=0.98;  vatom_crystal[i]="bcc";  vatom_pettifor_scale[i]=0.45; vatom_xray_scatt[i]=3.00145; vatom_pearson_coefficient[i]=0.0014588232;  // Li
  i++; vatom_symbol[i]="Be"; vatom_name[i]="Beryllium";    vatom_mass[i]*=9.0122;    vatom_volume[i]= 7.83290; vatom_valence_std[i]=2; vatom_valence_iupac[i]=2;  vatom_miedema_phi_star[i]=4.20;  vatom_miedema_nws[i]=1.60;  vatom_miedema_Vm[i]=2.9;  vatom_miedema_gamma_s[i]=1900;   vatom_miedema_BVm[i]=4.9;   vatom_radius[i]=0.114;  vatom_radius_covalent[i]=0.96;   vatom_electronegativity[i]=1.57;  vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=1.50; vatom_pearson_coefficient[i]=0.0;  // Be
  // p-electron systems
  i++; vatom_symbol[i]="B";  vatom_name[i]="Boron";        vatom_mass[i]*=10.81;     vatom_volume[i]= 5.88420; vatom_valence_std[i]=3; vatom_valence_iupac[i]=3;  vatom_miedema_phi_star[i]=4.75;  vatom_miedema_nws[i]=1.55;  vatom_miedema_Vm[i]=2.8;  vatom_miedema_gamma_s[i]= NNN;   vatom_miedema_BVm[i]=NNN;   vatom_radius[i]=0.097;  vatom_radius_covalent[i]=0.84;   vatom_electronegativity[i]=2.04;  vatom_crystal[i]="tet";  vatom_pettifor_scale[i]=2.00; vatom_pearson_coefficient[i]=0.00135391428;  // B
  i++; vatom_symbol[i]="C";  vatom_name[i]="Carbon";       vatom_mass[i]*=12.011;    vatom_volume[i]= 5.59490; vatom_valence_std[i]=4; vatom_valence_iupac[i]=4;  vatom_miedema_phi_star[i]=6.20;  vatom_miedema_nws[i]=1.90;  vatom_miedema_Vm[i]=1.8;  vatom_miedema_gamma_s[i]= NNN;   vatom_miedema_BVm[i]=NNN;   vatom_radius[i]=0.077;  vatom_radius_covalent[i]=0.76;   vatom_electronegativity[i]=2.55;  vatom_crystal[i]="dia";  vatom_pettifor_scale[i]=2.50; vatom_xray_scatt[i]=6.019; vatom_pearson_coefficient[i]=0.00007387218;   // C  //DX+CO20170904 vatom_radius_covalent uses sp3 hybridization (most common)
  i++; vatom_symbol[i]="N";  vatom_name[i]="Nitrogen";     vatom_mass[i]*=14.0067;   vatom_volume[i]= 7.59940; vatom_valence_std[i]=5; vatom_valence_iupac[i]=5;  vatom_miedema_phi_star[i]=7.00;  vatom_miedema_nws[i]=1.60;  vatom_miedema_Vm[i]=2.2;  vatom_miedema_gamma_s[i]= NNN;   vatom_miedema_BVm[i]=NNN;   vatom_radius[i]=0.071;  vatom_radius_covalent[i]=0.71;   vatom_electronegativity[i]=3.04;  vatom_crystal[i]="hex";  vatom_pettifor_scale[i]=3.00; vatom_pearson_coefficient[i]=0.00001857771; // N JX CHANGED VALENCE
  i++; vatom_symbol[i]="O";  vatom_name[i]="Oxygen";       vatom_mass[i]*=15.9994;   vatom_volume[i]= 7.78230; vatom_valence_std[i]=6; vatom_valence_iupac[i]=2;  vatom_miedema_phi_star[i]=6.97;  vatom_miedema_nws[i]=1.70;  vatom_miedema_Vm[i]=2.656;vatom_miedema_gamma_s[i]= NNN;   vatom_miedema_BVm[i]=NNN;   vatom_radius[i]=0.060;  vatom_radius_covalent[i]=0.66;   vatom_electronegativity[i]=3.44;  vatom_crystal[i]="cub";  vatom_pettifor_scale[i]=3.50; vatom_xray_scatt[i]=8.052; vatom_pearson_coefficient[i]=0.00003358805;  // O Table 27 of JX
  i++; vatom_symbol[i]="F";  vatom_name[i]="Fluorine";     vatom_mass[i]*=18.9984;   vatom_volume[i]= 9.99090; vatom_valence_std[i]=7; vatom_valence_iupac[i]=1;  vatom_miedema_phi_star[i]=NNN;   vatom_miedema_nws[i]=NNN;   vatom_miedema_Vm[i]=NNN;  vatom_miedema_gamma_s[i]= NNN;   vatom_miedema_BVm[i]=NNN;   vatom_radius[i]=NNN;    vatom_radius_covalent[i]=0.57;   vatom_electronegativity[i]=3.98;  vatom_crystal[i]="mcl";  vatom_pettifor_scale[i]=4.00; vatom_pearson_coefficient[i]=0.0;  // F
  i++; vatom_symbol[i]="Ne";vatom_name[i]="Neon";        vatom_mass[i]*=20.179;   vatom_volume[i]=19.9052; vatom_valence_std[i]=0; vatom_valence_iupac[i]=0; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.160; vatom_radius_covalent[i]=0.58;   vatom_electronegativity[i]=NNN; vatom_crystal[i]="fcc";  vatom_pettifor_scale[i]=0; vatom_pearson_coefficient[i]=0.00082783369; // Ne volume calculated with fcc-pawpbe

  // ROW3
  // s-electron systems
  i++; vatom_symbol[i]="Na";vatom_name[i]="Sodium";      vatom_mass[i]*=22.9898;  vatom_volume[i]=36.9135; vatom_valence_std[i]=1; vatom_valence_iupac[i]=1; vatom_miedema_phi_star[i]=2.70; vatom_miedema_nws[i]=0.82; vatom_miedema_Vm[i]=8.3; vatom_miedema_gamma_s[i]= 260;  vatom_miedema_BVm[i]=1.6;  vatom_radius[i]=0.186; vatom_radius_covalent[i]=1.66;   vatom_electronegativity[i]=0.93; vatom_crystal[i]="bcc";  vatom_pettifor_scale[i]=0.40; vatom_pearson_coefficient[i]=0.0;  // Na
  i++; vatom_symbol[i]="Mg";vatom_name[i]="Magnesium";   vatom_mass[i]*=24.305;   vatom_volume[i]=22.8178; vatom_valence_std[i]=2; vatom_valence_iupac[i]=2; vatom_miedema_phi_star[i]=3.45; vatom_miedema_nws[i]=1.17; vatom_miedema_Vm[i]=5.8; vatom_miedema_gamma_s[i]= 790;  vatom_miedema_BVm[i]=5.0;  vatom_radius[i]=0.160; vatom_radius_covalent[i]=1.41;   vatom_electronegativity[i]=1.31; vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=1.28; vatom_pearson_coefficient[i]=0.00073988271;  // Mg
  // p-electron systems
  i++; vatom_symbol[i]="Al";vatom_name[i]="Aluminium";   vatom_mass[i]*=26.9815;  vatom_volume[i]=16.4000; vatom_valence_std[i]=3; vatom_valence_iupac[i]=3; vatom_miedema_phi_star[i]=4.20; vatom_miedema_nws[i]=1.39; vatom_miedema_Vm[i]=4.6; vatom_miedema_gamma_s[i]=1200;  vatom_miedema_BVm[i]=7.2;  vatom_radius[i]=0.143; vatom_radius_covalent[i]=1.21;   vatom_electronegativity[i]=1.61; vatom_crystal[i]="fcc";  vatom_pettifor_scale[i]=1.66; vatom_pearson_coefficient[i]=0.0;  // Al
  i++; vatom_symbol[i]="Si";vatom_name[i]="Silicon";     vatom_mass[i]*=28.0855;  vatom_volume[i]=14.3536; vatom_valence_std[i]=4; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=4.70; vatom_miedema_nws[i]=1.50; vatom_miedema_Vm[i]=4.2; vatom_miedema_gamma_s[i]=1290;  vatom_miedema_BVm[i]=11.9; vatom_radius[i]=0.117; vatom_radius_covalent[i]=1.11;   vatom_electronegativity[i]=1.90; vatom_crystal[i]="dia";  vatom_pettifor_scale[i]=1.92; vatom_xray_scatt[i]=14.43; vatom_pearson_coefficient[i]=0.00020046752; // Si ???
  i++; vatom_symbol[i]="P"; vatom_name[i]="Phosphorus";  vatom_mass[i]*=30.9738;  vatom_volume[i]=14.1995; vatom_valence_std[i]=5; vatom_valence_iupac[i]=5; vatom_miedema_phi_star[i]=5.5;  vatom_miedema_nws[i]=1.65; vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.109; vatom_radius_covalent[i]=1.07;   vatom_electronegativity[i]=2.19; vatom_crystal[i]="cub";  vatom_pettifor_scale[i]=2.18; vatom_xray_scatt[i]=15.3133; vatom_pearson_coefficient[i]=0.0; // P   MIEDEMA = PAUL VAN DER PUT book
  i++; vatom_symbol[i]="S"; vatom_name[i]="Sulphur";     vatom_mass[i]*=32.06;    vatom_volume[i]=15.7301; vatom_valence_std[i]=6; vatom_valence_iupac[i]=6; vatom_miedema_phi_star[i]=5.6;  vatom_miedema_nws[i]=1.46; vatom_miedema_Vm[i]=4.376;vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.106; vatom_radius_covalent[i]=1.05;   vatom_electronegativity[i]=2.58; vatom_crystal[i]="orc";  vatom_pettifor_scale[i]=2.44; vatom_pearson_coefficient[i]=0.00016807795;  // S Table 27 of JX
  i++; vatom_symbol[i]="Cl";vatom_name[i]="Chlorine";    vatom_mass[i]*=35.453;   vatom_volume[i]=21.2947; vatom_valence_std[i]=7; vatom_valence_iupac[i]=7; vatom_miedema_phi_star[i]=5.32;  vatom_miedema_nws[i]=0.34;  vatom_miedema_Vm[i]=6.71; vatom_miedema_gamma_s[i]= 1013;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.107; vatom_radius_covalent[i]=1.02;   vatom_electronegativity[i]=3.16; vatom_crystal[i]="orc";  vatom_pettifor_scale[i]=2.70; vatom_pearson_coefficient[i]=0.00058238731;  // Cl interpolation phi_star, nws, Vm, gamma  JX CHANGED VALENCE    
  i++; vatom_symbol[i]="Ar";vatom_name[i]="Argon";       vatom_mass[i]*=39.948;   vatom_volume[i]=22.000; vatom_valence_std[i]=0; vatom_valence_iupac[i]=2; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.192; vatom_radius_covalent[i]=1.06;   vatom_electronegativity[i]=NNN;  vatom_crystal[i]="fcc";  vatom_pettifor_scale[i]=0; vatom_pearson_coefficient[i]=0.00003509919;  // Ar guessed volume, must double check from results  JX CHANGED VALENCE    

  // ROW4
  // s-electron systems
  i++; vatom_symbol[i]="K"; vatom_name[i]="Potassium";   vatom_mass[i]*=39.0983;  vatom_volume[i]=73.9091; vatom_valence_std[i]=1; vatom_valence_iupac[i]=1; vatom_miedema_phi_star[i]=2.25; vatom_miedema_nws[i]=0.65; vatom_miedema_Vm[i]=12.8;vatom_miedema_gamma_s[i]= 150;  vatom_miedema_BVm[i]=1.5;  vatom_radius[i]=0.231; vatom_radius_covalent[i]=2.03;   vatom_electronegativity[i]=0.82; vatom_crystal[i]="fcc";  vatom_pettifor_scale[i]=0.35; vatom_pearson_coefficient[i]=0.000164;  // K
  i++; vatom_symbol[i]="Ca";vatom_name[i]="Calcium";     vatom_mass[i]*=40.08;    vatom_volume[i]=42.1927; vatom_valence_std[i]=2; vatom_valence_iupac[i]=2; vatom_miedema_phi_star[i]=2.55; vatom_miedema_nws[i]=0.91; vatom_miedema_Vm[i]=8.8; vatom_miedema_gamma_s[i]= 490;  vatom_miedema_BVm[i]=4.0;  vatom_radius[i]=0.197; vatom_radius_covalent[i]=1.76;   vatom_electronegativity[i]=1.00; vatom_crystal[i]="bcc";  vatom_pettifor_scale[i]=0.60; vatom_pearson_coefficient[i]=0.000297564;  // Ca
  // d-electron systems: transition metals
  i++; vatom_symbol[i]="Sc";vatom_name[i]="Scandium";    vatom_mass[i]*=44.9559;  vatom_volume[i]=24.6739; vatom_valence_std[i]=3; vatom_valence_iupac[i]=3; vatom_miedema_phi_star[i]=3.25; vatom_miedema_nws[i]=1.27; vatom_miedema_Vm[i]=6.1; vatom_miedema_gamma_s[i]=1200;  vatom_miedema_BVm[i]=6.6;  vatom_radius[i]=0.160; vatom_radius_covalent[i]=1.70;   vatom_electronegativity[i]=1.36; vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=0.74; vatom_xray_scatt[i]=21.34; vatom_pearson_coefficient[i]=0.0;  // Sc
  i++; vatom_symbol[i]="Ti";vatom_name[i]="Titanium";    vatom_mass[i]*=47.9;     vatom_volume[i]=17.1035; vatom_valence_std[i]=4; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=3.65; vatom_miedema_nws[i]=1.47; vatom_miedema_Vm[i]=4.8; vatom_miedema_gamma_s[i]=2050;  vatom_miedema_BVm[i]=11.0; vatom_radius[i]=0.147; vatom_radius_covalent[i]=1.60;   vatom_electronegativity[i]=1.54; vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=0.79; vatom_xray_scatt[i]=22.24; vatom_pearson_coefficient[i]=0.000286456;  // Ti
  i++; vatom_symbol[i]="V"; vatom_name[i]="Vanadium";    vatom_mass[i]*=50.9415;  vatom_volume[i]=13.2086; vatom_valence_std[i]=5; vatom_valence_iupac[i]=5; vatom_miedema_phi_star[i]=4.25; vatom_miedema_nws[i]=1.64; vatom_miedema_Vm[i]=4.1; vatom_miedema_gamma_s[i]=2600;  vatom_miedema_BVm[i]=14.0; vatom_radius[i]=0.132; vatom_radius_covalent[i]=1.53;   vatom_electronegativity[i]=1.63; vatom_crystal[i]="bcc";  vatom_pettifor_scale[i]=0.84; vatom_pearson_coefficient[i]=9.54831E-07;  // V
  i++; vatom_symbol[i]="Cr";vatom_name[i]="Chromium";    vatom_mass[i]*=51.996;   vatom_volume[i]=11.4136; vatom_valence_std[i]=6; vatom_valence_iupac[i]=6; vatom_miedema_phi_star[i]=4.65; vatom_miedema_nws[i]=1.74; vatom_miedema_Vm[i]=3.7; vatom_miedema_gamma_s[i]=2400;  vatom_miedema_BVm[i]=14.0; vatom_radius[i]=0.125; vatom_radius_covalent[i]=1.39;   vatom_electronegativity[i]=1.66; vatom_crystal[i]="bcc";  vatom_pettifor_scale[i]=0.89; vatom_xray_scatt[i]=23.84; vatom_pearson_coefficient[i]=0.00013287;  // Cr
  i++; vatom_symbol[i]="Mn";vatom_name[i]="Manganese";   vatom_mass[i]*=54.93805; vatom_volume[i]=10.6487; vatom_valence_std[i]=7; vatom_valence_iupac[i]=7; vatom_miedema_phi_star[i]=4.45; vatom_miedema_nws[i]=1.61; vatom_miedema_Vm[i]=3.8; vatom_miedema_gamma_s[i]=1600;  vatom_miedema_BVm[i]=4.4;  vatom_radius[i]=0.112; vatom_radius_covalent[i]=1.61;   vatom_electronegativity[i]=1.55; vatom_crystal[i]="cub";  vatom_pettifor_scale[i]=0.94; vatom_xray_scatt[i]=24.46; vatom_pearson_coefficient[i]=1.67276E-32;  // vatom_xray_scatt[i]=24.3589; Mn JX CHANGED VALENCE //DX+CO20170904 vatom_radius_covalent[i] uses high spin configuration (most frequent)   
  i++; vatom_symbol[i]="Fe";vatom_name[i]="Iron";        vatom_mass[i]*=55.847;   vatom_volume[i]=10.2315; vatom_valence_std[i]=8; vatom_valence_iupac[i]=6; vatom_miedema_phi_star[i]=4.93; vatom_miedema_nws[i]=1.77; vatom_miedema_Vm[i]=3.7; vatom_miedema_gamma_s[i]=2550;  vatom_miedema_BVm[i]=12.0; vatom_radius[i]=0.124; vatom_radius_covalent[i]=1.52;   vatom_electronegativity[i]=1.83; vatom_crystal[i]="bcc";  vatom_pettifor_scale[i]=0.99; vatom_xray_scatt[i]=24.85; vatom_pearson_coefficient[i]=9.17912E-05;  // vatom_xray_scatt[i]=24.6830; Fe JX CHANGED VALENCE //DX+CO20170904 vatom_radius_covalent[i] uses high spin configuration (most frequent) 
  i++; vatom_symbol[i]="Co";vatom_name[i]="Cobalt";      vatom_mass[i]*=58.9332;  vatom_volume[i]=10.3205; vatom_valence_std[i]=9; vatom_valence_iupac[i]=5; vatom_miedema_phi_star[i]=5.10; vatom_miedema_nws[i]=1.75; vatom_miedema_Vm[i]=3.5; vatom_miedema_gamma_s[i]=2550;  vatom_miedema_BVm[i]=13.0; vatom_radius[i]=0.125; vatom_radius_covalent[i]=1.26;   vatom_electronegativity[i]=1.88; vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=1.04; vatom_xray_scatt[i]=24.59; vatom_pearson_coefficient[i]=0.0;  // Co JX CHANGED VALENCE //DX+CO20170904 vatom_radius_covalent[i] uses low spin configuration (most frequent)   
  i++; vatom_symbol[i]="Ni";vatom_name[i]="Nickel";      vatom_mass[i]*=58.69;    vatom_volume[i]=10.8664; vatom_valence_std[i]=10; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=5.20; vatom_miedema_nws[i]=1.75; vatom_miedema_Vm[i]=3.5; vatom_miedema_gamma_s[i]=2450;  vatom_miedema_BVm[i]=12.0; vatom_radius[i]=0.125; vatom_radius_covalent[i]=1.24;   vatom_electronegativity[i]=1.91; vatom_crystal[i]="fcc";  vatom_pettifor_scale[i]=1.09; vatom_xray_scatt[i]=25.02; vatom_pearson_coefficient[i]=0.000430773;  // Ni
  i++; vatom_symbol[i]="Cu";vatom_name[i]="Copper";      vatom_mass[i]*=63.546;   vatom_volume[i]=12.0159; vatom_valence_std[i]=11; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=4.55; vatom_miedema_nws[i]=1.47; vatom_miedema_Vm[i]=3.7; vatom_miedema_gamma_s[i]=1850;  vatom_miedema_BVm[i]=9.3;  vatom_radius[i]=0.128; vatom_radius_covalent[i]=1.32;   vatom_electronegativity[i]=1.90; vatom_crystal[i]="fcc";  vatom_pettifor_scale[i]=1.20; vatom_xray_scatt[i]=27.03; vatom_pearson_coefficient[i]=0.00021086;  // Cu JX CHANGED VALENCE    
  i++; vatom_symbol[i]="Zn";vatom_name[i]="Zinc";        vatom_mass[i]*=65.38;    vatom_volume[i]=15.0827; vatom_valence_std[i]=12; vatom_valence_iupac[i]=2; vatom_miedema_phi_star[i]=4.10; vatom_miedema_nws[i]=1.32; vatom_miedema_Vm[i]=4.4; vatom_miedema_gamma_s[i]=1020;  vatom_miedema_BVm[i]=5.5;  vatom_radius[i]=0.133; vatom_radius_covalent[i]=1.22;   vatom_electronegativity[i]=1.65; vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=1.44; vatom_xray_scatt[i]=28.44; vatom_pearson_coefficient[i]=0.000595597;  // Zn
  // p-electron systems
  i++; vatom_symbol[i]="Ga";vatom_name[i]="Gallium";     vatom_mass[i]*=69.737;   vatom_volume[i]=18.9039; vatom_valence_std[i]=3; vatom_valence_iupac[i]=3; vatom_miedema_phi_star[i]=4.10; vatom_miedema_nws[i]=1.31; vatom_miedema_Vm[i]=5.2; vatom_miedema_gamma_s[i]= 830;  vatom_miedema_BVm[i]=6.7;  vatom_radius[i]=0.135; vatom_radius_covalent[i]=1.22;   vatom_electronegativity[i]=1.81; vatom_crystal[i]="orc";  vatom_pettifor_scale[i]=1.68; vatom_pearson_coefficient[i]=0.000197588;  // Ga
  i++; vatom_symbol[i]="Ge";vatom_name[i]="Germanium";   vatom_mass[i]*=72.59;    vatom_volume[i]=19.2948; vatom_valence_std[i]=4; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=4.55; vatom_miedema_nws[i]=1.37; vatom_miedema_Vm[i]=4.6; vatom_miedema_gamma_s[i]=1030;  vatom_miedema_BVm[i]=10.5; vatom_radius[i]=0.122; vatom_radius_covalent[i]=1.20;   vatom_electronegativity[i]=2.01; vatom_crystal[i]="dia";  vatom_pettifor_scale[i]=1.92; vatom_pearson_coefficient[i]=0.00058782;  // Ge
  i++; vatom_symbol[i]="As";vatom_name[i]="Arsenic";     vatom_mass[i]*=74.9216;  vatom_volume[i]=19.0677; vatom_valence_std[i]=5; vatom_valence_iupac[i]=5; vatom_miedema_phi_star[i]=4.80; vatom_miedema_nws[i]=1.44; vatom_miedema_Vm[i]=5.2; vatom_miedema_gamma_s[i]=1000;  vatom_miedema_BVm[i]=5.1;  vatom_radius[i]=0.125; vatom_radius_covalent[i]=1.19;   vatom_electronegativity[i]=2.18; vatom_crystal[i]="rhl";  vatom_pettifor_scale[i]=2.16; vatom_pearson_coefficient[i]=0.0;  // As
  i++; vatom_symbol[i]="Se";vatom_name[i]="Selenium";    vatom_mass[i]*=78.96;    vatom_volume[i]=20.3733; vatom_valence_std[i]=6; vatom_valence_iupac[i]=6; vatom_miedema_phi_star[i]=5.17; vatom_miedema_nws[i]=1.40; vatom_miedema_Vm[i]=5.172;vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.116; vatom_radius_covalent[i]=1.20;   vatom_electronegativity[i]=2.55; vatom_crystal[i]="hex";  vatom_pettifor_scale[i]=2.40; vatom_pearson_coefficient[i]=0.00046279;  // Se Table 27 of JX
  i++; vatom_symbol[i]="Br";vatom_name[i]="Bromine";     vatom_mass[i]*=79.904;   vatom_volume[i]=26.3292; vatom_valence_std[i]=7; vatom_valence_iupac[i]=7; vatom_miedema_phi_star[i]=5.20; vatom_miedema_nws[i]=1.35; vatom_miedema_Vm[i]=7.31; vatom_miedema_gamma_s[i]= 943;  vatom_miedema_BVm[i]=3.4;  vatom_radius[i]=0.119; vatom_radius_covalent[i]=1.20;   vatom_electronegativity[i]=2.96; vatom_crystal[i]="orc";  vatom_pettifor_scale[i]=2.64; vatom_pearson_coefficient[i]=0.000156277;  // Br interpolation phi_star, nws, Vm, gamma, BVm JX CHANGED VALENCE    
  i++; vatom_symbol[i]="Kr";vatom_name[i]="Krypton";     vatom_mass[i]*=83.8;     vatom_volume[i]=-1.0000; vatom_valence_std[i]=0; vatom_valence_iupac[i]=2; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.197; vatom_radius_covalent[i]=1.16;   vatom_electronegativity[i]=3; vatom_crystal[i]="fcc";  vatom_pettifor_scale[i]=0; vatom_pearson_coefficient[i]=0.000248482;  // Kr

  // ROW5
  // s-electron systems
  i++; vatom_symbol[i]="Rb";vatom_name[i]="Rubidium";    vatom_mass[i]*=85.4678;  vatom_volume[i]=91.2738; vatom_valence_std[i]=1; vatom_valence_iupac[i]=1; vatom_miedema_phi_star[i]=2.10; vatom_miedema_nws[i]=0.60; vatom_miedema_Vm[i]=14.6;vatom_miedema_gamma_s[i]= 120;  vatom_miedema_BVm[i]=1.8;  vatom_radius[i]=0.251; vatom_radius_covalent[i]=2.20;   vatom_electronegativity[i]=0.82; vatom_crystal[i]="bcc";  vatom_pettifor_scale[i]=0.30; vatom_pearson_coefficient[i]=0.000109697;  // Rb
  i++; vatom_symbol[i]="Sr";vatom_name[i]="Strontium";   vatom_mass[i]*=87.62;    vatom_volume[i]=55.4105; vatom_valence_std[i]=2; vatom_valence_iupac[i]=2; vatom_miedema_phi_star[i]=2.40; vatom_miedema_nws[i]=0.84; vatom_miedema_Vm[i]=10.2;vatom_miedema_gamma_s[i]= 430;  vatom_miedema_BVm[i]=3.9;  vatom_radius[i]=0.215; vatom_radius_covalent[i]=1.95;   vatom_electronegativity[i]=0.95; vatom_crystal[i]="fcc";  vatom_pettifor_scale[i]=0.55; vatom_pearson_coefficient[i]=6.09969E-05;  // Sr
  // d-electron systems: transition metals
  i++; vatom_symbol[i]="Y"; vatom_name[i]="Yttrium";     vatom_mass[i]*=88.9059;  vatom_volume[i]=32.4546; vatom_valence_std[i]=3; vatom_valence_iupac[i]=3; vatom_miedema_phi_star[i]=3.20; vatom_miedema_nws[i]=1.21; vatom_miedema_Vm[i]=7.3; vatom_miedema_gamma_s[i]=1100;  vatom_miedema_BVm[i]=7.2;  vatom_radius[i]=0.181; vatom_radius_covalent[i]=1.90;   vatom_electronegativity[i]=1.22; vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=0.70; vatom_pearson_coefficient[i]=0.0;  // Y
  i++; vatom_symbol[i]="Zr";vatom_name[i]="Zirconium";   vatom_mass[i]*=91.22;    vatom_volume[i]=23.2561; vatom_valence_std[i]=4; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=3.40; vatom_miedema_nws[i]=1.39; vatom_miedema_Vm[i]=5.8; vatom_miedema_gamma_s[i]=1950;  vatom_miedema_BVm[i]=12.0; vatom_radius[i]=0.158; vatom_radius_covalent[i]=1.75;   vatom_electronegativity[i]=1.33; vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=0.76; vatom_pearson_coefficient[i]=0.000342629;  // Zr
  i++; vatom_symbol[i]="Nb";vatom_name[i]="Niobium";     vatom_mass[i]*=92.9064;  vatom_volume[i]=18.3132; vatom_valence_std[i]=5; vatom_valence_iupac[i]=5; vatom_miedema_phi_star[i]=4.00; vatom_miedema_nws[i]=1.62; vatom_miedema_Vm[i]=4.9; vatom_miedema_gamma_s[i]=2700;  vatom_miedema_BVm[i]=18.0; vatom_radius[i]=0.143; vatom_radius_covalent[i]=1.64;   vatom_electronegativity[i]=1.60; vatom_crystal[i]="bcc";  vatom_pettifor_scale[i]=0.82; vatom_pearson_coefficient[i]=0.0;  // Nb
  i++; vatom_symbol[i]="Mo";vatom_name[i]="Molybdenum";  vatom_mass[i]*=95.94;    vatom_volume[i]=15.6175; vatom_valence_std[i]=6; vatom_valence_iupac[i]=6; vatom_miedema_phi_star[i]=4.65; vatom_miedema_nws[i]=1.77; vatom_miedema_Vm[i]=4.4; vatom_miedema_gamma_s[i]=2950;  vatom_miedema_BVm[i]=26.0; vatom_radius[i]=0.136; vatom_radius_covalent[i]=1.54;   vatom_electronegativity[i]=2.16; vatom_crystal[i]="bcc";  vatom_pettifor_scale[i]=0.88; vatom_pearson_coefficient[i]=0.000598128;  // Mo
  i++; vatom_symbol[i]="Tc";vatom_name[i]="Technetium";  vatom_mass[i]*=98.9062;  vatom_volume[i]=14.4670; vatom_valence_std[i]=7; vatom_valence_iupac[i]=7; vatom_miedema_phi_star[i]=5.30; vatom_miedema_nws[i]=1.81; vatom_miedema_Vm[i]=4.2; vatom_miedema_gamma_s[i]=3050;  vatom_miedema_BVm[i]=26.0; vatom_radius[i]=NNN;   vatom_radius_covalent[i]=1.47;   vatom_electronegativity[i]=1.90; vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=0.94; vatom_pearson_coefficient[i]=0.0;  //Tc JX CHANGED VALENCE    
  i++; vatom_symbol[i]="Ru";vatom_name[i]="Ruthenium";   vatom_mass[i]*=101.07;   vatom_volume[i]=13.8390; vatom_valence_std[i]=8; vatom_valence_iupac[i]=8; vatom_miedema_phi_star[i]=5.40; vatom_miedema_nws[i]=1.83; vatom_miedema_Vm[i]=4.1; vatom_miedema_gamma_s[i]=3050;  vatom_miedema_BVm[i]=26.0; vatom_radius[i]=0.134; vatom_radius_covalent[i]=1.46;   vatom_electronegativity[i]=2.20; vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=1.00; vatom_pearson_coefficient[i]=0.000406665;  //Ru JX CHANGED VALENCE    
  i++; vatom_symbol[i]="Rh";vatom_name[i]="Rhodium";     vatom_mass[i]*=102.9055; vatom_volume[i]=14.1731; vatom_valence_std[i]=9; vatom_valence_iupac[i]=6; vatom_miedema_phi_star[i]=5.40; vatom_miedema_nws[i]=1.76; vatom_miedema_Vm[i]=4.1; vatom_miedema_gamma_s[i]=2750;  vatom_miedema_BVm[i]=23.0; vatom_radius[i]=0.134; vatom_radius_covalent[i]=1.42;   vatom_electronegativity[i]=2.28; vatom_crystal[i]="fcc";  vatom_pettifor_scale[i]=1.06; vatom_pearson_coefficient[i]=1.90706E-32;  // Rh
  i++; vatom_symbol[i]="Pd";vatom_name[i]="Palladium";   vatom_mass[i]*=106.4;    vatom_volume[i]=15.4596; vatom_valence_std[i]=10; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=5.45; vatom_miedema_nws[i]=1.67; vatom_miedema_Vm[i]=4.3; vatom_miedema_gamma_s[i]=2100;  vatom_miedema_BVm[i]=16.0; vatom_radius[i]=0.137; vatom_radius_covalent[i]=1.39;   vatom_electronegativity[i]=2.20; vatom_crystal[i]="fcc";  vatom_pettifor_scale[i]=1.12; vatom_pearson_coefficient[i]=0.000309478;  // Pd
  i++; vatom_symbol[i]="Ag";vatom_name[i]="Silver";      vatom_mass[i]*=107.8682; vatom_volume[i]=18.0678; vatom_valence_std[i]=11; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=4.45; vatom_miedema_nws[i]=1.39; vatom_miedema_Vm[i]=4.7; vatom_miedema_gamma_s[i]=1250;  vatom_miedema_BVm[i]=10.0; vatom_radius[i]=0.144; vatom_radius_covalent[i]=1.45;   vatom_electronegativity[i]=1.93; vatom_crystal[i]="fcc";  vatom_pettifor_scale[i]=1.18; vatom_xray_scatt[i]=47.18; vatom_pearson_coefficient[i]=8.57985E-05;  // Ag JX CHANGED VALENCE    
  i++; vatom_symbol[i]="Cd";vatom_name[i]="Cadmium";     vatom_mass[i]*=112.41;   vatom_volume[i]=22.0408; vatom_valence_std[i]=12; vatom_valence_iupac[i]=2; vatom_miedema_phi_star[i]=4.05; vatom_miedema_nws[i]=1.24; vatom_miedema_Vm[i]=5.5; vatom_miedema_gamma_s[i]= 780;  vatom_miedema_BVm[i]=6.10; vatom_radius[i]=0.150; vatom_radius_covalent[i]=1.44;   vatom_electronegativity[i]=1.69; vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=1.36; vatom_pearson_coefficient[i]=0.000271603;  // Cd
  // p-electron systems
  i++; vatom_symbol[i]="In";vatom_name[i]="Indium";      vatom_mass[i]*=114.82;   vatom_volume[i]=27.5233; vatom_valence_std[i]=3; vatom_valence_iupac[i]=3; vatom_miedema_phi_star[i]=3.90; vatom_miedema_nws[i]=1.17; vatom_miedema_Vm[i]=6.3; vatom_miedema_gamma_s[i]= 690;  vatom_miedema_BVm[i]=6.4;  vatom_radius[i]=0.157; vatom_radius_covalent[i]=1.42;   vatom_electronegativity[i]=1.78; vatom_crystal[i]="fct";  vatom_pettifor_scale[i]=1.60; vatom_pearson_coefficient[i]=1.24494E-05;  // In
  i++; vatom_symbol[i]="Sn";vatom_name[i]="Tin";         vatom_mass[i]*=118.69;   vatom_volume[i]=27.5555; vatom_valence_std[i]=4; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=4.15; vatom_miedema_nws[i]=1.24; vatom_miedema_Vm[i]=6.4; vatom_miedema_gamma_s[i]= 710;  vatom_miedema_BVm[i]=8.8;  vatom_radius[i]=0.158; vatom_radius_covalent[i]=1.39;   vatom_electronegativity[i]=1.96; vatom_crystal[i]="bct";  vatom_pettifor_scale[i]=1.84; vatom_pearson_coefficient[i]=0.000334085;  // Sn
  i++; vatom_symbol[i]="Sb";vatom_name[i]="Antimony";    vatom_mass[i]*=121.75;   vatom_volume[i]=27.1823; vatom_valence_std[i]=5; vatom_valence_iupac[i]=5; vatom_miedema_phi_star[i]=4.40; vatom_miedema_nws[i]=1.26; vatom_miedema_Vm[i]=6.6; vatom_miedema_gamma_s[i]= 680;  vatom_miedema_BVm[i]=7.0;  vatom_radius[i]=0.161; vatom_radius_covalent[i]=1.39;   vatom_electronegativity[i]=2.05; vatom_crystal[i]="rhl";  vatom_pettifor_scale[i]=2.08; vatom_pearson_coefficient[i]=6.60751E-05;  // Sb
  i++; vatom_symbol[i]="Te";vatom_name[i]="Tellurium";   vatom_mass[i]*=127.6;    vatom_volume[i]=28.1993; vatom_valence_std[i]=6; vatom_valence_iupac[i]=6; vatom_miedema_phi_star[i]=4.72; vatom_miedema_nws[i]=1.31;vatom_miedema_Vm[i]=6.439;vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.143; vatom_radius_covalent[i]=1.38;   vatom_electronegativity[i]=2.10; vatom_crystal[i]="hex";  vatom_pettifor_scale[i]=2.32; vatom_pearson_coefficient[i]=0.000283934;  // Te Table 27 of JX
  i++; vatom_symbol[i]="I"; vatom_name[i]="Iodine";     vatom_mass[i]*=126.9045; vatom_volume[i]=34.9784;  vatom_valence_std[i]=7; vatom_valence_iupac[i]=7; vatom_miedema_phi_star[i]=5.33; vatom_miedema_nws[i]=0.17;  vatom_miedema_Vm[i]=8.72;vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;vatom_radius[i]=0.136; vatom_radius_covalent[i]=1.39;   vatom_electronegativity[i]=2.66; vatom_crystal[i]="orc";  vatom_pettifor_scale[i]=2.56; vatom_pearson_coefficient[i]=0.0;  // I interpolation phi_star, nws, Vm,
  i++; vatom_symbol[i]="Xe";vatom_name[i]="Xenon";       vatom_mass[i]*=131.3;    vatom_volume[i]=-1.0000; vatom_valence_std[i]=0; vatom_valence_iupac[i]=8; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.218; vatom_radius_covalent[i]=1.40;   vatom_electronegativity[i]=2.60; vatom_crystal[i]="fcc";  vatom_pettifor_scale[i]=0; vatom_pearson_coefficient[i]=0.000267781;  //Xe JX CHANGED VALENCE    

  // ROW6
  // s-electron systems
  i++; vatom_symbol[i]="Cs";vatom_name[i]="Cesium";      vatom_mass[i]*=132.9054; vatom_volume[i]=117.281; vatom_valence_std[i]=1; vatom_valence_iupac[i]=1; vatom_miedema_phi_star[i]=1.95; vatom_miedema_nws[i]=0.55; vatom_miedema_Vm[i]=16.8;vatom_miedema_gamma_s[i]=  95;  vatom_miedema_BVm[i]=1.4;  vatom_radius[i]=0.265; vatom_radius_covalent[i]=2.44;   vatom_electronegativity[i]=0.79; vatom_crystal[i]="bcc";  vatom_pettifor_scale[i]=0.25; vatom_pearson_coefficient[i]=0.0;  // Cs
  i++; vatom_symbol[i]="Ba";vatom_name[i]="Barium";      vatom_mass[i]*=137.33;   vatom_volume[i]=62.6649; vatom_valence_std[i]=2; vatom_valence_iupac[i]=2; vatom_miedema_phi_star[i]=2.32; vatom_miedema_nws[i]=0.81; vatom_miedema_Vm[i]=11.3;vatom_miedema_gamma_s[i]= 370;  vatom_miedema_BVm[i]=3.9;  vatom_radius[i]=0.217; vatom_radius_covalent[i]=2.15;   vatom_electronegativity[i]=0.89; vatom_crystal[i]="bcc";  vatom_pettifor_scale[i]=0.50; vatom_pearson_coefficient[i]=6.23705E-05;  // Ba
  // d-electron systems: transition metals
  i++; vatom_symbol[i]="La";vatom_name[i]="Lanthanium";  vatom_mass[i]*=138.9055; vatom_volume[i]=36.8495; vatom_valence_std[i]=3; vatom_valence_iupac[i]=3; vatom_miedema_phi_star[i]=3.05;/*3.17*/ vatom_miedema_nws[i]=1.09;/*1.18*/ vatom_miedema_Vm[i]=8.0;/*7.98*/ vatom_miedema_gamma_s[i]= 900;  vatom_miedema_BVm[i]=5.5;  vatom_radius[i]=0.187; vatom_radius_covalent[i]=2.07;   vatom_electronegativity[i]=1.10; vatom_crystal[i]="hex";  vatom_pettifor_scale[i]=0.7480; vatom_pearson_coefficient[i]=4.65323E-08;  // La
  // lantanidies
  i++; vatom_symbol[i]="Ce";vatom_name[i]="Cerium";      vatom_mass[i]*=140.12;   vatom_volume[i]=26.4729; vatom_valence_std[i]=4; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=3.18;  vatom_miedema_nws[i]=1.19;  vatom_miedema_Vm[i]=7.76; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.182; vatom_radius_covalent[i]=2.04;   vatom_electronegativity[i]=1.12; vatom_crystal[i]="fcc";  vatom_pettifor_scale[i]=0.7460; vatom_pearson_coefficient[i]=2.24956E-05;  // Ce pettifor linear interpolation // miedema from Alonso-March.
  i++; vatom_symbol[i]="Pr";vatom_name[i]="Praseodymium";vatom_mass[i]*=140.9077; vatom_volume[i]=36.4987; vatom_valence_std[i]=5; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=3.19;  vatom_miedema_nws[i]=1.20;  vatom_miedema_Vm[i]=7.56; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.183; vatom_radius_covalent[i]=2.03;   vatom_electronegativity[i]=1.13; vatom_crystal[i]="hex";  vatom_pettifor_scale[i]=0.7440; vatom_pearson_coefficient[i]=0.0;  // Pr pettifor linear interpolation
  i++; vatom_symbol[i]="Nd";vatom_name[i]="Neodymium";   vatom_mass[i]*=144.24;   vatom_volume[i]=29.6719; vatom_valence_std[i]=6; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=3.19;  vatom_miedema_nws[i]=1.20;  vatom_miedema_Vm[i]=7.51; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.182; vatom_radius_covalent[i]=2.01;   vatom_electronegativity[i]=1.14; vatom_crystal[i]="hex";  vatom_pettifor_scale[i]=0.7420; vatom_pearson_coefficient[i]=0.000231599;  // Nd pettifor linear interpolation JX CHANGED VALENCE    
  i++; vatom_symbol[i]="Pm";vatom_name[i]="Promethium";  vatom_mass[i]*=146.92;   vatom_volume[i]=34.6133; vatom_valence_std[i]=7; vatom_valence_iupac[i]=3; vatom_miedema_phi_star[i]=3.19;  vatom_miedema_nws[i]=1.21;  vatom_miedema_Vm[i]=7.43; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=NNN;   vatom_radius_covalent[i]=1.99;   vatom_electronegativity[i]=1.13; vatom_crystal[i]="hex";  vatom_pettifor_scale[i]=0.7400; vatom_pearson_coefficient[i]=0.0;  // Pm pettifor linear interpolation
  i++; vatom_symbol[i]="Sm";vatom_name[i]="Samarium";    vatom_mass[i]*=150.4;    vatom_volume[i]=33.9484; vatom_valence_std[i]=8; vatom_valence_iupac[i]=3; vatom_miedema_phi_star[i]=3.20;  vatom_miedema_nws[i]=1.21;  vatom_miedema_Vm[i]=7.37; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.181; vatom_radius_covalent[i]=1.98;   vatom_electronegativity[i]=1.17; vatom_crystal[i]="rhl";  vatom_pettifor_scale[i]=0.7380; vatom_pearson_coefficient[i]=0.000334686;  // Sm pettifor linear interpolation
  i++; vatom_symbol[i]="Eu";vatom_name[i]="Europium";    vatom_mass[i]*=151.96;   vatom_volume[i]=43.1719; vatom_valence_std[i]=9; vatom_valence_iupac[i]=3; vatom_miedema_phi_star[i]=3.20;  vatom_miedema_nws[i]=1.21;  vatom_miedema_Vm[i]=7.36; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.204; vatom_radius_covalent[i]=1.98;   vatom_electronegativity[i]=1.20; vatom_crystal[i]="bcc";  vatom_pettifor_scale[i]=0.7360; vatom_pearson_coefficient[i]=4.32857E-05;  // Eu pettifor linear interpolation
  i++; vatom_symbol[i]="Gd";vatom_name[i]="Gadolinium";  vatom_mass[i]*=157.25;   vatom_volume[i]=32.5777; vatom_valence_std[i]=10; vatom_valence_iupac[i]=3; vatom_miedema_phi_star[i]=3.20;  vatom_miedema_nws[i]=1.21;  vatom_miedema_Vm[i]=7.34; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.180; vatom_radius_covalent[i]=1.96;   vatom_electronegativity[i]=1.20; vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=0.7340; vatom_pearson_coefficient[i]=0.000127674;  // Gd pettifor linear interpolation
  i++; vatom_symbol[i]="Tb";vatom_name[i]="Terbium";     vatom_mass[i]*=158.9254; vatom_volume[i]=32.0200; vatom_valence_std[i]=11; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=3.21;  vatom_miedema_nws[i]=1.22;  vatom_miedema_Vm[i]=7.20; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.177; vatom_radius_covalent[i]=1.94;   vatom_electronegativity[i]=1.10; vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=0.7320; vatom_pearson_coefficient[i]=0.0;  // Tb pettifor linear interpolation
  i++; vatom_symbol[i]="Dy";vatom_name[i]="Dysprosium";  vatom_mass[i]*=162.5;    vatom_volume[i]=31.5096; vatom_valence_std[i]=12; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=3.21;  vatom_miedema_nws[i]=1.22;  vatom_miedema_Vm[i]=7.12; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.177; vatom_radius_covalent[i]=1.92;   vatom_electronegativity[i]=1.22; vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=0.7300; vatom_pearson_coefficient[i]=5.20771E-05;  // Dy pettifor linear interpolation JX CHANGED VALENCE    
  i++; vatom_symbol[i]="Ho";vatom_name[i]="Holmium";     vatom_mass[i]*=164.9304; vatom_volume[i]=31.0155; vatom_valence_std[i]=13; vatom_valence_iupac[i]=3; vatom_miedema_phi_star[i]=3.22;  vatom_miedema_nws[i]=1.22;  vatom_miedema_Vm[i]=7.06; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.176; vatom_radius_covalent[i]=1.92;   vatom_electronegativity[i]=1.23; vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=0.7280; vatom_pearson_coefficient[i]=2.96961E-32;  // Ho pettifor linear interpolation
  i++; vatom_symbol[i]="Er";vatom_name[i]="Erbium";      vatom_mass[i]*=167.26;   vatom_volume[i]=30.5431; vatom_valence_std[i]=14; vatom_valence_iupac[i]=3; vatom_miedema_phi_star[i]=3.22;  vatom_miedema_nws[i]=1.23;  vatom_miedema_Vm[i]=6.98; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.175; vatom_radius_covalent[i]=1.89;   vatom_electronegativity[i]=1.24; vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=0.7260; vatom_pearson_coefficient[i]=7.24618E-05;  // Er pettifor linear interpolation
  i++; vatom_symbol[i]="Tm";vatom_name[i]="Thulium";     vatom_mass[i]*=168.9342; vatom_volume[i]=30.0016; vatom_valence_std[i]=15; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=3.22;  vatom_miedema_nws[i]=1.23;  vatom_miedema_Vm[i]=6.90; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.174; vatom_radius_covalent[i]=1.90;   vatom_electronegativity[i]=1.25; vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=0.7240; vatom_pearson_coefficient[i]=0.0;  // Tm pettifor linear interpolation JX CHANGED VALENCE    
  i++; vatom_symbol[i]="Yb";vatom_name[i]="Ytterbium";   vatom_mass[i]*=173.04;   vatom_volume[i]=39.4395; vatom_valence_std[i]=16; vatom_valence_iupac[i]=3; vatom_miedema_phi_star[i]=3.22;  vatom_miedema_nws[i]=1.23;  vatom_miedema_Vm[i]=6.86; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.193; vatom_radius_covalent[i]=1.87;   vatom_electronegativity[i]=1.10; vatom_crystal[i]="fcc";  vatom_pettifor_scale[i]=0.7220; vatom_pearson_coefficient[i]=8.54557E-05;  // Yb pettifor linear interpolation
  i++; vatom_symbol[i]="Lu";vatom_name[i]="Lutetium";    vatom_mass[i]*=174.967;  vatom_volume[i]=29.3515; vatom_valence_std[i]=17; vatom_valence_iupac[i]=3; vatom_miedema_phi_star[i]=3.22;  vatom_miedema_nws[i]=1.24;  vatom_miedema_Vm[i]=6.81; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.173; vatom_radius_covalent[i]=1.87;   vatom_electronegativity[i]=1.27; vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=0.7200; vatom_pearson_coefficient[i]=8.27273E-07;  // Lu
  // d-electron systems: transition metalsnnn";
  i++; vatom_symbol[i]="Hf";vatom_name[i]="Hafnium";     vatom_mass[i]*=178.49;   vatom_volume[i]=22.0408; vatom_valence_std[i]=4; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=3.55; vatom_miedema_nws[i]=1.43; vatom_miedema_Vm[i]=5.6; vatom_miedema_gamma_s[i]=2200;  vatom_miedema_BVm[i]=15.0; vatom_radius[i]=0.159; vatom_radius_covalent[i]=1.75;   vatom_electronegativity[i]=1.30; vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=0.775; vatom_pearson_coefficient[i]=5.25384E-05;  // Hf
  i++; vatom_symbol[i]="Ta";vatom_name[i]="Tantalum";    vatom_mass[i]*=180.9479; vatom_volume[i]=18.1100; vatom_valence_std[i]=5; vatom_valence_iupac[i]=5; vatom_miedema_phi_star[i]=4.05; vatom_miedema_nws[i]=1.63; vatom_miedema_Vm[i]=4.9; vatom_miedema_gamma_s[i]=3050;  vatom_miedema_BVm[i]=22.0; vatom_radius[i]=0.147; vatom_radius_covalent[i]=1.70;   vatom_electronegativity[i]=1.50; vatom_crystal[i]="bcc";  vatom_pettifor_scale[i]=0.83; vatom_pearson_coefficient[i]=3.66845E-09;  // Ta
  i++; vatom_symbol[i]="W"; vatom_name[i]="Tungsten";    vatom_mass[i]*=183.85;   vatom_volume[i]=15.9387; vatom_valence_std[i]=6; vatom_valence_iupac[i]=6; vatom_miedema_phi_star[i]=4.80; vatom_miedema_nws[i]=1.81; vatom_miedema_Vm[i]=4.5; vatom_miedema_gamma_s[i]=3300;  vatom_miedema_BVm[i]=31.0; vatom_radius[i]=0.137; vatom_radius_covalent[i]=1.62;   vatom_electronegativity[i]=2.36; vatom_crystal[i]="bcc";  vatom_pettifor_scale[i]=0.885; vatom_pearson_coefficient[i]=6.96679E-05;  // W
  i++; vatom_symbol[i]="Re";vatom_name[i]="Rhenium";     vatom_mass[i]*=186.2;    vatom_volume[i]=14.8941; vatom_valence_std[i]=7; vatom_valence_iupac[i]=7; vatom_miedema_phi_star[i]=5.40; vatom_miedema_nws[i]=1.86; vatom_miedema_Vm[i]=4.3; vatom_miedema_gamma_s[i]=3650;  vatom_miedema_BVm[i]=33.0; vatom_radius[i]=0.138; vatom_radius_covalent[i]=1.51;   vatom_electronegativity[i]=1.90; vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=0.94; vatom_pearson_coefficient[i]=2.70849E-05;  // Re
  i++; vatom_symbol[i]="Os";vatom_name[i]="Osmium";      vatom_mass[i]*=190.2;    vatom_volume[i]=14.2403; vatom_valence_std[i]=8; vatom_valence_iupac[i]=8; vatom_miedema_phi_star[i]=5.40; vatom_miedema_nws[i]=1.85; vatom_miedema_Vm[i]=4.2; vatom_miedema_gamma_s[i]=3500;  vatom_miedema_BVm[i]=35.0; vatom_radius[i]=0.135; vatom_radius_covalent[i]=1.44;   vatom_electronegativity[i]=2.20; vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=0.995; vatom_pearson_coefficient[i]=7.45234E-05;  // Os JX CHANGED VALENCE    
  i++; vatom_symbol[i]="Ir";vatom_name[i]="Iridium";     vatom_mass[i]*=192.22;   vatom_volume[i]=14.5561; vatom_valence_std[i]=9; vatom_valence_iupac[i]=8; vatom_miedema_phi_star[i]=5.55; vatom_miedema_nws[i]=1.83; vatom_miedema_Vm[i]=4.2; vatom_miedema_gamma_s[i]=3100;  vatom_miedema_BVm[i]=25.0; vatom_radius[i]=0.135; vatom_radius_covalent[i]=1.41;   vatom_electronegativity[i]=2.20; vatom_crystal[i]="fcc";  vatom_pettifor_scale[i]=1.05; vatom_pearson_coefficient[i]=2.53787E-05;  // Ir JX CHANGED VALENCE    
  i++; vatom_symbol[i]="Pt";vatom_name[i]="Platinum";    vatom_mass[i]*=195.09;   vatom_volume[i]=15.7298; vatom_valence_std[i]=10; vatom_valence_iupac[i]=6; vatom_miedema_phi_star[i]=5.65; vatom_miedema_nws[i]=1.78; vatom_miedema_Vm[i]=4.4; vatom_miedema_gamma_s[i]=2550;  vatom_miedema_BVm[i]=18.0; vatom_radius[i]=0.138; vatom_radius_covalent[i]=1.36;   vatom_electronegativity[i]=2.28; vatom_crystal[i]="fcc";  vatom_pettifor_scale[i]=1.105; vatom_pearson_coefficient[i]=3.39206E-05;  // Pt
  i++; vatom_symbol[i]="Au";vatom_name[i]="Gold";        vatom_mass[i]*=196.9665; vatom_volume[i]=18.1904; vatom_valence_std[i]=11; vatom_valence_iupac[i]=5; vatom_miedema_phi_star[i]=5.15; vatom_miedema_nws[i]=1.57; vatom_miedema_Vm[i]=4.7; vatom_miedema_gamma_s[i]=1550;  vatom_miedema_BVm[i]=18.0; vatom_radius[i]=0.144; vatom_radius_covalent[i]=1.36;   vatom_electronegativity[i]=2.54; vatom_crystal[i]="fcc";  vatom_pettifor_scale[i]=1.16; vatom_xray_scatt[i]=74.99; vatom_pearson_coefficient[i]=2.08217E-32;  // Au
  i++; vatom_symbol[i]="Hg";vatom_name[i]="Mercury";     vatom_mass[i]*=200.59;   vatom_volume[i]=29.7156; vatom_valence_std[i]=12; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=4.20; vatom_miedema_nws[i]=1.24; vatom_miedema_Vm[i]=5.8; vatom_miedema_gamma_s[i]= 610;  vatom_miedema_BVm[i]=4.0;  vatom_radius[i]=0.150; vatom_radius_covalent[i]=1.32;   vatom_electronegativity[i]=2.00; vatom_crystal[i]="rhl";  vatom_pettifor_scale[i]=1.32; vatom_pearson_coefficient[i]=6.52519E-05;  // Hg
  // p-electron systems
  i++; vatom_symbol[i]="Tl";vatom_name[i]="Thallium";    vatom_mass[i]*=204.37;   vatom_volume[i]=31.0721; vatom_valence_std[i]=3; vatom_valence_iupac[i]=3; vatom_miedema_phi_star[i]=3.90; vatom_miedema_nws[i]=1.12; vatom_miedema_Vm[i]=6.6; vatom_miedema_gamma_s[i]= 610;  vatom_miedema_BVm[i]=6.2;  vatom_radius[i]=0.171; vatom_radius_covalent[i]=1.45;   vatom_electronegativity[i]=1.62; vatom_crystal[i]="hcp";  vatom_pettifor_scale[i]=1.56; vatom_pearson_coefficient[i]=1.99659E-05;  // Tl electronegativity  2.04=>1.62
  i++; vatom_symbol[i]="Pb";vatom_name[i]="Lead";        vatom_mass[i]*=207.2;    vatom_volume[i]=31.6649; vatom_valence_std[i]=4; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=4.10; vatom_miedema_nws[i]=1.15; vatom_miedema_Vm[i]=6.9; vatom_miedema_gamma_s[i]= 610;  vatom_miedema_BVm[i]=7.9;  vatom_radius[i]=0.175; vatom_radius_covalent[i]=1.46;   vatom_electronegativity[i]=2.33; vatom_crystal[i]="fcc";  vatom_pettifor_scale[i]=1.80; vatom_pearson_coefficient[i]=1.94378E-05;  // Pb
  i++; vatom_symbol[i]="Bi";vatom_name[i]="Bismuth";     vatom_mass[i]*=208.9804; vatom_volume[i]=31.5691; vatom_valence_std[i]=5; vatom_valence_iupac[i]=5; vatom_miedema_phi_star[i]=4.15; vatom_miedema_nws[i]=1.16; vatom_miedema_Vm[i]=7.2; vatom_miedema_gamma_s[i]= 550;  vatom_miedema_BVm[i]=6.7;  vatom_radius[i]=0.182; vatom_radius_covalent[i]=1.48;   vatom_electronegativity[i]=2.02; vatom_crystal[i]="rhl";  vatom_pettifor_scale[i]=2.04; vatom_pearson_coefficient[i]=0.0;  // Bi
  i++; vatom_symbol[i]="Po";vatom_name[i]="Polonium";    vatom_mass[i]*=209.98;   vatom_volume[i]=NNN;     vatom_valence_std[i]=6; vatom_valence_iupac[i]=6; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.140; vatom_radius_covalent[i]=1.40;   vatom_electronegativity[i]=2.00; vatom_crystal[i]="sc";  vatom_pettifor_scale[i]=2.28; vatom_pearson_coefficient[i]=0.0;  // Po
  i++; vatom_symbol[i]="At";vatom_name[i]="Astatine";    vatom_mass[i]*=210;      vatom_volume[i]=NNN;     vatom_valence_std[i]=7; vatom_valence_iupac[i]=7; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=NNN;   vatom_radius_covalent[i]=1.50;   vatom_electronegativity[i]=2.20; vatom_crystal[i]="nnn";  vatom_pettifor_scale[i]=2.52; vatom_pearson_coefficient[i]=0.0;  // At
  i++; vatom_symbol[i]="Rn";vatom_name[i]="Radon";       vatom_mass[i]*=222;      vatom_volume[i]=NNN;     vatom_valence_std[i]=0; vatom_valence_iupac[i]=6; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=NNN;   vatom_radius_covalent[i]=1.50;   vatom_electronegativity[i]=2.2; vatom_crystal[i]="fcc";  vatom_pettifor_scale[i]=0; vatom_pearson_coefficient[i]=0.0;  // Rn

  // ROW7
  // s-electron systems
  i++; vatom_symbol[i]="Fr";vatom_name[i]="Francium";    vatom_mass[i]*=223.02;   vatom_volume[i]=NNN;     vatom_valence_std[i]=1; vatom_valence_iupac[i]=1; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=NNN;   vatom_radius_covalent[i]=2.60;   vatom_electronegativity[i]=0.70; vatom_crystal[i]="bcc";  vatom_pettifor_scale[i]=0; vatom_pearson_coefficient[i]=0.0;  // Fr
  i++; vatom_symbol[i]="Ra";vatom_name[i]="Radium";      vatom_mass[i]*=226.0254; vatom_volume[i]=-1.0000; vatom_valence_std[i]=2; vatom_valence_iupac[i]=2; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=NNN;   vatom_radius_covalent[i]=2.21;   vatom_electronegativity[i]=0.89; vatom_crystal[i]="bct";  vatom_pettifor_scale[i]=0; vatom_pearson_coefficient[i]=0.0;  // Ra
  // d-electron systems: transition metals
  i++; vatom_symbol[i]="Ac";vatom_name[i]="Actinium";    vatom_mass[i]*=227.03;   vatom_volume[i]=45.2437; vatom_valence_std[i]=3; vatom_valence_iupac[i]=3; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=NNN;   vatom_radius_covalent[i]=2.15;   vatom_electronegativity[i]=1.10; vatom_crystal[i]="fcc";  vatom_pettifor_scale[i]=0; vatom_pearson_coefficient[i]=0.0;  // Ac
  // actinidies
  i++; vatom_symbol[i]="Th";vatom_name[i]="Thorium";     vatom_mass[i]*=232.0381; vatom_volume[i]=31.9586; vatom_valence_std[i]=4; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=3.30; vatom_miedema_nws[i]=1.28; vatom_miedema_Vm[i]=7.3; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.180; vatom_radius_covalent[i]=2.06;   vatom_electronegativity[i]=1.30; vatom_crystal[i]="fcc";  vatom_pettifor_scale[i]=0; vatom_xray_scatt[i]=86.64; vatom_pearson_coefficient[i]=0.0; // Th
  i++; vatom_symbol[i]="Pa";vatom_name[i]="Protoactinium";vatom_mass[i]*=231.04;  vatom_volume[i]=NNN;     vatom_valence_std[i]=5; vatom_valence_iupac[i]=5; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=NNN;   vatom_radius_covalent[i]=2.00;   vatom_electronegativity[i]=1.50; vatom_crystal[i]="bct";  vatom_pettifor_scale[i]=0; vatom_pearson_coefficient[i]=0.0;  // Pa
  i++; vatom_symbol[i]="U"; vatom_name[i]="Uranium";     vatom_mass[i]*=238.03;   vatom_volume[i]=NNN;     vatom_valence_std[i]=6; vatom_valence_iupac[i]=6; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=0.138; vatom_radius_covalent[i]=1.96;   vatom_electronegativity[i]=1.38; vatom_crystal[i]="orc";  vatom_pettifor_scale[i]=0; vatom_pearson_coefficient[i]=1.15611E-06;  // U
  i++; vatom_symbol[i]="Np";vatom_name[i]="Neptunium";   vatom_mass[i]*=237.05;   vatom_volume[i]=NNN;     vatom_valence_std[i]=7; vatom_valence_iupac[i]=7; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=NNN;   vatom_radius_covalent[i]=1.90;   vatom_electronegativity[i]=NNN; vatom_crystal[i]="nnn";  vatom_pettifor_scale[i]=0;  vatom_pearson_coefficient[i]=0.0;  // Np
  i++; vatom_symbol[i]="Pu"; vatom_name[i]="Plutonium";  vatom_mass[i]*=244.06;   vatom_volume[i]=NNN;     vatom_valence_std[i]=8; vatom_valence_iupac[i]=7; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=NNN; vatom_radius_covalent[i]=1.87;   vatom_electronegativity[i]=NNN; vatom_crystal[i]="nnn";  vatom_pettifor_scale[i]=0; vatom_pearson_coefficient[i]=0.0;  // Pu
  i++; vatom_symbol[i]="Am";vatom_name[i]="Americium";   vatom_mass[i]*=243.06;   vatom_volume[i]=NNN;     vatom_valence_std[i]=9; vatom_valence_iupac[i]=7; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=NNN;   vatom_radius_covalent[i]=1.80;   vatom_electronegativity[i]=NNN; vatom_crystal[i]="nnn";  vatom_pettifor_scale[i]=0; vatom_pearson_coefficient[i]=0.0;  // Am
  i++; vatom_symbol[i]="Cm"; vatom_name[i]="Curium";     vatom_mass[i]*=247.07;   vatom_volume[i]=NNN;     vatom_valence_std[i]=10; vatom_valence_iupac[i]=8; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=NNN; vatom_radius_covalent[i]=1.69;   vatom_electronegativity[i]=NNN; vatom_crystal[i]="nnn";  vatom_pettifor_scale[i]=0; vatom_pearson_coefficient[i]=0.0;  // Cm

  i++; vatom_symbol[i]="Bk"; vatom_name[i]="Berkelium";  vatom_mass[i]*=247.07;   vatom_volume[i]=NNN;     vatom_valence_std[i]=11; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=NNN; vatom_radius_covalent[i]=NNN;   vatom_electronegativity[i]=NNN; vatom_crystal[i]="nnn";  vatom_pettifor_scale[i]=0; vatom_pearson_coefficient[i]=0.0;  // Bk
  i++; vatom_symbol[i]="Cf"; vatom_name[i]="Californium";vatom_mass[i]*=251.08;   vatom_volume[i]=NNN;     vatom_valence_std[i]=12; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=NNN; vatom_radius_covalent[i]=NNN;   vatom_electronegativity[i]=NNN; vatom_crystal[i]="nnn";  vatom_pettifor_scale[i]=0; vatom_pearson_coefficient[i]=0.0;  // Cf
  i++; vatom_symbol[i]="Es"; vatom_name[i]="Einsteinium";vatom_mass[i]*=252.08;   vatom_volume[i]=NNN;     vatom_valence_std[i]=13; vatom_valence_iupac[i]=4; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=NNN; vatom_radius_covalent[i]=NNN;   vatom_electronegativity[i]=NNN; vatom_crystal[i]="nnn";  vatom_pettifor_scale[i]=0; vatom_pearson_coefficient[i]=0.0;  // Es
  i++; vatom_symbol[i]="Fm"; vatom_name[i]="Fermium";    vatom_mass[i]*=257.1;    vatom_volume[i]=NNN;     vatom_valence_std[i]=14; vatom_valence_iupac[i]=3; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=NNN; vatom_radius_covalent[i]=NNN;   vatom_electronegativity[i]=NNN; vatom_crystal[i]="nnn";  vatom_pettifor_scale[i]=0; vatom_pearson_coefficient[i]=0.0;  // Fm
  i++; vatom_symbol[i]="Md"; vatom_name[i]="Mendelevium";vatom_mass[i]*=258.1;    vatom_volume[i]=NNN;     vatom_valence_std[i]=15; vatom_valence_iupac[i]=3; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=NNN; vatom_radius_covalent[i]=NNN;   vatom_electronegativity[i]=NNN; vatom_crystal[i]="nnn";  vatom_pettifor_scale[i]=0; vatom_pearson_coefficient[i]=0.0;  // Md
  i++; vatom_symbol[i]="No"; vatom_name[i]="Nobelium";   vatom_mass[i]*=259.1;    vatom_volume[i]=NNN;     vatom_valence_std[i]=16; vatom_valence_iupac[i]=3; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=NNN; vatom_radius_covalent[i]=NNN;   vatom_electronegativity[i]=NNN; vatom_crystal[i]="nnn";  vatom_pettifor_scale[i]=0; vatom_pearson_coefficient[i]=0.0;  // No
  i++; vatom_symbol[i]="Lr"; vatom_name[i]="Lawrencium"; vatom_mass[i]*=262.11;   vatom_volume[i]=NNN;     vatom_valence_std[i]=17; vatom_valence_iupac[i]=3; vatom_miedema_phi_star[i]=NNN;  vatom_miedema_nws[i]=NNN;  vatom_miedema_Vm[i]=NNN; vatom_miedema_gamma_s[i]= NNN;  vatom_miedema_BVm[i]=NNN;  vatom_radius[i]=NNN; vatom_radius_covalent[i]=NNN;   vatom_electronegativity[i]=NNN; vatom_crystal[i]="nnn";  vatom_pettifor_scale[i]=0; vatom_pearson_coefficient[i]=0.0;  // Lr

  //   int valence_WSETYAWAN[]={9999,1,0,1,2,3,4,-3,-2,-1,0,1,2,3,4,-3,-2,-1,0,1,2,3,4,3,3,2,3,3,2,1,2,3,4,-3,-2,-1,0,1,2,3,4,2,4,7,4,3,2,1,2,3,4,-3,-2,-1,0,1,2,3,3,3,3,3,3,2,3,3,3,3,3,3,3,3,4,5,3,4,4,4,2,1,2,3,4,3,0,0,0,1,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
  //   for(uint j=1;j<93;j++) {
  //   int test=vatom_valence_std[j];
  //   if(j<=9 && test>=5) {test=test-8;} // first row
  //   if(j<=17 && test>=5) {test=test-8;} // second row
  ////   if(test!=valence_WSETYAWAN[j]) 
  //cerr << vatom_symbol[j] << " " << vatom_valence_std[j] << " " <<  test  << " " << valence_WSETYAWAN[j];
  //if(test!=valence_WSETYAWAN[j]) cerr << "   ****   ";
  //cerr << endl;
  //}
  //throw aurostd::xerror(__AFLOW_FILE__,XPID+"atoms_initialize():","Throw for debugging purposes.",_GENERIC_ERROR_);


  //   for(int i=1;i<90;i++) {
  //     cerr << i << " " <<  vatom_symbol[i] << " " << vatom_pettifor_scale[i] << endl;
  //   }
  // cerr << i << " " << vatom_symbol[i] << endl; throw aurostd::xerror(__AFLOW_FILE__,XPID+"atoms_initialize","Throw for debugging purposes.",_GENERIC_ERROR_);

  // aconvasp stuff
  // All data collected from the NIST online tables:
  // http://physics.nist.gov/PhysRefData/FFast/html/form.html
  // All data are ideally for f1 values for Cu-alpha (wavelength=1.5418A, E=8.0416keV).
  // These are for E=7.9026keV (Cu-alpha is wavelength=1.5418A, E=8.0416keV).
  vatom_xray_scatt[1+2]=3.00145E+00; // Li
  vatom_xray_scatt[1+14]=1.53133E+01; // P
  vatom_xray_scatt[1+24]=2.43589E+01; // Mn
  vatom_xray_scatt[1+25]=2.46830E+01; // Fe

  // All data collected from the online tables:
  // http://www-cxro.lbl.gov/optical_constants/pert_form.html
  // All data are f1 values for Cu-alpha (wavelength=1.5418A, E=8.0416keV).
  vatom_xray_scatt[1+0]=1.000; // H
  vatom_xray_scatt[1+1]=2.000; // He
  vatom_xray_scatt[1+2]=3.001; // Li
  vatom_xray_scatt[1+6]=6.019; // C
  vatom_xray_scatt[1+7]=8.052; // O
  vatom_xray_scatt[1+13]=14.43; // P
  vatom_xray_scatt[1+14]=15.30; // P
  vatom_xray_scatt[1+20]=21.34; // Sc
  vatom_xray_scatt[1+21]=22.24; // Ti
  vatom_xray_scatt[1+23]=23.84; // Cr
  vatom_xray_scatt[1+24]=24.46; // Mn
  vatom_xray_scatt[1+25]=24.85; // Fe
  vatom_xray_scatt[1+26]=24.59; // Co
  vatom_xray_scatt[1+27]=25.02; // Ni
  vatom_xray_scatt[1+28]=27.03; // Cu
  vatom_xray_scatt[1+29]=28.44; // Zn
  vatom_xray_scatt[1+46]=47.18; // Ag
  vatom_xray_scatt[1+78]=74.99; // Au
  vatom_xray_scatt[1+89]=86.64; // Th

  // Atomic masses
  // All indices are the atomic number shifted back by one.
  // All masses are in kilograms

  // not useful anymore all the masses are declared
  //  vatom_mass=vector<double> (NUM_ELEMENTS,0.0);
  //  for(i=0;i<NUM_ELEMENTS;i++) {
  //    vatom_mass[i]=(double)(2*i)*AMU2KILOGRAM;
  //  }
  //  vatom_mass[1+0]=1.0079*AMU2KILOGRAM; //H
  //  vatom_mass[1+7]=15.9994*AMU2KILOGRAM; //O
  //  vatom_mass[1+24]=54.93805*AMU2KILOGRAM; //Mn

  // finish and copy


}

// **************************************************************************
// Function GetAtomNumber
// **************************************************************************
uint GetAtomNumber(const string& symbol) {
  for(uint iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==vatom_symbol.at(iat) || symbol==vatom_name.at(iat))
      return (uint) iat;
  cerr << "GetAtomNumber symbol not found symbol=" << symbol << "  vatom_name.size()=" << vatom_name.size() << endl;
  return 0; // no symbol found
}

// **************************************************************************
// Function GetAtomName
// **************************************************************************
std::string GetAtomName(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==vatom_symbol.at(iat) || symbol==vatom_name.at(iat))
      return vatom_name.at(iat);
  return symbol;
}
std::string GetAtomName(const uint& atnum) {
  if(atnum>=vatom_name.size() || atnum<=0) {
    cerr << "GetAtomName out of boundary  atnum=" << atnum << "  vatom_name.size()=" << vatom_name.size() << endl;
    return "not found";
  }
  return vatom_name.at(atnum);
}

// **************************************************************************
// Function GetAtomSymbol
// **************************************************************************
std::string GetAtomSymbol(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==vatom_symbol.at(iat) || symbol==vatom_symbol.at(iat))
      return vatom_symbol.at(iat);
  return symbol;
}
std::string GetAtomSymbol(const uint& atnum) {
  if(atnum>=vatom_symbol.size() || atnum<=0) {
    cerr << "GetAtomSymbol out of boundary  atnum=" << atnum << "  vatom_symbol.size()=" << vatom_symbol.size() << endl;
    return "not found";
  }
  return vatom_symbol.at(atnum);
}

// **************************************************************************
// Function GetAtomMass
// **************************************************************************
double GetAtomMass(const string& _symbol,bool clean) { //CO20181128
  string symbol=_symbol; //CO20181128
  if(clean){symbol=KBIN::VASP_PseudoPotential_CleanName(symbol);} //CO20181128
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==vatom_symbol.at(iat) || symbol==vatom_name.at(iat))
      return vatom_mass.at(iat);
  return (double) NNN;
}

double GetAtomMass(const uint& atnum) {
  if(atnum>=vatom_mass.size() || atnum<=0) {
    cerr << "GetAtomMass out of boundary  atnum=" << atnum << "  vatom_mass.size()=" << vatom_mass.size() << endl;
    return (double) NNN;
  }
  return vatom_mass.at(atnum);
}

// **************************************************************************
// Function GetAtomComptonCrossSection
// **************************************************************************
double GetAtomComptonCrossSection(const string& symbol) {
  //  cout << "symbol=" << symbol << endl;
  //  cout << "GetAtomNumber(symbol)=" << GetAtomNumber(symbol) << endl;
  return (double) GetAtomComptonCrossSection(GetAtomNumber(symbol));
}

double GetAtomComptonCrossSection(const uint& atnum) { //sigma_c at 662 KeV, compton cross section in barn (1 barn = 1e-28 m^2) [ref: Ortiz, Comp Mat Sci 44, 1042 (2009)]
  double sigma_c[]={0.29,0.57,0.86,1.15,1.43,1.72,2.00,2.29,2.58,2.86,3.15,3.43,3.72,4.00,4.29,4.57,4.86,5.14,5.43,5.71,5.99,6.28,6.56,6.84,7.13,7.41,7.69,7.98,8.26,8.54,8.82,9.11,9.39,9.67,9.95,10.24,10.52,10.79,11.08,11.36,11.64,11.92,12.20,12.48,12.76,13.04,13.32,13.60,13.88,14.16,14.44,14.72,15.00,15.27,15.55,15.83,16.11,16.38,16.66,16.94,17.22,17.50,17.77,18.05,18.33,18.60,18.88,19.16,19.44,19.71,19.98,20.26,20.54,20.81,21.08,21.36,21.64,21.91,22.19,22.46,22.73,23.01,23.28,23.55,23.82,24.10,24.38,24.64,24.92,-1.000};
  vector<double> vsigma_c;vsigma_c.push_back(NNN);
  for(uint i=0;sigma_c[i]>0;i++) vsigma_c.push_back(sigma_c[i]);
  if(atnum>=vsigma_c.size() || atnum<=0) {
    cerr << "GetAtomComptonCrossSection out of boundary  atnum=" << atnum << "  vsigma_c.size()=" << vsigma_c.size() << endl;
    return (double) NNN;
  }
  return vsigma_c.at(atnum);
}

// **************************************************************************
// Function GetAtomPhotoelectricCrossSection
// **************************************************************************
double GetAtomPhotoelectricCrossSection(const string& symbol) {
  //  cout << "symbol=" << symbol << endl;
  //  cout << "GetAtomNumber(symbol)=" << GetAtomNumber(symbol) << endl;
  return (double) GetAtomPhotoelectricCrossSection(GetAtomNumber(symbol));
}

double GetAtomPhotoelectricCrossSection(const uint& atnum) { //sigma_pe photoelectric cross section in barn (1 barn = 1e-28 m^2) %ref: Ortiz, Comp Mat Sci 44, 1042 (2009)
  double sigma_pe[]={8.79e-09,2.70e-07,2.62e-06,1.33e-05,4.59e-05,1.22e-04,2.75e-04,5.51e-04,9.85e-04,1.69e-03,2.62e-03,3.97e-03,6.10e-03,8.46e-03,1.18e-02,1.69e-02,2.19e-02,2.87e-02,3.68e-02,4.91e-02,5.94e-02,7.74e-02,9.47e-02,1.16e-01,1.40e-01,1.59e-01,1.94e-01,2.21e-01,2.65e-01,3.10e-01,3.56e-01,4.14e-01,4.68e-01,5.49e-01,6.24e-01,7.09e-01,7.98e-01,9.16e-01,1.02,1.10,1.27,1.49,1.57,1.74,1.92,2.15,2.36,2.60,2.87,3.06,3.30,3.66,4.06,4.31,4.68,5.03,5.47,5.86,6.41,6.95,7.43,7.96,8.56,9.02,9.78,1.05e+01,1.11e+01,1.19e+01,1.27e+01,1.37e+01,1.44e+01,1.48e+01,1.57e+01,1.65e+01,1.78e+01,1.94e+01,2.06e+01,2.11e+01,2.23e+01,2.36e+01,2.56e+01,2.65e+01,2.80e+01,3.00e+01,3.15e+01,3.33e+01,3.47e+01,3.65e+01,3.82e+01,-1.000};
  vector<double> vsigma_pe;vsigma_pe.push_back(NNN);
  for(uint i=0;sigma_pe[i]>0;i++) vsigma_pe.push_back(sigma_pe[i]);
  if(atnum>=vsigma_pe.size() || atnum<=0) {
    cerr << "GetAtomPhotoelectricCrossSection out of boundary  atnum=" << atnum << "  vsigma_pe.size()=" << vsigma_pe.size() << endl;
    return (double) NNN;
  }
  return vsigma_pe.at(atnum);
}

// **************************************************************************
// Function GetAtomVolume
// **************************************************************************
double GetAtomVolume(const string& _symbol,bool clean) { //CO20181128
  string symbol=_symbol; //CO20181128
  if(clean){symbol=KBIN::VASP_PseudoPotential_CleanName(symbol);} //CO20181128
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==vatom_symbol.at(iat) || symbol==vatom_name.at(iat))
      return vatom_volume.at(iat);
  return (double) NNN;
}
double GetAtomVolume(const uint& atnum) {
  if(atnum>=vatom_volume.size() || atnum<=0) {
    cerr << "GetAtomVolume out of boundary  atnum=" << atnum << "  vatom_volume.size()=" << vatom_volume.size() << endl;
    return (double) NNN;
  }
  return vatom_volume.at(atnum);
}

// **************************************************************************
// Function GetAtomValenceIupac
// **************************************************************************
int GetAtomValenceIupac(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==vatom_symbol.at(iat) || symbol==vatom_name.at(iat))
      return vatom_valence_iupac.at(iat);
  return (int) NNN;
}

int GetAtomValenceIupac(const uint& atnum) {
  if(atnum>=vatom_valence_iupac.size() || atnum<=0) {
    cerr << "GetAtomValenceIupac out of boundary  atnum=" << atnum << "  vatom_valence_iupac.size()=" << vatom_valence_iupac.size() << endl;
    return (int) NNN;
  }
  return vatom_valence_iupac.at(atnum);
}

// **************************************************************************
// Function GetAtomValenceStd
// **************************************************************************
int GetAtomValenceStd(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==vatom_symbol.at(iat) || symbol==vatom_name.at(iat))
      return vatom_valence_std.at(iat);
  return (int) NNN;
}

int GetAtomValenceStd(const uint& atnum) {
  if(atnum>=vatom_valence_std.size() || atnum<=0) {
    cerr << "GetAtomValenceStd out of boundary  atnum=" << atnum << "  vatom_valence_std.size()=" << vatom_valence_std.size() << endl;
    return (int) NNN;
  }
  return vatom_valence_std.at(atnum);
}

// **************************************************************************
// Function GetAtomRadius
// **************************************************************************
double GetAtomRadius(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==vatom_symbol.at(iat) || symbol==vatom_name.at(iat))
      return vatom_radius.at(iat);
  return (double) NNN;
}
double GetAtomRadius(const uint& atnum) {
  if(atnum>=vatom_radius.size() || atnum<=0) {
    cerr << "GetAtomRadius out of boundary  atnum=" << atnum << "  vatom_radius.size()=" << vatom_radius.size() << endl;
    return (double) NNN;
  }
  return vatom_radius.at(atnum);
}

//DX+CO20170904 START 
// **************************************************************************
// Function GetAtomRadiusCovalent
// **************************************************************************
double GetAtomRadiusCovalent(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==vatom_symbol.at(iat) || symbol==vatom_name.at(iat))
      return vatom_radius_covalent.at(iat);
  return (double) NNN;
}
double GetAtomRadiusCovalent(const uint& atnum) {
  if(atnum>=vatom_radius_covalent.size() || atnum<=0) {
    cerr << "GetAtomRadiusCovalent out of boundary  atnum=" << atnum << "  vatom_radius_covalent.size()=" << vatom_radius_covalent.size() << endl;
    return (double) NNN;
  }
  return vatom_radius_covalent.at(atnum);
}
//DX+CO20170904 END 

// **************************************************************************
// Function GetAtomElectronegativity
// **************************************************************************
double GetAtomElectronegativity(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==vatom_symbol.at(iat) || symbol==vatom_name.at(iat))
      return vatom_electronegativity.at(iat);
  return (double) NNN;
}
double GetAtomElectronegativity(const uint& atnum) {
  if(atnum>=vatom_electronegativity.size() || atnum<=0) {
    cerr << "GetAtomElectronegativity out of boundary  atnum=" << atnum << "  vatom_electronegativity.size()=" << vatom_electronegativity.size() << endl;
    return (double) NNN;
  }
  return vatom_electronegativity.at(atnum);
}

// **************************************************************************
// Function GetAtomCrystal
// **************************************************************************
string GetAtomCrystal(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==vatom_symbol.at(iat) || symbol==vatom_name.at(iat))
      return vatom_crystal.at(iat);
  return (string) "nnn";
}
string GetAtomCrystal(const uint& atnum) {
  if(atnum>=vatom_crystal.size() || atnum<=0) {
    cerr << "GetAtomCrystal out of boundary  atnum=" << atnum << "  vatom_crystal.size()=" << vatom_crystal.size() << endl;
    return (string) "nnn";
  }
  return vatom_crystal.at(atnum);
}

// **************************************************************************
// Function GetAtomPettiforScale
// **************************************************************************
double GetAtomPettiforScale(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==vatom_symbol.at(iat) || symbol==vatom_name.at(iat))
      return vatom_pettifor_scale.at(iat);
  return 0.0;
}

double GetAtomPettiforScale(const uint& atnum) {
  if(atnum>=vatom_crystal.size() || atnum<=0) {
    cerr << "GetAtomPettiforScale out of boundary  atnum=" << atnum << "  vatom_crystal.size()=" << vatom_crystal.size() << endl;
    return 0.0;
  }
  return vatom_pettifor_scale.at(atnum);
}

bool GetAtomPettiforScale(const vector<string>& vsymbol,vector<double>& vvalue) {
  vvalue.clear(); // delete
  for(uint i=0;i<vsymbol.size();i++)
    vvalue.push_back(GetAtomPettiforScale(vsymbol.at(i)));
  return TRUE;
}

bool GetAtomPettiforScale(const vector<uint>& vatnum,vector<double>& vvalue) {
  vvalue.clear(); // delete
  for(uint i=0;i<vatnum.size();i++)
    vvalue.push_back(GetAtomPettiforScale(vatnum.at(i)));
  return TRUE;
}

bool GetAtomPettiforScale(const vector<string>& vsymbol,xvector<double>& vvalue) {
  if(vvalue.rows!=(int) vsymbol.size()) return FALSE; // nothing to be ordered
  if(vvalue.lrows!=1) return FALSE; // start from 1
  for(uint i=0;i<vsymbol.size();i++)
    vvalue[i+1]=GetAtomPettiforScale(vsymbol.at(i));
  return TRUE;
}

bool GetAtomPettiforScale(const vector<uint>& vatnum,xvector<double>& vvalue) {
  if(vvalue.rows!=(int) vatnum.size()) return FALSE; // nothing to be ordered
  if(vvalue.lrows!=1) return FALSE; // start from 1
  for(uint i=0;i<vatnum.size();i++)
    vvalue[i+1]=GetAtomPettiforScale(vatnum.at(i));
  return TRUE;
}

bool SortAtomsPettiforScale(vector<string> &vsymbol,xvector<int> &vorder,xvector<double> &vvalue) {
  if(vorder.rows!=(int) vsymbol.size()) return FALSE; // nothing to be ordered
  if(vvalue.rows!=(int) vsymbol.size()) return FALSE; // nothing to be ordered
  if(vorder.lrows!=1) return FALSE; // start from 1 .. and contains order from 1
  if(vvalue.lrows!=1) return FALSE; // start from 1
  // build
  for(uint i=1;i<=vsymbol.size();i++) {
    vvalue[i]=GetAtomPettiforScale(vsymbol.at(i-1));
    vorder[i]=i;
  }
  aurostd::sort2(vsymbol.size(),vvalue,vorder);
  vector<string> vsymbol_tmp(vsymbol);
  vsymbol.clear();
  for(uint i=1;i<=vsymbol_tmp.size();i++)
    vsymbol.push_back(vsymbol_tmp.at(vorder[i]-1));
  return TRUE;
}

bool SortAtomsPettiforScale(vector<string> &vsymbol,vector<int> &vorder,vector<double> &vvalue) {
  //   vorder.clear();vvalue.clear();
  //   xvector<int> xvorder(vsymbol.size());
  //   xvector<double> xvvalue(vsymbol.size());
  //   SortAtomsPettiforScale(vsymbol,xvorder,xvvalue);
  //   for(uint i=0;i<vsymbol.size();i++) {
  //     vvalue.push_back(xvvalue(i+1));
  //     vorder.push_back(xvorder(i+1));
  //   }
  vector<string> vsymbol_tmp(vsymbol);
  vorder.clear();vvalue.clear();
  for(uint i=0;i<vsymbol.size();i++) {
    vvalue.push_back(GetAtomPettiforScale(vsymbol.at(i)));
    vorder.push_back(i);
  }
  aurostd::sort(vvalue,vorder);
  vsymbol.clear();
  for(uint i=0;i<vsymbol_tmp.size();i++)
    vsymbol.push_back(vsymbol_tmp.at(vorder.at(i)));
  return TRUE;
}

bool SortAtomsPettiforScale(vector<string> &vsymbol,vector<int> &vorder) {
  vector<double> vvalue;
  SortAtomsPettiforScale(vsymbol,vorder,vvalue);
  return TRUE;
}

bool SortAtomsPettiforScale(vector<string> &vsymbol,vector<double> &vvalue) {
  vector<int> vorder;
  SortAtomsPettiforScale(vsymbol,vorder,vvalue);
  return TRUE;
}

bool SortAtomsPettiforScale(vector<string> &vsymbol) {
  vector<double> vvalue;
  vector<int> vorder;
  SortAtomsPettiforScale(vsymbol,vorder,vvalue);
  return TRUE;
}

// **************************************************************************
// Function GetAtomXrayScatt
// **************************************************************************
double GetAtomXrayScatt(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==vatom_symbol.at(iat) || symbol==vatom_name.at(iat))
      return vatom_xray_scatt.at(iat);
  return (double) NNN;
}
double GetAtomXrayScatt(const uint& atnum) {
  if(atnum>=vatom_xray_scatt.size() || atnum<=0) {
    cerr << "GetAtomXrayScatt out of boundary  atnum=" << atnum << "  vatom_xray_scatt.size()=" << vatom_xray_scatt.size() << endl;
    return (double) NNN;
  }
  return vatom_xray_scatt.at(atnum);
}

//DX20181220 - get group of atoms - START 
// **************************************************************************
// Function GetGroupOfAtoms
// **************************************************************************
vector<string> GetGroupOfAtoms(string& group_name){
  vector<string> element_list;
  if(aurostd::tolower(group_name) == "metals"){
    element_list.push_back("Li"); element_list.push_back("Be"); element_list.push_back("Na"); element_list.push_back("Mg");
    element_list.push_back("Al"); element_list.push_back("K"); element_list.push_back("Ca"); element_list.push_back("Sc"); 
    element_list.push_back("Ti"); element_list.push_back("V"); element_list.push_back("Cr"); element_list.push_back("Mn"); 
    element_list.push_back("Fe"); element_list.push_back("Co"); element_list.push_back("Ni"); element_list.push_back("Cu"); 
    element_list.push_back("Zn"); element_list.push_back("Ga"); element_list.push_back("Rb"); element_list.push_back("Sr"); 
    element_list.push_back("Y"); element_list.push_back("Zr"); element_list.push_back("Nb"); element_list.push_back("Mo"); 
    element_list.push_back("Tc"); element_list.push_back("Ru"); element_list.push_back("Rh"); element_list.push_back("Pd"); 
    element_list.push_back("Ag"); element_list.push_back("Cd"); element_list.push_back("In"); element_list.push_back("Sn"); 
    element_list.push_back("Cs"); element_list.push_back("Ba"); element_list.push_back("La"); element_list.push_back("Ce"); 
    element_list.push_back("Pr"); element_list.push_back("Nd"); element_list.push_back("Pm"); element_list.push_back("Sm"); 
    element_list.push_back("Eu"); element_list.push_back("Gd"); element_list.push_back("Tb"); element_list.push_back("Dy"); 
    element_list.push_back("Ho"); element_list.push_back("Er"); element_list.push_back("Tm"); element_list.push_back("Yb"); 
    element_list.push_back("Lu"); element_list.push_back("Hf"); element_list.push_back("Ta"); element_list.push_back("W"); 
    element_list.push_back("Re"); element_list.push_back("Os"); element_list.push_back("Ir"); element_list.push_back("Pt"); 
    element_list.push_back("Au"); element_list.push_back("Hg"); element_list.push_back("Tl"); element_list.push_back("Pb"); 
    element_list.push_back("Bi"); element_list.push_back("Po");
  }
  if(aurostd::tolower(group_name) == "alkali"){
    element_list.push_back("Li"); element_list.push_back("Na"); element_list.push_back("K"); element_list.push_back("Rb");
    element_list.push_back("Cs");
  }
  if(aurostd::tolower(group_name) == "alkaline"){
    element_list.push_back("Be"); element_list.push_back("Mg"); element_list.push_back("Ca"); element_list.push_back("Sr");
    element_list.push_back("Ba");
  }
  if(aurostd::tolower(group_name) == "transition_metals" || aurostd::tolower(group_name) == "transition-metals"){
    element_list.push_back("Sc"); element_list.push_back("Ti");
    element_list.push_back("V"); element_list.push_back("Cr"); element_list.push_back("Mn"); element_list.push_back("Fe");
    element_list.push_back("Co"); element_list.push_back("Ni"); element_list.push_back("Cu"); element_list.push_back("Zn");
    element_list.push_back("Y");
    element_list.push_back("Zr"); element_list.push_back("Nb"); element_list.push_back("Mo"); element_list.push_back("Tc");
    element_list.push_back("Ru"); element_list.push_back("Rh"); element_list.push_back("Pd"); element_list.push_back("Ag");
    element_list.push_back("Cd");
    element_list.push_back("Hf");
    element_list.push_back("Ta"); element_list.push_back("W"); element_list.push_back("Re"); element_list.push_back("Os");
    element_list.push_back("Ir"); element_list.push_back("Pt"); element_list.push_back("Au"); element_list.push_back("Hg");
    element_list.push_back("Tl");
  }
  if(aurostd::tolower(group_name) == "lanthanides"){
    element_list.push_back("La"); element_list.push_back("Ce"); element_list.push_back("Pr"); element_list.push_back("Nd");
    element_list.push_back("Pm"); element_list.push_back("Sm"); element_list.push_back("Eu"); element_list.push_back("Gd");
    element_list.push_back("Tb"); element_list.push_back("Dy"); element_list.push_back("Ho"); element_list.push_back("Er");
    element_list.push_back("Tm"); element_list.push_back("Yb"); element_list.push_back("Lu");
  }
  if(aurostd::tolower(group_name) == "nonmetals" || aurostd::tolower(group_name) == "non-metals" || aurostd::tolower(group_name) == "non_metals"){
    cerr << "in non metals" << endl;
    element_list.push_back("He"); element_list.push_back("B"); element_list.push_back("C"); element_list.push_back("N");
    element_list.push_back("O"); element_list.push_back("F"); element_list.push_back("Ne"); element_list.push_back("Si");
    element_list.push_back("P"); element_list.push_back("S"); element_list.push_back("Cl"); element_list.push_back("Ar");
    element_list.push_back("As"); element_list.push_back("Se"); element_list.push_back("Br"); element_list.push_back("Kr");
    element_list.push_back("Te"); element_list.push_back("I"); element_list.push_back("Xe");
  }
  print(element_list);
  return element_list;
}
//DX20181220 - get group of atoms - END 

// **************************************************************************
// Function GetPearsonCoefficient
// **************************************************************************
double GetPearsonCoefficient(const string& symbol) {
  for (int iat = 0; iat < NUM_ELEMENTS; iat++) {
    if (symbol == vatom_symbol[iat] || symbol == vatom_name[iat]) {
      return GetPearsonCoefficient(iat);
    }
  }
  // If not found throw xerror
  string message = symbol + " is not a valid element name or symbol.";
  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
}

double GetPearsonCoefficient(const int& iat) {
  return vatom_pearson_coefficient.at(iat);
}

// **************************************************************************
// Function GetCompoundAttenuationLength
// **************************************************************************
double GetCompoundAttenuationLength(const vector<string>& species,const vector<double>& composition,const double& density) { // density in g/cm^3, return in cm
  if(species.size()!=composition.size()) {
    stringstream message; //CO20190629
    message << "species.size()[" << species.size() << "]!=composition.size()[" << composition.size() << "]";  //CO20190629
    throw aurostd::xerror(__AFLOW_FILE__,"GetCompoundAttenuationLength():",message,_INDEX_MISMATCH_);  //CO20190629
  }
  // cout << "Density=" << density << "<br>" << endl;
  double numerator=0.0,denominator=0.0;
  for(uint i=0;i<species.size();i++) {
    numerator+=composition.at(i)*GetAtomMass(species.at(i))*1000;// from Kg to grams  
    denominator+=density*(composition.at(i)*(GetAtomComptonCrossSection(species.at(i))+GetAtomPhotoelectricCrossSection(species.at(i)))*1e-24);
  }
  //  cout << "numerator=" << numerator << endl;
  //  cout << "denominator=" << denominator << endl;
  return (double) numerator/denominator; // in cm
}

double GetCompoundAttenuationLength(const deque<string>& _species,const deque<int>& _composition,const double& density) { // density in g/cm^3, return in cm
  vector<double> composition;
  for(uint i=0;i<_composition.size();i++) 
    composition.push_back(double(_composition.at(i)));
  vector<string> species;
  for(uint i=0;i<_species.size();i++) 
    species.push_back(string(_species.at(i)));
  return GetCompoundAttenuationLength(species,composition,density);
}

// **************************************************************************
// Function XATOM_AlphabetizationSpecies & XATOM_AlphabetizationCompound
// **************************************************************************
string XATOM_AlphabetizationSpecies(const string& speciesA,const string& speciesB) {
  string system="";
  if(speciesA<=speciesB) system=speciesA+speciesB; else system=speciesB+speciesA;
  return system;
}

string XATOM_AlphabetizationSpecies(const vector<string>& vspecies_in) {
  string system="";
  vector<string> vspecies(vspecies_in);std::sort(vspecies.begin(),vspecies.end());
  return aurostd::joinWDelimiter(vspecies,"");
  //[CO20200624 - OBSOLETE]aurostd::sort(vspecies);
  //[CO20200624 - OBSOLETE]for(uint i=0;i<vspecies.size();i++)
  //[CO20200624 - OBSOLETE]  system << vspecies.at(i);
  //[CO20200624 - OBSOLETE]return system.str();
}

string XATOM_AlphabetizationSpecies(const vector<string>& vspecies_in,const vector<double>& vnumbers_in) {
  stringstream system;
  vector<string> vspecies(vspecies_in);
  vector<double> vnumbers(vnumbers_in);
  aurostd::sort(vspecies,vnumbers);
  for(uint i=0;i<vspecies.size();i++)
    system << vspecies[i] << vnumbers[i];
  return system.str();
}

void XATOM_AlphabetizationSpecies(string& system, vector<string>& vspecies,vector<double>& vnumbers) {
  vspecies.clear();vnumbers.clear();
  KBIN::VASP_SplitAlloySpecies(KBIN::VASP_PseudoPotential_CleanName(system),vspecies,vnumbers);
  stringstream systemstream;
  aurostd::sort(vspecies,vnumbers);
  for(uint i=0;i<vspecies.size();i++)
    systemstream << vspecies.at(i);
  system=systemstream.str();
}

void XATOM_AlphabetizationCompound(string& system, vector<string>& vspecies,vector<double>& vnumbers) {
  vspecies.clear();vnumbers.clear();
  KBIN::VASP_SplitAlloySpecies(KBIN::VASP_PseudoPotential_CleanName(system),vspecies,vnumbers);
  stringstream systemstream;
  aurostd::sort(vspecies,vnumbers);
  for(uint i=0;i<vspecies.size();i++)
    systemstream << vspecies.at(i) << vnumbers.at(i);
  system=systemstream.str();
}

void XATOM_AlphabetizationSpecies(string& system, vector<string>& vspecies) {
  vector<double> vnumbers(vspecies.size());
  XATOM_AlphabetizationSpecies(system,vspecies,vnumbers);
}

void XATOM_AlphabetizationSpecies(string& system) {
  vector<string> vspecies;
  vector<double> vnumbers;
  XATOM_AlphabetizationSpecies(system,vspecies,vnumbers);
}

void XATOM_AlphabetizationCompound(string& system) {
  vector<string> vspecies;
  vector<double> vnumbers;
  XATOM_AlphabetizationCompound(system,vspecies,vnumbers);
}

// **************************************************************************
// Function XATOM_SplitAlloySpecies
// **************************************************************************
uint XATOM_SplitAlloySpecies(const string& alloy_in, vector<string> &speciesX) {
  speciesX=aurostd::getElements(alloy_in); //CO20200624
  //[CO20200624 - OBSOLETE]string alloy=alloy_in;
  //[CO20200624 - OBSOLETE]alloy=KBIN::VASP_PseudoPotential_CleanName(alloy); // always some cleaning is good
  //[CO20200624 - OBSOLETE]alloy=KBIN::VASP_PseudoPotential_CleanName(alloy); // always some cleaning is good
  //[CO20200624 - OBSOLETE]alloy=aurostd::RemoveNumbers(alloy);              // remove composition
  //[CO20200624 - OBSOLETE]speciesX.clear();
  //[CO20200624 - OBSOLETE]for(uint i=0;i<alloy.length();i++) {
  //[CO20200624 - OBSOLETE]  if(alloy[i]>='A' && alloy[i]<='Z') speciesX.push_back("");
  //[CO20200624 - OBSOLETE]  speciesX.at(speciesX.size()-1)+=alloy[i];
  //[CO20200624 - OBSOLETE]}
  //[CO20200624 - OBSOLETE]for(uint i=0;i<speciesX.size();i++)
  //[CO20200624 - OBSOLETE]  speciesX.at(i)=aurostd::CleanStringASCII(speciesX.at(i));
  return speciesX.size();
}

//[CO20200624 - OBSOLETE]uint new_XATOM_SplitAlloySpecies(string alloy_in, vector<string> &speciesX, vector<double> &natomsX) {
//[CO20200624 - OBSOLETE]  string alloy=alloy_in,alloyn;
//[CO20200624 - OBSOLETE]  string letters="QWERTYUIOPASDFGHJKLZXCVBNMqwertyuiopasdfghjklzxcvbnm";
//[CO20200624 - OBSOLETE]  string numbers="0123456789";
//[CO20200624 - OBSOLETE]  cerr << alloy << endl;// throw aurostd::xerror(__AFLOW_FILE__,XPID+"atoms_initialize","Throw for debugging purposes.",_GENERIC_ERROR_);
//[CO20200624 - OBSOLETE]
//[CO20200624 - OBSOLETE]  alloy=KBIN::VASP_PseudoPotential_CleanName(alloy); // always some cleaning is good
//[CO20200624 - OBSOLETE]  cerr << alloy << endl;// throw aurostd::xerror(__AFLOW_FILE__,XPID+"atoms_initialize","Throw for debugging purposes.",_GENERIC_ERROR_);
//[CO20200624 - OBSOLETE]  alloy=KBIN::VASP_PseudoPotential_CleanName(alloy); // always some cleaning is good
//[CO20200624 - OBSOLETE]  cerr << alloy << endl;// throw aurostd::xerror(__AFLOW_FILE__,XPID+"atoms_initialize","Throw for debugging purposes.",_GENERIC_ERROR_);
//[CO20200624 - OBSOLETE]  alloy=aurostd::CleanStringASCII(alloy);
//[CO20200624 - OBSOLETE]  cerr << alloy << endl;// throw aurostd::xerror(__AFLOW_FILE__,XPID+"atoms_initialize","Throw for debugging purposes.",_GENERIC_ERROR_);
//[CO20200624 - OBSOLETE]  for(uint i=0;i<alloy.length();i++)
//[CO20200624 - OBSOLETE]    for(uint j=0;j<letters.length();j++)
//[CO20200624 - OBSOLETE]      if(alloy[i]==letters[j] && alloy[i]!=0) {cerr << alloy[i] << endl; alloy[i]='_';}
//[CO20200624 - OBSOLETE]
//[CO20200624 - OBSOLETE]  cerr << alloy << endl; throw aurostd::xerror(__AFLOW_FILE__,XPID+"atoms_initialize","Throw for debugging purposes.",_GENERIC_ERROR_);
//[CO20200624 - OBSOLETE]
//[CO20200624 - OBSOLETE]  speciesX.clear();
//[CO20200624 - OBSOLETE]  for(uint i=0;i<alloy.length();i++) {
//[CO20200624 - OBSOLETE]    if(alloy[i]>='A' && alloy[i]<='Z') speciesX.push_back("");
//[CO20200624 - OBSOLETE]    speciesX.at(speciesX.size()-1)+=alloy[i];
//[CO20200624 - OBSOLETE]  }
//[CO20200624 - OBSOLETE]  for(uint i=0;i<speciesX.size();i++)
//[CO20200624 - OBSOLETE]    speciesX.at(i)=aurostd::CleanStringASCII(speciesX.at(i));   // clean it up so it does not have problems inside only letters_numbers
//[CO20200624 - OBSOLETE]  // now the atoms
//[CO20200624 - OBSOLETE]  natomsX.clear();
//[CO20200624 - OBSOLETE]  for(uint i=0;i<speciesX.size();i++) {
//[CO20200624 - OBSOLETE]    if(i<speciesX.size()-1)
//[CO20200624 - OBSOLETE]      natomsX.push_back(aurostd::string2utype<double>(alloyn.substr(alloyn.find(speciesX.at(i))+speciesX.at(i).length(),alloyn.find(speciesX.at(i+1))-alloyn.find(speciesX.at(i))-speciesX.at(i).length())));
//[CO20200624 - OBSOLETE]    else
//[CO20200624 - OBSOLETE]      natomsX.push_back(aurostd::string2utype<double>(alloyn.substr(alloyn.find(speciesX.at(i))+speciesX.at(i).length())));
//[CO20200624 - OBSOLETE]    if(abs(natomsX.at(natomsX.size()-1))<=0.00001) natomsX.at(natomsX.size()-1)=1.0;  // fix the no number = 1
//[CO20200624 - OBSOLETE]  }
//[CO20200624 - OBSOLETE]  //  for(uint i=0;i<natomsX.size();i++)
//[CO20200624 - OBSOLETE]  //  cerr << natomsX.at(i) << endl;
//[CO20200624 - OBSOLETE]  return speciesX.size();
//[CO20200624 - OBSOLETE]}

uint XATOM_SplitAlloySpecies(const string& alloy_in, vector<string> &speciesX, vector<double> &natomsX) {
  speciesX=aurostd::getElements(alloy_in,natomsX);
  //[CO20200624 - OBSOLETE]string alloy=alloy_in,alloyn;
  //[CO20200624 - OBSOLETE]alloy=KBIN::VASP_PseudoPotential_CleanName(alloy); // always some cleaning is good
  //[CO20200624 - OBSOLETE]alloy=KBIN::VASP_PseudoPotential_CleanName(alloy); // always some cleaning is good
  //[CO20200624 - OBSOLETE]alloyn=aurostd::CleanStringASCII(alloy);
  //[CO20200624 - OBSOLETE]alloy=aurostd::RemoveNumbers(alloy);              // remove composition
  //[CO20200624 - OBSOLETE]speciesX.clear();
  //[CO20200624 - OBSOLETE]for(uint i=0;i<alloy.length();i++) {
  //[CO20200624 - OBSOLETE]  if(alloy[i]>='A' && alloy[i]<='Z') speciesX.push_back("");
  //[CO20200624 - OBSOLETE]  speciesX.at(speciesX.size()-1)+=alloy[i];
  //[CO20200624 - OBSOLETE]}
  //[CO20200624 - OBSOLETE]for(uint i=0;i<speciesX.size();i++)
  //[CO20200624 - OBSOLETE]  speciesX.at(i)=aurostd::CleanStringASCII(speciesX.at(i));   // clean it up so it does not have problems inside only letters_numbers
  //[CO20200624 - OBSOLETE]// now the atoms
  //[CO20200624 - OBSOLETE]natomsX.clear();
  //[CO20200624 - OBSOLETE]for(uint i=0;i<speciesX.size();i++) {
  //[CO20200624 - OBSOLETE]  if(i<speciesX.size()-1)
  //[CO20200624 - OBSOLETE]    natomsX.push_back(aurostd::string2utype<double>(alloyn.substr(alloyn.find(speciesX.at(i))+speciesX.at(i).length(),alloyn.find(speciesX.at(i+1))-alloyn.find(speciesX.at(i))-speciesX.at(i).length())));
  //[CO20200624 - OBSOLETE]  else
  //[CO20200624 - OBSOLETE]    natomsX.push_back(aurostd::string2utype<double>(alloyn.substr(alloyn.find(speciesX.at(i))+speciesX.at(i).length())));
  //[CO20200624 - OBSOLETE]  if(abs(natomsX.at(natomsX.size()-1))<=0.00001) natomsX.at(natomsX.size()-1)=1.0;  // fix the no number = 1
  //[CO20200624 - OBSOLETE]}
  //[CO20200624 - OBSOLETE]//  for(uint i=0;i<natomsX.size();i++)
  //[CO20200624 - OBSOLETE]//  cerr << natomsX.at(i) << endl;
  return speciesX.size();
}


// **************************************************************************
// Function XATOM_SplitAlloyPseudoPotentials
// **************************************************************************
uint XATOM_SplitAlloyPseudoPotentials(const string& alloy_in, vector<string> &species_ppX) {
  species_ppX=aurostd::getElements(alloy_in,pp_string,false,false,true); //CO20200624 - no clean or sort, but do keep_pp
  //[CO20200624 - OBSOLETE]string alloy=alloy_in;
  //[CO20200624 - OBSOLETE]alloy=aurostd::RemoveNumbers(alloy);              // remove composition
  //[CO20200624 - OBSOLETE]species_ppX.clear();
  //[CO20200624 - OBSOLETE]for(uint i=0;i<alloy.length();i++) {
  //[CO20200624 - OBSOLETE]  if(alloy[i]>='A' && alloy[i]<='Z') species_ppX.push_back("");
  //[CO20200624 - OBSOLETE]  species_ppX.at(species_ppX.size()-1)+=alloy[i];
  //[CO20200624 - OBSOLETE]}
  //[CO20200624 - OBSOLETE]for(uint i=0;i<species_ppX.size();i++)
  //[CO20200624 - OBSOLETE]  species_ppX.at(i)=aurostd::CleanStringASCII(species_ppX.at(i));
  return species_ppX.size();
}

uint XATOM_SplitAlloyPseudoPotentials(const string& alloy_in, vector<string> &species_ppX, vector<double> &natomsX) {
  species_ppX=aurostd::getElements(alloy_in,natomsX,pp_string,false,false,true); //CO20200624 - no clean or sort, but do keep_pp - will return natomX to be all 1's
  //[CO20200624 - OBSOLETE]string alloy=alloy_in,alloyn;
  //[CO20200624 - OBSOLETE]alloyn=aurostd::CleanStringASCII(alloy);
  //[CO20200624 - OBSOLETE]alloy=aurostd::RemoveNumbers(alloy);              // remove composition
  //[CO20200624 - OBSOLETE]species_ppX.clear();
  //[CO20200624 - OBSOLETE]for(uint i=0;i<alloy.length();i++) {
  //[CO20200624 - OBSOLETE]  if(alloy[i]>='A' && alloy[i]<='Z') species_ppX.push_back("");
  //[CO20200624 - OBSOLETE]  species_ppX.at(species_ppX.size()-1)+=alloy[i];
  //[CO20200624 - OBSOLETE]}
  //[CO20200624 - OBSOLETE]for(uint i=0;i<species_ppX.size();i++)
  //[CO20200624 - OBSOLETE]  species_ppX.at(i)=aurostd::CleanStringASCII(species_ppX.at(i));   // clean it up so it does not have problems inside only letters_numbers
  //[CO20200624 - OBSOLETE]// now the atoms
  //[CO20200624 - OBSOLETE]natomsX.clear();
  //[CO20200624 - OBSOLETE]for(uint i=0;i<species_ppX.size();i++) {
  //[CO20200624 - OBSOLETE]  if(i<species_ppX.size()-1)
  //[CO20200624 - OBSOLETE]    natomsX.push_back(aurostd::string2utype<double>(alloyn.substr(alloyn.find(species_ppX.at(i))+species_ppX.at(i).length(),alloyn.find(species_ppX.at(i+1))-alloyn.find(species_ppX.at(i))-species_ppX.at(i).length())));
  //[CO20200624 - OBSOLETE]  else
  //[CO20200624 - OBSOLETE]    natomsX.push_back(aurostd::string2utype<double>(alloyn.substr(alloyn.find(species_ppX.at(i))+species_ppX.at(i).length())));
  //[CO20200624 - OBSOLETE]  if(abs(natomsX.at(natomsX.size()-1))<=0.00001) natomsX.at(natomsX.size()-1)=1.0;  // fix the no number = 1
  //[CO20200624 - OBSOLETE]}
  //[CO20200624 - OBSOLETE]//  for(uint i=0;i<natomsX.size();i++)
  //[CO20200624 - OBSOLETE]//  cerr << natomsX.at(i) << endl;
  return species_ppX.size();
}

//DX20200724 [OBSOLETE] //DX composition2stoichiometry - 20181009 - START
//DX20200724 [OBSOLETE] vector<uint> composition2stoichiometry(string& composition){
//DX20200724 [OBSOLETE]   vector<uint> stoichiometry;
//DX20200724 [OBSOLETE]   bool is_previous_alpha = false;
//DX20200724 [OBSOLETE]   bool is_previous_digit = false;
//DX20200724 [OBSOLETE]   string number_string = "";
//DX20200724 [OBSOLETE]   for(uint i=0;i<composition.size();i++){
//DX20200724 [OBSOLETE]     if(isalpha(composition[i])){
//DX20200724 [OBSOLETE]       if(is_previous_alpha){
//DX20200724 [OBSOLETE]         stoichiometry.push_back(1);
//DX20200724 [OBSOLETE]       }
//DX20200724 [OBSOLETE]       else if(is_previous_digit){
//DX20200724 [OBSOLETE]         stoichiometry.push_back(aurostd::string2utype<uint>(number_string));
//DX20200724 [OBSOLETE]       }
//DX20200724 [OBSOLETE]     }
//DX20200724 [OBSOLETE]     else if(isdigit(composition[i])){
//DX20200724 [OBSOLETE]       if(is_previous_digit){
//DX20200724 [OBSOLETE]         stringstream tmp; tmp << number_string << composition[i];
//DX20200724 [OBSOLETE]         number_string = tmp.str();
//DX20200724 [OBSOLETE]       }
//DX20200724 [OBSOLETE]       else {
//DX20200724 [OBSOLETE]         stringstream tmp; tmp << composition[i];
//DX20200724 [OBSOLETE]         number_string = tmp.str();
//DX20200724 [OBSOLETE]       }
//DX20200724 [OBSOLETE]     }
//DX20200724 [OBSOLETE]     is_previous_alpha = isalpha(composition[i]);
//DX20200724 [OBSOLETE]     is_previous_digit = isdigit(composition[i]);
//DX20200724 [OBSOLETE]   }
//DX20200724 [OBSOLETE]   if(is_previous_alpha){
//DX20200724 [OBSOLETE]     stoichiometry.push_back(1);
//DX20200724 [OBSOLETE]   }
//DX20200724 [OBSOLETE]   else if(is_previous_digit){
//DX20200724 [OBSOLETE]     stoichiometry.push_back(aurostd::string2utype<uint>(number_string));
//DX20200724 [OBSOLETE]   }
//DX20200724 [OBSOLETE]   return stoichiometry;
//DX20200724 [OBSOLETE] }
//DX20200724 [OBSOLETE] //DX composition2stoichiometry - 20181009 - END

// ***************************************************************************
// getAtomIndicesByType() //DX20210322
// ***************************************************************************
vector<uint> getAtomIndicesByType(const xstructure& xstr, int type) {

  // Get the atom indices of a given type from an xstructure

  stringstream message;

  uint natoms = xstr.atoms.size();

  vector<uint> indices_atoms_subset;
  for(uint i=0; i<natoms; i++) {
    if(xstr.atoms[i].type == type){ indices_atoms_subset.push_back(i); }
  }

  if(indices_atoms_subset.size() == 0){
    message << "No atoms found with type = " << type << ". Check structure.";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
  }

  return indices_atoms_subset;
}

// ***************************************************************************
// getAtomIndicesByName() //DX20210322
// ***************************************************************************
vector<uint> getAtomIndicesByName(const xstructure& xstr, const string& name) {

  // Get the atom indices of a given name/species from an xstructure

  stringstream message;

  uint natoms = xstr.atoms.size();
  string name_clean = KBIN::VASP_PseudoPotential_CleanName(name);

  vector<uint> indices_atoms_subset;
  for(uint i=0; i<natoms; i++) {
    if(xstr.atoms[i].cleanname == name_clean){ indices_atoms_subset.push_back(i); }
  }

  if(indices_atoms_subset.size() == 0){
    message << "No atoms found with name = " << name << " (note, using name_clean=" << name_clean << "). Check structure.";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
  }

  return indices_atoms_subset;
}


// ***************************************************************************
// getLeastFrequentAtomTypes() //DX20210322
// ***************************************************************************
vector<uint> getLeastFrequentAtomTypes(const xstructure& xstr) {

  // The least frequent atom type is the atom type with the smallest
  // concentration in the crystal. The atoms of this type are the minimal set
  // that exhibit the crystal periodicity (useful for finding
  // alternative lattices and translations).
  // Returns ALL LFA types (could be more than one)

  vector<uint> lfa_types; // lfa = least frequent atom

  // find minimum type count
  int type_count_min = aurostd::min(xstr.num_each_type);

  // find all species with this atom count
  for(uint i=0;i<xstr.num_each_type.size();i++){
    if(xstr.num_each_type[i] == type_count_min){ lfa_types.push_back(i); }
  }

  if(lfa_types.size() == 0){
    throw aurostd::xerror(__AFLOW_FILE__,"getLeastFrequentAtomTypes():","Least frequent atom type not found. Bad xstructure.",_INPUT_ERROR_);
  }

  return lfa_types;
}

// ***************************************************************************
// getLeastFrequentAtomSpecies() //DX20201230 - moved from XtalFinder
// ***************************************************************************
vector<string> getLeastFrequentAtomSpecies(const xstructure& xstr, bool clean) {

  // The least frequent atom species is the species with the smallest
  // concentration in the crystal. These atoms are the minimal set
  // that exhibit the crystal periodicity (useful for finding
  // alternative lattices and translations).
  // clean: cleans atom name (removes pseudopotential)
  // Returns ALL LFA species (could be more than one)

  vector<string> lfa_species; // lfa = least frequent atom

  // find minimum type count
  int type_count_min = aurostd::min(xstr.num_each_type);

  // find the first species with this atom count
  for(uint i=0;i<xstr.num_each_type.size();i++){
    if(xstr.num_each_type[i] == type_count_min){
      if(clean){ lfa_species.push_back(KBIN::VASP_PseudoPotential_CleanName(xstr.species[i])); }
      else{ lfa_species.push_back(xstr.species[i]); }
    }
  }

  if(lfa_species.size() == 0){
    throw aurostd::xerror(__AFLOW_FILE__,"getLeastFrequentAtomTypes():","Least frequent atom species not found. Bad xstructure.",_INPUT_ERROR_);
  }

  return lfa_species;
}

// **************************************************************************
// Function xstructure::GetElements() //DX20200728
// **************************************************************************
vector<string> xstructure::GetElements(bool clean_name, bool fake_names) const{ // Made function const //SD20220221

  // Returns the elements in the xstructure

  bool LDEBUG=(FALSE || XHOST.DEBUG);

  // ---------------------------------------------------------------------------
  // 1) try xstructure.species
  if(species.size() != 0){ //DX20210315 
    if(clean_name){
      vector<string> vspecies;
      for(uint i=0;i<species.size();i++){
        vspecies.push_back(KBIN::VASP_PseudoPotential_CleanName(species[i]));
      }
      return vspecies;
    }
    else{ return aurostd::deque2vector((*this).species); }
  }
  // ---------------------------------------------------------------------------
  // 2) try element names (check if first atom name is given)
  else if (!atoms[0].name.empty()){
    return GetElementsFromAtomNames(clean_name);
  }
  // ---------------------------------------------------------------------------
  // 3) if all are empty, return fake elements (optional) // Does not change xstructure //SD20220221 
  else if (atoms[0].name.empty() && fake_names){
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " WARNING: Atoms are not labeled" << endl;}
    return pflow::getFakeElements(num_each_type.size());
  }

  throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "There are no element names in the structure.",_RUNTIME_ERROR_);
}

// **************************************************************************
// Function xstructure::GetElements() //DX20200728
// **************************************************************************
vector<string> xstructure::GetElementsFromAtomNames(bool clean_name) const{ // Made function const //SD20220221

  // Extracts the species from the atom names


  vector<string> species_found;
  if(atoms.size()==0){ return species_found; }
  if(!atoms[0].name_is_given) { return species_found; }

  uint iat=0;
  string species_tmp = "";
  for(uint i=0;i<num_each_type.size();i++){
    species_tmp = atoms[iat].name; //always the first in the species set
    if(species_tmp.empty()){
      stringstream message;
      message << "Found 1 non-specified species within structure:" << endl;
      message << (*this) << endl;
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ERROR_);
    }
    for(int j=0;j<num_each_type[i];j++){
      // check all atoms of the same type have the same name
      if(atoms[iat].name!=species_tmp){
        stringstream message;
        message << "The number of each type and atom names do not agree." << endl;
        message << (*this) << endl;
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ERROR_);
      }
      iat++;
    }
    if(clean_name){ species_found.push_back(KBIN::VASP_PseudoPotential_CleanName(species_tmp)); }
    else{ species_found.push_back(species_tmp); }
  }
  return species_found;
}

// **************************************************************************
// Function xstructure::GetReducedComposition() //DX20200728
// **************************************************************************
vector<uint> xstructure::GetReducedComposition(bool numerical_sort){


  vector<uint> composition;
  for(uint i=0;i<num_each_type.size();i++){composition.push_back((uint)num_each_type[i]);}
  vector<uint> reduced_composition;

  // sort in numerical order (useful for prototypes)
  if(numerical_sort){ std::stable_sort(composition.begin(),composition.end()); }

  // reduce by GCD
  aurostd::reduceByGCD(composition, reduced_composition);

  return reduced_composition;
}

// ***************************************************************************
// Function getCentroidOfStructure() //DX20200728
// ***************************************************************************
xvector<double> getCentroidOfStructure(const xstructure& xstr, bool use_cpos, bool use_atom_mass){
  return getCentroidOfStructure(xstr.atoms,use_cpos,use_atom_mass);
}

// useful for molecules/non-periodic structures too
xvector<double> getCentroidOfStructure(const deque<_atom>& atoms, bool use_cpos, bool use_atom_mass){

  // Calculate centroid in a structure
  // This overload is useful for non-periodic structures as well
  // (e.g., needed for symmetry analysis of molecules)
  // Uses aurostd::getCentroid()

  // ---------------------------------------------------------------------------
  // coordinate type
  vector<xvector<double> > coordinates;
  // cpos
  if(use_cpos){
    for(uint i=0;i<atoms.size();i++){ coordinates.push_back(atoms[i].cpos); }
  }
  // fpos
  else{
    for(uint i=0;i<atoms.size();i++){ coordinates.push_back(atoms[i].fpos); }
  }

  // ---------------------------------------------------------------------------
  // coordinate weights
  vector<double> weights;
  // uses mass of atomic species
  if(use_atom_mass){
    for(uint i=0;i<atoms.size();i++){ weights.push_back(atoms[i].mass); }
  }
  // consider as geometric points
  else{
    for(uint i=0;i<atoms.size();i++){ weights.push_back(1.0); }
  }

  return getCentroid(coordinates,weights);
}

// ***************************************************************************
// Function getCentroidOfStructurePBC() //DX20200728
// ***************************************************************************
xvector<double> getCentroidOfStructurePBC(const xstructure& xstr, bool use_cpos, bool use_atom_mass){
  return getCentroidOfStructurePBC(xstr.atoms,xstr.lattice,use_cpos,use_atom_mass);
}

xvector<double> getCentroidOfStructurePBC(const deque<_atom>& atoms,
    xmatrix<double> lattice,
    bool use_cpos,
    bool use_atom_mass){

  // Calculate centroid in a structure with periodic boundary conditions
  // Uses aurostd::getCentroidPBC()

  xmatrix<double> cell;

  // ---------------------------------------------------------------------------
  // coordinate type
  vector<xvector<double> > coordinates;
  // cpos
  if(use_cpos){
    for(uint i=0;i<atoms.size();i++){ coordinates.push_back(atoms[i].cpos); }
    cell = lattice; // use lattice
  }
  // fpos
  else{
    for(uint i=0;i<atoms.size();i++){ coordinates.push_back(atoms[i].fpos); }
    cell = aurostd::eye<double>(3,3); // use unit cube
  }

  // ---------------------------------------------------------------------------
  // coordinate weights
  vector<double> weights;
  // uses mass of atomic species
  if(use_atom_mass){
    for(uint i=0;i<atoms.size();i++){ weights.push_back(atoms[i].mass); }
  }
  // consider as geometric points
  else{
    for(uint i=0;i<atoms.size();i++){ weights.push_back(1.0); }
  }

  return getCentroidPBC(coordinates,weights,cell);
}

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// _SYM_OP
// look into aflow.h for the definitions

// constructors
_sym_op::_sym_op() {free();}
_sym_op::_sym_op(const _sym_op& b){copy(b);}

// destructor
_sym_op::~_sym_op() {free();}

void _sym_op::free() {
  Uc.clear();Uf.clear();           // clear stuff
  generator.clear();               // clear stuff
  generator_coefficients.clear();  // clear stuff       //DX20171206 - generator coefficients
  SU2_matrix.clear();              // clear stuff	      //DX20180115 - 2x2 complex SU(2) matrix
  su2_coefficients.clear();        // clear stuff       //DX20180115 - su(2) coefficients on Pauli matrices
  angle=0.0;                       // clear stuff
  axis.clear();                    // clear stuff
  quaternion_vector.clear();       // clear stuff	//GG
  quaternion_matrix.clear();       // clear stuff	//GG
  str_type="";                     // clear stuff
  str_Hermann_Mauguin="";          // clear stuff
  str_Schoenflies="";              // clear stuff
  flag_inversion=FALSE;            // clear stuff
  is_pgroup=FALSE;                 // clear stuff
  is_pgroup_xtal=FALSE;            // clear stuff
  is_pgroupk_Patterson=FALSE;       // clear stuff      //DX20200129
  is_pgroupk=FALSE;                // clear stuff
  is_pgroupk_xtal=FALSE;           // clear stuff       //DX20171205
  ctau.clear();ftau.clear();       // clear stuff
  basis_atoms_map.clear();         // clear stuff
  basis_types_map.clear();         // clear stuff
  basis_map_calculated=FALSE;      // clear stuff
  is_fgroup=FALSE;                 // clear stuff
  ctrasl.clear();ftrasl.clear();   // clear stuff
  is_sgroup=FALSE;                 // clear stuff
  site=0;                          // clear stuff       //DX20170803
  is_agroup=FALSE;                 // clear stuff
}

void _sym_op::copy(const _sym_op& b){
  Uc=b.Uc;
  Uf=b.Uf;
  generator=b.generator;
  generator_coefficients=b.generator_coefficients; //DX20171206 - generator coefficients
  SU2_matrix=b.SU2_matrix;                         //DX20180115 - 2x2 complex SU(2) matrix 
  su2_coefficients=b.su2_coefficients;             //DX20180115 - su(2) coefficients on Pauli matrices
  angle=b.angle;
  axis=b.axis;
  quaternion_vector=b.quaternion_vector;	//GG
  quaternion_matrix=b.quaternion_matrix;	//GG
  str_type=b.str_type;
  str_Hermann_Mauguin=b.str_Hermann_Mauguin;
  str_Schoenflies=b.str_Schoenflies;
  flag_inversion=b.flag_inversion;
  is_pgroup=b.is_pgroup;
  is_pgroup_xtal=b.is_pgroup_xtal;
  is_pgroupk_Patterson=b.is_pgroupk_Patterson; //DX20200129
  is_pgroupk=b.is_pgroupk;
  is_pgroupk_xtal=b.is_pgroupk_xtal;  //DX20171205
  ctau=b.ctau;
  ftau=b.ftau;
  basis_atoms_map.clear();for(uint i=0;i<b.basis_atoms_map.size();i++){basis_atoms_map.push_back(b.basis_atoms_map.at(i));}
  basis_types_map.clear();for(uint i=0;i<b.basis_types_map.size();i++){basis_types_map.push_back(b.basis_types_map.at(i));}
  basis_map_calculated=b.basis_map_calculated;
  is_fgroup=b.is_fgroup;
  ctrasl=b.ctrasl;
  ftrasl=b.ftrasl;
  is_sgroup=b.is_sgroup;
  site=b.site;                        //DX20170803
  is_agroup=b.is_agroup;
}

const _sym_op& _sym_op::operator=(const _sym_op& b) {       // operator=
  if(this != &b) {free();copy(b);}
  return *this;
}

ostream& operator<<(ostream& oss,const _sym_op& symop) {
  xmatrix<double> Uexp(3,3);
  // oss.setf(std::ios::fixed,std::ios::floatfield);
  // oss.precision(10);
  if(symop.is_pgroup==TRUE) oss << " pgroup" << endl;
  if(symop.is_pgroup_xtal==TRUE) oss << " pgroup_xtal" << endl;
  if(symop.is_pgroupk_Patterson==TRUE) oss << " pgroupk_Patterson" << endl; //DX20200129
  if(symop.is_fgroup==TRUE) oss << " fgroup" << endl;
  if(symop.is_sgroup==TRUE) oss << " sgroup" << endl;
  if(symop.is_agroup==TRUE) oss << " agroup" << endl;
  if(symop.is_pgroupk==TRUE) oss << " pgroupk" << endl;
  if(symop.is_pgroupk_xtal==TRUE) oss << " pgroupk_xtal" << endl;  //DX20171205
  if(symop.is_agroup==TRUE) oss << " site:" << symop.site << endl; //DX20170803
  oss << " type: "<< symop.str_type << endl;
  oss << " Hermann_Mauguin: "<< symop.str_Hermann_Mauguin << endl;
  oss << " Schoenflies: "<< symop.str_Schoenflies << endl;
  //DX+CO START
  oss << "" << roundoff(symop.Uc,1e-8) << " Uc "<< endl; //CO roundoff for printing
  //oss << "" << symop.Uc << " Uc "<< endl;
  //DX+CO END
  Uexp=exp(symop.generator);
  // if(! a.symop_inversion[k] ) oss << "" << max(aurostd::abs(Uexp-Uc))
  // << " error " << endl; else oss << "" <<  max(aurostd::abs(Uexp+Uc)) << " error " << endl;
  //oss << "" << xint(symop.Uf) << " Uf "<< endl; //CO, WILL ZERO OUT 0.999999999999999999
  //CO START
  xmatrix<int> Uf_int(symop.Uf.lrows,symop.Uf.lcols,symop.Uf.urows,symop.Uf.ucols);
  for (int i=symop.Uf.lrows;i<=symop.Uf.urows;i++){
    for (int j=symop.Uf.lcols;j<=symop.Uf.ucols;j++){
      Uf_int[i][j]=nint(symop.Uf[i][j]);
    }
  }
  oss << "" << Uf_int << " Uf "<< endl;
  //CO END
  oss << "" << symop.generator << " A=generator U=+-exp(A) [not Uc and -1 if inversion]" << endl;
  //DX20171206
  oss << "" << symop.generator_coefficients << "  so(3) expansion coefficients on Lx, Ly, and Lz basis" << endl;
  //DX20171206
  //DX20180115 - adding SU(2) and su(2); specific xcomplex printing - START
  char buf11_re[80],buf11_im[80],buf12_re[80],buf12_im[80],buf21_re[80],buf21_im[80],buf22_re[80],buf22_im[80];
  string iobuf="%11.4le";
  //HE20221102 switching from deprecated sprintf to snprintf (eliminates the chance of buffer overflows)
  snprintf(buf11_re,80,iobuf.c_str(),symop.SU2_matrix(1,1).re);snprintf(buf11_im,80,iobuf.c_str(),symop.SU2_matrix(1,1).im);
  snprintf(buf12_re,80,iobuf.c_str(),symop.SU2_matrix(1,2).re);snprintf(buf12_im,80,iobuf.c_str(),symop.SU2_matrix(1,2).im);
  snprintf(buf21_re,80,iobuf.c_str(),symop.SU2_matrix(2,1).re);snprintf(buf21_im,80,iobuf.c_str(),symop.SU2_matrix(2,1).im);
  snprintf(buf22_re,80,iobuf.c_str(),symop.SU2_matrix(2,2).re);snprintf(buf22_im,80,iobuf.c_str(),symop.SU2_matrix(2,2).im);
  oss << " (" << buf11_re << "," << buf11_im << ") ";
  oss << " (" << buf12_re << "," << buf12_im << ")" << endl; 
  oss << " (" << buf21_re << "," << buf21_im << ") "; 
  oss << " (" << buf22_re << "," << buf22_im << ")" << "  SU(2) complex matrix [(real,imaginary)]" << endl; 
  //DX - formatting issues with xcomplex: oss << " " << symop.SU2_matrix << "  SU(2) complex matrix [(real,imaginary)]" << endl;
  char buf1_re[80],buf1_im[80],buf2_re[80],buf2_im[80],buf3_re[80],buf3_im[80];
  snprintf(buf1_re,80,iobuf.c_str(),symop.su2_coefficients(1).re);snprintf(buf1_im,80,iobuf.c_str(),symop.su2_coefficients(1).im);
  snprintf(buf2_re,80,iobuf.c_str(),symop.su2_coefficients(2).re);snprintf(buf2_im,80,iobuf.c_str(),symop.su2_coefficients(2).im);
  snprintf(buf3_re,80,iobuf.c_str(),symop.su2_coefficients(3).re);snprintf(buf3_im,80,iobuf.c_str(),symop.su2_coefficients(3).im);
  oss << " (" << buf1_re << "," << buf1_im << ") ";
  oss << " (" << buf2_re << "," << buf2_im << ") ";
  oss << " (" << buf3_re << "," << buf3_im << ")" << "  su(2) expansion coefficients on Pauli matrices [(real,imaginary)]" << endl;
  //DX - formatting issues with complex: oss << " " << symop.su2_coefficients << "  su(2) expansion coefficients on Pauli matrices [(real,imaginary)]" << endl;
  //DX20180115 - adding SU(2) and su(2); specific xcomplex printing - END
  oss << " "<< symop.angle << "  angle " << endl;
  oss << "" << symop.axis <<  "  axis " << endl;
  //GG START HERE
  // Quaternion Output
  // Roundoff used to print zeros instead of values e-17
  oss << "" << roundoff(symop.quaternion_vector, 1e-8) << " quaternion_vector " << endl;
  oss << "" << roundoff(symop.quaternion_matrix, 1e-8) << " quaternion_matrix " << endl;
  //GG STOP HERE
  oss << " "<< symop.flag_inversion << "   inversion " << endl;
  //   oss << " "<< symop.str_type << endl;
  //   oss << " "<< symop.str_Hermann_Mauguin << endl;
  //   oss << " "<< symop.str_Schoenflies << endl;
  if(symop.is_fgroup || symop.is_sgroup) oss << "" << symop.ctau << "    ctau " << endl;
  if(symop.is_fgroup || symop.is_sgroup) oss << "" << symop.ftau << "    ftau " << endl;
  if(symop.is_sgroup) oss << "" << symop.ctrasl << "    ctrasl " << endl;
  if(symop.is_sgroup) oss << "" << symop.ftrasl << "    ftrasl " << endl;

  //DX+CO START
  //if(symop.is_fgroup==TRUE||symop.is_agroup==TRUE) {
  if(symop.basis_map_calculated) {
    //DX+CO END
    oss << " - ";
    for(uint n=0;n<symop.basis_atoms_map.size();n++)
      oss << symop.basis_atoms_map[n] << " ";
    oss << "   basis_atoms_map" << " ";
    oss << " - ";
    for(uint n=0;n<symop.basis_types_map.size();n++)
      oss << symop.basis_types_map[n] << " ";
    oss << "   basis_types_map" << " ";
    oss << endl;
  }

  return oss;
}

void _sym_op::setUc(const xmatrix<double>& _Uc,const xmatrix<double>& lattice){ //CO20190520
  Uc=_Uc;
  xmatrix<double> f2c=trasp(lattice);xmatrix<double> c2f=inverse(trasp(lattice));
  Uf=c2f*Uc*f2c;
}
void _sym_op::setUf(const xmatrix<double>& _Uf,const xmatrix<double>& lattice){ //CO20190520
  Uf=_Uf;
  xmatrix<double> f2c=trasp(lattice);xmatrix<double> c2f=inverse(trasp(lattice));
  Uc=f2c*Uf*c2f;
}

//DX201801107 - add _kpoint class - START
// ***************************************************************************
// ***************************************************************************
// _kpoint

// constructors
_kpoint::_kpoint() {
  iomode=IOAFLOW_AUTO;
  klattice.clear();                // clear stuff
  fpos.clear();cpos.clear();       // clear stuff
  label="";                        // clear stuff
  is_transformed=FALSE;            // clear stuff
}

// destructor
_kpoint::~_kpoint() {
  free();
}

void _kpoint::free() {
}

// assignment operator
const _kpoint& _kpoint::operator=(const _kpoint& b) {       // operator=
  if(this != &b) {
    free();
    iomode=b.iomode;
    fpos=b.fpos;
    cpos=b.cpos;
    label=b.label;
    is_transformed=b.is_transformed;
  }
  return *this;
}

// print kpoint string
string _kpoint::str() const { //CO20220611 - negative sign formatting
  //[CO20220611 - adding spacing for minus]string tmp = "   " + aurostd::joinWDelimiter(xvecDouble2vecString(fpos,4,true,1e-4,FIXED_STREAM),"   ") + "   ! " + label;
  stringstream oss;
  oss.precision(4);
  oss << "   ";
  for(int i=fpos.lrows;i<=fpos.urows;i++){
    if(fpos[i]<10) oss << " ";  //will never happen
    if(!std::signbit(fpos[i])) oss << " ";
    oss << aurostd::utype2string(fpos[i],4,true,1e-4,FIXED_STREAM);
    if(i<fpos.urows) oss << "   ";
  }
  oss << "   ! " << label;
  if(is_transformed){ oss << "\'"; } //add prime
  return oss.str();
}

// operator<<
ostream& operator<<(ostream& oss,const _kpoint& kpt) {
  oss << kpt.str();
  return oss;
}

// transform kpoint
void _kpoint::TransformKpoint(const xmatrix<double>& P){
  fpos=fpos*P;
  klattice=aurostd::inverse(P)*klattice; //i.e., klattice'=Q*klattice
  is_transformed=true;
}
//DX201801107 - add _kpoint class - END

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// wyckoffsite_ITC
// look into aflow.h for the definitions
// constructors

wyckoffsite_ITC::wyckoffsite_ITC() {
  coord.clear();
  index=0; //DX20200427
  type="";
  wyckoffSymbol="";
  letter=""; //DX20180128 - add Wyckoff letter
  site_symmetry=""; //DX20180128 - add Wyckoff site symmetry
  multiplicity=0; //DX20180128 - add Wyckoff multiplicity
  site_occupation=1.0; //DX20191010 - add site occupation (default: 1.0)
  equations.clear(); //DX20180128 - add Wyckoff multiplicity
  parameter_index=0; //DX20200513
}

// destructor
wyckoffsite_ITC::~wyckoffsite_ITC() {
  free();
}

// empty free
void wyckoffsite_ITC::free() {
}

// operator=
const wyckoffsite_ITC& wyckoffsite_ITC::operator=(const wyckoffsite_ITC& b) {       // operator=
  if(this != &b) {
    free();
    coord=b.coord;
    index=b.index; //DX20200501
    type=b.type;
    wyckoffSymbol=b.wyckoffSymbol;
    letter=b.letter; //DX20180128 - add Wyckoff letter
    site_symmetry=b.site_symmetry; //DX20180128 - add Wyckoff site symmetry
    multiplicity=b.multiplicity; //DX20180128 - add Wyckoff multiplicity
    site_occupation=b.site_occupation; //DX20191010 - add site occupation
    equations=b.equations; //DX20180128 - add Wyckoff multiplicity
    parameter_index=b.parameter_index; //DX20200513
  }
  return *this;
}

//DX20190130 - add comparison operator so we can sort by Wyckoff letter, then species - START
// operator<
bool wyckoffsite_ITC::operator<(const wyckoffsite_ITC& b) const {       // operator<
  if(letter<b.letter){
    return true;
  }
  else if(letter==b.letter){
    if(type<b.type){
      return true;
    }
    else {
      return false;
    }
  }
  return false;
}
//DX20190130 - add comparison operator so we can sort by Wyckoff letter, then species - END

// copy 
wyckoffsite_ITC::wyckoffsite_ITC(const wyckoffsite_ITC& b) {
  free();
  coord=b.coord;
  index=b.index; //DX20200501
  type=b.type;
  wyckoffSymbol=b.wyckoffSymbol;
  letter=b.letter; //DX20180128 - add Wyckoff letter
  site_symmetry=b.site_symmetry; //DX20180128 - add Wyckoff site symmetry
  multiplicity=b.multiplicity; //DX20180128 - add Wyckoff multiplicity
  site_occupation=b.site_occupation; //DX20191010 - add site occupation
  equations=b.equations; //DX20180128 - add Wyckoff multiplicity
  parameter_index=b.parameter_index; //DX20200513
}

// operator <<
ostream& operator<<(ostream& oss,const wyckoffsite_ITC& site) {
  //DX20171212 [OBSOLETE] oss << "wyckoffsite_ITC operator<< " << endl;
  oss << " coord: "<< site.coord << endl;
  oss << " index: "<< site.index << endl; //DX20200501
  oss << " type: "<< site.type << endl;
  oss << " letter: "<< site.letter << endl;
  oss << " site_symmetry: "<< site.site_symmetry << endl;
  oss << " multiplicity: "<< site.multiplicity << endl;
  oss << " wyckoffSymbol: "<< site.wyckoffSymbol << endl;
  oss << " site_occupation: " << site.site_occupation << endl; //DX20191010 - add site occupation
  oss << " parameter_index: " << site.parameter_index << endl; //DX20200513
  oss << " equations: " << endl; //DX20191010 - add site occupation
  for(uint i=0;i<site.equations.size();i++){
    oss << "  " << aurostd::joinWDelimiter(site.equations[i],",") << endl;
  }
  return oss;
}

// ***************************************************************************
// AtomEnvironment Class - DX20191122
// ***************************************************************************
// ---------------------------------------------------------------------------
// AtomEnvironment (constructor)
AtomEnvironment::AtomEnvironment(){
  free();
}

// ---------------------------------------------------------------------------
// AtomEnvironment::free
void AtomEnvironment::free(){
  mode=0;
  element_center="";
  type_center=0;
  num_types=0;
  num_neighbors=0;
  elements_neighbor.clear();
  types_neighbor.clear();
  distances_neighbor.clear();
  coordinations_neighbor.clear();
  coordinates_neighbor.clear();
  facets.clear();
  facet_area.clear();
  area=0;
  volume=0;
  has_hull=false;
  facet_order.clear();
  facet_order.resize(8, 0);
}

// ---------------------------------------------------------------------------
// AtomEnvironment (destructor)
AtomEnvironment::~AtomEnvironment(){
  free();
}

// ---------------------------------------------------------------------------
// AtomEnvironment (copy constructor)
AtomEnvironment::AtomEnvironment(const AtomEnvironment& b){
  copy(b);
}

// ---------------------------------------------------------------------------
// AtomEnvironment::copy
void AtomEnvironment::copy(const AtomEnvironment& b) {
  mode=b.mode;
  element_center=b.element_center;
  type_center=b.type_center;
  num_neighbors=b.num_neighbors;
  num_types=b.num_types;
  elements_neighbor=b.elements_neighbor;
  types_neighbor=b.types_neighbor;
  distances_neighbor=b.distances_neighbor;
  coordinations_neighbor=b.coordinations_neighbor;
  coordinates_neighbor=b.coordinates_neighbor;
  facets=b.facets;
  facet_order=b.facet_order;
  facet_area=b.facet_area;
  area=b.area;
  volume=b.volume;
  has_hull=b.has_hull;
}

// ---------------------------------------------------------------------------
// AtomEnvironment::operator=
const AtomEnvironment& AtomEnvironment::operator=(const AtomEnvironment& b){
  if(this!=&b){
    copy(b);
  }
  return *this;
}

// ---------------------------------------------------------------------------
// AtomEnvironment::operator<< 
ostream& operator<<(ostream& oss, const AtomEnvironment& AtomEnvironment){

  oss << AtomEnvironment.toJSON(false).toString();

  return oss;
}

// ***************************************************************************
// AtomEnvironment::constructAtomEnvironmentHull() - HE20210408
// ***************************************************************************

/// @brief constructed a convex hull around the atomic environment
void AtomEnvironment::constructAtomEnvironmentHull(void){
  bool LDEBUG=(false || XHOST.DEBUG);
  if (has_hull) {
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " AE hull is already set" << endl;
    return;
  }

  vector<xvector<double> > points;
  for (uint t = 0; t < num_neighbors; t++) {
    points.push_back(index2Point(t));
  }

  if (LDEBUG) cerr << __AFLOW_FUNC__ << " create AE hull around " << num_neighbors << " atoms" << endl;
  xoption hull_options;
  hull_options.flag("CHULL::FULL_HULL", true);
  hull_options.flag("CHULL::SKIP_N+1_ENTHALPY_GAIN_ANALYSIS", true);
  hull_options.flag("CHULL::SKIP_STABILITY_CRITERION_ANALYSIS", true);
  hull_options.flag("CHULL::INCLUDE_OUTLIERS", true);
  hull_options.flag("CHULL::SEE_NEGLECT", false);
  chull::ConvexHull AEhull;
  AEhull = chull::ConvexHull(hull_options, points);

  if (LDEBUG) cerr << __AFLOW_FUNC__ << " resulting hull has " << AEhull.m_facets.size() << " raw facets" << endl;

  vector<vector<uint> > facet_collection;
  AEhull.getJoinedFacets(facet_collection);
  for (std::vector<vector<uint> >::const_iterator f = facet_collection.begin(); f != facet_collection.end(); ++f) {
    vector<uint> nf;
    for (std::vector<uint>::const_iterator v = f->begin(); v != f->end(); ++v) nf.push_back(*v);
    facets.push_back(nf);
  }
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " after joining " << facets.size() << " facets are remaining" << endl;

  for (std::vector<vector<uint> >::const_iterator f = facets.begin(); f != facets.end(); ++f) {
    if (f->size()<10) facet_order[f->size()-3]++;
    else facet_order[7]++;
  }
  for (uint t = 0; t < facet_collection.size(); t++) {
    vector<xvector<double> > facet_coords;
    for (std::vector<uint>::const_iterator ind=facet_collection[t].begin(); ind != facet_collection[t].end(); ind++){
      facet_coords.push_back(points[*ind]);
    }
    facet_area.push_back(aurostd::areaPointsOnPlane(facet_coords));
  }
  volume = aurostd::volume(points, facets, true);
  area = aurostd::sum(facet_area);
  has_hull = true;
}

// ***************************************************************************
// AtomEnvironment::index2Point() - HE20210408
// ***************************************************************************

/// @brief lookup function to map flat neighbor index back into element sorted coordinates_neighbor list
/// @param index neighbor index
/// @return neighbor coordinates
xvector<double> AtomEnvironment::index2Point(uint index){
  for(uint i=0;i<coordinates_neighbor.size();i++){
    if (index < coordinations_neighbor[i]) return coordinates_neighbor[i][index];
    else index -= coordinations_neighbor[i];
  }
  throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "index out of bounds", _INDEX_BOUNDS_);
}

// ***************************************************************************
// AtomEnvironment::toJSON() - HE20210408
// ***************************************************************************

/// @brief serialize AtomEnvironment class to json
/// @return json string
aurostd::JSONwriter AtomEnvironment::toJSON(bool full) const{

  aurostd::JSONwriter ae_json;

  // this is the XtalFinder printing mode, we need to retain
  // this printing mode since the work has been published
  // this is flatter than the other printing method //DX20210624
  if(mode == ATOM_ENVIRONMENT_MODE_1){
    // element_center
    ae_json.addString("element_center", element_center);
    // type_center
    ae_json.addNumber("type_center", type_center);
    // elements_neighbor
    ae_json.addVector("elements_neighbor", elements_neighbor);
    // distances_neighbor
    ae_json.addVector("distances_neighbor", distances_neighbor, 8, true);
    // coordinations_neighbor
    ae_json.addVector("coordinations_neighbor", coordinations_neighbor);

    // coordinates_neighbor
    vector<string> coordinates_all;
    for(uint i=0;i<coordinates_neighbor.size();i++){
      vector<string> coordinates_set;
      for(uint j=0;j<coordinates_neighbor[i].size();j++){
        aurostd::JSONwriter neighbor;
        neighbor.mergeRawJSON("["+aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(coordinates_neighbor[i][j],8,true,1e-6),",")+"]");
        coordinates_set.push_back(neighbor.toString(false));
      }
      aurostd::JSONwriter json_tmp;
      json_tmp.mergeRawJSON("["+aurostd::joinWDelimiter(coordinates_set,",")+"]");
      coordinates_all.push_back(json_tmp.toString(false));
    }
    ae_json.addVector("coordinates_neighbor",coordinates_all,false);
  }

  else{
    ae_json.addNumber("ae_mode", mode);
    ae_json.addString("center_element", element_center);
    ae_json.addNumber("center_element_index", type_center);
    ae_json.addNumber("element_count", num_types);

    //    "elements_neighbor"
    if (full) {
      vector<aurostd::JSONwriter> distance_collection;
      for (uint i = 0; i < elements_neighbor.size(); i++) {
        aurostd::JSONwriter distance_element;
        distance_element.addNumber("index", i);
        distance_element.addString("name", elements_neighbor[i]);
        distance_element.addNumber("min_distance", distances_neighbor[i]);
        distance_element.addNumber("coordination", coordinations_neighbor[i]);
        distance_collection.push_back(distance_element);
      }
      ae_json.addVector("neighbor_elements", distance_collection);
    }
    if (has_hull){
      ae_json.addNumber("volume", volume);
      ae_json.addNumber("area", area);
    }
    vector<aurostd::JSONwriter> neighbors;
    uint index=0;
    for (uint i = 0; i < coordinations_neighbor.size(); i++) {
      for (uint k = 0; k < coordinations_neighbor[i]; k++) {
        aurostd::JSONwriter neighbor;
        if (full) {
          neighbor.addNumber("index", index);
          neighbor.addString("element", elements_neighbor[i]);
          neighbor.addNumber("element_index", types_neighbor[i]);
          neighbor.addVector("coordinate", coordinates_neighbor[i][k]);
        }
        else {
          neighbor.addString("element", elements_neighbor[i]);
          neighbor.addVector("coordinate", coordinates_neighbor[i][k]);
        }
        neighbors.push_back(neighbor);
        index++;

      }
    }
    ae_json.addVector("neighbors", neighbors);

    if (has_hull && full) {
      vector <aurostd::JSONwriter> facets_collection;
      for (uint i = 0; i < facets.size(); i++) {
        aurostd::JSONwriter facet_entry;
        facet_entry.addNumber("area", facet_area[i]);
        facet_entry.addVector("vertices", facets[i]);
        facets_collection.push_back(facet_entry);
      }
      ae_json.addVector("facets", facets_collection);
      ae_json.addVector("facet_order", facet_order);
    }
  }
  return ae_json;
}


// ***************************************************************************
// AtomEnvironment::getAtomEnvironment() - DX20191122 
// ***************************************************************************
// determines the atomic environment around a central atom 
// current functionality:
//   - calculates the nearest neighbors by element-type, i.e., minimum coordination shell for a given element-type
//     the neighbor elements, types, distance, coordination, and coordinates are stored in the object
//     only one distance is sto
// preliminary functionality, can/will be expanded in the future
void AtomEnvironment::getAtomEnvironment(const xstructure& xstr, uint center_index, uint ae_mode){
  vector<string> neighbor_elements;
  getAtomEnvironment(xstr, center_index, neighbor_elements, ae_mode);
}

void AtomEnvironment::getAtomEnvironment(const xstructure& xstr, uint center_index, const vector<string>& neighbor_elements, uint ae_mode){

  // ---------------------------------------------------------------------------
  // ATOM_ENVIRONMENT_MODE_1 : default minimum coordination shell
  // [FUTURE] ATOM_ENVIRONMENT_MODE_2 : out to a certain radius
  // [FUTURE] ATOM_ENVIRONMENT_MODE_3 : largest gap in radial distribution function

  // ---------------------------------------------------------------------------
  // get central atom info
  mode=ae_mode;
  for(uint i=0;i<xstr.atoms.size();i++){
    if(i==center_index){
      element_center = xstr.atoms[i].name;
      type_center = xstr.atoms[i].type;
    }
  }

  num_types = xstr.species.size();

  // ---------------------------------------------------------------------------
  // ATOM_ENVIRONMENT_MODE_1 : minimum coordination environment for each type 
  if(mode==ATOM_ENVIRONMENT_MODE_1){
    for(uint i=0;i<xstr.species.size();i++){
      // check if types are restricted, otherwise get closest neighbors by type
      if(aurostd::WithinList(neighbor_elements, xstr.species[i]) || neighbor_elements.empty()){
        uint frequency = 0;
        double min_dist = AUROSTD_MAX_DOUBLE;
        vector<xvector<double> > coordinates;

        // calculate minimum coordination shell to a particular element-type 
        minimumCoordinationShell(xstr, center_index, min_dist, frequency, coordinates, xstr.species[i]);

        // store 
        elements_neighbor.push_back(xstr.species[i]);
        types_neighbor.push_back(i);
        distances_neighbor.push_back(min_dist);
        coordinations_neighbor.push_back(frequency);
        coordinates_neighbor.push_back(coordinates);
        num_neighbors += frequency;
      }
    }
  }

  // ---------------------------------------------------------------------------
  // [FUTURE] ATOM_ENVIRONMENT_MODE_2 : environment out to a given radius 

  // ---------------------------------------------------------------------------
  // [FUTURE] ATOM_ENVIRONMENT_MODE_3 : environment out to largest gap in radial distribution function 
  // i.e., GFA convention

}

// ***************************************************************************
// getAtomEnvironments() - DX20191122
// ***************************************************************************
vector<AtomEnvironment> getAtomEnvironments(const xstructure& xstr, uint mode){

  // Calculate the atomic environments in the structure

  vector<AtomEnvironment> environments;

  for(uint i=0;i<xstr.atoms.size();i++){
    AtomEnvironment env; 
    env.getAtomEnvironment(xstr, i, mode);
    environments.push_back(env);
  }
  return environments;
}

// ***************************************************************************
// writeAtomEnvironments() - HE20210723
// ***************************************************************************

void writeAtomEnvironments(vector<AtomEnvironment> AE, const std::map<string, string> meta_data){

  bool LDEBUG = (false || XHOST.DEBUG);

  vector<aurostd::JSONwriter> ae_collection;
  aurostd::JSONwriter ae_json;
  string file_name = "atomic_environment.json";
  string directory_name = "";
  string file_path = "";
  string file_extension = ".json";

  // set filetype
  filetype ftype = json_ft;
  if(XHOST.vflag_control.flag("PRINT_MODE::JSON")) {
    ftype = json_ft;
    file_extension = ".json";
  }
  else if(XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
    ftype = txt_ft;
    file_extension = ".txt";
  }

  // for now just JSON is supported
  if(ftype != json_ft){
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Just JSON is supported at the moment.", _INPUT_ERROR_);
  }

  // construct file path
  if(XHOST.vflag_control.flag("FILE")) file_name = XHOST.vflag_control.getattachedscheme("FILE");
  else {
    if (meta_data.find("auid")!=meta_data.end()) {
      string auid = meta_data.at("auid");
      if (auid.find("aflow:") != std::string::npos) file_name = auid.substr(6);
      else file_name = auid;
    }
  }

  // ensure that filename has the appropriate extension
  if (!(file_name.size() >= file_extension.size() && 0 == file_name.compare(file_name.size()-file_extension.size(), file_extension.size(), file_extension))){
    file_name += file_extension;
  };

  if(XHOST.vflag_control.flag("DIRECTORY")){
    directory_name = XHOST.vflag_control.getattachedscheme("DIRECTORY");
    aurostd::DirectoryMake(directory_name);
  }

  if (!directory_name.empty()) file_path = directory_name + "/" + file_name;
  else file_path = file_name;

  if(LDEBUG) cerr << __AFLOW_FUNC__ << " Saving " << AE.size() << " atomic environments" << endl;

  for(uint i=0; i<AE.size(); i++) ae_collection.push_back(AE[i].toJSON());

  if (!meta_data.empty()) {
    for (std::map<string, string>::const_iterator meta_entry=meta_data.begin(); meta_entry!=meta_data.end(); ++meta_entry){
      ae_json.addString(meta_entry->first, meta_entry->second);
    }
  }

  ae_json.addVector("atomic_environments", ae_collection);
  aurostd::string2file(ae_json.toString(), file_path, "WRITE");
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " Written to " << file_path << endl;
}



// ***************************************************************************
// getLFAAtomEnvironments() - DX20191122
// ***************************************************************************
vector<AtomEnvironment> getLFAAtomEnvironments(const xstructure& xstr, const string& lfa, const vector<string>& LFAs, uint mode){

  // Calculate the LFA atomic environments
  // i.e., only environments comprised of the least frequent atom (LFA) type
  // use-case: quickly screen for potential isoconfigurational structures (see AFLOW-XtalFinder)

  vector<AtomEnvironment> environments_LFA;

  for(uint i=0;i<xstr.atoms.size();i++){
    if(xstr.atoms[i].name == lfa){
      AtomEnvironment LFA_env; 
      LFA_env.getAtomEnvironment(xstr, i, LFAs, mode);
      environments_LFA.push_back(LFA_env);
    }
  }
  return environments_LFA;
}

// ***************************************************************************
// Reset dims for RadiusSphereLattice() - DX20191122
// ***************************************************************************
void resetLatticeDimensions(const xmatrix<double>& lattice, 
    double radius, 
    xvector<int>& dims,
    vector<xvector<double> >& l1, 
    vector<xvector<double> >& l2, 
    vector<xvector<double> >& l3, 
    vector<int>& a_index, 
    vector<int>& b_index, 
    vector<int>& c_index){

  // resets the lattice dimensions (dims) based on radius
  // generates lattice vectors (l1,l2,l3) right away = speed increase
  // stores dimension indices (a_index,b_index,c_index)
  // new dims explore order : zeroth cell to max dims = speed increase 
  // (can break early if match is found)
  // Useful for finding the closest neihbors or minimum interatomic distances:
  // once we find a neighbor, update/reduce how far we need to search to find
  // a closer neighbor
  //DX create function date: 20190705

  // ---------------------------------------------------------------------------
  // get new dimensions based on radius
  if(radius<=_ZERO_TOL_){ dims[1]=1; dims[2]=1; dims[3]=1; }
  else{ dims=LatticeDimensionSphere(lattice,radius); }
  // ---------------------------------------------------------------------------
  // clear old 
  l1.clear(); l2.clear(); l3.clear();
  a_index.clear(); b_index.clear(); c_index.clear();

  // ---------------------------------------------------------------------------
  // [NEW] - go from zeroth cell out
  // more likely to find match close to origin, why start so far away

  // push back zeroth cell : dims[1]=dims[2]=dims[3]=0
  l1.push_back(0*lattice(1));a_index.push_back(0);
  l2.push_back(0*lattice(2));b_index.push_back(0);
  l3.push_back(0*lattice(3));c_index.push_back(0);

  // push back 1,-1,2,-2,...dims,-dims
  for(int a=1;a<=dims[1];a++){l1.push_back(a*lattice(1));a_index.push_back(a); l1.push_back(-a*lattice(1));a_index.push_back(-a);}
  for(int b=1;b<=dims[2];b++){l2.push_back(b*lattice(2));b_index.push_back(b); l2.push_back(-b*lattice(2));b_index.push_back(-b);}
  for(int c=1;c<=dims[3];c++){l3.push_back(c*lattice(3));c_index.push_back(c); l3.push_back(-c*lattice(3));c_index.push_back(-c);}

}

// ***************************************************************************
// minimumCoordinationShellLatticeOnly() - DX20191122 
// ***************************************************************************
void minimumCoordinationShellLatticeOnly(const xmatrix<double>& lattice,
    double& min_dist, 
    uint& frequency, 
    vector<xvector<double> >& coordinates){

  // determine the minimum coordination shell of the lattice
  // i.e., find the set of closest neighbors to the origin
  // (overload: uses lattice radius) 

  // ---------------------------------------------------------------------------
  // determine necessary search radius
  double radius=RadiusSphereLattice(lattice);

  minimumCoordinationShellLatticeOnly(lattice, min_dist, frequency, coordinates, radius);
}

// ***************************************************************************
// minimumCoordinationShellLatticeOnly() - DX20191122 
// ***************************************************************************
void minimumCoordinationShellLatticeOnly(const xmatrix<double>& lattice,
    double& min_dist, 
    uint& frequency, 
    vector<xvector<double> >& coordinates, 
    double radius){

  // determine the minimum coordination shell of the lattice
  // i.e., find the set of closest neighbors to the origin
  // (overload: instantiates lattice dimension information) 

  // ---------------------------------------------------------------------------
  // instantiate lattice vectors 
  vector<xvector<double> > l1, l2, l3; 
  vector<int> a_index, b_index, c_index;
  xvector<int> dims(3); dims[1]=dims[2]=dims[3]=0; // declare/reset
  resetLatticeDimensions(lattice,radius,dims,l1,l2,l3,a_index,b_index,c_index);

  minimumCoordinationShellLatticeOnly(lattice, dims, l1, l2, l3, 
      a_index, b_index, c_index, 
      min_dist, frequency, coordinates, radius);
}

// ***************************************************************************
// minimumCoordinationShellLatticeOnly() - DX20191122
// ***************************************************************************
void minimumCoordinationShellLatticeOnly(const xmatrix<double>& lattice, 
    xvector<int>& dims,
    vector<xvector<double> >& l1, 
    vector<xvector<double> >& l2, 
    vector<xvector<double> >& l3, 
    vector<int>& a_index, 
    vector<int>& b_index, 
    vector<int>& c_index, 
    double& min_dist, 
    uint& frequency, 
    vector<xvector<double> >& coordinates,
    double radius){

  // determine the minimum coordination shell environment of the lattice
  // i.e., find the set of closest neighbors to the origin
  // stores l1, l2, l3, a_index, b_index, and c_index for external use
  // optional "radius" as enables more control over search space 
  // (and potential speed up, may not need to search as far as the lattice radius)

  xvector<double> tmp;

  // ---------------------------------------------------------------------------
  // reset lattice dimensions 
  resetLatticeDimensions(lattice,radius,dims,l1,l2,l3,a_index,b_index,c_index);

  double relative_tolerance = 10.0; // coordination shell thickness is ten percent of minimum distance

  // ---------------------------------------------------------------------------
  // loop through lattice vectors (stored before-hand in l1,l2,l3)
  for(uint m=0;m<l1.size();m++){
    xvector<double> a_component = l1[m];                  //DX : coord1-coord2+a*lattice(1)
    for(uint n=0;n<l2.size();n++){
      xvector<double> ab_component = a_component + l2[n]; //DX : coord1-coord2+a*lattice(1) + (b*lattice(2))
      for(uint p=0;p<l3.size();p++){
        if(!(m==0 && n==0 && p==0)){
          tmp = ab_component + l3[p];                     //DX : coord1-coord2+a*lattice(1) + (b*lattice(2)) + (c*lattice(3))
          double tmp_mod = aurostd::modulus(tmp);
          // ---------------------------------------------------------------------------
          // if found a new minimum distance and update coordination/frequency and coordinate 
          if(tmp_mod<min_dist){
            // ---------------------------------------------------------------------------
            // if new distance is close to the original it is the same coordination shell (add to coordination)
            // otherwise, reset coordination shell
            //DX - FIXED TOL (bad for undecorated prototypes) - if(aurostd::isequal(tmp_mod,min_dist,0.5)){ frequency+=1; } // within half an Angstrom
            if(aurostd::isequal(tmp_mod,min_dist,(min_dist/relative_tolerance))){ frequency+=1; coordinates.push_back(tmp); } // tenth of min_dist
            else{ frequency=1; coordinates.clear(); coordinates.push_back(tmp); } //initialize
            min_dist=tmp_mod;
            // ---------------------------------------------------------------------------
            // diminishing dims: if minimum distance changed, then we may not need to search as far
            // reset loop and search again based on new minimum distance
            if(!(dims[1]==1 && dims[2]==1 && dims[3]==1)){
              resetLatticeDimensions(lattice,min_dist,dims,l1,l2,l3,a_index,b_index,c_index);
              m=n=p=0;
              frequency=0; //reset
              coordinates.clear(); //DX20210222
            }
          }
          //DX - FIXED TOL (bad for undecorated prototypes) - else if(aurostd::isequal(tmp_mod,min_dist,0.5)){ frequency+=1; } // within half an Angstrom
          else if(aurostd::isequal(tmp_mod,min_dist,(min_dist/relative_tolerance))){ frequency+=1; coordinates.push_back(tmp); } // tenth of min dist
        }
      }
    }
  }
}

// ***************************************************************************
// minimumCoordinationShell() - DX20191122 
// ***************************************************************************
void minimumCoordinationShell(const xstructure& xstr, 
    uint center_index, 
    double& min_dist, 
    uint& frequency, 
    vector<xvector<double> >& coordinates){

  string type = "";

  minimumCoordinationShell(xstr, center_index, min_dist, frequency, coordinates, type);
}

// ***************************************************************************
// minimumCoordinationShell() - DX20191122
// ***************************************************************************
void minimumCoordinationShell(const xstructure& xstr, 
    uint center_index, 
    double& min_dist, 
    uint& frequency, 
    vector<xvector<double> >& coordinates, 
    const string& type){

  // determine the minimum coordination shell environment
  // "type" enables the search of environments by certain elements/types only
  // (e.g., find the neighborhood of oxygen atoms surrounding a magnesium center)

  xvector<double> tmp_coord;

  // ---------------------------------------------------------------------------
  // instantiate lattice vectors 
  vector<xvector<double> > l1, l2, l3; 
  vector<int> a_index, b_index, c_index;
  xvector<int> dims(3); dims[1]=dims[2]=dims[3]=0; // declare/reset
  //resetLatticeDimensions(lattice,radius,dims,l1,l2,l3,a_index,b_index,c_index);

  double relative_tolerance = 10.0; // coordination shell thickness is ten percent of minimum distance

  for(uint ii=0; ii<xstr.atoms.size(); ii++){
    // ---------------------------------------------------------------------------
    // if atom ii is not environment center, find minimum distance between center atom ii's images 
    if(ii!=center_index && (xstr.atoms[ii].name == type || type == "")){ //DX20191105 - added type=="" 
      xvector<double> incell_dist = xstr.atoms[center_index].cpos-xstr.atoms[ii].cpos;
      double incell_mod = aurostd::modulus(incell_dist);
      if(!(dims[1]==1 && dims[2]==1 && dims[3]==1) && incell_mod!=1e9){
        resetLatticeDimensions(xstr.lattice,incell_mod,dims,l1,l2,l3,a_index,b_index,c_index);
      }
      //DX20180423 - running vector in each loop saves computations; fewer duplicate operations
      for(uint m=0;m<l1.size();m++){
        xvector<double> a_component = incell_dist + l1[m];    //DX : coord1-coord2+a*lattice(1)
        for(uint n=0;n<l2.size();n++){
          xvector<double> ab_component = a_component + l2[n]; //DX : coord1-coord2+a*lattice(1) + (b*lattice(2))
          for(uint p=0;p<l3.size();p++){
            tmp_coord = ab_component + l3[p];                 //DX : coord1-coord2+a*lattice(1) + (b*lattice(2)) + (c*lattice(3))
            double tmp_mod = aurostd::modulus(tmp_coord);
            if(tmp_mod<min_dist){
              //DX - FIXED TOL (bad for undecorated prototypes) - if(aurostd::isequal(tmp_mod,min_dist,0.5)){ frequency+=1; } // within half an Angstrom
              if(aurostd::isequal(tmp_mod,min_dist,(min_dist/relative_tolerance))){ frequency+=1; coordinates.push_back(tmp_coord); } // tenth of min_dist
              else{ frequency=1; coordinates.clear(); coordinates.push_back(tmp_coord); } //initialize
              min_dist=tmp_mod;
            }
            //DX - FIXED TOL (bad for undecorated prototypes) else if(aurostd::isequal(tmp_mod,min_dist,0.5)){ frequency+=1; } // within half an Angstrom
            else if(aurostd::isequal(tmp_mod,min_dist,(min_dist/relative_tolerance))){ frequency+=1; coordinates.push_back(tmp_coord); } // tenth of min_dist
          }
        }
      }
    }
    // ---------------------------------------------------------------------------
    // if atom is environment center check its images, but only need to search as 
    // far as min_dist or lattice_radius (whichever is smaller)
    else if(ii==center_index && (xstr.atoms[ii].name == type || type == "")){ //DX20191105 - added type==""
      double lattice_radius=RadiusSphereLattice(xstr.lattice);
      double search_radius=min(lattice_radius,min_dist);

      // ---------------------------------------------------------------------------
      // use variant that stores the lattice dimension information so it can be 
      // updated for the "minimumCoordinationShell" function
      minimumCoordinationShellLatticeOnly(xstr.lattice, dims, l1, l2, l3, 
          a_index, b_index, c_index, 
          min_dist, frequency, coordinates, search_radius);
    }
  }
}


// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// STRUCTURE
// look into aflow.h for the definitions

void xstructure::free() { //DX20191220 - moved all initializations from constuctor into free()
  iomode=IOAFLOW_AUTO;        // what else could we do right now...
  title=""; //DX20191210
  directory="";
  prototype="";
  info="";
  if(title.length()==0) title="NO_TITLE_GIVEN";
  // num_atoms=0;
  // num_types=0;
  scale=1.0;
  neg_scale=FALSE;
  scale_second=DEFAULT_POCC_SITE_TOL; //DEFAULT_PARTIAL_OCCUPATION_TOLERANCE; //CO20180409
  neg_scale_second=FALSE; //CO20180409
  scale_third.isentry=false;  //CO20170803 - site tol
  scale_third.content_double=DEFAULT_POCC_STOICH_TOL;//DEFAULT_PARTIAL_OCCUPATION_TOLERANCE;  //CO20170803 - site tol
  coord_type[0]=coord_type[1]=0;
  coord_flag=_COORDS_FRACTIONAL_; // _COORDS_FRACTIONAL_ (0) fractional, _COORDS_CARTESIAN_ (1) cartesian.
  isd=FALSE; // !=0 => Selective dynamics, =0 => no selective dynamics.
  lattice.clear();
  lattice(1,1)=lattice(2,2)=lattice(3,3)=1.0;
  a=b=c=1.0;alpha=beta=gamma=90.0;
  klattice=ReciprocalLattice(lattice,scale);
  f2c=trasp(lattice);
  c2f=inverse(trasp(lattice));
  symbolic_math_representation_only=FALSE;                   //DX20180618
  constrained_symmetry_calculation=FALSE;                    //DX20180618
  symbolic_math_lattice.clear();                             //DX20180618
  num_parameters=0;                                          // number of parameters ANRL 20180618
  num_lattice_parameters=0;                                  // number of lattice parameters ANRL 20180618
  prototype_parameter_list.clear();                          // prototype parameter list ANRL 20180618
  prototype_parameter_values.clear();                        // prototype parameters values ANRL 20180618
  origin.clear();
  // TOLERANCES ------------------------
  equiv_fpos_epsilon=_EQUIV_FPOS_EPS_; // standard but you can change
  // NUM_EACH_TYPE ---------------------
  //ClearSpecies(); //CO20180420
  num_each_type.clear();    //CO20180420 - ClearSpecies()
  comp_each_type.clear();   //CO20180420 - ClearSpecies()
  stoich_each_type.clear(); //CO20171025 //CO20180420 - ClearSpecies()
  // SPECIES ---------------------------
  species.clear();species_pp.clear();species_pp_type.clear();species_pp_version.clear();species_pp_ZVAL.clear();species_pp_vLDAU.clear();species_volume.clear();species_mass.clear(); //CO20180420 - ClearSpecies()
  is_vasp4_poscar_format=TRUE;
  is_vasp5_poscar_format=FALSE;
  // ATOMS -----------------------------
  atoms.clear();
  // FLAGS -----------------------------
  primitive_calculated=FALSE; //DX20201005
  Niggli_calculated=FALSE;
  Niggli_avoid=FALSE;
  Minkowski_calculated=FALSE;
  Minkowski_avoid=FALSE;
  LatticeReduction_calculated=FALSE;
  LatticeReduction_avoid=FALSE;
  // LATTICE stuff ---------------------
  Standard_Lattice_calculated=FALSE;
  Standard_Lattice_avoid=FALSE;
  Standard_Lattice_primitive=FALSE;
  Standard_Lattice_conventional=FALSE;
  Standard_Lattice_has_failed=FALSE;
  bravais_lattice_type="";
  bravais_lattice_variation_type="";//WSETYAWAN mod
  bravais_lattice_system="";
  bravais_lattice_lattice_type="";
  bravais_lattice_lattice_variation_type="";//WSETYAWAN mod
  bravais_lattice_lattice_system="";
  pearson_symbol="";
  reciprocal_lattice_type="";
  reciprocal_lattice_variation_type="";//WSETYAWAN mod
  volume_changed_original2new=FALSE; //DX20181024
  transform_coordinates_original2new.clear(); //DX20181024
  transform_coordinates_new2original.clear(); //DX20181024
  rotate_lattice_original2new.clear(); //DX20181024
  rotate_lattice_new2original.clear(); //DX20181024
  //reciprocal_conventional_lattice_type="";
  bravais_superlattice_lattice.clear(); //DX20210209
  bravais_superlattice_type="";
  bravais_superlattice_variation_type="";
  bravais_superlattice_system="";
  pearson_symbol_superlattice="";
  // GENERAL PURPOSE LABEL -------------
  label_uint=0;
  label_int=0;
  label_double=0;
  // ORDER PARAMETER -------------------
  order_parameter_structure=FALSE;
  order_parameter_atoms.clear();
  order_parameter_orbit=1; // always orbit of itself
  order_parameter_sum=0;
  // PARTIAL OCCUPATION -------------------
  partial_occupation_flag=FALSE;
  partial_occupation_site_tol=DEFAULT_POCC_SITE_TOL; //DEFAULT_PARTIAL_OCCUPATION_TOLERANCE;     // DEFAULT //CO20180409
  partial_occupation_stoich_tol=DEFAULT_POCC_STOICH_TOL; //DEFAULT_PARTIAL_OCCUPATION_TOLERANCE;   // DEFAULT //CO20180409
  partial_occupation_HNF=0;
  partial_occupation_sublattice.clear();
  // FORCES/POSITIONS ------------------
  qm_calculated=FALSE;
  qm_scale=1.0;
  qm_lattice.clear();
  qm_lattice(1,1)=qm_lattice(2,2)=qm_lattice(3,3)=1.0;
  qm_klattice=ReciprocalLattice(qm_lattice,qm_scale);
  qm_f2c=trasp(qm_lattice);
  qm_c2f=inverse(trasp(qm_lattice));
  qm_origin.clear();
  qm_atoms.clear();
  qm_forces.clear();
  qm_forces_write=FALSE;
  qm_positions.clear();
  qm_positions_write=FALSE;
  qm_E_cell=0.0;qm_dE_cell=0.0;qm_H_cell=0.0;qm_PV_cell=0.0;qm_mag_cell=0.0;qm_P=0.0;
  qm_E_atom=0.0;qm_dE_atom=0.0;qm_H_atom=0.0;qm_PV_atom=0.0;qm_mag_atom=0.0; 
  // KPOINTS ---------------------------
  kpoints_k1=0;kpoints_k2=0;kpoints_k3=0;
  kpoints_s1=0;kpoints_s2=0;kpoints_s3=0;
  kpoints_kmax=0;kpoints_kppra=0;
  kpoints_mode=0;
  kpoints_kscheme="";
  //DX+CO START
  dist_nn_min=AUROSTD_NAN;    //CO
  // SYMMETRY TOLERANCE ----------------------------
  sym_eps=AUROSTD_NAN; //DX
  sym_eps_calculated=false; //DX, this means that it was calculated and set by the symmetry routines
  sym_eps_change_count=0; //DX20180222 - added tolerance count specific to structure
  sym_eps_no_scan=false; //DX20210331 - added no scan specific to structure
  //DX+CO END
  // PGROUP ----------------------------
  pgroup.clear();            // just initialize
  pgroup_calculated=FALSE;
  // PGROUP_XTAL ----------------------------
  pgroup_xtal.clear();        // just initialize
  pgroup_xtal_calculated=FALSE;
  crystal_family="";crystal_system="";point_group_crystal_class="";
  point_group_Shoenflies="";point_group_Hermann_Mauguin="";point_group_orbifold="";
  point_group_type="";point_group_order="";point_group_structure="";
  // PGROUPK_PATTERSON ---------------------------- //DX20200129
  pgroupk_Patterson.clear();        // just initialize
  pgroupk_Patterson_calculated=FALSE;
  // PGROUPK ----------------------------
  pgroupk.clear();            // just initialize
  pgroupk_calculated=FALSE;
  // PGROUPK_XTAL ----------------------------
  pgroupk_xtal.clear();            // just initialize //DX20171205 - Added pgroupk_xtal
  pgroupk_xtal_calculated=FALSE;                      //DX20171205 - Added pgroupk_xtal
  // FGROUP ----------------------------
  fgroup.clear();            // just initialize
  fgroup_calculated=FALSE;
  // SGROUP ----------------------------
  sgroup_radius=-_calculate_symmetry_default_sgroup_radius_; // symmetry not calculated
  sgroup_radius_dims.clear();
  sgroup.clear();            // just initialize
  sgroup_calculated=FALSE;
  // SITE POINT GROUP ------------------
  agroup_calculated=FALSE;
  for(uint i=0;i<agroup.size();i++)
    agroup.at(i).clear();
  agroup.clear();
  // INEQUIVALENT ATOMS ----------------
  iatoms_calculated=FALSE;
  for(uint i=0;i<iatoms.size();i++)
    iatoms.at(i).clear();
  iatoms.clear();
  // SPACE GROUP WITH PLATON/FINDSYM -----------
  spacegroup="";
  spacegrouplabel="";
  spacegroupnumber=0;
  spacegroupnumberoption=0;
  spacegroupoption="";
  is_spacegroup_platon=FALSE;
  is_spacegroup_findsym=FALSE;
  is_spacegroup_aflow=FALSE;
  // SPACE GROUP ITC (RHT)
  crystal_system_ITC="";  //RHT
  point_group_ITC=""; //RHT
  //DX+CO START
  bravais_label_ITC='X';  //RHT
  lattice_label_ITC='X'; //RHT
  space_group_ITC=0; //RHT
  //DX+CO END
  wyckoff_library_entry_ITC=""; //RHT
  wyccar_ITC.clear(); //RHT
  standard_lattice_ITC.clear(); //RHT
  standard_basis_ITC.clear(); //RHT
  wyckoff_sites_ITC.clear();   //RHT
  wyckoff_symbols_ITC.clear(); //RHT
  setting_ITC = 0; //DX20170830 - SGDATA
  origin_ITC.clear(); //DX20170830 - SGDATA
  general_position_ITC.clear(); //DX20170830 - SGDATA
  // GRID ATOMS ------------------------
  grid_atoms_calculated=FALSE;
  grid_atoms_dimsL.clear();
  grid_atoms_dimsH.clear();
  grid_atoms.clear();        // just initialize
  grid_atoms_number=0;
  grid_atoms_sc2pcMap.clear(); //CO20171025
  grid_atoms_pc2scMap.clear(); //CO20171025
  // LIJK OBEJCTS ----------------------
  lijk_calculated=FALSE;
  lijk_table.clear();
  lijk_cpos.clear();
  lijk_fpos.clear();
  lijk_dims.clear();
  // NEIGHBOR ------------------------
  //DX20210122 [OBSOLETE] neighbors_calculated=FALSE;
  //DX20210122 [OBSOLETE] neighbors_radius=0.0;
  //DX20210122 [OBSOLETE] neighbors_dradius=0.0;
  //DX20210122 [OBSOLETE] for(uint i=0;i<neighbors_atoms_func_r_vs_nn.size();i++)
  //DX20210122 [OBSOLETE]   neighbors_atoms_func_r_vs_nn.at(i).clear();
  //DX20210122 [OBSOLETE] neighbors_atoms_func_r_vs_nn.clear();
  //DX20210122 [OBSOLETE] for(uint i=0;i<neighbors_atoms_func_num_vs_nn.size();i++)
  //DX20210122 [OBSOLETE]   neighbors_atoms_func_num_vs_nn.at(i).clear();
  //DX20210122 [OBSOLETE] neighbors_atoms_func_num_vs_nn.clear();
  //DX20210122 [OBSOLETE] neighbors_func_r_vs_nn.clear();                   // contains function distance vs neighbors (all atoms)
  //DX20210122 [OBSOLETE] neighbors_func_num_vs_nn.clear();                    // contains function number vs neighbors (all atoms)
  // OUTPUT/ERROR ----------------------
  Niggli_has_failed=FALSE;
  Minkowski_has_failed=FALSE;
  LatticeReduction_has_failed=FALSE;
  write_lattice_flag=FALSE;
  write_klattice_flag=FALSE;
  write_inequivalent_flag=FALSE;
  write_DEBUG_flag=FALSE;
  error_flag=FALSE;
  error_string="";
  // -----------------------------------
}

// Constructors
xstructure::xstructure(const string& structure_title) { //CO20211122
  (*this).initialize(structure_title);  //CO20211122
  //CO20211122 [OBSOLETE - done in initialize]free(); //DX20191220 - moved contents below into free()
  //CO20211122 [OBSOLETE - done in initialize]title=structure_title;
}

// ifstream/istream
xstructure::xstructure(istream& _input,int _iomode) {
  (*this).initialize(_input, _iomode);
  //DX20210129 [OBSOLETE - done in initialize] free(); //DX20191220 - added free to initialize
  //DX20210129 [OBSOLETE - done in initialize] (*this).iomode=_iomode;
  //DX20210129 [OBSOLETE - done in initialize] _input >> (*this);
}

xstructure::xstructure(ifstream& _input,int _iomode) {
  (*this).initialize(_input, _iomode);
  //DX20210129 [OBSOLETE - done in initialize] free(); //DX20191220 - added free to initialize
  //DX20210129 [OBSOLETE - done in initialize] (*this).iomode=_iomode;
  //DX20210129 [OBSOLETE - done in initialize] _input >> (*this);
}

xstructure::xstructure(const stringstream& __input,int _iomode) { //DX20210129 - added const
  (*this).initialize(__input, _iomode); //DX20210129
  //DX20210129 [OBSOLETE - done in initialize] free(); //DX20191220 - added free to initialize
  //DX20210129 [OBSOLETE - done in initialize] (*this).iomode=_iomode;
  //DX20210129 [OBSOLETE - done in initialize] stringstream _input(__input.str());
  //DX20210129 [OBSOLETE - done in initialize] _input >> (*this);
}

xstructure::xstructure(const string& _input,int _iomode) {  //CO20211122
  (*this).initialize(_input, _iomode);  //CO20211122
  //CO20211122 [OBSOLETE - done in initialize]free(); //DX20191220 - added free to initialize
  //CO20211122 [OBSOLETE - done in initialize]stringstream strstream;
  //CO20211122 [OBSOLETE - done in initialize]aurostd::efile2stringstream(_input,strstream); //CO20171025
  //CO20211122 [OBSOLETE - done in initialize](*this).iomode=_iomode;
  //CO20211122 [OBSOLETE - done in initialize](*this).directory = _input; //DX20180526 - location of xstructure
  //CO20211122 [OBSOLETE - done in initialize]strstream >> (*this);
}

xstructure::xstructure(const string& url,const string& file,int _iomode) {  //CO20211122
  (*this).initialize(url, file, _iomode);  //CO20211122
  //CO20211122 [OBSOLETE - done in initialize]free(); //DX20191220 - added free to initialize
  //CO20211122 [OBSOLETE - done in initialize]stringstream strstream;
  //CO20211122 [OBSOLETE - done in initialize]aurostd::url2stringstream(url+"/"+file,strstream);
  //CO20211122 [OBSOLETE - done in initialize](*this).iomode=_iomode;
  //CO20211122 [OBSOLETE - done in initialize](*this).directory = url+"/"+file; //DX20180526 - location of xstructure
  //CO20211122 [OBSOLETE - done in initialize]strstream >> (*this);
}

void xstructure::copy(const xstructure& bstr) {
  // All the other stuff not set in the constructor
  iomode=bstr.iomode;
  title=bstr.title;
  directory=bstr.directory;
  prototype=bstr.prototype;
  info=bstr.info;
  //  num_types=bstr.num_types;
  //  num_atoms=bstr.num_atoms;
  scale=bstr.scale;
  neg_scale=bstr.neg_scale;
  scale_second=bstr.scale_second; //CO20180409
  neg_scale_second=bstr.neg_scale_second; //CO20180409
  scale_third=bstr.scale_third; //CO20170803 - site tol //CO20180409
  strcpy(coord_type,bstr.coord_type);
  coord_flag=bstr.coord_flag;
  isd=bstr.isd;
  lattice=bstr.lattice;
  a=bstr.a;b=bstr.b;c=bstr.c;
  alpha=bstr.alpha;beta=bstr.beta;gamma=bstr.gamma;
  klattice=bstr.klattice;
  origin=bstr.origin;
  f2c=bstr.f2c;
  c2f=bstr.c2f;
  symbolic_math_representation_only=bstr.symbolic_math_representation_only;                  //DX20180618
  constrained_symmetry_calculation=bstr.constrained_symmetry_calculation;                    //DX20180618
  symbolic_math_lattice=bstr.symbolic_math_lattice; //DX20180618
  num_parameters=bstr.num_parameters;                                      // number of parameters ANRL 20180618
  num_lattice_parameters=bstr.num_lattice_parameters;                      // number of lattice parameters ANRL 20180618
  prototype_parameter_list=bstr.prototype_parameter_list;                  // prototype parameter list ANRL 20180618
  prototype_parameter_values=bstr.prototype_parameter_values;              // prototype parameters values ANRL 20180618
  // TOLERANCES ------------------------
  equiv_fpos_epsilon=bstr.equiv_fpos_epsilon;
  // NUM_EACH_TYPE ---------------------
  num_each_type.clear();
  for(uint i=0;i<bstr.num_each_type.size();i++)
    num_each_type.push_back(bstr.num_each_type.at(i));
  comp_each_type.clear();
  for(uint i=0;i<bstr.comp_each_type.size();i++)
    comp_each_type.push_back(bstr.comp_each_type.at(i));
  stoich_each_type.clear(); //CO20171025
  for(uint i=0;i<bstr.stoich_each_type.size();i++) //CO20171025
    stoich_each_type.push_back(bstr.stoich_each_type.at(i)); //CO20171025
  // SPECIES ---------------------------
  species.clear();
  for(uint i=0;i<bstr.species.size();i++)
    species.push_back(bstr.species.at(i));
  species_pp.clear();
  for(uint i=0;i<bstr.species_pp.size();i++)
    species_pp.push_back(bstr.species_pp.at(i));
  species_pp_type.clear();
  for(uint i=0;i<bstr.species_pp_type.size();i++)
    species_pp_type.push_back(bstr.species_pp_type.at(i));
  species_pp_version.clear();
  for(uint i=0;i<bstr.species_pp_version.size();i++)
    species_pp_version.push_back(bstr.species_pp_version.at(i));
  species_pp_ZVAL.clear();
  for(uint i=0;i<bstr.species_pp_ZVAL.size();i++)
    species_pp_ZVAL.push_back(bstr.species_pp_ZVAL.at(i));
  species_pp_vLDAU.clear();
  for(uint i=0;i<bstr.species_pp_vLDAU.size();i++)
    species_pp_vLDAU.push_back(bstr.species_pp_vLDAU.at(i));
  species_volume.clear();
  for(uint i=0;i<bstr.species_volume.size();i++)
    species_volume.push_back(bstr.species_volume.at(i));
  species_mass.clear();
  for(uint i=0;i<bstr.species_mass.size();i++)
    species_mass.push_back(bstr.species_mass.at(i));
  is_vasp4_poscar_format=bstr.is_vasp4_poscar_format;
  is_vasp5_poscar_format=bstr.is_vasp5_poscar_format;
  // FLAGS -----------------------------
  primitive_calculated=bstr.primitive_calculated; //DX20201007
  Niggli_calculated=bstr.Niggli_calculated;
  Niggli_avoid=bstr.Niggli_avoid;
  Minkowski_calculated=bstr.Minkowski_calculated;
  Minkowski_avoid=bstr.Minkowski_avoid;
  LatticeReduction_calculated=bstr.LatticeReduction_calculated;
  LatticeReduction_avoid=bstr.LatticeReduction_avoid;
  // LATTICE stuff ---------------------
  Standard_Lattice_calculated=bstr.Standard_Lattice_calculated;
  Standard_Lattice_avoid=bstr.Standard_Lattice_avoid;
  Standard_Lattice_primitive=bstr.Standard_Lattice_primitive;
  Standard_Lattice_conventional=bstr.Standard_Lattice_conventional;
  Standard_Lattice_has_failed=bstr.Standard_Lattice_has_failed;
  bravais_lattice_type=bstr.bravais_lattice_type;
  bravais_lattice_variation_type=bstr.bravais_lattice_variation_type;
  bravais_lattice_system=bstr.bravais_lattice_system;
  bravais_lattice_lattice_type=bstr.bravais_lattice_lattice_type;
  bravais_lattice_lattice_variation_type=bstr.bravais_lattice_lattice_variation_type;
  bravais_lattice_lattice_system=bstr.bravais_lattice_lattice_system;
  pearson_symbol=bstr.pearson_symbol;
  reciprocal_lattice_type=bstr.reciprocal_lattice_type;
  reciprocal_lattice_variation_type=bstr.reciprocal_lattice_variation_type;
  bravais_superlattice_lattice=bstr.bravais_superlattice_lattice; //DX20210209
  bravais_superlattice_type=bstr.bravais_superlattice_type;
  bravais_superlattice_variation_type=bstr.bravais_superlattice_variation_type;
  bravais_superlattice_system=bstr.bravais_superlattice_system;
  pearson_symbol_superlattice=bstr.pearson_symbol_superlattice;
  volume_changed_original2new=bstr.volume_changed_original2new; //DX20181024
  transform_coordinates_original2new=bstr.transform_coordinates_original2new; //DX20181024
  transform_coordinates_new2original=bstr.transform_coordinates_new2original; //DX20181024
  rotate_lattice_original2new=bstr.rotate_lattice_original2new; //DX20181024
  rotate_lattice_new2original=bstr.rotate_lattice_new2original; //DX20181024
  // ATOMS -----------------------------
  atoms.clear();
  for(uint i=0;i<bstr.atoms.size();i++)
    atoms.push_back(bstr.atoms.at(i));
  // GENERAL PURPOSE LABEL -------------
  label_uint=bstr.label_uint;
  label_int=bstr.label_int;
  label_double=bstr.label_double;
  // ORDER PARAMETER -------------------
  order_parameter_structure=bstr.order_parameter_structure;
  order_parameter_atoms.clear();
  for(uint i=0;i<bstr.order_parameter_atoms.size();i++)
    order_parameter_atoms.push_back(bstr.order_parameter_atoms.at(i));
  order_parameter_orbit=bstr.order_parameter_orbit;
  order_parameter_sum=bstr.order_parameter_sum;
  // PARTIAL OCCUPATION -------------------
  partial_occupation_flag=bstr.partial_occupation_flag;
  partial_occupation_site_tol=bstr.partial_occupation_site_tol;     //CO20180409
  partial_occupation_stoich_tol=bstr.partial_occupation_stoich_tol; //CO20180409
  partial_occupation_HNF=bstr.partial_occupation_HNF;
  partial_occupation_sublattice.clear();
  for(uint i=0;i<bstr.partial_occupation_sublattice.size();i++)
    partial_occupation_sublattice.push_back(bstr.partial_occupation_sublattice.at(i));
  // FORCES/POSITIONS ------------------
  qm_calculated=bstr.qm_calculated;
  qm_scale=bstr.qm_scale;
  qm_lattice=bstr.qm_lattice;
  qm_klattice=bstr.qm_klattice;
  qm_f2c=bstr.qm_f2c;
  qm_c2f=bstr.qm_c2f;
  qm_origin=bstr.qm_origin;
  qm_atoms.clear();
  for(uint i=0;i<bstr.qm_atoms.size();i++)
    qm_atoms.push_back(bstr.qm_atoms.at(i));
  qm_forces.clear();
  for(uint i=0;i<bstr.qm_forces.size();i++)
    qm_forces.push_back(bstr.qm_forces.at(i));
  qm_forces_write=bstr.qm_forces_write;
  qm_positions.clear();
  for(uint i=0;i<bstr.qm_positions.size();i++)
    qm_positions.push_back(bstr.qm_positions.at(i));
  qm_positions_write=bstr.qm_positions_write;
  qm_E_cell=bstr.qm_E_cell;qm_dE_cell=bstr.qm_dE_cell;qm_H_cell=bstr.qm_H_cell;qm_PV_cell=bstr.qm_PV_cell;qm_mag_cell=bstr.qm_mag_cell;qm_P=bstr.qm_P;
  qm_E_atom=bstr.qm_E_atom;qm_dE_atom=bstr.qm_dE_atom;qm_H_atom=bstr.qm_H_atom;qm_PV_atom=bstr.qm_PV_atom;qm_mag_atom=bstr.qm_mag_atom;
  // KPOINTS ---------------------------
  kpoints_k1=bstr.kpoints_k1;kpoints_k2=bstr.kpoints_k2;kpoints_k3=bstr.kpoints_k3;
  kpoints_s1=bstr.kpoints_s1;kpoints_s2=bstr.kpoints_s2;kpoints_s3=bstr.kpoints_s3;
  kpoints_kmax=bstr.kpoints_kmax;kpoints_kppra=bstr.kpoints_kppra;
  kpoints_mode=bstr.kpoints_mode;
  kpoints_kscheme=bstr.kpoints_kscheme;
  //DX+CO START
  dist_nn_min=bstr.dist_nn_min;    //CO
  // SYMMETRY TOLERANCE ----------------------------
  sym_eps=bstr.sym_eps; //DX
  sym_eps_calculated=bstr.sym_eps_calculated; //DX
  sym_eps_change_count=bstr.sym_eps_change_count; //DX20180222 - added tolerance count specific to structure
  sym_eps_no_scan=bstr.sym_eps_no_scan; //DX20210331 - added no scan specific to structure
  //DX+CO END
  // PGROUP ----------------------------
  pgroup.clear();
  for(uint i=0;i<bstr.pgroup.size();i++)
    pgroup.push_back(bstr.pgroup.at(i));
  pgroup_calculated=bstr.pgroup_calculated;
  // PGROUP_XTAL ----------------------------
  pgroup_xtal.clear();
  for(uint i=0;i<bstr.pgroup_xtal.size();i++)
    pgroup_xtal.push_back(bstr.pgroup_xtal.at(i));
  pgroup_xtal_calculated=bstr.pgroup_xtal_calculated;
  crystal_family=bstr.crystal_family;
  crystal_system=bstr.crystal_system;
  point_group_crystal_class=bstr.point_group_crystal_class;
  point_group_Shoenflies=bstr.point_group_Shoenflies;
  point_group_Hermann_Mauguin=bstr.point_group_Hermann_Mauguin;
  point_group_orbifold=bstr.point_group_orbifold;
  point_group_type=bstr.point_group_type;
  point_group_order=bstr.point_group_order;
  point_group_structure=bstr.point_group_structure;
  // PGROUPK_PATTERSON ---------------------------- //DX20200129
  pgroupk_Patterson.clear();
  for(uint i=0;i<bstr.pgroupk_Patterson.size();i++)
    pgroupk_Patterson.push_back(bstr.pgroupk_Patterson.at(i));
  pgroupk_Patterson_calculated=bstr.pgroupk_Patterson_calculated;
  // PGROUPK ----------------------------
  pgroupk.clear();
  for(uint i=0;i<bstr.pgroupk.size();i++)
    pgroupk.push_back(bstr.pgroupk.at(i));
  pgroupk_calculated=bstr.pgroupk_calculated;
  // PGROUPK_XTAL ----------------------------
  pgroupk_xtal.clear();                                    //DX20171205 - Added pgroupk_xtal
  for(uint i=0;i<bstr.pgroupk_xtal.size();i++)             //DX20171205 - Added pgroupk_xtal
    pgroupk_xtal.push_back(bstr.pgroupk_xtal.at(i));       //DX20171205 - Added pgroupk_xtal
  pgroupk_xtal_calculated=bstr.pgroupk_xtal_calculated;    //DX20171205 - Added pgroupk_xtal
  // FGROUP ----------------------------
  fgroup.clear();
  for(uint i=0;i<bstr.fgroup.size();i++)
    fgroup.push_back(bstr.fgroup.at(i));
  fgroup_calculated=bstr.fgroup_calculated;
  // SGROUP ----------------------------
  sgroup_radius=bstr.sgroup_radius;
  sgroup_radius_dims=bstr.sgroup_radius_dims;
  sgroup.clear();
  for(uint i=0;i<bstr.sgroup.size();i++)
    sgroup.push_back(bstr.sgroup.at(i));
  sgroup_calculated=bstr.sgroup_calculated;
  // SITE POINT GROUP ------------------
  agroup_calculated=bstr.agroup_calculated;
  for(uint i=0;i<agroup.size();i++) agroup.at(i).clear();
  agroup.clear();
  agroup=std::vector<std::vector<_sym_op> > (bstr.agroup.size());
  for(uint i=0;i<bstr.agroup.size();i++)
    for(uint j=0;j<bstr.agroup.at(i).size();j++)
      agroup.at(i).push_back(bstr.agroup.at(i).at(j));
  // INEQUIVALENT ATOMS ----------------
  iatoms_calculated=bstr.iatoms_calculated;
  for(uint i=0;i<iatoms.size();i++) iatoms.at(i).clear();
  iatoms.clear();
  for(uint i=0;i<bstr.iatoms.size();i++) {
    iatoms.push_back(std::vector<int>(0));
    for(uint j=0;j<bstr.iatoms.at(i).size();j++)
      iatoms.at(i).push_back(bstr.iatoms.at(i).at(j));
  }
  // SPACE GROUP WITH PLATON/FINDSYM/AFLOW -----------
  spacegroup=bstr.spacegroup;
  spacegrouplabel=bstr.spacegrouplabel;
  spacegroupoption=bstr.spacegroupoption;
  spacegroupnumber=bstr.spacegroupnumber;
  spacegroupnumberoption=bstr.spacegroupnumberoption;
  is_spacegroup_platon=bstr.is_spacegroup_platon;
  is_spacegroup_findsym=bstr.is_spacegroup_findsym;
  is_spacegroup_aflow=bstr.is_spacegroup_aflow;
  // SPACE GROUP ITC -----------
  crystal_system_ITC=bstr.crystal_system_ITC;  //RHT
  point_group_ITC=bstr.point_group_ITC; //RHT
  bravais_label_ITC=bstr.bravais_label_ITC; //RHT
  lattice_label_ITC=bstr.lattice_label_ITC; //RHT
  space_group_ITC=bstr.space_group_ITC; //RHT
  wyckoff_library_entry_ITC=bstr.wyckoff_library_entry_ITC; //RHT
  wyccar_ITC.clear(); for(uint i=0;i<bstr.wyccar_ITC.size();i++) wyccar_ITC.push_back(bstr.wyccar_ITC.at(i)); //RHT
  standard_lattice_ITC=bstr.standard_lattice_ITC; //RHT
  standard_basis_ITC.clear(); for(uint i=0;i<bstr.standard_basis_ITC.size();i++) standard_basis_ITC.push_back(bstr.standard_basis_ITC.at(i)); //RHT
  wyckoff_sites_ITC.clear(); for(uint i=0;i<bstr.wyckoff_sites_ITC.size();i++) wyckoff_sites_ITC.push_back(bstr.wyckoff_sites_ITC.at(i)); //RHT
  wyckoff_symbols_ITC.clear(); for(uint i=0;i<bstr.wyckoff_symbols_ITC.size();i++) wyckoff_symbols_ITC.push_back(bstr.wyckoff_symbols_ITC.at(i)); //RHT
  setting_ITC=bstr.setting_ITC; //DX20170830 - SGDATA
  origin_ITC=bstr.origin_ITC; //DX20170830 - SGDATA
  general_position_ITC=bstr.general_position_ITC; //DX20170830 - SGDATA
  // GRID ATOMS ------------------------
  grid_atoms_calculated=bstr.grid_atoms_calculated;
  grid_atoms_dimsL=bstr.grid_atoms_dimsL;
  grid_atoms_dimsH=bstr.grid_atoms_dimsH;
  grid_atoms.clear();
  for(uint i=0;i<bstr.grid_atoms.size();i++)
    grid_atoms.push_back(bstr.grid_atoms.at(i));
  grid_atoms_number=bstr.grid_atoms_number;
  grid_atoms_sc2pcMap.clear(); for(uint i=0;i<bstr.grid_atoms_sc2pcMap.size();i++){grid_atoms_sc2pcMap.push_back(bstr.grid_atoms_sc2pcMap[i]);} //CO20171025
  grid_atoms_pc2scMap.clear(); for(uint i=0;i<bstr.grid_atoms_pc2scMap.size();i++){grid_atoms_pc2scMap.push_back(bstr.grid_atoms_pc2scMap[i]);} //CO20171025
  // LIJK OBEJCTS ----------------------
  lijk_calculated=bstr.lijk_calculated;
  lijk_table.clear();
  lijk_cpos.clear();
  lijk_fpos.clear();
  for(uint i=0;i<bstr.lijk_table.size();i++) {
    lijk_table.push_back(bstr.lijk_table.at(i));
    lijk_cpos.push_back(bstr.lijk_cpos.at(i));
    lijk_fpos.push_back(bstr.lijk_fpos.at(i));
  }
  lijk_dims=bstr.lijk_dims;
  // NEIGHBOR ------------------------
  //DX20210122 [OBSOLETE] neighbors_calculated=bstr.neighbors_calculated;
  //DX20210122 [OBSOLETE] neighbors_radius=bstr.neighbors_radius;
  //DX20210122 [OBSOLETE] neighbors_dradius=bstr.neighbors_dradius;
  //  for(uint i=0;i<neighbors_atoms_func_r_vs_nn.size();i++)
  //   neighbors_atoms_func_r_vs_nn.at(i).clear();
  //DX20210122 [OBSOLETE] neighbors_atoms_func_r_vs_nn.clear();
  //DX20210122 [OBSOLETE] for(uint i=0;i<bstr.neighbors_atoms_func_r_vs_nn.size();i++)
  //DX20210122 [OBSOLETE]   neighbors_atoms_func_r_vs_nn.push_back(bstr.neighbors_atoms_func_r_vs_nn.at(i));
  //  for(uint i=0;i<neighbors_atoms_func_num_vs_nn.size();i++)
  //   neighbors_atoms_func_num_vs_nn.at(i).clear();
  //DX20210122 [OBSOLETE] neighbors_atoms_func_num_vs_nn.clear();
  //DX20210122 [OBSOLETE] for(uint i=0;i<bstr.neighbors_atoms_func_num_vs_nn.size();i++)
  //DX20210122 [OBSOLETE]   neighbors_atoms_func_num_vs_nn.push_back(bstr.neighbors_atoms_func_num_vs_nn.at(i));
  //DX20210122 [OBSOLETE] neighbors_func_r_vs_nn.clear();
  //DX20210122 [OBSOLETE] for(uint i=0;i<bstr.neighbors_func_r_vs_nn.size();i++)
  //DX20210122 [OBSOLETE]   neighbors_func_r_vs_nn.push_back(bstr.neighbors_func_r_vs_nn.at(i));
  //DX20210122 [OBSOLETE] neighbors_func_num_vs_nn.clear();
  //DX20210122 [OBSOLETE] for(uint i=0;i<bstr.neighbors_func_num_vs_nn.size();i++)
  //DX20210122 [OBSOLETE]   neighbors_func_num_vs_nn.push_back(bstr.neighbors_func_num_vs_nn.at(i));
  // OUTPUT/ERROR ----------------------
  Niggli_has_failed=bstr.Niggli_has_failed;
  Minkowski_has_failed=bstr.Minkowski_has_failed;
  LatticeReduction_has_failed=bstr.LatticeReduction_has_failed;
  write_lattice_flag=bstr.write_lattice_flag;
  write_klattice_flag=bstr.write_klattice_flag;
  write_inequivalent_flag=bstr.write_inequivalent_flag;
  write_DEBUG_flag=bstr.write_DEBUG_flag;
  error_flag=bstr.error_flag;
  error_string=bstr.error_string;
  // ----------------------------------
}

//ME20200220 - from CO's function in apl::Supercell
void LightCopy(const xstructure& a, xstructure& b) {
  b.clear();
  stringstream POSCAR;
  POSCAR.str("");
  POSCAR << a;
  POSCAR >> b;
  //enable inequivalent flag to work
  for (uint i = 0; i < b.atoms.size(); i++) {
    b.atoms[i].equivalent = a.atoms[i].equivalent;
    b.atoms[i].is_inequivalent = a.atoms[i].is_inequivalent;
    b.atoms[i].num_equivalents = a.atoms[i].num_equivalents;
  }
  //enable inequivalent flag to work
  b.write_inequivalent_flag = a.write_inequivalent_flag;
  b.info = a.info;
}

// copy
xstructure::xstructure(const xstructure& b) {
  //  free();
  // *this=b;
  copy(b);
}

// destructor
xstructure::~xstructure() {
  free(); //DX20191220 - added free and moved contents below into free
}

// copies xtructures: b=a
const xstructure& xstructure::operator=(const xstructure& b) {  // operator=
  if(this!=&b) {
    free();
    copy(b);
  }
  return *this;
}

void xstructure::clear() { //DX20191220 - uppercase to lowercase clear
  xstructure _tmp;
  (*this)=_tmp;
}

void xstructure::clean() { //DX20191220 - uppercase to lowercase clean
  stringstream ss_xstr;
  ss_xstr << (*this);
  (*this).clear(); //DX20191220 - uppercase to lowercase clear
  ss_xstr >> (*this);
  ss_xstr.str("");
}

void xstructure::ClearSpecies() { //CO20180420 - helps with pocc, match with AddAtom()
  num_each_type.clear();
  comp_each_type.clear();
  stoich_each_type.clear();
  species.clear();species_pp.clear();species_pp_type.clear();species_pp_version.clear();
  species_pp_ZVAL.clear();
  species_pp_vLDAU.clear();
  species_volume.clear();
  species_mass.clear();
}

//ME20211004 - from POCC
void xstructure::CleanStructure() {
  neg_scale = false;  //NO negative scale... doesn't really matter, scale is one variable
  ReScale(1.0);
  ShiftOriginToAtom(0);
  BringInCell();
  clean(); //DX20191220 - uppercase to lowercase clean
}

void xstructure::initialize(const string& structure_title) { //CO20211122 - initialize structure; avoid copying of xstructure
  free(); //DX20191220 - moved contents below into free()
  title=structure_title;
}

void xstructure::initialize(istream& _input,int _iomode) { //DX20210129 - initialize structure; avoid copying of xstructure
  free(); //DX20191220 - added free to initialize
  (*this).iomode=_iomode;
  _input >> (*this);
}

void xstructure::initialize(ifstream& _input,int _iomode) { //DX20210129 - initialize structure; avoid copying of xstructure
  free(); //DX20191220 - added free to initialize
  (*this).iomode=_iomode;
  _input >> (*this);
}

void xstructure::initialize(const stringstream& __input,int _iomode) { //DX20210129 - initialize structure; avoid copying of xstructure
  free(); //DX20191220 - added free to initialize
  (*this).iomode=_iomode;
  stringstream _input(__input.str());
  _input >> (*this);
}

void xstructure::initialize(const string& _input,int _iomode) { //CO20211122 - initialize structure; avoid copying of xstructure
  free(); //DX20191220 - added free to initialize
  stringstream strstream;
  aurostd::efile2stringstream(_input,strstream); //CO20171025
  (*this).iomode=_iomode;
  (*this).directory = _input; //DX20180526 - location of xstructure
  strstream >> (*this);
}

void xstructure::initialize(const string& url,const string& file,int _iomode) { //CO20211122 - initialize structure; avoid copying of xstructure
  free(); //DX20191220 - added free to initialize
  stringstream strstream;
  string _url = url;
  if(url.find(":AFLOW")!=string::npos){ //HE20220615 safeguard against the direct use of AURLs as suggested by CO
    _url = aurostd::StringSubst(url,":AFLOW","/AFLOW");
  }
  aurostd::url2stringstream(_url+"/"+file,strstream);
  (*this).iomode=_iomode;
  (*this).directory = _url+"/"+file; //DX20180526 - location of xstructure
  strstream >> (*this);
}

// **************************************************************************
// Xstructure operator<< OUTPUT_XSTRUCTURE_OUTPUT 
// **************************************************************************
// print an xstructure in a variety of forms
ostream& operator<<(ostream& oss,const xstructure& a) { // operator<<
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream message;
  int a_iomode=a.iomode;
  //  DEBUG=TRUE;
  // ----------------------------------------------------------------------
  // PUT DEFAULT
  if(a_iomode==IOAFLOW_AUTO) a_iomode=IOVASP_AUTO; // some default
  // ----------------------------------------------------------------------
  // VASP OUTPUT
  if(a_iomode==IOVASP_AUTO || a_iomode==IOVASP_POSCAR || a_iomode==IOVASP_ABCCAR || a_iomode==IOVASP_WYCKCAR) { // VASP POSCAR
    oss.setf(std::ios::fixed,std::ios::floatfield);
    uint _precision_=_DOUBLE_WRITE_PRECISION_MAX_; //14; //was 16 SC 10 DM //CO20180515
    oss.precision(_precision_);
    //DX20180618 - Check for symbolic representaion only - START
    if(a.symbolic_math_representation_only){
      xstructure aa(a);
      oss << aa.PrintSymbolicMathRepresentation();
      return oss;
    }
    //DX20180618 - Check for symbolic representation only - END
    if(a_iomode==IOAFLOW_AUTO)   oss << a.title <<endl; // << " (AUTO) " << endl;
    if(a_iomode==IOVASP_AUTO)    oss << a.title <<endl; // << " (AUTO) " << endl;
    if(a_iomode==IOVASP_POSCAR)  oss << a.title <<endl; // << " (POSCAR) " << endl;
    if(a_iomode==IOVASP_ABCCAR)  oss << a.title <<endl; // << " (ABCCAR) " << endl;
    //DX20221128 - START
    if(a_iomode==IOVASP_WYCKCAR) {
      oss << a.title << "| SG: ";
      // check if space group has been calculated, otherwise GetSpaceGroupName() could fail //DX20221128
      try {
        oss << GetSpaceGroupName(a.space_group_ITC,a.directory) << " " << a.space_group_ITC << " PG: " << a.point_group_ITC << " BL: " << a.bravais_label_ITC << " | sym_eps: " << a.sym_eps << endl; //DX20210526 - extend title
      }
      catch(aurostd::xerror& re){ 
        try {
          oss << GetSpaceGroupName(a.spacegroupnumber,a.directory) << " " << a.spacegroupnumber << " PG: " << a.point_group_ITC << " BL: " << a.bravais_label_ITC << " | sym_eps: " << a.sym_eps << endl; //DX20210526 - extend title
        }
        catch(aurostd::xerror& re){ 
          if(LDEBUG){ message << "Cannot determine space group name" << endl; }
        }
      }
    }
    //DX20221128 - STOP
    if(a.neg_scale==FALSE || a_iomode==IOVASP_WYCKCAR) { //DX20210708 - wyccar should always use scale factor
      oss.precision(6);  //DM
      oss << a.scale; // << endl; //CO20170630
    } else {
      // oss << a.scale << endl;
      double s=a.scale;
      double vol=s*s*s*GetVol(a.lattice);
      oss.precision(6);  //DM
      oss << -1*vol;  // << endl; //CO20170630
    }
    oss.precision(_precision_);  //SC to cut/paste from matlab in format long
    //CO20170630, add pocc tol to make truly pocc readable by aflow
    if(a.partial_occupation_flag==TRUE) {
      oss.unsetf(ios_base::floatfield);
      oss << " "; //<< std::defaultfloat;
      if(a.neg_scale_second){oss << (-1)*a.partial_occupation_HNF;}
      else{oss << a.partial_occupation_site_tol;}
      if(1||a.scale_third.isentry){  //CO20170803 - stoich tol //always print
        oss << " "; 
        oss << a.partial_occupation_stoich_tol;
      }
      oss << std::fixed;
    }
    oss << endl;
    // ----------------------------------------------------------------------
    if(a_iomode==IOVASP_POSCAR || a_iomode==IOVASP_AUTO) {
      for(uint i=1;i<=3;i++) {
        for(uint j=1;j<=3;j++) {
          oss << " ";
          if(abs(a.lattice(i,j))<10.0) oss << " ";
          if(!std::signbit(a.lattice(i,j))) oss << " ";
          oss << a.lattice(i,j) << "";
        }
        oss << endl;
      }
    }
    // ----------------------------------------------------------------------
    if(a_iomode==IOVASP_ABCCAR){ //DX20210525 - separated a_iomode==IOVASP_WYCKCAR
      oss << " ";
      oss.precision(10);  //SC to cut/paste from matlab in format long
      if(abs(a.a)<10.0) oss << " ";
      if(!std::signbit(a.a)) oss << " ";
      oss << a.a << "";
      if(abs(a.b)<10.0) oss << " ";
      if(!std::signbit(a.b)) oss << " ";
      oss << a.b << "";
      if(abs(a.c)<10.0) oss << " ";
      if(!std::signbit(a.c)) oss << " ";
      oss << a.c << "";
      oss.precision(4);  //SC to cut/paste from matlab in format long
      if(abs(a.alpha)<10.0) oss << " ";
      if(!std::signbit(a.alpha)) oss << " ";
      oss << a.alpha << "";
      if(abs(a.beta)<10.0)  oss << " ";
      if(!std::signbit(a.beta))  oss << " ";
      oss << a.beta << "";
      if(abs(a.gamma)<10.0) oss << " ";
      if(!std::signbit(a.gamma)) oss << " ";
      oss << a.gamma << "";
      //DX20210621 [OBSOLETE] if(a_iomode==IOVASP_WYCKCAR) oss << "  " << a.spacegroupnumber << "";
      oss << endl;
      oss.precision(_precision_);  //SC to cut/paste from matlab in format long
    }
    else if(a_iomode==IOVASP_WYCKCAR){ //DX20210525 - note wyccar uses lattice parameters of the conventional cell
      xvector<double> data = Getabc_angles(a.standard_lattice_ITC,DEGREES);
      oss << " ";
      oss.precision(10);  //SC to cut/paste from matlab in format long
      if(abs(data(1))<10.0) oss << " ";
      if(!std::signbit(data(1))) oss << " ";
      oss << data(1) << "";
      if(abs(data(2))<10.0) oss << " ";
      if(!std::signbit(data(2))) oss << " ";
      oss << data(2) << "";
      if(abs(data(3))<10.0) oss << " ";
      if(!std::signbit(data(3))) oss << " ";
      oss << data(3) << "";
      oss.precision(4);  //SC to cut/paste from matlab in format long
      if(abs(data(4))<10.0) oss << " ";
      if(!std::signbit(data(4))) oss << " ";
      oss << data(4) << "";
      if(abs(data(5))<10.0)  oss << " ";
      if(!std::signbit(data(5)))  oss << " ";
      oss << data(5) << "";
      if(abs(a.gamma)<10.0) oss << " ";
      if(!std::signbit(data(6))) oss << " ";
      oss << data(6) << "";
      oss << " " << a.space_group_ITC << " " << a.setting_ITC;
      oss << endl;
      oss.precision(_precision_);  //SC to cut/paste from matlab in format long
    }
    // ----------------------------------------------------------------------
    if(a.is_vasp4_poscar_format==TRUE) {
    } // nothing to do

    if(a.is_vasp5_poscar_format==TRUE) {
      for(uint i=0;i<a.species_pp.size();i++) { //ME20190308 - species is empty when structure is based on vasp4 POSCAR
        //oss << a.species.at(i)<< " ";  OBSOLETE ME20190304 - species can contain PP information, which VASP doesn't recognize
        oss << KBIN::VASP_PseudoPotential_CleanName(a.species_pp.at(i)) << " ";  //ME20190308
      }
      oss << endl;
    }
    //DX20210526 - add WYCCAR format - START
    if(a_iomode==IOVASP_WYCKCAR) {
      oss << aurostd::joinWDelimiter(SYM::countWyckoffTypes(a.wyckoff_sites_ITC), " ") << endl; // gets the number of Wyckoff positions per type
      oss << "Direct(WYCCAR)" << endl; // wyccar is always in direct/fractional
      uint nWyckoff_sites = a.wyckoff_sites_ITC.size();
      double _coord = AUROSTD_MAX_DOUBLE;
      for(uint i=0;i<nWyckoff_sites;i++){
        _coord = AUROSTD_MAX_DOUBLE;
        oss << " ";
        for(uint j=1;j<=3;j++) {
          _coord=aurostd::roundoff(a.wyckoff_sites_ITC[i].coord(j),pow(10.0,-(double)_precision_));
          if(abs(_coord)<10.0) oss << " ";
          if(!std::signbit(_coord)) oss << " ";
          oss << _coord << " ";
        }
        oss << std::setprecision(_precision_) << std::left
          << " " << a.wyckoff_sites_ITC[i].type
          << " " << a.wyckoff_sites_ITC[i].multiplicity
          << " " << a.wyckoff_sites_ITC[i].letter
          << " " << a.wyckoff_sites_ITC[i].site_symmetry;
        oss << endl;
      }
    }
    //DX20210526 - add WYCCAR format - END
    // ----------------------------------------------------------------------
    //CO20170630 - fixing for POCC
    //[CO20180705 - we have const str&, so we can't modify atom arrangement, this MUST be done before structure is printed]a.MakeTypes();  //CO20180705 - repetita iuvant
    //[CO20180705 - we have const str&, so we can't modify atom arrangement, this MUST be done before structure is printed]std::stable_sort(a.atoms.begin(),a.atoms.end(),sortAtomsType);  //CO20180705 - this makes it necessary that atoms are properly typed
    if(a_iomode!=IOVASP_WYCKCAR) { //DX20210611 - do not do for Wyccar, atom count is not the same as number of Wyckoff positiosn
      if(a.partial_occupation_flag==TRUE) {
        //need to figure out the '+'
        uint iatom=0;
        vector<vector<uint> > vsame_pocc;
        double last_pocc=0.0;
        if(LDEBUG) {
          for(uint i=0;i<a.atoms.size();i++){
            cerr << __AFLOW_FUNC__ << " name=" << a.atoms[i].name << " type=" << a.atoms[i].type << " pocc=" << a.atoms[i].partial_occupation_value << endl;
          }
        }
        if(a.atoms.size()){
          for(uint i=0;i<a.num_each_type.size();i++){
            vsame_pocc.push_back(vector<uint>(0));  //for first atom
            vsame_pocc.back().push_back(0); //first atom
            last_pocc=a.atoms[iatom].partial_occupation_value;
            for(uint j=0;j<(uint)a.num_each_type[i];j++){
              if(aurostd::isequal(a.atoms[iatom].partial_occupation_value,last_pocc,_AFLOW_POCC_ZERO_TOL_)){vsame_pocc.back().back()+=1;} //same
              else { //'+'
                vsame_pocc.back().push_back(1);
                last_pocc=a.atoms[iatom].partial_occupation_value;
              }
              iatom++;
            }
          }
        }
        if(LDEBUG) {
          for(uint i=0;i<vsame_pocc.size();i++){
            for(uint j=0;j<vsame_pocc[i].size();j++){
              cerr << __AFLOW_FUNC__ << " vsame_pocc[" << i << "][" << j << "]=" << vsame_pocc[i][j] << endl;
            }
          }
        }
        iatom=0;
        for(uint i=0;i<vsame_pocc.size();i++){
          //if(vsame_pocc[i].size()==1){  //no '+'
          //  oss << vsame_pocc[i][0] << "*";
          //  //oss << std::defaultfloat;
          //  oss.unsetf(ios_base::floatfield);
          //  oss << a.atoms[iatom++].partial_occupation_value << std::fixed << " ";
          //} else {  //need '+'
          for(uint j=0;j<vsame_pocc[i].size();j++){
            oss << vsame_pocc[i][j] << "*";
            //oss << std::defaultfloat;
            oss.unsetf(ios_base::floatfield);
            oss << a.atoms[iatom].partial_occupation_value << std::fixed << (j!=vsame_pocc[i].size()-1?"+":" ");
            iatom+=vsame_pocc[i][j];
          }
          //}
        }

        //[OBSOLETE - CO20180705]for(uint i=0;i<a.num_each_type.size();i++){
        //[OBSOLETE - CO20180705]  oss << a.num_each_type.at(i) << "*";
        //[OBSOLETE - CO20180705]  //oss << std::defaultfloat;
        //[OBSOLETE - CO20180705]  oss.unsetf(ios_base::floatfield);
        //[OBSOLETE - CO20180705]  oss << a.atoms[i].partial_occupation_value << std::fixed << " ";
        //[OBSOLETE - CO20180705]}
      } else {for(uint i=0;i<a.num_each_type.size();i++){oss << a.num_each_type.at(i) << " ";}}
      oss << endl;
      if(a.isd) oss << "Selective Dynamics" << endl; // DONE YOYO BUG
      // oss << a.coord_type << endl;

      if(a.coord_flag==_COORDS_FRACTIONAL_) oss << "Direct(" << a.atoms.size() << ") ";
      if(a.coord_flag==_COORDS_CARTESIAN_)  oss << "Cartesian(" << a.atoms.size() << ") ";
      //  if(a.partial_occupation_flag==TRUE)  oss << "Pocc ";
      if(a.order_parameter_structure==TRUE)  oss << "OrderParameter(" << a.order_parameter_atoms.size() << ") ";
      if(1) { // write [A1B2C3D4]
        //      oss << "[";for(uint i=0;i<a.num_each_type.size();i++) {oss << char('A'+i) << a.num_each_type.at(i);}oss << "] ";
        //CO20170630, the original num_each_type doesn't work here, so we fix
        if(a.partial_occupation_flag==TRUE) {
          oss << "Partial ";
          //oss.precision(_pocc_precision_);  //CO20170630
          //oss << std::defaultfloat;
          int comp_prec=(int)ceil(log10(1.0/a.partial_occupation_stoich_tol));  //ceil ensures we round up above 1 //CO20181226
          oss.precision(comp_prec); //CO20181226
          oss.unsetf(ios_base::floatfield);
          oss << "[";for(uint i=0;i<a.comp_each_type.size();i++) {oss << char('A'+i) << a.comp_each_type.at(i);}oss << "] ";
          oss << std::fixed;
          oss.precision(_precision_);       //CO20170630 //CO20181226
        } else {
          oss << "["; for(uint i=0,k=0;i<a.num_each_type.size();k+=a.num_each_type.at(i),i++) { oss << char(a.atoms.at(k).type+65) << a.num_each_type.at(i);} oss << "] ";
        }
      }
      // done
      oss << endl;

      double _coord;  //CO20190322 - remove annoying -0.0000000
      for(uint iat=0;iat<a.atoms.size();iat++) {
        oss << " ";
        for(uint j=1;j<=3;j++) {
          //	oss << " ";
          if(a.coord_flag==_COORDS_FRACTIONAL_) {_coord=a.atoms.at(iat).fpos(j);} //CO20190322 - remove annoying -0.0000000
          if(a.coord_flag==_COORDS_CARTESIAN_)  {_coord=a.atoms.at(iat).cpos(j);} //CO20190322 - remove annoying -0.0000000

          _coord=aurostd::roundoff(_coord,pow(10.0,-(double)_precision_)); //CO20190322 - remove annoying -0.0000000
          if(abs(_coord)<10.0) oss << " "; //CO20190322 - remove annoying -0.0000000
          if(!std::signbit(_coord)) oss << " ";  //CO20190322 - remove annoying -0.0000000
          oss << _coord << " "; //CO20190322 - remove annoying -0.0000000

          //[CO20190322 OBSOLETE]if(a.coord_flag==_COORDS_FRACTIONAL_) {if(abs(a.atoms.at(iat).fpos(j))<10.0) oss << " ";if(!std::signbit(a.atoms.at(iat).fpos(j))) oss << " "; oss << a.atoms.at(iat).fpos(j) << " ";}
          //[CO20190322 OBSOLETE]if(a.coord_flag==_COORDS_CARTESIAN_)  {if(abs(a.atoms.at(iat).cpos(j))<10.0) oss << " ";if(!std::signbit(a.atoms.at(iat).cpos(j))) oss << " "; oss << a.atoms.at(iat).cpos(j) << " ";}
        }
        //  cout << aurostd::modulus(a.atoms.at(iat).cpos) << " ";
        if(a.isd==TRUE)
          oss << " " << a.atoms.at(iat).sd[0] << " " << a.atoms.at(iat).sd[1] << " " << a.atoms.at(iat).sd[2];
        if(a.atoms.at(iat).name_is_given==TRUE) {
          oss << " " << a.atoms.at(iat).name << " ";
          for(uint j=a.atoms.at(iat).name.length();j<5;j++) oss << " ";
        }
        if(a.partial_occupation_flag==TRUE) {
          //oss.precision(_pocc_precision_);  //CO20170630
          //if(a.atoms.at(iat).partial_occupation_flag==FALSE) oss << "-      ";
          //	if(a.atoms.at(iat).partial_occupation_flag==TRUE) oss << a.atoms.at(iat).partial_occupation_value << "  ";// << " (" << iat << "/" << a.partial_occupation_flags.size() << ")";
          //oss << std::defaultfloat;
          oss.unsetf(ios_base::floatfield);
          oss << "pocc=" << a.atoms.at(iat).partial_occupation_value << "  ";
          oss << std::fixed;
          //oss.precision(_precision_); //CO20170630
        }
        if(a.order_parameter_structure==TRUE) {
          if(a.atoms.at(iat).order_parameter_atom==FALSE) oss << "- ";
          if(a.atoms.at(iat).order_parameter_atom==TRUE) oss << a.atoms.at(iat).order_parameter_value << " ";// << " (" << iat << "/" << a.order_parameter_atoms.size() << ")";
        }
        if(a.write_inequivalent_flag==TRUE) {
          oss << " ";
          // ?	if(i<10) oss << "0";
          oss << iat << "[";
          if(a.atoms.at(iat).equivalent<10) oss << "0";
          oss << a.atoms.at(iat).equivalent << "]";
          if(a.atoms.at(iat).is_inequivalent) {
            oss <<"*";
            oss << "_(" << a.atoms.at(iat).num_equivalents << ") "; //<< "  index=" << a.atoms.at(iat).index_iatoms << "  ";
            //  " v" << a.iatoms.size() << "   burp ";
            // for(uint jat=0;jat<a.iatoms.size();jat++)  oss << a.iatoms.at(jat).size() << " ";
          }
        }
        if(a.qm_forces_write) {
          if(a.qm_calculated==TRUE)  oss << "F *(";
          if(a.qm_calculated==FALSE) oss << "F  (";
          for(uint j=1;j<=3;j++) {
            if(abs(a.qm_forces.at(iat)(j))<10.0) oss << " ";
            if(!std::signbit(a.qm_forces.at(iat)(j))) oss << " ";
            oss << a.qm_forces.at(iat)(j) << " ";
          }
          oss << ")_   ";
        }
        if(a.qm_positions_write) {
          if(a.qm_calculated==TRUE)  oss << "P *(";
          if(a.qm_calculated==FALSE) oss << "P  (";
          for(uint j=1;j<=3;j++) {
            if(abs(a.qm_positions.at(iat)(j))<10.0) oss << " ";
            if(!std::signbit(a.qm_positions.at(iat)(j))) oss << " ";
            oss << a.qm_positions.at(iat)(j) << " ";
          }
          oss << ")_   ";
        }
        if(a.write_DEBUG_flag) {
          oss << " s"<<a.atoms.at(iat).spin;
          //[CO20200130 - number->basis]oss << " n"<<a.atoms.at(iat).number;
          oss << " b"<<a.atoms.at(iat).basis;
          oss << " N("<<a.atoms.at(iat).cleanname;
          oss << " "<<a.atoms.at(iat).atomic_number<<" "<<" ["<<a.atoms.at(iat).type<<"] ";
          oss << " ijk("<<a.atoms.at(iat).ijk(1)<<","<<a.atoms.at(iat).ijk(2)<<","<<a.atoms.at(iat).ijk(3)<<")";
        }
        oss << endl;oss.flush();
      } // iat
    } //DX20210610 - end Wyccar if-statement
    if(a.write_lattice_flag) {
      oss << "DIRECT LATTICE per raw" << endl;
      for(uint i=1;i<=3;i++) {
        for(uint j=1;j<=3;j++) {
          oss << " ";
          if(!std::signbit(a.scale*a.lattice(i,j))) oss << " ";
          oss << a.scale*a.lattice(i,j) << " ";
        }
        oss << endl;
      }
    }
    if(a.write_klattice_flag) {
      oss << "RECIPROCAL LATTICE per raw" << endl;
      for(uint i=1;i<=3;i++) {
        for(uint j=1;j<=3;j++) {
          oss << " ";
          if(!std::signbit(a.klattice(i,j))) oss << " ";
          oss << a.klattice(i,j) << " ";
        }
        oss << endl;
      }
    }
    if(a.write_lattice_flag && a.write_klattice_flag) {
      oss << "ORTOGONALITY (a*b')/2pi=I" << endl;
      xmatrix<double> orto(3,3);
      orto=((a.scale*a.lattice))*trasp(a.klattice)/(2.0*pi);
      for(uint i=1;i<=3;i++) {
        for(uint j=1;j<=3;j++) {
          oss << " ";
          if(!std::signbit(orto(i,j))) oss << " ";
          oss << orto(i,j) << " ";
        }
        oss << endl;
      }
    }
    if(a.write_DEBUG_flag) {
      oss << "kpoints_k1,k2,k3 = " << a.kpoints_k1 << "," << a.kpoints_k2 << "," << a.kpoints_k3 << endl;
      oss << "kpoints_kmax     = " << a.kpoints_kmax << endl;
      oss << "kpoints_kppra    = " << a.kpoints_kppra << endl;
      oss << "kpoints_mode     = " << a.kpoints_mode << endl;
      oss << "kpoints_kscheme  = " << a.kpoints_kscheme << endl;
    }
    if(0 && a.partial_occupation_flag==TRUE) { //KY
      oss << "*******************************" << endl;
      for(int hnf=1;hnf<10;hnf++) {
        oss << hnf << " ";
        oss << endl;
      }
    }
    //DX20180618 - Check for symmetry constrained calculation - START
    if(a.constrained_symmetry_calculation){
      xstructure aa(a);
      oss << aa.PrintSymbolicMathRepresentation();
    }
    //DX20180618 - Check for symmetry constrained calculation - END
    return oss;
  } // END OF VASP
  // ----------------------------------------------------------------------
  //  QUANTUM ESPRESSO OUTPUT
  if(a_iomode==IOQE_AUTO || a_iomode==IOQE_GEOM) { // VASP POSCAR
    oss << "! AFLOW::QE BEGIN " << endl;
    uint depthQE=27;
    if(a_iomode==IOQE_AUTO) oss << "! " << a.title <<endl;//<< " (AUTO)" <<endl;
    if(a_iomode==IOQE_GEOM) oss << "! " << a.title <<endl;//<< " (GEOM)" <<endl;
    oss << aurostd::PaddedPOST("&system",depthQE," ") << " ! // aflow " << endl;
    oss << aurostd::PaddedPOST(" ibrav=0,",depthQE," ") << " ! // free " << endl;
    oss << aurostd::PaddedPOST(" nat="+aurostd::utype2string(a.atoms.size())+",",depthQE) << " ! // a.atoms.size() " << endl;
    oss << aurostd::PaddedPOST(" ntyp="+aurostd::utype2string(a.num_each_type.size()),depthQE) << " ! // a.num_each_type.size() " << endl;
    //oss << aurostd::PaddedPOST(" ntyp="+aurostd::utype2string(a.num_each_type.size())+",",depthQE) << " ! // a.num_each_type.size() " << endl;  //CO20171010
    //oss << aurostd::PaddedPOST(" ecutwfc=_AFLOW_ECUTWFC_,",depthQE," ") << " ! // fix these " << endl;  //CO20171010
    //oss << aurostd::PaddedPOST(" ecutrho=_AFLOW_ECUTRHO_",depthQE," ") << " ! // fix these " << endl;   //CO20171010
    oss << " /" << endl;
    oss.setf(std::ios::fixed,std::ios::floatfield);
    uint _precision_=_DOUBLE_WRITE_PRECISION_MAX_; //14; //was 16 SC 10 DM //CO20180515
    oss.precision(_precision_);
    if(a.coord_flag==_COORDS_FRACTIONAL_) oss << "ATOMIC_POSITIONS (crystal)" << endl;
    if(a.coord_flag==_COORDS_CARTESIAN_)  oss << "ATOMIC_POSITIONS (angstrom)" << endl;
    for(uint iat=0;iat<a.atoms.size();iat++) {
      oss << " ";
      if(a.atoms.at(iat).name_is_given==TRUE) {
        oss << " " << aurostd::PaddedPOST(KBIN::VASP_PseudoPotential_CleanName(a.atoms.at(iat).name),5," ") << " ";
      } else {
        message << "QE needs atoms species names";  //CO20190629
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_MISSING_); //CO20190629
      }
      for(uint j=1;j<=3;j++) {
        //  oss << " ";
        if(a.coord_flag==_COORDS_FRACTIONAL_) {if(abs(a.atoms.at(iat).fpos(j))<10.0) oss << " ";if(!std::signbit(a.atoms.at(iat).fpos(j))) oss << " "; oss << a.atoms.at(iat).fpos(j) << " ";}
        if(a.coord_flag==_COORDS_CARTESIAN_)  {if(abs(a.atoms.at(iat).cpos(j))<10.0) oss << " ";if(!std::signbit(a.atoms.at(iat).cpos(j))) oss << " "; oss << a.atoms.at(iat).cpos(j) << " ";}
      }
      oss << " ! // " << a.atoms.at(iat).cleanname << " ";
      if(a.write_inequivalent_flag==TRUE) {
        oss << " ";
        //	if(i<10) oss << "0";
        oss << iat << "[";
        if(a.atoms.at(iat).equivalent<10) oss << "0";
        oss << a.atoms.at(iat).equivalent << "]";
        if(a.atoms.at(iat).is_inequivalent) {
          oss <<"*";
          oss << "_(" << a.atoms.at(iat).num_equivalents << ") "; //<< "  index=" << a.atoms.at(iat).index_iatoms << "  ";
          //  " v" << a.iatoms.size() << "   burp ";
          // for(uint jat=0;jat<a.iatoms.size();jat++)  oss << a.iatoms.at(jat).size() << " ";
        }
      }
      oss << endl;
    }
    // ----------------------------------------------------------------------
    oss.precision(_precision_);  //SC to cut/paste from matlab in format long
    oss << "CELL_PARAMETERS (angstrom)" << endl ;
    {
      for(uint i=1;i<=3;i++) {
        for(uint j=1;j<=3;j++) {
          oss << " ";
          if(abs(a.lattice(i,j))<10.0) oss << " ";
          if(!std::signbit(a.lattice(i,j))) oss << " ";
          oss << a.lattice(i,j)*a.scale << ""; //DX20180215 - added scaling factor
        }
        oss << endl;
      }
    }
    oss << "# AFLOW::QE END " << endl;
    return oss;
  }
  // ----------------------------------------------------------------------
  //  CIF OUTPUT
  if(a.iomode==IOCIF) { // CIF
    pflow::PrintCIF(oss, a, a.spacegroupnumber, a.setting_ITC); //DX20210630 - add setting (otherwise, this will mess up rhl systems by mixing the hex and rhl setting)
    return oss;
  }

  // ----------------------------------------------------------------------
  //  ABINIT OUTPUT
  if(a_iomode==IOABINIT_AUTO || a_iomode==IOABINIT_GEOM) { // VASP POSCAR
    xstructure aa(a); //DX20210415 - need to make a copy to rescale
    aa.ReScale(1.0); //DX20210415
    oss << "# AFLOW::ABINIT BEGIN " << endl;
    uint _precision_=_DOUBLE_WRITE_PRECISION_MAX_; //14; //was 16 SC 10 DM //CO20180515
    oss.precision(_precision_);
    oss.setf(std::ios::fixed,std::ios::floatfield);
    if(a_iomode==IOABINIT_AUTO) oss << "# " << aa.title <<endl;//<< " (AUTO)" << endl;
    if(a_iomode==IOABINIT_GEOM) oss << "# " << aa.title <<endl;//<< " (GEOM)" << endl;
    oss << "acell   " << double(1) << "   " << double(1) << "   " << double(1) << "  ANGSTR" << endl; // scaling of the primitive vectors, in Bohr.
    for(uint j=1;j<=3;j++) { //CO20190908 - manual is misleading, it's row-based// each COLUMN of this array is one primitive translation
      if(j==1) oss << "rprim";
      if(j==2) oss << "     ";
      if(j==3) oss << "     ";
      for(uint i=1;i<=3;i++) {
        oss << " ";
        if(abs(aa.lattice(j,i))<10.0) oss << " ";  //CO20190908 - manual is misleading, it's row-based
        if(!std::signbit(aa.lattice(j,i))) oss << " "; //CO20190908 - manual is misleading, it's row-based
        oss << aa.lattice(j,i) << "";  //CO20190908 - manual is misleading, it's row-based
      }
      oss << endl;
    }
    oss << "natom " << aa.atoms.size() << endl;
    //DX20200313 - add atom type info via znucl - START
    oss << "znucl ";
    for(uint i=0;i<aa.species.size();i++){
      for(uint e=0;e<velement.size();e++){ //external variable (see aflow_xelement.h)
        if(velement[e].symbol == KBIN::VASP_PseudoPotential_CleanName(aa.species[i])){
          oss << e << " "; // index corresponds to Z value 
          break;
        }
      }
    }
    oss << endl;
    //DX20200313 - add atom type info via znucl - END
    oss << "typat ";
    //   for(uint i=0;i<aa.num_each_type.size();i++) oss << a.num_each_type.at(i) << " ";  oss << endl;
    for(uint i=0;i<aa.atoms.size();i++)
      oss << aa.atoms.at(i).type+1 << " ";
    oss << endl;
    if(aa.coord_flag==_COORDS_FRACTIONAL_) oss << "xred " << endl;
    if(aa.coord_flag==_COORDS_CARTESIAN_) oss << "xangst " << endl;
    for(uint iat=0;iat<aa.atoms.size();iat++) {
      oss << "      ";
      for(uint j=1;j<=3;j++) {
        if(aa.coord_flag==_COORDS_FRACTIONAL_) {
          if(abs(aa.atoms.at(iat).fpos(j))<10.0) oss << " ";
          if(!std::signbit(aa.atoms.at(iat).fpos(j))) oss << " ";
          oss << aa.atoms.at(iat).fpos(j) << " ";
        }
        if(aa.coord_flag==_COORDS_CARTESIAN_) {
          if(abs(aa.atoms.at(iat).cpos(j))<10.0) oss << " ";
          if(!std::signbit(aa.atoms.at(iat).cpos(j))) oss << " ";
          oss << aa.atoms.at(iat).cpos(j) << " ";
        }
      }
      oss << " # " << aa.atoms.at(iat).cleanname << " ";
      if(aa.write_inequivalent_flag==TRUE) {
        oss << " ";
        //	if(i<10) oss << "0";
        oss << iat << "[";
        if(aa.atoms.at(iat).equivalent<10) oss << "0";
        oss << aa.atoms.at(iat).equivalent << "]";
        if(aa.atoms.at(iat).is_inequivalent) {
          oss <<"*";
          oss << "_(" << aa.atoms.at(iat).num_equivalents << ") "; //<< "  index=" << a.atoms.at(iat).index_iatoms << "  ";
          //  " v" << aa.iatoms.size() << "   burp ";
          // for(uint jat=0;jat<aa.iatoms.size();jat++)  oss << a.iatoms.at(jat).size() << " ";
        }
      }
      oss << endl;
    }
    oss << "# AFLOW::ABINIT END " << endl;
    return oss;
  } 

  // ----------------------------------------------------------------------
  //  ELK OUTPUT //DX20200315
  if(a_iomode==IOELK_AUTO || a_iomode==IOELK_GEOM) { // ELK
    xstructure aa(a); //DX20210415 - need to make a copy to rescale
    aa.ReScale(1.0); //DX20210415
    oss << "# AFLOW::ELK BEGIN " << endl;
    uint _precision_=_DOUBLE_WRITE_PRECISION_MAX_; //14; //was 16 SC 10 DM //CO20180515
    oss.precision(_precision_);
    oss.setf(std::ios::fixed,std::ios::floatfield);
    oss << "# " << aa.title <<endl;
    oss << endl;
    //scaling factors
    oss << "scale" << endl << " " << aa.scale << endl << endl;
    oss << "scale1" << endl << " 1.0" << endl << endl; // returns unscaled (for now)
    oss << "scale2" << endl << " 1.0" << endl << endl; // returns unscaled (for now) 
    oss << "scale3" << endl << " 1.0" << endl << endl; // returns unscaled (for now)

    //lattice, note: convert to atomic units (Bohr)
    oss << "avec" << endl;
    oss << " " << aa.lattice(1)*angstrom2bohr << endl;
    oss << " " << aa.lattice(2)*angstrom2bohr << endl;
    oss << " " << aa.lattice(3)*angstrom2bohr << endl;
    oss << endl;

    // atom info
    oss << "atoms" << endl;
    oss << " " << setw(49) << std::left << aa.species.size();
    oss << ": nspecies" << endl;
    //DX20210409 [OBSOLETE] xvector<double> magnetic_field; // not currently supported; zero vector for now
    for(uint i=0;i<aa.num_each_type.size();i++){
      oss << setw(50) << std::left << "\'" + aa.species[i] + ".in\'";
      oss << ": spfname" << endl;
      oss << " " << setw(49) << std::left << aa.num_each_type[i];
      oss << ": natoms; atpos, bfcmt below" << endl;
      for(uint iat=0;iat<aa.atoms.size();iat++){
        if(aa.atoms[iat].name == aa.species[i]){
          // atom coordinates
          for(uint j=1;j<=3;j++) {
            if(aa.coord_flag==_COORDS_CARTESIAN_) {
              if(abs(aa.atoms.at(iat).fpos(j))<10.0) oss << " ";
              if(!std::signbit(aa.atoms.at(iat).fpos(j))) oss << " ";
              oss << aa.atoms.at(iat).fpos(j) << " ";
            }
            if(aa.coord_flag==_COORDS_FRACTIONAL_) {
              if(abs(aa.atoms.at(iat).fpos(j))<10.0) oss << " ";
              if(!std::signbit(aa.atoms.at(iat).fpos(j))) oss << " ";
              oss << aa.atoms.at(iat).fpos(j) << " ";
            }
          }
          // magnetic field //DX20210409 - updated with non-collinear spin
          for(uint j=1;j<=3;j++) {
            if(abs(aa.atoms[iat].noncoll_spin(j))<10.0) oss << " ";
            if(!std::signbit(aa.atoms[iat].noncoll_spin(j))) oss << " ";
            oss << aa.atoms[iat].noncoll_spin(j) << " ";
          }
          oss << endl;
        }
      }
    }
    oss << "# AFLOW::ELK END " << endl;
    return oss;
  }

  // ----------------------------------------------------------------------
  //  AIMS OUTPUT
  if(a_iomode==IOAIMS_AUTO || a_iomode==IOAIMS_GEOM) { // VASP POSCAR
    xstructure aa(a);
    aa.ReScale(1.0); //very important because there is NO scale factor in AIMS //CO20180420
    uint _precision_=_DOUBLE_WRITE_PRECISION_MAX_; //14; //was 16 SC 10 DM //CO20180515
    oss.precision(_precision_);
    oss.setf(std::ios::fixed,std::ios::floatfield);
    oss << "# " << aa.title <<endl;//<< " (AUTO)" << endl;
    oss << "# AFLOW::AIMS BEGIN " << endl;  //come after title
    //DX20180618 - Check for symbolic representaion only - START
    if(aa.symbolic_math_representation_only){
      oss << aa.PrintSymbolicMathRepresentation();
      oss << "# AFLOW::AIMS END " << endl;
      return oss;
    }
    //DX20180618 - Check for symbolic representation only - END
    for(uint i=1;i<=3;i++) { // each COLUMN of this array is one primitive translation
      oss << "lattice_vector ";
      for(uint j=1;j<=3;j++) {
        if(abs(aa.lattice(i,j))<100.0) oss << " ";
        if(abs(aa.lattice(i,j))<10.0) oss << " ";
        if(!std::signbit(aa.lattice(i,j))) oss << " ";
        oss << aa.lattice(i,j) << "";
      }
      oss << endl;
    }
    for(uint iat=0;iat<aa.atoms.size();iat++) {
      if(0||LDEBUG){
        cerr << "XSTRUCTURE<<: AIMS aa.coord_flag==" << aa.coord_flag << endl;
        cerr << "XSTRUCTURE<<: AIMS _COORDS_FRACTIONAL_==" << _COORDS_FRACTIONAL_ << endl;
        cerr << "XSTRUCTURE<<: AIMS _COORDS_CARTESIAN_==" << _COORDS_CARTESIAN_ << endl;
      }
      oss << (aa.coord_flag==_COORDS_FRACTIONAL_ ? "atom_frac" : "atom") << " ";
      for(uint j=1;j<=3;j++) {
        if(aa.coord_flag==_COORDS_FRACTIONAL_) {
          if(abs(aa.atoms.at(iat).fpos(j))<100.0) oss << " ";
          if(abs(aa.atoms.at(iat).fpos(j))<10.0) oss << " ";
          if(!std::signbit(aa.atoms.at(iat).fpos(j))) oss << " ";
          oss << aa.atoms.at(iat).fpos(j) << " ";
        }
        if(aa.coord_flag==_COORDS_CARTESIAN_)  {
          if(abs(aa.atoms.at(iat).cpos(j))<100.0) oss << " ";
          if(abs(aa.atoms.at(iat).cpos(j))<10.0) oss << " ";
          if(!std::signbit(aa.atoms.at(iat).cpos(j))) oss << " ";
          oss << aa.atoms.at(iat).cpos(j) << " ";
        }
      }
      oss << " " << aa.atoms.at(iat).cleanname << " ";
      if(aa.write_inequivalent_flag==TRUE) {
        oss << " # ";
        //	if(i<10) oss << "0";
        oss << iat << "[";
        if(aa.atoms.at(iat).equivalent<10) oss << "0";
        oss << aa.atoms.at(iat).equivalent << "]";
        if(aa.atoms.at(iat).is_inequivalent) {
          oss <<"*";
          oss << "_(" << aa.atoms.at(iat).num_equivalents << ") "; //<< "  index=" << aa.atoms.at(iat).index_iatoms << "  ";
          //  " v" << aa.iatoms.size() << "   burp ";
          // for(uint jat=0;jat<aa.iatoms.size();jat++)  oss << aa.iatoms.at(jat).size() << " ";
        }
      }
      oss << endl;
    }
    //DX20180618 - Check for symmetry constrained calculation - START
    if(aa.constrained_symmetry_calculation){
      oss << aa.PrintSymbolicMathRepresentation();
    }
    //DX20180618 - Check for symmetry constrained calculation - END
    oss << "# AFLOW::AIMS END " << endl;
    return oss;
  } 

  // ----------------------------------------------------------------------
  //  ATAT OUTPUT // SD20220123
  //  Alloy-Theoretic Automated Toolkit
  //  See: https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/manual.pdf
  if(a_iomode == IOATAT_STR) { // ATAT
    xstructure aa(a);
    uint _precision_=_DOUBLE_WRITE_PRECISION_MAX_; //14; //was 16 SC 10 DM //CO20180515
    oss.precision(_precision_);
    oss.setf(std::ios::fixed,std::ios::floatfield);
    xmatrix<double> axes = aurostd::eye<double>(3, 3); // set axes to idenity
    // write the axes
    for (uint i = 1; i <= 3; i++) {
      for (uint j = 1; j <= 3; j++) {
        oss << axes(i, j) << " ";
      }
      oss << endl;
    }
    // write the fractional cell vectors (== lattice)
    for (uint i = 1; i <= 3; i++) {
      for (uint j = 1; j <= 3; j++) {
        oss << " ";
        if (abs(aa.lattice(i, j)) < 10.0) {oss << " ";}
        if (!std::signbit(aa.lattice(i, j))) {oss << " ";}
        oss << aa.lattice(i, j) << "";
      }
      oss << endl;
    }
    // write the atoms
    xvector<double> coord(3);
    for (uint iat = 0; iat < aa.atoms.size(); iat++) {
      oss << " ";
      for (uint i=1; i <= 3; i++) {
        if (abs(aa.atoms[iat].cpos(i)) < 10.0) {oss << " ";}
        if (!std::signbit(aa.atoms[iat].cpos(i))) {oss << " ";}
        oss << aa.atoms[iat].cpos(i) << " ";
      }
      if (aa.atoms[iat].name_is_given == TRUE) {oss << aa.atoms[iat].cleanname;}
      oss << endl;
    }
    return oss;
  }
  // ----------------------------------------------------------------------

  oss << "NOT CODED YET" << endl;
  return oss;
}

// **************************************************************************
// PrintSymbolicMathRepresentation
// **************************************************************************
string xstructure::PrintSymbolicMathRepresentation(void){
  xstructure aa((*this));
  int a_iomode=aa.iomode;
  stringstream oss;
  //VASP
  if(a_iomode==IOAIMS_AUTO || a_iomode==IOVASP_AUTO || a_iomode==IOVASP_POSCAR) { // VASP POSCAR
    stringstream title;
    title << aa.title << " # parameters: " << aa.num_parameters << " # lattice parameters: " << aa.num_lattice_parameters << " # Wyckoff parameters: " << (aa.num_parameters-aa.num_lattice_parameters);
    oss << title.str() << endl; 
    oss << "1.0" << endl; //scaling factor
    for(uint i=0;i<aa.symbolic_math_lattice.size();i++){
      oss << "   " << aurostd::joinWDelimiter(aa.symbolic_math_lattice[i],"  ") << endl;
    }
    if(aa.is_vasp4_poscar_format==TRUE) {
    } // nothing to do

    if(aa.is_vasp5_poscar_format==TRUE) {
      for(uint i=0;i<aa.species.size();i++)
        oss << aa.species.at(i)<< " ";
      oss << endl;
    }
    if(aa.partial_occupation_flag==TRUE) {
      for(uint i=0;i<aa.num_each_type.size();i++){
        oss << aa.num_each_type.at(i) << "*";
        oss.unsetf(ios_base::floatfield);
        oss << aa.atoms[i].partial_occupation_value << std::fixed << " ";
      }
    } else {
      for(uint i=0;i<aa.num_each_type.size();i++){
        oss << aa.num_each_type.at(i) << " ";
      }
    }
    oss << endl;
    if(aa.coord_flag==_COORDS_FRACTIONAL_) {oss << "Direct(" << aa.atoms.size() << ") ";}
    if(aa.coord_flag==_COORDS_CARTESIAN_) {oss << "Cartesian(" << aa.atoms.size() << ") ";}
    if(aa.order_parameter_structure==TRUE) {oss << "OrderParameter(" << aa.order_parameter_atoms.size() << ") ";}
    if(1) { 
      if(aa.partial_occupation_flag==TRUE) {
        oss << "Partial ";
        //oss.precision(_pocc_precision_);  //CO20170630
        //oss << std::defaultfloat;
        oss.unsetf(ios_base::floatfield);
        oss << "[";for(uint i=0;i<aa.comp_each_type.size();i++) {oss << char('A'+i) << aa.comp_each_type.at(i);}oss << "] ";
        oss << std::fixed;
        //oss.precision(_precision_);       //CO20170630
      } else {
        oss << "["; for(uint i=0,k=0;i<aa.num_each_type.size();k+=aa.num_each_type.at(i),i++) { oss << char(aa.atoms.at(k).type+65) << aa.num_each_type.at(i);} oss << "] ";
      }
    }
    // done
    oss << endl;
    for(uint i=0;i<aa.atoms.size();i++){
      if(aa.coord_flag==_COORDS_FRACTIONAL_) {
        oss << "   " << aurostd::joinWDelimiter(aa.atoms[i].fpos_equation,"  ");
      }
      else if(aa.coord_flag==_COORDS_CARTESIAN_) {
        oss << "   " << aurostd::joinWDelimiter(aa.atoms[i].cpos_equation,"  ");
      }
      if(aa.atoms.at(i).name_is_given==TRUE) {
        oss << " " << aa.atoms.at(i).name << " ";
        for(uint j=aa.atoms.at(i).name.length();j<5;j++) oss << " ";
      }
      if(aa.partial_occupation_flag==TRUE) {
        oss.unsetf(ios_base::floatfield);
        oss << "pocc=" << aa.atoms.at(i).partial_occupation_value << "  ";
        oss << std::fixed;
      }
      oss << endl;
    }
  }
  //QE

  //ABINIT

  //AIMS
  if(a_iomode==IOAIMS_AUTO || a_iomode==IOAIMS_GEOM) { // AIMS GEOM 
    oss << "# format: symmetry_n_params [n n_lv n_fracpos]" << endl;
    oss << "symmetry_n_params " << aa.num_parameters << " " << aa.num_lattice_parameters << " " << (aa.num_parameters-aa.num_lattice_parameters) << endl; 
    // change lattice parameter ratio to separate parameters
    vector<string> parameter_list;
    for(uint i=0;i<aa.prototype_parameter_list.size();i++){
      if(aa.prototype_parameter_list[i] == "b/a"){
        parameter_list.push_back("b");
      }
      else if(aa.prototype_parameter_list[i] == "c/a"){
        parameter_list.push_back("c");
      }
      else {
        parameter_list.push_back(aa.prototype_parameter_list[i]);
      }
    }
    oss << "symmetry_params " << aurostd::joinWDelimiter(parameter_list," ") << endl;
    for(uint i=0;i<3;i++){
      oss << "symmetry_lv " << aurostd::joinWDelimiter(aa.symbolic_math_lattice[i]," , ") << endl;
    }
    if(aa.coord_flag==_COORDS_FRACTIONAL_) {
      for(uint i=0;i<aa.atoms.size();i++){
        oss << "symmetry_frac " << aurostd::joinWDelimiter(aa.atoms[i].fpos_equation," , ") << endl;
      }
    }
    if(aa.coord_flag==_COORDS_CARTESIAN_) {
      for(uint i=0;i<aa.atoms.size();i++){
        oss << "symmetry_cart " << aurostd::joinWDelimiter(aa.atoms[i].cpos_equation," , ") << endl;
      }
    }

    //for(uint i=1;i<=3;i++) { // each COLUMN of this array is one primitive translation
    //oss << "lattice_vector ";
    //for(uint j=1;j<=3;j++) {
    //if(abs(aa.lattice(i,j))<100.0) oss << " ";
    //if(abs(aa.lattice(i,j))<10.0) oss << " ";
    //if(!std::signbit(aa.lattice(i,j))) oss << " ";
    //oss << aa.lattice(i,j) << "";
    //}
    //oss << endl;
    //}
    //for(uint iat=0;iat<aa.atoms.size();iat++) {
    //if(0||LDEBUG){
    //cerr << "XSTRUCTURE<<: AIMS aa.coord_flag==" << aa.coord_flag << endl;
    //cerr << "XSTRUCTURE<<: AIMS _COORDS_FRACTIONAL_==" << _COORDS_FRACTIONAL_ << endl;
    //cerr << "XSTRUCTURE<<: AIMS _COORDS_CARTESIAN_==" << _COORDS_CARTESIAN_ << endl;
    //}
    //oss << (aa.coord_flag==_COORDS_FRACTIONAL_ ? "atom_frac" : "atom") << " ";
    //for(uint j=1;j<=3;j++) {
    //if(aa.coord_flag==_COORDS_FRACTIONAL_) {
    //if(abs(aa.atoms.at(iat).fpos(j))<100.0) oss << " ";
    //if(abs(aa.atoms.at(iat).fpos(j))<10.0) oss << " ";
    //if(!std::signbit(aa.atoms.at(iat).fpos(j))) oss << " ";
    //oss << aa.atoms.at(iat).fpos(j) << " ";
    //}
    //if(aa.coord_flag==_COORDS_CARTESIAN_)  {
    //if(abs(aa.atoms.at(iat).cpos(j))<100.0) oss << " ";
    //if(abs(aa.atoms.at(iat).cpos(j))<10.0) oss << " ";
    //if(!std::signbit(aa.atoms.at(iat).cpos(j))) oss << " ";
    //oss << aa.atoms.at(iat).cpos(j) << " ";
    //}
    //}
    //oss << " " << aa.atoms.at(iat).cleanname << " ";
    //if(aa.write_inequivalent_flag==TRUE) {
    //oss << " # ";
    ////	if(i<10) oss << "0";
    //oss << iat << "[";
    //if(aa.atoms.at(iat).equivalent<10) oss << "0";
    //oss << aa.atoms.at(iat).equivalent << "]";
    //if(aa.atoms.at(iat).is_inequivalent) {
    //oss <<"*";
    //oss << "_(" << aa.atoms.at(iat).num_equivalents << ") "; //<< "  index=" << aa.atoms.at(iat).index_iatoms << "  ";
    ////  " v" << aa.iatoms.size() << "   burp ";
    //// for(uint jat=0;jat<aa.iatoms.size();jat++)  oss << aa.iatoms.at(jat).size() << " ";
    //}
    //}
    //oss << endl;
    //}

  }
  return oss.str();
}

// **************************************************************************
// PrintUNCLE
// **************************************************************************
string xstructure::PrintUNCLE(void) {   // Print in uncle format
  stringstream oss;
  this->iomode=IOVASP_POSCAR;
  oss << "# Structure name:" << endl;
  oss << *this;
  return oss.str();
}

// // **************************************************************************
// // PrintADO
// // **************************************************************************
// string xstructure::PrintADO(string strin) {   // Print in ado format
//   stringstream oss;
//   uint _precision_=_DOUBLE_WRITE_PRECISION_MAX_; //14; //was 16 SC 10 DM //CO20180515
//   oss.setf(std::ios::fixed,std::ios::floatfield);
//   oss.precision(_precision_);
//   // oss << "# Structure name number of atoms in the cell:" << endl;
//   //  oss << "1=" << species.at(0) << " 2=" << species.at(1) << " 3=" << species.at(2) << " n=" << species.size() << endl;
//   oss << strin << "# atomic positions in cartesians" << endl;
//   for(uint iat=0;iat<atoms.size();iat++) {
//     oss << strin;
//     for(uint j=1;j<=3;j++) {
//       if(abs(atoms.at(iat).cpos(j))<10.0) oss << " ";if(atoms.at(iat).cpos(j)>=0.0) oss << " "; oss << atoms.at(iat).cpos(j) << " ";
//     }
//     if(species.size()>0) if(atoms.at(iat).name==species.at(0)) oss << " -1";
//     if(species.size()>1) if(atoms.at(iat).name==species.at(1)) oss << "  1";
//     oss << "  " << atoms.at(iat).name;// << " " << atoms.at(iat).type;
//     oss << endl;
//   }  
//   oss << strin << "# Unit cell vectors in units of a" << endl;
//   for(uint ilat=1;ilat<=3;ilat++) {
//     oss << strin;
//     for(uint j=1;j<=3;j++) {
//       if(abs(lattice[ilat][j])<10.0) oss << " ";if(lattice[ilat][j]>=0.0) oss << " "; oss << lattice[ilat][j] << " ";
//     }
//     oss << endl;
//   }
//   return oss.str();
// }


// **************************************************************************
// PAULING DETECTOR
// **************************************************************************
bool PAULING_WyckoffDetector(vector<string> &vinput);

string DeStupidizer(string &strin) {
  aurostd::StringSubst(strin,"\t"," ");
  string s=" ";for(char c=1;c<32;c++) {s[0]=c; aurostd::StringSubst(strin,s," ");}  //destupidization
  aurostd::StringSubst(strin,"  "," "); // cleaning
  return strin;
}

bool sortAtomsTypes(const _atom& a1,const _atom& a2) {  //CO20180705
  //CO20190218
  //sorting on DOUBLES is dangerous, we need to avoid flipping equivalent atoms
  //so we need to set a cutoff
  if(a1.type!=a2.type){return a1.type<a2.type;}
  if(!aurostd::isequal(a1.partial_occupation_value,a2.partial_occupation_value,_AFLOW_POCC_ZERO_TOL_)){return a1.partial_occupation_value>a2.partial_occupation_value;} //reverse, we actually want lowest type but highest occupation first
  return a1.basis<a2.basis; //maintain previous relative ordering
  //if(a1.type!=a2.type){return a1.type<a2.type;}
  //double dist1=aurostd::modulus(a1.fpos);
  //double dist2=aurostd::modulus(a2.fpos);
  //return dist1<dist2;
  //prettier (fpos standard)
  //if(a1.fpos.rows!=3 || a1.fpos.rows!=a2.fpos.rows){
  //  cerr << "XSTRUCTURE::sortAtomsNames:: bad cartesian coordinates" << endl;
  //  throw aurostd::xerror(__AFLOW_FILE__,XPID+"sortAtomsTypes():","Throw for debugging purposes.",_GENERIC_ERROR_);
  //}
  //for(uint i=1;i<=3;i++){
  //  if(a1.fpos[i]!=a2.fpos[i]){return a1.fpos[i]<a2.fpos[i];}
  //}
  //return false;
}

//ideal for AFLOW (alphabetized) + POCC (grouped by occupations)
bool sortAtomsNames(const _atom& a1,const _atom& a2) {  //CO20180705
  //CO20190218
  //sorting on DOUBLES is dangerous, we need to avoid flipping equivalent atoms
  //so we need to set a cutoff
  if(a1.name!=a2.name){return a1.name<a2.name;}
  if(!aurostd::isequal(a1.partial_occupation_value,a2.partial_occupation_value,_AFLOW_POCC_ZERO_TOL_)){return a1.partial_occupation_value>a2.partial_occupation_value;} //reverse, we actually want lowest type but highest occupation first
  return a1.basis<a2.basis; //maintain previous relative ordering
  //return true;
  //if(a1.name!=a2.name){return a1.name<a2.name;}
  //double dist1=aurostd::modulus(a1.fpos);
  //double dist2=aurostd::modulus(a2.fpos);
  //return dist1<dist2;
  //prettier (fpos standard)
  //if(a1.fpos.rows!=3 || a1.fpos.rows!=a2.fpos.rows){
  //  cerr << "XSTRUCTURE::sortAtomsNames:: bad cartesian coordinates" << endl;
  //  throw aurostd::xerror(__AFLOW_FILE__,XPID+"sortAtomsNames():","Throw for debugging purposes.",_GENERIC_ERROR_);
  //}
  //for(uint i=1;i<=3;i++){
  //  if(a1.fpos[i]!=a2.fpos[i]){return a1.fpos[i]<a2.fpos[i];}
  //}
  //return false;
}

bool sortAtomsDist(const _atom& a1,const _atom& a2) {
  //CO20190218
  //sorting on DOUBLES is dangerous, we need to avoid flipping equivalent atoms
  //so we need to set a cutoff
  if(a1.type!=a2.type){return a1.type<a2.type;}
  double dist1=aurostd::modulus(a1.cpos);
  double dist2=aurostd::modulus(a2.cpos);
  if(!aurostd::isequal(dist1,dist2,_ZERO_TOL_)){return dist1<dist2;} //CO20180705 - maybe we need to consider adding tol here
  return sortAtomsTypes(a1,a2); //CO20180705, pocc values!
}

//void sortAtomsDist() {std::stable_sort(atoms.begin(),atoms.end(),sortAtomsDist);}

bool sortAtomsEquiv(const _atom& a1,const _atom& a2){
  if(a1.type!=a2.type){return a1.type<a2.type;} //this is generally implied by equivalent, but not so for POCC, so keep
  if(a1.equivalent!=a2.equivalent){return a1.equivalent<a2.equivalent;}
  return sortAtomsTypes(a1,a2); //CO 180705, pocc values!
} //CO190101

// ---------------------------------------------------------------------------
// Wyckoff sorting function (by Wyckoff letter) //DX20200515
bool sortWyckoffByLetter(const wyckoffsite_ITC& a, const wyckoffsite_ITC& b) {
  // compare letter
  if(a.letter<b.letter){ return true; }
  // if letters the same, sort by type
  else if(a.letter==b.letter){
    if(a.type<b.type){ return true; }
    else { return false; }
  }
  return false;
}

// ---------------------------------------------------------------------------
// Wyckoff sorting function (by atom type) //DX20200515
bool sortWyckoffByType(const wyckoffsite_ITC& a, const wyckoffsite_ITC& b) {
  // compare type
  if(a.type<b.type){ return true; }
  // if types the same, sort by letter
  else if(a.type==b.type){
    if(a.letter<b.letter){ return true; }
    else if(a.letter==b.letter){ //DX20230129 - if the letter is the same
      if(a.parameter_index<b.parameter_index) { return true; } //DX20230130 - check the parameter index
      else { return false; } //DX20230130 - otherwise do not change order
    }
    else { return false; }
  }
  return false;
}

// **************************************************************************
// Xstructure operator>>  INPUT_XSTRUCTURE_INPUT
// **************************************************************************
// loads istream into a xstructure: istream >> xstructure
istream& operator>>(istream& cinput, xstructure& a) {
#define oss cout
  // this is also a constructor so everything should look well defined
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream message;

  if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
  if(LDEBUG){
    if(a.iomode==IOAFLOW_AUTO) cerr << __AFLOW_FUNC__ << " a.iomode = IOAFLOW_AUTO" << endl;
    if(a.iomode==IOVASP_AUTO) cerr << __AFLOW_FUNC__ << " a.iomode = IOVASP_AUTO" << endl;
    if(a.iomode==IOVASP_POSCAR) cerr << __AFLOW_FUNC__ << " a.iomode = IOVASP_POSCAR" << endl;
    if(a.iomode==IOVASP_ABCCAR) cerr << __AFLOW_FUNC__ << " a.iomode = IOVASP_ABCCAR" << endl;
    if(a.iomode==IOVASP_WYCKCAR) cerr << __AFLOW_FUNC__ << " a.iomode = IOVASP_WYCKCAR" << endl;
    if(a.iomode==IOQE_AUTO) cerr << __AFLOW_FUNC__ << " a.iomode = IOQE_AUTO" << endl;
    if(a.iomode==IOQE_GEOM) cerr << __AFLOW_FUNC__ << " a.iomode = IOQE_GEOM" << endl;
    if(a.iomode==IOAIMS_AUTO) cerr << __AFLOW_FUNC__ << " a.iomode = IOAIMS_AUTO" << endl;  //CO20171008
    if(a.iomode==IOAIMS_GEOM) cerr << __AFLOW_FUNC__ << " a.iomode = IOAIMS_GEOM" << endl;  //CO20171008
    if(a.iomode==IOABINIT_GEOM) cerr << __AFLOW_FUNC__ << " a.iomode = IOABINIT_GEOM" << endl;  //DX20200310
    if(a.iomode==IOELK_GEOM) cerr << __AFLOW_FUNC__ << " a.iomode = IOELK_GEOM" << endl;  //DX20200310
    if(a.iomode==IOCIF) cerr << __AFLOW_FUNC__ << " a.iomode = IOCIF" << endl;  //DX20180723
    if(a.iomode==IOATAT_STR) cerr << __AFLOW_FUNC__ << " a.iomode = IOATAT_STR" << endl;  //SD20220114
  }

  if(LDEBUG) cerr << __AFLOW_FUNC__ << " definitions" << endl;
  uint iline=0;
  vector<string> vinput,tokens;
  aurostd::stream2vectorstring(cinput,vinput);
  //CO20180702 - detect NO input
  string input_no_spaces=aurostd::joinWDelimiter(vinput,"");
  input_no_spaces=aurostd::RemoveWhiteSpaces(input_no_spaces);
  if(input_no_spaces.empty()){
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No input",_INPUT_MISSING_);  //CO20190629
  }  //CO20180702

  if(LDEBUG) cerr << __AFLOW_FUNC__ << " DeStupidizer" << endl;
  // now clean for comments, tabs, double spaces ... etc
  //CO20180409 - fixing for issues with # at the beginning of the line
  string::size_type loc;  //CO20180409
  for(uint i=1;i<vinput.size();i++) { // not 1st line
    //aurostd::string2tokens(vinput[i],tokens,"//");  //SC          //CO20180409 - not the best way, as it will screw up for lines that START with comment delimiter
    //if(i>0 && tokens.size()>0) vinput[i]=tokens.at(0);            //CO20180409 - not the best way, as it will screw up for lines that START with comment delimiter
    //aurostd::string2tokens(vinput[i],tokens,"#");  //SC           //CO20180409 - not the best way, as it will screw up for lines that START with comment delimiter
    //if(i>0 && tokens.size()>0) vinput[i]=tokens.at(0);            //CO20180409 - not the best way, as it will screw up for lines that START with comment delimiter
    //aurostd::string2tokens(vinput[i],tokens,"!");  // for QE      //CO20180409 - not the best way, as it will screw up for lines that START with comment delimiter
    //if(i>0 && tokens.size()>0) vinput[i]=tokens.at(0);            //CO20180409 - not the best way, as it will screw up for lines that START with comment delimiter
    if(i>0){  //CO20180409
      loc=vinput[i].find("//");vinput[i]=vinput[i].substr(0,loc);
      loc=vinput[i].find('#');vinput[i]=vinput[i].substr(0,loc);
      loc=vinput[i].find("!");vinput[i]=vinput[i].substr(0,loc);
    }
    DeStupidizer(vinput[i]);
  }
  if(vinput.size()==0) {
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No input",_INPUT_MISSING_);  //CO20190629
  }  //CO20180420
    
  if(LDEBUG){ //CO20230502
    for(uint i=0;i<vinput.size();i++){
      cerr << __AFLOW_FUNC__ << " vinput[i=" << i << "]=\"" << vinput[i] << "\"" << endl;
    }
  }

  //  for(uint i=0;i<vinput.size();i++) cerr << "[" << i << "] " <<  vinput[i] << " " << "[X]" << endl;   throw aurostd::xerror(__AFLOW_FILE__,XPID+"sortAtomsTypes():","Throw for debugging purposes.",_GENERIC_ERROR_);
  string sstring,stmp;
  bool IOMODE_found=FALSE;
  vector<double> poccaus; // partial occupation local host
  a.partial_occupation_sublattice.clear(); // partial occupation local host

  //CO20190219 - need to clear atoms
  //a.atoms.clear();  //CO20190219 - need to use RemoveAtom() (safe)
  //CO20190219 - really remove atoms
  for(uint i=a.atoms.size()-1;i<a.atoms.size();i--){
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " removing atom[" << i << "]" << endl;
    a.RemoveAtom(i);
  }

  // ----------------------------------------------------------------------
  // SEARTH FOR MODES 
  //  LDEBUG=TRUE;
  // PAULING PROTO DETECTOR
  if(!IOMODE_found) {
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " PAULING PROTO DETECTOR" << endl;
    aurostd::string2tokens(vinput.at(0),tokens);
    if(tokens.size() && (tokens.at(0)=="PAULING" || tokens.at(0)=="pauling" || tokens.at(0)=="Pauling")) {  //CO20180420 - ask for size() to print core file
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " PAULING PROTO DETECTOR = TRUE" << endl;
      IOMODE_found=TRUE;
      PAULING_WyckoffDetector(vinput);
      tokens.clear();
    }
  }

  // ----------------------------------------------------------------------
  // QUANTUM ESPRESSO FINDER
  if(!IOMODE_found) {
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " QUANTUM ESPRESSO DETECTOR" << endl;
    uint QE=0;
    bool QE_ERROR=FALSE;
    if(LDEBUG) for(uint i=0;i<vinput.size();i++) cerr << vinput[i] << endl;
    for(uint i=0;i<vinput.size();i++) QE+=aurostd::substring2bool(vinput[i],"&system",true)+aurostd::substring2bool(vinput[i],"&SYSTEM",true); //DX20180123 - added true to clean the spaces in string
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " QUANTUM ESPRESSO DETECTOR QE(&system)=" << QE << endl;
    for(uint i=0;i<vinput.size();i++) QE+=aurostd::substring2bool(vinput[i],"ibrav=",true)+aurostd::substring2bool(vinput[i],"IBRAV=",true); //DX20180123 - added true to clean the spaces in string
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " QUANTUM ESPRESSO DETECTOR QE(ibrav)=" << QE << endl;
    for(uint i=0;i<vinput.size();i++) QE+=aurostd::substring2bool(vinput[i],"nat=",true)+aurostd::substring2bool(vinput[i],"NAT=",true); //DX20180123 - added true to clean the spaces in string
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " QUANTUM ESPRESSO DETECTOR QE(nat)=" << QE << endl;
    for(uint i=0;i<vinput.size();i++) QE+=aurostd::substring2bool(vinput[i],"ntyp=",true)+aurostd::substring2bool(vinput[i],"NTYP=",true); //DX20180123 - added true to clean the spaces in string
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " QUANTUM ESPRESSO DETECTOR QE(ntyp)=" << QE << endl;
    for(uint i=0;i<vinput.size();i++) QE+=aurostd::substring2bool(vinput[i],"atomic_positions",true)+aurostd::substring2bool(vinput[i],"ATOMIC_POSITIONS",true); //DX20180123 - added true to clean the spaces in string
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " QUANTUM ESPRESSO DETECTOR QE(ATOMIC_POSITIONS)=" << QE << endl;

    for(uint i=0;i<vinput.size()&&QE==5;i++) {
      if(aurostd::substring2bool(vinput[i],"ATOMIC_POSITIONS") && 
          !aurostd::substring2bool(vinput[i],"crystal","CRYSTAL") && 
          !aurostd::substring2bool(vinput[i],"bohr","BOHR") && //DX added
          !aurostd::substring2bool(vinput[i],"angstrom","ANGSTROM")) {
        cerr << __AFLOW_FUNC__ << " QE input(1) not supported vinput.at(" << i << ")= \"" << vinput[i] << "\"" << endl;QE_ERROR=TRUE;
      }
      if( aurostd::substring2bool(vinput[i],"&system")) {
        //	cerr << __AFLOW_FUNC__ << " QE input(2) not supported vinput.at(" << i << ")= \"" << vinput[i] << "\"" << endl;QE_ERROR=TRUE;
      }
    }
    if(QE==5 && QE_ERROR) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"QE input errors",_INPUT_MISSING_); //CO20190629
    }
    if(QE==5 && !QE_ERROR) {
      a.iomode=IOQE_AUTO; // might need further discipline but for now it is ok.. 2013 May SC
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " QUANTUM ESPRESSO DETECTOR = TRUE" << endl;
      IOMODE_found=TRUE;
    }
  }

  // ----------------------------------------------------------------------
  //for CIF input //DX20180723 - add cif reader - START
  if(!IOMODE_found) {
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CIF DETECTOR" << endl;
    uint CIF=0;
    if(LDEBUG) for(uint i=0;i<vinput.size();i++) cerr << vinput[i] << endl;
    for(uint i=0;i<vinput.size();i++){ 
      if(aurostd::substring2bool(vinput[i],"loop_",true)){ CIF+=1; break;}
    }
    for(uint i=0;i<vinput.size();i++){ 
      if(aurostd::substring2bool(vinput[i],"_atom_site_",true)){ CIF+=1; break;}
    }
    for(uint i=0;i<vinput.size();i++){ 
      if(aurostd::substring2bool(vinput[i],"_cell_length_",true)){ CIF+=1; break;}
    }
    if(CIF==3){
      a.iomode=IOCIF; 
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " CIF DETECTOR = TRUE" << endl;
      IOMODE_found=TRUE;
    }
  }
  //DX20180723 - add cif reader - END

  // ----------------------------------------------------------------------
  // ABINIT input - START (DX20200310)
  // based on documentation found at https://docs.abinit.org/variables/basic/
  // NOTE: originally had acell and rprim as keywords, but they are not mandatory
  if(!IOMODE_found) {
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " ABINIT DETECTOR" << endl;
    uint ABINIT=0;
    // find number of atoms (natom)
    for(uint i=0;i<vinput.size();i++){ 
      if(aurostd::substring2bool(vinput[i],"natom",true)){ ABINIT+=1; break;} 
    }
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " ABINIT DETECTOR (natom)=" << ABINIT << endl;
    // find types of atoms (typat)
    for(uint i=0;i<vinput.size();i++){ 
      if(aurostd::substring2bool(vinput[i],"typat",true)){ ABINIT+=1; break;} 
    }
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " ABINIT DETECTOR (typat)=" << ABINIT << endl;
    // find atom positions (xred, xcart, xangst)
    for(uint i=0;i<vinput.size();i++){ 
      if(aurostd::substring2bool(vinput[i],"xred",true) || 
          aurostd::substring2bool(vinput[i],"xcart",true) || 
          aurostd::substring2bool(vinput[i],"xangst",true)){ ABINIT+=1; break;} 
    }
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " ABINIT DETECTOR (xred,xcart,xangst)=" << ABINIT << endl;

    if(ABINIT==3){
      a.iomode = IOABINIT_GEOM;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " ABINIT DETECTOR = TRUE" << endl;
      IOMODE_found = TRUE;
    }
  }
  // ABINIT input - END (DX20200310)

  // ----------------------------------------------------------------------
  // ELK input - START (DX20200310)
  if(!IOMODE_found) {
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " ELK DETECTOR" << endl;
    uint ELK=0;
    // find atoms keyword (atoms)
    for(uint i=0;i<vinput.size();i++){
      if(aurostd::substring2bool(vinput[i],"atoms",true)){ ELK+=1; break;}
    }
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " ELK DETECTOR (atoms)=" << ELK << endl;
    // find lattice keyword (avec)
    for(uint i=0;i<vinput.size();i++){
      if(aurostd::substring2bool(vinput[i],"avec",true)){ ELK+=1; break;}
    }
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " ELK DETECTOR (avec)=" << ELK << endl;

    if(ELK==2){
      a.iomode = IOELK_GEOM;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " ELK DETECTOR = TRUE" << endl;
      IOMODE_found = TRUE;
    }
  }
  // ELK input - END (DX20200310)

  // ----------------------------------------------------------------------
  // ATAT input - START (SD20220117)
  if(!IOMODE_found) {
    if (LDEBUG) cerr << __AFLOW_FUNC__ << " ATAT DETECTOR" << endl; 
    uint ATAT = 1;
    // count number of entries for axes and fractional cell vectors, check for correct type
    uint line = 0;
    for (; line < vinput.size() - 1 && line < 6 && ATAT; line++) {
      aurostd::string2tokens(vinput[line], tokens, " ");
      ATAT = (tokens.size() == 3 && aurostd::isfloat(tokens[0]) && aurostd::isfloat(tokens[1]) && aurostd::isfloat(tokens[2]));
    }
    // count number of entries for atom positions, check for correct type
    for (; line < vinput.size() - 1 && ATAT; line++) {
      aurostd::string2tokens(vinput[line],tokens," ");
      ATAT = (tokens.size() == 4 && aurostd::isfloat(tokens[0]) && aurostd::isfloat(tokens[1]) && aurostd::isfloat(tokens[2]));
      if (ATAT && !xelement::xelement::isElement(tokens[3])) {ATAT = 0;}
    }
    if (ATAT == 1) {
      a.iomode = IOATAT_STR;
      if (LDEBUG) cerr << __AFLOW_FUNC__ << " ATAT DETECTOR = TRUE" << endl; 
      IOMODE_found = TRUE; 
    }
  }
  // ATAT input - END (SD20220117)

  // ----------------------------------------------------------------------
  //for AIMS input - unfortunately, it's very generic so leave for last
  if(!IOMODE_found) {
    vector<string> tokens_line;
    //since it's so generic, let's be super strict, only look at the first word in the line
    //needs to match atom, atom_frac, or lattice
    for(uint i=0;i<vinput.size()&&!IOMODE_found;i++){
      aurostd::string2tokens(vinput[i],tokens_line," ");
      if(tokens_line.size()){
        IOMODE_found=(IOMODE_found || 
            tokens_line[0]=="atom" || 
            tokens_line[0]=="atom_frac" || 
            tokens_line[0]=="lattice_vector");
      }
    }
    if(IOMODE_found){
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " AIMS GEOM DETECTOR = TRUE" << endl;}
      a.iomode=IOAIMS_GEOM;
      //if atom_frac, then lattice MUST be provided
      //if atom, no need for lattice
      bool lat_found,lat_found_anywhere,lat1_found,lat2_found,lat3_found;
      bool atom_found,atom_found_anywhere,frac_found_anywhere;
      lat_found=lat_found_anywhere=lat1_found=lat2_found=lat3_found=atom_found=atom_found_anywhere=frac_found_anywhere=false;
      for(uint i=0;i<vinput.size();i++){
        tokens_line.clear();
        aurostd::string2tokens(vinput[i],tokens_line," ");
        if(!tokens_line.size()){continue;}
        lat_found=(tokens_line[0]=="lattice_vector");
        lat_found_anywhere=(lat_found_anywhere || lat_found);
        frac_found_anywhere=(frac_found_anywhere || tokens_line[0]=="atom_frac");
        atom_found=(aurostd::substring2bool(tokens_line[0],"atom"));
        atom_found_anywhere=(atom_found_anywhere || atom_found || frac_found_anywhere);
        if(LDEBUG) {
          cerr << __AFLOW_FUNC__ << " AIMS line[" << i+1 << "]: " << vinput[i] << endl;
          cerr << __AFLOW_FUNC__ << " AIMS lat_found=" << lat_found << endl;
          cerr << __AFLOW_FUNC__ << " AIMS lat_found_anywhere=" << lat_found_anywhere << endl;
          cerr << __AFLOW_FUNC__ << " AIMS lat1_found=" << lat1_found << endl;
          cerr << __AFLOW_FUNC__ << " AIMS lat2_found=" << lat2_found << endl;
          cerr << __AFLOW_FUNC__ << " AIMS lat3_found=" << lat3_found << endl;
          cerr << __AFLOW_FUNC__ << " AIMS atom_found=" << atom_found << endl;
          cerr << __AFLOW_FUNC__ << " AIMS atom_found_anywhere=" << atom_found_anywhere << endl;
          cerr << __AFLOW_FUNC__ << " AIMS frac_found_anywhere=" << frac_found_anywhere << endl;
        }
        if(lat_found || atom_found){
          if(lat_found && tokens_line.size()<4){  //could be more, but not less
            message << " AIMS input error, ";  //CO20190629
            message << "lattice_vector "; //CO20190629
            message << "at line[" << i+1 << "] is ill-defined" << endl; //CO20190629
            message << "line: " << vinput[i] << endl; //CO20190629
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
          }
          if(atom_found && tokens_line.size()<5){  //could be more, but not less (need name/type in last column)
            message << " AIMS input error, "; //CO20190629
            message << "atom position "; //CO20190629
            message << "at line[" << i+1 << "] "; //CO20190629
            if(tokens_line.size()==4){message << "is missing the atom name" << endl;} //CO20190629
            else {message << "is ill-defined" << endl;} //CO20190629
            message << "line: " << vinput[i] << endl;  //CO20190629
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
          }
        }
      }
      if(!atom_found_anywhere){
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"AIMS input error, no atoms found",_INPUT_ERROR_); //CO20190629
      }
      a.coord_flag=_COORDS_CARTESIAN_;
      if(lat_found_anywhere || frac_found_anywhere){
        a.coord_flag=_COORDS_FRACTIONAL_;
        uint lat_count=0;
        vector<string> lat;
        for(uint i=0;i<vinput.size();i++){
          aurostd::string2tokens(vinput[i],tokens_line," ");
          if(tokens_line.size() && tokens_line[0]=="lattice_vector"){
            lat_count++;
            if(LDEBUG) {cerr << __AFLOW_FUNC__ << " AIMS lat_count=" << lat_count << " line[" << i << "," << vinput[i] << "]" << endl;}
            if(lat_count==1){lat1_found=true;}
            else if(lat_count==2){lat2_found=true;}
            else if(lat_count==3){lat3_found=true;}
            else {
              throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"AIMS input error, too many lattice vectors found",_INPUT_ERROR_);  //CO20190629
            }
          }
        }
        if(!lat1_found || !lat2_found || !lat3_found){
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"AIMS input error, incomplete lattice vector specification (needed if atom_frac found)",_INPUT_ERROR_);  //CO20190629
        }
      }
    }
    else {if(LDEBUG) {cerr << __AFLOW_FUNC__ << " AIMS GEOM DETECTOR = FALSE" << endl;}}
  }
  // DESPERATE FINDING => VASP
  if(!IOMODE_found) {
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " VASP DETECTOR" << endl;
    if(a.iomode==IOAFLOW_AUTO) {
      a.iomode=IOVASP_AUTO; // still not found, need something to eat.
      //  if(a.iomode==IOVASP_AUTO) a.iomode=IOVASP_POSCAR;
      IOMODE_found=TRUE;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " VASP DETECTOR = TRUE" << endl;
    }
  }

  // ----------------------------------------------------------------------
  // SOME EXTRA VERBOSE
  if(LDEBUG) if(a.iomode==IOAFLOW_AUTO) cerr << __AFLOW_FUNC__ << " a.iomode = IOAFLOW_AUTO" << endl;
  if(LDEBUG) if(a.iomode==IOVASP_AUTO) cerr << __AFLOW_FUNC__ << " a.iomode = IOVASP_AUTO" << endl;
  if(LDEBUG) if(a.iomode==IOVASP_POSCAR) cerr << __AFLOW_FUNC__ << " a.iomode = IOVASP_POSCAR" << endl;
  if(LDEBUG) if(a.iomode==IOVASP_ABCCAR) cerr << __AFLOW_FUNC__ << " a.iomode = IOVASP_ABCCAR" << endl;
  if(LDEBUG) if(a.iomode==IOVASP_WYCKCAR) cerr << __AFLOW_FUNC__ << " a.iomode = IOVASP_WYCKCAR" << endl;
  if(LDEBUG) if(a.iomode==IOQE_AUTO) cerr << __AFLOW_FUNC__ << " a.iomode = IOQE_AUTO" << endl;
  if(LDEBUG) if(a.iomode==IOQE_GEOM) cerr << __AFLOW_FUNC__ << " a.iomode = IOQE_GEOM" << endl;
  if(LDEBUG) if(a.iomode==IOAIMS_AUTO) cerr << __AFLOW_FUNC__ << " a.iomode = IOAIMS_AUTO" << endl;  //CO20171008
  if(LDEBUG) if(a.iomode==IOAIMS_GEOM) cerr << __AFLOW_FUNC__ << " a.iomode = IOAIMS_GEOM" << endl;  //CO20171008
  if(LDEBUG) if(a.iomode==IOATAT_STR) cerr << __AFLOW_FUNC__ << " a.iomode = IOATAT_STR" << endl;  //SD20220114
  // ----------------------------------------------------------------------
  // VASP INPUT
  if(a.iomode==IOVASP_AUTO || a.iomode==IOVASP_POSCAR || a.iomode==IOVASP_ABCCAR || a.iomode==IOVASP_WYCKCAR) { // VASP POSCAR
    // for variable number of items
    //bool scale_second_flag=FALSE;//,scale_third_flag=FALSE; //CO20180409
    //double scale_second_value=0.0;//,scale_third_value=0.0;; //CO20180409
    //
    a.is_vasp4_poscar_format=FALSE;
    a.is_vasp5_poscar_format=FALSE;
    // -------------- POSCAR EXISTS
    int num_atoms=0;
    // -------------- TITLE
    // input.getline(stmp,MAX_TITLE_SIZE);title=stmp;
    if(vinput.size()-1<iline) {
      message << "missing line[" << iline << "]" << endl; //CO20190629
      for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
    }  //CO20180420 - check for missing lines
    //DX+ME20210525 [OBSOLETE - need to remove control code characters from input, important for web] a.title=vinput.at(iline++);
    aurostd::RemoveControlCodeCharactersFromString(vinput[iline++],a.title); //DX+ME20210525 - need to remove control code characters from input, important for web
    // -------------- SCALE
    //    input >> a.scale;
    if(vinput.size()-1<iline) {
      message << "missing line[" << iline << "]" << endl; //CO20190629
      for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
    }  //CO20180420 - check for missing lines
    stmp=vinput.at(iline++);
    aurostd::StringSubst(stmp,"\t"," ");aurostd::StringSubst(stmp,"  "," ");aurostd::StringSubst(stmp,"  "," ");
    aurostd::string2tokens(stmp,tokens);
    if(tokens.size()==0) {
      message << "missing second line in poscar" << endl; //CO20190629
      for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
    }
    // oss << tokens.size() <<  " = " << tokens.at(0) << endl;throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_);
    a.scale=0.0;
    if(tokens.size()>0) {a.scale=aurostd::string2utype<double>(tokens.at(0));}
    if(tokens.size()>1) {/*a.neg_scale_second=TRUE;*/a.scale_second=aurostd::string2utype<double>(tokens.at(1));a.neg_scale_second=std::signbit(a.scale_second);} //CO20180409
    if(tokens.size()>2) {a.scale_third.isentry=TRUE;a.scale_third.content_double=aurostd::string2utype<double>(tokens.at(2));}  //CO20170803 - site tol
    //  oss << a.scale << endl;
    // if(a.neg_scale_second) oss << a.scale_second << endl; //CO20180409
    // if(scale_third_flag) oss << scale_third_value << endl;
    // oss << sstring << endl;
    // -------------- UNIT CELL
    if(vinput.size()-1<iline) {
      message << "missing line[" << iline << "]" << endl; //CO20190629
      for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
    }  //CO20180420 - check for missing lines
    stmp=vinput.at(iline++);
    aurostd::string2tokens(stmp,tokens);
    a.iomode=IOVASP_POSCAR;
    if(tokens.size()==3) a.iomode=IOVASP_POSCAR;
    if(tokens.size()==6) a.iomode=IOVASP_ABCCAR;
    if(tokens.size()==7 || tokens.size()==8) a.iomode=IOVASP_WYCKCAR;
    //    oss << " LDEBUG token.size()=" << tokens.size() << "" << endl;
    // ---------------------------------------------------------------
    if(a.iomode==IOVASP_POSCAR) {
      // oss << __AFLOW_FUNC__ << " Chosen IOVASP_POSCAR" << endl;
      //    input >> a.lattice(1,1) >> a.lattice(1,2) >> a.lattice(1,3);
      a.lattice(1,1)=aurostd::string2utype<double>(tokens[0]);
      a.lattice(1,2)=aurostd::string2utype<double>(tokens[1]);
      a.lattice(1,3)=aurostd::string2utype<double>(tokens[2]);
      stringstream input_tmp;
      if(vinput.size()-1<iline) {
        message << "missing line[" << iline << "]" << endl; //CO20190629
        for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
      }  //CO20180420 - check for missing lines
      input_tmp.clear();input_tmp.str(vinput.at(iline++));
      input_tmp >> a.lattice(2,1) >> a.lattice(2,2) >> a.lattice(2,3);
      if(vinput.size()-1<iline) {
        message << "missing line[" << iline << "]" << endl; //CO20190629
        for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
      }  //CO20180420 - check for missing lines
      input_tmp.clear();input_tmp.str(vinput.at(iline++));
      input_tmp >> a.lattice(3,1) >> a.lattice(3,2) >> a.lattice(3,3);
      xvector<double> data(6);
      data=Getabc_angles(a.lattice,DEGREES);
      a.a=data[1];a.b=data[2];a.c=data[3];
      a.alpha=data[4];a.beta=data[5];a.gamma=data[6];
    }
    // ---------------------------------------------------------------
    if(a.iomode==IOVASP_ABCCAR) {
      // cerr << __AFLOW_FUNC__ << " Chosen IOVASP_ABCCAR" << endl;
      a.a=aurostd::string2utype<double>(tokens[0]);
      a.b=aurostd::string2utype<double>(tokens[1]);
      a.c=aurostd::string2utype<double>(tokens[2]);
      a.alpha=aurostd::string2utype<double>(tokens[3]);
      a.beta=aurostd::string2utype<double>(tokens[4]);
      a.gamma=aurostd::string2utype<double>(tokens[5]);
      a.lattice=GetClat(a.a,a.b,a.c,a.alpha,a.beta,a.gamma);
    }
    // ---------------------------------------------------------------
    if(a.iomode==IOVASP_WYCKCAR) {
      //     cerr << __AFLOW_FUNC__ << " Chosen IOVASP_WYCKCAR" << endl;
      //DX20170905 - GET SYM_EPS
      vector<string> title_tokens;
      aurostd::string2tokens(a.title,title_tokens);
      if(aurostd::substring2bool(title_tokens,"sym_eps")){
        a.sym_eps=aurostd::string2utype<double>(title_tokens[title_tokens.size()-1]);     
      } 
      //DX20170905 - GET SYM_EPS
      a.a=aurostd::string2utype<double>(tokens[0]);
      a.b=aurostd::string2utype<double>(tokens[1]);
      a.c=aurostd::string2utype<double>(tokens[2]);
      a.alpha=aurostd::string2utype<double>(tokens[3]);
      a.beta=aurostd::string2utype<double>(tokens[4]);
      a.gamma=aurostd::string2utype<double>(tokens[5]);
      a.spacegroupnumber=aurostd::string2utype<int>(tokens[6]);
      a.spacegroupnumberoption=0;// 1; no option
      if(tokens.size()>=8) a.spacegroupnumberoption=aurostd::string2utype<int>(tokens[7]);
      a.spacegrouplabel=GetSpaceGroupLabel(a.spacegroupnumber);
      a.spacegroup=GetSpaceGroupName(a.spacegroupnumber, a.directory); //DX20180526 - add directory
      a.lattice=GetClat(a.a,a.b,a.c,a.alpha,a.beta,a.gamma);
    }
    // If scale < 0 then it should be treated as the volume.
    a.neg_scale=FALSE;
    if(a.scale<0.0) {
      //CO20211130 - the usual scaling factor is the lattice parameter
      //for fcc, GetVol(a.lattice)=0.25
      double nvol=-1.0*(a.scale);
      double ovol=GetVol(a.lattice);
      a.scale=std::pow((double) (nvol/ovol),(double) 1.0/3.0);
      a.neg_scale=TRUE;
    }
    a.FixLattices(); // Reciprocal/f2c/c2f
    a.kpoints_k1=0;a.kpoints_k2=0;a.kpoints_k3=0;
    a.kpoints_kmax=0;a.kpoints_kppra=0;a.kpoints_kscheme="";
    clear(a.origin);
    // ---------------------------------------------------------------
    // -------------- CHECK VASP4 VASP5 and CHECK DIRECT/CARTESIANS STUFF
    //CO20190629 - shift to line with Direct/Cartesian to find POCC specification
    //don't worry about setting Selective Dynamics yet (but be aware of it)
    if(a.iomode==IOVASP_POSCAR) stmp=vinput[6]; //true for VASP4, change if VASP5 later
    if(a.iomode==IOVASP_ABCCAR) stmp=vinput[4];
    if(a.iomode==IOVASP_WYCKCAR) stmp=vinput[4];
    aurostd::StringSubst(stmp,"\t"," ");aurostd::StringSubst(stmp,"  "," ");aurostd::StringSubst(stmp,"  "," ");
    aurostd::string2tokens(stmp,tokens);
    if(tokens.size()==0) {
      //[CO20190629 - OBSOLETE]message << "missing D/C/S line" << endl;  //CO20190629
      message << "Missing \"Selective Dynamics\"/\"Direct\"/\"Cartesian\" line" << endl;  //CO20190629
      for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_MISSING_); //CO20190629
    }

    // VASP 4
    a.is_vasp4_poscar_format=TRUE;
    a.is_vasp5_poscar_format=FALSE;

    if(!(tokens[0][0]=='S' || tokens[0][0]=='s' || tokens[0][0]=='D' || tokens[0][0]=='d' || tokens[0][0]=='C' || tokens[0][0]=='c')){
      // VASP 5
      a.is_vasp4_poscar_format=FALSE;
      a.is_vasp5_poscar_format=TRUE;
      //Kesong Yang identified this bug, we need to shift the line if vasp5, also if selective dynamics
      uint ijline=0;
      if(a.iomode==IOVASP_POSCAR) ijline=7;
      if(a.iomode==IOVASP_ABCCAR) ijline=5;
      if(a.iomode==IOVASP_WYCKCAR) ijline=5;
      if(tokens[0][0]=='S' || tokens[0][0]=='s'){ijline++;} //selective dynamics
      stmp=vinput.at(ijline);
      aurostd::StringSubst(stmp,"\t"," ");aurostd::StringSubst(stmp,"  "," ");aurostd::StringSubst(stmp,"  "," ");
      aurostd::string2tokens(stmp,tokens);  //now we should have D/C line
    }

    if(!(tokens[0][0]=='S' || tokens[0][0]=='s' || tokens[0][0]=='D' || tokens[0][0]=='d' || tokens[0][0]=='C' || tokens[0][0]=='c')){
      message << "Missing \"Selective Dynamics\"/\"Direct\"/\"Cartesian\" line" << endl;  //CO20190629
      for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_MISSING_); //CO20190629
    }

    //[CO20190629 - OBSOLETE]if(tokens[0][0]=='S' || tokens[0][0]=='s' || tokens[0][0]=='D' || tokens[0][0]=='d' || tokens[0][0]=='C' || tokens[0][0]=='c') {
    //[CO20190629 - OBSOLETE]  // VASP 4
    //[CO20190629 - OBSOLETE]  a.is_vasp4_poscar_format=TRUE;
    //[CO20190629 - OBSOLETE]  a.is_vasp5_poscar_format=FALSE;
    //[CO20190629 - OBSOLETE]} else {
    //[CO20190629 - OBSOLETE]  // VASP 5
    //[CO20190629 - OBSOLETE]  a.is_vasp4_poscar_format=FALSE;
    //[CO20190629 - OBSOLETE]  a.is_vasp5_poscar_format=TRUE;
    //[CO20190629 - OBSOLETE]}

    a.partial_occupation_flag=FALSE;

    for(uint j=1;j<tokens.size();j++) if(tokens[j][0]=='P' || tokens[j][0]=='p') a.partial_occupation_flag=TRUE;
    //if(a.scale_second==0.0) { //CO20180409
    a.partial_occupation_HNF=0; // nothing defined //CO20180409
    a.partial_occupation_site_tol=DEFAULT_POCC_SITE_TOL; //DEFAULT_PARTIAL_OCCUPATION_TOLERANCE; // nothing defined  //CO20170803 - site tol
    a.partial_occupation_stoich_tol=DEFAULT_POCC_STOICH_TOL; //DEFAULT_PARTIAL_OCCUPATION_TOLERANCE; // nothing defined //CO20180409
    //}
    if(a.partial_occupation_flag) //&& a.neg_scale_second) //CO20180409
    { //CO20200106 - patching for auto-indenting
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " a.neg_scale_second=" << a.neg_scale_second << ", a.scale_second=" << a.scale_second << endl;}
      if(a.neg_scale_second){a.partial_occupation_HNF=(int) (-a.scale_second);} //CO20180409
      else {
        a.partial_occupation_site_tol=a.partial_occupation_stoich_tol=a.scale_second; //CO20180409
        if(a.scale_third.isentry){a.partial_occupation_stoich_tol=a.scale_third.content_double;} //CO20170803 - site tol
      }
      if(LDEBUG) {
        cerr << __AFLOW_FUNC__ << " a.partial_occupation_HNF=" << a.partial_occupation_HNF << endl;
        cerr << __AFLOW_FUNC__ << " a.partial_occupation_site_tol=" << a.partial_occupation_site_tol << endl;
        cerr << __AFLOW_FUNC__ << " a.partial_occupation_stoich_tol=" << a.partial_occupation_stoich_tol << endl;
      }
      //if(a.scale_second==0.0) {//CO20180409
      //  a.partial_occupation_stoich_tol=DEFAULT_PARTIAL_OCCUPATION_TOLERANCE; // nothing defined
      //  a.partial_occupation_site_tol=DEFAULT_PARTIAL_OCCUPATION_TOLERANCE; // nothing defined  //CO20170803 - site tol
      //}
      //     if(abs(a.partial_occupation_stoich_tol)>1e-12) cerr << "a.partial_occupation_stoich_tol=" << a.partial_occupation_stoich_tol << endl; //CO20180409
      //   if(abs(a.partial_occupation_HNF)>0.0) cerr << "a.partial_occupation_HNF=" << a.partial_occupation_HNF << endl;
    }
    //if(a.partial_occupation_flag && a.scale_third.isentry){a.partial_occupation_site_tol=a.scale_third.content_double;} //CO20170803 - site tol

    //last line was last lattice vector/geometry line
    if(a.is_vasp5_poscar_format) {
      if(vinput.size()-1<iline) {
        message << "missing line[" << iline << "]" << endl; //CO20190629
        for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
      }  //CO20180420 - check for missing lines
      stmp=vinput.at(iline++);  // to skip toward vasp5
    }
    // ---------------------------------------------------------------
    // -------------- ATOMS
    // Number of atoms of each type and total number of atoms.
    if(vinput.size()-1<iline) {
      message << "missing line[" << iline << "]" << endl; //CO20190629
      for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
    }  //CO20180420 - check for missing lines
    stmp=vinput.at(iline++);
    // The following is necessary because if the last lattice parameter has
    // tabs/spaces after it then the getline just grabs those spaces.  This
    // second getline will then be called to actually read in the number
    // of each type of atom.
    string tmpns;
    tmpns=aurostd::RemoveWhiteSpaces(stmp);
    if(string(tmpns).size()==0) {
      if(vinput.size()-1<iline) {
        message << "missing line[" << iline << "]" << endl; //CO20190629
        for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
      }  //CO20180420 - check for missing lines
      stmp=vinput.at(iline++);
    }
    a.num_each_type.clear();num_atoms=0;
    aurostd::string2tokens(stmp,tokens);
    // -------------- check for partial occupation
    vector<string> tokens_i=tokens,tokens_j,tokens_k;
    int number;double dpocc;
    //    a.partial_occupation_flag=FALSE;
    for(uint i=0;i<tokens_i.size()&&!a.partial_occupation_flag;i++)
      if(aurostd::substring2bool(tokens_i.at(i),"*") || aurostd::substring2bool(tokens_i.at(i),"+"))
        a.partial_occupation_flag=TRUE;
    // -------------- no partial occupation
    if(a.partial_occupation_flag==FALSE) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " PARTIAL OCCUPATION = FALSE" << endl;
      for(uint i=0;i<tokens_i.size();i++) {
        number=aurostd::string2utype<int>(tokens_i.at(i));dpocc=1.0;
        a.num_each_type.push_back(number); num_atoms+=number;
        for(uint l=0;l<(uint) number;l++) poccaus.push_back(dpocc);
        for(uint l=0;l<(uint) number;l++) a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
      }
    }
    // -------------- yes partial occupation
    if(a.partial_occupation_flag==TRUE) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " PARTIAL OCCUPATION = TRUE" << endl;
      for(uint i=0;i<tokens_i.size();i++) {
        if(!aurostd::substring2bool(tokens_i.at(i),"*") && !aurostd::substring2bool(tokens_i.at(i),"+")) { // NO POCC KEYWORD
          number=aurostd::string2utype<int>(tokens_i.at(i));dpocc=1.0;
          a.num_each_type.push_back(number); num_atoms+=number;
          for(uint l=0;l<(uint) number;l++) poccaus.push_back(dpocc);
          for(uint l=0;l<(uint) number;l++) a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
        } else {
          aurostd::string2tokens(tokens_i.at(i),tokens_j,"+");if(LDEBUG) cerr << "tokens_j.size()=" << tokens_j.size() << endl;
          number=0;int nnumber;
          for(uint j=0;j<tokens_j.size();j++) {
            aurostd::string2tokens(tokens_j.at(j),tokens_k,"*");
            if(tokens_k.size()==0) {
              message << "PARTIAL OCCUPATION error [1] tokens_k.size()==0, no *" << endl; //CO20190629
              for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
              throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
            }
            nnumber=aurostd::string2utype<int>(tokens_k.at(0));
            if(tokens_k.size()==1) dpocc=1.0;
            if(tokens_k.size()==2) {dpocc=aurostd::string2utype<double>(tokens_k.at(1));} // SOMETHING FOR ROUNDING DPOCC UP TO FRACTIONS....
            number+=nnumber; num_atoms+=nnumber;
            for(uint l=0;l<(uint) nnumber;l++) poccaus.push_back(dpocc);
            for(uint l=0;l<(uint) nnumber;l++) {
              if(tokens_k.size()==1) a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_); // no specie number
              if(tokens_k.size()==2) a.partial_occupation_sublattice.push_back(i); // put specie number
            }
            if(tokens_k.size()>=3) {
              message << "PARTIAL OCCUPATION error [1] tokens_k.size()>=3, too many *" << endl; //CO20190629
              for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
              throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
            }
          } // loop on +
          a.num_each_type.push_back(number);
        }
      }
      if(LDEBUG) {cerr << "P(" << poccaus.size()<< ") = ";for(uint j=0;j<poccaus.size();j++) cerr << poccaus.at(j) << " ";cerr << endl;}
    }
    //  cerr << "num_atoms=" << num_atoms << endl; throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_); // num_atoms is the SUM of the atoms in the numbers
    // -------------- COORDINATE TYPES
    // Type of coordinates (Fractional or Cartesian) - only 1st character matters (D/d or C/c).
    // This line might also be the Selective Dynamics line so you must check for that (S/s).
    a.isd=FALSE;
    string stmp;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " DEBUG [1]" << endl;
    if(vinput.size()-1<iline) {
      message << "missing line[" << iline << "]" << endl; //CO20190629
      for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
    }  //CO20180420 - check for missing lines
    stmp=vinput.at(iline++);
    aurostd::StringSubst(stmp,"\t"," ");
    std::vector<string> stmp_tokens;
    aurostd::string2tokens(stmp,stmp_tokens);
    // Note that if there are spaces at the beginning of the line we have to remove them.
    if(stmp_tokens.size()==0) {
      message << "Found blank line on line 7. This line should give coordinate type or selective dynamics." << endl;  //CO20190629
      for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
    } else {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " DEBUG [2]" << endl;
      string sstmp=stmp_tokens.at(0);
      a.order_parameter_structure=FALSE;
      if(sstmp[0]=='S' || sstmp[0]=='s') {
        a.isd=TRUE;
        if(vinput.size()-1<iline){
          message << "missing line[" << iline << "]" << endl; //CO20190629
          for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
        }  //CO20180420 - check for missing lines
        stmp=vinput.at(iline++);
        sstmp=aurostd::RemoveSpaces(stmp);
      }
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " DEBUG [3]" << endl;
      if(sstmp[0]=='D' || sstmp[0]=='d') {
        //	cerr << "FRAC" << endl;
        a.coord_type[0]=sstmp[0];
        a.coord_flag=_COORDS_FRACTIONAL_;
      } else {
        if(sstmp[0]=='C' || sstmp[0]=='c') {
          //	  cerr << "CART" << endl;
          a.coord_type[0]=sstmp[0];
          a.coord_flag=_COORDS_CARTESIAN_;
          if(a.iomode==IOVASP_WYCKCAR) {
            message << "WYCKOFF mode requires FRACTIONAL coordinates (DIRECT)." << endl;  //CO20190629
            for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
          }
        } else {
          message << "Did not find coordinate type D/d or C/c." << endl;  //CO20190629
          for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
        }
      }
      a.coord_type[1]='\0';
      if(stmp_tokens.size()>1) {
        for(uint ipos=1;ipos<stmp_tokens.size();ipos++) {
          if(stmp_tokens.at(ipos)[0]=='O' || stmp_tokens.at(ipos)[0]=='o') a.order_parameter_structure=TRUE;
          //	  if(stmp_tokens.at(ipos)[0]=='P' || stmp_tokens.at(ipos)[0]=='p') a.partial_occupation_flag=TRUE;
        }
      }
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " DEBUG a.order_parameter_structure=" << a.order_parameter_structure << endl;
    }
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " DEBUG [4]" << endl;

    // --------------  Basis atom positions and names
    //   int cnt;
    xvector<double> v(3);clear(v);
    xvector<double> vzero(3);clear(vzero);
    // Set default values: names, sd.
    string def_name="H";
    string def_sd = "TTT";
    // NOTHING TO CLEAR     clear(*cpos);clear(*fpos);clear(*types);
    //  Read in all the atom data.
    // The number of fields on a line can vary but must mean something, which depends on isd.
    // For isd=0:
    //   If nf=3 then there are only atomic positions.
    //   If nf>=4 then there are atomic positions and a name.
    // For isd!=0:
    //   If nf=3 then there are only atomic positions.
    //   If nf=4 then there are atomic positions and a name.
    //   If nf=6 then there are atomic positions and sd.
    //   If nf>=7 then there are atomic positions and sd and a name.

    // clear up a little
    // NOTHING TO CLEAR     for(int i=0;i<=MAX_ATOM_UCELL;i++) (*name)[i].clear();          // vector(0,MAX_ATOM_UCELL)
    // NOTHING TO CLEAR     for(int i=0;i<=MAX_ATOM_UCELL;i++) (*sd)[i].clear();             // vector(0,MAX_ATOM_UCELL)
    // NOTHING TO CLEAR     for(int i=1;i<=MAX_ATOM_UCELL;i++) (*name_is_given)(i)=FALSE; // xvector(1,MAX_ATOM_UCELL)

    int iat=-1;uint itype=0,j;uint iline_ref=iline;
    for(itype=0;itype<a.num_each_type.size();itype++) {
      //  cerr << a.num_each_type.at(itype) << endl;
      for(j=0;j<(uint) a.num_each_type.at(itype);j++) {
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " DEBUG [5] itype,j=" << itype << "," << j << endl;
        //	bool plug_atom=TRUE;
        iat++; // it startf from -1 so the first is ZERO
        //  cerr << iat << " type=" << itype << endl;
        //  for(int iat=1;iat<=(num_atoms);iat++) { //[CO20200106 - close bracket for indenting]}
        // cerr << iat << " " << (num_atoms) << endl;
        _atom atom;                                        // create new atom
        string stmp;
        stmp=vinput.at(iline++);
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " " << iline << " " << vinput.size() << "," << iline-iline_ref << "," << num_atoms << "," << stmp << endl;
        if(iline==vinput.size() && (iline-iline_ref<(uint) num_atoms)) {
          message << "Insufficient number of atom lines." << endl;  //CO20190629
          for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
        }  
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " DEBUG [6]" << endl;
        stmp=aurostd::RemoveCharacter(stmp,'\t');
        aurostd::StringSubst(stmp,"\t"," ");
        std::vector<string> stmp_tokens;
        aurostd::string2tokens(stmp,stmp_tokens);
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " DEBUG [6b] stmp_tokens.size()=" << stmp_tokens.size() << endl;
        if(stmp_tokens.size()<3) {
          message << "Insufficient number of atom entries in atom=" << iline-iline_ref << "" << endl; //CO20190629
          for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
        }  

        int id=0;
        string name;
        // Read in the atom positions.
        v(1)=atof(stmp_tokens.at(id++).c_str()); // when not a matter of speed, vectors
        v(2)=atof(stmp_tokens.at(id++).c_str()); // should be indicized with () instead of [],
        v(3)=atof(stmp_tokens.at(id++).c_str()); // b/c () does boundary checks !
        if(a.coord_flag==_COORDS_FRACTIONAL_) {
          //    for(int ii=1;ii<=3;ii++)
          atom.fpos=v;
          atom.cpos=F2C(a.lattice,atom.fpos);
        }
        if(a.coord_flag==_COORDS_CARTESIAN_) {
          atom.cpos=v;
          atom.fpos=C2F(a.lattice,atom.cpos);
        }
        //[CO20200130 - number->basis]atom.number=iat;    // reference position for convasp
        atom.basis=iat;     // position in the basis
        atom.ijk(1)=0;atom.ijk(2)=0;atom.ijk(3)=0; // inside the zero cell...
        atom.corigin(1)=0.0;atom.corigin(2)=0.0;atom.corigin(3)=0.0; // inside the zero cell
        atom.coord(1)=0.0;atom.coord(2)=0.0;atom.coord(3)=0.0; // inside the zero cell
        atom.spin=0.0;
        atom.noncoll_spin.clear(); //DX20171205 - non-collinear spin
        atom.type=itype;                // CONVASP_MODE if I want type 0,1,2,3,...
        atom.order_parameter_value=0;
        atom.order_parameter_atom=FALSE;
        atom.partial_occupation_value=1.0;
        atom.partial_occupation_flag=FALSE;
        //	plug_atom=TRUE;   // not used
        // NO ORDER PARAMETER
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " DEBUG [7]" << endl;
        if(a.order_parameter_structure==FALSE) {
          // stmp_tokens.size() = 4 (plus possible comments).
          // Read in the names.
          if(stmp_tokens.size()==4 || (stmp_tokens.size()>=4 && a.isd==0)) {
            atom.name=stmp_tokens.at(id++);
            atom.CleanName();
            atom.CleanSpin();
            atom.name_is_given=TRUE;
          }
          // stmp_tokens.size()=6
          if(a.isd!=0 && stmp_tokens.size()==6) {
            string sdt;
            for(int ic=1;ic<=3;ic++)
              sdt=sdt+stmp_tokens.at(id++);
            atom.sd=sdt;
          }
          // stmp_tokens.size()=7
          if(a.isd==TRUE && stmp_tokens.size()>=7) {
            string sdt;
            for(int ic=1;ic<=3;ic++)
              sdt=sdt+stmp_tokens.at(id++);
            atom.sd=sdt;
            atom.name=stmp_tokens.at(id++);
            atom.CleanName();
            atom.CleanSpin();
            atom.name_is_given=TRUE;
          }
        }
        // ORDER PARAMETER
        if(a.order_parameter_structure==TRUE) {
          if(stmp_tokens.size()!=5 && stmp_tokens.size()!=6) {
            if(a.order_parameter_structure==TRUE) {
              message << "With the order parameter you must specify " << endl;  //CO20190629
              message << "  x y z Name OrderParameter "                     << endl;  //CO20190629
              message << " where x,y,z are the coordinates " << endl; //CO20190629
              message << " Name is the symbol of the atom  " << endl; //CO20190629
              message << " Order parameter is -=none, *=consider, -1,0,1 (integer) values " << endl; //CO20190629
              throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
            }
          }
          if(stmp_tokens.size()==5 || stmp_tokens.size()==6) {
            atom.name=stmp_tokens.at(id);
            atom.CleanName();atom.CleanSpin();
            atom.name_is_given=TRUE;
            id++;
            if(a.order_parameter_structure==TRUE) {
              if(stmp_tokens.at(id)=="-") {
                atom.order_parameter_value=0;
                atom.order_parameter_atom=FALSE;
              } else {
                if(stmp_tokens.at(id)=="*") {
                  atom.order_parameter_value=0;
                  atom.order_parameter_atom=TRUE;
                } else {
                  atom.order_parameter_value=aurostd::string2utype<int>(stmp_tokens.at(id).c_str());
                  atom.order_parameter_atom=TRUE;
                  a.order_parameter_sum+=atom.order_parameter_value;
                }
              }
              id++;
            }
          }
        }
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " DEBUG [8]" << endl;
        // now plug atom into the atomlist
        //      cerr << atom.name << endl;
        //      cerr << atom.cleanname << endl;
        //      cerr << atom.name_is_given << endl;
        a.atoms.push_back(atom);
      } // iat loop
    }
  } //  end VASP POSCAR MODE


  // ----------------------------------------------------------------------
  // QE INPUT
  if(a.iomode==IOQE_AUTO || a.iomode==IOQE_GEOM) {
    // LDEBUG=TRUE;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " QUANTUM ESPRESSO IOQE_AUTO/GEOM" << endl;
    // START FROM CELL
    a.scale=1.0; // standard
    a.neg_scale=FALSE; // standard
    int ibrav_value = 0; //DX20180123 - added ibrav
    int nat_value = 0;   //DX20180123 - added nat
    bool bohr_lat = false;   //DX20180215 - added bohr for lattice
    bool bohr_pos = false;   //DX20180215 - added bohr for positions
    bool alat = false;   //DX20180215 - added alat
    uint iline = 0, jline = 0;
    iline=vinput.size();
    //DX20180123 - added nat - START
    for(uint i=0;i<vinput.size();i++){ 
      if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"NAT=",true)){
        vector<string> comma_tokens; 
        aurostd::string2tokens(aurostd::RemoveWhiteSpaces(vinput[i]),comma_tokens,","); //DX - it is possible to have multiple fields per line
        for(uint j=0;j<comma_tokens.size();j++){
          if(aurostd::substring2bool(aurostd::toupper(comma_tokens.at(j)),"NAT=",true)){
            vector<string> tokens; 
            aurostd::string2tokens(aurostd::RemoveWhiteSpaces(comma_tokens.at(j)),tokens,"=");
            nat_value = aurostd::string2utype<int>(tokens.at(1));
          }
        }
      }                                       
    }
    //DX20180123 - added nat - END
    for(uint i=0;i<vinput.size();i++){
      //DX20180123 [OBSOLETE]  if(aurostd::substring2bool(vinput[i],"CELL_PARAMETERS") && 
      //DX20180123 [OBSOLETE]    aurostd::substring2bool(vinput[i],"angstrom","ANGSTROM")) {iline=i+1;jline=i;}
      if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"CELL_PARAMETERS")){ //DX20170124 - added aurostd::toupper
        if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"ANGSTROM")){
          iline=i+1;jline=i;
        }
        else if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"ALAT")){ //DX20180215 - added alat bool
          iline=i+1;jline=i;                                                     //DX20180215 - added alat bool
          alat = true;                                                           //DX20180215 - added alat bool
        }
        //DX20180123 - added bohr - START
        else if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"BOHR")){
          iline=i+1;jline=i;
          bohr_lat = true;
        }
        //DX20180123 - added bohr - END
      }
      //DX20180123 - added ibrav - START
      else if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"IBRAV=",true)){
        vector<string> comma_tokens; 
        aurostd::string2tokens(aurostd::RemoveWhiteSpaces(vinput[i]),comma_tokens,","); //DX - it is possible to have multiple fields per line
        for(uint j=0;j<comma_tokens.size();j++){
          if(aurostd::substring2bool(aurostd::toupper(comma_tokens.at(j)),"IBRAV=",true)){
            vector<string> tokens; 
            aurostd::string2tokens(aurostd::RemoveWhiteSpaces(comma_tokens.at(j)),tokens,"=");
            ibrav_value = aurostd::string2utype<int>(tokens.at(1));
          }
        }
      }
      //DX20180123 - added ibrav - END
    }
    xvector<double> parameters(6); 
    uint celldm_count = 0;
    bool isabc = false; // distinguish between celldm and a,b,c,cosAB,cosAC,cosBC

    // celldm variant
    for(uint i=0;i<vinput.size();i++){
      if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"CELLDM")){
        vector<string> comma_tokens; 
        aurostd::string2tokens(aurostd::RemoveWhiteSpaces(vinput[i]),comma_tokens,","); //DX - it is possible to have multiple fields per line
        for(uint j=0;j<comma_tokens.size();j++){
          if(aurostd::substring2bool(aurostd::toupper(comma_tokens.at(j)),"CELLDM",true)){ //DX20170124 - added aurostd::toupper
            celldm_count++;
            vector<string> tokens; 
            aurostd::string2tokens(aurostd::RemoveWhiteSpaces(comma_tokens.at(j)),tokens,"=");
            // In case the celldms are not in order, or not all are given; need to check number
            vector<string> parentheses_tokens; 
            aurostd::string2tokens(aurostd::RemoveWhiteSpaces(tokens.at(0)),parentheses_tokens,"(");
            uint index = aurostd::string2utype<uint>(aurostd::RemoveCharacter(parentheses_tokens.at(1),')'));  
            parameters[index] = aurostd::string2utype<double>(tokens.at(1));
          }
        }
      } 
    }
    // a,b,c,cosAB,cosAC,cosBC variant
    if(celldm_count == 0){
      isabc = true;
      for(uint i=0;i<vinput.size();i++){
        if(aurostd::substring2bool(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])),"A=")){
          vector<string> comma_tokens; 
          aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])),comma_tokens,","); //DX - it is possible to have multiple fields per line
          for(uint j=0;j<comma_tokens.size();j++){
            if(aurostd::substring2bool(aurostd::toupper(comma_tokens.at(j)),"A=",true)){ //DX20170124 - added aurostd::toupper
              vector<string> tokens; 
              aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(comma_tokens.at(j))),tokens,"A=");
              if(tokens.size()==1){
                parameters[1] = aurostd::string2utype<double>(tokens.at(0)); //A only
              }
            }
          }
        } 
        if(aurostd::substring2bool(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])),"B=")){
          vector<string> comma_tokens; 
          aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])),comma_tokens,","); //DX - it is possible to have multiple fields per line
          for(uint j=0;j<comma_tokens.size();j++){
            if(aurostd::substring2bool(aurostd::toupper(comma_tokens.at(j)),"B=",true)){ //DX20170124 - added aurostd::toupper
              vector<string> tokens; 
              aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(comma_tokens.at(j))),tokens,"B=");
              if(tokens.size()==1){
                parameters[2] = aurostd::string2utype<double>(tokens.at(0)); //B only
              }
            }
          }
        } 
        if(aurostd::substring2bool(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])),"C=")){
          vector<string> comma_tokens; 
          aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])),comma_tokens,","); //DX - it is possible to have multiple fields per line
          for(uint j=0;j<comma_tokens.size();j++){
            if(aurostd::substring2bool(aurostd::toupper(comma_tokens.at(j)),"C=",true)){ //DX20170124 - added aurostd::toupper
              vector<string> tokens; 
              aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(comma_tokens.at(j))),tokens,"C=");
              if(tokens.size()==1){
                parameters[3] = aurostd::string2utype<double>(tokens.at(0)); //C only
              }
            }
          }
        } 
        if(aurostd::substring2bool(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])),"COSAB=")){
          vector<string> comma_tokens; 
          aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])),comma_tokens,","); //DX - it is possible to have multiple fields per line
          for(uint j=0;j<comma_tokens.size();j++){
            if(aurostd::substring2bool(aurostd::toupper(comma_tokens.at(j)),"COSAB=",true)){ //DX20170124 - added aurostd::toupper
              vector<string> tokens; 
              aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(comma_tokens.at(j))),tokens,"COSAB=");
              if(tokens.size()==1){
                parameters[4] = aurostd::string2utype<double>(tokens.at(0)); // COSAB only
              }
            }
          }
        } 
        if(aurostd::substring2bool(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])),"COSAC=")){
          vector<string> comma_tokens; 
          aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])),comma_tokens,","); //DX - it is possible to have multiple fields per line
          for(uint j=0;j<comma_tokens.size();j++){
            if(aurostd::substring2bool(aurostd::toupper(comma_tokens.at(j)),"COSAC=",true)){ //DX20170124 - added aurostd::toupper
              vector<string> tokens; 
              aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(comma_tokens.at(j))),tokens,"COSAC=");
              if(tokens.size()==1){
                parameters[5] = aurostd::string2utype<double>(tokens.at(0)); // COSAC only
              }
            }
          }
        } 
        if(aurostd::substring2bool(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])),"COSBC=")){
          vector<string> comma_tokens; 
          aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])),comma_tokens,","); //DX - it is possible to have multiple fields per line
          for(uint j=0;j<comma_tokens.size();j++){
            if(aurostd::substring2bool(aurostd::toupper(comma_tokens.at(j)),"COSBC=",true)){ //DX20170124 - added aurostd::toupper
              vector<string> tokens; 
              aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(comma_tokens.at(j))),tokens,"COSBC=");
              if(tokens.size()==1){
                parameters[6] = aurostd::string2utype<double>(tokens.at(0)); // COSBC only
              }
            }
          }
        } 
      }
    }
    //DX20180123 - added ibrav/parameters - START
    if(ibrav_value!=0){
      a.lattice = pflow::QE_ibrav2lattice(ibrav_value, parameters, isabc);
    }
    //DX20180123 - added ibrav/parameters - END
    else {
      if(iline<vinput.size()-3) { // IN a1 / a2 / a3 rows version
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " QUANTUM ESPRESSO FOUND 3 extra lines after, trying a1/a2/a3 on raws" << endl;
        stringstream input_tmp;
        input_tmp.clear();input_tmp.str(vinput.at(iline++));
        input_tmp >> a.lattice(1,1) >> a.lattice(1,2) >> a.lattice(1,3);
        input_tmp.clear();input_tmp.str(vinput.at(iline++));
        input_tmp >> a.lattice(2,1) >> a.lattice(2,2) >> a.lattice(2,3);
        input_tmp.clear();input_tmp.str(vinput.at(iline++));
        input_tmp >> a.lattice(3,1) >> a.lattice(3,2) >> a.lattice(3,3);
        if(bohr_lat){                          //DX20180123 - added bohr 
          a.lattice = a.lattice*bohr2angstrom; //DX20180123 - added bohr
        }                                      //DX20180123 - added bohr
        if(alat){                              //DX20180215 - added alat 
          if(isabc){
            a.lattice = a.lattice*parameters[1]; //DX20180123 - added alat
          }
          else {
            a.lattice = a.lattice*parameters[1]*bohr2angstrom; //DX20180123 - added alat
          }
        }                             //DX20180215 - added alat
        xvector<double> data(6);
        data=Getabc_angles(a.lattice,DEGREES);
        a.a=data[1];a.b=data[2];a.c=data[3];
        a.alpha=data[4];a.beta=data[5];a.gamma=data[6];
      }
    }
    // ---------------------------------------------------------------
    // SOME FOR ABCCAR
    // ...
    // ---------------------------------------------------------------
    // SOME FOR WYCCAR
    // ...
    a.FixLattices(); // Reciprocal/f2c/c2f
    a.kpoints_k1=0;a.kpoints_k2=0;a.kpoints_k3=0;
    a.kpoints_kmax=0;a.kpoints_kppra=0;a.kpoints_kscheme="";
    clear(a.origin);

    // cerr << "iline=" << iline << endl;
    // cerr << "vinput.size()=" << vinput.size() << endl;
    // NOW ADD ATOMS
    for(uint i=0;i<vinput.size();i++) {
      //DX20180124 [OBSOLETE] if(aurostd::substring2bool(vinput[i],"ATOMIC_POSITIONS") && aurostd::substring2bool(vinput[i],"crystal","CRYSTAL"))
      if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"ATOMIC_POSITIONS") && aurostd::substring2bool(aurostd::toupper(vinput[i]),"CRYSTAL")) //DX20180124 - added aurostd::toupper
      { //CO20200106 - patching for auto-indenting
        iline=i+1;a.coord_flag=_COORDS_FRACTIONAL_;
        jline=iline+nat_value; //DX20180123 - added nat
      }
      //DX20180124 [OBSOLETE] if(aurostd::substring2bool(vinput[i],"ATOMIC_POSITIONS") && aurostd::substring2bool(vinput[i],"angstrom","ANGSTROM"))
      if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"ATOMIC_POSITIONS") && aurostd::substring2bool(aurostd::toupper(vinput[i]),"ANGSTROM")) //DX20180124 - added aurostd::toupper
      { //CO20200106 - patching for auto-indenting
        iline=i+1;a.coord_flag=_COORDS_CARTESIAN_;
        jline=iline+nat_value; //DX20180123 - added nat
      }
      //DX20180124 -- added Bohr case - START
      if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"ATOMIC_POSITIONS") && aurostd::substring2bool(aurostd::toupper(vinput[i]),"BOHR")) { //DX20180124 - added aurostd::toupper
        iline=i+1;a.coord_flag=_COORDS_CARTESIAN_;
        jline=iline+nat_value; //DX20180123 - added nat
        bohr_pos = true;
      }
      //DX20180124 -- added Bohr case - END
    }

    for(uint i=iline;i<vinput.size()&&i<jline;i++) {
      _atom atom;                                        // create new atom
      string stmp;
      xvector<double> v(3);clear(v);
      stmp=vinput[i];       stmp=aurostd::RemoveCharacter(vinput[i],'\t');aurostd::StringSubst(stmp,"\t"," ");
      std::vector<string> stmp_tokens;
      aurostd::string2tokens(stmp,stmp_tokens);
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " DEBUG [6b] stmp_tokens.size()=" << stmp_tokens.size() << endl;
      if(stmp_tokens.size()<4) {
        message << "Insufficient number of atom entries in atom=" << i << "" << endl; //CO20190629
        for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
      }  
      int id=0;
      string name;
      // Read in the atom positions.
      name=stmp_tokens.at(id++);
      v(1)=aurostd::string2utype<double>(stmp_tokens.at(id++)); // when not a matter of speed, vectors
      v(2)=aurostd::string2utype<double>(stmp_tokens.at(id++)); // should be indicized with () instead of [],
      v(3)=aurostd::string2utype<double>(stmp_tokens.at(id++)); // b/c () does boundary checks !
      if(a.coord_flag==_COORDS_FRACTIONAL_) {
        atom.fpos=v;
        atom.cpos=F2C(a.lattice,atom.fpos);
      }
      if(a.coord_flag==_COORDS_CARTESIAN_) {
        if(bohr_pos){ //DX20180124 - added bohr case
          atom.cpos=bohr2angstrom*v; //DX20180124 - added bohr case
          atom.fpos=C2F(a.lattice,atom.cpos);
        }
        else { //DX20180124 - added bohr case
          atom.cpos=v;
          atom.fpos=C2F(a.lattice,atom.cpos);
        }
      }
      atom.name=name;
      atom.CleanName();
      atom.CleanSpin();
      atom.name_is_given=TRUE;

      //[CO20200130 - number->basis]atom.number=a.atoms.size();    // reference position for convasp
      atom.basis=a.atoms.size();     // position in the basis
      atom.ijk(1)=0;atom.ijk(2)=0;atom.ijk(3)=0; // inside the zero cell...
      atom.corigin(1)=0.0;atom.corigin(2)=0.0;atom.corigin(3)=0.0; // inside the zero cell
      atom.coord(1)=0.0;atom.coord(2)=0.0;atom.coord(3)=0.0; // inside the zero cell
      atom.spin=0.0;
      atom.noncoll_spin.clear(); //DX20171205 - non-collinear spin
      // FIXED BELOW atom.type=itype;                // CONVASP_MODE if I want type 0,1,2,3,...
      atom.order_parameter_value=0;
      atom.order_parameter_atom=FALSE;
      atom.partial_occupation_value=1.0;
      atom.partial_occupation_flag=FALSE;
      // DONE
      a.AddAtom(atom,true); //CO20230319 - add by species
      // NO PARTIAL OCCUPATION
      a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
    }
    // FIX TITLE
    // FIX ITYPE
    // for(uint i=0;i<a.atoms.size();i++) cerr << "a.atoms.at(i).type=" << a.atoms.at(i).type << endl;
    a.SpeciesPutAlphabetic(); // fight analphabetization
    uint iat=0;
    string _title=""; //CO20200731 - NO_TITLE_GIVEN clear
    for(uint itype=0;itype<a.num_each_type.size();itype++) 
      for(uint j=0;j<(uint) a.num_each_type.at(itype);j++) {
        if(j==0) _title+=a.atoms.at(iat).name+aurostd::utype2string(a.num_each_type.at(itype));
        a.atoms.at(iat++).type=itype;
      }
    if(!_title.empty()){a.title=_title;}  //CO20200731 - safer for empty titles
    // for(uint i=0;i<a.atoms.size();i++) cerr << "a.atoms.at(i).type=" << a.atoms.at(i).type << endl;
    a.partial_occupation_flag=FALSE;
    //    a.xstructure2vasp();
    a.is_vasp4_poscar_format=FALSE;
    a.is_vasp5_poscar_format=FALSE;
    // DONE ?
  } // QE INPUT

  // ----------------------------------------------------------------------
  // ABINIT INPUT (DX20200310)
  if(a.iomode==IOABINIT_AUTO || a.iomode==IOABINIT_GEOM) { // ABINIT
    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ABINIT READER begin" << endl; }

    a.scale=1.0; // standard
    a.neg_scale=FALSE; // standard

    // ----------------------------------------------------------------------
    // get lattice scaling: acell (optional keyword)
    xvector<double> acell; acell(1)=1.0; acell(2)=1.0; acell(3)=1.0;
    for(uint i=0;i<vinput.size();i++){
      if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"ACELL",true)){
        if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ABINIT READER acell line=" << vinput[i] << endl; }
        string acell_line = aurostd::toupper(vinput[i]);
        aurostd::StringSubst(acell_line,"ACELL","");

        // get units first (then remove from temp string for easy parsing)
        string lattice_vec_unit = "Bohr"; // default is Bohr
        if(aurostd::substring2bool(acell_line,"ANGST",true)){ 
          lattice_vec_unit = "Angstrom"; 
          aurostd::StringSubst(acell_line,"ANGSTROM","");
          aurostd::StringSubst(acell_line,"ANGSTR","");
          aurostd::StringSubst(acell_line,"ANGST","");
        }
        if(aurostd::substring2bool(acell_line,"BOHR",true)){
          lattice_vec_unit = "Bohr";
          aurostd::StringSubst(acell_line,"BOHR","");
        }

        // clean-up
        acell_line = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(acell_line);

        // check if explicitly given or uses multiplication
        // explicit : e.g., "acell 1.0 1.0 1.0"
        // multiplication : e.g., "acell 3*1.0"
        bool multiplication_variant = false;
        if(aurostd::substring2bool(acell_line,"*",true)){ multiplication_variant=true; }

        if(multiplication_variant){
          vector<string> tokens;
          aurostd::string2tokens(acell_line,tokens,"*");
          double factor = aurostd::string2utype<double>(aurostd::RemoveWhiteSpaces(tokens[1]));
          acell= factor*acell;
        }
        else{
          vector<string> tokens;
          uint number_tokens = aurostd::string2tokens(acell_line,tokens," ");
          if(number_tokens==3){
            acell(1) = aurostd::string2utype<double>(tokens[0]);
            acell(2) = aurostd::string2utype<double>(tokens[1]);
            acell(3) = aurostd::string2utype<double>(tokens[2]);
          }
          else{
            message << "Unable to parse the acell line; unexpected format. acell_line = " << acell_line << " tokens: " << aurostd::joinWDelimiter(tokens,",");
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_FILE_WRONG_FORMAT_);
          }
        }
        if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ABINIT READER extracted acell = " << acell << endl; }
      }
    }

    // ----------------------------------------------------------------------
    // get lattice (rprim) : given column-major (optional)
    // default: identity (i.e., cube)
    // supported formats:
    //   1) single line 9 fields
    //   2) three lines 3 fields per line
    bool is_lattice_line = false;
    uint lattice_line_count = 0;
    a.lattice = aurostd::identity((double) 0,3); //abinit default //DX20200521 - new identity format
    for(uint i=0;i<vinput.size();i++){
      if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"RPRIM",true)){
        if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ABINIT READER rprim line found = " << vinput[i] << endl; }
        is_lattice_line = true;
      }
      if(is_lattice_line){
        string rprim_line = aurostd::toupper(vinput[i]);
        aurostd::StringSubst(rprim_line,"RPRIM","");
        vector<string> tokens;
        uint number_tokens = aurostd::string2tokens(rprim_line,tokens," ");
        if(number_tokens==9){
          //column-major
          a.lattice(1,1) = aurostd::frac2dbl(tokens[0])*acell(1);
          a.lattice(2,1) = aurostd::frac2dbl(tokens[1])*acell(1);
          a.lattice(3,1) = aurostd::frac2dbl(tokens[2])*acell(1);
          a.lattice(1,2) = aurostd::frac2dbl(tokens[3])*acell(2);
          a.lattice(2,2) = aurostd::frac2dbl(tokens[4])*acell(2);
          a.lattice(3,2) = aurostd::frac2dbl(tokens[5])*acell(2);
          a.lattice(1,3) = aurostd::frac2dbl(tokens[6])*acell(3);
          a.lattice(2,3) = aurostd::frac2dbl(tokens[7])*acell(3);
          a.lattice(3,3) = aurostd::frac2dbl(tokens[8])*acell(3);
          is_lattice_line = false;
          break;
        }
        else if(number_tokens==3){
          lattice_line_count += 1;
          //column-major
          a.lattice(1,lattice_line_count) = aurostd::frac2dbl(tokens[0])*acell(lattice_line_count);
          a.lattice(2,lattice_line_count) = aurostd::frac2dbl(tokens[1])*acell(lattice_line_count);
          a.lattice(3,lattice_line_count) = aurostd::frac2dbl(tokens[2])*acell(lattice_line_count);
          if(lattice_line_count == 3){ is_lattice_line=false; break; }
        }
        else if(number_tokens==0){ continue; }
        else{
          message << "Unable to parse the rprim line; unexpected format.";
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_FILE_WRONG_FORMAT_);
        }
      }
    }
    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ABINIT READER lattice (row-major; AFLOW convention) = " << endl << a.lattice << endl; }

    // ----------------------------------------------------------------------
    // get number of atoms
    // NOTE: sometimes multiple keywords can be in the same line
    // natom : number of atoms in unit cell
    // natrd : number of atoms to read in file (useful for symmetry reduced files, e.g., iatom representations)
    uint number_of_atoms = 0;
    for(uint i=0;i<vinput.size();i++){
      if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"NATOM",true)){
        if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ABINIT READER natom line = " << vinput[i] << endl; }
        string natom_line = aurostd::toupper(vinput[i]);
        vector<string> tokens; 
        aurostd::string2tokens(natom_line,tokens," ");
        // the loop method protects against multiple keywords per line
        for(uint t=0;t<tokens.size();t++){
          if(tokens[t]=="NATOM"){ number_of_atoms = aurostd::string2utype<uint>(tokens[t+1]); break; }
        }
        if(number_of_atoms == 0){ 
          message << "Unable to parse the natom line; unexpected format. vinput[i] = \"" << natom_line << "\".";
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_FILE_WRONG_FORMAT_);
        }
      }
    }
    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ABINIT READER number of atoms; natom = " << number_of_atoms << endl; }

    // natrd (optional)
    bool found_natrd = false;
    uint number_of_atoms_to_read = 0;
    for(uint i=0;i<vinput.size();i++){
      if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"NATRD",true)){
        if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ABINIT READER natrd line = " << vinput[i] << endl; }
        found_natrd = true;
        string natrd_line = aurostd::toupper(vinput[i]);
        vector<string> tokens; 
        aurostd::string2tokens(natrd_line,tokens," ");
        // the loop method protects against multiple keywords per line
        for(uint t=0;t<tokens.size();t++){
          if(tokens[t]=="NATRD"){ number_of_atoms_to_read = aurostd::string2utype<uint>(tokens[t+1]); break; }
        }
        if(number_of_atoms_to_read == 0){ 
          message << "Unable to parse the natrd line; unexpected format. natrd_line = \"" << natrd_line << "\".";
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_FILE_WRONG_FORMAT_);
        }
        if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ABINIT READER number of atoms to read; natrd = " << number_of_atoms_to_read << endl; }
      }
    }
    if(!found_natrd){ number_of_atoms_to_read = number_of_atoms; } // if natrd not specified set to natom

    if(number_of_atoms != number_of_atoms_to_read){
      message << "The natoms != natrd, i.e., unit cell requires atoms (iatoms) to be expanded by symmetry. Functionality not yet supported.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_FILE_WRONG_FORMAT_);
    }

    // ----------------------------------------------------------------------
    // get number of types (optional)
    bool found_ntypat = false;
    uint number_of_atom_types = 0;
    for(uint i=0;i<vinput.size();i++){
      if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"NTYPAT",true)){
        if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ABINIT READER ntypat line = " << vinput[i] << endl; }
        found_ntypat = true;
        string ntypat_line = aurostd::toupper(vinput[i]);
        vector<string> tokens; 
        aurostd::string2tokens(ntypat_line,tokens," ");
        // the loop method protects against multiple keywords per line
        for(uint t=0;t<tokens.size();t++){
          if(tokens[t]=="NTYPAT"){ number_of_atom_types = aurostd::string2utype<uint>(tokens[t+1]); break; }
        }
        if(number_of_atoms_to_read == 0){ 
          message << "Unable to parse the ntypat line; unexpected format. ntypat_line = \"" << ntypat_line << "\".";
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_FILE_WRONG_FORMAT_);
        }
        if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ABINIT READER number of atom types; ntypat = " << number_of_atom_types << endl; }
      }
    }

    // ----------------------------------------------------------------------
    // get atom type index
    // in same order as xred, xcart, or xangst
    vector<uint> atom_types; //follows order of atoms
    for(uint i=0;i<vinput.size();i++){
      if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"TYPAT",true)){
        if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ABINIT READER typat line = " << vinput[i] << endl; }
        string typat_line = aurostd::toupper(vinput[i]);
        vector<string> tokens; 
        aurostd::string2tokens(typat_line,tokens," ");
        // the loop method protects against multiple keywords per line
        bool found_typat = false;
        bool is_field_attribute = false;
        for(uint t=0;t<tokens.size();t++){
          if(tokens[t]=="TYPAT"){ found_typat = true; is_field_attribute = true; continue; }
          if(is_field_attribute){
            if(aurostd::substring2bool(tokens[t],"*")){
              vector<string> sub_tokens; aurostd::string2tokens(tokens[t],sub_tokens,"*");
              uint multiplier = aurostd::string2utype<uint>(sub_tokens[0]);
              uint type_index = aurostd::string2utype<uint>(sub_tokens[1]);
              for(uint m=0;m<multiplier;m++){ atom_types.push_back(type_index); }
            }
            else if(tokens[t][0]>='0' && tokens[t][0]<='9'){ //is digit
              atom_types.push_back(aurostd::string2utype<uint>(tokens[t]));
            }
            else{
              // signals a new keyword in the same line
              is_field_attribute = false;
              break;
            }
          }
        }
        if(!found_typat){ continue; } //found ntypat not typat, try another line
      }
    }
    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ABINIT READER atom type sequence (matches order of atom positions); typat = " << aurostd::joinWDelimiter(atom_types,",") << endl; }

    // check the number of unique types matches ntypat
    vector<uint> unique_types = atom_types;
    std::stable_sort(unique_types.begin(),unique_types.end());
    unique_types.erase(std::unique(unique_types.begin(),unique_types.end()),unique_types.end());
    if(found_ntypat && (unique_types.size() != number_of_atom_types)){
      message << "Number of atom types does not match ntypat variable. typat=" << aurostd::joinWDelimiter(unique_types,",") << " | ntypat=" << number_of_atom_types;
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
    }
    if(atom_types.size() != number_of_atoms_to_read){ 
      message << "Number of atom types does not match the number of atoms to read. typat=" << aurostd::joinWDelimiter(atom_types,",") << " | natom=" << number_of_atoms;
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
    }

    // ----------------------------------------------------------------------
    // get atom type; element (znucl)
    vector<uint> nuclear_charge;
    for(uint i=0;i<vinput.size();i++){
      if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"ZNUCL",true)){
        if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ABINIT READER znucl line = " << vinput[i] << endl; }
        string znucl_line = aurostd::toupper(vinput[i]);
        vector<string> tokens; 
        aurostd::string2tokens(znucl_line,tokens," ");
        // the loop method protects against multiple keywords per line
        bool is_field_attribute = false;
        for(uint t=0;t<tokens.size();t++){
          if(tokens[t]=="ZNUCL"){ is_field_attribute = true; continue; }
          if(is_field_attribute){
            if(tokens[t][0]>='0' && tokens[t][0]<='9'){ //is digit
              nuclear_charge.push_back(aurostd::string2utype<uint>(tokens[t]));
            }
            else{
              // signals a new keyword in the same line
              is_field_attribute = false;
              break;
            }
          }
        }
      }
    }
    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ABINIT READER nuclear charge of each atom type (order corresponds to typat index); znucl = " << aurostd::joinWDelimiter(nuclear_charge,",") << endl; }

    // ----------------------------------------------------------------------
    // normally a mandatory keyword, but the AFLOW-ABINIT writer has not been printing this keyword
    // for backwards compatability, we will throw a warning (for now)
    if(nuclear_charge.size()==0){ 
      message << "The atom elements/types (znucl) are not specified. Using fictious atoms.";
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,std::cerr,_LOGGER_WARNING_);
    }

    // ----------------------------------------------------------------------
    // get atom positions
    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ABINIT READER reading atom positions" << endl; }
    deque<_atom> atoms_temp;
    bool is_atom_line = false;
    bool is_Bohr_units = false;
    uint atom_line_count = 0;
    for(uint i=0;i<vinput.size();i++){
      if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"XRED",true)){
        if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ABINIT READER xred (fractional) line = " << vinput[i] << endl; }
        a.coord_flag = _COORDS_FRACTIONAL_;
        is_atom_line = true;
      }
      else if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"XCART",true)){
        if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ABINIT READER xcart (Cartesian, unit=Bohr) line = " << vinput[i] << endl; }
        a.coord_flag = _COORDS_CARTESIAN_;
        is_Bohr_units = true;
        is_atom_line = true;
      }
      else if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"XANGST",true)){
        if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ABINIT READER xangst (Cartesian, unit=Angstrom) line = " << vinput[i] << endl; }
        a.coord_flag = _COORDS_CARTESIAN_;
        is_atom_line = true;
      }
      if(is_atom_line){
        string atom_line = aurostd::toupper(vinput[i]);
        aurostd::StringSubst(atom_line,"XRED",""); //for easier parsing
        aurostd::StringSubst(atom_line,"XCART",""); //for easier parsing 
        aurostd::StringSubst(atom_line,"XANGST",""); //for easier parsing
        vector<string> tokens;
        uint number_tokens = aurostd::string2tokens(atom_line,tokens," ");
        if(number_tokens==3){
          atom_line_count += 1;
          xvector<double> coordinate;
          coordinate(1) = aurostd::frac2dbl(tokens[0]);
          coordinate(2) = aurostd::frac2dbl(tokens[1]);
          coordinate(3) = aurostd::frac2dbl(tokens[2]);

          _atom atom;
          if(a.coord_flag == _COORDS_FRACTIONAL_){ 
            atom.fpos = coordinate; 
            atom.cpos = F2C(a.lattice,atom.fpos); 
          }
          else if(a.coord_flag == _COORDS_CARTESIAN_){ 
            atom.cpos = coordinate; 
            if(is_Bohr_units){ atom.cpos *= bohr2angstrom; } // AFLOW expects Angstroms
            atom.fpos = C2F(a.lattice,atom.cpos); 
          }
          atoms_temp.push_back(atom);
        }
        else if(number_tokens==(3*number_of_atoms_to_read)){
          for(uint t=0;t<tokens.size();t+=3){
            atom_line_count += 1;
            xvector<double> coordinate;
            coordinate(1) = aurostd::frac2dbl(tokens[t]);
            coordinate(2) = aurostd::frac2dbl(tokens[t+1]);
            coordinate(3) = aurostd::frac2dbl(tokens[t+2]);

            _atom atom;
            if(a.coord_flag == _COORDS_FRACTIONAL_){ 
              atom.fpos = coordinate; 
              atom.cpos = F2C(a.lattice,atom.fpos); 
            }
            else if(a.coord_flag == _COORDS_CARTESIAN_){ 
              atom.cpos = coordinate; 
              if(is_Bohr_units){ atom.cpos *= bohr2angstrom; } // AFLOW expects Angstroms
              atom.fpos = C2F(a.lattice,atom.cpos); 
            }
            atoms_temp.push_back(atom);
          }
        }
        else if(number_tokens==0){ continue; }
        else{
          message << "Unable to parse the atom line; unexpected format. vinput[i] = " << atom_line;
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_FILE_WRONG_FORMAT_);
        }
        if(atom_line_count == number_of_atoms_to_read){ break;}
      }
    }

    if(LDEBUG){ 
      for(uint iat=0;iat<atoms_temp.size();iat++){
        cerr << __AFLOW_FUNC__ << " ABINIT READER atom position [" << iat << "] = " << atoms_temp[iat] << endl;
      }
    }

    // ----------------------------------------------------------------------
    // ensure correct number of atoms
    if(atoms_temp.size() != number_of_atoms_to_read){ 
      message << "The number of atoms is not commensurate: number of atom positions found = " << atoms_temp.size() << " vs natom/natrd = " << number_of_atoms_to_read;
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,std::cerr,_RUNTIME_ERROR_);
    }

    // ----------------------------------------------------------------------
    // add atoms to xstructure 
    // add name/type info as well
    for(uint i=0;i<atoms_temp.size();i++){
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " ABINIT READER finding element/type for atom [" << a << "]" << endl; }
      xelement::xelement element;
      // use nuclear charge associated atom index
      if(nuclear_charge.size()){ element = xelement::xelement(nuclear_charge[atom_types[i]-1]); } //i-1 since index doesn't start at zero
      // use arbitrary element names (in order of increasing Z)
      else{ element = xelement::xelement(atom_types[i]); }
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " ABINIT READER element/type for atom [" << a << "] = " << element.symbol << endl; }
      atoms_temp[i].name=element.symbol;
      atoms_temp[i].CleanName();
      atoms_temp[i].CleanSpin();
      atoms_temp[i].name_is_given=TRUE;

      //[CO20200130 - number->basis]atoms_temp[i].number=atoms_temp.size();    // reference position for convasp
      atoms_temp[i].basis=atoms_temp.size();     // position in the basis
      atoms_temp[i].ijk(1)=0;atoms_temp[i].ijk(2)=0;atoms_temp[i].ijk(3)=0; // inside the zero cell...
      atoms_temp[i].corigin(1)=0.0;atoms_temp[i].corigin(2)=0.0;atoms_temp[i].corigin(3)=0.0; // inside the zero cell
      atoms_temp[i].coord(1)=0.0;atoms_temp[i].coord(2)=0.0;atoms_temp[i].coord(3)=0.0; // inside the zero cell
      atoms_temp[i].spin=0.0;
      atoms_temp[i].noncoll_spin.clear(); //DX20171205 - non-collinear spin
      // FIXED BELOW atom.type=itype;                // CONVASP_MODE if I want type 0,1,2,3,...
      atoms_temp[i].order_parameter_value=0;
      atoms_temp[i].order_parameter_atom=FALSE;
      atoms_temp[i].partial_occupation_value=1.0;
      atoms_temp[i].partial_occupation_flag=FALSE;
      // DONE
      a.AddAtom(atoms_temp[i],true);  //CO20230319 - add by species
      // NO PARTIAL OCCUPATION
      a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);

    }
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " ABINIT READER fixing atom information (alphabetize, make basis, set species, etc.)" << endl; }
    a.SpeciesPutAlphabetic();
    std::stable_sort(a.atoms.begin(),a.atoms.end(),sortAtomsNames);
    a.MakeBasis();
    a.MakeTypes(); //DX20190508 - otherwise types are not created
    // add system name to title
    a.buildGenericTitle(); //DX20210211
    //DX20210211 [OBSOLETE] uint iat=0;
    //DX20210211 [OBSOLETE - abinit does not generally have a title line] a.title += " "; // add space between existing title
    //DX20210211 [OBSOLETE] for(uint itype=0;itype<a.num_each_type.size();itype++)
    //DX20210211 [OBSOLETE]   for(uint j=0;j<(uint) a.num_each_type.at(itype);j++) {
    //DX20210211 [OBSOLETE]     if(j==0) a.title+=a.atoms.at(iat).name+aurostd::utype2string(a.num_each_type.at(itype));
    //DX20210211 [OBSOLETE]     a.atoms.at(iat++).type=itype;
    //DX20210211 [OBSOLETE]   }
    a.partial_occupation_flag=FALSE;
    a.is_vasp4_poscar_format=FALSE;
    a.is_vasp5_poscar_format=FALSE;

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " ABINIT READER - Finished" << endl; }

  }

  // ----------------------------------------------------------------------
  // ELK INPUT - START (DX20200310)
  // based on information from http://elk.sourceforge.net/elk.pdf
  if(a.iomode==IOELK_AUTO || a.iomode==IOELK_GEOM) { // ELK
    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ELK READER begin" << endl; }

    a.scale=1.0; // standard
    a.neg_scale=FALSE; // standard
    a.coord_flag = _COORDS_FRACTIONAL_; // always for ELK

    // ----------------------------------------------------------------------
    // get isotropic lattice scaling; two types
    //   1) scale : isotropic scaling
    //   2) scale1/2/3 : anisotropic scaling
    double isotropic_scaling = 1.0;
    xvector<double> anisotropic_scaling;
    anisotropic_scaling(1) = 1.0; anisotropic_scaling(2) = 1.0; anisotropic_scaling(3) = 1.0;

    for(uint i=0;i<vinput.size();i++){
      if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"SCALE",true)){

        string scaling_title_line = aurostd::toupper(vinput[i]);
        // ----------------------------------------------------------------------
        // anisotropic
        if(aurostd::substring2bool(scaling_title_line,"SCALE1",true) ||
            aurostd::substring2bool(scaling_title_line,"SCALE2",true) ||
            aurostd::substring2bool(scaling_title_line,"SCALE3",true)){
          if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ELK READER isotropic scaling line found = " << vinput[i] << endl; }
          uint vector_index = aurostd::string2utype<uint>(aurostd::RemoveWhiteSpacesFromTheFrontAndBack(aurostd::StringSubst(scaling_title_line,"SCALE","")));
          string scaling_line = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vinput[i+1]); // scaling value is in the next line
          anisotropic_scaling(vector_index) = aurostd::string2utype<double>(scaling_line);

          if(anisotropic_scaling(vector_index)<_ZERO_TOL_){
            message << "Anisotropic scaling factor is zero or negative for " << vector_index << " component;; unable to parse the scale line; unexpected format. scaling_line = " << scaling_line;
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_FILE_WRONG_FORMAT_);
          }
          if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ELK READER anisotropic scaling for vector " << vector_index << " = " << anisotropic_scaling(vector_index) << endl; }
        }
        // ----------------------------------------------------------------------
        // isotropic
        else{
          if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ELK READER isotropic scaling line found = " << vinput[i] << endl; }
          string scaling_line = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vinput[i+1]); // scaling is the next line
          isotropic_scaling = aurostd::string2utype<double>(scaling_line);
          if(isotropic_scaling<_ZERO_TOL_){
            message << "Isotropic scaling factor is zero or negative; unable to parse the scale line; unexpected format. scaling_line = " << scaling_line;
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_FILE_WRONG_FORMAT_);
          }
          if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ELK READER isotropic scaling = " << isotropic_scaling << endl; }
        }
      }
    }

    // ----------------------------------------------------------------------
    // get lattice (avec) : row-major
    bool is_lattice_line = false;
    uint lattice_line_count = 0;
    for(uint i=0;i<vinput.size();i++){
      if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"AVEC",true)){
        if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ELK READER avec line found = " << vinput[i] << endl; }
        is_lattice_line = true;
      }
      if(is_lattice_line){
        string avec_line = aurostd::toupper(vinput[i]);
        aurostd::StringSubst(avec_line,"AVEC","");
        vector<string> tokens;
        uint number_tokens = aurostd::string2tokens(avec_line,tokens," ");
        if(number_tokens==9){
          //row-major
          a.lattice(1,1) = aurostd::frac2dbl(tokens[0])*isotropic_scaling*anisotropic_scaling(1);
          a.lattice(1,2) = aurostd::frac2dbl(tokens[1])*isotropic_scaling*anisotropic_scaling(1);
          a.lattice(1,3) = aurostd::frac2dbl(tokens[2])*isotropic_scaling*anisotropic_scaling(1);
          a.lattice(2,1) = aurostd::frac2dbl(tokens[3])*isotropic_scaling*anisotropic_scaling(2);
          a.lattice(2,2) = aurostd::frac2dbl(tokens[4])*isotropic_scaling*anisotropic_scaling(2);
          a.lattice(2,3) = aurostd::frac2dbl(tokens[5])*isotropic_scaling*anisotropic_scaling(2);
          a.lattice(3,1) = aurostd::frac2dbl(tokens[6])*isotropic_scaling*anisotropic_scaling(3);
          a.lattice(3,2) = aurostd::frac2dbl(tokens[7])*isotropic_scaling*anisotropic_scaling(3);
          a.lattice(3,3) = aurostd::frac2dbl(tokens[8])*isotropic_scaling*anisotropic_scaling(3);
          is_lattice_line = false;
          break;
        }
        else if(number_tokens==3){
          lattice_line_count += 1;
          //row-major
          a.lattice(lattice_line_count,1) = aurostd::frac2dbl(tokens[0])*isotropic_scaling*anisotropic_scaling(lattice_line_count);
          a.lattice(lattice_line_count,2) = aurostd::frac2dbl(tokens[1])*isotropic_scaling*anisotropic_scaling(lattice_line_count);
          a.lattice(lattice_line_count,3) = aurostd::frac2dbl(tokens[2])*isotropic_scaling*anisotropic_scaling(lattice_line_count);
          if(lattice_line_count == 3){ is_lattice_line=false; break; }
        }
        else if(number_tokens==0){ continue; }
        else{
          message << "Unable to parse the avec line; unexpected format.";
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_FILE_WRONG_FORMAT_);
        }
      }
    }
    a.lattice = a.lattice*bohr2angstrom; // convert to Angstroms
    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ELK READER lattice (row-major; AFLOW convention) = " << endl << a.lattice << endl; }

    // ----------------------------------------------------------------------
    // get number of atom types (atoms)
    double number_of_atom_types = 0;

    for(uint i=0;i<vinput.size();i++){
      if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"ATOMS",true) &&
          !aurostd::substring2bool(aurostd::toupper(vinput[i]),"NATOMS",true)){
        vector<string> tokens;
        aurostd::string2tokens(vinput[i+1],tokens," "); // scaling value is in the next line
        number_of_atom_types = aurostd::string2utype<uint>(aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[0]));

        if(number_of_atom_types == 0){
          message << "Unable to parse the atoms line; unexpected format. atoms_line = \"" << tokens[0] << "\".";
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_FILE_WRONG_FORMAT_);
        }
        if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ELK READER number of atom types atoms = " << number_of_atom_types << endl; }
        break;
      }
    }

    // ----------------------------------------------------------------------
    // get atoms (elements, positions, etc.)
    uint number_of_element_dot_in_files = 0;
    vector<uint> num_each_type;

    for(uint i=0;i<vinput.size();i++){
      if(aurostd::substring2bool(aurostd::toupper(vinput[i]),".IN",true)){
        number_of_element_dot_in_files++;
        if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ELK READER \"ELEMENT.in\" file line found = " << vinput[i] << endl; }

        vector<string> tokens;
        // ----------------------------------------------------------------------
        // get element
        aurostd::string2tokens(vinput[i],tokens,".");
        string element_symbol = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[0]); //clean //DX20210409 - remove whitespace before as well
        element_symbol = aurostd::RemoveCharacterFromTheFrontAndBack(element_symbol,'\''); //clean
        element_symbol = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(element_symbol); //clean
        tokens.clear();
        if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ELK READER element extracted = " << element_symbol << endl; }

        // ----------------------------------------------------------------------
        // get number of that element (+1 line down)
        aurostd::string2tokens(vinput[i+1],tokens," "); // note: +1 line down
        uint num_of_this_type = aurostd::string2utype<uint>(aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[0])); //clean
        num_each_type.push_back(num_of_this_type);
        if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ELK READER number of " << element_symbol << " atoms = " << num_of_this_type << endl; }
        tokens.clear();

        // ----------------------------------------------------------------------
        // get subsequent atom positions
        for(uint iat=1;iat<=num_of_this_type;iat++){
          uint line_index = i+1+iat; // keep track of line index
          _atom atom;

          string atom_line = vinput[line_index];
          uint number_of_tokens = aurostd::string2tokens(atom_line,tokens," ");

          // first three fields are atom positions (fractional)
          if(number_of_tokens >= 3){
            atom.fpos(1) = aurostd::string2utype<double>(tokens[0]);
            atom.fpos(2) = aurostd::string2utype<double>(tokens[1]);
            atom.fpos(3) = aurostd::string2utype<double>(tokens[2]);

            if(number_of_tokens == 6){
              // subsequent are external magnetic field in Cartesian coords
              // not stored beyond here (for now)
              xvector<double> magnetic_field;
              magnetic_field(1) = aurostd::string2utype<double>(tokens[3]);
              magnetic_field(2) = aurostd::string2utype<double>(tokens[4]);
              magnetic_field(3) = aurostd::string2utype<double>(tokens[5]);
              atom.noncoll_spin_is_given = TRUE; //DX20210409
              atom.noncoll_spin = magnetic_field; //DX20210409
              if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ELK READER magnetic field/spin found : " << atom.noncoll_spin << endl; }
            }
          }
          // F2C
          atom.cpos = F2C(a.lattice,atom.fpos);

          atom.name=element_symbol;
          atom.CleanName();
          //DX20210409 [OBSOLETE] atom.CleanSpin();
          atom.name_is_given=TRUE;

          //[CO20200130 - number->basis]atom.number=a.atoms.size();    // reference position for convasp
          atom.basis=a.atoms.size();     // position in the basis
          atom.ijk(1)=0;atom.ijk(2)=0;atom.ijk(3)=0; // inside the zero cell...
          atom.corigin(1)=0.0;atom.corigin(2)=0.0;atom.corigin(3)=0.0; // inside the zero cell
          atom.coord(1)=0.0;atom.coord(2)=0.0;atom.coord(3)=0.0; // inside the zero cell
          //DX20210409 [OBSOLETE] atom.spin=0.0;
          //DX20210409 [OBSOLETE] atom.noncoll_spin.clear();
          atom.order_parameter_value=0;
          atom.order_parameter_atom=FALSE;
          atom.partial_occupation_value=1.0;
          atom.partial_occupation_flag=FALSE;

          if(LDEBUG){ cerr << __AFLOW_FUNC__ << " ELK READER atom added = " << atom << endl; }
          // DONE
          a.AddAtom(atom,true); //CO20230319 - add by species
          // NO PARTIAL OCCUPATION
          a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
        }
        i += num_of_this_type+1; // no need to read over the same lines
      }
    }

    // ----------------------------------------------------------------------
    // check number of element types is consistent
    if(number_of_atom_types == 0 || number_of_atom_types != number_of_element_dot_in_files){
      message << "The number of atom types do not match the number of \"ELEMENT.in\" files: atom = " << number_of_atom_types << " vs # .in files = " << number_of_element_dot_in_files;
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_FILE_WRONG_FORMAT_);
    }
    // ----------------------------------------------------------------------
    // check number of total atoms is consistent
    uint sum_each_type = aurostd::sum(num_each_type);
    if(sum_each_type != a.atoms.size()){
      message << "The total number of atoms does not match the sum of each atom type: a.atoms.size() = " << a.atoms.size() << " vs sum(num_each_type) = " << sum_each_type;;
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_FILE_WRONG_FORMAT_);
    }

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " ELK READER fixing atom information (alphabetize, make basis, set species, etc.)" << endl; }
    a.SpeciesPutAlphabetic();
    std::stable_sort(a.atoms.begin(),a.atoms.end(),sortAtomsNames);
    a.MakeBasis();
    a.MakeTypes();
    // add system name to title
    a.buildGenericTitle(); //DX20210211
    //DX20210211 [OBSOLETE] uint iat=0;
    //DX20210211 [OBSOLETE - Elk does not generally have a title line] a.title += " "; // add space between existing title
    //DX20210211 [OBSOLETE] for(uint itype=0;itype<a.num_each_type.size();itype++)
    //DX20210211 [OBSOLETE]   for(uint j=0;j<(uint) a.num_each_type.at(itype);j++) {
    //DX20210211 [OBSOLETE]     if(j==0) a.title+=a.atoms.at(iat).name+aurostd::utype2string(a.num_each_type.at(itype));
    //DX20210211 [OBSOLETE]    a.atoms.at(iat++).type=itype;
    //DX20210211 [OBSOLETE]  }
    // NO PARTIAL OCCUPATION
    a.partial_occupation_flag=FALSE;
    a.is_vasp4_poscar_format=FALSE;
    a.is_vasp5_poscar_format=FALSE;

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " ELK READER - Finished" << endl; }
  }

  // ----------------------------------------------------------------------
  // CIF INPUT
  if(a.iomode==IOCIF) { // CIF
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CIF [0]" << endl;
    a.scale=1.0; 
    a.neg_scale=FALSE; 
    a.lattice=aurostd::eye<double>(3,3); //CO20190520

    a.spacegroupnumber=0;
    a.spacegroupnumberoption=0;
    // get space group
    for(uint i=0;i<vinput.size();i++) {
      if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"_SYMMETRY_INT_TABLES_NUMBER")){ // converted to upper to be case insensitive
        vector<string> tokens; 
        aurostd::string2tokens(aurostd::RemoveWhiteSpaces(aurostd::toupper(vinput[i])),tokens,"_SYMMETRY_INT_TABLES_NUMBER"); // converted to upper to be case insensitive
        if(tokens.size()==1){
          a.spacegroupnumber = aurostd::string2utype<uint>(tokens.at(0));
        }
      }
      else if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"_SPACE_GROUP_IT_NUMBER")){ // converted to upper to be case insensitive
        vector<string> tokens; 
        aurostd::string2tokens(aurostd::RemoveWhiteSpaces(aurostd::toupper(vinput[i])),tokens,"_SPACE_GROUP_IT_NUMBER"); // converted to upper to be case insensitive
        if(tokens.size()==1){
          a.spacegroupnumber = aurostd::string2utype<uint>(tokens.at(0));
        }
      }
      //DX20190708 - added another space group variant - START
      else if(aurostd::substring2bool(aurostd::toupper(vinput.at(i)),"_SYMMETRY_SPACE_GROUP_NAME_H-M")){ // converted to upper to be case insensitive
        vector<string> tokens; 
        aurostd::string2tokens(vinput.at(i),tokens);
        if(aurostd::toupper(tokens.at(0))=="_SYMMETRY_SPACE_GROUP_NAME_H-M"){
          tokens.erase(tokens.begin());
          string spacegroupsymbol = aurostd::joinWDelimiter(tokens,"");
          spacegroupsymbol = aurostd::RemoveCharacterFromTheFrontAndBack(spacegroupsymbol,'\''); //clean
          spacegroupsymbol = aurostd::RemoveCharacterFromTheFrontAndBack(spacegroupsymbol,'\"'); //clean
          try{ 
            a.spacegroupnumber = GetSpaceGroupNumber(spacegroupsymbol); 
          } //DX20191029 - added try/catch sequence
          catch(aurostd::xerror& re){ 
            if(LDEBUG){ message << "Cannot determine space group setting from the Hermann-Mauguin symbol; non-standard setting."; pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, std::cerr, _LOGGER_WARNING_); }
          } //DX20191029 - added try/catch sequence
        }
      }
      //DX20190708 - added another space group variant - END
    }
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CIF [1]" << endl;
    // get space group setting
    string spacegroup_Hall="";
    for(uint i=0;i<vinput.size();i++) {
      if(aurostd::substring2bool(aurostd::toupper(vinput[i]),"_SYMMETRY_SPACE_GROUP_NAME_HALL")){ // converted to upper to be case insensitive //DX20200521 - missing "_SYMMETRY_" prefix
        vector<string> tokens; 
        aurostd::string2tokens(vinput.at(i),tokens); //DX20200521 - this line was missing
        //DX20190708 - fix Hall reader - START
        if(aurostd::toupper(tokens.at(0))=="_SYMMETRY_SPACE_GROUP_NAME_HALL"){ //DX20200521 - changed H-M to HALL
          tokens.erase(tokens.begin());
          spacegroup_Hall = aurostd::joinWDelimiter(tokens," "); //need a space here for Hall designation
          spacegroup_Hall = aurostd::RemoveCharacterFromTheFrontAndBack(spacegroup_Hall,'\''); //clean
          spacegroup_Hall = aurostd::RemoveCharacterFromTheFrontAndBack(spacegroup_Hall,'\"'); //clean
        }
        //DX20190708 - fix Hall reader - END
      }
    }
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CIF [2]" << endl;
    //DX20191029 - check if space group number is found - START
    if(a.spacegroupnumber==0){
      message << "Either space group number was not given or it was given in a non-standard setting.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__,message,_VALUE_ERROR_);
    }
    //DX20191029 - check if space group number is found - END
    //ME20220124 - Read symmetry operations without consistency checks first
    bool found_setting = false;
    vector<string> spacegroup_symop_xyz;
    bool found_symops=FALSE;
    bool found_symop_id=FALSE; //DX20190708
    for(uint i=0;i<vinput.size();i++) {
      if(aurostd::substring2bool(vinput[i],"_space_group_symop_operation_xyz") || aurostd::substring2bool(vinput[i],"_symmetry_equiv_pos_as_xyz")){
        found_symops=TRUE;
      }
      else if(aurostd::substring2bool(vinput[i],"_space_group_symop_id") || aurostd::substring2bool(vinput[i],"_symmetry_equiv_pos_site_id")){ //DX20190708
        found_symop_id=TRUE;
      }
      else if(found_symops){
        if (vinput[i].find("loop_") != string::npos) break; // End of symops
        vector<string> tokens;
        aurostd::string2tokens(vinput[i],tokens," ");
        //DX20181210 - account for many formats (i.e., x,y,z or 'x, y, z') - START
        if (tokens.size() == 0) continue;
        if (found_symop_id && !isdigit(tokens[0][0])) break;  // End of symops
        if (found_symop_id){ tokens.erase(tokens.begin()); } //erase symop index, not needed //DX20190708 - enclose in if-statement
        string symop = aurostd::joinWDelimiter(tokens,"");
        symop = aurostd::RemoveCharacter(symop,'\''); // remove '
        symop = aurostd::RemoveCharacter(symop,'\"'); // remove "
        symop = aurostd::RemoveWhiteSpaces(symop); // remove spaces
        try {
          symop = SYM::reorderWyckoffPosition(symop); //DX20190708 - standardize order of equation (variable first, then number)
        } catch (aurostd::xerror& e) {
          break;  // End of symops
        }
        spacegroup_symop_xyz.push_back(symop);
        //DX20181210 - account for many formats (i.e., x,y,z or 'x, y, z') - END
      }
    }
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CIF [3]" << endl;
    for(uint setting_number=1;setting_number<=2;setting_number++){
      string setting_string = aurostd::utype2string<uint>(setting_number);
      int general_wyckoff_multiplicity=0; // general Wyckoff position multiplicity
      vector<string> general_wyckoff_position; // general Wyckoff position equations
      // get general Wyckoff multiplicity and position saved in aflow
      SYM::getGeneralWyckoffMultiplicityAndPosition(a.spacegroupnumber, setting_string, general_wyckoff_multiplicity, general_wyckoff_position);
      // ME20220124 - moved up
      //      SYM::initsgs(setting_string);
      //      using SYM::gl_sgs;
      //      cerr << "find the spacegroupstring" << endl;
      //      string spacegroupstring = gl_sgs[a.spacegroupnumber - 1];
      //      cerr << "spacegroupstring: " << spacegroupstring << endl;
      //      vector<int> wyckoff_multiplicities = SYM::get_multiplicities(spacegroupstring);
      //      cerr << "wyckoff_multiplicity: " << wyckoff_multiplicities[1] << endl;
      //      // get symops from cif
      //      vector<string> spacegroup_symop_xyz;  
      //      int multiplicity_count=0;
      //      bool found_symops=FALSE;
      //      cerr << "wyckoff_multiplicities[1]: " << wyckoff_multiplicities[1] << endl;
      //      for(uint i=0;i<vinput.size();i++) {
      //	if(aurostd::substring2bool(vinput[i],"_space_group_symop_operation_xyz") || aurostd::substring2bool(vinput[i],"_symmetry_equiv_pos_as_xyz")){ // _space_group_symop_operation_xyz supersedes all
      //	  found_symops=TRUE;
      //	}
      //	else if(found_symops && multiplicity_count<wyckoff_multiplicities[1]){
      //	  multiplicity_count+=1;
      //	  vector<string> tokens; 
      //	  aurostd::string2tokens(vinput[i],tokens," ");
      //	  if(tokens.size()==2){
      //	    spacegroup_symop_xyz.push_back(tokens[1]);
      //	  }
      //	
      //	}
      //      }
      //      // check symops against general wyckoff position to determine setting
      //      //uint setting_number = 1; //default is first setting
      //      vector<string> general_wyckoff_position = SYM::findGeneralWyckoffPosition(spacegroupstring, wyckoff_multiplicities[1]);
      //      cerr << "general position" << endl;
      //      print(general_wyckoff_position);
      //      cerr << "====================================" << endl;
      //      cerr << "cif symop" << endl;
      //      print(spacegroup_symop_xyz);
      //      cerr << "general_wyckoff_position.size(): " << general_wyckoff_position.size() << endl;
      //      cerr << "spacegroup_symop_xyz.size(): " << spacegroup_symop_xyz.size() << endl;
      // check equations in cif
      //vector<string> spacegroup_symop_xyz;
      //int multiplicity_count=0;
      //bool found_symops=FALSE;
      //bool found_symop_id=FALSE; //DX20190708
      //for(uint i=0;i<vinput.size();i++) {
      //  if(aurostd::substring2bool(vinput[i],"_space_group_symop_operation_xyz") || aurostd::substring2bool(vinput[i],"_symmetry_equiv_pos_as_xyz")){
      //    found_symops=TRUE;
      //  }
      //  else if(aurostd::substring2bool(vinput[i],"_space_group_symop_id") || aurostd::substring2bool(vinput[i],"_symmetry_equiv_pos_site_id")){ //DX20190708
      //    found_symop_id=TRUE;
      //  }
      //  else if(found_symops && multiplicity_count<general_wyckoff_multiplicity){
      //    multiplicity_count+=1;
      //    vector<string> tokens;
      //    aurostd::string2tokens(vinput[i],tokens," ");
      //    //DX20181210 - account for many formats (i.e., x,y,z or 'x, y, z') - START
      //    if(found_symop_id){ tokens.erase(tokens.begin()); } //erase symop index, not needed //DX20190708 - enclose in if-statement
      //    string symop = aurostd::joinWDelimiter(tokens,"");
      //    symop = aurostd::RemoveCharacter(symop,'\''); // remove '
      //    symop = aurostd::RemoveCharacter(symop,'\"'); // remove "
      //    symop = aurostd::RemoveWhiteSpaces(symop); // remove spaces
      //    symop = SYM::reorderWyckoffPosition(symop); //DX20190708 - standardize order of equation (variable first, then number)
      //    spacegroup_symop_xyz.push_back(symop);
      //    //if(tokens.size()==2){
      //    //  spacegroup_symop_xyz.push_back(tokens[1]);
      //    //}
      //    //DX20181210 - account for many formats (i.e., x,y,z or 'x, y, z') - END
      //  }
      //}
      // compare cif and aflow's general position
      uint match_count=0;
      if(general_wyckoff_position.size()==spacegroup_symop_xyz.size()){
        for(uint i=0;i<general_wyckoff_position.size();i++){
          bool match_found = false;
          for(uint j=0;j<spacegroup_symop_xyz.size();j++){
            if(general_wyckoff_position[i]==spacegroup_symop_xyz[j]){
              match_found = true;
              match_count+=1;
              break;
            }
          }
          if(!match_found){
            if(LDEBUG) {
              cerr << "WARNING: Could not match " << i << ": " << general_wyckoff_position[i] << " to any symop in CIF. Trying other setting (if exists)." << endl;
            }
            break;
          }
        }
        if(match_count!=general_wyckoff_position.size()){
          //oss << "ERROR - xstructure::operator>>: Symmetry operations do not match between input operations and space group number/option." << endl; 
          //print(general_wyckoff_position);
          //print(spacegroup_symop_xyz);
          //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_);
          continue; //try a different setting
        }
        else {
          a.spacegroupnumberoption = setting_number;
          a.spacegroupoption = setting_string; //DX20191029
          found_setting = true;
          break; //found setting
        }
      }
    }
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CIF [4]" << endl;
    if(!found_setting){
      // ME20220124 - Changed to warning
      message << "Symmetry operations do not match between input operations and space group number/option.";  //CO20190629
      message << " Building structure using symmetry operations in CIF file with space group P1.";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, std::cerr, _LOGGER_WARNING_);
      a.spacegroupnumber = 1;
      if (LDEBUG) {
        for(uint i=0;i<vinput.size();i++) std::cerr << vinput[i] << endl;  //CO20190629 // ME20220124 - moved to LDEBUG because outputting the entire CIF makes the error message unreadable
      }
      //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
    }
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CIF [5]" << endl;
    // get lattice
    for(uint i=0;i<vinput.size();i++) {
      vector<string> tokens;
      aurostd::string2tokens(vinput[i],tokens," ");
      if(aurostd::substring2bool(vinput[i],"_cell_length_a") && tokens.size()==2){
        a.a = aurostd::string2utype<double>(tokens[1]); 
      }
      else if(aurostd::substring2bool(vinput[i],"_cell_length_b") && tokens.size()==2){
        a.b = aurostd::string2utype<double>(tokens[1]); 
      }
      else if(aurostd::substring2bool(vinput[i],"_cell_length_c") && tokens.size()==2){
        a.c = aurostd::string2utype<double>(tokens[1]); 
      }
      else if(aurostd::substring2bool(vinput[i],"_cell_angle_alpha") && tokens.size()==2){
        a.alpha = aurostd::string2utype<double>(tokens[1]); 
      }
      else if(aurostd::substring2bool(vinput[i],"_cell_angle_beta") && tokens.size()==2){
        a.beta = aurostd::string2utype<double>(tokens[1]); 
      }
      else if(aurostd::substring2bool(vinput[i],"_cell_angle_gamma") && tokens.size()==2){
        a.gamma = aurostd::string2utype<double>(tokens[1]); 
      }
    }
    a.lattice = GetClat(a.a, a.b, a.c, a.alpha, a.beta, a.gamma);
    a.FixLattices();
    a.partial_occupation_flag = FALSE;

    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CIF [6]" << endl;
    // ME20220124 - Convert xyz to symbolic representation
    vector<symbolic::Symbolic> spacegroup_symop_symbolic;
    if (!found_setting) {
      uint nsym = spacegroup_symop_xyz.size();
      vector<vector<string> > vsymops(nsym);
      for (uint i = 0; i < nsym; ++i) {
        aurostd::string2tokens(spacegroup_symop_xyz[i], vsymops[i], ",");
      }
      spacegroup_symop_symbolic = anrl::equations2SymbolicEquations(vsymops);
    }

    // get atoms
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CIF [7]" << endl;
    vector<string> atom_site_fields;
    bool found_atom_site_labels=FALSE;
    bool already_storing_atoms=FALSE;
    for(uint i=0;i<vinput.size();i++) {
      //DX20181210 - CIF can have multiple loops, reset after each loop and neglect aniso loop - START
      if(aurostd::substring2bool(vinput[i],"loop")){
        found_atom_site_labels = FALSE;
      }
      if(aurostd::substring2bool(vinput[i],"_atom_site") && !aurostd::substring2bool(vinput[i],"aniso")){
        if(already_storing_atoms){ break; } //DX20190718 - to handle format of Springer Materials cifs (adds extra fields at the end; do not read them)
        found_atom_site_labels=TRUE;
        atom_site_fields.push_back(vinput[i]);
      }
      else if(found_atom_site_labels && aurostd::RemoveWhiteSpaces(vinput[i]).size()!=0 && 
          vinput[i][0]=='_' && !aurostd::substring2bool(vinput[i],"aniso")){ //DX20190718 - to account for non-standard fields, e.g., _sm_ (Springer Materials) fields; not standard, but abundant
        if(already_storing_atoms){ break; } //DX20190718 - to handle format of Springer Materials cifs (adds extra fields at the end; do not read them)
        atom_site_fields.push_back(vinput[i]);
      }
      //DX20181210 - CIF can have multiple loops, reset after each loop and neglect aniso_U loop - END
      else if(found_atom_site_labels == TRUE && aurostd::RemoveWhiteSpaces(vinput[i]).size()!=0){
        vector<string> tokens;
        aurostd::string2tokens(vinput[i],tokens," ");
        //DX20190718 - check tokens first; Springer Materials has extra spaces - START
        for(uint t=0;t<tokens.size();t++){
          if(aurostd::RemoveWhiteSpaces(tokens[t])[0]=='\'' && aurostd::RemoveWhiteSpaces(tokens[t])[tokens[t].size()-1]!='\''){
            if(t+1<tokens.size()){
              if(aurostd::RemoveWhiteSpaces(tokens[t+1])[0]!='\'' && aurostd::RemoveWhiteSpaces(tokens[t+1])[tokens[t+1].size()-1]=='\''){
                tokens[t]+="_"+tokens[t+1];
                tokens.erase(tokens.begin()+t+1);
              }
            }
          }
        }
        //DX20190718 - check tokens first; Springer Materials has extra spaces - END
        if(tokens.size()==atom_site_fields.size()){
          _atom atom_tmp;
          wyckoffsite_ITC wyckoff_tmp; //DX20191029
          //ME20220124 - prepare for P1 if setting not found
          if (!found_setting) {
            wyckoff_tmp.letter = "a";
            wyckoff_tmp.multiplicity = 1;
            wyckoff_tmp.site_symmetry = "1";
            string eq_str = "x,y,z";
            vector<string> eq;
            aurostd::string2tokens(eq_str, eq, ",");
            wyckoff_tmp.equations.push_back(eq);
          }
          for(uint t=0;t<tokens.size();t++){
            if(aurostd::substring2bool(atom_site_fields.at(t),"_atom_site_type_symbol")){
              string name = aurostd::RemoveCharacterFromTheFrontAndBack(tokens[t],'\''); //DX20190718 - remove surrounding '' (common in Springer Materials cifs)
              // ME20220113 - name could have oxidation states (found in newer ICSD CIFs)
              if (name.length() > 2) name = name.substr(0, isdigit(name[1])?1:2);
              atom_tmp.name = name;
              atom_tmp.name_is_given = TRUE;
            }
            if(aurostd::substring2bool(atom_site_fields.at(t),"_atom_site_fract_x")){ atom_tmp.fpos[1] = aurostd::string2utype<double>(tokens[t]); }
            if(aurostd::substring2bool(atom_site_fields.at(t),"_atom_site_fract_y")){ atom_tmp.fpos[2] = aurostd::string2utype<double>(tokens[t]); }
            if(aurostd::substring2bool(atom_site_fields.at(t),"_atom_site_fract_z")){ atom_tmp.fpos[3] = aurostd::string2utype<double>(tokens[t]); }
            if(aurostd::substring2bool(atom_site_fields.at(t),"_atom_site_occupancy")){ 
              atom_tmp.partial_occupation_value = aurostd::string2utype<double>(tokens[t]); 
              wyckoff_tmp.site_occupation = atom_tmp.partial_occupation_value; //DX20191029 
              if(aurostd::abs(atom_tmp.partial_occupation_value-1.0)<1e-6){ atom_tmp.partial_occupation_flag = FALSE; }
              else { atom_tmp.partial_occupation_flag = TRUE; a.partial_occupation_flag = TRUE; }
            }
            if(aurostd::substring2bool(atom_site_fields.at(t),"_atom_site_symmetry_multiplicity")){ wyckoff_tmp.multiplicity=aurostd::string2utype<double>(tokens[t]); }
            if(aurostd::substring2bool(atom_site_fields.at(t),"_atom_site_Wyckoff_label") || aurostd::substring2bool(atom_site_fields.at(t),"_atom_site_Wyckoff_symbol")){ wyckoff_tmp.letter=tokens[t]; }  //ME20220125 - Added _atoms_site_Wyckoff_symbol check
          }
          wyckoff_tmp.type = atom_tmp.name; //DX20191029
          wyckoff_tmp.coord = atom_tmp.fpos; //DX20191029
          if(found_setting && wyckoff_tmp.multiplicity!=0 && wyckoff_tmp.letter!=""){
            SYM::getWyckoffInformation(a.spacegroupnumber, a.spacegroupoption, wyckoff_tmp.letter, wyckoff_tmp.multiplicity, wyckoff_tmp.site_symmetry, wyckoff_tmp.equations);
          }
          a.wyckoff_sites_ITC.push_back(wyckoff_tmp);
          atom_tmp.cpos=a.f2c*atom_tmp.fpos;
          if(LDEBUG) cerr << __AFLOW_FUNC__ << " pre AddAtom() [i=" << i << "]" << endl;
          a.AddAtom(atom_tmp,true); //CO20230319 - add by species
          if(LDEBUG) cerr << __AFLOW_FUNC__ << " post AddAtom() [i=" << i << "]" << endl;
          //ME20220124 - Use symmetry operations when setting unknown
          if (!found_setting) {
            uint natoms = a.atoms.size();  // For Wyckoff positions
            _atom at;
            deque<_atom> atoms_symop;
            symbolic::Symbolic x("x"), y("y"), z("z"), result;
            for (uint i = 0; i < spacegroup_symop_symbolic.size(); i++) {
              const symbolic::Symbolic& eq = spacegroup_symop_symbolic[i];
              const xvector<double>& fpos = atom_tmp.fpos;
              result = eq[x == fpos[1], y == fpos[2], z == fpos[3]];
              at = atom_tmp;
              at.fpos[1] = (double) result(0);
              at.fpos[2] = (double) result(1);
              at.fpos[3] = (double) result(2);
              at.cpos = a.f2c * at.fpos;
              atoms_symop.push_back(at);
            }
            a.AddAtom(atoms_symop,true);  //CO20230319 - add by species
            // Add Wyckoff positions
            uint natoms_added = a.atoms.size() - natoms;
            for (uint i = 0; i < natoms_added; ++i) a.wyckoff_sites_ITC.push_back(wyckoff_tmp);
          }
          already_storing_atoms=TRUE; //DX20190718 - to handle format of Springer Materials cifs (adds extra fields at the end; do not read them)
        } else {
          message << "Unexpected number of input fields based on _atom_site_[] information (tokens=" << tokens.size() << ", atom_sites_[]=" << atom_site_fields.size() << ")." <<  endl;  //CO20190629
          for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
        }
      }
    }
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CIF [8]" << endl;
    a=WyckoffPOSITIONS(a.spacegroupnumber,a.spacegroupnumberoption,a);
    a.isd=FALSE; // set Selective Dynamics to false
    //DX20191010 - moved loop that used to be here after re-alphabetizing
    a.SpeciesPutAlphabetic(); //DX20190508 - put alphabetic, needed for many AFLOW functions to work properly
    std::stable_sort(a.atoms.begin(),a.atoms.end(),sortAtomsNames); //DX20200312
    std::sort(a.wyckoff_sites_ITC.begin(), a.wyckoff_sites_ITC.end(), sortWyckoffByType); //DX20201014 - sort the Wyckoff positions too
    a.MakeBasis(); //DX20200803 - must be after alphabetic sort
    a.MakeTypes(); //DX20190508 - otherwise types are not created //DX20200803 - must be after alphabetic sort
    //DX20191010 - moved this loop - START
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CIF [9]" << endl;
    for(uint i=0;i<a.atoms.size();i++){
      if(a.atoms[i].partial_occupation_flag==TRUE){
        poccaus.push_back(a.atoms[i].partial_occupation_value);
        a.partial_occupation_sublattice.push_back(a.atoms[i].type);
      }
      else {
        poccaus.push_back(1.0);
        a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
      }
    }
    //DX20191010 - moved this loop - END
    a.is_vasp4_poscar_format=FALSE; //DX20190308 - needed or SPECIES section breaks
    a.is_vasp5_poscar_format=FALSE; //DX20190308 - needed or SPECIES section breaks

    // add title, CIFs do not generally have a canonical "title" line, so make one
    a.BringInCell();  //ME20220124
    a.buildGenericTitle(); //DX20210211
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CIF [DONE]" << endl;
  } // CIF INPUT

  // ----------------------------------------------------------------------
  // AIMS INPUT
  if(a.iomode==IOAIMS_AUTO || a.iomode==IOAIMS_GEOM) { // AIMS GEOM
    //CO20171008
    //remember, we already did all the debugging above
    //if we get here, we can assume geometry.in is solid!
    //if more debugging needed, test above where we detect iomode
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " AIMS IOAIMS_AUTO/IOAIMS_GEOM" << endl;
    // START FROM CELL
    a.scale=1.0; // standard
    a.neg_scale=FALSE; // standard
    a.lattice=aurostd::eye<double>(3,3);//CO20190520

    uint lat_count=0;
    //get lattice first, if available (c2f,f2c)
    for(uint i=0;i<vinput.size();i++){
      aurostd::string2tokens(vinput[i],tokens," ");
      if(tokens[0]=="lattice_vector" && tokens.size()>3){
        lat_count++;
        a.lattice(lat_count,1)=aurostd::string2utype<double>(tokens[1]);
        a.lattice(lat_count,2)=aurostd::string2utype<double>(tokens[2]);
        a.lattice(lat_count,3)=aurostd::string2utype<double>(tokens[3]);
      }
      if(lat_count==3){break;}
    }

    a.FixLattices();
    //a.f2c=trasp(a.lattice), a.c2f=inverse(a.f2c);

    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " AIMS lattice" << endl;
      cerr << a.lattice << endl;
      cerr << __AFLOW_FUNC__ << " AIMS f2c" << endl;
      cerr << a.f2c << endl;
      cerr << __AFLOW_FUNC__ << " AIMS c2f" << endl;
      cerr << a.c2f << endl;
    }

    deque<_atom> atoms;
    //now get atoms
    bool atom_found;
    _atom atom;
    bool atom_props_search=false; //AIMS stores more atom properties UNDER "atom" and before next "atom"
    for(uint i=0;i<vinput.size();i++){
      atom_found=false;
      aurostd::string2tokens(vinput[i],tokens," ");
      if(!tokens.size()){continue;}
      if(aurostd::substring2bool(tokens[0],"atom") && tokens.size()>4){
        atom_props_search=false;
        if(tokens[0]=="atom"){
          atom_found=true;atom.clear();
          atom.cpos[1]=aurostd::string2utype<double>(tokens[1]);
          atom.cpos[2]=aurostd::string2utype<double>(tokens[2]);
          atom.cpos[3]=aurostd::string2utype<double>(tokens[3]);
          atom.fpos=a.c2f*atom.cpos;
          atom.name=atom.cleanname=tokens[4];
          if(LDEBUG) {
            cerr << __AFLOW_FUNC__ << " AIMS atom[" << atom.name <<"] found (cartesian):" << endl;
            cerr << "    cpos" << atom.cpos << endl;
            cerr << "    fpos" << atom.fpos << endl;
          }
        } else if(tokens[0]=="atom_frac"){
          atom_found=true;atom.clear();
          atom.fpos[1]=aurostd::string2utype<double>(tokens[1]);
          atom.fpos[2]=aurostd::string2utype<double>(tokens[2]);
          atom.fpos[3]=aurostd::string2utype<double>(tokens[3]);
          atom.cpos=a.f2c*atom.fpos;
          atom.name=atom.cleanname=tokens[4];
          if(LDEBUG) {
            cerr << __AFLOW_FUNC__ << " AIMS atom[" << atom.name <<"] found (fractional):" << endl;
            cerr << "    fpos" << atom.fpos << endl;
            cerr << "    cpos" << atom.cpos << endl;
          }
        }//ignore else, garbage
        if(atom_found){
          atom_props_search=true; //so we can look for other props on the next line

          //set some defaults! - START
          atom.name_is_given=true;
          //atom.CleanName();
          if(LDEBUG) {
            cerr << __AFLOW_FUNC__ << " AIMS line=" << vinput[i] << endl;
            cerr << __AFLOW_FUNC__ << " AIMS atom.cleanname=" << atom.cleanname << endl;
          }
          atom.CleanSpin();
          //FIXED BELOW atom.number=atoms.size();    // reference position for convasp
          //FIXED BELOW atom.basis=atoms.size();     // position in the basis //MAKE SURE NOT TO SET BASIS HERE, WILL SCREW UP sortAtomsNames() BELOW
          atom.ijk(1)=0;atom.ijk(2)=0;atom.ijk(3)=0; // inside the zero cell...
          atom.corigin(1)=0.0;atom.corigin(2)=0.0;atom.corigin(3)=0.0; // inside the zero cell
          atom.coord(1)=0.0;atom.coord(2)=0.0;atom.coord(3)=0.0; // inside the zero cell
          atom.spin=0.0;
          atom.noncoll_spin.clear(); //DX20171205 - non-collinear spin
          // FIXED BELOW atom.type=itype;                // CONVASP_MODE if I want type 0,1,2,3,...
          atom.order_parameter_value=0;
          atom.order_parameter_atom=FALSE;
          atom.partial_occupation_value=1.0;
          atom.partial_occupation_flag=FALSE;
          //set some defaults! - STOP

          atoms.push_back(atom);
          //wait till later to add to structure, sort first that we don't screw up species
          //a.AddAtom(atom,true); //CO20230319 - add by species
          //a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);  //default
        }
      } else if(atom_props_search && atoms.size()){
        //look and store additional properties here to atoms.back();
      } //ignore else
    }

    //sort first, then assign types
    //MUST BE STABLE SORT, absolutely critical
    //otherwise the order of equivalent names will be mixed EVERY TIME
    //this will screw up forces in APL
    std::stable_sort(atoms.begin(),atoms.end(),sortAtomsNames); //safe because we do AddAtom() below

    //CO20230319 - DX, I believe this code is obsolete now that I patch AddAtom() with add_species, keep for safety
    //assign types
    uint itype=0;
    atoms[0].type=itype;
    for(uint i=1;i<atoms.size();i++){
      if(atoms[i].name!=atoms[i-1].name){itype++;}
      atoms[i].type=itype;
    }

    //a.atoms.clear();  //we clear atoms earlier with RemoveAtom()
    //now add atoms in right order (for species, etc.)
    for(uint i=0;i<atoms.size();i++){
      a.AddAtom(atoms[i],true);  //does num_each_type and comp_each_type  //CO20230319 - add by species
      a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
    }

    a.MakeBasis();

    if(LDEBUG) {
      for(uint i=0;i<a.num_each_type.size();i++){
        cerr << __AFLOW_FUNC__ << " AIMS num_each_type[" << i <<"]=" <<a.num_each_type[i] << ", ";
        cerr << "comp_each_type[" << i <<"]=" <<a.comp_each_type[i] << endl;
      }
    }

    //grab title last
    a.title.clear();
    if(vinput.size()){
      std::size_t pos=vinput[0].find_first_of('#');
      if(pos!=std::string::npos){
        a.title=vinput[0];
        a.title.erase(pos,1);
        a.title=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(a.title);
      }
    }
    if(a.title.empty()){a.buildGenericTitle();}

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " AIMS title=" << a.title << endl;}

    a.partial_occupation_flag=FALSE;
    a.is_vasp4_poscar_format=FALSE;
    a.is_vasp5_poscar_format=FALSE;
    // DONE ?

  } // AIMS INPUT
  // ----------------------------------------------------------------------

  // ----------------------------------------------------------------------
  // ATAT INPUT //SD20220114
  //  Alloy-Theoretic Automated Toolkit
  //  See: https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/manual.pdf
  if(a.iomode==IOATAT_STR) {
    if (LDEBUG) cerr << __AFLOW_FUNC__ << " ATAT IOATAT_STR" << endl;
    a.scale = 1.0;
    a.neg_scale = FALSE;
    xmatrix<double> axes(3, 3), frac_cell(3, 3);

    // read the axes
    uint line = 0;
    uint vec_count = 1;
    for (; line < vinput.size() && vec_count < 4; line++) {
      aurostd::string2tokens(vinput[line], tokens, " ");
      axes(vec_count, 1) = aurostd::string2utype<double>(tokens[0]);
      axes(vec_count, 2) = aurostd::string2utype<double>(tokens[1]);
      axes(vec_count, 3) = aurostd::string2utype<double>(tokens[2]);
      vec_count++;
    }
    // read the fractional cell vectors
    vec_count = 1;
    for (; line < vinput.size() && vec_count < 4; line++) {
      aurostd::string2tokens(vinput[line], tokens, " ");
      frac_cell(vec_count, 1) = aurostd::string2utype<double>(tokens[0]);
      frac_cell(vec_count, 2) = aurostd::string2utype<double>(tokens[1]);
      frac_cell(vec_count, 3) = aurostd::string2utype<double>(tokens[2]);
      vec_count++;
    }
    a.lattice = axes * frac_cell;
    a.FixLattices();
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ATAT lattice" << endl;
      cerr << a.lattice << endl;
      cerr << __AFLOW_FUNC__ << " ATAT f2c" << endl;
      cerr << a.f2c << endl;
      cerr << __AFLOW_FUNC__ << " ATAT c2f" << endl;
      cerr << a.c2f << endl;
    }

    // read atoms
    deque<_atom> atoms;
    _atom atom;
    xvector<double> avec(3);
    for (; line < vinput.size() - 1; line++) {
      atom.clear();
      aurostd::string2tokens(vinput[line], tokens, " ");
      avec(1) = aurostd::string2utype<double>(tokens[0]);
      avec(2) = aurostd::string2utype<double>(tokens[1]);
      avec(3) = aurostd::string2utype<double>(tokens[2]);
      atom.name = atom.cleanname = tokens[3];
      atom.cpos = trasp(axes) * avec;
      atom.fpos = a.c2f * atom.cpos;
      atom.name_is_given = (!atom.name.empty());
      atoms.push_back(atom);
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " ATAT atom[" << atom.name <<"] found:" << endl;
        cerr << "    fpos" << atom.fpos << endl;
        cerr << "    cpos" << atom.cpos << endl;
      }
    }
    std::stable_sort(atoms.begin(),atoms.end(),sortAtomsNames);

    //CO20230319 - DX, I believe this code is obsolete now that I patch AddAtom() with add_species, keep for safety
    // assign types
    uint itype = 0;
    atoms[0].type = itype;
    for (uint i = 1;i < atoms.size(); i++) {
      if (atoms[i].name != atoms[i - 1].name) {itype++;}
      atoms[i].type = itype;
    }

    // add atoms
    for (uint i=0; i < atoms.size(); i++) {
      a.AddAtom(atoms[i],true); //CO20230319 - add by species
      a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
    }
    a.SpeciesPutAlphabetic();

    // add additional attributes
    a.MakeBasis();
    a.MakeTypes();
    a.partial_occupation_flag = FALSE;
    a.is_vasp4_poscar_format = FALSE;
    a.is_vasp5_poscar_format = FALSE;
    a.buildGenericTitle();

  } // ATAT INPUT
  // ----------------------------------------------------------------------

  // COMMON CODE (NEED TO BE PATCHED, THOUGH).
  // FIX NORMAL AND PARTIAL OCCUPAITON
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " COMMON CODE [9]" << endl;
  if(a.partial_occupation_flag==FALSE) {
    // have partial
    a.comp_each_type.clear();
    for(uint i=0;i<a.num_each_type.size();i++)
      a.comp_each_type.push_back((double) a.num_each_type.at(i));
  } else {
    // have partial
    if(poccaus.size()!=a.atoms.size()) {
      message << "poccaus.size()=" << poccaus.size() << " a.atoms.size()=" << a.atoms.size() << " " << endl;  //CO20190629
      for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
    }
    // create list (empty)
    a.comp_each_type.clear();
    for(uint i=0;i<a.num_each_type.size();i++)
      a.comp_each_type.push_back(0.0);
    // screen each one
    for(uint i=0;i<a.atoms.size();i++) {
      if(poccaus.at(i)<1.0-1e-5) {
        a.atoms.at(i).partial_occupation_flag=TRUE;
        a.atoms.at(i).partial_occupation_value=poccaus.at(i);
        a.comp_each_type.at(a.atoms.at(i).type)+=a.atoms.at(i).partial_occupation_value;
      } else {
        a.atoms.at(i).partial_occupation_flag=FALSE;
        a.atoms.at(i).partial_occupation_value=1.0;
        a.comp_each_type.at(a.atoms.at(i).type)+=a.atoms.at(i).partial_occupation_value;
      }
    }
  }
  a.GetStoich();  //CO20170724
  // ---------------------------------------------------------------
  // -------------- CREATE PARTIAL OCCUPATION STUFF
  // cerr << a.atoms.size() << endl;
  // cerr << a.partial_occupation_sublattice.size() << endl;
  // for(uint i=0;i<a.partial_occupation_sublattice.size();i++) cerr << a.partial_occupation_sublattice.at(i) << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " PARTIAL OCCUPATION STUFF [10]" << endl;
  for(uint iatom=0;iatom<a.atoms.size();iatom++) {
    vector<uint> partial_occupation_sublattice_iatom;
    if(a.partial_occupation_sublattice.at(iatom)== _pocc_no_sublattice_) {
      partial_occupation_sublattice_iatom.push_back(iatom);
    } else {
      for(uint jatom=0;jatom<a.atoms.size();jatom++) {
        if(a.atoms.at(jatom).type==a.partial_occupation_sublattice.at(iatom))
          partial_occupation_sublattice_iatom.push_back(iatom);	  
      }
    }
    if(0) {    cerr << a.partial_occupation_sublattice.at(iatom) << ", ";
      for(uint jatom=0;jatom<partial_occupation_sublattice_iatom.size();jatom++)
        cerr << partial_occupation_sublattice_iatom.at(jatom) << " ";
      cerr << endl;
    }
  }
  // ---------------------------------------------------------------
  // REMOVE ORDER PARAMETER
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " REMOVE ORDER PARAMETER [11]" << endl;
  a.order_parameter_atoms.clear();
  for(uint i=0;i<a.atoms.size();i++)
    if(a.atoms.at(i).order_parameter_atom==TRUE)
      a.order_parameter_atoms.push_back(i);
  a.order_parameter_structure=(a.order_parameter_atoms.size()>0);
  // ---------------------------------------------------------------
  // -------------- SPECIES
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " SPECIES [12]" << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " a.is_vasp4_poscar_format=" << a.is_vasp4_poscar_format << endl;
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " a.is_vasp5_poscar_format=" << a.is_vasp5_poscar_format << endl;
  if(a.is_vasp4_poscar_format==TRUE) {
    a.species.clear();a.species_pp.clear();a.species_pp_type.clear();a.species_pp_version.clear();a.species_pp_ZVAL.clear();a.species_pp_vLDAU.clear();a.species_volume.clear();a.species_mass.clear();
    // Plug the species
    for(uint i=0,j=0;i<a.num_each_type.size();j+=a.num_each_type.at(i),i++) {
      a.species.push_back(a.atoms.at(j).name);
      a.species_pp.push_back(a.atoms.at(j).name);
      a.species_pp_type.push_back("");
      a.species_pp_version.push_back("");
      a.species_pp_ZVAL.push_back(0.0);
      a.species_pp_vLDAU.push_back(deque<double>());
      a.species_volume.push_back(0.0);
      a.species_mass.push_back(0.0);
    }
  }
  if(a.is_vasp5_poscar_format==TRUE) {
    a.species.clear();a.species_pp.clear();a.species_pp_type.clear();a.species_pp_version.clear();a.species_pp_ZVAL.clear();a.species_pp_vLDAU.clear();a.species_volume.clear();a.species_mass.clear();
    if(a.iomode==IOVASP_POSCAR) aurostd::string2tokens((vinput[5]),tokens);
    if(a.iomode==IOVASP_ABCCAR) aurostd::string2tokens((vinput[3]),tokens);
    if(a.iomode==IOVASP_WYCKCAR) aurostd::string2tokens((vinput[3]),tokens);
    if(a.num_each_type.size()!=tokens.size()) {
      message << "You need to specify the same number of species and atoms types" << endl;  //CO20190629
      message << "      a.num_each_type.size()=" << a.num_each_type.size() << endl; //CO20190629
      message << "      tokens.size()=" << tokens.size() << endl; //CO20190629
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
    }
    for(uint i=0;i<tokens.size();i++) {
      a.species.push_back(tokens.at(i));
      a.species_pp.push_back(tokens.at(i));
      a.species_pp_type.push_back("");
      a.species_pp_version.push_back("");
      a.species_pp_ZVAL.push_back(0.0);
      a.species_pp_vLDAU.push_back(deque<double>());
      a.species_volume.push_back(0.0);
      a.species_mass.push_back(0.0);
    }
    uint k=0;
    for(uint i=0;i<a.num_each_type.size();i++)
      for(uint j=0;j<(uint) a.num_each_type.at(i);j++) {
        if(a.atoms.at(k).name_is_given==FALSE) {
          a.atoms.at(k).name=a.species.at(i);
          a.atoms.at(k).name_is_given=TRUE;
        }
        k++;
      }
  }
  // ---------------------------------------------------------------
  // ALL atoms have been added. Now add the wyckoff ones
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " WYCKCAR [13]" << endl;
  if(a.iomode==IOVASP_WYCKCAR) {
    // cerr << "[" << a.atoms.size() << "]" << endl;
    a=WyckoffPOSITIONS(a.spacegroupnumber,a.spacegroupnumberoption,a);
    a.title=a.title+" (WYCKOFF "+a.spacegroup+" "+a.spacegrouplabel+")";
    // cerr << "[" << a.atoms.size() << "]" << endl;
  }

  // ---------------------------------------------------------------
  // Make spaces and links inside the qm part
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " QM_SPACE [14]" << endl;
  for(uint i=0;i<a.atoms.size();i++) {
    xvector<double> v(3); v.clear();
    _atom atom;
    atom=a.atoms.at(i);
    atom.cpos.clear();atom.fpos.clear();
    a.qm_forces.push_back(v);
    a.qm_positions.push_back(v);
    a.qm_atoms.push_back(atom);
  }

  // ---------------------------------------------------------------
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " WRAPPING [13]" << endl;
  // TOLERANCES ------------------------
  a.equiv_fpos_epsilon=_EQUIV_FPOS_EPS_; // standard but you can change
  // SORT ATOMS (FALSE) -----------------------------
  // AFLOW prefers alphabetic ordering; HOWEVER, sorting by default can cause
  // issues for functions/processes outside of AFLOW (e.g., settings in INCAR).
  // For now, it is safer to sort the atoms inside the particular AFLOW function
  // where it is needed (e.g., symmetry and prototype functions).
  // For readers/writers other than the VASP geometry file, we will always sort
  // alphabetically. //DX+CO20210706
  bool force_alphabetic_sorting=false;
  if(force_alphabetic_sorting){
    a.SpeciesPutAlphabetic();
    std::stable_sort(a.atoms.begin(),a.atoms.end(),sortAtomsNames);
  }
  // MAKE BASIS
  a.MakeBasis();
  if(force_alphabetic_sorting){ a.MakeTypes(); } //DX+CO20210706
  // FLAGS -----------------------------
  a.Niggli_calculated=FALSE;
  a.Niggli_avoid=FALSE;
  a.Minkowski_calculated=FALSE;
  a.Minkowski_avoid=FALSE;
  a.LatticeReduction_calculated=FALSE;
  a.LatticeReduction_avoid=FALSE;
  // LATTICE FLAGS -----------------------------
  a.Standard_Lattice_calculated=FALSE;
  a.Standard_Lattice_avoid=FALSE;
  a.Standard_Lattice_primitive=FALSE;
  a.Standard_Lattice_conventional=FALSE;
  a.Standard_Lattice_has_failed=FALSE;
  a.bravais_lattice_type="";
  a.bravais_lattice_variation_type="";
  a.bravais_lattice_system="";
  a.bravais_lattice_lattice_type="";
  a.bravais_lattice_lattice_variation_type="";
  a.bravais_lattice_lattice_system="";
  a.pearson_symbol="";
  a.reciprocal_lattice_type="";
  a.reciprocal_lattice_variation_type="";
  a.bravais_superlattice_type="";
  a.bravais_superlattice_variation_type="";
  a.bravais_superlattice_system="";
  a.pearson_symbol_superlattice="";
  // QM CALCULATED STUFF
  a.qm_origin.clear();
  a.qm_scale=0.0;
  a.qm_lattice.clear();
  a.qm_klattice.clear();
  a.qm_f2c.clear();
  a.qm_c2f.clear();
  a.qm_calculated=FALSE;
  a.qm_forces_write=FALSE;
  a.qm_positions_write=FALSE;
  // CALCULATED STUFF
  //DX+CO START
  a.dist_nn_min=AUROSTD_NAN;    //CO
  a.sym_eps=AUROSTD_NAN; //DX
  a.sym_eps_calculated=FALSE; //DX
  a.sym_eps_change_count=0; //DX20180222 - added tolerance count specific to structure
  a.sym_eps_no_scan=false; //DX20210331 - added no scan specific to structure
  //DX+CO END
  a.iatoms_calculated=FALSE;
  a.pgroup_calculated=FALSE;
  a.pgroup_xtal_calculated=FALSE;
  a.pgroupk_Patterson_calculated=FALSE; //DX20200129
  a.pgroupk_calculated=FALSE;
  a.pgroupk_xtal_calculated=FALSE; //DX20171205 - Added pgroupk_xtal
  a.fgroup_calculated=FALSE;
  a.sgroup_calculated=FALSE;
  a.grid_atoms_calculated=FALSE;
  a.lijk_calculated=FALSE;
  //DX20210122 [OBSOLETE] a.neighbors_calculated=FALSE;
  //DX20180712 START
  // ANRL SYMBOLIC MATH
  a.symbolic_math_representation_only=FALSE;
  a.constrained_symmetry_calculation=FALSE;
  //DX20180712 END
  // OUTPUT STUFF
  a.error_flag=FALSE;
  a.error_string="";
  a.write_lattice_flag=FALSE;
  a.write_inequivalent_flag=FALSE;
  a.write_klattice_flag=FALSE;
  a.write_DEBUG_flag=FALSE;
#ifdef XSTR_DEBUG
  a.write_lattice_flag=TRUE;
  a.write_inequivalent_flag=TRUE;
  a.write_klattice_flag=TRUE;
  a.write_DEBUG_flag=TRUE;
#endif
  //RF20200310 BEGIN
  for(uint i=0;i<a.atoms.size();i++) {  //CO20200624
    if (a.atoms.at(i).cleanname.empty()) a.atoms.at(i).CleanName(); //(KBIN::VASP_PseudoPotential_CleanName(a.atoms.at(i).name)); //CO20200624 - fixed CleanName()
  }
  //RF20200310 END
  // CHECKS
  if(a.atoms.size()!=a.qm_atoms.size())     {
    message << "a.atoms.size()!=a.qm_atoms.size() " << endl;  //CO20190629
    for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
  }
  if(a.atoms.size()!=a.qm_forces.size())    {
    message << "a.atoms.size()!=a.qm_forces.size() " << endl; //CO20190629
    for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
  }
  if(a.atoms.size()!=a.qm_positions.size()) {
    message << "a.atoms.size()!=a.qm_positions.size() " << endl;  //CO20190629
    for(uint i=0;i<vinput.size();i++) message << vinput[i] << endl;  //CO20190629
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_); //CO20190629
  }
  if(det(a.lattice)<0.0){ //CO20200201
    message << "Found negative determinant for lattice (det()=" << det(a.lattice) << "). Flip your basis."; //CO20200201
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_); //CO20200201
  } //CO20200201

  // ---------------------------------------------------------------
  if(LDEBUG) cerr << __AFLOW_FUNC__ << " DONE [99]" << endl;
  // DONE
  return cinput;
}

//DX20180124 - added ibrav to lattice - START
// **************************************************************************
// pflow::QE_ibrav2lattice
// **************************************************************************
namespace pflow {
  xmatrix<double> QE_ibrav2lattice(const int& ibrav, const xvector<double>& parameters, const bool& isabc){
    // Lattice based on ibrav and lattice parameters (see http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PW.html#ibrav)
    xmatrix<double> lattice;

    double a,b,c,gamma,beta,alpha=0;

    if(isabc){
      a = parameters(1);               // parameters(1) = a
      b = parameters(2);               // parameters(2) = b
      c = parameters(3);               // parameters(3) = c
    }
    else {
      a = parameters(1)*bohr2angstrom;               // parameters(1) = a   //DX20180215 - added bohr2angstrom (celldm is in Bohr)
      b = parameters(2)*parameters(1)*bohr2angstrom; // parameters(2) = b/a //DX20180215 - added bohr2angstrom (celldm is in Bohr)
      c = parameters(3)*parameters(1)*bohr2angstrom; // parameters(3) = c/a //DX20180215 - added bohr2angstrom (celldm is in Bohr)
    }
    gamma = acos(parameters(4));       // parameters(4) = cos(gamma) NOTE: Different than AFLOW's order convention of a,b,c,alpha,beta,gamma
    beta = acos(parameters(5));        // parameters(5) = cos(beta)  NOTE: Different than AFLOW's order convention of a,b,c,alpha,beta,gamma
    alpha = acos(parameters(6));       // parameters(6) = cos(alpha) NOTE: Different than AFLOW's order convention of a,b,c,alpha,beta,gamma

    xvector<double> xn(3);   xn(1)=1.0;xn(2)=0.0;xn(3)=0.0;
    xvector<double> yn(3);   yn(1)=0.0;yn(2)=1.0;yn(3)=0.0;
    xvector<double> zn(3);   zn(1)=0.0;zn(2)=0.0;zn(3)=1.0;
    xvector<double> a1(3),a2(3),a3(3);

    if(ibrav == 1){ //CUB
      a1 = a*xn; 
      a2 = a*yn; 
      a3 = a*zn; 
    }
    else if(ibrav == 2){ //FCC
      a1=-(1.0/2.0)*a*xn+(1.0/2.0)*a*zn;
      a2=(1.0/2.0)*a*yn+(1.0/2.0)*a*zn;
      a3=-(1.0/2.0)*a*xn+(1.0/2.0)*a*yn;
    }
    else if(ibrav == 3){ //BCC
      a1=(1.0/2.0)*a*xn+(1.0/2.0)*a*yn+(1.0/2.0)*a*zn;
      a2=-(1.0/2.0)*a*xn+(1.0/2.0)*a*yn+(1.0/2.0)*a*zn;
      a3=-(1.0/2.0)*a*xn-(1.0/2.0)*a*yn+(1.0/2.0)*a*zn;
    }
    else if(ibrav == -3){ //BCC (symmetric axis)
      a1=-(1.0/2.0)*a*xn+(1.0/2.0)*a*yn+(1.0/2.0)*a*zn;
      a2=(1.0/2.0)*a*xn-(1.0/2.0)*a*yn+(1.0/2.0)*a*zn;
      a3=(1.0/2.0)*a*xn+(1.0/2.0)*a*yn-(1.0/2.0)*a*zn;
    }
    else if(ibrav == 4){ //Hexagonal/Trigonal
      a1=a*xn;
      a2=-(1.0/2.0)*a*xn+(sqrt(3.0)/2.0)*a*yn;
      a3=c*zn;
    }
    else if(ibrav == 5){ //Trigonal R
      a1=sqrt((1.0-cos(gamma))/2.0)*a*xn-sqrt((1.0-cos(gamma))/6.0)*a*yn+sqrt((1.0+(2.0*cos(gamma)))/3.0)*a*zn;
      a2=2.0*sqrt((1.0-cos(gamma))/6.0)*a*yn+sqrt((1.0+(2.0*cos(gamma)))/3.0)*a*zn;
      a3=-sqrt((1.0-cos(gamma))/2.0)*a*xn-sqrt((1.0-cos(gamma))/6.0)*a*yn+sqrt((1.0+(2.0*cos(gamma)))/3.0)*a*zn;
    }
    else if(ibrav == -5){ //Trigonal R (3-fold axis)
      a1 = (a/sqrt(3.0))*((sqrt((1.0+(2.0*cos(gamma)))/3.0))-2.0*sqrt(2.0)*(sqrt((1.0-cos(gamma))/6.0)))*xn+(a/sqrt(3.0))*((sqrt((1.0+(2.0*cos(gamma)))/3.0))+sqrt(2.0)*(sqrt((1.0-cos(gamma))/6.0)))*yn+(a/sqrt(3.0))*((sqrt((1.0+(2.0*cos(gamma)))/3.0))+sqrt(2.0)*(sqrt((1.0-cos(gamma))/6.0)))*zn; 
      a2 = (a/sqrt(3.0))*((sqrt((1.0+(2.0*cos(gamma)))/3.0))+sqrt(2.0)*(sqrt((1.0-cos(gamma))/6.0)))*xn+(a/sqrt(3.0))*((sqrt((1.0+(2.0*cos(gamma)))/3.0))-2.0*sqrt(2.0)*(sqrt((1.0-cos(gamma))/6.0)))*yn+(a/sqrt(3.0))*((sqrt((1.0+(2.0*cos(gamma)))/3.0))+sqrt(2.0)*(sqrt((1.0-cos(gamma))/6.0)))*zn; 
      a3 = (a/sqrt(3.0))*((sqrt((1.0+(2.0*cos(gamma)))/3.0))+sqrt(2.0)*(sqrt((1.0-cos(gamma))/6.0)))*xn+(a/sqrt(3.0))*((sqrt((1.0+(2.0*cos(gamma)))/3.0))+sqrt(2.0)*(sqrt((1.0-cos(gamma))/6.0)))*yn+(a/sqrt(3.0))*((sqrt((1.0+(2.0*cos(gamma)))/3.0))-2.0*sqrt(2.0)*(sqrt((1.0-cos(gamma))/6.0)))*zn; 
    }
    else if(ibrav == 6){ //TET
      a1=a*xn;
      a2=a*yn;
      a3=c*zn;
    }
    else if(ibrav == 7){ //BCT
      a1=(1.0/2.0)*a*xn-(1.0/2.0)*a*yn+(1.0/2.0)*c*zn;
      a2=(1.0/2.0)*a*xn+(1.0/2.0)*a*yn+(1.0/2.0)*c*zn;
      a3=-(1.0/2.0)*a*xn-(1.0/2.0)*a*yn+(1.0/2.0)*c*zn;
    }
    else if(ibrav == 8){ //ORC
      a1=a*xn;
      a2=b*yn;
      a3=c*zn;
    }
    else if(ibrav == 9){ //ORCC
      a1=(1.0/2.0)*a*xn+(1.0/2.0)*b*yn;
      a2=-(1.0/2.0)*a*xn+(1.0/2.0)*b*yn;
      a3=c*zn;
    }
    else if(ibrav == -9){ //ORCC - 2
      a1=(1.0/2.0)*a*xn-(1.0/2.0)*b*yn;
      a2=(1.0/2.0)*a*xn+(1.0/2.0)*b*yn;
      a3=c*zn;
    }
    else if(ibrav == 10){ //ORCF
      a1=(1.0/2.0)*a*xn+(1.0/2.0)*c*zn;
      a2=(1.0/2.0)*a*xn+(1.0/2.0)*b*yn;
      a3=(1.0/2.0)*b*yn+(1.0/2.0)*c*zn;
    }
    else if(ibrav == 11){ //ORCI
      a1=(1.0/2.0)*a*xn+(1.0/2.0)*b*yn+(1.0/2.0)*c*zn;
      a2=-(1.0/2.0)*a*xn+(1.0/2.0)*b*yn+(1.0/2.0)*c*zn;
      a3=-(1.0/2.0)*a*xn-(1.0/2.0)*b*yn+(1.0/2.0)*c*zn;
    }
    else if(ibrav == 12){ //MCL (unique axis c)
      a1=a*xn;
      a2=cos(gamma)*b*xn+sin(gamma)*b*yn;
      a3=c*zn;
    }
    else if(ibrav == -12){ //MCL (unique axis b)
      a1=a*xn;
      a2=b*yn;
      a3=cos(beta)*c*xn+sin(beta)*c*zn;
    }
    else if(ibrav == 13){ //MCLC
      a1=(1.0/2.0)*a*xn-(1.0/2.0)*c*zn;
      a2=cos(gamma)*b*xn+sin(gamma)*b*yn;
      a3=(1.0/2.0)*a*xn+(1.0/2.0)*c*zn;
    }
    else if(ibrav == 14){ //TRI
      double cx=c*cos(beta);
      double cy=c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma);
      double cz=sqrt(pow(c,2.0)-pow(cx,2.0)-pow(cy,2.0));

      a1=a*xn;
      a2=b*cos(gamma)*xn+b*sin(gamma)*yn;
      a3=cx*xn+cy*yn+cz*zn;
    } 

    lattice(1,1)=a1(1);lattice(1,2)=a1(2);lattice(1,3)=a1(3);
    lattice(2,1)=a2(1);lattice(2,2)=a2(2);lattice(2,3)=a2(3);
    lattice(3,1)=a3(1);lattice(3,2)=a3(2);lattice(3,3)=a3(3);

    return lattice;
  }
}
//DX20180124 - added ibrav to lattice - START

// **************************************************************************
// xstructure::GetStoich
// **************************************************************************
// Get stoichiometries
bool xstructure::GetStoich(void) { //CO20171025
  double total_comp=0.0;
  for(uint i=0;i<comp_each_type.size();i++){
    total_comp+=comp_each_type[i];
  }
  stoich_each_type.clear();
  for(uint i=0;i<comp_each_type.size();i++){
    stoich_each_type.push_back( comp_each_type[i]/total_comp );
  }
  //CO20210916 - round-off printing errors can be a big challenge here
  //with PARTCAR showing: pocc=1, pocc=0.333, pocc=0.333, pocc=0.333 (LIB4/LIB/CNb_svTa_pvTi_sv:PAW_PBE/AB_cF8_225_a_b.AB:POCC_P0-1xA_P1-0.333xB-0.333xC-0.333xD)
  //we get this:
  //stoichiometry=0.50025013,0.16658329,0.16658329,0.16658329
  //instead of this:
  //stoichiometry=0.5,0.166666667,0.166666667,0.166666667
  //the problem is NOT the partial_occupation_values, but the sum and division for stoich
  //try to fix with aurostd::double2fraction()
  if(partial_occupation_flag){
    try{
      bool LDEBUG=(FALSE || XHOST.DEBUG); 
      total_comp=0.0;
      int numerator=0,denominator=0;
      double stoich=0.0;
      deque<double> vstoich;
      for(uint i=0;i<comp_each_type.size();i++){
        aurostd::double2fraction(comp_each_type[i],numerator,denominator,partial_occupation_stoich_tol);  //need the right tolerance here
        stoich=(double)numerator/(double)denominator;
        total_comp+=stoich;
        vstoich.push_back(stoich);
      }
      for(uint i=0;i<comp_each_type.size();i++){vstoich[i]/=total_comp;}
      stoich_each_type=vstoich;
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " stoich_each_type=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(stoich_each_type),",") << endl;}
    }
    catch(aurostd::xerror& re){;} //do nothing
  }
  return TRUE;
}

// **************************************************************************
// xstructure::sortAtomsEquivalent()
// **************************************************************************
// cluster together atoms by equivalent atoms
bool xstructure::sortAtomsEquivalent(void) {
  bool LDEBUG=(FALSE || XHOST.DEBUG); 
  if(partial_occupation_flag==false){
    if(!(*this).iatoms_calculated){pflow::PerformFullSymmetry(*this);}
    if(!(*this).iatoms_calculated){return false;}
  }
  if(LDEBUG) {
    bool write_inequivalent_flag=(*this).write_inequivalent_flag;
    (*this).write_inequivalent_flag=true;
    cerr << __AFLOW_FUNC__ << " structure before iatoms sorting" << endl;
    cerr << (*this) << endl;
    (*this).write_inequivalent_flag=write_inequivalent_flag;
  }
  bool sort_needed=true;
  if(partial_occupation_flag==false){
    sort_needed=false;
    int iatom=0;
    for(uint i=0;i<(*this).iatoms.size()&&!sort_needed;i++){
      for(uint j=0;j<(*this).iatoms[i].size()&&!sort_needed;j++){
        if((*this).iatoms[i][j]!=iatom){sort_needed=true;}
        iatom++;
      }
    }
  }
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " sort " << (sort_needed?string(""):string("NOT ")) << "needed for this structure!" << endl;}
  (*this).MakeBasis();  //make sure the basis is set, sortAtomsEquiv tries not to mess with relative order using basis
  if(!sort_needed){return true;}
  deque<_atom> atoms=(*this).atoms;
  std::stable_sort(atoms.begin(),atoms.end(),sortAtomsEquiv); //safe because we do AddAtom() below
  if(LDEBUG) {
    //check order before AddAtom()
    cerr << __AFLOW_FUNC__ << " newly sorted atoms pre-AddAtom()" << endl;
    bool print_RHT=false;
    bool verbose=false;
    bool print_cartesian=false;
    for(uint i=0;i<atoms.size();i++){
      print_RHT=atoms[i].print_RHT;atoms[i].print_RHT=false;
      verbose=atoms[i].verbose;atoms[i].verbose=false;
      print_cartesian=atoms[i].print_cartesian;atoms[i].print_cartesian=false;
      cerr << atoms[i] << endl;
      atoms[i].print_RHT=print_RHT;
      atoms[i].verbose=verbose;
      atoms[i].print_cartesian=print_cartesian;
    }
  }
  uint atoms_size=(*this).atoms.size();
  for(uint i=atoms_size-1;i<atoms_size;i--){
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " removing atom[" << i << "]" << endl;
    (*this).RemoveAtom(i);
  }
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " all atoms removed" << endl;}
  for(uint i=0;i<atoms.size();i++){
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " adding atom[" << i << "]" << endl;
    (*this).AddAtom(atoms[i],false);  //CO20230319 - add by type
  }
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " all atoms added back" << endl;}
  (*this).ClearSymmetry();  //new structure, clear symmetry
  if(LDEBUG) {
    cerr << __AFLOW_FUNC__ << " newly sorted structure" << endl;
    cerr << (*this) << endl;
  }
  return true;
}

// **************************************************************************
// xstructure::FixLattices
// **************************************************************************
// Fix all the lattices (you can do by hand, but here everything is done)
bool xstructure::FixLattices(void) {
  klattice=ReciprocalLattice(lattice,scale);
  f2c=trasp(lattice);
  c2f=inverse(trasp(lattice));
  // lattice is the reference... everything else depends
  //  if(iomode==IOVASP_POSCAR || iomode==IOVASP_AUTO) {
  xvector<double> data(6);
  data=Getabc_angles(lattice,DEGREES);
  a=data[1];b=data[2];c=data[3];
  alpha=data[4];beta=data[5];gamma=data[6];
  //  }
  //  if(iomode==IOVASP_ABCCAR || iomode==IOVASP_WYCKCAR) {
  //  lattice=GetClat(a,b,c,alpha,beta,gamma);
  // }
  return TRUE;
}


// **************************************************************************
// GetStructure
// **************************************************************************
// get a structure from a directory (POSCAR for VASP
xstructure GetStructure(const int& iomode,ifstream& input) {
  xstructure out;
  if(iomode==IOVASP_POSCAR) { // VASP POSCAR
    out.iomode=IOVASP_POSCAR;
    if(!input) {
      cerr << "EEEEE   File not found" << endl;
      input.clear();input.close();
      out.error_flag=TRUE;
      out.error_string="FILE NOT FOUND";
      return out;
    }
    input >> out;
  }
  return out;
}

// **************************************************************************
// GetStructure
// **************************************************************************
xstructure GetStructure(const int& iomode,const string& Directory) {
  xstructure out;
  if(iomode==IOVASP_POSCAR) { // VASP POSCAR
    out.iomode=IOVASP_POSCAR;
    ifstream input;
    string File;
    File=Directory+"/POSCAR";
    input.open(File.c_str(),std::ios::in);
    out=GetStructure(iomode,input);
    input.clear();input.close();
    return out;
  }
  return out;
}

// **************************************************************************
// xstructure::SetCoordinates
// **************************************************************************
// change coordinates type
void xstructure::SetCoordinates(int mode)  {
  switch(mode) {
    case _UPDATE_LATTICE_VECTORS_TO_ABCANGLES_ : {
                                                   throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"[1] mode="+aurostd::utype2string(mode),_INPUT_ERROR_);
                                                   break;
                                                 }
    case _UPDATE_LATTICE_ABCANGLES_TO_VECTORS_ : {
                                                   throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"[2] mode="+aurostd::utype2string(mode),_INPUT_ERROR_);
                                                   break;
                                                 }
    case _COORDS_CARTESIAN_ : {
                                coord_flag=_COORDS_CARTESIAN_;
                                strcpy(coord_type,"C");
                                break;
                              }
    case _COORDS_FRACTIONAL_ : {
                                 coord_flag=_COORDS_FRACTIONAL_;
                                 strcpy(coord_type,"D");
                                 break;
                               }
    default: {
               cerr << __AFLOW_FUNC__ << " NOTHING TO DO  mode=" << mode << endl;
             }
  }
}


// **************************************************************************
// xstructure::MakeBasis
// **************************************************************************
// This fixes basis and number
void xstructure::MakeBasis(void) {
  // need to update NUMBER and BASIS, number seems restricted to convasp and largely OBSOLETE (CO20190226)
  for(uint iatom=0;iatom<atoms.size();iatom++) {
    atoms.at(iatom).basis=iatom;
    //[CO20200130 - number->basis]atoms.at(iatom).number=iatom;
  }
}

// **************************************************************************
// xstructure::MakeTypes
// **************************************************************************
//CO20180420
void xstructure::MakeTypes(void) {
  // need to update TYPES based on num_each_type
  // type is usually used as an index for species
  // if we take a subset of atoms from another structure (POCC), need to reset first iatom to 0
  stringstream message;
  uint sum_atoms=0;
  for(uint itype=0;itype<num_each_type.size();itype++){sum_atoms+=num_each_type[itype];}
  if(sum_atoms!=atoms.size()){
    message << "num_each_type does not match atom count (sum_atoms=" << sum_atoms << " vs. atoms.size()=" << atoms.size() << ").";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ERROR_);
  }

  uint iat=0;
  for(uint itype=0;itype<num_each_type.size();itype++){
    for(uint j=0;j<(uint)num_each_type[itype];j++){
      atoms[iat++].type=itype;
    }
  }
}

// **************************************************************************
// xstructure::AddAtom() //DX20210202
// **************************************************************************
// This adds a deque<_atom> to the structure.
// More efficient than adding one atom at a time (AddAtom): use more
// efficient for-loop for atoms (upper-triangular) and update species/basis
// info once at the end

void xstructure::AddAtom(const deque<_atom>& _atoms_in, bool add_species, bool check_present) { //DX20210129  //CO20230220 - patching names vs. types bug //CO20230319 - adding add_species which assumes the species information of the incoming atoms is correct, otherwise add by type
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    cerr << __AFLOW_FUNC__ << " add_species=" << add_species << endl;
  }
  
  deque<_atom> atoms_in=_atoms_in;  //make a copy so we can fix types

  //CO20230319 = add_species fixes the "source of truth" problem below, we either add by species (and patch type), or we add by type
  //CO20230220 - found a "source of truth" issue: we cannot trust the incoming types of atoms_in, 
  //as the type many not match with species of (*this)
  //therefore, we need to rely on the name, and if it's not given, (DO NOT) throw an error, we have to pray that type has been set correctly
  //do NOT modify (*this), as we need to check if the atom sits on top of an existing atom before we really add it to (*this)
  //only modify the type of the incoming atoms vector
  //NOTE: MapAtom() below DEPENDS on having correct type implicitly, so this part of the routine MUST come first
  uint iat=0;
  if(add_species){  //CO20230319 - fix atom's type
    deque<string> species_new=species;
    //
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " ORIG species=" << aurostd::joinWDelimiter(species_new,",") << endl;
      cerr << __AFLOW_FUNC__ << " NAMES ORIG atoms_in=";
      for(iat=0;iat<atoms_in.size();iat++){
        cerr << atoms_in[iat].name << (iat<atoms_in.size()-1?",":"");
      }
      cerr << endl;
      cerr << __AFLOW_FUNC__ << " TYPES ORIG atoms_in=";
      for(iat=0;iat<atoms_in.size();iat++){
        cerr << atoms_in[iat].type << (iat<atoms_in.size()-1?",":"");
      }
      cerr << endl;
    }
    //
    bool FOUND_SPECIES=false;
    uint isp=0,species_position=0;
    for(iat=0;iat<atoms_in.size();iat++){
      _atom& atom=atoms_in[iat];
      if(atom.name.empty()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"atom.name.empty()",_INPUT_MISSING_);}
      FOUND_SPECIES=false;
      for(isp=0;isp<species_new.size()&&FOUND_SPECIES==false;isp++){
        if(atom.name==species_new[isp]){FOUND_SPECIES=TRUE;species_position=isp;}
      }
      if(!FOUND_SPECIES){species_new.push_back(atom.name);species_position=species_new.size()-1;}
      atom.type=species_position;
    }
    //
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " TYPES NEW atoms_in=";
      for(iat=0;iat<atoms_in.size();iat++){
        cerr << atoms_in[iat].type << (iat<atoms_in.size()-1?",":"");
      }
      cerr << endl;
      cerr << __AFLOW_FUNC__ << " NEW species=" << aurostd::joinWDelimiter(species_new,",") << endl;
    }
    //
  }
  for(iat=0;iat<atoms_in.size();iat++){
    _atom& atom=atoms_in[iat];
    atom.CleanName(); //moved up from below
    //[DX20170921 - Need to keep spin info]atom.CleanSpin();
    atom.basis=atoms.size()+iat;  //also update basis for sorting below
  }

  uint natoms_xstr = atoms.size();
  const deque<_atom>* ptr_atoms = &atoms_in;
  deque<_atom> atoms_unique;

  if(check_present){
    // use sym_eps if available; if not, use tenth of an Angstrom
    // (since this function adds atoms iteratively, we cannot use minimumDistance,
    // because it would change as we add new atoms) //DX20210202
    double tol=(*this).sym_eps;
    if(tol>=AUROSTD_NAN || tol<_ZERO_TOL_){ tol = 0.1; } // tenth of Angstrom

    // first check if the input atoms are unique
    // it is more efficient to use a double for-loop (upper-triangular)
    // as opposed to MapAtom(deque<_atom>, _atom); otherwise you end up
    // checking atoms twice //DX20210202
    bool FOUND_POSITION=FALSE;
    uint jat=0;
    for(iat=0;iat<atoms_in.size();iat++){
      FOUND_POSITION=FALSE;
      for(jat=iat+1;jat<atoms_in.size()&&FOUND_POSITION==FALSE;jat++){
        if(SYM::MapAtom(atoms_in[iat], atoms_in[jat], true, (*this).lattice, false, tol)){ FOUND_POSITION=TRUE; }
      }
      if(FOUND_POSITION){ continue; }
      // now check if any atoms in the xstructure are duplicates with the input atoms
      else if(natoms_xstr != 0){
        // ME20220120 - changed atoms_unique to atom_in since we are looping over atoms_in
        if(!SYM::MapAtom(atoms, atoms_in[iat], true, (*this).lattice, false, tol)){
          atoms_unique.push_back(atoms_in[iat]);
        }
      }
      // if no atoms in the xstructure, just add to the unique list
      else{ atoms_unique.push_back(atoms_in[iat]); }
    }
    ptr_atoms=&atoms_unique;
  }

  //DX20210202 [OBSOLETE - use MapAtom] if(check_present){
  //DX20210202 [OBSOLETE - use MapAtom]   // check that this atom is not already present
  //DX20210202 [OBSOLETE - use MapAtom]   xvector<double> a1(3),a2(3),a3(3),aijk(3);
  //DX20210202 [OBSOLETE - use MapAtom]   a1=lattice(1);a2=lattice(2);a3=lattice(3);
  //DX20210202 [OBSOLETE - use MapAtom]   for(uint iat=0;iat<atoms_in.size();iat++){
  //DX20210202 [OBSOLETE - use MapAtom]     bool FOUND_POSITION=FALSE;
  //DX20210202 [OBSOLETE - use MapAtom]     for(uint jat=iat+1;jat<atoms_in.size()&&FOUND_POSITION==FALSE;jat++){
  //DX20210202 [OBSOLETE - use MapAtom]       if(atoms_in[iat].type==atoms_in[jat].type && atoms_in[iat].name==atoms_in[jat].name){
  //DX20210202 [OBSOLETE - use MapAtom]         for(int i=-1;i<=1&&FOUND_POSITION==FALSE;i++){
  //DX20210202 [OBSOLETE - use MapAtom]           for(int j=-1;j<=1&&FOUND_POSITION==FALSE;j++){
  //DX20210202 [OBSOLETE - use MapAtom]             for(int k=-1;k<=1&&FOUND_POSITION==FALSE;k++) {
  //DX20210202 [OBSOLETE - use MapAtom]               aijk[1]=i;aijk[2]=j;aijk[3]=k;
  //DX20210202 [OBSOLETE - use MapAtom]               //	if(aurostd::modulus(atoms_in[iat].cpos-(((double)i)*a1+((double)j)*a2+((double)k)*a3+atoms_in[jat].cpos))<0.1) FOUND_POSITION=TRUE;
  //DX20210202 [OBSOLETE - use MapAtom]               //DX+CO START
  //DX20210202 [OBSOLETE - use MapAtom]               //DX if(aurostd::modulus(atoms_in[iat].fpos-(aijk+atoms_in[jat].fpos))<0.01) FOUND_POSITION=TRUE;
  //DX20210202 [OBSOLETE - use MapAtom]               if((*this).sym_eps!=AUROSTD_NAN && (*this).sym_eps<AUROSTD_NAN && (*this).sym_eps>1e-10){ //DX20171201 - Added (*this).sym_eps>1e-10 //DX20180215 - added (*this).sym_eps<AUROSTD_NAN (needed) 
  //DX20210202 [OBSOLETE - use MapAtom]                 if(aurostd::modulus((*this).f2c*(atoms_in[iat].fpos-(aijk+atoms_in[jat].fpos)))<=(*this).sym_eps) FOUND_POSITION=TRUE; //DX
  //DX20210202 [OBSOLETE - use MapAtom]               }
  //DX20210202 [OBSOLETE - use MapAtom]               else {
  //DX20210202 [OBSOLETE - use MapAtom]                 //if(aurostd::modulus(atoms[iat].fpos-(aijk+atom.fpos))<1e-10) FOUND_POSITION=TRUE; //DX
  //DX20210202 [OBSOLETE - use MapAtom]                 if(aurostd::modulus(atoms_in[iat].cpos-(((double)i)*a1+((double)j)*a2+((double)k)*a3+atoms_in[jat].cpos))<0.1) FOUND_POSITION=TRUE; //DX20171201
  //DX20210202 [OBSOLETE - use MapAtom]               }
  //DX20210202 [OBSOLETE - use MapAtom]               //DX+CO END
  //DX20210202 [OBSOLETE - use MapAtom]             }
  //DX20210202 [OBSOLETE - use MapAtom]           }
  //DX20210202 [OBSOLETE - use MapAtom]         }
  //DX20210202 [OBSOLETE - use MapAtom]       }
  //DX20210202 [OBSOLETE - use MapAtom]     }
  //DX20210202 [OBSOLETE - use MapAtom]     if(FOUND_POSITION){ continue; } // found no need to add it further
  //DX20210202 [OBSOLETE - use MapAtom]     else{ atoms.push_back(atoms_in[iat]); }
  //DX20210202 [OBSOLETE - use MapAtom]   }
  //DX20210202 [OBSOLETE - use MapAtom] }
  //DX20210202 [OBSOLETE - use MapAtom] else{
  //DX20210202 [OBSOLETE - use MapAtom]   atoms = atoms_in;
  //DX20210202 [OBSOLETE - use MapAtom] }

  // update the species: update num/comp each type or add new species
  for(iat=0;iat<ptr_atoms->size();iat++){
    UpdateSpecies(ptr_atoms->at(iat),add_species); //DX20210202 - consolidated code below into function
  }

  // add atoms to xstructure
  if(natoms_xstr == 0){ atoms = *ptr_atoms; } //if possible, do assignment instead of push_back (faster)
  else{ 
    for(iat=0;iat<ptr_atoms->size();iat++){
      atoms.push_back(ptr_atoms->at(iat));
    }
  }

  GetStoich();  //CO20170724
  //CO20230220 - BEWARE: this stable sort is NOT meant to sort species
  //it is meant to sort atoms added to the end of the list to the right species-set within already constructed species vector
  //do NOT use sortAtomNames
  //[CO20230220 - does not work, will rearrange differently than UpdateSpecies() above, so species order does not match atom order]//ME20220612 - Originally had sortAtomsTypes, but other parts of xstructure use sortAtomsNames,
  //[CO20230220 - does not work, will rearrange differently than UpdateSpecies() above, so species order does not match atom order]//leading to inconsistencies when the input structure is not alphabetic.
  std::stable_sort(atoms.begin(), atoms.end(), sortAtomsTypes); //[CO20230220 - does not work, will rearrange differently than UpdateSpecies() above, so species order does not match atom order] //sortAtomsNames
  MakeBasis(); // need to update NUMBER and BASIS
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " updated xstr=" << endl << (*this) << endl;}
}

// **************************************************************************
// xstructure::AddAtom
// **************************************************************************
// This adds an atom to the structure.

void xstructure::AddAtom(const _atom& atom, bool add_species, bool check_present) { //CO20230220 - patching names vs. types bug //CO20230319 - adding add_species which assumes the species information of the incoming atoms is correct, otherwise add by type
  bool LDEBUG=(FALSE || XHOST.DEBUG); 
  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    cerr << __AFLOW_FUNC__ << " add_species=" << add_species << endl;
  }
  
  _atom btom=atom; //DX20210202
  
  //CO20230319 = add_species fixes the "source of truth" problem below, we either add by species (and patch type), or we add by type
  //CO20230220 - found a "source of truth" issue: we cannot trust the incoming types of atoms_in, 
  //as the type many not match with species of (*this)
  //therefore, we need to rely on the name, and if it's not given, (DO NOT) throw an error, we have to pray that type has been set correctly
  //do NOT modify (*this), as we need to check if the atom sits on top of an existing atom before we really add it to (*this)
  //only modify the type of the incoming atoms vector
  //NOTE: MapAtom() below DEPENDS on having correct type implicitly, so this part of the routine MUST come first
  if(add_species){  //CO20230319 - fix atom's type
    if(btom.name.empty()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"atom.name.empty()",_INPUT_MISSING_);}
    deque<string> species_new=species;
    //
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " ORIG species=" << aurostd::joinWDelimiter(species_new,",") << endl;
      cerr << __AFLOW_FUNC__ << " TYPES ORIG atom=" << btom.type << endl;
    }
    //
    bool FOUND_SPECIES=false;
    uint species_position=0;
    uint isp=0;
    for(isp=0;isp<species_new.size()&&FOUND_SPECIES==false;isp++){
      if(btom.name==species_new[isp]){FOUND_SPECIES=TRUE;species_position=isp;}
    }
    if(!FOUND_SPECIES){species_new.push_back(btom.name);species_position=species_new.size()-1;}
    btom.type=species_position;
    //
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " TYPES NEW atom=" << btom.type << endl;
      cerr << __AFLOW_FUNC__ << " NEW species=" << aurostd::joinWDelimiter(species_new,",") << endl;
    }
    //
  }
  btom.CleanName(); //moved up from below
  //[DX20170921 - Need to keep spin info]btom.CleanSpin();
  btom.basis=atoms.size();  //also update basis for sorting below

  if(check_present){ //CO20210116 - AddCorners() should NOT check
    // use sym_eps if available; if not, use tenth of an Angstrom
    // (since this function adds atoms iteratively, we cannot use minimumDistance,
    // because it would change as we add new atoms) //DX20210202
    double tol=(*this).sym_eps;
    if(tol>=AUROSTD_NAN || tol<_ZERO_TOL_){ tol = 0.1; } // tenth of Angstrom
    if(SYM::MapAtom((*this).atoms, btom, true, (*this).lattice, false, tol)){ return; }
  }

  //DX20210202 [OBSOLETE - use MapAtom] if(check_present){ //CO20210116 - AddCorners() should NOT check
  //DX20210202 [OBSOLETE - use MapAtom]   // check that this atom is not already present
  //DX20210202 [OBSOLETE - use MapAtom]   xvector<double> a1(3),a2(3),a3(3),aijk(3);
  //DX20210202 [OBSOLETE - use MapAtom]   a1=lattice(1);a2=lattice(2);a3=lattice(3);
  //DX20210202 [OBSOLETE - use MapAtom]   bool FOUND_POSITION=FALSE;
  //DX20210202 [OBSOLETE - use MapAtom]   for(uint iat=0;iat<atoms.size()&&FOUND_POSITION==FALSE;iat++)
  //DX20210202 [OBSOLETE - use MapAtom]     if(atoms[iat].type==atom.type && atoms[iat].name==atom.name)
  //DX20210202 [OBSOLETE - use MapAtom]       for(int i=-1;i<=1&&FOUND_POSITION==FALSE;i++)
  //DX20210202 [OBSOLETE - use MapAtom]         for(int j=-1;j<=1&&FOUND_POSITION==FALSE;j++)
  //DX20210202 [OBSOLETE - use MapAtom]           for(int k=-1;k<=1&&FOUND_POSITION==FALSE;k++) {
  //DX20210202 [OBSOLETE - use MapAtom]             aijk[1]=i;aijk[2]=j;aijk[3]=k;
  //DX20210202 [OBSOLETE - use MapAtom]             //	if(aurostd::modulus(atoms[iat].cpos-(((double)i)*a1+((double)j)*a2+((double)k)*a3+atom.cpos))<0.1) FOUND_POSITION=TRUE;
  //DX20210202 [OBSOLETE - use MapAtom]             //DX+CO START
  //DX20210202 [OBSOLETE - use MapAtom]             //DX if(aurostd::modulus(atoms[iat].fpos-(aijk+atom.fpos))<0.01) FOUND_POSITION=TRUE;
  //DX20210202 [OBSOLETE - use MapAtom]             if((*this).sym_eps!=AUROSTD_NAN && (*this).sym_eps<AUROSTD_NAN && (*this).sym_eps>1e-10){ //DX20171201 - Added (*this).sym_eps>1e-10 //DX20180215 - added (*this).sym_eps<AUROSTD_NAN (needed) 
  //DX20210202 [OBSOLETE - use MapAtom]               if(aurostd::modulus((*this).f2c*(atoms[iat].fpos-(aijk+atom.fpos)))<=(*this).sym_eps) FOUND_POSITION=TRUE; //DX
  //DX20210202 [OBSOLETE - use MapAtom]             }
  //DX20210202 [OBSOLETE - use MapAtom]             else {
  //DX20210202 [OBSOLETE - use MapAtom]               //if(aurostd::modulus(atoms[iat].fpos-(aijk+atom.fpos))<1e-10) FOUND_POSITION=TRUE; //DX
  //DX20210202 [OBSOLETE - use MapAtom]               if(aurostd::modulus(atoms[iat].cpos-(((double)i)*a1+((double)j)*a2+((double)k)*a3+atom.cpos))<0.1) FOUND_POSITION=TRUE; //DX20171201
  //DX20210202 [OBSOLETE - use MapAtom]             }
  //DX20210202 [OBSOLETE - use MapAtom]             //DX+CO END
  //DX20210202 [OBSOLETE - use MapAtom]           }
  //DX20210202 [OBSOLETE - use MapAtom]   if(FOUND_POSITION==TRUE) return; // found no need to add it further
  //DX20210202 [OBSOLETE - use MapAtom] }

  // update the species: update num/comp each type or add new species
  UpdateSpecies(btom,add_species); //DX20210202 - consolidated code below into function
  atoms.push_back(btom);
  //[CO20230220 - moved up]if(btom.name_is_given) {
  //[CO20230220 - moved up]  btom.CleanName();
  //[CO20230220 - moved up]  //DX20170921 - Need to keep spin info  btom.CleanSpin();
  //[CO20230220 - moved up]}
  GetStoich();  //CO20170724
  //CO20230220 - BEWARE: this stable sort is NOT meant to sort species
  //it is meant to sort atoms added to the end of the list to the right species-set within already constructed species vector
  //do NOT use sortAtomNames
  //[CO20230220 - does not work, will rearrange differently than UpdateSpecies() above, so species order does not match atom order]//ME20220612 - Originally had sortAtomsTypes, but other parts of xstructure use sortAtomsNames,
  //[CO20230220 - does not work, will rearrange differently than UpdateSpecies() above, so species order does not match atom order]//leading to inconsistencies when the input structure is not alphabetic.
  std::stable_sort(atoms.begin(), atoms.end(), sortAtomsTypes); //[CO20230220 - does not work, will rearrange differently than UpdateSpecies() above, so species order does not match atom order] //sortAtomsNames

  //[CO20230220 - keep consistent with other AddAtom() function]// OLD STYLE
  //[CO20230220 - keep consistent with other AddAtom() function]//  atoms.push_back(btom);  MakeBasis(); return;
  //[CO20230220 - keep consistent with other AddAtom() function]// NEW STYLE
  //[CO20230220 - keep consistent with other AddAtom() function]bool found=FALSE;
  //[CO20230220 - keep consistent with other AddAtom() function]if(0)  for(uint iat=0;iat<atoms.size()&&!found;iat++) {
  //[CO20230220 - keep consistent with other AddAtom() function]  if(iat<atoms.size()-1) {
  //[CO20230220 - keep consistent with other AddAtom() function]    if(atoms[iat].type==btom.type && atoms.at(iat+1).type!=btom.type) {
  //[CO20230220 - keep consistent with other AddAtom() function]      //	if(LDEBUG)
  //[CO20230220 - keep consistent with other AddAtom() function]      cerr << "HERE1 iat=" << iat << "  atoms[iat].type=" << atoms[iat].type << "  btom.type=" << btom.type << endl;//" atoms.begin()=" <<  long(atoms.begin()) << endl;
  //[CO20230220 - keep consistent with other AddAtom() function]      atoms.insert(iat+atoms.begin()+1,btom); // potential problem  with CAST
  //[CO20230220 - keep consistent with other AddAtom() function]      found=TRUE;
  //[CO20230220 - keep consistent with other AddAtom() function]    }
  //[CO20230220 - keep consistent with other AddAtom() function]  }
  //[CO20230220 - keep consistent with other AddAtom() function]}
  //[CO20230220 - keep consistent with other AddAtom() function]if(1) {
  //[CO20230220 - keep consistent with other AddAtom() function]  // sort by types and partial occupation (highest occupation first)
  //[CO20230220 - keep consistent with other AddAtom() function]  std::deque<_atom>::iterator it=atoms.begin();
  //[CO20230220 - keep consistent with other AddAtom() function]  for(uint iat=0;iat<atoms.size()&&!found;iat++,it++) {
  //[CO20230220 - keep consistent with other AddAtom() function]    if(iat<atoms.size()-1) {
  //[CO20230220 - keep consistent with other AddAtom() function]      //	cerr << "HERE0 iat=" << iat << "  atoms[iat].type=" << atoms[iat].type << "  btom.type=" << btom.type << endl;
  //[CO20230220 - keep consistent with other AddAtom() function]      if((atoms[iat].type==btom.type && atoms.at(iat+1).type!=btom.type) || 
  //[CO20230220 - keep consistent with other AddAtom() function]          (atoms[iat].type==btom.type && atoms.at(iat+1).partial_occupation_value<btom.partial_occupation_value)) {  //CO20180705 - for pocc sorting, larger pocc ahead of smaller pocc
  //[CO20230220 - keep consistent with other AddAtom() function]        //	if(LDEBUG)
  //[CO20230220 - keep consistent with other AddAtom() function]        //	  cerr << "HERE1 iat=" << iat << "  atoms[iat].type=" << atoms[iat].type << "  btom.type=" << btom.type << endl;//" atoms.begin()=" <<  long(atoms.begin()) << endl;
  //[CO20230220 - keep consistent with other AddAtom() function]        atoms.insert(it+1,btom);  // it is iterator, fine for insert.
  //[CO20230220 - keep consistent with other AddAtom() function]        found=TRUE;
  //[CO20230220 - keep consistent with other AddAtom() function]      }
  //[CO20230220 - keep consistent with other AddAtom() function]    }
  //[CO20230220 - keep consistent with other AddAtom() function]  }
  //[CO20230220 - keep consistent with other AddAtom() function]}
  //[CO20230220 - keep consistent with other AddAtom() function]// if never found add at the end
  //[CO20230220 - keep consistent with other AddAtom() function]if(!found) atoms.push_back(btom);
  //[CO20230220 - keep consistent with other AddAtom() function]// 
  //[CO20230220 - keep consistent with other AddAtom() function]//  atoms.push_back(btom);
  // done
  MakeBasis(); // need to update NUMBER and BASIS
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " updated xstr=" << endl << (*this) << endl;}
}

// **************************************************************************
// xstructure::RemoveAtom
// **************************************************************************
// This removes an atom from the structure.
//CO20170721 - added some safety checks to make sure we weren't deleting an entry
// of a vector/deque that didn't exist (see if size() > itype)
// this is really important because if we build an xstructure on the fly and only occupy
// some of these attributes, the code will seg fault badly, and the error will not
// show up until the xstructure is deconstructed (super confusing)
// make sure to code safely always!
void xstructure::RemoveAtom(const uint& iatom) {
  if(iatom<atoms.size()) {
    uint itype=atoms.at(iatom).type;
    if(num_each_type.size()>itype) num_each_type.at(itype)--;
    if(comp_each_type.size()>itype) comp_each_type.at(itype)-=atoms.at(iatom).partial_occupation_value;
    if(num_each_type.size()>itype && num_each_type.at(itype)==0) {
      if(num_each_type.size()>itype)      num_each_type.erase(num_each_type.begin()+itype);    // erase num_each_type
      if(comp_each_type.size()>itype)     comp_each_type.erase(comp_each_type.begin()+itype);  // erase comp_each_type
      if(species.size()>itype)            species.erase(species.begin()+itype);  // erase species
      if(species_pp.size()>itype)         species_pp.erase(species_pp.begin()+itype);  // erase species_pp
      if(species_pp_type.size()>itype)    species_pp_type.erase(species_pp_type.begin()+itype);  // erase species_pp_type
      if(species_pp_version.size()>itype) species_pp_version.erase(species_pp_version.begin()+itype);  // erase species_pp_version
      if(species_pp_ZVAL.size()>itype)    species_pp_ZVAL.erase(species_pp_ZVAL.begin()+itype);  // erase species_pp_ZVAL
      if(species_pp_vLDAU.size()>itype)   species_pp_vLDAU.erase(species_pp_vLDAU.begin()+itype);  // erase species_pp_vLDAU
      if(species_volume.size()>itype)     species_volume.erase(species_volume.begin()+itype);  // erase species_volume
      if(species_mass.size()>itype)       species_mass.erase(species_mass.begin()+itype);  // erase species_mass
      //CO20170721 - might be better if we wrote function like MakeBasis() for types and did 
      //at the end, but it is not unsafe (seg fault) in the way it is written
      for(uint i=0;i<atoms.size();i++)
        if(i!=iatom && atoms.at(i).type>(int)itype)
          atoms.at(i).type--;
    }
    //CO20170721 - this is obsolete with MakeBasis() below!
    // do atoms
    //for(uint i=0;i<atoms.size();i++) {
    //  if(i!=iatom && atoms.at(i).number>atoms.at(iatom).number)
    //    atoms.at(i).number--;
    //  if(i!=iatom && atoms.at(i).basis>atoms.at(iatom).basis)
    //    atoms.at(i).basis--;
    //}
    atoms.erase(atoms.begin()+iatom);
    //    // do partial_occupation_flags
    //    for(uint i=0;i<partial_occupation_flags.size();i++)
    //      if(iatom==partial_occupation_flags.at(i))
    //	partial_occupation_flags.erase(partial_occupation_flags.begin()+i);
    // do order_parameter_atoms
    //CO20170721 - this is okay as it won't seg fault (doesn't delete anything bigger than .size() )
    for(uint i=0;i<order_parameter_atoms.size();i++)
      if(iatom==order_parameter_atoms.at(i))
        order_parameter_atoms.erase(order_parameter_atoms.begin()+i);
  }
  GetStoich();  //CO20170724
  // done
  MakeBasis(); // need to update NUMBER and BASIS
}


void xstructure::RemoveAtom(vector<uint>& v_atoms_to_remove) { //CO20181226
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  std::sort(v_atoms_to_remove.begin(),v_atoms_to_remove.end());v_atoms_to_remove.erase( std::unique( v_atoms_to_remove.begin(), v_atoms_to_remove.end() ), v_atoms_to_remove.end() ); //remove duplicates //CO20181226
  std::sort(v_atoms_to_remove.rbegin(),v_atoms_to_remove.rend()); //NOTE the r, reverse sort, that way when we remove, it doesn't affect other indices
  for(uint atom=0;atom<v_atoms_to_remove.size();atom++){
    if(LDEBUG) {cerr << "Removing Atom " <<  v_atoms_to_remove[atom] << endl;}
    RemoveAtom(v_atoms_to_remove[atom]);
  }
}

void xstructure::RemoveAtom(void) { //DX20210129
  // Removes all atoms from an xstructure and clears the related
  // atom/species variable; faster than removing one at a time

  atoms.clear();
  ClearSpecies(); //clears everything species related
  order_parameter_atoms.clear();
}

void xstructure::ReplaceAtoms(const deque<_atom>& new_atoms, bool check_present, bool sort_species){ //CO20190520 //DX20210129 - added check_present  //CO20230319 - add bool for sort_species
  //this is the SAFEST/CLEANEST way to replace atoms in an xstructure
  //it takes care of num_each_type, species, etc.
  bool LDEBUG=(FALSE || XHOST.DEBUG);

  //DX20210129 [OBSOLETE - remove all at once] for(uint i=atoms.size()-1;i<atoms.size();i--){  //removing atoms
  //DX20210129 [OBSOLETE - remove all at once]   if(LDEBUG) cerr << __AFLOW_FUNC__ << " removing atom[" << i << "]" << endl;
  //DX20210129 [OBSOLETE - remove all at once]   RemoveAtom(i);
  //DX20210129 [OBSOLETE - remove all at once] }
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " removing all atoms" << endl;}
  RemoveAtom(); //DX20210129 - remove all atoms and clear species variables

  if(LDEBUG){cerr << __AFLOW_FUNC__ << " adding new atoms" << endl;}
  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " new_atoms=";
    for(uint i=0;i<new_atoms.size();i++){
      cerr << new_atoms[i].name << (i<new_atoms.size()-1?",":" ");
    }
    cerr << endl;
  }
  //DX20210202 [OBSOLETE] for(uint i=0;i<new_atoms.size();i++){AddAtom(new_atoms[i], true);}  //adding atoms  //CO20230319 - add by species
  AddAtom(new_atoms, true, check_present);  //adding atoms  //CO20230319 - add by species

  if(sort_species){
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " PRE-SPECIES-SORT" << endl;
      cerr << __AFLOW_FUNC__ << " xstr=" << endl << (*this) << endl;
      cerr << __AFLOW_FUNC__ << " xstr.species=" << aurostd::joinWDelimiter(species,",") << endl;
    }

    SpeciesPutAlphabetic(); //DX20210129

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " POST-SPECIES-SORT" << endl;
      cerr << __AFLOW_FUNC__ << " xstr=" << endl << (*this) << endl;
      cerr << __AFLOW_FUNC__ << " xstr.species=" << aurostd::joinWDelimiter(species,",") << endl;
    }
  }

  //CO20230220 - test that species matches names
  uint itype=0,i=0,iatom=0;
  for(itype=0;itype<num_each_type.size();itype++){
    for(i=0;i<(uint)num_each_type[itype];i++){
      if(!atoms[iatom].name.empty() && !species[itype].empty() && atoms[iatom].name!=species[itype]){
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"atoms and species were incorrectly (re)sorted",_RUNTIME_ERROR_);
      }
      iatom++;
    }
  }
  //
}

// **************************************************************************
// xstructure::RemoveCopies xstructure::RemoveFractionalCopies xstructure::RemoveCartesianCopies
// **************************************************************************
// This removes atoms too close than tol
void xstructure::RemoveCopies(double tol) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  bool flag_isprimitive=FALSE,flag_remove=FALSE;
  uint iat1=0,iat2=0,irm=0;
  while(flag_isprimitive==FALSE) {
    if(iat1==atoms.size()) {
      flag_isprimitive=TRUE;
    } else {
      flag_remove=false;
      for(iat2=iat1+1;iat2<atoms.size()&&flag_remove==FALSE;iat2++)
        if((aurostd::modulus(atoms.at(iat1).fpos-atoms.at(iat2).fpos)<tol || aurostd::modulus(atoms.at(iat1).cpos-atoms.at(iat2).cpos)<tol) && flag_remove==FALSE) {
          flag_remove=TRUE;
          irm=iat2;
        }
      if(flag_remove==TRUE) RemoveAtom(irm);
      else iat1++;
    }
    if(LDEBUG) cout << "DEBUG (RemoveCopies) iat1=" << iat1 << endl;
  }
}
void xstructure::RemoveFractionalCopies(double tol) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  bool flag_isprimitive=FALSE,flag_remove=FALSE;
  uint iat1=0,iat2=0,irm=0;
  while(flag_isprimitive==FALSE) {
    if(iat1==atoms.size()) {
      flag_isprimitive=TRUE;
    } else {
      flag_remove=false;
      for(iat2=iat1+1;iat2<atoms.size()&&flag_remove==FALSE;iat2++){
        xvector<double> fdiff = atoms.at(iat1).fpos-atoms.at(iat2).fpos; //DX20180503 - account for PBC
        fdiff = SYM::minimizeDistanceFractionalMethod(fdiff); //DX20190613
        //DX20190613 [OBSOLETE] SYM::PBC(fdiff); //DX20180503 - account for PBC
        if(aurostd::modulus((*this).f2c*fdiff)<tol && flag_remove==FALSE) //DX20180503 - account for PBC and perform in Cartesian space
          //DX20180503 [OBSOLETE] if(aurostd::modulus(atoms.at(iat1).fpos-atoms.at(iat2).fpos)<tol && flag_remove==FALSE)
        { //CO20200106 - patching for auto-indenting
          flag_remove=TRUE;
          irm=iat2;
        }
      }
      if(flag_remove==TRUE) RemoveAtom(irm);
      else iat1++;
    }
    if(LDEBUG) cout << "DEBUG (RemoveFractionalCopies) iat1=" << iat1 << endl;
  }
}
void xstructure::RemoveCartesianCopies(double tol) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  bool flag_isprimitive=FALSE,flag_remove=FALSE;
  uint iat1=0,iat2=0,irm=0;
  while(flag_isprimitive==FALSE) {
    if(iat1==atoms.size()) {
      flag_isprimitive=TRUE;
    } else {
      flag_remove=false;
      for(iat2=iat1+1;iat2<atoms.size()&&flag_remove==FALSE;iat2++)
        if(aurostd::modulus(atoms.at(iat1).cpos-atoms.at(iat2).cpos)<tol && flag_remove==FALSE) {
          flag_remove=TRUE;
          irm=iat2;
        }
      if(flag_remove==TRUE) RemoveAtom(irm);
      else iat1++;
    }
    if(LDEBUG) cout << "DEBUG (RemoveCartesianCopies) iat1=" << iat1 << endl;
  }
}

// **************************************************************************
// xstructure::AddCorners
// **************************************************************************
void xstructure::AddCorners(void) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  xstructure str;
  BringInCell();str=*this;
  while(atoms.size()) RemoveAtom(0);   
  for(uint iat=0;iat<str.atoms.size();iat++) {
    for(double i=0;i<=1;i+=0.99) {
      for(double j=0;j<=1;j+=0.99) {
        for(double k=0;k<=1;k+=0.99) {
          _atom atom=str.atoms.at(iat);
          atom.fpos[1]+=i;atom.fpos[2]+=j;atom.fpos[3]+=k;
          atom.cpos=F2C(lattice,atom.fpos);
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " atom.fpos=" << atom.fpos;}
          if(atom.fpos[1]<=1.0 && atom.fpos[2]<=1.0 && atom.fpos[3]<=1.0){
            if(LDEBUG){cerr << " : adding atom";}
            AddAtom(atom,false,false);  //CO20210116 - do NOT check if atom is already there  //CO20230319 - add by type
          }
          if(LDEBUG){cerr << endl;}
        }
      }
    }
  }
  title=title+" with_corners";
}

// **************************************************************************
// xstructure::ShiftOriginToAtom
// **************************************************************************
// // Shift the origin to atom(iat)
void xstructure::ShiftOriginToAtom(const int& iat) {
  stringstream message;
  //DX+CO START
  if(iat<0 || iat>=(int)atoms.size()) {
    message << "iat=" << iat << " out of boundaries (0," << atoms.size()-1 << ").";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
  }
  xvector<double> frigin(3);
  origin=atoms[iat].cpos;
  frigin=atoms[iat].fpos;
  for(uint i=0;i<atoms.size();i++) {
    atoms[i].fpos=atoms[i].fpos-frigin;
    atoms[i].cpos=atoms[i].cpos-origin;
  }
  //CO TOO SLOW
  //origin=atoms.at(iat).cpos;
  //frigin=atoms.at(iat).fpos;
  //for(uint i=0;i<atoms.size();i++) {
  //  atoms.at(i).fpos=atoms.at(i).fpos-frigin;
  //  atoms.at(i).cpos=atoms.at(i).cpos-origin;
  //}
  //DX+CO END
}

// **************************************************************************
// SetSDNumbers
// **************************************************************************
xstructure SetSDNumbers(const xstructure& a, const vector<string>& in_sd) {
  // Note that in_sd has one name for each type, not one name for each atom.
  xstructure b(a);
  //  int size=in_sd.size();
  for(uint cnt=0;cnt<b.atoms.size();cnt++) {
    b.atoms[cnt].sd=in_sd[cnt];
    if(in_sd[cnt].size()<2) {
      cerr << "WARNING:  Atom=" << cnt << " you must specify SD strings 3 characters long (switching to TTT)" << endl;
      b.atoms.at(cnt).sd = "TTT";
    }
  }
  return b;
}


// **************************************************************************
// SetSDTypes
// **************************************************************************
xstructure SetSDTypes(const xstructure& a, const vector<string>& in_sd) {
  // Note that in_sd has one name for each type, not one name for each atom.
  xstructure b(a);
  int cnt=-1;
  int size=in_sd.size();
  for(int i=0;i<(int)size;i++) {
    if(i<(int)b.num_each_type.size()) {
      for(int j=0;j<(int)b.num_each_type.at(i);j++) {
        cnt++;
        b.atoms.at(cnt).sd=in_sd[i];
        if(in_sd[i].size()<2) {
          cerr << "WARNING:  Atom=" << cnt << " you must specify SD strings 3 characters long (switching to TTT)" << endl;
          b.atoms.at(cnt).sd="TTT";
        }
      }
    }
  }
  return b;
}

// **************************************************************************
// GetTypes
// **************************************************************************
vector<int> GetTypes(const xstructure& a) {
  vector<int> out_type;
  for(int i=0;i<(int)a.atoms.size();i++)
    out_type.push_back(a.atoms.at(i).type);
  return out_type;
}

// **************************************************************************
// GetNames
// **************************************************************************
vector<string> GetNames(const xstructure& a) {
  vector<string> out_name;
  for(int i=0;i<(int)a.atoms.size();i++)
    out_name.push_back(a.atoms.at(i).name);
  return out_name;
}

// **************************************************************************
// GetCleanNames
// **************************************************************************
vector<string> GetCleanNames(const xstructure& a) {
  vector<string> out_cleanname;
  for(int i=0;i<(int)a.atoms.size();i++)
    out_cleanname.push_back(a.atoms.at(i).cleanname);
  return out_cleanname;
}

// **************************************************************************
// GetSpins
// **************************************************************************
vector<double> GetSpins(const xstructure& a) {
  vector<double> out_spin;
  for(int i=0;i<(int)a.atoms.size();i++)
    out_spin.push_back(a.atoms.at(i).spin);
  return out_spin;
}

// ***************************************************************************
// GetElementName
// ***************************************************************************
string GetElementName(string stringin) {
  // need to clean up the _pv stuff of VASP
  for(uint i=0;i<NUM_ELEMENTS;i++)
    if(stringin==vatom_symbol.at(i))
      return vatom_name.at(i);
  return "NotFound";
}

// ***************************************************************************
// GetSpaceGroupName
// ***************************************************************************
string GetSpaceGroupName(int spacegroupnumber, string directory) {
  if(spacegroupnumber < 1 || spacegroupnumber > 230) { //DX20190708 - for xerror
    stringstream message; //DX20190708 - for xerror
    message << "routine: space group specified invalid (1-230): "; //DX20190708 - for xerror
    message << spacegroupnumber << " [dir=" << directory << "]." << endl; //DX20190708 - for xerror
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_); //DX20190708 - for xerror
  }
  // OK
  //DX+ME20190708 - changed subsequent "if" to "else if" -> efficiency
  //ME20210210 - converted to switch statement (faster)
  switch(spacegroupnumber) {
    case 1: // ------------------- 1  P1 #1
      return "P1";
    case 2: // ------------------- 2  P-1 #2
      return "P-1";
    case 3: // ------------------- 3  P2 #3
      return "P2";
    case 4: // ------------------- 4  P2_{1} #4
      return "P2_{1}";
    case 5: // ------------------- 5  C2 #5
      return "C2";
    case 6: // ------------------- 6  Pm #6
      return "Pm";
    case 7: // ------------------- 7  Pc #7
      return "Pc";
    case 8: // ------------------- 8  Cm #8
      return "Cm";
    case 9: // ------------------- 9  Cc #9
      return "Cc";
    case 10: // ------------------- 10  P2/m #10
      return "P2/m";
    case 11: // ------------------- 11  P2_{1}/m #11
      return "P2_{1}/m";
    case 12: // ------------------- 12  C2/m #12
      return "C2/m";
    case 13: // ------------------- 13  P2/c #13
      return "P2/c";
    case 14: // ------------------- 14  P2_{1}/c #14
      return "P2_{1}/c";
    case 15: // ------------------- 15  C2/c #15
      return "C2/c";
    case 16: // ------------------- 16  P222 #16
      return "P222";
    case 17: // ------------------- 17  P222_{1} #17
      return "P222_{1}";
    case 18: // ------------------- 18  P2_{1}2_{1}2 #18
      return "P2_{1}2_{1}2";
    case 19: // ------------------- 19  P2_{1}2_{1}2_{1} #19
      return "P2_{1}2_{1}2_{1}";
    case 20: // ------------------- 20  C222_{1} #20
      return "C222_{1}";
    case 21: // ------------------- 21  C222 #21
      return "C222";
    case 22: // ------------------- 22  F222 #22
      return "F222";
    case 23: // ------------------- 23  I222 #23
      return "I222";
    case 24: // ------------------- 24  I2_{1}2_{1}2_{1} #24
      return "I2_{1}2_{1}2_{1}";
    case 25: // ------------------- 25  Pmm2 #25
      return "Pmm2";
    case 26: // ------------------- 26  Pmc2_{1} #26
      return "Pmc2_{1}";
    case 27: // ------------------- 27  Pcc2 #27
      return "Pcc2";
    case 28: // ------------------- 28  Pma2 #28
      return "Pma2";
    case 29: // ------------------- 29  Pca2_{1} #29
      return "Pca2_{1}";
    case 30: // ------------------- 30  Pnc2 #30
      return "Pnc2";
    case 31: // ------------------- 31  Pmn2_{1} #31
      return "Pmn2_{1}";
    case 32: // ------------------- 32  Pba2 #32
      return "Pba2";
    case 33: // ------------------- 33  Pna2_{1} #33
      return "Pna2_{1}";
    case 34: // ------------------- 34  Pnn2 #34
      return "Pnn2";
    case 35: // ------------------- 35  Cmm2 #35
      return "Cmm2";
    case 36: // ------------------- 36  Cmc2_{1} #36
      return "Cmc2_{1}";
    case 37: // ------------------- 37  Ccc2 #37
      return "Ccc2";
    case 38: // ------------------- 38  Amm2 #38
      return "Amm2";
    case 39: // ------------------- 39  Aem2 #39
      return "Aem2";
    case 40: // ------------------- 40  Ama2 #40
      return "Ama2";
    case 41: // ------------------- 41  Aea2 #41
      return "Aea2";
    case 42: // ------------------- 42  Fmm2 #42
      return "Fmm2";
    case 43: // ------------------- 43  Fdd2 #43
      return "Fdd2";
    case 44: // ------------------- 44  Imm2 #44
      return "Imm2";
    case 45: // ------------------- 45  Iba2 #45
      return "Iba2";
    case 46: // ------------------- 46  Ima2 #46
      return "Ima2";
    case 47: // ------------------- 47  Pmmm #47
      return "Pmmm";
    case 48: // ------------------- 48  Pnnn #48
      return "Pnnn";
    case 49: // ------------------- 49  Pccm #49
      return "Pccm";
    case 50: // ------------------- 50  Pban #50
      return "Pban";
    case 51: // ------------------- 51  Pmma #51
      return "Pmma";
    case 52: // ------------------- 52  Pnna #52
      return "Pnna";
    case 53: // ------------------- 53  Pmna #53
      return "Pmna";
    case 54: // ------------------- 54  Pcca #54
      return "Pcca";
    case 55: // ------------------- 55  Pbam #55
      return "Pbam";
    case 56: // ------------------- 56  Pccn #56
      return "Pccn";
    case 57: // ------------------- 57  Pbcm #57
      return "Pbcm";
    case 58: // ------------------- 58  Pnnm #58
      return "Pnnm";
    case 59: // ------------------- 59  Pmmn #59
      return "Pmmn";
    case 60: // ------------------- 60  Pbcn #60
      return "Pbcn";
    case 61: // ------------------- 61  Pbca #61
      return "Pbca";
    case 62: // ------------------- 62  Pnma #62
      return "Pnma";
    case 63: // ------------------- 63  Cmcm #63
      return "Cmcm";
    case 64: // ------------------- 64  Cmce #64
      return "Cmce";
    case 65: // ------------------- 65  Cmmm #65
      return "Cmmm";
    case 66: // ------------------- 66  Cccm #66
      return "Cccm";
    case 67: // ------------------- 67  Cmme #67
      return "Cmme";
    case 68: // ------------------- 68  Ccce #68
      return "Ccce";
    case 69: // ------------------- 69  Fmmm #69
      return "Fmmm";
    case 70: // ------------------- 70  Fddd #70
      return "Fddd";
    case 71: // ------------------- 71  Immm #71
      return "Immm";
    case 72: // ------------------- 72  Ibam #72
      return "Ibam";
    case 73: // ------------------- 73  Ibca #73
      return "Ibca";
    case 74: // ------------------- 74  Imma #74
      return "Imma";
    case 75: // ------------------- 75  P4 #75
      return "P4";
    case 76: // ------------------- 76  P4_{1} #76
      return "P4_{1}";
    case 77: // ------------------- 77  P4_{2} #77
      return "P4_{2}";
    case 78: // ------------------- 78  P4_{3} #78
      return "P4_{3}";
    case 79: // ------------------- 79  I4 #79
      return "I4";
    case 80: // ------------------- 80  I4_{1} #80
      return "I4_{1}";
    case 81: // ------------------- 81  P-4 #81
      return "P-4";
    case 82: // ------------------- 82  I-4 #82
      return "I-4";
    case 83: // ------------------- 83  P4/m #83
      return "P4/m";
    case 84: // ------------------- 84  P4_{2}/m #84
      return "P4_{2}/m";
    case 85: // ------------------- 85  P4/n #85
      return "P4/n";
    case 86: // ------------------- 86  P4_{2}/n #86
      return "P4_{2}/n";
    case 87: // ------------------- 87  I4/m #87
      return "I4/m";
    case 88: // ------------------- 88  I4_{1}/a #88
      return "I4_{1}/a";
    case 89: // ------------------- 89  P422 #89
      return "P422";
    case 90: // ------------------- 90  P42_{1}2 #90
      return "P42_{1}2";
    case 91: // ------------------- 91  P4_{1}22 #91
      return "P4_{1}22";
    case 92: // ------------------- 92  P4_{1}2_{1}2 #92
      return "P4_{1}2_{1}2";
    case 93: // ------------------- 93  P4_{2}22 #93
      return "P4_{2}22";
    case 94: // ------------------- 94  P4_{2}2_{1}2 #94
      return "P4_{2}2_{1}2";
    case 95: // ------------------- 95  P4_{3}22 #95
      return "P4_{3}22";
    case 96: // ------------------- 96  P4_{3}2_{1}2 #96
      return "P4_{3}2_{1}2";
    case 97: // ------------------- 97  I422 #97
      return "I422";
    case 98: // ------------------- 98  I4_{1}22 #98
      return "I4_{1}22";
    case 99: // ------------------- 99  P4mm #99
      return "P4mm";
    case 100: // ------------------- 100  P4bm #100
      return "P4bm";
    case 101: // ------------------- 101  P4_{2}cm #101
      return "P4_{2}cm";
    case 102: // ------------------- 102  P4_{2}nm #102
      return "P4_{2}nm";
    case 103: // ------------------- 103  P4cc #103
      return "P4cc";
    case 104: // ------------------- 104  P4nc #104
      return "P4nc";
    case 105: // ------------------- 105  P4_{2}mc #105
      return "P4_{2}mc";
    case 106: // ------------------- 106  P4_{2}bc #106
      return "P4_{2}bc";
    case 107: // ------------------- 107  I4mm #107
      return "I4mm";
    case 108: // ------------------- 108  I4cm #108
      return "I4cm";
    case 109: // ------------------- 109  I4_{1}md #109
      return "I4_{1}md";
    case 110: // ------------------- 110  I4_{1}cd #110
      return "I4_{1}cd";
    case 111: // ------------------- 111  P-42m #111
      return "P-42m";
    case 112: // ------------------- 112  P-42c #112
      return "P-42c";
    case 113: // ------------------- 113  P-42_{1}m #113
      return "P-42_{1}m";
    case 114: // ------------------- 114  P-42_{1}c #114
      return "P-42_{1}c";
    case 115: // ------------------- 115  P-4m2 #115
      return "P-4m2";
    case 116: // ------------------- 116  P-4c2 #116
      return "P-4c2";
    case 117: // ------------------- 117  P-4b2 #117
      return "P-4b2";
    case 118: // ------------------- 118  P-4n2 #118
      return "P-4n2";
    case 119: // ------------------- 119  I-4m2 #119
      return "I-4m2";
    case 120: // ------------------- 120  I-4c2 #120
      return "I-4c2";
    case 121: // ------------------- 121  I-42m #121
      return "I-42m";
    case 122: // ------------------- 122  I-42d #122
      return "I-42d";
    case 123: // ------------------- 123  P4/mmm #123
      return "P4/mmm";
    case 124: // ------------------- 124  P4/mcc #124
      return "P4/mcc";
    case 125: // ------------------- 125  P4/nbm #125
      return "P4/nbm";
    case 126: // ------------------- 126  P4/nnc #126
      return "P4/nnc";
    case 127: // ------------------- 127  P4/mbm #127
      return "P4/mbm";
    case 128: // ------------------- 128  P4/mnc #128
      return "P4/mnc";
    case 129: // ------------------- 129  P4/nmm #129
      return "P4/nmm";
    case 130: // ------------------- 130  P4/ncc #130
      return "P4/ncc";
    case 131: // ------------------- 131  P4_{2}/mmc #131
      return "P4_{2}/mmc";
    case 132: // ------------------- 132  P4_{2}/mcm #132
      return "P4_{2}/mcm";
    case 133: // ------------------- 133  P4_{2}/nbc #133
      return "P4_{2}/nbc";
    case 134: // ------------------- 134  P4_{2}/nnm #134
      return "P4_{2}/nnm";
    case 135: // ------------------- 135  P4_{2}/mbc #135
      return "P4_{2}/mbc";
    case 136: // ------------------- 136  P4_{2}/mnm #136
      return "P4_{2}/mnm";
    case 137: // ------------------- 137  P4_{2}/nmc #137
      return "P4_{2}/nmc";
    case 138: // ------------------- 138  P4_{2}/ncm #138
      return "P4_{2}/ncm";
    case 139: // ------------------- 139  I4/mmm #139
      return "I4/mmm";
    case 140: // ------------------- 140  I4/mcm #140
      return "I4/mcm";
    case 141: // ------------------- 141  I4_{1}/amd #141
      return "I4_{1}/amd";
    case 142: // ------------------- 142  I4_{1}/acd #142
      return "I4_{1}/acd";
    case 143: // ------------------- 143  P3 #143
      return "P3";
    case 144: // ------------------- 144  P3_{1} #144
      return "P3_{1}";
    case 145: // ------------------- 145  P3_{2} #145
      return "P3_{2}";
    case 146: // ------------------- 146  R3 #146
      return "R3";
    case 147: // ------------------- 147  P-3 #147
      return "P-3";
    case 148: // ------------------- 148  R-3 #148
      return "R-3";
    case 149: // ------------------- 149  P312 #149
      return "P312";
    case 150: // ------------------- 150  P321 #150
      return "P321";
    case 151: // ------------------- 151  P3_{1}12 #151
      return "P3_{1}12";
    case 152: // ------------------- 152  P3_{1}21 #152
      return "P3_{1}21";
    case 153: // ------------------- 153  P3_{2}12 #153
      return "P3_{2}12";
    case 154: // ------------------- 154  P3_{2}21 #154
      return "P3_{2}21";
    case 155: // ------------------- 155  R32 #155
      return "R32";
    case 156: // ------------------- 156  P3m1 #156
      return "P3m1";
    case 157: // ------------------- 157  P31m #157
      return "P31m";
    case 158: // ------------------- 158  P3c1 #158
      return "P3c1";
    case 159: // ------------------- 159  P31c #159
      return "P31c";
    case 160: // ------------------- 160  R3m #160
      return "R3m";
    case 161: // ------------------- 161  R3c #161
      return "R3c";
    case 162: // ------------------- 162  P-31m #162
      return "P-31m";
    case 163: // ------------------- 163  P-31c #163
      return "P-31c";
    case 164: // ------------------- 164  P-3m1 #164
      return "P-3m1";
    case 165: // ------------------- 165  P-3c1 #165
      return "P-3c1";
    case 166: // ------------------- 166  R-3m #166
      return "R-3m";
    case 167: // ------------------- 167  R-3c #167
      return "R-3c";
    case 168: // ------------------- 168  P6 #168
      return "P6";
    case 169: // ------------------- 169  P6_{1} #169
      return "P6_{1}";
    case 170: // ------------------- 170  P6_{5} #170
      return "P6_{5}";
    case 171: // ------------------- 171  P6_{2} #171
      return "P6_{2}";
    case 172: // ------------------- 172  P6_{4} #172
      return "P6_{4}";
    case 173: // ------------------- 173  P6_{3} #173
      return "P6_{3}";
    case 174: // ------------------- 174  P-6 #174
      return "P-6";
    case 175: // ------------------- 175  P6/m #175
      return "P6/m";
    case 176: // ------------------- 176  P6_{3}/m #176
      return "P6_{3}/m";
    case 177: // ------------------- 177  P622 #177
      return "P622";
    case 178: // ------------------- 178  P6_{1}22 #178
      return "P6_{1}22";
    case 179: // ------------------- 179  P6_{5}22 #179
      return "P6_{5}22";
    case 180: // ------------------- 180  P6_{2}22 #180
      return "P6_{2}22";
    case 181: // ------------------- 181  P6_{4}22 #181
      return "P6_{4}22";
    case 182: // ------------------- 182  P6_{3}22 #182
      return "P6_{3}22";
    case 183: // ------------------- 183  P6mm #183
      return "P6mm";
    case 184: // ------------------- 184  P6cc #184
      return "P6cc";
    case 185: // ------------------- 185  P6_{3}cm #185
      return "P6_{3}cm";
    case 186: // ------------------- 186  P6_{3}mc #186
      return "P6_{3}mc";
    case 187: // ------------------- 187  P-6m2 #187
      return "P-6m2";
    case 188: // ------------------- 188  P-6c2 #188
      return "P-6c2";
    case 189: // ------------------- 189  P-62m #189
      return "P-62m";
    case 190: // ------------------- 190  P-62c #190
      return "P-62c";
    case 191: // ------------------- 191  P6/mmm #191
      return "P6/mmm";
    case 192: // ------------------- 192  P6/mcc #192
      return "P6/mcc";
    case 193: // ------------------- 193  P6_{3}/mcm #193
      return "P6_{3}/mcm";
    case 194: // ------------------- 194  P6_{3}/mmc #194
      return "P6_{3}/mmc";
    case 195: // ------------------- 195  P23 #195
      return "P23";
    case 196: // ------------------- 196  F23 #196
      return "F23";
    case 197: // ------------------- 197  I23 #197
      return "I23";
    case 198: // ------------------- 198  P2_{1}3 #198
      return "P2_{1}3";
    case 199: // ------------------- 199  I2_{1}3 #199
      return "I2_{1}3";
    case 200: // ------------------- 200  Pm-3 #200
      return "Pm-3";
    case 201: // ------------------- 201  Pn-3 #201
      return "Pn-3";
    case 202: // ------------------- 202  Fm-3 #202
      return "Fm-3";
    case 203: // ------------------- 203  Fd-3 #203
      return "Fd-3";
    case 204: // ------------------- 204  Im-3 #204
      return "Im-3";
    case 205: // ------------------- 205  Pa-3 #205
      return "Pa-3";
    case 206: // ------------------- 206  Ia-3 #206
      return "Ia-3";
    case 207: // ------------------- 207  P432 #207
      return "P432";
    case 208: // ------------------- 208  P4_{2}32 #208
      return "P4_{2}32";
    case 209: // ------------------- 209  F432 #209
      return "F432";
    case 210: // ------------------- 210  F4_{1}32 #210
      return "F4_{1}32";
    case 211: // ------------------- 211  I432 #211
      return "I432";
    case 212: // ------------------- 212  P4_{3}32 #212
      return "P4_{3}32";
    case 213: // ------------------- 213  P4_{1}32 #213
      return "P4_{1}32";
    case 214: // ------------------- 214  I4_{1}32 #214
      return "I4_{1}32";
    case 215: // ------------------- 215  P-43m #215
      return "P-43m";
    case 216: // ------------------- 216  F-43m #216
      return "F-43m";
    case 217: // ------------------- 217  I-43m #217
      return "I-43m";
    case 218: // ------------------- 218  P-43n #218
      return "P-43n";
    case 219: // ------------------- 219  F-43c #219
      return "F-43c";
    case 220: // ------------------- 220  I-43d #220
      return "I-43d";
    case 221: // ------------------- 221  Pm-3m #221
      return "Pm-3m";
    case 222: // ------------------- 222  Pn-3n #222
      return "Pn-3n";
    case 223: // ------------------- 223  Pm-3n #223
      return "Pm-3n";
    case 224: // ------------------- 224  Pn-3m #224
      return "Pn-3m";
    case 225: // ------------------- 225  Fm-3m #225
      return "Fm-3m";
    case 226: // ------------------- 226  Fm-3c #226
      return "Fm-3c";
    case 227: // ------------------- 227  Fd-3m #227
      return "Fd-3m";
    case 228: // ------------------- 228  Fd-3c #228
      return "Fd-3c";
    case 229: // ------------------- 229  Im-3m #229
      return "Im-3m";
    case 230: // ------------------- 230  Ia-3d #230
      return "Ia-3d";
  }
  // No spacegroup found
  return "";
}

// ***************************************************************************
// GetSpaceGroupNumber
// ***************************************************************************
int GetSpaceGroupNumber(const string& spacegroupsymbol, string directory) {
  //DX20190708
  stringstream message;
  int spacegroupnumber=0;
  if(spacegroupsymbol[0] != 'P' && spacegroupsymbol[0] != 'I' && spacegroupsymbol[0] != 'F' &&
      spacegroupsymbol[0] != 'R' && spacegroupsymbol[0] != 'C' && spacegroupsymbol[0] != 'A') {
    message << "routine: space group specified invalid (lattice centering not identified: P,I,F,R,C,A): ";
    message << "input symbol=" << spacegroupsymbol << " [dir=" << directory << "]." << endl;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
  }
  // OK
  //DX+ME20190708 - changed subsequent "if" to "else if" -> efficiency
  if(spacegroupsymbol=="P1") {  // ------------------- 1  P1 #1                                                                                 
    spacegroupnumber=1;}
  else if(spacegroupsymbol=="P-1") {  // ------------------- 2  P-1 #2
    spacegroupnumber=2;}
  else if(spacegroupsymbol=="P2") {  // ------------------- 3  P2 #3
    spacegroupnumber=3;}
  else if(spacegroupsymbol=="P2_{1}" || spacegroupsymbol=="P2_1" || spacegroupsymbol=="P21") {  // ------------------- 4  P2_{1}#4
    spacegroupnumber=4;}
  else if(spacegroupsymbol=="C2") {  // ------------------- 5  C2 #5
    spacegroupnumber=5;}
  else if(spacegroupsymbol=="Pm") {  // ------------------- 6  Pm #6
    spacegroupnumber=6;}
  else if(spacegroupsymbol=="Pc") {  // ------------------- 7  Pc #7
    spacegroupnumber=7;}
  else if(spacegroupsymbol=="Cm") {  // ------------------- 8  Cm #8
    spacegroupnumber=8;}
  else if(spacegroupsymbol=="Cc") {  // ------------------- 9  Cc #9
    spacegroupnumber=9;}
  else if(spacegroupsymbol=="P2/m") {  // ------------------- 10  P2/m #10
    spacegroupnumber=10;}
  else if(spacegroupsymbol=="P2_{1}/m" || spacegroupsymbol=="P2_1/m" || spacegroupsymbol=="P21/m") {  // ------------------- 11  P2_{1}/m #11
    spacegroupnumber=11;}
  else if(spacegroupsymbol=="C2/m") {  // ------------------- 12  C2/m #12
    spacegroupnumber=12;}
  else if(spacegroupsymbol=="P2/c") {  // ------------------- 13  P2/c #13
    spacegroupnumber=13;}
  else if(spacegroupsymbol=="P2_{1}/c" || spacegroupsymbol=="P2_1/c" || spacegroupsymbol=="P21/c") {  // ------------------- 14  P2_{1}/c #14
    spacegroupnumber=14;}
  else if(spacegroupsymbol=="C2/c") {  // ------------------- 15  C2/c #15
    spacegroupnumber=15;}
  else if(spacegroupsymbol=="P222") {  // ------------------- 16  P222 #16
    spacegroupnumber=16;}
  else if(spacegroupsymbol=="P222_{1}" || spacegroupsymbol=="P222_1" || spacegroupsymbol=="P2221") {  // ------------------- 17  P222_{1}#17
    spacegroupnumber=17;}
  else if(spacegroupsymbol=="P2_{1}2_{1}2" || spacegroupsymbol=="P2_12_12" || spacegroupsymbol=="P21212") {  // ------------------- 18  P2_{1}2_{1}2 #18
    spacegroupnumber=18;}
  else if(spacegroupsymbol=="P2_{1}2_{1}2_{1}" || spacegroupsymbol=="P2_12_12_1" || spacegroupsymbol=="P212121") {  // ------------------- 19  P2_{1}2_{1}2_{1}#19
    spacegroupnumber=19;}
  else if(spacegroupsymbol=="C222_{1}" || spacegroupsymbol=="C222_1" || spacegroupsymbol=="C2221") {  // ------------------- 20  C222_{1}#20
    spacegroupnumber=20;}
  else if(spacegroupsymbol=="C222") {  // ------------------- 21  C222 #21
    spacegroupnumber=21;}
  else if(spacegroupsymbol=="F222") {  // ------------------- 22  F222 #22
    spacegroupnumber=22;}
  else if(spacegroupsymbol=="I222") {  // ------------------- 23  I222 #23
    spacegroupnumber=23;}
  else if(spacegroupsymbol=="I2_{1}2_{1}2_{1}" || spacegroupsymbol=="I2_12_12_1" || spacegroupsymbol=="I212121") {  // ------------------- 24  I2_{1}2_{1}2_{1}#24
    spacegroupnumber=24;}
  else if(spacegroupsymbol=="Pmm2") {  // ------------------- 25  Pmm2 #25
    spacegroupnumber=25;}
  else if(spacegroupsymbol=="Pmc2_{1}" || spacegroupsymbol=="Pmc2_1" || spacegroupsymbol=="Pmc21") {  // ------------------- 26  Pmc2_{1}#26
    spacegroupnumber=26;}
  else if(spacegroupsymbol=="Pcc2") {  // ------------------- 27  Pcc2 #27
    spacegroupnumber=27;}
  else if(spacegroupsymbol=="Pma2") {  // ------------------- 28  Pma2 #28
    spacegroupnumber=28;}
  else if(spacegroupsymbol=="Pca2_{1}" || spacegroupsymbol=="Pca2_1" || spacegroupsymbol=="Pca21") {  // ------------------- 29  Pca2_{1}#29
    spacegroupnumber=29;}
  else if(spacegroupsymbol=="Pnc2") {  // ------------------- 30  Pnc2 #30
    spacegroupnumber=30;}
  else if(spacegroupsymbol=="Pmn2_{1}" || spacegroupsymbol=="Pmn2_1" || spacegroupsymbol=="Pmn21") {  // ------------------- 31  Pmn2_{1}#31
    spacegroupnumber=31;}
  else if(spacegroupsymbol=="Pba2") {  // ------------------- 32  Pba2 #32
    spacegroupnumber=32;}
  else if(spacegroupsymbol=="Pna2_{1}" || spacegroupsymbol=="Pna2_1" || spacegroupsymbol=="Pna21") {  // ------------------- 33  Pna2_{1}#33
    spacegroupnumber=33;}
  else if(spacegroupsymbol=="Pnn2") {  // ------------------- 34  Pnn2 #34
    spacegroupnumber=34;}
  else if(spacegroupsymbol=="Cmm2") {  // ------------------- 35  Cmm2 #35
    spacegroupnumber=35;}
  else if(spacegroupsymbol=="Cmc2_{1}" || spacegroupsymbol=="Cmc2_1" || spacegroupsymbol=="Cmc21") {  // ------------------- 36  Cmc2_{1}#36
    spacegroupnumber=36;}
  else if(spacegroupsymbol=="Ccc2") {  // ------------------- 37  Ccc2 #37
    spacegroupnumber=37;}
  else if(spacegroupsymbol=="Amm2") {  // ------------------- 38  Amm2 #38
    spacegroupnumber=38;}
  else if(spacegroupsymbol=="Aem2") {  // ------------------- 39  Aem2 #39
    spacegroupnumber=39;}
  else if(spacegroupsymbol=="Ama2") {  // ------------------- 40  Ama2 #40
    spacegroupnumber=40;}
  else if(spacegroupsymbol=="Aea2") {  // ------------------- 41  Aea2 #41
    spacegroupnumber=41;}
  else if(spacegroupsymbol=="Fmm2") {  // ------------------- 42  Fmm2 #42
    spacegroupnumber=42;}
  else if(spacegroupsymbol=="Fdd2") {  // ------------------- 43  Fdd2 #43
    spacegroupnumber=43;}
  else if(spacegroupsymbol=="Imm2") {  // ------------------- 44  Imm2 #44
    spacegroupnumber=44;}
  else if(spacegroupsymbol=="Iba2") {  // ------------------- 45  Iba2 #45
    spacegroupnumber=45;}
  else if(spacegroupsymbol=="Ima2") {  // ------------------- 46  Ima2 #46
    spacegroupnumber=46;}
  else if(spacegroupsymbol=="Pmmm") {  // ------------------- 47  Pmmm #47
    spacegroupnumber=47;}
  else if(spacegroupsymbol=="Pnnn") {  // ------------------- 48  Pnnn #48
    spacegroupnumber=48;}
  else if(spacegroupsymbol=="Pccm") {  // ------------------- 49  Pccm #49
    spacegroupnumber=49;}
  else if(spacegroupsymbol=="Pban") {  // ------------------- 50  Pban #50
    spacegroupnumber=50;}
  else if(spacegroupsymbol=="Pmma") {  // ------------------- 51  Pmma #51
    spacegroupnumber=51;}
  else if(spacegroupsymbol=="Pnna") {  // ------------------- 52  Pnna #52
    spacegroupnumber=52;}
  else if(spacegroupsymbol=="Pmna") {  // ------------------- 53  Pmna #53
    spacegroupnumber=53;}
  else if(spacegroupsymbol=="Pcca") {  // ------------------- 54  Pcca #54
    spacegroupnumber=54;}
  else if(spacegroupsymbol=="Pbam") {  // ------------------- 55  Pbam #55
    spacegroupnumber=55;}
  else if(spacegroupsymbol=="Pccn") {  // ------------------- 56  Pccn #56
    spacegroupnumber=56;}
  else if(spacegroupsymbol=="Pbcm") {  // ------------------- 57  Pbcm #57
    spacegroupnumber=57;}
  else if(spacegroupsymbol=="Pnnm") {  // ------------------- 58  Pnnm #58
    spacegroupnumber=58;}
  else if(spacegroupsymbol=="Pmmn") {  // ------------------- 59  Pmmn #59
    spacegroupnumber=59;}
  else if(spacegroupsymbol=="Pbcn") {  // ------------------- 60  Pbcn #60
    spacegroupnumber=60;}
  else if(spacegroupsymbol=="Pbca") {  // ------------------- 61  Pbca #61
    spacegroupnumber=61;}
  else if(spacegroupsymbol=="Pnma") {  // ------------------- 62  Pnma #62
    spacegroupnumber=62;}
  else if(spacegroupsymbol=="Cmcm") {  // ------------------- 63  Cmcm #63
    spacegroupnumber=63;}
  else if(spacegroupsymbol=="Cmce") {  // ------------------- 64  Cmce #64
    spacegroupnumber=64;}
  else if(spacegroupsymbol=="Cmmm") {  // ------------------- 65  Cmmm #65
    spacegroupnumber=65;}
  else if(spacegroupsymbol=="Cccm") {  // ------------------- 66  Cccm #66
    spacegroupnumber=66;}
  else if(spacegroupsymbol=="Cmme") {  // ------------------- 67  Cmme #67
    spacegroupnumber=67;}
  else if(spacegroupsymbol=="Ccce") {  // ------------------- 68  Ccce #68
    spacegroupnumber=68;}
  else if(spacegroupsymbol=="Fmmm") {  // ------------------- 69  Fmmm #69
    spacegroupnumber=69;}
  else if(spacegroupsymbol=="Fddd") {  // ------------------- 70  Fddd #70
    spacegroupnumber=70;}
  else if(spacegroupsymbol=="Immm") {  // ------------------- 71  Immm #71
    spacegroupnumber=71;}
  else if(spacegroupsymbol=="Ibam") {  // ------------------- 72  Ibam #72
    spacegroupnumber=72;}
  else if(spacegroupsymbol=="Ibca") {  // ------------------- 73  Ibca #73
    spacegroupnumber=73;}
  else if(spacegroupsymbol=="Imma") {  // ------------------- 74  Imma #74
    spacegroupnumber=74;}
  else if(spacegroupsymbol=="P4") {  // ------------------- 75  P4 #75
    spacegroupnumber=75;}
  else if(spacegroupsymbol=="P4_{1}" || spacegroupsymbol=="P4_1" || spacegroupsymbol=="P41") {  // ------------------- 76  P4_{1}#76
    spacegroupnumber=76;}
  else if(spacegroupsymbol=="P4_{2}" || spacegroupsymbol=="P4_2" || spacegroupsymbol=="P42") {  // ------------------- 77  P4_{2}#77
    spacegroupnumber=77;}
  else if(spacegroupsymbol=="P4_{3}" || spacegroupsymbol=="P4_3" || spacegroupsymbol=="P43") {  // ------------------- 78  P4_{3}#78
    spacegroupnumber=78;}
  else if(spacegroupsymbol=="I4") {  // ------------------- 79  I4 #79
    spacegroupnumber=79;}
  else if(spacegroupsymbol=="I4_{1}" || spacegroupsymbol=="I4_1" || spacegroupsymbol=="I41") {  // ------------------- 80  I4_{1}#80
    spacegroupnumber=80;}
  else if(spacegroupsymbol=="P-4") {  // ------------------- 81  P-4 #81
    spacegroupnumber=81;}
  else if(spacegroupsymbol=="I-4") {  // ------------------- 82  I-4 #82
    spacegroupnumber=82;}
  else if(spacegroupsymbol=="P4/m") {  // ------------------- 83  P4/m #83
    spacegroupnumber=83;}
  else if(spacegroupsymbol=="P4_{2}/m" || spacegroupsymbol=="P4_2/m" || spacegroupsymbol=="P42/m") {  // ------------------- 84  P4_{2}/m #84
    spacegroupnumber=84;}
  else if(spacegroupsymbol=="P4/n") {  // ------------------- 85  P4/n #85
    spacegroupnumber=85;}
  else if(spacegroupsymbol=="P4_{2}/n" || spacegroupsymbol=="P4_2/n" || spacegroupsymbol=="P42/n") {  // ------------------- 86  P4_{2}/n #86
    spacegroupnumber=86;}
  else if(spacegroupsymbol=="I4/m") {  // ------------------- 87  I4/m #87
    spacegroupnumber=87;}
  else if(spacegroupsymbol=="I4_{1}/a" || spacegroupsymbol=="I4_1/a" || spacegroupsymbol=="I41/a") {  // ------------------- 88  I4_{1}/a #88
    spacegroupnumber=88;}
  else if(spacegroupsymbol=="P422") {  // ------------------- 89  P422 #89
    spacegroupnumber=89;}
  else if(spacegroupsymbol=="P42_{1}2" || spacegroupsymbol=="P42_12" || spacegroupsymbol=="P4212") {  // ------------------- 90  P42_{1}2 #90
    spacegroupnumber=90;}
  else if(spacegroupsymbol=="P4_{1}22" || spacegroupsymbol=="P4_122" || spacegroupsymbol=="P4122") {  // ------------------- 91  P4_{1}22 #91
    spacegroupnumber=91;}
  else if(spacegroupsymbol=="P4_{1}2_{1}2" || spacegroupsymbol=="P4_12_12" || spacegroupsymbol=="P41212") {  // ------------------- 92  P4_{1}2_{1}2 #92
    spacegroupnumber=92;}
  else if(spacegroupsymbol=="P4_{2}22" || spacegroupsymbol=="P4_222" || spacegroupsymbol=="P4222") {  // ------------------- 93  P4_{2}22 #93
    spacegroupnumber=93;}
  else if(spacegroupsymbol=="P4_{2}2_{1}2" || spacegroupsymbol=="P4_22_12" || spacegroupsymbol=="P42212") {  // ------------------- 94  P4_{2}2_{1}2 #94
    spacegroupnumber=94;}
  else if(spacegroupsymbol=="P4_{3}22" || spacegroupsymbol=="P4_322" || spacegroupsymbol=="P4322") {  // ------------------- 95  P4_{3}22 #95
    spacegroupnumber=95;}
  else if(spacegroupsymbol=="P4_{3}2_{1}2" || spacegroupsymbol=="P4_32_12" || spacegroupsymbol=="P43212") {  // ------------------- 96  P4_{3}2_{1}2 #96
    spacegroupnumber=96;}
  else if(spacegroupsymbol=="I422") {  // ------------------- 97  I422 #97
    spacegroupnumber=97;}
  else if(spacegroupsymbol=="I4_{1}22" || spacegroupsymbol=="I4_122" || spacegroupsymbol=="I4122") {  // ------------------- 98  I4_{1}22 #98
    spacegroupnumber=98;}
  else if(spacegroupsymbol=="P4mm") {  // ------------------- 99  P4mm #99
    spacegroupnumber=99;}
  else if(spacegroupsymbol=="P4bm") {  // ------------------- 100  P4bm #100
    spacegroupnumber=100;}
  else if(spacegroupsymbol=="P4_{2}cm" || spacegroupsymbol=="P4_2cm" || spacegroupsymbol=="P42cm") {  // ------------------- 101  P4_{2}cm #101
    spacegroupnumber=101;}
  else if(spacegroupsymbol=="P4_{2}nm" || spacegroupsymbol=="P4_2nm" || spacegroupsymbol=="P42nm") {  // ------------------- 102  P4_{2}nm #102
    spacegroupnumber=102;}
  else if(spacegroupsymbol=="P4cc") {  // ------------------- 103  P4cc #103
    spacegroupnumber=103;}
  else if(spacegroupsymbol=="P4nc") {  // ------------------- 104  P4nc #104
    spacegroupnumber=104;}
  else if(spacegroupsymbol=="P4_{2}mc" || spacegroupsymbol=="P4_2mc" || spacegroupsymbol=="P42mc") {  // ------------------- 105  P4_{2}mc #105
    spacegroupnumber=105;}
  else if(spacegroupsymbol=="P4_{2}bc" || spacegroupsymbol=="P4_2bc" || spacegroupsymbol=="P42bc") {  // ------------------- 106  P4_{2}bc #106
    spacegroupnumber=106;}
  else if(spacegroupsymbol=="I4mm") {  // ------------------- 107  I4mm #107
    spacegroupnumber=107;}
  else if(spacegroupsymbol=="I4cm") {  // ------------------- 108  I4cm #108
    spacegroupnumber=108;}
  else if(spacegroupsymbol=="I4_{1}md" || spacegroupsymbol=="I4_1md" || spacegroupsymbol=="I41md") {  // ------------------- 109  I4_{1}md #109
    spacegroupnumber=109;}
  else if(spacegroupsymbol=="I4_{1}cd" || spacegroupsymbol=="I4_1cd" || spacegroupsymbol=="I41cd") {  // ------------------- 110  I4_{1}cd #110
    spacegroupnumber=110;}
  else if(spacegroupsymbol=="P-42m") {  // ------------------- 111  P-42m #111
    spacegroupnumber=111;}
  else if(spacegroupsymbol=="P-42c") {  // ------------------- 112  P-42c #112
    spacegroupnumber=112;}
  else if(spacegroupsymbol=="P-42_{1}m" || spacegroupsymbol=="P-42_1m" || spacegroupsymbol=="P-421m") {  // ------------------- 113  P-42_{1}m #113
    spacegroupnumber=113;}
  else if(spacegroupsymbol=="P-42_{1}c" || spacegroupsymbol=="P-42_1c" || spacegroupsymbol=="P-421c") {  // ------------------- 114  P-42_{1}c #114
    spacegroupnumber=114;}
  else if(spacegroupsymbol=="P-4m2") {  // ------------------- 115  P-4m2 #115
    spacegroupnumber=115;}
  else if(spacegroupsymbol=="P-4c2") {  // ------------------- 116  P-4c2 #116
    spacegroupnumber=116;}
  else if(spacegroupsymbol=="P-4b2") {  // ------------------- 117  P-4b2 #117
    spacegroupnumber=117;}
  else if(spacegroupsymbol=="P-4n2") {  // ------------------- 118  P-4n2 #118
    spacegroupnumber=118;}
  else if(spacegroupsymbol=="I-4m2") {  // ------------------- 119  I-4m2 #119
    spacegroupnumber=119;}
  else if(spacegroupsymbol=="I-4c2") {  // ------------------- 120  I-4c2 #120
    spacegroupnumber=120;}
  else if(spacegroupsymbol=="I-42m") {  // ------------------- 121  I-42m #121
    spacegroupnumber=121;}
  else if(spacegroupsymbol=="I-42d") {  // ------------------- 122  I-42d #122
    spacegroupnumber=122;}
  else if(spacegroupsymbol=="P4/mmm") {  // ------------------- 123  P4/mmm #123
    spacegroupnumber=123;}
  else if(spacegroupsymbol=="P4/mcc") {  // ------------------- 124  P4/mcc #124
    spacegroupnumber=124;}
  else if(spacegroupsymbol=="P4/nbm") {  // ------------------- 125  P4/nbm #125
    spacegroupnumber=125;}
  else if(spacegroupsymbol=="P4/nnc") {  // ------------------- 126  P4/nnc #126
    spacegroupnumber=126;}
  else if(spacegroupsymbol=="P4/mbm") {  // ------------------- 127  P4/mbm #127
    spacegroupnumber=127;}
  else if(spacegroupsymbol=="P4/mnc") {  // ------------------- 128  P4/mnc #128
    spacegroupnumber=128;}
  else if(spacegroupsymbol=="P4/nmm") {  // ------------------- 129  P4/nmm #129
    spacegroupnumber=129;}
  else if(spacegroupsymbol=="P4/ncc") {  // ------------------- 130  P4/ncc #130
    spacegroupnumber=130;}
  else if(spacegroupsymbol=="P4_{2}/mmc" || spacegroupsymbol=="P4_2/mmc" || spacegroupsymbol=="P42/mmc") {  // ------------------- 131  P4_{2}/mmc #131
    spacegroupnumber=131;}
  else if(spacegroupsymbol=="P4_{2}/mcm" || spacegroupsymbol=="P4_2/mcm" || spacegroupsymbol=="P42/mcm") {  // ------------------- 132  P4_{2}/mcm #132
    spacegroupnumber=132;}
  else if(spacegroupsymbol=="P4_{2}/nbc" || spacegroupsymbol=="P4_2/nbc" || spacegroupsymbol=="P42/nbc") {  // ------------------- 133  P4_{2}/nbc #133
    spacegroupnumber=133;}
  else if(spacegroupsymbol=="P4_{2}/nnm" || spacegroupsymbol=="P4_2/nnm" || spacegroupsymbol=="P42/nnm") {  // ------------------- 134  P4_{2}/nnm #134
    spacegroupnumber=134;}
  else if(spacegroupsymbol=="P4_{2}/mbc" || spacegroupsymbol=="P4_2/mbc" || spacegroupsymbol=="P42/mbc") {  // ------------------- 135  P4_{2}/mbc #135
    spacegroupnumber=135;}
  else if(spacegroupsymbol=="P4_{2}/mnm" || spacegroupsymbol=="P4_2/mnm" || spacegroupsymbol=="P42/mnm") {  // ------------------- 136  P4_{2}/mnm #136
    spacegroupnumber=136;}
  else if(spacegroupsymbol=="P4_{2}/nmc" || spacegroupsymbol=="P4_2/nmc" || spacegroupsymbol=="P42/nmc") {  // ------------------- 137  P4_{2}/nmc #137
    spacegroupnumber=137;}
  else if(spacegroupsymbol=="P4_{2}/ncm" || spacegroupsymbol=="P4_2/ncm" || spacegroupsymbol=="P42/ncm") {  // ------------------- 138  P4_{2}/ncm #138
    spacegroupnumber=138;}
  else if(spacegroupsymbol=="I4/mmm") {  // ------------------- 139  I4/mmm #139
    spacegroupnumber=139;}
  else if(spacegroupsymbol=="I4/mcm") {  // ------------------- 140  I4/mcm #140
    spacegroupnumber=140;}
  else if(spacegroupsymbol=="I4_{1}/amd" || spacegroupsymbol=="I4_1/amd" || spacegroupsymbol=="I41/amd") {  // ------------------- 141  I4_{1}/amd #141
    spacegroupnumber=141;}
  else if(spacegroupsymbol=="I4_{1}/acd" || spacegroupsymbol=="I4_1/acd" || spacegroupsymbol=="I41/acd") {  // ------------------- 142  I4_{1}/acd #142
    spacegroupnumber=142;}
  else if(spacegroupsymbol=="P3") {  // ------------------- 143  P3 #143
    spacegroupnumber=143;}
  else if(spacegroupsymbol=="P3_{1}" || spacegroupsymbol=="P3_1" || spacegroupsymbol=="P31") {  // ------------------- 144  P3_{1}#144
    spacegroupnumber=144;}
  else if(spacegroupsymbol=="P3_{2}" || spacegroupsymbol=="P3_2" || spacegroupsymbol=="P32") {  // ------------------- 145  P3_{2}#145
    spacegroupnumber=145;}
  else if(spacegroupsymbol=="R3") {  // ------------------- 146  R3 #146
    spacegroupnumber=146;}
  else if(spacegroupsymbol=="P-3") {  // ------------------- 147  P-3 #147
    spacegroupnumber=147;}
  else if(spacegroupsymbol=="R-3") {  // ------------------- 148  R-3 #148
    spacegroupnumber=148;}
  else if(spacegroupsymbol=="P312") {  // ------------------- 149  P312 #149
    spacegroupnumber=149;}
  else if(spacegroupsymbol=="P321") {  // ------------------- 150  P321 #150
    spacegroupnumber=150;}
  else if(spacegroupsymbol=="P3_{1}12" || spacegroupsymbol=="P3_112" || spacegroupsymbol=="P3112") {  // ------------------- 151  P3_{1}12 #151
    spacegroupnumber=151;}
  else if(spacegroupsymbol=="P3_{1}21" || spacegroupsymbol=="P3_121" || spacegroupsymbol=="P3121") {  // ------------------- 152  P3_{1}21 #152
    spacegroupnumber=152;}
  else if(spacegroupsymbol=="P3_{2}12" || spacegroupsymbol=="P3_212" || spacegroupsymbol=="P3212") {  // ------------------- 153  P3_{2}12 #153
    spacegroupnumber=153;}
  else if(spacegroupsymbol=="P3_{2}21" || spacegroupsymbol=="P3_221" || spacegroupsymbol=="P3221") {  // ------------------- 154  P3_{2}21 #154
    spacegroupnumber=154;}
  else if(spacegroupsymbol=="R32") {  // ------------------- 155  R32 #155
    spacegroupnumber=155;}
  else if(spacegroupsymbol=="P3m1") {  // ------------------- 156  P3m1 #156
    spacegroupnumber=156;}
  else if(spacegroupsymbol=="P31m") {  // ------------------- 157  P31m #157
    spacegroupnumber=157;}
  else if(spacegroupsymbol=="P3c1") {  // ------------------- 158  P3c1 #158
    spacegroupnumber=158;}
  else if(spacegroupsymbol=="P31c") {  // ------------------- 159  P31c #159
    spacegroupnumber=159;}
  else if(spacegroupsymbol=="R3m") {  // ------------------- 160  R3m #160
    spacegroupnumber=160;}
  else if(spacegroupsymbol=="R3c") {  // ------------------- 161  R3c #161
    spacegroupnumber=161;}
  else if(spacegroupsymbol=="P-31m") {  // ------------------- 162  P-31m #162
    spacegroupnumber=162;}
  else if(spacegroupsymbol=="P-31c") {  // ------------------- 163  P-31c #163
    spacegroupnumber=163;}
  else if(spacegroupsymbol=="P-3m1") {  // ------------------- 164  P-3m1 #164
    spacegroupnumber=164;}
  else if(spacegroupsymbol=="P-3c1") {  // ------------------- 165  P-3c1 #165
    spacegroupnumber=165;}
  else if(spacegroupsymbol=="R-3m") {  // ------------------- 166  R-3m #166
    spacegroupnumber=166;}
  else if(spacegroupsymbol=="R-3c") {  // ------------------- 167  R-3c #167
    spacegroupnumber=167;}
  else if(spacegroupsymbol=="P6") {  // ------------------- 168  P6 #168
    spacegroupnumber=168;}
  else if(spacegroupsymbol=="P6_{1}" || spacegroupsymbol=="P6_1" || spacegroupsymbol=="P61") {  // ------------------- 169  P6_{1}#169
    spacegroupnumber=169;}
  else if(spacegroupsymbol=="P6_{5}" || spacegroupsymbol=="P6_5" || spacegroupsymbol=="P65") {  // ------------------- 170  P6_{5}#170
    spacegroupnumber=170;}
  else if(spacegroupsymbol=="P6_{2}" || spacegroupsymbol=="P6_2" || spacegroupsymbol=="P62") {  // ------------------- 171  P6_{2}#171
    spacegroupnumber=171;}
  else if(spacegroupsymbol=="P6_{4}" || spacegroupsymbol=="P6_4" || spacegroupsymbol=="P64") {  // ------------------- 172  P6_{4}#172
    spacegroupnumber=172;}
  else if(spacegroupsymbol=="P6_{3}" || spacegroupsymbol=="P6_3" || spacegroupsymbol=="P63") {  // ------------------- 173  P6_{3}#173
    spacegroupnumber=173;}
  else if(spacegroupsymbol=="P-6") {  // ------------------- 174  P-6 #174
    spacegroupnumber=174;}
  else if(spacegroupsymbol=="P6/m") {  // ------------------- 175  P6/m #175
    spacegroupnumber=175;}
  else if(spacegroupsymbol=="P6_{3}/m" || spacegroupsymbol=="P6_3/m" || spacegroupsymbol=="P63/m") {  // ------------------- 176  P6_{3}/m #176
    spacegroupnumber=176;}
  else if(spacegroupsymbol=="P622") {  // ------------------- 177  P622 #177
    spacegroupnumber=177;}
  else if(spacegroupsymbol=="P6_{1}22" || spacegroupsymbol=="P6_122" || spacegroupsymbol=="P6122") {  // ------------------- 178  P6_{1}22 #178
    spacegroupnumber=178;}
  else if(spacegroupsymbol=="P6_{5}22" || spacegroupsymbol=="P6_522" || spacegroupsymbol=="P6522") {  // ------------------- 179  P6_{5}22 #179
    spacegroupnumber=179;}
  else if(spacegroupsymbol=="P6_{2}22" || spacegroupsymbol=="P6_222" || spacegroupsymbol=="P6222") {  // ------------------- 180  P6_{2}22 #180
    spacegroupnumber=180;}
  else if(spacegroupsymbol=="P6_{4}22" || spacegroupsymbol=="P6_422" || spacegroupsymbol=="P6422") {  // ------------------- 181  P6_{4}22 #181
    spacegroupnumber=181;}
  else if(spacegroupsymbol=="P6_{3}22" || spacegroupsymbol=="P6_322" || spacegroupsymbol=="P6322") {  // ------------------- 182  P6_{3}22 #182
    spacegroupnumber=182;}
  else if(spacegroupsymbol=="P6mm") {  // ------------------- 183  P6mm #183
    spacegroupnumber=183;}
  else if(spacegroupsymbol=="P6cc") {  // ------------------- 184  P6cc #184
    spacegroupnumber=184;}
  else if(spacegroupsymbol=="P6_{3}cm" || spacegroupsymbol=="P6_3cm" || spacegroupsymbol=="P63cm") {  // ------------------- 185  P6_{3}cm #185
    spacegroupnumber=185;}
  else if(spacegroupsymbol=="P6_{3}mc" || spacegroupsymbol=="P6_3mc" || spacegroupsymbol=="P63mc") {  // ------------------- 186  P6_{3}mc #186
    spacegroupnumber=186;}
  else if(spacegroupsymbol=="P-6m2") {  // ------------------- 187  P-6m2 #187
    spacegroupnumber=187;}
  else if(spacegroupsymbol=="P-6c2") {  // ------------------- 188  P-6c2 #188
    spacegroupnumber=188;}
  else if(spacegroupsymbol=="P-62m") {  // ------------------- 189  P-62m #189
    spacegroupnumber=189;}
  else if(spacegroupsymbol=="P-62c") {  // ------------------- 190  P-62c #190
    spacegroupnumber=190;}
  else if(spacegroupsymbol=="P6/mmm") {  // ------------------- 191  P6/mmm #191
    spacegroupnumber=191;}
  else if(spacegroupsymbol=="P6/mcc") {  // ------------------- 192  P6/mcc #192
    spacegroupnumber=192;}
  else if(spacegroupsymbol=="P6_{3}/mcm" || spacegroupsymbol=="P6_3/mcm" || spacegroupsymbol=="P63/mcm") {  // ------------------- 193  P6_{3}/mcm #193
    spacegroupnumber=193;}
  else if(spacegroupsymbol=="P6_{3}/mmc" || spacegroupsymbol=="P6_3/mmc" || spacegroupsymbol=="P63/mmc") {  // ------------------- 194  P6_{3}/mmc #194
    spacegroupnumber=194;}
  else if(spacegroupsymbol=="P23") {  // ------------------- 195  P23 #195
    spacegroupnumber=195;}
  else if(spacegroupsymbol=="F23") {  // ------------------- 196  F23 #196
    spacegroupnumber=196;}
  else if(spacegroupsymbol=="I23") {  // ------------------- 197  I23 #197
    spacegroupnumber=197;}
  else if(spacegroupsymbol=="P2_{1}3" || spacegroupsymbol=="P2_13" || spacegroupsymbol=="P213") {  // ------------------- 198  P2_{1}3 #198
    spacegroupnumber=198;}
  else if(spacegroupsymbol=="I2_{1}3" || spacegroupsymbol=="I2_13" || spacegroupsymbol=="I213") {  // ------------------- 199  I2_{1}3 #199
    spacegroupnumber=199;}
  else if(spacegroupsymbol=="Pm-3") {  // ------------------- 200  Pm-3 #200
    spacegroupnumber=200;}
  else if(spacegroupsymbol=="Pn-3") {  // ------------------- 201  Pn-3 #201
    spacegroupnumber=201;}
  else if(spacegroupsymbol=="Fm-3") {  // ------------------- 202  Fm-3 #202
    spacegroupnumber=202;}
  else if(spacegroupsymbol=="Fd-3") {  // ------------------- 203  Fd-3 #203
    spacegroupnumber=203;}
  else if(spacegroupsymbol=="Im-3") {  // ------------------- 204  Im-3 #204
    spacegroupnumber=204;}
  else if(spacegroupsymbol=="Pa-3") {  // ------------------- 205  Pa-3 #205
    spacegroupnumber=205;}
  else if(spacegroupsymbol=="Ia-3") {  // ------------------- 206  Ia-3 #206
    spacegroupnumber=206;}
  else if(spacegroupsymbol=="P432") {  // ------------------- 207  P432 #207
    spacegroupnumber=207;}
  else if(spacegroupsymbol=="P4_{2}32" || spacegroupsymbol=="P4_232" || spacegroupsymbol=="P4232") {  // ------------------- 208  P4_{2}32 #208
    spacegroupnumber=208;}
  else if(spacegroupsymbol=="F432") {  // ------------------- 209  F432 #209
    spacegroupnumber=209;}
  else if(spacegroupsymbol=="F4_{1}32" || spacegroupsymbol=="F4_132" || spacegroupsymbol=="F4132") {  // ------------------- 210  F4_{1}32 #210
    spacegroupnumber=210;}
  else if(spacegroupsymbol=="I432") {  // ------------------- 211  I432 #211
    spacegroupnumber=211;}
  else if(spacegroupsymbol=="P4_{3}32" || spacegroupsymbol=="P4_332" || spacegroupsymbol=="P4332") {  // ------------------- 212  P4_{3}32 #212
    spacegroupnumber=212;}
  else if(spacegroupsymbol=="P4_{1}32" || spacegroupsymbol=="P4_132" || spacegroupsymbol=="P4132") {  // ------------------- 213  P4_{1}32 #213
    spacegroupnumber=213;}
  else if(spacegroupsymbol=="I4_{1}32" || spacegroupsymbol=="I4_132" || spacegroupsymbol=="I4132") {  // ------------------- 214  I4_{1}32 #214
    spacegroupnumber=214;}
  else if(spacegroupsymbol=="P-43m") {  // ------------------- 215  P-43m #215
    spacegroupnumber=215;}
  else if(spacegroupsymbol=="F-43m") {  // ------------------- 216  F-43m #216
    spacegroupnumber=216;}
  else if(spacegroupsymbol=="I-43m") {  // ------------------- 217  I-43m #217
    spacegroupnumber=217;}
  else if(spacegroupsymbol=="P-43n") {  // ------------------- 218  P-43n #218
    spacegroupnumber=218;}
  else if(spacegroupsymbol=="F-43c") {  // ------------------- 219  F-43c #219
    spacegroupnumber=219;}
  else if(spacegroupsymbol=="I-43d") {  // ------------------- 220  I-43d #220
    spacegroupnumber=220;}
  else if(spacegroupsymbol=="Pm-3m") {  // ------------------- 221  Pm-3m #221
    spacegroupnumber=221;}
  else if(spacegroupsymbol=="Pn-3n") {  // ------------------- 222  Pn-3n #222
    spacegroupnumber=222;}
  else if(spacegroupsymbol=="Pm-3n") {  // ------------------- 223  Pm-3n #223
    spacegroupnumber=223;}
  else if(spacegroupsymbol=="Pn-3m") {  // ------------------- 224  Pn-3m #224
    spacegroupnumber=224;}
  else if(spacegroupsymbol=="Fm-3m") {  // ------------------- 225  Fm-3m #225
    spacegroupnumber=225;}
  else if(spacegroupsymbol=="Fm-3c") {  // ------------------- 226  Fm-3c #226
    spacegroupnumber=226;}
  else if(spacegroupsymbol=="Fd-3m") {  // ------------------- 227  Fd-3m #227
    spacegroupnumber=227;}
  else if(spacegroupsymbol=="Fd-3c") {  // ------------------- 228  Fd-3c #228
    spacegroupnumber=228;}
  else if(spacegroupsymbol=="Im-3m") {  // ------------------- 229  Im-3m #229
    spacegroupnumber=229;}
  else if(spacegroupsymbol=="Ia-3d") {  // ------------------- 230  Ia-3d #230
    spacegroupnumber=230;}
  // done
  else{
    message << "routine: space group specified invalid; perhaps non-ITC setting: ";
    message << "space group symbol=" << spacegroupsymbol << " [dir=" << directory << "].";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
  }
  return spacegroupnumber;
}

// ***************************************************************************
// GetSpaceGroupSchoenflies
// ***************************************************************************
string GetSpaceGroupSchoenflies(int spacegroupnumber, string directory) {
  if(spacegroupnumber < 1 || spacegroupnumber > 230) { //DX20190708 - for xerror
    stringstream message; //DX20190708 - for xerror
    message << "routine: space group specified invalid (1-230): "; //DX20190708 - for xerror
    message << spacegroupnumber << " [dir=" << directory << "]." << endl; //DX20190708 - for xerror
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_); //DX20190708 - for xerror
  }
  // OK
  //ME20210210 - converted to switch statement (faster)
  //DX+ME20190708 - changed subsequent "if" to "else if" -> efficiency
  switch(spacegroupnumber) {
    case 1: // ------------------- 1  C_{1}^{1} #1
      return "C_{1}^{1}";
    case 2: // ------------------- 2  C_{i}^{1} #2
      return "C_{i}^{1}";
    case 3: // ------------------- 3  C_{2}^{1} #3
      return "C_{2}^{1}";
    case 4: // ------------------- 4  C_{2}^{2} #4
      return "C_{2}^{2}";
    case 5: // ------------------- 5  C_{2}^{3} #5
      return "C_{2}^{3}";
    case 6: // ------------------- 6  C_{s}^{1} #6
      return "C_{s}^{1}";
    case 7: // ------------------- 7  C_{s}^{2} #7
      return "C_{s}^{2}";
    case 8: // ------------------- 8  C_{s}^{3} #8
      return "C_{s}^{3}";
    case 9: // ------------------- 9  C_{s}^{4} #9
      return "C_{s}^{4}";
    case 10: // ------------------- 10  C_{2h}^{1} #10
      return "C_{2h}^{1}";
    case 11: // ------------------- 11  C_{2h}^{2} #11
      return "C_{2h}^{2}";
    case 12: // ------------------- 12  C_{2h}^{3} #12
      return "C_{2h}^{3}";
    case 13: // ------------------- 13  C_{2h}^{4} #13
      return "C_{2h}^{4}";
    case 14: // ------------------- 14  C_{2h}^{5} #14
      return "C_{2h}^{5}";
    case 15: // ------------------- 15  C_{2h}^{6} #15
      return "C_{2h}^{6}";
    case 16: // ------------------- 16  D_{2}^{1} #16
      return "D_{2}^{1}";
    case 17: // ------------------- 17  D_{2}^{2} #17
      return "D_{2}^{2}";
    case 18: // ------------------- 18  D_{2}^{3} #18
      return "D_{2}^{3}";
    case 19: // ------------------- 19  D_{2}^{4} #19
      return "D_{2}^{4}";
    case 20: // ------------------- 20  D_{2}^{5} #20
      return "D_{2}^{5}";
    case 21: // ------------------- 21  D_{2}^{6} #21
      return "D_{2}^{6}";
    case 22: // ------------------- 22  D_{2}^{7} #22
      return "D_{2}^{7}";
    case 23: // ------------------- 23  D_{2}^{8} #23
      return "D_{2}^{8}";
    case 24: // ------------------- 24  D_{2}^{9} #24
      return "D_{2}^{9}";
    case 25: // ------------------- 25  C_{2v}^{1} #25
      return "C_{2v}^{1}";
    case 26: // ------------------- 26  C_{2v}^{2} #26
      return "C_{2v}^{2}";
    case 27: // ------------------- 27  C_{2v}^{3} #27
      return "C_{2v}^{3}";
    case 28: // ------------------- 28  C_{2v}^{4} #28
      return "C_{2v}^{4}";
    case 29: // ------------------- 29  C_{2v}^{5} #29
      return "C_{2v}^{5}";
    case 30: // ------------------- 30  C_{2v}^{6} #30
      return "C_{2v}^{6}";
    case 31: // ------------------- 31  C_{2v}^{7} #31
      return "C_{2v}^{7}";
    case 32: // ------------------- 32  C_{2v}^{8} #32
      return "C_{2v}^{8}";
    case 33: // ------------------- 33  C_{2v}^{9} #33
      return "C_{2v}^{9}";
    case 34: // ------------------- 34  C_{2v}^{10} #34
      return "C_{2v}^{10}";
    case 35: // ------------------- 35  C_{2v}^{11} #35
      return "C_{2v}^{11}";
    case 36: // ------------------- 36  C_{2v}^{12} #36
      return "C_{2v}^{12}";
    case 37: // ------------------- 37  C_{2v}^{13} #37
      return "C_{2v}^{13}";
    case 38: // ------------------- 38  C_{2v}^{14} #38
      return "C_{2v}^{14}";
    case 39: // ------------------- 39  C_{2v}^{15} #39
      return "C_{2v}^{15}";
    case 40: // ------------------- 40  C_{2v}^{16} #40
      return "C_{2v}^{16}";
    case 41: // ------------------- 41  C_{2v}^{17} #41
      return "C_{2v}^{17}";
    case 42: // ------------------- 42  C_{2v}^{18} #42
      return "C_{2v}^{18}";
    case 43: // ------------------- 43  C_{2v}^{19} #43
      return "C_{2v}^{19}";
    case 44: // ------------------- 44  C_{2v}^{20} #44
      return "C_{2v}^{20}";
    case 45: // ------------------- 45  C_{2v}^{21} #45
      return "C_{2v}^{21}";
    case 46: // ------------------- 46  C_{2v}^{22} #46
      return "C_{2v}^{22}";
    case 47: // ------------------- 47  D_{2h}^{1} #47
      return "D_{2h}^{1}";
    case 48: // ------------------- 48  D_{2h}^{2} #48
      return "D_{2h}^{2}";
    case 49: // ------------------- 49  D_{2h}^{3} #49
      return "D_{2h}^{3}";
    case 50: // ------------------- 50  D_{2h}^{4} #50
      return "D_{2h}^{4}";
    case 51: // ------------------- 51  D_{2h}^{5} #51
      return "D_{2h}^{5}";
    case 52: // ------------------- 52  D_{2h}^{6} #52
      return "D_{2h}^{6}";
    case 53: // ------------------- 53  D_{2h}^{7} #53
      return "D_{2h}^{7}";
    case 54: // ------------------- 54  D_{2h}^{8} #54
      return "D_{2h}^{8}";
    case 55: // ------------------- 55  D_{2h}^{9} #55
      return "D_{2h}^{9}";
    case 56: // ------------------- 56  D_{2h}^{10} #56
      return "D_{2h}^{10}";
    case 57: // ------------------- 57  D_{2h}^{11} #57
      return "D_{2h}^{11}";
    case 58: // ------------------- 58  D_{2h}^{12} #58
      return "D_{2h}^{12}";
    case 59: // ------------------- 59  D_{2h}^{13} #59
      return "D_{2h}^{13}";
    case 60: // ------------------- 60  D_{2h}^{14} #60
      return "D_{2h}^{14}";
    case 61: // ------------------- 61  D_{2h}^{15} #61
      return "D_{2h}^{15}";
    case 62: // ------------------- 62  D_{2h}^{16} #62
      return "D_{2h}^{16}";
    case 63: // ------------------- 63  D_{2h}^{17} #63
      return "D_{2h}^{17}";
    case 64: // ------------------- 64  D_{2h}^{18} #64
      return "D_{2h}^{18}";
    case 65: // ------------------- 65  D_{2h}^{19} #65
      return "D_{2h}^{19}";
    case 66: // ------------------- 66  D_{2h}^{20} #66
      return "D_{2h}^{20}";
    case 67: // ------------------- 67  D_{2h}^{21} #67
      return "D_{2h}^{21}";
    case 68: // ------------------- 68  D_{2h}^{22} #68
      return "D_{2h}^{22}";
    case 69: // ------------------- 69  D_{2h}^{23} #69
      return "D_{2h}^{23}";
    case 70: // ------------------- 70  D_{2h}^{24} #70
      return "D_{2h}^{24}";
    case 71: // ------------------- 71  D_{2h}^{25} #71
      return "D_{2h}^{25}";
    case 72: // ------------------- 72  D_{2h}^{26} #72
      return "D_{2h}^{26}";
    case 73: // ------------------- 73  D_{2h}^{27} #73
      return "D_{2h}^{27}";
    case 74: // ------------------- 74  D_{2h}^{28} #74
      return "D_{2h}^{28}";
    case 75: // ------------------- 75  C_{4}^{1} #75
      return "C_{4}^{1}";
    case 76: // ------------------- 76  C_{4}^{2} #76
      return "C_{4}^{2}";
    case 77: // ------------------- 77  C_{4}^{3} #77
      return "C_{4}^{3}";
    case 78: // ------------------- 78  C_{4}^{4} #78
      return "C_{4}^{4}";
    case 79: // ------------------- 79  C_{4}^{5} #79
      return "C_{4}^{5}";
    case 80: // ------------------- 80  C_{4}^{6} #80
      return "C_{4}^{6}";
    case 81: // ------------------- 81  S_{4}^{1} #81
      return "S_{4}^{1}";
    case 82: // ------------------- 82  S_{4}^{2} #82
      return "S_{4}^{2}";
    case 83: // ------------------- 83  C_{4h}^{1} #83
      return "C_{4h}^{1}";
    case 84: // ------------------- 84  C_{4h}^{2} #84
      return "C_{4h}^{2}";
    case 85: // ------------------- 85  C_{4h}^{3} #85
      return "C_{4h}^{3}";
    case 86: // ------------------- 86  C_{4h}^{4} #86
      return "C_{4h}^{4}";
    case 87: // ------------------- 87  C_{4h}^{5} #87
      return "C_{4h}^{5}";
    case 88: // ------------------- 88  C_{4h}^{6} #88
      return "C_{4h}^{6}";
    case 89: // ------------------- 89  D_{4}^{1} #89
      return "D_{4}^{1}";
    case 90: // ------------------- 90  D_{4}^{2} #90
      return "D_{4}^{2}";
    case 91: // ------------------- 91  D_{4}^{3} #91
      return "D_{4}^{3}";
    case 92: // ------------------- 92  D_{4}^{4} #92
      return "D_{4}^{4}";
    case 93: // ------------------- 93  D_{4}^{5} #93
      return "D_{4}^{5}";
    case 94: // ------------------- 94  D_{4}^{6} #94
      return "D_{4}^{6}";
    case 95: // ------------------- 95  D_{4}^{7} #95
      return "D_{4}^{7}";
    case 96: // ------------------- 96  D_{4}^{8} #96
      return "D_{4}^{8}";
    case 97: // ------------------- 97  D_{4}^{9} #97
      return "D_{4}^{9}";
    case 98: // ------------------- 98  D_{4}^{10} #98
      return "D_{4}^{10}";
    case 99: // ------------------- 99  C_{4v}^{1} #99
      return "C_{4v}^{1}";
    case 100: // ------------------- 100  C_{4v}^{2} #100
      return "C_{4v}^{2}";
    case 101: // ------------------- 101  C_{4v}^{3} #101
      return "C_{4v}^{3}";
    case 102: // ------------------- 102  C_{4v}^{4} #102
      return "C_{4v}^{4}";
    case 103: // ------------------- 103  C_{4v}^{5} #103
      return "C_{4v}^{5}";
    case 104: // ------------------- 104  C_{4v}^{6} #104
      return "C_{4v}^{6}";
    case 105: // ------------------- 105  C_{4v}^{7} #105
      return "C_{4v}^{7}";
    case 106: // ------------------- 106  C_{4v}^{8} #106
      return "C_{4v}^{8}";
    case 107: // ------------------- 107  C_{4v}^{9} #107
      return "C_{4v}^{9}";
    case 108: // ------------------- 108  C_{4v}^{10} #108
      return "C_{4v}^{10}";
    case 109: // ------------------- 109  C_{4v}^{11} #109
      return "C_{4v}^{11}";
    case 110: // ------------------- 110  C_{4v}^{12} #110
      return "C_{4v}^{12}";
    case 111: // ------------------- 111  D_{2d}^{1} #111
      return "D_{2d}^{1}";
    case 112: // ------------------- 112  D_{2d}^{2} #112
      return "D_{2d}^{2}";
    case 113: // ------------------- 113  D_{2d}^{3} #113
      return "D_{2d}^{3}";
    case 114: // ------------------- 114  D_{2d}^{4} #114
      return "D_{2d}^{4}";
    case 115: // ------------------- 115  D_{2d}^{5} #115
      return "D_{2d}^{5}";
    case 116: // ------------------- 116  D_{2d}^{6} #116
      return "D_{2d}^{6}";
    case 117: // ------------------- 117  D_{2d}^{7} #117
      return "D_{2d}^{7}";
    case 118: // ------------------- 118  D_{2d}^{8} #118
      return "D_{2d}^{8}";
    case 119: // ------------------- 119  D_{2d}^{9} #119
      return "D_{2d}^{9}";
    case 120: // ------------------- 120  D_{2d}^{10} #120
      return "D_{2d}^{10}";
    case 121: // ------------------- 121  D_{2d}^{11} #121
      return "D_{2d}^{11}";
    case 122: // ------------------- 122  D_{2d}^{12} #122
      return "D_{2d}^{12}";
    case 123: // ------------------- 123  D_{4h}^{1} #123
      return "D_{4h}^{1}";
    case 124: // ------------------- 124  D_{4h}^{2} #124
      return "D_{4h}^{2}";
    case 125: // ------------------- 125  D_{4h}^{3} #125
      return "D_{4h}^{3}";
    case 126: // ------------------- 126  D_{4h}^{4} #126
      return "D_{4h}^{4}";
    case 127: // ------------------- 127  D_{4h}^{5} #127
      return "D_{4h}^{5}";
    case 128: // ------------------- 128  D_{4h}^{6} #128
      return "D_{4h}^{6}";
    case 129: // ------------------- 129  D_{4h}^{7} #129
      return "D_{4h}^{7}";
    case 130: // ------------------- 130  D_{4h}^{8} #130
      return "D_{4h}^{8}";
    case 131: // ------------------- 131  D_{4h}^{9} #131
      return "D_{4h}^{9}";
    case 132: // ------------------- 132  D_{4h}^{10} #132
      return "D_{4h}^{10}";
    case 133: // ------------------- 133  D_{4h}^{11} #133
      return "D_{4h}^{11}";
    case 134: // ------------------- 134  D_{4h}^{12} #134
      return "D_{4h}^{12}";
    case 135: // ------------------- 135  D_{4h}^{13} #135
      return "D_{4h}^{13}";
    case 136: // ------------------- 136  D_{4h}^{14} #136
      return "D_{4h}^{14}";
    case 137: // ------------------- 137  D_{4h}^{15} #137
      return "D_{4h}^{15}";
    case 138: // ------------------- 138  D_{4h}^{16} #138
      return "D_{4h}^{16}";
    case 139: // ------------------- 139  D_{4h}^{17} #139
      return "D_{4h}^{17}";
    case 140: // ------------------- 140  D_{4h}^{18} #140
      return "D_{4h}^{18}";
    case 141: // ------------------- 141  D_{4h}^{19} #141
      return "D_{4h}^{19}";
    case 142: // ------------------- 142  D_{4h}^{20} #142
      return "D_{4h}^{20}";
    case 143: // ------------------- 143  C_{3}^{1} #143
      return "C_{3}^{1}";
    case 144: // ------------------- 144  C_{3}^{2} #144
      return "C_{3}^{2}";
    case 145: // ------------------- 145  C_{3}^{3} #145
      return "C_{3}^{3}";
    case 146: // ------------------- 146  C_{3}^{4} #146
      return "C_{3}^{4}";
    case 147: // ------------------- 147  C_{3i}^{1} #147
      return "C_{3i}^{1}";
    case 148: // ------------------- 148  C_{3i}^{2} #148
      return "C_{3i}^{2}";
    case 149: // ------------------- 149  D_{3}^{1} #149
      return "D_{3}^{1}";
    case 150: // ------------------- 150  D_{3}^{2} #150
      return "D_{3}^{2}";
    case 151: // ------------------- 151  D_{3}^{3} #151
      return "D_{3}^{3}";
    case 152: // ------------------- 152  D_{3}^{4} #152
      return "D_{3}^{4}";
    case 153: // ------------------- 153  D_{3}^{5} #153
      return "D_{3}^{5}";
    case 154: // ------------------- 154  D_{3}^{6} #154
      return "D_{3}^{6}";
    case 155: // ------------------- 155  D_{3}^{7} #155
      return "D_{3}^{7}";
    case 156: // ------------------- 156  C_{3v}^{1} #156
      return "C_{3v}^{1}";
    case 157: // ------------------- 157  C_{3v}^{2} #157
      return "C_{3v}^{2}";
    case 158: // ------------------- 158  C_{3v}^{3} #158
      return "C_{3v}^{3}";
    case 159: // ------------------- 159  C_{3v}^{4} #159
      return "C_{3v}^{4}";
    case 160: // ------------------- 160  C_{3v}^{5} #160
      return "C_{3v}^{5}";
    case 161: // ------------------- 161  C_{3v}^{6} #161
      return "C_{3v}^{6}";
    case 162: // ------------------- 162  D_{3d}^{1} #162
      return "D_{3d}^{1}";
    case 163: // ------------------- 163  D_{3d}^{2} #163
      return "D_{3d}^{2}";
    case 164: // ------------------- 164  D_{3d}^{3} #164
      return "D_{3d}^{3}";
    case 165: // ------------------- 165  D_{3d}^{4} #165
      return "D_{3d}^{4}";
    case 166: // ------------------- 166  D_{3d}^{5} #166
      return "D_{3d}^{5}";
    case 167: // ------------------- 167  D_{3d}^{6} #167
      return "D_{3d}^{6}";
    case 168: // ------------------- 168  C_{6}^{1} #168
      return "C_{6}^{1}";
    case 169: // ------------------- 169  C_{6}^{2} #169
      return "C_{6}^{2}";
    case 170: // ------------------- 170  C_{6}^{3} #170
      return "C_{6}^{3}";
    case 171: // ------------------- 171  C_{6}^{4} #171
      return "C_{6}^{4}";
    case 172: // ------------------- 172  C_{6}^{5} #172
      return "C_{6}^{5}";
    case 173: // ------------------- 173  C_{6}^{6} #173
      return "C_{6}^{6}";
    case 174: // ------------------- 174  C_{3h}^{1} #174
      return "C_{3h}^{1}";
    case 175: // ------------------- 175  C_{6h}^{1} #175
      return "C_{6h}^{1}";
    case 176: // ------------------- 176  C_{6h}^{2} #176
      return "C_{6h}^{2}";
    case 177: // ------------------- 177  D_{6}^{1} #177
      return "D_{6}^{1}";
    case 178: // ------------------- 178  D_{6}^{2} #178
      return "D_{6}^{2}";
    case 179: // ------------------- 179  D_{6}^{3} #179
      return "D_{6}^{3}";
    case 180: // ------------------- 180  D_{6}^{4} #180
      return "D_{6}^{4}";
    case 181: // ------------------- 181  D_{6}^{5} #181
      return "D_{6}^{5}";
    case 182: // ------------------- 182  D_{6}^{6} #182
      return "D_{6}^{6}";
    case 183: // ------------------- 183  C_{6v}^{1} #183
      return "C_{6v}^{1}";
    case 184: // ------------------- 184  C_{6v}^{2} #184
      return "C_{6v}^{2}";
    case 185: // ------------------- 185  C_{6v}^{3} #185
      return "C_{6v}^{3}";
    case 186: // ------------------- 186  C_{6v}^{4} #186
      return "C_{6v}^{4}";
    case 187: // ------------------- 187  D_{3h}^{1} #187
      return "D_{3h}^{1}";
    case 188: // ------------------- 188  D_{3h}^{2} #188
      return "D_{3h}^{2}";
    case 189: // ------------------- 189  D_{3h}^{3} #189
      return "D_{3h}^{3}";
    case 190: // ------------------- 190  D_{3h}^{4} #190
      return "D_{3h}^{4}";
    case 191: // ------------------- 191  D_{6h}^{1} #191
      return "D_{6h}^{1}";
    case 192: // ------------------- 192  D_{6h}^{2} #192
      return "D_{6h}^{2}";
    case 193: // ------------------- 193  D_{6h}^{3} #193
      return "D_{6h}^{3}";
    case 194: // ------------------- 194  D_{6h}^{4} #194
      return "D_{6h}^{4}";
    case 195: // ------------------- 195  T^{1} #195
      return "T^{1}";
    case 196: // ------------------- 196  T^{2} #196
      return "T^{2}";
    case 197: // ------------------- 197  T^{3} #197
      return "T^{3}";
    case 198: // ------------------- 198  T^{4} #198
      return "T^{4}";
    case 199: // ------------------- 199  T^{5} #199
      return "T^{5}";
    case 200: // ------------------- 200  T_{h}^{1} #200
      return "T_{h}^{1}";
    case 201: // ------------------- 201  T_{h}^{2} #201
      return "T_{h}^{2}";
    case 202: // ------------------- 202  T_{h}^{3} #202
      return "T_{h}^{3}";
    case 203: // ------------------- 203  T_{h}^{4} #203
      return "T_{h}^{4}";
    case 204: // ------------------- 204  T_{h}^{5} #204
      return "T_{h}^{5}";
    case 205: // ------------------- 205  T_{h}^{6} #205
      return "T_{h}^{6}";
    case 206: // ------------------- 206  T_{h}^{7} #206
      return "T_{h}^{7}";
    case 207: // ------------------- 207  O^{1} #207
      return "O^{1}";
    case 208: // ------------------- 208  O^{2} #208
      return "O^{2}";
    case 209: // ------------------- 209  O^{3} #209
      return "O^{3}";
    case 210: // ------------------- 210  O^{4} #210
      return "O^{4}";
    case 211: // ------------------- 211  O^{5} #211
      return "O^{5}";
    case 212: // ------------------- 212  O^{6} #212
      return "O^{6}";
    case 213: // ------------------- 213  O^{7} #213
      return "O^{7}";
    case 214: // ------------------- 214  O^{8} #214
      return "O^{8}";
    case 215: // ------------------- 215  T_{d}^{1} #215
      return "T_{d}^{1}";
    case 216: // ------------------- 216  T_{d}^{2} #216
      return "T_{d}^{2}";
    case 217: // ------------------- 217  T_{d}^{3} #217
      return "T_{d}^{3}";
    case 218: // ------------------- 218  T_{d}^{4} #218
      return "T_{d}^{4}";
    case 219: // ------------------- 219  T_{d}^{5} #219
      return "T_{d}^{5}";
    case 220: // ------------------- 220  T_{d}^{6} #220
      return "T_{d}^{6}";
    case 221: // ------------------- 221  O_{h}^{1} #221
      return "O_{h}^{1}";
    case 222: // ------------------- 222  O_{h}^{2} #222
      return "O_{h}^{2}";
    case 223: // ------------------- 223  O_{h}^{3} #223
      return "O_{h}^{3}";
    case 224: // ------------------- 224  O_{h}^{4} #224
      return "O_{h}^{4}";
    case 225: // ------------------- 225  O_{h}^{5} #225
      return "O_{h}^{5}";
    case 226: // ------------------- 226  O_{h}^{6} #226
      return "O_{h}^{6}";
    case 227: // ------------------- 227  O_{h}^{7} #227
      return "O_{h}^{7}";
    case 228: // ------------------- 228  O_{h}^{8} #228
      return "O_{h}^{8}";
    case 229: // ------------------- 229  O_{h}^{9} #229
      return "O_{h}^{9}";
    case 230: // ------------------- 230  O_{h}^{10} #230
      return "O_{h}^{10}";
  }
  // No space group found
  return "";
}

// ***************************************************************************
// GetSpaceGroupHall
// ***************************************************************************
string GetSpaceGroupHall(int spacegroupnumber, int setting, string directory) {
  //DX - Hall distinguishes space group setting.  This table assumes the first 
  //      setting that appears in the ITC.
  //      For more settings, they need to be hard-coded here.
  stringstream message; //DX20190708 - for xerror
  string spacegroup=""; //DX20190708 - for xerror
  if(spacegroupnumber < 1 || spacegroupnumber > 230) { //DX20190708 - for xerror
    message << "routine: space group specified invalid (1-230): "; //DX20190708 - for xerror
    message << spacegroupnumber << " [dir=" << directory << "]." << endl; //DX20190708 - for xerror
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_); //DX20190708 - for xerror
  }
  // OK
  if(setting==SG_SETTING_ANRL){ setting = anrl::getANRLSettingChoice(spacegroupnumber); } //DX20210420
  else if (setting==0){ //signals default //DX20180807
    // if RHL, AFLOW prefers hexagonal setting (i.e., setting=2)
    if(spacegroupnumber==146 || spacegroupnumber==148 || spacegroupnumber==155 || spacegroupnumber==160 || 
        spacegroupnumber==161 || spacegroupnumber==166 || spacegroupnumber==167){
      setting=2;
    }
    // else setting==SG_SETTING_1
    else {
      setting=SG_SETTING_1;
    }
  }
  if(setting < 1 || setting > 2) {
    message << "routine: setting choice is invalid (1 or 2 only): " << setting << " [dir=" << directory << "]."; //DX20190708 - for xerror
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_); //DX20190708 - for xerror
  }
  // OK
  //DX+ME20190708 - changed subsequent "if" to "else if" -> efficiency
  //ME20210210 - converted to switch statement (faster)
  switch (spacegroupnumber) {
    case 1: // ------------------- 1  P 1 #1
      return "P 1";
    case 2: // ------------------- 2  -P 1 #2
      return "-P 1";
    case 3: // ------------------- 3  P 2y #3
      return "P 2y";
    case 4: // ------------------- 4  P 2yb #4
      return "P 2yb";
    case 5: // ------------------- 5  C 2y #5
      return "C 2y";
    case 6: // ------------------- 6  P -2y #6
      return "P -2y";
    case 7: // ------------------- 7  P -2yc #7
      return "P -2yc";
    case 8: // ------------------- 8  C -2y #8
      return "C -2y";
    case 9: // ------------------- 9  C -2yc #9
      return "C -2yc";
    case 10: // ------------------- 10  -P 2y #10
      return "-P 2y";
    case 11: // ------------------- 11  -P 2yb #11
      return "-P 2yb";
    case 12: // ------------------- 12  -C 2y #12
      return "-C 2y";
    case 13: // ------------------- 13  -P 2yc #13
      return "-P 2yc";
    case 14: // ------------------- 14  -P 2ybc #14
      return "-P 2ybc";
    case 15: // ------------------- 15  -C 2yc #15
      return "-C 2yc";
    case 16: // ------------------- 16  P 2 2 #16
      return "P 2 2";
    case 17: // ------------------- 17  P 2c 2 #17
      return "P 2c 2";
    case 18: // ------------------- 18  P 2 2ab #18
      return "P 2 2ab";
    case 19: // ------------------- 19  P 2ac 2ab #19
      return "P 2ac 2ab";
    case 20: // ------------------- 20  C 2c 2 #20
      return "C 2c 2";
    case 21: // ------------------- 21  C 2 2 #21
      return "C 2 2";
    case 22: // ------------------- 22  F 2 2 #22
      return "F 2 2";
    case 23: // ------------------- 23  I 2 2 #23
      return "I 2 2";
    case 24: // ------------------- 24  I 2b 2c #24
      return "I 2b 2c";
    case 25: // ------------------- 25  P 2 -2 #25
      return "P 2 -2";
    case 26: // ------------------- 26  P 2c -2 #26
      return "P 2c -2";
    case 27: // ------------------- 27  P 2 -2c #27
      return "P 2 -2c";
    case 28: // ------------------- 28  P 2 -2a #28
      return "P 2 -2a";
    case 29: // ------------------- 29  P 2c -2ac #29
      return "P 2c -2ac";
    case 30: // ------------------- 30  P 2 -2bc #30
      return "P 2 -2bc";
    case 31: // ------------------- 31  P 2ac -2 #31
      return "P 2ac -2";
    case 32: // ------------------- 32  P 2 -2ab #32
      return "P 2 -2ab";
    case 33: // ------------------- 33  P 2c -2n #33
      return "P 2c -2n";
    case 34: // ------------------- 34  P 2 -2n #34
      return "P 2 -2n";
    case 35: // ------------------- 35  C 2 -2 #35
      return "C 2 -2";
    case 36: // ------------------- 36  C 2c -2 #36
      return "C 2c -2";
    case 37: // ------------------- 37  C 2 -2c #37
      return "C 2 -2c";
    case 38: // ------------------- 38  A 2 -2 #38
      return "A 2 -2";
    case 39: // ------------------- 39  A 2 -2c #39
      return "A 2 -2c";
    case 40: // ------------------- 40  A 2 -2a #40
      return "A 2 -2a";
    case 41: // ------------------- 41  A 2 -2ac #41
      return "A 2 -2ac";
    case 42: // ------------------- 42  F 2 -2 #42
      return "F 2 -2";
    case 43: // ------------------- 43  F 2 -2d #43
      return "F 2 -2d";
    case 44: // ------------------- 44  I 2 -2 #44
      return "I 2 -2";
    case 45: // ------------------- 45  I 2 -2c #45
      return "I 2 -2c";
    case 46: // ------------------- 46  I 2 -2a #46
      return "I 2 -2a";
    case 47: // ------------------- 47  -P 2 2 #47
      return "-P 2 2";
    case 48:
      if (setting==SG_SETTING_1) return "P 2 2 -1n";  // ------------------- 48  -P P 2 2 -1n #48 (setting 1)
      return "-P 2ab 2bc";                            // ------------------- 48  -P 2ab 2bc #48 (setting 2)
    case 49: // ------------------- 49  -P 2 2c #49
      return "-P 2 2c";
    case 50:
      if (setting==SG_SETTING_1) return "P 2 2 -1ab"; // ------------------- 50  -P P 2 2 -1ab #50 (setting 1)
      else return "-P 2ab 2b";             // ------------------- 50  -P 2ab 2b #50 (setting 2)
    case 51: // ------------------- 51  -P 2a 2a #51
      return "-P 2a 2a";
    case 52: // ------------------- 52  -P 2a 2bc #52
      return "-P 2a 2bc";
    case 53: // ------------------- 53  -P 2ac 2 #53
      return "-P 2ac 2";
    case 54: // ------------------- 54  -P 2a 2ac #54
      return "-P 2a 2ac";
    case 55: // ------------------- 55  -P 2 2ab #55
      return "-P 2 2ab";
    case 56: // ------------------- 56  -P 2ab 2ac #56
      return "-P 2ab 2ac";
    case 57: // ------------------- 57  -P 2c 2b #57
      return "-P 2c 2b";
    case 58: // ------------------- 58  -P 2 2n #58
      return "-P 2 2n";
    case 59:
      if (setting==SG_SETTING_1) return "P 2 2ab -1ab"; // ------------------- 59  P 2 2ab -1ab #59 (setting 1)
      else return "-P 2ab 2a";                          // ------------------- 59  -P 2ab 2a #59 (setting 2)
    case 60: // ------------------- 60  -P 2n 2ab #60
      return "-P 2n 2ab";
    case 61: // ------------------- 61  -P 2ac 2ab #61
      return "-P 2ac 2ab";
    case 62: // ------------------- 62  -P 2ac 2n #62
      return "-P 2ac 2n";
    case 63: // ------------------- 63  -C 2c 2 #63
      return "-C 2c 2";
    case 64: // ------------------- 64  -C 2bc 2 #64
      return "-C 2bc 2";
    case 65: // ------------------- 65  -C 2 2 #65
      return "-C 2 2";
    case 66: // ------------------- 66  -C 2 2c #66
      return "-C 2 2c";
    case 67: // ------------------- 67  -C 2b 2 #67
      return "-C 2b 2";
    case 68:
      if (setting==SG_SETTING_1) return "C 2 2 -1bc"; // ------------------- 68  C 2 2 -1bc #68 (setting 1)
      else return "-C 2b 2bc";                        // ------------------- 68  -C 2b 2bc #68 (setting 2)
    case 69: // ------------------- 69  -F 2 2 #69
      return "-F 2 2";
    case 70:
      if (setting==SG_SETTING_1) return "F 2 2 -1d"; // ------------------- 70  F 2 2 -1d #70 (setting 1)
      else return "-F 2uv 2vw";                      // ------------------- 70  -F 2uv 2vw #70 (setting 2)
    case 71: // ------------------- 71  -I 2 2 #71
      return "-I 2 2";
    case 72: // ------------------- 72  -I 2 2c #72
      return "-I 2 2c";
    case 73: // ------------------- 73  -I 2b 2c #73
      return "-I 2b 2c";
    case 74: // ------------------- 74  -I 2b 2 #74
      return "-I 2b 2";
    case 75: // ------------------- 75  P 4 #75
      return "P 4";
    case 76: // ------------------- 76  P 4w #76
      return "P 4w";
    case 77: // ------------------- 77  P 4c #77
      return "P 4c";
    case 78: // ------------------- 78  P 4cw #78
      return "P 4cw";
    case 79: // ------------------- 79  I 4 #79
      return "I 4";
    case 80: // ------------------- 80  I 4bw #80
      return "I 4bw";
    case 81: // ------------------- 81  P -4 #81
      return "P -4";
    case 82: // ------------------- 82  I -4 #82
      return "I -4";
    case 83: // ------------------- 83  -P 4 #83
      return "-P 4";
    case 84: // ------------------- 84  -P 4c #84
      return "-P 4c";
    case 85:
      if (setting==SG_SETTING_1) return "P 4ab -1ab"; // ------------------- 85  P 4ab -1ab #85 (setting 1)
      else return "-P 4a";                            // ------------------- 85  -P 4a #85 (setting 2)
    case 86:
      if (setting==SG_SETTING_1) return "P 4n -1n"; // ------------------- 86  P 4n -1n #86 (setting 1)
      else return "-P 4bc";                         // ------------------- 86  -P 4bc #86 (setting 2)
    case 87: // ------------------- 87  -I 4 #87
      return "-I 4";
    case 88:
      if (setting==SG_SETTING_1) return "I 4bw -1bw"; // ------------------- 88  I 4bw -1bw #88 (setting 1)
      else return "-I 4ad";                           // ------------------- 88  -I 4ad #88 (setting 2)
    case 89: // ------------------- 89  P 4 2 #89
      return "P 4 2";
    case 90: // ------------------- 90  P 4ab 2ab #90
      return "P 4ab 2ab";
    case 91: // ------------------- 91  P 4w 2c #91
      return "P 4w 2c";
    case 92: // ------------------- 92  P 4abw 2nw #92
      return "P 4abw 2nw";
    case 93: // ------------------- 93  P 4c 2 #93
      return "P 4c 2";
    case 94: // ------------------- 94  P 4n 2n #94
      return "P 4n 2n";
    case 95: // ------------------- 95  P 4cw 2c #95
      return "P 4cw 2c";
    case 96: // ------------------- 96  P 4nw 2abw #96
      return "P 4nw 2abw";
    case 97: // ------------------- 97  I 4 2 #97
      return "I 4 2";
    case 98: // ------------------- 98  I 4bw 2bw #98
      return "I 4bw 2bw";
    case 99: // ------------------- 99  P 4 -2 #99
      return "P 4 -2";
    case 100: // ------------------- 100  P 4 -2ab #100
      return "P 4 -2ab";
    case 101: // ------------------- 101  P 4c -2c #101
      return "P 4c -2c";
    case 102: // ------------------- 102  P 4n -2n #102
      return "P 4n -2n";
    case 103: // ------------------- 103  P 4 -2c #103
      return "P 4 -2c";
    case 104: // ------------------- 104  P 4 -2n #104
      return "P 4 -2n";
    case 105: // ------------------- 105  P 4c -2 #105
      return "P 4c -2";
    case 106: // ------------------- 106  P 4c -2ab #106
      return "P 4c -2ab";
    case 107: // ------------------- 107  I 4 -2 #107
      return "I 4 -2";
    case 108: // ------------------- 108  I 4 -2c #108
      return "I 4 -2c";
    case 109: // ------------------- 109  I 4bw -2 #109
      return "I 4bw -2";
    case 110: // ------------------- 110  I 4bw -2c #110
      return "I 4bw -2c";
    case 111: // ------------------- 111  P -4 2 #111
      return "P -4 2";
    case 112: // ------------------- 112  P -4 2c #112
      return "P -4 2c";
    case 113: // ------------------- 113  P -4 2ab #113
      return "P -4 2ab";
    case 114: // ------------------- 114  P -4 2n #114
      return "P -4 2n";
    case 115: // ------------------- 115  P -4 -2 #115
      return "P -4 -2";
    case 116: // ------------------- 116  P -4 -2c #116
      return "P -4 -2c";
    case 117: // ------------------- 117  P -4 -2ab #117
      return "P -4 -2ab";
    case 118: // ------------------- 118  P -4 -2n #118
      return "P -4 -2n";
    case 119: // ------------------- 119  I -4 -2 #119
      return "I -4 -2";
    case 120: // ------------------- 120  I -4 -2c #120
      return "I -4 -2c";
    case 121: // ------------------- 121  I -4 2 #121
      return "I -4 2";
    case 122: // ------------------- 122  I -4 2bw #122
      return "I -4 2bw";
    case 123: // ------------------- 123  -P 4 2 #123
      return "-P 4 2";
    case 124: // ------------------- 124  -P 4 2c #124
      return "-P 4 2c";
    case 125:
      if (setting==SG_SETTING_1) return "P 4 2 -1ab"; // ------------------- 125  P 4 2 -1ab #125 (setting 1)
      else return "-P 4a 2b";                         // ------------------- 125  -P 4a 2b #125 (setting 2)
    case 126:
      if (setting==SG_SETTING_1) return "P 4 2 -1n"; // ------------------- 126  P 4 2 -1n #126 (setting 1)
      else return "-P 4a 2bc";                       // ------------------- 126  -P 4a 2bc #126 (setting 2)
    case 127: // ------------------- 127  -P 4 2ab #127
      return "-P 4 2ab";
    case 128: // ------------------- 128  -P 4 2n #128
      return "-P 4 2n";
    case 129:
      if (setting==SG_SETTING_1) return "P 4ab 2ab -1ab"; // ------------------- 129  P 4ab 2ab -1ab #129 (setting 1)
      else return "-P 4a 2a";                             // ------------------- 129  -P 4a 2a #129 (setting 2)
    case 130:
      if (setting==SG_SETTING_1) return "P 4ab 2n -1ab"; // ------------------- 130  P 4ab 2n -1ab #130 (setting 1)
      else return "-P 4a 2ac";                           // ------------------- 130  -P 4a 2ac #130 (setting 2)
    case 131: // ------------------- 131  -P 4c 2 #131
      return "-P 4c 2";
    case 132: // ------------------- 132  -P 4c 2c #132
      return "-P 4c 2c";
    case 133:
      if (setting==SG_SETTING_1) return "P 4n 2c -1n"; // ------------------- 133  P 4n 2c -1n #133 (setting 1)
      else return "-P 4ac 2b";                         // ------------------- 133  -P 4ac 2b #133 (setting 2)
    case 134:
      if (setting==SG_SETTING_1) return "P 4n 2 -1n"; // ------------------- 134  -P 4ac 2bc #134 (setting 1)
      else return "-P 4ac 2bc";                       // ------------------- 134  -P 4ac 2bc #134 (setting 2)
    case 135: // ------------------- 135  -P 4c 2ab #135
      return "-P 4c 2ab";
    case 136: // ------------------- 136  -P 4n 2n #136
      return "-P 4n 2n";
    case 137:
      if (setting==SG_SETTING_1) return "P 4n 2n -1n"; // ------------------- 137  P 4n 2n -1n #137 (setting 1)
      else return "-P 4ac 2a";                         // ------------------- 137  -P 4ac 2a #137 (setting 2)
    case 138:
      if (setting==SG_SETTING_1) return "P 4n 2ab -1n"; // ------------------- 138  P 4n 2ab -1n #138 (setting 1)
      else return "-P 4ac 2ac";                         // ------------------- 138  -P 4ac 2ac #138 (setting 2)
    case 139: // ------------------- 139  -I 4 2 #139
      return "-I 4 2";
    case 140: // ------------------- 140  -I 4 2c #140
      return "-I 4 2c";
    case 141:
      if (setting==SG_SETTING_1) return "I 4bw 2bw -1bw"; // ------------------- 141  I 4bw 2bw -1bw #141 (setting 1)
      else return "-I 4bd 2";                             // ------------------- 141  -I 4bd 2 #141 (setting 2)
    case 142:
      if (setting==SG_SETTING_1) return "I 4bw 2aw -1bw"; // ------------------- 142  I 4bw 2aw -1bw #142 (setting 1)
      else return "-I 4bd 2c"; // ------------------- 142  -I 4bd 2c #142 (setting 2)
    case 143: // ------------------- 143  P 3 #143
      return "P 3";
    case 144: // ------------------- 144  P 31 #144
      return "P 31";
    case 145: // ------------------- 145  P 32 #145
      return "P 32";
    case 146:
      if (setting==SG_SETTING_1) return "P 3*"; // ------------------- 146  P 3* #146 (setting 1, rhl)
      else return "R 3";                        // ------------------- 146  R 3 #146 (setting 2, hex)
    case 147: // ------------------- 147  -P 3 #147
      return "-P 3";
    case 148:
      if (setting==SG_SETTING_1) return "-P 3*"; // ------------------- 148  -P 3* #148 (setting 1, rhl)
      else return "-R 3";                        // ------------------- 148  -R 3 #148 (setting 2, hex)
    case 149: // ------------------- 149  P 3 2 #149
      return "P 3 2";
    case 150: // ------------------- 150  P 3 2'' #150
      return "P 3 2''";
    case 151: // ------------------- 151  P 31 2c (0 0 1) #151
      return "P 31 2c (0 0 1)";
    case 152: // ------------------- 152  P 31 2'' #152
      return "P 31 2''";
    case 153: // ------------------- 153  P 32 2c (0 0 -1) #153
      return "P 32 2c (0 0 -1)";
    case 154: // ------------------- 154  P 32 2'' #154
      return "P 32 2''";
    case 155:
      if (setting==SG_SETTING_1) return "P 3* 2"; // ------------------- 155  P 3* 2 #155 (setting 1, rhl)
      else return "R 3 2''";                      // ------------------- 155  R 3 2'' #155 (setting 2, hex)
    case 156: // ------------------- 156  P 3 -2'' #156
      return "P 3 -2''";
    case 157: // ------------------- 157  P 3 -2 #157
      return "P 3 -2";
    case 158: // ------------------- 158  P 3 -2''c #158
      return "P 3 -2''c";
    case 159: // ------------------- 159  P 3 -2c #159
      return "P 3 -2c";
    case 160:
      if (setting==SG_SETTING_1) return "P 3* -2"; // ------------------- 160  P 3* -2 #160 (setting 1, rhl)
      else return "R 3 -2''";                      // ------------------- 160  R 3 -2'' #160 (setting 2, hex)
    case 161:
      if (setting==SG_SETTING_1) return "P 3* -2n"; // ------------------- 161  P 3* -2n #161 (setting 1, rhl)
      else return "R 3 -2''c";                      // ------------------- 161  R 3 -2''c #161 (setting 2, hex)
    case 162: // ------------------- 162  -P 3 2 #162
      return "-P 3 2";
    case 163: // ------------------- 163  -P 3 2c #163
      return "-P 3 2c";
    case 164: // ------------------- 164  -P 3 2'' #164
      return "-P 3 2''";
    case 165: // ------------------- 165  -P 3 2''c #165
      return "-P 3 2''c";
    case 166:
      if (setting==SG_SETTING_1) return "-P 3* 2"; // ------------------- 166  -P 3* 2 #166 (setting 1, rhl)
      else return "-R 3 2''";                      // ------------------- 166  -R 3 2'' #166 (setting 2, hex)
    case 167:
      if (setting==SG_SETTING_1) return "-P 3* 2n"; // ------------------- 167  -P 3* 2n #167 (setting 1, rhl)
      else return "-R 3 2''c";                      // ------------------- 167  -R 3 2''c #167 (setting 2, hex)
    case 168: // ------------------- 168  P 6 #168
      return "P 6";
    case 169: // ------------------- 169  P 61 #169
      return "P 61";
    case 170: // ------------------- 170  P 65 #170
      return "P 65";
    case 171: // ------------------- 171  P 62 #171
      return "P 62";
    case 172: // ------------------- 172  P 64 #172
      return "P 64";
    case 173: // ------------------- 173  P 6c #173
      return "P 6c";
    case 174: // ------------------- 174  P -6 #174
      return "P -6";
    case 175: // ------------------- 175  -P 6 #175
      return "-P 6";
    case 176: // ------------------- 176  -P 6c #176
      return "-P 6c";
    case 177: // ------------------- 177  P 6 2 #177
      return "P 6 2";
    case 178: // ------------------- 178  P 61 2 (0 0 -1) #178
      return "P 61 2 (0 0 -1)";
    case 179: // ------------------- 179  P 65 2 (0 0 1) #179
      return "P 65 2 (0 0 1)";
    case 180: // ------------------- 180  P 62 2c (0 0 1) #180
      return "P 62 2c (0 0 1)";
    case 181: // ------------------- 181  P 64 2c (0 0 -1) #181
      return "P 64 2c (0 0 -1)";
    case 182: // ------------------- 182  P 6c 2c #182
      return "P 6c 2c";
    case 183: // ------------------- 183  P 6 -2 #183
      return "P 6 -2";
    case 184: // ------------------- 184  P 6 -2c #184
      return "P 6 -2c";
    case 185: // ------------------- 185  P 6c -2 #185
      return "P 6c -2";
    case 186: // ------------------- 186  P 6c -2c #186
      return "P 6c -2c";
    case 187: // ------------------- 187  P -6 2 #187
      return "P -6 2";
    case 188: // ------------------- 188  P -6c 2 #188
      return "P -6c 2";
    case 189: // ------------------- 189  P -6 -2 #189
      return "P -6 -2";
    case 190: // ------------------- 190  P -6c -2c #190
      return "P -6c -2c";
    case 191: // ------------------- 191  -P 6 2 #191
      return "-P 6 2";
    case 192: // ------------------- 192  -P 6 2c #192
      return "-P 6 2c";
    case 193: // ------------------- 193  -P 6c 2 #193
      return "-P 6c 2";
    case 194: // ------------------- 194  -P 6c 2c #194
      return "-P 6c 2c";
    case 195: // ------------------- 195  P 2 2 3 #195
      return "P 2 2 3";
    case 196: // ------------------- 196  F 2 2 3 #196
      return "F 2 2 3";
    case 197: // ------------------- 197  I 2 2 3 #197
      return "I 2 2 3";
    case 198: // ------------------- 198  P 2ac 2ab 3 #198
      return "P 2ac 2ab 3";
    case 199: // ------------------- 199  I 2b 2c 3 #199
      return "I 2b 2c 3";
    case 200: // ------------------- 200  -P 2 2 3 #200
      return "-P 2 2 3";
    case 201:
      if (setting==SG_SETTING_1) return "P 2 2 3 -1n"; // ------------------- 201  P 2 2 3 -1n #201 (setting 1)
      else return "-P 2ab 2bc 3";                      // ------------------- 201  -P 2ab 2bc 3 #201 (setting 2)
    case 202: // ------------------- 202  -F 2 2 3 #202
      return "-F 2 2 3";
    case 203:
      if (setting==SG_SETTING_1) return "F 2 2 3 -1d"; // ------------------- 203  F 2 2 3 -1d #203 (setting 1)
      else return "-F 2uv 2vw 3";                      // ------------------- 203  -F 2uv 2vw 3 #203 (setting 2)
    case 204: // ------------------- 204  -I 2 2 3 #204
      return "-I 2 2 3";
    case 205: // ------------------- 205  -P 2ac 2ab 3 #205
      return "-P 2ac 2ab 3";
    case 206: // ------------------- 206  -I 2b 2c 3 #206
      return "-I 2b 2c 3";
    case 207: // ------------------- 207  P 4 2 3 #207
      return "P 4 2 3";
    case 208: // ------------------- 208  P 4n 2 3 #208
      return "P 4n 2 3";
    case 209: // ------------------- 209  F 4 2 3 #209
      return "F 4 2 3";
    case 210: // ------------------- 210  F 4d 2 3 #210
      return "F 4d 2 3";
    case 211: // ------------------- 211  I 4 2 3 #211
      return "I 4 2 3";
    case 212: // ------------------- 212  P 4acd 2ab 3 #212
      return "P 4acd 2ab 3";
    case 213: // ------------------- 213  P 4bd 2ab 3 #213
      return "P 4bd 2ab 3";
    case 214: // ------------------- 214  I 4bd 2c 3 #214
      return "I 4bd 2c 3";
    case 215: // ------------------- 215  P -4 2 3 #215
      return "P -4 2 3";
    case 216: // ------------------- 216  F -4 2 3 #216
      return "F -4 2 3";
    case 217: // ------------------- 217  I -4 2 3 #217
      return "I -4 2 3";
    case 218: // ------------------- 218  P -4n 2 3 #218
      return "P -4n 2 3";
    case 219: // ------------------- 219  F -4c 2 3 #219
      return "F -4c 2 3";
    case 220: // ------------------- 220  I -4bd 2c 3 #220
      return "I -4bd 2c 3";
    case 221: // ------------------- 221  -P 4 2 3 #221
      return "-P 4 2 3";
    case 222:
      if (setting==SG_SETTING_1) return "P 4 2 3 -1n"; // ------------------- 222  P 4 2 3 -1n #222 (setting 1)
      else return "-P 4a 2bc 3";                       // ------------------- 222  -P 4a 2bc 3 #222 (setting 2)
    case 223: // ------------------- 223  -P 4n 2 3 #223
      return "-P 4n 2 3";
    case 224:
      if (setting==SG_SETTING_1) return "P 4n 2 3 -1n"; // ------------------- 224  P 4n 2 3 -1n #224 (setting 1)
      else return "-P 4bc 2bc 3";                       // ------------------- 224  -P 4bc 2bc 3 #224 (setting 2)
    case 225: // ------------------- 225  -F 4 2 3 #225
      return "-F 4 2 3";
    case 226: // ------------------- 226  -F 4c 2 3 #226
      return "-F 4c 2 3";
    case 227:
      if (setting==SG_SETTING_1) return "F 4d 2 3 -1d"; // ------------------- 227  -F 4vw 2vw 3 #227 (setting 1)
      else return "-F 4vw 2vw 3";                       // ------------------- 227  -F 4vw 2vw 3 #227 (setting 2)
    case 228:
      if (setting==SG_SETTING_1) return "F 4d 2 3 -1cd"; // ------------------- 228  F 4d 2 3 -1cd #228 (setting 1)
      else return "-F 4cvw 2vw 3";                       // ------------------- 228  -F 4cvw 2vw 3 #228 (setting 2)
    case 229: // ------------------- 229  -I 4 2 3 #229
      return "-I 4 2 3";
    case 230: // ------------------- 230  -I 4bd 2c 3 #230
      return "-I 4bd 2c 3";
  }
  // No space group found
  return "";
}

// ***************************************************************************
// GetLaueLabel
// ***************************************************************************
string GetLaueLabel(string& point_group) {
  string laue = "";
  //DX+ME20190708 - changed subsequent "if" to "else if" -> efficiency
  // -1 
  if (point_group=="1" || point_group=="-1"){
    laue = "-1";
  }
  // 2/m
  else if(point_group=="2" || point_group=="m" || point_group=="2/m"){
    laue = "2/m";
  }
  // mmm
  else if(point_group=="222" || point_group=="mm2" || point_group=="mmm"){
    laue = "mmm";
  }
  // 4/m
  else if(point_group=="4" || point_group=="-4" || point_group=="4/m"){
    laue = "4/m";
  }
  // 4/mmmm
  else if(point_group=="422" || point_group=="4mm" || point_group=="-42m" || point_group=="-4m2" || point_group=="4/mmm"){
    laue = "4/mmm";
  }
  // -3
  else if(point_group=="3" || point_group=="-3"){
    laue = "-3";
  }
  // -3m
  else if(point_group=="312" || point_group=="321" || point_group=="32" || point_group=="31m" || point_group=="3m1" || point_group=="-31m" || point_group=="-3m1" || point_group=="-3m" || point_group=="3m"){
    laue = "-3m";
  }
  // 6/m
  else if(point_group=="6" || point_group=="-6" || point_group=="6/m"){
    laue = "6/m";
  }
  // 6/mmm
  else if(point_group=="622" || point_group=="6mm" || point_group=="-6m2" || point_group=="-62m" || point_group=="6/mmm"){
    laue = "6/mmm";
  }
  // m-3
  else if(point_group=="23" || point_group=="m-3"){
    laue = "m-3";
  }
  // m-3m
  else if(point_group=="432" || point_group=="-43m" || point_group=="m-3m"){
    laue = "m-3m";
  }
  return laue;
}

// ***************************************************************************
// GetSpaceGroupLabel
// ***************************************************************************
string GetSpaceGroupLabel(int spacegroupnumber) {
  string spacegrouplabel;
  spacegrouplabel="#"+aurostd::utype2string(spacegroupnumber);
  return spacegrouplabel;
}

// **************************************************************************
// Function MetricTensor
// **************************************************************************
// this function returns the metric tensor
// Corey Oses
xmatrix<double> MetricTensor(const xstructure& a) {return MetricTensor(a.lattice,a.scale);}

xmatrix<double> MetricTensor(const xmatrix<double>& lattice,double scale) {
  if(lattice.rows!=lattice.cols){
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Dimension mismatch, should be square lattice matrix.",_VALUE_ILLEGAL_);
  }
  xmatrix<double> metric_tensor(lattice.rows,lattice.cols);
  for(int i=lattice.lrows;i<=lattice.urows;i++){ //CO20190520
    for(int j=lattice.lcols;j<=lattice.ucols;j++){ //CO20190520
      metric_tensor(i,j)=aurostd::scalar_product(lattice(i),lattice(j));
    }
  }
  return scale*metric_tensor;
}

// **************************************************************************
// Function ReciprocalLattice
// **************************************************************************
// this function returns the reciprocal lattice
// from the original one... it takes the scale !
// Stefano Curtarolo
xmatrix<double> ReciprocalLattice(const xstructure& a){return ReciprocalLattice(a.lattice,a.scale);}

xmatrix<double> ReciprocalLattice(const xmatrix<double>& rlattice,double scale) {        // AFLOW_FUNCTION_IMPLEMENTATION
  xvector<double> a1(3),a2(3),a3(3),b1(3),b2(3),b3(3);
  xmatrix<double> klattice(3,3);  // kvectors are RAWS
  double norm;
  a1[1]=rlattice[1][1];a1[2]=rlattice[1][2];a1[3]=rlattice[1][3];
  a2[1]=rlattice[2][1];a2[2]=rlattice[2][2];a2[3]=rlattice[2][3];
  a3[1]=rlattice[3][1];a3[2]=rlattice[3][2];a3[3]=rlattice[3][3];
  norm=2.0*pi/(det(rlattice)*scale);
  b1=norm*vector_product(a2,a3);
  b2=norm*vector_product(a3,a1);
  b3=norm*vector_product(a1,a2);
  klattice[1][1]=b1[1];klattice[1][2]=b1[2];klattice[1][3]=b1[3];
  klattice[2][1]=b2[1];klattice[2][2]=b2[2];klattice[2][3]=b2[3];
  klattice[3][1]=b3[1];klattice[3][2]=b3[2];klattice[3][3]=b3[3];
  return klattice;
}

//xmatrix<double> ReciprocalLattice(const xmatrix<double>& rlattice) {        // AFLOW_FUNCTION_IMPLEMENTATION
//  return ReciprocalLattice(rlattice,1.0);
//}

// **************************************************************************
// Function KPPRA
// **************************************************************************
// This function calculates k1,k2,k3 starting from real lattice and NK total
// the function does not normalize with number of atoms so the calculation
// must be done somewhere else
string KPPRA(int& k1,int& k2,int& k3,const xmatrix<double>& rlattice,const int& NK) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream aus("");
  aus.precision(5);
  xmatrix<double> klattice(3,3);
  klattice=ReciprocalLattice(rlattice,1.0);
  xvector<double> b1(3),b2(3),b3(3),db1(3),db2(3),db3(3);
  double nb1,nb2,nb3;
  k1=1;k2=1;k3=1;
  int kk1,kk2,kk3,kk;
  b1=klattice(1);nb1=aurostd::modulus(b1);
  b2=klattice(2);nb2=aurostd::modulus(b2);
  b3=klattice(3);nb3=aurostd::modulus(b3);
  if(LDEBUG) aus << "KPPRA LDEBUG:  " << endl << rlattice << endl << endl << klattice << endl << b1 << endl << b2 << endl << b3 << endl << nb1 << endl << nb2 << endl << nb3 << endl;
  if(LDEBUG) aus << "KPPRA LDEBUG:  " << nb1 << " " << nb2 << " " << nb3 << " " << endl;
  if(NK>1) {
    bool found=FALSE;
    double dkdelta,dk;
    dkdelta=0.999;
    dk=aurostd::min(nb1,nb2,nb3);
    kk1=0;kk2=0;kk3=0;
    kk=0;
    int iverbose=0;
    while(!found) {
      kk++;
      if(dk<=1e-5) {
        k1=1;k2=1;k3=1;
        aus << "00000  MESSAGE KPOINTS KPPRA minimum not found k=[" << k1 << "," << k2 << "," << k3 << "]=" << k1*k2*k3 << endl;
      }
      kk1=(int) floor((double) nb1/dk);db1=b1/((double) kk1);
      kk2=(int) floor((double) nb2/dk);db2=b2/((double) kk2);
      kk3=(int) floor((double) nb3/dk);db3=b3/((double) kk3);
      if(kk1+kk2+kk3>iverbose) {
        //  if(!mod(kk,50) || kk1*kk2*kk3>=NK)
        aus << "00000  MESSAGE KPOINTS KPPRA minimizing k=[" << kk1 << "," << kk2 << "," << kk3 << "]=" << kk1*kk2*kk3 << " =[" << aurostd::modulus(db1) << "," << aurostd::modulus(db2) << "," << aurostd::modulus(db3) << "]   dk=" << dk << endl;
        iverbose=kk1+kk2+kk3;
      }
      if(kk1*kk2*kk3>=NK) {
        k1=kk1;k2=kk2;k3=kk3;
        found=TRUE;
      }
      dk=dk*dkdelta;
    }
  } else { // force 1 1 1 for Gamma
    k1=1;k2=1;k3=1;
  }
  db1=b1/((double) k1); db2=b2/((double) k2); db3=b3/((double) k3);
  aus << "00000  MESSAGE KPOINTS KPPRA routine [" << k1 << "," << k2 << "," << k3 << "]=" << k1*k2*k3 << "=[" << aurostd::modulus(db1) << "," << aurostd::modulus(db2) << "," << aurostd::modulus(db3) << "]   "  << endl;
  return aus.str();
}

string KPPRA_LAT(int& k1,int& k2,int& k3,const xmatrix<double>& rlattice,const int& NK) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream aus("");
  aus.precision(5);
  xvector<double> kdata(6),kdatagrid(6);
  xmatrix<double> klattice(3,3);
  klattice=ReciprocalLattice(rlattice,1.0);
  klattice=GetStandardPrimitive(klattice);
  kdata=Getabc_angles(klattice,DEGREES);

  xvector<double> b1(3),b2(3),b3(3),db1(3),db2(3),db3(3);
  double nb1,nb2,nb3;
  k1=1;k2=1;k3=1;
  int kk1,kk2,kk3,kk;
  b1=klattice(1);nb1=aurostd::modulus(b1);
  b2=klattice(2);nb2=aurostd::modulus(b2);
  b3=klattice(3);nb3=aurostd::modulus(b3);
  if(LDEBUG) aus << "KPPRA LDEBUG:  " << endl << rlattice << endl << endl << klattice << endl << b1 << endl << b2 << endl << b3 << endl << nb1 << endl << nb2 << endl << nb3 << endl;
  if(LDEBUG) aus << "KPPRA LDEBUG:  " << nb1 << " " << nb2 << " " << nb3 << " " << endl;
  if(NK>1) {
    bool found=FALSE;
    double dkdelta,dk;
    dkdelta=0.999;
    dk=aurostd::min(nb1,nb2,nb3);
    kk1=0;kk2=0;kk3=0;
    kk=0;
    int iverbose=0;
    while(!found) {
      kk++;
      if(dk<=1e-5) {
        k1=1;k2=1;k3=1;
        aus << "00000  MESSAGE KPOINTS KPPRA minimum not found k=[" << k1 << "," << k2 << "," << k3 << "]=" << k1*k2*k3 << endl;
      }
      kk1=(int) floor((double) nb1/dk);db1=b1/((double) kk1);
      kk2=(int) floor((double) nb2/dk);db2=b2/((double) kk2);
      kk3=(int) floor((double) nb3/dk);db3=b3/((double) kk3);
      kdatagrid=kdata;kdatagrid[1]=1.0;kdatagrid[2]=1.0*aurostd::modulus(db2)/aurostd::modulus(db1);kdatagrid[3]=1.0*aurostd::modulus(db3)/aurostd::modulus(db1);

      if(kk1+kk2+kk3>iverbose) {
        //  if(!mod(kk,50) || kk1*kk2*kk3>=NK)
        aus << "00000  MESSAGE KPOINTS KPPRA minimizing k=[" << kk1 << "," << kk2 << "," << kk3 << "," << GetLatticeType(kdatagrid) << "]=" << kk1*kk2*kk3 << " = [" << aurostd::modulus(db1) << "," << aurostd::modulus(db2) << "," << aurostd::modulus(db3) << "]   dk=" << dk << " " << endl;
        iverbose=kk1+kk2+kk3;
      }
      if(kk1*kk2*kk3>=NK) {
        k1=kk1;k2=kk2;k3=kk3;
        found=TRUE;
      }
      dk=dk*dkdelta;
    }
  } else { // force 1 1 1 for Gamma
    k1=1;k2=1;k3=1;
  }
  db1=b1/((double) k1); db2=b2/((double) k2); db3=b3/((double) k3);
  aus << "00000  MESSAGE KPOINTS KPPRA routine [" << k1 << "," << k2 << "," << k3 << "]=" << k1*k2*k3 << "=[" << aurostd::modulus(db1) << "," << aurostd::modulus(db2) << "," << aurostd::modulus(db3) << "]   "  << endl;
  return aus.str();
}


string KPPRA(xstructure& str,const int& _NK) {
  //  cerr << "KPPRA" << endl;
  int NK=1;
  NK= (int) ((double) _NK/str.atoms.size()+0.5);if(NK<1) NK=1;
  int k1=1,k2=1,k3=1;
  xmatrix<double> rlattice=str.lattice;
  rlattice=str.scale*rlattice;
  // string stringKPPRA=KPPRA_LAT(k1,k2,k3,rlattice,NK);
  string stringKPPRA=KPPRA(k1,k2,k3,rlattice,NK);
  str.kpoints_k1=k1;str.kpoints_k2=k2;str.kpoints_k3=k3;
  str.kpoints_kmax=max(str.kpoints_k1,str.kpoints_k2,str.kpoints_k3);
  str.kpoints_kppra=str.kpoints_k1*str.kpoints_k2*str.kpoints_k3*str.atoms.size();
  return stringKPPRA;
}

// **************************************************************************
// Function KPPRA_DELTA
// **************************************************************************
// This function calculates k1,k2,k3 starting from real lattice and DK
string KPPRA_DELTA(int& k1,int& k2,int& k3,const xmatrix<double>& rlattice,const double& DK) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream aus("");
  aus.precision(5);
  xmatrix<double> klattice(3,3);
  klattice=ReciprocalLattice(rlattice,1.0);
  xvector<double> b1(3),b2(3),b3(3),db1(3),db2(3),db3(3);
  double nb1,nb2,nb3;
  k1=1;k2=1;k3=1;
  int kk1,kk2,kk3,kk;
  b1=klattice(1);nb1=aurostd::modulus(b1);
  b2=klattice(2);nb2=aurostd::modulus(b2);
  b3=klattice(3);nb3=aurostd::modulus(b3);
  if(LDEBUG) aus << "KPPRA LDEBUG:  " << endl << rlattice << endl << endl << klattice << endl << b1 << endl << b2 << endl << b3 << endl << nb1 << endl << nb2 << endl << nb3 << endl;
  if(LDEBUG) aus << "KPPRA LDEBUG:  " << nb1 << " " << nb2 << " " << nb3 << " " << endl;
  if(DK>1.0e-6) {
    bool found=FALSE;
    double dkdelta,dk;
    dkdelta=0.9999;
    dk=aurostd::min(nb1,nb2,nb3);
    kk1=0;kk2=0;kk3=0;
    kk=0;
    int iverbose=0;
    while(!found) {
      kk++;
      kk1=(int) floor((double) nb1/dk);db1=b1/((double) kk1);
      kk2=(int) floor((double) nb2/dk);db2=b2/((double) kk2);
      kk3=(int) floor((double) nb3/dk);db3=b3/((double) kk3);
      if(kk1+kk2+kk3>iverbose) {
        //  if(!mod(kk,50) || kk1*kk2*kk3>=DK)
        aus << "00000  MESSAGE KPOINTS KPPRA minimizing k=[" << kk1 << "," << kk2 << "," << kk3 << "]=" << kk1*kk2*kk3 << " =[" << aurostd::modulus(db1) << "," << aurostd::modulus(db2) << "," << aurostd::modulus(db3) << "]   dk=" << dk << endl;
        iverbose=kk1+kk2+kk3;
      }
      if((aurostd::modulus(db1)<DK) && (aurostd::modulus(db2)<DK) && (aurostd::modulus(db3)<DK)) {
        found=TRUE;k1=kk1;k2=kk2;k3=kk3;}
      dk=dk*dkdelta;
    }
  } else { // force 1 1 1 for Gamma
    k1=1;k2=1;k3=1;
  }
  db1=b1/((double) k1); db2=b2/((double) k2); db3=b3/((double) k3);
  aus << "00000  MESSAGE KPOINTS KPPRA routine [" << k1 << "," << k2 << "," << k3 << "]=" << k1*k2*k3 << "=[" << aurostd::modulus(db1) << "," << aurostd::modulus(db2) << "," << aurostd::modulus(db3) << "]   "  << endl;
  return aus.str();
}

string KPPRA_DELTA(xstructure& str,const double& DK) {
  //  cerr << "KPPRA_DELTA" << endl;
  int k1=1,k2=1,k3=1;
  xmatrix<double> rlattice=str.lattice;
  rlattice=str.scale*rlattice;
  string stringKPPRA=KPPRA_DELTA(k1,k2,k3,rlattice,DK);
  str.kpoints_k1=k1;str.kpoints_k2=k2;str.kpoints_k3=k3;
  str.kpoints_kmax=max(str.kpoints_k1,str.kpoints_k2,str.kpoints_k3);
  str.kpoints_kppra=str.kpoints_k1*str.kpoints_k2*str.kpoints_k3*str.atoms.size();
  return stringKPPRA;
}

// **************************************************************************
// Function GetNBAND
// **************************************************************************
// returns estimated version of NBANDS starting from
// electrons, ions, spin and ispin
int GetNBANDS(int electrons,int nions,int spineach,bool ispin,int NPAR) {
  double out=0.0;
  out=max(ceil((electrons+4.0)/1.75)+max(nions/1.75,6.0),ceil(0.80*electrons)); // from VASP
  if(ispin) out+=(nions*spineach+1)/2;
  //  out*=1.2;  // safety from vasp
  out*=1.3;      // safety more
  out=out*1.1+5;    // Thu Jun 11 12:08:42 EDT 2009 // METAL PROJECT
  out*=1.075;    // Tue Oct 13 07:59:43 EDT 2009 // ICSD PROJECT
  out+=5;        // Sun Nov  1 10:41:20 EDT 2009 // ICSD PROJECT ORC
  out*=1.03;     // Tue Feb 26 15:15:36 EST 2013 // HELPS dielectric CALS
  out*=1.05;     // Mon Apr 23 13:40:02 EST 2018 // HELPS SCAN
  // cerr << "GetNBANDS=" << out << endl;
  if (nions < 100) {
    out*=std::pow((double) nions,(double) 0.025);  // rescale so for big numbers of ions you get extra bands // Wed Jun 23 12:29:01 EDT 2010
  } else {
    out*=std::pow((double) nions,(double) 0.06);  //ME20191028 - prior scaling factor not sufficient for supercells
  }
  //  cerr << "GetNBANDS=" << out << endl;
  // throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_);
  int nbands=(int)ceil(out);
  //CO20210315 START - adjust for NPAR
  if(NPAR>0){
    int increment=1; //(increase?+1:-1);
    while((nbands%NPAR)!=0){nbands+=increment;}
  }
  //CO20210315 END - adjust for NPAR
  return nbands;
}

// **************************************************************************
// Function GetZVAL from *CAR
// ***************************************************************************
double GetZVAL(const stringstream& sss,vector<double>& vZVAL) {
  xPOTCAR potcar;
  potcar.GetProperties(sss);
  vZVAL.clear();
  for(uint i=0;i<potcar.vZVAL.size();i++) 
    vZVAL.push_back(potcar.vZVAL.at(i));
  return potcar.ZVAL_sum;
}

double GetZVAL(const _xvasp& xvasp,vector<double>& vZVAL) {
  return GetZVAL(xvasp.POTCAR,vZVAL);
}

double GetZVAL(const string& directory,vector<double>& vZVAL) {
  stringstream sss("");vector<string> vfile;
  string file2search="POTCAR,OUTCAR,POTCAR.relax1,POTCAR.relax2,POTCAR.static,POTCAR.bands,OUTCAR.relax1,OUTCAR.relax2,OUTCAR.static,OUTCAR.bands";
  aurostd::string2tokens(file2search,vfile,",");
  for(uint i=0;i<vfile.size()&&!sss.str().length();i++) {
    if(!sss.str().length() && aurostd::FileExist(directory+"/"+vfile.at(i))) aurostd::file2stringstream(directory+"/"+vfile.at(i),sss);  
    if(!sss.str().length() && aurostd::EFileExist(directory+"/"+vfile.at(i))) aurostd::efile2stringstream(directory+"/"+vfile.at(i),sss);
  } //  cerr << sss.str() << endl;
  return GetZVAL(sss,vZVAL);
}

double GetCellAtomZVAL(const stringstream& sss,vector<double>& vZVAL,const stringstream& sstr,vector<double>& sZVAL,string mode) {
  vZVAL.clear();sZVAL.clear();
  GetZVAL(sss,vZVAL);
  stringstream aus(sstr.str()); aus << sstr.str();
  xstructure xstr(aus,IOAFLOW_AUTO);
  if(mode=="CELL" || mode=="") {
    double CellZVAL=xstr.GetZVAL(vZVAL);
    for(uint i=0;i<vZVAL.size();i++) 
      sZVAL.push_back(double(vZVAL.at(i)*double(xstr.num_each_type.at(i))));
    return CellZVAL;
  }
  if(mode=="ATOM") {
    double CellZVAL=xstr.GetZVAL(vZVAL)/double(xstr.atoms.size());
    for(uint i=0;i<vZVAL.size();i++) 
      sZVAL.push_back(double(vZVAL.at(i)*double(xstr.num_each_type.at(i)))/double(xstr.atoms.size()));
    return CellZVAL;
  }
  return 0.0;
}  

double GetCellAtomZVAL(const string& directory,vector<double>& vZVAL,vector<double>& sZVAL,string mode) {  // from directory POT/POS returns total ZVAL cell, vZVAL and sZVAL
  vector<string> vfile;
  // search for data
  stringstream sss("");
  aurostd::string2tokens("POTCAR,OUTCAR,POTCAR.relax1,POTCAR.relax2,POTCAR.static,POTCAR.bands,OUTCAR.relax1,OUTCAR.relax2,OUTCAR.static,OUTCAR.bands",vfile,",");
  for(uint i=0;i<vfile.size()&&!sss.str().length();i++) {
    if(!sss.str().length() && aurostd::FileExist(directory+"/"+vfile.at(i))) aurostd::file2stringstream(directory+"/"+vfile.at(i),sss);  
    if(!sss.str().length() && aurostd::EFileExist(directory+"/"+vfile.at(i))) aurostd::efile2stringstream(directory+"/"+vfile.at(i),sss);
  } //   cerr << sss.str() << endl;

  // search for xstructure
  stringstream sstr("");
  aurostd::string2tokens("POSCAR,CONTCAR,POSCAR.relax1,POSCAR.relax2,POSCAR.static,POSCAR.bands,CONTCAR.relax1,CONTCAR.relax2,CONTCAR.static,CONTCAR.bands",vfile,",");
  for(uint i=0;i<vfile.size()&&!sstr.str().length();i++) {
    //    cerr << "file=" << directory+"/"+vfile.at(i) << endl;
    if(!sstr.str().length() && aurostd::FileExist(directory+"/"+vfile.at(i))) aurostd::file2stringstream(directory+"/"+vfile.at(i),sstr);  
    if(!sstr.str().length() && aurostd::EFileExist(directory+"/"+vfile.at(i))) aurostd::efile2stringstream(directory+"/"+vfile.at(i),sstr);
  } //  cerr << sstr.str() << endl;
  // done
  return GetCellAtomZVAL(sss,vZVAL,sstr,sZVAL,mode);
}

// ***************************************************************************
// Function GetZVAL
// ***************************************************************************
// Given the ZVAL of each species, it returns total ZVAL of cell
double xstructure::GetZVAL(const vector<double>& vZVAL) {
  stringstream message;
  if(num_each_type.size()!=vZVAL.size()) {
    message << "num_each_type.size()=" << num_each_type.size() << endl;
    message << "vZVAL.size()=" << vZVAL.size() << endl;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
  }
  double CellZVAL=0.0;
  for(uint i=0;i<vZVAL.size();i++)  
    CellZVAL+=double(vZVAL.at(i)*num_each_type.at(i));
  return CellZVAL;
}

// **************************************************************************
// Function GetPOMASS from *CAR
// ***************************************************************************
double GetPOMASS(const stringstream& sss,vector<double>& vPOMASS) {
  xPOTCAR potcar;
  potcar.GetProperties(sss);
  vPOMASS.clear();
  for(uint i=0;i<potcar.vPOMASS.size();i++) 
    vPOMASS.push_back(potcar.vPOMASS.at(i));
  return potcar.POMASS_sum;
}

double GetPOMASS(const _xvasp& xvasp,vector<double>& vPOMASS) {
  return GetPOMASS(xvasp.POTCAR,vPOMASS);
}

double GetPOMASS(const string& directory,vector<double>& vPOMASS) {
  stringstream sss("");vector<string> vfile,_aus;
  aurostd::string2tokens("POTCAR,OUTCAR,POTCAR.relax1,POTCAR.relax2,POTCAR.static,POTCAR.bands,OUTCAR.relax1,OUTCAR.relax2,OUTCAR.static,OUTCAR.bands",vfile,",");
  for(uint i=0;i<vfile.size()&&!sss.str().length();i++) {
    if(!sss.str().length() && aurostd::FileExist(directory+"/"+vfile.at(i))) aurostd::file2stringstream(directory+"/"+vfile.at(i),sss);  
    if(!sss.str().length() && aurostd::EFileExist(directory+"/"+vfile.at(i))) aurostd::efile2stringstream(directory+"/"+vfile.at(i),sss);
  } //  cerr << sss.str() << endl;
  return GetPOMASS(sss,vPOMASS);
}

double GetCellAtomPOMASS(const stringstream& sss,vector<double>& vPOMASS,const stringstream& sstr,vector<double>& sPOMASS,string mode) {
  vPOMASS.clear();sPOMASS.clear();
  GetPOMASS(sss,vPOMASS);
  stringstream aus(sstr.str()); aus << sstr.str();
  xstructure xstr(aus,IOAFLOW_AUTO);
  if(mode=="CELL" || mode=="") {
    double CellPOMASS=xstr.GetPOMASS(vPOMASS);
    for(uint i=0;i<vPOMASS.size();i++) 
      sPOMASS.push_back(double(vPOMASS.at(i)*double(xstr.num_each_type.at(i))));
    return CellPOMASS;
  }
  if(mode=="ATOM") {
    double CellPOMASS=xstr.GetPOMASS(vPOMASS)/double(xstr.atoms.size());
    for(uint i=0;i<vPOMASS.size();i++) 
      sPOMASS.push_back(double(vPOMASS.at(i)*double(xstr.num_each_type.at(i)))/double(xstr.atoms.size()));
    return CellPOMASS;
  }
  return 0.0;
}  

double GetCellAtomPOMASS(const string& directory,vector<double>& vPOMASS,vector<double>& sPOMASS,string mode) {  // from directory POT/POS returns total POMASS cell, vPOMASS and sPOMASS
  vector<string> vfile,_aus;
  // search for data
  stringstream sss("");
  aurostd::string2tokens("POTCAR,OUTCAR,POTCAR.relax1,POTCAR.relax2,POTCAR.static,POTCAR.bands,OUTCAR.relax1,OUTCAR.relax2,OUTCAR.static,OUTCAR.bands",vfile,",");

  for(uint i=0;i<vfile.size()&&!sss.str().length();i++) {
    if(!sss.str().length() && aurostd::FileExist(directory+"/"+vfile.at(i))) aurostd::file2stringstream(directory+"/"+vfile.at(i),sss);  
    if(!sss.str().length() && aurostd::EFileExist(directory+"/"+vfile.at(i))) aurostd::efile2stringstream(directory+"/"+vfile.at(i),sss);
  } //  cerr << sss.str() << endl;

  // search for xstructure
  stringstream sstr("");
  aurostd::string2tokens("POSCAR,CONTCAR,POSCAR.relax1,POSCAR.relax2,POSCAR.static,POSCAR.bands,CONTCAR.relax1,CONTCAR.relax2,CONTCAR.static,CONTCAR.bands",vfile,",");

  for(uint i=0;i<vfile.size()&&!sstr.str().length();i++) {
    if(!sstr.str().length() && aurostd::FileExist(directory+"/"+vfile.at(i))) aurostd::file2stringstream(directory+"/"+vfile.at(i),sstr);  
    if(!sstr.str().length() && aurostd::EFileExist(directory+"/"+vfile.at(i))) aurostd::efile2stringstream(directory+"/"+vfile.at(i),sstr);
  } //  cerr << sstr.str() << endl;

  // done
  return GetCellAtomPOMASS(sss,vPOMASS,sstr,sPOMASS,mode);
}

// ***************************************************************************
// Function GetPOMASS
// ***************************************************************************
// Given the POMASS of each species, it returns total POMASS of cell
double xstructure::GetPOMASS(const vector<double>& vPOMASS) {
  stringstream message;
  if(num_each_type.size()!=vPOMASS.size()) {
    message << "num_each_type.size()=" << num_each_type.size() << endl;
    message << "vPOMASS.size()=" << vPOMASS.size() << endl;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
  }
  double CellPOMASS=0.0;
  for(uint i=0;i<vPOMASS.size();i++)  
    CellPOMASS+=double(vPOMASS.at(i)*num_each_type.at(i));
  return CellPOMASS;
}


// **************************************************************************
// Function GetVols, det,
// **************************************************************************
// Wrap up for GetVol functions ... useful to bring
// convasp framework to aflow.

double GetVol(const xmatrix<double>& lat) {
  return abs(det(lat));
}        // AFLOW_FUNCTION_IMPLEMENTATION
double det(const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3) {
  return scalar_product(v1,vector_product(v2,v3));
}
double det(const double& a11,const double& a12,const double& a13,
    const double& a21,const double& a22,const double& a23,
    const double& a31,const double& a32,const double& a33) {
  return a11*a22*a33+a12*a23*a31+a13*a21*a32-a13*a22*a31-a12*a21*a33-a11*a23*a32;
}
double GetVol(const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3) {
  return abs(det(v1,v2,v3));
}        // AFLOW_FUNCTION_IMPLEMENTATION


// ***************************************************************************
// Function Getabc_angles
// ***************************************************************************
// This function returns a,b,c,alpha,beta,gamma for the
// cell given the lattice xvectors. Dane Morgan, adjusted by SC

//#define _Getabc_angles Getabc_angles
//#define _Getabc_angles __NO_USE_Sortabc_angles

xvector<double> Getabc_angles(const xmatrix<double>& lat,int mode) {  // AFLOW_FUNCTION_IMPLEMENTATION
  xvector<double> data(6);
  data(1)=aurostd::modulus(lat(1));
  data(2)=aurostd::modulus(lat(2));
  data(3)=aurostd::modulus(lat(3));
  data(4)=angle(lat(2),lat(3));
  data(5)=angle(lat(3),lat(1));
  data(6)=angle(lat(1),lat(2));
  if(mode==DEGREES) {
    data(4)*=rad2deg;
    data(5)*=rad2deg;
    data(6)*=rad2deg;
  }
  return data;
}

xvector<double> Getabc_angles(const xmatrix<double>& lat,const xvector<int>& permut,int mode) {        // AFLOW_FUNCTION_IMPLEMENTATION
  xvector<double> data(6);
  data(1)=aurostd::modulus(lat(1));
  data(2)=aurostd::modulus(lat(2));
  data(3)=aurostd::modulus(lat(3));
  data(4)=angle(lat(2),lat(3));
  data(5)=angle(lat(3),lat(1));
  data(6)=angle(lat(1),lat(2));
  if(mode==DEGREES) {
    data(4)*=rad2deg;
    data(5)*=rad2deg;
    data(6)*=rad2deg;
  }
  // with permutation - from AVDV
  int i,imin,imax,imid;
  double dmin,dmax;
  dmax=dmin=data(1);
  imax=imin=imid=1;
  for(i=1;i<=3;i++) {
    if(dmax<=data(i)) {dmax=data(i);imax=i; }
    if(dmin>data(i)) {dmin=data(i);imin=i; }}
  // set imid
  for(i=1;i<=3;i++)if(i!=imin && i!=imax) imid=i;
  // if all lattice parameters are equal length, numerical noise may cause imin=imax
  if(imin==imax) for(i=1;i<=3;i++)if(i!=imin && i!=imid) imax=i;
  permut[1]=imin;permut[2]=imid;permut[3]=imax;
  // done permutation - from AVDV
  return data;
}

xvector<double> Getabc_angles(const xvector<double>& r1,      // AFLOW_FUNCTION_IMPLEMENTATION
    const xvector<double>& r2,      // AFLOW_FUNCTION_IMPLEMENTATION
    const xvector<double>& r3,      // AFLOW_FUNCTION_IMPLEMENTATION
    int mode) {              // AFLOW_FUNCTION_IMPLEMENTATION
  xmatrix<double> lat(3,3);
  lat(1,1)=r1(1);lat(1,2)=r1(2);lat(1,3)=r1(3);
  lat(2,1)=r2(1);lat(2,2)=r2(2);lat(2,3)=r2(3);
  lat(3,1)=r3(1);lat(3,2)=r3(2);lat(3,3)=r3(3);
  return Getabc_angles(lat,mode);
}

xvector<double> Getabc_angles(const xvector<double>& r1,      // AFLOW_FUNCTION_IMPLEMENTATION
    const xvector<double>& r2,      // AFLOW_FUNCTION_IMPLEMENTATION
    const xvector<double>& r3,      // AFLOW_FUNCTION_IMPLEMENTATION
    const xvector<int>& permut,     // AFLOW_FUNCTION_IMPLEMENTATION
    int mode) {              // AFLOW_FUNCTION_IMPLEMENTATION
  xmatrix<double> lat(3,3);
  lat(1,1)=r1(1);lat(1,2)=r1(2);lat(1,3)=r1(3);
  lat(2,1)=r2(1);lat(2,2)=r2(2);lat(2,3)=r2(3);
  lat(3,1)=r3(1);lat(3,2)=r3(2);lat(3,3)=r3(3);
  return Getabc_angles(lat,permut,mode);
}

xvector<double> Sortabc_angles(const xmatrix<double>& lat,int mode) {        // AFLOW_FUNCTION_IMPLEMENTATION
  // with permutation - from AVDV
  int i,imin,imax,imid;
  double dmin,dmax;
  dmax=dmin=aurostd::modulus(lat(1));
  imax=imin=imid=1;
  for(i=1;i<=3;i++) {
    if(dmax<=aurostd::modulus(lat(i))) {dmax=aurostd::modulus(lat(i));imax=i; }
    if(dmin>aurostd::modulus(lat(i))) {dmin=aurostd::modulus(lat(i));imin=i; }}
  // set imid
  for(i=1;i<=3;i++) if(i!=imin && i!=imax) imid=i;
  // if all lattice parameters are equal length, numerical noise may cause imin=imax
  if(imin==imax) for(i=1;i<=3;i++) if(i!=imin && i!=imid) imax=i;

  xvector<double> data(6);
  data(1)=aurostd::modulus(lat(imin));
  data(2)=aurostd::modulus(lat(imid));
  data(3)=aurostd::modulus(lat(imax));
  data(4)=angle(lat(imid),lat(imax));
  data(5)=angle(lat(imin),lat(imax));
  data(6)=angle(lat(imin),lat(imid));
  if(mode==DEGREES) {
    data(4)*=rad2deg;
    data(5)*=rad2deg;
    data(6)*=rad2deg;
  }
  return data;
}

// **************************************************************************
// Function GetClat
// **************************************************************************
// This function returns cartesian lattice xvectors for the cell given
// a,b,c,alpha,beta,gamma. Assumes angles are in radians.
// This routine aligns a along +X, b in XY plane with
// angle gamma from a in +Y direction (bx=b.a_unit=b*cos[gamma],
// by=sqrt(b^2-bx^2)=b*sin[gamma]), and c is then given by
// cx=c.a_unit=c*cos[beta],
// cy=c.by_unit=c*(cos[alpha]-cos[gamma]*cos[beta])/sin[gamma],
// cz=sqrt(c^2-cx^2-cy^2)
// Dane Morgan, adjusted by SC

xmatrix<double> GetClat(const xvector<double>& abc_angles) {   // AFLOW_FUNCTION_IMPLEMENTATION
  stringstream message;
  xmatrix<double> clattice(3,3);
  double a=abc_angles[1];
  double b=abc_angles[2];
  double c=abc_angles[3];
  double bc= abc_angles[4]*deg2rad; // angle from b to c (remove a)
  double ca= abc_angles[5]*deg2rad; // angle from c to a (remove b)
  double ab= abc_angles[6]*deg2rad; // angle from a to b (remove c)
  //  if(abs(bc)>6.3||abs(ca)>6.3||abs(ab)>6.3) { cerr << _AUROSTD_XLIBS_ERROR_ << "GetClat: angles must be in RADIANT " << endl;throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_);}
  clattice(1,1)=a;
  clattice(2,1)=b*cos(ab);
  clattice(2,2)=b*sin(ab);
  clattice(3,1)=c*cos(ca);
  if(ab<0.00000001) {
    message <<"The angle gamma from a to b is too small" << endl;
    message <<"gamma = " << ab << endl;
    message <<"STOPPING "<< endl;
    message << _AUROSTD_XLIBS_ERROR_ << "ERROR: STOPPING " << endl;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
  }
  clattice(3,2)=c*(cos(bc)-cos(ab)*cos(ca))/sin(ab);
  clattice(3,3)=sqrt(abs(c*c-clattice(3,2)*clattice(3,2)-clattice(3,1)*clattice(3,1)));
  return clattice;
}

xmatrix<double> GetClat(const double &a,const double &b,const double &c,const double &alpha,const double &beta,const double &gamma) {
  stringstream message;
  xmatrix<double> clattice(3,3);
  double bc= alpha*deg2rad; // angle from b to c (remove a)
  double ca= beta*deg2rad; // angle from c to a (remove b)
  double ab= gamma*deg2rad; // angle from a to b (remove c)
  //  if(abs(bc)>6.3||abs(ca)>6.3||abs(ab)>6.3) { cerr << _AUROSTD_XLIBS_ERROR_ << "GetClat: angles must be in RADIANT " << endl;throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_);}
  clattice(1,1)=a;
  clattice(2,1)=b*cos(ab);
  clattice(2,2)=b*sin(ab);
  clattice(3,1)=c*cos(ca);
  if(ab<0.00000001) {
    message <<"The angle gamma from a to b is too small" << endl;
    message <<"gamma = " << ab << endl;
    message <<"STOPPING "<< endl;
    message << _AUROSTD_XLIBS_ERROR_ << "ERROR: STOPPING " << endl;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
  }
  clattice(3,2)=c*(cos(bc)-cos(ab)*cos(ca))/sin(ab);
  clattice(3,3)=sqrt(abs(c*c-clattice(3,2)*clattice(3,2)-clattice(3,1)*clattice(3,1)));
  return clattice;
}

// **************************************************************************
// Function GetIntpolStr
// **************************************************************************
// This function gets the structure that is linearly interpolated
// a fraction f of the way between strA and strB.  Interpolation
// nincludes lattice parameters.  The scale factors of the interpolated
// structures are all set to 1.

xstructure GetIntpolStr(xstructure strA, xstructure strB, const double& f,const string& path_flag) {
  strA=ReScale(strA,1.0);
  strB=ReScale(strB,1.0);
  // Get new lattice params.
  xmatrix<double> lati(3,3);
  xmatrix<double> latA(3,3);latA=strA.lattice;
  xmatrix<double> latB(3,3);latB=strB.lattice;
  xvector<double> dl(3);
  for(int ic=1;ic<=3;ic++) {
    dl=f*(latB(ic)-latA(ic));
    for(int jc=1;jc<=3;jc++) lati(ic,jc)=latA(ic,jc)+dl(jc);
  }
  // Get new cart. coords.
  if(strA.atoms.size()!=strB.atoms.size()) {
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,_AUROSTD_XLIBS_ERROR_+" number of atoms must be the same in both structures!!",_INPUT_ILLEGAL_);
  }
  int size=strA.atoms.size();
  vector<xvector<double> > cposi(size,3);
  for(int iat=0;iat<size;iat++) {
    xvector<double> dp(3);
    dp=strB.atoms.at(iat).cpos-strA.atoms.at(iat).cpos;
    // If path_flag is n/N then the path is taken between
    // nearest images.  Otherwise path is between the atoms given.
    if(path_flag=="n" || path_flag=="N") {
      xvector<double> ddp(3);ddp=C2F(lati,dp);
      for(int ic=1;ic<=3;ic++) {
        ddp(ic)=ddp(ic)-nint(ddp(ic));
      }
      dp=F2C(lati,ddp);
    }
    dp=f*dp;
    cposi.at(iat)=strA.atoms.at(iat).cpos+dp;
  }
  xstructure stri=strA;
  stri.lattice=lati;
  for(int iat=0;iat<size;iat++) {
    stri.atoms.at(iat).cpos=cposi.at(iat);
    stri.atoms.at(iat).fpos=C2F(stri.lattice,stri.atoms.at(iat).cpos);
  }
  return stri;
}

// **************************************************************************
// Function RadiusSphereLattice
// **************************************************************************
// This function returns the radius of the sphere encompassing the whole cell
double RadiusSphereLattice(const xmatrix<double>& lattice,double scale) {
  double radius=0;
  for(int i=-1;i<=1;i++)
    for(int j=-1;j<=1;j++)
      for(int k=-1;k<=1;k++)
        if(aurostd::modulus((double)i*lattice(1)+(double)j*lattice(2)+(double)k*lattice(3))>radius)
          radius=aurostd::modulus((double)i*lattice(1)+(double)j*lattice(2)+(double)k*lattice(3));
  return scale*radius;
}

// **************************************************************************
// Function LatticeDimensionSphere
// **************************************************************************
// This function, given a lattice with vectors lat[3][3], finds the
// dimensions along the unit cell vectors such that a sphere of given radius
// fits within a uniform grid of 2dim[1]x2dim[2]x2dim[3] lattice points
// centered at the origin.
//
// The algorithm works by getting the normal (e.g. n1) to each pair of lattice
// vectors (e.g. a2, a3), scaling this normal to have length radius and
// then projecting this normal parallel to the a2,a3 plane onto the
// remaining lattice vector a1. This will tell us the number of a1 vectors
// needed to make a grid to encompass the sphere.
//
// Written as lat_dimension by Anton Van Der Ven
// Adjusted by Stefano Curtarolo

xvector<int> LatticeDimensionSphere(const xmatrix<double>& _lattice, double radius,double scale) {
  // Adapted from AVDV routine
#if DEBUG_LATTICE_DIMENSIONS
  bool LDEBUG=(FALSE || XHOST.DEBUG);
#endif
  xmatrix<double> lattice; lattice=scale*_lattice;
  int i,j,k;
  xmatrix<double> invlattice(3,3),normals(3,3),frac_normals(3,3);
  xvector<int> dim(3);
  double length;
  normals.clear();
  //get the normals to pairs of lattice vectors of length radius
  if(0) { // without eijk
    for(i=1;i<=3;i++)
      for(j=1;j<=3;j++)
        //     normals(i,j)=lattice((i+1)%3+1,(j+1)%3+1)*lattice((i+2)%3+1,(j+2)%3+1)-  // original
        //	lattice((i+1)%3+1,(j+2)%3+1)*lattice((i+2)%3+1,(j+1)%3+1);            // original
        normals(i,j)=lattice((i+1-1)%3+1,(j+1-1)%3+1)*lattice((i+2-1)%3+1,(j+2-1)%3+1)-
          lattice((i+1-1)%3+1,(j+2-1)%3+1)*lattice((i+2-1)%3+1,(j+1-1)%3+1);
    //  for(i=1;i<=3;i++)
    //   for(j=1;j<=3;j++)
    //    normals(i,j)=vector_product(lattice((i+1-1)%3+1),lattice((i+2-1)%3+1))[j];
  }
  if(1) { // with eijk
    for(int m=1;m<=3;m++)
      for(int n=1;n<=3;n++)
        for(int l=1;l<=3;l++) {
          normals(1,l)+=aurostd::eijk(l,m,n)*lattice(2,m)*lattice(3,n);
          normals(2,l)+=aurostd::eijk(l,m,n)*lattice(3,m)*lattice(1,n);
          normals(3,l)+=aurostd::eijk(l,m,n)*lattice(1,m)*lattice(2,n);
        }
  }
#if DEBUG_LATTICE_DIMENSIONS
  if(LDEBUG) {
    for(uint i=1;i<(uint)normals.rows+1;i++){
      cerr << __AFLOW_FUNC__ << " normals(" << i << ")=" << normals(i) << endl;
    }
  }
#endif
  if(0) { // with eijk and estarijk
    for(int m=1;m<=3;m++)
      for(int n=1;n<=3;n++)
        for(int l=1;l<=3;l++)
          for(int i=1;i<=3;i++)
            for(int j=1;j<=3;j++)
              for(int k=1;k<=3;k++)
                normals(i,l)+=aurostd::eijk(l,m,n)*lattice(j,m)*lattice(k,n)*aurostd::estarijk(i,j,k);
  }
  for(i=1;i<=3;i++) {
    length=aurostd::modulus(normals(i));
    for(j=1;j<=3;j++) normals(i,j)=radius*normals(i,j)/length;
  }
  //get the normals in the coordinates system of the lattice vectors
  invlattice=aurostd::inverse(lattice);
#if DEBUG_LATTICE_DIMENSIONS
  if(LDEBUG) { //CO20190520
    cerr << __AFLOW_FUNC__ << " normals="<< endl;cerr << normals << endl; //CO20190520
    cerr << __AFLOW_FUNC__ << " lattice="<< endl;cerr << lattice << endl; //CO20190520
    cerr << __AFLOW_FUNC__ << " inverse(lattice)="<< endl;cerr << invlattice << endl; //CO20190520
  } //CO20190520
#endif

  for(i=1;i<=3;i++) {
    for(j=1;j<=3;j++) {
      frac_normals(i,j)=0.0;
      for(k=1;k<=3;k++)
        frac_normals(i,j)=frac_normals(i,j)+normals(i,k)*invlattice(k,j);
    }
  }
  //the diagonals of frac_normal contain the dimensions of the lattice grid that
  //encompasses a sphere of radius = radius
  for(i=1;i<=3;i++){
#if DEBUG_LATTICE_DIMENSIONS
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " abs(frac_normals(i,i))=" << abs(frac_normals(i,i)) << endl;}
#endif
    dim(i)=(int)ceil(abs(frac_normals(i,i)));
  }
  if(max(dim)==0) { 
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,_AUROSTD_XLIBS_ERROR_+" dim=0!!",_INPUT_ILLEGAL_);
  }
  return dim;
}

xvector<int> LatticeDimensionSphere(const xstructure& str, double radius) {
  return LatticeDimensionSphere(str.lattice,radius,str.scale); //str.scale*str.lattice,radius);  //CO20171024
}

// **************************************************************************
// F2C and C2F exchange transformations
// **************************************************************************
// Stefano Curtarolo
// xf=(CF)*xc and xc=(FC)*xf   (FC)=inv(CF)
//  A is a property so it is a 2 indices tensor
//  x'=A*x and xf'=Af*xf and xc'=Ac*xc
//  xf'=(CF)*xc'=Af*xf=Af*(CF)*xc   mult inv(CF) *
//  xc'=inv(CF)*Af*(CF)*xc HENCE
//  Ac=inv(CF)*Af*(CF);
//  Af=inv(FC)*Ac*(FC)
//  in another way
//  Ac_ij=inv(CF)_il Af_lm (CF)_mj
//  Af_ij=inv(FC)_il Ac_lm (FC)_mj
//  same for operations, although here FC=inv(CF)
//  see notes on operation of force constant tensors

// **************************************************************************
// -------------------------------------------------------- for COLUMN vectors
xvector<double> F2C(const double& scale,const xmatrix<double>& lattice,const xvector<double>& fpos) {
  xmatrix<double> f2c(3,3);f2c=scale*trasp(lattice);
  return f2c*fpos;                                 // fpos are F components per COLUMS !
}

xvector<double> F2C(const xmatrix<double>& lattice,const xvector<double>& fpos) {
  return F2C(1.0,lattice,fpos); // VASP cartesian coordinates are intended normalized on the scale
}

xvector<double> C2F(const double& scale,const xmatrix<double>& lattice,const xvector<double>& cpos) {
  xmatrix<double> f2c(3,3);f2c=scale*trasp(lattice);
  return inverse(f2c)*cpos;                         // cpos are C components per COLUMS !
}

xvector<double> C2F(const xmatrix<double>& lattice,const xvector<double>& cpos) {
  return C2F(1.0,lattice,cpos); // VASP cartesian coordinates are intended normalized on the scale
}

// **************************************************************************
// ---------------------------------------------- for a set of COLUM matrices
xmatrix<double> F2C(const double& scale,const xmatrix<double>& lattice,const xmatrix<double>& fpos) {
  xmatrix<double> f2c(3,3);f2c=scale*trasp(lattice);
  return f2c*fpos;                                 // fpos are F components per COLUMS !
}

xmatrix<double> F2C(const xmatrix<double>& lattice,const xmatrix<double>& fpos) {
  return F2C(1.0,lattice,fpos); // VASP cartesian coordinates are intended normalized on the scale
}

xmatrix<double> C2F(const double& scale,const xmatrix<double>& lattice,const xmatrix<double>& cpos) {
  xmatrix<double> f2c(3,3);f2c=scale*trasp(lattice);
  return inverse(f2c)*cpos;                         // cpos are C components per COLUMS !
}

xmatrix<double> C2F(const xmatrix<double>& lattice,const xmatrix<double>& cpos) {
  return C2F(1.0,lattice,cpos); // VASP cartesian coordinates are intended normalized on the scale
}

// **************************************************************************
// -------------------------------------------------------- for _atoms with lattice
_atom F2C(const double& scale,const xmatrix<double>& lattice,const _atom& iatom) {
  xmatrix<double> f2c(3,3);f2c=scale*trasp(lattice);
  _atom oatom(iatom); oatom.cpos=f2c*iatom.fpos; return oatom;
}

_atom F2C(const xmatrix<double>& lattice,const _atom& iatom) {
  return F2C(1.0,lattice,iatom); // VASP cartesian coordinates are intended normalized on the scale
}

_atom C2F(const double& scale,const xmatrix<double>& lattice,const _atom& iatom) {
  xmatrix<double> f2c(3,3);f2c=scale*trasp(lattice);
  _atom oatom(iatom); oatom.fpos=inverse(f2c)*iatom.cpos; return oatom;
}

_atom C2F(const xmatrix<double>& lattice,const _atom& iatom) {
  return C2F(1.0,lattice,iatom); // VASP cartesian coordinates are intended normalized on the scale
}

// **************************************************************************
// -------------------------------------------------------- for _atoms with structure
_atom F2C(const double& scale,const xstructure& str,const _atom& iatom) {
  xmatrix<double> f2c(3,3);f2c=scale*trasp(str.lattice);
  _atom oatom(iatom); oatom.cpos=f2c*iatom.fpos; return oatom;
}

_atom F2C(const xstructure& str,const _atom& iatom) {
  return F2C(1.0,str.lattice,iatom); // VASP cartesian coordinates are intended normalized on the scale
}

_atom C2F(const double& scale,const xstructure& str,const _atom& iatom) {
  xmatrix<double> f2c(3,3);f2c=scale*trasp(str.lattice);
  _atom oatom(iatom); oatom.fpos=inverse(f2c)*iatom.cpos; return oatom;
}

_atom C2F(const xstructure& str,const _atom& iatom) {
  return C2F(1.0,str.lattice,iatom); // VASP cartesian coordinates are intended normalized on the scale
}

// **************************************************************************
// **************************************************************************
// -------------------------------------------------------- for OPERATORS (matrices)
xmatrix<double> FF2CC(const double& scale,const xmatrix<double>& lattice,const xmatrix<double>& fmat) {
  xmatrix<double> f2c(3,3);f2c=scale*trasp(lattice);
  return f2c*fmat*inverse(f2c);                    // fmat is an operation in F coordinates
}

xmatrix<double> FF2CC(const xmatrix<double>& lattice,const xmatrix<double>& fmat) {
  return FF2CC(1.0,lattice,fmat); // VASP cartesian coordinates are intended normalized on the scale
}

xmatrix<double> CC2FF(const double& scale,const xmatrix<double>& lattice,const xmatrix<double>& cmat) {
  xmatrix<double> f2c(3,3);f2c=scale*trasp(lattice);
  return inverse(f2c)*cmat*f2c;                    // cmat is an operation in C coordinates
}

xmatrix<double> CC2FF(const xmatrix<double>& lattice,const xmatrix<double>& cmat) {
  return CC2FF(1.0,lattice,cmat);  // VASP cartesian coordinates are intended normalized on the scale
}

// ***************************************************************************
// Function IdenticalAtoms
// ***************************************************************************
// Makes all atoms the same. Stefano Curtarolo Nov 2008
void xstructure::IdenticalAtoms(void) {
  xvector<double> fpos(3),cpos(3);
  // fix atoms
  for(uint i=0;i<atoms.size();i++) {
    fpos=atoms.at(i).fpos;    // save fpos
    cpos=atoms.at(i).cpos;    // save cpos
    atoms.at(i)=atoms.at(0);  // use this so to copy all info
    atoms.at(i).fpos=fpos;    // fix back fpos
    atoms.at(i).cpos=cpos;    // fix back cpos
    atoms.at(i).partial_occupation_flag=FALSE;
    atoms.at(i).partial_occupation_value=1.0;
  }
  // fix numbers
  num_each_type.clear();num_each_type.push_back(atoms.size());
  comp_each_type.clear();comp_each_type.push_back((double) atoms.size());
  GetStoich();  //CO20170724
  // fix species
  for(uint i=1;i<species.size();i++) {
    species.at(i)="";species_pp.at(i)="";species_volume.at(i)=0.0;species_mass.at(i)=0.0;
  }
}

xstructure IdenticalAtoms(const xstructure& str) {
  xstructure sstr(str);
  sstr.IdenticalAtoms();
  return sstr;
}

// ***************************************************************************
// Function SwapCoordinates
// ***************************************************************************
// Permute Coordinates i with j Stefano Curtarolo Oct 2009
// Wahyu Setyawan: keep the same right-handness of the original
void xstructure::SwapCoordinates(const uint& ii,const uint& jj) {            // Permute Coordinates i with j
  // is obvious
  uint kk=0;
  if((ii==1 && jj==2) || (ii==2 && jj==1)) kk=3;   // permutation
  if((ii==2 && jj==3) || (ii==3 && jj==2)) kk=1;   // permutation
  if((ii==3 && jj==1) || (ii==1 && jj==3)) kk=2;   // permutation
  if((ii==jj) || (ii<1) || (ii>3) || (jj<1) || (jj>3)) return; // nothing to do

  // do the lattice
  double tmp;
  for(uint i=1;i<=3;i++) {
    tmp=lattice[ii][i];lattice[ii][i]=lattice[jj][i];lattice[jj][i]=tmp;   // swap RAW_ii with RAW_jj
    lattice[kk][i]=-lattice[kk][i]; // keep the right-handness of the original
  }
  FixLattices(); // fix the lattices ...
  // cartesian atoms are the same (crystal was not really changed, only definitions of abc were
  // so now I can get the new fractional
  for(uint i=0;i<atoms.size();i++)
    atoms.at(i).fpos=C2F(lattice,atoms.at(i).cpos); // DONE, now it has new fpos with the same cpos !
  // done bring all in cell now
  BringInCell();
}

xstructure SwapCoordinates(const xstructure& str,const uint& ii,const uint& jj) {
  xstructure sstr(str);
  sstr.SwapCoordinates(ii,jj);
  return sstr;
}

// ***************************************************************************
// Function GetLatticeType()
// ***************************************************************************
// Determine the real, reciprocal, and superlattice symmetry information
// Stefano Curtarolo
// Modified by David Hicks (DX)
// Includes self-consistency loop to ensure descriptions are commensurate
// DX20210225 - cleaned/consolidated function
void xstructure::GetLatticeType(double sym_eps, bool no_scan) {
  xstructure str_sp,str_sc;
  GetLatticeType(str_sp,str_sc,sym_eps,no_scan);
}

void xstructure::GetLatticeType(xstructure& str_sp,xstructure& str_sc, double sym_eps, bool no_scan) {

  bool LDEBUG=(FALSE || XHOST.DEBUG);

  // ---------------------------------------------------------------------------
  // set symmetry tolerance based on the following sequence
  // 1) use input, 2) use sym_eps in xstructure, 3) calculate default
  double tolerance = sym_eps;
  if(tolerance==AUROSTD_MAX_DOUBLE){
    if((*this).sym_eps_calculated){ tolerance = (*this).sym_eps; }
    else{ tolerance=SYM::defaultTolerance((*this)); }
  }
  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [1] Set symmetry tolerance (starting sym_eps=" << tolerance << ")" << endl; }

  // keep track of self-consistent tolerance
  bool same_eps = false;
  uint count = 0;
  uint count_max = 100; // safety for while loop, don't calculate forever

  // update tolerance info in *this
  (*this).sym_eps=tolerance;
  (*this).sym_eps_calculated=true;
  double tolerance_orig=tolerance; //DX20210623

  // ---------------------------------------------------------------------------
  // loop over the real, reciprocal, and superlattice analysis until all
  // symmetries are commensurate at a common tolerance value
  while(!same_eps && count++ < count_max){

    // ---------------------------------------------------------------------------
    // clear to start
    (*this).ClearSymmetry();

    // ---------------------------------------------------------------------------
    // update the tolerance, it may have change during loop
    tolerance = (*this).sym_eps;
    no_scan = (*this).sym_eps_no_scan;

    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [2] Top of self-consistent lattice-type loop (calculating real, reciprocal, and superlattice types) (sym_eps=" << tolerance << ", sym_eps_change_count=" << (*this).sym_eps_change_count << ")" << endl; }

    // ---------------------------------------------------------------------------
    // check if consistency checks failed (maxed while loop iteration)
    // turn of scan
    if(count==count_max){
      no_scan=(*this).sym_eps_no_scan=true;
      tolerance = tolerance_orig; // set to original tolerance //DX20210623 - originally sym_eps, but this could be AUROSTD_MAX_DOUBLE;
      cerr << __AFLOW_FUNC__ << " Unable to calculate consistent symmetry. Calculating at original tolerance (sym_eps=" << sym_eps << ") and ignoring consistency checks." << endl;
    }

    // ---------------------------------------------------------------------------
    // REAL - pass in str_sp and str_sc to keep primitive and conventional info
    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [3] Calculate real lattice type (sym_eps=" << tolerance << ", sym_eps_change_count=" << (*this).sym_eps_change_count << ")" << endl; }
    (*this).GetRealLatticeType(str_sp, str_sc, tolerance);
    tolerance = (*this).sym_eps; // update the tolerance

    // ---------------------------------------------------------------------------
    // RECIPROCAL
    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [4] Calculate reciprocal lattice type (sym_eps=" << tolerance << ", sym_eps_change_count=" << (*this).sym_eps_change_count << ")" << endl; }
    (*this).GetReciprocalLatticeType(tolerance);
    if(!no_scan && (*this).sym_eps != tolerance){ continue; } // if tolerance changed, recalc

    // ---------------------------------------------------------------------------
    // SUPERLATTICE
    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [5] Calculate superlattice type (sym_eps=" << tolerance << ", sym_eps_change_count=" << (*this).sym_eps_change_count << ")" << endl; }
    (*this).GetSuperlatticeType(tolerance);
    if(!no_scan && (*this).sym_eps != tolerance){ continue; } // if tolerance changed, recalc

    // made it to the end with same sym_eps
    same_eps = true;
  }

  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [6] Lattice types calculation finished! (sym_eps=" << tolerance << ", sym_eps_change_count=" << (*this).sym_eps_change_count << ")" << endl; }

}

//DX20210302 [OBSOLETE] // ***************************************************************************
//DX20210302 [OBSOLETE] // Function GetLatticeType
//DX20210302 [OBSOLETE] // ***************************************************************************
//DX20210302 [OBSOLETE] void xstructure::GetLatticeType(void) {
//DX20210302 [OBSOLETE]   xstructure str_sp,str_sc;
//DX20210302 [OBSOLETE]   GetLatticeType(str_sp,str_sc);
//DX20210302 [OBSOLETE] }
//DX20210302 [OBSOLETE] 
//DX20210302 [OBSOLETE] void xstructure::GetLatticeType(xstructure& str_sp,xstructure& str_sc) {
//DX20210302 [OBSOLETE]   //  bool VERBOSE=TRUE;
//DX20210302 [OBSOLETE]   bool VERBOSE=FALSE;
//DX20210302 [OBSOLETE]   //DX double eps=0.002,epsang=0.02;
//DX20210302 [OBSOLETE]   // double eps=0.02,epsang=0.02;  //JX
//DX20210302 [OBSOLETE]   // DIRECT
//DX20210302 [OBSOLETE]   xstructure str_in;//,str_sp,str_sc;
//DX20210302 [OBSOLETE]   // start
//DX20210302 [OBSOLETE]   str_in=*this;
//DX20210302 [OBSOLETE]   str_in.title="NO_RECURSION";
//DX20210302 [OBSOLETE]   // str_in.GetPrimitive();
//DX20210302 [OBSOLETE]   // str_in.MinkowskiBasisReduction();
//DX20210302 [OBSOLETE]   // cerr << str_in << endl;
//DX20210302 [OBSOLETE]   // str_in.CalculateSymmetryPointGroup(TRUE);
//DX20210302 [OBSOLETE]   // str_in.CalculateSymmetryFactorGroup(TRUE);
//DX20210302 [OBSOLETE]   // str_in.CalculateSymmetryPointGroupCrystal(TRUE);
//DX20210302 [OBSOLETE]   //DX START
//DX20210302 [OBSOLETE]   bool same_eps = false;
//DX20210302 [OBSOLETE]   uint count = 0;
//DX20210302 [OBSOLETE]   while(same_eps == false && count++ < 100){
//DX20210302 [OBSOLETE]     //DX END
//DX20210302 [OBSOLETE]     if(0) {
//DX20210302 [OBSOLETE]       LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc); // STD tolerance  // ONLY BRAVAIS_CRYSTAL
//DX20210302 [OBSOLETE]       // LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,0.0001,0.001);
//DX20210302 [OBSOLETE]     }
//DX20210302 [OBSOLETE]     if(1) {
//DX20210302 [OBSOLETE]       //    cerr << XPID << "LATTICE::Bravais_Lattice_StructureDefault IN" << endl;
//DX20210302 [OBSOLETE]       LATTICE::Bravais_Lattice_StructureDefault(str_in,str_sp,str_sc); // STD tolerance  // ONLY BRAVAIS_CRYSTAL
//DX20210302 [OBSOLETE]       //  cerr << XPID << "LATTICE::Bravais_Lattice_StructureDefault OUT" << endl;
//DX20210302 [OBSOLETE]     }
//DX20210302 [OBSOLETE]     if(VERBOSE) cerr << "xstructure::GetLatticeType: [4]" << endl;
//DX20210302 [OBSOLETE]     if(str_sp.pgroup_calculated==FALSE) str_sp.CalculateSymmetryPointGroup(FALSE);// cerr << "POINT GROUP" << endl;
//DX20210302 [OBSOLETE]     if(str_sp.fgroup_calculated==FALSE) str_sp.CalculateSymmetryFactorGroup(FALSE); //cerr << "FACTOR GROUP" << endl;
//DX20210302 [OBSOLETE]     if(str_sp.pgroup_xtal_calculated==FALSE) str_sp.CalculateSymmetryPointGroupCrystal(FALSE); //cerr << "POINT GROUP XTAL" << endl;
//DX20210302 [OBSOLETE]     //  *this=str_sp; // more obvious but will mess up the structures.... we only want to take the properties
//DX20210302 [OBSOLETE]     this->bravais_lattice_type=str_sp.bravais_lattice_type;
//DX20210302 [OBSOLETE]     this->bravais_lattice_variation_type=str_sp.bravais_lattice_variation_type;
//DX20210302 [OBSOLETE]     this->bravais_lattice_system=str_sp.bravais_lattice_system;
//DX20210302 [OBSOLETE]     this->bravais_lattice_lattice_type=str_sp.bravais_lattice_lattice_type;
//DX20210302 [OBSOLETE]     this->bravais_lattice_lattice_variation_type=str_sp.bravais_lattice_lattice_variation_type;
//DX20210302 [OBSOLETE]     this->bravais_lattice_lattice_system=str_sp.bravais_lattice_lattice_system;
//DX20210302 [OBSOLETE]     this->volume_changed_original2new=str_sp.volume_changed_original2new; //DX20181024
//DX20210302 [OBSOLETE]     this->transform_coordinates_original2new=str_sp.transform_coordinates_original2new; //DX20181024
//DX20210302 [OBSOLETE]     this->transform_coordinates_new2original=str_sp.transform_coordinates_new2original; //DX20181024
//DX20210302 [OBSOLETE]     this->rotate_lattice_original2new=str_sp.rotate_lattice_original2new; //DX20181024
//DX20210302 [OBSOLETE]     this->rotate_lattice_new2original=str_sp.rotate_lattice_new2original; //DX20181024
//DX20210302 [OBSOLETE]     this->pearson_symbol=str_sp.pearson_symbol;
//DX20210302 [OBSOLETE]     this->crystal_family=str_sp.crystal_family;
//DX20210302 [OBSOLETE]     this->crystal_system=str_sp.crystal_system;
//DX20210302 [OBSOLETE]     this->point_group_crystal_class=str_sp.point_group_crystal_class;
//DX20210302 [OBSOLETE]     this->point_group_Shoenflies=str_sp.point_group_Shoenflies;
//DX20210302 [OBSOLETE]     this->point_group_Hermann_Mauguin=str_sp.point_group_Hermann_Mauguin;
//DX20210302 [OBSOLETE]     this->point_group_orbifold=str_sp.point_group_orbifold;
//DX20210302 [OBSOLETE]     this->point_group_type=str_sp.point_group_type;
//DX20210302 [OBSOLETE]     this->point_group_order=str_sp.point_group_order;
//DX20210302 [OBSOLETE]     this->point_group_structure=str_sp.point_group_structure;
//DX20210302 [OBSOLETE]     // RECIPROCAL
//DX20210302 [OBSOLETE]     xstructure str_reciprocal_in,str_reciprocal_sp,str_reciprocal_sc;
//DX20210302 [OBSOLETE]     str_reciprocal_in.lattice=this->klattice;str_reciprocal_in.FixLattices();
//DX20210302 [OBSOLETE]     str_reciprocal_in.title="NO_RECURSION";
//DX20210302 [OBSOLETE]     //DX+CO START
//DX20210302 [OBSOLETE]     this->sym_eps=str_reciprocal_in.sym_eps=str_reciprocal_sp.sym_eps=str_reciprocal_sc.sym_eps=str_sp.sym_eps; //DX
//DX20210302 [OBSOLETE]     this->sym_eps_calculated=str_reciprocal_in.sym_eps_calculated=str_reciprocal_sp.sym_eps_calculated=str_reciprocal_sc.sym_eps_calculated=str_sp.sym_eps_calculated; //DX
//DX20210302 [OBSOLETE]     this->sym_eps_change_count=str_reciprocal_in.sym_eps_change_count=str_reciprocal_sp.sym_eps_change_count=str_reciprocal_sc.sym_eps_change_count=str_sp.sym_eps_change_count; //DX20180222 - added sym_eps change count
//DX20210302 [OBSOLETE]     //DX+CO END
//DX20210302 [OBSOLETE]     _atom atom;str_reciprocal_in.AddAtom(atom,false); //CO20230319 - add by type
//DX20210302 [OBSOLETE]     //LATTICE::Standard_Lattice_Structure(str_reciprocal_in,str_reciprocal_sp,str_reciprocal_sc,eps,epsang); //SC OLD VERSION
//DX20210302 [OBSOLETE]     if(VERBOSE) cerr << "xstructure::GetLatticeType: [5]" << endl;
//DX20210302 [OBSOLETE]     //DX int ss=0; //JX
//DX20210302 [OBSOLETE]     //DX LATTICE::Standard_Lattice_Structure(str_reciprocal_in,str_reciprocal_sp,str_reciprocal_sc,eps,epsang,ss,_EPS_); //JX
//DX20210302 [OBSOLETE]     //DX20170814 START - Use real pgroup to calculate pgroupk and then set pgroupk from str_sp to the pgroup and pgroup_xtal of str_reciprocal_in
//DX20210302 [OBSOLETE]     //DX20170814 The pgroup and pgroup_xtal are the same for the str_reciprocal structure because there is only one atom at the origin
//DX20210302 [OBSOLETE]     //DX20170814 (i.e. lattice and crystal symmetry are the same for the reciprocal space crystal)
//DX20210302 [OBSOLETE]     //DX20170829 [OBSOLETE] -since performing full symmetry analysis by default - str_sp.CalculateSymmetryPointGroupKLattice(FALSE);
//DX20210302 [OBSOLETE]     //DX20180426 - possible that lattice exhibits lower symmetry than crystal (i.e., from str_sp); would need to pass lattice symmetry from Standard_Lattice, but that information is not stored out of scope, commenting out 5 lines below 
//DX20210302 [OBSOLETE]     //DX20180426 [OBSOLETE] - possible that lattice exhibits lower symmetry than crystal (i.e., from str_sp) - str_reciprocal_in.pgroup=str_reciprocal_sp.pgroup=str_reciprocal_sc.pgroup=str_sp.pgroupk;
//DX20210302 [OBSOLETE]     //DX20180426 [OBSOLETE] - possible that lattice exhibits lower symmetry than crystal (i.e., from str_sp) - str_reciprocal_in.pgroup_calculated=str_reciprocal_sp.pgroup_calculated=str_reciprocal_sc.pgroup_calculated=str_sp.pgroupk_calculated;
//DX20210302 [OBSOLETE]     //DX20180426 [OBSOLETE] - possible that lattice exhibits lower symmetry than crystal (i.e., from str_sp) - str_reciprocal_in.pgroup_xtal=str_reciprocal_sp.pgroup_xtal=str_reciprocal_sc.pgroup_xtal=str_sp.pgroupk;
//DX20210302 [OBSOLETE]     //DX20180426 [OBSOLETE] - possible that lattice exhibits lower symmetry than crystal (i.e., from str_sp) - str_reciprocal_in.pgroup_xtal_calculated=str_reciprocal_sp.pgroup_xtal_calculated=str_reciprocal_sc.pgroup_xtal_calculated=str_sp.pgroup_calculated;
//DX20210302 [OBSOLETE]     //DX20180426 [OBSOLETE] - possible that lattice exhibits lower symmetry than crystal (i.e., from str_sp) - str_reciprocal_in.pgroup_xtal_calculated=str_reciprocal_sp.pgroup_xtal_calculated=str_reciprocal_sc.pgroup_xtal_calculated=str_sp.pgroup_calculated;
//DX20210302 [OBSOLETE]     //DX20170814 END
//DX20210302 [OBSOLETE]     LATTICE::Standard_Lattice_StructureDefault(str_reciprocal_in,str_reciprocal_sp,str_reciprocal_sc,false); //DX //DX20180226 - do not need to do full sym for recip
//DX20210302 [OBSOLETE]     //DX START
//DX20210302 [OBSOLETE]     if(str_sp.sym_eps == str_reciprocal_sp.sym_eps){
//DX20210302 [OBSOLETE]       same_eps = true;
//DX20210302 [OBSOLETE]     }
//DX20210302 [OBSOLETE]     else {
//DX20210302 [OBSOLETE]       str_in.sym_eps = str_sp.sym_eps = str_sc.sym_eps = str_reciprocal_sp.sym_eps;
//DX20210302 [OBSOLETE]       str_in.sym_eps_change_count = str_sp.sym_eps_change_count = str_sc.sym_eps_change_count = str_reciprocal_sp.sym_eps_change_count; //DX20180222 - added sym_eps change count
//DX20210302 [OBSOLETE]     }
//DX20210302 [OBSOLETE]     //DX END
//DX20210302 [OBSOLETE]     this->reciprocal_lattice_type=str_reciprocal_sp.bravais_lattice_type;
//DX20210302 [OBSOLETE]     this->reciprocal_lattice_variation_type=str_reciprocal_sp.bravais_lattice_variation_type;
//DX20210302 [OBSOLETE]     if(VERBOSE) cerr << "xstructure::GetLatticeType: [6]" << endl;
//DX20210302 [OBSOLETE]     // SUPERLATTICE
//DX20210302 [OBSOLETE]     xstructure str_superlattice_in,str_superlattice_sp,str_superlattice_sc;
//DX20210302 [OBSOLETE]     str_superlattice_in=*this;
//DX20210302 [OBSOLETE]     str_superlattice_in.ClearSymmetry();  //DX20170814 - It wasn't cleared, so nothing was being calculated
//DX20210302 [OBSOLETE]     str_superlattice_in.title="NO_RECURSION";
//DX20210302 [OBSOLETE]     //str_superlattice_in.GetPrimitive(0.01);
//DX20210302 [OBSOLETE]     if(VERBOSE) cerr << str_superlattice_in << endl;
//DX20210302 [OBSOLETE]     if(VERBOSE) cerr << "xstructure::GetLatticeType: [7]" << endl;
//DX20210302 [OBSOLETE]     str_superlattice_in.IdenticalAtoms();  // make superlattice
//DX20210302 [OBSOLETE]     if(VERBOSE) cerr << "xstructure::GetLatticeType: [8]" << endl;
//DX20210302 [OBSOLETE]     if(VERBOSE) cerr << str_superlattice_in << endl;
//DX20210302 [OBSOLETE]     str_superlattice_in.GetPrimitive(0.005);
//DX20210302 [OBSOLETE]     if(VERBOSE) cerr << str_superlattice_in << endl;
//DX20210302 [OBSOLETE]     if(VERBOSE) cerr << "xstructure::GetLatticeType: [9]" << endl;
//DX20210302 [OBSOLETE]     str_superlattice_in.Minkowski_calculated=FALSE;
//DX20210302 [OBSOLETE]     if(VERBOSE) cerr << "xstructure::GetLatticeType: [10]" << endl;
//DX20210302 [OBSOLETE]     str_superlattice_in.MinkowskiBasisReduction();
//DX20210302 [OBSOLETE]     if(VERBOSE) cerr << "xstructure::GetLatticeType: [11]" << endl;
//DX20210302 [OBSOLETE]     if(VERBOSE) cerr << str_superlattice_in << endl;
//DX20210302 [OBSOLETE]     //  LATTICE::Standard_Lattice_Structure(str_superlattice_in,str_superlattice_sp,str_superlattice_sc); //SC OLD VERSION
//DX20210302 [OBSOLETE]     //DX ss=0; //JX
//DX20210302 [OBSOLETE]     //DX+CO START
//DX20210302 [OBSOLETE]     str_superlattice_in.sym_eps=str_superlattice_sp.sym_eps=str_superlattice_sc.sym_eps=str_sp.sym_eps; //DX
//DX20210302 [OBSOLETE]     str_superlattice_in.sym_eps_calculated=str_superlattice_sp.sym_eps_calculated=str_superlattice_sc.sym_eps_calculated=str_sp.sym_eps_calculated; //DX
//DX20210302 [OBSOLETE]     str_superlattice_in.sym_eps_change_count=str_superlattice_sp.sym_eps_change_count=str_superlattice_sc.sym_eps_change_count=str_sp.sym_eps_change_count; //DX20180222 - added sym_eps change count
//DX20210302 [OBSOLETE]     //DX+CO END
//DX20210302 [OBSOLETE]     //DX LATTICE::Standard_Lattice_Structure(str_superlattice_in,str_superlattice_sp,str_superlattice_sc,eps,epsang,ss,_EPS_); //JX
//DX20210302 [OBSOLETE]     LATTICE::Standard_Lattice_StructureDefault(str_superlattice_in,str_superlattice_sp,str_superlattice_sc,false); //DX //DX20180226 - do not need to do full sym for superlattice
//DX20210302 [OBSOLETE]     //DX START
//DX20210302 [OBSOLETE]     if(str_sp.sym_eps == str_superlattice_sp.sym_eps){
//DX20210302 [OBSOLETE]       same_eps = true;
//DX20210302 [OBSOLETE]     }
//DX20210302 [OBSOLETE]     else {
//DX20210302 [OBSOLETE]       str_sp.sym_eps = str_superlattice_sp.sym_eps;
//DX20210302 [OBSOLETE]       str_sp.sym_eps_change_count = str_superlattice_sp.sym_eps_change_count; //DX20180222 - added sym_eps change count
//DX20210302 [OBSOLETE]     }
//DX20210302 [OBSOLETE]     //DX END
//DX20210302 [OBSOLETE]     if(VERBOSE) cerr << "xstructure::GetLatticeType: [12]" << endl;
//DX20210302 [OBSOLETE]     this->bravais_superlattice_lattice=str_superlattice_sp.lattice; //DX20210209
//DX20210302 [OBSOLETE]     this->bravais_superlattice_type=str_superlattice_sp.bravais_lattice_type;
//DX20210302 [OBSOLETE]     this->bravais_superlattice_variation_type=str_superlattice_sp.bravais_lattice_variation_type;
//DX20210302 [OBSOLETE]     this->bravais_superlattice_system=str_superlattice_sp.bravais_lattice_system;
//DX20210302 [OBSOLETE]     this->pearson_symbol_superlattice=str_superlattice_sp.pearson_symbol;
//DX20210302 [OBSOLETE]     if(count==100){
//DX20210302 [OBSOLETE]       cerr << "ERROR in Bravais_Lattice_StructureDefault(): Unable to find reliable sym_eps." << endl;
//DX20210302 [OBSOLETE]       break;
//DX20210302 [OBSOLETE]     }
//DX20210302 [OBSOLETE]   } //DX while loop
//DX20210302 [OBSOLETE] }

// ***************************************************************************
// Function GetExtendedCrystallographicData() //DX20210302
// ***************************************************************************
// Determine the real, reciprocal, superlattice, and space group symmetries
// Includes self-consistency loop to ensure descriptions are commensurate
void xstructure::GetExtendedCrystallographicData(double sym_eps,
    bool no_scan,
    int setting) {
  xstructure str_sp,str_sc;
  GetExtendedCrystallographicData(str_sp,str_sc,sym_eps,no_scan,setting);
}

void xstructure::GetExtendedCrystallographicData(xstructure& str_sp,
    xstructure& str_sc,
    double sym_eps,
    bool no_scan,
    int setting) {

  bool LDEBUG=(FALSE || XHOST.DEBUG);

  // ---------------------------------------------------------------------------
  // set symmetry tolerance based on the following sequence
  // 1) use input, 2) use sym_eps in xstructure, 3) calculate default
  double tolerance = sym_eps;
  if(tolerance==AUROSTD_MAX_DOUBLE){
    if((*this).sym_eps_calculated){ tolerance = (*this).sym_eps; }
    else{ tolerance=SYM::defaultTolerance((*this)); }
  }
  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [1] Set symmetry tolerance (starting sym_eps=" << tolerance << ")" << endl; }

  // keep track of self-consistent tolerance
  bool force_perform = true;
  bool same_eps = false;
  uint count = 0;
  uint count_max = 100; // safety for while loop, don't calculate forever

  // update tolerance info in *this
  (*this).sym_eps=tolerance;
  (*this).sym_eps_calculated=true;
  double tolerance_orig=tolerance; //DX20210623

  // ---------------------------------------------------------------------------
  // loop over the real, reciprocal, and superlattice analysis until all
  // symmetries are commensurate with a common tolerance value
  while(!same_eps && count++ < count_max){

    // ---------------------------------------------------------------------------
    // clear to start
    (*this).ClearSymmetry();

    // ---------------------------------------------------------------------------
    // update the tolerance, it may have change during loop
    tolerance = (*this).sym_eps;
    no_scan = (*this).sym_eps_no_scan;

    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [2] Top of self-consistent extended crystallographic data loop (calculating lattice type and space group data) (sym_eps=" << tolerance << ", sym_eps_change_count=" << (*this).sym_eps_change_count << ")" << endl; }

    // ---------------------------------------------------------------------------
    // check if consistency checks failed (maxed while loop iteration)
    // turn off scan
    if(count==count_max){
      no_scan=(*this).sym_eps_no_scan=true;
      tolerance = tolerance_orig; // set to original tolerance //DX20210623 - originally sym_eps, but this could be AUROSTD_MAX_DOUBLE;
      cerr << __AFLOW_FUNC__ << " Unable to calculate consistent symmetry. Calculating at original tolerance (sym_eps=" << sym_eps << ") and ignoring consistency checks." << endl;
    }

    // ---------------------------------------------------------------------------
    // REAL, RECIPROCAL, and SUPERLATTICE data
    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [3] Calculate real, reciprocal, and superlattice information (sym_eps=" << tolerance << ", sym_eps_change_count=" << (*this).sym_eps_change_count << ")" << endl; }
    (*this).GetLatticeType(str_sp, str_sc, (*this).sym_eps);
    tolerance = (*this).sym_eps; // update the tolerance

    // ---------------------------------------------------------------------------
    // space group data
    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [4] Calculate the space group symmetry information (sym_eps=" << tolerance << ", sym_eps_change_count=" << (*this).sym_eps_change_count << ")" << endl; }
    (*this).SpaceGroup_ITC(tolerance, -1, setting, no_scan);
    if(!no_scan && (*this).sym_eps != tolerance){ continue; } // if tolerance changed, recalc

    // ---------------------------------------------------------------------------
    // check GetLatticeType vs SpaceGroup_ITC results
    if(!(*this).sym_eps_no_scan) {
      int multiplicity_of_primitive=str_sp.fgroup.size()/str_sp.pgroup_xtal.size();
      bool derivative_structure=false;
      string lattice_and_centering = LATTICE::Lattice2TypeAndCentering((*this).bravais_lattice_type); //DX20210412 - check centering
      string lattice_and_centering_from_sg = SYM::spacegroup2latticeAndCentering((*this).space_group_ITC); //DX20210412 - check centering
      if(!(lattice_and_centering == lattice_and_centering_from_sg && SYM::ComparePointGroupAndSpaceGroupString((*this),multiplicity_of_primitive,derivative_structure))){
        if(LDEBUG) {
          cerr << __AFLOW_FUNC__ << " WARNING: Space group symbol and point group symbol do not match. (sg=" << GetSpaceGroupName((*this).space_group_ITC,(*this).directory) << ", centering_sg=" << lattice_and_centering_from_sg << " | pg=" << (*this).point_group_Hermann_Mauguin << ", centering=" << lattice_and_centering << ") [dir=" << (*this).directory << "]" << endl;
        }
        if(!SYM::change_tolerance((*this),(*this).sym_eps,(*this).dist_nn_min,(*this).sym_eps_no_scan)){
          if(force_perform){
            if(LDEBUG) {
              cerr << __AFLOW_FUNC__ << " WARNING: Scan failed. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << (*this).directory << endl;
            }
            no_scan = (*this).sym_eps_no_scan = true; //DX20210331
            (*this).sym_eps = tolerance_orig; // set to original tolerance //DX20210623 - originally sym_eps, but this could be AUROSTD_MAX_DOUBLE;
          }
        }
        continue;
      }
    }

    // made it to the end with same sym_eps
    same_eps = true;
  }

  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [5] Extended crystallographic data calculation finished! (sym_eps=" << tolerance << ", sym_eps_change_count=" << count << ")" << endl; }

}

// ***************************************************************************
// Function GetRealLatticeType //DX20210209
// ***************************************************************************
// Determine the real space symmetry information
// Includes self-consistency loop to ensure descriptions are commensurate
void xstructure::GetRealLatticeType(double sym_eps) {
  xstructure str_sp,str_sc;
  GetRealLatticeType(str_sp,str_sc,sym_eps);
}

void xstructure::GetRealLatticeType(xstructure& str_sp,xstructure& str_sc, double sym_eps) {

  bool LDEBUG=(FALSE || XHOST.DEBUG);

  // ---------------------------------------------------------------------------
  // set symmetry tolerance based on the following sequence
  // 1) use input, 2) use sym_eps in xstructure, 3) calculate default
  double tolerance = sym_eps;
  if(tolerance==AUROSTD_MAX_DOUBLE){
    if((*this).sym_eps_calculated){ tolerance = (*this).sym_eps; }
    else{ tolerance=SYM::defaultTolerance((*this)); }
  }
  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [1] Set symmetry tolerance (starting sym_eps=" << tolerance << ")" << endl; }

  // need to create a copy
  xstructure str_in=*this;

  // update tolerance info for all xstructures
  str_in.sym_eps=str_sp.sym_eps=str_sc.sym_eps=tolerance;
  str_in.sym_eps_calculated=str_sp.sym_eps_calculated=str_sc.sym_eps_calculated=(*this).sym_eps_calculated;
  str_in.sym_eps_change_count=str_sp.sym_eps_change_count=str_sc.sym_eps_change_count=(*this).sym_eps_change_count;
  str_in.sym_eps_no_scan=str_sp.sym_eps_no_scan=str_sc.sym_eps_no_scan=(*this).sym_eps_no_scan; //DX20210430

  // calculate
  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [2]" << endl; }
  LATTICE::Bravais_Lattice_StructureDefault(str_in,str_sp,str_sc); // STD tolerance  // ONLY BRAVAIS_CRYSTAL

  // set properties
  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [3]" << endl; }
  if(str_sp.pgroup_calculated==FALSE) str_sp.CalculateSymmetryPointGroup(FALSE);// cerr << "POINT GROUP" << endl;
  if(str_sp.fgroup_calculated==FALSE) str_sp.CalculateSymmetryFactorGroup(FALSE); //cerr << "FACTOR GROUP" << endl;
  if(str_sp.pgroup_xtal_calculated==FALSE) str_sp.CalculateSymmetryPointGroupCrystal(FALSE); //cerr << "POINT GROUP XTAL" << endl;
  //  *this=str_sp; // more obvious but will mess up the structures.... we only want to take the properties
  this->bravais_lattice_type=str_sp.bravais_lattice_type;
  this->bravais_lattice_variation_type=str_sp.bravais_lattice_variation_type;
  this->bravais_lattice_system=str_sp.bravais_lattice_system;
  this->bravais_lattice_lattice_type=str_sp.bravais_lattice_lattice_type;
  this->bravais_lattice_lattice_variation_type=str_sp.bravais_lattice_lattice_variation_type;
  this->bravais_lattice_lattice_system=str_sp.bravais_lattice_lattice_system;
  this->volume_changed_original2new=str_sp.volume_changed_original2new; //DX20181024
  this->transform_coordinates_original2new=str_sp.transform_coordinates_original2new; //DX20181024
  this->transform_coordinates_new2original=str_sp.transform_coordinates_new2original; //DX20181024
  this->rotate_lattice_original2new=str_sp.rotate_lattice_original2new; //DX20181024
  this->rotate_lattice_new2original=str_sp.rotate_lattice_new2original; //DX20181024
  this->pearson_symbol=str_sp.pearson_symbol;
  this->crystal_family=str_sp.crystal_family;
  this->crystal_system=str_sp.crystal_system;
  this->point_group_crystal_class=str_sp.point_group_crystal_class;
  this->point_group_Shoenflies=str_sp.point_group_Shoenflies;
  this->point_group_Hermann_Mauguin=str_sp.point_group_Hermann_Mauguin;
  this->point_group_orbifold=str_sp.point_group_orbifold;
  this->point_group_type=str_sp.point_group_type;
  this->point_group_order=str_sp.point_group_order;
  this->point_group_structure=str_sp.point_group_structure;
  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [4] DONE" << endl; }

  // update sym_eps
  this->sym_eps=str_in.sym_eps=str_sc.sym_eps=str_sp.sym_eps; //DX
  this->sym_eps_calculated=str_in.sym_eps_calculated=str_sc.sym_eps_calculated=str_sp.sym_eps_calculated; //DX
  this->sym_eps_change_count=str_in.sym_eps_change_count=str_sc.sym_eps_change_count=str_sp.sym_eps_change_count; //DX20180222 - added sym_eps change count
  this->sym_eps_no_scan=str_in.sym_eps_no_scan=str_sc.sym_eps_no_scan=str_sp.sym_eps_no_scan; //DX20210430 - added no_scan

}

// ***************************************************************************
// Function GetReciprocalLatticeType //DX20210209
// ***************************************************************************
// Determine the reciprocal space symmetry information
// Includes self-consistency loop to ensure descriptions are commensurate
void xstructure::GetReciprocalLatticeType(double sym_eps) {
  xstructure str_sp,str_sc;
  GetReciprocalLatticeType(str_sp,str_sc,sym_eps);
}

void xstructure::GetReciprocalLatticeType(xstructure& str_sp,xstructure& str_sc, double sym_eps) {

  bool LDEBUG=(FALSE || XHOST.DEBUG);

  // ---------------------------------------------------------------------------
  // set symmetry tolerance based on the following sequence
  // 1) use input, 2) use sym_eps in xstructure, 3) calculate default
  double tolerance = sym_eps;
  if(tolerance==AUROSTD_MAX_DOUBLE){
    if((*this).sym_eps_calculated){ tolerance = (*this).sym_eps; }
    else{ tolerance=SYM::defaultTolerance((*this)); }
  }
  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [1] Set symmetry tolerance (starting sym_eps=" << tolerance << ")" << endl; }

  // ---------------------------------------------------------------------------
  // RECIPROCAL - use klattice an one atom (at the origin)
  xstructure str_in;
  str_in.lattice=this->klattice;str_in.FixLattices();
  _atom atom;str_in.AddAtom(atom,false);  //CO20230319 - add by type

  // update tolerance info for all xstructures
  str_in.sym_eps=str_sp.sym_eps=str_sc.sym_eps=tolerance;
  str_in.sym_eps_calculated=str_sp.sym_eps_calculated=str_sc.sym_eps_calculated=(*this).sym_eps_calculated;
  str_in.sym_eps_change_count=str_sp.sym_eps_change_count=str_sc.sym_eps_change_count=(*this).sym_eps_change_count;
  str_in.sym_eps_no_scan=str_sp.sym_eps_no_scan=str_sc.sym_eps_no_scan=(*this).sym_eps_no_scan; //DX20210430 - added no_scan

  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [1]" << endl; }
  //DX20170814 START - Use real pgroup to calculate pgroupk and then set pgroupk from str_sp to the pgroup and pgroup_xtal of str_reciprocal_in
  //DX20170814 The pgroup and pgroup_xtal are the same for the str_reciprocal structure because there is only one atom at the origin
  //DX20170814 (i.e. lattice and crystal symmetry are the same for the reciprocal space crystal)
  //DX20170829 [OBSOLETE] -since performing full symmetry analysis by default - str_sp.CalculateSymmetryPointGroupKLattice(FALSE);
  //DX20180426 - possible that lattice exhibits lower symmetry than crystal (i.e., from str_sp); would need to pass lattice symmetry from Standard_Lattice, but that information is not stored out of scope, commenting out 5 lines below
  //DX20180426 [OBSOLETE] - possible that lattice exhibits lower symmetry than crystal (i.e., from str_sp) - str_reciprocal_in.pgroup=str_reciprocal_sp.pgroup=str_reciprocal_sc.pgroup=str_sp.pgroupk;
  //DX20180426 [OBSOLETE] - possible that lattice exhibits lower symmetry than crystal (i.e., from str_sp) - str_reciprocal_in.pgroup_calculated=str_reciprocal_sp.pgroup_calculated=str_reciprocal_sc.pgroup_calculated=str_sp.pgroupk_calculated;
  //DX20180426 [OBSOLETE] - possible that lattice exhibits lower symmetry than crystal (i.e., from str_sp) - str_reciprocal_in.pgroup_xtal=str_reciprocal_sp.pgroup_xtal=str_reciprocal_sc.pgroup_xtal=str_sp.pgroupk;
  //DX20180426 [OBSOLETE] - possible that lattice exhibits lower symmetry than crystal (i.e., from str_sp) - str_reciprocal_in.pgroup_xtal_calculated=str_reciprocal_sp.pgroup_xtal_calculated=str_reciprocal_sc.pgroup_xtal_calculated=str_sp.pgroup_calculated;
  //DX20180426 [OBSOLETE] - possible that lattice exhibits lower symmetry than crystal (i.e., from str_sp) - str_reciprocal_in.pgroup_xtal_calculated=str_reciprocal_sp.pgroup_xtal_calculated=str_reciprocal_sc.pgroup_xtal_calculated=str_sp.pgroup_calculated;
  //DX20170814 END
  LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc,false); //DX //DX20180226 - do not need to do full sym for recip

  this->reciprocal_lattice_type=str_sp.bravais_lattice_type;
  this->reciprocal_lattice_variation_type=str_sp.bravais_lattice_variation_type;
  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [2]" << endl; }

  // update sym_eps
  this->sym_eps=str_in.sym_eps=str_sc.sym_eps=str_sp.sym_eps; //DX
  this->sym_eps_calculated=str_in.sym_eps_calculated=str_sc.sym_eps_calculated=str_sp.sym_eps_calculated; //DX
  this->sym_eps_change_count=str_in.sym_eps_change_count=str_sc.sym_eps_change_count=str_sp.sym_eps_change_count; //DX20180222 - added sym_eps change count
  this->sym_eps_no_scan=str_in.sym_eps_no_scan=str_sc.sym_eps_no_scan=str_sp.sym_eps_no_scan; //DX20210430 - added no_scan
}

// ***************************************************************************
// Function GetSuperlatticeType() //DX20210302
// ***************************************************************************
// Determine the superlattice symmetry information
// Includes self-consistency loop to ensure descriptions are commensurate
void xstructure::GetSuperlatticeType(double sym_eps) {
  xstructure str_sp,str_sc;
  GetSuperlatticeType(str_sp,str_sc,sym_eps);
}

void xstructure::GetSuperlatticeType(xstructure& str_sp,xstructure& str_sc, double sym_eps) {

  bool LDEBUG=(FALSE || XHOST.DEBUG);

  // ---------------------------------------------------------------------------
  // set symmetry tolerance based on the following sequence
  // 1) use input, 2) use sym_eps in xstructure, 3) calculate default
  double tolerance = sym_eps;
  if(tolerance==AUROSTD_MAX_DOUBLE){
    if((*this).sym_eps_calculated){ tolerance = (*this).sym_eps; }
    else{ tolerance=SYM::defaultTolerance((*this)); }
  }
  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [1] Set symmetry tolerance (starting sym_eps=" << tolerance << ")" << endl; }

  // ---------------------------------------------------------------------------
  // SUPERLATTICE - decorate with single atom type
  xstructure str_in=*this;
  str_in.ClearSymmetry();  // need to clear symmetry; otherwise, nothing is calculated
  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " [1]" << endl;
    cerr << str_in << endl;
  }
  // decorate with single atom type
  str_in.IdenticalAtoms();  // make superlattice
  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " [2]" << endl;
    cerr << str_in << endl;
  }
  // primitivize
  str_in.GetPrimitive(); //DX20210430 - remove obsolete eps=0.005
  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " [3]" << endl;
    cerr << str_in << endl;
  }
  // Minkowski
  str_in.Minkowski_calculated=FALSE;
  str_in.MinkowskiBasisReduction();
  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " [4]" << endl;
    cerr << str_in << endl;
  }

  // update tolerance info for all xstructures
  str_in.sym_eps=str_sp.sym_eps=str_sc.sym_eps=tolerance;
  str_in.sym_eps_calculated=str_sp.sym_eps_calculated=str_sc.sym_eps_calculated=(*this).sym_eps_calculated;
  str_in.sym_eps_change_count=str_sp.sym_eps_change_count=str_sc.sym_eps_change_count=(*this).sym_eps_change_count;
  str_in.sym_eps_no_scan=str_sp.sym_eps_no_scan=str_sc.sym_eps_no_scan=(*this).sym_eps_no_scan; //DX20210430 - added no_scan

  // main lattice function
  LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc,false); //DX //DX20180226 - do not need to do full sym for superlattice
  str_sp.ReScale(1.0); //DX20210211 - need to rescale to 1 since we aren't propagating the superlattice scaling factor
  this->bravais_superlattice_lattice=str_sp.lattice;
  this->bravais_superlattice_type=str_sp.bravais_lattice_type;
  this->bravais_superlattice_variation_type=str_sp.bravais_lattice_variation_type;
  this->bravais_superlattice_system=str_sp.bravais_lattice_system;
  this->pearson_symbol_superlattice=str_sp.pearson_symbol;

  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [6] DONE" << endl; }

  // update sym_eps
  this->sym_eps=str_in.sym_eps=str_sc.sym_eps=str_sp.sym_eps; //DX
  this->sym_eps_calculated=str_in.sym_eps_calculated=str_sc.sym_eps_calculated=str_sp.sym_eps_calculated; //DX
  this->sym_eps_change_count=str_in.sym_eps_change_count=str_sc.sym_eps_change_count=str_sp.sym_eps_change_count; //DX20180222 - added sym_eps change count
  this->sym_eps_no_scan=str_in.sym_eps_no_scan=str_sc.sym_eps_no_scan=str_sp.sym_eps_no_scan; //DX20210430 - added no_scan
}

string GetLatticeType(xmatrix<double> lattice) {
  //DX double eps=0.00010,epsang=0.01;
  xstructure str_in,str_sp,str_sc;
  str_in.lattice=lattice;str_in.FixLattices();
  str_in.title="NO_RECURSION";
  _atom atom;str_in.AddAtom(atom,false);  //CO20230319 - add by type
  // LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,eps,epsang); //SC OLD VERSION
  //DX int ss=0; //JX
  //DX LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,eps,epsang,ss,_EPS_); //JX
  LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc); //DX
  return str_sp.bravais_lattice_type;
  //  return str_sp.bravais_lattice_variation_type;
}

string GetLatticeType(xvector<double> data) {
  xmatrix<double> lattice(GetClat(data));
  return GetLatticeType(lattice);
}

// ***************************************************************************
// Function Standard_Primitive
// ***************************************************************************
// Lattice Reduction to Max Orthogonality (MINK) and then Niggly Form
xstructure Standard_Primitive_UnitCellForm(const xstructure& a) {
  xstructure str_in(a),str_sp,str_sc;
  if(str_in.Standard_Lattice_avoid==TRUE) return str_in;     // Nothing to do
  if(str_in.Standard_Lattice_primitive==TRUE) return str_in; // already Primitive
  LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc);
  return str_sp;
}

xstructure GetStandardPrimitive(const xstructure& a) {
  // cerr << "GetStandardPrimitive(const xstructure& a)" << endl;
  return Standard_Primitive_UnitCellForm(a);
}

void xstructure::Standard_Primitive_UnitCellForm(void) {
  xstructure str_sp,str_sc;
  if(Standard_Lattice_avoid==TRUE) return;     // Nothing to do
  if(Standard_Lattice_primitive==TRUE) return; // already Primitive
  // cerr << bravais_lattice_type << endl;
  LATTICE::Standard_Lattice_StructureDefault(*this,str_sp,str_sc);
  *this=str_sp;
}

void xstructure::GetStandardPrimitive(void) {
  Standard_Primitive_UnitCellForm();
}

xmatrix<double> GetStandardPrimitive(xmatrix<double> lattice) {
  xstructure str;
  str.lattice=lattice;str.scale=1.0;str.FixLattices();
  str.title="NO_RECURSION";
  _atom atom;str.AddAtom(atom,false); //CO20230319 - add by type
  str=GetStandardPrimitive(str);
  return str.lattice;
}

xvector<double> GetStandardPrimitive(xvector<double> data) {
  xmatrix<double> lattice(GetClat(data));
  return Getabc_angles(GetStandardPrimitive(lattice),DEGREES);
}

// ***************************************************************************
// Function Standard_Conventional
// ***************************************************************************
// Lattice Reduction to Max Orthogonality (MINK) and then Niggly Form
xstructure Standard_Conventional_UnitCellForm(const xstructure& a) {
  xstructure str_in(a),str_sp,str_sc;
  if(str_in.Standard_Lattice_avoid==TRUE) return str_in;        // Nothing to do
  if(str_in.Standard_Lattice_conventional==TRUE) return str_in; // already Conventional
  LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc);
  return str_sc;
}

xstructure GetStandardConventional(const xstructure& a) {
  return Standard_Conventional_UnitCellForm(a);
}

void xstructure::Standard_Conventional_UnitCellForm(void) {
  xstructure str_sp,str_sc;
  if(Standard_Lattice_avoid==TRUE) return;        // Nothing to do
  if(Standard_Lattice_conventional==TRUE) return; // already Conventional
  LATTICE::Standard_Lattice_StructureDefault(*this,str_sp,str_sc);
  *this=str_sc;
}

void xstructure::GetStandardConventional(void) {
  Standard_Conventional_UnitCellForm();
}

xmatrix<double> GetStandardConventional(xmatrix<double> lattice) {
  xstructure str;
  str.lattice=lattice;str.scale=1.0;str.FixLattices();
  str.title="NO_RECURSION";
  _atom atom;str.AddAtom(atom,false); //CO20230319 - add by type
  str=GetStandardConventional(str);
  return str.lattice;
}

xvector<double> GetStandardConventional(xvector<double> data) {
  xmatrix<double> lattice(GetClat(data));
  return Getabc_angles(GetStandardConventional(lattice),DEGREES);
}

// ***************************************************************************
// Function SpeciesLabel
// ***************************************************************************
// Returns the name of the specie (the name of the 1st atom of the specie
// Stefano Curtarolo - nov 2008
string xstructure::SpeciesLabel(const uint& A) {
  if(A>num_each_type.size()) return "xxx (outside boundaries)"; // outside boundaries
  string name;vector<string> tokens;

  uint i=0,A_start=0;//,A_stop=0;
  for(i=0;i<num_each_type.size();i++) {
    if(i<A)  {A_start+=num_each_type.at(i);}
    // [UNNECESSARY] if(i==A) {A_start=A_start;}
    //    if(i==A) {A_stop=A_start+num_each_type.at(i);}
  }
  if(atoms.at(A_start).name_is_given==FALSE) return "xxx (name not given)"; // name not given
  name=atoms.at(A_start).name;
  aurostd::StringSubst(name,"+"," ");aurostd::StringSubst(name,"-"," ");
  aurostd::string2tokens(name,tokens," ");
  name=tokens[0];
  return name;
}

string SpeciesLabel(const xstructure& str,const uint& A) {
  xstructure sstr(str);
  return sstr.SpeciesLabel(A);
}

// ***************************************************************************
// Function SpeciesSwap
// ***************************************************************************
// Permute Species A with B (safe for species C) Stefano Curtarolo Nov 2008
void xstructure::SpeciesSwap(const uint& specieA,const uint& specieB) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  // some useful checks
  if(num_each_type.size()==0) return; // empty structures have nothing to swap
  if(num_each_type.size()==1) return; // pure structures have nothing to swap
  // check done
  deque<_atom> atoms_buf;
  uint i=0,j=0,k=0;
  uint specieAA=0,specieAA_start=0,specieAA_stop=0;
  uint specieBB=0,specieBB_start=0,specieBB_stop=0;
  // PUT ALPHABETIC
  if(specieA<specieB) {specieAA=specieA;specieBB=specieB;}
  if(specieA>specieB) {specieAA=specieB;specieBB=specieA;}
  // CHECK
  if(specieA==specieB) return;                    // NO SWAP
  if(specieA>num_each_type.size()) return;  // NO SWAP nothing to swap
  if(specieB>num_each_type.size()) return;  // NO SWAP nothing to swap
  // FIND BOUNDARIES
  for(i=0;i<num_each_type.size();i++) {
    if(i<specieAA)  {specieAA_start+=num_each_type.at(i);}
    if(i==specieAA) {specieAA_stop=specieAA_start+num_each_type.at(i);}
    if(i<specieBB)  {specieBB_start+=num_each_type.at(i);}
    if(i==specieBB) {specieBB_stop=specieBB_start+num_each_type.at(i);}
  }
  if(LDEBUG) cerr << "DEBUG  atoms.size()=" << atoms.size() << endl;
  if(LDEBUG) cerr << "DEBUG  specieAA=" << specieAA << endl;
  if(LDEBUG) cerr << "DEBUG  specieBB=" << specieBB << endl;
  if(LDEBUG) cerr << "DEBUG  specieAA_start=" << specieAA_start << endl;
  if(LDEBUG) cerr << "DEBUG  specieAA_stop="  << specieAA_stop << endl;
  if(LDEBUG) cerr << "DEBUG  specieBB_start=" << specieBB_start << endl;
  if(LDEBUG) cerr << "DEBUG  specieBB_stop="  << specieBB_stop << endl;
  // ATOMS  -----------------------------
  atoms_buf.clear();
  for(i=0;i<atoms.size();i++) atoms_buf.push_back(atoms.at(i)); // create buffer
  atoms.clear();   // RECONSTRUCTING
  for(i=0;i<specieAA_start;i++) atoms.push_back(atoms_buf.at(i));                // before specieA preserve
  for(i=specieBB_start;i<specieBB_stop;i++) atoms.push_back(atoms_buf.at(i));    // during specieA swap with specieB
  for(i=specieAA_stop;i<specieBB_start;i++) atoms.push_back(atoms_buf.at(i));    // between specieA and specieB preserve
  for(i=specieAA_start;i<specieAA_stop;i++) atoms.push_back(atoms_buf.at(i));    // during specieB swap with specieA
  for(i=specieBB_stop;i<atoms_buf.size();i++) atoms.push_back(atoms_buf.at(i));  // after specieB preserve
  // swap numbers
  uint iaus;
  iaus=num_each_type.at(specieAA);num_each_type.at(specieAA)=num_each_type.at(specieBB);num_each_type.at(specieBB)=iaus;
  // now fix the types
  k=0;
  for(i=0;i<num_each_type.size();i++)
    for(j=0;j<(uint) num_each_type.at(i);j++)
      atoms.at(k++).type=i; // done
  // swap species
  string saus;double daus;
  if(specieAA<species.size() && specieBB<species.size()) {saus=species.at(specieAA);species.at(specieAA)=species.at(specieBB);species.at(specieBB)=saus;}
  if(specieAA<species_pp.size() && specieBB<species_pp.size()) {saus=species_pp.at(specieAA);species_pp.at(specieAA)=species_pp.at(specieBB);species_pp.at(specieBB)=saus;}
  if(specieAA<comp_each_type.size() && specieBB<comp_each_type.size()) {daus=comp_each_type.at(specieAA);comp_each_type.at(specieAA)=comp_each_type.at(specieBB);comp_each_type.at(specieBB)=daus;} //CO20180705
  if(specieAA<stoich_each_type.size() && specieBB<stoich_each_type.size()) {daus=stoich_each_type.at(specieAA);stoich_each_type.at(specieAA)=stoich_each_type.at(specieBB);stoich_each_type.at(specieBB)=daus;} //CO20180705
  if(specieAA<species_volume.size() && specieBB<species_volume.size()) {daus=species_volume.at(specieAA);species_volume.at(specieAA)=species_volume.at(specieBB);species_volume.at(specieBB)=daus;}
  if(specieAA<species_mass.size() && specieBB<species_mass.size()) {daus=species_mass.at(specieAA);species_mass.at(specieAA)=species_mass.at(specieBB);species_mass.at(specieBB)=daus;}
}

// ***************************************************************************
// Function SpeciesGetAlphabetic
// ***************************************************************************
// Tell if two species are alphabetic!  Stefano Curtarolo Nov 2008
bool xstructure::SpeciesGetAlphabetic(void) {
  stringstream message;
  if(num_each_type.size()!=species.size()) {
    message << "num_each_type.size()!=species.size()   ("<<num_each_type.size()<<","<<species.size()<<")" << endl;
    message << "num_each_type.size()="<<num_each_type.size()<<": "; for(uint i=0;i<num_each_type.size();i++) message << num_each_type.at(i) << " "; message << endl;
    message << "species.size()="<<species.size()<< ": "; for(uint i=0;i<species.size();i++) message << species.at(i) << " "; message << endl;
    message << "species_pp.size()="<<species_pp.size()<< ": "; for(uint i=0;i<species_pp.size();i++) message << species_pp.at(i) << " "; message << endl;
    message << "species_volume.size()="<<species_volume.size()<< ": "; for(uint i=0;i<species_volume.size();i++) message << species_volume.at(i) << " "; message << endl;
    message << "species_mass.size()="<<species_mass.size()<< ": "; for(uint i=0;i<species_mass.size();i++) message << species_mass.at(i) << " "; message << endl;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_RANGE_);
  }
  // some useful checks
  if(species.size()==0) return TRUE; // empty structures are always alphabetic
  if(species.size()==1) return TRUE; // pure structures are always alphabetic
  // check done
  vector<string> sspecies;
  for(uint isp=0;isp<species.size();isp++)
    sspecies.push_back(aurostd::RemoveNumbers(KBIN::VASP_PseudoPotential_CleanName(species.at(isp))));

  for(uint isp=0;isp<sspecies.size()-1;isp++)
    if(sspecies.at(isp)>sspecies.at(isp+1)) return FALSE;
  // otherwise return TRUE;
  return TRUE;
}

// ***************************************************************************
// Function SpeciesPutAlphabetic
// ***************************************************************************
// Tell if two species are alphabetic!  Stefano Curtarolo Nov 2008
bool xstructure::SpeciesPutAlphabetic(void) {
  stringstream message;
  if(num_each_type.size()!=species.size()) {
    message << "num_each_type.size()!=species.size()   ("<<num_each_type.size()<<","<<species.size()<<")";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_RANGE_);
  }
  // some useful checks
  if(species.size()==0) return TRUE; // empty structures are always alphabetic
  if(species.size()==1) return TRUE; // pure structures are always alphabetic
  // check done
  vector<string> sspecies;

  while(SpeciesGetAlphabetic()==FALSE) {
    sspecies.clear();
    for(uint isp=0;isp<species.size();isp++)
      sspecies.push_back(aurostd::RemoveNumbers(KBIN::VASP_PseudoPotential_CleanName(species.at(isp))));
    for(uint isp=0;isp<sspecies.size()-1;isp++)
      if(sspecies.at(isp)>sspecies.at(isp+1)) SpeciesSwap(isp,isp+1);
  }
  MakeBasis(); // repetita iuvant
  // otherwise return TRUE;
  return TRUE;
}

// ***************************************************************************
// Function SpeciesString
// ***************************************************************************
// Returns a string with the species.  Stefano Curtarolo Nov 2008
string xstructure::SpeciesString(void) {
  stringstream strstream;
  strstream.clear();strstream.str(std::string());
  for(uint i=0;i<num_each_type.size();i++) {
    strstream << SpeciesLabel(i);
    if(i!=num_each_type.size()-1)
      strstream << " ";
  }
  return strstream.str();
}

// ***************************************************************************
// Function SetSpecies
// ***************************************************************************
// Set the species  Stefano Curtarolo Nov 2014
// CO20221112 - do not set species explicitly, always use add/remove atoms, it updates EVERYTHING all at once and avoids seg faults
uint xstructure::SetSpecies(const deque<string>& species_in,bool sort_species) {  //CO20230319 - sort species
  bool LDEBUG=(FALSE || XHOST.DEBUG); 
  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " species_in=" << aurostd::joinWDelimiter(species_in,",") << endl;
    cerr << __AFLOW_FUNC__ << " xstr=" << endl << (*this) << endl;
    cerr << __AFLOW_FUNC__ << " xstr.species=" << aurostd::joinWDelimiter(species,",") << endl;
    cerr << __AFLOW_FUNC__ << " xstr.num_each_type=" << aurostd::joinWDelimiter(num_each_type,",") << endl;
  }
  if(species_in.size()!=species.size()){aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"species_in.size()!=species.size()",_VALUE_RANGE_);} //CO20190317
  if(species_in.size()!=num_each_type.size()){aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"species_in.size()!=num_each_type.size()",_VALUE_RANGE_);} //CO20190317
  uint iatom=0;
  deque<_atom> atoms_new;
  for(uint itype=0;itype<num_each_type.size();itype++){
    for(int j=0;j<num_each_type[itype];j++){
      atoms_new.push_back(atoms[iatom++]);
      atoms_new.back().name=species_in[itype];
      atoms_new.back().CleanName();
      //[DX20170921 - Need to keep spin info]atoms_new.back().CleanSpin();
      atoms_new.back().name_is_given=(!atoms_new.back().name.empty());
    }
  }
  ReplaceAtoms(atoms_new,true,sort_species); //CO20230319 - sort species
  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " xstr_updated=" << endl << (*this) << endl;
    cerr << __AFLOW_FUNC__ << " xstr_updated.species=" << aurostd::joinWDelimiter(species,",") << endl;
    cerr << __AFLOW_FUNC__ << " xstr_updated.num_each_type=" << aurostd::joinWDelimiter(num_each_type,",") << endl;
  }
  return species_in.size();
}

// ***************************************************************************
// Function UpdateSpecies() //DX20210202 [from AddAtom, consolidate to function]
// ***************************************************************************
void xstructure::UpdateSpecies(const _atom& atom,bool add_species){return xstructure::UpdateSpecies_20230319(atom,add_species);}  //CO20230319 - fixing species vs. types issue
void xstructure::UpdateSpecies_20230319(const _atom& atom,bool add_species){  //CO20230319 - fixing species vs. types issue

  // Update the species info based on the atom input
  // If the species is already in xstructure, update the number of types
  // and composition of each type, otherwise, add the new species info
  // This code was copied from AddAtom

  bool LDEBUG=(FALSE || XHOST.DEBUG); 
  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    cerr << __AFLOW_FUNC__ << " add_species=" << add_species << endl;
  }

  bool FOUND_SPECIES=FALSE;
  uint species_position=0;

  if(add_species){  //CO20230319 - if we want to decide where to add the atom based on atom.name
    if(atom.name.empty()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"atom.name.empty()",_INPUT_MISSING_);} //DX20210324 - check if name is empty
    for(uint isp=0;isp<species.size()&&FOUND_SPECIES==FALSE;isp++)
      if(atom.name==species[isp]) {FOUND_SPECIES=TRUE;species_position=isp;}
    if(!FOUND_SPECIES){
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " new_species=" << atom.name << endl;
      species_position=species.size();
    }
  }
  else{ //DX20210324 - if name is empty, then use types
    for(uint isp=0;isp<species.size()&&FOUND_SPECIES==FALSE;isp++)
      if(atom.type==(int)isp) {FOUND_SPECIES=TRUE;species_position=isp;}
    if(!FOUND_SPECIES){
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " new_type=" << atom.type << endl;
      species_position=atom.type;
    }
  }

  if(FOUND_SPECIES==FALSE) {
    while(species.size()<=species_position){
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " increasing species.size()=" << species.size() << endl;
      num_each_type.push_back(0);
      comp_each_type.push_back(0.0);
      species.push_back(""); // cerr << "UpdateSpecies=" << atom.name << endl;
      species_pp.push_back(""); // cerr << "UpdateSpecies=" << atom.name << endl;
      species_pp_type.push_back(""); // cerr << "UpdateSpecies=" << atom.name << endl;
      species_pp_version.push_back(""); // cerr << "UpdateSpecies=" << atom.name << endl;
      species_pp_ZVAL.push_back(0.0); // cerr << "UpdateSpecies=" << atom.name << endl;
      species_pp_vLDAU.push_back(deque<double>()); // cerr << "UpdateSpecies=" << atom.name << endl;
      species_volume.push_back(GetAtomVolume("")); // cerr << "UpdateSpecies=" << atom.name << endl;
      species_mass.push_back(GetAtomMass("")); // cerr << "UpdateSpecies=" << atom.name << endl;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " species.size()=" << species.size() << endl;
    }
  }
  //
  if(num_each_type[species_position]==0){
    species[species_position]=atom.name; // cerr << "UpdateSpecies=" << atom.name << endl;
    species_pp[species_position]=atom.name; // cerr << "UpdateSpecies=" << atom.name << endl;
    //species_pp_type[species_position]=""; // cerr << "UpdateSpecies=" << atom.name << endl;
    //species_pp_version[species_position]=""; // cerr << "UpdateSpecies=" << atom.name << endl;
    //species_pp_ZVAL[species_position]=0.0; // cerr << "UpdateSpecies=" << atom.name << endl;
    //species_pp_vLDAU[species_position]=deque<double>(); // cerr << "UpdateSpecies=" << atom.name << endl;
    species_volume[species_position]=GetAtomVolume(atom.name); // cerr << "UpdateSpecies=" << atom.name << endl;
    species_mass[species_position]=GetAtomMass(atom.name); // cerr << "UpdateSpecies=" << atom.name << endl;
  }
  num_each_type[species_position]++;
  comp_each_type[species_position]+=atom.partial_occupation_value;
}
void xstructure::UpdateSpecies_20230101(const _atom& atom,bool add_species){

  // Update the species info based on the atom input
  // If the species is already in xstructure, update the number of types
  // and composition of each type, otherwise, add the new species info
  // This code was copied from AddAtom

  bool LDEBUG=(FALSE || XHOST.DEBUG); 
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

  bool FOUND_SPECIES=FALSE;
  uint species_position=0;

  if(add_species){  //CO20230319 - if we want to decide where to add the atom based on atom.name
    if(atom.name.empty()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"atom.name.empty()",_INPUT_MISSING_);} //DX20210324 - check if name is empty
    for(uint isp=0;isp<species.size()&&FOUND_SPECIES==FALSE;isp++)
      if(atom.name==species[isp]) {FOUND_SPECIES=TRUE;species_position=isp;}
  }
  else{ //DX20210324 - if name is empty, then use types
    for(uint isp=0;isp<species.size()&&FOUND_SPECIES==FALSE;isp++)
      if(atom.type==(int)isp) {FOUND_SPECIES=TRUE;species_position=isp;}
  }

  if(FOUND_SPECIES==FALSE) {
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " new_species=" << atom.name << endl;
    num_each_type.push_back(1);
    comp_each_type.push_back(atom.partial_occupation_value);
    species.push_back(atom.name); // cerr << "UpdateSpecies=" << atom.name << endl;
    species_pp.push_back(atom.name); // cerr << "UpdateSpecies=" << atom.name << endl;
    species_pp_type.push_back(""); // cerr << "UpdateSpecies=" << atom.name << endl;
    species_pp_version.push_back(""); // cerr << "UpdateSpecies=" << atom.name << endl;
    species_pp_ZVAL.push_back(0.0); // cerr << "UpdateSpecies=" << atom.name << endl;
    species_pp_vLDAU.push_back(deque<double>()); // cerr << "UpdateSpecies=" << atom.name << endl;
    species_volume.push_back(GetAtomVolume(atom.name)); // cerr << "UpdateSpecies=" << atom.name << endl;
    species_mass.push_back(GetAtomMass(atom.name)); // cerr << "UpdateSpecies=" << atom.name << endl;
  } else {
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " increasing species_position " << species_position << endl;
    num_each_type[species_position]++;
    comp_each_type[species_position]+=atom.partial_occupation_value;
  }
}

// ***************************************************************************
// NIGGLI NIGGLI NIGGLI NIGGLI NIGGLI NIGGLI NIGGLI NIGGLI NIGGLI NIGGLI NIGGL
// ***************************************************************************
// GetNiggliCell
// ***************************************************************************
// Calculates the reduced Niggli cell.  It is based on
// a program of Eric Wu's.  Here is his + my documentation.
//FUNCTION:
//Subroutine calculates the reduced cell (Niggli cell) for a given
//primitive cell.  It also gives the transformation matrix to go from
//the primitive cell to the Niggli cell.
//
//REFERENCES:
//Y. Lepage, J. Appl. Cryst 15 255-229 (1982) and
//Y. Lepage, J. Appl. Cryst 20 264-269 (1987)
//International tables for crystallography Ch 5,9
//I Krivy, B. Gruber. Acta Cryst 1976 A32 297 (1975)
//W. Clegg., Acta Cryst A37 913-915 (1981)
//A.D. Mighell J. Appl Cryst 9 491-498 (1976)
//
//METHODOLOGY:
//The Niggli cell is UNIQUE and reduced. A Niggli cell is a Buerger cell
//(shortest lattice vectors) but  not necessarily the other way around,
//since the Buerger cell is not unique.  This info is useful because with
//the Niggli cell one knows instantly the type of bravais lattice
//and how to go to the conventional cell if you apply the correct
//algorithm.  There are 2 ways to do this.  The first is the straight
//way -  see International Tables Ch. 9 .  The other involves counting
//the # of 2fold axes and going from there, see the Lepage and Clegg
//references.  This second method is more useful if you want to find
//more than just bravais lattice (symmetry elements or space group
//for example).  The subroutine assumes you input a primitive cell
//and it spits out the Niggli cell and the transformation matrix
//to go from the primitive cell to Niggli cell.
//NB! There is a typo in step 3 of Krivy and Gruber that I corrected.
//[ (ksi*eta*ksi) should be (ksi*eta*zeta) in first part of Eq. 3. ]
//
//INPUT:
//in_lat: a primitive set of lattice vectors,
//in_lat(1,1)=a_x
//in_lat(1,2)=a_y
//in_lat(2,1)=b_x  ...
//
//OUTPUT:
//niggli_lat - reduced (Niggli) cell, same format
//P - Matrix to go from Primitive to Niggli cell
//Q - P^-1;
//
//For more info of the definitions here see Int. Tables Cryst. ch.5
//Let nlat be the niggli lattice, where the first column is a,
//second b, third c.  Let lat be the original cell equivalent.
//The P matrix to transform cells is generally defined by
//nlat = lat * P
//We will follow this convention for defining P and Q.
//Note that is the rest of the code, we define lattices with
//rows equal to a, b, and c.  Therefore, under our conventional
//definitions of nlat, call it nlat_c=nlat^T, and lat, call
//it lat_c=lat^T, we have
//nlat_c = P^T * lat_c
//where ^T means the transpose.

void _sdebug_GetNiggliCell(int,int,double,double,double,double,double,double,ostream&);

bool GetNiggliCell(const xmatrix<double>& in_lat,xmatrix<double>& niggli_lat,xmatrix<double>& P,xmatrix<double>& Q) {
  return GetNiggliCell_20180213(in_lat,niggli_lat,P,Q);
}

//DX20180213 - Niggli algorithm, fixed tolerances and added goto loop in step 6 - START
bool GetNiggliCell_20180213(const xmatrix<double>& in_lat,xmatrix<double>& niggli_lat,xmatrix<double>& P,xmatrix<double>& Q) {
  // return FALSE if failed
  // rerurn TRUE if ok

  //  double RENORM=1E+6;  //DM not used
  int MAXITER=100000;
  double TOL=_ZERO_TOL_; //DX20180212 - from 1e-15 to 1e-10

  // Initialize matrices for tranformations (3x3).
  xmatrix<double> m1(3,3);
  m1(1,2)=1.0; m1(2,1)=1.0; m1(3,3)=-1.0;
  xmatrix<double> m2(3,3);
  m2(1,1)=-1.0; m2(2,3)=1.0; m2(3,2)=1.0;
  xmatrix<double> m3(3,3);
  // Els 1,1 2,2 3,3 are changed later for m3
  xmatrix<double> m4(3,3);
  // Els 1,1 2,2 3,3 are changed later for m4
  xmatrix<double> m5(3,3);
  m5(1,1)=1.0; m5(2,2)=1.0; m5(3,3)=1.0;
  // El 3,2 is changed later for m5
  xmatrix<double> m6(3,3);
  m6(1,1)=1.0; m6(2,2)=1.0; m6(3,3)=1.0;
  // El 3,1 is changed later for m6
  xmatrix<double> m7(3,3);
  m7(1,1)=1.0; m7(2,2)=1.0; m7(3,3)=1.0;
  // El 2,1 is changed later for m7
  xmatrix<double> m8(3,3);
  m8(1,1)=1.0; m8(2,2)=1.0; m8(1,3)=1.0; m8(2,3)=1.0; m8(3,3)=1.0;

  // Initialize a, b, c, ksi, eta, zeta
  xvector<double> indat(6);
  indat=Getabc_angles(in_lat,RADIANS);
  double a=indat(1)*indat(1);
  double b=indat(2)*indat(2);
  double c=indat(3)*indat(3);
  double ksi=2.0*indat(2)*indat(3)*cos(indat(4));
  double eta=2.0*indat(1)*indat(3)*cos(indat(5));
  double zeta=2.0*indat(1)*indat(2)*cos(indat(6));
  // Round to RENORM decimal place to eliminate numerical errors
  //  double a=Nint(indat(1)*indat(1)*RENORM);
  //  double b=Nint(indat(2)*indat(2)*RENORM);
  //  double c=Nint(indat(3)*indat(3)*RENORM);
  //  double ksi=Nint(2.0*indat(2)*indat(3)*cos(indat(4))*RENORM);
  //  double eta=Nint(2.0*indat(1)*indat(3)*cos(indat(5))*RENORM);
  //  double zeta=Nint(2.0*indat(1)*indat(2)*cos(indat(6))*RENORM);

  // Dummy variables
  double temp,temp1,temp2,temp3;

  // Initialize tranformation matrix
  //  xmatrix<double> P(3,3);
  P(1,1)=1; P(2,2)=1; P(3,3)=1;

  int cnt=0;
  // tpx
  // _sdebug_GetNiggliCell(0,0,a,b,c,ksi,eta,zeta,cout);

  // Start loop

LoopHead:

  cnt++;
  //  cerr << "DEBUG: Niggli() cnt=" << cnt << endl;
  if(cnt>MAXITER) {
    //     stringstream oss;
    //     // -------
    //     oss << "EEEEE  aflow " << VERSION << endl;
    //     oss << "EEEEE  ERROR: CellReduceFuncs/GetNiggliCell" << endl;
    //     oss << "EEEEE  ERROR: Number of interations greater the MAXITER = " << MAXITER << endl;
    //     oss << "EEEEE  ERROR: This seems like too many - there is probably some problem." << endl;
    //     oss << "EEEEE  ERROR:." << endl;
    //     //  oss << endl;
    //     // -------
    //     cerr << oss.str();
    //     //  cout << oss.str();
    //     // -------
    return FALSE;
  }


  // Step 1
  //DX20180209 [OBSOLETE] if(((a-b)>TOL) || ((abs(a-b)<TOL) && (abs(ksi)>abs(eta)) ) )
  if(((a-b)>TOL) || ((abs(a-b)<TOL) && ((abs(ksi)-abs(eta))>TOL) ) ) //DX20180209 - more robust; precision
  { //CO20200106 - patching for auto-indenting
    temp=a;a=b;b=temp;
    temp=-ksi;ksi=-eta;eta=temp;
    P=P*m1;
    // tpx
    // _sdebug_GetNiggliCell(1,cnt,a,b,c,ksi,eta,zeta,cout);
  }

  // Step 2
  //DX20180209 [OBSOLETE] if(((b-c)>TOL) || ((abs(b-c)<TOL) && (abs(eta)>abs(zeta)) ) )
  if(((b-c)>TOL) || ((abs(b-c)<TOL) && ((abs(eta)-abs(zeta))>TOL) ) ) //DX20180209 - more robust; precision
  { //CO20200106 - patching for auto-indenting
    temp=b;b=c;c=temp;
    temp=-eta;eta=-zeta;zeta=temp;
    P=P*m2;
    // tpx
    // _sdebug_GetNiggliCell(2,cnt,a,b,c,ksi,eta,zeta,cout);
    goto LoopHead;
  }

  // Step 3
  //DX20180209 [OBSOLETE] if((ksi*eta*zeta)>0.0)
  if((ksi*eta*zeta)>TOL) //DX20180209 - more robust; precision
  { //CO20200106 - patching for auto-indenting
    m3(1,1)=SignNoZero(ksi);
    m3(2,2)=SignNoZero(eta);
    m3(3,3)=SignNoZero(zeta);
    P=P*m3;
    ksi=abs(ksi);
    eta=abs(eta);
    zeta=abs(zeta);
    // tpx
    // _sdebug_GetNiggliCell(3,cnt,a,b,c,ksi,eta,zeta,cout);
  }

  // Step 4
  //EG's original if((abs(ksi)<TOL)||(abs(eta)<TOL)||(abs(zeta)<TOL)||(ksi*eta*zeta<=0))
  if((abs(ksi)<TOL)||(abs(eta)<TOL)||(abs(zeta)<TOL)||(ksi*eta*zeta<=-TOL))  //DX20180212 - if any are zero, then, ksi*eta*zeta is zero
    //if(ksi*eta*zeta<=TOL) //DX20180209 - more robust; precision
  { //CO20200106 - patching for auto-indenting
    m4(1,1)=-(SignNoZero(ksi));
    m4(2,2)=-(SignNoZero(eta));
    m4(3,3)=-(SignNoZero(zeta));
    if(abs(ksi)<TOL) m4(1,1)=m4(2,2)*m4(3,3);
    if(abs(eta)<TOL) m4(2,2)=m4(1,1)*m4(3,3);
    if(abs(zeta)<TOL) m4(3,3)=m4(1,1)*m4(2,2);
    P=P*m4;
    ksi=-abs(ksi);
    eta=-abs(eta);
    zeta=-abs(zeta);
    // tpx
    // _sdebug_GetNiggliCell(4,cnt,a,b,c,ksi,eta,zeta,cout);
  }

  // Step 5
  if(((abs(ksi)-b)>TOL) ||
      ((abs(ksi-b)<TOL) && (2.0*eta-zeta)<-TOL) || //DX20180209 - more robust; precision
      ((abs(ksi+b)<TOL) && (zeta<-TOL)) ) //DX20180209 - more robust; precision
    //DX20180209 [OBSOLETE] if(((abs(ksi)-b)>=TOL) ||
    //DX20180209 [OBSOLETE] ((abs(ksi-b)<TOL) && (2.0*eta<zeta)) ||
    //DX20180209 [OBSOLETE] ((abs(ksi+b)<TOL) && (zeta<0)) )
  { //CO20200106 - patching for auto-indenting
    m5(2,3)=-(SignNoZero(ksi));
    P=P*m5;
    temp1=b+c-ksi*SignNoZero(ksi);
    temp2=eta-zeta*SignNoZero(ksi);
    temp3=ksi-2.0*b*SignNoZero(ksi);
    c=temp1;
    eta=temp2;
    ksi=temp3;
    // tpx
    // _sdebug_GetNiggliCell(5,cnt,a,b,c,ksi,eta,zeta,cout);
    goto LoopHead;
  }

  // Step 6
  if(((abs(eta)-a)>TOL) || 
      ((abs(eta-a)<TOL) && (2.0*ksi-zeta)<-TOL) || //DX20180209 - more robust; precision
      ((abs(eta+a)<TOL) && (zeta<-TOL))) //DX20180209 - more robust; precision
    //DX20180209 [OBSOLETE] ((abs(eta-a)<TOL) && (2.0*ksi<zeta)) ||
    //DX20180209 [OBSOLETE] ((abs(eta+a)<TOL) && (zeta<0)))
  { //CO20200106 - patching for auto-indenting
    m6(1,3)=-SignNoZero(eta);
    P=P*m6;
    temp1=a+c-eta*SignNoZero(eta);
    temp2=ksi-zeta*SignNoZero(eta);
    temp3=eta-2.0*a*SignNoZero(eta);
    c=temp1;
    ksi=temp2;
    eta=temp3;
    // tpx
    // _sdebug_GetNiggliCell(6,cnt,a,b,c,ksi,eta,zeta,cout);
    goto LoopHead; //DX20180212 - this was missing
  }

  // Step 7
  if(((abs(zeta)-a)>TOL) || 
      ((abs(zeta-a)<TOL) && (2.0*ksi-eta)<-TOL) || //DX20180209 - more robust; precision
      ((abs(zeta+a)<TOL) && (eta<-TOL))) //DX20180209 - more robust; precision
    //DX20180209 [OBSOLETE] ((abs(zeta-a)<TOL) && (2.0*ksi<eta)) ||
    //DX20180209 [OBSOLETE] ((abs(zeta+a)<TOL) && (eta<0)))
  { //CO20200106 - patching for auto-indenting
    m7(1,2)=-SignNoZero(zeta);
    P=P*m7;
    temp1=a+b-zeta*SignNoZero(zeta);
    temp2=ksi-eta*SignNoZero(zeta);
    temp3=zeta-2.0*a*SignNoZero(zeta);
    b=temp1;
    ksi=temp2;
    zeta=temp3;
    // tpx
    // _sdebug_GetNiggliCell(7,cnt,a,b,c,ksi,eta,zeta,cout);
    goto LoopHead;
  }

  // Step 8
  //DX20180209 [OBSOLETE] if((ksi+eta+zeta+a+b<0) ||
  //DX20180209 [OBOSLETE]   ((ksi+eta+zeta+a+b<0) && (2.0*(a+eta)+zeta>0)))
  if((ksi+eta+zeta+a+b<-TOL) || //DX20180209 - more robust; precision
      (abs(ksi+eta+zeta+a+b)<TOL && (2.0*(a+eta)+zeta>TOL))) //DX20180209 - more robust; precision
  { //CO20200106 - patching for auto-indenting
    P=P*m8;
    temp1=a+b+c+ksi+eta+zeta;
    temp2=2.0*b+ksi+zeta;
    temp3=2.0*a+eta+zeta;
    c=temp1;
    ksi=temp2;
    eta=temp3;
    // tpx
    // _sdebug_GetNiggliCell(8,cnt,a,b,c,ksi,eta,zeta,cout);
    goto LoopHead;
  }

  // Renormalize back to regular cell (divide by RENORM)
  xvector<double> outdat(6);
  //      outdat(1)=sqrt(a/RENORM);
  //      outdat(2)=sqrt(b/RENORM);
  //      outdat(3)=sqrt(c/RENORM);
  //      outdat(4)=acos(ksi/RENORM/2.0/outdat(2)/outdat(3));
  //      outdat(5)=acos(eta/RENORM/2.0/outdat(1)/outdat(3));
  //      outdat(6)=acos(zeta/RENORM/2.0/outdat(1)/outdat(2));
  outdat(1)=sqrt(a);
  outdat(2)=sqrt(b);
  outdat(3)=sqrt(c);
  outdat(4)=acos(ksi/2.0/outdat(2)/outdat(3));
  outdat(5)=acos(eta/2.0/outdat(1)/outdat(3));
  outdat(6)=acos(zeta/2.0/outdat(1)/outdat(2));

  // Get Niggli cell
  niggli_lat=trasp(P)*in_lat;

  // Get Q
  //  Q = P;
  // xmatrix<double> tmat(3,3);
  //  GaussJordan(Q,tmat); // Returns Q=P^-1.
  Q=inverse(P);

  // Checks
  // Make sure that a,b,c,alpha,beta,gamma are the same from
  // direct calculation and from using P to get niggli_lat.
  xvector<double> poutdat(6);
  poutdat=Getabc_angles(niggli_lat,RADIANS);
  int flag=0;
  for(int i=1;i<=6;i++) {
    if(abs(poutdat(i)-outdat(i))>2*TOL) {flag=1;}
  }
  if(flag) {
    stringstream oss;
    oss << "ERROR: CellReduceFuncs/GetNiggliCell" << endl;
    oss << "ERROR: Lattice parameters/angles as calculated" << endl;
    oss << "ERROR: with Niggli algorithm and P do not match." << endl;
    oss << "ERROR: a,b,c,alpha,beta,gamma from direct algorithm: ";
    oss << outdat << endl;
    oss << "ERROR: a,b,c,alpha,beta,gamma from P matrix: ";
    oss << poutdat << endl;
    oss << "ERROR: Returning." << endl;
    // output
    cerr << oss.str();
    //  cout << oss.str();
    // done
    return FALSE;
  }
  return TRUE; // perfect
}
//DX20180213 - Niggli algorithm, fixed tolerances and added goto loop in step 6 - END

//DX20180213 - Original Niggli routine below
bool GetNiggliCell_20180101(const xmatrix<double>& in_lat,xmatrix<double>& niggli_lat,xmatrix<double>& P,xmatrix<double>& Q) {
  // return FALSE if failed
  // rerurn TRUE if ok

  //  double RENORM=1E+6;  //DM not used
  int MAXITER=100000;
  double TOL=1e-14;

  // Initialize matrices for tranformations (3x3).
  xmatrix<double> m1(3,3);
  m1(1,2)=1.0; m1(2,1)=1.0; m1(3,3)=-1.0;
  xmatrix<double> m2(3,3);
  m2(1,1)=-1.0; m2(2,3)=1.0; m2(3,2)=1.0;
  xmatrix<double> m3(3,3);
  // Els 1,1 2,2 3,3 are changed later for m3
  xmatrix<double> m4(3,3);
  // Els 1,1 2,2 3,3 are changed later for m4
  xmatrix<double> m5(3,3);
  m5(1,1)=1.0; m5(2,2)=1.0; m5(3,3)=1.0;
  // El 3,2 is changed later for m5
  xmatrix<double> m6(3,3);
  m6(1,1)=1.0; m6(2,2)=1.0; m6(3,3)=1.0;
  // El 3,1 is changed later for m6
  xmatrix<double> m7(3,3);
  m7(1,1)=1.0; m7(2,2)=1.0; m7(3,3)=1.0;
  // El 2,1 is changed later for m7
  xmatrix<double> m8(3,3);
  m8(1,1)=1.0; m8(2,2)=1.0; m8(1,3)=1.0; m8(2,3)=1.0; m8(3,3)=1.0;

  // Initialize a, b, c, ksi, eta, zeta
  xvector<double> indat(6);
  indat=Getabc_angles(in_lat,RADIANS);
  double a=indat(1)*indat(1);
  double b=indat(2)*indat(2);
  double c=indat(3)*indat(3);
  double ksi=2.0*indat(2)*indat(3)*cos(indat(4));
  double eta=2.0*indat(1)*indat(3)*cos(indat(5));
  double zeta=2.0*indat(1)*indat(2)*cos(indat(6));
  // Round to RENORM decimal place to eliminate numerical errors
  //  double a=Nint(indat(1)*indat(1)*RENORM);
  //  double b=Nint(indat(2)*indat(2)*RENORM);
  //  double c=Nint(indat(3)*indat(3)*RENORM);
  //  double ksi=Nint(2.0*indat(2)*indat(3)*cos(indat(4))*RENORM);
  //  double eta=Nint(2.0*indat(1)*indat(3)*cos(indat(5))*RENORM);
  //  double zeta=Nint(2.0*indat(1)*indat(2)*cos(indat(6))*RENORM);

  // Dummy variables
  double temp,temp1,temp2,temp3;

  // Initialize tranformation matrix
  //  xmatrix<double> P(3,3);
  P(1,1)=1; P(2,2)=1; P(3,3)=1;

  int cnt=0;
  // tpx
  // _sdebug_GetNiggliCell(0,0,a,b,c,ksi,eta,zeta,cout);

  // Start loop

LoopHead:

  cnt++;
  //  cerr << "DEBUG: Niggli() cnt=" << cnt << endl;
  if(cnt>MAXITER) {
    //     stringstream oss;
    //     // -------
    //     oss << "EEEEE  aflow " << VERSION << endl;
    //     oss << "EEEEE  ERROR: CellReduceFuncs/GetNiggliCell" << endl;
    //     oss << "EEEEE  ERROR: Number of interations greater the MAXITER = " << MAXITER << endl;
    //     oss << "EEEEE  ERROR: This seems like too many - there is probably some problem." << endl;
    //     oss << "EEEEE  ERROR:." << endl;
    //     //  oss << endl;
    //     // -------
    //     cerr << oss.str();
    //     //  cout << oss.str();
    //     // -------
    return FALSE;
  }

  // Step 1
  if(((a-b)>TOL) || ((abs(a-b)<TOL) && (abs(ksi)>abs(eta)) ) ) {
    temp=a;a=b;b=temp;
    temp=-ksi;ksi=-eta;eta=temp;
    P=P*m1;
    // tpx
    // _sdebug_GetNiggliCell(1,cnt,a,b,c,ksi,eta,zeta,cout);
  }

  // Step 2
  if(((b-c)>TOL) || ((abs(b-c)<TOL) && (abs(eta)>abs(zeta)) ) ) {
    temp=b;b=c;c=temp;
    temp=-eta;eta=-zeta;zeta=temp;
    P=P*m2;
    // tpx
    // _sdebug_GetNiggliCell(2,cnt,a,b,c,ksi,eta,zeta,cout);
    goto LoopHead;
  }

  // Step 3
  if((ksi*eta*zeta)>0.0) {
    m3(1,1)=SignNoZero(ksi);
    m3(2,2)=SignNoZero(eta);
    m3(3,3)=SignNoZero(zeta);
    P=P*m3;
    ksi=abs(ksi);
    eta=abs(eta);
    zeta=abs(zeta);
    // tpx
    // _sdebug_GetNiggliCell(3,cnt,a,b,c,ksi,eta,zeta,cout);
  }

  // Step 4
  //EW's original if((abs(ksi)<TOL)||(abs(eta)<TOL)||(abs(zeta)<TOL)||(ksi*eta*zeta<=0))
  if(ksi*eta*zeta<=0)
  { //CO20200106 - patching for auto-indenting
    m4(1,1)=-(SignNoZero(ksi));
    m4(2,2)=-(SignNoZero(eta));
    m4(3,3)=-(SignNoZero(zeta));
    if(abs(ksi)<TOL) m4(1,1)=m4(2,2)*m4(3,3);
    if(abs(eta)<TOL) m4(2,2)=m4(1,1)*m4(3,3);
    if(abs(zeta)<TOL) m4(3,3)=m4(1,1)*m4(2,2);
    P=P*m4;
    ksi=-abs(ksi);
    eta=-abs(eta);
    zeta=-abs(zeta);
    // tpx
    // _sdebug_GetNiggliCell(4,cnt,a,b,c,ksi,eta,zeta,cout);
  }

  // Step 5
  if(((abs(ksi)-b)>=TOL) ||
      ((abs(ksi-b)<TOL) && (2.0*eta<zeta)) ||
      ((abs(ksi+b)<TOL) && (zeta<0)) ) {
    m5(2,3)=-(SignNoZero(ksi));
    P=P*m5;
    temp1=b+c-ksi*SignNoZero(ksi);
    temp2=eta-zeta*SignNoZero(ksi);
    temp3=ksi-2.0*b*SignNoZero(ksi);
    c=temp1;
    eta=temp2;
    ksi=temp3;
    // tpx
    // _sdebug_GetNiggliCell(5,cnt,a,b,c,ksi,eta,zeta,cout);
    goto LoopHead;
  }

  // Step 6
  if(((abs(eta)-a)>TOL) || ((abs(eta-a)<TOL) && (2.0*ksi<zeta)) ||
      ((abs(eta+a)<TOL) && (zeta<0))) {
    m6(1,3)=-SignNoZero(eta);
    P=P*m6;
    temp1=a+c-eta*SignNoZero(eta);
    temp2=ksi-zeta*SignNoZero(eta);
    temp3=eta-2.0*a*SignNoZero(eta);
    c=temp1;
    ksi=temp2;
    eta=temp3;
    // tpx
    // _sdebug_GetNiggliCell(6,cnt,a,b,c,ksi,eta,zeta,cout);
  }

  // Step 7
  if(((abs(zeta)-a)>TOL) || ((abs(zeta-a)<TOL) && (2.0*ksi<eta)) ||
      ((abs(zeta+a)<TOL) && (eta<0))) {
    m7(1,2)=-SignNoZero(zeta);
    P=P*m7;
    temp1=a+b-zeta*SignNoZero(zeta);
    temp2=ksi-eta*SignNoZero(zeta);
    temp3=zeta-2.0*a*SignNoZero(zeta);
    b=temp1;
    ksi=temp2;
    zeta=temp3;
    // tpx
    // _sdebug_GetNiggliCell(7,cnt,a,b,c,ksi,eta,zeta,cout);
    goto LoopHead;
  }

  // Step 8
  if((ksi+eta+zeta+a+b<0) ||
      ((ksi+eta+zeta+a+b<0) && (2.0*(a+eta)+zeta>0))) {
    P=P*m8;
    temp1=a+b+c+ksi+eta+zeta;
    temp2=2.0*b+ksi+zeta;
    temp3=2.0*a+eta+zeta;
    c=temp1;
    ksi=temp2;
    eta=temp3;
    // tpx
    // _sdebug_GetNiggliCell(8,cnt,a,b,c,ksi,eta,zeta,cout);
    goto LoopHead;
  }

  // Renormalize back to regular cell (divide by RENORM)
  xvector<double> outdat(6);
  //      outdat(1)=sqrt(a/RENORM);
  //      outdat(2)=sqrt(b/RENORM);
  //      outdat(3)=sqrt(c/RENORM);
  //      outdat(4)=acos(ksi/RENORM/2.0/outdat(2)/outdat(3));
  //      outdat(5)=acos(eta/RENORM/2.0/outdat(1)/outdat(3));
  //      outdat(6)=acos(zeta/RENORM/2.0/outdat(1)/outdat(2));
  outdat(1)=sqrt(a);
  outdat(2)=sqrt(b);
  outdat(3)=sqrt(c);
  outdat(4)=acos(ksi/2.0/outdat(2)/outdat(3));
  outdat(5)=acos(eta/2.0/outdat(1)/outdat(3));
  outdat(6)=acos(zeta/2.0/outdat(1)/outdat(2));

  // Get Niggli cell
  niggli_lat=trasp(P)*in_lat;

  // Get Q
  //  Q = P;
  // xmatrix<double> tmat(3,3);
  //  GaussJordan(Q,tmat); // Returns Q=P^-1.
  Q=inverse(P);

  // Checks
  // Make sure that a,b,c,alpha,beta,gamma are the same from
  // direct calculation and from using P to get niggli_lat.
  xvector<double> poutdat(6);
  poutdat=Getabc_angles(niggli_lat,RADIANS);
  int flag=0;
  for(int i=1;i<=6;i++) {
    if(abs(poutdat(i)-outdat(i))>2*TOL) {flag=1;}
  }
  if(flag) {
    stringstream oss;
    oss << "ERROR: CellReduceFuncs/GetNiggliCell" << endl;
    oss << "ERROR: Lattice parameters/angles as calculated" << endl;
    oss << "ERROR: with Niggli algorithm and P do not match." << endl;
    oss << "ERROR: a,b,c,alpha,beta,gamma from direct algorithm: ";
    oss << outdat << endl;
    oss << "ERROR: a,b,c,alpha,beta,gamma from P matrix: ";
    oss << poutdat << endl;
    oss << "ERROR: Returning." << endl;
    // output
    cerr << oss.str();
    //  cout << oss.str();
    // done
    return FALSE;
  }
  return TRUE; // perfect
}

xstructure GetNiggliStr(const xstructure& in_str) {
  xmatrix<double> niggli_lat(3,3);
  xmatrix<double> P(3,3);
  xmatrix<double> Q(3,3);
  xstructure sstr=in_str;
  double scale=sstr.scale;
  sstr=ReScale(sstr,1.0);
  bool is_Niggli_ok=TRUE;
  is_Niggli_ok=GetNiggliCell(sstr.lattice,niggli_lat,P,Q);
  if(is_Niggli_ok==FALSE) {
    sstr=in_str;
    sstr.Niggli_has_failed=TRUE;
    return sstr; // got a problem
  } else {
    // Create new str same as before but with Niggli cell params.
    // Transform cell parameters by xyz_niggli = Q * xyz_original
    // (see ITC, ch.5).
    xstructure niggli_str=sstr;
    niggli_str.Niggli_has_failed=FALSE;
    niggli_str.lattice=niggli_lat;
    niggli_str.FixLattices();
    niggli_str.atoms.clear();

    _atom atom;
    for(int ia=0;ia<(int)sstr.atoms.size();ia++) {
      atom=sstr.atoms.at(ia);
      atom.fpos=Q*sstr.atoms.at(ia).fpos;
      atom.cpos=F2C(niggli_str.lattice,atom.fpos);
      niggli_str.atoms.push_back(atom);
    }
    // Set lattice params, atom positions, put in the 0th cell, and reset scale.
    niggli_str=BringInCell(niggli_str);
    for(int ia=0;ia<(int)sstr.atoms.size();ia++)  // NEED TO CLEAR lattice positions of atoms
      clear(sstr.atoms.at(ia).ijk);           // because the niggly reproduces the same structure
    // DONE and fix it back
    niggli_str=ReScale(niggli_str,scale);
    // DONE
    niggli_str.Niggli_calculated=TRUE;
    return niggli_str;
  }
}

xmatrix<double> GetNiggliStr(const xmatrix<double>& lattice) {
  xmatrix<double> niggli_lat(3,3);
  xmatrix<double> P(3,3);
  xmatrix<double> Q(3,3);
  bool is_Niggli_ok=TRUE;
  is_Niggli_ok=GetNiggliCell(lattice,niggli_lat,P,Q);
  if(is_Niggli_ok==FALSE) return lattice;
  return niggli_lat;
}

void _sdebug_GetNiggliCell(int step,int iter,double a,double b,double c,
    double ksi,double eta,double zeta,ostream& sout) {
  if(0) {
    // tpx
    sout << "Step: " << step << endl;
    cout << "Iter: " << iter << endl;
    cout << "a,b,c,ksi,eta,zeta: " <<a<<" "<<b<<" "<<c <<" "<<ksi<<" "<<eta<<" "<<zeta<<endl;
    cout << "sum a+b+c: " << a+b+c << endl;
  }
}


// ***************************************************************************
// Function NiggliUnitCellForm
// ***************************************************************************
// Converts the unit cell to the standardized Niggli form.
// It is based on the GetNiggliStr() function written by Eric Wu and
// Dane Morgan, adapted by Stefano Curtarolo and present in aflow_pflow_funcs.cpp
xstructure NiggliUnitCellForm(const xstructure& a) {
  xstructure b(a);
  b.NiggliUnitCellForm();
  return b;
}

void xstructure::NiggliUnitCellForm(void) {
  if(Niggli_avoid==TRUE) return;
  if(Niggli_calculated==FALSE) {
    *this=GetNiggliStr(*this);
    Niggli_calculated=TRUE;
    FixLattices();
    Niggli_has_failed=FALSE;
  } else {
    if(!XHOST.QUIET) cout << XPID << "00000  MESSAGE NIGGLI Form already Calculated: skipping" << endl;
  }
}

xmatrix<double> NiggliUnitCellForm(const xmatrix<double>& lattice) {
  return GetNiggliStr(lattice);
}

// ***************************************************************************
// Function GetNiggliStructures() //DX20201006
// ***************************************************************************
void GetNiggliStructures(vector<xstructure>& structures, uint start_index, uint end_index){

  // Converts a set of xstructures to their Niggli representation
  // Optional indices can be included; useful for pre-distributed
  // threading schemes
  // Default: run over entire range


  // if end index is greater than structures.size(), then compute Niggli cell for all structures
  if(end_index > structures.size()){ end_index=structures.size(); }

  for(uint i=start_index;i<end_index;i++){
    structures[i].NiggliUnitCellForm();
  }
}

// **************************************************************************
// Function MinkowskiBasisReduction
// **************************************************************************
// This routine takes a set of basis vectors (that form a lattice)
// and reduces them so that they form the shortest possible basis.
// The reduction is performed so that each vector "a_i" is a close
// as possible to the origin while remaining in the affine plane which
// is defined by "a_j", "a_k" but shifted by "a_i", for any choice
// of even permutations of i,j,k in 1,2,3.
// See Lecture notes in computer science, ISSN 0302-974, ANTS - VI :
// algorithmic number theory, 2004, vol. 3076, pp. 338-357
// ISBN 3-540-22156-5
// Written by Gus Hart in F90, recoded by SC in C++ (Sep/08).
xstructure MinkowskiBasisReduction(const xstructure& a) {
  xstructure b(a);
  b.MinkowskiBasisReduction();
  return b;
}

void xstructure::MinkowskiBasisReduction(void) {
  if(LatticeReduction_avoid==TRUE) return;
  if(Minkowski_avoid==TRUE) return;
  if(Minkowski_calculated==FALSE) {
    xmatrix<double> basis(3,3);
    double old_scale=scale;
    ReScale(1.0);
    basis=trasp(lattice);
    basis=aurostd::reduce_to_shortest_basis(basis);
    lattice=trasp(basis);
    FixLattices(); // get f2c c2f and klattice
    for(uint i=0;i<atoms.size();i++)
      atoms[i].fpos=C2F(lattice,atoms[i].cpos);
    // atoms[i]=C2F(lattice,atoms[i]);  // it works, same just curiosity
    // atoms[i]=C2F(atoms[i]);          // it works, same just curiosity
    ReScale(old_scale);
    FixLattices();
    Minkowski_calculated=TRUE;
    Minkowski_has_failed=FALSE;
  } else {
    //    if(!QUIET) cout << XPID << "00000  MESSAGE MINKOWSKI Basis Reduction already Calculated: skipping" << endl;
  }
}

xmatrix<double> MinkowskiBasisReduction(const xmatrix<double>& lattice) {
  xmatrix<double> basis(3,3);
  basis=lattice;
  basis=trasp(basis);
  basis=aurostd::reduce_to_shortest_basis(basis);
  basis=trasp(basis);
  return basis;
}

// ***************************************************************************
// Function GetMinkowskiStructures() //DX20201006
// ***************************************************************************
void GetMinkowskiStructures(vector<xstructure>& structures, uint start_index, uint end_index){

  // Converts a set of xstructures to their Minkowski representation
  // Optional indices can be included; useful for pre-distributed
  // threading schemes
  // Default: run over entire range

  // if end index is greater than structures.size(), then compute Minkowski cell for all structures
  if(end_index > structures.size()){ end_index=structures.size(); }

  for(uint i=start_index;i<end_index;i++){
    structures[i].MinkowskiBasisReduction();
  }
}

// ***************************************************************************
// Function LatticeReduction
// ***************************************************************************
// Lattice Reduction to Max Orthogonality (MINK) and then Niggly Form
xstructure LatticeReduction(const xstructure& a) {
  xstructure b(a);
  b.LatticeReduction();
  return b;
}

void xstructure::LatticeReduction(void) {
  if(LatticeReduction_avoid==TRUE) return;
  if(LatticeReduction_calculated==FALSE) {
    MinkowskiBasisReduction();            // Minkowski first
    NiggliUnitCellForm();                 // Niggli Second
    FixLattices();                        // fix f2c c2f and K-space
    BringInCell();                        // bring atoms in new basis
    LatticeReduction_calculated=TRUE;
    LatticeReduction_has_failed=FALSE;
  } else {
    //   if(!QUIET) cout << XPID << "00000  MESSAGE LATTICE Basis Reduction already Calculated: skipping" << endl;
  }
}

xmatrix<double> LatticeReduction(const xmatrix<double>& lattice) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  xmatrix<double> basis(3,3);
  basis=lattice;
  basis=MinkowskiBasisReduction(basis); // Minkowski first
  basis=NiggliUnitCellForm(basis);      // Niggli Second
  if(LDEBUG) cerr << "WARNING: remeber to FixLattices for f2c,c2f and reciprocal and bring atoms in cell " << endl;
  return basis;
}

// ******************************************************************************
// Fold atoms into cell
// ******************************************************************************
// Folds atoms into the cell
//DX20190214 [OBSOLETE] deque<_atom> foldAtomsInCell(deque<_atom>& atoms, xmatrix<double>& c2f_new, xmatrix<double>& f2c_new, bool skew) { //CO20190520 - removed pointers for bools and doubles, added const where possible
//DX20190214 [OBSOLETE]   double tol = _SYM_TOL_;
//DX20190214 [OBSOLETE]   return foldAtomsInCell(atoms, c2f_new, f2c_new, skew, tol);
//DX20190214 [OBSOLETE]}

// xstructure::foldAtomsInCell()
// modify xstructure in-place
void xstructure::foldAtomsInCell(const xmatrix<double>& lattice_new, bool skew, double tol, bool check_min_dists) { //DX20210104
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " xstr(ORIG)=" << endl << (*this) << endl;}

  deque<_atom> atoms_new = ::foldAtomsInCell((*this), lattice_new, skew, tol, check_min_dists); // fold atoms in //DX20210118 - specify global namespace

  // update xstructure info
  (*this).lattice=lattice_new;
  // sort and update atom counts/order/types/basis/etc.
  std::stable_sort(atoms_new.begin(),atoms_new.end(),sortAtomsNames); //DX20210129
  (*this).ReplaceAtoms(atoms_new);
  (*this).BringInCell();
  (*this).FixLattices();
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " xstr(NEW)=" << endl << (*this) << endl;}
}

deque<_atom> foldAtomsInCell(const xstructure& a,const xmatrix<double>& lattice_new, bool skew, double tol, bool check_min_dists) { //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 = added check_min_dists bool
  bool LDEBUG=(FALSE || XHOST.DEBUG);

  double volume_original=abs(aurostd::det(a.lattice));
  double volume_new=abs(aurostd::det(lattice_new));
  bool fold_in_only=( (volume_new < volume_original) || (aurostd::isequal(volume_original,volume_new)) );

  deque<_atom> atoms_orig=a.atoms;  //need to make a copy for the pointer
  deque<_atom>* ptr_atoms=&atoms_orig;
  deque<_atom> atoms=*ptr_atoms; //DX20210129 - this need to be done before if-statement since we put atomic_grid inside if-statement
  if(!fold_in_only){
    xstructure atomic_grid; //stays empty if not needed //DX+ME20210111 - added inside if-statement
    double radius=RadiusSphereLattice(lattice_new);
    xvector<int> dims=LatticeDimensionSphere(a.lattice,radius);//int dim=max(dims)+1; //dim=3;  //CO20190520
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " a.lattice=" << endl;cerr << a.lattice << endl;
      cerr << __AFLOW_FUNC__ << " lattice_new=" << endl;cerr << lattice_new << endl;
      cerr << __AFLOW_FUNC__ << " vol(a.lattice)=" << abs(aurostd::det(a.lattice)) << endl;
      cerr << __AFLOW_FUNC__ << " vol(lattice_new)=" << abs(aurostd::det(lattice_new)) << endl;
      cerr << __AFLOW_FUNC__ << " radius(a.lattice)=" << RadiusSphereLattice(a.lattice) << endl;
      cerr << __AFLOW_FUNC__ << " radius(lattice_new)=" << radius << endl;
      cerr << __AFLOW_FUNC__ << " dims=" << dims << endl;
    }
    //[CO20190520 - excessive, too large of an exploration radius]xmatrix<double> supercell; supercell(1,1)=dims(1); supercell(2,2)=dims(2); supercell(3,3)=dims(3);  //NO NEED, function ensures radius is encompassed //be safe and go +1 out
    //xmatrix<double> supercell; supercell(1,1)=dim; supercell(2,2)=dim; supercell(3,3)=dim;  //NO NEED, function ensures radius is encompassed //be safe and go +1 out
    //vector<int> sc2pcMap, pc2scMap; //dummy
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " building atomic grid with dims=[" << dims << "]" << endl;}
    atomic_grid=a;atomic_grid.clean(); //DX20191220 - uppercase to lowercase clean
    atomic_grid.GenerateGridAtoms(dims[1],dims[2],dims[3]); //much faster than supercell
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " atomic grid built" << endl;}
    ptr_atoms=&atomic_grid.grid_atoms;  //CO20190808 - GenerateGridAtoms() populates grid_atoms, not atoms
    atoms=*ptr_atoms; //DX20210129 - set inside if-statement, otherwise, grid atoms goes out of scope
  }
  //DX20210129 [OBOSLETE] const deque<_atom> atoms=*ptr_atoms;

  return foldAtomsInCell(atoms,a.lattice,lattice_new,skew,tol,check_min_dists); //DX20190619 = added check_min_dists bool
}

deque<_atom> foldAtomsInCell(const deque<_atom>& atoms,const xmatrix<double>& lattice_orig,const xmatrix<double>& lattice_new,bool skew, double tol, bool check_min_dists) {  //DX20190619 - added check_min_dists bool
  bool LDEBUG=(FALSE || XHOST.DEBUG);

  deque<_atom> atoms_in_cell;

  xmatrix<double> f2c_new=trasp(lattice_new);
  xmatrix<double> c2f_new=inverse(f2c_new);

  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " f2c_new=" << endl;cerr << f2c_new << endl;
    cerr << __AFLOW_FUNC__ << " c2f_new=" << endl;cerr << c2f_new << endl;
  }

  _atom tmp;
  for (uint j = 0; j < atoms.size(); j++) {
    //[CO20190520 - this case is not needed]if(atoms_in_cell.size() == 0) {
    //[CO20190520 - this case is not needed]  //[OBSOLETE]atoms[j].fpos = c2f_new * atoms[j].cpos;
    //[CO20190520 - this case is not needed]  //[OBSOLETE]atoms[j].fpos = BringInCell(atoms[j].fpos);
    //[CO20190520 - this case is not needed]  //[OBSOLETE]atoms[j].cpos = f2c_new * atoms[j].fpos;
    //[CO20190520 - this case is not needed]  atoms_in_cell.push_back(atoms[j]);
    //[CO20190520 - this case is not needed]  atoms_in_cell.back().fpos = BringInCell(c2f_new * atoms[j].cpos);
    //[CO20190520 - this case is not needed]  atoms_in_cell.back().cpos = f2c_new * atoms_in_cell.back().fpos;
    //[CO20190520 - this case is not needed]  atoms_in_cell.back().ijk(1)=0; atoms_in_cell.back().ijk(2)=0; atoms_in_cell.back().ijk(3)=0;
    //[CO20190520 - this case is not needed]} else {  //[CO20200106 - close bracket for indenting]}
    //bool duplicate_atom = false;
    tmp.fpos = BringInCell(c2f_new * atoms[j].cpos);
    tmp.cpos = f2c_new * tmp.fpos;
    //[OBSOLETE]for (uint a = 0; a < atoms_in_cell.size(); a++) {
    //[OBSOLETE]  if(MapAtomsInNewCell(atoms_in_cell[a], tmp, c2f_orig, f2c_new, skew, tol))
    //[OBSOLETE]  if(MapAtoms(atoms_in_cell[a], tmp, c2f_orig, f2c_new, skew, tol))
    //[OBSOLETE]  { //CO20200106 - patching for auto-indenting
    //[OBSOLETE]    duplicate_atom = true;
    //[OBSOLETE]    break;
    //[OBSOLETE]  }
    //[OBSOLETE]}
    //[OBSOLETE]if(duplicate_atom == false) {
    if(!SYM::MapAtom(atoms_in_cell,tmp,false,lattice_new,f2c_new,skew,tol)){ //DX20190619 - lattice_new and f2c_new as input
      //[OBSOLETE]atoms[j].fpos = tmp.fpos; //BringInCell(tmp.fpos);
      //[OBSOLETE]atoms[j].cpos = tmp.cpos; //f2c_new * atoms[j].fpos;
      atoms_in_cell.push_back(atoms[j]);
      atoms_in_cell.back().fpos = tmp.fpos;
      atoms_in_cell.back().cpos = tmp.cpos;
      atoms_in_cell.back().ijk(1)=0; atoms_in_cell.back().ijk(2)=0; atoms_in_cell.back().ijk(3)=0;
    }
    //[CO20190520 - this case is not needed]}
  }

  if(check_min_dists){ //DX20190613
    double min_dist_orig=SYM::minimumDistance(atoms);  //lattice_orig //this does NOT work if we use GenerateGridAtoms (no longer periodic with lattice), so simply compare distances between atoms. NOTE: this is no longer the TRUE minimumDistance(), which requires knowledge of the lattice vectors
    double min_dist_new=SYM::minimumDistance(atoms_in_cell);  //lattice_new
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " lattice_orig=" << endl;cerr << lattice_orig << endl;
      cerr << __AFLOW_FUNC__ << " lattice_new=" << endl;cerr << lattice_new << endl;
      cerr << __AFLOW_FUNC__ << " atoms_orig=" << endl;for(uint i=0;i<atoms.size();i++){cerr << atoms[i] << endl;}
      cerr << __AFLOW_FUNC__ << " min_dist_orig=" << endl;cerr << min_dist_orig << endl;
      cerr << __AFLOW_FUNC__ << " min_dist_new=" << endl;cerr << min_dist_new << endl;
    }
    if(!aurostd::isequal(min_dist_orig,min_dist_new,0.1)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Minimum distance changed, check that atoms are not rotated",_RUNTIME_ERROR_);}
  }

  // sort atoms //DX+CO20210119
  std::stable_sort(atoms_in_cell.begin(),atoms_in_cell.end(),sortAtomsNames);

  return atoms_in_cell;
}

// ***************************************************************************
// Function GetPrimitiveVASP()
// ***************************************************************************
// get reduced lattice in reciprocal space, converts to real and folds atoms into this lattice
// this should be the fastest lattice for VASP
//CO20180409 - refer to standard primitive instead, this function is simply a test, probably not optimal to standard primitive
xstructure GetPrimitiveVASP(const xstructure& a) {
  double tol=a.sym_eps;
  if(tol==AUROSTD_NAN){tol=SYM::defaultTolerance(a);}
  return GetPrimitiveVASP(a,tol);
}

xstructure GetPrimitiveVASP(const xstructure& a,double tol) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  xstructure b(a);
  b.ReScale(1.0);
  xmatrix<double> klattice=ReciprocalLattice(b.lattice,b.scale);        //get reciprocal lattice
  klattice=MinkowskiBasisReduction(klattice);                           //get min of lattice
  xmatrix<double> lattice_new=ReciprocalLattice(klattice,b.scale);      //get direct lattice
  //get skew - START
  if(b.dist_nn_min==AUROSTD_NAN){b.dist_nn_min=SYM::minimumDistance(b.atoms,b.lattice);}
  bool skew=SYM::isLatticeSkewed(lattice_new,b.dist_nn_min,tol);
  //get skew - STOP
  b.lattice=lattice_new; //b.f2c=trasp(b.lattice); b.c2f=inverse(b.f2c);  //set new lattice, f2c, c2f
  //b.atoms=foldAtomsInCell(b.atoms,b.c2f,b.f2c,skew,tol);                //fold atoms into new lattice
  //[CO20190520 - ReplaceAtoms()]b.atoms=foldAtomsInCell(a,b.lattice,skew,tol,true);                   //fold atoms into new lattice, don't bother folding out (slow)
  deque<_atom> atoms=foldAtomsInCell(a,b.lattice,skew,tol);  //CO20190520 - ReplaceAtoms()
  b.ReplaceAtoms(atoms);  //CO20190520 - ReplaceAtoms()
  if((0||LDEBUG)&&!isequal(a.lattice,b.lattice,1e-3)){
    cerr << "-----------------------------------------------------------------------" << endl;
    cerr << "ORIG STRUCTURE " << endl;
    cerr << a;
    b.ReScale(1.0); b.ShiftOriginToAtom(0); b.BringInCell(); //fast clean for comparison
    cerr << "NEW STRUCTURE " << endl;
    cerr << b;
    cerr << "STRUCTURES IDENTICAL = " << (compare::structuresMatch(a,b,true) ? "TRUE" : "FALSE") << endl;
    cerr << "-----------------------------------------------------------------------" << endl;
  }
  return b;
}

// **************************************************************************
// Function BringInCell
// **************************************************************************
// these routines take atoms or structures and brings them
// to the unit cell (incell). There is some overloading for
// structures. SC Aug2007
// EDITED BY CO+DX to include a tolerance that converts from Cartesian
// space to direct space via covariant and contravariant transformations.
#define _EPS_roundoff_ 0.001
#define _incellcutoff_ (1.0-_EPS_roundoff_)
//DX+CO, IF EPSILON IS PROVIDED, THEN WE ASSUME YOU WANT US TO MOVE THE ATOMS, OTHERWISE IT'S A HARD CUT OFF

// BringInCell() - ROBUST (DX+ME+CO20190905)
// this function brings components/xvectors/_atoms into a unit cell
// the upper bound and lower bound of the cell can be adjusted
// (e.g., standard unit cell : 0.0 to 1.0 ; unit cell centered on origin: -0.5 to 0.5)
// the tolerance can be tuned and "shifts" the bounds, favoring the lower bound
// (e.g., for bounds 0.0 to 1.0, we favor the origin, bringing values 
// inside the cell if they are between "lower_bound-tolerance" and "upper_bound-tolerance")
// the AFLOW developers suggest a hard cutoff, e.g., _ZERO_TOL_ (DX+ME+CO)

// **************************************************************************
// BringInCellInPlace() (change value/object in place)

// -------------------------------------------------------------------
// double (change in place)
void BringInCellInPlace(double& component, double tolerance, double upper_bound, double lower_bound) {
  if (component == INFINITY || component != component || component == -INFINITY) {
    stringstream message; // Moving the stringstream outside the if-statement would add a lot to the run time (~1 sec). 
    message << "Value of component is invalid: (+-) INF or NAN value (component=" << component << ").";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ERROR_); //DX20190905 - replaced cerr with throw
  }
  if (std::signbit(tolerance)) { //DX20191115 
    stringstream message; // Moving the stringstream outside the if-statement would add a lot to the run time (~1 sec). 
    message << "Sign of tolerance is negative (tolerance=" << tolerance << ").";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_);
  }
  while (component - upper_bound >= -tolerance){ component -= 1.0; } //note: non-symmetric, favors values closer to lower bound
  while (component - lower_bound < -tolerance){ component += 1.0; }
}

// -------------------------------------------------------------------
// xvector (change in place)
void BringInCellInPlace(xvector<double>& fpos, double tolerance, double upper_bound, double lower_bound) {
  for (int i = fpos.lrows; i <= fpos.urows; i++) {
    BringInCellInPlace(fpos[i], tolerance, upper_bound, lower_bound);
  }
}

// -------------------------------------------------------------------
// _atom (change in place, updates both fpos and pos) 
void BringInCellInPlace(_atom& atom, const xmatrix<double>& lattice, double tolerance, double upper_bound, double lower_bound) { //DX20190904
  BringInCellInPlaceFPOS(atom, tolerance, upper_bound, lower_bound); // update fpos first

  // update cpos
  atom.cpos=F2C(lattice,atom.fpos); // update cpos next
  atom.isincell=TRUE;
}

// -------------------------------------------------------------------
// xstructure (change in place) 
void BringInCellInPlace(xstructure& xstr, double tolerance, double upper_bound, double lower_bound) { //DX20190904
  for(uint i=0;i<xstr.atoms.size();i++){
    BringInCellInPlace(xstr.atoms[i], xstr.lattice, tolerance, upper_bound, lower_bound);
  }
}

// **************************************************************************
// BringInCell() (return new value/object)

// -------------------------------------------------------------------
// double (return new double)
double BringInCell(double component_in, double tolerance, double upper_bound, double lower_bound) {
  double component_out = component_in;
  if (component_out == INFINITY || component_out != component_out || component_out == -INFINITY) {
    stringstream message; // Moving the stringstream outside the if-statement would add a lot to the run time (~1 sec). 
    message << "Value of component is invalid: (+-) INF or NAN value (component=" << component_out << ").";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ERROR_); //DX20190905 - replaced cerr with throw
  }
  if (std::signbit(tolerance)) { //DX20191115 
    stringstream message; // Moving the stringstream outside the if-statement would add a lot to the run time (~1 sec). 
    message << "Sign of tolerance is negative (tolerance=" << tolerance << ").";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_);
  }
  while (component_out - upper_bound >= -tolerance) { component_out -= 1.0; } //note: non-symmetric, favors values closer to lower bound
  while (component_out - lower_bound < -tolerance) { component_out += 1.0; }
  return component_out;
}

// -------------------------------------------------------------------
// xvector (return new xvector)
xvector<double> BringInCell(const xvector<double>& fpos_in, double tolerance, double upper_bound, double lower_bound) { //DX20190904
  xvector<double> fpos_out = fpos_in;
  for (int i = fpos_out.lrows; i <= fpos_out.urows; i++) {
    BringInCellInPlace(fpos_out[i], tolerance, upper_bound, lower_bound);
  }
  return fpos_out;
}

// -------------------------------------------------------------------
// _atom (return new _atom, update fpos only) 
_atom BringInCellFPOS(const _atom& atom_in, double tolerance, double upper_bound, double lower_bound) { //DX20190904
  _atom atom_out = atom_in;
  BringInCellInPlace(atom_out.fpos, tolerance, upper_bound, lower_bound);

  // update ijk
  for(int i=atom_out.fpos.lrows; i<=atom_out.fpos.urows; i++){
    atom_out.ijk(i) = (int)(atom_out.fpos[i]-atom_in.fpos[i]);
  }
  return atom_out;
}


// -------------------------------------------------------------------
// xstructure (return xstructure) 
xstructure BringInCell(const xstructure& xstr_in, double tolerance, double upper_bound, double lower_bound) { //DX20190904
  xstructure xstr_out = xstr_in;
  for(uint i=0;i<xstr_out.atoms.size();i++){
    BringInCellInPlace(xstr_out.atoms[i], xstr_out.lattice, tolerance, upper_bound, lower_bound);
  }
  return xstr_out;
}

// **************************************************************************
// BringInCellFPOS() (updates atom.fpos only)

// -------------------------------------------------------------------
// _atom (change in place, update fpos only) 
void BringInCellInPlaceFPOS(_atom& atom, double tolerance, double upper_bound, double lower_bound) { //DX20190904
  xvector<double> orig_fpos = atom.fpos; //DX - needed for ijk later
  BringInCellInPlace(atom.fpos, tolerance, upper_bound, lower_bound);

  // update ijk
  for(int i=atom.fpos.lrows; i<=atom.fpos.urows; i++){
    atom.ijk(i) = (int)(atom.fpos[i]-orig_fpos[i]);
  }
}

// -------------------------------------------------------------------
// _atom (return new _atom, update fpos and cpos) 
_atom BringInCell(const _atom& atom_in, const xmatrix<double>& lattice, double tolerance, double upper_bound, double lower_bound) { //DX20190904
  _atom atom_out = BringInCellFPOS(atom_in, tolerance, upper_bound, lower_bound);

  // update cpos
  atom_out.cpos=F2C(lattice,atom_out.fpos);
  atom_out.isincell=TRUE;

  return atom_out;
}

// **************************************************************************
// BringInCell() (method for xstructure)

// -------------------------------------------------------------------
// xstructure (xstructure method) 
void xstructure::BringInCell(double tolerance, double upper_bound, double lower_bound) { //DX20190904
  for(uint i=0;i<atoms.size();i++){
    BringInCellInPlace(atoms[i], lattice, tolerance, upper_bound, lower_bound); //DX "::" to access outside of xstructure class
  }
}

//DX20190905 [OBSOLETE] //CO20190114 - DO NOT USE OVERLOADS OF BRINGINCELL() WITH EPSILON UNLESS YOU KNOW WHAT YOU ARE DOING
//DX20190905 [OBSOLETE] //DEFAULT TO OVERLOADS WITHOUT EPSILON WHICH USE HARD CUTOFF OF _ZERO_TOL_ = 1e-10
//DX20190905 [OBSOLETE] xvector<double> BringInCell(const xvector<double>& v_in,double epsilon) {
//DX20190905 [OBSOLETE]   //DX+CO START
//DX20190905 [OBSOLETE]   return BringInCell_20160101(v_in,epsilon);
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] //CO20190114 - DO NOT USE OVERLOADS OF BRINGINCELL() WITH EPSILON UNLESS YOU KNOW WHAT YOU ARE DOING
//DX20190905 [OBSOLETE] //DEFAULT TO OVERLOADS WITHOUT EPSILON WHICH USE HARD CUTOFF OF _ZERO_TOL_ = 1e-10
//DX20190905 [OBSOLETE] xvector<double> BringInCell_20161115(const xvector<double>& v_in,double epsilon) {
//DX20190905 [OBSOLETE]   return BringInCell_20160101(v_in,epsilon); //SYM::mod_one_xvec(v_in);
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] //CO20190114 - DO NOT USE OVERLOADS OF BRINGINCELL() WITH EPSILON UNLESS YOU KNOW WHAT YOU ARE DOING
//DX20190905 [OBSOLETE] //DEFAULT TO OVERLOADS WITHOUT EPSILON WHICH USE HARD CUTOFF OF _ZERO_TOL_ = 1e-10
//DX20190905 [OBSOLETE] xvector<double> BringInCell_20160101(const xvector<double>& v_in,double epsilon) {
//DX20190905 [OBSOLETE]   double incelleps=1.0-epsilon;
//DX20190905 [OBSOLETE]   xvector<double> v_out(v_in.urows,v_in.lrows);
//DX20190905 [OBSOLETE]   for(int i=v_out.lrows;i<=v_out.urows;i++) {
//DX20190905 [OBSOLETE]     v_out(i)=v_in(i);
//DX20190905 [OBSOLETE]     while(v_out(i)> incelleps) v_out(i)-=1.0;
//DX20190905 [OBSOLETE]     while(v_out(i)< 0.0)       v_out(i)+=1.0;
//DX20190905 [OBSOLETE]     if(abs(v_out(i))<epsilon)  v_out(i)=0.0;
//DX20190905 [OBSOLETE]     if(v_out(i)> incelleps)    v_out(i)=0.0;
//DX20190905 [OBSOLETE]   }
//DX20190905 [OBSOLETE]   //v_out=roundoff(v_out,_EPS_sym_);
//DX20190905 [OBSOLETE]   //DX+CO END
//DX20190905 [OBSOLETE]   return v_out;
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] xvector<double> BringInCell(const xvector<double>& v_in) {
//DX20190905 [OBSOLETE]   //DX+CO START
//DX20190905 [OBSOLETE]   return SYM::mod_one_xvec(v_in); //hard cutoff, _ZERO_TOL_ = 1e-10
//DX20190905 [OBSOLETE]   //_EPS_sym_ is the universal tolerance here
//DX20190905 [OBSOLETE]   //return BringInCell(v_in,_EPS_sym_);
//DX20190905 [OBSOLETE]   //DX+CO END
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] //DX+CO START
//DX20190905 [OBSOLETE] //CO20190114 - DO NOT USE OVERLOADS OF BRINGINCELL() WITH EPSILON UNLESS YOU KNOW WHAT YOU ARE DOING
//DX20190905 [OBSOLETE] //DEFAULT TO OVERLOADS WITHOUT EPSILON WHICH USE HARD CUTOFF OF _ZERO_TOL_ = 1e-10
//DX20190905 [OBSOLETE] _atom BringInCell(const _atom& atom_in,const xmatrix<double>& lattice,double epsilon) {
//DX20190905 [OBSOLETE]  
//DX20190905 [OBSOLETE]   return BringInCell_20160101(atom_in,lattice,epsilon);
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] //CO20190114 - DO NOT USE OVERLOADS OF BRINGINCELL() WITH EPSILON UNLESS YOU KNOW WHAT YOU ARE DOING
//DX20190905 [OBSOLETE] //DEFAULT TO OVERLOADS WITHOUT EPSILON WHICH USE HARD CUTOFF OF _ZERO_TOL_ = 1e-10
//DX20190905 [OBSOLETE] _atom BringInCell_20161115(const _atom& atom_in,const xmatrix<double>& lattice,double epsilon) {
//DX20190905 [OBSOLETE]   return BringInCell_20160101(atom_in,lattice,epsilon);//BringInCell(atom_in,lattice);
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] //DX+CO END
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] //CO20190114 - DO NOT USE OVERLOADS OF BRINGINCELL() WITH EPSILON UNLESS YOU KNOW WHAT YOU ARE DOING
//DX20190905 [OBSOLETE] //DEFAULT TO OVERLOADS WITHOUT EPSILON WHICH USE HARD CUTOFF OF _ZERO_TOL_ = 1e-10
//DX20190905 [OBSOLETE] _atom BringInCell_20160101(const _atom& atom_in,const xmatrix<double>& lattice,double epsilon) {
//DX20190905 [OBSOLETE]   double _incelleps=1.0;
//DX20190905 [OBSOLETE]   if(epsilon>0.0) _incelleps=1.0-epsilon;else _incelleps=1.0;
//DX20190905 [OBSOLETE]   _atom atom;
//DX20190905 [OBSOLETE]   atom=atom_in;
//DX20190905 [OBSOLETE]   atom.ijk=atom_in.ijk;
//DX20190905 [OBSOLETE]   for(int i=1;i<=3;i++) {
//DX20190905 [OBSOLETE]     while(atom.fpos(i)> _incelleps) {
//DX20190905 [OBSOLETE]       atom.fpos(i)-=1.0;
//DX20190905 [OBSOLETE]       atom.ijk(i)++;
//DX20190905 [OBSOLETE]     }
//DX20190905 [OBSOLETE]     while(atom.fpos(i)<0.0) {
//DX20190905 [OBSOLETE]       atom.fpos(i)+=1.0;
//DX20190905 [OBSOLETE]       atom.ijk(i)--;
//DX20190905 [OBSOLETE]     }
//DX20190905 [OBSOLETE]     if(epsilon>0.0) { // roundoff only if epsilon>0.0
//DX20190905 [OBSOLETE]       if(abs(atom.fpos(i))<epsilon) atom.fpos(i)=0.0;
//DX20190905 [OBSOLETE]       if(atom.fpos(i)>_incelleps) atom.fpos(i)=0.0;
//DX20190905 [OBSOLETE]     }
//DX20190905 [OBSOLETE]   }
//DX20190905 [OBSOLETE]   atom.cpos=F2C(lattice,atom.fpos);
//DX20190905 [OBSOLETE]   atom.isincell=TRUE;
//DX20190905 [OBSOLETE]   return atom;
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] //DX+CO START
//DX20190905 [OBSOLETE] _atom BringInCell(const _atom& atom_in,const xmatrix<double>& lattice) {
//DX20190905 [OBSOLETE]   return BringInCell_20161115(atom_in,lattice);
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] _atom BringInCell_20161115(const _atom& atom_in,const xmatrix<double>& lattice) {
//DX20190905 [OBSOLETE]   //_atom atom;
//DX20190905 [OBSOLETE]   //atom=atom_in;
//DX20190905 [OBSOLETE]   //atom.ijk=atom_in.ijk; //just to be sure //it's part of the assignment operator now
//DX20190905 [OBSOLETE]   _atom atom=SYM::mod_one_atom(atom_in);
//DX20190905 [OBSOLETE]   atom.cpos=F2C(lattice,atom.fpos);
//DX20190905 [OBSOLETE]   atom.isincell=TRUE;
//DX20190905 [OBSOLETE]   return atom;
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] //DX+CO END
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] _atom BringInCell_20160101(const _atom& atom_in,const xmatrix<double>& lattice) {
//DX20190905 [OBSOLETE]   return BringInCell(atom_in,lattice,_EPS_roundoff_);
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] //_atom _BringInCell(const _atom& atom_in,const xmatrix<double>& lattice) {
//DX20190905 [OBSOLETE] //  return BringInCell(atom_in,lattice);
//DX20190905 [OBSOLETE] //}
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] //DX+CO START
//DX20190905 [OBSOLETE] //CO20190114 - DO NOT USE OVERLOADS OF BRINGINCELL() WITH EPSILON UNLESS YOU KNOW WHAT YOU ARE DOING
//DX20190905 [OBSOLETE] //DEFAULT TO OVERLOADS WITHOUT EPSILON WHICH USE HARD CUTOFF OF _ZERO_TOL_ = 1e-10
//DX20190905 [OBSOLETE] xstructure BringInCell(const xstructure& a,double epsilon) {
//DX20190905 [OBSOLETE]   return BringInCell_20160101(a,epsilon);
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] //CO20190114 - DO NOT USE OVERLOADS OF BRINGINCELL() WITH EPSILON UNLESS YOU KNOW WHAT YOU ARE DOING
//DX20190905 [OBSOLETE] //DEFAULT TO OVERLOADS WITHOUT EPSILON WHICH USE HARD CUTOFF OF _ZERO_TOL_ = 1e-10
//DX20190905 [OBSOLETE] xstructure BringInCell_20161115(const xstructure& a,double epsilon) {
//DX20190905 [OBSOLETE]   xstructure b(a); // copies everything
//DX20190905 [OBSOLETE]   for(int i=0;i<(int)a.atoms.size();i++){
//DX20190905 [OBSOLETE]     b.atoms[i]=BringInCell(a.atoms[i],a.lattice,epsilon);//BringInCell(a.atoms[i],a.lattice);
//DX20190905 [OBSOLETE]   }
//DX20190905 [OBSOLETE]   return b;
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] //CO20190114 - DO NOT USE OVERLOADS OF BRINGINCELL() WITH EPSILON UNLESS YOU KNOW WHAT YOU ARE DOING
//DX20190905 [OBSOLETE] //DEFAULT TO OVERLOADS WITHOUT EPSILON WHICH USE HARD CUTOFF OF _ZERO_TOL_ = 1e-10
//DX20190905 [OBSOLETE] xstructure BringInCell_20160101(const xstructure& a,double epsilon) {
//DX20190905 [OBSOLETE]   xstructure b(a); // copies everything
//DX20190905 [OBSOLETE]   for(int i=0;i<(int)a.atoms.size();i++) {
//DX20190905 [OBSOLETE]     b.atoms[i]=BringInCell(a.atoms[i],a.lattice,epsilon);//BringInCell(a.atoms[i],a.lattice);
//DX20190905 [OBSOLETE]   }
//DX20190905 [OBSOLETE]   return b;
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] xstructure BringInCell(const xstructure& a) {
//DX20190905 [OBSOLETE]   //DX+CO START
//DX20190905 [OBSOLETE]   return BringInCell_20161115(a);
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] //SUPER SLOW, AVOID!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//DX20190905 [OBSOLETE] //USE THE ONE THAT ONLY COPIES ATOMS, CO
//DX20190905 [OBSOLETE] //a.sym_eps is STRONGLY PREFERRED
//DX20190905 [OBSOLETE] //however, it requires knowledge (and calculating of the full symmetry, i.e. slow factor group)
//DX20190905 [OBSOLETE] //therefore, if you call the default without having a.sym_eps calculated, the assumption
//DX20190905 [OBSOLETE] //is that a constant tolerance is good enough (speed)
//DX20190905 [OBSOLETE] //if you want to do this correctly, calculate the full symmetry of the structure first (CalculateFullSymmetry() in
//DX20190905 [OBSOLETE] //aflow_aconvasp_main.cpp)
//DX20190905 [OBSOLETE] xstructure BringInCell_20161115(const xstructure& a) {
//DX20190905 [OBSOLETE]   xstructure b=a;
//DX20190905 [OBSOLETE]   //  double _eps_;
//DX20190905 [OBSOLETE]   //  if(a.sym_eps!=AUROSTD_NAN){ //Tolerance came from user or was calculated
//DX20190905 [OBSOLETE]   //    _eps_=a.sym_eps;
//DX20190905 [OBSOLETE]   //  }
//DX20190905 [OBSOLETE]   //  else {
//DX20190905 [OBSOLETE]   //    _eps_=_EPS_sym_;
//DX20190905 [OBSOLETE]   //    // Calculate point group/space group to find correct tolerance for system
//DX20190905 [OBSOLETE]   //    //SYM::CalculatePointGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss);       
//DX20190905 [OBSOLETE]   //    //_eps_=a.sym_eps;
//DX20190905 [OBSOLETE]   //  }
//DX20190905 [OBSOLETE]   //  return BringInCell(a,_eps_);
//DX20190905 [OBSOLETE]   for(uint i=0; i<b.atoms.size(); i++){
//DX20190905 [OBSOLETE]     b.atoms[i].fpos = SYM::mod_one_xvec(b.atoms[i].fpos);
//DX20190905 [OBSOLETE]   } 
//DX20190905 [OBSOLETE]   return b;
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] xstructure BringInCell_20160101(const xstructure& a) {
//DX20190905 [OBSOLETE]   return BringInCell(a,_EPS_roundoff_);
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] //DX+CO END
//DX20190905 [OBSOLETE] 
//DX20190905 [OBSOLETE] //CO20190114 - DO NOT USE OVERLOADS OF BRINGINCELL() WITH EPSILON UNLESS YOU KNOW WHAT YOU ARE DOING
//DX20190905 [OBSOLETE] //DEFAULT TO OVERLOADS WITHOUT EPSILON WHICH USE HARD CUTOFF OF _ZERO_TOL_ = 1e-10
//DX20190905 [OBSOLETE] void xstructure::BringInCell(double epsilon) {
//DX20190905 [OBSOLETE]   _atom BringInCell(const _atom& atom_in,const xmatrix<double>& lattice,double epsilon);
//DX20190905 [OBSOLETE]   for(int i=0;i<(int)atoms.size();i++){
//DX20190905 [OBSOLETE]     atoms[i]=BringInCell(atoms[i],lattice,epsilon);
//DX20190905 [OBSOLETE]   }
//DX20190905 [OBSOLETE] }
//DX20190905 [OBSOLETE] void xstructure::BringInCell(void) {
//DX20190905 [OBSOLETE]   //DX+CO START
//DX20190905 [OBSOLETE]   _atom BringInCell(const _atom& atom_in,const xmatrix<double>& lattice);
//DX20190905 [OBSOLETE]   for(int i=0;i<(int)atoms.size();i++){
//DX20190905 [OBSOLETE]     atoms[i]=BringInCell(atoms[i],lattice);
//DX20190905 [OBSOLETE]   }
//DX20190905 [OBSOLETE]   //TOO SLOW
//DX20190905 [OBSOLETE]   //xstructure BringInCell(const xstructure& a);
//DX20190905 [OBSOLETE]   //*this=BringInCell(*this);
//DX20190905 [OBSOLETE]   //BringInCell(_EPS_sym_);  //CO
//DX20190905 [OBSOLETE]   //DX+CO END
//DX20190905 [OBSOLETE] }

// ***************************************************************************
// Function BringInCompact
// ***************************************************************************
// Make a structure where all the atoms are all the
// atoms are mapped through the unit and neighbors cells
// to minimize the shortest possible bond with an adjacent atom
// This option is very useful if you run big and complicate
// molecules where atoms exit of the unit cell and you have
// problems understanding where they are because visualization
// packages do not show bonds anymore ...
// Anyway, it is easier to test than to describe. (SC 6 Aug 04)
xstructure BringInCompact(const xstructure& a) {
  xstructure sstr=a;
  sstr.BringInCompact();
  return sstr;
}

xstructure _old_BringInCompact(const xstructure& a) {
  xstructure sstr=a;
  // double TOL=1e-15;
  double min_bond,bond;
  // For direct coordinates
  for(uint iat=0;iat<sstr.atoms.size();iat++)
    sstr.atoms.at(iat)=BringInCell(sstr.atoms.at(iat),sstr.lattice);
  xvector<double> adref1pos(3),adtstpos(3),adtrgpos(3);
  xvector<double> acref1pos(3),actstpos(3),actrgpos(3);
  for(uint i=1;i<sstr.atoms.size();i++) {        // scan all the atoms to move except the first
    adtrgpos=sstr.atoms.at(i).fpos;
    actrgpos=sstr.atoms.at(i).cpos;
    min_bond=1.0e6;
    for(uint ii=0;ii<i;ii++) {                      // scan over all the reference atoms
      adref1pos=sstr.atoms.at(ii).fpos;
      for(int ic=1;ic<=3;ic++) {                     // scan over all the reference atoms
        acref1pos(ic)=0.0;                          // scan over all the reference atoms
        for(int jc=1;jc<=3;jc++)
          acref1pos(ic)=acref1pos(ic)+adref1pos(jc)*sstr.lattice(jc,ic);
        // scan over all the reference atoms
      }                                             // scan over all the reference atoms

      for(int i1=-1;i1<=1;i1++)                     // roll over first neighbor cells
        for(int j1=-1;j1<=1;j1++)                   // roll over first neighbor cells
          for(int k1=-1;k1<=1;k1++) {               // roll over first neighbor cells
            adtstpos(1)=sstr.atoms.at(i).fpos(1)+i1;
            adtstpos(2)=sstr.atoms.at(i).fpos(2)+j1;
            adtstpos(3)=sstr.atoms.at(i).fpos(3)+k1;      // test the atom
            for(int ic=1;ic<=3;ic++) {                // test the atom
              actstpos(ic)=0.0;                      // test the atom
              for(int jc=1;jc<=3;jc++)               // test the atom
                actstpos(ic)=actstpos(ic)+adtstpos(jc)*sstr.lattice(jc,ic);   // test the atom
            }                                        // test the atom  
            bond=aurostd::modulus(actstpos-acref1pos);    // test the bond distance
            // if((bond<min_bond)                  // if it is OK then DO IT !
            if((bond<1.03*min_bond && abs((int) i-(int) ii)<10) || (bond<0.98*min_bond) )
            { //CO20200106 - patching for auto-indenting
              // if it is OK then DO IT !
              min_bond=bond;                                      // update
              adtrgpos=adtstpos;
              actrgpos=actstpos;
            }
          }
    }
    // now the TRG atom has all the information about the bst position of the atom to be moved.
    // than EAT IT !
    for(int j=1;j<=3;j++) {
      sstr.atoms.at(i).fpos(j)=adtrgpos(j);
      sstr.atoms.at(i).cpos(j)=actrgpos(j);
    }
  }
  // For Cartesian coordinates (get from direct coords);
  for(uint iat=0;iat<sstr.atoms.size();iat++)
    sstr.atoms.at(iat).cpos=F2C(sstr.lattice,sstr.atoms.at(iat).fpos);
  return sstr;
}

void xstructure::BringInCompact(void) {
  double min_bond,bond;
  // For direct coordinates
  BringInCell();
  xvector<double> adref1pos(3),adtstpos(3),adtrgpos(3);
  xvector<double> acref1pos(3),actstpos(3),actrgpos(3);
  for(uint i=1;i<atoms.size();i++) {        // scan all the atoms to move except the first
    adtrgpos=atoms.at(i).fpos;
    actrgpos=atoms.at(i).cpos;
    min_bond=1.0e6;
    for(uint ii=0;ii<i;ii++) {                      // scan over all the reference atoms
      adref1pos=atoms.at(ii).fpos;
      for(int ic=1;ic<=3;ic++) {                     // scan over all the reference atoms
        acref1pos(ic)=0.0;                          // scan over all the reference atoms
        for(int jc=1;jc<=3;jc++)
          acref1pos(ic)=acref1pos(ic)+adref1pos(jc)*lattice(jc,ic);
        // scan over all the reference atoms
      }                                             // scan over all the reference atoms

      for(int i1=-1;i1<=1;i1++)                     // roll over first neighbor cells
        for(int j1=-1;j1<=1;j1++)                   // roll over first neighbor cells
          for(int k1=-1;k1<=1;k1++) {               // roll over first neighbor cells
            adtstpos(1)=atoms.at(i).fpos(1)+i1;
            adtstpos(2)=atoms.at(i).fpos(2)+j1;
            adtstpos(3)=atoms.at(i).fpos(3)+k1;      // test the atom
            for(int ic=1;ic<=3;ic++) {                // test the atom
              actstpos(ic)=0.0;                      // test the atom        
              for(int jc=1;jc<=3;jc++)               // test the atom
                actstpos(ic)=actstpos(ic)+adtstpos(jc)*lattice(jc,ic);   // test the atom
            }                                        // test the atom  
            bond=aurostd::modulus(actstpos-acref1pos);    // test the bond distance
            // if((bond<min_bond)                  // if it is OK then DO IT !
            if((bond<1.03*min_bond && abs((int) i-(int) ii)<10) || (bond<0.98*min_bond) )
            { //CO20200106 - patching for auto-indenting
              // if it is OK then DO IT !
              min_bond=bond;                                      // update
              adtrgpos=adtstpos;
              actrgpos=actstpos;
            }
          }
    }
    // now the TRG atom has all the information about the bst position of the atom to be moved.
    // than EAT IT !
    for(int j=1;j<=3;j++) {
      atoms.at(i).fpos(j)=adtrgpos(j);
      atoms.at(i).cpos(j)=actrgpos(j);
    }
  }
  // For Cartesian coordinates (get from direct coords);
  for(uint iat=0;iat<atoms.size();iat++)
    atoms.at(iat).cpos=F2C(lattice,atoms.at(iat).fpos);
}

// ***************************************************************************
// Function BringInWignerSeitz
// ***************************************************************************
// Make a structure where all the atoms are
// mapped to their images in the Wigner-Seitz cell.(SC 10Jan04)
void xstructure::BringInWignerSeitz(void) {
  xvector<double> rat(3),a1(3),a2(3),a3(3),a12(3),a31(3),a23(3),a123(3);
  double na1,na2,na3,na12,na31,na23,na123;
  a1=lattice(1);na1=aurostd::modulus(a1);
  a2=lattice(2);na2=aurostd::modulus(a2);
  a3=lattice(3);na3=aurostd::modulus(a3);
  a12=lattice(1)+lattice(2);na12=aurostd::modulus(a12);
  a31=lattice(1)+lattice(3);na31=aurostd::modulus(a31);
  a23=lattice(2)+lattice(3);na23=aurostd::modulus(a23);
  a123=lattice(1)+lattice(2)+lattice(3);na123=aurostd::modulus(a123);
  double proj_a1,proj_a2,proj_a3,proj_a12,proj_a31,proj_a23,proj_a123;

  BringInCell();

  for(uint iat=0;iat<atoms.size();iat++) {
    for(int i=-1;i<=1;i++)
      for(int j=-1;j<=1;j++)
        for(int k=-1;k<=1;k++) {
          rat=atoms.at(iat).cpos+((double)i)*a1+((double)j)*a2+((double)k)*a3;
          proj_a1=scalar_product(rat,a1)/na1/na1;
          proj_a2=scalar_product(rat,a2)/na2/na2;
          proj_a3=scalar_product(rat,a3)/na3/na3;
          proj_a12=scalar_product(rat,a12)/na12/na12;
          proj_a31=scalar_product(rat,a31)/na31/na31;
          proj_a23=scalar_product(rat,a23)/na23/na23;
          proj_a123=scalar_product(rat,a123)/na123/na123;
          if((proj_a1>-0.5 && proj_a1<=0.5) &&
              (proj_a2>-0.5 && proj_a2<=0.5) &&
              (proj_a3>-0.5 && proj_a3<=0.5) &&
              (proj_a12>-0.5 && proj_a12<=0.5) &&
              (proj_a31>-0.5 && proj_a31<=0.5) &&
              (proj_a23>-0.5 && proj_a23<=0.5) &&
              (proj_a123>-0.5 && proj_a123<=0.5)) {
            atoms.at(iat).cpos(1)=rat(1);
            atoms.at(iat).cpos(2)=rat(2);
            atoms.at(iat).cpos(3)=rat(3);
            i=10;j=10;k=10;
          }
        }
    atoms.at(iat).fpos=C2F(lattice,atoms.at(iat).cpos);
  }
}

xstructure BringInWignerSeitz(const xstructure& a) {
  xstructure str=a;
  str.BringInWignerSeitz();
  return str;
}


// **************************************************************************
// Function GetPrimitive
// Function IsTranslationFVector& IsTranslationCVector
// ***************************************************************************
// This function returns 1 if the given vector is a translation
// vector of the given structure and 0 otherwise. Input translation
// vector is expectd to be in fractional coordinates.
// written by Stefano Curtarolo (superfast edition)
//#define IsTranslationFVector IsTranslationFVectorFAST
//#define IsTranslationFVector IsTranslationFVectorORIGINAL
//#define IsTranslationFVector IsTranslationFVectorFAST_2011
#define IsTranslationFVector IsTranslationFVectorORIGINAL_2011

// isTranslationVector() //DX20210316 - moved from XtalFinder
// faster than subsequent variants below and more robust
bool isTranslationVector(const xstructure& xstr, const xvector<double>& vec, double tolerance, bool is_frac){

  // Check if input vector is a translation vector (i.e. preserves periodicity)
  // tolerance: default is half an Angstrom
  // (tolerance example: need at least tol=0.1 for As1_ICSD_158474 == As1_ICSD_162840 via XtalFinder)

  if(tolerance<_ZERO_TOL_) {
    stringstream message;
    message << "Zero tolerance: " << tolerance;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
  }

  xvector<double> cvec, fvec;
  if(is_frac){
    fvec = vec;
    cvec = xstr.f2c*vec;
  }
  else{
    fvec = xstr.c2f*vec;
    cvec = vec;
  }

  uint natoms = xstr.atoms.size();
  bool skew = false;
  uint count=0;

  // ---------------------------------------------------------------------------
  // check if applying the translation maps to another atom
  // use MapAtom to match type/name/spin/occupation/etc.
  for(uint d=0;d<natoms;d++){
    _atom tmp_atom = xstr.atoms[d];
    tmp_atom.cpos+=cvec;
    tmp_atom.fpos+=fvec;
    if(SYM::MapAtom(xstr.atoms,tmp_atom,true,xstr.lattice,xstr.f2c,skew,tolerance)){
      count++;
    }
    // match not found, violates periodicity, return immediately
    else { return false; }
  }
  return (count == natoms);
}

bool IsTranslationFVectorFAST(const xstructure& a, const xvector<double>& ftvec) {
  stringstream message;
  if(a.equiv_fpos_epsilon<1.0e-12) {
    message << "Zero tolerance: " << a.equiv_fpos_epsilon;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
  }
  double tolerance=a.equiv_fpos_epsilon;
  if(aurostd::modulus(ftvec)<=tolerance) return TRUE;

  uint CntGoodTrans=0;
  xvector<double> ftpos(3),diff(3);
  bool found_atom=FALSE,found_tvec=TRUE;
  xvector<int> types(0,a.atoms.size());       // mirror to enhance speed !
  for(uint iat=0;iat<a.atoms.size();iat++)    // mirror to enhance speed !
    types(iat)=a.atoms.at(iat).type;          // mirror to enhance speed !
  for(uint iat=0;iat<a.atoms.size()&&found_tvec;iat++) {
    // Get translated position and shift back into unit cell.
    //  ftpos=a.atoms.at(iat).fpos+ftvec;
    // ftpos=BringInCell(ftpos,tolerance);  // no because it will be needed later...
    // Find closest atom in unit cell to translated position (should only be one atom less than tolerance).
    found_atom=FALSE;
    for(uint jat=0;jat<a.atoms.size()&&!found_atom&&found_tvec;jat++) {
      if(types[iat]==types[jat]) {
        //        diff=ftpos-a.atoms.at(jat).fpos;
        diff=BringInCell(a.atoms.at(iat).fpos-a.atoms.at(jat).fpos+ftvec,tolerance);
        // If the translated atom maps onto another and its type is the same as the atom mapped onto then increment.
        if(aurostd::modulus(diff)<=tolerance) {CntGoodTrans++;found_atom=TRUE;}
      }
    } // For jat
    if(found_atom==FALSE) found_tvec=FALSE; // at least one atom does not have a copy...
  } // For iat
  // If we counted a good translation for every atoms then we have a translation vector.
  // cout << "GOODs=" << CntGoodTrans << endl;
  if(CntGoodTrans==a.atoms.size()) return TRUE;
  return FALSE;
}

bool IsTranslationFVectorORIGINAL(const xstructure& a, const xvector<double>& ftvec) {
  stringstream message;
  if(a.equiv_fpos_epsilon<1.0e-12) {
    message << "Zero tolerance: " << a.equiv_fpos_epsilon;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
  }
  double tolerance=a.equiv_fpos_epsilon;
  if(aurostd::modulus(ftvec)<=tolerance) return TRUE;
  uint CntGoodTrans=0;
  xvector<double> ftpos(3),diff(3);
  xvector<int> types(0,a.atoms.size());       // mirror to enhance speed !
  for(uint iat=0;iat<a.atoms.size();iat++)    // mirror to enhance speed !
    types(iat)=a.atoms.at(iat).type;          // mirror to enhance speed !
  for(uint iat=0;iat<a.atoms.size();iat++) {
    // Get translated position and shift back into unit cell.
    // ftpos=a.atoms.at(iat).fpos+ftvec;
    // ftpos=BringInCell(ftpos,tolerance); // no because it will be needed later...
    // Find closest atom in unit cell to translated position (should only be one atom less than tolerance).
    for(uint jat=0;jat<a.atoms.size();jat++) {
      if(types[iat]==types[jat]) {
        //        diff=ftpos-a.atoms.at(jat).fpos;
        diff=BringInCell(a.atoms.at(iat).fpos-a.atoms.at(jat).fpos+ftvec,tolerance);
        // If the translated atom maps onto another and its type is the same as the atom mapped onto then increment.
        if(aurostd::modulus(diff)<=tolerance) CntGoodTrans++;
      }
    } // For jat
  } // For iat
  // If we counted a good translation for every atoms then we have a translation vector.
  //  cout << "GOODs=" << CntGoodTrans << endl;
  if(CntGoodTrans==a.atoms.size()) return TRUE;
  return FALSE;
}

bool IsTranslationFVectorFAST_2011(const xstructure& a, const xvector<double>& ftvec) {
  stringstream message;
  if(a.equiv_fpos_epsilon<1.0e-12) {
    message << "Zero tolerance: " << a.equiv_fpos_epsilon;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
  }
  double tolerance=a.equiv_fpos_epsilon;
  if(aurostd::modulus(ftvec)<=tolerance) return TRUE;
  uint CntGoodTrans=0;
  xvector<double> ftpos(3),diff(3);
  bool found_atom=FALSE,found_tvec=TRUE;
  xvector<int> types(0,a.atoms.size());       // mirror to enhance speed !
  for(uint iat=0;iat<a.atoms.size();iat++)    // mirror to enhance speed !
    types(iat)=a.atoms.at(iat).type;          // mirror to enhance speed !
  for(uint iat=0;iat<a.atoms.size()&&found_tvec;iat++) {
    // Get translated position and shift back into unit cell.
    ftpos=a.atoms.at(iat).fpos+ftvec;
    ftpos=BringInCell(ftpos);
    // Find closest atom in unit cell to translated position (should only be one atom less than tolerance).
    found_atom=FALSE;
    for(uint jat=0;jat<a.atoms.size()&&!found_atom&&found_tvec;jat++) {
      if(types[iat]==types[jat]) {
        diff=ftpos-a.atoms.at(jat).fpos;
        // If the translated atom maps onto another and its type is the same as the atom mapped onto then increment.
        if(aurostd::modulus(diff)< tolerance ) {CntGoodTrans++;found_atom=TRUE;}
      }
    } // For jat
    if(found_atom==FALSE) found_tvec=FALSE; // at least one atom does not have a copy...
  } // For iat
  // If we counted a good translation for every atoms then we have a translation vector.
  if(CntGoodTrans==a.atoms.size()) return TRUE;
  return FALSE;
}

bool IsTranslationFVectorORIGINAL_2011(const xstructure& a, const xvector<double>& ftvec) {
  stringstream message;
  if(a.equiv_fpos_epsilon<1.0e-12) {
    message << "Zero tolerance: " << a.equiv_fpos_epsilon;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
  }
  double tolerance=a.equiv_fpos_epsilon;
  if(aurostd::modulus(ftvec)<=tolerance) return TRUE;

  uint CntGoodTrans=0;
  xvector<double> ftpos(3),diff(3);
  xvector<int> types(0,a.atoms.size());       // mirror to enhance speed !
  for(uint iat=0;iat<a.atoms.size();iat++)    // mirror to enhance speed !
    types(iat)=a.atoms.at(iat).type;          // mirror to enhance speed !
  for(uint iat=0;iat<a.atoms.size();iat++) {
    // Get translated position and shift back into unit cell.
    ftpos=a.atoms.at(iat).fpos+ftvec;
    ftpos=BringInCell(ftpos);
    // Find closest atom in unit cell to translated position (should only be one atom less than tolerance).
    for(uint jat=0;jat<a.atoms.size();jat++) {
      if(types[iat]==types[jat]) {
        diff=ftpos-a.atoms.at(jat).fpos;
        diff = SYM::minimizeDistanceFractionalMethod(diff); //DX20190613
        //DX20190613 [OBSOLETE] SYM::PBC(diff); //DX20190317 - need to use PBC for comparing difference vectors
        // If the translated atom maps onto another and its type is the same as the atom mapped onto then increment.
        if(aurostd::modulus(diff)< tolerance ) CntGoodTrans++;
      }
    } // For jat
  } // For iat
  // If we counted a good translation for every atoms then we have a translation vector.
  if(CntGoodTrans==a.atoms.size()) return TRUE;
  return FALSE;
}


bool IsTranslationCVector(const xstructure& a, const xvector<double>& ctvec) {
  stringstream message;
  if(a.equiv_fpos_epsilon<1.0e-12) {
    message << "Zero tolerance: " << a.equiv_fpos_epsilon;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
  }
  /* Input translation vector is expectd to be in cartesian coordinates. */
  return IsTranslationFVector(a,C2F(a.lattice,ctvec));
}

// **************************************************************************
// Function GetPrimitive
// **************************************************************************
// This funtion returns a structure object where everything
// has been transformed into a primitive lattice.
// Algorithm
// Construct a list of all candidate primitive lattice vectors.
// This list consists of the present lattice vectors and all basis
// vectors where translation by that basis vector leaves all the
// basis atom types unchanged. Loop over all possible triads in this list.
// For each triad find the volume. Store a list of all triads with the
// minimum volume. Then select the triad for which the first triad vector
// has the maximal projection onto the first lattice vector.  If there
// are more than one of these do the same for the second lattice vector,
// and then the third if needed. Then just pick one if there are still
// multiple candidates.

xstructure GetPrimitiveSINGLE(const xstructure& a,double tolerance);
xstructure GetPrimitiveMULTITHREAD(const xstructure& a,double tolerance);
xstructure GetPrimitive(const xstructure& a);

xstructure GetPrimitive(const xstructure& _a,double tolerance) {
  // return GetPrimitiveMULTITHREAD(a);

  xstructure b;
  b=GetPrimitive_20210322(_a, tolerance); //DX20210406 - new/fast variant

  //[DX TESTING] xstructure c;
  //[DX TESTING] c=GetPrimitiveMULTITHREAD(a,tolerance);
  //[DX TESTING] cerr << "STRUCTURES MATCH (ORIG vs OLD): " << compare::structuresMatch(a,b,false) << " (ORIG vs NEW): " << compare::structuresMatch(a,c,false) << endl;
  //[DX TESTING] cerr << "VOLUMES: " << b.atoms.size()/c.atoms.size() << " (mod=" << b.atoms.size()%c.atoms.size() << ")" << endl;

  // OLD VERSION - the version above is fast and does not need to be multithreaded
  //DX20210406 [OBSOLETE] bool LDEBUG=(FALSE || XHOST.DEBUG);
  //DX20210406 [OBSOLETE] if(LDEBUG) cerr << "GetPrimitive(const xstructure& _a,double tolerance): BEGIN" << endl;

  //DX20210406 [OBSOLETE] bool bool_backup=AFLOW_PTHREADS::FLAG;
  //DX20210406 [OBSOLETE] uint uint_backup=AFLOW_PTHREADS::MAX_PTHREADS;
  //DX20210406 [OBSOLETE] xstructure a(_a);
  //DX20210406 [OBSOLETE] if(tolerance<=0.0) a.equiv_fpos_epsilon=_EQUIV_FPOS_EPS_; else a.equiv_fpos_epsilon=tolerance;

  //DX20210406 [OBSOLETE] if(AFLOW_PTHREADS::RUNNING==0 && !AFLOW_PTHREADS::FLAG) { // force to save time
  //DX20210406 [OBSOLETE]   AFLOW_PTHREADS::FLAG=TRUE;   // put multithread by default in
  //DX20210406 [OBSOLETE]   AFLOW_PTHREADS::MAX_PTHREADS=XHOST.CPU_Cores;
  //DX20210406 [OBSOLETE] }
  //DX20210406 [OBSOLETE] // calc
  //DX20210406 [OBSOLETE] xstructure b;
  //DX20210406 [OBSOLETE] b=GetPrimitiveMULTITHREAD(a,tolerance);

  //DX20210406 [OBSOLETE] //if(AFLOW_PTHREADS::FLAG==TRUE) {
  //DX20210406 [OBSOLETE] //  //    cerr << "GetPrimitiveMULTITHREAD a.equiv_fpos_epsilon=" << a.equiv_fpos_epsilon << endl;
  //DX20210406 [OBSOLETE] //  b=GetPrimitiveMULTITHREAD(a,tolerance);
  //DX20210406 [OBSOLETE] //} else {
  //DX20210406 [OBSOLETE] //  //    cerr << "GetPrimitiveSINGLE a.equiv_fpos_epsilon=" << a.equiv_fpos_epsilon << endl;
  //DX20210406 [OBSOLETE] //  b=GetPrimitiveSINGLE(a,tolerance);
  //DX20210406 [OBSOLETE] //}
  //DX20210406 [OBSOLETE] // restore and return
  //DX20210406 [OBSOLETE] AFLOW_PTHREADS::FLAG=bool_backup;
  //DX20210406 [OBSOLETE] AFLOW_PTHREADS::MAX_PTHREADS=uint_backup;

  return b;
}

xstructure GetPrimitive(const xstructure& a) {
  //DX20210401 [OBSOLETE] return GetPrimitive(a,_EQUIV_FPOS_EPS_);
  double tolerance=AUROSTD_MAX_DOUBLE;
  if(a.sym_eps!=AUROSTD_NAN){ tolerance = a.sym_eps; }
  else{ tolerance = SYM::defaultTolerance(a); }
  return GetPrimitive(a,tolerance);
}

//ME20220208 [OBSOLETE - not used] pthread_mutex_t mutex_XATOM=PTHREAD_MUTEX_INITIALIZER;
//ME20220208 [OBSOLETE - not used] #define _PTHREAD_FLUSH_TIME_ 1
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used] typedef struct {
//ME20220208 [OBSOLETE - not used]   int      ITHREAD;
//ME20220208 [OBSOLETE - not used]   int      THREADS_MAX;
//ME20220208 [OBSOLETE - not used]   xstructure *pstr;
//ME20220208 [OBSOLETE - not used]   int      ispecie_min;
//ME20220208 [OBSOLETE - not used]   vector<xvector<double> > *ptvector;
//ME20220208 [OBSOLETE - not used]   xmatrix<double> *polattice;
//ME20220208 [OBSOLETE - not used]   int      itbusy;
//ME20220208 [OBSOLETE - not used]   int      step;
//ME20220208 [OBSOLETE - not used] } _threaded_GETTVECTORS_params;
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used] void *_threaded_GetTvectors(void *ptr) {
//ME20220208 [OBSOLETE - not used]   bool LDEBUG=(FALSE || XHOST.DEBUG);
//ME20220208 [OBSOLETE - not used]   _threaded_GETTVECTORS_params* pparams;
//ME20220208 [OBSOLETE - not used]   xmatrix<double> plattice(3,3);
//ME20220208 [OBSOLETE - not used]   xvector<double> fdisp(3),cdisp(3);
//ME20220208 [OBSOLETE - not used]   pparams = (_threaded_GETTVECTORS_params*) ptr;
//ME20220208 [OBSOLETE - not used]   (pparams->itbusy)=TRUE;
//ME20220208 [OBSOLETE - not used]   // pthread_mutex_lock(&mutex_XATOM);
//ME20220208 [OBSOLETE - not used]   AFLOW_PTHREADS::RUNNING++;
//ME20220208 [OBSOLETE - not used]   // cerr << "AFLOW_PTHREADS::RUNNING=" << AFLOW_PTHREADS::RUNNING << endl;
//ME20220208 [OBSOLETE - not used]   // pthread_mutex_unlock(&mutex_XATOM);
//ME20220208 [OBSOLETE - not used]   // CODE BEGIN
//ME20220208 [OBSOLETE - not used]   //  cerr << "debug " << (pparams->ITHREAD) << "/" << (pparams->THREADS_MAX) << endl;
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used]   double sstr_volume=(*pparams->pstr).Volume();
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used]   // CODE STEP1 BEGIN
//ME20220208 [OBSOLETE - not used]   if(pparams->step==1) {
//ME20220208 [OBSOLETE - not used]     if(LDEBUG) pthread_mutex_lock(&mutex_XATOM);   // LOCK
//ME20220208 [OBSOLETE - not used]     if(LDEBUG) cerr << "*_threaded_GetTvectors STEP1: " << AFLOW_PTHREADS::RUNNING << " " << pparams->ITHREAD << endl;
//ME20220208 [OBSOLETE - not used]     if(LDEBUG) pthread_mutex_unlock(&mutex_XATOM); // UNLOCK
//ME20220208 [OBSOLETE - not used]     for(uint iat1=0+(pparams->ITHREAD);iat1<(*pparams->pstr).atoms.size();iat1+=(pparams->THREADS_MAX)) { // does modulus thread max
//ME20220208 [OBSOLETE - not used]       if((*pparams->pstr).atoms.at(iat1).type==(pparams->ispecie_min))
//ME20220208 [OBSOLETE - not used]         for(uint iat2=0;iat2<(*pparams->pstr).atoms.size();iat2++)
//ME20220208 [OBSOLETE - not used]           if((*pparams->pstr).atoms.at(iat2).type==(pparams->ispecie_min)) {
//ME20220208 [OBSOLETE - not used]             fdisp=(*pparams->pstr).atoms.at(iat2).fpos-(*pparams->pstr).atoms.at(iat1).fpos;
//ME20220208 [OBSOLETE - not used]             cdisp=(*pparams->pstr).atoms.at(iat2).cpos-(*pparams->pstr).atoms.at(iat1).cpos;
//ME20220208 [OBSOLETE - not used]             if(aurostd::modulus(fdisp)>0.01 && aurostd::modulus(cdisp)>0.01)
//ME20220208 [OBSOLETE - not used]               if(IsTranslationFVector((*pparams->pstr),fdisp)) {
//ME20220208 [OBSOLETE - not used]                 pthread_mutex_lock(&mutex_XATOM);
//ME20220208 [OBSOLETE - not used]                 (*pparams->ptvector).push_back(cdisp);
//ME20220208 [OBSOLETE - not used]                 pthread_mutex_unlock(&mutex_XATOM);
//ME20220208 [OBSOLETE - not used]               }
//ME20220208 [OBSOLETE - not used]           }
//ME20220208 [OBSOLETE - not used]     }
//ME20220208 [OBSOLETE - not used]   } // CODE STEP1 END
//ME20220208 [OBSOLETE - not used]   // CODE STEP2 BEGIN
//ME20220208 [OBSOLETE - not used]   if(pparams->step==2) {
//ME20220208 [OBSOLETE - not used]     if(LDEBUG) pthread_mutex_lock(&mutex_XATOM);   // LOCK
//ME20220208 [OBSOLETE - not used]     if(LDEBUG) cerr << "*_threaded_GetTvectors STEP2: " << AFLOW_PTHREADS::RUNNING << " " << pparams->ITHREAD << endl;
//ME20220208 [OBSOLETE - not used]     if(LDEBUG) pthread_mutex_unlock(&mutex_XATOM); // UNLOCK
//ME20220208 [OBSOLETE - not used]     for(uint iu=0+(pparams->ITHREAD);iu<(*pparams->ptvector).size();iu+=(pparams->THREADS_MAX)) { // does modulus thread max
//ME20220208 [OBSOLETE - not used]       for(uint i=1;i<=3;i++) plattice[1][i]=(*pparams->ptvector).at(iu)[i];
//ME20220208 [OBSOLETE - not used]       for(uint iv=0;iv<(*pparams->ptvector).size()&& iv!=iu;iv++) {
//ME20220208 [OBSOLETE - not used]         for(uint i=1;i<=3;i++) plattice[2][i]=(*pparams->ptvector).at(iv)[i];
//ME20220208 [OBSOLETE - not used]         for(uint iw=0;iw<(*pparams->ptvector).size()&& iw!=iv && iw!=iu;iw++) {
//ME20220208 [OBSOLETE - not used]           for(uint i=1;i<=3;i++) plattice[3][i]=(*pparams->ptvector).at(iw)[i];
//ME20220208 [OBSOLETE - not used]           if(det(plattice)>0.999 && det(plattice)<sstr_volume) {   // no coplanar and contain at least 1 atom and smaller than the original cell
//ME20220208 [OBSOLETE - not used]             if(aurostd::isinteger(sstr_volume/det(plattice))) {    // integer ratio of volumes
//ME20220208 [OBSOLETE - not used]               if(det(plattice)<det((*pparams->polattice))) {                    // better than before
//ME20220208 [OBSOLETE - not used]                 if(LDEBUG) cout << XPID << "DEBUG"<<iu<<","<<iv<<","<<iw<<" "<< sstr_volume<<" "<<det(plattice)<<" "<<sstr_volume/det(plattice)<<endl;
//ME20220208 [OBSOLETE - not used]                 if(isdifferent(plattice,(*pparams->polattice),0.0001)) {
//ME20220208 [OBSOLETE - not used]                   plattice=MinkowskiBasisReduction(plattice);      // Minkowski first
//ME20220208 [OBSOLETE - not used]                   plattice=NiggliUnitCellForm(plattice);           // Niggli Second
//ME20220208 [OBSOLETE - not used]                   if(isdifferent(plattice,(*pparams->polattice),0.0001)) {
//ME20220208 [OBSOLETE - not used]                     pthread_mutex_lock(&mutex_XATOM);   // LOCK
//ME20220208 [OBSOLETE - not used]                     (*pparams->polattice)=plattice;
//ME20220208 [OBSOLETE - not used]                     pthread_mutex_unlock(&mutex_XATOM); // UNLOCK
//ME20220208 [OBSOLETE - not used]                   }
//ME20220208 [OBSOLETE - not used]                 }
//ME20220208 [OBSOLETE - not used]               }
//ME20220208 [OBSOLETE - not used]             }
//ME20220208 [OBSOLETE - not used]           }
//ME20220208 [OBSOLETE - not used]         }
//ME20220208 [OBSOLETE - not used]       }
//ME20220208 [OBSOLETE - not used]     }
//ME20220208 [OBSOLETE - not used]   } // CODE STEP2 END
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used]   (pparams->itbusy)=FALSE;
//ME20220208 [OBSOLETE - not used]   AFLOW_PTHREADS::RUNNING--;
//ME20220208 [OBSOLETE - not used]   aurostd::Sleep(_PTHREAD_FLUSH_TIME_);
//ME20220208 [OBSOLETE - not used]   return NULL;
//ME20220208 [OBSOLETE - not used] }

// **************************************************************************
// GetPrimitive() //DX20210406
// **************************************************************************
// This version is faster than the previous AFLOW variants
// Speed ups:
//   - returns immediately if the number of atoms/types indicate the cell
//     cannot be reduced (i.e., only one atom or one atom of a given type)
//   - optimized lattice vector search (search over only the least-frequent
//     atom type)
//   - optimized lattice search (volume checks, mulitplicity checks, etc.)
//   - uses transfomation method to convert between original and primitive
//     cell (as opposed to removing atoms)
// It also uses sym_eps to reduce the cell (previous methods did not)
xstructure GetPrimitive_20210322(const xstructure& a,double eps) {  //DX20210406
  xstructure xstr = a;
  xstr.GetPrimitive_20210322(eps);
  return xstr;
}

void xstructure::GetPrimitive_20210322(double eps) { //DX20210406

  bool LDEBUG=(FALSE || XHOST.DEBUG);

  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " BEGIN " << endl; }

  uint natoms_orig = (*this).atoms.size();
  // ---------------------------------------------------------------------------
  // if only one atom or atom type in the unit cell, the structure is already
  // primitivized
  if(natoms_orig == 1 || aurostd::min((*this).num_each_type) == 1) { return; }

  (*this).ReScale(1.0);
  (*this).FixLattices(); //DX20210407 - since we use c2f/f2c, update for safety

  // ---------------------------------------------------------------------------
  // set tolerance
  double tolerance = eps;
  if((*this).dist_nn_min==AUROSTD_NAN){ (*this).MinDist(); }

  if(tolerance==AUROSTD_MAX_DOUBLE){
    if((*this).sym_eps!=AUROSTD_NAN){ tolerance = (*this).sym_eps; }
    else{ tolerance = SYM::defaultTolerance((*this)); }
  }
  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [1] Tolerance = " << tolerance << endl; }

  double volume_orig=(*this).Volume();

  // ---------------------------------------------------------------------------
  // get least frequent atom type and the corresponding set of atoms
  // to search for possible lattice vectors (minimal set of atoms perserving
  // periodicity)
  uint atom_type_min = getLeastFrequentAtomTypes((*this))[0]; // normally a vector, grabbing only first one (there will always be one)
  vector<uint> vindices_atoms_min = getAtomIndicesByType((*this),atom_type_min);
  uint natoms_min = vindices_atoms_min.size();

  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [2] Subset of atoms to find lattice vectors: " << natoms_min << endl; }

  // generate list of vectors
  vector<xvector<double> > candidate_lattice_vector;
  candidate_lattice_vector.push_back((*this).lattice(1));  // lattice is made of good vectors
  candidate_lattice_vector.push_back((*this).lattice(2));  // lattice is made of good vectors
  candidate_lattice_vector.push_back((*this).lattice(3));  // lattice is made of good vectors

  // ---------------------------------------------------------------------------
  // get all lattice vectors
  // only need to check difference between 0th and ith atom coordinates for the
  // subset of atom indices
  vector<xvector<double> > diff_vectors;
  for(uint i=1; i<natoms_min; i++) {
    diff_vectors.push_back(::BringInCell((*this).atoms[vindices_atoms_min[i]].fpos - (*this).atoms[vindices_atoms_min[0]].fpos)); //need to get BringInCell from outside xstructure scope
  }
  //Translate by difference vectors and check if equivalent
  bool is_frac = true;
  for(uint d=0; d < diff_vectors.size(); d++){
    if(isTranslationVector((*this), diff_vectors[d], tolerance, is_frac)){
      candidate_lattice_vector.push_back((*this).f2c*diff_vectors[d]);
    }
  }

  double nlattice_vectors = candidate_lattice_vector.size();
  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " [3] number of lattice vectors=" << nlattice_vectors << endl; }

  if(LDEBUG){
    for(uint i=0;i<nlattice_vectors;i++){
      cerr << __AFLOW_FUNC__ << i << ": mod=" << aurostd::modulus(candidate_lattice_vector[i]) << " vec=" << candidate_lattice_vector[i] << endl;
    }
  }

  // ---------------------------------------------------------------------------
  // if only the three original lattice vectors were found, then the cell
  // is already primitivized
  if(nlattice_vectors == 3){ return; }

  xmatrix<double> lattice_tmp(3,3), plattice(3,3), olattice(3,3);
  olattice=(*this).lattice;                           // the lattice is always a good lattice
  double volume_min = aurostd::det((*this).lattice);
  double volume_tmp = AUROSTD_MAX_DOUBLE;
  double amin = (*this).dist_nn_min;
  double tol_vol = tolerance * amin * amin * amin;    // scale tolerance for volume
  double atom_number_ratio = 1;                       // track atom ratio during reduction
  double integer_tol = 0.05;                          // use somewhat large tol to find lattices
  double moduli_sum=AUROSTD_MAX_DOUBLE, moduli_sum_tmp=AUROSTD_MAX_DOUBLE; // ensures shortest sum of lattice vectors (chooses minimum Minkowski lattice, some noise can get introduced)

  // ---------------------------------------------------------------------------
  // check possible lattices, upper triangular for-loop only
  // criteria for volume of new lattice:
  //  1) must not be zero
  //  2) must be smaller than original lattice volume
  //  3) must be an integer multiple of original lattice
  //  4) must reduce atom count consistent with integer multiple
  //  5) must not remove entire atom types
  //  6) must form shortest vectors
  // Minkowski/Niggli will fix left-handed (negative determinants)
  for(uint iu=0;iu<nlattice_vectors;iu++) { for(uint i=1;i<=3;i++) lattice_tmp[1][i]=candidate_lattice_vector.at(iu)[i];
    for(uint iv=iu+1;iv<nlattice_vectors;iv++) { for(uint i=1;i<=3;i++) lattice_tmp[2][i]=candidate_lattice_vector.at(iv)[i];
      for(uint iw=iv+1;iw<nlattice_vectors;iw++) { for(uint i=1;i<=3;i++) lattice_tmp[3][i]=candidate_lattice_vector.at(iw)[i];
        volume_tmp = abs(lattice_tmp[1][1]*lattice_tmp[2][2]*lattice_tmp[3][3]+lattice_tmp[1][2]*lattice_tmp[2][3]*lattice_tmp[3][1]+ // FAST
            lattice_tmp[1][3]*lattice_tmp[2][1]*lattice_tmp[3][2]-lattice_tmp[1][3]*lattice_tmp[2][2]*lattice_tmp[3][1]-              // FAST
            lattice_tmp[1][2]*lattice_tmp[2][1]*lattice_tmp[3][3]-lattice_tmp[1][1]*lattice_tmp[2][3]*lattice_tmp[3][2]);             // FAST
        if(aurostd::abs(volume_tmp) > tol_vol){
          plattice=::MinkowskiBasisReduction(lattice_tmp);    // Minkowski first, "::" is needed to access outside xstructure scope
          plattice=::NiggliUnitCellForm(plattice);            // Niggli Second, "::" is needed to access outside xstructure scope
          volume_tmp = (plattice[1][1]*plattice[2][2]*plattice[3][3]+plattice[1][2]*plattice[2][3]*plattice[3][1]+  // FAST
              plattice[1][3]*plattice[2][1]*plattice[3][2]-plattice[1][3]*plattice[2][2]*plattice[3][1]-              // FAST
              plattice[1][2]*plattice[2][1]*plattice[3][3]-plattice[1][1]*plattice[2][3]*plattice[3][2]);             // FAST
          atom_number_ratio = (double)natoms_min/aurostd::nint(volume_orig/volume_tmp);
          moduli_sum_tmp = aurostd::modulus(plattice(1)) + aurostd::modulus(plattice(2)) + aurostd::modulus(plattice(3));
          if(volume_tmp > tol_vol &&                                                 // 1) tol > 0
            volume_tmp < volume_orig &&                                            // 2) new_tol < orig_tol
              aurostd::isinteger(volume_orig / volume_tmp, integer_tol) &&           // 3) new and orig volumes related by integer multiple
              aurostd::isinteger(atom_number_ratio, 0.05) &&                         // 4) number of atoms is consistent with integer multiple
              (atom_number_ratio-1.0)>-_ZERO_TOL_ &&                                 // 5) didn't lose any atoms
              moduli_sum_tmp < moduli_sum){                                          // 6) ensure shortest sum of vectors (Minkowski should ensure this, but noise can be introduced)
                moduli_sum = moduli_sum_tmp;
                olattice=plattice;
                volume_min = (olattice[1][1]*olattice[2][2]*olattice[3][3]+olattice[1][2]*olattice[2][3]*olattice[3][1]+   // FAST
                    olattice[1][3]*olattice[2][1]*olattice[3][2]-olattice[1][3]*olattice[2][2]*olattice[3][1]-             // FAST
                    olattice[1][2]*olattice[2][1]*olattice[3][3]-olattice[1][1]*olattice[2][3]*olattice[3][2]);            // FAST
              }
        }
      }
    }
  }
  plattice=olattice;

  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " reduced lattice=" << plattice << endl; }

  // ---------------------------------------------------------------------------
  // if the lattice remains the same, do not change or update the atoms
  if(aurostd::isequal(volume_min, volume_orig) || aurostd::isequal(plattice,(*this).lattice)){ return; }

  // ---------------------------------------------------------------------------
  //DX20210316 - used transformation method (more efficient)
  // ---------------------------------------------------------------------------
  // update xstructure
  xstructure prim = (*this);
  //prim.FixLattices();
  xmatrix<double> transformation_matrix = GetBasisTransformation((*this).lattice, plattice);
  xmatrix<double> rotation_matrix = aurostd::eye<double>(3,3);
  // try to primitivize, it may fail, so return original structure
  try{
    // this checks volumes, number of atoms, etc. internally
    prim.TransformStructure(transformation_matrix, rotation_matrix);
  }
  catch(aurostd::xerror& err){  //CO20230220 - patching names vs. types bug
    if(aurostd::substring2bool(err.whereFileName(),"aflow_xatom.cpp") && 
        aurostd::substring2bool(err.whereFunction(),"xstructure::ReplaceAtoms():") && 
        aurostd::substring2bool(err.what(),"atoms and species were incorrectly (re)sorted")){
      throw err;
    }
    return;
  }

  // no messed up volume
  // this is checked in TransformStructure, but kept as a safety
  // DX20210623 - the check has been improved:
  // first, check if the reduction factor (inverse of fraction) is an integer
  // then, check if that factor is consistent with the number of atoms
  // this method is less sensitive to the tolerance threshold
  double fraction=aurostd::det(prim.lattice)/aurostd::det((*this).lattice);
  double reduction_factor = 1.0/fraction;
  if(!aurostd::isinteger(reduction_factor,0.1)){
    stringstream message;
    message << "The original volume is not an integer multiple of the new volume: " << reduction_factor;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
  }
  uint reduction_factor_integer = (uint)round(reduction_factor);
  if(abs(natoms_orig-(double)reduction_factor_integer*prim.atoms.size())>0.1) {
    stringstream message;
    message << "ERROR   " << __AFLOW_FUNC__ << endl;
    message << "        supercell has the wrong number of atoms" << endl;
    message << "        volume original    = " << (*this).Volume() << endl;
    message << "        volume prim        = " << prim.Volume() << endl;
    message << "        a.scale            = " << (*this).scale << endl;
    message << "        b.scale            = " << prim.scale << endl;
    message << "        a.atoms.size()     = " << (*this).atoms.size() << endl;
    message << "        b.atoms.size()     = " << prim.atoms.size() << endl;
    message << "        fraction           = " << fraction << endl;
    message << "        reduction_factor   = " << reduction_factor << endl; //DX20210623
    message << "        supercell atoms    = " << reduction_factor*prim.atoms.size() << endl; //DX20210623
    message << prim << endl;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
  }
  // everything ok
  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " END [ok]=" << fraction << endl; }

  // set primitive representation
  (*this) = prim;
}

//ME20220208 [OBSOLETE - not used] xstructure GetPrimitiveMULTITHREAD(const xstructure& _a,double tolerance) {  // APRIL 2009 JUNE 2012 added tolerance
//ME20220208 [OBSOLETE - not used]   bool LDEBUG=(FALSE || XHOST.DEBUG);
//ME20220208 [OBSOLETE - not used]   if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN " << endl; //CO20200201
//ME20220208 [OBSOLETE - not used]   if(LDEBUG){cerr << __AFLOW_FUNC__ << " _a=" << endl;cerr << _a << endl;} //CO20200201
//ME20220208 [OBSOLETE - not used]   cout.setf(std::ios::fixed,std::ios::floatfield);
//ME20220208 [OBSOLETE - not used]   cout.precision(10);
//ME20220208 [OBSOLETE - not used]   xstructure a(_a);
//ME20220208 [OBSOLETE - not used]   xstructure sstr=a;
//ME20220208 [OBSOLETE - not used]   if(LDEBUG){cerr << __AFLOW_FUNC__ << " sstr=" << endl;cerr << sstr << endl;} //CO20200201
//ME20220208 [OBSOLETE - not used]   if(LDEBUG){cerr << __AFLOW_FUNC__ << " sstr.atoms.size()=" << sstr.atoms.size() << endl;}  //CO20200201
//ME20220208 [OBSOLETE - not used]   sstr.SetVolume(sstr.atoms.size());
//ME20220208 [OBSOLETE - not used]   if(LDEBUG){cerr << __AFLOW_FUNC__ << " sstr(post Vol)=" << endl;cerr << sstr << endl;} //CO20200201
//ME20220208 [OBSOLETE - not used]   sstr=ReScale(sstr,1.0);
//ME20220208 [OBSOLETE - not used]   if(LDEBUG){cerr << __AFLOW_FUNC__ << " sstr(pre BringInCell)=" << endl;cerr << sstr << endl;} //CO20200201
//ME20220208 [OBSOLETE - not used]   sstr=BringInCell(sstr);
//ME20220208 [OBSOLETE - not used]   if(LDEBUG){cerr << __AFLOW_FUNC__ << " sstr(post BringInCell)=" << endl;cerr << sstr << endl;} //CO20200201
//ME20220208 [OBSOLETE - not used]   //  sstr.CalculateSymmetry();
//ME20220208 [OBSOLETE - not used]   if(tolerance<=0.0) a.equiv_fpos_epsilon=_EQUIV_FPOS_EPS_; else a.equiv_fpos_epsilon=tolerance;
//ME20220208 [OBSOLETE - not used]   double sstr_volume=sstr.Volume();
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used]   if(LDEBUG) cerr << __AFLOW_FUNC__ << " [1] " << endl; //CO20200201
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used]   _aflags aflags;
//ME20220208 [OBSOLETE - not used]   // identify the minimum set of atoms
//ME20220208 [OBSOLETE - not used]   // bool PGROUPWRITE=FALSE,PGROUPKWRITE=FALSE,FGROUPWRITE=FALSE,IATOMSWRITE=FALSE;
//ME20220208 [OBSOLETE - not used]   // bool OSSWRITE=TRUE; // to FileMESSAGE, does not matter as it is /dev/null
//ME20220208 [OBSOLETE - not used]   // ofstream FileMESSAGE("/dev/stderr");
//ME20220208 [OBSOLETE - not used]   // SYM::CalculatePointGroup(FileMESSAGE,sstr,aflags,PGROUPWRITE,OSSWRITE,cout);
//ME20220208 [OBSOLETE - not used]   // SYM::CalculatePointGroupKLattice(FileMESSAGE,sstr,aflags,PGROUPKWRITE,OSSWRITE,cout);
//ME20220208 [OBSOLETE - not used]   // SYM::CalculateFactorGroup(FileMESSAGE,sstr,aflags,FGROUPWRITE,OSSWRITE,cout);
//ME20220208 [OBSOLETE - not used]   // SYM::CalculateInequivalentAtoms(FileMESSAGE,sstr,aflags,IATOMSWRITE,OSSWRITE,cout);
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used]   xmatrix<double> plattice(3,3),olattice(3,3);
//ME20220208 [OBSOLETE - not used]   xvector<double> fdisp(3),cdisp(3);
//ME20220208 [OBSOLETE - not used]   std::vector<xvector<double> > candidate_lattice_vector;
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used]   if(LDEBUG) cerr << __AFLOW_FUNC__ << " [2] " << endl; //CO20200201
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used]   int specie_min=sstr.num_each_type.at(0),ispecie_min=0,specie_min_threshold=_PRIM_MULTITHREAD_MIN_ATOMS_THRESHOLD_; // seems a good threshold
//ME20220208 [OBSOLETE - not used]   for(uint ispecie=0;ispecie<sstr.num_each_type.size();ispecie++)
//ME20220208 [OBSOLETE - not used]     if(sstr.num_each_type.at(ispecie)<specie_min) {
//ME20220208 [OBSOLETE - not used]       specie_min=sstr.num_each_type.at(ispecie);
//ME20220208 [OBSOLETE - not used]       ispecie_min=ispecie;
//ME20220208 [OBSOLETE - not used]     }
//ME20220208 [OBSOLETE - not used]   // cerr << "DEBUG specie_min=" << specie_min << endl;
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used]   // generate list of vectors
//ME20220208 [OBSOLETE - not used]   candidate_lattice_vector.clear();
//ME20220208 [OBSOLETE - not used]   candidate_lattice_vector.push_back(sstr.lattice(1));  // lattice is made of good vectors
//ME20220208 [OBSOLETE - not used]   candidate_lattice_vector.push_back(sstr.lattice(2));  // lattice is made of good vectors
//ME20220208 [OBSOLETE - not used]   candidate_lattice_vector.push_back(sstr.lattice(3));  // lattice is made of good vectors
//ME20220208 [OBSOLETE - not used]   // no threads
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used]   if(LDEBUG) cerr << __AFLOW_FUNC__ << " [3] " << endl; //CO20200201
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used]   if(!AFLOW_PTHREADS::FLAG || specie_min<=specie_min_threshold) {
//ME20220208 [OBSOLETE - not used]     //   cerr << "NO PTHREADS" << endl;
//ME20220208 [OBSOLETE - not used]     for(uint iat1=0;iat1<sstr.atoms.size();iat1++) { //      if(sstr.atoms.at(iat1).type==ispecie_min)
//ME20220208 [OBSOLETE - not used]       for(uint iat2=0;iat2<sstr.atoms.size();iat2++) { //           if(sstr.atoms.at(iat2).type==ispecie_min)
//ME20220208 [OBSOLETE - not used]         if(iat1!=iat2) {
//ME20220208 [OBSOLETE - not used]           fdisp=sstr.atoms.at(iat2).fpos-sstr.atoms.at(iat1).fpos;
//ME20220208 [OBSOLETE - not used]           cdisp=sstr.atoms.at(iat2).cpos-sstr.atoms.at(iat1).cpos;
//ME20220208 [OBSOLETE - not used]           if(aurostd::modulus(fdisp)>0.01 && aurostd::modulus(cdisp)>0.01)
//ME20220208 [OBSOLETE - not used]             if(IsTranslationFVector(sstr,fdisp))
//ME20220208 [OBSOLETE - not used]               candidate_lattice_vector.push_back(cdisp);
//ME20220208 [OBSOLETE - not used]         }
//ME20220208 [OBSOLETE - not used]       }
//ME20220208 [OBSOLETE - not used]     }
//ME20220208 [OBSOLETE - not used]   } else { // multithread
//ME20220208 [OBSOLETE - not used]     cerr << "START THREADS [1] (GetPrimitiveMULTITHREAD=" << AFLOW_PTHREADS::MAX_PTHREADS << ") [" << specie_min << "]" << endl;
//ME20220208 [OBSOLETE - not used]     AFLOW_PTHREADS::Clean_Threads();                                              // multithread clean
//ME20220208 [OBSOLETE - not used]     _threaded_GETTVECTORS_params params[MAX_ALLOCATABLE_PTHREADS];                     // multithread
//ME20220208 [OBSOLETE - not used]     // _threaded_GETTVECTORS_params params[AFLOW_PTHREADS::MAX_PTHREADS];               // multithread
//ME20220208 [OBSOLETE - not used]     // vector<_threaded_GETTVECTORS_params> params(AFLOW_PTHREADS::MAX_PTHREADS);          // multithread
//ME20220208 [OBSOLETE - not used]     // prepare
//ME20220208 [OBSOLETE - not used]     for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) {              // multithread
//ME20220208 [OBSOLETE - not used]       // construction of params[i]                                         // multithread
//ME20220208 [OBSOLETE - not used]       params[ithread].ITHREAD=ithread;                                     // multithread
//ME20220208 [OBSOLETE - not used]       params[ithread].THREADS_MAX=AFLOW_PTHREADS::MAX_PTHREADS;                      // multithread
//ME20220208 [OBSOLETE - not used]       params[ithread].pstr=&sstr;                                          // multithread
//ME20220208 [OBSOLETE - not used]       params[ithread].ispecie_min=ispecie_min;                             // multithread
//ME20220208 [OBSOLETE - not used]       params[ithread].ptvector=&candidate_lattice_vector;                  // multithread
//ME20220208 [OBSOLETE - not used]       params[ithread].polattice=&olattice;                                 // multithread
//ME20220208 [OBSOLETE - not used]       params[ithread].itbusy=ithread;                                      // multithread
//ME20220208 [OBSOLETE - not used]       params[ithread].step=1;  // 1=gettvectors                            // multithread
//ME20220208 [OBSOLETE - not used]     }                                                                      // multithread
//ME20220208 [OBSOLETE - not used]     // run
//ME20220208 [OBSOLETE - not used]     for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) {
//ME20220208 [OBSOLETE - not used]       AFLOW_PTHREADS::viret[ithread]=pthread_create(&(AFLOW_PTHREADS::vpthread[ithread]),NULL,_threaded_GetTvectors,(void*)&params[ithread]);
//ME20220208 [OBSOLETE - not used]       //  aurostd::Sleep(10);
//ME20220208 [OBSOLETE - not used]     }
//ME20220208 [OBSOLETE - not used]     // collect
//ME20220208 [OBSOLETE - not used]     for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++)
//ME20220208 [OBSOLETE - not used]       pthread_join((AFLOW_PTHREADS::vpthread[ithread]),NULL);
//ME20220208 [OBSOLETE - not used]   }
//ME20220208 [OBSOLETE - not used]   if(LDEBUG) cout << "DEBUG: candidate_lattice_vector.size()=" << candidate_lattice_vector.size() << endl;
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used]   // now split
//ME20220208 [OBSOLETE - not used]   int cnt=0;
//ME20220208 [OBSOLETE - not used]   olattice=sstr.lattice;                      // the lattice is always a good lattice
//ME20220208 [OBSOLETE - not used]   if(!AFLOW_PTHREADS::FLAG || specie_min<=specie_min_threshold) {
//ME20220208 [OBSOLETE - not used]     // now generate triplets
//ME20220208 [OBSOLETE - not used]     for(uint iu=0;iu<candidate_lattice_vector.size();iu++) {
//ME20220208 [OBSOLETE - not used]       for(uint i=1;i<=3;i++) plattice[1][i]=candidate_lattice_vector.at(iu)[i];
//ME20220208 [OBSOLETE - not used]       for(uint iv=0;iv<candidate_lattice_vector.size()&& iv!=iu;iv++) {
//ME20220208 [OBSOLETE - not used]         for(uint i=1;i<=3;i++) plattice[2][i]=candidate_lattice_vector.at(iv)[i];
//ME20220208 [OBSOLETE - not used]         for(uint iw=0;iw<candidate_lattice_vector.size()&& iw!=iv && iw!=iu;iw++) {
//ME20220208 [OBSOLETE - not used]           for(uint i=1;i<=3;i++) plattice[3][i]=candidate_lattice_vector.at(iw)[i];
//ME20220208 [OBSOLETE - not used]           if(det(plattice)>0.999 && det(plattice)<sstr_volume) {   // no coplanar and contain at least 1 atom and smaller than the original cell
//ME20220208 [OBSOLETE - not used]             if(aurostd::isinteger(sstr_volume/det(plattice),0.0001)) {    // integer ratio of volumes
//ME20220208 [OBSOLETE - not used]               //  cerr << sstr_volume/det(plattice) << " " << aurostd::isinteger(sstr_volume/det(plattice),0.001) << endl;
//ME20220208 [OBSOLETE - not used]               if(det(plattice)<det(olattice)) {                    // better than before
//ME20220208 [OBSOLETE - not used]                 if(LDEBUG) cout << XPID << "DEBUG"<<iu<<","<<iv<<","<<iw<<" "<< sstr_volume<<" "<<det(plattice)<<" "<<sstr_volume/det(plattice)<<endl;
//ME20220208 [OBSOLETE - not used]                 if(isdifferent(plattice,olattice,0.001)) {
//ME20220208 [OBSOLETE - not used]                   plattice=MinkowskiBasisReduction(plattice);      // Minkowski first
//ME20220208 [OBSOLETE - not used]                   plattice=NiggliUnitCellForm(plattice);           // Niggli Second
//ME20220208 [OBSOLETE - not used]                   if(isdifferent(plattice,olattice,0.001)) {
//ME20220208 [OBSOLETE - not used]                     olattice=plattice;
//ME20220208 [OBSOLETE - not used]                     cnt++;
//ME20220208 [OBSOLETE - not used]                   }
//ME20220208 [OBSOLETE - not used]                 }
//ME20220208 [OBSOLETE - not used]               }
//ME20220208 [OBSOLETE - not used]             }
//ME20220208 [OBSOLETE - not used]           }
//ME20220208 [OBSOLETE - not used]         }
//ME20220208 [OBSOLETE - not used]       }
//ME20220208 [OBSOLETE - not used]     }
//ME20220208 [OBSOLETE - not used]   } else { // multithread
//ME20220208 [OBSOLETE - not used]     //    if(LDEBUG)
//ME20220208 [OBSOLETE - not used]     cerr << "START THREADS [2] (GetPrimitiveMULTITHREAD=" << AFLOW_PTHREADS::MAX_PTHREADS << ") [" << specie_min << "]" << endl;
//ME20220208 [OBSOLETE - not used]     AFLOW_PTHREADS::Clean_Threads();                                              // multithread clean
//ME20220208 [OBSOLETE - not used]     _threaded_GETTVECTORS_params params[MAX_ALLOCATABLE_PTHREADS];                     // multithread
//ME20220208 [OBSOLETE - not used]     // _threaded_GETTVECTORS_params params[AFLOW_PTHREADS::MAX_PTHREADS];               // multithread
//ME20220208 [OBSOLETE - not used]     // vector<_threaded_GETTVECTORS_params> params(AFLOW_PTHREADS::MAX_PTHREADS);          // multithread
//ME20220208 [OBSOLETE - not used]     // prepare
//ME20220208 [OBSOLETE - not used]     for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) {              // multithread
//ME20220208 [OBSOLETE - not used]       // construction of params[i]                                         // multithread
//ME20220208 [OBSOLETE - not used]       params[ithread].ITHREAD=ithread;                                     // multithread
//ME20220208 [OBSOLETE - not used]       params[ithread].THREADS_MAX=AFLOW_PTHREADS::MAX_PTHREADS;                      // multithread
//ME20220208 [OBSOLETE - not used]       params[ithread].pstr=&sstr;                                          // multithread
//ME20220208 [OBSOLETE - not used]       params[ithread].ispecie_min=ispecie_min;                             // multithread
//ME20220208 [OBSOLETE - not used]       params[ithread].ptvector=&candidate_lattice_vector;                  // multithread
//ME20220208 [OBSOLETE - not used]       params[ithread].polattice=&olattice;                                 // multithread
//ME20220208 [OBSOLETE - not used]       params[ithread].itbusy=ithread;                                      // multithread
//ME20220208 [OBSOLETE - not used]       params[ithread].step=2;  // 1=gettvectors                            // multithread
//ME20220208 [OBSOLETE - not used]     }                                                                      // multithread
//ME20220208 [OBSOLETE - not used]     // run
//ME20220208 [OBSOLETE - not used]     for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) {
//ME20220208 [OBSOLETE - not used]       AFLOW_PTHREADS::viret[ithread]=pthread_create(&(AFLOW_PTHREADS::vpthread[ithread]),NULL,_threaded_GetTvectors,(void*)&params[ithread]);
//ME20220208 [OBSOLETE - not used]       //  aurostd::Sleep(10);
//ME20220208 [OBSOLETE - not used]     }
//ME20220208 [OBSOLETE - not used]     // collect
//ME20220208 [OBSOLETE - not used]     for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++)
//ME20220208 [OBSOLETE - not used]       pthread_join((AFLOW_PTHREADS::vpthread[ithread]),NULL);
//ME20220208 [OBSOLETE - not used]     if(LDEBUG) cerr << __AFLOW_FUNC__ << " END THREADS " << endl; //CO20200201
//ME20220208 [OBSOLETE - not used]   }
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used]   plattice=olattice;
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used]   // done
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used]   if(LDEBUG){cerr << __AFLOW_FUNC__ << " sstr=" << endl;cerr << sstr << endl;} //CO20200201
//ME20220208 [OBSOLETE - not used]   xstructure b=sstr;
//ME20220208 [OBSOLETE - not used]   b.lattice=plattice;//b.lattice=roundoff(b.lattice,_EPS_FPOS_EQUAL_);
//ME20220208 [OBSOLETE - not used]   b.FixLattices();
//ME20220208 [OBSOLETE - not used]   b.write_lattice_flag=FALSE;
//ME20220208 [OBSOLETE - not used]   b.write_klattice_flag=FALSE;
//ME20220208 [OBSOLETE - not used]   b.write_DEBUG_flag=FALSE;
//ME20220208 [OBSOLETE - not used]   // plug them all
//ME20220208 [OBSOLETE - not used]   for(uint iat=0;iat<sstr.atoms.size();iat++) {
//ME20220208 [OBSOLETE - not used]     b.atoms.at(iat).fpos=BringInCell(C2F(b.lattice,b.atoms.at(iat).cpos));
//ME20220208 [OBSOLETE - not used]     b.atoms.at(iat).cpos=F2C(b.lattice,b.atoms.at(iat).fpos);
//ME20220208 [OBSOLETE - not used]   }
//ME20220208 [OBSOLETE - not used]   // now remove them
//ME20220208 [OBSOLETE - not used]   b=BringInCell(b);b.RemoveFractionalCopies(tolerance);
//ME20220208 [OBSOLETE - not used]   b=BringInCell(b);b.RemoveCartesianCopies(0.01);
//ME20220208 [OBSOLETE - not used]   // rescale back to original scale.
//ME20220208 [OBSOLETE - not used]   b.SetVolume(Volume(a)*b.atoms.size()/a.atoms.size());
//ME20220208 [OBSOLETE - not used]   b=ReScale(b,a.scale);
//ME20220208 [OBSOLETE - not used]   //  // fix it up with the new Minkowsky and Niggli reductions // CANT DO AUTOMATICALLY
//ME20220208 [OBSOLETE - not used]   //  b=LatticeReduction(b);
//ME20220208 [OBSOLETE - not used]   // Put everything in new primitive cell.
//ME20220208 [OBSOLETE - not used]   b=BringInCell(b);
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used]   // check !
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used]   // no fractional ratio of atoms
//ME20220208 [OBSOLETE - not used]   double fraction_atoms=(double) a.atoms.size()/b.atoms.size();
//ME20220208 [OBSOLETE - not used]   if(!aurostd::isinteger(fraction_atoms,0.01)) return _a;
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used]   // nearest too close
//ME20220208 [OBSOLETE - not used]   if(NearestNeighbor(b)<5.0*tolerance)  return _a;
//ME20220208 [OBSOLETE - not used]   if(NearestNeighbor(b)<0.1)  return _a;
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used]   // no messed up volume
//ME20220208 [OBSOLETE - not used]   double fraction=Volume(a)/Volume(b);
//ME20220208 [OBSOLETE - not used]   if(abs(b.atoms.size()*fraction-a.atoms.size())>0.1) {
//ME20220208 [OBSOLETE - not used]     stringstream message;
//ME20220208 [OBSOLETE - not used]     message << "ERROR   xstructure xstructure::GetPrimitive(void)" << endl;
//ME20220208 [OBSOLETE - not used]     message << "        supercell has the wrong number of atoms" << endl;
//ME20220208 [OBSOLETE - not used]     message << "        volume original    = " << Volume(a) << endl;
//ME20220208 [OBSOLETE - not used]     message << "        volume prim        = " << Volume(b) << endl;
//ME20220208 [OBSOLETE - not used]     message << "        a.scale            = " << a.scale << endl;
//ME20220208 [OBSOLETE - not used]     message << "        b.scale            = " << b.scale << endl;
//ME20220208 [OBSOLETE - not used]     message << "        a.atoms.size()     = " << a.atoms.size() << endl;
//ME20220208 [OBSOLETE - not used]     message << "        b.atoms.size()     = " << b.atoms.size() << endl;
//ME20220208 [OBSOLETE - not used]     message << "        fraction           = " << fraction << endl;
//ME20220208 [OBSOLETE - not used]     message << "        supercell atoms    = " << fraction*b.atoms.size() << endl;
//ME20220208 [OBSOLETE - not used]     message << b << endl;
//ME20220208 [OBSOLETE - not used]     throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::GetPrimitive(void)",message,_RUNTIME_ERROR_);
//ME20220208 [OBSOLETE - not used]   }
//ME20220208 [OBSOLETE - not used]   // everything ok
//ME20220208 [OBSOLETE - not used]   if(LDEBUG) cerr << __AFLOW_FUNC__ << " END [ok]=" << fraction_atoms << endl;  //CO20200201
//ME20220208 [OBSOLETE - not used]   b.ClearSymmetry();  //CO20181226 - new structure, symmetry not calculated
//ME20220208 [OBSOLETE - not used]   b.primitive_calculated = TRUE; //DX20201007
//ME20220208 [OBSOLETE - not used]   return b;
//ME20220208 [OBSOLETE - not used]
//ME20220208 [OBSOLETE - not used] }


xstructure GetPrimitiveSINGLE(const xstructure& _a,double tolerance) {  // APRIL 2009JUNE 2012 added tolerance
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  cout.setf(std::ios::fixed,std::ios::floatfield);
  cout.precision(10);
  xstructure a(_a);
  xstructure sstr=a;
  sstr.SetVolume(sstr.atoms.size());
  sstr=ReScale(sstr,1.0);
  sstr=BringInCell(sstr);
  if(tolerance<=0.0) a.equiv_fpos_epsilon=_EQUIV_FPOS_EPS_; else a.equiv_fpos_epsilon=tolerance;
  double sstr_volume=sstr.Volume();

  _aflags aflags;

  xmatrix<double> plattice(3,3),olattice(3,3);
  xvector<double> fdisp(3),cdisp(3);
  std::vector<xvector<double> > candidate_lattice_vector;

  int specie_min=100000,ispecie_min=0;
  for(uint ispecie=0;ispecie<sstr.num_each_type.size();ispecie++)
    if(sstr.num_each_type.at(ispecie)<specie_min) {
      specie_min=sstr.num_each_type.at(ispecie);
      ispecie_min=ispecie;
    }

  // generate list of vectors
  candidate_lattice_vector.push_back(sstr.lattice(1));  // lattice is made of good vectors
  candidate_lattice_vector.push_back(sstr.lattice(2));  // lattice is made of good vectors
  candidate_lattice_vector.push_back(sstr.lattice(3));  // lattice is made of good vectors
  for(uint iat1=0;iat1<sstr.atoms.size();iat1++)
    if(sstr.atoms.at(iat1).type==ispecie_min)
      for(uint iat2=0;iat2<sstr.atoms.size();iat2++)
        if(sstr.atoms.at(iat2).type==ispecie_min) {
          fdisp=sstr.atoms.at(iat2).fpos-sstr.atoms.at(iat1).fpos;
          cdisp=sstr.atoms.at(iat2).cpos-sstr.atoms.at(iat1).cpos;
          if(aurostd::modulus(fdisp)>0.01 && aurostd::modulus(cdisp)>0.01)
            if(IsTranslationFVector(sstr,fdisp))
              candidate_lattice_vector.push_back(cdisp);
        }
  if(LDEBUG) cout << "DEBUG" << candidate_lattice_vector.size() << endl;
  int cnt=0;
  olattice=sstr.lattice;                      // the lattice is always a good lattice
  // now generate triplets
  for(uint iu=0;iu<candidate_lattice_vector.size();iu++) {
    for(uint i=1;i<=3;i++) plattice(1,i)=candidate_lattice_vector.at(iu)(i);
    for(uint iv=0;iv<candidate_lattice_vector.size()&& iv!=iu;iv++) {
      for(uint i=1;i<=3;i++) plattice(2,i)=candidate_lattice_vector.at(iv)(i);
      for(uint iw=0;iw<candidate_lattice_vector.size()&& iw!=iv && iw!=iu;iw++) {
        for(uint i=1;i<=3;i++) plattice(3,i)=candidate_lattice_vector.at(iw)(i);
        if(det(plattice)>0.999 && det(plattice)<sstr_volume) {   // no coplanar and contain at least 1 atom and smaller than the original cell
          if(aurostd::isinteger(sstr_volume/det(plattice))) {    // integer ratio of volumes
            if(det(plattice)<det(olattice)) {                    // better than before
              if(LDEBUG) cout << XPID << "DEBUG"<<iu<<","<<iv<<","<<iw<<" "<< sstr_volume<<" "<<det(plattice)<<" "<<sstr_volume/det(plattice)<<endl;
              if(isdifferent(plattice,olattice,0.0001)) {
                plattice=MinkowskiBasisReduction(plattice);      // Minkowski first
                plattice=NiggliUnitCellForm(plattice);           // Niggli Second
                if(isdifferent(plattice,olattice,0.0001)) {
                  olattice=plattice;
                  cnt++;
                }
              }
            }
          }
        }
      }
    }
  }
  plattice=olattice;
  // done

  xstructure b=sstr;
  b.lattice=plattice;//b.lattice=roundoff(b.lattice,_EPS_FPOS_EQUAL_);
  b.FixLattices();
  b.write_lattice_flag=FALSE;
  b.write_klattice_flag=FALSE;
  b.write_DEBUG_flag=FALSE;
  // plug them all
  for(uint iat=0;iat<sstr.atoms.size();iat++) {
    b.atoms.at(iat).fpos=BringInCell(C2F(b.lattice,b.atoms.at(iat).cpos));
    b.atoms.at(iat).cpos=F2C(b.lattice,b.atoms.at(iat).fpos);
  }
  // now remove them
  b.RemoveFractionalCopies();
  // rescale back to original scale.
  b.SetVolume(Volume(a)*b.atoms.size()/a.atoms.size());
  b=ReScale(b,a.scale);
  //  // fix it up with the new Minkowsky and Niggli reductions // CANT DO AUTOMATICALLY
  //  b=LatticeReduction(b);
  // Put everything in new primitive cell.
  b=BringInCell(b);
  // check !
  double fraction=Volume(a)/Volume(b);
  if(abs(b.atoms.size()*fraction-a.atoms.size())>0.1) {
    stringstream message;
    message << "ERROR   xstructure xstructure::GetPrimitive(void)" << endl;
    message << "        supercell has the wrong number of atoms" << endl;
    message << "        volume original    = " << Volume(a) << endl;
    message << "        volume prim        = " << Volume(b) << endl;
    message << "        a.scale            = " << a.scale << endl;
    message << "        b.scale            = " << b.scale << endl;
    message << "        a.atoms.size()     = " << a.atoms.size() << endl;
    message << "        b.atoms.size()     = " << b.atoms.size() << endl;
    message << "        fraction           = " << fraction << endl;
    message << "        supercell atoms    = " << fraction*b.atoms.size() << endl;
    message << b << endl;
    throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::GetPrimitiveSINGLE()",message,_RUNTIME_ERROR_);
  }
  // everything ok
  b.ClearSymmetry();  //CO20181226 - new structure, symmetry not calculated
  b.primitive_calculated = TRUE; //DX20201007
  return b;
}


xstructure GetPrimitive1(const xstructure& a) {  // MARCH 2009
  stringstream message;
  if(a.equiv_fpos_epsilon<1.0e-12) {
    message << "Zero tolerance: " << a.equiv_fpos_epsilon;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
  }
  double tolerance=a.equiv_fpos_epsilon;

  cout.setf(std::ios::fixed,std::ios::floatfield);
  cout.precision(10);
  xstructure sstr=a;
  sstr.SetVolume(sstr.atoms.size());
  sstr=ReScale(sstr,1.0);
  sstr=BringInCell(sstr);
  double sstr_volume=sstr.Volume();
  _aflags aflags;
  // identify the minimum set of atoms

  xmatrix<double> plattice(3,3),olattice(3,3);
  xvector<double> fdisp(3),cdisp(3);
  std::vector<xvector<double> > candidate_lattice_vector;

  int specie_min=100000,ispecie_min=0;
  for(uint ispecie=0;ispecie<sstr.num_each_type.size();ispecie++)
    if(sstr.num_each_type.at(ispecie)<specie_min) {
      specie_min=sstr.num_each_type.at(ispecie);
      ispecie_min=ispecie;
    }

  // generate list of vectors
  candidate_lattice_vector.push_back(sstr.lattice(1));  // lattice is made of good vectors
  candidate_lattice_vector.push_back(sstr.lattice(2));  // lattice is made of good vectors
  candidate_lattice_vector.push_back(sstr.lattice(3));  // lattice is made of good vectors
  olattice=sstr.lattice;                      // the lattice is always a good lattice
  for(uint iat1=0;iat1<sstr.atoms.size();iat1++) {
    for(uint iat2=0;iat2<sstr.atoms.size();iat2++) {
      if(iat2!=iat1 &&
          sstr.atoms.at(iat1).type==ispecie_min &&
          sstr.atoms.at(iat2).type==ispecie_min) {
        fdisp=sstr.atoms.at(iat2).fpos-sstr.atoms.at(iat1).fpos;
        cdisp=sstr.atoms.at(iat2).cpos-sstr.atoms.at(iat1).cpos;
        if(IsTranslationFVector(sstr,fdisp)) candidate_lattice_vector.push_back(cdisp);
      }
    }
  }
  //  cerr << candidate_lattice_vector.size() << endl;
  int cnt=0;
  // now generate triplets
  for(uint iu=0;iu<candidate_lattice_vector.size();iu++) {
    for(uint i=1;i<=3;i++) plattice(1,i)=candidate_lattice_vector.at(iu)(i);
    for(uint iv=0;iv<candidate_lattice_vector.size()&& iv!=iu;iv++) {
      for(uint i=1;i<=3;i++) plattice(2,i)=candidate_lattice_vector.at(iv)(i);
      for(uint iw=0;iw<candidate_lattice_vector.size()&& iw!=iv && iw!=iu;iw++) {
        for(uint i=1;i<=3;i++) plattice(3,i)=candidate_lattice_vector.at(iw)(i);
        if(det(plattice)>tolerance && det(plattice)<sstr_volume && det(plattice)<det(olattice)) { // well defined
          if(isdifferent(plattice,olattice,0.0001)) {
            plattice=MinkowskiBasisReduction(plattice); // Minkowski first
            plattice=NiggliUnitCellForm(plattice);      // Niggli Second
            if(isdifferent(plattice,olattice,0.0001)) {
              olattice=plattice;
              //      cerr << det(olattice) << " " << cnt<< endl;
              cnt++;
            }
          }
        }
      }
    }
  }
  plattice=olattice;
  // done

  _atom atom;
  xstructure b=sstr;
  b.atoms.clear();
  b.lattice=plattice;b.lattice=roundoff(b.lattice,tolerance);
  b.FixLattices();
  b.write_lattice_flag=FALSE;
  b.write_klattice_flag=FALSE;
  b.write_DEBUG_flag=FALSE;
  bool atom_found=FALSE;
  // for scanning around
  //  double radius=aurostd::modulus(sstr.lattice(1))+aurostd::modulus(sstr.lattice(2))+aurostd::modulus(sstr.lattice(3));
  //  int dims=max(LatticeDimensionSphere(sstr.lattice,radius));
  int dims=max(LatticeDimensionSphere(sstr.lattice,1.5*RadiusSphereLattice(sstr.lattice)));
  for(uint iat=0;iat<b.num_each_type.size();iat++) {
    b.num_each_type.at(iat)=0; // create enough space
    b.comp_each_type.at(iat)=0; // create enough space
  }

  for(uint iat=0;iat<sstr.atoms.size();iat++) {
    for(int i=-dims;i<=dims;i++)
      for(int j=-dims;j<=dims;j++)
        for(int k=-dims;k<=dims;k++) {
          //    atom=BringInCell(sstr.atoms.at(i),sstr.lattice);
          atom=sstr.atoms.at(iat);
          atom.cpos=atom.cpos+((double)i)*sstr.lattice(1)+((double)j)*sstr.lattice(2)+((double)k)*sstr.lattice(3);
          atom.fpos=C2F(b.lattice,atom.cpos);
          if(atom.fpos(1)>=-tolerance && atom.fpos(1)<1.0-tolerance &&
              atom.fpos(2)>=-tolerance && atom.fpos(2)<1.0-tolerance &&
              atom.fpos(3)>=-tolerance && atom.fpos(3)<1.0-tolerance) {      // found something inside
            for(uint ii=0;ii<b.atoms.size()&&!atom_found;ii++)
              atom_found=identical(atom.cpos,b.atoms.at(ii).cpos,0.1);       // look in all the list of operations
            // atom_found=FALSE;
            if(!atom_found) {
              atom.fpos=roundoff(atom.fpos,tolerance);
              atom.cpos=roundoff(atom.cpos,tolerance);
              b.atoms.push_back(atom);
              b.num_each_type.at(atom.type)++;         // CONVASP_MODE
              b.comp_each_type.at(atom.type)+=atom.partial_occupation_value;         // CONVASP_MODE
            }
          }
        }
  }
  b.GetStoich();  //CO20170724
  // rescale back to original scale.
  b.SetVolume(Volume(a)*b.atoms.size()/a.atoms.size());
  b=ReScale(b,a.scale);
  // fix it up with the new Minkowsky and Niggli reductions
  b=LatticeReduction(b);
  // Put everything in new primitive cell.
  b=BringInCell(b);
  // check !
  double fraction=Volume(a)/Volume(b);
  if(abs(b.atoms.size()*fraction-a.atoms.size())>0.1) {
    message << "ERROR   xstructure xstructure::GetPrimitive(void)" << endl;
    message << "        supercell has the wrong number of atoms" << endl;
    message << "        volume original    = " << Volume(a) << endl;
    message << "        volume prim        = " << Volume(b) << endl;
    message << "        a.scale            = " << a.scale << endl;
    message << "        b.scale            = " << b.scale << endl;
    message << "        a.atoms.size()     = " << a.atoms.size() << endl;
    message << "        b.atoms.size()     = " << b.atoms.size() << endl;
    message << "        fraction           = " << fraction << endl;
    message << "        supercell atoms    = " << fraction*b.atoms.size() << endl;
    message << b << endl;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
  }
  // everything ok
  b.primitive_calculated = TRUE; //DX20201007
  return b;
}



// second try
xstructure GetPrimitive2(const xstructure& a) {
  stringstream message;
  if(a.equiv_fpos_epsilon<1.0e-12) {
    message << "Zero tolerance: " << a.equiv_fpos_epsilon;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
  }
  double tolerance=a.equiv_fpos_epsilon;

  cout.setf(std::ios::fixed,std::ios::floatfield);
  cout.precision(10);

  xstructure sstr=a;
  //  sstr.lattice=lattice;
  // sstr.scale=scale;
  std::vector<xvector<double> > candidate_lattice_vector;
  xmatrix<double> plattice(3,3);
  vector<xmatrix<double> > plattice_list;
  xvector<double> fdisp(3),cdisp(3);

  // Get all the data from the structure for easy use.
  // Set scale to 1 so you don't need to rescale coordinates.
  sstr=ReScale(sstr,1.0);
  // Put everything in the unit cell.
  //  cerr << sstr.scale << endl;
  sstr=BringInCell(sstr);
  // cerr << sstr.scale << endl;
  // throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_);
  string title=sstr.title;
  int i;
  for(uint iat=0;iat<sstr.num_each_type.size();iat++) {
    sstr.num_each_type.at(iat)=0; // create enough space
    sstr.comp_each_type.at(iat)=0; // create enough space
  }
  // Create list of candidate plvs.

  candidate_lattice_vector.push_back(sstr.lattice(1));
  candidate_lattice_vector.push_back(sstr.lattice(2));
  candidate_lattice_vector.push_back(sstr.lattice(3));

  for(uint ia=1;ia<sstr.atoms.size();ia++) {
    // Here we calculate all displacements from the first atom, i.e., the
    // first atom is considered the origin.
    cdisp=sstr.atoms.at(ia).cpos-sstr.atoms.at(0).cpos;
    if(IsTranslationCVector(sstr,cdisp)) candidate_lattice_vector.push_back(cdisp);
  }
  // Now the candidate_lattice_vector have been found and we must take all possible
  // traids of *distinct* vectors and get the smallest volume.
  int num_cand=candidate_lattice_vector.size();
  double min_vol,new_vol;
  min_vol=GetVol(sstr.lattice);
  for(int iu=0;iu<num_cand;iu++) {
    for(i=1;i<=3;i++) plattice(1,i)=candidate_lattice_vector.at(iu)(i);
    for(int iv=0;iv<num_cand;iv++) {
      for(i=1;i<=3;i++) plattice(2,i)=candidate_lattice_vector.at(iv)(i);
      for(int iw=0;iw<num_cand;iw++) {
        if(iu!=iw && iu!=iv && iv!=iw) {
          for(i=1;i<=3;i++) plattice(3,i)=candidate_lattice_vector.at(iw)(i);
          new_vol=GetVol(plattice);
          // if(new_vol<=min_vol && new_vol>tolerance) min_vol=new_vol;
          if(new_vol<=min_vol && new_vol>1.0e-5) min_vol=new_vol;
        }
      } //iw loop
    } //iv loop
  } //iu loop
  // Now that we know the min_vol we must go back through all the triads
  // and store the set that have the min_vol.
  for(int iu=0;iu<num_cand;iu++) {
    for(i=1;i<=3;i++) plattice(1,i)=candidate_lattice_vector.at(iu)(i);
    for(int iv=0;iv<num_cand;iv++) {
      for(i=1;i<=3;i++) plattice(2,i)=candidate_lattice_vector.at(iv)(i);
      for(int iw=0;iw<num_cand;iw++) {
        if(iu!=iw && iu!=iv && iv!=iw) {
          for(i=1;i<=3;i++) plattice(3,i)=candidate_lattice_vector.at(iw)(i);
          new_vol=GetVol(plattice);
          if(aurostd::abs(new_vol-min_vol)<tolerance && new_vol>1.0e-5)
            plattice_list.push_back(plattice); // Add to plattice_list
          //if new_vol=min_vol
        } //if iu!=iv!=iw
      } //iw loop
    } //iv loop
  } //iu loop
  // Now we have a set of all possible primitive lattices
  // and we want to pick one that looks as much like the
  // original lattice as possible.
  double paw,pawopt=1e6;
  for(uint i=0;i<plattice_list.size();i++) {
    paw=0;
    for(int j=1;j<=3;j++) paw+=angle(plattice_list.at(i)(j),sstr.lattice(j));
    if(paw<pawopt) {pawopt=paw;plattice=plattice_list.at(i);}; // minimization
  }
  // For now just pick the first one.
  // if(plattice_list.size()>0) plattice=plattice_list.at(1);
  // If volume did not reduce then just keep original lattice
  if(aurostd::abs(GetVol(sstr.lattice)-GetVol(plattice))<tolerance) {plattice=sstr.lattice;}

  // Now we have the plvs for the new structure.  Now we must
  // construct a new structure object with the plvs and return it.
  // This new structure object will be just like the old one but
  // with new lattice vectors and new basis atoms (which will be a
  // subset of the original ones, chosen to have direct coordinates
  // less than 1).

  // To create new structure we must set all the parameters.
  // To get new basis atoms, loop over all the original basis atoms,
  // get new direct coordinates, make all less than one, then take unique
  // values.

  _atom atom;
  xstructure b=sstr;
  b.atoms.clear();
  b.lattice=plattice;b.lattice=roundoff(b.lattice,tolerance);
  b.FixLattices();
  b.write_lattice_flag=TRUE;
  b.write_klattice_flag=FALSE;
  b.write_DEBUG_flag=FALSE;
  bool atom_found=FALSE;
  // for scanning around
  double radius;
  radius=aurostd::modulus(sstr.lattice(1))+aurostd::modulus(sstr.lattice(2))+aurostd::modulus(sstr.lattice(3));
  int dims=max(LatticeDimensionSphere(sstr.lattice,radius));
  for(uint iat=0;iat<b.num_each_type.size();iat++) {
    b.num_each_type.at(iat)=0; // create enough space
    b.comp_each_type.at(iat)=0; // create enough space
  }

  for(uint iat=0;iat<sstr.atoms.size();iat++) {
    for(int i=-dims;i<=dims;i++)
      for(int j=-dims;j<=dims;j++)
        for(int k=-dims;k<=dims;k++) {
          //    atom=BringInCell(sstr.atoms.at(i),sstr.lattice);
          atom=sstr.atoms.at(iat);
          atom.cpos=atom.cpos+((double)i)*sstr.lattice(1)+((double)j)*sstr.lattice(2)+((double)k)*sstr.lattice(3);
          atom.fpos=C2F(b.lattice,atom.cpos);
          if(atom.fpos(1)>=-tolerance && atom.fpos(1)<1.0-tolerance &&
              atom.fpos(2)>=-tolerance && atom.fpos(2)<1.0-tolerance &&
              atom.fpos(3)>=-tolerance && atom.fpos(3)<1.0-tolerance) {      // found something inside
            for(uint ii=0;ii<b.atoms.size()&&!atom_found;ii++)
              atom_found=identical(atom.cpos,b.atoms.at(ii).cpos,0.1);       // look in all the list of operations
            // atom_found=FALSE;
            if(!atom_found) {
              atom.fpos=roundoff(atom.fpos,tolerance);
              atom.cpos=roundoff(atom.cpos,tolerance);
              b.atoms.push_back(atom);
              b.num_each_type.at(atom.type)++;         // CONVASP_MODE
              b.comp_each_type.at(atom.type)+=atom.partial_occupation_value;         // CONVASP_MODE
            }
          }
        }
  }
  b.GetStoich();  //CO20170724
  // some check
  double fraction=((double) GetVol(sstr.lattice)/GetVol(b.lattice));
  if(abs(b.atoms.size()*fraction-a.atoms.size())>0.1) {
    cerr << "ERROR   xstructure xstructure::GetPrimitive(void)" << endl;
    cerr << "        supercell has the wrong number of atoms" << endl;
    cerr << "        volume original    = " << GetVol(sstr.lattice) << endl;
    cerr << "        volume prim        = " << GetVol(b.lattice) << endl;
    //cerr << "        sstr.scale            = " << sstr.scale << endl;
    //cerr << "        b.scale            = " << b.scale << endl;
    cerr << "        sstr.atoms.size()     = " << a.atoms.size() << endl;
    cerr << "        b.atoms.size()     = " << b.atoms.size() << endl;
    cerr << "        fraction           = " << fraction << endl;
    cerr << "        supercell atoms    = " << fraction*b.atoms.size() << endl;
    cerr << GetVol(b.lattice) << endl;
    cerr << GetVol(a.lattice) << endl;
    //   throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::GetPrimitive():","Throw for debugging purposes.",_GENERIC_ERROR_);
  }
  // rescale back to original scale.
  b=ReScale(b,a.scale);
  // fix it up with the new Minkowsky and Niggli reductions
  b=LatticeReduction(b);
  // Put everything in new primitive cell.
  b=BringInCell(b);
  // everything ok
  b.primitive_calculated = TRUE; //DX20201007
  return b;
}

// third try
xstructure GetPrimitive3(const xstructure& a) {
  stringstream message;
  if(a.equiv_fpos_epsilon<1.0e-12) {
    message << "Zero tolerance: " << a.equiv_fpos_epsilon;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
  }
  double tolerance=a.equiv_fpos_epsilon;

  cout.setf(std::ios::fixed,std::ios::floatfield);
  cout.precision(10);
  xstructure sstr=a;
  sstr.SetVolume(sstr.atoms.size());
  sstr=ReScale(sstr,1.0);
  sstr=BringInCell(sstr);

  double sstr_volume=sstr.Volume();
  _aflags aflags;aflags.Directory="./";
  // identify the minimum set of atoms
  bool PGROUPWRITE=TRUE,FGROUPWRITE=TRUE;//,IATOMSWRITE=TRUE;
  bool OSSWRITE=FALSE; // to FileMESSAGE, does not matter as it is /dev/null
  ofstream FileMESSAGE("/dev/null");
  bool _QUIET_=XHOST.QUIET;XHOST.QUIET=TRUE;
  SYM::CalculatePointGroup(FileMESSAGE,sstr,aflags,PGROUPWRITE,OSSWRITE,cout);
  SYM::CalculateFactorGroup(FileMESSAGE,sstr,aflags,FGROUPWRITE,OSSWRITE,cout);
  //  SYM::CalculateInequivalentAtoms(FileMESSAGE,sstr,aflags,IATOMSWRITE,OSSWRITE,cout);

  xmatrix<double> plattice(3,3),olattice(3,3);
  xvector<double> fdisp(3),cdisp(3);
  std::vector<xvector<double> > candidate_lattice_vector;

  // generate list of vectors
  candidate_lattice_vector.push_back(sstr.lattice(1));  // lattice is made of good vectors
  candidate_lattice_vector.push_back(sstr.lattice(2));  // lattice is made of good vectors
  candidate_lattice_vector.push_back(sstr.lattice(3));  // lattice is made of good vectors
  olattice=sstr.lattice;                      // the lattice is always a good lattice

  bool sym_found;
  for(uint i=0;i<sstr.fgroup.size();i++) {
    if(aurostd::modulus(sstr.fgroup.at(i).ctau)>0.01) {
      //    cerr << i << " " << sstr.fgroup.at(i).ctau << endl;
      sym_found=FALSE;
      for(uint ii=0;ii<candidate_lattice_vector.size()&&!sym_found;ii++)
        sym_found=identical(sstr.fgroup.at(i).ctau,candidate_lattice_vector[ii],tolerance);    // look in all the list of operations
      if(sym_found==FALSE) {                                          // new operation, generate and save it
        candidate_lattice_vector.push_back(sstr.fgroup.at(i).ctau);
        cerr << i << " " << sstr.fgroup.at(i).ctau << endl;
      }
    }
  }

  cerr << candidate_lattice_vector.size() << endl;
  int cnt=0;
  // now generate triplets
  for(uint iu=0;iu<candidate_lattice_vector.size();iu++) {
    for(uint i=1;i<=3;i++) plattice(1,i)=candidate_lattice_vector.at(iu)(i);
    for(uint iv=0;iv<candidate_lattice_vector.size()&& iv!=iu;iv++) {
      for(uint i=1;i<=3;i++) plattice(2,i)=candidate_lattice_vector.at(iv)(i);
      for(uint iw=0;iw<candidate_lattice_vector.size()&& iw!=iv && iw!=iu;iw++) {
        for(uint i=1;i<=3;i++) plattice(3,i)=candidate_lattice_vector.at(iw)(i);
        if(det(plattice)>tolerance && det(plattice)<sstr_volume && det(plattice)<det(olattice)) { // well defined
          if(isdifferent(plattice,olattice,0.0001)) {
            plattice=MinkowskiBasisReduction(plattice); // Minkowski first
            plattice=NiggliUnitCellForm(plattice);      // Niggli Second
            if(isdifferent(plattice,olattice,0.0001)) {
              olattice=plattice;
              //      cerr << det(olattice) << " " << cnt<< endl;
              cnt++;
            }
          }
        }
      }
    }
  }
  plattice=olattice;
  // done

  _atom atom;
  xstructure b=sstr;
  b.atoms.clear();
  b.lattice=plattice;b.lattice=roundoff(b.lattice,tolerance);
  b.FixLattices();
  b.write_lattice_flag=FALSE;
  b.write_klattice_flag=FALSE;
  b.write_DEBUG_flag=FALSE;
  bool atom_found=FALSE;
  // for scanning around
  int dims=max(LatticeDimensionSphere(sstr.lattice,RadiusSphereLattice(sstr.lattice)));
  for(uint iat=0;iat<b.num_each_type.size();iat++) {
    b.num_each_type.at(iat)=0; // create enough space
    b.comp_each_type.at(iat)=0; // create enough space
  }

  for(uint iat=0;iat<sstr.atoms.size();iat++) {
    for(int i=-dims;i<=dims;i++)
      for(int j=-dims;j<=dims;j++)
        for(int k=-dims;k<=dims;k++) {
          //    atom=BringInCell(sstr.atoms.at(i),sstr.lattice);
          atom=sstr.atoms.at(iat);
          atom.cpos=atom.cpos+((double)i)*sstr.lattice(1)+((double)j)*sstr.lattice(2)+((double)k)*sstr.lattice(3);
          atom.fpos=C2F(b.lattice,atom.cpos);
          if(atom.fpos(1)>=-tolerance && atom.fpos(1)<1.0-tolerance &&
              atom.fpos(2)>=-tolerance && atom.fpos(2)<1.0-tolerance &&
              atom.fpos(3)>=-tolerance && atom.fpos(3)<1.0-tolerance) {      // found something inside
            for(uint ii=0;ii<b.atoms.size()&&!atom_found;ii++)
              atom_found=identical(atom.cpos,b.atoms.at(ii).cpos,0.1);       // look in all the list of operations
            // atom_found=FALSE;
            if(!atom_found) {
              atom.fpos=roundoff(atom.fpos,tolerance);
              atom.cpos=roundoff(atom.cpos,tolerance);
              b.atoms.push_back(atom);
              b.num_each_type.at(atom.type)++;         // CONVASP_MODE
              b.comp_each_type.at(atom.type)+=atom.partial_occupation_value;         // CONVASP_MODE
            }
          }
        }
  }
  b.GetStoich();  //CO20170724
  // rescale back to original scale.
  b.SetVolume(Volume(a)*b.atoms.size()/a.atoms.size());
  b=ReScale(b,a.scale);
  // fix it up with the new Minkowsky and Niggli reductions
  b=LatticeReduction(b);
  XHOST.QUIET=_QUIET_;

  // Put everything in new primitive cell.
  b=BringInCell(b);
  // check !
  double fraction=Volume(a)/Volume(b);
  if(abs(b.atoms.size()*fraction-a.atoms.size())>0.1) {
    message << "ERROR   xstructure xstructure::GetPrimitive(void)" << endl;
    message << "        supercell has the wrong number of atoms" << endl;
    message << "        volume original    = " << Volume(a) << endl;
    message << "        volume prim        = " << Volume(b) << endl;
    message << "        a.scale            = " << a.scale << endl;
    message << "        b.scale            = " << b.scale << endl;
    message << "        a.atoms.size()     = " << a.atoms.size() << endl;
    message << "        b.atoms.size()     = " << b.atoms.size() << endl;
    message << "        fraction           = " << fraction << endl;
    message << "        supercell atoms    = " << fraction*b.atoms.size() << endl;
    message << b << endl;
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
  }
  // everything ok
  b.primitive_calculated = TRUE; //DX20201007
  return b;
}




// ***************************************************************************
// Operator GetPrimitive
// ***************************************************************************

void xstructure::GetPrimitive(void) {
  extern xstructure GetPrimitive(const xstructure& a); // so it does not recurse
  xstructure a(*this);
  *this=GetPrimitive(a);
}

void xstructure::GetPrimitive(double tolerance) {
  extern xstructure GetPrimitive(const xstructure& a,double tolerance); // so it does not recurse
  xstructure a(*this);
  *this=GetPrimitive(a,tolerance);
}

void xstructure::GetPrimitive2(void) {
  extern xstructure GetPrimitive2(const xstructure& a); // so it does not recurse
  xstructure a(*this);
  *this=GetPrimitive2(a);
}

void xstructure::GetPrimitive3(void) {
  extern xstructure GetPrimitive3(const xstructure& a); // so it does not recurse
  xstructure a(*this);
  *this=GetPrimitive3(a);
}

// ***************************************************************************
// Function GetPrimitiveStructures() //DX20201006
// ***************************************************************************
void GetPrimitiveStructures(vector<xstructure>& structures, uint start_index, uint end_index){

  // Converts a set of xstructures to their primitive representation
  // Optional indices can be included; useful for pre-distributed
  // threading schemes
  // Default: run over entire range

  // if end index is greater than structures.size(), then compute primitive cell for all structures
  if(end_index > structures.size()){ end_index=structures.size(); }

  for(uint i=start_index;i<end_index;i++){
    structures[i].GetPrimitive();
  }
}

// ***************************************************************************
// Function MinDist
// ***************************************************************************
double xstructure::MinDist(void) {
  dist_nn_min=SYM::minimumDistance(*this);
  return dist_nn_min;
}

// ***************************************************************************
// Function NearestNeighbor() // moved from aflow_xproto.cpp
// ***************************************************************************
double NearestNeighbor(const xstructure &str_in) {
  return SYM::minimumDistance(str_in);
  //[CO20171024 OBSOLETE]xstructure str(str_in);
  //[CO20171024 OBSOLETE]// if(LDEBUG) { cerr << "NearestNeighbor 1" << endl; }
  //[CO20171024 OBSOLETE]// if(LDEBUG) { cerr << str.scale << endl; }
  //[CO20171024 OBSOLETE]str.ReScale(1.0);
  //[CO20171024 OBSOLETE]xvector<int> ndims(3);
  //[CO20171024 OBSOLETE]// str.neighbors_radius=RadiusSphereLattice(str.lattice);
  //[CO20171024 OBSOLETE]// str.neighbors_radius=max(modulus(str.lattice(1)),modulus(str.lattice(2)),modulus(str.lattice(3)));
  //[CO20171024 OBSOLETE]// ndims=LatticeDimensionSphere(str.lattice,str.neighbors_radius);
  //[CO20171024 OBSOLETE]ndims[1]=ndims[2]=ndims[3]=1;
  //[CO20171024 OBSOLETE]deque<_atom> vatoms;
  //[CO20171024 OBSOLETE]_atom atom;
  //[CO20171024 OBSOLETE]// if(LDEBUG) { cerr << "NearestNeighbor 2" << endl; }
  //[CO20171024 OBSOLETE]for(int i=-ndims[1];i<=ndims[1];i++) {
  //[CO20171024 OBSOLETE]  for(int j=-ndims[2];j<=ndims[2];j++) {
  //[CO20171024 OBSOLETE]    for(int k=-ndims[3];k<=ndims[3];k++) {
  //[CO20171024 OBSOLETE]for(uint iat=0;iat<str.atoms.size();iat++) {
  //[CO20171024 OBSOLETE]  atom=str.atoms.at(iat);
  //[CO20171024 OBSOLETE]  atom.cpos=atom.cpos+i*str.lattice(1)+j*str.lattice(2)+k*str.lattice(3);
  //[CO20171024 OBSOLETE]  vatoms.push_back(atom);
  //[CO20171024 OBSOLETE]}
  //[CO20171024 OBSOLETE]    }
  //[CO20171024 OBSOLETE]  }
  //[CO20171024 OBSOLETE]}
  //[CO20171024 OBSOLETE]double nndist=RadiusSphereLattice(str.lattice);
  //[CO20171024 OBSOLETE]for(uint i=0;i<vatoms.size();i++) {
  //[CO20171024 OBSOLETE]  for(uint j=0;j<vatoms.size();j++) {
  //[CO20171024 OBSOLETE]    if(i!=j)
  //[CO20171024 OBSOLETE]if(modulus(vatoms.at(i).cpos-vatoms.at(j).cpos) < nndist) nndist=modulus(vatoms.at(i).cpos-vatoms.at(j).cpos);
  //[CO20171024 OBSOLETE]  }
  //[CO20171024 OBSOLETE]}
  //[CO20171024 OBSOLETE]return nndist;
}

// ***************************************************************************
// Function NearestNeighbors() //DX20201230 - moved from XtalFinder 
// ***************************************************************************
vector<double> NearestNeighbors(const xstructure& xstr){

  // Determine the nearest neighbor distances centered on each atom
  // of the structure (needed for XtalFinder)

  vector<double> all_nn_distances;
  double nn = AUROSTD_MAX_DOUBLE;

  for(uint i=0;i<xstr.atoms.size();i++){
    nn = NearestNeighborToAtom(xstr,i);
    all_nn_distances.push_back(nn);
  }
  return all_nn_distances;
}

// ***************************************************************************
// Function NearestNeighborToAtom() //DX20201230 - moved from XtalFinder
// ***************************************************************************
double NearestNeighborToAtom(const xstructure& xstr, uint k) {

  // Find the minimum interatomic distance in the structure to atom k
  // Different than SYM::minimumDistance(): only considers one atom index
  // in minimization routine, as opposed to global minimimum
  // (considering one atom only affords speed ups)
  // Use resetLatticeDimension() to update search radius for nearest
  // neighbors: once we find a neighbor, update/reduce how far we
  // need to search to find a closer neighbor

  double min_dist=AUROSTD_MAX_DOUBLE;
  double prev_min_dist=0; //DX20190716
  xmatrix<double> lattice = xstr.lattice; //NEW

  //DX speed increase
  //perhaps can speed up even more, since the lattice doesn't change for the xstr...
  vector<xvector<double> > l1, l2, l3;
  vector<int> a_index, b_index, c_index;
  xvector<int> dims(3); //DX20190710 - use robust method
  dims[1]=dims[2]=dims[3]=0; //reset

  xvector<double> tmp_coord, incell_dist, a_component, ab_component; //DX20200329
  double incell_mod=AUROSTD_MAX_DOUBLE;

  uint ii=0, m=0, n=0, p=0, m_size=0, n_size=0, p_size=0;

  for(ii=0; ii<xstr.atoms.size(); ii++){
    if(ii!=k){
      if(min_dist<prev_min_dist){
        if(!(dims[1]==1 && dims[2]==1 && dims[3]==1)){
          // update the dimensions based on new search radius (min_dist)
          resetLatticeDimensions(lattice,min_dist,dims,l1,l2,l3,a_index,b_index,c_index);
          prev_min_dist=min_dist;
          m_size = l1.size(); n_size = l2.size(); p_size = l3.size();
        }
      }
      incell_dist = xstr.atoms[k].cpos-xstr.atoms[ii].cpos;
      incell_mod = aurostd::modulus(incell_dist);
      if(incell_mod<min_dist){
        if(!(dims[1]==1 && dims[2]==1 && dims[3]==1)){
          // update the dimensions based on new search radius (incell_mod)
          resetLatticeDimensions(lattice,incell_mod,dims,l1,l2,l3,a_index,b_index,c_index);
          m_size = l1.size(); n_size = l2.size(); p_size = l3.size();
        }
        prev_min_dist=incell_mod;
      }
      //DX20180423 - running vector in each loop saves computations; fewer duplicate operations
      for(m=0;m<m_size;m++){
        a_component = incell_dist + l1[m];    //DX : coord1-coord2+a*lattice(1)
        for(n=0;n<n_size;n++){
          ab_component = a_component + l2[n]; //DX : coord1-coord2+a*lattice(1) + (b*lattice(2))
          for(p=0;p<p_size;p++){
            tmp_coord = ab_component + l3[p]; //DX : coord1-coord2+a*lattice(1) + (b*lattice(2)) + (c*lattice(3))
            min_dist=aurostd::min(min_dist,aurostd::modulus(tmp_coord));
          }
        }
      }
    }
  }

  return min_dist;
}

// ***************************************************************************
// Function ReScale
// ***************************************************************************
xstructure ReScale(const xstructure& a, const double &in_scale) {
  // This resets scale and changes the cell parameters and coordinates
  // appropriately.  Keeps volume fixed.
  if(in_scale==0.0) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"ReScale()","in_scale must be non zero",_INPUT_ILLEGAL_);}
  xstructure b(a);
  if(aurostd::identical(b.scale,in_scale,_ZERO_TOL_)){return b;}  //try hard not to introduce precision errors, currently we print scale with precision 6
  b.lattice=b.lattice*b.scale/in_scale;
  b.origin=b.origin*b.scale/in_scale;
  b.f2c=trasp(b.lattice);
  b.c2f=inverse(trasp(b.lattice));
  // klattice already contained the scale so it does not need to be fixed
  // b.klattice=b.klattice*in_scale/b.scale;
  for(int i=0;i<(int)b.atoms.size();i++) {
    b.atoms.at(i).fpos=a.atoms.at(i).fpos;
    b.atoms.at(i).cpos=a.atoms.at(i).cpos*b.scale/in_scale;
  }
  if(b.fgroup_calculated) {
    for(int fg=0;fg<(int)b.fgroup.size();fg++) {
      b.fgroup[fg].ctau=b.fgroup[fg].ctau*b.scale/in_scale;
    }
  }
  if(b.sgroup_calculated) {
    for(int sg=0;sg<(int)b.sgroup.size();sg++) {
      b.sgroup[sg].ctau=b.sgroup[sg].ctau*b.scale/in_scale;
      b.sgroup[sg].ctrasl=b.sgroup[sg].ctrasl*b.scale/in_scale;
    }
  }
  b.scale=in_scale;
  b.FixLattices();  // touched scale, then fix the lattices
  return b;
}

void xstructure::ReScale(const double &in_scale) {
  if(in_scale==0.0) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"ReScale()","in_scale must be non zero",_INPUT_ILLEGAL_);}
  if(aurostd::identical(scale,in_scale,_ZERO_TOL_)){return;}  //try hard not to introduce precision errors, currently we print scale with precision 6
  lattice=lattice*scale/in_scale;
  origin=origin*scale/in_scale;
  f2c=trasp(lattice);
  c2f=inverse(trasp(lattice));
  // klattice already contained the scale so it does not need to be fixed
  // klattice=klattice*in_scale/scale;
  for(int i=0;i<(int)atoms.size();i++){
    atoms.at(i).cpos=atoms.at(i).cpos*scale/in_scale;
  }
  if(fgroup_calculated) {
    for(int fg=0;fg<(int)fgroup.size();fg++) {
      fgroup[fg].ctau=fgroup[fg].ctau*scale/in_scale;
    }
  }
  if(sgroup_calculated) {
    for(int sg=0;sg<(int)sgroup.size();sg++) {
      sgroup[sg].ctau=sgroup[sg].ctau*scale/in_scale;
      sgroup[sg].ctrasl=sgroup[sg].ctrasl*scale/in_scale;
    }
  }
  scale=in_scale;
  FixLattices(); // touched scale, then fix the lattices
}

// ***************************************************************************
// Function SetScale
// ***************************************************************************
void xstructure::SetScale(const double &in_scale) {
  if(in_scale==0.0) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"SetScale()","in_scale must be non zero",_INPUT_ILLEGAL_);}
  scale=in_scale;
  FixLattices();  // touched scale, then fix the lattices
}

xstructure SetScale(const xstructure& a,const double &in_scale) {
  // This resets scale.  Keeps volume fixed.
  if(in_scale==0.0) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"SetScale()","in_scale must be non zero",_INPUT_ILLEGAL_);}
  xstructure b;b=a;
  b.scale=in_scale;
  b.FixLattices();  // touched scale, then fix the lattices
  return b;
}

//AS20200514 START
// ***************************************************************************
// Function UpdateCartesianCoordinates
// ***************************************************************************
void xstructure::UpdateCartesianCoordinates()
{
  for (uint at=0; at<atoms.size(); at++) atoms[at].cpos = scale*f2c*atoms[at].fpos;
}
//AS20200514 END

// ***************************************************************************
// Function SetVolume
// ***************************************************************************
void xstructure::SetVolume(const double &in_volume) {
  if(in_volume==0.0){throw aurostd::xerror(__AFLOW_FILE__,"SetVolume()","in_scale must be non zero",_INPUT_ILLEGAL_);} //CO20200201
  if(det(lattice)<0.0){ //CO20200201
    stringstream message; //CO20200201
    message << "Found negative determinant for lattice (det()=" << det(lattice) << "). Flip your basis."; //CO20200201
    throw aurostd::xerror(__AFLOW_FILE__,"SetVolume()",message,_INPUT_ILLEGAL_); //CO20200201
  } //CO20200201
  scale=std::pow((double) in_volume/det(lattice),(double) 1.0/3.0);
  FixLattices();  // touched scale, then fix the lattices
}

xstructure SetVolume(const xstructure& a,const double &in_volume) {
  xstructure b(a);
  b.SetVolume(in_volume);
  return b;
  //[CO20200201 OBSOLETE]if(in_volume==0.0){throw aurostd::xerror(__AFLOW_FILE__,"SetVolume()","in_scale must be non zero",_INPUT_ILLEGAL_);} //CO20200201
  //[CO20200201 OBSOLETE]xstructure b;b=a;
  //[CO20200201 OBSOLETE]if(det(b.lattice)<0.0){ //CO20200201
  //[CO20200201 OBSOLETE]  stringstream message; //CO20200201
  //[CO20200201 OBSOLETE]  message << "Found negative determinant for lattice (det()=" << det(b.lattice) << "). Flip your basis."; //CO20200201
  //[CO20200201 OBSOLETE]  throw aurostd::xerror(__AFLOW_FILE__,"SetVolume()",message,_INPUT_ILLEGAL_); //CO20200201
  //[CO20200201 OBSOLETE]} //CO20200201
  //[CO20200201 OBSOLETE]b.scale=std::pow((double) in_volume/det(b.lattice),(double) 1.0/3.0);
  //[CO20200201 OBSOLETE]b.FixLattices(); // touched scale, need to fix the lattices
  //[CO20200201 OBSOLETE]return b;
}

// ***************************************************************************
// Function SetAutoVolume
// ***************************************************************************
void xstructure::SetAutoVolume(bool use_AFLOW_defaults_in) {  //CO20191010
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream message;

  if(LDEBUG) {
    cerr << __AFLOW_FUNC__ << " xstr=" << endl << (*this) << endl;
    cerr << __AFLOW_FUNC__ << " xstr.species=" << aurostd::joinWDelimiter((*this).species,",") << endl;
    cerr << __AFLOW_FUNC__ << " fixing volume" << endl;
    cerr << __AFLOW_FUNC__ << " volume_orig=" << GetVolume() << endl;
  }
  double volume=0; //,voli=0;
  bool use_AFLOW_defaults=use_AFLOW_defaults_in;
  //try and pull from species_volume first
  for(uint i=0;i<atoms.size()&&!use_AFLOW_defaults;i++){
    for(uint j=0;j<num_each_type.size()&&!use_AFLOW_defaults;j++){
      if(atoms[i].name==species[j]){
        const double& voli=species_volume[j];
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " atoms[i].name=" << atoms[i].name << " atoms[i].vol=" << voli << endl;}
        if(voli==NNN || aurostd::isequal(voli,0.0,_ZERO_TOL_)){use_AFLOW_defaults=true;}
        if(aurostd::isequal(atoms[i].partial_occupation_value,0.0,_ZERO_TOL_)){
          message << "partial_occupation_value==0.0 for atom=" << i;
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
        }
        volume+=atoms[i].partial_occupation_value*voli;
      }
    }
  }
  //otherwise get defaults from AFLOW
  if(use_AFLOW_defaults || abs(volume)<_XPROTO_ZERO_VOL_){
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " using automatic volumes" << endl;}
    volume=0;
    double voli=0.0;
    for(uint i=0;i<atoms.size();i++){
      for(uint j=0;j<num_each_type.size();j++){
        if(atoms[i].name==species[j]){
          voli=GetAtomVolume(atoms[i].name);
          if(voli==NNN || aurostd::isequal(voli,0.0,_ZERO_TOL_)){
            message << "No volume found for " << atoms[i].name << " (auto volumes)";
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
          }
          if(aurostd::isequal(atoms[i].partial_occupation_value,0.0,_ZERO_TOL_)){
            message << "partial_occupation_value==0.0 for atom=" << i;
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
          }
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " atoms[i=" << i << "].name=" << atoms[i].name << " pocc=" << atoms[i].partial_occupation_value << " vol=" << voli << endl;}
          volume+=atoms[i].partial_occupation_value*voli;
        }
      }
    }
  }
  if(abs(volume)<_XPROTO_ZERO_VOL_){
    message << "Final volume==0, check species default volumes";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
  }
  SetVolume(volume);
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " volume_new=" << GetVolume() << endl;}
}

// ***************************************************************************
// Function InflateLattice
// ***************************************************************************
void xstructure::InflateLattice(const double &coefficient) {
  if(coefficient==0.0) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::InflateLattice()","coefficient must be non zero.",_INPUT_ILLEGAL_);}
  //  scale=coefficient*scale;
  lattice=coefficient*lattice;
  FixLattices();  // touched scale/lattice, then fix the lattices
}

xstructure InflateLattice(const xstructure& a,const double &coefficient) {
  // This resets scale.  Keeps volume fixed.
  if(coefficient==0.0) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::InflateLattice()","coefficient must be non zero.",_INPUT_ILLEGAL_);}
  xstructure b;b=a;
  //  b.scale=coefficient*b.scale;
  b.lattice=coefficient*b.lattice;
  b.FixLattices();  // touched scale/lattice, then fix the lattices
  return b;
}

// ***************************************************************************
// Function InflateVolume
// ***************************************************************************
void xstructure::InflateVolume(const double &coefficient) {
  if(coefficient==0.0) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::InflateVolume()","coefficient must be non zero",_INPUT_ILLEGAL_);}
  // scale=std::pow((double) coefficient,(double) 1/3)*scale;
  lattice=std::pow((double) coefficient,(double) 1/3)*lattice;
  FixLattices();  // touched scale/lattice, then fix the lattices
  UpdateCartesianCoordinates(); //AS20200514
}

xstructure InflateVolume(const xstructure& a,const double &coefficient) {
  if(coefficient==0.0) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::InflateVolume()","coefficient must be non zero",_INPUT_ILLEGAL_);}
  xstructure b;b=a;
  //  b.scale=std::pow((double) coefficient,(double) 1/3)*b.scale;
  b.lattice=std::pow((double) coefficient,(double) 1/3)*b.lattice;
  b.FixLattices(); // touched scale/lattice, need to fix the lattices
  b.UpdateCartesianCoordinates(); //AS20200514
  return b;
}


// ***************************************************************************
// Function GetVolume
// ***************************************************************************
double xstructure::GetVolume(void) const {  //CO20200201
  return scale*scale*scale*det(lattice);
}

double GetVolume(const xstructure& a) {
  return a.GetVolume(); //[CO20200520 - OBSOLETE]a.scale*a.scale*a.scale*det(a.lattice);
}

// ***************************************************************************
// Function Volume
// ***************************************************************************
double xstructure::Volume(void) const { //CO20200201
  return GetVolume(); //CO20200201
  //[CO20200201]return scale*scale*scale*det(lattice);
}

double Volume(const xstructure& a) {
  return a.GetVolume(); //CO20200201
  //[CO20200201]return a.scale*a.scale*a.scale*det(a.lattice);
}

_atom BringCloseToOrigin(_atom& atom, xmatrix<double>& f2c){
  _atom atom_out = atom;
  xvector<double> v_in = atom.fpos;
  xvector<int> ijk = atom.ijk;
  for(uint i=1; i<4; i++){
    while(v_in(i)>(1.0-_ZERO_TOL_)){
      v_in(i) = v_in(i)-1;
      ijk(i)-=1;
    }
    while(v_in(i)<=-_ZERO_TOL_){  // fixed DX
      v_in(i) = v_in(i)+1;
      ijk(i)+=1;
    }
  }
  atom_out.fpos = v_in;
  atom_out.cpos = f2c*atom_out.fpos;
  atom_out.ijk = ijk;
  return atom_out;
}

bool uniqueAtomInCell(_atom& atom, deque<_atom>& atoms){
  if(inCell(atom.fpos)){
    if(!alreadyInCell(atom, atoms)){
      return TRUE;
    }
    return FALSE;
  }
  return FALSE;
}

//DX START
//DX20191125 [OBSOLETE - GENERALIZED BELOW] bool inCell(xvector<double>& pos_vec){
//DX20191125 [OBSOLETE - GENERALIZED BELOW]   //DX20180726 [OBSOLETE] if(pos_vec(1)>=-_ZERO_TOL_ && pos_vec(1)<1.0-_ZERO_TOL_ &&
//DX20191125 [OBSOLETE - GENERALIZED BELOW]   //DX20180726 [OBSOLETE]    pos_vec(2)>=-_ZERO_TOL_ && pos_vec(2)<1.0-_ZERO_TOL_ &&
////DX20180726 [OBSOLETE]    pos_vec(3)>=-_ZERO_TOL_ && pos_vec(3)<1.0-_ZERO_TOL_) {      // found something inside
//DX20191125 [OBSOLETE - GENERALIZED BELOW]   if(pos_vec(1)>=-_ZERO_TOL_ && pos_vec(1)<1.0+_ZERO_TOL_ && //DX20180726 - changed from 1.0-_ZERO_TOL_ to 1.0+_ZERO_TOL_
//DX20191125 [OBSOLETE - GENERALIZED BELOW]      pos_vec(2)>=-_ZERO_TOL_ && pos_vec(2)<1.0+_ZERO_TOL_ && //DX20180726 - changed from 1.0-_ZERO_TOL_ to 1.0+_ZERO_TOL_
//DX20191125 [OBSOLETE - GENERALIZED BELOW]      pos_vec(3)>=-_ZERO_TOL_ && pos_vec(3)<1.0+_ZERO_TOL_) {      // found something inside //DX20180726 - changed from 1.0-_ZERO_TOL_ to 1.0+_ZERO_TOL_
//DX20191125 [OBSOLETE - GENERALIZED BELOW]     return TRUE;
//DX20191125 [OBSOLETE - GENERALIZED BELOW]   }
//DX20191125 [OBSOLETE - GENERALIZED BELOW]   return FALSE;
//DX20191125 [OBSOLETE - GENERALIZED BELOW] }
//DX20191125 [OBSOLETE - GENERALIZED BELOW] }
//DX END

// ***************************************************************************
// atomInCell() 
// ***************************************************************************
bool atomInCell(const _atom& atom, double tolerance, double upper_bound, double lower_bound) { //ME+DX20210203 - added bounds

  // check if the atom is in the unit cell based on fractional coordinates
  // if you use the non-default tolerance (i.e., _ZERO_TOL_), this alone is not robust 
  // Note: check over each component and returning false immediately (faster)

  return inCell(atom.fpos, tolerance, upper_bound, lower_bound);

}

// ***************************************************************************
// inCell() 
// ***************************************************************************
// ME20210128 - Added bounds
bool inCell(const xvector<double>& pos_vec, double tolerance, double upper_bound, double lower_bound) {

  // check if the position is in the unit cell based on fractional coordinates
  // if you use the non-default tolerance (i.e., _ZERO_TOL_), this alone is not robust 
  // Note: check over each component and returning false immediately (faster)

  for(uint f=1;f<4;f++){
    // ME20210128: Used to be pos_vec[f] > 1.0 + tolerance.
    // Adjusted to use the same cut-off criterion as bringInCell
    if((pos_vec[f] - upper_bound) >= -tolerance || (pos_vec[f] - lower_bound) < -tolerance){ //allows tunable cutoff
      return false;
    }
  }
  return true;
}


//DX20180726 - check if already in cell - START
bool alreadyInCell(_atom& atom, deque<_atom> atoms){
  for(uint i=0;i<atoms.size();i++){ 
    xvector<double> fdiff = atom.fpos-atoms[i].fpos;
    fdiff = SYM::minimizeDistanceFractionalMethod(fdiff); //DX20190613
    //DX20190613 [OBSOLETE] SYM::PBC(fdiff);
    if(aurostd::abs(fdiff(1))<_ZERO_TOL_ && 
        aurostd::abs(fdiff(2))<_ZERO_TOL_ && 
        aurostd::abs(fdiff(3))<_ZERO_TOL_){
      return TRUE;
    }
  }
  return FALSE;
}
//DX20180726 - check if already in cell - END

// ***************************************************************************
// Function GetSuperCell
// ***************************************************************************
// This funtion creates SuperCells with 9 or 3 elements...
// the old routine by Dane Morgan does not work, so I rewrote this
// one from scratch.
// The algorithm is simple: make the bigger cell
// b.lattice=supercell*a.lattice;
// and generate a bunch of atoms around the old cell and check
// if they are inside the new cell... There is a checksum with ERROR
// if the numbers do not match. Stefano Curtarolo (aug07).
//CO (2017) adding maps to/from supercell/primitive structure, get_symmetry and get_full_basis flags
// get_symmetry will propagate symmetry of primitive cell
// WARNING: if the structure is a derivative structure (i.e., non-uniform expansion, see derivative_structure) then
// NOT all symmetry operations will work, as the symmetry is reduced
// YOU need to check which of these don't work, unless you also calculate get_full_basis,
// which checks the validity of the symmetry operation by finding the full basis map

//CO START
xstructure GetSuperCell(const xstructure& aa, const xmatrix<double> &supercell,vector<int>& sc2pcMap,vector<int>& pc2scMap,
    bool get_symmetry, bool get_full_basis, bool force_supercell_matrix,bool force_strict_pc2scMap) //DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
  //xstructure GetSuperCell(const xstructure& aa, const xmatrix<double> &supercell)
{ //CO20200106 - patching for auto-indenting
  //#define _eps_scell_ 1.0e-5
  //#define _eps_scell_ 0.0
  // check for error
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream message;
  double vol_supercell=det(supercell);
  if(abs(vol_supercell)<0.001){message << "Singular supercell matrix";throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_);}
  xstructure a(aa); a.ReScale(1.0); //the nuclear option, the only way not to mess around with scale EVERYWHERE
  //DO NOT MODIFY STRUCTURE IN HERE, WE WANT TO PROPAGATE SYMMETRY FROM PRIMITIVE STRUCTURE!
  //a.BringInCell();
  xstructure b(a); //actual supercell, need to copy!
  //pflow::CalculateFullSymmetry(a);
  _atom atom,atom2;
  _sym_op pSymOp,fSymOp,aSymOp;  //CO
  int i,j,k; //,dim;
  uint iat;//,at_indx; //CO
  xvector<int> dims(3);
  xvector<double> cshift(3);  //CO
  vector<xvector<double> > cshifts; //CO
  //double zeroTol=1e-10;
  double radius;
  b.lattice=supercell*a.lattice;
  // the scale is kept the same....it is always saved as positive
  b.FixLattices();
  for(i=0;i<(int) b.num_each_type.size();i++) b.num_each_type[i]=0;
  for(i=0;i<(int) b.comp_each_type.size();i++) b.comp_each_type[i]=0;
  b.atoms.clear();
  sc2pcMap.clear(); pc2scMap.clear(); //CO20170722 - clear this everytime

  bool skew = false; //DX20190319 - declared outside loop

  if(get_symmetry){ //DX20190319 - added if statement; don't calc unless necessary
    if(b.dist_nn_min==AUROSTD_NAN){b.dist_nn_min=SYM::minimumDistance(a);} //calculate dist_nn_min
    if(b.sym_eps==AUROSTD_NAN){b.sym_eps=SYM::defaultTolerance(b);}
    skew = SYM::isLatticeSkewed(b.lattice,b.dist_nn_min,b.sym_eps); //DX20190319 - declared above
  } //DX //DX20190319 - added if statement; don't calc unless necessary

  double nx=supercell(1,1),ny=supercell(2,2),nz=supercell(3,3);
  //can only propagate if nx==ny==nz and diagonal
  //FALSE means we have a derivative structure
  //bool derivative_structure=!(abs(supercell(1,1)-supercell(2,2))<zeroTol && abs(supercell(1,1)-supercell(3,3))<zeroTol &&
  //                            abs(supercell(1,2))<zeroTol && abs(supercell(1,3))<zeroTol &&
  //                            abs(supercell(2,1))<zeroTol && abs(supercell(2,3))<zeroTol &&
  //                            abs(supercell(3,1))<zeroTol && abs(supercell(3,2))<zeroTol);
  bool derivative_structure=!(aurostd::isdiagonal(supercell) && 
      aurostd::isequal(nx,ny) && 
      aurostd::isequal(ny,nz));
  //CO
  //VERY IMPORTANT, as we need to reduce the symmetry of the derivative structure
  //getting the basis serves as a validation for the symmetry operator
  //get_full_basis = get_full_basis || derivative_structure;

  //CO START
  // symmetry stuff
  b.ClearSymmetry();  //clear first
  //CO END

  //DX20190319 - added option to expand strictly by uniform supercell matrix - START
  //CO20190409 - note, force_supercell_matrix is PURELY for speed up purposes, the other approach should yield the same results (just slower)
  if(force_supercell_matrix && aurostd::isdiagonal(supercell)){dims[1]=nx;dims[2]=ny;dims[3]=nz;} //[CO20190520 - fixing for dims] && !derivative_structure){dim=nx;} // if here, then nx=ny=nz
  else{
    radius=RadiusSphereLattice(b.lattice);
    dims=LatticeDimensionSphere(a.lattice,radius);//[CO20190520 - EXCESSIVE]dim=max(dims)+1;
  } 

  if(LDEBUG){cerr << __AFLOW_FUNC__ << " dims=" << dims << endl;}
  //DX20190319 - added option to expand strictly by uniform supercell matrix - END
  // if(LDEBUG) cerr << "DEBUG  dims=" << dims << " " << " radius=" << radius << endl;  // DEBUG

  bool match=false;

  //CO20190409 - this pc2scMap issue is more complicated...
  //this is how we resolve: force_strict_pc2scMap
  //if force_strict_pc2scMap==false (default), then the pc2scMap returns the first of the equivalent atoms
  //this is good because of how the algorithm enumerates equivalent atoms: 
  //if you know the supercell size (POCC N_HNF), then the equivalent atoms are next n_hnf atoms
  //however, for APL, we want pc2scMap to return the atom of the original primitive cell (not an equivalent atom)
  //therefore, use force_strict_pc2scMap=true
  //for force_strict_pc2scMap==true, the supercell matrix should be diagonal, otherwise you are not guaranteed to get the i==0 && j==0 && k==0 atom
  //if force_strict_pc2scMap==true and it does not find the i==0 && j==0 && k==0, it throws an error
  //CO20190114 - pc2scMap only makes sense for true supercell expansions
  //there are cases where we use GetSuperCell to convert between representations (see aflow_lattice.cpp)
  //in this case, a mapping is not really possible/useful, so clear out pc2scMap
  bool ignore_pcmap=false;
  bool pcmap=false;

  //CO20181226 - do a check if iatoms really calculated for this cell
  uint atoms_size_check=0;
  for(uint i=0;i<a.iatoms.size();i++){atoms_size_check+=a.iatoms[i].size();}
  if(atoms_size_check!=a.atoms.size()){a.ClearSymmetry();}  //it's possible that the cell was transformed somewhere in aflow, and the symmetry was not cleared

  //CO START
  if(a.iatoms_calculated){
    //CO START
    //create bins
    for(uint ii=0;ii<a.iatoms.size();ii++){
      b.iatoms.push_back(vector<int>(0));
    }
    for(uint ia=0;ia<a.iatoms.size();ia++){
      for(uint iia=0;iia<a.iatoms[ia].size();iia++){
        pcmap=false;
        //[CO20190520 - EXCESSIVE]for(i=-dim;i<=dim;i++) {  //[CO20200106 - close bracket for indenting]}
        //[CO20190520 - EXCESSIVE]  for(j=-dim;j<=dim;j++) {  //[CO20200106 - close bracket for indenting]}
        //[CO20190520 - EXCESSIVE]    for(k=-dim;k<=dim;k++) {  //[CO20200106 - close bracket for indenting]}
        for(i=-dims[1];i<=dims[1];i++) {
          for(j=-dims[2];j<=dims[2];j++) {
            for(k=-dims[3];k<=dims[3];k++) {
              atom=a.atoms[a.iatoms[ia][iia]];
              //cerr << "atom " << a.iatoms[ia][iia] << " fpos_UNrot " << atom.fpos << endl;
              cshift=((double)i)*a.lattice(1)+((double)j)*a.lattice(2)+((double)k)*a.lattice(3);
              atom.cpos=atom.cpos+cshift;
              atom.fpos=b.c2f*atom.cpos; //C2F(b.lattice,atom.cpos);               // put in fractional of new basis
              //  atom.fpos=roundoff(atom.fpos);
              //cerr << "atom " << a.iatoms[ia][iia] << " fpos_rot   " << atom.fpos << endl;
              //cerr << "atom_fpos: " << atom.fpos << endl;
              //DX20180726 - if(inCell(atom.fpos)) //hard cut off
              if(uniqueAtomInCell(atom, b.atoms)) //soft cut off; then check images later //DX20180726
              { //CO20200106 - patching for auto-indenting
                //if(atom.fpos(1)>=-_eps_scell_ && atom.fpos(1)<1.0-_eps_scell_ &&
                //  atom.fpos(2)>=-_eps_scell_ && atom.fpos(2)<1.0-_eps_scell_ &&
                //  atom.fpos(3)>=-_eps_scell_ && atom.fpos(3)<1.0-_eps_scell_)      // found something inside
                // atom=BringInCell(atom,b.lattice);
                b.num_each_type[atom.type]++;                  // CONVASP_MODE
                b.comp_each_type[atom.type]+=atom.partial_occupation_value;                  // CONVASP_MODE

                //we found a new iatom
                if(b.iatoms[ia].empty()){
                  atom.equivalent=b.atoms.size(); //reference self
                  atom.is_inequivalent=TRUE;  //iatom
                } else {
                  //eq atom
                  atom.equivalent=b.iatoms[ia][0];  //reference first iatom
                  atom.is_inequivalent=FALSE; //equivalent atom
                }

                //ijk
                atom.ijk(1)=i;atom.ijk(2)=j;atom.ijk(3)=k;

                //DX20180726 - bring close to origin 
                atom = BringCloseToOrigin(atom, b.f2c);  //updates fpos/cpos/ijk

                b.atoms.push_back(atom);  //do NOT use AddAtom(), AddAtom() rearranges per species and we need to know the mapping
                b.iatoms[ia].push_back(b.atoms.size()-1);
                //save cshifts for fgroups later...
                match=false;
                for(uint cf=0;cf<cshifts.size()&&!match;cf++){
                  if(aurostd::isequal(cshift,cshifts[cf],_ZERO_TOL_)){
                    match=true;
                  }
                }
                if(!match){
                  cshifts.push_back(cshift);
                }

                //mapping
                sc2pcMap.push_back(a.iatoms[ia][iia]);
                // ME20210506 - Strict mapping is done outside to account for non-diagonal supercells
                if(!ignore_pcmap && !pcmap && !force_strict_pc2scMap){
                  //if(force_strict_pc2scMap==true){  //only if i==0 && j==0 && k==0 atom
                  //  if(i==0 && j==0 && k==0){pc2scMap.push_back(b.atoms.size()-1);pcmap=true;}
                  //} else {pc2scMap.push_back(b.atoms.size()-1);pcmap=true;}
                  pc2scMap.push_back(b.atoms.size()-1);pcmap=true;
                }
                //[CO20190116 - OBSOLETE]if(ignore_pcmap==false && i==0 && j==0 && k==0){
                //[CO20190116 - OBSOLETE]  pc2scMap.push_back(b.atoms.size()-1);
                //[CO20190116 - OBSOLETE]  pcmap=true;
                //[CO20190116 - OBSOLETE]}
                //[CO20190116 - OBSOLETE]//pc2scMap is sort of irrelevant, we just need to pick ONE equivalent atom (there are many)
                //[CO20190116 - OBSOLETE]//so just pick the first one
                //[CO20190116 - OBSOLETE]if(!pcmap){
                //[CO20190116 - OBSOLETE]  pc2scMap.push_back(b.atoms.size()-1);
                //[CO20190116 - OBSOLETE]  pcmap=true;
                //[CO20190116 - OBSOLETE]}
                //matching cpos does not work!
                //matching by index (i,j,k) does not work because of inCell()
                //if(!i&&!j&&!k) pc2scMap.push_back(b.atoms.size()-1);
                //if(aurostd::identical(a.atoms[a.iatoms[ia][iia]].cpos,atom.cpos,1e-6)){
                //  pc2scMap[a.iatoms[ia][iia]]=b.atoms.size()-1;
                //}
              }
            }
          }
        }
        // ME20210506 - Strict mapping is done outside to account for non-diagonal supercells
        if(!ignore_pcmap && !pcmap && !force_strict_pc2scMap){
          // ME20210506 - Strict mapping is done outside to account for non-diagonal supercells
          //if(force_strict_pc2scMap){
          //  message << "pc2scMap not found for atom[i=" << a.iatoms[ia][iia] << "]";
          //  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INDEX_MISMATCH_);
          //}
          ignore_pcmap=true;
          pc2scMap.clear();
        }
      }
    }

    // save the number of equivalents
    uint iequivalent=0;
    for(uint iat=0;iat<b.atoms.size();iat++) {
      if(b.atoms[iat].is_inequivalent) {
        b.atoms[iat].num_equivalents=b.iatoms[iequivalent].size();
        b.atoms[iat].index_iatoms=iequivalent;
        iequivalent++;
      }
    }

    b.iatoms_calculated=TRUE;
    //CO END
  } else {
    for(uint ia=0;ia<a.atoms.size();ia++) {
      pcmap=false;
      //[CO20190520 - EXCESSIVE]for(i=-dim;i<=dim;i++){ //[CO20200106 - close bracket for indenting]}
      //[CO20190520 - EXCESSIVE]  for(j=-dim;j<=dim;j++){ //[CO20200106 - close bracket for indenting]}
      //[CO20190520 - EXCESSIVE]    for(k=-dim;k<=dim;k++){ //[CO20200106 - close bracket for indenting]}
      for(i=-dims[1];i<=dims[1];i++){
        for(j=-dims[2];j<=dims[2];j++){
          for(k=-dims[3];k<=dims[3];k++){
            atom=a.atoms[ia];
            //atom.cpos=atom.cpos+(((double)i)*a.lattice(1)+((double)j)*a.lattice(2)+((double)k)*a.lattice(3));
            cshift=((double)i)*a.lattice(1)+((double)j)*a.lattice(2)+((double)k)*a.lattice(3);
            atom.cpos=atom.cpos+cshift;
            atom.fpos=b.c2f*atom.cpos; //C2F(b.lattice,atom.cpos);               // put in fractional of new basis
            //  atom.fpos=roundoff(atom.fpos);
            //DX20180726 - if(inCell(atom.fpos)) //hard cut off
            if(uniqueAtomInCell(atom, b.atoms)) //soft cut off; then check images later //DX20180726
              //if(atom.fpos(1)>=-_eps_scell_ && atom.fpos(1)<1.0-_eps_scell_ &&
              //  atom.fpos(2)>=-_eps_scell_ && atom.fpos(2)<1.0-_eps_scell_ &&
              //  atom.fpos(3)>=-_eps_scell_ && atom.fpos(3)<1.0-_eps_scell_)      // found something inside
              // atom=BringInCell(atom,b.lattice);
            { //CO20200106 - patching for auto-indenting
              b.num_each_type[atom.type]++;                  // CONVASP_MODE
              b.comp_each_type[atom.type]+=atom.partial_occupation_value;                  // CONVASP_MODE
              //ijk
              atom.ijk(1)=i;atom.ijk(2)=j;atom.ijk(3)=k;
              //DX20180726 - bring close to origin 
              atom = BringCloseToOrigin(atom, b.f2c);  //updates fpos/cpos/ijk
              b.atoms.push_back(atom);  //do NOT use AddAtom(), AddAtom() rearranges per species and we need to know the mapping
              //save cshifts for fgroups later...
              match=false;
              for(uint cf=0;cf<cshifts.size()&&!match;cf++){
                if(aurostd::isequal(cshift,cshifts[cf],_ZERO_TOL_)){
                  match=true;
                }
              }
              if(!match){
                cshifts.push_back(cshift);
              }

              //mapping
              sc2pcMap.push_back(ia);
              // ME20210506 - Strict mapping is done outside to account for non-diagonal supercells
              if(!ignore_pcmap && !pcmap && !force_strict_pc2scMap){
                //if(force_strict_pc2scMap==true){  //only if i==0 && j==0 && k==0 atom
                //  if(i==0 && j==0 && k==0){pc2scMap.push_back(b.atoms.size()-1);pcmap=true;}
                //} else {pc2scMap.push_back(b.atoms.size()-1);pcmap=true;}
                pc2scMap.push_back(b.atoms.size()-1);pcmap=true;
              }
              //[CO20190116 - OBSOLETE]if(ignore_pcmap==false && i==0 && j==0 && k==0){
              //[CO20190116 - OBSOLETE]  pc2scMap.push_back(b.atoms.size()-1);
              //[CO20190116 - OBSOLETE]  pcmap=true;
              //[CO20190116 - OBSOLETE]}
              //[CO20190116 - OBSOLETE]//pc2scMap is sort of irrelevant, we just need to pick ONE equivalent atom (there are many)
              //[CO20190116 - OBSOLETE]//so just pick the first one
              //[CO20190116 - OBSOLETE]if(!pcmap){
              //[CO20190116 - OBSOLETE]  pc2scMap.push_back(b.atoms.size()-1);
              //[CO20190116 - OBSOLETE]  pcmap=true;
              //[CO20190116 - OBSOLETE]}
              //matching cpos does not work!
              //matching by index (i,j,k) does not work because of inCell()
              //if(!i&&!j&&!k) pc2scMap.push_back(b.atoms.size()-1);
              //if(aurostd::identical(a.atoms[ia].cpos,atom.cpos,1e-6)){
              //  pc2scMap[a.iatoms[ia][iia]]=b.atoms.size()-1;
              //}
            }
          }
        }
      }
      // ME20210506 - Strict mapping is done outside to account for non-diagonal supercells
      if(!ignore_pcmap && !pcmap && !force_strict_pc2scMap){
        //if(force_strict_pc2scMap){
        //  message << "pc2scMap not found for atom[i=" << ia << "]";
        //  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _INDEX_MISMATCH_);
        //}
        ignore_pcmap=true;
        pc2scMap.clear();
      }
    }
  }
  //CO END
  // ME20210506 - The old method for force_strict_pc2scMap only works for diagonal
  // supercells. This method is brute-force but should work for most non-diagonal
  // cells.
  if (!ignore_pcmap && force_strict_pc2scMap) {
    pc2scMap.clear();
    uint pcatoms = a.atoms.size();
    pc2scMap.resize(pcatoms);
    uint scatoms = b.atoms.size();
    uint nshifts = cshifts.size();
    uint s = 0, ipc = 0, isc = 0;
    for (s = 0; s < nshifts; s++) {
      for (ipc = 0; ipc < pcatoms; ipc++) {
        for (isc = 0; isc < scatoms; isc++) {
          if(aurostd::identical(a.atoms[ipc].cpos + cshifts[s],b.atoms[isc].cpos,_FLOAT_TOL_)) {
            pc2scMap[ipc] = isc;
            break;
          }
        }
        if (isc == scatoms) break;
      }
      if (ipc == pcatoms) break;
    }
    if (s == nshifts) {
      message << "pc2scMap not found";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
  }

  b.GetStoich();  //CO20170724
  b.MakeBasis(); // need to update NUMBER and BASIS

  if(LDEBUG) {
    cerr << __AFLOW_FUNC__ << " sc2pcMap=" << aurostd::joinWDelimiter(sc2pcMap," ") << endl;
    cerr << __AFLOW_FUNC__ << " pc2scMap=" << aurostd::joinWDelimiter(pc2scMap," ") << endl;
  }

  // some check
  double fraction=((double) GetVol(b.lattice)/GetVol(a.lattice));
  double density_a=a.scale*a.scale*a.scale*abs(det(a.lattice))/a.atoms.size();
  double density_b=b.scale*b.scale*b.scale*abs(det(b.lattice))/b.atoms.size();
  if(abs(b.atoms.size()-a.atoms.size()*fraction)>0.1 || abs(density_a-density_b)>0.001) {
    //[CO20190520 - OBSOLETE]cerr << "ERROR   xstructure xstructure::GetSuperCell     " << endl;
    //[CO20190520 - OBSOLETE]cerr << "         supercell has the wrong number of atoms" << endl;
    //[CO20190520 - OBSOLETE]cerr << "         b.atoms.size()     = " << b.atoms.size() << endl;
    //[CO20190520 - OBSOLETE]cerr << "         a.atoms.size()     = " << a.atoms.size() << endl;
    //[CO20190520 - OBSOLETE]cerr << "         b.scale            = " << b.scale << endl;
    //[CO20190520 - OBSOLETE]cerr << "         a.scale            = " << a.scale << endl;
    //[CO20190520 - OBSOLETE]cerr << "         fraction           = " << fraction << endl;
    //[CO20190520 - OBSOLETE]cerr << "         supercell atoms    = " << fraction*a.atoms.size() << endl;
    //[CO20190520 - OBSOLETE]cerr << "         b.density          = " << density_b << endl;
    //[CO20190520 - OBSOLETE]cerr << "         a.density          = " << density_a << endl;
    //[CO20190520 - OBSOLETE]cerr << "         b.lattice          = " << endl;
    //[CO20190520 - OBSOLETE]cerr << b.lattice << endl;
    //[CO20190520 - OBSOLETE]cerr << "         a.lattice          = " << endl;
    //[CO20190520 - OBSOLETE]cerr << a.lattice << endl;

    message << "Supercell has the wrong number of atoms" << endl;
    message << "b.atoms.size()     = " << b.atoms.size() << endl;
    message << "a.atoms.size()     = " << a.atoms.size() << endl;
    message << "b.scale            = " << b.scale << endl;
    message << "a.scale            = " << a.scale << endl;
    message << "fraction           = " << fraction << endl;
    message << "supercell atoms    = " << fraction*a.atoms.size() << endl;
    message << "b.density          = " << density_b << endl;
    message << "a.density          = " << density_a << endl;
    message << "b.lattice          = " << endl;
    message << b.lattice << endl;
    message << "a.lattice          = " << endl;
    message << a.lattice << endl;

    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
  }

  bool pretend_uniform=false;//true;  //CO TEST, REMOVE ME
  get_full_basis = !pretend_uniform && get_full_basis;

  //CO START
  //for now, we focus on pgroup, fgroup, iatoms, and site symmetry (agroup), add as you need
  if(get_symmetry){
    ofstream FileMESSAGE;
    _aflags aflags;
    _kflags kflags;
    bool _write_=false;    //CO no verbose, annoys JJPR
    bool osswrite=false;   //CO no verbose, annoys JJPR
    bool same_pgroups=true;
    bool calculated_pgroups=false;
    bool CALCULATE_FULL_SYMMETRY_ROBUSTLY = false;  //for testing APL
    bool KRUN=true;  //FORCE FULL CALC//false; //CO20181226
    //ostream& oss=cout;  //defined in macro at top of file
    if(derivative_structure){
      KRUN = KRUN && SYM::CalculatePointGroup(FileMESSAGE,b,aflags,_write_,osswrite,oss);
      if(LDEBUG) {
        if(!KRUN){
          cerr << "Symmetry propagation FAILED with derivative structure at point group" << endl;
        } else {
          cerr << "Symmetry propagation PASSED with derivative structure at point group" << endl;
        }
      }
      same_pgroups=(KRUN && SYM::PointGroupsIdentical(a.pgroup,b.pgroup,b.sym_eps,false));  //DX20171207 - added is_same_lattice
      calculated_pgroups=KRUN;
    }
    if((KRUN && !pretend_uniform && !same_pgroups) || CALCULATE_FULL_SYMMETRY_ROBUSTLY){
      if(CALCULATE_FULL_SYMMETRY_ROBUSTLY){cerr << "Calculating symmetry of supercell robustly" << endl;}
      //we calculate earlier to see if there's a mismatch
      //SYM::CalculatePointGroup(FileMESSAGE,b,aflags,_write_,osswrite,oss);
      //do all at same sym_eps as primitive cell
      KRUN = KRUN && SYM::CalculateFactorGroup(FileMESSAGE,b,aflags,_write_,osswrite,oss);
      if(LDEBUG) {
        if(!KRUN){
          cerr << "Symmetry propagation FAILED with derivative structure at factor group" << endl;
        } else {
          cerr << "Symmetry propagation PASSED with derivative structure at factor group" << endl;
        }
      }
      KRUN = KRUN && SYM::CalculatePointGroupCrystal(FileMESSAGE,b,aflags,_write_,osswrite,oss);
      if(LDEBUG) {
        if(!KRUN){
          cerr << "Symmetry propagation FAILED with derivative structure at point group crystal" << endl;
        } else {
          cerr << "Symmetry propagation PASSED with derivative structure at point group crystal" << endl;
        }
      }
      //if(!a.iatoms_calculated){
      KRUN = KRUN && SYM::CalculateInequivalentAtoms(FileMESSAGE,b,aflags,_write_,osswrite,oss); //100% necessary, new pgroups means different symmetry, different iatoms
      if(LDEBUG) {
        if(!KRUN){
          cerr << "Symmetry propagation FAILED with derivative structure at iatoms" << endl;
        } else {
          cerr << "Symmetry propagation PASSED with derivative structure at iatoms" << endl;
        }
      }
      int agroup_calculation_mode=(get_full_basis ? 0 : 1 );
      if(CALCULATE_FULL_SYMMETRY_ROBUSTLY){agroup_calculation_mode=2;}
      //}
      //AGAIN, many fgroups, but not many pgroups for derivative structures, let's see if this is faster...
      //if(!SYM::CalculateSitePointGroup(FileMESSAGE,b,true,aflags,_write_,osswrite,oss)){  //iatoms only  
      KRUN = KRUN && SYM::CalculateSitePointGroup(FileMESSAGE,b,agroup_calculation_mode,aflags,_write_,osswrite,oss);  //don't waste time calculating basis_map for eatoms, really never use them anyway
      if(LDEBUG) {
        if(!KRUN){
          cerr << "Symmetry propagation FAILED with derivative structure at agroup" << endl;
        } else {
          cerr << "Symmetry propagation PASSED with derivative structure at agroup" << endl;
        }
      }
      //}
      //validate that we have good symmetry here
    } else if(KRUN){
      //////////////////////////////////////////////////////////////////////////
      //PGROUP
      if(KRUN && a.pgroup_calculated && !calculated_pgroups){
        for(uint i=0;i<a.pgroup.size() && KRUN;i++){
          pSymOp=a.pgroup[i];
          pSymOp.Uf=b.c2f * pSymOp.Uc * b.f2c;
          pSymOp.basis_atoms_map.clear();
          pSymOp.basis_types_map.clear();
          pSymOp.basis_map_calculated=FALSE;
          //CO+DX, getFullSymBasis does not make sense for just rotations
          //if(get_full_basis){
          //  KRUN = KRUN && SYM::getFullSymBasis(b.atoms,b.lattice,b.c2f,b.f2c,pSymOp,FALSE,skew,b.sym_eps,pSymOp.basis_atoms_map,pSymOp.basis_types_map);
          //  if(LDEBUG) {
          //    if(!KRUN){
          //      cerr << "Symmetry propagation FAILED with uniform supercell structure at point group" << endl;
          //    } else {
          //      cerr << "Symmetry propagation PASSED with uniform supercell structure at point group" << endl;
          //}
          //  }
          //  pSymOp.basis_map_calculated=KRUN;
          //}
          //if(KRUN){b.pgroup.push_back(pSymOp);}
          if(KRUN){SYM::AddSymmetryToStructure(b,pSymOp.Uc,pSymOp.Uf,pSymOp.ctau,pSymOp.ftau,pSymOp.ctrasl,pSymOp.ftrasl,pSymOp.basis_atoms_map,pSymOp.basis_types_map,pSymOp.basis_map_calculated,_PGROUP_,FALSE);}  //CO20170706 - make sure quaternion is updated
        }
        //if(KRUN && b.pgroup.size()){
        //b.pgroup_calculated=true;
        b.pgroup_calculated=(KRUN && b.pgroup.size());
        //}
      }
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      //FGROUP
      if(KRUN && a.fgroup_calculated) { 
        for(uint fg=0;fg<a.fgroup.size() && KRUN;fg++){
          for(uint cs=0;cs<cshifts.size() && KRUN;cs++){
            fSymOp=a.fgroup[fg];
            fSymOp.ctau=fSymOp.ctau+cshifts[cs];
            fSymOp.ftau=b.c2f*fSymOp.ctau;
            if(inCell(fSymOp.ftau)){ //DX CHANGE HERE; NO MORE TOL_ABC_RES
              // We have to correct the Uf for each symop since we have changed the lattice...
              fSymOp.Uf=b.c2f * fSymOp.Uc * b.f2c;
              fSymOp.basis_atoms_map.clear();
              fSymOp.basis_types_map.clear();
              fSymOp.basis_map_calculated=FALSE;
              for(uint iii=0;iii<b.atoms.size();iii++){
                fSymOp.basis_atoms_map.push_back(0);
                fSymOp.basis_types_map.push_back(0);
              }
              fSymOp.basis_map_calculated=false;
              //calculate basis_atoms_map and basis_types_map
              if(get_full_basis) {
                //NOPE, we will calculate
                KRUN = KRUN && SYM::getFullSymBasis(b.atoms,b.lattice,b.c2f,b.f2c,fSymOp,TRUE,skew,b.sym_eps,fSymOp.basis_atoms_map,fSymOp.basis_types_map);
                if(LDEBUG) {
                  if(!KRUN){
                    cerr << "Symmetry propagation FAILED with uniform supercell structure at factor group" << endl;
                  } else {
                    cerr << "Symmetry propagation PASSED with uniform supercell structure at factor group" << endl;
                  }
                }
                //if(!SYM::getFullSymBasis(b.atoms,b.lattice,b.c2f,b.f2c,fSymOp,TRUE,skew,b.sym_eps,fSymOp.basis_atoms_map,fSymOp.basis_types_map)) {
                //cerr << "Unable to find atom/types basis for fgroup" << endl;
                //cerr << fSymOp << endl;
                //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_);
                //KRUN = FALSE;
                //}
                fSymOp.basis_map_calculated=KRUN;
              }
              //if(KRUN){b.fgroup.push_back(fSymOp);}
              if(KRUN){SYM::AddSymmetryToStructure(b,fSymOp.Uc,fSymOp.Uf,fSymOp.ctau,fSymOp.ftau,fSymOp.ctrasl,fSymOp.ftrasl,fSymOp.basis_atoms_map,fSymOp.basis_types_map,fSymOp.basis_map_calculated,_FGROUP_,FALSE);}  //CO20170706 - make sure quaternion is updated
            }
          }
        }
        //if(b.fgroup.size()){
        //  b.fgroup_calculated=TRUE;   //there are more, but these are the important ones
        b.fgroup_calculated=(KRUN && b.fgroup.size());   //there are more, but these are the important ones
        //}
      }
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      //PGROUP_XTAL
      if(KRUN && a.pgroup_xtal_calculated){
        for(uint i=0;i<a.pgroup_xtal.size() && KRUN;i++){
          pSymOp=a.pgroup_xtal[i];
          pSymOp.Uf=b.c2f * pSymOp.Uc * b.f2c;
          pSymOp.basis_atoms_map.clear();
          pSymOp.basis_types_map.clear();
          pSymOp.basis_map_calculated=FALSE;
          //CO+DX, getFullSymBasis does not make sense for just rotations
          //if(get_full_basis){
          //  KRUN = KRUN && SYM::getFullSymBasis(b.atoms,b.lattice,b.c2f,b.f2c,pSymOp,FALSE,skew,b.sym_eps,pSymOp.basis_atoms_map,pSymOp.basis_types_map);
          //  if(LDEBUG) {
          //    if(!KRUN){
          //      cerr << "Symmetry propagation FAILED with uniform supercell structure at point group crystal" << endl;
          //    } else {
          //      cerr << "Symmetry propagation PASSED with uniform supercell structure at point group crystal" << endl;
          //    }
          //  }
          //  pSymOp.basis_map_calculated=KRUN;
          //}
          //if(KRUN){b.pgroup_xtal.push_back(pSymOp);}
          if(KRUN){SYM::AddSymmetryToStructure(b,pSymOp.Uc,pSymOp.Uf,pSymOp.ctau,pSymOp.ftau,pSymOp.ctrasl,pSymOp.ftrasl,pSymOp.basis_atoms_map,pSymOp.basis_types_map,pSymOp.basis_map_calculated,_PGROUP_XTAL_,FALSE);}  //CO20170706 - make sure quaternion is updated
        }
        //if(KRUN && b.pgroup_xtal.size()){
        //b.pgroup_xtal_calculated=true;
        b.pgroup_xtal_calculated=(KRUN && b.pgroup_xtal.size());
        //}
      }
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      //IF A.IATOMS, then we already propagated, otherwise, just run full routine for iatoms
      if(KRUN && !a.iatoms_calculated){
        KRUN = KRUN && SYM::CalculateInequivalentAtoms(FileMESSAGE,b,aflags,_write_,osswrite,oss);  //do all atoms
        if(LDEBUG) {
          if(!KRUN){
            cerr << "Symmetry propagation FAILED with uniform supercell structure at agroup" << endl;
          } else {
            cerr << "Symmetry propagation PASSED with uniform supercell structure at agroup" << endl;
          }
        }
      }
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      //AGROUP
      if(KRUN && a.agroup_calculated){
        if(b.iatoms_calculated && b.fgroup_calculated && get_full_basis){ //the only way we get a speed up is if we can use fgroup.basis_atoms_map
          //xstructure bb=b;
          deque<_atom> b_atoms=b.atoms;
          xvector<double> origin(3),frigin(3);
          //create space for agroups
          for(uint iia=0;iia<b.atoms.size();iia++){
            b.agroup.push_back(vector<_sym_op>(0));
          }
          for(uint ia=0;ia<b.iatoms.size() && KRUN;ia++){
            iat=b.iatoms[ia][0];
            //let's recycle what we have to improve speed
            origin=b.atoms[iat].cpos;
            frigin=b.atoms[iat].fpos;
            for(uint ii=0;ii<b.atoms.size();ii++){
              //go back to original first, then subtract new origin
              b_atoms[ii].cpos=b.atoms[ii].cpos-origin;
              b_atoms[ii].fpos=b.atoms[ii].fpos-frigin;
              //now bring in cell
              b_atoms[ii]=BringInCell(b_atoms[ii],b.lattice);
            }
            //bb.ShiftOriginToAtom(iat);
            //bb.BringInCell();

            //IATOMS ONLY
            for(uint iia=0;iia<a.agroup[sc2pcMap[iat]].size() && KRUN;iia++){
              aSymOp=a.agroup[sc2pcMap[iat]][iia];
              // We have to correct the Uf for each symop since we have changed the lattice...
              aSymOp.Uf=b.c2f * aSymOp.Uc * b.f2c;
              //no longer necessary since we force a basis map calculation
              //aSymOp.basis_atoms_map.clear();
              //aSymOp.basis_types_map.clear();
              //aSymOp.basis_map_calculated=FALSE;
              //for(uint iii=0;iii<b.atoms.size();iii++){
              //aSymOp.basis_atoms_map.push_back(0);
              //aSymOp.basis_types_map.push_back(0);
              //}
              //fSymOp.basis_map_calculated=false;
              KRUN = KRUN && SYM::getFullSymBasis(b_atoms,b.lattice,b.c2f,b.f2c,aSymOp,TRUE,skew,b.sym_eps,aSymOp.basis_atoms_map,aSymOp.basis_types_map);
              if(LDEBUG) {
                if(!KRUN){
                  cerr << "Symmetry propagation FAILED with uniform supercell structure at agroup" << endl;
                } else {
                  cerr << "Symmetry propagation PASSED with uniform supercell structure at agroup" << endl;
                }
              }
              //if(!SYM::getFullSymBasis(b_atoms,b.lattice,b.c2f,b.f2c,aSymOp,TRUE,skew,b.sym_eps,aSymOp.basis_atoms_map,aSymOp.basis_types_map)) {
              //cerr << "Unable to find atom/types basis for agroup" << endl;
              //cerr << aSymOp << endl;
              //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_);
              //KRUN = FALSE;
              //}
              aSymOp.basis_map_calculated=KRUN;
              //if(KRUN){b.agroup[iat].push_back(aSymOp);}
              if(KRUN){SYM::AddSymmetryToStructure(b,iat,aSymOp.Uc,aSymOp.Uf,aSymOp.ctau,aSymOp.ftau,aSymOp.ctrasl,aSymOp.ftrasl,aSymOp.basis_atoms_map,aSymOp.basis_types_map,aSymOp.basis_map_calculated,_AGROUP_,FALSE);}  //CO20170706 - make sure quaternion is updated
            }
          }

          //EATOMS FOLLOW
          KRUN = KRUN && SYM::CalculateSitePointGroup_EquivalentSites(b,get_full_basis,b.sym_eps);
          if(LDEBUG) {
            if(!KRUN){
              cerr << "Symmetry propagation FAILED with uniform supercell structure at agroup equivalent" << endl;
            } else {
              cerr << "Symmetry propagation PASSED with uniform supercell structure at agroup equivalent" << endl;
            }
          }
          //KRUN = FALSE;
          //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unable to propagate site symmetry to equivalent atoms.",_RUNTIME_ERROR_);
        } else {
          //can be faster than procedure above because there are MANY fgroups
          KRUN = KRUN && SYM::CalculateSitePointGroup(FileMESSAGE,b,1,aflags,_write_,osswrite,oss);  //we already know get_full_basis==FALSE, so don't waste time calculating for eatoms
          if(LDEBUG) {
            if(!KRUN){
              cerr << "Symmetry propagation FAILED with uniform supercell structure at agroup" << endl;
            } else {
              cerr << "Symmetry propagation PASSED with uniform supercell structure at agroup" << endl;
            }
          }
        }
        //if(b.agroup.size() && b.agroup[0].size()){      //just a fast check to see we have agroups somewhere (we should always get identity)
        b.agroup_calculated=(KRUN && b.agroup[0].size()); //just a fast check to see we have agroups somewhere (we should always get identity)
        //  b.agroup_calculated=TRUE;
        //}
      }
      //////////////////////////////////////////////////////////////////////////
    }
    if(!KRUN){
      oss << (aflags.QUIET?"":"00000  MESSAGE ") << "SUPERCELL Symmetry propagation FAILED" << Message(__AFLOW_FILE__,aflags) << endl;
      oss << (aflags.QUIET?"":"00000  MESSAGE ") << "SUPERCELL Symmetry retrying with symmetry scan" << Message(__AFLOW_FILE__,aflags) << endl;
      b.ClearSymmetry(); //CO20181226
      pflow::PerformFullSymmetry(b,FileMESSAGE,aflags,kflags,osswrite,oss);
      //FOOLPROOF!!!!!!!!!
      //no need for krun here with force_perform: if it fails with full scan, it will calculate at default tolerance and keep going!
    }
  }

  //cerr << "prim: " << a.pgroup.size() << ", sup: " << b.pgroup.size() << endl;    //CO REMOVE

  //CO END

  b.ReScale(aa.scale);  //the nuclear option, the only way not to mess around with scale EVERYWHERE
  return b;
}
//CO END

xstructure GetSuperCell(const xstructure& a, const xvector<double>& supercell,vector<int>& sc2pcMap,vector<int>& pc2scMap,
    bool get_symmetry, bool get_full_basis, bool force_supercell_matrix,bool force_strict_pc2scMap) //DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
  //xstructure GetSuperCell(const xstructure& a, const xvector<double>& supercell)
{ //CO20200106 - patching for auto-indenting
  xmatrix<double> _supercell(3,3);
  if(supercell.rows==9) {
    _supercell(1,1)=supercell(1);_supercell(1,2)=supercell(2);_supercell(1,3)=supercell(3);
    _supercell(2,1)=supercell(4);_supercell(2,2)=supercell(5);_supercell(2,3)=supercell(6);
    _supercell(3,1)=supercell(7);_supercell(3,2)=supercell(8);_supercell(3,3)=supercell(9);
    return GetSuperCell(a,_supercell,sc2pcMap,pc2scMap,get_symmetry,get_full_basis,force_supercell_matrix,force_strict_pc2scMap); //DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
  }
  if(supercell.rows==3) {
    _supercell(1,1)=supercell(1);_supercell(2,2)=supercell(2);_supercell(3,3)=supercell(3);
    return GetSuperCell(a,_supercell,sc2pcMap,pc2scMap,get_symmetry,get_full_basis,force_supercell_matrix,force_strict_pc2scMap); //DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
  }
  stringstream message;
  message << "Matrix must have 9 or 3 elements";
  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_);
}

xstructure GetSuperCell(const xstructure& a, const xvector<int>& supercell,vector<int>& sc2pcMap,vector<int>& pc2scMap,
    bool get_symmetry, bool get_full_basis, bool force_supercell_matrix,bool force_strict_pc2scMap) //DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
  //xstructure GetSuperCell(const xstructure& a, const xvector<int>& supercell)
{ //CO20200106 - patching for auto-indenting
  xmatrix<double> _supercell(3,3);
  if(supercell.rows==9) {
    _supercell(1,1)=supercell(1);_supercell(1,2)=supercell(2);_supercell(1,3)=supercell(3);
    _supercell(2,1)=supercell(4);_supercell(2,2)=supercell(5);_supercell(2,3)=supercell(6);
    _supercell(3,1)=supercell(7);_supercell(3,2)=supercell(8);_supercell(3,3)=supercell(9);
    return GetSuperCell(a,_supercell,sc2pcMap,pc2scMap,get_symmetry,get_full_basis,force_supercell_matrix,force_strict_pc2scMap); //DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
  }
  if(supercell.rows==3) {
    _supercell(1,1)=supercell(1);_supercell(2,2)=supercell(2);_supercell(3,3)=supercell(3);
    return GetSuperCell(a,_supercell,sc2pcMap,pc2scMap,get_symmetry,get_full_basis,force_supercell_matrix,force_strict_pc2scMap); //DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
  }
  stringstream message;
  message << "Matrix must have 9 or 3 elements";
  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_);
}

xstructure GetSuperCell(const xstructure& a, const int& sc11,const int& sc12,const int& sc13, const int& sc21,const int& sc22,const int& sc23, const int& sc31,const int& sc32,const int& sc33,vector<int>& sc2pcMap,vector<int>& pc2scMap,
    bool get_symmetry, bool get_full_basis, bool force_supercell_matrix,bool force_strict_pc2scMap) //DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
  //xstructure GetSuperCell(const xstructure& a, const int& sc11,const int& sc12,const int& sc13, const int& sc21,const int& sc22,const int& sc23, const int& sc31,const int& sc32,const int& sc33)
{ //CO20200106 - patching for auto-indenting
  xmatrix<double> _supercell(3,3);
  _supercell.clear();
  _supercell(1,1)=(double) sc11;_supercell(1,2)=(double) sc12;_supercell(1,3)=(double) sc13;
  _supercell(2,1)=(double) sc21;_supercell(2,2)=(double) sc22;_supercell(2,3)=(double) sc23;
  _supercell(3,1)=(double) sc31;_supercell(3,2)=(double) sc32;_supercell(3,3)=(double) sc33;
  return GetSuperCell(a,_supercell,sc2pcMap,pc2scMap,get_symmetry,get_full_basis,force_supercell_matrix,force_strict_pc2scMap); //DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
}

xstructure GetSuperCell(const xstructure& a, const int& sc1,const int& sc2,const int& sc3,vector<int>& sc2pcMap,vector<int>& pc2scMap,
    bool get_symmetry, bool get_full_basis, bool force_supercell_matrix,bool force_strict_pc2scMap) //DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
  //xstructure GetSuperCell(const xstructure& a, const int& sc1,const int& sc2,const int& sc3)
{ //CO20200106 - patching for auto-indenting
  xmatrix<double> _supercell(3,3);
  _supercell.clear();
  _supercell(1,1)=(double) sc1;
  _supercell(2,2)=(double) sc2;
  _supercell(3,3)=(double) sc3;
  return GetSuperCell(a,_supercell,sc2pcMap,pc2scMap,get_symmetry,get_full_basis,force_supercell_matrix,force_strict_pc2scMap); //DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
}

//CO START
xstructure GetSuperCell(const xstructure& aa, const xmatrix<double> &supercell) {
  vector<int> sc2pcMap;
  vector<int> pc2scMap;
  bool get_symmetry=false;
  bool get_full_basis=false;
  bool force_supercell_matrix=false; //DX20190319 - added force_supercell_matrix
  bool force_strict_pc2scMap=false; //CO20190409 - added force_strict_pc2scMap
  return GetSuperCell(aa,supercell,sc2pcMap,pc2scMap,get_symmetry,get_full_basis,force_supercell_matrix,force_strict_pc2scMap); //DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
}

xstructure GetSuperCell(const xstructure& a, const xvector<double>& supercell) {
  vector<int> sc2pcMap;
  vector<int> pc2scMap;
  bool get_symmetry=false;
  bool get_full_basis=false;
  bool force_supercell_matrix=false; //DX20190319 - added force_supercell_matrix
  bool force_strict_pc2scMap=false; //CO20190409 - added force_strict_pc2scMap
  return GetSuperCell(a,supercell,sc2pcMap,pc2scMap,get_symmetry,get_full_basis,force_supercell_matrix,force_strict_pc2scMap); //DX20190319 - added force_supercell_matrix //CO20190409 - added force_strict_pc2scMap
}

xstructure GetSuperCell(const xstructure& a, const xvector<int>& supercell) {
  vector<int> sc2pcMap;
  vector<int> pc2scMap;
  bool get_symmetry=false;
  bool get_full_basis=false;
  bool force_supercell_matrix=false; //DX20190319 - added force_supercell_matrix
  bool force_strict_pc2scMap=false; //CO20190409 - added force_strict_pc2scMap
  return GetSuperCell(a,supercell,sc2pcMap,pc2scMap,get_symmetry,get_full_basis,force_supercell_matrix,force_strict_pc2scMap);  //DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
}

xstructure GetSuperCell(const xstructure& a, const int& sc11,const int& sc12,const int& sc13, const int& sc21,const int& sc22,const int& sc23, const int& sc31,const int& sc32,const int& sc33) {
  vector<int> sc2pcMap;
  vector<int> pc2scMap;
  bool get_symmetry=false;
  bool get_full_basis=false;
  bool force_supercell_matrix=false; //DX20190319 - added force_supercell_matrix
  bool force_strict_pc2scMap=false; //CO20190409 - added force_strict_pc2scMap
  return GetSuperCell(a,sc11,sc12,sc13,sc21,sc22,sc23,sc31,sc32,sc33,sc2pcMap,pc2scMap,get_symmetry,get_full_basis,force_supercell_matrix,force_strict_pc2scMap); //DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
}

xstructure GetSuperCell(const xstructure& a, const int& sc1,const int& sc2,const int& sc3) {
  vector<int> sc2pcMap;
  vector<int> pc2scMap;
  bool get_symmetry=false;
  bool get_full_basis=false;
  bool force_supercell_matrix=false; //DX20190319 - added force_supercell_matrix
  bool force_strict_pc2scMap=false; //CO20190409 - added force_strict_pc2scMap
  return GetSuperCell(a,sc1,sc2,sc3,sc2pcMap,pc2scMap,get_symmetry,get_full_basis,force_supercell_matrix,force_strict_pc2scMap); //DX20190319 - added force_supercell_matrix //CO20190409 - added force_strict_pc2scMap
}
//CO END


// ***************************************************************************
// Function ClearSymmetry
// ***************************************************************************
void xstructure::ClearSymmetry(void) {

  bool LDEBUG = (FALSE || XHOST.DEBUG); //DX20210406

  // PGROUP ----------------------------
  pgroup.clear();            // just initialize
  pgroup_calculated=FALSE;
  bravais_lattice_lattice_type="";bravais_lattice_lattice_variation_type="";bravais_lattice_lattice_system=""; //DX20210430 - missing
  // PGROUP_XTAL ----------------------------
  pgroup_xtal.clear();        // just initialize
  pgroup_xtal_calculated=FALSE;
  crystal_family="";crystal_system="";point_group_crystal_class="";
  point_group_Shoenflies="";point_group_Hermann_Mauguin="";point_group_orbifold="";
  point_group_type="";point_group_order="";point_group_structure="";
  pearson_symbol=""; //DX20210430 - missing
  bravais_lattice_type="";bravais_lattice_variation_type="";bravais_lattice_system=""; //DX20210430 - missing
  // PGROUPK_PATTERSON ---------------------------- //DX20200129
  pgroupk_Patterson.clear();        // just initialize
  pgroupk_Patterson_calculated=FALSE;
  // PGROUPK ----------------------------
  pgroupk.clear();            // just initialize
  pgroupk_calculated=FALSE;
  reciprocal_lattice_type="";reciprocal_lattice_variation_type=""; //DX20210430 - missing
  // PGROUPK_XTAL ----------------------------
  pgroupk_xtal.clear();            // just initialize //DX20171205 - Added pgroupk_xtal
  pgroupk_xtal_calculated=FALSE;                      //DX20171205 - Added pgroupk_xtal
  // FGROUP ----------------------------
  fgroup.clear();            // just initialize
  fgroup_calculated=FALSE;
  // SGROUP ----------------------------
  sgroup_radius=-_calculate_symmetry_default_sgroup_radius_; // symmetry not calculated
  sgroup_radius_dims.clear();
  sgroup.clear();            // just initialize
  sgroup_calculated=FALSE;
  // SITE POINT GROUP ------------------
  agroup_calculated=FALSE;
  for(uint i=0;i<agroup.size();i++)
    agroup.at(i).clear();
  agroup.clear();
  // INEQUIVALENT ATOMS ----------------
  iatoms_calculated=FALSE;
  for(uint i=0;i<iatoms.size();i++)
    iatoms.at(i).clear();
  iatoms.clear();
  for(uint i=0;i<atoms.size();i++){atoms[i].ClearSymmetry();} //CO20190219
  // SUPERLATTICE //DX20210430 - missing
  bravais_superlattice_lattice.clear();
  bravais_superlattice_type="";bravais_superlattice_variation_type="";
  bravais_superlattice_system="";pearson_symbol_superlattice="";

  if(LDEBUG){ cerr << XPID << "xstructure::ClearSymmetry(): All symmetry attributes have been cleared." << endl; } //DX20210406

}

//DX - Consider using pflow::CalculateFullSymmetry in aflow_aconvasp_main.cpp.
//      It contains consistency checks for the symmetry analysis.
// ***************************************************************************
// Function CalculateSymmetry
// ***************************************************************************
bool xstructure::CalculateSymmetry(bool ossverbose,double radius) {
  ofstream FileDevNull("/dev/null");

  _aflags aflags; aflags.Directory="./"; aflags.QUIET=TRUE;
  _kflags kflags; pflow::defaultKFlags4SymWrite(kflags,false); pflow::defaultKFlags4SymCalc(kflags,true);

  (*this).LatticeReduction_avoid=TRUE;
  (*this).sgroup_radius=radius;

  return pflow::PerformFullSymmetry(*this,FileDevNull,aflags,kflags,ossverbose,oss);

  //SYM::CalculatePointGroup(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
  //SYM::CalculatePointGroupKLattice(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
  //SYM::CalculateSitePointGroup(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
  //SYM::CalculateFactorGroup(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
  //SYM::CalculateInequivalentAtoms(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
  //SYM::CalculateSpaceGroup(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
}

bool xstructure::CalculateSymmetry(void) {
  LatticeReduction_avoid=TRUE;
  return CalculateSymmetry(FALSE,_calculate_symmetry_default_sgroup_radius_*MaxStructureLattice(*this));
}

bool CalculateSymmetry(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose,double radius) {
  ofstream FileDevNull("/dev/null");

  _aflags aflags; aflags.Directory="./"; aflags.QUIET=TRUE;
  _kflags kflags; pflow::defaultKFlags4SymWrite(kflags,fffverbose); pflow::defaultKFlags4SymCalc(kflags,true);

  str.LatticeReduction_avoid=TRUE;
  str.sgroup_radius=radius;

  return pflow::PerformFullSymmetry(str,FileDevNull,aflags,kflags,ossverbose,oss);
  //SYM::CalculatePointGroup(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
  //SYM::CalculatePointGroupKLattice(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
  //SYM::CalculateSitePointGroup(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
  //SYM::CalculateFactorGroup(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
  //SYM::CalculateInequivalentAtoms(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
  //SYM::CalculateSpaceGroup(FileDevNull,str,aflags,TRUE,ossverbose,oss);
}

bool CalculateSymmetry(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose) {
  str.LatticeReduction_avoid=TRUE;
  return CalculateSymmetry(str,ossverbose,oss,fffverbose,_calculate_symmetry_default_sgroup_radius_*MaxStructureLattice(str));
}

bool CalculateSymmetry(xstructure& str,bool ossverbose,ostream& oss,double radius) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  return CalculateSymmetry(str,ossverbose,oss,FALSE,radius);
}

bool CalculateSymmetry(xstructure& str,bool ossverbose,double radius) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  return CalculateSymmetry(str,ossverbose,cout,FALSE,radius);
}

bool CalculateSymmetry(xstructure& str,double radius) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  return CalculateSymmetry(str,FALSE,cout,FALSE,radius);
}

bool CalculateSymmetry(xstructure& str,bool ossverbose) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  return CalculateSymmetry(str,ossverbose,cout,FALSE,_calculate_symmetry_default_sgroup_radius_*MaxStructureLattice(str));
}

bool CalculateSymmetry(xstructure& str) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  return CalculateSymmetry(str,FALSE,cout,FALSE,_calculate_symmetry_default_sgroup_radius_*MaxStructureLattice(str));
}

//DX - Consider using pflow::CalculateFullSymmetry in aflow_aconvasp_main.cpp.
//      It contains consistency checks for the symmetry analysis.
// ***************************************************************************
// Function CalculateSymmetryPointGroup
// ***************************************************************************
void xstructure::CalculateSymmetryPointGroup(bool ossverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  LatticeReduction_avoid=TRUE;
  SYM::CalculatePointGroup(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
}

void xstructure::CalculateSymmetryPointGroup(void) {
  LatticeReduction_avoid=TRUE;
  CalculateSymmetryPointGroup(FALSE);
}

void CalculateSymmetryPointGroup(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  str.LatticeReduction_avoid=TRUE;
  SYM::CalculatePointGroup(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
}

void CalculateSymmetryPointGroup(xstructure& str,bool ossverbose,ostream& oss) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroup(str,ossverbose,oss,FALSE);
}

void CalculateSymmetryPointGroup(xstructure& str,bool ossverbose) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroup(str,ossverbose,cout,FALSE);
}

void CalculateSymmetryPointGroup(xstructure& str) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroup(str,FALSE,cout,FALSE);
}

//DX - Consider using pflow::CalculateFullSymmetry in aflow_aconvasp_main.cpp.
//      It contains consistency checks for the symmetry analysis.
// ***************************************************************************
// Function CalculateSymmetryPointGroupCrystal
// ***************************************************************************
void xstructure::CalculateSymmetryPointGroupCrystal(bool ossverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  LatticeReduction_avoid=TRUE;
  if(pgroup_calculated==FALSE) SYM::CalculatePointGroup(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
  if(fgroup_calculated==FALSE) SYM::CalculateFactorGroup(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
  SYM::CalculatePointGroupCrystal(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
}

void xstructure::CalculateSymmetryPointGroupCrystal(void) {
  LatticeReduction_avoid=TRUE;
  if(pgroup_calculated==FALSE) CalculateSymmetryPointGroup(FALSE);
  if(fgroup_calculated==FALSE) CalculateSymmetryFactorGroup(FALSE);
  CalculateSymmetryPointGroupCrystal(FALSE);
}

void CalculateSymmetryPointGroupCrystal(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  str.LatticeReduction_avoid=TRUE;
  if(str.pgroup_calculated==FALSE) SYM::CalculatePointGroup(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
  if(str.fgroup_calculated==FALSE) SYM::CalculateFactorGroup(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
  SYM::CalculatePointGroupCrystal(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
}

void CalculateSymmetryPointGroupCrystal(xstructure& str,bool ossverbose,ostream& oss) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroupCrystal(str,ossverbose,oss,FALSE);
}

void CalculateSymmetryPointGroupCrystal(xstructure& str,bool ossverbose) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroupCrystal(str,ossverbose,cout,FALSE);
}

void CalculateSymmetryPointGroupCrystal(xstructure& str) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  if(str.pgroup_calculated==FALSE) CalculateSymmetryPointGroup(str,FALSE,cout,FALSE);
  if(str.fgroup_calculated==FALSE) CalculateSymmetryFactorGroup(str,FALSE,cout,FALSE);
  CalculateSymmetryPointGroupCrystal(str,FALSE,cout,FALSE);
}

//DX - Consider using pflow::CalculateFullSymmetry in aflow_aconvasp_main.cpp.
//      It contains consistency checks for the symmetry analysis.
// ***************************************************************************
// Function CalculateSymmetryFactorGroup
// ***************************************************************************
void xstructure::CalculateSymmetryFactorGroup(bool ossverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  LatticeReduction_avoid=TRUE;
  SYM::CalculatePointGroup(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
  SYM::CalculateFactorGroup(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
}

void xstructure::CalculateSymmetryFactorGroup(void) {
  LatticeReduction_avoid=TRUE;
  CalculateSymmetryFactorGroup(FALSE);
}

void CalculateSymmetryFactorGroup(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  str.LatticeReduction_avoid=TRUE;
  SYM::CalculatePointGroup(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
  SYM::CalculateFactorGroup(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
}

void CalculateSymmetryFactorGroup(xstructure& str,bool ossverbose,ostream& oss) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryFactorGroup(str,ossverbose,oss,FALSE);
}

void CalculateSymmetryFactorGroup(xstructure& str,bool ossverbose) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryFactorGroup(str,ossverbose,cout,FALSE);
}

void CalculateSymmetryFactorGroup(xstructure& str) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryFactorGroup(str,FALSE,cout,FALSE);
}

//DX - Consider using pflow::CalculateFullSymmetry in aflow_aconvasp_main.cpp.
//      It contains consistency checks for the symmetry analysis.
//ME20200114 - made capitalization more consistent with other functions
// ***************************************************************************
// Function CalculateSymmetryPointGroupKLattice
// ***************************************************************************
void xstructure::CalculateSymmetryPointGroupKLattice(bool ossverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  // LatticeReduction_avoid=TRUE;  // not necssary for klattice
  SYM::CalculatePointGroupKLattice(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
}

void xstructure::CalculateSymmetryPointGroupKLattice(void) {
  // LatticeReduction_avoid=TRUE;  // not necssary for klattice
  CalculateSymmetryPointGroupKLattice(FALSE);
}

void CalculateSymmetryPointGroupKLattice(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  // str.LatticeReduction_avoid=TRUE;  // not necssary for klattice
  SYM::CalculatePointGroupKLattice(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
}

void CalculateSymmetryPointGroupKLattice(xstructure& str,bool ossverbose,ostream& oss) {
  // str.LatticeReduction_avoid=TRUE;  // not necssary for klattice
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroupKLattice(str,ossverbose,oss,FALSE);
}

void CalculateSymmetryPointGroupKLattice(xstructure& str,bool ossverbose) {
  //  str.LatticeReduction_avoid=TRUE;  // not necssary for klattice
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroupKLattice(str,ossverbose,cout,FALSE);
}

void CalculateSymmetryPointGroupKLattice(xstructure& str) {
  //  str.LatticeReduction_avoid=TRUE;  // not necssary for klattice
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroupKLattice(str,FALSE,cout,FALSE);
}

//ME20200114 - added missing function
// ***************************************************************************
// Function CalculateSymmetryPointGroupKCrystal
// ***************************************************************************
void xstructure::CalculateSymmetryPointGroupKCrystal(bool ossverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  SYM::CalculatePointGroupKCrystal(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
}

void xstructure::CalculateSymmetryPointGroupKCrystal(void) {
  CalculateSymmetryPointGroupKCrystal(FALSE);
}

void CalculateSymmetryPointGroupKCrystal(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  SYM::CalculatePointGroupKCrystal(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
}

void CalculateSymmetryPointGroupKCrystal(xstructure& str,bool ossverbose,ostream& oss) {
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroupKCrystal(str,ossverbose,oss,FALSE);
}

void CalculateSymmetryPointGroupKCrystal(xstructure& str,bool ossverbose) {
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroupKCrystal(str,ossverbose,cout,FALSE);
}

void CalculateSymmetryPointGroupKCrystal(xstructure& str) {
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroupKCrystal(str,FALSE,cout,FALSE);
}

//ME20200129
// ***************************************************************************
// Function CalculateSymmetryPointGroupKPatterson
// ***************************************************************************

void xstructure::CalculateSymmetryPointGroupKPatterson(bool ossverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  SYM::CalculatePointGroupKPatterson(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
}

void xstructure::CalculateSymmetryPointGroupKPatterson(void) {
  CalculateSymmetryPointGroupKPatterson(FALSE);
}

void CalculateSymmetryPointGroupKPatterson(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  SYM::CalculatePointGroupKPatterson(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
}

void CalculateSymmetryPointGroupKPatterson(xstructure& str,bool ossverbose,ostream& oss) {
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroupKPatterson(str,ossverbose,oss,FALSE);
}

void CalculateSymmetryPointGroupKPatterson(xstructure& str,bool ossverbose) {
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroupKPatterson(str,ossverbose,cout,FALSE);
}

void CalculateSymmetryPointGroupKPatterson(xstructure& str) {
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroupKPatterson(str,FALSE,cout,FALSE);
}

// ***************************************************************************
// Function fixEmptyAtomNames()
// ***************************************************************************
void xstructure::fixEmptyAtomNames(bool force_fix){
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  if(species.size()==species_pp.size()) { //CO20190218
    for(uint itype=0;itype<species.size();itype++) {
      if((force_fix || species.at(itype)=="") && species_pp.at(itype)!=""){
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " species_pp.at(" << itype << ")=" << species_pp.at(itype) << endl;}
        species.at(itype)=species_pp.at(itype); //KBIN::VASP_PseudoPotential_CleanName(species_pp.at(itype));  //CO20181226 KEEP PP INFO if available (auto aflow.in)
      }
    }
  }  // cormac I`ll write a short pflow for this stuff
  if(species.size()==num_each_type.size()){
    int iatom=0;
    for(uint itype=0;itype<num_each_type.size();itype++) {
      string s=string(species.at(itype));
      species.at(itype)=s;
      for(int j=0;j<num_each_type.at(itype);j++) {
        atoms.at(iatom).name=s;    // CONVASP_MODE
        atoms.at(iatom).CleanName();
        atoms.at(iatom).CleanSpin();
        atoms.at(iatom).name_is_given=TRUE;
        iatom++;
      }
    }
  }
}

// ***************************************************************************
// Function buildGenericTitle()
// ***************************************************************************
void xstructure::buildGenericTitle(bool vasp_input,bool force_fix){
  if(vasp_input){pflow::fixEmptyAtomNames(*this,force_fix);}
  title.clear();
  title=getGenericTitleXStructure(*this,false); //no latex
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]uint iat=0;
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]//if any names missing from atoms, lets use generic names
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]bool atom_names=true;
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]for(uint i=0;i<atoms.size()&&atom_names;i++){if(atoms[i].name.empty()){atom_names=false;}} //CO20180316 - use pp names
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]for(uint itype=0;itype<num_each_type.size();itype++){
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]  for(uint j=0;j<(uint) num_each_type.at(itype);j++) {
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]    if(j==0){
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]      if(atom_names){title+=atoms.at(iat).name;} //CO20180316 - use pp names
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]      else {title+=char('A'+itype);}
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]      title+=aurostd::utype2string(num_each_type.at(itype));
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]    }
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]    iat++;
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]  }
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]}
}

// ***************************************************************************
// Function xstructure2qe
// ***************************************************************************
void xstructure::xstructure2qe(void) {
  ReScale(1.0);
  neg_scale=FALSE;
  coord_flag=_COORDS_FRACTIONAL_; 
  iomode=IOQE_GEOM;
  if(title.empty()) {buildGenericTitle();}  //CO20171008 - pushed all of this to a function
  return;
}

// ***************************************************************************
// Function xstructure2vasp
// ***************************************************************************
void xstructure::xstructure2vasp(void) {
  //DX20180425 [OBSOLETE] scale=1.0;
  ReScale(1.0); //DX20180425 - needs to be generic and rescale to 1.0; see other xstructure converters
  neg_scale=FALSE;
  coord_flag=_COORDS_FRACTIONAL_; 
  iomode=IOVASP_AUTO;
  //  cerr << "title=\"" << title << "\"" << endl;
  if(title.empty()) {buildGenericTitle();}  //CO20171008 - pushed all of this to a function
  //  uint iat=0;
  //  for(uint itype=0;itype<num_each_type.size();itype++) 
  //    for(uint j=0;j<(uint) num_each_type.at(itype);j++) {
  //      if(j==0) title+=atoms.at(iat).name+aurostd::utype2string(num_each_type.at(itype));
  //      //      atoms.at(iat++).type=itype;
  //    }
  return;
}

// ***************************************************************************
// Function xstructure2itc
// ***************************************************************************
void xstructure::xstructure2itc(void) { //CO20220613
  ReScale(1.0);
  neg_scale=FALSE;
  coord_flag=_COORDS_FRACTIONAL_;
  char iomode_orig=iomode;  //save
  iomode=IOVASP_WYCKCAR;
  if(title.empty()) {buildGenericTitle();}  //CO20171008 - pushed all of this to a function
  if(!partial_occupation_flag){ //CO20220715
    (*this).spacegroupnumber = (*this).SpaceGroup_ITC();
    (*this).lattice = (*this).standard_lattice_ITC; // need to update the lattice; may have rotated
  }
  stringstream ss;
  ss << (*this);
  (*this).clear();
  ss >> (*this);
  iomode=iomode_orig;
  return;
}

// ***************************************************************************
// Function xstructure2aims
// ***************************************************************************
void xstructure::xstructure2aims(void) {
  ReScale(1.0);
  neg_scale=FALSE;
  coord_flag=_COORDS_FRACTIONAL_; 
  iomode=IOAIMS_GEOM;
  if(title.empty()) {buildGenericTitle();}  //CO20171008
  return;
}

// ***************************************************************************
// Function xstructure2abinit
// ***************************************************************************
void xstructure::xstructure2abinit(void) {
  ReScale(1.0);
  neg_scale=FALSE;
  coord_flag=_COORDS_FRACTIONAL_; 
  iomode=IOABINIT_GEOM;
  if(title.empty()) {buildGenericTitle();}  //CO20171008 - pushed all of this to a function
  return;
}

// ***************************************************************************
// Function xstructure2cif
// ***************************************************************************
void xstructure::xstructure2cif(void) { //DX20190131
  ReScale(1.0);
  neg_scale=FALSE;
  coord_flag=_COORDS_FRACTIONAL_; 
  iomode=IOCIF;
  if(!partial_occupation_flag){ //CO20220715
    (*this).spacegroupnumber = (*this).SpaceGroup_ITC();
    (*this).lattice = (*this).standard_lattice_ITC; // need to update the lattice; may have rotated
  }
  if(title.empty()) {buildGenericTitle();}  //CO20171008 - pushed all of this to a function
  return;
}

// ***************************************************************************
// Function xstructure2abccar
// ***************************************************************************
void xstructure::xstructure2abccar(void) { //DX20190131
  ReScale(1.0);
  neg_scale=FALSE;
  coord_flag=_COORDS_FRACTIONAL_; 
  iomode=IOVASP_ABCCAR;
  if(title.empty()) {buildGenericTitle();}  //CO20171008 - pushed all of this to a function
  return;
}

// ***************************************************************************
// Function xstructure2elk
// ***************************************************************************
void xstructure::xstructure2elk(void) { //DX20200313
  ReScale(1.0);
  neg_scale=FALSE;
  coord_flag=_COORDS_FRACTIONAL_; 
  iomode=IOELK_GEOM;
  if(title.empty()) {buildGenericTitle();}  //CO20171008 - pushed all of this to a function
  return;
}

// ***************************************************************************
// Function xstructure2atat
// ***************************************************************************
void xstructure::xstructure2atat(void) { //SD20220123
  ReScale(1.0);
  neg_scale=FALSE;
  coord_flag=_COORDS_FRACTIONAL_;
  iomode=IOATAT_STR;
  if(title.empty()) {buildGenericTitle();}  //CO20171008 - pushed all of this to a function
  return;
}

// ***************************************************************************
// Function platon2print
// ***************************************************************************
string xstructure::platon2print(bool P_EQUAL,bool P_EXACT,double P_ang,double P_d1,double P_d2,double P_d3) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream oss;
  oss.setf(std::ios::fixed,std::ios::floatfield);
  xstructure str=*this;

  // Deal with too small volume problem
  if(str.GetVolume() < (double) PLATON_MIN_VOLUME_PER_ATOM*str.atoms.size()) {
    str.SetVolume(PLATON_MIN_VOLUME_PER_ATOM*str.atoms.size()) ;
    if(LDEBUG) cerr << "platon2print:  PLATON FIXED VOLUME PER ATOM = " << PLATON_MIN_VOLUME_PER_ATOM << endl;
  }
  if(LDEBUG) cerr << "platon2print: volume=" << str.GetVolume() << endl;
  // Set scale to 1 so you don't need to rescale coordinates.
  str.ReScale(1.0);
  std::vector<string> el_names(7);
  int k;
  el_names[0]="Ag";
  el_names[1]="Zr";
  el_names[2]="Cd";
  el_names[3]="Mo";
  el_names[4]="Fe";
  el_names[5]="W";
  el_names[6]="O";
  // Print out data in Platon format

  oss << "TITL " << str.title << endl;
  oss.precision(8);
  oss << "CELL "
    << str.scale*aurostd::modulus(str.lattice(1,1),str.lattice(1,2),str.lattice(1,3)) << " "
    << str.scale*aurostd::modulus(str.lattice(2,1),str.lattice(2,2),str.lattice(2,3)) << " "
    << str.scale*aurostd::modulus(str.lattice(3,1),str.lattice(3,2),str.lattice(3,3)) << " "
    << aurostd::angle(str.lattice(2,1),str.lattice(2,2),str.lattice(2,3),str.lattice(3,1),str.lattice(3,2),str.lattice(3,3)) << " "
    << aurostd::angle(str.lattice(1,1),str.lattice(1,2),str.lattice(1,3),str.lattice(3,1),str.lattice(3,2),str.lattice(3,3)) << " "
    << aurostd::angle(str.lattice(1,1),str.lattice(1,2),str.lattice(1,3),str.lattice(2,1),str.lattice(2,2),str.lattice(2,3))
    << endl;
  // oss << str.atoms.size() << endl;
  oss.precision(10);
  k=0;
  if(str.num_each_type.size()!=2) {
    for(uint i=0;i<str.num_each_type.size();i++)
      for(int j=1;j<=str.num_each_type.at(i);j++) {
        //      printf("%c%i ",i+66,j);
        if(str.atoms.at(k).name_is_given) {
          oss << str.atoms.at(k).cleanname << j << " ";
        }
        else { // Must make up a name
          if(i>(el_names.size()-1)) { // No more default names so make everything W
            oss << "W" << j << " ";
          }
          else { // Use default names
            oss << el_names[i] << j << " ";
          }
        }
        oss << str.atoms.at(k).fpos(1)<< " " << str.atoms.at(k).fpos(2) << " " << str.atoms.at(k).fpos(3) << endl;
        k++;
      }
  }
  if(str.num_each_type.size()==2) {
    for(uint i=0;i<str.num_each_type.size();i++)
      for(int j=1;j<=str.num_each_type.at(i);j++) {
        if(str.num_each_type.at(0)>str.num_each_type.at(1))
          oss << el_names[1-i] << j << " " << str.atoms.at(k).fpos(1)<< " " << str.atoms.at(k).fpos(2) << " " << str.atoms.at(k).fpos(3) << endl;
        else
          oss << el_names[i] << j << " " << str.atoms.at(k).fpos(1)<< " " << str.atoms.at(k).fpos(2) << " " << str.atoms.at(k).fpos(3) << endl;
        k++;
      }
  }
  // oss <<  num_each_type.at(i) << " ";
  //  int str.atoms.size()=cpos.size();
  // for(i=0;i<str.atoms.size();i++)
  //  if(str.coord_flag==0)
  // oss << " " << str.atoms.at(i).fpos(1)<< " " << str.atoms.at(i).fpos(2) << " " << stry.atoms.at(i).fpos(3) << endl;
  //  printf("CALC ADDSYM ");
  oss << "CALC ADDSYM ";// << endl;
  if(P_EQUAL) oss << "EQUAL ";
  if(P_EXACT) oss << "EXACT ";
  if(!aurostd::isequal(P_ang,PLATON_TOLERANCE_ANGLE,1.0e-6) ||
      !aurostd::isequal(P_d1,PLATON_TOLERANCE_D1,1.0e-6) ||
      !aurostd::isequal(P_d1,PLATON_TOLERANCE_D2,1.0e-6) ||
      !aurostd::isequal(P_d1,PLATON_TOLERANCE_D3,1.0e-6))
    //    printf("%14.10f %14.10f %14.10f %14.10f ",P_ang,P_d1,P_d2,P_d3);  printf("\n");
    oss << P_ang << " " << P_d1 << " " << P_d2 << " " << P_d3 << " ";//
  oss << endl;
  return oss.str();
}

// ***************************************************************************
// Function DecorateWithElements()
// ***************************************************************************
void xstructure::DecorateWithElements(void) {

  // Apply an element to each atom type.
  // Elements are first alphabetized to follow the AFLOW convention


  // elements need to be alphabetic for AFLOW
  deque<string> elements;
  for(uint i=0;i<velement.size();i++){ elements.push_back(velement[i].symbol); } // from xelement
  std::stable_sort(elements.begin(), elements.end());

  // update species and atom names;
  SetSpecies(elements);

  //DX20200727 [OBSOLETE] int iatom=0;
  //DX20200727 [OBSOLETE] for(uint itype=0;itype<num_each_type.size();itype++)
  //DX20200727 [OBSOLETE]   for(int j=0;j<num_each_type.at(itype);j++) {
  //DX20200727 [OBSOLETE]     //    if(atoms.at(iatom).name_is_given==FALSE)
  //DX20200727 [OBSOLETE]     {
  //DX20200727 [OBSOLETE]       atoms.at(iatom).name=vatom_symbol[2+atoms.at(iatom).type];
  //DX20200727 [OBSOLETE]       if(atoms.at(iatom).type==0) atoms.at(iatom).name="Ag"; // works....
  //DX20200727 [OBSOLETE]       //	if(atoms.at(iatom).type==1) atoms.at(iatom).name="Au"; // works....
  //DX20200727 [OBSOLETE]       atoms.at(iatom).CleanName();
  //DX20200727 [OBSOLETE]       //DX20170921 - Need to keep spin info atoms.at(iatom).CleanSpin();
  //DX20200727 [OBSOLETE]       atoms.at(iatom).name_is_given=TRUE;
  //DX20200727 [OBSOLETE]     }
  //DX20200727 [OBSOLETE]     iatom++;
  //DX20200727 [OBSOLETE]   }
}

// ***************************************************************************
// xstructure::DecorateWithFakeElements() //DX20200724
// ***************************************************************************
void xstructure::DecorateWithFakeElements(){

  // Apply a fake letter to each atom type.
  // Using letters (not elements) to avoid confusion with real materials
  // i.e., prototype vs material
  // In the case of compounds with more
  // than 26 species it is necessary to add more characters to this string


  // get fake elements
  vector<string> fake_elements = pflow::getFakeElements(num_each_type.size());

  // update species atom names;
  SetSpecies(aurostd::vector2deque(fake_elements));
}

// ***************************************************************************
// Function SpaceGroup
// ***************************************************************************
string xstructure::platon2sg(bool P_EQUAL,bool P_EXACT,double P_ang,double P_d1,double P_d2,double P_d3) {
  xstructure str=*this;
  stringstream aus;
  // string directory="/tmp/_aflow_platon_"+XHOST.ostrPID.str()+"_"+XHOST.ostrTID.str();  // dont change what works //CO20200502 - threadID
  string directory=XHOST.tmpfs+"/_aflow_platon_"+XHOST.ostrPID.str()+"_"+XHOST.ostrTID.str();  // dont change what works  //CO20200502 - threadID
  string file=directory+"/aflow_platon_";
  string file_spf=file+".spf";
  string file_out=file+".out";
  string output;
  vector<string> space_group;
  aurostd::DirectoryMake(directory);
  str.DecorateWithElements(); //DX20200727 - FakeNames() -> DecorateWithElements();
  aus << str.platon2print(P_EQUAL,P_EXACT,P_ang,P_d1,P_d2,P_d3);
  aurostd::stringstream2file(aus,file_spf);aus.clear();aus.str(std::string());
  aus << "cd " << directory << endl;
  aus << "platon -o -c " << file_spf << " | grep \"Space Group\" | head -1 | sed \"s/Space Group //g\"  | sed \"s/No:/#/g\" | sed \"s/,/#/g\" | sed \"s/ //g\" > " << file_out << endl;
  aurostd::execute(aus);
  aurostd::string2tokens(aurostd::file2string(file_out),space_group,"#");
  aus.clear();aus.str(std::string());
  aus << "rm -rf " << directory << endl;
  aurostd::execute(aus);
  spacegroup=space_group[0]+" #"+space_group[1];
  spacegrouplabel=space_group[0];
  spacegroupnumber=aurostd::string2utype<int>(space_group[1]);
  spacegroupoption="FIX: to extract from platon";
  is_spacegroup_platon=TRUE;
  return spacegroup;
}

string xstructure::findsym2sg(double tolerance) {
  xstructure str=*this;
  // Read in input file.
  stringstream oss;
  string findsym_in=aurostd::TmpFileCreate("findsym.in");
  oss << str.findsym2print(tolerance);
  aurostd::stringstream2file(oss,findsym_in);
  vector<string> vlines,tokens;
  FINDSYM::Write("data_space.txt","./");
  FINDSYM::Write("data_wyckoff.txt","./");
  aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("findsym")+" < "+findsym_in),vlines);
  FROZSL::Delete("data_space.txt","./");
  FROZSL::Delete("data_wyckoff.txt","./");
  aurostd::RemoveFile(findsym_in);
  aurostd::RemoveFile("./findsym.log");
  for(uint i=0;i<vlines.size();i++) {
    if(aurostd::substring2bool(vlines.at(i),"Space Group"))
      aurostd::string2tokens(vlines.at(i),tokens);
  }
  //    for(uint i=0;i<tokens.size();i++) cerr << tokens.at(i) << endl;
  spacegroup=tokens.at(4)+" #"+tokens.at(2);
  spacegrouplabel=tokens.at(4);
  spacegroupnumber=aurostd::string2utype<int>(tokens.at(2));
  spacegroupoption="FIX: to extract from findsym";
  is_spacegroup_findsym=TRUE;
  return spacegroup;
}

// ***************************************************************************
// xstructure::findsym2execute
// ***************************************************************************
// This funtion executes findsys
string xstructure::findsym2execute(double tolerance) {
  xstructure str=*this;
  // Read in input file.
  stringstream oss;
  string findsym_in=aurostd::TmpFileCreate("findsym.in");
  oss << str.findsym2print(tolerance);
  aurostd::stringstream2file(oss,findsym_in);
  vector<string> vlines,tokens;
  FINDSYM::Write("data_space.txt","./");
  FINDSYM::Write("data_wyckoff.txt","./");
  string out=aurostd::execute2string(XHOST.command("findsym")+" < "+findsym_in);
  FROZSL::Delete("data_space.txt","./");
  FROZSL::Delete("data_wyckoff.txt","./");
  aurostd::RemoveFile(findsym_in);
  aurostd::RemoveFile("./findsym.log");
  return out;
}

// ***************************************************************************
// xstructure::findsym2print
// ***************************************************************************
// This funtion prints out structural data in a format for findsys
string xstructure::findsym2print(double tolerance) {
  xstructure sstr=*this;
  stringstream oss;
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(10);
  // Set scale to 1 so you don't need to rescale coordinates.
  sstr.ReScale(1.0);
  // Print out data in Findsym format
  oss << sstr.title << endl;
  oss << tolerance << "                               accuracy (0=highest=1e-6) " << endl;
  if(0) { // old
    oss << "1                                          form lattice: (1=3 vectors) (2=lengths,angles) " << endl;
    oss << sstr.lattice(1,1) << " " << sstr.lattice(1,2) << " " << sstr.lattice(1,3) << "    Bravais lattice" << endl;
    oss << sstr.lattice(2,1) << " " << sstr.lattice(2,2) << " " << sstr.lattice(2,3) << "    Bravais lattice" << endl;
    oss << sstr.lattice(3,1) << " " << sstr.lattice(3,2) << " " << sstr.lattice(3,3) << "    Bravais lattice" << endl;
    oss << "1.0 0.0 0.0                   unit cell in function of lattice parameter" << endl;
    oss << "0.0 1.0 0.0                   unit cell in function of lattice parameter" << endl;
    oss << "0.0 0.0 1.0                   unit cell in function of lattice parameter" << endl;
    oss << sstr.atoms.size() << "                             number of atoms in the primitive unit cell" << endl;
    for(uint i=0;i<sstr.num_each_type.size();i++)
      for(int j=1;j<=sstr.num_each_type.at(i);j++)
        oss <<  sstr.num_each_type.at(i) << " ";
  }
  if(1) { //new
    oss << "2 Form of lattice parameters: to be entered as lengths and angles " << endl;
    oss << sstr.a << " " << sstr.b << " " << sstr.c << " " << sstr.alpha << " " << sstr.beta << " " << sstr.gamma << " a, b, c, alpha, beta, gamma " << endl;
    oss << "1 Vectors defining the unit cell" << endl;
    oss << "1.0 0.0 0.0                   unit cell in function of lattice parameter" << endl;
    oss << "0.0 1.0 0.0                   unit cell in function of lattice parameter" << endl;
    oss << "0.0 0.0 1.0                   unit cell in function of lattice parameter" << endl;
    oss << sstr.atoms.size() << "                            number of atoms in the primitive unit cell" << endl;
    for(uint i=0;i<sstr.atoms.size();i++)
      oss <<  sstr.atoms.at(i).type+1 << " ";
  }  
  oss << endl;
  for(uint i=0;i<sstr.atoms.size();i++) {
    //  if(coord_flag==0)
    oss << " " << sstr.atoms.at(i).fpos(1) << " " << sstr.atoms.at(i).fpos(2) << " " << sstr.atoms.at(i).fpos(3) << endl;
    //  if(coord_flag==1) oss << " " << sstr.atoms.at(i).cpos(1) << " " << sstr.atoms.at(i).cpos(2) << " " << sstr.atoms.at(i).cpos(3) << endl;
  }
  //  oss << endl;
  return oss.str();
}
// Stefano Curtarolo MIT 2002 - DUKE 2013

// ***************************************************************************
// Rotate
// ***************************************************************************
// Rotation is done around the origin of the structure.
// Define origin q, a point p, and a rotation of p around
// q (R_q(p)).  Then R_q(p)=R_0(p-q)+q=R_0(p)-R_0(q)+q.
// We can evalutate the R_0 terms by simply mulitplying by
// a matrix.  Then we must add q to all the final cartesian
// positions.

// ---------------------------------------------------------------------------
// returns new xstructure (makes a copy)
xstructure Rotate(const xstructure&a, const xmatrix<double>& rm) {
  xstructure xstr_rotated = a;
  xstr_rotated.Rotate(rm);
  return xstr_rotated;
}

// ---------------------------------------------------------------------------
// modifies in-place (efficient)
void xstructure::Rotate(const xmatrix<double>& rm) {
  bool LDEBUG=(FALSE || XHOST.DEBUG); //CO20190520
  if(LDEBUG) { //CO20190520
    cerr << __AFLOW_FUNC__ << " (*this)=" << endl;cerr << (*this) << endl; //CO20190520
    cerr << __AFLOW_FUNC__ << " (*this).origin=" << (*this).origin << endl; //CO20190520
    cerr << __AFLOW_FUNC__ << " rm=" << endl;cerr << rm << endl; //ME20200204
  }
  if (aurostd::isidentity(rm)) return;  //ME20200204 - no need to go through all the motions for identity matrix
  // Get R_0(p) for all cartesian positions.
  xmatrix<double> nlattice(3,3);
  nlattice=trasp((*this).lattice);
  (*this).lattice=trasp(rm*nlattice);
  (*this).FixLattices();  //CO20190409 - so we don't need to keep redefining f2c/c2f
  const xmatrix<double>& f2c=(*this).f2c; //CO20190520
  const xmatrix<double>& c2f=(*this).c2f; //CO20190520
  uint natoms = (*this).atoms.size(); //DX+ME20210111 - set variable to optimize for-loops
  //DX+ME20210111 [OBSOLETE - consolidate into one for-loop, below] for(uint ia=0;ia<natoms;ia++){(*this).atoms[ia].cpos=f2c*(*this).atoms[ia].fpos;}  //CO20190409 - so we don't need to keep redefining f2c/c2f
  //[CO20190409 - OBSOLETE]for(int ia=0;ia<(int)b.atoms.size();ia++)
  //[CO20190409 - OBSOLETE]  b.atoms.at(ia).cpos=F2C(b.lattice,b.atoms.at(ia).fpos);
  //Get R_0(q)
  xvector<double> r_orig(3);
  r_orig=rm*(*this).origin;
  // Assign new cartesian positions
  //DX+ME20210111 [OBSOLETE - consolidate into one for-loop, below] for(uint ia=0;ia<natoms;ia++){(*this).atoms[ia].cpos+=(-r_orig+(*this).origin);}
  // Get all the direct coords.
  //DX+ME20210111 [OBSOLETE - consolidate into one for-loop, below] for(uint ia=0;ia<natoms;ia++){(*this).atoms[ia].fpos=c2f*(*this).atoms[ia].cpos;}  //CO20190409 - so we don't need to keep redefining f2c/c2f
  for(uint ia=0;ia<natoms;ia++){ //DX+ME20210111 consolidate into one for-loop
    (*this).atoms[ia].cpos=f2c*(*this).atoms[ia].fpos;
    (*this).atoms[ia].cpos+=(-r_orig+(*this).origin);
    (*this).atoms[ia].fpos=c2f*(*this).atoms[ia].cpos;
  }
  //[CO20190409 - OBSOLETE]for(int ia=0;ia<(int)b.atoms.size();ia++)
  //[CO20190409 - OBSOLETE]  b.atoms.at(ia).fpos=C2F(b.lattice,b.atoms.at(ia).cpos);
  return;
}

// ***************************************************************************
// GetRotationMatrix // moved from pflow
// ***************************************************************************
// This gets a rotation matrix from 3 angles assumed
// to represent a rotation around x, then y, then z.
// Angles are assumed to be in radians.
aurostd::matrix<double> GetRotationMatrix(const vector<double>& angles) {  //CO20200404 pflow::matrix()->aurostd::matrix()
  // Sin and cos.
  vector<double> sn(3,0.0);
  vector<double> cs(3,0.0);
  for(int ic=0;ic<3;ic++) {
    sn[ic]=sin(angles[ic]);
    cs[ic]=cos(angles[ic]);
  }
  // Set rotation matrix (do x, then y, then z rotation).
  aurostd::matrix<double> xm(3,3);pflow::VVset(xm,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
  aurostd::matrix<double> ym(3,3);pflow::VVset(ym,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
  aurostd::matrix<double> zm(3,3);pflow::VVset(zm,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()

  xm[0][0]=1;
  xm[1][1]=cs[0];
  xm[1][2]=-sn[0];
  xm[2][1]=sn[0];
  xm[2][2]=cs[0];

  ym[0][0]=cs[1];
  ym[0][2]=sn[1];
  ym[1][1]=1;
  ym[2][0]=-sn[1];
  ym[2][2]=cs[1];

  zm[0][0]=cs[2];
  zm[0][1]=-sn[2];
  zm[1][0]=sn[2];
  zm[1][1]=cs[2];
  zm[2][2]=1;

  aurostd::matrix<double> rm;  //CO20200404 pflow::matrix()->aurostd::matrix()
  rm=pflow::MMmult(zm,pflow::MMmult(ym,xm));
  return rm;
}

// ***************************************************************************
// RotateStrVec // moved from pflow
// ***************************************************************************
// This rotates each structure in the vstr.
// The rotation goes from an initial to a final set of angles
// in steps of the primitive rotation.
// The primitive rotation is the rotation around x, y, then z
// by the amount of change given in rot divided by the number
// of structures - 1.  All frames are initially rotated
// according to the initial rotation angles.  The rotation is
// then done as one primitive rotation per frame.
void RotateStrVec(vector<xstructure>& vstr, const vector<double>& rot) {
  // Get initial rotation matrix.
  vector<double> angles(3);
  for(int ic=0;ic<3;ic++) {
    angles[ic]=rot[2*ic];
    angles[ic]=angles[ic]*TWOPI/360.0;
  }
  aurostd::matrix<double> irm=GetRotationMatrix(angles); //CO20200404 pflow::matrix()->aurostd::matrix()

  // get primitive rotation matrix.
  int s=vstr.size()-1;
  if(s<1) s=1;
  for(int ic=0;ic<3;ic++) {
    angles[ic]=(rot[2*ic+1]-rot[2*ic])/(s);
    angles[ic]=angles[ic]*TWOPI/360.0;
  }
  aurostd::matrix<double> prm=GetRotationMatrix(angles); //CO20200404 pflow::matrix()->aurostd::matrix()
  aurostd::matrix<double> rm=irm;  //CO20200404 pflow::matrix()->aurostd::matrix()
  xmatrix<double> xprm(3,3); xprm=aurostd::matrix2xmatrix(prm);  //CO20200404 pflow::matrix()->aurostd::matrix()
  xmatrix<double> xrm(3,3);  xrm=aurostd::matrix2xmatrix(rm);  //CO20200404 pflow::matrix()->aurostd::matrix()
  for(int is=0;is<(int)vstr.size();is++) {
    //    xrm=aurostd::matrix2xmatrix(rm); //CO20200404 pflow::matrix()->aurostd::matrix()
    //   vstr[is]=Rotate(vstr[is],xrm);
    //DX20210127 [OBSOLETE, avoid copying xstructure] vstr[is]=Rotate(vstr[is],aurostd::matrix2xmatrix(rm)); //CO20200404 pflow::matrix()->aurostd::matrix()
    vstr[is].Rotate(aurostd::matrix2xmatrix(rm)); //CO20200404 pflow::matrix()->aurostd::matrix() //DX20210127 - do not make copy of xstructure
    rm=pflow::MMmult(prm,rm);
  }
}

// ***************************************************************************
//  Function GetLTCell
// ***************************************************************************
// This function operates a linear transform on a cell.
xstructure GetLTCell(const xmatrix<double>& lt, const xstructure& str) {
  xstructure nstr(str);
  nstr.lattice=str.lattice*trasp(lt);
  for(int iat=0;iat<(int)nstr.atoms.size();iat++) {
    //    nstr.atoms.at(iat).cpos=str.atoms.at(iat).cpos*trasp(lt);
    //  nstr.atoms.at(iat).fpos=C2F(nstr.lattice,nstr.atoms.at(iat).cpos);
    nstr.atoms.at(iat)=F2C(nstr.lattice,nstr.atoms.at(iat));
  }
  return nstr;
}

// ***************************************************************************
//  Function GetLTFVCell
// ***************************************************************************
// This function operates a linear transform on a cell.  Rotates all vectors
// by phi around nvec.  Formula from Goldstein, 1980 (Actually used version on
// Mathematic web site).

xstructure GetLTFVCell(const xvector<double>& nvec, const double phi, const xstructure& str) {
  //  xmatrix<double> fpos=str.fpos;
  xmatrix<double> nlat(3,3);
  xvector<double> nhat(3),v1(3),v2(3),v3(3),r(3);
  nhat=nvec/aurostd::modulus(nvec);
  double tmp;
  double cphi=std::cos(phi);
  double sphi=std::sin(phi);
  for(int i=1;i<=3;i++) {
    r=str.lattice(i);
    double nhatpr=scalar_product(nhat,r);
    v1=cphi*r;
    tmp=nhatpr*(1.0-cphi);
    v2=tmp*nhat;
    v3=vector_product(r,nhat);
    v3=sphi*v3;
    for(int j=1;j<=3;j++) nlat(i,j)=nlat(i,j)+v1(j);
    for(int j=1;j<=3;j++) nlat(i,j)=nlat(i,j)+v2(j);
    for(int j=1;j<=3;j++) nlat(i,j)=nlat(i,j)+v3(j);
  }
  xstructure nstr=str;
  nstr.lattice=nlat;
  for(int iat=0;iat<(int)nstr.atoms.size();iat++) {
    nstr.atoms.at(iat).fpos=str.atoms.at(iat).fpos;// *trasp(lt);
    nstr.atoms.at(iat).cpos=F2C(nstr.lattice,nstr.atoms.at(iat).fpos);
  }
  return nstr;
}

// ***************************************************************************
//  Function ShiftPos ShiftCpos ShiftFpos
// ***************************************************************************
// This function shifts the position of the atoms !
// DX20201215 - added in-place methods

// make a copy
xstructure ShiftPos(const xstructure& a,const xvector<double>& shift, bool is_frac) { //DX20210113 - "int flag" to "bool is_frac"
  xstructure b(a);
  b.ShiftPos(shift, is_frac);
  return b;
}

// modify in-place //DX2021215
void xstructure::ShiftPos(const xvector<double>& shift, bool is_frac) {
  // Cartesian shift
  if(is_frac==FALSE) { (*this).ShiftCPos(shift); } //DX20210113 - use function, reduce code
  // Direct coords shift
  else { (*this).ShiftFPos(shift); } //DX20210113 - use function, reduce code
}

// make a copy
xstructure ShiftCPos(const xstructure& a,const xvector<double>& shift) {
  xstructure b(a);
  b.ShiftCPos(shift);
  return b;
}

// modify in-place //DX2021215
void xstructure::ShiftCPos(const xvector<double>& shift) {
  uint natoms = (*this).atoms.size(); //DX20210111 - make variable to optimize for-loops
  xmatrix<double> c2f=(*this).scale*inverse(trasp((*this).lattice)); //DX+ME20210111
  (*this).coord_flag=TRUE;
  for(uint ia=0;ia<natoms;ia++) {
    (*this).atoms[ia].cpos=(*this).atoms[ia].cpos+shift;
    (*this).atoms[ia].fpos=c2f*(*this).atoms[ia].cpos; //DX+ME20210111 - use pre-calculated c2f, optimize
  }
}

// make a copy
xstructure ShiftFPos(const xstructure& a,const xvector<double>& shift) {
  xstructure b(a);
  b.ShiftCPos(shift);
  return b;
}

// modify in-place //DX2021215
void xstructure::ShiftFPos(const xvector<double>& shift) {
  uint natoms = (*this).atoms.size(); //DX20210111 - make variable to optimize for-loops
  xmatrix<double> f2c=(*this).scale*trasp((*this).lattice); //DX+ME20210111
  (*this).coord_flag=FALSE;
  for(uint ia=0;ia<natoms;ia++) {
    (*this).atoms[ia].fpos=(*this).atoms[ia].fpos+shift;
    (*this).atoms[ia].cpos=f2c*(*this).atoms[ia].fpos; //DX+ME20210111 - use pre-calculated f2c, optimize
  }
}

// **************************************************************************
// MaxStructureLattice and MinStructureLattice
// **************************************************************************
double MaxStructureLattice(const xstructure& str) {
  return max(aurostd::modulus(str.lattice(1)),aurostd::modulus(str.lattice(2)),aurostd::modulus(str.lattice(3)));
}

double MinStructureLattice(const xstructure& str) {
  return min(aurostd::modulus(str.lattice(1)),aurostd::modulus(str.lattice(2)),aurostd::modulus(str.lattice(3)));
}

// **************************************************************************
// AtomCDisp
// **************************************************************************
// This function returns the Cartesian displacement vector at2-at1 between two atoms.
xvector<double> AtomCDisp(const _atom& at1, const _atom& at2) {
  xvector<double> diff(3);
  diff=at2.cpos-at1.cpos;
  return diff;
}

// **************************************************************************
// AtomDist
// **************************************************************************
// This function returns the distance between two atoms.
double AtomDist(const xstructure& str,const _atom& atom1,const _atom& atom2) {  // with structure
  xvector<double> pos1(3),pos2(3);
  // contains shift WITH the ijk lattice and origin of structure
  pos1=str.scale*(atom1.cpos+atom1.ijk(1)*str.lattice(1)+atom1.ijk(2)*str.lattice(2)+atom1.ijk(3)*str.lattice(3))+str.origin;
  pos2=str.scale*(atom2.cpos+atom2.ijk(1)*str.lattice(1)+atom2.ijk(2)*str.lattice(2)+atom2.ijk(3)*str.lattice(3))+str.origin;
  return aurostd::modulus(pos1-pos2);
}
double AtomDist(const _atom& at1, const _atom& at2) {  // without structure
  return aurostd::modulus(AtomCDisp(at1,at2));
}

// **************************************************************************
// SameAtom
// **************************************************************************
#define _atom_eps_tol_       0.01
bool SameAtom(const xstructure& str,const _atom& atom1,const _atom& atom2) {
  if(AtomDist(str,atom1,atom2)<_atom_eps_tol_ && atom1.type==atom2.type) return TRUE;
  return FALSE;
}
bool SameAtom(const _atom& atom1,const _atom& atom2) {
  if(AtomDist(atom1,atom2)<_atom_eps_tol_ && atom1.type==atom2.type) return TRUE;
  return FALSE;
}

// **************************************************************************
// DifferentAtom
// **************************************************************************
bool DifferentAtom(const xstructure& str,const _atom& atom1,const _atom& atom2) {
  if(AtomDist(str,atom1,atom2)>_atom_eps_tol_ || atom1.type!=atom2.type) return TRUE;
  return FALSE;
}

// **************************************************************************
// GetDistMatrix
// **************************************************************************
//CO20171024
xmatrix<double> GetDistMatrix(const xstructure& aa){
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  xstructure a(aa); a.ReScale(1.0);
  xstructure xstr_cluster=a;
  //int dim=aurostd::max(LatticeDimensionSphere(a.lattice,RadiusSphereLattice(a.lattice))); //OVERKILL
  GenerateGridAtoms(xstr_cluster,RadiusSphereLattice(a.lattice));//dim,dim,dim);            //OVERKILL 
  vector<vector<int> > atom_types;
  atom_types.push_back(vector<int>(0)); atom_types.back().push_back(0);
  bool found;
  uint indx_found;
  for(uint i=1;i<a.atoms.size();i++){
    found=false;
    for(uint j=0;j<atom_types.size()&&!found;j++){
      if(!atom_types[j].size()){
        cerr << "GetDistMatrix: atom_types indices are compromised!" << endl;
        xmatrix<double> dummy;  //dummy
        return dummy;
      }  //should never come up
      if(a.atoms[i].type==a.atoms[atom_types[j][0]].type){
        found=true;
        indx_found=j;
      }
    }
    if(!found){
      atom_types.push_back(vector<int>(0));
      indx_found=atom_types.size()-1;
    }
    atom_types[indx_found].push_back(i);
  }
  if(a.num_each_type.size()){
    if(atom_types.size()!=a.num_each_type.size()){
      cerr << "GetDistMatrix: odd count of atom types, it's probably poorly sorted" << endl;
      xmatrix<double> dummy;  //dummy
      return dummy;
    }
  }
  if(LDEBUG) {
    for(uint it1=0;it1<atom_types.size();it1++){cerr << "GetDistMatrix: atom_types[" << it1 <<"][0]=" << atom_types[it1][0] << endl;}
  }
  xmatrix<double> distsij(1,1,atom_types.size(),atom_types.size());
  for(uint it1=0;it1<atom_types.size();it1++){
    for(uint it2=0;it2<atom_types.size();it2++){
      distsij(it1+1,it2+1)=distsij(it2+1,it1+1)=AUROSTD_MAX_DOUBLE;
    }
  }
  uint atom1=0,atom2=0;
  uint it1=0,it2=0,ia1=0,ia2=0;
  double distij,min_dist;
  for(it1=0;it1<atom_types.size();it1++){      //type 1 (must be in this order)
    for(it2=it1;it2<atom_types.size();it2++){  //type 2 (must be in this order)
      if(LDEBUG) {cerr << "GetDistMatrix: finding min dist between itype=" << it1 << " and itype=" << it2 << endl;}
      min_dist=AUROSTD_MAX_DOUBLE;                                      //reset min_dist
      for(ia1=0;ia1<atom_types[it1].size();ia1++){                 //must go through all atoms of same type
        atom1=xstr_cluster.grid_atoms_pc2scMap[atom_types[it1][ia1]];   //get respective index in cluster
        for(ia2=0;ia2<atom_types[it2].size();ia2++){               //must go through all atoms of the same type
          for(atom2=0;atom2<(uint)xstr_cluster.grid_atoms_number;atom2++){ //go through all atoms of the cluster
            if(atom1!=atom2 && a.atoms[atom_types[it2][ia2]].type==xstr_cluster.grid_atoms[atom2].type){  //cannot be same index (dist=0), and must be the types we want
              distij=AtomDist(xstr_cluster.grid_atoms[atom1],xstr_cluster.grid_atoms[atom2]); //distance
              if(0&&LDEBUG){
                cerr << "GetDistMatrix: trying dist=" << distij << " between atom[" << atom1 << ",type="<<xstr_cluster.grid_atoms[atom1].type<<"] and atom[" << atom2 << ",type="<<xstr_cluster.grid_atoms[atom2].type<<"]" << endl;
              }
              if(distij<min_dist){  //grab min_dist
                distsij(it1+1,it2+1)=distsij(it2+1,it1+1)=distij;
                min_dist=distij;
                if(LDEBUG) {
                  cerr << "GetDistMatrix: found new min dist=" << min_dist << " between atom[" << atom1 << ",type="<<xstr_cluster.grid_atoms[atom1].type<<"] and atom[" << atom2 << ",type="<<xstr_cluster.grid_atoms[atom2].type<<"]" << endl;
                }
              }
            }
          }
        }
      }
    }
  }
  return distsij;
}

// **************************************************************************
// GetNBONDXX
// **************************************************************************
//CO20171024
vector<double> GetNBONDXX(const xstructure& a){
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  xmatrix<double> distsij=GetDistMatrix(a);
  vector<double> dists;
  for(uint it1=1;it1<(uint)distsij.rows+1;it1++){
    for(uint it2=it1;it2<(uint)distsij.cols+1;it2++){
      dists.push_back(distsij(it1,it2));
    }
  }
  if(LDEBUG) {
    cerr << "GetNBONDXX: DIST_MATRIX" << endl;
    cerr << distsij << endl;
    cerr << "GetNBONDXX: MIN_DIST=" << SYM::minimumDistance(a) << endl;
  }
  return dists;
}

// **************************************************************************
// GenerateGridAtoms
// **************************************************************************
// make grid of atoms !
int xstructure::GenerateGridAtoms(double radius) {return GenerateGridAtoms(LatticeDimensionSphere((*this),radius));} // radius is not normalized over the scale
int xstructure::GenerateGridAtoms(int d) {return GenerateGridAtoms(-d,d,-d,d,-d,d);}
int xstructure::GenerateGridAtoms(int d1,int d2,int d3) {return GenerateGridAtoms(-d1,d1,-d2,d2,-d3,d3);}
int xstructure::GenerateGridAtoms(const xvector<int>& dims) {return GenerateGridAtoms(-dims(1),dims(1),-dims(2),dims(2),-dims(3),dims(3));} //CO20200912
int xstructure::GenerateGridAtoms(int i1,int i2,int j1,int j2,int k1,int k2) {return GenerateGridAtoms_20191218(i1,i2,j1,j2,k1,k2);}
int xstructure::GenerateGridAtoms_20190520(int i1,int i2,int j1,int j2,int k1,int k2) { //DX20191218 - added date [ORIG]
  bool LDEBUG=(FALSE || XHOST.DEBUG); //CO20190520
  if(LDEBUG) { //CO20190520
    cerr << __AFLOW_FUNC__ << " str=" << endl;cerr << (*this) << endl; //CO20190520
    cerr << __AFLOW_FUNC__ << " i=" << i1 << ":" << i2 << endl; //CO20190520
    cerr << __AFLOW_FUNC__ << " j=" << j1 << ":" << j2 << endl; //CO20190520
    cerr << __AFLOW_FUNC__ << " k=" << k1 << ":" << k2 << endl; //CO20190520
  } //CO20190520
  // same scale as before
  grid_atoms.clear();
  grid_atoms_sc2pcMap.clear(); grid_atoms_pc2scMap.clear();
  _atom atom;
  BringInCell();  // are INCELL.
  //xvector<double> a1(3),a2(3),a3(3);                     // a1,a2,a3 are the rows of the lattice matrix
  //a1=lattice(1);a2=lattice(2);a3=lattice(3); // a1,a2,a3 are the rows of the lattice matrix
  const xvector<double>& a1=lattice(1);  //CO20190520 - no need to make copies
  const xvector<double>& a2=lattice(2);  //CO20190520 - no need to make copies
  const xvector<double>& a3=lattice(3);  //CO20190520 - no need to make copies

  for(uint iat=0;iat<atoms.size();iat++){
    grid_atoms.push_back(atoms[iat]);  // put first the unit cell ! //DX20190709 - at to [] = speed increase
    grid_atoms_pc2scMap.push_back(grid_atoms.size()-1); //CO20171025 
    grid_atoms_sc2pcMap.push_back(iat); //CO20171025
  }
  for(int i=i1;i<=i2;i++) {
    for(int j=j1;j<=j2;j++) {
      for(int k=k1;k<=k2;k++) {
        if(i!=0 || j!=0 || k!=0) {
          for(uint iat=0;iat<atoms.size();iat++) {
            atom=atoms[iat]; //DX20190709 - at to [] = speed increase
            atom.isincell=FALSE; // these are OUT OF CELL
            atom.cpos=((double)i)*a1+((double)j)*a2+((double)k)*a3+atoms[iat].cpos; //DX20190709 - at to [] = speed increase
            atom.fpos[1]=i+atoms[iat].fpos[1]; //DX20190709 - at to [] = speed increase
            atom.fpos[2]=j+atoms[iat].fpos[2]; //DX20190709 - at to [] = speed increase
            atom.fpos[3]=k+atoms[iat].fpos[3]; //DX20190709 - at to [] = speed increase
            grid_atoms.push_back(atom);
            grid_atoms_sc2pcMap.push_back(iat); //CO20171025
            if(LDEBUG) { //CO20190520
              cerr << __AFLOW_FUNC__ << " grid_atoms[" << grid_atoms.size()-1 << "].cpos=" << grid_atoms.back().cpos << endl; //CO20190520
              cerr << __AFLOW_FUNC__ << " grid_atoms[" << grid_atoms.size()-1 << "].fpos=" << grid_atoms.back().fpos << endl; //CO20190520
              cerr << __AFLOW_FUNC__ << " grid_atoms[" << grid_atoms.size()-1 << "]=" << grid_atoms.back() << endl; //DX20191218
              cerr << __AFLOW_FUNC__ << " grid_atoms_sc2pcMap[" << grid_atoms.size()-1 << "]=" << grid_atoms_sc2pcMap.back() << endl; //DX20191218
            } //CO20190520
          }
        }
      }
    }
  }
  if(0){  //CO20190808 - quick check of mindist
    double min_dist_local=AUROSTD_MAX_DOUBLE,min_dist=AUROSTD_MAX_DOUBLE;
    for(uint i=0;i<grid_atoms.size()-1;i++){
      for(uint j=i+1;j<grid_atoms.size();j++){
        min_dist_local=aurostd::modulus(grid_atoms[i].cpos-grid_atoms[j].cpos);
        if(min_dist_local<min_dist){
          min_dist=min_dist_local;
        }
      }
    }
    if(!aurostd::isequal(min_dist,SYM::minimumDistance((*this)),0.1)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Minimum distance changed, check that atoms are not rotated",_RUNTIME_ERROR_);}
  }
  grid_atoms_calculated=TRUE;
  grid_atoms_dimsL[1]=i1;grid_atoms_dimsL[2]=j1;grid_atoms_dimsL[3]=k1;
  grid_atoms_dimsH[1]=i2;grid_atoms_dimsH[2]=j2;grid_atoms_dimsH[3]=k2;
  grid_atoms_number=grid_atoms.size();
  //  for(uint i=0;i<grid_atoms.size();i++)
  //   cerr << grid_atoms.at(i).cpos << endl;
  // cerr << grid_atoms.size() << endl;
  return grid_atoms.size();
}

int xstructure::GenerateGridAtoms_20191218(int i1,int i2,int j1,int j2,int k1,int k2) { //DX20191218 - [NEW]
  bool LDEBUG=(FALSE || XHOST.DEBUG); //CO20190520
  if(LDEBUG) { //CO20190520
    cerr << __AFLOW_FUNC__ << " str=" << endl;cerr << (*this) << endl; //CO20190520
    cerr << __AFLOW_FUNC__ << " i=" << i1 << ":" << i2 << endl; //CO20190520
    cerr << __AFLOW_FUNC__ << " j=" << j1 << ":" << j2 << endl; //CO20190520
    cerr << __AFLOW_FUNC__ << " k=" << k1 << ":" << k2 << endl; //CO20190520
  } //CO20190520
  // same scale as before
  grid_atoms.clear();
  grid_atoms_sc2pcMap.clear(); grid_atoms_pc2scMap.clear();
  _atom atom;
  BringInCell();  // are INCELL.
  //xvector<double> a1(3),a2(3),a3(3);                     // a1,a2,a3 are the rows of the lattice matrix
  //a1=lattice(1);a2=lattice(2);a3=lattice(3); // a1,a2,a3 are the rows of the lattice matrix
  const xvector<double>& a1=lattice(1);  //CO20190520 - no need to make copies
  const xvector<double>& a2=lattice(2);  //CO20190520 - no need to make copies
  const xvector<double>& a3=lattice(3);  //CO20190520 - no need to make copies
  //DX20190709 - calculate and store once = speed - START
  vector<xvector<double> > l1, l2, l3;
  vector<int> a_index, b_index, c_index;
  for(int i=i1;i<=i2;i++){l1.push_back(i*a1);a_index.push_back(i);}
  for(int j=j1;j<=j2;j++){l2.push_back(j*a2);b_index.push_back(j);}
  for(int k=k1;k<=k2;k++){l3.push_back(k*a3);c_index.push_back(k);}
  //DX20191218 - calculate and store once = speed - END

  // resize vectors - DX20191122
  uint num_grid_atoms = atoms.size()*l1.size()*l2.size()*l3.size();
  grid_atoms.resize(num_grid_atoms);
  grid_atoms_pc2scMap.resize(atoms.size()); //DX20191218 - should be the size of the the primitive cell, not grid
  grid_atoms_sc2pcMap.resize(num_grid_atoms);

  uint grid_atom_count = 0; // keep track of index - DX20191122

  for(uint iat=0;iat<atoms.size();iat++){
    //grid_atoms.push_back(atoms[iat]);  // put first the unit cell ! //DX20190709 - at to [] = speed increase
    //grid_atoms_pc2scMap.push_back(grid_atoms.size()-1); //CO20171025 
    //grid_atoms_sc2pcMap.push_back(iat); //CO20171025
    grid_atoms[grid_atom_count] = atoms[iat];  // put first the unit cell ! //DX20190709 - at to [] = speed increase
    grid_atoms_pc2scMap[grid_atom_count] = iat; //DX20191218 - use the index of the primitive cell not the running count of grid_atoms since it was resized
    grid_atoms_sc2pcMap[grid_atom_count] = iat; //CO20171025
    grid_atom_count++; //DX20191122
  }
  //for(int i=i1;i<=i2;i++) {
  //for(int j=j1;j<=j2;j++) {
  //for(int k=k1;k<=k2;k++) {
  //if(i!=0 || j!=0 || k!=0) {
  //for(uint iat=0;iat<atoms.size();iat++) {
  //atom=atoms[iat]; //DX20190709 - at to [] = speed increase
  //atom.isincell=FALSE; // these are OUT OF CELL
  //atom.cpos=((double)i)*a1+((double)j)*a2+((double)k)*a3+atoms[iat].cpos; //DX20190709 - at to [] = speed increase
  //atom.fpos[1]=i+atoms[iat].fpos[1]; //DX20190709 - at to [] = speed increase
  //atom.fpos[2]=j+atoms[iat].fpos[2]; //DX20190709 - at to [] = speed increase
  //atom.fpos[3]=k+atoms[iat].fpos[3]; //DX20190709 - at to [] = speed increase
  //grid_atoms.push_back(atom);
  //grid_atoms_sc2pcMap.push_back(iat); //CO20171025
  //if(LDEBUG) { //CO20190520
  //cerr << __AFLOW_FUNC__ << " grid_atoms[" << grid_atoms.size()-1 << "].cpos=" << grid_atoms.back().cpos << endl; //CO20190520
  //cerr << __AFLOW_FUNC__ << " grid_atoms[" << grid_atoms.size()-1 << "].fpos=" << grid_atoms.back().fpos << endl; //CO20190520
  //} //CO20190520
  //}
  //}
  //}
  //}
  //}
  xvector<double> a_component(3), ab_component(3), abc_component(3); //DX+ME20191107 - define outside loop (speed increase)
  uint natoms = atoms.size(); //DX20191107 - initialize natoms outside loop (speed increase)
  for(uint i=0;i<l1.size();i++) {
    a_component = l1[i];                           //DX : i*lattice(1)
    for(uint j=0;j<l2.size();j++) {
      ab_component = a_component + l2[j];          //DX : i*lattice(1) + j*lattice(2)
      for(uint k=0;k<l3.size();k++) {
        //DX20191218 [WRONG INDICES] if(i!=0 || j!=0 || k!=0)
        if(a_index[i]!=0 || b_index[j]!=0 || c_index[k]!=0) //DX20191218
        { //CO20200106 - patching for auto-indenting
          abc_component = ab_component + l3[k];    //DX : i*lattice(1) + j*lattice(2) + k*lattice(3)
          for(uint iat=0;iat<natoms;iat++) {       //DX20191107 - replace atoms.size() with natoms
            atom=atoms[iat];                   //DX20190709 - at to [] = speed increase
            atom.isincell=FALSE;                   // these are OUT OF CELL
            //DX20191127 [OBOSLETE] atom.cpos=abc_component+atoms[iat].cpos; //DX20190709 - at to [] = speed increase
            atom.cpos+=abc_component;              //DX20190709 - at to [] = speed increase //CO20191127 
            //DX20191127 [OBOSLETE] atom.fpos[1]=a_index[i]+atoms[iat].fpos[1]; //DX20190709 - at to [] = speed increase
            //DX20191127 [OBOSLETE] atom.fpos[2]=b_index[j]+atoms[iat].fpos[2]; //DX20190709 - at to [] = speed increase
            //DX20191127 [OBOSLETE] atom.fpos[3]=c_index[k]+atoms[iat].fpos[3]; //DX20190709 - at to [] = speed increase
            atom.fpos[1]+=a_index[i];              //DX20190709 - at to [] = speed increase //CO20191127
            atom.fpos[2]+=b_index[j];              //DX20190709 - at to [] = speed increase //CO20191127
            atom.fpos[3]+=c_index[k];              //DX20190709 - at to [] = speed increase //CO20191127
            //DX20191122 [OBSOLETE-PUSH_BACK] grid_atoms.push_back(atom);
            //DX20191122 [OBSOLETE-PUSH_BACK] grid_atoms_sc2pcMap.push_back(iat); //CO20171025
            grid_atoms[grid_atom_count] = atom;
            grid_atoms_sc2pcMap[grid_atom_count] = iat; //CO20171025
            //DX20200320 [OBSOLETE-moved outside of loop for speed] if(LDEBUG) { //CO20190520
            //DX20200320 [OBSOLETE-moved outside of loop for speed]   //DX20191122 [OBSOLETE-PUSH_BACK] cerr << __AFLOW_FUNC__ << " grid_atoms[" << grid_atoms.size()-1 << "].cpos=" << grid_atoms.back().cpos << endl; //CO20190520
            //DX20200320 [OBSOLETE-moved outside of loop for speed]   //DX20191122 [OBSOLETE-PUSH_BACK] cerr << __AFLOW_FUNC__ << " grid_atoms[" << grid_atoms.size()-1 << "].fpos=" << grid_atoms.back().fpos << endl; //CO20190520
            //DX20200320 [OBSOLETE-moved outside of loop for speed]   cerr << __AFLOW_FUNC__ << " grid_atoms[" << grid_atom_count << "].cpos=" << grid_atoms[grid_atom_count].cpos << endl; //CO20190520
            //DX20200320 [OBSOLETE-moved outside of loop for speed]   cerr << __AFLOW_FUNC__ << " grid_atoms[" << grid_atom_count << "].fpos=" << grid_atoms[grid_atom_count].fpos << endl; //CO20190520
            //DX20200320 [OBSOLETE-moved outside of loop for speed]   cerr << __AFLOW_FUNC__ << " grid_atoms[" << grid_atom_count << "]=" << grid_atoms[grid_atom_count] << endl; //DX20191218
            //DX20200320 [OBSOLETE-moved outside of loop for speed]   cerr << __AFLOW_FUNC__ << " grid_atoms_sc2pcMap[" << grid_atom_count << "]=" << grid_atoms_sc2pcMap[grid_atom_count] << endl; //DX20191218
            //DX20200320 [OBSOLETE-moved outside of loop for speed] } //CO20190520
            grid_atom_count++; //DX20191122
          }
        }
      }
    }
  }
  //DX20200320 - moved outside of loop so that LDEBUG boolean is not checked every time (speed increase when grid atoms is large)
  if(LDEBUG) { //CO20190520
    for(uint i=0;i<grid_atoms.size();i++){
      //DX20191122 [OBSOLETE-PUSH_BACK] cerr << __AFLOW_FUNC__ << " grid_atoms[" << grid_atoms.size()-1 << "].cpos=" << grid_atoms.back().cpos << endl; //CO20190520
      //DX20191122 [OBSOLETE-PUSH_BACK] cerr << __AFLOW_FUNC__ << " grid_atoms[" << grid_atoms.size()-1 << "].fpos=" << grid_atoms.back().fpos << endl; //CO20190520
      cerr << __AFLOW_FUNC__ << " grid_atoms[" << i << "].cpos=" << grid_atoms[i].cpos << endl; //CO20190520
      cerr << __AFLOW_FUNC__ << " grid_atoms[" << i << "].fpos=" << grid_atoms[i].fpos << endl; //CO20190520
      cerr << __AFLOW_FUNC__ << " grid_atoms[" << i << "]=" << grid_atoms[i] << endl; //DX20191218
      cerr << __AFLOW_FUNC__ << " grid_atoms_sc2pcMap[" << i << "]=" << grid_atoms_sc2pcMap[i] << endl; //DX20191218
    } //CO20190520
  }
  if(0){  //CO20190808 - quick check of mindist
    double min_dist_local=AUROSTD_MAX_DOUBLE,min_dist=AUROSTD_MAX_DOUBLE;
    for(uint i=0;i<grid_atoms.size()-1;i++){
      for(uint j=i+1;j<grid_atoms.size();j++){
        min_dist_local=aurostd::modulus(grid_atoms[i].cpos-grid_atoms[j].cpos);
        if(min_dist_local<min_dist){
          min_dist=min_dist_local;
        }
      }
    }
    if(!aurostd::isequal(min_dist,SYM::minimumDistance((*this)),0.1)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Minimum distance changed, check that atoms are not rotated",_RUNTIME_ERROR_);}
  }
  grid_atoms_calculated=TRUE;
  grid_atoms_dimsL[1]=i1;grid_atoms_dimsL[2]=j1;grid_atoms_dimsL[3]=k1;
  grid_atoms_dimsH[1]=i2;grid_atoms_dimsH[2]=j2;grid_atoms_dimsH[3]=k2;
  grid_atoms_number=grid_atoms.size();
  //  for(uint i=0;i<grid_atoms.size();i++)
  //   cerr << grid_atoms.at(i).cpos << endl;
  // cerr << grid_atoms.size() << endl;
  return grid_atoms.size();
}

int GenerateGridAtoms(xstructure& str) {return str.GenerateGridAtoms(-1,1,-1,1,-1,1);}
int GenerateGridAtoms(xstructure& str,double radius) {return str.GenerateGridAtoms(radius);}  //CO20200912 - double
int GenerateGridAtoms(xstructure& str,int d) {return str.GenerateGridAtoms(d);}
int GenerateGridAtoms(xstructure& str,int d1,int d2,int d3) {return str.GenerateGridAtoms(d1,d2,d3);}
int GenerateGridAtoms(xstructure& str,int i1,int i2,int j1,int j2,int k1,int k2) {return str.GenerateGridAtoms(i1,i2,j1,j2,k1,k2);} //DX20191218
int GenerateGridAtoms(xstructure& str,const xvector<int>& dims) {return str.GenerateGridAtoms(dims);} //CO20200912

// **************************************************************************
// GenerateLIJK table stuff
// **************************************************************************
namespace aurostd {   // INT
  class _ssort_int_value0123 {                    // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const {
        if(v1[0]<v2[0]) return TRUE;
        if(v1[0]==v2[0]) {
          if(v1[1]>v2[1]) return TRUE;
          if(v1[1]==v2[1]) {
            if(v1[2]>v2[2]) return TRUE;
            if(v1[2]==v2[2]) {
              if(v1[3]>v2[3]) return TRUE;
            }
          }
        }
        return FALSE;
      }
  };
}

int xstructure::GenerateLIJK(double radius) {
  lijk_table.clear();
  xvector<double> a1(3),a2(3),a3(3); //    a1,a2,a3 are the rows of the lattice matrix
  a1=lattice(1);a2=lattice(2);a3=lattice(3);    // a1,a2,a3 are the rows of the lattice matrix
  vector<int> _lijk(4);
  xvector<int> int_ijk(3);
  xvector<double> cpos(3),fpos(3);
  xvector<int> dims(3);
  dims=LatticeDimensionSphere(lattice,radius);  // for(int i=0;i<=3;i++) dims[i]+=1;
  lijk_dims[1]=dims[1];
  lijk_dims[2]=dims[2];
  lijk_dims[3]=dims[3];
  vector<vector<int> > _temp_lijk_table;
  _temp_lijk_table.clear();

  for(int i=-lijk_dims[1];i<=lijk_dims[1];i++) {
    for(int j=-lijk_dims[2];j<=lijk_dims[2];j++) {
      for(int k=-lijk_dims[3];k<=lijk_dims[3];k++) {
        cpos=((double)i)*a1+((double)j)*a2+((double)k)*a3;
        // 	if(aurostd::modulus(cpos)<=radius)  // I put them all since I will not be scanning through.
        {
          _lijk[0]=(int) (10000000.0*sqrt(aurostd::modulus(cpos))); // for the next sorting
          _lijk[1]=i;_lijk[2]=j;_lijk[3]=k;
          _temp_lijk_table.push_back(_lijk);
        }
      }
    }
  }
  sort(_temp_lijk_table.begin(),_temp_lijk_table.end(),aurostd::_ssort_int_value0123());

  for(uint i=0;i<_temp_lijk_table.size();i++) {
    _temp_lijk_table[i][0]=i; // relabel
    int_ijk[1]=_temp_lijk_table.at(i)[1];
    int_ijk[2]=_temp_lijk_table.at(i)[2];
    int_ijk[3]=_temp_lijk_table.at(i)[3];
    lijk_table.push_back(int_ijk);
  }

  for(uint i=0;i<lijk_table.size();i++) {
    cpos=((double)lijk_table.at(i)[1])*a1+((double)lijk_table.at(i)[2])*a2+((double)lijk_table.at(i)[3])*a3;
    fpos[1]=lijk_table.at(i)[1];fpos[2]=lijk_table.at(i)[2];fpos[3]=lijk_table.at(i)[3];
    lijk_cpos.push_back(cpos);
    lijk_fpos.push_back(fpos);
  }

  lijk_calculated=TRUE;
  return lijk_table.size();
}


// **************************************************************************
// LIJK table: L => IJK
// **************************************************************************
void l2ijk(const xstructure& str,const int &l,int &i,int &j,int &k) {
  if(l<0) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"l2ijk():","l2ijk error: l<0",_VALUE_RANGE_);}
  if(l>(int) str.lijk_table.size()) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"l2ijk():","l>str.lijk_table.size()",_VALUE_RANGE_);}
  i=str.lijk_table[l][1];j=str.lijk_table[l][2];k=str.lijk_table[l][3];
}

void l2ijk(const xstructure& str,const int &l,xvector<int>& ijk) {
  if(l<0) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"l2ijk():","l2ijk error: l<0",_VALUE_RANGE_);}
  if(l>(int) str.lijk_table.size()) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"l2ijk():","l>str.lijk_table.size()",_VALUE_RANGE_);}
  ijk[1]=str.lijk_table[l][1];ijk[2]=str.lijk_table[l][2];ijk[3]=str.lijk_table[l][3];
}

xvector<int> l2ijk(const xstructure& str,const int &l) {
  xvector<int> ijk(3);
  int i,j,k;
  l2ijk(str,l,i,j,k);
  ijk(1)=i;ijk(2)=j;ijk(3)=k;
  return ijk;
}


// **************************************************************************
// LIJK table: IJK => L
// **************************************************************************
void ijk2l(const xstructure& str,int &l,const int &i,const int &j,const int &k) {
  l=-1;
  if(i<-str.lijk_dims(1) || i>str.lijk_dims(1)) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"ijk2l():","i out of boundary",_VALUE_RANGE_);}
  if(j<-str.lijk_dims(2) || j>str.lijk_dims(2)) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"ijk2l():","j out of boundary",_VALUE_RANGE_);}
  if(k<-str.lijk_dims(3) || k>str.lijk_dims(3)) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"ijk2l():","k out of boundary",_VALUE_RANGE_);}
  for(uint ll=0;ll<str.lijk_table.size();ll++) // start search
    if(str.lijk_table.at(ll)[1]==i)         // faster comparison one at a time
      if(str.lijk_table.at(ll)[2]==j)       // faster comparison one at a time
        if(str.lijk_table.at(ll)[3]==k)     // faster comparison one at a time
          l=ll;
  if(l<0) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"ijk2l():","l not found",_VALUE_RANGE_);}
}

void ijk2l(const xstructure& str,int &l,const xvector<int>& ijk) {
  ijk2l(str,l,ijk(1),ijk(2),ijk(3));
}


int ijk2l(const xstructure& str,const int &i,const int &j,const int &k) {
  int l;
  ijk2l(str,l,i,j,k);
  return l;
}

int ijk2l(const xstructure& str,const xvector<int>& ijk) {
  int l;
  ijk2l(str,l,ijk(1),ijk(2),ijk(3));
  return l;
}

// **************************************************************************
// input2AIMSxstr input2ABINITxstr input2QExstr input2VASPxstr
// ***************************************************************************
xstructure input2AIMSxstr(istream& input) {
  xstructure a(input,IOAFLOW_AUTO);
  //  if(a.iomode==IOVASP_AUTO || a.iomode==IOVASP_POSCAR || a.iomode==IOVASP_ABCCAR || a.iomode==IOVASP_WYCKCAR) 
  a.xstructure2aims();
  return a;
}

xstructure input2ABINITxstr(istream& input) {
  xstructure a(input,IOAFLOW_AUTO);
  //  if(a.iomode==IOVASP_AUTO || a.iomode==IOVASP_POSCAR || a.iomode==IOVASP_ABCCAR || a.iomode==IOVASP_WYCKCAR) 
  a.xstructure2abinit();
  return a;
}

xstructure input2QExstr(istream& input) {
  xstructure a(input,IOAFLOW_AUTO);
  //  if(a.iomode==IOVASP_AUTO || a.iomode==IOVASP_POSCAR || a.iomode==IOVASP_ABCCAR || a.iomode==IOVASP_WYCKCAR) 
  a.xstructure2qe();
  return a;
}

xstructure input2VASPxstr(istream& input,bool vasp5) {  //CO20210119 - added vasp5
  xstructure a(input,IOAFLOW_AUTO);
  //  if(a.iomode==IOQE_AUTO || a.iomode==IOQE_GEOM)
  a.xstructure2vasp();
  if(vasp5){a.is_vasp4_poscar_format=false;a.is_vasp5_poscar_format=true;}  //CO20210119
  //  cerr << a.title << endl;
  return a;
}

xstructure input2ITCxstr(istream& input) {  //CO20220613
  xstructure a(input,IOAFLOW_AUTO);
  a.xstructure2itc();
  return a;
}

xstructure input2ELKxstr(istream& input) { //DX20200313
  xstructure a(input,IOAFLOW_AUTO);
  a.xstructure2elk();
  return a;
}

xstructure input2ATATxstr(istream& input) { //SD20220123
  xstructure a(input,IOAFLOW_AUTO);
  a.xstructure2atat();
  return a;
}

// **************************************************************************
// r_lattice
// ***************************************************************************
xvector<double> r_lattice(const xstructure& str,const int &l) {
  xvector<double> rrr(3);
  int i,j,k;
  l2ijk(str,l,i,j,k);
  rrr=((double)i)*str.lattice(1)+((double)j)*str.lattice(2)+((double)k)*str.lattice(3);
  return rrr;
}

xvector<double> r_lattice(const xstructure& str,const int &i,const int &j,const int &k) {
  xvector<double> rrr(3);
  rrr=((double)i)*str.lattice(1)+((double)j)*str.lattice(2)+((double)k)*str.lattice(3);
  return rrr;
}

xvector<double> r_lattice(const xstructure& str,const xvector<int>& ijk) {
  xvector<double> rrr(3);
  rrr=((double)ijk(1))*str.lattice(1)+((double)ijk(2))*str.lattice(2)+((double)ijk(3))*str.lattice(3);
  return rrr;
}

// **************************************************************************
// Function GetUnitCellRep
// **************************************************************************
// Dane Morgan / Stefano Curtarolo
// This function calculated the representation of a point in terms of a position in cell 0 and its true unit cell.
// Output position in cell 0 is given with same coordinate type (Cart or Direct) as input position.
// Note that there is always an ambiguity of being in a given cell with position 0, or one cell over with position 1.  This
// is broken by forcing all cell positions to be at 0 if they are within TOL of 1.  This gives consistent cell image locations.
//CO20170717 - Compare this with PBC() function, this function brings in to 0th cell, PBC brings in to cell between -0.5 and 0.5
// PBC() is good for minimizing overall fpos
void GetUnitCellRep(const xvector<double>& ppos,xvector<double>& p_cell0,
    xvector<int>& ijk,const xmatrix<double>& lattice,
    const bool coord_flag) {
  double TOL=1e-11;
  xvector<double> cpos(3),fpos(3);
  // p_cell0=xvector<double> (3);
  // ijk=xvector<int> (3);

  clear(ijk);clear(p_cell0);

  if(coord_flag==TRUE) // ppos is in cartesian coords
    fpos=C2F(lattice,ppos);

  for(int ic=1;ic<=3;ic++) {
    int s;
    s=SignNoZero(fpos(ic));
    if(s>0) ijk(ic)=int(fpos(ic));
    else    ijk(ic)=int(fpos(ic))-1;
    p_cell0(ic)=fpos(ic)-ijk(ic);
    // Boundary Correction.  p_cell0 is >=0 and <1.  
    // If p_cell0 is within TOL of 1
    // then set it to 0 and shift ijk up by 1.
    if(abs(1-p_cell0(ic))<TOL) {
      p_cell0(ic)=0;
      ijk(ic)=ijk(ic)+1;
    }
  }
  if(coord_flag==TRUE) // ppos is in cartesian coords
    p_cell0=F2C(lattice,p_cell0);
}

// **************************************************************************
// Function COMPARE_GetNeighData
// **************************************************************************
// For sort algorithm in GetNeighData
class compare_GetNeighData {
  public:  
    int operator()(const _atom& a, const _atom& b) {
      bool LDEBUG=(FALSE || XHOST.DEBUG);
      if(LDEBUG) cerr << "operator()(const _atom& a, const _atom& b)" << endl;
      double tol=1e-15;
      // Sort by distance
      if(aurostd::isequal(GetDistFromOrigin(a),GetDistFromOrigin(b),tol)) {
        if(LDEBUG) cerr << "operator()(const _atom& a, const _atom& b)  isequal(GetDistFromOrigin(a),GetDistFromOrigin(b),tol)" << endl;
        // Sort by unit cell values
        if(a.name==b.name) {
          if(LDEBUG) cerr << "operator()(const _atom& a, const _atom& b)  a.name==b.name" << endl;
          xvector<int> ijka=a.ijk;
          xvector<int> ijkb=b.ijk;
          int va=100*ijka(1)+10*ijka(2)+1*ijka(3);
          int vb=100*ijkb(1)+10*ijkb(2)+1*ijkb(3);
          return va<vb;
        } else {
          if(LDEBUG) cerr << "operator()(const _atom& a, const _atom& b)  a.name!=b.name" << endl;
          return a.name<b.name;
        }
      } else {
        if(LDEBUG) cerr << "operator()(const _atom& a, const _atom& b)  !isequal(GetDistFromOrigin(a),GetDistFromOrigin(b),tol)" << endl;
        return GetDistFromOrigin(a)<GetDistFromOrigin(b);
      }
      // Sort by name
      // if(a.name==b.name) {
      //  return GetDistFromOrigin(a)<GetDistFromOrigin(b);
      // }
      // else {
      //  return a.name<b.name;
      // }
    }
};

//RF20200831 - checkStructure - START
// **************************************************************************
// checkStructure
// **************************************************************************
// rescale structure to 1 and check whether e.g. species and atoms are present
void xstructure::checkStructure(){
  bool LDEBUG = (FALSE || XHOST.DEBUG);
  stringstream message;
  (*this).ReScale(1.0); // rescales scaling factor in second line of POSCAR to 1, needed for correct distances
  //throw some general information such as input structure
  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << endl << "INPUT STRUCTURE:" << endl;
    cerr << __AFLOW_FUNC__ << (*this) << endl;
  }
  // check whether the species vector is populated, otherwise throw error
  if ((*this).species.size() == 0){
    message << " BAD NEWS: It seems there are no species in the structure. Please adjust the structure and rerun.";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
  }
  // check whether there are any atoms in the structure
  if ((*this).atoms.size() == 0){
    message << " BAD NEWS: It seems there are no atoms in the structure. Please adjust the structure and rerun.";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
  }
  // if species of atoms are not known like in VASP4 format, throw error
  for(uint k=0,ksize=(*this).atoms.size();k<ksize;k++){
    if ((*this).atoms[k].cleanname == ""){
      message << " BAD NEWS: It seems you are providing a structure without complete species information as input. This implementation requires a structure with the species information included. For a VASP4 POSCAR, the species must be written on the right side next to the coordinates for each atom. Please adjust the structure and rerun.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
    }
  }
}
//RF20200831 - checkStructure - END


// **************************************************************************
// Function GetNeighbors
// **************************************************************************
// rewrite of GetNeighData()
// atoms_cell is the atoms for which neighbors are found (could be iatoms), hence it determines the sizes of i_neighbors and distances
// i_neighbors are indices to grid_atoms (xstructure attribute)
// CO20200912
void xstructure::GetNeighbors(deque<deque<uint> >& i_neighbors,deque<deque<double> >& distances,double rmin,bool prim,bool unique_only){
  deque<_atom> atoms_cell;
  return GetNeighbors(atoms_cell,i_neighbors,distances,rmin,prim,unique_only);
}
void xstructure::GetNeighbors(deque<_atom>& atoms_cell,deque<deque<uint> >& i_neighbors,deque<deque<double> >& distances,double rmin,bool prim,bool unique_only){
  double rmax=RadiusSphereLattice((*this).scale*(*this).lattice);
  return GetNeighbors(atoms_cell,i_neighbors,distances,rmax,rmin,prim,unique_only);
}
void xstructure::GetNeighbors(deque<deque<uint> >& i_neighbors,deque<deque<double> >& distances,double rmax,double rmin,bool prim,bool unique_only){
  deque<_atom> atoms_cell;
  return GetNeighbors(atoms_cell,i_neighbors,distances,rmax,rmin,prim,unique_only);
}
void xstructure::GetNeighbors(deque<_atom>& atoms_cell,deque<deque<uint> >& i_neighbors,deque<deque<double> >& distances,double rmax,double rmin,bool prim,bool unique_only){
  bool LDEBUG=(FALSE || XHOST.DEBUG);

  uint i=0,k=0;

  if(prim){GetPrimitive();}
  ReScale(1.0);
  //get atoms_cell
  atoms_cell.clear(); //clear
  for(i=0;i<i_neighbors.size();i++){i_neighbors[i].clear();} i_neighbors.clear(); //clear
  for(i=0;i<distances.size();i++){distances[i].clear();} distances.clear(); //clear
  vector<uint> atomscell2atoms_mapping;
  if(unique_only){
    if(iatoms_calculated==false){CalculateSymmetry();}
    for(i=0;i<iatoms.size();i++){
      atoms_cell.push_back(atoms[iatoms[i][0]]);
      atomscell2atoms_mapping.push_back(iatoms[i][0]);
      i_neighbors.push_back(deque<uint>(0));
      distances.push_back(deque<double>(0));
    }
  }else{
    for(uint i=0;i<atoms.size();i++){
      atoms_cell.push_back(atoms[i]);
      atomscell2atoms_mapping.push_back(i);
      i_neighbors.push_back(deque<uint>(0));
      distances.push_back(deque<double>(0));
    }
  }

  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " atoms_cell.size()=" << atoms_cell.size() << endl;
    cerr << __AFLOW_FUNC__ << " atomscell2atoms_mapping=" << aurostd::joinWDelimiter(atomscell2atoms_mapping,",") << endl;
  }

  if(LDEBUG){cerr << __AFLOW_FUNC__ << " GenerateGridAtoms(): starting" << endl;}
  GenerateGridAtoms(rmax);
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " GenerateGridAtoms(): done" << endl;}
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " grid_atoms_number=" << grid_atoms_number << endl;}

  if(LDEBUG){cerr << __AFLOW_FUNC__ << " generating neighbors" << endl;}
  //[STILL SLOWER THAN AUROSTD::SORT()]//[SHAVES 1 SECOND]vector<double> v_mindist;
  double dist=0.0;
  uint atom1=0,atom2=0;
  for(i=0;i<atoms_cell.size();i++){
    //[STILL SLOWER THAN AUROSTD::SORT()]//[SHAVES 1 SECOND]v_mindist.push_back(AUROSTD_MAX_DOUBLE);
    atom1=grid_atoms_pc2scMap[atomscell2atoms_mapping[i]];
    for(atom2=0;atom2<(uint)grid_atoms_number;atom2++){
      if(atom1==atom2){continue;} //skip self
      dist=AtomDist(grid_atoms[atom1],grid_atoms[atom2]); //distance
      if(dist>rmin && dist<=rmax){
        i_neighbors[i].push_back(atom2);
        distances[i].push_back(dist);
        //[STILL SLOWER THAN AUROSTD::SORT()]//[SHAVES 1 SECOND]if(dist<v_mindist.back()){v_mindist.back()=dist;}
      }
    }
  }

  if(0){
    //an attempt to shave off more runtime: only shaves off a second AT BEST for "make check_cce"
    //the slowest part of the algorithm WAS sorting (using insertion sort), so reducing vector size would help
    //now using quicksort, reducing the vector size does not help much (very efficient)
    //to fully realize the 1 second shave off, replace min(distances[k]) with v_mindist (see above)
    vector<uint> i2remove;
    double cutoff=0.0;
    for(k=0;k<distances.size();k++){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " START distances[k=" << k << "].size()=" << distances[k].size() << endl;}
      i2remove.clear();
      cutoff=1.25*min(distances[k]); //v_mindist[k];  //focus on first shell: does not require an input tol
      for(i=0;i<distances[k].size();i++){
        if(distances[k][i]>cutoff){i2remove.push_back(i);}
      }
      for(i=i2remove.size()-1;i<i2remove.size();i--){
        i_neighbors[k].erase(i_neighbors[k].begin()+i2remove[i]);
        distances[k].erase(distances[k].begin()+i2remove[i]);
      }
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " STOP  distances[k=" << k << "].size()=" << distances[k].size() << endl;}
    }
  }

  //now sort
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " sorting" << endl;}
  //[TOO SLOW FOR LARGE VECTORS]uint uint_tmp=0;
  //[TOO SLOW FOR LARGE VECTORS]double double_tmp=0;
  //[TOO SLOW FOR LARGE VECTORS]for(k=0;k<i_neighbors.size();k++){
  //[TOO SLOW FOR LARGE VECTORS]  for(i=0;i<i_neighbors[k].size();i++){
  //[TOO SLOW FOR LARGE VECTORS]    for(j=i+1;j<i_neighbors[k].size();j++){
  //[TOO SLOW FOR LARGE VECTORS]      if(distances[k][j]<distances[k][i]){
  //[TOO SLOW FOR LARGE VECTORS]        //std::swap(i_neighbors[k][i],i_neighbors[k][j]);
  //[TOO SLOW FOR LARGE VECTORS]        uint_tmp=i_neighbors[k][i];
  //[TOO SLOW FOR LARGE VECTORS]        i_neighbors[k][i]=i_neighbors[k][j];
  //[TOO SLOW FOR LARGE VECTORS]        i_neighbors[k][j]=uint_tmp;
  //[TOO SLOW FOR LARGE VECTORS]        //std::swap(distances[k][i],distances[k][j]);
  //[TOO SLOW FOR LARGE VECTORS]        double_tmp=distances[k][i];
  //[TOO SLOW FOR LARGE VECTORS]        distances[k][i]=distances[k][j];
  //[TOO SLOW FOR LARGE VECTORS]        distances[k][j]=double_tmp;
  //[TOO SLOW FOR LARGE VECTORS]      }
  //[TOO SLOW FOR LARGE VECTORS]    }
  //[TOO SLOW FOR LARGE VECTORS]  }
  //[TOO SLOW FOR LARGE VECTORS]}
  //[STILL SLOWER THAN AUROSTD::SORT()]vector<uint> i_sorting;
  //[STILL SLOWER THAN AUROSTD::SORT()]deque<uint> i_neighborsk; //copy
  //[STILL SLOWER THAN AUROSTD::SORT()]deque<double> distancesk; //copy
  //[STILL SLOWER THAN AUROSTD::SORT()]for(k=0;k<i_neighbors.size();k++){
  //[STILL SLOWER THAN AUROSTD::SORT()]  i_neighborsk=i_neighbors[k];  //copy
  //[STILL SLOWER THAN AUROSTD::SORT()]  distancesk=distances[k];  //copy
  //[STILL SLOWER THAN AUROSTD::SORT()]  i_sorting.clear();for(i=0;i<i_neighbors[k].size();i++){i_sorting.push_back(i);}
  //[STILL SLOWER THAN AUROSTD::SORT()]  std::sort(i_sorting.begin(),i_sorting.end(),[&](uint a, uint b) {return distances[k][a]<distances[k][b];});
  //[STILL SLOWER THAN AUROSTD::SORT()]  //reorder
  //[STILL SLOWER THAN AUROSTD::SORT()]  //https://stackoverflow.com/questions/838384/reorder-vector-using-a-vector-of-indices
  //[STILL SLOWER THAN AUROSTD::SORT()]  //do NOT use a reordering algorithm, it takes a LONG time
  //[STILL SLOWER THAN AUROSTD::SORT()]  //instead make a copy, memory is not a problem here, as the order of the memory used does not change (2N vs. 4N)
  //[STILL SLOWER THAN AUROSTD::SORT()]  for(i=0;i<i_sorting.size();i++){
  //[STILL SLOWER THAN AUROSTD::SORT()]    i_neighbors[k][i]=i_neighborsk[i_sorting[i]];
  //[STILL SLOWER THAN AUROSTD::SORT()]    distances[k][i]=distancesk[i_sorting[i]];
  //[STILL SLOWER THAN AUROSTD::SORT()]  }
  //[STILL SLOWER THAN AUROSTD::SORT()]}
  for(k=0;k<i_neighbors.size();k++){aurostd::sort(distances[k],i_neighbors[k]);}

  if(LDEBUG){
    for(k=0;k<i_neighbors.size();k++){
      cerr << __AFLOW_FUNC__ << " ATOMS_CELL[k=" << k << "]: " << atoms_cell[k].name << endl;
      for(i=0;i<i_neighbors[k].size();i++){
        cerr << "  neighbor[i=" << i << "]: " << grid_atoms[i_neighbors[k][i]].name << " dist=" << distances[k][i] << endl;
      }
    }
  }

}


void GetNeighbors(const xstructure& xstr_in,deque<deque<uint> >& i_neighbors,deque<deque<double> >& distances,double rmin,bool prim,bool unique_only){
  xstructure xstr(xstr_in);
  return xstr.GetNeighbors(i_neighbors,distances,rmin,prim,unique_only);
}
void GetNeighbors(const xstructure& xstr_in,deque<_atom>& atoms_cell,deque<deque<uint> >& i_neighbors,deque<deque<double> >& distances,double rmin,bool prim,bool unique_only){
  xstructure xstr(xstr_in);
  return xstr.GetNeighbors(atoms_cell,i_neighbors,distances,rmin,prim,unique_only);
}
void GetNeighbors(const xstructure& xstr_in,deque<deque<uint> >& i_neighbors,deque<deque<double> >& distances,double rmax,double rmin,bool prim,bool unique_only){
  xstructure xstr(xstr_in);
  return xstr.GetNeighbors(i_neighbors,distances,rmax,rmin,prim,unique_only);
}
void GetNeighbors(const xstructure& xstr_in,deque<_atom>& atoms_cell,deque<deque<uint> >& i_neighbors,deque<deque<double> >& distances,double rmax,double rmin,bool prim,bool unique_only){
  xstructure xstr(xstr_in);
  return xstr.GetNeighbors(atoms_cell,i_neighbors,distances,rmax,rmin,prim,unique_only);
}


// **************************************************************************
// Function GetCoordinations
// **************************************************************************
// CO20200912
void xstructure::GetCoordinations(deque<deque<uint> >& coordinations,double rmin,double tol,bool prim,bool unique_only){
  deque<deque<uint> > i_neighbors;
  deque<deque<double> > distances;
  return GetCoordinations(i_neighbors,distances,coordinations,tol,rmin,prim,unique_only);
}
void xstructure::GetCoordinations(deque<_atom>& atoms_cell,deque<deque<uint> >& coordinations,double rmin,double tol,bool prim,bool unique_only){
  deque<deque<uint> > i_neighbors;
  deque<deque<double> > distances;
  return GetCoordinations(atoms_cell,i_neighbors,distances,coordinations,tol,rmin,prim,unique_only);
}
void xstructure::GetCoordinations(deque<deque<uint> >& coordinations,double rmax,double rmin,double tol,bool prim,bool unique_only){
  deque<deque<uint> > i_neighbors;
  deque<deque<double> > distances;
  return GetCoordinations(i_neighbors,distances,coordinations,rmax,rmin,tol,prim,unique_only);
}
void xstructure::GetCoordinations(deque<_atom>& atoms_cell,deque<deque<uint> >& coordinations,double rmax,double rmin,double tol,bool prim,bool unique_only){
  deque<deque<uint> > i_neighbors;
  deque<deque<double> > distances;
  return GetCoordinations(atoms_cell,i_neighbors,distances,coordinations,rmax,rmin,tol,prim,unique_only);
}
void xstructure::GetCoordinations(deque<deque<uint> >& i_neighbors,deque<deque<double> >& distances,deque<deque<uint> >& coordinations,double rmin,double tol,bool prim,bool unique_only){
  deque<_atom> atoms_cell;
  return GetCoordinations(atoms_cell,i_neighbors,distances,coordinations,tol,rmin,prim,unique_only);
}
void xstructure::GetCoordinations(deque<_atom>& atoms_cell,deque<deque<uint> >& i_neighbors,deque<deque<double> >& distances,deque<deque<uint> >& coordinations,double rmin,double tol,bool prim,bool unique_only){
  double rmax=RadiusSphereLattice((*this).scale*(*this).lattice);
  return GetCoordinations(atoms_cell,i_neighbors,distances,coordinations,rmax,rmin,tol,prim,unique_only);
}
void xstructure::GetCoordinations(deque<deque<uint> >& i_neighbors,deque<deque<double> >& distances,deque<deque<uint> >& coordinations,double rmax,double rmin,double tol,bool prim,bool unique_only){
  deque<_atom> atoms_cell;
  return GetCoordinations(atoms_cell,i_neighbors,distances,coordinations,rmax,rmin,tol,prim,unique_only);
}
void xstructure::GetCoordinations(deque<_atom>& atoms_cell,deque<deque<uint> >& i_neighbors,deque<deque<double> >& distances,deque<deque<uint> >& coordinations,double rmax,double rmin,double tol,bool prim,bool unique_only){
  bool LDEBUG=(FALSE || XHOST.DEBUG);

  GetNeighbors(atoms_cell,i_neighbors,distances,rmax,rmin,prim,unique_only);

  uint i=0,j=0,k=0;
  for(i=0;i<coordinations.size();i++){coordinations[i].clear();} coordinations.clear(); //clear

  for(k=0;k<i_neighbors.size();k++){
    coordinations.push_back(deque<uint>(0));
    j=0;  //index to start at
    coordinations.back().push_back(0);
    for(i=0;i<i_neighbors[k].size();i++){
      if(aurostd::isequal(distances[k][j],distances[k][i],tol)){coordinations.back().back()+=1;}
      else{
        j=i;
        coordinations.back().push_back(1);
      }
    }
  }

  if(LDEBUG){
    for(k=0;k<i_neighbors.size();k++){
      cerr << __AFLOW_FUNC__ << " ATOMS_CELL[k=" << k << "]: " << atoms_cell[k].name << endl;
      for(i=0;i<coordinations[k].size();i++){
        cerr << "  coordination[i=" << i << "]: " << coordinations[k][i] << endl;
      }
    }
  }

}

void GetCoordinations(const xstructure& xstr_in,deque<deque<uint> >& coordinations,double rmin,double tol,bool prim,bool unique_only){
  xstructure xstr(xstr_in);
  return xstr.GetCoordinations(coordinations,tol,rmin,prim,unique_only);
}
void GetCoordinations(const xstructure& xstr_in,deque<_atom>& atoms_cell,deque<deque<uint> >& coordinations,double rmin,double tol,bool prim,bool unique_only){
  xstructure xstr(xstr_in);
  return xstr.GetCoordinations(atoms_cell,coordinations,tol,rmin,prim,unique_only);
}
void GetCoordinations(const xstructure& xstr_in,deque<deque<uint> >& coordinations,double rmax,double rmin,double tol,bool prim,bool unique_only){
  xstructure xstr(xstr_in);
  return xstr.GetCoordinations(coordinations,rmax,rmin,tol,prim,unique_only);
}
void GetCoordinations(const xstructure& xstr_in,deque<_atom>& atoms_cell,deque<deque<uint> >& coordinations,double rmax,double rmin,double tol,bool prim,bool unique_only){
  xstructure xstr(xstr_in);
  return xstr.GetCoordinations(atoms_cell,coordinations,rmax,rmin,tol,prim,unique_only);
}
void GetCoordinations(const xstructure& xstr_in,deque<deque<uint> >& i_neighbors,deque<deque<double> >& distances,deque<deque<uint> >& coordinations,double rmin,double tol,bool prim,bool unique_only){
  xstructure xstr(xstr_in);
  return xstr.GetCoordinations(i_neighbors,distances,coordinations,tol,rmin,prim,unique_only);
}
void GetCoordinations(const xstructure& xstr_in,deque<_atom>& atoms_cell,deque<deque<uint> >& i_neighbors,deque<deque<double> >& distances,deque<deque<uint> >& coordinations,double rmin,double tol,bool prim,bool unique_only){
  xstructure xstr(xstr_in);
  return xstr.GetCoordinations(atoms_cell,i_neighbors,distances,coordinations,rmin,tol,prim,unique_only);
}
void GetCoordinations(const xstructure& xstr_in,deque<deque<uint> >& i_neighbors,deque<deque<double> >& distances,deque<deque<uint> >& coordinations,double rmax,double rmin,double tol,bool prim,bool unique_only){
  xstructure xstr(xstr_in);
  return xstr.GetCoordinations(i_neighbors,distances,coordinations,rmax,rmin,tol,prim,unique_only);
}
void GetCoordinations(const xstructure& xstr_in,deque<_atom>& atoms_cell,deque<deque<uint> >& i_neighbors,deque<deque<double> >& distances,deque<deque<uint> >& coordinations,double rmax,double rmin,double tol,bool prim,bool unique_only){
  xstructure xstr(xstr_in);
  return xstr.GetCoordinations(atoms_cell,i_neighbors,distances,coordinations,rmax,rmin,tol,prim,unique_only);
}

// **************************************************************************
// Function GetNeighData
// **************************************************************************
// This function collects all the neighbor data between rmin and rmax and stores 
// it for each atom in a vector of atom objects in order of increasing distance.  

void xstructure::GetNeighData(const double rmax,deque<deque<_atom> >& neigh_mat,const double rmin) const {  //CO20220623
  deque<_atom> atoms_cell;
  return GetNeighData(atoms_cell,rmax,neigh_mat,rmin);
}
void xstructure::GetNeighData(deque<_atom>& atoms_cell,const double rmax,deque<deque<_atom> >& neigh_mat,const double rmin) const {
  //CO20220623 - rewritten using GetNeighbors() (MUCH faster)
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  xstructure xstr_tmp(*this); //make copy so we can modify xstr
  deque<deque<uint> > i_neighbors;
  deque<deque<double> > distances;
  bool prim=true;bool unique_only=false;
  xstr_tmp.GetNeighbors(atoms_cell,i_neighbors,distances,rmax,rmin,prim,unique_only);
  //now build neigh_mat
  uint i=0,j=0;
  for(i=0;i<neigh_mat.size();i++){neigh_mat[i].clear();}neigh_mat.clear();  //clear
  for(i=0;i<i_neighbors.size();i++){
    neigh_mat.push_back(deque<_atom>(0));
    neigh_mat[i].push_back(atoms_cell[i]);  //first atom is self
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " neigh_mat[i=" << i << "][j=0].cpos=" << neigh_mat[i].back().cpos << endl;}
    for(j=0;j<i_neighbors[i].size();j++){
      neigh_mat[i].push_back(xstr_tmp.grid_atoms[i_neighbors[i][j]]);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " neigh_mat[i=" << i << "][j=" << j << "].cpos=" << neigh_mat[i].back().cpos << "; dist=" << distances[i][j] << endl;}
    }
  }
}
void xstructure::GetNeighData_20220101(const deque<_atom>& in_atom_vec,		
    const double rmin, const double rmax,
    deque<deque<_atom> >& neigh_mat) const {  
  //CO20220623 - AVOID USING THIS FUNCTION, use the one above
  //this one relies on an input of in_atom_vec which is usually the atoms of the structure
  //there is only one exception: GetGoodShellPoints() in aflow_pflow_print.cpp; this function should be rewritten in the future
  //considering the name of this function, it should only operate on the atoms of the structure
  //no point of having two sources of truth, it's confusing
  //the rewrite for an input of atoms should avoid ad-hoc, a posteriori mappings introducing noise in positions: GetUnitCellRep() and ConvertAtomToLat()
  double epsilon=1.0e-3;  // sometimes you have wrong images due to roundoff
  deque<_atom> neigh_vec;
  // Get data from str.
  xstructure sstr(*this);
  // Set scale to 1 so you don't need to rescale coordinates.
  // [OBSOLETE] sstr=ReScale(sstr,1.0);
  sstr.ReScale(1.0);
  xmatrix<double> original_lattice(3,3);
  original_lattice=sstr.lattice;
  // Use Niggli cell to avoid missing neighbors in angular cells.
  // This involves converting your atoms into the 0th Niggli cell,
  // keeping track of the unit cell you came from, and then shifting back
  // at the end.  This shifting must be done for both atom_vec (at the
  // beginning) and neigh_mat (at the end).
  //sstr.NiggliUnitCellForm();  
  sstr.MinkowskiBasisReduction();
  //    cerr << sstr << endl;

  xmatrix<double> lat(3,3);
  lat=sstr.lattice;

  //  Convert input atom_vec to Niggli lattice
  deque<_atom> atom_vec=in_atom_vec;
  xmatrix<int>  atom_shifts((int) atom_vec.size()-1,3,0,1);
  for(uint ia=0;ia<atom_vec.size();ia++) {
    _atom a=atom_vec.at(ia);
    xvector<double> p_cell0(3);
    xvector<int>    atom_shifts_temp(3);
    GetUnitCellRep(a.cpos,p_cell0,atom_shifts_temp,lat,1);
    atom_vec.at(ia).cpos=p_cell0;
    atom_vec.at(ia)=ConvertAtomToLat(atom_vec.at(ia),lat);
    atom_shifts(ia,1)=atom_shifts_temp(1);
    atom_shifts(ia,2)=atom_shifts_temp(2);
    atom_shifts(ia,3)=atom_shifts_temp(3);
  }    

  // Create vector of all atoms in all periodic images needed so that for atoms in the unit cell every possible neighbor within rmax
  // is included.  This can be done by looping over all atoms in all unit cells such that if the unit cells are given by
  // n0*v0,n1*v1,n2*v2, then ni_max*vi>rmax for all i=0->2.  This is not rigorous but seems to work quite well.  It failed for very
  // angular cells, but by using the Niggli reduced cell this pitfall seems to be avoided (but no proof, so be careful).

  int imax,jmax,kmax;
  // Find imax
  // (algorithm 1) approximate
  imax=(int)max(rmax/aurostd::modulus(lat(1)),rmax/aurostd::modulus(lat(2)),rmax/aurostd::modulus(lat(3)))+2;
  jmax=imax;
  kmax=imax;
  // (algorithm 2) exact, the requirement is that atoms must in incell.
  // we should find nlat1,nlat2,nlat3
  // where nlat1 is the lat1 component that is perpendicular to the plane made by lat2 & lat3.
  // then imax=ceil(rmax/nlat1)
  // This algorithm is implemented in function LatticeDimensionSphere
  xvector<int> dims(3);
  dims=LatticeDimensionSphere(lat,rmax);
  imax=dims(1);
  jmax=dims(2);
  kmax=dims(3);
  xvector<int> ijk(3);
  deque<_atom> all_atom_vec;
  // latt maybe a rotated version of POSCAR, so need to incell-ized the fpos and cpos
  // sstr.BringInCell(); // with roundoff
  sstr.BringInCell(); //DX+RF20200121 - negative input is no longer accepted // no roundoff
  //DX+RF20200121 - negative input is no longer accepted : sstr.BringInCell(-1.0); // no roundoff
  for(ijk(1)=-imax;ijk(1)<=imax;ijk(1)++) {
    for(ijk(2)=-jmax;ijk(2)<=jmax;ijk(2)++) {
      for(ijk(3)=-kmax;ijk(3)<=kmax;ijk(3)++) {
        xvector<double> ctpos(3);
        xvector<double> ftpos(3);
        for(uint iat=0;iat<sstr.atoms.size();iat++) {
          _atom a;
          a=sstr.atoms.at(iat);
          a.name=sstr.atoms.at(iat).name;
          a.basis=iat; //[CO20200130 - number->basis]a.number=iat;
          a.ijk=ijk;
          for(int ic=1;ic<=3;ic++)
            ctpos(ic)=sstr.atoms.at(iat).cpos(ic)+ijk(1)*lat(1,ic)+ijk(2)*lat(2,ic)+ijk(3)*lat(3,ic);
          ftpos=C2F(lat,ctpos);
          a.cpos=ctpos;
          a.fpos=ftpos;
          a.type=sstr.atoms.at(iat).type;
          all_atom_vec.push_back(a);
        } // iat
      } // ijk
    } // ijk
  } // ijk

  // Now build neighbors list for each atom on atom list.  Each
  // neighbor list will be row in the neigh_mat matrix.
  for(uint ia1=0;ia1<atom_vec.size();ia1++) {
    xvector<double> corigin(3);
    _atom at=atom_vec.at(ia1);
    corigin = at.cpos;
    at.corigin=corigin;
    deque<_atom> neigh_vec;
    neigh_vec.push_back(at);
    for(uint ia2=0;ia2<all_atom_vec.size();ia2++) {
      double dist = AtomDist(at,all_atom_vec.at(ia2));
      if(dist<=rmax && dist>=rmin && dist>=epsilon) {
        all_atom_vec.at(ia2).corigin=corigin;
        neigh_vec.push_back(all_atom_vec.at(ia2));
      } // if
    } // ia2
    sort(neigh_vec.begin()+1,neigh_vec.end(),compare_GetNeighData());
    neigh_mat.push_back(neigh_vec);
  } // ia1

  //  Convert neigh_mat from Niggli cell to original lattice
  for(uint ia=0;ia<neigh_mat.size();ia++) {
    for(uint ja=0;ja<neigh_mat.at(ia).size();ja++) {
      xvector<double> fpos(3);
      xvector<double> cpos(3);
      _atom a;
      a=neigh_mat.at(ia).at(ja);
      fpos=a.fpos;
      for(int ic=1;ic<=3;ic++) {
        fpos(ic)=fpos(ic)+atom_shifts(ia,ic);
      }
      cpos=F2C(lat,fpos);
      a.cpos=cpos;
      neigh_mat.at(ia).at(ja)=ConvertAtomToLat(a,original_lattice);
    }
  }
}

//[CO20220623 - OBSOLETE]// **************************************************************************
//[CO20220623 - OBSOLETE]// Function GetStrNeighData
//[CO20220623 - OBSOLETE]// **************************************************************************
//[CO20220623 - OBSOLETE]// This function collects all the neighbor data out to some
//[CO20220623 - OBSOLETE]// cutoff and stores it for each atom in the structure.
//[CO20220623 - OBSOLETE]void xstructure::GetStrNeighData(const double cutoff,deque<deque<_atom> >& neigh_mat) const { //RF+CO20200513
//[CO20220623 - OBSOLETE]  deque<_atom> atom_vec;
//[CO20220623 - OBSOLETE]  neigh_mat.clear();
//[CO20220623 - OBSOLETE]  // Get data from str.
//[CO20220623 - OBSOLETE]  // Set scale to 1 so you don't need to rescale coordinates.
//[CO20220623 - OBSOLETE]  xstructure sstr(*this);
//[CO20220623 - OBSOLETE]  // [OBSOLETE]  sstr=ReScale(sstr,1.0);
//[CO20220623 - OBSOLETE]  sstr.ReScale(1.0);
//[CO20220623 - OBSOLETE]  // Create atom objects for each atom in structure.
//[CO20220623 - OBSOLETE]  xvector<int> ijk(3);ijk.clear();
//[CO20220623 - OBSOLETE]  for(uint iat=0;iat<sstr.atoms.size();iat++) {
//[CO20220623 - OBSOLETE]    _atom a=sstr.atoms.at(iat);
//[CO20220623 - OBSOLETE]    a.name=sstr.atoms.at(iat).name;
//[CO20220623 - OBSOLETE]    a.basis=iat; //[CO20200130 - number->basis]a.number=iat;
//[CO20220623 - OBSOLETE]    a.ijk=sstr.atoms.at(iat).ijk;
//[CO20220623 - OBSOLETE]    a.cpos=sstr.atoms.at(iat).cpos;
//[CO20220623 - OBSOLETE]    a.fpos=sstr.atoms.at(iat).fpos;//cerr << sstr.atoms.at(iat).fpos << endl;
//[CO20220623 - OBSOLETE]    a.type=sstr.atoms.at(iat).type;
//[CO20220623 - OBSOLETE]    atom_vec.push_back(a);
//[CO20220623 - OBSOLETE]  }
//[CO20220623 - OBSOLETE]  double rmin=1e-6;
//[CO20220623 - OBSOLETE]  // [OBSOLETE] GetNeighData(atom_vec,sstr,rmin,cutoff,neigh_mat);
//[CO20220623 - OBSOLETE]  sstr.GetNeighData_20220101(atom_vec,rmin,cutoff,neigh_mat);
//[CO20220623 - OBSOLETE]}

// **************************************************************************
// Function GetBasisTransformation //DX20201015
// **************************************************************************
xmatrix<double> GetBasisTransformation(const xmatrix<double>& lattice_original, const xmatrix<double>& lattice_new) {
  return lattice_new*inverse(lattice_original);
}

// **************************************************************************
// Function GetBasisTransformationInternalTranslations //DX20201124
// **************************************************************************
vector<xvector<double> > GetBasisTransformationInternalTranslations(const xmatrix<double>& basis_transformation) {

  // Given a basis transformation matrix, determine any internal lattice
  // translation(s). This is necessary if the basis transformation increases
  // the volume of the cell, otherwise, there are no internal translations
  // (return immediately).
  // Another way to think of this: if you expand your lattice/cell, this
  // function finds all the lattice points in the new cell

  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream message;

  vector<xvector<double> > translations;

  //DX20210520 [OBSOLETE] double cell_volume_change = aurostd::abs(aurostd::det(basis_transformation));

  //DX20210520 [OBSOLETE] if(LDEBUG){ cerr << __AFLOW_FUNC__ << " changed in cell volume from basis transformation: " << cell_volume_change << endl; }

  // ---------------------------------------------------------------------------
  // check if the basis transformation makes the cell larger and find
  // corresponding internal translations
  // DX20210520 - DO NOT EXCLUDE BASED ON DETERMINANT
  // It is possible to stretch in one direction and compress in another and still
  // get determinant=1 (e.g., POCC structures)
  //DX20210520 [OBSOELTE] if(cell_volume_change-1.0>_AUROSTD_XSCALAR_TOLERANCE_INTEGER_){}

  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " cell size increases. Finding internal translations." << endl; }

  // ---------------------------------------------------------------------------
  // get inverse matrix (Q)
  xmatrix<double> inverse_transform = aurostd::inverse(basis_transformation);

  // ---------------------------------------------------------------------------
  // to get translations take the "larger cell" in fractional coordinates
  // and perform the inverse operation (Q) to see how small it gets,
  // then these are the internal translations
  xmatrix<double> lattice_frac = aurostd::eye<double>(3,3);
  xmatrix<double> lattice_shrink = inverse_transform*lattice_frac;

  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " shrunken lattice: " << lattice_shrink << endl; }

  // ---------------------------------------------------------------------------
  // Now that we have the shortest internal translations from lattice shrink
  // (forms a basis), we need to find all the internal translations inside this
  // cell via linear combinations of this basis.
  // To determine how many combinations we need (i.e. how far to expand), we can
  // use LatticeDimensionSphere(). Since lattice_shrink is in fractional
  // coordinates, we need to find the necessary dimensions in each direction
  // to fill the cell (i.e., the unit box). //DX20210111
  xvector<int> dims=LatticeDimensionSphere(lattice_shrink,1.0);
  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " number of times to apply each internal translation: " << dims[1] << "," << dims[2] << "," << dims[3] << endl; }

  // ---------------------------------------------------------------------------
  // create all linear combinations of translations, filter out duplicates later
  xvector<double> a_vec=lattice_shrink(1);
  xvector<double> b_vec=lattice_shrink(2);
  xvector<double> c_vec=lattice_shrink(3);
  xvector<double> a_vec_scaled, b_vec_scaled, c_vec_scaled;
  for(int a=0;a<=dims[1];a++){ //DX20210506 - need <=
    a_vec_scaled = (double)a*a_vec;
    translations.push_back(a_vec_scaled);
    for(int b=0;b<=dims[2];b++){ //DX20210506 - need <=
      b_vec_scaled = (double)b*b_vec;
      translations.push_back(b_vec_scaled);
      translations.push_back(a_vec_scaled+b_vec_scaled);
      for(int c=0;c<=dims[3];c++){ //DX20210506 - need <=
        c_vec_scaled = (double)c*c_vec;
        translations.push_back(c_vec_scaled);
        translations.push_back(a_vec_scaled+c_vec_scaled);
        translations.push_back(b_vec_scaled+c_vec_scaled);
        translations.push_back(a_vec_scaled+b_vec_scaled+c_vec_scaled);
      }
    }
  }

  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " # translations:" << translations.size() << endl;
    for(uint t=0;t<translations.size();t++){
      cerr << __AFLOW_FUNC__ << " translations:" << translations[t] << endl;
    }
  }

  // ---------------------------------------------------------------------------
  // filter out unique translations 
  vector<xvector<double> > unique_translations;
  bool unique = true;
  for(uint t=0;t<translations.size();t++){
    xvector<double> translation_incell = BringInCell(translations[t]);
    unique = true;
    for(uint u=0;u<unique_translations.size() && unique ;u++){
      unique = !(aurostd::isequal(translation_incell,unique_translations[u]));
    }
    if(unique){ unique_translations.push_back(translation_incell); }
  }

  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " # unique_translations:" << unique_translations.size() << endl;
    for(uint t=0;t<unique_translations.size();t++){
      cerr << __AFLOW_FUNC__ << " unique_translations:" << unique_translations[t] << endl;
    }
  }
  translations = unique_translations;

  // ---------------------------------------------------------------------------
  // if the cell size remains the same or shrinks, no internal translations
  //DX20210520 [OBSOELTE] else{}
  //DX20210520 [OBSOELTE]  // use null vector
  //DX20210520 [OBSOELTE]  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " cell size remains the same or reduced. No internal translations." << endl; }
  //DX20210520 [OBSOELTE]  xvector<double> zero_xvector;
  //DX20210520 [OBSOELTE]  translations.push_back(zero_xvector);
  return translations;
}

// **************************************************************************
// Function GetRotation //DX20201015
// **************************************************************************
xmatrix<double> GetRotation(const xmatrix<double>& lattice_original, const xmatrix<double>& lattice_new) {
  return aurostd::inverse(lattice_original)*lattice_new; 
}

// **************************************************************************
// Function ChangeBasis() //DX20201015
// **************************************************************************
// Convert a structure (lattice and atom positions) into a new representation
// based on the input transformation matrix.
// The transformation matrix is generally NOT a unitary transformation -
// it can change the volume of the cell - otherwise it would be a rotation
// (use Rotate() instead).
// The procedure is generalized for transformations that enlarge (supercell)
// or reduce (primitivize) the structure.
// Enlarging the cell: search for unique internal translations based on 
// transformation matrix.
// Reducing the cell: remove duplicate atom positions.
// Example transformation matrix (4x1x1 supercell expansion):
//   -4.0000e+00  0.0000e+00  0.0000e+00
//   -1.0000e+00  0.0000e+00  1.0000e+00
//   -1.0000e+00  1.0000e+00  0.0000e+00

// ---------------------------------------------------------------------------
// returns new xstructure (makes a copy)
xstructure ChangeBasis(const xstructure& xstr, const xmatrix<double>& transformation_matrix) {
  xstructure xstr_transformed = xstr;
  xstr_transformed.ChangeBasis(transformation_matrix);
  return xstr_transformed;
}

// ---------------------------------------------------------------------------
// modifies in-place (efficient)
void xstructure::ChangeBasis(const xmatrix<double>& transformation_matrix) {

  // if the transformation matrix is the identity, don't do anything
  if(aurostd::isidentity(transformation_matrix)){ return; }

  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream message;

  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " structure BEFORE basis transformation:" << endl;
    cerr << (*this) << endl;
  }

  uint natoms_orig = (*this).atoms.size();
  uint natoms_transformed = 0;
  bool is_integer_multiple_transformation = true;
  if((*this).dist_nn_min==AUROSTD_NAN || (*this).dist_nn_min==AUROSTD_MAX_DOUBLE){ (*this).dist_nn_min=SYM::minimumDistance((*this)); }

  // if the cell is being reduced in a symmetry routine, it is best to set the same-atom tolerance to
  // the minimum interatomic distance minus the sym_eps (e.g., dist_nn_min-sym_eps), since this is the "new resolution"
  // of the atom postions; 90% of dist_nn_min is used so that the nearest-neighbor is not removed erroneously //DX20210716
  double tol=0.9*(*this).dist_nn_min; //DX20210316 - 0.1 is not robust, used fraction of min_dist //DX20210716 - default tol of 90% of dist_nn_min should be sufficient
  if((*this).sym_eps!=AUROSTD_NAN && (*this).sym_eps!=AUROSTD_MAX_DOUBLE){ tol -= (*this).sym_eps; } //DX20210716 - if sym_eps is valid, then subtract it from the tol, since this is the new "resolution" for the system

  // ---------------------------------------------------------------------------
  // transform the lattice 
  xmatrix<double> lattice_orig = (*this).lattice;
  (*this).lattice = transformation_matrix*(*this).lattice;
  (*this).FixLattices();

  // ---------------------------------------------------------------------------
  // get internal translations from basis transformation (i.e. transforming
  // to larger cells)
  vector<xvector<double> > translations = GetBasisTransformationInternalTranslations(transformation_matrix);

  // ---------------------------------------------------------------------------
  // transform the atom positions
  deque<_atom> atom_basis;
  _atom atom_tmp;
  xmatrix<double> forig2fnew=inverse(trasp(transformation_matrix)); // Q*pos , but need to transpose Q for AFLOW xmatrix convention
  xmatrix<double> f2c=trasp((*this).lattice); // Q*pos , but need to transpose Q for AFLOW xmatrix convention
  for(uint i=0;i<(*this).atoms.size();i++){
    for(uint t=0;t<translations.size();t++){
      atom_tmp = (*this).atoms[i];
      atom_tmp.fpos=forig2fnew*((*this).atoms[i].fpos);
      atom_tmp.fpos=::BringInCell(atom_tmp.fpos+translations[t]);
      atom_tmp.cpos=f2c*atom_tmp.fpos;
      if(!SYM::MapAtom(atom_basis, atom_tmp, true, (*this).lattice, false, tol)){ //DX20210324 - to account for duplicate translations
        atom_basis.push_back(atom_tmp);
      }
    }
  }

  // ---------------------------------------------------------------------------
  // calculate change in basis transformation determinant
  // (shift to 1 to easily see if reduces or expands)
  double basis_transformation_det_change = aurostd::abs(aurostd::det(transformation_matrix))-1.0;

  // ---------------------------------------------------------------------------
  // reduce the cell: remove any duplicate atoms
  // use _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_ to be consistent with AUROSTD's
  // isinteger tolerance
  if(basis_transformation_det_change < -_AUROSTD_XSCALAR_TOLERANCE_IDENTITY_){
    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " removing duplicate atoms (cell has been reduced)." << endl; }

    bool skew = false;
    deque<_atom> new_basis = ::foldAtomsInCell(atom_basis, lattice_orig, (*this).lattice, skew, tol, false); //false: don't check atom mappings (slow) //DX20210118 - add global namespace
    atom_basis = new_basis;

    // check atom count
    natoms_transformed = atom_basis.size();
    is_integer_multiple_transformation = ((*this).num_each_type.size() == 1 || natoms_orig%natoms_transformed==0); //DX20210316 - integer multiple does not apply to unaries
  }
  // ---------------------------------------------------------------------------
  // enlarge the cell: update the atom count information
  else if(basis_transformation_det_change > _AUROSTD_XSCALAR_TOLERANCE_INTEGER_){
    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " cell size has increased." << endl; }
    // check atom count
    natoms_transformed = atom_basis.size();
    is_integer_multiple_transformation = ((*this).num_each_type.size() == 1 || natoms_transformed%natoms_orig==0); //DX20210316 - integer multiple does not apply to unaries
  }

  // ---------------------------------------------------------------------------
  // is integer multiple transformation (reduce or enlarge)
  if(!is_integer_multiple_transformation){
    message << "Number of atoms is no longer an integer multiple with respect to the input structure"
      << " original: " << natoms_orig
      << " transformed: " << natoms_transformed
      << "; check the transformation matrix or same-atom tolerance.";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
  }

  // ---------------------------------------------------------------------------
  // if the number of atoms changed (i.e., change in determinant is zero),
  // update the atom counts/order/types/etc.
  if(!aurostd::isequal(aurostd::abs(basis_transformation_det_change), _ZERO_TOL_, _AUROSTD_XSCALAR_TOLERANCE_INTEGER_)){
    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " updating atom count information." << endl; }
    std::stable_sort(atom_basis.begin(),atom_basis.end(),sortAtomsNames); //DX20210129
    (*this).ReplaceAtoms(atom_basis, false); //false: check_atom_overlap
  }
  // ---------------------------------------------------------------------------
  // if the transformation preserves the volume, one-to-one mappings
  // no need to update species/types/etc. (i.e., ReplaceAtoms() is not needed)
  else{
    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " cell size remains the same (updating atom positions)." << endl; }
    (*this).atoms = atom_basis;
  }

  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " structure AFTER basis transformation:" << endl;
    cerr << (*this) << endl;
  }
}

// **************************************************************************
// Function TransformStructure() //DX20201125
// **************************************************************************
// ---------------------------------------------------------------------------
// returns new xstructure (makes a copy)
xstructure TransformStructure(const xstructure& xstr,
    const xmatrix<double>& transformation_matrix,
    const xmatrix<double>& rotation) {

  xvector<double> origin_shift;
  return TransformStructure(xstr, transformation_matrix, rotation, origin_shift);
}

xstructure TransformStructure(const xstructure& xstr,
    const xmatrix<double>& transformation_matrix,
    const xmatrix<double>& rotation,
    const xvector<double>& origin_shift,
    bool is_shift_frac){
  xstructure xstr_transformed = xstr;
  xstr_transformed.TransformStructure(transformation_matrix, rotation, origin_shift, is_shift_frac);
  return xstr_transformed;
}

// ---------------------------------------------------------------------------
// modifies in-place (efficient)
void xstructure::TransformStructure(
    const xmatrix<double>& transformation_matrix,
    const xmatrix<double>& rotation) {

  xvector<double> origin_shift;
  (*this).TransformStructure(transformation_matrix, rotation, origin_shift);
}

void xstructure::TransformStructure(
    const xmatrix<double>& transformation_matrix,
    const xmatrix<double>& rotation,
    const xvector<double>& origin_shift,
    bool is_shift_frac) {

  bool LDEBUG=(FALSE || XHOST.DEBUG);

  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " basis transformation: " << transformation_matrix << endl;
    cerr << __AFLOW_FUNC__ << " rotation (R): " << rotation << endl;
  }

  // ---------------------------------------------------------------------------
  // changed basis
  (*this).ChangeBasis(transformation_matrix);
  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " structure after CHANGING BASIS: " << (*this) << endl; }

  // ---------------------------------------------------------------------------
  // rotate
  (*this).Rotate(rotation);
  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " structure after ROTATING: " << (*this) << endl; }

  // ---------------------------------------------------------------------------
  // rotate
  bool coordinate_flag = (*this).coord_flag; // store original coordinate-type
  (*this).ShiftPos(origin_shift,is_shift_frac);
  (*this).coord_flag=coordinate_flag; // set back to original coordinate-type
  (*this).BringInCell(); //DX20210116
  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " structure after shifting origin: " << (*this) << endl; }
}

// **************************************************************************
// GenerateXXXXXXXXXXx
// **************************************************************************
void xstructure::qm_clear(void) {
  qm_calculated=FALSE;
  qm_scale=0.0;
  qm_lattice.clear();
  qm_origin.clear();
  qm_klattice.clear();
  qm_f2c.clear();
  qm_c2f.clear();
  qm_origin.clear();
  qm_forces_write=FALSE;
  qm_positions_write=FALSE;
  qm_atoms.clear();
  qm_forces.clear();
  qm_positions.clear();
  for(uint iat=0;iat<atoms.size();iat++) {
    qm_atoms.push_back(atoms.at(iat));          // just plug something, then I`ll clean
    qm_atoms.at(iat).cpos.clear();
    qm_atoms.at(iat).fpos.clear();
    qm_forces.push_back(atoms.at(iat).cpos);    // just plug something, then I`ll clean
    qm_forces.at(iat).clear();
    qm_positions.push_back(atoms.at(iat).cpos); // just plug something, then I`ll clean
    qm_positions.at(iat).clear();
  }
  qm_E_cell=0.0;qm_dE_cell=0.0;qm_H_cell=0.0;qm_PV_cell=0.0;qm_mag_cell=0.0;qm_P=0.0;
  qm_E_atom=0.0;qm_dE_atom=0.0;qm_H_atom=0.0;qm_PV_atom=0.0;qm_mag_atom=0.0; 

  if(atoms.size()!=qm_atoms.size())     {throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::qm_clear():","[1] atoms.size()!=qm_atoms.size().",_VALUE_ERROR_);}
  if(atoms.size()!=qm_forces.size())    {throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::qm_clear():","[2] atoms.size()!=qm_forces.size().",_VALUE_ERROR_);}
  if(atoms.size()!=qm_positions.size()) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::qm_clear():","[3] atoms.size()!=qm_positions.size().",_VALUE_ERROR_);}

}

void xstructure::qm_recycle(void) {
  if(atoms.size()!=qm_atoms.size())     {throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::qm_recycle():","[1] atoms.size()!=qm_atoms.size().",_VALUE_ERROR_);}
  if(atoms.size()!=qm_forces.size())    {throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::qm_recycle():","[2] atoms.size()!=qm_forces.size().",_VALUE_ERROR_);}
  if(atoms.size()!=qm_positions.size()) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::qm_recycle():","[3] atoms.size()!=qm_positions.size().",_VALUE_ERROR_);}
  scale=qm_scale;
  lattice=qm_lattice;
  klattice=qm_klattice;
  f2c=qm_f2c;
  c2f=qm_c2f;
  origin=qm_origin;
  FixLattices();
  for(uint i=0;i<atoms.size();i++) {                 // copy cpos/fpos
    atoms.at(i).cpos=qm_atoms.at(i).cpos;     // get cpos
    atoms.at(i).fpos=qm_atoms.at(i).fpos;     // get fpos
  }
  qm_clear();
}

void xstructure::qm_load(const string& Directory,const string& suffix,int iomode) {
  double data_natoms=double(atoms.size());
  if(iomode!=IOVASP_POSCAR) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::qm_load():","Only IOVASP_POSCAR is supported.",_FILE_WRONG_FORMAT_);};
  if(iomode==IOVASP_POSCAR) {
    xOUTCAR outcar;
    if(aurostd::FileEmpty(Directory+"/OUTCAR"+suffix))   {throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::qm_load():","Empty OUTCAR.",_FILE_CORRUPT_);}  //PN+JJPR FIXED BUG

    outcar.GetPropertiesFile(Directory+"/OUTCAR"+suffix); //SD20221019 - don't perform natoms check yet

    if(abs(data_natoms-(double) outcar.natoms)>0.1) {
      stringstream message;
      message << "ERROR void xstructure::qm_load: data_natoms(" << data_natoms << ")!= (int) outcar.natoms(" << outcar.natoms << ") ..." << endl;
      message << "      Directory=" << Directory << endl;
      message << "      suffix=" << suffix << endl;
      message << "      iomode=" << iomode << endl;
      throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::qm_load():",message,_FILE_WRONG_FORMAT_);;
    }

    //    cerr << atoms.size() << endl;
    qm_clear();
    if(atoms.size()!=qm_atoms.size())     {throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::qm_load():","[1] atoms.size()!=qm_atoms.size().",_VALUE_ERROR_);}
    if(atoms.size()!=qm_forces.size())    {throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::qm_load():","[2] atoms.size()!=qm_forces.size().",_VALUE_ERROR_);}
    if(atoms.size()!=qm_positions.size()) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::qm_load():","[3] atoms.size()!=qm_positions.size().",_VALUE_ERROR_);}

    // NEW WITH xOUTCAR
    qm_forces.clear(); for(uint i=0;i<outcar.vforces.size();i++)  qm_forces.push_back(outcar.vforces.at(i)); 
    qm_positions.clear(); for(uint i=0;i<outcar.vpositions_cartesian.size();i++)  qm_positions.push_back(outcar.vpositions_cartesian.at(i)); 
    if(atoms.size()!=qm_forces.size())    {throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::qm_load():","[4] atoms.size()!=qm_forces.size().",_VALUE_ERROR_);}
    if(atoms.size()!=qm_positions.size()) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::qm_load():","[5] atoms.size()!=qm_positions.size().",_VALUE_ERROR_);}

    // NEW WITH xVASPRUNXML
    xVASPRUNXML vasprunxml;
    if(aurostd::FileEmpty(Directory+"/vasprun.xml"+suffix))   {throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::qm_load():","Empty vasprun.xml.",_FILE_CORRUPT_);} //PN+JJPR FIXED BUG
    //vasprunxml.GetPropertiesFile(Directory+"/vasprun.xml");
    vasprunxml.GetPropertiesFile(Directory+"/vasprun.xml"+suffix); //PN+JJPR FIXED BUG
    qm_forces.clear(); for(uint i=0;i<vasprunxml.vforces.size();i++)  qm_forces.push_back(vasprunxml.vforces.at(i)); 


    // OLD
    stringstream CONTCAR_stringstream;
    CONTCAR_stringstream.str(std::string()); CONTCAR_stringstream << aurostd::FileToString(Directory+"/CONTCAR"+suffix);

    qm_E_cell=outcar.energy_cell;qm_E_atom=outcar.energy_atom;
    qm_H_cell=outcar.enthalpy_cell;qm_H_atom=outcar.enthalpy_atom;
    qm_PV_cell=outcar.PV_cell;qm_PV_atom=outcar.PV_atom;
    qm_P=outcar.pressure;
    qm_mag_cell=outcar.mag_cell;qm_mag_atom=outcar.mag_atom;

    // CONTCAR OPERATIONS ---------------------------------------------------------------
    // cerr << "DEBUG: xstructure::qm_clear(void)" << endl;
    // cerr << CONTCAR_stringstream.str() << endl;
    // cerr << "DEBUG: xstructure::qm_clear(void)" << endl;
    // string CONTCAR_string="";
    xstructure b(CONTCAR_stringstream);                               // LOAD it in
    qm_scale=b.scale;
    qm_lattice=b.lattice;
    qm_klattice=b.klattice;
    qm_f2c=b.f2c;
    qm_c2f=b.c2f;
    qm_origin=b.origin;
    for(uint i=0;i<atoms.size();i++) {                 // copy cpos/fpos
      qm_atoms.at(i)=atoms.at(i);
      qm_atoms.at(i).cpos=b.atoms.at(i).cpos;     // get cpos
      qm_atoms.at(i).fpos=b.atoms.at(i).fpos;     // get fpos
    }
    qm_calculated=TRUE;
  }
}

// **************************************************************************
// should be in xlibs
// SEARCH SORT WORLD "sort for vector/deque of double_int"

// namespace aurostd {
//   void sort(vector<double>& varg1,vector<int>& varg2) {
//     vector<aurostd::_double_int_> vv(varg1.size());
//     for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);}
//     sort(vv.begin(),vv.end(),_sort_double_int_());
//     for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;}
//   }
//   void sort(deque<double>& varg1,deque<int>& varg2) {
//     deque<aurostd::_double_int_> vv(varg1.size());
//     for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);}
//     sort(vv.begin(),vv.end(),_sort_double_int_());
//     for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;}
//   }
// }

// namespace aurostd {
//   void sort(vector<double>& varg1,vector<double>& varg2) {
//     vector<aurostd::_double_double_> vv(varg1.size());
//     for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);}
//     sort(vv.begin(),vv.end(),_sort_double_double_());
//     for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;}
//   }
//   void sort(deque<double>& varg1,deque<double>& varg2) {
//     deque<aurostd::_double_double_> vv(varg1.size());
//     for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);}
//     sort(vv.begin(),vv.end(),_sort_double_double_());
//     for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;}
//   }
// }

// **************************************************************************
// PAULING DETECTOR
// **************************************************************************
bool PAULING_WyckoffDetector(vector<string> &vinput) {
  vector<string> tokens,elements;
  aurostd::string2tokens(vinput.at(0),tokens);
  for(uint i=1;i<tokens.size();i++)
    elements.push_back(tokens.at(i));
  for(uint i=0;i<elements.size();i++) {
    //  cout << "scanning " << elements.at(i) << endl;
    for(uint j=1;j<vinput.size();j++) {
      aurostd::string2tokens(vinput[j],tokens);
      if(tokens.size()==7) {
        if(tokens.at(1)==elements.at(i)) {
          cout << "   ";
          for(uint k=0;k<3;k++) cout << aurostd::PaddedPOST(tokens.at(4+k),7," ") << "  ";
          cout << char(i+65) << "   "; // shift to ascii
          cout << aurostd::PaddedPRE(tokens.at(0),7," ") << "  ";
          for(uint k=1;k<4;k++) cout << aurostd::PaddedPOST(tokens.at(k),4," ") << "  ";
          cout << endl;
        }
      }
    }
  }

  throw aurostd::xerror(__AFLOW_FILE__,XPID+"PAULING_WyckoffDetector():","Throw for debugging purposes.",_GENERIC_ERROR_);
}

//DX20170831 - xstructure2json - START
// **************************************************************************
// xstructure2json
// **************************************************************************
string xstructure2json(xstructure& xstr) {
  string eendl="";
  bool roff=true; //round off
  stringstream sss;
  stringstream sscontent_json;
  vector<string> vcontent_json;

  // TITLE 
  if(xstr.title.size()){
    sscontent_json << "\"title\":\"" << aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xstr.title) << "\"" << eendl; //CO20220531 - removing whitespace
  } else {
    if(PRINT_NULL_JSON){ sscontent_json << "\"title\":null" << eendl;}
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // SCALE/VOL 
  if(xstr.scale){
    sscontent_json << "\"scale\":" << xstr.scale << "" << eendl; //DX20180306 - number not string (removed quotations)
  } else {
    if(PRINT_NULL_JSON){ sscontent_json << "\"scale\":null" << eendl;}
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // LATTICE
  if(xstr.lattice.rows){
    sscontent_json << "\"lattice\":[" << aurostd::xmatDouble2String(xstr.lattice,_DOUBLE_WRITE_PRECISION_MAX_,roff) << "]" << eendl; //CO20180515
  } else {
    if(PRINT_NULL_JSON){ sscontent_json << "\"lattice\":null" << eendl;}
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // SPECIES 
  if(xstr.species.size()){
    deque<string> cleaned_species; //DX20190612 - cleaned species names
    for(uint i=0;i<xstr.species.size(); i++) { cleaned_species.push_back(KBIN::VASP_PseudoPotential_CleanName(xstr.species[i])); } //DX20190612 - cleaned species names
    sscontent_json << "\"species\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(cleaned_species,"\""),",") << "]" << eendl; //DX20190612 - cleaned species names
  } else {
    if(PRINT_NULL_JSON){ sscontent_json << "\"species\":null" << eendl;}
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // NUMBER OF TYPES
  if(xstr.num_each_type.size()){
    sscontent_json << "\"number_each_type\":[" << aurostd::joinWDelimiter(xstr.num_each_type,",") << "]" << eendl;
  } else {
    if(PRINT_NULL_JSON){ sscontent_json << "\"number_each_type\":null" << eendl;}
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // COORDINATES TYPE
  if(xstr.coord_flag==_COORDS_FRACTIONAL_){
    sscontent_json << "\"coordinates_type\":\"direct\"" << eendl; //CO20171025
  } else if(xstr.coord_flag==_COORDS_CARTESIAN_){
    sscontent_json << "\"coordinates_type\":\"Cartesian\"" << eendl;
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // ATOMS
  if(xstr.atoms.size()){
    sscontent_json << "\"atoms\":[" << eendl;
    for(uint i=0;i<xstr.atoms.size();i++){
      sscontent_json << atom2json(xstr.atoms[i],xstr.coord_flag,xstr.partial_occupation_flag) << eendl;
      if(i != xstr.atoms.size()-1){
        sscontent_json << ",";
      }
    }
    sscontent_json << "]" << eendl;
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  sss << "{" << aurostd::joinWDelimiter(vcontent_json,",")  << "}" << eendl;
  return sss.str();
}
//DX20170831 - xstructure2json - END

//DX20170831 - atom2json - START
// **************************************************************************
// atom2json
// **************************************************************************
string atom2json(_atom& atom, int coord_flag, int poccupation) {
  string eendl="";
  bool roff=true; //round off
  stringstream sss;
  stringstream sscontent_json;
  vector<string> vcontent_json;

  // NAME 
  if(atom.name.size()){
    sscontent_json << "\"name\":\"" << KBIN::VASP_PseudoPotential_CleanName(atom.name) << "\"" << eendl; //DX20190612 - added function to clean names
  } else {
    if(PRINT_NULL_JSON){ sscontent_json << "\"name\":null" << eendl;}
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // POSITION
  if(coord_flag==_COORDS_FRACTIONAL_){
    sscontent_json << "\"position\":[" << aurostd::joinWDelimiter(xvecDouble2vecString(atom.fpos,_DOUBLE_WRITE_PRECISION_MAX_,roff),",") << "]" << eendl; //CO20180515
  } else if(coord_flag==_COORDS_CARTESIAN_){
    sscontent_json << "\"position\":[" << aurostd::joinWDelimiter(xvecDouble2vecString(atom.cpos,_DOUBLE_WRITE_PRECISION_MAX_,roff),",") << "]" << eendl; //CO20180515
  } else {
    if(PRINT_NULL_JSON){ sscontent_json << "\"position\":null" << eendl;}
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // OCCUPANCY
  if(poccupation==TRUE){
    sscontent_json << "\"occupancy\":" << atom.partial_occupation_value << eendl;
  } else if(poccupation==FALSE){
    sscontent_json << "\"occupancy\":1.0"<< eendl;
  } else {
    if(PRINT_NULL_JSON){ sscontent_json << "\"occupancy\":null" << eendl;}
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  sss << "{" << aurostd::joinWDelimiter(vcontent_json,",")  << "}";
  return sss.str();
} 
//DX20170831 - atom2json

#endif  // _AFLOW_XATOM_CPP

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2023              *
// *                                                                        *
// **************************************************************************
