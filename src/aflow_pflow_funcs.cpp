// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo
// Dane Morgan

#ifndef _AFLOW_PFLOW_FUNCS_CPP_
#define _AFLOW_PFLOW_FUNCS_CPP_

#include "aflow_pflow.h"

// ***************************************************************************
// XRAY XRAY XRAY XRAY XRAY XRAY XRAY XRAY XRAY XRAY XRAY XRAY XRAY XRAY XRAY
// ***************************************************************************
// Function DebyeWallerFactor
// ***************************************************************************
// All data must be sent in SI units.
// This formula is from B.E. Warren, X-ray Diffraction, eq. 5.9.
// It is based on the Debye approximation and only holds for T>>TDebye.
// Note that this function returns the sqrt of what we usually call the
// DW factor.  This is because this term multiplies the scattering
// factors f, not f^2.  DW^2 would be the appropriate term to modulate
// the intensity, but DW is appropriate to modulate the scattering factors
// (see Warren,eq.3.24).
double DebyeWallerFactor(double theta,
    double temp, double debye_temp,
    double mass,double lambda) {
  double st=sin(theta);
  double h=PLANCKSCONSTANT_h; //ME20181020
  double twoB=h*h*temp*12.0/(mass*KBOLTZ*debye_temp*debye_temp);
  double twoM=twoB*st*st/(lambda*lambda);
  double DWfact=exp(-1.0*twoM/2.0); // Use e^-M, not e^-2M.
  return DWfact;
}

// ***************************************************************************
// getGenericTitleXStructure()
// ***************************************************************************
string getGenericTitleXStructure(const xstructure& xstr,bool latex){ //CO20190520
  //CO20200624 - used to be num_each_type, now it's comp_each_type (works for both pocc and non-pocc)
  //use pocc default for precision
  string title="";
  uint iat=0;
  int comp_prec=(int)ceil(log10(1.0/xstr.partial_occupation_stoich_tol));  //ceil ensures we round up above 1 //CO20181226
  bool atom_names=true;

  for(uint i=0;i<xstr.atoms.size()&&atom_names;i++){if(xstr.atoms[i].cleanname.empty()){atom_names=false;}}
  for(uint itype=0;itype<xstr.num_each_type.size();itype++){
    for(uint j=0;j<(uint)xstr.num_each_type[itype];j++) {
      if(j==0){
        if(atom_names){title+=xstr.atoms[iat].cleanname;} //CO20200624 - never use pp names, never mix pp with composition
        else {title+=char('A'+itype);}
        if(latex){title+="$_{"+aurostd::utype2string(xstr.comp_each_type[itype],comp_prec)+"}$";}
        else {title+=aurostd::utype2string(xstr.comp_each_type[itype],comp_prec);}
      }
      iat++;
    }
  }
  return title;
}

// ***************************************************************************
// balanceEquation()
// ***************************************************************************
xvector<double> balanceChemicalEquation(const vector<xvector<double> >& lhs,const vector<xvector<double> >& rhs,
    bool normalize,double tol){ //CO20180801
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  xvector<double> dummy;
  if(lhs.size()==0 || rhs.size()==0){return dummy;}
  int dim=lhs[0].rows;
  for(uint i=1;i<lhs.size();i++){if(lhs[i].rows!=dim){return dummy;}}
  for(uint i=0;i<rhs.size();i++){if(rhs[i].rows!=dim){return dummy;}}

  if(LDEBUG) {
    cerr << __AFLOW_FUNC__ << " lhs=";
    for(uint i=0;i<lhs.size();i++){
      cerr << lhs[i] << (i!=lhs.size()-1?", ":"");
    }
    cerr << endl;
    cerr << __AFLOW_FUNC__ << " rhs=";
    for(uint i=0;i<rhs.size();i++){
      cerr << rhs[i] << (i!=rhs.size()-1?", ":"");
    }
    cerr << endl;
  }

  //[OBSOLETE CO20180801]vector<xvector<double> > lhs;
  //[OBSOLETE CO20180801]vector<xvector<double> > rhs;
  //[OBSOLETE CO20180801]if(reduce){
  //[OBSOLETE CO20180801]  for(uint i=0;i<_lhs.size();i++){lhs.push_back(aurostd::reduceByGCD(_lhs[i],tol));}
  //[OBSOLETE CO20180801]  for(uint i=0;i<_rhs.size();i++){rhs.push_back(aurostd::reduceByGCD(_rhs[i],tol));}
  //[OBSOLETE CO20180801]  if(LDEBUG) {
  //[OBSOLETE CO20180801]    cerr << __AFLOW_FUNC__ << " lhs REDUCED=";
  //[OBSOLETE CO20180801]    for(uint i=0;i<lhs.size();i++){
  //[OBSOLETE CO20180801]      cerr << lhs[i] << (i!=lhs.size()-1?", ":"");
  //[OBSOLETE CO20180801]    }
  //[OBSOLETE CO20180801]    cerr << endl;
  //[OBSOLETE CO20180801]    cerr << __AFLOW_FUNC__ << " rhs REDUCED=";
  //[OBSOLETE CO20180801]    for(uint i=0;i<rhs.size();i++){
  //[OBSOLETE CO20180801]      cerr << rhs[i] << (i!=rhs.size()-1?", ":"");
  //[OBSOLETE CO20180801]    }
  //[OBSOLETE CO20180801]    cerr << endl;
  //[OBSOLETE CO20180801]  }
  //[OBSOLETE CO20180801]} else {lhs=_lhs;rhs=_rhs;}

  // needs to organized in the following way
  // [[Mn=2,Cu=1,Fe=0],   //left_hand_side
  // [Mn=0,Cu=5,Fe=3],    //left_hand_side
  // [Mn=-1,Cu=-1,Fe=0],  //right_hand_side
  // [Mn=-2,Cu=-1,Fe=0]]  //right_hand_side
  // i.e. compounds in rows, elements in cols
  // left_hand_side of reaction is POSITIVE
  // right_hand_side of reaction is NEGATIVE

  xmatrix<double> composition_matrix(lhs.size()+rhs.size(),dim);
  int counter;
  for(uint i=0;i<lhs.size();i++){counter=1;for(int j=lhs[i].lrows;j<=lhs[i].urows;j++){composition_matrix(i+1,counter++)=lhs[i][j];}}             //lhs
  for(uint i=0;i<rhs.size();i++){counter=1;for(int j=rhs[i].lrows;j<=rhs[i].urows;j++){composition_matrix(i+1+lhs.size(),counter++)=-rhs[i][j];}} //rhs

  if(LDEBUG) {
    cerr << __AFLOW_FUNC__ << " composition matrix:" << endl;
    cerr << composition_matrix << endl;
  }

  return balanceChemicalEquation(composition_matrix,normalize,tol);
}

//normalize === set first coefficient to 1
xvector<double> balanceChemicalEquation(const xmatrix<double>& composition_matrix,bool normalize,double tol){ //CO20191110
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream message;
  if(composition_matrix.rows<=composition_matrix.cols){
    message << "Composition matrix (m<=n) will NOT yield a viable null space for this analysis";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_); //CO20190226
  }
  //[CO20191110 - OBSOLETE]xmatrix<double> composition_matrix=_composition_matrix;
  //[CO20191110 - OBSOLETE]xmatrix<double> Q=aurostd::generalHouseHolderQRDecomposition(composition_matrix);
  xmatrix<double> Q,R;
  aurostd::QRDecomposition_HouseHolder(composition_matrix,Q,R);
  if(LDEBUG) {
    cerr << __AFLOW_FUNC__ << " Q:" << endl;
    cerr << Q << endl;
  }
  xvector<double> coef(R.rows);
  for(int i=1;i<Q.rows+1;i++){coef[i]=Q(i,Q.ucols);}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " PRE-coefficients =" << coef << endl;}
  if(normalize){
    if(abs(coef[1])<tol){return coef;} //avoid divide by 0
    coef/=coef[1];
  }
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " POST-coefficients=" << coef << endl;}
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " checking: all coef's are >0" << endl;}
  //check that none are negative, that means we did NOT provide a valid equation (check products)
  for(int i=1;i<coef.rows+1;i++){
    if(coef[i]<-_ZERO_TOL_){
      message << "Found negative coef[" << i << "], invalid chemical equation (run with --debug to see): coef=" << coef;
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_RANGE_); //CO20190226
    }
  }
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " checking: all scalar products should yield 0" << endl;}
  double sum;
  for(int i=1;i<composition_matrix.cols+1;i++){
    if(LDEBUG) {cerr << "component-" << i << ": sum=";}
    sum=0.0;
    for(int j=1;j<composition_matrix.rows+1;j++){sum+=coef(j)*composition_matrix(j,i);}
    if(LDEBUG) {cerr << sum << endl;}
    if(abs(sum)>_ZERO_TOL_){
      message << "Chemical equation was not balanced (run with --debug to see): coef=" << coef;
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_RANGE_); //CO20190226
    }
  }
  return coef;
}

namespace pflow {
  //CO20190321
  //follows procedure outlined in: De Leon et al., PRL 114, 165502 (2015) (supp info)
#define DEFAULT_SURFACE_LAYERS 3  //See W. Sun and G. Ceder, Surface Science 617 (2013) 53-59, fix/relax these layers of the slab
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow,istream& input,ostream& oss){
    xstructure xstr_in(input,IOAFLOW_AUTO);
    return GeneralizedStackingFaultEnergyCalculation(vpflow,xstr_in,oss);
  }
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow,const xstructure& xstr_in,ostream& oss){
    _aflags aflags; aflags.Directory=".";
    _kflags kflags;
    _vflags vflags;
    return GeneralizedStackingFaultEnergyCalculation(vpflow,xstr_in,aflags,kflags,vflags,oss);
  }
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow,istream& input,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ostream& oss){
    xstructure xstr_in(input,IOAFLOW_AUTO);
    return GeneralizedStackingFaultEnergyCalculation(vpflow,xstr_in,aflags,kflags,vflags,oss);
  }
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow,const xstructure& xstr_in,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ostream& oss){
    ofstream FileMESSAGE;
    return GeneralizedStackingFaultEnergyCalculation(vpflow,xstr_in,aflags,kflags,vflags,FileMESSAGE,oss);
  }
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow,istream& input,ofstream& FileMESSAGE,ostream& oss){
    xstructure xstr_in(input,IOAFLOW_AUTO);
    return GeneralizedStackingFaultEnergyCalculation(vpflow,xstr_in,FileMESSAGE,oss);
  }
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow,const xstructure& xstr_in,ofstream& FileMESSAGE,ostream& oss){
    _aflags aflags; aflags.Directory=".";
    _kflags kflags;
    _vflags vflags;
    return GeneralizedStackingFaultEnergyCalculation(vpflow,xstr_in,aflags,kflags,vflags,FileMESSAGE,oss);
  }
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow,istream& input,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss){
    xstructure xstr_in(input,IOAFLOW_AUTO);
    return GeneralizedStackingFaultEnergyCalculation(vpflow,xstr_in,aflags,kflags,vflags,FileMESSAGE,oss);
  }
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow,const xstructure& xstr_in,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss){
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;
    std::streamsize prec = 8;
    bool check_min_dist=true; //turn off if it gets too slow
    int count_check_min_dist=0;

    message << aflow::Banner("BANNER_NORMAL");
    pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_RAW_);  //first to screen (not logged, file not opened)
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " starting" << endl;}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - get conventional cell
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " get conventional cell" << endl;}

    xstructure xstr_bulk(xstr_in);
    double min_dist=xstr_bulk.dist_nn_min;
    if(min_dist==AUROSTD_NAN){min_dist=SYM::minimumDistance(xstr_bulk);}
    double min_dist_orig=min_dist;
    if(check_min_dist){ //sanity check as we rotate structure/atoms
      min_dist=xstr_bulk.MinDist();
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " mindist[" << count_check_min_dist++ << "]=" << min_dist << endl;}
      if(!aurostd::isequal(min_dist_orig,min_dist)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Minimum distance changed",_VALUE_ERROR_);}
    }

    bool convert_sconv=true;
    if(convert_sconv){xstr_bulk=Standard_Conventional_UnitCellForm(xstr_bulk);xstr_bulk.clean();}  //best to work with standard conventional unitcell //DX20191220 - uppercase to lowercase clean
    if(check_min_dist){ //sanity check as we rotate structure/atoms
      min_dist=xstr_bulk.MinDist();
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " mindist[" << count_check_min_dist++ << "]=" << min_dist << endl;}
      if(!aurostd::isequal(min_dist_orig,min_dist)){
        //throw a warning here instead, minimum distance MIGHT change with sconv
        //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Minimum distance changed",_VALUE_ERROR_);
        message << "Minimum distance changed (sprim -> sconv)";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_WARNING_);
        min_dist_orig=min_dist;
      }
    }
    xstr_bulk.ReScale(1.0);
    xstr_bulk.ShiftOriginToAtom(0);xstr_bulk.origin=0.0; //reset origin
    xstr_bulk.BringInCell();
    xstr_bulk.clean();  //clear origin! //DX20191220 - uppercase to lowercase clean
    if(LDEBUG) {xstr_bulk.write_DEBUG_flag=TRUE;}
    //xstr_bulk.coord_flag=_COORDS_CARTESIAN_;  //much more accurate for this type of calculation

    message << "structure(standard conventional)=";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
    message << xstr_bulk << endl;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_RAW_);

    if(check_min_dist){ //sanity check as we rotate structure/atoms
      min_dist=xstr_bulk.MinDist();
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " mindist[" << count_check_min_dist++ << "]=" << min_dist << endl;}
      if(!aurostd::isequal(min_dist_orig,min_dist)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Minimum distance changed",_VALUE_ERROR_);}
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - get conventional cell
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - get sym info for structure (mostly for FPOSMatch)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    xstructure xstr_sym(xstr_bulk);  //make a copy so we don't carry around all the symmetry in memory as we make copies of a
    pflow::PerformFullSymmetry(xstr_sym);
    xstr_bulk.dist_nn_min=xstr_sym.dist_nn_min;
    xstr_bulk.sym_eps=xstr_sym.sym_eps;
    bool skew=SYM::isLatticeSkewed(xstr_bulk.lattice,xstr_bulk.dist_nn_min,xstr_bulk.sym_eps);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " xstr_bulk.dist_nn_min=" << xstr_bulk.dist_nn_min << endl;
      cerr << __AFLOW_FUNC__ << " xstr_bulk.sym_eps=" << xstr_bulk.sym_eps << endl;
      cerr << __AFLOW_FUNC__ << " skew=" << skew << endl;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - get sym info for structure (mostly for FPOSMatch)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - read flags
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " reading flags" << endl;}

    int h_s=1,k_s=1,l_s=0;  //hkl of shear
    xvector<int> hkl_s,hkl_s_ORIG;
    double step_size=0.2;   //step size
    int fixed_layers=DEFAULT_SURFACE_LAYERS;     //number of fixed layers
    bool spin_off=false;
    bool partial_dissociation=false;

    vector<string> tokens;
    aurostd::string2tokens(vpflow.getattachedscheme("GENERALIZED_STACKING_FAULT_ENERGY::SHEAR_DIRECTION"),tokens,",");
    if(tokens.size()==3){
      h_s=aurostd::string2utype<int>(tokens[0]);
      k_s=aurostd::string2utype<int>(tokens[1]);
      l_s=aurostd::string2utype<int>(tokens[2]);
    }
    hkl_s[1]=h_s;hkl_s[2]=k_s;hkl_s[3]=l_s;
    if(hkl_s[1]==0 && hkl_s[2]==0 && hkl_s[3]==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"hkl_s=(0,0,0)",_INPUT_ERROR_);}
    hkl_s_ORIG=hkl_s;
    string step_size_string=vpflow.getattachedscheme("GENERALIZED_STACKING_FAULT_ENERGY::STEP_SIZE"); //step size
    if(aurostd::isfloat(step_size_string)){
      double _step_size=aurostd::string2utype<double>(step_size_string);
      if(_step_size>0.0 && _step_size<1.0){step_size=_step_size;}
    }
    string steps_string=vpflow.getattachedscheme("GENERALIZED_STACKING_FAULT_ENERGY::STEPS"); //step size
    if(aurostd::isfloat(steps_string)){
      int _steps=aurostd::string2utype<int>(steps_string);
      if(_steps>0){step_size=1.0/_steps;}
    }
    string fixed_layers_string=vpflow.getattachedscheme("GENERALIZED_STACKING_FAULT_ENERGY::FIXED_LAYERS");
    if(aurostd::isfloat(fixed_layers_string)){
      int _fixed_layers=aurostd::string2utype<int>(fixed_layers_string);
      if(_fixed_layers>0){fixed_layers=_fixed_layers;}
    }
    spin_off=vpflow.flag("GENERALIZED_STACKING_FAULT_ENERGY::SPIN_OFF");
    partial_dissociation=vpflow.flag("GENERALIZED_STACKING_FAULT_ENERGY::PARTIAL_DISSOCIATION");

    std::streamsize prec_original = message.precision(); //original
    std::ios_base::fmtflags ff_original = message.flags();  //original
    message.precision(prec);
    message.unsetf(std::ios_base::floatfield);

    message << "shear_direction" << hkl_s;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
    message << "step_size=" << step_size;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
    message << "fixed_layers=" << fixed_layers;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
    message << "spin=" << (spin_off?"OFF":"ON");pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
    message << "partial_dissociation=" << (partial_dissociation?"ON":"OFF");pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);

    message.precision(prec_original); //set back
    message.flags(ff_original); //set back

    aurostd::xoption slab_flags=vpflow; //may modify based on how many layers in unit cell

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - read flags
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - create slab
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    message << "Creating slab";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);

    xvector<int> hkl_i;
    int total_layers;
    xstructure xstr_slab_newbasis;  //xstr_rotated
    vector<int> sc2pcMap_slab,pc2scMap_slab;

    //create define for rigidrotation vs. slabbasis
    //convert n_s (fractional) to cartesian then fraction of new basis (S=CB)

    xstructure xstr_slab;
    xmatrix<double> rotation; //=xstr_slab_newbasis.lattice*inverse(xstr_bulk.lattice); //xstr_slab_newbasis.lattice = rotation * xstr_bulk.lattice
    //[CO20190805 - this does NOT work]if(0){
    //[CO20190805 - this does NOT work]  xstr_slab=slab::CreateSlab_RigidRotation(slab_flags,xstr_bulk,hkl_i,total_layers,rotation,xstr_slab_newbasis,sc2pcMap_slab,pc2scMap_slab,aflags,FileMESSAGE,oss);
    //[CO20190805 - this does NOT work]}else{
    xstr_slab=slab::CreateSlab_SurfaceLattice(slab_flags,xstr_bulk,hkl_i,total_layers,rotation,xstr_slab_newbasis,sc2pcMap_slab,pc2scMap_slab,aflags,FileMESSAGE,AUROSTD_MAX_DOUBLE,oss);
    //[CO20190805 - this does NOT work]}
    //cerr << rotation << endl << endl;
    //rotation=inverse(xstr_bulk.lattice)*xstr_slab_newbasis.lattice;
    //cerr << rotation << endl << endl;

    //xmatrix<double> rotation=xstr_slab_newbasis.c2f*xstr_bulk.f2c;  //VERY POWERFUL
    //xmatrix<double> rotation=trasp(trasp(xstr_slab_newbasis.lattice)*inverse(trasp(xstr_bulk.lattice)));
    //xmatrix<double> rotation=aurostd::inverse(trasp(xstr_slab_newbasis.lattice*aurostd::inverse(xstr_bulk.lattice)));
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " xstr_bulk.lattice=" << endl;cerr << xstr_bulk.lattice << endl;
      cerr << __AFLOW_FUNC__ << " xstr_slab_newbasis.lattice=" << endl;cerr << xstr_slab_newbasis.lattice << endl;
      cerr << __AFLOW_FUNC__ << " rotation=" << endl;cerr << rotation << endl;
    }

    if(total_layers%2!=0){message << "total_layers is odd, it is better to pick an even number (top vs. bottom)";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_WARNING_);}

    //rotate n_s too
    //[CO20190408 - do rotation last!]bool rotate_shear=false;
    //[CO20190408 - do rotation last!]if(rotate_shear){
    //[CO20190408 - do rotation last!]  n_s=R*n_s;  //CO20190322 - don't think this needs to be rotated
    //[CO20190408 - do rotation last!]  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " rotated n_s=" << n_s << endl;}  //CO20190322 - don't think this needs to be rotated
    //[CO20190408 - do rotation last!]}
    //[CO20190408 - do rotation last!]//create rotated primitive cell: faster to rotate smaller cell than bigger one
    //[CO20190408 - do rotation last!]xstructure xstr_slab_newbasis(a);
    //[CO20190408 - do rotation last!]if(0){  //try rotating atoms, DOES NOT WORK
    //[CO20190408 - do rotation last!]  _sym_op rot; rot.setUc(R,xstr_bulk.lattice); rot.is_pgroup=true;  //just rotation
    //[CO20190408 - do rotation last!]  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " sym_op=" << endl;cerr << rot << endl;}
    //[CO20190408 - do rotation last!]  if(1){  //try ApplyXstructure(), same as ApplyAtomValidate on all atoms individually
    //[CO20190408 - do rotation last!]    xstr_slab_newbasis=SYM::ApplyXstructure(rot,xstr_bulk,false);  //no incell, we'll figure this out later
    //[CO20190408 - do rotation last!]  } else {  //try ApplyAtomValidate
    //[CO20190408 - do rotation last!]    for(uint i=0;i<xstr_bulk.atoms.size();i++){
    //[CO20190408 - do rotation last!]      if(!SYM::ApplyAtomValidate(xstr_bulk.atoms[i],xstr_slab_newbasis.atoms[i],rot,xstr_bulk,true,false)){
    //[CO20190408 - do rotation last!]        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"ApplyAtomValidate() failed",_VALUE_ERROR_);
    //[CO20190408 - do rotation last!]      }
    //[CO20190408 - do rotation last!]    }
    //[CO20190408 - do rotation last!]  }
    //[CO20190408 - do rotation last!]} else {  //try rotating lattice, seems to work
    //[CO20190408 - do rotation last!]  //xstr_slab_newbasis.lattice=R * xstr_slab_newbasis.lattice * inverse(R);
    //[CO20190408 - do rotation last!]  xstr_slab_newbasis.lattice=R * xstr_slab_newbasis.lattice;
    //[CO20190408 - do rotation last!]  xstr_slab_newbasis.FixLattices();
    //[CO20190408 - do rotation last!]  const xmatrix<double>& f2c=xstr_slab_newbasis.f2c;
    //[CO20190408 - do rotation last!]  for(uint i=0;i<xstr_slab_newbasis.atoms.size();i++){xstr_slab_newbasis.atoms[i].cpos=f2c*xstr_slab_newbasis.atoms[i].fpos;}
    //[CO20190408 - do rotation last!]  //xstr_slab_newbasis=Standard_Conventional_UnitCellForm(a);xstr_slab_newbasis.clean(); //DX20191220 - uppercase to lowercase clean
    //[CO20190408 - do rotation last!]}
    //[CO20190408 - do rotation last!]if(LDEBUG) {cerr << __AFLOW_FUNC__ << " xstr_slab_newbasis=" << endl;cerr << xstr_slab_newbasis << endl;}
    //[CO20190408 - do rotation last!]if(check_min_dist){ //sanity check as we rotate structure/atoms
    //[CO20190408 - do rotation last!]  min_dist=xstr_slab_newbasis.MinDist();
    //[CO20190408 - do rotation last!]  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " mindist[" << count_check_min_dist++ << "]=" << min_dist << endl;}
    //[CO20190408 - do rotation last!]  if(!aurostd::isequal(min_dist_orig,min_dist)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Minimum distance changed",_VALUE_ERROR_);}
    //[CO20190408 - do rotation last!]}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - create slab
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - defining hkl normals
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " defining HKL normals" << endl;}

    //we need UN-ROTATED lattice here so we can get the right distance
    xvector<double> n_i=slab::HKLPlane2Normal(xstr_bulk.lattice,hkl_i);
    xvector<double> n_i_ORIG=n_i;
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " n_i[hkl=" << hkl_i << "]=" << n_i << endl;}

    //[CO20190515 - not necessarily true that n_i and n_s are perpendicular, n_s must have no 0 component AFTER rotation]//n_i and n_s must be perpendicular
    //[CO20190515 - not necessarily true that n_i and n_s are perpendicular, n_s must have no 0 component AFTER rotation]if(!aurostd::isequal(aurostd::scalar_product(n_i,n_s),0.0,0.1)){
    //[CO20190515 - not necessarily true that n_i and n_s are perpendicular, n_s must have no 0 component AFTER rotation]  message << "n_i(" << n_i << ") is not perpendicular to n_s(" << n_s << ")" << endl;
    //[CO20190515 - not necessarily true that n_i and n_s are perpendicular, n_s must have no 0 component AFTER rotation]  message << "n_i DOT n_s = " << aurostd::scalar_product(n_i,n_s) << endl;
    //[CO20190515 - not necessarily true that n_i and n_s are perpendicular, n_s must have no 0 component AFTER rotation]  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ERROR_);
    //[CO20190515 - not necessarily true that n_i and n_s are perpendicular, n_s must have no 0 component AFTER rotation]}

    //[CO20190515 - WRONG, [hkl] WITH brackets is ALREADY in direct space]xvector<double> n_s=HKLPlane2Normal(xstr_bulk.lattice,hkl_s);  //we need UN-ROTATED lattice here so we can get the right distance
    //[CO20190515 - WRONG, [hkl] WITH brackets is ALREADY in direct space]double d_layers=getDistanceBetweenImages(xstr_bulk.lattice,h_s,k_s,l_s); //this depends on UN-ROTATED lattice
    xvector<double> n_s=xstr_bulk.f2c*aurostd::xvector2utype<int,double>(hkl_s);n_s/=aurostd::modulus(n_s); //f2c=trasp(xstr_bulk.lattice)
    xvector<double> n_s_ORIG=n_s;
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " n_s[hkl=" << hkl_s << "]=" << n_s << endl;}

    double angle_ns_unrotated=aurostd::angle(n_i,n_s);
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " angle_ns(unrotated)=" << rad2deg*angle_ns_unrotated << endl;}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - defining hkl normals
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - get n_i_rotated and n_s_rotated
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //rotate n_i and n_s
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " n_i[hkl=" << hkl_i << "](unrotated)=" << n_i << endl;}
    xvector<double> n_i_rotated=rotation*n_i;n_i_rotated/=aurostd::modulus(n_i_rotated);  //rotate to match structure rotation
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " n_i[hkl=" << hkl_i << "](rotated)  =" << n_i_rotated << endl;}

    if(!aurostd::isequal(abs(n_i_rotated[3]),1.0)){
      message << "Rotation unsuccessful, n_i should be aligned along the z-axis" << endl;
      message << "n_i(rotated)  =" << n_i_rotated << endl;
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
    }

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " n_s[hkl=" << hkl_s << "](unrotated)=" << n_s << endl;}
    xvector<double> n_s_rotated=rotation*n_s;n_s_rotated/=aurostd::modulus(n_s_rotated);  //rotate to match structure rotation
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " n_s[hkl=" << hkl_s << "](rotated)  =" << n_s_rotated << endl;}

    double angle_ns_rotated=aurostd::angle(n_i_rotated,n_s_rotated);
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " angle_ns(unrotated)=" << rad2deg*angle_ns_unrotated << endl;
      cerr << __AFLOW_FUNC__ << " angle_ns(rotated)  =" << rad2deg*angle_ns_rotated << endl;
    }

    if(!aurostd::isequal(angle_ns_unrotated,angle_ns_rotated)){
      message << "Rotation unsuccessful" << endl;
      message << "angle_ns(unrotated)=" << rad2deg*angle_ns_unrotated << endl;
      message << "angle_ns(rotated)  =" << rad2deg*angle_ns_rotated << endl;
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - get n_i_rotated and n_s_rotated
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - find symmetrically equivalent n_s if necessary
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(!partial_dissociation){
      if(!aurostd::isequal(rad2deg*angle_ns_unrotated,90.0) || aurostd::isequal(n_s_rotated[3],0.0)){ //not sure about partial dissociation
        message << "Need to find symmetrically equivalent n_s (angle_ns_unrotated(degrees)=" << rad2deg*angle_ns_unrotated << ",n_s(rotated)=" << n_s_rotated << ")";
        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_WARNING_);

        const vector<_sym_op>& v_pgroups=xstr_sym.pgroup_xtal;
        vector<xvector<double> > v_n_s_sym,v_n_s_sym_rotated;
        vector<uint> v_pgs;
        xvector<double> n_s_sym,n_s_sym_rotated;
        double angle_nssym_unrotated=0.0;
        for(uint pg=0;pg<v_pgroups.size();pg++){
          n_s_sym=v_pgroups[pg].Uc*n_s;
          angle_nssym_unrotated=aurostd::angle(n_i,n_s_sym);
          if(aurostd::isequal(rad2deg*angle_nssym_unrotated,90.0)){
            n_s_sym_rotated=rotation*n_s_sym;
            if(aurostd::isequal(n_s_sym_rotated[3],0.0)){
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " found viable symmetrically equivalent n_s=" << n_s_sym << ", n_s_sym_rotated=" << n_s_sym_rotated << ", pg=" << pg << endl;}
              v_n_s_sym.push_back(n_s_sym);
              v_n_s_sym_rotated.push_back(n_s_sym_rotated);
              v_pgs.push_back(pg);
            }
          }
        }

        if(v_n_s_sym.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot find any viable symmetrically equivalent n_s",_RUNTIME_ERROR_);}

        //sort the two lists (symmetrically equivalent and rotated variants) by the rotated variants
        //look for 1 0 0 first, then -1 0 0, then
        //not absolutely necessary, but it will prevent us from rerunning symmetrically equivalent in the future (directory name)
        //it's also cleaner
        xvector<double> xvtmp,xvtmp_rotated;
        uint uitmp=0;
        bool swap=false;
        bool sing_dir_vec_i=false;      //if direction i is single direction vector (prefer these)
        bool sing_dir_vec_j=false;      //if direction j is single direction vector (prefer these)
        bool equ_comp_vec_i=false;      //if component 1 and 2 are equivalent (component 3 is 0) (prefer these)
        bool equ_comp_vec_j=false;      //if component 1 and 2 are equivalent (component 3 is 0) (prefer these)
        bool components_equ=false;      //if components are equal
        bool components_equ_abs=false;  //if abs() of components are equal
        bool j_gt_i_abs=false;          //if abs() of j k-component is greater than abs() of i k-component
        bool j_gt_i=false;              //if j k-component is greater than i k-component
        int lrows=n_s.lrows;
        for(uint i=0;i<v_n_s_sym_rotated.size()-1;i++){
          for(uint j=i+1;j<v_n_s_sym_rotated.size();j++){
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " comparing " << v_n_s_sym_rotated[i] << " with " << v_n_s_sym_rotated[j] << endl;}

            swap=false;

            //first look for single direction vectors (1.0000e+00  0.0000e+00  0.0000e+00 vs 5.0000e-01 -8.6603e-01  0.0000e+00), we prefer these
            //second look for vectors with equal components (7.0711e-01 -7.0711e-01  0.0000e+00)

            if(swap==false){
              sing_dir_vec_i=false;
              sing_dir_vec_j=false;
              for(int k=v_n_s_sym_rotated[i].lrows;k<=v_n_s_sym_rotated[i].urows;k++){
                if(aurostd::isequal(abs(v_n_s_sym_rotated[i][k]),1.0)){
                  if(LDEBUG){cerr << __AFLOW_FUNC__ << " [i] is single direction" << endl;}
                  sing_dir_vec_i=true;
                }
                if(aurostd::isequal(abs(v_n_s_sym_rotated[j][k]),1.0)){
                  if(LDEBUG){cerr << __AFLOW_FUNC__ << " [j] is single direction" << endl;}
                  sing_dir_vec_j=true;
                }
              }
              if(sing_dir_vec_i != sing_dir_vec_j){
                if(sing_dir_vec_j==true && sing_dir_vec_i==false){swap=true;}
                else{continue;}
              }
            }

            if(swap==false){
              equ_comp_vec_i=false;
              equ_comp_vec_j=false;
              if(aurostd::isequal(abs(v_n_s_sym_rotated[i][lrows]),abs(v_n_s_sym_rotated[i][lrows+1]))){equ_comp_vec_i=true;}
              if(aurostd::isequal(abs(v_n_s_sym_rotated[j][lrows]),abs(v_n_s_sym_rotated[j][lrows+1]))){equ_comp_vec_j=true;}
              if(equ_comp_vec_i != equ_comp_vec_j){
                if(equ_comp_vec_j==true && equ_comp_vec_i==false){swap=true;}
                //otherwise compare by elements
              }
            }

            //finally compare individual components
            if(swap==false){
              for(int k=v_n_s_sym_rotated[i].lrows;k<=v_n_s_sym_rotated[i].urows && swap==false;k++){
                components_equ=aurostd::isequal(v_n_s_sym_rotated[i][k],v_n_s_sym_rotated[j][k]);
                if(LDEBUG){cerr << __AFLOW_FUNC__ << " [i][k=" << k << "] ?== [j][k=" << k << "] == " << components_equ << endl;}
                if(components_equ==false){
                  components_equ_abs=aurostd::isequal(abs(v_n_s_sym_rotated[i][k]),abs(v_n_s_sym_rotated[j][k]));
                  j_gt_i_abs=bool(abs(v_n_s_sym_rotated[j][k])>abs(v_n_s_sym_rotated[i][k]));
                  j_gt_i=bool(v_n_s_sym_rotated[j][k]>v_n_s_sym_rotated[i][k]);
                  if(LDEBUG){
                    cerr << __AFLOW_FUNC__ << " abs([i][k=" << k << "]) ?== abs([j][k=" << k << "]) == " << components_equ_abs << endl;
                    cerr << __AFLOW_FUNC__ << " abs([j][k=" << k << "]) ?> abs([i][k=" << k << "]) == " << j_gt_i_abs << endl;
                    cerr << __AFLOW_FUNC__ << " [j][k=" << k << "] ?> [i][k=" << k << "] == " << j_gt_i << endl;
                  }
                  if( (components_equ_abs==false && j_gt_i_abs==true) || 
                      (components_equ_abs==true && components_equ==false && j_gt_i==true) ){swap=true;}
                  break;  //very important, comparing 5 with -5
                }
              }
            }
            if(swap){
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " swapping " << v_n_s_sym_rotated[i] << " and " << v_n_s_sym_rotated[j] << endl;}
              xvtmp=v_n_s_sym[i];
              xvtmp_rotated=v_n_s_sym_rotated[i];
              uitmp=v_pgs[i];
              v_n_s_sym[i]=v_n_s_sym[j];
              v_n_s_sym_rotated[i]=v_n_s_sym_rotated[j];
              v_pgs[i]=v_pgs[j];
              v_n_s_sym[j]=xvtmp;
              v_n_s_sym_rotated[j]=xvtmp_rotated;
              v_pgs[j]=uitmp;
            }
          }
        }

        if(LDEBUG){
          for(uint i=0;i<v_n_s_sym.size();i++){
            cerr << __AFLOW_FUNC__ << " v_n_s_sym[i=" << i << "]=" << v_n_s_sym[i] << ", v_n_s_sym_rotated[i=" << i << "]=" << v_n_s_sym_rotated[i] << ", v_pgs[i=" << i << "]="<< v_pgs[i] << endl;
          }
        }

        n_s=v_n_s_sym[0];
        n_s_rotated=v_n_s_sym_rotated[0];


        //fix hkl_s as well for directory labeling
        //cannot use n_s, as it is normalized
        //rotated hkl_s by v_pgroups
        //rotate in Cartesian coordinates, then convert to fractional
        xvector<double> n_s_tmp=xstr_bulk.f2c*aurostd::xvector2utype<int,double>(hkl_s);  //do not normalize
        n_s_tmp=v_pgroups[v_pgs[0]].Uc*n_s_tmp; //rotate
        hkl_s=aurostd::xvector2utype<double,int>(xstr_bulk.c2f*n_s_tmp);  //convert to fractional

        message << "Selecting symmetrically equivalent hkl_s=" << hkl_s << endl;
        message << "n_s=" << n_s << endl;
        message << "n_s(rotated)=" << n_s_rotated << endl;
        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
      }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - find symmetrically equivalent n_s if necessary
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - resolve layers count
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //d_spacing is based on lattice ONLY, no crystal information
    //d_layers,d_cells are based on CRYSTAL (lattice + basis)
    //d_layers finds distance in direction of n_s to next structure image
    //d_cells makes sure you loop outside cell at least once
    //d_spacing <= d_layers <= d_cells
    double d_spacing_i=slab::getSpacingHKLPlane(xstr_bulk,hkl_i); //aurostd::modulus(xstr_slab.lattice(1))/sqrt(h_i*h_i+k_i*k_i+l_i*l_i);
    double d_layers_i=slab::getDistanceBetweenImages(xstr_bulk,n_i,false); //this depends on UN-ROTATED lattice
    double d_cells_i=slab::getDistanceBetweenImages(xstr_bulk,n_i,true); //go outside cell
    int layers_per_cell_i=(int)(d_cells_i/d_layers_i);  //floor
    int supercell_layers_i=(total_layers+layers_per_cell_i-1)/layers_per_cell_i;  //ceil //(double)total_layers;
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " n_i[hkl=" << hkl_i << "]=" << n_i << endl;
      cerr << __AFLOW_FUNC__ << " d_spacing_i=" << d_spacing_i << endl;
      cerr << __AFLOW_FUNC__ << " d_layers_i=" << d_layers_i << endl;
      cerr << __AFLOW_FUNC__ << " d_cells_i=" << d_cells_i << endl;
      cerr << __AFLOW_FUNC__ << " layers_per_cell_i=" << layers_per_cell_i << endl;
      cerr << __AFLOW_FUNC__ << " supercell_layers_i=" << supercell_layers_i << endl;
    }
    double d_spacing_s=slab::getSpacingHKLPlane(xstr_bulk,hkl_s); //aurostd::modulus(xstr_slab.lattice(1))/sqrt(h_s*h_s+k_s*k_s+l_s*l_s);
    double d_layers_s=slab::getDistanceBetweenImages(xstr_bulk,n_s,false); //this depends on UN-ROTATED lattice
    double d_cells_s=slab::getDistanceBetweenImages(xstr_bulk,n_s,true); //go outside cell
    int layers_per_cell_s=(int)(d_cells_s/d_layers_s);  //floor
    int supercell_layers_s=(total_layers+layers_per_cell_s-1)/layers_per_cell_s;  //ceil //(double)total_layers;
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " n_s[hkl=" << hkl_s << "]=" << n_s << endl;
      cerr << __AFLOW_FUNC__ << " d_spacing_s=" << d_spacing_s << endl;
      cerr << __AFLOW_FUNC__ << " d_layers_s=" << d_layers_s << endl;
      cerr << __AFLOW_FUNC__ << " d_cells_s=" << d_cells_s << endl;
      cerr << __AFLOW_FUNC__ << " layers_per_cell_s=" << layers_per_cell_s << endl;
      cerr << __AFLOW_FUNC__ << " supercell_layers_s=" << supercell_layers_s << endl;
    }
    //[CO20190520 - do NOT reduce total_layers based on shear direction]if(layers_per_cell_s<1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"layers_per_cell_s<1 (="+aurostd::utype2string(layers_per_cell_s)+")",_INPUT_ERROR_);}
    //[CO20190520 - do NOT reduce total_layers based on shear direction]if(layers_per_cell_s>1){  //let's adjust total_layers input to CreateSlab_RigidRotation() which builds supercell assuming layers_per_cell_s==1
    //[CO20190520 - do NOT reduce total_layers based on shear direction]  int total_layers=DEFAULT_TOTAL_LAYERS;
    //[CO20190520 - do NOT reduce total_layers based on shear direction]  string total_layers_string=vpflow.getattachedscheme("CREATE_SLAB::TOTAL_LAYERS");
    //[CO20190520 - do NOT reduce total_layers based on shear direction]  if(aurostd::isfloat(total_layers_string)){
    //[CO20190520 - do NOT reduce total_layers based on shear direction]    int _total_layers=aurostd::string2utype<int>(total_layers_string);
    //[CO20190520 - do NOT reduce total_layers based on shear direction]    if(_total_layers>0){total_layers=_total_layers;}
    //[CO20190520 - do NOT reduce total_layers based on shear direction]  }
    //[CO20190520 - do NOT reduce total_layers based on shear direction]  if(LDEBUG){cerr << __AFLOW_FUNC__ << " total_layers(pre)=" << total_layers << endl;}
    //[CO20190520 - do NOT reduce total_layers based on shear direction]  total_layers=(total_layers+layers_per_cell_s-1)/layers_per_cell_s;  //ceil
    //[CO20190520 - do NOT reduce total_layers based on shear direction]  if(LDEBUG){cerr << __AFLOW_FUNC__ << " total_layers(post)=" << total_layers << endl;}
    //[CO20190520 - do NOT reduce total_layers based on shear direction]  slab_flags.pop_attached("CREATE_SLAB::TOTAL_LAYERS");
    //[CO20190520 - do NOT reduce total_layers based on shear direction]  slab_flags.push_attached("CREATE_SLAB::TOTAL_LAYERS",aurostd::utype2string(total_layers));
    //[CO20190520 - do NOT reduce total_layers based on shear direction]}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - resolve layers count
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - rotate n_i and n_s
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //swap for rotated variants, we saved originals already
    n_i=n_i_rotated;
    n_s=n_s_rotated;

    if(!partial_dissociation){
      //check that there is NO z component to n_s after rotation
      if(!aurostd::isequal(n_s[3],0.0)){
        message << "The shear plane and plane of interest are incommensurate" << endl;
        message << "rotation=" << endl; message << rotation << endl;
        message << "hkl_i=" << hkl_i << ", n_i(rotated)=" << n_i << endl;
        message << "hkl_s=" << hkl_s << ", n_s(rotated)=" << n_s << endl;
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
      }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - rotate n_i and n_s
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - define half plane
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //[CO20190423 - too much work, rely on atom.ijk instead!]if(0){  //too much work, rely on atom.ijk instead!
    //[CO20190423 - too much work, rely on atom.ijk instead!]if(LDEBUG) {cerr << __AFLOW_FUNC__ << " defining half plane" << endl;}
    //[CO20190423 - too much work, rely on atom.ijk instead!]
    //[CO20190423 - too much work, rely on atom.ijk instead!]int count_total,count_above,count_below;
    //[CO20190423 - too much work, rely on atom.ijk instead!]count_total=xstr_slab.atoms.size();
    //[CO20190423 - too much work, rely on atom.ijk instead!]count_above=count_total/2;count_below=count_total-count_above;
    //[CO20190423 - too much work, rely on atom.ijk instead!]//count_above=count_below=count_total/2;
    //[CO20190423 - too much work, rely on atom.ijk instead!]//int count_above=0,count_below=0,count_total=xstr_slab.atoms.size();
    //[CO20190423 - too much work, rely on atom.ijk instead!]
    //[CO20190423 - too much work, rely on atom.ijk instead!]if(LDEBUG) {
    //[CO20190423 - too much work, rely on atom.ijk instead!]  cerr << __AFLOW_FUNC__ << " count_above=" << count_above << endl;
    //[CO20190423 - too much work, rely on atom.ijk instead!]  cerr << __AFLOW_FUNC__ << " count_below=" << count_below << endl;
    //[CO20190423 - too much work, rely on atom.ijk instead!]  cerr << __AFLOW_FUNC__ << " count_total=" << count_total << endl;
    //[CO20190423 - too much work, rely on atom.ijk instead!]}
    //[CO20190423 - too much work, rely on atom.ijk instead!]if(count_above!=count_below){
    //[CO20190423 - too much work, rely on atom.ijk instead!]  message << "count_above!=count_below";
    //[CO20190423 - too much work, rely on atom.ijk instead!]  pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_WARNING_);
    //[CO20190423 - too much work, rely on atom.ijk instead!]}
    //[CO20190423 - too much work, rely on atom.ijk instead!]if(count_above+count_below!=count_total){
    //[CO20190423 - too much work, rely on atom.ijk instead!]  message << "count_above+count_below!=count_total";
    //[CO20190423 - too much work, rely on atom.ijk instead!]  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ERROR_); //CO20190226
    //[CO20190423 - too much work, rely on atom.ijk instead!]}
    //[CO20190423 - too much work, rely on atom.ijk instead!]
    //[CO20190423 - too much work, rely on atom.ijk instead!]//half-plane is bad idea, too much room for error
    //[CO20190423 - too much work, rely on atom.ijk instead!]////define half plane along c-axis
    //[CO20190423 - too much work, rely on atom.ijk instead!]////normal vector is simply c-axis
    //[CO20190423 - too much work, rely on atom.ijk instead!]//xvector<double> n_hp=xstr_slab_newbasis_supercell.lattice(3);n_hp/=aurostd::modulus(n_hp); //unit vector
    //[CO20190423 - too much work, rely on atom.ijk instead!]//xvector<double> p_hp=aurostd::modulus(xstr_slab_newbasis_supercell.lattice(3))/2.0*n_hp; //half plane value
    //[CO20190423 - too much work, rely on atom.ijk instead!]//double D_hp=-aurostd::scalar_product(n_hp,p_hp);
    //[CO20190423 - too much work, rely on atom.ijk instead!]//if(LDEBUG) {
    //[CO20190423 - too much work, rely on atom.ijk instead!]//  cerr << __AFLOW_FUNC__ << " n_hp=" << n_hp << endl;
    //[CO20190423 - too much work, rely on atom.ijk instead!]//  cerr << __AFLOW_FUNC__ << " p_hp=" << p_hp << endl;
    //[CO20190423 - too much work, rely on atom.ijk instead!]//  cerr << __AFLOW_FUNC__ << " D_hp=" << D_hp << endl;
    //[CO20190423 - too much work, rely on atom.ijk instead!]//}

    //[CO20190423 - too much work, rely on atom.ijk instead!]//define zero plane along c-axis
    //[CO20190423 - too much work, rely on atom.ijk instead!]//normal vector is simply c-axis
    //[CO20190423 - too much work, rely on atom.ijk instead!]xvector<double> n_zp=xstr_slab.lattice(3);n_zp/=aurostd::modulus(n_zp); //unit vector
    //[CO20190423 - too much work, rely on atom.ijk instead!]xvector<double> p_zp; //0,0,0
    //[CO20190423 - too much work, rely on atom.ijk instead!]double D_zp=-aurostd::scalar_product(n_zp,p_zp);
    //[CO20190423 - too much work, rely on atom.ijk instead!]if(LDEBUG) {
    //[CO20190423 - too much work, rely on atom.ijk instead!]  cerr << __AFLOW_FUNC__ << " n_zp=" << n_zp << endl;
    //[CO20190423 - too much work, rely on atom.ijk instead!]  cerr << __AFLOW_FUNC__ << " p_zp=" << p_zp << endl;
    //[CO20190423 - too much work, rely on atom.ijk instead!]  cerr << __AFLOW_FUNC__ << " D_zp=" << D_zp << endl;
    //[CO20190423 - too much work, rely on atom.ijk instead!]}
    //[CO20190423 - too much work, rely on atom.ijk instead!]}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - define half plane
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - identify selective dynamics
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " identifying selective dynamics" << endl;}

    //k in atom.ijk will go from 0 to supercell_layers_i-1

    //sometimes k_min==0, sometimes k_min==1 (look at GetSuperCell())
    //get k_min and k_max first
    int k=0,k_min=xstr_slab.atoms.size(),k_max=-xstr_slab.atoms.size();
    for(uint i=0;i<xstr_slab.atoms.size();i++){
      const _atom& atom=xstr_slab.atoms[i];
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " atom.fpos=" << atom.fpos << ", atom.ijk=" << atom.ijk << endl;}
      k=atom.ijk[3];
      if(k<k_min){k_min=k;}
      if(k>k_max){k_max=k;}
    }
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " k_min=" << k_min << ", k_max=" << k_max << endl;}
    if(k_max-k_min+1!=supercell_layers_i){  //test of stupidity
      message << "k_max-k_min+1!=supercell_layers_i";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ERROR_);
    }

    uint count_total_fixed=0,count_bottom_fixed=0,count_top_fixed=0; //tests of stupidity
    bool fixed_bottom=false,fixed_top=false;

    xstr_slab.isd=true; //set selective dynamics
    for(uint i=0;i<xstr_slab.atoms.size();i++){
      _atom& atom=xstr_slab.atoms[i];
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " atom.fpos=" << atom.fpos << ", atom.ijk=" << atom.ijk << endl;}
      k=atom.ijk[3];
      fixed_bottom=(k<(k_min+fixed_layers));
      fixed_top=(k>(k_min+(supercell_layers_i-fixed_layers-1)));
      //keeping top/bottom fixed provides shielding to the effect of vacuum
      if(fixed_bottom || fixed_top){  //keep fixed
        atom.sd="FFF";  //keep fixed
        count_total_fixed++;
        if(fixed_bottom){count_bottom_fixed++;}
        else if(fixed_top){count_top_fixed++;}
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " atom.ijk=" << atom.ijk << " is FIXED" << endl;}
      } else {
        //we only allow relaxation in z-direction because we want to allow planes of atoms to come together/go apart
        //atom.sd="TTT";  //default, allow relaxation
        atom.sd="FFT";  //default, allow relaxation //CORRECTION: only relax in z direction: Vitek Phil Mag 18, 773-786 (1968), also http://theory.cm.utexas.edu/forum/viewtopic.php?t=3301
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " atom.ijk=" << atom.ijk << " will RELAX" << endl;}
      }
    }
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " count_bottom_fixed=" << count_bottom_fixed << endl;
      cerr << __AFLOW_FUNC__ << " count_top_fixed=" << count_top_fixed << endl;
      cerr << __AFLOW_FUNC__ << " count_bottom_fixed+count_top_fixed=" << count_bottom_fixed+count_top_fixed << endl;
      cerr << __AFLOW_FUNC__ << " count_total_fixed=" << count_total_fixed << endl;
    }
    if(count_total_fixed!=count_bottom_fixed+count_top_fixed){
      message << "count_total_fixed!=count_bottom_fixed+count_top_fixed (" << count_total_fixed << "!=" << count_bottom_fixed+count_top_fixed << ")";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ERROR_); //CO20190226
    }

    //[CO20190423 - too much work, rely on atom.ijk instead!]if(0){  //too much work, rely on atom.ijk instead!
    //[CO20190423 - too much work, rely on atom.ijk instead!]//http://mathworld.wolfram.com/Plane.html
    //[CO20190423 - too much work, rely on atom.ijk instead!]xstr_slab.isd=true; //set selective dynamics
    //[CO20190423 - too much work, rely on atom.ijk instead!]double signed_point_plane_distance=0.0;
    //[CO20190423 - too much work, rely on atom.ijk instead!]vector<atom_plane_dist> v_apd;
    //[CO20190423 - too much work, rely on atom.ijk instead!]for(uint i=0;i<xstr_slab.atoms.size();i++){
    //[CO20190423 - too much work, rely on atom.ijk instead!]  xstr_slab.atoms[i].sd="TTT";  //default, change to FFF later
    //[CO20190423 - too much work, rely on atom.ijk instead!]  signed_point_plane_distance=aurostd::scalar_product(n_zp,xstr_slab.atoms[i].cpos)+D_zp; //n_zp is already normalized to 1.0
    //[CO20190423 - too much work, rely on atom.ijk instead!]  if(LDEBUG) {
    //[CO20190423 - too much work, rely on atom.ijk instead!]    cerr << __AFLOW_FUNC__ << " atom.ijk=" << xstr_slab.atoms[i].ijk << endl;
    //[CO20190423 - too much work, rely on atom.ijk instead!]    cerr << __AFLOW_FUNC__ << " signed_point_plane_distance[atom=" << i << "]=" << signed_point_plane_distance << endl;
    //[CO20190423 - too much work, rely on atom.ijk instead!]  }
    //[CO20190423 - too much work, rely on atom.ijk instead!]  //(std::signbit(signed_point_plane_distance + _ZERO_TOL_ ) ? count_below++ : count_above++);  //not great about atoms right on plane
    //[CO20190423 - too much work, rely on atom.ijk instead!]  //[do later]((signed_point_plane_distance<-_ZERO_TOL_)? count_below++ : count_above++);
    //[CO20190423 - too much work, rely on atom.ijk instead!]  v_apd.push_back(atom_plane_dist());
    //[CO20190423 - too much work, rely on atom.ijk instead!]  v_apd.back().index=i;
    //[CO20190423 - too much work, rely on atom.ijk instead!]  v_apd.back().distance=signed_point_plane_distance;
    //[CO20190423 - too much work, rely on atom.ijk instead!]}
    //[CO20190423 - too much work, rely on atom.ijk instead!]std::sort(v_apd.begin(),v_apd.end()); //sort by distance, we will take from bottom/top halves
    //[CO20190423 - too much work, rely on atom.ijk instead!]
    //[CO20190423 - too much work, rely on atom.ijk instead!]//get number to keep fixed, assume a layer is one unit cell
    //[CO20190423 - too much work, rely on atom.ijk instead!]//uint count_keep_fixed=xstr_bulk.atoms.size() * supercell_layers_i * supercell_layers_i * fixed_layers;
    //[CO20190423 - too much work, rely on atom.ijk instead!]uint count_keep_fixed=xstr_bulk.atoms.size() * xy_dims * xy_dims * fixed_layers;
    //[CO20190423 - too much work, rely on atom.ijk instead!]if(LDEBUG) {cerr << __AFLOW_FUNC__ << " count_keep_fixed=" << count_keep_fixed << endl;}
    //[CO20190423 - too much work, rely on atom.ijk instead!]if(2*count_keep_fixed>xstr_slab.atoms.size()){
    //[CO20190423 - too much work, rely on atom.ijk instead!]  message << "2*count_keep_fixed>xstr_slab.atoms.size()";
    //[CO20190423 - too much work, rely on atom.ijk instead!]  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ERROR_); //CO20190226
    //[CO20190423 - too much work, rely on atom.ijk instead!]}
    //[CO20190423 - too much work, rely on atom.ijk instead!]
    //[CO20190423 - too much work, rely on atom.ijk instead!]//fixed selective dynamics stuff at once
    //[CO20190423 - too much work, rely on atom.ijk instead!]uint count_check=0;
    //[CO20190423 - too much work, rely on atom.ijk instead!]for(uint i=0;i<count_keep_fixed;i++){xstr_slab.atoms[v_apd[i].index].sd="FFF";count_check++;} //bottom
    //[CO20190423 - too much work, rely on atom.ijk instead!]if(count_check!=count_keep_fixed){
    //[CO20190423 - too much work, rely on atom.ijk instead!]  message << "count_check!=count_keep_fixed [1]";
    //[CO20190423 - too much work, rely on atom.ijk instead!]  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ERROR_); //CO20190226
    //[CO20190423 - too much work, rely on atom.ijk instead!]}
    //[CO20190423 - too much work, rely on atom.ijk instead!]count_check=0;
    //[CO20190423 - too much work, rely on atom.ijk instead!]for(uint i=v_apd.size()-1;i>v_apd.size()-1-count_keep_fixed;i--){xstr_slab.atoms[v_apd[i].index].sd="FFF";count_check++;} //top
    //[CO20190423 - too much work, rely on atom.ijk instead!]if(count_check!=count_keep_fixed){
    //[CO20190423 - too much work, rely on atom.ijk instead!]  message << "count_check!=count_keep_fixed [2]";
    //[CO20190423 - too much work, rely on atom.ijk instead!]  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ERROR_); //CO20190226
    //[CO20190423 - too much work, rely on atom.ijk instead!]}
    //[CO20190423 - too much work, rely on atom.ijk instead!]}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - identify selective dynamics
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - create shear sub-directories
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //I was going to create a check to compare fpos' of atoms from shear_fraction==0 and shear_fraction==1
    //however, they will NOT be necessarily the same
    //if we rotate along 111 and shear along 112, the atoms will go into a new layer (z-axis)
    //need to check this thoroughly

    int dir_count=0,dir_count_total=0;
    bool create_final_duplicate_directory=(true || LDEBUG); //debug really, we don't need the 1.0 shear fraction directory
    double shear_fraction=0.0;
    double shear_fraction_final=(1.0+(create_final_duplicate_directory?_ZERO_TOL_:-_ZERO_TOL_));  //it's just a LITTLE bigger than 1.0 for the while loop

    //for partial dissociation
    double z_shift_fpos=0.0;  //for partial dissociation
    xvector<double> n_z=xstr_slab_newbasis.lattice(3);
    if(partial_dissociation){
      xvector<double> total_shift_cpos=(d_layers_s) * n_s;
      xvector<double> total_shift_fpos=xstr_slab_newbasis.c2f*total_shift_cpos;
      z_shift_fpos=total_shift_fpos[3];
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " total_shift_cpos=" << total_shift_cpos << endl;
        cerr << __AFLOW_FUNC__ << " total_shift_fpos=" << total_shift_fpos << endl;
        cerr << __AFLOW_FUNC__ << " z_shift_fpos=" << z_shift_fpos << endl;
      }
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " z_shift_fpos=" << z_shift_fpos << endl;}

    //stringstream POSCAR;
    stringstream new_AflowIn_ss;
    stringstream arun_dirname;
    _xvasp xvasp;
    //_aflags aflags; aflags.Directory=".";
    //_kflags kflags; 
    //_vflags vflags;
    //[CO20190405 - does not work, could be multiple layers within this spacing]double d_layers_s=getSpacingHKLPlane(xstr_slab.lattice,h_s,k_s,l_s); //aurostd::modulus(xstr_slab.lattice(1))/sqrt(h_s*h_s+k_s*k_s+l_s*l_s);
    //get dir_count_total
    dir_count=0;shear_fraction=0.0;
    while(shear_fraction<shear_fraction_final){shear_fraction+=step_size;dir_count_total++;}
    dir_count=0;shear_fraction=0.0;

    //check if vasp is already done
    bool all_vasp_done=true;
    string FileName;
    while(shear_fraction<shear_fraction_final){
      xvasp.clear();
      //create directory name
      xvasp.Directory=aflags.Directory;
      xvasp.AVASP_arun=true;
      xvasp.AVASP_arun_mode="GSFE"; //generalized stacking fault energy
      arun_dirname.str("");
      arun_dirname << std::setfill('0') << std::setw(aurostd::getZeroPadding(dir_count_total)) << dir_count+1 << "_";
      //arun_dirname << "iHKL" << aurostd::utype2string(hkl_i[1]) << aurostd::utype2string(hkl_i[2]) << aurostd::utype2string(hkl_i[3]) << "-";
      //arun_dirname << "sHKL" << aurostd::utype2string(hkl_s[1]) << aurostd::utype2string(hkl_s[2]) << aurostd::utype2string(hkl_s[3]) << "-";
      //arun_dirname << "sf" << aurostd::utype2string(shear_fraction,8);
      arun_dirname << "PL=" << aurostd::joinWDelimiter(hkl_i,":") << "-"; //cannot have - before numbers as there can be negative directions, use : to help readability
      arun_dirname << "DIR=" << aurostd::joinWDelimiter(hkl_s,":") << "-";//cannot have - before numbers as there can be negative directions, use : to help readability
      arun_dirname << "FRAC=" << aurostd::utype2string(shear_fraction,3);  //fix precision here if you want
      xvasp.AVASP_arun_runname=arun_dirname.str();
      AVASP_populateXVASP(aflags,kflags,vflags,xvasp);

      FileName=xvasp.Directory + "/" + DEFAULT_AFLOW_QMVASP_OUT;
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking for: " << FileName << endl;}
      if(!aurostd::EFileExist(FileName)){all_vasp_done=false;break;}
      shear_fraction+=step_size;
      dir_count++;
    }
    dir_count=0;shear_fraction=0.0;

    if(!all_vasp_done){message << "Creating sheared VASP runs";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);}

    //create directories if vasp is not done
    int half_k=supercell_layers_i/2;  //floor
    //half_k=(supercell_layers_i+2-1)/2;  //ceil
    uint count_total=0,count_bottom=0,count_top=0; //tests of stupidity
    while(shear_fraction<shear_fraction_final){
      if(all_vasp_done==false){
        xvasp.clear();
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " shear_fraction=" << shear_fraction << endl;}

        //apply shear
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " applying shear" << endl;}
        xstructure xstr_shear(xstr_slab);
        xstr_shear.write_DEBUG_flag=FALSE;  //FORCE
        const xmatrix<double>& c2f=xstr_slab.c2f;

        count_total=0;count_bottom=0;count_top=0;
        for(uint i=0;i<xstr_shear.atoms.size();i++){
          _atom& atom=xstr_shear.atoms[i];
          if(LDEBUG) {cerr << __AFLOW_FUNC__ << " atom.fpos=" << atom.fpos << endl;}
          k=atom.ijk[3];
          if(k<(k_min+half_k)){ //not shearing
            if(LDEBUG) {cerr << __AFLOW_FUNC__ << " atom.ijk=" << atom.ijk << " NOT shearing" << endl;}
            count_bottom++;
          } else {  //shearing
            if(LDEBUG) {cerr << __AFLOW_FUNC__ << " atom.ijk=" << atom.ijk << " shearing" << endl;}
            atom.cpos += (shear_fraction * d_layers_s) * n_s;
            if(partial_dissociation){atom.cpos -= (shear_fraction * z_shift_fpos) * n_z;}
            atom.fpos = c2f*atom.cpos; //C2F(xstr_shear.lattice,atom.cpos);
            count_top++;
          }
          count_total++;
        }
        if(LDEBUG) {
          cerr << __AFLOW_FUNC__ << " count_bottom=" << count_bottom << endl;
          cerr << __AFLOW_FUNC__ << " count_top=" << count_top << endl;
          cerr << __AFLOW_FUNC__ << " count_bottom+count_top=" << count_bottom+count_top << endl;
          cerr << __AFLOW_FUNC__ << " count_total=" << count_total << endl;
          cerr << __AFLOW_FUNC__ << " xstr_shear.atoms.size()=" << xstr_shear.atoms.size() << endl;
        }
        if(count_total!=xstr_shear.atoms.size()){
          message << "count_total!=xstr_shear.atoms.size()";
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ERROR_);
        }
        if(count_total!=count_bottom+count_top){
          message << "count_total!=count_bottom+count_top";
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ERROR_);
        }

        //[CO20190423 - too much work, rely on atom.ijk instead!]if(0){  //too much work, rely on atom.ijk instead!
        //[CO20190423 - too much work, rely on atom.ijk instead!]  for(uint i=v_apd.size()-1;i>v_apd.size()-1-count_above;i--){  //top (furthest away)
        //[CO20190423 - too much work, rely on atom.ijk instead!]    _atom& atom=xstr_shear.atoms[v_apd[i].index];
        //[CO20190423 - too much work, rely on atom.ijk instead!]    atom.cpos += (shear_fraction * d_layers_s) * n_s;
        //[CO20190423 - too much work, rely on atom.ijk instead!]    atom.fpos = c2f*atom.cpos; //C2F(xstr_shear.lattice,atom.cpos);
        //[CO20190423 - too much work, rely on atom.ijk instead!]  }
        //[CO20190423 - too much work, rely on atom.ijk instead!]}
        //for(uint i=0;i<v_apd.size();i++)
        //  //if(!std::signbit(v_apd[i].distance)) //not great about atoms right on plane
        //  { //CO20200106 - patching for auto-indenting
        //  if(signed_point_plane_distance>=-_ZERO_TOL_){ //above plane, opposite of what is defined above for count_below/count_above
        //    xstr_shear.atoms[v_apd[i].index].cpos += shear_fraction * d_layers_s;
        //    xstr_shear.atoms[v_apd[i].index].fpos = C2F(xstr_shear.lattice,xstr_shear.atoms[v_apd[i].index].cpos);
        //  }
        //}
        xstr_shear.BringInCell(); //wrap around
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " xstr_shear=" << endl;cerr << xstr_shear << endl;}

        if(check_min_dist){ //sanity check as we rotate structure/atoms
          min_dist=xstr_shear.MinDist();
          if(LDEBUG) {cerr << __AFLOW_FUNC__ << " mindist[" << (count_check_min_dist++)+dir_count << "]=" << min_dist << endl;}
          if(!aurostd::isequal(min_dist_orig,min_dist)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Minimum distance changed",_VALUE_ERROR_);}
        }

        //if(0){
        //  string destination=aurostd::utype2string(h_i)+aurostd::utype2string(k_i)+aurostd::utype2string(l_i) + "/" + aurostd::utype2string(dir_count);
        //  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " creating " << destination << endl;}
        //  aurostd::DirectoryMake(destination);
        //  destination+="/POSCAR";
        //  stringstream POSCAR; POSCAR.str("");
        //  POSCAR << xstr;
        //  aurostd::stringstream2file(POSCAR,destination);
        //}

        //load in xstructure
        xvasp.str.clear(); //DX20191220 - uppercase to lowercase clear
        xvasp.str=xstr_shear;

        //create directory name
        xvasp.Directory=aflags.Directory;
        xvasp.AVASP_arun=true;
        xvasp.AVASP_arun_mode="GSFE"; //generalized stacking fault energy
        arun_dirname.str("");
        arun_dirname << std::setfill('0') << std::setw(aurostd::getZeroPadding(dir_count_total)) << dir_count+1 << "_";
        //arun_dirname << "iHKL" << aurostd::utype2string(hkl_i[1]) << aurostd::utype2string(hkl_i[2]) << aurostd::utype2string(hkl_i[3]) << "-";
        //arun_dirname << "sHKL" << aurostd::utype2string(hkl_s[1]) << aurostd::utype2string(hkl_s[2]) << aurostd::utype2string(hkl_s[3]) << "-";
        //arun_dirname << "sf" << aurostd::utype2string(shear_fraction,8);
        arun_dirname << "PL=" << aurostd::joinWDelimiter(hkl_i,":") << "-"; //cannot have - before numbers as there can be negative directions, use : to help readability
        arun_dirname << "DIR=" << aurostd::joinWDelimiter(hkl_s,":") << "-";//cannot have - before numbers as there can be negative directions, use : to help readability
        arun_dirname << "FRAC=" << aurostd::utype2string(shear_fraction,3);  //fix precision here if you want
        xvasp.AVASP_arun_runname=arun_dirname.str();
        AVASP_populateXVASP(aflags,kflags,vflags,xvasp);

        setPreserveUnitCell(xvasp);
        //set k-points to 11x11x1 as per 10.1103/PhysRevLett.114.165502, this is just a quick fix for now
        if(1){
          xvasp.aopts.flag("FLAG::KPOINTS_IMPLICIT",FALSE);
          //xvasp.aopts.flag("FLAG::KPOINTS_EXPLICIT",TRUE);
          xvasp.aopts.flag("FLAG::KPOINTS_EXPLICIT_START_STOP",TRUE);
          xvasp.AVASP_KPOINTS_EXPLICIT_START_STOP.str(""); //clear
          xvasp.AVASP_KPOINTS_EXPLICIT_START_STOP << "KPOINTS" << endl;
          xvasp.AVASP_KPOINTS_EXPLICIT_START_STOP << "0" << endl;
          xvasp.AVASP_KPOINTS_EXPLICIT_START_STOP << "Gamma" << endl;
          xvasp.AVASP_KPOINTS_EXPLICIT_START_STOP << "11 11 1" << endl;
          xvasp.AVASP_KPOINTS_EXPLICIT_START_STOP << "0 0 0" << endl;
        }

        //relax ions ONLY
        xvasp.aopts.flag("AFLOWIN_FLAG::RELAX_TYPE",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::RELAX_TYPE","IONS");
        xvasp.aopts.flag("FLAG::VOLUME_PRESERVED",TRUE);  //no volume changes +/*
        if(spin_off){xvasp.aopts.flag("FLAG::AVASP_SPIN",FALSE);}

        //make aflow.in
        AVASP_MakeSingleAFLOWIN(xvasp,new_AflowIn_ss,true); //false,-1,false);  //don't write/print and hence don't pthread
      }
      shear_fraction+=step_size;
      dir_count++;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - create shear sub-directories
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(!all_vasp_done){
      message << "Now waiting for sheared VASP runs";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_COMPLETE_);
      return;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - postprocessing
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - postprocessing
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  }

  //CO20190321
  //follows procedure outlined in: W. Sun and G. Ceder, Surface Science 617 (2013) 53-59
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow,istream& input,ostream& oss){
    xstructure xstr_in(input,IOAFLOW_AUTO);
    return CleavageEnergyCalculation(vpflow,xstr_in,oss);
  }
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow,const xstructure& xstr_in,ostream& oss){
    _aflags aflags; aflags.Directory=".";
    _kflags kflags;
    _vflags vflags;
    return CleavageEnergyCalculation(vpflow,xstr_in,aflags,kflags,vflags,oss);
  }
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow,istream& input,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ostream& oss){
    xstructure xstr_in(input,IOAFLOW_AUTO);
    return CleavageEnergyCalculation(vpflow,xstr_in,aflags,kflags,vflags,oss);
  }
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow,const xstructure& xstr_in,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ostream& oss){
    ofstream FileMESSAGE;
    return CleavageEnergyCalculation(vpflow,xstr_in,aflags,kflags,vflags,FileMESSAGE,oss);
  }
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow,istream& input,ofstream& FileMESSAGE,ostream& oss){
    xstructure xstr_in(input,IOAFLOW_AUTO);
    return CleavageEnergyCalculation(vpflow,xstr_in,FileMESSAGE,oss);
  }
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow,const xstructure& xstr_in,ofstream& FileMESSAGE,ostream& oss){
    _aflags aflags; aflags.Directory=".";
    _kflags kflags;
    _vflags vflags;
    return CleavageEnergyCalculation(vpflow,xstr_in,aflags,kflags,vflags,FileMESSAGE,oss);
  }
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow,istream& input,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss){
    xstructure xstr_in(input,IOAFLOW_AUTO);
    return CleavageEnergyCalculation(vpflow,xstr_in,aflags,kflags,vflags,FileMESSAGE,oss);
  }
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow,const xstructure& xstr_in,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss){
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;
    std::streamsize prec = 8;
    bool check_min_dist=true; //turn off if it gets too slow
    int count_check_min_dist=0;

    message << aflow::Banner("BANNER_NORMAL");
    pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_RAW_);  //first to screen (not logged, file not opened)
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " starting" << endl;}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - get conventional cell
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " get conventional cell" << endl;}

    xstructure xstr_bulk(xstr_in);
    double min_dist=xstr_bulk.dist_nn_min;
    if(min_dist==AUROSTD_NAN){min_dist=SYM::minimumDistance(xstr_bulk);}
    double min_dist_orig=min_dist;
    if(check_min_dist){ //sanity check as we rotate structure/atoms
      min_dist=xstr_bulk.MinDist();
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " mindist[" << count_check_min_dist++ << "]=" << min_dist << endl;}
      if(!aurostd::isequal(min_dist_orig,min_dist)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Minimum distance changed",_VALUE_ERROR_);}
    }

    bool convert_sconv=true;
    if(convert_sconv){xstr_bulk=Standard_Conventional_UnitCellForm(xstr_bulk);xstr_bulk.clean();}  //best to work with standard conventional unitcell //DX20191220 - uppercase to lowercase clean
    if(check_min_dist){ //sanity check as we rotate structure/atoms
      min_dist=xstr_bulk.MinDist();
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " mindist[" << count_check_min_dist++ << "]=" << min_dist << endl;}
      if(!aurostd::isequal(min_dist_orig,min_dist)){
        //throw a warning here instead, minimum distance MIGHT change with sconv
        //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Minimum distance changed",_VALUE_ERROR_);
        message << "Minimum distance changed (sprim -> sconv)";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_WARNING_);
        min_dist_orig=min_dist;
      }
    }
    xstr_bulk.ReScale(1.0);
    xstr_bulk.ShiftOriginToAtom(0);xstr_bulk.origin=0.0; //reset origin
    xstr_bulk.BringInCell();
    xstr_bulk.clean();  //clear origin! //DX20191220 - uppercase to lowercase clean
    if(LDEBUG) {xstr_bulk.write_DEBUG_flag=TRUE;}
    //xstr_bulk.coord_flag=_COORDS_CARTESIAN_;  //much more accurate for this type of calculation

    message << "structure(standard conventional)=";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
    message << xstr_bulk << endl;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_RAW_);

    if(check_min_dist){ //sanity check as we rotate structure/atoms
      min_dist=xstr_bulk.MinDist();
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " mindist[" << count_check_min_dist++ << "]=" << min_dist << endl;}
      if(!aurostd::isequal(min_dist_orig,min_dist)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Minimum distance changed",_VALUE_ERROR_);}
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - get conventional cell
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - get sym info for structure (mostly for FPOSMatch)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    xstructure xstr_sym(xstr_bulk);  //make a copy so we don't carry around all the symmetry in memory as we make copies of a
    pflow::PerformFullSymmetry(xstr_sym);
    xstr_bulk.dist_nn_min=xstr_sym.dist_nn_min;
    xstr_bulk.sym_eps=xstr_sym.sym_eps;
    bool skew=SYM::isLatticeSkewed(xstr_bulk.lattice,xstr_bulk.dist_nn_min,xstr_bulk.sym_eps);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " xstr_bulk.dist_nn_min=" << xstr_bulk.dist_nn_min << endl;
      cerr << __AFLOW_FUNC__ << " xstr_bulk.sym_eps=" << xstr_bulk.sym_eps << endl;
      cerr << __AFLOW_FUNC__ << " skew=" << skew << endl;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - get sym info for structure (mostly for FPOSMatch)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - read flags
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " reading flags" << endl;}

    int h_i=1,k_i=1,l_i=1;  //hkl of interest //needed for dirname
    xvector<int> hkl_i;     //needed for dirname
    int relaxation_layers=DEFAULT_SURFACE_LAYERS;      //number of relaxation layers
    bool spin_off=false;

    vector<string> tokens;
    aurostd::string2tokens(vpflow.getattachedscheme("CREATE_SLAB::PLANE_INTEREST"),tokens,",");
    if(tokens.size()==3){
      h_i=aurostd::string2utype<int>(tokens[0]);
      k_i=aurostd::string2utype<int>(tokens[1]);
      l_i=aurostd::string2utype<int>(tokens[2]);
    }
    hkl_i[1]=h_i;hkl_i[2]=k_i;hkl_i[3]=l_i;
    if(hkl_i[1]==0 && hkl_i[2]==0 && hkl_i[3]==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"hkl_i=(0,0,0)",_INPUT_ERROR_);}
    string relaxation_layers_string=vpflow.getattachedscheme("CLEAVAGE_ENERGY::RELAXATION_LAYERS"); //step size
    if(aurostd::isfloat(relaxation_layers_string)){
      int _relaxation_layers=aurostd::string2utype<int>(relaxation_layers_string);
      if(_relaxation_layers>0){relaxation_layers=_relaxation_layers;}
    }
    spin_off=vpflow.flag("CLEAVAGE_ENERGY::SPIN_OFF");

    std::streamsize prec_original = message.precision(); //original
    std::ios_base::fmtflags ff_original = message.flags();  //original
    message.precision(prec);
    message.unsetf(std::ios_base::floatfield);

    message << "relaxation_layers=" << relaxation_layers;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
    message << "spin=" << (spin_off?"OFF":"ON");pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);

    message.precision(prec_original); //set back
    message.flags(ff_original); //set back

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - read flags
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - defining hkl normals
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " defining HKL normals" << endl;}

    xvector<double> n_i=slab::HKLPlane2Normal(xstr_bulk.lattice,hkl_i);
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " n_i[h=" << hkl_i << "]=" << n_i << endl;}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - defining hkl normals
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - run bulk (with our parameters)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int dir_count=0,dir_count_total=2;  //two total dirs, bulk and slab
    bool all_vasp_done=true;
    string FileName;

    stringstream new_AflowIn_ss;
    stringstream arun_dirname;
    _xvasp xvasp;

    xvasp.clear();
    //create directory name
    xvasp.Directory=aflags.Directory;
    xvasp.AVASP_arun=true;
    xvasp.AVASP_arun_mode="CLEAVENG"; //generalized stacking fault energy
    arun_dirname.str("");
    arun_dirname << std::setfill('0') << std::setw(aurostd::getZeroPadding(dir_count_total)) << dir_count+1 << "_";
    arun_dirname << "iHKL" << aurostd::utype2string(hkl_i[1]) << aurostd::utype2string(hkl_i[2]) << aurostd::utype2string(hkl_i[3]) << "-";
    arun_dirname << "bulk";
    xvasp.AVASP_arun_runname=arun_dirname.str();
    AVASP_populateXVASP(aflags,kflags,vflags,xvasp);

    all_vasp_done=true;
    FileName=xvasp.Directory + "/" + DEFAULT_AFLOW_QMVASP_OUT;
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking for: " << FileName << endl;}
    if(!aurostd::EFileExist(FileName)){all_vasp_done=false;}

    bool convert_sprim=true;
    if(!all_vasp_done){
      message << "Creating bulk VASP run";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);

      xvasp.clear();

      //load in xstructure
      xvasp.str.clear(); //DX20191220 - uppercase to lowercase clear
      xvasp.str=xstr_bulk;
      if(convert_sprim){
        xvasp.str=Standard_Primitive_UnitCellForm(xvasp.str);
        xvasp.str.ReScale(1.0);
        xvasp.str.ShiftOriginToAtom(0);xvasp.str.origin=0.0; //reset origin
        xvasp.str.BringInCell();
        xvasp.str.clean();  //clear origin! //DX20191220 - uppercase to lowercase clean
      }
      xvasp.str.write_DEBUG_flag=FALSE;  //FORCE

      //create directory name
      xvasp.Directory=aflags.Directory;
      xvasp.AVASP_arun=true;
      xvasp.AVASP_arun_mode="CLEAVENG"; //generalized stacking fault energy
      arun_dirname.str("");
      arun_dirname << std::setfill('0') << std::setw(aurostd::getZeroPadding(dir_count_total)) << dir_count+1 << "_";
      arun_dirname << "iHKL" << aurostd::utype2string(hkl_i[1]) << aurostd::utype2string(hkl_i[2]) << aurostd::utype2string(hkl_i[3]) << "-";
      arun_dirname << "bulk";
      xvasp.AVASP_arun_runname=arun_dirname.str();
      AVASP_populateXVASP(aflags,kflags,vflags,xvasp);

      setPreserveUnitCell(xvasp);
      //set k-points to 11x11x1 as per 10.1103/PhysRevLett.114.165502, this is just a quick fix for now
      if(1){
        xvasp.aopts.flag("FLAG::KPOINTS_IMPLICIT",FALSE);
        //xvasp.aopts.flag("FLAG::KPOINTS_EXPLICIT",TRUE);
        xvasp.aopts.flag("FLAG::KPOINTS_EXPLICIT_START_STOP",TRUE);
        xvasp.AVASP_KPOINTS_EXPLICIT_START_STOP.str(""); //clear
        xvasp.AVASP_KPOINTS_EXPLICIT_START_STOP << "KPOINTS" << endl;
        xvasp.AVASP_KPOINTS_EXPLICIT_START_STOP << "0" << endl;
        xvasp.AVASP_KPOINTS_EXPLICIT_START_STOP << "Gamma" << endl;
        xvasp.AVASP_KPOINTS_EXPLICIT_START_STOP << "11 11 11" << endl;
        xvasp.AVASP_KPOINTS_EXPLICIT_START_STOP << "0 0 0" << endl;
      }

      //relax ions ONLY
      xvasp.aopts.flag("AFLOWIN_FLAG::RELAX_TYPE",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::RELAX_TYPE","ALL");
      if(spin_off){xvasp.aopts.flag("FLAG::AVASP_SPIN",FALSE);}

      //make aflow.in
      AVASP_MakeSingleAFLOWIN(xvasp,new_AflowIn_ss,true); //false,-1,false);  //don't write/print and hence don't pthread
    }
    dir_count++;

    if(!all_vasp_done){
      message << "Now waiting for the bulk VASP run";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_COMPLETE_);
      return;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - run bulk (with our parameters)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - grab relaxed bulk
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    xstructure xstr_relaxed=KBIN::GetMostRelaxedStructure(xvasp.Directory);
    message << "structure(bulk_relaxed)=";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
    message << xstr_relaxed << endl;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_RAW_);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - grab relaxed bulk
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - create slab
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    message << "Creating slab";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);

    //xvector<int> hkl_i;
    int total_layers = 0;
    xmatrix<double> rotation;
    xstructure xstr_slablattice;
    vector<int> sc2pcMap_slab,pc2scMap_slab;

    xstructure xstr_slab=slab::CreateSlab_SurfaceLattice(vpflow,xstr_relaxed,hkl_i,total_layers,rotation,xstr_slablattice,sc2pcMap_slab,pc2scMap_slab,aflags,FileMESSAGE,AUROSTD_MAX_DOUBLE,oss);

    if(total_layers%2!=0){message << "total_layers is odd, it is better to pick an even number (top vs. bottom)";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_WARNING_);}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - create slab
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - resolve layers count
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " resolving layers count" << endl;}

    double d_spacing=slab::getSpacingHKLPlane(xstr_bulk,hkl_i); //aurostd::modulus(xstr_slab.lattice(1))/sqrt(h_s*h_s+k_s*k_s+l_s*l_s);
    double d_layers=slab::getDistanceBetweenImages(xstr_bulk,n_i,false); //this depends on UN-ROTATED lattice
    double d_cells=slab::getDistanceBetweenImages(xstr_bulk,n_i,true); //go outside cell
    int layers_per_cell=(int)(d_cells/d_layers);  //floor
    int supercell_layers=(total_layers+layers_per_cell-1)/layers_per_cell;  //ceil //(double)total_layers;
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " n_i[h=" << hkl_i << "]=" << n_i << endl;
      cerr << __AFLOW_FUNC__ << " d_spacing=" << d_spacing << endl;
      cerr << __AFLOW_FUNC__ << " d_layers=" << d_layers << endl;
      cerr << __AFLOW_FUNC__ << " d_cells=" << d_cells << endl;
      cerr << __AFLOW_FUNC__ << " abs(d_layers-d_cells)=" << abs(d_layers-d_cells) << endl;
      cerr << __AFLOW_FUNC__ << " layers_per_cell=" << layers_per_cell << endl;
      cerr << __AFLOW_FUNC__ << " supercell_layers=" << supercell_layers << endl;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - resolve layers count
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - identify selective dynamics
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " identifying selective dynamics" << endl;}

    //k in atom.ijk will go from 0 to supercell_layers-1

    //sometimes k_min==0, sometimes k_min==1 (look at GetSuperCell())
    //get k_min and k_max first
    int k=0,k_min=xstr_slab.atoms.size(),k_max=-xstr_slab.atoms.size();
    for(uint i=0;i<xstr_slab.atoms.size();i++){
      const _atom& atom=xstr_slab.atoms[i];
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " atom.fpos=" << atom.fpos << ", atom.ijk=" << atom.ijk << endl;}
      k=atom.ijk[3];
      if(k<k_min){k_min=k;}
      if(k>k_max){k_max=k;}
    }
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " k_min=" << k_min << ", k_max=" << k_max << endl;}
    if(k_max-k_min+1!=supercell_layers){  //test of stupidity
      message << "k_max-k_min!=supercell_layers";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ERROR_);
    }

    uint count_total_relaxing=0,count_bottom_relaxing=0,count_top_relaxing=0; //tests of stupidity
    bool relax_bottom=false,relax_top=false;

    xstr_slab.isd=true; //set selective dynamics
    for(uint i=0;i<xstr_slab.atoms.size();i++){
      _atom& atom=xstr_slab.atoms[i];
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " atom.fpos=" << atom.fpos << ", atom.ijk=" << atom.ijk << endl;}
      k=atom.ijk[3];
      relax_bottom=(k<(k_min+relaxation_layers));
      relax_top=(k>(k_min+(supercell_layers-relaxation_layers-1)));
      if(relax_bottom || relax_top){  //relax
        //atom.sd="TTT";  //default, allow relaxation
        atom.sd="FFT";  //default, allow relaxation //CORRECTION: only relax in z direction: Vitek Phil Mag 18, 773-786 (1968), also http://theory.cm.utexas.edu/forum/viewtopic.php?t=3301
        count_total_relaxing++;
        if(relax_bottom){count_bottom_relaxing++;}
        else if(relax_top){count_top_relaxing++;}
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " atom.ijk=" << atom.ijk << " will RELAX" << endl;}
      } else {
        atom.sd="FFF";  //keep fixed
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " atom.ijk=" << atom.ijk << " is FIXED" << endl;}
      }
    }
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " count_bottom_relaxing=" << count_bottom_relaxing << endl;
      cerr << __AFLOW_FUNC__ << " count_top_relaxing=" << count_top_relaxing << endl;
      cerr << __AFLOW_FUNC__ << " count_bottom_relaxing+count_top_relaxing=" << count_bottom_relaxing+count_top_relaxing << endl;
      cerr << __AFLOW_FUNC__ << " count_total_relaxing=" << count_total_relaxing << endl;
    }
    if(count_total_relaxing!=count_bottom_relaxing+count_top_relaxing){
      message << "count_total_relaxing!=count_bottom_relaxing+count_top_relaxing (" << count_total_relaxing << "!=" << count_bottom_relaxing+count_top_relaxing << ")";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ERROR_); //CO20190226
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - identify selective dynamics
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - run slab
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    xvasp.clear();
    //create directory name
    xvasp.Directory=aflags.Directory;
    xvasp.AVASP_arun=true;
    xvasp.AVASP_arun_mode="CLEAVENG"; //generalized stacking fault energy
    arun_dirname.str("");
    arun_dirname << std::setfill('0') << std::setw(aurostd::getZeroPadding(dir_count_total)) << dir_count+1 << "_";
    arun_dirname << "iHKL" << aurostd::utype2string(hkl_i[1]) << aurostd::utype2string(hkl_i[2]) << aurostd::utype2string(hkl_i[3]) << "-";
    arun_dirname << "slab";
    xvasp.AVASP_arun_runname=arun_dirname.str();
    AVASP_populateXVASP(aflags,kflags,vflags,xvasp);

    all_vasp_done=true;
    FileName=xvasp.Directory + "/" + DEFAULT_AFLOW_QMVASP_OUT;
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking for: " << FileName << endl;}
    if(!aurostd::EFileExist(FileName)){all_vasp_done=false;}

    if(!all_vasp_done){
      message << "Creating slab VASP run";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);

      xvasp.clear();
      //load in xstructure
      xvasp.str.clear(); //DX20191220 - uppercase to lowercase clear
      xvasp.str=xstr_slab;
      xvasp.str.write_DEBUG_flag=FALSE;  //FORCE

      //create directory name
      xvasp.Directory=aflags.Directory;
      xvasp.AVASP_arun=true;
      xvasp.AVASP_arun_mode="CLEAVENG"; //generalized stacking fault energy
      arun_dirname.str("");
      arun_dirname << std::setfill('0') << std::setw(aurostd::getZeroPadding(dir_count_total)) << dir_count+1 << "_";
      arun_dirname << "iHKL" << aurostd::utype2string(hkl_i[1]) << aurostd::utype2string(hkl_i[2]) << aurostd::utype2string(hkl_i[3]) << "-";
      arun_dirname << "slab";
      xvasp.AVASP_arun_runname=arun_dirname.str();
      AVASP_populateXVASP(aflags,kflags,vflags,xvasp);

      setPreserveUnitCell(xvasp);
      //set k-points to 11x11x1 as per 10.1103/PhysRevLett.114.165502, this is just a quick fix for now
      if(1){
        xvasp.aopts.flag("FLAG::KPOINTS_IMPLICIT",FALSE);
        //xvasp.aopts.flag("FLAG::KPOINTS_EXPLICIT",TRUE);
        xvasp.aopts.flag("FLAG::KPOINTS_EXPLICIT_START_STOP",TRUE);
        xvasp.AVASP_KPOINTS_EXPLICIT_START_STOP.str(""); //clear
        xvasp.AVASP_KPOINTS_EXPLICIT_START_STOP << "KPOINTS" << endl;
        xvasp.AVASP_KPOINTS_EXPLICIT_START_STOP << "0" << endl;
        xvasp.AVASP_KPOINTS_EXPLICIT_START_STOP << "Gamma" << endl;
        xvasp.AVASP_KPOINTS_EXPLICIT_START_STOP << "11 11 1" << endl;
        xvasp.AVASP_KPOINTS_EXPLICIT_START_STOP << "0 0 0" << endl;
      }

      //relax ions ONLY
      xvasp.aopts.flag("AFLOWIN_FLAG::RELAX_TYPE",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::RELAX_TYPE","IONS");
      if(spin_off){xvasp.aopts.flag("FLAG::AVASP_SPIN",FALSE);}

      //make aflow.in
      AVASP_MakeSingleAFLOWIN(xvasp,new_AflowIn_ss,true); //false,-1,false);  //don't write/print and hence don't pthread

    }
    dir_count++;

    if(!all_vasp_done){
      message << "Now waiting for the slab VASP run";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_COMPLETE_);
      return;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - run slab
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - postprocessing
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - postprocessing
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  }
} // namespace pflow

namespace pflow {
  bool findClosedPackingPlane(istream& input){xstructure xstr_in(input,IOAFLOW_AUTO);return findClosedPackingPlane(xstr_in);} //CO20191110
  bool findClosedPackingPlane(const xstructure& xstr_in){ //CO20191110
    bool LDEBUG=(FALSE || XHOST.DEBUG);

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " xstr_in=" << endl;cerr << xstr_in << endl;}

    xstructure xstr_grid(xstr_in);
    xstr_grid.GenerateGridAtoms(3); //-3:3 in every direction
    const deque<_atom>& grid_atoms=xstr_grid.grid_atoms;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " grid_atoms=" << endl;for(uint i=0;i<grid_atoms.size();i++){cerr << "i=" << i << ": cpos=" << grid_atoms[i].cpos << endl;}}

    aurostd::xcombos xc(grid_atoms.size(),4,'C');  //need to check coplanarity: requires 4 points
    vector<int> combo;

    while(xc.increment()){
      combo=xc.getCombo();
      const _atom& a1=grid_atoms[combo[0]];
      const _atom& a2=grid_atoms[combo[1]];
      const _atom& a3=grid_atoms[combo[2]];
      const _atom& a4=grid_atoms[combo[3]];

      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " checking coplanarity of:" << endl;
        cerr << a1.cpos << endl;
        cerr << a2.cpos << endl;
        cerr << a3.cpos << endl;
        cerr << a4.cpos << endl;
      }


      //check if they are coplanar
      //use formula from https://en.wikipedia.org/wiki/Coplanarity
      if(aurostd::isequal( aurostd::scalar_product( aurostd::vector_product( a2.cpos-a1.cpos , a4.cpos-a1.cpos ) , a3.cpos-a1.cpos ) , 0.0)){
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " found coplanar set" << endl;


        }
      }
    }

    return true;
  }
} // namespace pflow

// ***************************************************************************
// ParseChemFormula
// ***************************************************************************
void ParseChemFormula(string& ChemFormula,vector<string>& ChemNameRef,vector<float>& ChemConcRef) {
  //ChemFormula=MgB2.3 -> ChemNameRef=["Mg","B"] and ChemConcRef=[1,2.3]

  //input:
  //ChemFormula

  //output:
  //ChemNameRef, ChemConcRef

  //wahyu@alumni.duke.edu

  uint nchar;
  //  int i,j,k,L,itmp;
  float AtomConc;
  string AtomSymbol;

  ChemNameRef.clear(); ChemConcRef.clear();
  while(ChemFormula.size()>0) {
    nchar=3;
    if(ChemFormula.size()>nchar-1) {
      //check 3-character atom symbol
      ParseChemFormulaIndividual(3,ChemFormula,AtomSymbol,AtomConc);
      if(AtomSymbol!="") {
        ChemNameRef.push_back(AtomSymbol);
        ChemConcRef.push_back(AtomConc);
      }//if valid 3-character atomic symbol found
      else {
        //check 2-character atom symbol
        ParseChemFormulaIndividual(2,ChemFormula,AtomSymbol,AtomConc);
        if(AtomSymbol!="") {
          ChemNameRef.push_back(AtomSymbol);
          ChemConcRef.push_back(AtomConc);
        }//if valid 2-character atomic symbol found
        else {
          //check 1-character atom symbol
          ParseChemFormulaIndividual(1,ChemFormula,AtomSymbol,AtomConc);
          if(AtomSymbol!="") {
            ChemNameRef.push_back(AtomSymbol);
            ChemConcRef.push_back(AtomConc);
          }//if valid 1-character atomic symbol found
        }//else nchar 1
      }//else nchar 2
    }//if ChemFormula.size()>2
    else {
      nchar=2;
      if(ChemFormula.size()>nchar-1) {
        //check 2-character atom symbol
        ParseChemFormulaIndividual(2,ChemFormula,AtomSymbol,AtomConc);
        if(AtomSymbol!="") {
          ChemNameRef.push_back(AtomSymbol);
          ChemConcRef.push_back(AtomConc);
        }//if valid 2-character atomic symbol found
        else {
          //check 1-character atom symbol
          ParseChemFormulaIndividual(1,ChemFormula,AtomSymbol,AtomConc);
          if(AtomSymbol!="") {
            ChemNameRef.push_back(AtomSymbol);
            ChemConcRef.push_back(AtomConc);
          }//if valid 1-character atomic symbol found
        }//else nchar 1
      }//if ChemFormula.size()>1
      else {
        nchar=1;
        ParseChemFormulaIndividual(1,ChemFormula,AtomSymbol,AtomConc);
        if(AtomSymbol!="") {
          ChemNameRef.push_back(AtomSymbol);
          ChemConcRef.push_back(AtomConc);
        }//if valid 1-character atomic symbol found    
      }//if ChemFormula.size()==1
    }
  }//while not empty

}
// ***************************************************************************
// ParseChemFormulaIndividual
// ***************************************************************************
void ParseChemFormulaIndividual(uint nchar, string& ChemFormula, string& AtomSymbol, float& AtomConc) {
  //get 1 valid pair of atom symbol and its concentration from ChemFormula.
  //start from the beginning of ChemFormula. If found, return the AtomSymbol and AtomConc
  //and truncate it from the ChemFormula.

  //input:
  //nchar     : the number of characters in the atom symbol

  //wahyu@alumni.duke.edu

  uint iret;
  int i;
  string sN,sC;

  AtomSymbol="";
  AtomConc=0.0;
  sN=ChemFormula.substr(0,nchar);
  iret=GetAtomNumber(sN);
  if(iret>0) {//valid atomic symbol
    AtomSymbol=sN;
    if(ChemFormula.size()==nchar) {//last atom with no concentration
      ChemFormula="";
      AtomConc=1.0;
    }
    else {
      ChemFormula=ChemFormula.substr(nchar,ChemFormula.size());
      if((ChemFormula[0]>'0' && ChemFormula[0]<'9') || ChemFormula[0]=='.') {
        //parse the concentration
        for(i=0;i<(int)ChemFormula.size();i++) {
          if(!((ChemFormula[i]>'0' && ChemFormula[i]<'9') || ChemFormula[i]=='.')) break;
        }
        sC=ChemFormula.substr(0,i);
        if(sC[0]=='.') sC="0"+sC;//append 0 in case the first character is '.'
        if(sC[sC.size()-1]=='.') sC=sC+"0";
        AtomConc=aurostd::string2utype<float>(sC);
        if(i==(int)ChemFormula.size()) ChemFormula="";
        else ChemFormula=ChemFormula.substr(i,ChemFormula.size());
      }
      else  AtomConc=1.0;//set concentration to 1.0 if atom symbol is valid but concentration is not given
    }
  }//if valid atomic symbol found  
}

// ***************************************************************************
// GetXRAY HELPER FUNCTIONS
// ***************************************************************************
// Corey Oses 20190620
namespace pflow {
  void GetXray2ThetaIntensity(const xstructure& str,
      vector<double>& v_twotheta,
      vector<double>& v_intensity,
      double lambda) { //CO20190520 //CO20190620 - v_amplitude can be calculated later
    int num_atoms=str.atoms.size();
    //  vector<string> names=str.GetNames();
    vector<double> dist,sf;
    vector<double> scatt_fact(num_atoms,0.0);
    vector<double> mass(num_atoms,0.0);
    vector<double> twoB_vec(num_atoms,0.0);
    //pflow::GetXray(str,dist,sf,lambda,scatt_fact,mass,twoB_vec);
    vector<vector<double> > ids;
    aurostd::matrix<double> data;  //CO20200404 pflow::matrix()->aurostd::matrix()
    pflow::GetXrayData(str,dist,sf,scatt_fact,mass,twoB_vec,ids,data,lambda); //CO20190620 - intmax can be grabbed later

    v_twotheta.clear();

    double tol=XRAY_THETA_TOL;

    for(uint i=0;i<data.size();i++) {
      double theta=0;
      if(data[i][3]>0) {
        double term=lambda/(2.0*data[i][3]);
        if(term<=1) {
          theta=asin(term);
          theta=theta*360.0/TWOPI; // rad->degrees
          if(theta>tol){
            v_twotheta.push_back(2.0*theta);
            v_intensity.push_back(data[i][4]);
          } // if theta<tol
        } // if term<=1
      } // if dist>0
    } // for
  }

  vector<uint> GetXrayPeaks(const xstructure& str,
      vector<double>& v_twotheta,
      vector<double>& v_intensity,
      vector<double>& v_intensity_smooth,
      double lambda) { //CO20190520 //CO20190620 - v_peaks_amplitude not needed
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " input str=" << endl;cerr << str << endl;}
    GetXray2ThetaIntensity(str,v_twotheta,v_intensity,lambda);  //v_amplitude can be grabbed later
    return GetXrayPeaks(v_twotheta,v_intensity,v_intensity_smooth);
  }
  vector<uint> GetXrayPeaks(const vector<double>& v_twotheta,const vector<double>& v_intensity,vector<double>& v_intensity_smooth) { //CO20190520 //CO20190620 - v_peaks_amplitude not needed
    bool LDEBUG=(FALSE || XHOST.DEBUG);

    if(v_twotheta.size()<2){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"v_twotheta.size()<2",_VALUE_ILLEGAL_);}
    if(v_twotheta.size()!=v_intensity.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"v_twotheta.size()!=v_intensity.size()",_VALUE_ILLEGAL_);}

    xvector<double> xv_intensity=aurostd::vector2xvector<double>(v_intensity),xv_intensity_smooth;
    uint smoothing_iterations=4,avg_window=4;int width_maximum=1;double significance_multiplier=1.0;  //defaults

    //fix avg_window to derive from degree separation, not points
    double twotheta_smoothing=5.0;  //avg over 5 degrees
    double twotheta_diff=(v_twotheta[1]-v_twotheta[0]);
    double davg_window=twotheta_smoothing/twotheta_diff;
    avg_window=max(avg_window,(uint)davg_window);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " twotheta_smoothing=" << twotheta_smoothing << endl;
      cerr << __AFLOW_FUNC__ << " twotheta_diff=" << twotheta_diff << endl;
      cerr << __AFLOW_FUNC__ << " davg_window=" << davg_window << endl;
      cerr << __AFLOW_FUNC__ << " avg_window=" << avg_window << endl;
    }

    //CO COME BACK, plug in 20-120 range
    vector<int> _peak_indices=aurostd::getPeaks(xv_intensity,xv_intensity_smooth,smoothing_iterations,avg_window,width_maximum,significance_multiplier);  //indices wrt xvector (starts at 1), need to convert
    v_intensity_smooth=aurostd::xvector2vector<double>(xv_intensity_smooth);
    vector<uint> peak_indices;for(uint i=0;i<_peak_indices.size();i++){peak_indices.push_back(_peak_indices[i]-xv_intensity.lrows);} //shift indices for xvector -> vector conversion

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " _peak_indices=" << aurostd::joinWDelimiter(_peak_indices,",") << endl;
      cerr << __AFLOW_FUNC__ << "  peak_indices=" << aurostd::joinWDelimiter(peak_indices,",") << endl;
    }

    return peak_indices;

    //[CO20190620 - moved to aurostd::getPeaks()]v_intensity_smooth.clear();
    //[CO20190620 - moved to aurostd::getPeaks()]v_peaks_twotheta.clear();
    //[CO20190620 - moved to aurostd::getPeaks()]v_peaks_intensity.clear();
    //[CO20190620 - moved to aurostd::getPeaks()]
    //[CO20190620 - moved to aurostd::getPeaks()]//using method outlined here: https://dsp.stackexchange.com/questions/1302/peak-detection-approach
    //[CO20190620 - moved to aurostd::getPeaks()]//raw data is X
    //[CO20190620 - moved to aurostd::getPeaks()]//smooth data via moving average to get Y
    //[CO20190620 - moved to aurostd::getPeaks()]//stddev(X-Y) is sigma
    //[CO20190620 - moved to aurostd::getPeaks()]//detect peaks when (X-Y)>multiplier*sigma
    //[CO20190620 - moved to aurostd::getPeaks()]
    //[CO20190620 - moved to aurostd::getPeaks()]//smooth data
    //[CO20190620 - moved to aurostd::getPeaks()]uint smoothing_iterations=4;int avg_window=4;
    //[CO20190620 - moved to aurostd::getPeaks()]xvector<double> xv_intensity=aurostd::vector2xvector<double>(v_intensity);
    //[CO20190620 - moved to aurostd::getPeaks()]xvector<double> xv_intensity_smooth=xv_intensity;
    //[CO20190620 - moved to aurostd::getPeaks()]for(uint i=0;i<smoothing_iterations;i++){xv_intensity_smooth=aurostd::moving_average(xv_intensity_smooth,avg_window);}
    //[CO20190620 - moved to aurostd::getPeaks()]xvector<double> diff=xv_intensity-xv_intensity_smooth;
    //[CO20190620 - moved to aurostd::getPeaks()]double sigma=aurostd::stddev(diff);
    //[CO20190620 - moved to aurostd::getPeaks()]double multiplier=1;
    //[CO20190620 - moved to aurostd::getPeaks()]
    //[CO20190620 - moved to aurostd::getPeaks()]bool local_maximum=false;
    //[CO20190620 - moved to aurostd::getPeaks()]bool significant=false;
    //[CO20190620 - moved to aurostd::getPeaks()]for(int i=xv_intensity.lrows;i<=xv_intensity.urows;i++){
    //[CO20190620 - moved to aurostd::getPeaks()]  local_maximum=(i>xv_intensity.lrows && i<xv_intensity.urows && xv_intensity[i]>xv_intensity[i-1] && xv_intensity[i]>xv_intensity[i+1]);
    //[CO20190620 - moved to aurostd::getPeaks()]  significant=(diff[i]>multiplier*sigma);
    //[CO20190620 - moved to aurostd::getPeaks()]  if(local_maximum && significant){
    //[CO20190620 - moved to aurostd::getPeaks()]    v_peaks_twotheta.push_back(v_twotheta[i-xv_intensity.lrows]);
    //[CO20190620 - moved to aurostd::getPeaks()]    v_peaks_intensity.push_back(v_intensity[i-xv_intensity.lrows]);
    //[CO20190620 - moved to aurostd::getPeaks()]    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " PEAK[two-theta=" << v_twotheta[i-xv_intensity.lrows] << "]=" << xv_intensity[i] << endl;}
    //[CO20190620 - moved to aurostd::getPeaks()]  }
    //[CO20190620 - moved to aurostd::getPeaks()]}
    //[CO20190620 - moved to aurostd::getPeaks()]
    //[CO20190620 - moved to aurostd::getPeaks()]v_intensity_smooth=aurostd::xvector2vector(xv_intensity_smooth);
    //[CO20190620 - moved to aurostd::getPeaks()]
    //[CO20190620 - moved to aurostd::getPeaks()]if(0){  //don't bother sorting
    //[CO20190620 - moved to aurostd::getPeaks()]  for(uint i=0;i<v_peaks_twotheta.size();i++){
    //[CO20190620 - moved to aurostd::getPeaks()]    for(uint j=i+1;j<v_peaks_twotheta.size();j++){
    //[CO20190620 - moved to aurostd::getPeaks()]      if(v_peaks_intensity[j]>v_peaks_intensity[i]){  //CO20190620 - use intensity vs. amplitude
    //[CO20190620 - moved to aurostd::getPeaks()]        double tmp_twotheta=v_peaks_twotheta[i];
    //[CO20190620 - moved to aurostd::getPeaks()]        double tmp_intensity=v_peaks_intensity[i];
    //[CO20190620 - moved to aurostd::getPeaks()]        
    //[CO20190620 - moved to aurostd::getPeaks()]        v_peaks_twotheta[i]=v_peaks_twotheta[j];
    //[CO20190620 - moved to aurostd::getPeaks()]        v_peaks_intensity[i]=v_peaks_intensity[j];
    //[CO20190620 - moved to aurostd::getPeaks()]        
    //[CO20190620 - moved to aurostd::getPeaks()]        v_peaks_twotheta[j]=tmp_twotheta;
    //[CO20190620 - moved to aurostd::getPeaks()]        v_peaks_intensity[j]=tmp_intensity;
    //[CO20190620 - moved to aurostd::getPeaks()]      }
    //[CO20190620 - moved to aurostd::getPeaks()]    }
    //[CO20190620 - moved to aurostd::getPeaks()]  }
    //[CO20190620 - moved to aurostd::getPeaks()]}
  }
} // namespace pflow

//[CO20190629 - replaced with aurostd::compareVecElements<double>]// ***************************************************************************
//[CO20190629 - replaced with aurostd::compareVecElements<double>]// PrintXRAY ids_cmp
//[CO20190629 - replaced with aurostd::compareVecElements<double>]// ***************************************************************************
//[CO20190629 - replaced with aurostd::compareVecElements<double>]// This function sorts by theta (reverse sort by distance)
//[CO20190629 - replaced with aurostd::compareVecElements<double>]class ids_cmp{
//[CO20190629 - replaced with aurostd::compareVecElements<double>]public:
//[CO20190629 - replaced with aurostd::compareVecElements<double>]  int operator()(const vector<double>& a, const vector<double>& b)
//[CO20190629 - replaced with aurostd::compareVecElements<double>]  {return a[0]>b[0];} // Sorts in increasing order.
//[CO20190629 - replaced with aurostd::compareVecElements<double>]};

//[CO20190629 - replaced with aurostd::compareVecElements<int>]// ***************************************************************************
//[CO20190629 - replaced with aurostd::compareVecElements<int>]// PrintXRAY hkl_cmp
//[CO20190629 - replaced with aurostd::compareVecElements<int>]// ***************************************************************************
//[CO20190629 - replaced with aurostd::compareVecElements<int>]// This function sorts hkl
//[CO20190629 - replaced with aurostd::compareVecElements<int>]class hkl_cmp{
//[CO20190629 - replaced with aurostd::compareVecElements<int>]public:
//[CO20190629 - replaced with aurostd::compareVecElements<int>]  int operator()(const vector<int>& a, const vector<int>& b)
//[CO20190629 - replaced with aurostd::compareVecElements<int>]  {
//[CO20190629 - replaced with aurostd::compareVecElements<int>]    if(a[0]!=b[0]){return a[0]>b[0];}
//[CO20190629 - replaced with aurostd::compareVecElements<int>]    if(a[1]!=b[1]){return a[1]>b[1];}
//[CO20190629 - replaced with aurostd::compareVecElements<int>]    //if(a[2]!=b[2]){return a[2]>b[2];}
//[CO20190629 - replaced with aurostd::compareVecElements<int>]    return a[2]>b[2];   //return something...
//[CO20190629 - replaced with aurostd::compareVecElements<int>]    //[CO20190620 - what if h, k,or l is bigger than 10?]int na=a[0]*100+a[1]*10+a[2];
//[CO20190629 - replaced with aurostd::compareVecElements<int>]    //[CO20190620 - what if h, k,or l is bigger than 10?]int nb=b[0]*100+b[1]*10+b[2];
//[CO20190629 - replaced with aurostd::compareVecElements<int>]    //[CO20190620 - what if h, k,or l is bigger than 10?]return na>nb;
//[CO20190629 - replaced with aurostd::compareVecElements<int>]  }
//[CO20190629 - replaced with aurostd::compareVecElements<int>]};
//[CO20190629 - replaced with aurostd::compareVecElements<int>]//      int t=1,f=1;
//[CO20190629 - replaced with aurostd::compareVecElements<int>]//      if(a[0]<b[0]) {
//[CO20190629 - replaced with aurostd::compareVecElements<int>]//	return false;
//[CO20190629 - replaced with aurostd::compareVecElements<int>]//      }
//[CO20190629 - replaced with aurostd::compareVecElements<int>]//      else {
//[CO20190629 - replaced with aurostd::compareVecElements<int>]//	if(a[1]<b[1]) {
//[CO20190629 - replaced with aurostd::compareVecElements<int>]//	  return false;
//[CO20190629 - replaced with aurostd::compareVecElements<int>]//	}
//[CO20190629 - replaced with aurostd::compareVecElements<int>]//	else {
//[CO20190629 - replaced with aurostd::compareVecElements<int>]//	  if(a[2]<b[2]) {
//[CO20190629 - replaced with aurostd::compareVecElements<int>]//	    return false;
//[CO20190629 - replaced with aurostd::compareVecElements<int>]//	  }
//[CO20190629 - replaced with aurostd::compareVecElements<int>]//	  else {
//[CO20190629 - replaced with aurostd::compareVecElements<int>]//	    return true;
//[CO20190629 - replaced with aurostd::compareVecElements<int>]//	  }
//[CO20190629 - replaced with aurostd::compareVecElements<int>]//	}
//[CO20190629 - replaced with aurostd::compareVecElements<int>]//      }
//[CO20190629 - replaced with aurostd::compareVecElements<int>]//    }

namespace pflow {
  void GetXrayData(const xstructure& str,
      vector<double>& dist,
      vector<double>& sf,
      vector<double>& scatt_fact,
      vector<double>& mass,vector<double>& twoB_vec,
      vector<vector<double> >& ids,
      aurostd::matrix<double>& data, //CO20200404 pflow::matrix()->aurostd::matrix()
      double lambda) { //CO20190520  //CO20190620 - intmax can be grabbed later
    //int num_atoms=str.atoms.size();
    //  vector<string> names=str.GetNames();
    //vector<double> dist,sf;
    //vector<double> scatt_fact(num_atoms,0.0);
    //vector<double> mass(num_atoms,0.0);
    //vector<double> twoB_vec(num_atoms,0.0);
    pflow::GetXray(str,dist,sf,scatt_fact,mass,twoB_vec,lambda);

    double tol=1.0E-5;
    //int w1=4; // int
    //int w2=12; // some doubles
    //int w3=20; // Integrated intensities
    int tlen=sf.size();
    int len=Nint(std::pow((double) tlen,(double) 1.0/3.0));
    int kmx=(len-1)/2; // len should be odd.

    // Sort by theta (reverse sort by distance).
    // Define an id pointer to sort.
    vector<double> v(5);
    // aurostd::matrix<double> ids(tlen,v);  //CO20200404 pflow::matrix()->aurostd::matrix()
    //vector<vector<double> > ids(tlen,v);  //CO20190409
    ids.resize(tlen,v);
    for(int i0=-kmx;i0<=kmx;i0++) {
      for(int i1=-kmx;i1<=kmx;i1++) {
        for(int i2=-kmx;i2<=kmx;i2++) {
          int ii0=i0+kmx;
          int ii1=i1+kmx;
          int ii2=i2+kmx;
          int id=ii2+ii1*len+ii0*len*len;
          ids[id][0]=dist[id];
          ids[id][1]=id;
          ids[id][2]=i0;
          ids[id][3]=i1;
          ids[id][4]=i2;
        } // i2
      } // i1
    } // i0

    //[CO20190629 - waste of a class]sort(ids.begin(),ids.end(),ids_cmp());
    //[CO20190629 - does NOT work, sort ONLY by 0th index]sort(ids.rbegin(),ids.rend(),aurostd::compareVecElements<double>);  //CO20190629 - note that it is in descending order by distance (greater go first)
    //[CO20190629 - rbegin()/rend() != descending sort, this WILL change results]sort(ids.rbegin(),ids.rend(),aurostd::compareVecElement<double>(0));  //CO20190629 - note that it is in descending order by distance (greater go first)
    sort(ids.begin(),ids.end(),aurostd::compareVecElement<double>(0,false));  //CO20190629 - note that it is in descending order by distance (greater go first)

    // Add corrections to all the amplitudes.
    // Get max amplitude for normalizing and percentages.
    double ampmax=1E-8;
    for(int i0=-kmx;i0<=kmx;i0++) {
      for(int i1=-kmx;i1<=kmx;i1++) {
        for(int i2=-kmx;i2<=kmx;i2++) {
          int ii0=i0+kmx;
          int ii1=i1+kmx;
          int ii2=i2+kmx;
          int id1=ii2+ii1*len+ii0*len*len;
          int id=(int) ids[id1][1];
          ii0=(int) ids[id1][2];
          ii1=(int) ids[id1][3];
          ii2=(int) ids[id1][4];
          double theta=0;
          if(dist[id]>0) {
            double term=lambda/(2.0*dist[id]);
            if(term<=1) {
              theta=asin(term);
              if(theta>tol) {
                sf[id]=sf[id]*CorrectionFactor(theta);
                if(sf[id]>ampmax) {ampmax=sf[id];}
              } // if theta>tol
            } // if term<=1
          } // if dist>0
        } // i2
      } // i1
    } // i0

    //[CO20190520 - printing moved to PrintXray()]// Print out all data.
    //[CO20190520 - printing moved to PrintXray()]oss << "Wavelength (Ang) = " << lambda << endl;
    //[CO20190520 - printing moved to PrintXray()]oss << "Atom_Name  ScattFact   Mass(amu)   B(Ang)(DW=exp(-B*sin(theta)^2/lambda^2))" << endl; //CO20190329
    //[CO20190520 - printing moved to PrintXray()]for(uint iat=0;iat<(uint) num_atoms;iat++) {
    //[CO20190520 - printing moved to PrintXray()]  oss <<setw(4)<<iat+1<<" "<<setw(4)<<str.atoms.at(iat).cleanname << " ";
    //[CO20190520 - printing moved to PrintXray()]  if(str.atoms.at(iat).cleanname.length()>1) oss << " ";
    //[CO20190520 - printing moved to PrintXray()]  oss <<setw(w2)<<scatt_fact[iat]<<setw(w2)<<KILOGRAM2AMU*mass[iat]<<setw(w2)<<1E+20*twoB_vec[iat]/2.0<<endl;
    //[CO20190520 - printing moved to PrintXray()]}
    //[CO20190520 - printing moved to PrintXray()]oss << "******************** All data ********************" << endl;
    //[CO20190520 - printing moved to PrintXray()]oss << "2*theta      Intensity            h    k    l    dist         keyword " << endl;
    //[CO20190520 - printing moved to PrintXray()]for(int i0=-kmx;i0<=kmx;i0++) {
    //[CO20190520 - printing moved to PrintXray()]  for(int i1=-kmx;i1<=kmx;i1++) {
    //[CO20190520 - printing moved to PrintXray()]    for(int i2=-kmx;i2<=kmx;i2++) {
    //[CO20190520 - printing moved to PrintXray()]      int ii0=i0+kmx;
    //[CO20190520 - printing moved to PrintXray()]      int ii1=i1+kmx;
    //[CO20190520 - printing moved to PrintXray()]      int ii2=i2+kmx;
    //[CO20190520 - printing moved to PrintXray()]      int id1=ii2+ii1*len+ii0*len*len;
    //[CO20190520 - printing moved to PrintXray()]      int id=(int) ids[id1][1];
    //[CO20190520 - printing moved to PrintXray()]      ii0=(int) ids[id1][2];
    //[CO20190520 - printing moved to PrintXray()]      ii1=(int) ids[id1][3];
    //[CO20190520 - printing moved to PrintXray()]      ii2=(int) ids[id1][4];
    //[CO20190520 - printing moved to PrintXray()]      double theta=0;
    //[CO20190520 - printing moved to PrintXray()]      if(dist[id]>0) {
    //[CO20190520 - printing moved to PrintXray()]        double term=lambda/(2.0*dist[id]);
    //[CO20190520 - printing moved to PrintXray()]        if(term<=1) {
    //[CO20190520 - printing moved to PrintXray()]          theta=asin(term);
    //[CO20190520 - printing moved to PrintXray()]          theta=theta*360.0/TWOPI; // rad->degrees
    //[CO20190520 - printing moved to PrintXray()]          if(theta>tol) oss
    //[CO20190520 - printing moved to PrintXray()]                          <<setw(w2)<<2.0*theta<<" " // angle
    //[CO20190520 - printing moved to PrintXray()]                          <<setw(w3)<<sf[id]<<" " // sf
    //[CO20190520 - printing moved to PrintXray()]                          <<setw(w1)<<ii0<<" " // h
    //[CO20190520 - printing moved to PrintXray()]                          <<setw(w1)<<ii1<<" " // k
    //[CO20190520 - printing moved to PrintXray()]                          <<setw(w1)<<ii2<<" " // l
    //[CO20190520 - printing moved to PrintXray()]                          <<setw(w2)<<dist[id]<<" " // dist
    //[CO20190520 - printing moved to PrintXray()]                          << "SINGLE"
    //[CO20190520 - printing moved to PrintXray()]                          << endl;
    //[CO20190520 - printing moved to PrintXray()]        } // if term<=1
    //[CO20190520 - printing moved to PrintXray()]      } // if dist>0
    //[CO20190520 - printing moved to PrintXray()]    } // i2
    //[CO20190520 - printing moved to PrintXray()]  } // i1
    //[CO20190520 - printing moved to PrintXray()]} // i0

    // Now group everything at the same distance together and only store one entry for each distance.
    // Choose hkl such that h is the largest, k the second, and l the third.
    // Multiply the sf by the number of degenerate points at that distance.
    // Get max integrated intensity for normalizations and percentages.
    //[CO20190620 - not needed here]double intmax=1E-8;
    double odist=dist[(int)ids[0][1]]; // Initialize odsit to first distance.
    double osf=sf[(int)ids[0][1]]; // Initialize osf to first distance.
    vector<vector<int> > hkl_list;
    vector<int> hkl(3);
    //aurostd::matrix<double> data;  //CO20200404 pflow::matrix()->aurostd::matrix()
    for(int i0=-kmx;i0<=kmx;i0++) {
      for(int i1=-kmx;i1<=kmx;i1++) {
        for(int i2=-kmx;i2<=kmx;i2++) {
          int ii0=i0+kmx;
          int ii1=i1+kmx;
          int ii2=i2+kmx;
          int id1=ii2+ii1*len+ii0*len*len;
          int id=(int) ids[id1][1];
          ii0=(int) ids[id1][2];
          ii1=(int) ids[id1][3];
          ii2=(int) ids[id1][4];
          double pdist=dist[id]; // Get present distance.
          // Create vector of all the hkl values with the same distance
          if(aurostd::abs(pdist-odist)<tol) { // Add present h,k,l to hkl_list.
            hkl[0]=ii0;		    
            hkl[1]=ii1;		    
            hkl[2]=ii2;		    
            hkl_list.push_back(hkl);
          }
          else { // Store one hkl, dist, total sf, multiplicity of point in data vector and then reset hkl_list vector to new hkl.
            vector<double> datav(6);
            // Sort hkl
            //[CO20190629 - waste of a class]sort(hkl_list.begin(),hkl_list.end(),hkl_cmp());
            sort(hkl_list.rbegin(),hkl_list.rend(),aurostd::compareVecElements<int>);  //CO20190629 - note that it is in descending order by hkl (greater go first)
            datav[0]=(double) hkl_list[0][0];  
            datav[1]=(double) hkl_list[0][1];  
            datav[2]=(double) hkl_list[0][2];  
            datav[3]=odist;
            datav[4]=osf*hkl_list.size();
            datav[5]=hkl_list.size();
            data.push_back(datav);
            //[CO20190620 - not needed here]if(datav[4]>intmax) {intmax=datav[4];}
            vector<int> v(0);
            hkl_list= vector<vector<int> > (0,v);
            hkl[0]=ii0;		    
            hkl[1]=ii1;		    
            hkl[2]=ii2;		    
            hkl_list.push_back(hkl);
          }
          odist=pdist;
          osf=sf[id];
        } // i2
      } // i1
    } // i0

    //[CO20190520 - printing moved to PrintXray()]// Output grouped data
    //[CO20190520 - printing moved to PrintXray()]oss << "******************** Grouped data ********************" << endl;
    //[CO20190520 - printing moved to PrintXray()]oss << "2*theta      IntIntensity         %ofMaxInt    h    k    l    dist         mult. correction    keyword " << endl;
    //[CO20190520 - printing moved to PrintXray()]for(uint i=0;i<data.size();i++) {
    //[CO20190520 - printing moved to PrintXray()]  double theta=0;
    //[CO20190520 - printing moved to PrintXray()]  if(data[i][3]>0) {
    //[CO20190520 - printing moved to PrintXray()]    double term=lambda/(2.0*data[i][3]);
    //[CO20190520 - printing moved to PrintXray()]    if(term<=1) {
    //[CO20190520 - printing moved to PrintXray()]      theta=asin(term);
    //[CO20190520 - printing moved to PrintXray()]      theta=theta*360.0/TWOPI; // rad->degrees
    //[CO20190520 - printing moved to PrintXray()]      if(theta>tol) oss
    //[CO20190520 - printing moved to PrintXray()]                      <<setw(w2)<<2.0*theta<<" " // angle
    //[CO20190520 - printing moved to PrintXray()]                      <<setw(w3)<<setprecision(2)<<data[i][4]<<setprecision(PREC_DEFAULT)<<" " // sf
    //[CO20190520 - printing moved to PrintXray()]                      <<setw(w2)<<setprecision(2)<<100*data[i][4]/intmax<<setprecision(PREC_DEFAULT)<<" " // % max sf
    //[CO20190520 - printing moved to PrintXray()]                      <<setw(w1)<<(int)data[i][0]<<" " // h
    //[CO20190520 - printing moved to PrintXray()]                      <<setw(w1)<<(int)data[i][1]<<" " // k
    //[CO20190520 - printing moved to PrintXray()]                      <<setw(w1)<<(int)data[i][2]<<" " // l
    //[CO20190520 - printing moved to PrintXray()]                      <<setw(w2)<<data[i][3]<<" " // dist
    //[CO20190520 - printing moved to PrintXray()]                      <<setw(5)<<(int)data[i][5]<<" " // mult.
    //[CO20190520 - printing moved to PrintXray()]                      <<setw(w2)<<CorrectionFactor(theta*TWOPI/360.0)<<" " // correction.
    //[CO20190520 - printing moved to PrintXray()]                      << " GROUP"
    //[CO20190520 - printing moved to PrintXray()]                      << endl;
    //[CO20190520 - printing moved to PrintXray()]    } // if term<=1
    //[CO20190520 - printing moved to PrintXray()]  } // if dist>0
    //[CO20190520 - printing moved to PrintXray()]} // for

    //[CO20190520 - printing moved to PrintXray()]// Output data to plot
    //[CO20190520 - printing moved to PrintXray()]oss << "******************** To Plot data ********************" << endl;
    //[CO20190520 - printing moved to PrintXray()]oss << "2*theta      Amplitude    keyword " << endl;
    //[CO20190520 - printing moved to PrintXray()]for(uint i=0;i<data.size();i++) {
    //[CO20190520 - printing moved to PrintXray()]  double theta=0;
    //[CO20190520 - printing moved to PrintXray()]  if(data[i][3]>0) {
    //[CO20190520 - printing moved to PrintXray()]    double term=lambda/(2.0*data[i][3]);
    //[CO20190520 - printing moved to PrintXray()]    if(term<=1) {
    //[CO20190520 - printing moved to PrintXray()]      theta=asin(term);
    //[CO20190520 - printing moved to PrintXray()]      theta=theta*360.0/TWOPI; // rad->degrees
    //[CO20190520 - printing moved to PrintXray()]      if(theta>tol) {
    //[CO20190520 - printing moved to PrintXray()]        // initial 0.
    //[CO20190520 - printing moved to PrintXray()]        oss <<setw(w2)<<2.0*theta<<" " // angle
    //[CO20190520 - printing moved to PrintXray()]            <<setw(w2)<<"0"<<" " // sf
    //[CO20190520 - printing moved to PrintXray()]            << "TOPLOT"
    //[CO20190520 - printing moved to PrintXray()]            << endl;
    //[CO20190520 - printing moved to PrintXray()]        // true value of sf/intmax.
    //[CO20190520 - printing moved to PrintXray()]        oss <<setw(w2)<<2.0*theta<<" " // angle
    //[CO20190520 - printing moved to PrintXray()]            <<setw(w2)<<setprecision(2)<<100*data[i][4]/intmax<<setprecision(PREC_DEFAULT)<<" " // sf
    //[CO20190520 - printing moved to PrintXray()]            << "TOPLOT"
    //[CO20190520 - printing moved to PrintXray()]            << endl;
    //[CO20190520 - printing moved to PrintXray()]        // final 0.
    //[CO20190520 - printing moved to PrintXray()]        oss <<setw(w2)<<2.0*theta<<" " // angle
    //[CO20190520 - printing moved to PrintXray()]            <<setw(w2)<<"0"<<" " // sf
    //[CO20190520 - printing moved to PrintXray()]            << "TOPLOT"
    //[CO20190520 - printing moved to PrintXray()]            << endl;
    //[CO20190520 - printing moved to PrintXray()]        // tpx
    //[CO20190520 - printing moved to PrintXray()]      } // if theta<tol
    //[CO20190520 - printing moved to PrintXray()]    } // if term<=1
    //[CO20190520 - printing moved to PrintXray()]  } // if dist>0
    //[CO20190520 - printing moved to PrintXray()]} // for
  } // end routine
} // namespace pflow

// ***************************************************************************
// GetXRAY
// ***************************************************************************
// This function gets XRAY following the convasp framework.
// Dane Morgan
// updated by Corey Oses 20190520
namespace pflow {
  void GetXray(const xstructure& str, vector<double>& dist,vector<double>& sf,
      vector<double>& scatt_fact, //CO20190520
      vector<double>& mass, vector<double>& twoB_vec,double lambda) {
    stringstream message; //CO20190322
    // Get data from str.
    // Set scale to 1 so you don't need to rescale coordinates.
    xstructure sstr=str;
    sstr=ReScale(sstr,1.0);

    xmatrix<double> rlat(3,3);rlat=sstr.klattice;
    int num_atoms=sstr.atoms.size();
    //  aurostd::matrix<double> fpos=GetFpos(sstr);  //CO20200404 pflow::matrix()->aurostd::matrix()

    // Set parameters for Debye Waller factors.
    double temp=300;
    double debye_temp=300;
    double dw,theta;

    // Get scattering factors and masses and Debye-Waller 2B values.
    double h=PLANCKSCONSTANT_h;  //ME20181020
    for(int i=0;i<num_atoms;i++) {
      //_atom atom=sstr.atoms.at(i); atom.CleanName();  //CO20190322 - why make another copy of atom if we already made copy of xstructure??
      _atom& atom=sstr.atoms[i]; atom.CleanName();
      if(atom.name_is_given==FALSE || atom.cleanname.empty()){  //CO20190322
        message << "Need to provide atom names";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
      }
      scatt_fact[i]=GetXrayScattFactor(atom.cleanname,lambda);
      mass[i]=GetAtomMass(atom.cleanname);
      twoB_vec[i]=h*h*temp*12.0/(mass[i]*KBOLTZ*debye_temp*debye_temp);
    }

    //  Get max h,k,l value that gives any scattering (theta<=90).
    int kmx0=(int) (4.0*PI/(modulus(rlat(1))*lambda)+1);
    int kmx1=(int) (4.0*PI/(modulus(rlat(2))*lambda)+1);
    int kmx2=(int) (4.0*PI/(modulus(rlat(3))*lambda)+1);
    int kmx=max(kmx0,kmx1);kmx=max(kmx2,kmx);
    int len=2*kmx+1;  //simply -kmx to kmx in for loop
    int tlen=len*len*len;  //all combos of hkl, convert 3D HKL search to 1D
    double sfr; //real part of structure factor
    double sfi; //imaginary part of structure factor
    dist = vector<double> (tlen,0.0);
    sf = vector<double> (tlen,0.0);

    //initialization
    int ii0=0;
    int ii1=0;
    int ii2=0;
    int id=0;
    xvector<double> rv(3);
    double rvnorm=0.0;
    double term=0.0;
    double gdotr;

    for(int i0=-kmx;i0<=kmx;i0++) {
      for(int i1=-kmx;i1<=kmx;i1++) {
        for(int i2=-kmx;i2<=kmx;i2++) {
          ii0=i0+kmx;
          ii1=i1+kmx;
          ii2=i2+kmx;
          id=ii2+ii1*len+ii0*len*len; //this is index that converts 3D HKL search to 1D
          sfr=0.0;
          sfi=0.0;
          for(int ic=1;ic<=3;ic++){rv(ic)=i0*rlat(1,ic)+i1*rlat(2,ic)+i2*rlat(3,ic);}
          rvnorm=modulus(rv);
          if(rvnorm>0.0) dist[id]=TWOPI/rvnorm;
          for(int iat=0;iat<num_atoms;iat++) {
            // Get Debye Waller factor.
            dw=1;
            if(dist[id]>0) {
              term=lambda/(2.0*dist[id]);
              if(term<=1) { //angle, must be less than equal to 1 or domain error occurs
                theta=std::asin(term);  //radians
                theta*=360.0/TWOPI; // rad->degrees
                dw=DebyeWallerFactor(theta,temp,debye_temp,mass[iat],lambda*1.0E-10);
              }
            }
            const _atom& atom=sstr.atoms[iat];
            //  double gdotr=i0*fpos[iat][0]+i1*fpos[iat][1]+i2*fpos[iat][2];
            gdotr=i0*atom.fpos(1)+i1*atom.fpos(2)+i2*atom.fpos(3);
            gdotr*=TWOPI;
            sfr+=dw*scatt_fact[iat]*std::cos(gdotr);
            sfi-=dw*scatt_fact[iat]*std::sin(gdotr);
          }
          sf[id]=(sfr*sfr)+(sfi*sfi); //sum the squares, vector addition
        } // i2
      } // i1
    } // i0
  }
}

// ***************************************************************************
// RDF RDF RDF RDF RDF RDF RDF RDF RDF RDF RDF RDF RDF RDF RDF RDF RDF RDF RDF
// ***************************************************************************
// Function GetRDF
// ***************************************************************************
namespace pflow {
  vector<vector<int> > getvvitypesRDF(const xstructure& xstr){  //CO20220627
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    vector<vector<int> > _vvitypes,vvitypes;
    vector<int> vitypes,vsizes;
    uint i=0,j=0;
    for(i=0;i<2;i++){
      _vvitypes.push_back(vector<int>(0));
      for(j=0;j<xstr.num_each_type.size();j++){_vvitypes.back().push_back(j);}  //initial set, not all of these are viable options
      vsizes.push_back(xstr.num_each_type.size());
    }
    aurostd::xcombos xc(vsizes,'E');
    while(xc.increment()){
      vitypes=xc.applyCombo(_vvitypes);
      if(vitypes[1]<vitypes[0]){continue;}  //filter so we don't get duplicates
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " vitypes=" << aurostd::joinWDelimiter(vitypes,",") << endl;}
      vvitypes.push_back(vitypes);
    }
    return vvitypes;
  }
  void GetRDF(const xstructure& xstr_in,aurostd::xmatrix<double>& rdf_all,const double rmax,const int nbins,bool raw_counts,const double sigma,const int _window_gaussian) { //CO20220624 - new procedure
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}
    xstructure xstr(xstr_in);
    if(xstr.species.size()==0){xstr.species=aurostd::vector2deque(pflow::getFakeElements(xstr.num_each_type.size()));}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " species=" << aurostd::joinWDelimiter(xstr.species,",") << endl;}
    deque<_atom> atoms_cell;
    deque<deque<uint> > i_neighbors;
    deque<deque<double> > distances;
    double rmin=0;
    //SC+CO20220711 - unique_only can only be true if we scale properly by the Wyckoff multiplicity
    //easier to brute force count
    bool prim=true,unique_only=false;
    xstr.GetNeighbors(atoms_cell,i_neighbors,distances,rmax,rmin,prim,unique_only);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " atoms_cell.size()=" << atoms_cell.size() << endl;}
    //for compound ABC, need counts for AA, AB, AC, BB, BC, CC
    //this is bit-enumeration for sets of two (i,j) ensuring i<j
    vector<vector<int> > vvitypes=getvvitypesRDF(xstr);
    rdf_all=aurostd::xmatrix<double>(nbins,vvitypes.size()+1); //resize, +1 because last column will be totals
    uint i=0;
    //SC+CO20220711
    //if the unit cell has more than one atom of each type, we need to average their contributions
    //this is equivalent to adding partial contribution based on num_each_type
    vector<double> vincrements;
    for(i=0;i<atoms_cell.size();i++){vincrements.push_back( 1.0/xstr.num_each_type[ atoms_cell[i].type ] );}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " vincrements=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vincrements),",") << endl;}
    int itype=0,ibin=0;
    double dist=0.0,rad=0.0;
    double drad=rmax/(double)nbins; //sphere shell volume=4*pi*rad^2*drad
    vector<int> vitypes;
    uint j=0,k=0;
    for(i=0;i<i_neighbors.size();i++){
      for(j=0;j<i_neighbors[i].size();j++){
        dist=distances[i][j];
        if(dist<rmax){
          ibin=int((dist/rmax)*nbins)+rdf_all.lrows;
          //need to figure out which bin to increment
          //first create vitypes of two atoms, sort just in case
          vitypes.clear();
          vitypes.push_back(atoms_cell[i].type);vitypes.push_back(xstr.grid_atoms[i_neighbors[i][j]].type);
          std::sort(vitypes.begin(),vitypes.end());
          //search vvitypes for match
          itype=AUROSTD_MAX_INT;
          for(k=0;k<vvitypes.size()&&itype==AUROSTD_MAX_INT;k++){
            if(vitypes==vvitypes[k]){itype=k;}
          }
          if(itype==AUROSTD_MAX_INT){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"could not find itype",_INDEX_MISMATCH_);}
          itype+=rdf_all.lcols; //+1 because xmatrix starts with 1
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << " ibin=" << ibin << " (dist=" << dist << ")" << endl;
            cerr << __AFLOW_FUNC__ << " itype=" << itype << endl;
          }
          rdf_all[ibin][itype]+=vincrements[i];
          rdf_all[ibin][rdf_all.ucols]+=vincrements[i];  //last row is totals
        }
      }
    }
    if(raw_counts==false){
      for(ibin=rdf_all.lrows;ibin<=rdf_all.urows;ibin++){
        rad=drad*(ibin-rdf_all.lrows);
        if(ibin==rdf_all.lrows && !aurostd::iszero(rdf_all[ibin][rdf_all.ucols],_ZERO_TOL_)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"first bin cannot be normalized (rad==0), increase nbins",_VALUE_ILLEGAL_);}
        for(itype=rdf_all.lcols;itype<=rdf_all.ucols;itype++){
          if(aurostd::nonZeroWithinTol(rad,_ZERO_TOL_)){rdf_all[ibin][itype]/=(4.0*PI*std::pow(rad,2.0)*drad);}
        }
      }
    }
    if(!aurostd::iszero(sigma)){  //smooth
      aurostd::xvector<double> _xcol,xcol,filter_gaussian;
      int window_gaussian=_window_gaussian;
      if(window_gaussian==0){window_gaussian=aurostd::gaussian_filter_get_window(sigma);}
      filter_gaussian=aurostd::gaussian_filter_xv(sigma,window_gaussian,rdf_all.lrows);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " filter_gaussian[sigma=" << sigma << ",window_gaussian=" << window_gaussian << "]=" << filter_gaussian << endl;}
      for(itype=rdf_all.lcols;itype<=rdf_all.ucols;itype++){
        _xcol=rdf_all.getcol(itype);  //extract column
        xcol=aurostd::convolution(_xcol,filter_gaussian,CONV_SHAPE_SAME);
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " xcol(pre)[itype=" << itype << "]=" << _xcol << endl;
          cerr << __AFLOW_FUNC__ << " xcol(post)[itype=" << itype << "]=" << xcol << endl;
        }
        for(ibin=rdf_all.lrows;ibin<=rdf_all.urows;ibin++){rdf_all[ibin][itype]=xcol[ibin];}  //plug back in
      }
    }
    if(LDEBUG){PrintRDF(xstr,rmax,nbins,rdf_all,cerr);}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " END" << endl;}
  }
// This function gets the radial distribution functions (RDF).
// The RDF for an atom is just a binned histogram of the
// number of atoms (possibly of a given type) a distance r away.  
// The RDF for a type is just the average RDF's for all the atoms
// of that type.
// The rdf is stored as follows: I=(0,natoms-1),J,K=(0,ntypes)
// row: J+(ntypes+1)*I = atom I / Type J RDF.
//    (if row: J=ntypes then the rdf is the atom I / All types RDF.)
// I don't think the following ever got coded.
// lines (ntypes+1)*natoms+K+(ntypes+1)*J = type J / type K average RDF.
//    (if row: K=ntypes then the rdf is the type I / All types RDF.)
// Dane Morgan, Modified by Stefano Curtarolo
// CO20220627 - not working, tethered to a few old routines, will need to be rewritten (later)
// use new procedure above GetRDF()
  void GetRDF_20220101(const xstructure& str, const double rmax,
      const int nbins, aurostd::matrix<double>& rdf_all) {  //CO20200404 pflow::matrix()->aurostd::matrix()
    int natoms=str.atoms.size();
    std::deque<int> num_each_type=str.num_each_type;
    int ntyp=num_each_type.size();
    rdf_all=aurostd::matrix<double> ((ntyp+1)*natoms,nbins); //CO20200404 pflow::matrix()->aurostd::matrix()
    deque<deque<_atom> > nmat;
    // [OBSOLETE]    pflow::GetNeighData(str,rmax,nmat);  //CO20220623 - using new GetNeighData()
    str.GetNeighData(rmax,nmat);   // once GetRD goes in xstructure I can remove the copy //CO20220623 - using new GetNeighData()
    // for(int i=0;i<nmat[0].size();i++) cout << AtomDist(nmat[0][0],nmat[0][i]) << " "; cout << endl;
    for(int I1=0;I1<(int)nmat.size();I1++) { // Each atom for which we find RDF.
      int I2=1;
      double dist=0;
      int s2=(int)nmat.at(I1).size();
      _atom a1=nmat.at(I1).at(0);
      // Check each atom neighbor for I1 in proper range.
      dist=AtomDist(a1,nmat.at(I1).at(I2));
      while (dist<rmax && I2<s2) {
        int ib=int((dist/rmax)*nbins); // first bin is [0,rmax/nbins).
        int J2=nmat.at(I1).at(I2).type;     // CONVASP_MODE
        //      cout << I1 << " " << I2 << " " <<  J2 << " " << s2 << " " << dist << " " << ib << endl;
        rdf_all[(ntyp+1)*I1+J2][ib]++; // Does all binning for atom/type RDFs.
        I2++;
        if(I2<s2) dist=AtomDist(a1,nmat[I1][I2]);
      }  
      // Get sum over all types
      for(int ib=0;ib<nbins;ib++) {
        for(int it=0;it<ntyp;it++) {
          rdf_all[(ntyp+1)*I1+ntyp][ib]+=rdf_all[(ntyp+1)*I1+it][ib];
        }
      }
    } // I1 loop
  }
}
// ***************************************************************************
// Function GetRDFShells
// ***************************************************************************
//   This function gets the neighbor shells from the
//   radial distribution functions.  Shells are found
//   by the following method.  The code calculates
//   the change in rdf from drdf=rdf[i+1]-rdf[i].  We
//   step through drdf starting adding up rdf into shell1
//   while drdf>=0 and then while drdf<=0.  When drdf gets
//   >=0 again we start a new shell until drdf passes
//   through <=0 and is >=0 again, at which point we again
//   start a new shell.  In other words, we increment the
//   shell we are adding to when drdf goes from >=0 to <0 to
//  >=0 again.
namespace pflow {
  void GetRDFShells(const xstructure& str,const double& rmax,const int& nbins,
      const int& smooth_width,const aurostd::matrix<double>& rdf,  //CO20200404 pflow::matrix()->aurostd::matrix()
      aurostd::matrix<double>& rdfsh,aurostd::matrix<double>& rdfsh_loc) { //CO20200404 pflow::matrix()->aurostd::matrix()
    // double TOL=1e-5; //DM not used
    // int natom=(int)str.atoms.size(); //DM not used
    if(smooth_width) {;} // phony just to keep smooth_width busy
    if(str.atoms.size()) {;} // phony just to keep str busy

    int _rdi=(int)rdf.size();
    rdfsh = aurostd::matrix<double>  (_rdi,0); pflow::VVset(rdfsh,0);  //CO20200404 pflow::matrix()->aurostd::matrix()
    rdfsh_loc =  aurostd::matrix<double> (_rdi,0);  pflow::VVset(rdfsh_loc,0); //CO20200404 pflow::matrix()->aurostd::matrix()
    double dr=(rmax/(double)nbins);
    for(int i=0;i<_rdi;i++) {

      /*
      // get smoothed rdf_sm
      vector<double> rdf_sm(nbins,0.0);
      int cnt=0;
      for(int ib=0;ib<nbins;ib++) {
      for(int ism=-(smooth_width-1)/2;ism<=(smooth_width-1)/2;ism++) {
      int id=ib+ism;
      if(id>0 && id<nbins) {
      cnt++;
      rdf_sm[ib]+=rdf[i][ib+ism];
      }
      }
      if(cnt>0) rdf_sm[ib]=rdf_sm[ib]/cnt;
      cnt=0;
      }
      rdf_sm=SmoothFunc(rdf[i],smooth_width);
      */

      // get drdf,ddrdf of rdf
      vector<double> drdf(nbins,0.0);
      vector<double> ddrdf(nbins,0.0);
      for(int ib=1;ib<nbins-1;ib++) {
        drdf[ib]=(rdf[i][ib+1]-rdf[i][ib-1]);
      }
      for(int ib=1;ib<nbins-1;ib++) {
        //      ddrdf[ib]=(rdf[i][ib+1]+rdf[i][ib-1]-2*rdf[i][ib]);
        ddrdf[ib]=(drdf[ib+1]-drdf[ib-1]);
      }

      // get all zeros of drdf with ddrf!=0
      // -999 means nothing,-1 means ddrdf<0, 0 means ddrdf=0, +1 means ddrdf>0.
      vector<int> drdf_zeros(nbins-1,-999);
      for(int ib=1;ib<nbins-1;ib++) {
        // If drdf==0
        if(drdf[ib]==0 && ddrdf[ib]<0) drdf_zeros[ib]=-1;
        if(drdf[ib]==0 && ddrdf[ib]==0) drdf_zeros[ib]=0;
        if(drdf[ib]==0 && ddrdf[ib]==1) drdf_zeros[ib]=1;
        // If drdf== is passing through 0 from +->- or -->+
        if(drdf[ib]!=0) {
          if(drdf[ib-1]>0 && drdf[ib+1]<0) drdf_zeros[ib]=-1;
          if(drdf[ib-1]<0 && drdf[ib+1]>0) drdf_zeros[ib]=+1;
          // If drdf== is becoming or changinf from 0 through
          // 0->+/- or +/-->0.
          if(drdf[ib-1]==0 && ddrdf[ib]<0) drdf_zeros[ib]=-1;
          if(drdf[ib-1]==0 && ddrdf[ib]==0) drdf_zeros[ib]=0;
          if(drdf[ib-1]==0 && ddrdf[ib]>0) drdf_zeros[ib]=1;
          if(drdf[ib+1]==0 && ddrdf[ib]<0) drdf_zeros[ib]=-1;
          if(drdf[ib+1]==0 && ddrdf[ib]==0) drdf_zeros[ib]=0;
          if(drdf[ib+1]==0 && ddrdf[ib]>0) drdf_zeros[ib]=1;
        }

        // tpx
        //      cout << "drdf_zeros " << ib << " " << rdf[i][ib] << " " << drdf[ib] << " " << ddrdf[ib] << " " << drdf_zeros[ib] << endl;
      }

      // get the actual shell atoms counts and avg radius.
      double shell=0;
      double rsh_avg=0;
      double state_dn=0;
      for(int ib=0;ib<nbins-2;ib++) {
        //tpx
        //           cout << i << " " << ib << " "
        //	         << rdf[i][ib] << " " << drdf[ib] << " " << shell << " " << endl;
        shell=shell+rdf[i][ib];
        double rad=(dr*(double)ib)+dr/2.0;
        rsh_avg=rsh_avg+rdf[i][ib]*rad;
        if(drdf_zeros[ib]==-1) { // Set state_dn
          state_dn=1;
        }
        if(state_dn && drdf_zeros[ib]==1) { // New shell
          if(shell>0) rsh_avg=rsh_avg/shell;
          rdfsh[i].push_back(shell);
          rdfsh_loc[i].push_back(rsh_avg);
          shell=0;
          rsh_avg=0;
          state_dn=0;
        }
      } // for ib
    } // for i
  }
}

// ***************************************************************************
// Function RdfSh_RMS
// ***************************************************************************
// This function compares the rdf shells of two
// atoms for each type and returns the rms.  
namespace pflow {
  double RdfSh_RMS(const int iaA, const int iaB, const int nsh_max,const int nt,
      const aurostd::matrix<double>& rdfsh_all_A
      ,const aurostd::matrix<double>& rdfsh_all_B) { //CO20200404 pflow::matrix()->aurostd::matrix()
    double rms=0;
    int cnt=0;
    for(int it=0;it<nt;it++) {
      int idA=(nt+1)*iaA+it;
      int idB=(nt+1)*iaB+it;
      int tempsh=min((int) rdfsh_all_A[idA].size(),(int) rdfsh_all_B[idB].size());
      int nsh = min(tempsh,nsh_max);
      for(int ish=0;ish<nsh;ish++) {
        cnt++;
        rms+=(rdfsh_all_A[idA][ish]-rdfsh_all_B[idB][ish])
          *(rdfsh_all_A[idA][ish]-rdfsh_all_B[idB][ish]);
        // tpx
        //      cout << "iaA,iaB,it,ish,cnt,rms " <<iaA<<" "<<iaB<<" "<<it<<" "<<ish<<" "<< cnt << " " << rms << endl;
      }
    }
    if(cnt>0) rms=sqrt(rms/cnt);
    return rms;
  }
}

// ***************************************************************************
// Function CmpRDFShell
// ***************************************************************************
// This function compares the rdf shells of two structures.
// Two atoms similarity are based on the rms error between
// their first nsh nn shells (for every type).  Atoms must
// be the same type to even be compared.
// The results are given in:
// best_match: For each atom in str_A gives the best matching atom
//             in str_B of the smae type.
// rms_mat:  Gives the rms error between every pair of atoms.
// For the best matches I start with the first atom of str_A and
// compare to all of B.  Each successive match for an atom of str_A
// is performed excluding the previously matched atoms of str_B.
// str_A and str_B must have the same number of each type of atom.
namespace pflow {
  void CmpRDFShells(const xstructure& str_A, const xstructure& str_B,
      const aurostd::matrix<double>& rdfsh_all_A,
      const aurostd::matrix<double>& rdfsh_all_B,
      const int nsh, vector<int>& best_match,
      aurostd::matrix<double>& rms_mat) {  //CO20200404 pflow::matrix()->aurostd::matrix()
    double TOL=1e-15;
    std::deque<int> netype_A=str_A.num_each_type;
    std::deque<int> netype_B=str_B.num_each_type;
    // Exit if A and B have different numbers of any types of atoms.
    if(!pflow::VVequal(netype_A,netype_B) || netype_A.size()!=netype_B.size()) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"structures A and B do not have the same number of each type of atom",_INPUT_ILLEGAL_); //CO20200624
    }
    int nat=str_A.atoms.size();
    int nt=netype_A.size();

    // Get rms_mat
    rms_mat = aurostd::matrix<double> (nat,nat); pflow::VVset(rms_mat,-1.0); //CO20200404 pflow::matrix()->aurostd::matrix()

    for(int iaA=0;iaA<nat;iaA++) {
      for(int iaB=0;iaB<nat;iaB++) {
        rms_mat[iaA][iaB]=RdfSh_RMS(iaA,iaB,nsh,nt,rdfsh_all_A,rdfsh_all_B);
      }
    }
    // Get best_match matrix.
    // For each atom in A check the rms error for all atoms of that
    // type in B.  Then take the best one and then remove it from
    // the list for further comparisons.
    best_match=vector<int> (nat,-1);
    vector<int> cand_Batoms(nat,0);
    for(int i=0;i<nat;i++) cand_Batoms[i]=i;
    for(int iaA=0;iaA<(int)rms_mat.size();iaA++) {
      double rms=-10;
      double trms=-10;
      int best_at_id=-1;
      int typeA=str_A.atoms.at(iaA).type;
      // Check each B atom and if it is the right type and lowest rms
      // then track its id.
      for(int ipB=0;ipB<(int)cand_Batoms.size();ipB++) {
        int iaB=cand_Batoms[ipB];
        int typeB=str_B.atoms.at(iaB).type;
        if(typeA==typeB) {
          trms=rms_mat[iaA][iaB];
          if(rms<-1) rms=trms+1.0;
          if(trms<rms-TOL) {
            best_at_id=iaB;
            rms=trms;
          }
        }
      }
      best_match[iaA]=best_at_id;
      // Remove best_at_id from future comparisons.
      cand_Batoms.erase(find(cand_Batoms.begin(),cand_Batoms.end(),best_at_id));
    }
  }
}

// ***************************************************************************
// Function GetSmoothRDF
// ***************************************************************************
// This function gets a smoothed RDF.
namespace pflow {
  aurostd::matrix<double> GetSmoothRDF(const aurostd::matrix<double>& rdf, //CO20200404 pflow::matrix()->aurostd::matrix()
      const double& sigma) {
    aurostd::matrix<double> rdf_sm=rdf;  //CO20200404 pflow::matrix()->aurostd::matrix()
    for(int i=0;i<(int)rdf.size();i++) {
      rdf_sm[i]=pflow::SmoothFunc(rdf[i],sigma);
    }
    return rdf_sm;
  }
}
// ***************************************************************************
// COMPARE STRUCTURES COMPARE STRUCTURES COMPARE STRUCTURES COMPARE STRUCTURES
// ***************************************************************************
// ***************************************************************************
// CmpStrDist
// ***************************************************************************
// Compares the distances within rcut between str1 and str2.
// Assigns matrix of distances.  Each row is for a different
// pair type. Pair types are identified by pair i,j in row
// j+ntypes*i.  1<= i,j <=ntypes. I always take j>=i.
// Original by Dane Morgan, modified by STefano Curtarolo for type shift !
namespace pflow {
  void CmpStrDist(xstructure str1, xstructure str2,const double& cutoff,
      aurostd::matrix<double>& dist1, aurostd::matrix<double>& dist2,
      aurostd::matrix<double>& dist_diff,aurostd::matrix<double>& dist_diff_n) { //CO20200404 pflow::matrix()->aurostd::matrix()

    deque<deque<_atom> > neigh_mat1;
    deque<deque<_atom> > neigh_mat2;
    // [OBSOLETE] pflow::GetNeighData(str1,cutoff,neigh_mat1);  //CO20220623 - using new GetNeighData()
    // [OBSOLETE] pflow::GetNeighData(str2,cutoff,neigh_mat2);  //CO20220623 - using new GetNeighData()
    str1.GetNeighData(cutoff,neigh_mat1); //CO20220623 - using new GetNeighData()
    str2.GetNeighData(cutoff,neigh_mat2); //CO20220623 - using new GetNeighData()

    int ntypes1=(int)str1.num_each_type.size();
    int ntypes2=(int)str2.num_each_type.size();
    int nprs1=(ntypes1*(ntypes1+1))/2;
    int nprs2=(ntypes2*(ntypes2+1))/2;
    vector<double> tmpvec;
    dist1=aurostd::matrix<double> (nprs1,tmpvec);  //CO20200404 pflow::matrix()->aurostd::matrix()
    dist2=aurostd::matrix<double> (nprs2,tmpvec);  //CO20200404 pflow::matrix()->aurostd::matrix()

    // Collect distances.
    // dist1
    for(int ia=0;ia<(int)neigh_mat1.size();ia++) {
      _atom a = neigh_mat1[ia][0];
      for(int in=1;in<(int)neigh_mat1[ia].size();in++) {
        _atom an = neigh_mat1[ia][in];
        int ta=a.type;              // CONVASP_MODE
        int tna=an.type;            // CONVASP_MODE
        int i=std::min(ta,tna);
        int j=std::max(ta,tna);
        int id=j-i+ntypes1*i-max(0,i*(i-1)/2);
        dist1[id].push_back(AtomDist(a,an));
      } // in
    } // ia
    // dist2
    for(int ia=0;ia<(int)neigh_mat2.size();ia++) {
      _atom a = neigh_mat2[ia][0];
      for(int in=1;in<(int)neigh_mat2[ia].size();in++) {
        _atom an = neigh_mat2[ia][in];
        int ta=a.type;            // CONVASP_MODE
        int tna=an.type;          // CONVASP_MODE
        int i=std::min(ta,tna);
        int j=std::max(ta,tna);
        int id=j-i+ntypes2*i-max(0,i*(i-1)/2);
        dist2[id].push_back(AtomDist(a,an));
      } // in
    } // ia

    // Sort dist vectors.
    for(int ip=0;ip<nprs1;ip++) {
      sort(dist1[ip].begin(),dist1[ip].end());  
    }
    for(int ip=0;ip<nprs2;ip++) {
      sort(dist2[ip].begin(),dist2[ip].end());  
    }

    // tpx
    //  pflow::Vout(dist1[0],cout);
    //  pflow::Vout(dist2[0],cout);
    // Assign distance difference matrix
    int ntypes_min=std::min(ntypes1,ntypes2);
    int nprs_min=ntypes_min*(ntypes_min+1)/2;
    dist_diff = aurostd::matrix<double> (nprs_min,tmpvec); //CO20200404 pflow::matrix()->aurostd::matrix()
    dist_diff_n = aurostd::matrix<double> (nprs_min,tmpvec); //CO20200404 pflow::matrix()->aurostd::matrix()
    for(int i=0;i<ntypes_min;i++) {
      for(int j=i;j<ntypes_min;j++) {
        int id1=j-i+ntypes1*i-max(0,i*(i-1)/2);
        int id2=j-i+ntypes2*i-max(0,i*(i-1)/2);
        int idmin=j-i+ntypes_min*i-max(0,i*(i-1)/2);
        int num_dist=std::min((int)dist1[id1].size(),(int)dist2[id2].size());
        for(int ip=0;ip<num_dist;ip++) {
          double d1=dist1[id1][ip];
          double d2=dist2[id2][ip];
          dist_diff[idmin].push_back(d2-d1);
          dist_diff_n[idmin].push_back(2*(d2-d1)/(d1+d2));
        }// ip
      }// j
    }// i
  } // end routine
}

// ****************************************************************************************************
// PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA P
// ****************************************************************************************************
namespace pflow {
  // Constructors
  pdosdata::pdosdata() {
  }
  void pdosdata::PrintParams(ostream& outf, const std::vector<string>& Ltotnames) {
    outf << "EMIN = " << emin << endl;
    outf << "EMAX = " << emax << endl;
    outf << "EFERMI = " << efermi << endl;
    outf << "NBINS = " << nbins << endl;
    outf << "SPIN = " << spin << endl;
    outf << "NLM = " << nlm << endl;
    outf << "NATOMS = " << natoms << endl;
    outf << "SMOOTH_SIGMA = " << smooth_sigma << endl;
    outf << "PRINT_PARAMS = " << print_params << endl;
    for(int i=0;i<(int)pdos_at.size();i++) {
      outf << "# CASE = " << i+1 << endl;
      outf << "  ATOMS = ";
      Vout(pdos_at[i],outf);
      if(pdos_k.size()>0) {
        outf << "  KPOINTS = ";
        Vout(pdos_k[i],outf);
      }
      if(pdos_b.size()>0) {
        outf << "  BANDS = ";
        Vout(pdos_b[i],outf);
      }
      outf << "  LMVALUES = ";
      Vout(pdos_lm[i],outf);
      outf << "  # ( LMVALUES = ";
      for(int j=0;j<(int)pdos_lm[i].size();j++) {
        outf << " " << Ltotnames[pdos_lm[i][j]-1];
      }
      outf << " )" << endl;
    }
  }

  void pdosdata::PrintPDOS(ostream& outf, const int& sp) {
    outf.precision(5);
    outf.setf(std::ios::fixed,std::ios::floatfield);
    outf.setf(std::ios::left,std::ios::adjustfield);

    if(sp==0) outf << "# Energy         Up            Cumulative_Up" << endl;
    if(sp==1) outf << "# Energy         Up             Dn             Up-Dn        Cumulative_Up   Cumulative_Dn   Cumulative_Up-Dn" << endl;
    for(int ib=0;ib<(int)pdos.size();ib++) {
      for(int i=0;i<(int)pdos[ib].size();i++) {
        outf << setw(10) << pdos[ib][i] << "     ";
      }
      outf << endl;
    }
  }
}

// ****************************************************************************************************
// RAY TRACING RAY TRACING RAY TRACING RAY TRACING RAY TRACING RAY TRACING RAY TRACING RAY TRACING RAY
// ****************************************************************************************************
namespace pflow {
  void rtparams::free() {
  }

  void rtparams::copy(const rtparams& b) {
    resx=b.resx;
    resy=b.resy;
    viewdir=b.viewdir;
    viewdir_s=b.viewdir_s;
    updir=b.updir;
    updir_s=b.updir_s;
    zoom=b.zoom;
    aspectratio=b.aspectratio;
    antialiasing=b.antialiasing;
    raydepth=b.raydepth;
    center=b.center;
    center_guide=b.center_guide;
    center_s=b.center_s;
    background=b.background;
    lightcenter=b.lightcenter;
    lightrad=b.lightrad;
    lightcolor=b.lightcolor;
    sphtex_tex=b.sphtex_tex;
    sphtex_tex_def=b.sphtex_tex_def;	//CO20190329
    sphtex_color=b.sphtex_color;
    sphtex_color_def=b.sphtex_color_def;
    sphtex_phong=b.sphtex_phong;
    sphtex_phong_def=b.sphtex_phong_def;
    sphtex_names=b.sphtex_names;
    sph_rad=b.sph_rad;
    sph_rad_def=b.sph_rad_def;
    shading=b.shading;
    outfile=b.outfile;
    sc=b.sc;
    sc_s=b.sc_s;
    calc_type=b.calc_type;
    input_files=b.input_files;
    first_set=b.first_set;
    insert_file=b.insert_file;
    rotation=b.rotation;
    // Plane variables
    plane=b.plane;
    plane_s=b.plane_s;
    plane_center=b.plane_center;
    plane_center_def=b.plane_center_def;
    plane_center_s=b.plane_center_s;
    plane_normal=b.plane_normal;
    plane_normal_def=b.plane_normal_def;
    plane_normal_s=b.plane_normal_s;
    plane_color=b.plane_color;
    plane_color_def=b.plane_color_def;
    plane_color_s=b.plane_color_s;
    planetex_tex=b.planetex_tex;
    planetex_tex_def=b.planetex_tex_def;
    planetex_tex_s=b.planetex_tex_s;
  }

  // Constructors
  rtparams::rtparams() {
    calc_type = 0;
    resx=600;
    resy=600;
    viewdir = vector<double> (3,0.0);
    viewdir[1]=1;
    viewdir_s = 0;
    updir = vector<double> (3,0.0);
    updir[2]=1;
    updir_s=0;
    zoom=1;
    aspectratio=1;
    antialiasing=0;
    raydepth=12;
    center = vector<double> (3,0.0);
    center_guide = vector<double> (6,0.0);
    center_s = 0;
    background = vector<double> (3,1.0);
    lightcenter = aurostd::matrix<double> (1,center);  //CO20200404 pflow::matrix()->aurostd::matrix()
    lightrad = vector<double> (1,0.001);
    vector<double> color (3,0.0);
    lightcolor = aurostd::matrix<double> (1,color);  //CO20200404 pflow::matrix()->aurostd::matrix()
    sphtex_tex_def = vector<double> (4);
    sphtex_tex_def[0] = 0.3; // ambient
    sphtex_tex_def[1] = 0.6; // diffuse
    sphtex_tex_def[2] = 0.2; // specular
    sphtex_tex_def[3] = 1.0; // opacity
    sphtex_color_def = vector<double> (3,0.5);
    sphtex_phong_def = vector<double> (2,0.0);
    sphtex_phong_def[1] = 10000;
    sph_rad_def = 1;
    shading = "fullshade";
    outfile = "POSCAR_RT";
    sc=aurostd::matrix<double> (3,3); pflow::VVset(sc,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
    sc[0][0]=1;
    sc[1][1]=1;
    sc[2][2]=1;
    sc_s=0;
    first_set=1;
    insert_file="NO_INSERT_FILE";
    rotation = vector<double> (6,0.0);
    struct_origin = vector<double> (3,0.0);
    //  int struct_origin_s=0;
    // Plane variables
    plane=0;
    plane_s=0;
    plane_center = vector<double> (3);
    plane_normal = vector<double> (3);
    plane_color = vector<double> (3);
    planetex_tex = vector<double> (4);
    plane_center_def = vector<double> (3);
    plane_normal_def = vector<double> (3);
    plane_color_def = vector<double> (3);
    planetex_tex_def = vector<double> (4);
    plane_center_s=0;
    plane_normal_s=0;
    plane_color_s=0;
    planetex_tex_s=0;
  }

  rtparams::rtparams(const rtparams& b) {copy(b);}

  const rtparams& rtparams::operator=(const rtparams& b) {
    if(this != &b) {
      free();
      copy(b);
    }
    return *this;
  }

  void rtparams::Print() const {
    cout << "CENTER = " << center_guide[0] << " " << center_guide[1] << " " << center_guide[2] << " " << center_guide[3] << " " << center_guide[4] << " " << center_guide[5] << endl;
    cout << "VIEWDIR = " << viewdir[0] << " " << viewdir[1] << " " << viewdir[2] << endl;
    cout << "UPDIR = " << updir[0] << " " << updir[1] << " " << updir[2] << endl;
    cout << "STRUCTURE_ORIGIN = " << struct_origin[0] << " " << struct_origin[1] << " " << struct_origin[2] << endl;
    cout << "ROTATION = " << rotation[0] << " " << rotation[1]
      << " " << rotation[2] << " " << rotation[3]
      << " " << rotation[4] << " " << rotation[5] << endl;
  }
}

// ***************************************************************************
// GetDatFromOUTCAR
// ***************************************************************************
// This gets the lattice vectors and num_each_type from an
// OUTCAR file.
// Dane Morgan style
namespace pflow {
  void GetDatFromOutcar(vector<aurostd::matrix<double> >& lat_vec, //CO20200404 pflow::matrix()->aurostd::matrix()
      deque<int>& num_each_type, ifstream& outcar_inf) {
    int first_lat=1;
    vector<string> a(3,"Z");
    lat_vec.clear();
    string sdum;
    while (outcar_inf >> a[2]) {
      string key= (a[0]+" "+a[1]+" "+a[2]);
      if(key=="reciprocal lattice vectors") {
        aurostd::matrix<double> lat(3,3);  //CO20200404 pflow::matrix()->aurostd::matrix()
        for(int ic=0;ic<3;ic++) {
          outcar_inf >> lat[ic][0] >> lat[ic][1] >> lat[ic][2];
          outcar_inf >>sdum>>sdum>>sdum;
        }
        if(!first_lat) { // skip first lat.
          lat_vec.push_back(lat);
        }
        first_lat=0;
      }
      if(key=="ions per type") {
        string s;
        int id=0;
        getline(outcar_inf,s);
        aurostd::GetNextVal(s,id).c_str(); // remove an "=" sign.
        while (id<(int)s.size()) {
          int n=atoi(aurostd::GetNextVal(s,id).c_str());
          num_each_type.push_back(n);
        }
      }
      a[0]=a[1];
      a[1]=a[2];
    }
  }

}

// ***************************************************************************
// GetDatFromXDATCAR
// ***************************************************************************
// This gets the fpos from XDATCAR.
namespace pflow {
  void GetDatFromXdatcar(vector<aurostd::matrix<double> >& fpos_vec, //CO20200404 pflow::matrix()->aurostd::matrix()
      ifstream& xdatcar_inf)
  {
    fpos_vec.clear();
    string s;
    aurostd::matrix<double> dpmat(0);  //CO20200404 pflow::matrix()->aurostd::matrix()
    vector<double> dpvec(3);
    // Read in header lines
    for(int il=0;il<6;il++) {
      getline(xdatcar_inf,s);
    }
    // Read in all sets of fpos.
    int keep_reading=1;
    while (keep_reading) {
      getline(xdatcar_inf,s);
      int id;
      string key;
      while (s.size()>1 && key!="Konfig=") {
        id=0;
        for(int ic=0;ic<3;ic++) {
          dpvec[ic]=atof(aurostd::GetNextVal(s,id).c_str());
        }
        dpmat.push_back(dpvec);
        if(!getline(xdatcar_inf,s)) keep_reading=0;
        id=0;
        key=aurostd::GetNextVal(s,id);
      }
      fpos_vec.push_back(dpmat);
      dpmat=aurostd::matrix<double> (0); //CO20200404 pflow::matrix()->aurostd::matrix()
    }
  }
}

// ***************************************************************************
// GetStrVecFromOUTCAR_XDATCAR
// ***************************************************************************
namespace pflow {
  vector<xstructure> GetStrVecFromOUTCAR_XDATCAR(ifstream& outcar_inf, ifstream& xdatcar_inf) {
    vector<aurostd::matrix<double> > fpos_vec; //CO20200404 pflow::matrix()->aurostd::matrix()
    vector<aurostd::matrix<double> > lat_vec;  //CO20200404 pflow::matrix()->aurostd::matrix()
    deque<int> num_each_type;
    pflow::GetDatFromOutcar(lat_vec,num_each_type,outcar_inf);
    pflow::GetDatFromXdatcar(fpos_vec,xdatcar_inf);
    if(lat_vec.size()!=fpos_vec.size()) {
      cout << endl;
      cerr << "WARNING: RayTraceFuncs.cc/GetStrVecFromOUTCAR_XDATCAR" << endl;
      cerr << "WARNING: Number of lattice vector and positions of images are not equal." << endl;
      cerr << "WARNING: Number of lattices: " << lat_vec.size()<< endl;
      cerr << "WARNING: Number of positions: " << fpos_vec.size()<< endl;
      cerr << "WARNING: Resizing number of lattices to match number of positions and using last lattice to fill out vector if needed. " << endl;
      cout << endl;
    }
    if(lat_vec.size()>fpos_vec.size()) {lat_vec.resize(fpos_vec.size());}
    if(lat_vec.size()<fpos_vec.size()) {
      aurostd::matrix<double> tlat=lat_vec[lat_vec.size()-1];  //CO20200404 pflow::matrix()->aurostd::matrix()
      for(int i=(int)lat_vec.size();i<(int)fpos_vec.size();i++) {
        lat_vec.push_back(tlat);
      }
    }
    vector<xstructure> vstr;
    xstructure str;
    int nstr=(int)lat_vec.size();
    for(int is=0;is<nstr;is++) {
      int nat=fpos_vec[0].size();
      vector<string> names(nat,"H");
      vector<int> names_were_given(nat,FALSE);
      str=pflow::SetLat(str,lat_vec[is]);
      //DX20210118 [OBSOLETE] str=pflow::SetNumEachType(str,num_each_type);
      str.num_each_type = num_each_type; //DX20210202 - replace SetNumEachType
      str=pflow::AddAllAtomPos(str,fpos_vec[is],0);
      str=pflow::SetAllAtomNames(str,names);
      str=pflow::SetNamesWereGiven(str,names_were_given);
      vstr.push_back(str);
    }
    return vstr;
  }
}

// ***************************************************************************
// PrintStrVec
// ***************************************************************************
namespace pflow {
  void PrintStrVec(const vector<xstructure>& vstr, ostream& outf) {
    if(vstr.size()==0) return;
    for(int ist=0;ist<(int)vstr.size()-1;ist++) {
      outf << vstr[ist];
      outf << endl;
    }
    outf << vstr[vstr.size()-1];
  }
}

// ***************************************************************************
// ReadInRTParams
// ***************************************************************************
// Reads the ray tracing parameters and stores them in the rtparam object.
// Input data is stored in rtinfile using tokens of the form
// TOKEN=value.  The input can have arbitrary spaces and blank lines.
// The input can also have arbitrary stuff after value as long as it
// is separated by a space.  The value can sometimes be more than one
// field.  Comments are any sequence of text on a single line that
// follows a "#".
namespace pflow {
  void ReadInRTParams(ifstream& rtinfile, pflow::rtparams& rtp) {

    // Read in all the tokens.
    string s,s_ns;
    vector<string> token_vec;
    vector<string> val_vec;
    while (!rtinfile.eof()) {
      getline(rtinfile,s);
      s_ns=aurostd::RemoveSpaces(s); // Get string with no spaces.
      if(s_ns.size()>0) { // Make sure line was not blank
        if(s_ns[0]!='#') { // Exclude comment lines
          string token,sval;
          int id=s.find('=');
          if(id>=(int)s.length()) {
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"the following token is incorrectly formatted: "+s,_INPUT_ILLEGAL_); //CO20200624
          }
          token=s.substr(0,id);
          token=aurostd::RemoveSpaces(token);
          int i_f=std::min(s.length(),s.find('#')); // End of string
          sval=s.substr(id+1,i_f-id-1);
          token_vec.push_back(token);
          val_vec.push_back(sval);
        } //if s_ns[0]!=#
      } // if s_ns.size>0

    } // while !rtinfile.eof

    // Read in values for variables in the rtparams object.
    int first_light = 1;
    vector<double> tmp;
    for(int i=0;i<(int)token_vec.size();i++) {
      int id=0;
      string tok=token_vec[i];
      string sval;
      double val;
      int found_token=0;
      if(tok=="RESX") {
        id=0;
        val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
        rtp.resx=val;
        found_token=1;
      }// RESX
      if(tok=="RESY") {
        id=0;
        val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
        rtp.resy=val;
        found_token=1;
      }// RESY
      if(tok=="VIEWDIR") {
        // Note that there are 3 values here.
        int id=0;
        for(int ic=0;ic<3;ic++) {
          val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
          rtp.viewdir[ic]=val;
        }
        found_token=1;
        rtp.viewdir_s=1;
      }// VIEWDIR
      if(tok=="UPDIR") {
        // Note that there are 3 values here.
        int id=0;
        for(int ic=0;ic<3;ic++) {
          val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
          rtp.updir[ic]=val;
        }
        found_token=1;
        rtp.updir_s=1;
      }// UPDIR
      if(tok=="ZOOM") {
        id=0;
        val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
        rtp.zoom=val;
        found_token=1;
      }// ZOOM
      if(tok=="ASPECTRATIO") {
        id=0;
        val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
        rtp.aspectratio=val;
        found_token=1;
      }// ASPECTRATIO
      if(tok=="ANTIALIASING") {
        id=0;
        val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
        rtp.antialiasing=val;
        found_token=1;
      }// ANTIALIASING
      if(tok=="RAYDEPTH") {
        id=0;
        val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
        rtp.raydepth=val;
        found_token=1;
      }// RAYDEPTH
      if(tok=="CENTER") {
        // Note that there are 6 values here.
        int id=0;
        for(int ic=0;ic<6;ic++) {
          val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
          rtp.center_guide[ic]=val;
        }
        rtp.center_s=1;
        found_token=1;
      }// CENTER
      if(tok=="LIGHT") {
        /* This is a little confusing.  There can be multiple
           LIGHT tokens.  Each time one is found a new LIGHT is
           added to the rtparams object.  Note that there are
           7 values here.  3 for center, 1 for rad, 3 for color.
           These values must be pushed back onto the correct vectors
           and matrices.  However, for the first LIGHT token, one must
           overwrite the default rather than add a new light.  
           */
        vector<double> lcen(3);
        double lrad;
        vector<double> lcolor(3);
        int id=0;
        for(int ic=0;ic<3;ic++) {
          val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
          lcen[ic]=val;
        }
        val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
        lrad=val;
        for(int ic=0;ic<3;ic++) {
          val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
          lcolor[ic]=val;
        }
        if(first_light) {
          first_light=0;
          rtp.lightcenter[0]=lcen;
          rtp.lightrad[0]=lrad;
          rtp.lightcolor[0]=lcolor;
        }
        else {
          rtp.lightcenter.push_back(lcen);
          rtp.lightrad.push_back(lrad);
          rtp.lightcolor.push_back(lcolor);
        }
        found_token=1;
      }// LIGHT
      if(tok=="BACKGROUND") {
        // Note that there are 3 values here.
        int id=0;
        for(int ic=0;ic<3;ic++) {
          val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
          rtp.background[ic]=val;
        }
        found_token=1;
      }// BACKGROUND
      if(tok=="ATOMCOLOR") {
        // Note that there are 4 values here.
        int id=0;
        tmp = vector<double> (4);
        for(int ic=0;ic<(int)tmp.size();ic++) {
          val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
          tmp[ic]=val;
        }
        rtp.sphtex_color.push_back(tmp);
        found_token=1;
      }// ATOMCOLOR
      if(tok=="ATOMTEXTURE") {
        // Note that there are 5 values here.
        int id=0;
        tmp = vector<double> (5);
        for(int ic=0;ic<(int)tmp.size();ic++) {
          val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
          tmp[ic]=val;
        }
        rtp.sphtex_tex.push_back(tmp);
        found_token=1;
      }// ATOMTEXTURE
      if(tok=="ATOMRAD") {
        // Note that there are 2 values here.
        int id=0;
        tmp = vector<double> (2);
        for(int ic=0;ic<(int)tmp.size();ic++) {
          val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
          rtp.sph_rad.push_back(val);
        }
        found_token=1;
      }//ATOMRAD
      if(tok=="SHADING") {
        int id=0;
        rtp.shading = aurostd::GetNextVal(val_vec[i],id);
        found_token=1;
      }// SHADING
      if(tok=="OUTFILE") {
        int id=0;
        rtp.outfile = aurostd::GetNextVal(val_vec[i],id);
        found_token=1;
      }// OUTFILE
      if(tok=="SUPERCELL") {
        // Note that there are 9 values here.
        int id=0;
        for(int ic=0;ic<3;ic++) {
          for(int jc=0;jc<3;jc++) {
            val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
            rtp.sc[ic][jc]=val;
          }
        }
        found_token=1;
        rtp.sc_s=1;
      }//SUPERCELL
      if(tok=="CALCTYPE") {
        int id=0;
        rtp.calc_type=atoi(aurostd::GetNextVal(val_vec[i],id).c_str());
        found_token=1;
      }//CALCTYPE
      if(tok=="INFILE") {
        int id=0;
        while (id<(int)val_vec[i].size()) {
          string s;
          s=aurostd::GetNextVal(val_vec[i],id);
          if(s.size()>0) rtp.input_files.push_back(s);
        }
        found_token=1;
      }//INSERT_FILE
      if(tok=="INSERT_FILE") {
        int id=0;
        string s;
        s=aurostd::GetNextVal(val_vec[i],id);
        rtp.insert_file=s;
        found_token=1;
      }//INSERT_FILE
      if(tok=="ROTATION") {
        int id=0;
        for(int ic=0;ic<6;ic++) {
          string s=aurostd::GetNextVal(val_vec[i],id);
          rtp.rotation[ic] = atof(s.c_str());
          //rtp.rotation[ic] = atof(aurostd::GetNextVal(val_vec[ic],id).c_str());
        }
        found_token=1;
      }//ROTATION
      if(tok=="STRUCTURE_ORIGIN") {
        int id=0;
        for(int ic=0;ic<3;ic++) {
          string s=aurostd::GetNextVal(val_vec[i],id);
          rtp.struct_origin[ic] = atof(s.c_str());
        }
        rtp.struct_origin_s=1;
        found_token=1;
      }//STRUCTURE_ORIGIN
      if(tok=="PLANE") {
        int id=0;
        rtp.plane=atoi(aurostd::GetNextVal(val_vec[i],id).c_str());
        found_token=1;
      }//PLANE
      if(tok=="PLANECENTER") {
        // Note that there are 3 values here.
        int id=0;
        for(int ic=0;ic<(int)rtp.plane_center.size();ic++) {
          val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
          rtp.plane_center[ic]=val;
        }
        found_token=1;
        rtp.updir_s=1;
      }// PLANECENTER
      if(tok=="PLANENORMAL") {
        // Note that there are 3 values here.
        int id=0;
        for(int ic=0;ic<(int)rtp.plane_normal.size();ic++) {
          val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
          rtp.plane_normal[ic]=val;
        }
        found_token=1;
      }// PLANENORMAL
      if(tok=="PLANETEXTURE") {
        // Note that there are 4 values here.
        int id=0;
        for(int ic=0;ic<(int)rtp.planetex_tex.size();ic++) {
          val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
          rtp.planetex_tex[ic]=val;
        }
        found_token=1;
      }// PLANETEXTURE
      if(tok=="PLANECOLOR") {
        // Note that there are 3 values here.
        int id=0;
        for(int ic=0;ic<(int)rtp.plane_color.size();ic++) {
          val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
          rtp.plane_color[ic]=val;
        }
        found_token=1;
      }// PLANECOLOR

      if(!found_token) {
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"token is not recognized: "+tok,_INPUT_ILLEGAL_); //CO20200624
      }

    } // for i
  }// end routine
}

// ***************************************************************************
// ReadInStrVec
// ***************************************************************************
namespace pflow {
  void ReadInStrVec(vector<xstructure>& vstr, ifstream& strlist_inf)  {
    vstr.clear();
    stringstream sss;
    vector<string> vline,vtmp;
    aurostd::stream2vectorstring(strlist_inf,vline);
    uint iline=0;
    while (iline<vline.size()) {
      if(aurostd::RemoveWhiteSpaces(vline.at(iline))=="") {
        xstructure str(sss);
        vstr.push_back(str);
        sss.clear();sss.str("");
      } else {
        sss << vline.at(iline) << endl;
      }
      iline++;
    } 
    // [OBSOLETE]    //  int cnt=0;
    // [OBSOLETE]       xstructure str;
    // [OBSOLETE]    // Read in first structure.
    // [OBSOLETE]    strlist_inf >> str;
    // [OBSOLETE]    //  cnt++;
    // [OBSOLETE]    //  cout << "Read in structure: " << cnt << endl;
    // [OBSOLETE]    vstr.push_back(str);
    // [OBSOLETE]    // Read in all remaining structures.
    // [OBSOLETE]    while (getline(strlist_inf,dum)) {
    // [OBSOLETE]      xstructure str;
    // [OBSOLETE]      strlist_inf >> str;
    // [OBSOLETE]      //    cnt++;
    // [OBSOLETE]      //    cout << "Read in structure: " << cnt << endl;
    // [OBSOLETE]      vstr.push_back(str);
    // [OBSOLETE]    }
  }
}

// ***************************************************************************
// JoinStrVec
// ***************************************************************************
// Note that this can be made more memory efficient
// quite easily if that is needed.
namespace pflow {
  void JoinStrVec(vector<xstructure> vstr1,
      vector<xstructure> vstr2,
      vector<xstructure>& vstrtot) {
    // cout << "In JOIN" << endl;
    // cout << "strlist1 size " << vstr1.size() << endl;
    // cout << "strlist2 size " << vstr2.size() << endl;
    int size1=vstr1.size();
    int size2=vstr2.size();
    if(size1<size2) {// Pad vstr1
      xstructure str=vstr1[size1-1];
      for(int is=size1;is<size2;is++) {
        vstr1.push_back(str);
      }
    }
    if(size2<size1) {// Pad vstr2
      xstructure str=vstr2[size2-1];
      for(int is=size2;is<size1;is++) {
        vstr2.push_back(str);
      }
    }
    for(int is=0;is<size1;is++) {
      xstructure str1=vstr1[is];
      xstructure str2=vstr2[is];
      aurostd::matrix<double> lat1=pflow::GetLat(str1);  //CO20200404 pflow::matrix()->aurostd::matrix()
      xmatrix<double> xlat1(3,3);xlat1=str1.lattice;
      aurostd::matrix<double> cpos2=pflow::GetCpos(str2);  //CO20200404 pflow::matrix()->aurostd::matrix()
      deque<int> num_each_type_1=str1.num_each_type; //DX20210119 pflow::GetNumEachType() -> xstructure.num_each_type
      int num_types_1=num_each_type_1.size();
      deque<int> num_each_type_2=str2.num_each_type; //DX20210119 pflow::GetNumEachType() -> xstructure.num_each_type
      int num_types_2=num_each_type_2.size();
      // Loop over str2 atoms, assign them a type and c/d positions in an atom.
      int cnt=0;
      for(int it=0;it<num_types_2;it++) {
        for(int ia=0;ia<num_each_type_2.at(it);ia++) {
          _atom a;
          a=pflow::SetCpos(a,cpos2[cnt]);
          a.fpos=C2F(xlat1,a.cpos);
          //	a=pflow::SetFpos(a,C2F(lat1,cpos2[cnt]));
          a=pflow::SetType(a,num_types_1+it);
          str1.AddAtom(a,false);  //CO20230319 - add by type
          cnt++;
        }
      }
      //      cout << str1.atoms.size() << endl;
      vstrtot.push_back(str1);
    }//is
  }//end routine
}

// ***************************************************************************
// SetStrFromRTParams
// ***************************************************************************
// This assumes you have read in the rtparams with
// SetRTParams.  It uses both a structure and the rtparams
// and alters the structure based on rtparams.  Here the
// structure is made into a supercell if needed and the
// struc_origin is set.
namespace pflow {
  void SetStrFromRTParams(xstructure& str, pflow::rtparams& rtp) {
    // Make structral adjustments
    if(rtp.sc_s) {
      xmatrix<double> _supercell(3,3);
      _supercell=aurostd::matrix2xmatrix(rtp.sc);  //CO20200404 pflow::matrix()->aurostd::matrix()
      // Make a supercell if it was set by user.  Note that this moves all atoms into the unit cell.
      str=GetSuperCell(str,_supercell);
      // Reset rtparams now that you have a supercell.
      pflow::SetRTParams(str,rtp);
    }
    str=pflow::SetOrigin(str,rtp.struct_origin);
  }
}
// ***************************************************************************
// SuperCellStrVec
// ***************************************************************************
// Gets supercells for every strucutre in the structure vector.
namespace pflow {
  void SuperCellStrVec(vector<xstructure>& vstr, const aurostd::matrix<double>& sc) {  //CO20200404 pflow::matrix()->aurostd::matrix()
    if(vstr.size()==0) return;
    for(int is=0;is<(int)vstr.size();is++) {
      xmatrix<double> _supercell(3,3);
      _supercell=aurostd::matrix2xmatrix(sc);  //CO20200404 pflow::matrix()->aurostd::matrix()
      vstr[is]=GetSuperCell(vstr[is],_supercell); // Makes a supercell.
    }
  }
}
// ***************************************************************************
// UpDateRTParams
// *************************************************
// This simply makes changes to rtparams as the program steps
// through different frames.  This is use for values that must
// evolve during the movie (e.g., the center).
namespace pflow {
  void UpDateRTParams(pflow::rtparams& rtp, const int& istr, int nstr) {
    // Update center
    // tpx tpx ??? (cut time to 1/2 movie for moving center by /2).
    //nstr=(int) (nstr/1.5);
    double r=0;
    if(nstr>1) {
      r = (double)istr/(double)(nstr-1);
    }
    if(r>1) r=1;
    if(rtp.center_s) {
      for(int ic=0;ic<3;ic++) {
        rtp.center[ic]=rtp.center_guide[2*ic]+r*(rtp.center_guide[2*ic+1]-rtp.center_guide[2*ic]);
      }
    }
    //tpx
    cout << "CENTER " ;
    pflow::Vout(rtp.center,cout);
  }
}

// ***************************************************************************
// SetRTParams
// ***************************************************************************
// This assumes you have read in the rtparams with
// SetRTParams.  It uses both a structure and the rtparams
// and alters rtparams based on what has been read in and the
// structure.  The rtparam defaults are set based on the
// ntypes and other structural characteristics.
namespace pflow {
  void SetRTParams(xstructure& str, pflow::rtparams& rtp) {
    // Now set all defaults for rtp that
    // use structural information.
    str=ReScale(str,1.0);
    aurostd::matrix<double> lat=pflow::GetLat(str);  //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> cpos=pflow::GetCpos(str);  //CO20200404 pflow::matrix()->aurostd::matrix()
    int nat=cpos.size();

    // For all items that must be set for each type we do that here.
    // These include sphtex_* and sph_rad.
    // Note that we protect against reading in data for atom types that do not exist.
    // Also note that this makes of unique formats for the rtparam
    // data that only exist the when they are read in.  Therefore,
    // we must only do this the first time this set function is
    // called.  This is controlled by the first_set parameter.
    if(rtp.first_set) {
      int ntypes = str.num_each_type.size(); //DX20210119 pflow::GetNumEachType() -> xstructure.num_each_type
      aurostd::matrix<double> tmp; //CO20200404 pflow::matrix()->aurostd::matrix()
      int size;
      // Set sphtex_tex
      size=rtp.sphtex_tex_def.size();
      tmp=aurostd::matrix<double> (ntypes,rtp.sphtex_tex_def); //CO20200404 pflow::matrix()->aurostd::matrix()
      for(int i=0;i<(int)rtp.sphtex_tex.size();i++) {
        int id=(int)rtp.sphtex_tex[i][0]-1;
        for(int j=1;j<size+1;j++) {
          if(id>=0 && id<ntypes) tmp[id][j-1]=rtp.sphtex_tex[i][j];
        }
      }
      rtp.sphtex_tex=tmp;
      // Set sphtex_color
      size=rtp.sphtex_color_def.size();
      tmp=aurostd::matrix<double> (ntypes,rtp.sphtex_color_def); //CO20200404 pflow::matrix()->aurostd::matrix()
      for(int i=0;i<(int)rtp.sphtex_color.size();i++) {
        int id=(int)rtp.sphtex_color[i][0]-1;
        for(int j=1;j<size+1;j++) {
          if(id>=0 && id<ntypes) tmp[id][j-1]=rtp.sphtex_color[i][j];
        }
      }
      rtp.sphtex_color=tmp;
      // Set sphtex_phong
      size=rtp.sphtex_phong_def.size();
      tmp=aurostd::matrix<double> (ntypes,rtp.sphtex_phong_def); //CO20200404 pflow::matrix()->aurostd::matrix()
      for(int i=0;i<(int)rtp.sphtex_phong.size();i++) {
        int id=(int)rtp.sphtex_phong[i][0]-1;
        for(int j=1;j<size+1;j++) {
          if(id>=0 && id<ntypes) tmp[id][j-1]=rtp.sphtex_phong[i][j];
        }
      }
      rtp.sphtex_phong=tmp;
      // Set sph_rad
      size=1;
      vector<double> vtmp(ntypes,rtp.sph_rad_def);
      for(int i=0;i<(int)rtp.sph_rad.size()/2;i++) {
        int id=(int)rtp.sph_rad[2*i]-1;
        if(id>=0 && id<ntypes) vtmp[id]=rtp.sph_rad[2*i+1];
      }
      rtp.sph_rad=vtmp;
      // Set sphtex_names
      rtp.sphtex_names = vector<string> (ntypes,"txt_atom_type_");
      for(int it=0;it<ntypes;it++) {
        ostringstream os;
        os << it+1 << ends;
        string s(os.str());
        rtp.sphtex_names[it]=rtp.sphtex_names[it]+s;
      }
      rtp.first_set=0;
    }// if first_set    

    //  Set view direction along the 0 axis.
    //  if(!rtp.viewdir_s) {rtp.viewdir=lat[0];}

    //  Set up direction along the 1 axis.
    //  if(!rtp.updir_s) {rtp.updir=lat[0];}

    /*
       Set up struct_origin, which is the point around which rotation
       will take place.  Here we default to the first moment of the atom
       positions.
       */
    if(! rtp.struct_origin_s) {// If struct_orig was not set by user
      // rtp.struct_origin=SVprod(0.5,VVsum(lat[0],VVsum(lat[1],lat[2])));
      vector<double> mom1=xvector2vector(GetMom1(str));
      rtp.struct_origin=mom1;
    }

    /*
       Set up center (where camera is located).
       Define the center to be displaced from the centroid along the view
       direction. Define center to be located a distance d=-2.0*R/tan(45)=-2.0*R
       from the centroid, where R is the the largest projected distance
       of an atom from the centroid (first moment) in the image plane (the
       plane defined to be perpendicular to the view direction and containing the
       centroid).  So define
       mom1=first moment
       vd=unit vector along view direction
       up=unit vector along up direction
       upp=unit vector consisting of projection of up vector into the image plane.
       vdp=unit vector in image plane perpendicular to vd and updir.

       Actually, instead of mom1 being the centroid just set it to the
       struct_origin.  All else is the same.
       */
    if(! rtp.center_s) {// If center was not set by user
      vector<double> mom1=rtp.struct_origin;
      vector<double> vd=pflow::SVprod(1.0/pflow::norm(rtp.viewdir),rtp.viewdir);
      vector<double> up=pflow::SVprod(1.0/pflow::norm(rtp.updir),rtp.updir);
      vector<double> upp;
      upp=pflow::VVsum(up,pflow::SVprod(-pflow::VVprod(up,vd),vd));
      upp=pflow::SVprod(1.0/pflow::norm(upp),upp);
      vector<double> vdp=pflow::VVcross(vd,upp);
      double proj=0;
      //      double max_iat=0;
      for(int iat=0;iat<nat;iat++) {
        double p1=pflow::VVprod(pflow::VVdiff(cpos[iat],mom1),vdp);
        double p2=pflow::VVprod(pflow::VVdiff(cpos[iat],mom1),upp);
        double nproj=sqrt(p1*p1+p2*p2);
        if(nproj>proj) {
          proj=nproj;
          //	  max_iat=iat;
        }
      }
      double displ=-2.0*proj;
      rtp.center=pflow::VVsum(mom1,pflow::SVprod(displ,vd));
    }// if ! center_s
    else { // if center is set by user
      for(int ic=0;ic<3;ic++) {
        rtp.center[ic]=rtp.center_guide[2*ic];
      }
    }
  }
}

// ***************************************************************************
// GetRTDatFile
// ***************************************************************************
namespace pflow {
  void GetRTDatFile(xstructure str, const pflow::rtparams& rtp,
      ostringstream& rtdat_file) {
    rtdat_file << "BEGIN_SCENE" << endl;
    rtdat_file << "  RESOLUTION " << rtp.resx << " " << rtp.resy << endl;

    //string tga="tga";
    //string outfile;
    //outfile=rtp.outfile+"_"+tga;
    //rtdat_file << "  OUTFILE " << outfile << endl;

    rtdat_file << endl;
    rtdat_file << "CAMERA" << endl;
    rtdat_file << "  ZOOM " << rtp.zoom << endl;
    rtdat_file << "  ASPECTRATIO " << rtp.aspectratio <<  endl;
    rtdat_file << "  ANTIALIASING " << rtp.antialiasing << endl;
    rtdat_file << "  RAYDEPTH " << rtp.raydepth << endl;
    rtdat_file << "  CENTER " << rtp.center[0] << " " << rtp.center[1] << " " << rtp.center[2] << endl;
    rtdat_file << "  VIEWDIR " << rtp.viewdir[0] << " " << rtp.viewdir[1] << " " << rtp.viewdir[2] << endl;
    rtdat_file << "  UPDIR " << rtp.updir[0] << " " << rtp.updir[1] << " " << rtp.updir[2] << endl;
    rtdat_file <<  endl;
    rtdat_file << "END_CAMERA" << endl;
    rtdat_file <<  endl;
    rtdat_file << "BACKGROUND " << rtp.background[0] << " " << rtp.background[1] << " " << rtp.background[2] << endl;
    rtdat_file <<  endl;
    for(int il=0;il<(int)rtp.lightcenter.size();il++) {
      rtdat_file << "LIGHT ";
      rtdat_file << "CENTER " << rtp.lightcenter[il][0] << " " << rtp.lightcenter[il][1] << " " << rtp.lightcenter[il][2];
      rtdat_file << " RAD " << rtp.lightrad[il];
      rtdat_file << " COLOR " << rtp.lightcolor[il][0] << " " << rtp.lightcolor[il][1] << " " << rtp.lightcolor[il][2];
      rtdat_file <<  endl << endl;
    }
    for(int iat=0;iat<(int) rtp.sphtex_names.size();iat++) {
      rtdat_file << "TEXDEF " << rtp.sphtex_names[iat] << " ";
      rtdat_file << "AMBIENT " << rtp.sphtex_tex[iat][0] << " ";
      rtdat_file << "DIFFUSE " << rtp.sphtex_tex[iat][1] << " ";
      rtdat_file << "SPECULAR " << rtp.sphtex_tex[iat][2] << " ";
      rtdat_file << "OPACITY " << rtp.sphtex_tex[iat][3] << " ";
      rtdat_file << endl;
      rtdat_file << "PHONG PLASTIC " << rtp.sphtex_phong[iat][0] << " " << "PHONG_SIZE " << rtp.sphtex_phong[iat][1] << endl;
      rtdat_file << "  COLOR " << rtp.sphtex_color[iat][0] << " " << rtp.sphtex_color[iat][1] << " " << rtp.sphtex_color[iat][2] << endl;
      rtdat_file << "  TEXFUNC 0" << endl;
      rtdat_file << endl;
    }
    // Do atoms
    str=ReScale(str,1.0);
    aurostd::matrix<double> cpos=pflow::GetCpos(str);  //CO20200404 pflow::matrix()->aurostd::matrix()
    deque<int> num_each_type=str.num_each_type; //DX20210119 pflow::GetNumEachType() -> xstructure.num_each_type
    //int nat=(int)cpos.size();
    int ntype=(int)num_each_type.size();
    int cnt=0;
    for(int it=0;it<ntype;it++) {
      for(int iat=0;iat<num_each_type.at(it);iat++) {
        rtdat_file << "SPHERE CENTER " << cpos[cnt][0] << " " << cpos[cnt][1] << " " << cpos[cnt][2];
        rtdat_file << " RAD " << rtp.sph_rad[it];
        rtdat_file << " " << rtp.sphtex_names[it] << endl;
        rtdat_file << endl;
        cnt++;
      }
    }
    // Do plane
    if(rtp.plane==1) {
      rtdat_file << "PLANE" << endl;
      rtdat_file << "  CENTER " << rtp.plane_center[0] << " " << rtp.plane_center[1] << " " << rtp.plane_center[2] << endl;
      rtdat_file << "  NORMAL " << rtp.plane_center[0] << " " << rtp.plane_center[1] << " " << rtp.plane_center[2] << endl;
      rtdat_file << "  TEXTURE " << endl;
      rtdat_file << "    AMBIENT " << rtp.planetex_tex[0] << " ";
      rtdat_file << "DIFFUSE " << rtp.planetex_tex[1] << " ";
      rtdat_file << "SPECULAR " << rtp.planetex_tex[2] << " ";
      rtdat_file << "OPACITY " << rtp.planetex_tex[3] << " ";
      rtdat_file << "  COLOR " << rtp.plane_color[0] << " " << rtp.plane_color[1] << " " << rtp.plane_color[2] << endl;
      rtdat_file << "  TEXFUNC 0" << endl;
      cout << endl;
    }

    rtdat_file << "END_SCENE" << endl;
    rtdat_file << ends;
  }
}

// ***************************************************************************
// PrintRDatFile
// ***************************************************************************
namespace pflow {
  string PrintRTDatFile(ostringstream& rtdat_file, const pflow::rtparams& rt_params) {  
    string filename = rt_params.outfile+".dat";
    ostringstream outss;
    outss << filename << ends;
    ofstream outf(outss.str().c_str());
    outf << rtdat_file.str();
    return filename;
  }
}

// ***************************************************************************
// CreateRTtgaFile
// ***************************************************************************
namespace pflow {
  string CreateRTtgaFile(const string& datfile, const pflow::rtparams& rt_params) {  
    // Run tachyon on datfile.
    ostringstream tachyon_stream;
    tachyon_stream << "tachyon " << "-" << rt_params.shading << " " << datfile << ends;
    // system(tachyon_stream.str().c_str());
    aurostd::execute(tachyon_stream);
    // Move outfile.tga to user defined name.
    ostringstream mv_stream;
    string tgafile = rt_params.outfile+".tga";
    mv_stream << "mv outfile.tga " << tgafile << ends;
    aurostd::execute(mv_stream);
    //  system(mv_stream.str().c_str());
    return tgafile;
  }
}

// ***************************************************************************
// CreateRTjpgFile
// ***************************************************************************
namespace pflow {
  string CreateRTjpgFile(const string& tgafile, const pflow::rtparams& rt_params) {  
    // Run convert.
    ostringstream convert_stream;
    string jpgfile = rt_params.outfile+".jpg";
    convert_stream << XHOST.command("convert") << " " << tgafile << " " << jpgfile << ends;
    aurostd::execute(convert_stream);
    // system(convert_stream.str().c_str());
    return jpgfile;
  }
}

// ***************************************************************************
// GetRTencFile
// ***************************************************************************
namespace pflow {
  void GetRTencFile(const pflow::rtparams& rtp, const int nim,
      ostringstream& os) {
    os << "PATTERN         IBBPBBPBBPBBPBB" << endl;
    os << "OUTPUT          "<< rtp.outfile << ".mpg" << endl;
    os << "GOP_SIZE        16" << endl;
    os << "SLICES_PER_FRAME        1" << endl;
    os << "BASE_FILE_FORMAT        PPM" << endl;
    os << endl;
    os << "INPUT_CONVERT   djpeg *" << endl;
    os << "INPUT_DIR       ." << endl;
    os << "INPUT" << endl;
    os << rtp.outfile << "_*.jpg       [" << aurostd::PaddedNumString(0,5) << "-" << nim-1 << "+1]" << endl;
    os << "END_INPUT" << endl;
    os << endl;
    os << "IQSCALE         8" << endl;
    os << "PQSCALE         10" << endl;
    os << "BQSCALE         25" << endl;
    os << "PIXEL           HALF" << endl;
    os << "RANGE           10" << endl;
    os << "PSEARCH_ALG     LOGARITHMIC" << endl;
    os << "BSEARCH_ALG     CROSS2" << endl;
    os << "REFERENCE_FRAME DECODED" << endl;
    os << ends;
  }
}

// ***************************************************************************
// PrintRTencFile
// ***************************************************************************
namespace pflow {
  string PrintRTencFile(const pflow::rtparams& rt_params, ostringstream& rtenc_file) {
    string filename = rt_params.outfile+".enc";
    ostringstream outss;
    outss << filename << ends;
    ofstream outf(outss.str().c_str());
    outf << rtenc_file.str();
    return filename;
  }
}

// ***************************************************************************
// CreateRTmgpFile
// ***************************************************************************
namespace pflow {
  string CreateRTmpgFile(const pflow::rtparams& rt_params, const string& encfile) {
    // Run mpeg_encoder on enc_file.
    ostringstream encoder_stream;
    encoder_stream << "mpeg_encode " << " " << encfile << ends;
    /*  char* encoder_char = encoder_stream.str().c_str();
        system(encoder_char);*/
    system(encoder_stream.str().c_str());
    // Get outfile.mpg name
    string file;
    file = rt_params.outfile+".mpg";
    return file;
  }
}

// ***************************************************************************
// RayTraceManager
// ***************************************************************************
namespace pflow {
  void RayTraceManager(vector<string> argv) {
    ifstream rtinfile(argv.at(2).c_str()); // File where RT params are input.
    aurostd::InFileExistCheck("RayTraceFuncs.cc/RayTraceManager",argv.at(2),rtinfile);
    pflow::rtparams rtp; // Object that stores RT params.
    ReadInRTParams(rtinfile,rtp); // Sets RT params object from input.
    switch (rtp.calc_type) {
      case 0:{ // A single file of structures.
               // Read in structure list.
               if(rtp.input_files.size()<1) {
                 throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"must specify a structure list file to open with the token INFILE (e.g., INFILE = STRLIST)",_FILE_CORRUPT_); //CO20200624
               }
               ifstream strlist_inf(rtp.input_files[0].c_str());
               aurostd::InFileExistCheck("RayTraceFuncs.cc/RayTraceManager",rtp.input_files[0].c_str(),strlist_inf);
               vector<xstructure> vstr;
               cout << endl;
               cout << "Reading in structure list." << endl;
               cout << endl;
               pflow::ReadInStrVec(vstr,strlist_inf);
               int nstr=vstr.size();
               cout << "We are working with this many structures: " << nstr << endl;
               // Set rtparams.
               pflow::SetRTParams(vstr[0],rtp);
               // Set structures from rtparams.
               for(int is=0;is<nstr;is++) {
                 pflow::SetStrFromRTParams(vstr[is],rtp);
               }
               cout << endl;
               cout << "Here are the values of the ray trace parameters." << endl;
               rtp.Print();
               cout << endl;
               // Rotate vstr
               RotateStrVec(vstr,rtp.rotation);
               // Create dat,tga,jpg files.
               string datfile,tgafile,jpgfile;
               for(int is=0;is<nstr;is++) {
                 cout << endl;
                 cout << "Creating dat/tga/jpg file num " << is+1 << " out of a total of " << nstr << endl;
                 cout << endl;
                 // Update rtparams
                 pflow::UpDateRTParams(rtp,is,nstr);
                 // Get dat file stringstream.
                 ostringstream rtdat_file;
                 xstructure str=vstr[is];
                 pflow::GetRTDatFile(str,rtp,rtdat_file);       // Puts formatted .dat file into a stringstream.
                 datfile=pflow::PrintRTDatFile(rtdat_file,rtp); // Prints outfile.dat file and returns name.
                 tgafile=pflow::CreateRTtgaFile(datfile,rtp);   // Creates outfile.tga file and returns name.
                 jpgfile=pflow::CreateRTjpgFile(tgafile,rtp);   // Creates outfile.jpg file and returns name.
                 // Copy files to proper names (outfile.num.suffix).
                 ostringstream mvcmd_dat;
                 ostringstream mvcmd_tga;
                 ostringstream mvcmd_jpg;
                 string name;
                 string num = aurostd::PaddedNumString(is,5);
                 string datname = (rtp.outfile+"_"+num+".dat");
                 mvcmd_dat << "mv " << datfile << " " << datname << ends;
                 system(mvcmd_dat.str().c_str());
                 string tganame = (rtp.outfile+"_"+num+".tga");
                 mvcmd_tga << "mv " << tgafile << " " << tganame << ends;
                 system(mvcmd_tga.str().c_str());
                 string jpgname = (rtp.outfile+"_"+num+".jpg");
                 mvcmd_jpg << "mv " << jpgfile << " " << jpgname << ends;
                 system(mvcmd_jpg.str().c_str());
                 if(nstr>1) { // Remove dat,tga files if there are a lot of them.
                   ostringstream rmcmd_dat;
                   ostringstream rmcmd_tga;
                   rmcmd_dat << "/bin/rm " << datname << ends;
                   system(rmcmd_dat.str().c_str());
                   //	rmcmd_tga << "/bin/rm " << tganame << ends;
                   //	system(rmcmd_tga.str());
                 }
               }
               // Get and print enc file for mpeg encoder.
               ostringstream rtenc_file;
               GetRTencFile(rtp,nstr,rtenc_file);
               string encfile;
               encfile=PrintRTencFile(rtp,rtenc_file);
               // Create mpg file with mpeg encoder.
               cout << endl;
               cout << "Creating mpg file " << rtp.outfile << ".mpg" << endl;
               cout << endl;
               string mpgfile;
               mpgfile=CreateRTmpgFile(rtp,encfile);
               break;
             }
      default:{
                throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"invalid CALCTYPE = "+aurostd::utype2string(rtp.calc_type),_INPUT_ILLEGAL_); //CO20200624
              }
    }//switch
  }
}

// ***************************************************************************
// GetRotationMatrix [OBSOLETE - moved to xatom]
// ***************************************************************************

// ***************************************************************************
// RotateStrVec //DX20210127 [OBSOLETE - moved to xatom]
// ***************************************************************************

// ***************************************************************************
// PROJDATA PROJDATA PROJDATA PROJDATA PROJDATA PROJDATA PROJDATA PROJDATA PRO
// ***************************************************************************
namespace pflow {
  // Constructors
  projdata::projdata()
  {
    sp=0;
    rspin=2;
    nl_max=4; // for s,p,d,f orbitals
    nlm_max=16; // for s,p,d,f orbitals
    nlmtot_max=20; // for s,p,d,f orbitals + p,d,f,all totals
    nl=nl_max; // Default to max values
    nlm=nlm_max; // Default to max values
    nlmtot=nlmtot_max; // Default to max values
    LMnames = vector<string> (nlm);
    LMnames[0]="S     ";
    LMnames[1]="Py    ";
    LMnames[2]="Pz    ";
    LMnames[3]="Px    ";
    LMnames[4]="Dxy   ";
    LMnames[5]="Dyz   ";
    LMnames[6]="Dz2   ";
    LMnames[7]="Dxz   ";
    LMnames[8]="Dx2-y2";
    LMnames[9]="F1    ";
    LMnames[10]="F2    ";
    LMnames[11]="F3    ";
    LMnames[12]="F4    ";
    LMnames[13]="F5    ";
    LMnames[14]="F6    ";
    LMnames[15]="F7    ";
    Lnames = vector<string> (nl);
    Lnames[0]="S     ";
    Lnames[1]="P     ";
    Lnames[2]="D     ";
    Lnames[3]="F     ";
    LLMnames = vector<string> (nlm+4);
    LLMnames[0]="S     ";
    LLMnames[1]="Py    ";
    LLMnames[2]="Pz    ";
    LLMnames[3]="Px    ";
    LLMnames[4]="Ptot  ";
    LLMnames[5]="Dxy   ";
    LLMnames[6]="Dyz   ";
    LLMnames[7]="Dz2   ";
    LLMnames[8]="Dxz   ";
    LLMnames[9]="Dx2-y2";
    LLMnames[10]="Dtot  ";
    LLMnames[11]="F1    ";
    LLMnames[12]="F2    ";
    LLMnames[13]="F3    ";
    LLMnames[14]="F4    ";
    LLMnames[15]="F5    ";
    LLMnames[16]="F6    ";
    LLMnames[17]="F7    ";
    LLMnames[18]="Ftot  ";
    LLMnames[19]="Tot   ";
  }

  void projdata::Print(ostream& outf) {
    outf.precision(3);
    outf.setf(std::ios::fixed,std::ios::floatfield);
    outf.setf(std::ios::left,std::ios::adjustfield);

    int ioncnt;
    int w1=5;
    int w2=6;
    string key;

    outf << endl;
    outf << "**************************************************" << endl;
    outf << endl;
    outf << "WARNING:  There is some confusing business about these projections." << endl;
    outf << "WARNING:  To get results that match the vasp summations you may need" << endl;
    outf << "WARNING:  to alter the vasp source code.  See the convasp help by" << endl;
    outf << "WARNING:  typing convasp -h .  Read the section on -pocc." << endl;
    outf << endl;
    outf << "The basic data" << endl;
    outf << "nlm " << nlm << endl;
    outf << "nkpts " << nkpts << endl;
    outf << "nbands " << nbands << endl;
    outf << "nions " << nions << endl;
    outf << "ntypes " << ntypes << endl;
    outf << "type  " << "num of that type " << endl;
    for(int it=0;it<ntypes;it++) {
      outf << it << "         " << num_each_type.at(it) << endl;
    }
    outf << "Spin polarization (0=no,1=yes) " << sp << endl;
    outf << "rspin (spin degeneracy) " << rspin << endl;
    outf << "Fermi Weights (nbands X nkpts)" << endl;
    outf << "Band  kpt1  kpt2 ..." << endl;
    outf << "Up spin" << endl;
    for(int ib=0;ib<(int)wfermi_u.size();ib++) {
      outf << " " << ib;
      for(int ik=0;ik<(int)wfermi_u[ib].size();ik++) {
        outf << "   " << wfermi_u[ib][ik];
      }
      outf << endl;
    }
    outf << "Down spin" << endl;
    for(int ib=0;ib<(int)wfermi_d.size();ib++) {
      outf << " " << ib;
      for(int ik=0;ik<(int)wfermi_d[ib].size();ik++) {
        outf << "   " << wfermi_d[ib][ik];
      }
      outf << endl;
    }

    outf << "kpt values and weights" << endl;
    for(int ik=0;ik<(int)wkpt.size();ik++) {
      outf << "  " << ik+1 << "  "
        << kpts[ik][0] << "  "
        << kpts[ik][1] << "  "
        << kpts[ik][2] << "      "
        << wkpt[ik] << endl;
    }

    outf << "Band energies" << endl;
    outf << "   " << "kpt   " << "band  " << "energy (up/dn) " << endl;
    for(int ik=0;ik<(int)ener_k_b_u.size();ik++) {
      for(int ib=0;ib<(int)ener_k_b_u[ik].size();ib++) {
        outf << "   " << setw(5) << ik+1 << " " << setw(5) << ib+1 << " " << setprecision(5) << ener_k_b_u[ik][ib];
        if(sp) outf << "  " << setprecision(5) << ener_k_b_d[ik][ib];
        outf << endl;
      }
    }

    outf.precision(3);

    outf << "Projections" << endl;
    outf << "   " << "type  " << "ion   " << "kpt   " << "band  " << "lm    " << "projection (up/dn)" << endl;
    ioncnt=-1;
    for(int it=0;it<(int)pdat_u.size();it++) {
      for(int iat=0;iat<num_each_type.at(it);iat++) {
        ioncnt++;
        for(int ik=0;ik<(int)pdat_u[it].size();ik++) {
          for(int ib=0;ib<(int)pdat_u[it][ik][iat].size();ib++) {
            for(int ilm=0;ilm<(int)pdat_u[it][ik][iat][ib].size();ilm++) {
              outf << "   ";
              outf << setw(w2) << it+1;
              outf << setw(w2) << ioncnt+1;
              outf << setw(w2) << ik+1;
              outf << setw(w2) << ib+1;
              outf << "   " << LMnames[ilm];
              outf << "  " << setw(w2) << pdat_u[it][ik][iat][ib][ilm];
              if(sp) outf << "  " << setw(w2) << pdat_d[it][ik][iat][ib][ilm];
              outf << endl;
            }
          }
        }
      }
    }

    // Occupations ION:KPT:BND:LM
    outf << endl;
    outf << "**************************************************" << endl;
    outf << endl;
    outf << "Occupations vs. ION:KPT:BND:LM" << endl;
    outf << endl;
    key="ION:KPT:BND:LM";
    outf << "Up spin"<< endl;
    outf << "   " << "ion   " << "type  " << "kpt   " << "band  ";
    for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
      outf << LLMnames[ilm] << "  ";
    }
    outf << endl;
    ioncnt=-1;
    for(int it=0;it<ntypes;it++) {
      for(int iat=0;iat<num_each_type.at(it);iat++) {
        ioncnt++;
        for(int ik=0;ik<nkpts;ik++) {
          for(int ib=0;ib<nbands;ib++) {
            outf << "   " << setw(w1) << ioncnt+1;
            outf << " " << setw(w1) << it+1;
            outf << " " << setw(w1) << ik+1;
            outf << " " << setw(w1) << ib+1;
            outf << " " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][0];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][1];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][2];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][3];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][1];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][4];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][5];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][6];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][7];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][8];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][2];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][9];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][10];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][11];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][12];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][13];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][14];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][15];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][3];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][4];
            outf << " " << key;
            outf << endl;
          }
        }
      }
    }
    if(sp) {
      outf << "Down spin"<< endl;
      outf << "   " << "ion   " << "type  " << "kpt   " << "band  ";
      for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
        outf << LLMnames[ilm] << "  ";
      }
      outf << endl;
      ioncnt=-1;
      for(int it=0;it<ntypes;it++) {
        for(int iat=0;iat<num_each_type.at(it);iat++) {
          ioncnt++;
          for(int ik=0;ik<nkpts;ik++) {
            for(int ib=0;ib<nbands;ib++) {
              outf << "   " << setw(w1) << ioncnt+1;
              outf << " " << setw(w1) << it+1;
              outf << " " << setw(w1) << ik+1;
              outf << " " << setw(w1) << ib+1;
              outf << " " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][0];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][1];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][2];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][3];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][1];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][4];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][5];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][6];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][7];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][8];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][2];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][9];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][10];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][11];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][12];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][13];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][14];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][15];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][3];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][4];
              outf << " " << key;
              outf << endl;
            }
          }
        }
      }
      outf << "Up-Down spin"<< endl;
      outf << "   " << "ion   " << "type  " << "kpt   " << "band  ";
      for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
        outf << LLMnames[ilm] << "  ";
      }
      outf << endl;
      ioncnt=-1;
      for(int it=0;it<ntypes;it++) {
        for(int iat=0;iat<num_each_type.at(it);iat++) {
          ioncnt++;
          for(int ik=0;ik<nkpts;ik++) {
            for(int ib=0;ib<nbands;ib++) {
              outf << "   " << setw(w1) << ioncnt+1;
              outf << " " << setw(w1) << it+1;
              outf << " " << setw(w1) << ik+1;
              outf << " " << setw(w1) << ib+1;
              outf << " " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][0]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][0];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][1]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][1];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][2]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][2];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][3]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][3];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][1]-occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][1];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][4]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][4];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][5]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][5];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][6]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][6];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][7]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][7];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][8]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][8];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][2]-occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][2];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][9]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][9];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][10]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][10];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][11]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][11];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][12]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][12];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][13]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][13];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][14]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][14];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][15]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][15];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][3]-occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][3];
              outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][4]-occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][4];
              outf << " " << key;
              outf << endl;
            }
          }
        }
      }
    }

    // Occupations ION:BND:LM
    outf << endl;
    outf << "**************************************************" << endl;
    outf << endl;
    outf << "Occupations vs. ION:BND:LM" << endl;
    outf << endl;
    key="ION:BND:LM";
    outf << "Up spin"<< endl;
    outf << "   " << "ion   " << "type  " << "band  ";
    for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
      outf << LLMnames[ilm] << "  ";
    }
    outf << endl;
    ioncnt=-1;
    for(int it=0;it<ntypes;it++) {
      for(int iat=0;iat<num_each_type.at(it);iat++) {
        ioncnt++;
        for(int ib=0;ib<nbands;ib++) {
          outf << "   " << setw(w1) << ioncnt+1;
          outf << " " << setw(w1) << it+1;
          outf << " " << setw(w1) << ib+1;
          outf << " " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][0];
          outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][1];
          outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][2];
          outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][3];
          outf << "  " << setw(w2) << occ_vs_ion_bnd_l_u[ioncnt][ib][1];
          outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][4];
          outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][5];
          outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][6];
          outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][7];
          outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][8];
          outf << "  " << setw(w2) << occ_vs_ion_bnd_l_u[ioncnt][ib][2];
          outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][9];
          outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][10];
          outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][11];
          outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][12];
          outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][13];
          outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][14];
          outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][15];
          outf << "  " << setw(w2) << occ_vs_ion_bnd_l_u[ioncnt][ib][3];
          outf << "  " << setw(w2) << occ_vs_ion_bnd_l_u[ioncnt][ib][4];
          outf << " " << key;
          outf << endl;
        }
      }
    }
    if(sp) {
      outf << "Down spin"<< endl;
      outf << "   " << "ion   " << "type  " << "band  ";
      for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
        outf << LLMnames[ilm] << "  ";
      }
      outf << endl;
      ioncnt=-1;
      for(int it=0;it<ntypes;it++) {
        for(int iat=0;iat<num_each_type.at(it);iat++) {
          ioncnt++;
          for(int ib=0;ib<nbands;ib++) {
            outf << "   " << setw(w1) << ioncnt+1;
            outf << " " << setw(w1) << it+1;
            outf << " " << setw(w1) << ib+1;
            outf << " " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][0];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][1];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][2];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][3];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_l_d[ioncnt][ib][1];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][4];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][5];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][6];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][7];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][8];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_l_d[ioncnt][ib][2];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][9];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][10];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][11];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][12];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][13];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][14];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][15];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_l_d[ioncnt][ib][3];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_l_d[ioncnt][ib][4];
            outf << " " << key;
            outf << endl;
          }
        }
      }
      outf << "Up-Down spin"<< endl;
      outf << "   " << "ion   " << "type  " << "band  ";
      for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
        outf << LLMnames[ilm] << "  ";
      }
      outf << endl;
      ioncnt=-1;
      for(int it=0;it<ntypes;it++) {
        for(int iat=0;iat<num_each_type.at(it);iat++) {
          ioncnt++;
          for(int ib=0;ib<nbands;ib++) {
            outf << "   " << setw(w1) << ioncnt+1;
            outf << " " << setw(w1) << it+1;
            outf << " " << setw(w1) << ib+1;
            outf << " " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][0]-occ_vs_ion_bnd_lm_d[ioncnt][ib][0];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][1]-occ_vs_ion_bnd_lm_d[ioncnt][ib][1];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][2]-occ_vs_ion_bnd_lm_d[ioncnt][ib][2];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][3]-occ_vs_ion_bnd_lm_d[ioncnt][ib][3];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_l_u[ioncnt][ib][1]-occ_vs_ion_bnd_l_d[ioncnt][ib][1];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][4]-occ_vs_ion_bnd_lm_d[ioncnt][ib][4];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][5]-occ_vs_ion_bnd_lm_d[ioncnt][ib][5];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][6]-occ_vs_ion_bnd_lm_d[ioncnt][ib][6];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][7]-occ_vs_ion_bnd_lm_d[ioncnt][ib][7];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][8]-occ_vs_ion_bnd_lm_d[ioncnt][ib][8];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_l_u[ioncnt][ib][2]-occ_vs_ion_bnd_l_d[ioncnt][ib][2];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][9]-occ_vs_ion_bnd_lm_d[ioncnt][ib][9];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][10]-occ_vs_ion_bnd_lm_d[ioncnt][ib][10];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][11]-occ_vs_ion_bnd_lm_d[ioncnt][ib][11];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][12]-occ_vs_ion_bnd_lm_d[ioncnt][ib][12];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][13]-occ_vs_ion_bnd_lm_d[ioncnt][ib][13];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][14]-occ_vs_ion_bnd_lm_d[ioncnt][ib][14];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][15]-occ_vs_ion_bnd_lm_d[ioncnt][ib][15];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_l_u[ioncnt][ib][3]-occ_vs_ion_bnd_l_d[ioncnt][ib][3];
            outf << "  " << setw(w2) << occ_vs_ion_bnd_l_u[ioncnt][ib][4]-occ_vs_ion_bnd_l_d[ioncnt][ib][4];
            outf << " " << key;
            outf << endl;
          }
        }
      }
    }// if sp

    // Occupations ION:KPT:LM
    outf << endl;
    outf << "**************************************************" << endl;
    outf << endl;
    outf << "Occupations vs. ION:KPT:LM" << endl;
    outf << endl;
    key="ION:KPT:LM";
    outf << "Up spin"<< endl;
    outf << "   " << "ion   " << "type  " << "kpt   ";
    for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
      outf << LLMnames[ilm] << "  ";
    }
    outf << endl;
    ioncnt=-1;
    for(int it=0;it<ntypes;it++) {
      for(int iat=0;iat<num_each_type.at(it);iat++) {
        ioncnt++;
        for(int ik=0;ik<nkpts;ik++) {
          outf << "   " << setw(w1) << ioncnt+1;
          outf << " " << setw(w1) << it+1;
          outf << " " << setw(w1) << ik+1;
          outf << " " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][0];
          outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][1];
          outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][2];
          outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][3];
          outf << "  " << setw(w2) << occ_vs_ion_kpt_l_u[ioncnt][ik][1];
          outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][4];
          outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][5];
          outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][6];
          outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][7];
          outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][8];
          outf << "  " << setw(w2) << occ_vs_ion_kpt_l_u[ioncnt][ik][2];
          outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][9];
          outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][10];
          outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][11];
          outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][12];
          outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][13];
          outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][14];
          outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][15];
          outf << "  " << setw(w2) << occ_vs_ion_kpt_l_u[ioncnt][ik][3];
          outf << "  " << setw(w2) << occ_vs_ion_kpt_l_u[ioncnt][ik][4];
          outf << " " << key;
          outf << endl;
        }
      }
    }
    if(sp) {
      outf << "Down spin"<< endl;
      outf << "   " << "ion   " << "type  " << "kpt   ";
      for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
        outf << LLMnames[ilm] << "  ";
      }
      outf << endl;
      ioncnt=-1;
      for(int it=0;it<ntypes;it++) {
        for(int iat=0;iat<num_each_type.at(it);iat++) {
          ioncnt++;
          for(int ik=0;ik<nkpts;ik++) {
            outf << "   " << setw(w1) << ioncnt+1;
            outf << " " << setw(w1) << it+1;
            outf << " " << setw(w1) << ik+1;
            outf << " " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][0];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][1];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][2];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][3];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_l_d[ioncnt][ik][1];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][4];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][5];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][6];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][7];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][8];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_l_d[ioncnt][ik][2];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][9];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][10];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][11];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][12];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][13];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][14];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][15];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_l_d[ioncnt][ik][3];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_l_d[ioncnt][ik][4];
            outf << " " << key;
            outf << endl;
          }
        }
      }
      outf << "Up-Down spin"<< endl;
      outf << "   " << "ion   " << "type  " << "kpt   ";
      for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
        outf << LLMnames[ilm] << "  ";
      }
      outf << endl;
      ioncnt=-1;
      for(int it=0;it<ntypes;it++) {
        for(int iat=0;iat<num_each_type.at(it);iat++) {
          ioncnt++;
          for(int ik=0;ik<nkpts;ik++) {
            outf << "   " << setw(5) << ioncnt+1;
            outf << " " << setw(5) << it+1;
            outf << " " << setw(5) << ik+1;
            outf << " " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][0]-occ_vs_ion_kpt_lm_d[ioncnt][ik][0];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][1]-occ_vs_ion_kpt_lm_d[ioncnt][ik][1];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][2]-occ_vs_ion_kpt_lm_d[ioncnt][ik][2];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][3]-occ_vs_ion_kpt_lm_d[ioncnt][ik][3];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_l_u[ioncnt][ik][1]-occ_vs_ion_kpt_l_d[ioncnt][ik][1];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][4]-occ_vs_ion_kpt_lm_d[ioncnt][ik][4];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][5]-occ_vs_ion_kpt_lm_d[ioncnt][ik][5];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][6]-occ_vs_ion_kpt_lm_d[ioncnt][ik][6];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][7]-occ_vs_ion_kpt_lm_d[ioncnt][ik][7];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][8]-occ_vs_ion_kpt_lm_d[ioncnt][ik][8];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_l_u[ioncnt][ik][2]-occ_vs_ion_kpt_l_d[ioncnt][ik][2];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][9]-occ_vs_ion_kpt_lm_d[ioncnt][ik][9];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][10]-occ_vs_ion_kpt_lm_d[ioncnt][ik][10];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][11]-occ_vs_ion_kpt_lm_d[ioncnt][ik][11];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][12]-occ_vs_ion_kpt_lm_d[ioncnt][ik][12];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][13]-occ_vs_ion_kpt_lm_d[ioncnt][ik][13];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][14]-occ_vs_ion_kpt_lm_d[ioncnt][ik][14];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][15]-occ_vs_ion_kpt_lm_d[ioncnt][ik][15];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_l_u[ioncnt][ik][3]-occ_vs_ion_kpt_l_d[ioncnt][ik][3];
            outf << "  " << setw(w2) << occ_vs_ion_kpt_l_u[ioncnt][ik][4]-occ_vs_ion_kpt_l_d[ioncnt][ik][4];
            outf << " " << key;
            outf << endl;
          }
        }
      }
    }

    // Occupations ION:LM
    outf << endl;
    outf << "**************************************************" << endl;
    outf << endl;
    outf << "Occupations vs. ION:LM" << endl;
    outf << endl;
    key="ION:LM";
    outf << "Up spin"<< endl;
    outf << "   " << "ion   " << "type  ";
    for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
      outf << LLMnames[ilm] << "  ";
    }
    outf << endl;
    ioncnt=-1;
    for(int it=0;it<ntypes;it++) {
      for(int iat=0;iat<num_each_type.at(it);iat++) {
        ioncnt++;
        outf << "   " << setw(w1) << ioncnt+1;
        outf << " " << setw(w1) << it+1;
        outf << " " << setw(w2) << occ_vs_ion_lm_u[ioncnt][0];
        outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][1];
        outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][2];
        outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][3];
        outf << "  " << setw(w2) << occ_vs_ion_l_u[ioncnt][1];
        outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][4];
        outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][5];
        outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][6];
        outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][7];
        outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][8];
        outf << "  " << setw(w2) << occ_vs_ion_l_u[ioncnt][2];
        outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][9];
        outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][10];
        outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][11];
        outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][12];
        outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][13];
        outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][14];
        outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][15];
        outf << "  " << setw(w2) << occ_vs_ion_l_u[ioncnt][3];
        outf << "  " << setw(w2) << occ_vs_ion_l_u[ioncnt][4];
        outf << " " << key;
        outf << endl;
      }
    }
    if(sp) {
      outf << "Down spin"<< endl;
      outf << "   " << "ion   " << "type  ";
      for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
        outf << LLMnames[ilm] << "  ";
      }
      outf << endl;
      ioncnt=-1;
      for(int it=0;it<ntypes;it++) {
        for(int iat=0;iat<num_each_type.at(it);iat++) {
          ioncnt++;
          outf << "   " << setw(w1) << ioncnt+1;
          outf << " " << setw(w1) << it+1;
          outf << " " << setw(w2) << occ_vs_ion_lm_d[ioncnt][0];
          outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][1];
          outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][2];
          outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][3];
          outf << "  " << setw(w2) << occ_vs_ion_l_d[ioncnt][1];
          outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][4];
          outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][5];
          outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][6];
          outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][7];
          outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][8];
          outf << "  " << setw(w2) << occ_vs_ion_l_d[ioncnt][2];
          outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][9];
          outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][10];
          outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][11];
          outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][12];
          outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][13];
          outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][14];
          outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][15];
          outf << "  " << setw(w2) << occ_vs_ion_l_d[ioncnt][3];
          outf << "  " << setw(w2) << occ_vs_ion_l_d[ioncnt][4];
          outf << " " << key;
          outf << endl;
        }
      }
      outf << "Up-Down spin"<< endl;
      outf << "   " << "ion   " << "type  ";
      for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
        outf << LLMnames[ilm] << "  ";
      }
      outf << endl;
      ioncnt=-1;
      for(int it=0;it<ntypes;it++) {
        for(int iat=0;iat<num_each_type.at(it);iat++) {
          ioncnt++;
          outf << "   " << setw(5) << ioncnt+1;
          outf << " " << setw(5) << it+1;
          outf << " " << setw(w2) << occ_vs_ion_lm_u[ioncnt][0]-occ_vs_ion_lm_d[ioncnt][0];
          outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][1]-occ_vs_ion_lm_d[ioncnt][1];
          outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][2]-occ_vs_ion_lm_d[ioncnt][2];
          outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][3]-occ_vs_ion_lm_d[ioncnt][3];
          outf << "  " << setw(w2) << occ_vs_ion_l_u[ioncnt][1]-occ_vs_ion_l_d[ioncnt][1];
          outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][4]-occ_vs_ion_lm_d[ioncnt][4];
          outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][5]-occ_vs_ion_lm_d[ioncnt][5];
          outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][6]-occ_vs_ion_lm_d[ioncnt][6];
          outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][7]-occ_vs_ion_lm_d[ioncnt][7];
          outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][8]-occ_vs_ion_lm_d[ioncnt][8];
          outf << "  " << setw(w2) << occ_vs_ion_l_u[ioncnt][2]-occ_vs_ion_l_d[ioncnt][2];
          outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][9]-occ_vs_ion_lm_d[ioncnt][9];
          outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][10]-occ_vs_ion_lm_d[ioncnt][10];
          outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][11]-occ_vs_ion_lm_d[ioncnt][11];
          outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][12]-occ_vs_ion_lm_d[ioncnt][12];
          outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][13]-occ_vs_ion_lm_d[ioncnt][13];
          outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][14]-occ_vs_ion_lm_d[ioncnt][14];
          outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][15]-occ_vs_ion_lm_d[ioncnt][15];
          outf << "  " << setw(w2) << occ_vs_ion_l_u[ioncnt][3]-occ_vs_ion_l_d[ioncnt][3];
          outf << "  " << setw(w2) << occ_vs_ion_l_u[ioncnt][4]-occ_vs_ion_l_d[ioncnt][4];
          outf << " " << key;
          outf << endl;
        }
      }
    }// sp
  }
}

// ***************************************************************************
// PROJFUNCS PROJFUNCS PROJFUNCS PROJFUNCS PROJFUNCS PROJFUNCS PROJFUNCS PROJF
// ***************************************************************************

// ***************************************************************************
// ProcessProjection
// ***************************************************************************
// The projections might be a complex amplitude, a true
// real probability, or a probability with a complex
// phase.  Depending on which, you need to process
// the projection differently.  Hence this routine, so
// all processing can be changed at once.
namespace pflow {
  std::complex<double> ProcessProjection(const std::complex<double>& proj) {
    std::complex<double> p;
    p=proj;
    // For proj=probability with phase.
    // This allows for only positive probabilities (does not exactly match vasp output).
    // return sqrt(p*conj(p));
    // For proj=real probability with complex part always zero.
    // This allows for negative probabilities (does exactly match vasp output).
    return p;
  }
}
// ***************************************************************************
// ReadInProj
// ***************************************************************************
// Need to read everything in twice if it is spin polarized.
// rspin=1 for spin polarized, rspin=2 for non-spin polarized.
namespace pflow {
  void ReadInProj(projdata& pd) {
    int have_aug;
    string s;
    char c;
    ifstream infile(pd.PROOUTinfile.c_str());
    aurostd::InFileExistCheck("ReadInProj",pd.PROOUTinfile.c_str(),infile);
    string sdum;
    // Initial data
    infile >> sdum; // Title
    infile >> sdum >> sdum >> sdum >> pd.nkpts;
    infile >> sdum >> sdum >> sdum >> pd.nbands;
    infile >> sdum >> sdum >> sdum >> pd.nions;
    infile >> pd.ntypes >> pd.ntypes;
    pd.num_each_type=vector<int> (pd.ntypes);
    for(int it=0;it<pd.ntypes;it++) {
      infile >> pd.num_each_type.at(it);
    }
    // For getting kpoint energies.
    int skip=2;
    if(pd.nions==1) skip=1;

    // Fermi weights (kpts (inner loop) and bands (outer loop)).
    pd.wfermi_u = aurostd::matrix<double> (pd.nbands,pd.nkpts,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
    for(int ib=0;ib<pd.nbands;ib++) {
      for(int ik=0;ik<pd.nkpts;ik++) {
        infile >> pd.wfermi_u[ib][ik];
        // tpx
        // cout<<"ib,ik,wfermi_u "<<ib<<" "<<ik<<" "<<pd.wfermi_u[ib][ik]<<endl;
      }
    }

    // Check to see if there are spd or spdf electrons (nl=3 or nl=4).
    if(aurostd::FindIfStringInStream("s      p      d      f",infile)) {
      pd.nl=4;
    }
    else {
      pd.nl=3;
    }
    if(pd.nl==3) {
      pd.nlm=9; // for s,p,d orbitals
      pd.nlmtot=12; // for s,p,d orbitals + p,d,all totals
    }

    // Projections
    aurostd::matrix<std::complex<double> > m(pd.nbands,pd.nlm,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<aurostd::matrix<std::complex<double> > > mm(pd.nkpts,pd.nions,m);  //CO20200404 pflow::matrix()->aurostd::matrix()
    pd.pdat_u=vector<aurostd::matrix<aurostd::matrix<std::complex<double> > > > (pd.ntypes,mm);  //CO20200404 pflow::matrix()->aurostd::matrix()
    vector<double> rval(pd.nlm);
    vector<double> ival(pd.nlm);
    for(int it=0;it<pd.ntypes;it++) {
      for(int ik=0;ik<pd.nkpts;ik++) {
        for(int iat=0;iat<pd.num_each_type.at(it);iat++) {
          for(int ib=0;ib<pd.nbands;ib++) {
            for(int ilm=0;ilm<pd.nlm;ilm++) {
              infile >> rval[ilm] >> ival[ilm];
              std::complex<double> ctmp (rval[ilm],ival[ilm]);
              pd.pdat_u[it][ik][iat][ib][ilm]=ctmp;
              // tpx
              //  cout << "PROJ ib,ilm " << ib << " " << ilm << " " <<  rval[ilm] << " " <<  ival[ilm] << " " << pd.pdat_u[it][ik][iat][ib][ilm] << endl;
            }
          }
        }
      }
    }

    // Augmented Projections
    getline(infile,s); // Gets final carriage return from previous line.
    getline(infile,s); // Gets "augmentation part"
    c = infile.peek();
    if(c=='#') {
      have_aug=0; // There are no augmented projections.
    }
    else {
      have_aug=1; // There are augmented projections.
    }
    if(have_aug) {
      for(int ik=0;ik<pd.nkpts;ik++) {
        for(int ib=0;ib<pd.nbands;ib++) {
          for(int it=0;it<pd.ntypes;it++) {
            for(int iat=0;iat<pd.num_each_type.at(it);iat++) {
              for(int ilm=0;ilm<pd.nlm;ilm++) {
                infile >> rval[ilm];
                //tpx
                // cout << "Rval " << rval[ilm] << " " << ilm  << endl;
              }
              for(int ilm=0;ilm<pd.nlm;ilm++) {
                pd.pdat_u[it][ik][iat][ib][ilm]=
                  pd.pdat_u[it][ik][iat][ib][ilm]+rval[ilm];
              }
            }
          }
        }
      }
      getline(infile,s); // Gets final carriage return last line of augmentation charges.
    }//if have_aug

    // Kpt weights and values and energy vs. kpt,bnd
    pd.wkpt = vector<double> (pd.nkpts,0.0);
    pd.kpts = aurostd::matrix<double> (pd.nkpts,3,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
    pd.ener_k_b_u = aurostd::matrix<double> (pd.nkpts,pd.nbands,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
    getline(infile,s); // Gets line "# of k-points: ..."
    for(int ik=0;ik<pd.nkpts;ik++) {
      getline(infile,s); // Gets blank line.
      infile >>sdum>>sdum>>sdum>>pd.kpts[ik][0]>>pd.kpts[ik][1]>>pd.kpts[ik][2]>>sdum>>sdum>>pd.wkpt[ik];
      getline(infile,s); // Gets final carriage return from previous line.
      getline(infile,s); // Gets line.
      for(int ib=0;ib<pd.nbands;ib++) {
        infile >>sdum>>sdum>>sdum>>sdum>>pd.ener_k_b_u[ik][ib]>>sdum>>sdum>>sdum;
        for(int ii=0;ii<pd.nions+skip;ii++) {
          // Gets lines (nbands*(nions+skip)).  Note that
          // the numer of blank lines varies in different vasp versions
          // so this loop does not count blank lines.
          getline(infile,s);
          if(aurostd::RemoveSpaces(s).size()==0) {
            ii=ii-1; // Don't count blank lines.
          }
        }
      }
    }

    // Do again if spin polarized.  Look for another PROOUT title line.
    // At this point lines could be blank, end of file, or PROOUT.  We want
    // to keep reading if they are blank, stop if they are EOF, and do
    // spin-polarized if they are PROOUT.
    infile >> s; // Skip all whitespace and get PROOUT title if it exists.
    if(aurostd::RemoveSpaces(s)=="PROOUT") {
      // Set spin polarization dependent variables.
      pd.sp=1;
      pd.rspin=1;

      getline(infile,s); // Gets endline after PROOUT.
      getline(infile,s); // get line.
      getline(infile,s); // get line.

      // Fermi weights (kpts (inner loop) and bands (outer loop)).
      pd.wfermi_d = aurostd::matrix<double> (pd.nbands,pd.nkpts,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
      for(int ib=0;ib<pd.nbands;ib++) {
        for(int ik=0;ik<pd.nkpts;ik++) {
          infile >> pd.wfermi_d[ib][ik];
        }
      }

      // Projections
      aurostd::matrix<std::complex<double> > m(pd.nbands,pd.nlm,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
      aurostd::matrix<aurostd::matrix<std::complex<double> > > mm(pd.nkpts,pd.nions,m);  //CO20200404 pflow::matrix()->aurostd::matrix()
      pd.pdat_d=vector<aurostd::matrix<aurostd::matrix<std::complex<double> > > > (pd.ntypes,mm);  //CO20200404 pflow::matrix()->aurostd::matrix()
      vector<double> rval(pd.nlm);
      vector<double> ival(pd.nlm);
      for(int it=0;it<pd.ntypes;it++) {
        for(int ik=0;ik<pd.nkpts;ik++) {
          for(int iat=0;iat<pd.num_each_type.at(it);iat++) {
            for(int ib=0;ib<pd.nbands;ib++) {
              for(int ilm=0;ilm<pd.nlm;ilm++) {
                infile >> rval[ilm] >> ival[ilm];
              }
              for(int ilm=0;ilm<pd.nlm;ilm++) {
                std::complex<double> ctmp (rval[ilm],ival[ilm]);
                pd.pdat_d[it][ik][iat][ib][ilm]=ctmp;
              }
            }
          }
        }
      }

      // Augmented Projections
      string s;
      getline(infile,s); // Gets final carriage return from previous line.
      getline(infile,s); // Gets "augmentation part"
      c = infile.peek();
      if(c=='#') {
        have_aug=0; // There are no augmented projections.
      }
      else {
        have_aug=1; // There are augmented projections.
      }
      if(have_aug) {
        for(int ik=0;ik<pd.nkpts;ik++) {
          for(int ib=0;ib<pd.nbands;ib++) {
            for(int it=0;it<pd.ntypes;it++) {
              for(int iat=0;iat<pd.num_each_type.at(it);iat++) {
                for(int ilm=0;ilm<pd.nlm;ilm++) {
                  infile >> rval[ilm];
                }
                for(int ilm=0;ilm<pd.nlm;ilm++) {
                  pd.pdat_d[it][ik][iat][ib][ilm]=
                    pd.pdat_d[it][ik][iat][ib][ilm]+rval[ilm];
                }
              }
            }
          }
        }
        getline(infile,s); // Gets final carriage return last line of augmentation charges.
      }//if have_aug

      // Kpt weights and values and energy vs. kpt,bnd
      pd.ener_k_b_d = aurostd::matrix<double> (pd.nkpts,pd.nbands,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
      getline(infile,s); // Gets line "# of k-points: ..."
      for(int ik=0;ik<pd.nkpts;ik++) {
        getline(infile,s); // Gets blank line.
        infile >>sdum>>sdum>>sdum>>pd.kpts[ik][0]>>pd.kpts[ik][1]>>pd.kpts[ik][2]>>sdum>>sdum>>pd.wkpt[ik];
        getline(infile,s); // Gets final carriage return from previous line.
        getline(infile,s); // Gets line.
        for(int ib=0;ib<pd.nbands;ib++) {
          infile >>sdum>>sdum>>sdum>>sdum>>pd.ener_k_b_d[ik][ib]>>sdum>>sdum>>sdum;
          for(int ii=0;ii<pd.nions+skip;ii++) {
            // Gets lines (nbands*(nions+skip)).  Note that
            // the numer of blank lines varies in different vasp versions
            // so this loop does not count blank lines.
            getline(infile,s);
            if(aurostd::RemoveSpaces(s).size()==0) {
              ii=ii-1; // Don't count blank lines.
            }
          }
        }
      }

    }// if spin polarized

  }// end
}

// ***************************************************************************
// CalcNeatProj
// ***************************************************************************
// Calculates all the basic occupancies.  If only_occ=true
// then it does this only for occupied states.  If only_occ=false
// then it does this for all states, which does not really
// give occupations, but sets things up for PDOS calculations.
namespace pflow {
  void CalcNeatProj(projdata& pd, int only_occ) {
    int ioncnt;

    // Get occ_vs_ion_kpt_bnd_lm_(u,d)
    // Get occ_vs_ion_kpt_lm_(u,d)
    // Get occ_vs_ion_bnd_lm_(u,d)
    // Get occ_vs_ion_lm_(u,d)
    aurostd::matrix<double> m1 (pd.nbands,pd.nlm_max,0.0); //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> m2 (pd.nkpts,pd.nlm_max,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> m3 (pd.nbands,pd.nlm_max,0.0); //CO20200404 pflow::matrix()->aurostd::matrix()
    pd.occ_vs_ion_kpt_bnd_lm_u = aurostd::matrix<aurostd::matrix<double> > (pd.nions,pd.nkpts,m1); //CO20200404 pflow::matrix()->aurostd::matrix()
    pd.occ_vs_ion_kpt_lm_u = vector<aurostd::matrix<double> > (pd.nions,m2); //CO20200404 pflow::matrix()->aurostd::matrix()
    pd.occ_vs_ion_bnd_lm_u = vector<aurostd::matrix<double> > (pd.nions,m3); //CO20200404 pflow::matrix()->aurostd::matrix()
    pd.occ_vs_ion_lm_u = aurostd::matrix<double> (pd.nions,pd.nlm_max,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
    if(pd.sp) {
      pd.occ_vs_ion_kpt_bnd_lm_d = aurostd::matrix<aurostd::matrix<double> > (pd.nions,pd.nkpts,m1); //CO20200404 pflow::matrix()->aurostd::matrix()
      pd.occ_vs_ion_kpt_lm_d = vector<aurostd::matrix<double> > (pd.nions,m2); //CO20200404 pflow::matrix()->aurostd::matrix()
      pd.occ_vs_ion_bnd_lm_d = vector<aurostd::matrix<double> > (pd.nions,m3); //CO20200404 pflow::matrix()->aurostd::matrix()
      pd.occ_vs_ion_lm_d = aurostd::matrix<double> (pd.nions,pd.nlm_max,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
    }
    ioncnt=-1;
    for(int it=0;it<pd.ntypes;it++) {
      for(int iat=0;iat<pd.num_each_type.at(it);iat++) {
        ioncnt++;
        for(int ik=0;ik<pd.nkpts;ik++) {
          for(int ib=0;ib<pd.nbands;ib++) {
            for(int ilm=0;ilm<pd.nlm;ilm++) {
              double temp;
              temp=ProcessProjection(pd.pdat_u[it][ik][iat][ib][ilm]).real();
              if(only_occ) {
                pd.occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][ilm]=pd.rspin*pd.wfermi_u[ib][ik]*(temp);
              }
              else {
                pd.occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][ilm]=pd.rspin*(temp);
              }
              pd.occ_vs_ion_kpt_lm_u[ioncnt][ik][ilm]=
                pd.occ_vs_ion_kpt_lm_u[ioncnt][ik][ilm]+
                pd.occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][ilm];
              pd.occ_vs_ion_bnd_lm_u[ioncnt][ib][ilm]=
                pd.occ_vs_ion_bnd_lm_u[ioncnt][ib][ilm]+
                pd.wkpt[ik]*pd.occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][ilm];
              pd.occ_vs_ion_lm_u[ioncnt][ilm]=
                pd.occ_vs_ion_lm_u[ioncnt][ilm]+
                pd.wkpt[ik]*pd.occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][ilm];
              if(pd.sp) {
                temp=ProcessProjection(pd.pdat_d[it][ik][iat][ib][ilm]).real();
                if(only_occ) {
                  pd.occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][ilm]=pd.rspin*pd.wfermi_d[ib][ik]*(temp);
                }
                else {
                  pd.occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][ilm]=pd.rspin*(temp);
                }
                pd.occ_vs_ion_kpt_lm_d[ioncnt][ik][ilm]=
                  pd.occ_vs_ion_kpt_lm_d[ioncnt][ik][ilm]+
                  pd.occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][ilm];
                pd.occ_vs_ion_bnd_lm_d[ioncnt][ib][ilm]=
                  pd.occ_vs_ion_bnd_lm_d[ioncnt][ib][ilm]+
                  pd.wkpt[ik]*pd.occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][ilm];
                pd.occ_vs_ion_lm_d[ioncnt][ilm]=
                  pd.occ_vs_ion_lm_d[ioncnt][ilm]+
                  pd.wkpt[ik]*pd.occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][ilm];
              }
            }
          }
        }
      }
    }

    // Get sums over m.
    // Get occ_vs_ion_l_(u,d)
    // Get occ_vs_ion_kpt_l_(u,d)
    // Get occ_vs_ion_bnd_l_(u,d)
    // Get occ_vs_ion_kpt_bnd_l_(u,d)
    // This assumes the usual spd orbitals.
    aurostd::matrix<double> n1 (pd.nbands,pd.nl_max+1,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> n2 (pd.nkpts,pd.nl_max+1,0.0); //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> n3 (pd.nbands,pd.nl_max+1,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
    pd.occ_vs_ion_kpt_bnd_l_u = aurostd::matrix<aurostd::matrix<double> > (pd.nions,pd.nkpts,n1);  //CO20200404 pflow::matrix()->aurostd::matrix()
    pd.occ_vs_ion_kpt_l_u = vector<aurostd::matrix<double> > (pd.nions,n2);  //CO20200404 pflow::matrix()->aurostd::matrix()
    pd.occ_vs_ion_bnd_l_u = vector<aurostd::matrix<double> > (pd.nions,n3);  //CO20200404 pflow::matrix()->aurostd::matrix()
    pd.occ_vs_ion_l_u = aurostd::matrix<double> (pd.nions,pd.nl_max+1,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
    if(pd.sp) {
      pd.occ_vs_ion_kpt_bnd_l_d = aurostd::matrix<aurostd::matrix<double> > (pd.nions,pd.nkpts,n1);  //CO20200404 pflow::matrix()->aurostd::matrix()
      pd.occ_vs_ion_kpt_l_d = vector<aurostd::matrix<double> > (pd.nions,n2);  //CO20200404 pflow::matrix()->aurostd::matrix()
      pd.occ_vs_ion_bnd_l_d = vector<aurostd::matrix<double> > (pd.nions,n3);  //CO20200404 pflow::matrix()->aurostd::matrix()
      pd.occ_vs_ion_l_d = aurostd::matrix<double> (pd.nions,pd.nl_max+1,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
    }

    ioncnt=-1;
    for(int it=0;it<pd.ntypes;it++) {
      for(int iat=0;iat<pd.num_each_type.at(it);iat++) {
        ioncnt++;
        for(int ik=0;ik<pd.nkpts;ik++) {
          for(int ib=0;ib<pd.nbands;ib++) {

            // S
            pd.occ_vs_ion_l_u[ioncnt][0]+=pd.occ_vs_ion_lm_u[ioncnt][0];
            pd.occ_vs_ion_kpt_l_u[ioncnt][ik][0]+=pd.occ_vs_ion_kpt_lm_u[ioncnt][ik][0];
            pd.occ_vs_ion_bnd_l_u[ioncnt][ib][0]+=pd.occ_vs_ion_bnd_lm_u[ioncnt][ib][0];
            pd.occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][0]+=pd.occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][0];
            if(pd.sp) {
              pd.occ_vs_ion_l_d[ioncnt][0]+=pd.occ_vs_ion_lm_d[ioncnt][0];
              pd.occ_vs_ion_kpt_l_d[ioncnt][ik][0]+=pd.occ_vs_ion_kpt_lm_d[ioncnt][ik][0];
              pd.occ_vs_ion_bnd_l_d[ioncnt][ib][0]+=pd.occ_vs_ion_bnd_lm_d[ioncnt][ib][0];
              pd.occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][0]+=pd.occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][0];
            }

            // P
            for(int ip=0;ip<3;ip++) {
              pd.occ_vs_ion_l_u[ioncnt][1]+=pd.occ_vs_ion_lm_u[ioncnt][1+ip];
              pd.occ_vs_ion_kpt_l_u[ioncnt][ik][1]+=pd.occ_vs_ion_kpt_lm_u[ioncnt][ik][1+ip];
              pd.occ_vs_ion_bnd_l_u[ioncnt][ib][1]+=pd.occ_vs_ion_bnd_lm_u[ioncnt][ib][1+ip];
              pd.occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][1]+=pd.occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][1+ip];
              if(pd.sp) {
                pd.occ_vs_ion_l_d[ioncnt][1]+=pd.occ_vs_ion_lm_d[ioncnt][1+ip];
                pd.occ_vs_ion_kpt_l_d[ioncnt][ik][1]+=pd.occ_vs_ion_kpt_lm_d[ioncnt][ik][1+ip];
                pd.occ_vs_ion_bnd_l_d[ioncnt][ib][1]+=pd.occ_vs_ion_bnd_lm_d[ioncnt][ib][1+ip];
                pd.occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][1]+=pd.occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][1+ip];
              }
            }

            // D
            for(int id=0;id<5;id++) {
              pd.occ_vs_ion_l_u[ioncnt][2]+=pd.occ_vs_ion_lm_u[ioncnt][4+id];
              pd.occ_vs_ion_kpt_l_u[ioncnt][ik][2]+=pd.occ_vs_ion_kpt_lm_u[ioncnt][ik][4+id];
              pd.occ_vs_ion_bnd_l_u[ioncnt][ib][2]+=pd.occ_vs_ion_bnd_lm_u[ioncnt][ib][4+id];
              pd.occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][2]+=pd.occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][4+id];
              if(pd.sp) {
                pd.occ_vs_ion_l_d[ioncnt][2]+=pd.occ_vs_ion_lm_d[ioncnt][4+id];
                pd.occ_vs_ion_kpt_l_d[ioncnt][ik][2]+=pd.occ_vs_ion_kpt_lm_d[ioncnt][ik][4+id];
                pd.occ_vs_ion_bnd_l_d[ioncnt][ib][2]+=pd.occ_vs_ion_bnd_lm_d[ioncnt][ib][4+id];
                pd.occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][2]+=pd.occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][4+id];
              }
            }

            // F
            for(int id=0;id<7;id++) {
              pd.occ_vs_ion_l_u[ioncnt][3]+=pd.occ_vs_ion_lm_u[ioncnt][9+id];
              pd.occ_vs_ion_kpt_l_u[ioncnt][ik][3]+=pd.occ_vs_ion_kpt_lm_u[ioncnt][ik][9+id];
              pd.occ_vs_ion_bnd_l_u[ioncnt][ib][3]+=pd.occ_vs_ion_bnd_lm_u[ioncnt][ib][9+id];
              pd.occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][3]+=pd.occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][9+id];
              if(pd.sp) {
                pd.occ_vs_ion_l_d[ioncnt][3]+=pd.occ_vs_ion_lm_d[ioncnt][9+id];
                pd.occ_vs_ion_kpt_l_d[ioncnt][ik][3]+=pd.occ_vs_ion_kpt_lm_d[ioncnt][ik][9+id];
                pd.occ_vs_ion_bnd_l_d[ioncnt][ib][3]+=pd.occ_vs_ion_bnd_lm_d[ioncnt][ib][9+id];
                pd.occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][3]+=pd.occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][9+id];
              }
            }

          }//ib
        }//ik

        // Normalize sums over m for multiple counting.
        for(int il=0;il<pd.nl;il++) {
          pd.occ_vs_ion_l_u[ioncnt][il]/=(pd.nbands*pd.nkpts);
          if(pd.sp)  pd.occ_vs_ion_l_d[ioncnt][il]/=(pd.nbands*pd.nkpts);
          for(int ik=0;ik<pd.nkpts;ik++) {
            pd.occ_vs_ion_kpt_l_u[ioncnt][ik][il]/=(pd.nbands);

            if(pd.sp) pd.occ_vs_ion_kpt_l_d[ioncnt][ik][il]/=(pd.nbands);
          }
          for(int ib=0;ib<pd.nbands;ib++) {
            pd.occ_vs_ion_bnd_l_u[ioncnt][ib][il]/=(pd.nkpts);
            if(pd.sp) pd.occ_vs_ion_bnd_l_d[ioncnt][ib][il]/=(pd.nkpts);
          }
        }

        // Calculate total sums over L.

        int iid=pd.nl_max;
        for(int id=0;id<pd.nl_max;id++) {
          pd.occ_vs_ion_l_u[ioncnt][iid]+=pd.occ_vs_ion_l_u[ioncnt][id];	
          if(pd.sp) pd.occ_vs_ion_l_d[ioncnt][iid]+=pd.occ_vs_ion_l_d[ioncnt][id];
          for(int ik=0;ik<pd.nkpts;ik++) {
            pd.occ_vs_ion_kpt_l_u[ioncnt][ik][iid]+=pd.occ_vs_ion_kpt_l_u[ioncnt][ik][id];
            if(pd.sp) pd.occ_vs_ion_kpt_l_d[ioncnt][ik][iid]+=pd.occ_vs_ion_kpt_l_d[ioncnt][ik][id];
          }
          for(int ib=0;ib<pd.nbands;ib++) {
            pd.occ_vs_ion_bnd_l_u[ioncnt][ib][iid]+=pd.occ_vs_ion_bnd_l_u[ioncnt][ib][id];
            if(pd.sp) pd.occ_vs_ion_bnd_l_d[ioncnt][ib][iid]+=pd.occ_vs_ion_bnd_l_d[ioncnt][ib][id];
          }
          for(int ik=0;ik<pd.nkpts;ik++) {
            for(int ib=0;ib<pd.nbands;ib++) {
              pd.occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][iid]+=pd.occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][id];
              if(pd.sp) pd.occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][iid]+=pd.occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][id];
            }
          }
        }

      }// at
    }// type

    // Create single occ variable for all the l,m and sums over m together (for lmtot).
    // This way you can just loop over lm from 0 to lmtot-1.
    m1 = aurostd::matrix<double> (pd.nbands,pd.nlmtot,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
    pd.occ_vs_ion_kpt_bnd_lmtot_u = aurostd::matrix<aurostd::matrix<double> > (pd.nions,pd.nkpts,m1);  //CO20200404 pflow::matrix()->aurostd::matrix()
    if(pd.sp) pd.occ_vs_ion_kpt_bnd_lmtot_d = aurostd::matrix<aurostd::matrix<double> > (pd.nions,pd.nkpts,m1);  //CO20200404 pflow::matrix()->aurostd::matrix()
    for(int ii=0;ii<pd.nions;ii++) {
      for(int ik=0;ik<pd.nkpts;ik++) {
        for(int ib=0;ib<pd.nbands;ib++) {
          pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][0]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][0];
          pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][1]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][1];
          pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][2]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][2];
          pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][3]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][3];
          pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][4]=pd.occ_vs_ion_kpt_bnd_l_u[ii][ik][ib][1];
          pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][5]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][4];
          pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][6]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][5];
          pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][7]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][6];
          pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][8]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][7];
          pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][9]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][8];
          pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][10]=pd.occ_vs_ion_kpt_bnd_l_u[ii][ik][ib][2];
          pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][11]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][9];
          pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][12]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][10];
          pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][13]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][11];
          pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][14]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][12];
          pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][15]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][13];
          pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][16]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][14];
          pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][17]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][15];
          pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][18]=pd.occ_vs_ion_kpt_bnd_l_u[ii][ik][ib][3];
          pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][19]=pd.occ_vs_ion_kpt_bnd_l_u[ii][ik][ib][4];
          if(pd.sp) {
            pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][0]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][0];
            pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][1]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][1];
            pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][2]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][2];
            pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][3]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][3];
            pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][4]=pd.occ_vs_ion_kpt_bnd_l_d[ii][ik][ib][1];
            pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][5]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][4];
            pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][6]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][5];
            pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][7]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][6];
            pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][8]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][7];
            pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][9]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][8];
            pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][10]=pd.occ_vs_ion_kpt_bnd_l_d[ii][ik][ib][2];
            pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][11]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][9];
            pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][12]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][10];
            pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][13]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][11];
            pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][14]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][12];
            pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][15]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][13];
            pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][16]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][14];
            pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][17]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][15];
            pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][18]=pd.occ_vs_ion_kpt_bnd_l_d[ii][ik][ib][3];
            pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][19]=pd.occ_vs_ion_kpt_bnd_l_d[ii][ik][ib][4];
          }
        }
      }
    }
  }
}

// ***************************************************************************
// ReadInPDOSData
// ***************************************************************************
// Reads in parameters that determine what partial DOS
//   will be calculated.
namespace pflow {
  void ReadInPDOSData(const pflow::projdata& prd, pflow::pdosdata& pdd) {
    stringstream message;

    ifstream infile(pdd.PDOSinfile.c_str());
    aurostd::InFileExistCheck("ReadInPDOSData",pdd.PDOSinfile.c_str(),infile);

    // Defaults (emin,emax,nbins,smooth_sigma,print_params)
    double emn=prd.ener_k_b_u[0][0];
    double emx=emn;
    double e;
    for(int ik=0;ik<prd.nkpts;ik++) {
      for(int ib=0;ib<prd.nbands;ib++) {
        e = prd.ener_k_b_u[ik][ib];
        if(emn>e) emn=e;
        if(emx<e) emx=e;
        if(prd.sp) e = prd.ener_k_b_d[ik][ib];
        if(emn>e) emn=e;
        if(emx<e) emx=e;
      }
    }
    pdd.emin=emn-0.5;
    pdd.emax=emx+0.5;
    pdd.nbins=300;
    pdd.smooth_sigma=(pdd.emax-pdd.emin)/pdd.nbins;
    pdd.print_params=0;

    // Read in all the tokens.

    string s,s_ns;
    vector<string> token_vec;
    vector<string> val_vec;
    while (!infile.eof()) {
      getline(infile,s);
      s_ns=aurostd::RemoveSpaces(s); // Get string with no spaces.
      if(s_ns.size()>0) { // Make sure line was not blank
        if(s_ns[0]!='#') { // Exclude comment lines
          string token;
          string sval;
          int id=s.find('=');
          if(id>=(int)s.length()) {
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"the following token is incorrectly formatted: "+s,_INPUT_ILLEGAL_); //CO20200624
          }
          token=s.substr(0,id);
          token=aurostd::RemoveSpaces(token);
          int i_f=std::min((int)s.length(),(int)s.find('#')); // End of string
          sval=s.substr(id+1,i_f-id-1);
          token_vec.push_back(token);
          val_vec.push_back(sval);
        } //if s_ns[0]!=#
      } // if s_ns.size>0

    } // while !infile.eof

    // Read in values for variables in the pdosdata object.
    vector<int> tmp;
    int atom_cnt=0;
    for(int i=0;i<(int)token_vec.size();i++) {
      int id=0;
      string tok=token_vec[i];
      string sval;
      // double dval;
      int ival;
      int found_token=0;
      if(tok=="ATOMS") {
        id=0;
        atom_cnt++;
        tmp.clear();
        // set up atoms matrix
        while (id<(int)val_vec[i].size()) {
          string s;
          s=aurostd::GetNextVal(val_vec[i],id);
          if(s.size()>0) {
            ival=atoi(s.c_str());
            tmp.push_back(ival);
          }
        }
        pdd.pdos_at.push_back(tmp);
        // add vector to pdos_k,pdos_b,pdos_lm
        tmp.clear();
        pdd.pdos_k.push_back(tmp);
        pdd.pdos_b.push_back(tmp);
        pdd.pdos_lm.push_back(tmp);
        found_token=1;
      }// ATOMS
      if(tok=="KPOINTS") {
        if((int)pdd.pdos_k.size()!=atom_cnt || atom_cnt==0) AtomCntError(tok,(int)pdd.pdos_k.size(),atom_cnt);
        id=0;
        tmp.clear();
        // set up atoms matrix
        while (id<(int)val_vec[i].size()) {
          string s;
          s=aurostd::GetNextVal(val_vec[i],id);
          if(s.size()>0) {
            ival=atoi(s.c_str());
            pdd.pdos_k[pdd.pdos_k.size()-1].push_back(ival);
          }
        }
        found_token=1;
      }// KPOINTS
      if(tok=="BANDS") {
        id=0;
        tmp.clear();
        // set up atoms matrix
        while (id<(int)val_vec[i].size()) {
          string s;
          s=aurostd::GetNextVal(val_vec[i],id);
          if(s.size()>0) {
            ival=atoi(s.c_str());
            pdd.pdos_b[pdd.pdos_b.size()-1].push_back(ival);
          }
        }
        found_token=1;
      }// LMVALUES
      if(tok=="LMVALUES") {
        if((int)pdd.pdos_lm.size()!=atom_cnt || atom_cnt==0) AtomCntError(tok,(int)pdd.pdos_lm.size(),atom_cnt);
        id=0;
        tmp.clear();
        // set up atoms matrix
        while (id<(int)val_vec[i].size()) {
          string s;
          s=aurostd::GetNextVal(val_vec[i],id);
          if(s.size()>0) {
            ival=atoi(s.c_str());
            pdd.pdos_lm[pdd.pdos_lm.size()-1].push_back(ival);
          }
        }
        found_token=1;
      }// LMVALUES
      if(tok=="EMIN") {
        id=0;
        pdd.emin=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
        found_token=1;
      }// EMIN
      if(tok=="EMAX") {
        id=0;
        pdd.emax=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
        found_token=1;
      }// EMAX
      if(tok=="NBINS") {
        id=0;
        pdd.nbins=atoi(aurostd::GetNextVal(val_vec[i],id).c_str());
        found_token=1;
      }// NBIN
      if(tok=="SMOOTH_SIGMA") {
        id=0;
        pdd.smooth_sigma=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
        found_token=1;
      }// SMOOTH_SIGMA
      if(tok=="PRINT_PARAMS") {
        id=0;
        pdd.print_params=atoi(aurostd::GetNextVal(val_vec[i],id).c_str());
        found_token=1;
      }// PRINT_PARAMS

      if(!found_token) {
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"token is not recognized: "+tok,_INPUT_ILLEGAL_); //CO20200624
      }

    } // for i

    // Defaults (all atoms, all kpts, all bands, lm = nlmtot = all l+m)
    vector<int> vv(4);
    vv[0]=prd.nions;
    vv[1]=prd.nkpts;
    vv[2]=prd.nbands;
    vv[3]=prd.nlmtot;
    int mx = *max_element(vv.begin(),vv.end());
    vector<int> cnt(mx);
    for(int i=1;i<=mx;i++) {
      cnt[i-1]=i;
    }
    //  int mx = max(prd.nions,prd.nkpts,prd.nbands,prd.nlmtot);
    vector<int> atv(prd.nions);
    for(int i=0;i<(int)pdd.pdos_at.size();i++) {
      if(pdd.pdos_at[i].size()==0) pdd.pdos_at[i] = vector<int> (cnt.begin(),cnt.begin()+prd.nions);
      if(pdd.pdos_k[i].size()==0) pdd.pdos_k[i] = vector<int> (cnt.begin(),cnt.begin()+prd.nkpts);
      if(pdd.pdos_b[i].size()==0) pdd.pdos_b[i] = vector<int> (cnt.begin(),cnt.begin()+prd.nbands);
      if(pdd.pdos_lm[i].size()==0) pdd.pdos_lm[i] = vector<int> (1,prd.nlmtot);
    }

    // Check that all data is within bounds
    for(int i=0;i<(int)pdd.pdos_at.size();i++) { // cases
      for(int ia=0;ia<(int)pdd.pdos_at[i].size();ia++) {
        int atp=pdd.pdos_at[i][ia];
        if(atp<1 || atp>prd.nions) {
          message << "Error in case: " << i+1 << endl;
          message << "too low or high an atom number has been entered in entry: " << ia+1 << endl;
          message << "Bad value is: " << atp << endl;
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_); //CO20200624
        }
      }
      for(int ik=0;ik<(int)pdd.pdos_k[i].size();ik++) {
        int kp=pdd.pdos_k[i][ik];
        if(kp<1 || kp>prd.nkpts) {
          message << "Error in case: " << i+1 << endl;
          message << "too low or high a kpt number has been entered in entry: " << ik+1 << endl;
          message << "Bad value is: " << kp << endl;
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_); //CO20200624
        }
      }
      for(int ib=0;ib<(int)pdd.pdos_b[i].size();ib++) {
        int bp=pdd.pdos_b[i][ib];
        if(bp<1 || bp>prd.nbands) {
          message << "Error in case: " << i+1 << endl;
          message << "too low or high a bnd number has been entered in entry: " << ib+1 << endl;
          message << "Bad value is: " << bp << endl;
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_); //CO20200624
        }
      }
      for(int ilm=0;ilm<(int)pdd.pdos_lm[i].size();ilm++) {
        int lmp=pdd.pdos_lm[i][ilm];
        if(lmp<1 || lmp>prd.nlmtot) {
          message << "Error in case: " << i+1 << endl;
          message << "too low or high an lm number has been entered in entry: " << ilm+1 << endl;
          message << "Bad value is: " << lmp << endl;
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_); //CO20200624
        }
      }
    }
  }// end routine
}

// ***************************************************************************
// CalcPDOS
// ***************************************************************************
// Calculates the total PDOS based on the settings in pdd.
namespace pflow {
  void CalcPDOS(const projdata& prd, pdosdata& pdd) {
    double de=(pdd.emax-pdd.emin)/pdd.nbins;
    double e;
    int ibin;
    int sp_size=1;
    if(prd.sp) sp_size=3;
    pdd.pdos = aurostd::matrix<double> (pdd.nbins,2*sp_size+1,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
    for(int ib=0;ib<pdd.nbins;ib++) {
      pdd.pdos[ib][0]=pdd.emin+(double)ib*de;
    }
    for(int i=0;i<(int)pdd.pdos_at.size();i++) { // cases
      for(int ia=0;ia<(int)pdd.pdos_at[i].size();ia++) {
        int at=pdd.pdos_at[i][ia]-1;
        for(int ik=0;ik<(int)pdd.pdos_k[i].size();ik++) {
          int kp=pdd.pdos_k[i][ik]-1;
          for(int ib=0;ib<(int)pdd.pdos_b[i].size();ib++) {
            int bp=pdd.pdos_b[i][ib]-1;
            for(int ilm=0;ilm<(int)pdd.pdos_lm[i].size();ilm++) {
              int lmp=pdd.pdos_lm[i][ilm]-1;
              e = prd.ener_k_b_u[kp][bp];
              ibin=int((e-pdd.emin)/de);
              if(ibin>=0 && ibin<pdd.nbins) {
                pdd.pdos[ibin][1]+=prd.occ_vs_ion_kpt_bnd_lmtot_u[at][kp][bp][lmp]*prd.wkpt[kp];
                if(prd.sp) {
                  e = prd.ener_k_b_d[kp][bp];
                  ibin=int((e-pdd.emin)/de);
                  pdd.pdos[ibin][2]+=prd.occ_vs_ion_kpt_bnd_lmtot_d[at][kp][bp][lmp]*prd.wkpt[kp];
                }
              }// ibin in range
            }//ilm
          }//ik
        }//ib
      }//iat
    }//cases
    // Now do Up-Dn if spin polarized.
    if(prd.sp) {
      for(int ibin=0;ibin<pdd.nbins;ibin++) {
        pdd.pdos[ibin][3]=pdd.pdos[ibin][1]-pdd.pdos[ibin][2];
      }
    }

    // Normalize DOS to a density
    for(int ib=0;ib<pdd.nbins;ib++) {
      for(int i=1;i<sp_size+1;i++) {
        pdd.pdos[ib][i]/=de;
      }
    }

    // Smooth DOS
    SmoothPDOS(prd,pdd);

    // Get integrated DOS
    for(int i=1;i<sp_size+1;i++) {
      pdd.pdos[0][i+sp_size]=pdd.pdos[0][i]*de;
    }
    for(int ib=1;ib<pdd.nbins;ib++) {
      for(int i=1;i<sp_size+1;i++) {
        pdd.pdos[ib][i+sp_size]=pdd.pdos[ib-1][i+sp_size]+pdd.pdos[ib][i]*de;
      }
    }
  }//end routine
}

// ***************************************************************************
// SmoothPDOS
// ***************************************************************************
// Smooths the total PDOS based on Gaussian smoothing.
namespace pflow {
  void SmoothPDOS(const projdata& prd, pdosdata& pdd) {
    double de=(pdd.emax-pdd.emin)/pdd.nbins;
    //   double e; //DM not used
    //  int ibin; //DM not used
    int sp_size=1;
    if(prd.sp) sp_size=3;

    // Get smoothed pdos.
    for(int is=1;is<=sp_size;is++) {
      vector<double> tdos(pdd.nbins);
      for(int ib=0;ib<pdd.nbins;ib++) {
        tdos[ib]=pdd.pdos[ib][is];
      }
      double sig=pdd.smooth_sigma/de;
      tdos=SmoothFunc(tdos,sig);
      for(int ib=0;ib<pdd.nbins;ib++) {
        pdd.pdos[ib][is]=tdos[ib];
      }
    }

    //Old smooth DOS calc - does not seem to work right!
    //int range=(int)(5.0*pdd.smooth_sigma/de); // Uses gaussian out to 5 sigma.
    //matrix<double> wt(pdd.nbins,2*range+1,0.0);

    //// Get Gaussian weights
    //vector<double> norm(pdd.nbins,0.0);
    //for(int ib=0;ib<pdd.nbins;ib++) {
    //for(int i=-range;i<=range;i++) {
    //double x=(double)(i)*de+de/2.0;
    //if((ib+i)>=0 && (ib+i)<pdd.nbins) {
    //wt[ib][i+range]=Normal(x,de/2,pdd.smooth_sigma);
    //norm[ib]=norm[ib]+wt[ib][i+range];
    //}
    //}
    //}
    //// Normalize to one
    //for(int ib=0;ib<pdd.nbins;ib++) {
    //for(int i=-range;i<=range;i++) {
    //wt[ib][i+range]=wt[ib][i+range]/norm[ib];
    //}
    //}

    //// Average in weighted nearby bins.
    //for(int is=1;is<=sp_size;is++) {
    //for(int ib=0;ib<pdd.nbins;ib++) {
    //double tdos=0;
    //for(int i=-range;i<=range;i++) {
    //if((ib+i)>0 && (ib+i)<pdd.nbins) {
    //tdos=tdos+wt[ib][i+range]*pdd.pdos[ib+i][is];
    //}
    //}
    //pdd.pdos[ib][is]=tdos;
    //}
    //}
  }
}

// ***************************************************************************
// AtomCntError
// ***************************************************************************
// Smooths the total PDOS based on Gaussian smoothing.
namespace pflow {
  void AtomCntError(const string& tok, const int tokcnt, const int atom_cnt) {
    stringstream message;
    message << "in AtomCntError " <<endl;
    message << "The token " << tok << " has been defined (possibly by default) this many times: " << tokcnt << endl;
    message << "You have used the token ATOMS this many times: " << atom_cnt << endl;
    message << "These must be equal and >0 - the token ATOMS must precede each case and must occurr at least once. "<< endl;
    message << "This error is probably due to your having forgot to use the ATOMS token. Please type ATOMS=*** as the first line for each case."<< endl;
    throw aurostd::xerror(__AFLOW_FILE__,"pflow::AtomCntError():",message,_INPUT_ILLEGAL_); //CO20200624
  }
}

// ***************************************************************************
// SUMPDOSFUNCS SUMPDOSFUNCS SUMPDOSFUNCS SUMPDOSFUNCS SUMPDOSFUNCS SUMPDOSFUN
// ***************************************************************************

// ***************************************************************************
// ReadSumDOSParams
// ***************************************************************************
// Reads in parameters from PDOSParams_infile
// that decide what summations to perform over pdos.
// This routine makes use of the pdosdata object which
// contains some parameters that are not used (parameters
// for kpoints, bands, etc.)

namespace pflow {
  void ReadSumDOSParams(ifstream& infile, pflow::pdosdata& pdd) {
    // Defaults (emin,emax,nbins,smooth_sigma,print_params,efermi)
    stringstream message;
    pdd.emin=0.0;
    pdd.emax=0.0;
    pdd.nbins=301;
    pdd.smooth_sigma=(pdd.emax-pdd.emin)/pdd.nbins;
    pdd.print_params=0;
    pdd.efermi=-999;
    pdd.spin=1; // 1: non spin polarized, 2 spin polarized.          
    pdd.nlm=9; // all the possible l and m:  9 for s,p,d,  16 for s,p,d,f

    // Read in all the tokens.
    string s,s_ns;
    vector<string> token_vec;
    vector<string> val_vec;
    while (!infile.eof()) {
      getline(infile,s);
      s_ns=aurostd::RemoveSpaces(s); // Get string with no spaces.
      if(s_ns.size()>0) { // Make sure line was not blank
        if(s_ns[0]!='#') { // Exclude comment lines
          string token;
          string sval;
          int id=s.find('=');
          if(id>=(int)s.length()) {
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"the following token is incorrectly formatted: "+s,_INPUT_ILLEGAL_); //CO20200624
          }
          token=s.substr(0,id);
          token=aurostd::RemoveSpaces(token);
          int i_f=std::min((int)s.length(),(int)s.find('#')); // End of string
          sval=s.substr(id+1,i_f-id-1);
          token_vec.push_back(token);
          val_vec.push_back(sval);
        } //if s_ns[0]!=#
      } // if s_ns.size>0
    } // while !infile.eof

    // Read in values for variables in the pdosdata object.
    vector<int> tmp;
    int atom_cnt=0;
    for(int i=0;i<(int)token_vec.size();i++) {
      int id=0;
      string tok=token_vec[i];
      string sval;
      // double dval; //DM not used
      int ival;
      int found_token=0;
      if(tok=="ATOMS") {
        id=0;
        atom_cnt++;
        tmp.clear();
        // set up atoms matrix
        while (id<(int)val_vec[i].size()) {
          string s;
          s=aurostd::GetNextVal(val_vec[i],id);
          if(s.size()>0) {
            ival=atoi(s.c_str());
            tmp.push_back(ival);
          }
        }
        pdd.pdos_at.push_back(tmp);
        // add vector to pdos_lm
        tmp.clear();
        pdd.pdos_lm.push_back(tmp);
        found_token=1;
      }// ATOMS
      if(tok=="LMVALUES") {
        if((int)pdd.pdos_lm.size()!=atom_cnt || atom_cnt==0) AtomCntError(tok,(int)pdd.pdos_lm.size(),atom_cnt);
        id=0;
        tmp.clear();
        // set up atoms matrix
        while (id<(int)val_vec[i].size()) {
          string s;
          s=aurostd::GetNextVal(val_vec[i],id);
          if(s.size()>0) {
            ival=atoi(s.c_str());
            pdd.pdos_lm[pdd.pdos_lm.size()-1].push_back(ival);
          }
        }
        found_token=1;
      }// LMVALUES
      if(tok=="SMOOTH_SIGMA") {
        id=0;
        pdd.smooth_sigma=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
        found_token=1;
      }// SMOOTH_SIGMA
      if(tok=="PRINT_PARAMS") {
        id=0;
        pdd.print_params=atoi(aurostd::GetNextVal(val_vec[i],id).c_str());
        found_token=1;
      }// PRINT_PARAMS
      if(tok=="EFERMI") {
        id=0;
        pdd.efermi=atoi(aurostd::GetNextVal(val_vec[i],id).c_str());
        found_token=1;
      }// EFERMI
      if(tok=="SPIN") {
        id=0;
        pdd.spin=atoi(aurostd::GetNextVal(val_vec[i],id).c_str());
        found_token=1;
      }// SPIN
      if(tok=="NLM") {
        id=0;
        pdd.nlm=atoi(aurostd::GetNextVal(val_vec[i],id).c_str());
        found_token=1;
      }// NLM
      if(!found_token) {
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"token is not recognized: "+tok,_INPUT_ILLEGAL_); //CO20200624
      }
    } // for i over tokens

    // Check that all data is within bounds
    for(int i=0;i<(int)pdd.pdos_at.size();i++) { // cases
      for(int ilm=0;ilm<(int)pdd.pdos_lm[i].size();ilm++) {
        int lmp=pdd.pdos_lm[i][ilm];
        if(lmp<1 || lmp>pdd.nlm) {
          message << "Error in case: " << i+1 << endl;
          message << "You have entered too low or high an lm number in entry: " << ilm+1 << endl;
          message << "Bad value is: " << lmp << endl;
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_); //CO20200624
        }
      }
    }
  }// end routine
}

// ***************************************************************************
// ReadSumDOSParams
// ***************************************************************************
// Reads in pdos data from PDOSinfile.
namespace pflow {
  void ReadInPDOSData(aurostd::matrix<aurostd::matrix<double> >& allpdos, pflow::pdosdata& pdd,  //CO20200404 pflow::matrix()->aurostd::matrix()
      ifstream& infile) {
    string sdum;
    // Get natoms
    infile >> pdd.natoms;
    getline(infile,sdum);
    // Skip lines
    for(int i=0;i<4;i++) {
      getline(infile,sdum);
    }
    // Get number of bins and E Fermi
    double tmpe;
    infile >> sdum >> sdum >> pdd.nbins >> tmpe;
    getline(infile,sdum);
    if(abs(pdd.efermi+999)<1e-6) { // Reset efermi if it has special -999 value
      pdd.efermi=tmpe;
    }
    // Skip total DOS
    for(int i=0;i<pdd.nbins;i++) {
      getline(infile,sdum);
    }
    // Initialize size of allpdos[atoms][spin][lm][bin]
    aurostd::matrix<double> tmp(pdd.nlm+1,pdd.nbins);  //CO20200404 pflow::matrix()->aurostd::matrix()
    allpdos = aurostd::matrix<aurostd::matrix<double> > (pdd.natoms,pdd.spin,tmp); //CO20200404 pflow::matrix()->aurostd::matrix()
    // tpx
    //cout << spin << " " << pdd.nlm << " " << pdd.natoms << " " << pdd.nbins << endl;
    // Loop over each atom
    for(int ia=0;ia<pdd.natoms;ia++) {
      // Throw away first line since it is not DOS info
      getline(infile,sdum);
      // This is needed since getline just gets the end of the previous line when
      // the previous line has been read in by >>.  Only for atom 0 is the
      // previous line read in by getline, so only one more getline is needed.
      if(ia>0) {
        getline(infile,sdum);
      }
      // Loop over each bin
      for(int ib=0;ib<pdd.nbins;ib++) {
        // Read in energy into first spin and nlm array elements.
        infile >> allpdos[ia][0][0][ib];
        // Loop over each lm
        for(int ilm=1;ilm<=pdd.nlm;ilm++) {
          // Loop over each spin
          for(int is=0;is<pdd.spin;is++) {
            infile >> allpdos[ia][is][ilm][ib];
            // tpx
            /*
               cout << "ia ib is ilm allpdos " << ia << " "
               << ib << " "
               << is << " "
               << ilm << " "
               << allpdos[ia][is][ilm][ib] << endl;
               */
          } // spin
        } // lm
      } // bin
    } // at
  } // end routine
}

// ***************************************************************************
// SumPDOS
// ***************************************************************************
// Sums the PDOS accoring to the parameters specified.
namespace pflow {
  void SumPDOS(const aurostd::matrix<aurostd::matrix<double> >& allpdos, pflow::pdosdata& pdd) { //CO20200404 pflow::matrix()->aurostd::matrix()
    stringstream message;
    if(pdd.spin==1) {
      pdd.pdos = aurostd::matrix<double> (pdd.nbins,3);  //CO20200404 pflow::matrix()->aurostd::matrix()
    }
    else {
      pdd.pdos = aurostd::matrix<double> (pdd.nbins,7);  //CO20200404 pflow::matrix()->aurostd::matrix()
    }
    // Set energies
    for(int ib=0;ib<pdd.nbins;ib++) {
      pdd.pdos[ib][0]=allpdos[0][0][0][ib]-pdd.efermi;
    }
    // Loop over flagged atoms and lm and sum
    for(int ic=0;ic<(int)pdd.pdos_at.size();ic++) {
      for(int ida=0;ida<(int)pdd.pdos_at[ic].size();ida++) {
        int ia=pdd.pdos_at[ic][ida]-1;
        // Error check
        if(ia<0 || ia>=pdd.natoms) {
          message << "Error in case: " << ic+1 << endl;
          message << "You have entered too low or high an atom number in entry: " << ida+1 << endl;
          message << "Min/Max allowed values are: " << "1 / " << pdd.natoms << endl;
          message << "Bad value is: " << ia+1 << endl;
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_); //CO20200624
        }
        // Loop over each lm
        for(int idlm=0;idlm<(int)pdd.pdos_lm[ic].size();idlm++) {
          int ilm=pdd.pdos_lm[ic][idlm];
          // Error check
          if(ilm<1 || ilm>pdd.nlm) {
            message << "Error in case: " << ic+1 << endl;
            message << "You have entered too low or high an lm number in entry: " << idlm+1 << endl;
            message << "Min/Max allowed values are: " << "1 / " << pdd.nlm << endl;
            message << "Bad value is: " << ilm << endl;
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_); //CO20200624
          }
          // Loop over each bin
          for(int ib=0;ib<pdd.nbins;ib++) {
            // Do spin sum
            //  double val; //DM not used
            // Note: allpdos[atoms][spin][lm][bin]
            switch (pdd.spin) {
              case 1:{ // non spin polarized
                       pdd.pdos[ib][1]=pdd.pdos[ib][1]+allpdos[ia][0][ilm][ib];
                       break;
                     }
              case 2:{ // spin polarized
                       pdd.pdos[ib][1]=pdd.pdos[ib][1]+allpdos[ia][0][ilm][ib];
                       pdd.pdos[ib][2]=pdd.pdos[ib][2]+allpdos[ia][1][ilm][ib];
                       break;
                     }
              default:
                     message << "ERROR: Error in case: " << ic+1 << endl;
                     message << "ERROR: You have entered too low or high an spin number" << endl;
                     message << "ERROR: Min/Max allowed values are: " << "1 / 2" << endl;
                     message << "ERROR: Bad value is: " << pdd.spin << endl;
                     throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_); //CO20200624
            } // switch spin
          } // bin
        } // lm
      } // at
    } // case
    // Get cumulative and difference values
    double dE=pdd.pdos[1][0]-pdd.pdos[0][0];
    if(pdd.spin==1) {
      pdd.pdos[0][2]=0;
      for(int ib=1;ib<pdd.nbins;ib++) {
        pdd.pdos[ib][2]=pdd.pdos[ib-1][2]+pdd.pdos[ib][1]*dE;
      }
    }
    else {
      pdd.pdos[0][3]=pdd.pdos[0][1]-pdd.pdos[0][2];
      pdd.pdos[0][4]=0;
      pdd.pdos[0][5]=0;
      pdd.pdos[0][6]=0;
      for(int ib=1;ib<pdd.nbins;ib++) {
        pdd.pdos[ib][3]=pdd.pdos[ib][1]-pdd.pdos[ib][2];
        pdd.pdos[ib][4]=pdd.pdos[ib-1][4]+pdd.pdos[ib][1]*dE;
        pdd.pdos[ib][5]=pdd.pdos[ib-1][5]+pdd.pdos[ib][2]*dE;
        pdd.pdos[ib][6]=pdd.pdos[ib-1][6]+pdd.pdos[ib][3]*dE;
      }
    }
  }
}

// ***************************************************************************
// RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBF
// ***************************************************************************

// ***************************************************************************
// Function TotalAtomDist
// ***************************************************************************
// NOTE: THIS FUNCTION IS OUTDATED BY RBPOSCAR DISP FUNC BELOW.
// This function gets the total distances between
// two structures by taking the sqrt of the sum of
// the squared distances between each atom.  
// If path_flag = N/n then the nearest
// images are used to measure distance.
namespace pflow {
  double TotalAtomDist(xstructure str, xstructure str00, const string& path_flag) {
    str=ReScale(str,1);
    str00=ReScale(str00,1);
    aurostd::matrix<double> cpos=pflow::GetCpos(str);  //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> cpos00=pflow::GetCpos(str00);  //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> lat=pflow::GetLat(str);  //CO20200404 pflow::matrix()->aurostd::matrix()
    int nat=std::min(cpos.size(),cpos00.size());
    double dtot=0;
    for(int iat=0;iat<nat;iat++) {
      vector<double> dp=pflow::VVdiff(cpos[iat],cpos00[iat]);
      // If path_flag is n/N then the path is taken between
      // nearest images.  Otherwise path is between the atoms given.
      if(path_flag=="n" || path_flag=="N") {
        vector<double> ddp=pflow::vecC2F(lat,dp);
        for(int ic=0;ic<3;ic++) {
          ddp[ic]=ddp[ic]-Nint(ddp[ic]);
        }
        dp=pflow::vecF2C(lat,ddp);
      }
      dtot=dtot+pflow::norm(dp)*pflow::norm(dp);
    }
    return sqrt(dtot);
  }
}

// ***************************************************************************
// Function GetRBDir
// ***************************************************************************
// This function gets the directories of the rubber
// band run.  The directories are assumed to
// go from 00 to 0(nim+1) (or 00 to (nim+1) if nim>8).
namespace pflow {
  vector<string> GetRBDir(const int& nim) {
    vector<string> rbdir(nim+2);
    for(int im=0;im<nim+2;im++) {
      ostringstream tmp;
      if(im<10) {
        tmp << "0" << im << ends;
      }
      else {
        tmp << im << ends;
      }
      rbdir[im]=tmp.str();
    }
    return rbdir;
  }
}

// ***************************************************************************
// Function GetRBEner
// ***************************************************************************
// This function gets the energy from the directories
// of a rubber band run. The energies are looked for
// in the OSZICAR file in the line before the last line.
namespace pflow {
  vector<double> GetRBEner(const int& nim) {
    vector<double> ener((nim+2),0.0);
    std::string tmpfile=aurostd::TmpFileCreate();
    vector<string> rbdir=GetRBDir(nim);
    for(int im=0;im<nim+2;im++) {
      ostringstream cmd;
      cmd << "cd " << rbdir[im] << ";";
      cmd << "tail -2 OSZICAR | head -1 | awk '{if($1==\"CG\") {print $4;}} {if($1==\"DAV:\") {print $3;}} {if($1==\"RMM:\") {print $3;}}' >> " << tmpfile << ";" 	<< "cd ..;" << ends;
      system(cmd.str().c_str());
      ifstream ifaustream(tmpfile.c_str());
      ifaustream >> ener[im];
      aurostd::RemoveFile(tmpfile);
    }
    return ener;
  }
}

// ***************************************************************************
// Function GetRBStruct
// ***************************************************************************
// This function gets the structures from a rubber
// band run.  The directories are assumed to go from 00 to 0(nim+1)
// (or 00 to (nim+1) if nim>8).  
namespace pflow {
  vector<xstructure> GetRBStruct(const int& nim) {
    vector<xstructure> vstr(nim+2);
    vector<string> rbdir=GetRBDir(nim);
    for(int im=0;im<nim+2;im++) {
      ostringstream file;
      if(im==0 || im==(nim+1)) { // First and last files are POSCAR
        file << rbdir[im] << "/POSCAR" << ends;
      }
      else {
        file << rbdir[im] << "/CONTCAR" << ends;
      }
      ifstream ifaustream(file.str().c_str());
      ifaustream >> vstr[im];
      ifaustream.close();
      ifaustream.clear();
    }
    return vstr;
  }
}

// ***************************************************************************
// Function GetRBDistCum
// ***************************************************************************
// This function gets the distances between images
// from the directories of a rubber band run. The initial distance
// for image i is taken by calculating the sum of the Euclidean
// distances between all atoms from image i to image i-1.
// The distances are then summed, so that d[i] is the cumulative
// distance to structure i.
namespace pflow {
  vector<double> GetRBDistCum(const vector<xstructure>& vstr, const string& path_flag) {
    int nim = vstr.size()-2;
    vector<double> d((nim+2),0.0);
    xstructure str=vstr[0];
    xstructure last_str=vstr[0];
    xstructure diffstr;
    aurostd::matrix<double> cm;  //CO20200404 pflow::matrix()->aurostd::matrix()

    for(int im=0;im<nim+2;im++) {
      str=vstr[im];
      double dist=0.0;
      RBPoscarDisp(last_str,str,diffstr,dist,cm,path_flag);
      d[im]=dist;
      last_str=str;
    }

    // Turn d into an cumulative distance.
    for(int i=1;i<(int) d.size();i++) {
      d[i]=d[i]+d[i-1];
    }
    return d;
  }
}

// ***************************************************************************
// Function GetRBDistFromStrI
// ***************************************************************************
// This function gets the distances between images
// from the directories of a rubber band run and a given
// structure, I (dist[j] = strI - str[j]).
namespace pflow {
  vector<double> GetRBDistFromStrI(const vector<xstructure>& vstr,
      const xstructure& strI,
      const string& path_flag) {
    int nim = vstr.size()-2;
    vector<double> d((nim+2),0.0);
    xstructure diffstr;
    aurostd::matrix<double> cm;  //CO20200404 pflow::matrix()->aurostd::matrix()
    for(int im=0;im<nim+2;im++) {
      RBPoscarDisp(vstr[im],strI,diffstr,d[im],cm,path_flag);
    }
    return d;
  }
}

// ***************************************************************************
//  Function RBPoscarDisp
// ***************************************************************************
// This function gets the displacement between two POSCAR files (str2-str1).  
namespace pflow {
  void RBPoscarDisp(const xstructure& str1in, const xstructure& str2in,
      xstructure& diffstr, double& totdist, aurostd::matrix<double>& cm, //CO20200404 pflow::matrix()->aurostd::matrix()
      const string& path_flag) {
    diffstr=str1in;
    xstructure str1=str1in;
    xstructure str2=str2in;
    str1=ReScale(str1,1);
    str2=ReScale(str2,1);
    aurostd::matrix<double> cpos1=pflow::GetCpos(str1);  //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> cpos2=pflow::GetCpos(str2);  //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> lat1=pflow::GetLat(str1);  //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> lat2=pflow::GetLat(str1);  //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> latdiff(3,3);  //CO20200404 pflow::matrix()->aurostd::matrix()
    cm = aurostd::matrix<double> (2);  //CO20200404 pflow::matrix()->aurostd::matrix()
    int nat=min(cpos1.size(),cpos2.size());
    aurostd::matrix<double> cposdiff(nat,3); //CO20200404 pflow::matrix()->aurostd::matrix()
    totdist=0;
    for(int ic=0;ic<3;ic++) {
      latdiff[ic]=pflow::VVdiff(lat2[ic],lat1[ic]);
    }
    for(int iat=0;iat<nat;iat++) {
      vector<double> dp=pflow::VVdiff(cpos2[iat],cpos1[iat]);
      // If path_flag is n/N then the path is taken between
      // nearest images.  Otherwise path is between the atoms given.
      if(path_flag=="n" || path_flag=="N") {
        vector<double> ddp=pflow::vecC2F(lat1,dp);
        for(int ic=0;ic<3;ic++) {
          ddp[ic]=ddp[ic]-Nint(ddp[ic]);
        }
        dp=pflow::vecF2C(lat1,ddp);
      }
      totdist=totdist+pflow::norm(dp)*pflow::norm(dp);
      cposdiff[iat]=dp;
      // tpx I think
      //    Vout(cposdiff[iat],cout);
      //    cout << "  " << iat << endl;
    }
    totdist=sqrt(totdist);
    diffstr=pflow::SetAllAtomPos(diffstr,cposdiff,1);
    //  diffstr.SetLat(latdiff); // Just keep lat=lat1.
    cm[0]=xvector2vector(GetMom1(str1));
    cm[1]=xvector2vector(GetMom1(str2));
  }
}

// ***************************************************************************
// CHARGE FUNCS CHARGE FUNCS CHARGE FUNCS CHARGE FUNCS CHARGE FUNCS CHARGE FUN
// ***************************************************************************

// **************************************************
// Function CompAperpB
// **************************************************
// This function returns the component of A perpendicular to B.

namespace pflow {
  vector<double> CompAperpB(const vector<double>& a, const vector<double>& b) {
    double dp;
    dp=pflow::VVprod(a,b);
    double nb=pflow::norm(b);
    return pflow::VVdiff(a,pflow::SVprod(dp/(nb*nb),b));
  }
}

namespace pflow {
  vector<double> GetDispToAtomDir(const xstructure& a, const vector<double>& in_fpos, const int at_num) {
    // returns in_fpos-(fpos of atom at_num) in direct coordinates.  in_fpos should be in direct coords.
    vector<double> disp(3,0.0);
    for(int ic=0;ic<3;ic++) {
      //   disp[ic]=in_fpos[ic]-fpos[at_num][ic]; // CONVASP
      disp[ic]=in_fpos[ic]-a.atoms.at(at_num).fpos(ic+1); // AFLOW STRUCTURE
      //  if(disp[ic]>0.5) disp[ic]=disp[ic]-1;
      //  if(disp[ic]<-0.5) disp[ic]=disp[ic]+1;
    }
    return disp;
  }
  double GetDistToAtomDir(const xstructure& a, const vector<double>& in_fpos, const int at_num) {
    // returns norm(in_fpos-(pos of atom at_num)).  in_fpos should be in direct coords.
    vector<double> disp(3);
    disp=GetDispToAtomDir(a,in_fpos,at_num);
    disp=pflow::vecF2C(xmatrix2matrix(a.lattice),disp);
    return a.scale*pflow::norm(disp);
  }
  vector<double> GetDispToAtomImageDir(const xstructure& a,const vector<double>& in_fpos, const int at_num) {
    // returns {in_fpos}-{fpos of atom at_num} in direct coordinates, where the difference is taken
    // between the nearest images (braces denote the fact that each set point is really a set of points
    // in different image cells.  in_fpos should be in direct coords.
    vector<double> disp(3,0.0);
    for(int ic=0;ic<3;ic++) {
      //      disp[ic]=in_fpos[ic]-fpos[at_num][ic]; // CONVASP
      disp[ic]=in_fpos[ic]-a.atoms.at(at_num).fpos(ic+1); // AFLOW
      disp[ic]=disp[ic]-(int)disp[ic];
      if(disp[ic]>0.5) disp[ic]=disp[ic]-1;
      if(disp[ic]<-0.5) disp[ic]=disp[ic]+1;
    }
    return disp;
  }
  double GetDistToAtomImageDir(const xstructure& a,const vector<double>& in_fpos, const int at_num) {
    // returns {in_fpos}-{fpos of atom at_num} in direct coordinates, where the difference is taken
    // between the nearest images (braces denote the fact that each set point is really a set of points
    // in different image cells.  in_fpos should be in direct coords.
    vector<double> disp(3);
    disp=pflow::GetDispToAtomImageDir(a,in_fpos,at_num);
    disp=pflow::vecF2C(xmatrix2matrix(a.lattice),disp);
    return a.scale*norm(disp);
  }
}

namespace pflow {
  void GetChgInt(vector<aurostd::matrix<double> >& rad_chg_int, aurostd::matrix<double>& vor_chg_int) {  //CO20200404 pflow::matrix()->aurostd::matrix()
    // Read in CHGCAR
    // Read in intial POSCAR format at beginning of file
    xstructure str;
    cin >> str;
    // Read in the grid for total charge
    vector<int> ngrid(3);
    int npts;
    cin >> ngrid[0] >> ngrid[1] >> ngrid[2];
    npts=ngrid[0]*ngrid[1]*ngrid[2];
    // Read in the total charge
    vector<double> chg_tot(npts,0.0);
    for(int i=0;i<npts;i++) {
      cin >> chg_tot[i];
    }
    // Read in dummy fields
    // If first field is augmentation then read 4*natoms+natoms-1 more fields. For vasp 4.4.5.
    // If first field is not augmentation then read natoms-1 more fields. For vasp 4.4.1.
    string ds;
    int natoms=pflow::GetNumAtoms(str);
    cin >> ds;
    string keyword = "augmentation";
    if(ds==keyword) {
      for(int i=0;i<(5*natoms-1);i++) {
        cin >> ds;
      }
    }
    else {
      for(int i=0;i<(natoms-1);i++) {
        cin >> ds;
      }
    }

    // Read in the grid for diff charge
    cin >> ngrid[0] >> ngrid[1] >> ngrid[2];
    npts=ngrid[0]*ngrid[1]*ngrid[2];
    // Read in the diff charge
    vector<double> chg_diff(npts,0.0);
    for(int i=0;i<npts;i++) {
      cin >> chg_diff[i];
    }
    // Read in 4*natoms=4*natoms dummy fields.  This is
    // not necessary and the lines only exist in vasp 4.4.5.
    /*
       for(int i=0;i<4*natoms;i++) {
       cin >> ds;
       }
       */

    // Loop over grid and bin charge densities.

    // Radial charge integration.
    // Integral over a sphere will be done by setting up NRBIN bins, each
    // of witdh DR, our to RMAX.  Then each charge at each point will be put is the
    // appropriate bin for each atom.  Note that the points and all possible relevant
    // images must be looped over.
    double DR=0.05; // width of radial bins.
    double RMAX=3; // Max distance to bin to.
    int NRBIN=(int)(RMAX/DR); // First bin is 0->DR, last bin is RMAX-DR->RMAX
    aurostd::matrix<double> dum_mat(NRBIN,5,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
    rad_chg_int = vector<aurostd::matrix<double> > (natoms,dum_mat); //CO20200404 pflow::matrix()->aurostd::matrix()
    vector<int> ig(3);
    // Determine the number of images along each lattice vector which must be considered.
    // You will center the charge density on each atom so the max number of images along a
    // given lattice param should be no more than however many lattice params are needed
    // to make sure that you are more than RMAX from an atom along that lattice param direction.
    vector<int> imax(3);
    double scale=pflow::GetScale(str);
    aurostd::matrix<double> lat=pflow::GetLat(str);  //CO20200404 pflow::matrix()->aurostd::matrix()
    for(int ic=0;ic<3;ic++) {
      imax[ic]=(int)(RMAX/(scale*pflow::norm(lat[ic]))+1)*ngrid[ic];
    }
    //    double sumchg=0;  //DM not used
    // Atom loop
    aurostd::matrix<double> at_fpos=pflow::GetFpos(str); //CO20200404 pflow::matrix()->aurostd::matrix()
    for(int iat=0;iat<natoms;iat++) {
      // Initialize radial positions in rad_chg_int    
      for(int ib=0;ib<NRBIN;ib++) {
        rad_chg_int[iat][ib][0]=(ib+1)*DR;
      }
      // Find location of atom on grid.
      vector<int> at_grid_pos(3);
      for(int ic=0;ic<3;ic++) {
        at_grid_pos[ic]=Nint(at_fpos[iat][ic]*ngrid[ic]);
      }
      // Grid loop: This loops over grid points shifted to center around grid point approximation to
      // position of iat.  The ig values therefore need to be modified to be true grid points
      // (true_ig) by undoing the shift (adding back the at_grid_pos).
      for(ig[0]=-imax[0];ig[0]<imax[0];ig[0]++) {
        for(ig[1]=-imax[1];ig[1]<imax[1];ig[1]++) {
          for(ig[2]=-imax[2];ig[2]<imax[2];ig[2]++) {
            // Get true grid points
            vector<int> true_ig(3);
            vector<double> chg_fpos(3);
            for(int ic=0;ic<3;ic++) {
              true_ig[ic]=ig[ic]+at_grid_pos[ic];
              // Get true position of this grid point in direct coords.
              chg_fpos[ic]=((float)true_ig[ic]/(float)ngrid[ic]);
              // Shift true grid point back into the 000 cell.
              true_ig[ic]=true_ig[ic]-(int)((float)true_ig[ic]/(float)ngrid[ic])*ngrid[ic];
              if(true_ig[ic]<0) true_ig[ic]=true_ig[ic]+ngrid[ic];
            } // ic
            // Get charge values for each grid position
            int id=true_ig[0]+true_ig[1]*ngrid[0]+true_ig[2]*ngrid[0]*ngrid[1];
            double chgtot=chg_tot[id];
            double chgdiff=chg_diff[id];
            double chgup=0.5*(chg_tot[id]+chg_diff[id]);
            double chgdn=0.5*(chg_tot[id]-chg_diff[id]);
            vector<double> fpos(3);
            // Get distance from iat to true grid point position.
            double dist=pflow::GetDistToAtomDir(str,chg_fpos,iat);
            // Get bin for this atom and increment radial charge density.
            int ibin=(int)(dist/DR);	  
            if(ibin<NRBIN) {
              rad_chg_int[iat][ibin][1]=rad_chg_int[iat][ibin][1]+chgtot;
              rad_chg_int[iat][ibin][2]=rad_chg_int[iat][ibin][2]+chgdiff;
              rad_chg_int[iat][ibin][3]=rad_chg_int[iat][ibin][3]+chgup;
              rad_chg_int[iat][ibin][4]=rad_chg_int[iat][ibin][4]+chgdn;
            } // if ibin<NRBIN
          } // ig[2]
        } // ig[1]
      } // ig[0]
    } // iat

    // Voronoi charge integration
    vor_chg_int = aurostd::matrix<double> (natoms,4,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
    for(ig[0]=0;ig[0]<ngrid[0];ig[0]++) {
      for(ig[1]=0;ig[1]<ngrid[1];ig[1]++) {
        for(ig[2]=0;ig[2]<ngrid[2];ig[2]++) {
          // Get charge values for this grid position
          int id=ig[0]+ig[1]*ngrid[0]+ig[2]*ngrid[0]*ngrid[1];
          double chgtot=chg_tot[id];
          double chgdiff=chg_diff[id];
          double chgup=0.5*(chg_tot[id]+chg_diff[id]);
          double chgdn=0.5*(chg_tot[id]-chg_diff[id]);
          vector<double> fpos(3);
          // Set up grid positions in direct coordinates.
          for(int ic=0;ic<3;ic++) {
            fpos[ic]=((float)ig[ic])/(float)ngrid[ic];
          } // ic
          // Loop over atoms
          double dist;
          double min_dist=pflow::GetDistToAtomImageDir(str,fpos,0);
          int min_dist_id=0;
          for(int iat=0;iat<natoms;iat++) {
            // Get "image" distance to each atom (distance to closest image)
            dist=pflow::GetDistToAtomImageDir(str,fpos,iat);
            // Keep track of which atom is the minimum "image" distance away
            if(dist<min_dist) {
              min_dist=dist;
              min_dist_id=iat;
            } // if dist<min_dist
          } // iat
          vor_chg_int[min_dist_id][0]=vor_chg_int[min_dist_id][0]+chgtot;
          vor_chg_int[min_dist_id][1]=vor_chg_int[min_dist_id][1]+chgdiff;
          vor_chg_int[min_dist_id][2]=vor_chg_int[min_dist_id][2]+chgup;
          vor_chg_int[min_dist_id][3]=vor_chg_int[min_dist_id][3]+chgdn;
        } // ig[2]
      } // ig[1]
    } // ig[0]

    // Normalize and turn rad_chg_int from chg in each radial shell to cumulative integral.
    for(int iat=0;iat<natoms;iat++) {
      for(int i=0;i<4;i++) {
        vor_chg_int[iat][i]=vor_chg_int[iat][i]/npts;
      }
      for(int i=1;i<5;i++) {
        for(int ib=0;ib<NRBIN;ib++) {
          rad_chg_int[iat][ib][i]=rad_chg_int[iat][ib][i]/npts;
          if(ib>0) rad_chg_int[iat][ib][i]=(rad_chg_int[iat][ib][i]+rad_chg_int[iat][ib-1][i]);
        } // ib
      } // for loop over tot,diff,up,dn
    } // iat
  }
}

namespace pflow {
  void pd_params::Print(ostream& outf) const{
    outf << type << " # Type " << endl;
    outf << scale << " # scale " << endl;
    for(int ic=0;ic<3;ic++) {
      for(int jc=0;jc<3;jc++) {
        outf << pts[ic][jc] << " ";
      }
      outf << " # pts " << endl;
    }
    for(int ic=0;ic<3;ic++) {
      for(int jc=0;jc<3;jc++) {
        outf << dpts[ic][jc] << " ";
      }
      outf << " # dpts " << endl;
    }
    outf << orig_loc << " # origin_loc " << endl;
    outf << ortho << " # ortho " << endl;
  }
}

// ***************************************************************************
// pflow::ReadCHGCAR
// ***************************************************************************
namespace pflow {
  bool ReadCHGCAR(xstructure& str,
      stringstream& chgcar_header,
      vector<int>& ngrid,
      vector<int>& format_dim,    
      vector<double>& chg_tot,
      vector<double>& chg_diff,
      stringstream& chgcar_ss,
      ostream& oss) {
    // This funtion reads in the CHGCAR or AECCAR file from a file.
    // Operates under assumption that it was formatted by v46s (updated by CO)
    // v46s is standard for AFLOWLIB
    // Format:
    //      Starts with POSCAR
    //      newline
    //      NX NY NZ
    //      NX*NY*NZ entries with a newline after each one (grouped by 5 usually, but 
    //          this code can read any grouping) (Tot dens = pUP+pDn)
    //      PAW augmentation occupancies, one for each atom
    //      If spin-polarized:
    //          0's line
    //          NX NY NZ
    //          NX*NY*NZ entries with a newline after each one (grouped by 5) 
    //              (Diff dens = pUP-pDN)
    //          PAW augmentation occupancies, one for each atom
    // AECCAR's will stop after total density
    //
    // PAW augmentation occupancies are not used at all in AFLOW, and all other ReadCHGCAR 
    // scripts simply ignore it
    // there also seems to be some major discrepancies between how different versions of 
    // VASP format the CHGCAR
    // OBSOLETE:
    // here's text for previous versions (OBSOLETE):
    //      If first field is augmentation then read 4*natoms+natoms-1 more fields. For vasp 4.4.5.
    //      If first field is not augmentation then read natoms-1 more fields. For vasp 4.4.1.
    // previous version text (OBSOLETE):
    //      Here is the format of a CHGCAR file (for v445 with
    //      spin polarization from DEC linux).
    //      Note that for Intel linux vasp445 the newlines after
    //      each chg density line seem to missing, however, the
    //      code still reads this in correctly.
    //      newline
    //      Starts with a POSCAR file
    //      newline
    //      NX NY NZ
    //      NX*NY*NZ entries with a newline after each one (grouped by 5) (Tot dens = pUP+pDn)
    //      Natom lines, each with 4 fields, with a newline after each one
    //      Natom entries (grouped by 5).
    //      NX NY NZ
    //      NX*NY*NZ entries with a newline after each one (grouped by 5) (Diff dens = pUP-pDN)
    //      Natom lines, each with 4 fields, with a newline after each one
    //      Natom entries (grouped by 5). // I don't think these are actually in the CHGCAR.


    // DEBUG
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) oss << __AFLOW_FUNC__ << "BEGIN" << endl;

    string content,line, keyword="augmentation";
    vector<string> vcontent, sum_tokens, line_tokens;
    vector<int> ngrid_compare;
    content=chgcar_ss.str();
    aurostd::string2vectorstring(content,vcontent,true,true); //consecutive=true, trim_edges=true, CO20170613, we have issues with "\n \n" vs. "\n\n", we want both to be separated
    stringstream poscar;
    chgcar_header.str("");
    int grid_point=0;
    uint natoms=0, linecount=0, npts=1, npts_compare=1, numcolumns=0, num_skip=0;
    bool got_chgtot=false;      //so we know that we're now reading chgdiff
    if(LDEBUG) oss << __AFLOW_FUNC__ << "CHECKING FORMAT OF SYSTEM NAME AND SCALING FACTOR" << endl;
    //do some quick checks to make sure it's format correctly, otherwise kill the script
    //first line is comment line, skip check
    //next line contains scaling factor
    aurostd::string2tokens(vcontent.at(1),line_tokens," ");
    if(line_tokens.size()!=1) {
      oss << endl;
      oss << __AFLOW_FUNC__ << "ERROR: Unrecognized CHGCAR header format (scaling factor)." << endl;
      oss << endl;
      return FALSE;
    }
    line_tokens.clear();
    if(LDEBUG) oss << __AFLOW_FUNC__ << "FORMAT FOR SCALING FACTOR LOOKS OK" << endl;

    //next three lines should contain lattice vectors
    if(LDEBUG) oss << __AFLOW_FUNC__ << "CHECKING FORMAT FOR LATTICE VECTORS" << endl;
    for(uint i=2;i<5;i++) {
      aurostd::string2tokens(vcontent.at(i),line_tokens," ");
      if(line_tokens.size()!=3) {
        oss << endl;
        oss << __AFLOW_FUNC__ << "ERROR: Unrecognized CHGCAR header format (lattice vectors)." << endl;
        oss << endl;
        return FALSE;
      }
    }
    line_tokens.clear();
    if(LDEBUG) oss << __AFLOW_FUNC__ << "FORMAT FOR LATTICE VECTORS LOOKS OK" << endl;

    //DX20190618 - check for VASP5 - START
    //mirrors check in aflow_xatom.cpp
    bool is_vasp_4_poscar = true; //VASP4 (default)
    uint poscar_header_count = 7 ; // VASP4 (default)
    uint species_line_number = 5 ; // VASP4 (default)
    string stmp = vcontent.at(6);
    aurostd::StringSubst(stmp,"\t"," ");aurostd::StringSubst(stmp,"  "," ");aurostd::StringSubst(stmp,"  "," ");
    aurostd::string2tokens(stmp,line_tokens);
    if(line_tokens.at(0)[0]!='S' && line_tokens.at(0)[0]!='s' && line_tokens.at(0)[0]!='D' && line_tokens.at(0)[0]!='d' && line_tokens.at(0)[0]!='C' && line_tokens.at(0)[0]!='c') {
      is_vasp_4_poscar = false; // then VASP5
      poscar_header_count = 8; // then VASP5
      species_line_number = 6; // then VASP5
    }
    line_tokens.clear();
    if(LDEBUG) oss << __AFLOW_FUNC__ << "IS VASP4 FORMAT = " << is_vasp_4_poscar << endl;
    //DX20190618 - check for VASP5 - END

    //if we get here, assume CHGCAR is formatted correctly for now
    //read in header, scaling factor, lattice vectors, number of atoms, coordinate type
    if(LDEBUG) oss << __AFLOW_FUNC__ << "READING POSCAR" << endl;
    //DX20190618 [OBSOLETE] for(uint i=0;i<7;i++)
    for(uint i=0;i<poscar_header_count;i++) //DX20190618
    { //CO20200106 - patching for auto-indenting
      poscar << vcontent.at(i) << endl;
      chgcar_header << vcontent.at(i) << endl;
      linecount=i;
    }
    //read in number of atoms
    //DX20190618 aurostd::string2tokens(vcontent.at(5),sum_tokens," ");
    aurostd::string2tokens(vcontent.at(species_line_number),sum_tokens," "); //DX20190618
    for(uint i=0;i<sum_tokens.size();i++) {
      natoms+=aurostd::string2utype<uint>(sum_tokens.at(i));
    }
    linecount++;   //start on coordinate line
    //read in coordinates
    for(uint i=0;i<natoms;i++) {
      poscar << vcontent.at(linecount) << endl;
      chgcar_header << vcontent.at(linecount) << endl;
      linecount++;
    }
    //define xstruture
    if(LDEBUG) oss << __AFLOW_FUNC__ << "DEFINE XSTRUCTURE" << endl;
    str=poscar;
    if(LDEBUG) oss << str << endl;
    //skip empty line
    linecount++;
    //now we get chg_tot and chg_diff (if it exists)
    if(LDEBUG) oss << __AFLOW_FUNC__ << "READ CHG VALUES" << endl;
    while(vcontent.size()>linecount) {
      //see if it's augmentation occupancies that we're reading, and skip appropriate number of lines
      if(vcontent.at(linecount).compare(0,12,keyword)==0) {
        if(LDEBUG) oss << __AFLOW_FUNC__ << "GET FORMAT FOR AUGMENTATION OCCUPANCIES" << endl;
        for(uint i=0;i<natoms;i++) {
          //check formatting of line first, 'augmentation occupancies natom npts'
          aurostd::string2tokens(vcontent.at(linecount),line_tokens," ");
          if(line_tokens.size()!=4) {
            oss << endl;
            oss << __AFLOW_FUNC__ << "ERROR: Unrecognized CHGCAR format (augmentation occupancies)." << endl;
            oss << endl;
            return FALSE;
          }
          npts=aurostd::string2utype<int>(line_tokens.at(3));
          line_tokens.clear();
          aurostd::string2tokens(vcontent.at(linecount+1),line_tokens," ");
          numcolumns=line_tokens.size();
          // this npts (augmentation occupancies) will be different than chgtot, grab it for comparsion
          format_dim.push_back(npts);
          format_dim.push_back(numcolumns);
          num_skip=(uint)(std::ceil((double)npts/(double)numcolumns));    //skip right number of liens
          linecount+=num_skip+1;  //+1 to get to next augmentation line
          if(LDEBUG) oss << __AFLOW_FUNC__ << "SUCCESSFULLY GATHERED FORMAT FOR AUGMENTATION OCCUPANCIES" << endl;
        }
      } else {    //must be chg_tot or chg_diff
        npts=1; ngrid.clear();      //clear for next round
        line_tokens.clear();
        //if gathering chg_diff, we have to ignore extra line of 0's
        if(got_chgtot) {
          linecount++;
        }
        //get npts
        sum_tokens.clear();
        aurostd::string2tokens(vcontent.at(linecount),sum_tokens, " ");
        //check formatting, this line should contain npts
        if(sum_tokens.size()!=3) {
          oss << endl;
          oss << __AFLOW_FUNC__ << "ERROR: Unrecognized CHGCAR format (number of grid points)." << endl;
          oss << endl;
          return FALSE;
        }
        for(uint i=0;i<sum_tokens.size();i++) {
          grid_point=aurostd::string2utype<int>(sum_tokens.at(i));   
          npts=npts*grid_point;
          ngrid.push_back(grid_point);
        }
        //initialize once
        if(!got_chgtot) {
          if(LDEBUG) oss << __AFLOW_FUNC__ << "READING CHG_DIFF" << endl;
          //now with npts, initialize chg_tot and chg_diff
          chg_tot=vector<double>(npts,0.0);
          ngrid_compare=ngrid;
          npts_compare=npts;  //save for comparison later with chg_diff
        } else {
          if(LDEBUG) oss << __AFLOW_FUNC__ << "READING CHG_TOT" << endl;
          chg_diff=vector<double>(npts,1.0);
          if(npts!=npts_compare) {    //chg_tot and chg_diff have different number of points
            oss << endl;
            oss << __AFLOW_FUNC__ << "ERROR: Number of grid points for CHG_tot and CHG_diff do not match. " << endl;
            oss << __AFLOW_FUNC__ << "ERROR: npts_total: " << npts << endl;
            oss << __AFLOW_FUNC__ << "ERROR: npts_diff: " << npts_compare << endl;
            oss << __AFLOW_FUNC__ << "ERROR: This will give nonsense!! " << endl;
            oss << endl;
            return FALSE;
          }
          if(ngrid!=ngrid_compare) {  //chg_tot and chg_diff have different grid
            oss << endl;
            oss << __AFLOW_FUNC__ << "ERROR: Grids for CHG_tot and CHG_diff do not match. " << endl;
            oss << __AFLOW_FUNC__ << "ERROR: ngrid_total: " << ngrid.at(0) << " " << ngrid.at(1) << " " << ngrid.at(2) << endl;
            oss << __AFLOW_FUNC__ << "ERROR: ngrid_diff: " << ngrid_compare.at(0) << " " << ngrid_compare.at(1) << " " << ngrid_compare.at(2) << endl;
            oss << __AFLOW_FUNC__ << "ERROR: This will give nonsense!! " << endl;
            oss << endl;
            return FALSE;
          }
        }
        //assume VASP CHGCAR format (number of columns) is constant, grab number
        linecount++;
        aurostd::string2tokens(vcontent.at(linecount),line_tokens," ");
        numcolumns=line_tokens.size();
        //already grab npts for comparsion, we just need numcolumns
        format_dim.push_back(numcolumns);

        for(uint i=0;i<(uint) npts;i+=numcolumns) {
          line_tokens.clear();
          aurostd::string2tokens(vcontent.at(linecount),line_tokens," ");
          numcolumns=line_tokens.size();  //there may be lines at end with only 1 or 2 entries
          for(uint j=0;j<numcolumns;j++) {
            if(!got_chgtot) {
              chg_tot.at(i+j)=aurostd::string2utype<double>(line_tokens.at(j));
            } else {
              chg_diff.at(i+j)=aurostd::string2utype<double>(line_tokens.at(j));
            }
          }
          linecount++;
        }
        if(!got_chgtot) {
          got_chgtot=true;
        }
        line_tokens.clear();
        if(LDEBUG) oss << __AFLOW_FUNC__ << "SUCCESSFULLY READ CHG VALUES" << endl;
      }
    }
    if(LDEBUG) oss << __AFLOW_FUNC__ << "DONE" << endl;
    return TRUE;
  }
} // namespace pflow

// ***************************************************************************
// pflow::ReadCHGCAR
// ***************************************************************************
namespace pflow {
  bool ReadChg(xstructure& str,
      vector<int>& ngrid,
      vector<double>& chg_tot,
      vector<double>& chg_diff,
      istream& chgfile) {
    //OLD READCHG FORMAT, just keeping it to not screw up other functions that depend on it
    stringstream chgcar_header,chgcar_ss;
    vector<int> format_dim;
    chgcar_ss << chgfile.rdbuf();
    ostringstream oss;
    return ReadCHGCAR(str,chgcar_header,ngrid,format_dim,chg_tot,chg_diff,chgcar_ss,oss);
  }
} // namespace pflow

// ***************************************************************************
// GetChgInt
// ***************************************************************************
//  This funtion reads in the CHGCAR file from standard input
//  and returns a matrix of charge integrals in the Voronoi volume around
//  each atom and a vector of matrices which gives, for each atom,
//  the integrated charges in a sphere vs. the radius of the sphere.
namespace pflow {
  void GetChgInt(vector<aurostd::matrix<double> >& rad_chg_int,  //CO20200404 pflow::matrix()->aurostd::matrix()
      aurostd::matrix<double>& vor_chg_int,  //CO20200404 pflow::matrix()->aurostd::matrix()
      xstructure& str,
      vector<int>& ngrid,
      vector<double>& chg_tot,
      vector<double>& chg_diff) {
    // Loop over grid and bin charge densities.

    // Radial charge integration.
    // Integral over a sphere will be done by setting up NRBIN bins, each
    // of witdh DR, our to RMAX.  Then each charge at each point will be put is the
    // appropriate bin for each atom.  Note that the points and all possible relevant
    // images must be looped over.
    double DR=0.05; // width of radial bins.
    double RMAX=3; // Max distance to bin to.
    int NRBIN=(int)(RMAX/DR); // First bin is 0->DR, last bin is RMAX-DR->RMAX
    aurostd::matrix<double> dum_mat(NRBIN,5,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
    int natoms=pflow::GetNumAtoms(str);
    int npts=ngrid[0]*ngrid[1]*ngrid[2];
    rad_chg_int = vector<aurostd::matrix<double> > (natoms,dum_mat); //CO20200404 pflow::matrix()->aurostd::matrix()
    vector<int> ig(3);
    // Determine the number of images along each lattice vector which must be considered.
    // You will center the charge density on each atom so the max number of images along a
    // given lattice param should be no more than however many lattice params are needed
    // to make sure that you are more than RMAX from an atom along that lattice param direction.
    vector<int> imax(3);
    double scale=pflow::GetScale(str);
    aurostd::matrix<double> lat=pflow::GetLat(str);  //CO20200404 pflow::matrix()->aurostd::matrix()
    for(int ic=0;ic<3;ic++) {
      imax[ic]=(int)(RMAX/(scale*norm(lat[ic]))+1)*ngrid[ic];
    }
    //   double sumchg=0;   //DM not used
    // Atom loop
    aurostd::matrix<double> at_fpos=pflow::GetFpos(str); //CO20200404 pflow::matrix()->aurostd::matrix()
    for(int iat=0;iat<natoms;iat++) {
      // Initialize radial positions in rad_chg_int    
      for(int ib=0;ib<NRBIN;ib++) {
        rad_chg_int[iat][ib][0]=(ib+1)*DR;
      }
      // Find location of atom on grid.
      vector<int> at_grid_pos(3);
      for(int ic=0;ic<3;ic++) {
        at_grid_pos[ic]=Nint(at_fpos[iat][ic]*ngrid[ic]);
      }
      // Grid loop: This loops over grid points shifted to center around grid point approximation to
      // position of iat.  The ig values therefore need to be modified to be true grid points
      // (true_ig) by undoing the shift (adding back the at_grid_pos).
      for(ig[0]=-imax[0];ig[0]<imax[0];ig[0]++) {
        for(ig[1]=-imax[1];ig[1]<imax[1];ig[1]++) {
          for(ig[2]=-imax[2];ig[2]<imax[2];ig[2]++) {
            // Get true grid points
            vector<int> true_ig(3);
            vector<double> chg_fpos(3);
            for(int ic=0;ic<3;ic++) {
              true_ig[ic]=ig[ic]+at_grid_pos[ic];
              // Get true position of this grid point in direct coords.
              chg_fpos[ic]=((float)true_ig[ic]/(float)ngrid[ic]);
              // Shift true grid point back into the 000 cell.
              true_ig[ic]=true_ig[ic]-(int)((float)true_ig[ic]/(float)ngrid[ic])*ngrid[ic];
              if(true_ig[ic]<0) true_ig[ic]=true_ig[ic]+ngrid[ic];
            } // ic
            // Get charge values for each grid position
            int id=true_ig[0]+true_ig[1]*ngrid[0]+true_ig[2]*ngrid[0]*ngrid[1];
            double chgtot=chg_tot[id];
            double chgdiff=chg_diff[id];
            double chgup=0.5*(chg_tot[id]+chg_diff[id]);
            double chgdn=0.5*(chg_tot[id]-chg_diff[id]);
            vector<double> fpos(3);
            // Get distance from iat to true grid point position.
            double dist=pflow::GetDistToAtomDir(str,chg_fpos,iat);
            // Get bin for this atom and increment radial charge density.
            int ibin=(int)(dist/DR);	  
            if(ibin<NRBIN) {
              rad_chg_int[iat][ibin][1]=rad_chg_int[iat][ibin][1]+chgtot;
              rad_chg_int[iat][ibin][2]=rad_chg_int[iat][ibin][2]+chgdiff;
              rad_chg_int[iat][ibin][3]=rad_chg_int[iat][ibin][3]+chgup;
              rad_chg_int[iat][ibin][4]=rad_chg_int[iat][ibin][4]+chgdn;
            } // if ibin<NRBIN
          } // ig[2]
        } // ig[1]
      } // ig[0]
    } // iat

    // Voronoi charge integration
    vor_chg_int = aurostd::matrix<double> (natoms,4,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
    for(ig[0]=0;ig[0]<ngrid[0];ig[0]++) {
      for(ig[1]=0;ig[1]<ngrid[1];ig[1]++) {
        for(ig[2]=0;ig[2]<ngrid[2];ig[2]++) {
          // Get charge values for this grid position
          int id=ig[0]+ig[1]*ngrid[0]+ig[2]*ngrid[0]*ngrid[1];
          double chgtot=chg_tot[id];
          double chgdiff=chg_diff[id];
          double chgup=0.5*(chg_tot[id]+chg_diff[id]);
          double chgdn=0.5*(chg_tot[id]-chg_diff[id]);
          vector<double> fpos(3);
          // Set up grid positions in direct coordinates.
          for(int ic=0;ic<3;ic++) {
            fpos[ic]=((float)ig[ic])/(float)ngrid[ic];
          } // ic
          // Loop over atoms
          double dist;
          double min_dist=pflow::GetDistToAtomImageDir(str,fpos,0);
          int min_dist_id=0;
          for(int iat=0;iat<natoms;iat++) {
            // Get "image" distance to each atom (distance to closest image)
            dist=pflow::GetDistToAtomImageDir(str,fpos,iat);
            // Keep track of which atom is the minimum "image" distance away
            if(dist<min_dist) {
              min_dist=dist;
              min_dist_id=iat;
            } // if dist<min_dist
          } // iat
          vor_chg_int[min_dist_id][0]=vor_chg_int[min_dist_id][0]+chgtot;
          vor_chg_int[min_dist_id][1]=vor_chg_int[min_dist_id][1]+chgdiff;
          vor_chg_int[min_dist_id][2]=vor_chg_int[min_dist_id][2]+chgup;
          vor_chg_int[min_dist_id][3]=vor_chg_int[min_dist_id][3]+chgdn;
        } // ig[2]
      } // ig[1]
    } // ig[0]

    // Normalize and turn rad_chg_int from chg in each radial shell to cumulative integral.
    for(int iat=0;iat<natoms;iat++) {
      for(int i=0;i<4;i++) {
        vor_chg_int[iat][i]=vor_chg_int[iat][i]/npts;
      }
      for(int i=1;i<5;i++) {
        for(int ib=0;ib<NRBIN;ib++) {
          rad_chg_int[iat][ib][i]=rad_chg_int[iat][ib][i]/npts;
          if(ib>0) rad_chg_int[iat][ib][i]=(rad_chg_int[iat][ib][i]+rad_chg_int[iat][ib-1][i]);
        } // ib
      } // for loop over tot,diff,up,dn
    } // iat
  }
}

// ***************************************************************************
// ReadPlaneDensParams
// ***************************************************************************
// Reads in parameters.
// Input file format:

// D # Coordinates for following points (Direct/Cartesian)
// scale # Scale factor - edges of plane get mult. by this (but not origin).
// x y z # origin point
// x y z # "X" axis
// x y z # "Y" axis
// Nx Ny # Number of X and Y grid points
// Middle # Location for origin (Middle/Corner).
// Ortho # Whether to use Y orthogonal to X (Ortho/Strict).
namespace pflow {

  void ReadPlaneDensParams(const xstructure& str, pd_params& pdp, istream& infile) {
    // Note that this converts everything to direct coordinates.
    string s;
    double TOL=1e-7;
    pdp.pts = aurostd::matrix<double> (3,3); //CO20200404 pflow::matrix()->aurostd::matrix()
    pdp.dpts = aurostd::matrix<double> (3,3);  //CO20200404 pflow::matrix()->aurostd::matrix()
    infile >> pdp.type;
    getline(infile,s);
    infile >> pdp.scale;
    getline(infile,s);
    for(int ic=0;ic<3;ic++) {
      for(int jc=0;jc<3;jc++) {
        infile >> pdp.pts[ic][jc];
      }
      getline(infile,s);
      // transform points to cart coordinates if necesary.
      if(pdp.type[0]=='D' || pdp.type[0]=='d') pdp.pts[ic]=pflow::vecF2C(pflow::GetLat(str),pdp.pts[ic]);
    }
    infile >> pdp.Nx >> pdp.Ny;
    getline(infile,s);
    infile >> pdp.orig_loc;
    getline(infile,s);
    infile >> pdp.ortho;
    getline(infile,s);
    pdp.dpts[0]=pdp.pts[0];
    pdp.dpts[1]=VVdiff(pdp.pts[1],pdp.pts[0]);
    pdp.dpts[2]=VVdiff(pdp.pts[2],pdp.pts[0]);
    // Orthogonalize if necessary
    // Finds vector of same length as dpts[2] but perp. to dpts[1].
    if(pdp.ortho[0]=='O' || pdp.ortho[0]=='o') {
      double norm_orig=norm(pdp.dpts[2]);
      pdp.dpts[2]=pflow::CompAperpB(pdp.dpts[2],pdp.dpts[1]);
      double norm_new=norm(pdp.dpts[2]);
      if(norm_new<TOL) {
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"two axis are parallel and no plane can be determined",_INPUT_ILLEGAL_); //CO20200624
      }
      pdp.dpts[2]=SVprod(norm_orig/norm_new,pdp.dpts[2]);
    }
    // Shift origin to middle if necessary.
    // I always place dpts[0] at the corner.  So to put origin
    // at the middle shift dpts[0] by -0.5*(dpts[1]+dpts[2]).
    if(pdp.orig_loc[0]=='M' || pdp.orig_loc[0]=='m') {
      vector<double> shift = pflow::SVprod(0.5,VVsum(pdp.dpts[1],pdp.dpts[2]));
      pdp.dpts[0]=pflow::VVdiff(pdp.dpts[0],shift);
    }
    // Convert everything to direct coordinates
    for(int ic=0;ic<3;ic++) {
      pdp.pts[ic]=pflow::vecC2F(pflow::GetLat(str),pdp.pts[ic]);
      pdp.dpts[ic]=pflow::vecC2F(pflow::GetLat(str),pdp.dpts[ic]);
      // rescale coordinates
      if(ic>0) pdp.dpts[ic]=pflow::SVprod(pdp.scale,pdp.dpts[ic]);
    }
  }
}

// ***************************************************************************
// GetPlaneDens
// ***************************************************************************
// Gets density in plane.  Work always within direct coordinates.
// The original densities are stored in chg_tot, chg_diff such that
// if the direct coordinates of a grid point are i/Nx,j/Ny,k/Nz then
// the index for the charge at that point is id=i+j*Nx+k*Nx*Ny.  The
// charge in the plane will be stored the same way.  Let Npx,Npy =
// plane grid.  Then if the coordinates of a grid point in the plane are
// i/Npx, j/Npy then the index for the charge at that point in id=i+j*Npx.
// You need to be able to access charge densities and structural
// parameters for this routine.
namespace pflow {
  void GetPlaneDens(const pd_params& pdp,
      vector<double>& dens2d_tot,
      vector<double>& dens2d_diff,
      const xstructure& str,
      const vector<int>& ngrid,
      const vector<double>& chg_tot,
      const vector<double>& chg_diff) {
    if(str.atoms.size()) {;} // phony just to keep str busy
    int Nx=pdp.Nx;
    int Ny=pdp.Ny;
    int npts=Nx*Ny;
    dens2d_tot=vector<double> (npts,0.0);
    dens2d_diff=vector<double> (npts,0.0);
    // Loop over each point in the 2d plane
    for(int i=0;i<Nx;i++) {
      for(int j=0;j<Ny;j++) {
        // Get id2d for this point.
        int id2d =i+Nx*j;
        // Get direct coordinates for this point.
        vector<double> pt(3,0.0);
        vector<double> Xpt(3,0.0);
        vector<double> Ypt(3,0.0);
        Xpt=SVprod(float(i)/float(Nx-1),pdp.dpts[1]);
        Ypt=SVprod(float(j)/float(Ny-1),pdp.dpts[2]);
        pt=VVsum(Xpt,Ypt);
        pt=VVsum(pdp.dpts[0],pt);
        // Get charge at this point.
        // Simply use value of charge at nearest point
        // on the 3d charge grid.
        int ii=Nint(ngrid[0]*pt[0])%ngrid[0];
        int jj=Nint(ngrid[1]*pt[1])%ngrid[1];
        int kk=Nint(ngrid[2]*pt[2])%ngrid[2];
        if(ii<0) ii=ii+ngrid[0];
        if(jj<0) jj=jj+ngrid[1];
        if(kk<0) kk=kk+ngrid[2];
        int id3d = ii+ngrid[0]*jj+ngrid[0]*ngrid[1]*kk;
        dens2d_tot[id2d]=chg_tot[id3d];
        dens2d_diff[id2d]=chg_diff[id3d];
      }
    }  
  }
}

// ***************************************************************************
// PrintPlaneDens
// ***************************************************************************
namespace pflow {
  void PrintPlaneDens(const pd_params& pdp,
      const vector<double>& dens2d_tot,
      const vector<double>& dens2d_diff,
      const xstructure& str) {  
    if(str.atoms.size()) {;} // phony just to keep str busy

    ofstream outf_tot("dens2d.tot.out");
    ofstream outf_up("dens2d.up.out");
    ofstream outf_dn("dens2d.dn.out");
    ofstream outf_diff("dens2d.diff.out");

    int Nx=pdp.Nx;
    int Ny=pdp.Ny;
    for(int i=0;i<Nx;i++) {
      for(int j=0;j<Ny;j++) {
        int id=i+Nx*j;
        outf_tot << dens2d_tot[id] << " ";
        outf_diff << dens2d_diff[id] << " ";
        outf_up << (dens2d_tot[id]+dens2d_diff[id])/2.0 << " ";
        outf_dn << (dens2d_tot[id]-dens2d_diff[id])/2.0 << " ";
      }
      outf_tot << endl;
      outf_diff << endl;
      outf_up << endl;
      outf_dn << endl;
    }
  }
}

// ***************************************************************************
// EWALD FUNCS EWALD FUNCS EWALD FUNCS EWALD FUNCS EWALD FUNCS EWALD FUNCS EWA
// ***************************************************************************
double CONV_FACT=(1.0E+10*E_ELECTRON/(4.0*PI*EPS_VACUUM));

// ***************************************************************************
// Ewald
// ***************************************************************************
// Calculates Ewald Sum.  Based on routine from Eric Wu.

// ******************************************************************
// **                   Subroutine ewald                           **
// **                                                              **
// **   Subroutine ewald calculates the ewald sum in eV            **
// **   Also, this calculates the ewald gradients in ev   since    **
// **         it is just a small calculation to do- ewald scales   **
// **         as N^2, extra force calculation  adds as  N          **
// **                                                              **
// **   NB if you want to convert forces to ev/A you need the      **
// **   inverse of the lattice vectors matrix (A) then do a        **
// **    Fx=ainvx*fx+ainvy*fy+ainvz*fz                             **
// **                 where Fx is in ev/A fx is reduced            **
// **   to get the components along each cooredinate (I think)     **
// **                                                              **
// **                                                              **
// **   References:                                                **
// **     www.dl.ac.uk/TCS/Software/DL_POLY/USRMAN/node62.html     **
// **     www.ee.duke.edu/~ayt/ewaldpaper/node8.html               **
// **     Moldy (molecular dynamics) manual                        **
// **   http://www.earth.ox.ac.uk/~keith/moldy-manual/node11.html  **
// **     Allen+tildesly p159 (for cubic)                          **
// **     Kittel Appendix B                                        **
// **                                                              **
// **   Formula:                                                   **
// **                                                              **
// **   For the energy:                                            **
// **      Etotal = Erecip + Ereal + Epoint                        **
// **                                                              **
// **      Epoint = [1/(4*pi*eo)]                                  **
// **                   * sqrt(eta)/sqrt(pi) * sum(i) qi^2         **
// **              +[1/(4*pi*e0)] * [pi/(2*vol*eta)]               **
// **                   * [sum(i) qi]^2   (for charged cell)       **
// **      Ereal = [1/(4*pi*e0)]* 0.5 *                            **
// **                      sum(ij) (qi*qj*erfc(sqrt(eta)*rij)/rij  **      
// **      Erecip = [1/(4*pi*e0)] *0.5* [4*pi/vol] *               **
// **               sum(G) 1/G^2 * exp[-G^2/(4*eta)] *             **
// **               |sum(i) (qi*exp(-iGri)|^2                      **
// **                                                              **
// **                                                              **
// **  For the forces:                                             **
// **      procedure: take the neg. derivatives.                   **
// **                 The point is const,                          **
// **                 and disappears.  So                          **
// **                                                              **
// **      Fj=  Fjrecip + Fjreal                                   **
// **                                                              **
// **      Fjxrecip = +[1/(4*pi*e0)] * [4*pi/vol] *                **
// **                 qj*  sum(G) Gx * G*1/G^2 *                   **
// **                  exp[-G^2/(4*eta)]*                          **
// **                             [-sin(Grj)*sum(i) qicos(Gri)]    **
// **                            +[cos(Grj)]*sum(i) qisin(Gri)]    **
// **       to get this, recall exp(-iGr)=cos(Gr)+isin(Gr)         **
// **       and put in Erecip.  Then take the derivatives of       **
// **       cos^2(Gr) and sin^2(Gr)                                **
// **                                                              **
// **      Fjxreal =   -[1/(4*pi*e0)]  * qj *                      **
// **                  Sum(i) qi/rij^3*                            **
// **                         [2*sqrt(eta)*rij/sqrt(pi)            **
// **                                  *exp(-eta*rij^2)            **
// **                          +erfc(sqrt(eta)*rij)]*rxij          **
// **                recall                                        **
// **        erf(x)=2/sqrt(pi)integral(0,x)exp(-t^2)dt             **
// **                                                              **
// **                                                              **
// **   input variables                                            **
// **                                                              **
// **     lat(3,3) real space lattice vectors input (ax ay az      **
// **                                              bx by bz        **
// **                                              cx cy cz)       **
// **                                               but stored     **
// **                                       ax bx cx...         .  **
// **                                       as per fort convention **
// **     rlat(3,3) recip space lattice vectors (see above)        **
// **     natoms  number of atoms                                  **
// **     ntype   # of types of atoms                              **
// **     vol     cell volume (A^3)                                **
// **     atfpos(natoms,3)   position of atoms in fractional       **
// **     atchg(natoms)      nominal charge of atoms by TYPE       **
// **     attyp(natoms)        type of atom (type of at 1, etc.    **
// **                                                              **
// **   output variables                                           **
// **                                                              **
// **     eewalde   ewald energy in   eV                           **
// **     ewaldf(3,natoms)   negative derivative of ewald          **
// **                        energy                                **
// **                          aka ewald forces ev/A               **
// ******************************************************************

namespace pflow {
  void Ewald(const xstructure& in_str,double& epoint,double& ereal,
      double& erecip,double& eewald,double& eta,const double& SUMTOL) {
    xstructure str=in_str;
    str=PutInCell(str);
    aurostd::matrix<double> lat = pflow::GetScaledLat(str);  //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> rlat = pflow::RecipLat(lat); //CO20200404 pflow::matrix()->aurostd::matrix()
    int natoms=pflow::GetNumAtoms(str);
    double vol=pflow::GetVol(lat);
    aurostd::matrix<double> atfpos=pflow::GetFpos(str);  //CO20200404 pflow::matrix()->aurostd::matrix()
    vector<int> attyp=GetTypes(str);
    vector<string> names=GetNames(str);
    // Charges should be stored in names - turn into vector of doubles.
    vector<double> atchg(natoms);
    double totchg=0;
    for(int i=0;i<(int) names.size();i++) {
      atchg[i]=atof(names[i].c_str());
      totchg=totchg+atchg[i];
    }

    // Constants
    double TOL=1e-6;

    // Check chg neutral
    if(abs(totchg)>TOL) {
      cerr << "WARNING: Cell is not neutral" << endl;
      cerr << "WARNING: Total charge = " << totchg << endl;
      cerr << "WARNING: Convasp will use jellium background to force neutrality."<< endl;
    }

    // Get eta if it was not input (it is <0 if it needs to be found).
    if(eta<=0) {eta=GetEta(natoms,vol);}
    // Constants
    double rteta=sqrt(eta);
    // Point energy
    epoint=GetPointEner(rteta,atchg,vol);
    // Recipricol energy
    erecip=GetRecipEner(eta,atchg,vol,rlat,atfpos,SUMTOL);
    // Real space energy
    ereal=GetRealEner(eta,atchg,vol,lat,atfpos,SUMTOL);
    // Sum all terms
    eewald=(ereal+erecip-epoint);
  }
}

// ***************************************************************************
// GetEta
// ***************************************************************************
namespace pflow {
  double GetEta(const int& natoms,const double& vol) {
    double term=std::pow((double) (5.5*natoms/vol/vol),(double) 1.0/6.0);
    double eta=sqrt(PI)/1.2*term;
    eta=eta*eta;
    return eta;
  }
}

// ***************************************************************************
// GetPointEner
// ***************************************************************************
// Gets point energy.
namespace pflow {
  double GetPointEner(const double& rteta, const vector<double>& atchg,
      const double& vol) {
    double ept=0;
    double chg_cell_corr=0;
    for(int ia=0;ia<(int) atchg.size();ia++) {
      ept=ept+atchg[ia]*atchg[ia];
      chg_cell_corr=chg_cell_corr+atchg[ia];
    }
    ept=ept*rteta/RTPI;
    chg_cell_corr=chg_cell_corr*chg_cell_corr*PI/(2.0*vol*rteta*rteta); //chg_cell_corr*PI/(2.0*vol*rteta*rteta); //Wei Xie 20180316
    ept=ept+chg_cell_corr;
    ept=CONV_FACT*ept;
    return ept;
  }
}

// ***************************************************************************
// GetRecipEner
// ***************************************************************************
// Gets recipricol energy.
namespace pflow {
  double GetRecipEner(const double& eta,const vector<double>& atchg,
      const double& vol,const aurostd::matrix<double>& rlat, //CO20200404 pflow::matrix()->aurostd::matrix()
      const aurostd::matrix<double>& atfpos,const double& SUMTOL) {  //CO20200404 pflow::matrix()->aurostd::matrix()
    double TOL=1e-16;
    double log_eps=-30; // This assumes sf<
    double arg=0;
    double erecip=0;
    int gcnt=1;
    int gcont=1; // whether or not to add more shells.
    while(gcont) {
      double maxterm=0;
      vector<int> ig(3);
      for(ig[0]=-gcnt;ig[0]<=gcnt;ig[0]++) {
        for(ig[1]=-gcnt;ig[1]<=gcnt;ig[1]++) {
          for(ig[2]=-gcnt;ig[2]<=gcnt;ig[2]++) {
            if( (abs(ig[0])==gcnt) || (abs(ig[1])==gcnt) ||
                (abs(ig[2])==gcnt) || gcnt==1) { // Just does new shells
              vector<double> gvec=pflow::VMmult(ig,rlat);
              double gsq=pflow::VVdot(gvec,gvec);
              if(gsq>TOL) { // Avoid gsq=0.
                // get structure factor
                double sfr=0;
                double sfi=0;
                for(int ia=0;ia<(int) atfpos.size();ia++) {
                  double exparg=TWOPI*pflow::VVprod(atfpos[ia],ig);
                  sfr=sfr+atchg[ia]*cos(exparg);
                  sfi=sfi+atchg[ia]*sin(exparg);
                }
                if(aurostd::abs(sfr)<1e-16) {sfr=0.0;}
                if(aurostd::abs(sfi)<1e-16) {sfi=0.0;}
                double sf=sfr*sfr+sfi*sfi;
                // Get expval.  In order not to get floating point
                // errors due to manipulating small numbers we must
                // not evaluate the exp when it produces irrelevantly
                // small arguments.  This happens when
                // expval*sf << eps (where safe eps=1e-30), or equivalently,
                // when ln(sf)-ln(g^2)-g^2/(4eta) < ln(eps)
                arg=gsq/(4.0*eta);
                double expval;
                double term1=0;
                if(sf>0) {
                  term1=log(sf)-log(gsq)-gsq/(4*eta);
                }
                if(term1<log_eps) { // Avoids some floating point exceptions.
                  expval=0;
                }
                else {
                  expval=exp(-arg)/gsq;
                }
                double term=expval*sf;
                //	      if(term<1e-16) {term=0;}
                erecip=erecip+term;
                // Set maxterm to max of abs(term) and expval.
                // The expval term is needed since in some shells
                // term may be very small due to coincidently small sf,
                // even before convergence.
                if(aurostd::abs(term)>maxterm) {maxterm=aurostd::abs(term);}
                if(expval>maxterm) {maxterm=expval;}
                // ADD FORCES ???
              } // If gsq>TOL
            } // If to do only new shells (ig==gcnt)
          } // ig0
        } // ig1
      } // ig2
      gcnt++; // Add new shell
      if(maxterm<SUMTOL) {gcont=0;} // Do not add another shell
    } // while gcont
    erecip=erecip*0.5*4.0*PI/vol; // CGS units, convert to eV at end.
    erecip=CONV_FACT*erecip;
    // FORCES STUFF HERE ???				  
    return erecip;				  
  } // End GetRecipEner
}

// ***************************************************************************
// GetRealEner
// ***************************************************************************
// Gets real space energy.
namespace pflow {
  double GetRealEner(const double& eta,const vector<double>& atchg,
      const double& vol,const aurostd::matrix<double>& lat,  //CO20200404 pflow::matrix()->aurostd::matrix()
      const aurostd::matrix<double>& atfpos,const double& SUMTOL) {  //CO20200404 pflow::matrix()->aurostd::matrix()
    if(vol) {;} // phony just to keep vol busy

    double ereal=0;
    double TOL=1e-5;
    double log_eps=-30;
    double erfcarg=0;
    double rteta=sqrt(eta);
    int nat=atfpos.size();
    int rcnt=1;
    int rcont=1;
    while(rcont) {
      double maxterm=0;
      vector<int> ir(3);
      for(ir[0]=-rcnt;ir[0]<=rcnt;ir[0]++) {
        for(ir[1]=-rcnt;ir[1]<=rcnt;ir[1]++) {
          for(ir[2]=-rcnt;ir[2]<=rcnt;ir[2]++) {
            if( (abs(ir[0])==rcnt) || (abs(ir[1])==rcnt) ||
                (abs(ir[2])==rcnt) || rcnt==1) { // Just does new shells
              for(int ia=0;ia<nat;ia++) { // All atoms in unit cell
                for(int ja=0;ja<nat;ja++) { // All neighbors in cell given by ir
                  // Get displacement between atoms ia(cell 0) and ja(cell ir).
                  vector<double> disp=pflow::VVdiff(atfpos[ja],atfpos[ia]);
                  disp=pflow::VVsum(disp,ir);
                  disp=pflow::vecF2C(lat,disp);
                  double dist=pflow::norm(disp);
                  if(dist>TOL) { // Avoid atom dist to itself.
                    erfcarg=rteta*dist;
                    // Get erfcval.  In order not to get floating point
                    // errors due to manipulating small numbers we must
                    // not evaluate the erfc when it produces irrelevantly
                    // small arguments.  This happens when
                    // erfc(x)*abs(q) << eps (where safe eps=1e-30), or equivalently,
                    // when ln(abs(q))-x^2-ln(x)-0.5*ln(pi) < ln(eps) (where I have
                    // used the asymptotic form erfc(x)~e^-x^2/(x*rt(pi))
                    // for large x.
                    double erfcval;
                    double term1=atchg[ia]*atchg[ja]/dist;
                    double term2=0;
                    if(abs(term1)>0) {
                      term2=log(abs(term1))-erfcarg*erfcarg-log(erfcarg)-0.5*log(PI);
                    }
                    if(term2<log_eps) { // Avoids some floating point exceptions.
                      erfcval=0;
                    }
                    else {
                      erfcval=erfc(erfcarg);
                    }
                    // tpx
                    // cout << "term1 erfcval " << term1 << " " << erfcval << endl;
                    double term=term1*erfcval;
                    ereal=ereal+term;
                    if(aurostd::abs(term)>maxterm) {maxterm=aurostd::abs(term);}
                    // FORCES ???
                  } // rdist>TOL
                } // ja
              } // ia
            } // if to do only new shells (ir==rcnt)
          } // rg0
        } // rg1
      } // rg2
      rcnt++; // Add new shell
      if(maxterm<SUMTOL) {rcont=0;} // Do not add another shell
    } // while rcont
    ereal=ereal*0.5; // CGS units, convert to eV at end.
    ereal=CONV_FACT*ereal;
    return ereal;				  
  }  
}

// ***************************************************************************
// GetSreenedESEner
// ***************************************************************************
// Gets electrostatic energy with screening (does sum in real space).
namespace pflow {
  double ScreenedESEner(const xstructure& in_str,const double& Ks,const double& SUMTOL) {
    xstructure str=in_str;
    str=PutInCell(str);
    aurostd::matrix<double> lat = pflow::GetScaledLat(str);  //CO20200404 pflow::matrix()->aurostd::matrix()
    int natoms=pflow::GetNumAtoms(str);
    // double vol=pflow::GetVol(lat);  //DM not used
    aurostd::matrix<double> atfpos=pflow::GetFpos(str);  //CO20200404 pflow::matrix()->aurostd::matrix()
    vector<int> attyp=GetTypes(str);
    vector<string> names=GetNames(str);
    // Charges should be stored in names - turn into vector of doubles.
    vector<double> atchg(natoms);
    double totchg=0;
    for(int i=0;i<(int)names.size();i++) {
      atchg[i]=atof(names[i].c_str());
      totchg=totchg+atchg[i];
    }
    double ereal=0;
    double TOL=1e-5;
    int nat=atfpos.size();
    int rcnt=1;
    int rcont=1;
    while(rcont) {
      double maxterm=0;
      vector<int> ir(3);
      for(ir[0]=-rcnt;ir[0]<=rcnt;ir[0]++) {
        for(ir[1]=-rcnt;ir[1]<=rcnt;ir[1]++) {
          for(ir[2]=-rcnt;ir[2]<=rcnt;ir[2]++) {
            if( (abs(ir[0])==rcnt) || (abs(ir[1])==rcnt) ||
                (abs(ir[2])==rcnt) || rcnt==1) { // Just does new shells
              for(int ia=0;ia<nat;ia++) { // All atoms in unit cell
                for(int ja=0;ja<nat;ja++) { // All neighbors in cell given by ir
                  // Get displacement between atoms ia(cell 0) and ja(cell ir).
                  vector<double> disp=pflow::VVdiff(atfpos[ja],atfpos[ia]);
                  disp=pflow::VVsum(disp,ir);
                  disp=pflow::vecF2C(lat,disp);
                  double dist=pflow::norm(disp);
                  if(dist>TOL) { // Avoid atom dist to itself.
                    double term=exp(-Ks*dist)*atchg[ia]*atchg[ja]/dist;
                    ereal=ereal+term;
                    if(aurostd::abs(term)>maxterm) {maxterm=aurostd::abs(term);}
                    // FORCES ???
                  } // rdist>TOL
                } // ja
              } // ia
            } // if to do only new shells (ir==rcnt)
          } // rg0
        } // rg1
      } // rg2
      rcnt++; // Add new shell
      if(maxterm<SUMTOL) {rcont=0;} // Do not add another shell
    } // while rcont
    ereal=ereal*0.5; // CGS units, convert to eV at end.
    ereal=CONV_FACT*ereal;
    return ereal;				  
  }  
}

namespace pflow {
  void BZMAX(istream & input) {
    //Usage: aflow --BZmax < POSCAR

    xstructure str_in(input,IOAFLOW_AUTO);
    xstructure str_sp,str_sc;
    LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc); //This functions format the str_in, and stores data in str_sp and str_sc
    double grid=16.0;
    string lattice, stmp, line, kpath_string;
    stringstream KPATH, strline;
    bool foundBZ;
    xvector<double> b1(3); xvector<double> b2(3); xvector<double> b3(3); //reciprocal lattice
    xmatrix<double> klattice;    //Declare a matrix to store reciprocal lattice

    lattice=str_sp.bravais_lattice_variation_type;
    kpath_string=LATTICE::KPOINTS_Directions(lattice, str_sp.lattice, grid,str_sp.iomode,foundBZ);
    KPATH.str(kpath_string);  //converting kpath from string type to stringstream
    klattice=str_sp.klattice;

    b1[1]=str_sp.klattice(1,1); b1[2]=str_sp.klattice(1,2); b1[3]= str_sp.klattice(1,3);
    b2[1]=str_sp.klattice(2,1); b2[2]=str_sp.klattice(2,2); b2[3]= str_sp.klattice(2,3);
    b3[1]=str_sp.klattice(3,1); b3[2]=str_sp.klattice(3,2); b3[3]= str_sp.klattice(3,3);

    //////////////READING KPOINTS//////////////////////////////////////////////////////////
    for (int i=0; i<4; i++) {getline(KPATH, stmp);} //Read the first 4 lines
    vector<vector<double> > kpoints;
    vector<string> kpointslabel;
    int count=0, j=0; //count is the number of rows of kpoints
    while (getline(KPATH, line)) {
      if(aurostd::CountWordsinString(line)>=3) {
        vector<double> kpt_tmp;
        double a1, a2, a3;
        string a5;
        strline.clear();
        strline.str(line);
        strline >> a1 >> a2 >> a3 >> stmp >> a5;
        kpt_tmp.push_back(a1);
        kpt_tmp.push_back(a2);
        kpt_tmp.push_back(a3);
        kpoints.push_back(kpt_tmp);
        kpointslabel.push_back(a5);
        j++;
      }
      count =j;
    }

    //////////////READING KPOINTS//////////////////////////////////////////////////////////
    //Converting kpoints into Cartesian
    int NKPOINTS=count;
    vector<vector<double> > kpcart;
    vector<double> klinedistance;
    //vector<double> klinedistance(NKPOINTS);
    for (int i=0; i<NKPOINTS; i++) {
      vector<double> kpcart_tmp;
      double c1=0, c2=0, c3=0, distance;
      c1=kpoints[i][0]*b1[1]+kpoints[i][1]*b2[1]+kpoints[i][2]*b3[1];
      c2=kpoints[i][0]*b1[2]+kpoints[i][1]*b2[2]+kpoints[i][2]*b3[2];
      c3=kpoints[i][0]*b1[3]+kpoints[i][1]*b2[3]+kpoints[i][2]*b3[3];
      kpcart_tmp.push_back(c1);
      kpcart_tmp.push_back(c2);
      kpcart_tmp.push_back(c3);
      kpcart.push_back(kpcart_tmp);
      distance=sqrt(c1*c1+c2*c2+c3*c3);
      klinedistance.push_back(distance);
    }
    //Sorting
    for (int i=0; i<(NKPOINTS-1); i++) {
      for (int j=i+1; j<NKPOINTS; j++) {
        if(klinedistance[j] > klinedistance[i]) {
          //Sorting kline distance
          double tmp;
          tmp=klinedistance[i];
          klinedistance[i]=klinedistance[j];
          klinedistance[j]=tmp;
          //Meanwhile change the order of the kpoints label
          string kplabel_tmp;
          kplabel_tmp=kpointslabel[i];
          kpointslabel[i]=kpointslabel[j];
          kpointslabel[j]=kplabel_tmp;
          //Meanwhile change the order of the kpoints label
          double k0, k1, k2;
          k0=kpoints[i][0];
          k1=kpoints[i][1];
          k2=kpoints[i][3];
          kpoints[i][0]=kpoints[j][0];
          kpoints[i][1]=kpoints[j][1];
          kpoints[i][2]=kpoints[j][2];
          kpoints[j][0]=k0;
          kpoints[j][1]=k1;
          kpoints[j][2]=k2;
        }
      }
    }
    //Formating output
    cout.precision(8);
    cout << "Bravais Lattice:   " << lattice << endl;
    cout << "              KPOINTS              DISTANCE   TYPE" << endl;
    cout << kpoints[0][0] <<" ";
    cout << kpoints[0][1] <<" ";
    cout << kpoints[0][2]<<"  ";
    cout << klinedistance[0] <<"   ";
    cout << kpointslabel[0] << endl;

    for (int i=1; i<(NKPOINTS-1); i++) {
      if(kpointslabel[i].compare(kpointslabel[i-1])!=0) {
        cout << kpoints[i][0] <<" ";
        cout << kpoints[i][1] <<" ";
        cout << kpoints[i][2]<<"  ";
        cout << klinedistance[i] <<"   ";
        cout << kpointslabel[i] << endl;
      }
    }
  }
}


//ME20190628 BEGIN
// prettyPrintCompound
namespace pflow {

  static const double ZERO_TOL = 1E-8;  // from CHULL

  string prettyPrintCompound(const string& compound, vector_reduction_type vred, bool exclude1, filetype ftype) {  //char mode  //CO20190629
    vector<double> vcomposition;
    vector<string> vspecies =  aurostd::getElements(compound, vcomposition);
    return prettyPrintCompound(vspecies, vcomposition, vred, exclude1, ftype);  //mode  //CO20190629
  }

  string prettyPrintCompound(const vector<string>& vspecies,const vector<uint>& vcomposition,vector_reduction_type vred,bool exclude1,filetype ftype) {  // overload //char mode //DX20200727
    vector<double> vcomposition_dbl;
    for(uint i=0;i<vcomposition.size();i++){vcomposition_dbl.push_back((double)vcomposition[i]);}
    return prettyPrintCompound(vspecies,aurostd::vector2xvector<double>(vcomposition_dbl),vred,exclude1,ftype); //mode //CO20190629
  }

  // Moved here from the ConvexHull class
  string prettyPrintCompound(const vector<string>& vspecies,const vector<double>& vcomposition,vector_reduction_type vred,bool exclude1,filetype ftype) {  // overload //char mode //CO20190629
    return prettyPrintCompound(vspecies,aurostd::vector2xvector<double>(vcomposition),vred,exclude1,ftype); //mode //CO20190629
  }

  string prettyPrintCompound(const vector<string>& vspecies,const xvector<uint>& vcomposition,vector_reduction_type vred,bool exclude1,filetype ftype) {  // overload //char mode //DX20200727
    xvector<double> vcomposition_dbl(vcomposition.rows);
    for(int i=1;i<=vcomposition.rows;i++){vcomposition_dbl(i)=(double)vcomposition[i];}
    return prettyPrintCompound(vspecies,vcomposition_dbl,vred,exclude1,ftype); //mode //CO20190629
  }

  string prettyPrintCompound(const vector<string>& vspecies,const xvector<double>& vcomposition,vector_reduction_type vred,bool exclude1,filetype ftype) {  // main function //char mode //CO20190629
    // 2-D, we usually want vred=gcd_vrt true for convex points, and no_vrt elsewhere
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;
    uint precision=COEF_PRECISION;
    stringstream output;output.precision(precision);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " vspecies=" << aurostd::joinWDelimiter(vspecies,",") << endl;
      cerr << __AFLOW_FUNC__ << " vcomposition=" << aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(vcomposition),",") << endl;
    }
    if(vspecies.size()!=(uint)vcomposition.rows) {
      message << "vspecies.size() != vcomposition.rows" << endl;
      message << "vspecies=" << aurostd::joinWDelimiter(vspecies,",") << endl;
      message << "vcomposition=" << aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(vcomposition),",") << endl;
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message, _INDEX_MISMATCH_);
    }
    // special case, unary
    if(vspecies.size() == 1) {
      output << vspecies[0];
      if(!exclude1) {output << (vred==gcd_vrt?1:vcomposition[vcomposition.lrows]);}
      return output.str();
    }
    xvector<double> comp=vcomposition;
    xvector<double> final_comp=comp; //DX20191125
    //DX20191125 [OBSOLETE] if(vred==gcd_vrt){comp=aurostd::reduceByGCD(comp,ZERO_TOL);}
    if(vred==gcd_vrt){aurostd::reduceByGCD(comp,final_comp,ZERO_TOL);} //DX20191125 - new function form
    else if(vred==frac_vrt){final_comp=aurostd::normalizeSumToOne(comp,ZERO_TOL);}
    else if(vred==no_vrt){;}
    else {throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__,"Unknown reduce mode",_INPUT_UNKNOWN_);}
    if(std::abs(aurostd::sum(final_comp)) < ZERO_TOL){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Empty composition");}
    for(uint i=0,fl_size_i=vspecies.size();i<fl_size_i;i++) {
      output << vspecies[i];
      if(!(exclude1 && aurostd::identical(final_comp[i+final_comp.lrows],1.0,ZERO_TOL))) {
        if(ftype==latex_ft) {output << "$_{"; //mode==_latex_ //CO20190629
        } else if(ftype==gnuplot_ft){output<< "_{";}  //mode==_gnuplot_ //CO20190629
        output << final_comp[i+final_comp.lrows];
        if(ftype==latex_ft) {output << "}$";} //mode==_latex_ //CO20190629
        else if(ftype==gnuplot_ft){output<< "}";} //mode==_gnuplot_ //CO20190629
      }
    }
    return output.str();
  }
}

//[CO20200526 - EASY TEMPLATE CLASS]namespace pflow {
//[CO20200526 - EASY TEMPLATE CLASS]  AQueue::AQueue(ostream& oss) : xStream(oss),m_initialized(false) {initialize();}
//[CO20200526 - EASY TEMPLATE CLASS]  AQueue::AQueue(ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize();}
//[CO20200526 - EASY TEMPLATE CLASS]  AQueue::AQueue(const AQueue& b) : xStream(*b.getOFStream(),*b.getOSS()) {copy(b);} // copy PUBLIC
//[CO20200526 - EASY TEMPLATE CLASS]  
//[CO20200526 - EASY TEMPLATE CLASS]  AQueue::~AQueue() {xStream::free();free();}
//[CO20200526 - EASY TEMPLATE CLASS]
//[CO20200526 - EASY TEMPLATE CLASS]  const AQueue& AQueue::operator=(const AQueue& other) {
//[CO20200526 - EASY TEMPLATE CLASS]    if(this!=&other) {copy(other);}
//[CO20200526 - EASY TEMPLATE CLASS]    return *this;
//[CO20200526 - EASY TEMPLATE CLASS]  }
//[CO20200526 - EASY TEMPLATE CLASS]
//[CO20200526 - EASY TEMPLATE CLASS]  void AQueue::clear() {free();}  //clear PUBLIC
//[CO20200526 - EASY TEMPLATE CLASS]  void AQueue::free() {
//[CO20200526 - EASY TEMPLATE CLASS]    m_initialized=false;
//[CO20200526 - EASY TEMPLATE CLASS]  }
//[CO20200526 - EASY TEMPLATE CLASS]  
//[CO20200526 - EASY TEMPLATE CLASS]  void AQueue::copy(const AQueue& b) {  //copy PRIVATE
//[CO20200526 - EASY TEMPLATE CLASS]    m_initialized=b.m_initialized;
//[CO20200526 - EASY TEMPLATE CLASS]  }
//[CO20200526 - EASY TEMPLATE CLASS]  
//[CO20200526 - EASY TEMPLATE CLASS]  bool AQueue::initialize(ostream& oss) {
//[CO20200526 - EASY TEMPLATE CLASS]    xStream::initialize(oss);
//[CO20200526 - EASY TEMPLATE CLASS]    return initialize();
//[CO20200526 - EASY TEMPLATE CLASS]  }
//[CO20200526 - EASY TEMPLATE CLASS]
//[CO20200526 - EASY TEMPLATE CLASS]  bool AQueue::initialize(ofstream& FileMESSAGE,ostream& oss) {
//[CO20200526 - EASY TEMPLATE CLASS]    xStream::initialize(FileMESSAGE,oss);
//[CO20200526 - EASY TEMPLATE CLASS]    return initialize();
//[CO20200526 - EASY TEMPLATE CLASS]  }
//[CO20200526 - EASY TEMPLATE CLASS]
//[CO20200526 - EASY TEMPLATE CLASS]  bool AQueue::initialize() {
//[CO20200526 - EASY TEMPLATE CLASS]    free();
//[CO20200526 - EASY TEMPLATE CLASS]    return false;
//[CO20200526 - EASY TEMPLATE CLASS]  }
//[CO20200526 - EASY TEMPLATE CLASS]}

#define _AQUEUE_DEBUG_ false

//CO20200526 - queueing class
namespace pflow {
  void AJob::free(){
    m_index=AUROSTD_MAX_UINT;
    m_id=AUROSTD_MAX_UINT;
    m_user="";
    m_status=JOB_QUEUED;
    m_ncpus=0;  //total ncpus specified for the job
    m_vinodes.clear();
    m_vncpus.clear(); //ncpus allocated to each node (sum of this should be m_ncpus)
    m_vipartitions.clear();
  }
}

namespace pflow {
  void ANode::free(){
    m_index=AUROSTD_MAX_UINT;
    m_name="";
    m_status=NODE_OFFLINE;
    m_ncpus=0;
    m_ncpus_occupied=0;
    m_properties="";
    m_vijobs.clear();
    m_vipartitions.clear();
  }
  bool ANode::isStatus(const node_status& status) const {
    //[keep silent]bool LDEBUG=(FALSE || _AQUEUE_DEBUG_ || XHOST.DEBUG);
    //[keep silent]if(LDEBUG){
    //[keep silent]  //status
    //[keep silent]  if(status==NODE_FREE){cerr << __AFLOW_FUNC__ << " status==NODE_FREE" << endl;}
    //[keep silent]  else if(status==NODE_OCCUPIED){cerr << __AFLOW_FUNC__ << " status==NODE_OCCUPIED" << endl;}
    //[keep silent]  else if(status==NODE_FULL){cerr << __AFLOW_FUNC__ << " status==NODE_FULL" << endl;}
    //[keep silent]  else if(status==NODE_DOWN){cerr << __AFLOW_FUNC__ << " status==NODE_DOWN" << endl;}
    //[keep silent]  else if(status==NODE_OFFLINE){cerr << __AFLOW_FUNC__ << " status==NODE_OFFLINE" << endl;}
    //[keep silent]  else if(status==NODE_OPERATIONAL){cerr << __AFLOW_FUNC__ << " status==NODE_OPERATIONAL" << endl;}
    //[keep silent]  else if(status==NODE_NONOPERATIONAL){cerr << __AFLOW_FUNC__ << " status==NODE_NONOPERATIONAL" << endl;}
    //[keep silent]  //m_status
    //[keep silent]  if(m_status==NODE_FREE){cerr << __AFLOW_FUNC__ << " m_status==NODE_FREE" << endl;}
    //[keep silent]  else if(m_status==NODE_OCCUPIED){cerr << __AFLOW_FUNC__ << " m_status==NODE_OCCUPIED" << endl;}
    //[keep silent]  else if(m_status==NODE_FULL){cerr << __AFLOW_FUNC__ << " m_status==NODE_FULL" << endl;}
    //[keep silent]  else if(m_status==NODE_DOWN){cerr << __AFLOW_FUNC__ << " m_status==NODE_DOWN" << endl;}
    //[keep silent]  else if(m_status==NODE_OFFLINE){cerr << __AFLOW_FUNC__ << " m_status==NODE_OFFLINE" << endl;}
    //[keep silent]}
    if(status==NODE_OPERATIONAL && (m_status==NODE_FREE||m_status==NODE_OCCUPIED||m_status==NODE_FULL)){return true;}
    else if(status==NODE_NONOPERATIONAL && (m_status==NODE_DOWN||m_status==NODE_OFFLINE)){return true;}
    return m_status==status;
  }
}

namespace pflow {
  void APartition::free(){
    m_index=AUROSTD_MAX_UINT;
    m_name="";
    m_properties_node="";
    m_inodes.clear();
    m_vijobs.clear();
  }
}

namespace pflow {
  AQueue::AQueue(ostream& oss) : xStream(oss),m_initialized(false) {initialize();}
  AQueue::AQueue(ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize();}
  AQueue::AQueue(const aurostd::xoption& vpflow,ostream& oss) : xStream(oss),m_initialized(false) {initialize(vpflow);}
  AQueue::AQueue(const aurostd::xoption& vpflow,ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize(vpflow);}
  AQueue::AQueue(const AQueue& b) : xStream(*b.getOFStream(),*b.getOSS()) {copy(b);} // copy PUBLIC

  AQueue::~AQueue() {xStream::free();free();}

  const AQueue& AQueue::operator=(const AQueue& other) {
    if(this!=&other) {copy(other);}
    return *this;
  }

  void AQueue::clear() {free();}  //clear PUBLIC
  void AQueue::freeQueue() {
    m_qsys=QUEUE_SLURM;
    m_partitions.clear();
    m_nodes.clear();
    m_jobs.clear();
  }
  void AQueue::free() {
    m_initialized=false;
    m_flags.clear();
    freeQueue();
  }

  void AQueue::copy(const AQueue& b) {  //copy PRIVATE
    m_initialized=b.m_initialized;
    m_flags=b.m_flags;
    m_qsys=b.m_qsys;
    m_partitions.clear();for(uint i=0;i<b.m_partitions.size();i++){m_partitions.push_back(b.m_partitions[i]);}
    m_nodes.clear();for(uint i=0;i<b.m_nodes.size();i++){m_nodes.push_back(b.m_nodes[i]);}
    m_jobs.clear();for(uint i=0;i<b.m_jobs.size();i++){m_jobs.push_back(b.m_jobs[i]);}
  }

  bool AQueue::initialize(ostream& oss) {
    xStream::initialize(oss);
    return initialize();
  }
  bool AQueue::initialize(ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize();
  }
  bool AQueue::initialize(const aurostd::xoption& vpflow,ostream& oss) {
    xStream::initialize(oss);
    return initialize(vpflow);
  }
  bool AQueue::initialize(const aurostd::xoption& vpflow,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(vpflow);
  }

  bool AQueue::initialize() {
    free();
    return false;
  }
  bool AQueue::initialize(const aurostd::xoption& vpflow) {
    free();
    setFlags(vpflow);
    return true;
  }

  void AQueue::setFlags(const aurostd::xoption& vpflow){m_flags=vpflow;}

  uint AQueue::getNNodes() const {return m_nodes.size();}
  uint AQueue::getNCPUS() const {
    uint ncpus=0;
    for(uint i=0;i<m_nodes.size();i++){ncpus+=m_nodes[i].m_ncpus;}
    return ncpus;
  }
  uint AQueue::getNNodes(const APartition& partition) const {return partition.m_inodes.size();}
  uint AQueue::getNCPUS(const APartition& partition) const {
    uint ncpus=0;
    for(uint i=0;i<partition.m_inodes.size();i++){ncpus+=m_nodes[partition.m_inodes[i]].m_ncpus;}
    return ncpus;
  }
  uint AQueue::getNNodes(const APartition& partition,const node_status& status) const {
    uint nnodes=0;
    for(uint i=0;i<partition.m_inodes.size();i++){
      if(m_nodes[partition.m_inodes[i]].isStatus(status)){nnodes+=1;}
    }
    return nnodes;
  }
  uint AQueue::getNCPUS(const APartition& partition,const node_status& status_node,const cpus_status& status_cpus) const {
    uint ncpus=0;
    if(status_cpus==CPUS_TOTAL){
      for(uint i=0;i<partition.m_inodes.size();i++){
        if(m_nodes[partition.m_inodes[i]].isStatus(status_node)){
          ncpus+=m_nodes[partition.m_inodes[i]].m_ncpus;
        }
      }
    }
    else if(status_cpus==CPUS_FREE){
      for(uint i=0;i<partition.m_inodes.size();i++){
        if(m_nodes[partition.m_inodes[i]].isStatus(status_node)){
          ncpus+=m_nodes[partition.m_inodes[i]].m_ncpus-m_nodes[partition.m_inodes[i]].m_ncpus_occupied;
        }
      }
    }
    else if(status_cpus==CPUS_OCCUPIED){
      for(uint i=0;i<partition.m_inodes.size();i++){
        if(m_nodes[partition.m_inodes[i]].isStatus(status_node)){
          ncpus+=m_nodes[partition.m_inodes[i]].m_ncpus_occupied;
        }
      }
    }
    return ncpus;
  }
  uint AQueue::getNCPUS(const string& user,const string& partition,const job_status& status) const {
    uint ipartition=partitionName2Index(partition);
    if(ipartition>m_partitions.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"ipartition>m_partitions.size()-1",_INDEX_BOUNDS_);}
    return getNCPUS(user,m_partitions[ipartition],status);
  }
  uint AQueue::getNCPUS(const string& user,const APartition& partition,const job_status& status) const {
    bool LDEBUG=(FALSE || _AQUEUE_DEBUG_ || XHOST.DEBUG);
    uint ncpus=0;
    uint ijob=0;
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " partition=" << partition.m_name << endl;}
    for(uint i=0;i<partition.m_vijobs.size();i++){
      ijob=partition.m_vijobs[i];
      if(ijob>m_jobs.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"ijob>m_jobs.size()-1",_INDEX_BOUNDS_);}
      const AJob& job=m_jobs[ijob];
      if(job.m_user==user && job.m_status==status){
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " job=" << job.m_id << " ncpus=" << job.m_ncpus << endl;}
        ncpus+=job.m_ncpus;
      }
    }
    return ncpus;
  }
  double AQueue::getPercentage(const string& user,const string& partition,const job_status& status) const {
    uint ipartition=partitionName2Index(partition);
    if(ipartition>m_partitions.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"ipartition>m_partitions.size()-1",_INDEX_BOUNDS_);}
    return getPercentage(user,m_partitions[ipartition],status);
  }
  double AQueue::getPercentage(const string& user,const APartition& partition,const job_status& status) const {
    uint denom=getNCPUS(partition,NODE_OPERATIONAL);
    if(denom==0){return 0.0;}
    return (double)getNCPUS(user,partition,status)/(double)denom;
  }

  uint AQueue::nodeName2Index(const string& name) const {
    for(uint inode=0;inode<m_nodes.size();inode++){
      if(m_nodes[inode].m_name==name){
        return inode;
      }
    }
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No node found with name="+name,_INPUT_UNKNOWN_);
  }
  uint AQueue::partitionName2Index(const string& name) const {
    for(uint ipartition=0;ipartition<m_partitions.size();ipartition++){
      if(m_partitions[ipartition].m_name==name){
        return ipartition;
      }
    }
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No partition found with name="+name,_INPUT_UNKNOWN_);
  }

  uint getTORQUEIDFromString(const string& torqueid_str) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string tmp="";
    vector<string> tokens;
    aurostd::string2tokens(torqueid_str,tokens,".");
    if(tokens.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"tokens.size()==0 [1]",_RUNTIME_ERROR_);}
    tmp=tokens[0];
    aurostd::string2tokens(tmp,tokens,"/");
    if(tokens.size()==1){tmp=tokens[0];}  //no '/' found
    else if(tokens.size()==2){tmp=tokens[1];}
    else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"tokens.size()!=1||2",_RUNTIME_ERROR_);}
    if(!aurostd::isfloat(tmp)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"!aurostd::isfloat(tmp)",_RUNTIME_ERROR_);}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " id=" << tmp << endl;}
    return aurostd::string2utype<uint>(tmp);
  }

  bool AQueue::addJob(const AJob& _job){  //add job to m_jobs
    bool LDEBUG=(FALSE || _AQUEUE_DEBUG_ || XHOST.DEBUG);
    bool job_added=false;
    bool found_job=false,found_node=false;
    uint i=0,j=0;
    for(i=0;i<m_jobs.size()&&found_job==false;i++){
      if(m_jobs[i].m_id==_job.m_id){
        //job found
        AJob& job=m_jobs[i];
        if(job.m_status!=_job.m_status){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"job.m_status!=_job.m_status",_RUNTIME_ERROR_);}  //mismatch in statuses
        const job_status& jstat=job.m_status;
        if(jstat==JOB_RUNNING && _job.m_vinodes.size()==0){ //then where is it running? we should have node information (m_vinodes and m_vncpus)
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"RUNNING job="+aurostd::utype2string(_job.m_id)+" has not been assigned to a node",_RUNTIME_ERROR_);
        }
        if(_job.m_vinodes.size()>0){  //node information is available for mapping to a node
          //to add a job, its record should have been taken from one node (at a time)
          if(_job.m_vinodes.size()!=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"_job.m_vinodes.size()!=1",_RUNTIME_ERROR_);}
          if(_job.m_vncpus.size()!=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"_job.m_vncpus.size()!=1",_RUNTIME_ERROR_);}
          if(_job.m_vinodes.size()!=_job.m_vncpus.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"_job.m_vinodes.size()!=_job.m_vncpus.size()",_RUNTIME_ERROR_);}
          const uint inode=_job.m_vinodes[0];
          const uint ncpus=_job.m_vncpus[0];
          found_node=false;
          for(j=0;j<job.m_vinodes.size()&&found_node==false;j++){
            if(job.m_vinodes[j]==inode){
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " adding " << ncpus << " cpus to job=" << job.m_id << ",inode=" << job.m_vinodes[j] << endl;}
              job.m_ncpus+=ncpus;
              job.m_vncpus[j]+=ncpus;
              found_node=true;
            }
          }
          if(!found_node){
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " job=" << job.m_id << " is also running on inode=" << inode << " with " << ncpus << " cpus" << endl;}
            job.m_ncpus+=ncpus;
            job.m_vinodes.push_back(inode);
            job.m_vncpus.push_back(ncpus);
          }
        }
        found_job=true;
      }
    }
    if(!found_job){
      m_jobs.push_back(_job);m_jobs.back().m_index=m_jobs.size()-1;job_added=true;
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " adding NEW job=" << m_jobs.back().m_id << endl;}
    }
    return job_added;
  }

  bool AQueue::addPartition(const APartition& _partition){
    bool LDEBUG=(FALSE || _AQUEUE_DEBUG_ || XHOST.DEBUG);
    uint i=0;
    for(i=0;i<m_partitions.size();i++){
      if(m_partitions[i].m_name==_partition.m_name){return false;}
    }
    bool partition_added=false;
    m_partitions.push_back(_partition);m_partitions.back().m_index=m_partitions.size()-1;partition_added=true;
    APartition& partition=m_partitions.back();
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " adding partition=" << partition.m_name << endl;}
    if(partition.m_properties_node.empty()){  //only available if user is root
      partition.m_properties_node=partition.m_name; //default string to search
      //[resort to default]throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Not sure how to define partitions for HOST="+XHOST.hostname,_RUNTIME_ERROR_);
      if(XHOST.hostname=="qrats.materials.duke.edu"){ //protection, ideally the aflow user would be "root"
        //hack
        if(partition.m_name=="priority"){partition.m_properties_node="sharedCompute";}
        else if(partition.m_name=="debug"){partition.m_properties_node="debug";}
        else if(partition.m_name=="batch"){partition.m_properties_node="sharedCompute";}
        else if(partition.m_name=="research"){partition.m_properties_node="sharedCompute";}
        else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown queue="+partition.m_name,_RUNTIME_ERROR_);}
      }
    }
    return partition_added;
  }

#define DELIM_PROPERTIES_NODE ","

  bool AQueue::addNode(const ANode& _node){
    bool LDEBUG=(FALSE || _AQUEUE_DEBUG_ || XHOST.DEBUG);
    uint inode=0;
    for(inode=0;inode<m_nodes.size();inode++){
      if(m_nodes[inode].m_name==_node.m_name){
        ANode& node=m_nodes[inode];
        //the find()!=string::npos could be problematic if we have queues like 'fast' and 'faster'
        if(node.m_properties!=_node.m_properties&&node.m_properties.find(_node.m_properties)==string::npos){  //QUSER patch: a comma-separated string is important for quser where each node is listed N times for each partition it is associated with
          //from http://docs.adaptivecomputing.com/torque/4-2-8/Content/topics/4-serverPolicies/mappingQueueToRes.htm
          //"TORQUE does not currently provide a simple mechanism for mapping queues to nodes. However, schedulers such as Moab and Maui can provide this functionality."
          //since properties can be anything, we add additional matching queues as a comma-separated string
          //[sort]if(!node.m_properties.empty()){node.m_properties+=DELIM_PROPERTIES_NODE;}
          //[sort]if(!_node.m_properties.empty()){node.m_properties+=_node.m_properties;}
          vector<string> tokens;
          aurostd::string2tokens(node.m_properties,tokens,DELIM_PROPERTIES_NODE);
          tokens.push_back(_node.m_properties);
          std::sort(tokens.begin(),tokens.end());
          node.m_properties=aurostd::joinWDelimiter(tokens,DELIM_PROPERTIES_NODE);
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " adding properties=\"" << _node.m_properties << "\" to node=" << node.m_name << " properties (NEW properties=" << node.m_properties << ")" << endl;}
        }
        return false;
      }
    }
    bool node_added=false;
    m_nodes.push_back(_node);m_nodes.back().m_index=m_nodes.size()-1;node_added=true;
    ANode& node=m_nodes.back();
    //determine OCCUPIED vs. FULL
    if(node.m_status==NODE_FREE && node.m_ncpus_occupied>0){node.m_status=NODE_OCCUPIED;}
    if(node.m_status==NODE_OCCUPIED && node.m_ncpus_occupied>=node.m_ncpus){node.m_status=NODE_FULL;}
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " adding ";
      if(node.m_status==NODE_FREE){cerr << "FREE";}
      else if(node.m_status==NODE_OCCUPIED){cerr << "OCCUPIED";}
      else if(node.m_status==NODE_FULL){cerr << "FULL";}
      else if(node.m_status==NODE_DOWN){cerr << "DOWN";}
      else if(node.m_status==NODE_OFFLINE){cerr << "OFFLINE";}
      else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown node.status",_RUNTIME_ERROR_);}
      cerr << " node=" << node.m_name << endl;
    }
    return node_added;
  }

  void AQueue::nodePartitionMapping(ANode& node){
    bool LDEBUG=(FALSE || _AQUEUE_DEBUG_ || XHOST.DEBUG);
    //get partition
    bool found=false;
    node.m_vipartitions.clear();
    uint ipartition=0,iproperty=0;
    //parse node.m_properties by DELIM_PROPERTIES_NODE
    vector<string> vproperties_node;
    aurostd::string2tokens(node.m_properties,vproperties_node,DELIM_PROPERTIES_NODE);
    //search for node-partition matches
    for(ipartition=0;ipartition<m_partitions.size();ipartition++){
      APartition& partition=m_partitions[ipartition];
      for(iproperty=0;iproperty<vproperties_node.size();iproperty++){
        if(vproperties_node[iproperty]==partition.m_properties_node){
          partition.m_inodes.push_back(node.m_index);
          node.m_vipartitions.push_back(ipartition);
          found=true;
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " mapping node=" << node.m_name << " to partition=" << partition.m_name << endl;}
        }
      }
    }
    if(!found){
      APartition _partition;_partition.free();
      _partition.m_name=node.m_properties;
      //set m_properties_node inside addPartition()
      if(!addPartition(_partition)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Partition already exists",_RUNTIME_ERROR_);}
      ipartition=m_partitions.size()-1;
      APartition& partition=m_partitions[ipartition];
      partition.m_inodes.push_back(node.m_index);
      node.m_vipartitions.push_back(ipartition);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " mapping node=" << node.m_name << " to partition=" << partition.m_name << endl;}
    }
  }

  void AQueue::jobMapping(AJob& job){
    bool LDEBUG=(FALSE || _AQUEUE_DEBUG_ || XHOST.DEBUG);
    uint i=0,ipartition=0;
    //map job to partition
    //mapping MUST have happened upon adding the job, just checking and verbose
    if(job.m_vipartitions.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"job="+aurostd::utype2string(job.m_id)+" has not been assigned a partition",_INPUT_ERROR_);}
    for(i=0;i<job.m_vipartitions.size();i++){
      ipartition=job.m_vipartitions[i];
      if(ipartition>m_partitions.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"ipartition>m_partitions.size()-1",_INDEX_BOUNDS_);}
      APartition& partition=m_partitions[ipartition];
      partition.m_vijobs.push_back(job.m_index);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " mapping job=" << job.m_id << " to partition=" << partition.m_name << endl;}
    }
    //map RUNNING job to node
    if(job.m_status==JOB_RUNNING && job.m_vinodes.size()==0){ //then where is it running? we should have node information (m_vinodes and m_vncpus)
      //it's possible that, in the non-zero amount of time between processing pbsnodes and qstat, a job status changes from 'Q' or 'H' to 'R'
      //this "transition state' is always possible when the commands are not run at the exact same time
      //therefore, it is better to sleep for a few seconds, and re-reun the whole processing script again
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"RUNNING job="+aurostd::utype2string(job.m_id)+" has not been assigned a node",_RUNTIME_EXTERNAL_FAIL_); //this will issue a re-run
    }
    if(job.m_vinodes.size()>0){
      uint inode=0;
      for(i=0;i<job.m_vinodes.size();i++){
        inode=job.m_vinodes[i];
        if(inode>m_nodes.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"inode>m_nodes.size()-1",_INDEX_BOUNDS_);}
        ANode& node=m_nodes[inode];
        node.m_vijobs.push_back(job.m_index);
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " mapping job=" << job.m_id << " to node=" << node.m_name << endl;}
        //[WRONG - job maps to node and node to partitions, but job maps to 1 (usually) partition]if(node.m_vipartitions.size()==0){nodePartitionMapping(node);}  //safety
        //[WRONG - job maps to node and node to partitions, but job maps to 1 (usually) partition]//we could check that ALL nodes belong to a partition, but we might have a separate node...
        //[WRONG - job maps to node and node to partitions, but job maps to 1 (usually) partition]for(j=0;j<node.m_vipartitions.size();j++){
        //[WRONG - job maps to node and node to partitions, but job maps to 1 (usually) partition]  ipartition=node.m_vipartitions[j];
        //[WRONG - job maps to node and node to partitions, but job maps to 1 (usually) partition]  if(ipartition>m_partitions.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"ipartition>m_partitions.size()-1",_INDEX_BOUNDS_);}
        //[WRONG - job maps to node and node to partitions, but job maps to 1 (usually) partition]  APartition& partition=m_partitions[ipartition];
        //[WRONG - job maps to node and node to partitions, but job maps to 1 (usually) partition]  partition.m_vijobs.push_back(job.m_index);
        //[WRONG - job maps to node and node to partitions, but job maps to 1 (usually) partition]  job.m_vipartitions.push_back(ipartition);
        //[WRONG - job maps to node and node to partitions, but job maps to 1 (usually) partition]  if(LDEBUG){cerr << __AFLOW_FUNC__ << " mapping job=" << job.m_id << " to partition=" << partition.m_name << endl;}
        //[WRONG - job maps to node and node to partitions, but job maps to 1 (usually) partition]}
      }
    }
  }

#define ATTEMPTS_GETQUEUE_MAX 10
  void AQueue::getQueue() {
    bool LDEBUG=(TRUE || _AQUEUE_DEBUG_ || XHOST.DEBUG);

    bool SUCCESS=false;
    uint attempts=0;
    while(attempts++<ATTEMPTS_GETQUEUE_MAX && SUCCESS==false){
      try {freeQueue();processQueue();SUCCESS=true;}
      catch(aurostd::xerror& err) {
        SUCCESS=false;
        if(err.whatCode()!=_RUNTIME_EXTERNAL_FAIL_){attempts=ATTEMPTS_GETQUEUE_MAX;}
        if(attempts<ATTEMPTS_GETQUEUE_MAX){
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " failed attempt, trying to process queue again" << endl;}
          aurostd::Sleep(2);
        }
        else{
          if(SUCCESS==false){
            pflow::logger(err.whereFileName(),err.whereFunction(),err.what(),std::cout,_LOGGER_ERROR_);
          }
        }
      }
    }
  }

  void AQueue::readNodesPartitionsSLURM(){
    bool LDEBUG=(FALSE || _AQUEUE_DEBUG_ || XHOST.DEBUG);

    vector<string> lines,tokens,tokens2;
    uint iline=0;
    string line="",key="",value="",tmp="";

    //run sinfo
    //https://slurm.schedmd.com/sinfo.html
    if(!aurostd::IsCommandAvailable("sinfo")){
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"no sinfo command found for SLURM configuration",_RUNTIME_EXTERNAL_MISS_);
    }
    lines=aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("sinfo")+" -N -l"));
    APartition _partition;_partition.free();
    ANode _node;_node.free();
    for(iline=0;iline<lines.size();iline++){
      if(lines[iline].find("NODELIST")!=string::npos){continue;}  //skip date and header
      line=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(lines[iline]);
      if(line.empty()){continue;}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " line=\"" << line << "\"" << endl;}
      _node.free();
      aurostd::string2tokens(line,tokens," ");
      if(tokens.size()<11){continue;}   //skip date and header //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"tokens.size()<11",_RUNTIME_ERROR_);
      _node.m_name=tokens[0];
      //tokens[2] is partition, addPartition() takes care of checking that it's not already in the list
      _partition.free();
      _partition.m_name=tokens[2];aurostd::StringSubst(_partition.m_name,"*","");  //this signifies default queue, we don't care
      addPartition(_partition);
      _node.m_properties=_partition.m_name; //for node-partition mapping() later
      //tokens[3] is status
      if(tokens[3].find('*')!=string::npos){_node.m_status=NODE_DOWN;}  //catch first, doesn't matter what status it is in, IT'S NOT RESPONDING
      else if(tokens[3]=="allocated+"||tokens[3]=="allocated"||tokens[3]=="alloc"||
          tokens[3]=="completing"||tokens[3]=="comp"||
          FALSE){_node.m_status=NODE_FULL;}  //most likely to appear first, quicker to appear at the top
      else if(tokens[3]=="idle"){_node.m_status=NODE_FREE;}
      else if(tokens[3]=="mixed"||tokens[3]=="mix"){_node.m_status=NODE_OCCUPIED;}
      else if(tokens[3]=="reserved"||tokens[3]=="resv"||tokens[3]=="unknown"||tokens[3]=="unk"){_node.m_status=NODE_OFFLINE;}
      else if(tokens[3]=="down"||
          tokens[3]=="drained"||tokens[3]=="drain"||tokens[3]=="draining"||tokens[3]=="drng"||
          tokens[3]=="fail"||tokens[3]=="failing"||tokens[3]=="failg"||
          tokens[3]=="future"||tokens[3]=="futr"||
          tokens[3]=="maint"||tokens[3]=="reboot"||
          tokens[3]=="perfctrs"||tokens[3]=="npc"||
          tokens[3]=="power_down"||tokens[3]=="powering_down"||tokens[3]=="pow_dn"||
          tokens[3]=="power_up"||tokens[3]=="powering_up"||tokens[3]=="pow_up"||
          tokens[3]=="no_respond"||
          FALSE){_node.m_status=NODE_DOWN;} //leave for last, many variants to consider
      //tokens[4] is ncpus
      if(!aurostd::isfloat(tokens[4])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"tokens[4] is NOT a float (ncpus)",_RUNTIME_ERROR_);}
      _node.m_ncpus=aurostd::string2utype<double>(tokens[4]);
      addNode(_node);
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " found " << m_partitions.size() << " partitions" << endl;}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " found " << m_nodes.size() << " nodes" << endl;}
  }

  void AQueue::readJobsSLURM(){
    bool LDEBUG=(FALSE || _AQUEUE_DEBUG_ || XHOST.DEBUG);

    vector<string> lines,tokens,tokens2;
    uint iline=0,i=0;
    string line="",key="",value="",tmp="";

    //run squeue
    if(!aurostd::IsCommandAvailable("squeue")){
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"no squeue command found for SLURM configuration",_RUNTIME_EXTERNAL_MISS_);
    }
    lines=aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("squeue")+" --format=\"\%.18i %.9P %.8j \%.8u %.2t %.10M \%.6D \%C \%R\"")); //man squeue
    AJob _job;_job.free();
    for(iline=0;iline<lines.size();iline++){
      if(lines[iline].find("JOBID")!=string::npos){continue;}  //skip date and header
      line=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(lines[iline]);
      if(line.empty()){continue;}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " line=\"" << line << "\"" << endl;}
      aurostd::string2tokens(line,tokens," ");
      //tokens[0] is job id
      if(!aurostd::isfloat(tokens[0])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"!aurostd::isfloat(tokens[0])",_RUNTIME_ERROR_);}
      _job.free();
      _job.m_id=aurostd::string2utype<uint>(tokens[0]);
      //tokens[1] is partition
      aurostd::string2tokens(tokens[1],tokens2,",");
      for(i=0;i<tokens2.size();i++){
        _job.m_vipartitions.push_back(partitionName2Index(tokens2[i]));
      }
      //tokens[4] is status
      //https://curc.readthedocs.io/en/latest/running-jobs/squeue-status-codes.html
      if(tokens[4]=="PD"){_job.m_status=JOB_QUEUED;}
      else if(tokens[4]=="R"||tokens[4]=="ST"){_job.m_status=JOB_RUNNING;}  //ST because it retains cores
      else if(tokens[4]=="S"){_job.m_status=JOB_HELD;} //S because cores are given up for other jobs (transition state)
      else if(tokens[4]=="CD"||tokens[4]=="CG"||
          tokens[4]=="F"||
          tokens[4]=="PR"||
          FALSE){_job.m_status=JOB_DONE;}
      else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown job status="+tokens[4],_RUNTIME_ERROR_);}
      //tokens[7] is ncpus
      if(!aurostd::isfloat(tokens[7])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"!aurostd::isfloat(tokens[7])",_RUNTIME_ERROR_);}
      _job.m_ncpus=aurostd::string2utype<uint>(tokens[7]);
      //tokens[8] is node if running
      if(_job.m_status==JOB_RUNNING){
        _job.m_vinodes.push_back(nodeName2Index(tokens[8]));
        _job.m_vncpus.push_back(_job.m_ncpus);
      }
      addJob(_job);
    }
  }

  void AQueue::readPartitionsTORQUE(){
    bool LDEBUG=(FALSE || _AQUEUE_DEBUG_ || XHOST.DEBUG);

    vector<string> lines,tokens,tokens2;
    uint iline=0;
    string line="",key="",value="",tmp="";

    //run qstat -f -Q
    if(!aurostd::IsCommandAvailable("qstat")){
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"no qstat command found for TORQUE configuration",_RUNTIME_EXTERNAL_MISS_);
    }
    //'qstat -f' gives FULL job information, but takes a long time
    lines=aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("qstat")+" -f -Q"));
    APartition _partition;_partition.free();
    _partition.m_name="";
    _partition.m_properties_node="";
    for(iline=0;iline<lines.size();iline++){
      line=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(lines[iline]);
      if(line.empty()){continue;}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " line=\"" << line << "\"" << endl;}
      if(line.find("Queue:")!=string::npos){
        if(!_partition.m_name.empty()){addPartition(_partition);}
        aurostd::string2tokens(line,tokens,":");
        if(tokens.size()!=2){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"tokens.size()!=2",_RUNTIME_ERROR_);}
        _partition.free(); //reset
        _partition.m_name=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[1]);
      }
      else{
        aurostd::string2tokens(line,tokens,"=");
        if(tokens.size()!=2){continue;}   //sometimes we get garbage lines like 'lete:0' //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"!tokens.size()",_RUNTIME_ERROR_);
        key=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[0]);
        value=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[1]);
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " key=" << key << endl;}
        if(key=="resources_default.neednodes"){_partition.m_properties_node=value;}  //really only works for root, but we leave here for now
      }
    }
    if(!_partition.m_name.empty()){addPartition(_partition);}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " found " << m_partitions.size() << " partitions" << endl;}
  }

  void AQueue::readNodesJobsTORQUE(){
    bool LDEBUG=(FALSE || _AQUEUE_DEBUG_ || XHOST.DEBUG);

    vector<string> lines,tokens,tokens2;
    uint iline=0,i=0;
    string line="",key="",value="",tmp="";

    //run pbsnodes
    if(!aurostd::IsCommandAvailable("pbsnodes")){
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"no pbsnodes command found for TORQUE configuration",_RUNTIME_EXTERNAL_MISS_);
    }
    lines=aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("pbsnodes")));
    ANode _node;_node.free();
    AJob _job;_job.free();
    for(iline=0;iline<lines.size();iline++){
      line=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(lines[iline]);
      if(line.empty()){continue;}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " line=\"" << line << "\"" << endl;}
      if(line.find('=')==string::npos){
        if(!_node.m_properties.empty()){addNode(_node);} //we have node information loaded up
        _node.free();  //reset node
        _node.m_name=line;
        continue;
      }
      else{
        aurostd::string2tokens(line,tokens,"=");
        if(tokens.size()!=2){continue;}   //sometimes we get garbage lines like 'lete:0' //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"!tokens.size()",_RUNTIME_ERROR_);
        key=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[0]);
        value=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[1]);
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " key=" << key << endl;}
        if(key=="state"){
          //possible states: http://docs.adaptivecomputing.com/torque/4-0-2/Content/topics/commands/pbsnodes.htm
          //"active", "all", "busy", "down", "free", "job-exclusive", "job-sharing", "offline", "reserve", "state-unknown", "time-shared", and "up"
          if(value.find("down")!=string::npos){_node.m_status=NODE_DOWN;} //down first so we can rectify
          else if(value.find("offline")!=string::npos){_node.m_status=NODE_OFFLINE;} //offline next
          else if(value.find("job-exclusive")!=string::npos){_node.m_status=NODE_OCCUPIED;}
          else if(value.find("free")!=string::npos){_node.m_status=NODE_FREE;}
          else{
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"unknown state from pbsnodes: state=\""+value+"\"",_RUNTIME_ERROR_);  //more states to account for
          }
        }
        else if(key=="np"){_node.m_ncpus=aurostd::string2utype<uint>(value);}
        else if(key=="properties"){_node.m_properties=value;}  //needed to match with queue
        else if(key=="jobs"){
          aurostd::string2tokens(value,tokens,",");
          _node.m_ncpus_occupied=tokens.size();  //easy, this is all we need for now
          //parse 0/709843.qrats.materials.duke.edu
          for(i=0;i<tokens.size();i++){
            _job.free();
            _job.m_id=getTORQUEIDFromString(tokens[i]);
            _job.m_status=JOB_RUNNING; //if it's here, it's running
            _job.m_ncpus=1;
            _job.m_vinodes.push_back(m_nodes.size());  //we haven't added the node yet, so this will be 0 first...
            _job.m_vncpus.push_back(_job.m_ncpus);
            addJob(_job);
          }
        }
        //no else, there are other keys that we don't process
      }
    }
    if(!_node.m_properties.empty()){addNode(_node);} //we have node information loaded up
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " found " << m_nodes.size() << " nodes" << endl;}
  }

  void AQueue::readJobsTORQUE(){
    bool LDEBUG=(FALSE || _AQUEUE_DEBUG_ || XHOST.DEBUG);

    vector<string> lines,tokens,tokens2;
    uint iline=0,i=0;
    string line="",key="",value="",tmp="";
    uint id=AUROSTD_MAX_UINT;
    bool found=false;
    AJob _job;_job.free();

    //get jobs user
    if(!aurostd::IsCommandAvailable("qstat")){
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"no qstat command found for TORQUE configuration",_RUNTIME_EXTERNAL_MISS_);
    }
    //'qstat -f' gives FULL job information, but takes a long time
    lines=aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("qstat")+" -a"));
    for(iline=0;iline<lines.size();iline++){
      if(lines[iline].find("Job ID")!=string::npos){continue;}  //skip header lines
      line=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(lines[iline]);
      if(line.empty()){continue;}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " line=\"" << line << "\"" << endl;}
      aurostd::string2tokens(line,tokens," ");
      if(tokens.size()<11){continue;}   //skip date and header //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"tokens.size()<11",_RUNTIME_ERROR_);
      //tokens[0] has job id
      tmp=tokens[0];
      if(tmp.find("---")!=string::npos){continue;}  //skip header //must be longer than 2 --, see 'Req'd Memory'
      //709759.qrats.materials
      id=getTORQUEIDFromString(tmp);
      found=false;
      for(i=0;i<m_jobs.size()&&found==false;i++){
        if(m_jobs[i].m_id==id){
          //get user: tokens[1]
          m_jobs[i].m_user=tokens[1];
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " adding user=" << m_jobs[i].m_user << " to job=" << m_jobs[i].m_id << endl;}
          //get partition: tokens[2]
          m_jobs[i].m_vipartitions.push_back(partitionName2Index(tokens[2]));
          found=true;
        }
      }
      if(found==false){ //add job
        _job.free();
        _job.m_id=id;
        _job.m_user=tokens[1];
        //tokens[2] has queue, very important for non-running jobs
        _job.m_vipartitions.push_back(partitionName2Index(tokens[2]));
        //tokens[6] has ncpus
        if(!aurostd::isfloat(tokens[6])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"!aurostd::isfloat(tokens[6])",_RUNTIME_ERROR_);}
        _job.m_ncpus=aurostd::string2utype<uint>(tokens[6]);
        //tokens[9] has job status
        //http://docs.adaptivecomputing.com/torque/4-1-4/Content/topics/commands/qstat.htm
        if(tokens[9]=="Q"){_job.m_status=JOB_QUEUED;}
        else if(tokens[9]=="R"){_job.m_status=JOB_RUNNING;}
        else if(tokens[9]=="H"||tokens[9]=="S"){_job.m_status=JOB_HELD;}
        else if(tokens[9]=="C"||tokens[9]=="E"||tokens[9]=="W"){_job.m_status=JOB_DONE;} //not sure about W...
        else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown job status="+tokens[9],_RUNTIME_ERROR_);} //T?
        addJob(_job);
      }
    }
  }

  void AQueue::processQueue() {
    bool LDEBUG=(FALSE || _AQUEUE_DEBUG_ || XHOST.DEBUG);

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    if(aurostd::IsCommandAvailable("squeue")){m_qsys=QUEUE_SLURM;}
    else if(aurostd::IsCommandAvailable("pbsnodes")){m_qsys=QUEUE_TORQUE;}
    else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"cannot identify queuing system",_RUNTIME_ERROR_);}

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " m_qsys=";
      if(m_qsys==QUEUE_SLURM){cerr << "SLURM";}
      if(m_qsys==QUEUE_TORQUE){cerr << "TORQUE";}
      cerr << endl;
    }

    if(m_qsys==QUEUE_SLURM){
      readNodesPartitionsSLURM();
      readJobsSLURM();
    }
    else if(m_qsys==QUEUE_TORQUE){
      readPartitionsTORQUE();
      readNodesJobsTORQUE();  //get job numbers and association to nodes here
      readJobsTORQUE();  //get job users here, validate information from readNodesJobsTORQUE()

      //may need showq later for node-specific submission
      //if(LDEBUG){cerr << __AFLOW_FUNC__ << " XHOST.command(\"showq\")=" << XHOST.command("showq") << endl;}
    }

    //get mappings now - clear first START
    for(uint inode=0;inode<m_nodes.size();inode++){
      ANode& node=m_nodes[inode];
      node.m_vijobs.clear();
      node.m_vipartitions.clear();
    }
    for(uint ipartition=0;ipartition<m_partitions.size();ipartition++){
      APartition& partition=m_partitions[ipartition];
      partition.m_inodes.clear();
      partition.m_vijobs.clear();
    }
    //get mappings now - clear first END

    //get node mappings
    for(uint inode=0;inode<m_nodes.size();inode++){
      ANode& node=m_nodes[inode];
      nodePartitionMapping(node);
    }
    //get job mappings
    for(uint ijob=0;ijob<m_jobs.size();ijob++){
      AJob& job=m_jobs[ijob];
      jobMapping(job);
    }

    if(LDEBUG){
      //see node to job mapping
      uint ijob=0;
      for(uint inode=0;inode<m_nodes.size();inode++){
        ANode& node=m_nodes[inode];
        cerr << __AFLOW_FUNC__ << " node=" << node.m_name << " has jobs=";
        for(uint i=0;i<node.m_vijobs.size();i++){
          ijob=node.m_vijobs[i];
          if(ijob>m_jobs.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"ijob>m_jobs.size()-1",_INDEX_BOUNDS_);}
          AJob& job=m_jobs[ijob];
          cerr << job.m_id << (i<node.m_vijobs.size()-1?",":"");
        }
        cerr << endl;
      }
    }

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " " << XHOST.hostname << " has " << getNNodes() << " nodes (" << getNCPUS() << " cpus)" << endl;
      bool add_comma=false;
      for(uint i=0;i<m_partitions.size();i++){
        add_comma=false;
        cerr << __AFLOW_FUNC__ << " " << m_partitions[i].m_name << ":";
        if(getNCPUS(m_partitions[i],NODE_FREE)>0){cerr << (add_comma?",":"") << " " << getNNodes(m_partitions[i],NODE_FREE) << " free nodes (" << getNCPUS(m_partitions[i],NODE_FREE) << " cpus)";if(add_comma==false){add_comma=true;}}
        if(getNCPUS(m_partitions[i],NODE_OCCUPIED)>0){cerr << (add_comma?",":"") << " " << getNNodes(m_partitions[i],NODE_OCCUPIED) << " occupied nodes (" << getNCPUS(m_partitions[i],NODE_OCCUPIED,CPUS_FREE) << " cpus free, " << getNCPUS(m_partitions[i],NODE_OCCUPIED,CPUS_OCCUPIED) << " cpus occupied)";if(add_comma==false){add_comma=true;}}
        if(getNCPUS(m_partitions[i],NODE_FULL)>0){cerr << (add_comma?",":"") << " " << getNNodes(m_partitions[i],NODE_FULL) << " full nodes (" << getNCPUS(m_partitions[i],NODE_FULL) << " cpus)";if(add_comma==false){add_comma=true;}}
        if(getNCPUS(m_partitions[i],NODE_DOWN)>0){cerr << (add_comma?",":"") << " " << getNNodes(m_partitions[i],NODE_DOWN) << " down nodes (" << getNCPUS(m_partitions[i],NODE_DOWN) << " cpus)";if(add_comma==false){add_comma=true;}}
        if(getNCPUS(m_partitions[i],NODE_OFFLINE)>0){cerr << (add_comma?",":"") << " " << getNNodes(m_partitions[i],NODE_OFFLINE) << " offline nodes (" << getNCPUS(m_partitions[i],NODE_OFFLINE) << " cpus)";if(add_comma==false){add_comma=true;}}
        if(getNCPUS(m_partitions[i])>0){cerr << (add_comma?",":"") << " " << getNNodes(m_partitions[i]) << " TOTAL nodes (" << getNCPUS(m_partitions[i]) << " cpus)";if(add_comma==false){add_comma=true;}}
        cerr << endl;
      }
    }

  }
}

//CO20200526 - queueing class
namespace pflow {
  string getQueueStatus(const aurostd::xoption& vpflow){  //CO20200526
    AQueue aqueue(vpflow);
    aqueue.getQueue();
    //[only test on qrats]cerr << aqueue.getPercentage("aflow","batch",JOB_RUNNING) << endl;
    return "";
  }
}

// ***************************************************************************
// pflow::getFakeElements()
// ***************************************************************************
namespace pflow{
  vector<string> getFakeElements(uint nspecies){

    // Return vector of fake elements
    // Useful for determining "elements" for prototypes


    if(nspecies>26){
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"There are more than 26 species, this function must be modified to include more fake elements.",_RUNTIME_ERROR_);
    }

    vector<string> elements;
    for(uint i=0;i<nspecies;i++){
      stringstream ss_letter; ss_letter << char('A'+i); // cannot type cast char to string directly //DX20200907 - use ASCII
      elements.push_back(ss_letter.str());
    }

    return elements;
  }
}

// ***************************************************************************
// pflow::hasRealElements()
// ***************************************************************************
namespace pflow{
  bool hasRealElements(const xstructure& xstr){

    // Determine if elements in the xstructure are real/physical.
    // Uses xelement

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;

    if(xstr.species.size() > 0){
      xelement::xelement element;
      uint nspecies = xstr.species.size();
      for(uint i=0;i<nspecies;i++){
        try{
          element = xelement::xelement(KBIN::VASP_PseudoPotential_CleanName(xstr.species[i]));
        }
        catch(aurostd::xerror& re){
          if(LDEBUG){ cerr << __AFLOW_FUNC__ << xstr.species[i] << " is not a real element." << endl; }
          return false;
        }
      }
    }
    else{
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"The species are empty.",_INPUT_NUMBER_);
    }

    return true;
  }
}

// ***************************************************************************
// pflow::getSymmetryTolerance() //DX20200820
// ***************************************************************************
namespace pflow{
  double getSymmetryTolerance(const xstructure& xstr, const string& tolerance_string){

    // Return the symmetry tolerance
    // options:
    //  1) tight = min_nn_dist/100
    //  2) loose = min_nn_dist/10
    //  3) number = user defined (Angstroms)

    stringstream message;

    double default_tolerance=SYM::defaultTolerance(xstr);
    double tolerance = AUROSTD_NAN;
    if(!tolerance_string.empty()){
      if(aurostd::toupper(tolerance_string[0]) == 'T'){ //Tight
        tolerance=default_tolerance;
      }
      else if(aurostd::toupper(tolerance_string[0]) == 'L'){ //Loose
        tolerance=default_tolerance*10.0;
      }
      else {
        tolerance=aurostd::string2utype<double>(tolerance_string);
      }
    }
    else {
      tolerance = default_tolerance;
    }
    if(tolerance < 1e-10){
      message << "Tolerance cannot be zero (i.e. less than 1e-10): tol=" << tolerance << ".";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_RANGE_);
    }

    return tolerance;
  }
}

// ***************************************************************************
// pflow::getSymmetryToleranceSpectrum() //DX20200820
// ***************************************************************************
namespace pflow{
  vector<double> getSymmetryToleranceSpectrum(const string& tolerance_range_string){

    // Return the symmetry tolerance spectrum
    // Expected input: "start:end:step"

    stringstream message;

    vector<double> tolerance_spectrum;

    vector<string> tokens;
    if(aurostd::string2tokens(tolerance_range_string,tokens,":") == 3){
      double start = aurostd::string2utype<double>(tokens[0]);
      double end = aurostd::string2utype<double>(tokens[1]);
      uint nsteps = aurostd::string2utype<uint>(tokens[2])-1;
      if(end<start){
        message << "End of the range cannot be less than the beginning of the range: start=" << start << ", end=" << end;
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
      }
      double interval = (end-start)/(double)nsteps;
      for(uint i=0;i<=nsteps;i++){ tolerance_spectrum.push_back(start+((double)i*interval)); }
    }
    else{
      message << "Expected three inputs: first=range_start:second=range_end:third=nsteps.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
    }

    return tolerance_spectrum;
  }
}



// ***************************************************************************
// pflow::writeAtomicEnvironment() //HE20210617
// ***************************************************************************

namespace pflow {

  /// @brief write collection of atomic environments into json files
  /// @param auid AFLOW ID
  /// @param aeMode enviroment definition (see ATOM_ENVIRONMENT_MODE_X in aflow.h)
  /// @param radius change the search radius [FUTURE]
  ///
  /// Feature request: analyse a structure file directly
  void outputAtomicEnvironment(const string &auid, uint aeMode, double radius) {
    bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) cerr << __AFLOW_FUNC__ << " Start" << endl;

    // for FUTURE use
    if (radius == 0){}

    string aurl = "";
    aflowlib::_aflowlib_entry entry;
    xstructure str;
    std::map<string, string> meta_data;

    // Quickly match auid to aurl (speedup from 15s to 0.5s compared to single use of AflowlibLocator)
    if (!auid.empty())
      aurl = aurostd::execute2string(XHOST.command("aflow_data") + " vLIBS | grep -B1 \"" + auid + "\" | head -n 1");
    // Fallback if quick search failed
    if (!auid.empty() && aurl.empty()) {
      cerr << __AFLOW_FUNC__ << " Quick search failed! Trying standard methode." << endl;
      aflowlib::AflowlibLocator(auid, aurl, "AFLOWLIB_AUID2AURL");
    }

    if (!aurl.empty()) {
      entry.aurl = aurl;
      entry.auid = auid;
      loadXstructures(entry);
      meta_data["aurl"] = aurl;
      meta_data["auid"] = auid;
      str = entry.vstr.back();
      if (LDEBUG) cerr << __AFLOW_FUNC__ << " AUID: " << auid << endl;
      if (LDEBUG) cerr << __AFLOW_FUNC__ << " AURL: " << aurl << endl;
    } else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "could not load a structure", _INPUT_ERROR_);
    }

    vector <AtomEnvironment> AE = getAtomEnvironments(str, aeMode);
    for(uint i=0; i<AE.size(); i++) AE[i].constructAtomEnvironmentHull();
    writeAtomEnvironments(AE, meta_data);
  }
}

// ***************************************************************************
// pflow::getSpaceGroupSetting() //DX20210420
// ***************************************************************************
namespace pflow{
  uint getSpaceGroupSetting(const string& setting_string, uint mode_default){

    // Return the space group setting
    // options:
    //  1) 1 (SG_SETTING_1)    = rhombohedral: rhl setting, monoclinic: unique-axis b, centrosymmetric: origin on high-symmetry site
    //  2) 2 (SG_SETTING_2)    = rhombohedral: hex setting, monoclinic: unique-axis c, centrosymmetric: origin on inversion site
    //  3) 3 (SG_SETTING_ANRL) = rhomobhedral: rhl setting, monoclinic: unique-axis b, centrosymmetric: origin on inversion site 
    // The mode_default variable is an optional input: --prototype command defaults to the AFLOW setting, while everything else defaults to SG_SETTING_1

    stringstream message;

    uint setting = mode_default; //default

    // space group setting
    if(!setting_string.empty()){
      if(aurostd::tolower(setting_string) == "aflow" || aurostd::tolower(setting_string) == "anrl"){
        setting=SG_SETTING_ANRL;
      }
      else {
        int setting_num=aurostd::string2utype<int>(setting_string);
        if(setting_num==1){setting=SG_SETTING_1;}
        else if(setting_num==2){setting=SG_SETTING_2;}
        if(setting_num!=SG_SETTING_1 && setting_num!=SG_SETTING_2){
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Setting must be 1, 2, or \"aflow\" (for rhombohedral systems: 1=rhl setting and 2=hex setting; for monoclinic systems: 1=unique axis-b and 2=unique axis-c).",_INPUT_ILLEGAL_);
        }
      }
    }

    return setting;
  }
}

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
#endif
