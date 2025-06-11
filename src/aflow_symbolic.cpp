// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *           Aflow DAVID HICKS - Duke University 2014-2021                 *
// *                                                                         *
// ***************************************************************************
// Written by David Hicks (DX) - 2020

#ifndef _AFLOW_SYMBOLIC_CPP
#define _AFLOW_SYMBOLIC_CPP
#include "aflow.h"
#include "aflow_symmetry_spacegroup.h"
#include "aflow_compare_structure.h"
#include "aflow_pflow.h"
#include "aflow_symbolic.h"

#define _DEBUG_SYMBOLIC_ false //DX20200625


#if USE_SYMBOLIC_SOURCE //DX20200831 - defined in aflow.h
// *************************************************************************** 
// *************************************************************************** 
// Symbolic math functions leveraging the SymbolicC++ source code
// 
// SymbolicC++ is written by Y. Hardy, W.-H. Steeb, and T.K. Shi 
// SymbolicC++ is open-source via the GNU-GPL licensce (version 2)
// see http://symboliccpp.sourceforge.net/ for more information
//
// AFLOW implementation:
// The SymbolicC++ files are included with AFLOW 
// (compiled locally, see SYMBOLICCPLUSPLUS/ directory).
// Note: The source has been modified to allow inclusion of the 
// functions into multiple places in AFLOW (function defintions moved from 
// header files to CPP files; to overcome the "One Definition Rule")
// Below are generic helper functions, e.g., convert string to Symbolic 
// representation and ANRL prototype specific functions
//
// Generic functions are in the "symbolic" namespace, and 
// ANRL specific functions are in the "anrl" namespace
//
// For more information on the AFLOW implementation, contact 
// David Hicks (DX) david.hicks@duke.edu
// *************************************************************************** 
// *************************************************************************** 

// ---------------------------------------------------------------------------
// GENERIC SYMBOLIC MATH FUNCTIONS

// ***************************************************************************
// template <typename T> char const* get_type()
// ***************************************************************************
template <typename T> char const* get_type(T const& object) {

  // get the type of the object
  // this function is needed to avoid compiler warning so that
  // the template type can be recognized at compile time

  return typeid(*object).name();
}

// *************************************************************************** 
// SymbolicWyckoffSite initializeSymbolicWyckoffSite() 
// *************************************************************************** 
SymbolicWyckoffSite initializeSymbolicWyckoffSite(const wyckoffsite_ITC& Wyckoff){

  // initialize symbolic Wyckoff site struct 

  SymbolicWyckoffSite Wyckoff_symbolic;
  Wyckoff_symbolic.coord = Wyckoff.coord;
  Wyckoff_symbolic.index = Wyckoff.index;
  Wyckoff_symbolic.type = Wyckoff.type;
  Wyckoff_symbolic.wyckoffSymbol = Wyckoff.wyckoffSymbol;
  Wyckoff_symbolic.letter = Wyckoff.letter;
  Wyckoff_symbolic.site_symmetry = Wyckoff.site_symmetry;
  Wyckoff_symbolic.multiplicity = Wyckoff.multiplicity;
  Wyckoff_symbolic.site_occupation = Wyckoff.site_occupation;
  Wyckoff_symbolic.equations = anrl::equations2SymbolicEquations(Wyckoff.equations);
  Wyckoff_symbolic.parameter_index = Wyckoff.parameter_index;
  return Wyckoff_symbolic;
}

// *************************************************************************** 
// substituteVariableWithParameterDesignation() 
// *************************************************************************** 
// vector<SymbolicWyckoffSite>
void substituteVariableWithParameterDesignation(vector<SymbolicWyckoffSite>& Wyckoff_sites_symbolic){
  for(uint i=0;i<Wyckoff_sites_symbolic.size();i++){
    substituteVariableWithParameterDesignation(Wyckoff_sites_symbolic[i]);
  }
}

// SymbolicWyckoffSite
void substituteVariableWithParameterDesignation(SymbolicWyckoffSite& Wyckoff_symbolic){

  // convert a generic variable to the parameter designation, e.g., x -> x2

  for(uint i=0;i<Wyckoff_symbolic.equations.size();i++){
    Wyckoff_symbolic.equations[i]=Wyckoff_symbolic.equations[i].subst("x","x"+aurostd::utype2string<uint>(Wyckoff_symbolic.parameter_index));
    Wyckoff_symbolic.equations[i]=Wyckoff_symbolic.equations[i].subst("y","y"+aurostd::utype2string<uint>(Wyckoff_symbolic.parameter_index));
    Wyckoff_symbolic.equations[i]=Wyckoff_symbolic.equations[i].subst("z","z"+aurostd::utype2string<uint>(Wyckoff_symbolic.parameter_index));
  }
}

// *************************************************************************** 
// symbolic::string2symbolic
// *************************************************************************** 
namespace symbolic {
  Symbolic string2symbolic(const string& str){

    // Convert a string into symbolic math notation
    // The SYMBOLICC++ library cannot convert a string into a symbol so
    // we must do it ourselves
    // NOTE: this function is rudimentary and NOT generalized, 
    // if this is to be used beyond Wyckoff positions, we need to add more 
    // operators/functions (e.g., sin, cos, exponentials, etc.)

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_SYMBOLIC_);
    stringstream message;

    // ---------------------------------------------------------------------------
    // break up string into doubles and characters
    vector<SYM::sdouble> sdouble_temp = SYM::simplify(str);

    // ---------------------------------------------------------------------------
    // loop over characters and build symbolic expression 
    Symbolic out;
    for(uint i=0;i<sdouble_temp.size();i++){
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " sdouble_temp: (dbl) " << sdouble_temp[i].dbl 
          << " (chr) " << sdouble_temp[i].chr << endl;
      }
      // ---------------------------------------------------------------------------
      // if contains char 
      if(sdouble_temp[i].chr != '\0'){
        stringstream ss_char; ss_char << sdouble_temp[i].chr;
        out += Symbolic(sdouble_temp[i].dbl)*Symbolic(ss_char.str()); // '*' is a big assumption, may need to generalize in the future
      }
      // ---------------------------------------------------------------------------
      // else contains only double
      else{
        out += Symbolic(sdouble_temp[i].dbl);
      }
    }

    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " symbolic expression: out=" << out << endl; }
    return out;
  }
}

// *************************************************************************** 
// symbolic::isEqualSymbolic()
// *************************************************************************** 
namespace symbolic {
  bool isEqual(const Symbolic& a, const Symbolic& b){

    // relies on SymbolicC++'s implementation of "=="
    // this is sufficient for now, e.g., (1e-9)*x+(1e-6)*y+0 = 0

    bool VERBOSE=FALSE; // VERBOSE INSTEAD OF LDEBUG SINCE A NESTED FUNCTION

    Symbolic diff = a - b;

    if(VERBOSE){ 
      cerr << __AFLOW_FUNC__ << " a-b=" << diff << endl;
    }

    return (diff==_SYMBOLIC_ZERO_);
  }
}

// *************************************************************************** 
// symbolic::isEqualSymbolicVector()
// *************************************************************************** 
namespace symbolic {
  bool isEqualVector(const Symbolic& a_vec, const Symbolic& b_vec){

    // check if Symbolic vector inputs are equal

    // ---------------------------------------------------------------------------
    // check that the inputs' types are SymbolicMatrix (vector)
    const char* a_vec_info = get_type(a_vec); //DX20200901
    const char* b_vec_info = get_type(b_vec); //DX20200901
    if(a_vec_info != typeid(SymbolicMatrix).name() ||
        b_vec_info != typeid(SymbolicMatrix).name()){
      stringstream message;
      message << "One or both of the inputs are not a SymbolicMatrix (i.e., typeids are different):"
        << " a_vec (input) id: " << a_vec_info
        << " b_vec (input) id: " << b_vec_info
        << " SymbolicMatrix id: " << typeid(SymbolicMatrix).name();
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
    }

    bool VERBOSE=FALSE; // VERBOSE INSTEAD OF LDEBUG SINCE FUNCITON IS NESTED

    if(VERBOSE){ 
      cerr << __AFLOW_FUNC__ << " a_vec-b_vec=" << (a_vec-b_vec) << endl;
    }

    for(uint i=0;i<3;i++){
      if(!isEqual(a_vec(i),b_vec(i))){ return false; }
    }

    return true;
  }
}

// *************************************************************************** 
// symbolic::matrix2VectorVectorString() 
// *************************************************************************** 
namespace symbolic {
  vector<vector<string> > matrix2VectorVectorString(const Symbolic& lattice){

    // convert symbolic matrix to vector<vector<string> >
    // perhaps check if type is matrix?

    // ---------------------------------------------------------------------------
    // check that input type is a SymbolicMatrix 
    const char* lattice_info = get_type(lattice); //DX20200901
    if(lattice_info != typeid(SymbolicMatrix).name()){
      stringstream message;
      message << "The input is not a SymbolicMatrix (i.e., typeids are different): lattice (input) id: "
        << lattice_info << " SymbolicMatrix id: " << typeid(SymbolicMatrix).name();
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
    }

    vector<vector<string> > vvstring;

    for(uint i=0;i<3;i++){
      vector<string> row;
      for(uint j=0;j<3;j++){
        stringstream ss_component; ss_component << lattice.row(i)(j);
        row.push_back(ss_component.str());
      }
      vvstring.push_back(row);
    }

    return vvstring;
  }
}

// *************************************************************************** 
// symbolic::BringInCell() 
// *************************************************************************** 
namespace symbolic {
  Symbolic BringInCell(const Symbolic& vec_in, double tolerance, double upper_bound, double lower_bound){

    // bring symbolic math inside the cell
    // it is impossible to know what the variable (x, y, or z) will be 
    // to truly bring the coordinate in the cell, but this function brings 
    // any constant term(s) in cell
    // here, we use the Symbolic.coeff(<variable>,<order>) function, where
    // <variable>: is the variable of interest
    // <order>: is the order/degree of the variable
    // TRICK: define a constant and provide order 1, i.e.,
    // Symbolic constant = Symbolic(1); equation.coeff(constant,1);
    // tune the upper and lower bound (e.g, 0 to 1 or -0.5 to 0.5)
    // NOTE: preference to lower bound (e.g., 0) vs upper bound (e.g., 1)

    Symbolic vec_out = vec_in;

    // ---------------------------------------------------------------------------
    // define a constant so we can determine the "coefficients of the constant"
    Symbolic constant = Symbolic(1); // 1 is a great constant

    for(uint i=0;i<3;i++){
      // ---------------------------------------------------------------------------
      // bring in cell based on constant term
      while(double(vec_out(i).coeff(constant,1))-upper_bound>=-tolerance){ vec_out(i) -= 1; }
      while(double(vec_out(i).coeff(constant,1)-lower_bound)<-tolerance){ vec_out(i) += 1; }
    }
    return vec_out;
  }
}


// ---------------------------------------------------------------------------
// ANRL SPECIFIC SYMBOLIC MATH FUNCTIONS

// *************************************************************************** 
// anrl::SymbolicANRLPrimitiveLattices()
// *************************************************************************** 
namespace anrl {
  symbolic::Symbolic SymbolicANRLPrimitiveLattices(const string& lattice_and_centering, const char& space_group_letter){

    // Grab symbolic representation of primitive lattice in the ANRL convention.

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_SYMBOLIC_);

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " lattice and centering: " << lattice_and_centering << endl;
      cerr << __AFLOW_FUNC__ << " first character of space group symbol: " << space_group_letter << endl;
    }

    symbolic::Symbolic lattice("L",3,3);

    // ---------------------------------------------------------------------------
    // create symbolic list of characters
    symbolic::Symbolic a("a");
    symbolic::Symbolic b("b");
    symbolic::Symbolic c("c");
    symbolic::Symbolic cx("cx");
    symbolic::Symbolic cy("cy");
    symbolic::Symbolic cz("cz");
    symbolic::Symbolic beta("beta");
    symbolic::Symbolic gamma("gamma");

    // ---------------------------------------------------------------------------
    // triclinic (aP)
    if(lattice_and_centering == "aP"){
      lattice = ((a, _SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_),
          (b*cos(gamma), b*sin(gamma), _SYMBOLIC_ZERO_),
          (cx, cy, cz));
    }

    // ---------------------------------------------------------------------------
    // simple monoclinic (mP)
    if(lattice_and_centering == "mP"){
      lattice = ((a, _SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_),
          (_SYMBOLIC_ZERO_, b, _SYMBOLIC_ZERO_),
          (c*cos(beta), _SYMBOLIC_ZERO_, c*sin(beta)));
    }

    // ---------------------------------------------------------------------------
    // base-centered monoclinic (mC)
    if(lattice_and_centering == "mC"){
      lattice = ((0.5*a, -0.5*b, _SYMBOLIC_ZERO_), //DX20210106 - need parentheses and floats
          (0.5*a, 0.5*b, _SYMBOLIC_ZERO_),
          (c*cos(beta), _SYMBOLIC_ZERO_ , c*sin(beta)));
    }

    // ---------------------------------------------------------------------------
    // simple orthorhombic (oP)
    if(lattice_and_centering == "oP"){
      lattice = ((a, _SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_),
          (_SYMBOLIC_ZERO_, b, _SYMBOLIC_ZERO_),
          (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, c));
    }

    // ---------------------------------------------------------------------------
    // base-centered orthorhombic (oC)
    if(lattice_and_centering == "oC"){
      if(space_group_letter == 'A'){
        lattice = ((a, _SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_),
            (_SYMBOLIC_ZERO_, 0.5*b, -0.5*c),
            (_SYMBOLIC_ZERO_, 0.5*b, 0.5*c));
      }
      if(space_group_letter == 'C'){
        lattice = ((0.5*a, -0.5*b, _SYMBOLIC_ZERO_),
            (0.5*a, 0.5*b, _SYMBOLIC_ZERO_),
            (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, c));
      }
    }

    // ---------------------------------------------------------------------------
    // body-centered orthorhombic (oI)
    if(lattice_and_centering == "oI"){
      lattice = ((-0.5*a, 0.5*b, 0.5*c),
          (0.5*a, -0.5*b, 0.5*c),
          (0.5*a, 0.5*b, -0.5*c));
    }

    // ---------------------------------------------------------------------------
    // face-centered orthorhombic (oF)
    if(lattice_and_centering == "oF"){
      lattice = ((_SYMBOLIC_ZERO_, 0.5*b, 0.5*c),
          (0.5*a, _SYMBOLIC_ZERO_, 0.5*c),
          (0.5*a, 0.5*b, _SYMBOLIC_ZERO_));
    }

    // ---------------------------------------------------------------------------
    // simple tetragonal (tP)
    if(lattice_and_centering == "tP"){
      lattice = ((a, _SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_),
          (_SYMBOLIC_ZERO_, a, _SYMBOLIC_ZERO_),
          (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, c));
    }

    // ---------------------------------------------------------------------------
    // body-centered tetragonal (tI)
    if(lattice_and_centering == "tI"){
      lattice = ((-0.5*a, 0.5*a, 0.5*c),
          (0.5*a, -0.5*a, 0.5*c),
          (0.5*a, 0.5*a, -0.5*c));
    }

    // ---------------------------------------------------------------------------
    // hexagonal/trigonal (hP)
    if(lattice_and_centering == "hP"){
      lattice = ((0.5*a, -(sqrt(3.0)/2.0)*a, _SYMBOLIC_ZERO_),
          (0.5*a, (sqrt(3.0)/2.0)*a, _SYMBOLIC_ZERO_),
          (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, c));
    }

    // ---------------------------------------------------------------------------
    // rhombohedral (hR)
    if(lattice_and_centering == "hR"){
      lattice = ((0.5*a, -(1.0/(2.0*sqrt(3.0)))*a, (1.0/3.0)*c),
          (_SYMBOLIC_ZERO_, (1.0/sqrt(3.0))*a, (1.0/3.0)*c),
          (-0.5*a, -(1.0/(2.0*sqrt(3.0)))*a, (1.0/3.0)*c));
    }

    // ---------------------------------------------------------------------------
    // simple cubic (cP)
    if(lattice_and_centering == "cP"){
      lattice = ((a, _SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_),
          (_SYMBOLIC_ZERO_, a, _SYMBOLIC_ZERO_),
          (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, a));
    }

    // ---------------------------------------------------------------------------
    // body-centered cubic (cI)
    if(lattice_and_centering == "cI"){
      lattice = ((-0.5*a, 0.5*a, 0.5*a),
          (0.5*a, -0.5*a, 0.5*a),
          (0.5*a, 0.5*a, -0.5*a));
    }

    // ---------------------------------------------------------------------------
    // face-centered cubic (cF)
    if(lattice_and_centering == "cF"){
      lattice = ((_SYMBOLIC_ZERO_, 0.5*a, 0.5*a),
          (0.5*a, _SYMBOLIC_ZERO_, 0.5*a),
          (0.5*a, 0.5*a, _SYMBOLIC_ZERO_));
    }

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " lattice: " << lattice << endl;
    }

    return lattice;
  }
}

// *************************************************************************** 
// anrl::equations2SymbolicEquations()
// *************************************************************************** 
namespace anrl {
  vector<symbolic::Symbolic> equations2SymbolicEquations(const vector<vector<string> >& equations){

    // Convert equations (vector<vector<string> >) to symbolic notation.

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_SYMBOLIC_);
    stringstream message;

    vector<symbolic::Symbolic> symbolic_equations;

    for(uint i=0;i<equations.size();i++){
      symbolic::Symbolic position("pos", 3);
      if(equations[i].size()!=3){
        message << "Equation " << i << " does not have 3 coordinates (problem with ITC library coordinates): " << aurostd::joinWDelimiter(equations[i],",");
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_GENERIC_ERROR_);
      }
      for(uint j=0;j<equations[i].size();j++){
        if(equations[i][j] == "0"){
          position(j) = _SYMBOLIC_ZERO_;
        }
        else{
          position(j) = symbolic::string2symbolic(equations[i][j]);
        }
      }
      symbolic_equations.push_back(position);
    }

    if(LDEBUG){
      for(uint i=0;i<symbolic_equations.size();i++){
        cerr << __AFLOW_FUNC__ << " equations (string): " << aurostd::joinWDelimiter(equations[i],",") << " --> equations (Symbolic): " << symbolic_equations[i] << endl;
      }
    }

    return symbolic_equations;
  }
}

// *************************************************************************** 
// anrl::cartesian2lattice()
// *************************************************************************** 
namespace anrl {
  symbolic::Symbolic cartesian2lattice(const symbolic::Symbolic& lattice, const symbolic::Symbolic& cartesian_coordinate){

    // converts Cartesian coordinates to lattice coordinates
    // this is the symbolic math equivalent to C2F()

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_SYMBOLIC_);

    // ---------------------------------------------------------------------------
    // calculate volume (symbolic) 
    symbolic::Symbolic volume = lattice.row(0)|((lattice.row(1)%(lattice.row(2))).transpose()); // | is dot product, % is cross product
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " symbolic volume=" << volume << endl;
    }

    // ---------------------------------------------------------------------------
    // determine transformation matrix (vectors) 
    symbolic::Symbolic b1 = (lattice.row(1)%(lattice.row(2)))/volume;
    symbolic::Symbolic b2 = (lattice.row(2)%(lattice.row(0)))/volume;
    symbolic::Symbolic b3 = (lattice.row(0)%(lattice.row(1)))/volume;

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " transformation vectors b1=" << b1 << endl;
      cerr << __AFLOW_FUNC__ << " transformation vectors b2=" << b2 << endl;
      cerr << __AFLOW_FUNC__ << " transformation vectors b3=" << b3 << endl;
    }

    // ---------------------------------------------------------------------------
    // convert to lattice coordinate 
    symbolic::Symbolic lattice_coordinate("latt_coord",3);
    lattice_coordinate(0)=(cartesian_coordinate|b1).simplify();
    lattice_coordinate(1)=(cartesian_coordinate|b2).simplify();
    lattice_coordinate(2)=(cartesian_coordinate|b3).simplify();

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " cartesian coord=" << cartesian_coordinate << " --> " << lattice_coordinate << endl;
    }

    return lattice_coordinate;
  }
}

// *************************************************************************** 
// anrl::getXYZ2LatticeTransformation() 
// *************************************************************************** 
namespace anrl{
  symbolic::Symbolic getXYZ2LatticeTransformation(const string& lattice_and_centering){

    // gets transformation matrix from XYZ coordinates (ITC equations) to 
    // the Cartesian coordinates with respect to the lattice vectors.

    symbolic::Symbolic xyz2lattice("xyz2lattice",3,3);

    // ---------------------------------------------------------------------------
    // create symbolic list of characters
    symbolic::Symbolic a("a");
    symbolic::Symbolic b("b");
    symbolic::Symbolic c("c");
    symbolic::Symbolic cx("cx");
    symbolic::Symbolic cy("cy");
    symbolic::Symbolic cz("cz");
    symbolic::Symbolic beta("beta");
    symbolic::Symbolic gamma("gamma");

    // ---------------------------------------------------------------------------
    if(lattice_and_centering == "aP"){
      xyz2lattice = ((a, b*cos(gamma), cx),
          (_SYMBOLIC_ZERO_, b*sin(gamma), cy),
          (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, cz));
    } 

    else if(lattice_and_centering == "mP" || lattice_and_centering == "mC"){
      xyz2lattice = ((a, _SYMBOLIC_ZERO_, c*cos(beta)),
          (_SYMBOLIC_ZERO_, b, _SYMBOLIC_ZERO_),
          (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, c*sin(beta)));
    }
    else if(lattice_and_centering == "oP" || lattice_and_centering == "oC" || 
        lattice_and_centering == "oI" || lattice_and_centering == "oF"){
      xyz2lattice = ((a, _SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_),
          (_SYMBOLIC_ZERO_, b, _SYMBOLIC_ZERO_),
          (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, c));
    }
    else if(lattice_and_centering == "tP" || lattice_and_centering == "tI"){
      xyz2lattice = ((a, _SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_),
          (_SYMBOLIC_ZERO_, a, _SYMBOLIC_ZERO_),
          (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, c));
    }
    else if(lattice_and_centering == "hP" || lattice_and_centering == "hR"){
      xyz2lattice = ((a, _SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_),
          (_SYMBOLIC_ZERO_, a, _SYMBOLIC_ZERO_),
          (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, c));
    }
    else if(lattice_and_centering == "cP" || lattice_and_centering == "cF" || lattice_and_centering == "cI"){
      xyz2lattice = ((a, _SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_),
          (_SYMBOLIC_ZERO_, a, _SYMBOLIC_ZERO_),
          (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, a));
    }

    return xyz2lattice;
  }
}

// *************************************************************************** 
// anrl::getEquationsForCenteredLattices() 
// *************************************************************************** 
namespace anrl {
  vector<symbolic::Symbolic> getEquationsForCenteredLattices(const string& lattice_and_centering,
      const symbolic::Symbolic& lattice,
      const vector<symbolic::Symbolic>& conventional_equations){

    // convert equations to lattice equations for centered lattice (C, I, F).

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_SYMBOLIC_);
    stringstream message;

    vector<symbolic::Symbolic> lattice_equations;

    // ---------------------------------------------------------------------------
    // get symbolic transformation matrix for a particular lattice 
    symbolic::Symbolic xyz2lattice = getXYZ2LatticeTransformation(lattice_and_centering);
    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " symbolic transformation from xyz to lattice: " << xyz2lattice << endl; }

    // ---------------------------------------------------------------------------
    // transform symbolic coordinates 
    for(uint i=0;i<conventional_equations.size();i++){
      symbolic::Symbolic vec = xyz2lattice*conventional_equations[i];
      lattice_equations.push_back(cartesian2lattice(lattice, vec).simplify());
    }

    // ---------------------------------------------------------------------------
    // [LDEBUG] print conventional vs centered lattice equations
    if(LDEBUG){
      for(uint i=0;i<conventional_equations.size();i++){
        cerr << "conventional: " << conventional_equations[i] << " to centered lattice: " << lattice_equations[i] << endl;
      }
    }

    // ---------------------------------------------------------------------------
    // bring lattice equations in cell after transformation 
    for(uint i=0;i<lattice_equations.size();i++){
      lattice_equations[i] = symbolic::BringInCell(lattice_equations[i]);
    }

    // ---------------------------------------------------------------------------
    // [LDEBUG] print conventional vs centered lattice equations AFTER bring in cell
    if(LDEBUG){
      for(uint i=0;i<conventional_equations.size();i++){
        cerr << "conventional: " << conventional_equations[i] << " to centered lattice: " << lattice_equations[i] << endl;
      }
    }

    // ---------------------------------------------------------------------------
    // remove duplicates after bring in cell
    vector<symbolic::Symbolic> primitive_lattice_equations;
    for(uint i=0;i<lattice_equations.size();i++){
      bool is_unique_equation = true;
      for(uint j=0;j<primitive_lattice_equations.size();j++){
        if(symbolic::isEqualVector(primitive_lattice_equations[j],lattice_equations[i])){
          is_unique_equation = false;
          break;
        }
      }
      if(is_unique_equation){ primitive_lattice_equations.push_back(lattice_equations[i]); }
    }

    if(LDEBUG){
      cerr << "# equations before: " << lattice_equations.size() << " vs # equations after: " << primitive_lattice_equations.size() << endl;
    }

    return primitive_lattice_equations;
  }
}

// *************************************************************************** 
// anrl::getEquationsForCenteredLattices() 
// *************************************************************************** 
namespace anrl {
  vector<symbolic::Symbolic> convertEquations2FractionalEquations(const string& lattice_and_centering,
      const symbolic::Symbolic& lattice,
      const vector<symbolic::Symbolic> conventional_equations){

    // convert equations to lattice equations (fractional).

    vector<symbolic::Symbolic> transformed_equations;

    // ---------------------------------------------------------------------------
    // if C-, I-, F-centered, then the equations from ITC are xyz coordinates
    if(lattice_and_centering == "mC" || lattice_and_centering == "oC" || lattice_and_centering == "oI" ||
        lattice_and_centering == "oF" || lattice_and_centering == "tI" || lattice_and_centering == "cI" ||
        lattice_and_centering == "cF"){
      transformed_equations = getEquationsForCenteredLattices(lattice_and_centering, lattice, conventional_equations);  
    }

    // ---------------------------------------------------------------------------
    // if hexagonal/trigonal/rhombohedral, then 
    else if(lattice_and_centering == "hP" || lattice_and_centering == "hR"){
      transformed_equations = conventional_equations;
    }
    // ---------------------------------------------------------------------------
    // else 
    else{
      transformed_equations = conventional_equations;
    }

    return transformed_equations;
  }
}

// *************************************************************************** 
// symbolic::addSymbolicEquation2Atoms() 
// *************************************************************************** 
namespace anrl {
  void addSymbolicEquation2Atoms(const vector<symbolic::Symbolic>& equations, deque<_atom>& atoms, bool isfpos){

    // add symbolic equations to atom.fpos_equation or atom.cpos_equation
    // converts from variables from Symbolic to string
    // DEFAULT: update atom.fpos_equation

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_SYMBOLIC_);
    stringstream message;

    // ---------------------------------------------------------------------------
    // ensure sizes of atoms and symbolic equations match 
    if(equations.size() != atoms.size()){
      message << "The number of equations and atoms do not match. Check tolerances. #equations=" << equations.size() << ", #atoms=" << atoms.size();
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
    }

    // ---------------------------------------------------------------------------
    // ensure sizes of atoms and symbolic equations match 
    for(uint i=0;i<atoms.size();i++){
      vector<string> coordinate;
      for(uint j=0;j<3;j++){
        stringstream ss_pos; ss_pos << equations[i].row(j);
        coordinate.push_back(ss_pos.str());
      }
      if(isfpos){ atoms[i].fpos_equation = coordinate; }
      else{ atoms[i].cpos_equation = coordinate; }
    }

    // ---------------------------------------------------------------------------
    // [LDEBUG] print string (fpos/cpos) version of symbolic equation
    if(LDEBUG){
      if(isfpos){
        for(uint i=0;i<atoms.size();i++){
          cerr << __AFLOW_FUNC__ << " fpos_equation=" << aurostd::joinWDelimiter(atoms[i].fpos_equation,",") << endl;
        } 
      }
      else{
        for(uint i=0;i<atoms.size();i++){
          cerr << __AFLOW_FUNC__ << " cpos_equation=" << aurostd::joinWDelimiter(atoms[i].cpos_equation,",") << endl;
        } 
      }
    }
  }
}
#endif // USE_SYMBOLIC_SOURCE

#endif // _AFLOW_SYMBOLIC_CPP
