// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - David Hicks - 2018
// FILE "ANRL/aflow_anrl_A12B7C2_hP21_174_2j2k_ajk_cf.cpp"

#ifndef _AFLOW_ANRL_A12B7C2_hP21_174_2j2k_ajk_cf_CPP // AFLOW_REMOVE_GREP
#define _AFLOW_ANRL_A12B7C2_hP21_174_2j2k_ajk_cf_CPP // AFLOW_REMOVE_GREP
#include "../aflow.h" // AFLOW_REMOVE_GREP

namespace anrl {
  uint WebANRL_A12B7C2_hP21_174_2j2k_ajk_cf(stringstream &web,bool LDEBUG);
}

namespace anrl {
  uint PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf(stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG) {
    // system A12B7C2_hP21_174_2j2k_ajk_cf

    if(XHOST.vflag_control.flag("WWW")) {
      WebANRL_A12B7C2_hP21_174_2j2k_ajk_cf(web,LDEBUG); // PLUG WEB STUFF
#ifdef _ANRL_NOWEB_
      web << "no web";
      cout << web.str() << endl;
#else
      cout << web.str() << endl;
#endif
      return 0; //DX20200727
    }

    vector<double> vparameters;
    aurostd::string2tokens(parameters,vparameters,",");

    uint nspecies,natoms,spacegroup,nunderscores,nparameters;
    string label,Pearson_symbol,params,Strukturbericht,prototype,dialect;

    anrl::vproto2tokens(proto_line,label,nspecies,natoms,spacegroup,nunderscores,nparameters,Pearson_symbol,params,Strukturbericht,prototype,dialect);

    anrl::PrototypeANRL_Consistency(vparameters.size(),nparameters,prototype,label,
        Strukturbericht,Pearson_symbol,spacegroup,params,print_mode);    

    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: FOUND" << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: label=" << label << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: nspecies=" << nspecies << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: natoms=" << natoms << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: spacegroup=" << spacegroup << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: nunderscores=" << nunderscores << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: nparameters=" <<  nparameters << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: Pearson_symbol=" << Pearson_symbol << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: params=" << params << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: Strukturbericht=" << Strukturbericht << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: prototype=" << prototype << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: dialect=" << dialect << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: vparameters.size()=" << vparameters.size() << endl;}

    xvector<double> xn(3);   xn(1)=1.0;xn(2)=0.0;xn(3)=0.0;
    xvector<double> yn(3);   yn(1)=0.0;yn(2)=1.0;yn(3)=0.0;
    xvector<double> zn(3);   zn(1)=0.0;zn(2)=0.0;zn(3)=1.0;
    xvector<double> a1(3),a2(3),a3(3);

    if(print_mode==1 && vparameters.size()==0){
      for(uint n=0;n<nparameters;n++){
        vparameters.push_back(0);
      }
    }

    uint i=0;
    double a=vparameters.at(i++);                  if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: a=" << a << endl;}
    double covera=vparameters.at(i++),c=covera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: c=" << c << " (c/a=" << covera << ")" << endl;}

    double x4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: x4=" << x4 << endl;}
    double y4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: y4=" << y4 << endl;}
    double x5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: x5=" << x5 << endl;}
    double y5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: y5=" << y5 << endl;}
    double x6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: x6=" << x6 << endl;}
    double y6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: y6=" << y6 << endl;}
    double x7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: x7=" << x7 << endl;}
    double y7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: y7=" << y7 << endl;}
    double x8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: x8=" << x8 << endl;}
    double y8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: y8=" << y8 << endl;}
    double x9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: x9=" << x9 << endl;}
    double y9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf: y9=" << y9 << endl;}

    str.iomode=IOVASP_AUTO;
    str.title=label+" params="+parameters+" SG="+aurostd::utype2string(spacegroup)+DOI_ANRL; //CO20190520
    str.scale=1.0;

    a1=(1.0/2.0)*a*xn-(sqrt(3.0)/2.0)*a*yn;
    a2=(1.0/2.0)*a*xn+(sqrt(3.0)/2.0)*a*yn;
    a3=c*zn;

    str.lattice(1,1)=a1(1);str.lattice(1,2)=a1(2);str.lattice(1,3)=a1(3);
    str.lattice(2,1)=a2(1);str.lattice(2,2)=a2(2);str.lattice(2,3)=a2(3);
    str.lattice(3,1)=a3(1);str.lattice(3,2)=a3(2);str.lattice(3,3)=a3(3);

    // symbolic representation of lattice vectors
    vector<string> a1_equation, a2_equation, a3_equation;
    a1_equation.push_back("(1.0/2.0)*a");a1_equation.push_back("-(sqrt(3.0)/2.0)*a");a1_equation.push_back("0");
    a2_equation.push_back("(1.0/2.0)*a");a2_equation.push_back("(sqrt(3.0)/2.0)*a");a2_equation.push_back("0");
    a3_equation.push_back("0");a3_equation.push_back("0");a3_equation.push_back("c");
    str.symbolic_math_lattice.push_back(a1_equation);
    str.symbolic_math_lattice.push_back(a2_equation);
    str.symbolic_math_lattice.push_back(a3_equation);

    str.num_lattice_parameters = 2;

    str.num_parameters = vparameters.size();
    vector<string> parameter_list; aurostd::string2tokens(params,parameter_list,",");
    str.prototype_parameter_list = parameter_list;
    str.prototype_parameter_values = vparameters;

    if(print_mode!=1){
      str.FixLattices(); // Reciprocal/f2c/c2f
    }

    _atom atom;

    atom.name="A"; atom.type=0;                                       // atom B4
    atom.fpos(1)=x4;atom.fpos(2)=y4;atom.fpos(3)=0.0;                     // atom B4
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x4");atom.fpos_equation.push_back("y4");atom.fpos_equation.push_back("0.0");// atom B4 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B4 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B4

    atom.name="A"; atom.type=0;                                       // atom B5
    atom.fpos(1)=-y4;atom.fpos(2)=(x4-y4);atom.fpos(3)=0.0;                     // atom B5
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y4");atom.fpos_equation.push_back("(x4-y4)");atom.fpos_equation.push_back("0.0");// atom B5 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B5 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B5

    atom.name="A"; atom.type=0;                                       // atom B6
    atom.fpos(1)=(-x4+y4);atom.fpos(2)=-x4;atom.fpos(3)=0.0;                     // atom B6
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x4+y4)");atom.fpos_equation.push_back("-x4");atom.fpos_equation.push_back("0.0");// atom B6 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B6 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B6

    atom.name="A"; atom.type=0;                                       // atom B7
    atom.fpos(1)=x5;atom.fpos(2)=y5;atom.fpos(3)=0.0;                     // atom B7
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("y5");atom.fpos_equation.push_back("0.0");// atom B7 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B7 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B7

    atom.name="A"; atom.type=0;                                       // atom B8
    atom.fpos(1)=-y5;atom.fpos(2)=(x5-y5);atom.fpos(3)=0.0;                     // atom B8
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y5");atom.fpos_equation.push_back("(x5-y5)");atom.fpos_equation.push_back("0.0");// atom B8 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B8 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B8

    atom.name="A"; atom.type=0;                                       // atom B9
    atom.fpos(1)=(-x5+y5);atom.fpos(2)=-x5;atom.fpos(3)=0.0;                     // atom B9
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x5+y5)");atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("0.0");// atom B9 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B9 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B9

    atom.name="A"; atom.type=0;                                       // atom B13
    atom.fpos(1)=x7;atom.fpos(2)=y7;atom.fpos(3)=(1.0/2.0);                     // atom B13
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x7");atom.fpos_equation.push_back("y7");atom.fpos_equation.push_back("(1.0/2.0)");// atom B13 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B13 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B13

    atom.name="A"; atom.type=0;                                       // atom B14
    atom.fpos(1)=-y7;atom.fpos(2)=(x7-y7);atom.fpos(3)=(1.0/2.0);                     // atom B14
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y7");atom.fpos_equation.push_back("(x7-y7)");atom.fpos_equation.push_back("(1.0/2.0)");// atom B14 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B14 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B14

    atom.name="A"; atom.type=0;                                       // atom B15
    atom.fpos(1)=(-x7+y7);atom.fpos(2)=-x7;atom.fpos(3)=(1.0/2.0);                     // atom B15
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x7+y7)");atom.fpos_equation.push_back("-x7");atom.fpos_equation.push_back("(1.0/2.0)");// atom B15 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B15 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B15

    atom.name="A"; atom.type=0;                                       // atom B16
    atom.fpos(1)=x8;atom.fpos(2)=y8;atom.fpos(3)=(1.0/2.0);                     // atom B16
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x8");atom.fpos_equation.push_back("y8");atom.fpos_equation.push_back("(1.0/2.0)");// atom B16 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B16 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B16

    atom.name="A"; atom.type=0;                                       // atom B17
    atom.fpos(1)=-y8;atom.fpos(2)=(x8-y8);atom.fpos(3)=(1.0/2.0);                     // atom B17
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y8");atom.fpos_equation.push_back("(x8-y8)");atom.fpos_equation.push_back("(1.0/2.0)");// atom B17 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B17 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B17

    atom.name="A"; atom.type=0;                                       // atom B18
    atom.fpos(1)=(-x8+y8);atom.fpos(2)=-x8;atom.fpos(3)=(1.0/2.0);                     // atom B18
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x8+y8)");atom.fpos_equation.push_back("-x8");atom.fpos_equation.push_back("(1.0/2.0)");// atom B18 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B18 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B18

    atom.name="B"; atom.type=1;                                       // atom B1
    atom.fpos(1)=0.0;atom.fpos(2)=0.0;atom.fpos(3)=0.0;                     // atom B1
    atom.fpos_equation.clear();atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("0.0");// atom B1 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B1 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B1

    atom.name="B"; atom.type=1;                                       // atom B10
    atom.fpos(1)=x6;atom.fpos(2)=y6;atom.fpos(3)=0.0;                     // atom B10
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("y6");atom.fpos_equation.push_back("0.0");// atom B10 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B10 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B10

    atom.name="B"; atom.type=1;                                       // atom B11
    atom.fpos(1)=-y6;atom.fpos(2)=(x6-y6);atom.fpos(3)=0.0;                     // atom B11
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y6");atom.fpos_equation.push_back("(x6-y6)");atom.fpos_equation.push_back("0.0");// atom B11 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B11 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B11

    atom.name="B"; atom.type=1;                                       // atom B12
    atom.fpos(1)=(-x6+y6);atom.fpos(2)=-x6;atom.fpos(3)=0.0;                     // atom B12
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x6+y6)");atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("0.0");// atom B12 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B12 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B12

    atom.name="B"; atom.type=1;                                       // atom B19
    atom.fpos(1)=x9;atom.fpos(2)=y9;atom.fpos(3)=(1.0/2.0);                     // atom B19
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x9");atom.fpos_equation.push_back("y9");atom.fpos_equation.push_back("(1.0/2.0)");// atom B19 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B19 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B19

    atom.name="B"; atom.type=1;                                       // atom B20
    atom.fpos(1)=-y9;atom.fpos(2)=(x9-y9);atom.fpos(3)=(1.0/2.0);                     // atom B20
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y9");atom.fpos_equation.push_back("(x9-y9)");atom.fpos_equation.push_back("(1.0/2.0)");// atom B20 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B20 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B20

    atom.name="B"; atom.type=1;                                       // atom B21
    atom.fpos(1)=(-x9+y9);atom.fpos(2)=-x9;atom.fpos(3)=(1.0/2.0);                     // atom B21
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x9+y9)");atom.fpos_equation.push_back("-x9");atom.fpos_equation.push_back("(1.0/2.0)");// atom B21 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B21 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B21

    atom.name="C"; atom.type=2;                                       // atom B2
    atom.fpos(1)=(1.0/3.0);atom.fpos(2)=(2.0/3.0);atom.fpos(3)=0.0;                     // atom B2
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/3.0)");atom.fpos_equation.push_back("(2.0/3.0)");atom.fpos_equation.push_back("0.0");// atom B2 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B2 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B2

    atom.name="C"; atom.type=2;                                       // atom B3
    atom.fpos(1)=(2.0/3.0);atom.fpos(2)=(1.0/3.0);atom.fpos(3)=(1.0/2.0);                     // atom B3
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(2.0/3.0)");atom.fpos_equation.push_back("(1.0/3.0)");atom.fpos_equation.push_back("(1.0/2.0)");// atom B3 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B3 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B3


    return str.atoms.size();  
  }
} // namespace anrl

namespace anrl {
  uint WebANRL_A12B7C2_hP21_174_2j2k_ajk_cf(stringstream& web,bool LDEBUG) {
#ifndef _ANRL_NOWEB_
#endif

    if(LDEBUG) {cerr << "anrl:: WebANRL_A12B7C2_hP21_174_2j2k_ajk_cf: web.str().size()=" << web.str().size() << endl;}

    return web.str().size();
  }
} // namespace anrl

#endif // AFLOW_REMOVE_GREP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

