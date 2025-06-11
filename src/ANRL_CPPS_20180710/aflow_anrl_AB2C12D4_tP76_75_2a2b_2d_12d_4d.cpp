// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - David Hicks - 2018
// FILE "ANRL/aflow_anrl_AB2C12D4_tP76_75_2a2b_2d_12d_4d.cpp"

#ifndef _AFLOW_ANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d_CPP // AFLOW_REMOVE_GREP
#define _AFLOW_ANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d_CPP // AFLOW_REMOVE_GREP
#include "../aflow.h" // AFLOW_REMOVE_GREP

namespace anrl {
  uint WebANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d(stringstream &web,bool LDEBUG);
}

namespace anrl {
  uint PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d(stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG) {
    // system AB2C12D4_tP76_75_2a2b_2d_12d_4d

    if(XHOST.vflag_control.flag("WWW")) {
      WebANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d(web,LDEBUG); // PLUG WEB STUFF
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

    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: FOUND" << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: label=" << label << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: nspecies=" << nspecies << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: natoms=" << natoms << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: spacegroup=" << spacegroup << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: nunderscores=" << nunderscores << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: nparameters=" <<  nparameters << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: Pearson_symbol=" << Pearson_symbol << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: params=" << params << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: Strukturbericht=" << Strukturbericht << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: prototype=" << prototype << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: dialect=" << dialect << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: vparameters.size()=" << vparameters.size() << endl;}

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
    double a=vparameters.at(i++);                  if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: a=" << a << endl;}
    double covera=vparameters.at(i++),c=covera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: c=" << c << " (c/a=" << covera << ")" << endl;}

    double z1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z1=" << z1 << endl;}
    double z2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z2=" << z2 << endl;}
    double z3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z3=" << z3 << endl;}
    double z4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z4=" << z4 << endl;}
    double x5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: x5=" << x5 << endl;}
    double y5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: y5=" << y5 << endl;}
    double z5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z5=" << z5 << endl;}
    double x6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: x6=" << x6 << endl;}
    double y6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: y6=" << y6 << endl;}
    double z6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z6=" << z6 << endl;}
    double x7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: x7=" << x7 << endl;}
    double y7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: y7=" << y7 << endl;}
    double z7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z7=" << z7 << endl;}
    double x8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: x8=" << x8 << endl;}
    double y8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: y8=" << y8 << endl;}
    double z8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z8=" << z8 << endl;}
    double x9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: x9=" << x9 << endl;}
    double y9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: y9=" << y9 << endl;}
    double z9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z9=" << z9 << endl;}
    double x10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: x10=" << x10 << endl;}
    double y10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: y10=" << y10 << endl;}
    double z10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z10=" << z10 << endl;}
    double x11=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: x11=" << x11 << endl;}
    double y11=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: y11=" << y11 << endl;}
    double z11=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z11=" << z11 << endl;}
    double x12=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: x12=" << x12 << endl;}
    double y12=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: y12=" << y12 << endl;}
    double z12=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z12=" << z12 << endl;}
    double x13=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: x13=" << x13 << endl;}
    double y13=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: y13=" << y13 << endl;}
    double z13=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z13=" << z13 << endl;}
    double x14=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: x14=" << x14 << endl;}
    double y14=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: y14=" << y14 << endl;}
    double z14=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z14=" << z14 << endl;}
    double x15=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: x15=" << x15 << endl;}
    double y15=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: y15=" << y15 << endl;}
    double z15=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z15=" << z15 << endl;}
    double x16=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: x16=" << x16 << endl;}
    double y16=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: y16=" << y16 << endl;}
    double z16=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z16=" << z16 << endl;}
    double x17=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: x17=" << x17 << endl;}
    double y17=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: y17=" << y17 << endl;}
    double z17=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z17=" << z17 << endl;}
    double x18=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: x18=" << x18 << endl;}
    double y18=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: y18=" << y18 << endl;}
    double z18=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z18=" << z18 << endl;}
    double x19=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: x19=" << x19 << endl;}
    double y19=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: y19=" << y19 << endl;}
    double z19=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z19=" << z19 << endl;}
    double x20=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: x20=" << x20 << endl;}
    double y20=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: y20=" << y20 << endl;}
    double z20=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z20=" << z20 << endl;}
    double x21=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: x21=" << x21 << endl;}
    double y21=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: y21=" << y21 << endl;}
    double z21=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z21=" << z21 << endl;}
    double x22=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: x22=" << x22 << endl;}
    double y22=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: y22=" << y22 << endl;}
    double z22=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: z22=" << z22 << endl;}

    str.iomode=IOVASP_AUTO;
    str.title=label+" params="+parameters+" SG="+aurostd::utype2string(spacegroup)+DOI_ANRL; //CO20190520
    str.scale=1.0;

    a1=a*xn;
    a2=a*yn;
    a3=c*zn;

    str.lattice(1,1)=a1(1);str.lattice(1,2)=a1(2);str.lattice(1,3)=a1(3);
    str.lattice(2,1)=a2(1);str.lattice(2,2)=a2(2);str.lattice(2,3)=a2(3);
    str.lattice(3,1)=a3(1);str.lattice(3,2)=a3(2);str.lattice(3,3)=a3(3);

    // symbolic representation of lattice vectors
    vector<string> a1_equation, a2_equation, a3_equation;
    a1_equation.push_back("a");a1_equation.push_back("0");a1_equation.push_back("0");
    a2_equation.push_back("0");a2_equation.push_back("a");a2_equation.push_back("0");
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

    atom.name="A"; atom.type=0;                                       // atom B1
    atom.fpos(1)=0.0;atom.fpos(2)=0.0;atom.fpos(3)=z1;                     // atom B1
    atom.fpos_equation.clear();atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("z1");// atom B1 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B1 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B1

    atom.name="A"; atom.type=0;                                       // atom B2
    atom.fpos(1)=0.0;atom.fpos(2)=0.0;atom.fpos(3)=z2;                     // atom B2
    atom.fpos_equation.clear();atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("z2");// atom B2 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B2 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B2

    atom.name="A"; atom.type=0;                                       // atom B3
    atom.fpos(1)=(1.0/2.0);atom.fpos(2)=(1.0/2.0);atom.fpos(3)=z3;                     // atom B3
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("z3");// atom B3 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B3 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B3

    atom.name="A"; atom.type=0;                                       // atom B4
    atom.fpos(1)=(1.0/2.0);atom.fpos(2)=(1.0/2.0);atom.fpos(3)=z4;                     // atom B4
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("z4");// atom B4 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B4 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B4

    atom.name="B"; atom.type=1;                                       // atom B5
    atom.fpos(1)=x5;atom.fpos(2)=y5;atom.fpos(3)=z5;                     // atom B5
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("y5");atom.fpos_equation.push_back("z5");// atom B5 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B5 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B5

    atom.name="B"; atom.type=1;                                       // atom B6
    atom.fpos(1)=-x5;atom.fpos(2)=-y5;atom.fpos(3)=z5;                     // atom B6
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("-y5");atom.fpos_equation.push_back("z5");// atom B6 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B6 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B6

    atom.name="B"; atom.type=1;                                       // atom B7
    atom.fpos(1)=-y5;atom.fpos(2)=x5;atom.fpos(3)=z5;                     // atom B7
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y5");atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("z5");// atom B7 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B7 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B7

    atom.name="B"; atom.type=1;                                       // atom B8
    atom.fpos(1)=y5;atom.fpos(2)=-x5;atom.fpos(3)=z5;                     // atom B8
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y5");atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("z5");// atom B8 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B8 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B8

    atom.name="B"; atom.type=1;                                       // atom B9
    atom.fpos(1)=x6;atom.fpos(2)=y6;atom.fpos(3)=z6;                     // atom B9
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("y6");atom.fpos_equation.push_back("z6");// atom B9 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B9 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B9

    atom.name="B"; atom.type=1;                                       // atom B10
    atom.fpos(1)=-x6;atom.fpos(2)=-y6;atom.fpos(3)=z6;                     // atom B10
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("-y6");atom.fpos_equation.push_back("z6");// atom B10 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B10 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B10

    atom.name="B"; atom.type=1;                                       // atom B11
    atom.fpos(1)=-y6;atom.fpos(2)=x6;atom.fpos(3)=z6;                     // atom B11
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y6");atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("z6");// atom B11 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B11 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B11

    atom.name="B"; atom.type=1;                                       // atom B12
    atom.fpos(1)=y6;atom.fpos(2)=-x6;atom.fpos(3)=z6;                     // atom B12
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y6");atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("z6");// atom B12 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B12 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B12

    atom.name="C"; atom.type=2;                                       // atom B13
    atom.fpos(1)=x7;atom.fpos(2)=y7;atom.fpos(3)=z7;                     // atom B13
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x7");atom.fpos_equation.push_back("y7");atom.fpos_equation.push_back("z7");// atom B13 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B13 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B13

    atom.name="C"; atom.type=2;                                       // atom B14
    atom.fpos(1)=-x7;atom.fpos(2)=-y7;atom.fpos(3)=z7;                     // atom B14
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x7");atom.fpos_equation.push_back("-y7");atom.fpos_equation.push_back("z7");// atom B14 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B14 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B14

    atom.name="C"; atom.type=2;                                       // atom B15
    atom.fpos(1)=-y7;atom.fpos(2)=x7;atom.fpos(3)=z7;                     // atom B15
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y7");atom.fpos_equation.push_back("x7");atom.fpos_equation.push_back("z7");// atom B15 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B15 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B15

    atom.name="C"; atom.type=2;                                       // atom B16
    atom.fpos(1)=y7;atom.fpos(2)=-x7;atom.fpos(3)=z7;                     // atom B16
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y7");atom.fpos_equation.push_back("-x7");atom.fpos_equation.push_back("z7");// atom B16 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B16 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B16

    atom.name="C"; atom.type=2;                                       // atom B17
    atom.fpos(1)=x8;atom.fpos(2)=y8;atom.fpos(3)=z8;                     // atom B17
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x8");atom.fpos_equation.push_back("y8");atom.fpos_equation.push_back("z8");// atom B17 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B17 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B17

    atom.name="C"; atom.type=2;                                       // atom B18
    atom.fpos(1)=-x8;atom.fpos(2)=-y8;atom.fpos(3)=z8;                     // atom B18
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x8");atom.fpos_equation.push_back("-y8");atom.fpos_equation.push_back("z8");// atom B18 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B18 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B18

    atom.name="C"; atom.type=2;                                       // atom B19
    atom.fpos(1)=-y8;atom.fpos(2)=x8;atom.fpos(3)=z8;                     // atom B19
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y8");atom.fpos_equation.push_back("x8");atom.fpos_equation.push_back("z8");// atom B19 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B19 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B19

    atom.name="C"; atom.type=2;                                       // atom B20
    atom.fpos(1)=y8;atom.fpos(2)=-x8;atom.fpos(3)=z8;                     // atom B20
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y8");atom.fpos_equation.push_back("-x8");atom.fpos_equation.push_back("z8");// atom B20 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B20 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B20

    atom.name="C"; atom.type=2;                                       // atom B21
    atom.fpos(1)=x9;atom.fpos(2)=y9;atom.fpos(3)=z9;                     // atom B21
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x9");atom.fpos_equation.push_back("y9");atom.fpos_equation.push_back("z9");// atom B21 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B21 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B21

    atom.name="C"; atom.type=2;                                       // atom B22
    atom.fpos(1)=-x9;atom.fpos(2)=-y9;atom.fpos(3)=z9;                     // atom B22
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x9");atom.fpos_equation.push_back("-y9");atom.fpos_equation.push_back("z9");// atom B22 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B22 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B22

    atom.name="C"; atom.type=2;                                       // atom B23
    atom.fpos(1)=-y9;atom.fpos(2)=x9;atom.fpos(3)=z9;                     // atom B23
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y9");atom.fpos_equation.push_back("x9");atom.fpos_equation.push_back("z9");// atom B23 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B23 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B23

    atom.name="C"; atom.type=2;                                       // atom B24
    atom.fpos(1)=y9;atom.fpos(2)=-x9;atom.fpos(3)=z9;                     // atom B24
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y9");atom.fpos_equation.push_back("-x9");atom.fpos_equation.push_back("z9");// atom B24 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B24 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B24

    atom.name="C"; atom.type=2;                                       // atom B25
    atom.fpos(1)=x10;atom.fpos(2)=y10;atom.fpos(3)=z10;                     // atom B25
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x10");atom.fpos_equation.push_back("y10");atom.fpos_equation.push_back("z10");// atom B25 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B25 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B25

    atom.name="C"; atom.type=2;                                       // atom B26
    atom.fpos(1)=-x10;atom.fpos(2)=-y10;atom.fpos(3)=z10;                     // atom B26
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x10");atom.fpos_equation.push_back("-y10");atom.fpos_equation.push_back("z10");// atom B26 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B26 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B26

    atom.name="C"; atom.type=2;                                       // atom B27
    atom.fpos(1)=-y10;atom.fpos(2)=x10;atom.fpos(3)=z10;                     // atom B27
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y10");atom.fpos_equation.push_back("x10");atom.fpos_equation.push_back("z10");// atom B27 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B27 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B27

    atom.name="C"; atom.type=2;                                       // atom B28
    atom.fpos(1)=y10;atom.fpos(2)=-x10;atom.fpos(3)=z10;                     // atom B28
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y10");atom.fpos_equation.push_back("-x10");atom.fpos_equation.push_back("z10");// atom B28 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B28 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B28

    atom.name="C"; atom.type=2;                                       // atom B29
    atom.fpos(1)=x11;atom.fpos(2)=y11;atom.fpos(3)=z11;                     // atom B29
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x11");atom.fpos_equation.push_back("y11");atom.fpos_equation.push_back("z11");// atom B29 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B29 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B29

    atom.name="C"; atom.type=2;                                       // atom B30
    atom.fpos(1)=-x11;atom.fpos(2)=-y11;atom.fpos(3)=z11;                     // atom B30
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x11");atom.fpos_equation.push_back("-y11");atom.fpos_equation.push_back("z11");// atom B30 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B30 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B30

    atom.name="C"; atom.type=2;                                       // atom B31
    atom.fpos(1)=-y11;atom.fpos(2)=x11;atom.fpos(3)=z11;                     // atom B31
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y11");atom.fpos_equation.push_back("x11");atom.fpos_equation.push_back("z11");// atom B31 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B31 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B31

    atom.name="C"; atom.type=2;                                       // atom B32
    atom.fpos(1)=y11;atom.fpos(2)=-x11;atom.fpos(3)=z11;                     // atom B32
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y11");atom.fpos_equation.push_back("-x11");atom.fpos_equation.push_back("z11");// atom B32 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B32 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B32

    atom.name="C"; atom.type=2;                                       // atom B33
    atom.fpos(1)=x12;atom.fpos(2)=y12;atom.fpos(3)=z12;                     // atom B33
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x12");atom.fpos_equation.push_back("y12");atom.fpos_equation.push_back("z12");// atom B33 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B33 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B33

    atom.name="C"; atom.type=2;                                       // atom B34
    atom.fpos(1)=-x12;atom.fpos(2)=-y12;atom.fpos(3)=z12;                     // atom B34
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x12");atom.fpos_equation.push_back("-y12");atom.fpos_equation.push_back("z12");// atom B34 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B34 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B34

    atom.name="C"; atom.type=2;                                       // atom B35
    atom.fpos(1)=-y12;atom.fpos(2)=x12;atom.fpos(3)=z12;                     // atom B35
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y12");atom.fpos_equation.push_back("x12");atom.fpos_equation.push_back("z12");// atom B35 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B35 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B35

    atom.name="C"; atom.type=2;                                       // atom B36
    atom.fpos(1)=y12;atom.fpos(2)=-x12;atom.fpos(3)=z12;                     // atom B36
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y12");atom.fpos_equation.push_back("-x12");atom.fpos_equation.push_back("z12");// atom B36 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B36 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B36

    atom.name="C"; atom.type=2;                                       // atom B37
    atom.fpos(1)=x13;atom.fpos(2)=y13;atom.fpos(3)=z13;                     // atom B37
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x13");atom.fpos_equation.push_back("y13");atom.fpos_equation.push_back("z13");// atom B37 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B37 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B37

    atom.name="C"; atom.type=2;                                       // atom B38
    atom.fpos(1)=-x13;atom.fpos(2)=-y13;atom.fpos(3)=z13;                     // atom B38
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x13");atom.fpos_equation.push_back("-y13");atom.fpos_equation.push_back("z13");// atom B38 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B38 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B38

    atom.name="C"; atom.type=2;                                       // atom B39
    atom.fpos(1)=-y13;atom.fpos(2)=x13;atom.fpos(3)=z13;                     // atom B39
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y13");atom.fpos_equation.push_back("x13");atom.fpos_equation.push_back("z13");// atom B39 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B39 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B39

    atom.name="C"; atom.type=2;                                       // atom B40
    atom.fpos(1)=y13;atom.fpos(2)=-x13;atom.fpos(3)=z13;                     // atom B40
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y13");atom.fpos_equation.push_back("-x13");atom.fpos_equation.push_back("z13");// atom B40 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B40 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B40

    atom.name="C"; atom.type=2;                                       // atom B41
    atom.fpos(1)=x14;atom.fpos(2)=y14;atom.fpos(3)=z14;                     // atom B41
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x14");atom.fpos_equation.push_back("y14");atom.fpos_equation.push_back("z14");// atom B41 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B41 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B41

    atom.name="C"; atom.type=2;                                       // atom B42
    atom.fpos(1)=-x14;atom.fpos(2)=-y14;atom.fpos(3)=z14;                     // atom B42
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x14");atom.fpos_equation.push_back("-y14");atom.fpos_equation.push_back("z14");// atom B42 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B42 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B42

    atom.name="C"; atom.type=2;                                       // atom B43
    atom.fpos(1)=-y14;atom.fpos(2)=x14;atom.fpos(3)=z14;                     // atom B43
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y14");atom.fpos_equation.push_back("x14");atom.fpos_equation.push_back("z14");// atom B43 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B43 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B43

    atom.name="C"; atom.type=2;                                       // atom B44
    atom.fpos(1)=y14;atom.fpos(2)=-x14;atom.fpos(3)=z14;                     // atom B44
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y14");atom.fpos_equation.push_back("-x14");atom.fpos_equation.push_back("z14");// atom B44 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B44 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B44

    atom.name="C"; atom.type=2;                                       // atom B45
    atom.fpos(1)=x15;atom.fpos(2)=y15;atom.fpos(3)=z15;                     // atom B45
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x15");atom.fpos_equation.push_back("y15");atom.fpos_equation.push_back("z15");// atom B45 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B45 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B45

    atom.name="C"; atom.type=2;                                       // atom B46
    atom.fpos(1)=-x15;atom.fpos(2)=-y15;atom.fpos(3)=z15;                     // atom B46
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x15");atom.fpos_equation.push_back("-y15");atom.fpos_equation.push_back("z15");// atom B46 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B46 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B46

    atom.name="C"; atom.type=2;                                       // atom B47
    atom.fpos(1)=-y15;atom.fpos(2)=x15;atom.fpos(3)=z15;                     // atom B47
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y15");atom.fpos_equation.push_back("x15");atom.fpos_equation.push_back("z15");// atom B47 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B47 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B47

    atom.name="C"; atom.type=2;                                       // atom B48
    atom.fpos(1)=y15;atom.fpos(2)=-x15;atom.fpos(3)=z15;                     // atom B48
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y15");atom.fpos_equation.push_back("-x15");atom.fpos_equation.push_back("z15");// atom B48 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B48 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B48

    atom.name="C"; atom.type=2;                                       // atom B49
    atom.fpos(1)=x16;atom.fpos(2)=y16;atom.fpos(3)=z16;                     // atom B49
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x16");atom.fpos_equation.push_back("y16");atom.fpos_equation.push_back("z16");// atom B49 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B49 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B49

    atom.name="C"; atom.type=2;                                       // atom B50
    atom.fpos(1)=-x16;atom.fpos(2)=-y16;atom.fpos(3)=z16;                     // atom B50
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x16");atom.fpos_equation.push_back("-y16");atom.fpos_equation.push_back("z16");// atom B50 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B50 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B50

    atom.name="C"; atom.type=2;                                       // atom B51
    atom.fpos(1)=-y16;atom.fpos(2)=x16;atom.fpos(3)=z16;                     // atom B51
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y16");atom.fpos_equation.push_back("x16");atom.fpos_equation.push_back("z16");// atom B51 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B51 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B51

    atom.name="C"; atom.type=2;                                       // atom B52
    atom.fpos(1)=y16;atom.fpos(2)=-x16;atom.fpos(3)=z16;                     // atom B52
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y16");atom.fpos_equation.push_back("-x16");atom.fpos_equation.push_back("z16");// atom B52 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B52 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B52

    atom.name="C"; atom.type=2;                                       // atom B53
    atom.fpos(1)=x17;atom.fpos(2)=y17;atom.fpos(3)=z17;                     // atom B53
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x17");atom.fpos_equation.push_back("y17");atom.fpos_equation.push_back("z17");// atom B53 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B53 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B53

    atom.name="C"; atom.type=2;                                       // atom B54
    atom.fpos(1)=-x17;atom.fpos(2)=-y17;atom.fpos(3)=z17;                     // atom B54
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x17");atom.fpos_equation.push_back("-y17");atom.fpos_equation.push_back("z17");// atom B54 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B54 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B54

    atom.name="C"; atom.type=2;                                       // atom B55
    atom.fpos(1)=-y17;atom.fpos(2)=x17;atom.fpos(3)=z17;                     // atom B55
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y17");atom.fpos_equation.push_back("x17");atom.fpos_equation.push_back("z17");// atom B55 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B55 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B55

    atom.name="C"; atom.type=2;                                       // atom B56
    atom.fpos(1)=y17;atom.fpos(2)=-x17;atom.fpos(3)=z17;                     // atom B56
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y17");atom.fpos_equation.push_back("-x17");atom.fpos_equation.push_back("z17");// atom B56 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B56 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B56

    atom.name="C"; atom.type=2;                                       // atom B57
    atom.fpos(1)=x18;atom.fpos(2)=y18;atom.fpos(3)=z18;                     // atom B57
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x18");atom.fpos_equation.push_back("y18");atom.fpos_equation.push_back("z18");// atom B57 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B57 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B57

    atom.name="C"; atom.type=2;                                       // atom B58
    atom.fpos(1)=-x18;atom.fpos(2)=-y18;atom.fpos(3)=z18;                     // atom B58
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x18");atom.fpos_equation.push_back("-y18");atom.fpos_equation.push_back("z18");// atom B58 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B58 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B58

    atom.name="C"; atom.type=2;                                       // atom B59
    atom.fpos(1)=-y18;atom.fpos(2)=x18;atom.fpos(3)=z18;                     // atom B59
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y18");atom.fpos_equation.push_back("x18");atom.fpos_equation.push_back("z18");// atom B59 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B59 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B59

    atom.name="C"; atom.type=2;                                       // atom B60
    atom.fpos(1)=y18;atom.fpos(2)=-x18;atom.fpos(3)=z18;                     // atom B60
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y18");atom.fpos_equation.push_back("-x18");atom.fpos_equation.push_back("z18");// atom B60 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B60 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B60

    atom.name="D"; atom.type=3;                                       // atom B61
    atom.fpos(1)=x19;atom.fpos(2)=y19;atom.fpos(3)=z19;                     // atom B61
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x19");atom.fpos_equation.push_back("y19");atom.fpos_equation.push_back("z19");// atom B61 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B61 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B61

    atom.name="D"; atom.type=3;                                       // atom B62
    atom.fpos(1)=-x19;atom.fpos(2)=-y19;atom.fpos(3)=z19;                     // atom B62
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x19");atom.fpos_equation.push_back("-y19");atom.fpos_equation.push_back("z19");// atom B62 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B62 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B62

    atom.name="D"; atom.type=3;                                       // atom B63
    atom.fpos(1)=-y19;atom.fpos(2)=x19;atom.fpos(3)=z19;                     // atom B63
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y19");atom.fpos_equation.push_back("x19");atom.fpos_equation.push_back("z19");// atom B63 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B63 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B63

    atom.name="D"; atom.type=3;                                       // atom B64
    atom.fpos(1)=y19;atom.fpos(2)=-x19;atom.fpos(3)=z19;                     // atom B64
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y19");atom.fpos_equation.push_back("-x19");atom.fpos_equation.push_back("z19");// atom B64 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B64 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B64

    atom.name="D"; atom.type=3;                                       // atom B65
    atom.fpos(1)=x20;atom.fpos(2)=y20;atom.fpos(3)=z20;                     // atom B65
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x20");atom.fpos_equation.push_back("y20");atom.fpos_equation.push_back("z20");// atom B65 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B65 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B65

    atom.name="D"; atom.type=3;                                       // atom B66
    atom.fpos(1)=-x20;atom.fpos(2)=-y20;atom.fpos(3)=z20;                     // atom B66
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x20");atom.fpos_equation.push_back("-y20");atom.fpos_equation.push_back("z20");// atom B66 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B66 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B66

    atom.name="D"; atom.type=3;                                       // atom B67
    atom.fpos(1)=-y20;atom.fpos(2)=x20;atom.fpos(3)=z20;                     // atom B67
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y20");atom.fpos_equation.push_back("x20");atom.fpos_equation.push_back("z20");// atom B67 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B67 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B67

    atom.name="D"; atom.type=3;                                       // atom B68
    atom.fpos(1)=y20;atom.fpos(2)=-x20;atom.fpos(3)=z20;                     // atom B68
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y20");atom.fpos_equation.push_back("-x20");atom.fpos_equation.push_back("z20");// atom B68 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B68 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B68

    atom.name="D"; atom.type=3;                                       // atom B69
    atom.fpos(1)=x21;atom.fpos(2)=y21;atom.fpos(3)=z21;                     // atom B69
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x21");atom.fpos_equation.push_back("y21");atom.fpos_equation.push_back("z21");// atom B69 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B69 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B69

    atom.name="D"; atom.type=3;                                       // atom B70
    atom.fpos(1)=-x21;atom.fpos(2)=-y21;atom.fpos(3)=z21;                     // atom B70
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x21");atom.fpos_equation.push_back("-y21");atom.fpos_equation.push_back("z21");// atom B70 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B70 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B70

    atom.name="D"; atom.type=3;                                       // atom B71
    atom.fpos(1)=-y21;atom.fpos(2)=x21;atom.fpos(3)=z21;                     // atom B71
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y21");atom.fpos_equation.push_back("x21");atom.fpos_equation.push_back("z21");// atom B71 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B71 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B71

    atom.name="D"; atom.type=3;                                       // atom B72
    atom.fpos(1)=y21;atom.fpos(2)=-x21;atom.fpos(3)=z21;                     // atom B72
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y21");atom.fpos_equation.push_back("-x21");atom.fpos_equation.push_back("z21");// atom B72 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B72 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B72

    atom.name="D"; atom.type=3;                                       // atom B73
    atom.fpos(1)=x22;atom.fpos(2)=y22;atom.fpos(3)=z22;                     // atom B73
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x22");atom.fpos_equation.push_back("y22");atom.fpos_equation.push_back("z22");// atom B73 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B73 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B73

    atom.name="D"; atom.type=3;                                       // atom B74
    atom.fpos(1)=-x22;atom.fpos(2)=-y22;atom.fpos(3)=z22;                     // atom B74
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x22");atom.fpos_equation.push_back("-y22");atom.fpos_equation.push_back("z22");// atom B74 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B74 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B74

    atom.name="D"; atom.type=3;                                       // atom B75
    atom.fpos(1)=-y22;atom.fpos(2)=x22;atom.fpos(3)=z22;                     // atom B75
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y22");atom.fpos_equation.push_back("x22");atom.fpos_equation.push_back("z22");// atom B75 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B75 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B75

    atom.name="D"; atom.type=3;                                       // atom B76
    atom.fpos(1)=y22;atom.fpos(2)=-x22;atom.fpos(3)=z22;                     // atom B76
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y22");atom.fpos_equation.push_back("-x22");atom.fpos_equation.push_back("z22");// atom B76 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B76 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B76


    return str.atoms.size();  
  }
} // namespace anrl

namespace anrl {
  uint WebANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d(stringstream& web,bool LDEBUG) {
#ifndef _ANRL_NOWEB_
#endif

    if(LDEBUG) {cerr << "anrl:: WebANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d: web.str().size()=" << web.str().size() << endl;}

    return web.str().size();
  }
} // namespace anrl

#endif // AFLOW_REMOVE_GREP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

