// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - David Hicks - 2020
// FILE "ANRL/aflow_anrl_ABC4_mC48_12_gi_hi_2i3j.cpp"

#ifndef _AFLOW_ANRL_ABC4_mC48_12_gi_hi_2i3j_CPP // AFLOW_REMOVE_GREP
#define _AFLOW_ANRL_ABC4_mC48_12_gi_hi_2i3j_CPP // AFLOW_REMOVE_GREP
#include "../aflow.h" // AFLOW_REMOVE_GREP

namespace anrl {
  uint WebANRL_ABC4_mC48_12_gi_hi_2i3j(stringstream &web,bool LDEBUG);
}

namespace anrl {
  uint PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j(stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG) {
    // system ABC4_mC48_12_gi_hi_2i3j

    if(XHOST.vflag_control.flag("WWW")) {
      WebANRL_ABC4_mC48_12_gi_hi_2i3j(web,LDEBUG); // PLUG WEB STUFF
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

    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: FOUND" << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: label=" << label << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: nspecies=" << nspecies << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: natoms=" << natoms << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: spacegroup=" << spacegroup << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: nunderscores=" << nunderscores << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: nparameters=" <<  nparameters << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: Pearson_symbol=" << Pearson_symbol << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: params=" << params << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: Strukturbericht=" << Strukturbericht << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: prototype=" << prototype << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: dialect=" << dialect << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: vparameters.size()=" << vparameters.size() << endl;}

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
    double a=vparameters.at(i++);                  if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: a=" << a << endl;}
    double bovera=vparameters.at(i++),b=bovera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: b=" << b << " (b/a=" << bovera << ")" << endl;}
    double covera=vparameters.at(i++),c=covera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: c=" << c << " (c/a=" << covera << ")" << endl;}
    double beta=vparameters.at(i++);               if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: beta=" << beta << endl;}

    double y1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: y1=" << y1 << endl;}
    double y2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: y2=" << y2 << endl;}
    double x3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: x3=" << x3 << endl;}
    double z3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: z3=" << z3 << endl;}
    double x4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: x4=" << x4 << endl;}
    double z4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: z4=" << z4 << endl;}
    double x5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: x5=" << x5 << endl;}
    double z5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: z5=" << z5 << endl;}
    double x6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: x6=" << x6 << endl;}
    double z6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: z6=" << z6 << endl;}
    double x7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: x7=" << x7 << endl;}
    double y7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: y7=" << y7 << endl;}
    double z7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: z7=" << z7 << endl;}
    double x8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: x8=" << x8 << endl;}
    double y8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: y8=" << y8 << endl;}
    double z8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: z8=" << z8 << endl;}
    double x9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: x9=" << x9 << endl;}
    double y9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: y9=" << y9 << endl;}
    double z9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: z9=" << z9 << endl;}

    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: cos(beta)=" << cos(deg2rad*beta)  << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_mC48_12_gi_hi_2i3j: sin(beta)=" << sin(deg2rad*beta)  << endl;}

    str.iomode=IOVASP_AUTO;
    str.title=label+" params="+parameters+" SG="+aurostd::utype2string(spacegroup)+DOI_ANRL;
    str.scale=1.0;

    a1=(1.0/2.0)*a*xn-(1.0/2.0)*b*yn;
    a2=(1.0/2.0)*a*xn+(1.0/2.0)*b*yn;
    a3=c*cos(deg2rad*beta)*xn+c*sin(deg2rad*beta)*zn;

    str.lattice(1,1)=a1(1);str.lattice(1,2)=a1(2);str.lattice(1,3)=a1(3);
    str.lattice(2,1)=a2(1);str.lattice(2,2)=a2(2);str.lattice(2,3)=a2(3);
    str.lattice(3,1)=a3(1);str.lattice(3,2)=a3(2);str.lattice(3,3)=a3(3);

    // symbolic representation of lattice vectors
    vector<string> a1_equation, a2_equation, a3_equation;
    a1_equation.push_back("(1.0/2.0)*a");a1_equation.push_back("-(1.0/2.0)*b");a1_equation.push_back("0");
    a2_equation.push_back("(1.0/2.0)*a");a2_equation.push_back("(1.0/2.0)*b");a2_equation.push_back("0");
    a3_equation.push_back("c*cos(beta)");a3_equation.push_back("0");a3_equation.push_back("c*sin(beta)");
    str.symbolic_math_lattice.push_back(a1_equation);
    str.symbolic_math_lattice.push_back(a2_equation);
    str.symbolic_math_lattice.push_back(a3_equation);

    str.num_lattice_parameters = 4;

    str.num_parameters = vparameters.size();
    vector<string> parameter_list; aurostd::string2tokens(params,parameter_list,",");
    str.prototype_parameter_list = parameter_list;
    str.prototype_parameter_values = vparameters;

    if(print_mode!=1){
      str.FixLattices(); // Reciprocal/f2c/c2f
    }

    _atom atom;


    atom.name="A"; atom.type=0;                                       // atom B1
    atom.fpos(1)=-y1;atom.fpos(2)=y1;atom.fpos(3)=0.0;                     // atom B1
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y1");atom.fpos_equation.push_back("y1");atom.fpos_equation.push_back("0.0");// atom B1 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B1 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B1

    atom.name="A"; atom.type=0;                                       // atom B2
    atom.fpos(1)=y1;atom.fpos(2)=-y1;atom.fpos(3)=0.0;                     // atom B2
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y1");atom.fpos_equation.push_back("-y1");atom.fpos_equation.push_back("0.0");// atom B2 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B2 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B2

    atom.name="A"; atom.type=0;                                       // atom B5
    atom.fpos(1)=x3;atom.fpos(2)=x3;atom.fpos(3)=z3;                     // atom B5
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("z3");// atom B5 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B5 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B5

    atom.name="A"; atom.type=0;                                       // atom B6
    atom.fpos(1)=-x3;atom.fpos(2)=-x3;atom.fpos(3)=-z3;                     // atom B6
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("-z3");// atom B6 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B6 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B6

    atom.name="B"; atom.type=1;                                       // atom B3
    atom.fpos(1)=-y2;atom.fpos(2)=y2;atom.fpos(3)=(1.0/2.0);                     // atom B3
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y2");atom.fpos_equation.push_back("y2");atom.fpos_equation.push_back("(1.0/2.0)");// atom B3 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B3 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B3

    atom.name="B"; atom.type=1;                                       // atom B4
    atom.fpos(1)=y2;atom.fpos(2)=-y2;atom.fpos(3)=(1.0/2.0);                     // atom B4
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y2");atom.fpos_equation.push_back("-y2");atom.fpos_equation.push_back("(1.0/2.0)");// atom B4 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B4 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B4

    atom.name="B"; atom.type=1;                                       // atom B7
    atom.fpos(1)=x4;atom.fpos(2)=x4;atom.fpos(3)=z4;                     // atom B7
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x4");atom.fpos_equation.push_back("x4");atom.fpos_equation.push_back("z4");// atom B7 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B7 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B7

    atom.name="B"; atom.type=1;                                       // atom B8
    atom.fpos(1)=-x4;atom.fpos(2)=-x4;atom.fpos(3)=-z4;                     // atom B8
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x4");atom.fpos_equation.push_back("-x4");atom.fpos_equation.push_back("-z4");// atom B8 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B8 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B8

    atom.name="C"; atom.type=2;                                       // atom B9
    atom.fpos(1)=x5;atom.fpos(2)=x5;atom.fpos(3)=z5;                     // atom B9
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("z5");// atom B9 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B9 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B9

    atom.name="C"; atom.type=2;                                       // atom B10
    atom.fpos(1)=-x5;atom.fpos(2)=-x5;atom.fpos(3)=-z5;                     // atom B10
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("-z5");// atom B10 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B10 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B10

    atom.name="C"; atom.type=2;                                       // atom B11
    atom.fpos(1)=x6;atom.fpos(2)=x6;atom.fpos(3)=z6;                     // atom B11
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("z6");// atom B11 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B11 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B11

    atom.name="C"; atom.type=2;                                       // atom B12
    atom.fpos(1)=-x6;atom.fpos(2)=-x6;atom.fpos(3)=-z6;                     // atom B12
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("-z6");// atom B12 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B12 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B12

    atom.name="C"; atom.type=2;                                       // atom B13
    atom.fpos(1)=(x7-y7);atom.fpos(2)=(x7+y7);atom.fpos(3)=z7;                     // atom B13
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x7-y7)");atom.fpos_equation.push_back("(x7+y7)");atom.fpos_equation.push_back("z7");// atom B13 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B13 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B13

    atom.name="C"; atom.type=2;                                       // atom B14
    atom.fpos(1)=(-x7-y7);atom.fpos(2)=(-x7+y7);atom.fpos(3)=-z7;                     // atom B14
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x7-y7)");atom.fpos_equation.push_back("(-x7+y7)");atom.fpos_equation.push_back("-z7");// atom B14 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B14 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B14

    atom.name="C"; atom.type=2;                                       // atom B15
    atom.fpos(1)=(-x7+y7);atom.fpos(2)=(-x7-y7);atom.fpos(3)=-z7;                     // atom B15
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x7+y7)");atom.fpos_equation.push_back("(-x7-y7)");atom.fpos_equation.push_back("-z7");// atom B15 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B15 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B15

    atom.name="C"; atom.type=2;                                       // atom B16
    atom.fpos(1)=(x7+y7);atom.fpos(2)=(x7-y7);atom.fpos(3)=z7;                     // atom B16
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x7+y7)");atom.fpos_equation.push_back("(x7-y7)");atom.fpos_equation.push_back("z7");// atom B16 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B16 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B16

    atom.name="C"; atom.type=2;                                       // atom B17
    atom.fpos(1)=(x8-y8);atom.fpos(2)=(x8+y8);atom.fpos(3)=z8;                     // atom B17
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x8-y8)");atom.fpos_equation.push_back("(x8+y8)");atom.fpos_equation.push_back("z8");// atom B17 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B17 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B17

    atom.name="C"; atom.type=2;                                       // atom B18
    atom.fpos(1)=(-x8-y8);atom.fpos(2)=(-x8+y8);atom.fpos(3)=-z8;                     // atom B18
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x8-y8)");atom.fpos_equation.push_back("(-x8+y8)");atom.fpos_equation.push_back("-z8");// atom B18 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B18 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B18

    atom.name="C"; atom.type=2;                                       // atom B19
    atom.fpos(1)=(-x8+y8);atom.fpos(2)=(-x8-y8);atom.fpos(3)=-z8;                     // atom B19
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x8+y8)");atom.fpos_equation.push_back("(-x8-y8)");atom.fpos_equation.push_back("-z8");// atom B19 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B19 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B19

    atom.name="C"; atom.type=2;                                       // atom B20
    atom.fpos(1)=(x8+y8);atom.fpos(2)=(x8-y8);atom.fpos(3)=z8;                     // atom B20
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x8+y8)");atom.fpos_equation.push_back("(x8-y8)");atom.fpos_equation.push_back("z8");// atom B20 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B20 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B20

    atom.name="C"; atom.type=2;                                       // atom B21
    atom.fpos(1)=(x9-y9);atom.fpos(2)=(x9+y9);atom.fpos(3)=z9;                     // atom B21
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x9-y9)");atom.fpos_equation.push_back("(x9+y9)");atom.fpos_equation.push_back("z9");// atom B21 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B21 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B21

    atom.name="C"; atom.type=2;                                       // atom B22
    atom.fpos(1)=(-x9-y9);atom.fpos(2)=(-x9+y9);atom.fpos(3)=-z9;                     // atom B22
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x9-y9)");atom.fpos_equation.push_back("(-x9+y9)");atom.fpos_equation.push_back("-z9");// atom B22 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B22 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B22

    atom.name="C"; atom.type=2;                                       // atom B23
    atom.fpos(1)=(-x9+y9);atom.fpos(2)=(-x9-y9);atom.fpos(3)=-z9;                     // atom B23
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x9+y9)");atom.fpos_equation.push_back("(-x9-y9)");atom.fpos_equation.push_back("-z9");// atom B23 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B23 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B23

    atom.name="C"; atom.type=2;                                       // atom B24
    atom.fpos(1)=(x9+y9);atom.fpos(2)=(x9-y9);atom.fpos(3)=z9;                     // atom B24
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x9+y9)");atom.fpos_equation.push_back("(x9-y9)");atom.fpos_equation.push_back("z9");// atom B24 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B24 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B24


    return str.atoms.size();  
  }
} // namespace anrl

namespace anrl {
  uint WebANRL_ABC4_mC48_12_gi_hi_2i3j(stringstream& web,bool LDEBUG) {
#ifndef _ANRL_NOWEB_
#endif

    if(LDEBUG) {cerr << "anrl:: WebANRL_ABC4_mC48_12_gi_hi_2i3j: web.str().size()=" << web.str().size() << endl;}

    return web.str().size();
  }
} // namespace anrl

#endif // AFLOW_REMOVE_GREP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
