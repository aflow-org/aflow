// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - David Hicks - 2018
// FILE "ANRL/aflow_anrl_ABC4_tI96_142_e_ab_2g.cpp"

#ifndef _AFLOW_ANRL_ABC4_tI96_142_e_ab_2g_CPP // AFLOW_REMOVE_GREP
#define _AFLOW_ANRL_ABC4_tI96_142_e_ab_2g_CPP // AFLOW_REMOVE_GREP
#include "../aflow.h" // AFLOW_REMOVE_GREP

namespace anrl {
  uint WebANRL_ABC4_tI96_142_e_ab_2g(stringstream &web,bool LDEBUG);
}

namespace anrl {
  uint PrototypeANRL_ABC4_tI96_142_e_ab_2g(stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG) {
    // system ABC4_tI96_142_e_ab_2g

    if(XHOST.vflag_control.flag("WWW")) {
      WebANRL_ABC4_tI96_142_e_ab_2g(web,LDEBUG); // PLUG WEB STUFF
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

    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: FOUND" << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: label=" << label << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: nspecies=" << nspecies << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: natoms=" << natoms << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: spacegroup=" << spacegroup << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: nunderscores=" << nunderscores << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: nparameters=" <<  nparameters << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: Pearson_symbol=" << Pearson_symbol << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: params=" << params << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: Strukturbericht=" << Strukturbericht << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: prototype=" << prototype << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: dialect=" << dialect << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: vparameters.size()=" << vparameters.size() << endl;}

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
    double a=vparameters.at(i++);                  if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: a=" << a << endl;}
    double covera=vparameters.at(i++),c=covera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: c=" << c << " (c/a=" << covera << ")" << endl;}

    double x3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: x3=" << x3 << endl;}
    double x4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: x4=" << x4 << endl;}
    double y4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: y4=" << y4 << endl;}
    double z4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: z4=" << z4 << endl;}
    double x5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: x5=" << x5 << endl;}
    double y5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: y5=" << y5 << endl;}
    double z5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC4_tI96_142_e_ab_2g: z5=" << z5 << endl;}

    str.iomode=IOVASP_AUTO;
    str.title=label+" params="+parameters+" SG="+aurostd::utype2string(spacegroup)+DOI_ANRL; //CO20190520
    str.scale=1.0;

    a1=-(1.0/2.0)*a*xn+(1.0/2.0)*a*yn+(1.0/2.0)*c*zn;
    a2=(1.0/2.0)*a*xn-(1.0/2.0)*a*yn+(1.0/2.0)*c*zn;
    a3=(1.0/2.0)*a*xn+(1.0/2.0)*a*yn-(1.0/2.0)*c*zn;

    str.lattice(1,1)=a1(1);str.lattice(1,2)=a1(2);str.lattice(1,3)=a1(3);
    str.lattice(2,1)=a2(1);str.lattice(2,2)=a2(2);str.lattice(2,3)=a2(3);
    str.lattice(3,1)=a3(1);str.lattice(3,2)=a3(2);str.lattice(3,3)=a3(3);

    // symbolic representation of lattice vectors
    vector<string> a1_equation, a2_equation, a3_equation;
    a1_equation.push_back("-(1.0/2.0)*a");a1_equation.push_back("(1.0/2.0)*a");a1_equation.push_back("(1.0/2.0)*c");
    a2_equation.push_back("(1.0/2.0)*a");a2_equation.push_back("-(1.0/2.0)*a");a2_equation.push_back("(1.0/2.0)*c");
    a3_equation.push_back("(1.0/2.0)*a");a3_equation.push_back("(1.0/2.0)*a");a3_equation.push_back("-(1.0/2.0)*c");
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

    atom.name="A"; atom.type=0;                                       // atom B9
    atom.fpos(1)=(1.0/4.0);atom.fpos(2)=((1.0/4.0)+x3);atom.fpos(3)=x3;                     // atom B9
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("((1.0/4.0)+x3)");atom.fpos_equation.push_back("x3");// atom B9 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B9 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B9

    atom.name="A"; atom.type=0;                                       // atom B10
    atom.fpos(1)=(3.0/4.0);atom.fpos(2)=((1.0/4.0)-x3);atom.fpos(3)=((1.0/2.0)-x3);                     // atom B10
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("((1.0/4.0)-x3)");atom.fpos_equation.push_back("((1.0/2.0)-x3)");// atom B10 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B10 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B10

    atom.name="A"; atom.type=0;                                       // atom B11
    atom.fpos(1)=((1.0/4.0)+x3);atom.fpos(2)=(3.0/4.0);atom.fpos(3)=x3;                     // atom B11
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/4.0)+x3)");atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("x3");// atom B11 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B11 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B11

    atom.name="A"; atom.type=0;                                       // atom B12
    atom.fpos(1)=((1.0/4.0)-x3);atom.fpos(2)=(1.0/4.0);atom.fpos(3)=((1.0/2.0)-x3);                     // atom B12
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/4.0)-x3)");atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("((1.0/2.0)-x3)");// atom B12 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B12 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B12

    atom.name="A"; atom.type=0;                                       // atom B13
    atom.fpos(1)=(3.0/4.0);atom.fpos(2)=((3.0/4.0)-x3);atom.fpos(3)=-x3;                     // atom B13
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("((3.0/4.0)-x3)");atom.fpos_equation.push_back("-x3");// atom B13 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B13 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B13

    atom.name="A"; atom.type=0;                                       // atom B14
    atom.fpos(1)=(1.0/4.0);atom.fpos(2)=((3.0/4.0)+x3);atom.fpos(3)=((1.0/2.0)+x3);                     // atom B14
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("((3.0/4.0)+x3)");atom.fpos_equation.push_back("((1.0/2.0)+x3)");// atom B14 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B14 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B14

    atom.name="A"; atom.type=0;                                       // atom B15
    atom.fpos(1)=((3.0/4.0)-x3);atom.fpos(2)=(1.0/4.0);atom.fpos(3)=-x3;                     // atom B15
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((3.0/4.0)-x3)");atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("-x3");// atom B15 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B15 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B15

    atom.name="A"; atom.type=0;                                       // atom B16
    atom.fpos(1)=((3.0/4.0)+x3);atom.fpos(2)=(3.0/4.0);atom.fpos(3)=((1.0/2.0)+x3);                     // atom B16
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((3.0/4.0)+x3)");atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("((1.0/2.0)+x3)");// atom B16 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B16 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B16

    atom.name="B"; atom.type=1;                                       // atom B1
    atom.fpos(1)=(5.0/8.0);atom.fpos(2)=(3.0/8.0);atom.fpos(3)=(1.0/4.0);                     // atom B1
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(5.0/8.0)");atom.fpos_equation.push_back("(3.0/8.0)");atom.fpos_equation.push_back("(1.0/4.0)");// atom B1 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B1 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B1

    atom.name="B"; atom.type=1;                                       // atom B2
    atom.fpos(1)=(3.0/8.0);atom.fpos(2)=(5.0/8.0);atom.fpos(3)=(3.0/4.0);                     // atom B2
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(3.0/8.0)");atom.fpos_equation.push_back("(5.0/8.0)");atom.fpos_equation.push_back("(3.0/4.0)");// atom B2 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B2 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B2

    atom.name="B"; atom.type=1;                                       // atom B3
    atom.fpos(1)=(7.0/8.0);atom.fpos(2)=(1.0/8.0);atom.fpos(3)=(3.0/4.0);                     // atom B3
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(7.0/8.0)");atom.fpos_equation.push_back("(1.0/8.0)");atom.fpos_equation.push_back("(3.0/4.0)");// atom B3 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B3 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B3

    atom.name="B"; atom.type=1;                                       // atom B4
    atom.fpos(1)=(1.0/8.0);atom.fpos(2)=(7.0/8.0);atom.fpos(3)=(1.0/4.0);                     // atom B4
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/8.0)");atom.fpos_equation.push_back("(7.0/8.0)");atom.fpos_equation.push_back("(1.0/4.0)");// atom B4 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B4 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B4

    atom.name="B"; atom.type=1;                                       // atom B5
    atom.fpos(1)=(3.0/8.0);atom.fpos(2)=(1.0/8.0);atom.fpos(3)=(1.0/4.0);                     // atom B5
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(3.0/8.0)");atom.fpos_equation.push_back("(1.0/8.0)");atom.fpos_equation.push_back("(1.0/4.0)");// atom B5 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B5 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B5

    atom.name="B"; atom.type=1;                                       // atom B6
    atom.fpos(1)=(1.0/8.0);atom.fpos(2)=(3.0/8.0);atom.fpos(3)=(3.0/4.0);                     // atom B6
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/8.0)");atom.fpos_equation.push_back("(3.0/8.0)");atom.fpos_equation.push_back("(3.0/4.0)");// atom B6 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B6 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B6

    atom.name="B"; atom.type=1;                                       // atom B7
    atom.fpos(1)=(5.0/8.0);atom.fpos(2)=(7.0/8.0);atom.fpos(3)=(3.0/4.0);                     // atom B7
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(5.0/8.0)");atom.fpos_equation.push_back("(7.0/8.0)");atom.fpos_equation.push_back("(3.0/4.0)");// atom B7 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B7 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B7

    atom.name="B"; atom.type=1;                                       // atom B8
    atom.fpos(1)=(7.0/8.0);atom.fpos(2)=(5.0/8.0);atom.fpos(3)=(1.0/4.0);                     // atom B8
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(7.0/8.0)");atom.fpos_equation.push_back("(5.0/8.0)");atom.fpos_equation.push_back("(1.0/4.0)");// atom B8 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B8 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B8

    atom.name="C"; atom.type=2;                                       // atom B17
    atom.fpos(1)=(y4+z4);atom.fpos(2)=(x4+z4);atom.fpos(3)=(x4+y4);                     // atom B17
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(y4+z4)");atom.fpos_equation.push_back("(x4+z4)");atom.fpos_equation.push_back("(x4+y4)");// atom B17 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B17 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B17

    atom.name="C"; atom.type=2;                                       // atom B18
    atom.fpos(1)=((1.0/2.0)-y4+z4);atom.fpos(2)=(z4-x4);atom.fpos(3)=((1.0/2.0)-x4-y4);                     // atom B18
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y4+z4)");atom.fpos_equation.push_back("(z4-x4)");atom.fpos_equation.push_back("((1.0/2.0)-x4-y4)");// atom B18 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B18 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B18

    atom.name="C"; atom.type=2;                                       // atom B19
    atom.fpos(1)=(x4+z4);atom.fpos(2)=((1.0/2.0)-y4+z4);atom.fpos(3)=(x4-y4);                     // atom B19
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x4+z4)");atom.fpos_equation.push_back("((1.0/2.0)-y4+z4)");atom.fpos_equation.push_back("(x4-y4)");// atom B19 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B19 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B19

    atom.name="C"; atom.type=2;                                       // atom B20
    atom.fpos(1)=(z4-x4);atom.fpos(2)=(y4+z4);atom.fpos(3)=((1.0/2.0)-x4+y4);                     // atom B20
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(z4-x4)");atom.fpos_equation.push_back("(y4+z4)");atom.fpos_equation.push_back("((1.0/2.0)-x4+y4)");// atom B20 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B20 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B20

    atom.name="C"; atom.type=2;                                       // atom B21
    atom.fpos(1)=(y4-z4);atom.fpos(2)=((1.0/2.0)-x4-z4);atom.fpos(3)=((1.0/2.0)-x4+y4);                     // atom B21
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(y4-z4)");atom.fpos_equation.push_back("((1.0/2.0)-x4-z4)");atom.fpos_equation.push_back("((1.0/2.0)-x4+y4)");// atom B21 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B21 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B21

    atom.name="C"; atom.type=2;                                       // atom B22
    atom.fpos(1)=((1.0/2.0)-y4-z4);atom.fpos(2)=((1.0/2.0)+x4-z4);atom.fpos(3)=(x4-y4);                     // atom B22
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y4-z4)");atom.fpos_equation.push_back("((1.0/2.0)+x4-z4)");atom.fpos_equation.push_back("(x4-y4)");// atom B22 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B22 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B22

    atom.name="C"; atom.type=2;                                       // atom B23
    atom.fpos(1)=((1.0/2.0)+x4-z4);atom.fpos(2)=(y4-z4);atom.fpos(3)=(x4+y4);                     // atom B23
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x4-z4)");atom.fpos_equation.push_back("(y4-z4)");atom.fpos_equation.push_back("(x4+y4)");// atom B23 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B23 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B23

    atom.name="C"; atom.type=2;                                       // atom B24
    atom.fpos(1)=((1.0/2.0)-x4-z4);atom.fpos(2)=((1.0/2.0)-y4-z4);atom.fpos(3)=((1.0/2.0)-x4-y4);                     // atom B24
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x4-z4)");atom.fpos_equation.push_back("((1.0/2.0)-y4-z4)");atom.fpos_equation.push_back("((1.0/2.0)-x4-y4)");// atom B24 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B24 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B24

    atom.name="C"; atom.type=2;                                       // atom B25
    atom.fpos(1)=-(y4+z4);atom.fpos(2)=-(x4+z4);atom.fpos(3)=-(x4+y4);                     // atom B25
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-(y4+z4)");atom.fpos_equation.push_back("-(x4+z4)");atom.fpos_equation.push_back("-(x4+y4)");// atom B25 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B25 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B25

    atom.name="C"; atom.type=2;                                       // atom B26
    atom.fpos(1)=((1.0/2.0)+y4-z4);atom.fpos(2)=(x4-z4);atom.fpos(3)=((1.0/2.0)+x4+y4);                     // atom B26
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y4-z4)");atom.fpos_equation.push_back("(x4-z4)");atom.fpos_equation.push_back("((1.0/2.0)+x4+y4)");// atom B26 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B26 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B26

    atom.name="C"; atom.type=2;                                       // atom B27
    atom.fpos(1)=-(x4+z4);atom.fpos(2)=((1.0/2.0)+y4-z4);atom.fpos(3)=(y4-x4);                     // atom B27
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-(x4+z4)");atom.fpos_equation.push_back("((1.0/2.0)+y4-z4)");atom.fpos_equation.push_back("(y4-x4)");// atom B27 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B27 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B27

    atom.name="C"; atom.type=2;                                       // atom B28
    atom.fpos(1)=(x4-z4);atom.fpos(2)=-(y4+z4);atom.fpos(3)=((1.0/2.0)+x4-y4);                     // atom B28
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x4-z4)");atom.fpos_equation.push_back("-(y4+z4)");atom.fpos_equation.push_back("((1.0/2.0)+x4-y4)");// atom B28 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B28 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B28

    atom.name="C"; atom.type=2;                                       // atom B29
    atom.fpos(1)=(z4-y4);atom.fpos(2)=((1.0/2.0)+x4+z4);atom.fpos(3)=((1.0/2.0)+x4-y4);                     // atom B29
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(z4-y4)");atom.fpos_equation.push_back("((1.0/2.0)+x4+z4)");atom.fpos_equation.push_back("((1.0/2.0)+x4-y4)");// atom B29 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B29 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B29

    atom.name="C"; atom.type=2;                                       // atom B30
    atom.fpos(1)=((1.0/2.0)+y4+z4);atom.fpos(2)=((1.0/2.0)-x4+z4);atom.fpos(3)=(y4-x4);                     // atom B30
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y4+z4)");atom.fpos_equation.push_back("((1.0/2.0)-x4+z4)");atom.fpos_equation.push_back("(y4-x4)");// atom B30 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B30 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B30

    atom.name="C"; atom.type=2;                                       // atom B31
    atom.fpos(1)=((1.0/2.0)-x4+z4);atom.fpos(2)=(z4-y4);atom.fpos(3)=-(x4+y4);                     // atom B31
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x4+z4)");atom.fpos_equation.push_back("(z4-y4)");atom.fpos_equation.push_back("-(x4+y4)");// atom B31 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B31 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B31

    atom.name="C"; atom.type=2;                                       // atom B32
    atom.fpos(1)=((1.0/2.0)+x4+z4);atom.fpos(2)=((1.0/2.0)+y4+z4);atom.fpos(3)=((1.0/2.0)+x4+y4);                     // atom B32
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x4+z4)");atom.fpos_equation.push_back("((1.0/2.0)+y4+z4)");atom.fpos_equation.push_back("((1.0/2.0)+x4+y4)");// atom B32 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B32 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B32

    atom.name="C"; atom.type=2;                                       // atom B33
    atom.fpos(1)=(y5+z5);atom.fpos(2)=(x5+z5);atom.fpos(3)=(x5+y5);                     // atom B33
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(y5+z5)");atom.fpos_equation.push_back("(x5+z5)");atom.fpos_equation.push_back("(x5+y5)");// atom B33 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B33 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B33

    atom.name="C"; atom.type=2;                                       // atom B34
    atom.fpos(1)=((1.0/2.0)-y5+z5);atom.fpos(2)=(z5-x5);atom.fpos(3)=((1.0/2.0)-x5-y5);                     // atom B34
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y5+z5)");atom.fpos_equation.push_back("(z5-x5)");atom.fpos_equation.push_back("((1.0/2.0)-x5-y5)");// atom B34 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B34 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B34

    atom.name="C"; atom.type=2;                                       // atom B35
    atom.fpos(1)=(x5+z5);atom.fpos(2)=((1.0/2.0)-y5+z5);atom.fpos(3)=(x5-y5);                     // atom B35
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x5+z5)");atom.fpos_equation.push_back("((1.0/2.0)-y5+z5)");atom.fpos_equation.push_back("(x5-y5)");// atom B35 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B35 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B35

    atom.name="C"; atom.type=2;                                       // atom B36
    atom.fpos(1)=(z5-x5);atom.fpos(2)=(y5+z5);atom.fpos(3)=((1.0/2.0)-x5+y5);                     // atom B36
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(z5-x5)");atom.fpos_equation.push_back("(y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5+y5)");// atom B36 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B36 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B36

    atom.name="C"; atom.type=2;                                       // atom B37
    atom.fpos(1)=(y5-z5);atom.fpos(2)=((1.0/2.0)-x5-z5);atom.fpos(3)=((1.0/2.0)-x5+y5);                     // atom B37
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5-z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5+y5)");// atom B37 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B37 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B37

    atom.name="C"; atom.type=2;                                       // atom B38
    atom.fpos(1)=((1.0/2.0)-y5-z5);atom.fpos(2)=((1.0/2.0)+x5-z5);atom.fpos(3)=(x5-y5);                     // atom B38
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5-z5)");atom.fpos_equation.push_back("(x5-y5)");// atom B38 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B38 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B38

    atom.name="C"; atom.type=2;                                       // atom B39
    atom.fpos(1)=((1.0/2.0)+x5-z5);atom.fpos(2)=(y5-z5);atom.fpos(3)=(x5+y5);                     // atom B39
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x5-z5)");atom.fpos_equation.push_back("(y5-z5)");atom.fpos_equation.push_back("(x5+y5)");// atom B39 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B39 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B39

    atom.name="C"; atom.type=2;                                       // atom B40
    atom.fpos(1)=((1.0/2.0)-x5-z5);atom.fpos(2)=((1.0/2.0)-y5-z5);atom.fpos(3)=((1.0/2.0)-x5-y5);                     // atom B40
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x5-z5)");atom.fpos_equation.push_back("((1.0/2.0)-y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5-y5)");// atom B40 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B40 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B40

    atom.name="C"; atom.type=2;                                       // atom B41
    atom.fpos(1)=-(y5+z5);atom.fpos(2)=-(x5+z5);atom.fpos(3)=-(x5+y5);                     // atom B41
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-(y5+z5)");atom.fpos_equation.push_back("-(x5+z5)");atom.fpos_equation.push_back("-(x5+y5)");// atom B41 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B41 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B41

    atom.name="C"; atom.type=2;                                       // atom B42
    atom.fpos(1)=((1.0/2.0)+y5-z5);atom.fpos(2)=(x5-z5);atom.fpos(3)=((1.0/2.0)+x5+y5);                     // atom B42
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y5-z5)");atom.fpos_equation.push_back("(x5-z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5+y5)");// atom B42 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B42 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B42

    atom.name="C"; atom.type=2;                                       // atom B43
    atom.fpos(1)=-(x5+z5);atom.fpos(2)=((1.0/2.0)+y5-z5);atom.fpos(3)=(y5-x5);                     // atom B43
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-(x5+z5)");atom.fpos_equation.push_back("((1.0/2.0)+y5-z5)");atom.fpos_equation.push_back("(y5-x5)");// atom B43 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B43 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B43

    atom.name="C"; atom.type=2;                                       // atom B44
    atom.fpos(1)=(x5-z5);atom.fpos(2)=-(y5+z5);atom.fpos(3)=((1.0/2.0)+x5-y5);                     // atom B44
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x5-z5)");atom.fpos_equation.push_back("-(y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5-y5)");// atom B44 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B44 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B44

    atom.name="C"; atom.type=2;                                       // atom B45
    atom.fpos(1)=(z5-y5);atom.fpos(2)=((1.0/2.0)+x5+z5);atom.fpos(3)=((1.0/2.0)+x5-y5);                     // atom B45
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(z5-y5)");atom.fpos_equation.push_back("((1.0/2.0)+x5+z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5-y5)");// atom B45 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B45 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B45

    atom.name="C"; atom.type=2;                                       // atom B46
    atom.fpos(1)=((1.0/2.0)+y5+z5);atom.fpos(2)=((1.0/2.0)-x5+z5);atom.fpos(3)=(y5-x5);                     // atom B46
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5+z5)");atom.fpos_equation.push_back("(y5-x5)");// atom B46 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B46 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B46

    atom.name="C"; atom.type=2;                                       // atom B47
    atom.fpos(1)=((1.0/2.0)-x5+z5);atom.fpos(2)=(z5-y5);atom.fpos(3)=-(x5+y5);                     // atom B47
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x5+z5)");atom.fpos_equation.push_back("(z5-y5)");atom.fpos_equation.push_back("-(x5+y5)");// atom B47 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B47 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B47

    atom.name="C"; atom.type=2;                                       // atom B48
    atom.fpos(1)=((1.0/2.0)+x5+z5);atom.fpos(2)=((1.0/2.0)+y5+z5);atom.fpos(3)=((1.0/2.0)+x5+y5);                     // atom B48
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x5+z5)");atom.fpos_equation.push_back("((1.0/2.0)+y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5+y5)");// atom B48 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B48 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B48


    return str.atoms.size();  
  }
} // namespace anrl

namespace anrl {
  uint WebANRL_ABC4_tI96_142_e_ab_2g(stringstream& web,bool LDEBUG) {
#ifndef _ANRL_NOWEB_
#endif

    if(LDEBUG) {cerr << "anrl:: WebANRL_ABC4_tI96_142_e_ab_2g: web.str().size()=" << web.str().size() << endl;}

    return web.str().size();
  }
} // namespace anrl

#endif // AFLOW_REMOVE_GREP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

