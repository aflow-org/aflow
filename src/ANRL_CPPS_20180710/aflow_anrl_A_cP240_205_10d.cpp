// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - David Hicks - 2018
// FILE "ANRL/aflow_anrl_A_cP240_205_10d.cpp"

#ifndef _AFLOW_ANRL_A_cP240_205_10d_CPP // AFLOW_REMOVE_GREP
#define _AFLOW_ANRL_A_cP240_205_10d_CPP // AFLOW_REMOVE_GREP
#include "../aflow.h" // AFLOW_REMOVE_GREP

namespace anrl {
  uint WebANRL_A_cP240_205_10d(stringstream &web,bool LDEBUG);
}

namespace anrl {
  uint PrototypeANRL_A_cP240_205_10d(stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG) {
    // system A_cP240_205_10d

    if(XHOST.vflag_control.flag("WWW")) {
      WebANRL_A_cP240_205_10d(web,LDEBUG); // PLUG WEB STUFF
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

    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: FOUND" << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: label=" << label << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: nspecies=" << nspecies << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: natoms=" << natoms << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: spacegroup=" << spacegroup << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: nunderscores=" << nunderscores << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: nparameters=" <<  nparameters << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: Pearson_symbol=" << Pearson_symbol << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: params=" << params << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: Strukturbericht=" << Strukturbericht << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: prototype=" << prototype << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: dialect=" << dialect << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: vparameters.size()=" << vparameters.size() << endl;}

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
    double a=vparameters.at(i++);                  if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: a=" << a << endl;}

    double x1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: x1=" << x1 << endl;}
    double y1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: y1=" << y1 << endl;}
    double z1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: z1=" << z1 << endl;}
    double x2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: x2=" << x2 << endl;}
    double y2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: y2=" << y2 << endl;}
    double z2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: z2=" << z2 << endl;}
    double x3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: x3=" << x3 << endl;}
    double y3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: y3=" << y3 << endl;}
    double z3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: z3=" << z3 << endl;}
    double x4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: x4=" << x4 << endl;}
    double y4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: y4=" << y4 << endl;}
    double z4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: z4=" << z4 << endl;}
    double x5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: x5=" << x5 << endl;}
    double y5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: y5=" << y5 << endl;}
    double z5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: z5=" << z5 << endl;}
    double x6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: x6=" << x6 << endl;}
    double y6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: y6=" << y6 << endl;}
    double z6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: z6=" << z6 << endl;}
    double x7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: x7=" << x7 << endl;}
    double y7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: y7=" << y7 << endl;}
    double z7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: z7=" << z7 << endl;}
    double x8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: x8=" << x8 << endl;}
    double y8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: y8=" << y8 << endl;}
    double z8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: z8=" << z8 << endl;}
    double x9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: x9=" << x9 << endl;}
    double y9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: y9=" << y9 << endl;}
    double z9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: z9=" << z9 << endl;}
    double x10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: x10=" << x10 << endl;}
    double y10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: y10=" << y10 << endl;}
    double z10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A_cP240_205_10d: z10=" << z10 << endl;}

    str.iomode=IOVASP_AUTO;
    str.title=label+" params="+parameters+" SG="+aurostd::utype2string(spacegroup)+DOI_ANRL; //CO20190520
    str.scale=1.0;

    a1=a*xn;
    a2=a*yn;
    a3=a*zn;

    str.lattice(1,1)=a1(1);str.lattice(1,2)=a1(2);str.lattice(1,3)=a1(3);
    str.lattice(2,1)=a2(1);str.lattice(2,2)=a2(2);str.lattice(2,3)=a2(3);
    str.lattice(3,1)=a3(1);str.lattice(3,2)=a3(2);str.lattice(3,3)=a3(3);

    // symbolic representation of lattice vectors
    vector<string> a1_equation, a2_equation, a3_equation;
    a1_equation.push_back("a");a1_equation.push_back("0");a1_equation.push_back("0");
    a2_equation.push_back("0");a2_equation.push_back("a");a2_equation.push_back("0");
    a3_equation.push_back("0");a3_equation.push_back("0");a3_equation.push_back("a");
    str.symbolic_math_lattice.push_back(a1_equation);
    str.symbolic_math_lattice.push_back(a2_equation);
    str.symbolic_math_lattice.push_back(a3_equation);

    str.num_lattice_parameters = 1;

    str.num_parameters = vparameters.size();
    vector<string> parameter_list; aurostd::string2tokens(params,parameter_list,",");
    str.prototype_parameter_list = parameter_list;
    str.prototype_parameter_values = vparameters;

    if(print_mode!=1){
      str.FixLattices(); // Reciprocal/f2c/c2f
    }

    _atom atom;

    atom.name="A"; atom.type=0;                                       // atom B1
    atom.fpos(1)=x1;atom.fpos(2)=y1;atom.fpos(3)=z1;                     // atom B1
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x1");atom.fpos_equation.push_back("y1");atom.fpos_equation.push_back("z1");// atom B1 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B1 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B1

    atom.name="A"; atom.type=0;                                       // atom B2
    atom.fpos(1)=((1.0/2.0)-x1);atom.fpos(2)=-y1;atom.fpos(3)=((1.0/2.0)+z1);                     // atom B2
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x1)");atom.fpos_equation.push_back("-y1");atom.fpos_equation.push_back("((1.0/2.0)+z1)");// atom B2 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B2 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B2

    atom.name="A"; atom.type=0;                                       // atom B3
    atom.fpos(1)=-x1;atom.fpos(2)=((1.0/2.0)+y1);atom.fpos(3)=((1.0/2.0)-z1);                     // atom B3
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x1");atom.fpos_equation.push_back("((1.0/2.0)+y1)");atom.fpos_equation.push_back("((1.0/2.0)-z1)");// atom B3 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B3 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B3

    atom.name="A"; atom.type=0;                                       // atom B4
    atom.fpos(1)=((1.0/2.0)+x1);atom.fpos(2)=((1.0/2.0)-y1);atom.fpos(3)=-z1;                     // atom B4
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x1)");atom.fpos_equation.push_back("((1.0/2.0)-y1)");atom.fpos_equation.push_back("-z1");// atom B4 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B4 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B4

    atom.name="A"; atom.type=0;                                       // atom B5
    atom.fpos(1)=z1;atom.fpos(2)=x1;atom.fpos(3)=y1;                     // atom B5
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z1");atom.fpos_equation.push_back("x1");atom.fpos_equation.push_back("y1");// atom B5 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B5 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B5

    atom.name="A"; atom.type=0;                                       // atom B6
    atom.fpos(1)=((1.0/2.0)+z1);atom.fpos(2)=((1.0/2.0)-x1);atom.fpos(3)=-y1;                     // atom B6
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z1)");atom.fpos_equation.push_back("((1.0/2.0)-x1)");atom.fpos_equation.push_back("-y1");// atom B6 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B6 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B6

    atom.name="A"; atom.type=0;                                       // atom B7
    atom.fpos(1)=((1.0/2.0)-z1);atom.fpos(2)=-x1;atom.fpos(3)=((1.0/2.0)+y1);                     // atom B7
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z1)");atom.fpos_equation.push_back("-x1");atom.fpos_equation.push_back("((1.0/2.0)+y1)");// atom B7 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B7 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B7

    atom.name="A"; atom.type=0;                                       // atom B8
    atom.fpos(1)=-z1;atom.fpos(2)=((1.0/2.0)+x1);atom.fpos(3)=((1.0/2.0)-y1);                     // atom B8
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z1");atom.fpos_equation.push_back("((1.0/2.0)+x1)");atom.fpos_equation.push_back("((1.0/2.0)-y1)");// atom B8 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B8 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B8

    atom.name="A"; atom.type=0;                                       // atom B9
    atom.fpos(1)=y1;atom.fpos(2)=z1;atom.fpos(3)=x1;                     // atom B9
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y1");atom.fpos_equation.push_back("z1");atom.fpos_equation.push_back("x1");// atom B9 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B9 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B9

    atom.name="A"; atom.type=0;                                       // atom B10
    atom.fpos(1)=-y1;atom.fpos(2)=((1.0/2.0)+z1);atom.fpos(3)=((1.0/2.0)-x1);                     // atom B10
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y1");atom.fpos_equation.push_back("((1.0/2.0)+z1)");atom.fpos_equation.push_back("((1.0/2.0)-x1)");// atom B10 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B10 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B10

    atom.name="A"; atom.type=0;                                       // atom B11
    atom.fpos(1)=((1.0/2.0)+y1);atom.fpos(2)=((1.0/2.0)-z1);atom.fpos(3)=-x1;                     // atom B11
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y1)");atom.fpos_equation.push_back("((1.0/2.0)-z1)");atom.fpos_equation.push_back("-x1");// atom B11 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B11 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B11

    atom.name="A"; atom.type=0;                                       // atom B12
    atom.fpos(1)=((1.0/2.0)-y1);atom.fpos(2)=-z1;atom.fpos(3)=((1.0/2.0)+x1);                     // atom B12
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y1)");atom.fpos_equation.push_back("-z1");atom.fpos_equation.push_back("((1.0/2.0)+x1)");// atom B12 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B12 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B12

    atom.name="A"; atom.type=0;                                       // atom B13
    atom.fpos(1)=-x1;atom.fpos(2)=-y1;atom.fpos(3)=-z1;                     // atom B13
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x1");atom.fpos_equation.push_back("-y1");atom.fpos_equation.push_back("-z1");// atom B13 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B13 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B13

    atom.name="A"; atom.type=0;                                       // atom B14
    atom.fpos(1)=((1.0/2.0)+x1);atom.fpos(2)=y1;atom.fpos(3)=((1.0/2.0)-z1);                     // atom B14
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x1)");atom.fpos_equation.push_back("y1");atom.fpos_equation.push_back("((1.0/2.0)-z1)");// atom B14 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B14 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B14

    atom.name="A"; atom.type=0;                                       // atom B15
    atom.fpos(1)=x1;atom.fpos(2)=((1.0/2.0)-y1);atom.fpos(3)=((1.0/2.0)+z1);                     // atom B15
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x1");atom.fpos_equation.push_back("((1.0/2.0)-y1)");atom.fpos_equation.push_back("((1.0/2.0)+z1)");// atom B15 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B15 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B15

    atom.name="A"; atom.type=0;                                       // atom B16
    atom.fpos(1)=((1.0/2.0)-x1);atom.fpos(2)=((1.0/2.0)+y1);atom.fpos(3)=z1;                     // atom B16
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x1)");atom.fpos_equation.push_back("((1.0/2.0)+y1)");atom.fpos_equation.push_back("z1");// atom B16 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B16 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B16

    atom.name="A"; atom.type=0;                                       // atom B17
    atom.fpos(1)=-z1;atom.fpos(2)=-x1;atom.fpos(3)=-y1;                     // atom B17
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z1");atom.fpos_equation.push_back("-x1");atom.fpos_equation.push_back("-y1");// atom B17 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B17 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B17

    atom.name="A"; atom.type=0;                                       // atom B18
    atom.fpos(1)=((1.0/2.0)-z1);atom.fpos(2)=((1.0/2.0)+x1);atom.fpos(3)=y1;                     // atom B18
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z1)");atom.fpos_equation.push_back("((1.0/2.0)+x1)");atom.fpos_equation.push_back("y1");// atom B18 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B18 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B18

    atom.name="A"; atom.type=0;                                       // atom B19
    atom.fpos(1)=((1.0/2.0)+z1);atom.fpos(2)=x1;atom.fpos(3)=((1.0/2.0)-y1);                     // atom B19
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z1)");atom.fpos_equation.push_back("x1");atom.fpos_equation.push_back("((1.0/2.0)-y1)");// atom B19 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B19 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B19

    atom.name="A"; atom.type=0;                                       // atom B20
    atom.fpos(1)=z1;atom.fpos(2)=((1.0/2.0)-x1);atom.fpos(3)=((1.0/2.0)+y1);                     // atom B20
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z1");atom.fpos_equation.push_back("((1.0/2.0)-x1)");atom.fpos_equation.push_back("((1.0/2.0)+y1)");// atom B20 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B20 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B20

    atom.name="A"; atom.type=0;                                       // atom B21
    atom.fpos(1)=-y1;atom.fpos(2)=-z1;atom.fpos(3)=-x1;                     // atom B21
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y1");atom.fpos_equation.push_back("-z1");atom.fpos_equation.push_back("-x1");// atom B21 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B21 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B21

    atom.name="A"; atom.type=0;                                       // atom B22
    atom.fpos(1)=y1;atom.fpos(2)=((1.0/2.0)-z1);atom.fpos(3)=((1.0/2.0)+x1);                     // atom B22
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y1");atom.fpos_equation.push_back("((1.0/2.0)-z1)");atom.fpos_equation.push_back("((1.0/2.0)+x1)");// atom B22 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B22 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B22

    atom.name="A"; atom.type=0;                                       // atom B23
    atom.fpos(1)=((1.0/2.0)-y1);atom.fpos(2)=((1.0/2.0)+z1);atom.fpos(3)=x1;                     // atom B23
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y1)");atom.fpos_equation.push_back("((1.0/2.0)+z1)");atom.fpos_equation.push_back("x1");// atom B23 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B23 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B23

    atom.name="A"; atom.type=0;                                       // atom B24
    atom.fpos(1)=((1.0/2.0)+y1);atom.fpos(2)=z1;atom.fpos(3)=((1.0/2.0)-x1);                     // atom B24
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y1)");atom.fpos_equation.push_back("z1");atom.fpos_equation.push_back("((1.0/2.0)-x1)");// atom B24 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B24 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B24

    atom.name="A"; atom.type=0;                                       // atom B25
    atom.fpos(1)=x2;atom.fpos(2)=y2;atom.fpos(3)=z2;                     // atom B25
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x2");atom.fpos_equation.push_back("y2");atom.fpos_equation.push_back("z2");// atom B25 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B25 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B25

    atom.name="A"; atom.type=0;                                       // atom B26
    atom.fpos(1)=((1.0/2.0)-x2);atom.fpos(2)=-y2;atom.fpos(3)=((1.0/2.0)+z2);                     // atom B26
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x2)");atom.fpos_equation.push_back("-y2");atom.fpos_equation.push_back("((1.0/2.0)+z2)");// atom B26 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B26 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B26

    atom.name="A"; atom.type=0;                                       // atom B27
    atom.fpos(1)=-x2;atom.fpos(2)=((1.0/2.0)+y2);atom.fpos(3)=((1.0/2.0)-z2);                     // atom B27
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x2");atom.fpos_equation.push_back("((1.0/2.0)+y2)");atom.fpos_equation.push_back("((1.0/2.0)-z2)");// atom B27 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B27 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B27

    atom.name="A"; atom.type=0;                                       // atom B28
    atom.fpos(1)=((1.0/2.0)+x2);atom.fpos(2)=((1.0/2.0)-y2);atom.fpos(3)=-z2;                     // atom B28
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x2)");atom.fpos_equation.push_back("((1.0/2.0)-y2)");atom.fpos_equation.push_back("-z2");// atom B28 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B28 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B28

    atom.name="A"; atom.type=0;                                       // atom B29
    atom.fpos(1)=z2;atom.fpos(2)=x2;atom.fpos(3)=y2;                     // atom B29
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z2");atom.fpos_equation.push_back("x2");atom.fpos_equation.push_back("y2");// atom B29 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B29 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B29

    atom.name="A"; atom.type=0;                                       // atom B30
    atom.fpos(1)=((1.0/2.0)+z2);atom.fpos(2)=((1.0/2.0)-x2);atom.fpos(3)=-y2;                     // atom B30
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z2)");atom.fpos_equation.push_back("((1.0/2.0)-x2)");atom.fpos_equation.push_back("-y2");// atom B30 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B30 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B30

    atom.name="A"; atom.type=0;                                       // atom B31
    atom.fpos(1)=((1.0/2.0)-z2);atom.fpos(2)=-x2;atom.fpos(3)=((1.0/2.0)+y2);                     // atom B31
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z2)");atom.fpos_equation.push_back("-x2");atom.fpos_equation.push_back("((1.0/2.0)+y2)");// atom B31 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B31 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B31

    atom.name="A"; atom.type=0;                                       // atom B32
    atom.fpos(1)=-z2;atom.fpos(2)=((1.0/2.0)+x2);atom.fpos(3)=((1.0/2.0)-y2);                     // atom B32
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z2");atom.fpos_equation.push_back("((1.0/2.0)+x2)");atom.fpos_equation.push_back("((1.0/2.0)-y2)");// atom B32 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B32 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B32

    atom.name="A"; atom.type=0;                                       // atom B33
    atom.fpos(1)=y2;atom.fpos(2)=z2;atom.fpos(3)=x2;                     // atom B33
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y2");atom.fpos_equation.push_back("z2");atom.fpos_equation.push_back("x2");// atom B33 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B33 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B33

    atom.name="A"; atom.type=0;                                       // atom B34
    atom.fpos(1)=-y2;atom.fpos(2)=((1.0/2.0)+z2);atom.fpos(3)=((1.0/2.0)-x2);                     // atom B34
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y2");atom.fpos_equation.push_back("((1.0/2.0)+z2)");atom.fpos_equation.push_back("((1.0/2.0)-x2)");// atom B34 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B34 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B34

    atom.name="A"; atom.type=0;                                       // atom B35
    atom.fpos(1)=((1.0/2.0)+y2);atom.fpos(2)=((1.0/2.0)-z2);atom.fpos(3)=-x2;                     // atom B35
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y2)");atom.fpos_equation.push_back("((1.0/2.0)-z2)");atom.fpos_equation.push_back("-x2");// atom B35 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B35 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B35

    atom.name="A"; atom.type=0;                                       // atom B36
    atom.fpos(1)=((1.0/2.0)-y2);atom.fpos(2)=-z2;atom.fpos(3)=((1.0/2.0)+x2);                     // atom B36
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y2)");atom.fpos_equation.push_back("-z2");atom.fpos_equation.push_back("((1.0/2.0)+x2)");// atom B36 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B36 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B36

    atom.name="A"; atom.type=0;                                       // atom B37
    atom.fpos(1)=-x2;atom.fpos(2)=-y2;atom.fpos(3)=-z2;                     // atom B37
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x2");atom.fpos_equation.push_back("-y2");atom.fpos_equation.push_back("-z2");// atom B37 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B37 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B37

    atom.name="A"; atom.type=0;                                       // atom B38
    atom.fpos(1)=((1.0/2.0)+x2);atom.fpos(2)=y2;atom.fpos(3)=((1.0/2.0)-z2);                     // atom B38
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x2)");atom.fpos_equation.push_back("y2");atom.fpos_equation.push_back("((1.0/2.0)-z2)");// atom B38 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B38 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B38

    atom.name="A"; atom.type=0;                                       // atom B39
    atom.fpos(1)=x2;atom.fpos(2)=((1.0/2.0)-y2);atom.fpos(3)=((1.0/2.0)+z2);                     // atom B39
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x2");atom.fpos_equation.push_back("((1.0/2.0)-y2)");atom.fpos_equation.push_back("((1.0/2.0)+z2)");// atom B39 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B39 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B39

    atom.name="A"; atom.type=0;                                       // atom B40
    atom.fpos(1)=((1.0/2.0)-x2);atom.fpos(2)=((1.0/2.0)+y2);atom.fpos(3)=z2;                     // atom B40
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x2)");atom.fpos_equation.push_back("((1.0/2.0)+y2)");atom.fpos_equation.push_back("z2");// atom B40 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B40 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B40

    atom.name="A"; atom.type=0;                                       // atom B41
    atom.fpos(1)=-z2;atom.fpos(2)=-x2;atom.fpos(3)=-y2;                     // atom B41
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z2");atom.fpos_equation.push_back("-x2");atom.fpos_equation.push_back("-y2");// atom B41 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B41 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B41

    atom.name="A"; atom.type=0;                                       // atom B42
    atom.fpos(1)=((1.0/2.0)-z2);atom.fpos(2)=((1.0/2.0)+x2);atom.fpos(3)=y2;                     // atom B42
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z2)");atom.fpos_equation.push_back("((1.0/2.0)+x2)");atom.fpos_equation.push_back("y2");// atom B42 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B42 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B42

    atom.name="A"; atom.type=0;                                       // atom B43
    atom.fpos(1)=((1.0/2.0)+z2);atom.fpos(2)=x2;atom.fpos(3)=((1.0/2.0)-y2);                     // atom B43
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z2)");atom.fpos_equation.push_back("x2");atom.fpos_equation.push_back("((1.0/2.0)-y2)");// atom B43 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B43 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B43

    atom.name="A"; atom.type=0;                                       // atom B44
    atom.fpos(1)=z2;atom.fpos(2)=((1.0/2.0)-x2);atom.fpos(3)=((1.0/2.0)+y2);                     // atom B44
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z2");atom.fpos_equation.push_back("((1.0/2.0)-x2)");atom.fpos_equation.push_back("((1.0/2.0)+y2)");// atom B44 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B44 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B44

    atom.name="A"; atom.type=0;                                       // atom B45
    atom.fpos(1)=-y2;atom.fpos(2)=-z2;atom.fpos(3)=-x2;                     // atom B45
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y2");atom.fpos_equation.push_back("-z2");atom.fpos_equation.push_back("-x2");// atom B45 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B45 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B45

    atom.name="A"; atom.type=0;                                       // atom B46
    atom.fpos(1)=y2;atom.fpos(2)=((1.0/2.0)-z2);atom.fpos(3)=((1.0/2.0)+x2);                     // atom B46
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y2");atom.fpos_equation.push_back("((1.0/2.0)-z2)");atom.fpos_equation.push_back("((1.0/2.0)+x2)");// atom B46 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B46 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B46

    atom.name="A"; atom.type=0;                                       // atom B47
    atom.fpos(1)=((1.0/2.0)-y2);atom.fpos(2)=((1.0/2.0)+z2);atom.fpos(3)=x2;                     // atom B47
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y2)");atom.fpos_equation.push_back("((1.0/2.0)+z2)");atom.fpos_equation.push_back("x2");// atom B47 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B47 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B47

    atom.name="A"; atom.type=0;                                       // atom B48
    atom.fpos(1)=((1.0/2.0)+y2);atom.fpos(2)=z2;atom.fpos(3)=((1.0/2.0)-x2);                     // atom B48
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y2)");atom.fpos_equation.push_back("z2");atom.fpos_equation.push_back("((1.0/2.0)-x2)");// atom B48 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B48 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B48

    atom.name="A"; atom.type=0;                                       // atom B49
    atom.fpos(1)=x3;atom.fpos(2)=y3;atom.fpos(3)=z3;                     // atom B49
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("y3");atom.fpos_equation.push_back("z3");// atom B49 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B49 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B49

    atom.name="A"; atom.type=0;                                       // atom B50
    atom.fpos(1)=((1.0/2.0)-x3);atom.fpos(2)=-y3;atom.fpos(3)=((1.0/2.0)+z3);                     // atom B50
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x3)");atom.fpos_equation.push_back("-y3");atom.fpos_equation.push_back("((1.0/2.0)+z3)");// atom B50 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B50 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B50

    atom.name="A"; atom.type=0;                                       // atom B51
    atom.fpos(1)=-x3;atom.fpos(2)=((1.0/2.0)+y3);atom.fpos(3)=((1.0/2.0)-z3);                     // atom B51
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("((1.0/2.0)+y3)");atom.fpos_equation.push_back("((1.0/2.0)-z3)");// atom B51 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B51 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B51

    atom.name="A"; atom.type=0;                                       // atom B52
    atom.fpos(1)=((1.0/2.0)+x3);atom.fpos(2)=((1.0/2.0)-y3);atom.fpos(3)=-z3;                     // atom B52
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x3)");atom.fpos_equation.push_back("((1.0/2.0)-y3)");atom.fpos_equation.push_back("-z3");// atom B52 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B52 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B52

    atom.name="A"; atom.type=0;                                       // atom B53
    atom.fpos(1)=z3;atom.fpos(2)=x3;atom.fpos(3)=y3;                     // atom B53
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z3");atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("y3");// atom B53 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B53 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B53

    atom.name="A"; atom.type=0;                                       // atom B54
    atom.fpos(1)=((1.0/2.0)+z3);atom.fpos(2)=((1.0/2.0)-x3);atom.fpos(3)=-y3;                     // atom B54
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z3)");atom.fpos_equation.push_back("((1.0/2.0)-x3)");atom.fpos_equation.push_back("-y3");// atom B54 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B54 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B54

    atom.name="A"; atom.type=0;                                       // atom B55
    atom.fpos(1)=((1.0/2.0)-z3);atom.fpos(2)=-x3;atom.fpos(3)=((1.0/2.0)+y3);                     // atom B55
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z3)");atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("((1.0/2.0)+y3)");// atom B55 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B55 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B55

    atom.name="A"; atom.type=0;                                       // atom B56
    atom.fpos(1)=-z3;atom.fpos(2)=((1.0/2.0)+x3);atom.fpos(3)=((1.0/2.0)-y3);                     // atom B56
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z3");atom.fpos_equation.push_back("((1.0/2.0)+x3)");atom.fpos_equation.push_back("((1.0/2.0)-y3)");// atom B56 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B56 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B56

    atom.name="A"; atom.type=0;                                       // atom B57
    atom.fpos(1)=y3;atom.fpos(2)=z3;atom.fpos(3)=x3;                     // atom B57
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y3");atom.fpos_equation.push_back("z3");atom.fpos_equation.push_back("x3");// atom B57 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B57 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B57

    atom.name="A"; atom.type=0;                                       // atom B58
    atom.fpos(1)=-y3;atom.fpos(2)=((1.0/2.0)+z3);atom.fpos(3)=((1.0/2.0)-x3);                     // atom B58
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y3");atom.fpos_equation.push_back("((1.0/2.0)+z3)");atom.fpos_equation.push_back("((1.0/2.0)-x3)");// atom B58 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B58 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B58

    atom.name="A"; atom.type=0;                                       // atom B59
    atom.fpos(1)=((1.0/2.0)+y3);atom.fpos(2)=((1.0/2.0)-z3);atom.fpos(3)=-x3;                     // atom B59
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y3)");atom.fpos_equation.push_back("((1.0/2.0)-z3)");atom.fpos_equation.push_back("-x3");// atom B59 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B59 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B59

    atom.name="A"; atom.type=0;                                       // atom B60
    atom.fpos(1)=((1.0/2.0)-y3);atom.fpos(2)=-z3;atom.fpos(3)=((1.0/2.0)+x3);                     // atom B60
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y3)");atom.fpos_equation.push_back("-z3");atom.fpos_equation.push_back("((1.0/2.0)+x3)");// atom B60 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B60 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B60

    atom.name="A"; atom.type=0;                                       // atom B61
    atom.fpos(1)=-x3;atom.fpos(2)=-y3;atom.fpos(3)=-z3;                     // atom B61
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("-y3");atom.fpos_equation.push_back("-z3");// atom B61 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B61 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B61

    atom.name="A"; atom.type=0;                                       // atom B62
    atom.fpos(1)=((1.0/2.0)+x3);atom.fpos(2)=y3;atom.fpos(3)=((1.0/2.0)-z3);                     // atom B62
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x3)");atom.fpos_equation.push_back("y3");atom.fpos_equation.push_back("((1.0/2.0)-z3)");// atom B62 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B62 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B62

    atom.name="A"; atom.type=0;                                       // atom B63
    atom.fpos(1)=x3;atom.fpos(2)=((1.0/2.0)-y3);atom.fpos(3)=((1.0/2.0)+z3);                     // atom B63
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("((1.0/2.0)-y3)");atom.fpos_equation.push_back("((1.0/2.0)+z3)");// atom B63 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B63 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B63

    atom.name="A"; atom.type=0;                                       // atom B64
    atom.fpos(1)=((1.0/2.0)-x3);atom.fpos(2)=((1.0/2.0)+y3);atom.fpos(3)=z3;                     // atom B64
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x3)");atom.fpos_equation.push_back("((1.0/2.0)+y3)");atom.fpos_equation.push_back("z3");// atom B64 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B64 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B64

    atom.name="A"; atom.type=0;                                       // atom B65
    atom.fpos(1)=-z3;atom.fpos(2)=-x3;atom.fpos(3)=-y3;                     // atom B65
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z3");atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("-y3");// atom B65 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B65 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B65

    atom.name="A"; atom.type=0;                                       // atom B66
    atom.fpos(1)=((1.0/2.0)-z3);atom.fpos(2)=((1.0/2.0)+x3);atom.fpos(3)=y3;                     // atom B66
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z3)");atom.fpos_equation.push_back("((1.0/2.0)+x3)");atom.fpos_equation.push_back("y3");// atom B66 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B66 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B66

    atom.name="A"; atom.type=0;                                       // atom B67
    atom.fpos(1)=((1.0/2.0)+z3);atom.fpos(2)=x3;atom.fpos(3)=((1.0/2.0)-y3);                     // atom B67
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z3)");atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("((1.0/2.0)-y3)");// atom B67 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B67 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B67

    atom.name="A"; atom.type=0;                                       // atom B68
    atom.fpos(1)=z3;atom.fpos(2)=((1.0/2.0)-x3);atom.fpos(3)=((1.0/2.0)+y3);                     // atom B68
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z3");atom.fpos_equation.push_back("((1.0/2.0)-x3)");atom.fpos_equation.push_back("((1.0/2.0)+y3)");// atom B68 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B68 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B68

    atom.name="A"; atom.type=0;                                       // atom B69
    atom.fpos(1)=-y3;atom.fpos(2)=-z3;atom.fpos(3)=-x3;                     // atom B69
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y3");atom.fpos_equation.push_back("-z3");atom.fpos_equation.push_back("-x3");// atom B69 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B69 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B69

    atom.name="A"; atom.type=0;                                       // atom B70
    atom.fpos(1)=y3;atom.fpos(2)=((1.0/2.0)-z3);atom.fpos(3)=((1.0/2.0)+x3);                     // atom B70
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y3");atom.fpos_equation.push_back("((1.0/2.0)-z3)");atom.fpos_equation.push_back("((1.0/2.0)+x3)");// atom B70 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B70 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B70

    atom.name="A"; atom.type=0;                                       // atom B71
    atom.fpos(1)=((1.0/2.0)-y3);atom.fpos(2)=((1.0/2.0)+z3);atom.fpos(3)=x3;                     // atom B71
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y3)");atom.fpos_equation.push_back("((1.0/2.0)+z3)");atom.fpos_equation.push_back("x3");// atom B71 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B71 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B71

    atom.name="A"; atom.type=0;                                       // atom B72
    atom.fpos(1)=((1.0/2.0)+y3);atom.fpos(2)=z3;atom.fpos(3)=((1.0/2.0)-x3);                     // atom B72
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y3)");atom.fpos_equation.push_back("z3");atom.fpos_equation.push_back("((1.0/2.0)-x3)");// atom B72 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B72 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B72

    atom.name="A"; atom.type=0;                                       // atom B73
    atom.fpos(1)=x4;atom.fpos(2)=y4;atom.fpos(3)=z4;                     // atom B73
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x4");atom.fpos_equation.push_back("y4");atom.fpos_equation.push_back("z4");// atom B73 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B73 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B73

    atom.name="A"; atom.type=0;                                       // atom B74
    atom.fpos(1)=((1.0/2.0)-x4);atom.fpos(2)=-y4;atom.fpos(3)=((1.0/2.0)+z4);                     // atom B74
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x4)");atom.fpos_equation.push_back("-y4");atom.fpos_equation.push_back("((1.0/2.0)+z4)");// atom B74 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B74 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B74

    atom.name="A"; atom.type=0;                                       // atom B75
    atom.fpos(1)=-x4;atom.fpos(2)=((1.0/2.0)+y4);atom.fpos(3)=((1.0/2.0)-z4);                     // atom B75
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x4");atom.fpos_equation.push_back("((1.0/2.0)+y4)");atom.fpos_equation.push_back("((1.0/2.0)-z4)");// atom B75 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B75 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B75

    atom.name="A"; atom.type=0;                                       // atom B76
    atom.fpos(1)=((1.0/2.0)+x4);atom.fpos(2)=((1.0/2.0)-y4);atom.fpos(3)=-z4;                     // atom B76
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x4)");atom.fpos_equation.push_back("((1.0/2.0)-y4)");atom.fpos_equation.push_back("-z4");// atom B76 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B76 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B76

    atom.name="A"; atom.type=0;                                       // atom B77
    atom.fpos(1)=z4;atom.fpos(2)=x4;atom.fpos(3)=y4;                     // atom B77
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z4");atom.fpos_equation.push_back("x4");atom.fpos_equation.push_back("y4");// atom B77 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B77 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B77

    atom.name="A"; atom.type=0;                                       // atom B78
    atom.fpos(1)=((1.0/2.0)+z4);atom.fpos(2)=((1.0/2.0)-x4);atom.fpos(3)=-y4;                     // atom B78
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z4)");atom.fpos_equation.push_back("((1.0/2.0)-x4)");atom.fpos_equation.push_back("-y4");// atom B78 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B78 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B78

    atom.name="A"; atom.type=0;                                       // atom B79
    atom.fpos(1)=((1.0/2.0)-z4);atom.fpos(2)=-x4;atom.fpos(3)=((1.0/2.0)+y4);                     // atom B79
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z4)");atom.fpos_equation.push_back("-x4");atom.fpos_equation.push_back("((1.0/2.0)+y4)");// atom B79 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B79 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B79

    atom.name="A"; atom.type=0;                                       // atom B80
    atom.fpos(1)=-z4;atom.fpos(2)=((1.0/2.0)+x4);atom.fpos(3)=((1.0/2.0)-y4);                     // atom B80
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z4");atom.fpos_equation.push_back("((1.0/2.0)+x4)");atom.fpos_equation.push_back("((1.0/2.0)-y4)");// atom B80 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B80 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B80

    atom.name="A"; atom.type=0;                                       // atom B81
    atom.fpos(1)=y4;atom.fpos(2)=z4;atom.fpos(3)=x4;                     // atom B81
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y4");atom.fpos_equation.push_back("z4");atom.fpos_equation.push_back("x4");// atom B81 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B81 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B81

    atom.name="A"; atom.type=0;                                       // atom B82
    atom.fpos(1)=-y4;atom.fpos(2)=((1.0/2.0)+z4);atom.fpos(3)=((1.0/2.0)-x4);                     // atom B82
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y4");atom.fpos_equation.push_back("((1.0/2.0)+z4)");atom.fpos_equation.push_back("((1.0/2.0)-x4)");// atom B82 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B82 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B82

    atom.name="A"; atom.type=0;                                       // atom B83
    atom.fpos(1)=((1.0/2.0)+y4);atom.fpos(2)=((1.0/2.0)-z4);atom.fpos(3)=-x4;                     // atom B83
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y4)");atom.fpos_equation.push_back("((1.0/2.0)-z4)");atom.fpos_equation.push_back("-x4");// atom B83 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B83 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B83

    atom.name="A"; atom.type=0;                                       // atom B84
    atom.fpos(1)=((1.0/2.0)-y4);atom.fpos(2)=-z4;atom.fpos(3)=((1.0/2.0)+x4);                     // atom B84
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y4)");atom.fpos_equation.push_back("-z4");atom.fpos_equation.push_back("((1.0/2.0)+x4)");// atom B84 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B84 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B84

    atom.name="A"; atom.type=0;                                       // atom B85
    atom.fpos(1)=-x4;atom.fpos(2)=-y4;atom.fpos(3)=-z4;                     // atom B85
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x4");atom.fpos_equation.push_back("-y4");atom.fpos_equation.push_back("-z4");// atom B85 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B85 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B85

    atom.name="A"; atom.type=0;                                       // atom B86
    atom.fpos(1)=((1.0/2.0)+x4);atom.fpos(2)=y4;atom.fpos(3)=((1.0/2.0)-z4);                     // atom B86
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x4)");atom.fpos_equation.push_back("y4");atom.fpos_equation.push_back("((1.0/2.0)-z4)");// atom B86 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B86 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B86

    atom.name="A"; atom.type=0;                                       // atom B87
    atom.fpos(1)=x4;atom.fpos(2)=((1.0/2.0)-y4);atom.fpos(3)=((1.0/2.0)+z4);                     // atom B87
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x4");atom.fpos_equation.push_back("((1.0/2.0)-y4)");atom.fpos_equation.push_back("((1.0/2.0)+z4)");// atom B87 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B87 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B87

    atom.name="A"; atom.type=0;                                       // atom B88
    atom.fpos(1)=((1.0/2.0)-x4);atom.fpos(2)=((1.0/2.0)+y4);atom.fpos(3)=z4;                     // atom B88
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x4)");atom.fpos_equation.push_back("((1.0/2.0)+y4)");atom.fpos_equation.push_back("z4");// atom B88 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B88 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B88

    atom.name="A"; atom.type=0;                                       // atom B89
    atom.fpos(1)=-z4;atom.fpos(2)=-x4;atom.fpos(3)=-y4;                     // atom B89
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z4");atom.fpos_equation.push_back("-x4");atom.fpos_equation.push_back("-y4");// atom B89 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B89 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B89

    atom.name="A"; atom.type=0;                                       // atom B90
    atom.fpos(1)=((1.0/2.0)-z4);atom.fpos(2)=((1.0/2.0)+x4);atom.fpos(3)=y4;                     // atom B90
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z4)");atom.fpos_equation.push_back("((1.0/2.0)+x4)");atom.fpos_equation.push_back("y4");// atom B90 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B90 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B90

    atom.name="A"; atom.type=0;                                       // atom B91
    atom.fpos(1)=((1.0/2.0)+z4);atom.fpos(2)=x4;atom.fpos(3)=((1.0/2.0)-y4);                     // atom B91
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z4)");atom.fpos_equation.push_back("x4");atom.fpos_equation.push_back("((1.0/2.0)-y4)");// atom B91 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B91 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B91

    atom.name="A"; atom.type=0;                                       // atom B92
    atom.fpos(1)=z4;atom.fpos(2)=((1.0/2.0)-x4);atom.fpos(3)=((1.0/2.0)+y4);                     // atom B92
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z4");atom.fpos_equation.push_back("((1.0/2.0)-x4)");atom.fpos_equation.push_back("((1.0/2.0)+y4)");// atom B92 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B92 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B92

    atom.name="A"; atom.type=0;                                       // atom B93
    atom.fpos(1)=-y4;atom.fpos(2)=-z4;atom.fpos(3)=-x4;                     // atom B93
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y4");atom.fpos_equation.push_back("-z4");atom.fpos_equation.push_back("-x4");// atom B93 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B93 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B93

    atom.name="A"; atom.type=0;                                       // atom B94
    atom.fpos(1)=y4;atom.fpos(2)=((1.0/2.0)-z4);atom.fpos(3)=((1.0/2.0)+x4);                     // atom B94
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y4");atom.fpos_equation.push_back("((1.0/2.0)-z4)");atom.fpos_equation.push_back("((1.0/2.0)+x4)");// atom B94 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B94 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B94

    atom.name="A"; atom.type=0;                                       // atom B95
    atom.fpos(1)=((1.0/2.0)-y4);atom.fpos(2)=((1.0/2.0)+z4);atom.fpos(3)=x4;                     // atom B95
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y4)");atom.fpos_equation.push_back("((1.0/2.0)+z4)");atom.fpos_equation.push_back("x4");// atom B95 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B95 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B95

    atom.name="A"; atom.type=0;                                       // atom B96
    atom.fpos(1)=((1.0/2.0)+y4);atom.fpos(2)=z4;atom.fpos(3)=((1.0/2.0)-x4);                     // atom B96
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y4)");atom.fpos_equation.push_back("z4");atom.fpos_equation.push_back("((1.0/2.0)-x4)");// atom B96 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B96 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B96

    atom.name="A"; atom.type=0;                                       // atom B97
    atom.fpos(1)=x5;atom.fpos(2)=y5;atom.fpos(3)=z5;                     // atom B97
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("y5");atom.fpos_equation.push_back("z5");// atom B97 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B97 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B97

    atom.name="A"; atom.type=0;                                       // atom B98
    atom.fpos(1)=((1.0/2.0)-x5);atom.fpos(2)=-y5;atom.fpos(3)=((1.0/2.0)+z5);                     // atom B98
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x5)");atom.fpos_equation.push_back("-y5");atom.fpos_equation.push_back("((1.0/2.0)+z5)");// atom B98 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B98 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B98

    atom.name="A"; atom.type=0;                                       // atom B99
    atom.fpos(1)=-x5;atom.fpos(2)=((1.0/2.0)+y5);atom.fpos(3)=((1.0/2.0)-z5);                     // atom B99
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("((1.0/2.0)+y5)");atom.fpos_equation.push_back("((1.0/2.0)-z5)");// atom B99 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B99 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B99

    atom.name="A"; atom.type=0;                                       // atom B100
    atom.fpos(1)=((1.0/2.0)+x5);atom.fpos(2)=((1.0/2.0)-y5);atom.fpos(3)=-z5;                     // atom B100
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x5)");atom.fpos_equation.push_back("((1.0/2.0)-y5)");atom.fpos_equation.push_back("-z5");// atom B100 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B100 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B100

    atom.name="A"; atom.type=0;                                       // atom B101
    atom.fpos(1)=z5;atom.fpos(2)=x5;atom.fpos(3)=y5;                     // atom B101
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z5");atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("y5");// atom B101 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B101 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B101

    atom.name="A"; atom.type=0;                                       // atom B102
    atom.fpos(1)=((1.0/2.0)+z5);atom.fpos(2)=((1.0/2.0)-x5);atom.fpos(3)=-y5;                     // atom B102
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5)");atom.fpos_equation.push_back("-y5");// atom B102 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B102 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B102

    atom.name="A"; atom.type=0;                                       // atom B103
    atom.fpos(1)=((1.0/2.0)-z5);atom.fpos(2)=-x5;atom.fpos(3)=((1.0/2.0)+y5);                     // atom B103
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z5)");atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("((1.0/2.0)+y5)");// atom B103 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B103 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B103

    atom.name="A"; atom.type=0;                                       // atom B104
    atom.fpos(1)=-z5;atom.fpos(2)=((1.0/2.0)+x5);atom.fpos(3)=((1.0/2.0)-y5);                     // atom B104
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z5");atom.fpos_equation.push_back("((1.0/2.0)+x5)");atom.fpos_equation.push_back("((1.0/2.0)-y5)");// atom B104 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B104 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B104

    atom.name="A"; atom.type=0;                                       // atom B105
    atom.fpos(1)=y5;atom.fpos(2)=z5;atom.fpos(3)=x5;                     // atom B105
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y5");atom.fpos_equation.push_back("z5");atom.fpos_equation.push_back("x5");// atom B105 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B105 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B105

    atom.name="A"; atom.type=0;                                       // atom B106
    atom.fpos(1)=-y5;atom.fpos(2)=((1.0/2.0)+z5);atom.fpos(3)=((1.0/2.0)-x5);                     // atom B106
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y5");atom.fpos_equation.push_back("((1.0/2.0)+z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5)");// atom B106 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B106 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B106

    atom.name="A"; atom.type=0;                                       // atom B107
    atom.fpos(1)=((1.0/2.0)+y5);atom.fpos(2)=((1.0/2.0)-z5);atom.fpos(3)=-x5;                     // atom B107
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y5)");atom.fpos_equation.push_back("((1.0/2.0)-z5)");atom.fpos_equation.push_back("-x5");// atom B107 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B107 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B107

    atom.name="A"; atom.type=0;                                       // atom B108
    atom.fpos(1)=((1.0/2.0)-y5);atom.fpos(2)=-z5;atom.fpos(3)=((1.0/2.0)+x5);                     // atom B108
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y5)");atom.fpos_equation.push_back("-z5");atom.fpos_equation.push_back("((1.0/2.0)+x5)");// atom B108 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B108 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B108

    atom.name="A"; atom.type=0;                                       // atom B109
    atom.fpos(1)=-x5;atom.fpos(2)=-y5;atom.fpos(3)=-z5;                     // atom B109
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("-y5");atom.fpos_equation.push_back("-z5");// atom B109 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B109 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B109

    atom.name="A"; atom.type=0;                                       // atom B110
    atom.fpos(1)=((1.0/2.0)+x5);atom.fpos(2)=y5;atom.fpos(3)=((1.0/2.0)-z5);                     // atom B110
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x5)");atom.fpos_equation.push_back("y5");atom.fpos_equation.push_back("((1.0/2.0)-z5)");// atom B110 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B110 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B110

    atom.name="A"; atom.type=0;                                       // atom B111
    atom.fpos(1)=x5;atom.fpos(2)=((1.0/2.0)-y5);atom.fpos(3)=((1.0/2.0)+z5);                     // atom B111
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("((1.0/2.0)-y5)");atom.fpos_equation.push_back("((1.0/2.0)+z5)");// atom B111 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B111 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B111

    atom.name="A"; atom.type=0;                                       // atom B112
    atom.fpos(1)=((1.0/2.0)-x5);atom.fpos(2)=((1.0/2.0)+y5);atom.fpos(3)=z5;                     // atom B112
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x5)");atom.fpos_equation.push_back("((1.0/2.0)+y5)");atom.fpos_equation.push_back("z5");// atom B112 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B112 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B112

    atom.name="A"; atom.type=0;                                       // atom B113
    atom.fpos(1)=-z5;atom.fpos(2)=-x5;atom.fpos(3)=-y5;                     // atom B113
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z5");atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("-y5");// atom B113 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B113 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B113

    atom.name="A"; atom.type=0;                                       // atom B114
    atom.fpos(1)=((1.0/2.0)-z5);atom.fpos(2)=((1.0/2.0)+x5);atom.fpos(3)=y5;                     // atom B114
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5)");atom.fpos_equation.push_back("y5");// atom B114 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B114 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B114

    atom.name="A"; atom.type=0;                                       // atom B115
    atom.fpos(1)=((1.0/2.0)+z5);atom.fpos(2)=x5;atom.fpos(3)=((1.0/2.0)-y5);                     // atom B115
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z5)");atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("((1.0/2.0)-y5)");// atom B115 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B115 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B115

    atom.name="A"; atom.type=0;                                       // atom B116
    atom.fpos(1)=z5;atom.fpos(2)=((1.0/2.0)-x5);atom.fpos(3)=((1.0/2.0)+y5);                     // atom B116
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z5");atom.fpos_equation.push_back("((1.0/2.0)-x5)");atom.fpos_equation.push_back("((1.0/2.0)+y5)");// atom B116 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B116 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B116

    atom.name="A"; atom.type=0;                                       // atom B117
    atom.fpos(1)=-y5;atom.fpos(2)=-z5;atom.fpos(3)=-x5;                     // atom B117
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y5");atom.fpos_equation.push_back("-z5");atom.fpos_equation.push_back("-x5");// atom B117 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B117 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B117

    atom.name="A"; atom.type=0;                                       // atom B118
    atom.fpos(1)=y5;atom.fpos(2)=((1.0/2.0)-z5);atom.fpos(3)=((1.0/2.0)+x5);                     // atom B118
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y5");atom.fpos_equation.push_back("((1.0/2.0)-z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5)");// atom B118 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B118 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B118

    atom.name="A"; atom.type=0;                                       // atom B119
    atom.fpos(1)=((1.0/2.0)-y5);atom.fpos(2)=((1.0/2.0)+z5);atom.fpos(3)=x5;                     // atom B119
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y5)");atom.fpos_equation.push_back("((1.0/2.0)+z5)");atom.fpos_equation.push_back("x5");// atom B119 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B119 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B119

    atom.name="A"; atom.type=0;                                       // atom B120
    atom.fpos(1)=((1.0/2.0)+y5);atom.fpos(2)=z5;atom.fpos(3)=((1.0/2.0)-x5);                     // atom B120
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y5)");atom.fpos_equation.push_back("z5");atom.fpos_equation.push_back("((1.0/2.0)-x5)");// atom B120 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B120 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B120

    atom.name="A"; atom.type=0;                                       // atom B121
    atom.fpos(1)=x6;atom.fpos(2)=y6;atom.fpos(3)=z6;                     // atom B121
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("y6");atom.fpos_equation.push_back("z6");// atom B121 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B121 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B121

    atom.name="A"; atom.type=0;                                       // atom B122
    atom.fpos(1)=((1.0/2.0)-x6);atom.fpos(2)=-y6;atom.fpos(3)=((1.0/2.0)+z6);                     // atom B122
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x6)");atom.fpos_equation.push_back("-y6");atom.fpos_equation.push_back("((1.0/2.0)+z6)");// atom B122 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B122 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B122

    atom.name="A"; atom.type=0;                                       // atom B123
    atom.fpos(1)=-x6;atom.fpos(2)=((1.0/2.0)+y6);atom.fpos(3)=((1.0/2.0)-z6);                     // atom B123
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("((1.0/2.0)+y6)");atom.fpos_equation.push_back("((1.0/2.0)-z6)");// atom B123 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B123 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B123

    atom.name="A"; atom.type=0;                                       // atom B124
    atom.fpos(1)=((1.0/2.0)+x6);atom.fpos(2)=((1.0/2.0)-y6);atom.fpos(3)=-z6;                     // atom B124
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x6)");atom.fpos_equation.push_back("((1.0/2.0)-y6)");atom.fpos_equation.push_back("-z6");// atom B124 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B124 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B124

    atom.name="A"; atom.type=0;                                       // atom B125
    atom.fpos(1)=z6;atom.fpos(2)=x6;atom.fpos(3)=y6;                     // atom B125
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z6");atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("y6");// atom B125 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B125 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B125

    atom.name="A"; atom.type=0;                                       // atom B126
    atom.fpos(1)=((1.0/2.0)+z6);atom.fpos(2)=((1.0/2.0)-x6);atom.fpos(3)=-y6;                     // atom B126
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z6)");atom.fpos_equation.push_back("((1.0/2.0)-x6)");atom.fpos_equation.push_back("-y6");// atom B126 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B126 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B126

    atom.name="A"; atom.type=0;                                       // atom B127
    atom.fpos(1)=((1.0/2.0)-z6);atom.fpos(2)=-x6;atom.fpos(3)=((1.0/2.0)+y6);                     // atom B127
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z6)");atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("((1.0/2.0)+y6)");// atom B127 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B127 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B127

    atom.name="A"; atom.type=0;                                       // atom B128
    atom.fpos(1)=-z6;atom.fpos(2)=((1.0/2.0)+x6);atom.fpos(3)=((1.0/2.0)-y6);                     // atom B128
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z6");atom.fpos_equation.push_back("((1.0/2.0)+x6)");atom.fpos_equation.push_back("((1.0/2.0)-y6)");// atom B128 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B128 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B128

    atom.name="A"; atom.type=0;                                       // atom B129
    atom.fpos(1)=y6;atom.fpos(2)=z6;atom.fpos(3)=x6;                     // atom B129
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y6");atom.fpos_equation.push_back("z6");atom.fpos_equation.push_back("x6");// atom B129 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B129 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B129

    atom.name="A"; atom.type=0;                                       // atom B130
    atom.fpos(1)=-y6;atom.fpos(2)=((1.0/2.0)+z6);atom.fpos(3)=((1.0/2.0)-x6);                     // atom B130
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y6");atom.fpos_equation.push_back("((1.0/2.0)+z6)");atom.fpos_equation.push_back("((1.0/2.0)-x6)");// atom B130 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B130 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B130

    atom.name="A"; atom.type=0;                                       // atom B131
    atom.fpos(1)=((1.0/2.0)+y6);atom.fpos(2)=((1.0/2.0)-z6);atom.fpos(3)=-x6;                     // atom B131
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y6)");atom.fpos_equation.push_back("((1.0/2.0)-z6)");atom.fpos_equation.push_back("-x6");// atom B131 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B131 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B131

    atom.name="A"; atom.type=0;                                       // atom B132
    atom.fpos(1)=((1.0/2.0)-y6);atom.fpos(2)=-z6;atom.fpos(3)=((1.0/2.0)+x6);                     // atom B132
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y6)");atom.fpos_equation.push_back("-z6");atom.fpos_equation.push_back("((1.0/2.0)+x6)");// atom B132 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B132 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B132

    atom.name="A"; atom.type=0;                                       // atom B133
    atom.fpos(1)=-x6;atom.fpos(2)=-y6;atom.fpos(3)=-z6;                     // atom B133
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("-y6");atom.fpos_equation.push_back("-z6");// atom B133 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B133 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B133

    atom.name="A"; atom.type=0;                                       // atom B134
    atom.fpos(1)=((1.0/2.0)+x6);atom.fpos(2)=y6;atom.fpos(3)=((1.0/2.0)-z6);                     // atom B134
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x6)");atom.fpos_equation.push_back("y6");atom.fpos_equation.push_back("((1.0/2.0)-z6)");// atom B134 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B134 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B134

    atom.name="A"; atom.type=0;                                       // atom B135
    atom.fpos(1)=x6;atom.fpos(2)=((1.0/2.0)-y6);atom.fpos(3)=((1.0/2.0)+z6);                     // atom B135
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("((1.0/2.0)-y6)");atom.fpos_equation.push_back("((1.0/2.0)+z6)");// atom B135 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B135 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B135

    atom.name="A"; atom.type=0;                                       // atom B136
    atom.fpos(1)=((1.0/2.0)-x6);atom.fpos(2)=((1.0/2.0)+y6);atom.fpos(3)=z6;                     // atom B136
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x6)");atom.fpos_equation.push_back("((1.0/2.0)+y6)");atom.fpos_equation.push_back("z6");// atom B136 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B136 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B136

    atom.name="A"; atom.type=0;                                       // atom B137
    atom.fpos(1)=-z6;atom.fpos(2)=-x6;atom.fpos(3)=-y6;                     // atom B137
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z6");atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("-y6");// atom B137 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B137 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B137

    atom.name="A"; atom.type=0;                                       // atom B138
    atom.fpos(1)=((1.0/2.0)-z6);atom.fpos(2)=((1.0/2.0)+x6);atom.fpos(3)=y6;                     // atom B138
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z6)");atom.fpos_equation.push_back("((1.0/2.0)+x6)");atom.fpos_equation.push_back("y6");// atom B138 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B138 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B138

    atom.name="A"; atom.type=0;                                       // atom B139
    atom.fpos(1)=((1.0/2.0)+z6);atom.fpos(2)=x6;atom.fpos(3)=((1.0/2.0)-y6);                     // atom B139
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z6)");atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("((1.0/2.0)-y6)");// atom B139 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B139 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B139

    atom.name="A"; atom.type=0;                                       // atom B140
    atom.fpos(1)=z6;atom.fpos(2)=((1.0/2.0)-x6);atom.fpos(3)=((1.0/2.0)+y6);                     // atom B140
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z6");atom.fpos_equation.push_back("((1.0/2.0)-x6)");atom.fpos_equation.push_back("((1.0/2.0)+y6)");// atom B140 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B140 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B140

    atom.name="A"; atom.type=0;                                       // atom B141
    atom.fpos(1)=-y6;atom.fpos(2)=-z6;atom.fpos(3)=-x6;                     // atom B141
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y6");atom.fpos_equation.push_back("-z6");atom.fpos_equation.push_back("-x6");// atom B141 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B141 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B141

    atom.name="A"; atom.type=0;                                       // atom B142
    atom.fpos(1)=y6;atom.fpos(2)=((1.0/2.0)-z6);atom.fpos(3)=((1.0/2.0)+x6);                     // atom B142
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y6");atom.fpos_equation.push_back("((1.0/2.0)-z6)");atom.fpos_equation.push_back("((1.0/2.0)+x6)");// atom B142 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B142 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B142

    atom.name="A"; atom.type=0;                                       // atom B143
    atom.fpos(1)=((1.0/2.0)-y6);atom.fpos(2)=((1.0/2.0)+z6);atom.fpos(3)=x6;                     // atom B143
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y6)");atom.fpos_equation.push_back("((1.0/2.0)+z6)");atom.fpos_equation.push_back("x6");// atom B143 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B143 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B143

    atom.name="A"; atom.type=0;                                       // atom B144
    atom.fpos(1)=((1.0/2.0)+y6);atom.fpos(2)=z6;atom.fpos(3)=((1.0/2.0)-x6);                     // atom B144
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y6)");atom.fpos_equation.push_back("z6");atom.fpos_equation.push_back("((1.0/2.0)-x6)");// atom B144 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B144 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B144

    atom.name="A"; atom.type=0;                                       // atom B145
    atom.fpos(1)=x7;atom.fpos(2)=y7;atom.fpos(3)=z7;                     // atom B145
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x7");atom.fpos_equation.push_back("y7");atom.fpos_equation.push_back("z7");// atom B145 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B145 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B145

    atom.name="A"; atom.type=0;                                       // atom B146
    atom.fpos(1)=((1.0/2.0)-x7);atom.fpos(2)=-y7;atom.fpos(3)=((1.0/2.0)+z7);                     // atom B146
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x7)");atom.fpos_equation.push_back("-y7");atom.fpos_equation.push_back("((1.0/2.0)+z7)");// atom B146 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B146 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B146

    atom.name="A"; atom.type=0;                                       // atom B147
    atom.fpos(1)=-x7;atom.fpos(2)=((1.0/2.0)+y7);atom.fpos(3)=((1.0/2.0)-z7);                     // atom B147
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x7");atom.fpos_equation.push_back("((1.0/2.0)+y7)");atom.fpos_equation.push_back("((1.0/2.0)-z7)");// atom B147 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B147 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B147

    atom.name="A"; atom.type=0;                                       // atom B148
    atom.fpos(1)=((1.0/2.0)+x7);atom.fpos(2)=((1.0/2.0)-y7);atom.fpos(3)=-z7;                     // atom B148
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x7)");atom.fpos_equation.push_back("((1.0/2.0)-y7)");atom.fpos_equation.push_back("-z7");// atom B148 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B148 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B148

    atom.name="A"; atom.type=0;                                       // atom B149
    atom.fpos(1)=z7;atom.fpos(2)=x7;atom.fpos(3)=y7;                     // atom B149
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z7");atom.fpos_equation.push_back("x7");atom.fpos_equation.push_back("y7");// atom B149 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B149 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B149

    atom.name="A"; atom.type=0;                                       // atom B150
    atom.fpos(1)=((1.0/2.0)+z7);atom.fpos(2)=((1.0/2.0)-x7);atom.fpos(3)=-y7;                     // atom B150
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z7)");atom.fpos_equation.push_back("((1.0/2.0)-x7)");atom.fpos_equation.push_back("-y7");// atom B150 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B150 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B150

    atom.name="A"; atom.type=0;                                       // atom B151
    atom.fpos(1)=((1.0/2.0)-z7);atom.fpos(2)=-x7;atom.fpos(3)=((1.0/2.0)+y7);                     // atom B151
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z7)");atom.fpos_equation.push_back("-x7");atom.fpos_equation.push_back("((1.0/2.0)+y7)");// atom B151 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B151 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B151

    atom.name="A"; atom.type=0;                                       // atom B152
    atom.fpos(1)=-z7;atom.fpos(2)=((1.0/2.0)+x7);atom.fpos(3)=((1.0/2.0)-y7);                     // atom B152
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z7");atom.fpos_equation.push_back("((1.0/2.0)+x7)");atom.fpos_equation.push_back("((1.0/2.0)-y7)");// atom B152 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B152 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B152

    atom.name="A"; atom.type=0;                                       // atom B153
    atom.fpos(1)=y7;atom.fpos(2)=z7;atom.fpos(3)=x7;                     // atom B153
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y7");atom.fpos_equation.push_back("z7");atom.fpos_equation.push_back("x7");// atom B153 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B153 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B153

    atom.name="A"; atom.type=0;                                       // atom B154
    atom.fpos(1)=-y7;atom.fpos(2)=((1.0/2.0)+z7);atom.fpos(3)=((1.0/2.0)-x7);                     // atom B154
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y7");atom.fpos_equation.push_back("((1.0/2.0)+z7)");atom.fpos_equation.push_back("((1.0/2.0)-x7)");// atom B154 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B154 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B154

    atom.name="A"; atom.type=0;                                       // atom B155
    atom.fpos(1)=((1.0/2.0)+y7);atom.fpos(2)=((1.0/2.0)-z7);atom.fpos(3)=-x7;                     // atom B155
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y7)");atom.fpos_equation.push_back("((1.0/2.0)-z7)");atom.fpos_equation.push_back("-x7");// atom B155 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B155 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B155

    atom.name="A"; atom.type=0;                                       // atom B156
    atom.fpos(1)=((1.0/2.0)-y7);atom.fpos(2)=-z7;atom.fpos(3)=((1.0/2.0)+x7);                     // atom B156
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y7)");atom.fpos_equation.push_back("-z7");atom.fpos_equation.push_back("((1.0/2.0)+x7)");// atom B156 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B156 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B156

    atom.name="A"; atom.type=0;                                       // atom B157
    atom.fpos(1)=-x7;atom.fpos(2)=-y7;atom.fpos(3)=-z7;                     // atom B157
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x7");atom.fpos_equation.push_back("-y7");atom.fpos_equation.push_back("-z7");// atom B157 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B157 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B157

    atom.name="A"; atom.type=0;                                       // atom B158
    atom.fpos(1)=((1.0/2.0)+x7);atom.fpos(2)=y7;atom.fpos(3)=((1.0/2.0)-z7);                     // atom B158
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x7)");atom.fpos_equation.push_back("y7");atom.fpos_equation.push_back("((1.0/2.0)-z7)");// atom B158 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B158 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B158

    atom.name="A"; atom.type=0;                                       // atom B159
    atom.fpos(1)=x7;atom.fpos(2)=((1.0/2.0)-y7);atom.fpos(3)=((1.0/2.0)+z7);                     // atom B159
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x7");atom.fpos_equation.push_back("((1.0/2.0)-y7)");atom.fpos_equation.push_back("((1.0/2.0)+z7)");// atom B159 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B159 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B159

    atom.name="A"; atom.type=0;                                       // atom B160
    atom.fpos(1)=((1.0/2.0)-x7);atom.fpos(2)=((1.0/2.0)+y7);atom.fpos(3)=z7;                     // atom B160
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x7)");atom.fpos_equation.push_back("((1.0/2.0)+y7)");atom.fpos_equation.push_back("z7");// atom B160 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B160 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B160

    atom.name="A"; atom.type=0;                                       // atom B161
    atom.fpos(1)=-z7;atom.fpos(2)=-x7;atom.fpos(3)=-y7;                     // atom B161
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z7");atom.fpos_equation.push_back("-x7");atom.fpos_equation.push_back("-y7");// atom B161 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B161 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B161

    atom.name="A"; atom.type=0;                                       // atom B162
    atom.fpos(1)=((1.0/2.0)-z7);atom.fpos(2)=((1.0/2.0)+x7);atom.fpos(3)=y7;                     // atom B162
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z7)");atom.fpos_equation.push_back("((1.0/2.0)+x7)");atom.fpos_equation.push_back("y7");// atom B162 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B162 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B162

    atom.name="A"; atom.type=0;                                       // atom B163
    atom.fpos(1)=((1.0/2.0)+z7);atom.fpos(2)=x7;atom.fpos(3)=((1.0/2.0)-y7);                     // atom B163
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z7)");atom.fpos_equation.push_back("x7");atom.fpos_equation.push_back("((1.0/2.0)-y7)");// atom B163 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B163 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B163

    atom.name="A"; atom.type=0;                                       // atom B164
    atom.fpos(1)=z7;atom.fpos(2)=((1.0/2.0)-x7);atom.fpos(3)=((1.0/2.0)+y7);                     // atom B164
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z7");atom.fpos_equation.push_back("((1.0/2.0)-x7)");atom.fpos_equation.push_back("((1.0/2.0)+y7)");// atom B164 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B164 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B164

    atom.name="A"; atom.type=0;                                       // atom B165
    atom.fpos(1)=-y7;atom.fpos(2)=-z7;atom.fpos(3)=-x7;                     // atom B165
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y7");atom.fpos_equation.push_back("-z7");atom.fpos_equation.push_back("-x7");// atom B165 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B165 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B165

    atom.name="A"; atom.type=0;                                       // atom B166
    atom.fpos(1)=y7;atom.fpos(2)=((1.0/2.0)-z7);atom.fpos(3)=((1.0/2.0)+x7);                     // atom B166
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y7");atom.fpos_equation.push_back("((1.0/2.0)-z7)");atom.fpos_equation.push_back("((1.0/2.0)+x7)");// atom B166 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B166 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B166

    atom.name="A"; atom.type=0;                                       // atom B167
    atom.fpos(1)=((1.0/2.0)-y7);atom.fpos(2)=((1.0/2.0)+z7);atom.fpos(3)=x7;                     // atom B167
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y7)");atom.fpos_equation.push_back("((1.0/2.0)+z7)");atom.fpos_equation.push_back("x7");// atom B167 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B167 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B167

    atom.name="A"; atom.type=0;                                       // atom B168
    atom.fpos(1)=((1.0/2.0)+y7);atom.fpos(2)=z7;atom.fpos(3)=((1.0/2.0)-x7);                     // atom B168
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y7)");atom.fpos_equation.push_back("z7");atom.fpos_equation.push_back("((1.0/2.0)-x7)");// atom B168 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B168 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B168

    atom.name="A"; atom.type=0;                                       // atom B169
    atom.fpos(1)=x8;atom.fpos(2)=y8;atom.fpos(3)=z8;                     // atom B169
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x8");atom.fpos_equation.push_back("y8");atom.fpos_equation.push_back("z8");// atom B169 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B169 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B169

    atom.name="A"; atom.type=0;                                       // atom B170
    atom.fpos(1)=((1.0/2.0)-x8);atom.fpos(2)=-y8;atom.fpos(3)=((1.0/2.0)+z8);                     // atom B170
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x8)");atom.fpos_equation.push_back("-y8");atom.fpos_equation.push_back("((1.0/2.0)+z8)");// atom B170 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B170 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B170

    atom.name="A"; atom.type=0;                                       // atom B171
    atom.fpos(1)=-x8;atom.fpos(2)=((1.0/2.0)+y8);atom.fpos(3)=((1.0/2.0)-z8);                     // atom B171
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x8");atom.fpos_equation.push_back("((1.0/2.0)+y8)");atom.fpos_equation.push_back("((1.0/2.0)-z8)");// atom B171 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B171 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B171

    atom.name="A"; atom.type=0;                                       // atom B172
    atom.fpos(1)=((1.0/2.0)+x8);atom.fpos(2)=((1.0/2.0)-y8);atom.fpos(3)=-z8;                     // atom B172
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x8)");atom.fpos_equation.push_back("((1.0/2.0)-y8)");atom.fpos_equation.push_back("-z8");// atom B172 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B172 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B172

    atom.name="A"; atom.type=0;                                       // atom B173
    atom.fpos(1)=z8;atom.fpos(2)=x8;atom.fpos(3)=y8;                     // atom B173
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z8");atom.fpos_equation.push_back("x8");atom.fpos_equation.push_back("y8");// atom B173 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B173 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B173

    atom.name="A"; atom.type=0;                                       // atom B174
    atom.fpos(1)=((1.0/2.0)+z8);atom.fpos(2)=((1.0/2.0)-x8);atom.fpos(3)=-y8;                     // atom B174
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z8)");atom.fpos_equation.push_back("((1.0/2.0)-x8)");atom.fpos_equation.push_back("-y8");// atom B174 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B174 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B174

    atom.name="A"; atom.type=0;                                       // atom B175
    atom.fpos(1)=((1.0/2.0)-z8);atom.fpos(2)=-x8;atom.fpos(3)=((1.0/2.0)+y8);                     // atom B175
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z8)");atom.fpos_equation.push_back("-x8");atom.fpos_equation.push_back("((1.0/2.0)+y8)");// atom B175 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B175 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B175

    atom.name="A"; atom.type=0;                                       // atom B176
    atom.fpos(1)=-z8;atom.fpos(2)=((1.0/2.0)+x8);atom.fpos(3)=((1.0/2.0)-y8);                     // atom B176
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z8");atom.fpos_equation.push_back("((1.0/2.0)+x8)");atom.fpos_equation.push_back("((1.0/2.0)-y8)");// atom B176 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B176 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B176

    atom.name="A"; atom.type=0;                                       // atom B177
    atom.fpos(1)=y8;atom.fpos(2)=z8;atom.fpos(3)=x8;                     // atom B177
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y8");atom.fpos_equation.push_back("z8");atom.fpos_equation.push_back("x8");// atom B177 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B177 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B177

    atom.name="A"; atom.type=0;                                       // atom B178
    atom.fpos(1)=-y8;atom.fpos(2)=((1.0/2.0)+z8);atom.fpos(3)=((1.0/2.0)-x8);                     // atom B178
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y8");atom.fpos_equation.push_back("((1.0/2.0)+z8)");atom.fpos_equation.push_back("((1.0/2.0)-x8)");// atom B178 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B178 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B178

    atom.name="A"; atom.type=0;                                       // atom B179
    atom.fpos(1)=((1.0/2.0)+y8);atom.fpos(2)=((1.0/2.0)-z8);atom.fpos(3)=-x8;                     // atom B179
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y8)");atom.fpos_equation.push_back("((1.0/2.0)-z8)");atom.fpos_equation.push_back("-x8");// atom B179 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B179 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B179

    atom.name="A"; atom.type=0;                                       // atom B180
    atom.fpos(1)=((1.0/2.0)-y8);atom.fpos(2)=-z8;atom.fpos(3)=((1.0/2.0)+x8);                     // atom B180
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y8)");atom.fpos_equation.push_back("-z8");atom.fpos_equation.push_back("((1.0/2.0)+x8)");// atom B180 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B180 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B180

    atom.name="A"; atom.type=0;                                       // atom B181
    atom.fpos(1)=-x8;atom.fpos(2)=-y8;atom.fpos(3)=-z8;                     // atom B181
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x8");atom.fpos_equation.push_back("-y8");atom.fpos_equation.push_back("-z8");// atom B181 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B181 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B181

    atom.name="A"; atom.type=0;                                       // atom B182
    atom.fpos(1)=((1.0/2.0)+x8);atom.fpos(2)=y8;atom.fpos(3)=((1.0/2.0)-z8);                     // atom B182
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x8)");atom.fpos_equation.push_back("y8");atom.fpos_equation.push_back("((1.0/2.0)-z8)");// atom B182 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B182 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B182

    atom.name="A"; atom.type=0;                                       // atom B183
    atom.fpos(1)=x8;atom.fpos(2)=((1.0/2.0)-y8);atom.fpos(3)=((1.0/2.0)+z8);                     // atom B183
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x8");atom.fpos_equation.push_back("((1.0/2.0)-y8)");atom.fpos_equation.push_back("((1.0/2.0)+z8)");// atom B183 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B183 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B183

    atom.name="A"; atom.type=0;                                       // atom B184
    atom.fpos(1)=((1.0/2.0)-x8);atom.fpos(2)=((1.0/2.0)+y8);atom.fpos(3)=z8;                     // atom B184
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x8)");atom.fpos_equation.push_back("((1.0/2.0)+y8)");atom.fpos_equation.push_back("z8");// atom B184 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B184 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B184

    atom.name="A"; atom.type=0;                                       // atom B185
    atom.fpos(1)=-z8;atom.fpos(2)=-x8;atom.fpos(3)=-y8;                     // atom B185
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z8");atom.fpos_equation.push_back("-x8");atom.fpos_equation.push_back("-y8");// atom B185 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B185 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B185

    atom.name="A"; atom.type=0;                                       // atom B186
    atom.fpos(1)=((1.0/2.0)-z8);atom.fpos(2)=((1.0/2.0)+x8);atom.fpos(3)=y8;                     // atom B186
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z8)");atom.fpos_equation.push_back("((1.0/2.0)+x8)");atom.fpos_equation.push_back("y8");// atom B186 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B186 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B186

    atom.name="A"; atom.type=0;                                       // atom B187
    atom.fpos(1)=((1.0/2.0)+z8);atom.fpos(2)=x8;atom.fpos(3)=((1.0/2.0)-y8);                     // atom B187
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z8)");atom.fpos_equation.push_back("x8");atom.fpos_equation.push_back("((1.0/2.0)-y8)");// atom B187 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B187 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B187

    atom.name="A"; atom.type=0;                                       // atom B188
    atom.fpos(1)=z8;atom.fpos(2)=((1.0/2.0)-x8);atom.fpos(3)=((1.0/2.0)+y8);                     // atom B188
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z8");atom.fpos_equation.push_back("((1.0/2.0)-x8)");atom.fpos_equation.push_back("((1.0/2.0)+y8)");// atom B188 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B188 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B188

    atom.name="A"; atom.type=0;                                       // atom B189
    atom.fpos(1)=-y8;atom.fpos(2)=-z8;atom.fpos(3)=-x8;                     // atom B189
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y8");atom.fpos_equation.push_back("-z8");atom.fpos_equation.push_back("-x8");// atom B189 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B189 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B189

    atom.name="A"; atom.type=0;                                       // atom B190
    atom.fpos(1)=y8;atom.fpos(2)=((1.0/2.0)-z8);atom.fpos(3)=((1.0/2.0)+x8);                     // atom B190
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y8");atom.fpos_equation.push_back("((1.0/2.0)-z8)");atom.fpos_equation.push_back("((1.0/2.0)+x8)");// atom B190 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B190 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B190

    atom.name="A"; atom.type=0;                                       // atom B191
    atom.fpos(1)=((1.0/2.0)-y8);atom.fpos(2)=((1.0/2.0)+z8);atom.fpos(3)=x8;                     // atom B191
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y8)");atom.fpos_equation.push_back("((1.0/2.0)+z8)");atom.fpos_equation.push_back("x8");// atom B191 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B191 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B191

    atom.name="A"; atom.type=0;                                       // atom B192
    atom.fpos(1)=((1.0/2.0)+y8);atom.fpos(2)=z8;atom.fpos(3)=((1.0/2.0)-x8);                     // atom B192
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y8)");atom.fpos_equation.push_back("z8");atom.fpos_equation.push_back("((1.0/2.0)-x8)");// atom B192 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B192 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B192

    atom.name="A"; atom.type=0;                                       // atom B193
    atom.fpos(1)=x9;atom.fpos(2)=y9;atom.fpos(3)=z9;                     // atom B193
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x9");atom.fpos_equation.push_back("y9");atom.fpos_equation.push_back("z9");// atom B193 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B193 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B193

    atom.name="A"; atom.type=0;                                       // atom B194
    atom.fpos(1)=((1.0/2.0)-x9);atom.fpos(2)=-y9;atom.fpos(3)=((1.0/2.0)+z9);                     // atom B194
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x9)");atom.fpos_equation.push_back("-y9");atom.fpos_equation.push_back("((1.0/2.0)+z9)");// atom B194 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B194 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B194

    atom.name="A"; atom.type=0;                                       // atom B195
    atom.fpos(1)=-x9;atom.fpos(2)=((1.0/2.0)+y9);atom.fpos(3)=((1.0/2.0)-z9);                     // atom B195
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x9");atom.fpos_equation.push_back("((1.0/2.0)+y9)");atom.fpos_equation.push_back("((1.0/2.0)-z9)");// atom B195 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B195 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B195

    atom.name="A"; atom.type=0;                                       // atom B196
    atom.fpos(1)=((1.0/2.0)+x9);atom.fpos(2)=((1.0/2.0)-y9);atom.fpos(3)=-z9;                     // atom B196
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x9)");atom.fpos_equation.push_back("((1.0/2.0)-y9)");atom.fpos_equation.push_back("-z9");// atom B196 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B196 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B196

    atom.name="A"; atom.type=0;                                       // atom B197
    atom.fpos(1)=z9;atom.fpos(2)=x9;atom.fpos(3)=y9;                     // atom B197
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z9");atom.fpos_equation.push_back("x9");atom.fpos_equation.push_back("y9");// atom B197 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B197 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B197

    atom.name="A"; atom.type=0;                                       // atom B198
    atom.fpos(1)=((1.0/2.0)+z9);atom.fpos(2)=((1.0/2.0)-x9);atom.fpos(3)=-y9;                     // atom B198
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z9)");atom.fpos_equation.push_back("((1.0/2.0)-x9)");atom.fpos_equation.push_back("-y9");// atom B198 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B198 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B198

    atom.name="A"; atom.type=0;                                       // atom B199
    atom.fpos(1)=((1.0/2.0)-z9);atom.fpos(2)=-x9;atom.fpos(3)=((1.0/2.0)+y9);                     // atom B199
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z9)");atom.fpos_equation.push_back("-x9");atom.fpos_equation.push_back("((1.0/2.0)+y9)");// atom B199 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B199 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B199

    atom.name="A"; atom.type=0;                                       // atom B200
    atom.fpos(1)=-z9;atom.fpos(2)=((1.0/2.0)+x9);atom.fpos(3)=((1.0/2.0)-y9);                     // atom B200
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z9");atom.fpos_equation.push_back("((1.0/2.0)+x9)");atom.fpos_equation.push_back("((1.0/2.0)-y9)");// atom B200 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B200 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B200

    atom.name="A"; atom.type=0;                                       // atom B201
    atom.fpos(1)=y9;atom.fpos(2)=z9;atom.fpos(3)=x9;                     // atom B201
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y9");atom.fpos_equation.push_back("z9");atom.fpos_equation.push_back("x9");// atom B201 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B201 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B201

    atom.name="A"; atom.type=0;                                       // atom B202
    atom.fpos(1)=-y9;atom.fpos(2)=((1.0/2.0)+z9);atom.fpos(3)=((1.0/2.0)-x9);                     // atom B202
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y9");atom.fpos_equation.push_back("((1.0/2.0)+z9)");atom.fpos_equation.push_back("((1.0/2.0)-x9)");// atom B202 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B202 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B202

    atom.name="A"; atom.type=0;                                       // atom B203
    atom.fpos(1)=((1.0/2.0)+y9);atom.fpos(2)=((1.0/2.0)-z9);atom.fpos(3)=-x9;                     // atom B203
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y9)");atom.fpos_equation.push_back("((1.0/2.0)-z9)");atom.fpos_equation.push_back("-x9");// atom B203 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B203 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B203

    atom.name="A"; atom.type=0;                                       // atom B204
    atom.fpos(1)=((1.0/2.0)-y9);atom.fpos(2)=-z9;atom.fpos(3)=((1.0/2.0)+x9);                     // atom B204
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y9)");atom.fpos_equation.push_back("-z9");atom.fpos_equation.push_back("((1.0/2.0)+x9)");// atom B204 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B204 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B204

    atom.name="A"; atom.type=0;                                       // atom B205
    atom.fpos(1)=-x9;atom.fpos(2)=-y9;atom.fpos(3)=-z9;                     // atom B205
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x9");atom.fpos_equation.push_back("-y9");atom.fpos_equation.push_back("-z9");// atom B205 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B205 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B205

    atom.name="A"; atom.type=0;                                       // atom B206
    atom.fpos(1)=((1.0/2.0)+x9);atom.fpos(2)=y9;atom.fpos(3)=((1.0/2.0)-z9);                     // atom B206
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x9)");atom.fpos_equation.push_back("y9");atom.fpos_equation.push_back("((1.0/2.0)-z9)");// atom B206 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B206 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B206

    atom.name="A"; atom.type=0;                                       // atom B207
    atom.fpos(1)=x9;atom.fpos(2)=((1.0/2.0)-y9);atom.fpos(3)=((1.0/2.0)+z9);                     // atom B207
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x9");atom.fpos_equation.push_back("((1.0/2.0)-y9)");atom.fpos_equation.push_back("((1.0/2.0)+z9)");// atom B207 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B207 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B207

    atom.name="A"; atom.type=0;                                       // atom B208
    atom.fpos(1)=((1.0/2.0)-x9);atom.fpos(2)=((1.0/2.0)+y9);atom.fpos(3)=z9;                     // atom B208
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x9)");atom.fpos_equation.push_back("((1.0/2.0)+y9)");atom.fpos_equation.push_back("z9");// atom B208 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B208 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B208

    atom.name="A"; atom.type=0;                                       // atom B209
    atom.fpos(1)=-z9;atom.fpos(2)=-x9;atom.fpos(3)=-y9;                     // atom B209
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z9");atom.fpos_equation.push_back("-x9");atom.fpos_equation.push_back("-y9");// atom B209 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B209 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B209

    atom.name="A"; atom.type=0;                                       // atom B210
    atom.fpos(1)=((1.0/2.0)-z9);atom.fpos(2)=((1.0/2.0)+x9);atom.fpos(3)=y9;                     // atom B210
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z9)");atom.fpos_equation.push_back("((1.0/2.0)+x9)");atom.fpos_equation.push_back("y9");// atom B210 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B210 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B210

    atom.name="A"; atom.type=0;                                       // atom B211
    atom.fpos(1)=((1.0/2.0)+z9);atom.fpos(2)=x9;atom.fpos(3)=((1.0/2.0)-y9);                     // atom B211
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z9)");atom.fpos_equation.push_back("x9");atom.fpos_equation.push_back("((1.0/2.0)-y9)");// atom B211 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B211 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B211

    atom.name="A"; atom.type=0;                                       // atom B212
    atom.fpos(1)=z9;atom.fpos(2)=((1.0/2.0)-x9);atom.fpos(3)=((1.0/2.0)+y9);                     // atom B212
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z9");atom.fpos_equation.push_back("((1.0/2.0)-x9)");atom.fpos_equation.push_back("((1.0/2.0)+y9)");// atom B212 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B212 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B212

    atom.name="A"; atom.type=0;                                       // atom B213
    atom.fpos(1)=-y9;atom.fpos(2)=-z9;atom.fpos(3)=-x9;                     // atom B213
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y9");atom.fpos_equation.push_back("-z9");atom.fpos_equation.push_back("-x9");// atom B213 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B213 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B213

    atom.name="A"; atom.type=0;                                       // atom B214
    atom.fpos(1)=y9;atom.fpos(2)=((1.0/2.0)-z9);atom.fpos(3)=((1.0/2.0)+x9);                     // atom B214
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y9");atom.fpos_equation.push_back("((1.0/2.0)-z9)");atom.fpos_equation.push_back("((1.0/2.0)+x9)");// atom B214 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B214 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B214

    atom.name="A"; atom.type=0;                                       // atom B215
    atom.fpos(1)=((1.0/2.0)-y9);atom.fpos(2)=((1.0/2.0)+z9);atom.fpos(3)=x9;                     // atom B215
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y9)");atom.fpos_equation.push_back("((1.0/2.0)+z9)");atom.fpos_equation.push_back("x9");// atom B215 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B215 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B215

    atom.name="A"; atom.type=0;                                       // atom B216
    atom.fpos(1)=((1.0/2.0)+y9);atom.fpos(2)=z9;atom.fpos(3)=((1.0/2.0)-x9);                     // atom B216
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y9)");atom.fpos_equation.push_back("z9");atom.fpos_equation.push_back("((1.0/2.0)-x9)");// atom B216 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B216 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B216

    atom.name="A"; atom.type=0;                                       // atom B217
    atom.fpos(1)=x10;atom.fpos(2)=y10;atom.fpos(3)=z10;                     // atom B217
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x10");atom.fpos_equation.push_back("y10");atom.fpos_equation.push_back("z10");// atom B217 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B217 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B217

    atom.name="A"; atom.type=0;                                       // atom B218
    atom.fpos(1)=((1.0/2.0)-x10);atom.fpos(2)=-y10;atom.fpos(3)=((1.0/2.0)+z10);                     // atom B218
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x10)");atom.fpos_equation.push_back("-y10");atom.fpos_equation.push_back("((1.0/2.0)+z10)");// atom B218 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B218 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B218

    atom.name="A"; atom.type=0;                                       // atom B219
    atom.fpos(1)=-x10;atom.fpos(2)=((1.0/2.0)+y10);atom.fpos(3)=((1.0/2.0)-z10);                     // atom B219
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x10");atom.fpos_equation.push_back("((1.0/2.0)+y10)");atom.fpos_equation.push_back("((1.0/2.0)-z10)");// atom B219 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B219 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B219

    atom.name="A"; atom.type=0;                                       // atom B220
    atom.fpos(1)=((1.0/2.0)+x10);atom.fpos(2)=((1.0/2.0)-y10);atom.fpos(3)=-z10;                     // atom B220
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x10)");atom.fpos_equation.push_back("((1.0/2.0)-y10)");atom.fpos_equation.push_back("-z10");// atom B220 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B220 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B220

    atom.name="A"; atom.type=0;                                       // atom B221
    atom.fpos(1)=z10;atom.fpos(2)=x10;atom.fpos(3)=y10;                     // atom B221
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z10");atom.fpos_equation.push_back("x10");atom.fpos_equation.push_back("y10");// atom B221 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B221 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B221

    atom.name="A"; atom.type=0;                                       // atom B222
    atom.fpos(1)=((1.0/2.0)+z10);atom.fpos(2)=((1.0/2.0)-x10);atom.fpos(3)=-y10;                     // atom B222
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z10)");atom.fpos_equation.push_back("((1.0/2.0)-x10)");atom.fpos_equation.push_back("-y10");// atom B222 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B222 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B222

    atom.name="A"; atom.type=0;                                       // atom B223
    atom.fpos(1)=((1.0/2.0)-z10);atom.fpos(2)=-x10;atom.fpos(3)=((1.0/2.0)+y10);                     // atom B223
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z10)");atom.fpos_equation.push_back("-x10");atom.fpos_equation.push_back("((1.0/2.0)+y10)");// atom B223 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B223 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B223

    atom.name="A"; atom.type=0;                                       // atom B224
    atom.fpos(1)=-z10;atom.fpos(2)=((1.0/2.0)+x10);atom.fpos(3)=((1.0/2.0)-y10);                     // atom B224
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z10");atom.fpos_equation.push_back("((1.0/2.0)+x10)");atom.fpos_equation.push_back("((1.0/2.0)-y10)");// atom B224 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B224 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B224

    atom.name="A"; atom.type=0;                                       // atom B225
    atom.fpos(1)=y10;atom.fpos(2)=z10;atom.fpos(3)=x10;                     // atom B225
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y10");atom.fpos_equation.push_back("z10");atom.fpos_equation.push_back("x10");// atom B225 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B225 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B225

    atom.name="A"; atom.type=0;                                       // atom B226
    atom.fpos(1)=-y10;atom.fpos(2)=((1.0/2.0)+z10);atom.fpos(3)=((1.0/2.0)-x10);                     // atom B226
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y10");atom.fpos_equation.push_back("((1.0/2.0)+z10)");atom.fpos_equation.push_back("((1.0/2.0)-x10)");// atom B226 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B226 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B226

    atom.name="A"; atom.type=0;                                       // atom B227
    atom.fpos(1)=((1.0/2.0)+y10);atom.fpos(2)=((1.0/2.0)-z10);atom.fpos(3)=-x10;                     // atom B227
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y10)");atom.fpos_equation.push_back("((1.0/2.0)-z10)");atom.fpos_equation.push_back("-x10");// atom B227 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B227 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B227

    atom.name="A"; atom.type=0;                                       // atom B228
    atom.fpos(1)=((1.0/2.0)-y10);atom.fpos(2)=-z10;atom.fpos(3)=((1.0/2.0)+x10);                     // atom B228
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y10)");atom.fpos_equation.push_back("-z10");atom.fpos_equation.push_back("((1.0/2.0)+x10)");// atom B228 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B228 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B228

    atom.name="A"; atom.type=0;                                       // atom B229
    atom.fpos(1)=-x10;atom.fpos(2)=-y10;atom.fpos(3)=-z10;                     // atom B229
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x10");atom.fpos_equation.push_back("-y10");atom.fpos_equation.push_back("-z10");// atom B229 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B229 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B229

    atom.name="A"; atom.type=0;                                       // atom B230
    atom.fpos(1)=((1.0/2.0)+x10);atom.fpos(2)=y10;atom.fpos(3)=((1.0/2.0)-z10);                     // atom B230
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x10)");atom.fpos_equation.push_back("y10");atom.fpos_equation.push_back("((1.0/2.0)-z10)");// atom B230 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B230 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B230

    atom.name="A"; atom.type=0;                                       // atom B231
    atom.fpos(1)=x10;atom.fpos(2)=((1.0/2.0)-y10);atom.fpos(3)=((1.0/2.0)+z10);                     // atom B231
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x10");atom.fpos_equation.push_back("((1.0/2.0)-y10)");atom.fpos_equation.push_back("((1.0/2.0)+z10)");// atom B231 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B231 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B231

    atom.name="A"; atom.type=0;                                       // atom B232
    atom.fpos(1)=((1.0/2.0)-x10);atom.fpos(2)=((1.0/2.0)+y10);atom.fpos(3)=z10;                     // atom B232
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x10)");atom.fpos_equation.push_back("((1.0/2.0)+y10)");atom.fpos_equation.push_back("z10");// atom B232 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B232 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B232

    atom.name="A"; atom.type=0;                                       // atom B233
    atom.fpos(1)=-z10;atom.fpos(2)=-x10;atom.fpos(3)=-y10;                     // atom B233
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z10");atom.fpos_equation.push_back("-x10");atom.fpos_equation.push_back("-y10");// atom B233 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B233 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B233

    atom.name="A"; atom.type=0;                                       // atom B234
    atom.fpos(1)=((1.0/2.0)-z10);atom.fpos(2)=((1.0/2.0)+x10);atom.fpos(3)=y10;                     // atom B234
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z10)");atom.fpos_equation.push_back("((1.0/2.0)+x10)");atom.fpos_equation.push_back("y10");// atom B234 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B234 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B234

    atom.name="A"; atom.type=0;                                       // atom B235
    atom.fpos(1)=((1.0/2.0)+z10);atom.fpos(2)=x10;atom.fpos(3)=((1.0/2.0)-y10);                     // atom B235
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z10)");atom.fpos_equation.push_back("x10");atom.fpos_equation.push_back("((1.0/2.0)-y10)");// atom B235 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B235 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B235

    atom.name="A"; atom.type=0;                                       // atom B236
    atom.fpos(1)=z10;atom.fpos(2)=((1.0/2.0)-x10);atom.fpos(3)=((1.0/2.0)+y10);                     // atom B236
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z10");atom.fpos_equation.push_back("((1.0/2.0)-x10)");atom.fpos_equation.push_back("((1.0/2.0)+y10)");// atom B236 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B236 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B236

    atom.name="A"; atom.type=0;                                       // atom B237
    atom.fpos(1)=-y10;atom.fpos(2)=-z10;atom.fpos(3)=-x10;                     // atom B237
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y10");atom.fpos_equation.push_back("-z10");atom.fpos_equation.push_back("-x10");// atom B237 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B237 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B237

    atom.name="A"; atom.type=0;                                       // atom B238
    atom.fpos(1)=y10;atom.fpos(2)=((1.0/2.0)-z10);atom.fpos(3)=((1.0/2.0)+x10);                     // atom B238
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y10");atom.fpos_equation.push_back("((1.0/2.0)-z10)");atom.fpos_equation.push_back("((1.0/2.0)+x10)");// atom B238 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B238 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B238

    atom.name="A"; atom.type=0;                                       // atom B239
    atom.fpos(1)=((1.0/2.0)-y10);atom.fpos(2)=((1.0/2.0)+z10);atom.fpos(3)=x10;                     // atom B239
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y10)");atom.fpos_equation.push_back("((1.0/2.0)+z10)");atom.fpos_equation.push_back("x10");// atom B239 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B239 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B239

    atom.name="A"; atom.type=0;                                       // atom B240
    atom.fpos(1)=((1.0/2.0)+y10);atom.fpos(2)=z10;atom.fpos(3)=((1.0/2.0)-x10);                     // atom B240
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y10)");atom.fpos_equation.push_back("z10");atom.fpos_equation.push_back("((1.0/2.0)-x10)");// atom B240 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B240 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B240


    return str.atoms.size();  
  }
} // namespace anrl

namespace anrl {
  uint WebANRL_A_cP240_205_10d(stringstream& web,bool LDEBUG) {
#ifndef _ANRL_NOWEB_
#endif

    if(LDEBUG) {cerr << "anrl:: WebANRL_A_cP240_205_10d: web.str().size()=" << web.str().size() << endl;}

    return web.str().size();
  }
} // namespace anrl

#endif // AFLOW_REMOVE_GREP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

