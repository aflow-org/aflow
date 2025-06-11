// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - David Hicks - 2018
// FILE "ANRL/aflow_anrl_ABC3_oP20_30_2a_c_3c.cpp"

#ifndef _AFLOW_ANRL_ABC3_oP20_30_2a_c_3c_CPP // AFLOW_REMOVE_GREP
#define _AFLOW_ANRL_ABC3_oP20_30_2a_c_3c_CPP // AFLOW_REMOVE_GREP
#include "../aflow.h" // AFLOW_REMOVE_GREP

namespace anrl {
  uint WebANRL_ABC3_oP20_30_2a_c_3c(stringstream &web,bool LDEBUG);
}

namespace anrl {
  uint PrototypeANRL_ABC3_oP20_30_2a_c_3c(stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG) {
    // system ABC3_oP20_30_2a_c_3c

    if(XHOST.vflag_control.flag("WWW")) {
      WebANRL_ABC3_oP20_30_2a_c_3c(web,LDEBUG); // PLUG WEB STUFF
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

    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: FOUND" << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: label=" << label << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: nspecies=" << nspecies << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: natoms=" << natoms << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: spacegroup=" << spacegroup << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: nunderscores=" << nunderscores << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: nparameters=" <<  nparameters << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: Pearson_symbol=" << Pearson_symbol << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: params=" << params << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: Strukturbericht=" << Strukturbericht << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: prototype=" << prototype << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: dialect=" << dialect << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: vparameters.size()=" << vparameters.size() << endl;}

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
    double a=vparameters.at(i++);                  if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: a=" << a << endl;}
    double bovera=vparameters.at(i++),b=bovera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: b=" << b << " (b/a=" << bovera << ")" << endl;}
    double covera=vparameters.at(i++),c=covera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: c=" << c << " (c/a=" << covera << ")" << endl;}

    double z1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: z1=" << z1 << endl;}
    double z2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: z2=" << z2 << endl;}
    double x3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: x3=" << x3 << endl;}
    double y3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: y3=" << y3 << endl;}
    double z3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: z3=" << z3 << endl;}
    double x4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: x4=" << x4 << endl;}
    double y4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: y4=" << y4 << endl;}
    double z4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: z4=" << z4 << endl;}
    double x5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: x5=" << x5 << endl;}
    double y5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: y5=" << y5 << endl;}
    double z5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: z5=" << z5 << endl;}
    double x6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: x6=" << x6 << endl;}
    double y6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: y6=" << y6 << endl;}
    double z6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC3_oP20_30_2a_c_3c: z6=" << z6 << endl;}

    str.iomode=IOVASP_AUTO;
    str.title=label+" params="+parameters+" SG="+aurostd::utype2string(spacegroup)+DOI_ANRL; //CO20190520
    str.scale=1.0;

    a1=a*xn;
    a2=b*yn;
    a3=c*zn;  

    str.lattice(1,1)=a1(1);str.lattice(1,2)=a1(2);str.lattice(1,3)=a1(3);
    str.lattice(2,1)=a2(1);str.lattice(2,2)=a2(2);str.lattice(2,3)=a2(3);
    str.lattice(3,1)=a3(1);str.lattice(3,2)=a3(2);str.lattice(3,3)=a3(3);

    // symbolic representation of lattice vectors
    vector<string> a1_equation, a2_equation, a3_equation;
    a1_equation.push_back("a");a1_equation.push_back("0");a1_equation.push_back("0");
    a2_equation.push_back("0");a2_equation.push_back("b");a2_equation.push_back("0");
    a3_equation.push_back("0");a3_equation.push_back("0");a3_equation.push_back("c");
    str.symbolic_math_lattice.push_back(a1_equation);
    str.symbolic_math_lattice.push_back(a2_equation);
    str.symbolic_math_lattice.push_back(a3_equation);

    str.num_lattice_parameters = 3;

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
    atom.fpos(1)=0.0;atom.fpos(2)=(1.0/2.0);atom.fpos(3)=((1.0/2.0)+z1);                     // atom B2
    atom.fpos_equation.clear();atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("((1.0/2.0)+z1)");// atom B2 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B2 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B2

    atom.name="A"; atom.type=0;                                       // atom B3
    atom.fpos(1)=0.0;atom.fpos(2)=0.0;atom.fpos(3)=z2;                     // atom B3
    atom.fpos_equation.clear();atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("z2");// atom B3 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B3 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B3

    atom.name="A"; atom.type=0;                                       // atom B4
    atom.fpos(1)=0.0;atom.fpos(2)=(1.0/2.0);atom.fpos(3)=((1.0/2.0)+z2);                     // atom B4
    atom.fpos_equation.clear();atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("((1.0/2.0)+z2)");// atom B4 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B4 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B4

    atom.name="B"; atom.type=1;                                       // atom B5
    atom.fpos(1)=x3;atom.fpos(2)=y3;atom.fpos(3)=z3;                     // atom B5
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("y3");atom.fpos_equation.push_back("z3");// atom B5 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B5 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B5

    atom.name="B"; atom.type=1;                                       // atom B6
    atom.fpos(1)=-x3;atom.fpos(2)=-y3;atom.fpos(3)=z3;                     // atom B6
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("-y3");atom.fpos_equation.push_back("z3");// atom B6 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B6 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B6

    atom.name="B"; atom.type=1;                                       // atom B7
    atom.fpos(1)=x3;atom.fpos(2)=((1.0/2.0)-y3);atom.fpos(3)=((1.0/2.0)+z3);                     // atom B7
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("((1.0/2.0)-y3)");atom.fpos_equation.push_back("((1.0/2.0)+z3)");// atom B7 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B7 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B7

    atom.name="B"; atom.type=1;                                       // atom B8
    atom.fpos(1)=-x3;atom.fpos(2)=((1.0/2.0)+y3);atom.fpos(3)=((1.0/2.0)+z3);                     // atom B8
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("((1.0/2.0)+y3)");atom.fpos_equation.push_back("((1.0/2.0)+z3)");// atom B8 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B8 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B8

    atom.name="C"; atom.type=2;                                       // atom B9
    atom.fpos(1)=x4;atom.fpos(2)=y4;atom.fpos(3)=z4;                     // atom B9
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x4");atom.fpos_equation.push_back("y4");atom.fpos_equation.push_back("z4");// atom B9 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B9 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B9

    atom.name="C"; atom.type=2;                                       // atom B10
    atom.fpos(1)=-x4;atom.fpos(2)=-y4;atom.fpos(3)=z4;                     // atom B10
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x4");atom.fpos_equation.push_back("-y4");atom.fpos_equation.push_back("z4");// atom B10 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B10 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B10

    atom.name="C"; atom.type=2;                                       // atom B11
    atom.fpos(1)=x4;atom.fpos(2)=((1.0/2.0)-y4);atom.fpos(3)=((1.0/2.0)+z4);                     // atom B11
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x4");atom.fpos_equation.push_back("((1.0/2.0)-y4)");atom.fpos_equation.push_back("((1.0/2.0)+z4)");// atom B11 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B11 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B11

    atom.name="C"; atom.type=2;                                       // atom B12
    atom.fpos(1)=-x4;atom.fpos(2)=((1.0/2.0)+y4);atom.fpos(3)=((1.0/2.0)+z4);                     // atom B12
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x4");atom.fpos_equation.push_back("((1.0/2.0)+y4)");atom.fpos_equation.push_back("((1.0/2.0)+z4)");// atom B12 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B12 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B12

    atom.name="C"; atom.type=2;                                       // atom B13
    atom.fpos(1)=x5;atom.fpos(2)=y5;atom.fpos(3)=z5;                     // atom B13
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("y5");atom.fpos_equation.push_back("z5");// atom B13 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B13 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B13

    atom.name="C"; atom.type=2;                                       // atom B14
    atom.fpos(1)=-x5;atom.fpos(2)=-y5;atom.fpos(3)=z5;                     // atom B14
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("-y5");atom.fpos_equation.push_back("z5");// atom B14 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B14 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B14

    atom.name="C"; atom.type=2;                                       // atom B15
    atom.fpos(1)=x5;atom.fpos(2)=((1.0/2.0)-y5);atom.fpos(3)=((1.0/2.0)+z5);                     // atom B15
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("((1.0/2.0)-y5)");atom.fpos_equation.push_back("((1.0/2.0)+z5)");// atom B15 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B15 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B15

    atom.name="C"; atom.type=2;                                       // atom B16
    atom.fpos(1)=-x5;atom.fpos(2)=((1.0/2.0)+y5);atom.fpos(3)=((1.0/2.0)+z5);                     // atom B16
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("((1.0/2.0)+y5)");atom.fpos_equation.push_back("((1.0/2.0)+z5)");// atom B16 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B16 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B16

    atom.name="C"; atom.type=2;                                       // atom B17
    atom.fpos(1)=x6;atom.fpos(2)=y6;atom.fpos(3)=z6;                     // atom B17
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("y6");atom.fpos_equation.push_back("z6");// atom B17 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B17 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B17

    atom.name="C"; atom.type=2;                                       // atom B18
    atom.fpos(1)=-x6;atom.fpos(2)=-y6;atom.fpos(3)=z6;                     // atom B18
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("-y6");atom.fpos_equation.push_back("z6");// atom B18 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B18 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B18

    atom.name="C"; atom.type=2;                                       // atom B19
    atom.fpos(1)=x6;atom.fpos(2)=((1.0/2.0)-y6);atom.fpos(3)=((1.0/2.0)+z6);                     // atom B19
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("((1.0/2.0)-y6)");atom.fpos_equation.push_back("((1.0/2.0)+z6)");// atom B19 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B19 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B19

    atom.name="C"; atom.type=2;                                       // atom B20
    atom.fpos(1)=-x6;atom.fpos(2)=((1.0/2.0)+y6);atom.fpos(3)=((1.0/2.0)+z6);                     // atom B20
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("((1.0/2.0)+y6)");atom.fpos_equation.push_back("((1.0/2.0)+z6)");// atom B20 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B20 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B20


    return str.atoms.size();  
  }
} // namespace anrl

namespace anrl {
  uint WebANRL_ABC3_oP20_30_2a_c_3c(stringstream& web,bool LDEBUG) {
#ifndef _ANRL_NOWEB_
#endif

    if(LDEBUG) {cerr << "anrl:: WebANRL_ABC3_oP20_30_2a_c_3c: web.str().size()=" << web.str().size() << endl;}

    return web.str().size();
  }
} // namespace anrl

#endif // AFLOW_REMOVE_GREP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

