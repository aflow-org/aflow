// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - David Hicks - 2018
// FILE "ANRL/aflow_anrl_A5BCD6_cF416_228_eg_c_b_h.cpp"

#ifndef _AFLOW_ANRL_A5BCD6_cF416_228_eg_c_b_h_CPP // AFLOW_REMOVE_GREP
#define _AFLOW_ANRL_A5BCD6_cF416_228_eg_c_b_h_CPP // AFLOW_REMOVE_GREP
#include "../aflow.h" // AFLOW_REMOVE_GREP

namespace anrl {
  uint WebANRL_A5BCD6_cF416_228_eg_c_b_h(stringstream &web,bool LDEBUG);
}

namespace anrl {
  uint PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h(stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG) {
    // system A5BCD6_cF416_228_eg_c_b_h

    if(XHOST.vflag_control.flag("WWW")) {
      WebANRL_A5BCD6_cF416_228_eg_c_b_h(web,LDEBUG); // PLUG WEB STUFF
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

    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h: FOUND" << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h: label=" << label << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h: nspecies=" << nspecies << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h: natoms=" << natoms << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h: spacegroup=" << spacegroup << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h: nunderscores=" << nunderscores << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h: nparameters=" <<  nparameters << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h: Pearson_symbol=" << Pearson_symbol << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h: params=" << params << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h: Strukturbericht=" << Strukturbericht << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h: prototype=" << prototype << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h: dialect=" << dialect << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h: vparameters.size()=" << vparameters.size() << endl;}

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
    double a=vparameters.at(i++);                  if(LDEBUG) { cerr << "anrl::PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h: a=" << a << endl;}

    double x3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h: x3=" << x3 << endl;}
    double y4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h: y4=" << y4 << endl;}
    double x5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h: x5=" << x5 << endl;}
    double y5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h: y5=" << y5 << endl;}
    double z5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h: z5=" << z5 << endl;}

    str.iomode=IOVASP_AUTO;
    str.title=label+" params="+parameters+" SG="+aurostd::utype2string(spacegroup)+DOI_ANRL; //CO20190520
    str.scale=1.0;

    a1=(1.0/2.0)*a*yn+(1.0/2.0)*a*zn;
    a2=(1.0/2.0)*a*xn+(1.0/2.0)*a*zn;
    a3=(1.0/2.0)*a*xn+(1.0/2.0)*a*yn;

    str.lattice(1,1)=a1(1);str.lattice(1,2)=a1(2);str.lattice(1,3)=a1(3);
    str.lattice(2,1)=a2(1);str.lattice(2,2)=a2(2);str.lattice(2,3)=a2(3);
    str.lattice(3,1)=a3(1);str.lattice(3,2)=a3(2);str.lattice(3,3)=a3(3);

    // symbolic representation of lattice vectors
    vector<string> a1_equation, a2_equation, a3_equation;
    a1_equation.push_back("0");a1_equation.push_back("(1.0/2.0)*a");a1_equation.push_back("(1.0/2.0)*a");
    a2_equation.push_back("(1.0/2.0)*a");a2_equation.push_back("0");a2_equation.push_back("(1.0/2.0)*a");
    a3_equation.push_back("(1.0/2.0)*a");a3_equation.push_back("(1.0/2.0)*a");a3_equation.push_back("0");
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

    atom.name="A"; atom.type=0;                                       // atom B17
    atom.fpos(1)=x3;atom.fpos(2)=x3;atom.fpos(3)=x3;                     // atom B17
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("x3");// atom B17 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B17 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B17

    atom.name="A"; atom.type=0;                                       // atom B18
    atom.fpos(1)=x3;atom.fpos(2)=x3;atom.fpos(3)=((1.0/2.0)-3.0*x3);                     // atom B18
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("((1.0/2.0)-3.0*x3)");// atom B18 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B18 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B18

    atom.name="A"; atom.type=0;                                       // atom B19
    atom.fpos(1)=x3;atom.fpos(2)=((1.0/2.0)-3.0*x3);atom.fpos(3)=x3;                     // atom B19
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("((1.0/2.0)-3.0*x3)");atom.fpos_equation.push_back("x3");// atom B19 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B19 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B19

    atom.name="A"; atom.type=0;                                       // atom B20
    atom.fpos(1)=((1.0/2.0)-3.0*x3);atom.fpos(2)=x3;atom.fpos(3)=x3;                     // atom B20
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-3.0*x3)");atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("x3");// atom B20 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B20 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B20

    atom.name="A"; atom.type=0;                                       // atom B21
    atom.fpos(1)=((1.0/2.0)-x3);atom.fpos(2)=((1.0/2.0)-x3);atom.fpos(3)=3.0*x3;                     // atom B21
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x3)");atom.fpos_equation.push_back("((1.0/2.0)-x3)");atom.fpos_equation.push_back("3.0*x3");// atom B21 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B21 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B21

    atom.name="A"; atom.type=0;                                       // atom B22
    atom.fpos(1)=((1.0/2.0)-x3);atom.fpos(2)=((1.0/2.0)-x3);atom.fpos(3)=((1.0/2.0)-x3);                     // atom B22
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x3)");atom.fpos_equation.push_back("((1.0/2.0)-x3)");atom.fpos_equation.push_back("((1.0/2.0)-x3)");// atom B22 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B22 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B22

    atom.name="A"; atom.type=0;                                       // atom B23
    atom.fpos(1)=((1.0/2.0)-x3);atom.fpos(2)=3.0*x3;atom.fpos(3)=((1.0/2.0)-x3);                     // atom B23
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x3)");atom.fpos_equation.push_back("3.0*x3");atom.fpos_equation.push_back("((1.0/2.0)-x3)");// atom B23 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B23 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B23

    atom.name="A"; atom.type=0;                                       // atom B24
    atom.fpos(1)=3.0*x3;atom.fpos(2)=((1.0/2.0)-x3);atom.fpos(3)=((1.0/2.0)-x3);                     // atom B24
    atom.fpos_equation.clear();atom.fpos_equation.push_back("3.0*x3");atom.fpos_equation.push_back("((1.0/2.0)-x3)");atom.fpos_equation.push_back("((1.0/2.0)-x3)");// atom B24 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B24 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B24

    atom.name="A"; atom.type=0;                                       // atom B25
    atom.fpos(1)=-x3;atom.fpos(2)=-x3;atom.fpos(3)=-x3;                     // atom B25
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("-x3");// atom B25 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B25 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B25

    atom.name="A"; atom.type=0;                                       // atom B26
    atom.fpos(1)=-x3;atom.fpos(2)=-x3;atom.fpos(3)=((1.0/2.0)+3.0*x3);                     // atom B26
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("((1.0/2.0)+3.0*x3)");// atom B26 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B26 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B26

    atom.name="A"; atom.type=0;                                       // atom B27
    atom.fpos(1)=-x3;atom.fpos(2)=((1.0/2.0)+3.0*x3);atom.fpos(3)=-x3;                     // atom B27
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("((1.0/2.0)+3.0*x3)");atom.fpos_equation.push_back("-x3");// atom B27 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B27 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B27

    atom.name="A"; atom.type=0;                                       // atom B28
    atom.fpos(1)=((1.0/2.0)+3.0*x3);atom.fpos(2)=-x3;atom.fpos(3)=-x3;                     // atom B28
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+3.0*x3)");atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("-x3");// atom B28 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B28 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B28

    atom.name="A"; atom.type=0;                                       // atom B29
    atom.fpos(1)=((1.0/2.0)+x3);atom.fpos(2)=((1.0/2.0)+x3);atom.fpos(3)=-3.0*x3;                     // atom B29
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x3)");atom.fpos_equation.push_back("((1.0/2.0)+x3)");atom.fpos_equation.push_back("-3.0*x3");// atom B29 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B29 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B29

    atom.name="A"; atom.type=0;                                       // atom B30
    atom.fpos(1)=((1.0/2.0)+x3);atom.fpos(2)=((1.0/2.0)+x3);atom.fpos(3)=((1.0/2.0)+x3);                     // atom B30
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x3)");atom.fpos_equation.push_back("((1.0/2.0)+x3)");atom.fpos_equation.push_back("((1.0/2.0)+x3)");// atom B30 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B30 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B30

    atom.name="A"; atom.type=0;                                       // atom B31
    atom.fpos(1)=((1.0/2.0)+x3);atom.fpos(2)=-3.0*x3;atom.fpos(3)=((1.0/2.0)+x3);                     // atom B31
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x3)");atom.fpos_equation.push_back("-3.0*x3");atom.fpos_equation.push_back("((1.0/2.0)+x3)");// atom B31 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B31 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B31

    atom.name="A"; atom.type=0;                                       // atom B32
    atom.fpos(1)=-3.0*x3;atom.fpos(2)=((1.0/2.0)+x3);atom.fpos(3)=((1.0/2.0)+x3);                     // atom B32
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-3.0*x3");atom.fpos_equation.push_back("((1.0/2.0)+x3)");atom.fpos_equation.push_back("((1.0/2.0)+x3)");// atom B32 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B32 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B32

    atom.name="A"; atom.type=0;                                       // atom B33
    atom.fpos(1)=(3.0/4.0);atom.fpos(2)=((1.0/4.0)-2.0*y4);atom.fpos(3)=((1.0/4.0)+2.0*y4);                     // atom B33
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("((1.0/4.0)-2.0*y4)");atom.fpos_equation.push_back("((1.0/4.0)+2.0*y4)");// atom B33 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B33 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B33

    atom.name="A"; atom.type=0;                                       // atom B34
    atom.fpos(1)=((1.0/4.0)-2.0*y4);atom.fpos(2)=(3.0/4.0);atom.fpos(3)=(1.0/4.0);                     // atom B34
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/4.0)-2.0*y4)");atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("(1.0/4.0)");// atom B34 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B34 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B34

    atom.name="A"; atom.type=0;                                       // atom B35
    atom.fpos(1)=((1.0/4.0)+2.0*y4);atom.fpos(2)=(1.0/4.0);atom.fpos(3)=(3.0/4.0);                     // atom B35
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/4.0)+2.0*y4)");atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("(3.0/4.0)");// atom B35 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B35 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B35

    atom.name="A"; atom.type=0;                                       // atom B36
    atom.fpos(1)=(1.0/4.0);atom.fpos(2)=((1.0/4.0)+2.0*y4);atom.fpos(3)=((1.0/4.0)-2.0*y4);                     // atom B36
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("((1.0/4.0)+2.0*y4)");atom.fpos_equation.push_back("((1.0/4.0)-2.0*y4)");// atom B36 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B36 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B36

    atom.name="A"; atom.type=0;                                       // atom B37
    atom.fpos(1)=((1.0/4.0)+2.0*y4);atom.fpos(2)=(3.0/4.0);atom.fpos(3)=((1.0/4.0)-2.0*y4);                     // atom B37
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/4.0)+2.0*y4)");atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("((1.0/4.0)-2.0*y4)");// atom B37 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B37 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B37

    atom.name="A"; atom.type=0;                                       // atom B38
    atom.fpos(1)=(1.0/4.0);atom.fpos(2)=((1.0/4.0)-2.0*y4);atom.fpos(3)=(3.0/4.0);                     // atom B38
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("((1.0/4.0)-2.0*y4)");atom.fpos_equation.push_back("(3.0/4.0)");// atom B38 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B38 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B38

    atom.name="A"; atom.type=0;                                       // atom B39
    atom.fpos(1)=(3.0/4.0);atom.fpos(2)=((1.0/4.0)+2.0*y4);atom.fpos(3)=(1.0/4.0);                     // atom B39
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("((1.0/4.0)+2.0*y4)");atom.fpos_equation.push_back("(1.0/4.0)");// atom B39 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B39 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B39

    atom.name="A"; atom.type=0;                                       // atom B40
    atom.fpos(1)=((1.0/4.0)-2.0*y4);atom.fpos(2)=(1.0/4.0);atom.fpos(3)=((1.0/4.0)+2.0*y4);                     // atom B40
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/4.0)-2.0*y4)");atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("((1.0/4.0)+2.0*y4)");// atom B40 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B40 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B40

    atom.name="A"; atom.type=0;                                       // atom B41
    atom.fpos(1)=((1.0/4.0)-2.0*y4);atom.fpos(2)=((1.0/4.0)+2.0*y4);atom.fpos(3)=(3.0/4.0);                     // atom B41
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/4.0)-2.0*y4)");atom.fpos_equation.push_back("((1.0/4.0)+2.0*y4)");atom.fpos_equation.push_back("(3.0/4.0)");// atom B41 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B41 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B41

    atom.name="A"; atom.type=0;                                       // atom B42
    atom.fpos(1)=(3.0/4.0);atom.fpos(2)=(1.0/4.0);atom.fpos(3)=((1.0/4.0)-2.0*y4);                     // atom B42
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("((1.0/4.0)-2.0*y4)");// atom B42 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B42 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B42

    atom.name="A"; atom.type=0;                                       // atom B43
    atom.fpos(1)=(1.0/4.0);atom.fpos(2)=(3.0/4.0);atom.fpos(3)=((1.0/4.0)+2.0*y4);                     // atom B43
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("((1.0/4.0)+2.0*y4)");// atom B43 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B43 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B43

    atom.name="A"; atom.type=0;                                       // atom B44
    atom.fpos(1)=((1.0/4.0)+2.0*y4);atom.fpos(2)=((1.0/4.0)-2.0*y4);atom.fpos(3)=(1.0/4.0);                     // atom B44
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/4.0)+2.0*y4)");atom.fpos_equation.push_back("((1.0/4.0)-2.0*y4)");atom.fpos_equation.push_back("(1.0/4.0)");// atom B44 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B44 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B44

    atom.name="A"; atom.type=0;                                       // atom B45
    atom.fpos(1)=(1.0/4.0);atom.fpos(2)=((3.0/4.0)+2.0*y4);atom.fpos(3)=((3.0/4.0)-2.0*y4);                     // atom B45
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("((3.0/4.0)+2.0*y4)");atom.fpos_equation.push_back("((3.0/4.0)-2.0*y4)");// atom B45 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B45 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B45

    atom.name="A"; atom.type=0;                                       // atom B46
    atom.fpos(1)=((3.0/4.0)+2.0*y4);atom.fpos(2)=(1.0/4.0);atom.fpos(3)=(3.0/4.0);                     // atom B46
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((3.0/4.0)+2.0*y4)");atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("(3.0/4.0)");// atom B46 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B46 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B46

    atom.name="A"; atom.type=0;                                       // atom B47
    atom.fpos(1)=((3.0/4.0)-2.0*y4);atom.fpos(2)=(3.0/4.0);atom.fpos(3)=(1.0/4.0);                     // atom B47
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((3.0/4.0)-2.0*y4)");atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("(1.0/4.0)");// atom B47 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B47 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B47

    atom.name="A"; atom.type=0;                                       // atom B48
    atom.fpos(1)=(3.0/4.0);atom.fpos(2)=((3.0/4.0)-2.0*y4);atom.fpos(3)=((3.0/4.0)+2.0*y4);                     // atom B48
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("((3.0/4.0)-2.0*y4)");atom.fpos_equation.push_back("((3.0/4.0)+2.0*y4)");// atom B48 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B48 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B48

    atom.name="A"; atom.type=0;                                       // atom B49
    atom.fpos(1)=((3.0/4.0)-2.0*y4);atom.fpos(2)=(1.0/4.0);atom.fpos(3)=((3.0/4.0)+2.0*y4);                     // atom B49
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((3.0/4.0)-2.0*y4)");atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("((3.0/4.0)+2.0*y4)");// atom B49 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B49 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B49

    atom.name="A"; atom.type=0;                                       // atom B50
    atom.fpos(1)=(3.0/4.0);atom.fpos(2)=((3.0/4.0)+2.0*y4);atom.fpos(3)=(1.0/4.0);                     // atom B50
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("((3.0/4.0)+2.0*y4)");atom.fpos_equation.push_back("(1.0/4.0)");// atom B50 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B50 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B50

    atom.name="A"; atom.type=0;                                       // atom B51
    atom.fpos(1)=(1.0/4.0);atom.fpos(2)=((3.0/4.0)-2.0*y4);atom.fpos(3)=(3.0/4.0);                     // atom B51
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("((3.0/4.0)-2.0*y4)");atom.fpos_equation.push_back("(3.0/4.0)");// atom B51 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B51 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B51

    atom.name="A"; atom.type=0;                                       // atom B52
    atom.fpos(1)=((3.0/4.0)+2.0*y4);atom.fpos(2)=(3.0/4.0);atom.fpos(3)=((3.0/4.0)-2.0*y4);                     // atom B52
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((3.0/4.0)+2.0*y4)");atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("((3.0/4.0)-2.0*y4)");// atom B52 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B52 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B52

    atom.name="A"; atom.type=0;                                       // atom B53
    atom.fpos(1)=((3.0/4.0)+2.0*y4);atom.fpos(2)=((3.0/4.0)-2.0*y4);atom.fpos(3)=(1.0/4.0);                     // atom B53
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((3.0/4.0)+2.0*y4)");atom.fpos_equation.push_back("((3.0/4.0)-2.0*y4)");atom.fpos_equation.push_back("(1.0/4.0)");// atom B53 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B53 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B53

    atom.name="A"; atom.type=0;                                       // atom B54
    atom.fpos(1)=(1.0/4.0);atom.fpos(2)=(3.0/4.0);atom.fpos(3)=((3.0/4.0)+2.0*y4);                     // atom B54
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("((3.0/4.0)+2.0*y4)");// atom B54 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B54 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B54

    atom.name="A"; atom.type=0;                                       // atom B55
    atom.fpos(1)=(3.0/4.0);atom.fpos(2)=(1.0/4.0);atom.fpos(3)=((3.0/4.0)-2.0*y4);                     // atom B55
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("((3.0/4.0)-2.0*y4)");// atom B55 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B55 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B55

    atom.name="A"; atom.type=0;                                       // atom B56
    atom.fpos(1)=((3.0/4.0)-2.0*y4);atom.fpos(2)=((3.0/4.0)+2.0*y4);atom.fpos(3)=(3.0/4.0);                     // atom B56
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((3.0/4.0)-2.0*y4)");atom.fpos_equation.push_back("((3.0/4.0)+2.0*y4)");atom.fpos_equation.push_back("(3.0/4.0)");// atom B56 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B56 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B56

    atom.name="B"; atom.type=1;                                       // atom B9
    atom.fpos(1)=0.0;atom.fpos(2)=0.0;atom.fpos(3)=0.0;                     // atom B9
    atom.fpos_equation.clear();atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("0.0");// atom B9 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B9 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B9

    atom.name="B"; atom.type=1;                                       // atom B10
    atom.fpos(1)=0.0;atom.fpos(2)=0.0;atom.fpos(3)=(1.0/2.0);                     // atom B10
    atom.fpos_equation.clear();atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("(1.0/2.0)");// atom B10 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B10 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B10

    atom.name="B"; atom.type=1;                                       // atom B11
    atom.fpos(1)=0.0;atom.fpos(2)=(1.0/2.0);atom.fpos(3)=0.0;                     // atom B11
    atom.fpos_equation.clear();atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("0.0");// atom B11 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B11 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B11

    atom.name="B"; atom.type=1;                                       // atom B12
    atom.fpos(1)=(1.0/2.0);atom.fpos(2)=0.0;atom.fpos(3)=0.0;                     // atom B12
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("0.0");// atom B12 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B12 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B12

    atom.name="B"; atom.type=1;                                       // atom B13
    atom.fpos(1)=(1.0/2.0);atom.fpos(2)=(1.0/2.0);atom.fpos(3)=0.0;                     // atom B13
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("0.0");// atom B13 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B13 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B13

    atom.name="B"; atom.type=1;                                       // atom B14
    atom.fpos(1)=(1.0/2.0);atom.fpos(2)=(1.0/2.0);atom.fpos(3)=(1.0/2.0);                     // atom B14
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("(1.0/2.0)");// atom B14 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B14 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B14

    atom.name="B"; atom.type=1;                                       // atom B15
    atom.fpos(1)=(1.0/2.0);atom.fpos(2)=0.0;atom.fpos(3)=(1.0/2.0);                     // atom B15
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("(1.0/2.0)");// atom B15 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B15 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B15

    atom.name="B"; atom.type=1;                                       // atom B16
    atom.fpos(1)=0.0;atom.fpos(2)=(1.0/2.0);atom.fpos(3)=(1.0/2.0);                     // atom B16
    atom.fpos_equation.clear();atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("(1.0/2.0)");// atom B16 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B16 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B16

    atom.name="C"; atom.type=2;                                       // atom B1
    atom.fpos(1)=(1.0/4.0);atom.fpos(2)=(1.0/4.0);atom.fpos(3)=(1.0/4.0);                     // atom B1
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("(1.0/4.0)");// atom B1 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B1 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B1

    atom.name="C"; atom.type=2;                                       // atom B2
    atom.fpos(1)=(1.0/4.0);atom.fpos(2)=(1.0/4.0);atom.fpos(3)=(3.0/4.0);                     // atom B2
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("(3.0/4.0)");// atom B2 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B2 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B2

    atom.name="C"; atom.type=2;                                       // atom B3
    atom.fpos(1)=(1.0/4.0);atom.fpos(2)=(3.0/4.0);atom.fpos(3)=(1.0/4.0);                     // atom B3
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("(1.0/4.0)");// atom B3 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B3 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B3

    atom.name="C"; atom.type=2;                                       // atom B4
    atom.fpos(1)=(3.0/4.0);atom.fpos(2)=(1.0/4.0);atom.fpos(3)=(1.0/4.0);                     // atom B4
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("(1.0/4.0)");// atom B4 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B4 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B4

    atom.name="C"; atom.type=2;                                       // atom B5
    atom.fpos(1)=(3.0/4.0);atom.fpos(2)=(3.0/4.0);atom.fpos(3)=(3.0/4.0);                     // atom B5
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("(3.0/4.0)");// atom B5 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B5 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B5

    atom.name="C"; atom.type=2;                                       // atom B6
    atom.fpos(1)=(3.0/4.0);atom.fpos(2)=(3.0/4.0);atom.fpos(3)=(1.0/4.0);                     // atom B6
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("(1.0/4.0)");// atom B6 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B6 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B6

    atom.name="C"; atom.type=2;                                       // atom B7
    atom.fpos(1)=(3.0/4.0);atom.fpos(2)=(1.0/4.0);atom.fpos(3)=(3.0/4.0);                     // atom B7
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("(3.0/4.0)");// atom B7 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B7 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B7

    atom.name="C"; atom.type=2;                                       // atom B8
    atom.fpos(1)=(1.0/4.0);atom.fpos(2)=(3.0/4.0);atom.fpos(3)=(3.0/4.0);                     // atom B8
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("(3.0/4.0)");// atom B8 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B8 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B8

    atom.name="D"; atom.type=3;                                       // atom B57
    atom.fpos(1)=(-x5+y5+z5);atom.fpos(2)=(x5-y5+z5);atom.fpos(3)=(x5+y5-z5);                     // atom B57
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x5+y5+z5)");atom.fpos_equation.push_back("(x5-y5+z5)");atom.fpos_equation.push_back("(x5+y5-z5)");// atom B57 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B57 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B57

    atom.name="D"; atom.type=3;                                       // atom B58
    atom.fpos(1)=(x5-y5+z5);atom.fpos(2)=(-x5+y5+z5);atom.fpos(3)=((1.0/2.0)-x5-y5-z5);                     // atom B58
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x5-y5+z5)");atom.fpos_equation.push_back("(-x5+y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5-y5-z5)");// atom B58 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B58 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B58

    atom.name="D"; atom.type=3;                                       // atom B59
    atom.fpos(1)=(x5+y5-z5);atom.fpos(2)=((1.0/2.0)-x5-y5-z5);atom.fpos(3)=(-x5+y5+z5);                     // atom B59
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x5+y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5-y5-z5)");atom.fpos_equation.push_back("(-x5+y5+z5)");// atom B59 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B59 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B59

    atom.name="D"; atom.type=3;                                       // atom B60
    atom.fpos(1)=((1.0/2.0)-x5-y5-z5);atom.fpos(2)=(x5+y5-z5);atom.fpos(3)=(x5-y5+z5);                     // atom B60
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x5-y5-z5)");atom.fpos_equation.push_back("(x5+y5-z5)");atom.fpos_equation.push_back("(x5-y5+z5)");// atom B60 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B60 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B60

    atom.name="D"; atom.type=3;                                       // atom B61
    atom.fpos(1)=(x5+y5-z5);atom.fpos(2)=(-x5+y5+z5);atom.fpos(3)=(x5-y5+z5);                     // atom B61
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x5+y5-z5)");atom.fpos_equation.push_back("(-x5+y5+z5)");atom.fpos_equation.push_back("(x5-y5+z5)");// atom B61 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B61 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B61

    atom.name="D"; atom.type=3;                                       // atom B62
    atom.fpos(1)=((1.0/2.0)-x5-y5-z5);atom.fpos(2)=(x5-y5+z5);atom.fpos(3)=(-x5+y5+z5);                     // atom B62
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x5-y5-z5)");atom.fpos_equation.push_back("(x5-y5+z5)");atom.fpos_equation.push_back("(-x5+y5+z5)");// atom B62 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B62 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B62

    atom.name="D"; atom.type=3;                                       // atom B63
    atom.fpos(1)=(-x5+y5+z5);atom.fpos(2)=(x5+y5-z5);atom.fpos(3)=((1.0/2.0)-x5-y5-z5);                     // atom B63
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x5+y5+z5)");atom.fpos_equation.push_back("(x5+y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5-y5-z5)");// atom B63 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B63 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B63

    atom.name="D"; atom.type=3;                                       // atom B64
    atom.fpos(1)=(x5-y5+z5);atom.fpos(2)=((1.0/2.0)-x5-y5-z5);atom.fpos(3)=(x5+y5-z5);                     // atom B64
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x5-y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5-y5-z5)");atom.fpos_equation.push_back("(x5+y5-z5)");// atom B64 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B64 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B64

    atom.name="D"; atom.type=3;                                       // atom B65
    atom.fpos(1)=(x5-y5+z5);atom.fpos(2)=(x5+y5-z5);atom.fpos(3)=(-x5+y5+z5);                     // atom B65
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x5-y5+z5)");atom.fpos_equation.push_back("(x5+y5-z5)");atom.fpos_equation.push_back("(-x5+y5+z5)");// atom B65 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B65 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B65

    atom.name="D"; atom.type=3;                                       // atom B66
    atom.fpos(1)=(-x5+y5+z5);atom.fpos(2)=((1.0/2.0)-x5-y5-z5);atom.fpos(3)=(x5-y5+z5);                     // atom B66
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x5+y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5-y5-z5)");atom.fpos_equation.push_back("(x5-y5+z5)");// atom B66 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B66 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B66

    atom.name="D"; atom.type=3;                                       // atom B67
    atom.fpos(1)=((1.0/2.0)-x5-y5-z5);atom.fpos(2)=(-x5+y5+z5);atom.fpos(3)=(x5+y5-z5);                     // atom B67
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x5-y5-z5)");atom.fpos_equation.push_back("(-x5+y5+z5)");atom.fpos_equation.push_back("(x5+y5-z5)");// atom B67 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B67 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B67

    atom.name="D"; atom.type=3;                                       // atom B68
    atom.fpos(1)=(x5+y5-z5);atom.fpos(2)=(x5-y5+z5);atom.fpos(3)=((1.0/2.0)-x5-y5-z5);                     // atom B68
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x5+y5-z5)");atom.fpos_equation.push_back("(x5-y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5-y5-z5)");// atom B68 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B68 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B68

    atom.name="D"; atom.type=3;                                       // atom B69
    atom.fpos(1)=((1.0/2.0)+x5-y5-z5);atom.fpos(2)=((1.0/2.0)-x5+y5-z5);atom.fpos(3)=(x5+y5+z5);                     // atom B69
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x5-y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5+y5-z5)");atom.fpos_equation.push_back("(x5+y5+z5)");// atom B69 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B69 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B69

    atom.name="D"; atom.type=3;                                       // atom B70
    atom.fpos(1)=((1.0/2.0)-x5+y5-z5);atom.fpos(2)=((1.0/2.0)+x5-y5-z5);atom.fpos(3)=((1.0/2.0)-x5-y5+z5);                     // atom B70
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x5+y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5-y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5-y5+z5)");// atom B70 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B70 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B70

    atom.name="D"; atom.type=3;                                       // atom B71
    atom.fpos(1)=((1.0/2.0)-x5-y5+z5);atom.fpos(2)=(x5+y5+z5);atom.fpos(3)=((1.0/2.0)-x5+y5-z5);                     // atom B71
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x5-y5+z5)");atom.fpos_equation.push_back("(x5+y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5+y5-z5)");// atom B71 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B71 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B71

    atom.name="D"; atom.type=3;                                       // atom B72
    atom.fpos(1)=(x5+y5+z5);atom.fpos(2)=((1.0/2.0)-x5-y5+z5);atom.fpos(3)=((1.0/2.0)+x5-y5-z5);                     // atom B72
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x5+y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5-y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5-y5-z5)");// atom B72 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B72 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B72

    atom.name="D"; atom.type=3;                                       // atom B73
    atom.fpos(1)=((1.0/2.0)-x5-y5+z5);atom.fpos(2)=((1.0/2.0)+x5-y5-z5);atom.fpos(3)=(x5+y5+z5);                     // atom B73
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x5-y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5-y5-z5)");atom.fpos_equation.push_back("(x5+y5+z5)");// atom B73 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B73 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B73

    atom.name="D"; atom.type=3;                                       // atom B74
    atom.fpos(1)=(x5+y5+z5);atom.fpos(2)=((1.0/2.0)-x5+y5-z5);atom.fpos(3)=((1.0/2.0)-x5-y5+z5);                     // atom B74
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x5+y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5+y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5-y5+z5)");// atom B74 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B74 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B74

    atom.name="D"; atom.type=3;                                       // atom B75
    atom.fpos(1)=((1.0/2.0)+x5-y5-z5);atom.fpos(2)=((1.0/2.0)-x5-y5+z5);atom.fpos(3)=((1.0/2.0)-x5+y5-z5);                     // atom B75
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x5-y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5-y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5+y5-z5)");// atom B75 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B75 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B75

    atom.name="D"; atom.type=3;                                       // atom B76
    atom.fpos(1)=((1.0/2.0)-x5+y5-z5);atom.fpos(2)=(x5+y5+z5);atom.fpos(3)=((1.0/2.0)+x5-y5-z5);                     // atom B76
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x5+y5-z5)");atom.fpos_equation.push_back("(x5+y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5-y5-z5)");// atom B76 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B76 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B76

    atom.name="D"; atom.type=3;                                       // atom B77
    atom.fpos(1)=((1.0/2.0)-x5+y5-z5);atom.fpos(2)=((1.0/2.0)-x5-y5+z5);atom.fpos(3)=(x5+y5+z5);                     // atom B77
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x5+y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5-y5+z5)");atom.fpos_equation.push_back("(x5+y5+z5)");// atom B77 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B77 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B77

    atom.name="D"; atom.type=3;                                       // atom B78
    atom.fpos(1)=((1.0/2.0)+x5-y5-z5);atom.fpos(2)=(x5+y5+z5);atom.fpos(3)=((1.0/2.0)-x5-y5+z5);                     // atom B78
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x5-y5-z5)");atom.fpos_equation.push_back("(x5+y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5-y5+z5)");// atom B78 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B78 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B78

    atom.name="D"; atom.type=3;                                       // atom B79
    atom.fpos(1)=(x5+y5+z5);atom.fpos(2)=((1.0/2.0)+x5-y5-z5);atom.fpos(3)=((1.0/2.0)-x5+y5-z5);                     // atom B79
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x5+y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5-y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5+y5-z5)");// atom B79 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B79 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B79

    atom.name="D"; atom.type=3;                                       // atom B80
    atom.fpos(1)=((1.0/2.0)-x5-y5+z5);atom.fpos(2)=((1.0/2.0)-x5+y5-z5);atom.fpos(3)=((1.0/2.0)+x5-y5-z5);                     // atom B80
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x5-y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5+y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5-y5-z5)");// atom B80 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B80 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B80

    atom.name="D"; atom.type=3;                                       // atom B81
    atom.fpos(1)=(x5-y5-z5);atom.fpos(2)=(-x5+y5-z5);atom.fpos(3)=(-x5-y5+z5);                     // atom B81
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x5-y5-z5)");atom.fpos_equation.push_back("(-x5+y5-z5)");atom.fpos_equation.push_back("(-x5-y5+z5)");// atom B81 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B81 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B81

    atom.name="D"; atom.type=3;                                       // atom B82
    atom.fpos(1)=(-x5+y5-z5);atom.fpos(2)=(x5-y5-z5);atom.fpos(3)=((1.0/2.0)+x5+y5+z5);                     // atom B82
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x5+y5-z5)");atom.fpos_equation.push_back("(x5-y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5+y5+z5)");// atom B82 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B82 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B82

    atom.name="D"; atom.type=3;                                       // atom B83
    atom.fpos(1)=(-x5-y5+z5);atom.fpos(2)=((1.0/2.0)+x5+y5+z5);atom.fpos(3)=(x5-y5-z5);                     // atom B83
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x5-y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5+y5+z5)");atom.fpos_equation.push_back("(x5-y5-z5)");// atom B83 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B83 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B83

    atom.name="D"; atom.type=3;                                       // atom B84
    atom.fpos(1)=((1.0/2.0)+x5+y5+z5);atom.fpos(2)=(-x5-y5+z5);atom.fpos(3)=(-x5+y5-z5);                     // atom B84
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x5+y5+z5)");atom.fpos_equation.push_back("(-x5-y5+z5)");atom.fpos_equation.push_back("(-x5+y5-z5)");// atom B84 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B84 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B84

    atom.name="D"; atom.type=3;                                       // atom B85
    atom.fpos(1)=(-x5-y5+z5);atom.fpos(2)=(x5-y5-z5);atom.fpos(3)=(-x5+y5-z5);                     // atom B85
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x5-y5+z5)");atom.fpos_equation.push_back("(x5-y5-z5)");atom.fpos_equation.push_back("(-x5+y5-z5)");// atom B85 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B85 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B85

    atom.name="D"; atom.type=3;                                       // atom B86
    atom.fpos(1)=((1.0/2.0)+x5+y5+z5);atom.fpos(2)=(-x5+y5-z5);atom.fpos(3)=(x5-y5-z5);                     // atom B86
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x5+y5+z5)");atom.fpos_equation.push_back("(-x5+y5-z5)");atom.fpos_equation.push_back("(x5-y5-z5)");// atom B86 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B86 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B86

    atom.name="D"; atom.type=3;                                       // atom B87
    atom.fpos(1)=(x5-y5-z5);atom.fpos(2)=(-x5-y5+z5);atom.fpos(3)=((1.0/2.0)+x5+y5+z5);                     // atom B87
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x5-y5-z5)");atom.fpos_equation.push_back("(-x5-y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5+y5+z5)");// atom B87 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B87 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B87

    atom.name="D"; atom.type=3;                                       // atom B88
    atom.fpos(1)=(-x5+y5-z5);atom.fpos(2)=((1.0/2.0)+x5+y5+z5);atom.fpos(3)=(-x5-y5+z5);                     // atom B88
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x5+y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5+y5+z5)");atom.fpos_equation.push_back("(-x5-y5+z5)");// atom B88 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B88 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B88

    atom.name="D"; atom.type=3;                                       // atom B89
    atom.fpos(1)=(-x5+y5-z5);atom.fpos(2)=(-x5-y5+z5);atom.fpos(3)=(x5-y5-z5);                     // atom B89
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x5+y5-z5)");atom.fpos_equation.push_back("(-x5-y5+z5)");atom.fpos_equation.push_back("(x5-y5-z5)");// atom B89 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B89 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B89

    atom.name="D"; atom.type=3;                                       // atom B90
    atom.fpos(1)=(x5-y5-z5);atom.fpos(2)=((1.0/2.0)+x5+y5+z5);atom.fpos(3)=(-x5+y5-z5);                     // atom B90
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x5-y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5+y5+z5)");atom.fpos_equation.push_back("(-x5+y5-z5)");// atom B90 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B90 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B90

    atom.name="D"; atom.type=3;                                       // atom B91
    atom.fpos(1)=((1.0/2.0)+x5+y5+z5);atom.fpos(2)=(x5-y5-z5);atom.fpos(3)=(-x5-y5+z5);                     // atom B91
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x5+y5+z5)");atom.fpos_equation.push_back("(x5-y5-z5)");atom.fpos_equation.push_back("(-x5-y5+z5)");// atom B91 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B91 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B91

    atom.name="D"; atom.type=3;                                       // atom B92
    atom.fpos(1)=(-x5-y5+z5);atom.fpos(2)=(-x5+y5-z5);atom.fpos(3)=((1.0/2.0)+x5+y5+z5);                     // atom B92
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x5-y5+z5)");atom.fpos_equation.push_back("(-x5+y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5+y5+z5)");// atom B92 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B92 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B92

    atom.name="D"; atom.type=3;                                       // atom B93
    atom.fpos(1)=((1.0/2.0)-x5+y5+z5);atom.fpos(2)=((1.0/2.0)+x5-y5+z5);atom.fpos(3)=(-x5-y5-z5);                     // atom B93
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x5+y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5-y5+z5)");atom.fpos_equation.push_back("(-x5-y5-z5)");// atom B93 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B93 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B93

    atom.name="D"; atom.type=3;                                       // atom B94
    atom.fpos(1)=((1.0/2.0)+x5-y5+z5);atom.fpos(2)=((1.0/2.0)-x5+y5+z5);atom.fpos(3)=((1.0/2.0)+x5+y5-z5);                     // atom B94
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x5-y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5+y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5+y5-z5)");// atom B94 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B94 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B94

    atom.name="D"; atom.type=3;                                       // atom B95
    atom.fpos(1)=((1.0/2.0)+x5+y5-z5);atom.fpos(2)=(-x5-y5-z5);atom.fpos(3)=((1.0/2.0)+x5-y5+z5);                     // atom B95
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x5+y5-z5)");atom.fpos_equation.push_back("(-x5-y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5-y5+z5)");// atom B95 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B95 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B95

    atom.name="D"; atom.type=3;                                       // atom B96
    atom.fpos(1)=(-x5-y5-z5);atom.fpos(2)=((1.0/2.0)+x5+y5-z5);atom.fpos(3)=((1.0/2.0)-x5+y5+z5);                     // atom B96
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x5-y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5+y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5+y5+z5)");// atom B96 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B96 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B96

    atom.name="D"; atom.type=3;                                       // atom B97
    atom.fpos(1)=((1.0/2.0)+x5+y5-z5);atom.fpos(2)=((1.0/2.0)-x5+y5+z5);atom.fpos(3)=(-x5-y5-z5);                     // atom B97
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x5+y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5+y5+z5)");atom.fpos_equation.push_back("(-x5-y5-z5)");// atom B97 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B97 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B97

    atom.name="D"; atom.type=3;                                       // atom B98
    atom.fpos(1)=(-x5-y5-z5);atom.fpos(2)=((1.0/2.0)+x5-y5+z5);atom.fpos(3)=((1.0/2.0)+x5+y5-z5);                     // atom B98
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x5-y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5-y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5+y5-z5)");// atom B98 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B98 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B98

    atom.name="D"; atom.type=3;                                       // atom B99
    atom.fpos(1)=((1.0/2.0)-x5+y5+z5);atom.fpos(2)=((1.0/2.0)+x5+y5-z5);atom.fpos(3)=((1.0/2.0)+x5-y5+z5);                     // atom B99
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x5+y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5+y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5-y5+z5)");// atom B99 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B99 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B99

    atom.name="D"; atom.type=3;                                       // atom B100
    atom.fpos(1)=((1.0/2.0)+x5-y5+z5);atom.fpos(2)=(-x5-y5-z5);atom.fpos(3)=((1.0/2.0)-x5+y5+z5);                     // atom B100
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x5-y5+z5)");atom.fpos_equation.push_back("(-x5-y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5+y5+z5)");// atom B100 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B100 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B100

    atom.name="D"; atom.type=3;                                       // atom B101
    atom.fpos(1)=((1.0/2.0)+x5-y5+z5);atom.fpos(2)=((1.0/2.0)+x5+y5-z5);atom.fpos(3)=(-x5-y5-z5);                     // atom B101
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x5-y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5+y5-z5)");atom.fpos_equation.push_back("(-x5-y5-z5)");// atom B101 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B101 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B101

    atom.name="D"; atom.type=3;                                       // atom B102
    atom.fpos(1)=((1.0/2.0)-x5+y5+z5);atom.fpos(2)=(-x5-y5-z5);atom.fpos(3)=((1.0/2.0)+x5+y5-z5);                     // atom B102
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x5+y5+z5)");atom.fpos_equation.push_back("(-x5-y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5+y5-z5)");// atom B102 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B102 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B102

    atom.name="D"; atom.type=3;                                       // atom B103
    atom.fpos(1)=(-x5-y5-z5);atom.fpos(2)=((1.0/2.0)-x5+y5+z5);atom.fpos(3)=((1.0/2.0)+x5-y5+z5);                     // atom B103
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x5-y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5+y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5-y5+z5)");// atom B103 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B103 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B103

    atom.name="D"; atom.type=3;                                       // atom B104
    atom.fpos(1)=((1.0/2.0)+x5+y5-z5);atom.fpos(2)=((1.0/2.0)+x5-y5+z5);atom.fpos(3)=((1.0/2.0)-x5+y5+z5);                     // atom B104
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x5+y5-z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5-y5+z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5+y5+z5)");// atom B104 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B104 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B104


    return str.atoms.size();  
  }
} // namespace anrl

namespace anrl {
  uint WebANRL_A5BCD6_cF416_228_eg_c_b_h(stringstream& web,bool LDEBUG) {
#ifndef _ANRL_NOWEB_
#endif

    if(LDEBUG) {cerr << "anrl:: WebANRL_A5BCD6_cF416_228_eg_c_b_h: web.str().size()=" << web.str().size() << endl;}

    return web.str().size();
  }
} // namespace anrl

#endif // AFLOW_REMOVE_GREP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

