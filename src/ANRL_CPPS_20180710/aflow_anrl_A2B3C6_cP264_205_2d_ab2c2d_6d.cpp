// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - David Hicks - 2018
// FILE "ANRL/aflow_anrl_A2B3C6_cP264_205_2d_ab2c2d_6d.cpp"

#ifndef _AFLOW_ANRL_A2B3C6_cP264_205_2d_ab2c2d_6d_CPP // AFLOW_REMOVE_GREP
#define _AFLOW_ANRL_A2B3C6_cP264_205_2d_ab2c2d_6d_CPP // AFLOW_REMOVE_GREP
#include "../aflow.h" // AFLOW_REMOVE_GREP

namespace anrl {
  uint WebANRL_A2B3C6_cP264_205_2d_ab2c2d_6d(stringstream &web,bool LDEBUG);
}

namespace anrl {
  uint PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d(stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG) {
    // system A2B3C6_cP264_205_2d_ab2c2d_6d

    if(XHOST.vflag_control.flag("WWW")) {
      WebANRL_A2B3C6_cP264_205_2d_ab2c2d_6d(web,LDEBUG); // PLUG WEB STUFF
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

    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: FOUND" << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: label=" << label << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: nspecies=" << nspecies << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: natoms=" << natoms << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: spacegroup=" << spacegroup << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: nunderscores=" << nunderscores << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: nparameters=" <<  nparameters << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: Pearson_symbol=" << Pearson_symbol << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: params=" << params << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: Strukturbericht=" << Strukturbericht << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: prototype=" << prototype << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: dialect=" << dialect << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: vparameters.size()=" << vparameters.size() << endl;}

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
    double a=vparameters.at(i++);                  if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: a=" << a << endl;}

    double x3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: x3=" << x3 << endl;}
    double x4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: x4=" << x4 << endl;}
    double x5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: x5=" << x5 << endl;}
    double y5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: y5=" << y5 << endl;}
    double z5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: z5=" << z5 << endl;}
    double x6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: x6=" << x6 << endl;}
    double y6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: y6=" << y6 << endl;}
    double z6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: z6=" << z6 << endl;}
    double x7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: x7=" << x7 << endl;}
    double y7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: y7=" << y7 << endl;}
    double z7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: z7=" << z7 << endl;}
    double x8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: x8=" << x8 << endl;}
    double y8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: y8=" << y8 << endl;}
    double z8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: z8=" << z8 << endl;}
    double x9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: x9=" << x9 << endl;}
    double y9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: y9=" << y9 << endl;}
    double z9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: z9=" << z9 << endl;}
    double x10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: x10=" << x10 << endl;}
    double y10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: y10=" << y10 << endl;}
    double z10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: z10=" << z10 << endl;}
    double x11=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: x11=" << x11 << endl;}
    double y11=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: y11=" << y11 << endl;}
    double z11=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: z11=" << z11 << endl;}
    double x12=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: x12=" << x12 << endl;}
    double y12=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: y12=" << y12 << endl;}
    double z12=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: z12=" << z12 << endl;}
    double x13=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: x13=" << x13 << endl;}
    double y13=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: y13=" << y13 << endl;}
    double z13=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: z13=" << z13 << endl;}
    double x14=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: x14=" << x14 << endl;}
    double y14=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: y14=" << y14 << endl;}
    double z14=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: z14=" << z14 << endl;}

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

    atom.name="A"; atom.type=0;                                       // atom B25
    atom.fpos(1)=x5;atom.fpos(2)=y5;atom.fpos(3)=z5;                     // atom B25
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("y5");atom.fpos_equation.push_back("z5");// atom B25 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B25 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B25

    atom.name="A"; atom.type=0;                                       // atom B26
    atom.fpos(1)=((1.0/2.0)-x5);atom.fpos(2)=-y5;atom.fpos(3)=((1.0/2.0)+z5);                     // atom B26
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x5)");atom.fpos_equation.push_back("-y5");atom.fpos_equation.push_back("((1.0/2.0)+z5)");// atom B26 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B26 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B26

    atom.name="A"; atom.type=0;                                       // atom B27
    atom.fpos(1)=-x5;atom.fpos(2)=((1.0/2.0)+y5);atom.fpos(3)=((1.0/2.0)-z5);                     // atom B27
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("((1.0/2.0)+y5)");atom.fpos_equation.push_back("((1.0/2.0)-z5)");// atom B27 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B27 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B27

    atom.name="A"; atom.type=0;                                       // atom B28
    atom.fpos(1)=((1.0/2.0)+x5);atom.fpos(2)=((1.0/2.0)-y5);atom.fpos(3)=-z5;                     // atom B28
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x5)");atom.fpos_equation.push_back("((1.0/2.0)-y5)");atom.fpos_equation.push_back("-z5");// atom B28 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B28 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B28

    atom.name="A"; atom.type=0;                                       // atom B29
    atom.fpos(1)=z5;atom.fpos(2)=x5;atom.fpos(3)=y5;                     // atom B29
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z5");atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("y5");// atom B29 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B29 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B29

    atom.name="A"; atom.type=0;                                       // atom B30
    atom.fpos(1)=((1.0/2.0)+z5);atom.fpos(2)=((1.0/2.0)-x5);atom.fpos(3)=-y5;                     // atom B30
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5)");atom.fpos_equation.push_back("-y5");// atom B30 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B30 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B30

    atom.name="A"; atom.type=0;                                       // atom B31
    atom.fpos(1)=((1.0/2.0)-z5);atom.fpos(2)=-x5;atom.fpos(3)=((1.0/2.0)+y5);                     // atom B31
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z5)");atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("((1.0/2.0)+y5)");// atom B31 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B31 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B31

    atom.name="A"; atom.type=0;                                       // atom B32
    atom.fpos(1)=-z5;atom.fpos(2)=((1.0/2.0)+x5);atom.fpos(3)=((1.0/2.0)-y5);                     // atom B32
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z5");atom.fpos_equation.push_back("((1.0/2.0)+x5)");atom.fpos_equation.push_back("((1.0/2.0)-y5)");// atom B32 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B32 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B32

    atom.name="A"; atom.type=0;                                       // atom B33
    atom.fpos(1)=y5;atom.fpos(2)=z5;atom.fpos(3)=x5;                     // atom B33
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y5");atom.fpos_equation.push_back("z5");atom.fpos_equation.push_back("x5");// atom B33 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B33 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B33

    atom.name="A"; atom.type=0;                                       // atom B34
    atom.fpos(1)=-y5;atom.fpos(2)=((1.0/2.0)+z5);atom.fpos(3)=((1.0/2.0)-x5);                     // atom B34
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y5");atom.fpos_equation.push_back("((1.0/2.0)+z5)");atom.fpos_equation.push_back("((1.0/2.0)-x5)");// atom B34 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B34 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B34

    atom.name="A"; atom.type=0;                                       // atom B35
    atom.fpos(1)=((1.0/2.0)+y5);atom.fpos(2)=((1.0/2.0)-z5);atom.fpos(3)=-x5;                     // atom B35
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y5)");atom.fpos_equation.push_back("((1.0/2.0)-z5)");atom.fpos_equation.push_back("-x5");// atom B35 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B35 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B35

    atom.name="A"; atom.type=0;                                       // atom B36
    atom.fpos(1)=((1.0/2.0)-y5);atom.fpos(2)=-z5;atom.fpos(3)=((1.0/2.0)+x5);                     // atom B36
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y5)");atom.fpos_equation.push_back("-z5");atom.fpos_equation.push_back("((1.0/2.0)+x5)");// atom B36 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B36 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B36

    atom.name="A"; atom.type=0;                                       // atom B37
    atom.fpos(1)=-x5;atom.fpos(2)=-y5;atom.fpos(3)=-z5;                     // atom B37
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("-y5");atom.fpos_equation.push_back("-z5");// atom B37 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B37 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B37

    atom.name="A"; atom.type=0;                                       // atom B38
    atom.fpos(1)=((1.0/2.0)+x5);atom.fpos(2)=y5;atom.fpos(3)=((1.0/2.0)-z5);                     // atom B38
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x5)");atom.fpos_equation.push_back("y5");atom.fpos_equation.push_back("((1.0/2.0)-z5)");// atom B38 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B38 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B38

    atom.name="A"; atom.type=0;                                       // atom B39
    atom.fpos(1)=x5;atom.fpos(2)=((1.0/2.0)-y5);atom.fpos(3)=((1.0/2.0)+z5);                     // atom B39
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("((1.0/2.0)-y5)");atom.fpos_equation.push_back("((1.0/2.0)+z5)");// atom B39 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B39 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B39

    atom.name="A"; atom.type=0;                                       // atom B40
    atom.fpos(1)=((1.0/2.0)-x5);atom.fpos(2)=((1.0/2.0)+y5);atom.fpos(3)=z5;                     // atom B40
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x5)");atom.fpos_equation.push_back("((1.0/2.0)+y5)");atom.fpos_equation.push_back("z5");// atom B40 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B40 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B40

    atom.name="A"; atom.type=0;                                       // atom B41
    atom.fpos(1)=-z5;atom.fpos(2)=-x5;atom.fpos(3)=-y5;                     // atom B41
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z5");atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("-y5");// atom B41 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B41 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B41

    atom.name="A"; atom.type=0;                                       // atom B42
    atom.fpos(1)=((1.0/2.0)-z5);atom.fpos(2)=((1.0/2.0)+x5);atom.fpos(3)=y5;                     // atom B42
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5)");atom.fpos_equation.push_back("y5");// atom B42 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B42 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B42

    atom.name="A"; atom.type=0;                                       // atom B43
    atom.fpos(1)=((1.0/2.0)+z5);atom.fpos(2)=x5;atom.fpos(3)=((1.0/2.0)-y5);                     // atom B43
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z5)");atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("((1.0/2.0)-y5)");// atom B43 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B43 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B43

    atom.name="A"; atom.type=0;                                       // atom B44
    atom.fpos(1)=z5;atom.fpos(2)=((1.0/2.0)-x5);atom.fpos(3)=((1.0/2.0)+y5);                     // atom B44
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z5");atom.fpos_equation.push_back("((1.0/2.0)-x5)");atom.fpos_equation.push_back("((1.0/2.0)+y5)");// atom B44 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B44 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B44

    atom.name="A"; atom.type=0;                                       // atom B45
    atom.fpos(1)=-y5;atom.fpos(2)=-z5;atom.fpos(3)=-x5;                     // atom B45
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y5");atom.fpos_equation.push_back("-z5");atom.fpos_equation.push_back("-x5");// atom B45 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B45 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B45

    atom.name="A"; atom.type=0;                                       // atom B46
    atom.fpos(1)=y5;atom.fpos(2)=((1.0/2.0)-z5);atom.fpos(3)=((1.0/2.0)+x5);                     // atom B46
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y5");atom.fpos_equation.push_back("((1.0/2.0)-z5)");atom.fpos_equation.push_back("((1.0/2.0)+x5)");// atom B46 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B46 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B46

    atom.name="A"; atom.type=0;                                       // atom B47
    atom.fpos(1)=((1.0/2.0)-y5);atom.fpos(2)=((1.0/2.0)+z5);atom.fpos(3)=x5;                     // atom B47
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y5)");atom.fpos_equation.push_back("((1.0/2.0)+z5)");atom.fpos_equation.push_back("x5");// atom B47 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B47 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B47

    atom.name="A"; atom.type=0;                                       // atom B48
    atom.fpos(1)=((1.0/2.0)+y5);atom.fpos(2)=z5;atom.fpos(3)=((1.0/2.0)-x5);                     // atom B48
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y5)");atom.fpos_equation.push_back("z5");atom.fpos_equation.push_back("((1.0/2.0)-x5)");// atom B48 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B48 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B48

    atom.name="A"; atom.type=0;                                       // atom B49
    atom.fpos(1)=x6;atom.fpos(2)=y6;atom.fpos(3)=z6;                     // atom B49
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("y6");atom.fpos_equation.push_back("z6");// atom B49 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B49 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B49

    atom.name="A"; atom.type=0;                                       // atom B50
    atom.fpos(1)=((1.0/2.0)-x6);atom.fpos(2)=-y6;atom.fpos(3)=((1.0/2.0)+z6);                     // atom B50
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x6)");atom.fpos_equation.push_back("-y6");atom.fpos_equation.push_back("((1.0/2.0)+z6)");// atom B50 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B50 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B50

    atom.name="A"; atom.type=0;                                       // atom B51
    atom.fpos(1)=-x6;atom.fpos(2)=((1.0/2.0)+y6);atom.fpos(3)=((1.0/2.0)-z6);                     // atom B51
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("((1.0/2.0)+y6)");atom.fpos_equation.push_back("((1.0/2.0)-z6)");// atom B51 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B51 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B51

    atom.name="A"; atom.type=0;                                       // atom B52
    atom.fpos(1)=((1.0/2.0)+x6);atom.fpos(2)=((1.0/2.0)-y6);atom.fpos(3)=-z6;                     // atom B52
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x6)");atom.fpos_equation.push_back("((1.0/2.0)-y6)");atom.fpos_equation.push_back("-z6");// atom B52 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B52 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B52

    atom.name="A"; atom.type=0;                                       // atom B53
    atom.fpos(1)=z6;atom.fpos(2)=x6;atom.fpos(3)=y6;                     // atom B53
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z6");atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("y6");// atom B53 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B53 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B53

    atom.name="A"; atom.type=0;                                       // atom B54
    atom.fpos(1)=((1.0/2.0)+z6);atom.fpos(2)=((1.0/2.0)-x6);atom.fpos(3)=-y6;                     // atom B54
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z6)");atom.fpos_equation.push_back("((1.0/2.0)-x6)");atom.fpos_equation.push_back("-y6");// atom B54 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B54 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B54

    atom.name="A"; atom.type=0;                                       // atom B55
    atom.fpos(1)=((1.0/2.0)-z6);atom.fpos(2)=-x6;atom.fpos(3)=((1.0/2.0)+y6);                     // atom B55
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z6)");atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("((1.0/2.0)+y6)");// atom B55 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B55 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B55

    atom.name="A"; atom.type=0;                                       // atom B56
    atom.fpos(1)=-z6;atom.fpos(2)=((1.0/2.0)+x6);atom.fpos(3)=((1.0/2.0)-y6);                     // atom B56
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z6");atom.fpos_equation.push_back("((1.0/2.0)+x6)");atom.fpos_equation.push_back("((1.0/2.0)-y6)");// atom B56 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B56 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B56

    atom.name="A"; atom.type=0;                                       // atom B57
    atom.fpos(1)=y6;atom.fpos(2)=z6;atom.fpos(3)=x6;                     // atom B57
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y6");atom.fpos_equation.push_back("z6");atom.fpos_equation.push_back("x6");// atom B57 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B57 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B57

    atom.name="A"; atom.type=0;                                       // atom B58
    atom.fpos(1)=-y6;atom.fpos(2)=((1.0/2.0)+z6);atom.fpos(3)=((1.0/2.0)-x6);                     // atom B58
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y6");atom.fpos_equation.push_back("((1.0/2.0)+z6)");atom.fpos_equation.push_back("((1.0/2.0)-x6)");// atom B58 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B58 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B58

    atom.name="A"; atom.type=0;                                       // atom B59
    atom.fpos(1)=((1.0/2.0)+y6);atom.fpos(2)=((1.0/2.0)-z6);atom.fpos(3)=-x6;                     // atom B59
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y6)");atom.fpos_equation.push_back("((1.0/2.0)-z6)");atom.fpos_equation.push_back("-x6");// atom B59 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B59 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B59

    atom.name="A"; atom.type=0;                                       // atom B60
    atom.fpos(1)=((1.0/2.0)-y6);atom.fpos(2)=-z6;atom.fpos(3)=((1.0/2.0)+x6);                     // atom B60
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y6)");atom.fpos_equation.push_back("-z6");atom.fpos_equation.push_back("((1.0/2.0)+x6)");// atom B60 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B60 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B60

    atom.name="A"; atom.type=0;                                       // atom B61
    atom.fpos(1)=-x6;atom.fpos(2)=-y6;atom.fpos(3)=-z6;                     // atom B61
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("-y6");atom.fpos_equation.push_back("-z6");// atom B61 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B61 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B61

    atom.name="A"; atom.type=0;                                       // atom B62
    atom.fpos(1)=((1.0/2.0)+x6);atom.fpos(2)=y6;atom.fpos(3)=((1.0/2.0)-z6);                     // atom B62
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x6)");atom.fpos_equation.push_back("y6");atom.fpos_equation.push_back("((1.0/2.0)-z6)");// atom B62 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B62 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B62

    atom.name="A"; atom.type=0;                                       // atom B63
    atom.fpos(1)=x6;atom.fpos(2)=((1.0/2.0)-y6);atom.fpos(3)=((1.0/2.0)+z6);                     // atom B63
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("((1.0/2.0)-y6)");atom.fpos_equation.push_back("((1.0/2.0)+z6)");// atom B63 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B63 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B63

    atom.name="A"; atom.type=0;                                       // atom B64
    atom.fpos(1)=((1.0/2.0)-x6);atom.fpos(2)=((1.0/2.0)+y6);atom.fpos(3)=z6;                     // atom B64
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x6)");atom.fpos_equation.push_back("((1.0/2.0)+y6)");atom.fpos_equation.push_back("z6");// atom B64 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B64 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B64

    atom.name="A"; atom.type=0;                                       // atom B65
    atom.fpos(1)=-z6;atom.fpos(2)=-x6;atom.fpos(3)=-y6;                     // atom B65
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z6");atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("-y6");// atom B65 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B65 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B65

    atom.name="A"; atom.type=0;                                       // atom B66
    atom.fpos(1)=((1.0/2.0)-z6);atom.fpos(2)=((1.0/2.0)+x6);atom.fpos(3)=y6;                     // atom B66
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z6)");atom.fpos_equation.push_back("((1.0/2.0)+x6)");atom.fpos_equation.push_back("y6");// atom B66 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B66 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B66

    atom.name="A"; atom.type=0;                                       // atom B67
    atom.fpos(1)=((1.0/2.0)+z6);atom.fpos(2)=x6;atom.fpos(3)=((1.0/2.0)-y6);                     // atom B67
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z6)");atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("((1.0/2.0)-y6)");// atom B67 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B67 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B67

    atom.name="A"; atom.type=0;                                       // atom B68
    atom.fpos(1)=z6;atom.fpos(2)=((1.0/2.0)-x6);atom.fpos(3)=((1.0/2.0)+y6);                     // atom B68
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z6");atom.fpos_equation.push_back("((1.0/2.0)-x6)");atom.fpos_equation.push_back("((1.0/2.0)+y6)");// atom B68 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B68 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B68

    atom.name="A"; atom.type=0;                                       // atom B69
    atom.fpos(1)=-y6;atom.fpos(2)=-z6;atom.fpos(3)=-x6;                     // atom B69
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y6");atom.fpos_equation.push_back("-z6");atom.fpos_equation.push_back("-x6");// atom B69 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B69 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B69

    atom.name="A"; atom.type=0;                                       // atom B70
    atom.fpos(1)=y6;atom.fpos(2)=((1.0/2.0)-z6);atom.fpos(3)=((1.0/2.0)+x6);                     // atom B70
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y6");atom.fpos_equation.push_back("((1.0/2.0)-z6)");atom.fpos_equation.push_back("((1.0/2.0)+x6)");// atom B70 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B70 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B70

    atom.name="A"; atom.type=0;                                       // atom B71
    atom.fpos(1)=((1.0/2.0)-y6);atom.fpos(2)=((1.0/2.0)+z6);atom.fpos(3)=x6;                     // atom B71
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y6)");atom.fpos_equation.push_back("((1.0/2.0)+z6)");atom.fpos_equation.push_back("x6");// atom B71 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B71 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B71

    atom.name="A"; atom.type=0;                                       // atom B72
    atom.fpos(1)=((1.0/2.0)+y6);atom.fpos(2)=z6;atom.fpos(3)=((1.0/2.0)-x6);                     // atom B72
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y6)");atom.fpos_equation.push_back("z6");atom.fpos_equation.push_back("((1.0/2.0)-x6)");// atom B72 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B72 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B72

    atom.name="B"; atom.type=1;                                       // atom B1
    atom.fpos(1)=0.0;atom.fpos(2)=0.0;atom.fpos(3)=0.0;                     // atom B1
    atom.fpos_equation.clear();atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("0.0");// atom B1 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B1 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B1

    atom.name="B"; atom.type=1;                                       // atom B2
    atom.fpos(1)=(1.0/2.0);atom.fpos(2)=0.0;atom.fpos(3)=(1.0/2.0);                     // atom B2
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("(1.0/2.0)");// atom B2 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B2 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B2

    atom.name="B"; atom.type=1;                                       // atom B3
    atom.fpos(1)=0.0;atom.fpos(2)=(1.0/2.0);atom.fpos(3)=(1.0/2.0);                     // atom B3
    atom.fpos_equation.clear();atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("(1.0/2.0)");// atom B3 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B3 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B3

    atom.name="B"; atom.type=1;                                       // atom B4
    atom.fpos(1)=(1.0/2.0);atom.fpos(2)=(1.0/2.0);atom.fpos(3)=0.0;                     // atom B4
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("0.0");// atom B4 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B4 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B4

    atom.name="B"; atom.type=1;                                       // atom B5
    atom.fpos(1)=(1.0/2.0);atom.fpos(2)=(1.0/2.0);atom.fpos(3)=(1.0/2.0);                     // atom B5
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("(1.0/2.0)");// atom B5 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B5 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B5

    atom.name="B"; atom.type=1;                                       // atom B6
    atom.fpos(1)=0.0;atom.fpos(2)=(1.0/2.0);atom.fpos(3)=0.0;                     // atom B6
    atom.fpos_equation.clear();atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("0.0");// atom B6 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B6 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B6

    atom.name="B"; atom.type=1;                                       // atom B7
    atom.fpos(1)=(1.0/2.0);atom.fpos(2)=0.0;atom.fpos(3)=0.0;                     // atom B7
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("0.0");// atom B7 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B7 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B7

    atom.name="B"; atom.type=1;                                       // atom B8
    atom.fpos(1)=0.0;atom.fpos(2)=0.0;atom.fpos(3)=(1.0/2.0);                     // atom B8
    atom.fpos_equation.clear();atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("(1.0/2.0)");// atom B8 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B8 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B8

    atom.name="B"; atom.type=1;                                       // atom B9
    atom.fpos(1)=x3;atom.fpos(2)=x3;atom.fpos(3)=x3;                     // atom B9
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("x3");// atom B9 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B9 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B9

    atom.name="B"; atom.type=1;                                       // atom B10
    atom.fpos(1)=((1.0/2.0)-x3);atom.fpos(2)=-x3;atom.fpos(3)=((1.0/2.0)+x3);                     // atom B10
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x3)");atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("((1.0/2.0)+x3)");// atom B10 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B10 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B10

    atom.name="B"; atom.type=1;                                       // atom B11
    atom.fpos(1)=-x3;atom.fpos(2)=((1.0/2.0)+x3);atom.fpos(3)=((1.0/2.0)-x3);                     // atom B11
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("((1.0/2.0)+x3)");atom.fpos_equation.push_back("((1.0/2.0)-x3)");// atom B11 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B11 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B11

    atom.name="B"; atom.type=1;                                       // atom B12
    atom.fpos(1)=((1.0/2.0)+x3);atom.fpos(2)=((1.0/2.0)-x3);atom.fpos(3)=-x3;                     // atom B12
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x3)");atom.fpos_equation.push_back("((1.0/2.0)-x3)");atom.fpos_equation.push_back("-x3");// atom B12 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B12 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B12

    atom.name="B"; atom.type=1;                                       // atom B13
    atom.fpos(1)=-x3;atom.fpos(2)=-x3;atom.fpos(3)=-x3;                     // atom B13
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("-x3");// atom B13 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B13 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B13

    atom.name="B"; atom.type=1;                                       // atom B14
    atom.fpos(1)=((1.0/2.0)+x3);atom.fpos(2)=x3;atom.fpos(3)=((1.0/2.0)-x3);                     // atom B14
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x3)");atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("((1.0/2.0)-x3)");// atom B14 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B14 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B14

    atom.name="B"; atom.type=1;                                       // atom B15
    atom.fpos(1)=x3;atom.fpos(2)=((1.0/2.0)-x3);atom.fpos(3)=((1.0/2.0)+x3);                     // atom B15
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("((1.0/2.0)-x3)");atom.fpos_equation.push_back("((1.0/2.0)+x3)");// atom B15 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B15 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B15

    atom.name="B"; atom.type=1;                                       // atom B16
    atom.fpos(1)=((1.0/2.0)-x3);atom.fpos(2)=((1.0/2.0)+x3);atom.fpos(3)=x3;                     // atom B16
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x3)");atom.fpos_equation.push_back("((1.0/2.0)+x3)");atom.fpos_equation.push_back("x3");// atom B16 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B16 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B16

    atom.name="B"; atom.type=1;                                       // atom B17
    atom.fpos(1)=x4;atom.fpos(2)=x4;atom.fpos(3)=x4;                     // atom B17
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x4");atom.fpos_equation.push_back("x4");atom.fpos_equation.push_back("x4");// atom B17 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B17 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B17

    atom.name="B"; atom.type=1;                                       // atom B18
    atom.fpos(1)=((1.0/2.0)-x4);atom.fpos(2)=-x4;atom.fpos(3)=((1.0/2.0)+x4);                     // atom B18
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x4)");atom.fpos_equation.push_back("-x4");atom.fpos_equation.push_back("((1.0/2.0)+x4)");// atom B18 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B18 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B18

    atom.name="B"; atom.type=1;                                       // atom B19
    atom.fpos(1)=-x4;atom.fpos(2)=((1.0/2.0)+x4);atom.fpos(3)=((1.0/2.0)-x4);                     // atom B19
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x4");atom.fpos_equation.push_back("((1.0/2.0)+x4)");atom.fpos_equation.push_back("((1.0/2.0)-x4)");// atom B19 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B19 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B19

    atom.name="B"; atom.type=1;                                       // atom B20
    atom.fpos(1)=((1.0/2.0)+x4);atom.fpos(2)=((1.0/2.0)-x4);atom.fpos(3)=-x4;                     // atom B20
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x4)");atom.fpos_equation.push_back("((1.0/2.0)-x4)");atom.fpos_equation.push_back("-x4");// atom B20 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B20 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B20

    atom.name="B"; atom.type=1;                                       // atom B21
    atom.fpos(1)=-x4;atom.fpos(2)=-x4;atom.fpos(3)=-x4;                     // atom B21
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x4");atom.fpos_equation.push_back("-x4");atom.fpos_equation.push_back("-x4");// atom B21 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B21 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B21

    atom.name="B"; atom.type=1;                                       // atom B22
    atom.fpos(1)=((1.0/2.0)+x4);atom.fpos(2)=x4;atom.fpos(3)=((1.0/2.0)-x4);                     // atom B22
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x4)");atom.fpos_equation.push_back("x4");atom.fpos_equation.push_back("((1.0/2.0)-x4)");// atom B22 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B22 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B22

    atom.name="B"; atom.type=1;                                       // atom B23
    atom.fpos(1)=x4;atom.fpos(2)=((1.0/2.0)-x4);atom.fpos(3)=((1.0/2.0)+x4);                     // atom B23
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x4");atom.fpos_equation.push_back("((1.0/2.0)-x4)");atom.fpos_equation.push_back("((1.0/2.0)+x4)");// atom B23 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B23 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B23

    atom.name="B"; atom.type=1;                                       // atom B24
    atom.fpos(1)=((1.0/2.0)-x4);atom.fpos(2)=((1.0/2.0)+x4);atom.fpos(3)=x4;                     // atom B24
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x4)");atom.fpos_equation.push_back("((1.0/2.0)+x4)");atom.fpos_equation.push_back("x4");// atom B24 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B24 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B24

    atom.name="B"; atom.type=1;                                       // atom B73
    atom.fpos(1)=x7;atom.fpos(2)=y7;atom.fpos(3)=z7;                     // atom B73
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x7");atom.fpos_equation.push_back("y7");atom.fpos_equation.push_back("z7");// atom B73 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B73 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B73

    atom.name="B"; atom.type=1;                                       // atom B74
    atom.fpos(1)=((1.0/2.0)-x7);atom.fpos(2)=-y7;atom.fpos(3)=((1.0/2.0)+z7);                     // atom B74
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x7)");atom.fpos_equation.push_back("-y7");atom.fpos_equation.push_back("((1.0/2.0)+z7)");// atom B74 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B74 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B74

    atom.name="B"; atom.type=1;                                       // atom B75
    atom.fpos(1)=-x7;atom.fpos(2)=((1.0/2.0)+y7);atom.fpos(3)=((1.0/2.0)-z7);                     // atom B75
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x7");atom.fpos_equation.push_back("((1.0/2.0)+y7)");atom.fpos_equation.push_back("((1.0/2.0)-z7)");// atom B75 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B75 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B75

    atom.name="B"; atom.type=1;                                       // atom B76
    atom.fpos(1)=((1.0/2.0)+x7);atom.fpos(2)=((1.0/2.0)-y7);atom.fpos(3)=-z7;                     // atom B76
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x7)");atom.fpos_equation.push_back("((1.0/2.0)-y7)");atom.fpos_equation.push_back("-z7");// atom B76 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B76 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B76

    atom.name="B"; atom.type=1;                                       // atom B77
    atom.fpos(1)=z7;atom.fpos(2)=x7;atom.fpos(3)=y7;                     // atom B77
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z7");atom.fpos_equation.push_back("x7");atom.fpos_equation.push_back("y7");// atom B77 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B77 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B77

    atom.name="B"; atom.type=1;                                       // atom B78
    atom.fpos(1)=((1.0/2.0)+z7);atom.fpos(2)=((1.0/2.0)-x7);atom.fpos(3)=-y7;                     // atom B78
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z7)");atom.fpos_equation.push_back("((1.0/2.0)-x7)");atom.fpos_equation.push_back("-y7");// atom B78 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B78 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B78

    atom.name="B"; atom.type=1;                                       // atom B79
    atom.fpos(1)=((1.0/2.0)-z7);atom.fpos(2)=-x7;atom.fpos(3)=((1.0/2.0)+y7);                     // atom B79
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z7)");atom.fpos_equation.push_back("-x7");atom.fpos_equation.push_back("((1.0/2.0)+y7)");// atom B79 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B79 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B79

    atom.name="B"; atom.type=1;                                       // atom B80
    atom.fpos(1)=-z7;atom.fpos(2)=((1.0/2.0)+x7);atom.fpos(3)=((1.0/2.0)-y7);                     // atom B80
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z7");atom.fpos_equation.push_back("((1.0/2.0)+x7)");atom.fpos_equation.push_back("((1.0/2.0)-y7)");// atom B80 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B80 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B80

    atom.name="B"; atom.type=1;                                       // atom B81
    atom.fpos(1)=y7;atom.fpos(2)=z7;atom.fpos(3)=x7;                     // atom B81
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y7");atom.fpos_equation.push_back("z7");atom.fpos_equation.push_back("x7");// atom B81 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B81 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B81

    atom.name="B"; atom.type=1;                                       // atom B82
    atom.fpos(1)=-y7;atom.fpos(2)=((1.0/2.0)+z7);atom.fpos(3)=((1.0/2.0)-x7);                     // atom B82
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y7");atom.fpos_equation.push_back("((1.0/2.0)+z7)");atom.fpos_equation.push_back("((1.0/2.0)-x7)");// atom B82 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B82 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B82

    atom.name="B"; atom.type=1;                                       // atom B83
    atom.fpos(1)=((1.0/2.0)+y7);atom.fpos(2)=((1.0/2.0)-z7);atom.fpos(3)=-x7;                     // atom B83
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y7)");atom.fpos_equation.push_back("((1.0/2.0)-z7)");atom.fpos_equation.push_back("-x7");// atom B83 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B83 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B83

    atom.name="B"; atom.type=1;                                       // atom B84
    atom.fpos(1)=((1.0/2.0)-y7);atom.fpos(2)=-z7;atom.fpos(3)=((1.0/2.0)+x7);                     // atom B84
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y7)");atom.fpos_equation.push_back("-z7");atom.fpos_equation.push_back("((1.0/2.0)+x7)");// atom B84 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B84 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B84

    atom.name="B"; atom.type=1;                                       // atom B85
    atom.fpos(1)=-x7;atom.fpos(2)=-y7;atom.fpos(3)=-z7;                     // atom B85
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x7");atom.fpos_equation.push_back("-y7");atom.fpos_equation.push_back("-z7");// atom B85 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B85 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B85

    atom.name="B"; atom.type=1;                                       // atom B86
    atom.fpos(1)=((1.0/2.0)+x7);atom.fpos(2)=y7;atom.fpos(3)=((1.0/2.0)-z7);                     // atom B86
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x7)");atom.fpos_equation.push_back("y7");atom.fpos_equation.push_back("((1.0/2.0)-z7)");// atom B86 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B86 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B86

    atom.name="B"; atom.type=1;                                       // atom B87
    atom.fpos(1)=x7;atom.fpos(2)=((1.0/2.0)-y7);atom.fpos(3)=((1.0/2.0)+z7);                     // atom B87
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x7");atom.fpos_equation.push_back("((1.0/2.0)-y7)");atom.fpos_equation.push_back("((1.0/2.0)+z7)");// atom B87 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B87 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B87

    atom.name="B"; atom.type=1;                                       // atom B88
    atom.fpos(1)=((1.0/2.0)-x7);atom.fpos(2)=((1.0/2.0)+y7);atom.fpos(3)=z7;                     // atom B88
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x7)");atom.fpos_equation.push_back("((1.0/2.0)+y7)");atom.fpos_equation.push_back("z7");// atom B88 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B88 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B88

    atom.name="B"; atom.type=1;                                       // atom B89
    atom.fpos(1)=-z7;atom.fpos(2)=-x7;atom.fpos(3)=-y7;                     // atom B89
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z7");atom.fpos_equation.push_back("-x7");atom.fpos_equation.push_back("-y7");// atom B89 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B89 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B89

    atom.name="B"; atom.type=1;                                       // atom B90
    atom.fpos(1)=((1.0/2.0)-z7);atom.fpos(2)=((1.0/2.0)+x7);atom.fpos(3)=y7;                     // atom B90
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z7)");atom.fpos_equation.push_back("((1.0/2.0)+x7)");atom.fpos_equation.push_back("y7");// atom B90 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B90 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B90

    atom.name="B"; atom.type=1;                                       // atom B91
    atom.fpos(1)=((1.0/2.0)+z7);atom.fpos(2)=x7;atom.fpos(3)=((1.0/2.0)-y7);                     // atom B91
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z7)");atom.fpos_equation.push_back("x7");atom.fpos_equation.push_back("((1.0/2.0)-y7)");// atom B91 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B91 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B91

    atom.name="B"; atom.type=1;                                       // atom B92
    atom.fpos(1)=z7;atom.fpos(2)=((1.0/2.0)-x7);atom.fpos(3)=((1.0/2.0)+y7);                     // atom B92
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z7");atom.fpos_equation.push_back("((1.0/2.0)-x7)");atom.fpos_equation.push_back("((1.0/2.0)+y7)");// atom B92 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B92 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B92

    atom.name="B"; atom.type=1;                                       // atom B93
    atom.fpos(1)=-y7;atom.fpos(2)=-z7;atom.fpos(3)=-x7;                     // atom B93
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y7");atom.fpos_equation.push_back("-z7");atom.fpos_equation.push_back("-x7");// atom B93 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B93 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B93

    atom.name="B"; atom.type=1;                                       // atom B94
    atom.fpos(1)=y7;atom.fpos(2)=((1.0/2.0)-z7);atom.fpos(3)=((1.0/2.0)+x7);                     // atom B94
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y7");atom.fpos_equation.push_back("((1.0/2.0)-z7)");atom.fpos_equation.push_back("((1.0/2.0)+x7)");// atom B94 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B94 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B94

    atom.name="B"; atom.type=1;                                       // atom B95
    atom.fpos(1)=((1.0/2.0)-y7);atom.fpos(2)=((1.0/2.0)+z7);atom.fpos(3)=x7;                     // atom B95
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y7)");atom.fpos_equation.push_back("((1.0/2.0)+z7)");atom.fpos_equation.push_back("x7");// atom B95 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B95 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B95

    atom.name="B"; atom.type=1;                                       // atom B96
    atom.fpos(1)=((1.0/2.0)+y7);atom.fpos(2)=z7;atom.fpos(3)=((1.0/2.0)-x7);                     // atom B96
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y7)");atom.fpos_equation.push_back("z7");atom.fpos_equation.push_back("((1.0/2.0)-x7)");// atom B96 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B96 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B96

    atom.name="B"; atom.type=1;                                       // atom B97
    atom.fpos(1)=x8;atom.fpos(2)=y8;atom.fpos(3)=z8;                     // atom B97
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x8");atom.fpos_equation.push_back("y8");atom.fpos_equation.push_back("z8");// atom B97 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B97 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B97

    atom.name="B"; atom.type=1;                                       // atom B98
    atom.fpos(1)=((1.0/2.0)-x8);atom.fpos(2)=-y8;atom.fpos(3)=((1.0/2.0)+z8);                     // atom B98
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x8)");atom.fpos_equation.push_back("-y8");atom.fpos_equation.push_back("((1.0/2.0)+z8)");// atom B98 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B98 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B98

    atom.name="B"; atom.type=1;                                       // atom B99
    atom.fpos(1)=-x8;atom.fpos(2)=((1.0/2.0)+y8);atom.fpos(3)=((1.0/2.0)-z8);                     // atom B99
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x8");atom.fpos_equation.push_back("((1.0/2.0)+y8)");atom.fpos_equation.push_back("((1.0/2.0)-z8)");// atom B99 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B99 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B99

    atom.name="B"; atom.type=1;                                       // atom B100
    atom.fpos(1)=((1.0/2.0)+x8);atom.fpos(2)=((1.0/2.0)-y8);atom.fpos(3)=-z8;                     // atom B100
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x8)");atom.fpos_equation.push_back("((1.0/2.0)-y8)");atom.fpos_equation.push_back("-z8");// atom B100 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B100 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B100

    atom.name="B"; atom.type=1;                                       // atom B101
    atom.fpos(1)=z8;atom.fpos(2)=x8;atom.fpos(3)=y8;                     // atom B101
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z8");atom.fpos_equation.push_back("x8");atom.fpos_equation.push_back("y8");// atom B101 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B101 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B101

    atom.name="B"; atom.type=1;                                       // atom B102
    atom.fpos(1)=((1.0/2.0)+z8);atom.fpos(2)=((1.0/2.0)-x8);atom.fpos(3)=-y8;                     // atom B102
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z8)");atom.fpos_equation.push_back("((1.0/2.0)-x8)");atom.fpos_equation.push_back("-y8");// atom B102 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B102 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B102

    atom.name="B"; atom.type=1;                                       // atom B103
    atom.fpos(1)=((1.0/2.0)-z8);atom.fpos(2)=-x8;atom.fpos(3)=((1.0/2.0)+y8);                     // atom B103
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z8)");atom.fpos_equation.push_back("-x8");atom.fpos_equation.push_back("((1.0/2.0)+y8)");// atom B103 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B103 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B103

    atom.name="B"; atom.type=1;                                       // atom B104
    atom.fpos(1)=-z8;atom.fpos(2)=((1.0/2.0)+x8);atom.fpos(3)=((1.0/2.0)-y8);                     // atom B104
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z8");atom.fpos_equation.push_back("((1.0/2.0)+x8)");atom.fpos_equation.push_back("((1.0/2.0)-y8)");// atom B104 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B104 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B104

    atom.name="B"; atom.type=1;                                       // atom B105
    atom.fpos(1)=y8;atom.fpos(2)=z8;atom.fpos(3)=x8;                     // atom B105
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y8");atom.fpos_equation.push_back("z8");atom.fpos_equation.push_back("x8");// atom B105 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B105 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B105

    atom.name="B"; atom.type=1;                                       // atom B106
    atom.fpos(1)=-y8;atom.fpos(2)=((1.0/2.0)+z8);atom.fpos(3)=((1.0/2.0)-x8);                     // atom B106
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y8");atom.fpos_equation.push_back("((1.0/2.0)+z8)");atom.fpos_equation.push_back("((1.0/2.0)-x8)");// atom B106 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B106 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B106

    atom.name="B"; atom.type=1;                                       // atom B107
    atom.fpos(1)=((1.0/2.0)+y8);atom.fpos(2)=((1.0/2.0)-z8);atom.fpos(3)=-x8;                     // atom B107
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y8)");atom.fpos_equation.push_back("((1.0/2.0)-z8)");atom.fpos_equation.push_back("-x8");// atom B107 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B107 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B107

    atom.name="B"; atom.type=1;                                       // atom B108
    atom.fpos(1)=((1.0/2.0)-y8);atom.fpos(2)=-z8;atom.fpos(3)=((1.0/2.0)+x8);                     // atom B108
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y8)");atom.fpos_equation.push_back("-z8");atom.fpos_equation.push_back("((1.0/2.0)+x8)");// atom B108 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B108 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B108

    atom.name="B"; atom.type=1;                                       // atom B109
    atom.fpos(1)=-x8;atom.fpos(2)=-y8;atom.fpos(3)=-z8;                     // atom B109
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x8");atom.fpos_equation.push_back("-y8");atom.fpos_equation.push_back("-z8");// atom B109 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B109 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B109

    atom.name="B"; atom.type=1;                                       // atom B110
    atom.fpos(1)=((1.0/2.0)+x8);atom.fpos(2)=y8;atom.fpos(3)=((1.0/2.0)-z8);                     // atom B110
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x8)");atom.fpos_equation.push_back("y8");atom.fpos_equation.push_back("((1.0/2.0)-z8)");// atom B110 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B110 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B110

    atom.name="B"; atom.type=1;                                       // atom B111
    atom.fpos(1)=x8;atom.fpos(2)=((1.0/2.0)-y8);atom.fpos(3)=((1.0/2.0)+z8);                     // atom B111
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x8");atom.fpos_equation.push_back("((1.0/2.0)-y8)");atom.fpos_equation.push_back("((1.0/2.0)+z8)");// atom B111 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B111 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B111

    atom.name="B"; atom.type=1;                                       // atom B112
    atom.fpos(1)=((1.0/2.0)-x8);atom.fpos(2)=((1.0/2.0)+y8);atom.fpos(3)=z8;                     // atom B112
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x8)");atom.fpos_equation.push_back("((1.0/2.0)+y8)");atom.fpos_equation.push_back("z8");// atom B112 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B112 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B112

    atom.name="B"; atom.type=1;                                       // atom B113
    atom.fpos(1)=-z8;atom.fpos(2)=-x8;atom.fpos(3)=-y8;                     // atom B113
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z8");atom.fpos_equation.push_back("-x8");atom.fpos_equation.push_back("-y8");// atom B113 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B113 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B113

    atom.name="B"; atom.type=1;                                       // atom B114
    atom.fpos(1)=((1.0/2.0)-z8);atom.fpos(2)=((1.0/2.0)+x8);atom.fpos(3)=y8;                     // atom B114
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z8)");atom.fpos_equation.push_back("((1.0/2.0)+x8)");atom.fpos_equation.push_back("y8");// atom B114 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B114 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B114

    atom.name="B"; atom.type=1;                                       // atom B115
    atom.fpos(1)=((1.0/2.0)+z8);atom.fpos(2)=x8;atom.fpos(3)=((1.0/2.0)-y8);                     // atom B115
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z8)");atom.fpos_equation.push_back("x8");atom.fpos_equation.push_back("((1.0/2.0)-y8)");// atom B115 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B115 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B115

    atom.name="B"; atom.type=1;                                       // atom B116
    atom.fpos(1)=z8;atom.fpos(2)=((1.0/2.0)-x8);atom.fpos(3)=((1.0/2.0)+y8);                     // atom B116
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z8");atom.fpos_equation.push_back("((1.0/2.0)-x8)");atom.fpos_equation.push_back("((1.0/2.0)+y8)");// atom B116 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B116 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B116

    atom.name="B"; atom.type=1;                                       // atom B117
    atom.fpos(1)=-y8;atom.fpos(2)=-z8;atom.fpos(3)=-x8;                     // atom B117
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y8");atom.fpos_equation.push_back("-z8");atom.fpos_equation.push_back("-x8");// atom B117 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B117 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B117

    atom.name="B"; atom.type=1;                                       // atom B118
    atom.fpos(1)=y8;atom.fpos(2)=((1.0/2.0)-z8);atom.fpos(3)=((1.0/2.0)+x8);                     // atom B118
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y8");atom.fpos_equation.push_back("((1.0/2.0)-z8)");atom.fpos_equation.push_back("((1.0/2.0)+x8)");// atom B118 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B118 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B118

    atom.name="B"; atom.type=1;                                       // atom B119
    atom.fpos(1)=((1.0/2.0)-y8);atom.fpos(2)=((1.0/2.0)+z8);atom.fpos(3)=x8;                     // atom B119
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y8)");atom.fpos_equation.push_back("((1.0/2.0)+z8)");atom.fpos_equation.push_back("x8");// atom B119 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B119 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B119

    atom.name="B"; atom.type=1;                                       // atom B120
    atom.fpos(1)=((1.0/2.0)+y8);atom.fpos(2)=z8;atom.fpos(3)=((1.0/2.0)-x8);                     // atom B120
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y8)");atom.fpos_equation.push_back("z8");atom.fpos_equation.push_back("((1.0/2.0)-x8)");// atom B120 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B120 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B120

    atom.name="C"; atom.type=2;                                       // atom B121
    atom.fpos(1)=x9;atom.fpos(2)=y9;atom.fpos(3)=z9;                     // atom B121
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x9");atom.fpos_equation.push_back("y9");atom.fpos_equation.push_back("z9");// atom B121 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B121 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B121

    atom.name="C"; atom.type=2;                                       // atom B122
    atom.fpos(1)=((1.0/2.0)-x9);atom.fpos(2)=-y9;atom.fpos(3)=((1.0/2.0)+z9);                     // atom B122
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x9)");atom.fpos_equation.push_back("-y9");atom.fpos_equation.push_back("((1.0/2.0)+z9)");// atom B122 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B122 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B122

    atom.name="C"; atom.type=2;                                       // atom B123
    atom.fpos(1)=-x9;atom.fpos(2)=((1.0/2.0)+y9);atom.fpos(3)=((1.0/2.0)-z9);                     // atom B123
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x9");atom.fpos_equation.push_back("((1.0/2.0)+y9)");atom.fpos_equation.push_back("((1.0/2.0)-z9)");// atom B123 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B123 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B123

    atom.name="C"; atom.type=2;                                       // atom B124
    atom.fpos(1)=((1.0/2.0)+x9);atom.fpos(2)=((1.0/2.0)-y9);atom.fpos(3)=-z9;                     // atom B124
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x9)");atom.fpos_equation.push_back("((1.0/2.0)-y9)");atom.fpos_equation.push_back("-z9");// atom B124 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B124 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B124

    atom.name="C"; atom.type=2;                                       // atom B125
    atom.fpos(1)=z9;atom.fpos(2)=x9;atom.fpos(3)=y9;                     // atom B125
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z9");atom.fpos_equation.push_back("x9");atom.fpos_equation.push_back("y9");// atom B125 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B125 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B125

    atom.name="C"; atom.type=2;                                       // atom B126
    atom.fpos(1)=((1.0/2.0)+z9);atom.fpos(2)=((1.0/2.0)-x9);atom.fpos(3)=-y9;                     // atom B126
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z9)");atom.fpos_equation.push_back("((1.0/2.0)-x9)");atom.fpos_equation.push_back("-y9");// atom B126 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B126 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B126

    atom.name="C"; atom.type=2;                                       // atom B127
    atom.fpos(1)=((1.0/2.0)-z9);atom.fpos(2)=-x9;atom.fpos(3)=((1.0/2.0)+y9);                     // atom B127
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z9)");atom.fpos_equation.push_back("-x9");atom.fpos_equation.push_back("((1.0/2.0)+y9)");// atom B127 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B127 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B127

    atom.name="C"; atom.type=2;                                       // atom B128
    atom.fpos(1)=-z9;atom.fpos(2)=((1.0/2.0)+x9);atom.fpos(3)=((1.0/2.0)-y9);                     // atom B128
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z9");atom.fpos_equation.push_back("((1.0/2.0)+x9)");atom.fpos_equation.push_back("((1.0/2.0)-y9)");// atom B128 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B128 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B128

    atom.name="C"; atom.type=2;                                       // atom B129
    atom.fpos(1)=y9;atom.fpos(2)=z9;atom.fpos(3)=x9;                     // atom B129
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y9");atom.fpos_equation.push_back("z9");atom.fpos_equation.push_back("x9");// atom B129 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B129 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B129

    atom.name="C"; atom.type=2;                                       // atom B130
    atom.fpos(1)=-y9;atom.fpos(2)=((1.0/2.0)+z9);atom.fpos(3)=((1.0/2.0)-x9);                     // atom B130
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y9");atom.fpos_equation.push_back("((1.0/2.0)+z9)");atom.fpos_equation.push_back("((1.0/2.0)-x9)");// atom B130 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B130 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B130

    atom.name="C"; atom.type=2;                                       // atom B131
    atom.fpos(1)=((1.0/2.0)+y9);atom.fpos(2)=((1.0/2.0)-z9);atom.fpos(3)=-x9;                     // atom B131
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y9)");atom.fpos_equation.push_back("((1.0/2.0)-z9)");atom.fpos_equation.push_back("-x9");// atom B131 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B131 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B131

    atom.name="C"; atom.type=2;                                       // atom B132
    atom.fpos(1)=((1.0/2.0)-y9);atom.fpos(2)=-z9;atom.fpos(3)=((1.0/2.0)+x9);                     // atom B132
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y9)");atom.fpos_equation.push_back("-z9");atom.fpos_equation.push_back("((1.0/2.0)+x9)");// atom B132 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B132 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B132

    atom.name="C"; atom.type=2;                                       // atom B133
    atom.fpos(1)=-x9;atom.fpos(2)=-y9;atom.fpos(3)=-z9;                     // atom B133
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x9");atom.fpos_equation.push_back("-y9");atom.fpos_equation.push_back("-z9");// atom B133 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B133 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B133

    atom.name="C"; atom.type=2;                                       // atom B134
    atom.fpos(1)=((1.0/2.0)+x9);atom.fpos(2)=y9;atom.fpos(3)=((1.0/2.0)-z9);                     // atom B134
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x9)");atom.fpos_equation.push_back("y9");atom.fpos_equation.push_back("((1.0/2.0)-z9)");// atom B134 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B134 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B134

    atom.name="C"; atom.type=2;                                       // atom B135
    atom.fpos(1)=x9;atom.fpos(2)=((1.0/2.0)-y9);atom.fpos(3)=((1.0/2.0)+z9);                     // atom B135
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x9");atom.fpos_equation.push_back("((1.0/2.0)-y9)");atom.fpos_equation.push_back("((1.0/2.0)+z9)");// atom B135 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B135 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B135

    atom.name="C"; atom.type=2;                                       // atom B136
    atom.fpos(1)=((1.0/2.0)-x9);atom.fpos(2)=((1.0/2.0)+y9);atom.fpos(3)=z9;                     // atom B136
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x9)");atom.fpos_equation.push_back("((1.0/2.0)+y9)");atom.fpos_equation.push_back("z9");// atom B136 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B136 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B136

    atom.name="C"; atom.type=2;                                       // atom B137
    atom.fpos(1)=-z9;atom.fpos(2)=-x9;atom.fpos(3)=-y9;                     // atom B137
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z9");atom.fpos_equation.push_back("-x9");atom.fpos_equation.push_back("-y9");// atom B137 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B137 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B137

    atom.name="C"; atom.type=2;                                       // atom B138
    atom.fpos(1)=((1.0/2.0)-z9);atom.fpos(2)=((1.0/2.0)+x9);atom.fpos(3)=y9;                     // atom B138
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z9)");atom.fpos_equation.push_back("((1.0/2.0)+x9)");atom.fpos_equation.push_back("y9");// atom B138 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B138 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B138

    atom.name="C"; atom.type=2;                                       // atom B139
    atom.fpos(1)=((1.0/2.0)+z9);atom.fpos(2)=x9;atom.fpos(3)=((1.0/2.0)-y9);                     // atom B139
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z9)");atom.fpos_equation.push_back("x9");atom.fpos_equation.push_back("((1.0/2.0)-y9)");// atom B139 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B139 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B139

    atom.name="C"; atom.type=2;                                       // atom B140
    atom.fpos(1)=z9;atom.fpos(2)=((1.0/2.0)-x9);atom.fpos(3)=((1.0/2.0)+y9);                     // atom B140
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z9");atom.fpos_equation.push_back("((1.0/2.0)-x9)");atom.fpos_equation.push_back("((1.0/2.0)+y9)");// atom B140 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B140 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B140

    atom.name="C"; atom.type=2;                                       // atom B141
    atom.fpos(1)=-y9;atom.fpos(2)=-z9;atom.fpos(3)=-x9;                     // atom B141
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y9");atom.fpos_equation.push_back("-z9");atom.fpos_equation.push_back("-x9");// atom B141 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B141 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B141

    atom.name="C"; atom.type=2;                                       // atom B142
    atom.fpos(1)=y9;atom.fpos(2)=((1.0/2.0)-z9);atom.fpos(3)=((1.0/2.0)+x9);                     // atom B142
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y9");atom.fpos_equation.push_back("((1.0/2.0)-z9)");atom.fpos_equation.push_back("((1.0/2.0)+x9)");// atom B142 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B142 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B142

    atom.name="C"; atom.type=2;                                       // atom B143
    atom.fpos(1)=((1.0/2.0)-y9);atom.fpos(2)=((1.0/2.0)+z9);atom.fpos(3)=x9;                     // atom B143
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y9)");atom.fpos_equation.push_back("((1.0/2.0)+z9)");atom.fpos_equation.push_back("x9");// atom B143 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B143 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B143

    atom.name="C"; atom.type=2;                                       // atom B144
    atom.fpos(1)=((1.0/2.0)+y9);atom.fpos(2)=z9;atom.fpos(3)=((1.0/2.0)-x9);                     // atom B144
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y9)");atom.fpos_equation.push_back("z9");atom.fpos_equation.push_back("((1.0/2.0)-x9)");// atom B144 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B144 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B144

    atom.name="C"; atom.type=2;                                       // atom B145
    atom.fpos(1)=x10;atom.fpos(2)=y10;atom.fpos(3)=z10;                     // atom B145
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x10");atom.fpos_equation.push_back("y10");atom.fpos_equation.push_back("z10");// atom B145 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B145 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B145

    atom.name="C"; atom.type=2;                                       // atom B146
    atom.fpos(1)=((1.0/2.0)-x10);atom.fpos(2)=-y10;atom.fpos(3)=((1.0/2.0)+z10);                     // atom B146
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x10)");atom.fpos_equation.push_back("-y10");atom.fpos_equation.push_back("((1.0/2.0)+z10)");// atom B146 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B146 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B146

    atom.name="C"; atom.type=2;                                       // atom B147
    atom.fpos(1)=-x10;atom.fpos(2)=((1.0/2.0)+y10);atom.fpos(3)=((1.0/2.0)-z10);                     // atom B147
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x10");atom.fpos_equation.push_back("((1.0/2.0)+y10)");atom.fpos_equation.push_back("((1.0/2.0)-z10)");// atom B147 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B147 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B147

    atom.name="C"; atom.type=2;                                       // atom B148
    atom.fpos(1)=((1.0/2.0)+x10);atom.fpos(2)=((1.0/2.0)-y10);atom.fpos(3)=-z10;                     // atom B148
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x10)");atom.fpos_equation.push_back("((1.0/2.0)-y10)");atom.fpos_equation.push_back("-z10");// atom B148 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B148 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B148

    atom.name="C"; atom.type=2;                                       // atom B149
    atom.fpos(1)=z10;atom.fpos(2)=x10;atom.fpos(3)=y10;                     // atom B149
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z10");atom.fpos_equation.push_back("x10");atom.fpos_equation.push_back("y10");// atom B149 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B149 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B149

    atom.name="C"; atom.type=2;                                       // atom B150
    atom.fpos(1)=((1.0/2.0)+z10);atom.fpos(2)=((1.0/2.0)-x10);atom.fpos(3)=-y10;                     // atom B150
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z10)");atom.fpos_equation.push_back("((1.0/2.0)-x10)");atom.fpos_equation.push_back("-y10");// atom B150 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B150 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B150

    atom.name="C"; atom.type=2;                                       // atom B151
    atom.fpos(1)=((1.0/2.0)-z10);atom.fpos(2)=-x10;atom.fpos(3)=((1.0/2.0)+y10);                     // atom B151
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z10)");atom.fpos_equation.push_back("-x10");atom.fpos_equation.push_back("((1.0/2.0)+y10)");// atom B151 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B151 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B151

    atom.name="C"; atom.type=2;                                       // atom B152
    atom.fpos(1)=-z10;atom.fpos(2)=((1.0/2.0)+x10);atom.fpos(3)=((1.0/2.0)-y10);                     // atom B152
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z10");atom.fpos_equation.push_back("((1.0/2.0)+x10)");atom.fpos_equation.push_back("((1.0/2.0)-y10)");// atom B152 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B152 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B152

    atom.name="C"; atom.type=2;                                       // atom B153
    atom.fpos(1)=y10;atom.fpos(2)=z10;atom.fpos(3)=x10;                     // atom B153
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y10");atom.fpos_equation.push_back("z10");atom.fpos_equation.push_back("x10");// atom B153 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B153 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B153

    atom.name="C"; atom.type=2;                                       // atom B154
    atom.fpos(1)=-y10;atom.fpos(2)=((1.0/2.0)+z10);atom.fpos(3)=((1.0/2.0)-x10);                     // atom B154
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y10");atom.fpos_equation.push_back("((1.0/2.0)+z10)");atom.fpos_equation.push_back("((1.0/2.0)-x10)");// atom B154 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B154 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B154

    atom.name="C"; atom.type=2;                                       // atom B155
    atom.fpos(1)=((1.0/2.0)+y10);atom.fpos(2)=((1.0/2.0)-z10);atom.fpos(3)=-x10;                     // atom B155
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y10)");atom.fpos_equation.push_back("((1.0/2.0)-z10)");atom.fpos_equation.push_back("-x10");// atom B155 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B155 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B155

    atom.name="C"; atom.type=2;                                       // atom B156
    atom.fpos(1)=((1.0/2.0)-y10);atom.fpos(2)=-z10;atom.fpos(3)=((1.0/2.0)+x10);                     // atom B156
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y10)");atom.fpos_equation.push_back("-z10");atom.fpos_equation.push_back("((1.0/2.0)+x10)");// atom B156 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B156 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B156

    atom.name="C"; atom.type=2;                                       // atom B157
    atom.fpos(1)=-x10;atom.fpos(2)=-y10;atom.fpos(3)=-z10;                     // atom B157
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x10");atom.fpos_equation.push_back("-y10");atom.fpos_equation.push_back("-z10");// atom B157 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B157 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B157

    atom.name="C"; atom.type=2;                                       // atom B158
    atom.fpos(1)=((1.0/2.0)+x10);atom.fpos(2)=y10;atom.fpos(3)=((1.0/2.0)-z10);                     // atom B158
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x10)");atom.fpos_equation.push_back("y10");atom.fpos_equation.push_back("((1.0/2.0)-z10)");// atom B158 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B158 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B158

    atom.name="C"; atom.type=2;                                       // atom B159
    atom.fpos(1)=x10;atom.fpos(2)=((1.0/2.0)-y10);atom.fpos(3)=((1.0/2.0)+z10);                     // atom B159
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x10");atom.fpos_equation.push_back("((1.0/2.0)-y10)");atom.fpos_equation.push_back("((1.0/2.0)+z10)");// atom B159 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B159 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B159

    atom.name="C"; atom.type=2;                                       // atom B160
    atom.fpos(1)=((1.0/2.0)-x10);atom.fpos(2)=((1.0/2.0)+y10);atom.fpos(3)=z10;                     // atom B160
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x10)");atom.fpos_equation.push_back("((1.0/2.0)+y10)");atom.fpos_equation.push_back("z10");// atom B160 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B160 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B160

    atom.name="C"; atom.type=2;                                       // atom B161
    atom.fpos(1)=-z10;atom.fpos(2)=-x10;atom.fpos(3)=-y10;                     // atom B161
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z10");atom.fpos_equation.push_back("-x10");atom.fpos_equation.push_back("-y10");// atom B161 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B161 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B161

    atom.name="C"; atom.type=2;                                       // atom B162
    atom.fpos(1)=((1.0/2.0)-z10);atom.fpos(2)=((1.0/2.0)+x10);atom.fpos(3)=y10;                     // atom B162
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z10)");atom.fpos_equation.push_back("((1.0/2.0)+x10)");atom.fpos_equation.push_back("y10");// atom B162 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B162 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B162

    atom.name="C"; atom.type=2;                                       // atom B163
    atom.fpos(1)=((1.0/2.0)+z10);atom.fpos(2)=x10;atom.fpos(3)=((1.0/2.0)-y10);                     // atom B163
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z10)");atom.fpos_equation.push_back("x10");atom.fpos_equation.push_back("((1.0/2.0)-y10)");// atom B163 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B163 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B163

    atom.name="C"; atom.type=2;                                       // atom B164
    atom.fpos(1)=z10;atom.fpos(2)=((1.0/2.0)-x10);atom.fpos(3)=((1.0/2.0)+y10);                     // atom B164
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z10");atom.fpos_equation.push_back("((1.0/2.0)-x10)");atom.fpos_equation.push_back("((1.0/2.0)+y10)");// atom B164 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B164 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B164

    atom.name="C"; atom.type=2;                                       // atom B165
    atom.fpos(1)=-y10;atom.fpos(2)=-z10;atom.fpos(3)=-x10;                     // atom B165
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y10");atom.fpos_equation.push_back("-z10");atom.fpos_equation.push_back("-x10");// atom B165 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B165 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B165

    atom.name="C"; atom.type=2;                                       // atom B166
    atom.fpos(1)=y10;atom.fpos(2)=((1.0/2.0)-z10);atom.fpos(3)=((1.0/2.0)+x10);                     // atom B166
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y10");atom.fpos_equation.push_back("((1.0/2.0)-z10)");atom.fpos_equation.push_back("((1.0/2.0)+x10)");// atom B166 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B166 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B166

    atom.name="C"; atom.type=2;                                       // atom B167
    atom.fpos(1)=((1.0/2.0)-y10);atom.fpos(2)=((1.0/2.0)+z10);atom.fpos(3)=x10;                     // atom B167
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y10)");atom.fpos_equation.push_back("((1.0/2.0)+z10)");atom.fpos_equation.push_back("x10");// atom B167 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B167 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B167

    atom.name="C"; atom.type=2;                                       // atom B168
    atom.fpos(1)=((1.0/2.0)+y10);atom.fpos(2)=z10;atom.fpos(3)=((1.0/2.0)-x10);                     // atom B168
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y10)");atom.fpos_equation.push_back("z10");atom.fpos_equation.push_back("((1.0/2.0)-x10)");// atom B168 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B168 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B168

    atom.name="C"; atom.type=2;                                       // atom B169
    atom.fpos(1)=x11;atom.fpos(2)=y11;atom.fpos(3)=z11;                     // atom B169
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x11");atom.fpos_equation.push_back("y11");atom.fpos_equation.push_back("z11");// atom B169 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B169 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B169

    atom.name="C"; atom.type=2;                                       // atom B170
    atom.fpos(1)=((1.0/2.0)-x11);atom.fpos(2)=-y11;atom.fpos(3)=((1.0/2.0)+z11);                     // atom B170
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x11)");atom.fpos_equation.push_back("-y11");atom.fpos_equation.push_back("((1.0/2.0)+z11)");// atom B170 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B170 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B170

    atom.name="C"; atom.type=2;                                       // atom B171
    atom.fpos(1)=-x11;atom.fpos(2)=((1.0/2.0)+y11);atom.fpos(3)=((1.0/2.0)-z11);                     // atom B171
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x11");atom.fpos_equation.push_back("((1.0/2.0)+y11)");atom.fpos_equation.push_back("((1.0/2.0)-z11)");// atom B171 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B171 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B171

    atom.name="C"; atom.type=2;                                       // atom B172
    atom.fpos(1)=((1.0/2.0)+x11);atom.fpos(2)=((1.0/2.0)-y11);atom.fpos(3)=-z11;                     // atom B172
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x11)");atom.fpos_equation.push_back("((1.0/2.0)-y11)");atom.fpos_equation.push_back("-z11");// atom B172 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B172 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B172

    atom.name="C"; atom.type=2;                                       // atom B173
    atom.fpos(1)=z11;atom.fpos(2)=x11;atom.fpos(3)=y11;                     // atom B173
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z11");atom.fpos_equation.push_back("x11");atom.fpos_equation.push_back("y11");// atom B173 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B173 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B173

    atom.name="C"; atom.type=2;                                       // atom B174
    atom.fpos(1)=((1.0/2.0)+z11);atom.fpos(2)=((1.0/2.0)-x11);atom.fpos(3)=-y11;                     // atom B174
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z11)");atom.fpos_equation.push_back("((1.0/2.0)-x11)");atom.fpos_equation.push_back("-y11");// atom B174 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B174 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B174

    atom.name="C"; atom.type=2;                                       // atom B175
    atom.fpos(1)=((1.0/2.0)-z11);atom.fpos(2)=-x11;atom.fpos(3)=((1.0/2.0)+y11);                     // atom B175
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z11)");atom.fpos_equation.push_back("-x11");atom.fpos_equation.push_back("((1.0/2.0)+y11)");// atom B175 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B175 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B175

    atom.name="C"; atom.type=2;                                       // atom B176
    atom.fpos(1)=-z11;atom.fpos(2)=((1.0/2.0)+x11);atom.fpos(3)=((1.0/2.0)-y11);                     // atom B176
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z11");atom.fpos_equation.push_back("((1.0/2.0)+x11)");atom.fpos_equation.push_back("((1.0/2.0)-y11)");// atom B176 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B176 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B176

    atom.name="C"; atom.type=2;                                       // atom B177
    atom.fpos(1)=y11;atom.fpos(2)=z11;atom.fpos(3)=x11;                     // atom B177
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y11");atom.fpos_equation.push_back("z11");atom.fpos_equation.push_back("x11");// atom B177 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B177 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B177

    atom.name="C"; atom.type=2;                                       // atom B178
    atom.fpos(1)=-y11;atom.fpos(2)=((1.0/2.0)+z11);atom.fpos(3)=((1.0/2.0)-x11);                     // atom B178
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y11");atom.fpos_equation.push_back("((1.0/2.0)+z11)");atom.fpos_equation.push_back("((1.0/2.0)-x11)");// atom B178 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B178 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B178

    atom.name="C"; atom.type=2;                                       // atom B179
    atom.fpos(1)=((1.0/2.0)+y11);atom.fpos(2)=((1.0/2.0)-z11);atom.fpos(3)=-x11;                     // atom B179
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y11)");atom.fpos_equation.push_back("((1.0/2.0)-z11)");atom.fpos_equation.push_back("-x11");// atom B179 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B179 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B179

    atom.name="C"; atom.type=2;                                       // atom B180
    atom.fpos(1)=((1.0/2.0)-y11);atom.fpos(2)=-z11;atom.fpos(3)=((1.0/2.0)+x11);                     // atom B180
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y11)");atom.fpos_equation.push_back("-z11");atom.fpos_equation.push_back("((1.0/2.0)+x11)");// atom B180 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B180 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B180

    atom.name="C"; atom.type=2;                                       // atom B181
    atom.fpos(1)=-x11;atom.fpos(2)=-y11;atom.fpos(3)=-z11;                     // atom B181
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x11");atom.fpos_equation.push_back("-y11");atom.fpos_equation.push_back("-z11");// atom B181 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B181 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B181

    atom.name="C"; atom.type=2;                                       // atom B182
    atom.fpos(1)=((1.0/2.0)+x11);atom.fpos(2)=y11;atom.fpos(3)=((1.0/2.0)-z11);                     // atom B182
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x11)");atom.fpos_equation.push_back("y11");atom.fpos_equation.push_back("((1.0/2.0)-z11)");// atom B182 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B182 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B182

    atom.name="C"; atom.type=2;                                       // atom B183
    atom.fpos(1)=x11;atom.fpos(2)=((1.0/2.0)-y11);atom.fpos(3)=((1.0/2.0)+z11);                     // atom B183
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x11");atom.fpos_equation.push_back("((1.0/2.0)-y11)");atom.fpos_equation.push_back("((1.0/2.0)+z11)");// atom B183 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B183 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B183

    atom.name="C"; atom.type=2;                                       // atom B184
    atom.fpos(1)=((1.0/2.0)-x11);atom.fpos(2)=((1.0/2.0)+y11);atom.fpos(3)=z11;                     // atom B184
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x11)");atom.fpos_equation.push_back("((1.0/2.0)+y11)");atom.fpos_equation.push_back("z11");// atom B184 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B184 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B184

    atom.name="C"; atom.type=2;                                       // atom B185
    atom.fpos(1)=-z11;atom.fpos(2)=-x11;atom.fpos(3)=-y11;                     // atom B185
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z11");atom.fpos_equation.push_back("-x11");atom.fpos_equation.push_back("-y11");// atom B185 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B185 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B185

    atom.name="C"; atom.type=2;                                       // atom B186
    atom.fpos(1)=((1.0/2.0)-z11);atom.fpos(2)=((1.0/2.0)+x11);atom.fpos(3)=y11;                     // atom B186
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z11)");atom.fpos_equation.push_back("((1.0/2.0)+x11)");atom.fpos_equation.push_back("y11");// atom B186 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B186 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B186

    atom.name="C"; atom.type=2;                                       // atom B187
    atom.fpos(1)=((1.0/2.0)+z11);atom.fpos(2)=x11;atom.fpos(3)=((1.0/2.0)-y11);                     // atom B187
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z11)");atom.fpos_equation.push_back("x11");atom.fpos_equation.push_back("((1.0/2.0)-y11)");// atom B187 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B187 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B187

    atom.name="C"; atom.type=2;                                       // atom B188
    atom.fpos(1)=z11;atom.fpos(2)=((1.0/2.0)-x11);atom.fpos(3)=((1.0/2.0)+y11);                     // atom B188
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z11");atom.fpos_equation.push_back("((1.0/2.0)-x11)");atom.fpos_equation.push_back("((1.0/2.0)+y11)");// atom B188 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B188 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B188

    atom.name="C"; atom.type=2;                                       // atom B189
    atom.fpos(1)=-y11;atom.fpos(2)=-z11;atom.fpos(3)=-x11;                     // atom B189
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y11");atom.fpos_equation.push_back("-z11");atom.fpos_equation.push_back("-x11");// atom B189 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B189 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B189

    atom.name="C"; atom.type=2;                                       // atom B190
    atom.fpos(1)=y11;atom.fpos(2)=((1.0/2.0)-z11);atom.fpos(3)=((1.0/2.0)+x11);                     // atom B190
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y11");atom.fpos_equation.push_back("((1.0/2.0)-z11)");atom.fpos_equation.push_back("((1.0/2.0)+x11)");// atom B190 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B190 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B190

    atom.name="C"; atom.type=2;                                       // atom B191
    atom.fpos(1)=((1.0/2.0)-y11);atom.fpos(2)=((1.0/2.0)+z11);atom.fpos(3)=x11;                     // atom B191
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y11)");atom.fpos_equation.push_back("((1.0/2.0)+z11)");atom.fpos_equation.push_back("x11");// atom B191 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B191 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B191

    atom.name="C"; atom.type=2;                                       // atom B192
    atom.fpos(1)=((1.0/2.0)+y11);atom.fpos(2)=z11;atom.fpos(3)=((1.0/2.0)-x11);                     // atom B192
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y11)");atom.fpos_equation.push_back("z11");atom.fpos_equation.push_back("((1.0/2.0)-x11)");// atom B192 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B192 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B192

    atom.name="C"; atom.type=2;                                       // atom B193
    atom.fpos(1)=x12;atom.fpos(2)=y12;atom.fpos(3)=z12;                     // atom B193
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x12");atom.fpos_equation.push_back("y12");atom.fpos_equation.push_back("z12");// atom B193 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B193 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B193

    atom.name="C"; atom.type=2;                                       // atom B194
    atom.fpos(1)=((1.0/2.0)-x12);atom.fpos(2)=-y12;atom.fpos(3)=((1.0/2.0)+z12);                     // atom B194
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x12)");atom.fpos_equation.push_back("-y12");atom.fpos_equation.push_back("((1.0/2.0)+z12)");// atom B194 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B194 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B194

    atom.name="C"; atom.type=2;                                       // atom B195
    atom.fpos(1)=-x12;atom.fpos(2)=((1.0/2.0)+y12);atom.fpos(3)=((1.0/2.0)-z12);                     // atom B195
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x12");atom.fpos_equation.push_back("((1.0/2.0)+y12)");atom.fpos_equation.push_back("((1.0/2.0)-z12)");// atom B195 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B195 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B195

    atom.name="C"; atom.type=2;                                       // atom B196
    atom.fpos(1)=((1.0/2.0)+x12);atom.fpos(2)=((1.0/2.0)-y12);atom.fpos(3)=-z12;                     // atom B196
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x12)");atom.fpos_equation.push_back("((1.0/2.0)-y12)");atom.fpos_equation.push_back("-z12");// atom B196 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B196 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B196

    atom.name="C"; atom.type=2;                                       // atom B197
    atom.fpos(1)=z12;atom.fpos(2)=x12;atom.fpos(3)=y12;                     // atom B197
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z12");atom.fpos_equation.push_back("x12");atom.fpos_equation.push_back("y12");// atom B197 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B197 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B197

    atom.name="C"; atom.type=2;                                       // atom B198
    atom.fpos(1)=((1.0/2.0)+z12);atom.fpos(2)=((1.0/2.0)-x12);atom.fpos(3)=-y12;                     // atom B198
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z12)");atom.fpos_equation.push_back("((1.0/2.0)-x12)");atom.fpos_equation.push_back("-y12");// atom B198 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B198 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B198

    atom.name="C"; atom.type=2;                                       // atom B199
    atom.fpos(1)=((1.0/2.0)-z12);atom.fpos(2)=-x12;atom.fpos(3)=((1.0/2.0)+y12);                     // atom B199
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z12)");atom.fpos_equation.push_back("-x12");atom.fpos_equation.push_back("((1.0/2.0)+y12)");// atom B199 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B199 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B199

    atom.name="C"; atom.type=2;                                       // atom B200
    atom.fpos(1)=-z12;atom.fpos(2)=((1.0/2.0)+x12);atom.fpos(3)=((1.0/2.0)-y12);                     // atom B200
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z12");atom.fpos_equation.push_back("((1.0/2.0)+x12)");atom.fpos_equation.push_back("((1.0/2.0)-y12)");// atom B200 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B200 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B200

    atom.name="C"; atom.type=2;                                       // atom B201
    atom.fpos(1)=y12;atom.fpos(2)=z12;atom.fpos(3)=x12;                     // atom B201
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y12");atom.fpos_equation.push_back("z12");atom.fpos_equation.push_back("x12");// atom B201 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B201 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B201

    atom.name="C"; atom.type=2;                                       // atom B202
    atom.fpos(1)=-y12;atom.fpos(2)=((1.0/2.0)+z12);atom.fpos(3)=((1.0/2.0)-x12);                     // atom B202
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y12");atom.fpos_equation.push_back("((1.0/2.0)+z12)");atom.fpos_equation.push_back("((1.0/2.0)-x12)");// atom B202 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B202 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B202

    atom.name="C"; atom.type=2;                                       // atom B203
    atom.fpos(1)=((1.0/2.0)+y12);atom.fpos(2)=((1.0/2.0)-z12);atom.fpos(3)=-x12;                     // atom B203
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y12)");atom.fpos_equation.push_back("((1.0/2.0)-z12)");atom.fpos_equation.push_back("-x12");// atom B203 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B203 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B203

    atom.name="C"; atom.type=2;                                       // atom B204
    atom.fpos(1)=((1.0/2.0)-y12);atom.fpos(2)=-z12;atom.fpos(3)=((1.0/2.0)+x12);                     // atom B204
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y12)");atom.fpos_equation.push_back("-z12");atom.fpos_equation.push_back("((1.0/2.0)+x12)");// atom B204 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B204 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B204

    atom.name="C"; atom.type=2;                                       // atom B205
    atom.fpos(1)=-x12;atom.fpos(2)=-y12;atom.fpos(3)=-z12;                     // atom B205
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x12");atom.fpos_equation.push_back("-y12");atom.fpos_equation.push_back("-z12");// atom B205 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B205 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B205

    atom.name="C"; atom.type=2;                                       // atom B206
    atom.fpos(1)=((1.0/2.0)+x12);atom.fpos(2)=y12;atom.fpos(3)=((1.0/2.0)-z12);                     // atom B206
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x12)");atom.fpos_equation.push_back("y12");atom.fpos_equation.push_back("((1.0/2.0)-z12)");// atom B206 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B206 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B206

    atom.name="C"; atom.type=2;                                       // atom B207
    atom.fpos(1)=x12;atom.fpos(2)=((1.0/2.0)-y12);atom.fpos(3)=((1.0/2.0)+z12);                     // atom B207
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x12");atom.fpos_equation.push_back("((1.0/2.0)-y12)");atom.fpos_equation.push_back("((1.0/2.0)+z12)");// atom B207 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B207 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B207

    atom.name="C"; atom.type=2;                                       // atom B208
    atom.fpos(1)=((1.0/2.0)-x12);atom.fpos(2)=((1.0/2.0)+y12);atom.fpos(3)=z12;                     // atom B208
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x12)");atom.fpos_equation.push_back("((1.0/2.0)+y12)");atom.fpos_equation.push_back("z12");// atom B208 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B208 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B208

    atom.name="C"; atom.type=2;                                       // atom B209
    atom.fpos(1)=-z12;atom.fpos(2)=-x12;atom.fpos(3)=-y12;                     // atom B209
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z12");atom.fpos_equation.push_back("-x12");atom.fpos_equation.push_back("-y12");// atom B209 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B209 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B209

    atom.name="C"; atom.type=2;                                       // atom B210
    atom.fpos(1)=((1.0/2.0)-z12);atom.fpos(2)=((1.0/2.0)+x12);atom.fpos(3)=y12;                     // atom B210
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z12)");atom.fpos_equation.push_back("((1.0/2.0)+x12)");atom.fpos_equation.push_back("y12");// atom B210 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B210 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B210

    atom.name="C"; atom.type=2;                                       // atom B211
    atom.fpos(1)=((1.0/2.0)+z12);atom.fpos(2)=x12;atom.fpos(3)=((1.0/2.0)-y12);                     // atom B211
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z12)");atom.fpos_equation.push_back("x12");atom.fpos_equation.push_back("((1.0/2.0)-y12)");// atom B211 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B211 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B211

    atom.name="C"; atom.type=2;                                       // atom B212
    atom.fpos(1)=z12;atom.fpos(2)=((1.0/2.0)-x12);atom.fpos(3)=((1.0/2.0)+y12);                     // atom B212
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z12");atom.fpos_equation.push_back("((1.0/2.0)-x12)");atom.fpos_equation.push_back("((1.0/2.0)+y12)");// atom B212 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B212 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B212

    atom.name="C"; atom.type=2;                                       // atom B213
    atom.fpos(1)=-y12;atom.fpos(2)=-z12;atom.fpos(3)=-x12;                     // atom B213
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y12");atom.fpos_equation.push_back("-z12");atom.fpos_equation.push_back("-x12");// atom B213 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B213 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B213

    atom.name="C"; atom.type=2;                                       // atom B214
    atom.fpos(1)=y12;atom.fpos(2)=((1.0/2.0)-z12);atom.fpos(3)=((1.0/2.0)+x12);                     // atom B214
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y12");atom.fpos_equation.push_back("((1.0/2.0)-z12)");atom.fpos_equation.push_back("((1.0/2.0)+x12)");// atom B214 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B214 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B214

    atom.name="C"; atom.type=2;                                       // atom B215
    atom.fpos(1)=((1.0/2.0)-y12);atom.fpos(2)=((1.0/2.0)+z12);atom.fpos(3)=x12;                     // atom B215
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y12)");atom.fpos_equation.push_back("((1.0/2.0)+z12)");atom.fpos_equation.push_back("x12");// atom B215 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B215 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B215

    atom.name="C"; atom.type=2;                                       // atom B216
    atom.fpos(1)=((1.0/2.0)+y12);atom.fpos(2)=z12;atom.fpos(3)=((1.0/2.0)-x12);                     // atom B216
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y12)");atom.fpos_equation.push_back("z12");atom.fpos_equation.push_back("((1.0/2.0)-x12)");// atom B216 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B216 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B216

    atom.name="C"; atom.type=2;                                       // atom B217
    atom.fpos(1)=x13;atom.fpos(2)=y13;atom.fpos(3)=z13;                     // atom B217
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x13");atom.fpos_equation.push_back("y13");atom.fpos_equation.push_back("z13");// atom B217 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B217 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B217

    atom.name="C"; atom.type=2;                                       // atom B218
    atom.fpos(1)=((1.0/2.0)-x13);atom.fpos(2)=-y13;atom.fpos(3)=((1.0/2.0)+z13);                     // atom B218
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x13)");atom.fpos_equation.push_back("-y13");atom.fpos_equation.push_back("((1.0/2.0)+z13)");// atom B218 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B218 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B218

    atom.name="C"; atom.type=2;                                       // atom B219
    atom.fpos(1)=-x13;atom.fpos(2)=((1.0/2.0)+y13);atom.fpos(3)=((1.0/2.0)-z13);                     // atom B219
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x13");atom.fpos_equation.push_back("((1.0/2.0)+y13)");atom.fpos_equation.push_back("((1.0/2.0)-z13)");// atom B219 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B219 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B219

    atom.name="C"; atom.type=2;                                       // atom B220
    atom.fpos(1)=((1.0/2.0)+x13);atom.fpos(2)=((1.0/2.0)-y13);atom.fpos(3)=-z13;                     // atom B220
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x13)");atom.fpos_equation.push_back("((1.0/2.0)-y13)");atom.fpos_equation.push_back("-z13");// atom B220 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B220 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B220

    atom.name="C"; atom.type=2;                                       // atom B221
    atom.fpos(1)=z13;atom.fpos(2)=x13;atom.fpos(3)=y13;                     // atom B221
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z13");atom.fpos_equation.push_back("x13");atom.fpos_equation.push_back("y13");// atom B221 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B221 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B221

    atom.name="C"; atom.type=2;                                       // atom B222
    atom.fpos(1)=((1.0/2.0)+z13);atom.fpos(2)=((1.0/2.0)-x13);atom.fpos(3)=-y13;                     // atom B222
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z13)");atom.fpos_equation.push_back("((1.0/2.0)-x13)");atom.fpos_equation.push_back("-y13");// atom B222 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B222 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B222

    atom.name="C"; atom.type=2;                                       // atom B223
    atom.fpos(1)=((1.0/2.0)-z13);atom.fpos(2)=-x13;atom.fpos(3)=((1.0/2.0)+y13);                     // atom B223
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z13)");atom.fpos_equation.push_back("-x13");atom.fpos_equation.push_back("((1.0/2.0)+y13)");// atom B223 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B223 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B223

    atom.name="C"; atom.type=2;                                       // atom B224
    atom.fpos(1)=-z13;atom.fpos(2)=((1.0/2.0)+x13);atom.fpos(3)=((1.0/2.0)-y13);                     // atom B224
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z13");atom.fpos_equation.push_back("((1.0/2.0)+x13)");atom.fpos_equation.push_back("((1.0/2.0)-y13)");// atom B224 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B224 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B224

    atom.name="C"; atom.type=2;                                       // atom B225
    atom.fpos(1)=y13;atom.fpos(2)=z13;atom.fpos(3)=x13;                     // atom B225
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y13");atom.fpos_equation.push_back("z13");atom.fpos_equation.push_back("x13");// atom B225 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B225 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B225

    atom.name="C"; atom.type=2;                                       // atom B226
    atom.fpos(1)=-y13;atom.fpos(2)=((1.0/2.0)+z13);atom.fpos(3)=((1.0/2.0)-x13);                     // atom B226
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y13");atom.fpos_equation.push_back("((1.0/2.0)+z13)");atom.fpos_equation.push_back("((1.0/2.0)-x13)");// atom B226 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B226 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B226

    atom.name="C"; atom.type=2;                                       // atom B227
    atom.fpos(1)=((1.0/2.0)+y13);atom.fpos(2)=((1.0/2.0)-z13);atom.fpos(3)=-x13;                     // atom B227
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y13)");atom.fpos_equation.push_back("((1.0/2.0)-z13)");atom.fpos_equation.push_back("-x13");// atom B227 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B227 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B227

    atom.name="C"; atom.type=2;                                       // atom B228
    atom.fpos(1)=((1.0/2.0)-y13);atom.fpos(2)=-z13;atom.fpos(3)=((1.0/2.0)+x13);                     // atom B228
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y13)");atom.fpos_equation.push_back("-z13");atom.fpos_equation.push_back("((1.0/2.0)+x13)");// atom B228 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B228 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B228

    atom.name="C"; atom.type=2;                                       // atom B229
    atom.fpos(1)=-x13;atom.fpos(2)=-y13;atom.fpos(3)=-z13;                     // atom B229
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x13");atom.fpos_equation.push_back("-y13");atom.fpos_equation.push_back("-z13");// atom B229 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B229 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B229

    atom.name="C"; atom.type=2;                                       // atom B230
    atom.fpos(1)=((1.0/2.0)+x13);atom.fpos(2)=y13;atom.fpos(3)=((1.0/2.0)-z13);                     // atom B230
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x13)");atom.fpos_equation.push_back("y13");atom.fpos_equation.push_back("((1.0/2.0)-z13)");// atom B230 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B230 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B230

    atom.name="C"; atom.type=2;                                       // atom B231
    atom.fpos(1)=x13;atom.fpos(2)=((1.0/2.0)-y13);atom.fpos(3)=((1.0/2.0)+z13);                     // atom B231
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x13");atom.fpos_equation.push_back("((1.0/2.0)-y13)");atom.fpos_equation.push_back("((1.0/2.0)+z13)");// atom B231 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B231 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B231

    atom.name="C"; atom.type=2;                                       // atom B232
    atom.fpos(1)=((1.0/2.0)-x13);atom.fpos(2)=((1.0/2.0)+y13);atom.fpos(3)=z13;                     // atom B232
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x13)");atom.fpos_equation.push_back("((1.0/2.0)+y13)");atom.fpos_equation.push_back("z13");// atom B232 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B232 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B232

    atom.name="C"; atom.type=2;                                       // atom B233
    atom.fpos(1)=-z13;atom.fpos(2)=-x13;atom.fpos(3)=-y13;                     // atom B233
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z13");atom.fpos_equation.push_back("-x13");atom.fpos_equation.push_back("-y13");// atom B233 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B233 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B233

    atom.name="C"; atom.type=2;                                       // atom B234
    atom.fpos(1)=((1.0/2.0)-z13);atom.fpos(2)=((1.0/2.0)+x13);atom.fpos(3)=y13;                     // atom B234
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z13)");atom.fpos_equation.push_back("((1.0/2.0)+x13)");atom.fpos_equation.push_back("y13");// atom B234 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B234 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B234

    atom.name="C"; atom.type=2;                                       // atom B235
    atom.fpos(1)=((1.0/2.0)+z13);atom.fpos(2)=x13;atom.fpos(3)=((1.0/2.0)-y13);                     // atom B235
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z13)");atom.fpos_equation.push_back("x13");atom.fpos_equation.push_back("((1.0/2.0)-y13)");// atom B235 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B235 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B235

    atom.name="C"; atom.type=2;                                       // atom B236
    atom.fpos(1)=z13;atom.fpos(2)=((1.0/2.0)-x13);atom.fpos(3)=((1.0/2.0)+y13);                     // atom B236
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z13");atom.fpos_equation.push_back("((1.0/2.0)-x13)");atom.fpos_equation.push_back("((1.0/2.0)+y13)");// atom B236 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B236 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B236

    atom.name="C"; atom.type=2;                                       // atom B237
    atom.fpos(1)=-y13;atom.fpos(2)=-z13;atom.fpos(3)=-x13;                     // atom B237
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y13");atom.fpos_equation.push_back("-z13");atom.fpos_equation.push_back("-x13");// atom B237 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B237 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B237

    atom.name="C"; atom.type=2;                                       // atom B238
    atom.fpos(1)=y13;atom.fpos(2)=((1.0/2.0)-z13);atom.fpos(3)=((1.0/2.0)+x13);                     // atom B238
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y13");atom.fpos_equation.push_back("((1.0/2.0)-z13)");atom.fpos_equation.push_back("((1.0/2.0)+x13)");// atom B238 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B238 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B238

    atom.name="C"; atom.type=2;                                       // atom B239
    atom.fpos(1)=((1.0/2.0)-y13);atom.fpos(2)=((1.0/2.0)+z13);atom.fpos(3)=x13;                     // atom B239
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y13)");atom.fpos_equation.push_back("((1.0/2.0)+z13)");atom.fpos_equation.push_back("x13");// atom B239 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B239 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B239

    atom.name="C"; atom.type=2;                                       // atom B240
    atom.fpos(1)=((1.0/2.0)+y13);atom.fpos(2)=z13;atom.fpos(3)=((1.0/2.0)-x13);                     // atom B240
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y13)");atom.fpos_equation.push_back("z13");atom.fpos_equation.push_back("((1.0/2.0)-x13)");// atom B240 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B240 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B240

    atom.name="C"; atom.type=2;                                       // atom B241
    atom.fpos(1)=x14;atom.fpos(2)=y14;atom.fpos(3)=z14;                     // atom B241
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x14");atom.fpos_equation.push_back("y14");atom.fpos_equation.push_back("z14");// atom B241 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B241 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B241

    atom.name="C"; atom.type=2;                                       // atom B242
    atom.fpos(1)=((1.0/2.0)-x14);atom.fpos(2)=-y14;atom.fpos(3)=((1.0/2.0)+z14);                     // atom B242
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x14)");atom.fpos_equation.push_back("-y14");atom.fpos_equation.push_back("((1.0/2.0)+z14)");// atom B242 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B242 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B242

    atom.name="C"; atom.type=2;                                       // atom B243
    atom.fpos(1)=-x14;atom.fpos(2)=((1.0/2.0)+y14);atom.fpos(3)=((1.0/2.0)-z14);                     // atom B243
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x14");atom.fpos_equation.push_back("((1.0/2.0)+y14)");atom.fpos_equation.push_back("((1.0/2.0)-z14)");// atom B243 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B243 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B243

    atom.name="C"; atom.type=2;                                       // atom B244
    atom.fpos(1)=((1.0/2.0)+x14);atom.fpos(2)=((1.0/2.0)-y14);atom.fpos(3)=-z14;                     // atom B244
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x14)");atom.fpos_equation.push_back("((1.0/2.0)-y14)");atom.fpos_equation.push_back("-z14");// atom B244 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B244 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B244

    atom.name="C"; atom.type=2;                                       // atom B245
    atom.fpos(1)=z14;atom.fpos(2)=x14;atom.fpos(3)=y14;                     // atom B245
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z14");atom.fpos_equation.push_back("x14");atom.fpos_equation.push_back("y14");// atom B245 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B245 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B245

    atom.name="C"; atom.type=2;                                       // atom B246
    atom.fpos(1)=((1.0/2.0)+z14);atom.fpos(2)=((1.0/2.0)-x14);atom.fpos(3)=-y14;                     // atom B246
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z14)");atom.fpos_equation.push_back("((1.0/2.0)-x14)");atom.fpos_equation.push_back("-y14");// atom B246 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B246 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B246

    atom.name="C"; atom.type=2;                                       // atom B247
    atom.fpos(1)=((1.0/2.0)-z14);atom.fpos(2)=-x14;atom.fpos(3)=((1.0/2.0)+y14);                     // atom B247
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z14)");atom.fpos_equation.push_back("-x14");atom.fpos_equation.push_back("((1.0/2.0)+y14)");// atom B247 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B247 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B247

    atom.name="C"; atom.type=2;                                       // atom B248
    atom.fpos(1)=-z14;atom.fpos(2)=((1.0/2.0)+x14);atom.fpos(3)=((1.0/2.0)-y14);                     // atom B248
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z14");atom.fpos_equation.push_back("((1.0/2.0)+x14)");atom.fpos_equation.push_back("((1.0/2.0)-y14)");// atom B248 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B248 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B248

    atom.name="C"; atom.type=2;                                       // atom B249
    atom.fpos(1)=y14;atom.fpos(2)=z14;atom.fpos(3)=x14;                     // atom B249
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y14");atom.fpos_equation.push_back("z14");atom.fpos_equation.push_back("x14");// atom B249 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B249 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B249

    atom.name="C"; atom.type=2;                                       // atom B250
    atom.fpos(1)=-y14;atom.fpos(2)=((1.0/2.0)+z14);atom.fpos(3)=((1.0/2.0)-x14);                     // atom B250
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y14");atom.fpos_equation.push_back("((1.0/2.0)+z14)");atom.fpos_equation.push_back("((1.0/2.0)-x14)");// atom B250 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B250 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B250

    atom.name="C"; atom.type=2;                                       // atom B251
    atom.fpos(1)=((1.0/2.0)+y14);atom.fpos(2)=((1.0/2.0)-z14);atom.fpos(3)=-x14;                     // atom B251
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y14)");atom.fpos_equation.push_back("((1.0/2.0)-z14)");atom.fpos_equation.push_back("-x14");// atom B251 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B251 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B251

    atom.name="C"; atom.type=2;                                       // atom B252
    atom.fpos(1)=((1.0/2.0)-y14);atom.fpos(2)=-z14;atom.fpos(3)=((1.0/2.0)+x14);                     // atom B252
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y14)");atom.fpos_equation.push_back("-z14");atom.fpos_equation.push_back("((1.0/2.0)+x14)");// atom B252 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B252 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B252

    atom.name="C"; atom.type=2;                                       // atom B253
    atom.fpos(1)=-x14;atom.fpos(2)=-y14;atom.fpos(3)=-z14;                     // atom B253
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x14");atom.fpos_equation.push_back("-y14");atom.fpos_equation.push_back("-z14");// atom B253 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B253 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B253

    atom.name="C"; atom.type=2;                                       // atom B254
    atom.fpos(1)=((1.0/2.0)+x14);atom.fpos(2)=y14;atom.fpos(3)=((1.0/2.0)-z14);                     // atom B254
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x14)");atom.fpos_equation.push_back("y14");atom.fpos_equation.push_back("((1.0/2.0)-z14)");// atom B254 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B254 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B254

    atom.name="C"; atom.type=2;                                       // atom B255
    atom.fpos(1)=x14;atom.fpos(2)=((1.0/2.0)-y14);atom.fpos(3)=((1.0/2.0)+z14);                     // atom B255
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x14");atom.fpos_equation.push_back("((1.0/2.0)-y14)");atom.fpos_equation.push_back("((1.0/2.0)+z14)");// atom B255 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B255 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B255

    atom.name="C"; atom.type=2;                                       // atom B256
    atom.fpos(1)=((1.0/2.0)-x14);atom.fpos(2)=((1.0/2.0)+y14);atom.fpos(3)=z14;                     // atom B256
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x14)");atom.fpos_equation.push_back("((1.0/2.0)+y14)");atom.fpos_equation.push_back("z14");// atom B256 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B256 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B256

    atom.name="C"; atom.type=2;                                       // atom B257
    atom.fpos(1)=-z14;atom.fpos(2)=-x14;atom.fpos(3)=-y14;                     // atom B257
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-z14");atom.fpos_equation.push_back("-x14");atom.fpos_equation.push_back("-y14");// atom B257 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B257 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B257

    atom.name="C"; atom.type=2;                                       // atom B258
    atom.fpos(1)=((1.0/2.0)-z14);atom.fpos(2)=((1.0/2.0)+x14);atom.fpos(3)=y14;                     // atom B258
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-z14)");atom.fpos_equation.push_back("((1.0/2.0)+x14)");atom.fpos_equation.push_back("y14");// atom B258 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B258 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B258

    atom.name="C"; atom.type=2;                                       // atom B259
    atom.fpos(1)=((1.0/2.0)+z14);atom.fpos(2)=x14;atom.fpos(3)=((1.0/2.0)-y14);                     // atom B259
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+z14)");atom.fpos_equation.push_back("x14");atom.fpos_equation.push_back("((1.0/2.0)-y14)");// atom B259 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B259 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B259

    atom.name="C"; atom.type=2;                                       // atom B260
    atom.fpos(1)=z14;atom.fpos(2)=((1.0/2.0)-x14);atom.fpos(3)=((1.0/2.0)+y14);                     // atom B260
    atom.fpos_equation.clear();atom.fpos_equation.push_back("z14");atom.fpos_equation.push_back("((1.0/2.0)-x14)");atom.fpos_equation.push_back("((1.0/2.0)+y14)");// atom B260 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B260 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B260

    atom.name="C"; atom.type=2;                                       // atom B261
    atom.fpos(1)=-y14;atom.fpos(2)=-z14;atom.fpos(3)=-x14;                     // atom B261
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y14");atom.fpos_equation.push_back("-z14");atom.fpos_equation.push_back("-x14");// atom B261 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B261 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B261

    atom.name="C"; atom.type=2;                                       // atom B262
    atom.fpos(1)=y14;atom.fpos(2)=((1.0/2.0)-z14);atom.fpos(3)=((1.0/2.0)+x14);                     // atom B262
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y14");atom.fpos_equation.push_back("((1.0/2.0)-z14)");atom.fpos_equation.push_back("((1.0/2.0)+x14)");// atom B262 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B262 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B262

    atom.name="C"; atom.type=2;                                       // atom B263
    atom.fpos(1)=((1.0/2.0)-y14);atom.fpos(2)=((1.0/2.0)+z14);atom.fpos(3)=x14;                     // atom B263
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-y14)");atom.fpos_equation.push_back("((1.0/2.0)+z14)");atom.fpos_equation.push_back("x14");// atom B263 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B263 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B263

    atom.name="C"; atom.type=2;                                       // atom B264
    atom.fpos(1)=((1.0/2.0)+y14);atom.fpos(2)=z14;atom.fpos(3)=((1.0/2.0)-x14);                     // atom B264
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+y14)");atom.fpos_equation.push_back("z14");atom.fpos_equation.push_back("((1.0/2.0)-x14)");// atom B264 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B264 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B264


    return str.atoms.size();  
  }
} // namespace anrl

namespace anrl {
  uint WebANRL_A2B3C6_cP264_205_2d_ab2c2d_6d(stringstream& web,bool LDEBUG) {
#ifndef _ANRL_NOWEB_
#endif

    if(LDEBUG) {cerr << "anrl:: WebANRL_A2B3C6_cP264_205_2d_ab2c2d_6d: web.str().size()=" << web.str().size() << endl;}

    return web.str().size();
  }
} // namespace anrl

#endif // AFLOW_REMOVE_GREP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

