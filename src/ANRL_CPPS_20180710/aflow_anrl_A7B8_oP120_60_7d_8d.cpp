// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - David Hicks - 2018
// FILE "ANRL/aflow_anrl_A7B8_oP120_60_7d_8d.cpp"

#ifndef _AFLOW_ANRL_A7B8_oP120_60_7d_8d_CPP // AFLOW_REMOVE_GREP
#define _AFLOW_ANRL_A7B8_oP120_60_7d_8d_CPP // AFLOW_REMOVE_GREP
#include "../aflow.h" // AFLOW_REMOVE_GREP

namespace anrl {
  uint WebANRL_A7B8_oP120_60_7d_8d(stringstream &web,bool LDEBUG);
}

namespace anrl {
  uint PrototypeANRL_A7B8_oP120_60_7d_8d(stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG) {
    // system A7B8_oP120_60_7d_8d

    if(XHOST.vflag_control.flag("WWW")) {
      WebANRL_A7B8_oP120_60_7d_8d(web,LDEBUG); // PLUG WEB STUFF
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

    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: FOUND" << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: label=" << label << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: nspecies=" << nspecies << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: natoms=" << natoms << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: spacegroup=" << spacegroup << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: nunderscores=" << nunderscores << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: nparameters=" <<  nparameters << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: Pearson_symbol=" << Pearson_symbol << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: params=" << params << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: Strukturbericht=" << Strukturbericht << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: prototype=" << prototype << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: dialect=" << dialect << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: vparameters.size()=" << vparameters.size() << endl;}

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
    double a=vparameters.at(i++);                  if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: a=" << a << endl;}
    double bovera=vparameters.at(i++),b=bovera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: b=" << b << " (b/a=" << bovera << ")" << endl;}
    double covera=vparameters.at(i++),c=covera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: c=" << c << " (c/a=" << covera << ")" << endl;}

    double x1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: x1=" << x1 << endl;}
    double y1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: y1=" << y1 << endl;}
    double z1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: z1=" << z1 << endl;}
    double x2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: x2=" << x2 << endl;}
    double y2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: y2=" << y2 << endl;}
    double z2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: z2=" << z2 << endl;}
    double x3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: x3=" << x3 << endl;}
    double y3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: y3=" << y3 << endl;}
    double z3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: z3=" << z3 << endl;}
    double x4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: x4=" << x4 << endl;}
    double y4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: y4=" << y4 << endl;}
    double z4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: z4=" << z4 << endl;}
    double x5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: x5=" << x5 << endl;}
    double y5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: y5=" << y5 << endl;}
    double z5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: z5=" << z5 << endl;}
    double x6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: x6=" << x6 << endl;}
    double y6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: y6=" << y6 << endl;}
    double z6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: z6=" << z6 << endl;}
    double x7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: x7=" << x7 << endl;}
    double y7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: y7=" << y7 << endl;}
    double z7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: z7=" << z7 << endl;}
    double x8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: x8=" << x8 << endl;}
    double y8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: y8=" << y8 << endl;}
    double z8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: z8=" << z8 << endl;}
    double x9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: x9=" << x9 << endl;}
    double y9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: y9=" << y9 << endl;}
    double z9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: z9=" << z9 << endl;}
    double x10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: x10=" << x10 << endl;}
    double y10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: y10=" << y10 << endl;}
    double z10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: z10=" << z10 << endl;}
    double x11=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: x11=" << x11 << endl;}
    double y11=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: y11=" << y11 << endl;}
    double z11=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: z11=" << z11 << endl;}
    double x12=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: x12=" << x12 << endl;}
    double y12=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: y12=" << y12 << endl;}
    double z12=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: z12=" << z12 << endl;}
    double x13=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: x13=" << x13 << endl;}
    double y13=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: y13=" << y13 << endl;}
    double z13=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: z13=" << z13 << endl;}
    double x14=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: x14=" << x14 << endl;}
    double y14=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: y14=" << y14 << endl;}
    double z14=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: z14=" << z14 << endl;}
    double x15=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: x15=" << x15 << endl;}
    double y15=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: y15=" << y15 << endl;}
    double z15=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A7B8_oP120_60_7d_8d: z15=" << z15 << endl;}

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
    atom.fpos(1)=x1;atom.fpos(2)=y1;atom.fpos(3)=z1;                     // atom B1
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x1");atom.fpos_equation.push_back("y1");atom.fpos_equation.push_back("z1");// atom B1 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B1 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B1

    atom.name="A"; atom.type=0;                                       // atom B2
    atom.fpos(1)=((1.0/2.0)-x1);atom.fpos(2)=((1.0/2.0)-y1);atom.fpos(3)=((1.0/2.0)+z1);                     // atom B2
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x1)");atom.fpos_equation.push_back("((1.0/2.0)-y1)");atom.fpos_equation.push_back("((1.0/2.0)+z1)");// atom B2 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B2 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B2

    atom.name="A"; atom.type=0;                                       // atom B3
    atom.fpos(1)=-x1;atom.fpos(2)=y1;atom.fpos(3)=((1.0/2.0)-z1);                     // atom B3
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x1");atom.fpos_equation.push_back("y1");atom.fpos_equation.push_back("((1.0/2.0)-z1)");// atom B3 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B3 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B3

    atom.name="A"; atom.type=0;                                       // atom B4
    atom.fpos(1)=((1.0/2.0)+x1);atom.fpos(2)=((1.0/2.0)-y1);atom.fpos(3)=-z1;                     // atom B4
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x1)");atom.fpos_equation.push_back("((1.0/2.0)-y1)");atom.fpos_equation.push_back("-z1");// atom B4 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B4 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B4

    atom.name="A"; atom.type=0;                                       // atom B5
    atom.fpos(1)=-x1;atom.fpos(2)=-y1;atom.fpos(3)=-z1;                     // atom B5
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x1");atom.fpos_equation.push_back("-y1");atom.fpos_equation.push_back("-z1");// atom B5 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B5 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B5

    atom.name="A"; atom.type=0;                                       // atom B6
    atom.fpos(1)=((1.0/2.0)+x1);atom.fpos(2)=((1.0/2.0)+y1);atom.fpos(3)=((1.0/2.0)-z1);                     // atom B6
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x1)");atom.fpos_equation.push_back("((1.0/2.0)+y1)");atom.fpos_equation.push_back("((1.0/2.0)-z1)");// atom B6 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B6 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B6

    atom.name="A"; atom.type=0;                                       // atom B7
    atom.fpos(1)=x1;atom.fpos(2)=-y1;atom.fpos(3)=((1.0/2.0)+z1);                     // atom B7
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x1");atom.fpos_equation.push_back("-y1");atom.fpos_equation.push_back("((1.0/2.0)+z1)");// atom B7 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B7 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B7

    atom.name="A"; atom.type=0;                                       // atom B8
    atom.fpos(1)=((1.0/2.0)-x1);atom.fpos(2)=((1.0/2.0)+y1);atom.fpos(3)=z1;                     // atom B8
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x1)");atom.fpos_equation.push_back("((1.0/2.0)+y1)");atom.fpos_equation.push_back("z1");// atom B8 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B8 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B8

    atom.name="A"; atom.type=0;                                       // atom B9
    atom.fpos(1)=x2;atom.fpos(2)=y2;atom.fpos(3)=z2;                     // atom B9
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x2");atom.fpos_equation.push_back("y2");atom.fpos_equation.push_back("z2");// atom B9 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B9 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B9

    atom.name="A"; atom.type=0;                                       // atom B10
    atom.fpos(1)=((1.0/2.0)-x2);atom.fpos(2)=((1.0/2.0)-y2);atom.fpos(3)=((1.0/2.0)+z2);                     // atom B10
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x2)");atom.fpos_equation.push_back("((1.0/2.0)-y2)");atom.fpos_equation.push_back("((1.0/2.0)+z2)");// atom B10 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B10 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B10

    atom.name="A"; atom.type=0;                                       // atom B11
    atom.fpos(1)=-x2;atom.fpos(2)=y2;atom.fpos(3)=((1.0/2.0)-z2);                     // atom B11
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x2");atom.fpos_equation.push_back("y2");atom.fpos_equation.push_back("((1.0/2.0)-z2)");// atom B11 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B11 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B11

    atom.name="A"; atom.type=0;                                       // atom B12
    atom.fpos(1)=((1.0/2.0)+x2);atom.fpos(2)=((1.0/2.0)-y2);atom.fpos(3)=-z2;                     // atom B12
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x2)");atom.fpos_equation.push_back("((1.0/2.0)-y2)");atom.fpos_equation.push_back("-z2");// atom B12 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B12 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B12

    atom.name="A"; atom.type=0;                                       // atom B13
    atom.fpos(1)=-x2;atom.fpos(2)=-y2;atom.fpos(3)=-z2;                     // atom B13
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x2");atom.fpos_equation.push_back("-y2");atom.fpos_equation.push_back("-z2");// atom B13 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B13 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B13

    atom.name="A"; atom.type=0;                                       // atom B14
    atom.fpos(1)=((1.0/2.0)+x2);atom.fpos(2)=((1.0/2.0)+y2);atom.fpos(3)=((1.0/2.0)-z2);                     // atom B14
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x2)");atom.fpos_equation.push_back("((1.0/2.0)+y2)");atom.fpos_equation.push_back("((1.0/2.0)-z2)");// atom B14 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B14 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B14

    atom.name="A"; atom.type=0;                                       // atom B15
    atom.fpos(1)=x2;atom.fpos(2)=-y2;atom.fpos(3)=((1.0/2.0)+z2);                     // atom B15
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x2");atom.fpos_equation.push_back("-y2");atom.fpos_equation.push_back("((1.0/2.0)+z2)");// atom B15 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B15 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B15

    atom.name="A"; atom.type=0;                                       // atom B16
    atom.fpos(1)=((1.0/2.0)-x2);atom.fpos(2)=((1.0/2.0)+y2);atom.fpos(3)=z2;                     // atom B16
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x2)");atom.fpos_equation.push_back("((1.0/2.0)+y2)");atom.fpos_equation.push_back("z2");// atom B16 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B16 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B16

    atom.name="A"; atom.type=0;                                       // atom B17
    atom.fpos(1)=x3;atom.fpos(2)=y3;atom.fpos(3)=z3;                     // atom B17
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("y3");atom.fpos_equation.push_back("z3");// atom B17 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B17 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B17

    atom.name="A"; atom.type=0;                                       // atom B18
    atom.fpos(1)=((1.0/2.0)-x3);atom.fpos(2)=((1.0/2.0)-y3);atom.fpos(3)=((1.0/2.0)+z3);                     // atom B18
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x3)");atom.fpos_equation.push_back("((1.0/2.0)-y3)");atom.fpos_equation.push_back("((1.0/2.0)+z3)");// atom B18 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B18 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B18

    atom.name="A"; atom.type=0;                                       // atom B19
    atom.fpos(1)=-x3;atom.fpos(2)=y3;atom.fpos(3)=((1.0/2.0)-z3);                     // atom B19
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("y3");atom.fpos_equation.push_back("((1.0/2.0)-z3)");// atom B19 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B19 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B19

    atom.name="A"; atom.type=0;                                       // atom B20
    atom.fpos(1)=((1.0/2.0)+x3);atom.fpos(2)=((1.0/2.0)-y3);atom.fpos(3)=-z3;                     // atom B20
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x3)");atom.fpos_equation.push_back("((1.0/2.0)-y3)");atom.fpos_equation.push_back("-z3");// atom B20 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B20 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B20

    atom.name="A"; atom.type=0;                                       // atom B21
    atom.fpos(1)=-x3;atom.fpos(2)=-y3;atom.fpos(3)=-z3;                     // atom B21
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("-y3");atom.fpos_equation.push_back("-z3");// atom B21 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B21 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B21

    atom.name="A"; atom.type=0;                                       // atom B22
    atom.fpos(1)=((1.0/2.0)+x3);atom.fpos(2)=((1.0/2.0)+y3);atom.fpos(3)=((1.0/2.0)-z3);                     // atom B22
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x3)");atom.fpos_equation.push_back("((1.0/2.0)+y3)");atom.fpos_equation.push_back("((1.0/2.0)-z3)");// atom B22 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B22 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B22

    atom.name="A"; atom.type=0;                                       // atom B23
    atom.fpos(1)=x3;atom.fpos(2)=-y3;atom.fpos(3)=((1.0/2.0)+z3);                     // atom B23
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("-y3");atom.fpos_equation.push_back("((1.0/2.0)+z3)");// atom B23 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B23 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B23

    atom.name="A"; atom.type=0;                                       // atom B24
    atom.fpos(1)=((1.0/2.0)-x3);atom.fpos(2)=((1.0/2.0)+y3);atom.fpos(3)=z3;                     // atom B24
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x3)");atom.fpos_equation.push_back("((1.0/2.0)+y3)");atom.fpos_equation.push_back("z3");// atom B24 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B24 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B24

    atom.name="A"; atom.type=0;                                       // atom B25
    atom.fpos(1)=x4;atom.fpos(2)=y4;atom.fpos(3)=z4;                     // atom B25
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x4");atom.fpos_equation.push_back("y4");atom.fpos_equation.push_back("z4");// atom B25 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B25 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B25

    atom.name="A"; atom.type=0;                                       // atom B26
    atom.fpos(1)=((1.0/2.0)-x4);atom.fpos(2)=((1.0/2.0)-y4);atom.fpos(3)=((1.0/2.0)+z4);                     // atom B26
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x4)");atom.fpos_equation.push_back("((1.0/2.0)-y4)");atom.fpos_equation.push_back("((1.0/2.0)+z4)");// atom B26 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B26 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B26

    atom.name="A"; atom.type=0;                                       // atom B27
    atom.fpos(1)=-x4;atom.fpos(2)=y4;atom.fpos(3)=((1.0/2.0)-z4);                     // atom B27
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x4");atom.fpos_equation.push_back("y4");atom.fpos_equation.push_back("((1.0/2.0)-z4)");// atom B27 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B27 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B27

    atom.name="A"; atom.type=0;                                       // atom B28
    atom.fpos(1)=((1.0/2.0)+x4);atom.fpos(2)=((1.0/2.0)-y4);atom.fpos(3)=-z4;                     // atom B28
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x4)");atom.fpos_equation.push_back("((1.0/2.0)-y4)");atom.fpos_equation.push_back("-z4");// atom B28 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B28 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B28

    atom.name="A"; atom.type=0;                                       // atom B29
    atom.fpos(1)=-x4;atom.fpos(2)=-y4;atom.fpos(3)=-z4;                     // atom B29
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x4");atom.fpos_equation.push_back("-y4");atom.fpos_equation.push_back("-z4");// atom B29 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B29 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B29

    atom.name="A"; atom.type=0;                                       // atom B30
    atom.fpos(1)=((1.0/2.0)+x4);atom.fpos(2)=((1.0/2.0)+y4);atom.fpos(3)=((1.0/2.0)-z4);                     // atom B30
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x4)");atom.fpos_equation.push_back("((1.0/2.0)+y4)");atom.fpos_equation.push_back("((1.0/2.0)-z4)");// atom B30 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B30 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B30

    atom.name="A"; atom.type=0;                                       // atom B31
    atom.fpos(1)=x4;atom.fpos(2)=-y4;atom.fpos(3)=((1.0/2.0)+z4);                     // atom B31
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x4");atom.fpos_equation.push_back("-y4");atom.fpos_equation.push_back("((1.0/2.0)+z4)");// atom B31 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B31 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B31

    atom.name="A"; atom.type=0;                                       // atom B32
    atom.fpos(1)=((1.0/2.0)-x4);atom.fpos(2)=((1.0/2.0)+y4);atom.fpos(3)=z4;                     // atom B32
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x4)");atom.fpos_equation.push_back("((1.0/2.0)+y4)");atom.fpos_equation.push_back("z4");// atom B32 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B32 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B32

    atom.name="A"; atom.type=0;                                       // atom B33
    atom.fpos(1)=x5;atom.fpos(2)=y5;atom.fpos(3)=z5;                     // atom B33
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("y5");atom.fpos_equation.push_back("z5");// atom B33 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B33 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B33

    atom.name="A"; atom.type=0;                                       // atom B34
    atom.fpos(1)=((1.0/2.0)-x5);atom.fpos(2)=((1.0/2.0)-y5);atom.fpos(3)=((1.0/2.0)+z5);                     // atom B34
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x5)");atom.fpos_equation.push_back("((1.0/2.0)-y5)");atom.fpos_equation.push_back("((1.0/2.0)+z5)");// atom B34 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B34 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B34

    atom.name="A"; atom.type=0;                                       // atom B35
    atom.fpos(1)=-x5;atom.fpos(2)=y5;atom.fpos(3)=((1.0/2.0)-z5);                     // atom B35
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("y5");atom.fpos_equation.push_back("((1.0/2.0)-z5)");// atom B35 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B35 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B35

    atom.name="A"; atom.type=0;                                       // atom B36
    atom.fpos(1)=((1.0/2.0)+x5);atom.fpos(2)=((1.0/2.0)-y5);atom.fpos(3)=-z5;                     // atom B36
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x5)");atom.fpos_equation.push_back("((1.0/2.0)-y5)");atom.fpos_equation.push_back("-z5");// atom B36 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B36 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B36

    atom.name="A"; atom.type=0;                                       // atom B37
    atom.fpos(1)=-x5;atom.fpos(2)=-y5;atom.fpos(3)=-z5;                     // atom B37
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("-y5");atom.fpos_equation.push_back("-z5");// atom B37 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B37 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B37

    atom.name="A"; atom.type=0;                                       // atom B38
    atom.fpos(1)=((1.0/2.0)+x5);atom.fpos(2)=((1.0/2.0)+y5);atom.fpos(3)=((1.0/2.0)-z5);                     // atom B38
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x5)");atom.fpos_equation.push_back("((1.0/2.0)+y5)");atom.fpos_equation.push_back("((1.0/2.0)-z5)");// atom B38 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B38 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B38

    atom.name="A"; atom.type=0;                                       // atom B39
    atom.fpos(1)=x5;atom.fpos(2)=-y5;atom.fpos(3)=((1.0/2.0)+z5);                     // atom B39
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("-y5");atom.fpos_equation.push_back("((1.0/2.0)+z5)");// atom B39 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B39 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B39

    atom.name="A"; atom.type=0;                                       // atom B40
    atom.fpos(1)=((1.0/2.0)-x5);atom.fpos(2)=((1.0/2.0)+y5);atom.fpos(3)=z5;                     // atom B40
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x5)");atom.fpos_equation.push_back("((1.0/2.0)+y5)");atom.fpos_equation.push_back("z5");// atom B40 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B40 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B40

    atom.name="A"; atom.type=0;                                       // atom B41
    atom.fpos(1)=x6;atom.fpos(2)=y6;atom.fpos(3)=z6;                     // atom B41
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("y6");atom.fpos_equation.push_back("z6");// atom B41 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B41 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B41

    atom.name="A"; atom.type=0;                                       // atom B42
    atom.fpos(1)=((1.0/2.0)-x6);atom.fpos(2)=((1.0/2.0)-y6);atom.fpos(3)=((1.0/2.0)+z6);                     // atom B42
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x6)");atom.fpos_equation.push_back("((1.0/2.0)-y6)");atom.fpos_equation.push_back("((1.0/2.0)+z6)");// atom B42 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B42 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B42

    atom.name="A"; atom.type=0;                                       // atom B43
    atom.fpos(1)=-x6;atom.fpos(2)=y6;atom.fpos(3)=((1.0/2.0)-z6);                     // atom B43
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("y6");atom.fpos_equation.push_back("((1.0/2.0)-z6)");// atom B43 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B43 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B43

    atom.name="A"; atom.type=0;                                       // atom B44
    atom.fpos(1)=((1.0/2.0)+x6);atom.fpos(2)=((1.0/2.0)-y6);atom.fpos(3)=-z6;                     // atom B44
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x6)");atom.fpos_equation.push_back("((1.0/2.0)-y6)");atom.fpos_equation.push_back("-z6");// atom B44 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B44 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B44

    atom.name="A"; atom.type=0;                                       // atom B45
    atom.fpos(1)=-x6;atom.fpos(2)=-y6;atom.fpos(3)=-z6;                     // atom B45
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("-y6");atom.fpos_equation.push_back("-z6");// atom B45 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B45 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B45

    atom.name="A"; atom.type=0;                                       // atom B46
    atom.fpos(1)=((1.0/2.0)+x6);atom.fpos(2)=((1.0/2.0)+y6);atom.fpos(3)=((1.0/2.0)-z6);                     // atom B46
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x6)");atom.fpos_equation.push_back("((1.0/2.0)+y6)");atom.fpos_equation.push_back("((1.0/2.0)-z6)");// atom B46 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B46 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B46

    atom.name="A"; atom.type=0;                                       // atom B47
    atom.fpos(1)=x6;atom.fpos(2)=-y6;atom.fpos(3)=((1.0/2.0)+z6);                     // atom B47
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("-y6");atom.fpos_equation.push_back("((1.0/2.0)+z6)");// atom B47 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B47 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B47

    atom.name="A"; atom.type=0;                                       // atom B48
    atom.fpos(1)=((1.0/2.0)-x6);atom.fpos(2)=((1.0/2.0)+y6);atom.fpos(3)=z6;                     // atom B48
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x6)");atom.fpos_equation.push_back("((1.0/2.0)+y6)");atom.fpos_equation.push_back("z6");// atom B48 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B48 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B48

    atom.name="A"; atom.type=0;                                       // atom B49
    atom.fpos(1)=x7;atom.fpos(2)=y7;atom.fpos(3)=z7;                     // atom B49
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x7");atom.fpos_equation.push_back("y7");atom.fpos_equation.push_back("z7");// atom B49 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B49 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B49

    atom.name="A"; atom.type=0;                                       // atom B50
    atom.fpos(1)=((1.0/2.0)-x7);atom.fpos(2)=((1.0/2.0)-y7);atom.fpos(3)=((1.0/2.0)+z7);                     // atom B50
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x7)");atom.fpos_equation.push_back("((1.0/2.0)-y7)");atom.fpos_equation.push_back("((1.0/2.0)+z7)");// atom B50 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B50 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B50

    atom.name="A"; atom.type=0;                                       // atom B51
    atom.fpos(1)=-x7;atom.fpos(2)=y7;atom.fpos(3)=((1.0/2.0)-z7);                     // atom B51
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x7");atom.fpos_equation.push_back("y7");atom.fpos_equation.push_back("((1.0/2.0)-z7)");// atom B51 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B51 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B51

    atom.name="A"; atom.type=0;                                       // atom B52
    atom.fpos(1)=((1.0/2.0)+x7);atom.fpos(2)=((1.0/2.0)-y7);atom.fpos(3)=-z7;                     // atom B52
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x7)");atom.fpos_equation.push_back("((1.0/2.0)-y7)");atom.fpos_equation.push_back("-z7");// atom B52 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B52 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B52

    atom.name="A"; atom.type=0;                                       // atom B53
    atom.fpos(1)=-x7;atom.fpos(2)=-y7;atom.fpos(3)=-z7;                     // atom B53
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x7");atom.fpos_equation.push_back("-y7");atom.fpos_equation.push_back("-z7");// atom B53 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B53 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B53

    atom.name="A"; atom.type=0;                                       // atom B54
    atom.fpos(1)=((1.0/2.0)+x7);atom.fpos(2)=((1.0/2.0)+y7);atom.fpos(3)=((1.0/2.0)-z7);                     // atom B54
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x7)");atom.fpos_equation.push_back("((1.0/2.0)+y7)");atom.fpos_equation.push_back("((1.0/2.0)-z7)");// atom B54 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B54 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B54

    atom.name="A"; atom.type=0;                                       // atom B55
    atom.fpos(1)=x7;atom.fpos(2)=-y7;atom.fpos(3)=((1.0/2.0)+z7);                     // atom B55
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x7");atom.fpos_equation.push_back("-y7");atom.fpos_equation.push_back("((1.0/2.0)+z7)");// atom B55 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B55 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B55

    atom.name="A"; atom.type=0;                                       // atom B56
    atom.fpos(1)=((1.0/2.0)-x7);atom.fpos(2)=((1.0/2.0)+y7);atom.fpos(3)=z7;                     // atom B56
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x7)");atom.fpos_equation.push_back("((1.0/2.0)+y7)");atom.fpos_equation.push_back("z7");// atom B56 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B56 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B56

    atom.name="B"; atom.type=1;                                       // atom B57
    atom.fpos(1)=x8;atom.fpos(2)=y8;atom.fpos(3)=z8;                     // atom B57
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x8");atom.fpos_equation.push_back("y8");atom.fpos_equation.push_back("z8");// atom B57 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B57 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B57

    atom.name="B"; atom.type=1;                                       // atom B58
    atom.fpos(1)=((1.0/2.0)-x8);atom.fpos(2)=((1.0/2.0)-y8);atom.fpos(3)=((1.0/2.0)+z8);                     // atom B58
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x8)");atom.fpos_equation.push_back("((1.0/2.0)-y8)");atom.fpos_equation.push_back("((1.0/2.0)+z8)");// atom B58 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B58 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B58

    atom.name="B"; atom.type=1;                                       // atom B59
    atom.fpos(1)=-x8;atom.fpos(2)=y8;atom.fpos(3)=((1.0/2.0)-z8);                     // atom B59
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x8");atom.fpos_equation.push_back("y8");atom.fpos_equation.push_back("((1.0/2.0)-z8)");// atom B59 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B59 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B59

    atom.name="B"; atom.type=1;                                       // atom B60
    atom.fpos(1)=((1.0/2.0)+x8);atom.fpos(2)=((1.0/2.0)-y8);atom.fpos(3)=-z8;                     // atom B60
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x8)");atom.fpos_equation.push_back("((1.0/2.0)-y8)");atom.fpos_equation.push_back("-z8");// atom B60 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B60 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B60

    atom.name="B"; atom.type=1;                                       // atom B61
    atom.fpos(1)=-x8;atom.fpos(2)=-y8;atom.fpos(3)=-z8;                     // atom B61
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x8");atom.fpos_equation.push_back("-y8");atom.fpos_equation.push_back("-z8");// atom B61 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B61 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B61

    atom.name="B"; atom.type=1;                                       // atom B62
    atom.fpos(1)=((1.0/2.0)+x8);atom.fpos(2)=((1.0/2.0)+y8);atom.fpos(3)=((1.0/2.0)-z8);                     // atom B62
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x8)");atom.fpos_equation.push_back("((1.0/2.0)+y8)");atom.fpos_equation.push_back("((1.0/2.0)-z8)");// atom B62 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B62 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B62

    atom.name="B"; atom.type=1;                                       // atom B63
    atom.fpos(1)=x8;atom.fpos(2)=-y8;atom.fpos(3)=((1.0/2.0)+z8);                     // atom B63
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x8");atom.fpos_equation.push_back("-y8");atom.fpos_equation.push_back("((1.0/2.0)+z8)");// atom B63 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B63 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B63

    atom.name="B"; atom.type=1;                                       // atom B64
    atom.fpos(1)=((1.0/2.0)-x8);atom.fpos(2)=((1.0/2.0)+y8);atom.fpos(3)=z8;                     // atom B64
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x8)");atom.fpos_equation.push_back("((1.0/2.0)+y8)");atom.fpos_equation.push_back("z8");// atom B64 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B64 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B64

    atom.name="B"; atom.type=1;                                       // atom B65
    atom.fpos(1)=x9;atom.fpos(2)=y9;atom.fpos(3)=z9;                     // atom B65
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x9");atom.fpos_equation.push_back("y9");atom.fpos_equation.push_back("z9");// atom B65 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B65 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B65

    atom.name="B"; atom.type=1;                                       // atom B66
    atom.fpos(1)=((1.0/2.0)-x9);atom.fpos(2)=((1.0/2.0)-y9);atom.fpos(3)=((1.0/2.0)+z9);                     // atom B66
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x9)");atom.fpos_equation.push_back("((1.0/2.0)-y9)");atom.fpos_equation.push_back("((1.0/2.0)+z9)");// atom B66 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B66 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B66

    atom.name="B"; atom.type=1;                                       // atom B67
    atom.fpos(1)=-x9;atom.fpos(2)=y9;atom.fpos(3)=((1.0/2.0)-z9);                     // atom B67
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x9");atom.fpos_equation.push_back("y9");atom.fpos_equation.push_back("((1.0/2.0)-z9)");// atom B67 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B67 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B67

    atom.name="B"; atom.type=1;                                       // atom B68
    atom.fpos(1)=((1.0/2.0)+x9);atom.fpos(2)=((1.0/2.0)-y9);atom.fpos(3)=-z9;                     // atom B68
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x9)");atom.fpos_equation.push_back("((1.0/2.0)-y9)");atom.fpos_equation.push_back("-z9");// atom B68 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B68 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B68

    atom.name="B"; atom.type=1;                                       // atom B69
    atom.fpos(1)=-x9;atom.fpos(2)=-y9;atom.fpos(3)=-z9;                     // atom B69
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x9");atom.fpos_equation.push_back("-y9");atom.fpos_equation.push_back("-z9");// atom B69 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B69 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B69

    atom.name="B"; atom.type=1;                                       // atom B70
    atom.fpos(1)=((1.0/2.0)+x9);atom.fpos(2)=((1.0/2.0)+y9);atom.fpos(3)=((1.0/2.0)-z9);                     // atom B70
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x9)");atom.fpos_equation.push_back("((1.0/2.0)+y9)");atom.fpos_equation.push_back("((1.0/2.0)-z9)");// atom B70 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B70 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B70

    atom.name="B"; atom.type=1;                                       // atom B71
    atom.fpos(1)=x9;atom.fpos(2)=-y9;atom.fpos(3)=((1.0/2.0)+z9);                     // atom B71
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x9");atom.fpos_equation.push_back("-y9");atom.fpos_equation.push_back("((1.0/2.0)+z9)");// atom B71 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B71 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B71

    atom.name="B"; atom.type=1;                                       // atom B72
    atom.fpos(1)=((1.0/2.0)-x9);atom.fpos(2)=((1.0/2.0)+y9);atom.fpos(3)=z9;                     // atom B72
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x9)");atom.fpos_equation.push_back("((1.0/2.0)+y9)");atom.fpos_equation.push_back("z9");// atom B72 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B72 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B72

    atom.name="B"; atom.type=1;                                       // atom B73
    atom.fpos(1)=x10;atom.fpos(2)=y10;atom.fpos(3)=z10;                     // atom B73
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x10");atom.fpos_equation.push_back("y10");atom.fpos_equation.push_back("z10");// atom B73 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B73 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B73

    atom.name="B"; atom.type=1;                                       // atom B74
    atom.fpos(1)=((1.0/2.0)-x10);atom.fpos(2)=((1.0/2.0)-y10);atom.fpos(3)=((1.0/2.0)+z10);                     // atom B74
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x10)");atom.fpos_equation.push_back("((1.0/2.0)-y10)");atom.fpos_equation.push_back("((1.0/2.0)+z10)");// atom B74 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B74 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B74

    atom.name="B"; atom.type=1;                                       // atom B75
    atom.fpos(1)=-x10;atom.fpos(2)=y10;atom.fpos(3)=((1.0/2.0)-z10);                     // atom B75
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x10");atom.fpos_equation.push_back("y10");atom.fpos_equation.push_back("((1.0/2.0)-z10)");// atom B75 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B75 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B75

    atom.name="B"; atom.type=1;                                       // atom B76
    atom.fpos(1)=((1.0/2.0)+x10);atom.fpos(2)=((1.0/2.0)-y10);atom.fpos(3)=-z10;                     // atom B76
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x10)");atom.fpos_equation.push_back("((1.0/2.0)-y10)");atom.fpos_equation.push_back("-z10");// atom B76 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B76 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B76

    atom.name="B"; atom.type=1;                                       // atom B77
    atom.fpos(1)=-x10;atom.fpos(2)=-y10;atom.fpos(3)=-z10;                     // atom B77
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x10");atom.fpos_equation.push_back("-y10");atom.fpos_equation.push_back("-z10");// atom B77 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B77 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B77

    atom.name="B"; atom.type=1;                                       // atom B78
    atom.fpos(1)=((1.0/2.0)+x10);atom.fpos(2)=((1.0/2.0)+y10);atom.fpos(3)=((1.0/2.0)-z10);                     // atom B78
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x10)");atom.fpos_equation.push_back("((1.0/2.0)+y10)");atom.fpos_equation.push_back("((1.0/2.0)-z10)");// atom B78 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B78 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B78

    atom.name="B"; atom.type=1;                                       // atom B79
    atom.fpos(1)=x10;atom.fpos(2)=-y10;atom.fpos(3)=((1.0/2.0)+z10);                     // atom B79
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x10");atom.fpos_equation.push_back("-y10");atom.fpos_equation.push_back("((1.0/2.0)+z10)");// atom B79 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B79 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B79

    atom.name="B"; atom.type=1;                                       // atom B80
    atom.fpos(1)=((1.0/2.0)-x10);atom.fpos(2)=((1.0/2.0)+y10);atom.fpos(3)=z10;                     // atom B80
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x10)");atom.fpos_equation.push_back("((1.0/2.0)+y10)");atom.fpos_equation.push_back("z10");// atom B80 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B80 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B80

    atom.name="B"; atom.type=1;                                       // atom B81
    atom.fpos(1)=x11;atom.fpos(2)=y11;atom.fpos(3)=z11;                     // atom B81
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x11");atom.fpos_equation.push_back("y11");atom.fpos_equation.push_back("z11");// atom B81 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B81 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B81

    atom.name="B"; atom.type=1;                                       // atom B82
    atom.fpos(1)=((1.0/2.0)-x11);atom.fpos(2)=((1.0/2.0)-y11);atom.fpos(3)=((1.0/2.0)+z11);                     // atom B82
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x11)");atom.fpos_equation.push_back("((1.0/2.0)-y11)");atom.fpos_equation.push_back("((1.0/2.0)+z11)");// atom B82 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B82 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B82

    atom.name="B"; atom.type=1;                                       // atom B83
    atom.fpos(1)=-x11;atom.fpos(2)=y11;atom.fpos(3)=((1.0/2.0)-z11);                     // atom B83
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x11");atom.fpos_equation.push_back("y11");atom.fpos_equation.push_back("((1.0/2.0)-z11)");// atom B83 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B83 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B83

    atom.name="B"; atom.type=1;                                       // atom B84
    atom.fpos(1)=((1.0/2.0)+x11);atom.fpos(2)=((1.0/2.0)-y11);atom.fpos(3)=-z11;                     // atom B84
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x11)");atom.fpos_equation.push_back("((1.0/2.0)-y11)");atom.fpos_equation.push_back("-z11");// atom B84 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B84 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B84

    atom.name="B"; atom.type=1;                                       // atom B85
    atom.fpos(1)=-x11;atom.fpos(2)=-y11;atom.fpos(3)=-z11;                     // atom B85
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x11");atom.fpos_equation.push_back("-y11");atom.fpos_equation.push_back("-z11");// atom B85 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B85 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B85

    atom.name="B"; atom.type=1;                                       // atom B86
    atom.fpos(1)=((1.0/2.0)+x11);atom.fpos(2)=((1.0/2.0)+y11);atom.fpos(3)=((1.0/2.0)-z11);                     // atom B86
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x11)");atom.fpos_equation.push_back("((1.0/2.0)+y11)");atom.fpos_equation.push_back("((1.0/2.0)-z11)");// atom B86 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B86 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B86

    atom.name="B"; atom.type=1;                                       // atom B87
    atom.fpos(1)=x11;atom.fpos(2)=-y11;atom.fpos(3)=((1.0/2.0)+z11);                     // atom B87
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x11");atom.fpos_equation.push_back("-y11");atom.fpos_equation.push_back("((1.0/2.0)+z11)");// atom B87 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B87 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B87

    atom.name="B"; atom.type=1;                                       // atom B88
    atom.fpos(1)=((1.0/2.0)-x11);atom.fpos(2)=((1.0/2.0)+y11);atom.fpos(3)=z11;                     // atom B88
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x11)");atom.fpos_equation.push_back("((1.0/2.0)+y11)");atom.fpos_equation.push_back("z11");// atom B88 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B88 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B88

    atom.name="B"; atom.type=1;                                       // atom B89
    atom.fpos(1)=x12;atom.fpos(2)=y12;atom.fpos(3)=z12;                     // atom B89
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x12");atom.fpos_equation.push_back("y12");atom.fpos_equation.push_back("z12");// atom B89 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B89 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B89

    atom.name="B"; atom.type=1;                                       // atom B90
    atom.fpos(1)=((1.0/2.0)-x12);atom.fpos(2)=((1.0/2.0)-y12);atom.fpos(3)=((1.0/2.0)+z12);                     // atom B90
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x12)");atom.fpos_equation.push_back("((1.0/2.0)-y12)");atom.fpos_equation.push_back("((1.0/2.0)+z12)");// atom B90 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B90 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B90

    atom.name="B"; atom.type=1;                                       // atom B91
    atom.fpos(1)=-x12;atom.fpos(2)=y12;atom.fpos(3)=((1.0/2.0)-z12);                     // atom B91
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x12");atom.fpos_equation.push_back("y12");atom.fpos_equation.push_back("((1.0/2.0)-z12)");// atom B91 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B91 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B91

    atom.name="B"; atom.type=1;                                       // atom B92
    atom.fpos(1)=((1.0/2.0)+x12);atom.fpos(2)=((1.0/2.0)-y12);atom.fpos(3)=-z12;                     // atom B92
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x12)");atom.fpos_equation.push_back("((1.0/2.0)-y12)");atom.fpos_equation.push_back("-z12");// atom B92 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B92 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B92

    atom.name="B"; atom.type=1;                                       // atom B93
    atom.fpos(1)=-x12;atom.fpos(2)=-y12;atom.fpos(3)=-z12;                     // atom B93
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x12");atom.fpos_equation.push_back("-y12");atom.fpos_equation.push_back("-z12");// atom B93 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B93 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B93

    atom.name="B"; atom.type=1;                                       // atom B94
    atom.fpos(1)=((1.0/2.0)+x12);atom.fpos(2)=((1.0/2.0)+y12);atom.fpos(3)=((1.0/2.0)-z12);                     // atom B94
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x12)");atom.fpos_equation.push_back("((1.0/2.0)+y12)");atom.fpos_equation.push_back("((1.0/2.0)-z12)");// atom B94 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B94 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B94

    atom.name="B"; atom.type=1;                                       // atom B95
    atom.fpos(1)=x12;atom.fpos(2)=-y12;atom.fpos(3)=((1.0/2.0)+z12);                     // atom B95
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x12");atom.fpos_equation.push_back("-y12");atom.fpos_equation.push_back("((1.0/2.0)+z12)");// atom B95 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B95 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B95

    atom.name="B"; atom.type=1;                                       // atom B96
    atom.fpos(1)=((1.0/2.0)-x12);atom.fpos(2)=((1.0/2.0)+y12);atom.fpos(3)=z12;                     // atom B96
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x12)");atom.fpos_equation.push_back("((1.0/2.0)+y12)");atom.fpos_equation.push_back("z12");// atom B96 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B96 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B96

    atom.name="B"; atom.type=1;                                       // atom B97
    atom.fpos(1)=x13;atom.fpos(2)=y13;atom.fpos(3)=z13;                     // atom B97
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x13");atom.fpos_equation.push_back("y13");atom.fpos_equation.push_back("z13");// atom B97 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B97 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B97

    atom.name="B"; atom.type=1;                                       // atom B98
    atom.fpos(1)=((1.0/2.0)-x13);atom.fpos(2)=((1.0/2.0)-y13);atom.fpos(3)=((1.0/2.0)+z13);                     // atom B98
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x13)");atom.fpos_equation.push_back("((1.0/2.0)-y13)");atom.fpos_equation.push_back("((1.0/2.0)+z13)");// atom B98 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B98 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B98

    atom.name="B"; atom.type=1;                                       // atom B99
    atom.fpos(1)=-x13;atom.fpos(2)=y13;atom.fpos(3)=((1.0/2.0)-z13);                     // atom B99
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x13");atom.fpos_equation.push_back("y13");atom.fpos_equation.push_back("((1.0/2.0)-z13)");// atom B99 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B99 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B99

    atom.name="B"; atom.type=1;                                       // atom B100
    atom.fpos(1)=((1.0/2.0)+x13);atom.fpos(2)=((1.0/2.0)-y13);atom.fpos(3)=-z13;                     // atom B100
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x13)");atom.fpos_equation.push_back("((1.0/2.0)-y13)");atom.fpos_equation.push_back("-z13");// atom B100 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B100 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B100

    atom.name="B"; atom.type=1;                                       // atom B101
    atom.fpos(1)=-x13;atom.fpos(2)=-y13;atom.fpos(3)=-z13;                     // atom B101
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x13");atom.fpos_equation.push_back("-y13");atom.fpos_equation.push_back("-z13");// atom B101 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B101 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B101

    atom.name="B"; atom.type=1;                                       // atom B102
    atom.fpos(1)=((1.0/2.0)+x13);atom.fpos(2)=((1.0/2.0)+y13);atom.fpos(3)=((1.0/2.0)-z13);                     // atom B102
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x13)");atom.fpos_equation.push_back("((1.0/2.0)+y13)");atom.fpos_equation.push_back("((1.0/2.0)-z13)");// atom B102 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B102 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B102

    atom.name="B"; atom.type=1;                                       // atom B103
    atom.fpos(1)=x13;atom.fpos(2)=-y13;atom.fpos(3)=((1.0/2.0)+z13);                     // atom B103
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x13");atom.fpos_equation.push_back("-y13");atom.fpos_equation.push_back("((1.0/2.0)+z13)");// atom B103 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B103 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B103

    atom.name="B"; atom.type=1;                                       // atom B104
    atom.fpos(1)=((1.0/2.0)-x13);atom.fpos(2)=((1.0/2.0)+y13);atom.fpos(3)=z13;                     // atom B104
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x13)");atom.fpos_equation.push_back("((1.0/2.0)+y13)");atom.fpos_equation.push_back("z13");// atom B104 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B104 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B104

    atom.name="B"; atom.type=1;                                       // atom B105
    atom.fpos(1)=x14;atom.fpos(2)=y14;atom.fpos(3)=z14;                     // atom B105
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x14");atom.fpos_equation.push_back("y14");atom.fpos_equation.push_back("z14");// atom B105 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B105 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B105

    atom.name="B"; atom.type=1;                                       // atom B106
    atom.fpos(1)=((1.0/2.0)-x14);atom.fpos(2)=((1.0/2.0)-y14);atom.fpos(3)=((1.0/2.0)+z14);                     // atom B106
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x14)");atom.fpos_equation.push_back("((1.0/2.0)-y14)");atom.fpos_equation.push_back("((1.0/2.0)+z14)");// atom B106 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B106 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B106

    atom.name="B"; atom.type=1;                                       // atom B107
    atom.fpos(1)=-x14;atom.fpos(2)=y14;atom.fpos(3)=((1.0/2.0)-z14);                     // atom B107
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x14");atom.fpos_equation.push_back("y14");atom.fpos_equation.push_back("((1.0/2.0)-z14)");// atom B107 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B107 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B107

    atom.name="B"; atom.type=1;                                       // atom B108
    atom.fpos(1)=((1.0/2.0)+x14);atom.fpos(2)=((1.0/2.0)-y14);atom.fpos(3)=-z14;                     // atom B108
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x14)");atom.fpos_equation.push_back("((1.0/2.0)-y14)");atom.fpos_equation.push_back("-z14");// atom B108 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B108 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B108

    atom.name="B"; atom.type=1;                                       // atom B109
    atom.fpos(1)=-x14;atom.fpos(2)=-y14;atom.fpos(3)=-z14;                     // atom B109
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x14");atom.fpos_equation.push_back("-y14");atom.fpos_equation.push_back("-z14");// atom B109 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B109 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B109

    atom.name="B"; atom.type=1;                                       // atom B110
    atom.fpos(1)=((1.0/2.0)+x14);atom.fpos(2)=((1.0/2.0)+y14);atom.fpos(3)=((1.0/2.0)-z14);                     // atom B110
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x14)");atom.fpos_equation.push_back("((1.0/2.0)+y14)");atom.fpos_equation.push_back("((1.0/2.0)-z14)");// atom B110 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B110 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B110

    atom.name="B"; atom.type=1;                                       // atom B111
    atom.fpos(1)=x14;atom.fpos(2)=-y14;atom.fpos(3)=((1.0/2.0)+z14);                     // atom B111
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x14");atom.fpos_equation.push_back("-y14");atom.fpos_equation.push_back("((1.0/2.0)+z14)");// atom B111 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B111 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B111

    atom.name="B"; atom.type=1;                                       // atom B112
    atom.fpos(1)=((1.0/2.0)-x14);atom.fpos(2)=((1.0/2.0)+y14);atom.fpos(3)=z14;                     // atom B112
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x14)");atom.fpos_equation.push_back("((1.0/2.0)+y14)");atom.fpos_equation.push_back("z14");// atom B112 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B112 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B112

    atom.name="B"; atom.type=1;                                       // atom B113
    atom.fpos(1)=x15;atom.fpos(2)=y15;atom.fpos(3)=z15;                     // atom B113
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x15");atom.fpos_equation.push_back("y15");atom.fpos_equation.push_back("z15");// atom B113 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B113 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B113

    atom.name="B"; atom.type=1;                                       // atom B114
    atom.fpos(1)=((1.0/2.0)-x15);atom.fpos(2)=((1.0/2.0)-y15);atom.fpos(3)=((1.0/2.0)+z15);                     // atom B114
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x15)");atom.fpos_equation.push_back("((1.0/2.0)-y15)");atom.fpos_equation.push_back("((1.0/2.0)+z15)");// atom B114 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B114 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B114

    atom.name="B"; atom.type=1;                                       // atom B115
    atom.fpos(1)=-x15;atom.fpos(2)=y15;atom.fpos(3)=((1.0/2.0)-z15);                     // atom B115
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x15");atom.fpos_equation.push_back("y15");atom.fpos_equation.push_back("((1.0/2.0)-z15)");// atom B115 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B115 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B115

    atom.name="B"; atom.type=1;                                       // atom B116
    atom.fpos(1)=((1.0/2.0)+x15);atom.fpos(2)=((1.0/2.0)-y15);atom.fpos(3)=-z15;                     // atom B116
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x15)");atom.fpos_equation.push_back("((1.0/2.0)-y15)");atom.fpos_equation.push_back("-z15");// atom B116 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B116 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B116

    atom.name="B"; atom.type=1;                                       // atom B117
    atom.fpos(1)=-x15;atom.fpos(2)=-y15;atom.fpos(3)=-z15;                     // atom B117
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x15");atom.fpos_equation.push_back("-y15");atom.fpos_equation.push_back("-z15");// atom B117 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B117 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B117

    atom.name="B"; atom.type=1;                                       // atom B118
    atom.fpos(1)=((1.0/2.0)+x15);atom.fpos(2)=((1.0/2.0)+y15);atom.fpos(3)=((1.0/2.0)-z15);                     // atom B118
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)+x15)");atom.fpos_equation.push_back("((1.0/2.0)+y15)");atom.fpos_equation.push_back("((1.0/2.0)-z15)");// atom B118 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B118 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B118

    atom.name="B"; atom.type=1;                                       // atom B119
    atom.fpos(1)=x15;atom.fpos(2)=-y15;atom.fpos(3)=((1.0/2.0)+z15);                     // atom B119
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x15");atom.fpos_equation.push_back("-y15");atom.fpos_equation.push_back("((1.0/2.0)+z15)");// atom B119 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B119 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B119

    atom.name="B"; atom.type=1;                                       // atom B120
    atom.fpos(1)=((1.0/2.0)-x15);atom.fpos(2)=((1.0/2.0)+y15);atom.fpos(3)=z15;                     // atom B120
    atom.fpos_equation.clear();atom.fpos_equation.push_back("((1.0/2.0)-x15)");atom.fpos_equation.push_back("((1.0/2.0)+y15)");atom.fpos_equation.push_back("z15");// atom B120 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B120 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B120


    return str.atoms.size();  
  }
} // namespace anrl

namespace anrl {
  uint WebANRL_A7B8_oP120_60_7d_8d(stringstream& web,bool LDEBUG) {
#ifndef _ANRL_NOWEB_
#endif

    if(LDEBUG) {cerr << "anrl:: WebANRL_A7B8_oP120_60_7d_8d: web.str().size()=" << web.str().size() << endl;}

    return web.str().size();
  }
} // namespace anrl

#endif // AFLOW_REMOVE_GREP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

