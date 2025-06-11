// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - David Hicks - 2020
// FILE "ANRL/aflow_anrl_A3B5C_mC54_8_3a3b_9a3b_3a.cpp"

#ifndef _AFLOW_ANRL_A3B5C_mC54_8_3a3b_9a3b_3a_CPP // AFLOW_REMOVE_GREP
#define _AFLOW_ANRL_A3B5C_mC54_8_3a3b_9a3b_3a_CPP // AFLOW_REMOVE_GREP
#include "../aflow.h" // AFLOW_REMOVE_GREP

namespace anrl {
  uint WebANRL_A3B5C_mC54_8_3a3b_9a3b_3a(stringstream &web,bool LDEBUG);
}

namespace anrl {
  uint PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a(stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG) {
    // system A3B5C_mC54_8_3a3b_9a3b_3a

    if(XHOST.vflag_control.flag("WWW")) {
      WebANRL_A3B5C_mC54_8_3a3b_9a3b_3a(web,LDEBUG); // PLUG WEB STUFF
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

    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: FOUND" << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: label=" << label << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: nspecies=" << nspecies << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: natoms=" << natoms << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: spacegroup=" << spacegroup << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: nunderscores=" << nunderscores << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: nparameters=" <<  nparameters << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: Pearson_symbol=" << Pearson_symbol << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: params=" << params << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: Strukturbericht=" << Strukturbericht << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: prototype=" << prototype << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: dialect=" << dialect << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: vparameters.size()=" << vparameters.size() << endl;}

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
    double a=vparameters.at(i++);                  if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: a=" << a << endl;}
    double bovera=vparameters.at(i++),b=bovera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: b=" << b << " (b/a=" << bovera << ")" << endl;}
    double covera=vparameters.at(i++),c=covera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: c=" << c << " (c/a=" << covera << ")" << endl;}
    double beta=vparameters.at(i++);               if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: beta=" << beta << endl;}

    double x1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: x1=" << x1 << endl;}
    double z1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: z1=" << z1 << endl;}
    double x2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: x2=" << x2 << endl;}
    double z2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: z2=" << z2 << endl;}
    double x3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: x3=" << x3 << endl;}
    double z3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: z3=" << z3 << endl;}
    double x4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: x4=" << x4 << endl;}
    double z4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: z4=" << z4 << endl;}
    double x5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: x5=" << x5 << endl;}
    double z5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: z5=" << z5 << endl;}
    double x6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: x6=" << x6 << endl;}
    double z6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: z6=" << z6 << endl;}
    double x7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: x7=" << x7 << endl;}
    double z7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: z7=" << z7 << endl;}
    double x8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: x8=" << x8 << endl;}
    double z8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: z8=" << z8 << endl;}
    double x9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: x9=" << x9 << endl;}
    double z9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: z9=" << z9 << endl;}
    double x10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: x10=" << x10 << endl;}
    double z10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: z10=" << z10 << endl;}
    double x11=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: x11=" << x11 << endl;}
    double z11=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: z11=" << z11 << endl;}
    double x12=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: x12=" << x12 << endl;}
    double z12=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: z12=" << z12 << endl;}
    double x13=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: x13=" << x13 << endl;}
    double z13=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: z13=" << z13 << endl;}
    double x14=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: x14=" << x14 << endl;}
    double z14=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: z14=" << z14 << endl;}
    double x15=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: x15=" << x15 << endl;}
    double z15=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: z15=" << z15 << endl;}
    double x16=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: x16=" << x16 << endl;}
    double y16=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: y16=" << y16 << endl;}
    double z16=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: z16=" << z16 << endl;}
    double x17=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: x17=" << x17 << endl;}
    double y17=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: y17=" << y17 << endl;}
    double z17=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: z17=" << z17 << endl;}
    double x18=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: x18=" << x18 << endl;}
    double y18=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: y18=" << y18 << endl;}
    double z18=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: z18=" << z18 << endl;}
    double x19=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: x19=" << x19 << endl;}
    double y19=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: y19=" << y19 << endl;}
    double z19=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: z19=" << z19 << endl;}
    double x20=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: x20=" << x20 << endl;}
    double y20=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: y20=" << y20 << endl;}
    double z20=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: z20=" << z20 << endl;}
    double x21=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: x21=" << x21 << endl;}
    double y21=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: y21=" << y21 << endl;}
    double z21=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: z21=" << z21 << endl;}

    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: cos(beta)=" << cos(deg2rad*beta)  << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3B5C_mC54_8_3a3b_9a3b_3a: sin(beta)=" << sin(deg2rad*beta)  << endl;}

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
    atom.fpos(1)=x1;atom.fpos(2)=x1;atom.fpos(3)=z1;                     // atom B1
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x1");atom.fpos_equation.push_back("x1");atom.fpos_equation.push_back("z1");// atom B1 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B1 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B1

    atom.name="A"; atom.type=0;                                       // atom B2
    atom.fpos(1)=x2;atom.fpos(2)=x2;atom.fpos(3)=z2;                     // atom B2
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x2");atom.fpos_equation.push_back("x2");atom.fpos_equation.push_back("z2");// atom B2 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B2 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B2

    atom.name="A"; atom.type=0;                                       // atom B3
    atom.fpos(1)=x3;atom.fpos(2)=x3;atom.fpos(3)=z3;                     // atom B3
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("z3");// atom B3 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B3 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B3

    atom.name="A"; atom.type=0;                                       // atom B16
    atom.fpos(1)=(x16-y16);atom.fpos(2)=(x16+y16);atom.fpos(3)=z16;                     // atom B16
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x16-y16)");atom.fpos_equation.push_back("(x16+y16)");atom.fpos_equation.push_back("z16");// atom B16 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B16 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B16

    atom.name="A"; atom.type=0;                                       // atom B17
    atom.fpos(1)=(x16+y16);atom.fpos(2)=(x16-y16);atom.fpos(3)=z16;                     // atom B17
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x16+y16)");atom.fpos_equation.push_back("(x16-y16)");atom.fpos_equation.push_back("z16");// atom B17 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B17 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B17

    atom.name="A"; atom.type=0;                                       // atom B18
    atom.fpos(1)=(x17-y17);atom.fpos(2)=(x17+y17);atom.fpos(3)=z17;                     // atom B18
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x17-y17)");atom.fpos_equation.push_back("(x17+y17)");atom.fpos_equation.push_back("z17");// atom B18 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B18 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B18

    atom.name="A"; atom.type=0;                                       // atom B19
    atom.fpos(1)=(x17+y17);atom.fpos(2)=(x17-y17);atom.fpos(3)=z17;                     // atom B19
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x17+y17)");atom.fpos_equation.push_back("(x17-y17)");atom.fpos_equation.push_back("z17");// atom B19 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B19 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B19

    atom.name="A"; atom.type=0;                                       // atom B20
    atom.fpos(1)=(x18-y18);atom.fpos(2)=(x18+y18);atom.fpos(3)=z18;                     // atom B20
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x18-y18)");atom.fpos_equation.push_back("(x18+y18)");atom.fpos_equation.push_back("z18");// atom B20 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B20 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B20

    atom.name="A"; atom.type=0;                                       // atom B21
    atom.fpos(1)=(x18+y18);atom.fpos(2)=(x18-y18);atom.fpos(3)=z18;                     // atom B21
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x18+y18)");atom.fpos_equation.push_back("(x18-y18)");atom.fpos_equation.push_back("z18");// atom B21 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B21 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B21

    atom.name="B"; atom.type=1;                                       // atom B4
    atom.fpos(1)=x4;atom.fpos(2)=x4;atom.fpos(3)=z4;                     // atom B4
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x4");atom.fpos_equation.push_back("x4");atom.fpos_equation.push_back("z4");// atom B4 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B4 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B4

    atom.name="B"; atom.type=1;                                       // atom B5
    atom.fpos(1)=x5;atom.fpos(2)=x5;atom.fpos(3)=z5;                     // atom B5
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("z5");// atom B5 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B5 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B5

    atom.name="B"; atom.type=1;                                       // atom B6
    atom.fpos(1)=x6;atom.fpos(2)=x6;atom.fpos(3)=z6;                     // atom B6
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("z6");// atom B6 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B6 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B6

    atom.name="B"; atom.type=1;                                       // atom B7
    atom.fpos(1)=x7;atom.fpos(2)=x7;atom.fpos(3)=z7;                     // atom B7
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x7");atom.fpos_equation.push_back("x7");atom.fpos_equation.push_back("z7");// atom B7 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B7 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B7

    atom.name="B"; atom.type=1;                                       // atom B8
    atom.fpos(1)=x8;atom.fpos(2)=x8;atom.fpos(3)=z8;                     // atom B8
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x8");atom.fpos_equation.push_back("x8");atom.fpos_equation.push_back("z8");// atom B8 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B8 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B8

    atom.name="B"; atom.type=1;                                       // atom B9
    atom.fpos(1)=x9;atom.fpos(2)=x9;atom.fpos(3)=z9;                     // atom B9
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x9");atom.fpos_equation.push_back("x9");atom.fpos_equation.push_back("z9");// atom B9 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B9 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B9

    atom.name="B"; atom.type=1;                                       // atom B10
    atom.fpos(1)=x10;atom.fpos(2)=x10;atom.fpos(3)=z10;                     // atom B10
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x10");atom.fpos_equation.push_back("x10");atom.fpos_equation.push_back("z10");// atom B10 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B10 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B10

    atom.name="B"; atom.type=1;                                       // atom B11
    atom.fpos(1)=x11;atom.fpos(2)=x11;atom.fpos(3)=z11;                     // atom B11
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x11");atom.fpos_equation.push_back("x11");atom.fpos_equation.push_back("z11");// atom B11 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B11 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B11

    atom.name="B"; atom.type=1;                                       // atom B12
    atom.fpos(1)=x12;atom.fpos(2)=x12;atom.fpos(3)=z12;                     // atom B12
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x12");atom.fpos_equation.push_back("x12");atom.fpos_equation.push_back("z12");// atom B12 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B12 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B12

    atom.name="B"; atom.type=1;                                       // atom B22
    atom.fpos(1)=(x19-y19);atom.fpos(2)=(x19+y19);atom.fpos(3)=z19;                     // atom B22
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x19-y19)");atom.fpos_equation.push_back("(x19+y19)");atom.fpos_equation.push_back("z19");// atom B22 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B22 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B22

    atom.name="B"; atom.type=1;                                       // atom B23
    atom.fpos(1)=(x19+y19);atom.fpos(2)=(x19-y19);atom.fpos(3)=z19;                     // atom B23
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x19+y19)");atom.fpos_equation.push_back("(x19-y19)");atom.fpos_equation.push_back("z19");// atom B23 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B23 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B23

    atom.name="B"; atom.type=1;                                       // atom B24
    atom.fpos(1)=(x20-y20);atom.fpos(2)=(x20+y20);atom.fpos(3)=z20;                     // atom B24
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x20-y20)");atom.fpos_equation.push_back("(x20+y20)");atom.fpos_equation.push_back("z20");// atom B24 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B24 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B24

    atom.name="B"; atom.type=1;                                       // atom B25
    atom.fpos(1)=(x20+y20);atom.fpos(2)=(x20-y20);atom.fpos(3)=z20;                     // atom B25
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x20+y20)");atom.fpos_equation.push_back("(x20-y20)");atom.fpos_equation.push_back("z20");// atom B25 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B25 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B25

    atom.name="B"; atom.type=1;                                       // atom B26
    atom.fpos(1)=(x21-y21);atom.fpos(2)=(x21+y21);atom.fpos(3)=z21;                     // atom B26
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x21-y21)");atom.fpos_equation.push_back("(x21+y21)");atom.fpos_equation.push_back("z21");// atom B26 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B26 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B26

    atom.name="B"; atom.type=1;                                       // atom B27
    atom.fpos(1)=(x21+y21);atom.fpos(2)=(x21-y21);atom.fpos(3)=z21;                     // atom B27
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x21+y21)");atom.fpos_equation.push_back("(x21-y21)");atom.fpos_equation.push_back("z21");// atom B27 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B27 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B27

    atom.name="C"; atom.type=2;                                       // atom B13
    atom.fpos(1)=x13;atom.fpos(2)=x13;atom.fpos(3)=z13;                     // atom B13
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x13");atom.fpos_equation.push_back("x13");atom.fpos_equation.push_back("z13");// atom B13 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B13 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B13

    atom.name="C"; atom.type=2;                                       // atom B14
    atom.fpos(1)=x14;atom.fpos(2)=x14;atom.fpos(3)=z14;                     // atom B14
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x14");atom.fpos_equation.push_back("x14");atom.fpos_equation.push_back("z14");// atom B14 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B14 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B14

    atom.name="C"; atom.type=2;                                       // atom B15
    atom.fpos(1)=x15;atom.fpos(2)=x15;atom.fpos(3)=z15;                     // atom B15
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x15");atom.fpos_equation.push_back("x15");atom.fpos_equation.push_back("z15");// atom B15 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B15 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B15


    return str.atoms.size();  
  }
} // namespace anrl

namespace anrl {
  uint WebANRL_A3B5C_mC54_8_3a3b_9a3b_3a(stringstream& web,bool LDEBUG) {
#ifndef _ANRL_NOWEB_
#endif

    if(LDEBUG) {cerr << "anrl:: WebANRL_A3B5C_mC54_8_3a3b_9a3b_3a: web.str().size()=" << web.str().size() << endl;}

    return web.str().size();
  }
} // namespace anrl

#endif // AFLOW_REMOVE_GREP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
