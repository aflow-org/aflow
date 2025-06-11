// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - David Hicks - 2020
// FILE "ANRL/aflow_anrl_AB3C_aP30_2_3i_9i_3i.cpp"

#ifndef _AFLOW_ANRL_AB3C_aP30_2_3i_9i_3i_CPP // AFLOW_REMOVE_GREP
#define _AFLOW_ANRL_AB3C_aP30_2_3i_9i_3i_CPP // AFLOW_REMOVE_GREP
#include "../aflow.h" // AFLOW_REMOVE_GREP

namespace anrl {
  uint WebANRL_AB3C_aP30_2_3i_9i_3i(stringstream &web,bool LDEBUG);
}

namespace anrl {
  uint PrototypeANRL_AB3C_aP30_2_3i_9i_3i(stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG) {
    // system AB3C_aP30_2_3i_9i_3i

    if(XHOST.vflag_control.flag("WWW")) {
      WebANRL_AB3C_aP30_2_3i_9i_3i(web,LDEBUG); // PLUG WEB STUFF
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

    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: FOUND" << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: label=" << label << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: nspecies=" << nspecies << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: natoms=" << natoms << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: spacegroup=" << spacegroup << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: nunderscores=" << nunderscores << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: nparameters=" <<  nparameters << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: Pearson_symbol=" << Pearson_symbol << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: params=" << params << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: Strukturbericht=" << Strukturbericht << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: prototype=" << prototype << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: dialect=" << dialect << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: vparameters.size()=" << vparameters.size() << endl;}

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
    double a=vparameters.at(i++);                  if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: a=" << a << endl;}
    double bovera=vparameters.at(i++),b=bovera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: b=" << b << " (b/a=" << bovera << ")" << endl;}
    double covera=vparameters.at(i++),c=covera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: c=" << c << " (c/a=" << covera << ")" << endl;}
    double alpha=vparameters.at(i++);              if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: alpha=" << alpha << endl;}
    double beta=vparameters.at(i++);               if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: beta=" << beta << endl;}
    double gamma=vparameters.at(i++);              if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: gamma=" << gamma << endl;}

    double x1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: x1=" << x1 << endl;}
    double y1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: y1=" << y1 << endl;}
    double z1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: z1=" << z1 << endl;}
    double x2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: x2=" << x2 << endl;}
    double y2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: y2=" << y2 << endl;}
    double z2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: z2=" << z2 << endl;}
    double x3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: x3=" << x3 << endl;}
    double y3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: y3=" << y3 << endl;}
    double z3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: z3=" << z3 << endl;}
    double x4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: x4=" << x4 << endl;}
    double y4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: y4=" << y4 << endl;}
    double z4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: z4=" << z4 << endl;}
    double x5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: x5=" << x5 << endl;}
    double y5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: y5=" << y5 << endl;}
    double z5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: z5=" << z5 << endl;}
    double x6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: x6=" << x6 << endl;}
    double y6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: y6=" << y6 << endl;}
    double z6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: z6=" << z6 << endl;}
    double x7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: x7=" << x7 << endl;}
    double y7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: y7=" << y7 << endl;}
    double z7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: z7=" << z7 << endl;}
    double x8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: x8=" << x8 << endl;}
    double y8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: y8=" << y8 << endl;}
    double z8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: z8=" << z8 << endl;}
    double x9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: x9=" << x9 << endl;}
    double y9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: y9=" << y9 << endl;}
    double z9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: z9=" << z9 << endl;}
    double x10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: x10=" << x10 << endl;}
    double y10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: y10=" << y10 << endl;}
    double z10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: z10=" << z10 << endl;}
    double x11=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: x11=" << x11 << endl;}
    double y11=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: y11=" << y11 << endl;}
    double z11=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: z11=" << z11 << endl;}
    double x12=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: x12=" << x12 << endl;}
    double y12=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: y12=" << y12 << endl;}
    double z12=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: z12=" << z12 << endl;}
    double x13=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: x13=" << x13 << endl;}
    double y13=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: y13=" << y13 << endl;}
    double z13=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: z13=" << z13 << endl;}
    double x14=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: x14=" << x14 << endl;}
    double y14=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: y14=" << y14 << endl;}
    double z14=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: z14=" << z14 << endl;}
    double x15=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: x15=" << x15 << endl;}
    double y15=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: y15=" << y15 << endl;}
    double z15=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: z15=" << z15 << endl;}

    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: cos(alpha)=" << cos(deg2rad*alpha)  << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: sin(alpha)=" << sin(deg2rad*alpha)  << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: cos(beta)=" << cos(deg2rad*beta)  << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: sin(beta)=" << sin(deg2rad*beta)  << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: cos(gamma)=" << cos(deg2rad*gamma)  << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB3C_aP30_2_3i_9i_3i: sin(gamma)=" << sin(deg2rad*gamma)  << endl;}

    str.iomode=IOVASP_AUTO;
    str.title=label+" params="+parameters+" SG="+aurostd::utype2string(spacegroup)+DOI_ANRL;
    str.scale=1.0;

    double cx=c*cos(deg2rad*beta);
    double cy=c*(cos(deg2rad*alpha)-cos(deg2rad*beta)*cos(deg2rad*gamma))/sin(deg2rad*gamma);
    double cz=sqrt(pow(c,2.0)-pow(cx,2.0)-pow(cy,2.0));

    a1=a*xn;
    a2=b*cos(deg2rad*gamma)*xn+b*sin(deg2rad*gamma)*yn;
    a3=cx*xn+cy*yn+cz*zn;

    str.lattice(1,1)=a1(1);str.lattice(1,2)=a1(2);str.lattice(1,3)=a1(3);
    str.lattice(2,1)=a2(1);str.lattice(2,2)=a2(2);str.lattice(2,3)=a2(3);
    str.lattice(3,1)=a3(1);str.lattice(3,2)=a3(2);str.lattice(3,3)=a3(3);

    // symbolic representation of lattice vectors
    vector<string> a1_equation, a2_equation, a3_equation;
    a1_equation.push_back("a");a1_equation.push_back("0");a1_equation.push_back("0");
    a2_equation.push_back("b*cos(gamma)");a2_equation.push_back("b*sin(gamma)");a2_equation.push_back("0");
    a3_equation.push_back("cx");a3_equation.push_back("cy");a3_equation.push_back("cz");
    str.symbolic_math_lattice.push_back(a1_equation);
    str.symbolic_math_lattice.push_back(a2_equation);
    str.symbolic_math_lattice.push_back(a3_equation);

    str.num_lattice_parameters = 6;

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
    atom.fpos(1)=-x1;atom.fpos(2)=-y1;atom.fpos(3)=-z1;                     // atom B2
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x1");atom.fpos_equation.push_back("-y1");atom.fpos_equation.push_back("-z1");// atom B2 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B2 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B2

    atom.name="A"; atom.type=0;                                       // atom B3
    atom.fpos(1)=x2;atom.fpos(2)=y2;atom.fpos(3)=z2;                     // atom B3
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x2");atom.fpos_equation.push_back("y2");atom.fpos_equation.push_back("z2");// atom B3 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B3 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B3

    atom.name="A"; atom.type=0;                                       // atom B4
    atom.fpos(1)=-x2;atom.fpos(2)=-y2;atom.fpos(3)=-z2;                     // atom B4
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x2");atom.fpos_equation.push_back("-y2");atom.fpos_equation.push_back("-z2");// atom B4 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B4 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B4

    atom.name="A"; atom.type=0;                                       // atom B5
    atom.fpos(1)=x3;atom.fpos(2)=y3;atom.fpos(3)=z3;                     // atom B5
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("y3");atom.fpos_equation.push_back("z3");// atom B5 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B5 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B5

    atom.name="A"; atom.type=0;                                       // atom B6
    atom.fpos(1)=-x3;atom.fpos(2)=-y3;atom.fpos(3)=-z3;                     // atom B6
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("-y3");atom.fpos_equation.push_back("-z3");// atom B6 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B6 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B6

    atom.name="B"; atom.type=1;                                       // atom B7
    atom.fpos(1)=x4;atom.fpos(2)=y4;atom.fpos(3)=z4;                     // atom B7
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x4");atom.fpos_equation.push_back("y4");atom.fpos_equation.push_back("z4");// atom B7 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B7 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B7

    atom.name="B"; atom.type=1;                                       // atom B8
    atom.fpos(1)=-x4;atom.fpos(2)=-y4;atom.fpos(3)=-z4;                     // atom B8
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x4");atom.fpos_equation.push_back("-y4");atom.fpos_equation.push_back("-z4");// atom B8 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B8 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B8

    atom.name="B"; atom.type=1;                                       // atom B9
    atom.fpos(1)=x5;atom.fpos(2)=y5;atom.fpos(3)=z5;                     // atom B9
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("y5");atom.fpos_equation.push_back("z5");// atom B9 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B9 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B9

    atom.name="B"; atom.type=1;                                       // atom B10
    atom.fpos(1)=-x5;atom.fpos(2)=-y5;atom.fpos(3)=-z5;                     // atom B10
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("-y5");atom.fpos_equation.push_back("-z5");// atom B10 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B10 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B10

    atom.name="B"; atom.type=1;                                       // atom B11
    atom.fpos(1)=x6;atom.fpos(2)=y6;atom.fpos(3)=z6;                     // atom B11
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("y6");atom.fpos_equation.push_back("z6");// atom B11 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B11 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B11

    atom.name="B"; atom.type=1;                                       // atom B12
    atom.fpos(1)=-x6;atom.fpos(2)=-y6;atom.fpos(3)=-z6;                     // atom B12
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("-y6");atom.fpos_equation.push_back("-z6");// atom B12 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B12 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B12

    atom.name="B"; atom.type=1;                                       // atom B13
    atom.fpos(1)=x7;atom.fpos(2)=y7;atom.fpos(3)=z7;                     // atom B13
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x7");atom.fpos_equation.push_back("y7");atom.fpos_equation.push_back("z7");// atom B13 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B13 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B13

    atom.name="B"; atom.type=1;                                       // atom B14
    atom.fpos(1)=-x7;atom.fpos(2)=-y7;atom.fpos(3)=-z7;                     // atom B14
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x7");atom.fpos_equation.push_back("-y7");atom.fpos_equation.push_back("-z7");// atom B14 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B14 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B14

    atom.name="B"; atom.type=1;                                       // atom B15
    atom.fpos(1)=x8;atom.fpos(2)=y8;atom.fpos(3)=z8;                     // atom B15
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x8");atom.fpos_equation.push_back("y8");atom.fpos_equation.push_back("z8");// atom B15 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B15 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B15

    atom.name="B"; atom.type=1;                                       // atom B16
    atom.fpos(1)=-x8;atom.fpos(2)=-y8;atom.fpos(3)=-z8;                     // atom B16
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x8");atom.fpos_equation.push_back("-y8");atom.fpos_equation.push_back("-z8");// atom B16 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B16 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B16

    atom.name="B"; atom.type=1;                                       // atom B17
    atom.fpos(1)=x9;atom.fpos(2)=y9;atom.fpos(3)=z9;                     // atom B17
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x9");atom.fpos_equation.push_back("y9");atom.fpos_equation.push_back("z9");// atom B17 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B17 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B17

    atom.name="B"; atom.type=1;                                       // atom B18
    atom.fpos(1)=-x9;atom.fpos(2)=-y9;atom.fpos(3)=-z9;                     // atom B18
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x9");atom.fpos_equation.push_back("-y9");atom.fpos_equation.push_back("-z9");// atom B18 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B18 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B18

    atom.name="B"; atom.type=1;                                       // atom B19
    atom.fpos(1)=x10;atom.fpos(2)=y10;atom.fpos(3)=z10;                     // atom B19
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x10");atom.fpos_equation.push_back("y10");atom.fpos_equation.push_back("z10");// atom B19 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B19 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B19

    atom.name="B"; atom.type=1;                                       // atom B20
    atom.fpos(1)=-x10;atom.fpos(2)=-y10;atom.fpos(3)=-z10;                     // atom B20
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x10");atom.fpos_equation.push_back("-y10");atom.fpos_equation.push_back("-z10");// atom B20 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B20 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B20

    atom.name="B"; atom.type=1;                                       // atom B21
    atom.fpos(1)=x11;atom.fpos(2)=y11;atom.fpos(3)=z11;                     // atom B21
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x11");atom.fpos_equation.push_back("y11");atom.fpos_equation.push_back("z11");// atom B21 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B21 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B21

    atom.name="B"; atom.type=1;                                       // atom B22
    atom.fpos(1)=-x11;atom.fpos(2)=-y11;atom.fpos(3)=-z11;                     // atom B22
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x11");atom.fpos_equation.push_back("-y11");atom.fpos_equation.push_back("-z11");// atom B22 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B22 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B22

    atom.name="B"; atom.type=1;                                       // atom B23
    atom.fpos(1)=x12;atom.fpos(2)=y12;atom.fpos(3)=z12;                     // atom B23
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x12");atom.fpos_equation.push_back("y12");atom.fpos_equation.push_back("z12");// atom B23 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B23 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B23

    atom.name="B"; atom.type=1;                                       // atom B24
    atom.fpos(1)=-x12;atom.fpos(2)=-y12;atom.fpos(3)=-z12;                     // atom B24
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x12");atom.fpos_equation.push_back("-y12");atom.fpos_equation.push_back("-z12");// atom B24 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B24 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B24

    atom.name="C"; atom.type=2;                                       // atom B25
    atom.fpos(1)=x13;atom.fpos(2)=y13;atom.fpos(3)=z13;                     // atom B25
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x13");atom.fpos_equation.push_back("y13");atom.fpos_equation.push_back("z13");// atom B25 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B25 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B25

    atom.name="C"; atom.type=2;                                       // atom B26
    atom.fpos(1)=-x13;atom.fpos(2)=-y13;atom.fpos(3)=-z13;                     // atom B26
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x13");atom.fpos_equation.push_back("-y13");atom.fpos_equation.push_back("-z13");// atom B26 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B26 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B26

    atom.name="C"; atom.type=2;                                       // atom B27
    atom.fpos(1)=x14;atom.fpos(2)=y14;atom.fpos(3)=z14;                     // atom B27
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x14");atom.fpos_equation.push_back("y14");atom.fpos_equation.push_back("z14");// atom B27 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B27 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B27

    atom.name="C"; atom.type=2;                                       // atom B28
    atom.fpos(1)=-x14;atom.fpos(2)=-y14;atom.fpos(3)=-z14;                     // atom B28
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x14");atom.fpos_equation.push_back("-y14");atom.fpos_equation.push_back("-z14");// atom B28 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B28 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B28

    atom.name="C"; atom.type=2;                                       // atom B29
    atom.fpos(1)=x15;atom.fpos(2)=y15;atom.fpos(3)=z15;                     // atom B29
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x15");atom.fpos_equation.push_back("y15");atom.fpos_equation.push_back("z15");// atom B29 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B29 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B29

    atom.name="C"; atom.type=2;                                       // atom B30
    atom.fpos(1)=-x15;atom.fpos(2)=-y15;atom.fpos(3)=-z15;                     // atom B30
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x15");atom.fpos_equation.push_back("-y15");atom.fpos_equation.push_back("-z15");// atom B30 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B30 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B30


    return str.atoms.size();  
  }
} // namespace anrl

namespace anrl {
  uint WebANRL_AB3C_aP30_2_3i_9i_3i(stringstream& web,bool LDEBUG) {
#ifndef _ANRL_NOWEB_
#endif

    if(LDEBUG) {cerr << "anrl:: WebANRL_AB3C_aP30_2_3i_9i_3i: web.str().size()=" << web.str().size() << endl;}

    return web.str().size();
  }
} // namespace anrl

#endif // AFLOW_REMOVE_GREP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
