// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - David Hicks - 2018
// FILE "ANRL/aflow_anrl_A6B4C16D_oP108_27_abcd4e_4e_16e_e.cpp"

#ifndef _AFLOW_ANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e_CPP // AFLOW_REMOVE_GREP
#define _AFLOW_ANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e_CPP // AFLOW_REMOVE_GREP
#include "../aflow.h" // AFLOW_REMOVE_GREP

namespace anrl {
  uint WebANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e(stringstream &web,bool LDEBUG);
}

namespace anrl {
  uint PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e(stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG) {
    // system A6B4C16D_oP108_27_abcd4e_4e_16e_e

    if(XHOST.vflag_control.flag("WWW")) {
      WebANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e(web,LDEBUG); // PLUG WEB STUFF
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

    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: FOUND" << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: label=" << label << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: nspecies=" << nspecies << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: natoms=" << natoms << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: spacegroup=" << spacegroup << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: nunderscores=" << nunderscores << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: nparameters=" <<  nparameters << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: Pearson_symbol=" << Pearson_symbol << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: params=" << params << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: Strukturbericht=" << Strukturbericht << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: prototype=" << prototype << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: dialect=" << dialect << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: vparameters.size()=" << vparameters.size() << endl;}

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
    double a=vparameters.at(i++);                  if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: a=" << a << endl;}
    double bovera=vparameters.at(i++),b=bovera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: b=" << b << " (b/a=" << bovera << ")" << endl;}
    double covera=vparameters.at(i++),c=covera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: c=" << c << " (c/a=" << covera << ")" << endl;}

    double z1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z1=" << z1 << endl;}
    double z2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z2=" << z2 << endl;}
    double z3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z3=" << z3 << endl;}
    double z4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z4=" << z4 << endl;}
    double x5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x5=" << x5 << endl;}
    double y5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y5=" << y5 << endl;}
    double z5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z5=" << z5 << endl;}
    double x6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x6=" << x6 << endl;}
    double y6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y6=" << y6 << endl;}
    double z6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z6=" << z6 << endl;}
    double x7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x7=" << x7 << endl;}
    double y7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y7=" << y7 << endl;}
    double z7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z7=" << z7 << endl;}
    double x8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x8=" << x8 << endl;}
    double y8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y8=" << y8 << endl;}
    double z8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z8=" << z8 << endl;}
    double x9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x9=" << x9 << endl;}
    double y9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y9=" << y9 << endl;}
    double z9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z9=" << z9 << endl;}
    double x10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x10=" << x10 << endl;}
    double y10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y10=" << y10 << endl;}
    double z10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z10=" << z10 << endl;}
    double x11=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x11=" << x11 << endl;}
    double y11=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y11=" << y11 << endl;}
    double z11=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z11=" << z11 << endl;}
    double x12=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x12=" << x12 << endl;}
    double y12=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y12=" << y12 << endl;}
    double z12=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z12=" << z12 << endl;}
    double x13=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x13=" << x13 << endl;}
    double y13=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y13=" << y13 << endl;}
    double z13=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z13=" << z13 << endl;}
    double x14=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x14=" << x14 << endl;}
    double y14=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y14=" << y14 << endl;}
    double z14=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z14=" << z14 << endl;}
    double x15=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x15=" << x15 << endl;}
    double y15=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y15=" << y15 << endl;}
    double z15=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z15=" << z15 << endl;}
    double x16=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x16=" << x16 << endl;}
    double y16=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y16=" << y16 << endl;}
    double z16=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z16=" << z16 << endl;}
    double x17=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x17=" << x17 << endl;}
    double y17=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y17=" << y17 << endl;}
    double z17=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z17=" << z17 << endl;}
    double x18=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x18=" << x18 << endl;}
    double y18=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y18=" << y18 << endl;}
    double z18=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z18=" << z18 << endl;}
    double x19=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x19=" << x19 << endl;}
    double y19=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y19=" << y19 << endl;}
    double z19=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z19=" << z19 << endl;}
    double x20=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x20=" << x20 << endl;}
    double y20=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y20=" << y20 << endl;}
    double z20=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z20=" << z20 << endl;}
    double x21=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x21=" << x21 << endl;}
    double y21=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y21=" << y21 << endl;}
    double z21=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z21=" << z21 << endl;}
    double x22=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x22=" << x22 << endl;}
    double y22=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y22=" << y22 << endl;}
    double z22=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z22=" << z22 << endl;}
    double x23=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x23=" << x23 << endl;}
    double y23=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y23=" << y23 << endl;}
    double z23=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z23=" << z23 << endl;}
    double x24=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x24=" << x24 << endl;}
    double y24=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y24=" << y24 << endl;}
    double z24=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z24=" << z24 << endl;}
    double x25=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x25=" << x25 << endl;}
    double y25=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y25=" << y25 << endl;}
    double z25=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z25=" << z25 << endl;}
    double x26=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x26=" << x26 << endl;}
    double y26=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y26=" << y26 << endl;}
    double z26=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z26=" << z26 << endl;}
    double x27=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x27=" << x27 << endl;}
    double y27=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y27=" << y27 << endl;}
    double z27=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z27=" << z27 << endl;}
    double x28=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x28=" << x28 << endl;}
    double y28=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y28=" << y28 << endl;}
    double z28=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z28=" << z28 << endl;}
    double x29=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: x29=" << x29 << endl;}
    double y29=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: y29=" << y29 << endl;}
    double z29=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: z29=" << z29 << endl;}

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
    atom.fpos(1)=0.0;atom.fpos(2)=0.0;atom.fpos(3)=((1.0/2.0)+z1);                     // atom B2
    atom.fpos_equation.clear();atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("((1.0/2.0)+z1)");// atom B2 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B2 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B2

    atom.name="A"; atom.type=0;                                       // atom B3
    atom.fpos(1)=0.0;atom.fpos(2)=(1.0/2.0);atom.fpos(3)=z2;                     // atom B3
    atom.fpos_equation.clear();atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("z2");// atom B3 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B3 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B3

    atom.name="A"; atom.type=0;                                       // atom B4
    atom.fpos(1)=0.0;atom.fpos(2)=(1.0/2.0);atom.fpos(3)=((1.0/2.0)+z2);                     // atom B4
    atom.fpos_equation.clear();atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("((1.0/2.0)+z2)");// atom B4 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B4 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B4

    atom.name="A"; atom.type=0;                                       // atom B5
    atom.fpos(1)=(1.0/2.0);atom.fpos(2)=0.0;atom.fpos(3)=z3;                     // atom B5
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("z3");// atom B5 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B5 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B5

    atom.name="A"; atom.type=0;                                       // atom B6
    atom.fpos(1)=(1.0/2.0);atom.fpos(2)=0.0;atom.fpos(3)=((1.0/2.0)+z3);                     // atom B6
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("((1.0/2.0)+z3)");// atom B6 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B6 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B6

    atom.name="A"; atom.type=0;                                       // atom B7
    atom.fpos(1)=(1.0/2.0);atom.fpos(2)=(1.0/2.0);atom.fpos(3)=z4;                     // atom B7
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("z4");// atom B7 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B7 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B7

    atom.name="A"; atom.type=0;                                       // atom B8
    atom.fpos(1)=(1.0/2.0);atom.fpos(2)=(1.0/2.0);atom.fpos(3)=((1.0/2.0)+z4);                     // atom B8
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("((1.0/2.0)+z4)");// atom B8 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B8 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B8

    atom.name="A"; atom.type=0;                                       // atom B9
    atom.fpos(1)=x5;atom.fpos(2)=y5;atom.fpos(3)=z5;                     // atom B9
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("y5");atom.fpos_equation.push_back("z5");// atom B9 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B9 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B9

    atom.name="A"; atom.type=0;                                       // atom B10
    atom.fpos(1)=-x5;atom.fpos(2)=-y5;atom.fpos(3)=z5;                     // atom B10
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("-y5");atom.fpos_equation.push_back("z5");// atom B10 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B10 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B10

    atom.name="A"; atom.type=0;                                       // atom B11
    atom.fpos(1)=x5;atom.fpos(2)=-y5;atom.fpos(3)=((1.0/2.0)+z5);                     // atom B11
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("-y5");atom.fpos_equation.push_back("((1.0/2.0)+z5)");// atom B11 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B11 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B11

    atom.name="A"; atom.type=0;                                       // atom B12
    atom.fpos(1)=-x5;atom.fpos(2)=y5;atom.fpos(3)=((1.0/2.0)+z5);                     // atom B12
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("y5");atom.fpos_equation.push_back("((1.0/2.0)+z5)");// atom B12 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B12 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B12

    atom.name="A"; atom.type=0;                                       // atom B13
    atom.fpos(1)=x6;atom.fpos(2)=y6;atom.fpos(3)=z6;                     // atom B13
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("y6");atom.fpos_equation.push_back("z6");// atom B13 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B13 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B13

    atom.name="A"; atom.type=0;                                       // atom B14
    atom.fpos(1)=-x6;atom.fpos(2)=-y6;atom.fpos(3)=z6;                     // atom B14
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("-y6");atom.fpos_equation.push_back("z6");// atom B14 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B14 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B14

    atom.name="A"; atom.type=0;                                       // atom B15
    atom.fpos(1)=x6;atom.fpos(2)=-y6;atom.fpos(3)=((1.0/2.0)+z6);                     // atom B15
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("-y6");atom.fpos_equation.push_back("((1.0/2.0)+z6)");// atom B15 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B15 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B15

    atom.name="A"; atom.type=0;                                       // atom B16
    atom.fpos(1)=-x6;atom.fpos(2)=y6;atom.fpos(3)=((1.0/2.0)+z6);                     // atom B16
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("y6");atom.fpos_equation.push_back("((1.0/2.0)+z6)");// atom B16 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B16 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B16

    atom.name="A"; atom.type=0;                                       // atom B17
    atom.fpos(1)=x7;atom.fpos(2)=y7;atom.fpos(3)=z7;                     // atom B17
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x7");atom.fpos_equation.push_back("y7");atom.fpos_equation.push_back("z7");// atom B17 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B17 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B17

    atom.name="A"; atom.type=0;                                       // atom B18
    atom.fpos(1)=-x7;atom.fpos(2)=-y7;atom.fpos(3)=z7;                     // atom B18
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x7");atom.fpos_equation.push_back("-y7");atom.fpos_equation.push_back("z7");// atom B18 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B18 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B18

    atom.name="A"; atom.type=0;                                       // atom B19
    atom.fpos(1)=x7;atom.fpos(2)=-y7;atom.fpos(3)=((1.0/2.0)+z7);                     // atom B19
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x7");atom.fpos_equation.push_back("-y7");atom.fpos_equation.push_back("((1.0/2.0)+z7)");// atom B19 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B19 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B19

    atom.name="A"; atom.type=0;                                       // atom B20
    atom.fpos(1)=-x7;atom.fpos(2)=y7;atom.fpos(3)=((1.0/2.0)+z7);                     // atom B20
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x7");atom.fpos_equation.push_back("y7");atom.fpos_equation.push_back("((1.0/2.0)+z7)");// atom B20 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B20 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B20

    atom.name="A"; atom.type=0;                                       // atom B21
    atom.fpos(1)=x8;atom.fpos(2)=y8;atom.fpos(3)=z8;                     // atom B21
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x8");atom.fpos_equation.push_back("y8");atom.fpos_equation.push_back("z8");// atom B21 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B21 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B21

    atom.name="A"; atom.type=0;                                       // atom B22
    atom.fpos(1)=-x8;atom.fpos(2)=-y8;atom.fpos(3)=z8;                     // atom B22
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x8");atom.fpos_equation.push_back("-y8");atom.fpos_equation.push_back("z8");// atom B22 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B22 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B22

    atom.name="A"; atom.type=0;                                       // atom B23
    atom.fpos(1)=x8;atom.fpos(2)=-y8;atom.fpos(3)=((1.0/2.0)+z8);                     // atom B23
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x8");atom.fpos_equation.push_back("-y8");atom.fpos_equation.push_back("((1.0/2.0)+z8)");// atom B23 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B23 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B23

    atom.name="A"; atom.type=0;                                       // atom B24
    atom.fpos(1)=-x8;atom.fpos(2)=y8;atom.fpos(3)=((1.0/2.0)+z8);                     // atom B24
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x8");atom.fpos_equation.push_back("y8");atom.fpos_equation.push_back("((1.0/2.0)+z8)");// atom B24 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B24 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B24

    atom.name="B"; atom.type=1;                                       // atom B25
    atom.fpos(1)=x9;atom.fpos(2)=y9;atom.fpos(3)=z9;                     // atom B25
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x9");atom.fpos_equation.push_back("y9");atom.fpos_equation.push_back("z9");// atom B25 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B25 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B25

    atom.name="B"; atom.type=1;                                       // atom B26
    atom.fpos(1)=-x9;atom.fpos(2)=-y9;atom.fpos(3)=z9;                     // atom B26
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x9");atom.fpos_equation.push_back("-y9");atom.fpos_equation.push_back("z9");// atom B26 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B26 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B26

    atom.name="B"; atom.type=1;                                       // atom B27
    atom.fpos(1)=x9;atom.fpos(2)=-y9;atom.fpos(3)=((1.0/2.0)+z9);                     // atom B27
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x9");atom.fpos_equation.push_back("-y9");atom.fpos_equation.push_back("((1.0/2.0)+z9)");// atom B27 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B27 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B27

    atom.name="B"; atom.type=1;                                       // atom B28
    atom.fpos(1)=-x9;atom.fpos(2)=y9;atom.fpos(3)=((1.0/2.0)+z9);                     // atom B28
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x9");atom.fpos_equation.push_back("y9");atom.fpos_equation.push_back("((1.0/2.0)+z9)");// atom B28 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B28 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B28

    atom.name="B"; atom.type=1;                                       // atom B29
    atom.fpos(1)=x10;atom.fpos(2)=y10;atom.fpos(3)=z10;                     // atom B29
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x10");atom.fpos_equation.push_back("y10");atom.fpos_equation.push_back("z10");// atom B29 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B29 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B29

    atom.name="B"; atom.type=1;                                       // atom B30
    atom.fpos(1)=-x10;atom.fpos(2)=-y10;atom.fpos(3)=z10;                     // atom B30
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x10");atom.fpos_equation.push_back("-y10");atom.fpos_equation.push_back("z10");// atom B30 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B30 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B30

    atom.name="B"; atom.type=1;                                       // atom B31
    atom.fpos(1)=x10;atom.fpos(2)=-y10;atom.fpos(3)=((1.0/2.0)+z10);                     // atom B31
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x10");atom.fpos_equation.push_back("-y10");atom.fpos_equation.push_back("((1.0/2.0)+z10)");// atom B31 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B31 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B31

    atom.name="B"; atom.type=1;                                       // atom B32
    atom.fpos(1)=-x10;atom.fpos(2)=y10;atom.fpos(3)=((1.0/2.0)+z10);                     // atom B32
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x10");atom.fpos_equation.push_back("y10");atom.fpos_equation.push_back("((1.0/2.0)+z10)");// atom B32 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B32 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B32

    atom.name="B"; atom.type=1;                                       // atom B33
    atom.fpos(1)=x11;atom.fpos(2)=y11;atom.fpos(3)=z11;                     // atom B33
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x11");atom.fpos_equation.push_back("y11");atom.fpos_equation.push_back("z11");// atom B33 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B33 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B33

    atom.name="B"; atom.type=1;                                       // atom B34
    atom.fpos(1)=-x11;atom.fpos(2)=-y11;atom.fpos(3)=z11;                     // atom B34
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x11");atom.fpos_equation.push_back("-y11");atom.fpos_equation.push_back("z11");// atom B34 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B34 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B34

    atom.name="B"; atom.type=1;                                       // atom B35
    atom.fpos(1)=x11;atom.fpos(2)=-y11;atom.fpos(3)=((1.0/2.0)+z11);                     // atom B35
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x11");atom.fpos_equation.push_back("-y11");atom.fpos_equation.push_back("((1.0/2.0)+z11)");// atom B35 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B35 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B35

    atom.name="B"; atom.type=1;                                       // atom B36
    atom.fpos(1)=-x11;atom.fpos(2)=y11;atom.fpos(3)=((1.0/2.0)+z11);                     // atom B36
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x11");atom.fpos_equation.push_back("y11");atom.fpos_equation.push_back("((1.0/2.0)+z11)");// atom B36 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B36 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B36

    atom.name="B"; atom.type=1;                                       // atom B37
    atom.fpos(1)=x12;atom.fpos(2)=y12;atom.fpos(3)=z12;                     // atom B37
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x12");atom.fpos_equation.push_back("y12");atom.fpos_equation.push_back("z12");// atom B37 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B37 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B37

    atom.name="B"; atom.type=1;                                       // atom B38
    atom.fpos(1)=-x12;atom.fpos(2)=-y12;atom.fpos(3)=z12;                     // atom B38
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x12");atom.fpos_equation.push_back("-y12");atom.fpos_equation.push_back("z12");// atom B38 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B38 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B38

    atom.name="B"; atom.type=1;                                       // atom B39
    atom.fpos(1)=x12;atom.fpos(2)=-y12;atom.fpos(3)=((1.0/2.0)+z12);                     // atom B39
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x12");atom.fpos_equation.push_back("-y12");atom.fpos_equation.push_back("((1.0/2.0)+z12)");// atom B39 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B39 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B39

    atom.name="B"; atom.type=1;                                       // atom B40
    atom.fpos(1)=-x12;atom.fpos(2)=y12;atom.fpos(3)=((1.0/2.0)+z12);                     // atom B40
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x12");atom.fpos_equation.push_back("y12");atom.fpos_equation.push_back("((1.0/2.0)+z12)");// atom B40 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B40 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B40

    atom.name="C"; atom.type=2;                                       // atom B41
    atom.fpos(1)=x13;atom.fpos(2)=y13;atom.fpos(3)=z13;                     // atom B41
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x13");atom.fpos_equation.push_back("y13");atom.fpos_equation.push_back("z13");// atom B41 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B41 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B41

    atom.name="C"; atom.type=2;                                       // atom B42
    atom.fpos(1)=-x13;atom.fpos(2)=-y13;atom.fpos(3)=z13;                     // atom B42
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x13");atom.fpos_equation.push_back("-y13");atom.fpos_equation.push_back("z13");// atom B42 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B42 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B42

    atom.name="C"; atom.type=2;                                       // atom B43
    atom.fpos(1)=x13;atom.fpos(2)=-y13;atom.fpos(3)=((1.0/2.0)+z13);                     // atom B43
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x13");atom.fpos_equation.push_back("-y13");atom.fpos_equation.push_back("((1.0/2.0)+z13)");// atom B43 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B43 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B43

    atom.name="C"; atom.type=2;                                       // atom B44
    atom.fpos(1)=-x13;atom.fpos(2)=y13;atom.fpos(3)=((1.0/2.0)+z13);                     // atom B44
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x13");atom.fpos_equation.push_back("y13");atom.fpos_equation.push_back("((1.0/2.0)+z13)");// atom B44 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B44 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B44

    atom.name="C"; atom.type=2;                                       // atom B45
    atom.fpos(1)=x14;atom.fpos(2)=y14;atom.fpos(3)=z14;                     // atom B45
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x14");atom.fpos_equation.push_back("y14");atom.fpos_equation.push_back("z14");// atom B45 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B45 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B45

    atom.name="C"; atom.type=2;                                       // atom B46
    atom.fpos(1)=-x14;atom.fpos(2)=-y14;atom.fpos(3)=z14;                     // atom B46
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x14");atom.fpos_equation.push_back("-y14");atom.fpos_equation.push_back("z14");// atom B46 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B46 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B46

    atom.name="C"; atom.type=2;                                       // atom B47
    atom.fpos(1)=x14;atom.fpos(2)=-y14;atom.fpos(3)=((1.0/2.0)+z14);                     // atom B47
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x14");atom.fpos_equation.push_back("-y14");atom.fpos_equation.push_back("((1.0/2.0)+z14)");// atom B47 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B47 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B47

    atom.name="C"; atom.type=2;                                       // atom B48
    atom.fpos(1)=-x14;atom.fpos(2)=y14;atom.fpos(3)=((1.0/2.0)+z14);                     // atom B48
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x14");atom.fpos_equation.push_back("y14");atom.fpos_equation.push_back("((1.0/2.0)+z14)");// atom B48 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B48 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B48

    atom.name="C"; atom.type=2;                                       // atom B49
    atom.fpos(1)=x15;atom.fpos(2)=y15;atom.fpos(3)=z15;                     // atom B49
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x15");atom.fpos_equation.push_back("y15");atom.fpos_equation.push_back("z15");// atom B49 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B49 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B49

    atom.name="C"; atom.type=2;                                       // atom B50
    atom.fpos(1)=-x15;atom.fpos(2)=-y15;atom.fpos(3)=z15;                     // atom B50
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x15");atom.fpos_equation.push_back("-y15");atom.fpos_equation.push_back("z15");// atom B50 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B50 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B50

    atom.name="C"; atom.type=2;                                       // atom B51
    atom.fpos(1)=x15;atom.fpos(2)=-y15;atom.fpos(3)=((1.0/2.0)+z15);                     // atom B51
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x15");atom.fpos_equation.push_back("-y15");atom.fpos_equation.push_back("((1.0/2.0)+z15)");// atom B51 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B51 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B51

    atom.name="C"; atom.type=2;                                       // atom B52
    atom.fpos(1)=-x15;atom.fpos(2)=y15;atom.fpos(3)=((1.0/2.0)+z15);                     // atom B52
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x15");atom.fpos_equation.push_back("y15");atom.fpos_equation.push_back("((1.0/2.0)+z15)");// atom B52 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B52 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B52

    atom.name="C"; atom.type=2;                                       // atom B53
    atom.fpos(1)=x16;atom.fpos(2)=y16;atom.fpos(3)=z16;                     // atom B53
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x16");atom.fpos_equation.push_back("y16");atom.fpos_equation.push_back("z16");// atom B53 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B53 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B53

    atom.name="C"; atom.type=2;                                       // atom B54
    atom.fpos(1)=-x16;atom.fpos(2)=-y16;atom.fpos(3)=z16;                     // atom B54
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x16");atom.fpos_equation.push_back("-y16");atom.fpos_equation.push_back("z16");// atom B54 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B54 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B54

    atom.name="C"; atom.type=2;                                       // atom B55
    atom.fpos(1)=x16;atom.fpos(2)=-y16;atom.fpos(3)=((1.0/2.0)+z16);                     // atom B55
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x16");atom.fpos_equation.push_back("-y16");atom.fpos_equation.push_back("((1.0/2.0)+z16)");// atom B55 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B55 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B55

    atom.name="C"; atom.type=2;                                       // atom B56
    atom.fpos(1)=-x16;atom.fpos(2)=y16;atom.fpos(3)=((1.0/2.0)+z16);                     // atom B56
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x16");atom.fpos_equation.push_back("y16");atom.fpos_equation.push_back("((1.0/2.0)+z16)");// atom B56 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B56 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B56

    atom.name="C"; atom.type=2;                                       // atom B57
    atom.fpos(1)=x17;atom.fpos(2)=y17;atom.fpos(3)=z17;                     // atom B57
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x17");atom.fpos_equation.push_back("y17");atom.fpos_equation.push_back("z17");// atom B57 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B57 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B57

    atom.name="C"; atom.type=2;                                       // atom B58
    atom.fpos(1)=-x17;atom.fpos(2)=-y17;atom.fpos(3)=z17;                     // atom B58
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x17");atom.fpos_equation.push_back("-y17");atom.fpos_equation.push_back("z17");// atom B58 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B58 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B58

    atom.name="C"; atom.type=2;                                       // atom B59
    atom.fpos(1)=x17;atom.fpos(2)=-y17;atom.fpos(3)=((1.0/2.0)+z17);                     // atom B59
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x17");atom.fpos_equation.push_back("-y17");atom.fpos_equation.push_back("((1.0/2.0)+z17)");// atom B59 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B59 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B59

    atom.name="C"; atom.type=2;                                       // atom B60
    atom.fpos(1)=-x17;atom.fpos(2)=y17;atom.fpos(3)=((1.0/2.0)+z17);                     // atom B60
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x17");atom.fpos_equation.push_back("y17");atom.fpos_equation.push_back("((1.0/2.0)+z17)");// atom B60 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B60 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B60

    atom.name="C"; atom.type=2;                                       // atom B61
    atom.fpos(1)=x18;atom.fpos(2)=y18;atom.fpos(3)=z18;                     // atom B61
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x18");atom.fpos_equation.push_back("y18");atom.fpos_equation.push_back("z18");// atom B61 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B61 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B61

    atom.name="C"; atom.type=2;                                       // atom B62
    atom.fpos(1)=-x18;atom.fpos(2)=-y18;atom.fpos(3)=z18;                     // atom B62
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x18");atom.fpos_equation.push_back("-y18");atom.fpos_equation.push_back("z18");// atom B62 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B62 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B62

    atom.name="C"; atom.type=2;                                       // atom B63
    atom.fpos(1)=x18;atom.fpos(2)=-y18;atom.fpos(3)=((1.0/2.0)+z18);                     // atom B63
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x18");atom.fpos_equation.push_back("-y18");atom.fpos_equation.push_back("((1.0/2.0)+z18)");// atom B63 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B63 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B63

    atom.name="C"; atom.type=2;                                       // atom B64
    atom.fpos(1)=-x18;atom.fpos(2)=y18;atom.fpos(3)=((1.0/2.0)+z18);                     // atom B64
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x18");atom.fpos_equation.push_back("y18");atom.fpos_equation.push_back("((1.0/2.0)+z18)");// atom B64 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B64 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B64

    atom.name="C"; atom.type=2;                                       // atom B65
    atom.fpos(1)=x19;atom.fpos(2)=y19;atom.fpos(3)=z19;                     // atom B65
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x19");atom.fpos_equation.push_back("y19");atom.fpos_equation.push_back("z19");// atom B65 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B65 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B65

    atom.name="C"; atom.type=2;                                       // atom B66
    atom.fpos(1)=-x19;atom.fpos(2)=-y19;atom.fpos(3)=z19;                     // atom B66
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x19");atom.fpos_equation.push_back("-y19");atom.fpos_equation.push_back("z19");// atom B66 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B66 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B66

    atom.name="C"; atom.type=2;                                       // atom B67
    atom.fpos(1)=x19;atom.fpos(2)=-y19;atom.fpos(3)=((1.0/2.0)+z19);                     // atom B67
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x19");atom.fpos_equation.push_back("-y19");atom.fpos_equation.push_back("((1.0/2.0)+z19)");// atom B67 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B67 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B67

    atom.name="C"; atom.type=2;                                       // atom B68
    atom.fpos(1)=-x19;atom.fpos(2)=y19;atom.fpos(3)=((1.0/2.0)+z19);                     // atom B68
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x19");atom.fpos_equation.push_back("y19");atom.fpos_equation.push_back("((1.0/2.0)+z19)");// atom B68 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B68 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B68

    atom.name="C"; atom.type=2;                                       // atom B69
    atom.fpos(1)=x20;atom.fpos(2)=y20;atom.fpos(3)=z20;                     // atom B69
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x20");atom.fpos_equation.push_back("y20");atom.fpos_equation.push_back("z20");// atom B69 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B69 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B69

    atom.name="C"; atom.type=2;                                       // atom B70
    atom.fpos(1)=-x20;atom.fpos(2)=-y20;atom.fpos(3)=z20;                     // atom B70
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x20");atom.fpos_equation.push_back("-y20");atom.fpos_equation.push_back("z20");// atom B70 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B70 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B70

    atom.name="C"; atom.type=2;                                       // atom B71
    atom.fpos(1)=x20;atom.fpos(2)=-y20;atom.fpos(3)=((1.0/2.0)+z20);                     // atom B71
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x20");atom.fpos_equation.push_back("-y20");atom.fpos_equation.push_back("((1.0/2.0)+z20)");// atom B71 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B71 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B71

    atom.name="C"; atom.type=2;                                       // atom B72
    atom.fpos(1)=-x20;atom.fpos(2)=y20;atom.fpos(3)=((1.0/2.0)+z20);                     // atom B72
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x20");atom.fpos_equation.push_back("y20");atom.fpos_equation.push_back("((1.0/2.0)+z20)");// atom B72 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B72 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B72

    atom.name="C"; atom.type=2;                                       // atom B73
    atom.fpos(1)=x21;atom.fpos(2)=y21;atom.fpos(3)=z21;                     // atom B73
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x21");atom.fpos_equation.push_back("y21");atom.fpos_equation.push_back("z21");// atom B73 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B73 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B73

    atom.name="C"; atom.type=2;                                       // atom B74
    atom.fpos(1)=-x21;atom.fpos(2)=-y21;atom.fpos(3)=z21;                     // atom B74
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x21");atom.fpos_equation.push_back("-y21");atom.fpos_equation.push_back("z21");// atom B74 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B74 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B74

    atom.name="C"; atom.type=2;                                       // atom B75
    atom.fpos(1)=x21;atom.fpos(2)=-y21;atom.fpos(3)=((1.0/2.0)+z21);                     // atom B75
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x21");atom.fpos_equation.push_back("-y21");atom.fpos_equation.push_back("((1.0/2.0)+z21)");// atom B75 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B75 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B75

    atom.name="C"; atom.type=2;                                       // atom B76
    atom.fpos(1)=-x21;atom.fpos(2)=y21;atom.fpos(3)=((1.0/2.0)+z21);                     // atom B76
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x21");atom.fpos_equation.push_back("y21");atom.fpos_equation.push_back("((1.0/2.0)+z21)");// atom B76 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B76 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B76

    atom.name="C"; atom.type=2;                                       // atom B77
    atom.fpos(1)=x22;atom.fpos(2)=y22;atom.fpos(3)=z22;                     // atom B77
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x22");atom.fpos_equation.push_back("y22");atom.fpos_equation.push_back("z22");// atom B77 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B77 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B77

    atom.name="C"; atom.type=2;                                       // atom B78
    atom.fpos(1)=-x22;atom.fpos(2)=-y22;atom.fpos(3)=z22;                     // atom B78
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x22");atom.fpos_equation.push_back("-y22");atom.fpos_equation.push_back("z22");// atom B78 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B78 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B78

    atom.name="C"; atom.type=2;                                       // atom B79
    atom.fpos(1)=x22;atom.fpos(2)=-y22;atom.fpos(3)=((1.0/2.0)+z22);                     // atom B79
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x22");atom.fpos_equation.push_back("-y22");atom.fpos_equation.push_back("((1.0/2.0)+z22)");// atom B79 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B79 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B79

    atom.name="C"; atom.type=2;                                       // atom B80
    atom.fpos(1)=-x22;atom.fpos(2)=y22;atom.fpos(3)=((1.0/2.0)+z22);                     // atom B80
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x22");atom.fpos_equation.push_back("y22");atom.fpos_equation.push_back("((1.0/2.0)+z22)");// atom B80 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B80 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B80

    atom.name="C"; atom.type=2;                                       // atom B81
    atom.fpos(1)=x23;atom.fpos(2)=y23;atom.fpos(3)=z23;                     // atom B81
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x23");atom.fpos_equation.push_back("y23");atom.fpos_equation.push_back("z23");// atom B81 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B81 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B81

    atom.name="C"; atom.type=2;                                       // atom B82
    atom.fpos(1)=-x23;atom.fpos(2)=-y23;atom.fpos(3)=z23;                     // atom B82
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x23");atom.fpos_equation.push_back("-y23");atom.fpos_equation.push_back("z23");// atom B82 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B82 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B82

    atom.name="C"; atom.type=2;                                       // atom B83
    atom.fpos(1)=x23;atom.fpos(2)=-y23;atom.fpos(3)=((1.0/2.0)+z23);                     // atom B83
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x23");atom.fpos_equation.push_back("-y23");atom.fpos_equation.push_back("((1.0/2.0)+z23)");// atom B83 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B83 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B83

    atom.name="C"; atom.type=2;                                       // atom B84
    atom.fpos(1)=-x23;atom.fpos(2)=y23;atom.fpos(3)=((1.0/2.0)+z23);                     // atom B84
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x23");atom.fpos_equation.push_back("y23");atom.fpos_equation.push_back("((1.0/2.0)+z23)");// atom B84 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B84 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B84

    atom.name="C"; atom.type=2;                                       // atom B85
    atom.fpos(1)=x24;atom.fpos(2)=y24;atom.fpos(3)=z24;                     // atom B85
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x24");atom.fpos_equation.push_back("y24");atom.fpos_equation.push_back("z24");// atom B85 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B85 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B85

    atom.name="C"; atom.type=2;                                       // atom B86
    atom.fpos(1)=-x24;atom.fpos(2)=-y24;atom.fpos(3)=z24;                     // atom B86
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x24");atom.fpos_equation.push_back("-y24");atom.fpos_equation.push_back("z24");// atom B86 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B86 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B86

    atom.name="C"; atom.type=2;                                       // atom B87
    atom.fpos(1)=x24;atom.fpos(2)=-y24;atom.fpos(3)=((1.0/2.0)+z24);                     // atom B87
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x24");atom.fpos_equation.push_back("-y24");atom.fpos_equation.push_back("((1.0/2.0)+z24)");// atom B87 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B87 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B87

    atom.name="C"; atom.type=2;                                       // atom B88
    atom.fpos(1)=-x24;atom.fpos(2)=y24;atom.fpos(3)=((1.0/2.0)+z24);                     // atom B88
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x24");atom.fpos_equation.push_back("y24");atom.fpos_equation.push_back("((1.0/2.0)+z24)");// atom B88 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B88 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B88

    atom.name="C"; atom.type=2;                                       // atom B89
    atom.fpos(1)=x25;atom.fpos(2)=y25;atom.fpos(3)=z25;                     // atom B89
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x25");atom.fpos_equation.push_back("y25");atom.fpos_equation.push_back("z25");// atom B89 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B89 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B89

    atom.name="C"; atom.type=2;                                       // atom B90
    atom.fpos(1)=-x25;atom.fpos(2)=-y25;atom.fpos(3)=z25;                     // atom B90
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x25");atom.fpos_equation.push_back("-y25");atom.fpos_equation.push_back("z25");// atom B90 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B90 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B90

    atom.name="C"; atom.type=2;                                       // atom B91
    atom.fpos(1)=x25;atom.fpos(2)=-y25;atom.fpos(3)=((1.0/2.0)+z25);                     // atom B91
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x25");atom.fpos_equation.push_back("-y25");atom.fpos_equation.push_back("((1.0/2.0)+z25)");// atom B91 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B91 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B91

    atom.name="C"; atom.type=2;                                       // atom B92
    atom.fpos(1)=-x25;atom.fpos(2)=y25;atom.fpos(3)=((1.0/2.0)+z25);                     // atom B92
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x25");atom.fpos_equation.push_back("y25");atom.fpos_equation.push_back("((1.0/2.0)+z25)");// atom B92 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B92 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B92

    atom.name="C"; atom.type=2;                                       // atom B93
    atom.fpos(1)=x26;atom.fpos(2)=y26;atom.fpos(3)=z26;                     // atom B93
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x26");atom.fpos_equation.push_back("y26");atom.fpos_equation.push_back("z26");// atom B93 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B93 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B93

    atom.name="C"; atom.type=2;                                       // atom B94
    atom.fpos(1)=-x26;atom.fpos(2)=-y26;atom.fpos(3)=z26;                     // atom B94
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x26");atom.fpos_equation.push_back("-y26");atom.fpos_equation.push_back("z26");// atom B94 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B94 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B94

    atom.name="C"; atom.type=2;                                       // atom B95
    atom.fpos(1)=x26;atom.fpos(2)=-y26;atom.fpos(3)=((1.0/2.0)+z26);                     // atom B95
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x26");atom.fpos_equation.push_back("-y26");atom.fpos_equation.push_back("((1.0/2.0)+z26)");// atom B95 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B95 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B95

    atom.name="C"; atom.type=2;                                       // atom B96
    atom.fpos(1)=-x26;atom.fpos(2)=y26;atom.fpos(3)=((1.0/2.0)+z26);                     // atom B96
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x26");atom.fpos_equation.push_back("y26");atom.fpos_equation.push_back("((1.0/2.0)+z26)");// atom B96 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B96 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B96

    atom.name="C"; atom.type=2;                                       // atom B97
    atom.fpos(1)=x27;atom.fpos(2)=y27;atom.fpos(3)=z27;                     // atom B97
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x27");atom.fpos_equation.push_back("y27");atom.fpos_equation.push_back("z27");// atom B97 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B97 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B97

    atom.name="C"; atom.type=2;                                       // atom B98
    atom.fpos(1)=-x27;atom.fpos(2)=-y27;atom.fpos(3)=z27;                     // atom B98
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x27");atom.fpos_equation.push_back("-y27");atom.fpos_equation.push_back("z27");// atom B98 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B98 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B98

    atom.name="C"; atom.type=2;                                       // atom B99
    atom.fpos(1)=x27;atom.fpos(2)=-y27;atom.fpos(3)=((1.0/2.0)+z27);                     // atom B99
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x27");atom.fpos_equation.push_back("-y27");atom.fpos_equation.push_back("((1.0/2.0)+z27)");// atom B99 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B99 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B99

    atom.name="C"; atom.type=2;                                       // atom B100
    atom.fpos(1)=-x27;atom.fpos(2)=y27;atom.fpos(3)=((1.0/2.0)+z27);                     // atom B100
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x27");atom.fpos_equation.push_back("y27");atom.fpos_equation.push_back("((1.0/2.0)+z27)");// atom B100 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B100 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B100

    atom.name="C"; atom.type=2;                                       // atom B101
    atom.fpos(1)=x28;atom.fpos(2)=y28;atom.fpos(3)=z28;                     // atom B101
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x28");atom.fpos_equation.push_back("y28");atom.fpos_equation.push_back("z28");// atom B101 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B101 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B101

    atom.name="C"; atom.type=2;                                       // atom B102
    atom.fpos(1)=-x28;atom.fpos(2)=-y28;atom.fpos(3)=z28;                     // atom B102
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x28");atom.fpos_equation.push_back("-y28");atom.fpos_equation.push_back("z28");// atom B102 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B102 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B102

    atom.name="C"; atom.type=2;                                       // atom B103
    atom.fpos(1)=x28;atom.fpos(2)=-y28;atom.fpos(3)=((1.0/2.0)+z28);                     // atom B103
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x28");atom.fpos_equation.push_back("-y28");atom.fpos_equation.push_back("((1.0/2.0)+z28)");// atom B103 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B103 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B103

    atom.name="C"; atom.type=2;                                       // atom B104
    atom.fpos(1)=-x28;atom.fpos(2)=y28;atom.fpos(3)=((1.0/2.0)+z28);                     // atom B104
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x28");atom.fpos_equation.push_back("y28");atom.fpos_equation.push_back("((1.0/2.0)+z28)");// atom B104 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B104 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B104

    atom.name="D"; atom.type=3;                                       // atom B105
    atom.fpos(1)=x29;atom.fpos(2)=y29;atom.fpos(3)=z29;                     // atom B105
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x29");atom.fpos_equation.push_back("y29");atom.fpos_equation.push_back("z29");// atom B105 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B105 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B105

    atom.name="D"; atom.type=3;                                       // atom B106
    atom.fpos(1)=-x29;atom.fpos(2)=-y29;atom.fpos(3)=z29;                     // atom B106
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x29");atom.fpos_equation.push_back("-y29");atom.fpos_equation.push_back("z29");// atom B106 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B106 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B106

    atom.name="D"; atom.type=3;                                       // atom B107
    atom.fpos(1)=x29;atom.fpos(2)=-y29;atom.fpos(3)=((1.0/2.0)+z29);                     // atom B107
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x29");atom.fpos_equation.push_back("-y29");atom.fpos_equation.push_back("((1.0/2.0)+z29)");// atom B107 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B107 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B107

    atom.name="D"; atom.type=3;                                       // atom B108
    atom.fpos(1)=-x29;atom.fpos(2)=y29;atom.fpos(3)=((1.0/2.0)+z29);                     // atom B108
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x29");atom.fpos_equation.push_back("y29");atom.fpos_equation.push_back("((1.0/2.0)+z29)");// atom B108 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B108 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B108


    return str.atoms.size();  
  }
} // namespace anrl

namespace anrl {
  uint WebANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e(stringstream& web,bool LDEBUG) {
#ifndef _ANRL_NOWEB_
#endif

    if(LDEBUG) {cerr << "anrl:: WebANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e: web.str().size()=" << web.str().size() << endl;}

    return web.str().size();
  }
} // namespace anrl

#endif // AFLOW_REMOVE_GREP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

