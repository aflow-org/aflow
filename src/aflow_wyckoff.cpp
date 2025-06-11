// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 1994-2011
// Positions taken from the Bilbao Crystallographic Database

#ifndef _WYCKOFF_CPP_
#define _WYCKOFF_CPP_

#include "aflow.h"

// all vector./matrix must be loaded before

void SpaceGroupOptionRequired(uint &spacegroup, uint &option) {
  cerr << __AFLOW_FUNC__ << " Wyckoff Spacegroup " << spacegroup << " requires option 1 or 2 (" << option << ")" << endl;
  if(option==0 || option>3) {
    option=1;
    cerr << __AFLOW_FUNC__ << " Wyckoff Spacegroup " << spacegroup << " taking option=" << option << " (let`s hope it is the right one, check the concentrations and space-group)" << endl;
  }
  if(option==3) {
    option=2;
    cerr << __AFLOW_FUNC__ << " Wyckoff Spacegroup " << spacegroup << " taking option=" << option << " (let`s hope it is the right one, check the concentrations and space-group)" << endl;
  }
}

bool SpaceGroupOptionRequired(uint sg) {
  if(sg==3   || sg==4   || sg==5   || sg==6   || sg==7   ||   sg==8) return TRUE;
  if(sg==9   || sg==10  || sg==11  || sg==12  || sg==13  ||  sg==14) return TRUE;
  if(sg==15  || sg==48  || sg==50  || sg==59  || sg==68  ||  sg==70) return TRUE;
  if(sg==85  || sg==86  || sg==88  || sg==125 || sg==126 || sg==129) return TRUE;
  if(sg==130 || sg==133 || sg==134 || sg==137 || sg==138 || sg==141) return TRUE;
  if(sg==142 || sg==146 || sg==148 || sg==155 || sg==160 || sg==161) return TRUE;
  if(sg==166 || sg==167 || sg==201 || sg==203 || sg==222 || sg==224) return TRUE;
  if(sg==227 || sg==228) return TRUE;
  return FALSE;
}


// ----------------------------------------------------------------------------
xvector<double> wv(const double &x,const double &y,const double &z) {
  xvector<double> out(3);
  out[1]=x;out[2]=y;out[3]=z;
  return out;
}

// ----------------------------------------------------------------------------
void wa(_atom& a,xstructure &str) {
  a=BringInCell(a,str.lattice);
  a.cpos=F2C(str.lattice,a.fpos);
  str.AddAtom(a,false); //CO20230319 - add by type
}

// ----------------------------------------------------------------------------
xstructure WyckoffPOSITIONS(uint spacegroup, xstructure strin) {
  return  WyckoffPOSITIONS(spacegroup,(uint)1, strin);
}

xstructure WyckoffPOSITIONS(uint spacegroup_in, uint option_in, xstructure strin) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  if (spacegroup_in < 1 || spacegroup_in > 230) {
    string message = "Invalid space group " + aurostd::utype2string<uint>(spacegroup_in);
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message,_VALUE_RANGE_);
  }

  if(LDEBUG) cerr << "WyckoffPOSITIONS [0]" << endl;

  xvector<double> o(3);
  xstructure str(strin);
  double x=0,y=0,z=0;
  _atom a;
  uint spacegroup=spacegroup_in;
  uint option=option_in;

  while (str.atoms.size()>0) { str.RemoveAtom(0); } // strip all atoms
  // str.species.clear();str.species_pp.clear();str.species_pp_type.clear();str.species_pp_version.clear();str.species_pp_ZVAL.clear();species_pp_vLDAU.clear();str.species_volume.clear();str.species_mass.clear(); // patch for RemoveAtom

  if(LDEBUG) cerr << "WyckoffPOSITIONS [1]" << endl;

  for(uint i=0;i<strin.atoms.size();i++) {
    x=strin.atoms[i].fpos[1];
    y=strin.atoms[i].fpos[2];
    z=strin.atoms[i].fpos[3];
    a=strin.atoms[i];

    //     //(0,0,0)(1./2,1./2,0)
    //     for(uint j=1;j<=2;j++) {
    //       if(j==1) o=wv(0,0,0);
    //       if(j==2) o=wv(1./2,1./2,0);
    //     }
    //     //(0,0,0)(0,1./2,1./2)
    //     for(uint j=1;j<=2;j++) {
    //       if(j==1) o=wv(0,0,0);
    //       if(j==2) o=wv(0,1./2,1./2);
    //     }
    //     //(0,0,0)(1./2,1./2,1./2)
    //     for(uint j=1;j<=2;j++) {
    //       if(j==1) o=wv(0,0,0);
    //       if(j==2) o=wv(1./2,1./2,1./2);
    //     }
    //     //(0,0,0)(2./3,1./3,1./3)(1./3,2./3,2./3)
    //     for(uint j=1;j<=3;j++) {
    //       if(j==1) o=wv(0,0,0);
    //       if(j==2) o=wv(2./3,1./3,1./3);
    //       if(j==3) o=wv(1./3,2./3,2./3);
    //     }
    //     //(0,0,0)(0,1./2,1./2)(1./2,0,1./2)(1./2,1./2,0)
    //     //(0,0,0)(1./2,0,1./2)(0,1./2,1./2)(1./2,1./2,0)
    //     for(uint j=1;j<=4;j++) {
    //       if(j==1) o=wv(0,0,0);
    //       if(j==2) o=wv(0,1./2,1./2);
    //       if(j==3) o=wv(1./2,0,1./2);
    //       if(j==4) o=wv(1./2,1./2,0);
    //     }

    // for all spacegroups
    o=wv(0,0,0);

    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    if(option!=1 && option!=2) {
      if(spacegroup::SpaceGroupOptionRequired(spacegroup)) SpaceGroupOptionRequired(spacegroup,option);
      //    cerr << spacegroup << " " << spacegroup::SpaceGroupOptionRequired(spacegroup) << endl;
      // if(SpaceGroupOptionRequired(spacegroup)) SpaceGroupOptionRequired(spacegroup,option);
    }
    if(spacegroup::SpaceGroupOptionRequired(spacegroup)) str.spacegroupnumberoption=option;    // plug option inside   
    switch (spacegroup) {
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 1  P1 #1
      case 1:
        str.spacegroup="P1";
        str.spacegrouplabel="#1";
        str.spacegroupnumber=1;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 2 P-1 #2
      case 2:
        str.spacegroup="P-1";
        str.spacegrouplabel="#2";
        str.spacegroupnumber=2;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 3  P2 #3
      case 3:
        if(option==1) {
          str.spacegroup="P2";
          str.spacegrouplabel="#3";
          str.spacegroupnumber=3;
          str.spacegroupoption="unique axis b";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
        } else {
          str.spacegroup="P2";
          str.spacegrouplabel="#3";
          str.spacegroupnumber=3;
          str.spacegroupoption="unique axis c";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 4  P2_{1} #4
      case 4:
        if(option==1) {
          str.spacegroup="P2_{1}";
          str.spacegrouplabel="#4";
          str.spacegroupnumber=4;
          str.spacegroupoption="unique axis b";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,-z);wa(a,str);
        } else {
          str.spacegroup="P2_{1}";
          str.spacegrouplabel="#4";
          str.spacegroupnumber=4;
          str.spacegroupoption="unique axis c";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 5  C2 #5
      case 5:
        if(option==1) {
          str.spacegroup="C2";
          str.spacegrouplabel="#5";
          str.spacegroupnumber=5;
          str.spacegroupoption="unique axis b";
          for(uint j=1;j<=2;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(1./2,1./2,0);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x,y,-z);wa(a,str);
          }
        } else {
          str.spacegroup="C2";
          str.spacegrouplabel="#5";
          str.spacegroupnumber=5;
          str.spacegroupoption="unique axis c";
          for(uint j=1;j<=2;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(0,1./2,1./2);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x,-y,z);wa(a,str);
          }
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 6  Pm #6
      case 6:
        if(option==1) {
          str.spacegroup="Pm";
          str.spacegrouplabel="#6";
          str.spacegroupnumber=6;
          str.spacegroupoption="unique axis b";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(x,-y,z);wa(a,str);
        } else {
          str.spacegroup="Pm";
          str.spacegrouplabel="#6";
          str.spacegroupnumber=6;
          str.spacegroupoption="unique axis c";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(x,y,-z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 7  Pc #7
      case 7:
        if(option==1) {
          str.spacegroup="Pc";
          str.spacegrouplabel="#7";
          str.spacegroupnumber=7;
          str.spacegroupoption="unique axis b";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
        } else {
          str.spacegroup="Pc";
          str.spacegrouplabel="#7";
          str.spacegroupnumber=7;
          str.spacegroupoption="unique axis c";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(x+1./2,y,-z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 8  Cm #8
      case 8:
        if(option==1) {
          str.spacegroup="Cm";
          str.spacegrouplabel="#8";
          str.spacegroupnumber=8;
          str.spacegroupoption="unique axis b";
          for(uint j=1;j<=2;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(1./2,1./2,0);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(x,-y,z);wa(a,str);
          }
        } else {
          str.spacegroup="Cm";
          str.spacegrouplabel="#8";
          str.spacegroupnumber=8;
          str.spacegroupoption="unique axis c";
          for(uint j=1;j<=2;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(0,1./2,1./2);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(x,y,-z);wa(a,str);
          }
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 9  Cc #9
      case 9:
        if(option==1) {
          str.spacegroup="Cc";
          str.spacegrouplabel="#9";
          str.spacegroupnumber=9;
          str.spacegroupoption="unique axis b";
          for(uint j=1;j<=2;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(1./2,1./2,0);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
          }
        } else {
          str.spacegroup="Cc";
          str.spacegrouplabel="#9";
          str.spacegroupnumber=9;
          str.spacegroupoption="unique axis c";
          for(uint j=1;j<=2;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(0,1./2,1./2);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(x+1./2,y,-z);wa(a,str);
          }
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 10  P2./m #10
      case 10:
        if(option==1) {
          str.spacegroup="P2./m";
          str.spacegrouplabel="#10";
          str.spacegroupnumber=10;
          str.spacegroupoption="unique axis b";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,z);wa(a,str);
        } else {
          str.spacegroup="P2./m";
          str.spacegrouplabel="#10";
          str.spacegroupnumber=10;
          str.spacegroupoption="unique axis c";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,y,-z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 11  P2_{1}./m #11
      case 11:
        if(option==1) {
          str.spacegroup="P2_{1}./m";
          str.spacegrouplabel="#11";
          str.spacegroupnumber=11;
          str.spacegroupoption="unique axis b";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,z);wa(a,str);
        } else {
          str.spacegroup="P2_{1}./m";
          str.spacegrouplabel="#11";
          str.spacegroupnumber=11;
          str.spacegroupoption="unique axis c";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,y,-z+1./2);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 12  C2./m #12
      case 12:
        if(option==1) {
          str.spacegroup="C2./m";
          str.spacegrouplabel="#12";
          str.spacegroupnumber=12;
          str.spacegroupoption="unique axis b";
          for(uint j=1;j<=2;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(1./2,1./2,0);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x,y,-z);wa(a,str);
            a.fpos=o+wv(-x,-y,-z);wa(a,str);
            a.fpos=o+wv(x,-y,z);wa(a,str);
          }
        } else {
          str.spacegroup="C2./m";
          str.spacegrouplabel="#12";
          str.spacegroupnumber=12;
          str.spacegroupoption="unique axis c";
          for(uint j=1;j<=2;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(0,1./2,1./2);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x,-y,z);wa(a,str);
            a.fpos=o+wv(-x,-y,-z);wa(a,str);
            a.fpos=o+wv(x,y,-z);wa(a,str);
          }
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 13  P2./c #13
      case 13:
        if(option==1) {
          str.spacegroup="P2./c";
          str.spacegrouplabel="#13";
          str.spacegroupnumber=13;
          str.spacegroupoption="unique axis b";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
        } else {
          str.spacegroup="P2./c";
          str.spacegrouplabel="#13";
          str.spacegroupnumber=13;
          str.spacegroupoption="unique axis c";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y,-z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 14  P2_{1}./c #14
      case 14:
        if(option==1) {
          str.spacegroup="P2_{1}./c";
          str.spacegrouplabel="#14";
          str.spacegroupnumber=14;
          str.spacegroupoption="unique axis b";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
        } else {
          str.spacegroup="P2_{1}./c";
          str.spacegrouplabel="#14";
          str.spacegroupnumber=14;
          str.spacegroupoption="unique axis c";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 15  C2./c #15
      case 15:
        if(option==1) {
          str.spacegroup="C2./c";
          str.spacegrouplabel="#15";
          str.spacegroupnumber=15;
          str.spacegroupoption="unique axis b";
          for(uint j=1;j<=2;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(1./2,1./2,0);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
            a.fpos=o+wv(-x,-y,-z);wa(a,str);
            a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
          }
        } else {
          str.spacegroup="C2./c";
          str.spacegrouplabel="#15";
          str.spacegroupnumber=15;
          str.spacegroupoption="unique axis c";
          for(uint j=1;j<=2;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(0,1./2,1./2);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x+1./2,-y,z);wa(a,str);
            a.fpos=o+wv(-x,-y,-z);wa(a,str);
            a.fpos=o+wv(x+1./2,y,-z);wa(a,str);
          }
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 16  P222 #16
      case 16:
        str.spacegroup="P222";
        str.spacegrouplabel="#16";
        str.spacegroupnumber=16;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 17  P222_{1} #17
      case 17:
        str.spacegroup="P222_{1}";
        str.spacegrouplabel="#17";
        str.spacegroupnumber=17;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(x,-y,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 18  P2_{1}2_{1}2 #18
      case 18:
        str.spacegroup="P2_{1}2_{1}2";
        str.spacegrouplabel="#18";
        str.spacegroupnumber=18;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 19  P2_{1}2_{1}2_{1} #19
      case 19:
        str.spacegroup="P2_{1}2_{1}2_{1}";
        str.spacegrouplabel="#19";
        str.spacegroupnumber=19;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 20  C222_{1} #20
      case 20:
        str.spacegroup="C222_{1}";
        str.spacegrouplabel="#20";
        str.spacegroupnumber=20;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 21  C222 #21
      case 21:
        str.spacegroup="C222";
        str.spacegrouplabel="#21";
        str.spacegroupnumber=21;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 22  F222 #22
      case 22:
        str.spacegroup="F222";
        str.spacegrouplabel="#22";
        str.spacegroupnumber=22;
        str.spacegroupoption="";
        for(uint j=1;j<=4;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(0,1./2,1./2);
          if(j==3) o=wv(1./2,0,1./2);
          if(j==4) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 23  I222 #23
      case 23:
        str.spacegroup="I222";
        str.spacegrouplabel="#23";
        str.spacegroupnumber=23;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 24  I2_{1}2_{1}2_{1} #24
      case 24:
        str.spacegroup="I2_{1}2_{1}2_{1}";
        str.spacegrouplabel="#24";
        str.spacegroupnumber=24;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 25  Pmm2 #25
      case 25:
        str.spacegroup="Pmm2";
        str.spacegrouplabel="#25";
        str.spacegroupnumber=25;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(x,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 26  Pmc2_{1} #26
      case 26:
        str.spacegroup="Pmc2_{1}";
        str.spacegrouplabel="#26";
        str.spacegroupnumber=26;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 27  Pcc2 #27
      case 27:
        str.spacegroup="Pcc2";
        str.spacegrouplabel="#27";
        str.spacegroupnumber=27;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 28  Pma2 #28
      case 28:
        str.spacegroup="Pma2";
        str.spacegrouplabel="#28";
        str.spacegroupnumber=28;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(x+1./2,-y,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,y,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 29  Pca2_{1} #29
      case 29:
        str.spacegroup="Pca2_{1}";
        str.spacegrouplabel="#29";
        str.spacegroupnumber=29;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,y,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 30  Pnc2 #30
      case 30:
        str.spacegroup="Pnc2";
        str.spacegrouplabel="#30";
        str.spacegroupnumber=30;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y+1./2,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 31  Pmn2_{1} #31
      case 31:
        str.spacegroup="Pmn2_{1}";
        str.spacegrouplabel="#31";
        str.spacegroupnumber=31;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 32  Pba2 #32
      case 32:
        str.spacegroup="Pba2";
        str.spacegrouplabel="#32";
        str.spacegroupnumber=32;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 33  Pna2_{1} #33
      case 33:
        str.spacegroup="Pna2_{1}";
        str.spacegrouplabel="#33";
        str.spacegroupnumber=33;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 34  Pnn2 #34
      case 34:
        str.spacegroup="Pnn2";
        str.spacegrouplabel="#34";
        str.spacegroupnumber=34;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 35  Cmm2 #35
      case 35:
        str.spacegroup="Cmm2";
        str.spacegrouplabel="#35";
        str.spacegroupnumber=35;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 36  Cmc2_{1} #36
      case 36:
        str.spacegroup="Cmc2_{1}";
        str.spacegrouplabel="#36";
        str.spacegroupnumber=36;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 37  Ccc2 #37
      case 37:
        str.spacegroup="Ccc2";
        str.spacegrouplabel="#37";
        str.spacegroupnumber=37;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 38  Amm2 #38
      case 38:
        str.spacegroup="Amm2";
        str.spacegrouplabel="#38";
        str.spacegroupnumber=38;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(0,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 39  Aem2 #39
      case 39:
        str.spacegroup="Aem2";
        str.spacegrouplabel="#39";
        str.spacegroupnumber=39;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(0,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 40  Ama2 #40
      case 40:
        str.spacegroup="Ama2";
        str.spacegrouplabel="#40";
        str.spacegroupnumber=40;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(0,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 41  Aea2 #41
      case 41:
        str.spacegroup="Aea2";
        str.spacegrouplabel="#41";
        str.spacegroupnumber=41;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(0,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 42  Fmm2 #42
      case 42:
        str.spacegroup="Fmm2";
        str.spacegrouplabel="#42";
        str.spacegroupnumber=42;
        str.spacegroupoption="";
        for(uint j=1;j<=4;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(0,1./2,1./2);
          if(j==3) o=wv(1./2,0,1./2);
          if(j==4) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 43  Fdd2 #43
      case 43:
        str.spacegroup="Fdd2";
        str.spacegrouplabel="#43";
        str.spacegroupnumber=43;
        str.spacegroupoption="";
        for(uint j=1;j<=4;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(0,1./2,1./2);
          if(j==3) o=wv(1./2,0,1./2);
          if(j==4) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(x+1./4,-y+1./4,z+1./4);wa(a,str);
          a.fpos=o+wv(-x+1./4,y+1./4,z+1./4);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 44  Imm2 #44
      case 44:
        str.spacegroup="Imm2";
        str.spacegrouplabel="#44";
        str.spacegroupnumber=44;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 45  Iba2 #45
      case 45:
        str.spacegroup="Iba2";
        str.spacegrouplabel="#45";
        str.spacegroupnumber=45;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 46  Ima2 #46
      case 46:
        str.spacegroup="Ima2";
        str.spacegrouplabel="#46";
        str.spacegroupnumber=46;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 47  Pmmm #47
      case 47:
        str.spacegroup="Pmmm";
        str.spacegrouplabel="#47";
        str.spacegroupnumber=47;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,-z);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 48  Pnnn #48
      case 48:
        if(option==1) {
          str.spacegroup="Pnnn";
          str.spacegrouplabel="#48";
          str.spacegroupnumber=48;
          str.spacegroupoption="origin choice 1";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
        } else {
          str.spacegroup="Pnnn";
          str.spacegrouplabel="#48";
          str.spacegroupnumber=48;
          str.spacegroupoption="origin choice 2";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y,-z+1./2);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,z+1./2);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 49  Pccm #49
      case 49:
        str.spacegroup="Pccm";
        str.spacegrouplabel="#49";
        str.spacegroupnumber=49;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(x,-y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 50  Pban #50
      case 50:
        if(option==1) {
          str.spacegroup="Pban";
          str.spacegrouplabel="#50";
          str.spacegroupnumber=50;
          str.spacegroupoption="origin choice 1";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
        } else {
          str.spacegroup="Pban";
          str.spacegrouplabel="#50";
          str.spacegroupnumber=50;
          str.spacegroupoption="origin choice 2";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,-z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 51  Pmma #51
      case 51:
        str.spacegroup="Pmma";
        str.spacegrouplabel="#51";
        str.spacegroupnumber=51;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,-y,-z);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,y,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 52  Pnna #52
      case 52:
        str.spacegroup="Pnna";
        str.spacegrouplabel="#52";
        str.spacegroupnumber=52;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,-y,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(x,-y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,y,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y+1./2,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 53  Pmna #53
      case 53:
        str.spacegroup="Pmna";
        str.spacegrouplabel="#53";
        str.spacegroupnumber=53;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(x,-y,-z);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 54  Pcca #54
      case 54:
        str.spacegroup="Pcca";
        str.spacegrouplabel="#54";
        str.spacegroupnumber=54;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,y,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 55  Pbam #55
      case 55:
        str.spacegroup="Pbam";
        str.spacegrouplabel="#55";
        str.spacegroupnumber=55;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 56  Pccn #56
      case 56:
        str.spacegroup="Pccn";
        str.spacegrouplabel="#56";
        str.spacegroupnumber=56;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
        a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
        a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,y,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 57  Pbcm #57
      case 57:
        str.spacegroup="Pbcm";
        str.spacegrouplabel="#57";
        str.spacegroupnumber=57;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(x,-y+1./2,-z);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y+1./2,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 58  Pnnm #58
      case 58:
        str.spacegroup="Pnnm";
        str.spacegrouplabel="#58";
        str.spacegroupnumber=58;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 59  Pmmn #59
      case 59:
        if(option==1) {
          str.spacegroup="Pmmn";
          str.spacegrouplabel="#59";
          str.spacegroupnumber=59;
          str.spacegroupoption="origin choice 1";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
        } else {
          str.spacegroup="Pmmn";
          str.spacegrouplabel="#59";
          str.spacegroupnumber=59;
          str.spacegroupoption="origin choice 2";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y,-z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 60  Pbcn #60
      case 60:
        str.spacegroup="Pbcn";
        str.spacegrouplabel="#60";
        str.spacegroupnumber=60;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,-y+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 61  Pbca #61
      case 61:
        str.spacegroup="Pbca";
        str.spacegrouplabel="#61";
        str.spacegroupnumber=61;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 62  Pnma #62
      case 62:
        str.spacegroup="Pnma";
        str.spacegrouplabel="#62";
        str.spacegroupnumber=62;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y+1./2,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(x,-y+1./2,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 63  Cmcm #63
      case 63:
        str.spacegroup="Cmcm";
        str.spacegrouplabel="#63";
        str.spacegroupnumber=63;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,y,-z+1./2);wa(a,str);
          a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 64  Cmce #64
      case 64:
        str.spacegroup="Cmce";
        str.spacegrouplabel="#64";
        str.spacegroupnumber=64;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 65  Cmmm #65
      case 65:
        str.spacegroup="Cmmm";
        str.spacegrouplabel="#65";
        str.spacegroupnumber=65;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 66  Cccm #66
      case 66:
        str.spacegroup="Cccm";
        str.spacegrouplabel="#66";
        str.spacegroupnumber=66;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
          a.fpos=o+wv(x,-y,-z+1./2);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 67  Cmme #67
      case 67:
        str.spacegroup="Cmme";
        str.spacegrouplabel="#67";
        str.spacegroupnumber=67;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 68  Ccce #68
      case 68:
        if(option==1) {
          str.spacegroup="Ccce";
          str.spacegrouplabel="#68";
          str.spacegroupnumber=68;
          str.spacegroupoption="origin choice 1";
          for(uint j=1;j<=2;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(1./2,1./2,0);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
            a.fpos=o+wv(-x,y,-z);wa(a,str);
            a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
            a.fpos=o+wv(-x,-y+1./2,-z+1./2);wa(a,str);
            a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
            a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
            a.fpos=o+wv(-x+1./2,y,z+1./2);wa(a,str);
          }
        } else {
          str.spacegroup="Ccce";
          str.spacegrouplabel="#68";
          str.spacegroupnumber=68;
          str.spacegroupoption="origin choice 2";
          for(uint j=1;j<=2;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(1./2,1./2,0);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x+1./2,-y,z);wa(a,str);
            a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
            a.fpos=o+wv(x+1./2,-y,-z+1./2);wa(a,str);
            a.fpos=o+wv(-x,-y,-z);wa(a,str);
            a.fpos=o+wv(x+1./2,y,-z);wa(a,str);
            a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
            a.fpos=o+wv(-x+1./2,y,z+1./2);wa(a,str);
          }
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 69  Fmmm #69
      case 69:
        str.spacegroup="Fmmm";
        str.spacegrouplabel="#69";
        str.spacegroupnumber=69;
        str.spacegroupoption="";
        for(uint j=1;j<=4;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(0,1./2,1./2);
          if(j==3) o=wv(1./2,0,1./2);
          if(j==4) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 70  Fddd #70
      case 70:
        if(option==1) {
          str.spacegroup="Fddd";
          str.spacegrouplabel="#70";
          str.spacegroupnumber=70;
          str.spacegroupoption="origin choice 1";
          for(uint j=1;j<=4;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(0,1./2,1./2);
            if(j==3) o=wv(1./2,0,1./2);
            if(j==4) o=wv(1./2,1./2,0);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x,-y,z);wa(a,str);
            a.fpos=o+wv(-x,y,-z);wa(a,str);
            a.fpos=o+wv(x,-y,-z);wa(a,str);
            a.fpos=o+wv(-x+1./4,-y+1./4,-z+1./4);wa(a,str);
            a.fpos=o+wv(x+1./4,y+1./4,-z+1./4);wa(a,str);
            a.fpos=o+wv(x+1./4,-y+1./4,z+1./4);wa(a,str);
            a.fpos=o+wv(-x+1./4,y+1./4,z+1./4);wa(a,str);
          }
        } else {
          str.spacegroup="Fddd";
          str.spacegrouplabel="#70";
          str.spacegroupnumber=70;
          str.spacegroupoption="origin choice 2";
          for(uint j=1;j<=4;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(0,1./2,1./2);
            if(j==3) o=wv(1./2,0,1./2);
            if(j==4) o=wv(1./2,1./2,0);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x+3./4,-y+3./4,z);wa(a,str);
            a.fpos=o+wv(-x+3./4,y,-z+3./4);wa(a,str);
            a.fpos=o+wv(x,-y+3./4,-z+3./4);wa(a,str);
            a.fpos=o+wv(-x,-y,-z);wa(a,str);
            a.fpos=o+wv(x+1./4,y+1./4,-z);wa(a,str);
            a.fpos=o+wv(x+1./4,-y,z+1./4);wa(a,str);
            a.fpos=o+wv(-x,y+1./4,z+1./4);wa(a,str);
          }
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 71  Immm #71
      case 71:
        str.spacegroup="Immm";
        str.spacegrouplabel="#71";
        str.spacegroupnumber=71;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 72  Ibam #72
      case 72:
        str.spacegroup="Ibam";
        str.spacegrouplabel="#72";
        str.spacegroupnumber=72;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 73  Ibca #73
      case 73:
        str.spacegroup="Ibca";
        str.spacegrouplabel="#73";
        str.spacegroupnumber=73;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 74  Imma #74
      case 74:
        str.spacegroup="Imma";
        str.spacegrouplabel="#74";
        str.spacegroupnumber=74;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 75  P4 #75
      case 75:
        str.spacegroup="P4";
        str.spacegrouplabel="#75";
        str.spacegroupnumber=75;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y,x,z);wa(a,str);
        a.fpos=o+wv(y,-x,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 76  P4_{1} #76
      case 76:
        str.spacegroup="P4_{1}";
        str.spacegrouplabel="#76";
        str.spacegroupnumber=76;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-y,x,z+1./4);wa(a,str);
        a.fpos=o+wv(y,-x,z+3./4);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 77  P4_{2} #77
      case 77:
        str.spacegroup="P4_{2}";
        str.spacegrouplabel="#77";
        str.spacegroupnumber=77;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(y,-x,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 78  P4_{3} #78
      case 78:
        str.spacegroup="P4_{3}";
        str.spacegrouplabel="#78";
        str.spacegroupnumber=78;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-y,x,z+3./4);wa(a,str);
        a.fpos=o+wv(y,-x,z+1./4);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 79  I4 #79
      case 79:
        str.spacegroup="I4";
        str.spacegrouplabel="#79";
        str.spacegroupnumber=79;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-y,x,z);wa(a,str);
          a.fpos=o+wv(y,-x,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 80  I4_{1} #80
      case 80:
        str.spacegroup="I4_{1}";
        str.spacegrouplabel="#80";
        str.spacegroupnumber=80;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-y,x+1./2,z+1./4);wa(a,str);
          a.fpos=o+wv(y+1./2,-x,z+3./4);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 81  P-4 #81
      case 81:
        str.spacegroup="P-4";
        str.spacegrouplabel="#81";
        str.spacegroupnumber=81;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(y,-x,-z);wa(a,str);
        a.fpos=o+wv(-y,x,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 82  I-4 #82
      case 82:
        str.spacegroup="I-4";
        str.spacegrouplabel="#82";
        str.spacegroupnumber=82;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(y,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x,-z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 83  P4./m #83
      case 83:
        str.spacegroup="P4./m";
        str.spacegrouplabel="#83";
        str.spacegroupnumber=83;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y,x,z);wa(a,str);
        a.fpos=o+wv(y,-x,z);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(y,-x,-z);wa(a,str);
        a.fpos=o+wv(-y,x,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 84  P4_{2}./m #84
      case 84:
        str.spacegroup="P4_{2}./m";
        str.spacegrouplabel="#84";
        str.spacegroupnumber=84;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(y,-x,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(y,-x,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y,x,-z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 85  P4./n #85
      case 85:
        if(option==1) {
          str.spacegroup="P4./n";
          str.spacegrouplabel="#85";
          str.spacegroupnumber=85;
          str.spacegroupoption="origin choice 1";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,x+1./2,z);wa(a,str);
          a.fpos=o+wv(y+1./2,-x+1./2,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(y,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x,-z);wa(a,str);
        } else {
          str.spacegroup="P4./n";
          str.spacegrouplabel="#85";
          str.spacegroupnumber=85;
          str.spacegroupoption="origin choice 2";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,x,z);wa(a,str);
          a.fpos=o+wv(y,-x+1./2,z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(y+1./2,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x+1./2,-z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 86  P4_{2}./n #86
      case 86:
        if(option==1) {
          str.spacegroup="P4_{2}./n";
          str.spacegrouplabel="#86";
          str.spacegroupnumber=86;
          str.spacegroupoption="origin choice 1";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(y,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x,-z);wa(a,str);
        } else {
          str.spacegroup="P4_{2}./n";
          str.spacegrouplabel="#86";
          str.spacegroupnumber=86;
          str.spacegroupoption="origin choice 2";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-y,x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,-x,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(y,-x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,x,-z+1./2);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 87  I4./m #87
      case 87:
        str.spacegroup="I4./m";
        str.spacegrouplabel="#87";
        str.spacegroupnumber=87;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-y,x,z);wa(a,str);
          a.fpos=o+wv(y,-x,z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,y,-z);wa(a,str);
          a.fpos=o+wv(y,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x,-z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 88  I4_{1}./a #88
      case 88:
        if(option==1) {
          str.spacegroup="I4_{1}./a";
          str.spacegrouplabel="#88";
          str.spacegroupnumber=88;
          str.spacegroupoption="origin choice 1";
          for(uint j=1;j<=2;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(1./2,1./2,1./2);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x+1./2,-y+1./2,z+1./2);wa(a,str);
            a.fpos=o+wv(-y,x+1./2,z+1./4);wa(a,str);
            a.fpos=o+wv(y+1./2,-x,z+3./4);wa(a,str);
            a.fpos=o+wv(-x,-y+1./2,-z+1./4);wa(a,str);
            a.fpos=o+wv(x+1./2,y,-z+3./4);wa(a,str);
            a.fpos=o+wv(y,-x,-z);wa(a,str);
            a.fpos=o+wv(-y+1./2,x+1./2,-z+1./2);wa(a,str);
          }
        } else {
          str.spacegroup="I4_{1}./a";
          str.spacegrouplabel="#88";
          str.spacegroupnumber=88;
          str.spacegroupoption="origin choice 2";
          for(uint j=1;j<=2;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(1./2,1./2,1./2);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
            a.fpos=o+wv(-y+3./4,x+1./4,z+1./4);wa(a,str);
            a.fpos=o+wv(y+3./4,-x+3./4,z+3./4);wa(a,str);
            a.fpos=o+wv(-x,-y,-z);wa(a,str);
            a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
            a.fpos=o+wv(y+1./4,-x+3./4,-z+3./4);wa(a,str);
            a.fpos=o+wv(-y+1./4,x+1./4,-z+1./4);wa(a,str);
          }
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 89  P422 #89
      case 89:
        str.spacegroup="P422";
        str.spacegrouplabel="#89";
        str.spacegroupnumber=89;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y,x,z);wa(a,str);
        a.fpos=o+wv(y,-x,z);wa(a,str);
        a.fpos=o+wv(-x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,-z);wa(a,str);
        a.fpos=o+wv(y,x,-z);wa(a,str);
        a.fpos=o+wv(-y,-x,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 90  P42_{1}2 #90
      case 90:
        str.spacegroup="P42_{1}2";
        str.spacegrouplabel="#90";
        str.spacegroupnumber=90;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y+1./2,x+1./2,z);wa(a,str);
        a.fpos=o+wv(y+1./2,-x+1./2,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
        a.fpos=o+wv(y,x,-z);wa(a,str);
        a.fpos=o+wv(-y,-x,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 91  P4_{1}22 #91
      case 91:
        str.spacegroup="P4_{1}22";
        str.spacegrouplabel="#91";
        str.spacegroupnumber=91;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-y,x,z+1./4);wa(a,str);
        a.fpos=o+wv(y,-x,z+3./4);wa(a,str);
        a.fpos=o+wv(-x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,-z+1./2);wa(a,str);
        a.fpos=o+wv(y,x,-z+3./4);wa(a,str);
        a.fpos=o+wv(-y,-x,-z+1./4);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 92  P4_{1}2_{1}2 #92
      case 92:
        str.spacegroup="P4_{1}2_{1}2";
        str.spacegrouplabel="#92";
        str.spacegroupnumber=92;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-y+1./2,x+1./2,z+1./4);wa(a,str);
        a.fpos=o+wv(y+1./2,-x+1./2,z+3./4);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,-z+1./4);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,-z+3./4);wa(a,str);
        a.fpos=o+wv(y,x,-z);wa(a,str);
        a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 93  P4_{2}22 #93
      case 93:
        str.spacegroup="P4_{2}22";
        str.spacegrouplabel="#93";
        str.spacegroupnumber=93;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(y,-x,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,-z);wa(a,str);
        a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 94  P4_{2}2_{1}2 #94
      case 94:
        str.spacegroup="P4_{2}2_{1}2";
        str.spacegrouplabel="#94";
        str.spacegroupnumber=94;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(y,x,-z);wa(a,str);
        a.fpos=o+wv(-y,-x,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 95  P4_{3}22 #95
      case 95:
        str.spacegroup="P4_{3}22";
        str.spacegrouplabel="#95";
        str.spacegroupnumber=95;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-y,x,z+3./4);wa(a,str);
        a.fpos=o+wv(y,-x,z+1./4);wa(a,str);
        a.fpos=o+wv(-x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,-z+1./2);wa(a,str);
        a.fpos=o+wv(y,x,-z+1./4);wa(a,str);
        a.fpos=o+wv(-y,-x,-z+3./4);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 96  P4_{3}2_{1}2 #96
      case 96:
        str.spacegroup="P4_{3}2_{1}2";
        str.spacegrouplabel="#96";
        str.spacegroupnumber=96;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-y+1./2,x+1./2,z+3./4);wa(a,str);
        a.fpos=o+wv(y+1./2,-x+1./2,z+1./4);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,-z+3./4);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,-z+1./4);wa(a,str);
        a.fpos=o+wv(y,x,-z);wa(a,str);
        a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 97  I422 #97
      case 97:
        str.spacegroup="I422";
        str.spacegrouplabel="#97";
        str.spacegroupnumber=97;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-y,x,z);wa(a,str);
          a.fpos=o+wv(y,-x,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(y,x,-z);wa(a,str);
          a.fpos=o+wv(-y,-x,-z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 98  I4_{1}22 #98
      case 98:
        str.spacegroup="I4_{1}22";
        str.spacegrouplabel="#98";
        str.spacegroupnumber=98;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-y,x+1./2,z+1./4);wa(a,str);
          a.fpos=o+wv(y+1./2,-x,z+3./4);wa(a,str);
          a.fpos=o+wv(-x+1./2,y,-z+3./4);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,-z+1./4);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(-y,-x,-z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 99  P4mm #99
      case 99:
        str.spacegroup="P4mm";
        str.spacegrouplabel="#99";
        str.spacegroupnumber=99;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y,x,z);wa(a,str);
        a.fpos=o+wv(y,-x,z);wa(a,str);
        a.fpos=o+wv(x,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,z);wa(a,str);
        a.fpos=o+wv(-y,-x,z);wa(a,str);
        a.fpos=o+wv(y,x,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 100  P4bm #100
      case 100:
        str.spacegroup="P4bm";
        str.spacegrouplabel="#100";
        str.spacegroupnumber=100;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y,x,z);wa(a,str);
        a.fpos=o+wv(y,-x,z);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
        a.fpos=o+wv(-y+1./2,-x+1./2,z);wa(a,str);
        a.fpos=o+wv(y+1./2,x+1./2,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 101  P4_{2}cm #101
      case 101:
        str.spacegroup="P4_{2}cm";
        str.spacegrouplabel="#101";
        str.spacegroupnumber=101;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(y,-x,z+1./2);wa(a,str);
        a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
        a.fpos=o+wv(-y,-x,z);wa(a,str);
        a.fpos=o+wv(y,x,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 102  P4_{2}nm #102
      case 102:
        str.spacegroup="P4_{2}nm";
        str.spacegrouplabel="#102";
        str.spacegroupnumber=102;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-y,-x,z);wa(a,str);
        a.fpos=o+wv(y,x,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 103  P4cc #103
      case 103:
        str.spacegroup="P4cc";
        str.spacegrouplabel="#103";
        str.spacegroupnumber=103;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y,x,z);wa(a,str);
        a.fpos=o+wv(y,-x,z);wa(a,str);
        a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
        a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
        a.fpos=o+wv(y,x,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 104  P4nc #104
      case 104:
        str.spacegroup="P4nc";
        str.spacegrouplabel="#104";
        str.spacegroupnumber=104;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y,x,z);wa(a,str);
        a.fpos=o+wv(y,-x,z);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 105  P4_{2}mc #105
      case 105:
        str.spacegroup="P4_{2}mc";
        str.spacegrouplabel="#105";
        str.spacegroupnumber=105;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(y,-x,z+1./2);wa(a,str);
        a.fpos=o+wv(x,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,z);wa(a,str);
        a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
        a.fpos=o+wv(y,x,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 106  P4_{2}bc #106
      case 106:
        str.spacegroup="P4_{2}bc";
        str.spacegrouplabel="#106";
        str.spacegroupnumber=106;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(y,-x,z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
        a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 107  I4mm #107
      case 107:
        str.spacegroup="I4mm";
        str.spacegrouplabel="#107";
        str.spacegroupnumber=107;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-y,x,z);wa(a,str);
          a.fpos=o+wv(y,-x,z);wa(a,str);
          a.fpos=o+wv(x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
          a.fpos=o+wv(-y,-x,z);wa(a,str);
          a.fpos=o+wv(y,x,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 108  I4cm #108
      case 108:
        str.spacegroup="I4cm";
        str.spacegrouplabel="#108";
        str.spacegroupnumber=108;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-y,x,z);wa(a,str);
          a.fpos=o+wv(y,-x,z);wa(a,str);
          a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
          a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
          a.fpos=o+wv(y,x,z+1./2);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 109  I4_{1}md #109
      case 109:
        str.spacegroup="I4_{1}md";
        str.spacegrouplabel="#109";
        str.spacegroupnumber=109;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-y,x+1./2,z+1./4);wa(a,str);
          a.fpos=o+wv(y+1./2,-x,z+3./4);wa(a,str);
          a.fpos=o+wv(x,-y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-y,-x+1./2,z+1./4);wa(a,str);
          a.fpos=o+wv(y+1./2,x,z+3./4);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 110  I4_{1}cd #110
      case 110:
        str.spacegroup="I4_{1}cd";
        str.spacegrouplabel="#110";
        str.spacegroupnumber=110;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-y,x+1./2,z+1./4);wa(a,str);
          a.fpos=o+wv(y+1./2,-x,z+3./4);wa(a,str);
          a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
          a.fpos=o+wv(-y,-x+1./2,z+3./4);wa(a,str);
          a.fpos=o+wv(y+1./2,x,z+1./4);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 111  P-42m #111
      case 111:
        str.spacegroup="P-42m";
        str.spacegrouplabel="#111";
        str.spacegroupnumber=111;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(y,-x,-z);wa(a,str);
        a.fpos=o+wv(-y,x,-z);wa(a,str);
        a.fpos=o+wv(-x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,-z);wa(a,str);
        a.fpos=o+wv(-y,-x,z);wa(a,str);
        a.fpos=o+wv(y,x,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 112  P-42c #112
      case 112:
        str.spacegroup="P-42c";
        str.spacegrouplabel="#112";
        str.spacegroupnumber=112;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(y,-x,-z);wa(a,str);
        a.fpos=o+wv(-y,x,-z);wa(a,str);
        a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(x,-y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
        a.fpos=o+wv(y,x,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 113  P-42_{1}m #113
      case 113:
        str.spacegroup="P-42_{1}m";
        str.spacegrouplabel="#113";
        str.spacegroupnumber=113;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(y,-x,-z);wa(a,str);
        a.fpos=o+wv(-y,x,-z);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
        a.fpos=o+wv(-y+1./2,-x+1./2,z);wa(a,str);
        a.fpos=o+wv(y+1./2,x+1./2,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 114  P-42_{1}c #114
      case 114:
        str.spacegroup="P-42_{1}c";
        str.spacegrouplabel="#114";
        str.spacegroupnumber=114;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(y,-x,-z);wa(a,str);
        a.fpos=o+wv(-y,x,-z);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 115  P-4m2 #115
      case 115:
        str.spacegroup="P-4m2";
        str.spacegrouplabel="#115";
        str.spacegroupnumber=115;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(y,-x,-z);wa(a,str);
        a.fpos=o+wv(-y,x,-z);wa(a,str);
        a.fpos=o+wv(x,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,z);wa(a,str);
        a.fpos=o+wv(y,x,-z);wa(a,str);
        a.fpos=o+wv(-y,-x,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 116  P-4c2 #116
      case 116:
        str.spacegroup="P-4c2";
        str.spacegrouplabel="#116";
        str.spacegroupnumber=116;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(y,-x,-z);wa(a,str);
        a.fpos=o+wv(-y,x,-z);wa(a,str);
        a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
        a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 117  P-4b2 #117
      case 117:
        str.spacegroup="P-4b2";
        str.spacegrouplabel="#117";
        str.spacegroupnumber=117;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(y,-x,-z);wa(a,str);
        a.fpos=o+wv(-y,x,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
        a.fpos=o+wv(y+1./2,x+1./2,-z);wa(a,str);
        a.fpos=o+wv(-y+1./2,-x+1./2,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 118  P-4n2 #118
      case 118:
        str.spacegroup="P-4n2";
        str.spacegrouplabel="#118";
        str.spacegroupnumber=118;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(y,-x,-z);wa(a,str);
        a.fpos=o+wv(-y,x,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 119  I-4m2 #119
      case 119:
        str.spacegroup="I-4m2";
        str.spacegrouplabel="#119";
        str.spacegroupnumber=119;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(y,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x,-z);wa(a,str);
          a.fpos=o+wv(x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
          a.fpos=o+wv(y,x,-z);wa(a,str);
          a.fpos=o+wv(-y,-x,-z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 120  I-4c2 #120
      case 120:
        str.spacegroup="I-4c2";
        str.spacegrouplabel="#120";
        str.spacegroupnumber=120;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(y,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x,-z);wa(a,str);
          a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
          a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
          a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 121  I-42m #121
      case 121:
        str.spacegroup="I-42m";
        str.spacegrouplabel="#121";
        str.spacegroupnumber=121;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(y,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x,-z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(-y,-x,z);wa(a,str);
          a.fpos=o+wv(y,x,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 122  I-42d #122
      case 122:
        str.spacegroup="I-42d";
        str.spacegrouplabel="#122";
        str.spacegroupnumber=122;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(y,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x,-z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y,-z+3./4);wa(a,str);
          a.fpos=o+wv(x+1./2,-y,-z+3./4);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x,z+3./4);wa(a,str);
          a.fpos=o+wv(y+1./2,x,z+3./4);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 123  P4./mmm #123
      case 123:
        str.spacegroup="P4./mmm";
        str.spacegrouplabel="#123";
        str.spacegroupnumber=123;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y,x,z);wa(a,str);
        a.fpos=o+wv(y,-x,z);wa(a,str);
        a.fpos=o+wv(-x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,-z);wa(a,str);
        a.fpos=o+wv(y,x,-z);wa(a,str);
        a.fpos=o+wv(-y,-x,-z);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(y,-x,-z);wa(a,str);
        a.fpos=o+wv(-y,x,-z);wa(a,str);
        a.fpos=o+wv(x,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,z);wa(a,str);
        a.fpos=o+wv(-y,-x,z);wa(a,str);
        a.fpos=o+wv(y,x,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 124  P4./mcc #124
      case 124:
        str.spacegroup="P4./mcc";
        str.spacegrouplabel="#124";
        str.spacegroupnumber=124;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y,x,z);wa(a,str);
        a.fpos=o+wv(y,-x,z);wa(a,str);
        a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(x,-y,-z+1./2);wa(a,str);
        a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(y,-x,-z);wa(a,str);
        a.fpos=o+wv(-y,x,-z);wa(a,str);
        a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
        a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
        a.fpos=o+wv(y,x,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 125  P4./nbm #125
      case 125:
        if(option==1) {
          str.spacegroup="P4./nbm";
          str.spacegrouplabel="#125";
          str.spacegroupnumber=125;
          str.spacegroupoption="origin choice 1";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-y,x,z);wa(a,str);
          a.fpos=o+wv(y,-x,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(y,x,-z);wa(a,str);
          a.fpos=o+wv(-y,-x,-z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(y+1./2,-x+1./2,-z);wa(a,str);
          a.fpos=o+wv(-y+1./2,x+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,z);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,z);wa(a,str);
        } else {
          str.spacegroup="P4./nbm";
          str.spacegrouplabel="#125";
          str.spacegroupnumber=125;
          str.spacegroupoption="origin choice 2";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,x,z);wa(a,str);
          a.fpos=o+wv(y,-x+1./2,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,-z);wa(a,str);
          a.fpos=o+wv(y,x,-z);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,-z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(y+1./2,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,z);wa(a,str);
          a.fpos=o+wv(-y,-x,z);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 126  P4./nnc #126
      case 126:
        if(option==1) {
          str.spacegroup="P4./nnc";
          str.spacegrouplabel="#126";
          str.spacegroupnumber=126;
          str.spacegroupoption="origin choice 1";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-y,x,z);wa(a,str);
          a.fpos=o+wv(y,-x,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(y,x,-z);wa(a,str);
          a.fpos=o+wv(-y,-x,-z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,-x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
        } else {
          str.spacegroup="P4./nnc";
          str.spacegrouplabel="#126";
          str.spacegroupnumber=126;
          str.spacegroupoption="origin choice 2";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,x,z);wa(a,str);
          a.fpos=o+wv(y,-x+1./2,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y,-z+1./2);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(y+1./2,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 127  P4./mbm #127
      case 127:
        str.spacegroup="P4./mbm";
        str.spacegrouplabel="#127";
        str.spacegroupnumber=127;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y,x,z);wa(a,str);
        a.fpos=o+wv(y,-x,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
        a.fpos=o+wv(y+1./2,x+1./2,-z);wa(a,str);
        a.fpos=o+wv(-y+1./2,-x+1./2,-z);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(y,-x,-z);wa(a,str);
        a.fpos=o+wv(-y,x,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
        a.fpos=o+wv(-y+1./2,-x+1./2,z);wa(a,str);
        a.fpos=o+wv(y+1./2,x+1./2,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 128  P4./mnc #128
      case 128:
        str.spacegroup="P4./mnc";
        str.spacegrouplabel="#128";
        str.spacegroupnumber=128;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y,x,z);wa(a,str);
        a.fpos=o+wv(y,-x,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(y,-x,-z);wa(a,str);
        a.fpos=o+wv(-y,x,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 129  P4./nmm #129
      case 129:
        if(option==1) {
          str.spacegroup="P4./nmm";
          str.spacegrouplabel="#129";
          str.spacegroupnumber=129;
          str.spacegroupoption="origin choice 1";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,x+1./2,z);wa(a,str);
          a.fpos=o+wv(y+1./2,-x+1./2,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
          a.fpos=o+wv(y,x,-z);wa(a,str);
          a.fpos=o+wv(-y,-x,-z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(y,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x,-z);wa(a,str);
          a.fpos=o+wv(x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,z);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,z);wa(a,str);
        } else {
          str.spacegroup="P4./nmm";
          str.spacegrouplabel="#129";
          str.spacegroupnumber=129;
          str.spacegroupoption="origin choice 2";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,x,z);wa(a,str);
          a.fpos=o+wv(y,-x+1./2,z);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y,-z);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,-z);wa(a,str);
          a.fpos=o+wv(-y,-x,-z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(y+1./2,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x+1./2,-z);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,z);wa(a,str);
          a.fpos=o+wv(y,x,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 130  P4./ncc #130
      case 130:
        if(option==1) {
          str.spacegroup="P4./ncc";
          str.spacegrouplabel="#130";
          str.spacegroupnumber=130;
          str.spacegroupoption="origin choice 1";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,x+1./2,z);wa(a,str);
          a.fpos=o+wv(y+1./2,-x+1./2,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
          a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(y,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x,-z);wa(a,str);
          a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
        } else {
          str.spacegroup="P4./ncc";
          str.spacegrouplabel="#130";
          str.spacegroupnumber=130;
          str.spacegroupoption="origin choice 2";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,x,z);wa(a,str);
          a.fpos=o+wv(y,-x+1./2,z);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-y,-z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(y+1./2,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x+1./2,-z);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,y,z+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(y,x,z+1./2);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 131  P4_{2}./mmc #131
      case 131:
        str.spacegroup="P4_{2}./mmc";
        str.spacegrouplabel="#131";
        str.spacegroupnumber=131;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(y,-x,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,-z);wa(a,str);
        a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(y,-x,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y,x,-z+1./2);wa(a,str);
        a.fpos=o+wv(x,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,z);wa(a,str);
        a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
        a.fpos=o+wv(y,x,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 132  P4_{2}./mcm #132
      case 132:
        str.spacegroup="P4_{2}./mcm";
        str.spacegrouplabel="#132";
        str.spacegroupnumber=132;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(y,-x,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(x,-y,-z+1./2);wa(a,str);
        a.fpos=o+wv(y,x,-z);wa(a,str);
        a.fpos=o+wv(-y,-x,-z);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(y,-x,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y,x,-z+1./2);wa(a,str);
        a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
        a.fpos=o+wv(-y,-x,z);wa(a,str);
        a.fpos=o+wv(y,x,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 133  P4_{2}./nbc #133
      case 133:
        if(option==1) {
          str.spacegroup="P4_{2}./nbc";
          str.spacegrouplabel="#133";
          str.spacegroupnumber=133;
          str.spacegroupoption="origin choice 1";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
          a.fpos=o+wv(x,-y,-z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,-z);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,-z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(y,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
          a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
          a.fpos=o+wv(y,x,z+1./2);wa(a,str);
        } else {
          str.spacegroup="P4_{2}./nbc";
          str.spacegrouplabel="#133";
          str.spacegroupnumber=133;
          str.spacegroupoption="origin choice 2";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,x,z+1./2);wa(a,str);
          a.fpos=o+wv(y,-x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,-z);wa(a,str);
          a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(y+1./2,-x,-z+1./2);wa(a,str);
          a.fpos=o+wv(-y,x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,z);wa(a,str);
          a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 134  P4_{2}./nnm #134
      case 134:
        if(option==1) {
          str.spacegroup="P4_{2}./nnm";
          str.spacegrouplabel="#134";
          str.spacegroupnumber=134;
          str.spacegroupoption="origin choice 1";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(y,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-y,-x,z);wa(a,str);
          a.fpos=o+wv(y,x,z);wa(a,str);
        } else {
          str.spacegroup="P4_{2}./nnm";
          str.spacegrouplabel="#134";
          str.spacegroupnumber=134;
          str.spacegroupoption="origin choice 2";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,x,z+1./2);wa(a,str);
          a.fpos=o+wv(y,-x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,y,-z+1./2);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(y,x,-z);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,-z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(y+1./2,-x,-z+1./2);wa(a,str);
          a.fpos=o+wv(-y,x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-y,-x,z);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 135  P4_{2}./mbc #135
      case 135:
        str.spacegroup="P4_{2}./mbc";
        str.spacegrouplabel="#135";
        str.spacegroupnumber=135;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(y,-x,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
        a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(y,-x,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y,x,-z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
        a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 136  P4_{2}./mnm #136
      case 136:
        str.spacegroup="P4_{2}./mnm";
        str.spacegrouplabel="#136";
        str.spacegroupnumber=136;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(y,x,-z);wa(a,str);
        a.fpos=o+wv(-y,-x,-z);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(y+1./2,-x+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y+1./2,x+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-y,-x,z);wa(a,str);
        a.fpos=o+wv(y,x,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 137  P4_{2}./nmc #137
      case 137:
        if(option==1) {
          str.spacegroup="P4_{2}./nmc";
          str.spacegrouplabel="#137";
          str.spacegroupnumber=137;
          str.spacegroupoption="origin choice 1";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(y,x,-z);wa(a,str);
          a.fpos=o+wv(-y,-x,-z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(y,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x,-z);wa(a,str);
          a.fpos=o+wv(x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
        } else {
          str.spacegroup="P4_{2}./nmc";
          str.spacegrouplabel="#137";
          str.spacegroupnumber=137;
          str.spacegroupoption="origin choice 2";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,x,z+1./2);wa(a,str);
          a.fpos=o+wv(y,-x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y,-z);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(y+1./2,-x,-z+1./2);wa(a,str);
          a.fpos=o+wv(-y,x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(y,x,z+1./2);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 138  P4_{2}./ncm #138
      case 138:
        if(option==1) {
          str.spacegroup="P4_{2}./ncm";
          str.spacegrouplabel="#138";
          str.spacegroupnumber=138;
          str.spacegroupoption="origin choice 1";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
          a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
          a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(y,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x,-z);wa(a,str);
          a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,z);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,z);wa(a,str);
        } else {
          str.spacegroup="P4_{2}./ncm";
          str.spacegrouplabel="#138";
          str.spacegroupnumber=138;
          str.spacegroupoption="origin choice 2";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,x,z+1./2);wa(a,str);
          a.fpos=o+wv(y,-x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-y,-z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,-z);wa(a,str);
          a.fpos=o+wv(-y,-x,-z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(y+1./2,-x,-z+1./2);wa(a,str);
          a.fpos=o+wv(-y,x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,y,z+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,z);wa(a,str);
          a.fpos=o+wv(y,x,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 139  I4./mmm #139
      case 139:
        str.spacegroup="I4./mmm";
        str.spacegrouplabel="#139";
        str.spacegroupnumber=139;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-y,x,z);wa(a,str);
          a.fpos=o+wv(y,-x,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(y,x,-z);wa(a,str);
          a.fpos=o+wv(-y,-x,-z);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,y,-z);wa(a,str);
          a.fpos=o+wv(y,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x,-z);wa(a,str);
          a.fpos=o+wv(x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
          a.fpos=o+wv(-y,-x,z);wa(a,str);
          a.fpos=o+wv(y,x,z);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 140  I4./mcm #140
      case 140:
        str.spacegroup="I4./mcm";
        str.spacegrouplabel="#140";
        str.spacegroupnumber=140;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-y,x,z);wa(a,str);
          a.fpos=o+wv(y,-x,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
          a.fpos=o+wv(x,-y,-z+1./2);wa(a,str);
          a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
          a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,y,-z);wa(a,str);
          a.fpos=o+wv(y,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x,-z);wa(a,str);
          a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
          a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
          a.fpos=o+wv(y,x,z+1./2);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 141  I4_{1}./amd #141
      case 141:
        if(option==1) {
          str.spacegroup="I4_{1}./amd";
          str.spacegrouplabel="#141";
          str.spacegroupnumber=141;
          str.spacegroupoption="origin choice 1";
          for(uint j=1;j<=2;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(1./2,1./2,1./2);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x+1./2,-y+1./2,z+1./2);wa(a,str);
            a.fpos=o+wv(-y,x+1./2,z+1./4);wa(a,str);
            a.fpos=o+wv(y+1./2,-x,z+3./4);wa(a,str);
            a.fpos=o+wv(-x+1./2,y,-z+3./4);wa(a,str);
            a.fpos=o+wv(x,-y+1./2,-z+1./4);wa(a,str);
            a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
            a.fpos=o+wv(-y,-x,-z);wa(a,str);
            a.fpos=o+wv(-x,-y+1./2,-z+1./4);wa(a,str);
            a.fpos=o+wv(x+1./2,y,-z+3./4);wa(a,str);
            a.fpos=o+wv(y,-x,-z);wa(a,str);
            a.fpos=o+wv(-y+1./2,x+1./2,-z+1./2);wa(a,str);
            a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
            a.fpos=o+wv(-x,y,z);wa(a,str);
            a.fpos=o+wv(-y+1./2,-x,z+3./4);wa(a,str);
            a.fpos=o+wv(y,x+1./2,z+1./4);wa(a,str);
          }
        } else {
          str.spacegroup="I4_{1}./amd";
          str.spacegrouplabel="#141";
          str.spacegroupnumber=141;
          str.spacegroupoption="origin choice 2";
          for(uint j=1;j<=2;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(1./2,1./2,1./2);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
            a.fpos=o+wv(-y+1./4,x+3./4,z+1./4);wa(a,str);
            a.fpos=o+wv(y+1./4,-x+1./4,z+3./4);wa(a,str);
            a.fpos=o+wv(-x+1./2,y,-z+1./2);wa(a,str);
            a.fpos=o+wv(x,-y,-z);wa(a,str);
            a.fpos=o+wv(y+1./4,x+3./4,-z+1./4);wa(a,str);
            a.fpos=o+wv(-y+1./4,-x+1./4,-z+3./4);wa(a,str);
            a.fpos=o+wv(-x,-y,-z);wa(a,str);
            a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
            a.fpos=o+wv(y+3./4,-x+1./4,-z+3./4);wa(a,str);
            a.fpos=o+wv(-y+3./4,x+3./4,-z+1./4);wa(a,str);
            a.fpos=o+wv(x+1./2,-y,z+1./2);wa(a,str);
            a.fpos=o+wv(-x,y,z);wa(a,str);
            a.fpos=o+wv(-y+3./4,-x+1./4,z+3./4);wa(a,str);
            a.fpos=o+wv(y+3./4,x+3./4,z+1./4);wa(a,str);
          }
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 142  I4_{1}./acd #142
      case 142:
        if(option==1) {
          str.spacegroup="I4_{1}./acd";
          str.spacegrouplabel="#142";
          str.spacegroupnumber=142;
          str.spacegroupoption="origin choice 1";
          for(uint j=1;j<=2;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(1./2,1./2,1./2);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x+1./2,-y+1./2,z+1./2);wa(a,str);
            a.fpos=o+wv(-y,x+1./2,z+1./4);wa(a,str);
            a.fpos=o+wv(y+1./2,-x,z+3./4);wa(a,str);
            a.fpos=o+wv(-x+1./2,y,-z+1./4);wa(a,str);
            a.fpos=o+wv(x,-y+1./2,-z+3./4);wa(a,str);
            a.fpos=o+wv(y+1./2,x+1./2,-z);wa(a,str);
            a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
            a.fpos=o+wv(-x,-y+1./2,-z+1./4);wa(a,str);
            a.fpos=o+wv(x+1./2,y,-z+3./4);wa(a,str);
            a.fpos=o+wv(y,-x,-z);wa(a,str);
            a.fpos=o+wv(-y+1./2,x+1./2,-z+1./2);wa(a,str);
            a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
            a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
            a.fpos=o+wv(-y+1./2,-x,z+1./4);wa(a,str);
            a.fpos=o+wv(y,x+1./2,z+3./4);wa(a,str);
          }
        } else {
          str.spacegroup="I4_{1}./acd";
          str.spacegrouplabel="#142";
          str.spacegroupnumber=142;
          str.spacegroupoption="origin choice 2";
          for(uint j=1;j<=2;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(1./2,1./2,1./2);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
            a.fpos=o+wv(-y+1./4,x+3./4,z+1./4);wa(a,str);
            a.fpos=o+wv(y+1./4,-x+1./4,z+3./4);wa(a,str);
            a.fpos=o+wv(-x+1./2,y,-z);wa(a,str);
            a.fpos=o+wv(x,-y,-z+1./2);wa(a,str);
            a.fpos=o+wv(y+1./4,x+3./4,-z+3./4);wa(a,str);
            a.fpos=o+wv(-y+1./4,-x+1./4,-z+1./4);wa(a,str);
            a.fpos=o+wv(-x,-y,-z);wa(a,str);
            a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
            a.fpos=o+wv(y+3./4,-x+1./4,-z+3./4);wa(a,str);
            a.fpos=o+wv(-y+3./4,x+3./4,-z+1./4);wa(a,str);
            a.fpos=o+wv(x+1./2,-y,z);wa(a,str);
            a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
            a.fpos=o+wv(-y+3./4,-x+1./4,z+1./4);wa(a,str);
            a.fpos=o+wv(y+3./4,x+3./4,z+3./4);wa(a,str);
          }
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 143  P3 #143
      case 143:
        str.spacegroup="P3";
        str.spacegrouplabel="#143";
        str.spacegroupnumber=143;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 144  P3_{1} #144
      case 144:
        str.spacegroup="P3_{1}";
        str.spacegrouplabel="#144";
        str.spacegroupnumber=144;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z+1./3);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z+2./3);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 145  P3_{2} #145
      case 145:
        str.spacegroup="P3_{2}";
        str.spacegrouplabel="#145";
        str.spacegroupnumber=145;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z+2./3);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z+1./3);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 146  R3 #146
      case 146:
        if(option==1) {
          str.spacegroup="R3";
          str.spacegrouplabel="#146";
          str.spacegroupnumber=146;
          str.spacegroupoption="rhombohedral axes";
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
        } else {
          str.spacegroup="R3";
          str.spacegrouplabel="#146";
          str.spacegroupnumber=146;
          str.spacegroupoption="hexagonal axes";
          for(uint j=1;j<=3;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(2./3,1./3,1./3);
            if(j==3) o=wv(1./3,2./3,2./3);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-y,x-y,z);wa(a,str);
            a.fpos=o+wv(-x+y,-x,z);wa(a,str);
          }
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 147  P-3 #147
      case 147:
        str.spacegroup="P-3";
        str.spacegrouplabel="#147";
        str.spacegroupnumber=147;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(y,-x+y,-z);wa(a,str);
        a.fpos=o+wv(x-y,x,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 148  R-3 #148
      case 148:
        if(option==1) {
          str.spacegroup="R-3";
          str.spacegrouplabel="#148";
          str.spacegroupnumber=148;
          str.spacegroupoption="rhombohedral axes";
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(-z,-x,-y);wa(a,str);
          a.fpos=o+wv(-y,-z,-x);wa(a,str);
        } else {
          str.spacegroup="R-3";
          str.spacegrouplabel="#148";
          str.spacegroupnumber=148;
          str.spacegroupoption="hexagonal axes";
          for(uint j=1;j<=3;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(2./3,1./3,1./3);
            if(j==3) o=wv(1./3,2./3,2./3);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-y,x-y,z);wa(a,str);
            a.fpos=o+wv(-x+y,-x,z);wa(a,str);
            a.fpos=o+wv(-x,-y,-z);wa(a,str);
            a.fpos=o+wv(y,-x+y,-z);wa(a,str);
            a.fpos=o+wv(x-y,x,-z);wa(a,str);
          }
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 149  P312 #149
      case 149:
        str.spacegroup="P312";
        str.spacegrouplabel="#149";
        str.spacegroupnumber=149;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(-y,-x,-z);wa(a,str);
        a.fpos=o+wv(-x+y,y,-z);wa(a,str);
        a.fpos=o+wv(x,x-y,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 150  P321 #150
      case 150:
        str.spacegroup="P321";
        str.spacegrouplabel="#150";
        str.spacegroupnumber=150;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(y,x,-z);wa(a,str);
        a.fpos=o+wv(x-y,-y,-z);wa(a,str);
        a.fpos=o+wv(-x,-x+y,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 151  P3_{1}12 #151
      case 151:
        str.spacegroup="P3_{1}12";
        str.spacegrouplabel="#151";
        str.spacegroupnumber=151;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z+1./3);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z+2./3);wa(a,str);
        a.fpos=o+wv(-y,-x,-z+2./3);wa(a,str);
        a.fpos=o+wv(-x+y,y,-z+1./3);wa(a,str);
        a.fpos=o+wv(x,x-y,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 152  P3_{1}21 #152
      case 152:
        str.spacegroup="P3_{1}21";
        str.spacegrouplabel="#152";
        str.spacegroupnumber=152;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z+1./3);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z+2./3);wa(a,str);
        a.fpos=o+wv(y,x,-z);wa(a,str);
        a.fpos=o+wv(x-y,-y,-z+2./3);wa(a,str);
        a.fpos=o+wv(-x,-x+y,-z+1./3);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 153  P3_{2}12 #153
      case 153:
        str.spacegroup="P3_{2}12";
        str.spacegrouplabel="#153";
        str.spacegroupnumber=153;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z+2./3);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z+1./3);wa(a,str);
        a.fpos=o+wv(-y,-x,-z+1./3);wa(a,str);
        a.fpos=o+wv(-x+y,y,-z+2./3);wa(a,str);
        a.fpos=o+wv(x,x-y,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 154  P3_{2}21 #154
      case 154:
        str.spacegroup="P3_{2}21";
        str.spacegrouplabel="#154";
        str.spacegroupnumber=154;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z+2./3);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z+1./3);wa(a,str);
        a.fpos=o+wv(y,x,-z);wa(a,str);
        a.fpos=o+wv(x-y,-y,-z+1./3);wa(a,str);
        a.fpos=o+wv(-x,-x+y,-z+2./3);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 155  R32 #155
      case 155:
        if(option==1) {
          str.spacegroup="R32";
          str.spacegrouplabel="#155";
          str.spacegroupnumber=155;
          str.spacegroupoption="rhombohedral axes";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-z,-y,-x);wa(a,str);
          a.fpos=o+wv(-y,-x,-z);wa(a,str);
          a.fpos=o+wv(-x,-z,-y);wa(a,str);
        } else {
          str.spacegroup="R32";
          str.spacegrouplabel="#155";
          str.spacegroupnumber=155;
          str.spacegroupoption="hexagonal axes";
          for(uint j=1;j<=3;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(2./3,1./3,1./3);
            if(j==3) o=wv(1./3,2./3,2./3);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-y,x-y,z);wa(a,str);
            a.fpos=o+wv(-x+y,-x,z);wa(a,str);
            a.fpos=o+wv(y,x,-z);wa(a,str);
            a.fpos=o+wv(x-y,-y,-z);wa(a,str);
            a.fpos=o+wv(-x,-x+y,-z);wa(a,str);
          }
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 156  P3m1 #156
      case 156:
        str.spacegroup="P3m1";
        str.spacegrouplabel="#156";
        str.spacegroupnumber=156;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(-y,-x,z);wa(a,str);
        a.fpos=o+wv(-x+y,y,z);wa(a,str);
        a.fpos=o+wv(x,x-y,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 157  P31m #157
      case 157:
        str.spacegroup="P31m";
        str.spacegrouplabel="#157";
        str.spacegroupnumber=157;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(y,x,z);wa(a,str);
        a.fpos=o+wv(x-y,-y,z);wa(a,str);
        a.fpos=o+wv(-x,-x+y,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 158  P3c1 #158
      case 158:
        str.spacegroup="P3c1";
        str.spacegrouplabel="#158";
        str.spacegroupnumber=158;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+y,y,z+1./2);wa(a,str);
        a.fpos=o+wv(x,x-y,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 159  P31c #159
      case 159:
        str.spacegroup="P31c";
        str.spacegrouplabel="#159";
        str.spacegroupnumber=159;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(x-y,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-x+y,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 160  R3m #160
      case 160:
        if(option==1) {
          str.spacegroup="R3m";
          str.spacegrouplabel="#160";
          str.spacegroupnumber=160;
          str.spacegroupoption="rhombohedral axes";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(z,y,x);wa(a,str);
          a.fpos=o+wv(y,x,z);wa(a,str);
          a.fpos=o+wv(x,z,y);wa(a,str);
        } else {
          str.spacegroup="R3m";
          str.spacegrouplabel="#160";
          str.spacegroupnumber=160;
          str.spacegroupoption="hexagonal axes";
          for(uint j=1;j<=3;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(2./3,1./3,1./3);
            if(j==3) o=wv(1./3,2./3,2./3);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-y,x-y,z);wa(a,str);
            a.fpos=o+wv(-x+y,-x,z);wa(a,str);
            a.fpos=o+wv(-y,-x,z);wa(a,str);
            a.fpos=o+wv(-x+y,y,z);wa(a,str);
            a.fpos=o+wv(x,x-y,z);wa(a,str);
          }
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 161  R3c #161
      case 161:
        if(option==1) {
          str.spacegroup="R3c";
          str.spacegrouplabel="#161";
          str.spacegroupnumber=161;
          str.spacegroupoption="rhombohedral axes";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(z+1./2,y+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,z+1./2,y+1./2);wa(a,str);
        } else {
          str.spacegroup="R3c";
          str.spacegrouplabel="#161";
          str.spacegroupnumber=161;
          str.spacegroupoption="hexagonal axes";
          for(uint j=1;j<=3;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(2./3,1./3,1./3);
            if(j==3) o=wv(1./3,2./3,2./3);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-y,x-y,z);wa(a,str);
            a.fpos=o+wv(-x+y,-x,z);wa(a,str);
            a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
            a.fpos=o+wv(-x+y,y,z+1./2);wa(a,str);
            a.fpos=o+wv(x,x-y,z+1./2);wa(a,str);
          }
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 162  P-31m #162
      case 162:
        str.spacegroup="P-31m";
        str.spacegrouplabel="#162";
        str.spacegroupnumber=162;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(-y,-x,-z);wa(a,str);
        a.fpos=o+wv(-x+y,y,-z);wa(a,str);
        a.fpos=o+wv(x,x-y,-z);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(y,-x+y,-z);wa(a,str);
        a.fpos=o+wv(x-y,x,-z);wa(a,str);
        a.fpos=o+wv(y,x,z);wa(a,str);
        a.fpos=o+wv(x-y,-y,z);wa(a,str);
        a.fpos=o+wv(-x,-x+y,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 163  P-31c #163
      case 163:
        str.spacegroup="P-31c";
        str.spacegrouplabel="#163";
        str.spacegroupnumber=163;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x+y,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(x,x-y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(y,-x+y,-z);wa(a,str);
        a.fpos=o+wv(x-y,x,-z);wa(a,str);
        a.fpos=o+wv(y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(x-y,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-x+y,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 164  P-3m1 #164
      case 164:
        str.spacegroup="P-3m1";
        str.spacegrouplabel="#164";
        str.spacegroupnumber=164;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(y,x,-z);wa(a,str);
        a.fpos=o+wv(x-y,-y,-z);wa(a,str);
        a.fpos=o+wv(-x,-x+y,-z);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(y,-x+y,-z);wa(a,str);
        a.fpos=o+wv(x-y,x,-z);wa(a,str);
        a.fpos=o+wv(-y,-x,z);wa(a,str);
        a.fpos=o+wv(-x+y,y,z);wa(a,str);
        a.fpos=o+wv(x,x-y,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 165  P-3c1 #165
      case 165:
        str.spacegroup="P-3c1";
        str.spacegrouplabel="#165";
        str.spacegroupnumber=165;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
        a.fpos=o+wv(x-y,-y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-x+y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(y,-x+y,-z);wa(a,str);
        a.fpos=o+wv(x-y,x,-z);wa(a,str);
        a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+y,y,z+1./2);wa(a,str);
        a.fpos=o+wv(x,x-y,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 166  R-3m #166
      case 166:
        if(option==1) {
          str.spacegroup="R-3m";
          str.spacegrouplabel="#166";
          str.spacegroupnumber=166;
          str.spacegroupoption="rhombohedral axes";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-z,-y,-x);wa(a,str);
          a.fpos=o+wv(-y,-x,-z);wa(a,str);
          a.fpos=o+wv(-x,-z,-y);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(-z,-x,-y);wa(a,str);
          a.fpos=o+wv(-y,-z,-x);wa(a,str);
          a.fpos=o+wv(z,y,x);wa(a,str);
          a.fpos=o+wv(y,x,z);wa(a,str);
          a.fpos=o+wv(x,z,y);wa(a,str);
        } else {
          str.spacegroup="R-3m";
          str.spacegrouplabel="#166";
          str.spacegroupnumber=166;
          str.spacegroupoption="hexagonal axes";
          for(uint j=1;j<=3;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(2./3,1./3,1./3);
            if(j==3) o=wv(1./3,2./3,2./3);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-y,x-y,z);wa(a,str);
            a.fpos=o+wv(-x+y,-x,z);wa(a,str);
            a.fpos=o+wv(y,x,-z);wa(a,str);
            a.fpos=o+wv(x-y,-y,-z);wa(a,str);
            a.fpos=o+wv(-x,-x+y,-z);wa(a,str);
            a.fpos=o+wv(-x,-y,-z);wa(a,str);
            a.fpos=o+wv(y,-x+y,-z);wa(a,str);
            a.fpos=o+wv(x-y,x,-z);wa(a,str);
            a.fpos=o+wv(-y,-x,z);wa(a,str);
            a.fpos=o+wv(-x+y,y,z);wa(a,str);
            a.fpos=o+wv(x,x-y,z);wa(a,str);
          }
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 167  R-3c #167
      case 167:
        if(option==1) {
          str.spacegroup="R-3c";
          str.spacegrouplabel="#167";
          str.spacegroupnumber=167;
          str.spacegroupoption="rhombohedral axes";
          o=wv(0,0,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-z+1./2,-y+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,-z+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(-z,-x,-y);wa(a,str);
          a.fpos=o+wv(-y,-z,-x);wa(a,str);
          a.fpos=o+wv(z+1./2,y+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,z+1./2,y+1./2);wa(a,str);
        } else {
          str.spacegroup="R-3c";
          str.spacegrouplabel="#167";
          str.spacegroupnumber=167;
          str.spacegroupoption="hexagonal axes";
          for(uint j=1;j<=3;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(2./3,1./3,1./3);
            if(j==3) o=wv(1./3,2./3,2./3);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-y,x-y,z);wa(a,str);
            a.fpos=o+wv(-x+y,-x,z);wa(a,str);
            a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
            a.fpos=o+wv(x-y,-y,-z+1./2);wa(a,str);
            a.fpos=o+wv(-x,-x+y,-z+1./2);wa(a,str);
            a.fpos=o+wv(-x,-y,-z);wa(a,str);
            a.fpos=o+wv(y,-x+y,-z);wa(a,str);
            a.fpos=o+wv(x-y,x,-z);wa(a,str);
            a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
            a.fpos=o+wv(-x+y,y,z+1./2);wa(a,str);
            a.fpos=o+wv(x,x-y,z+1./2);wa(a,str);
          }
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 168  P6 #168
      case 168:
        str.spacegroup="P6";
        str.spacegrouplabel="#168";
        str.spacegroupnumber=168;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(y,-x+y,z);wa(a,str);
        a.fpos=o+wv(x-y,x,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 169  P6_{1} #169
      case 169:
        str.spacegroup="P6_{1}";
        str.spacegrouplabel="#169";
        str.spacegroupnumber=169;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z+1./3);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z+2./3);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(y,-x+y,z+5./6);wa(a,str);
        a.fpos=o+wv(x-y,x,z+1./6);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 170  P6_{5} #170
      case 170:
        str.spacegroup="P6_{5}";
        str.spacegrouplabel="#170";
        str.spacegroupnumber=170;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z+2./3);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z+1./3);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(y,-x+y,z+1./6);wa(a,str);
        a.fpos=o+wv(x-y,x,z+5./6);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 171  P6_{2} #171
      case 171:
        str.spacegroup="P6_{2}";
        str.spacegrouplabel="#171";
        str.spacegroupnumber=171;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z+2./3);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z+1./3);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(y,-x+y,z+2./3);wa(a,str);
        a.fpos=o+wv(x-y,x,z+1./3);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 172  P6_{4} #172
      case 172:
        str.spacegroup="P6_{4}";
        str.spacegrouplabel="#172";
        str.spacegroupnumber=172;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z+1./3);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z+2./3);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(y,-x+y,z+1./3);wa(a,str);
        a.fpos=o+wv(x-y,x,z+2./3);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 173  P6_{3} #173
      case 173:
        str.spacegroup="P6_{3}";
        str.spacegrouplabel="#173";
        str.spacegroupnumber=173;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(y,-x+y,z+1./2);wa(a,str);
        a.fpos=o+wv(x-y,x,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 174  P-6 #174
      case 174:
        str.spacegroup="P-6";
        str.spacegrouplabel="#174";
        str.spacegroupnumber=174;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(-y,x-y,-z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 175  P6./m #175
      case 175:
        str.spacegroup="P6./m";
        str.spacegrouplabel="#175";
        str.spacegroupnumber=175;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(y,-x+y,z);wa(a,str);
        a.fpos=o+wv(x-y,x,z);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(y,-x+y,-z);wa(a,str);
        a.fpos=o+wv(x-y,x,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(-y,x-y,-z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 176  P6_{3}./m #176
      case 176:
        str.spacegroup="P6_{3}./m";
        str.spacegrouplabel="#176";
        str.spacegroupnumber=176;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(y,-x+y,z+1./2);wa(a,str);
        a.fpos=o+wv(x-y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(y,-x+y,-z);wa(a,str);
        a.fpos=o+wv(x-y,x,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y,x-y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x+y,-x,-z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 177  P622 #177
      case 177:
        str.spacegroup="P622";
        str.spacegrouplabel="#177";
        str.spacegroupnumber=177;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(y,-x+y,z);wa(a,str);
        a.fpos=o+wv(x-y,x,z);wa(a,str);
        a.fpos=o+wv(y,x,-z);wa(a,str);
        a.fpos=o+wv(x-y,-y,-z);wa(a,str);
        a.fpos=o+wv(-x,-x+y,-z);wa(a,str);
        a.fpos=o+wv(-y,-x,-z);wa(a,str);
        a.fpos=o+wv(-x+y,y,-z);wa(a,str);
        a.fpos=o+wv(x,x-y,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 178  P6_{1}22 #178
      case 178:
        str.spacegroup="P6_{1}22";
        str.spacegrouplabel="#178";
        str.spacegroupnumber=178;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z+1./3);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z+2./3);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(y,-x+y,z+5./6);wa(a,str);
        a.fpos=o+wv(x-y,x,z+1./6);wa(a,str);
        a.fpos=o+wv(y,x,-z+1./3);wa(a,str);
        a.fpos=o+wv(x-y,-y,-z);wa(a,str);
        a.fpos=o+wv(-x,-x+y,-z+2./3);wa(a,str);
        a.fpos=o+wv(-y,-x,-z+5./6);wa(a,str);
        a.fpos=o+wv(-x+y,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(x,x-y,-z+1./6);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 179  P6_{5}22 #179
      case 179:
        str.spacegroup="P6_{5}22";
        str.spacegrouplabel="#179";
        str.spacegroupnumber=179;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z+2./3);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z+1./3);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(y,-x+y,z+1./6);wa(a,str);
        a.fpos=o+wv(x-y,x,z+5./6);wa(a,str);
        a.fpos=o+wv(y,x,-z+2./3);wa(a,str);
        a.fpos=o+wv(x-y,-y,-z);wa(a,str);
        a.fpos=o+wv(-x,-x+y,-z+1./3);wa(a,str);
        a.fpos=o+wv(-y,-x,-z+1./6);wa(a,str);
        a.fpos=o+wv(-x+y,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(x,x-y,-z+5./6);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 180  P6_{2}22 #180
      case 180:
        str.spacegroup="P6_{2}22";
        str.spacegrouplabel="#180";
        str.spacegroupnumber=180;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z+2./3);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z+1./3);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(y,-x+y,z+2./3);wa(a,str);
        a.fpos=o+wv(x-y,x,z+1./3);wa(a,str);
        a.fpos=o+wv(y,x,-z+2./3);wa(a,str);
        a.fpos=o+wv(x-y,-y,-z);wa(a,str);
        a.fpos=o+wv(-x,-x+y,-z+1./3);wa(a,str);
        a.fpos=o+wv(-y,-x,-z+2./3);wa(a,str);
        a.fpos=o+wv(-x+y,y,-z);wa(a,str);
        a.fpos=o+wv(x,x-y,-z+1./3);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 181  P6_{4}22 #181
      case 181:
        str.spacegroup="P6_{4}22";
        str.spacegrouplabel="#181";
        str.spacegroupnumber=181;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z+1./3);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z+2./3);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(y,-x+y,z+1./3);wa(a,str);
        a.fpos=o+wv(x-y,x,z+2./3);wa(a,str);
        a.fpos=o+wv(y,x,-z+1./3);wa(a,str);
        a.fpos=o+wv(x-y,-y,-z);wa(a,str);
        a.fpos=o+wv(-x,-x+y,-z+2./3);wa(a,str);
        a.fpos=o+wv(-y,-x,-z+1./3);wa(a,str);
        a.fpos=o+wv(-x+y,y,-z);wa(a,str);
        a.fpos=o+wv(x,x-y,-z+2./3);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 182  P6_{3}22 #182
      case 182:
        str.spacegroup="P6_{3}22";
        str.spacegrouplabel="#182";
        str.spacegroupnumber=182;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(y,-x+y,z+1./2);wa(a,str);
        a.fpos=o+wv(x-y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(y,x,-z);wa(a,str);
        a.fpos=o+wv(x-y,-y,-z);wa(a,str);
        a.fpos=o+wv(-x,-x+y,-z);wa(a,str);
        a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x+y,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(x,x-y,-z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 183  P6mm #183
      case 183:
        str.spacegroup="P6mm";
        str.spacegrouplabel="#183";
        str.spacegroupnumber=183;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(y,-x+y,z);wa(a,str);
        a.fpos=o+wv(x-y,x,z);wa(a,str);
        a.fpos=o+wv(-y,-x,z);wa(a,str);
        a.fpos=o+wv(-x+y,y,z);wa(a,str);
        a.fpos=o+wv(x,x-y,z);wa(a,str);
        a.fpos=o+wv(y,x,z);wa(a,str);
        a.fpos=o+wv(x-y,-y,z);wa(a,str);
        a.fpos=o+wv(-x,-x+y,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 184  P6cc #184
      case 184:
        str.spacegroup="P6cc";
        str.spacegrouplabel="#184";
        str.spacegroupnumber=184;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(y,-x+y,z);wa(a,str);
        a.fpos=o+wv(x-y,x,z);wa(a,str);
        a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+y,y,z+1./2);wa(a,str);
        a.fpos=o+wv(x,x-y,z+1./2);wa(a,str);
        a.fpos=o+wv(y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(x-y,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-x+y,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 185  P6_{3}cm #185
      case 185:
        str.spacegroup="P6_{3}cm";
        str.spacegrouplabel="#185";
        str.spacegroupnumber=185;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(y,-x+y,z+1./2);wa(a,str);
        a.fpos=o+wv(x-y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+y,y,z+1./2);wa(a,str);
        a.fpos=o+wv(x,x-y,z+1./2);wa(a,str);
        a.fpos=o+wv(y,x,z);wa(a,str);
        a.fpos=o+wv(x-y,-y,z);wa(a,str);
        a.fpos=o+wv(-x,-x+y,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 186  P6_{3}mc #186
      case 186:
        str.spacegroup="P6_{3}mc";
        str.spacegrouplabel="#186";
        str.spacegroupnumber=186;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(y,-x+y,z+1./2);wa(a,str);
        a.fpos=o+wv(x-y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(-y,-x,z);wa(a,str);
        a.fpos=o+wv(-x+y,y,z);wa(a,str);
        a.fpos=o+wv(x,x-y,z);wa(a,str);
        a.fpos=o+wv(y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(x-y,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-x+y,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 187  P-6m2 #187
      case 187:
        str.spacegroup="P-6m2";
        str.spacegrouplabel="#187";
        str.spacegroupnumber=187;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(-y,x-y,-z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,-z);wa(a,str);
        a.fpos=o+wv(-y,-x,z);wa(a,str);
        a.fpos=o+wv(-x+y,y,z);wa(a,str);
        a.fpos=o+wv(x,x-y,z);wa(a,str);
        a.fpos=o+wv(-y,-x,-z);wa(a,str);
        a.fpos=o+wv(-x+y,y,-z);wa(a,str);
        a.fpos=o+wv(x,x-y,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 188  P-6c2 #188
      case 188:
        str.spacegroup="P-6c2";
        str.spacegrouplabel="#188";
        str.spacegroupnumber=188;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(x,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y,x-y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x+y,-x,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+y,y,z+1./2);wa(a,str);
        a.fpos=o+wv(x,x-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-y,-x,-z);wa(a,str);
        a.fpos=o+wv(-x+y,y,-z);wa(a,str);
        a.fpos=o+wv(x,x-y,-z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 189  P-62m #189
      case 189:
        str.spacegroup="P-62m";
        str.spacegrouplabel="#189";
        str.spacegroupnumber=189;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(-y,x-y,-z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,-z);wa(a,str);
        a.fpos=o+wv(y,x,-z);wa(a,str);
        a.fpos=o+wv(x-y,-y,-z);wa(a,str);
        a.fpos=o+wv(-x,-x+y,-z);wa(a,str);
        a.fpos=o+wv(y,x,z);wa(a,str);
        a.fpos=o+wv(x-y,-y,z);wa(a,str);
        a.fpos=o+wv(-x,-x+y,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 190  P-62c #190
      case 190:
        str.spacegroup="P-62c";
        str.spacegrouplabel="#190";
        str.spacegroupnumber=190;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(x,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y,x-y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x+y,-x,-z+1./2);wa(a,str);
        a.fpos=o+wv(y,x,-z);wa(a,str);
        a.fpos=o+wv(x-y,-y,-z);wa(a,str);
        a.fpos=o+wv(-x,-x+y,-z);wa(a,str);
        a.fpos=o+wv(y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(x-y,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-x+y,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 191  P6./mmm #191
      case 191:
        str.spacegroup="P6./mmm";
        str.spacegrouplabel="#191";
        str.spacegroupnumber=191;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(y,-x+y,z);wa(a,str);
        a.fpos=o+wv(x-y,x,z);wa(a,str);
        a.fpos=o+wv(y,x,-z);wa(a,str);
        a.fpos=o+wv(x-y,-y,-z);wa(a,str);
        a.fpos=o+wv(-x,-x+y,-z);wa(a,str);
        a.fpos=o+wv(-y,-x,-z);wa(a,str);
        a.fpos=o+wv(-x+y,y,-z);wa(a,str);
        a.fpos=o+wv(x,x-y,-z);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(y,-x+y,-z);wa(a,str);
        a.fpos=o+wv(x-y,x,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(-y,x-y,-z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,-z);wa(a,str);
        a.fpos=o+wv(-y,-x,z);wa(a,str);
        a.fpos=o+wv(-x+y,y,z);wa(a,str);
        a.fpos=o+wv(x,x-y,z);wa(a,str);
        a.fpos=o+wv(y,x,z);wa(a,str);
        a.fpos=o+wv(x-y,-y,z);wa(a,str);
        a.fpos=o+wv(-x,-x+y,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 192  P6./mcc #192
      case 192:
        str.spacegroup="P6./mcc";
        str.spacegrouplabel="#192";
        str.spacegroupnumber=192;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(y,-x+y,z);wa(a,str);
        a.fpos=o+wv(x-y,x,z);wa(a,str);
        a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
        a.fpos=o+wv(x-y,-y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-x+y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x+y,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(x,x-y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(y,-x+y,-z);wa(a,str);
        a.fpos=o+wv(x-y,x,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(-y,x-y,-z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,-z);wa(a,str);
        a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+y,y,z+1./2);wa(a,str);
        a.fpos=o+wv(x,x-y,z+1./2);wa(a,str);
        a.fpos=o+wv(y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(x-y,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-x+y,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 193  P6_{3}./mcm #193
      case 193:
        str.spacegroup="P6_{3}./mcm";
        str.spacegrouplabel="#193";
        str.spacegroupnumber=193;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(y,-x+y,z+1./2);wa(a,str);
        a.fpos=o+wv(x-y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
        a.fpos=o+wv(x-y,-y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-x+y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y,-x,-z);wa(a,str);
        a.fpos=o+wv(-x+y,y,-z);wa(a,str);
        a.fpos=o+wv(x,x-y,-z);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(y,-x+y,-z);wa(a,str);
        a.fpos=o+wv(x-y,x,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y,x-y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x+y,-x,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+y,y,z+1./2);wa(a,str);
        a.fpos=o+wv(x,x-y,z+1./2);wa(a,str);
        a.fpos=o+wv(y,x,z);wa(a,str);
        a.fpos=o+wv(x-y,-y,z);wa(a,str);
        a.fpos=o+wv(-x,-x+y,z);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 194  P6_{3}./mmc #194
      case 194:
        str.spacegroup="P6_{3}./mmc";
        str.spacegrouplabel="#194";
        str.spacegroupnumber=194;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-y,x-y,z);wa(a,str);
        a.fpos=o+wv(-x+y,-x,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(y,-x+y,z+1./2);wa(a,str);
        a.fpos=o+wv(x-y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(y,x,-z);wa(a,str);
        a.fpos=o+wv(x-y,-y,-z);wa(a,str);
        a.fpos=o+wv(-x,-x+y,-z);wa(a,str);
        a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x+y,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(x,x-y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(y,-x+y,-z);wa(a,str);
        a.fpos=o+wv(x-y,x,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y,x-y,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x+y,-x,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y,-x,z);wa(a,str);
        a.fpos=o+wv(-x+y,y,z);wa(a,str);
        a.fpos=o+wv(x,x-y,z);wa(a,str);
        a.fpos=o+wv(y,x,z+1./2);wa(a,str);
        a.fpos=o+wv(x-y,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,-x+y,z+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 195  P23 #195
      case 195:
        str.spacegroup="P23";
        str.spacegrouplabel="#195";
        str.spacegroupnumber=195;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,-z);wa(a,str);
        a.fpos=o+wv(z,x,y);wa(a,str);
        a.fpos=o+wv(z,-x,-y);wa(a,str);
        a.fpos=o+wv(-z,-x,y);wa(a,str);
        a.fpos=o+wv(-z,x,-y);wa(a,str);
        a.fpos=o+wv(y,z,x);wa(a,str);
        a.fpos=o+wv(-y,z,-x);wa(a,str);
        a.fpos=o+wv(y,-z,-x);wa(a,str);
        a.fpos=o+wv(-y,-z,x);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 196  F23 #196
      case 196:
        str.spacegroup="F23";
        str.spacegrouplabel="#196";
        str.spacegroupnumber=196;
        str.spacegroupoption="";
        for(uint j=1;j<=4;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(0,1./2,1./2);
          if(j==3) o=wv(1./2,0,1./2);
          if(j==4) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,-x,y);wa(a,str);
          a.fpos=o+wv(-z,x,-y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y,z,-x);wa(a,str);
          a.fpos=o+wv(y,-z,-x);wa(a,str);
          a.fpos=o+wv(-y,-z,x);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 197  I23 #197
      case 197:
        str.spacegroup="I23";
        str.spacegrouplabel="#197";
        str.spacegroupnumber=197;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,-x,y);wa(a,str);
          a.fpos=o+wv(-z,x,-y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y,z,-x);wa(a,str);
          a.fpos=o+wv(y,-z,-x);wa(a,str);
          a.fpos=o+wv(-y,-z,x);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 198  P2_{1}3 #198
      case 198:
        str.spacegroup="P2_{1}3";
        str.spacegrouplabel="#198";
        str.spacegroupnumber=198;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
        a.fpos=o+wv(z,x,y);wa(a,str);
        a.fpos=o+wv(z+1./2,-x+1./2,-y);wa(a,str);
        a.fpos=o+wv(-z+1./2,-x,y+1./2);wa(a,str);
        a.fpos=o+wv(-z,x+1./2,-y+1./2);wa(a,str);
        a.fpos=o+wv(y,z,x);wa(a,str);
        a.fpos=o+wv(-y,z+1./2,-x+1./2);wa(a,str);
        a.fpos=o+wv(y+1./2,-z+1./2,-x);wa(a,str);
        a.fpos=o+wv(-y+1./2,-z,x+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 199  I2_{1}3 #199
      case 199:
        str.spacegroup="I2_{1}3";
        str.spacegrouplabel="#199";
        str.spacegroupnumber=199;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z+1./2,-x+1./2,-y);wa(a,str);
          a.fpos=o+wv(-z+1./2,-x,y+1./2);wa(a,str);
          a.fpos=o+wv(-z,x+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y,z+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,-z+1./2,-x);wa(a,str);
          a.fpos=o+wv(-y+1./2,-z,x+1./2);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 200  Pm-3 #200
      case 200:
        str.spacegroup="Pm-3";
        str.spacegrouplabel="#200";
        str.spacegroupnumber=200;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,-z);wa(a,str);
        a.fpos=o+wv(z,x,y);wa(a,str);
        a.fpos=o+wv(z,-x,-y);wa(a,str);
        a.fpos=o+wv(-z,-x,y);wa(a,str);
        a.fpos=o+wv(-z,x,-y);wa(a,str);
        a.fpos=o+wv(y,z,x);wa(a,str);
        a.fpos=o+wv(-y,z,-x);wa(a,str);
        a.fpos=o+wv(y,-z,-x);wa(a,str);
        a.fpos=o+wv(-y,-z,x);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,z);wa(a,str);
        a.fpos=o+wv(-z,-x,-y);wa(a,str);
        a.fpos=o+wv(-z,x,y);wa(a,str);
        a.fpos=o+wv(z,x,-y);wa(a,str);
        a.fpos=o+wv(z,-x,y);wa(a,str);
        a.fpos=o+wv(-y,-z,-x);wa(a,str);
        a.fpos=o+wv(y,-z,x);wa(a,str);
        a.fpos=o+wv(-y,z,x);wa(a,str);
        a.fpos=o+wv(y,z,-x);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 201  Pn-3 #201
      case 201:
        if(option==1) {
          str.spacegroup="Pn-3";
          str.spacegrouplabel="#201";
          str.spacegroupnumber=201;
          str.spacegroupoption="origin choice 1";
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,-x,y);wa(a,str);
          a.fpos=o+wv(-z,x,-y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y,z,-x);wa(a,str);
          a.fpos=o+wv(y,-z,-x);wa(a,str);
          a.fpos=o+wv(-y,-z,x);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-z+1./2,-x+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(-z+1./2,x+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(z+1./2,x+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(z+1./2,-x+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,-z+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,-z+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,z+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,z+1./2,-x+1./2);wa(a,str);
        } else {
          str.spacegroup="Pn-3";
          str.spacegrouplabel="#201";
          str.spacegroupnumber=201;
          str.spacegroupoption="origin choice 2";
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y,-z+1./2);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z,-x+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(-z+1./2,-x+1./2,y);wa(a,str);
          a.fpos=o+wv(-z+1./2,x,-y+1./2);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y+1./2,z,-x+1./2);wa(a,str);
          a.fpos=o+wv(y,-z+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,-z+1./2,x);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,x+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(z+1./2,x+1./2,-y);wa(a,str);
          a.fpos=o+wv(z+1./2,-x,y+1./2);wa(a,str);
          a.fpos=o+wv(-y,-z,-x);wa(a,str);
          a.fpos=o+wv(y+1./2,-z,x+1./2);wa(a,str);
          a.fpos=o+wv(-y,z+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,z+1./2,-x);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 202  Fm-3 #202
      case 202:
        str.spacegroup="Fm-3";
        str.spacegrouplabel="#202";
        str.spacegroupnumber=202;
        str.spacegroupoption="";
        for(uint j=1;j<=4;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(0,1./2,1./2);
          if(j==3) o=wv(1./2,0,1./2);
          if(j==4) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,-x,y);wa(a,str);
          a.fpos=o+wv(-z,x,-y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y,z,-x);wa(a,str);
          a.fpos=o+wv(y,-z,-x);wa(a,str);
          a.fpos=o+wv(-y,-z,x);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
          a.fpos=o+wv(-z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,x,y);wa(a,str);
          a.fpos=o+wv(z,x,-y);wa(a,str);
          a.fpos=o+wv(z,-x,y);wa(a,str);
          a.fpos=o+wv(-y,-z,-x);wa(a,str);
          a.fpos=o+wv(y,-z,x);wa(a,str);
          a.fpos=o+wv(-y,z,x);wa(a,str);
          a.fpos=o+wv(y,z,-x);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 203  Fd-3 #203
      case 203:
        if(option==1) {
          str.spacegroup="Fd-3";
          str.spacegrouplabel="#203";
          str.spacegroupnumber=203;
          str.spacegroupoption="origin choice 1";
          for(uint j=1;j<=4;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(0,1./2,1./2);
            if(j==3) o=wv(1./2,0,1./2);
            if(j==4) o=wv(1./2,1./2,0);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x,-y,z);wa(a,str);
            a.fpos=o+wv(-x,y,-z);wa(a,str);
            a.fpos=o+wv(x,-y,-z);wa(a,str);
            a.fpos=o+wv(z,x,y);wa(a,str);
            a.fpos=o+wv(z,-x,-y);wa(a,str);
            a.fpos=o+wv(-z,-x,y);wa(a,str);
            a.fpos=o+wv(-z,x,-y);wa(a,str);
            a.fpos=o+wv(y,z,x);wa(a,str);
            a.fpos=o+wv(-y,z,-x);wa(a,str);
            a.fpos=o+wv(y,-z,-x);wa(a,str);
            a.fpos=o+wv(-y,-z,x);wa(a,str);
            a.fpos=o+wv(-x+1./4,-y+1./4,-z+1./4);wa(a,str);
            a.fpos=o+wv(x+1./4,y+1./4,-z+1./4);wa(a,str);
            a.fpos=o+wv(x+1./4,-y+1./4,z+1./4);wa(a,str);
            a.fpos=o+wv(-x+1./4,y+1./4,z+1./4);wa(a,str);
            a.fpos=o+wv(-z+1./4,-x+1./4,-y+1./4);wa(a,str);
            a.fpos=o+wv(-z+1./4,x+1./4,y+1./4);wa(a,str);
            a.fpos=o+wv(z+1./4,x+1./4,-y+1./4);wa(a,str);
            a.fpos=o+wv(z+1./4,-x+1./4,y+1./4);wa(a,str);
            a.fpos=o+wv(-y+1./4,-z+1./4,-x+1./4);wa(a,str);
            a.fpos=o+wv(y+1./4,-z+1./4,x+1./4);wa(a,str);
            a.fpos=o+wv(-y+1./4,z+1./4,x+1./4);wa(a,str);
            a.fpos=o+wv(y+1./4,z+1./4,-x+1./4);wa(a,str);
          }
        } else {
          str.spacegroup="Fd-3";
          str.spacegrouplabel="#203";
          str.spacegroupnumber=203;
          str.spacegroupoption="origin choice 2";
          for(uint j=1;j<=4;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(0,1./2,1./2);
            if(j==3) o=wv(1./2,0,1./2);
            if(j==4) o=wv(1./2,1./2,0);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x+3./4,-y+3./4,z);wa(a,str);
            a.fpos=o+wv(-x+3./4,y,-z+3./4);wa(a,str);
            a.fpos=o+wv(x,-y+3./4,-z+3./4);wa(a,str);
            a.fpos=o+wv(z,x,y);wa(a,str);
            a.fpos=o+wv(z,-x+3./4,-y+3./4);wa(a,str);
            a.fpos=o+wv(-z+3./4,-x+3./4,y);wa(a,str);
            a.fpos=o+wv(-z+3./4,x,-y+3./4);wa(a,str);
            a.fpos=o+wv(y,z,x);wa(a,str);
            a.fpos=o+wv(-y+3./4,z,-x+3./4);wa(a,str);
            a.fpos=o+wv(y,-z+3./4,-x+3./4);wa(a,str);
            a.fpos=o+wv(-y+3./4,-z+3./4,x);wa(a,str);
            a.fpos=o+wv(-x,-y,-z);wa(a,str);
            a.fpos=o+wv(x+1./4,y+1./4,-z);wa(a,str);
            a.fpos=o+wv(x+1./4,-y,z+1./4);wa(a,str);
            a.fpos=o+wv(-x,y+1./4,z+1./4);wa(a,str);
            a.fpos=o+wv(-z,-x,-y);wa(a,str);
            a.fpos=o+wv(-z,x+1./4,y+1./4);wa(a,str);
            a.fpos=o+wv(z+1./4,x+1./4,-y);wa(a,str);
            a.fpos=o+wv(z+1./4,-x,y+1./4);wa(a,str);
            a.fpos=o+wv(-y,-z,-x);wa(a,str);
            a.fpos=o+wv(y+1./4,-z,x+1./4);wa(a,str);
            a.fpos=o+wv(-y,z+1./4,x+1./4);wa(a,str);
            a.fpos=o+wv(y+1./4,z+1./4,-x);wa(a,str);
          }
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 204  Im-3 #204
      case 204:
        str.spacegroup="Im-3";
        str.spacegrouplabel="#204";
        str.spacegroupnumber=204;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,-x,y);wa(a,str);
          a.fpos=o+wv(-z,x,-y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y,z,-x);wa(a,str);
          a.fpos=o+wv(y,-z,-x);wa(a,str);
          a.fpos=o+wv(-y,-z,x);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
          a.fpos=o+wv(-z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,x,y);wa(a,str);
          a.fpos=o+wv(z,x,-y);wa(a,str);
          a.fpos=o+wv(z,-x,y);wa(a,str);
          a.fpos=o+wv(-y,-z,-x);wa(a,str);
          a.fpos=o+wv(y,-z,x);wa(a,str);
          a.fpos=o+wv(-y,z,x);wa(a,str);
          a.fpos=o+wv(y,z,-x);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 205  Pa-3 #205
      case 205:
        str.spacegroup="Pa-3";
        str.spacegrouplabel="#205";
        str.spacegroupnumber=205;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
        a.fpos=o+wv(z,x,y);wa(a,str);
        a.fpos=o+wv(z+1./2,-x+1./2,-y);wa(a,str);
        a.fpos=o+wv(-z+1./2,-x,y+1./2);wa(a,str);
        a.fpos=o+wv(-z,x+1./2,-y+1./2);wa(a,str);
        a.fpos=o+wv(y,z,x);wa(a,str);
        a.fpos=o+wv(-y,z+1./2,-x+1./2);wa(a,str);
        a.fpos=o+wv(y+1./2,-z+1./2,-x);wa(a,str);
        a.fpos=o+wv(-y+1./2,-z,x+1./2);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
        a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
        a.fpos=o+wv(-z,-x,-y);wa(a,str);
        a.fpos=o+wv(-z+1./2,x+1./2,y);wa(a,str);
        a.fpos=o+wv(z+1./2,x,-y+1./2);wa(a,str);
        a.fpos=o+wv(z,-x+1./2,y+1./2);wa(a,str);
        a.fpos=o+wv(-y,-z,-x);wa(a,str);
        a.fpos=o+wv(y,-z+1./2,x+1./2);wa(a,str);
        a.fpos=o+wv(-y+1./2,z+1./2,x);wa(a,str);
        a.fpos=o+wv(y+1./2,z,-x+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 206  Ia-3 #206
      case 206:
        str.spacegroup="Ia-3";
        str.spacegrouplabel="#206";
        str.spacegroupnumber=206;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z+1./2,-x+1./2,-y);wa(a,str);
          a.fpos=o+wv(-z+1./2,-x,y+1./2);wa(a,str);
          a.fpos=o+wv(-z,x+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y,z+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,-z+1./2,-x);wa(a,str);
          a.fpos=o+wv(-y+1./2,-z,x+1./2);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
          a.fpos=o+wv(-z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z+1./2,x+1./2,y);wa(a,str);
          a.fpos=o+wv(z+1./2,x,-y+1./2);wa(a,str);
          a.fpos=o+wv(z,-x+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(-y,-z,-x);wa(a,str);
          a.fpos=o+wv(y,-z+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,z+1./2,x);wa(a,str);
          a.fpos=o+wv(y+1./2,z,-x+1./2);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 207  P432 #207
      case 207:
        str.spacegroup="P432";
        str.spacegrouplabel="#207";
        str.spacegroupnumber=207;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,-z);wa(a,str);
        a.fpos=o+wv(z,x,y);wa(a,str);
        a.fpos=o+wv(z,-x,-y);wa(a,str);
        a.fpos=o+wv(-z,-x,y);wa(a,str);
        a.fpos=o+wv(-z,x,-y);wa(a,str);
        a.fpos=o+wv(y,z,x);wa(a,str);
        a.fpos=o+wv(-y,z,-x);wa(a,str);
        a.fpos=o+wv(y,-z,-x);wa(a,str);
        a.fpos=o+wv(-y,-z,x);wa(a,str);
        a.fpos=o+wv(y,x,-z);wa(a,str);
        a.fpos=o+wv(-y,-x,-z);wa(a,str);
        a.fpos=o+wv(y,-x,z);wa(a,str);
        a.fpos=o+wv(-y,x,z);wa(a,str);
        a.fpos=o+wv(x,z,-y);wa(a,str);
        a.fpos=o+wv(-x,z,y);wa(a,str);
        a.fpos=o+wv(-x,-z,-y);wa(a,str);
        a.fpos=o+wv(x,-z,y);wa(a,str);
        a.fpos=o+wv(z,y,-x);wa(a,str);
        a.fpos=o+wv(z,-y,x);wa(a,str);
        a.fpos=o+wv(-z,y,x);wa(a,str);
        a.fpos=o+wv(-z,-y,-x);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 208  P4_{2}32 #208
      case 208:
        str.spacegroup="P4_{2}32";
        str.spacegrouplabel="#208";
        str.spacegroupnumber=208;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,-z);wa(a,str);
        a.fpos=o+wv(z,x,y);wa(a,str);
        a.fpos=o+wv(z,-x,-y);wa(a,str);
        a.fpos=o+wv(-z,-x,y);wa(a,str);
        a.fpos=o+wv(-z,x,-y);wa(a,str);
        a.fpos=o+wv(y,z,x);wa(a,str);
        a.fpos=o+wv(-y,z,-x);wa(a,str);
        a.fpos=o+wv(y,-z,-x);wa(a,str);
        a.fpos=o+wv(-y,-z,x);wa(a,str);
        a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,z+1./2,-y+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,z+1./2,y+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,-z+1./2,-y+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-z+1./2,y+1./2);wa(a,str);
        a.fpos=o+wv(z+1./2,y+1./2,-x+1./2);wa(a,str);
        a.fpos=o+wv(z+1./2,-y+1./2,x+1./2);wa(a,str);
        a.fpos=o+wv(-z+1./2,y+1./2,x+1./2);wa(a,str);
        a.fpos=o+wv(-z+1./2,-y+1./2,-x+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 209  F432 #209
      case 209:
        str.spacegroup="F432";
        str.spacegrouplabel="#209";
        str.spacegroupnumber=209;
        str.spacegroupoption="";
        for(uint j=1;j<=4;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(0,1./2,1./2);
          if(j==3) o=wv(1./2,0,1./2);
          if(j==4) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,-x,y);wa(a,str);
          a.fpos=o+wv(-z,x,-y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y,z,-x);wa(a,str);
          a.fpos=o+wv(y,-z,-x);wa(a,str);
          a.fpos=o+wv(-y,-z,x);wa(a,str);
          a.fpos=o+wv(y,x,-z);wa(a,str);
          a.fpos=o+wv(-y,-x,-z);wa(a,str);
          a.fpos=o+wv(y,-x,z);wa(a,str);
          a.fpos=o+wv(-y,x,z);wa(a,str);
          a.fpos=o+wv(x,z,-y);wa(a,str);
          a.fpos=o+wv(-x,z,y);wa(a,str);
          a.fpos=o+wv(-x,-z,-y);wa(a,str);
          a.fpos=o+wv(x,-z,y);wa(a,str);
          a.fpos=o+wv(z,y,-x);wa(a,str);
          a.fpos=o+wv(z,-y,x);wa(a,str);
          a.fpos=o+wv(-z,y,x);wa(a,str);
          a.fpos=o+wv(-z,-y,-x);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 210  F4_{1}32 #210
      case 210:
        str.spacegroup="F4_{1}32";
        str.spacegrouplabel="#210";
        str.spacegroupnumber=210;
        str.spacegroupoption="";
        for(uint j=1;j<=4;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(0,1./2,1./2);
          if(j==3) o=wv(1./2,0,1./2);
          if(j==4) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y,-z+1./2);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z+1./2,-x,-y+1./2);wa(a,str);
          a.fpos=o+wv(-z,-x+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(-z+1./2,x+1./2,-y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y+1./2,z+1./2,-x);wa(a,str);
          a.fpos=o+wv(y+1./2,-z,-x+1./2);wa(a,str);
          a.fpos=o+wv(-y,-z+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(y+3./4,x+1./4,-z+3./4);wa(a,str);
          a.fpos=o+wv(-y+1./4,-x+1./4,-z+1./4);wa(a,str);
          a.fpos=o+wv(y+1./4,-x+3./4,z+3./4);wa(a,str);
          a.fpos=o+wv(-y+3./4,x+3./4,z+1./4);wa(a,str);
          a.fpos=o+wv(x+3./4,z+1./4,-y+3./4);wa(a,str);
          a.fpos=o+wv(-x+3./4,z+3./4,y+1./4);wa(a,str);
          a.fpos=o+wv(-x+1./4,-z+1./4,-y+1./4);wa(a,str);
          a.fpos=o+wv(x+1./4,-z+3./4,y+3./4);wa(a,str);
          a.fpos=o+wv(z+3./4,y+1./4,-x+3./4);wa(a,str);
          a.fpos=o+wv(z+1./4,-y+3./4,x+3./4);wa(a,str);
          a.fpos=o+wv(-z+3./4,y+3./4,x+1./4);wa(a,str);
          a.fpos=o+wv(-z+1./4,-y+1./4,-x+1./4);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 211  I432 #211
      case 211:
        str.spacegroup="I432";
        str.spacegrouplabel="#211";
        str.spacegroupnumber=211;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,-x,y);wa(a,str);
          a.fpos=o+wv(-z,x,-y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y,z,-x);wa(a,str);
          a.fpos=o+wv(y,-z,-x);wa(a,str);
          a.fpos=o+wv(-y,-z,x);wa(a,str);
          a.fpos=o+wv(y,x,-z);wa(a,str);
          a.fpos=o+wv(-y,-x,-z);wa(a,str);
          a.fpos=o+wv(y,-x,z);wa(a,str);
          a.fpos=o+wv(-y,x,z);wa(a,str);
          a.fpos=o+wv(x,z,-y);wa(a,str);
          a.fpos=o+wv(-x,z,y);wa(a,str);
          a.fpos=o+wv(-x,-z,-y);wa(a,str);
          a.fpos=o+wv(x,-z,y);wa(a,str);
          a.fpos=o+wv(z,y,-x);wa(a,str);
          a.fpos=o+wv(z,-y,x);wa(a,str);
          a.fpos=o+wv(-z,y,x);wa(a,str);
          a.fpos=o+wv(-z,-y,-x);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 212  P4_{3}32 #212
      case 212:
        str.spacegroup="P4_{3}32";
        str.spacegrouplabel="#212";
        str.spacegroupnumber=212;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
        a.fpos=o+wv(z,x,y);wa(a,str);
        a.fpos=o+wv(z+1./2,-x+1./2,-y);wa(a,str);
        a.fpos=o+wv(-z+1./2,-x,y+1./2);wa(a,str);
        a.fpos=o+wv(-z,x+1./2,-y+1./2);wa(a,str);
        a.fpos=o+wv(y,z,x);wa(a,str);
        a.fpos=o+wv(-y,z+1./2,-x+1./2);wa(a,str);
        a.fpos=o+wv(y+1./2,-z+1./2,-x);wa(a,str);
        a.fpos=o+wv(-y+1./2,-z,x+1./2);wa(a,str);
        a.fpos=o+wv(y+1./4,x+3./4,-z+3./4);wa(a,str);
        a.fpos=o+wv(-y+1./4,-x+1./4,-z+1./4);wa(a,str);
        a.fpos=o+wv(y+3./4,-x+3./4,z+1./4);wa(a,str);
        a.fpos=o+wv(-y+3./4,x+1./4,z+3./4);wa(a,str);
        a.fpos=o+wv(x+1./4,z+3./4,-y+3./4);wa(a,str);
        a.fpos=o+wv(-x+3./4,z+1./4,y+3./4);wa(a,str);
        a.fpos=o+wv(-x+1./4,-z+1./4,-y+1./4);wa(a,str);
        a.fpos=o+wv(x+3./4,-z+3./4,y+1./4);wa(a,str);
        a.fpos=o+wv(z+1./4,y+3./4,-x+3./4);wa(a,str);
        a.fpos=o+wv(z+3./4,-y+3./4,x+1./4);wa(a,str);
        a.fpos=o+wv(-z+3./4,y+1./4,x+3./4);wa(a,str);
        a.fpos=o+wv(-z+1./4,-y+1./4,-x+1./4);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 213  P4_{1}32 #213
      case 213:
        str.spacegroup="P4_{1}32";
        str.spacegrouplabel="#213";
        str.spacegroupnumber=213;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
        a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
        a.fpos=o+wv(z,x,y);wa(a,str);
        a.fpos=o+wv(z+1./2,-x+1./2,-y);wa(a,str);
        a.fpos=o+wv(-z+1./2,-x,y+1./2);wa(a,str);
        a.fpos=o+wv(-z,x+1./2,-y+1./2);wa(a,str);
        a.fpos=o+wv(y,z,x);wa(a,str);
        a.fpos=o+wv(-y,z+1./2,-x+1./2);wa(a,str);
        a.fpos=o+wv(y+1./2,-z+1./2,-x);wa(a,str);
        a.fpos=o+wv(-y+1./2,-z,x+1./2);wa(a,str);
        a.fpos=o+wv(y+3./4,x+1./4,-z+1./4);wa(a,str);
        a.fpos=o+wv(-y+3./4,-x+3./4,-z+3./4);wa(a,str);
        a.fpos=o+wv(y+1./4,-x+1./4,z+3./4);wa(a,str);
        a.fpos=o+wv(-y+1./4,x+3./4,z+1./4);wa(a,str);
        a.fpos=o+wv(x+3./4,z+1./4,-y+1./4);wa(a,str);
        a.fpos=o+wv(-x+1./4,z+3./4,y+1./4);wa(a,str);
        a.fpos=o+wv(-x+3./4,-z+3./4,-y+3./4);wa(a,str);
        a.fpos=o+wv(x+1./4,-z+1./4,y+3./4);wa(a,str);
        a.fpos=o+wv(z+3./4,y+1./4,-x+1./4);wa(a,str);
        a.fpos=o+wv(z+1./4,-y+1./4,x+3./4);wa(a,str);
        a.fpos=o+wv(-z+1./4,y+3./4,x+1./4);wa(a,str);
        a.fpos=o+wv(-z+3./4,-y+3./4,-x+3./4);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 214  I4_{1}32 #214
      case 214:
        str.spacegroup="I4_{1}32";
        str.spacegrouplabel="#214";
        str.spacegroupnumber=214;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z+1./2,-x+1./2,-y);wa(a,str);
          a.fpos=o+wv(-z+1./2,-x,y+1./2);wa(a,str);
          a.fpos=o+wv(-z,x+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y,z+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,-z+1./2,-x);wa(a,str);
          a.fpos=o+wv(-y+1./2,-z,x+1./2);wa(a,str);
          a.fpos=o+wv(y+3./4,x+1./4,-z+1./4);wa(a,str);
          a.fpos=o+wv(-y+3./4,-x+3./4,-z+3./4);wa(a,str);
          a.fpos=o+wv(y+1./4,-x+1./4,z+3./4);wa(a,str);
          a.fpos=o+wv(-y+1./4,x+3./4,z+1./4);wa(a,str);
          a.fpos=o+wv(x+3./4,z+1./4,-y+1./4);wa(a,str);
          a.fpos=o+wv(-x+1./4,z+3./4,y+1./4);wa(a,str);
          a.fpos=o+wv(-x+3./4,-z+3./4,-y+3./4);wa(a,str);
          a.fpos=o+wv(x+1./4,-z+1./4,y+3./4);wa(a,str);
          a.fpos=o+wv(z+3./4,y+1./4,-x+1./4);wa(a,str);
          a.fpos=o+wv(z+1./4,-y+1./4,x+3./4);wa(a,str);
          a.fpos=o+wv(-z+1./4,y+3./4,x+1./4);wa(a,str);
          a.fpos=o+wv(-z+3./4,-y+3./4,-x+3./4);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 215  P-43m #215
      case 215:
        str.spacegroup="P-43m";
        str.spacegrouplabel="#215";
        str.spacegroupnumber=215;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,-z);wa(a,str);
        a.fpos=o+wv(z,x,y);wa(a,str);
        a.fpos=o+wv(z,-x,-y);wa(a,str);
        a.fpos=o+wv(-z,-x,y);wa(a,str);
        a.fpos=o+wv(-z,x,-y);wa(a,str);
        a.fpos=o+wv(y,z,x);wa(a,str);
        a.fpos=o+wv(-y,z,-x);wa(a,str);
        a.fpos=o+wv(y,-z,-x);wa(a,str);
        a.fpos=o+wv(-y,-z,x);wa(a,str);
        a.fpos=o+wv(y,x,z);wa(a,str);
        a.fpos=o+wv(-y,-x,z);wa(a,str);
        a.fpos=o+wv(y,-x,-z);wa(a,str);
        a.fpos=o+wv(-y,x,-z);wa(a,str);
        a.fpos=o+wv(x,z,y);wa(a,str);
        a.fpos=o+wv(-x,z,-y);wa(a,str);
        a.fpos=o+wv(-x,-z,y);wa(a,str);
        a.fpos=o+wv(x,-z,-y);wa(a,str);
        a.fpos=o+wv(z,y,x);wa(a,str);
        a.fpos=o+wv(z,-y,-x);wa(a,str);
        a.fpos=o+wv(-z,y,-x);wa(a,str);
        a.fpos=o+wv(-z,-y,x);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 216  F-43m #216
      case 216:
        str.spacegroup="F-43m";
        str.spacegrouplabel="#216";
        str.spacegroupnumber=216;
        str.spacegroupoption="";
        for(uint j=1;j<=4;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(0,1./2,1./2);
          if(j==3) o=wv(1./2,0,1./2);
          if(j==4) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,-x,y);wa(a,str);
          a.fpos=o+wv(-z,x,-y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y,z,-x);wa(a,str);
          a.fpos=o+wv(y,-z,-x);wa(a,str);
          a.fpos=o+wv(-y,-z,x);wa(a,str);
          a.fpos=o+wv(y,x,z);wa(a,str);
          a.fpos=o+wv(-y,-x,z);wa(a,str);
          a.fpos=o+wv(y,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x,-z);wa(a,str);
          a.fpos=o+wv(x,z,y);wa(a,str);
          a.fpos=o+wv(-x,z,-y);wa(a,str);
          a.fpos=o+wv(-x,-z,y);wa(a,str);
          a.fpos=o+wv(x,-z,-y);wa(a,str);
          a.fpos=o+wv(z,y,x);wa(a,str);
          a.fpos=o+wv(z,-y,-x);wa(a,str);
          a.fpos=o+wv(-z,y,-x);wa(a,str);
          a.fpos=o+wv(-z,-y,x);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 217  I-43m #217
      case 217:
        str.spacegroup="I-43m";
        str.spacegrouplabel="#217";
        str.spacegroupnumber=217;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,-x,y);wa(a,str);
          a.fpos=o+wv(-z,x,-y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y,z,-x);wa(a,str);
          a.fpos=o+wv(y,-z,-x);wa(a,str);
          a.fpos=o+wv(-y,-z,x);wa(a,str);
          a.fpos=o+wv(y,x,z);wa(a,str);
          a.fpos=o+wv(-y,-x,z);wa(a,str);
          a.fpos=o+wv(y,-x,-z);wa(a,str);
          a.fpos=o+wv(-y,x,-z);wa(a,str);
          a.fpos=o+wv(x,z,y);wa(a,str);
          a.fpos=o+wv(-x,z,-y);wa(a,str);
          a.fpos=o+wv(-x,-z,y);wa(a,str);
          a.fpos=o+wv(x,-z,-y);wa(a,str);
          a.fpos=o+wv(z,y,x);wa(a,str);
          a.fpos=o+wv(z,-y,-x);wa(a,str);
          a.fpos=o+wv(-z,y,-x);wa(a,str);
          a.fpos=o+wv(-z,-y,x);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 218  P-43n #218
      case 218:
        str.spacegroup="P-43n";
        str.spacegrouplabel="#218";
        str.spacegroupnumber=218;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,-z);wa(a,str);
        a.fpos=o+wv(z,x,y);wa(a,str);
        a.fpos=o+wv(z,-x,-y);wa(a,str);
        a.fpos=o+wv(-z,-x,y);wa(a,str);
        a.fpos=o+wv(-z,x,-y);wa(a,str);
        a.fpos=o+wv(y,z,x);wa(a,str);
        a.fpos=o+wv(-y,z,-x);wa(a,str);
        a.fpos=o+wv(y,-z,-x);wa(a,str);
        a.fpos=o+wv(-y,-z,x);wa(a,str);
        a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(y+1./2,-x+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y+1./2,x+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,z+1./2,y+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,z+1./2,-y+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,-z+1./2,y+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-z+1./2,-y+1./2);wa(a,str);
        a.fpos=o+wv(z+1./2,y+1./2,x+1./2);wa(a,str);
        a.fpos=o+wv(z+1./2,-y+1./2,-x+1./2);wa(a,str);
        a.fpos=o+wv(-z+1./2,y+1./2,-x+1./2);wa(a,str);
        a.fpos=o+wv(-z+1./2,-y+1./2,x+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 219  F-43c #219
      case 219:
        str.spacegroup="F-43c";
        str.spacegrouplabel="#219";
        str.spacegroupnumber=219;
        str.spacegroupoption="";
        for(uint j=1;j<=4;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(0,1./2,1./2);
          if(j==3) o=wv(1./2,0,1./2);
          if(j==4) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,-x,y);wa(a,str);
          a.fpos=o+wv(-z,x,-y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y,z,-x);wa(a,str);
          a.fpos=o+wv(y,-z,-x);wa(a,str);
          a.fpos=o+wv(-y,-z,x);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,-x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,z+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,z+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,-z+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-z+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(z+1./2,y+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(z+1./2,-y+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(-z+1./2,y+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(-z+1./2,-y+1./2,x+1./2);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 220  I-43d #220
      case 220:
        str.spacegroup="I-43d";
        str.spacegrouplabel="#220";
        str.spacegroupnumber=220;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z+1./2,-x+1./2,-y);wa(a,str);
          a.fpos=o+wv(-z+1./2,-x,y+1./2);wa(a,str);
          a.fpos=o+wv(-z,x+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y,z+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,-z+1./2,-x);wa(a,str);
          a.fpos=o+wv(-y+1./2,-z,x+1./2);wa(a,str);
          a.fpos=o+wv(y+1./4,x+1./4,z+1./4);wa(a,str);
          a.fpos=o+wv(-y+1./4,-x+3./4,z+3./4);wa(a,str);
          a.fpos=o+wv(y+3./4,-x+1./4,-z+3./4);wa(a,str);
          a.fpos=o+wv(-y+3./4,x+3./4,-z+1./4);wa(a,str);
          a.fpos=o+wv(x+1./4,z+1./4,y+1./4);wa(a,str);
          a.fpos=o+wv(-x+3./4,z+3./4,-y+1./4);wa(a,str);
          a.fpos=o+wv(-x+1./4,-z+3./4,y+3./4);wa(a,str);
          a.fpos=o+wv(x+3./4,-z+1./4,-y+3./4);wa(a,str);
          a.fpos=o+wv(z+1./4,y+1./4,x+1./4);wa(a,str);
          a.fpos=o+wv(z+3./4,-y+1./4,-x+3./4);wa(a,str);
          a.fpos=o+wv(-z+3./4,y+3./4,-x+1./4);wa(a,str);
          a.fpos=o+wv(-z+1./4,-y+3./4,x+3./4);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 221  Pm-3m #221
      case 221:
        str.spacegroup="Pm-3m";
        str.spacegrouplabel="#221";
        str.spacegroupnumber=221;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,-z);wa(a,str);
        a.fpos=o+wv(z,x,y);wa(a,str);
        a.fpos=o+wv(z,-x,-y);wa(a,str);
        a.fpos=o+wv(-z,-x,y);wa(a,str);
        a.fpos=o+wv(-z,x,-y);wa(a,str);
        a.fpos=o+wv(y,z,x);wa(a,str);
        a.fpos=o+wv(-y,z,-x);wa(a,str);
        a.fpos=o+wv(y,-z,-x);wa(a,str);
        a.fpos=o+wv(-y,-z,x);wa(a,str);
        a.fpos=o+wv(y,x,-z);wa(a,str);
        a.fpos=o+wv(-y,-x,-z);wa(a,str);
        a.fpos=o+wv(y,-x,z);wa(a,str);
        a.fpos=o+wv(-y,x,z);wa(a,str);
        a.fpos=o+wv(x,z,-y);wa(a,str);
        a.fpos=o+wv(-x,z,y);wa(a,str);
        a.fpos=o+wv(-x,-z,-y);wa(a,str);
        a.fpos=o+wv(x,-z,y);wa(a,str);
        a.fpos=o+wv(z,y,-x);wa(a,str);
        a.fpos=o+wv(z,-y,x);wa(a,str);
        a.fpos=o+wv(-z,y,x);wa(a,str);
        a.fpos=o+wv(-z,-y,-x);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,z);wa(a,str);
        a.fpos=o+wv(-z,-x,-y);wa(a,str);
        a.fpos=o+wv(-z,x,y);wa(a,str);
        a.fpos=o+wv(z,x,-y);wa(a,str);
        a.fpos=o+wv(z,-x,y);wa(a,str);
        a.fpos=o+wv(-y,-z,-x);wa(a,str);
        a.fpos=o+wv(y,-z,x);wa(a,str);
        a.fpos=o+wv(-y,z,x);wa(a,str);
        a.fpos=o+wv(y,z,-x);wa(a,str);
        a.fpos=o+wv(-y,-x,z);wa(a,str);
        a.fpos=o+wv(y,x,z);wa(a,str);
        a.fpos=o+wv(-y,x,-z);wa(a,str);
        a.fpos=o+wv(y,-x,-z);wa(a,str);
        a.fpos=o+wv(-x,-z,y);wa(a,str);
        a.fpos=o+wv(x,-z,-y);wa(a,str);
        a.fpos=o+wv(x,z,y);wa(a,str);
        a.fpos=o+wv(-x,z,-y);wa(a,str);
        a.fpos=o+wv(-z,-y,x);wa(a,str);
        a.fpos=o+wv(-z,y,-x);wa(a,str);
        a.fpos=o+wv(z,-y,-x);wa(a,str);
        a.fpos=o+wv(z,y,x);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 222  Pn-3n #222
      case 222:
        if(option==1) {
          str.spacegroup="Pn-3n";
          str.spacegrouplabel="#222";
          str.spacegroupnumber=222;
          str.spacegroupoption="origin choice 1";
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,-x,y);wa(a,str);
          a.fpos=o+wv(-z,x,-y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y,z,-x);wa(a,str);
          a.fpos=o+wv(y,-z,-x);wa(a,str);
          a.fpos=o+wv(-y,-z,x);wa(a,str);
          a.fpos=o+wv(y,x,-z);wa(a,str);
          a.fpos=o+wv(-y,-x,-z);wa(a,str);
          a.fpos=o+wv(y,-x,z);wa(a,str);
          a.fpos=o+wv(-y,x,z);wa(a,str);
          a.fpos=o+wv(x,z,-y);wa(a,str);
          a.fpos=o+wv(-x,z,y);wa(a,str);
          a.fpos=o+wv(-x,-z,-y);wa(a,str);
          a.fpos=o+wv(x,-z,y);wa(a,str);
          a.fpos=o+wv(z,y,-x);wa(a,str);
          a.fpos=o+wv(z,-y,x);wa(a,str);
          a.fpos=o+wv(-z,y,x);wa(a,str);
          a.fpos=o+wv(-z,-y,-x);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-z+1./2,-x+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(-z+1./2,x+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(z+1./2,x+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(z+1./2,-x+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,-z+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,-z+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,z+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,z+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,-x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,-z+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-z+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,z+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,z+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(-z+1./2,-y+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(-z+1./2,y+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(z+1./2,-y+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(z+1./2,y+1./2,x+1./2);wa(a,str);
        } else {
          str.spacegroup="Pn-3n";
          str.spacegrouplabel="#222";
          str.spacegroupnumber=222;
          str.spacegroupoption="origin choice 2";
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y,-z+1./2);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z,-x+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(-z+1./2,-x+1./2,y);wa(a,str);
          a.fpos=o+wv(-z+1./2,x,-y+1./2);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y+1./2,z,-x+1./2);wa(a,str);
          a.fpos=o+wv(y,-z+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,-z+1./2,x);wa(a,str);
          a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(y,-x+1./2,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,x,z);wa(a,str);
          a.fpos=o+wv(x,z,-y+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,z,y);wa(a,str);
          a.fpos=o+wv(-x+1./2,-z+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(x,-z+1./2,y);wa(a,str);
          a.fpos=o+wv(z,y,-x+1./2);wa(a,str);
          a.fpos=o+wv(z,-y+1./2,x);wa(a,str);
          a.fpos=o+wv(-z+1./2,y,x);wa(a,str);
          a.fpos=o+wv(-z+1./2,-y+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,x+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(z+1./2,x+1./2,-y);wa(a,str);
          a.fpos=o+wv(z+1./2,-x,y+1./2);wa(a,str);
          a.fpos=o+wv(-y,-z,-x);wa(a,str);
          a.fpos=o+wv(y+1./2,-z,x+1./2);wa(a,str);
          a.fpos=o+wv(-y,z+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,z+1./2,-x);wa(a,str);
          a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-y,x+1./2,-z);wa(a,str);
          a.fpos=o+wv(y+1./2,-x,-z);wa(a,str);
          a.fpos=o+wv(-x,-z,y+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-z,-y);wa(a,str);
          a.fpos=o+wv(x+1./2,z+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(-x,z+1./2,-y);wa(a,str);
          a.fpos=o+wv(-z,-y,x+1./2);wa(a,str);
          a.fpos=o+wv(-z,y+1./2,-x);wa(a,str);
          a.fpos=o+wv(z+1./2,-y,-x);wa(a,str);
          a.fpos=o+wv(z+1./2,y+1./2,x+1./2);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 223  Pm-3n #223
      case 223:
        str.spacegroup="Pm-3n";
        str.spacegrouplabel="#223";
        str.spacegroupnumber=223;
        str.spacegroupoption="";
        o=wv(0,0,0);
        a.fpos=o+wv(x,y,z);wa(a,str);
        a.fpos=o+wv(-x,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,-z);wa(a,str);
        a.fpos=o+wv(z,x,y);wa(a,str);
        a.fpos=o+wv(z,-x,-y);wa(a,str);
        a.fpos=o+wv(-z,-x,y);wa(a,str);
        a.fpos=o+wv(-z,x,-y);wa(a,str);
        a.fpos=o+wv(y,z,x);wa(a,str);
        a.fpos=o+wv(-y,z,-x);wa(a,str);
        a.fpos=o+wv(y,-z,-x);wa(a,str);
        a.fpos=o+wv(-y,-z,x);wa(a,str);
        a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,z+1./2,-y+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,z+1./2,y+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,-z+1./2,-y+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-z+1./2,y+1./2);wa(a,str);
        a.fpos=o+wv(z+1./2,y+1./2,-x+1./2);wa(a,str);
        a.fpos=o+wv(z+1./2,-y+1./2,x+1./2);wa(a,str);
        a.fpos=o+wv(-z+1./2,y+1./2,x+1./2);wa(a,str);
        a.fpos=o+wv(-z+1./2,-y+1./2,-x+1./2);wa(a,str);
        a.fpos=o+wv(-x,-y,-z);wa(a,str);
        a.fpos=o+wv(x,y,-z);wa(a,str);
        a.fpos=o+wv(x,-y,z);wa(a,str);
        a.fpos=o+wv(-x,y,z);wa(a,str);
        a.fpos=o+wv(-z,-x,-y);wa(a,str);
        a.fpos=o+wv(-z,x,y);wa(a,str);
        a.fpos=o+wv(z,x,-y);wa(a,str);
        a.fpos=o+wv(z,-x,y);wa(a,str);
        a.fpos=o+wv(-y,-z,-x);wa(a,str);
        a.fpos=o+wv(y,-z,x);wa(a,str);
        a.fpos=o+wv(-y,z,x);wa(a,str);
        a.fpos=o+wv(y,z,-x);wa(a,str);
        a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
        a.fpos=o+wv(-y+1./2,x+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(y+1./2,-x+1./2,-z+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,-z+1./2,y+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,-z+1./2,-y+1./2);wa(a,str);
        a.fpos=o+wv(x+1./2,z+1./2,y+1./2);wa(a,str);
        a.fpos=o+wv(-x+1./2,z+1./2,-y+1./2);wa(a,str);
        a.fpos=o+wv(-z+1./2,-y+1./2,x+1./2);wa(a,str);
        a.fpos=o+wv(-z+1./2,y+1./2,-x+1./2);wa(a,str);
        a.fpos=o+wv(z+1./2,-y+1./2,-x+1./2);wa(a,str);
        a.fpos=o+wv(z+1./2,y+1./2,x+1./2);wa(a,str);
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 224  Pn-3m #224
      case 224:
        if(option==1) {
          str.spacegroup="Pn-3m";
          str.spacegrouplabel="#224";
          str.spacegroupnumber=224;
          str.spacegroupoption="origin choice 1";
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,-x,y);wa(a,str);
          a.fpos=o+wv(-z,x,-y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y,z,-x);wa(a,str);
          a.fpos=o+wv(y,-z,-x);wa(a,str);
          a.fpos=o+wv(-y,-z,x);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,z+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,z+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,-z+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-z+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(z+1./2,y+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(z+1./2,-y+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(-z+1./2,y+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(-z+1./2,-y+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-z+1./2,-x+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(-z+1./2,x+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(z+1./2,x+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(z+1./2,-x+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,-z+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,-z+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,z+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,z+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(-y,-x,z);wa(a,str);
          a.fpos=o+wv(y,x,z);wa(a,str);
          a.fpos=o+wv(-y,x,-z);wa(a,str);
          a.fpos=o+wv(y,-x,-z);wa(a,str);
          a.fpos=o+wv(-x,-z,y);wa(a,str);
          a.fpos=o+wv(x,-z,-y);wa(a,str);
          a.fpos=o+wv(x,z,y);wa(a,str);
          a.fpos=o+wv(-x,z,-y);wa(a,str);
          a.fpos=o+wv(-z,-y,x);wa(a,str);
          a.fpos=o+wv(-z,y,-x);wa(a,str);
          a.fpos=o+wv(z,-y,-x);wa(a,str);
          a.fpos=o+wv(z,y,x);wa(a,str);
        } else {
          str.spacegroup="Pn-3m";
          str.spacegrouplabel="#224";
          str.spacegroupnumber=224;
          str.spacegroupoption="origin choice 2";
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,y,-z+1./2);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z,-x+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(-z+1./2,-x+1./2,y);wa(a,str);
          a.fpos=o+wv(-z+1./2,x,-y+1./2);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y+1./2,z,-x+1./2);wa(a,str);
          a.fpos=o+wv(y,-z+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,-z+1./2,x);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,-z);wa(a,str);
          a.fpos=o+wv(-y,-x,-z);wa(a,str);
          a.fpos=o+wv(y+1./2,-x,z+1./2);wa(a,str);
          a.fpos=o+wv(-y,x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,z+1./2,-y);wa(a,str);
          a.fpos=o+wv(-x,z+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(-x,-z,-y);wa(a,str);
          a.fpos=o+wv(x+1./2,-z,y+1./2);wa(a,str);
          a.fpos=o+wv(z+1./2,y+1./2,-x);wa(a,str);
          a.fpos=o+wv(z+1./2,-y,x+1./2);wa(a,str);
          a.fpos=o+wv(-z,y+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(-z,-y,-x);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,x+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(z+1./2,x+1./2,-y);wa(a,str);
          a.fpos=o+wv(z+1./2,-x,y+1./2);wa(a,str);
          a.fpos=o+wv(-y,-z,-x);wa(a,str);
          a.fpos=o+wv(y+1./2,-z,x+1./2);wa(a,str);
          a.fpos=o+wv(-y,z+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,z+1./2,-x);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,z);wa(a,str);
          a.fpos=o+wv(y,x,z);wa(a,str);
          a.fpos=o+wv(-y+1./2,x,-z+1./2);wa(a,str);
          a.fpos=o+wv(y,-x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,-z+1./2,y);wa(a,str);
          a.fpos=o+wv(x,-z+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(x,z,y);wa(a,str);
          a.fpos=o+wv(-x+1./2,z,-y+1./2);wa(a,str);
          a.fpos=o+wv(-z+1./2,-y+1./2,x);wa(a,str);
          a.fpos=o+wv(-z+1./2,y,-x+1./2);wa(a,str);
          a.fpos=o+wv(z,-y+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(z,y,x);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 225  Fm-3m #225
      case 225:
        str.spacegroup="Fm-3m";
        str.spacegrouplabel="#225";
        str.spacegroupnumber=225;
        str.spacegroupoption="";
        for(uint j=1;j<=4;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(0,1./2,1./2);
          if(j==3) o=wv(1./2,0,1./2);
          if(j==4) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,-x,y);wa(a,str);
          a.fpos=o+wv(-z,x,-y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y,z,-x);wa(a,str);
          a.fpos=o+wv(y,-z,-x);wa(a,str);
          a.fpos=o+wv(-y,-z,x);wa(a,str);
          a.fpos=o+wv(y,x,-z);wa(a,str);
          a.fpos=o+wv(-y,-x,-z);wa(a,str);
          a.fpos=o+wv(y,-x,z);wa(a,str);
          a.fpos=o+wv(-y,x,z);wa(a,str);
          a.fpos=o+wv(x,z,-y);wa(a,str);
          a.fpos=o+wv(-x,z,y);wa(a,str);
          a.fpos=o+wv(-x,-z,-y);wa(a,str);
          a.fpos=o+wv(x,-z,y);wa(a,str);
          a.fpos=o+wv(z,y,-x);wa(a,str);
          a.fpos=o+wv(z,-y,x);wa(a,str);
          a.fpos=o+wv(-z,y,x);wa(a,str);
          a.fpos=o+wv(-z,-y,-x);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
          a.fpos=o+wv(-z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,x,y);wa(a,str);
          a.fpos=o+wv(z,x,-y);wa(a,str);
          a.fpos=o+wv(z,-x,y);wa(a,str);
          a.fpos=o+wv(-y,-z,-x);wa(a,str);
          a.fpos=o+wv(y,-z,x);wa(a,str);
          a.fpos=o+wv(-y,z,x);wa(a,str);
          a.fpos=o+wv(y,z,-x);wa(a,str);
          a.fpos=o+wv(-y,-x,z);wa(a,str);
          a.fpos=o+wv(y,x,z);wa(a,str);
          a.fpos=o+wv(-y,x,-z);wa(a,str);
          a.fpos=o+wv(y,-x,-z);wa(a,str);
          a.fpos=o+wv(-x,-z,y);wa(a,str);
          a.fpos=o+wv(x,-z,-y);wa(a,str);
          a.fpos=o+wv(x,z,y);wa(a,str);
          a.fpos=o+wv(-x,z,-y);wa(a,str);
          a.fpos=o+wv(-z,-y,x);wa(a,str);
          a.fpos=o+wv(-z,y,-x);wa(a,str);
          a.fpos=o+wv(z,-y,-x);wa(a,str);
          a.fpos=o+wv(z,y,x);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 226  Fm-3c #226
      case 226:
        str.spacegroup="Fm-3c";
        str.spacegrouplabel="#226";
        str.spacegroupnumber=226;
        str.spacegroupoption="";
        for(uint j=1;j<=4;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(0,1./2,1./2);
          if(j==3) o=wv(1./2,0,1./2);
          if(j==4) o=wv(1./2,1./2,0);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,-x,y);wa(a,str);
          a.fpos=o+wv(-z,x,-y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y,z,-x);wa(a,str);
          a.fpos=o+wv(y,-z,-x);wa(a,str);
          a.fpos=o+wv(-y,-z,x);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,z+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,z+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,-z+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-z+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(z+1./2,y+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(z+1./2,-y+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(-z+1./2,y+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(-z+1./2,-y+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
          a.fpos=o+wv(-z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,x,y);wa(a,str);
          a.fpos=o+wv(z,x,-y);wa(a,str);
          a.fpos=o+wv(z,-x,y);wa(a,str);
          a.fpos=o+wv(-y,-z,-x);wa(a,str);
          a.fpos=o+wv(y,-z,x);wa(a,str);
          a.fpos=o+wv(-y,z,x);wa(a,str);
          a.fpos=o+wv(y,z,-x);wa(a,str);
          a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,-x+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,-z+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-z+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,z+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,z+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(-z+1./2,-y+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(-z+1./2,y+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(z+1./2,-y+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(z+1./2,y+1./2,x+1./2);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 227  Fd-3m #227
      case 227:
        if(option==1) {
          str.spacegroup="Fd-3m";
          str.spacegrouplabel="#227";
          str.spacegroupnumber=227;
          str.spacegroupoption="origin choice 1";
          for(uint j=1;j<=4;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(0,1./2,1./2);
            if(j==3) o=wv(1./2,0,1./2);
            if(j==4) o=wv(1./2,1./2,0);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x,-y+1./2,z+1./2);wa(a,str);
            a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
            a.fpos=o+wv(x+1./2,-y,-z+1./2);wa(a,str);
            a.fpos=o+wv(z,x,y);wa(a,str);
            a.fpos=o+wv(z+1./2,-x,-y+1./2);wa(a,str);
            a.fpos=o+wv(-z,-x+1./2,y+1./2);wa(a,str);
            a.fpos=o+wv(-z+1./2,x+1./2,-y);wa(a,str);
            a.fpos=o+wv(y,z,x);wa(a,str);
            a.fpos=o+wv(-y+1./2,z+1./2,-x);wa(a,str);
            a.fpos=o+wv(y+1./2,-z,-x+1./2);wa(a,str);
            a.fpos=o+wv(-y,-z+1./2,x+1./2);wa(a,str);
            a.fpos=o+wv(y+3./4,x+1./4,-z+3./4);wa(a,str);
            a.fpos=o+wv(-y+1./4,-x+1./4,-z+1./4);wa(a,str);
            a.fpos=o+wv(y+1./4,-x+3./4,z+3./4);wa(a,str);
            a.fpos=o+wv(-y+3./4,x+3./4,z+1./4);wa(a,str);
            a.fpos=o+wv(x+3./4,z+1./4,-y+3./4);wa(a,str);
            a.fpos=o+wv(-x+3./4,z+3./4,y+1./4);wa(a,str);
            a.fpos=o+wv(-x+1./4,-z+1./4,-y+1./4);wa(a,str);
            a.fpos=o+wv(x+1./4,-z+3./4,y+3./4);wa(a,str);
            a.fpos=o+wv(z+3./4,y+1./4,-x+3./4);wa(a,str);
            a.fpos=o+wv(z+1./4,-y+3./4,x+3./4);wa(a,str);
            a.fpos=o+wv(-z+3./4,y+3./4,x+1./4);wa(a,str);
            a.fpos=o+wv(-z+1./4,-y+1./4,-x+1./4);wa(a,str);
            a.fpos=o+wv(-x+1./4,-y+1./4,-z+1./4);wa(a,str);
            a.fpos=o+wv(x+1./4,y+3./4,-z+3./4);wa(a,str);
            a.fpos=o+wv(x+3./4,-y+3./4,z+1./4);wa(a,str);
            a.fpos=o+wv(-x+3./4,y+1./4,z+3./4);wa(a,str);
            a.fpos=o+wv(-z+1./4,-x+1./4,-y+1./4);wa(a,str);
            a.fpos=o+wv(-z+3./4,x+1./4,y+3./4);wa(a,str);
            a.fpos=o+wv(z+1./4,x+3./4,-y+3./4);wa(a,str);
            a.fpos=o+wv(z+3./4,-x+3./4,y+1./4);wa(a,str);
            a.fpos=o+wv(-y+1./4,-z+1./4,-x+1./4);wa(a,str);
            a.fpos=o+wv(y+3./4,-z+3./4,x+1./4);wa(a,str);
            a.fpos=o+wv(-y+3./4,z+1./4,x+3./4);wa(a,str);
            a.fpos=o+wv(y+1./4,z+3./4,-x+3./4);wa(a,str);
            a.fpos=o+wv(-y+1./2,-x,z+1./2);wa(a,str);
            a.fpos=o+wv(y,x,z);wa(a,str);
            a.fpos=o+wv(-y,x+1./2,-z+1./2);wa(a,str);
            a.fpos=o+wv(y+1./2,-x+1./2,-z);wa(a,str);
            a.fpos=o+wv(-x+1./2,-z,y+1./2);wa(a,str);
            a.fpos=o+wv(x+1./2,-z+1./2,-y);wa(a,str);
            a.fpos=o+wv(x,z,y);wa(a,str);
            a.fpos=o+wv(-x,z+1./2,-y+1./2);wa(a,str);
            a.fpos=o+wv(-z+1./2,-y,x+1./2);wa(a,str);
            a.fpos=o+wv(-z,y+1./2,-x+1./2);wa(a,str);
            a.fpos=o+wv(z+1./2,-y+1./2,-x);wa(a,str);
            a.fpos=o+wv(z,y,x);wa(a,str);
          }
        } else {
        //   cerr << "aflow_wyckoff.cpp: spacegroup=" << spacegroup << " option=" << option << endl;
          str.spacegroup="Fd-3m";
          str.spacegrouplabel="#227";
          str.spacegroupnumber=227;
          str.spacegroupoption="origin choice 2";
          for(uint j=1;j<=4;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(0,1./2,1./2);
            if(j==3) o=wv(1./2,0,1./2);
            if(j==4) o=wv(1./2,1./2,0);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x+3./4,-y+1./4,z+1./2);wa(a,str);
            a.fpos=o+wv(-x+1./4,y+1./2,-z+3./4);wa(a,str);
            a.fpos=o+wv(x+1./2,-y+3./4,-z+1./4);wa(a,str);
            a.fpos=o+wv(z,x,y);wa(a,str);
            a.fpos=o+wv(z+1./2,-x+3./4,-y+1./4);wa(a,str);
            a.fpos=o+wv(-z+3./4,-x+1./4,y+1./2);wa(a,str);
            a.fpos=o+wv(-z+1./4,x+1./2,-y+3./4);wa(a,str);
            a.fpos=o+wv(y,z,x);wa(a,str);
            a.fpos=o+wv(-y+1./4,z+1./2,-x+3./4);wa(a,str);
            a.fpos=o+wv(y+1./2,-z+3./4,-x+1./4);wa(a,str);
            a.fpos=o+wv(-y+3./4,-z+1./4,x+1./2);wa(a,str);
            a.fpos=o+wv(y+3./4,x+1./4,-z+1./2);wa(a,str);
            a.fpos=o+wv(-y,-x,-z);wa(a,str);
            a.fpos=o+wv(y+1./4,-x+1./2,z+3./4);wa(a,str);
            a.fpos=o+wv(-y+1./2,x+3./4,z+1./4);wa(a,str);
            a.fpos=o+wv(x+3./4,z+1./4,-y+1./2);wa(a,str);
            a.fpos=o+wv(-x+1./2,z+3./4,y+1./4);wa(a,str);
            a.fpos=o+wv(-x,-z,-y);wa(a,str);
            a.fpos=o+wv(x+1./4,-z+1./2,y+3./4);wa(a,str);
            a.fpos=o+wv(z+3./4,y+1./4,-x+1./2);wa(a,str);
            a.fpos=o+wv(z+1./4,-y+1./2,x+3./4);wa(a,str);
            a.fpos=o+wv(-z+1./2,y+3./4,x+1./4);wa(a,str);
            a.fpos=o+wv(-z,-y,-x);wa(a,str);
            a.fpos=o+wv(-x,-y,-z);wa(a,str);
            a.fpos=o+wv(x+1./4,y+3./4,-z+1./2);wa(a,str);
            a.fpos=o+wv(x+3./4,-y+1./2,z+1./4);wa(a,str);
            a.fpos=o+wv(-x+1./2,y+1./4,z+3./4);wa(a,str);
            a.fpos=o+wv(-z,-x,-y);wa(a,str);
            a.fpos=o+wv(-z+1./2,x+1./4,y+3./4);wa(a,str);
            a.fpos=o+wv(z+1./4,x+3./4,-y+1./2);wa(a,str);
            a.fpos=o+wv(z+3./4,-x+1./2,y+1./4);wa(a,str);
            a.fpos=o+wv(-y,-z,-x);wa(a,str);
            a.fpos=o+wv(y+3./4,-z+1./2,x+1./4);wa(a,str);
            a.fpos=o+wv(-y+1./2,z+1./4,x+3./4);wa(a,str);
            a.fpos=o+wv(y+1./4,z+3./4,-x+1./2);wa(a,str);
            a.fpos=o+wv(-y+1./4,-x+3./4,z+1./2);wa(a,str);
            a.fpos=o+wv(y,x,z);wa(a,str);
            a.fpos=o+wv(-y+3./4,x+1./2,-z+1./4);wa(a,str);
            a.fpos=o+wv(y+1./2,-x+1./4,-z+3./4);wa(a,str);
            a.fpos=o+wv(-x+1./4,-z+3./4,y+1./2);wa(a,str);
            a.fpos=o+wv(x+1./2,-z+1./4,-y+3./4);wa(a,str);
            a.fpos=o+wv(x,z,y);wa(a,str);
            a.fpos=o+wv(-x+3./4,z+1./2,-y+1./4);wa(a,str);
            a.fpos=o+wv(-z+1./4,-y+3./4,x+1./2);wa(a,str);
            a.fpos=o+wv(-z+3./4,y+1./2,-x+1./4);wa(a,str);
            a.fpos=o+wv(z+1./2,-y+1./4,-x+3./4);wa(a,str);
            a.fpos=o+wv(z,y,x);wa(a,str);
          }
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 228  Fd-3c #228
      case 228:
        if(option==1) {
          str.spacegroup="Fd-3c";
          str.spacegrouplabel="#228";
          str.spacegroupnumber=228;
          str.spacegroupoption="origin choice 1";
          for(uint j=1;j<=4;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(0,1./2,1./2);
            if(j==3) o=wv(1./2,0,1./2);
            if(j==4) o=wv(1./2,1./2,0);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x,-y+1./2,z+1./2);wa(a,str);
            a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
            a.fpos=o+wv(x+1./2,-y,-z+1./2);wa(a,str);
            a.fpos=o+wv(z,x,y);wa(a,str);
            a.fpos=o+wv(z+1./2,-x,-y+1./2);wa(a,str);
            a.fpos=o+wv(-z,-x+1./2,y+1./2);wa(a,str);
            a.fpos=o+wv(-z+1./2,x+1./2,-y);wa(a,str);
            a.fpos=o+wv(y,z,x);wa(a,str);
            a.fpos=o+wv(-y+1./2,z+1./2,-x);wa(a,str);
            a.fpos=o+wv(y+1./2,-z,-x+1./2);wa(a,str);
            a.fpos=o+wv(-y,-z+1./2,x+1./2);wa(a,str);
            a.fpos=o+wv(y+3./4,x+1./4,-z+3./4);wa(a,str);
            a.fpos=o+wv(-y+1./4,-x+1./4,-z+1./4);wa(a,str);
            a.fpos=o+wv(y+1./4,-x+3./4,z+3./4);wa(a,str);
            a.fpos=o+wv(-y+3./4,x+3./4,z+1./4);wa(a,str);
            a.fpos=o+wv(x+3./4,z+1./4,-y+3./4);wa(a,str);
            a.fpos=o+wv(-x+3./4,z+3./4,y+1./4);wa(a,str);
            a.fpos=o+wv(-x+1./4,-z+1./4,-y+1./4);wa(a,str);
            a.fpos=o+wv(x+1./4,-z+3./4,y+3./4);wa(a,str);
            a.fpos=o+wv(z+3./4,y+1./4,-x+3./4);wa(a,str);
            a.fpos=o+wv(z+1./4,-y+3./4,x+3./4);wa(a,str);
            a.fpos=o+wv(-z+3./4,y+3./4,x+1./4);wa(a,str);
            a.fpos=o+wv(-z+1./4,-y+1./4,-x+1./4);wa(a,str);
            a.fpos=o+wv(-x+3./4,-y+3./4,-z+3./4);wa(a,str);
            a.fpos=o+wv(x+3./4,y+1./4,-z+1./4);wa(a,str);
            a.fpos=o+wv(x+1./4,-y+1./4,z+3./4);wa(a,str);
            a.fpos=o+wv(-x+1./4,y+3./4,z+1./4);wa(a,str);
            a.fpos=o+wv(-z+3./4,-x+3./4,-y+3./4);wa(a,str);
            a.fpos=o+wv(-z+1./4,x+3./4,y+1./4);wa(a,str);
            a.fpos=o+wv(z+3./4,x+1./4,-y+1./4);wa(a,str);
            a.fpos=o+wv(z+1./4,-x+1./4,y+3./4);wa(a,str);
            a.fpos=o+wv(-y+3./4,-z+3./4,-x+3./4);wa(a,str);
            a.fpos=o+wv(y+1./4,-z+1./4,x+3./4);wa(a,str);
            a.fpos=o+wv(-y+1./4,z+3./4,x+1./4);wa(a,str);
            a.fpos=o+wv(y+3./4,z+1./4,-x+1./4);wa(a,str);
            a.fpos=o+wv(-y,-x+1./2,z);wa(a,str);
            a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
            a.fpos=o+wv(-y+1./2,x,-z);wa(a,str);
            a.fpos=o+wv(y,-x,-z+1./2);wa(a,str);
            a.fpos=o+wv(-x,-z+1./2,y);wa(a,str);
            a.fpos=o+wv(x,-z,-y+1./2);wa(a,str);
            a.fpos=o+wv(x+1./2,z+1./2,y+1./2);wa(a,str);
            a.fpos=o+wv(-x+1./2,z,-y);wa(a,str);
            a.fpos=o+wv(-z,-y+1./2,x);wa(a,str);
            a.fpos=o+wv(-z+1./2,y,-x);wa(a,str);
            a.fpos=o+wv(z,-y,-x+1./2);wa(a,str);
            a.fpos=o+wv(z+1./2,y+1./2,x+1./2);wa(a,str);
          }
        } else {
          str.spacegroup="Fd-3c";
          str.spacegrouplabel="#228";
          str.spacegroupnumber=228;
          str.spacegroupoption="origin choice 2";
          for(uint j=1;j<=4;j++) {
            if(j==1) o=wv(0,0,0);
            if(j==2) o=wv(0,1./2,1./2);
            if(j==3) o=wv(1./2,0,1./2);
            if(j==4) o=wv(1./2,1./2,0);
            a.fpos=o+wv(x,y,z);wa(a,str);
            a.fpos=o+wv(-x+1./4,-y+3./4,z+1./2);wa(a,str);
            a.fpos=o+wv(-x+3./4,y+1./2,-z+1./4);wa(a,str);
            a.fpos=o+wv(x+1./2,-y+1./4,-z+3./4);wa(a,str);
            a.fpos=o+wv(z,x,y);wa(a,str);
            a.fpos=o+wv(z+1./2,-x+1./4,-y+3./4);wa(a,str);
            a.fpos=o+wv(-z+1./4,-x+3./4,y+1./2);wa(a,str);
            a.fpos=o+wv(-z+3./4,x+1./2,-y+1./4);wa(a,str);
            a.fpos=o+wv(y,z,x);wa(a,str);
            a.fpos=o+wv(-y+3./4,z+1./2,-x+1./4);wa(a,str);
            a.fpos=o+wv(y+1./2,-z+1./4,-x+3./4);wa(a,str);
            a.fpos=o+wv(-y+1./4,-z+3./4,x+1./2);wa(a,str);
            a.fpos=o+wv(y+3./4,x+1./4,-z);wa(a,str);
            a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
            a.fpos=o+wv(y+1./4,-x,z+3./4);wa(a,str);
            a.fpos=o+wv(-y,x+3./4,z+1./4);wa(a,str);
            a.fpos=o+wv(x+3./4,z+1./4,-y);wa(a,str);
            a.fpos=o+wv(-x,z+3./4,y+1./4);wa(a,str);
            a.fpos=o+wv(-x+1./2,-z+1./2,-y+1./2);wa(a,str);
            a.fpos=o+wv(x+1./4,-z,y+3./4);wa(a,str);
            a.fpos=o+wv(z+3./4,y+1./4,-x);wa(a,str);
            a.fpos=o+wv(z+1./4,-y,x+3./4);wa(a,str);
            a.fpos=o+wv(-z,y+3./4,x+1./4);wa(a,str);
            a.fpos=o+wv(-z+1./2,-y+1./2,-x+1./2);wa(a,str);
            a.fpos=o+wv(-x,-y,-z);wa(a,str);
            a.fpos=o+wv(x+3./4,y+1./4,-z+1./2);wa(a,str);
            a.fpos=o+wv(x+1./4,-y+1./2,z+3./4);wa(a,str);
            a.fpos=o+wv(-x+1./2,y+3./4,z+1./4);wa(a,str);
            a.fpos=o+wv(-z,-x,-y);wa(a,str);
            a.fpos=o+wv(-z+1./2,x+3./4,y+1./4);wa(a,str);
            a.fpos=o+wv(z+3./4,x+1./4,-y+1./2);wa(a,str);
            a.fpos=o+wv(z+1./4,-x+1./2,y+3./4);wa(a,str);
            a.fpos=o+wv(-y,-z,-x);wa(a,str);
            a.fpos=o+wv(y+1./4,-z+1./2,x+3./4);wa(a,str);
            a.fpos=o+wv(-y+1./2,z+3./4,x+1./4);wa(a,str);
            a.fpos=o+wv(y+3./4,z+1./4,-x+1./2);wa(a,str);
            a.fpos=o+wv(-y+1./4,-x+3./4,z);wa(a,str);
            a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
            a.fpos=o+wv(-y+3./4,x,-z+1./4);wa(a,str);
            a.fpos=o+wv(y,-x+1./4,-z+3./4);wa(a,str);
            a.fpos=o+wv(-x+1./4,-z+3./4,y);wa(a,str);
            a.fpos=o+wv(x,-z+1./4,-y+3./4);wa(a,str);
            a.fpos=o+wv(x+1./2,z+1./2,y+1./2);wa(a,str);
            a.fpos=o+wv(-x+3./4,z,-y+1./4);wa(a,str);
            a.fpos=o+wv(-z+1./4,-y+3./4,x);wa(a,str);
            a.fpos=o+wv(-z+3./4,y,-x+1./4);wa(a,str);
            a.fpos=o+wv(z,-y+1./4,-x+3./4);wa(a,str);
            a.fpos=o+wv(z+1./2,y+1./2,x+1./2);wa(a,str);
          }
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 229  Im-3m #229
      case 229:
        str.spacegroup="Im-3m";
        str.spacegrouplabel="#229";
        str.spacegroupnumber=229;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,-z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,-x,y);wa(a,str);
          a.fpos=o+wv(-z,x,-y);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y,z,-x);wa(a,str);
          a.fpos=o+wv(y,-z,-x);wa(a,str);
          a.fpos=o+wv(-y,-z,x);wa(a,str);
          a.fpos=o+wv(y,x,-z);wa(a,str);
          a.fpos=o+wv(-y,-x,-z);wa(a,str);
          a.fpos=o+wv(y,-x,z);wa(a,str);
          a.fpos=o+wv(-y,x,z);wa(a,str);
          a.fpos=o+wv(x,z,-y);wa(a,str);
          a.fpos=o+wv(-x,z,y);wa(a,str);
          a.fpos=o+wv(-x,-z,-y);wa(a,str);
          a.fpos=o+wv(x,-z,y);wa(a,str);
          a.fpos=o+wv(z,y,-x);wa(a,str);
          a.fpos=o+wv(z,-y,x);wa(a,str);
          a.fpos=o+wv(-z,y,x);wa(a,str);
          a.fpos=o+wv(-z,-y,-x);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x,y,-z);wa(a,str);
          a.fpos=o+wv(x,-y,z);wa(a,str);
          a.fpos=o+wv(-x,y,z);wa(a,str);
          a.fpos=o+wv(-z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z,x,y);wa(a,str);
          a.fpos=o+wv(z,x,-y);wa(a,str);
          a.fpos=o+wv(z,-x,y);wa(a,str);
          a.fpos=o+wv(-y,-z,-x);wa(a,str);
          a.fpos=o+wv(y,-z,x);wa(a,str);
          a.fpos=o+wv(-y,z,x);wa(a,str);
          a.fpos=o+wv(y,z,-x);wa(a,str);
          a.fpos=o+wv(-y,-x,z);wa(a,str);
          a.fpos=o+wv(y,x,z);wa(a,str);
          a.fpos=o+wv(-y,x,-z);wa(a,str);
          a.fpos=o+wv(y,-x,-z);wa(a,str);
          a.fpos=o+wv(-x,-z,y);wa(a,str);
          a.fpos=o+wv(x,-z,-y);wa(a,str);
          a.fpos=o+wv(x,z,y);wa(a,str);
          a.fpos=o+wv(-x,z,-y);wa(a,str);
          a.fpos=o+wv(-z,-y,x);wa(a,str);
          a.fpos=o+wv(-z,y,-x);wa(a,str);
          a.fpos=o+wv(z,-y,-x);wa(a,str);
          a.fpos=o+wv(z,y,x);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // 230  Ia-3d #230
      case 230:
        str.spacegroup="Ia-3d";
        str.spacegrouplabel="#230";
        str.spacegroupnumber=230;
        str.spacegroupoption="";
        for(uint j=1;j<=2;j++) {
          if(j==1) o=wv(0,0,0);
          if(j==2) o=wv(1./2,1./2,1./2);
          a.fpos=o+wv(x,y,z);wa(a,str);
          a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
          a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
          a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
          a.fpos=o+wv(z,x,y);wa(a,str);
          a.fpos=o+wv(z+1./2,-x+1./2,-y);wa(a,str);
          a.fpos=o+wv(-z+1./2,-x,y+1./2);wa(a,str);
          a.fpos=o+wv(-z,x+1./2,-y+1./2);wa(a,str);
          a.fpos=o+wv(y,z,x);wa(a,str);
          a.fpos=o+wv(-y,z+1./2,-x+1./2);wa(a,str);
          a.fpos=o+wv(y+1./2,-z+1./2,-x);wa(a,str);
          a.fpos=o+wv(-y+1./2,-z,x+1./2);wa(a,str);
          a.fpos=o+wv(y+3./4,x+1./4,-z+1./4);wa(a,str);
          a.fpos=o+wv(-y+3./4,-x+3./4,-z+3./4);wa(a,str);
          a.fpos=o+wv(y+1./4,-x+1./4,z+3./4);wa(a,str);
          a.fpos=o+wv(-y+1./4,x+3./4,z+1./4);wa(a,str);
          a.fpos=o+wv(x+3./4,z+1./4,-y+1./4);wa(a,str);
          a.fpos=o+wv(-x+1./4,z+3./4,y+1./4);wa(a,str);
          a.fpos=o+wv(-x+3./4,-z+3./4,-y+3./4);wa(a,str);
          a.fpos=o+wv(x+1./4,-z+1./4,y+3./4);wa(a,str);
          a.fpos=o+wv(z+3./4,y+1./4,-x+1./4);wa(a,str);
          a.fpos=o+wv(z+1./4,-y+1./4,x+3./4);wa(a,str);
          a.fpos=o+wv(-z+1./4,y+3./4,x+1./4);wa(a,str);
          a.fpos=o+wv(-z+3./4,-y+3./4,-x+3./4);wa(a,str);
          a.fpos=o+wv(-x,-y,-z);wa(a,str);
          a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
          a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
          a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
          a.fpos=o+wv(-z,-x,-y);wa(a,str);
          a.fpos=o+wv(-z+1./2,x+1./2,y);wa(a,str);
          a.fpos=o+wv(z+1./2,x,-y+1./2);wa(a,str);
          a.fpos=o+wv(z,-x+1./2,y+1./2);wa(a,str);
          a.fpos=o+wv(-y,-z,-x);wa(a,str);
          a.fpos=o+wv(y,-z+1./2,x+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./2,z+1./2,x);wa(a,str);
          a.fpos=o+wv(y+1./2,z,-x+1./2);wa(a,str);
          a.fpos=o+wv(-y+1./4,-x+3./4,z+3./4);wa(a,str);
          a.fpos=o+wv(y+1./4,x+1./4,z+1./4);wa(a,str);
          a.fpos=o+wv(-y+3./4,x+3./4,-z+1./4);wa(a,str);
          a.fpos=o+wv(y+3./4,-x+1./4,-z+3./4);wa(a,str);
          a.fpos=o+wv(-x+1./4,-z+3./4,y+3./4);wa(a,str);
          a.fpos=o+wv(x+3./4,-z+1./4,-y+3./4);wa(a,str);
          a.fpos=o+wv(x+1./4,z+1./4,y+1./4);wa(a,str);
          a.fpos=o+wv(-x+3./4,z+3./4,-y+1./4);wa(a,str);
          a.fpos=o+wv(-z+1./4,-y+3./4,x+3./4);wa(a,str);
          a.fpos=o+wv(-z+3./4,y+3./4,-x+1./4);wa(a,str);
          a.fpos=o+wv(z+3./4,-y+1./4,-x+3./4);wa(a,str);
          a.fpos=o+wv(z+1./4,y+1./4,x+1./4);wa(a,str);
        }
        break;
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
      // ----------------------------------------------------------------------------
    }
  }
  // cycl atoms
  // cerr << "DONE atoms=" << str.atoms.size() << endl;
  //DONE
  str.MakeBasis();
  return str;
}

//----------------------------------------------------------------------------

#endif//_WYCKOFF_IMPLEMENTATIONS_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************


