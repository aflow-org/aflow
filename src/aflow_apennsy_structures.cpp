// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

#ifndef _AFLOW_APENNSY_STRUCTURES_CPP_
#define _AFLOW_APENNSY_STRUCTURES_CPP_
#include "aflow.h"

// **************************************************************************
// Load Structures types HTQC or GUS
// **************************************************************************
namespace aflowlib {
  bool _aflowlib_entry::FixDescription(void) {
    pureA=FALSE,pureB=FALSE;
    fcc=FALSE;bcc=FALSE;hcp=FALSE;
    string alt_structure_description;
    if(structure_name=="0") {structure_description="nnn";alt_structure_description="nnn";return TRUE;}
  
    // LIB1 AND UNITARY
    if(structure_name=="A1") {structure_description="A1";alt_structure_description="A1 FCC Cu";fcc=TRUE;pureA=TRUE;return TRUE;}
    if(structure_name=="A2") {structure_description="A2";alt_structure_description="A2 BCC W" ;bcc=TRUE;pureA=TRUE;return TRUE;}
    if(structure_name=="A3") {structure_description="A3";alt_structure_description="A3 HCP Mg";hcp=TRUE;pureA=TRUE;return TRUE;}

    // LIB1 BUT BINARY
    if(structure_name=="1") {structure_description="A1";alt_structure_description="A1 FCC Cu";fcc=TRUE;pureB=TRUE;return TRUE;}
    if(structure_name=="f2") {structure_description="A1";alt_structure_description="A1 FCC Cu";fcc=TRUE;pureB=TRUE;return TRUE;}
    if(structure_name=="2") {structure_description="A1";alt_structure_description="A1 FCC Cu";fcc=TRUE;pureA=TRUE;return TRUE;}
    if(structure_name=="f1") {structure_description="A1";alt_structure_description="A1 FCC Cu";fcc=TRUE;pureA=TRUE;return TRUE;}
    if(structure_name=="1.cub") {structure_description="A1.cub";alt_structure_description="A1 FCC Cu";fcc=TRUE;pureB=TRUE;return TRUE;}
    if(structure_name=="2.cub") {structure_description="A1.cub";alt_structure_description="A1 FCC Cu";fcc=TRUE;pureA=TRUE;return TRUE;}
    if(structure_name=="58") {structure_description="A2";alt_structure_description="A2 BCC W";bcc=TRUE;pureB=TRUE;return TRUE;}
    if(structure_name=="b2") {structure_description="A2";alt_structure_description="A2 BCC W";bcc=TRUE;pureB=TRUE;return TRUE;}
    if(structure_name=="59") {structure_description="A2";alt_structure_description="A2 BCC W";bcc=TRUE;pureA=TRUE;return TRUE;}
    if(structure_name=="b1") {structure_description="A2";alt_structure_description="A2 BCC W";bcc=TRUE;pureA=TRUE;return TRUE;}
    if(structure_name=="58.cub") {structure_description="A2.cub";alt_structure_description="A2 BCC W";bcc=TRUE;pureB=TRUE;return TRUE;}
    if(structure_name=="59.cub") {structure_description="A2.cub";alt_structure_description="A2 BCC W";bcc=TRUE;pureA=TRUE;return TRUE;}
    if(structure_name=="115") {structure_description="A3";alt_structure_description="A3 HCP Mg";hcp=TRUE;pureB=TRUE; return TRUE;}
    if(structure_name=="h2") {structure_description="A3";alt_structure_description="A3 HCP Mg";hcp=TRUE;pureB=TRUE; return TRUE;}
    if(structure_name=="117") {structure_description="A3";alt_structure_description="A3 HCP Mg";hcp=TRUE;pureA=TRUE;}
    if(structure_name=="h1") {structure_description="A3";alt_structure_description="A3 HCP Mg";hcp=TRUE;pureA=TRUE;}

    if(structure_name=="301") {structure_description="A4 DIA C";pureA=TRUE;return TRUE;}
    if(structure_name=="302") {structure_description="A4 DIA C";pureB=TRUE;return TRUE;}

    if(structure_name=="303") {structure_description="A6 TET In";pureA=TRUE;return TRUE;}
    if(structure_name=="304") {structure_description="A6 TET In";pureB=TRUE;return TRUE;}

    if(structure_name=="305") {structure_description="A5 \\beta-Sn";pureA=TRUE;return TRUE;}
    if(structure_name=="306") {structure_description="A5 \\beta-Sn";pureB=TRUE;return TRUE;}
  
    if(structure_name=="307") {structure_description="A7 RHL \\alpha-As hP6";pureA=TRUE;return TRUE;}
    if(structure_name=="308") {structure_description="A7 RHL \\alpha-As hP6";pureB=TRUE;return TRUE;}
  
    if(structure_name=="A7.A") {structure_description="A7 RHL \\alpha-As hR2";pureA=TRUE;return TRUE;}
    if(structure_name=="A7.B") {structure_description="A7 RHL \\alpha-As hR2";pureB=TRUE;return TRUE;}
  
    if(structure_name=="323") {structure_description="A8 \\gamma-Se";pureA=TRUE;return TRUE;}
    if(structure_name=="324") {structure_description="A8 \\gamma-Se";pureB=TRUE;return TRUE;}
  
    if(structure_name=="325") {structure_description="A9 graphite";pureA=TRUE;return TRUE;}
    if(structure_name=="326") {structure_description="A9 graphite";pureB=TRUE;return TRUE;}
  
    if(structure_name=="317") {structure_description="A11 \\alpha-Ga";pureA=TRUE;return TRUE;}
    if(structure_name=="318") {structure_description="A11 \\alpha-Ga";pureB=TRUE;return TRUE;}
  
    if(structure_name=="319") {structure_description="A12 \\alpha-Mn";pureA=TRUE;return TRUE;}
    if(structure_name=="320") {structure_description="A12 \\alpha-Mn";pureB=TRUE;return TRUE;}
  
    if(structure_name=="321") {structure_description="A13 \\beta-Mn";pureA=TRUE;return TRUE;}
    if(structure_name=="322") {structure_description="A13 \\beta-Mn";pureB=TRUE;return TRUE;}
  
    if(structure_name=="327") {structure_description="Aa \\alpha-Pa";pureA=TRUE;return TRUE;}
    if(structure_name=="328") {structure_description="Aa \\alpha-Pa";pureB=TRUE;return TRUE;}
  
    if(structure_name=="329") {structure_description="Ab \\beta-U";pureA=TRUE;return TRUE;}
    if(structure_name=="330") {structure_description="Ab \\beta-U";pureB=TRUE;return TRUE;}
  
    if(structure_name=="331") {structure_description="Ac \\alpha-Np";pureA=TRUE;return TRUE;}
    if(structure_name=="332") {structure_description="Ac \\alpha-Np";pureB=TRUE;return TRUE;}
  
    if(structure_name=="333") {structure_description="Ad \\beta-Np";pureA=TRUE;return TRUE;}
    if(structure_name=="334") {structure_description="Ad \\beta-Np";pureB=TRUE;return TRUE;}
  
    // fcc
    if(structure_name=="3") {structure_description="L1_{0} AuCu";fcc=TRUE;return TRUE;}
    if(structure_name=="4") {structure_description="L1_{1} CuPt";fcc=TRUE;return TRUE;}
    if(structure_name=="5") {structure_description="\\beta_{1} FCC_{AB2}^{[100]}";fcc=TRUE;return TRUE;}
    if(structure_name=="6") {structure_description="\\beta_{2} FCC_{AB2}^{[100]}";fcc=TRUE;return TRUE;}
    if(structure_name=="7") {structure_description="MoPt_{2}";fcc=TRUE;return TRUE;}
    if(structure_name=="8") {structure_description="MoPt_{2}";fcc=TRUE;return TRUE;}
    if(structure_name=="9") {structure_description="\\alpha_{1} FCC_{AB2}^{[111]}";fcc=TRUE;return TRUE;}
    if(structure_name=="10") {structure_description="\\alpha_{2} FCC_{AB2}^{[111]}";fcc=TRUE;return TRUE;}
    if(structure_name=="11") {structure_description="M(CdPt)";fcc=TRUE;return TRUE;}
    if(structure_name=="12") {structure_description="M(CdPt)";fcc=TRUE;return TRUE;}
    if(structure_name=="13") {structure_description="M(PtTc) Z1 FCC_{AB3}^{[001]}";fcc=TRUE;return TRUE;}
    if(structure_name=="14") {structure_description="Z2 FCC_{A2B2}^{[001]}";fcc=TRUE;return TRUE;}
    if(structure_name=="15") {structure_description="M(PtTc) Z3 FCC_{AB3}^{[001]}";fcc=TRUE;return TRUE;}
    if(structure_name=="16") {structure_description="FCC_{AB3}^{[113]}";fcc=TRUE;return TRUE;}
    if(structure_name=="17") {structure_description="W2 FCC_{A2B2}^{[311]}";fcc=TRUE;return TRUE;}
    if(structure_name=="18") {structure_description="FCC_{AB3}^{[113]}";fcc=TRUE;return TRUE;}
    if(structure_name=="19") {structure_description="Y1 FCC_{AB3}^{[011]}";fcc=TRUE;return TRUE;}
    if(structure_name=="20") {structure_description="M(AgZr) Y2 FCC_{A2B2}^{[011]}";fcc=TRUE;return TRUE;}
    if(structure_name=="21") {structure_description="Y3 FCC_{AB3}^{[011]}";fcc=TRUE;return TRUE;}
    if(structure_name=="22") {structure_description="D0_{22} Al3Ti";fcc=TRUE;return TRUE;}
    if(structure_name=="23") {structure_description="CH 40 NbAs NbP";fcc=TRUE;return TRUE;}
    if(structure_name=="24") {structure_description="D0_{22} Al3Ti";fcc=TRUE;return TRUE;}
    if(structure_name=="25") {structure_description="L1_{2} AuCu3";fcc=TRUE;return TRUE;}
    if(structure_name=="26") {structure_description="L1_{2} AuCu3";fcc=TRUE;return TRUE;}
    if(structure_name=="27") {structure_description="V1 FCC_{AB3}^{[111]}";fcc=TRUE;return TRUE;}
    if(structure_name=="28") {structure_description="V2 FCC_{A2B2}^{[111]}";fcc=TRUE;return TRUE;}
    if(structure_name=="29") {structure_description="V3 FCC_{AB3}^{[111]}";fcc=TRUE;return TRUE;}
    //  bcc
    if(structure_name=="60") {structure_description="\\gamma-IrV";bcc=TRUE;return TRUE;}
    if(structure_name=="61") {structure_description="B2 CsCl";bcc=TRUE;return TRUE;}
    if(structure_name=="62") {structure_description="BCC_{AB2}^{[211]}";bcc=TRUE;return TRUE;}
    if(structure_name=="63") {structure_description="BCC_{AB2}^{[211]}";bcc=TRUE;return TRUE;}
    // if(structure_name=="64") {structure_description="BCC_{AB2}^{[011]}";bcc=TRUE;return TRUE;}
    // if(structure_name=="65") {structure_description="BCC_{AB2}^{[011]}";bcc=TRUE;return TRUE;}
    if(structure_name=="66") {structure_description="C11_b MoSi_{2} BCC_{AB2}^{[001]}";bcc=TRUE;return TRUE;}
    if(structure_name=="67") {structure_description="C11_b MoSi_{2} BCC_{AB2}^{[001]}";bcc=TRUE;return TRUE;}
    if(structure_name=="68") {bcc=TRUE;return TRUE;}
    if(structure_name=="69") {bcc=TRUE;return TRUE;}
    if(structure_name=="70") {bcc=TRUE;return TRUE;}
    if(structure_name=="71") {bcc=TRUE;return TRUE;}
    if(structure_name=="72") {bcc=TRUE;return TRUE;}
    if(structure_name=="73") {bcc=TRUE;return TRUE;}
    if(structure_name=="74") {bcc=TRUE;return TRUE;}
    if(structure_name=="75") {bcc=TRUE;return TRUE;}
    if(structure_name=="76") {structure_description="BCC_{AB3}^{[001]}";bcc=TRUE;return TRUE;}
    if(structure_name=="77") {structure_description="B11 \\gamma-CuTi BCC_{A2B2}^{[001]}";bcc=TRUE;return TRUE;}
    if(structure_name=="78") {structure_description="BCC_{AB3}^{[001]}";bcc=TRUE;return TRUE;}
    if(structure_name=="79") {bcc=TRUE;return TRUE;}
    if(structure_name=="80") {bcc=TRUE;return TRUE;}
    if(structure_name=="81") {bcc=TRUE;return TRUE;}
    if(structure_name=="82") {structure_description="L6_0";bcc=TRUE;return TRUE;}
    if(structure_name=="83") {structure_description="L6_0";bcc=TRUE;return TRUE;}
    if(structure_name=="84") {structure_description="D0_{3} AlFe_{3}";bcc=TRUE;return TRUE;}
    if(structure_name=="85") {structure_description="B32 NaTl";bcc=TRUE;return TRUE;}
    if(structure_name=="86") {structure_description="D0_{3} AlFe_{3}";bcc=TRUE;return TRUE;}
    // hcp
    if(structure_name=="116") {structure_description="B_h WC";hcp=TRUE;return TRUE;}
    if(structure_name=="118") {hcp=TRUE;return TRUE;}
    if(structure_name=="121") {hcp=TRUE;return TRUE;}
    if(structure_name=="119") {structure_description="oP4 CuTe";hcp=TRUE;return TRUE;}
    if(structure_name=="120") {structure_description="B19 AuCd";hcp=TRUE;return TRUE;}
    if(structure_name=="123") {hcp=TRUE;return TRUE;}
    if(structure_name=="122") {hcp=TRUE;return TRUE;}
    if(structure_name=="124") {hcp=TRUE;return TRUE;}
    if(structure_name=="125") {hcp=TRUE;return TRUE;}
    if(structure_name=="127") {hcp=TRUE;return TRUE;}
    if(structure_name=="126") {hcp=TRUE;return TRUE;}
    if(structure_name=="128") {hcp=TRUE;return TRUE;}
    if(structure_name=="132") {hcp=TRUE;return TRUE;}
    if(structure_name=="129") {hcp=TRUE;return TRUE;}
    if(structure_name=="134") {hcp=TRUE;return TRUE;}
    if(structure_name=="130") {hcp=TRUE;return TRUE;}
    if(structure_name=="137") {hcp=TRUE;return TRUE;}
    if(structure_name=="131") {hcp=TRUE;return TRUE;}
    if(structure_name=="135") {hcp=TRUE;return TRUE;}
    if(structure_name=="133") {hcp=TRUE;return TRUE;}
    if(structure_name=="140") {hcp=TRUE;return TRUE;}
    if(structure_name=="136") {hcp=TRUE;return TRUE;}
    if(structure_name=="138") {hcp=TRUE;return TRUE;}
    if(structure_name=="139") {hcp=TRUE;return TRUE;}
    if(structure_name=="141") {hcp=TRUE;return TRUE;}
    if(structure_name=="145") {hcp=TRUE;return TRUE;}
    if(structure_name=="142") {hcp=TRUE;return TRUE;}
    if(structure_name=="147") {hcp=TRUE;return TRUE;}
    if(structure_name=="143") {structure_description="ZrSi_{2} OC12";hcp=TRUE;return TRUE;}
    if(structure_name=="150") {structure_description="ZrSi_{2} OC12";hcp=TRUE;return TRUE;}
    //
    // LAVES OPMEGAS AND OTHERSreturn TRUE;}
    if(structure_name=="178") {structure_description="C14 MgZn_{2}";return TRUE;}
    if(structure_name=="179") {structure_description="C14 MgZn_{2}";return TRUE;}
    if(structure_name=="180") {structure_description="D0_{11} Fe_{3}C";return TRUE;}
    if(structure_name=="181") {structure_description="D0_{11} Fe_{3}C";return TRUE;}
    if(structure_name=="182") {structure_description="C15 Cu_{2}Mg";return TRUE;}
    if(structure_name=="183") {structure_description="C15 Cu_{2}Mg";return TRUE;}
    if(structure_name=="184") {structure_description="A15 Cr_{3}Si";return TRUE;}
    if(structure_name=="185") {structure_description="A15 Cr_{3}Si";return TRUE;}
    if(structure_name=="186") {structure_description="D0_{19} Ni_{3}Sn";return TRUE;}
    if(structure_name=="187") {structure_description="D0_{19} Ni_{3}Sn";return TRUE;}
    if(structure_name=="188") {structure_description="C49 ZrSi_{2}";return TRUE;}
    if(structure_name=="189") {structure_description="C49 ZrSi_{2}";return TRUE;}
    if(structure_name=="190") {structure_description="\\Omega Z=14";return TRUE;} // "C6 CdI2"
    if(structure_name=="191") {structure_description="\\Omega Z=14";return TRUE;} // "C6 CdI2"
    if(structure_name=="192") {structure_description="B33 CrB";return TRUE;}
    if(structure_name=="193") {structure_description="B33 CrB";return TRUE;}
    if(structure_name=="194") {structure_description="TlI";return TRUE;}
    if(structure_name=="195") {structure_description="B20 FeSi";return TRUE;}
    if(structure_name=="196") {structure_description="CdTi";return TRUE;}
    if(structure_name=="197") {structure_description="CdTi";return TRUE;}
    if(structure_name=="198") {structure_description="InTh";return TRUE;}
    if(structure_name=="199") {structure_description="SeTl";return TRUE;}
    if(structure_name=="201") {structure_description="B1 NaCl";return TRUE;}
    if(structure_name=="202") {structure_description="D1_{3} Al_{4}Ba";return TRUE;}
    if(structure_name=="203") {structure_description="D1_{3} Al_{4}Ba";return TRUE;}
    if(structure_name=="204") {structure_description="D2_d CaCu_{5}";return TRUE;}
    if(structure_name=="205") {structure_description="D2_d CaCu_{5}";return TRUE;}
    if(structure_name=="206") {structure_description="D2_b ThMn12";return TRUE;}
    if(structure_name=="207") {structure_description="D2_b ThMn12";return TRUE;}
    if(structure_name=="208") {structure_description="C22 Fe_{2}P";return TRUE;}
    if(structure_name=="209") {structure_description="C22 Fe_{2}P";return TRUE;}
    if(structure_name=="210") {structure_description="C37 Co_{2}Si";return TRUE;}
    if(structure_name=="211") {structure_description="C37 Co_{2}Si";return TRUE;}
    if(structure_name=="212") {structure_description="Th_{2}Zn_{17}";return TRUE;}
    if(structure_name=="213") {structure_description="Th_{2}Zn_{17}";return TRUE;}
    if(structure_name=="214") {structure_description="Th_{3}P_{4}";return TRUE;}
    if(structure_name=="215") {structure_description="Th_{3}P_{4}";return TRUE;}
    if(structure_name=="216") {structure_description="C32 AlB_{2}";return TRUE;}
    if(structure_name=="217") {structure_description="C32 AlB_{2}";return TRUE;}
    if(structure_name=="218") {structure_description="B3 ZnS zincblende";return TRUE;}
    if(structure_name=="219") {structure_description="B4 ZnS Wurtzite";return TRUE;}
    if(structure_name=="220") {structure_description="B8_1 NiAs";return TRUE;}
    if(structure_name=="221") {structure_description="B8_1 NiAs";return TRUE;}
    if(structure_name=="222") {structure_description="D8_{8} Mn_{5}Si_{3}";return TRUE;}
    if(structure_name=="223") {structure_description="D8_{8} Mn_{5}Si_{3}";return TRUE;}
    if(structure_name=="224") {structure_description="Th_{2}Ni_{17}";return TRUE;}
    if(structure_name=="225") {structure_description="Th_{2}Ni_{17}";return TRUE;}
    if(structure_name=="226") {structure_description="A15 Cr_{3}Si";return TRUE;}
    if(structure_name=="227") {structure_description="A15 Cr_{3}Si";return TRUE;}
    if(structure_name=="228") {structure_description="C38 Cu_{2}Sb";return TRUE;}
    if(structure_name=="229") {structure_description="C38 Cu_{2}Sb";return TRUE;}
    if(structure_name=="230") {structure_description="C16 Al_{2}Cu";return TRUE;}
    if(structure_name=="231") {structure_description="C16 Al_{2}Cu";return TRUE;}
    if(structure_name=="232") {structure_description="C11_b MoSi_{2} BCC_{AB2}^{[001]}";return TRUE;}
    if(structure_name=="233") {structure_description="C11_b MoSi_{2} BCC_{AB2}^{[001]}";return TRUE;}
    if(structure_name=="234") {structure_description="C16 Al_{2}Cu";return TRUE;}
    if(structure_name=="235") {structure_description="C16 Al_{2}Cu";return TRUE;}
    if(structure_name=="236") {structure_description="hP142 Cd58Y13";return TRUE;}
    if(structure_name=="237") {structure_description="hP142 Cd58Y13";return TRUE;}
    if(structure_name=="238") {structure_description="hp24 PuAl_{3} - Co_{3}V";return TRUE;}
    if(structure_name=="239") {structure_description="hp24 PuAl_{3} - Co_{3}V";return TRUE;}
    // 240,241 ??
    if(structure_name=="242") {structure_description="oP24 NbPd_{3}";return TRUE;}
    if(structure_name=="243") {structure_description="oP24 NbPd_{3}";return TRUE;}
    if(structure_name=="244") {structure_description="D0_{24} Ni_{3}Ti";return TRUE;}
    if(structure_name=="245") {structure_description="D0_{24} Ni_{3}Ti";return TRUE;}
    if(structure_name=="246") {structure_description="hP6 CaIn_{2}";return TRUE;}
    if(structure_name=="247") {structure_description="hP6 CaIn_{2}";return TRUE;}
    if(structure_name=="248") {structure_description="C_c ThSi_{2}";return TRUE;}
    if(structure_name=="249") {structure_description="C_c ThSi_{2}";return TRUE;}
    //if(structure_name=="250") {structure_description="C11_a CaC_{2}";return TRUE;}
    //if(structure_name=="251") {structure_description="C11_a CaC_{2}";return TRUE;}
    if(structure_name=="252") {structure_description="C15_b AuBe_{5}";return TRUE;}
    if(structure_name=="253") {structure_description="C15_b AuBe_{5}";return TRUE;}
    if(structure_name=="254") {structure_description="W_{5}Si_{3}";return TRUE;}
    if(structure_name=="255") {structure_description="W_{5}Si_{3}";return TRUE;}
    if(structure_name=="256") {structure_description="B27 BFe";return TRUE;}
    if(structure_name=="257") {structure_description="C18 marcasite FeS_{2}";return TRUE;}
    if(structure_name=="258") {structure_description="C18 marcasite FeS_{2}";return TRUE;}
    if(structure_name=="259") {structure_description="Bi_{2}Te_{3}";return TRUE;}
    if(structure_name=="260") {structure_description="Bi_{2}Te_{3}";return TRUE;}
    if(structure_name=="261") {structure_description="NiTi_{2}";return TRUE;}
    if(structure_name=="262") {structure_description="NiTi_{2}";return TRUE;}
    if(structure_name=="263") {structure_description="Cu_{4}Ti_{3}";return TRUE;}
    if(structure_name=="264") {structure_description="Cu_{4}Ti_{3}";return TRUE;}
    // 265,266
    if(structure_name=="267") {structure_description="C18 marcasite FeS_{2}";return TRUE;}
    if(structure_name=="268") {structure_description="C18 marcasite FeS_{2}";return TRUE;}
    if(structure_name=="269") {structure_description="C6 CdI_{2}";return TRUE;}
    if(structure_name=="270") {structure_description="C6 CdI_{2}";return TRUE;}
    if(structure_name=="271") {structure_description="Cd_{3}Y";return TRUE;}
    if(structure_name=="272") {structure_description="Cd_{3}Y";return TRUE;}
    if(structure_name=="273") {structure_description="B8_{2}";return TRUE;}
    if(structure_name=="274") {structure_description="B8_{2}";return TRUE;}
    if(structure_name=="273p") {structure_description="Fe_{1.5}Ge_{1}";return TRUE;}
    if(structure_name=="274p") {structure_description="Fe_{1.5}Ge_{1}";return TRUE;}
    if(structure_name=="275") {structure_description="Ni_{2}Si";return TRUE;}
    if(structure_name=="276") {structure_description="Ni_{2}Si";return TRUE;}
    if(structure_name=="277") {structure_description="D0_a \\beta-Cu_{3}Ti";return TRUE;}
    if(structure_name=="278") {structure_description="D0_a \\beta-Cu_{3}Ti";return TRUE;}
    if(structure_name=="279") {structure_description="D0_{23} Al_{3}Zr";return TRUE;}
    if(structure_name=="280") {structure_description="D0_{23} Al_{3}Zr";return TRUE;}
    if(structure_name=="281") {structure_description="C2 pyrite FeS_{2} ";return TRUE;}
    if(structure_name=="282") {structure_description="C2 pyrite FeS_{2} ";return TRUE;}
    if(structure_name=="283") {structure_description="GdSi-2(1.)";return TRUE;}
    if(structure_name=="284") {structure_description="GdSi-2(1.)";return TRUE;}
    if(structure_name=="285") {structure_description="D1_a MoNi_{4}";return TRUE;}
    if(structure_name=="286") {structure_description="D1_a MoNi_{4}";return TRUE;}
    if(structure_name=="287") {structure_description="CuZr_{2}";return TRUE;}
    if(structure_name=="288") {structure_description="CuZr_{2}";return TRUE;}
    if(structure_name=="289") {structure_description="D0_{9} \\alpha-ReO_{3}";return TRUE;}
    if(structure_name=="290") {structure_description="D0_{9} \\alpha-ReO_{3}";return TRUE;}
    if(structure_name=="291") {structure_description="B10 PbO PbS";return TRUE;}
    if(structure_name=="292") {structure_description="B10 PbO PbS";return TRUE;}
    // **************************************************************************
    // Al Mg CALPHAD paper
    //concentrations.pbZERO((double);if(structure_name=="BiMg.AL1.02_03") {structure_description="BiMg.AL1.02-03";return TRUE;}
    //concentrations.pbZERO((double);if(structure_name=="MgIn.AL1.06_03") {structure_description="MgIn.AL1.06-03";return TRUE;}
    //concentrations.pbZERO((double);if(structure_name=="MgIn.AL1.20_08") {structure_description="MgIn.AL1.20-08";return TRUE;}
    //concentrations.pbZERO((double);if(structure_name=="MgSb.AL1.03_02") {structure_description="MgSb.AL1.03-02";return TRUE;}
    //concentrations.pbZERO((double);if(structure_name=="BiIn.AL1.02_04") {structure_description="BiIn.AL1.02-04";return TRUE;}
    //concentrations.pbZERO((double);if(structure_name=="BiIn.AL1.12_20") {structure_description="BiIn.AL1.12-20";return TRUE;}
    // **************************************************************************
    // Check Nstructures
    if(structure_name=="309") {structure_description="Ca_{7}Ge";return TRUE;}
    if(structure_name=="310") {structure_description="Ca_{7}Ge";return TRUE;}
    if(structure_name=="311") {structure_description="Pt_{8}Ti";return TRUE;}
    if(structure_name=="312") {structure_description="Pt_{8}Ti";return TRUE;}
    if(structure_name=="313") {structure_description="V_{4}Zn_{5}";return TRUE;}
    if(structure_name=="314") {structure_description="V_{4}Zn_{5}";return TRUE;}
    if(structure_name=="315") {structure_description="C36 MgNi_{2}";return TRUE;}
    if(structure_name=="316") {structure_description="C36 MgNi_{2}";return TRUE;}
    // 327 ??
    if(structure_name=="359") {structure_description="Ga_{4}Ti_{5}";return TRUE;}
    if(structure_name=="360") {structure_description="Ga_{4}Ti_{5}";return TRUE;}
    if(structure_name=="361") {structure_description="Al_{2}Zr_{3}";return TRUE;}
    if(structure_name=="362") {structure_description="Al_{2}Zr_{3}";return TRUE;}
    if(structure_name=="363") {structure_description="Al_{3}Zr_{4}";return TRUE;}
    if(structure_name=="364") {structure_description="Al_{3}Zr_{4}";return TRUE;}
    if(structure_name=="365") {structure_description="Al_{3}Zr_{2}";return TRUE;}
    if(structure_name=="366") {structure_description="Al_{3}Zr_{2}";return TRUE;}
    if(structure_name=="367") {structure_description="NaZn_{13}";return TRUE;}
    if(structure_name=="368") {structure_description="NaZn_{13}";return TRUE;}
    if(structure_name=="369") {structure_description="CaB_{6}";return TRUE;}
    if(structure_name=="370") {structure_description="CaB_{6}";return TRUE;}
    if(structure_name=="371") {structure_description="D5_{19} Al_{3}Ni_{2}";return TRUE;}
    if(structure_name=="372") {structure_description="D5_{19} Al_{3}Ni_{2}";return TRUE;}
    if(structure_name=="373") {structure_description="Ga_{4}Ni";return TRUE;}
    if(structure_name=="374") {structure_description="Ga_{4}Ni";return TRUE;}
    if(structure_name=="375") {structure_description="Ga_{3}Pt_{5}";return TRUE;}
    if(structure_name=="376") {structure_description="Ga_{3}Pt_{5}";return TRUE;}
    if(structure_name=="377") {structure_description="Ga_{9}Ni_{13}";return TRUE;}
    if(structure_name=="378") {structure_description="Ga_{9}Ni_{13}";return TRUE;}
    if(structure_name=="379") {structure_description="Ga_{4}Ni_{3}";return TRUE;}
    if(structure_name=="380") {structure_description="Ga_{4}Ni_{3}";return TRUE;}
    if(structure_name=="381") {structure_description="Ga_{2}Hf";return TRUE;}
    if(structure_name=="382") {structure_description="Ga_{2}Hf";return TRUE;}
    if(structure_name=="383") {structure_description="C_{6}Mn_{23}";return TRUE;}
    if(structure_name=="384") {structure_description="C_{6}Mn_{23}";return TRUE;}
    if(structure_name=="385") {structure_description="Cu_{8}Hf_{3}";return TRUE;}
    if(structure_name=="386") {structure_description="Cu_{8}Hf_{3}";return TRUE;}
    if(structure_name=="387") {structure_description="Ni_{10}Zr_{7}";return TRUE;}
    if(structure_name=="388") {structure_description="Ni_{10}Zr_{7}";return TRUE;}
    if(structure_name=="389") {structure_description="Hg_{2}Pt";return TRUE;}
    if(structure_name=="390") {structure_description="Hg_{2}Pt";return TRUE;}
    if(structure_name=="391") {structure_description="MgB_{4}";return TRUE;}
    if(structure_name=="392") {structure_description="MgB_{4}";return TRUE;}
    if(structure_name=="393") {structure_description="MgB_{7}";return TRUE;}
    if(structure_name=="394") {structure_description="MgB_{7}";return TRUE;}
    if(structure_name=="395") {structure_description="Ge_{3}Rh_{5}";return TRUE;}
    if(structure_name=="396") {structure_description="Ge_{3}Rh_{5}";return TRUE;}
    if(structure_name=="397") {structure_description="In_{4}Ti_{3}";return TRUE;}
    if(structure_name=="398") {structure_description="In_{4}Ti_{3}";return TRUE;}
    if(structure_name=="399") {structure_description="V_{7.49}Sb_{9}";return TRUE;}
    if(structure_name=="400") {structure_description="V_{7.49}Sb_{9}";return TRUE;}
    if(structure_name=="401") {structure_description="Li_{3}B_{14}";return TRUE;}
    if(structure_name=="402") {structure_description="Li_{3}B_{14}";return TRUE;}
    if(structure_name=="403") {structure_description="LiB_{3}";return TRUE;}
    if(structure_name=="404") {structure_description="LiB_{3}";return TRUE;}
    if(structure_name=="405") {structure_description="LiB-MS1";return TRUE;}
    if(structure_name=="406") {structure_description="LiB-MS1";return TRUE;}
    if(structure_name=="407") {structure_description="LiB-MS2";return TRUE;}
    if(structure_name=="408") {structure_description="LiB-MS2";return TRUE;}
    if(structure_name=="411") {structure_description="Re_{24}Ti_{5}";return TRUE;}
    if(structure_name=="412") {structure_description="Re_{24}Ti_{5}";return TRUE;}
    if(structure_name=="413") {structure_description="Ni_{7}Zr_{2}";return TRUE;}
    if(structure_name=="414") {structure_description="Ni_{7}Zr_{2}";return TRUE;}
    if(structure_name=="415") {structure_description="Zr_{9}Pt11";return TRUE;}
    if(structure_name=="416") {structure_description="Zr_{9}Pt11";return TRUE;}
    if(structure_name=="417") {structure_description="Au_{4}Zr";return TRUE;}
    if(structure_name=="418") {structure_description="Au_{4}Zr";return TRUE;}
    if(structure_name=="419") {structure_description="Hf_{54}Os_{17}";return TRUE;}
    if(structure_name=="420") {structure_description="Hf_{54}Os_{17}";return TRUE;}
    if(structure_name=="421") {structure_description="C40 CrSi_{2}";return TRUE;}
    if(structure_name=="422") {structure_description="C40 CrSi_{2}";return TRUE;}
    if(structure_name=="423.t0") {structure_description="Ag_{51}Gd_{14}.t0";return TRUE;}
    if(structure_name=="424.t0") {structure_description="Ag_{51}Gd_{14}.t0";return TRUE;}
    if(structure_name=="423.t1") {structure_description="Ag_{51}Gd_{14}.t1";return TRUE;}
    if(structure_name=="424.t1") {structure_description="Ag_{51}Gd_{14}.t1";return TRUE;}
    if(structure_name=="423.t2") {structure_description="Ag_{51}Gd_{14}.t2";return TRUE;}
    if(structure_name=="424.t2") {structure_description="Ag_{51}Gd_{14}.t2";return TRUE;}
    if(structure_name=="425") {structure_description="Zn_{22}Zr";return TRUE;}
    if(structure_name=="426") {structure_description="Zn_{22}Zr";return TRUE;}
    if(structure_name=="427") {structure_description="As_{2}Ti";return TRUE;}
    if(structure_name=="428") {structure_description="As_{2}Ti";return TRUE;}
    if(structure_name=="429") {structure_description="Ge_{10}Ho_{11}";return TRUE;}
    if(structure_name=="430") {structure_description="Ge_{10}Ho_{11}";return TRUE;}
    if(structure_name=="431") {structure_description="Ir_{3}Zr_{5}";return TRUE;}
    if(structure_name=="432") {structure_description="Ir_{3}Zr_{5}";return TRUE;}
    if(structure_name=="433") {structure_description="CaCl_{2}";return TRUE;}
    if(structure_name=="434") {structure_description="CaCl_{2}";return TRUE;}
    if(structure_name=="435") {structure_description="CFe_{4}";return TRUE;}
    if(structure_name=="436") {structure_description="CFe_{4}";return TRUE;}
    if(structure_name=="437") {structure_description="C_{2}Mn_{5}";return TRUE;}
    if(structure_name=="438") {structure_description="C_{2}Mn_{5}";return TRUE;}
    if(structure_name=="439") {structure_description="C_{3}Mn_{7}";return TRUE;}
    if(structure_name=="440") {structure_description="C_{3}Mn_{7}";return TRUE;}
    if(structure_name=="441") {structure_description="Fe_{3}Th_{7}";return TRUE;}
    if(structure_name=="442") {structure_description="Fe_{3}Th_{7}";return TRUE;}
    if(structure_name=="443") {structure_description="F_{3}Fe-#167";return TRUE;}
    if(structure_name=="444") {structure_description="F_{3}Fe-#167";return TRUE;}
    if(structure_name=="445") {structure_description="F_{3}Fe-#150";return TRUE;}
    if(structure_name=="446") {structure_description="F_{3}Fe-#150";return TRUE;}
    if(structure_name=="447") {structure_description="NiTi";return TRUE;}
    if(structure_name=="448") {structure_description="NiTi";return TRUE;}
    if(structure_name=="449") {structure_description="C3 Ag_2O";return TRUE;}
    if(structure_name=="450") {structure_description="C3 Ag_2O";return TRUE;}
    if(structure_name=="451") {structure_description="Bb \\eta-Ag_{2}Zn";return TRUE;}
    if(structure_name=="452") {structure_description="Bb \\eta-Ag_{2}Zn";return TRUE;}
    if(structure_name=="453") {structure_description="D0_{15} AlCl_{3}";return TRUE;}
    if(structure_name=="454") {structure_description="D0_{15} AlCl_{3}";return TRUE;}
    if(structure_name=="455") {structure_description="D0_{14} AlF_{3}";return TRUE;}
    if(structure_name=="456") {structure_description="D0_{14} AlF_{3}";return TRUE;}
    if(structure_name=="457") {structure_description="D5_1 \\alpha-Al_{2}O_{3}";return TRUE;}
    if(structure_name=="458") {structure_description="D5_1 \\alpha-Al_{2}O_{3}";return TRUE;}
    if(structure_name=="459") {structure_description="Al_{12}W";return TRUE;}
    if(structure_name=="460") {structure_description="Al_{12}W";return TRUE;}
    if(structure_name=="461") {structure_description="Ag_{2}Se";return TRUE;}
    if(structure_name=="462") {structure_description="Ag_{2}Se";return TRUE;}
    if(structure_name=="463") {structure_description="AlPd";return TRUE;}
    if(structure_name=="464") {structure_description="AlPd";return TRUE;}
    if(structure_name=="465") {structure_description="AlSb";return TRUE;}
    if(structure_name=="466") {structure_description="AlSb";return TRUE;}
    if(structure_name=="467") {structure_description="Cu_{5}Zn_{8} \\gamma-Brass";return TRUE;}
    if(structure_name=="468") {structure_description="Cu_{5}Zn_{8} \\gamma-Brass";return TRUE;}
    if(structure_name=="469") {structure_description="BaPb_{3}";return TRUE;}
    if(structure_name=="470") {structure_description="BaPb_{3}";return TRUE;}
    if(structure_name=="471") {structure_description="Hf_{3}Sc*-h321";return TRUE;}
    if(structure_name=="472") {structure_description="Hf_{3}Sc*-h321";return TRUE;}
    if(structure_name=="473") {structure_description="Hf_{5}Sc*-h51";return TRUE;}
    if(structure_name=="474") {structure_description="Hf_{5}Sc*-h51";return TRUE;}
    if(structure_name=="475") {structure_description="BiHf_{2}*-134";return TRUE;}
    if(structure_name=="476") {structure_description="BiHf_{2}*-134";return TRUE;}
    if(structure_name=="477") {structure_description="Hf_{5}Pb*-f63";return TRUE;}
    if(structure_name=="478") {structure_description="Hf_{5}Pb*-f63";return TRUE;}
    if(structure_name=="479") {structure_description="HfPd_{5}";return TRUE;}
    if(structure_name=="480") {structure_description="HfPd_{5}";return TRUE;}
    if(structure_name=="481") {structure_description="Re_{25}Zr_{21}";return TRUE;}
    if(structure_name=="482") {structure_description="Re_{25}Zr_{21}";return TRUE;}
    if(structure_name=="483") {structure_description="B11_{3}";return TRUE;}
    if(structure_name=="484") {structure_description="B11_{3}";return TRUE;}
    if(structure_name=="485") {structure_description="B11_{3}'";return TRUE;}
    if(structure_name=="486") {structure_description="B11_{3}'";return TRUE;}
    if(structure_name=="487") {structure_description="Z3";return TRUE;}
    if(structure_name=="488") {structure_description="Z3";return TRUE;}
    if(structure_name=="489") {structure_description="Z3'";return TRUE;}
    if(structure_name=="490") {structure_description="Z3'";return TRUE;}
    if(structure_name=="491") {structure_description="Al_{13}Co_{4} oP102";return TRUE;}
    if(structure_name=="492") {structure_description="Al_{13}Co_{4} oP102";return TRUE;}
    if(structure_name=="493") {structure_description="Al_{13}Co_{4} mP102";return TRUE;}
    if(structure_name=="494") {structure_description="Al_{13}Co_{4} mp102";return TRUE;}
    if(structure_name=="495") {structure_description="Al_{5}Co_{2} D8_{11}";return TRUE;}
    if(structure_name=="496") {structure_description="Al_{5}Co_{2} D8_{11}";return TRUE;}
    if(structure_name=="497") {structure_description="Al_{9}Co_{2}";return TRUE;}
    if(structure_name=="498") {structure_description="Al_{9}Co_{2}";return TRUE;}
    if(structure_name=="499") {structure_description="Pd_{4}Se";return TRUE;}
    if(structure_name=="500") {structure_description="Pd_{4}Se";return TRUE;}
    if(structure_name=="501") {structure_description="NaCd_{2}-Samson";return TRUE;}
    if(structure_name=="502") {structure_description="NaCd_{2}-Samson";return TRUE;}
    if(structure_name=="503") {structure_description="Er_{3}Ru_{2}";return TRUE;}
    if(structure_name=="504") {structure_description="Er_{3}Ru_{2}";return TRUE;}
    if(structure_name=="505") {structure_description="Pt_{3}Sr_{7}";return TRUE;}
    if(structure_name=="506") {structure_description="Pt_{3}Sr_{7}";return TRUE;}
    if(structure_name=="507") {structure_description="Ir_{4}Sc_{11}";return TRUE;}
    if(structure_name=="508") {structure_description="Ir_{4}Sc_{11}";return TRUE;}
    if(structure_name=="509") {structure_description="Ru_{25}Y_{44}";return TRUE;}
    if(structure_name=="510") {structure_description="Ru_{25}Y_{44}";return TRUE;}
    if(structure_name=="511") {structure_description="RuZn_{6}";return TRUE;}
    if(structure_name=="512") {structure_description="RuZn_{6}";return TRUE;}
    if(structure_name=="513") {structure_description="NbRu-\\beta''";return TRUE;}
    if(structure_name=="514") {structure_description="NbRu-\\beta''";return TRUE;}
    if(structure_name=="515") {structure_description="\\delta-CdNi";return TRUE;}
    if(structure_name=="516") {structure_description="\\delta-CdNi";return TRUE;}
    if(structure_name=="517") {structure_description="Cd_{2}Ce";return TRUE;}
    if(structure_name=="518") {structure_description="Cd_{2}Ce";return TRUE;}
    if(structure_name=="519") {structure_description="Hg_{2}U";return TRUE;}
    if(structure_name=="520") {structure_description="Hg_{2}U";return TRUE;}
    if(structure_name=="521") {structure_description="InMg_{2}";return TRUE;}
    if(structure_name=="522") {structure_description="InMg_{2}";return TRUE;}
    if(structure_name=="523") {structure_description="Cd_{3}Er";return TRUE;}
    if(structure_name=="524") {structure_description="Cd_{3}Er";return TRUE;}
    if(structure_name=="525") {structure_description="CdMg_{3}";return TRUE;}
    if(structure_name=="526") {structure_description="CdMg_{3}";return TRUE;}
    if(structure_name=="527") {structure_description="CeNi_{3}";return TRUE;}
    if(structure_name=="528") {structure_description="CeNi_{3}";return TRUE;}
    if(structure_name=="529") {structure_description="CoSc_{3}";return TRUE;}
    if(structure_name=="530") {structure_description="CoSc_{3}";return TRUE;}
    if(structure_name=="531") {structure_description="Pb_{3}Sr";return TRUE;}
    if(structure_name=="532") {structure_description="Pb_{3}Sr";return TRUE;}
    if(structure_name=="533") {structure_description="PuNi_{3}";return TRUE;}
    if(structure_name=="534") {structure_description="PuNi_{3}";return TRUE;}
    if(structure_name=="535") {structure_description="YZn_{3}";return TRUE;}
    if(structure_name=="536") {structure_description="YZn_{3}";return TRUE;}
    if(structure_name=="537") {structure_description="Co_{2}Y_{2}*-28";return TRUE;}
    if(structure_name=="538") {structure_description="Co_{2}Y_{2}*-28";return TRUE;}
    if(structure_name=="539") {structure_description="Sc_{2}Zr*-130";return TRUE;}
    if(structure_name=="540") {structure_description="Sc_{2}Zr*-130";return TRUE;}
    if(structure_name=="541") {structure_description="Mo_{3}Ti*-81";return TRUE;}
    if(structure_name=="542") {structure_description="Mo_{3}Ti*-81";return TRUE;}
    if(structure_name=="543") {structure_description="MoTi*-80-AB";return TRUE;}
    if(structure_name=="544") {structure_description="MoTi*-80-BA";return TRUE;}
    if(structure_name=="545") {structure_description="ReTi_{2}*-81";return TRUE;}
    if(structure_name=="546") {structure_description="ReTi_{2}*-81";return TRUE;}
    if(structure_name=="547") {structure_description="Hf_{2}Tl*-6";return TRUE;}
    if(structure_name=="548") {structure_description="Hf_{2}Tl*-6";return TRUE;}
    // if(structure_name=="549") {structure_description="Be_{2}Zn*-65";return TRUE;}
    // if(structure_name=="550") {structure_description="Be_{2}Zn*-65";return TRUE;}
    if(structure_name=="551") {structure_description="Re_{3}Ru*-124";return TRUE;}
    if(structure_name=="552") {structure_description="Re_{3}Ru*-124";return TRUE;}
    if(structure_name=="553") {structure_description="D0_{21}-Cu_{2.82}P";return TRUE;}
    if(structure_name=="554") {structure_description="D0_{21}-Cu_{2.82}P";return TRUE;}
    if(structure_name=="555") {structure_description="Mg_{2}Au";return TRUE;}
    if(structure_name=="556") {structure_description="Mg_{2}Au";return TRUE;}
    if(structure_name=="557") {structure_description="BiF_{3}";return TRUE;}
    if(structure_name=="558") {structure_description="BiF_{3}";return TRUE;}
    if(structure_name=="559") {structure_description="Al_{12}Mg_{17}";return TRUE;}
    if(structure_name=="560") {structure_description="Al_{12}Mg_{17}";return TRUE;}
    if(structure_name=="561") {structure_description="MgAu_{3-x}";return TRUE;}
    if(structure_name=="562") {structure_description="MgAu_{3-x}";return TRUE;}
    if(structure_name=="563") {structure_description="MgAu_{3+x}";return TRUE;}
    if(structure_name=="564") {structure_description="MgAu_{3+x}";return TRUE;}
    if(structure_name=="565") {structure_description="Mg_{2}Cu";return TRUE;}
    if(structure_name=="566") {structure_description="Mg_{2}Cu";return TRUE;}
    if(structure_name=="567") {structure_description="Mg_{2}Ga";return TRUE;}
    if(structure_name=="568") {structure_description="Mg_{2}Ga";return TRUE;}
    if(structure_name=="569") {structure_description="MgGa_{2}";return TRUE;}
    if(structure_name=="570") {structure_description="MgGa_{2}";return TRUE;}
    if(structure_name=="571") {structure_description="Mg_{5}Ga_{2}";return TRUE;}
    if(structure_name=="572") {structure_description="Mg_{5}Ga_{2}";return TRUE;}
    if(structure_name=="573") {structure_description="Mg_{2}Ga_{5}";return TRUE;}
    if(structure_name=="574") {structure_description="Mg_{2}Ga_{5}";return TRUE;}
    if(structure_name=="575") {structure_description="MgGa";return TRUE;}
    if(structure_name=="576") {structure_description="MgGa";return TRUE;}
    if(structure_name=="577") {structure_description="Mg_{3}Hg";return TRUE;}
    if(structure_name=="578") {structure_description="Mg_{3}Hg";return TRUE;}
    if(structure_name=="579") {structure_description="Mg_{3}In";return TRUE;}
    if(structure_name=="580") {structure_description="Mg_{3}In";return TRUE;}
    if(structure_name=="581") {structure_description="Mg_{44}Rh_{7}";return TRUE;}
    if(structure_name=="582") {structure_description="Mg_{44}Rh_{7}";return TRUE;}
    if(structure_name=="581b") {structure_description="Mg_{44}Rh_{7}";return TRUE;}
    if(structure_name=="582b") {structure_description="Mg_{44}Rh_{7}";return TRUE;}
    if(structure_name=="581c") {structure_description="Mg_{44}Rh_{7}";return TRUE;}
    if(structure_name=="582c") {structure_description="Mg_{44}Rh_{7}";return TRUE;}
    if(structure_name=="583") {structure_description="CaF_{2}";return TRUE;}
    if(structure_name=="584") {structure_description="CaF_{2}";return TRUE;}
    if(structure_name=="585") {structure_description="D0_{18}-Al_{3}Ir";return TRUE;}
    if(structure_name=="586") {structure_description="D0_{18}-Al_{3}Ir";return TRUE;}
    if(structure_name=="587") {structure_description="Mn_{23}Th_{6}";return TRUE;}
    if(structure_name=="588") {structure_description="Mn_{23}Th_{6}";return TRUE;}
    if(structure_name=="589") {structure_description="Mg_{2}Zn_{11}";return TRUE;}
    if(structure_name=="590") {structure_description="Mg_{2}Zn_{11}";return TRUE;}
    ;if(structure_name=="591") {structure_description="Mg_{4}Zn_{7}";return TRUE;}
    if(structure_name=="592") {structure_description="Mg_{4}Zn_{7}";return TRUE;}
    if(structure_name=="593") {structure_description="Mg_{3}Ru_{2}";return TRUE;}
    if(structure_name=="594") {structure_description="Mg_{3}Ru_{2}";return TRUE;}
    if(structure_name=="595") {structure_description="Au_{2}V";return TRUE;}
    if(structure_name=="596") {structure_description="Au_{2}V";return TRUE;}
    if(structure_name=="597") {structure_description="Sr_{9}Mg_{38}";return TRUE;}
    if(structure_name=="598") {structure_description="Sr_{9}Mg_{38}";return TRUE;}
    if(structure_name=="600.AAAAA") {structure_description="\\sigma_{AAAAA}";pureA=TRUE;return TRUE;}
    if(structure_name=="600.AAAAB") {structure_description="\\sigma_{AAAAB}";return TRUE;}
    if(structure_name=="600.AAABA") {structure_description="\\sigma_{AAABA}";return TRUE;}
    if(structure_name=="600.AAABB") {structure_description="\\sigma_{AAABB}";return TRUE;}
    if(structure_name=="600.AABAA") {structure_description="\\sigma_{AABAA}";return TRUE;}
    if(structure_name=="600.AABAB") {structure_description="\\sigma_{AABAB}";return TRUE;}
    if(structure_name=="600.AABBA") {structure_description="\\sigma_{AABBA}";return TRUE;}
    if(structure_name=="600.AABBB") {structure_description="\\sigma_{AABBB}";return TRUE;}
    if(structure_name=="600.ABAAA") {structure_description="\\sigma_{ABAAA}";return TRUE;}
    if(structure_name=="600.ABAAB") {structure_description="\\sigma_{ABAAB}";return TRUE;}
    if(structure_name=="600.ABABA") {structure_description="\\sigma_{ABABA}";return TRUE;}
    if(structure_name=="600.ABABB") {structure_description="\\sigma_{ABABB}";return TRUE;}
    if(structure_name=="600.ABBAA") {structure_description="\\sigma_{ABBAA}";return TRUE;}
    if(structure_name=="600.ABBAB") {structure_description="\\sigma_{ABBAB}";return TRUE;}
    if(structure_name=="600.ABBBA") {structure_description="\\sigma_{ABBBA}";return TRUE;}
    if(structure_name=="600.ABBBB") {structure_description="\\sigma_{ABBBB}";return TRUE;}
    if(structure_name=="600.BAAAA") {structure_description="\\sigma_{BAAAA}";return TRUE;}
    if(structure_name=="600.BAAAB") {structure_description="\\sigma_{BAAAB}";return TRUE;}
    if(structure_name=="600.BAABA") {structure_description="\\sigma_{BAABA}";return TRUE;}
    if(structure_name=="600.BAABB") {structure_description="\\sigma_{BAABB}";return TRUE;}
    if(structure_name=="600.BABAA") {structure_description="\\sigma_{BABAA}";return TRUE;}
    if(structure_name=="600.BABAB") {structure_description="\\sigma_{BABAB}";return TRUE;}
    if(structure_name=="600.BABBA") {structure_description="\\sigma_{BABBA}";return TRUE;}
    if(structure_name=="600.BABBB") {structure_description="\\sigma_{BABBB}";return TRUE;}
    if(structure_name=="600.BBAAA") {structure_description="\\sigma_{BBAAA}";return TRUE;}
    if(structure_name=="600.BBAAB") {structure_description="\\sigma_{BBAAB}";return TRUE;}
    if(structure_name=="600.BBABA") {structure_description="\\sigma_{BBABA}";return TRUE;}
    if(structure_name=="600.BBABB") {structure_description="\\sigma_{BBABB}";return TRUE;}
    if(structure_name=="600.BBBAA") {structure_description="\\sigma_{BBBAA}";return TRUE;}
    if(structure_name=="600.BBBAB") {structure_description="\\sigma_{BBBAB}";return TRUE;}
    if(structure_name=="600.BBBBA") {structure_description="\\sigma_{BBBBA}";return TRUE;}
    if(structure_name=="600.BBBBB") {structure_description="\\sigma_{BBBBB}";pureB=TRUE;return TRUE;}
    if(structure_name=="611") {structure_description="Rh_{13}Sc_{57}";return TRUE;}
    if(structure_name=="612") {structure_description="Rh_{13}Sc_{57}";return TRUE;}
    if(structure_name=="613") {structure_description="Fe_{7}W_{6}";return TRUE;}
    if(structure_name=="614") {structure_description="Fe_{7}W_{6}";return TRUE;}
    if(structure_name=="615") {structure_description="Au_{1}Sn_{2}";return TRUE;}
    if(structure_name=="616") {structure_description="Au_{1}Sn_{2}";return TRUE;}
    if(structure_name=="617") {structure_description="Bi_{2}Pd_{1}";return TRUE;}
    if(structure_name=="618") {structure_description="Bi_{2}Pd_{1}";return TRUE;}
    if(structure_name=="619") {structure_description="Bi_{2}Pd_{5}";return TRUE;}
    if(structure_name=="620") {structure_description="Bi_{2}Pd_{5}";return TRUE;}
    if(structure_name=="621") {structure_description="Bi_{3}Ni_{1}";return TRUE;}
    if(structure_name=="622") {structure_description="Bi_{3}Ni_{1}";return TRUE;}
    if(structure_name=="623") {structure_description="Bi_{3}Y_{5}";return TRUE;}
    if(structure_name=="624") {structure_description="Bi_{3}Y_{5}";return TRUE;}
    if(structure_name=="625") {structure_description="Bi_{4}Rh_{1}";return TRUE;}
    if(structure_name=="626") {structure_description="Bi_{4}Rh_{1}";return TRUE;}
    if(structure_name=="627") {structure_description="Bi_{1}Mn_{3}";return TRUE;}
    if(structure_name=="628") {structure_description="Bi_{1}Mn_{3}";return TRUE;}
    if(structure_name=="629") {structure_description="Bi_{1}Pd_{3}";return TRUE;}
    if(structure_name=="630") {structure_description="Bi_{1}Pd_{3}";return TRUE;}
    if(structure_name=="631") {structure_description="Bi_{1}Ti_{2}";return TRUE;}
    if(structure_name=="632") {structure_description="Bi_{1}Ti_{2}";return TRUE;}
    if(structure_name=="633") {structure_description="Co_{1}Sb_{2}";return TRUE;}
    if(structure_name=="634") {structure_description="Co_{1}Sb_{2}";return TRUE;}
    if(structure_name=="635") {structure_description="Cu_{7}Hg_{6}";return TRUE;}
    if(structure_name=="636") {structure_description="Cu_{7}Hg_{6}";return TRUE;}
    if(structure_name=="637") {structure_description="Hg_{2}K_{1}";return TRUE;}
    if(structure_name=="638") {structure_description="Hg_{2}K_{1}";return TRUE;}
    if(structure_name=="639") {structure_description="Ir_{2}Zn_{11}";return TRUE;}
    if(structure_name=="640") {structure_description="Ir_{2}Zn_{11}";return TRUE;}
    if(structure_name=="641") {structure_description="Ni_{3}P_{1}";return TRUE;}
    if(structure_name=="642") {structure_description="Ni_{3}P_{1}";return TRUE;}
    if(structure_name=="643.AB") {structure_description="KP_{15} aP32";return TRUE;}
    if(structure_name=="643.BA") {structure_description="KP_{15} aP32";return TRUE;}
    if(structure_name=="644.AB") {structure_description="Pu_{28}Zr tI116";return TRUE;}
    if(structure_name=="644.BA") {structure_description="Pu_{28}Zr tI116";return TRUE;}
    if(structure_name=="645.AB") {structure_description="B_{12}U cF52";return TRUE;}
    if(structure_name=="645.BA") {structure_description="B_{12}U cF52";return TRUE;}
    if(structure_name=="646.AB") {structure_description="TiZn_{16} oS68";return TRUE;}
    if(structure_name=="646.BA") {structure_description="TiZn_{16} oS68";return TRUE;}
    if(structure_name=="647.AB") {structure_description="Rh_{2}Y_{3} tI140";return TRUE;}
    if(structure_name=="647.BA") {structure_description="Rh_{2}Y_{3} tI140";return TRUE;}
    if(structure_name=="648.AB") {structure_description="Rh_{3}Pu_{5} tP32";return TRUE;}
    if(structure_name=="648.BA") {structure_description="Rh_{3}Pu_{5} tP32";return TRUE;}
    if(structure_name=="649.AB") {structure_description="SV_{3} tI32";return TRUE;}
    if(structure_name=="649.BA") {structure_description="SV_{3} tI32";return TRUE;}
    if(structure_name=="650.AB") {structure_description="IrTa oP12";return TRUE;}
    if(structure_name=="650.BA") {structure_description="IrTa oP12";return TRUE;}
    if(structure_name=="651.AB") {structure_description="Ge_{4}Sm_{5} oP36";return TRUE;}
    if(structure_name=="651.BA") {structure_description="Ge_{4}Sm_{5} oP36";return TRUE;}
    if(structure_name=="652.AB") {structure_description="Cl_{2}Pb oP12";return TRUE;}
    if(structure_name=="652.BA") {structure_description="Cl_{2}Pb oP12";return TRUE;}
    if(structure_name=="653.AB") {structure_description="NbPt_{3} mP16 ";return TRUE;}
    if(structure_name=="653.BA") {structure_description="NbPt_{3} mP16 ";return TRUE;}
    if(structure_name=="654.AB") {structure_description="Hg_{4}Pt cI10";return TRUE;}
    if(structure_name=="654.BA") {structure_description="Hg_{4}Pt cI10";return TRUE;}
    if(structure_name=="655.AB") {structure_description="Pd_{4}Pu_{3} hR14";return TRUE;}
    if(structure_name=="655.BA") {structure_description="Pd_{4}Pu_{3} hR14";return TRUE;}
    if(structure_name=="656.AB") {structure_description="Er_{3}Ni_{2} hR15";return TRUE;}
    if(structure_name=="656.BA") {structure_description="Er_{3}Ni_{2} hR15";return TRUE;}
    if(structure_name=="657.AB") {structure_description="Pd_{2}Ti oI6";return TRUE;}
    if(structure_name=="657.BA") {structure_description="Pd_{2}Ti oI6";return TRUE;}
    if(structure_name=="658.AB") {structure_description="Pd_{5}Ti_{3} tP8";return TRUE;}
    if(structure_name=="658.BA") {structure_description="Pd_{5}Ti_{3} tP8";return TRUE;}
    if(structure_name=="659.AB") {structure_description="Pd_{3}Ti_{2} oS20";return TRUE;}
    if(structure_name=="659.BA") {structure_description="Pd_{3}Ti_{2} oS20";return TRUE;}
    if(structure_name=="660.AB") {structure_description="Hg_{5}Mn_{2} tP14";return TRUE;}
    if(structure_name=="660.BA") {structure_description="Hg_{5}Mn_{2} tP14";return TRUE;}
    if(structure_name=="661.AB") {structure_description="Pu_{3}Pd_{5} oS32";return TRUE;}
    if(structure_name=="661.BA") {structure_description="Pu_{3}Pd_{5} oS32";return TRUE;}
    if(structure_name=="662.AB") {fcc=TRUE;structure_description="Ga_{2}Zr oS12";return TRUE;}
    if(structure_name=="662.BA") {fcc=TRUE;structure_description="Ga_{2}Zr oS12";return TRUE;}
    if(structure_name=="663.AB") {structure_description="FeP_{4} mP30";return TRUE;}
    if(structure_name=="663.BA") {structure_description="FeP_{4} mP30";return TRUE;}
    if(structure_name=="664.AB") {structure_description="FeP_{4} mS40";return TRUE;}
    if(structure_name=="664.BA") {structure_description="FeP_{4} mS40";return TRUE;}
    if(structure_name=="665.AB") {structure_description="In_{3}Ir tP16";return TRUE;}
    if(structure_name=="665.BA") {structure_description="In_{3}Ir tP16";return TRUE;}
    if(structure_name=="666.AB") {structure_description="Fe_{6}Ge_{5} mS44";return TRUE;}
    if(structure_name=="666.BA") {structure_description="Fe_{6}Ge_{5} mS44";return TRUE;}
    if(structure_name=="667.AB") {structure_description="Fe_{3}Ga_{4} mS42";return TRUE;}
    if(structure_name=="667.BA") {structure_description="Fe_{3}Ga_{4} mS42";return TRUE;}
    if(structure_name=="668.AB") {structure_description="Al_{6}Mn oS28";return TRUE;}
    if(structure_name=="668.BA") {structure_description="Al_{6}Mn oS28";return TRUE;}
    if(structure_name=="669.AB") {structure_description="AsFe oP8";return TRUE;}
    if(structure_name=="669.BA") {structure_description="AsFe oP8";return TRUE;}
    if(structure_name=="670.AB") {structure_description="CoSn hP6";return TRUE;}
    if(structure_name=="670.BA") {structure_description="CoSn hP6";return TRUE;}
    if(structure_name=="671.AB") {structure_description="As_{3}Co cI32";return TRUE;}
    if(structure_name=="671.BA") {structure_description="As_{3}Co cI32";return TRUE;}
    if(structure_name=="672.AB") {fcc=TRUE;structure_description="A_{15}B_{1}-FCC cI32";return TRUE;}
    if(structure_name=="672.BA") {fcc=TRUE;structure_description="A_{15}B_{1}-FCC cI32";return TRUE;}
    if(structure_name=="673.AB") {bcc=TRUE;structure_description="A_{15}B_{1}-BCC cP16";return TRUE;}
    if(structure_name=="673.BA") {bcc=TRUE;structure_description="A_{15}B_{1}-BCC cP16";return TRUE;}
    if(structure_name=="674.A") {pureA=TRUE;structure_description="C32_{pure} hP3";return TRUE;}
    if(structure_name=="674.B") {pureB=TRUE;structure_description="C32_{pure} hP3";return TRUE;}
    if(structure_name=="675.A") {pureA=TRUE;structure_description="A15_{pure} cP8";return TRUE;}
    if(structure_name=="675.B") {pureB=TRUE;structure_description="A15_{pure} cP8";return TRUE;}
    if(structure_name=="676.AB") {structure_description="Fe_{6.5}Ge_{4} hP22 (ideal A14B)";return TRUE;}
    if(structure_name=="676.BA") {structure_description="Fe_{6.5}Ge_{4} hP22 (ideal A14B)";return TRUE;}
    if(structure_name=="676p.AB") {structure_description="Fe_{6.5}Ge_{4} oS42 (real A13B)";return TRUE;}
    if(structure_name=="676p.BA") {structure_description="Fe_{6.5}Ge_{4} oS42 (real A13B)";return TRUE;}
    if(structure_name=="677.AB") {structure_description="Al_{8}Cr_{5} hR26";return TRUE;}
    if(structure_name=="677.BA") {structure_description="Al_{8}Cr_{5} hR26";return TRUE;}
    if(structure_name=="678.AB") {structure_description="Co_{5}Ge_{7} tI24";return TRUE;}
    if(structure_name=="678.BA") {structure_description="Co_{5}Ge_{7} tI24";return TRUE;}
    if(structure_name=="679.AB") {structure_description="E2_{1} cP5";return TRUE;}
    if(structure_name=="679.BA") {structure_description="E2_{1} cP5";return TRUE;}
    if(structure_name=="680.AB") {structure_description="\\beta-NW_{2} tP12";return TRUE;}
    if(structure_name=="680.BA") {structure_description="\\beta-NW_{2} tP12";return TRUE;}
    if(structure_name=="681.AB") {structure_description="\\alpha-MoS_{2} hR3";return TRUE;}
    if(structure_name=="681.BA") {structure_description="\\alpha-MoS_{2} hR3";return TRUE;}
    if(structure_name=="682.AB") {structure_description="D8_{5r} hR13";return TRUE;}
    if(structure_name=="682.BA") {structure_description="D8_{5r} hR13";return TRUE;}
    if(structure_name=="683.AB") {structure_description="Na_{13}Cl_{12}v7 cI50";return TRUE;}
    if(structure_name=="683.BA") {structure_description="Na_{13}Cl_{12}v7 cI50";return TRUE;}
    if(structure_name=="684.A") {pureA=TRUE;structure_description="\\alpha-N_{2} cP8";return TRUE;}
    if(structure_name=="684.B") {pureB=TRUE;structure_description="\\alpha-N_{2} cP8";return TRUE;}  
    if(structure_name=="685.A") {pureA=TRUE;structure_description="\\alpha-N_{2}p cP8";return TRUE;}
    if(structure_name=="685.B") {pureB=TRUE;structure_description="\\alpha-N_{2}p cP8";return TRUE;}
    if(structure_name=="686.A") {pureA=TRUE;structure_description="\\beta-N_{2} hP4";return TRUE;}
    if(structure_name=="686.B") {pureB=TRUE;structure_description="\\beta-N_{2} hP4";return TRUE;}
    if(structure_name=="687.A") {pureA=TRUE;structure_description="\\epsilon-N_{2} hR16";return TRUE;}
    if(structure_name=="687.B") {pureB=TRUE;structure_description="\\epsilon-N_{2} hR16";return TRUE;}
    if(structure_name=="688.A") {pureA=TRUE;structure_description="\\gamma-N_{2} tP4";return TRUE;}
    if(structure_name=="688.B") {pureB=TRUE;structure_description="\\gamma-N_{2} tP4";return TRUE;}
    if(structure_name=="689.AB") {structure_description="WN_{2}_194af hP6";return TRUE;}
    if(structure_name=="689.BA") {structure_description="WN_{2}_194af hP6";return TRUE;}
    if(structure_name=="690.AB") {structure_description="Ag_{2}Te mP12";return TRUE;}
    if(structure_name=="690.BA") {structure_description="Ag_{2}Te mP12";return TRUE;}
    if(structure_name=="691.AB") {structure_description="C1 cF12";return TRUE;}
    if(structure_name=="691.BA") {structure_description="C1 cF12";return TRUE;}
    if(structure_name=="692.AB") {structure_description="C19 hR3";return TRUE;}
    if(structure_name=="692.BA") {structure_description="C19 hR3";return TRUE;}
    if(structure_name=="693.AB") {structure_description="C4 tP6";return TRUE;}
    if(structure_name=="693.BA") {structure_description="C4 tP6";return TRUE;}
    if(structure_name=="694.AB") {structure_description="C7 hP6";return TRUE;}
    if(structure_name=="694.BA") {structure_description="C7 hP6";return TRUE;}
    if(structure_name=="695.AB") {structure_description="WN_{2}_187fg hP3";return TRUE;}
    if(structure_name=="695.BA") {structure_description="WN_{2}_187fg hP3";return TRUE;}
    if(structure_name=="696.AB") {structure_description="WN_{2}_194de hP6";return TRUE;}
    if(structure_name=="696.BA") {structure_description="WN_{2}_194de hP6";return TRUE;}
    if(structure_name=="697.AB") {structure_description="Hagg_N_{4}W6 tP10";return TRUE;}
    if(structure_name=="697.BA") {structure_description="Hagg_N_{4}W6 tP10";return TRUE;}
    if(structure_name=="698.AB") {structure_description="MoO_{3} oP16";return TRUE;}
    if(structure_name=="698.BA") {structure_description="MoO_{3} oP16";return TRUE;}
    if(structure_name=="699.AB") {structure_description="O_{3}W hP12";return TRUE;}
    if(structure_name=="699.BA") {structure_description="O_{3}W hP12";return TRUE;}
    if(structure_name=="700.AB") {structure_description="P_{3}Tc oP16";return TRUE;}
    if(structure_name=="700.BA") {structure_description="P_{3}Tc oP16";return TRUE;}
    if(structure_name=="701.AB") {structure_description="WN_{4}_164ad hP5";return TRUE;}
    if(structure_name=="701.BA") {structure_description="WN_{4}_164ad hP5";return TRUE;}
    if(structure_name=="702.AB") {structure_description="B_{4}W hP20";return TRUE;}
    if(structure_name=="702.BA") {structure_description="B_{4}W hP20";return TRUE;}
    if(structure_name=="703.AB") {structure_description="P_{4}Re oP40";return TRUE;}
    if(structure_name=="703.BA") {structure_description="P_{4}Re oP40";return TRUE;}
    if(structure_name=="704.AB") {structure_description="B8_{1}v2 hP14";return TRUE;}
    if(structure_name=="704.BA") {structure_description="B8_{1}v2 hP14";return TRUE;}
    if(structure_name=="705.AB") {structure_description="B8_{1}v2p oS28";return TRUE;}
    if(structure_name=="705.BA") {structure_description="B8_{1}v2p oS28";return TRUE;}
    if(structure_name=="706.AB") {structure_description="Na_{4}Cl_{3}v1 oS14";return TRUE;}
    if(structure_name=="706.BA") {structure_description="Na_{4}Cl_{3}v1 oS14";return TRUE;}
    if(structure_name=="707.AB") {structure_description="C_{4}W_{3} cP7";return TRUE;}
    if(structure_name=="707.BA") {structure_description="C_{4}W_{3} cP7";return TRUE;}
    if(structure_name=="708.AB") {structure_description="S_{3}U_{4} cP7";return TRUE;}
    if(structure_name=="708.BA") {structure_description="S_{3}U_{4} cP7";return TRUE;}
    if(structure_name=="709.AB") {structure_description="O_{5}Nb_{2} mS14";return TRUE;}
    if(structure_name=="709.BA") {structure_description="O_{5}Nb_{2} mS14";return TRUE;}
    if(structure_name=="710.AB") {structure_description="B8_{1}v3 hP13";return TRUE;}
    if(structure_name=="710.BA") {structure_description="B8_{1}v3 hP13";return TRUE;}
    if(structure_name=="711.AB") {structure_description="Na_{7}Cl_{6}_v3 cF52";return TRUE;}
    if(structure_name=="711.BA") {structure_description="Na_{7}Cl_{6}_v3 cF52";return TRUE;}
    if(structure_name=="712.AB") {structure_description="Na_{7}Cl_{8}v1 cF60";return TRUE;}
    if(structure_name=="712.BA") {structure_description="Na_{7}Cl_{8}v1 cF60";return TRUE;}
    if(structure_name=="713.AB") {structure_description="B8_{1}v1 hP15";return TRUE;}
    if(structure_name=="713.BA") {structure_description="B8_{1}v1 hP15";return TRUE;}
    if(structure_name=="714.AB") {structure_description="B8_{1}v1p hP15";return TRUE;}
    if(structure_name=="714.BA") {structure_description="B8_{1}v1p hP15";return TRUE;}
    if(structure_name=="715.AB") {structure_description="Na_{1}Cl_{1}v2 tP2";return TRUE;}
    if(structure_name=="715.BA") {structure_description="Na_{1}Cl_{1}v2 tP2";return TRUE;}
    if(structure_name=="716.AB") {structure_description="Na_{3}Cl_{3}v2 tP6";return TRUE;}
    if(structure_name=="716.BA") {structure_description="Na_{3}Cl_{3}v2 tP6";return TRUE;}
    if(structure_name=="717.AB") {structure_description="B8_{1}v2pp hP14";return TRUE;}
    if(structure_name=="717.BA") {structure_description="B8_{1}v2pp hP14";return TRUE;}
    if(structure_name=="718.AB") {structure_description="B3v2 hR6";return TRUE;}
    if(structure_name=="718.BA") {structure_description="B3v2 hR6";return TRUE;}
    if(structure_name=="719.AB") {structure_description="Na_{7}Cl_{7}v2 cF56";return TRUE;}
    if(structure_name=="719.BA") {structure_description="Na_{7}Cl_{7}v2 cF56";return TRUE;}
    if(structure_name=="720.AB") {structure_description="B1 cF8";return TRUE;}
    if(structure_name=="720.BA") {structure_description="B1 cF8";return TRUE;}
    if(structure_name=="721.AB") {structure_description="NTa hP6";return TRUE;}
    if(structure_name=="721.BA") {structure_description="NTa hP6";return TRUE;}
    if(structure_name=="722.AB") {structure_description="NbO cP6";return TRUE;}
    if(structure_name=="722.BA") {structure_description="NbO cP6";return TRUE;}
    if(structure_name=="723.AB") {structure_description="NbS hP16";return TRUE;}
    if(structure_name=="723.BA") {structure_description="NbS hP16";return TRUE;}
    if(structure_name=="724.AB") {structure_description="Hagg_NW_{2} tP3";return TRUE;}
    if(structure_name=="724.BA") {structure_description="Hagg_NW_{2} tP3";return TRUE;}
    if(structure_name=="725.AB") {structure_description="Mo_{2}N tI12";return TRUE;}
    if(structure_name=="725.BA") {structure_description="Mo_{2}N tI12";return TRUE;}
    if(structure_name=="726.AB") {structure_description="NV_{2} hP9";return TRUE;}
    if(structure_name=="726.BA") {structure_description="NV_{2} hP9";return TRUE;}
    if(structure_name=="727.AB") {structure_description="NW_{2}_Pearson hP9";return TRUE;}
    if(structure_name=="727.BA") {structure_description="NW_{2}_Pearson hP9";return TRUE;}
    if(structure_name=="728.A") {pureA=TRUE;structure_description="beta-W cP8";return TRUE;}
    if(structure_name=="728.B") {pureB=TRUE;structure_description="beta-W cP8";return TRUE;}
    if(structure_name=="729.AB") {structure_description="Ge_{8}Mn_{11} oP76";return TRUE;}
    if(structure_name=="729.BA") {structure_description="Ge_{8}Mn_{11} oP76";return TRUE;}
    if(structure_name=="730.AB") {structure_description="NaTe_{3} hP48";return TRUE;}
    if(structure_name=="730.BA") {structure_description="NaTe_{3} hP48";return TRUE;}
    if(structure_name=="731.AB") {structure_description="NaTe oP48";return TRUE;}
    if(structure_name=="731.BA") {structure_description="NaTe oP48";return TRUE;}
    if(structure_name=="732.AB") {structure_description="Sb_{3}Sr_{2} mP40";return TRUE;}
    if(structure_name=="732.BA") {structure_description="Sb_{3}Sr_{2} mP40";return TRUE;}
    if(structure_name=="733.AB") {structure_description="SrAs_{3} mS16";return TRUE;}
    if(structure_name=="733.BA") {structure_description="SrAs_{3} mS16";return TRUE;}
    if(structure_name=="734.AB") {structure_description="La_{2}Sb tI12";return TRUE;}
    if(structure_name=="734.BA") {structure_description="La_{2}Sb tI12";return TRUE;}
    if(structure_name=="735.AB") {structure_description="Ga_{5}Tm_{3} oP32";return TRUE;}
    if(structure_name=="735.BA") {structure_description="Ga_{5}Tm_{3} oP32";return TRUE;}
    if(structure_name=="ICSD_10509.AB") {structure_description="ICSD_{10509} aP15";return TRUE;}
    if(structure_name=="ICSD_10509.BA") {structure_description="ICSD_{10509} aP15";return TRUE;}
    if(structure_name=="ICSD_14374.AB") {structure_description="ICSD_{14374} cP138";return TRUE;}
    if(structure_name=="ICSD_14374.BA") {structure_description="ICSD_{14374} cP138";return TRUE;}
    if(structure_name=="ICSD_43249.AB") {structure_description="ICSD_{43249} oP8";return TRUE;}
    if(structure_name=="ICSD_43249.BA") {structure_description="ICSD_{43249} oP8";return TRUE;}
    if(structure_name=="ICSD_56273.AB") {structure_description="ICSD_{56273} cP2";return TRUE;}
    if(structure_name=="ICSD_56273.BA") {structure_description="ICSD_{56273} cP2";return TRUE;}
    if(structure_name=="ICSD_104638.AB") {structure_description="ICSD_{104638} oP102";return TRUE;}
    if(structure_name=="ICSD_104638.BA") {structure_description="ICSD_{104638} oP102";return TRUE;}
    if(structure_name=="ICSD_155840.AB") {structure_description="ICSD_{155840} cP16";return TRUE;}
    if(structure_name=="ICSD_155840.BA") {structure_description="ICSD_{155840} cP16";return TRUE;}
    if(structure_name=="ICSD_155842.AB") {structure_description="ICSD_{155842} tP16";return TRUE;}
    if(structure_name=="ICSD_155842.BA") {structure_description="ICSD_{155842} tP16";return TRUE;}
    if(structure_name=="ICSD_155843.AB") {structure_description="ICSD_{155843} tP4";return TRUE;}
    if(structure_name=="ICSD_155843.BA") {structure_description="ICSD_{155843} tP4";return TRUE;}
    if(structure_name=="ICSD_155844.AB") {structure_description="ICSD_{155844} tP16";return TRUE;}
    if(structure_name=="ICSD_155844.BA") {structure_description="ICSD_{155844} tP16";return TRUE;}
    if(structure_name=="ICSD_155846.AB") {structure_description="ICSD_{155846} cP16";return TRUE;}
    if(structure_name=="ICSD_155846.BA") {structure_description="ICSD_{155846} cP16";return TRUE;}
    if(structure_name=="ICSD_169787.AB") {structure_description="ICSD_{169787} hP9";return TRUE;}
    if(structure_name=="ICSD_169787.BA") {structure_description="ICSD_{169787} hP9";return TRUE;}
    if(structure_name=="ICSD_249632.AB") {structure_description="ICSD_{249632} tP14";return TRUE;}
    if(structure_name=="ICSD_249632.BA") {structure_description="ICSD_{249632} tP14";return TRUE;}
    if(structure_name=="ICSD_607482.AB") {structure_description="ICSD_{607482} cF16";return TRUE;}
    if(structure_name=="ICSD_607482.BA") {structure_description="ICSD_{607482} cF16";return TRUE;}
    
    if(structure_name=="ICSD_105856.AB") {structure_description="Pt_{7}Zn_{12}";return TRUE;}
    if(structure_name=="ICSD_105856.BA") {structure_description="Pt_{7}Zn_{12}";return TRUE;}
    if(structure_name=="ICSD_107575.AB") {structure_description="CoZn_{13}";return TRUE;}
    if(structure_name=="ICSD_107575.BA") {structure_description="CoZn_{13}";return TRUE;}
    
    if(aurostd::substring2bool(structure_name,"ICSD_") && aurostd::substring2bool(structure_name,".")) {
      vector<string> tokens;
      aurostd::string2tokens(structure_name,tokens,".");
      structure_description=tokens.at(0);
      aurostd::StringSubst(structure_description,"ICSD_","ICSD_{");
      structure_description+="}";
      return TRUE;}

    // default
    structure_description=structure_name;
    if(structure_description.at(0)=='f') fcc=TRUE;
    if(structure_description.at(0)=='b') bcc=TRUE;
    if(structure_description.at(0)=='h') hcp=TRUE;

    return FALSE;
  }
}
#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
