// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
//

#ifndef _AFLOW_APENNSY_PHASES_CPP_
#define _AFLOW_APENNSY_PHASES_CPP_
#include "aflow.h"

// ***************************************************************************
//  (double) 2/3;62;"/62/")_STR62)
//  (double) 1/3;63;"/63/")_STR62)
//  (double) 2/3;64;"/64/")_STR64)
//  (double) 1/3;65;"/65/")_STR64)

#define _Ag_GS_ string("fcc")
#define _Al_GS_ string("fcc")
#define _Au_GS_ string("fcc")
#define _B_GS_  string("tet")
#define _Be_GS_ string("hex")
#define _Bi_GS_ string("rhl")
#define _Cd_GS_ string("hcp")
#define _Co_GS_ string("hcp")
#define _Cr_GS_ string("bcc")
#define _Cu_GS_ string("fcc")
#define _Fe_GS_ string("bcc")
#define _Ga_GS_ string("A11")
#define _Ge_GS_ string("dia")

#define _Hf_GS_ string("hcp")
#define _Hg_GS_ string("rhl")
#define _In_GS_ string("A6")
#define _Ir_GS_ string("fcc")
#define _K_GS_  string("hcp (bcc)")
#define _La_GS_ string("hcp/fcc")
#define _Li_GS_ string("fcc/bcc       ")
#define _Mg_GS_ string("hcp")
#define _Mn_GS_ string("A12 ")
#define _Mo_GS_ string("bcc")
#define _Na_GS_ string("fcc/hcp (bcc)")
#define _Nb_GS_ string("bcc")
#define _Ni_GS_ string("fcc")
#define _Os_GS_ string("hcp")
#define _Pb_GS_ string("fcc")
#define _Pd_GS_ string("fcc")
#define _Pt_GS_ string("fcc")
#define _Re_GS_ string("hcp")
#define _Rh_GS_ string("fcc")
#define _Ru_GS_ string("hcp")
#define _Sc_GS_ string("hcp")
#define _Si_GS_ string("dia")
#define _Sm_GS_ string("rhl")
#define _Sn_GS_ string("dia (tet)")
#define _Sr_GS_ string("fcc")
#define _Ta_GS_ string("bcc")
#define _Tc_GS_ string("hcp")
#define _Ti_GS_ string("hcp")
#define _Tl_GS_ string("hcp")
#define _V_GS_ string("bcc")
#define _W_GS_ string("bcc")
#define _Y_GS_ string("hcp")
#define _Zn_GS_ string("hcp")
#define _Zr_GS_ string("hcp")

using aurostd::isequal;

// ***************************************************************************

bool APENNSY_Parameters::FixGndStatesNamesConcentrations(bool verbose) {
  uint j,k;
  bool vPRE=verbose,vPOS=verbose;
  vPRE=FALSE;vPOS=FALSE;
  cerr << "LIB_MODE=" << LIB_MODE << endl;
  vector<string> ventries;
  vector<uint> vjindex;
  if(verbose) cerr << "***************************************************************" << endl;
  if(verbose) cerr << "FixGndStatesNamesConcentrations: start" << endl;
 
  for(k=0;k<alloys.size();k++) {
    vjindex.clear();
    for(j=1;j<ZConcentrations.at(k).size();j++) {
      if(ZConcentrations.at(k).at(j)>-0.001 && ZGNDlibrary.at(k).at(j)>=0)  { // check that it is a gndstate
	// cerr << "ZLibrary.at(k).size()=" << ZLibrary.at(k).size() << endl;
	// cerr << "ZConcentrations.at(k).size()=" << ZConcentrations.at(k).size() << endl;
	// cerr << "ZConcentrations.at(k).at(j)=" << ZConcentrations.at(k).at(j) << endl;
		
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  AgAu
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Ag","Au")) {
	    PhaseDiagramNameLeft.at(k)="Ag";PhaseDiagramNameRight.at(k)="Au";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ag_GS_; // Ag
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D1_a/tie";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="(see table)  ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C37/MoPt_2  ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_0";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="  C37/MoPt_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="            L1_2/D0_{22}/D0_{23}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Au_GS_; // Au
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  AgAl
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ag","Al")) {
	    PhaseDiagramNameLeft.at(k)="Ag";PhaseDiagramNameRight.at(k)="Al";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ag_GS_; // Ag
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Al_GS_; // Al
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  AgCd
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Ag","Cd")) {
	    PhaseDiagramNameLeft.at(k)="Ag";PhaseDiagramNameRight.at(k)="Cd";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ag_GS_; // Ag
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{19}/D0_{22}/D0_{24}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C37";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B2/B19/B27";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="ZrSi_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{19}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Cd_GS_; // Cd
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  AgHf
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ag","Hf")) {
	    PhaseDiagramNameLeft.at(k)="Ag";PhaseDiagramNameRight.at(k)="Hf";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ag_GS_; // Ag
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B11";
            if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b/CuZr_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  AgMg
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Ag","Mg")) {
	    PhaseDiagramNameLeft.at(k)="Ag";PhaseDiagramNameRight.at(k)="Mg";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ag_GS_; // Ag
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{23}/D0_{24}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B2";
	    //if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B8_2/Ni_2Si/tie";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="      D0_{19}/D0_{a}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Mg_GS_; // Mg
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  AgNa
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Ag","Na")) {
	    PhaseDiagramNameLeft.at(k)="Ag";PhaseDiagramNameRight.at(k)="Na";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ag_GS_; // Ag
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Na_GS_; // Na
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  AgPd
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ag","Pd")) {
	    PhaseDiagramNameLeft.at(k)="Ag";PhaseDiagramNameRight.at(k)="Pd";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ag_GS_);
	    ventries.push_back("CuPt_{7}");
	    ventries.push_back("HfPd_{5}^{*}  ");
	    ventries.push_back("D0_{23}");
	    ventries.push_back("C37");
	    ventries.push_back("                                 Ag_{2}Pd_{3}^{+}");
	    ventries.push_back("             L1_{1}");
	    ventries.push_back("          CdPt_3^{*}"); 
	    ventries.push_back(_Pd_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Ag","Pd")) {
	    PhaseDiagramNameLeft.at(k)="Pd";PhaseDiagramNameRight.at(k)="Ag";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pd_GS_; // Pd
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D1_a/tie";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_1";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C37/C49";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_2/D0_{22}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ag_GS_; // Ag
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** AgPt
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ag","Pt")) {
	    PhaseDiagramNameLeft.at(k)="Ag";PhaseDiagramNameRight.at(k)="Pt";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ag_GS_);
	    ventries.push_back("CuPt_{7}");
	    ventries.push_back("Ag_{3}Pt_{2}^{+}");
	    ventries.push_back("          L1_{1}");
	    ventries.push_back(_Pt_GS_);
	  }
	// ***************************************************************  AgTi
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Ag","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Ti";PhaseDiagramNameRight.at(k)="Ag";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B11";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ag_GS_; // Ag
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	if(LIB_MODE==LIB_MODEX || LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ag","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Ag";PhaseDiagramNameRight.at(k)="Ti";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ag_GS_; // Ag
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B11";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  AgY
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ag","Y")) {
	    PhaseDiagramNameLeft.at(k)="Ag";PhaseDiagramNameRight.at(k)="Y";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ag_GS_; // Ag
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Y_GS_; // Y
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Ag","Y")) {
	    PhaseDiagramNameLeft.at(k)="Y";PhaseDiagramNameRight.at(k)="Ag";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Y_GS_; // Y
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C37/tie";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b";
	    //if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="A15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="  D0_a";
	    //if(isequal(ZConcentrations.at(k).at(j),(double) 5/6,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C15_b/tie";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ag_GS_; // Ag
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  AgZr
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Ag","Zr")) {
	    PhaseDiagramNameLeft.at(k)="Zr";PhaseDiagramNameRight.at(k)="Ag";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zr_GS_; // Zr
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_STR15_Z3+"/tie";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B11";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C6/C32";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ag_GS_; // Ag
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  AlHf
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Al","Hf")) {
	    PhaseDiagramNameLeft.at(k)="Al";PhaseDiagramNameRight.at(k)="Hf";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Al_GS_; // Al
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{23}";
 	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C14";
       	    if(isequal(ZConcentrations.at(k).at(j),(double) 4/7,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="Al_3Zr_4";
 	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  AlSc
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Al","Sc")) {
	    PhaseDiagramNameLeft.at(k)="Sc";PhaseDiagramNameRight.at(k)="Al";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Sc_GS_; // Sc
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{19}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B8_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Al_GS_; // Al
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  AlTi
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Al","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Al";PhaseDiagramNameRight.at(k)="Ti";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Al_GS_; // Al
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C32"; //
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  AlZr
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Al","Zr")) {
	    PhaseDiagramNameLeft.at(k)="Al";PhaseDiagramNameRight.at(k)="Zr";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Al_GS_; // Al
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C32"; //
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zr_GS_; // Zr
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  AuCd
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Au","Cd")) {
	    PhaseDiagramNameLeft.at(k)="Au";PhaseDiagramNameRight.at(k)="Cd";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Au_GS_; // Au
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{24}/D0_{19}/Al_3Pu";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B19";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L6_0";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Cd_GS_; // Cd
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  AuHf
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Au","Hf")) {
	    PhaseDiagramNameLeft.at(k)="Au";PhaseDiagramNameRight.at(k)="Hf";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Au_GS_; // Au // Au
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="Au_4Zr";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="   D0_a";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="   C11_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/7,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="          Cu_4Ti_3";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="      B11";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=" C11_b/CuZr_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  AuNb
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Au","Nb")) {
	    PhaseDiagramNameLeft.at(k)="Au";PhaseDiagramNameRight.at(k)="Nb";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Au_GS_; // Au
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C32 (paw-gga)";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_STR64;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Nb_GS_; // Nb
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  AuPd
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Au","Pd")) {
	    PhaseDiagramNameLeft.at(k)="Au";PhaseDiagramNameRight.at(k)="Pd";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Au_GS_);
	    ventries.push_back("HfPd_{5}^{*}");
	    ventries.push_back("D0_{23}");
	    ventries.push_back("C49");
	    ventries.push_back("NbP");
	    ventries.push_back("         L1_{2}");
	    ventries.push_back(_Pd_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Au","Pd")) {
	    PhaseDiagramNameLeft.at(k)="Au";PhaseDiagramNameRight.at(k)="Pd";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Au_GS_; // Au
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D1_a/tie      ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{23}/D0_{22}/L1_2     ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C49/C37";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_STR23_CH;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="                    D0_{23}/D0_{22}/L1_2/tie";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pd_GS_; // Pd
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  AuRu
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Au","Ru")) {
	    PhaseDiagramNameLeft.at(k)="Au";PhaseDiagramNameRight.at(k)="Ru";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Au_GS_; // Au // Au
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ru_GS_; // Ru
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  AuSc
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Au","Sc")) {
	    PhaseDiagramNameLeft.at(k)="Au";PhaseDiagramNameRight.at(k)="Sc";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Au_GS_; // Au
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="  D1_a   ";
	    // if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="   D0_a";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="  C11_b/MoPt_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B2/B19";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C37  ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Sc_GS_; // Sc
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  AuTi
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Au","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Au";PhaseDiagramNameRight.at(k)="Ti";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Au_GS_; // Au
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="  D1_a   ";
	    //if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="  D0_a";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b/MoPt_2            ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/7,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="        Cu_4Ti_3";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="        B11";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="A15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	if(LIB_MODE==LIB_MODEX || LIB_MODE==LIB_MODEN)
	  if(strsub(alloys.at(k),"Au","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Au";PhaseDiagramNameRight.at(k)="Ti";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Au_GS_; // Au
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="  D1_a   ";
	    //if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="  D0_a";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b/MoPt_2            ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/7,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="        Cu_4Ti_3";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="        B11";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="A15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  AuY
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Au","Y")) {
	    PhaseDiagramNameLeft.at(k)="Au";PhaseDiagramNameRight.at(k)="Y";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Au_GS_; // Au
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="   D0_a";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="  B33";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C37";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Y_GS_; // Y
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  AuZr
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Au","Zr")) {
	    PhaseDiagramNameLeft.at(k)="Au";PhaseDiagramNameRight.at(k)="Zr";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Au_GS_; // Au
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="   D0_a";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/7,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="          Cu_4Ti_3";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="                           B11/"+_STR14_Z2;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="         CuZr_2/C11_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="A15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zr_GS_; // Zr
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  BNa
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"B","Na")) {
	    PhaseDiagramNameLeft.at(k)="B";PhaseDiagramNameRight.at(k)="Na";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_B_GS_; // B
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Na_GS_; // Na
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  BTi
	if(LIB_MODE==LIB_MODEX || LIB_MODE==LIB_MODEN)
	  if(strsub(alloys.at(k),"B","Ti")) {
	    PhaseDiagramNameLeft.at(k)="B";PhaseDiagramNameRight.at(k)="Ti";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_B_GS_; // B
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C32";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="        CrB/FeB-b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  BSm
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"B","Sm")) {
	    PhaseDiagramNameLeft.at(k)="B";PhaseDiagramNameRight.at(k)="Sm";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_B_GS_; // B
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/(6+1),CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="SmB_6 (cP7)         ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 4/(16+4),CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="    SmB_4 (tP20)";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 8/(20+8),CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="                           Sm_{2}B_{5} (mP28)";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Sm_GS_; // Sm
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  BeCr
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Be","Cr")) {
	    PhaseDiagramNameLeft.at(k)="Be";PhaseDiagramNameRight.at(k)="Cr";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Be_GS_; // Be
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Cr_GS_; // Cr
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  BeHf
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Be","Hf")) {
	    PhaseDiagramNameLeft.at(k)="Be";PhaseDiagramNameRight.at(k)="Hf";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Be_GS_;// Be
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/14,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="          D2_3";
            if(isequal(ZConcentrations.at(k).at(j),(double) 2/19,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="          Th_2Zn_{17}";
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/6,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="                      D2_d";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  BeNa
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Be","Na")) {
	    PhaseDiagramNameLeft.at(k)="Be";PhaseDiagramNameRight.at(k)="Na";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Be_GS_; // Be
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Na_GS_; // Na
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  BeTi
	if(LIB_MODE==LIB_MODEX || LIB_MODE==LIB_MODEN)
	  if(strsub(alloys.at(k),"Be","Ti")) {
	    cerr << "HERE" << endl;
	    PhaseDiagramNameLeft.at(k)="Be ";PhaseDiagramNameRight.at(k)="Ti";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Be_GS_; // B
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/6,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D2_d CaCu_5";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="        B2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  BiHf
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Bi","Hf")) {
	    PhaseDiagramNameLeft.at(k)="Bi";PhaseDiagramNameRight.at(k)="Hf";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Bi_GS_;// Bi
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C16";
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B11";
            if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="BiHf_2^*";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  CdHf
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Cd","Hf")) {
	    PhaseDiagramNameLeft.at(k)="Cd";PhaseDiagramNameRight.at(k)="Hf";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Cd_GS_;// Cd
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="CdTi/B11   ";
            if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="   "+_STR5_BETA1;	  
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  CdPd
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Cd","Pd")) {
	    PhaseDiagramNameLeft.at(k)="Cd";PhaseDiagramNameRight.at(k)="Pd";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Cd_GS_);
	    ventries.push_back("Ir_{2}Zn_{11}");
	    ventries.push_back("Hg_{2}Pt     ");
	    ventries.push_back("L1_{0}");
	    ventries.push_back("C37");
	    ventries.push_back("     D0_{22}");
	    ventries.push_back("     D1_{a}");
	    ventries.push_back("            HfPd_{5}^{*}");
	    ventries.push_back("     CuPt_{7}");
	    ventries.push_back(_Pd_GS_);
	  }
	// *************************************************************** CdPt
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Cd","Pt")) {
	    PhaseDiagramNameLeft.at(k)="Cd";PhaseDiagramNameRight.at(k)="Pt";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Cd_GS_);
	    ventries.push_back("Hg_{4}Pt");
	    ventries.push_back("D0_{11}  ");
	    ventries.push_back("Hg_{2}Pt     ");
	    ventries.push_back("         L1_{0}");
	    ventries.push_back("             CdPt_3^{*}");
	    ventries.push_back("   CuPt_{7}");
	    ventries.push_back(_Pt_GS_);
	  }
	// ***************************************************************  CdPd
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Cd","Pd")) {
	    PhaseDiagramNameLeft.at(k)="Pd";PhaseDiagramNameRight.at(k)="Cd";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Cd_GS_; // Cd
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{22}/NbPd_3";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C37";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_0";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{19}/NbPd_3/Al_3Pu/D0_{24}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pd_GS_; // Pd
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  CdPt
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Cd","Pt")) {
	    PhaseDiagramNameLeft.at(k)="Cd";PhaseDiagramNameRight.at(k)="Pt";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="_Cd_GS_";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{11}/D0_{a}/D0_{22}     ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C37/C16/tie     ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_0";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="     CdPt_3^{proto}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pt_GS_; // Pt
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** CdRh
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Cd","Rh")) {
	    PhaseDiagramNameLeft.at(k)="Cd";PhaseDiagramNameRight.at(k)="Rh";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Cd_GS_);
	    ventries.push_back("Hg_{4}Pt");
	    ventries.push_back("Hg_{2}Pt     ");
	    ventries.push_back(_Rh_GS_);
	  }
	if(0 && LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Cd","Rh")) {
	    PhaseDiagramNameLeft.at(k)="Cd";PhaseDiagramNameRight.at(k)="Rh";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="_Cd_GS_"; // Cd
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{24}/Al_3Pu";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C37";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Rh_GS_; // Rh
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Cd","Rh")) {
	    PhaseDiagramNameLeft.at(k)="Rh";PhaseDiagramNameRight.at(k)="Cd";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Rh_GS_; // Rh
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C37           ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{24}/Al_3Pu";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="_Cd_GS_"; // Cd
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  CdTi
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Cd","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Ti";PhaseDiagramNameRight.at(k)="Cd";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="_Cd_GS_";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/7,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="               Cu_4Ti_3/tie";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="       B11";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	if(LIB_MODE==LIB_MODEX || LIB_MODE==LIB_MODEN)
	  if(strsub(alloys.at(k),"Cd","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Cd";PhaseDiagramNameRight.at(k)="Ti";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="_Cd_GS_";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b/CuZr_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="       B11";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  CdY
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Cd","Y")) {
	    PhaseDiagramNameLeft.at(k)="Y";PhaseDiagramNameRight.at(k)="Cd";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Y_GS_; // Y
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C37/C49  ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C6";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="Cd_3Er";//D0_{19}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 5/6,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="";//YCd_5^{proto}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="_Cd_GS_";
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  CdZr
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Cd","Zr")) {
	    PhaseDiagramNameLeft.at(k)="Zr";PhaseDiagramNameRight.at(k)="Cd";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="_Cd_GS_";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="A15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B11 (paw-gga) ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zr_GS_; // Zr
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  CoHf
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Co","Hf")) {
	    PhaseDiagramNameLeft.at(k)="Co";PhaseDiagramNameRight.at(k)="Hf";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="HCP";// Co
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C14";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B33";
            if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C37";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  CoPd
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Co","Pd")) {
	    PhaseDiagramNameLeft.at(k)="Co";PhaseDiagramNameRight.at(k)="Pd";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Co_GS_);
	    //	    ventries.push_back("            "+_STR13_Z1);
	    //	            ventries.push_back("            tet-"+_STR13_Z1+" c/a=2.8"); 
	            ventries.push_back("            tet-"+_STR13_Z1+" c/a=2.8"); 
	    ventries.push_back(_Pd_GS_);
	  }
	// *************************************************************** CoPt
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Co","Pt")) {
	    PhaseDiagramNameLeft.at(k)="Co";PhaseDiagramNameRight.at(k)="Pt";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Co_GS_);
	    ventries.push_back("D0_{19}");
	    ventries.push_back("      CuZr_{2}");
	    ventries.push_back("          HfPd_{5}^{*}");
	    ventries.push_back(_Pt_GS_);
	  }
	// ***************************************************************  CoTi
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Co","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Co";PhaseDiagramNameRight.at(k)="Ti";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="HCP";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C14/C36";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="     C37/CuZr_2/NiTi_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  CrHf
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Cr","Hf")) {
	    PhaseDiagramNameLeft.at(k)="Cr";PhaseDiagramNameRight.at(k)="Hf";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Cr_GS_; // Cr
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="   C15"; //
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** CrIr
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Cr","Ir")) {
	    PhaseDiagramNameLeft.at(k)="Cr";PhaseDiagramNameRight.at(k)="Ir";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Cr_GS_);
	    ventries.push_back("B19");
	    ventries.push_back("C37");
	    ventries.push_back("     D0_{19}");
	    ventries.push_back(_Ir_GS_);
	  }
	// ***************************************************************  CrOs  
        if(LIB_MODE==LIB_MODEX)  
          if(strsub(alloys.at(k),"Cr","Os")) {  
            PhaseDiagramNameLeft.at(k)="Cr";PhaseDiagramNameRight.at(k)="Os";
            if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
            if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Cr_GS_; // Cr  
            if(isequal(ZConcentrations.at(k).at(j),(double) 0.75,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{19}";
            if(isequal(ZConcentrations.at(k).at(j),(double) 5.0/6.0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="   Hf_5Sc^{*}";
            if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Os_GS_; // Os  
            if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
          }  
	// ***************************************************************  CrPd
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Cr","Pd")) {
	    PhaseDiagramNameLeft.at(k)="Cr";PhaseDiagramNameRight.at(k)="Pd";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Cr_GS_);
	    ventries.push_back("    L1_{2}");
	    ventries.push_back("             HfPd_{5}^{*}");
	    ventries.push_back(_Pd_GS_);
	  }
	// *************************************************************** CrPt
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Cr","Pt")) {
	    PhaseDiagramNameLeft.at(k)="Cr";PhaseDiagramNameRight.at(k)="Pt";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Cr_GS_);
	    ventries.push_back("B19    ");
	    ventries.push_back("     L1_{2}");
	    ventries.push_back("    Pt_{8}Ti");
	    ventries.push_back(_Pt_GS_);
	  }
	// *************************************************************** CrRh
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Cr","Rh")) {
	    PhaseDiagramNameLeft.at(k)="Cr";PhaseDiagramNameRight.at(k)="Rh";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Cr_GS_);
	    ventries.push_back("C37              ");
	    ventries.push_back("      L1_{2}");
	    ventries.push_back("    CuPt_{7}");
	    ventries.push_back(_Rh_GS_);
	  }
	// ***************************************************************  CrTi
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Cr","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Cr";PhaseDiagramNameRight.at(k)="Ti";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Cr_GS_; // Cr
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  CuHf
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Cu","Hf")) {
	    PhaseDiagramNameLeft.at(k)="Cu";PhaseDiagramNameRight.at(k)="Hf";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Cu_GS_; // Cu
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/6,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C15_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/11,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="Cu_8Hf_3";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 7/17,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="         Ni_{10}Hf_7";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b/CuZr_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  CuPd
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Cu","Pd")) {
	    PhaseDiagramNameLeft.at(k)="Cu";PhaseDiagramNameRight.at(k)="Pd";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Cu_GS_);
	    ventries.push_back("HfPd_{5}^{*}");
	    ventries.push_back("D0_{23}");
	    ventries.push_back("Ga_{2}Hf    ");
	    ventries.push_back("B2");
	    ventries.push_back("   L1_{2}");
	    ventries.push_back("     CuPt_{7}");
	    ventries.push_back(_Pd_GS_);
	  }
	// *************************************************************** CuPt
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Cu","Pt")) {
	    PhaseDiagramNameLeft.at(k)="Cu";PhaseDiagramNameRight.at(k)="Pt";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Cu_GS_);
	    ventries.push_back("CuPt_{7}");
	    ventries.push_back("L1_{2}");
	    ventries.push_back("L1_{1}      ");
	    ventries.push_back("                CdPt_{3}^{*}");
	    ventries.push_back("    CuPt_{7}");
	    ventries.push_back(_Pt_GS_);
	  }
	// *************************************************************** CuRh
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Cu","Rh")) {
	    PhaseDiagramNameLeft.at(k)="Cu";PhaseDiagramNameRight.at(k)="Rh";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Cu_GS_);
	    ventries.push_back("CuPt_{7}");
	    ventries.push_back(_Rh_GS_);
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  FeHf
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Fe","Hf")) {
	    PhaseDiagramNameLeft.at(k)="Fe";PhaseDiagramNameRight.at(k)="Hf";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Fe_GS_; // Fe
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/6,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C15_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** FeIr
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Fe","Ir")) {
	    PhaseDiagramNameLeft.at(k)="Fe";PhaseDiagramNameRight.at(k)="Ir";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Fe_GS_);
	    ventries.push_back("L6_{0}");
	    ventries.push_back("NbP");
	    ventries.push_back("    D0_{22}");
	    ventries.push_back(_Ir_GS_);
	  }
	// ***************************************************************  FePd
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Fe","Pd")) {
	    PhaseDiagramNameLeft.at(k)="Fe";PhaseDiagramNameRight.at(k)="Pd";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Fe_GS_);
	    ventries.push_back("   CuZr_{2}");
	    ventries.push_back("    D0_{23}");
	    ventries.push_back("           HfPd_{5}^{*}");
	    ventries.push_back("     Pt_{8}Ti");
	    ventries.push_back(_Pd_GS_);
	  }
	// *************************************************************** FePt
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Fe","Pt")) {
	    PhaseDiagramNameLeft.at(k)="Fe";PhaseDiagramNameRight.at(k)="Pt";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Fe_GS_);
	    ventries.push_back("L1_{0}");
	    ventries.push_back("Ga_{2}Hf");
	    ventries.push_back("             tet-L1_{2} c/a=.992");
	    ventries.push_back("         HfPd_5^{*}");
	    ventries.push_back(_Pt_GS_);
	  }
	// *************************************************************** FeRh
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Fe","Rh")) {
	    PhaseDiagramNameLeft.at(k)="Fe";PhaseDiagramNameRight.at(k)="Rh";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Fe_GS_);
	    ventries.push_back(_STR78+"    ");
	    ventries.push_back("Fe_{2}Rh^{*}     ");
	    ventries.push_back("      D0_{24}");
	    ventries.push_back(_Rh_GS_);
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  GaHf
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ga","Hf")) {
	    PhaseDiagramNameLeft.at(k)="Ga";PhaseDiagramNameRight.at(k)="Hf";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ga_GS_; // Ga
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="     Al_3Zr_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 5/8,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D8_8";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="       C16";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  GeZr
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ge","Zr")) {
	    PhaseDiagramNameLeft.at(k)="Ge";PhaseDiagramNameRight.at(k)="Zr";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ge_GS_; // Ge
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zr_GS_; // Zr
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  HfHg
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","Hg")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="Hg";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="  C11_b/CuZr_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hg_GS_; // Hg
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  HfIn
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","In")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="In";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_2";
            if(isequal(ZConcentrations.at(k).at(j),(double) 4/7,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="In_4Ti_3";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_In_GS_; // In
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  HfIr
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","Ir")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="Ir";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Hf_GS_);
	    ventries.push_back("C37      ");
	    ventries.push_back("       Ir_{5}Zr_{3}");
	    ventries.push_back("B27");
	    ventries.push_back("     Ga_{2}Hf");
	    ventries.push_back("       L1_{2}");
	    ventries.push_back(_Ir_GS_);
	  }
	// ***************************************************************  HfMg
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","Mg")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="Mg";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="CdTi";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Mg_GS_; // Mg
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  HfMn
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","Mn")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="Mn";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
            if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C14";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Mn_GS_; // Mn
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  HfMo
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","Mo")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="Mo";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
            if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Mo_GS_; // Mo
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  HfNi
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","Ni")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="Ni";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="TlI";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{24}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ni_GS_; // Ni
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  HfOs
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","Os")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="Os";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Os_GS_; // Os
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  HfPb
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","Pb")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="Pb";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/6,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="Hf_5Pb^*";
            if(isequal(ZConcentrations.at(k).at(j),(double) 3/8,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="     W_5Si_3";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pb_GS_; // Pb
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** HfPd
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","Pd")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="Pd";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Hf_GS_);
	    ventries.push_back("C11_b/CuZr_2            ");
	    ventries.push_back("B33       ");
	    ventries.push_back("Pd_{3}Ti_{2}                        ");
	    ventries.push_back("Pd_{5}Ti_{3}       ");
	    ventries.push_back("     D0_{24}");
	    ventries.push_back("           HfPd_{5}^{*}");
	    ventries.push_back("    Pt_{8}Ti");
	    ventries.push_back(_Pd_GS_);
	  }
	// ***************************************************************  HfPt
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","Pt")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="Pt";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="NiTi_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B33/TlI";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{24}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 8/9,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="Pt_8Ti    ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pt_GS_; // Pt
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  HfRe
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","Re")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="Re";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="Mo_3Ti^*    ";
            if(isequal(ZConcentrations.at(k).at(j),(double) 3/7,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="Al_3Zr_4   "; // will not show up
            if(isequal(ZConcentrations.at(k).at(j),(double) 25/46,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="       Re_{25}Zr_{21}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C14";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 24/29,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="     Re_{24}Ti_5";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Re_GS_; // Re
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  HfRh
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","Rh")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="Rh";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="CuZr_2";
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B27    ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 5/8,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="Ge_3Rh_5";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Rh_GS_; // Rh
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** HfRu
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","Ru")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="Ru";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Hf_GS_);
	    ventries.push_back("B2");
	    ventries.push_back(_Ru_GS_);
	  }
	// ***************************************************************  HfSc
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","Sc")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="Sc";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/6,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="Hf_5Sc^*       "; //
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="            Hf_3Sc^*"; //
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Sc_GS_; // Sc
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  HfSn
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","Sn")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="Sn";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
            if(isequal(ZConcentrations.at(k).at(j),(double) 3/8,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D8_8      "; //
	    if(isequal(ZConcentrations.at(k).at(j),(double) 4/9,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="      Ga_4Ti_5"; //
            if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C40"; //
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Sn_GS_; // Sn
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  HfTc
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","Tc")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="Tc";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="Mo_3Ti^*     "; //
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C49"; //
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B2"; //
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C14"; //
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Tc_GS_; // Tc
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  HfTi
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="Ti";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C32    "; //
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  HfTl
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","Tl")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="Tl";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="           "+_STR6_BETA2; //
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Tl_GS_; // Tl
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  HfV
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","V")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="V";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_V_GS_; // V
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  HfW
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","W")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="W";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
            if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C15     "; //
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_W_GS_; // W
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  HfZn
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","Zn")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="Zn";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
            if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b/CuZr_2   "; //
            if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="CaIn_2     "; //
            if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="YCd_3  "; //
            if(isequal(ZConcentrations.at(k).at(j),(double) 22/23,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="Zn_{22}Zr         "; //
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zn_GS_; // Zn
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  HfZr
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hf","Zr")) {
	    PhaseDiagramNameLeft.at(k)="Hf";PhaseDiagramNameRight.at(k)="Zr";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Hf_GS_; // Hf
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C32     "; //
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zr_GS_; // Zr
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  HgPd
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hg","Pd")) {
	    PhaseDiagramNameLeft.at(k)="Hg";PhaseDiagramNameRight.at(k)="Pd";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Hg_GS_);
	    ventries.push_back("Hg_{4}Pt");
	    ventries.push_back("Hg_{2}Pt    ");
	    ventries.push_back("L1_{0}      ");
	    ventries.push_back("Ga_{3}Pt_{5}  ");
	    ventries.push_back("    C37");
	    ventries.push_back("     D0_{22}");
	    ventries.push_back("     D1_{a}");
	    ventries.push_back(_Pd_GS_);
	  }
	// *************************************************************** HgPt
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hg","Pt")) {
	    PhaseDiagramNameLeft.at(k)="Hg";PhaseDiagramNameRight.at(k)="Pt";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Hg_GS_);
	    ventries.push_back("Hg_{4}Pt");
	    ventries.push_back(_Pt_GS_);
	  }
	// *************************************************************** HgRh
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Hg","Rh")) {
	    PhaseDiagramNameLeft.at(k)="Hg";PhaseDiagramNameRight.at(k)="Rh";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Hg_GS_);
	    ventries.push_back("Hg_{4}Pt");
	    ventries.push_back(_Rh_GS_);
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  InRe
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"In","Re")) {
	    PhaseDiagramNameLeft.at(k)="In";PhaseDiagramNameRight.at(k)="Re";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_In_GS_; // In
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Re_GS_; // Re
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  InRh
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"In","Rh")) {
	    PhaseDiagramNameLeft.at(k)="In";PhaseDiagramNameRight.at(k)="Rh";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_In_GS_; // In
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Rh_GS_; // Rh
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  InRu
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"In","Ru")) {
	    PhaseDiagramNameLeft.at(k)="In";PhaseDiagramNameRight.at(k)="Ru";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_In_GS_; // In
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ru_GS_; // Ru
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  InSc
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"In","Sc")) {
	    PhaseDiagramNameLeft.at(k)="In";PhaseDiagramNameRight.at(k)="Sc";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_In_GS_; // In
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Sc_GS_; // Sc
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  InSi
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"In","Si")) {
	    PhaseDiagramNameLeft.at(k)="In";PhaseDiagramNameRight.at(k)="Si";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_In_GS_; // In
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Si_GS_; // Si
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  InSn
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"In","Sn")) {
	    PhaseDiagramNameLeft.at(k)="In";PhaseDiagramNameRight.at(k)="Sn";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_In_GS_; // In
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Sn_GS_; // Sn
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  InSr
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"In","Sr")) {
	    PhaseDiagramNameLeft.at(k)="In";PhaseDiagramNameRight.at(k)="Sr";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_In_GS_; // In
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Sr_GS_; // Sr
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  InTa
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"In","Ta")) {
	    PhaseDiagramNameLeft.at(k)="In";PhaseDiagramNameRight.at(k)="Ta";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_In_GS_; // In
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ta_GS_; // Ta
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  InTc
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"In","Tc")) {
	    PhaseDiagramNameLeft.at(k)="In";PhaseDiagramNameRight.at(k)="Tc";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_In_GS_; // In
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Tc_GS_; // Tc
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  InTi
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"In","Ti")) {
	    PhaseDiagramNameLeft.at(k)="In";PhaseDiagramNameRight.at(k)="Ti";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_In_GS_; // In
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// *************************************************************** IrMn
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ir","Mn")) {
	    PhaseDiagramNameLeft.at(k)="Ir";PhaseDiagramNameRight.at(k)="Mn";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ir_GS_);
	    ventries.push_back("L1_{2}");
	    ventries.push_back("B19");
	    ventries.push_back("C37");
	    ventries.push_back("      L6_{0}");
	    ventries.push_back(_Ru_GS_);
	  }
	
	// *************************************************************** IrMo
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ir","Mo")) {
	    PhaseDiagramNameLeft.at(k)="Ir";PhaseDiagramNameRight.at(k)="Mo";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ir_GS_);
	    ventries.push_back("D0_{19}");
	    ventries.push_back("C37");
	    ventries.push_back("B19");
	    ventries.push_back(_Mo_GS_);
	  }
	// *************************************************************** IrNb
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ir","Nb")) {
	    PhaseDiagramNameLeft.at(k)="Ir";PhaseDiagramNameRight.at(k)="Nb";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ir_GS_);
	    ventries.push_back("Co_{3}V   ");
	    ventries.push_back("             L1_{0}");
	    ventries.push_back("         \\sigma_{ABBAB}");
	    ventries.push_back("A15");
	    ventries.push_back(_Nb_GS_);
	  }
	// *************************************************************** IrNi
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ir","Ni")) {
	    PhaseDiagramNameLeft.at(k)="Ir";PhaseDiagramNameRight.at(k)="Ni";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ir_GS_);
	    ventries.push_back("NbP");
	    ventries.push_back(_Ni_GS_);
	  }
	// ***************************************************************  IrOs  
        if(LIB_MODE==LIB_MODEX)  
          if(strsub(alloys.at(k),"Ir","Os")) {  
            PhaseDiagramNameLeft.at(k)="Ir";PhaseDiagramNameRight.at(k)="Os";
            if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
            if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ir_GS_; // Ir  
            if(isequal(ZConcentrations.at(k).at(j),(double) 0.111111,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="Pt_8Ti";
            if(isequal(ZConcentrations.at(k).at(j),(double) 0.833333,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="Hf_5Sc^{*}";
            if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Os_GS_; // Os  
            if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
          }  
	// *************************************************************** IrRe
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ir","Re")) {
	    PhaseDiagramNameLeft.at(k)="Ir";PhaseDiagramNameRight.at(k)="Re";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ir_GS_);
	    ventries.push_back("Pt_{8}Ti");
	    ventries.push_back("Ir_{2}Tc^{*}");
	    ventries.push_back("B19");
	    ventries.push_back("D0_{19}");
	    ventries.push_back(_Re_GS_);
	  }
	// *************************************************************** IrRh
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ir","Rh")) {
	    PhaseDiagramNameLeft.at(k)="Ir";PhaseDiagramNameRight.at(k)="Rh";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ir_GS_);
	    ventries.push_back(_STR18_W3);
	    ventries.push_back("Pd_{2}Ti     ");
	    ventries.push_back("                 "+_STR17_W2);
	    ventries.push_back("    Pd_{2}Ti");
	    ventries.push_back(_Rh_GS_);
	  }
	// *************************************************************** IrRu
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ir","Ru")) {
	    PhaseDiagramNameLeft.at(k)="Ir";PhaseDiagramNameRight.at(k)="Ru";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ir_GS_);
	    ventries.push_back("Pt_{8}Ti");
	    ventries.push_back("L1_{2}");
	    ventries.push_back("B19");
	    ventries.push_back("Ir_{2}Tc^{*}             ");
	    ventries.push_back("     D0_{19}");
	    ventries.push_back("          Hf_{5}Sc^{*}");
	    ventries.push_back(_Ru_GS_);
	  }
	// *************************************************************** IrSc
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ir","Sc")) {
	    PhaseDiagramNameLeft.at(k)="Ir";PhaseDiagramNameRight.at(k)="Sc";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ir_GS_);
	    ventries.push_back("CuPt_{7}");
	    ventries.push_back("C14");
	    ventries.push_back("B2");
	    ventries.push_back("             Ir_{4}Sc_{11}");
	    ventries.push_back("        Mg_{44}Rh_{7}");
	    ventries.push_back(_Sc_GS_);
	  }
	// *************************************************************** IrTa
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ir","Ta")) {
	    PhaseDiagramNameLeft.at(k)="Ir";PhaseDiagramNameRight.at(k)="Ta";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ir_GS_);
	    ventries.push_back("Co_{3}V");
	    ventries.push_back("                  Ga_{2}Hf");
	    ventries.push_back("         L1_{0}");
	    ventries.push_back("         \\sigma_{ABBAB}");
	    ventries.push_back("A15");
	    ventries.push_back(_Ta_GS_);
	  }
	// *************************************************************** IrTc
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ir","Tc")) {
	    PhaseDiagramNameLeft.at(k)="Ir";PhaseDiagramNameRight.at(k)="Tc";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ir_GS_);
	    ventries.push_back("Pt_{8}Ti");
	    ventries.push_back("Ir_{2}Tc^{*}");
	    ventries.push_back("B19");
	    ventries.push_back("      D0_{19}");
	    ventries.push_back(_Tc_GS_);
	  }
	// *************************************************************** IrTi
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ir","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Ir";PhaseDiagramNameRight.at(k)="Ti";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ir_GS_);
	    ventries.push_back("CuPt_{7}");
	    ventries.push_back("L1_{2}");
	    ventries.push_back("Ga_{2}Hf            ");
	    ventries.push_back("    Ga_{3}Pt_{5}");
	    ventries.push_back("L1_{0}");
	    ventries.push_back("C11_{b}");
	    ventries.push_back("A15");
	    ventries.push_back(_Ti_GS_);
	  }
	// ***************************************************************  IrTi
	if(0 && LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ir","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Ir";PhaseDiagramNameRight.at(k)="Ti";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ir_GS_; // Ir
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/8,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="Ca_7Ge";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_0";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="  C11_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="  A15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** IrV
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ir","V")) {
	    PhaseDiagramNameLeft.at(k)="Ir";PhaseDiagramNameRight.at(k)="V";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ir_GS_);
	    ventries.push_back("D0_{19}");
	    ventries.push_back("    L1_{0}"); // no gamma-IrV /\\gamma-IrV");
	    ventries.push_back("A15");
	    ventries.push_back("Pt_{8}Ti");
	    ventries.push_back(_V_GS_);
	  }
	// *************************************************************** IrW
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ir","W")) {
	    PhaseDiagramNameLeft.at(k)="Ir";PhaseDiagramNameRight.at(k)="W";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ir_GS_);
	    ventries.push_back("Pt_{8}Ti");
	    ventries.push_back("D0_{19}");
	    ventries.push_back("    C37");
	    ventries.push_back("B19");
	    ventries.push_back(_W_GS_);
	  }
	// *************************************************************** IrY
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ir","Y")) {
	    PhaseDiagramNameLeft.at(k)="Ir";PhaseDiagramNameRight.at(k)="Y";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ir_GS_);
	    ventries.push_back("C15");
	    ventries.push_back("B2");
	    ventries.push_back("    Pu_{5}Rh_{3}");
	    ventries.push_back("        C_{2}Mn_{5}");
	    ventries.push_back("        D0_{11}");
	    ventries.push_back(_Y_GS_);
	  }
	// *************************************************************** IrZn
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ir","Zn")) {
	    PhaseDiagramNameLeft.at(k)="Ir";PhaseDiagramNameRight.at(k)="Zn";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ir_GS_);
	    ventries.push_back("IrZn^{+}");
	    ventries.push_back("C49");
	    ventries.push_back("      NbPd_{3}");
	    ventries.push_back("              Ir_{2}Zn_{11}");
	    ventries.push_back(_Zn_GS_);
	  }
	// *************************************************************** IrZr
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ir","Zr")) {
	    PhaseDiagramNameLeft.at(k)="Ir";PhaseDiagramNameRight.at(k)="Zr";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ir_GS_);
	    ventries.push_back("L1_{2}");
	    ventries.push_back("Ga_{2}Hf   ");
	    ventries.push_back("NiTi");
	    ventries.push_back("         Ir_{3}Zr_{5}");
	    ventries.push_back("C37");
	    ventries.push_back("      SV_{3}");
	    ventries.push_back(_Zr_GS_);
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  LaRe
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"La","Re")) {
	    PhaseDiagramNameLeft.at(k)="La";PhaseDiagramNameRight.at(k)="Re";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_La_GS_; // La
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Re_GS_; // Re
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  LiMg
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Li","Mg")) {
	    PhaseDiagramNameLeft.at(k)="Li";PhaseDiagramNameRight.at(k)="Mg";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Li_GS_; // Li
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D1_a";	    
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C49/C11_b";	    
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B2/MoTi^*";	    
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b";	    
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Mg_GS_; // Mg
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  MgRe
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Mg","Re")) {
	    PhaseDiagramNameLeft.at(k)="Mg";PhaseDiagramNameRight.at(k)="Re";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Mg_GS_; // Mg
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Re_GS_; // Re
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
        // ***************************************************************  MnOs  
        if(LIB_MODE==LIB_MODEX)  
          if(strsub(alloys.at(k),"Mn","Os")) {  
            PhaseDiagramNameLeft.at(k)="Mn";PhaseDiagramNameRight.at(k)="Os";
            if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
            if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=string(_Mn_GS_+string("      ")); // Mn  
            if(isequal(ZConcentrations.at(k).at(j),(double) 0.5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B19";
            if(isequal(ZConcentrations.at(k).at(j),(double) 0.75,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{19}";
            if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Os_GS_; // Os  
            if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
          }  
	// *************************************************************** MnPd
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Mn","Pd")) {
	    PhaseDiagramNameLeft.at(k)="Mn";PhaseDiagramNameRight.at(k)="Pd";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Mn_GS_);
	    ventries.push_back("Ga_{3}Pt_{5}                       ");
	    ventries.push_back("C37");
	    ventries.push_back("    L1_{2}");
	    ventries.push_back("           HfPd_{5}^{*}");
	    ventries.push_back("       Pt_{8}Ti");
	    ventries.push_back(_Pd_GS_);
	  }
	// *************************************************************** MnPt
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Mn","Pt")) {
	    PhaseDiagramNameLeft.at(k)="Mn";PhaseDiagramNameRight.at(k)="Pt";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Mn_GS_);
	    ventries.push_back("D0_{19}");
	    ventries.push_back("Ga_{3}Pt_{5}          ");
	    ventries.push_back("Ga_{2}Hf");
	    ventries.push_back("         L1_{2}");
	    ventries.push_back("    Pt_{8}Ti");
	    ventries.push_back(_Pt_GS_);
	  }
	// *************************************************************** MnRh
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Mn","Rh")) {
	    PhaseDiagramNameLeft.at(k)="Mn";PhaseDiagramNameRight.at(k)="Rh";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Mn_GS_);
	    ventries.push_back("B2");
	    ventries.push_back("     L1_{2}");
	    ventries.push_back("     CuPt_{7}");
	    ventries.push_back(_Rh_GS_);
	  }
	// *************************************************************** MnRu
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Mn","Ru")) {
	    PhaseDiagramNameLeft.at(k)="Mn";PhaseDiagramNameRight.at(k)="Ru";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Mn_GS_);
	    ventries.push_back("Re_{24}Ti_{5}   ");
	    ventries.push_back(_Ru_GS_);
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  MoNb
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Mo","Nb")) {
	    PhaseDiagramNameLeft.at(k)="Nb";PhaseDiagramNameRight.at(k)="Mo";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Mo_GS_; // Mo
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{3}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Nb_GS_; // Nb
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  MoOs  
        if(LIB_MODE==LIB_MODEX)  
          if(strsub(alloys.at(k),"Mo","Os")) {  
            PhaseDiagramNameLeft.at(k)="Mo";PhaseDiagramNameRight.at(k)="Os";
            if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
            if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Mo_GS_; // Mo  
            if(isequal(ZConcentrations.at(k).at(j),(double) 0.75,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{19}";
            if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Os_GS_; // Os  
            if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
          }  
	// ***************************************************************  MoPd
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Mo","Pd")) {
	    PhaseDiagramNameLeft.at(k)="Mo";PhaseDiagramNameRight.at(k)="Pd";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Mo_GS_);
	    ventries.push_back("    MoPt_{2}");
	    ventries.push_back("    D1_{a}");
	    ventries.push_back("       Pt_{8}Ti");
	    ventries.push_back(_Pd_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Mo","Pd")) {
	    PhaseDiagramNameLeft.at(k)="Mo";PhaseDiagramNameRight.at(k)="Pd";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Mo_GS_; // Mo
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="ZrSi_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 4/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="  D1_a";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 5/6,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="MoPd_5^{proto}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pd_GS_; // Pd
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** MoPt
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Mo","Pt")) {
	    PhaseDiagramNameLeft.at(k)="Mo";PhaseDiagramNameRight.at(k)="Pt";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Mo_GS_);
	    ventries.push_back("B19    ");
	    ventries.push_back("      MoPt_{2}");
	    ventries.push_back("  D1_{a}");
	    ventries.push_back("   Pt_{8}Ti");
	    ventries.push_back(_Pt_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Mo","Pt")) {
	    PhaseDiagramNameLeft.at(k)="Mo";PhaseDiagramNameRight.at(k)="Pt";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Mo_GS_; // Mo
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B19";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="MoPt_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{22}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 4/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="  D1_a";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pt_GS_; // Pt
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** MoRh
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Mo","Rh")) {
	    PhaseDiagramNameLeft.at(k)="Mo";PhaseDiagramNameRight.at(k)="Rh";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Mo_GS_);
	    ventries.push_back("B19     ");
	    ventries.push_back("C37           ");
	    ventries.push_back("    D0_{19}");
	    ventries.push_back("   Pt_{8}Ti");
	    ventries.push_back(_Rh_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Mo","Rh")) {
	    PhaseDiagramNameLeft.at(k)="Mo";PhaseDiagramNameRight.at(k)="Rh";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Mo_GS_; // Mo
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B19";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C37           ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="CdMg_3    ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Rh_GS_; // Rh
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	if(0 && LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Mo","Rh")) {
	    PhaseDiagramNameLeft.at(k)="Mo";PhaseDiagramNameRight.at(k)="Rh";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Mo_GS_; // Mo
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B19";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C37           ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="CdMg_3    ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Rh_GS_; // Rh
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** MoRu
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Mo","Ru")) {
	    PhaseDiagramNameLeft.at(k)="Mo";PhaseDiagramNameRight.at(k)="Ru";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Mo_GS_);
	    ventries.push_back("\\sigma_{AABAB}");
	    ventries.push_back(_Ru_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Mo","Ru")) {
	    PhaseDiagramNameLeft.at(k)="Mo";PhaseDiagramNameRight.at(k)="Ru";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Mo_GS_; // Mo
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{19}   ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ru_GS_; // Ru
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  MoTi
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Mo","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Ti";PhaseDiagramNameRight.at(k)="Mo";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="MoTi_3^{proto}/tie    ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_STR62;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="MoTi^{proto}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="  Mo_3Ti^{proto}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 4/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="      D1_a/tie";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 5/6,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="        Mo_5Ti^{proto}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Mo_GS_; // Mo
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	if(LIB_MODE==LIB_MODEX || LIB_MODE==LIB_MODEN)
	  if(strsub(alloys.at(k),"Mo","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Mo";PhaseDiagramNameRight.at(k)="Ti";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Mo_GS_; // Mo
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="MoTi_3^{proto}/tie    ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_STR62;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="MoTi^{proto}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="  Mo_3Ti^{proto}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="      D1_a/tie";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/6,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="        Mo_5Ti^{proto}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  MoZr
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Mo","Zr")) {
	    PhaseDiagramNameLeft.at(k)="Zr";PhaseDiagramNameRight.at(k)="Mo";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zr_GS_; // Zr
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/6,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="Mo_5Ti^{proto}/tie  ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="MoZr_3^{proto}/tie  ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="MoZr^{proto}/tie      ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Mo_GS_; // Mo
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
        // ***************************************************************  NbOs  
        if(LIB_MODE==LIB_MODEX)  
          if(strsub(alloys.at(k),"Nb","Os")) {  
            PhaseDiagramNameLeft.at(k)="Nb";PhaseDiagramNameRight.at(k)="Os";
       	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Nb_GS_);
	    ventries.push_back("HfPd_5^{*}    ");
            ventries.push_back("A15");
            ventries.push_back("                                  \\sigma_{BAABA}");
            ventries.push_back("                               Al_{12}Mg_{17}");
            ventries.push_back("D0_{24}");
            ventries.push_back(_Os_GS_); // Os  
           }  
	// *************************************************************** NbPd
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Nb","Pd")) {
	    PhaseDiagramNameLeft.at(k)="Nb";PhaseDiagramNameRight.at(k)="Pd";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Nb_GS_);
	    ventries.push_back("Nb_{3}Pd^{+}        ");
	    ventries.push_back("CuZr_{2}     ");
	    ventries.push_back("MoPt_{2}                ");
	    ventries.push_back("             NbPd_{3}");
	    ventries.push_back("             HfPd_{5}^{*}");
	    ventries.push_back("       Pt_{8}Ti");
	    ventries.push_back(_Pd_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Nb","Pd")) {
	    PhaseDiagramNameLeft.at(k)="Nb";PhaseDiagramNameRight.at(k)="Pd";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Nb_GS_; // Nb
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_STR64+"       ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="MoPt_2                 ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{22}  ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pd_GS_; // Pd
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** NbPt
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Nb","Pt")) {
	    PhaseDiagramNameLeft.at(k)="Nb";PhaseDiagramNameRight.at(k)="Pt";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Nb_GS_);
	    ventries.push_back("A15");
	    ventries.push_back("L1_{0}     ");
	    ventries.push_back("     MoPt_{2}");
	    ventries.push_back("      NbPt_{3}/D0_a");
	    ventries.push_back("    Pt_{8}Ti");
	    ventries.push_back(_Pt_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Nb","Pt")) {
	    PhaseDiagramNameLeft.at(k)="Nb";PhaseDiagramNameRight.at(k)="Pt";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Nb_GS_; // Nb
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="A15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_0";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="MoPt_2           ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_a";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pt_GS_; // Pt
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** NbRh
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Nb","Rh")) {
	    PhaseDiagramNameLeft.at(k)="Nb";PhaseDiagramNameRight.at(k)="Rh";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Nb_GS_);
	    ventries.push_back("Pt_{8}Ti");
	    ventries.push_back("A15");
	    ventries.push_back("\\sigma_{BAABA}");
	    ventries.push_back("L1_{0}     ");
	    ventries.push_back("       Co_{3}V");
	    ventries.push_back(_Rh_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Nb","Rh")) {
	    PhaseDiagramNameLeft.at(k)="Nb";PhaseDiagramNameRight.at(k)="Rh";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Nb_GS_; // Nb
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="A15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_0";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="Al_3Pu/L1_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Rh_GS_; // Rh hcp ??
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	if(0 && LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Nb","Rh")) {
	    PhaseDiagramNameLeft.at(k)="Nb";PhaseDiagramNameRight.at(k)="Rh";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Nb_GS_; // Nb
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="A15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_0";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="Al_3Pu/L1_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Rh_GS_; // Rh hcp ??
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** NbRu
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Nb","Ru")) {
	    PhaseDiagramNameLeft.at(k)="Nb";PhaseDiagramNameRight.at(k)="Ru";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Nb_GS_);
	    ventries.push_back("Pt_{8}Ti");
	    ventries.push_back("Nb_{5}Ru^{*}     ");
	    ventries.push_back("L6_{0}");
	    ventries.push_back("Ga_{3}Pt_{5}");
	    ventries.push_back("Ga_{3}Pt_{5}");
	    ventries.push_back(_Ru_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Nb","Ru")) {
	    PhaseDiagramNameLeft.at(k)="Nb";PhaseDiagramNameRight.at(k)="Ru";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Nb_GS_; // Nb
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/6,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="Mo_5Ti^{proto}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_3";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C37";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{24}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ru_GS_; // Ru
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  NbTc
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Nb","Tc")) {
	    PhaseDiagramNameLeft.at(k)="Nb";PhaseDiagramNameRight.at(k)="Tc";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Nb_GS_; // Nb
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="Nb_3Tc^{proto}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Tc_GS_; // Tc
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// *************************************************************** NiPd
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ni","Pd")) {
	    PhaseDiagramNameLeft.at(k)="Ni";PhaseDiagramNameRight.at(k)="Pd";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ni_GS_);
	    ventries.push_back("            tet-"+_STR13_Z1+" c/a=2.7"); 
	    // ventries.push_back("              "+_STR13_Z1);
	    ventries.push_back(_Pd_GS_);
	  }
	// *************************************************************** NiPt
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ni","Pt")) {
	    PhaseDiagramNameLeft.at(k)="Ni";PhaseDiagramNameRight.at(k)="Pt";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ni_GS_);
	    ventries.push_back("D0_{22}");
	    ventries.push_back("    L1_{0}");
	    ventries.push_back("      CuZr_{2}");
	    ventries.push_back("      D0_{23}");
	    ventries.push_back(_Pt_GS_);
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
        // ***************************************************************  OsRe  
        if(LIB_MODE==LIB_MODEX)  
          if(strsub(alloys.at(k),"Os","Re")) {  
            PhaseDiagramNameLeft.at(k)="Os";PhaseDiagramNameRight.at(k)="Re";
     	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Os_GS_);
	    ventries.push_back("D0_{19}");
            ventries.push_back("B19");
            ventries.push_back("    Sc_2Zr^*");
            ventries.push_back("     Re_3Ru^*");
            ventries.push_back(_Re_GS_);
 	  }  
        // ***************************************************************  OsRh  
        if(LIB_MODE==LIB_MODEX)  
          if(strsub(alloys.at(k),"Os","Rh")) {  
            PhaseDiagramNameLeft.at(k)="Os";PhaseDiagramNameRight.at(k)="Rh";
       	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Os_GS_);
            ventries.push_back("RhRu^{*}");
            ventries.push_back(_Rh_GS_);
  	  }  
        // ***************************************************************  OsRu  
        if(LIB_MODE==LIB_MODEX)  
          if(strsub(alloys.at(k),"Os","Ru")) {  
            PhaseDiagramNameLeft.at(k)="Os";PhaseDiagramNameRight.at(k)="Ru";
       	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Os_GS_);
	    ventries.push_back("D0_a");
            ventries.push_back("B19");
	    ventries.push_back("  D0_a");
            ventries.push_back("       Hf_5Sc^{*}");
            ventries.push_back(_Ru_GS_);
	  }  
	// ***************************************************************  OsSc
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Os","Sc")) {
	    PhaseDiagramNameLeft.at(k)="Os";PhaseDiagramNameRight.at(k)="Sc";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Os_GS_);
	    ventries.push_back("C14");
	    ventries.push_back("        "+_STR5_BETA1);
	    ventries.push_back("             Ir_4Sc_{11}");
	    ventries.push_back("        Mg_{44}Rh_{7}");
	    ventries.push_back(_Sc_GS_);
	  }
        // ***************************************************************  OsTa  
        if(LIB_MODE==LIB_MODEX)  
          if(strsub(alloys.at(k),"Os","Ta")) {  
            PhaseDiagramNameLeft.at(k)="Os";PhaseDiagramNameRight.at(k)="Ta";
    	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Os_GS_);
	    ventries.push_back("Ga_2Hf      ");
	    ventries.push_back("Al_{12}Mg_{17}                           ");
	    ventries.push_back("\\sigma_{ABBAB}       ");
	    ventries.push_back("A15    ");
            ventries.push_back(_Ta_GS_);
 	  }  
        // ***************************************************************  OsTc  
        if(LIB_MODE==LIB_MODEX)  
          if(strsub(alloys.at(k),"Os","Tc")) {  
            PhaseDiagramNameLeft.at(k)="Os";PhaseDiagramNameRight.at(k)="Tc";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Os_GS_);
	    ventries.push_back("Hf_5Sc^{*}         ");
	    ventries.push_back("D0_{19}");
	    ventries.push_back("B19");
	    ventries.push_back("    D0_{19}");
            ventries.push_back(_Tc_GS_);
	  }  
	// ***************************************************************  OsTi  
        if(LIB_MODE==LIB_MODEX)  
          if(strsub(alloys.at(k),"Os","Ti")) {  
            PhaseDiagramNameLeft.at(k)="Os";PhaseDiagramNameRight.at(k)="Ti";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Os_GS_);
	    ventries.push_back("B2");
	    ventries.push_back("C49");
	    ventries.push_back("       Mo_3Ti^{*}");
            ventries.push_back(_Ti_GS_);
  	  }  
        // ***************************************************************  OsV  
        if(LIB_MODE==LIB_MODEX)  
          if(strsub(alloys.at(k),"Os","V")) {  
            PhaseDiagramNameLeft.at(k)="Os";PhaseDiagramNameRight.at(k)="V";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Os_GS_);
	    ventries.push_back("Re_3Ru^{*}         ");
	    ventries.push_back("Ga_3Pt_5                      ");
	    ventries.push_back("C11_b         ");
	    ventries.push_back("D0_3   ");
	    ventries.push_back("        Mo_5Ti^{*}");
            ventries.push_back(_V_GS_);
 	  }  
	// ***************************************************************  OsY
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Os","Y")) {
	    PhaseDiagramNameLeft.at(k)="Os";PhaseDiagramNameRight.at(k)="Y";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Os_GS_);
	    ventries.push_back("C14");
	    ventries.push_back("    D0_{11}");
	    ventries.push_back(_Y_GS_);
	  }
	// ***************************************************************  OsW  
        if(LIB_MODE==LIB_MODEX)  
          if(strsub(alloys.at(k),"Os","W")) {  
            PhaseDiagramNameLeft.at(k)="Os";PhaseDiagramNameRight.at(k)="W";
   	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Os_GS_);
            ventries.push_back("D0_{19}");
            ventries.push_back(_W_GS_);
 	  } 
	// ***************************************************************  OsZr  
	if(LIB_MODE==LIB_MODEX)  
          if(strsub(alloys.at(k),"Os","Zr")) {  
            PhaseDiagramNameLeft.at(k)="Os";PhaseDiagramNameRight.at(k)="Zr";
     	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Os_GS_);
            ventries.push_back("C14");
            ventries.push_back("B2");
            ventries.push_back("       Ir_4Sc_{11}");
            ventries.push_back("D1_a");
            ventries.push_back(_Zr_GS_);
    	  }  
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  PdPt
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pd","Pt")) {
	    PhaseDiagramNameLeft.at(k)="Pd";PhaseDiagramNameRight.at(k)="Pt";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pd_GS_);
	    ventries.push_back("CuPt_{7}");
	    ventries.push_back("CdPt_3^{*}  ");
	    ventries.push_back("L1_{1}       ");
	    ventries.push_back("         L1_{2}");
	    ventries.push_back("         CuPt_{7}");
	    ventries.push_back(_Pt_GS_);
	  }
	// *************************************************************** PdRe
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pd","Re")) {
	    PhaseDiagramNameLeft.at(k)="Pd";PhaseDiagramNameRight.at(k)="Re";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pd_GS_);
	    ventries.push_back("D0_{19}");
	    ventries.push_back(_Re_GS_);
	  }
	// *************************************************************** PdSc
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pd","Sc")) {
	    PhaseDiagramNameLeft.at(k)="Pd";PhaseDiagramNameRight.at(k)="Sc";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pd_GS_);
	    ventries.push_back("Pt_{8}Ti");
            ventries.push_back("HfPd_{5}^{*}    "); // HfPd$_5^{\star}$
	    ventries.push_back("L1_{2}");
	    ventries.push_back("C37   ");
            ventries.push_back("      Pd_{4}Pu_{3}");
	    ventries.push_back("B2");
	    ventries.push_back("NiTi_{2}");
	    ventries.push_back(_Sc_GS_);
	  }
	// *************************************************************** PdTa
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pd","Ta")) {
	    PhaseDiagramNameLeft.at(k)="Pd";PhaseDiagramNameRight.at(k)="Ta";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pd_GS_);
	    ventries.push_back("Pt_{8}Ti");
	    ventries.push_back("HfPd_{5}^{*}    ");
	    ventries.push_back("D0_{22}");
	    ventries.push_back("                              MoPt_{2}");
	    ventries.push_back("B11");
	    ventries.push_back("       Pt_{8}Ti");
	    ventries.push_back(_Ta_GS_);
	  }
	// *************************************************************** PdTc
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pd","Tc")) {
	    PhaseDiagramNameLeft.at(k)="Pd";PhaseDiagramNameRight.at(k)="Tc";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pd_GS_);
	    ventries.push_back("RhRu^{*}");
	    ventries.push_back("    D0_{19}");
	    ventries.push_back(_Tc_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Pd","Tc")) {
	    PhaseDiagramNameLeft.at(k)="Tc";PhaseDiagramNameRight.at(k)="Pd";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Tc_GS_; // Tc
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{19}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pd_GS_; // Pd
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** PdTi
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pd","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Pd";PhaseDiagramNameRight.at(k)="Ti";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pd_GS_);
	    // [OBSOLETE] ventries.push_back("Pt_{8}Ti");
	    // [OBSOLETE] ventries.push_back("Pd_{5}Sc^{*}        ");
	    ventries.push_back("HfPd_{5}^{*}        ");
	    ventries.push_back("D0_{24}");
	    ventries.push_back("                 Pd_2Ti");
	    ventries.push_back("                                Pd_5Ti_3");
	    ventries.push_back("                                                 Pd_3Ti_2");
	    ventries.push_back("                C11_b/CuZr_{2}");
	    ventries.push_back("A15");
	    ventries.push_back(_Ti_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Pd","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Ti";PhaseDiagramNameRight.at(k)="Pd";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="A15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B19/L1_0    ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="MoPt_2                 ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{24}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pd_GS_; // Pd
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	if(0 && (LIB_MODE==LIB_MODEX || LIB_MODE==LIB_MODEN))
	  if(strsub(alloys.at(k),"Pd","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Pd";PhaseDiagramNameRight.at(k)="Ti";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="A15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B19/L1_0";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="       MoPt_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{24}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pd_GS_; // Pd
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** PdV
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pd","V")) {
	    PhaseDiagramNameLeft.at(k)="Pd";PhaseDiagramNameRight.at(k)="V";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pd_GS_);
	    ventries.push_back("Pt_{8}Ti");
	    ventries.push_back("NbPd_{3}      ");
	    ventries.push_back("                         MoPt_{2}");
	    ventries.push_back("         Mo_{5}Ti^{*}");
	    ventries.push_back("        Pt_{8}Ti");
	    ventries.push_back(_V_GS_);
	  }
	// *************************************************************** PdW
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pd","W")) {
	    PhaseDiagramNameLeft.at(k)="Pd";PhaseDiagramNameRight.at(k)="W";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pd_GS_);
	    ventries.push_back("Pt_{8}Ti");
	    ventries.push_back(_W_GS_);
	  }
	// *************************************************************** PdY
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pd","Y")) {
	    PhaseDiagramNameLeft.at(k)="Pd";PhaseDiagramNameRight.at(k)="Y";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pd_GS_);
	    ventries.push_back("CuPt_{7}");
	    ventries.push_back("L1_{2}");
	    ventries.push_back("Pd_{4}Pu_{3}");
	    ventries.push_back("B33");
	    ventries.push_back("C37");
	    ventries.push_back("             D0_{11}");
	    ventries.push_back(_Y_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Pd","Y")) {
	    PhaseDiagramNameLeft.at(k)="Y";PhaseDiagramNameRight.at(k)="Pd";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Y_GS_; // Y
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{11}  ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C37  ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B27/B33";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pd_GS_; // Pd
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** PdZn
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pd","Zn")) {
	    PhaseDiagramNameLeft.at(k)="Pd";PhaseDiagramNameRight.at(k)="Zn";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pd_GS_);
	    ventries.push_back("Pt_{8}Ti");
	    ventries.push_back("C37");
	    ventries.push_back("L1_{0}     ");
	    ventries.push_back("     D0_{22}");
	    ventries.push_back("              Ir_{2}Zn_{11}");
	    ventries.push_back(_Zn_GS_);
	  }
	// ***************************************************************  PdZr
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pd","Zr")) {
	    PhaseDiagramNameLeft.at(k)="Pd";PhaseDiagramNameRight.at(k)="Zr";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pd_GS_);
	    ventries.push_back("Pt_{8}Ti");
	    ventries.push_back("HfPd_{5}^{*}");
	    ventries.push_back("D0_{24}");
	    ventries.push_back("B33");
	    ventries.push_back("       C11_b/CuZr_{2}");
	    ventries.push_back(_Zr_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Pd","Zr")) {
	    PhaseDiagramNameLeft.at(k)="Zr";PhaseDiagramNameRight.at(k)="Pd";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zr_GS_; // Zr
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_STR15_Z3;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B27/B33       ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b/tie (paw-gga)                                 ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{24}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 5/6,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C2/m #12";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pd_GS_; // Pd
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// *************************************************************** PtRe
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pt","Re")) {
	    PhaseDiagramNameLeft.at(k)="Pt";PhaseDiagramNameRight.at(k)="Re";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pt_GS_);
	    ventries.push_back(_STR76+"    ");
	    ventries.push_back("    D0_{19}");
	    ventries.push_back(_Re_GS_);
	  }
	// *************************************************************** PtRh
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pt","Rh")) {
	    PhaseDiagramNameLeft.at(k)="Pt";PhaseDiagramNameRight.at(k)="Rh";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pt_GS_);
	    ventries.push_back("NbP");
	    ventries.push_back("    Pd_{2}Ti");
	    ventries.push_back("      D0_{22}");
	    ventries.push_back(_Rh_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Pt","Rh")) {
	    PhaseDiagramNameLeft.at(k)="Pt";PhaseDiagramNameRight.at(k)="Rh";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pt_GS_; // Pt
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D1_a";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{22}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_STR23_CH;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C49";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{22}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 4/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D1_a";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Rh_GS_; // Rh
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	if(0 && LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pt","Rh")) {
	    PhaseDiagramNameLeft.at(k)="Pt";PhaseDiagramNameRight.at(k)="Rh";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pt_GS_; // Pt
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D1_a";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{22}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_STR23_CH;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C49";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{22}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 4/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D1_a";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Rh_GS_; // Rh
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  PtRu
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Pt","Ru")) {
	    PhaseDiagramNameLeft.at(k)="Pt";PhaseDiagramNameRight.at(k)="Ru";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pt_GS_; // Pt
	    // if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_STR13_Z1;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_STR15_Z3+"/tie   ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_STR14_Z2;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ru_GS_; // Ru
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	if(0 && LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pt","Ru")) {
	    PhaseDiagramNameLeft.at(k)="Pt";PhaseDiagramNameRight.at(k)="Ru";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pt_GS_; // Pt
	    // if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_STR13_Z1;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_STR15_Z3+"/tie   ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_STR14_Z2;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ru_GS_; // Ru
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pt","Ru")) {
	    PhaseDiagramNameLeft.at(k)="Pt";PhaseDiagramNameRight.at(k)="Ru";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pt_GS_);
	    ventries.push_back("CdTi");
	    ventries.push_back(_Ru_GS_);
	  }
	// *************************************************************** PtSc
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pt","Sc")) {
	    PhaseDiagramNameLeft.at(k)="Pt";PhaseDiagramNameRight.at(k)="Sc";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pt_GS_);
	    ventries.push_back("Pt_{8}Ti");
	    ventries.push_back("L1_{2}");
	    ventries.push_back("Ga_{2}Hf     ");
            ventries.push_back("Pd_{4}Pu_{3}");
	    ventries.push_back("B2");
	    ventries.push_back("    Cl_{2}Pb");
	    ventries.push_back("           Rh_{13}Sc_{57}");
	    ventries.push_back(_Sc_GS_);
	  }
	// *************************************************************** PtTa
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pt","Ta")) {
	    PhaseDiagramNameLeft.at(k)="Pt";PhaseDiagramNameRight.at(k)="Ta";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pt_GS_);
	    ventries.push_back("Pt_{8}Ti");
	    ventries.push_back("NbPt_{3}");
	    ventries.push_back("            Au_{2}V");
	    ventries.push_back("          L1_{0}");
	    ventries.push_back("  {\\sigma_{BBBAB}}");
	    ventries.push_back("           {{{A15}}}");
	    ventries.push_back("   Pt_{8}Ti");
	    ventries.push_back(_Ta_GS_);
	  }
	// ***************************************************************  PtTc
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pt","Tc")) {
	    PhaseDiagramNameLeft.at(k)="Pt";PhaseDiagramNameRight.at(k)="Tc";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pt_GS_);
	    //  ventries.push_back(_STR13_Z1+"     ");
	    ventries.push_back(_STR76+"     ");
	    ventries.push_back("Ir_{2}Tc^{*}");
	    ventries.push_back("      D0_{19}");
	    ventries.push_back(_Tc_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Pt","Tc")) {
	    PhaseDiagramNameLeft.at(k)="Pt";PhaseDiagramNameRight.at(k)="Tc";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pt_GS_; // Pt
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_STR13_Z1;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{19}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Tc_GS_; // Tc
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** PtTi
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pt","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Pt";PhaseDiagramNameRight.at(k)="Ti";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pt_GS_);
	    ventries.push_back("Pt_{8}Ti");
	    ventries.push_back("HfPd_{5}^{*}");
	    ventries.push_back("PuAl_{3}");
	    //   ventries.push_back("C49     ");
	    //   ventries.push_back("     Pd_3Ti_2");
	    ventries.push_back("     Au_{5}GaZn_{2}");
	    ventries.push_back("NiTi");
	    ventries.push_back("A15");
	    ventries.push_back(_Ti_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Pt","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Pt";PhaseDiagramNameRight.at(k)="Ti";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pt_GS_; // Pt
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{24}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C49/C37";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B19";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="A15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	if(0 && (LIB_MODE==LIB_MODEX || LIB_MODE==LIB_MODEN))
	  if(strsub(alloys.at(k),"Pt","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Pt";PhaseDiagramNameRight.at(k)="Ti";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pt_GS_; // Pt
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{24}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C49/C37";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B19";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="A15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** PtV
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pt","V")) {
	    PhaseDiagramNameLeft.at(k)="Pt";PhaseDiagramNameRight.at(k)="V";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pt_GS_);
	    ventries.push_back("Pt_{8}Ti");
	    ventries.push_back("D0_{a}");
	    ventries.push_back("MoPt_{2}   ");
	    ventries.push_back("           L1_{0}");
	    ventries.push_back("A15");
	    ventries.push_back("   Pt_{8}Ti");
	    ventries.push_back(_V_GS_);
	  }

	// *************************************************************** PtW
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pt","W")) {
	    PhaseDiagramNameLeft.at(k)="Pt";PhaseDiagramNameRight.at(k)="W";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pt_GS_);
	    ventries.push_back("Pt_{8}Ti");
	    ventries.push_back("D1_{a}");
	    ventries.push_back("MoPt_{2}");
	    ventries.push_back(_W_GS_);
	  }
	// ***************************************************************  PtY
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pt","Y")) {
	    PhaseDiagramNameLeft.at(k)="Pt";PhaseDiagramNameRight.at(k)="Y";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pt_GS_);
	    ventries.push_back("D2_{d}");
	    ventries.push_back("L1_{2}");
	    ventries.push_back("C15");
            ventries.push_back("Pd_{4}Pu_{3}");
	    ventries.push_back("B33");
	    ventries.push_back("    Cl_{2}Pb");
	    ventries.push_back("    D0_{11}");
	    ventries.push_back(_Y_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Pt","Y")) {
	    PhaseDiagramNameLeft.at(k)="Pt";PhaseDiagramNameRight.at(k)="Y";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pt_GS_; // Pt
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B33";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 5/8,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=" D8_8";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C37";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{11}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Y_GS_; // Y
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** PtZn
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pt","Zn")) {
	    PhaseDiagramNameLeft.at(k)="Pt";PhaseDiagramNameRight.at(k)="Zn";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pt_GS_);
	    ventries.push_back("CuPt_{7}");
	    ventries.push_back("CdPt_3^{*}       ");
	    ventries.push_back("L1_{0}");
	    ventries.push_back("     Pt_{7}Zn_{12}");
	    ventries.push_back("      C49");
	    ventries.push_back("    D0_{22}");
	    ventries.push_back("           Ir_{2}Zn_{11}");
	    ventries.push_back(_Zn_GS_);
	  }

	// *************************************************************** PtZr
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Pt","Zr")) {
	    PhaseDiagramNameLeft.at(k)="Pt";PhaseDiagramNameRight.at(k)="Zr";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Pt_GS_);
	    ventries.push_back("Pt_{8}Ti");
	    ventries.push_back("D0_{24}");
	    ventries.push_back("B33");
	    ventries.push_back("C16");
	    ventries.push_back(_Zr_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Pt","Zr")) {
	    PhaseDiagramNameLeft.at(k)="Pt";PhaseDiagramNameRight.at(k)="Zr";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Pt_GS_; // Pt
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{24}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B33";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C16";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="A15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zr_GS_; // Zr
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// *************************************************************** ReRh
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Re","Rh")) {
	    PhaseDiagramNameLeft.at(k)="Re";PhaseDiagramNameRight.at(k)="Rh";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Re_GS_);
	    ventries.push_back("D0_{19}");
	    ventries.push_back("B19");
	    ventries.push_back("             Ir_{2}Tc^{*}");
	    ventries.push_back(_Rh_GS_);
	  }
	// *************************************************************** ReRu
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Re","Ru")) {
	    PhaseDiagramNameLeft.at(k)="Re";PhaseDiagramNameRight.at(k)="Ru";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Re_GS_);
	    ventries.push_back("Re_{3}Ru^{*}      ");
	    ventries.push_back("B19");
	    ventries.push_back("   D0_{19}");
	    ventries.push_back(_Ru_GS_);
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  RhRu
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Rh","Ru")) {
	    PhaseDiagramNameLeft.at(k)="Rh";PhaseDiagramNameRight.at(k)="Ru";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Rh_GS_);
	    ventries.push_back("Pt_{8}Ti");
	    ventries.push_back("RhRu^{*}");
	    ventries.push_back("           RhRu_{2}^{*}");
	    ventries.push_back("            RhRu_{5}^{*}");
	    ventries.push_back(_Ru_GS_);
	  }
	if(0 && LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Rh","Ru")) {
	    PhaseDiagramNameLeft.at(k)="Rh";PhaseDiagramNameRight.at(k)="Ru";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Rh_GS_; // Rh
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="RhRu^{proto}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="RhRu_2^{proto}";
 	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="127";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ru_GS_; // Ru
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Rh","Ru")) {
	    PhaseDiagramNameLeft.at(k)="Ru";PhaseDiagramNameRight.at(k)="Rh";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ru_GS_; // Ru
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="127";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="RhRu_2^{proto}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="RhRu^{proto}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Rh_GS_; // Rh
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** RhSc
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Rh","Sc")) {
	    PhaseDiagramNameLeft.at(k)="Rh";PhaseDiagramNameRight.at(k)="Sc";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Rh_GS_);
	    ventries.push_back("CuPt_{7}");
	    ventries.push_back("L1_{2}");
	    ventries.push_back("B2");
	    ventries.push_back("          Ir_{4}Sc_{11}");
	    ventries.push_back("          Rh_{13}Sc_{57}");
	    ventries.push_back("        Mg_{44}Rh_{7}");
	    ventries.push_back(_Sc_GS_);
	  }
	// *************************************************************** RhTa
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Rh","Ta")) {
	    PhaseDiagramNameLeft.at(k)="Rh";PhaseDiagramNameRight.at(k)="Ta";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Rh_GS_);
	    ventries.push_back("L1_{2}");
	    ventries.push_back("                         Ga_{2}Hf");
	    ventries.push_back("A15");
	    ventries.push_back("           RuTc_{5}^{*}");
	    ventries.push_back("       Pt_{8}Ti");
	    ventries.push_back(_Ta_GS_);
	  }
	// *************************************************************** RhTc
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Rh","Tc")) {
	    PhaseDiagramNameLeft.at(k)="Rh";PhaseDiagramNameRight.at(k)="Tc";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Rh_GS_);
	    ventries.push_back("Ir_{2}Tc^{*}");
	    ventries.push_back("B19");
	    ventries.push_back("     D0_{19}");
	    ventries.push_back(_Tc_GS_);
	  }
	if(0 && LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Rh","Tc")) {
	    PhaseDiagramNameLeft.at(k)="Rh";PhaseDiagramNameRight.at(k)="Tc";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Rh_GS_; // Rh
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="ZrSi_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B19";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{19}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Tc_GS_; // Tc
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Rh","Tc")) {
	    PhaseDiagramNameLeft.at(k)="Tc";PhaseDiagramNameRight.at(k)="Rh";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Tc_GS_; // Tc
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{19}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B19";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="ZrSi_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Rh_GS_; // Rh
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** RhTi
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Rh","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Rh";PhaseDiagramNameRight.at(k)="Ti";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Rh_GS_);
	    ventries.push_back("CuPt_{7}");
	    ventries.push_back("L1_{2}");
	    ventries.push_back("Ge_{3}Rh_{5}     ");
	    ventries.push_back("        L1_{0}");
	    ventries.push_back("        C11_{b}");
	    ventries.push_back(_Ti_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Rh","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Ti";PhaseDiagramNameRight.at(k)="Rh";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b/CuZr_2   ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_0";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C37";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Rh_GS_; // Rh
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** RhV
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Rh","V")) {
	    PhaseDiagramNameLeft.at(k)="Rh";PhaseDiagramNameRight.at(k)="V";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Rh_GS_);
	    ventries.push_back("HfPd_{5}^{*}  ");
	    ventries.push_back("D0_{19}");
	    ventries.push_back("          L1_{0}");
	    ventries.push_back("A15");
	    ventries.push_back("          RuTc_{5}^{*}");
	    ventries.push_back("    Pt_{8}Ti");
	    ventries.push_back(_V_GS_);
	  }
	// *************************************************************** RhW
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Rh","W")) {
	    PhaseDiagramNameLeft.at(k)="Rh";PhaseDiagramNameRight.at(k)="W";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Rh_GS_);
	    ventries.push_back("Pt_{8}Ti");
	    ventries.push_back("D0_{19}");
	    ventries.push_back("          C37");
	    ventries.push_back(_W_GS_);
	  }
	// *************************************************************** RhY
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Rh","Y")) {
	    PhaseDiagramNameLeft.at(k)="Rh";PhaseDiagramNameRight.at(k)="Y";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Rh_GS_);
	    ventries.push_back("CeNi_{3}");
	    ventries.push_back("C15");
	    ventries.push_back("B2");
	    ventries.push_back("          Pu_{5}Rh_{3}");
	    ventries.push_back("          Fe_{3}Th_{7}");
	    ventries.push_back("     D0_{11}");
	    ventries.push_back(_Y_GS_);
	  }
	if(0 && LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Rh","Y")) {
	    PhaseDiagramNameLeft.at(k)="Rh";PhaseDiagramNameRight.at(k)="Y";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Rh_GS_; // Rh
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C37";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{11}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Y_GS_; // Y
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Rh","Y")) {
	    PhaseDiagramNameLeft.at(k)="Y";PhaseDiagramNameRight.at(k)="Rh";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Y_GS_; // Y
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{11}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C37";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Rh_GS_; // Rh
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** RhZn
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Rh","Zn")) {
	    PhaseDiagramNameLeft.at(k)="Rh";PhaseDiagramNameRight.at(k)="Zn";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Rh_GS_);
	    ventries.push_back("B2");
	    ventries.push_back("Ga_{3}Pt_{5}        ");
	    ventries.push_back("       ZrSi_{2}");
	    //	    ventries.push_back("                 Ir_{2}Tc^{*}");
	    ventries.push_back("      D0_{23}");
	    ventries.push_back("          Rh_{2}Zn_{11}");
	     ventries.push_back("RhZn_{13}      ");
	    ventries.push_back(_Zn_GS_);
	  }
	// ***************************************************************  RhZr
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Rh","Zr")) {
	    PhaseDiagramNameLeft.at(k)="Rh";PhaseDiagramNameRight.at(k)="Zr";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Rh_GS_);
	    ventries.push_back("L1_{2}");
	    ventries.push_back("Pd_{5}Pu_{3}       ");
            ventries.push_back("                Pd_{4}Pu_{3}");
	    ventries.push_back("       B33");
	    ventries.push_back("     C11_{b}");
	    ventries.push_back("     SV_{3}");
	    ventries.push_back(_Zr_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Rh","Zr")) {
	    PhaseDiagramNameLeft.at(k)="Zr";PhaseDiagramNameRight.at(k)="Rh";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zr_GS_; // Zr
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D1_{a}/tie    ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_STR15_Z3+"/tie  ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C16";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B27";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C37/tie   ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Rh_GS_; // Rh
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// *************************************************************** RuSc
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ru","Sc")) {
	    PhaseDiagramNameLeft.at(k)="Ru";PhaseDiagramNameRight.at(k)="Sc";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ru_GS_);
	    ventries.push_back("C14");
	    ventries.push_back("B2");
	    ventries.push_back("   C11_b");
	    ventries.push_back("          Ir_4Sc_{11}");
	    ventries.push_back("         Mg_{44}Rh_{7}");
	    ventries.push_back(_Sc_GS_);
	  }
	// *************************************************************** RuTa
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ru","Ta")) {
	    PhaseDiagramNameLeft.at(k)="Ru";PhaseDiagramNameRight.at(k)="Ta";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ru_GS_);
	    ventries.push_back("Ga_{3}Pt_{5}");
	    ventries.push_back("          Ga_{3}Pt_{5}");
	    ventries.push_back("          "+_STR13_Z1);
	    ventries.push_back("           Nb_{5}Ru^{*}");
	    ventries.push_back(_Ta_GS_);
	  }
	// *************************************************************** RuTc
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ru","Tc")) {
	    PhaseDiagramNameLeft.at(k)="Ru";PhaseDiagramNameRight.at(k)="Tc";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ru_GS_);
	    ventries.push_back("D0_{19}");
	    ventries.push_back("B19");
	    ventries.push_back("     D0_{19}");
	    ventries.push_back("            RuTc_{5}^{*}");
	    ventries.push_back(_Tc_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Ru","Tc")) {
	    PhaseDiagramNameLeft.at(k)="Tc";PhaseDiagramNameRight.at(k)="Ru";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Tc_GS_; // Tc
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{19}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B19";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{19}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ru_GS_; // Ru
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** RuTi
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ru","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Ru";PhaseDiagramNameRight.at(k)="Ti";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ru_GS_);
	    ventries.push_back("B2");
	    ventries.push_back("C49");
	    ventries.push_back("             Mo_{3}Ti^{*}");
	    ventries.push_back(_Ti_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Ru","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Ti";PhaseDiagramNameRight.at(k)="Ru";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="RuTi_3^{proto}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C49";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ru_GS_; // Ru
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// *************************************************************** RuV
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ru","V")) {
	    PhaseDiagramNameLeft.at(k)="Ru";PhaseDiagramNameRight.at(k)="V";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ru_GS_);
	    ventries.push_back("Re_{3}Ru^{*}      ");
	    ventries.push_back("C37");
	    ventries.push_back("Ga_{3}Pt_{5}               ");
	    ventries.push_back("     C11_{b}");
	    ventries.push_back("          Mo_{3}Ti^{*}");
	    ventries.push_back("    D1_{a}");
	    ventries.push_back("           Nb_{5}Ru^{*}");
	    ventries.push_back("       Pt_{8}Ti");
	    ventries.push_back(_V_GS_);
	  }
	// *************************************************************** RuW
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ru","W")) {
	    PhaseDiagramNameLeft.at(k)="Ru";PhaseDiagramNameRight.at(k)="W";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ru_GS_);
	    ventries.push_back("D0_{19}");
	    ventries.push_back(_W_GS_);
	  }
	// ***************************************************************  RuY
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Ru","Y")) {
	    PhaseDiagramNameLeft.at(k)="Y";PhaseDiagramNameRight.at(k)="Ru";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Y_GS_; // Y
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{11}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C16";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C14";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ru_GS_; // Ru
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  RuY
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ru","Y")) {
	    PhaseDiagramNameLeft.at(k)="Ru";PhaseDiagramNameRight.at(k)="Y";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ru_GS_+"  ");
	    ventries.push_back("C14");
	    ventries.push_back("Ru_{25}Y_{44}   ");
	    ventries.push_back("          C_{2}Mn_{5}");
	    ventries.push_back("      D0_{11}");
	    ventries.push_back(_Y_GS_);
	  }
	// *************************************************************** RuZn
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ru","Zn")) {
	    PhaseDiagramNameLeft.at(k)="Ru";PhaseDiagramNameRight.at(k)="Zn";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ru_GS_);
	    ventries.push_back("    L1_{2}");
	    ventries.push_back("      RuZn_{6}");
	    ventries.push_back(_Zn_GS_);
	  }
	// ***************************************************************  RuZr
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ru","Zr")) {
	    PhaseDiagramNameLeft.at(k)="Ru";PhaseDiagramNameRight.at(k)="Zr";
	    ventries.clear();vjindex.push_back(j);
	    ventries.push_back(_Ru_GS_);
	    ventries.push_back("B2");
	    ventries.push_back(_Zr_GS_);
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Ru","Zr")) {
	    PhaseDiagramNameLeft.at(k)="Zr";PhaseDiagramNameRight.at(k)="Ru";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zr_GS_; // Zr
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D1_{a}   ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ru_GS_; // Ru
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  ScTa
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Sc","Ta")) {
	    PhaseDiagramNameLeft.at(k)="Sc";PhaseDiagramNameRight.at(k)="Ta";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Sc_GS_; // Sc
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ta_GS_; // Ta
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  ScTc
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Sc","Tc")) {
	    PhaseDiagramNameLeft.at(k)="Sc";PhaseDiagramNameRight.at(k)="Tc";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Sc_GS_; // Sc
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Tc_GS_; // Tc
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  ScTi
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Sc","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Sc";PhaseDiagramNameRight.at(k)="Ti";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Sc_GS_; // Sc
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  ScV
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Sc","V")) {
	    PhaseDiagramNameLeft.at(k)="Sc";PhaseDiagramNameRight.at(k)="V";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Sc_GS_; // Sc
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_V_GS_; // V
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  ScW
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Sc","W")) {
	    PhaseDiagramNameLeft.at(k)="Sc";PhaseDiagramNameRight.at(k)="W";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Sc_GS_; // Sc
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_W_GS_; // W
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  ScY
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Sc","Y")) {
	    PhaseDiagramNameLeft.at(k)="Sc";PhaseDiagramNameRight.at(k)="Y";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Sc_GS_; // Sc
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Y_GS_; // Y
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  ScZn
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Sc","Zn")) {
	    PhaseDiagramNameLeft.at(k)="Sc";PhaseDiagramNameRight.at(k)="Zn";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Sc_GS_; // Sc
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zn_GS_; // Zn
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  ScZr
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Sc","Zr")) {
	    PhaseDiagramNameLeft.at(k)="Sc";PhaseDiagramNameRight.at(k)="Zr";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Sc_GS_; // Sc
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zr_GS_; // Zr
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  TaTc
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ta","Tc")) {
	    PhaseDiagramNameLeft.at(k)="Ta";PhaseDiagramNameRight.at(k)="Tc";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ta_GS_; // Ta
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Tc_GS_; // Tc
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  TaTi
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ta","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Ta";PhaseDiagramNameRight.at(k)="Ti";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ta_GS_; // Ta
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  TaV
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ta","V")) {
	    PhaseDiagramNameLeft.at(k)="Ta";PhaseDiagramNameRight.at(k)="V";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ta_GS_; // Ta
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_V_GS_; // V
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  TaW
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ta","W")) {
	    PhaseDiagramNameLeft.at(k)="Ta";PhaseDiagramNameRight.at(k)="W";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ta_GS_; // Ta
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_W_GS_; // W
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  TaY
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ta","Y")) {
	    PhaseDiagramNameLeft.at(k)="Ta";PhaseDiagramNameRight.at(k)="Y";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ta_GS_; // Ta
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Y_GS_; // Y
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  TaZn
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ta","Zn")) {
	    PhaseDiagramNameLeft.at(k)="Ta";PhaseDiagramNameRight.at(k)="Zn";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ta_GS_; // Ta
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zn_GS_; // Zn
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  TaZr
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ta","Zr")) {
	    PhaseDiagramNameLeft.at(k)="Ta";PhaseDiagramNameRight.at(k)="Zr";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ta_GS_; // Ta
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zr_GS_; // Zr
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  TcTi
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Tc","Ti")) {
	    PhaseDiagramNameLeft.at(k)="Ti";PhaseDiagramNameRight.at(k)="Tc";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="TcTi_3^{proto}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C49";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C11_b";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Tc_GS_; // Tc
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  TcY
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Tc","Y")) {
	    PhaseDiagramNameLeft.at(k)="Y";PhaseDiagramNameRight.at(k)="Tc";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Y_GS_; // Y
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D0_{11}";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C14";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Tc_GS_; // Tc
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  TcZr
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Tc","Zr")) {
	    PhaseDiagramNameLeft.at(k)="Zr";PhaseDiagramNameRight.at(k)="Tc";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zr_GS_; // Zr
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D1_{a}   ";
	    //if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="TcZr_3^{proto}/tie    ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C49";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C14";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Tc_GS_; // Tc
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  TiW
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ti","W")) {
	    PhaseDiagramNameLeft.at(k)="Ti";PhaseDiagramNameRight.at(k)="W";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_STR62;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=" ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 4/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D1_a";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 5/6,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=" ";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_W_GS_; // W
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  TiZn
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Ti","Zn")) {
	    PhaseDiagramNameLeft.at(k)="Ti";PhaseDiagramNameRight.at(k)="Zn";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Ti_GS_; // Ti
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="A15";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_0/B2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="C14";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zn_GS_; // Zn
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  VZr
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"V","Zr")) {
	    PhaseDiagramNameLeft.at(k)="V";PhaseDiagramNameRight.at(k)="Zr";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_V_GS_; // V
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zr_GS_; // Zr
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  WZr
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"W","Zr")) {
	    PhaseDiagramNameLeft.at(k)="W";PhaseDiagramNameRight.at(k)="Zr";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_W_GS_; // W
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zr_GS_; // Zr
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  YZr
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Y","Zr")) {
	    PhaseDiagramNameLeft.at(k)="Y";PhaseDiagramNameRight.at(k)="Zr";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Y_GS_; // Y
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zr_GS_; // Zr
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  ZnZr
	if(LIB_MODE==LIB_MODEX)
	  if(strsub(alloys.at(k),"Zn","Zr")) {
	    PhaseDiagramNameLeft.at(k)="Zn";PhaseDiagramNameRight.at(k)="Zr";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zn_GS_; // Zn
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Zr_GS_; // Zr
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ******************************************************************************************************************************
	// *                                                                                                                            *
	// *                                                  Mg-Calphad Project                                                        *
	// *                                                                                                                            *
	// ******************************************************************************************************************************
	// ***************************************************************  Mg-Calphad Project
	// ***************************************************************  MgBi
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Mg","Bi")) {  // FIXED
	    PhaseDiagramNameLeft.at(k)="Mg";PhaseDiagramNameRight.at(k)="Bi";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Mg_GS_; // Mg
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D5_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="A7";
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  MgIn
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Mg","In")) { // FIXED
	    PhaseDiagramNameLeft.at(k)="Mg";PhaseDiagramNameRight.at(k)="In";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Mg_GS_; // Mg
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/3,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="Mg_2In";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_0";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 3/4,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="L1_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_In_GS_; // In
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  MgSb
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Mg","Sb")) { // FIXED
	    PhaseDiagramNameLeft.at(k)="Mg";PhaseDiagramNameRight.at(k)="Sb";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_Mg_GS_; // Mg
	    if(isequal(ZConcentrations.at(k).at(j),(double) 2/5,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="D5_2";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="A7";
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  InBi
	if(strsub(alloys.at(k),"In","Bi")) { // FIXED
	  PhaseDiagramNameLeft.at(k)="In";PhaseDiagramNameRight.at(k)="Bi";
	  if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_In_GS_; // In
	  if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B10";
	  if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="A7";
	  if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	}
	// ***************************************************************  InSb
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"In","Sb")) { // FIXED
	    PhaseDiagramNameLeft.at(k)="In";PhaseDiagramNameRight.at(k)="Sb";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=_In_GS_; // In
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1/2,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="B3";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="A7";
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************  BiSb
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys.at(k),"Bi","Sb")) { // FIXED
	    PhaseDiagramNameLeft.at(k)="Bi";PhaseDiagramNameRight.at(k)="Sb";
	    if(vPRE) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	    if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="A7";
	    if(isequal(ZConcentrations.at(k).at(j),(double) 1,CEPSILON)) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="A7";
	    if(vPOS) cerr << (alloys.at(k)) << " " << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name << endl;
	  }
	// ***************************************************************
	// ***************************************************************

	// fix spaces
	if(ZGNDlibrary.at(k).at(j)>=0) {
	  if(ZConcentrations.at(k).at(j)<0.5) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name=ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name+"  ";
	  if(ZConcentrations.at(k).at(j)>0.5 && ZConcentrations.at(k).at(j)<1.0) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="           "+ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name;
	  if(ZConcentrations.at(k).at(j)>0.99) ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name="  "+ZLibrary.at(k).at(ZGNDlibrary.at(k).at(j)).structure_name;
	}

        if(0) {
	  string speciesA,speciesB;
	  KBIN::VASP_SplitAlloySpecies(alloys.at(k),speciesA,speciesB);
	  KBIN::VASP_PseudoPotential_CleanName(speciesA);
	  KBIN::VASP_PseudoPotential_CleanName(speciesB);	  
	  if(isequal(ZConcentrations.at(k).at(j),(double) 0,CEPSILON)) {
	    cout << "        // ***************************************************************  " << speciesA << "" << speciesB << "  " << endl;
	    cout << "        if(LIB_MODE==LIB_MODEX)" << "  " << endl;
	    cout << "          if(strsub(alloys.at(k),\"" << speciesA << "\",\"" << speciesB << "\")) {" << "  " << endl;
	    cout << "            PhaseDiagramNameLeft.at(k)=\"" << speciesA << "\";PhaseDiagramNameRight.at(k)=\"" << speciesB << "\";" << "  " << endl;
	    cout << "	    ventries.clear();vjindex.push_back(j);" << endl;
	    cout << "	    ventries.push_back(_" << speciesA << "_GS_);" << endl;
	  }
	  if(ZConcentrations.at(k).at(j)>0.0 && ZConcentrations.at(k).at(j)<1.0) {
	    cout << "	    ventries.push_back(\"XXX\");" << endl;
	  }
	  if(ZConcentrations.at(k).at(j)>=1.0) {
	    cout << "	    ventries.push_back(_" << speciesB << "_GS_);" << endl;
	    cout << "          }"  << endl;
	  }
	}  
	// done with prelim
      }
    } // j-cycle 
    // cerr << "ventries.size()=" << ventries.size() << endl;
    // cerr << "vjindex.size()=" << vjindex.size() << endl;
    // cerr << "ZGNDlibrary.at(k).size()=" << ZGNDlibrary.at(k).size() << endl;
    // now fix one alloy
    if(ventries.size()==vjindex.size()) {
      for(uint i=0;i<vjindex.size();i++) {
	// cerr << "i=" << i << "  vjindex.at(i)=" << vjindex.at(i) << "   ZGNDlibrary.at(k).at(vjindex.at(i))=" << ZGNDlibrary.at(k).at(vjindex.at(i)) << endl;
	if(ZGNDlibrary.at(k).at(vjindex.at(i))>=0) 
	  ZLibrary.at(k).at(ZGNDlibrary.at(k).at(vjindex.at(i))).structure_name=ventries.at(i);
      }
    } else {
      cerr << "ventries.size()=" << ventries.size() << endl;
      cerr << "vjindex.size()=" << vjindex.size() << endl;
      cerr << "ZGNDlibrary.at(k).size()=" << ZGNDlibrary.at(k).size() << endl;
    }
  } // k-cycle
  // now fix them all
  
  if(verbose) cerr << "FixGndStatesNamesConcentrations: stop" << endl;
  if(verbose) cerr << "***************************************************************" << endl;
  return TRUE;
}

#endif
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
