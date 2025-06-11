// **************************************************************************
// *                                                                        *
// *  PHASES STEFANO CURTAROLO MASSACHUSETTS INSTITUTE OF TECHNOLOGY 2003   *
// *                                                                        *
// **************************************************************************

/*
  #define C0000000 0.00000000000000  // C[1]
  #define C0111111 0.11111111111111  // C[2]
  #define C0125000 0.12500000000000  // C[3]
  #define C0166666 0.16666666666666  // C[4]
  #define C0200000 0.20000000000000  // C[5]
  #define C0250000 0.25000000000000  // C[6]
  #define C0285714 0.28571428571429  // not used
  #define C0333333 0.33333333333333  // C[7]
  #define C0375000 0.37500000000000  // C[8]
  #define C0400000 0.40000000000000  // C[9]
  #define C0428571 0.42857142857143  // C[10]
  #define C4444444 0.44444444444444  // C[11]
  #define C0500000 0.50000000000000  // C[12]
  #define C5555555 0.55555555555555  // C[13]
  #define C0571428 0.57142857142857  // C[14]
  #define C0600000 0.60000000000000  // C[15]
  #define C0625000 0.62500000000000  // C[16]
  #define C0666666 0.66666666666666  // C[17]
  #define C0750000 0.75000000000000  // C[18]
  #define C0800000 0.80000000000000  // C[19]
  #define C0833333 0.83333333333333  // C[20]
  #define C0875000 0.87500000000000  // C[21]
  #define C0888888 0.88888888888888  // C[22]
*/

// ***************************************************************
// ***************************************************************
// ***************************************************************
// DECLARATIONS

bool strsub(string strzero, string subA);
bool strsub(string strzero, string subA, string subB);
bool strsub(string strzero, string subA, string subB, string subC);
string ElementSymbolToString(string stringin);

#ifdef ZUNGER_PROTO
#define _STR5_BETA1   "\\beta1"
#define _STR6_BETA2   "\\beta2"
#define _STR9_ALPHA1  "\\alpha1"
#define _STR10_ALPHA2 "\\alpha2"
#define _STR13_Z1     "Z1"
#define _STR14_Z2     "Z2"
#define _STR15_Z3     "Z3"
#define _STR17_W2     "W2"
#define _STR19_Y1     "Y1"
#define _STR20_Y2     "Y2"
#define _STR21_Y3     "Y3"
#define _STR23_CH     "CH"
#define _STR27_V1     "V1"
#define _STR28_V2     "V2"
#define _STR29_V3     "V3"
#else
#define _STR5_BETA1   "FCC_{AB2}^{100}"
#define _STR6_BETA2   "FCC_{AB2}^{100}"
#define _STR9_ALPHA1  "FCC_{AB2}^{[111]}"
#define _STR10_ALPH2  "FCC_{AB2}^{[111]}"
#define _STR13_Z1     "FCC_{AB3}^{[001]}"
#define _STR14_Z2     "FCC_{A2B2}^{[001]}"
#define _STR15_Z3     "FCC_{AB3}^{[001]}"
#define _STR17_W2     "FCC_{A2B2}^{[311]}"
#define _STR19_Y1     "FCC_{AB3}^{[011]}"
#define _STR20_Y2     "FCC_{A2B2}^{[011]}"
#define _STR21_Y3     "FCC_{AB3}^{[011]}"
//#define _STR23_CH   "NbAs"
#define _STR23_CH     "FCC_{A2B2}^{[201]}"
#define _STR27_V1     "FCC_{AB3}^{[111]}"
#define _STR28_V2     "FCC_{A2B2}^{[111]}"
#define _STR29_V3     "FCC_{AB3}^{[111]}"
#endif

//  C0666666;62;"/62/""BCC_{AB2}^{[211]}"
//  C0333333;63;"/63/""BCC_{AB2}^{[211]}"
//  C0666666;64;"/64/""BCC_{AB2}^{[011]}"
//  C0333333;65;"/65/""BCC_{AB2}^{[011]}"

// ***************************************************************
// ***************************************************************
// ***************************************************************


bool PENNSY_Parameters::FixGndStatesNamesConcentrations(bool verbose) {
  int j,k;
  bool vPRE=verbose,vPOS=verbose;
  if(verbose) cerr << "***************************************************************" << endl;
  if(verbose) cerr << "FixGndStatesNamesConcentrations: start" << endl;
  for(k=1;k<=Nalloys;k++)
    for(j=1;j<=Nconcentrations2;j++) {
      GndStatesConcentrations[j][k]=-1.0;
      GndStatesNames[j][k].clear(); // empty the string
    }

  for(k=1;k<=Nalloys;k++) {
    for(j=1;j<=Nconcentrations2;j++) {
      if(GNDlibrary[j][k]) {
	GndStatesConcentrations[j][k]=concentrations[GNDlibrary[j][k]];
	//	GndStatesNames[j][k],strclean(structures_description[GNDlibrary[j][k]].c_str()));
	GndStatesNames[j][k]=structures_description[GNDlibrary[j][k]];

	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  AgAu
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Ag","Au")) {
	    PhaseDiagramNameLeft[k]="Ag";PhaseDiagramNameRight[k]="Au";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0200000)) GndStatesNames[j][k]="D1_a/tie";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="(see table)  ";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C37/MoPt_2  ";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="L1_0";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="  C37/MoPt_2";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="            L1_2/D0_{22}/D0_{23}";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  AgCd
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Ag","Cd")) {
	    PhaseDiagramNameLeft[k]="Ag";PhaseDiagramNameRight[k]="Cd";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="D0_{19}/D0_{22}/D0_{24}";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C37";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B2/B19/B27";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="ZrSi_2";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="D0_{19}";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  AgMg
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Ag","Mg")) {
	    PhaseDiagramNameLeft[k]="Ag"; PhaseDiagramNameRight[k]="Mg";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="D0_{23}/D0_{24}";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B2";
	    //if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="B8_2/Ni_2Si/tie";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="      D0_{19}/D0_{a}";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  AgNa
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Ag","Na")) {
	    PhaseDiagramNameLeft[k]="Ag"; PhaseDiagramNameRight[k]="Na";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C15";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  AgPd
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Ag","Pd")) {
	    PhaseDiagramNameLeft[k]="Pd"; PhaseDiagramNameRight[k]="Ag";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0200000)) GndStatesNames[j][k]="D1_a/tie";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="L1_1";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C37/C49";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="L1_2/D0_{22}";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  AgPt
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Ag","Pt")) {
	    PhaseDiagramNameLeft[k]="Ag"; PhaseDiagramNameRight[k]="Pt";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="troubles";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="troubles";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  AgRe
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ag","Re")) {
	    PhaseDiagramNameLeft[k]="Ag"; PhaseDiagramNameRight[k]="Re";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Ag
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  AgRh
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ag","Rh")) {
	    PhaseDiagramNameLeft[k]="Ag"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Ag
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  AgRu
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ag","Ru")) {
	    PhaseDiagramNameLeft[k]="Ag"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Ag
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  AgTi
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Ag","Ti")) {
	    PhaseDiagramNameLeft[k]="Ti"; PhaseDiagramNameRight[k]="Ag";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B11";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODEN)
	  if(strsub(alloys[k],"Ag","Ti")) {
	    PhaseDiagramNameLeft[k]="Ag"; PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B11";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ag","Ti")) {
	    PhaseDiagramNameLeft[k]="Ag"; PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B11";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  AgY
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Ag","Y")) {
	    PhaseDiagramNameLeft[k]="Y"; PhaseDiagramNameRight[k]="Ag";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C37/tie";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B2";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C11_b";
	    //if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="A15";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="  D0_a";
	    //if(isequal(GndStatesConcentrations[j][k],C0833333)) GndStatesNames[j][k]="C15_b/tie";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  AgZr
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Ag","Zr")) {
	    PhaseDiagramNameLeft[k]="Zr"; PhaseDiagramNameRight[k]="Ag";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="FCC_{AB3}^{[001]}/tie";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B11";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C6/C32";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  AlSc
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Al","Sc")) {
	    PhaseDiagramNameLeft[k]="Sc"; PhaseDiagramNameRight[k]="Al";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="D0_{19}";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="B8_2";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B2";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C15";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="L1_2";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  AuCd
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Au","Cd")) {
	    PhaseDiagramNameLeft[k]="Au"; PhaseDiagramNameRight[k]="Cd";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="D0_{24}/D0_{19}/Al_3Pu";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B19";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="L6_0";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  AuNb
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Au","Nb")) {
	    PhaseDiagramNameLeft[k]="Au"; PhaseDiagramNameRight[k]="Nb";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C32 (paw-gga)";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="BCC_{AB2}^{[011]}";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  AuPd
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Au","Pd")) {
	    PhaseDiagramNameLeft[k]="Au"; PhaseDiagramNameRight[k]="Pd";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0200000)) GndStatesNames[j][k]="D1_a/tie      ";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="D0_{23}/D0_{22}/L1_2     ";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C49/C37";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]=_STR23_CH;
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="                    D0_{23}/D0_{22}/L1_2/tie";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  AuRe
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Au","Re")) {
	    PhaseDiagramNameLeft[k]="Au"; PhaseDiagramNameRight[k]="Re";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Au
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  AuRh
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Au","Rh")) {
	    PhaseDiagramNameLeft[k]="Au"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Au
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  AuRu
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Au","Ru")) {
	    PhaseDiagramNameLeft[k]="Au"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Au
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Au","Ru")) {
	    PhaseDiagramNameLeft[k]="Au"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="HOLD FOR ACCURACY";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="HOLD FOR ACCURACY";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="HOLD FOR ACCURACY";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="HOLD FOR ACCURACY";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  AuSc
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Au","Sc")) {
	    PhaseDiagramNameLeft[k]="Au"; PhaseDiagramNameRight[k]="Sc";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0200000)) GndStatesNames[j][k]="  D1_a   ";
	    // if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="   D0_a";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="  C11_b/MoPt_2";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B2/B19";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C37  ";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  AuTi
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Au","Ti")) {
	    PhaseDiagramNameLeft[k]="Au"; PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0200000)) GndStatesNames[j][k]="  D1_a   ";
	    //if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="  D0_a";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C11_b/MoPt_2            ";
	    if(isequal(GndStatesConcentrations[j][k],C0428571)) GndStatesNames[j][k]="        Cu_4Ti_3";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="        B11";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="A15";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODEN)
	  if(strsub(alloys[k],"Au","Ti")) {
	    PhaseDiagramNameLeft[k]="Au"; PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0200000)) GndStatesNames[j][k]="  D1_a   ";
	    //if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="  D0_a";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C11_b/MoPt_2            ";
	    if(isequal(GndStatesConcentrations[j][k],C0428571)) GndStatesNames[j][k]="        Cu_4Ti_3";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="        B11";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="A15";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Au","Ti")) {
	    PhaseDiagramNameLeft[k]="Au"; PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0200000)) GndStatesNames[j][k]="  D1_a   ";
	    //if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="  D0_a";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C11_b       ";
	    if(isequal(GndStatesConcentrations[j][k],C0428571)) GndStatesNames[j][k]="         Cu_4Ti_3";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="       B11";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="A15";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  AuY
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Au","Y")) {
	    PhaseDiagramNameLeft[k]="Au"; PhaseDiagramNameRight[k]="Y";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="   D0_a";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="  B33";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C37";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  AuZr
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Au","Zr")) {
	    PhaseDiagramNameLeft[k]="Au"; PhaseDiagramNameRight[k]="Zr";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="   D0_a";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C0428571)) GndStatesNames[j][k]="          Cu_4Ti_3";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="                           B11/FCC_{A2B2}^{[001]}";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="         CuZr_2/C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="A15";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  BTi
	if(LIB_MODE==LIB_MODEN)
	  if(strsub(alloys[k],"B","Ti")) {
	    PhaseDiagramNameLeft[k]="B"; PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C32";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="        CrB/FeB-b";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  BeTi
	if(LIB_MODE==LIB_MODEN)
	  if(strsub(alloys[k],"Be","Ti")) {
	    cerr << "HERE" << endl;
	    PhaseDiagramNameLeft[k]="Be "; PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0166666)) GndStatesNames[j][k]="D2_d CaCu_5";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="        B2";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  CdPd
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Cd","Pd")) {
	    PhaseDiagramNameLeft[k]="Pd"; PhaseDiagramNameRight[k]="Cd";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="D0_{22}/NbPd_3";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C37";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="L1_0";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="D0_{19}/NbPd_3/Al_3Pu/D0_{24}";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  CdPt
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Cd","Pt")) {
	    PhaseDiagramNameLeft[k]="Cd"; PhaseDiagramNameRight[k]="Pt";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="D0_{11}/D0_{a}/D0_{22}     ";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C37/C16/tie     ";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="L1_0";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="     CdPt_3^{proto}";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  CdRe
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Cd","Re")) {
	    PhaseDiagramNameLeft[k]="Cd"; PhaseDiagramNameRight[k]="Re";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Cd
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  CdRh
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Cd","Rh")) {
	    PhaseDiagramNameLeft[k]="Cd"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Cd
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Cd","Rh")) {
	    PhaseDiagramNameLeft[k]="Rh"; PhaseDiagramNameRight[k]="Cd";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C37           ";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="D0_{24}/Al_3Pu";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Cd
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  CdRu
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Cd","Ru")) {
	    PhaseDiagramNameLeft[k]="Cd"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Cd
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  CdTi
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Cd","Ti")) {
	    PhaseDiagramNameLeft[k]="Ti";PhaseDiagramNameRight[k]="Cd";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C0428571)) GndStatesNames[j][k]="               Cu_4Ti_3/tie";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="       B11";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODEN)
	  if(strsub(alloys[k],"Cd","Ti")) {
	    PhaseDiagramNameLeft[k]="Cd";PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C11_b/CuZr_2";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="       B11";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Cd","Ti")) {
	    PhaseDiagramNameLeft[k]="Cd";PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="       B11";
	    // if(isequal(GndStatesConcentrations[j][k],C06666666)) GndStatesNames[j][k]="C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  CdY
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Cd","Y")) {
	    PhaseDiagramNameLeft[k]="Y"; PhaseDiagramNameRight[k]="Cd";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C37/C49  ";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B2";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C6";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="Cd_3Er";//D0_{19}";
	    if(isequal(GndStatesConcentrations[j][k],C0833333)) GndStatesNames[j][k]="";//YCd_5^{proto}";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  CdZr
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Cd","Zr")) {
	    PhaseDiagramNameLeft[k]="Zr"; PhaseDiagramNameRight[k]="Cd";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="A15";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B11 (paw-gga) ";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="L1_2";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  CoRe
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Co","Re")) {
	    PhaseDiagramNameLeft[k]="Co"; PhaseDiagramNameRight[k]="Re";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Co
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  CoRh
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Co","Rh")) {
	    PhaseDiagramNameLeft[k]="Co"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Co
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  CoRu
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Co","Ru")) {
	    PhaseDiagramNameLeft[k]="Co"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Co
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  CoTi
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Co","Ti")) {
	    PhaseDiagramNameLeft[k]="Co";PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="L1_2";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C14/C36";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B2";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="     C37/CuZr_2/NiTi_2";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  CrRe
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Cr","Re")) {
	    PhaseDiagramNameLeft[k]="Cr"; PhaseDiagramNameRight[k]="Re";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Cr
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  CrRh
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Cr","Rh")) {
	    PhaseDiagramNameLeft[k]="Cr"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Cr
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  CrRu
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Cr","Ru")) {
	    PhaseDiagramNameLeft[k]="Cr"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Cr
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  CrTi
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Cr","Ti")) {
	    PhaseDiagramNameLeft[k]="Cr";PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C15";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  CuRe
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Cu","Re")) {
	    PhaseDiagramNameLeft[k]="Cu"; PhaseDiagramNameRight[k]="Re";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Cu
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  CuRh
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Cu","Rh")) {
	    PhaseDiagramNameLeft[k]="Cu"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Cu
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  CuRu
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Cu","Ru")) {
	    PhaseDiagramNameLeft[k]="Cu"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Cu
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  FeRe
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Fe","Re")) {
	    PhaseDiagramNameLeft[k]="Fe"; PhaseDiagramNameRight[k]="Re";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Fe
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  FeRh
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Fe","Rh")) {
	    PhaseDiagramNameLeft[k]="Fe"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Fe
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  FeRu
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Fe","Ru")) {
	    PhaseDiagramNameLeft[k]="Fe"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Fe
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  HfRe
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Hf","Re")) {
	    PhaseDiagramNameLeft[k]="Hf"; PhaseDiagramNameRight[k]="Re";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Hf
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  HfRh
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Hf","Rh")) {
	    PhaseDiagramNameLeft[k]="Hf"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Hf
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  HfRu
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Hf","Ru")) {
	    PhaseDiagramNameLeft[k]="Hf"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Hf
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  HgRe
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Hg","Re")) {
	    PhaseDiagramNameLeft[k]="Hg"; PhaseDiagramNameRight[k]="Re";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="RHL"; // Hg
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  HgRh
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Hg","Rh")) {
	    PhaseDiagramNameLeft[k]="Hg"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="RHL"; // Hg
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  HgRu
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Hg","Ru")) {
	    PhaseDiagramNameLeft[k]="Hg"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="RHL"; // Hg
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  IrRe
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ir","Re")) {
	    PhaseDiagramNameLeft[k]="Ir"; PhaseDiagramNameRight[k]="Re";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Ir
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  IrRh
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ir","Rh")) {
	    PhaseDiagramNameLeft[k]="Ir"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Ir
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  IrRu
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ir","Ru")) {
	    PhaseDiagramNameLeft[k]="Ir"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Ir
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  IrTi
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ir","Ti")) {
	    PhaseDiagramNameLeft[k]="Ir";PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0125000)) GndStatesNames[j][k]="Ca_7Ge";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="L1_2";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="L1_0";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="  C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="  A15";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  IrRe
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ir","Re")) {
	    PhaseDiagramNameLeft[k]="Ir"; PhaseDiagramNameRight[k]="Re";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Ir
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  IrRh
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ir","Rh")) {
	    PhaseDiagramNameLeft[k]="Ir"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Ir
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  IrRu
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ir","Ru")) {
	    PhaseDiagramNameLeft[k]="Ir"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Ir
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  LaRe
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"La","Re")) {
	    PhaseDiagramNameLeft[k]="La"; PhaseDiagramNameRight[k]="Re";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // La
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  LaRh
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"La","Rh")) {
	    PhaseDiagramNameLeft[k]="La"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // La
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  LaRu
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"La","Ru")) {
	    PhaseDiagramNameLeft[k]="La"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // La
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  MgTi
	if(LIB_MODE==LIB_MODEN)
	  if(strsub(alloys[k],"Mg","Ti")) {
	    PhaseDiagramNameLeft[k]="Mg";PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  MnRe
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Mn","Re")) {
	    PhaseDiagramNameLeft[k]="Mn"; PhaseDiagramNameRight[k]="Re";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="CUB"; // Mn
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  MnRh
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Mn","Rh")) {
	    PhaseDiagramNameLeft[k]="Mn"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="CUB"; // Mn
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  MnRu
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Mn","Ru")) {
	    PhaseDiagramNameLeft[k]="Mn"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="CUB"; // Mn
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  MoNb
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Mo","Nb")) {
	    PhaseDiagramNameLeft[k]="Nb"; PhaseDiagramNameRight[k]="Mo";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B2";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="D0_{3}";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  MoPd
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Mo","Pd")) {
	    PhaseDiagramNameLeft[k]="Mo"; PhaseDiagramNameRight[k]="Pd";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="ZrSi_2";
	    if(isequal(GndStatesConcentrations[j][k],C0800000)) GndStatesNames[j][k]="  D1_a";
	    if(isequal(GndStatesConcentrations[j][k],C0833333)) GndStatesNames[j][k]="MoPd_5^{proto}";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  MoPt
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Mo","Pt")) {
	    PhaseDiagramNameLeft[k]="Mo"; PhaseDiagramNameRight[k]="Pt";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B19";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="MoPt_2";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="D0_{22}";
	    if(isequal(GndStatesConcentrations[j][k],C0800000)) GndStatesNames[j][k]="  D1_a";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  MoRe
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Mo","Re")) {
	    PhaseDiagramNameLeft[k]="Mo";PhaseDiagramNameRight[k]="Re";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Mo
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  MoRh
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Mo","Rh")) {
	    PhaseDiagramNameLeft[k]="Mo";PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Mo
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Mo","Rh")) {
	    PhaseDiagramNameLeft[k]="Mo";PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Mo
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B19";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C37           ";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="CdMg_3    ";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  MoRu
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Mo","Ru")) {
	    PhaseDiagramNameLeft[k]="Mo"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Mo
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Mo","Ru")) {
	    PhaseDiagramNameLeft[k]="Mo"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Mo
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="D0_{19}   ";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  MgTi
	if(LIB_MODE==LIB_MODEN)
	  if(strsub(alloys[k],"Mg","Ti")) {
	    PhaseDiagramNameLeft[k]="Mg"; PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  MoTi
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Mo","Ti")) {
	    PhaseDiagramNameLeft[k]="Ti"; PhaseDiagramNameRight[k]="Mo";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="MoTi_3^{proto}/tie    ";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="BCC_{AB2}^{[211]}";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="MoTi^{proto}";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="  Mo_3Ti^{proto}";
	    if(isequal(GndStatesConcentrations[j][k],C0800000)) GndStatesNames[j][k]="      D1_a/tie";
	    if(isequal(GndStatesConcentrations[j][k],C0833333)) GndStatesNames[j][k]="        Mo_5Ti^{proto}";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODEN)
	  if(strsub(alloys[k],"Mo","Ti")) {
	    PhaseDiagramNameLeft[k]="Mo"; PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="MoTi_3^{proto}/tie    ";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="BCC_{AB2}^{[211]}";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="MoTi^{proto}";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="  Mo_3Ti^{proto}";
	    if(isequal(GndStatesConcentrations[j][k],C0200000)) GndStatesNames[j][k]="      D1_a/tie";
	    if(isequal(GndStatesConcentrations[j][k],C0166666)) GndStatesNames[j][k]="        Mo_5Ti^{proto}";
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  MoZr
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Mo","Zr")) {
	    PhaseDiagramNameLeft[k]="Zr";PhaseDiagramNameRight[k]="Mo";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0166666)) GndStatesNames[j][k]="Mo_5Ti^{proto}/tie  ";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="MoZr_3^{proto}/tie  ";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="MoZr^{proto}/tie      ";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C15";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// *************************************************************** NbPd
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Nb","Pd")) {
	    PhaseDiagramNameLeft[k]="Nb";PhaseDiagramNameRight[k]="Pd";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="BCC_{AB2}^{[011]}       ";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="MoPt_2                 ";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="D0_{22}  ";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// *************************************************************** NbPt
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Nb","Pt")) {
	    PhaseDiagramNameLeft[k]="Nb";PhaseDiagramNameRight[k]="Pt";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="A15";
	    if(isequal(GndStatesConcentrations[j][k],C0600000)) GndStatesNames[j][k]="L1_0";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="MoPt_2           ";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="D0_a";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// *************************************************************** NbRe
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Nb","Re")) {
	    PhaseDiagramNameLeft[k]="Nb";PhaseDiagramNameRight[k]="Re";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Nb
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// *************************************************************** NbRh
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Nb","Rh")) {
	    PhaseDiagramNameLeft[k]="Nb";PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Nb
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Nb","Rh")) {
	    PhaseDiagramNameLeft[k]="Nb";PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Nb
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="A15";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="L1_0";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="Al_3Pu/L1_2";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh hcp ??
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// *************************************************************** NbRu
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Nb","Ru")) {
	    PhaseDiagramNameLeft[k]="Nb";PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Nb
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Nb","Ru")) {
	    PhaseDiagramNameLeft[k]="Nb";PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Nb
	    if(isequal(GndStatesConcentrations[j][k],C0166666)) GndStatesNames[j][k]="Mo_5Ti^{proto}";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="D0_3";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C37";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="D0_{24}";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  NbTc
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Nb","Tc")) {
	    PhaseDiagramNameLeft[k]="Nb"; PhaseDiagramNameRight[k]="Tc";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="Nb_3Tc^{proto}";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B2";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  NiRe
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ni","Re")) {
	    PhaseDiagramNameLeft[k]="Ni"; PhaseDiagramNameRight[k]="Re";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Ni
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  NiRh
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ni","Rh")) {
	    PhaseDiagramNameLeft[k]="Ni"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Ni
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  NiRu
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ni","Ru")) {
	    PhaseDiagramNameLeft[k]="Ni"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Ni
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  OsRe
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Os","Re")) {
	    PhaseDiagramNameLeft[k]="Os"; PhaseDiagramNameRight[k]="Re";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Os
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  OsRh
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Os","Rh")) {
	    PhaseDiagramNameLeft[k]="Os"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Os
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  OsRu
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Os","Ru")) {
	    PhaseDiagramNameLeft[k]="Os"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Os
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  PdPt
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Pd","Pt")) {
	    PhaseDiagramNameLeft[k]="Pd"; PhaseDiagramNameRight[k]="Pt";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="Pd_3Pt^{proto}/tie    ";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="L1_1";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="PdPt_3^{proto}";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  PdRe
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Pd","Re")) {
	    PhaseDiagramNameLeft[k]="Pd"; PhaseDiagramNameRight[k]="Re";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Pd
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  PdRh
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Pd","Rh")) {
	    PhaseDiagramNameLeft[k]="Pd"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Pd
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  PdRu
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Pd","Ru")) {
	    PhaseDiagramNameLeft[k]="Pd"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Pd
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  PdY
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Pd","Y")) {
	    PhaseDiagramNameLeft[k]="Y"; PhaseDiagramNameRight[k]="Pd";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="D0_{11}  ";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C37  ";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B27/B33";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="L1_2";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  PdTc
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Pd","Tc")) {
	    PhaseDiagramNameLeft[k]="Tc"; PhaseDiagramNameRight[k]="Pd";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="D0_{19}";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// *************************************************************** PdTi
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Pd","Ti")) {
	    PhaseDiagramNameLeft[k]="Ti";PhaseDiagramNameRight[k]="Pd";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="A15";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B19/L1_0    ";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="MoPt_2                 ";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="D0_{24}";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODEN)
	  if(strsub(alloys[k],"Pd","Ti")) {
	    PhaseDiagramNameLeft[k]="Pd";PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="A15";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B19/L1_0";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="       MoPt_2";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="D0_{24}";
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  PdZr
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Pd","Zr")) {
	    PhaseDiagramNameLeft[k]="Zr";PhaseDiagramNameRight[k]="Pd";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]=_STR15_Z3;
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B27/B33       ";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C11_b/tie (paw-gga)                                 ";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="D0_{24}";
	    if(isequal(GndStatesConcentrations[j][k],C0833333)) GndStatesNames[j][k]="C2/m #12";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  PtRe
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Pt","Re")) {
	    PhaseDiagramNameLeft[k]="Pt"; PhaseDiagramNameRight[k]="Re";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Pt
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  PtRh
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Pt","Rh")) {
	    PhaseDiagramNameLeft[k]="Pt"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Pt
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Pt","Rh")) {
	    PhaseDiagramNameLeft[k]="Pt"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Pt
	    if(isequal(GndStatesConcentrations[j][k],C0200000)) GndStatesNames[j][k]="D1_a";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="D0_{22}";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]=_STR23_CH;
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C49";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="D0_{22}";
	    if(isequal(GndStatesConcentrations[j][k],C0800000)) GndStatesNames[j][k]="D1_a";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  PtRu
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Pt","Ru")) {
	    PhaseDiagramNameLeft[k]="Pt"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Pt
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Pt","Ru")) {
	    PhaseDiagramNameLeft[k]="Pt"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Pt
	    // if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]=_STR13_Z1;
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="FCC_{AB3}^{[001]}/tie   ";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]=_STR14_Z2;
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";// Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  PtTc
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Pt","Tc")) {
	    PhaseDiagramNameLeft[k]="Pt"; PhaseDiagramNameRight[k]="Tc";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]=_STR13_Z1;
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="D0_{19}";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// *************************************************************** PtTi
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Pt","Ti")) {
	    PhaseDiagramNameLeft[k]="Pt";PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="D0_{24}";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C49/C37";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B19";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="A15";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODEN)
	  if(strsub(alloys[k],"Pt","Ti")) {
	    PhaseDiagramNameLeft[k]="Pt";PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="D0_{24}";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C49/C37";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B19";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="A15";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  PtY
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Pt","Y")) {
	    PhaseDiagramNameLeft[k]="Pt"; PhaseDiagramNameRight[k]="Y";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="L1_2";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C15";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B33";
	    if(isequal(GndStatesConcentrations[j][k],C0625000)) GndStatesNames[j][k]=" D8_8";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C37";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="D0_{11}";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  PtZr
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Pt","Zr")) {
	    PhaseDiagramNameLeft[k]="Pt"; PhaseDiagramNameRight[k]="Zr";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="D0_{24}";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B33";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C16";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="A15";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  ReRh
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Re","Rh")) {
	    PhaseDiagramNameLeft[k]="Re"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  ReRu
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Re","Ru")) {
	    PhaseDiagramNameLeft[k]="Re"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  ReSc
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Re","Sc")) {
	    PhaseDiagramNameLeft[k]="Re"; PhaseDiagramNameRight[k]="Sc";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Sc
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  ReTa
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Re","Ta")) {
	    PhaseDiagramNameLeft[k]="Re"; PhaseDiagramNameRight[k]="Ta";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC"; // Ta
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  ReTc
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Re","Tc")) {
	    PhaseDiagramNameLeft[k]="Re"; PhaseDiagramNameRight[k]="Tc";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Tc
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  ReTi
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Re","Ti")) {
	    PhaseDiagramNameLeft[k]="Re"; PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ti
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  ReV
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Re","V")) {
	    PhaseDiagramNameLeft[k]="Re"; PhaseDiagramNameRight[k]="V";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC"; // V
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  ReW
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Re","W")) {
	    PhaseDiagramNameLeft[k]="Re"; PhaseDiagramNameRight[k]="W";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC"; // W
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  ReY
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Re","Y")) {
	    PhaseDiagramNameLeft[k]="Re"; PhaseDiagramNameRight[k]="Y";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Y
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  ReZn
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Re","Zn")) {
	    PhaseDiagramNameLeft[k]="Re"; PhaseDiagramNameRight[k]="Zn";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Zn
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  ReZr
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Re","Zr")) {
	    PhaseDiagramNameLeft[k]="Re"; PhaseDiagramNameRight[k]="Zr";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Re
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Tc
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  RhRu
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Rh","Ru")) {
	    PhaseDiagramNameLeft[k]="Rh"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Rh","Ru")) {
	    PhaseDiagramNameLeft[k]="Ru"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="127";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="RhRu_2^{proto}";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="RhRu^{proto}";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  RhSc
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Rh","Sc")) {
	    PhaseDiagramNameLeft[k]="Rh"; PhaseDiagramNameRight[k]="Sc";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Sc
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  RhTa
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Rh","Ta")) {
	    PhaseDiagramNameLeft[k]="Rh"; PhaseDiagramNameRight[k]="Ta";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC"; // Ta
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  RhTc
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Rh","Tc")) {
	    PhaseDiagramNameLeft[k]="Rh"; PhaseDiagramNameRight[k]="Tc";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Tc
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Rh","Tc")) {
	    PhaseDiagramNameLeft[k]="Tc"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Tc
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="D0_{19}";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B19";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="ZrSi_2";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  RhTi
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Rh","Ti")) {
	    PhaseDiagramNameLeft[k]="Rh"; PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ti
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Rh","Ti")) {
	    PhaseDiagramNameLeft[k]="Ti"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C11_b/CuZr_2   ";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="L1_0";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C37";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="L1_2";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODEN)
	  if(strsub(alloys[k],"Rh","Ti")) {
	    PhaseDiagramNameLeft[k]="Rh"; PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C11_b/CuZr_2   ";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="L1_0";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C37";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="L1_2";
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  RhV
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Rh","V")) {
	    PhaseDiagramNameLeft[k]="Rh"; PhaseDiagramNameRight[k]="V";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC"; // V
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  RhW
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Rh","W")) {
	    PhaseDiagramNameLeft[k]="Rh"; PhaseDiagramNameRight[k]="W";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC"; // W
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  RhY
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Rh","Y")) {
	    PhaseDiagramNameLeft[k]="Rh"; PhaseDiagramNameRight[k]="Y";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Y
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Rh","Y")) {
	    PhaseDiagramNameLeft[k]="Y"; PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Y
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="D0_{11}";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C37";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B2";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C15";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  RhZn
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Rh","Zn")) {
	    PhaseDiagramNameLeft[k]="Rh"; PhaseDiagramNameRight[k]="Zn";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Zn
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  RhZr
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Rh","Zr")) {
	    PhaseDiagramNameLeft[k]="Rh"; PhaseDiagramNameRight[k]="Zr";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Tc
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Rh","Zr")) {
	    PhaseDiagramNameLeft[k]="Zr";PhaseDiagramNameRight[k]="Rh";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Tc
	    if(isequal(GndStatesConcentrations[j][k],C0200000)) GndStatesNames[j][k]="D1_{a}/tie    ";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="FCC_{AB3}^{[001]}/tie  ";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C16";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B27";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C37/tie   ";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="L1_2";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="FCC"; // Rh
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// *************************************************************** RuSc
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ru","Sc")) {
	    PhaseDiagramNameLeft[k]="Ru";PhaseDiagramNameRight[k]="Sc";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Sc
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// *************************************************************** RuTa
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ru","Ta")) {
	    PhaseDiagramNameLeft[k]="Ru";PhaseDiagramNameRight[k]="Ta";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Ta
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// *************************************************************** RuTc
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ru","Tc")) {
	    PhaseDiagramNameLeft[k]="Ru";PhaseDiagramNameRight[k]="Tc";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Tc
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Ru","Tc")) {
	    PhaseDiagramNameLeft[k]="Tc"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Tc
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="D0_{19}";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B19";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="D0_{19}";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// *************************************************************** RuTi
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ru","Ti")) {
	    PhaseDiagramNameLeft[k]="Ru";PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="RuTi_3^{proto}";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C49";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B2";
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Ti
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Ru","Ti")) {
	    PhaseDiagramNameLeft[k]="Ti";PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Ti
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="RuTi_3^{proto}";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C49";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B2";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODEN)
	  if(strsub(alloys[k],"Ru","Ti")) {
	    PhaseDiagramNameLeft[k]="Ru";PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="RuTi_3^{proto}";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C49";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B2";
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Ti
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  RuV
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ru","V")) {
	    PhaseDiagramNameLeft[k]="Ru"; PhaseDiagramNameRight[k]="V";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC"; // V
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  RuW
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ru","W")) {
	    PhaseDiagramNameLeft[k]="Ru"; PhaseDiagramNameRight[k]="W";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC"; // W
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  RuY
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ru","Y")) {
	    PhaseDiagramNameLeft[k]="Ru"; PhaseDiagramNameRight[k]="Y";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Y
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Ru","Y")) {
	    PhaseDiagramNameLeft[k]="Y"; PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Y
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="D0_{11}";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C16";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C14";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  RuZn
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ru","Zn")) {
	    PhaseDiagramNameLeft[k]="Ru"; PhaseDiagramNameRight[k]="Zn";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Zn
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  RuZr
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ru","Zr")) {
	    PhaseDiagramNameLeft[k]="Ru"; PhaseDiagramNameRight[k]="Zr";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Zr
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Ru","Zr")) {
	    PhaseDiagramNameLeft[k]="Zr";PhaseDiagramNameRight[k]="Ru";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Zr
	    if(isequal(GndStatesConcentrations[j][k],C0200000)) GndStatesNames[j][k]="D1_{a}   ";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B2";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ru
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  ScTa
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Sc","Ta")) {
	    PhaseDiagramNameLeft[k]="Sc"; PhaseDiagramNameRight[k]="Ta";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Sc
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC"; // Ta
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  ScTc
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Sc","Tc")) {
	    PhaseDiagramNameLeft[k]="Sc"; PhaseDiagramNameRight[k]="Tc";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Sc
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Tc
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  ScTi
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Sc","Ti")) {
	    PhaseDiagramNameLeft[k]="Sc"; PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Sc
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ti
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  ScV
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Sc","V")) {
	    PhaseDiagramNameLeft[k]="Sc"; PhaseDiagramNameRight[k]="W";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Sc
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC"; // V
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  ScW
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Sc","W")) {
	    PhaseDiagramNameLeft[k]="Sc"; PhaseDiagramNameRight[k]="W";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Sc
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC"; // W
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  ScY
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Sc","Y")) {
	    PhaseDiagramNameLeft[k]="Sc"; PhaseDiagramNameRight[k]="Y";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Sc
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Y
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  ScZn
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Sc","Zn")) {
	    PhaseDiagramNameLeft[k]="Sc"; PhaseDiagramNameRight[k]="Zn";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Sc
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Zn
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  ScZr
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Sc","Zr")) {
	    PhaseDiagramNameLeft[k]="Sc"; PhaseDiagramNameRight[k]="Zr";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Sc
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Zr
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  TaTc
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ta","Tc")) {
	    PhaseDiagramNameLeft[k]="Ta"; PhaseDiagramNameRight[k]="Tc";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Ta
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Tc
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  TaTi
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ta","Ti")) {
	    PhaseDiagramNameLeft[k]="Ta"; PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Ta
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ti
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  TaV
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ta","V")) {
	    PhaseDiagramNameLeft[k]="Ta"; PhaseDiagramNameRight[k]="W";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Ta
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC"; // V
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  TaW
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ta","W")) {
	    PhaseDiagramNameLeft[k]="Ta"; PhaseDiagramNameRight[k]="W";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Ta
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC"; // W
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  TaY
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ta","Y")) {
	    PhaseDiagramNameLeft[k]="Ta"; PhaseDiagramNameRight[k]="Y";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Ta
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Y
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  TaZn
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ta","Zn")) {
	    PhaseDiagramNameLeft[k]="Ta"; PhaseDiagramNameRight[k]="Zn";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Ta
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Zn
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  TaZr
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ta","Zr")) {
	    PhaseDiagramNameLeft[k]="Ta"; PhaseDiagramNameRight[k]="Zr";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // Ta
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Zr
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  TcTi
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Tc","Ti")) {
	    PhaseDiagramNameLeft[k]="Tc"; PhaseDiagramNameRight[k]="Ti";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Tc
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Ti
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Tc","Ti")) {
	    PhaseDiagramNameLeft[k]="Ti"; PhaseDiagramNameRight[k]="Tc";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Ti
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="TcTi_3^{proto}";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C49";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B2";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C11_b";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Tc
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  TcV
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Tc","V")) {
	    PhaseDiagramNameLeft[k]="Tc"; PhaseDiagramNameRight[k]="V";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Tc
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC"; // V
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  TcW
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Tc","W")) {
	    PhaseDiagramNameLeft[k]="Tc"; PhaseDiagramNameRight[k]="W";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Tc
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC"; // W
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  TcY
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Tc","Y")) {
	    PhaseDiagramNameLeft[k]="Tc"; PhaseDiagramNameRight[k]="Y";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Tc
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Y
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Tc","Y")) {
	    PhaseDiagramNameLeft[k]="Y"; PhaseDiagramNameRight[k]="Tc";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Y
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="D0_{11}";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C14";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Tc
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  TcZn
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Tc","Zn")) {
	    PhaseDiagramNameLeft[k]="Tc"; PhaseDiagramNameRight[k]="Zn";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Tc
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Zn
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  TcZr
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Tc","Zr")) {
	    PhaseDiagramNameLeft[k]="Tc"; PhaseDiagramNameRight[k]="Zr";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Tc
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Zr
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Tc","Zr")) {
	    PhaseDiagramNameLeft[k]="Zr"; PhaseDiagramNameRight[k]="Tc";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Zr
	    if(isequal(GndStatesConcentrations[j][k],C0200000)) GndStatesNames[j][k]="D1_{a}   ";
	    //if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="TcZr_3^{proto}/tie    ";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="C49";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B2";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C14";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Tc
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  TiV
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ti","V")) {
	    PhaseDiagramNameLeft[k]="Ti"; PhaseDiagramNameRight[k]="W";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Ti
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC"; // V
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  TiW
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ti","W")) {
	    PhaseDiagramNameLeft[k]="Ti"; PhaseDiagramNameRight[k]="W";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Ti
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="BCC_{AB2}^{[211]}";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]=" ";
	    if(isequal(GndStatesConcentrations[j][k],C0800000)) GndStatesNames[j][k]="D1_a";
	    if(isequal(GndStatesConcentrations[j][k],C0833333)) GndStatesNames[j][k]=" ";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC"; // W
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  TiY
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ti","Y")) {
	    PhaseDiagramNameLeft[k]="Ti"; PhaseDiagramNameRight[k]="Y";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Ti
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Y
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  TiZn
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ti","Zn")) {
	    PhaseDiagramNameLeft[k]="Ti"; PhaseDiagramNameRight[k]="Zn";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Ti
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="A15";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="L1_0/B2";
	    if(isequal(GndStatesConcentrations[j][k],C0666666)) GndStatesNames[j][k]="C14";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="L1_2";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Zn
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  TiZr
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Ti","Zr")) {
	    PhaseDiagramNameLeft[k]="Ti"; PhaseDiagramNameRight[k]="Zr";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Ti
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Zr
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Ti","Zr")) {
	    PhaseDiagramNameLeft[k]="Ti"; PhaseDiagramNameRight[k]="Zr";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Ti
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Zr
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  VW
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"V","W")) {
	    PhaseDiagramNameLeft[k]="V"; PhaseDiagramNameRight[k]="W";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // V
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="BCC"; // W
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  VY
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"V","Y")) {
	    PhaseDiagramNameLeft[k]="V"; PhaseDiagramNameRight[k]="Y";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // V
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Y
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  VZn
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"V","Zn")) {
	    PhaseDiagramNameLeft[k]="V"; PhaseDiagramNameRight[k]="Zn";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // V
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Zr
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  VZr
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"V","Zr")) {
	    PhaseDiagramNameLeft[k]="V"; PhaseDiagramNameRight[k]="Zr";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // V
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Zr
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  WY
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"W","Y")) {
	    PhaseDiagramNameLeft[k]="W"; PhaseDiagramNameRight[k]="Y";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // W
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Y
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  WZn
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"W","Zn")) {
	    PhaseDiagramNameLeft[k]="W"; PhaseDiagramNameRight[k]="Zn";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // W
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Zr
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  WZr
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"W","Zr")) {
	    PhaseDiagramNameLeft[k]="W"; PhaseDiagramNameRight[k]="Zr";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="BCC"; // W
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Zr
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  YZn
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Y","Zn")) {
	    PhaseDiagramNameLeft[k]="Y"; PhaseDiagramNameRight[k]="Zn";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Y
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Zr
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  YZr
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Y","Zr")) {
	    PhaseDiagramNameLeft[k]="Y"; PhaseDiagramNameRight[k]="Zr";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Y
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Zr
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ******************************************************************************************************************************
	// ***************************************************************  ZnZr
	if(LIB_MODE==LIB_MODEU)
	  if(strsub(alloys[k],"Zn","Zr")) {
	    PhaseDiagramNameLeft[k]="Zn"; PhaseDiagramNameRight[k]="Zr";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP"; // Zn
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="HCP"; // Zr
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ******************************************************************************************************************************
	// *                                                                                                                            *
	// *                                                  Mg-Calphad Project                                                        *
	// *                                                                                                                            *
	// ******************************************************************************************************************************
	// ***************************************************************  Mg-Calphad Project
	// ***************************************************************  MgBi
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Mg","Bi")) {  // FIXED
	    PhaseDiagramNameLeft[k]="Mg"; PhaseDiagramNameRight[k]="Bi";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0400000)) GndStatesNames[j][k]="D5_2";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="A7";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  MgIn
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Mg","In")) { // FIXED
	    PhaseDiagramNameLeft[k]="Mg"; PhaseDiagramNameRight[k]="In";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0250000)) GndStatesNames[j][k]="L1_2";
	    if(isequal(GndStatesConcentrations[j][k],C0333333)) GndStatesNames[j][k]="Mg_2In";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="L1_0";
	    if(isequal(GndStatesConcentrations[j][k],C0750000)) GndStatesNames[j][k]="L1_2";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="A6 In";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  MgSb
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Mg","Sb")) { // FIXED
	    PhaseDiagramNameLeft[k]="Mg"; PhaseDiagramNameRight[k]="Sb";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="HCP";
	    if(isequal(GndStatesConcentrations[j][k],C0400000)) GndStatesNames[j][k]="D5_2";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="A7";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  InBi
	if(strsub(alloys[k],"In","Bi")) { // FIXED
	  PhaseDiagramNameLeft[k]="In"; PhaseDiagramNameRight[k]="Bi";
	  if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="A6";
	  if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B10";
	  if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="A7";
	  if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	}
	// ***************************************************************  InSb
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"In","Sb")) { // FIXED
	    PhaseDiagramNameLeft[k]="In"; PhaseDiagramNameRight[k]="Sb";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="A6 In";
	    if(isequal(GndStatesConcentrations[j][k],C0500000)) GndStatesNames[j][k]="B3";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="A7";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************  BiSb
	if(LIB_MODE==LIB_MODE2)
	  if(strsub(alloys[k],"Bi","Sb")) { // FIXED
	    PhaseDiagramNameLeft[k]="Bi"; PhaseDiagramNameRight[k]="Sb";
	    if(vPRE) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	    if(isequal(GndStatesConcentrations[j][k],C0000000)) GndStatesNames[j][k]="A7";
	    if(isequal(GndStatesConcentrations[j][k],C1000000)) GndStatesNames[j][k]="A7";
	    if(vPOS) cerr << (alloys[k]) << " " << GndStatesConcentrations[j][k] << " " << GndStatesNames[j][k] << endl;
	  }
	// ***************************************************************
      }
    }
  }
  if(verbose) cerr << "FixGndStatesNamesConcentrations: stop" << endl;
  if(verbose) cerr << "***************************************************************" << endl;
  return TRUE;
}

string PENNSY_Parameters::MatlabGndStatesNamesConcentrations(bool _verbose) {
  stringstream oss;
  int i,j,k;
  double x,y,z,minE,maxE,deltaE;
  bool xbool;
  if(_verbose) cerr << "***************************************************************" << endl;
  if(_verbose) cerr << "MatlabGndStatesNamesConcentrations: start" << endl;
  oss << "disp('**************************************************************************');   " << endl;
  oss << "disp('*                                                                        *');   " << endl;
  oss << "disp('*    GNDstate analyzer -  Stefano Curtarolo                              *');   " << endl;
  oss << "disp('*                                                                        *');   " << endl;
  oss << "disp('**************************************************************************');   " << endl;
  for(k=1;k<=Nalloys;k++) {
    oss << "%disp('**************************************************************************');   " << "  % "<<alloys[k] << endl;
    oss << "%disp('*                                                                        *');   " << "  % "<<alloys[k] << endl;
    oss << "%disp('*    GNDstate analyzer -  Stefano Curtarolo                              *');   " << "  % "<<alloys[k] << endl;
    oss << "%disp('*                                                                        *');   " << "  % "<<alloys[k] << endl;
    oss << "%disp('**************************************************************************');   " << "  % "<<alloys[k] << endl;
    oss << "%disp(' Generated = " << TODAY << "');   " << "  % "<<alloys[k]<<endl;
    oss << "disp('" << alloys[k] << "')" << "  % "<<alloys[k]<<endl;
    oss << "                                                                                      " << "  % "<<alloys[k]<<endl;
    oss << "clear;                                                                                " << "  % "<<alloys[k]<<endl;
    oss << "FONTDIM=12;                                                                           " << "  % "<<alloys[k]<<endl;
    oss << "COLOR=[0 0 0];                                                                        " << "  % "<<alloys[k]<<endl;
    oss << "LINEWIDTH=1;                                                                          " << "  % "<<alloys[k]<<endl;
    oss << "PAPERPOSITION=[1 1 5 4];                                                              " << "  % "<<alloys[k]<<endl;
    oss << "PAPERSIZE=[8.5 11];                                                                   " << "  % "<<alloys[k]<<endl;
    oss << "POSITION=[420 540 560 420];                                                           " << "  % "<<alloys[k]<<endl;
    oss << "AXISWIDTH=1;                                                                          " << "  % "<<alloys[k]<<endl;
    oss << "AXISPOSITION=[0.136646 0.164179 0.774845 0.78678];                                    " << "  % "<<alloys[k]<<endl;
    oss << "clf;                                                                                  " << "  % "<<alloys[k]<<endl;
    oss << "% ****************************************************************************** " << "  % "<<alloys[k]<<endl;
    oss << "% ****************************************************************************** " << "  % "<<alloys[k]<<endl;
    oss << "% CONVEX_HULL_MATLAB CODE START                                          " << "  % "<<alloys[k]<<endl;
    oss << "% " << alloys[k] << "  % "<<alloys[k]<<endl;
    oss << "cd /home/auro/work/AFLOW3" << "  % "<<alloys[k]<<endl;
    oss << "H_FIGURE=figure;" << endl;
    oss << "set(H_FIGURE,'PaperUnits','inches');" << endl;
    oss << "set(H_FIGURE,'PaperOrientation','portrait');" << endl;
    oss << "set(H_FIGURE,'PaperPosition',PAPERPOSITION); % left top width height" << endl;
    oss << "set(H_FIGURE,'PaperSize',PAPERSIZE);" << endl;
    oss << "set(H_FIGURE,'PaperType','usletter');" << endl;
    oss << "set(H_FIGURE,'Position',POSITION);" << endl;
    oss << "Cb=[";
    for(i=1;i<=Nconcentrations2;i++)
      if(GNDlibrary[i][k]>0)
	oss << "  " << 100*C[i];
    oss << "];" << endl;
    oss << "EF=[";
    for(i=1;i<=Nconcentrations2;i++)
      if(GNDlibrary[i][k]>0)
	if(Library[GNDlibrary[i][k]][k].energy_formation_per_atom>=0)
	  oss << "  " << Library[GNDlibrary[i][k]][k].energy_formation_per_atom;
	else
	  oss << " " << Library[GNDlibrary[i][k]][k].energy_formation_per_atom;
    oss << "];" << endl;
    oss << "L=[";
    for(i=1;i<=Nconcentrations2;i++)
      if(GNDlibrary[i][k]>0)
	oss << "    " << CleanNameStructure(structures_name.at(GNDlibrary[i][k]));
    oss << "];" << endl;
    minE=0.0;maxE=0.0;
    for(i=1;i<=Nconcentrations2;i++) {
      if(GNDlibrary[i][k]>0)
	if(minE>Library[GNDlibrary[i][k]][k].energy_formation_per_atom)
	  minE=Library[GNDlibrary[i][k]][k].energy_formation_per_atom;
    }
    maxE=-minE/10+0.01;
    deltaE=maxE-minE;minE-=deltaE/20;minE-=deltaE/20;minE-=deltaE/20;deltaE=maxE-minE;
    //oss << "plot(Cb ,EF ,'b+-','LineWidth',0.5*LINEWIDTH);hold on;axis([0 100 %f %f]);\n",minE,maxE);
    oss << "plot(Cb ,EF ,'LineWidth',0.5*LINEWIDTH);hold on;axis([0 100 " << minE << " " << maxE << "])" << endl;
    for(i=1;i<=Nstructures;i++) {
      xbool=TRUE; for(j=1;j<=Nconcentrations2;j++) if(i==GNDlibrary[j][k]) xbool=FALSE;
      if(xbool && Library[i][k].gnd_available)
	oss << "plot(" << 100*concentrations[i] << "," << Library[i][k].energy_formation_per_atom << ",'rx','LineWidth',0.2*LINEWIDTH);hold on;" << endl;
    }

    for(i=1;i<=Nconcentrations2;i++)
      if(GNDlibrary[i][k]>0) {
	z=x=C[i];y=Library[GNDlibrary[i][k]][k].energy_formation_per_atom;
	oss << "plot(" << 100*x << "," << y << ",'b+','LineWidth',LINEWIDTH);hold on;";
	if(x<0.5) x-=0.05;if(x>0.5) x+=0.03;if(x<0.1) x-=0.03;
	oss << "text(" << 100*x-0.5*strlen(GndStatesNames[i][k].c_str()) << "," << y-deltaE/15 << ",'" << GndStatesNames[i][k] << "','fontsize',0.75*FONTDIM);" << endl;
      }
    oss << "text(" <<  -5.0 << "," << minE-deltaE/10 << ",'" << PhaseDiagramNameLeft[k] << "','fontsize',FONTDIM);" << endl;
    oss << "text(" << 100.0 << "," << minE-deltaE/10 << ",'" << PhaseDiagramNameRight[k] << "','fontsize',FONTDIM);" << endl;
    oss << "drawnow;" << "  % "<<alloys[k]<<endl;
    oss << "xlabel('Atomic Percent " << ElementSymbolToString(PhaseDiagramNameRight[k]) << "','fontsize',FONTDIM);" << endl;
    oss << "ylabel('eV/atom','fontsize',FONTDIM);" << endl;
    oss << "H_AXIS=gca;" << endl;
    oss << "set(H_AXIS,'LineWidth',AXISWIDTH);" << endl;
    oss << "set(H_AXIS,'Position',AXISPOSITION);" << endl;
    oss << "set(H_AXIS,'FontSize',FONTDIM);" << endl;
    oss << "print -depsc " << speciesAB[k] << ".eps " << endl;
    oss << "% CONVEX_HULL_MATLAB CODE END                                         " << "  % "<<alloys[k]<<endl;
    oss << "% *********************************************************************************** " << "  % "<<alloys[k]<<endl;
  }
  if(_verbose) cerr << "MatlabGndStatesNamesConcentrations: stop" << endl;
  if(_verbose) cerr << "***************************************************************" << endl;
  oss << "exit" << endl;
  return oss.str();
}

// ***************************************************************
// ***************************************************************
// ***************************************************************

bool strsub(string strzero, string subA) {
  return (bool) (strstr(strzero.c_str(),subA.c_str()));
}

bool strsub(string strzero, string subA, string subB) {
  return (bool) (strstr(strzero.c_str(),subA.c_str()) && strstr(strzero.c_str(),subB.c_str()));
}

bool strsub(string strzero, string subA, string subB, string subC) {
  return (bool) (strstr(strzero.c_str(),subA.c_str()) && strstr(strzero.c_str(),subB.c_str()) && strstr(strzero.c_str(),subC.c_str()));
}

string ElementSymbolToString(string stringin) {
  // need to clean up the _pv stuff of VASP
#define ELEMENT_STRING_MAXLEN   32
  string output;
  output="NotFound";
  if(stringin=="Ag") return "Silver";
  if(stringin=="Al") return "Aluminum";
  if(stringin=="Ar") return "Argon";
  if(stringin=="As") return "Arsenic";
  if(stringin=="Au") return "Gold";
  if(stringin=="B")  return "Boron";
  if(stringin=="Bi") return "Bismuth";
  if(stringin=="Be") return "Beryllium";
  if(stringin=="Br") return "Bromine";
  if(stringin=="C")  return "Carbon";
  if(stringin=="Ca") return "Calcium";
  if(stringin=="Cd") return "Cadmium";
  if(stringin=="Cl") return "Clorine";
  if(stringin=="Co") return "Cobalt";
  if(stringin=="Cr") return "Chromium";
  if(stringin=="Cu") return "Copper";
  if(stringin=="F")  return "Fluorine";
  if(stringin=="Fe") return "Iron";
  if(stringin=="Ga") return "Gallium";
  if(stringin=="Ge") return "Germanium";
  if(stringin=="H")  return "Hydrogen";
  if(stringin=="He") return "Helium";
  if(stringin=="Hf") return "Hafnium";
  if(stringin=="Hg") return "Mercury";
  if(stringin=="In") return "Indium";
  if(stringin=="K")  return "Potassium";
  if(stringin=="Kr") return "Krypton";
  if(stringin=="Li") return "Litium";
  if(stringin=="Mg") return "Magnesium";
  if(stringin=="Mn") return "Manganese";
  if(stringin=="Mo") return "Molybdenum";
  if(stringin=="N")  return "Nitrogen";
  if(stringin=="Na") return "Sodium";
  if(stringin=="Nb") return "Niobium";
  if(stringin=="Ne") return "Neon";
  if(stringin=="Ni") return "Nickel";
  if(stringin=="O")  return "Oxigen";
  if(stringin=="P")  return "Phosphorus";
  if(stringin=="Pd") return "Palladium";
  if(stringin=="Po") return "Polonium";
  if(stringin=="Pt") return "Platinum";
  if(stringin=="Pu") return "Plutonium";
  if(stringin=="S")  return "Sulfur";
  if(stringin=="Sb") return "Antimony";
  if(stringin=="Sc") return "Scandium";
  if(stringin=="Se") return "Selenium";
  if(stringin=="Si") return "Silicon";
  if(stringin=="Re") return "Rhenium";
  if(stringin=="Rh") return "Rhodium";
  if(stringin=="Ru") return "Ruthenium";
  if(stringin=="Ta") return "Tantalum";
  if(stringin=="Tc") return "Technetium";
  if(stringin=="Ti") return "Titanium";
  if(stringin=="U")  return "Uranium";
  if(stringin=="V")  return "Vanadium";
  if(stringin=="W")  return "Tungstem";
  if(stringin=="Xe") return "Xenon";
  if(stringin=="Y")  return "Yttrium";
  if(stringin=="Zn") return "Zinc";
  if(stringin=="Zr") return "Zirconium";
  return output;
}

// **************************************************************************
// *                                                                        *
// *  PHASES STEFANO CURTAROLO MASSACHUSETTS INSTITUTE OF TECHNOLOGY 2003   *
// *                                                                        *
// **************************************************************************

