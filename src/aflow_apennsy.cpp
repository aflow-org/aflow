// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

#ifndef _AFLOW_APENNSY_CPP_
#define _AFLOW_APENNSY_CPP_
#include "aflow.h"
//    vector<string> vrefs;
//    SystemReferences(alloys[k],vrefs);
//    for(uint i=0;i<vrefs.size();i++) oss << strfile1 << "REFERENCE: " << vrefs.at(i) << endl;

// ---------------------------------------------------------------------------

// ****************************************************************************************************
// ****************************************************************************************************
// CLASS APENNSY_Parameters

// constructors
APENNSY_Parameters::APENNSY_Parameters() {
  PseudopotentialNoclean=FALSE;
  LIBstructures_calcs=0;
  LIBstructures_holes=0;
  LIB_MODE=0;
  _flag_ce_cswap_=FALSE;
  _flag_load_xstructures=FALSE;
  Load_ALL=TRUE;
  Load_FCC=FALSE;
  Load_BCC=FALSE;
  Load_HCP=FALSE;
  alloysmesg.clear();
  alloys.clear();
  alloysRAW.clear();
  speciesA.clear();speciesB.clear();speciesAB.clear();
  pseudosA.clear();pseudosB.clear();
  Alloy2MiscibilityHT.clear();
  Alloy2MiscibilityExperiments.clear();
  Alloy2MiscibilityMiedema.clear();
  Alloy2MiscibilityHumeRothery.clear();
  PhaseDiagramNameLeft.clear();
  PhaseDiagramNameRight.clear();
  ZGNDlibrary.clear();
  ZMINElibrary.clear();
  ZConcentrations.clear();
  AlloyStructureIdentity.clear();
  RankLib.clear();
  RankLibInt.clear();
  ZLibrary.clear();
  // complex systems
  enthalpy_formation_atom_max.clear();
  entropic_temperature_max.clear();
  alloy_order.clear();
  alloy_holes.clear();
  alloy_deltamuBArich.clear();
  alloy_deltamuABrich.clear();
}

// destructor
APENNSY_Parameters::~APENNSY_Parameters() {
  free();
}

void APENNSY_Parameters::free() {
}

void APENNSY_Parameters::copy(const APENNSY_Parameters& b) {
  PseudopotentialNoclean=b.PseudopotentialNoclean;
  LIBstructures_calcs=b.LIBstructures_calcs;
  LIBstructures_holes=b.LIBstructures_holes;
  LIB_MODE=b.LIB_MODE;
  _flag_ce_cswap_=b._flag_ce_cswap_;
  _flag_load_xstructures=b._flag_load_xstructures;
  Load_ALL=b.Load_ALL;
  Load_FCC=b.Load_FCC;
  Load_BCC=b.Load_BCC;
  Load_HCP=b.Load_HCP;
  alloysmesg << b.alloysmesg.str();
  alloys.clear(); for(uint i=0;i<b.alloys.size();i++) alloys.push_back(b.alloys.at(i));
  alloysRAW.clear(); for(uint i=0;i<b.alloysRAW.size();i++) alloysRAW.push_back(b.alloysRAW.at(i));
  speciesA.clear(); for(uint i=0;i<b.speciesA.size();i++) speciesA.push_back(b.speciesA.at(i));
  speciesB.clear(); for(uint i=0;i<b.speciesB.size();i++) speciesB.push_back(b.speciesB.at(i));
  speciesAB.clear(); for(uint i=0;i<b.speciesAB.size();i++) speciesAB.push_back(b.speciesAB.at(i));
  pseudosA.clear(); for(uint i=0;i<b.pseudosA.size();i++) pseudosA.push_back(b.pseudosA.at(i));
  pseudosB.clear(); for(uint i=0;i<b.pseudosB.size();i++) pseudosB.push_back(b.pseudosB.at(i));
  Alloy2MiscibilityHT.clear(); for(uint i=0;i<b.Alloy2MiscibilityHT.size();i++) Alloy2MiscibilityHT.push_back(b.Alloy2MiscibilityHT.at(i));
  Alloy2MiscibilityExperiments.clear(); for(uint i=0;i<b.Alloy2MiscibilityExperiments.size();i++) Alloy2MiscibilityExperiments.push_back(b.Alloy2MiscibilityExperiments.at(i));
  Alloy2MiscibilityMiedema.clear(); for(uint i=0;i<b.Alloy2MiscibilityMiedema.size();i++) Alloy2MiscibilityMiedema.push_back(b.Alloy2MiscibilityMiedema.at(i));
  Alloy2MiscibilityHumeRothery.clear(); for(uint i=0;i<b.Alloy2MiscibilityHumeRothery.size();i++) Alloy2MiscibilityHumeRothery.push_back(b.Alloy2MiscibilityHumeRothery.at(i));
  PhaseDiagramNameLeft.clear(); for(uint i=0;i<b.PhaseDiagramNameLeft.size();i++) PhaseDiagramNameLeft.push_back(b.PhaseDiagramNameLeft.at(i));
  PhaseDiagramNameRight.clear(); for(uint i=0;i<b.PhaseDiagramNameRight.size();i++) PhaseDiagramNameRight.push_back(b.PhaseDiagramNameRight.at(i));
  //
  for(uint i=0;i<ZLibrary.size();i++) { ZLibrary.at(i).clear(); } ZLibrary.clear();  // MORE
  for(uint i=0;i<b.ZLibrary.size();i++) ZLibrary.push_back(b.ZLibrary.at(i));
  //
  for(uint i=0;i<ZGNDlibrary.size();i++) { ZGNDlibrary.at(i).clear(); } ZGNDlibrary.clear();  // MORE
  for(uint i=0;i<b.ZGNDlibrary.size();i++) ZGNDlibrary.push_back(b.ZGNDlibrary.at(i));
  //
  for(uint i=0;i<ZMINElibrary.size();i++) { ZMINElibrary.at(i).clear(); } ZMINElibrary.clear();  // MORE
  for(uint i=0;i<b.ZMINElibrary.size();i++) ZMINElibrary.push_back(b.ZMINElibrary.at(i));
  //
  for(uint i=0;i<ZConcentrations.size();i++) { ZConcentrations.at(i).clear(); } ZConcentrations.clear();  // MORE
  for(uint i=0;i<b.ZConcentrations.size();i++) ZConcentrations.push_back(b.ZConcentrations.at(i));
  // complex systems
  enthalpy_formation_atom_max.clear(); for(uint i=0;i<b.enthalpy_formation_atom_max.size();i++) enthalpy_formation_atom_max.push_back(b.enthalpy_formation_atom_max.at(i));
  entropic_temperature_max.clear(); for(uint i=0;i<b.entropic_temperature_max.size();i++) entropic_temperature_max.push_back(b.entropic_temperature_max.at(i));
  alloy_order.clear(); for(uint i=0;i<b.alloy_order.size();i++) alloy_order.push_back(b.alloy_order.at(i));
  alloy_holes.clear(); for(uint i=0;i<b.alloy_holes.size();i++) alloy_holes.push_back(b.alloy_holes.at(i));
  alloy_deltamuBArich.clear(); for(uint i=0;i<b.alloy_deltamuBArich.size();i++) alloy_deltamuBArich.push_back(b.alloy_deltamuBArich.at(i));
  alloy_deltamuABrich.clear(); for(uint i=0;i<b.alloy_deltamuABrich.size();i++) alloy_deltamuABrich.push_back(b.alloy_deltamuABrich.at(i));
}

// copy
APENNSY_Parameters::APENNSY_Parameters(const APENNSY_Parameters& b) {
  //  free();
  // *this=b;
  copy(b);
}

// copies xstructures: b=a
const APENNSY_Parameters& APENNSY_Parameters::operator=(const APENNSY_Parameters& b) {  // operator=
  if(this!=&b) {
    free();
    copy(b);
  }
  return *this;
}

void APENNSY_Parameters::clear() {
  APENNSY_Parameters tmp;
  copy(tmp);
}

// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************

// ***************************************************************************
// This function loads the dynamic model to store all the data for the library
// ***************************************************************************
bool APENNSY_Parameters::LoadDynamicMemory(bool _verbose) {
  if(_verbose) {;} // dummy load
  // Nstructures
  for(uint i=0;i<alloys.size();i++) {
    // Nstructures
    // Nalloy*Nstructures
    ZConcentrations.push_back(*(new vector<double>(0)));
    ZLibrary.push_back(*(new vector<aflowlib::_aflowlib_entry>(0)));

    // alloys.size()*Nconcentrations[2]
    ZGNDlibrary.push_back(*(new vector<int>));   // generate only the alloys.size()
    ZMINElibrary.push_back(*(new vector<uint>));   // generate only the alloys.size()
    // alloys.size()*Nconcentrations[2]*X
    RankLib.push_back(*(new vector<vector<int> >));   // generate only the alloys.size()
    RankLibInt.push_back(*(new vector<vector<int> >));   // generate only the alloys.size()
    // alloys.size()*NstructuresHTQC*NstructuresHTQC
    AlloyStructureIdentity.push_back(*(new vector<vector<bool> >));   // generate only the alloys.size()
    // alloys.size()
    Alloy2MiscibilityHT.push_back(MISCIBILITY_SYSTEM_UNKNOWN);
    Alloy2MiscibilityExperiments.push_back(MISCIBILITY_SYSTEM_UNKNOWN);
    Alloy2MiscibilityMiedema.push_back(MISCIBILITY_SYSTEM_UNKNOWN);
    Alloy2MiscibilityHumeRothery.push_back(MISCIBILITY_SYSTEM_UNKNOWN);
    PhaseDiagramNameLeft.push_back("");
    PhaseDiagramNameRight.push_back("");
    // AlloyProperties* temp;temp=new AlloyProperties;ALibrary.push_back(*temp);
    enthalpy_formation_atom_max.push_back(double(0.0));
    entropic_temperature_max.push_back(double(0.0));
    alloy_order.push_back(double(0.0));
    alloy_holes.push_back(int(0));
    alloy_deltamuBArich.push_back(double(0.0));
    alloy_deltamuABrich.push_back(double(0.0));
    // alloys.size()*Nstructures
    //     AlloyStructureIdentity.push_back(temp); // moved to aflow_apennsy_gndstate.cpp
    //   // alloys.size()*Nconcentrations[2]
    //   ZMINElibrary.push_back(*(new vector<uint>(Nconcentrations[2]+1)));
    //   // alloys.size()*Nconcentrations[2]
  }

  // alloys.size()
  // RankLib.push_back(temp); // moved to aflow_apennsy_gndstate.cpp
  // RankLibInt.push_back(temp); // moved to aflow_apennsy_gndstate.cpp

  return TRUE;
}

// ***************************************************************************
// isequal with tolerance
// ***************************************************************************
//[CO20200404 - OBSOLETE]template <class utype> bool _isequal(const utype& c1,const utype& c2) {
//[CO20200404 - OBSOLETE]#define isequal_threshold 0.0001
//[CO20200404 - OBSOLETE]  if(abs(c1-c2)<isequal_threshold) return TRUE;
//[CO20200404 - OBSOLETE]  else return FALSE;
//[CO20200404 - OBSOLETE]}

// ***************************************************************************
// simple search string
// ***************************************************************************
bool strsub(string strzero, string subA) {
  return (bool) (strstr(strzero.c_str(),subA.c_str()));
}

// ***************************************************************************
// simple search string
// ***************************************************************************
bool strsub(string strzero, string subA, string subB) {
  return (bool) (strstr(strzero.c_str(),subA.c_str()) && strstr(strzero.c_str(),subB.c_str()));
}

// ***************************************************************************
// simple search string
// ***************************************************************************
bool strsub(string strzero, string subA, string subB, string subC) {
  return (bool) (strstr(strzero.c_str(),subA.c_str()) && strstr(strzero.c_str(),subB.c_str()) && strstr(strzero.c_str(),subC.c_str()));
}

// ***************************************************************************
// CleanNameStructur
// ***************************************************************************
string CleanNameStructure(string stringin) {
  string strout=stringin;
  strout=aurostd::RemoveSubStringFirst(strout,"/");
  strout=aurostd::RemoveSubStringFirst(strout,"/");
  return strout;
}

// ***************************************************************************
// APENNSY_Parameters::SplitAlloySpecies and APENNSY_Parameters::SplitAlloyPseudoPotentials
// ***************************************************************************
bool APENNSY_Parameters::SplitAlloySpecies(void) {
  for (uint k=0;k<alloys.size();k++) {
    speciesA.push_back("");speciesB.push_back("");
    KBIN::VASP_SplitAlloySpecies(alloys.at(k),speciesA.at(k),speciesB.at(k));
    if(PseudopotentialNoclean==TRUE) { speciesAB.push_back(alloys.at(k));
    } else { speciesAB.push_back(speciesA.at(k)+speciesB.at(k)); }
  }
  return TRUE;
}

bool APENNSY_Parameters::SplitAlloyPseudoPotentials(void) {
  for (uint k=0;k<alloys.size();k++) {
    pseudosA.push_back("");pseudosB.push_back("");
    KBIN::VASP_SplitAlloyPseudoPotentials(alloys.at(k),pseudosA.at(k),pseudosB.at(k));
  }
  return TRUE;
}

// ***************************************************************************
// Is_alphabetic - tells if an alloy is alphabetic or not
// ***************************************************************************
bool _NEW_is_alphabetic(string alloy_dir) {
  string in_file_name_aflowlib_entry=alloy_dir+string("/3/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT);
  vector<string> vaflowlib_entry,tokens;
  if(XHOST.APENNSY_USE_SERVER) { // FINDING ALPHABETIC FROM aflowlib.out - local
    if(aurostd::FileExist(in_file_name_aflowlib_entry)) {
      aurostd::string2tokens(aurostd::file2string(in_file_name_aflowlib_entry),vaflowlib_entry,"|");
    } else {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"file not found [1]: "+in_file_name_aflowlib_entry,_FILE_CORRUPT_); //CO20200624
    }
  }
  if(XHOST.APENNSY_SERVER_AFLOWLIB_ORG) { // FINDING ALPHABETIC FROM aflowlib.out - web
    aurostd::url2tokens(in_file_name_aflowlib_entry,vaflowlib_entry,"|"); 
    // [OBSOLETE]  aurostd::string2tokens(aurostd::url2string(in_file_name_aflowlib_entry),vaflowlib_entry,"|");
    if(!vaflowlib_entry.size()) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"file not found [2]: "+in_file_name_aflowlib_entry,_FILE_CORRUPT_); //CO20200624
    }
  }
  vector<string> vspeciesD,vspeciesP;
  for(uint i=0;i<vaflowlib_entry.size();i++) {
    if(aurostd::substring2bool(vaflowlib_entry.at(i),"aurl=")) {
      aurostd::string2tokens(aurostd::substring2string(vaflowlib_entry.at(i),"aurl=",1,TRUE),tokens,"/");
      for(uint j=0;j<tokens.size();j++) if(tokens.at(j)=="3") XATOM_SplitAlloySpecies(tokens.at(j-1),vspeciesD);
      if(vspeciesD.size()==0) {;}
      if(vspeciesD.size()==1) {;}
      if(vspeciesD.size()==2) {
        //	cout << "[5] " << "[" << vspeciesD.at(0) << "," << vspeciesD.at(1) << "]" << endl;
      }
      if(vspeciesD.size()>=3) {
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"unsupported number of species (vspecies) in "+in_file_name_aflowlib_entry,_INPUT_ILLEGAL_); //CO20200624
      }
    }
    if(aurostd::substring2bool(vaflowlib_entry.at(i),"species=") && !aurostd::substring2bool(vaflowlib_entry.at(i),"nspecies=")) {
      aurostd::string2tokens(aurostd::substring2string(vaflowlib_entry.at(i),"species=",1,TRUE),vspeciesP,",");
      if(vspeciesP.size()==0) {;}
      if(vspeciesP.size()==1) {;}
      if(vspeciesP.size()==2) {
        //	cout << "[5] " << "[" << vspeciesP.at(0) << "," << vspeciesP.at(1) << "]" << endl;
      }
      if(vspeciesP.size()>=3) {
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"unsupported number of species (vspeciesP) in "+in_file_name_aflowlib_entry,_INPUT_ILLEGAL_); //CO20200624
      }
    }
  }
  if(vspeciesD.size()==0 && vspeciesP.size()==0) return TRUE;
  if(vspeciesD.size()==1 && vspeciesP.size()==1 && vspeciesD.at(0)==vspeciesP.at(0)) return TRUE;
  if(vspeciesD.at(0)==vspeciesD.at(1) && vspeciesP.at(0)==vspeciesP.at(1) && vspeciesD.at(0)==vspeciesP.at(0) && vspeciesD.at(1)==vspeciesP.at(1)) return TRUE;
  if(vspeciesD.at(0) <vspeciesD.at(1) && vspeciesP.at(0) <vspeciesP.at(1) && vspeciesD.at(0)==vspeciesP.at(0) && vspeciesD.at(1)==vspeciesP.at(1)) return TRUE;
  if(vspeciesD.at(0) >vspeciesD.at(1) && vspeciesP.at(0) <vspeciesP.at(1) && vspeciesD.at(0)==vspeciesP.at(1) && vspeciesD.at(1)==vspeciesP.at(0)) return FALSE;

  stringstream message;
  message << vspeciesD.size() << " " << vspeciesP.size() << endl;
  message << vspeciesD.at(0) << " " << vspeciesD.at(1) << " " << vspeciesP.at(0) << " " << vspeciesP.at(1) << endl;
  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_); //CO20200624
  return FALSE;
}

bool _OLD_is_alphabetic(string alloy_dir) {
  FILE *in_file_pointer;
  string in_file_name;
  static char string_line[1024],*strptr;
  int i=0;
  char A1=0,A2=0,B1=0,B2=0;
  bool out=FALSE;

  // FINDING CONCENTRATIONS FROM POSCAR.relax1  // OBSOLETE
  in_file_name=alloy_dir+string("/3/POSCAR.relax1"); // OBSOLETE
  in_file_pointer=fopen(in_file_name.c_str(),"r");
  if(in_file_pointer==NULL) {
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"file not found: "+in_file_name,_FILE_CORRUPT_); //CO20200624
  }

  while(fgets(string_line,1024,in_file_pointer)) {
    i++;
    if(i==1) {
      A1=string_line[0];A2=string_line[1];
      if(A1<65 || A1>90) {
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Error alphabetic (1)",_INPUT_ILLEGAL_); //CO20200624
      }
      A1=A1-65;
      if(A2<97 || A2>122)  A2=0; else A2=A2-97;
    }
    if(i==8) {
      strtod(string_line,&strptr);
      strtod(strptr,&strptr);
      strtod(strptr,&strptr);
      while(strptr[0]==32) strptr++;
      B1=strptr[0];B2=strptr[1];
      if(B1<65 || B1>90) {
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Error alphabetic (2)",_INPUT_ILLEGAL_); //CO20200624
      }
      B1=B1-65;
      if(B2<97 || B2>122)  B2=0; else B2=B2-97;
      //  cerr << strptr << endl;
    }
    // if(i>=8) fprintf(stderr,"%i %i %i %i\n",A1,A2,B1,B2);
    // cerr << i << " " << string_line; cerr.flush();
  }
  fclose(in_file_pointer);
  //  cerr << endl << A1 << " " << A2 << " " << B1 << " " << B2 << endl;
  if(A1==B1 && A2==B2) out=TRUE;
  else out=FALSE;
  //  cerr << alloy_dir << " " << out << endl;
  return (bool) out;
}

bool _is_alphabetic(string alloy_dir) {
  return _NEW_is_alphabetic(alloy_dir);
}

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
