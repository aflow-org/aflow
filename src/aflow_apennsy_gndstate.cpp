// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

#ifndef _AFLOW_APENNSY_GNDSTATE_CPP_
#define _AFLOW_APENNSY_GNDSTATE_CPP_
#include "aflow.h"
#include "aflow_apennsy.h"
#include "aflowlib.h"

#define _APENNSY_STYLE_OLD_  FALSE
#define _APENNSY_STYLE_NEW_  TRUE
#define INF 1E9

extern vector<string> vLibrary_ALL;

// DEVIL DIRS: dir!="64" && dir!="65" && dir!="549" && dir!="550" && dir!="f8269" && dir!="f9083" && dir!="f8819") 

// Nconcentrations[2]=ZConcentrations.at(nalloy).size()-1
// NstructuresHTQC[2]=ZLibrary.at(nalloy).size()

// **************************************************************************
// LoadLibrary - this is the main function loading all the energies
// **************************************************************************
bool APENNSY_Parameters::LoadLibrary(_aflags &aflags) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream message;
  // START_NEW_STYLE 
  // LOAD LIBRARY
  if(XHOST.APENNSY_USE_LIBRARY) aflowlib::LOAD_Library_ALL(TRUE);
  // END_NEW_STYLE
  if(0)cerr << aflags.Directory << endl;  // warning: unused parameter ‘aflags’ [-Wunused-parameter]

  aflowlib::_aflowlib_entry struct_proprts;
  LIBstructures_calcs=0;
  LIBstructures_holes=0;

  string alloy;vector<string> tokens;
  bool alphabetic;
  // xvector<int> gndstates(Nconcentrations[2]);
#ifdef TMELTING
  FILE *in_file_pointer_TMELTING;
  char string_line_TMELTING[1024],alloy_TMELTING[1024],*string_TMELTING_ptr;
  in_file_pointer_TMELTING=fopen(TMELTING_FILE,"r");
  if(in_file_pointer_TMELTING==NULL) {
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"file not found [M]: "+MELTING_FILE,_FILE_CORRUPT_); //CO20200624
  }
  fgets(string_line_TMELTING,1024,in_file_pointer_TMELTING);
  cerr << string_line_TMELTING << endl;
#endif
  for(uint nalloy=0;nalloy<alloys.size();nalloy++)
    //  cerr << "DEBUG [0]  " << "Nalloy=" << nalloy << endl;
    // for(nalloy=1;nalloy<=1;nalloy++)
  { //CO20200106 - patching for auto-indenting
#ifdef TMELTING
    fgets(string_line_TMELTING,1024,in_file_pointer_TMELTING);
    if(!strstr(string_line_TMELTING,alloys.at(nalloy))) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"string "+alloys.at(nalloy)+" not found in "+string_line_TMELTING,_FILE_CORRUPT_); //CO20200624
    }
    for(uint i=1;i<=Ntemperatures;i++) {
      if(i==1) Tmelting(i,nalloy)=strtod(string_line_TMELTING+10,&string_TMELTING_ptr);
      else Tmelting(i,nalloy)=strtod(string_TMELTING_ptr,&string_TMELTING_ptr); // A .25Lb.25Hb .5Lb .5Hb.75Lb.75Hb   B
      // cerr << Tmelting(i,nalloy) << "\t";  // crazy but true
    }
    Ta=Tmelting(1,nalloy);Tb=Tmelting(8,nalloy);
    TFmelting(1,nalloy)=0;
    TFmelting(2,nalloy)=Tmelting(2,nalloy)-0.75*Ta-0.25*Tb;
    TFmelting(3,nalloy)=Tmelting(3,nalloy)-0.75*Ta-0.25*Tb;
    TFmelting(4,nalloy)=Tmelting(4,nalloy)-0.50*Ta-0.50*Tb;
    TFmelting(5,nalloy)=Tmelting(5,nalloy)-0.50*Ta-0.50*Tb;
    TFmelting(6,nalloy)=Tmelting(4,nalloy)-0.25*Ta-0.75*Tb;
    TFmelting(7,nalloy)=Tmelting(5,nalloy)-0.25*Ta-0.75*Tb;
    TFmelting(8,nalloy)=0;
    //      for(i=1;i<=8;i++) cerr << TFmelting(i,nalloy) << "\t";  // crazy but true
    for(i=1;i<=4;i++) {
      TmeltingMinMax(i,nalloy)=strtod(string_TMELTING_ptr,&string_TMELTING_ptr);
      //	 cerr << TmeltingMinMax(i,nalloy) << "\t";  // crazy but true
    }
#endif
    alloy=alloysRAW.at(nalloy);
    vector<string> tokens_tmp;
    aurostd::string2tokens(alloy,tokens_tmp,"/");
    Alloy2MiscibilityExperiments.at(nalloy)=MiscibilityExperimentsCheck(tokens_tmp.at(tokens_tmp.size()-1));
    string alloy_tmp=tokens_tmp.at(tokens_tmp.size()-1);
    string pseudoA,pseudoB,speciesA,speciesB;
    KBIN::VASP_SplitAlloyPseudoPotentials(alloy_tmp,pseudoA,pseudoB);
    KBIN::VASP_SplitAlloySpecies(alloy_tmp,speciesA,speciesB);

    // START_NEW_STYLE
    // grepping inside the LIRARY
    // uint aflowlib::GREP_Species_ALL(vector<string> vspecies,                   // IN   [0,nspecies[ the species Ag,Cd,...   nspecies=number of these items    nspecies=naries
    //                vector<string>& vspecies_pp,                      // IN   [0,nspecies[ the pseudopotentials Ag_pv, Cd_sv
    // 		      vector<vector<string> >& vList,                   // OUT  [0,naries[*[0,vList.size()[ returns the lines of the library containing A,B,C,AB,AC,BC,ABC....
    // 		      vector<vector<vector<string> > > &vList_tokens,   // OUT  [0,naries[*[0,vList.size()[*[0,size_tokens[ returns the tokens for each line of the vList
    // 		      vector<vector<vector<string> > > &vList_species,  // OUT  [0,naries[*[0,vList.size()[*nspecies  returns the species present
    // 		      vector<double> &vList_Hmin,                       // OUT  [0,nspecies[ returns the min enthalpy for reference
    // 		      vector<string> &vList_Pmin,                       // OUT  [0,nspecies[ returns the prototype for reference
    // 		      vector<uint> &vList_Imin,                         // OUT  [0,nspecies[ returns the line index of vList.at(ispecies), in which we have the min enthalpy for reference
    // 		      vector<vector<vector<double> > > &vList_concs,    // OUT  [0,naries[*[0,vList.size()[*[0.nspecies[ the concentrations AxAyCz... where x+y+z=1 and it contains also ZEROS so that 0 0.25 0.75 is allowed
    // 		      vector<vector<double> > &vList_Ef) {              // OUT  [0,naries[*[0,vList.size()[ returns the formation energy of the list  //[CO20200106 - close bracket for indenting]}
    //vector<string> vspecies;                         // IN   [0,nspecies[ the species Ag,Cd,...   nspecies=number of these items    nspecies=naries
    //vector<string> vspecies_pp;                      // IN   [0,nspecies[ the pseudopotentials Ag_pv, Cd_sv
    //vector<vector<string> > vList;                   // OUT  [0,naries[*[0,vList.size()[ returns the lines of the library containing A,B,C,AB,AC,BC,ABC....
    //vector<vector<vector<string> > > vList_tokens;   // OUT  [0,naries[*[0,vList.size()[*[0,size_tokens[ returns the tokens for each line of the vList
    //vector<vector<vector<string> > > vList_species;  // OUT  [0,naries[*[0,vList.size()[*nspecies  returns the species present
    //vector<double> vList_Hmin;                       // OUT  [0,nspecies[ returns the min enthalpy for reference
    //vector<string> vList_Pmin;                       // OUT  [0,nspecies[ returns the prototype for reference
    //vector<uint> vList_Imin;                         // OUT  [0,nspecies[ returns the line index of vList.at(ispecies), in which we have the min enthalpy for reference
    //vector<vector<vector<double> > > vList_concs;    // OUT  [0,naries[*[0,vList.size()[*[0.nspecies[ the concentrations AxAyCz... where x+y+z=1 and it contains also ZEROS so that 0 0.25 0.75 is allowed
    //vector<vector<double> > vList_Ef;                // OUT  [0,naries[*[0,vList.size()[ returns the formation energy of the list

    // START_NEW_STYLE
    // grepping inside the LIRARY
    // uint aflowlib::GREP_Species_ALL(vector<string> vspecies,                   // IN   [0,nspecies[ the species Ag,Cd,...   nspecies=number of these items    nspecies=naries
    //                vector<string>& vspecies_pp,                      // IN   [0,nspecies[ the pseudopotentials Ag_pv, Cd_sv
    // 		      vector<vector<string> >& vList);                   // OUT  [0,naries[*[0,vList.size()[ returns the lines of the library containing A,B,C,AB,AC,BC,ABC....
    vector<string> vspecies;                         // IN   [0,nspecies[ the species Ag,Cd,...   nspecies=number of these items    nspecies=naries
    vector<string> vspecies_pp;                      // IN   [0,nspecies[ the pseudopotentials Ag_pv, Cd_sv
    vector<vector<string> > vList;                   // OUT  [0,naries[*[0,vList.size()[ returns the lines of the library containing A,B,C,AB,AC,BC,ABC....

    vspecies.push_back(speciesA);vspecies.push_back(speciesB);
    vspecies_pp.push_back(pseudoA);vspecies_pp.push_back(pseudoB);

    if(LDEBUG) cerr << " Creating concentrations" << endl;
    ZConcentrations.at(nalloy).clear();ZConcentrations.at(nalloy).push_back(-1.0);
    vector<string> vaflowlibentry_url,vaflowlibentry_url_aus;   
    vaflowlibentry_url.clear();vaflowlibentry_url_aus.clear();
    vector<aflowlib::_aflowlib_entry> vaflowlib,vaflowlib_aus;
    vaflowlib.clear();vaflowlib_aus.clear();
    aflowlib::_aflowlib_entry _aflowlib_tmp;

    // LOAD DATA !
    if(XHOST.APENNSY_USE_LIBRARY) { // from library in LOCAL SERVER
      if(LDEBUG) cerr << " XHOST.APENNSY_USE_LIBRARY" << endl;
      // aflowlib::GREP_Species_ALL(vspecies,vspecies_pp,vList,vList_tokens,vList_species,vList_Hmin,vList_Pmin,vList_Imin,vList_concs,vList_Ef);
      aflowlib::GREP_Species_ALL(vspecies,vspecies_pp,vList);
      if(LDEBUG) cerr << " LIBS entries= " << vList.at(0).size()+vList.at(1).size() << endl;
      for(uint ilist=0;ilist<vList.at(0).size();ilist++) {_aflowlib_tmp.Load(vList.at(0).at(ilist),cout);vaflowlib_aus.push_back(_aflowlib_tmp);};
      for(uint ilist=0;ilist<vList.at(1).size();ilist++) {_aflowlib_tmp.Load(vList.at(1).at(ilist),cout);vaflowlib_aus.push_back(_aflowlib_tmp);};
    } 
    // from directory in LOCAL SERVER
    vector<string> vfiles;
    if(XHOST.APENNSY_USE_SERVER) {
      if(LDEBUG) cerr << " XHOST.APENNSY_USE_SERVER" << endl;
      aurostd::string2tokens(aurostd::execute2string(string("find "+string(alloy)+" -name "+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)),vfiles,"\n");
      for(uint ilist=0;ilist<vfiles.size();ilist++) 
        if(vfiles.at(ilist)!="") {_aflowlib_tmp.file2aflowlib(vfiles.at(ilist),cout); vaflowlib_aus.push_back(_aflowlib_tmp);}
      for(uint ilist=0;ilist<vfiles.size();ilist++)
        if(vfiles.at(ilist)!="") {vaflowlibentry_url_aus.push_back(vfiles.at(ilist));}
    }
    // from directory in AFLOWLIB WEB SERVER
    if(XHOST.APENNSY_SERVER_AFLOWLIB_ORG) {
      if(LDEBUG) cerr << " XHOST.APENNSY_SERVER_AFLOWLIB_ORG" << endl;
      // [OBSOLETE]   aurostd::url2tokens(string(alloy)+"/"+_XENTRY_+"?aflowlib_entries",vfiles,","); 
      aurostd::url2tokens(string(alloy)+"/?aflowlib_entries",vfiles,","); 
      for(uint ilist=0;ilist<vfiles.size();ilist++) 
        if(vfiles.at(ilist)!="") {
          string url=string(alloy)+"/"+vfiles.at(ilist)+"/";
          aflowlib::_aflowlib_entry _aflowlib_tmp;
          _aflowlib_tmp.url2aflowlib(url,cout);
          if(!_aflowlib_tmp.entry.empty()) {
            vaflowlib_aus.push_back(_aflowlib_tmp);
            vaflowlibentry_url_aus.push_back(url);
            //	    vaflowlib.push_back(_aflowlib_tmp);
          }
        }  
    }

    for(uint ilist=0;ilist<vaflowlib_aus.size();ilist++) {
      bool goodENTRY=TRUE;
      if(aurostd::substring2bool(vaflowlib_aus.at(ilist).entry,"LDAU")) goodENTRY=FALSE;
      if(aurostd::substring2bool(vaflowlib_aus.at(ilist).entry,"PAW_GGA") && !aurostd::substring2bool(vaflowlib_aus.at(ilist).entry,"B_hSm_3")) goodENTRY=FALSE;
      if(aurostd::substring2bool(vaflowlib_aus.at(ilist).entry,"NUPDOWN")) goodENTRY=FALSE;
      if(goodENTRY) vaflowlib.push_back(vaflowlib_aus.at(ilist));
      if(goodENTRY) vaflowlibentry_url.push_back(vaflowlibentry_url_aus.at(ilist));
    }
    if(LDEBUG) cerr << "vaflowlib_aus.size()=" << vaflowlib_aus.size() << endl;
    if(LDEBUG) cerr << "vaflowlib.size()=" << vaflowlib.size() << endl;
    if(LDEBUG) cerr << "vaflowlibentry_url_aus.size()=" << vaflowlibentry_url_aus.size() << endl;
    if(LDEBUG) cerr << "vaflowlibentry_url.size()=" << vaflowlibentry_url.size() << endl;

    // PREPARE CONCENTRATIONS
    for(uint ilist=0;ilist<vaflowlib.size();ilist++) {
      string prototype="";string species="";int natoms,nspecies;vector<int> composition;
      double Cb=0.0;
      //  [OBSOLETE]  aflowlib::TokenExtractAFLOWLIB(vaflowlibentry.at(ilist),"prototype=",prototype); // cerr << "(1) prototype=" << prototype << endl;
      // [OBSOLETE]   aflowlib::TokenExtractAFLOWLIB(vaflowlib.at(ilist).entry,"species=",species); // cerr << "(1) species=" << species << " " << speciesA << "," << speciesB << endl;
      // [OBSOLETE] aflowlib::TokenExtractAFLOWLIB(vaflowlib.at(ilist).entry,"nspecies=",nspecies); // cerr << "(1) nspecies=" << nspecies << endl;
      species=vaflowlib.at(ilist).species;
      nspecies=vaflowlib.at(ilist).nspecies;
      if(nspecies==1) {
        if(species==speciesB) Cb=1.0;
        if(species==speciesA) Cb=0.0;
      }
      if(nspecies==2) {
        // [OBSOLETE] aflowlib::TokenExtractAFLOWLIB(vaflowlib.at(ilist).entry,"composition=",composition); // cerr << "(1) composition=" << composition << endl;
        // [OBSOLETE] aflowlib::TokenExtractAFLOWLIB(vaflowlib.at(ilist).entry,"natoms=",natoms); // cerr << "(1) natoms=" << natoms << endl;
        for(uint i=0;i<vaflowlib.at(ilist).vcomposition.size();i++) composition.push_back(vaflowlib.at(ilist).vcomposition.at(i));
        natoms=vaflowlib.at(ilist).natoms;
        Cb=((double) composition.at(1))/((double) natoms);
      }
      bool found=FALSE;
      for(uint ii=0;ii<ZConcentrations.at(nalloy).size()&&!found;ii++) found=(aurostd::abs(ZConcentrations.at(nalloy).at(ii)-Cb)<apennsy_epsilon);
      // check if neglection will affect number of concentrations
      if(aflags.vflag.flag("APENNSY::NEGLECT_STRUCTURES")) 
        for(uint i=0;i<aflags.APENNSY_NEGLECT_STRUCTURES_vstrs.size();i++) 
          if(vaflowlib.at(ilist).prototype==aflags.APENNSY_NEGLECT_STRUCTURES_vstrs.at(i))
            found=TRUE; // it will be skipped
      // if everything is ok then:
      if(!found) ZConcentrations.at(nalloy).push_back(Cb);
    }
    std::sort(ZConcentrations.at(nalloy).begin(),ZConcentrations.at(nalloy).end());
    if(LDEBUG) cerr << "ZConcentrations.at(nalloy).size()=" << ZConcentrations.at(nalloy).size() << endl;

    // PREPARE DYNAMIC DATA
    xvector<int> gndstates(ZConcentrations.at(nalloy).size()-1);    
    // alloys.size()*ZConcentrations.size() AlloyStructureIdentity
    ZGNDlibrary.at(nalloy).clear(); for(uint ii=0;ii<ZConcentrations.at(nalloy).size()-1+1;ii++) ZGNDlibrary.at(nalloy).push_back((uint) 0);
    if(LDEBUG) cerr << "ZGNDlibrary.at(nalloy).size()=" << ZGNDlibrary.at(nalloy).size() << endl;
    ZMINElibrary.at(nalloy).clear(); for(uint ii=0;ii<ZConcentrations.at(nalloy).size()-1+1;ii++) ZMINElibrary.at(nalloy).push_back((uint) 0);
    if(LDEBUG) cerr << "ZMINElibrary.at(nalloy).size()=" << ZMINElibrary.at(nalloy).size() << endl;
    // alloys.size()*ZConcentrations.size()*X
    RankLib.at(nalloy).clear();
    RankLibInt.at(nalloy).clear();
    for(uint j=0;j<ZConcentrations.at(nalloy).size();j++) {
      RankLib.at(nalloy).push_back(*(new vector<int>)); // alloys.size()*ZConcentrations.size()*X // MaxNumberConcentrations[2]
      RankLibInt.at(nalloy).push_back(*(new vector<int>)); // alloys.size()*ZConcentrations.size()*X // MaxNumberConcentrations[2]
    }
    if(LDEBUG) cerr << "RankLib.size()=" << RankLib.size() << endl;
    if(LDEBUG) cerr << "RankLibInt.size()=" << RankLibInt.size() << endl;
    if(LDEBUG) cerr << "RankLib.at(nalloy).size()=" << RankLib.at(nalloy).size() << endl;
    if(LDEBUG) cerr << "RankLibInt.at(nalloy).size()=" << RankLibInt.at(nalloy).size() << endl;

    //   DEBUG=TRUE;
    // cerr << "*" << pseudoA << "*" << pseudoB << endl;
    // cerr << "*************************************** "<< endl;
    cerr << "GNDSTATE: doing " << (nalloy<100 ? " " : "") << (nalloy<10 ? " " : "") << nalloy << " "; // << alloy;// << endl;
    alphabetic=_is_alphabetic(alloy);
    // alphabetic=is_alphabetic(alloys.at(nalloy));
    if(alphabetic)  cerr << "ALPHABETIC     ";
    if(!alphabetic) cerr << "NON ALPHABETIC ";
    // -------------------------------------------------------------------------------------
    // LOAD CLUSTER EXPANSION DATA !
    alloy_holes.at(nalloy)=0;

    if(LDEBUG) cerr << "AlloyStructureIdentity.size()=" << AlloyStructureIdentity.size() << endl;
    if(LDEBUG) cerr << "AlloyStructureIdentity.at(nalloy).size()=" << AlloyStructureIdentity.at(nalloy).size() << endl;
    if(LDEBUG) cerr << "vaflowlib.size()=" << vaflowlib.size() << endl;

    ZLibrary.at(nalloy).clear();


    for(uint i=0;i<vaflowlib.size();i++) { // i in vaflowlib.size()
      if(LDEBUG) cerr << "i=" << i << endl;
      aflowlib::_aflowlib_entry plug,struct_proprts;
      plug=vaflowlib.at(i);
      // [OBSOLETE]   plug_entry=aurostd::CleanStringASCII(plug_entry);
      string prototype="";string natoms="";string composition="";string species="";
      // [OBSOLETE] aflowlib::TokenExtractAFLOWLIB(plug.entry,"prototype=",prototype);
      prototype=plug.prototype;
      if(LDEBUG)  cerr << i << "   " << alloy << "/" << prototype << endl;
      cerr.flush();
      // ENERGY PROPERTIES

      //    if(isavailable==TRUE && strnumber!=INF) {
      // 	if((Load_FCC==FALSE && Load_BCC==FALSE && Load_HCP==FALSE) || Load_ALL==TRUE || prototype=="3") struct_proprts.GetVASP(alloy,i,ishole,alphabetic,*this);
      // 	// check if FCC BCC HCP
      // 	if(Load_FCC==TRUE && (prototype.at(0)=='f' || (structures[2].at(strnumber).fcc==TRUE))) struct_proprts.GetVASP(alloy,i,ishole,alphabetic,*this);
      // 	if(Load_BCC==TRUE && (prototype.at(0)=='b' || (structures[2].at(strnumber).bcc==TRUE))) struct_proprts.GetVASP(alloy,i,ishole,alphabetic,*this);
      // 	if(Load_HCP==TRUE && (prototype.at(0)=='h' || (structures[2].at(strnumber).hcp==TRUE))) struct_proprts.GetVASP(alloy,i,ishole,alphabetic,*this);
      // 	// REMOVE
      //      }
      bool ishole=FALSE;
      // cerr << pseudoA << "," << pseudoB << "," << prototype << "," << endl;
      // if(pseudoB=="Ag") cerr << alloy << endl;
      if(prototype=="64" || prototype=="65" || prototype=="549" || prototype=="550" || prototype=="f8269" || prototype=="f9083" || prototype=="f8819") ishole=TRUE; // DEVIL
      if((pseudoA=="Ag" && prototype=="303") || (pseudoB=="Ag" && prototype=="304")) ishole=TRUE;  // bad Ag is a wrong relaxation
      if((pseudoA=="Ag" && prototype=="323") || (pseudoB=="Ag" && prototype=="324")) ishole=TRUE;  // bad Ag is a wrong relaxation
      if((pseudoA=="Au" && prototype=="323") || (pseudoB=="Au" && prototype=="324")) ishole=TRUE;  // bad Au is a wrong relaxation
      if((pseudoA=="Al_h" && prototype=="307") || (pseudoB=="Al_h" && prototype=="308")) ishole=TRUE;  // bad Al_h pseudopotential !
      if((pseudoA=="Al_h" && prototype=="A7.A") || (pseudoB=="Al_h" && prototype=="A7.B")) ishole=TRUE;  // bad Al_h pseudopotential !
      if((pseudoA=="Al_h" && prototype=="323") || (pseudoB=="Al_h" && prototype=="324")) ishole=TRUE;  // bad Al_h pseudopotential !
      if((pseudoA=="Ca_sv" && prototype=="303") || (pseudoB=="Ca_sv" && prototype=="304")) ishole=TRUE;  // bad Ca_sv is a wrong relaxation
      if((pseudoA=="Ca_sv" && prototype=="323") || (pseudoB=="Ca_sv" && prototype=="324")) ishole=TRUE;  // bad Ca_sv is a wrong relaxation
      if((pseudoA=="Cd" && prototype=="323") || (pseudoB=="Cd" && prototype=="324")) ishole=TRUE;  // bad Cd is a wrong relaxation
      if((pseudoA=="Cu_pv" && prototype=="303") || (pseudoB=="Cu_pv" && prototype=="304")) ishole=TRUE;  // bad Cu_pv is a wrong relaxation
      if((pseudoA=="Cu_pv" && prototype=="323") || (pseudoB=="Cu_pv" && prototype=="324")) ishole=TRUE;  // bad Cu_pv is a wrong relaxation
      if((pseudoA=="Fe_pv" && prototype=="307") || (pseudoB=="Fe_pv" && prototype=="308")) ishole=TRUE;  // bad Fe_pv is a wrong relaxation
      if((pseudoA=="Fe_pv" && prototype=="A7.A") || (pseudoB=="Fe_pv" && prototype=="A7.B")) ishole=TRUE;  // bad Fe_pv is a wrong relaxation
      if((pseudoA=="Ge_h" && prototype=="305") || (pseudoB=="Ge_h" && prototype=="306")) ishole=TRUE;  // bad Ge_h is a wrong relaxation
      if((pseudoA=="In_d" && prototype=="323") || (pseudoB=="In_d" && prototype=="324")) ishole=TRUE;  // bad In_d is a wrong relaxation
      if((pseudoA=="Ir" && prototype=="303") || (pseudoB=="Ir" && prototype=="304")) ishole=TRUE;  // bad Ir is a wrong relaxation
      if((pseudoA=="K_sv" && prototype=="307") || (pseudoB=="K_sv" && prototype=="308")) ishole=TRUE;  // bad K_sv is a wrong relaxation
      if((pseudoA=="K_sv" && prototype=="A7.A") || (pseudoB=="K_sv" && prototype=="A7.B")) ishole=TRUE;  // bad K_sv is a wrong relaxation
      if((pseudoA=="La" && prototype=="303") || (pseudoB=="La" && prototype=="304")) ishole=TRUE;  // bad La is a wrong relaxation
      if((pseudoA=="La" && prototype=="323") || (pseudoB=="La" && prototype=="324")) ishole=TRUE;  // bad La is a wrong relaxation
      if((pseudoA=="Li_sv" && prototype=="307") || (pseudoB=="Li_sv" && prototype=="308")) ishole=TRUE;  // bad Li_sv is a wrong relaxation
      if((pseudoA=="Li_sv" && prototype=="A7.A") || (pseudoB=="Li_sv" && prototype=="A7.B")) ishole=TRUE;  // bad Li_sv is a wrong relaxation
      if((pseudoA=="Na_pv" && prototype=="307") || (pseudoB=="Na_pv" && prototype=="308")) ishole=TRUE;  // bad Na_pv is a wrong relaxation
      if((pseudoA=="Na_pv" && prototype=="A7.A") || (pseudoB=="Na_pv" && prototype=="A7.B")) ishole=TRUE;  // bad Na_pv is a wrong relaxation
      if((pseudoA=="Ni_pv" && prototype=="303") || (pseudoB=="Ni_pv" && prototype=="304")) ishole=TRUE;  // bad Ni_pv is a wrong relaxation
      if((pseudoA=="Ni_pv" && prototype=="323") || (pseudoB=="Ni_pv" && prototype=="324")) ishole=TRUE;  // bad Ni_pv is a wrong relaxation
      if((pseudoA=="Pb_d" && prototype=="303") || (pseudoB=="Pb_d" && prototype=="304")) ishole=TRUE;  // bad Pb_d is a wrong relaxation
      if((pseudoA=="Pb_d" && prototype=="323") || (pseudoB=="Pb_d" && prototype=="324")) ishole=TRUE;  // bad Pb_d is a wrong relaxation
      if((pseudoA=="Pd_pv" && prototype=="303") || (pseudoB=="Pd_pv" && prototype=="304")) ishole=TRUE;  // bad Pd_pv is a wrong relaxation
      if((pseudoA=="Pd_pv" && prototype=="323") || (pseudoB=="Pd_pv" && prototype=="324")) ishole=TRUE;  // bad Pd_pv is a wrong relaxation
      if((pseudoA=="Pt" && prototype=="303") || (pseudoB=="Pt" && prototype=="304")) ishole=TRUE;  // bad Pt is a wrong relaxation
      if((pseudoA=="Pt" && prototype=="317") || (pseudoB=="Pt" && prototype=="318")) ishole=TRUE;  // bad Pt is a wrong relaxation

      if((pseudoA=="Rh_pv" && prototype=="303") || (pseudoB=="Rh_pv" && prototype=="304")) ishole=TRUE;  // bad Rh_pv is a wrong relaxation
      if((pseudoA=="Si_h" && prototype=="305") || (pseudoB=="Si_h" && prototype=="306")) ishole=TRUE;  // bad Si_h is a wrong relaxation
      if((pseudoA=="Si_h" && prototype=="307") || (pseudoB=="Si_h" && prototype=="308")) ishole=TRUE;  // bad Si_h is a wrong relaxation
      if((pseudoA=="Si_h" && prototype=="A7.A") || (pseudoB=="Si_h" && prototype=="A7.B")) ishole=TRUE;  // bad Si_h is a wrong relaxation
      if((pseudoA=="Si_h" && prototype=="323") || (pseudoB=="Si_h" && prototype=="324")) ishole=TRUE;  // bad Si_h is a wrong relaxation
      if((pseudoA=="Ta_pv" && prototype=="307") || (pseudoB=="Ta_pv" && prototype=="308")) ishole=TRUE;  // bad Ta_pv is a wrong relaxation
      if((pseudoA=="Ta_pv" && prototype=="A7.A") || (pseudoB=="Ta_pv" && prototype=="A7.B")) ishole=TRUE;  // bad Ta_pv is a wrong relaxation
      if((pseudoA=="B_h" && prototype=="317") || (pseudoB=="B_h" && prototype=="318")) ishole=TRUE;  // bad B_h is a wrong relaxation

      //   cerr << prototype << " " << struct_proprts.energy_atom << endl;
      if(ishole==FALSE) {
        LIBstructures_calcs++;
        plug.vstr.clear();
        //	if(struct_proprts.vstr.size()!=0 && struct_proprts.vstr.size()!=3) {
        //	  cerr << "ERROR: struct_proprts.vstr.size()=" << struct_proprts.vstr.size() << endl;
        //    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"struct_proprts.vstr.size()="+aurostd::utype2string(struct_proprts.vstr.size()),_INDEX_MISMATCH_); //CO20200624
        //	}
        xstructure xstr;
        if(XHOST.APENNSY_USE_SERVER) {
          string aurl;	
          // [OBSOLETE] aflowlib::TokenExtractAFLOWLIB(plug.entry,"aurl=",aurl);
          aurl=plug.aurl;
          aurl=aurostd::RemoveSubString(aurl,"nietzsche.mems.duke.edu:");aurl=aurostd::RemoveSubString(aurl,"nietzsche:");
          aurl=aurostd::RemoveSubString(aurl,"materials.duke.edu:");aurl=aurostd::RemoveSubString(aurl,"materials:");
          aurl=aurostd::RemoveSubString(aurl,"aflowlib.mems.duke.edu:");
          aurl=aurostd::RemoveSubString(aurl,"aflowlib.duke.edu:");aurl=aurostd::RemoveSubString(aurl,"aflowlib:");
          aurl=aurostd::StringSubst(aurl,"/LIB/","/RAW/");
          aurl=aurostd::StringSubst(aurl,"AFLOWDATA/LIB2_RAW/","/common/LIB2/RAW/");
          // PRE
          if(!aurostd::FileExist(aurl+"/POSCAR.relax1") && !aurostd::FileExist(aurl+"/POSCAR.orig")) {
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"aurl=\""+aurl+"\" needed POSCAR.orig/POSCAR.relax1",_FILE_CORRUPT_); //CO20200624
          }
          if(!xstr.atoms.size() && aurostd::FileExist(aurl+"/POSCAR.orig")) {xstr=xstructure(aurl+"/POSCAR.orig",IOVASP_AUTO);}
          if(!xstr.atoms.size() && aurostd::FileExist(aurl+"/POSCAR.relax1")) {xstr=xstructure(aurl+"/POSCAR.relax1",IOVASP_AUTO);}
          plug.vstr.push_back(xstr);xstr.clear(); //DX20191220 - uppercase to lowercase clear
          // MID
          if(!aurostd::FileExist(aurl+"/CONTCAR.relax1") && !aurostd::FileExist(aurl+"/POSCAR.relax2")) {
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"aurl=\""+aurl+"\" needed POSCAR.relax2/CONTCAR.relax1",_FILE_CORRUPT_); //CO20200624
          }
          if(!xstr.atoms.size() && aurostd::FileExist(aurl+"/POSCAR.relax2")) {xstr=xstructure(aurl+"/POSCAR.relax2",IOVASP_AUTO);}
          if(!xstr.atoms.size() && aurostd::FileExist(aurl+"/CONTCAR.relax1")) {xstr=xstructure(aurl+"/CONTCAR.relax1",IOVASP_AUTO);}
          plug.vstr.push_back(xstr);xstr.clear(); //DX20191220 - uppercase to lowercase clear
          // POST
          if(!aurostd::FileExist(aurl+"/CONTCAR.relax") && !aurostd::FileExist(aurl+"/CONTCAR.relax2")) {
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"aurl=\""+aurl+"\" needed CONTCAR.relax/CONTCAR.relax2",_FILE_CORRUPT_); //CO20200624
          }
          if(!xstr.atoms.size() && aurostd::FileExist(aurl+"/CONTCAR.relax")) {xstr=xstructure(aurl+"/CONTCAR.relax",IOVASP_AUTO);}
          if(!xstr.atoms.size() && aurostd::FileExist(aurl+"/CONTCAR.relax2")) {xstr=xstructure(aurl+"/CONTCAR.relax2",IOVASP_AUTO);}
          plug.vstr.push_back(xstr);xstr.clear(); //DX20191220 - uppercase to lowercase clear
          //	plug.vstr.clear(); for(uint ii=0;ii<struct_proprts.vstr.size();ii++) plug.vstr.push_back(struct_proprts.vstr.at(ii));
          //	plug.nrelaxations=struct_proprts.nrelaxations;
        }

        if(XHOST.APENNSY_SERVER_AFLOWLIB_ORG) { 	  // load files
          // PRE
          if(!aurostd::substring2bool(vaflowlib.at(i).vfiles,"POSCAR.relax1") && !aurostd::substring2bool(vaflowlib.at(i).vfiles,"POSCAR.orig")) {
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"vaflowlibentry_url.at(i)=\""+vaflowlibentry_url.at(i)+"\" needed POSCAR.orig/POSCAR.relax1",_FILE_CORRUPT_); //CO20200624
          }
          if(!xstr.atoms.size() && aurostd::substring2bool(vaflowlib.at(i).vfiles,"POSCAR.orig")) {
            stringstream aus(aurostd::url2string(vaflowlibentry_url.at(i)+"/POSCAR.orig"));xstr=xstructure(aus,IOVASP_AUTO);}
          if(!xstr.atoms.size() && aurostd::substring2bool(vaflowlib.at(i).vfiles,"POSCAR.relax1")) {
            stringstream aus(aurostd::url2string(vaflowlibentry_url.at(i)+"/POSCAR.relax1"));xstr=xstructure(aus,IOVASP_AUTO);}
          plug.vstr.push_back(xstr);xstr.clear(); //DX20191220 - uppercase to lowercase clear
          // MID
          if(!aurostd::substring2bool(vaflowlib.at(i).vfiles,"CONTCAR.relax1") && !aurostd::substring2bool(vaflowlib.at(i).vfiles,"POSCAR.relax2")) {
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"vaflowlibentry_url.at(i)=\""+vaflowlibentry_url.at(i)+"\" needed POSCAR.relax2/CONTCAR.relax1",_FILE_CORRUPT_); //CO20200624
          }
          if(!xstr.atoms.size() && aurostd::substring2bool(vaflowlib.at(i).vfiles,"CONTCAR.relax1")) {
            stringstream aus(aurostd::url2string(vaflowlibentry_url.at(i)+"/CONTCAR.relax1"));xstr=xstructure(aus,IOVASP_AUTO);}
          if(!xstr.atoms.size() && aurostd::substring2bool(vaflowlib.at(i).vfiles,"POSCAR.relax2")) {
            stringstream aus(aurostd::url2string(vaflowlibentry_url.at(i)+"/POSCAR.relax2"));xstr=xstructure(aus,IOVASP_AUTO);}
          plug.vstr.push_back(xstr);xstr.clear(); //DX20191220 - uppercase to lowercase clear
          // POST
          if(!aurostd::substring2bool(vaflowlib.at(i).vfiles,"CONTCAR.relax") && !aurostd::substring2bool(vaflowlib.at(i).vfiles,"CONTCAR.relax2")) {
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"vaflowlibentry_url.at(i)=\""+vaflowlibentry_url.at(i)+"\" needed CONTCAR.relax/CONTCAR.relax2",_FILE_CORRUPT_); //CO20200624
          }
          if(!xstr.atoms.size() && aurostd::substring2bool(vaflowlib.at(i).vfiles,"CONTCAR.relax")) {
            stringstream aus(aurostd::url2string(vaflowlibentry_url.at(i)+"/CONTCAR.relax"));xstr=xstructure(aus,IOVASP_AUTO);}
          if(!xstr.atoms.size() && aurostd::substring2bool(vaflowlib.at(i).vfiles,"CONTCAR.relax2")) {
            stringstream aus(aurostd::url2string(vaflowlibentry_url.at(i)+"/CONTCAR.relax2"));xstr=xstructure(aus,IOVASP_AUTO);}
          plug.vstr.push_back(xstr);xstr.clear(); //DX20191220 - uppercase to lowercase clear

        }


        // also holes must have something
        plug.vNsgroup.clear();plug.vNsgroup.push_back(1);plug.vNsgroup.push_back(1);plug.vNsgroup.push_back(1);   // also holes must have something
        plug.vsgroup.clear();plug.vsgroup.push_back("1");plug.vsgroup.push_back("1");plug.vsgroup.push_back("1");   // also holes must have something
        // SPACEGROUP PROPERTIES

        plug.energy_atom_relax1=plug.energy_atom;  // WRAPPER NEEDS TO BE FIXED

        // tiny corrections
        if(pseudoA=="Cd" && pseudoB=="Pt" && prototype=="181") { //gamma_IrV
          plug.enthalpy_formation_atom-=0.0013;plug.enthalpy_formation_cell=plug.natoms*plug.enthalpy_formation_atom;
        }
        if(pseudoA=="Ir" && pseudoB=="V_sv" && prototype=="291") { //gamma_IrV
          plug.enthalpy_formation_cell-=0.001;plug.enthalpy_formation_atom-=0.0005;plug.enthalpy_cell-=0.001;plug.enthalpy_atom-=0.005;
        }
        if(pseudoA=="Hf_pv" && pseudoB=="Pd_pv" && prototype=="192") { // HfPd
          plug.enthalpy_formation_atom-=0.003;plug.enthalpy_formation_cell=plug.natoms*plug.enthalpy_formation_atom;
          plug.enthalpy_atom=plug.enthalpy_formation_atom;plug.enthalpy_cell=plug.natoms*plug.enthalpy_atom;
        }
        if(pseudoA=="Ir" && pseudoB=="Nb_sv" && prototype=="600.ABBAB") { // sigma
          plug.enthalpy_formation_cell+=0.001;plug.enthalpy_formation_atom+=0.0005;plug.enthalpy_cell+=0.001;plug.enthalpy_atom+=0.005;
          plug.enthalpy_formation_cell+=0.001;plug.enthalpy_formation_atom+=0.0005;plug.enthalpy_cell+=0.001;plug.enthalpy_atom+=0.005;
        }
        if(pseudoA=="Os_pv" && pseudoB=="Re_pv" && prototype=="122") { // sigma
          plug.enthalpy_formation_atom-=0.001;plug.enthalpy_formation_cell=plug.natoms*plug.enthalpy_formation_atom;
          plug.enthalpy_atom=plug.enthalpy_formation_atom;plug.enthalpy_cell=plug.natoms*plug.enthalpy_atom;
        }
        if(pseudoA=="Os_pv" && pseudoB=="Re_pv" && prototype=="448") { // sigma
          ishole=TRUE; // bug
        }
        if(pseudoA=="Rh_pv" && pseudoB=="Zr_sv" && prototype=="381") { // wrong channel
          ishole=TRUE; // bug
        }


        if(LDEBUG) cerr << "plug.entry=" << plug.entry << endl;
        if(LDEBUG) cerr << "prototype=" << prototype << endl;

        struct_proprts.GetSGROUP(plug.entry);
        plug.vNsgroup.clear(); for(uint ii=0;ii<struct_proprts.vNsgroup.size();ii++) plug.vNsgroup.push_back(struct_proprts.vNsgroup.at(ii));
        plug.vsgroup.clear(); for(uint ii=0;ii<struct_proprts.vsgroup.size();ii++) plug.vsgroup.push_back(struct_proprts.vsgroup.at(ii));
        if(LDEBUG) cerr << "plug.vNsgroup.size()=" << plug.vNsgroup.size() << endl;
        if(LDEBUG) cerr << "plug.vsgroup.size()=" << plug.vsgroup.size() << endl;

        plug.structure_name=plug.prototype;
        plug.structure_description=plug.prototype;

        plug.FixDescription();  // fix description
        if(LDEBUG) cerr << "plug.structure_name=" << plug.structure_name << endl;
        if(LDEBUG) cerr << "plug.structure_description=" << plug.structure_description << endl;

        vector<string> tokens;

        aurostd::string2tokens(plug.composition,tokens,",");
        plug.pureA=FALSE;plug.pureB=FALSE;
        if(tokens.size()==2) {
          plug.stoich_a=aurostd::string2utype<double>(tokens.at(0))/((double) plug.natoms);                 // STOICH_A
          plug.stoich_b=aurostd::string2utype<double>(tokens.at(1))/((double) plug.natoms);                 // STOICH_B
        }
        if(tokens.size()==1 && plug.species==speciesA) {plug.stoich_a=1.0;plug.stoich_b=0.0;plug.pureA=TRUE;plug.pureB=FALSE;}
        if(tokens.size()==1 && plug.species==speciesB) {plug.stoich_a=0.0;plug.stoich_b=1.0;plug.pureA=FALSE;plug.pureB=TRUE;}

        // aurostd::string2tokens(plug.geometry,tokens,",");
        // //	cerr << "tokens.size()=" << tokens.size() << endl;
        // if(tokens.size()!=6) {
        //  message << "(get_dimensions_properties) abc-tokens " << endl;message << plug.entry << endl;
        //  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_NUMBER_); //CO20200624
        // }

        // plug.unit_cell_a=(double) aurostd::string2utype<double>(tokens.at(0));         // UCELLD unit cell: a (max)
        // plug.unit_cell_b=(double) aurostd::string2utype<double>(tokens.at(1));         // UCELLD unit cell: b
        // plug.unit_cell_c=(double) aurostd::string2utype<double>(tokens.at(2));         // UCELLD unit cell: b
        // plug.unit_cell_alpha=(double) aurostd::string2utype<double>(tokens.at(3));     // UCELLD unit cell: alpha (b,c)
        // plug.unit_cell_beta=(double) aurostd::string2utype<double>(tokens.at(4));      // UCELLD unit cell: beta  (a,c)
        // plug.unit_cell_gamma=(double) aurostd::string2utype<double>(tokens.at(5));     // UCELLD unit cell: gamma (a,b)

        aurostd::StringSubst(plug.entry,"NNN #0","");
        //  cerr << plug.entry << endl;

        aurostd::string2tokens(plug.nbondxx,tokens,",");
        //	cerr << "plug.nbondxx=" << plug.nbondxx << endl;
        //	cerr << "tokens.size()=" << tokens.size() << endl;
        if(tokens.size()==1) {plug.bond_aa=aurostd::string2utype<double>(tokens.at(0));plug.bond_ab=0.0;plug.bond_bb=0.0;}
        if(tokens.size()==3) {plug.bond_aa=aurostd::string2utype<double>(tokens.at(0));plug.bond_ab=aurostd::string2utype<double>(tokens.at(1));plug.bond_bb=aurostd::string2utype<double>(tokens.at(2));}
        if(tokens.size()==0) {
          message << "(get_dimensions_properties) nbondxx not enough tokens " << endl;message << plug.entry << endl;
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_NUMBER_); //CO20200624
        }
        // done

        if(!ishole && aflags.vflag.flag("APENNSY::NEGLECT_STRUCTURES")) 
          for(uint i=0;i<aflags.APENNSY_NEGLECT_STRUCTURES_vstrs.size()&&!ishole;i++) 
            if(prototype==aflags.APENNSY_NEGLECT_STRUCTURES_vstrs.at(i)) ishole=TRUE;

        if(!ishole) ZLibrary.at(nalloy).push_back(plug);

      } // i as Nstructures

      // fix A6 becoming FCC: A6 => FCC sg(FCC=#225) sg(A6)=#129
      //      for(uint i1=1;i1<structures[2].size();i1++)
      // 	if(structures[2].at(i1).pureA==TRUE && structures[2].at(i1).name=="303") // A6 pure A
      // 	  for(uint i2=1;i2<structures[2].size();i2++)
      // 	    if(structures[2].at(i2).pureA==TRUE && structures[2].at(i2).name=="2") { // FCC pure A
      // 	      if(ZLibrary.at(nalloy).at(i1).vNsgroup.back()==225 && ZLibrary.at(nalloy).at(i1).enthalpy_atom<=ZLibrary.at(nalloy).at(i2).enthalpy_atom) {
      // 		ZLibrary.at(nalloy).at(i1).enthalpy_atom=ZLibrary.at(nalloy).at(i2).enthalpy_atom+0.001999;}} // A6 => FCC sg(FCC=#225) sg(A6)=#129
      //      for(uint i1=1;i1<structures[2].size();i1++)
      // 	if(structures[2].at(i1).pureB==TRUE && structures[2].at(i1).name=="304") // A6 pureB
      // 	  for(uint i2=1;i2<structures[2].size();i2++)
      // 	    if(structures[2].at(i2).pureB==TRUE && structures[2].at(i2).name=="1") { // FCC pureB
      // 	      if(ZLibrary.at(nalloy).at(i1).vNsgroup.back()==225 && ZLibrary.at(nalloy).at(i1).enthalpy_atom<=ZLibrary.at(nalloy).at(i2).enthalpy_atom) {
      // 		ZLibrary.at(nalloy).at(i1).enthalpy_atom=ZLibrary.at(nalloy).at(i2).enthalpy_atom+0.001999;}} // A6 => FCC sg(FCC=#225) sg(A6)=#129

    }
    if(LDEBUG)
      cerr << "ZLibrary.at(nalloy).size()=" << ZLibrary.at(nalloy).size() << endl;

    //   for(uint j=0;j<ZLibrary.at(nalloy).size();j++) cerr << "ZLibrary.at(nalloy).at(" << j << ").vNsgroup.size()=" << ZLibrary.at(nalloy).at(j).vsgroup.at(1) << endl;

    AlloyStructureIdentity.at(nalloy).clear();
    for(uint j=0;j<ZLibrary.at(nalloy).size();j++)
      AlloyStructureIdentity.at(nalloy).push_back(*(new vector<bool>(ZLibrary.at(nalloy).size()+1))); // alloys.size()*ZLibrary.at(nalloy).size()*ZLibrary.at(nalloy).size()

    //   // from FORM ENERGY TO ORDER
    //    // ORDER PER STRUCTURES
    //    double xa=0.0,xb=0.0;
    //    for(i=1;i<=ZLibrary.at(nalloy).size();i++)  { // make formation energy
    //      xa=ZLibrary.at(nalloy).at(i).stoich_a;
    //      xb=ZLibrary.at(nalloy).at(i).stoich_b;
    //      if(xa>0.0 && xb>0.0) {
    //  	ZLibrary.at(nalloy).at(i).order_formation_atom=(ZLibrary.at(nalloy).at(i).enthalpy_formation_atom/(xa*log(xa)+xb*log(xb)))*1000*log(2.0);//eV2K;
    //      } else {
    // 	ZLibrary.at(nalloy).at(i).order_formation_atom=0.0;
    //      }
    //    }
    //    // ORDER PER ALLOY
    //    alloy_order.at(nalloy)=0.0;
    //    for(i=1;i<=ZLibrary.at(nalloy).size();i++)
    //      if(ZLibrary.at(nalloy).at(i).gnd_available==TRUE)
    // 	//  alloy_order.at(nalloy)=(double) max(alloy_order.at(nalloy),ZLibrary.at(nalloy).at(i).order_formation_atom*log(2.0)/1000/KBOLTZEV);
    // 	//  alloy_order.at(nalloy)=(double) max(alloy_order.at(nalloy),ZLibrary.at(nalloy).at(i).order_formation_atom*log(2.0)*1000);
    // 	alloy_order.at(nalloy)=(double) max(alloy_order.at(nalloy),ZLibrary.at(nalloy).at(i).order_formation_atom);

    // determine miscibility from HT
    gndstates=xconvexhull(nalloy,*this);
    for(uint ii=1;ii<=ZConcentrations.at(nalloy).size()-1;ii++) 
      ZGNDlibrary.at(nalloy).at(ii)=gndstates(ii);
    gndstates=xminenergy(nalloy,*this); 
    for(uint ii=1;ii<=ZConcentrations.at(nalloy).size()-1;ii++) 
      ZMINElibrary.at(nalloy).at(ii)=gndstates(ii);

    // calculate max Hf and max Ts
    enthalpy_formation_atom_max.at(nalloy)=0.0;
    entropic_temperature_max.at(nalloy)=0.0;
    for(uint i=0;i<ZLibrary.at(nalloy).size();i++) {
      if(ZLibrary.at(nalloy).at(i).enthalpy_formation_atom<enthalpy_formation_atom_max.at(nalloy))
        enthalpy_formation_atom_max.at(nalloy)=ZLibrary.at(nalloy).at(i).enthalpy_formation_atom;
      if(ZLibrary.at(nalloy).at(i).entropic_temperature>entropic_temperature_max.at(nalloy))
        entropic_temperature_max.at(nalloy)=ZLibrary.at(nalloy).at(i).entropic_temperature;
    }

    // fprintf(stderr," [");for(uint ii=1;ii<=ZConcentrations.at(nalloy).size()-1;ii++) fprintf(stderr," %3i",structures_number[ZGNDlibrary.at(nalloy).at(ii)]); fprintf(stderr," ]  ");
    // fprintf(stderr," [");for(uint ii=1;ii<=ZConcentrations.at(nalloy).size()-1;ii++) fprintf(stderr," %3i",structures_number[MINElibrary[ii].at(nalloy)]); fprintf(stderr," ]  ");
    // Determine MISCIBILITY
    gndstates=xconvexhull(nalloy,*this);
    if(xphaseseparation(gndstates,nalloy,*this)==FALSE) {
      Alloy2MiscibilityHT.at(nalloy)=MISCIBILITY_SYSTEM_MISCIBLE;
    } else {
      if(ZLibrary.at(nalloy).size()>= MISCIBILITY_SYSTEM_CUTOFF) {
        Alloy2MiscibilityHT.at(nalloy)=MISCIBILITY_SYSTEM_NOMIX;
      } else {
        Alloy2MiscibilityHT.at(nalloy)=MISCIBILITY_SYSTEM_UNKNOWN;
      }
    }
    cerr << " ";
    //    // determine chemical potentials
    //    uint point1,point0;
    //    point1=DGNDlibrary.at(nalloy).at(1);point0=DGNDlibrary.at(nalloy).at(0);
    //    alloy_deltamuBArich.at(nalloy)=
    //      (ZLibrary.at(nalloy).at(point1).enthalpy_formation_atom-ZLibrary.at(nalloy).at(point0).enthalpy_formation_atom)/
    //      (structures[2].at(point1).vconcentrations.at(1)-structures[2].at(point0).vconcentrations.at(1));
    //    point1=DGNDlibrary.at(nalloy).at(DGNDlibrary.at(nalloy).size()-2);point0=DGNDlibrary.at(nalloy).at(DGNDlibrary.at(nalloy).size()-1);
    //    alloy_deltamuABrich.at(nalloy)=
    //      (ZLibrary.at(nalloy).at(point1).enthalpy_formation_atom-ZLibrary.at(nalloy).at(point0).enthalpy_formation_atom)/
    //      (structures[2].at(point1).vconcentrations.at(0)-structures[2].at(point0).vconcentrations.at(0));


    // Miscibility write
    cerr << APENNSY_MiscibilityString(Alloy2MiscibilityHT.at(nalloy)) << " ";           // WRITE MISCIBILITYAlloy2MiscibilityHT.at(nalloy)
    // PAULING: determine and write miscibility
    Alloy2MiscibilityExperiments.at(nalloy)=MiscibilityExperimentsCheck(alloy_tmp);
    cerr << APENNSY_MiscibilityString(Alloy2MiscibilityExperiments.at(nalloy)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityExperiments.at(nalloy)
    // MIEDEMA: determine and write miscibility
    Alloy2MiscibilityMiedema.at(nalloy)=MiscibilityMiedemaCheck(alloy_tmp);
    cerr << APENNSY_MiscibilityString(Alloy2MiscibilityMiedema.at(nalloy)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityMiedema.at(nalloy)
    // HUMEROTHERY: determine and write miscibility
    Alloy2MiscibilityHumeRothery.at(nalloy)=MiscibilityHumeRotheryCheck(alloy_tmp);
    cerr << APENNSY_MiscibilityString(Alloy2MiscibilityHumeRothery.at(nalloy)) << " ";  // WRITE MISCIBILITY Alloy2MiscibilityHumeRothery.at(nalloy)

    cerr << "  " << aurostd::PaddedPOST(alloy_tmp,11," ");

    if(xphaseseparation(gndstates,nalloy,*this)==FALSE) {
      gndstates=xminenergy(nalloy,*this);
      //  cout << "XXXXXX " << alloy << "   structure2run=" << structures_number[gndstates(5)] << endl;
    }
    for (uint i=0;i<ZLibrary.at(nalloy).size();i++) ZLibrary.at(nalloy).at(i).distance_gnd=GNDdistance(i,nalloy,*this);     // calcualte distance from CONVEX HULL !
    for (uint i=0;i<ZLibrary.at(nalloy).size();i++) ZLibrary.at(nalloy).at(i).distance_tie=TIELINEdistance(i,nalloy,*this); // calcualte distance from CONVEX HULL without that concentration!

    // PRINT the holes
    //   if(holes>1) cerr << "     holes=" << (holes<100 ? " " : "") << (holes<10 ? " " : "") << holes << " ";
    //   if(holes>1) cerr << "     " << (int) 100*(Nstructures-holes)/Nstructures << "% ";
    //  if(holes>1)
    cerr << "  " << (int) ZLibrary.at(nalloy).size(); // << Hlibrary(nalloy) << " ";
    // plot ORDER TENDENCY
    cerr << "   ";
    cerr.setf(std::ios::fixed,std::ios::floatfield);
    cerr.precision(8);
    cerr << int(alloy_order.at(nalloy));
    // print MISCIBILITY
    //   int miscibility=MiscibilityCheck(string(alloy+30));
    aurostd::string2tokens(alloy,tokens,"/");int miscibility=MiscibilityCheck(tokens.at(tokens.size()-1));
    if(miscibility==MISCIBILITY_SYSTEM_UNKNOWN) cerr << "    MISCIBILITY_SYSTEM_UNKNOWN";
    for(uint i=0;i<ZGNDlibrary.at(nalloy).size();i++) if(ZGNDlibrary.at(nalloy).at(i)>=0) cerr << "    " << ZLibrary.at(nalloy).at(ZGNDlibrary.at(nalloy).at(i)).structure_name;

    //    cerr << " " << alloy_deltamuBArich.at(nalloy) << " " << alloy_deltamuABrich.at(nalloy) << " "; // UNCALCULATED
    cerr << endl;
  } // nalloy<alloys.size()
  cerr << "Total structures calculated = " << LIBstructures_calcs  << endl;
  // FINISHED
  return TRUE;
}

// ***************************************************************************
// APENNSY_MiscibilityString
// ***************************************************************************
string APENNSY_MiscibilityString(int flag) {
  if(flag==MISCIBILITY_SYSTEM_SOLUTION) return "M";
  if(flag==MISCIBILITY_SYSTEM_MISCIBLE) return "C";
  if(flag==MISCIBILITY_SYSTEM_NOMIX) return "S";
  if(flag==MISCIBILITY_SYSTEM_UNKNOWN) return "?";
  return  "?";
}

// **************************************************************************
// **************************************************************************
// **************************************************************************
// CONVEX HULL FUNCTIONS

xvector<int> xconvexhull(int nalloy,const APENNSY_Parameters &params) {
  // given colums of alloys and concentration,
  // gives formation energy and convexhull
  // cerr << "params.ZConcentrations.at(nalloy).size()-1=" << params.ZConcentrations.at(nalloy).size()-1 << endl;
  // cerr << "params.ZLibrary.at(nalloy).size()=" << params.ZLibrary.at(nalloy).size() << endl;
  xvector<double> gndEFlibrary(params.ZConcentrations.at(nalloy).size()-1);
  xvector<int> gndEFindex(params.ZConcentrations.at(nalloy).size()-1);
  double c1,c2,c3,alpha,beta,energy;
  //  for(uint i=0;i<params.ZConcentrations.at(nalloy).size();i++) 
  //  cerr << "i=" << i << " params.ZConcentrations.at(nalloy).at(i)=" << params.ZConcentrations.at(nalloy).at(i) << endl;

  for(uint i=1;i<params.ZConcentrations.at(nalloy).size();i++) {
    gndEFlibrary(i)=APENNSY_INF;
    gndEFindex(i)=-1;  // default
    for (uint j=0;j<params.ZLibrary.at(nalloy).size();j++) {
      //cout << "j=" << j << " " << params.ZLibrary.at(nalloy).at(j).structure_name << " " << params.ZLibrary.at(nalloy).at(j).stoich_b << " " << endl;
      if(aurostd::isequal(params.ZLibrary.at(nalloy).at(j).stoich_b,params.ZConcentrations.at(nalloy).at(i),CEPSILON) && params.ZLibrary.at(nalloy).at(j).enthalpy_formation_atom<gndEFlibrary(i)) {
        gndEFlibrary(i)=params.ZLibrary.at(nalloy).at(j).enthalpy_formation_atom;
        gndEFindex(i)=j;
      }
    }
    //    cout << "i=" << i <<  "   gndEFindex(i)=" << gndEFindex(i) << " " << params.ZLibrary.at(nalloy).at(gndEFindex(i)).structure_name << endl;
  }
  // for(uint i=1;i<=params.ZConcentrations.at(nalloy).size()-1;i++)
  //  if(gndEFindex(i)>0) cerr << params.ZLibrary.at(nalloy).at(gndEFindex(i)).structure_name << endl;

  for(uint j=1;j<params.ZConcentrations.at(nalloy).size();j++) {
    for(uint k=1;k<params.ZConcentrations.at(nalloy).size();k++) {
      for(uint i=1;i<params.ZConcentrations.at(nalloy).size();i++) {
        if(gndEFindex(i)>=0 && gndEFindex(j)>=0 && gndEFindex(k)>=0 && j<i && i<k) {
          c1=params.ZConcentrations.at(nalloy).at(j);c2=params.ZConcentrations.at(nalloy).at(k);c3=params.ZConcentrations.at(nalloy).at(i);alpha=(c3-c1)/(c2-c1);beta=(c2-c3)/(c2-c1);
          energy=beta*params.ZLibrary.at(nalloy).at(gndEFindex(j)).enthalpy_formation_atom+alpha*params.ZLibrary.at(nalloy).at(gndEFindex(k)).enthalpy_formation_atom;
          energy=beta*params.ZLibrary.at(nalloy).at(gndEFindex(j)).enthalpy_formation_atom+alpha*params.ZLibrary.at(nalloy).at(gndEFindex(k)).enthalpy_formation_atom;
          energy=beta*params.ZLibrary.at(nalloy).at(gndEFindex(j)).enthalpy_formation_atom+alpha*params.ZLibrary.at(nalloy).at(gndEFindex(k)).enthalpy_formation_atom;
          if(energy<params.ZLibrary.at(nalloy).at(gndEFindex(i)).enthalpy_formation_atom) gndEFindex(i)=-1;
        }
      }
    }
  }
  //  for(uint i=1;i<=params.ZConcentrations.at(nalloy).size()-1;i++)
  //   if(gndEFindex(i)>=0) cerr << params.ZLibrary.at(nalloy).at(gndEFindex(i)).structure_name << endl;

  //   cerr << gndEFindex << endl;
  return gndEFindex;
}

// **************************************************************************
// **************************************************************************
// **************************************************************************
// MIMIMUN ENERGY FUNCTIONS

xvector<int> xminenergy(int nalloy,const APENNSY_Parameters &params) {
  // given colums of alloys and concentration,
  // gives formation energy and minenergy
  uint i,j;
  xvector<double> gndEFlibrary(params.ZConcentrations.at(nalloy).size()-1);
  xvector<int> gndEFindex(params.ZConcentrations.at(nalloy).size()-1);
  for(i=1;i<=params.ZConcentrations.at(nalloy).size()-1;i++) {
    gndEFlibrary(i)=APENNSY_INF;
    gndEFindex(i)=0;
    for (j=0;j<params.ZLibrary.at(nalloy).size();j++) {
      if(aurostd::isequal(params.ZLibrary.at(nalloy).at(j).stoich_b,params.ZConcentrations.at(nalloy).at(i),CEPSILON) && params.ZLibrary.at(nalloy).at(j).enthalpy_formation_atom<gndEFlibrary(i)) {
        gndEFlibrary(i)=params.ZLibrary.at(nalloy).at(j).enthalpy_formation_atom;
        gndEFindex(i)=j;
      }
    }
  }
  return gndEFindex;
}

// **************************************************************************
// **************************************************************************
// **************************************************************************
// PHASE SEPARATION

bool xphaseseparation(const xvector<int>& converhull,const int& nalloy,const APENNSY_Parameters &params) {
  bool out=TRUE;
  for (uint i=2;i<=params.ZConcentrations.at(nalloy).size()-1-1;i++)
    if(converhull(i)>0) out=FALSE;
  return out;
}

// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************
// DISTANCES STUFF

double GNDdistance(const int& i, const int& nalloy,const APENNSY_Parameters &params) {  // i = structure nalloy=alloy
  // cerr << "GNDdistance BEGIN" << endl;
  double distance=0.0,Eaa,Ebb,Ecc,caa,cbb,ccc;
  int kk_mid=0,kk_low=0,kk_high=0;
  bool done=FALSE;
  for(uint kk=1;kk<=params.ZConcentrations.at(nalloy).size()-1;kk++)
    if(aurostd::isequal(params.ZConcentrations.at(nalloy).at(kk),params.ZLibrary.at(nalloy).at(i).stoich_b,CEPSILON)) {
      kk_mid=kk;
      if(params.ZGNDlibrary.at(nalloy).at(kk)>=0) {
        distance=params.ZLibrary.at(nalloy).at(i).enthalpy_formation_atom-params.ZLibrary.at(nalloy).at(params.ZGNDlibrary.at(nalloy).at(kk)).enthalpy_formation_atom;
        done=TRUE;
        // questo se il punto e` un GNDSTATE
      }
    }
  // se non e` un GNDSTATE, allora calcola distanza
  if(done==TRUE) return distance;
  for(int kk=1;kk<kk_mid;kk++) if(params.ZGNDlibrary.at(nalloy).at(kk)>=0) kk_low=kk;
  for(int kk=params.ZConcentrations.at(nalloy).size()-1;kk>kk_mid;kk--) if(params.ZGNDlibrary.at(nalloy).at(kk)>=0) kk_high=kk;
  //for(int kk=1;kk<=params.ZConcentrations.at(nalloy).size()-1;kk++) cout << params.GNDlibrary(kk,nalloy) << " ";cout << endl;
  //for(int kk=1;kk<=params.ZConcentrations.at(nalloy).size()-1;kk++) cout << params.MINElibrary(kk,nalloy) << " ";cout << endl;
  Eaa=params.ZLibrary.at(nalloy).at(params.ZGNDlibrary.at(nalloy).at(kk_low)).enthalpy_formation_atom;
  caa=params.ZLibrary.at(nalloy).at(params.ZGNDlibrary.at(nalloy).at(kk_low)).stoich_b;
  Ebb=params.ZLibrary.at(nalloy).at(i).enthalpy_formation_atom;
  cbb=params.ZLibrary.at(nalloy).at(i).stoich_b;
  Ecc=params.ZLibrary.at(nalloy).at(params.ZGNDlibrary.at(nalloy).at(kk_high)).enthalpy_formation_atom;
  ccc=params.ZLibrary.at(nalloy).at(params.ZGNDlibrary.at(nalloy).at(kk_high)).stoich_b;

  distance=Ebb-(Ecc+(Eaa-Ecc)*(cbb-ccc)/(caa-ccc));
  //cout << i << " " << nalloy << " " << params.ZLibrary.at(nalloy).at(i).stoich_b << " " << kk_low << " " << kk_high << " " << Eaa << " " << Ebb << "  " << Ecc << " " << caa << " " << cbb << "  " << ccc << "  " << distance << endl;
  // cerr << "GNDdistance END" << endl;

  return distance;
}

xvector<double> GNDdistance(const int& nalloy,const APENNSY_Parameters &params) {  // nalloy=alloy
  xvector<double> distance((params.ZLibrary.at(nalloy).size()));
  for (uint i=1;i<=params.ZLibrary.at(nalloy).size();i++)
    distance(i)=GNDdistance(i,nalloy,params);
  return distance;
}

double TIELINEdistance(const int& i, const int& nalloy,const APENNSY_Parameters &params) {  // i = structure nalloy=alloy
  double distance=0.0,Eaa,Ebb,Ecc,caa,cbb,ccc;
  int kk_mid=0,kk_low,kk_high;
  kk_mid=1;kk_low=1;kk_high=params.ZConcentrations.at(nalloy).size()-1;
  for(uint kk=1;kk<=params.ZConcentrations.at(nalloy).size()-1;kk++)
    if(aurostd::isequal(params.ZConcentrations.at(nalloy).at(kk),params.ZLibrary.at(nalloy).at(i).stoich_b,CEPSILON)) kk_mid=kk;    // GET THE MID
  for(int kk=1;kk<kk_mid;kk++) if(params.ZGNDlibrary.at(nalloy).at(kk)>=0) kk_low=kk;                              // GET THE LOW
  for(int kk=params.ZConcentrations.at(nalloy).size()-1;kk>kk_mid;kk--) if(params.ZGNDlibrary.at(nalloy).at(kk)>=0) kk_high=kk;               // GET THE HIGH
  //for(int kk=1;kk<=params.ZConcentrations.at(nalloy).size()-1;kk++) cout << params.GNDlibrary(kk,nalloy) << " ";cout << endl;
  //for(int kk=1;kk<=params.ZConcentrations.at(nalloy).size()-1;kk++) cout << params.MINElibrary(kk,nalloy) << " ";cout << endl;
  Eaa=params.ZLibrary.at(nalloy).at(params.ZGNDlibrary.at(nalloy).at(kk_low)).enthalpy_formation_atom;
  caa=params.ZLibrary.at(nalloy).at(params.ZGNDlibrary.at(nalloy).at(kk_low)).stoich_b;
  Ebb=params.ZLibrary.at(nalloy).at(i).enthalpy_formation_atom;
  cbb=params.ZLibrary.at(nalloy).at(i).stoich_b;
  Ecc=params.ZLibrary.at(nalloy).at(params.ZGNDlibrary.at(nalloy).at(kk_high)).enthalpy_formation_atom;
  ccc=params.ZLibrary.at(nalloy).at(params.ZGNDlibrary.at(nalloy).at(kk_high)).stoich_b;
  distance=Ebb-(Ecc+(Eaa-Ecc)*(cbb-ccc)/(caa-ccc));
  //cout << i << " " << nalloy << " " << params.ZLibrary.at(nalloy).at(i).stoich_b << " " << kk_low << " " << kk_high << " " << Eaa << " " << Ebb << "  " << Ecc << " " << caa << " " << cbb << "  " << ccc << "  " << distance << endl;
  return distance;
}

xvector<double> TIELINEdistance(const int& nalloy,const APENNSY_Parameters &params) {  // nalloy=alloy
  xvector<double> distance((params.ZLibrary.at(nalloy).size()));
  for (uint i=1;i<=(params.ZLibrary.at(nalloy).size());i++)
    distance(i)=TIELINEdistance(i,nalloy,params);
  return distance;
}

double MINEdistance(const int& i, const int& nalloy,const APENNSY_Parameters &params) {  // i = structure nalloy=alloy
  double distance=0.0;
  for (uint kk=1;kk<=params.ZConcentrations.at(nalloy).size()-1;kk++) {
    if(aurostd::isequal(params.ZConcentrations.at(nalloy).at(kk),params.ZLibrary.at(nalloy).at(i).stoich_b,CEPSILON))
      distance=params.ZLibrary.at(nalloy).at(i).enthalpy_formation_atom-params.ZLibrary.at(nalloy).at(params.ZMINElibrary.at(nalloy).at(kk)).enthalpy_formation_atom;
  }
  return distance;
}

xvector<double> MINEdistance(const int& nalloy,const APENNSY_Parameters &params) {  // nalloy=alloy
  xvector<double> distance((params.ZLibrary.at(nalloy).size()));
  for (uint i=1;i<=(params.ZLibrary.at(nalloy).size());i++)
    distance(i)=MINEdistance(i,nalloy,params);
  return distance;
}

// **************************************************************************
// MakeRankLib - these are the functions making the Ranking (groundstates)
// **************************************************************************
void APENNSY_Parameters::MakeRankLib(void) {
  vector<double> ausE;
  vector<int> ausI;
  bool VERBOSE=FALSE,VERBOSE1=FALSE;
  // cerr << "MakeRankLib IN" << endl;
  for(uint nalloy=0;nalloy<alloys.size();nalloy++) {  // nalloy alloys
    for(uint j=1;j<ZConcentrations.at(nalloy).size();j++) {
      RankLib.at(nalloy).at(j).clear();
      for(uint k=0;k<ZLibrary.at(nalloy).size();k++)
        if(aurostd::isequal(ZLibrary.at(nalloy).at(k).stoich_b,ZConcentrations.at(nalloy).at(j),CEPSILON))
          RankLib.at(nalloy).at(j).push_back(k);
    }
    if(VERBOSE) {
      cout << "****************************************************************************" << endl;
      cout << "RANKING Before/AFTER     alloy=" << nalloy << endl;
    }
    for(uint j=1;j<ZConcentrations.at(nalloy).size();j++) {
      // cerr << "MakeRankLib MED ZConcentrations.at(nalloy).size()=" << ZConcentrations.at(nalloy).size()<< endl;
      if(ZConcentrations.at(nalloy).at(j)>-CEPSILON && ZConcentrations.at(nalloy).at(j)<1.0+CEPSILON) { // [0,1]
        // NOW PRINT NUMBERS
        if(VERBOSE) {
          printf("  %3.3f ",ZConcentrations.at(nalloy).at(j));
          for(uint k=0;k<RankLib.at(nalloy).at(j).size();k++) printf("%6i ",RankLib.at(nalloy).at(j).at(k));
          cout << endl;
        }
        // NOW PRINT ENERGIES
        if(VERBOSE1) {
          printf("  %3.3f ",ZConcentrations.at(nalloy).at(j));
          for(uint k=0;k<RankLib.at(nalloy).at(j).size();k++) printf("%+3.3f ",ZLibrary.at(nalloy).at(RankLib.at(nalloy).at(j).at(k)).enthalpy_formation_atom);
          cout << endl;
        }
        // NOW ORDER STUFF !

        // UP TO HERE
        ausE.clear();ausI.clear();
        for(uint k=0;k<RankLib.at(nalloy).at(j).size();k++) {
          ausE.push_back(ZLibrary.at(nalloy).at(RankLib.at(nalloy).at(j).at(k)).enthalpy_formation_atom);
          ausI.push_back(RankLib.at(nalloy).at(j).at(k));
        }
        aurostd::sort(ausE,ausI);
        for(uint k=0;k<RankLib.at(nalloy).at(j).size();k++)
          RankLib.at(nalloy).at(j).at(k)=ausI.at(k);

        // NOW PRINT NUMBERS
        if(VERBOSE) {
          printf("  %3.3f ",ZConcentrations.at(nalloy).at(j));
          for(uint k=0;k<RankLib.at(nalloy).at(j).size();k++) printf("%6i ",RankLib.at(nalloy).at(j).at(k));
          cout << endl;
        }
        // NOW PRINT ENERGIES
        if(VERBOSE1) {
          printf("  %3.3f ",ZConcentrations.at(nalloy).at(j));
          for(uint k=0;k<RankLib.at(nalloy).at(j).size();k++) printf("%+3.3f ",ZLibrary.at(nalloy).at(RankLib.at(nalloy).at(j).at(k)).enthalpy_formation_atom);
          cout << endl;
        }
      }
    }
  }  // RankLib has been Built !
  // fix reference noise
  for(uint k=0;k<alloys.size();k++) { // if(Alloy2MiscibilityHT.at(k)!=MISCIBILITY_SYSTEM_MISCIBLE)
    for(uint j=0;j<ZConcentrations.at(k).size();j++) {
      for(int i=RankLib.at(k).at(j).size()-1;i>=0;i--) {
        if(ZConcentrations.at(k).at(j)<CEPSILON || ZConcentrations.at(k).at(j)>1-CEPSILON) {
          ZLibrary.at(k).at(RankLib.at(k).at(j).at((uint) i)).enthalpy_formation_atom -= ZLibrary.at(k).at(RankLib.at(k).at(j).at(0)).enthalpy_formation_atom;
          // cerr << "MakeRankLib:" << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(RankLib.at(k).at(j).at((uint) i)).enthalpy_formation_atom << " " << (uint) i << endl;
        }
      }
    }
  }
  // cerr << "MakeRankLib OUT" << endl;
}

// **************************************************************************
// CheckAlloyRelaxationStructure - Check if structures have relaxed
// **************************************************************************
void APENNSY_Parameters::CheckAlloyRelaxationStructures(void) {   // check if 2 structures are the same
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  if(LDEBUG) cerr << "CheckAlloyRelaxationStructures IN" << endl;
  double ausERR,normAA1,normAB1,normBB1,normAA2,normAB2,normBB2;
  for(uint nalloy=0;nalloy<alloys.size();nalloy++)
    for(uint j1=0;j1<ZLibrary.at(nalloy).size();j1++)
      for(uint j2=0;j2<ZLibrary.at(nalloy).size();j2++)
        AlloyStructureIdentity.at(nalloy).at(j1).at(j2)=FALSE;
  if(LDEBUG) cerr << "CheckAlloyRelaxationStructures MED" << endl;
  for(uint nalloy=0;nalloy<alloys.size();nalloy++) {
    for(uint j1=0;j1<ZLibrary.at(nalloy).size();j1++) {
      AlloyStructureIdentity.at(nalloy).at(j1).at(j1)=TRUE;
      if(LDEBUG) cerr << "CheckAlloyRelaxationStructures MED nalloy=" << nalloy << " j1=" << j1 << endl;
      if(LDEBUG) cerr << "ZLibrary.at(nalloy).at(j1).structure_name=" << ZLibrary.at(nalloy).at(j1).structure_name << endl;
      if(LDEBUG) cerr << "ZLibrary.at(nalloy).at(j1).volume_cell=" << ZLibrary.at(nalloy).at(j1).volume_cell << endl;
      if(LDEBUG) cerr << "ZLibrary.at(nalloy).at(j1).vgeometry.at(0)=" << ZLibrary.at(nalloy).at(j1).vgeometry.at(0) << endl;
      if(LDEBUG) cerr << "ZLibrary.at(nalloy).at(j1).vgeometry.at(1)=" << ZLibrary.at(nalloy).at(j1).vgeometry.at(1) << endl;
      normAA1=std::pow((double) ZLibrary.at(nalloy).at(j1).volume_cell,(double) 1.0/3.0);
      normAB1=std::pow((double) ZLibrary.at(nalloy).at(j1).vgeometry.at(0),(double) 1.0/3.0);
      normBB1=std::pow((double) ZLibrary.at(nalloy).at(j1).vgeometry.at(1),(double) 1.0/3.0);
      if(LDEBUG) cerr << ZLibrary.at(nalloy).at(j1).vNsgroup.back() << endl;
      for(uint j2=j1;j2<ZLibrary.at(nalloy).size();j2++) {
        if(LDEBUG) cerr << "CheckAlloyRelaxationStructures MED nalloy=" << nalloy << " j1=" << j1  << " j2=" << j2 << endl;
        if(LDEBUG) cerr << "ZLibrary.at(nalloy).at(j1).structure_name=" << ZLibrary.at(nalloy).at(j1).structure_name << endl;
        if(LDEBUG) cerr << "ZLibrary.at(nalloy).at(j1).volume_cell=" << ZLibrary.at(nalloy).at(j1).volume_cell << endl;
        if(LDEBUG) cerr << "ZLibrary.at(nalloy).at(j1).vgeometry.at(0)=" << ZLibrary.at(nalloy).at(j1).vgeometry.at(0) << endl;
        if(LDEBUG) cerr << "ZLibrary.at(nalloy).at(j1).vgeometry.at(1)=" << ZLibrary.at(nalloy).at(j1).vgeometry.at(1) << endl;
        if(LDEBUG) cerr << "ZLibrary.at(nalloy).at(j1).vNsgroup.back()=" << ZLibrary.at(nalloy).at(j1).vNsgroup.back() << endl;
        if(LDEBUG) cerr << "ZLibrary.at(nalloy).at(j2).structure_name=" << ZLibrary.at(nalloy).at(j2).structure_name << endl;
        if(LDEBUG) cerr << "ZLibrary.at(nalloy).at(j2).volume_cell=" << ZLibrary.at(nalloy).at(j2).volume_cell << endl;
        if(LDEBUG) cerr << "ZLibrary.at(nalloy).at(j2).vgeometry.at(0)=" << ZLibrary.at(nalloy).at(j2).vgeometry.at(0) << endl;
        if(LDEBUG) cerr << "ZLibrary.at(nalloy).at(j2).vgeometry.at(1)=" << ZLibrary.at(nalloy).at(j2).vgeometry.at(1) << endl;
        if(LDEBUG) cerr << "ZLibrary.at(nalloy).at(j2).vNsgroup.back()=" << ZLibrary.at(nalloy).at(j2).vNsgroup.back() << endl;
        normAA2=std::pow((double) ZLibrary.at(nalloy).at(j2).volume_cell,(double) 1.0/3.0);  
        normAB2=std::pow((double) ZLibrary.at(nalloy).at(j2).vgeometry.at(0),(double) 1.0/3.0);
        normBB2=std::pow((double) ZLibrary.at(nalloy).at(j2).vgeometry.at(1),(double) 1.0/3.0);
        if(aurostd::isequal(ZLibrary.at(nalloy).at(j1).stoich_b,ZLibrary.at(nalloy).at(j2).stoich_b,CEPSILON) && ZLibrary.at(nalloy).at(j1).vNsgroup.back()==ZLibrary.at(nalloy).at(j2).vNsgroup.back()) {
          ausERR=abs(ZLibrary.at(nalloy).at(j1).enthalpy_formation_atom-ZLibrary.at(nalloy).at(j2).enthalpy_formation_atom);
          if(ausERR<dENERGY_CUTOFF) {
            ausERR=abs(2*(ZLibrary.at(nalloy).at(j1).volume_cell-ZLibrary.at(nalloy).at(j2).volume_cell)/(ZLibrary.at(nalloy).at(j1).volume_cell+ZLibrary.at(nalloy).at(j2).volume_cell)); // volume
            if(ausERR<dVOLUME_CUTOFF) {
              if(2*abs((ZLibrary.at(nalloy).at(j1).bond_aa*normAA1-ZLibrary.at(nalloy).at(j2).bond_aa*normAA2)/(ZLibrary.at(nalloy).at(j1).bond_aa*normAA1+ZLibrary.at(nalloy).at(j2).bond_aa*normAA2)) < dBONDSD_CUTOFF &&
                  2*abs((ZLibrary.at(nalloy).at(j1).bond_ab*normAB1-ZLibrary.at(nalloy).at(j2).bond_ab*normAB2)/(ZLibrary.at(nalloy).at(j1).bond_ab*normAB1+ZLibrary.at(nalloy).at(j2).bond_ab*normAB2)) < dBONDSD_CUTOFF &&
                  2*abs((ZLibrary.at(nalloy).at(j1).bond_bb*normBB1-ZLibrary.at(nalloy).at(j2).bond_bb*normBB2)/(ZLibrary.at(nalloy).at(j1).bond_bb*normBB1+ZLibrary.at(nalloy).at(j2).bond_bb*normBB2)) < dBONDSD_CUTOFF) {
                //   cerr << ausERR << endl;
                AlloyStructureIdentity.at(nalloy).at(j1).at(j2)=TRUE;
                AlloyStructureIdentity.at(nalloy).at(j2).at(j1)=TRUE;
              }
            }
          }
        }
      }
    }
  }
  if(LDEBUG) cerr << "CheckAlloyRelaxationStructures OUT" << endl;
}

// **************************************************************************
// **************************************************************************
// **************************************************************************

// **************************************************************************
// LIBRARY X
bool APENNSY_Parameters::LibLoadAlloysLIB2(_aflags &aflags) {
  // ------------------------------------------------------------------------
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag"))
    cerr << "APENNSY_Parameters::LibLoadAlloysLIB2(_aflags &aflags)" << endl;
  // ------------------------------------------------------------------------
  // LOAD
  alloys.clear();alloysRAW.clear();
  string command;
  vector<string> list;
  if(XHOST.APENNSY_USE_SERVER) { 
    bool TEST=FALSE;//TRUE;
    if(!TEST) command="ls "+init::AFLOW_Projects_Directories("LIB2")+"/RAW/ "; // NORMAL
    else command="ls "+init::AFLOW_Projects_Directories("LIB2")+"/RAW/ | grep -v Cs|  grep -E 'Os|Ru' | grep -v -E 'Al|As|B_h|Ba|Be|Bi|Br|Ca|Cl|Cs|Ga|Ge|In|K_|Li|Mg|Sb|Se|Si|Sn|Sr'"; // DEBUG
    //  cerr << aurostd::execute2string(command) << endl;
    aurostd::string2tokens(aurostd::execute2string(command),list,"\n");
  }
  if(XHOST.APENNSY_SERVER_AFLOWLIB_ORG) {
    // [OBSOLETE]   aurostd::url2tokens(AFLOWLIB_PROJECT_GNDSTATE+string("/"+_XENTRY_+"?aflowlib_entries"),list,",");
    aurostd::url2tokens(AFLOWLIB_PROJECT_GNDSTATE+string("/?aflowlib_entries"),list,",");
  }

  for(uint i=0;i<list.size();i++) {
    bool ok=TRUE;
    if(ok) ok=ok && !aurostd::substring2bool(list.at(i),"core");
    if(ok) ok=ok && !aurostd::substring2bool(list.at(i),"bin");
    if(ok) ok=ok && !aurostd::substring2bool(list.at(i),"nohup");
    if(ok) ok=ok && !aurostd::substring2bool(list.at(i),"x");
    if(ok) ok=ok && !aurostd::substring2bool(list.at(i),"Cs_svPd_pv");
    if(ok) ok=ok && !aurostd::substring2bool(list.at(i),"Ce");
    if(ok) ok=ok && !aurostd::substring2bool(list.at(i),"Gd");
    if(ok) ok=ok && !aurostd::substring2bool(list.at(i),"AlMg");
    if(ok) ok=ok && !aurostd::substring2bool(list.at(i),"GaMg");
    if(ok) ok=ok && !aurostd::substring2bool(list.at(i),"GeMg");
    if(ok) ok=ok && !aurostd::substring2bool(list.at(i),"Mg_pvSi");
    if(ok) alloys.push_back(list.at(i));
    //    if(ok) cerr << list.at(i) << endl;
  }
  aurostd::sort(alloys);
  alloysmesg.clear();alloysmesg.str(std::string());alloysmesg << "Loading LIBRARY X: Nalloys=" << alloys.size() << endl;
  // DONE
  if(XHOST.APENNSY_USE_SERVER) {  
    for(uint i=0;i<alloys.size();i++) alloysRAW.push_back(init::AFLOW_Projects_Directories("LIB2")+"/RAW/"+alloys.at(i));
  }
  if(XHOST.APENNSY_SERVER_AFLOWLIB_ORG) {
    for(uint i=0;i<alloys.size();i++) alloysRAW.push_back(AFLOWLIB_PROJECT_GNDSTATE+"/"+alloys.at(i));
  }
  // for(uint i=0;i<alloys.size();i++) cerr << i << " " << alloys.at(i) << " " << alloysRAW.at(i) << endl;
  LIB_MODE=LIB_MODEX;
  return TRUE;
}

// **************************************************************************
// LIBRARY U
bool APENNSY_Parameters::LibLoadAlloysLIB2U(_aflags &aflags) {
  // ------------------------------------------------------------------------
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag"))
    cerr << "APENNSY_Parameters::LibLoadAlloysLIB2U(_aflags &aflags)" << endl;
  // ------------------------------------------------------------------------
  // LOAD
  alloys.clear();alloysRAW.clear();
  string command;
  vector<string> list;
  command="AgAu,AgCd,AgCo,AgCr_pv,AgCu_pv,AgFe_pv,AgHf_pv,AgHg,AgIr,AgLa,AgMn_pv,AgMo_pv,AgNb_sv,AgNi_pv,AgOs_pv,AgPd_pv,AgPt,AgRe_pv,AgRh_pv,AgRu_pv,AgSc_sv,AgTa_pv,AgTc_pv,AgTi_sv,AgV_sv,AgW_pv,AgY_sv,AgZn,AgZr_sv,AuCd,AuCo,AuCr_pv,AuCu_pv,AuFe_pv,AuHf_pv,AuHg,AuIr,AuLa,AuMn_pv,AuMo_pv,AuNb_sv,AuNi_pv,AuOs_pv,AuPd_pv,AuPt,AuRe_pv,AuRh_pv,AuRu_pv,AuSc_sv,AuTa_pv,AuTc_pv,AuTi_sv,AuV_sv,AuW_pv,AuY_sv,AuZn,AuZr_sv,CdCo,CdCr_pv,CdCu_pv,CdFe_pv,CdHf_pv,CdHg,CdIr,CdLa,CdMn_pv,CdMo_pv,CdNb_sv,CdNi_pv,CdOs_pv,CdPd_pv,CdPt,CdRe_pv,CdRh_pv,CdRu_pv,CdSc_sv,CdTa_pv,CdTc_pv,CdTi_sv,CdV_sv,CdW_pv,CdY_sv,CdZn,CdZr_sv,CoCr_pv,CoCu_pv,CoFe_pv,CoHf_pv,CoHg,CoIr,CoLa,CoMn_pv,CoMo_pv,CoNb_sv,CoNi_pv,CoOs_pv,CoPd_pv,CoPt,CoRe_pv,CoRh_pv,CoRu_pv,CoSc_sv,CoTa_pv,CoTc_pv,CoTi_sv,CoV_sv,CoW_pv,CoY_sv,CoZn,CoZr_sv,Cr_pvCu_pv,Cr_pvFe_pv,Cr_pvHf_pv,Cr_pvHg,Cr_pvIr,Cr_pvLa,Cr_pvMn_pv,Cr_pvMo_pv,Cr_pvNb_sv,Cr_pvNi_pv,Cr_pvOs_pv,Cr_pvPd_pv,Cr_pvPt,Cr_pvRe_pv,Cr_pvRh_pv,Cr_pvRu_pv,Cr_pvSc_sv,Cr_pvTa_pv,Cr_pvTc_pv,Cr_pvTi_sv,Cr_pvV_sv,Cr_pvW_pv,Cr_pvY_sv,Cr_pvZn,Cr_pvZr_sv,Cu_pvFe_pv,Cu_pvHf_pv,Cu_pvHg,Cu_pvIr,Cu_pvLa,Cu_pvMn_pv,Cu_pvMo_pv,Cu_pvNb_sv,Cu_pvNi_pv,Cu_pvOs_pv,Cu_pvPd_pv,Cu_pvPt,Cu_pvRe_pv,Cu_pvRh_pv,Cu_pvRu_pv,Cu_pvSc_sv,Cu_pvTa_pv,Cu_pvTc_pv,Cu_pvTi_sv,Cu_pvV_sv,Cu_pvW_pv,Cu_pvY_sv,Cu_pvZn,Cu_pvZr_sv,Fe_pvHf_pv,Fe_pvHg,Fe_pvIr,Fe_pvLa,Fe_pvMn_pv,Fe_pvMo_pv,Fe_pvNb_sv,Fe_pvNi_pv,Fe_pvOs_pv,Fe_pvPd_pv,Fe_pvPt,Fe_pvRe_pv,Fe_pvRh_pv,Fe_pvRu_pv,Fe_pvSc_sv,Fe_pvTa_pv,Fe_pvTc_pv,Fe_pvTi_sv,Fe_pvV_sv,Fe_pvW_pv,Fe_pvY_sv,Fe_pvZn,Fe_pvZr_sv,Hf_pvHg,Hf_pvIr,Hf_pvLa,Hf_pvMn_pv,Hf_pvMo_pv,Hf_pvNb_sv,Hf_pvNi_pv,Hf_pvOs_pv,Hf_pvPd_pv,Hf_pvPt,Hf_pvRe_pv,Hf_pvRh_pv,Hf_pvRu_pv,Hf_pvSc_sv,Hf_pvTa_pv,Hf_pvTc_pv,Hf_pvTi_sv,Hf_pvV_sv,Hf_pvW_pv,Hf_pvY_sv,Hf_pvZn,Hf_pvZr_sv,HgIr,HgLa,HgMn_pv,HgMo_pv,HgNb_sv,HgNi_pv,HgOs_pv,HgPd_pv,HgPt,HgRe_pv,HgRh_pv,HgRu_pv,HgSc_sv,HgTa_pv,HgTc_pv,HgTi_sv,HgV_sv,HgW_pv,HgY_sv,HgZn,HgZr_sv,IrLa,IrMn_pv,IrMo_pv,IrNb_sv,IrNi_pv,IrOs_pv,IrPd_pv,IrPt,IrRe_pv,IrRh_pv,IrRu_pv,IrSc_sv,IrTa_pv,IrTc_pv,IrTi_sv,IrV_sv,IrW_pv,IrY_sv,IrZn,IrZr_sv,LaMn_pv,LaMo_pv,LaNb_sv,LaNi_pv,LaOs_pv,LaPd_pv,LaPt,LaRe_pv,LaRh_pv,LaRu_pv,LaSc_sv,LaTa_pv,LaTc_pv,LaTi_sv,LaV_sv,LaW_pv,LaY_sv,LaZn,LaZr_sv,Mn_pvMo_pv,Mn_pvNb_sv,Mn_pvNi_pv,Mn_pvOs_pv,Mn_pvPd_pv,Mn_pvPt,Mn_pvRe_pv,Mn_pvRh_pv,Mn_pvRu_pv,Mn_pvSc_sv,Mn_pvTa_pv,Mn_pvTc_pv,Mn_pvTi_sv,Mn_pvV_sv,Mn_pvW_pv,Mn_pvY_sv,Mn_pvZn,Mn_pvZr_sv,Mo_pvNb_sv,Mo_pvNi_pv,Mo_pvOs_pv,Mo_pvPd_pv,Mo_pvPt,Mo_pvRe_pv,Mo_pvRh_pv,Mo_pvRu_pv,Mo_pvSc_sv,Mo_pvTa_pv,Mo_pvTc_pv,Mo_pvTi_sv,Mo_pvV_sv,Mo_pvW_pv,Mo_pvY_sv,Mo_pvZn,Mo_pvZr_sv,Nb_svNi_pv,Nb_svOs_pv,Nb_svPd_pv,Nb_svPt,Nb_svRe_pv,Nb_svRh_pv,Nb_svRu_pv,Nb_svSc_sv,Nb_svTa_pv,Nb_svTc_pv,Nb_svTi_sv,Nb_svV_sv,Nb_svW_pv,Nb_svY_sv,Nb_svZn,Nb_svZr_sv,Ni_pvOs_pv,Ni_pvPd_pv,Ni_pvPt,Ni_pvRe_pv,Ni_pvRh_pv,Ni_pvRu_pv,Ni_pvSc_sv,Ni_pvTa_pv,Ni_pvTc_pv,Ni_pvTi_sv,Ni_pvV_sv,Ni_pvW_pv,Ni_pvY_sv,Ni_pvZn,Ni_pvZr_sv,Os_pvPd_pv,Os_pvPt,Os_pvRe_pv,Os_pvRh_pv,Os_pvRu_pv,Os_pvSc_sv,Os_pvTa_pv,Os_pvTc_pv,Os_pvTi_sv,Os_pvV_sv,Os_pvW_pv,Os_pvY_sv,Os_pvZn,Os_pvZr_sv,Pd_pvPt,Pd_pvRe_pv,Pd_pvRh_pv,Pd_pvRu_pv,Pd_pvSc_sv,Pd_pvTa_pv,Pd_pvTc_pv,Pd_pvTi_sv,Pd_pvV_sv,Pd_pvW_pv,Pd_pvY_sv,Pd_pvZn,Pd_pvZr_sv,PtRe_pv,PtRh_pv,PtRu_pv,PtSc_sv,PtTa_pv,PtTc_pv,PtTi_sv,PtV_sv,PtW_pv,PtY_sv,PtZn,PtZr_sv,Re_pvRh_pv,Re_pvRu_pv,Re_pvSc_sv,Re_pvTa_pv,Re_pvTc_pv,Re_pvTi_sv,Re_pvV_sv,Re_pvW_pv,Re_pvY_sv,Re_pvZn,Re_pvZr_sv,Rh_pvRu_pv,Rh_pvSc_sv,Rh_pvTa_pv,Rh_pvTc_pv,Rh_pvTi_sv,Rh_pvV_sv,Rh_pvW_pv,Rh_pvY_sv,Rh_pvZn,Rh_pvZr_sv,Ru_pvSc_sv,Ru_pvTa_pv,Ru_pvTc_pv,Ru_pvTi_sv,Ru_pvV_sv,Ru_pvW_pv,Ru_pvY_sv,Ru_pvZn,Ru_pvZr_sv,Sc_svTa_pv,Sc_svTc_pv,Sc_svTi_sv,Sc_svV_sv,Sc_svW_pv,Sc_svY_sv,Sc_svZn,Sc_svZr_sv,Ta_pvTc_pv,Ta_pvTi_sv,Ta_pvV_sv,Ta_pvW_pv,Ta_pvY_sv,Ta_pvZn,Ta_pvZr_sv,Tc_pvTi_sv,Tc_pvV_sv,Tc_pvW_pv,Tc_pvY_sv,Tc_pvZn,Tc_pvZr_sv,Ti_svV_sv,Ti_svW_pv,Ti_svY_sv,Ti_svZn,Ti_svZr_sv,V_svW_pv,V_svY_sv,V_svZn,V_svZr_sv,W_pvY_sv,W_pvZn,W_pvZr_sv,Y_svZn,Y_svZr_sv,ZnZr_sv";
  aurostd::string2tokens(command,list,",");
  for(uint i=0;i<list.size();i++) alloys.push_back(list.at(i));
  alloysmesg.clear();alloysmesg.str(std::string());alloysmesg << "Loading LIBRARY U: Nalloys=" << alloys.size() << endl;
  // DONE
  if(XHOST.APENNSY_USE_SERVER) {  
    for(uint i=0;i<alloys.size();i++) alloysRAW.push_back(init::AFLOW_Projects_Directories("LIB2")+"/RAW/"+alloys.at(i));
  }
  if(XHOST.APENNSY_SERVER_AFLOWLIB_ORG) {
    for(uint i=0;i<alloys.size();i++) alloysRAW.push_back(AFLOWLIB_PROJECT_GNDSTATE+"/"+alloys.at(i));
  }
  // for(uint i=0;i<alloys.size();i++) cerr << i << " " << alloys.at(i) << " " << alloysRAW.at(i) << endl;
  LIB_MODE=LIB_MODEX;
  return TRUE;
}

// **************************************************************************
// LIBRARY PGM
bool APENNSY_Parameters::LibLoadAlloysLIB2PGM(_aflags &aflags) {
  // ------------------------------------------------------------------------
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag"))
    cerr << "APENNSY_Parameters::LibLoadAlloysLIB2PGM(_aflags &aflags)" << endl;
  // ------------------------------------------------------------------------
  // LOAD
  alloys.clear();alloysRAW.clear();
  string command;
  vector<string> list;
  command="AgIr,AgOs_pv,AgPd_pv,AgPt,AgRh_pv,AgRu_pv,AuIr,AuOs_pv,AuPd_pv,AuPt,AuRh_pv,AuRu_pv,CdIr,CdOs_pv,CdPd_pv,CdPt,CdRh_pv,CdRu_pv,CoIr,CoOs_pv,CoPd_pv,CoPt,CoRh_pv,CoRu_pv,Cr_pvIr,Cr_pvOs_pv,Cr_pvPd_pv,Cr_pvPt,Cr_pvRh_pv,Cr_pvRu_pv,Cu_pvIr,Cu_pvOs_pv,Cu_pvPd_pv,Cu_pvPt,Cu_pvRh_pv,Cu_pvRu_pv,Fe_pvIr,Fe_pvOs_pv,Fe_pvPd_pv,Fe_pvPt,Fe_pvRh_pv,Fe_pvRu_pv,Hf_pvIr,Hf_pvOs_pv,Hf_pvPd_pv,Hf_pvPt,Hf_pvRh_pv,Hf_pvRu_pv,HgIr,HgOs_pv,HgPd_pv,HgPt,HgRh_pv,HgRu_pv,IrMn_pv,IrMo_pv,IrNb_sv,IrNi_pv,IrOs_pv,IrPd_pv,IrPt,IrRe_pv,IrRh_pv,IrRu_pv,IrSc_sv,IrTa_pv,IrTc_pv,IrTi_sv,IrV_sv,IrW_pv,IrY_sv,IrZn,IrZr_sv,Mn_pvOs_pv,Mn_pvPd_pv,Mn_pvPt,Mn_pvRh_pv,Mn_pvRu_pv,Mo_pvOs_pv,Mo_pvPd_pv,Mo_pvPt,Mo_pvRh_pv,Mo_pvRu_pv,Nb_svOs_pv,Nb_svPd_pv,Nb_svPt,Nb_svRh_pv,Nb_svRu_pv,Ni_pvOs_pv,Ni_pvPd_pv,Ni_pvPt,Ni_pvRh_pv,Ni_pvRu_pv,Os_pvPd_pv,Os_pvPt,Os_pvRe_pv,Os_pvRh_pv,Os_pvRu_pv,Os_pvSc_sv,Os_pvTa_pv,Os_pvTc_pv,Os_pvTi_sv,Os_pvV_sv,Os_pvW_pv,Os_pvY_sv,Os_pvZn,Os_pvZr_sv,Pd_pvPt,Pd_pvRe_pv,Pd_pvRh_pv,Pd_pvRu_pv,Pd_pvSc_sv,Pd_pvTa_pv,Pd_pvTc_pv,Pd_pvTi_sv,Pd_pvV_sv,Pd_pvW_pv,Pd_pvY_sv,Pd_pvZn,Pd_pvZr_sv,PtRe_pv,PtRh_pv,PtRu_pv,PtSc_sv,PtTa_pv,PtTc_pv,PtTi_sv,PtV_sv,PtW_pv,PtY_sv,PtZn,PtZr_sv,Re_pvRh_pv,Re_pvRu_pv,Rh_pvRu_pv,Rh_pvSc_sv,Rh_pvTa_pv,Rh_pvTc_pv,Rh_pvTi_sv,Rh_pvV_sv,Rh_pvW_pv,Rh_pvY_sv,Rh_pvZn,Rh_pvZr_sv,Ru_pvSc_sv,Ru_pvTa_pv,Ru_pvTc_pv,Ru_pvTi_sv,Ru_pvV_sv,Ru_pvW_pv,Ru_pvY_sv,Ru_pvZn,Ru_pvZr_sv";
  aurostd::string2tokens(command,list,",");
  for(uint i=0;i<list.size();i++) alloys.push_back(list.at(i));
  alloysmesg.clear();alloysmesg.str(std::string());alloysmesg << "Loading LIBRARY PGM: Nalloys=" << alloys.size() << endl;
  // DONE
  if(XHOST.APENNSY_USE_SERVER) {  
    for(uint i=0;i<alloys.size();i++) alloysRAW.push_back(init::AFLOW_Projects_Directories("LIB2")+"/RAW/"+alloys.at(i));
  }
  if(XHOST.APENNSY_SERVER_AFLOWLIB_ORG) {
    for(uint i=0;i<alloys.size();i++) alloysRAW.push_back(AFLOWLIB_PROJECT_GNDSTATE+"/"+alloys.at(i));
  }
  // for(uint i=0;i<alloys.size();i++) cerr << i << " " << alloys.at(i) << " " << alloysRAW.at(i) << endl;
  LIB_MODE=LIB_MODEX;
  return TRUE;
}

// **************************************************************************
// LIBRARY generic
bool APENNSY_Parameters::LibLoadAlloysALLOY(string alloy_name,_aflags &aflags) {
  // ------------------------------------------------------------------------
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag"))
    cerr << "APENNSY_Parameters::LibLoadAlloysALLOY(_aflags &aflags)" << endl;
  // ------------------------------------------------------------------------
  string names=alloy_name,specieA,specieB,specieC;
  if(PseudopotentialNoclean==FALSE) names=KBIN::VASP_PseudoPotential_CleanName(KBIN::VASP_PseudoPotential_CleanName(names));
  // LOAD
  alloys.clear();alloysRAW.clear();
  // Library X
  APENNSY_Parameters paramsX_tmp;
  paramsX_tmp.PseudopotentialNoclean=PseudopotentialNoclean; // inherit cleaniness

  paramsX_tmp.LibLoadAlloysLIB2(aflags);
  if(PseudopotentialNoclean==FALSE) {
    for(uint i=0;i<paramsX_tmp.alloys.size();i++) {
      //  cerr << "[" << names << "] [" << KBIN::VASP_PseudoPotential_CleanName(KBIN::VASP_PseudoPotential_CleanName(paramsX_tmp.alloys.at(i))) << "]" << endl;
      // if(aurostd::substring2bool(KBIN::VASP_PseudoPotential_CleanName(KBIN::VASP_PseudoPotential_CleanName(paramsX_tmp.alloys.at(i))),names))
      if(names==KBIN::VASP_PseudoPotential_CleanName(KBIN::VASP_PseudoPotential_CleanName(paramsX_tmp.alloys.at(i)))) // only one
      { //CO20200106 - patching for auto-indenting
        // if(names==paramsX_tmp.alloys.at(i)) { // only one  //[CO20200106 - close bracket for indenting]}
        cerr << paramsX_tmp.alloys.at(i) << endl;
        alloys.push_back(paramsX_tmp.alloys.at(i));
        alloysRAW.push_back(paramsX_tmp.alloysRAW.at(i));
      }
    }
  } else {
    for(uint i=0;i<paramsX_tmp.alloys.size();i++) {
      //     cerr << names << " * " << paramsX_tmp.alloys.at(i) << endl;
      if(names==paramsX_tmp.alloys.at(i)) { // only one
        alloys.push_back(paramsX_tmp.alloys.at(i));
        alloysRAW.push_back(paramsX_tmp.alloysRAW.at(i));
      }
    }
  }
  // fix verbose
  alloysmesg.clear();alloysmesg.str(std::string());alloysmesg << "Loading LIBRARY " <<  names << ": Nalloys=" << alloys.size() << endl;
  // DONE
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag"))
    for(uint i=0;i<alloys.size();i++)
      cerr << "Loading " << " " << alloys.at(i) << " " << alloysRAW.at(i) << " " << endl;
  LIB_MODE=LIB_MODEX;
  return TRUE;
}

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************



