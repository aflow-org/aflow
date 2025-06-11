// ***************************************************************************
// *                                                                         *
// *              Aflow COREY OSES - Duke University 2003-2021               *
// *                                                                         *
// ***************************************************************************
// Written by Corey Oses 2020

#ifndef _AUROSTD_XPARSER_CPP_
#define _AUROSTD_XPARSER_CPP_

#ifndef _AUROSTD_XPARSER_H_
#include "aurostd_xparser.h"
#endif

namespace aurostd {

  void VASP_PseudoPotential_CleanName_InPlace(string& species,bool capital_letters_only,bool remove_floats) { //CO20190712  //CO20210623 - added remove_floats
    //WARNING: to anyone adding to this list, BE CAREFUL to avoid adding entries that contain capital letters
    //they must be added to CAPITAL_LETTERS_PP_LIST in aurostd.h
    //these pp suffixes cause problems when parsing compounds (capital letters)

    vector<string> vCAPITAL_LETTERS_PP;
    aurostd::string2tokens(CAPITAL_LETTERS_PP_LIST,vCAPITAL_LETTERS_PP,",");
    for(uint i=0;i<vCAPITAL_LETTERS_PP.size();i++){//capital letter ones to watch out for when parsing compounds
      aurostd::RemoveSubStringInPlace(species,vCAPITAL_LETTERS_PP[i]);
    }

    if(capital_letters_only==false){
      //from AFLOW.org database //CO20210315 - must come before .5 (removed below)
      aurostd::RemoveSubStringInPlace(species,"pot_LDA/");
      aurostd::RemoveSubStringInPlace(species,"pot_GGA/");
      aurostd::RemoveSubStringInPlace(species,"pot_PBE/");
      aurostd::RemoveSubStringInPlace(species,"potpaw_LDA/");
      aurostd::RemoveSubStringInPlace(species,"potpaw_GGA/");
      aurostd::RemoveSubStringInPlace(species,"potpaw_PBE/");
      aurostd::RemoveSubStringInPlace(species,"potpaw_LDA.54/");
      aurostd::RemoveSubStringInPlace(species,"potpaw_PBE.54/");

      //general database  //CO20210315 - must come before .5 (removed below)
      aurostd::RemoveSubStringInPlace(species,DEFAULT_VASP_POTCAR_DIR_POT_LDA+"/");
      aurostd::RemoveSubStringInPlace(species,DEFAULT_VASP_POTCAR_DIR_POT_GGA+"/");
      aurostd::RemoveSubStringInPlace(species,DEFAULT_VASP_POTCAR_DIR_POT_PBE+"/");
      aurostd::RemoveSubStringInPlace(species,DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA+"/");
      aurostd::RemoveSubStringInPlace(species,DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA+"/");
      aurostd::RemoveSubStringInPlace(species,DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE+"/");
      aurostd::RemoveSubStringInPlace(species,DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN+"/");
      aurostd::RemoveSubStringInPlace(species,DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN+"/");

      aurostd::RemoveSubStringInPlace(species,"_old");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Si_h_old
      aurostd::RemoveSubStringInPlace(species,".old");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Mg_pv.old
      aurostd::RemoveSubStringInPlace(species,"_vnew");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Pd_vnew
      aurostd::RemoveSubStringInPlace(species,"_new2");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Ti_sv_new2
      aurostd::RemoveSubStringInPlace(species,"_new");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Au_new

      aurostd::RemoveSubStringInPlace(species,"_pvf");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Cu_pvf
      aurostd::RemoveSubStringInPlace(species,"_rel");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Pb_d_rel
      aurostd::RemoveSubStringInPlace(species,"_ref");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Ge_d_GW_ref
      aurostd::RemoveSubStringInPlace(species,"_local");  //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/C_local
      aurostd::RemoveSubStringInPlace(species,"_nopc");  //CO20190712 - potpaw_LDA/potpaw_PBE.20100505/Si_nopc
      aurostd::RemoveSubStringInPlace(species,".nrel");  //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/Ga_pv_GW.nrel
      aurostd::RemoveSubStringInPlace(species,"_nr");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/C_h_nr
      aurostd::RemoveSubStringInPlace(species,"_nc");  //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/H_nc_GW
      aurostd::RemoveSubStringInPlace(species,"_n");  //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/As_GW_n
      aurostd::RemoveSubStringInPlace(species,"_parsv");  //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/Mg_pv_parsv_GW
      aurostd::RemoveSubStringInPlace(species,"_sv2");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Li_sv2
      aurostd::RemoveSubStringInPlace(species,"_sv");
      aurostd::RemoveSubStringInPlace(species,"_vs"); //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/N_vs
      aurostd::RemoveSubStringInPlace(species,"_pv");
      aurostd::RemoveSubStringInPlace(species,"_dr");  //CO20190712 - BEFORE _d //potpaw_LDA/potpaw_LDA.20100505/Pb_dr
      aurostd::RemoveSubStringInPlace(species,"_d3");  //CO20190712 - BEFORE _d //potpaw_PBE/potpaw_PBE.20100506/Ge_d3
      aurostd::RemoveSubStringInPlace(species,"_d2");  //CO20190712 - BEFORE _d //potpaw_LDA/potpaw_LDA.05May2010/As_d2_GW
      aurostd::RemoveSubStringInPlace(species,"_d");
      aurostd::RemoveSubStringInPlace(species,"_soft");  //CO20190712 - BEFORE _s
      aurostd::RemoveSubStringInPlace(species,"_s");
      //[CO20190712 - OBSOLETE really _n and _2]aurostd::RemoveSubStringInPlace(species,"_2_n");
      aurostd::RemoveSubStringInPlace(species,"_h");
      aurostd::RemoveSubStringInPlace(species,"_f");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Cu_f
      aurostd::RemoveSubStringInPlace(species,"_af"); //CO20191110 - SHACHAR aflow pp 

      aurostd::RemoveSubStringInPlace(species,"_1");
      aurostd::RemoveSubStringInPlace(species,"_2");
      aurostd::RemoveSubStringInPlace(species,"_3");

      //CO20210623 - selectively remove floats, this might interfere with extracting composition
      if(remove_floats){
        aurostd::RemoveSubStringInPlace(species,"1.75"); //CO20190712 - potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H1.75
        aurostd::RemoveSubStringInPlace(species,"1.66"); //CO20190712 - potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H1.66
        aurostd::RemoveSubStringInPlace(species,"1.33"); //CO20190712 - potpaw_PBE.54/potpaw_PBE.54.04Sep2015/H1.33
        aurostd::RemoveSubStringInPlace(species,"1.25"); //CO20190712 - before all other decimal numbers
        aurostd::RemoveSubStringInPlace(species,"1.5"); //CO20190712 - potpaw_PBE/potpaw_PBE.06May2010/H1.5
        aurostd::RemoveSubStringInPlace(species,".75");  //CO20190712 - before 0.5
        aurostd::RemoveSubStringInPlace(species,".25");  //CO20190712 - potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H.25
        aurostd::RemoveSubStringInPlace(species,".66"); //CO20190712 - potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H.66
        aurostd::RemoveSubStringInPlace(species,".33"); //CO20190712 - potpaw_PBE.54/potpaw_PBE.54.04Sep2015/H.33
        aurostd::RemoveSubStringInPlace(species,".42"); //CO20190712 - potpaw_PBE.54/potpaw_PBE.54.04Sep2015/H.42
        aurostd::RemoveSubStringInPlace(species,".58"); //CO20190712 - before 0.5 //potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H.58
        aurostd::RemoveSubStringInPlace(species,".5");
      }

      aurostd::RemoveSubStringInPlace(species,"+1");
      aurostd::RemoveSubStringInPlace(species,"+3");
      aurostd::RemoveSubStringInPlace(species,"+5");
      aurostd::RemoveSubStringInPlace(species,"+7");
      aurostd::RemoveSubStringInPlace(species,"-1");
      aurostd::RemoveSubStringInPlace(species,"-3");
      aurostd::RemoveSubStringInPlace(species,"-5");
      aurostd::RemoveSubStringInPlace(species,"-7");

      aurostd::RemoveSubStringInPlace(species,"__"); //CO20190712 - BEFORE _ - potpaw_LDA/potpaw_LDA.05May2010/Si_sv_GW__
      aurostd::RemoveSubStringInPlace(species,"_");  //CO20190712  //potpaw_LDA/potpaw_LDA.05May2010/Si_sv_GW_
    }
  }

  //use only as a supplement for getElements(), do NOT use outside
  //this assumes a very simple Mn2Pd5, non-stoich is ok (e.g., Mn2.5Pd5)
  //no pseudo potential specification
  //no junk at the end (_ICSD_, :LDAU2, :PAW_PBE, .OLD, etc.), pre-process before
  //this is FASTER than getElements(), but not as robust for general input (specialized)
  void elementsFromCompositionString(const string& input,vector<string>& velements){vector<double> vcomposition;return elementsFromCompositionString(input,velements,vcomposition);}  //CO20190712
  template<class utype> void elementsFromCompositionString(const string& input,vector<string>& velements,vector<utype>& vcomposition){ //CO20190712
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    velements.clear();
    vcomposition.clear();  //ME20190628

    //////////////////////////////////////////////////////////////////////////////
    // START Checks for correct input by counting number of uppercase letters
    //////////////////////////////////////////////////////////////////////////////

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " original input=" << input << endl;}

    //CO20180409 - running through input twice, no need, simply check at the end
    //uint numberOfElements = 0;
    //for (uint i = 0; i < input.size(); i++) {
    //  if(isupper(input[i])) {
    //    numberOfElements++;
    //  }
    //}
    //if(numberOfElements == 0) {
    //  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Elements must be properly capitalized", FileMESSAGE, oss, _LOGGER_ERROR_);
    //  return velements;
    //}

    //////////////////////////////////////////////////////////////////////////////
    // END Checks for correct input by counting number of uppercase letters
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // START Parsing input
    //////////////////////////////////////////////////////////////////////////////

    //CO20180316 - fixed this function to be simpler, too complicated before
    string auxstr;
    for (uint i = 0; i < input.size(); i++) {
      if(isupper(input[i])) {
        auxstr.clear();
        auxstr+=input[i++];
        while (((i < input.size()) && (input[i]>='a' && input[i]<='z') )){auxstr+=input[i++];}
        i--;  //very important since we increase at the top of the loop (equivalent to i+j-1)
        //while (((i < input.size()) && isalpha(input[i]) && !isupper(input[i]))){auxstr+=input[i++];}
        //while(!(clean && !isalpha(input[i]))) //(input[i]=='_' || input[i]==':' || input[i]=='.' || isdigit(input[i]))))
        //isalpha() saves us from all issues with VASP_PseudoPotential_CleanName() except, e.g., potpaw_PBE/Na, we took care of that above
        //if(clean)
        //{ //CO20200106 - patching for auto-indenting
        //  auxstr = KBIN::VASP_PseudoPotential_CleanName(auxstr);  //fix vasp pp
        //  //CO20180409 - again, no need to run through essentially a third time, we already cut at these characters
        //  //look for bad characters and cut the string
        //  //for(uint j=1;j<auxstr.size();j++){
        //  //  if(auxstr[j]=='_' || auxstr[j]==':' || isdigit(auxstr[j])){auxstr=auxstr.substr(0,j);break;}  //fix aflow stuff like ':'
        //  //}
        //}
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " element found: " << auxstr << endl;}
        velements.push_back(auxstr);
        //ME20190628 - get composition, too
      } else if ( (input[i]>='0' && input[i]<='9') || (input[i] == '.')) {  //CO20190712 - just in case we have H.25 (not good form but try to catch anyway, never produced by aflow automatically)
        auxstr.clear();
        auxstr += input[i++];
        while ((i < input.size()) && ( (input[i]>='0' && input[i]<='9') || (input[i] == '.'))) {auxstr += input[i++];}
        i--;
        if (LDEBUG) {
          std::cerr << __AFLOW_FUNC__ << " found element count: " << auxstr << " of element " << (velements.size() - 1) << ".";
          if (vcomposition.size() != velements.size()) {
            std::cerr << " Will add ones to elements " << vcomposition.size() << " to " << (velements.size() - 2) << ".";
          }
          std::cerr << std::endl;
        }
        // Add implicit ones
        for (uint i = vcomposition.size(); i < velements.size() - 1; i++){vcomposition.push_back((utype)1.0);}
        vcomposition.push_back(aurostd::string2utype<utype>(auxstr));
      }
    }
    // Add implicit ones
    for (uint i = vcomposition.size(); i < velements.size(); i++) vcomposition.push_back((utype)1.0);
  }

  //use only as a supplement for getElements(), do NOT use outside
  //this assumes Mn_pvPt
  //no composition information
  //no junk at the end (_ICSD_, :LDAU2, :PAW_PBE, .OLD, etc.), pre-process before
  void elementsFromPPString(const string& input,vector<string>& velements,bool keep_pp){ //CO20190712
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    velements=getElements(input);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " velements=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(velements,"\""),",") << endl;}
    if(keep_pp==false){return;}

    //copy info into vspecies and clear velements
    vector<string> vspecies;
    for(uint i=0;i<velements.size();i++){vspecies.push_back(velements[i]);}
    velements.clear();

    //parse string around these elements
    string::size_type loc1=0,loc2=string::npos;
    vector<string> vCAPITAL_LETTERS_PP;
    aurostd::string2tokens(CAPITAL_LETTERS_PP_LIST,vCAPITAL_LETTERS_PP,",");
    bool found_CAPITAL_LETTERS_PP=false;
    bool found_CAPITAL_LETTERS=false;
    for(uint i=0;i<vspecies.size();i++){
      if((i+1)>=vspecies.size()){loc2=string::npos;}
      else{loc2=input.find(vspecies[i+1],loc1);}
      while(loc2!=string::npos){
        found_CAPITAL_LETTERS_PP=false;
        for(uint j=0;j<vCAPITAL_LETTERS_PP.size()&&found_CAPITAL_LETTERS_PP==false;j++){
          if((loc2-(vCAPITAL_LETTERS_PP[j].size()-1))<input.size()){continue;}
          found_CAPITAL_LETTERS=true;
          for(uint k=0;k<vCAPITAL_LETTERS_PP[j].size()&&found_CAPITAL_LETTERS==true;k++){
            if(input[loc2-k]!=vCAPITAL_LETTERS_PP[j][vCAPITAL_LETTERS_PP[j].size()-k-1]){found_CAPITAL_LETTERS=false;}
          }
          if(found_CAPITAL_LETTERS){found_CAPITAL_LETTERS_PP=true;}
        }
        if(found_CAPITAL_LETTERS_PP==false){break;}
        //[OBSOLETE] Do not pick W from _GW (tungsten)
        //[OBSOLETE]if (!( (loc2-2)<input.size() && (input[loc2] == 'W') && (input[loc2-1] == 'G') && (input[loc2-2] == '_') )){break;} //(loc2-2)<input.size() because loc2 is utype, it's always >0, loc2-2 can wrap around to a big number though
        loc2=input.find(vspecies[i],loc2+1);
      }
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " loc1=" << loc1 << ", loc2=" << loc2 << endl;}
      if(loc2==string::npos){ //CO20210315
        velements.push_back(input.substr(loc1));
        break;
      }else{
        velements.push_back(input.substr(loc1,loc2-loc1));  //loc2-loc1 because it is the distance
        loc1=loc2;
      }
    }

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " velements=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(velements,"\""),",") << endl;}

  }

  // ***************************************************************************
  // aurostd::getElements(string input,ostream&
  // oss,ofstream& FileMESSAGE)
  // ***************************************************************************
  // returns UNSORTED vector<string> from string
  vector<string> getElements(const string& input){ //CO20190712 //borrowed from XATOM_SplitAlloySpecies() //slow since we create many strings, but definitely works
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " original input=\"" << input << "\"" << endl;}
    string alloy=input;
    //[CO20190712 - no need for multiple passes anymore]for(uint i=1;i<=2;i++){alloy=KBIN::VASP_PseudoPotential_CleanName(alloy);} //be certain you clean everything, especially _GW (worst offender)
    aurostd::VASP_PseudoPotential_CleanName_InPlace(alloy); //be certain you clean everything, especially _GW (worst offender)
    aurostd::RemoveNumbersInPlace(alloy);              // remove composition
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " cleaned input=\"" << alloy << "\"" << endl;}
    vector<string> vspecies;
    for(uint i=0;i<alloy.length();i++) {
      if(alloy[i]>='A' && alloy[i]<='Z') vspecies.push_back("");
      vspecies.back()+=alloy[i];
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " vspecies pre ASCII clean=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vspecies,"\""),",") << endl;}
    for(uint i=0;i<vspecies.size();i++){aurostd::CleanStringASCII_InPlace(vspecies[i]);}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " vspecies post ASCII clean=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vspecies,"\""),",") << endl;}
    return vspecies;
  }
  vector<string> getElements(const string& input,elements_string_type e_str_type,bool clean,bool sort_elements,bool keep_pp,ostream& oss) {  // overload
    ofstream FileMESSAGE;
    return getElements(input,e_str_type,FileMESSAGE,clean,sort_elements,keep_pp,oss);
  }
  //ME20190628 - added variant that also determines the composition
  template<class utype> vector<string> getElements(const string& input,vector<utype>& vcomposition,bool clean,bool sort_elements,bool keep_pp,ostream& oss) {
    ofstream FileMESSAGE;
    return getElements(input,vcomposition,composition_string,FileMESSAGE,clean,sort_elements,keep_pp,oss);  //this gets composition_string by default, pp_string has no composition
  }
  //cannot deduce utype from this construction
  vector<string> getElements(const string& input,elements_string_type e_str_type,ofstream& FileMESSAGE,bool clean,bool sort_elements,bool keep_pp,ostream& oss) {  // overload
    vector<double> vcomposition;
    return getElements(input,vcomposition,e_str_type,FileMESSAGE,clean,sort_elements,keep_pp,oss);
  }
  template<class utype> vector<string> getElements(const string& input,vector<utype>& vcomposition,elements_string_type e_str_type,bool clean,bool sort_elements,bool keep_pp,ostream& oss) { // overload
    ofstream FileMESSAGE;
    return getElements(input,vcomposition,e_str_type,FileMESSAGE,clean,sort_elements,keep_pp,oss);  //this gets composition_string by default, pp_string has no composition
  }
  template<class utype> vector<string> getElements(const string& _input,vector<utype>& vcomposition,elements_string_type e_str_type,ofstream& FileMESSAGE,bool clean,bool sort_elements,bool keep_pp,ostream& oss) { // main function
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    vector<string> velements;
    vcomposition.clear();  //ME20190628

    //////////////////////////////////////////////////////////////////////////////
    // START Checks for correct input by counting number of uppercase letters
    //////////////////////////////////////////////////////////////////////////////

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " original input=" << _input << endl;}

    if(_input.empty()) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Empty input", FileMESSAGE, oss, _LOGGER_ERROR_);
      return velements;
    }

    string input=_input;

    if(clean && (e_str_type==composition_string || (e_str_type==pp_string && keep_pp==false))){
      //in case we run into potpaw_PBE/Na, but only works for single elements, must be before check for isupper(input[0]) 
      bool capital_letters_only=false;  //default
      bool remove_floats=false; //CO20210623 - explicitly KEEP floats for composition
      aurostd::VASP_PseudoPotential_CleanName_InPlace(input,capital_letters_only,remove_floats);
    }

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " checking input [1] =" << input << endl;}

    if(!isupper(input[0])) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Elements must be properly capitalized (input="+input+")", FileMESSAGE, oss, _LOGGER_ERROR_);
      return velements;
    }

    //we have a LIB1 problem... grab first everything before :
    //this is safe, as aflow generally introduces : in prototype, e.g., :LDAU2
    //this is safe anyway because elements would be BEFORE :
    if(clean){
      //FAST
      string::size_type loc;
      //:
      loc=input.find(':');input=input.substr(0,loc);
      //_ICSD_
      loc=input.find("_ICSD_");input=input.substr(0,loc);
      //SLOW
      //vector<string> tokens;
      //aurostd::string2tokens(input,tokens,":");
      //input=tokens[0];
    }

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " checking input [2] =" << input << endl;}

    //CO20180409 - running through input twice, no need, simply check at the end
    //uint numberOfElements = 0;
    //for (uint i = 0; i < input.size(); i++) {
    //  if(isupper(input[i])) {
    //    numberOfElements++;
    //  }
    //}
    //if(numberOfElements == 0) {
    //  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Elements must be properly capitalized", FileMESSAGE, oss, _LOGGER_ERROR_);
    //  return velements;
    //}

    //////////////////////////////////////////////////////////////////////////////
    // END Checks for correct input by counting number of uppercase letters
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // START Parsing input
    //////////////////////////////////////////////////////////////////////////////

    if(e_str_type==composition_string){elementsFromCompositionString(input,velements,vcomposition);}
    else if(e_str_type==pp_string){elementsFromPPString(input,velements,keep_pp);}
    else{
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown compound designation",_INPUT_ILLEGAL_);
    }

    if(clean){
      for(uint i=0;i<velements.size();i++){aurostd::CleanStringASCII_InPlace(velements[i]);}  //CO20190712 - extra cleaning from XATOM_SplitAlloySpecies
    }

    // Add implicit ones
    for (uint i = vcomposition.size(); i < velements.size(); i++) vcomposition.push_back((utype)1.0);

    //////////////////////////////////////////////////////////////////////////////
    // END Parsing input
    //////////////////////////////////////////////////////////////////////////////

    if(velements.size()==0){pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "No elements found", FileMESSAGE, oss, _LOGGER_ERROR_);}

    if(sort_elements && velements.size()>1){
      //this is MORE efficient that std::swap which has a copy constructor inside
      //http://www.cplusplus.com/reference/algorithm/swap/
      string etmp="";
      utype ctmp=(utype)0.0;
      for(uint i=0;i<velements.size()-1;i++){
        for(uint j=i+1;j<velements.size();j++){
          if(velements[i]>velements[j]){
            etmp=velements[j];  //fix old j
            velements[j]=velements[i];  //swap
            velements[i]=etmp;  //set i to old j
            ctmp=vcomposition[j]; //fix old j
            vcomposition[j]=vcomposition[i];  //swap
            vcomposition[i]=ctmp; //set i to old j
          }
        }
      }
    }

    return velements;
  }
} // namespace aurostd

//AS20201214 BEGIN JSONwriter
namespace aurostd {
  //***************************************************************************
  void JSONwriter::free() { content.clear(); }
  void JSONwriter::clear() { free(); }
  void JSONwriter::copy(const JSONwriter &jw){
    content = jw.content;
  }

  JSONwriter::JSONwriter() { free(); }
  JSONwriter::JSONwriter(const JSONwriter &jw)
  {
    free();
    copy(jw);
  }
  JSONwriter::~JSONwriter() { free(); }

  const JSONwriter& JSONwriter::operator=(const JSONwriter &jw){
    copy(jw);
    return *this;
  }
  //***************************************************************************
  //  addNumber block
  template <typename utype> void JSONwriter::addNumber(const string &key,
      const utype value)
  {
    content.push_back("\"" + key + "\":" + aurostd::utype2string(value));
  }

  //***************************************************************************
  //  addVector block
  template <typename utype> void JSONwriter::addVector(const string &key,
      const utype &value)
  {
    content.push_back("\"" + key + "\":[" + aurostd::joinWDelimiter(value, ",") + "]");
  }

  void JSONwriter::addVector(const string &key, const vector<double> &value,
      int precision, bool roundoff, double tol)
  {
    content.push_back("\"" + key + "\":[");
    content.back() += aurostd::joinWDelimiter(
        aurostd::vecDouble2vecString(value, precision, roundoff, tol),",");
    content.back() += "]";
  }

  void JSONwriter::addVector(const string &key, const deque<double> &value, int precision,
      bool roundoff, double tol)
  {
    content.push_back("\"" + key + "\":[");
    content.back() += aurostd::joinWDelimiter(
        aurostd::vecDouble2vecString(value, precision, roundoff, tol),",");
    content.back() += "]";
  }

  void JSONwriter::addVector(const string &key, const xvector<double> &value, int precision,
      bool roundoff, double tol) //DX20210308 - added xvector variant
  {
    content.push_back("\"" + key + "\":[");
    content.back() += aurostd::joinWDelimiter(
        aurostd::xvecDouble2vecString(value, precision, roundoff, tol),",");
    content.back() += "]";
  }

  void JSONwriter::addVector(const string &key, const vector<string> &value, bool wrap) //DX20210301 - added wrap
  {
    content.push_back("\"" + key + "\":[");
    if(wrap){content.back() += aurostd::joinWDelimiter(aurostd::wrapVecEntries(value, "\""), ","); }
    else{ content.back() += aurostd::joinWDelimiter(value, ","); }
    content.back() += "]";
  }

  void JSONwriter::addVector(const string &key, const deque<string> &value, bool wrap) //DX20210301 - added wrap
  {
    content.push_back("\"" + key + "\":[");
    if(wrap){content.back() += aurostd::joinWDelimiter(aurostd::wrapVecEntries(value, "\""), ","); }
    else{ content.back() += aurostd::joinWDelimiter(value, ","); }
    content.back() += "]";
  }

  void JSONwriter::addVector(const string &key, vector<JSONwriter> &value) //AS20210309
  {
    content.push_back("\"" + key + "\":[");
    uint size = value.size();
    if (size){
      for (uint i=0; i<size-1; i++){
        content.back() += value[i].toString() + ",";
      }
      content.back() += value[size-1].toString();
    }
    content.back() += "]";
  }

  //***************************************************************************
  //  addMatrix block
  template <typename utype> void JSONwriter::addMatrix(const string &key,
      const utype &value)
  {
    content.push_back("\"" + key + "\":[");
    vector<string> matrix;
    for (uint i=0; i<value.size(); i++){
      matrix.push_back("[" + aurostd::joinWDelimiter(value[i], ",") + "]");
    }
    content.back() += aurostd::joinWDelimiter(matrix, ",");
    content.back() += "]";
  }

  void JSONwriter::addMatrix(const string &key, const vector<vector<double> > &value,
      int precision, bool roundoff, double tol)
  {
    content.push_back("\"" + key + "\":[");
    vector<string> matrix;
    for (uint i=0; i<value.size(); i++){
      matrix.push_back("[" + aurostd::joinWDelimiter(
            aurostd::vecDouble2vecString(value[i], precision, roundoff, tol), ",") + "]");
    }
    content.back() += aurostd::joinWDelimiter(matrix, ",");
    content.back() += "]";
  }

  void JSONwriter::addMatrix(const string &key, const deque<deque<double> > &value,
      int precision, bool roundoff, double tol)
  {
    content.push_back("\"" + key + "\":[");
    vector<string> matrix;
    for (uint i=0; i<value.size(); i++){
      matrix.push_back("[" + aurostd::joinWDelimiter(
            aurostd::vecDouble2vecString(value[i], precision, roundoff, tol), ",") + "]");
    }
    content.back() += aurostd::joinWDelimiter(matrix, ",");
    content.back() += "]";
  }

  void JSONwriter::addMatrix(const string &key, const xmatrix<double> &value,
      int precision, bool roundoff, double tol) //DX20210308 - added xmatrix variant
  {
    content.push_back("\"" + key + "\":[");
    content.back() += xmatDouble2String(value, precision, roundoff, tol);
    content.back() += "]";
  }

  void JSONwriter::addMatrix(const string &key, const vector<vector<string> > &value) //DX20210211
  {
    content.push_back("\"" + key + "\":[");
    vector<string> matrix;
    for (uint i=0; i<value.size(); i++){
      matrix.push_back("[" + aurostd::joinWDelimiter(
            aurostd::wrapVecEntries(value[i],"\""), ",") + "]");
    }
    content.back() += aurostd::joinWDelimiter(matrix, ",");
    content.back() += "]";
  }

  void JSONwriter::addMatrix(const string &key, const deque<deque<string> > &value) //DX20210211
  {
    content.push_back("\"" + key + "\":[");
    vector<string> matrix;
    for (uint i=0; i<value.size(); i++){
      matrix.push_back("[" + aurostd::joinWDelimiter(
            aurostd::wrapVecEntries(value[i],"\""), ",") + "]");
    }
    content.back() += aurostd::joinWDelimiter(matrix, ",");
    content.back() += "]";
  }

  //***************************************************************************
  void JSONwriter::addString(const string &key, const string &value)
  {
    content.push_back("\"" + key + "\":" + "\"" + value + "\"");
  }

  //***************************************************************************
  void JSONwriter::addBool(const string &key, bool value)
  {
    content.push_back("\"" + key + "\":" + (value ? "true" : "false"));
  }

  //***************************************************************************
  void JSONwriter::mergeRawJSON(const string &value) //DX20210304 - changed name from addRaw to mergeRawJSON
  {
    content.push_back(value);
  }

  //***************************************************************************
  /// Add null value to a key //DX20210301
  void JSONwriter::addNull(const string &key)
  {
    content.push_back("\"" + key + "\":null");
  }

  //***************************************************************************
  /// Add "raw" value for a particular key //DX20210301
  /// Enables JSONs to be passed in as strings, e.g., xstructure2json()
  void JSONwriter::addRawJSON(const string &key, const string& value)
  {
    content.push_back("\"" + key + "\":" + value);
  }

  //***************************************************************************
  /// Adds JSON object as a new value (i.e. with curly braces)
  void JSONwriter::addJSON(const string &key, JSONwriter &value)
  {
    content.push_back("\"" + key + "\":" + value.toString(true));
  }

  //***************************************************************************
  /// Merges/inserts a given JSON object (i.e. with no curly braces)
  void JSONwriter::mergeJSON(JSONwriter &value)
  {
    content.push_back(value.toString(false));
  }

  //***************************************************************************
  string JSONwriter::toString(bool wrap)
  {
    if (wrap){
      return "{" + aurostd::joinWDelimiter(content, ",") + "}";
    }
    else{
      return aurostd::joinWDelimiter(content, ",");
    }
  }
}
//AS20201214 END

//***************************************************************************
// BENCHMARK FOR JSONWriter //AS20210309
//***************************************************************************
// Comparing vector of strings vs vector of JSONs
// Results from 20210309:
//   vector of strings 474.436 ms
//   vector of jsons   481.034 ms
// Keeping benchmark to test any performance degredation with new objects
// Need header<random> to run:

//AS20210309 [BENCHMARK] #include <random>
//AS20210309 [BENCHMARK] namespace aurostd {
//AS20210309 [BENCHMARK] 	void test_vector_string(int niterations)
//AS20210309 [BENCHMARK] 	{
//AS20210309 [BENCHMARK] 		std::random_device rd;
//AS20210309 [BENCHMARK] 		std::mt19937 mt(rd());
//AS20210309 [BENCHMARK] 		std::uniform_int_distribution<int> dist(0,100);
//AS20210309 [BENCHMARK]
//AS20210309 [BENCHMARK] 		aurostd::JSONwriter json;
//AS20210309 [BENCHMARK]
//AS20210309 [BENCHMARK] 		vector<string> set;
//AS20210309 [BENCHMARK] 		xvector<double> position;
//AS20210309 [BENCHMARK] 		for(int i=0;i<niterations;i++){
//AS20210309 [BENCHMARK] 			json.clear();
//AS20210309 [BENCHMARK] 			position.clear();
//AS20210309 [BENCHMARK]
//AS20210309 [BENCHMARK] 			position(1) = dist(mt);
//AS20210309 [BENCHMARK] 			position(2) = dist(mt);
//AS20210309 [BENCHMARK] 			position(3) = dist(mt);
//AS20210309 [BENCHMARK] 			json.addVector("position", position, _AFLOWLIB_DATA_GEOMETRY_PREC_);
//AS20210309 [BENCHMARK] 			json.addString("name", "Al");
//AS20210309 [BENCHMARK] 			json.addNumber("mulitiplicity", dist(mt));
//AS20210309 [BENCHMARK] 			json.addString("Wyckoff_letter", "a");
//AS20210309 [BENCHMARK] 			json.addString("site_symmetry",  "1");
//AS20210309 [BENCHMARK] 			set.push_back(json.toString(true));
//AS20210309 [BENCHMARK] 		}
//AS20210309 [BENCHMARK]
//AS20210309 [BENCHMARK] 		aurostd::JSONwriter json_out;
//AS20210309 [BENCHMARK] 		json_out.addVector("test", set, false);
//AS20210309 [BENCHMARK]
//AS20210309 [BENCHMARK] 		//  cout << json_out.toString() << endl;
//AS20210309 [BENCHMARK] 	}
//AS20210309 [BENCHMARK]
//AS20210309 [BENCHMARK] 	void test_vector_json(int niterations)
//AS20210309 [BENCHMARK] 	{
//AS20210309 [BENCHMARK] 		std::random_device rd;
//AS20210309 [BENCHMARK] 		std::mt19937 mt(rd());
//AS20210309 [BENCHMARK] 		std::uniform_int_distribution<int> dist(0,100);
//AS20210309 [BENCHMARK]
//AS20210309 [BENCHMARK] 		aurostd::JSONwriter json;
//AS20210309 [BENCHMARK] 		vector<aurostd::JSONwriter> set;
//AS20210309 [BENCHMARK]
//AS20210309 [BENCHMARK] 		xvector<double> position;
//AS20210309 [BENCHMARK] 		for(int i=0;i<niterations;i++){
//AS20210309 [BENCHMARK] 			json.clear();
//AS20210309 [BENCHMARK]
//AS20210309 [BENCHMARK] 			position(1) = dist(mt);
//AS20210309 [BENCHMARK] 			position(2) = dist(mt);
//AS20210309 [BENCHMARK] 			position(3) = dist(mt);
//AS20210309 [BENCHMARK] 			json.addVector("position", position, _AFLOWLIB_DATA_GEOMETRY_PREC_);
//AS20210309 [BENCHMARK] 			json.addString("name", "Al");
//AS20210309 [BENCHMARK] 			json.addNumber("mulitiplicity", dist(mt));
//AS20210309 [BENCHMARK] 			json.addString("Wyckoff_letter", "a");
//AS20210309 [BENCHMARK] 			json.addString("site_symmetry",  "1");
//AS20210309 [BENCHMARK] 			set.push_back(json);
//AS20210309 [BENCHMARK] 		}
//AS20210309 [BENCHMARK]
//AS20210309 [BENCHMARK] 		aurostd::JSONwriter json_out;
//AS20210309 [BENCHMARK] 		json_out.addVector("test", set);
//AS20210309 [BENCHMARK]
//AS20210309 [BENCHMARK] 		//  cout << json_out.toString() << endl;
//AS20210309 [BENCHMARK] 	}
//AS20210309 [BENCHMARK]
//AS20210309 [BENCHMARK] 	void run_vector_string_vs_json_test(void)
//AS20210309 [BENCHMARK] 	{
//AS20210309 [BENCHMARK] 		int num = 10000;
//AS20210309 [BENCHMARK] 		int n_avg = 20;
//AS20210309 [BENCHMARK]
//AS20210309 [BENCHMARK] 		std::clock_t start, duration;
//AS20210309 [BENCHMARK]
//AS20210309 [BENCHMARK] 		for (int i=0; i<n_avg; i++){
//AS20210309 [BENCHMARK] 			start = std::clock();
//AS20210309 [BENCHMARK] 			test_vector_string(num);
//AS20210309 [BENCHMARK] 			duration += std::clock() - start;
//AS20210309 [BENCHMARK] 		}
//AS20210309 [BENCHMARK] 		cout << "vector of strings " << 1000.0 * duration/CLOCKS_PER_SEC/n_avg << " ms\n";
//AS20210309 [BENCHMARK]
//AS20210309 [BENCHMARK] 		duration = 0;
//AS20210309 [BENCHMARK] 		for (int i=0; i<n_avg; i++){
//AS20210309 [BENCHMARK] 			start = std::clock();
//AS20210309 [BENCHMARK] 			test_vector_json(num);
//AS20210309 [BENCHMARK] 			duration += std::clock() - start;
//AS20210309 [BENCHMARK] 		}
//AS20210309 [BENCHMARK] 		cout << "vector of jsons " << 1000.0 * duration/CLOCKS_PER_SEC/n_avg << " ms\n";
//AS20210309 [BENCHMARK] 	}
//AS20210309 [BENCHMARK] }

//ME2020408 - JSON reader
//Moved from the AflowDB class
namespace aurostd {

  //extractJsonKeysAflow//////////////////////////////////////////////////////
  // This function extracts keys from an aflowlib.json file. It is much
  // faster than using SQLite's JSON extension, but was designed to only
  // work for the aflowlib.json. It cannot handle nested JSONs!
  vector<string> extractJsonKeysAflow(const string& json) {
    vector<string> keys;
    string substring = "";
    string::size_type pos = 0, lastPos = 0, dpos = 0, quote1  = 0, quote2 = 0, colon = 0;

    // Find the first comma - this is either the end of the key-value pair
    // or part of an array. Either way, the key is inside.
    pos = json.find(",");
    lastPos = 1;  // First character is a curly brace, so skip
    dpos = pos - lastPos;
    while ((pos != string::npos) || (lastPos != string::npos)) {
      // A comma could be separating a key-value pair an array
      // or numbers or strings
      substring = json.substr(lastPos, dpos);

      // Find the colon - if there is no colon, it cannot be a key-value pair
      colon = substring.find(":");
      if (colon != string::npos) {
        // A key is enclosed in quotes, so there must be at least two of them
        quote1 = substring.find("\"");
        if (quote1 != string::npos) {
          quote2 = substring.find("\"", quote1 + 1);
          // Most non-keys are filtered out by now. There could still be array
          // elements left. In that case, however, the colon is between the quotes,
          // so make sure that the first two quotes appear before the colon and
          // take everything in-between as the key. This breaks if quotes, colons,
          // and commas are inside a string in the right sequence, but should not
          // be the case in AFLOW's JSON files.
          if ((quote2 != string::npos) && (quote1 < colon) && (quote2 < colon)) {
            substring = substring.substr(quote1 + 1, quote2 - quote1 - 1);
            if (!substring.empty()) keys.push_back(substring);
          }
        }
      }
      // Move on to the next comma
      lastPos = json.find_first_not_of(",", pos);
      pos = json.find(",", lastPos);
      dpos = pos - lastPos;
    }
    return keys;
  }

  //extractJsonValueAflow/////////////////////////////////////////////////////
  // This function extracts values from an aflowlib.json file. It is much
  // faster than using SQLite's JSON extension, but has was designed to only
  // work for the aflowlib.json. It cannot handle nested JSONs!
  string extractJsonValueAflow(const string& json, string key) {
    string value = "";
    key = "\"" + key + "\":";
    string::size_type start = 0, end = 0;
    start = json.find(key);
    if (start != string::npos) {
      start += key.length();
      end = json.find("\":", start);
      if (end != string::npos) {
        // In case there is any white space between key and value
        value = aurostd::RemoveWhiteSpacesFromTheFront(json.substr(start, end - start));
        // If we have a nested object, "value" should only be '{' + white space by now.
        if (value[0] == '{') {
          string message = "JSON parser cannot read nested objects.";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
        }
        end = value.find_last_of(",");
        value = value.substr(0, end);
      } else {
        end = json.find("}", start);
        // In case there is any white space between key and value
        value = aurostd::RemoveWhiteSpacesFromTheFront(json.substr(start, end - start));
        // If we have a nested object, it should start with '{'
        if (value[0] == '{') {
          string message = "JSON parser cannot read nested objects.";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
        }
      }
    }
    return value;
  }

  //extractJsonVectorAflow/////////////////////////////////////////////////////
  vector<string> extractJsonVectorAflow(const string& json, string key) { //SD20220504
    string value = extractJsonValueAflow(json, key);
    string::size_type start = value.find("["), stop = value.rfind("]");
    vector<string> vec;
    aurostd::string2tokens(value.substr(start + 1, stop - 1), vec, ",");
    return vec;
  }

  //extractJsonMatrixAflow/////////////////////////////////////////////////////
  vector<vector<string>> extractJsonMatrixAflow(const string& json, string key) { //SD20220504
    string value = extractJsonValueAflow(json, key);
    string::size_type start = value.find("[["), stop = value.rfind("]]");
    vector<vector<string>> mat;
    vector<string> vec, tokens;
    aurostd::string2tokensByDelimiter(value.substr(start + 2, stop - 2), tokens, "],[");
    for (uint i = 0; i < tokens.size(); i++) {
      aurostd::string2tokens(tokens[i], vec, ",");
      mat.push_back(vec);
    }
    return mat;
  }

}

#endif // _AUROSTD_XPARSER_CPP_d

// **************************************************************************
// *                                                                        *
// *              Aflow COREY OSES - Duke University 2003-2021              *
// *                                                                        *
// **************************************************************************
