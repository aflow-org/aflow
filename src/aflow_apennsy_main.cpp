// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

#ifndef _AFLOW_APENNSY_MAIN_CPP_
#define _AFLOW_APENNSY_MAIN_CPP_
#include "aflow.h"
#include "aflow_apennsy.h"
// [OBSOLETE] #include "aflow_contrib_junkai_basic.h"
// [OBSOLETE] #include "aflow_contrib_junkai_protocheck.h"

//#include "aflow_xvaspin.h"   // TO REMOVE IN THE FUTURE
//#define _APENNSY_COPYRIGHT_  string("Curtarolo AFLOW http://"+AFLOWLIB_MATERIALS_SERVER)
//#define _APENNSY_COPYRIGHT_  string("http://+AFLOWLIB_MATERIALS_SERVER +" -  http://msg.byu.edu")
#define _APENNSY_COPYRIGHT_  string("http://"+AFLOWLIB_MATERIALS_SERVER)
#define _APENNSY_MAX_STRUCTURES_NUMBER_ 1000

#define GNUPLOT_BITMAP_SIZE string("size 1200,800")
#define GNUPLOT_VECTOR_SIZE string("size 50cm,30cm")
#define LATEX_TEXT_HEIGHT   string("\\textheight 268mm %LATEX")

// vector<string> vrefs;
// SystemReferences(alloys.at(k),vrefs);
// for(uint i=0;i<vrefs.size();i++) oss << strfile1 << "REFERENCE: " << vrefs.at(i) << endl;

// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
uint ApennsyARGs(vector<string> &argv,vector<string> &cmds,aurostd::xoption &vflag) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  // GENERAL STUFF

  vflag.flag("APENNSY::HELP",aurostd::args2flag(argv,cmds,"--HELP|--help"));
  vflag.flag("APENNSY::VERBOSE_flag",aurostd::args2flag(argv,cmds,"--VERBOSE|--verbose"));
  vflag.flag("APENNSY_LIST_flag",aurostd::args2flag(argv,cmds,"--LIST|--list"));

  //  XHOST.APENNSY_USE_SERVER=FALSE;XHOST.APENNSY_USE_LIBRARY=FALSE;XHOST.APENNSY_SERVER_AFLOWLIB_ORG=FALSE;
  vflag.flag("APENNSY::SERVER_AFLOWLIB_ORG_flag",aurostd::args2flag(argv,cmds,"--server=aflow.org"));
  if(vflag.flag("APENNSY::SERVER_AFLOWLIB_ORG_flag")) {XHOST.APENNSY_USE_SERVER=FALSE;XHOST.APENNSY_USE_LIBRARY=FALSE;XHOST.APENNSY_SERVER_AFLOWLIB_ORG=TRUE;}

  vflag.flag("APENNSY::LOAD_LIB2",aurostd::args2flag(argv,cmds,"--lib2|--LIB2|--libraryX|--LIBRARYX"));
  vflag.flag("APENNSY::LOAD_LIB2U",aurostd::args2flag(argv,cmds,"--lib2u|--LIB2U|--libraryU|--LIBRARYU"));
  vflag.flag("APENNSY::LOAD_LIB2PGM",aurostd::args2flag(argv,cmds,"--lib2pgm|--LIB2PGM|--libraryPGM|--LIBRARYPGM"));
  vflag.flag("APENNSY::LOAD_LIB2X",aurostd::args2flag(argv,cmds,"--lib2x|--LIB2X|--libraryXX|--LIBRARYXX"));
  vflag.flag("APENNSY::APOOL_PUBLIC",aurostd::args2flag(argv,cmds,"--apool|--apool_public"));
  vflag.flag("APENNSY::APOOL_PRIVATE",aurostd::args2flag(argv,cmds,"--apool_private"));
  vflag.flag("APENNSY::APOOL_TEST",aurostd::args2flag(argv,cmds,"--apool_test"));
  vflag.flag("APENNSY::LOAD_ALLOY",aurostd::args2flag(argv,cmds,"--alloys|--ALLOYS|--alloy|--ALLOY"));
  // XHOST.vflag_control.flag("PRINT_MODE::EPS")=aurostd::args2flag(argv,cmds,"--EPS|--eps");

  vflag.flag("APENNSY::ENTHALPY_LIST",aurostd::args2flag(argv,cmds,"--ENERGY|--energy|--ENTHALPY|--enthalpy"));
  vflag.flag("APENNSY::DATA",aurostd::args2flag(argv,cmds,"--DATA|--data"));
  vflag.flag("APENNSY::UNCLE",aurostd::args2flag(argv,cmds,"--UNCLE|--uncle"));
  vflag.flag("APENNSY::ENTHALPY_TOT",aurostd::args2flag(argv,cmds,"--Htot|-Htot|--enthalpy_total"));
  vflag.flag("APENNSY::ENTHALPY_ATOM",aurostd::args2flag(argv,cmds,"--Hat|-Hat|--enthalpy_atom"));
  vflag.flag("APENNSY::ENTHALPY_FORMATION_ATOM",aurostd::args2flag(argv,cmds,"--Hfat|-Hfat|--enthalpy_formation_atom"));
  if(!vflag.flag("APENNSY::ENTHALPY_TOT") && 
      !vflag.flag("APENNSY::ENTHALPY_ATOM") && 
      !vflag.flag("APENNSY::ENTHALPY_FORMATION_ATOM")) 
    vflag.flag("APENNSY::ENTHALPY_TOT",TRUE);

  vflag.flag("APENNSY::WEB",aurostd::args2flag(argv,cmds,"--WEB|--web"));
  vflag.flag("APENNSY::ALL",aurostd::args2flag(argv,cmds,"--ALL|--all"));
  vflag.flag("APENNSY::FCC",aurostd::args2flag(argv,cmds,"--FCC|--fcc"));
  vflag.flag("APENNSY::BCC",aurostd::args2flag(argv,cmds,"--BCC|--bcc"));
  vflag.flag("APENNSY::HCP",aurostd::args2flag(argv,cmds,"--HCP|--hcp"));

  vflag.flag("APENNSY::PS_ENERGY_LIST",aurostd::args2flag(argv,cmds,"--PSENERGY|--psenergy"));
  vflag.flag("APENNSY::CONVEX_HULL",aurostd::args2flag(argv,cmds,"--HULLS|--hulls|--HULL|--hull"));
  vflag.flag("APENNSY::MATLAB",aurostd::args2flag(argv,cmds,"--MATLAB|--matlab"));
  vflag.flag("APENNSY::GNUPLOT",aurostd::args2flag(argv,cmds,"--GNUPLOT|--gnuplot"));
  vflag.flag("APENNSY::SMALL_CONVEX_HULL_MATLAB",aurostd::args2flag(argv,cmds,"--SHULLS|--shulls|--SHULL|--shull"));
  vflag.flag("APENNSY::HISTOGRAM_LIST",aurostd::args2flag(argv,cmds,"--HISTOGRAM|--histogram"));
  vflag.flag("APENNSY::RULES",aurostd::args2flag(argv,cmds,"--RULES|--rules"));
  vflag.flag("APENNSY::STRUCTURES",aurostd::args2flag(argv,cmds,"--STR|--str"));
  vflag.flag("APENNSY::VASPIN",aurostd::args2flag(argv,cmds,"--VASPIN|--vaspin"));
  vflag.flag("APENNSY::ORDER",aurostd::args2flag(argv,cmds,"--order|--order"));
  vflag.flag("APENNSY::INFO",aurostd::args2flag(argv,cmds,"--INFORMATION|--information|--INFO|--info"));
  vflag.flag("APENNSY::MISCIBILITY",aurostd::args2flag(argv,cmds,"--MIX|--mix|--MISCIBILITY|--miscibility"));
  vflag.flag("APENNSY::MISCIBILITY_EXPERIMENTS",aurostd::args2flag(argv,cmds,"--experiments|--EXPERIMENTS"));
  vflag.flag("APENNSY::MISCIBILITY_MIEDEMA",aurostd::args2flag(argv,cmds,"--miedema|--MIEDEMA"));
  vflag.flag("APENNSY::MISCIBILITY_HUMEROTHERY",aurostd::args2flag(argv,cmds,"--humerothery|--HUMEROTHERY"));
  vflag.flag("APENNSY::MISCIBILITY_TABLE",aurostd::args2flag(argv,cmds,"--table|--TABLE"));
  vflag.flag("APENNSY::MISCIBILITY_STATISTICS",aurostd::args2flag(argv,cmds,"--statistics|--STATISTICS|--stat|--STAT"));
  vflag.flag("APENNSY::STRUCTURE_VOLUMES",aurostd::args2flag(argv,cmds,"--volumes|--vol|--VOLUMES"));
  vflag.flag("APENNSY::REFERENCE",aurostd::args2flag(argv,cmds,"--reference|--ref|--REFERENCE|--REF|--references|--refs|--REFERENCES|--REFS"));

  vflag.flag("APENNSY::UPDATE",aurostd::args2flag(argv,cmds,"--UPDATE|--update"));
  if(vflag.flag("APENNSY::UPDATE")) {
    vflag.flag("APENNSY::ENTHALPY_LIST",TRUE);
    vflag.flag("APENNSY::CONVEX_HULL",TRUE);
    vflag.flag("APENNSY::GNUPLOT",TRUE);
  }
  vflag.flag("APENNSY::NEGLECT_STRUCTURES",aurostd::args2attachedflag(argv,cmds,"--neglect="));
  vflag.flag("APENNSY::CSWAP",aurostd::args2flag(argv,cmds,"--CSWAP|--cswap"));
  vflag.flag("APENNSY::NOCLEAN",aurostd::args2flag(argv,cmds,"--noclean"));

  // PRIORITIES FOR SNAPSHOT
  if(XHOST.vflag_control.flag("APENNSY::LATEX_SNAPSHOT")) {
    XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT",TRUE);
    XHOST.vflag_control.flag("APENNSY::LATEX_CITE",TRUE);
    XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS",TRUE);
    vflag.flag("APENNSY::ENTHALPY_LIST",TRUE);
    vflag.flag("APENNSY::CONVEX_HULL",TRUE);
    vflag.flag("APENNSY::GNUPLOT",TRUE);
  } 

  if(LDEBUG) cout << "ApennsyARGs: xscheme=" << vflag.xscheme << endl;
  if(LDEBUG) cout << "ApennsyARGs: vxscheme.size()=" << vflag.vxscheme.size() << endl;
  if(LDEBUG) cout << "ApennsyARGs: argv.size()=" << argv.size() << endl;

  return vflag.vxscheme.size();
}


// ---------------------------------------------------------------------------
int Apennsymain(vector<string> &argv,vector<string> &cmds) {
  stringstream message;
  //  bool LDEBUG=(FALSE || XHOST.DEBUG);
  if(cmds.size()) {;} //dummy load
  // std::vector<string> cmds;
  string systemS;
  // INITIALIZE ******************************************************************
  ostringstream aus;
  _aflags aflags;
  string label="",date=aflow_get_time_string_short();
  APENNSY_Parameters params;

  // *****************************************************************************
  // Introduction
  cerr << aflow::Banner("BANNER_BIG");
  // *****************************************************************************
  if(!XHOST.is_command("gnuplot")) {
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"gnuplot is missing",_RUNTIME_INIT_); //CO20200624
  } 
  if(!XHOST.is_command("latex")) {
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"latex is missing",_RUNTIME_INIT_); //CO20200624
  }
  if(!XHOST.is_command("pdflatex")) {
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"pdflatex is missing",_RUNTIME_INIT_); //CO20200624
  }
  if(!XHOST.is_command("dvips")) {
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"dvips is missing",_RUNTIME_INIT_); //CO20200624
  }
  if(!XHOST.is_command("dvipdf")) {
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"dvipdf is missing",_RUNTIME_INIT_); //CO20200624
  }
  if(!XHOST.is_command("ps2pdf")) {
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"ps2pdf is missing",_RUNTIME_INIT_); //CO20200624
  }

  // [OBSOLETE] XHOST.vflag_apennsy.clear();                                 // inside init::Init::InitMachine
  // [OBSOLETE] ApennsyARGs(argv,cmds,XHOST.vflag_apennsy);         // inside init::Init::InitMachine
  aflags.vflag=XHOST.vflag_apennsy; 

  if(aflags.vflag.flag("APENNSY::HELP")) {
    cout << init::InitGlobalObject("README_AFLOW_APENNSY_TXT");cout << aflow::Banner("BANNER_BIG");return 1;  //expressed correctly
  }
  if(aflags.vflag.vxscheme.size()==0 || argv.size()==1) {
    cout << init::InitGlobalObject("README_AFLOW_APENNSY_TXT");cout << aflow::Banner("BANNER_BIG");return 0;  //error
  }

  aflowlib::_aflowlib_entry aflowlib_server;
  vector<string> aflowlib_server_content;

  //   if(aflags.vflag.flag("APENNSY::SERVER_AFLOWLIB_ORG_flag")) {
  //     aflowlib_server.url2aflowlib("aflow.org",cout);
  //     for(uint i=0;i<aflowlib_server.vserver.size();i++) {
  //       aurostd::StringSubst(aflowlib_server.vserver.at(i),":","/");aurostd::StringSubst(aflowlib_server.vserver.at(i),"//","/");
  //       aurostd::string2tokens(aurostd::url2string(aflowlib_server.vserver.at(i)+"?aflowlib_entries"),aflowlib_server.vserverdir.at(i),",");
  //     }
  //     if(aflags.APENNSY_LIST_flag) {
  //       for(uint i=0;i<aflowlib_server.vserver.size();i++) {
  // 	cerr << "SERVER(" << i << ")=" << aflowlib_server.vserver.at(i) << endl; 
  // 	cerr << "aflowlib_server.vserverdir.at(i).size()=" << aflowlib_server.vserverdir.at(i).size() << endl;
  //       }    
  //     }
  //   }

  params.Load_ALL=TRUE;
  params.PseudopotentialNoclean=aflags.vflag.flag("APENNSY::NOCLEAN");

  if(aflags.vflag.flag("APENNSY::ALL")) aflags.APENNSY_LATTICE_flag="ALL";
  if(aflags.vflag.flag("APENNSY::FCC")) aflags.APENNSY_LATTICE_flag="FCC";
  if(aflags.vflag.flag("APENNSY::BCC")) aflags.APENNSY_LATTICE_flag="BCC";
  if(aflags.vflag.flag("APENNSY::HCP")) aflags.APENNSY_LATTICE_flag="HCP";
  if(aflags.vflag.flag("APENNSY::WEB") && aflags.APENNSY_LATTICE_flag!="FCC" && 
      aflags.APENNSY_LATTICE_flag!="BCC" && aflags.APENNSY_LATTICE_flag!="HCP")  aflags.APENNSY_LATTICE_flag="ALL";

  if(aflags.vflag.flag("APENNSY::WEB") && aflags.APENNSY_LATTICE_flag=="ALL") {params.Load_ALL=TRUE;}
  if(aflags.vflag.flag("APENNSY::WEB") && aflags.APENNSY_LATTICE_flag=="FCC") {params.Load_FCC=TRUE;params.Load_ALL=FALSE;aflags.vflag.flag("APENNSY::ALL",FALSE);}
  if(aflags.vflag.flag("APENNSY::WEB") && aflags.APENNSY_LATTICE_flag=="BCC") {params.Load_BCC=TRUE;params.Load_ALL=FALSE;aflags.vflag.flag("APENNSY::ALL",FALSE);}
  if(aflags.vflag.flag("APENNSY::WEB") && aflags.APENNSY_LATTICE_flag=="HCP") {params.Load_HCP=TRUE;params.Load_ALL=FALSE;aflags.vflag.flag("APENNSY::ALL",FALSE);}


  // PRIORITIES
  if(aflags.vflag.flag("APENNSY::WEB")) XHOST.vflag_control.flag("PRINT_MODE::PNG",TRUE);
  if(XHOST.vflag_control.flag("PRINT_MODE::PNG")) {
    XHOST.vflag_control.flag("PRINT_MODE::PDF",FALSE);
    XHOST.vflag_control.flag("PRINT_MODE::GIF",FALSE);
    XHOST.vflag_control.flag("PRINT_MODE::JPG",FALSE);
    XHOST.vflag_control.flag("PRINT_MODE::EPS",FALSE);
  }
  if(XHOST.vflag_control.flag("PRINT_MODE::PDF")) {
    XHOST.vflag_control.flag("PRINT_MODE::GIF",FALSE);
    XHOST.vflag_control.flag("PRINT_MODE::JPG",FALSE);
    XHOST.vflag_control.flag("PRINT_MODE::EPS",FALSE);
  }    
  if(XHOST.vflag_control.flag("PRINT_MODE::EPS") || XHOST.vflag_control.flag("PRINT_MODE::PDF")) {  // VECTOR
    aflags.APENNSY_GNUPLOT_FONT_str=DEFAULT_GNUPLOT_EPS_FONT;
    aflags.APENNSY_GNUPLOT_FONT_BOLD_str=DEFAULT_GNUPLOT_EPS_FONT_BOLD;
    aflags.APENNSY_GNUPLOT_FONT_ITALICS_str=DEFAULT_GNUPLOT_EPS_FONT_ITALICS;
  }
  if(XHOST.vflag_control.flag("PRINT_MODE::PNG") || XHOST.vflag_control.flag("PRINT_MODE::JPG") || XHOST.vflag_control.flag("PRINT_MODE::GIF")) {  // BITMAP
    aflags.APENNSY_GNUPLOT_FONT_str=DEFAULT_GNUPLOT_PNG_FONT;
    aflags.APENNSY_GNUPLOT_FONT_BOLD_str=DEFAULT_GNUPLOT_PNG_FONT_BOLD;
    aflags.APENNSY_GNUPLOT_FONT_ITALICS_str=DEFAULT_GNUPLOT_PNG_FONT_ITALICS;
  }

  // START
  if(aflags.vflag.flag("APENNSY::NEGLECT_STRUCTURES")) cerr << "aflags.vflag.flag(\"APENNSY::NEGLECT_STRUCTURES\")" << endl;
  aflags.APENNSY_NEGLECT_STRUCTURES_vstrs.clear();
  if(aflags.vflag.flag("APENNSY::NEGLECT_STRUCTURES")) {
    aurostd::get_itemized_vector_string_from_input(argv,"--neglect",aflags.APENNSY_NEGLECT_STRUCTURES_vstrs,",");
    for(uint i=0;i<aflags.APENNSY_NEGLECT_STRUCTURES_vstrs.size();i++)
      cerr << "APENNSY skipping structure name=\"" << aflags.APENNSY_NEGLECT_STRUCTURES_vstrs.at(i) << "\"" << endl;
  }

  // ---------------------------------------------------------------
  // Priorities
  // aflags.vflag.flag("APENNSY::ENTHALPY_LIST",TRUE);
  // if(aflags.vflag.flag("APENNSY::CONVEX_HULL")       or aflags.APENNSY_STRUCTURES or
  // aflags.vflag.flag("APENNSY::SMALL_CONVEX_HULL_MATLAB") or
  // or aflags.vflag.flag("APENNSY::VASPIN") or
  // aflags.vflag.flag("APENNSY::PS_ENERGY_LIST") or aflags.APENNSY_HISTOGRAM_LIST or
  // XHOST.vflag_control.flag("PRINT_MODE::HTML"))
  // aflags.vflag.flag("APENNSY::ENTHALPY_LIST",FALSE);
  if(aflags.vflag.flag("APENNSY::CONVEX_HULL") && !aflags.vflag.flag("APENNSY::MATLAB") && !aflags.vflag.flag("APENNSY::GNUPLOT"))
    // aflags.vflag.flag("APENNSY::MATLAB",TRUE);
    aflags.vflag.flag("APENNSY::GNUPLOT",TRUE);  //RHT
  // ---------------------------------------------------------------
  // Introductions
  // cerr << aflow::Banner("BANNER_BIG");
  // ---------------------------------------------------------------
  // Some Verbosity
  if(aflags.vflag.flag("APENNSY::CONVEX_HULL"))             cerr << "APENNSY::CONVEX_HULL" << endl;
  if(aflags.vflag.flag("APENNSY::MATLAB"))                  cerr << "APENNSY::MATLAB" << endl;
  if(aflags.vflag.flag("APENNSY::GNUPLOT"))                 cerr << "APENNSY::GNUPLOT TEST" << endl;
  if(XHOST.vflag_control.flag("APENNSY::LATEX_SNAPSHOT"))   cerr << "XHOST.vflag_control.flag(\"APENNSY::LATEX_SNAPSHOT\")" << endl;
  if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT"))     cerr << "XHOST.vflag_control.flag(\"APENNSY::LATEX_OUTPUT\")" << endl;
  if(XHOST.vflag_control.flag("APENNSY::LATEX_CITE"))              cerr << "XHOST.vflag_control.flag(\"AAPENNSY::LATEX_CITE\")" << endl;
  if(XHOST.vflag_control.flag("PRINT_MODE::HTML"))          cerr << "XHOST.vflag_control.flag(\"PRINT_MODE::HTML\"))" << endl;
  if(XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS"))    cerr << "XHOST.vflag_control.flag(\"PRINT_MODE::HYPERLINKS\"))" << endl;
  if(aflags.vflag.flag("APENNSY::ENTHALPY_LIST"))           cerr << "APENNSY::ENTHALPY_LIST: tex file" << endl;
  if(aflags.vflag.flag("APENNSY::DATA"))                    cerr << "APENNSY::DATA: file" << endl;
  if(aflags.vflag.flag("APENNSY::UNCLE"))                   cerr << "APENNSY::UNCLE: input files for UNCLE" << endl;
  if(aflags.vflag.flag("APENNSY::WEB"))                     cerr << "APENNSY::WEB: input files for WEB" << endl;
  if(aflags.APENNSY_LATTICE_flag=="ALL")                    cerr << "APENNSY_LATTICE_flag: " << aflags.APENNSY_LATTICE_flag << " structures" << endl;
  if(aflags.vflag.flag("APENNSY::HISTOGRAM_LIST"))          cerr << "APENNSY::HISTOGRAM_LIST: tex file" << endl;
  if(aflags.vflag.flag("APENNSY::PS_ENERGY_LIST"))          cerr << "PHASE SEPARATING ENERGY_LIST: tex file" << endl;
  if(aflags.vflag.flag("APENNSY::RULES"))                   cerr << "APENNSY::RULES: generating ruels for LIB2" << endl;
  if(aflags.vflag.flag("APENNSY::STRUCTURES"))              cerr << "APENNSY::STRUCTURES STRUCTURE RECOGNITION" << endl;
  if(aflags.vflag.flag("APENNSY::VASPIN"))                  apennsy_std::MakeVaspIn(FALSE);
  if(aflags.vflag.flag("APENNSY::ORDER"))                   cerr << "APENNSY::ORDER: generating order table" << endl;
  if(aflags.vflag.flag("APENNSY::INFO"))                    cerr << "APENNSY::INFO: generating info.m" << endl;
  if(aflags.vflag.flag("APENNSY::MISCIBILITY"))                  cerr << "APENNSY::MISCIBILITY: generating aflow_nomix.cpp" << endl;
  if(aflags.vflag.flag("APENNSY::MISCIBILITY_EXPERIMENTS")) cerr << "APENNSY::MISCIBILITY_EXPERIMENTS: generating aflow_mix_experiments.cpp" << endl;
  if(aflags.vflag.flag("APENNSY::MISCIBILITY_MIEDEMA"))     cerr << "APENNSY::MISCIBILITY_MIEDEMA: generating Miedema solubility info" << endl;
  if(aflags.vflag.flag("APENNSY::MISCIBILITY_HUMEROTHERY")) cerr << "APENNSY::MISCIBILITY_HUMEROTHERY: generating Hume-Rothery solid solution info" << endl;
  if(aflags.vflag.flag("APENNSY::MISCIBILITY_TABLE"))       cerr << "APENNSY::MISCIBILITY_TABLE: generating table.tex" << endl;
  if(aflags.vflag.flag("APENNSY::MISCIBILITY_STATISTICS"))  cerr << "APENNSY::MISCIBILITY_STATISTICS: generating statistics" << endl;
  if(aflags.vflag.flag("APENNSY::REFERENCE"))               cerr << "APENNSY::REFERENCE" << endl;
  if(aflags.vflag.flag("APENNSY::LOAD_ALLOY")) {
    systemS=aurostd::args2string(argv,"--alloy|--ALLOY|--alloys|--ALLOYS","xxxx");

    if(systemS=="MgMo") systemS="MoMg"; // fix for NON ALPHABETIC
    if(systemS=="MgNa") systemS="NaMg"; // fix for NON ALPHABETIC
    if(systemS=="MgNb") systemS="NbMg"; // fix for NON ALPHABETIC
    if(systemS=="MgOs") systemS="OsMg"; // fix for NON ALPHABETIC
    if(systemS=="MgPb") systemS="PbMg"; // fix for NON ALPHABETIC
    if(systemS=="MgRb") systemS="RbMg"; // fix for NON ALPHABETIC
    if(systemS=="MgRe") systemS="ReMg"; // fix for NON ALPHABETIC
    if(systemS=="MgRh") systemS="RhMg"; // fix for NON ALPHABETIC
    if(systemS=="MgRu") systemS="RuMg"; // fix for NON ALPHABETIC
    if(systemS=="MgSc") systemS="ScMg"; // fix for NON ALPHABETIC
    if(systemS=="MgSn") systemS="SnMg"; // fix for NON ALPHABETIC
    if(systemS=="MgSr") systemS="SrMg"; // fix for NON ALPHABETIC
    if(systemS=="MgTa") systemS="TaMg"; // fix for NON ALPHABETIC
    if(systemS=="MgTi") systemS="TiMg"; // fix for NON ALPHABETIC
    if(systemS=="MgV")  systemS="VMg";  // fix for NON ALPHABETIC
    if(systemS=="MgW")  systemS="WMg";  // fix for NON ALPHABETIC
    if(systemS=="MgY")  systemS="YMg";  // fix for NON ALPHABETIC
    if(systemS=="MgZn") systemS="ZnMg"; // fix for NON ALPHABETIC
    if(systemS=="MgZr") systemS="ZrMg"; // fix for NON ALPHABETIC

    // DONE, move on
    cerr << "[" << systemS << "]" << endl;
  }
  if(aflags.vflag.flag("APENNSY::RULES")) aflags.vflag.flag("APENNSY::LOAD_LIB2",TRUE); // save some typing
  // if(aflags.vflag.flag("APENNSY::ORDER")) aflags.vflag.flag("APENNSY::LOAD_LIB2",TRUE); // sove some typing
  if(aflags.vflag.flag("APENNSY::ORDER")) aflags.vflag.flag("APENNSY::LOAD_LIB2U",TRUE); // sove some typing
  // *****************************************************************************
  // cerr << argv.size() << endl;
  // FILE *in_file_pointer;static char in_file_name[1024];static char string_line[1024];
  if(!aflags.vflag.flag("APENNSY::LOAD_LIB2X") && !aflags.vflag.flag("APENNSY::APOOL_PUBLIC") && !aflags.vflag.flag("APENNSY::APOOL_PRIVATE") && !aflags.vflag.flag("APENNSY::APOOL_TEST")) {
    // ---------------------------------------------------------------
    // PARAMS Pseudo Global Parameter
    // APENNSY_Parameters params;
    params._flag_ce_cswap_=aflags.vflag.flag("APENNSY::CSWAP");
    if(aflags.vflag.flag("APENNSY::WEB")) params._flag_load_xstructures=TRUE; else params._flag_load_xstructures=FALSE;
    // ---------------------------------------------------------------
    // LOAD STUFF
    if(aflags.vflag.flag("APENNSY::LOAD_LIB2"))  params.LibLoadAlloysLIB2(aflags);  // Load LIB2
    if(aflags.vflag.flag("APENNSY::LOAD_LIB2U"))  params.LibLoadAlloysLIB2U(aflags);  // Load LIB2
    if(aflags.vflag.flag("APENNSY::LOAD_LIB2PGM"))  params.LibLoadAlloysLIB2PGM(aflags);  // Load LIB2
    if(aflags.vflag.flag("APENNSY::LOAD_ALLOY"))     params.LibLoadAlloysALLOY(systemS,aflags); // Load ALLOY*
    if(params.alloys.size()==0) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No alloys to load",_INPUT_MISSING_); //CO20200624
    }
    // Loaded now start
    cerr << params.alloysmesg.str();            // Done some verbose
    // params.LibLoadAlloys();                  // Make the alloys name vector
    params.SplitAlloySpecies();                 // Divide the species
    // for(int k=0;k<params.alloys.size();k++) cerr << "params.speciesAB.at(k)=" << params.speciesAB.at(k) <<  endl;
    params.SplitAlloyPseudoPotentials();        // Divide the pseudopotentials
    params.LoadDynamicMemory(TRUE);             // Load Dynamic Memory
    params.LoadLibrary(aflags);                 // Load the whole Library
    params.MakeRankLib();                       // bottom
    params.CheckAlloyRelaxationStructures();    // Check Relaxation of Structures in Alloys
    // for(i=1;i<=Nconcentrations[2];i++) cerr << params.C.at(i) << " " <<  NumberConcentrations[2].at(i) << endl;
    // for(int k=0;k<params.alloys.size();k++) cerr << params.alloys.at(k) << " A=" << params.speciesA.at(k) << " B=" << params.speciesB.at(k) << " pA=" << pseudosA.at(k) << " pB=" << pseudosB.at(k) << endl;

    if(0)    if(XHOST.vflag_control.flag("OSS::COUT") && aflags.vflag.flag("APENNSY::SMALL_CONVEX_HULL_MATLAB")) {
      params.FixGndStatesNamesConcentrations(TRUE);
      cout << params.MatlabGndStatesNamesConcentrations(TRUE);
    }

    if(aflags.vflag.flag("APENNSY::HISTOGRAM_LIST")) cout << params.APENNSY_HistogramList(aflags);
    if(aflags.vflag.flag("APENNSY::DATA")) {cout << params.APENNSY_Data(TRUE,aflags);cout.flush();}
    if(aflags.vflag.flag("APENNSY::UNCLE")) {cout << params.APENNSY_UNCLE(TRUE,aflags);cout.flush();}
    if(aflags.vflag.flag("APENNSY::WEB") && !aflags.vflag.flag("APENNSY::REFERENCE")) {cout << params.APENNSY_Web(TRUE,aflags);cout.flush();}
    if(XHOST.vflag_control.flag("OSS::COUT") && aflags.vflag.flag("APENNSY::ENTHALPY_LIST")) cout << params.APENNSY_EnergyList(TRUE,aflags);
    if(aflags.vflag.flag("APENNSY::REFERENCE")) {cout << params.APENNSY_Reference(TRUE,aflags);cout.flush();}
    // if(XHOST.vflag_control.flag("OSS::COUT") && aflags.vflag.flag("APENNSY::CONVEX_HULL") && aflags.vflag.flag("APENNSY::GNUPLOT"))cout << params.APENNSY_ConvexHull(TRUE,aflags,_GNUPLOT_OUTPUT_MODE_);
    // if(XHOST.vflag_control.flag("OSS::COUT") && aflags.vflag.flag("APENNSY::CONVEX_HULL") && aflags.vflag.flag("APENNSY::MATLAB"))  cout << params.APENNSY_ConvexHull(TRUE,aflags,_MATLAB_OUTPUT_MODE_);
    // if(aflags.vflag.flag("APENNSY::LOAD_ALLOY") && aflags.vflag.flag("APENNSY::CONVEX_HULL") && aflags.vflag.flag("APENNSY::GNUPLOT")) { params.APENNSY_ConvexHull(TRUE,aflags,_GNUPLOT_OUTPUT_MODE_); }
    if(XHOST.vflag_control.flag("OSS::COUT") && aflags.vflag.flag("APENNSY::STRUCTURE_VOLUMES")) {cout << params.APENNSY_StructureVolumes(TRUE,aflags);}
    if(aflags.vflag.flag("APENNSY::LOAD_ALLOY") && aflags.vflag.flag("APENNSY::STRUCTURE_VOLUMES")) { params.APENNSY_StructureVolumes(TRUE,aflags);}

    if(aflags.vflag.flag("APENNSY::PS_ENERGY_LIST")) cout << params.APENNSY_PS_EnergyList(aflags);
    if(aflags.vflag.flag("APENNSY::RULES")) {cout << params.APENNSY_Rules(aflags);cout.flush();}
    if(aflags.vflag.flag("APENNSY::STRUCTURES")) cout << params.APENNSY_Structures(aflags);
    if(aflags.vflag.flag("APENNSY::ORDER")) {cout << params.APENNSY_Order(aflags);cout.flush();}
    if(aflags.vflag.flag("APENNSY::INFO")) {cout << params.APENNSY_SimulationInformation(aflags,"TIME,CORES,TIMECORES,MEM");cout.flush();}
    if(aflags.vflag.flag("APENNSY::MISCIBILITY")) {cout << params.APENNSY_Miscibility(aflags);cout.flush();}
    if(aflags.vflag.flag("APENNSY::MISCIBILITY_EXPERIMENTS")) {cout << params.APENNSY_Miscibility_Experiments(aflags);cout.flush();}
    if(aflags.vflag.flag("APENNSY::MISCIBILITY_MIEDEMA")) {cout << params.APENNSY_Miscibility_Miedema(aflags);cout.flush();}
    if(aflags.vflag.flag("APENNSY::MISCIBILITY_HUMEROTHERY")) {cout << params.APENNSY_Miscibility_HumeRothery(aflags);cout.flush();}
    if(aflags.vflag.flag("APENNSY::MISCIBILITY_TABLE")) {cout << params.APENNSY_Miscibility_Table(aflags);cout.flush();}
    if(aflags.vflag.flag("APENNSY::MISCIBILITY_STATISTICS")) {cout << params.APENNSY_Miscibility_Statistics(aflags);cout.flush();}

    if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) {
      for(uint k=0;k<params.alloys.size();k++) {
        if(XHOST.vflag_control.flag("PRINT_MODE::PDF"))
          cout << "<A href=http://" << XHOST.AFLOW_MATERIALS_SERVER << "/AFLOW/" << params.speciesAB.at(k) << ".pdf>" << params.speciesAB.at(k) << "</A> |" << endl;
        if(XHOST.vflag_control.flag("PRINT_MODE::GIF"))
          cout << "<A href=http://" << XHOST.AFLOW_MATERIALS_SERVER << "/AFLOW/" << params.speciesAB.at(k) << ".gif>" << params.speciesAB.at(k) << "</A> |" << endl;
        if(XHOST.vflag_control.flag("PRINT_MODE::JPG"))
          cout << "<A href=http://" << XHOST.AFLOW_MATERIALS_SERVER << "/AFLOW/" << params.speciesAB.at(k) << ".jpg>" << params.speciesAB.at(k) << "</A> |" << endl;
        if(XHOST.vflag_control.flag("PRINT_MODE::PNG"))
          cout << "<A href=http://" << XHOST.AFLOW_MATERIALS_SERVER << "/AFLOW/" << params.speciesAB.at(k) << ".png>" << params.speciesAB.at(k) << "</A> |" << endl;
      }
    }

    if(!XHOST.vflag_control.flag("OSS::COUT") && aflags.vflag.flag("APENNSY::SMALL_CONVEX_HULL_MATLAB")) {
      // check for MATLAB
      aurostd::CommandRequired(DEFAULT_KBIN_MATLAB_BIN); // MATLAB MUST BE AVAILABLE
      stringstream aus;aus.clear();aus.str(std::string());
      ofstream FileMATLAB(string("apennsy.m").c_str());
      params.FixGndStatesNamesConcentrations(TRUE);
      aus << params.MatlabGndStatesNamesConcentrations(TRUE);
      FileMATLAB << aus.str() << endl;
      FileMATLAB.flush();FileMATLAB.close();
      aus.clear();aus.str(std::string());

      aus << "export DISPLAY=:0.0" << endl << DEFAULT_KBIN_MATLAB_BIN << " -r " << string("apennsy") << endl;
      aurostd::execute(aus);

      //JX aus << "export DISPLAY=:0.0" << endl << DEFAULT_KBIN_MATLAB_BIN << " -r " << string("apennsy") << endl;
      //JX aurostd::execute(aus);
      //JX aus.clear();aus.str(std::string());
      //JX // aus << "mv *eps " << init::AFLOW_Projects_Directories("LIB2") << "/EPS/HULL/" << endl;
      //JX aus << "rm -f apennsy.m" << endl;
      //JX aurostd::execute(aus);

    }

    if(!XHOST.vflag_control.flag("OSS::COUT") && aflags.vflag.flag("APENNSY::CONVEX_HULL") && aflags.vflag.flag("APENNSY::MATLAB")) {
      // check for MATLAB
      aurostd::CommandRequired(DEFAULT_KBIN_MATLAB_BIN); // MATLAB MUST BE AVAILABLE
      stringstream aus;aus.clear();aus.str(std::string());
      ofstream FileMATLAB(string("apennsy.m").c_str());
      aus << params.APENNSY_ConvexHull(TRUE,aflags,_MATLAB_OUTPUT_MODE_);
      FileMATLAB << aus.str() << endl;
      FileMATLAB.flush();FileMATLAB.close();
      aus.clear();aus.str(std::string());
      aus << "export DISPLAY=:0.0" << endl << DEFAULT_KBIN_MATLAB_BIN << " -r " << string("apennsy") << endl;
      aurostd::execute(aus);
      aus.clear();aus.str(std::string());
      // aus << "mv *eps " << init::AFLOW_Projects_Directories("LIB2") << "/EPS/HULL/" << endl;
      if(!XHOST.vflag_control.flag("KEEP::MAT")) aus << "rm -f apennsy.m" << endl;
      aurostd::execute(aus);
    }

    if(!XHOST.vflag_control.flag("OSS::COUT") && aflags.vflag.flag("APENNSY::CONVEX_HULL") && aflags.vflag.flag("APENNSY::GNUPLOT")) {
      // for(int k=0;k<params.alloys.size();k++) cerr << "params.speciesAB.at(k)=" << params.speciesAB.at(k) <<  endl;
      // check for GNUPLOT
      // aurostd::CommandRequired("gnuplot"); // GNUPLOT MUST BE AVAILABLE
      stringstream aus;aus.clear();aus.str(std::string());
      string FileGNUPLOT_title=aurostd::TmpFileCreate("gsplot.gp");
      //  string FileGNUPLOT_title="gsplot.gp";
      ofstream FileGNUPLOT(string(FileGNUPLOT_title).c_str());
      aus << params.APENNSY_ConvexHull(TRUE,aflags,_GNUPLOT_OUTPUT_MODE_);
      FileGNUPLOT << aus.str() << endl;
      FileGNUPLOT.flush();FileGNUPLOT.close();
      aus.clear();aus.str(std::string());
      aus << "export DISPLAY=:0.0" << endl << XHOST.command("gnuplot") << " " << " " << FileGNUPLOT_title << endl;
      // [OBSOLETE] aus << "rm -f str_info.dat* hull.dat*" << endl;
      // [OBSOLETE] aus << "rm -f " << FileGNUPLOT_title  << endl;
      aurostd::execute(aus);
      if(!XHOST.vflag_control.flag("KEEP::GPL")) aurostd::RemoveFile(FileGNUPLOT_title);
    }

    // [OBSOLETE]   if(!XHOST.vflag_control.flag("OSS::COUT") && aflags.vflag.flag("APENNSY_PROTOCHECK")) { // never called
    // [OBSOLETE]   stringstream aus;aus.clear();aus.str(std::string());
    // [OBSOLETE]   aus << params.APENNSY_ProtoCheck(TRUE,aflags);
    // [OBSOLETE]   aurostd::execute(aus);
    // [OBSOLETE] }

    if(!XHOST.vflag_control.flag("OSS::COUT") && aflags.vflag.flag("APENNSY::ENTHALPY_LIST")) {

      stringstream aus;aus.clear();aus.str(std::string());
      string namefile="apennsy";
      label=date;
      uint calculations_total=0;
      for(uint k=0;k<params.alloys.size();k++) calculations_total+=params.ZLibrary.at(k).size();


      if(calculations_total>999) label=aurostd::PaddedPRE(aurostd::utype2string(calculations_total),5,"0");
      else label=aurostd::PaddedPRE(aurostd::utype2string(calculations_total),3,"0");
      if(aflags.vflag.flag("APENNSY::LOAD_LIB2")) namefile=namefile+".X."+label;
      if(aflags.vflag.flag("APENNSY::LOAD_LIB2U")) namefile=namefile+".U."+label;
      if(aflags.vflag.flag("APENNSY::LOAD_LIB2PGM")) namefile=namefile+".PGM."+label;
      if(aflags.vflag.flag("APENNSY::LOAD_ALLOY")) namefile=namefile+"."+systemS+"."+label;

      ofstream FileTEX(string(namefile+".tex").c_str());
      aus << params.APENNSY_EnergyList(TRUE,aflags);
      FileTEX << aus.str() << endl;
      FileTEX.flush();FileTEX.close();
      aus.clear();aus.str(std::string());
      aus << XHOST.command("latex") << " " << namefile << ".tex | grep -v Underfull" << endl;
      if(XHOST.vflag_control.flag("APENNSY::LATEX_CITE")) aus << XHOST.command("latex") << " " << namefile << ".tex | grep -v Underfull" << endl;
      if(XHOST.vflag_control.flag("APENNSY::LATEX_CITE")) aus << XHOST.command("latex") << " " << namefile << ".tex | grep -v Underfull" << endl;
      // [OBSOLETE]  aus << XHOST.command("dvips") << " " << namefile << ".dvi" << " -o " << namefile << ".eps" << endl;
      // [OBSOLETE]  aus << XHOST.command("ps2pdf") << " " << namefile << ".eps " << namefile << ".pdf" << endl;
      aus << XHOST.command("dvipdf") << " " << namefile << ".dvi" << endl;
      aus << "rm -f ";
      if(!(XHOST.vflag_control.flag("flag_KEEP::TEX") || XHOST.vflag_control.flag("APENNSY::LATEX_SNAPSHOT"))) aus << namefile << ".tex ";
      if(!XHOST.vflag_control.flag("KEEP::EPS")) aus << namefile << ".eps ";
      if(!XHOST.vflag_control.flag("KEEP::DVI")) aus << namefile << ".dvi ";
      if(!XHOST.vflag_control.flag("KEEP::TOC")) aus << namefile << ".toc ";
      aus << namefile << ".log " << namefile << ".out " << namefile << ".aux " << endl;
      aurostd::execute(aus);
    }

  }
  if(aflags.vflag.flag("APENNSY::LOAD_LIB2X") && !aflags.vflag.flag("APENNSY::APOOL_PUBLIC") && !aflags.vflag.flag("APENNSY::APOOL_PRIVATE") && !aflags.vflag.flag("APENNSY::APOOL_TEST")) {
    // ---------------------------------------------------------------
    // PARAMS Pseudo Global Parameter
    // APENNSY_Parameters params;
    params._flag_ce_cswap_=aflags.vflag.flag("APENNSY::CSWAP");
    if(aflags.vflag.flag("APENNSY::WEB")) params._flag_load_xstructures=TRUE; else params._flag_load_xstructures=FALSE;
    stringstream aus;aus.clear();aus.str(std::string());
    // string label=aflow_get_time_string_short();
    // ---------------------------------------------------------------
    // LOAD STUFF
    // APENNSY_Parameters params;
    // params.LibLoadAlloysALLOY("Ti",aflags);
    params.LibLoadAlloysLIB2(aflags);       // it will also sink all the structures which are not part of the analysis
    cerr << "APENNSY: Loading " << params.alloys.size() << " alloys" << endl;
    cerr << params.alloysmesg.str();            // Done some verbose
    params.SplitAlloySpecies();                 // Divide the species
    params.SplitAlloyPseudoPotentials();        // Divide the pseudopotentials
    params.LoadDynamicMemory(TRUE);             // Load Dynamic Memory
    params.LoadLibrary(aflags);                       // Load the whole Library
    params.MakeRankLib();                       // bottom
    params.CheckAlloyRelaxationStructures();    // Check Relaxation of Structures in Alloys

    if(!aflags.vflag.flag("APENNSY::MATLAB") && !aflags.vflag.flag("APENNSY::GNUPLOT")) aflags.vflag.flag("APENNSY::GNUPLOT",TRUE); //default output mode
    // if(aflags.vflag.flag("APENNSY::LOAD_ALLOY") && aflags.vflag.flag("APENNSY::CONVEX_HULL") && aflags.vflag.flag("APENNSY::GNUPLOT")) { params.APENNSY_ConvexHull(TRUE,aflags,_GNUPLOT_OUTPUT_MODE_); }

    if(aflags.vflag.flag("APENNSY::GNUPLOT")) {
      // check for GNUPLOT
      // aurostd::CommandRequired("gnuplot"); // GNUPLOT MUST BE AVAILABLE
      aus << params.APENNSY_ConvexHull(TRUE,aflags,_GNUPLOT_OUTPUT_MODE_);
      stringstream aus;aus.clear();aus.str(std::string());
      string FileGNUPLOT_title=aurostd::TmpFileCreate("gsplot.gp");
      //  string FileGNUPLOT_title="gsplot.gp";
      ofstream FileGNUPLOT(string(FileGNUPLOT_title).c_str());
      aus << params.APENNSY_ConvexHull(TRUE,aflags,_GNUPLOT_OUTPUT_MODE_);
      FileGNUPLOT << aus.str() << endl;
      FileGNUPLOT.flush();FileGNUPLOT.close();
      aus.clear();aus.str(std::string());
      aus << "export DISPLAY=:0.0" << endl << XHOST.command("gnuplot") << " " << " " << FileGNUPLOT_title << endl;
      // [OBSOLETE] aus << "rm -f str_info.dat* hull.dat*" << endl;
      // [OBSOLETE] aus << "rm -f " << FileGNUPLOT_title  << endl;
      aurostd::execute(aus);
      if(!XHOST.vflag_control.flag("KEEP::GPL")) aurostd::RemoveFile(FileGNUPLOT_title);
    }

    if(aflags.vflag.flag("APENNSY::MATLAB")) {
      // check for MATLAB
      aurostd::CommandRequired(DEFAULT_KBIN_MATLAB_BIN); // MATLAB MUST BE AVAILABLE
      aus << params.APENNSY_ConvexHull(TRUE,aflags,_MATLAB_OUTPUT_MODE_);
      ofstream FileMATLAB_LIB2X(string("apennsyXX.m").c_str());
      FileMATLAB_LIB2X << aus.str() << endl;
      FileMATLAB_LIB2X.close();
      aus.clear();aus.str(std::string());
      aus << "export DISPLAY=:0.0" << endl << DEFAULT_KBIN_MATLAB_BIN << " -r " << string("apennsyXX") << endl;
      aurostd::execute(aus);
    }

    vector<APENNSY_Parameters> vparams(params.alloys.size()+1);  // let`s see if it explodes

    for(uint k=0;k<vparams.size();k++) vparams.at(k).PseudopotentialNoclean=params.PseudopotentialNoclean; // copy for cleaning names.

    for(uint k=0;k<params.alloys.size();k++) {
      vparams.at(k)._flag_ce_cswap_=aflags.vflag.flag("APENNSY::CSWAP");
      if(aflags.vflag.flag("APENNSY::WEB")) vparams.at(k)._flag_load_xstructures=TRUE; else vparams.at(k)._flag_load_xstructures=FALSE;

      label=aurostd::PaddedPRE(aurostd::utype2string(params.ZLibrary.at(k).size()),3,"0");
      string namefile="";
      if(vparams.at(k).PseudopotentialNoclean==FALSE) namefile=string("apennsy."+KBIN::VASP_PseudoPotential_CleanName(params.alloys.at(k))+"."+label);
      if(vparams.at(k).PseudopotentialNoclean==TRUE) namefile=string("apennsy."+params.alloys.at(k)+"."+label);
      // cerr << "params.alloys.at(k)=" << params.alloys.at(k) << endl;
      // cerr << "(vparams.at(k).PseudopotentialNoclean=" << vparams.at(k).PseudopotentialNoclean << endl;
      // cerr << "NAMEFILE=" << namefile << endl;

      ofstream FileENERGY_LIB2X((namefile+".tex").c_str());
      aus.clear();aus.str(std::string());
      // cerr << params.alloys.at(k) << "  A=" << params.speciesA.at(k) << "  B=" << params.speciesB.at(k) << "  pA=" << params.pseudosA.at(k) << "  pB=" << params.pseudosB.at(k) << endl;

      if(vparams.at(k).PseudopotentialNoclean==FALSE) vparams.at(k).LibLoadAlloysALLOY(KBIN::VASP_PseudoPotential_CleanName(params.alloys.at(k)),aflags);  // load the particular alloy
      if(vparams.at(k).PseudopotentialNoclean==TRUE) vparams.at(k).LibLoadAlloysALLOY(params.alloys.at(k),aflags);  // load the particular alloy

      if(vparams.at(k).alloys.size()==0) {
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No alloys to load",_INPUT_MISSING_); //CO20200624
      }
      // Loaded now start
      cerr << vparams.at(k).alloysmesg.str();            // Done some verbose
      vparams.at(k).SplitAlloySpecies();                 // Divide the species
      vparams.at(k).SplitAlloyPseudoPotentials();        // Divide the pseudopotentials
      vparams.at(k).LoadDynamicMemory(TRUE);             // Load Dynamic Memory
      vparams.at(k).LoadLibrary(aflags);                 // Load the whole Library
      vparams.at(k).MakeRankLib();                       // bottom
      vparams.at(k).CheckAlloyRelaxationStructures();    // Check Relaxation of Structures in Alloys
      aus << vparams.at(k).APENNSY_EnergyList(TRUE,aflags);
      FileENERGY_LIB2X << aus.str() << endl;
      FileENERGY_LIB2X.close();

      aus.clear();aus.str(std::string());
      // aus << "mv " << KBIN::VASP_PseudoPotential_CleanName(params.alloys.at(k)) << ".eps " << init::AFLOW_Projects_Directories("LIB2") << "/EPS/HULL/" << endl;

      aus << XHOST.command("latex") << " " << namefile << ".tex  | grep -v Overfull | grep -v cmtt | grep -v type | grep -v BCC | grep -v share " << endl;
      // [OBSOLETE] aus << XHOST.command("dvips") << " " << namefile << ".dvi" << " -o " << namefile << ".eps" << endl;
      // [OBSOLETE] aus << XHOST.command("ps2pdf") << " " << namefile << ".eps " << namefile << ".pdf" << endl;
      aus << XHOST.command("dvipdf") << " " << namefile << ".dvi" << endl;
      aus << "rm -f ";
      if(!(XHOST.vflag_control.flag("flag_KEEP::TEX") || XHOST.vflag_control.flag("APENNSY::LATEX_SNAPSHOT"))) aus << namefile << ".tex ";
      if(!XHOST.vflag_control.flag("KEEP::EPS")) aus << namefile << ".eps ";
      if(!XHOST.vflag_control.flag("KEEP::DVI")) aus << namefile << ".dvi ";
      if(!XHOST.vflag_control.flag("KEEP::TOC")) aus << namefile << ".toc ";
      aus << namefile << ".log " << namefile << ".out " << namefile << ".aux " << endl;

      aurostd::execute(aus);
    }
    aus.clear();aus.str(std::string());
    aus << "rm -f " << string("apennsyXX") << endl;
    aurostd::execute(aus);
  }
  // finished
  // ----------------------------------------------------------------------------------------------------------
  // APOOL, do the online thing PUBLIC or PRIVATE
  if(aflags.vflag.flag("APENNSY::APOOL_PUBLIC") || aflags.vflag.flag("APENNSY::APOOL_PRIVATE") || aflags.vflag.flag("APENNSY::APOOL_TEST")) {
    stringstream oss;
    // APENNSY_Parameters params;
    // params.LibLoadAlloysALLOY("Ti",aflags);
    params.LibLoadAlloysLIB2(aflags);       // it will also sink all the structures which are not part of the analysis
    cerr << "APENNSY: Loading " << params.alloys.size() << " alloys" << endl;
    cerr << params.alloysmesg.str();            // Done some verbose
    params.SplitAlloySpecies();                 // Divide the species
    params.SplitAlloyPseudoPotentials();        // Divide the pseudopotentials

    // DISCLAIMER
    oss << "<b>Library of online structural <i>ab initio</i> calculations.</b>" << endl;
    oss << "<A href=http://" << XHOST.AFLOW_WEB_SERVER << "/aflow.html><B>HELP</B></A><!P><br>" << endl;
    oss << "<font size=1pt><font COLOR=black>" << "AFLOW V" << string(AFLOW_VERSION) << " </font><br><br>" << endl;
    if(0) {
      // oss << "<br> " << endl;
      oss << "<b>DISCLAIMER for the \"binary alloy library\"</b><br>" << endl;
      oss << "Data obtained from the <b>\"<i>aflow ab-initio binary alloy library</i>\"</b> is to be used for lawful purposes only. You are welcome to download structural information from our site at no charge only for scientific use. For commercial use of the data or to access other parts of the database, please contact Prof. S. Curtarolo (SC) regarding licensing policies. SC and his collaborators make every reasonable effort to ensure the accuracy and validity of the information provided. All data and information contained herein are provided \"as is\", without any express or implied warranty. Duke University, SC and collaborators accept no liability or responsibility for any errors or omissions in the content on the site or for damages as a result of relying on information contained within this site. Data may be updated, changed, suspended or withdrawn at any time for any reason at the sole discretion of SC. <br>" << endl;
      oss << "<b>Fair Use Acknowledgement</b><br>" << endl;
      oss << "In research publications, patents, books and other published media, the use of \"<i>aflow ab-initio binary alloy library</i>\" data <b>requires</b> acknowledgement to the reference: [S. Curtarolo, D. Morgan, K. Persson, J. Rodgers, and G. Ceder, <i>Predicting Crystal Structures with Data Mining of Quantum Calculations</i>, Phys. Rev. Lett. 91, 135503 (2003)] <b>in addition</b> to the reference(s) specified inside the \"raw\" data downloaded and listed below (check the lines REFERENCE). If no REFERENCE is provided, data must be considered incomplete or not accurate and it shall be used with extreme care." << endl;
    }
    oss << "<br>" << endl;
    oss << "<!*********************************************************************************************>" << endl;
    oss << "<img border=0 width=100% height=1 src=http://" << XHOST.AFLOW_WEB_SERVER << "/auro/images/line.gif>" << endl;
    oss << "<!Note: Avaliable alloys (other alloys will follow soon)> <!br>" << endl;
    oss << "<br>" << endl;

    // PRE
    oss << "<br>" << endl;
    if(aflags.vflag.flag("APENNSY::APOOL_TEST")) oss << "<font color=green>GREEN-ALLOYS are available.</font>";
    if(aflags.vflag.flag("APENNSY::APOOL_TEST")) oss << "<font color=red>RED-ALLOYS are in progress (inquire for updates).</font>";
    if(!aflags.vflag.flag("APENNSY::APOOL_TEST")) oss << " Alloys are <font color=green>Aa</font><font color=red>Bb</font>. ";
    oss << " <i>(last update: " << TODAY << ")</i><br>" << endl;
    oss << endl;

    // if(aflags.vflag.flag("APENNSY::APOOL_PUBLIC"))  oss << "<FORM ACTION=\"http://" << XHOST.AFLOW_MATERIALS_SERVER << "/cgi-bin/awrapper.cgi\" METHOD=\"POST\" TARGET=\"_parent\">" << endl;
    // if(aflags.vflag.flag("APENNSY::APOOL_TEST"))           oss << "<FORM ACTION=\"http://" << XHOST.AFLOW_MATERIALS_SERVER << "/cgi-bin/awrapper.cgi\" METHOD=\"POST\" TARGET=\"_parent\">" << endl;
    // if(aflags.vflag.flag("APENNSY::APOOL_PRIVATE")) oss << "<FORM ACTION=\"http://" << XHOST.AFLOW_MATERIALS_SERVER << "/cgi-bin/awrapper2.cgi\" METHOD=\"POST\" TARGET=\"_parent\">" << endl;
    if(aflags.vflag.flag("APENNSY::APOOL_PUBLIC"))  oss << "<FORM ACTION=\"http://" << XHOST.AFLOW_MATERIALS_SERVER << "/php/apool.php\" METHOD=\"POST\" TARGET=\"_parent\">" << endl;
    oss << "<input type=\"hidden\" name=\"job\" value=\"awrapper_apool\">" << endl;
    oss << "<font style='font-family:courier; font-size:9pt; color: #000000'>" << endl;
    oss << "<input type=\"radio\" name=\"lattice\" value=\"all\" checked>all&nbsp;" << endl;
    oss << "<input type=\"radio\" name=\"lattice\" value=\"bcc\">bcc&nbsp;" << endl;
    oss << "<input type=\"radio\" name=\"lattice\" value=\"fcc\">fcc&nbsp;" << endl;
    oss << "<input type=\"radio\" name=\"lattice\" value=\"hcp\">hcp&nbsp;" << endl;
    // oss << "<input type=\"radio\" name=\"lattice\" value=\"bcc_uncle\">bcc (uncle) &nbsp;" << endl;
    // oss << "<input type=\"radio\" name=\"lattice\" value=\"fcc_uncle\">fcc (uncle) &nbsp;" << endl;
    // oss << "<input type=\"radio\" name=\"lattice\" value=\"hcp_uncle\">hcp (uncle) &nbsp;" << endl;
    oss << "<input type=\"radio\" name=\"lattice\" value=\"phasediagram\">phase-diagram&nbsp;" << endl;
    oss << "<br> " << endl;

    // MAKE LIST
    if(aflags.vflag.flag("APENNSY::APOOL_TEST")) {
      int i=0,split=17;
      string alloyname;
      vector<string> vrefs;
      for(uint k=0;k<params.alloys.size();k++) {
        i++;
        alloyname=params.speciesAB.at(k);
        if(alloyname.length()==3) alloyname+="&nbsp;"; // &nbsp;
        if(alloyname.length()==2) alloyname+="&nbsp;&nbsp;"; // &nbsp;
        vrefs.clear();
        SystemReferences(params.speciesAB.at(k),vrefs);
        if((vrefs.size()>0 && aflags.vflag.flag("APENNSY::APOOL_PUBLIC")) || aflags.vflag.flag("APENNSY::APOOL_PRIVATE")) {
          oss << "<input type=\"radio\" name=\"alloy\" value=\"" << params.speciesAB.at(k) << "\"";
          if(alloyname=="NbRh") oss << " checked";
          oss << "><font color=green>" << alloyname << "</font><!available>" << endl;
        } else {
          oss << "<input type=\"radio\" disabled name=\"alloy\" value=\"" << params.speciesAB.at(k) << "\" onClick=\"this.checked=false; alert('n!')\"><font color=red>" << alloyname << "</font>" << endl;
        }
        if(mod(i,split)==0) oss << "<br>" << endl;
      }
    }

    // MAKE SQUARE
    if(!aflags.vflag.flag("APENNSY::APOOL_TEST")) {
      if(aflags.vflag.flag("APENNSY::APOOL_PUBLIC") || aflags.vflag.flag("APENNSY::APOOL_PRIVATE")) {
        // int i=0;
        string alloynameAB,alloynameBA;
        vector<string> vrefsAB,vrefsBA,vspecies;
        for(uint k=0;k<params.alloys.size();k++) {
          if(params.speciesA.at(k)!="B" && params.speciesA.at(k)!="C" && params.speciesA.at(k)!="K" && params.speciesA.at(k)!="Be" && params.speciesA.at(k)!="Ba") vspecies.push_back(params.speciesA.at(k));
          if(params.speciesB.at(k)!="B" && params.speciesB.at(k)!="C" && params.speciesB.at(k)!="K" && params.speciesB.at(k)!="Be" && params.speciesB.at(k)!="Ba") vspecies.push_back(params.speciesB.at(k));
        }	
        aurostd::sort_remove_duplicates(vspecies);
        oss << "<br>" << endl;
        for(uint kB=0;kB<vspecies.size();kB++) {
          if(kB==0) {
            oss << "&nbsp;&nbsp;&nbsp;";
          } else {
            oss << "<font color=red>" << vspecies.at(kB) << "</font>" << (vspecies.at(kB).length()==1?"&nbsp;":"") << " ";//&nbsp;";  // VERTICAL COLUMS
          }
          for(uint kA=0;kA<=kB;kA++) {
            alloynameAB=vspecies.at(kA)+vspecies.at(kB);
            alloynameBA=vspecies.at(kB)+vspecies.at(kA);
            // if(alloynameAB.length()==3) alloynameAB+="&nbsp;"; // &nbsp;
            // if(alloynameAB.length()==2) alloynameAB+="&nbsp;&nbsp;"; // &nbsp;
            if(kB==kA) {
              oss << " <font color=green>" << vspecies.at(kA) << "</font>";
            } else {
              vrefsAB.clear();SystemReferences(alloynameAB,vrefsAB); //
              vrefsBA.clear();SystemReferences(alloynameBA,vrefsBA); // in the event that the alloy is not tabulated in alphabetic
              string alloyname="";
              if(AlloyAlphabeticLIBRARY(alloynameAB)==TRUE) {alloyname=alloynameAB;} else {alloyname=alloynameBA;}
              if((aflags.vflag.flag("APENNSY::APOOL_PUBLIC") && vrefsAB.size()+vrefsBA.size()>0) || aflags.vflag.flag("APENNSY::APOOL_PRIVATE")) {
                oss << "<input type=\"radio\" name=\"alloy\" value=\"" << alloyname << "\"";
                if(alloyname=="NbRh") oss << " checked";
                oss << " >";// << endl;
              } else {
                oss << "<input type=\"radio\" disabled name=\"alloy\" value=\"" << alloyname << "\" onClick=\"this.checked=false; alert('n!')\">";// << endl;
              }
            }
          }
          oss << "<br>" << endl;
        }
      }

    }
    // POST
    oss << "</font>" << endl;
    oss << " " << endl;
    oss << "<br> <br>" << endl;
    oss << "<INPUT TYPE=\"submit\" VALUE=\"Start\">" << endl;
    oss << "<INPUT TYPE=\"reset\" VALUE=\"Reset\">" << endl;
    oss << "(be patient)<p>" << endl;
    oss << "</FORM>" << endl;
    oss << " " << endl;

    oss << " " << endl;
    oss << "<img border=0 width=100% height=1 src=http://" << XHOST.AFLOW_WEB_SERVER << "/auro/images/line.gif><br><br>" << endl;

    if(1) { // REFERENCES
      oss << "<x>" << endl;
      vector<_outreach> voutreach;
      // check aflow.h:#define AFLOW_PHP_APOOL_REFERENCES
      vector<uint> vwnumber;aurostd::string2tokens(AFLOW_PHP_APOOL_REFERENCES,vwnumber,",");
      voutreach_load(voutreach,"PUBLICATIONS");
      // oss << "LOADED " << voutreach.size()-3 << " voutreach <br>" << endl;
      oss << "<b><i>References</i></b><br>" << endl;
      for(uint iarticle=0;iarticle<vwnumber.size();iarticle++) {
        for(uint i=0;i<voutreach.size();i++) {
          if(vwnumber.at(iarticle)==voutreach.at(i).wnumber) {
            // GOT THE ARTICLE
            if(iarticle+1<10)  oss << "<sup>&nbsp;&nbsp;" << iarticle+1 << "</sup> "; // LABEL
            if(iarticle+1>=10) oss << "<sup>" << iarticle+1 << "</sup> "; // LABEL
            // AUTHORS
            string authors;
            for(uint iauth=0;iauth<voutreach.at(i).vauthor.size();iauth++) {
              authors+=voutreach.at(i).vauthor.at(iauth);
              if(voutreach.at(i).vauthor.size()==2 && iauth==voutreach.at(i).vauthor.size()-2) authors+=" and ";
              if(voutreach.at(i).vauthor.size()!=2 && iauth==voutreach.at(i).vauthor.size()-2) authors+=", and ";
              if(iauth!=voutreach.at(i).vauthor.size()-2) authors+=", ";
            }
            aurostd::StringSubst(authors,"Csányi","Csanyi");
            oss << authors << endl;
            // TITLE
            string title=voutreach.at(i).title;
            title="<i>"+title+"</i>, ";
            oss << title;
            // JOURNAL with year
            string journal=voutreach.at(i).journal;
            oss << "" << journal << ". ";
            // DONE
            oss << " <br>" << endl;
          }
        }
      }
      oss << "</x>" << endl;
      oss << "<br>" << endl;
      // oss << "<img border=0 width=100% height=5 src=http://" << XHOST.AFLOW_WEB_SERVER << "/auro/images/line.gif>" << endl;
      oss << "<br>" << endl;
      // oss << "<br>" << endl;
      cout << oss.str();
    }
  }

  return 1;
}

// **********************************************************************************************
// CLASSES PRINT STUFF
// **********************************************************************************************
string Verbatim(bool VERBATIM) {
  if(VERBATIM==TRUE) {
    if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) return "\\verb|";
    if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) return "<tt>";
  }
  if(VERBATIM==FALSE) {
    if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) return "|";
    if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) return "</tt>";
  }
  return "";
}

string StrWithLink(string alloy, string structure_name) {
  stringstream oss;
  string _structure_name=structure_name;
  aurostd::StringSubst(_structure_name,"_","\\_");

  //cerr << "XHOST.vflag_control.flag(\"APENNSY::LATEX_OUTPUT\")=" << XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT") << endl;
  //cerr << "XHOST.vflag_control.flag(\"PRINT_MODE::HYPERLINKS\")=" << XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS") << endl;
  //cerr << "XHOST.vflag_control.flag(\"PRINT_MODE::HTML\")=" << XHOST.vflag_control.flag("PRINT_MODE::HTML") << endl;

  int found=0;
  if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT") && !XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS")) {
    found=1;
    oss << "" << aurostd::PaddedPRE(CleanNameStructure(structure_name),APENNSY_STR_DEPTH);
  }
  if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT") && XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS")) {
    found=2;
    oss << Verbatim(FALSE) << "{\\href{http://" << AFLOW_MATERIALS_SERVER_DEFAULT << "/AFLOWDATA/LIB2_RAW/" << alloy << "/" << structure_name << "/" << _XENTRY_ << "}{" 
      << CleanNameStructure(_structure_name) << "}}" << Verbatim(TRUE) << aurostd::PaddedPOST("",APENNSY_STR_DEPTH-1-CleanNameStructure(structure_name).length());
    //	<< aurostd::PaddedPRE(CleanNameStructure(structure_name),APENNSY_STR_DEPTH-1) << " }}" << Verbatim(TRUE);
  }
  if(XHOST.vflag_control.flag("PRINT_MODE::HTML") && !XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS")) {
    found=3;
    oss << "" << aurostd::PaddedPRE(CleanNameStructure(structure_name),APENNSY_STR_DEPTH);
  }
  if(XHOST.vflag_control.flag("PRINT_MODE::HTML") && XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS")) {
    found=4;
    oss << "" << aurostd::PaddedPRE("",APENNSY_STR_DEPTH-(CleanNameStructure(structure_name).length())) 
      << "<A href=http://" << AFLOW_MATERIALS_SERVER_DEFAULT << "/AFLOWDATA/LIB2_RAW/"
      << alloy << "/" << structure_name << "/" << _XENTRY_ << ">" << CleanNameStructure(structure_name) << "</A>";
  }
  if(found==0) {
    cerr << "StrWithLink found=" << found << endl;
    cerr << "XHOST.vflag_control.flag(\"APENNSY::LATEX_OUTPUT\")=" << XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT") << endl;
    cerr << "XHOST.vflag_control.flag(\"PRINT_MODE::HYPERLINKS\")=" << XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS") << endl;
    cerr << "XHOST.vflag_control.flag(\"PRINT_MODE::HTML\")=" << XHOST.vflag_control.flag("PRINT_MODE::HTML") << endl;
  }
  return oss.str();
}

// **********************************************************************************************
// APENNSY_Parameters::APENNSY_EnergyList
// **********************************************************************************************
string APENNSY_Parameters::APENNSY_EnergyList(bool _verbose,_aflags &aflags) {
  stringstream oss;
  double dHf,dHfd,dHft,Ts,normAA1,normAB1,normBB1,normAA2,normAB2,normBB2,errorV,errorAA,errorAB,errorBB;
  string BANNER=Verbatim(TRUE)+"************************************************************************************************************** "+Verbatim(FALSE);
  if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) BANNER+="<br>";
  if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) BANNER+="\\\\";

  vector<uint> vwnumber_global;


  uint calculations_total=0;
  for(uint k=0;k<alloys.size();k++) calculations_total+=ZLibrary.at(k).size();
  if(XHOST.vflag_control.flag("APENNSY::LATEX_SNAPSHOT")) {
    cerr << "SNAPSHOT WITH calculations_total=" << calculations_total << endl;
  }

  // ************************************************************************
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag")) {;}; // phony
  // ************************************************************************
  oss.setf(std::ios::fixed,std::ios::floatfield);
  if(_verbose) cerr << "**********************************************************" << endl;
  if(_verbose) cerr << "* MAKING ENERGY LIST                                     *" << endl;
  if(_verbose) cerr << "**********************************************************" << endl;
  if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) {
    if(!XHOST.vflag_control.flag("APENNSY::LATEX_CITE")) oss << "\\documentclass" << APENNSY_LATEX_DOC_TYPE_TWOCOLUMS << " %LATEX" << endl;
    if(XHOST.vflag_control.flag("APENNSY::LATEX_CITE"))  oss << "\\documentclass[aps,amssymb,twocolumn]{article} %LATEX" << endl;
    oss << "\\usepackage{epsfig} %LATEX" << endl;
    oss << "\\voffset -40mm %LATEX" << endl;
    if(!XHOST.vflag_control.flag("APENNSY::LATEX_CITE")) oss << "\\hoffset -40mm %LATEX" << endl;
    if(XHOST.vflag_control.flag("APENNSY::LATEX_CITE")) oss << "\\hoffset -15mm %LATEX" << endl;
    if(!XHOST.vflag_control.flag("APENNSY::LATEX_CITE")) oss << LATEX_TEXT_HEIGHT << endl; // oss << "\\textheight 288mm %LATEX" << endl;
    if(XHOST.vflag_control.flag("APENNSY::LATEX_CITE")) oss << "\\textheight 266mm %LATEX" << endl;
    oss << "\\textwidth 190mm %LATEX" << endl;
    oss << "\\usepackage{fancyvrb}     % HYPERLINKS" << endl;

    if(XHOST.vflag_control.flag("APENNSY::LATEX_CITE")) {
      oss << "\\usepackage{tocloft} %LATEX for table of contents " << endl;
      oss << "\\renewcommand{\\cftsecleader}{\\cftdotfill{\\cftdotsep}} %LATEX" << endl;
      oss << "\\renewcommand\\cftsecfont{\\small} %LATEX for tableofcontent" << endl;
      oss << "\\renewcommand\\cftsecpagefont{\\small} %LATEX for tableofcontent" << endl;
      oss << "\\usepackage{fancyhdr} %LATEX" << endl;
    }

    // if(XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS")) oss << "\\usepackage{ifthen}     % HYPERLINKS" << endl;
    if(XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS")) oss << "\\usepackage{hyperref}      % HYPERLINKS" << endl;
    if(XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS")) oss << "%\\hypersetup{linkcolor=blue} % HYPERLINKS" << endl;
    if(XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS")) oss << "\\hypersetup{colorlinks=true,linkcolor=blue,urlcolor=blue} % HYPERLINKS" << endl;
    // if(XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS")) oss << "\\def\\hyperlinks{true} % HYPERLINKS" << endl;
    // if(XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS")) oss << "%\\def\\hyperlinks{false} % HYPERLINKS" << endl;
    if(XHOST.vflag_control.flag("APENNSY::LATEX_SNAPSHOT")) {
      oss << "\\usepackage{authblk} %LATEX %SNAPSHOT" << endl;
    }

    oss << "\\begin{document} %LATEX" << endl;
    // oss << "\\title{" << alloysRAW[1] << ", Project} %LATEX" << endl;
    // oss << "\\author{S. Curtarolo:  http://" << AFLOWLIB_MATERIALS_SERVER << "} %LATEX" << endl;
    // oss << "\\date{\\today} %LATEX" << endl;
    // oss << "\\maketitle %LATEX" << endl;
    // oss << "" << APENNSY_LATEX_FONT_SMALL << " %LATEX" << endl;

    if(XHOST.vflag_control.flag("APENNSY::LATEX_SNAPSHOT")) {
      oss << "\\onecolumn %LATEX %SNAPSHOT" << endl;
      oss << "\\title{\\ \\\\ \\ \\\\ AFLOW.ORG: Snapshot of 3d-4d-5d binary intermetallics\\\\ (" << TODAY ")} %LATEX %SNAPSHOT" << endl;
      oss << "\\author[1]{Camilo E. Calderon} %LATEX %SNAPSHOT" << endl;
      oss << "\\author[1]{Ohad Levy} %LATEX %SNAPSHOT" << endl;
      oss << "\\author[2]{Gus L.W. Hart} %LATEX %SNAPSHOT" << endl;
      oss << "\\author[3,4]{Stefano Curtarolo} %LATEX %SNAPSHOT" << endl;
      oss << "\\affil[1]{\\small \\textit{Department of Mechanical Engineering and Materials Science, Duke University, Durham, NC 27708, USA}} %LATEX %SNAPSHOT" << endl;
      oss << "\\affil[2]{\\small \\textit{Department of Physics and Astronomy, Brigham Young University, Provo UT 84602, USA}} %LATEX %SNAPSHOT" << endl;
      oss << "\\affil[3]{\\small \\textit{Materials Science, Electrical Engineering, Chemistry and Physics, Duke University, Durham, NC 27708, USA}} %LATEX %SNAPSHOT" << endl;
      oss << "\\affil[4]{\\small \\textit{email: stefano@duke.edu}} %LATEX %SNAPSHOT" << endl;
      oss << "\\maketitle %LATEX %SNAPSHOT" << endl;
      oss << "\\begin{abstract} %LATEX %SNAPSHOT" << endl;
      oss << "The document includes thermodynamic information for 3d-4d-5d binary intermetallics from the \textsf{aflow.org} repository %LATEX %SNAPSHOT" << endl;
      oss << "(data retrieved on " << TODAY << " and comprising " << calculations_total << " quantum mechanical entries). %LATEX %SNAPSHOT" << endl;
      oss << "The following binary intermetallics are included:";
      for(uint k=0;k<alloys.size();k++)
        oss << " " << "\\hyperref[" << speciesAB.at(k) << "]{" << KBIN::VASP_PseudoPotential_CleanName(speciesAB.at(k)) << "}" << (k<alloys.size()-2?",":", and"); 
      oss << ". %LATEX %SNAPSHOT" << endl;
      oss << "The reader is referred to Refs. \\cite{curtarolo:art75,curtarolo:art65} for the notation. %LATEX %SNAPSHOT" << endl;
      oss << "If appropriate, each alloy contains references to published articles where accurate post-processed interpretation of the data is presented. %LATEX %SNAPSHOT" << endl;
      oss << "\\end{abstract} %LATEX %SNAPSHOT" << endl;
      oss << "\\twocolumn %LATEX %SNAPSHOT" << endl;
      oss << "\\newpage %LATEX %SNAPSHOT" << endl;    
    }
    oss << "\\renewcommand{\\baselinestretch}{0.70} %LATEX" << endl;

    if(XHOST.vflag_control.flag("APENNSY::LATEX_CITE")) {
      oss << "\\pagenumbering{arabic} %LATEX %CITE" << endl;
      //  oss << "\\pagestyle{plain} %LATEX %CITE" << endl;
      oss << "\\pagestyle{fancy} %LATEX %CITE" << endl;
      oss << "\\lhead{} \\chead{} \\rhead{} %LATEX %CITE" << endl;
      oss << "\\lfoot{{\\tiny \\textsf{\\hyperref[referencetoc]{\\underline{Contents}}}}} \\rfoot{{\\tiny \\textsf{\\href{http://www.aflow.org}{www.aflow.org}}}} \\cfoot{\\thepage} %LATEX %CITE" << endl;
      oss << "\\footskip 14pt %LATEX %CITE" << endl;
      oss << "\\label{referencetoc} %LATEX %CITE" << endl;
      oss << "\\renewcommand{\\baselinestretch}{0.50} %LATEX" << endl;
      oss << "\\tableofcontents %LATEX %CITE" << endl;
      oss << "\\newpage %LATEX %CITE" << endl;
      oss << "\\renewcommand{\\baselinestretch}{0.70} %LATEX" << endl;    
      oss << "\\onecolumn %LATEX %CITE" << endl;
    }
    oss << "" << APENNSY_LATEX_FONT_FOOTNOTESIZE << " %LATEX" << endl;
    // oss << "{\\large [" << alloysRAW[1] << "] " << Nstructures << "*" << alloys.size() << "=" << Nstructures*alloys.size() << " structures}" << endl;
    // oss << "\\cleardoublepage %LATEX" << endl;
    // oss << " %LATEX" << endl;
  }

  for(uint k=0;k<alloys.size();k++)  // if(Alloy2MiscibilityHT.at(k)!=MISCIBILITY_SYSTEM_MISCIBLE)
  {    
    if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) {
      if(!XHOST.vflag_control.flag("APENNSY::LATEX_CITE")) oss << "{\\LARGE " << Verbatim(TRUE) << speciesAB.at(k) << Verbatim(FALSE) << "}\\ \\\\" << endl;
      //	if(XHOST.vflag_control.flag("APENNSY::LATEX_CITE")) oss << "\\section{" << speciesAB.at(k) << "\\ \\cite{eqweqweq}}" << endl << "\\vskip -3mm" << endl;
      if(XHOST.vflag_control.flag("APENNSY::LATEX_CITE")) {
        vector<uint> vwnumber_local,vwnumber;
        SystemReferences(speciesAB.at(k),vwnumber_local);
        for(uint iart=0;iart<vwnumber_local.size();iart++)                   // remove aflow and aflowlib papers
          if(vwnumber_local.at(iart)!=65 && vwnumber_local.at(iart)!=75)     // remove aflow and aflowlib papers
            vwnumber.push_back(vwnumber_local.at(iart));                     // remove aflow and aflowlib papers
        oss << "\\section{\\hyperref[" << speciesAB.at(k) << "]{" << speciesAB.at(k) << "} \\ ";
        if(vwnumber.size()>0) {
          oss << "\\cite{";
          for(uint iart=0;iart<vwnumber.size();iart++)
            oss << "curtarolo:art" << vwnumber.at(iart) << (iart<vwnumber.size()-1?",":"");
          for(uint iart1=0;iart1<vwnumber.size();iart1++) {
            bool found=FALSE;
            for(uint iart2=0;iart2<vwnumber_global.size()&&!found;iart2++)
              if(vwnumber.at(iart1)==vwnumber_global.at(iart2)) found=TRUE;
            if(!found) vwnumber_global.push_back(vwnumber.at(iart1));
          }
          oss << "}"; // end cite
        }
        oss << "}";
        oss << "\\label{" << speciesAB.at(k) << "}" << endl << "\\vskip -3mm" << endl;
      }

      //	if(XHOST.vflag_control.flag("APENNSY::LATEX_CITE")) oss << "\\cite{eqweqweq}" << endl;
      if(XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS")) 
        oss << "{ \\href{http://" << AFLOW_MATERIALS_SERVER_DEFAULT << "/AFLOWDATA/LIB2_RAW/" << alloys.at(k) << "/" << _XENTRY_ << "}{[alloy in aflow.org]}}\\ " << endl;
      //	if(!XHOST.vflag_control.flag("APENNSY::LATEX_CITE")) 
      oss << "{\\verb|[Hf_max=" << enthalpy_formation_atom_max.at(k) << " (eV), ";
      oss << "Ts_max=" << aurostd::utype2string(round(entropic_temperature_max.at(k)),5) << " (K), ";
      oss << "pseudopotentials=" << alloys.at(k) << "] |} \\\\  %LATEX" << endl;
      oss << "{" << APENNSY_LATEX_FONT << " \\verb|[calcs=" << ZLibrary.at(k).size() << "] [date=" << TODAY << "]";
      oss << "[" << _APENNSY_COPYRIGHT_ << "]";
      oss << "[Aflow " << string(AFLOW_VERSION) << "  (C) " << XHOST.Copyright_Years << " Stefano Curtarolo]   |} \\\\" << endl;
      // oss << "\\section{" << speciesAB.at(k) << "}  %LATEX" << endl;
      // oss << "\\label{alloy" << speciesAB.at(k) << "}"  %LATEX << endl;
      // oss << "\\begin{verbatim}" << endl;
      oss << BANNER << endl; // " *********************************************************************************************
      oss << BANNER << endl; // " *********************************************************************************************
      oss << "\\verb|% " << alloys.at(k) << "|\\\\"; //" << alloysmesg.str() << "|";
    }
    if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) {
      //	oss << "<font face=\"Verdana, Arial, Helvetica, sans-serif\" size=12 color=black>";
      oss << "<font size=12 color=black>";
      if(!XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS")) oss << speciesAB.at(k);
      if(XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS")) oss << "<A href=http://" << AFLOW_MATERIALS_SERVER_DEFAULT << "/AFLOWDATA/LIB2_RAW/"<< alloys.at(k) << "/" << _XENTRY_ << ">" << speciesAB.at(k) << "</A>" << endl;
      oss << "</font>" << endl;
      oss << "[Hf_max=" << enthalpy_formation_atom_max.at(k) << " (eV), ";
      oss << "Ts_max=" << aurostd::utype2string(round(entropic_temperature_max.at(k)),5) << " (K), ";
      oss << "pseudopotentials=" << alloys.at(k) << "] <!br>" << endl;
      oss << "[calcs=" << ZLibrary.at(k).size() << "] [date=" << TODAY << "] <!br>";
      oss << "[" << _APENNSY_COPYRIGHT_ << "] <br>";
      oss << "[Aflow " << string(AFLOW_VERSION) << "  (C) " << XHOST.Copyright_Years << " Stefano Curtarolo] <br> " << endl;
      oss << BANNER; // " *********************************************************************************************
      oss << BANNER; // " *********************************************************************************************
      oss << Verbatim(TRUE) 
        << alloys.at(k)
        << Verbatim(FALSE) << "<br>";
      //	  oss << " " << alloys.at(k) << " " << alloysmesg.str() << "|";
    }
    for(uint j=0;j<ZConcentrations.at(k).size();j++) {
      if(ZConcentrations.at(k).at(j)>-CEPSILON && ZConcentrations.at(k).at(j)<1+0+CEPSILON) {  // [0,1]
        oss.precision(7);
        oss << Verbatim(TRUE);
        if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) oss << "<font color=#228b22>";
        oss << "1-CONCENTRATION=" << ZConcentrations.at(k).at(j);
        if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) oss << "</font>";
        oss << Verbatim(FALSE);
        if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) oss << "<br>";
        if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) oss << "\\\\";
        uint i1=0;
        for(uint i=0;i<RankLib.at(k).at(j).size();i++) {
          uint ii=RankLib.at(k).at(j).at(i);
          if(aurostd::isequal(ZLibrary.at(k).at(ii).stoich_b,ZConcentrations.at(k).at(j),CEPSILON)) {
            dHf=ZLibrary.at(k).at(ii).enthalpy_formation_atom-ZLibrary.at(k).at(RankLib.at(k).at(j).at(0)).enthalpy_formation_atom;
            // if(ZConcentrations.at(k).at(j)<CEPSILON || ZConcentrations.at(k).at(j)>1-CEPSILON) cerr << ZConcentrations.at(k).at(j) << " " << ZLibrary.at(k).at(RankLib.at(k).at(j).at(0)).enthalpy_formation_atom << endl;
            if(dHf<10.0) {
              Ts=round(ZLibrary.at(k).at(ii).entropic_temperature); if(abs(Ts)<0.01) Ts=0.0; if(Ts<0.0) Ts=0.0;
              if(i<=_APENNSY_MAX_STRUCTURES_NUMBER_) {
                oss << Verbatim(TRUE);
                if(ii==ZMINElibrary.at(k).at(j)) i1=ii;
                if(i1>0 && AlloyStructureIdentity.at(k).at(i1).at(ii)) oss << "*"; else oss << " ";
                oss << "str=" << StrWithLink(alloys.at(k),ZLibrary.at(k).at(ii).structure_name);
                //	  cerr << "k=" << k << " ii=" << ii << " alloys.at(k),ZLibrary.at(k).at(ii).structure_name=" << ZLibrary.at(k).at(ii).structure_name << endl;
                // << " (" << aurostd::PaddedPRE((int)ii,APENNSY_INT_DEPTH) << ")"; // (NNN) not needed anymore
                oss.precision(4);
                // Hf
                oss << "  Hf=";
                if(ZLibrary.at(k).at(ii).enthalpy_formation_atom>=0.0) oss << "+";
                oss << ZLibrary.at(k).at(ii).enthalpy_formation_atom;
                // dHf=
                oss.precision(4);
                oss << "  dHf=";
                oss << dHf;
                // Ts=
                oss << "  Ts=";
                oss << aurostd::PaddedPRE(aurostd::utype2string(Ts,5),5);//aurostd::PaddedPRE(aurostd::utype2string(Ts,5),6);
                // sg
                oss << "  sg=[" << aurostd::PaddedPRE(ZLibrary.at(k).at(ii).vsgroup.front(),APENNSY_LSTR_DEPTH);
                if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) oss << Verbatim(FALSE) << " $|$" << Verbatim(TRUE);
                if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) oss << " | ";
                oss << aurostd::PaddedPRE(ZLibrary.at(k).at(ii).vsgroup.back(),APENNSY_LSTR_DEPTH) << "]";
                if(ZLibrary.at(k).at(ii).vNsgroup.front()!=ZLibrary.at(k).at(ii).vNsgroup.back()) oss << "XX"; else oss << "  ";
                // oss << " type=%s",structure_name.at(ii));
                if(ZLibrary.at(k).at(ii).structure_description[1]=='n') { // oss << " type="
                  //		  oss << " " << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH);
                  if(XHOST.vflag_control.flag("PRINT_MODE::HTML"))  { oss << " " << aurostd::latex2html(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name));
                  } else { oss << " " << CleanNameStructure(ZLibrary.at(k).at(ii).structure_name); }
                } else {  // oss << " type="
                  if(XHOST.vflag_control.flag("PRINT_MODE::HTML"))  { oss << " " << aurostd::latex2html(ZLibrary.at(k).at(ii).structure_description);
                  } else { oss << " " << ZLibrary.at(k).at(ii).structure_description; }
                }
                oss.precision(1);
                oss << " (" << abs(ZLibrary.at(k).at(ii).spin_atom) << ")";
                // if(EFlibrary(i,k)<0.001) oss << "             ";
                if(ii==ZMINElibrary.at(k).at(j) && (int) ii!=ZGNDlibrary.at(k).at(j)) oss << " Hmin";
                if(ii==ZMINElibrary.at(k).at(j) && (int) ii==ZGNDlibrary.at(k).at(j)) {
                  if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) oss << "<font color=blue>";
                  oss << " Hmin-GND";	
                  if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) oss << "</font>";
                }
                if(ZLibrary.at(k).at(ii).energy_atom_relax1 < ZLibrary.at(k).at(ii).energy_atom - 0.005) {
                  oss << " [E1- E2 (" << 1000*(ZLibrary.at(k).at(ii).energy_atom_relax1-ZLibrary.at(k).at(ii).energy_atom) << ")]"; // diff
                  if(_verbose) cerr << " ** WARNING ********* ElibraryOSZICAR1 < Elibrary *********** " << alloys.at(k) << "/" << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH) << " " << ZLibrary.at(k).at(ii).energy_atom_relax1-ZLibrary.at(k).at(ii).energy_atom << endl;
                }
                //	oss << "" << "\\ifthenelse {\\equal{\\hyperlinks}{true}}{{\\href{http://" << AFLOW_MATERIALS_SERVER_DEFAULT << "}{str}}}{{str}}" << "\\begin{verbatim}";
                oss << Verbatim(FALSE) << endl;
                // if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) oss << "<br>";
                if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) oss << "\\\\";
              }	
              if(0) { // OSZICAR1 > OSZICAR2
                if((ZLibrary.at(k).at(ii).energy_atom_relax1-ZLibrary.at(k).at(ii).energy_atom)<-0.0) { // OLD STYLE
                  oss << " " << Verbatim(TRUE)
                    << "str=" << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH) << " (" << aurostd::PaddedPRE((int)ii,APENNSY_INT_DEPTH) << ")"
                    << Verbatim(FALSE) << endl;
                  oss << " " << Verbatim(TRUE)
                    << "** WARNING ********* ElibraryOSZICAR1 < Elibrary *********** "
                    << alloys.at(k) << "/" << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH) << " " 
                    << ZLibrary.at(k).at(ii).energy_atom_relax1-ZLibrary.at(k).at(ii).energy_atom << Verbatim(FALSE) << endl;
                  if(_verbose) cerr << " ** WARNING ********* ElibraryOSZICAR1 < Elibrary *********** "
                    << alloys.at(k) << "/" << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH) << " "
                      << ZLibrary.at(k).at(ii).energy_atom_relax1-ZLibrary.at(k).at(ii).energy_atom << endl;
                }
              }
              if(dHf>4.0 && dHf<900.0) {
                oss << " " << Verbatim(TRUE)
                  << "str=" << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH) << " (" << aurostd::PaddedPRE((int)ii,APENNSY_INT_DEPTH) << ")"
                  << Verbatim(FALSE) << endl;
                oss << " " << Verbatim(TRUE)
                  << alloys.at(k) << "/" << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH) << " " << dHf 
                  << Verbatim(FALSE) << endl;
                if(_verbose) cerr << " ** ERROR ********* ENERGY TOO BIG *********** "
                  << alloys.at(k) << "/" << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH) << " " << dHf << endl;
              }
            }
          }
        }
      }
    }
    oss << BANNER; // " *********************************************************************************************
    oss << BANNER; // " *********************************************************************************************
    if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) oss << Verbatim(TRUE);
    //      if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) oss << "<br>" << endl;
    oss << "RELAXATION STRUCTURE TEST STEFANO                                         " << alloys.at(k) << "  ";
    if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) oss << Verbatim(FALSE);
    if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) oss << "\\\\";
    oss << endl;
    if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) oss << Verbatim(TRUE);
    oss << "  " << aurostd::PaddedPRE("STR1",APENNSY_STR_DEPTH);
    oss << "  " << aurostd::PaddedPRE("STR2",APENNSY_STR_DEPTH);
    oss << "   " << aurostd::PaddedPOST("conc",APENNSY_DBL_DEPTH);
    oss << " " << aurostd::PaddedPRE("space_group",APENNSY_LSTR_DEPTH) << "";
    oss << "    dHf(eV)  dV%  dAA%  dAB%  dBB%    " << alloys.at(k) << "  ";
    if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) oss << Verbatim(FALSE);
    if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) oss << "\\\\";
    oss << endl;
    // for(i1=0;i1<this->ZLibrary.at(k).size();i1++)
    for(uint i1=0;i1<this->ZLibrary.at(k).size();i1++)
    { //CO20200106 - patching for auto-indenting
      //	  for(i2=i1;i2<this->ZLibrary.at(k).size();i2++)
      for(uint i2=i1+1;i2<this->ZLibrary.at(k).size();i2++)
      { //CO20200106 - patching for auto-indenting
        if(AlloyStructureIdentity.at(k).at(i1).at(i2)==TRUE) {
          if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) oss << Verbatim(TRUE);
          normAA1=std::pow((double) ZLibrary.at(k).at(i1).volume_cell,(double) 1.0/3.0);
          normAB1=std::pow((double) ZLibrary.at(k).at(i1).vgeometry.at(0),(double) 1.0/3.0);
          normBB1=std::pow((double) ZLibrary.at(k).at(i1).vgeometry.at(1),(double) 1.0/3.0);
          normAA2=std::pow((double) ZLibrary.at(k).at(i2).volume_cell,(double) 1.0/3.0);
          normAB2=std::pow((double) ZLibrary.at(k).at(i2).vgeometry.at(0),(double) 1.0/3.0);
          normBB2=std::pow((double) ZLibrary.at(k).at(i2).vgeometry.at(1),(double) 1.0/3.0);
          errorV=abs(ZLibrary.at(k).at(i1).volume_cell-ZLibrary.at(k).at(i2).volume_cell)/(ZLibrary.at(k).at(i1).volume_cell+ZLibrary.at(k).at(i2).volume_cell);
          errorAA=abs(ZLibrary.at(k).at(i1).bond_aa*normAA1-ZLibrary.at(k).at(i2).bond_aa*normAA2)/(ZLibrary.at(k).at(i1).bond_aa*normAA1+ZLibrary.at(k).at(i2).bond_aa*normAA2);
          errorAB=abs(ZLibrary.at(k).at(i1).bond_ab*normAB1-ZLibrary.at(k).at(i2).bond_ab*normAB2)/(ZLibrary.at(k).at(i1).bond_ab*normAB1+ZLibrary.at(k).at(i2).bond_ab*normAB2);
          errorBB=abs(ZLibrary.at(k).at(i1).bond_bb*normBB1-ZLibrary.at(k).at(i2).bond_bb*normBB2)/(ZLibrary.at(k).at(i1).bond_bb*normBB1+ZLibrary.at(k).at(i2).bond_bb*normBB2);
          //	if(100.0*2.0*max(errorV,errorAA,errorAB,errorBB)<1.0)
          {
            oss << "  " << StrWithLink(alloys.at(k),ZLibrary.at(k).at(i1).structure_name) // aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(i1).structure_name),APENNSY_STR_DEPTH)
              << "  " << StrWithLink(alloys.at(k),ZLibrary.at(k).at(i2).structure_name) // aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(i2).structure_name),APENNSY_STR_DEPTH)
              //  << "   " << aurostd::PaddedPOST(ZLibrary.at(k).at(i1).stoich_b,APENNSY_DBL_DEPTH);
              << "   " << aurostd::PaddedPOST(ZLibrary.at(k).at(i1).stoich_b,APENNSY_DBL_DEPTH);
            // oss << "  %3i  %3i   %5.5f",structures_number[i1],structures_number[i2],ZLibrary.at(k).at(i1).stoich_b);
            oss << " [" << aurostd::PaddedPRE(ZLibrary.at(k).at(i2).vsgroup.back(),APENNSY_LSTR_DEPTH) << "]";
            oss.precision(4);
            oss << "  " << abs(ZLibrary.at(k).at(i1).enthalpy_formation_atom-ZLibrary.at(k).at(i2).enthalpy_formation_atom);
            oss.precision(2);
            oss << "  " << 100*2*errorV;
            oss.precision(2);
            oss << "  " << 100*2*errorAA;
            oss.precision(2);
            oss << "  " << 100*2*errorAB;
            oss.precision(2);
            oss << "  " << 100*2*errorBB;
            oss << " " << aurostd::PaddedPRE(alloys.at(k),APENNSY_LSTR_DEPTH);
            if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) oss << Verbatim(FALSE);
            if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) oss << "\\\\";
            oss << endl;
          }	
        }
      }
    }
    //      if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) oss << "<br>" << endl;
    oss << BANNER; // " *********************************************************************************************
    oss << BANNER; // " *********************************************************************************************
    for(uint j=0;j<ZConcentrations.at(k).size();j++) {
      if(ZConcentrations.at(k).at(j)>CEPSILON && ZConcentrations.at(k).at(j)<1.0-CEPSILON) {  // ]0,1[
        oss << Verbatim(TRUE);
        if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) oss << "<font color=#228b22>";
        oss << "2-CONCENTRATION=" << ZConcentrations.at(k).at(j) << Verbatim(FALSE) << endl;
        if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) oss << "</font>";
        // if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) oss << "<br>";
        if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) oss << "\\\\";
        oss << Verbatim(TRUE);
        uint ii=RankLib.at(k).at(j).at(0);
        dHfd=ZLibrary.at(k).at(ii).distance_gnd;
        dHft=ZLibrary.at(k).at(ii).distance_tie;
        Ts=ZLibrary.at(k).at(ii).entropic_temperature;if(abs(Ts)<0.01) Ts=0.0; if(Ts<0.0) Ts=0.0;
        oss << " " << "str=" << StrWithLink(alloys.at(k),ZLibrary.at(k).at(ii).structure_name);
        oss << " (" << aurostd::PaddedPRE((int)ii,APENNSY_INT_DEPTH) << ")";
        oss.precision(4);
        oss << " dHfd=" << dHfd;
        oss.precision(4);
        oss << " dHft=" << dHft;
        oss << " Ts=" << aurostd::utype2string(Ts,5);//aurostd::PaddedPRE(aurostd::utype2string(Ts,5),6);
        if(dHfd>0.003) oss << "  -  2 PHASE REGION ";
        else {
          if(dHft<-0.003) oss << "  -  GNDSTATE ";
          else oss << "  -  GNDSTATE/TIE_LINE ";
        }
        oss << " " << speciesAB.at(k);
        oss << Verbatim(FALSE) << endl;
        // if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) oss << "<br>";
        if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) oss << "\\\\";
        // now the RANKS
        for(uint i=0;i<RankLib.at(k).at(j).size();i++) {
          uint ii=RankLib.at(k).at(j).at(i);
          // FIX  if(ZLibrary.at(k).at(ii).stoich_b==ZConcentrations.at(k).at(j))
          if(aurostd::isequal(ZLibrary.at(k).at(ii).stoich_b,ZConcentrations.at(k).at(j),CEPSILON))
          { //CO20200106 - patching for auto-indenting
            dHf=ZLibrary.at(k).at(ii).enthalpy_formation_atom-ZLibrary.at(k).at(RankLib.at(k).at(j).at(0)).enthalpy_formation_atom;
            if(dHf<0.015) {
              oss << Verbatim(TRUE);
              oss << " " << "str=" << StrWithLink(alloys.at(k),ZLibrary.at(k).at(ii).structure_name);	 
              oss.precision(4);
              oss << " dHf=" << dHf;
              oss.precision(3);
              oss << " Vat=" << ZLibrary.at(k).at(ii).volume_atom;
              oss.precision(4);
              oss << " bnd=[" << ZLibrary.at(k).at(ii).bond_aa;
              if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) oss << Verbatim(FALSE) << " $|$" << Verbatim(TRUE);
              if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) oss << " | ";
              oss << ZLibrary.at(k).at(ii).bond_ab;
              if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) oss << Verbatim(FALSE) << " $|$" << Verbatim(TRUE);
              if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) oss << " | ";
              oss << ZLibrary.at(k).at(ii).bond_bb << "]";
              oss << " sgF=[" << aurostd::PaddedPRE(ZLibrary.at(k).at(ii).vsgroup.back(),APENNSY_LSTR_DEPTH) << "]";
              if(ZLibrary.at(k).at(ii).vNsgroup.front()!=ZLibrary.at(k).at(ii).vNsgroup.back()) oss << "XX";else oss << "  ";
              if(ZLibrary.at(k).at(ii).structure_description[1]=='n') {
                if(XHOST.vflag_control.flag("PRINT_MODE::HTML"))  { oss << " type=" << aurostd::latex2html(ZLibrary.at(k).at(ii).structure_description);
                } else { oss << " type=" << aurostd::latex2html(aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH));}
              } else {
                if(XHOST.vflag_control.flag("PRINT_MODE::HTML"))  { oss << " type=" << aurostd::latex2html(ZLibrary.at(k).at(ii).structure_description);
                } else { oss << " type=" << ZLibrary.at(k).at(ii).structure_description; }
              }
              oss << Verbatim(FALSE) << endl;
              //		if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) oss << "<br>";
              if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) oss << "\\\\";
            }
          }
        }
      }
    }
    oss << BANNER; // " *********************************************************************************************
    oss << BANNER; // " *********************************************************************************************
    if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) {
      string label=aurostd::PaddedPRE(aurostd::utype2string<uint>(ZLibrary.at(k).size()),3,"0");
      oss << "\\newpage  %LATEX" << endl; 
      oss << "\\cleardoublepage  %LATEX" << endl;
      // oss << "\\psfig{file=" << "" << init::AFLOW_Projects_Directories("LIB2") << "/EPS/HULL/" << alloysRAW.at(k) << "/" << speciesAB.at(k) << ".eps}  %LATEX" << endl;
      // oss << "\\psfig{file=" << "" << init::AFLOW_Projects_Directories("LIB2") << "/EPS/HULL/" << speciesAB.at(k) << ".eps}  %LATEX" << endl;
      // oss << "\\psfig{file=" << speciesAB.at(k) << "."+label << ".eps}  %LATEX" << endl;
      oss << "\\psfig{file=" << speciesAB.at(k) << ".eps}  %LATEX" << endl;
      oss << "\\clearpage  %LATEX" << endl;
      oss << "\\newpage  %LATEX" << endl;
    }      
    if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) {
      oss << Verbatim(TRUE) << endl;
      oss << "str =  structurstre/prototype label in the aflow database (aflow --protos)." << endl;
      oss << "<br>" << endl; 
      oss << "Hf  =  formation enthalpy in eV with respect to the most stable pure elements configuration present in the database \"LIB1\"." << endl;
      oss << "       <i>H<sub>F</sub>(A<sub>x</sub>B<sub>1-x</sub>)</i> represents the ordering-strength of a mixture <i>A<sub>x</sub>B<sub>1-x</sub></i> against decomposition into" << endl;
      oss << "       its pure constituents, at the appropriate fractions, <i>xA</i> and <i>(1-x)B<i>." << endl;
      oss << "<br>" << endl;
      oss << "Ts =   entropic temperature in K, a qualitative description of resilience against ideal disorder. For a binary compond, <i>A<sub>x</sub>B<sub>1-x</sub></i>, " << endl;
      oss << "       T<sub>s</sub>[A<sub>x</sub>B<sub>1-x</sub>]=-&Delta;H(A<sub>x</sub>B<sub>1-x</sub>)/[k<sub>B</sub>(x&middot;log(x)+(1-x)&middot;log(1-x))]. ";
      oss <<        "For an alloy having <it>i-</it>compounds, T<sub>s</sub>=max<sub>i</sub>{T<sub>s</sub>[A<sub>x<sub>i</sub></sub>B<sub>1-x<sub>i</sub></sub>]}." << endl;
      oss << "       See <a href=http://arxiv.org/abs/1308.4357>http://arxiv.org/abs/1308.4357</a>." << endl;
      oss << "<br>" << endl;
      oss << "sg =   [ orig | relax ]  indicates the spacegroup of the starting (\"orig\") and final relaxed geometries (\"relax\"). Inside the " << endl;
      oss << "       structure directory, the \"orig\" file is POSCAR.orig or POSCAR.relax1 while the \"relax\" file is CONTCAR.relax or the " << endl;
      oss << "       CONTCAR.relaxN where N is the highest possible integer." << endl;
      oss << "<br>" << endl;
      oss << "XX     after the space group symbols, indicates that the space group changed during the calculation. This is important so that " << endl;
      oss << "       the real ground state can be found. The user should search for phases with similar energy (~few meV) to find the correct " << endl;
      oss << "       ground state and/or using geometrical means (e.g. visualization, radial distibution function, atomic environments..) to " << endl;
      oss << "       understand the true nature of the most stable configuration. DO NOT TRUST only the minimum energy criterium <it>a-priori</i>." << endl;
      oss << "<br>" << endl;
      oss << Verbatim(FALSE) << endl;
    }
  }
  if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT") && XHOST.vflag_control.flag("APENNSY::LATEX_CITE")) {
    vector<_outreach> voutreach;
    voutreach_load(voutreach,"PUBLICATIONS");
    //    oss << "\\newpage  %LATEX" << endl;
    oss << "\\twocolumn" << endl;
    oss << "%\\bibliographystyle{plain}" << endl;
    oss << "%\\bibliography{\\jobname}" << endl;
    oss << "\\begin{thebibliography}{10}  %LATEX" << endl;
    // oss << "\\addcontentsline{toc}{section}{References}" << endl;
    // oss << "  \\expandafter\\ifx\\csname urlstyle\\endcsname\\relax  %LATEX" << endl;
    // oss << "  \\providecommand{\\doi}[1]{doi:\\discretionary{}{}{}#1}\\else  %LATEX" << endl;
    // oss << "  \\providecommand{\\doi}{doi:\\discretionary{}{}{}\\begingroup  %LATEX" << endl;
    // oss << "    \\urlstyle{rm}\\Url}\\fi  %LATEX" << endl;
    // oss << "  \\providecommand{\\selectlanguage}[1]{\\relax}  %LATEX" << endl;
    // oss << "  \\providecommand{\\bibAnnote}[2]{%  %LATEX" << endl;
    // oss << "    \\begin{quotation}\\noindent\\textsc{Key:} #1\\\\  %LATEX" << endl;
    // oss << "     \\textsc{Annotation:}\\ #2\\end{quotation}}  %LATEX" << endl;
    oss << "\\bibitem{curtarolo:art75} %LATEX" << endl;
    oss << " S. Curtarolo, W. Setyawan, S. Wang, J. Xue, K. Yang, R. H. Taylor, L. J. Nelson, G. L. W. Hart, S. Sanvito, M. Buongiorno Nardelli, N. Mingo, and O. Levy, % ART75 - aflow 30767 %LATEX" << endl;
    oss << " \\textit{AFLOW.ORG: a distributed materials properties repository from high-throughput \\textit{ab initio} calculations}, % ART75 - aflow 30767 %LATEX" << endl;
    oss << " Comp. Mat. Sci. \\textbf{58}, 227-235 (2012).  % ART75 - aflow 30767 %LATEX" << endl;
    oss << " %LATEX" << endl;
    oss << "\\bibitem{curtarolo:art65} %LATEX" << endl;
    oss << " S. Curtarolo, W. Setyawan, G. L. W. Hart, M. Jahnatek, R. V. Chepulskii, R. H. Taylor, S. Wang, J. Xue, K. Yang, O. Levy, M. Mehl, H. T. Stokes, D. O. Demchenko, and D. Morgan, % ART65 - aflow 30767 %LATEX" << endl;
    oss << " \\textit{AFLOW: an automatic framework for high-throughput materials discovery},  % ART65 - aflow 30767 %LATEX" << endl;
    oss << " Comp. Mat. Sci. \\textbf{58}, 218-226 (2012).  % ART65 - aflow 30767 %LATEX" << endl;
    oss << " %LATEX" << endl;
    for(uint iart=0;iart<vwnumber_global.size();iart++) {
      oss << "\\bibitem{curtarolo:art" << vwnumber_global.at(iart) << "}" << endl;
      bool found=false;
      for(uint iart2=0;iart2<voutreach.size()&&!found;iart2++)
        if(vwnumber_global.at(iart)==voutreach.at(iart2).wnumber) {
          // [OBSOLETE]	  voutreach.at(iart2).print_mode="LATEX";
          XHOST.vflag_control.flag("PRINT_MODE::LATEX",TRUE);
          voutreach.at(iart2).vextra_html.clear();voutreach.at(iart2).vextra_latex.clear();
          oss << voutreach.at(iart2) << endl;
          found=true;
        }
    }
    oss << "\\end{thebibliography}  %LATEX" << endl;
  }
  if(XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT")) {
    oss << "\\end{document}  %LATEX" << endl;
  }
  if(_verbose) cerr << "**********************************************************" << endl;
  if(_verbose) cerr << "* MAKING ENERGY LIST                                     *" << endl;
  if(_verbose) cerr << "**********************************************************" << endl;
  // FINISHED
  // CLEAR UP THINGS
  // if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) return aurostd::latex2html(oss.str());
  return oss.str();
}

// **********************************************************************************************
// APENNSY_Parameters::APENNSY_Data
// **********************************************************************************************
string APENNSY_Parameters::APENNSY_Data(bool _verbose,_aflags &aflags) {
  stringstream oss;
  double H,Hf,Hf0,dHf;
  // ************************************************************************
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag")) {;}; // phony
  // ************************************************************************
  oss.setf(std::ios::fixed,std::ios::floatfield);
  if(_verbose) cerr << "**********************************************************" << endl;
  if(_verbose) cerr << "* MAKING DATA                                            *" << endl;
  if(_verbose) cerr << "**********************************************************" << endl;
  for(uint k=0;k<alloys.size();k++)
    // if(Alloy2MiscibilityHT.at(k)!=MISCIBILITY_SYSTEM_MISCIBLE)
  {
    oss << "# " << speciesAB.at(k) << " (" << alloys.at(k) << ")" ;
    oss << " [calcs=" << ZLibrary.at(k).size() << "] [date=" << TODAY << "]";
    // oss << " [" << _APENNSY_COPYRIGHT_ << "]";
    oss << " [aflow " << string(AFLOW_VERSION) << " (C) " << XHOST.Copyright_Years << " Stefano Curtarolo]}" << endl;
    oss << "#  alloy      str       C_B    E enthal.      Hf         dHf    spin_atom  sg0 sg2    V/atom      bondAA     bondAB     bondBB    desc " << endl;
    for(uint j=0;j<ZConcentrations.at(k).size();j++) {
      if(ZConcentrations.at(k).at(j)>-CEPSILON && ZConcentrations.at(k).at(j)<1.0+CEPSILON) { // [0,1]
        oss.precision(7);
        //	  i1=0;
        for(uint i=0;i<RankLib.at(k).at(j).size();i++) {
          uint ii=RankLib.at(k).at(j).at(i);
          if(aurostd::isequal(ZLibrary.at(k).at(ii).stoich_b,ZConcentrations.at(k).at(j),CEPSILON)) {
            H=ZLibrary.at(k).at(ii).enthalpy_atom;
            Hf=ZLibrary.at(k).at(ii).enthalpy_formation_atom;
            Hf0=ZLibrary.at(k).at(RankLib.at(k).at(j).at(0)).enthalpy_formation_atom;
            dHf=Hf-Hf0;
            oss << "  " << aurostd::PaddedPRE("("+speciesAB.at(k)+")",6);
            oss << "   " << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH);
            oss.precision(7);oss << "  ";oss << ZLibrary.at(k).at(ii).stoich_b;
            oss.precision(7);oss << "  ";if(H>=0.0) oss << " " << H; else oss << H;
            oss.precision(7);oss << " ";if(Hf>=0.0) oss << " " << Hf; else oss << Hf;
            oss.precision(7);oss << " ";if(dHf>=0.0) oss << " " << dHf; else oss << dHf;
            oss.precision(3);oss << "  " << abs(ZLibrary.at(k).at(ii).spin_atom) << " ";
            oss.precision(3);
            oss << " " << aurostd::PaddedPRE(aurostd::utype2string(ZLibrary.at(k).at(ii).vNsgroup.front()),3) << " " << aurostd::PaddedPRE(aurostd::utype2string(ZLibrary.at(k).at(ii).vNsgroup.back()),3) << " "; // SPACE GROUP
            oss.precision(7);oss << "  ";oss << ZLibrary.at(k).at(ii).volume_cell/ZLibrary.at(k).at(ii).natoms;
            oss.precision(7);oss << "  ";oss << ZLibrary.at(k).at(ii).bond_aa;
            oss.precision(7);oss << "  ";oss << ZLibrary.at(k).at(ii).bond_ab;
            oss.precision(7);oss << "  ";oss << ZLibrary.at(k).at(ii).bond_bb;
            if(ZLibrary.at(k).at(ii).structure_description[1]=='n') {
              oss << " " << aurostd::PaddedPOST("\""+CleanNameStructure(ZLibrary.at(k).at(ii).structure_name+"\""),33);
            } else {
              oss << " " << aurostd::PaddedPOST("\""+ZLibrary.at(k).at(ii).structure_description+"\"",33);
            }
            //		oss << " # " << speciesAB.at(k) << "  (data)";
            //if(ii==ZMINElibrary.at(k).at(j)) oss << " Hmin";
            //if((int) ii==ZGNDlibrary.at(k).at(j))  oss << "-GND";	
            oss << endl;
            if(dHf>4.0 && dHf<900.0) {
              oss << "str=" << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH) << " (" << aurostd::PaddedPRE((int)ii,APENNSY_INT_DEPTH) << ")" << endl;
              oss << " ** ERROR ********* ENERGY TOO BIG *********** "
                << alloys.at(k) << "/" << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH) << " " << dHf << endl;
              if(_verbose) cerr << " ** ERROR ********* ENERGY TOO BIG *********** "
                << alloys.at(k) << "/" << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH) << " " << dHf << endl;
            }
          }
        }
      }
    }
  }
  if(_verbose) cerr << "**********************************************************" << endl;
  if(_verbose) cerr << "* MAKING DATA                                            *" << endl;
  if(_verbose) cerr << "**********************************************************" << endl;
  // FINISHED
  return oss.str();
}

// **********************************************************************************************
// APENNSY_Parameters::APENNSY_UNCLE
// **********************************************************************************************
string APENNSY_Parameters::APENNSY_UNCLE(bool _verbose,_aflags &aflags) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);//TRUE;
  stringstream oss;
  xstructure str;
  // ************************************************************************
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag")) {;}; // phony
  // ************************************************************************
  // oss.setf(std::ios::fixed,std::ios::floatfield);
  // oss.setf(std::ios::fixed,std::ios::floatfield);
  if(_verbose) cerr << "**********************************************************" << endl;
  if(_verbose) cerr << "* " << aurostd::PaddedPOST("MAKING UNCLE",54," ") << " *" << endl;
  if(_verbose) cerr << "**********************************************************" << endl;

  if(aflags.APENNSY_LATTICE_flag=="ALL")                oss << "# Printing ALL lattices (--all)" << endl;
  if(aflags.APENNSY_LATTICE_flag=="FCC")                oss << "# Printing FCC lattices (--fcc)" << endl;
  if(aflags.APENNSY_LATTICE_flag=="BCC")                oss << "# Printing BCC lattices (--bcc)" << endl;
  if(aflags.APENNSY_LATTICE_flag=="HCP")                oss << "# Printing HCP lattices (--hcp)" << endl;
  if(aflags.vflag.flag("APENNSY::ENTHALPY_TOT")==TRUE)            oss << "# Printing Enthalpy_Total [eV] (--enthalpy_total)" << endl;
  if(aflags.vflag.flag("APENNSY::ENTHALPY_ATOM")==TRUE)           oss << "# Printing Enthalpy_Atom [eV] (--enthalpy_atom)" << endl;
  if(aflags.vflag.flag("APENNSY::ENTHALPY_FORMATION_ATOM")==TRUE) oss << "# Printing Enthalpy_Formation_Atom [eV] (--enthalpy_formation_atom)" << endl;

  for(uint k=0;k<alloys.size();k++)
    // if(Alloy2MiscibilityHT.at(k)!=MISCIBILITY_SYSTEM_MISCIBLE)
  {
    vector<uint> structures_lattice;
    structures_lattice.clear();
    if(aflags.APENNSY_LATTICE_flag=="ALL") {structures_lattice.clear();for(uint i=0;i<ZLibrary.at(k).size();i++) structures_lattice.push_back(i);}
    if(aflags.APENNSY_LATTICE_flag=="FCC") {structures_lattice.clear();for(uint i=0;i<ZLibrary.at(k).size();i++) if(ZLibrary.at(k).at(i).fcc==TRUE) structures_lattice.push_back(i);}
    if(aflags.APENNSY_LATTICE_flag=="BCC") {structures_lattice.clear();for(uint i=0;i<ZLibrary.at(k).size();i++) if(ZLibrary.at(k).at(i).bcc==TRUE) structures_lattice.push_back(i);}
    if(aflags.APENNSY_LATTICE_flag=="HCP") {structures_lattice.clear();for(uint i=0;i<ZLibrary.at(k).size();i++) if(ZLibrary.at(k).at(i).hcp==TRUE) structures_lattice.push_back(i);}  

    oss << "# " << speciesAB.at(k) << " (" << alloys.at(k) << ")" ;
    oss << " [calcs=" << ZLibrary.at(k).size() << "] [date=" << TODAY << "]";
    // oss << " [" << _APENNSY_COPYRIGHT_ << "]";
    oss << " [aflow " << string(AFLOW_VERSION) << " (C) " << XHOST.Copyright_Years << " Stefano Curtarolo]}" << endl;
    if(LDEBUG) cerr << "structures_lattice.size()=" << structures_lattice.size() << endl;
    for(uint j=0;j<structures_lattice.size();j++) {
      if(LDEBUG) cerr << "j=" << j << endl;
      uint ii=structures_lattice.at(j);
      if(LDEBUG) cerr << "ii=" << ii << endl;  
      if(LDEBUG) cerr << "ZLibrary.at(k).at(ii).vstr.size()=" << ZLibrary.at(k).at(ii).vstr.size() << endl;
      {// && ii<40) { //[CO20200106 - close bracket for indenting]}
        cerr << "Loading " << alloys.at(k) << "/" << ZLibrary.at(k).at(ii).structure_name << endl;
        deque<string> _vspecies;_vspecies.push_back(speciesA.at(k));_vspecies.push_back(speciesB.at(k));
        str=aflowlib::PrototypeLibraries(cerr,aurostd::RemoveCharacter(ZLibrary.at(k).at(ii).structure_name,'/'),"",_vspecies,LIBRARY_MODE_HTQC);
        str.SetCoordinates(_COORDS_CARTESIAN_);
        str.SetVolume(str.atoms.size());
        if(str.num_each_type.size()==1) {
          // if(LDEBUG)
          cerr << "DEBUG pure system correction" << endl;
          //	    cerr << str.species.size() << endl;
          if(ZLibrary.at(k).at(ii).pureA==TRUE) {str.num_each_type.push_back(0);str.comp_each_type.push_back(0);str.species.push_back("X");}
          if(ZLibrary.at(k).at(ii).pureB==TRUE) {str.num_each_type.push_front(0);str.comp_each_type.push_front(0);str.species.push_front("X");}
        }
        oss << str.PrintUNCLE();
        oss.precision(16);
        // if(aflags.vflag.flag("APENNSY::ENTHALPY_TOT")) oss << "# Energy Total" << endl << E*str.atoms.size() << endl;
        // if(aflags.vflag.flag("APENNSY::ENTHALPY_ATOM"))  oss << "# Energy Atom" << endl << E << endl;
        // if(aflags.vflag.flag("APENNSY::ENTHALPY_FORMATION_ATOM"))  oss << "# Energy Atom" << endl << E << endl;
        if(aflags.vflag.flag("APENNSY::ENTHALPY_TOT"))
          oss << "# Enthalpy_Total [eV] (--enthalpy_total)" << endl << ZLibrary.at(k).at(ii).enthalpy_cell << endl;
        if(aflags.vflag.flag("APENNSY::ENTHALPY_ATOM"))  
          oss << "# Enthalpy_Atom [eV] (--enthalpy_atom)" << endl << ZLibrary.at(k).at(ii).enthalpy_atom << endl;
        if(aflags.vflag.flag("APENNSY::ENTHALPY_FORMATION_ATOM")) 
          oss << "# Enthalpy_Formation_Atom [eV] (--enthalpy_formation_atom)" << endl << ZLibrary.at(k).at(ii).enthalpy_formation_atom << endl;
        //	  oss << "  " << aurostd::PaddedPRE("("+speciesAB.at(k)+")",6);
        // oss << "   " << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH);
        // oss.precision(7);oss << "  ";oss << ZLibrary.at(k).at(ii).stoich_b;
        // oss.precision(7);oss << "  ";if(E>=0.0) oss << " " << E; else oss << E;
        // oss.precision(7);oss << " ";if(Hf>=0.0) oss << " " << Hf; else oss << Hf;
        // oss.precision(3);oss << "  " << abs(ZLibrary.at(k).at(ii).spin_atom) << " ";
        // oss << endl;
        if(ZLibrary.at(k).at(ii).enthalpy_formation_atom>100.0) {
          cerr << "str=" << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH) << " (" << aurostd::PaddedPRE((int)ii,APENNSY_INT_DEPTH) << ")" << endl;
          cerr << " ** ERROR ********* ENERGY TOO BIG *********** "
            << alloys.at(k) << "/" << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH) << " " << ZLibrary.at(k).at(ii).enthalpy_formation_atom << endl;
          if(_verbose) cerr << " ** ERROR ********* ENERGY TOO BIG *********** "
            << alloys.at(k) << "/" << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH) << " " << ZLibrary.at(k).at(ii).enthalpy_formation_atom << endl;
        }
      } // end of HTQC
    }
  }
  if(_verbose) cerr << "**********************************************************" << endl;
  if(_verbose) cerr << "* " << aurostd::PaddedPOST("MAKING UNCLE",54," ") << " *" << endl;
  if(_verbose) cerr << "**********************************************************" << endl;
  // FINISHED
  return oss.str();
}

// **********************************************************************************************
// APENNSY_Parameters::APENNSY_Web
// **********************************************************************************************
string APENNSY_Parameters::APENNSY_Web(bool _verbose,_aflags &aflags) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream oss;
  uint _precision_=14; //was 16 SC 10 DM
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(_precision_);
  double H_cell,H_atom,Hf_atom,spin_atom,volume_atom,stoich_a,stoich_b;
  // ************************************************************************
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag")) {;}; // phony
  // ************************************************************************
  // oss.setf(std::ios::fixed,std::ios::floatfield);
  // oss.setf(std::ios::fixed,std::ios::floatfield);
  if(_verbose) cerr << "**********************************************************" << endl;
  if(_verbose) cerr << "* " << aurostd::PaddedPOST("MAKING WEB",54," ") << " *" << endl;
  if(_verbose) cerr << "**********************************************************" << endl;

  for(uint k=0;k<alloys.size();k++) {
    xstructure strPRE,strMID,strPOST;
    vector<uint> structures_lattice;
    structures_lattice.clear();
    if(aflags.APENNSY_LATTICE_flag=="ALL") {structures_lattice.clear();for(uint i=0;i<ZLibrary.at(k).size();i++) structures_lattice.push_back(i);}
    if(aflags.APENNSY_LATTICE_flag=="FCC") {structures_lattice.clear();for(uint i=0;i<ZLibrary.at(k).size();i++) if(ZLibrary.at(k).at(i).fcc==TRUE) structures_lattice.push_back(i);}
    if(aflags.APENNSY_LATTICE_flag=="BCC") {structures_lattice.clear();for(uint i=0;i<ZLibrary.at(k).size();i++) if(ZLibrary.at(k).at(i).bcc==TRUE) structures_lattice.push_back(i);}
    if(aflags.APENNSY_LATTICE_flag=="HCP") {structures_lattice.clear();for(uint i=0;i<ZLibrary.at(k).size();i++) if(ZLibrary.at(k).at(i).hcp==TRUE) structures_lattice.push_back(i);}  
    // if(aflags.vflag.flag("APENNSY::ENTHALPY_TOT")) oss << "# Printing Energy_Total (--energy_total)" << endl;
    // if(aflags.vflag.flag("APENNSY::ENTHALPY_ATOM"))  oss << "# Printing Enthalpy_Atom (--enthalpy_atom)" << endl;

    // if(Alloy2MiscibilityHT.at(k)!=MISCIBILITY_SYSTEM_MISCIBLE)
    // oss << "[aflow " << string(AFLOW_VERSION) << " (C) " << XHOST.Copyright_Years << " Stefano Curtarolo]" << endl;
    string strfile1="";
    if(PseudopotentialNoclean==FALSE) strfile1="["+KBIN::VASP_PseudoPotential_CleanName(alloys.at(k))+"] ";
    if(PseudopotentialNoclean==TRUE) strfile1="["+(alloys.at(k))+"] ";
    // INTRO
    if(aflags.APENNSY_LATTICE_flag=="ALL") oss << strfile1 << "# WEB file for ALL " << speciesAB.at(k) << " (" << alloys.at(k) << ")";
    if(aflags.APENNSY_LATTICE_flag=="FCC") oss << strfile1 << "# WEB file for FCC " << speciesAB.at(k) << " (" << alloys.at(k) << ")";
    if(aflags.APENNSY_LATTICE_flag=="BCC") oss << strfile1 << "# WEB file for BCC " << speciesAB.at(k) << " (" << alloys.at(k) << ")";
    if(aflags.APENNSY_LATTICE_flag=="HCP") oss << strfile1 << "# WEB file for HCP " << speciesAB.at(k) << " (" << alloys.at(k) << ")";
    oss << " [calcs=" << ZLibrary.at(k).size()-1 << "] [date=" << TODAY << "]";  // must -1 to remove "/3/" !!
    // oss << " [" << _APENNSY_COPYRIGHT_ << "]";
    oss << " [aflow " << string(AFLOW_VERSION) << " (C) " << XHOST.Copyright_Years << " Stefano Curtarolo]" << endl;
    // MAKE REFERENCES
    vector<string> vrefs;
    SystemReferences(alloys.at(k),vrefs);
    for(uint i=0;i<vrefs.size();i++) oss << strfile1 << "REFERENCE: " << vrefs.at(i) << endl;
    // PROTOTYPES
    oss << strfile1 << "PROTOTYPES: http://" << XHOST.AFLOW_MATERIALS_SERVER << "/AFLOW/proto.pdf" << endl;
    // MAKE STRUCTURES
    if(LDEBUG) cerr << "structures_lattice.size()=" << structures_lattice.size() << endl;
    for(uint j=0;j<structures_lattice.size();j++) {
      if(LDEBUG) cerr << "j=" << j << endl;
      uint ii=structures_lattice.at(j);
      if(LDEBUG) cerr << "ii=" << ii << endl;  
      if(LDEBUG) cerr << "ZLibrary.at(k).at(ii).vstr.size()=" << ZLibrary.at(k).at(ii).vstr.size() << endl;
      // HTQC
      {// && ii<40) { //[CO20200106 - close bracket for indenting]}
        strPRE=ZLibrary.at(k).at(ii).vstr.front();strPRE.BringInCell();
        strMID=ZLibrary.at(k).at(ii).vstr.at(1);strMID.BringInCell();
        strPOST=ZLibrary.at(k).at(ii).vstr.back();strPOST.BringInCell();
        //	cerr << "Loading " << alloys.at(k) << "/" << ZLibrary.at(k).at(ii).structure_name << endl;
        if(strPRE.num_each_type.size()==1 && 0) {
          // if(LDEBUG)
          cerr << "DEBUG pure system correction" << endl;
          if(ZLibrary.at(k).at(ii).pureA==TRUE) {strPRE.num_each_type.push_back(0);strPRE.comp_each_type.push_back(0);strPRE.species.push_back("X");}
          if(ZLibrary.at(k).at(ii).pureB==TRUE) {strPRE.num_each_type.push_front(0);strPRE.comp_each_type.push_front(0);strPRE.species.push_front("X");}
        }
        if(LDEBUG) cerr << "APENNSY_Web: 1" << endl;
        string strfile2="";
        stringstream straus;vector<string> tokens;
        if(PseudopotentialNoclean==FALSE) strfile2="["+KBIN::VASP_PseudoPotential_CleanName(alloys.at(k))+"/"+CleanNameStructure(ZLibrary.at(k).at(ii).structure_name)+"] ";
        if(PseudopotentialNoclean==TRUE) strfile2="["+(alloys.at(k))+"/"+CleanNameStructure(ZLibrary.at(k).at(ii).structure_name)+"] ";
        oss << strfile2 << "# **************************************************************************************************************" << endl;
        oss << strfile2 << "# ----------------------------------------------------------------------- " << endl;
        if(PseudopotentialNoclean==FALSE) oss << strfile2 << KBIN::VASP_PseudoPotential_CleanName(alloys.at(k)) << "/" << CleanNameStructure(ZLibrary.at(k).at(ii).structure_name) << endl;
        if(PseudopotentialNoclean==TRUE) oss << strfile2 << (alloys.at(k)) << "/" << CleanNameStructure(ZLibrary.at(k).at(ii).structure_name) << endl;
        oss << strfile2 << "# ---- Structure PRE ---------------------------------------------------- " << endl;
        straus.str(std::string());straus << strPRE;aurostd::string2tokens(straus.str(),tokens,"\n");for(uint i=0;i<tokens.size();i++) oss << strfile2 << tokens.at(i) << endl;
        // 	oss << strfile2 << "# Structure MID ---------------------------------------------------- " << endl;
        // 	straus.str(std::string());straus << strMID;aurostd::string2tokens(straus.str(),tokens,"\n");for(uint i=0;i<tokens.size();i++) oss << strfile2 << tokens.at(i) << endl;
        oss << strfile2 << "# ---- Structure POST --------------------------------------------------- " << endl;
        straus.str(std::string());straus << strPOST;aurostd::string2tokens(straus.str(),tokens,"\n");for(uint i=0;i<tokens.size();i++) oss << strfile2 << tokens.at(i) << endl;
        oss << strfile2 << "# ---- DATA ------------------------------------------------------------- " << endl;
        oss.precision(14);
        if(LDEBUG) cerr << "APENNSY_Web: 2" << endl;
        H_cell=ZLibrary.at(k).at(ii).enthalpy_cell;
        H_atom=ZLibrary.at(k).at(ii).enthalpy_atom;
        Hf_atom=ZLibrary.at(k).at(ii).enthalpy_formation_atom;
        spin_atom=ZLibrary.at(k).at(ii).spin_atom;
        volume_atom=ZLibrary.at(k).at(ii).volume_atom;
        stoich_a=ZLibrary.at(k).at(ii).stoich_a;
        stoich_b=ZLibrary.at(k).at(ii).stoich_b;
        oss << strfile2 << " " << (H_cell  >=0.0?" ":"") << H_cell  << "  # H [eV] (VASP) " << endl;
        oss << strfile2 << " " << (H_atom  >=0.0?" ":"") << H_atom  << "   # H/at [eV] (VASP) " << endl;
        oss << strfile2 << " " << (Hf_atom  >=0.0?" ":"") << Hf_atom  << "   # Hf_atom [eV] (VASP) " << endl;
        oss << strfile2 << " " << (spin_atom>=0.0?" ":"") << spin_atom << "   # Mom/at " << endl;
        oss << strfile2 << " " << (volume_atom<10.0?" ":"") <<  volume_atom << "   # Volume/at " << endl;
        oss << strfile2 << " " << " "              << stoich_a << "   # Ca " << endl;
        oss << strfile2 << " " << " "              << stoich_b << "   # Cb " << endl;
        oss << strfile2 << " " << aurostd::PaddedPRE(ZLibrary.at(k).at(ii).vsgroup.front(),15) << "     # space group PRE " << endl;
        //oss << strfile2 << " " << aurostd::PaddedPRE(ZLibrary.at(k).at(ii).vsgroup.at(1),15) << "     # space group MID " << endl;
        oss << strfile2 << " " << aurostd::PaddedPRE(ZLibrary.at(k).at(ii).vsgroup.back(),15) << "     # space group POST ";
        if(ZLibrary.at(k).at(ii).vsgroup.front()!=ZLibrary.at(k).at(ii).vsgroup.back()) oss << " XX ";
        oss << endl;
        oss << strfile2 << "# ---- URL -------------------------------------------------------------- " << endl;
        oss << strfile2 << " " << "http://" << AFLOW_MATERIALS_SERVER_DEFAULT << "/AFLOWDATA/LIB2_RAW/" << alloys.at(k) << "/" << ZLibrary.at(k).at(ii).structure_name << "/" << _XENTRY_;
        oss << endl;
        if(XHOST.vflag_control.flag("PRINT_MODE::HTML") && XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS")) {
          oss << strfile2 << "# ---- HYPERLINK -------------------------------------------------------- " << endl;
          oss << strfile2 << " link=<A href=" << "http://" << AFLOW_MATERIALS_SERVER_DEFAULT << "/AFLOWDATA/LIB2_RAW/" << alloys.at(k) << "/" << ZLibrary.at(k).at(ii).structure_name << "/" << _XENTRY_ << ">" << alloys.at(k) << "/" << ZLibrary.at(k).at(ii).structure_name << "</A>" << endl;
        }
        oss << strfile2 << "# ---- aflowlib.out ----------------------------------------------- " << endl;
        oss << strfile2 << " " << ZLibrary.at(k).at(ii).entry;
        oss << endl;
        if(LDEBUG) cerr << "APENNSY_Web: 3" << endl;
        // <<  "|" << aurostd::PaddedPRE(ZLibrary.at(k).at(ii).sgroupPOST,APENNSY_LSTR_DEPTH) << "]";
        //	    if(ZLibrary.at(k).at(ii).NsgroupPRE!=ZLibrary.at(k).at(ii).NsgroupPOST) oss << "XX"; else oss << "  ";

        if(Hf_atom>100) {
          cerr << "str=" << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH) << " (" << aurostd::PaddedPRE((int)ii,APENNSY_INT_DEPTH) << ")" << endl;
          cerr << " ** ERROR ********* ENERGY TOO BIG *********** "
            << alloys.at(k) << "/" << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH) << " " << Hf_atom << endl;
          if(_verbose) cerr << " ** ERROR ********* ENERGY TOO BIG *********** "
            << alloys.at(k) << "/" << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH) << " " << Hf_atom << endl;
        }
        oss << strfile2 << "# **************************************************************************************************************" << endl;

      } // end of HTQC

    } // structure
    // PROTOTYPES
    oss << strfile1 << "PROTOTYPES: http://" << XHOST.AFLOW_MATERIALS_SERVER << "/awrapper.html" << endl;
    // MAKE REFERENCES
    for(uint i=0;i<vrefs.size();i++) oss << strfile1 << "REFERENCE: " << vrefs.at(i) << endl;
    // INTRO
    if(aflags.APENNSY_LATTICE_flag=="ALL") oss << strfile1 << "# WEB file for ALL " << speciesAB.at(k) << " (" << alloys.at(k) << ")";
    if(aflags.APENNSY_LATTICE_flag=="FCC") oss << strfile1 << "# WEB file for FCC " << speciesAB.at(k) << " (" << alloys.at(k) << ")";
    if(aflags.APENNSY_LATTICE_flag=="BCC") oss << strfile1 << "# WEB file for BCC " << speciesAB.at(k) << " (" << alloys.at(k) << ")";
    if(aflags.APENNSY_LATTICE_flag=="HCP") oss << strfile1 << "# WEB file for HCP " << speciesAB.at(k) << " (" << alloys.at(k) << ")";
    oss << " [calcs=" << ZLibrary.at(k).size()-1 << "] [date=" << TODAY << "]";  // must -1 to remove "/3/" !!
    // oss << " [" << _APENNSY_COPYRIGHT_ << "]";
    oss << " [aflow " << string(AFLOW_VERSION) << " (C) " << XHOST.Copyright_Years << " Stefano Curtarolo]" << endl;
    // done
  } // alloy
  if(_verbose) cerr << "**********************************************************" << endl;
  if(_verbose) cerr << "* " << aurostd::PaddedPOST("MAKING WEB",54," ") << " *" << endl;
  if(_verbose) cerr << "**********************************************************" << endl;
  // FINISHED
  return oss.str();
}

// **********************************************************************************************
// APENNSY_Parameters::APENNSY_Reference
// **********************************************************************************************
string APENNSY_Parameters::APENNSY_Reference(bool _verbose,_aflags &aflags) {
  stringstream oss;
  // ************************************************************************
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag")) {;}; // phony
  // ************************************************************************
  // oss.setf(std::ios::fixed,std::ios::floatfield);
  // oss.setf(std::ios::fixed,std::ios::floatfield);
  if(_verbose) cerr << "**********************************************************" << endl;
  if(_verbose) cerr << "* " << aurostd::PaddedPOST("MAKING REFERENCE",54," ") << " *" << endl;
  if(_verbose) cerr << "**********************************************************" << endl;

  for(uint k=0;k<alloys.size();k++) {
    string strfile1="",strbr="";
    if(PseudopotentialNoclean==FALSE) strfile1="["+KBIN::VASP_PseudoPotential_CleanName(alloys.at(k))+"] ";
    if(PseudopotentialNoclean==TRUE) strfile1="["+(alloys.at(k))+"] ";
    if(aflags.vflag.flag("APENNSY::WEB")) strbr="<br>";
    oss << "[aflow " << string(AFLOW_VERSION) << " (C) " << XHOST.Copyright_Years << " Stefano Curtarolo]" << strbr << endl;
    // MAKE REFERENCES
    vector<string> vrefs;
    SystemReferences(alloys.at(k),vrefs,TRUE);
    for(uint i=0;i<vrefs.size();i++) {
      //      vrefs.at(i).AUTHORS_ETAL=TRUE;
      oss << strfile1 << "REF: " << vrefs.at(i) << strbr << endl;
    }
  } // alloy
  if(_verbose) cerr << "**********************************************************" << endl;
  if(_verbose) cerr << "* " << aurostd::PaddedPOST("MAKING WEB",54," ") << " *" << endl;
  if(_verbose) cerr << "**********************************************************" << endl;
  // FINISHED
  return oss.str();
}

// **********************************************************************************************
// APENNSY_Parameters::APENNSY_StructureVolumes
// **********************************************************************************************
string APENNSY_Parameters::APENNSY_StructureVolumes(bool _verbose,_aflags &aflags) {
  stringstream oss;
  if(_verbose) cerr << "*******CREATING VOLUMES PLOT********" << endl;
  double minE=0,maxE=0,minV=0,minC=0,V_A=0,V_B=0;
  bool xbool;
  // ************************************************************************
  oss.setf(std::ios::fixed,std::ios::floatfield);
  for(uint k=0;k<alloys.size();k++) {
    // common
    minE=0.0;maxE=0.0;minV=0.0;minC=0.0;
    for(uint i=1;i<ZConcentrations.at(k).size();i++) {
      if(ZGNDlibrary.at(k).at(i)>=0) {

        if(minE>ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom)
        {
          minE=ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom;
          minV=ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).volume_atom;//for Vergard analysis
          minC=ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).stoich_b;
        }
        if(maxE<ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom) maxE=ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom;
      }
    }
    V_A=ZLibrary.at(k)[ZGNDlibrary.at(k).at(1)].volume_atom; //for Vergard analysis
    V_B=ZLibrary.at(k)[ZGNDlibrary.at(k).at(ZConcentrations.at(k).size()-1)].volume_atom;
    // cout << "minimum energy" << minE << endl << "minimum volume" << minV << endl << minC << endl;

    minE-=0.02;maxE=0.05;
    ofstream gnu;
    if(_verbose) cerr << "compiling hull data...";
    // [OBSOLETE] string hull_dat_file="./hull.dat."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+"."+aurostd::utype2string(k);  //CO20200502 - threadID
    string hull_dat_file=aurostd::TmpFileCreate("hull_dat_"+aurostd::utype2string(k));
    gnu.open(hull_dat_file.c_str());
    gnu.setf(ios_base::fixed, ios_base::floatfield);
    gnu.precision(2);

    //Ground state concentrations
    gnu << "GROUNDSTATE INFORMATION:" << endl;
    gnu << setw(18) << "conc." << setw(18) << "H_" << setw(18) << "vol." << setw(18) << "S_time" << setw(18) << "Vergard_bi" << setw(18) << "Vergard_ideal" << setw(8) << "HT_name" << setw(18) << "name" << endl;
    gnu << "*********************************************************************************************" << endl;
    double AvVOL=0;int count=0;
    for(uint i=1;i<ZConcentrations.at(k).size();i++) {
      if(ZGNDlibrary.at(k).at(i)>=0) {
        gnu << setw(18) << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).stoich_b;
        if(ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom>=0) {
          gnu << setw(18) << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom;
          AvVOL=AvVOL+ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).volume_atom;
          count++;
        }
        else {
          gnu << setw(18) << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom;
          AvVOL=AvVOL+ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).volume_atom;
          count++;
        }
        gnu << setw(18) << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).volume_atom;
        //HERE
        gnu << setw(18) << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).calculation_time << " ";


        if(ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).stoich_b<minC)//Vergard analysis
          gnu << setw(18) << ((minV-V_A)/(minC))*(ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).stoich_b)+V_A;// << endl;
        else
          gnu << setw(18) << ((minV-V_B)/(minC-1))*(ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).stoich_b-1)+V_B;// << endl;

        gnu << setw(18) << (V_B-V_A)*(ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).stoich_b)+V_A;

        gnu << setw(8) << CleanNameStructure(ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).structure_name);      	
        if((ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).structure_description)=="nnn") {
          gnu << setw(18) << CleanNameStructure(ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).structure_name) << endl;
        }
        else {
          gnu << setw(18) << (ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).structure_description) << endl;
        }

        //gnu << "      " << (ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).structure_description) << endl;

      }
    }

    AvVOL=abs(AvVOL/count);

    gnu.close();
    if(_verbose) {cerr << "done" << endl;
      cerr << "compiling structure data...";}
    ofstream gnu1;
    //    [OBSOLETE] string str_info_dat_file="./str_info.dat."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+"."+aurostd::utype2string(k); //CO20200502 - threadID
    string str_info_dat_file=aurostd::TmpFileCreate("str_info_dat_"+aurostd::utype2string(k));

    gnu1.open(str_info_dat_file.c_str());
    gnu1 << setw(18) << "conc." << setw(18) << "H_" << setw(18) << "vol." << setw(18) << "S_time" << setw(18) << "HT_name" << "      " << "name" << endl;
    gnu1 << "*********************************************************************************************" << endl;
    // cerr << "compiling structure data...2";
    // for(i=1;i<ZConcentrations.at(k).size();i++)
    for(int i=this->ZLibrary.at(k).size()-1;i>=0;i--) {
      // cout << i << endl;
      xbool=TRUE;
      for(uint j=1;j<ZConcentrations.at(k).size();j++)
        if(i==(int) ZGNDlibrary.at(k).at(j))
          xbool=FALSE;
      if(xbool) {
        gnu1 << setw(18) << ZLibrary.at(k).at(i).stoich_b;
        gnu1 << setw(18) << ZLibrary.at(k).at(i).enthalpy_formation_atom;
        gnu1 << setw(18) << ZLibrary.at(k).at(i).volume_atom;
        //HERE
        gnu1 << setw(18) << ZLibrary.at(k).at(i).calculation_time;
        gnu1 << setw(18) << CleanNameStructure(ZLibrary.at(k).at(i).structure_name);
        if((ZLibrary.at(k).at(i).structure_description)=="nnn") {
          gnu1 << "      " << CleanNameStructure(ZLibrary.at(k).at(i).structure_name) << endl;
        } else {
          gnu1 << "      " << (ZLibrary.at(k).at(i).structure_description) << endl;
        }
      }
    }
    gnu1.close();
    if(_verbose) cerr << "done" << endl;
    //GNUPLOT SCRIPT EPS/PDF
    cerr << "preparing gnuplot script..." << endl;
    string syscommand1="head -4 "+hull_dat_file+ "> temp1"+XHOST.ostrPID.str()+XHOST.ostrTID.str()+" ;tail -1 "+ hull_dat_file + " > temp2"+XHOST.ostrPID.str()+XHOST.ostrTID.str()+";cat "+ "temp1"+XHOST.ostrPID.str()+XHOST.ostrTID.str()+ " temp2"+XHOST.ostrPID.str()+XHOST.ostrTID.str()+">vergard.dat."+XHOST.ostrPID.str()+XHOST.ostrTID.str()+";rm -rf temp1"+XHOST.ostrPID.str()+XHOST.ostrTID.str()+";rm -rf temp2"+XHOST.ostrPID.str()+XHOST.ostrTID.str(); //CO20200502 - threadID
    aurostd::execute(syscommand1);
    oss << "set terminal postscript landscape color enhanced '" << aflags.APENNSY_GNUPLOT_FONT_str << "' 11" << endl;  // EPS
    // oss << "set terminal pdf landscape color enhanced '" << aflags.APENNSY_GNUPLOT_FONT_str << "' 11" << endl;  // PDF
    oss << "set output 'gfig_vol_" << speciesAB.at(k) << ".eps'" << endl;//make name automatic  // EPS 
    //  oss << "set output 'gfig_vol_" << speciesAB.at(k) << ".pdf'" << endl;//make name automatic  // PDF
    oss << "set title '[" << speciesAB.at(k) << "]" << "   opar=" << alloy_order.at(k)
      << "   calcs=" << ZLibrary.at(k).size()
      << "   date: " << TODAY
      << "   [" << _APENNSY_COPYRIGHT_ << "]' font '" << aflags.APENNSY_GNUPLOT_FONT_str << ", 11'" << endl;
    oss << "set xlabel '" << speciesB.at(k) << "' font '" << aflags.APENNSY_GNUPLOT_FONT_str << ",13'" << endl;
    oss << "set ylabel 'vol./atom' font '" << aflags.APENNSY_GNUPLOT_FONT_str << ", 13'" << endl;
    oss << "#set label '" << speciesA.at(k) << "' at screen .1, screen .02" << endl;//make label automatic
    oss << "#set label '" << speciesB.at(k) << "' at screen .96, screen .02" << endl;//make label automatic
    oss << "set xrange[0:1]" << endl;
    oss << "set yrange[" << minV*.8 << ":" << AvVOL*1.5 << "]" << endl;
    oss << "yoffset=-" << AvVOL/20 << endl;//horizontal and vertical offsets to make plot nice
    oss << "xoffset=0" << endl;
    //oss << "plot '" << hull_dat_file << "' u 1:3 w linespoints pt 1 lw 1.3 ps 2 lt rgb 'blue' t '','' u ($1+($1/2-.5)*xoffset):($3+yoffset/1.5):6 with labels font '" << aflags.APENNSY_GNUPLOT_FONT_BOLD_str << ",12' lt 2 t '','" << str_info_dat_file << "' u 1:3 w points pt 2 lt rgb 'red' ps .8 t '','' u ($1+($1/2-.5)*xoffset):($3-yoffset/20):6 with labels font '" << aflags.APENNSY_GNUPLOT_FONT_str << ",4' lt 2 t '','vergard.dat' u 1:3 w linespoints pt 3 lw 1.3 ps 2 lt rgb 'green' t ''" << endl;
    oss << "plot '" << hull_dat_file << "' u 1:3 w linespoints pt 1 lw 1.3 ps 2 lt rgb 'blue' t '','' u ($1+($1/2-.5)*xoffset):($3+yoffset/1.5):8 t '' with labels font '" << aflags.APENNSY_GNUPLOT_FONT_BOLD_str << ",12' lt 2 ,'' u 1:5 w linespoints pt 3 lw 2 ps 1 lt rgb 'brown' t '','' u 1:6 w linespoints pt 3 lw 2 ps 0 lt rgb 'orange' t '','" << str_info_dat_file << "' u 1:3 w points pt 2 lt rgb 'black' ps .5 t ''" << endl;
    if(!XHOST.vflag_control.flag("KEEP::GPL")) oss << "system `" << "rm -f " << hull_dat_file << " " <<  str_info_dat_file << "`" << endl;  // DELETE hull_dat_file,str_info_dat_file

    //echo -e "plot '-' w lines lw 5 \n 6 6 \n 7 7 \n e" | gnuplot -persist
    // CANT DO IN THIS WAY SINCE IT IS NOT EXTE$NSIVE LOOK CONVEX HULL AND HOW THE OSS IS COUT-ED AND USED OUTSIDE THE ROUTINE
    // if(_verbose) cerr << "running...";
    // stringstream aus;aus.clear();aus.str(std::string());
    // aus << XHOST.command("gnuplot") << " " << " gsplot.gp." << XHOST.ostrPID.str() << "." << XHOST.ostrTID.str() << ";rm -rf " << hull_dat_file << ";rm -rf " << str_info_dat_file << ";rm -rf gsplot.gp." << XHOST.ostrPID.str() << "." << XHOST.ostrTID.str(); //CO20200502 - threadID
    // aurostd::execute(aus);
    // if(_verbose) cerr << "done" << endl;
  }

  return oss.str();
}

// **********************************************************************************************
// APENNSY_Parameters::APENNSY_ConvexHull
// **********************************************************************************************
string APENNSY_Parameters::APENNSY_ConvexHull(bool _verbose,_aflags &aflags,uint mode) {
  stringstream oss;
  double x,y,minE=0,maxE=0;
  bool xbool;
  // ************************************************************************
  oss.setf(std::ios::fixed,std::ios::floatfield);
  if(mode==_MATLAB_OUTPUT_MODE_) oss << "% *********************************************************************************** " << endl;
  if(mode==_MATLAB_OUTPUT_MODE_) oss << "% *********************************************************************************** " << endl;
  if(_verbose) cerr << "**********************************************************" << endl;
  if(_verbose) cerr << "* MAKING GROUNDSTATES PLOTS                              *" << endl;
  if(_verbose) cerr << "**********************************************************" << endl;

  for(uint k=0;k<alloys.size();k++) {
    // common
    minE=0.0;maxE=0.0;
    for(uint i=1;i<ZConcentrations.at(k).size();i++) {
      if(ZGNDlibrary.at(k).at(i)>=0) {
        if(minE>ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom) minE=ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom;
        if(maxE<ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom) maxE=ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom;
      }
    }
    minE-=0.02;maxE=0.05;

    if(mode==_MATLAB_OUTPUT_MODE_) {
      // oss << "disp('*********************************************************************');   " << "  % " << alloys.at(k) << endl;
      // oss << "disp('*    GNDstate analyzer -  Stefano Curtarolo                         *');   " << "  % " << alloys.at(k) << endl;
      // oss << "disp('*********************************************************************');   " << "  % " << alloys.at(k) << endl;
      // oss << "disp(' Generated = " << TODAY << "');   " << "  % " << alloys.at(k) << endl;
      oss << "disp('" << alloys.at(k) << "')" << "  % " << alloys.at(k) << endl;
      oss << "                                                                                 " << "  % " << alloys.at(k) << endl;
      oss << "FONTDIM=16;                                                                      " << "  % " << alloys.at(k) << endl;
      oss << "COLOR=[0 0 0];                                                                   " << "  % " << alloys.at(k) << endl;
      oss << "LINEWIDTH=0.5;                                                                   " << "  % " << alloys.at(k) << endl;
      oss << "PAPERPOSITION=[1 1 5 4];                                                         " << "  % " << alloys.at(k) << endl;
      oss << "PAPERSIZE=[8.5 11];                                                              " << "  % " << alloys.at(k) << endl;
      oss << "POSITION=[420 540 660 320];                                                      " << "  % " << alloys.at(k) << endl;
      oss << "AXISWIDTH=1;                                                                     " << "  % " << alloys.at(k) << endl;
      oss << "%AXISPOSITION=[0.136646 0.164179 0.774845 0.78678];                              " << "  % " << alloys.at(k) << endl;
      oss << "AXISPOSITION=[0.136646 0.164179 0.754845 0.72678];                               " << "  % " << alloys.at(k) << endl;
      oss << "%AXISPOSITION=[0.136646 0.164179 0.4557455 0.4754795];                           " << "  % " << alloys.at(k) << endl;
      oss << "                                                                                 " << "  % " << alloys.at(k) << endl;
      oss << "% ****************************************************************************** " << "  % " << alloys.at(k) << endl;
      oss << "% ****************************************************************************** " << "  % " << alloys.at(k) << endl;
      oss << "% CONVEX_HULL_MATLAB CODE START                                                  " << "  % " << alloys.at(k) << endl;
      oss << "% " << alloys.at(k) << "  % " << alloys.at(k) << endl;
      // oss << "cd /home/auro/work/AFLOW3" << "  % " << alloys.at(k) << endl;
      // oss << "clear;clf;%figure;" << endl;
      ///////////////////////// oss << "figure;" << "  % " << alloys.at(k) << endl;
      oss << "Cb=[";
      for(uint i=1;i<ZConcentrations.at(k).size();i++)
        if(ZGNDlibrary.at(k).at(i)>=0) {
          oss << "  " << ZConcentrations.at(k).at(i);
          // if(_verbose) cerr  << "%f  %4i    %f",ZConcentrations.at(k).at(i),structures_number.at(ZGNDlibrary.at(k).at(i)),EFlibrary(ZGNDlibrary.at(k).at(i),k));
          // oss << "  % " << alloys.at(k) << endl;
        }
      oss << "];";oss << "  % " << alloys.at(k) << endl;
      oss << "EF=[";
      for(uint i=1;i<ZConcentrations.at(k).size();i++)
        if(ZGNDlibrary.at(k).at(i)>=0) {
          if(ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom>=0)
            oss << "  " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom;
          else
            oss << " " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom;
          //	if(_verbose) cerr  << "%f  %4i    %f ",ZConcentrations.at(k).at(i),structures_number.at(ZGNDlibrary.at(k).at(i)),ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom);
          // oss << "  % " << alloys.at(k) << endl;	
        }
      oss << "];";oss << "  % " << alloys.at(k) << endl;
      oss << "%L=[";
      for(uint i=1;i<ZConcentrations.at(k).size();i++)
        if(ZGNDlibrary.at(k).at(i)>=0) {
          oss << "    " << CleanNameStructure(ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).structure_name);
          // if(_verbose) cerr  << "%f  %4i    %f ",ZConcentrations.at(k).at(i),structures_number.at(ZGNDlibrary.at(k).at(i)),ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom);
          // oss << "  % " << alloys.at(k) << endl;	
        }
      oss << "];";oss << "  % " << alloys.at(k) << endl;
      // for(uint i=1;i<ZConcentrations.at(k).size();i++)
      // if(ZGNDlibrary.at(k).at(i)>=0) {
      // // aus_string << "text(%f,%f,'%i');",ZConcentrations.at(k).at(i),ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom,structures_number.at(ZGNDlibrary.at(k).at(i)));
      // oss << "  % " << alloys.at(k) << endl;	
      // }

      oss << "subplot(1,1,1);";oss << "  % " << alloys.at(k) << endl;
      oss << "plot(Cb ,EF ,'b+-','LineWidth',0.5*LINEWIDTH);hold on;axis([0 1 " << minE << " " << maxE << "]);";oss << "  % " << alloys.at(k) << endl;
      // for(uint i=0;i<this->ZLibrary.at(k).size();i++)
      for(int i=this->ZLibrary.at(k).size()-1;i>=0;i--)
      { //CO20200106 - patching for auto-indenting
        xbool=TRUE;
        for(uint j=1;j<ZConcentrations.at(k).size();j++)
          if(i==(int) ZGNDlibrary.at(k).at(j)) xbool=FALSE;
        if(xbool) {
          oss << "plot(" << ZLibrary.at(k).at(i).stoich_b << "," << ZLibrary.at(k).at(i).enthalpy_formation_atom << ",'rx','LineWidth',0.2*LINEWIDTH);hold on;";
          if(!XHOST.QUIET) {
            if(ZLibrary.at(k).at(i).structure_description[1]=='n') {
              oss << "text(" << ZLibrary.at(k).at(i).stoich_b+0.005-0.02*_isodd(i) << "," << ZLibrary.at(k).at(i).enthalpy_formation_atom << ",'" << CleanNameStructure(ZLibrary.at(k).at(i).structure_name) << "','fontsize',0.2*FONTDIM);";
            } else {
              oss << "text(" << ZLibrary.at(k).at(i).stoich_b+0.005-0.02*_isodd(i) << "," << ZLibrary.at(k).at(i).enthalpy_formation_atom << ",'" << ZLibrary.at(k).at(i).structure_description << "','fontsize',0.2*FONTDIM);";
            }
          }
          oss << "  % " << alloys.at(k) << endl;
        }
      }
      // GROUNDTSTATES
      for(uint i=1;i<ZConcentrations.at(k).size();i++)
        if(ZGNDlibrary.at(k).at(i)>=0) {
          x=ZConcentrations.at(k).at(i);y=ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom;
          // oss_positions_symbols_colors(x,y,z);
          oss << "plot(" << x << "," << y << ",'b+','LineWidth',LINEWIDTH);hold on;";
          if(!XHOST.QUIET || aurostd::isequal(ZConcentrations.at(k).at(i),0.0,0.001) || aurostd::isequal(ZConcentrations.at(k).at(i),1.0,0.001)) {
            if(ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).structure_description[1]=='n') {
              oss << "text(" << x+0.005 << "," << y-0.007 << ",'" << CleanNameStructure(ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).structure_name) << "','fontsize',0.5*FONTDIM);"; // structures_number.at(ZGNDlibrary.at(k).at(i)));
            } else {
              oss << "text(" << x+0.005 << "," << y-0.007 << ",'" << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).structure_description << "','fontsize',0.5*FONTDIM);";
            }
          }
          oss << "  % " << alloys.at(k) << endl;
        }
      if(XHOST.QUIET) {
        //	oss << "text(0.17," << (1.3*maxE+minE)/2 << ",'AFLOW  stefano.curtarolo@duke.edu','fontsize',1.4*FONTDIM);hold on;" << "  % " << alloys.at(k) << endl;
        oss << "text(0.40," << (1.3*maxE+minE)/2 << ",'" << speciesAB.at(k) << " (aflow)','fontsize',1.8*FONTDIM);hold on;" << "  % " << alloys.at(k) << endl;
      }

      oss << "drawnow;";oss << "  % " << alloys.at(k) << endl;

      // oss << "title('" << speciesAB.at(k) << "       calcs=" << ZLibrary.at(k).size() << "/" << this->ZLibrary.at(k).size() << "    date:" << TODAY << "   [" << _APENNSY_COPYRIGHT_ << "]','fontsize',0.8*FONTDIM);";oss << "  % " << alloys.at(k) << endl;
      oss.setf(std::ios::fixed,std::ios::floatfield);oss.precision(6);
      string title=aurostd::string2latex(speciesAB.at(k));
      oss << "title('" << title
        << "       opar=" << alloy_order.at(k)
        << "   calcs=" << ZLibrary.at(k).size()
        << "   date:" << TODAY
        << "   [" << _APENNSY_COPYRIGHT_ << "]','fontsize',0.8*FONTDIM);";
      oss << "  % " << alloys.at(k) << endl;

      oss << "xlabel('" << speciesB.at(k) << "','fontsize',FONTDIM);";
      oss << "ylabel('eV/atom','fontsize',FONTDIM);";oss << "  % " << alloys.at(k) << endl;

      oss << "orient landscape;";oss << "  % " << alloys.at(k) << endl;
      oss << "print -depsc " << speciesAB.at(k) << ".eps " << "  % " << alloys.at(k) << endl;
      if(XHOST.vflag_control.flag("PRINT_MODE::PDF")) oss << "system('ps2pdf " << speciesAB.at(k) << ".eps " << speciesAB.at(k) << ".pdf " << "');   % " << alloys.at(k) << endl;
      if(XHOST.vflag_control.flag("PRINT_MODE::GIF")) oss << "system('convert -rotate 90 " << speciesAB.at(k) << ".eps " << speciesAB.at(k) << ".gif " << "');   % " << alloys.at(k) << endl;
      if(XHOST.vflag_control.flag("PRINT_MODE::JPG")) oss << "system('convert -rotate 90 " << speciesAB.at(k) << ".eps " << speciesAB.at(k) << ".jpg " << "');   % " << alloys.at(k) << endl;
      if(XHOST.vflag_control.flag("PRINT_MODE::PNG")) oss << "system('convert -rotate 90 " << speciesAB.at(k) << ".eps " << speciesAB.at(k) << ".png " << "');   % " << alloys.at(k) << endl;
      oss << "drawnow;" << "  % " << alloys.at(k) << endl;
      ////////////////// oss << "close;" << "  % " << alloys.at(k) << endl;
      oss << "clf('reset');" << "  % " << alloys.at(k) << endl;
      oss << "% CONVEX_HULL_MATLAB CODE END                                         " << "  % " << alloys.at(k) << endl;
      oss << "% *********************************************************************************** " << "  % " << alloys.at(k) << endl;
    }

    if(mode==_GNUPLOT_OUTPUT_MODE_) {
      //WORKING HERE (RHT)

      ofstream gnu;
      if(_verbose) cerr << speciesAB.at(k) << endl << "compiling hull data...";
      //  [OBSOLETE] string hull_dat_file="./hull.dat."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+"."+aurostd::utype2string(k); //CO20200502 - threadID
      string hull_dat_file=aurostd::TmpFileCreate("hull_dat_"+aurostd::utype2string(k));
      gnu.open(hull_dat_file.c_str());
      gnu.setf(ios_base::fixed, ios_base::floatfield);
      gnu.precision(8);

      //Ground state concentrations
      gnu << "GROUNDSTATE INFORMATION:" << endl;
      gnu << setw(18) << "conc." << setw(18) << "H_" << setw(18) << "vol." << setw(18) << "S_time" << setw(18) << "HT_name" << "      " << "name" << endl;
      gnu << "*********************************************************************************************" << endl;
      double AvEn=0;int count=0;
      for(uint i=1;i<ZConcentrations.at(k).size();i++) {
        if(ZGNDlibrary.at(k).at(i)>=0) {
          gnu << " " << setw(18) << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).stoich_b;
          if(ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom>=0) {
            gnu << " " <<  setw(18) << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom;
            AvEn=AvEn+ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom;
            count++;
          } else {
            gnu << " " <<  setw(18) << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom;
            AvEn=AvEn+ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom;
            count++;
          }
          gnu << " " << setw(18) << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).volume_atom << " ";
          gnu << " " << setw(18) << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).calculation_time << " ";
          gnu << " " <<  setw(8) << CleanNameStructure(ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).structure_name) << " ";      	
          if((ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).structure_description)=="nnn") {
            gnu << " " << setw(18) << CleanNameStructure(ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).structure_name) << endl;
          } else {
            gnu << " " <<   setw(18) << (ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).structure_description) << endl;
          }

          //gnu << "      " << (ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).structure_description) << endl;

        }
      }

      AvEn=abs(AvEn/count);

      gnu.close();
      if(_verbose) cerr << "done" << endl << "compiling structure data...";
      ofstream gnu1;
      //  [OBSOLETE] string str_info_dat_file="./str_info.dat."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+"."+aurostd::utype2string(k); //CO20200502 - threadID
      string str_info_dat_file=aurostd::TmpFileCreate("str_info_dat_"+aurostd::utype2string(k));
      gnu1.open(str_info_dat_file.c_str());
      gnu1 << setw(18) << "conc." << setw(18) << "H_" << setw(18) << "vol." << setw(18) << "S_time" << setw(18) << "HT_name" << "      " << "name" << endl;
      gnu1 << "*********************************************************************************************" << endl;
      // for(uint i=1;i<ZConcentrations.at(k).size();i++)
      for(int i=this->ZLibrary.at(k).size()-1;i>=0;i--) {
        xbool=TRUE;
        for(uint j=1;j<ZConcentrations.at(k).size();j++)
          if(i==(int) ZGNDlibrary.at(k).at(j)) xbool=FALSE;
        //	cerr << "i=" << i << "  xbool=" << xbool << endl;
        if(xbool)  {
          gnu1 << " " << setw(18) << ZLibrary.at(k).at(i).stoich_b;
          gnu1 << " " << setw(18) << ZLibrary.at(k).at(i).enthalpy_formation_atom;
          gnu1 << " " << setw(18) << ZLibrary.at(k).at(i).volume_atom;
          gnu1 << " " << setw(18) << ZLibrary.at(k).at(i).calculation_time;
          gnu1 << " " << setw(18) << CleanNameStructure(ZLibrary.at(k).at(i).structure_name);
          if((ZLibrary.at(k).at(i).structure_description)=="nnn") {
            gnu1 << "      " << CleanNameStructure(ZLibrary.at(k).at(i).structure_name) << endl;
          } else {
            gnu1 << "      " << (ZLibrary.at(k).at(i).structure_description) << endl;
          }
        }
      }
      gnu1.close();
      bool immi=false;
      count=0;
      // if(abs(AvEn)<1.0e-12) // RICHARD a double can not be ==0
      if(abs(AvEn)<0.001) // RICHARD a double can not be ==0  we might have better
      { //CO20200106 - patching for auto-indenting
        immi=true;
        for(int i=this->ZLibrary.at(k).size()-1;i>=0;i--) {
          xbool=TRUE;
          for(uint j=1;j<ZConcentrations.at(k).size();j++)
            if(i==(int) ZGNDlibrary.at(k).at(j))
              xbool=FALSE;
          AvEn+=ZLibrary.at(k).at(i).enthalpy_formation_atom;
          count++;
        }
        AvEn=AvEn/count;
      }

      immi=TRUE;
      for(uint i=0;i<ZLibrary.at(k).size()&&immi;i++) if(ZLibrary.at(k).at(i).enthalpy_formation_atom<-0.001) immi=FALSE;  

      if(_verbose) cerr << "done" << endl;
      // GNUPLOT SCRIPT EPS and PNG
      cerr << "preparing gnuplot script..." << endl;
      oss << "reset" << endl;  // richard says so
      /*
         cerr << "XHOST.vflag_control.flag("PRINT_MODE::PNG")=" << XHOST.vflag_control.flag("PRINT_MODE::PNG") << endl;
         cerr << "XHOST.vflag_control.flag("PRINT_MODE::PDF")=" << XHOST.vflag_control.flag("PRINT_MODE::PDF") << endl;
         cerr << "XHOST.vflag_control.flag("PRINT_MODE::JPG")=" << XHOST.vflag_control.flag("PRINT_MODE::JPG") << endl;
         cerr << "XHOST.vflag_control.flag("PRINT_MODE::GIF")=" << XHOST.vflag_control.flag("PRINT_MODE::GIF") << endl;
         cerr << "XHOST.vflag_control.flag("PRINT_MODE::EPS")=" << XHOST.vflag_control.flag("PRINT_MODE::EPS") << endl;
         cerr << "aflags.vflag.flag("APENNSY::WEB")=" << aflags.vflag.flag("APENNSY::WEB") << endl;
         */
      if(XHOST.vflag_control.flag("PRINT_MODE::PNG")) {
        oss << "set terminal png truecolor nocrop enhanced font '" << aflags.APENNSY_GNUPLOT_FONT_str << ",11' " << GNUPLOT_BITMAP_SIZE << " " << endl;
        oss << "set output '" << speciesAB.at(k) << ".png'" << endl;//make name automatic  // drop MATLAB
      }
      if(XHOST.vflag_control.flag("PRINT_MODE::JPG")) {
        oss << "set terminal jpeg nointerlace nocrop font '" << aflags.APENNSY_GNUPLOT_FONT_str << ",11' " << GNUPLOT_BITMAP_SIZE << " " << endl;
        oss << "set output '" << speciesAB.at(k) << ".jpg'" << endl;//make name automatic  // drop MATLAB
      }
      if(XHOST.vflag_control.flag("PRINT_MODE::GIF")) {
        oss << "set terminal gif nointerlace nocrop font '" << aflags.APENNSY_GNUPLOT_FONT_str << ",11' " << GNUPLOT_BITMAP_SIZE << " " << endl;
        oss << "set output '" << speciesAB.at(k) << ".gif'" << endl;//make name automatic  // drop MATLAB
      }
      if(XHOST.vflag_control.flag("PRINT_MODE::EPS")) {
        oss << "set terminal postscript landscape color enhanced '" << aflags.APENNSY_GNUPLOT_FONT_str << "' 11" << endl;
        oss << "set output '" << speciesAB.at(k) << ".eps'" << endl;//make name automatic  // drop MATLAB
      }
      if(XHOST.vflag_control.flag("PRINT_MODE::PDF")) {
        oss << "set terminal pdf color enhanced font '" << aflags.APENNSY_GNUPLOT_FONT_str << ", 11' " << GNUPLOT_VECTOR_SIZE << " " << endl;
        oss << "set output '" << speciesAB.at(k) << ".pdf'" << endl;//make name automatic  // drop MATLAB
      }
      string title=aurostd::string2latex(speciesAB.at(k));
      // cerr << title << endl;
      oss << "set title '[" << title << "]" << "   opar=" << alloy_order.at(k)
        << "   calcs=" << ZLibrary.at(k).size()
        << "   date: " << TODAY
        << "   [" << _APENNSY_COPYRIGHT_ << "]' font '" << aflags.APENNSY_GNUPLOT_FONT_str << ", 11'" << endl;
      oss << "set xlabel '" << speciesB.at(k) << "' font '" << aflags.APENNSY_GNUPLOT_FONT_str << ",13'" << endl;
      oss << "set ylabel 'eV/atom' font '" << aflags.APENNSY_GNUPLOT_FONT_str << ", 13'" << endl;
      oss << "#set label '" << speciesA.at(k) << "' at screen .1, screen .02" << endl;//make label automatic
      oss << "#set label '" << speciesB.at(k) << "' at screen .96, screen .02" << endl;//make label automatic
      oss << "set xrange[0:1]" << endl;
      if(1) { // RICHARD THIS DOES NOT WORK
        if(immi==true) {
          //	  cerr << "immi=1 AvEn=" << AvEn << endl;
          oss << "set yrange[" << -AvEn/1.618033 << ":" << AvEn+1 << "]" << endl;//abs(ZLibrary.at(1).at(5).enthalpy_formation_atom) << "]" << endl;
        } else {
          //	  cerr << "immi=0 AvEn=" << AvEn << endl;
          oss << "set yrange[:" << AvEn << "]" << endl;//abs(ZLibrary.at(1).at(5).enthalpy_formation_atom) << "]" << endl;
        }
      } else {
        oss << "set yrange[" << minE*1.075 << ":" << maxE << "]" << endl;   // RICHARD WHY NOT USING THIS
      }
      oss << "yoffset=-" << AvEn/10 << endl;//horizontal and vertical offsets to make plot nice
      oss << "xoffset=0" << endl;
      if(aflags.vflag.flag("APENNSY::WEB")) {
        oss << "plot '" << hull_dat_file << "'u 1:2 w linespoints pt 1 lw 1.3 ps 2 lt rgb 'blue' t '','" << str_info_dat_file << "' u 1:2 w points pt 2 lt rgb 'red' ps .8 t ''" << endl;
      } else {
        if(XHOST.vflag_control.flag("PRINT_MODE::EPS") || XHOST.vflag_control.flag("PRINT_MODE::PDF")) {	
          oss << "plot '" << hull_dat_file << "'u 1:2 w linespoints pt 1 lw 1.3 ps 2 lt rgb 'blue' t '','' u ($1+($1/2-.5)*xoffset):($2+yoffset):6 t '' with labels font '" << aflags.APENNSY_GNUPLOT_FONT_BOLD_str << ",12' lt 2 ,'" << str_info_dat_file << "' u 1:2 w points pt 2 lt rgb 'red' ps .8 t '','' u ($1+($1/2-.5)*xoffset):($2-yoffset/2):6 t '' with labels font '" << aflags.APENNSY_GNUPLOT_FONT_str << ",4' lt 2  " << endl;
        }
        if((XHOST.vflag_control.flag("PRINT_MODE::PNG") || XHOST.vflag_control.flag("PRINT_MODE::JPG") || XHOST.vflag_control.flag("PRINT_MODE::GIF"))&& aflags.vflag.flag("APENNSY::WEB")==FALSE) {	
          oss << "plot '" << hull_dat_file << "'u 1:2 w linespoints pt 1 lw 1.3 ps 2 lt rgb 'blue' t '','' u ($1+($1/2-.5)*xoffset):($2+yoffset):6 t '' with labels font '" << aflags.APENNSY_GNUPLOT_FONT_BOLD_str << ",12' lt 2 ,'" << str_info_dat_file << "' u 1:2 w points pt 2 lt rgb 'red' ps .8 t '','' u ($1+($1/2-.5)*xoffset):($2-yoffset/2):6 t '' with labels font '" << aflags.APENNSY_GNUPLOT_FONT_str << ",4' lt 2   " << endl;
        }
      }
      if(!XHOST.vflag_control.flag("KEEP::GPL")) oss << "system `" << "rm -f " << hull_dat_file << " " <<  str_info_dat_file << "`" << endl;  // DELETE hull_dat_file,str_info_dat_file      
      //echo -e "plot '-' w lines lw 5 \n 6 6 \n 7 7 \n e" | gnuplot -persist
    }
  }
  // system("ls -alt");

  if(mode==_MATLAB_OUTPUT_MODE_) oss << "exit" << endl;
  if(mode==_GNUPLOT_OUTPUT_MODE_) {;};////oss << "exit" << endl;
  if(_verbose) cerr << "**********************************************************" << endl;
  if(_verbose) cerr << "* MAKING GROUNDSTATES PLOTS                              *" << endl;
  if(_verbose) cerr << "**********************************************************" << endl;
  // FINISHED
  return oss.str();
}

// [OBSOLETE] // **********************************************************************************************
// [OBSOLETE] // APENNSY_Parameters::APENNSY_ProtoCheck
// [OBSOLETE] // **********************************************************************************************
// [OBSOLETE] string APENNSY_Parameters::APENNSY_ProtoCheck(bool _verbose,_aflags& aflags) {
// [OBSOLETE]   stringstream oss;
// [OBSOLETE]   int LEV=0;
// [OBSOLETE]   vector<string> s;
// [OBSOLETE]   vector<vector<aflowlib::_aflowlib_entry> > Mylib;
// [OBSOLETE]   if(_verbose) cerr << "******CHECKING PROTOTYPES******" << endl;
// [OBSOLETE]   //static char string_line[1024],*string_line_ptr1;                                                                                                                                                                  
// [OBSOLETE]   //int natomsA=0,natomsB=0;                                                                                                                                                                                          
// [OBSOLETE]   uint i,j,k;
// [OBSOLETE]   //bool alphabetic;                                                                                                                                                                                                  
// [OBSOLETE]   char alloy[128];
// [OBSOLETE]   Mylib=ZLibrary;
// [OBSOLETE]   string tmp;
// [OBSOLETE] 
// [OBSOLETE]   cerr << alloysRAW.at(0) << " " << alloysRAW.at(1) << endl;
// [OBSOLETE]   xstructure strtmp;
// [OBSOLETE]   for(k=0;k<alloys.size();k++) {
// [OBSOLETE]     cerr << alloys.at(k) << endl << endl;
// [OBSOLETE]     vector<int> checkbox(ZLibrary.at(k).size());
// [OBSOLETE]     for(i=0;i<ZLibrary.at(k).size();i++) checkbox.at(i)=1;  // CHECKOUT
// [OBSOLETE]     for(i=1;i<ZConcentrations.at(k).size();i++) {
// [OBSOLETE]       if(ZGNDlibrary.at(k).at(i)>=0) {
// [OBSOLETE]         double ca,uenthalpy;
// [OBSOLETE]         vector<aflowlib::_aflowlib_entry> range;
// [OBSOLETE]         strcpy(alloy,alloysRAW.at(k).c_str());
// [OBSOLETE]         ca=ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).stoich_a;
// [OBSOLETE]         uenthalpy=ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_atom;
// [OBSOLETE]         cerr << "ENERGY  " << endl;
// [OBSOLETE]         cerr << "GND:  " << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_atom << endl << CleanNameStructure(ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).structure_name) << endl;
// [OBSOLETE]         cerr << "Array:  " << endl;
// [OBSOLETE]         for(j=0;j<ZLibrary.at(k).size();j++) {
// [OBSOLETE]           if(ZLibrary.at(k).at(j).stoich_a==ca&&checkbox.at(j)) {
// [OBSOLETE]             if(((ZLibrary.at(k).at(j).enthalpy_atom-uenthalpy)>-0.01)&&((ZLibrary.at(k).at(j).enthalpy_atom-uenthalpy)<0.01)) {
// [OBSOLETE]               checkbox.at(j)=0;
// [OBSOLETE]               cerr << ZLibrary.at(k).at(j).enthalpy_atom << "    " << ZLibrary.at(k).at(j).structure_name;
// [OBSOLETE] 	      // FILE *in_file_pointer;
// [OBSOLETE]               string in_file_name;
// [OBSOLETE]               in_file_name=alloy+ZLibrary.at(k).at(j).structure_name;
// [OBSOLETE]               //cout << in_file_name << endl;                                                                                                                                                                        
// [OBSOLETE]               in_file_name+=string("CONTCAR.relax2");
// [OBSOLETE]               // cout << in_file_name << endl;                                                                                                                                                                        
// [OBSOLETE] 	      // in_file_pointer=fopen(in_file_name.c_str(),"r");
// [OBSOLETE]               strtmp=xstructure(in_file_name,IOVASP_AUTO);
// [OBSOLETE] 	      // s=CheckProtoType(strtmp,LEV,cerr);
// [OBSOLETE]               CheckProtoType(strtmp,LEV,cerr,s);
// [OBSOLETE]               if(s.at(0)=="NEW") {
// [OBSOLETE]                 cerr << endl;
// [OBSOLETE]                 if(ZLibrary.at(k).at(j).structure_description=="nnn") {
// [OBSOLETE]                   ZLibrary.at(k).at(j).structure_description="*"+CleanNameStructure(ZLibrary.at(k).at(j).structure_name);
// [OBSOLETE]                   cerr << ZLibrary.at(k).at(j).structure_description << endl;
// [OBSOLETE]                 }
// [OBSOLETE]                 else {
// [OBSOLETE]                   ZLibrary.at(k).at(j).structure_description="*"+ZLibrary.at(k).at(j).structure_description;
// [OBSOLETE]                   cerr << ZLibrary.at(k).at(j).structure_description << endl;
// [OBSOLETE]                 }
// [OBSOLETE]               }
// [OBSOLETE]               else {
// [OBSOLETE]                 ZLibrary.at(k).at(j).structure_description=s.at(1);
// [OBSOLETE]                 cerr << ZLibrary.at(k).at(j).structure_description << endl;
// [OBSOLETE] 
// [OBSOLETE]               }
// [OBSOLETE]             }
// [OBSOLETE]           }
// [OBSOLETE]         }
// [OBSOLETE]         cerr << endl;
// [OBSOLETE]       }
// [OBSOLETE]     }
// [OBSOLETE]   }
// [OBSOLETE] 
// [OBSOLETE]   stringstream aus;aus.clear();aus.str(std::string());
// [OBSOLETE]   string FileGNUPLOT_title=aurostd::TmpFileCreate("gsplot.gp");
// [OBSOLETE]   //  string FileGNUPLOT_title="gsplot.gp";
// [OBSOLETE]   ofstream FileGNUPLOT(string(FileGNUPLOT_title).c_str());
// [OBSOLETE]   aus  << APENNSY_ConvexHull(TRUE,aflags,_GNUPLOT_OUTPUT_MODE_);
// [OBSOLETE]   FileGNUPLOT << aus.str() << endl;
// [OBSOLETE]   FileGNUPLOT.flush();FileGNUPLOT.close();
// [OBSOLETE]   aus.clear();aus.str(std::string());
// [OBSOLETE]   aus << "export DISPLAY=:0.0" << endl << XHOST.command("gnuplot") << " " << " " << FileGNUPLOT_title << endl;
// [OBSOLETE]   // [OBSOLETE] aus << "rm -f str_info.dat* hull.dat*" << endl;
// [OBSOLETE]   // [OBSOLETE] aus << "rm -f " << FileGNUPLOT_title  << endl;
// [OBSOLETE]   aurostd::execute(aus);
// [OBSOLETE]   if(!XHOST.vflag_control.flag("KEEP::GPL")) aurostd::RemoveFile(FileGNUPLOT_title);
// [OBSOLETE] 
// [OBSOLETE]   return oss.str();
// [OBSOLETE] }

// **********************************************************************************************
// APENNSY_Parameters::APENNSY_PS_EnergyList
// **********************************************************************************************
string APENNSY_Parameters::APENNSY_PS_EnergyList(_aflags &aflags) {
  stringstream oss;
  double dHf,dHfd;
  // ************************************************************************
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag")) {;}; // phony
  // ************************************************************************
  oss.setf(std::ios::fixed,std::ios::floatfield);

  cerr << "**********************************************************" << endl;
  cerr << "* MAKING PHASE SEPARATING SLOPPY ENERGY LIST             *" << endl;
  cerr << "**********************************************************" << endl;
  oss << "\\documentclass" << APENNSY_LATEX_DOC_TYPE << endl;
  oss << "\\usepackage{epsfig}" << endl;
  oss << "\\begin{document}" << endl;
  oss << "\\voffset -40mm" << endl;
  oss << "\\hoffset -5mm" << endl;
  oss << LATEX_TEXT_HEIGHT << endl; // oss << "\\textheight 240mm %LATEX" << endl;
  oss << "\\textwidth 180mm" << endl;
  oss << "" << APENNSY_LATEX_FONT_SMALL << "" << endl;

  if(0) {
    oss << endl;
    oss << "\\begin{center}" << endl;
    oss << "\\begin{tabular}{||c||} \\hline" << endl;
    oss << "Systems without intermetallic compounds \\\\ \\hline" << endl;
    oss << "\\begin{tabular}{c|c|c} " << endl;
    oss << "System & Lower Energy & E$_f$ \\\\ " << endl;
    oss << "       & Structure    & (meV/atom) \\\\ \\hline" << endl;
    for(uint k=0;k<alloys.size();k++) {
      xvector<double> FindMinimun(ZConcentrations.at(k).size()-1);
      if(Alloy2MiscibilityHT.at(k)==MISCIBILITY_SYSTEM_MISCIBLE) {
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Compound forming alloy "+alloys.at(k),_VALUE_ILLEGAL_); //CO20200624
      }
      { //find mimumum distance
        for(uint j=1;j<=ZConcentrations.at(k).size()-1;j++) {
          if(j==1 || j==ZConcentrations.at(k).size()-1) { FindMinimun(j)=APENNSY_INF; }
          else FindMinimun(j)=ZLibrary.at(k).at(RankLib.at(k).at(j).at(0)).distance_gnd;
        }
        uint j=mini(FindMinimun);
        uint ii=RankLib.at(k).at(j).at(0);dHfd=ZLibrary.at(k).at(ii).distance_gnd;
        cerr << " " << alloys.at(k);
        cerr << " " << ZConcentrations.at(k).at(j);
        cerr << " str=" << aurostd::string2utype<int>(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name)) << " (" << ii << ")";
        cerr << " dHfd=" << dHfd;
        if(ZLibrary.at(k).at(ii).vNsgroup.front()!=ZLibrary.at(k).at(ii).vNsgroup.back()) oss << "XX";else oss << "  ";
        //	  if(ZLibrary.at(k).at(ii).structure_description[1]=='n') {
        cerr << " type=" << (int) aurostd::string2utype<int>(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name));
        oss << alloys.at(k) << "& ";
        oss << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH) << "& ";
        oss.precision(5);
        oss << dHfd << "\\\\ \\hline" << endl;
        //	  } else {
        // cerr << " type=%s",ZLibrary.at(k).at(ii).structure_description);
        // oss << alloys.at(k) << "& ";
        // oss << ZLibrary.at(k).at(ii).structure_description << "& ";
        // oss.precision(5);
        // oss << dHfd << "\\\\ \\hline" << endl;
        //}
        cerr << endl;
      }
    }
    oss << "\\end{tabular} \\\\ \\hline" << endl;
    oss << "\\end{tabular}" << endl;
    oss << "\\end{center}" << endl;
  }
  for(uint k=0;k<alloys.size();k++) {
    oss << "\\begin{center}{\\LARGE \\begin{verbatim}" << speciesAB.at(k) << "\\end{verbatim}}\\end{center}" << endl;
    oss << "\\begin{verbatim}" << endl;
    oss << "% ********************************************************************************************* ";
    oss << "% ********************************************************************************************* " << endl;
    oss << "% ********************************************************************************************* ";
    oss << "% ********************************************************************************************* " << endl;
    oss << "% " << alloys.at(k) << "  " << alloysmesg.str();
    for(uint j=0;j<=ZConcentrations.at(k).size();j++) {
      if(ZConcentrations.at(k).at(j)>-CEPSILON && ZConcentrations.at(k).at(j)<1.0+CEPSILON) { // [0,1]
        oss.setf(std::ios::fixed,std::ios::floatfield);
        oss.precision(5);
        oss << "3-CONCENTRATION=" << ZConcentrations.at(k).at(j) << endl;	
        uint ii=RankLib.at(k).at(j).at(0);dHfd=ZLibrary.at(k).at(ii).distance_gnd;
        oss << "  str=" << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH) << " (" << aurostd::PaddedPRE((int) ii,APENNSY_INT_DEPTH) << ")";
        oss.precision(4);
        oss << " dHfd=" << dHfd << endl;
        uint i1=0;
        for(uint i=0;i<RankLib.at(k).at(j).size();i++) {
          ii=RankLib.at(k).at(j).at(i);
          if(aurostd::isequal(ZLibrary.at(k).at(ii).stoich_b,ZConcentrations.at(k).at(j),CEPSILON)) {
            dHf=ZLibrary.at(k).at(ii).enthalpy_formation_atom-ZLibrary.at(k).at(RankLib.at(k).at(j).at(0)).enthalpy_formation_atom;
            dHfd=ZLibrary.at(k).at(ii).distance_gnd;
            if(i<=10) {
              if(ii==ZMINElibrary.at(k).at(j)) i1=ii;
              if(i1>0 && AlloyStructureIdentity.at(k).at(i1).at(ii)) oss << " *"; else oss << "  ";	
              //oss << "str=%3i",structures_number.at(ii));
              oss << "str=" << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH) << " (" << aurostd::PaddedPRE((int) ii,APENNSY_INT_DEPTH) << ")";
              oss.precision(4);
              oss << " Hf=" << ZLibrary.at(k).at(ii).enthalpy_formation_atom;
              oss.precision(4);
              oss << " dHf=" << dHf;
              oss.precision(4);
              oss << " dHfd=f" << dHfd;
              oss.precision(3);
              oss << " Vat=" << ZLibrary.at(k).at(ii).volume_atom;
              oss.precision(4);
              oss << " bnd=[" << ZLibrary.at(k).at(ii).bond_aa << "|" << ZLibrary.at(k).at(ii).bond_ab << "|" << ZLibrary.at(k).at(ii).bond_bb << "]";
              oss << " sgF=[" << aurostd::PaddedPRE(ZLibrary.at(k).at(ii).vsgroup.back(),APENNSY_STR_DEPTH) << "]";
              oss << " sg=[" << aurostd::PaddedPRE(ZLibrary.at(k).at(ii).vsgroup.front(),APENNSY_STR_DEPTH) << "|" << aurostd::PaddedPRE(ZLibrary.at(k).at(ii).vsgroup.back(),APENNSY_STR_DEPTH) << "]";
              if(ZLibrary.at(k).at(ii).vNsgroup.front()!=ZLibrary.at(k).at(ii).vNsgroup.back()) oss << "XX";else oss << "  ";
              //oss << " type=" << structure_name.at(ii);
              if(ZLibrary.at(k).at(ii).structure_description[1]=='n') {
                oss << " type=" << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH);
              } else {
                oss << " type=" << ZLibrary.at(k).at(ii).structure_description;
              }
              // if(ZLibrary.at(k).at(i).enthalpy_formation_atom<0.001) oss << "             ";
              if(ii==ZMINElibrary.at(k).at(j)) oss << " Hmin";
              if((int) ii==ZGNDlibrary.at(k).at(j))  oss << "-GND";	
              oss << endl;
            }
          }
        }
      }
    }
    oss << "% ********************************************************************************************* ";
    oss << "% ********************************************************************************************* " << endl;
    oss << "% ********************************************************************************************* ";
    oss << "% ********************************************************************************************* " << endl;
    oss << "\\end{verbatim}" << endl;
    oss << "\\newpage" << endl;
    oss << "\\ \\ \\\\" << endl;
    oss << "\\newpage" << endl;
  }
  oss << "\\end{document}" << endl;
  cerr << "**********************************************************" << endl;
  cerr << "* MAKING ENERGY LIST                                     *" << endl;
  cerr << "**********************************************************" << endl;
  return oss.str();
}

// **********************************************************************************************
// APENNSY_Parameters::MatlabGndStatesNamesConcentrations
// **********************************************************************************************
string APENNSY_Parameters::MatlabGndStatesNamesConcentrations(bool _verbose) {  // SMALL HULLS SHULL
  stringstream oss;
  double x,y,minE,maxE,deltaE;
  bool xbool;
  bool meV=FALSE;
  double Eratio=1.0;
  if(_verbose) cerr << "***************************************************************" << endl;
  if(_verbose) cerr << "MatlabGndStatesNamesConcentrations: start" << endl;
  oss << "disp('**************************************************************************');   " << endl;
  oss << "disp('*                                                                        *');   " << endl;
  oss << "disp('*    GNDstate analyzer -  Stefano Curtarolo                              *');   " << endl;
  oss << "disp('*                                                                        *');   " << endl;
  oss << "disp('**************************************************************************');   " << endl;
  for(uint k=0;k<alloys.size();k++) {
    oss << "%disp('**************************************************************************');   " << "  % " << alloys.at(k) << endl;
    oss << "%disp('*                                                                        *');   " << "  % " << alloys.at(k) << endl;
    oss << "%disp('*    GNDstate analyzer -  Stefano Curtarolo                              *');   " << "  % " << alloys.at(k) << endl;
    oss << "%disp('*                                                                        *');   " << "  % " << alloys.at(k) << endl;
    oss << "%disp('**************************************************************************');   " << "  % " << alloys.at(k) << endl;
    oss << "%disp(' Generated = " << TODAY << "');   " << "  % " << alloys.at(k) << endl;
    oss << "disp('" << alloys.at(k) << "')" << "  % " << alloys.at(k) << endl;
    oss << "                                                                                      " << "  % " << alloys.at(k) << endl;
    oss << "clear;                                                                                " << "  % " << alloys.at(k) << endl;
    oss << "FONTDIM=12;                                                                           " << "  % " << alloys.at(k) << endl;
    oss << "COLOR=[0 0 0];                                                                        " << "  % " << alloys.at(k) << endl;
    oss << "LINEWIDTH=1;                                                                          " << "  % " << alloys.at(k) << endl;
    oss << "PAPERPOSITION=[1 1 5 4];                                                              " << "  % " << alloys.at(k) << endl;
    oss << "PAPERSIZE=[8.5 11];                                                                   " << "  % " << alloys.at(k) << endl;
    oss << "POSITION=[420 540 560 420];                                                           " << "  % " << alloys.at(k) << endl;
    oss << "AXISWIDTH=1;                                                                          " << "  % " << alloys.at(k) << endl;
    oss << "AXISPOSITION=[0.136646 0.164179 0.774845 0.78678];                                    " << "  % " << alloys.at(k) << endl;
    oss << "clf;                                                                                  " << "  % " << alloys.at(k) << endl;
    oss << "% ****************************************************************************** " << "  % " << alloys.at(k) << endl;
    oss << "% ****************************************************************************** " << "  % " << alloys.at(k) << endl;
    oss << "% CONVEX_HULL_MATLAB CODE START                                          " << "  % " << alloys.at(k) << endl;
    oss << "% " << alloys.at(k) << "  % " << alloys.at(k) << endl;
    // oss << "cd /home/auro/work/AFLOW3" << "  % " << alloys.at(k) << endl;
    oss << "H_FIGURE=figure;" << endl;
    oss << "set(H_FIGURE,'PaperUnits','inches');" << endl;
    oss << "set(H_FIGURE,'PaperOrientation','portrait');" << endl;
    oss << "set(H_FIGURE,'PaperPosition',PAPERPOSITION); % left top width height" << endl;
    oss << "set(H_FIGURE,'PaperSize',PAPERSIZE);" << endl;
    oss << "set(H_FIGURE,'PaperType','usletter');" << endl;
    oss << "set(H_FIGURE,'Position',POSITION);" << endl;
    // get boundaries
    minE=0.0;
    for(uint i=1;i<=ZConcentrations.at(k).size()-1;i++) {
      if(ZGNDlibrary.at(k).at(i)>=0)
        if(minE>ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom)
          minE=ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom;
    }
    maxE=-minE/10+0.01; // looks good here (graphically)
    if(minE<-0.100) {Eratio=1.0;meV=FALSE;} else {Eratio=1000.0;meV=TRUE;} minE*=Eratio;maxE*=Eratio;
    deltaE=maxE-minE;minE-=deltaE/20;minE-=deltaE/20;minE-=deltaE/20;deltaE=maxE-minE;
    // done

    oss << "Cb=[";
    for(uint i=1;i<=ZConcentrations.at(k).size()-1;i++)
      if(ZGNDlibrary.at(k).at(i)>=0)
        oss << "  " << 100*ZConcentrations.at(k).at(i);
    oss << "];" << endl;
    oss << "EF=[";
    for(uint i=1;i<=ZConcentrations.at(k).size()-1;i++)
      if(ZGNDlibrary.at(k).at(i)>=0) {
        if(Eratio*ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom>=0) {
          oss << "  " << Eratio*ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom;
        } else {
          oss << " " << Eratio*ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom;
        }
      }
    oss << "];" << endl;
    oss << "%L=[";
    for(uint i=1;i<=ZConcentrations.at(k).size()-1;i++)
      if(ZGNDlibrary.at(k).at(i)>=0)
        oss << "    " << CleanNameStructure(ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).structure_name);
    oss << "];" << endl;
    //
    //oss << "plot(Cb ,EF ,'b+-','LineWidth',0.5*LINEWIDTH);hold on;axis([0 100 %f %f]);\n",minE,maxE);
    oss << "plot(Cb ,EF ,'LineWidth',0.5*LINEWIDTH);hold on;axis([0 100 " << minE << " " << maxE << "])" << endl;
    for(uint i=0;i<(this->ZLibrary.at(k).size());i++) {
      xbool=TRUE; for(uint j=1;j<=ZConcentrations.at(k).size()-1;j++) if((int) i==ZGNDlibrary.at(k).at(j)) xbool=FALSE;
      if(xbool)
        oss << "plot(" << 100*ZLibrary.at(k).at(i).stoich_b << "," << Eratio*ZLibrary.at(k).at(i).enthalpy_formation_atom << ",'rx','LineWidth',0.2*LINEWIDTH);hold on;" << endl;
    }

    for(uint i=1;i<=ZConcentrations.at(k).size()-1;i++)
      if(ZGNDlibrary.at(k).at(i)>=0) {
        x=ZConcentrations.at(k).at(i);
        y=Eratio*ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).enthalpy_formation_atom;
        oss << "plot(" << 100*x << "," << y << ",'b+','LineWidth',LINEWIDTH);hold on;";
        y-=deltaE/100;
        if(x>=0.2 && x<0.5) x-=0.04;
        if(x>=0.1 && x<0.2) {x-=0.03;}
        if(x>0.5) x+=0.03;//if(x>0.5) x+=0.03;
        if(x<0.1) x-=0.02;
        if(x<=0.0) {x-=0.04;y+=deltaE/15;}
        if(x>=1.0) {x+=0.01;y+=deltaE/15;}
        if(x>0.0 && x<1.0) y+=deltaE/40;
        // cerr << x << "," << y << endl;
        if(x<=0.0 || x>=1.0) y+=deltaE/15; //RHT FIX THIS
        oss << "text(" << 100*x-0.5*strlen(ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).structure_name.c_str()) << "," << y-deltaE/15 << ",'" << ZLibrary.at(k).at(ZGNDlibrary.at(k).at(i)).structure_name << "','fontsize',0.75*FONTDIM);" << endl;
      }
    oss << "text(" << -5.0 << "," << minE-deltaE/10 << ",'" << PhaseDiagramNameLeft.at(k) << "','fontsize',FONTDIM);" << endl;
    oss << "text(" << 100.0 << "," << minE-deltaE/10 << ",'" << PhaseDiagramNameRight.at(k) << "','fontsize',FONTDIM);" << endl;
    oss << "drawnow;" << "  % " << alloys.at(k) << endl;
    oss << "xlabel('Atomic Percent " << GetElementName(PhaseDiagramNameRight.at(k)) << "','fontsize',FONTDIM);" << endl;
    if(meV==FALSE) oss << "ylabel('eV/atom','fontsize',FONTDIM);" << endl;
    if(meV==TRUE) oss << "ylabel('meV/atom','fontsize',FONTDIM);" << endl;
    oss << "H_AXIS=gca;" << endl;
    oss << "set(H_AXIS,'LineWidth',AXISWIDTH);" << endl;
    oss << "set(H_AXIS,'Position',AXISPOSITION);" << endl;
    oss << "set(H_AXIS,'FontSize',FONTDIM);" << endl;
    oss << "print -depsc fig_" << speciesAB.at(k) << ".eps " << endl;
    // oss << "print -dpdf fig_" << speciesAB.at(k) << ".pdf " << endl;
    //if(XHOST.vflag_control.flag("PRINT_MODE::PDF")) 
    // oss << "system('ps2pdf " << speciesAB.at(k) << ".eps " << speciesAB.at(k) << ".pdf " << "');   % " << alloys.at(k) << endl;
    oss << "% CONVEX_HULL_MATLAB CODE END                                         " << "  % " << alloys.at(k) << endl;
    oss << "% *********************************************************************************** " << "  % " << alloys.at(k) << endl;
  }
  if(_verbose) cerr << "MatlabGndStatesNamesConcentrations: stop" << endl;
  if(_verbose) cerr << "***************************************************************" << endl;
  oss << "exit" << endl;
  // done
  return oss.str();
}

// **********************************************************************************************
// APENNSY_Parameters::APENNSY_HistogramList
// **********************************************************************************************
string APENNSY_Parameters::APENNSY_HistogramList(_aflags &aflags) {
  stringstream oss;
  // ************************************************************************
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag")) {;}; // phony
  // ************************************************************************
  cerr << "**********************************************************" << endl;
  cerr << "* MAKING HISTOGRAM LIST                                  *" << endl;
  cerr << "**********************************************************" << endl;
  // oss << max(EFDlibrary) << endl;
  xvector<int> EFDhistogram(0,1000);
  EFDhistogram.reset(); //CO20190808 reset sets all to 0
  uint i,j,k;

  for(j=0;j<alloys.size();j++) {
    for(i=0;i<ZLibrary.at(j).size();i++) {
      k=(int) ((double) ZLibrary.at(j).at(i).distance_gnd*1000+0.999);
      // oss << k << " " << ((double) EFDlibrary(i,j)*1000) << endl;
      if(ZLibrary.at(j).at(i).stoich_b>apennsy_epsilon || ZLibrary.at(j).at(i).stoich_b<1-apennsy_epsilon) {
        if(k<=1000) EFDhistogram(k)++;
      } } }
  oss << "Ehisto=[";
  for(i=0;i<=200;i++) {
    oss << EFDhistogram(i);// << endl;
    if(i<200) oss << ",";
  }
  oss << "]" << endl;
  for(i=0;i<=1000;i++)
    oss << EFDhistogram(i) << "   " << "]" << i-1 << "," << i << "] meV/atom" << endl;
  // done
  return oss.str();
}

// **********************************************************************************************
// APENNSY_Parameters::APENNSY_Rules
// **********************************************************************************************
string APENNSY_Parameters::APENNSY_Rules(_aflags &aflags) {
  stringstream oss;
  // ************************************************************************
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag")) {;}; // phony
  // ************************************************************************
  string elements_char[]={"Ag","Au","Cd","Co","Cr","Cu","Fe","Hf","Hg","Ir","La","Mn","Mo","Nb","Ni","Os","Pd","Pt","Re","Rh","Ru","Sc","Ta","Tc","Ti","V","W","Y","Zn","Zr","XXX"};

  vector<string> elements;
  for(uint i=0;elements_char[i]!="XXX";i++)
    elements.push_back(elements_char[i]);
  SortAtomsPettiforScale(elements);
  // some DEBUG
  uint i,j,k;
  cerr << "RULES" << endl;
  vector<uint> alloy_miscible_HT,alloy_miscible_Experiments,alloy_miscible_Miedema,structures_available;
  vector<vector<int> > rules;
  uint alloy_miscible_HT_size;
  uint structures_available_size;
  for(k=0;k<alloys.size();k++) {
    if(Alloy2MiscibilityHT.at(k)==MISCIBILITY_SYSTEM_MISCIBLE) alloy_miscible_HT.push_back(k);
    if(Alloy2MiscibilityExperiments.at(k)==MISCIBILITY_SYSTEM_MISCIBLE) alloy_miscible_Experiments.push_back(k);
    if(Alloy2MiscibilityMiedema.at(k)==MISCIBILITY_SYSTEM_MISCIBLE) alloy_miscible_Miedema.push_back(k);
  } // done
  // FIX RULES // cerr << " " << alloy_miscible_HT.size() << " " << alloy_miscible_Experiments.size() << " " << alloy_miscible_Miedema.size() << " " << endl;
  // FIX RULES   for(j=1;j<structures[2].size();j++) {
  // FIX RULES    bool found=TRUE;
  // FIX RULES     for(i=0;i<alloy_miscible_HT.size()&&found;i++) {
  // FIX RULES        k=alloy_miscible_HT.at(i);
  // FIX RULES       found=TRUE;
  // FIX RULES     }
  // FIX RULES     if(found && structures[2].at(j).vconcentrations.at(0)>epsilon && structures[2].at(j).stoich_b>epsilon) structures_available.push_back(j);
  // FIX RULES     // else cerr << structures[2].at(j).name << " " << structures[2].at(j).vconcentrations.at(0) << "  " << structures[2].at(j).stoich_b << endl;
  // FIX RULES   }
  // cache
  bool CACHE_WRITE=FALSE;
  alloy_miscible_HT_size=alloy_miscible_HT.size();
  structures_available_size=structures_available.size();
  if(CACHE_WRITE) {
    // alloy_miscible_HT
    cout << alloy_miscible_HT.size() << " ";
    // structures_available
    cout << structures_available.size() << " ";
    for(j=0;j<structures_available.size();j++) cout << structures_available.at(j) << " ";
    cout << endl;
  }

  // FIX RULES  if(LDEBUG) {
  // FIX RULES    for(j=0;j<structures_available_size;j++) {
  // FIX RULES      cerr << structures[2].at(structures_available.at(j)).name << endl;
  // FIX RULES    }
  // FIX RULES  }
  cerr << " " << alloy_miscible_HT_size << " " << structures_available_size << " " << endl;
  // now make matrix of rules search
  for(i=0;i<structures_available_size;i++)
    rules.push_back(*(new vector<int>(alloy_miscible_HT_size)));
  // fill it up
  for(uint ii=0;ii<rules.size();ii++) { // rows = structures
    for(uint jj=0;jj<rules.at(ii).size();jj++) { // colums = alloys
      rules.at(ii).at(jj)=0;
      k=alloy_miscible_HT.at(jj);
      j=structures_available.at(ii);
      if(ZLibrary.at(k).at(j).enthalpy_formation_atom<0.0) rules.at(ii).at(jj)=1;
    }
  }
  // structures_available_size=rules.size()
  uint cnt,cntmax;
  // --------------------------------------------------------------------------------------------------------------------
  // try 1 body
  if(aurostd::args2flag(XHOST.argv,"-1")) {
    cerr << "RULES 1 body" << endl;
    cntmax=0;
    for(uint ii=0;ii<structures_available_size;ii++) { // rows = structures
      cnt=0;
      for(uint jj=0;jj<alloy_miscible_HT_size;jj++) // colums = alloys
        if(rules.at(ii).at(jj)==1) {cnt++;}
      if(cnt>cntmax) cntmax=cnt;
      // FIX RULES     if(cnt==alloy_miscible_HT_size) cerr << structures[2].at(structures_available.at(ii)).name << " " << cnt << endl;
    }
    cerr << "RULES 1 body max = " << cntmax << endl;
  }
  // --------------------------------------------------------------------------------------------------------------------
  // try 2 bodies
  if(aurostd::args2flag(XHOST.argv,"-2")) {
    cerr << "RULES 2 bodies" << endl;
    cntmax=0;
    for(uint ii1=0;ii1<structures_available_size;ii1++) { // rows = structures
      for(uint ii2=ii1+1;ii2<structures_available_size;ii2++) { // rows = structures
        cnt=0;
        for(uint jj=0;jj<alloy_miscible_HT_size;jj++) // colums = alloys
          if(rules[ii1].at(jj)==1 || rules[ii2].at(jj)==1) {cnt++;}
        if(cnt>cntmax) cntmax=cnt;
        // FIX RULES	if(cnt==alloy_miscible_HT_size)	  cerr << structures[2].at(structures_available.at(ii1)).name << " "	      << structures[2].at(structures_available.at(ii2)).name << " "	      << cnt << endl;
      }
    }
    cerr << "RULES 2 body max = " << cntmax << endl;
  }
  // --------------------------------------------------------------------------------------------------------------------
  // try 3 bodies
  if(aurostd::args2flag(XHOST.argv,"-3")) {
    cerr << "RULES 3 bodies" << endl;
    cntmax=0;
    for(uint ii1=0;ii1<structures_available_size;ii1++) { // rows = structures
      for(uint ii2=ii1+1;ii2<structures_available_size;ii2++) { // rows = structures
        for(uint ii3=ii2+1;ii3<structures_available_size;ii3++) { // rows = structures
          cnt=0;
          for(uint jj=0;jj<alloy_miscible_HT_size;jj++) // colums = alloys
            if(rules[ii1].at(jj)==1 || rules[ii2].at(jj)==1 || rules[ii3].at(jj)==1) {cnt++;}
          if(cnt>cntmax) {
            cntmax=cnt;
            // FIX RULES	  cerr << structures[2].at(structures_available.at(ii1)).name << " " << structures[2].at(structures_available.at(ii2)).name << " " << structures[2].at(structures_available.at(ii3)).name << " "  << cntmax << endl;
          }
        }
      }
    }
    cerr << "RULES 3 body max = " << cntmax << endl;
  }
  // --------------------------------------------------------------------------------------------------------------------
  // try 4 bodies
  if(aurostd::args2flag(XHOST.argv,"-4")) {
    cerr << "RULES 4 bodies" << endl;
    cntmax=0;
    for(uint ii1=0;ii1<structures_available_size;ii1++) { // rows = structures
      cerr << ii1 << endl;
      for(uint ii2=ii1+1;ii2<structures_available_size;ii2++) { // rows = structures
        for(uint ii3=ii2+1;ii3<structures_available_size;ii3++) { // rows = structures
          for(uint ii4=ii3+1;ii4<structures_available_size;ii4++) { // rows = structures
            cnt=0;
            for(uint jj=0;jj<alloy_miscible_HT_size;jj++) // colums = alloys
              if(rules[ii1].at(jj)==1 || rules[ii2].at(jj)==1 || rules[ii3].at(jj)==1 || rules[ii4].at(jj)==1) {cnt++;}
            if(cnt>cntmax) {
              cntmax=cnt;
              // FIX RULES	      cerr << structures[2].at(structures_available.at(ii1)).name << " "
              // FIX RULES		  << structures[2].at(structures_available.at(ii2)).name << " "
              // FIX RULES		  << structures[2].at(structures_available.at(ii3)).name << " "
              // FIX RULES		  << structures[2].at(structures_available.at(ii4)).name << " "
              // FIX RULES		  << cnt << endl;
            }
          }
        }
      }
    }
    cerr << "RULES 4 body max = " << cntmax << endl;
  }
  // --------------------------------------------------------------------------------------------------------------------
  // try 5 bodies
  if(aurostd::args2flag(XHOST.argv,"-5")) {
    cerr << "RULES 5 bodies" << endl;
    cntmax=0;
    for(uint ii1=0;ii1<structures_available_size;ii1++) { // rows = structures
      cerr << ii1 << endl;
      for(uint ii2=ii1+1;ii2<structures_available_size;ii2++) { // rows = structures
        for(uint ii3=ii2+1;ii3<structures_available_size;ii3++) { // rows = structures
          for(uint ii4=ii3+1;ii4<structures_available_size;ii4++) { // rows = structures
            for(uint ii5=ii4+1;ii5<structures_available_size;ii5++) { // rows = structures
              cnt=0;
              for(uint jj=0;jj<alloy_miscible_HT_size;jj++) // colums = alloys
                if(rules[ii1].at(jj)==1 || rules[ii2].at(jj)==1 || rules[ii3].at(jj)==1 || rules[ii4].at(jj)==1 ||
                    rules[ii5].at(jj)==1)
                {cnt++;}
              if(cnt>cntmax) {
                cntmax=cnt;
                // FIX RULES		cerr << structures[2].at(structures_available.at(ii1)).name << " "
                // FIX RULES		    << structures[2].at(structures_available.at(ii2)).name << " "
                // FIX RULES		    << structures[2].at(structures_available.at(ii3)).name << " "
                // FIX RULES		    << structures[2].at(structures_available.at(ii4)).name << " "
                // FIX RULES		    << structures[2].at(structures_available.at(ii5)).name << " "
                // FIX RULES		    << cntmax << endl;
              }
            }
          }
        }
      }
    }
    cerr << "RULES 5 body max = " << cntmax << endl;
  }
  // --------------------------------------------------------------------------------------------------------------------
  // try 6 bodies
  if(aurostd::args2flag(XHOST.argv,"-6")) {
    cerr << "RULES 6 bodies" << endl;
    cntmax=0;
    for(uint ii1=0;ii1<structures_available_size;ii1++) { // rows = structures
      cerr << ii1 << endl;
      for(uint ii2=ii1+1;ii2<structures_available_size;ii2++) { // rows = structures
        for(uint ii3=ii2+1;ii3<structures_available_size;ii3++) { // rows = structures
          for(uint ii4=ii3+1;ii4<structures_available_size;ii4++) { // rows = structures
            for(uint ii5=ii4+1;ii5<structures_available_size;ii5++) { // rows = structures
              for(uint ii6=ii5+1;ii6<structures_available_size;ii6++) { // rows = structures
                cnt=0;
                for(uint jj=0;jj<alloy_miscible_HT_size;jj++) // colums = alloys
                  if(rules[ii1].at(jj)==1 || rules[ii2].at(jj)==1 || rules[ii3].at(jj)==1 || rules[ii4].at(jj)==1 ||
                      rules[ii5].at(jj)==1 || rules[ii6].at(jj)==1)
                  {cnt++;} else {jj=alloy_miscible_HT_size;};
                if(cnt>cntmax) {
                  cntmax=cnt;
                  // FIX RULES		  cerr << structures[2].at(structures_available.at(ii1)).name << " "
                  // FIX RULES		      << structures[2].at(structures_available.at(ii2)).name << " "
                  // FIX RULES		      << structures[2].at(structures_available.at(ii3)).name << " "
                  // FIX RULES		      << structures[2].at(structures_available.at(ii4)).name << " "
                  // FIX RULES		      << structures[2].at(structures_available.at(ii5)).name << " "
                  // FIX RULES		      << structures[2].at(structures_available.at(ii6)).name << " "
                  // FIX RULES		      << cntmax << endl;
                }
              }
            }
          }
        }
      }
    }
    cerr << "RULES 6 body max = " << cntmax << endl;
  }
  // --------------------------------------------------------------------------------------------------------------------
  // try 7 bodies
  if(aurostd::args2flag(XHOST.argv,"-7")) {
    cerr << "RULES 7 bodies" << endl;
    cntmax=0;
    for(uint ii1=0;ii1<structures_available_size;ii1++) { // rows = structures
      cerr << ii1 << endl;
      for(uint ii2=ii1+1;ii2<structures_available_size;ii2++) { // rows = structures
        for(uint ii3=ii2+1;ii3<structures_available_size;ii3++) { // rows = structures
          for(uint ii4=ii3+1;ii4<structures_available_size;ii4++) { // rows = structures
            for(uint ii5=ii4+1;ii5<structures_available_size;ii5++) { // rows = structures
              for(uint ii6=ii5+1;ii6<structures_available_size;ii6++) { // rows = structures
                for(uint ii7=ii6+1;ii7<structures_available_size;ii7++) { // rows = structures
                  cnt=0;
                  for(uint jj=0;jj<alloy_miscible_HT_size;jj++) // colums = alloys
                    if(rules[ii1].at(jj)==1 || rules[ii2].at(jj)==1 || rules[ii3].at(jj)==1 || rules[ii4].at(jj)==1 ||
                        rules[ii5].at(jj)==1 || rules[ii6].at(jj)==1 || rules[ii7].at(jj)==1)
                    {cnt++;} else {jj=alloy_miscible_HT_size;};
                  if(cnt>cntmax) {
                    cntmax=cnt;
                    // FIX RULES		    cerr << structures[2].at(structures_available.at(ii1)).name << " "
                    // FIX RULES			 << structures[2].at(structures_available.at(ii2)).name << " "
                    // FIX RULES			 << structures[2].at(structures_available.at(ii3)).name << " "
                    // FIX RULES			 << structures[2].at(structures_available.at(ii4)).name << " "
                    // FIX RULES			 << structures[2].at(structures_available.at(ii5)).name << " "
                    // FIX RULES			 << structures[2].at(structures_available.at(ii6)).name << " "
                    // FIX RULES			 << structures[2].at(structures_available.at(ii7)).name << " "
                    // FIX RULES			 << cntmax << endl;
                  }
                }
              }
            }
          }
        }
      }
    }
    cerr << "RULES 7 body max = " << cntmax << endl;
  }
  // --------------------------------------------------------------------------------------------------------------------
  // try 8 bodies
  if(aurostd::args2flag(XHOST.argv,"-8")) {
    cerr << "RULES 8 bodies" << endl;
    cntmax=0;
    for(uint ii1=0;ii1<structures_available_size;ii1++) { // rows = structures
      cerr << ii1 << endl;
      for(uint ii2=ii1+1;ii2<structures_available_size;ii2++) { // rows = structures
        for(uint ii3=ii2+1;ii3<structures_available_size;ii3++) { // rows = structures
          for(uint ii4=ii3+1;ii4<structures_available_size;ii4++) { // rows = structures
            for(uint ii5=ii4+1;ii5<structures_available_size;ii5++) { // rows = structures
              for(uint ii6=ii5+1;ii6<structures_available_size;ii6++) { // rows = structures
                for(uint ii7=ii6+1;ii7<structures_available_size;ii7++) { // rows = structures
                  for(uint ii8=ii7+1;ii8<structures_available_size;ii8++) { // rows = structures
                    cnt=0;
                    for(uint jj=0;jj<alloy_miscible_HT_size;jj++) // colums = alloys
                      if(rules[ii1].at(jj)==1 || rules[ii2].at(jj)==1 || rules[ii3].at(jj)==1 || rules[ii4].at(jj)==1 ||
                          rules[ii5].at(jj)==1 || rules[ii6].at(jj)==1 || rules[ii7].at(jj)==1 || rules[ii8].at(jj)==1)
                      {cnt++;} else {jj=alloy_miscible_HT_size;};
                    if(cnt>cntmax) {
                      cntmax=cnt;
                      // FIX RULES		      cerr << structures[2].at(structures_available.at(ii1)).name << " "
                      // FIX RULES			  << structures[2].at(structures_available.at(ii2)).name << " "
                      // FIX RULES			  << structures[2].at(structures_available.at(ii3)).name << " "
                      // FIX RULES			  << structures[2].at(structures_available.at(ii4)).name << " "
                      // FIX RULES			  << structures[2].at(structures_available.at(ii5)).name << " "
                      // FIX RULES			  << structures[2].at(structures_available.at(ii6)).name << " "
                      // FIX RULES			  << structures[2].at(structures_available.at(ii7)).name << " "
                      // FIX RULES			  << structures[2].at(structures_available.at(ii8)).name << " "
                      // FIX RULES			  << cntmax << endl;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    cerr << "RULES 8 body max = " << cntmax << endl;
  }
  return oss.str();
}

// **********************************************************************************************
// APENNSY_Parameters::APENNSY_Structures
// **********************************************************************************************
string APENNSY_Parameters::APENNSY_Structures(_aflags &aflags) {
  stringstream oss;
  int i,j,k,ii;
  // ************************************************************************
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag")) {;}; // phony
  // ************************************************************************
  k=1;
  for(i=1;i<=(int)(ZConcentrations.at(k).size()-1)/2;i++) {
    for(j=0;j<(int)(this->ZLibrary.at(k).size());j++) {
      ii=j;
      if((aurostd::string2utype<int>(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name))>100 &&
            aurostd::string2utype<int>(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name))<180) ||
          aurostd::string2utype<int>(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name)) > 178)
        if(ZLibrary.at(k).at(ii).stoich_b==ZConcentrations.at(k).at(i) || ZLibrary.at(k).at(ii).stoich_b==ZConcentrations.at(k).at(ZConcentrations.at(k).size()-1+1-i)) {
          oss << " str=" << aurostd::PaddedPRE(CleanNameStructure(ZLibrary.at(k).at(ii).structure_name),APENNSY_STR_DEPTH) << " (" << aurostd::PaddedPRE(ii,APENNSY_INT_DEPTH) << ")";
          oss.precision(4);
          oss << " (check) Ca=" << ZLibrary.at(k).at(ii).stoich_b;
          oss.precision(4);
          oss << " (check) Cb= " << 1.0-ZLibrary.at(k).at(ii).stoich_b;
          oss.precision(2);
          oss << " At=" << ZLibrary.at(k).at(ii).natoms;
          oss.precision(2);
          oss << " Na=" << (int) ((double) 0.01+ZLibrary.at(k).at(ii).natoms*ZLibrary.at(k).at(ii).stoich_b);
          oss.precision(2);
          oss << " Nb=" << (int) ((double) 0.01+ZLibrary.at(k).at(ii).natoms*(1.0-ZLibrary.at(k).at(ii).stoich_b));
          oss << " sg=[" << aurostd::PaddedPRE(ZLibrary.at(k).at(ii).vsgroup.front(),APENNSY_STR_DEPTH) << "]";
          oss << " type=" << ZLibrary.at(k).at(ii).structure_description;
          if(0) {if(ZLibrary.at(k).at(ii).vNsgroup.front()!=ZLibrary.at(k).at(ii).vNsgroup.back()) oss << "XX";else oss << "  ";}
          oss << endl;
        }
    }
  }
  // done
  return oss.str();
}

// **********************************************************************************************
// APENNSY_Parameters::APENNSY_Order
// **********************************************************************************************
string APENNSY_Parameters::APENNSY_Order(_aflags &aflags) {
  stringstream oss;
  // ************************************************************************
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag")) {;}; // phony
  // ************************************************************************
  string elements_char[]={"Ag","Au","Cd","Co","Cr","Cu","Fe","Hf","Hg","Ir","La","Mn","Mo","Nb","Ni","Os","Pd","Pt","Re","Rh","Ru","Sc","Ta","Tc","Ti","V","W","Y","Zn","Zr","XXX"};
  vector<string> elements;
  for(uint i=0;elements_char[i]!="XXX";i++)
    elements.push_back(elements_char[i]);
  SortAtomsPettiforScale(elements);
  // some DEBUG
  if(0) {
    vector<double> vvalue;
    vector<int> vorder;
    SortAtomsPettiforScale(elements,vorder,vvalue);
    for(uint i=0;i<elements.size();i++)
      cerr << i << " " << elements.at(i) << " " << GetAtomPettiforScale(elements.at(i)) << endl;
  }

  string systemAB,systemBA;
  int k;

  // first line
  uint space=5;

  bool COMPOUND=FALSE,NOMIX=FALSE,UNKNOWN=FALSE,NOTSTUDIED=FALSE;
  int mixAB,mixBA;

  oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  for(uint j=0;j<elements.size();j++) {
    if(j==0) oss << "elements={";
    oss << "'" << elements.at(j) << "'";
    if(j<elements.size()-1) oss << ",";
    if(j==elements.size()-1) oss << "};" << endl;
  }

  oss << "%  ";// oss << aurostd::PaddedPRE(" ",space);
  for(uint i=0;i<elements.size();i++)
    oss << aurostd::PaddedPRE(elements.at(i),space);
  oss << endl;

  // %%%%%%%%% T DISORDER %%%%%%%%%%%%%%
  oss << "%%%%%%%%% T DISORDER %%%%%%%%%%%%%%" << endl;
  double order;
  for(uint i=0;i<elements.size();i++) {
    if(i==0) oss << "T=[" << endl;
    // oss << aurostd::PaddedPRE(elements.at(i),space);
    oss << "  ";
    for(uint j=0;j<elements.size();j++) {
      if(i==j) { order=0;
      } else {
        systemAB=elements.at(i)+elements.at(j);
        systemBA=elements.at(j)+elements.at(i);
        k=-1;
        for(uint kk=0;kk<alloys.size()&&k==-1;kk++) {
          if(speciesAB[kk]==systemAB || speciesAB[kk]==systemBA) k=kk;
        }
        // check if k==-1
        order=alloy_order.at((uint) k);
      }
      oss << aurostd::PaddedPRE(aurostd::utype2string(int(order)),space);
    }
    oss << endl;
    if(i==elements.size()-1) oss << "];" << endl;
  }
  // %%%%%%%%% X EXPERIMENTS %%%%%%%%%%%%%%
  oss << "%%%%%%%%% X EXPERIMENTS %%%%%%%%%%%%%%" << endl;
  for(uint i=0;i<elements.size();i++) {
    mixAB=0,mixBA=0;
    if(i==0) oss << "X=[" << endl;
    // oss << aurostd::PaddedPRE(elements.at(i),space);
    oss << "  ";
    for(uint j=0;j<elements.size();j++) {
      COMPOUND=FALSE,NOMIX=FALSE,UNKNOWN=FALSE,NOTSTUDIED=FALSE;
      if(i==j) { NOMIX=TRUE;
      } else {
        systemAB=elements.at(i)+elements.at(j);
        systemBA=elements.at(j)+elements.at(i);
        k=-1;
        for(uint kk=0;kk<alloys.size()&&k==-1;kk++) {
          if(speciesAB[kk]==systemAB || speciesAB[kk]==systemBA) k=kk;
        }
        // check if k==-1
        mixAB=MiscibilityExperimentsCheck(systemAB);
        mixBA=MiscibilityExperimentsCheck(systemBA);
        COMPOUND=(mixAB==MISCIBILITY_SYSTEM_MISCIBLE || mixBA==MISCIBILITY_SYSTEM_MISCIBLE);
        NOMIX=!COMPOUND && (mixAB==MISCIBILITY_SYSTEM_NOMIX || mixBA==MISCIBILITY_SYSTEM_NOMIX);
        UNKNOWN=!COMPOUND && !NOMIX;
      }
      if(COMPOUND) oss << aurostd::PaddedPRE("1",space);
      if(NOMIX) oss << aurostd::PaddedPRE("0",space);
      if(UNKNOWN || NOTSTUDIED) oss << aurostd::PaddedPRE("-1",space);
    }
    oss << endl;
    if(i==elements.size()-1) oss << "];" << endl;
  }

  // %%%%%%%%% M MIEDEMA %%%%%%%%%%%%%%
  oss << "%%%%%%%%% M MIEDEMA %%%%%%%%%%%%%%" << endl;
  mixAB=0,mixBA=0;
  for(uint i=0;i<elements.size();i++) {
    if(i==0) oss << "M=[" << endl;
    // oss << aurostd::PaddedPRE(elements.at(i),space);
    oss << "  ";
    for(uint j=0;j<elements.size();j++) {
      COMPOUND=FALSE,NOMIX=FALSE,UNKNOWN=FALSE,NOTSTUDIED=FALSE;
      if(i==j) { NOMIX=TRUE;
      } else {
        systemAB=elements.at(i)+elements.at(j);
        systemBA=elements.at(j)+elements.at(i);
        k=-1;
        for(uint kk=0;kk<alloys.size()&&k==-1;kk++) {
          if(speciesAB[kk]==systemAB || speciesAB[kk]==systemBA) k=kk;
        }
        // check if k==-1
        mixAB=MiscibilityMiedemaCheck(systemAB);
        mixBA=MiscibilityMiedemaCheck(systemBA);
        COMPOUND=(mixAB==MISCIBILITY_SYSTEM_MISCIBLE || mixBA==MISCIBILITY_SYSTEM_MISCIBLE);
        NOMIX=!COMPOUND && (mixAB==MISCIBILITY_SYSTEM_NOMIX || mixBA==MISCIBILITY_SYSTEM_NOMIX);
        UNKNOWN=!COMPOUND && !NOMIX;
      }
      if(COMPOUND) oss << aurostd::PaddedPRE("1",space);
      if(NOMIX) oss << aurostd::PaddedPRE("0",space);
      if(UNKNOWN || NOTSTUDIED) oss << aurostd::PaddedPRE("-1",space);
    }
    oss << endl;
    if(i==elements.size()-1) oss << "];" << endl;
  }

  // %%%%%%%%% M HUME-ROTHERY %%%%%%%%%%%%%%
  oss << "%%%%%%%%% H HUME-ROTHERY %%%%%%%%%%%%%%" << endl;
  mixAB=0,mixBA=0;
  for(uint i=0;i<elements.size();i++) {
    if(i==0) oss << "H=[" << endl;
    // oss << aurostd::PaddedPRE(elements.at(i),space);
    oss << "  ";
    for(uint j=0;j<elements.size();j++) {
      COMPOUND=FALSE,NOMIX=FALSE,UNKNOWN=FALSE,NOTSTUDIED=FALSE;
      if(i==j) { NOMIX=TRUE;
      } else {
        systemAB=elements.at(i)+elements.at(j);
        systemBA=elements.at(j)+elements.at(i);
        k=-1;
        for(uint kk=0;kk<alloys.size()&&k==-1;kk++) {
          if(speciesAB[kk]==systemAB || speciesAB[kk]==systemBA) k=kk;
        }
        // check if k==-1
        mixAB=MiscibilityHumeRotheryCheck(systemAB);
        mixBA=MiscibilityHumeRotheryCheck(systemBA);
        COMPOUND=(mixAB==MISCIBILITY_SYSTEM_MISCIBLE || mixBA==MISCIBILITY_SYSTEM_MISCIBLE);
        NOMIX=!COMPOUND && (mixAB==MISCIBILITY_SYSTEM_NOMIX || mixBA==MISCIBILITY_SYSTEM_NOMIX);
        UNKNOWN=!COMPOUND && !NOMIX;
      }
      if(COMPOUND) oss << aurostd::PaddedPRE("1",space);
      if(NOMIX) oss << aurostd::PaddedPRE("0",space);
      if(UNKNOWN || NOTSTUDIED) oss << aurostd::PaddedPRE("-1",space);
    }
    oss << endl;
    if(i==elements.size()-1) oss << "];" << endl;
  }
  return oss.str();
}

// **********************************************************************************************
// APENNSY_Parameters::APENNSY_SimulationInformation
// **********************************************************************************************
string APENNSY_Parameters::APENNSY_SimulationInformation(_aflags &aflags,string mode) {
  stringstream oss;
  // ************************************************************************
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag")) {;}; // phony
  // ************************************************************************
  if(mode=="TIME,CORES,TIMECORES,MEM") oss << aurostd::PaddedPRE("time(sec)",10) << " "
    << aurostd::PaddedPRE("cores",3) << " " 
      << aurostd::PaddedPRE("time*cores",10) << ""
      << aurostd::PaddedPRE("mem(MB)",8) << " " << endl;

  for(uint k=0;k<alloys.size();k++) {  // oss << " // " << alloys.at(k) << endl;
    for(uint j=0;j<ZLibrary.at(k).size();j++) {
      // string systemstructure=alloys.at(k)+"/"+ZLibrary.at(k1).at(j).structure_name;
      if(mode=="TIME") oss << ZLibrary.at(k).at(j).calculation_time << endl;
      if(mode=="CORES") oss << ZLibrary.at(k).at(j).calculation_cores << endl;
      if(mode=="TIMECORES") oss << ZLibrary.at(k).at(j).calculation_time*double(ZLibrary.at(k).at(j).calculation_cores) << endl;
      if(mode=="MEM") oss << ZLibrary.at(k).at(j).calculation_memory << endl;

      if(mode=="TIME,CORES,TIMECORES,MEM") oss << aurostd::PaddedPRE(ZLibrary.at(k).at(j).calculation_time,10) << " "
        << aurostd::PaddedPRE(ZLibrary.at(k).at(j).calculation_cores,3) << " " 
          << aurostd::PaddedPRE(ZLibrary.at(k).at(j).calculation_time*double(ZLibrary.at(k).at(j).calculation_cores),10) << " "
          << aurostd::PaddedPRE(ZLibrary.at(k).at(j).calculation_memory,8) << " " << endl;

      // 	oss << aurostd::PaddedPRE(aurostd::utype2string(ZLibrary.at(k).at(j).calculation_time),10) << " % " << aurostd::PaddedPOST(systemstructure,13) << "  proto = " << structures_description.at(j) << endl;
      // else
      //	oss << k << " " << j << " UN-AVAILABLE" << endl;	
    }
  }
  // done
  return oss.str();
}

// **********************************************************************************************
// APENNSY_Parameters::APENNSY_Miscibility
// **********************************************************************************************
string APENNSY_Parameters::APENNSY_Miscibility(_aflags &aflags) {
  stringstream oss;
  // ************************************************************************
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag")) {;}; // phony
  // ************************************************************************
  oss << "// ***************************************************************************" << endl;
  oss << "// *                                                                         *" << endl;
  oss << "// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *" << endl;
  oss << "// *                                                                         *" << endl;
  oss << "// ***************************************************************************" << endl;
  oss << "// Stefano Curtarolo - 2009 Duke" << endl;
  oss << "" << endl;
  oss << "// ***************************************************************************" << endl;
  oss << "// Total structures holes      = " << LIBstructures_holes << endl;
  oss << "// Total structures calculated = " << LIBstructures_calcs << endl;
  oss << "// Extracted in date  [date=" << TODAY << "]" << endl;
  oss << "// ***************************************************************************" << endl;
  oss << "" << endl;
  oss << "#ifndef _AFLOW_NOMIX_CPP" << endl;
  oss << "#define _AFLOW_NOMIX_CPP" << endl;
  oss << "#include \"aflow.h\"" << endl;
  oss << "" << endl;
  oss << "// ***************************************************************************" << endl;
  oss << "int MiscibilityCheck(string system_in) {  // (nomix,unknown,mix)" << endl;
  oss << "string system=system_in;vector<string> vspecies;" << endl;
  oss << " // XATOM_SplitAlloySpecies(KBIN::VASP_PseudoPotential_CleanName(system_in),vspecies);" << endl;
  oss << "XATOM_AlphabetizationSpecies(system,vspecies);" << endl;
  oss << "if(system.find(\"Cs\")!=string::npos) return MISCIBILITY_SYSTEM_NOMIX; // Cs bug !" << endl;
  oss << "if(vspecies.size()==1) return MISCIBILITY_SYSTEM_MISCIBLE; // pure element is miscible with itself" << endl;
  // saved in alphabetic order
  // oss << "  string system=KBIN::VASP_PseudoPotential_CleanName(system_in);" << endl;
  // oss << "  if(system.length()<2) return MISCIBILITY_SYSTEM_MISCIBLE; // pure element is miscible with itself" << endl;
  // oss << "  // based on Total structures calculated=" << LIBstructures_calcs << endl;
  // oss << "  // extracted in date  [date=" << TODAY << "]";
  for(uint k=0;k<alloys.size();k++) {
    string specieA,specieB;
    // oss << "# " << speciesAB.at(k) << " (" << alloys.at(k) << ")" ;
    // oss << " [calcs=" << ZLibrary.at(k).size() << "] [date=" << TODAY << "]";
    // oss << "  if(system==\"" << speciesAB.at(k) << "\")";
    KBIN::VASP_SplitAlloySpecies(speciesAB.at(k),specieA,specieB);
    // put them always in alphabetic order
    if(specieA <= specieB) oss << "  if(system==\"" << specieA << specieB << "\")";
    if(specieA >  specieB) oss << "  if(system==\"" << specieB << specieA << "\")";
    if(Alloy2MiscibilityHT.at(k)==MISCIBILITY_SYSTEM_MISCIBLE) {
      oss << " return MISCIBILITY_SYSTEM_MISCIBLE;";
    }
    if(Alloy2MiscibilityHT.at(k)==MISCIBILITY_SYSTEM_NOMIX) {
      oss << " return MISCIBILITY_SYSTEM_NOMIX;";
    }
    if(Alloy2MiscibilityHT.at(k)==MISCIBILITY_SYSTEM_UNKNOWN) {
      oss << " return MISCIBILITY_SYSTEM_UNKNOWN;";
    }
    oss << " // " << alloys.at(k) << " calcs=" << ZLibrary.at(k).size() << "/" << MISCIBILITY_SYSTEM_CUTOFF;
    oss << " elements=" << GetAtomNumber(specieA) << "," << GetAtomNumber(specieB);
    oss << " [" << TODAY << "]";
    oss << endl;
  }
  oss << "  // If not found then not known " << endl;
  oss << "  return MISCIBILITY_SYSTEM_UNKNOWN; // unknown " << endl;
  oss << "} " << endl;
  oss << "// ***************************************************************************" << endl;
  oss << "// Total structures holes      = " << LIBstructures_holes << endl;
  oss << "// Total structures calculated = " << LIBstructures_calcs << endl;
  oss << "// Extracted in date  [date=" << TODAY << "]" << endl;
  oss << "// ***************************************************************************" << endl;
  oss << "#endif " << endl;
  oss << " " << endl;
  oss << "// ***************************************************************************" << endl;
  oss << "// *                                                                         *" << endl;
  oss << "// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *" << endl;
  oss << "// *                                                                         *" << endl;
  oss << "// ***************************************************************************" << endl;
  // done
  return oss.str();
}

// **********************************************************************************************
// APENNSY_Parameters::APENNSY_Miscibility_Experiments
// **********************************************************************************************
string APENNSY_Parameters::APENNSY_Miscibility_Experiments(_aflags &aflags) {
  stringstream oss;
  // ************************************************************************
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag")) {;}; // phony
  // ************************************************************************
  oss << "// ***************************************************************************" << endl;
  oss << "// *                                                                         *" << endl;
  oss << "// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *" << endl;
  oss << "// *                                                                         *" << endl;
  oss << "// ***************************************************************************" << endl;
  oss << "// Ohad Levy - 2009 Duke" << endl;
  oss << "" << endl;
  oss << "// based on Massalski and the Pauling File " << endl;
  oss << "" << endl;
  oss << "#ifndef _AFLOW_MIX_EXPERIMENTS_CPP" << endl;
  oss << "#define _AFLOW_MIX_EXPERIMENTS_CPP" << endl;
  oss << "#include \"aflow.h\"" << endl;
  oss << "" << endl;
  oss << "// ***************************************************************************" << endl;
  oss << "int MiscibilityExperimentsCheck(string system_in) {  // (nomix,unknown,mix)" << endl;
  oss << "  string system=KBIN::VASP_PseudoPotential_CleanName(system_in);" << endl;
  oss << "  if(system.length()<2) return MISCIBILITY_SYSTEM_MISCIBLE; // pure element is miscible with itself" << endl;
  oss << "  // int C=MISCIBILITY_SYSTEM_MISCIBLE; " << endl;
  oss << "  int S=MISCIBILITY_SYSTEM_NOMIX; " << endl;
  oss << "  int U=MISCIBILITY_SYSTEM_UNKNOWN; " << endl;
  oss << "  int N=MISCIBILITY_SYSTEM_NOT_STUDIED; " << endl;
  // oss << "  // based on Total structures calculated=" << LIBstructures_calcs << endl;
  // oss << "  // extracted in date  [date=" << TODAY << "]";
  for(uint k=0;k<alloys.size();k++) {
    // oss << "# " << speciesAB.at(k) << " (" << alloys.at(k) << ")" ;
    // oss << " [calcs=" << ZLibrary.at(k).size() << "] [date=" << TODAY << "]";
    oss << "  if(system==\"" << speciesAB.at(k) << "\")";
    if(Alloy2MiscibilityExperiments.at(k)==MISCIBILITY_SYSTEM_MISCIBLE) {
      oss << " return return MISCIBILITY_SYSTEM_MISCIBLE;";
    }
    if(Alloy2MiscibilityExperiments.at(k)==MISCIBILITY_SYSTEM_NOMIX) {
      oss << " return S;";
    }
    if(Alloy2MiscibilityExperiments.at(k)==MISCIBILITY_SYSTEM_UNKNOWN) {
      oss << " return U;";
    }
    oss << " // Massalski and Pauling " << endl;
  }
  oss << "  // If not found then UNKNOWN " << endl;
  oss << "  return MISCIBILITY_SYSTEM_UNKNOWN; // unknown " << endl;
  oss << "} " << endl;
  oss << "// ***************************************************************************" << endl;
  oss << "#endif " << endl;
  oss << " " << endl;
  oss << "// ***************************************************************************" << endl;
  oss << "// *                                                                         *" << endl;
  oss << "// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *" << endl;
  oss << "// *                                                                         *" << endl;
  oss << "// ***************************************************************************" << endl;
  // done
  return oss.str();
}

// **********************************************************************************************
// APENNSY_Parameters::APENNSY_Miscibility_Miedema
// **********************************************************************************************
string APENNSY_Parameters::APENNSY_Miscibility_Miedema(_aflags &aflags) {
  stringstream oss;
  // ************************************************************************
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag")) {;}; // phony
  // ************************************************************************
  oss << "// ***************************************************************************" << endl;
  oss << "// *                                                                         *" << endl;
  oss << "// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *" << endl;
  oss << "// *                                                                         *" << endl;
  oss << "// ***************************************************************************" << endl;
  // HERE
  for(uint k=0;k<alloys.size();k++) {
    oss << "  " << aurostd::PaddedPOST(speciesAB.at(k),4) << " = " ;
    oss << APENNSY_MiscibilityString(Alloy2MiscibilityHT.at(k)) << " ";           // WRITE MISCIBILITY Alloy2MiscibilityHT.at(k)
    oss << APENNSY_MiscibilityString(Alloy2MiscibilityExperiments.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityExperiments.at(k)
    oss << APENNSY_MiscibilityString(Alloy2MiscibilityMiedema.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityMiedema.at(k)
    oss << APENNSY_MiscibilityString(Alloy2MiscibilityHumeRothery.at(k)) << " ";  // WRITE MISCIBILITY Alloy2MiscibilityHumeRothery.at(k)
    // oss << " // from the Miedema file " << endl;
    oss << endl;
  }
  // done
  return oss.str();
}

// **********************************************************************************************
// APENNSY_Parameters::APENNSY_Miscibility_HumeRothery
// **********************************************************************************************
string APENNSY_Parameters::APENNSY_Miscibility_HumeRothery(_aflags &aflags) {
  stringstream oss;
  // ************************************************************************
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag")) {;}; // phony
  // ************************************************************************
  oss << "// ***************************************************************************" << endl;
  oss << "// *                                                                         *" << endl;
  oss << "// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *" << endl;
  oss << "// *                                                                         *" << endl;
  oss << "// ***************************************************************************" << endl;
  for(uint k=0;k<alloys.size();k++) {
    oss << "  " << speciesAB.at(k) << "=" ;
    oss << APENNSY_MiscibilityString(Alloy2MiscibilityHT.at(k)) << " ";           // WRITE MISCIBILITY Alloy2MiscibilityHT.at(k)
    oss << APENNSY_MiscibilityString(Alloy2MiscibilityExperiments.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityExperiments.at(k)
    oss << APENNSY_MiscibilityString(Alloy2MiscibilityMiedema.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityMiedema.at(k)
    oss << APENNSY_MiscibilityString(Alloy2MiscibilityHumeRothery.at(k)) << " ";  // WRITE MISCIBILITY Alloy2MiscibilityHumeRothery.at(k)
    oss << " // from the Miedema file " << endl;
  }
  // done
  return oss.str();
}

// **********************************************************************************************
// APENNSY_Parameters::APENNSY_Miscibility_Table
// **********************************************************************************************
string APENNSY_Parameters::APENNSY_Miscibility_Table(_aflags &aflags) {
  stringstream oss;
  // ************************************************************************
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag")) {;}; // phony
  // ************************************************************************
  oss << "%% ***************************************************************************" << endl;
  oss << "%% *                                                                         *" << endl;
  oss << "%% *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *" << endl;
  oss << "%% *                                                                         *" << endl;
  oss << "%% ***************************************************************************" << endl;
  oss << "\\begin{table}[htb]" << endl;
  oss << "  \\caption{Comparison chart (Ohad, we put all of them and then I focus on the disagreement) } \\label{table1}" << endl;
  oss << "  \\scriptsize" << endl;
  oss << "  \\begin{tabular}{c|c|c|c}" << endl;
  oss << "     System & Exper.\\cite{Pauling,Massalski} C/S & HT & Agreement \\\\ \\hline" << endl;
  for(uint k=0;k<alloys.size();k++) {
    if(k==100 || k==201) {
      oss << "  \\scriptsize" << endl;
      oss << "  \\begin{tabular}{c|c|c|c}" << endl;
      oss << "     System & Exper.\\cite{Pauling,Massalski} C/S & HT & Agreement \\\\ \\hline" << endl;
    }
    oss << "      " << speciesAB.at(k) << " ";
    if(Alloy2MiscibilityExperiments.at(k)==MISCIBILITY_SYSTEM_MISCIBLE) oss << " & C";
    if(Alloy2MiscibilityExperiments.at(k)==MISCIBILITY_SYSTEM_NOMIX) oss << " & S";
    if(Alloy2MiscibilityExperiments.at(k)==MISCIBILITY_SYSTEM_UNKNOWN) oss << " & ?";
    if(Alloy2MiscibilityHT.at(k)==MISCIBILITY_SYSTEM_MISCIBLE) oss << " & C";
    if(Alloy2MiscibilityHT.at(k)==MISCIBILITY_SYSTEM_NOMIX) oss << " & S";
    if(Alloy2MiscibilityHT.at(k)==MISCIBILITY_SYSTEM_UNKNOWN) oss << " & ?";
    oss << "  \\\\ \\hline" << endl;
    if(k==200) {
      oss << "   \\end{tabular}" << endl;
      oss << "\\end{table} " << endl;
    }
  }
  oss << "   \\end{tabular}" << endl;
  oss << "\\end{table} " << endl;
  oss << "%% ***************************************************************************" << endl;
  oss << "%% *                                                                         *" << endl;
  oss << "%% *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *" << endl;
  oss << "%% *                                                                         *" << endl;
  oss << "%% ***************************************************************************" << endl;
  // done
  return oss.str();
}

// **********************************************************************************************
// APENNSY_Parameters::APENNSY_Miscibility_Statistics
// **********************************************************************************************
string APENNSY_Parameters::APENNSY_Miscibility_Statistics(_aflags &aflags) {
  stringstream oss;
  bool PRINT_SYSTEMS=FALSE;
  // ************************************************************************
  if(aflags.vflag.flag("APENNSY::VERBOSE_flag")) {;}; // phony
  oss << endl;
  // ************************************************************************
  // oss << "%% ***************************************************************************" << endl;
  // oss << "%% *                                                                         *" << endl;
  // oss << "%% *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *" << endl;
  // oss << "%% *                                                                         *" << endl;
  // oss << "%% ***************************************************************************" << endl;
  //#define MISCIBILITY_SYSTEM_NOT_STUDIED
  //#define MISCIBILITY_SYSTEM_SOLUTION
  //#define MISCIBILITY_SYSTEM_MISCIBLE
  //#define MISCIBILITY_SYSTEM_UNKNOWN
  //#define MISCIBILITY_SYSTEM_NOMIX
  //#define MISCIBILITY_SYSTEM_CUTOFF
  // Done vector ----------------------------------
  vector<bool> done(alloys.size()+1);
  // MIX-MIX-MIX  ----------------------------------
  oss << "MIX-MIX-MIX" << endl;
  uint cnt=0,cnt1=0,cnt2=0,cnt3=0;
  for(uint k=0;k<alloys.size();k++) {
    if(Alloy2MiscibilityHT.at(k)==MISCIBILITY_SYSTEM_MISCIBLE) cnt1++;
    if(Alloy2MiscibilityExperiments.at(k)==MISCIBILITY_SYSTEM_MISCIBLE) cnt2++;
    if(Alloy2MiscibilityMiedema.at(k)==MISCIBILITY_SYSTEM_MISCIBLE) cnt3++;
  } // done
  oss << " " << cnt1 << " " << cnt2 << " " << cnt3 << " " << endl;
  // -----------------------------------------------------------------------------------------
  oss << "***************************************************************************" << endl;
  oss << "HT VS EXPERIMENTS" << endl;
  oss << "***************************************************************************" << endl;
  // -------------------------------------------
  // clear
  for(uint k=0;k<alloys.size();k++) done.at(k)=FALSE;
  // MIX-MIX  ----------------------------------
  oss << "MIX-MIX-\?\?\?" << endl;
  for(uint k=0;k<alloys.size();k++) { if(k==0) cnt=0;
    if(done.at(k)==FALSE && Alloy2MiscibilityHT.at(k)==MISCIBILITY_SYSTEM_MISCIBLE && Alloy2MiscibilityExperiments.at(k)==MISCIBILITY_SYSTEM_MISCIBLE) {
      done.at(k)=TRUE;cnt++;
      // MISCIBILITY
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHT.at(k)) << " ";           // WRITE MISCIBILITY Alloy2MiscibilityHT.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityExperiments.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityExperiments.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityMiedema.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityMiedema.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHumeRothery.at(k)) << " ";  // WRITE MISCIBILITY Alloy2MiscibilityHumeRothery.at(k)
      // if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(alloys.at(k),11," ");
      if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(speciesAB.at(k),6," ");
      if(PRINT_SYSTEMS) oss << "  " << (int) cnt;
      if(PRINT_SYSTEMS) oss << "  " << (int) ZLibrary.at(k).size(); // << Hlibrary(k) << " ";
      if(PRINT_SYSTEMS) oss << endl;
    }
  } // done
  oss << "  " << (int) cnt << endl;
  // NOMIX-NOMIX  ----------------------------------
  oss << "NOMIX-NOMIX-\?\?\?" << endl;
  for(uint k=0;k<alloys.size();k++) { if(k==0) cnt=0;
    if(done.at(k)==FALSE && Alloy2MiscibilityHT.at(k)==MISCIBILITY_SYSTEM_NOMIX && Alloy2MiscibilityExperiments.at(k)==MISCIBILITY_SYSTEM_NOMIX) {
      done.at(k)=TRUE;cnt++;
      // MISCIBILITY
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHT.at(k)) << " ";           // WRITE MISCIBILITY Alloy2MiscibilityHT.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityExperiments.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityExperiments.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityMiedema.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityMiedema.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHumeRothery.at(k)) << " ";  // WRITE MISCIBILITY Alloy2MiscibilityHumeRothery.at(k)
      // if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(alloys.at(k),11," ");
      if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(speciesAB.at(k),6," ");
      if(PRINT_SYSTEMS) oss << "  " << (int) cnt;
      if(PRINT_SYSTEMS) oss << "  " << (int) ZLibrary.at(k).size(); // << Hlibrary(k) << " ";
      if(PRINT_SYSTEMS) oss << endl;
    }
  } // done
  oss << "  " << (int) cnt << endl;
  // MIX-NOMIX  ----------------------------------
  oss << "MIX-NOMIX-\?\?\?" << endl;
  for(uint k=0;k<alloys.size();k++) { if(k==0) cnt=0;
    if(done.at(k)==FALSE && Alloy2MiscibilityHT.at(k)==MISCIBILITY_SYSTEM_MISCIBLE && Alloy2MiscibilityExperiments.at(k)==MISCIBILITY_SYSTEM_NOMIX) {
      done.at(k)=TRUE;cnt++;
      // MISCIBILITY
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHT.at(k)) << " ";           // WRITE MISCIBILITY Alloy2MiscibilityHT.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityExperiments.at(k)) << " ";  // WRITE MISCIBILITY Alloy2MiscibilityExperiments.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityMiedema.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityMiedema.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHumeRothery.at(k)) << " ";  // WRITE MISCIBILITY Alloy2MiscibilityHumeRothery.at(k)
      // if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(alloys.at(k),11," ");
      if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(speciesAB.at(k),6," ");
      if(PRINT_SYSTEMS) oss << "  " << (int) cnt;
      if(PRINT_SYSTEMS) oss << "  " << (int) ZLibrary.at(k).size(); // << Hlibrary(k) << " ";
      if(PRINT_SYSTEMS) oss << endl;
    }
  } // done
  oss << "  " << (int) cnt << endl;
  // NOMIX-NOMIX  ----------------------------------
  oss << "NOMIX-MIX-\?\?\?" << endl;
  for(uint k=0;k<alloys.size();k++) { if(k==0) cnt=0;
    if(done.at(k)==FALSE && Alloy2MiscibilityHT.at(k)==MISCIBILITY_SYSTEM_NOMIX && Alloy2MiscibilityExperiments.at(k)==MISCIBILITY_SYSTEM_MISCIBLE) {
      done.at(k)=TRUE;cnt++;
      // MISCIBILITY
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHT.at(k)) << " ";           // WRITE MISCIBILITY Alloy2MiscibilityHT.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityExperiments.at(k)) << " ";  // WRITE MISCIBILITY Alloy2MiscibilityExperiments.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityMiedema.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityMiedema.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHumeRothery.at(k)) << " ";  // WRITE MISCIBILITY Alloy2MiscibilityHumeRothery.at(k)
      // if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(alloys.at(k),11," ");
      if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(speciesAB.at(k),6," ");
      if(PRINT_SYSTEMS) oss << "  " << (int) cnt;
      if(PRINT_SYSTEMS) oss << "  " << (int) ZLibrary.at(k).size(); // << Hlibrary(k) << " ";
      if(PRINT_SYSTEMS) oss << endl;
    }
  } // done
  oss << "  " << (int) cnt << endl;
  // LEFT OVER  ----------------------------------
  oss << "LEFT OVER" << endl;
  for(uint k=0;k<alloys.size();k++) { if(k==0) cnt=0;
    if(done.at(k)==FALSE) {
      done.at(k)=TRUE;cnt++;
      // MISCIBILITY
      oss << APENNSY_MiscibilityString(Alloy2MiscibilityHT.at(k)) << " ";           // WRITE MISCIBILITY Alloy2MiscibilityHT.at(k)
      oss << APENNSY_MiscibilityString(Alloy2MiscibilityExperiments.at(k)) << " ";  // WRITE MISCIBILITY Alloy2MiscibilityExperiments.at(k)
      oss << APENNSY_MiscibilityString(Alloy2MiscibilityMiedema.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityMiedema.at(k)
      oss << APENNSY_MiscibilityString(Alloy2MiscibilityHumeRothery.at(k)) << " ";  // WRITE MISCIBILITY Alloy2MiscibilityHumeRothery.at(k)
      // oss << "  " << aurostd::PaddedPOST(alloys.at(k),11," ");
      oss << "  " << aurostd::PaddedPOST(speciesAB.at(k),6," ");
      oss << "  " << (int) cnt;
      oss << "  " << (int) ZLibrary.at(k).size(); // << Hlibrary(k) << " ";
      oss << endl;
    }
  } // done
  oss << "  " << (int) cnt << endl;
  // -----------------------------------------------------------------------------------------
  oss << "***************************************************************************" << endl;
  oss << "HT VS MIEDEMA" << endl;
  oss << "***************************************************************************" << endl;
  // clear ----------------------------------------
  for(uint k=0;k<alloys.size();k++) done.at(k)=FALSE;
  // MIX-MIX  ----------------------------------
  oss << "MIX-\?\?\?-MIX" << endl;
  for(uint k=0;k<alloys.size();k++) { if(k==0) cnt=0;
    if(done.at(k)==FALSE && Alloy2MiscibilityHT.at(k)==MISCIBILITY_SYSTEM_MISCIBLE && Alloy2MiscibilityMiedema.at(k)==MISCIBILITY_SYSTEM_MISCIBLE) {
      done.at(k)=TRUE;cnt++;
      // MISCIBILITY
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHT.at(k)) << " ";           // WRITE MISCIBILITY Alloy2MiscibilityHT.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityExperiments.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityExperiments.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityMiedema.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityMiedema.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHumeRothery.at(k)) << " ";  // WRITE MISCIBILITY Alloy2MiscibilityHumeRothery.at(k)
      // if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(alloys.at(k),11," ");
      if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(speciesAB.at(k),6," ");
      if(PRINT_SYSTEMS) oss << "  " << (int) cnt;
      if(PRINT_SYSTEMS) oss << "  " << (int) ZLibrary.at(k).size(); // << Hlibrary(k) << " ";
      if(PRINT_SYSTEMS) oss << endl;
    }
  } // done
  oss << "  " << (int) cnt << endl;
  // NOMIX-NOMIX  ----------------------------------
  oss << "NOMIX-\?\?\?-NOMIX" << endl;
  for(uint k=0;k<alloys.size();k++) { if(k==0) cnt=0;
    if(done.at(k)==FALSE && Alloy2MiscibilityHT.at(k)==MISCIBILITY_SYSTEM_NOMIX && Alloy2MiscibilityMiedema.at(k)==MISCIBILITY_SYSTEM_NOMIX) {
      done.at(k)=TRUE;cnt++;
      // MISCIBILITY
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHT.at(k)) << " ";           // WRITE MISCIBILITY Alloy2MiscibilityHT.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityExperiments.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityExperiments.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityMiedema.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityMiedema.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHumeRothery.at(k)) << " ";  // WRITE MISCIBILITY Alloy2MiscibilityHumeRothery.at(k)
      // if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(alloys.at(k),11," ");
      if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(speciesAB.at(k),6," ");
      if(PRINT_SYSTEMS) oss << "  " << (int) cnt;
      if(PRINT_SYSTEMS) oss << "  " << (int) ZLibrary.at(k).size(); // << Hlibrary(k) << " ";
      if(PRINT_SYSTEMS) oss << endl;
    }
  } // done
  oss << "  " << (int) cnt << endl;
  // MIX-NOMIX  ----------------------------------
  oss << "MIX-\?\?\?-NOMIX" << endl;
  for(uint k=0;k<alloys.size();k++) { if(k==0) cnt=0;
    if(done.at(k)==FALSE && Alloy2MiscibilityHT.at(k)==MISCIBILITY_SYSTEM_MISCIBLE && Alloy2MiscibilityMiedema.at(k)==MISCIBILITY_SYSTEM_NOMIX) {
      done.at(k)=TRUE;cnt++;
      // MISCIBILITY
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHT.at(k)) << " ";           // WRITE MISCIBILITY Alloy2MiscibilityHT.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityExperiments.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityExperiments.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityMiedema.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityMiedema.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHumeRothery.at(k)) << " ";  // WRITE MISCIBILITY Alloy2MiscibilityHumeRothery.at(k)
      // if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(alloys.at(k),11," ");
      if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(speciesAB.at(k),6," ");
      if(PRINT_SYSTEMS) oss << "  " << (int) cnt;
      if(PRINT_SYSTEMS) oss << "  " << (int) ZLibrary.at(k).size(); // << Hlibrary(k) << " ";
      if(PRINT_SYSTEMS) oss << endl;
    }
  } // done
  oss << "  " << (int) cnt << endl;
  // NOMIX-NOMIX  ----------------------------------
  oss << "NOMIX-\?\?\?-MIX" << endl;
  for(uint k=0;k<alloys.size();k++) { if(k==0) cnt=0;
    if(done.at(k)==FALSE && Alloy2MiscibilityHT.at(k)==MISCIBILITY_SYSTEM_NOMIX && Alloy2MiscibilityMiedema.at(k)==MISCIBILITY_SYSTEM_MISCIBLE) {
      done.at(k)=TRUE;cnt++;
      // MISCIBILITY
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHT.at(k)) << " ";           // WRITE MISCIBILITY Alloy2MiscibilityHT.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityExperiments.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityExperiments.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityMiedema.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityMiedema.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHumeRothery.at(k)) << " ";  // WRITE MISCIBILITY Alloy2MiscibilityHumeRothery.at(k)
      // if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(alloys.at(k),11," ");
      if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(speciesAB.at(k),6," ");
      if(PRINT_SYSTEMS) oss << "  " << (int) cnt;
      if(PRINT_SYSTEMS) oss << "  " << (int) ZLibrary.at(k).size(); // << Hlibrary(k) << " ";
      if(PRINT_SYSTEMS) oss << endl;
    }
  } // done
  oss << "  " << (int) cnt << endl;
  // LEFT OVER  ----------------------------------
  oss << "LEFT OVER" << endl;
  for(uint k=0;k<alloys.size();k++) { if(k==0) cnt=0;
    if(done.at(k)==FALSE) {
      done.at(k)=TRUE;cnt++;
      // MISCIBILITY
      oss << APENNSY_MiscibilityString(Alloy2MiscibilityHT.at(k)) << " ";           // WRITE MISCIBILITY Alloy2MiscibilityHT.at(k)
      oss << APENNSY_MiscibilityString(Alloy2MiscibilityExperiments.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityExperiments.at(k)
      oss << APENNSY_MiscibilityString(Alloy2MiscibilityMiedema.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityMiedema.at(k)
      oss << APENNSY_MiscibilityString(Alloy2MiscibilityHumeRothery.at(k)) << " ";  // WRITE MISCIBILITY Alloy2MiscibilityHumeRothery.at(k)
      // oss << "  " << aurostd::PaddedPOST(alloys.at(k),11," ");
      oss << "  " << aurostd::PaddedPOST(speciesAB.at(k),6," ");
      oss << "  " << (int) cnt;
      oss << "  " << (int) ZLibrary.at(k).size(); // << Hlibrary(k) << " ";
      oss << endl;
    }
  } // done
  oss << "  " << (int) cnt << endl;
  // done
  // -----------------------------------------------------------------------------------------
  oss << "***************************************************************************" << endl;
  oss << "PAULING VS MIEDEMA" << endl;
  oss << "***************************************************************************" << endl;
  // clear ----------------------------------------
  for(uint k=0;k<alloys.size();k++) done.at(k)=FALSE;
  // MIX-MIX  ----------------------------------
  oss << "\?\?\?-MIX-MIX" << endl;
  for(uint k=0;k<alloys.size();k++) { if(k==0) cnt=0;
    if(done.at(k)==FALSE && Alloy2MiscibilityExperiments.at(k)==MISCIBILITY_SYSTEM_MISCIBLE && Alloy2MiscibilityMiedema.at(k)==MISCIBILITY_SYSTEM_MISCIBLE) {
      done.at(k)=TRUE;cnt++;
      // MISCIBILITY
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHT.at(k)) << " ";           // WRITE MISCIBILITY Alloy2MiscibilityHT.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityExperiments.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityExperiments.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityMiedema.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityMiedema.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHumeRothery.at(k)) << " ";  // WRITE MISCIBILITY Alloy2MiscibilityHumeRothery.at(k)
      // if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(alloys.at(k),11," ");
      if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(speciesAB.at(k),6," ");
      if(PRINT_SYSTEMS) oss << "  " << (int) cnt;
      if(PRINT_SYSTEMS) oss << "  " << (int) ZLibrary.at(k).size(); // << Hlibrary(k) << " ";
      if(PRINT_SYSTEMS) oss << endl;
    }
  } // done
  oss << "  " << (int) cnt << endl;
  // NOMIX-NOMIX  ----------------------------------
  oss << "\?\?\?-NOMIX--NOMIX" << endl;
  for(uint k=0;k<alloys.size();k++) { if(k==0) cnt=0;
    if(done.at(k)==FALSE && Alloy2MiscibilityExperiments.at(k)==MISCIBILITY_SYSTEM_NOMIX && Alloy2MiscibilityMiedema.at(k)==MISCIBILITY_SYSTEM_NOMIX) {
      done.at(k)=TRUE;cnt++;
      // MISCIBILITY
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHT.at(k)) << " ";           // WRITE MISCIBILITY Alloy2MiscibilityHT.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityExperiments.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityExperiments.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityMiedema.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityMiedema.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHumeRothery.at(k)) << " ";  // WRITE MISCIBILITY Alloy2MiscibilityHumeRothery.at(k)
      // if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(alloys.at(k),11," ");
      if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(speciesAB.at(k),6," ");
      if(PRINT_SYSTEMS) oss << "  " << (int) cnt;
      if(PRINT_SYSTEMS) oss << "  " << (int) ZLibrary.at(k).size(); // << Hlibrary(k) << " ";
      if(PRINT_SYSTEMS) oss << endl;
    }
  } // done
  oss << "  " << (int) cnt << endl;
  // MIX-NOMIX  ----------------------------------
  oss << "\?\?\?-MIX--NOMIX" << endl;
  for(uint k=0;k<alloys.size();k++) { if(k==0) cnt=0;
    if(done.at(k)==FALSE && Alloy2MiscibilityExperiments.at(k)==MISCIBILITY_SYSTEM_MISCIBLE && Alloy2MiscibilityMiedema.at(k)==MISCIBILITY_SYSTEM_NOMIX) {
      done.at(k)=TRUE;cnt++;
      // MISCIBILITY
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHT.at(k)) << " ";           // WRITE MISCIBILITY Alloy2MiscibilityHT.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityExperiments.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityExperiments.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityMiedema.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityMiedema.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHumeRothery.at(k)) << " ";  // WRITE MISCIBILITY Alloy2MiscibilityHumeRothery.at(k)
      // if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(alloys.at(k),11," ");
      if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(speciesAB.at(k),6," ");
      if(PRINT_SYSTEMS) oss << "  " << (int) cnt;
      if(PRINT_SYSTEMS) oss << "  " << (int) ZLibrary.at(k).size(); // << Hlibrary(k) << " ";
      if(PRINT_SYSTEMS) oss << endl;
    }
  } // done
  oss << "  " << (int) cnt << endl;
  // NOMIX-NOMIX  ----------------------------------
  oss << "\?\?\?-NOMIX--MIX" << endl;
  for(uint k=0;k<alloys.size();k++) { if(k==0) cnt=0;
    if(done.at(k)==FALSE && Alloy2MiscibilityExperiments.at(k)==MISCIBILITY_SYSTEM_NOMIX && Alloy2MiscibilityMiedema.at(k)==MISCIBILITY_SYSTEM_MISCIBLE) {
      done.at(k)=TRUE;cnt++;
      // MISCIBILITY
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHT.at(k)) << " ";           // WRITE MISCIBILITY Alloy2MiscibilityHT.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityExperiments.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityExperiments.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityMiedema.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityMiedema.at(k)
      if(PRINT_SYSTEMS) oss << APENNSY_MiscibilityString(Alloy2MiscibilityHumeRothery.at(k)) << " ";  // WRITE MISCIBILITY Alloy2MiscibilityHumeRothery.at(k)
      // if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(alloys.at(k),11," ");
      if(PRINT_SYSTEMS) oss << "  " << aurostd::PaddedPOST(speciesAB.at(k),6," ");
      if(PRINT_SYSTEMS) oss << "  " << (int) cnt;
      if(PRINT_SYSTEMS) oss << "  " << (int) ZLibrary.at(k).size(); // << Hlibrary(k) << " ";
      if(PRINT_SYSTEMS) oss << endl;
    }
  } // done
  oss << "  " << (int) cnt << endl;
  // done
  // LEFT OVER  ----------------------------------
  oss << "LEFT OVER" << endl;
  for(uint k=0;k<alloys.size();k++) { if(k==0) cnt=0;
    if(done.at(k)==FALSE) {
      done.at(k)=TRUE;cnt++;
      // MISCIBILITY
      oss << APENNSY_MiscibilityString(Alloy2MiscibilityHT.at(k)) << " ";           // WRITE MISCIBILITY Alloy2MiscibilityHT.at(k)
      oss << APENNSY_MiscibilityString(Alloy2MiscibilityExperiments.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityExperiments.at(k)
      oss << APENNSY_MiscibilityString(Alloy2MiscibilityMiedema.at(k)) << " ";      // WRITE MISCIBILITY Alloy2MiscibilityMiedema.at(k)
      oss << APENNSY_MiscibilityString(Alloy2MiscibilityHumeRothery.at(k)) << " ";  // WRITE MISCIBILITY Alloy2MiscibilityHumeRothery.at(k)
      // oss << "  " << aurostd::PaddedPOST(alloys.at(k),11," ");
      oss << "  " << aurostd::PaddedPOST(speciesAB.at(k),6," ");
      oss << "  " << (int) cnt;
      oss << "  " << (int) ZLibrary.at(k).size(); // << Hlibrary(k) << " ";
      oss << endl;
    }
  } // done
  oss << "  " << (int) cnt << endl;
  // -------------------------------------------------------------------------
  return oss.str();
}

// ***************************************************************************
// ***************************************************************************

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
