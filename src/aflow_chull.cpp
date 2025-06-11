// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *           Aflow COREY OSES - Duke University 2013-2021                  *
// *                                                                         *
// ***************************************************************************
// Written by Corey Oses
// corey.oses@duke.edu
// Previous versions also written by Eric Perim and Eric Gossett

#ifndef _AFLOW_CHULL_CPP_
#define _AFLOW_CHULL_CPP_

#include "aflow.h"
#include "aflow_chull.h"
#include "aflow_compare_structure.h"

//[CO20220630 - OBSOLETE]// Some parts are written within the C++0x support in GCC, especially std::thread,
//[CO20220630 - OBSOLETE]// which is implemented in gcc 4.4 and higher. For multithreads with std::thread see:
//[CO20220630 - OBSOLETE]// http://www.justsoftwaresolutions.co.uk/threading/multithreading-in-c++0x-part-1-starting-threads.html
//[CO20220630 - OBSOLETE]#if GCC_VERSION >= 40400
//[CO20220630 - OBSOLETE]#define AFLOW_CHULL_MULTITHREADS_ENABLE
//[CO20220630 - OBSOLETE]#include <thread>
//[CO20220630 - OBSOLETE]#include <mutex>
//[CO20220630 - OBSOLETE]static std::mutex m;
//[CO20220630 - OBSOLETE]#else
//[CO20220630 - OBSOLETE]#warning "The multithread parts of APL will not be included, since they need gcc 4.4 and higher (C++0x support)."
//[CO20220630 - OBSOLETE]#endif

#define _DEBUG_CHULL_ false  //CO20190116

#define _AFLOW_CHULL_PRINT_LOGO_1 TRUE
#define _AFLOW_CHULL_PRINT_LOGO_2 FALSE

// DEFINITIONS
const std::string AFLOW_WEB = string("http://" + AFLOWLIB_MATERIALS_SERVER);
const std::string NOMAD_WEB = string("http://www.nomad-coe.eu/");
const std::string ENTRY_PAGE_URL_PREFIX = string(AFLOW_WEB + "/material/?id=");  //"/material.php?id="
const std::string AFLOW_HULL_ENDPOINT_STRING = "aflow_hull_endpoint";
const std::string LATEX_DEFAULT_COLORS = "blue,green,red,brown,cyan,lime,magenta,olive,orange,pink,purple,teal,violet,darkgray,gray,lightgray,black,white,yellow"; //putting favorites near beginning, not so good ones in the back
const std::string LATEX_DVIPS_COLORS = "Apricot,Aquamarine,Bittersweet,Black,Blue,BlueGreen,BlueViolet,BrickRed,Brown,BurntOrange,CadetBlue,CarnationPink,Cerulean,CornflowerBlue,Cyan,Dandelion,DarkOrchid,Emerald,ForestGreen,Fuchsia,Goldenrod,Gray,GreenGreen,Yellow,JungleGreen,Lavender,LimeGreen,Magenta,Mahogany,Maroon,Melon,MidnightBlue,Mulberry,NavyBlue,OliveGreen,Orange,OrangeRed,Orchid,Peach,Periwinkle,PineGreen,Plum,ProcessBlue,Purple,RawSienna,Red,RedOrange,RedViolet,Rhodamine,RoyalBlue,RoyalPurple,RubineRed,Salmon,SeaGreen,Sepia,SkyBlue,SpringGreen,Tan,TealBlue,Thistle,Turquoise,Violet,VioletRed,White,WildStrawberry,Yellow,YellowGreen,YellowOrange";  //contains defaults and more
const std::string LATEX_COLORS_TO_AVOID = "black,white,yellow,darkgray,gray,lightgray";

// USAGE FLAGS
//[CO20180819 - MOVED TO AFLOWRC]const bool IGNORE_BAD_DATABASE = true;        //skip bad entries
const bool CORRECT_BAD_DATABASE = true;                                        //make minor corrections, carried over from apennsy (SC)
const bool PRINT_DIST2HULL_COL_TEX = false;                                    //print Dist2hull column in tex, there's no need because it's not used for anything in the image
const bool GET_DECOMPOSITION_POLYMORPHS = true;                                //print decomposition information for polymorphs

// LATEX PRINTING MODES
const char ADDPLOT_MODE_HULL_POINTS = 'P';
const char ADDPLOT_MODE_OFF_HULL_POINTS = 'O';
const char ADDPLOT_MODE_HULL_FACETS = 'F';
const char ADDPLOT_MODE_HULL_FACETSDROP_SHADOWS = 'D';
const char ADDPLOT_MODE_HEAT_MAPS = 'H';

// CITATION INFO
const std::string CHULL_CITE = "Use of this data welcomes reference to the following publication:"; //Use of this data welcomes reference to the following publication: //cite:
const std::string CHULL_AUTHORS = "C. Oses, E. Gossett, D. Hicks, F. Rose, M. J. Mehl, E. Perim, I. Takeuchi, S. Sanvito, M. Scheffler, Y. Lederer, O. Levy, C. Toher, and S. Curtarolo";
const std::string CHULL_TITLE = "AFLOW-CHULL: Cloud-Oriented Platform for Autonomous Phase Stability Analysis";
const std::string CHULL_JOURNAL_LATEX = "J. Chem. Inf. Model. \\textbf{58}(12), 2477-2490 (2018). \\href{https://doi.org/10.1021/acs.jcim.8b00393}{doi:10.1021/acs.jcim.8b00393}";
const std::string CHULL_JOURNAL_TXT = "J. Chem. Inf. Model. 58(12), 2477-2490 (2018). doi:10.1021/acs.jcim.8b00393";
//const std::string CHULL_JOURNAL_LATEX = "submitted \\href{https://arxiv.org/abs/1806.06901}{arXiv:1806.06901} (2018)";
//const std::string CHULL_JOURNAL_TXT = "submitted arXiv:1806.06901 (2018)";
//const std::string CHULL_JOURNAL = "submitted to Comput. Mater. Sci. (2018)";

// LATEX
const double LATEX_PT2INCH = 100.0/7227.0;  //https://tex.stackexchange.com/questions/8260/what-are-the-various-units-ex-em-in-pt-bp-dd-pc-expressed-in-mm
const double LATEX_WIDTH_LETTER_STD = 8.5;  //inches
const double LATEX_LEFT_MARGIN_LETTER_STD = 0.5;  //inches
const double LATEX_RIGHT_MARGIN_LETTER_STD = 0.5;  //inches
const double LATEX_TABCOLSEP_STD = 6;  //pt
const double LATEX_ARRAYRULEWIDTH_STD = .4;  //pt

//CO20180419 - moved to aurostd::xerror and aurostd::xerror
//namespace chull {
//  aurostd::xerror::aurostd::xerror(__AFLOW_FILE__,const std::string& function,const std::string& message) : std::runtime_error(message),f_name(function) {}  // I/O or computer type errors (no entries loaded)
//  aurostd::xerror::aurostd::xerror(__AFLOW_FILE__,const std::string& function,std::stringstream& message) : std::runtime_error(message.str()),f_name(function) {message.str("");}  // I/O or computer type errors (no entries loaded)
//  string aurostd::xerror::where(){return f_name;}
//  aurostd::xerror::aurostd::xerror(__AFLOW_FILE__,const std::string& function,const std::string& message) : std::logic_error(message),f_name(function) {}    //errors in logic, unintended (and insurmountable) use of functionality
//  aurostd::xerror::aurostd::xerror(__AFLOW_FILE__,const std::string& function,std::stringstream& message) : std::logic_error(message.str()),f_name(function) {message.str("");}    //errors in logic, unintended (and insurmountable) use of functionality
//  string aurostd::xerror::where(){return f_name;}
//} // namespace chull

namespace chull {
  bool convexHull(const aurostd::xoption& _vpflow) {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    ostream& oss = cout;
    ofstream FileMESSAGE;
    stringstream message;
    
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}
    
    aurostd::xoption vpflow(_vpflow);
    string directory=getPath(vpflow,FileMESSAGE,oss);
    vpflow.flag("CHULL::PATH",true);
    vpflow.pop_attached("CHULL::PATH");
    vpflow.push_attached("CHULL::PATH",directory);  //rectify directory in case it changes (to pwd)
    _aflags aflags; aflags.Directory=directory;

    //do this immediately, we don't want JSON-reader to bomb because there's no output
    if(vpflow.flag("CHULL::SCREEN_ONLY")){
      if(!(vpflow.flag("CHULL::TEXT_DOC")||vpflow.flag("CHULL::JSON_DOC"))){vpflow.flag("CHULL::JSON_DOC",true);}
      if(vpflow.flag("CHULL::TEXT_DOC")&&vpflow.flag("CHULL::JSON_DOC")){vpflow.flag("CHULL::TEXT_DOC",false);} //only one output allowed
    }

    //////////////////////////////////////////////////////////////////////////////
    // START Display usage, if requested
    //////////////////////////////////////////////////////////////////////////////

    string usage_usage="aflow --convex_hull=|--chull --alloy=MnPdPt[,AlCuZn,...] [--np=1] [chull_options] [--destination=[DIRECTORY]]";
    vector<string> usage_options;
    usage_options.push_back(usage_usage);
    usage_options.push_back(" ");
    usage_options.push_back("chull_options:");
    usage_options.push_back(" ");
    usage_options.push_back("GENERAL OPTIONS:");
    usage_options.push_back("--usage");
    usage_options.push_back("--print=|--p=|--output=|--o=latex|pdf|png|json|text|jupyter|jupyter2|jupyter3");
    usage_options.push_back("--keep=tex|--keep_tex|--keeptex|--tex");
    usage_options.push_back("--keep=log|--keep_log|--keeplog|--log");
    usage_options.push_back("--keep=tex,log");
    usage_options.push_back(" ");
    usage_options.push_back("LOADING OPTIONS:");
    usage_options.push_back("--load_library=|--loadlibrary=|--ll=icsd|lib1|lib2|lib3");
    usage_options.push_back("--load_API|--load_api|--loadapi|--lapi|--api");
    usage_options.push_back("--load_entries_entry_output|--loadentriesentryoutput|--leo");
    usage_options.push_back("--neglect=|--ban=aflow:bb0d45ab555bc208);aflow:fb9eaa58604ce774");
    usage_options.push_back("--see_neglect|--seeneglect|--sn");
    usage_options.push_back("--remove_extreme_points=|--removeextremepoints=|--remove_extrema=|--removeextrema=|--rep=-1000");
    usage_options.push_back("--entropic_temperature|--entropictemperature|--entroptemp");
    usage_options.push_back(" ");
    usage_options.push_back("ANALYSIS OPTIONS:");
    usage_options.push_back("--distance_to_hull=|--distancetohull=|--distance2hull=|--dist2hull=|--d2h=aflow:bb0d45ab555bc208,aflow:fb9eaa58604ce774");
    usage_options.push_back("--stability_criterion=|--stabilitycriterion=|--stable_criterion=|--scriterion=|--sc=aflow:bb0d45ab555bc208,aflow:fb9eaa58604ce774");
    usage_options.push_back("--n+1_enthalpy_gain=|--=|--n+1enthalpygain=|--n+1energygain=|--n+1egain=|--n1egain=|--n+1_enthalpygain=|--n+1+energygain=|--n+1_egain=|--nplus1=aflow:bb0d45ab555bc208,aflow:fb9eaa58604ce774");
    usage_options.push_back("--hull_formation_enthalpy=|--hull_enthalpy=|--hull_energy=0.25,0.25");
    usage_options.push_back("--skip_structure_comparison|--skipstructruecomparison|--skipstructcomp|--ssc");
    usage_options.push_back("--skip_stability_criterion_analysis|--skipstabilitycriterionanalysis|--skip_scriterion|--skipscriterion|--sscriterion");
    usage_options.push_back("--skip_n_plus_1_enthalpy_gain_analysis|--skip_n_plus_1_energy_gain_analysis|--skipnplus1enthalpygainanalysis|--skipnplus1energygainanalysis|--skip_nplus1|--skipnplus1|--snp1|--snpo");
    usage_options.push_back("--include_skewed_hulls|--include_skewed|--ish");
    usage_options.push_back("--include_unreliable_hulls|--include_unreliable|--iuh");
    usage_options.push_back("--include_outliers|--io");
    usage_options.push_back("--strict_outlier_analysis|--soa");
    usage_options.push_back("--include_ill_converged|--iic");
    usage_options.push_back("--force");
    usage_options.push_back(" ");
    usage_options.push_back("LATEX/PDF/PNG OPTIONS:");
    usage_options.push_back("--image_only|--imageonly|--image");
    usage_options.push_back("--no_document|--nodocument|--no_doc|--nodoc|--full_page_image|--fullpageimage");
    usage_options.push_back("--document_only|--documentonly|--doc_only|--doconly|--doc");
    usage_options.push_back("--latex_output|--latexoutput");
    usage_options.push_back("--latex_interactive|--latexinteractive");
    usage_options.push_back("--light_contrast|--lightcontrast|--lc");
    usage_options.push_back("--large_font|--largefont|--large|--lf");
    usage_options.push_back("--png_resolution=|--pngresolution=|--pngr=300");
    usage_options.push_back("--plot_iso_max_latent_heat|--plot_isomax|--iso_max|--isomax");
    usage_options.push_back(" ");
    // output usage
    if(vpflow.flag("CHULL::USAGE")) {
      if(!vpflow.flag("CHULL::SCREEN_ONLY")){init::MessageOption( "--usage", "CHULL()", usage_options);}
      if(vpflow.flag("CHULL::SCREEN_ONLY")&&vpflow.flag("CHULL::JSON_DOC")){oss << "{}";} //so JSON-reader doesn't bomb
      return TRUE;
    }

    //////////////////////////////////////////////////////////////////////////////
    // END Display usage, if requested
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // START Flag manipulation
    //////////////////////////////////////////////////////////////////////////////

    // address possible flag issues
    // get libraries to load from flags
    vector<string> vlibraries;
    if(!vpflow.flag("PFLOW::LOAD_LIBRARY")) {vpflow.push_attached("PFLOW::LOAD_LIBRARY","all");}
    aurostd::string2tokens(vpflow.getattachedscheme("PFLOW::LOAD_LIBRARY"),vlibraries,",");
    // special case, all
    if((vlibraries.size() == 1) && (vlibraries[0] == "all")) {pflow::defaultLoadEntriesFlags(vpflow, FileMESSAGE, oss, std::string("A"), false, true);}
    else{
      bool found=false;
      for(uint i=0,fl_size_i=vlibraries.size();i<fl_size_i;i++) {
        found = false;
        for(uint lib=1;lib<=_AFLOW_LIB_MAX_ && !found;lib++) {
          string LIB = aurostd::utype2string(lib);
          if(!found && (aurostd::toupper(vlibraries[i]) == "LIB" + LIB)) {
            pflow::defaultLoadEntriesFlags(vpflow, FileMESSAGE, oss, LIB, false, true);
            found = true;
          }
        }
        if(!found && (aurostd::toupper(vlibraries[i]) == "ICSD")) {
          pflow::defaultLoadEntriesFlags(vpflow, FileMESSAGE, oss, "ICSD", false, true);
          found = true;
        }
        //}else{
        if(!found) {
          message << "Incorrect input for loadlibraries \"" << vlibraries[i] << "\"";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
          if(vpflow.flag("CHULL::SCREEN_ONLY")&&vpflow.flag("CHULL::JSON_DOC")){oss << "{}";} //so JSON-reader doesn't bomb
          return FALSE;
        }
    }
    }
    // get desired output from flags
    // give an error if input is not as desired
    if(vpflow.flag("CHULL::OUTPUT")) {
      vector<string> out_forms;
      aurostd::string2tokens(vpflow.getattachedscheme("CHULL::OUTPUT"), out_forms, ",");
      for(uint i=0,fl_size_i=out_forms.size();i<fl_size_i;i++) {
        if(!((out_forms[i][0] == 'A' || out_forms[i][0] == 'a') ||
              (out_forms[i][0] == 'F' || out_forms[i][0] == 'f') ||
              (out_forms[i][0] == 'T' || out_forms[i][0] == 't') ||
              (out_forms[i][0] == 'J' || out_forms[i][0] == 'j') ||
              (out_forms[i][0] == 'W' || out_forms[i][0] == 'w') ||
              (out_forms[i][0] == 'L' || out_forms[i][0] == 'l' || out_forms[i][0] == 'P' || out_forms[i][0] == 'p'))) {
          message << "Incorrect input for output \"" << out_forms[i] << "\"";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
          if(vpflow.flag("CHULL::SCREEN_ONLY")&&vpflow.flag("CHULL::JSON_DOC")){oss << "{}";} //so JSON-reader doesn't bomb
          return FALSE;
        }
      }
    }

    //////////////////////////////////////////////////////////////////////////////
    // END Flag manipulation
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // START Gathering hull inputs
    //////////////////////////////////////////////////////////////////////////////

    message << aflow::Banner("BANNER_NORMAL");
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_RAW_);  //first to screen (not logged, file not opened)
    message << "Processing inputs. Eliminating degenerate (duplicate-element) inputs and duplicates.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);

    // get elements input
    vector<string> vinputs, velements;
    string inputs = vpflow.getattachedscheme("PFLOW::ALLOY");
    if(inputs.empty()) {
      message << "No input given for elements";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      //really drive the point home
      if(!vpflow.flag("CHULL::SCREEN_ONLY")){init::ErrorOption( "--alloy=" + vpflow.getattachedscheme("PFLOW::ALLOY"), "PFLOW()", usage_options);}
      if(vpflow.flag("CHULL::SCREEN_ONLY")&&vpflow.flag("CHULL::JSON_DOC")){oss << "{}";} //so JSON-reader doesn't bomb
      return FALSE;
    }
    inputs=aurostd::RemoveWhiteSpaces(inputs);  //CO20200531 - sometimes web injects spaces
    aurostd::StringSubst(inputs,"-","");;  //CO20200531 - removing '-' from web
    //
    vector<string> tokens_comma;
    uint nary;
    string original_input;
    bool verbose_elimination=true;
    if(aurostd::substring2bool(inputs,":")){  //enumeration mode
      vector<string> tokens_colon;
      aurostd::string2tokens(inputs,tokens_colon,":");
      nary=tokens_colon.size();
      vector<vector<string> > vcomponents;  //< <Ag,Au>,<Mn> >
      vector<int> vsizes;
      for(uint i=0,fl_size_i=tokens_colon.size();i<fl_size_i;i++){
        vcomponents.push_back(vector<string>(0));
        if(tokens_colon[i]=="LIB2"){tokens_colon[i]=SPECIE_RAW_LIB2;verbose_elimination=false;}
        if(tokens_colon[i]=="LIB2U"){tokens_colon[i]=SPECIE_RAW_LIB2U;verbose_elimination=false;}
        if(tokens_colon[i]=="LIB3"){tokens_colon[i]=SPECIE_RAW_LIB3;verbose_elimination=false;}
        aurostd::string2tokens(tokens_colon[i],tokens_comma,",");
        for(uint j=0,fl_size_j=tokens_comma.size();j<fl_size_j;j++){vcomponents[i].push_back(tokens_comma[j]);}
        vsizes.push_back(vcomponents[i].size());
      }
      aurostd::xcombos xc(vsizes,'E');
      while(xc.increment()){
        velements=xc.applyCombo(vcomponents);
        if(verbose_elimination){velements=aurostd::getElements(aurostd::joinWDelimiter(velements,""),pp_string,FileMESSAGE,true,true,false,oss);}  //fast way to eliminate pseudopotential information  //clean and sort, do not keep_pp  //CO20190712
        //[CO20190712 - OBSOLETE]if(verbose_elimination){velements=pflow::getAlphabeticVectorString(aurostd::joinWDelimiter(velements,""),FileMESSAGE,oss);}  //fast way to eliminate pseudopotential information
        std::sort(velements.begin(),velements.end());
        original_input=aurostd::joinWDelimiter(velements,"");
        std::sort(velements.begin(),velements.end());velements.erase( std::unique( velements.begin(), velements.end() ), velements.end() );  //MnMnPd is same as MnPd
        if(velements.size()==nary){vinputs.push_back(aurostd::joinWDelimiter(velements,""));}  //ignore MnPd if requested ternaries
        else{
          if(verbose_elimination){
            message << "Ignoring degenerate (duplicate-element) input (" << original_input << ")";
            pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
          }
        }
      }
    }else{  //simple list mode
      if(inputs=="LIB3HEUSLERS"&&aurostd::FileExist("/common/LIB3/LIB")){
        //[CO20221112 - might have ':', messing up combination input below]inputs=aurostd::execute2string("ls /common/LIB3/LIB | tr -s '\\r\\n' ',' | sed -e 's/,$/\\n/'");
        //[CO20221112 - too slow]inputs=aurostd::joinWDelimiter(aurostd::string2vectorstring(aurostd::execute2string("ls /common/LIB3/LIB")),",");
        inputs=aurostd::execute2string("find /common/LIB3/LIB -maxdepth 2 -name 'T000*' | grep 'T0001\\|T0002\\|T0003' | xargs dirname | sed 's/\\/common\\/LIB3\\/LIB\\///' | sort | uniq | tr -s '\\r\\n' ',' | sed -e 's/,$/\\n/'");
        //[CO20221112 - too slow]velements.clear();
        //[CO20221112 - too slow]for(uint i=0;i<vlib3.size();i++){
        //[CO20221112 - too slow]  velements=aurostd::getElements(vlib3[i],pp_string,FileMESSAGE,true,true,false,oss);
        //[CO20221112 - too slow]  vlib3[i]=aurostd::joinWDelimiter(velements,"");
        //[CO20221112 - too slow]}
        //[CO20221112 - too slow]inputs=aurostd::joinWDelimiter(vlib3,",");
        //[CO20221112 - too slow]velements.clear();
      }
      //
      aurostd::string2tokens(inputs, tokens_comma, ",");
      for(uint i=0,fl_size_i=tokens_comma.size();i<fl_size_i;i++){
        velements=aurostd::getElements(tokens_comma[i],pp_string,FileMESSAGE,true,true,false,oss);  //clean and sort, do not keep_pp  //CO20190712
        //[CO20190712 - OBSOLETE]velements=pflow::getAlphabeticVectorString(tokens_comma[i], FileMESSAGE,oss);
        nary=velements.size();
        original_input=aurostd::joinWDelimiter(velements,"");
        std::sort(velements.begin(),velements.end());velements.erase( std::unique( velements.begin(), velements.end() ), velements.end() );  //MnMnPd is same as MnPd
        if(velements.size()==nary){vinputs.push_back(aurostd::joinWDelimiter(velements,""));}  //ignore MnPd if requested ternaries
        else{
          if(verbose_elimination){
            message << "Ignoring degenerate (duplicate-element) input (" << original_input << ")";
            pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
          }
        }
      }
    }

    //remove duplicates from vinputs //SILENT
    if(verbose_elimination){
      vector<string> _vinputs=vinputs;vinputs.clear();
      if(_vinputs.size()){vinputs.push_back(_vinputs[0]);}
      for(uint i=1,fl_size_i=_vinputs.size();i<fl_size_i;i++){ //VERBOSE
        if(aurostd::WithinList(vinputs,_vinputs[i])){
          message << "Ignoring duplicate input (" << _vinputs[i] << ")";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
          continue;
        }
        vinputs.push_back(_vinputs[i]);
      }
    }else{std::sort(vinputs.begin(),vinputs.end());vinputs.erase( std::unique( vinputs.begin(), vinputs.end() ), vinputs.end() );}

    message << "Total convex hull inputs: " << vinputs.size();pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);

    //////////////////////////////////////////////////////////////////////////////
    // END Gathering hull inputs
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // START Adding --sc=XX to CHULL::NEGLECT if --output=web
    //////////////////////////////////////////////////////////////////////////////
    //[CO20210315 - using --fake_hull_sc]if(vpflow.flag("CHULL::STABILITY_CRITERION")&&vpflow.flag("CHULL::WEB_DOC")){
    //[CO20210315 - using --fake_hull_sc]  string sc_input=vpflow.getattachedscheme("CHULL::STABILITY_CRITERION");
    //[CO20210315 - using --fake_hull_sc]  if(vpflow.flag("CHULL::NEGLECT")){
    //[CO20210315 - using --fake_hull_sc]    string neglect=vpflow.getattachedscheme("CHULL::NEGLECT");
    //[CO20210315 - using --fake_hull_sc]    if(!neglect.empty()){neglect+=",";}
    //[CO20210315 - using --fake_hull_sc]    neglect+=sc_input;
    //[CO20210315 - using --fake_hull_sc]    vpflow.pop_attached("CHULL::NEGLECT");
    //[CO20210315 - using --fake_hull_sc]    vpflow.push_attached("CHULL::NEGLECT",neglect);
    //[CO20210315 - using --fake_hull_sc]    vpflow.flag("CHULL::NEGLECT",true); //repetita iuvant
    //[CO20210315 - using --fake_hull_sc]  }else{
    //[CO20210315 - using --fake_hull_sc]    vpflow.flag("CHULL::NEGLECT",true);
    //[CO20210315 - using --fake_hull_sc]    vpflow.push_attached("CHULL::NEGLECT",sc_input);
    //[CO20210315 - using --fake_hull_sc]  }
    //[CO20210315 - using --fake_hull_sc]  if(LDEBUG){cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"CHULL::NEGLECT\")=" << vpflow.getattachedscheme("CHULL::NEGLECT") << endl;}
    //[CO20210315 - using --fake_hull_sc]}
    //////////////////////////////////////////////////////////////////////////////
    // END Adding --sc=XX to CHULL::NEGLECT if --output=web
    //////////////////////////////////////////////////////////////////////////////
    
    //////////////////////////////////////////////////////////////////////////////
    // START Looping over hull inputs and creating desired output
    //////////////////////////////////////////////////////////////////////////////
    
    //do NOT pass FileMESSAGE
    //oss is attached to cout or cerr, so it doesn't matter

    uint ntasks=vinputs.size();
    bool Krun=true,KrunSingle=true;
#ifdef AFLOW_MULTITHREADS_ENABLE
    int ncpus=KBIN::get_NCPUS();
    if(ncpus>1){
      xthread::xThread xt(oss,ncpus, 1);
      vector<uint> vKrun(ntasks,1); //c++ doesn't like vector<bool>
      std::function<void(uint,const vector<string>&,const aurostd::xoption&,const _aflags&,vector<uint>&,ostream&)> fn = 
        [&] (uint i,const vector<string>& vinputs,const aurostd::xoption& xoptions,const _aflags& aflags,vector<uint>& vKrun,ostream& oss) {
          vKrun[i]=(convexHull(vinputs[i],xoptions,aflags,oss,true)?1:0);
        };
      //create dumb oss for threaded environment
      ofstream devnull("/dev/null");  //NULL
      ostream oss_dn(devnull.rdbuf());
      xt.run(ntasks,fn,vinputs,vpflow,aflags,vKrun,oss_dn);
      char logger_type=_LOGGER_COMPLETE_;
      for(uint i=0;i<ntasks;i++){
        Krun=(vKrun[i]==1?true:false);
        if(Krun){logger_type=_LOGGER_COMPLETE_;message << vinputs[i] << " completed successfully";}
        else{logger_type=_LOGGER_ERROR_;message << vinputs[i] << " did not complete successfully, run serially";}
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, logger_type);
        Krun=(Krun && KrunSingle);
      }
    }else{
      for(uint i=0;i<ntasks;i++){
        KrunSingle=convexHull(vinputs[i],vpflow,aflags,oss,(bool)i);
        Krun=(Krun && KrunSingle);
      }
    }
#else
    for(uint i=0;i<ntasks;i++){
      KrunSingle=convexHull(vinputs[i],vpflow,aflags,oss,(bool)i);
      Krun=(Krun && KrunSingle);
    }
#endif

    //////////////////////////////////////////////////////////////////////////////
    // END Looping over hull inputs and creating desired output
    //////////////////////////////////////////////////////////////////////////////
    return Krun;
  }

  bool convexHull(const string& input,const aurostd::xoption& _vpflow,const _aflags& aflags,ostream& oss,bool silence_flag_check) {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    stringstream message;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN with input=" << input << endl;}

    string directory=aflags.Directory;
    string log_name="";
    string alloy="";
    aurostd::xoption vpflow(_vpflow);
    ofstream FileMESSAGE;
    
    // go through each request
    // create log specific to that request
    vector<string> velements = aurostd::getElements(input,pp_string,FileMESSAGE,true,true,false,oss); //clean and sort, do not keep_pp  //CO20190712
    //[CO20190712 - OBSOLETE]velements = pflow::getAlphabeticVectorString(input, FileMESSAGE,oss);
    if(!velements.size()){
      message << "Invalid input (" << input << "), please capitalize element symbols";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      return false;/*Krun=false;*continue;*/
    }
    if(velements.size()<2){
      message << "Trivial input (" << input << "), enter binaries or higher";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      return false;/*Krun=false;*continue;*/
    }
    if(XHOST.vflag_control.flag("WWW")&&velements.size()>6){ //CO20200404 - new web flag
      message << velements.size() << "-dimensional hulls cannot be calculated through the web portal (max=6D), please download the AFLOW binary";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      return false;/*Krun=false;*continue;*/
    }
    alloy=aurostd::joinWDelimiter(velements,"");
    if(vpflow.flag("CHULL::LOG")) {
      log_name = "aflow_" + alloy + "_hull.log";
      string log_destination = directory + log_name;  // no output before banner //CO20180220
      FileMESSAGE.open(log_destination.c_str());
    }
    // spit out banner for only the first request
    message << aflow::Banner("BANNER_NORMAL");
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_RAW_, true);  //i //no screen, first to be logged
    message << "Starting " << alloy << " " << pflow::arity_string(velements.size(),false,false) << " convex hull";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);
    getPath(vpflow, FileMESSAGE, oss, false); //CO20180220 - directory stuff for logging
    chull::flagCheck(vpflow, velements, FileMESSAGE, oss, (bool)silence_flag_check);  // spit out all flag options

    ////////////////////////////////////////////////////////////////////////////
    // START Stability criterion calculation
    ////////////////////////////////////////////////////////////////////////////
    //if(vpflow.flag("CHULL::STABILITY_CRITERION")) {
    //  message << "Starting stable criterion calculation of " << vpflow.getattachedscheme("CHULL::STABILITY_CRITERION");
    //  message << " on " << input << " hull";
    //  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);
    //  vector<string> vauid;
    //  vector<double> vscriterion;
    //  aurostd::string2tokens(vpflow.getattachedscheme("CHULL::STABILITY_CRITERION"), vauid, ",");
    //  if(!calculateStabilityCriterion(vpflow,velements,vauid,vscriterion,FileMESSAGE,oss)) {return false;/*Krun=false;*continue;*/} //HAS to be thermal hull by virtue of input
    //  if(vpflow.flag("CHULL::SCREEN_ONLY")){
    //    if(vpflow.flag("CHULL::TEXT_DOC")){
    //      for(uint ia=0,fl_size_ia=vauid.size();ia<fl_size_ia;ia++) {
    //        message << vauid[ia] << ": " << chull::convertUnits(vscriterion[ia], (!vpflow.flag("CHULL::ENTROPIC_TEMPERATURE")?_m_:_std_)) << endl;
    //      }
    //      oss << message.str();
    //    }else if(vpflow.flag("CHULL::JSON_DOC")){
    //      vector<string> vmes;
    //      stringstream dummy;
    //      for(uint ia=0,fl_size_ia=vauid.size();ia<fl_size_ia;ia++) {
    //        dummy << "\"" <<vauid[ia] << "\":" << chull::convertUnits(vscriterion[ia], (!vpflow.flag("CHULL::ENTROPIC_TEMPERATURE")?_m_:_std_));
    //        vmes.push_back(dummy.str()); dummy.str("");
    //      }
    //      oss << aurostd::wrapString(aurostd::joinWDelimiter(vmes,","),"{","}");
    //    }else{ //.log only, but obsolete now anyway since it defaults to json
    //      message << "Unknown print option, only --print=text or --print=json available";
    //      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
    //      if(vpflow.flag("CHULL::LOG")) {FileMESSAGE.close();}
    //      return false;/*Krun=false;*continue;*/
    //    }
    //  }else{
    //    for(uint ia=0,fl_size_ia=vauid.size();ia<fl_size_ia;ia++) {
    //      if(!vpflow.flag("CHULL::ENTROPIC_TEMPERATURE")) {
    //        message << vauid[ia] << " criterion = " << chull::convertUnits(vscriterion[ia], _m_) << " (meV/atom)";
    //        if(std::signbit(vscriterion[ia])) {  //-4e-13 is still negative!
    //          message << ", may NOT be on the hull (negative value)";
    //          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
    //        }else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_COMPLETE_);}
    //      }else{
    //        message << vauid[ia] << " criterion = " << vscriterion[ia] << " (K)";
    //        if(!std::signbit(vscriterion[ia])) {
    //          message << ", may NOT be on the hull (positive value)";
    //          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
    //        }else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_COMPLETE_);}
    //      }
    //    }
    //  }
    //  continue;
    //}
    ////////////////////////////////////////////////////////////////////////////
    // END Stability criterion calculation
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // START Jupyter notebook writing files
    ////////////////////////////////////////////////////////////////////////////
    //MB20190301: For generating Jupyter notebook
    if(vpflow.flag("CHULL::WRITE_JUPYTER2") || vpflow.flag("CHULL::WRITE_JUPYTER3")){

      string aflow_chull_jupyter_json=aurostd::EmbData::get_content("aflow_chull_jupyter.json");
      string aflow_chull_jupyter_subdir = "AFLOW_CHULL_JUPYTER";
      string jupyter_directory=directory + aflow_chull_jupyter_subdir;
      aurostd::DirectoryMake(jupyter_directory);
      aurostd::StringSubst(aflow_chull_jupyter_json,"AFLOW_ALLOY_INSERT_HERE",alloy);
      string ver = (vpflow.flag("CHULL::WRITE_JUPYTER3")) ? "3" : "2";
      string full_ver = (vpflow.flag("CHULL::WRITE_JUPYTER3")) ? "3" : "2";
      aurostd::StringSubst(aflow_chull_jupyter_json,"<ver>",ver);
      aurostd::StringSubst(aflow_chull_jupyter_json,"<full ver>",full_ver);

      // write needed files
      stringstream output;
      output << aflow_chull_jupyter_json;
      aurostd::stringstream2file(output,jupyter_directory+'/'+"notebook.ipynb");
      aurostd::EmbData::save_to_file("aflow_chull_jupyter_plotter.py", jupyter_directory+'/'+"aflow_chull_plotter.py");
      aurostd::EmbData::save_to_file("aflow_chull_python.py", jupyter_directory+'/'+"aflow_chull.py");
      aurostd::EmbData::save_to_file("aflow_chull_jupyter_requirements.txt", jupyter_directory+'/'+"requirements.txt");

      message << "Created " << aflow_chull_jupyter_subdir << " directory for " << input << " hull";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_COMPLETE_);

      return true;
    }     

    ////////////////////////////////////////////////////////////////////////////
    // STOP Jupyter notebook writing files
    ////////////////////////////////////////////////////////////////////////////

    //speed ups for command line options
    if(vpflow.flag("CHULL::DIST2HULL") || vpflow.flag("CHULL::STABILITY_CRITERION") || vpflow.flag("CHULL::N+1_ENTHALPY_GAIN") || vpflow.flag("CHULL::HULL_FORMATION_ENTHALPY")) {
      vpflow.flag("CHULL::SKIP_THERMO_PROPERTIES_EXTRACTION",true);
      message << "Skipping thermodynamic properties extraction";pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);
    }

    //skip calculating all of the lower dimensional hulls
    //this works if you are finding the hull_energy somewhere in the middle of the hull
    //the edges are problematic (lower dimensional hulls)
    //so only skip calculating these hulls if you are calculating hull_energy of the ND hull
    if(0){  //not actually faster, lower dimensional hulls reduce number of points in ND hull, tested 7D-hull: 45 mins vs. 1 hr
      if(vpflow.flag("CHULL::HULL_FORMATION_ENTHALPY")){
        vector<double> _coords;
        aurostd::string2tokens<double>(vpflow.getattachedscheme("CHULL::HULL_FORMATION_ENTHALPY"), _coords, ",");
        if(_coords.size()==velements.size()-1||_coords.size()==velements.size()){
          bool at_edge=false;
          double sum=0.0;
          for(uint ia=0;ia<(velements.size()-1)&&ia<_coords.size();ia++){
            if(abs(_coords[ia])<ZERO_COEF_TOL||abs(1.0-_coords[ia])<ZERO_COEF_TOL){at_edge=true;}
            sum+=_coords[ia];
          }
          double coord_last=1.0-sum;
          at_edge=(at_edge || (abs(coord_last)<ZERO_COEF_TOL||abs(1.0-coord_last)<ZERO_COEF_TOL));
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << " coords=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(_coords,4),",") << endl;
            cerr << __AFLOW_FUNC__ << " sum=" << sum << endl;
            cerr << __AFLOW_FUNC__ << " coord_last=" << coord_last << endl;
            cerr << __AFLOW_FUNC__ << " at_edge=" << at_edge << endl;
          }
          if(at_edge==false){
            vpflow.flag("CHULL::CALCULATE_HIGHEST_DIMENSION_ONLY",true);
            message << "Calculating the highest dimensional hull ONLY";pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);
          }
        }
      }
    }

    ////////////////////////////////////////////////////////////////////////////
    // START Hull initialization
    ////////////////////////////////////////////////////////////////////////////

    ConvexHull hull(vpflow,velements,FileMESSAGE,oss);
    if(!hull.m_initialized) {
      message << "Hull was not created successfully";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      if(vpflow.flag("CHULL::LOG")) {FileMESSAGE.close();}
      if(vpflow.flag("CHULL::SCREEN_ONLY")&&vpflow.flag("CHULL::JSON_DOC")){oss << "{}";}
      return false;/*Krun=false;*continue;*/
    }
    uint dimension = hull.getDim();
    if(!dimension) {
      message << "Hull has no dimensions";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      if(vpflow.flag("CHULL::LOG")) {FileMESSAGE.close();}
      return false;/*Krun=false;*continue;*/
    }
    if(dimension < 2) {
      message << "Unable to calculate hulls with dimensions less than 2";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      if(vpflow.flag("CHULL::LOG")) {FileMESSAGE.close();}
      return false;/*Krun=false;*continue;*/
    }
    if(dimension != velements.size()) {
      message << "Dimension of hull does not reflect the number of elements";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      if(vpflow.flag("CHULL::LOG")) {FileMESSAGE.close();}
      return false;/*Krun=false;*continue;*/
    }

    ////////////////////////////////////////////////////////////////////////////
    // END Hull initialization
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // START Distance to hull calculation
    ////////////////////////////////////////////////////////////////////////////
    if(vpflow.flag("CHULL::DIST2HULL")) {
      message << "Starting distance to hull calculation of " << vpflow.getattachedscheme("CHULL::DIST2HULL");
      message << " on " << input << " hull";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);
      vector<string> vauid;
      aurostd::string2tokens(vpflow.getattachedscheme("CHULL::DIST2HULL"), vauid, ",");
      vector<double> vdist2hull;
      //NB: to anyone who is using the convex hull object
      //outside of declaration/initialization, all functions should be wrapped
      //in try/catch's to avoid hard exits
      //proceed otherwise at your own risk
      try{vdist2hull=hull.getDistancesToHull(vauid);}
      catch(aurostd::xerror& err){
        pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
        if(vpflow.flag("CHULL::LOG")) {FileMESSAGE.close();}
        if(vpflow.flag("CHULL::SCREEN_ONLY")&&vpflow.flag("CHULL::JSON_DOC")){oss << "{}";} //so JSON-reader doesn't bomb
        return true;
      }

      //set correct sign convention
      //since this SHOULD be a positive distance (inside hull), we will negate for formation_energy_hull
      //for (uint ia = 0; ia < vauid.size(); ia++) {
      //  if(!vpflow.flag("CHULL::ENTROPIC_TEMPERATURE")) {vdist2hull[ia]=abs(vdist2hull[ia]);}
      //}

      for(uint ia=0,fl_size_ia=vauid.size();ia<fl_size_ia;ia++) {
        if(!vpflow.flag("CHULL::ENTROPIC_TEMPERATURE")) {message << vauid[ia] << " dist2hull = " << chull::convertUnits(vdist2hull[ia], _m_) << " (meV/atom)";}
        else{message << vauid[ia] << " dist2hull = " << vdist2hull[ia] << " (K)";}
        if(aurostd::zeroWithinTol(vdist2hull[ia],ZERO_TOL)) {  //do not issue a warning either way, it's simply the distance
          message << ", may be on the hull";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_COMPLETE_);
        }else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_COMPLETE_);}
        //if(0&&!std::signbit(vdist2hull[ia])) {  //do not issue a warning either way, it's simply the distance
        //  message << ", may NOT be off the hull (positive value)";
        //  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
        //}else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_COMPLETE_);}
      }

      if(vpflow.flag("CHULL::SCREEN_ONLY")){
        if(vpflow.flag("CHULL::TEXT_DOC")){
          for(uint ia=0,fl_size_ia=vauid.size();ia<fl_size_ia;ia++) {
            message << vauid[ia] << ": " << chull::convertUnits(vdist2hull[ia], (!vpflow.flag("CHULL::ENTROPIC_TEMPERATURE")?_m_:_std_)) << endl;
          }
          oss << message.str();
        }else if(vpflow.flag("CHULL::JSON_DOC")){
          vector<string> vmes;
          stringstream dummy;
          for(uint ia=0,fl_size_ia=vauid.size();ia<fl_size_ia;ia++) {
            dummy << "\"" <<vauid[ia] << "\":" << chull::convertUnits(vdist2hull[ia], (!vpflow.flag("CHULL::ENTROPIC_TEMPERATURE")?_m_:_std_));
            vmes.push_back(dummy.str()); dummy.str("");
          }
          oss << aurostd::wrapString(aurostd::joinWDelimiter(vmes,","),"{","}");
        }else{ //.log only, but obsolete now anyway since it defaults to json
          message << "Unknown print option, only --print=text or --print=json available";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
          if(vpflow.flag("CHULL::LOG")) {FileMESSAGE.close();}
          return false;/*Krun=false;*continue;*/
        }
      }
      return true;
    }
    ////////////////////////////////////////////////////////////////////////////
    // END Distance to hull calculation
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // START Stability criterion calculation
    ////////////////////////////////////////////////////////////////////////////
    if(vpflow.flag("CHULL::STABILITY_CRITERION")) { //CO20210201 - chull-web SS plotter //&&(!vpflow.flag("CHULL::WEB_DOC"))
      message << "Starting stable criterion calculation of " << vpflow.getattachedscheme("CHULL::STABILITY_CRITERION");
      message << " on " << input << " hull";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);
      vector<string> vauid;
      aurostd::string2tokens(vpflow.getattachedscheme("CHULL::STABILITY_CRITERION"), vauid, ",");
      vector<double> vscriterion;
      //NB: to anyone who is using the convex hull object
      //outside of declaration/initialization, all functions should be wrapped
      //in try/catch's to avoid hard exits
      //proceed otherwise at your own risk
      try{vscriterion=hull.getStabilityCriterion(vauid);}
      catch(aurostd::xerror& err){
        pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
        if(vpflow.flag("CHULL::LOG")) {FileMESSAGE.close();}
        if(vpflow.flag("CHULL::SCREEN_ONLY")&&vpflow.flag("CHULL::JSON_DOC")){oss << "{}";} //so JSON-reader doesn't bomb
        return true;
      }

      ////set correct sign convention
      ////since this SHOULD be a negative distance (outside hull), we will negate for formation_energy_hull
      //for (uint ia = 0; ia < vauid.size(); ia++) {
      //  if(!vpflow.flag("CHULL::ENTROPIC_TEMPERATURE")) {vscriterion[ia]=-vscriterion[ia];}
      //}

      for(uint ia=0,fl_size_ia=vauid.size();ia<fl_size_ia;ia++) {
        if(!vpflow.flag("CHULL::ENTROPIC_TEMPERATURE")) {message << vauid[ia] << " criterion = " << chull::convertUnits(vscriterion[ia], _m_) << " (meV/atom)";}
        else{message << vauid[ia] << " criterion = " << vscriterion[ia] << " (K)";}
        if(std::signbit(vscriterion[ia])) {
          message << ", may NOT be on the hull (negative value)";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
        }else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_COMPLETE_);}
      }

      if(vpflow.flag("CHULL::SCREEN_ONLY")){
        if(vpflow.flag("CHULL::TEXT_DOC")){
          for(uint ia=0,fl_size_ia=vauid.size();ia<fl_size_ia;ia++) {
            message << vauid[ia] << ": " << chull::convertUnits(vscriterion[ia], (!vpflow.flag("CHULL::ENTROPIC_TEMPERATURE")?_m_:_std_)) << endl;
          }
          oss << message.str();
        }else if(vpflow.flag("CHULL::JSON_DOC")){
          vector<string> vmes;
          stringstream dummy;
          for(uint ia=0,fl_size_ia=vauid.size();ia<fl_size_ia;ia++) {
            dummy << "\"" <<vauid[ia] << "\":" << chull::convertUnits(vscriterion[ia], (!vpflow.flag("CHULL::ENTROPIC_TEMPERATURE")?_m_:_std_));
            vmes.push_back(dummy.str()); dummy.str("");
          }
          oss << aurostd::wrapString(aurostd::joinWDelimiter(vmes,","),"{","}");
        }else{ //.log only, but obsolete now anyway since it defaults to json
          message << "Unknown print option, only --print=text or --print=json available";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
          if(vpflow.flag("CHULL::LOG")) {FileMESSAGE.close();}
          return false;/*Krun=false;*continue;*/
        }
      }
      return true;
    }
    ////////////////////////////////////////////////////////////////////////////
    // END Stability criterion calculation
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // START N+1 enthalpy gain calculation
    ////////////////////////////////////////////////////////////////////////////
    if(vpflow.flag("CHULL::N+1_ENTHALPY_GAIN")) {
      message << "Starting N+1 enthalpy gain calculation of " << vpflow.getattachedscheme("CHULL::N+1_ENTHALPY_GAIN");
      message << " on " << input << " hull";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);
      vector<string> vauid;
      aurostd::string2tokens(vpflow.getattachedscheme("CHULL::N+1_ENTHALPY_GAIN"), vauid, ",");
      vector<double> vn1egain;
      //NB: to anyone who is using the convex hull object
      //outside of declaration/initialization, all functions should be wrapped
      //in try/catch's to avoid hard exits
      //proceed otherwise at your own risk
      try{vn1egain=hull.getNPlus1EnthalpyGain(vauid);}
      catch(aurostd::xerror& err){
        pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
        if(vpflow.flag("CHULL::LOG")) {FileMESSAGE.close();}
        if(vpflow.flag("CHULL::SCREEN_ONLY")&&vpflow.flag("CHULL::JSON_DOC")){oss << "{}";} //so JSON-reader doesn't bomb
        return true;
      }

      ////set correct sign convention
      ////since this SHOULD be a negative distance (outside hull), we will negate for formation_energy_hull
      //for (uint ia = 0; ia < vauid.size(); ia++) {
      //  if(!vpflow.flag("CHULL::ENTROPIC_TEMPERATURE")) {vn1egain[ia]=-vn1egain[ia];}
      //}

      for(uint ia=0,fl_size_ia=vauid.size();ia<fl_size_ia;ia++) {
        if(!vpflow.flag("CHULL::ENTROPIC_TEMPERATURE")) {message << vauid[ia] << " n+1_enthalpy_gain = " << chull::convertUnits(vn1egain[ia], _m_) << " (meV/atom)";}
        else{message << vauid[ia] << " n+1_enthalpy_gain = " << vn1egain[ia] << " (K)";}
        if(std::signbit(vn1egain[ia])) {
          message << ", may NOT be on the hull (negative value)";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
        }else{pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_COMPLETE_);}
      }

      if(vpflow.flag("CHULL::SCREEN_ONLY")){
        if(vpflow.flag("CHULL::TEXT_DOC")){
          for(uint ia=0,fl_size_ia=vauid.size();ia<fl_size_ia;ia++) {
            message << vauid[ia] << ": " << chull::convertUnits(vn1egain[ia], (!vpflow.flag("CHULL::ENTROPIC_TEMPERATURE")?_m_:_std_)) << endl;
          }
          oss << message.str();
        }else if(vpflow.flag("CHULL::JSON_DOC")){
          vector<string> vmes;
          stringstream dummy;
          for(uint ia=0,fl_size_ia=vauid.size();ia<fl_size_ia;ia++) {
            dummy << "\"" <<vauid[ia] << "\":" << chull::convertUnits(vn1egain[ia], (!vpflow.flag("CHULL::ENTROPIC_TEMPERATURE")?_m_:_std_));
            vmes.push_back(dummy.str()); dummy.str("");
          }
          oss << aurostd::wrapString(aurostd::joinWDelimiter(vmes,","),"{","}");
        }else{ //.log only, but obsolete now anyway since it defaults to json
          message << "Unknown print option, only --print=text or --print=json available";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
          if(vpflow.flag("CHULL::LOG")) {FileMESSAGE.close();}
          return false;/*Krun=false;*continue;*/
        }
      }
      return true;
    }
    ////////////////////////////////////////////////////////////////////////////
    // END N+1 enthalpy gain calculation
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // START Hull formation enthalpy calculation
    ////////////////////////////////////////////////////////////////////////////
    if(vpflow.flag("CHULL::HULL_FORMATION_ENTHALPY")) {
      message << "Starting calculation of the formation enthalpy at " << vpflow.getattachedscheme("CHULL::HULL_FORMATION_ENTHALPY");
      message << " on " << input << " hull";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);
      vector<double> _coords;
      xvector<double> coords(dimension);
      aurostd::string2tokens<double>(vpflow.getattachedscheme("CHULL::HULL_FORMATION_ENTHALPY"), _coords, ",");
      for(uint j=0,fl_size_j=_coords.size();j<fl_size_j&&j<dimension;j++){coords[j+coords.lrows]=_coords[j];}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " coords=" << coords << endl;}
      double dist2hull=0.0;
      //NB: to anyone who is using the convex hull object
      //outside of declaration/initialization, all functions should be wrapped
      //in try/catch's to avoid hard exits
      //proceed otherwise at your own risk
      try{
        ChullPoint cp(coords,FileMESSAGE,oss,hull.m_has_stoich_coords,hull.m_formation_energy_hull,false);  //not a real point
        dist2hull=hull.getDistanceToHull(cp,false,true);  //do not redo, get signed distance (this is energy)
        bool should_be_positive=!hull.m_lower_hull;
        bool correct_sign_vertical_distance=chull::correctSignVerticalDistance(dist2hull,should_be_positive);
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " dist2hull=" << dist2hull << endl;
          cerr << __AFLOW_FUNC__ << " correct_sign_vertical_distance=" << correct_sign_vertical_distance << endl;
        }
        if(!correct_sign_vertical_distance){dist2hull*=-1.0;}
      }
      catch(aurostd::xerror& err){
        pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
        if(vpflow.flag("CHULL::LOG")) {FileMESSAGE.close();}
        if(vpflow.flag("CHULL::SCREEN_ONLY")&&vpflow.flag("CHULL::JSON_DOC")){oss << "{}";} //so JSON-reader doesn't bomb
        return true;
      }

      if(vpflow.flag("CHULL::SCREEN_ONLY")){
        uint precision=COEF_PRECISION;
        double roundoff_tol=5.0*pow(10,-((int)precision)-1);
        stringstream hull_energy_ss;
        hull_energy_ss << "\"hull_energy";
        hull_energy_ss << aurostd::wrapString(aurostd::joinWDelimiter(xvecDouble2vecString(coords,precision,true,roundoff_tol,FIXED_STREAM),","),"[","]");
        hull_energy_ss << "\": ";
        hull_energy_ss << chull::convertUnits(dist2hull, (!vpflow.flag("CHULL::ENTROPIC_TEMPERATURE")?_m_:_std_));
        oss << aurostd::wrapString(hull_energy_ss.str(),"{","}");
      }else{
        //[CO20190801 - OBSOLETE]message << "hull_energy[coords=" << aurostd::joinWDelimiter(coords,",");
        //[CO20190801 - OBSOLETE]for(int j=coords.lrows;j<=coords.urows;j++){
        //[CO20190801 - OBSOLETE]  message << coords[j];
        //[CO20190801 - OBSOLETE]  if(j!=((int)dimension)-1){message << ",";}
        //[CO20190801 - OBSOLETE]}
        //[CO20190801 - OBSOLETE]message << "] = ";
        uint precision=COEF_PRECISION;
        double roundoff_tol=5.0*pow(10,-((int)precision)-1);
        message << "hull_energy" << aurostd::wrapString("coords="+aurostd::joinWDelimiter(xvecDouble2vecString(coords,precision,true,roundoff_tol,FIXED_STREAM),","),"[","]");
        message << " = ";
        if(!vpflow.flag("CHULL::ENTROPIC_TEMPERATURE")) {message << chull::convertUnits(dist2hull, _m_) << " (meV/atom)";}
        else{message << dist2hull << " (K)";}
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_COMPLETE_);
      }
      return true;
    }
    ////////////////////////////////////////////////////////////////////////////
    // END Hull formation enthalpy calculation
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // START outputs
    ////////////////////////////////////////////////////////////////////////////

    if(vpflow.flag("CHULL::TEXT_DOC")) {if(!hull.write(txt_ft)) {if(vpflow.flag("CHULL::LOG")) {FileMESSAGE.close();} return false;/*Krun=false;*continue;*/}} // text doc
    if(vpflow.flag("CHULL::JSON_DOC")) {if(!hull.write(json_ft)) {if(vpflow.flag("CHULL::LOG")) {FileMESSAGE.close();} return false;/*Krun=false;*continue;*/}} // json doc
    if(vpflow.flag("CHULL::WEB_DOC")) {if(!hull.write(chull_web_ft)) {if(vpflow.flag("CHULL::LOG")) {FileMESSAGE.close();} return false;/*Krun=false;*continue;*/}} // web-specific json doc
    if(vpflow.flag("CHULL::LATEX_DOC")||vpflow.flag("CHULL::PNG_IMAGE")) {if(!hull.write(latex_ft)) {if(vpflow.flag("CHULL::LOG")) {FileMESSAGE.close();} return false;/*Krun=false;*continue;*/}} // latex doc, keep last as it will trip on useless plots

    ////////////////////////////////////////////////////////////////////////////
    // END outputs
    ////////////////////////////////////////////////////////////////////////////

    // close input specific log
    if(vpflow.flag("CHULL::LOG")) {FileMESSAGE.close();}
    return true;
  }
} // namespace chull

namespace chull {
  //***************************************************************************//
  // chull::getPath(aurostd::xoption& vpflow,const bool& silent,ostringstream&
  // oss,ofstream& FileMESSAGE)
  //***************************************************************************//
  // gets path to redirect output
  string getPath(bool add_backslash) {return aurostd::getPWD() + (add_backslash?string("/"):string(""));} //[CO20191112 - OBSOLETE]aurostd::execute2string(XHOST.command("pwd"))
  string getPath(const aurostd::xoption& vpflow, ostream& oss, bool silent) {  // overload
    ofstream FileMESSAGE;
    return getPath(vpflow, FileMESSAGE, oss, silent);
  }
  string getPath(const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& oss, bool silent) {  // main function
    stringstream message;
    if(!vpflow.flag("CHULL::PATH")) {
      string pwd = getPath();
      if(!silent){
        message << "Directing output to current directory: " << pwd;
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_OPTION_); //, silent);  //CO20180220 - silent now means print AT ALL
      }
      return pwd;
    }
    return getPath(vpflow.getattachedscheme("CHULL::PATH"), FileMESSAGE, oss, silent);
  }
  //***************************************************************************//
  // chull::getPath(string _path,const bool& silent,ostream& oss,ofstream&
  // FileMESSAGE)
  //***************************************************************************//
  // gets path to redirect output
  string getPath(string _path, ostream& oss, bool silent) {  // overload
    ofstream FileMESSAGE;
    return getPath(_path, FileMESSAGE, oss, silent);
  }
  string getPath(string _path, ofstream& FileMESSAGE, ostream& oss, bool silent) {  // main function
    stringstream message;
    string pwd = getPath(false);
    string home = XHOST.home; //aurostd::execute2string(XHOST.command("echo") + " $HOME");
    string path="";
    // add '/' if _path doesn't already have it
    if(_path[_path.length() - 1] != '/') {
      _path += "/";
    }
    // remove ./ for relative paths
    if(_path[0] == '.') {
      _path = _path.substr(1, _path.length());
      if(!_path.empty() && _path[0] == '/') {
        _path = _path.substr(1, _path.length());
      }
    }
    // home doesn't have last '/'
    if(!_path.empty() && _path[0] == '~') {
      _path = _path.substr(1, _path.length());
      if(!_path.empty() && _path[0] == '/') {
        _path = _path.substr(1, _path.length());
      }
      _path = home + "/" + _path;
    }
    // if not root path (starting with '/'), it's a relative path, add it to pwd,
    // pwd doesn't have last '/'
    if(_path.empty() || (!_path.empty() && _path[0] != '/')) {
      path = pwd + "/";
      path += _path;
    }else{
      path = _path;
    }

    //test of stupidity
    aurostd::DirectoryMake(path);
    if(!aurostd::IsDirectory(path)){
      message << path << " does not seem to be a viable directory, changing to pwd=" << pwd;
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_WARNING_);
      path=pwd+"/";
    }

    if(!silent){
      message << "Directing output to " << path;
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_OPTION_); //, silent);  //CO20180220 - silent now means print AT AL
    }
    return path;
  }
} // namespace chull

namespace chull {
  //***************************************************************************//
  // chull::flagCheck(aurostd::xoption& vpflow,bool silent,ostringstream& oss,ofstream& FileMESSAGE)
  //***************************************************************************//
  // logs which flags are on
  void flagCheck(aurostd::xoption& vpflow, const vector<string>& velements, ostream& oss, bool silent) {  // overload
    ofstream FileMESSAGE;
    flagCheck(vpflow, velements, FileMESSAGE, oss, silent);
  }
  void flagCheck(aurostd::xoption& vpflow, const vector<string>& velements, ofstream& FileMESSAGE, ostream& oss, bool silent) {  // main function
    stringstream message;
    string directory=getPath(vpflow,FileMESSAGE,oss);
    _aflags aflags; aflags.Directory=directory;
    if(vpflow.flag("CHULL::TEXT_DOC")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::TEXT_DOC set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::JSON_DOC")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::JSON_DOC set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::WEB_DOC")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::WEB_DOC set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::LATEX_DOC")||vpflow.flag("CHULL::PNG_IMAGE")) {
      if(vpflow.flag("CHULL::LATEX_DOC")){
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::LATEX_DOC set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
      }
      if(vpflow.flag("CHULL::PNG_IMAGE")){
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::PNG_IMAGE set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
      }
      if(vpflow.flag("CHULL::IMAGE_ONLY")) {
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::IMAGE_ONLY set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
      }
      if(vpflow.flag("CHULL::NO_DOC")) {
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::NO_DOC set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
      }
      if(vpflow.flag("CHULL::DOC_ONLY")) {
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::DOC_ONLY set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
      }
      if(vpflow.flag("CHULL::KEEP_TEX")) {
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::KEEP_TEX set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
      }
      if(vpflow.flag("CHULL::LATEX_OUTPUT")) {
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::LATEX_OUTPUT set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
      }
      if(vpflow.flag("CHULL::LATEX_INTERACTIVE")) {
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::LATEX_INTERACTIVE set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
      }
      if(vpflow.flag("CHULL::LIGHT_CONTRAST")) {
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::LIGHT_CONTRAST set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
      }
      if(vpflow.flag("CHULL::LARGE_FONT")) {
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::LARGE_FONT set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
      }
      if(vpflow.flag("CHULL::PLOT_ISO_MAX_LATENT_HEAT")) {
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::PLOT_ISO_MAX_LATENT_HEAT set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
      }
    }
    if(vpflow.flag("CHULL::SCREEN_ONLY")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::SCREEN_ONLY set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::LOG")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::LOG set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    string lib_count_string,load_lib_flag_name;
    for(uint lib=1,fl_size_lib=velements.size();lib<=fl_size_lib && lib<=_AFLOW_LIB_MAX_;lib++) {
      lib_count_string=aurostd::utype2string(lib);
      load_lib_flag_name="PFLOW::LOAD_ENTRIES_LOAD_LIB"+lib_count_string;
      if(vpflow.flag(load_lib_flag_name)) {
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, load_lib_flag_name+" set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
      }
    }
    if(vpflow.flag("PFLOW::LOAD_ENTRIES_NARIES_MINUS_ONE")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "PFLOW::LOAD_ENTRIES_NARIES_MINUS_ONE set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("PFLOW::LOAD_ENTRIES_LOAD_ICSD")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "PFLOW::LOAD_ENTRIES_LOAD_ICSD set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("PFLOW::LOAD_API")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "PFLOW::LOAD_API set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("PFLOW::LOAD_ENTRIES_ENTRY_OUTPUT")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "PFLOW::LOAD_ENTRIES_ENTRY_OUTPUT set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::NEGLECT")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::NEGLECT set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::SEE_NEGLECT")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::SEE_NEGLECT set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::CALCULATE_FAKE_HULL_STABILITY_CRITERION")) { //CO20210315
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::CALCULATE_FAKE_HULL_STABILITY_CRITERION set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::CALCULATE_FAKE_HULL_N+1_ENTHALPY_GAIN")) { //SK20200327
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::CALCULATE_FAKE_HULL_N+1_ENTHALPY_GAIN set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::REMOVE_EXTREMA")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::REMOVE_EXTREMA set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::ENTROPIC_TEMPERATURE")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::ENTROPIC_TEMPERATURE set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::INCLUDE_PAW_GGA")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::INCLUDE_PAW_GGA set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::DIST2HULL")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::DIST2HULL set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::STABILITY_CRITERION")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::STABILITY_CRITERION set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::HULL_FORMATION_ENTHALPY")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::HULL_FORMATION_ENTHALPY set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::DIST_FOR_EQUIVALENCE_ANALYSIS")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::DIST_FOR_EQUIVALENCE_ANALYSIS set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::DIST_FOR_EQUIVALENCE_ANALYSIS_ND_ONLY")){
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::DIST_FOR_EQUIVALENCE_ANALYSIS_ND_ONLY set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::SKIP_STRUCTURE_COMPARISON")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::SKIP_STRUCTURE_COMPARISON set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::SKIP_STABILITY_CRITERION_ANALYSIS")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::SKIP_STABILITY_CRITERION_ANALYSIS set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::SKIP_N+1_ENTHALPY_GAIN_ANALYSIS")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::SKIP_N+1_ENTHALPY_GAIN_ANALYSIS set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::INCLUDE_SKEWED_HULLS")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::INCLUDE_SKEWED_HULLS set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::INCLUDE_UNRELIABLE_HULLS")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::INCLUDE_UNRELIABLE_HULLS set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::INCLUDE_OUTLIERS")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::INCLUDE_OUTLIERS set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::STRICT_OUTLIER_ANALYSIS")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::STRICT_OUTLIER_ANALYSIS set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::INCLUDE_ILL_CONVERGED")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::INCLUDE_ILL_CONVERGED set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("CHULL::CALCULATE_HIGHEST_DIMENSION_ONLY")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::CALCULATE_HIGHEST_DIMENSION_ONLY set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if(vpflow.flag("FORCE")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "CHULL::FORCE set to TRUE", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
  }
} // namespace chull

//namespace chull {
//  bool calculateStabilityCriterion(const vector<string>& velements,const string& auid,double& scriterion,ostream& oss){
//    ofstream FileMESSAGE;
//    return calculateStabilityCriterion(velements,auid,scriterion,FileMESSAGE,oss);
//  }
//  bool calculateStabilityCriterion(const vector<string>& velements,const string& auid,double& scriterion,ofstream& FileMESSAGE,ostream& oss){
//    vector<string> vauid;
//    vauid.push_back(auid);
//    vector<double> vscriterion;
//    if(!calculateStabilityCriterion(velements,vauid,vscriterion,FileMESSAGE,oss)) {return FALSE;}
//    scriterion = vscriterion[0];
//    return true;
//  }
//  bool calculateStabilityCriterion(const vector<string>& velements,const vector<string>& vauid,vector<double>& vscriterion,ostream& oss){
//    ofstream FileMESSAGE;
//    return calculateStabilityCriterion(velements,vauid,vscriterion,FileMESSAGE,oss);
//  }
//  bool calculateStabilityCriterion(const vector<string>& velements,const vector<string>& vauid,vector<double>& vscriterion,ofstream& FileMESSAGE,ostream& oss){
//    aurostd::xoption vpflow;
//    pflow::defaultLoadEntriesFlags(vpflow,FileMESSAGE,oss,std::string("A"),false,true);
//    return calculateStabilityCriterion(vpflow,velements,vauid,vscriterion,FileMESSAGE,oss);
//  }
//  bool calculateStabilityCriterion(const aurostd::xoption& vpflow,const vector<string>& velements,double& scriterion,ostream& oss){
//    ofstream FileMESSAGE;
//    return calculateStabilityCriterion(vpflow,velements,scriterion,FileMESSAGE,oss);
//  }
//  bool calculateStabilityCriterion(const aurostd::xoption& vpflow,const vector<string>& velements,double& scriterion,ofstream& FileMESSAGE,ostream& oss){
//    stringstream message;
//    if(!vpflow.flag("CHULL::NEGLECT")) {
//      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,"CHULL::NEGLECT not set",FileMESSAGE,oss,_LOGGER_ERROR_);
//      return false;
//    }
//    vector<string> points_neglect;
//    aurostd::string2tokens(vpflow.getattachedscheme("CHULL::NEGLECT"),points_neglect,",");
//    if(points_neglect.size()!=1) {
//      message << "Can only handle one AUID at a time, " << points_neglect.size() << " given";
//      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,FileMESSAGE,oss,_LOGGER_ERROR_);
//      return FALSE;
//    }
//    return calculateStabilityCriterion(vpflow,velements,points_neglect[0],scriterion,FileMESSAGE,oss);
//  }
//  bool calculateStabilityCriterion(const aurostd::xoption& vpflow,const vector<string>& velements,const string& auid,double& scriterion,ostream& oss){
//    ofstream FileMESSAGE;
//    return calculateStabilityCriterion(vpflow,velements,auid,scriterion,FileMESSAGE,oss);
//  }
//  bool calculateStabilityCriterion(const aurostd::xoption& vpflow,const vector<string>& velements,const string& auid,double& scriterion,ofstream& FileMESSAGE,ostream& oss){
//    vector<string> vauid;
//    vauid.push_back(auid);
//    vector<double> vscriterion;
//    if(!calculateStabilityCriterion(vpflow,velements,vauid,vscriterion,FileMESSAGE,oss)) {return FALSE;}
//    scriterion = vscriterion[0];
//    return true;
//  }
//  bool calculateStabilityCriterion(const aurostd::xoption& vpflow,const vector<string>& velements,vector<double>& vscriterion,ostream& oss){
//    ofstream FileMESSAGE;
//    return calculateStabilityCriterion(vpflow,velements,vscriterion,FileMESSAGE,oss);
//  }
//  bool calculateStabilityCriterion(const aurostd::xoption& vpflow,const vector<string>& velements,vector<double>& vscriterion,ofstream& FileMESSAGE,ostream& oss){
//    stringstream message;
//    if(!vpflow.flag("CHULL::NEGLECT")) {
//      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,"CHULL::NEGLECT not set",FileMESSAGE,oss,_LOGGER_ERROR_);
//      return FALSE;
//    }
//    vector<string> points_neglect;
//    aurostd::string2tokens(vpflow.getattachedscheme("CHULL::NEGLECT"),points_neglect,",");
//    return calculateStabilityCriterion(vpflow,velements,points_neglect,vscriterion,FileMESSAGE,oss);
//  }
//  bool calculateStabilityCriterion(const aurostd::xoption& vpflow,const vector<string>& velements,const vector<string>& vauid,vector<double>& vscriterion,ostream& oss) {
//    ofstream FileMESSAGE;
//    return calculateStabilityCriterion(vpflow,velements,vauid,vscriterion,FileMESSAGE,oss);
//  }
//  bool calculateStabilityCriterion(const aurostd::xoption& vpflow,const vector<string>& velements,const vector<string>& vauid,vector<double>& vscriterion,ofstream& FileMESSAGE,ostream& oss) {
//    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
//    stringstream message;
//    
//    //don't only remove vauid, but also equivalent gstates
//    //we need organized points, simply initialize dummy hull instead of calculating TWO full hulls
//    ConvexHull dummy(vpflow,FileMESSAGE,oss);
//    try{dummy.initializePoints(velements);}
//    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), FileMESSAGE, oss, _LOGGER_ERROR_);return false;}
//    dummy.m_initialized=true; //hack so we can get at the g-states
//    uint i_point,i_coord_group,g_state;
//    vector<uint> eq_gstates;
//    bool found=false;
//    vector<ChullPoint> points_to_neglect; //follows vauid
//    for(uint i=0,fl_size_i=vauid.size();i<fl_size_i;i++){
//      const string& auid=vauid[i];
//      if(auid.empty()){
//        message << "Empty auid found";
//        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,FileMESSAGE,oss,_LOGGER_ERROR_);
//        return false;
//      }
//      if(!dummy.findPoint(auid,i_point)){
//        message << "Specified auid not found on hull (auid=" << auid << ")";
//        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,FileMESSAGE,oss,_LOGGER_ERROR_);
//        return false;
//      }
//      const ChullPoint& point=dummy.m_points[i_point];  //this point may not be on the hull, it may be an equivalent structure, but coordgroup is
//      if(!dummy.getCoordGroupIndex(point,i_coord_group)){
//        message << "Coordgroup index not set (auid=" << auid << ")";
//        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,FileMESSAGE,oss,_LOGGER_ERROR_);
//        return false;
//      }
//      if(!dummy.m_coord_groups[i_coord_group].m_points.size()){
//        message << "No points found within coordgroup (auid=" << auid << ")";
//        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,FileMESSAGE,oss,_LOGGER_ERROR_);
//        return false;
//      }
//      //assume user knows what he's doing, we will check for sure later
//      dummy.m_coord_groups[i_coord_group].m_is_on_hull=true;
//      g_state=dummy.m_coord_groups[i_coord_group].m_ref_state;
//      dummy.m_points[g_state].m_is_g_state=true;
//      eq_gstates=dummy.getEquivalentEntries(g_state);
//      if(!eq_gstates.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No equivalent states found (not even self)");}
//      //check to make sure entry with specified auid among gstates
//      found=false;
//      for(uint j=0,fl_size_j=eq_gstates.size();j<fl_size_j;j++){
//        if(dummy.m_points[eq_gstates[j]].m_entry.auid==auid){found=true;}
//        points_to_neglect.push_back(dummy.m_points[eq_gstates[j]]);
//      }
//      if(!found){
//        message << "Point was not found to be an equivalent ground-state structure (auid=" << auid << ")";
//        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,FileMESSAGE,oss,_LOGGER_ERROR_);
//        return false;
//      }
//    }
//    
//    //////////////////////////////////////////////////////////////////////////////
//    // START Getting entries to be neglected
//    //////////////////////////////////////////////////////////////////////////////
//  
//    vector<ChullPoint> new_points;
//    const vector<ChullPoint>& points=dummy.m_points;
//    for(uint i=0,fl_size_i=points.size();i<fl_size_i;i++){
//      if(points[i].m_is_artificial){continue;}  //they will be added again
//      found=false;
//      for(uint j=0,fl_size_j=points_to_neglect.size();j<fl_size_j&&!found;j++){
//        if(points[i].m_entry.auid==points_to_neglect[j].m_entry.auid){
//          message << "Removing equivalent ground-state (auid=" << points_to_neglect[j].m_entry.auid << ")";
//          pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,FileMESSAGE,oss,_LOGGER_MESSAGE_);
//          found=true;
//        }
//      }
//      if(found){continue;}
//      new_points.push_back(points[i]);
//    }
//    
//    //////////////////////////////////////////////////////////////////////////////
//    // END Getting entries to be neglected
//    //////////////////////////////////////////////////////////////////////////////
//    
//    message << "Creating new hull without relevant g-states";
//    pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,FileMESSAGE,oss,_LOGGER_MESSAGE_);
//    aurostd::xoption cflags=vpflow;
//    cflags.flag("CHULL::SKIP_THERMO_PROPERTIES_EXTRACTION",true); //thermo properties NOT needed, just need hull
//    cflags.flag("CHULL::SKIP_STABILITY_CRITERION_ANALYSIS",true); //thermo properties NOT needed, just need hull
//    ConvexHull hull(cflags,new_points,velements,FileMESSAGE,oss,true,true);
//    if(!hull.m_initialized){return false;}
//    //since i_nary and i_alloy don't change, the getDistanceToHull function should work fine (getRelevantFacets())
//    uint i_nary,i_alloy;
//    for(uint i=0,fl_size_i=points_to_neglect.size();i<fl_size_i;i++) {
//      ChullPoint& point=points_to_neglect[i];
//      if(!dummy.getAlloyIndex(point,i_nary,i_alloy)){
//        message << "Alloy index not set (auid=" << point.m_entry.auid << ")";
//        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,FileMESSAGE,oss,_LOGGER_ERROR_);
//        return false;
//      }
//      if(point.isUnary()){point.setHullCoords();} //set to most general coords (m_coords), this reflects relevantFacets()
//      else{
//        xvector<int> elements_present=dummy.m_naries[i_nary].m_alloys[i_alloy].m_elements_present;
//        if(LDEBUG){cerr << __AFLOW_FUNC__ << " elements_present=" << elements_present << endl;}
//        point.setHullCoords(elements_present);  //just to be sure
//      }
//      try{vscriterion.push_back(hull.getDistanceToHull(point));}
//      catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), FileMESSAGE, oss, _LOGGER_ERROR_);return false;}
//    }
//
//    return true;
//  }
//} // namespace chull

namespace chull {
  // ***************************************************************************
  // chull::convertUnits(double value,char units)
  // ***************************************************************************
  // returns value in desired units
  double convertUnits(double value, char units) {
    if(units == _m_) {return value * 1e3;}
    return value;
  }
  //these functions are auxiliary - if entry could be passed instead of point, then call these functions
  //otherwise, keep as methods of the ChullPoint class
  double H_f_atom(const ChullPoint& point, char units) {
    if(point.m_has_entry){return H_f_atom(point.m_entry,units);}
    if(point.m_formation_energy_coord){return convertUnits(point.getLastCoord(),units);}
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No formation energy available for ChullPoint",_INPUT_ILLEGAL_);
    return AUROSTD_NAN;
  }
  double H_f_atom(const aflowlib::_aflowlib_entry& entry, char units){
    double d_tmp=entry.enthalpyFormationAtom(0);  //entry.enthalpy_formation_atom - ADDING CCE @ 0K
    if(d_tmp==AUROSTD_NAN){return AUROSTD_NAN;}
    return convertUnits(d_tmp,units);
  }
  double T_S(const ChullPoint& point){
    if(!point.m_has_entry){return AUROSTD_NAN;}
    return T_S(point.m_entry);
  }
  double T_S(const aflowlib::_aflowlib_entry& entry){
    double d_tmp=entry.entropic_temperature;
    if(d_tmp==AUROSTD_NAN){return AUROSTD_NAN;}
    return d_tmp;
  }
  double EFA(const ChullPoint& point, char units){
    if(!point.m_has_entry){return AUROSTD_NAN;}
    return EFA(point.m_entry,units);
  }
  double EFA(const aflowlib::_aflowlib_entry& entry, char units){
    double d_tmp=entry.entropy_forming_ability;
    if(d_tmp==AUROSTD_NAN){return AUROSTD_NAN;}
    return convertUnits(d_tmp,units);
  }
  double isoMaxLatentHeat(const ChullPoint& point, double x, char units){
    if(!point.m_has_entry){return AUROSTD_NAN;}
    return isoMaxLatentHeat(point.m_entry,x,units);
  }
  double isoMaxLatentHeat(const aflowlib::_aflowlib_entry& entry, double x, char units){
    if(T_S(entry)==AUROSTD_NAN){return AUROSTD_NAN;}
    if(x<=0||x>=1.0){return AUROSTD_NAN;} //protect log()
    double d_tmp=((double)KBOLTZEV)*T_S(entry)*(x*log(x)+(1.0-x)*log(1.0-x));
    return convertUnits(d_tmp,units);
  }


  //for stoichiometries and results of an algorithm (distances from hull), use soft cutoff
  //otherwise, for above/below hull (stability from energy), use hard cutoff (per artificial points)
  //bool aurostd::greaterEqualZero(double val,bool soft_cutoff){return ( soft_cutoff? val>=ZERO_TOL : val>=0.0 );}
  //bool aurostd::lessEqualZero(double val,bool soft_cutoff){return ( soft_cutoff? val<=-ZERO_TOL : val<=0.0 );}

  bool subspaceBelongs(const xvector<int>& space,const xvector<int>& subspace){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " space=" << space << endl;
      cerr << __AFLOW_FUNC__ << " subspace=" << subspace << endl;
    }
    if(subspace.rows!=space.rows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Dimension mismatch between space and subspace");}
    if(sum(subspace)>sum(space)){return false;}  //cannot be of higher dimension
    for(int i=subspace.lrows;i<=subspace.urows;i++){if(subspace[i]==1&&space[i]==0){return false;}}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " is relevant!" << endl;}
    return true;
  }

  bool correctSignVerticalDistance(double dist_2_hull,bool should_be_positive) {
    if( should_be_positive && aurostd::notPositive(dist_2_hull,true,ZERO_TOL)){return false;}
    if(!should_be_positive && aurostd::notNegative(dist_2_hull,true,ZERO_TOL)){return false;}
    return true;
  }

  xvector<double> getTruncatedCoords(const xvector<double>& coords,const xvector<int>& elements_present) {
    if(coords.rows!=elements_present.rows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Reduction invalid, coords mismatch");}
    uint h_dim=sum(elements_present);
    xvector<double> red_coords(coords.lrows,coords.lrows+h_dim-1);
    vector<uint> relevant_indices=getRelevantIndices(elements_present);
    for(uint i=0,fl_size_i=relevant_indices.size();i<fl_size_i;i++){red_coords[i+red_coords.lrows]=coords[relevant_indices[i]];}
    return red_coords;
  }

  vector<uint> getRelevantIndices(const xvector<int>& elements_present) {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    vector<uint> relevant_indices;
    for(int i=elements_present.lrows;i<=elements_present.urows;i++){
      if(elements_present[i]==1){relevant_indices.push_back(i);}
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " relevant indices=";
      for(uint i=0,fl_size_i=relevant_indices.size();i<fl_size_i;i++){cerr << relevant_indices[i] << (i!=relevant_indices.size()-1?",":"");}
      cerr << " (elements_present=" << elements_present << ")" << endl;
    }
    return relevant_indices;
  }

  bool coordsIdentical(const xvector<double>& coords1,const xvector<double>& coords2){return identical(coords1,coords2,ZERO_TOL);}
} // namespace chull

//CO20180420 - moved to xStream (xclasses.cpp)
//namespace chull {
////--------------------------------------------------------------------------------
//// class ChullClassTemplate
////--------------------------------------------------------------------------------
////--------------------------------------------------------------------------------
//// constructor
////--------------------------------------------------------------------------------
//ChullClassTemplate::ChullClassTemplate() : p_FileMESSAGE(NULL),f_new_ofstream(false) {} //{free();}
//ChullClassTemplate::~ChullClassTemplate() {freeAll();}
//void ChullClassTemplate::free() {}
//void ChullClassTemplate::xStream::free() {
//  //if(1){  //keep this OFF. technically, we will have a memory leak, but FacetPoint deletes p_FileMESSAGE prematurely. really, we need shared pointers.
//  p_oss=NULL;
//  if(f_new_ofstream){delete p_FileMESSAGE;}  //first delete, then set to null
//  p_FileMESSAGE=NULL;
//  f_new_ofstream=false;
//  //}
//}
//void ChullClassTemplate::freeAll(){free();xStream::free();}
//void ChullClassTemplate::setOFStream(ofstream& FileMESSAGE){p_FileMESSAGE=&FileMESSAGE;}
//void ChullClassTemplate::setOSS(ostream& oss) {p_oss=&oss;}
//} // namespace chull

namespace chull {
  //--------------------------------------------------------------------------------
  // class ChullPointLight
  //--------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------
  // constructor
  //--------------------------------------------------------------------------------
  ChullPointLight::ChullPointLight(ostream& oss) : xStream(oss),m_initialized(false) {;}
  ChullPointLight::ChullPointLight(ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {;}
  ChullPointLight::ChullPointLight(const ChullPointLight& b) : xStream(*b.getOFStream(),*b.getOSS()) {copy(b);} // copy PUBLIC  //upcasting is allowed, works for ChullPointLight and ChullPoint

  ChullPointLight::~ChullPointLight() {xStream::free();free();}

  const ChullPointLight& ChullPointLight::operator=(const ChullPointLight& other) { //upcasting is allowed, works for ChullPointLight and ChullPoint
    if(this!=&other) {copy(other);}
    return *this;
  }
  bool ChullPointLight::operator<(const ChullPointLight& other) const {
    //NB: this is ALWAYS sorted in descending order of stoich, no need to make options for ascending order
    //but, sorts in ascending order for energy
    if(m_coords.lrows!=other.m_coords.lrows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_coords.lrows!=other.m_coords.lrows");}
    if(m_coords.rows!=other.m_coords.rows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_coords.rows!=other.m_coords.rows");}
    for(int i=m_coords.lrows;i<=m_coords.urows;i++){if(m_coords[i]!=other.m_coords[i]){return (m_coords[i]<other.m_coords[i]);}}
    return false;
  }

  void ChullPointLight::clear() {free();}  //clear PUBLIC
  void ChullPointLight::free() {
    m_initialized=false;
    m_coords.clear();
    m_has_stoich_coords=false;
    m_formation_energy_coord=false;
    m_is_artificial=false;
    m_has_entry=false;
    m_i_nary=AUROSTD_MAX_UINT;
    m_i_alloy=AUROSTD_MAX_UINT;
    h_coords.clear();
  }

  void ChullPointLight::copy(const ChullPointLight& b) {  //copy PRIVATE  //upcasting is allowed, works for ChullPointLight and ChullPoint
    xStream::copy(b);
    m_initialized=b.m_initialized;
    m_coords=b.m_coords;
    m_has_stoich_coords=b.m_has_stoich_coords;
    m_formation_energy_coord=b.m_formation_energy_coord;
    m_is_artificial=b.m_is_artificial;
    m_has_entry=b.m_has_entry;
    m_i_nary=b.m_i_nary;
    m_i_alloy=b.m_i_alloy;
    h_coords=b.h_coords;
  }

  double ChullPointLight::getLastCoord() const {return m_coords[m_coords.urows];}

  //--------------------------------------------------------------------------------
  // class ChullPoint
  //--------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------
  // constructor
  //--------------------------------------------------------------------------------
  ChullPoint::ChullPoint(ostream& oss,bool has_stoich_coords,bool formation_energy_coord,bool is_artificial) : ChullPointLight(oss) {initialize(has_stoich_coords,formation_energy_coord,is_artificial);}
  ChullPoint::ChullPoint(const xvector<double>& coord,ostream& oss,bool has_stoich_coords,bool formation_energy_coord,bool is_artificial) : ChullPointLight(oss) {initialize(coord,has_stoich_coords,formation_energy_coord,is_artificial);}
  ChullPoint::ChullPoint(const vector<string>& velements,const aflowlib::_aflowlib_entry& entry,ostream& oss,bool formation_energy_coord) : ChullPointLight(oss) {initialize(velements,entry,formation_energy_coord);}
  ChullPoint::ChullPoint(ofstream& FileMESSAGE,ostream& oss,bool has_stoich_coords,bool formation_energy_coord,bool is_artificial) : ChullPointLight(FileMESSAGE,oss) {initialize(has_stoich_coords,formation_energy_coord,is_artificial);}
  ChullPoint::ChullPoint(const xvector<double>& coord,ofstream& FileMESSAGE,ostream& oss,bool has_stoich_coords,bool formation_energy_coord,bool is_artificial) : ChullPointLight(FileMESSAGE,oss) {initialize(coord,has_stoich_coords,formation_energy_coord,is_artificial);}
  ChullPoint::ChullPoint(const vector<string>& velements,const aflowlib::_aflowlib_entry& entry,ofstream& FileMESSAGE,ostream& oss,bool formation_energy_coord) : ChullPointLight(FileMESSAGE,oss) {initialize(velements,entry,formation_energy_coord);}
  ChullPoint::ChullPoint(const ChullPoint& b) : xStream(*b.getOFStream(),*b.getOSS()),ChullPointLight(b) {copy(b);} // copy PUBLIC

  ChullPoint::~ChullPoint() {xStream::free();free();}

  const ChullPoint& ChullPoint::operator=(const ChullPoint& other) {
    if(this!=&other) {copy(other);}
    return *this;
  }

  void ChullPoint::clear() {free();}  //clear PUBLIC
  void ChullPoint::free() {
    ChullPointLight::free();
    m_entry.clear(); if(m_entry.vsg.size()==0){m_entry.vsg.push_back(NOSG);} if(m_entry.vsg2.size()==0){m_entry.vsg2.push_back(NOSG);}  //hack so it doesn't break with front(),back(),[0]
    s_coords.clear();
    c_coords.clear();
    m_elements_present.clear();
    cleanPointForHullTransfer();
  }

  void ChullPoint::copy(const ChullPoint& b) {  //copy PRIVATE
    //xStream::copy(b); //done inside ChullPointLight::copy()
    ChullPointLight::copy(b);
    m_entry=b.m_entry; if(m_entry.vsg.size()==0){m_entry.vsg.push_back(NOSG);} if(m_entry.vsg2.size()==0){m_entry.vsg2.push_back(NOSG);}  //hack so it doesn't break with front(),back(),[0]
    m_i_coord_group=b.m_i_coord_group;
    m_i_icsd=b.m_i_icsd;
    m_found_icsd=b.m_found_icsd;
    s_coords=b.s_coords;
    c_coords=b.c_coords;
    m_elements_present=b.m_elements_present;
    m_calculated_equivalent_entries=b.m_calculated_equivalent_entries;
    m_equivalent_entries=b.m_equivalent_entries;
    m_is_on_hull=b.m_is_on_hull;
    m_is_g_state=b.m_is_g_state;
    m_is_equivalent_g_state=b.m_is_equivalent_g_state;
    m_is_sym_equivalent_g_state=b.m_is_sym_equivalent_g_state;
    m_dist_2_hull=b.m_dist_2_hull;
    //[OBSOLETE - reduce by frac_vrt always! so use coord_group values]m_decomp_coefs=b.m_decomp_coefs;
    m_stability_criterion=b.m_stability_criterion;
    m_n_plus_1_enthalpy_gain=b.m_n_plus_1_enthalpy_gain;
  }

  bool ChullPoint::initialize(ostream& oss,bool has_stoich_coords,bool formation_energy_coord,bool is_artificial) {
    xStream::initialize(oss);
    return initialize(has_stoich_coords,formation_energy_coord,is_artificial);
  }

  bool ChullPoint::initialize(const xvector<double>& coord,ostream& oss,bool has_stoich_coords,bool formation_energy_coord,bool is_artificial) {
    xStream::initialize(oss);
    return initialize(coord,has_stoich_coords,formation_energy_coord,is_artificial);
  }

  bool ChullPoint::initialize(const vector<string>& velements,const aflowlib::_aflowlib_entry& entry,ostream& oss,bool formation_energy_coord) {
    xStream::initialize(oss);
    return initialize(velements,entry,formation_energy_coord);
  }

  bool ChullPoint::initialize(ofstream& FileMESSAGE,ostream& oss,bool has_stoich_coords,bool formation_energy_coord,bool is_artificial) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(has_stoich_coords,formation_energy_coord,is_artificial);
  }

  bool ChullPoint::initialize(const xvector<double>& coord,ofstream& FileMESSAGE,ostream& oss,bool has_stoich_coords,bool formation_energy_coord,bool is_artificial) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(coord,has_stoich_coords,formation_energy_coord,is_artificial);
  }

  bool ChullPoint::initialize(const vector<string>& velements,const aflowlib::_aflowlib_entry& entry,ofstream& FileMESSAGE,ostream& oss,bool formation_energy_coord) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(velements,entry,formation_energy_coord);
  }

  bool ChullPoint::initialize(bool has_stoich_coords,bool formation_energy_coord,bool is_artificial) {
    free();
    m_has_stoich_coords=has_stoich_coords;
    m_formation_energy_coord=formation_energy_coord;
    m_is_artificial=is_artificial;
    m_initialized=false;  //no point
    return m_initialized;
  }

  bool ChullPoint::initialize(const xvector<double>& coord,bool has_stoich_coords,bool formation_energy_coord,bool is_artificial) {
    //we start by setting the most general coords, m_coords, relating to highest d-hull
    //then, if we can, we derive stoich and composition coords
    //it may seem "backwards", as stoich and composition are most accesible to entries, but
    //this code is GENERAL (any type of coords, not just energy/stoich)
    free();
    initializeCoords(coord,formation_energy_coord);
    m_has_stoich_coords=has_stoich_coords;
    if(m_has_stoich_coords){setStoichCoords();}
    m_is_artificial=is_artificial;
    m_initialized=true;
    return m_initialized;
  }

  bool ChullPoint::initialize(const vector<string>& velements,const aflowlib::_aflowlib_entry& entry,bool formation_energy_coord) {
    //we start by setting the most general coords, m_coords, relating to highest d-hull
    //then, if we can, we derive stoich and composition coords
    //it may seem "backwards", as stoich and composition are most accesible to entries, but
    //this code is GENERAL (any type of coords, not just energy/stoich)
    free();
    initializeCoords(velements,entry,formation_energy_coord);
    m_has_stoich_coords=true;
    setStoichCoords();
    m_is_artificial=false;
    m_initialized=true;
    return m_initialized;
  }

  //use hard cutoff here
  bool ChullPoint::isWithinHalfHull(bool lower_hull) const {return (lower_hull ? aurostd::lessEqualZero(getLastCoord()) : aurostd::greaterEqualZero(getLastCoord()) );}

  bool ChullPoint::isGState() const {
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized point");}
    return m_is_on_hull || m_is_g_state || m_is_equivalent_g_state;
  }

  xvector<double> ChullPoint::getStoichiometricCoords() const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}
    if(m_coords.rows==1){ //trivial unary hull
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_coords.rows==1: generating null xvector" << endl;}
      if(!isUnary()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Found non-unary point on unary hull",_INPUT_UNKNOWN_);}
      xvector<double> r_coords=aurostd::null_xv<double>(); //get null vector //(m_coords.lrows,m_coords.lrows);  //just return vector length 1 with 0 inside
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " r_coords=" << r_coords << endl;}
      return r_coords;
    }
    xvector<double> r_coords(m_coords.lrows,m_coords.urows-1);
    for(int i=m_coords.lrows;i<=m_coords.urows-1;i++){r_coords[i]=m_coords[i];}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " r_coords=" << r_coords << endl;}
    return r_coords;
  }

  xvector<double> ChullPoint::getTruncatedReducedCoords(const xvector<int>& elements_present,vector_reduction_type vred) const {
    if(vred==frac_vrt || vred==no_vrt){return getTruncatedSCoords(elements_present);}
    else if(vred==gcd_vrt){return getTruncatedCCoords(elements_present,true);}
    else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown reduce mode",_INPUT_UNKNOWN_);}
    xvector<double> out;return out;
  }

  xvector<double> ChullPoint::getTruncatedSCoords(const xvector<int>& elements_present) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(!m_has_stoich_coords){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Non-stoich coordinates");}
    xvector<double> coords=getTruncatedCoords(s_coords,elements_present);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " truncated stoich for " << s_coords << " is " << coords << " (elements_present=" << elements_present << ")" << endl;}
    return coords;
  }

  xvector<double> ChullPoint::getTruncatedCCoords(const xvector<int>& elements_present,bool reduce) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(!m_has_stoich_coords){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Non-stoich coordinates");}
    xvector<double> coords=getTruncatedCoords(c_coords,elements_present);
    //DX20191125 [OBSOLETE] if(reduce){coords=aurostd::reduceByGCD(coords);}
    //DX20191125 [OBSOLETE] if(LDEBUG){cerr << __AFLOW_FUNC__ << " truncated " << (reduce?"reduced":"") << " comp for " << c_coords << " is " << coords << " (elements_present=" << elements_present << ")" << endl;}
    //DX20191125 [OBSOLETE] return coords;
    xvector<double> final_coords=coords; //DX20191125
    if(reduce){ aurostd::reduceByGCD(coords, final_coords);} //DX20191125 - new form of function
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " truncated " << (reduce?"reduced":"") << " comp for " << c_coords << " is " << final_coords << " (elements_present=" << elements_present << ")" << endl;} //DX20191125 - changed coords to final_coords
    return final_coords; //DX20191125 - changed coords to final_coords
  }

  xvector<double> ChullPoint::getReducedCCoords() const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(!m_has_stoich_coords){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Non-stoich coordinates");}
    xvector<double> coords; aurostd::reduceByGCD(c_coords,coords,ZERO_TOL); //DX20191125 - new form of function
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " reduced comp for " << c_coords << " is " << coords << endl;}
    return coords;
  }

  uint ChullPoint::loadXstructures(bool relaxed_only) {  //load relaxed only
    if(m_entry.prototype.find("POCC")!=string::npos){return false;} //POCC entries have no composition
    if(!pflow::loadXstructures(m_entry,*p_FileMESSAGE,*p_oss,relaxed_only)){return false;}
    return m_entry.vstr.size();
  }

  bool ChullPoint::getMostRelaxedXstructure(xstructure& xstr) const { //this is const!
    if(m_entry.prototype.find("POCC")!=string::npos){return false;} //POCC entries have no composition
    if(pflow::loadXstructureLibEntry(m_entry,xstr)){return true;}
    aflowlib::_aflowlib_entry entry; entry.auid=m_entry.auid; entry.aurl=m_entry.aurl; entry.vspecies=m_entry.vspecies;  //fast copy  //CO20221110 - pass vspecies here to populate xstructure
    if(!pflow::loadXstructures(entry,*p_FileMESSAGE,*p_oss,true)){return false;}
    if(entry.vstr.size()==1){xstr=entry.vstr[0]; return true;}
    return false;
  }

  //small get()'s of fundamental types get copies, otherwise const&
  uint ChullPoint::getDim() const {return m_coords.rows;}
  bool ChullPoint::isUnary() const {return m_i_nary==0;}
  double ChullPoint::getFormationEnthalpy() const {return H_f_atom(*this,_std_);} //m_entry.enthalpy_formation_atom;
  double ChullPoint::getEntropicTemperature() const {return T_S(*this);} //m_entry.entropic_temperature;
  double ChullPoint::getEntropyFormingAbility() const {return EFA(*this,_std_);} //m_entry.entropic_temperature;
  const vector<string>& ChullPoint::getVSG() const {return m_entry.vsg2;}
  const string& ChullPoint::getSG() const {return getVSG().back();}  //tight tolerance fine!  //vsg === LOOSE //vsg2 === TIGHT // doesn't make sense
  double ChullPoint::getDist2Hull(char units) const {
    if(m_dist_2_hull==AUROSTD_MAX_DOUBLE){return AUROSTD_MAX_DOUBLE;}
    if(m_formation_energy_coord){return convertUnits(m_dist_2_hull,units);}
    else{return m_dist_2_hull;}  //no unit conversions coded yet here
  }
  double ChullPoint::getStabilityCriterion(char units) const {
    if(m_stability_criterion==AUROSTD_MAX_DOUBLE){return AUROSTD_MAX_DOUBLE;}
    if(m_formation_energy_coord){return convertUnits(m_stability_criterion,units);}
    else{return m_stability_criterion;}  //no unit conversions coded yet here
  }
  double ChullPoint::getRelativeStabilityCriterion() const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(m_stability_criterion==AUROSTD_MAX_DOUBLE){return AUROSTD_MAX_DOUBLE;}
    if(aurostd::identical(getLastCoord(),0.0,ZERO_TOL)){return 0.0;}
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " m_stability_criterion=" << m_stability_criterion << endl;
      cerr << __AFLOW_FUNC__ << " getLastCoord()=" << getLastCoord() << endl;
    }
    return abs(m_stability_criterion/getLastCoord()); //abs() because they are generally opposite signs //delivers as decimal, show as percentage
  }
  double ChullPoint::getNPlus1EnthalpyGain(char units) const {
    if(m_n_plus_1_enthalpy_gain==AUROSTD_MAX_DOUBLE){return AUROSTD_MAX_DOUBLE;}
    if(m_formation_energy_coord){return convertUnits(m_n_plus_1_enthalpy_gain,units);}
    else{return m_n_plus_1_enthalpy_gain;}  //no unit conversions coded yet here
  }
  double ChullPoint::getEntropyStabilizationCoefficient(char units) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(getDist2Hull()==AUROSTD_MAX_DOUBLE){return AUROSTD_MAX_DOUBLE;}
    if(getEntropyFormingAbility()==AUROSTD_NAN){return AUROSTD_MAX_DOUBLE;}
    if(aurostd::zeroWithinTol(getEntropyFormingAbility(),ZERO_TOL)){return 0.0;} //protect from division by zero
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " dist2hull=" << getDist2Hull() << endl;
      cerr << __AFLOW_FUNC__ << " EFA=" << getEntropyFormingAbility() << endl;
    }
    return convertUnits(sqrt(getDist2Hull() / getEntropyFormingAbility()),units);
  }


  /// \brief sort the vertex indexes based on their angle around a central normal vector
  /// \param facet list of vertex indexes forming the facet
  /// \param facet_id id of one of the original facets (to use its already calculated normal vector)
  ///
  /// Sorting the vertex indexes avoid crossing lines when drawing or calculating the area or volume later.
  void ConvexHull::sortFacetVertices(vector<uint> &facet, const uint &facet_id){//HE20210510
    xvector<double> center(3,1);
    const uint num_points = facet.size();
    xvector<uint> index_list(num_points,1);
    xvector<double> angle_list(num_points,1);
    for (std::vector<uint>::const_iterator p_id = facet.begin(); p_id != facet.end(); ++p_id) center += m_points[*p_id].m_coords;

    center /= num_points;
    const xvector<double> start_vector = m_points[facet[0]].m_coords - center;
    const xvector<double> normal = m_facets[facet_id].m_normal;

    // first index is used for the start_vector, therefore the angle is set to 0.0
    angle_list[1] = 0;
    index_list[1] = facet[0];
    for (uint i=1; i<num_points; i++){
      const xvector<double> next_vector = m_points[facet[i]].m_coords - center;
      const double dot = aurostd::scalar_product(start_vector, next_vector);
      const double det = start_vector[1]*next_vector[2]*normal[3] + next_vector[1]*normal[2]*start_vector[3] + normal[1]*start_vector[2]*next_vector[3]
        - start_vector[3]*next_vector[2]*normal[1] - next_vector[3]*normal[2]*start_vector[1] - normal[3]*start_vector[2]*next_vector[1];
      angle_list[i+1] = atan2(det, dot);
      index_list[i+1] = facet[i];
    }
    aurostd::quicksort2(num_points, angle_list, index_list);
    for (uint i=0; i<num_points; i++) facet[i] = index_list[i+1];
  }

  /// @brief generates list of facets, if two neighboring facets are coplanar join them
  /// @param facet_collection output vector containing lists of vertex indexes
  /// @param angle_threshold max angle between two facts in radian to be still considered coplanar
  /// @note the facets contain only vertices
  void ConvexHull::getJoinedFacets(vector<vector<uint> > &facet_collection, const double angle_threshold) {//HE20210510
    bool LDEBUG=(false || XHOST.DEBUG);

    if (m_dim != 3) throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "facet joining is just available in 3D", _VALUE_RANGE_);
    vector <vector<uint> > raw_facets;
    vector<xvector<double> > normals;
    std::map<uint, vector<uint> >  point_neighbors;
    std::list<std::pair<uint, uint> > raw_join_list;
    vector<vector<uint> > join_list;
    vector<uint> remove_facet;
    vector<vector<uint> > raw_facet_collection;
    std::map<uint, uint> corner_check;

    // Collect information on each facet
    for (std::vector<ChullFacet>::const_iterator facet = h_facets.begin(); facet != h_facets.end(); ++facet){

      vector <uint> vertices;

      // Represent facet based on the vertices indexes
      for(std::vector<chull::FacetPoint>::const_iterator vert = facet->m_vertices.begin(); vert != facet->m_vertices.end(); ++vert){
        vertices.push_back(vert->ch_index);
      }
      raw_facets.push_back(vertices);
      normals.push_back(facet->m_normal);
    }


    // Build lookup for neighboring points
    // (Base point is included to make a check easier)
    for (std::vector<vector<uint> >::const_iterator facet = raw_facets.begin(); facet != raw_facets.end(); ++facet){
      for (std::vector<uint>::const_iterator ind = facet->begin(); ind != facet->end(); ++ind){
        std::copy(facet->begin(),facet->end(), std::inserter(point_neighbors[*ind],point_neighbors[*ind].end()));
      }
    }

    // Ensure that point_neighbors is unique
    for (std::map<uint, vector<uint> >::iterator point_neighbors_entry = point_neighbors.begin();
        point_neighbors_entry != point_neighbors.end(); ++point_neighbors_entry){
      std::sort(point_neighbors_entry->second.begin(),point_neighbors_entry->second.end());
      point_neighbors_entry->second.erase( std::unique( point_neighbors_entry->second.begin(), point_neighbors_entry->second.end() ), point_neighbors_entry->second.end() );
    }

    uint raw_facets_size = raw_facets.size();
    // Check for each facet, if their neighbors have an equivalent normal vector
    if (LDEBUG) cerr << __AFLOW_FUNC__ << " coplanar | angle | facets | n1 | n2" << endl;
    for (uint i1=0; i1<raw_facets_size; i1++){
      const vector<uint> base_facet = raw_facets[i1];
      for (uint i2=i1+1; i2<raw_facets_size; i2++){
        uint check=0;
        const vector<uint> compare_facet = raw_facets[i2];
        for (uint k=0; k<3; k++) {
          for (uint j=0; j<3; j++) {
            if (base_facet[k]==compare_facet[j]) check++;
          }
        }
        // neighbor share two vertices
        if (check != 2) continue;
        double check_angle = aurostd::angle(normals[i1], normals[i2]);
        if (check_angle<angle_threshold){
          if(LDEBUG) cerr << __AFLOW_FUNC__ << " YES | ";
          raw_join_list.push_back(std::make_pair(i1, i2));
          remove_facet.push_back(i1);
          remove_facet.push_back(i2);
        }
        else{
          if (LDEBUG) cerr << __AFLOW_FUNC__ << " NO  | ";
        }
        if (LDEBUG){
          cerr << check_angle << " | ";
          cerr << i1 << ", " << i2 << " | ";
          for (uint k=1; k<4; k++) cerr << normals[i1][k] << ", ";
          cerr << "| ";
          for (uint k=1; k<4; k++) cerr << normals[i2][k] << ", ";
          cerr << endl;
        }
      }
    }

    // Combine the joined pairs into complete facets
    while (raw_join_list.size()){
      std::pair<uint, uint> start=*raw_join_list.begin();
      vector<uint> new_facet;
      new_facet.push_back(start.first); new_facet.push_back(start.second);
      raw_join_list.erase(std::find(raw_join_list.begin(), raw_join_list.end(), start));
      std::vector<std::pair<uint, uint> > to_delete;
      uint found = 1;
      while (found){
        found = 0;
        to_delete.clear();
        for (std::list< std::pair<uint, uint> >::const_iterator next_ptr = raw_join_list.begin(); next_ptr != raw_join_list.end(); ++next_ptr) {
          std::pair<uint, uint> next = *next_ptr;
          if (std::find(new_facet.begin(), new_facet.end(), next.first) != new_facet.end()) {
            found ++;
            new_facet.push_back(next.second);
            to_delete.push_back(next);
          }
          else if (std::find(new_facet.begin(), new_facet.end(), next.second) != new_facet.end())  {
            found ++;
            new_facet.push_back(next.first);
            to_delete.push_back(next);
          }
        }
        for (std::vector<std::pair<uint, uint> >::const_iterator next = to_delete.begin(); next != to_delete.end(); ++next) {
          raw_join_list.erase(std::find(raw_join_list.begin(), raw_join_list.end(), *next));
        }
      }
      join_list.push_back(new_facet);
    }

    // Build the raw facet collection
    //HE20220319 sorting needed for consistent 3D rendering (the order defines the facet normal)
    for (uint i=0; i<raw_facets_size; i++){
      sortFacetVertices(raw_facets[i], i);
      if (std::find(remove_facet.begin(), remove_facet.end(), i) == remove_facet.end()) raw_facet_collection.push_back(raw_facets[i]);
    }
    // Add joined facets
    for (std::vector<vector<uint> >::const_iterator to_join = join_list.begin(); to_join != join_list.end(); ++to_join) {
      vector<uint> vertices;
      for (vector<uint>::const_iterator f_id = to_join->begin(); f_id != to_join->end(); ++f_id) {
        for (std::vector<uint>::const_iterator p_id = raw_facets[*f_id].begin(); p_id != raw_facets[*f_id].end(); ++p_id) {
          if (std::find(vertices.begin(), vertices.end(), *p_id) == vertices.end()) vertices.push_back(*p_id); // ensure vertices vector is unique
        }
      }
      // sort the vertices for std::set_difference(); point_neighbors are already sorted
      std::sort(vertices.begin(), vertices.end());

      // If a vertex has no outside neighbor it is removed
      vector<uint> vertices_to_remove;
      for (vector<uint>::const_iterator p_id = vertices.begin(); p_id != vertices.end(); ++p_id){
        std::vector<uint> diff_result;
        std::set_difference(point_neighbors[*p_id].begin(), point_neighbors[*p_id].end(), vertices.begin(), vertices.end(), std::inserter(diff_result, diff_result.end()));
        if (!diff_result.size()) vertices_to_remove.push_back(*p_id);
      }
      for (std::vector<uint>::const_iterator p_id = vertices_to_remove.begin(); p_id != vertices_to_remove.end(); ++p_id) {
        vertices.erase(std::find(vertices.begin(), vertices.end(), *p_id));
      }
      //HE20220319 never add facet with less than three vertices
      if (vertices.size()>=3) {
        sortFacetVertices(vertices, *to_join->begin());
        raw_facet_collection.push_back(vertices);
      }
    }

    // remove points that are on an edge and not a corner
    for (uint i_facet = 0; i_facet < raw_facet_collection.size(); i_facet++) {
      for (uint i_vert = 0; i_vert < raw_facet_collection[i_facet].size(); i_vert++) {
        corner_check[raw_facet_collection[i_facet][i_vert]]++;
      }
    }

    for (uint i_facet = 0; i_facet < raw_facet_collection.size(); i_facet++) {
      vector<uint> vec_vertices;
      for (uint i_vert = 0; i_vert < raw_facet_collection[i_facet].size(); i_vert++) {
        if (corner_check[raw_facet_collection[i_facet][i_vert]]>2)  vec_vertices.push_back(raw_facet_collection[i_facet][i_vert]);
      }
      facet_collection.push_back(vec_vertices);
    }
    //HE20220413 START
    // Check for ghost facets with less than 3 vertices
    // (they can be created when joining facets that are not perfectly coplanar)
    // start with the largest index when removing, as vector will shrink and indexes would change
    for (uint i_facet=facet_collection.size()-1; i_facet > 0; i_facet--){
      if (facet_collection[i_facet].size() <= 2) {
        facet_collection.erase(facet_collection.begin() + i_facet);
      }
    }
    //HE20220413 END
  }
  //since we don't check ALL attributes of entry, then we weed out MORE
  //entries existing in different catalogs will not be strictly identical
  //avoid this by comparing only the most pertinent information
  bool ChullPoint::entryIdentical(const aflowlib::_aflowlib_entry& other) const {
    return (m_entry.compound==other.compound) && 
      (m_entry.prototype==other.prototype) && //believe it or not, compound + prototype is PROBABLY enough, but let's be sure
      (aurostd::identical(H_f_atom(m_entry,_std_),H_f_atom(other,_std_),ENERGY_TOL)); //&&
    //[UNRELIABLE BECAUSE OF TOLERANCE](m_entry.vsg2==other.vsg2);
  }

  void ChullPoint::initializeCoords(const xvector<double>& coord,bool formation_energy_coord) {
    setGenCoords(coord,formation_energy_coord);
    setHullCoords();  //default, will change later with increasing dims
  }

  void ChullPoint::initializeCoords(const vector<string>& velements,const aflowlib::_aflowlib_entry& entry,bool formation_energy_coord) {
    addEntry(entry);
    setGenCoords(velements,m_entry,formation_energy_coord);
    setHullCoords();  //default, will change later with increasing dims
    m_has_entry=true;
  }

  void ChullPoint::addEntry(const aflowlib::_aflowlib_entry& entry) {
    m_entry=entry;
    if(CORRECT_BAD_DATABASE){m_entry.correctBadDatabase(*p_FileMESSAGE,true,*p_oss);} //verbose //carried over from apennsy (SC)
  }

  void ChullPoint::setGenCoords(const xvector<double>& coord,bool formation_energy_coord) {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    m_coords=coord;
    m_formation_energy_coord=formation_energy_coord;
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_within_half_hull=" << isWithinHalfHull(m_formation_energy_coord) << endl;}  //m_formation_energy_coord===lower_hull
  }

  void ChullPoint::setGenCoords(const vector<string>& velements,const aflowlib::_aflowlib_entry& entry,bool formation_energy_coord) {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(entry.vcomposition.size()==0&&entry.vstoichiometry.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No vcomposition or vstoichiometry found for entry.auid="+entry.auid,_RUNTIME_ERROR_);}
    xvector<double> coord(velements.size());
    bool found=false; //not necessary, we already check in entryValid()
    double c_sum=0.0;
    if(entry.vcomposition.size()>0){
      for(uint i=0,fl_size_i=entry.vcomposition.size();i<fl_size_i;i++){c_sum+=entry.vcomposition[i];}  //derive stoich exactly!
      if(c_sum==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"c_sum==0 (entry.auid="+entry.auid+",entry.aurl="+entry.aurl+")",_RUNTIME_ERROR_);}
      for(uint i=0,fl_size_i=velements.size();i<fl_size_i-1;i++){
        found=false;
        for(uint j=0,fl_size_j=entry.vspecies.size();j<fl_size_j && !found;j++){
          if(velements[i]==entry.vspecies[j]){
            coord[i+coord.lrows]=entry.vcomposition[j]/c_sum;
            found=true;
          }
        }
        //[might be from lower hull]if(!found){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"element not found: "+velements[i],_RUNTIME_ERROR_);}
      }
    }else{  //pocc structures have no vcomposition, only vstoichiometry
      //get exact fraction
      int numerator=0,denominator=0;
      double stoich=0.0;
      vector<double> vstoich;
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " converting stoich[aurl=" << entry.aurl << "]=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(entry.vstoichiometry),",") << " to fractions" << endl;}
      for(uint i=0,fl_size_i=entry.vstoichiometry.size();i<fl_size_i;i++){  //derive stoich exactly!
        aurostd::double2fraction(entry.vstoichiometry[i],numerator,denominator,DEFAULT_POCC_STOICH_TOL);  //ZERO_TOL*10 is preferred, but due to buy in GetStoich(), use DEFAULT_POCC_STOICH_TOL instead (for now) //ZERO_TOL=1e-8, make slightly looser than what's written to aflowlib.out
        stoich=(double)numerator/(double)denominator;
        c_sum+=stoich;
        vstoich.push_back(stoich);
      }
      if(c_sum==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"c_sum==0 (entry.auid="+entry.auid+",entry.aurl="+entry.aurl+")",_RUNTIME_ERROR_);}
      if(!aurostd::identical(c_sum,1.0,ZERO_TOL)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"c_sum!=1 (entry.auid="+entry.auid+",entry.aurl="+entry.aurl+")",_RUNTIME_ERROR_);}
      for(uint i=0,fl_size_i=velements.size();i<fl_size_i-1;i++){
        found=false;
        for(uint j=0,fl_size_j=entry.vspecies.size();j<fl_size_j && !found;j++){
          if(velements[i]==entry.vspecies[j]){
            coord[i+coord.lrows]=vstoich[j]/c_sum;
            found=true;
          }
        }
        //[might be from lower hull]if(!found){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"element not found: "+velements[i],_RUNTIME_ERROR_);}
      }
      if(0){  //entry.vstoichiometry has write round-offs, better to derive fraction exactly
        for(uint i=0,fl_size_i=entry.vstoichiometry.size();i<fl_size_i;i++){c_sum+=entry.vstoichiometry[i];}  //derive stoich exactly!
        if(c_sum==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"c_sum==0 (entry.auid="+entry.auid+",entry.aurl="+entry.aurl+")",_RUNTIME_ERROR_);}
        for(uint i=0,fl_size_i=velements.size();i<fl_size_i-1;i++){
          found=false;
          for(uint j=0,fl_size_j=entry.vspecies.size();j<fl_size_j && !found;j++){
            if(velements[i]==entry.vspecies[j]){
              coord[i+coord.lrows]=entry.vstoichiometry[j]/c_sum;
              found=true;
            }
          }
          //[might be from lower hull]if(!found){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"element not found: "+velements[i],_RUNTIME_ERROR_);}
        }
      }
    }
    if(formation_energy_coord){coord[coord.urows]=H_f_atom(entry);} //entry.enthalpy_formation_atom
    else{coord[coord.urows]=entry.entropic_temperature;}  //entropic temperature is positive for stable compounds (we want upper half convex-hull)
    setGenCoords(coord,formation_energy_coord);
  }

  vector<uint> ChullPoint::getRelevantIndices(const xvector<int>& elements_present) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}
    if(!m_has_stoich_coords){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Non-stoich coordinates");}
    if(elements_present.rows!=s_coords.rows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Dimension mismatch between point and elements_present");}
    for(int i=elements_present.lrows;i<=elements_present.urows;i++){
      if(elements_present[i]==0 && aurostd::nonZeroWithinTol(s_coords[i],ZERO_TOL)){
        stringstream message;
        message << "Attempting to reduce non-zero coord (i=" << i << ",s_coords=" << s_coords << "), ";
        message << "elements_present=" << elements_present;
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
      }
    }
    return chull::getRelevantIndices(elements_present);
  }

  void ChullPoint::setStoichCoords() {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(!m_has_stoich_coords){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Non-stoich coordinates: aurl="+m_entry.aurl);}
    double c_sum=0.0; //concentration sum
    xvector<double> stoich(m_coords.urows,m_coords.lrows);
    xvector<int> elements_present(m_coords.urows,m_coords.lrows);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_coords=" << m_coords << endl;}
    for(int j=m_coords.lrows;j<=m_coords.urows-1;j++){
      if(std::signbit(m_coords[j])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Negative stoich coordinate found: aurl="+m_entry.aurl);} //no negative numbers in stoich coordinates, only energy
      stoich[j]=m_coords[j];
      if(aurostd::nonZeroWithinTol(m_coords[j],ZERO_TOL)){elements_present[j]=1;}
      c_sum+=m_coords[j];
    }
    stoich[stoich.urows]=(1.0-c_sum); //hidden dimension
    if(std::signbit(stoich[stoich.urows])){
      //necessary check now because POCC entries have no composition, only stoich, so write-out errors will be prevalent
      if(aurostd::zeroWithinTol(stoich[stoich.urows],ZERO_TOL)){stoich[stoich.urows]=0.0;}  //only zero out for the last coord, as it's derived from the subtraction of the others
      else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Negative stoich coordinate found: aurl="+m_entry.aurl);}
    }  //no negative numbers
    if(aurostd::nonZeroWithinTol(stoich[stoich.urows],ZERO_TOL)){elements_present[elements_present.urows]=1;}   //check if nary++
    s_coords=stoich;
    c_coords=s_coords;
    if(m_has_entry){
      c_sum=0.0;
      for(uint i=0,fl_size_i=m_entry.vcomposition.size();i<fl_size_i;i++){c_sum+=m_entry.vcomposition[i];}  //essentially UNDO s_coords calculation
      c_coords*=c_sum;
    }
    m_elements_present=elements_present;
    m_i_nary=sum(elements_present)-1;
    setHullCoords(elements_present);  //default, will change later with increasing dims
  }

  void ChullPoint::setHullCoords() {setHullCoords(m_coords);}  //default to m_coords
  void ChullPoint::setHullCoords(const xvector<double>& coords) {h_coords=coords;}
  void ChullPoint::setHullCoords(const xvector<int>& elements_present) {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " elements_present=" << elements_present << endl;
      cerr << __AFLOW_FUNC__ << " s_coords=" << s_coords << endl;
      cerr << __AFLOW_FUNC__ << " c_coords=" << c_coords << endl;
    }
    xvector<double> coords=chull::getTruncatedCoords(m_coords,elements_present);
    //vector<uint> relevant_indices=getRelevantIndices(elements_present);
    //for(uint i=0,fl_size_i=relevant_indices.size();i<fl_size_i;i++){coords[i+coords.lrows]=m_coords[relevant_indices[i]];}
    coords[coords.urows]=getLastCoord(); //overwrite last coord appropriately
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " setting h_coords=" << coords << endl;}
    setHullCoords(coords);
  }

  void ChullPoint::reduceCoords(const xvector<int>& elements_present) {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);

    //redoing a bit of ChullPoint::initialize()
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_coords(pre )=" << m_coords << endl;}
    xvector<double> coords=chull::getTruncatedCoords(m_coords,elements_present);
    coords[coords.urows]=getLastCoord(); //overwrite last coord appropriately
    initializeCoords(coords,m_formation_energy_coord);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_coords(post)=" << m_coords << endl;}
    if(m_has_stoich_coords){setStoichCoords();}
  }

  void ChullPoint::cleanPointForHullCalc() {h_coords.clear();}
  void ChullPoint::cleanPointForHullTransfer() {
    m_i_alloy=AUROSTD_MAX_UINT;
    m_i_coord_group=AUROSTD_MAX_UINT;
    m_i_icsd=AUROSTD_MAX_UINT;
    m_found_icsd=false;
    m_calculated_equivalent_entries=false;
    m_equivalent_entries.clear();
    m_is_on_hull=false;
    m_is_g_state=false;
    m_is_equivalent_g_state=false;
    m_is_sym_equivalent_g_state=false;
    m_dist_2_hull=AUROSTD_MAX_DOUBLE;
    //[OBSOLETE - reduce by frac_vrt always! so use coord_group values]m_decomp_coefs.clear();
    m_stability_criterion=AUROSTD_MAX_DOUBLE;
    m_n_plus_1_enthalpy_gain=AUROSTD_MAX_DOUBLE;
    cleanPointForHullCalc();
  }
} // namespace chull

namespace chull {
  //--------------------------------------------------------------------------------
  // constructor
  //--------------------------------------------------------------------------------
  FacetPoint::FacetPoint() {free();}
  FacetPoint::FacetPoint(const ChullPointLight& point,uint index){initialize(point,index);}  //need BOTH point and index, otherwise, just use point/index independently
  FacetPoint::FacetPoint(const FacetPoint& b) {copy(b);}  // copy PUBLIC
  FacetPoint::~FacetPoint() {free();}

  const FacetPoint& FacetPoint::operator=(const FacetPoint& other) {
    if(this!=&other) {copy(other);}
    return *this;
  }

  //simple sort, so ONLY sort facets from the same hull
  bool FacetPoint::operator<(const FacetPoint& other) const {return ch_index<other.ch_index;}
  void FacetPoint::clear() {FacetPoint a;copy(a);} //clear PUBLIC
  void FacetPoint::free() {
    m_initialized=false;
    ch_index=AUROSTD_MAX_UINT;
    ch_point.clear();
  }

  void FacetPoint::copy(const FacetPoint& b) { //copy PRIVATE
    m_initialized=b.m_initialized;
    ch_index=b.ch_index;
    ch_point=b.ch_point;
  }

  void FacetPoint::initialize(const ChullPointLight& point,uint index) {
    ch_point=point;
    ch_index=index;
    m_initialized=true;
  }
} // namespace chull

//nice sorters for points we know sit on a thermo hull (stoich coords + energy dimension)
namespace chull {
  bool sortThermoPoints::operator() (const FacetPoint& fpi,const FacetPoint& fpj) const {
    const ChullPointLight& ci=fpi.ch_point;
    const ChullPointLight& cj=fpj.ch_point;
    return (*this).operator()(ci,cj);
  }

  bool sortThermoPoints::operator() (const ChullPointLight& ci,const ChullPointLight& cj) const {  //upcasting is allowed, works for ChullPointLight and ChullPoint
    if(!(ci.m_initialized && cj.m_initialized)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Points not initialized");}
    //do not first sort binaries from ternaries, screws up facet sorting
    //keep sorting based on concentration of elements in relative order
    //if(ci.m_has_stoich_coords&&cj.m_has_stoich_coords){
    //  if(ci.m_i_nary!=cj.m_i_nary){return ci.m_i_nary<cj.m_i_nary;} //binaries before ternaries
    //}
    if(ci.m_coords.rows!=cj.m_coords.rows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Dimension mismatch among points");}
    for(int i=ci.m_coords.lrows;i<=ci.m_coords.urows-1;i++){
      if(ci.m_coords[i]!=cj.m_coords[i]){
        return m_sort_stoich_ascending ? (ci.m_coords[i]<cj.m_coords[i]) : (ci.m_coords[i]>cj.m_coords[i]);
      }
    }
    if(ci.getLastCoord()!=cj.getLastCoord()){
      return m_sort_energy_ascending ? (ci.getLastCoord()<cj.getLastCoord()) : (ci.getLastCoord()>cj.getLastCoord());
    }
    return false; //true; //breaks if true
  }
  //[CO20221111 - VERY SLOW]bool sortLIB2Entries::operator() (const aflowlib::_aflowlib_entry i_entry,const aflowlib::_aflowlib_entry j_entry) const {
  //[CO20221111 - VERY SLOW]  bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
  //[CO20221111 - VERY SLOW]  if(i_entry.catalog!=j_entry.catalog){return i_entry.catalog<j_entry.catalog;}
  //[CO20221111 - VERY SLOW]  if(i_entry.prototype!=j_entry.prototype){return i_entry.prototype<j_entry.prototype;}
  //[CO20221111 - VERY SLOW]  if(i_entry.catalog=="LIB2" && i_entry.vspecies.size()==1 && j_entry.vspecies.size()==1){
  //[CO20221111 - VERY SLOW]    vector<string> i_vspecies=i_entry.getSpeciesAURL(*p_FileMESSAGE,*p_oss); //in this case, species from AUID is NOT the real species
  //[CO20221111 - VERY SLOW]    vector<string> j_vspecies=j_entry.getSpeciesAURL(*p_FileMESSAGE,*p_oss); //in this case, species from AUID is NOT the real species
  //[CO20221111 - VERY SLOW]    if(i_vspecies.size()!=j_vspecies.size()){return i_vspecies.size()<j_vspecies.size();} //effectively sorting by catalog
  //[CO20221111 - VERY SLOW]    if(LDEBUG){
  //[CO20221111 - VERY SLOW]      cerr << __AFLOW_FUNC__ << " i_entry.auid=" << i_entry.auid << " i_entry.aurl=" << i_entry.aurl << " i_vspecies=" << aurostd::joinWDelimiter(i_vspecies,",") << endl;
  //[CO20221111 - VERY SLOW]      cerr << __AFLOW_FUNC__ << " j_entry.auid=" << j_entry.auid << " j_entry.aurl=" << j_entry.aurl << " j_vspecies=" << aurostd::joinWDelimiter(j_vspecies,",") << endl;
  //[CO20221111 - VERY SLOW]    }
  //[CO20221111 - VERY SLOW]    uint k=0;
  //[CO20221111 - VERY SLOW]    bool i_within_hull=true;for(k=0;k<i_vspecies.size()&&i_within_hull;k++){if(!aurostd::WithinList(m_velements,i_vspecies[k])){i_within_hull=false;}}
  //[CO20221111 - VERY SLOW]    bool j_within_hull=true;for(k=0;k<j_vspecies.size()&&j_within_hull;k++){if(!aurostd::WithinList(m_velements,j_vspecies[k])){j_within_hull=false;}}
  //[CO20221111 - VERY SLOW]    if(LDEBUG){cerr << __AFLOW_FUNC__ << " i_within_hull=" << i_within_hull << " j_within_hull=" << j_within_hull << endl;}
  //[CO20221111 - VERY SLOW]    if(i_within_hull!=j_within_hull){return i_within_hull;}
  //[CO20221111 - VERY SLOW]    for(k=0;k<i_vspecies.size();k++){if(i_vspecies[k]!=j_vspecies[k]){return i_vspecies[k]<j_vspecies[k];}}
  //[CO20221111 - VERY SLOW]  }
  //[CO20221111 - VERY SLOW]  //
  //[CO20221111 - VERY SLOW]  return i_entry.aurl<j_entry.aurl; //they cannot be identical as they are physical directories, sort like `ls`
  //[CO20221111 - VERY SLOW]}
  _aflowlib_entry_LIB2sorting::_aflowlib_entry_LIB2sorting(aflowlib::_aflowlib_entry& entry,uint index,vector<string>& velements_chull,ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_index(index) {
      //[CO20221112 - no pointers]m_entry=&entry;
      m_auid=entry.auid;
      m_aurl=entry.aurl;
      m_catalog=entry.catalog;
      m_prototype=entry.prototype;
      m_nspecies=entry.nspecies;
      m_velements_chull=velements_chull;
      m_species_AURL=entry.getSpeciesAURL(*p_FileMESSAGE,*p_oss);
    }
  _aflowlib_entry_LIB2sorting::~_aflowlib_entry_LIB2sorting(){xStream::free();m_species_AURL.clear();}
  bool _aflowlib_entry_LIB2sorting::operator<(const _aflowlib_entry_LIB2sorting& other) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(m_catalog!=other.m_catalog){return m_catalog<other.m_catalog;}
    if(m_prototype!=other.m_prototype){return m_prototype<other.m_prototype;}
    if(m_catalog=="LIB2" && m_nspecies==1 && other.m_nspecies==1){
      if(m_species_AURL.size()!=other.m_species_AURL.size()){return m_species_AURL.size()<other.m_species_AURL.size();} //effectively sorting by catalog
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " m_auid=" << m_auid << " m_aurl=" << m_aurl << " m_species_AURL=" << aurostd::joinWDelimiter(m_species_AURL,",") << endl;
        cerr << __AFLOW_FUNC__ << " other.m_auid=" << other.m_auid << " other.m_aurl=" << other.m_aurl << " other.m_species_AURL=" << aurostd::joinWDelimiter(other.m_species_AURL,",") << endl;
      }
      uint k=0;
      bool i_within_hull=true;for(k=0;k<m_species_AURL.size()&&i_within_hull;k++){if(!aurostd::WithinList(m_velements_chull,m_species_AURL[k])){i_within_hull=false;}}
      bool j_within_hull=true;for(k=0;k<other.m_species_AURL.size()&&j_within_hull;k++){if(!aurostd::WithinList(m_velements_chull,other.m_species_AURL[k])){j_within_hull=false;}}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " i_within_hull=" << i_within_hull << " j_within_hull=" << j_within_hull << endl;}
      if(i_within_hull!=j_within_hull){return i_within_hull;}
      for(k=0;k<m_species_AURL.size();k++){if(m_species_AURL[k]!=other.m_species_AURL[k]){return m_species_AURL[k]<other.m_species_AURL[k];}}
    }
    //
    return m_aurl<other.m_aurl; //they cannot be identical as they are physical directories, sort like `ls`
  }
} //namespace chull

namespace chull {
  //--------------------------------------------------------------------------------
  // class ChullFacet
  //--------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------
  // constructor
  //--------------------------------------------------------------------------------
  ChullFacet::ChullFacet(ostream& oss) : xStream(oss),m_initialized(false) {free();}
  ChullFacet::ChullFacet(ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {free();}
  ChullFacet::ChullFacet(const ChullFacet& b) : xStream(*b.getOFStream(),*b.getOSS()) {copy(b);} // copy PUBLIC
  ChullFacet::~ChullFacet() {xStream::free();free();}

  const ChullFacet& ChullFacet::operator=(const ChullFacet& other) {
    if(this!=&other) {copy(other);}
    return *this;
  }

  bool ChullFacet::operator<(const ChullFacet& other) const {
    if(m_vertices.size()!=other.m_vertices.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Dimension mismatch between facet points");}

    //simply sort by m_point (indices)
    for(uint i=0,fl_size_i=m_vertices.size();i<fl_size_i;i++){
      if(!(m_vertices[i].m_initialized && other.m_vertices[i].m_initialized)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized facetpoint");}
      if(m_vertices[i].ch_index!=other.m_vertices[i].ch_index){return m_vertices[i].ch_index<other.m_vertices[i].ch_index;}
    }
    return false;
  }

  void ChullFacet::clear() {ChullFacet a;copy(a);}  //clear PRIVATE
  void ChullFacet::free() {
    m_initialized=false;
    m_vertices.clear();
    m_dim=AUROSTD_MAX_UINT;
    m_has_stoich_coords=false;
    m_formation_energy_coord=false;
    m_content=AUROSTD_MAX_DOUBLE;
    m_directive_vectors.clear();
    m_normal.clear();
    m_offset=AUROSTD_MAX_DOUBLE;
    m_facet_centroid.clear();
    m_hull_reference.clear();
    m_is_hypercollinear=false;
    m_is_vertical=false;
    m_is_artificial=false;
    m_in_lower_hemisphere=false;
    m_ridges.clear();
    //see xStream::free()
    //p_oss=NULL;
    //p_FileMESSAGE=NULL;
    //f_new_ofstream=false;
    cleanFacet();
  }

  void ChullFacet::copy(const ChullFacet& b) {  //copy PRIVATE
    xStream::copy(b);
    m_initialized=b.m_initialized;
    m_vertices.clear(); for(uint i=0,fl_size_i=b.m_vertices.size();i<fl_size_i;i++){m_vertices.push_back(b.m_vertices[i]);}
    m_dim=b.m_dim;
    m_has_stoich_coords=b.m_has_stoich_coords;
    m_formation_energy_coord=b.m_formation_energy_coord;
    m_content=b.m_content;
    m_directive_vectors.clear(); for(uint i=0,fl_size_i=b.m_directive_vectors.size();i<fl_size_i;i++){m_directive_vectors.push_back(b.m_directive_vectors[i]);}
    m_normal=b.m_normal;
    m_offset=b.m_offset;
    m_facet_centroid=b.m_facet_centroid;
    m_hull_reference=b.m_hull_reference;
    m_is_hypercollinear=b.m_is_hypercollinear;
    m_is_vertical=b.m_is_vertical;
    m_is_artificial=b.m_is_artificial;
    m_in_lower_hemisphere=b.m_in_lower_hemisphere;
    f_visited=b.f_visited;
    f_outside_set.clear(); for(uint i=0,fl_size_i=b.f_outside_set.size();i<fl_size_i;i++){f_outside_set.push_back(b.f_outside_set[i]);}
    f_furthest_point=b.f_furthest_point;
    m_ridges.clear(); for(uint i=0,fl_size_i=b.m_ridges.size();i<fl_size_i;i++){m_ridges.push_back(b.m_ridges[i]);}
    f_neighbors.clear(); for(uint i=0,fl_size_i=b.f_neighbors.size();i<fl_size_i;i++){f_neighbors.push_back(b.f_neighbors[i]);}
  }

  bool ChullFacet::shareRidge(const ChullFacet& other) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    bool match;
    for(uint i=0,fl_size_i=m_ridges.size();i<fl_size_i;i++){
      const vector<uint>& ridge_indices1=m_ridges[i].getCHIndices();
      for(uint j=0,fl_size_j=other.m_ridges.size();j<fl_size_j;j++){
        const vector<uint>& ridge_indices2=other.m_ridges[j].getCHIndices();
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " comparing ";
          for(uint ri=0,fl_size_ri=ridge_indices1.size();ri<fl_size_ri;ri++){cerr << ridge_indices1[ri] << " ";}
          cerr << " vs. ";
          for(uint rj=0,fl_size_rj=ridge_indices1.size();rj<fl_size_rj;rj++){cerr << ridge_indices2[rj] << " ";}
        }
        match=(ridge_indices1==ridge_indices2);
        if(LDEBUG){cerr << (match?"MATCH":"NO") << endl;}
        if(match){return true;}
      }
    }
    return false;
  }

  bool ChullFacet::isPointOnFacet(const FacetPoint& fp) const {
    if(!fp.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized facetpoint");}
    return isPointOnFacet(fp.ch_index);
  }

  bool ChullFacet::isPointOnFacet(uint i_point) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(m_vertices.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Facet has no vertices");}
    //if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized facet");}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " checking if point[" << i_point << "] is on this facet" << endl;}
    for(uint i=0,fl_size_i=m_vertices.size();i<fl_size_i;i++){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_vertices[i].ch_index==" << m_vertices[i].ch_index << " ?= i_point==" << i_point << endl;}
      if(m_vertices[i].ch_index==i_point){return true;}
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " point[" << i_point << "] is not on this facet" << endl;}
    return false;
  }

  bool ChullFacet::isPointOutside(const FacetPoint& f_point) const {
    if(!f_point.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized facetpoint");}
    if(isPointOnFacet(f_point)){return false;}
    return isPointOutside(f_point.ch_point);
  }

  bool ChullFacet::isPointOutside(const ChullPointLight& point) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Facet not initialized");}
    double dist_point=getSignedPointPlaneDistance(point);
    bool is_outside=false;
    //special case no longer needed, but we keep here just in case
    //if(m_is_vertical){
    //  if(LDEBUG){cerr << __AFLOW_FUNC__ << " performing special VERTICAL facet check for outside-ness" << endl;}
    //  xvector<double> ref=m_hull_reference;
    //  for(int i=ref.lrows;i<=ref.urows-1;i++){ref[i]=1.0/m_dim;}
    //  if(LDEBUG){cerr << __AFLOW_FUNC__ << " ref=" << ref << endl;}
    //  double dist_ref=getSignedPointPlaneDistance(point);
    //  is_outside=(std::signbit(dist_point)!=std::signbit(dist_ref) && aurostd::nonZeroWithinTol(dist_point,ZERO_TOL));
    //  if(LDEBUG){cerr << __AFLOW_FUNC__ << " sign of distance indicates point is " << (is_outside?"OUTSIDE":"INSIDE") << " hull" << endl;}
    //  return is_outside;
    //}
    //if(LDEBUG){cerr << __AFLOW_FUNC__ << " performing normal check for outside-ness (not vertical facet)" << endl;}
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " looking at facet with vertices: ";
      for(uint i=0,fl_size_i=m_vertices.size();i<fl_size_i;i++){
        cerr << m_vertices[i].ch_point.m_coords << " | ";
      }
      cerr << "m_normal=" << m_normal << endl;
      cerr << __AFLOW_FUNC__ << " checking if point[h_coords=" << point.h_coords << "] is on this facet" << endl;
      cerr << __AFLOW_FUNC__ << " dist=" << dist_point << endl;
    }     
    //point is on facet
    if(aurostd::zeroWithinTol(dist_point,ZERO_TOL)){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " point appears to be right on top of facet" << endl;}
      return false;
    }
    is_outside=(std::signbit(dist_point));
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " sign of distance indicates point is " << (is_outside?"OUTSIDE":"INSIDE") << " hull" << endl;}
    return is_outside;                //shortcut with inward-aligned normal
  }

  //sign depends on normal, if normal vector and point are in the same half-space, point-plane distance is positive, negative otherwise
  double ChullFacet::getSignedPointPlaneDistance(const ChullPointLight& point) const {return getSignedPointPlaneDistance(point.h_coords);}
  double ChullFacet::getSignedPointPlaneDistance(const xvector<double>& point) const {
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Facet not initialized");}
    if(point.rows!=m_vertices[0].ch_point.h_coords.rows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Dimension mismatch between point and facet");}
    xvector<double> diff=point-m_vertices[0].ch_point.h_coords;
    return scalar_product(m_normal,diff);
  }

  double ChullFacet::getSignedVerticalDistanceToZero(const ChullPointLight& point) const {
    //get energy of facet at stoichiometry of input (point)
    return getSignedVerticalDistanceToZero(point.h_coords);
  }

  double ChullFacet::getSignedVerticalDistanceToZero(const xvector<double>& point) const {
    //get energy of facet at stoichiometry of input (point)
    //distance from 0 to point on facet (vertical projection)
    //negative for projections BELOW 0, positive for projections ABOVE 0
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized facet");}
    if(m_normal.rows!=point.rows){
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " normal=" << m_normal << endl;
        cerr << __AFLOW_FUNC__ << " point=" << point << endl;
      }
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Dimension mismatch between facet["+aurostd::utype2string(m_normal.rows)+"] and point["+aurostd::utype2string(point.rows)+"]");
    }
    double dist=-m_offset;
    for(int i=point.lrows;i<=point.urows-1;i++){dist-=m_normal[i]*point[i];} //note that we don't go to rows (because we assume it is zero)
    dist/=m_normal[m_normal.urows];
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " point=" << point << " projected onto facet with normal=" << m_normal << " and offset=" << m_offset << endl;
      cerr << __AFLOW_FUNC__ << " dist=" << dist << endl;
    }
    return dist;
  }

  double ChullFacet::getSignedVerticalDistance(const ChullPointLight& point) const {
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized facet");}
    if(!point.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized point");}
    return getSignedVerticalDistance(point.h_coords);
  }

  double ChullFacet::getSignedVerticalDistance(const xvector<double>& point) const {
    //NB: we approximate the facet to be a hyperplane
    //also, this is the vertical distance, not the shortest (vertical is chemically meaningful)
    //to find the true distance between point and facet, must use quadratic programming
    //https://www.mathworks.com/matlabcentral/answers/107595-how-can-i-find-the-minimum-distance-from-convex-boundary
    //https://stackoverflow.com/questions/18230259/computing-distance-from-a-point-to-a-triangulation-in-3d-with-matlab
    //this is H_hull - H_f, which is negative for points on top of the hull (lower hull)
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized facet");}
    double dist=getSignedVerticalDistanceToZero(point);
    //[OBSOLETE CO20180828]if(zero_point_projection_only){return dist;}
    //sign of distance:
    //independent of lower/upper hull:  above hull is negative, below hull is positive
    dist-=point[point.urows];
    return dist;
  }

  vector<uint> ChullFacet::getCHIndices() const {
    if(!m_vertices.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No vertices found");}
    vector<uint> vi;
    for(uint i=0,fl_size_i=m_vertices.size();i<fl_size_i;i++){vi.push_back(m_vertices[i].ch_index);}
    return vi;
  }

  //[CO20200508 - OBSOLETE]void ChullFacet::create(ostream& oss) { //this it NOT an initialization, as we do this piece by piece
  //[CO20200508 - OBSOLETE]  xStream::free();
  //[CO20200508 - OBSOLETE]  ofstream* _p_FileMESSAGE=new ofstream();f_new_ofstream=true;
  //[CO20200508 - OBSOLETE]  create(*_p_FileMESSAGE,oss);
  //[CO20200508 - OBSOLETE]  f_new_ofstream=true;  //override
  //[CO20200508 - OBSOLETE]}

  //[CO20200508 - OBSOLETE]void ChullFacet::create(ofstream& FileMESSAGE,ostream& oss) { //this it NOT an initialization, as we do this piece by piece
  //[CO20200508 - OBSOLETE]  xStream::free();
  //[CO20200508 - OBSOLETE]  setOFStream(FileMESSAGE); f_new_ofstream=false;
  //[CO20200508 - OBSOLETE]  setOSS(oss);
  //[CO20200508 - OBSOLETE]}

  //MOVED TO xStream
  //void ChullFacet::setOFStream(ofstream& FileMESSAGE){p_FileMESSAGE=&FileMESSAGE;}
  //void ChullFacet::setOSS(ostream& oss) {p_oss=&oss;}

  void ChullFacet::addVertex(const ChullPointLight& point,uint index) {FacetPoint fp(point,index);return addVertex(fp);} //no need for full copy
  void ChullFacet::addVertex(const FacetPoint& fp){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    stringstream message;
    if(!fp.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized facetpoint");}
    const ChullPointLight& point=fp.ch_point;
    if(m_vertices.size()==0){
      m_has_stoich_coords=point.m_has_stoich_coords;
      m_formation_energy_coord=point.m_formation_energy_coord;
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " setting has_stoich_coords=" << m_has_stoich_coords << endl;
        cerr << __AFLOW_FUNC__ << " setting formation_energy_coord=" << m_formation_energy_coord << endl;
      }
    }
    else{
      //if(m_has_stoich_coords && !point.m_has_stoich_coords)
      if(m_has_stoich_coords != point.m_has_stoich_coords) //spit warning for either mismatch
      { //CO20200106 - patching for auto-indenting
        message << "Mismatch among coord types (stoich vs. non-stoich coords), assuming non-stoich coords";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        m_has_stoich_coords=false; //(m_has_stoich_coords && point.m_has_stoich_coords);
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " facet.m_has_stoich_coords=" << m_has_stoich_coords << endl;
          for(uint i=0,fl_size_i=m_vertices.size();i<fl_size_i;i++){
            cerr << __AFLOW_FUNC__ << " m_vertices[" << i << "].ch_point.m_has_stoich_coords=" << m_vertices[i].ch_point.m_has_stoich_coords;
            cerr << ", h_coords=" << m_vertices[i].ch_point.h_coords << endl;
          }
          cerr << __AFLOW_FUNC__ << " point.m_has_stoich_coords=" << point.m_has_stoich_coords << endl;
        }
      }
      //if(m_formation_energy_coord && !point.m_formation_energy_coord)
      if(m_formation_energy_coord != point.m_formation_energy_coord) //spit warning for either mismatch
      { //CO20200106 - patching for auto-indenting
        message << "Mismatch among coord types (formation_energy vs. non-formation_energy), assuming non-formation_energy coords";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        m_formation_energy_coord=false; //(m_formation_energy_coord && point.m_formation_energy_coord);
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " facet.m_formation_energy_coord=" << m_formation_energy_coord << endl;
          for(uint i=0,fl_size_i=m_vertices.size();i<fl_size_i;i++){
            cerr << __AFLOW_FUNC__ << " m_vertices[" << i << "].ch_point.m_formation_energy_coord=" << m_vertices[i].ch_point.m_formation_energy_coord;
            cerr << ", h_coords=" << m_vertices[i].ch_point.h_coords << endl;
          }
          cerr << __AFLOW_FUNC__ << " point.m_formation_energy_coord=" << point.m_formation_energy_coord << endl;
        }
      }
    }
    m_vertices.push_back(fp);
  }

  void ChullFacet::initialize(const xvector<double>& ref,uint h_dim,bool check_validity){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " initialize facet with points: " << endl;
      for(uint i=0,fl_size_i=m_vertices.size();i<fl_size_i;i++){cerr << "    " << m_vertices[i].ch_point.h_coords << endl;}
    }
    m_dim=h_dim;
    m_hull_reference=ref;
    setDirectiveVectors(check_validity);
    setNormal(check_validity);
    setOffset();  //FIX after normal alignment
    setCentroid();
    alignNormalInward(); //hull centroid
    setVertical();
    setArtificial();
    setHemisphere();
    setRidges();
    m_initialized=true;
  }

  bool ChullFacet::hasValidPoints(string& error){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    error.clear();
    if(!m_vertices.size()){error="Facet has no defining points";return false;}
    for(uint i=0,fl_size_i=m_vertices.size();i<fl_size_i;i++){if(!m_vertices[i].m_initialized){error="Uninitialized facetpoint";}}
    if((uint)m_vertices[0].ch_point.h_coords.rows!=m_vertices.size()){
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " m_vertices[0].ch_point.h_coords.rows=" << m_vertices[0].ch_point.h_coords.rows;
        cerr << " vs. m_vertices.size()=" << m_vertices.size() << endl;
      }
      error="Dimension mismatch among points and coordinates";return false;
    }
    for(uint i=1,fl_size_i=m_vertices.size();i<fl_size_i;i++){
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " " << m_vertices[i].ch_index << "  h_coords=" << m_vertices[i].ch_point.h_coords << endl;
        cerr << __AFLOW_FUNC__ << " m_vertices[i].ch_point.h_coords.rows=" << m_vertices[i].ch_point.h_coords.rows << " vs. m_vertices[0].ch_point.h_coords.rows=" << m_vertices[0].ch_point.h_coords.rows << endl;
      }
      if(m_vertices[i].ch_point.h_coords.rows!=m_vertices[0].ch_point.h_coords.rows){error="Dimension mismatch among facet points";return false;}
    }
    for(uint i=0,fl_size_i=m_vertices.size();i<fl_size_i;i++) {
      for(uint j=i+1,fl_size_j=m_vertices.size();j<fl_size_j;j++) {
        if(m_vertices[i].ch_index==m_vertices[j].ch_index){error="Facet points are degenerate";return false;}
      }
    }
    if(m_vertices.size()!=(uint)m_vertices[0].ch_point.h_coords.rows){error="Facet has wrong number of defining points given dimension";return false;}
    return true;
  }

  void ChullFacet::setContent(){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    m_content=AUROSTD_MAX_DOUBLE;
    m_is_hypercollinear=true;
    string error; if(!hasValidPoints(error)){return;}
    m_content=0.0;  //if we don't set later, it's because it's really zero
    xmatrix<double> B(m_vertices.size(),m_vertices.size(),1,1); //for determinant
    for(uint i=0,fl_size_i=m_vertices.size();i<fl_size_i;i++){
      for(uint k=0,fl_size_k=m_vertices.size();k<fl_size_k;k++){
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << "B(" << i << "+1," << k << "+1)=";
          cerr << "aurostd::modulussquare([" << m_vertices[i].ch_point.h_coords << "] - [" << m_vertices[k].ch_point.h_coords << "])" << endl;
        }
        B(i+1,k+1)=aurostd::modulussquare(m_vertices[i].ch_point.h_coords-m_vertices[k].ch_point.h_coords);
      }
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " B=" << endl; cerr << B << endl;}
    double CMdetB=aurostd::CMdet(B);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_has_stoich_coords=" << std::boolalpha << m_has_stoich_coords << ", m_formation_energy_coord=" << std::boolalpha << m_formation_energy_coord << endl;}
    //tolerance is CRITICAL, because we could accidentally skip a facet and get the wrong decomposition coefficients (negative sign)
    //tolerance should scale with dimension, i.e., if we are measuring length of line, use TOL, area of triangle, use TOL^2, volume of tetrahedron, use TOL^3
    //it is not clear to me what TOL should be, but I believe it is related to ZERO_COEF_TOL
    //heuristically, ZERO_COEF_TOL works well for lengths and areas (chull paper)
    //since ZERO_COEF_TOL^2 is ZERO_TOL, make this a simple if/else
    //double tol=(m_has_stoich_coords&&m_formation_energy_coord ? ZERO_COEF_TOL : ZERO_TOL);  //1e-4 if stoich coords
    double tol=ZERO_TOL;
    if((m_dim<4) && m_has_stoich_coords && m_formation_energy_coord){tol=ZERO_COEF_TOL;}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " CMdet(B)=" << CMdetB << endl;}
    //the coef will ONLY decrease the content (fraction)
    //if we don't check CMdetB before multiplying the coef, we might get -nan
    //therefore, use tol for both CHdetB and content
    if(aurostd::zeroWithinTol(CMdetB,tol)){return;} //error="CMdet(B) is zero, shows hyper-collinearity";
    double j=m_vertices.size()-1; //two vertices == line segment == 1-simplex
    double coef=pow(-1.0,j+1.0)/(pow(2.0,j)*pow(aurostd::factorial(j),2));
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " j=" << j << endl;
      cerr << __AFLOW_FUNC__ << " coef=" << coef << endl;
    }
    m_content=sqrt(coef*CMdetB);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " content=" << m_content << endl;}
    if(aurostd::zeroWithinTol(m_content,tol)){return;} //error="simplex content is zero, shows hyper-collinearity";
    m_is_hypercollinear=false;
  }

  void ChullFacet::setDirectiveVectors(bool check_validity){  //perhaps we already checked...
    m_directive_vectors.clear();
    string error;
    if(check_validity && !hasValidPoints(error)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,error);}

    std::sort(m_vertices.begin(),m_vertices.end());
    for(uint i=1,fl_size_i=m_vertices.size();i<fl_size_i;i++){
      m_directive_vectors.push_back( m_vertices[i].ch_point.h_coords - m_vertices[0].ch_point.h_coords );
    }
  }

  bool ChullFacet::pointsMatchDirectiveVectors(string& error){
    error.clear();
    //now we can define 
    if(m_directive_vectors.size()!=m_vertices.size()-1){error="Dimension mismatch among directive vectors and points";return false;}
    if((uint)m_directive_vectors[0].rows!=m_vertices.size()){error="Dimension mismatch among directive vectors coordinates and points";return false;}
    return true;
  }

  bool ChullFacet::hasValidDirectiveVectors(string& error){
    error.clear();
    if(!m_directive_vectors.size()){error="No directive vectors found";return false;}
    for(uint i=1,fl_size_i=m_directive_vectors.size();i<fl_size_i;i++){
      if(m_directive_vectors[i].rows!=m_directive_vectors[0].rows){error="Dimension mismatch among directive vectors coordinates";return false;}
    }
    for(uint i=0,fl_size_i=m_directive_vectors.size();i<fl_size_i;i++){if(aurostd::zeroWithinTol(modulus(m_directive_vectors[i]),ZERO_TOL)){error="Ill-defined directive vectors";return false;}}
    return true;
  }

  bool ChullFacet::hasCollinearVectors(bool check_validity){  //perhaps we already checked...
    string error;
    if(check_validity && !pointsMatchDirectiveVectors(error)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,error);}
    if(check_validity && !hasValidDirectiveVectors(error)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,error);}
    for(uint i=0,fl_size_i=m_directive_vectors.size();i<(fl_size_i-1);i++) {
      for(uint j=i+1,fl_size_j=m_directive_vectors.size();j<fl_size_j;j++) {
        if(aurostd::isCollinear(m_directive_vectors[i],m_directive_vectors[j],ZERO_TOL)){return true;}
      }
    }
    return false;
  }

  bool ChullFacet::isValid(string& error) {
    //technically, facets that show hyper-collinearity, i.e., three points on a line, four points on a plane, etc., are not really facets
    //we detect these specifically by checking collinearity in directive_vectors (2D) or, more generally in N dimensions, by checking
    //content (volume) of simplex
    //but by removing these facets, we screw up the neighboring determination, and hence create gaps in the hull
    //so keep them
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    error.clear();
    if(!hasValidPoints(error)){return false;}
    setContent(); //we need to preserve neighbors, so do not kill these pseudo-facets
    setDirectiveVectors(false);
    if(!pointsMatchDirectiveVectors(error)){return false;}
    if(!hasValidDirectiveVectors(error)){return false;} //we need to preserve neighbors, so do not kill these pseudo-facets
    if(0&&hasCollinearVectors(false)){  //we need to preserve neighbors, so do not kill these pseudo-facets
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " shows collinarity" << endl;}
      error="Directive vectors are collinear";return false;
    }
    return true;
  }

  void ChullFacet::setNormal(bool check_validity){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    m_normal.clear();
    string error;
    if(check_validity && !isValid(error)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,error);}
    if(LDEBUG){
      for(uint i=0,fl_size_i=m_directive_vectors.size();i<fl_size_i;i++){
        cerr << __AFLOW_FUNC__ << " directive_vector[" << i << "]=" << m_directive_vectors[i] << endl;
      }
    }
    if(!m_directive_vectors.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No directive vectors calculated");}
    m_normal=aurostd::getGeneralNormal(m_directive_vectors);
    if(aurostd::zeroWithinTol(aurostd::modulus(m_normal),ZERO_TOL)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid normal calculated");}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " normal=" << m_normal << endl;}
  }

  void ChullFacet::setOffset(){
    if(!m_vertices.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Facet has not been defined");}
    const xvector<double>& plane_point=m_vertices[0].ch_point.h_coords;
    if(m_normal.rows!=plane_point.rows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Dimension mismatch between normal and point");}
    m_offset=-scalar_product(m_normal,plane_point);
  }

  void ChullFacet::setCentroid() {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    vector<xvector<double> > points;
    for(uint i=0,fl_size_i=m_vertices.size();i<fl_size_i;i++){points.push_back(m_vertices[i].ch_point.h_coords);}
    m_facet_centroid=aurostd::getCentroid(points);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " centroid: " << m_facet_centroid << endl;}
  }

  //_AFLOW_CHULL_VERTICAL_PLANE_TOLERANCE_ = 1e-4 for H_f_atom, else 1e-9 for T_S
  void ChullFacet::setVertical(){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking if vertical hull: normal=" << m_normal << endl;}
    double tol=(m_has_stoich_coords&&m_formation_energy_coord ? ZERO_COEF_TOL : ZERO_TOL);  //1e-4 if stoich coords
    m_is_vertical=aurostd::zeroWithinTol(m_normal[m_normal.urows],tol); //simple
  }

  void ChullFacet::setArtificial(){  //in half hulls, this finds the facet of all artificial points
    for(uint i=0,fl_size_i=m_vertices.size();i<fl_size_i;i++){if(!m_vertices[i].ch_point.m_is_artificial){m_is_artificial=false;return;}}
    m_is_artificial=true;return;
  }

  void ChullFacet::alignNormalInward() { 
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    //we want normal pointing inward
    double dist=scalar_product(m_normal,m_hull_reference)+m_offset;  //point-plane distance formula without normalization
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " m_normal=" << m_normal << endl;
      cerr << __AFLOW_FUNC__ << " h_reference=" << m_hull_reference << endl;
      cerr << __AFLOW_FUNC__ << " dist=" << dist << endl;
    }
    bool negate=std::signbit(dist);
    if(negate){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " NEGATING" << endl;}
      m_normal=-m_normal; 
      m_offset=-m_offset; //flips sign of offset
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " aligned normal=" << m_normal << endl;}
  }

  void ChullFacet::setHemisphere() {
    //vertical facets are NOT considered lower_hemisphere
    m_in_lower_hemisphere=(!m_is_vertical && !std::signbit(m_normal[m_normal.urows])); //not flat and upward pointed normal
  }

  void ChullFacet::setFurthestPoint(){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    f_furthest_point.clear();
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Facet not initialized");}
    if(!f_outside_set.size()){return;}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " getting furthest point for facet with normal " << m_normal << " (is_artificial=" << m_is_artificial << ")" << endl;}
    double dist=AUROSTD_MAX_DOUBLE,max_dist=0;
    for(uint i=0,fl_size_i=f_outside_set.size();i<fl_size_i;i++){
      dist=abs(getSignedPointPlaneDistance(f_outside_set[i].ch_point.h_coords)); //we only care about magnitude here
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " point=" << f_outside_set[i].ch_point.h_coords << " is " << dist << " from facet" << endl;}
      if(dist>max_dist){
        max_dist=dist;
        f_furthest_point=f_outside_set[i];
      }
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " furthest point from facet:" << endl;
      for(uint i=0,fl_size_i=m_vertices.size();i<fl_size_i;i++){cerr << "       point[" << m_vertices[i].ch_index << "]=" << m_vertices[i].ch_point.h_coords << endl;}
      cerr << "    is point[" << f_furthest_point.ch_index << "]=" << f_furthest_point.ch_point.h_coords << endl;
    }
  }

  void ChullFacet::setRidges(){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    m_ridges.clear();
    if(!m_vertices.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Facet has no vertices");}
    for(uint fl_size_i=m_vertices.size(),i=(fl_size_i-1);i<fl_size_i;i--){  //wrap around to HUGE number
      m_ridges.push_back(ChullFacet(*p_FileMESSAGE,*p_oss));  //CO20180305
      for(uint j=0,fl_size_j=m_vertices.size();j<fl_size_j;j++){
        if(i!=j){m_ridges.back().addVertex(m_vertices[j]);}
      }
      if(m_ridges.back().m_vertices.size()!=m_dim-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Ridge vertex count and facet dimension mismatch");}
      std::sort(m_ridges.back().m_vertices.begin(),m_ridges.back().m_vertices.end());
    }
    std::sort(m_ridges.begin(),m_ridges.end());

    if(LDEBUG){
      for(uint i=0,fl_size_i=m_ridges.size();i<fl_size_i;i++){
        cerr << __AFLOW_FUNC__ << " ridges[" << i << "]: ";
        for(uint j=0,fl_size_j=m_ridges[i].m_vertices.size();j<fl_size_j;j++){
          cerr << m_ridges[i].m_vertices[j].ch_index << " " ;
        }
        cerr << endl;
      }
    }
  }

  void ChullFacet::cleanFacet() {
    f_visited=false;
    f_outside_set.clear();
    f_furthest_point.clear();
    f_neighbors.clear();
  }
} // namespace chull

namespace chull {
  //--------------------------------------------------------------------------------
  // constructor
  //--------------------------------------------------------------------------------
  CoordGroup::CoordGroup() {free();}
  CoordGroup::CoordGroup(const xvector<double>& coord,bool has_stoich_coords) {initialize(coord,has_stoich_coords);}
  CoordGroup::CoordGroup(const CoordGroup& b) {copy(b);}  // copy PUBLIC
  CoordGroup::~CoordGroup() {free();}

  const CoordGroup& CoordGroup::operator=(const CoordGroup& other) {
    if(this!=&other) {copy(other);}
    return *this;
  }

  bool CoordGroup::operator<(const CoordGroup& other) const {
    // safety, so it doesn't break, but it's outside scope of function
    if(m_coords.rows!=other.m_coords.rows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Dimension mismatch between stoichiometries");} //return (m_coords.rows<other.m_coords.rows);
    for(int i=m_coords.lrows;i<=m_coords.urows;i++) {
      if(m_coords(i)!=other.m_coords(i)) {return (m_coords(i)<other.m_coords(i));}
    }
    return false;
  }

  void CoordGroup::clear() {CoordGroup a;copy(a);} //clear PUBLIC
  void CoordGroup::free() {
    m_initialized=false;
    m_coords.clear();
    m_points.clear();
    m_has_stoich_coords=false;
    m_has_artificial_unary=false;
    m_is_on_hull=false;
    m_hull_member=AUROSTD_MAX_UINT;
    m_ref_state=AUROSTD_MAX_UINT;
    m_candidate_hull_points.clear();
    m_i_nary=AUROSTD_MAX_UINT;
    m_i_alloy=AUROSTD_MAX_UINT;
    m_nearest_facet=AUROSTD_MAX_UINT;
    m_nearest_distance=AUROSTD_MAX_DOUBLE;
    m_decomp_phases.clear();
    m_decomp_coefs.clear();
    m_equilibrium_phases.clear();
    m_calculated_equivalent_g_states=false;
    m_equivalent_g_states.clear();
    m_sym_equivalent_g_states.clear();
    m_stability_criterion=AUROSTD_MAX_DOUBLE;
    m_n_plus_1_enthalpy_gain=AUROSTD_MAX_DOUBLE;
    m_found_icsd_g_state=false;
    m_i_icsd_g_state=AUROSTD_MAX_UINT;
  }

  void CoordGroup::copy(const CoordGroup& b) {
    m_initialized=b.m_initialized;
    m_coords=b.m_coords;
    m_points.clear(); for(uint i=0,fl_size_i=b.m_points.size();i<fl_size_i;i++){m_points.push_back(b.m_points[i]);}
    m_has_stoich_coords=b.m_has_stoich_coords;
    m_has_artificial_unary=b.m_has_artificial_unary;
    m_is_on_hull=b.m_is_on_hull;
    m_hull_member=b.m_hull_member;
    m_ref_state=b.m_ref_state;
    m_candidate_hull_points.clear(); for(uint i=0,fl_size_i=b.m_candidate_hull_points.size();i<fl_size_i;i++){m_candidate_hull_points.push_back(b.m_candidate_hull_points[i]);}
    m_i_nary=b.m_i_nary;
    m_i_alloy=b.m_i_alloy;
    m_nearest_facet=b.m_nearest_facet;
    m_nearest_distance=b.m_nearest_distance;
    m_decomp_phases.clear(); for(uint i=0,fl_size_i=b.m_decomp_phases.size();i<fl_size_i;i++){m_decomp_phases.push_back(b.m_decomp_phases[i]);}
    m_decomp_coefs=b.m_decomp_coefs;
    for(uint i=0,fl_size_i=m_equilibrium_phases.size();i<fl_size_i;i++){m_equilibrium_phases.clear();} m_equilibrium_phases.clear(); for(uint i=0;i<b.m_equilibrium_phases.size();i++){m_equilibrium_phases.push_back(b.m_equilibrium_phases[i]);}
    m_calculated_equivalent_g_states=b.m_calculated_equivalent_g_states;
    m_equivalent_g_states.clear(); for(uint i=0,fl_size_i=b.m_equivalent_g_states.size();i<fl_size_i;i++){m_equivalent_g_states.push_back(b.m_equivalent_g_states[i]);}
    m_sym_equivalent_g_states.clear(); for(uint i=0,fl_size_i=b.m_sym_equivalent_g_states.size();i<fl_size_i;i++){m_sym_equivalent_g_states.push_back(b.m_sym_equivalent_g_states[i]);}
    m_stability_criterion=b.m_stability_criterion;
    m_n_plus_1_enthalpy_gain=b.m_n_plus_1_enthalpy_gain;
    m_found_icsd_g_state=b.m_found_icsd_g_state;
    m_i_icsd_g_state=b.m_i_icsd_g_state;
  }

  void CoordGroup::initialize(const xvector<double>& coord,bool has_stoich_coords) {
    free();
    m_coords=coord;
    m_has_stoich_coords=has_stoich_coords;
    m_has_artificial_unary=false;
    m_initialized=true;
  }

  xvector<int> CoordGroup::getElementsPresent() const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_coords=" << m_coords << endl;}
    for(int i=m_coords.lrows;i<=m_coords.urows;i++){
      if(std::signbit(m_coords[i]) || m_coords[i]>1.0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coord("+aurostd::utype2string(i)+") is outside of [0,1] range of a generalized stoichiometry coordinate");}
    }
    double hid_dim=1.0-sum(m_coords);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " hid_dim=" << hid_dim << endl;}
    if(std::signbit(hid_dim) || hid_dim>1.0) {
      if(aurostd::zeroWithinTol(hid_dim,ZERO_TOL)){hid_dim=0.0;}  //only zero out for the last coord, as it's derived from the subtraction of the others
      else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coord("+aurostd::utype2string(m_coords.rows)+") is outside of [0,1] range of a generalized stoichiometry coordinate");}
    }

    xvector<int> elements_present(m_coords.lrows,m_coords.urows+1);
    for(int i=m_coords.lrows;i<=m_coords.urows;i++){if(aurostd::nonZeroWithinTol(m_coords[i],ZERO_TOL)){elements_present[i]=1;}}
    if(aurostd::nonZeroWithinTol(hid_dim,ZERO_TOL)){elements_present[elements_present.urows]=1;}

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " elements_present=" << elements_present << endl;}
    return elements_present;
  }

  uint CoordGroup::getDim() const {
    if(!m_initialized){
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");
    }
    return sum(getElementsPresent());
  }

} // namespace chull

namespace chull {
  //--------------------------------------------------------------------------------
  // constructor
  //--------------------------------------------------------------------------------
  Alloy::Alloy() {free();}
  Alloy::Alloy(const xvector<int>& elements_present) {initialize(elements_present);}
  Alloy::Alloy(const Alloy& b) {copy(b);}
  Alloy::~Alloy() {free();}

  const Alloy& Alloy::operator=(const Alloy& other) {
    if(this!=&other) {copy(other);}
    return *this;
  }

  bool Alloy::operator<(const Alloy& other) const {
    // safety, so it doesn't break, but it's outside scope of function
    if(m_elements_present.rows!=other.m_elements_present.rows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Dimension mismatch between alloys");} //return (m_elements_present.rows<other.m_elements_present.rows);
    for(int i=m_elements_present.lrows;i<=m_elements_present.urows;i++) {
      if(m_elements_present[i]!=other.m_elements_present[i]) {return (m_elements_present[i]<other.m_elements_present[i]);}
    }
    return false;
  }

  void Alloy::clear() {Alloy a;copy(a);} //clear PUBLIC
  void Alloy::free() { 
    m_initialized=false;
    m_elements_present.clear();
    m_dim=AUROSTD_MAX_UINT;
    m_coord_groups.clear();
    m_facets.clear();
  }

  void Alloy::copy(const Alloy& b) {
    m_initialized=b.m_initialized;
    m_elements_present=b.m_elements_present;
    m_dim=b.m_dim;
    m_coord_groups.clear(); for(uint i=0,fl_size_i=b.m_coord_groups.size();i<fl_size_i;i++){m_coord_groups.push_back(b.m_coord_groups[i]);}
    m_facets.clear(); for(uint i=0,fl_size_i=b.m_facets.size();i<fl_size_i;i++){m_facets.push_back(b.m_facets[i]);}
  }

  void Alloy::initialize(const xvector<int>& elements_present){
    free();
    m_elements_present=elements_present;
    m_dim=sum(elements_present);
    m_initialized=true;
  }

  bool Alloy::belongs2Hull(const xvector<int>& elements_present_hull) const {return subspaceBelongs(elements_present_hull,m_elements_present);}
} // namespace chull

namespace chull {
  //--------------------------------------------------------------------------------
  // constructor
  //--------------------------------------------------------------------------------
  Nary::Nary() {free();}
  Nary::Nary(uint dim) {initialize(dim);}
  Nary::Nary(const Nary& b) {copy(b);}
  Nary::~Nary() {free();}

  const Nary& Nary::operator=(const Nary& other) {
    if(this!=&other) {copy(other);}
    return *this;
  }

  bool Nary::operator<(const Nary& other) const {return (nary<other.nary);}

  void Nary::clear() {Nary a;copy(a);}
  void Nary::free() {
    m_initialized=false;
    nary=AUROSTD_MAX_UINT;
    m_alloys.clear();
  }

  void Nary::copy(const Nary& b) {
    m_initialized=b.m_initialized;
    nary=b.nary;
    m_alloys.clear(); for(uint i=0,fl_size_i=b.m_alloys.size();i<fl_size_i;i++){m_alloys.push_back(b.m_alloys[i]);}
  }

  void Nary::initialize(uint dim) {
    free();
    nary=dim;
    m_initialized=true;
  }
} // namespace chull

namespace chull {
  //--------------------------------------------------------------------------------
  // class ConvexHull
  //--------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------
  // constructor
  //--------------------------------------------------------------------------------
  ConvexHull::ConvexHull(ostream& oss) : xStream(oss),m_initialized(false) {initialize();}
  ConvexHull::ConvexHull(const string& alloy,ostream& oss) : xStream(oss),m_initialized(false) {initialize(alloy);}
  ConvexHull::ConvexHull(const vector<string>& velements,ostream& oss) : xStream(oss),m_initialized(false) {initialize(velements);}
  ConvexHull::ConvexHull(const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries,ostream& oss) : xStream(oss),m_initialized(false) {initialize(velements,entries);}
  ConvexHull::ConvexHull(const vector<xvector<double> >& vcoords,ostream& oss,bool has_stoich_coords,bool formation_energy_hull,bool add_artificial_unaries) : xStream(oss),m_initialized(false) {initialize(vcoords,has_stoich_coords,formation_energy_hull,add_artificial_unaries);}
  ConvexHull::ConvexHull(const vector<ChullPoint>& vpoints,ostream& oss,bool formation_energy_hull,bool add_artificial_unaries) : xStream(oss),m_initialized(false) {initialize(vpoints,formation_energy_hull,add_artificial_unaries);}
  ConvexHull::ConvexHull(const vector<ChullPoint>& vpoints,const vector<string>& velements,ostream& oss,bool formation_energy_hull,bool add_artificial_unaries) : xStream(oss),m_initialized(false) {initialize(vpoints,velements,formation_energy_hull,add_artificial_unaries);}
  ConvexHull::ConvexHull(ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize();}
  ConvexHull::ConvexHull(const string& alloy,ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize(alloy);}
  ConvexHull::ConvexHull(const vector<string>& velements,ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize(velements);}
  ConvexHull::ConvexHull(const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries,ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize(velements,entries);}
  ConvexHull::ConvexHull(const vector<xvector<double> >& vcoords,ofstream& FileMESSAGE,ostream& oss,bool has_stoich_coords,bool formation_energy_hull,bool add_artificial_unaries) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize(vcoords,has_stoich_coords,formation_energy_hull,add_artificial_unaries);}
  ConvexHull::ConvexHull(const vector<ChullPoint>& vpoints,ofstream& FileMESSAGE,ostream& oss,bool formation_energy_hull,bool add_artificial_unaries) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize(vpoints,formation_energy_hull,add_artificial_unaries);}
  ConvexHull::ConvexHull(const vector<ChullPoint>& vpoints,const vector<string>& velements,ofstream& FileMESSAGE,ostream& oss,bool formation_energy_hull,bool add_artificial_unaries) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize(vpoints,velements,formation_energy_hull,add_artificial_unaries);}
  ConvexHull::ConvexHull(const aurostd::xoption& vpflow,ostream& oss) : xStream(oss),m_initialized(false) {initialize(vpflow);}
  ConvexHull::ConvexHull(const aurostd::xoption& vpflow,const string& alloy,ostream& oss) : xStream(oss),m_initialized(false) {initialize(vpflow,alloy);}
  ConvexHull::ConvexHull(const aurostd::xoption& vpflow,const vector<string>& velements,ostream& oss) : xStream(oss),m_initialized(false) {initialize(vpflow,velements);}
  ConvexHull::ConvexHull(const aurostd::xoption& vpflow,const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries,ostream& oss) : xStream(oss),m_initialized(false) {initialize(vpflow,velements,entries);}
  ConvexHull::ConvexHull(const aurostd::xoption& vpflow,const vector<xvector<double> >& vcoords,ostream& oss,bool has_stoich_coords,bool formation_energy_hull,bool add_artificial_unaries) : xStream(oss),m_initialized(false) {initialize(vpflow,vcoords,has_stoich_coords,formation_energy_hull,add_artificial_unaries);}
  ConvexHull::ConvexHull(const aurostd::xoption& vpflow,const vector<ChullPoint>& vpoints,ostream& oss,bool formation_energy_hull,bool add_artificial_unaries) : xStream(oss),m_initialized(false) {initialize(vpflow,vpoints,formation_energy_hull,add_artificial_unaries);}
  ConvexHull::ConvexHull(const aurostd::xoption& vpflow,const vector<ChullPoint>& vpoints,const vector<string>& velements,ostream& oss,bool formation_energy_hull,bool add_artificial_unaries) : xStream(oss),m_initialized(false) {initialize(vpflow,vpoints,velements,formation_energy_hull,add_artificial_unaries);}
  ConvexHull::ConvexHull(const aurostd::xoption& vpflow,ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize(vpflow);}
  ConvexHull::ConvexHull(const aurostd::xoption& vpflow,const string& alloy,ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize(vpflow,alloy);}
  ConvexHull::ConvexHull(const aurostd::xoption& vpflow,const vector<string>& velements,ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize(vpflow,velements);}
  ConvexHull::ConvexHull(const aurostd::xoption& vpflow,const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries,ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize(vpflow,velements,entries);}
  ConvexHull::ConvexHull(const aurostd::xoption& vpflow,const vector<xvector<double> >& vcoords,ofstream& FileMESSAGE,ostream& oss,bool has_stoich_coords,bool formation_energy_hull,bool add_artificial_unaries) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize(vpflow,vcoords,has_stoich_coords,formation_energy_hull,add_artificial_unaries);}
  ConvexHull::ConvexHull(const aurostd::xoption& vpflow,const vector<ChullPoint>& vpoints,ofstream& FileMESSAGE,ostream& oss,bool formation_energy_hull,bool add_artificial_unaries) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize(vpflow,vpoints,formation_energy_hull,add_artificial_unaries);}
  ConvexHull::ConvexHull(const aurostd::xoption& vpflow,const vector<ChullPoint>& vpoints,const vector<string>& velements,ofstream& FileMESSAGE,ostream& oss,bool formation_energy_hull,bool add_artificial_unaries) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize(vpflow,vpoints,velements,formation_energy_hull,add_artificial_unaries);}
  ConvexHull::ConvexHull(const ConvexHull& b) : xStream(*b.getOFStream(),*b.getOSS()) {copy(b);}

  ConvexHull::~ConvexHull() {
    xStream::free();
    free();
    m_allowed_dft_types.clear();
  }

  const ConvexHull& ConvexHull::operator=(const ConvexHull& other) {
    if(this!=&other) {copy(other);}
    return *this;
  }

  void ConvexHull::clear() {ConvexHull a;copy(a);}  //clear PRIVATE
  void ConvexHull::free() {
    m_initialized=false;
    m_velements.clear();
    m_icsd_entries.clear();
    m_points.clear();
    m_naries.clear();
    m_coord_groups.clear();
    m_dim=AUROSTD_MAX_UINT;
    m_half_hull=false;
    m_lower_hull=false;
    m_has_stoich_coords=false;
    m_add_artificial_unaries=false;
    m_thermo_hull=false;
    m_formation_energy_hull=false;
    m_facets.clear();
    m_i_facets.clear();
    m_sort_energy_ascending=true;
    m_cflags.clear();
    m_aflags.clear();
    //see xStream::free()
    //p_oss=NULL;
    //p_FileMESSAGE=NULL;
    //f_new_ofstream=false;
    m_allow_all_formation_energies=DEFAULT_CHULL_ALLOW_ALL_FORMATION_ENERGIES;
    aurostd::string2tokens(DEFAULT_CHULL_ALLOWED_DFT_TYPES,m_allowed_dft_types,",");
    if(!m_allowed_dft_types.size()){m_allowed_dft_types.push_back("PAW_PBE");}
    cleanHull();
  }

  void ConvexHull::copy(const ConvexHull& b) {  //copy PRIVATE
    xStream::copy(b);
    m_initialized=b.m_initialized;
    m_velements.clear(); for(uint i=0,fl_size_i=b.m_velements.size();i<fl_size_i;i++){m_velements.push_back(b.m_velements[i]);}
    m_icsd_entries.clear(); for(uint i=0,fl_size_i=b.m_icsd_entries.size();i<fl_size_i;i++){m_icsd_entries.push_back(b.m_icsd_entries[i]);}
    m_points.clear(); for(uint i=0,fl_size_i=b.m_points.size();i<fl_size_i;i++){m_points.push_back(b.m_points[i]);}
    m_naries.clear(); for(uint i=0,fl_size_i=b.m_naries.size();i<fl_size_i;i++){m_naries.push_back(b.m_naries[i]);}
    m_coord_groups.clear(); for(uint i=0,fl_size_i=b.m_coord_groups.size();i<fl_size_i;i++){m_coord_groups.push_back(b.m_coord_groups[i]);}
    m_dim=b.m_dim;
    m_half_hull=b.m_half_hull;
    m_lower_hull=b.m_lower_hull;
    m_has_stoich_coords=b.m_has_stoich_coords;
    m_add_artificial_unaries=b.m_add_artificial_unaries;
    m_thermo_hull=b.m_thermo_hull;
    m_formation_energy_hull=b.m_formation_energy_hull;
    m_facets.clear(); for(uint i=0,fl_size_i=b.m_facets.size();i<fl_size_i;i++){m_facets.push_back(b.m_facets[i]);}
    m_i_facets.clear(); for(uint i=0,fl_size_i=b.m_i_facets.size();i<fl_size_i;i++){m_i_facets.push_back(b.m_i_facets[i]);}
    m_sort_energy_ascending=b.m_sort_energy_ascending;
    m_cflags=b.m_cflags;
    m_aflags=b.m_aflags;
    m_allow_all_formation_energies=b.m_allow_all_formation_energies;
    m_allowed_dft_types.clear(); for(uint i=0,fl_size_i=b.m_allowed_dft_types.size();i<fl_size_i;i++){m_allowed_dft_types.push_back(b.m_allowed_dft_types[i]);}
    h_dim=b.h_dim;
    m_elements_present=b.m_elements_present;
    h_points.clear(); for(uint i=0,fl_size_i=b.h_points.size();i<fl_size_i;i++){h_points.push_back(b.h_points[i]);}
    h_centroid=b.h_centroid;
    h_reference=b.h_reference;
    h_facets.clear(); for(uint i=0,fl_size_i=b.h_facets.size();i<fl_size_i;i++){h_facets.push_back(b.h_facets[i]);}
    h_visible_facets.clear(); for(uint i=0,fl_size_i=b.h_visible_facets.size();i<fl_size_i;i++){h_visible_facets.push_back(b.h_visible_facets[i]);}
    h_horizon_ridges.clear(); for(uint i=0,fl_size_i=b.h_horizon_ridges.size();i<fl_size_i;i++){h_horizon_ridges.push_back(b.h_horizon_ridges[i]);}
  }

  bool ConvexHull::initialize(ostream& oss) {
    xStream::initialize(oss);
    return initialize();
  }

  bool ConvexHull::initialize(const string& alloy,ostream& oss) {
    xStream::initialize(oss);
    return initialize(alloy);
  }

  bool ConvexHull::initialize(const vector<string>& velements,ostream& oss) {
    xStream::initialize(oss);
    return initialize(velements);
  }

  bool ConvexHull::initialize(const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries,ostream& oss) {
    xStream::initialize(oss);
    return initialize(velements,entries);
  }

  bool ConvexHull::initialize(const vector<xvector<double> >& vcoords,ostream& oss,bool has_stoich_coords,bool formation_enthalpy_hull,bool add_artificial_unaries) {
    xStream::initialize(oss);
    return initialize(vcoords,has_stoich_coords,formation_enthalpy_hull,add_artificial_unaries);
  }

  bool ConvexHull::initialize(const vector<ChullPoint>& vpoints,ostream& oss,bool formation_enthalpy_hull,bool add_artificial_unaries) {
    xStream::initialize(oss);
    return initialize(vpoints,formation_enthalpy_hull,add_artificial_unaries);
  }

  bool ConvexHull::initialize(const vector<ChullPoint>& vpoints,const vector<string>& velements,ostream& oss,bool formation_enthalpy_hull,bool add_artificial_unaries) {
    xStream::initialize(oss);
    return initialize(vpoints,velements,formation_enthalpy_hull,add_artificial_unaries);
  }

  bool ConvexHull::initialize(ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize();
  }

  bool ConvexHull::initialize(const string& alloy,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(alloy);
  }

  bool ConvexHull::initialize(const vector<string>& velements,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(velements);
  }

  bool ConvexHull::initialize(const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(velements,entries);
  }

  bool ConvexHull::initialize(const vector<xvector<double> >& vcoords,ofstream& FileMESSAGE,ostream& oss,bool has_stoich_coords,bool formation_enthalpy_hull,bool add_artificial_unaries) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(vcoords,has_stoich_coords,formation_enthalpy_hull,add_artificial_unaries);
  }

  bool ConvexHull::initialize(const vector<ChullPoint>& vpoints,ofstream& FileMESSAGE,ostream& oss,bool formation_enthalpy_hull,bool add_artificial_unaries) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(vpoints,formation_enthalpy_hull,add_artificial_unaries);
  }

  bool ConvexHull::initialize(const vector<ChullPoint>& vpoints,const vector<string>& velements,ofstream& FileMESSAGE,ostream& oss,bool formation_enthalpy_hull,bool add_artificial_unaries) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(vpoints,velements,formation_enthalpy_hull,add_artificial_unaries);
  }

  bool ConvexHull::initialize() {
    free();
    try{
      setDefaultCFlags();
      setDirectory();
      m_initialized=false;  //no points
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }

  bool ConvexHull::initialize(const string& alloy) {
    free();
    try{
      setDefaultCFlags();
      setDirectory();
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return createHull(alloy);
  }

  bool ConvexHull::initialize(const vector<string>& velements) {
    free();
    try{
      setDefaultCFlags();
      setDirectory();
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return createHull(velements);
  }

  bool ConvexHull::initialize(const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries) {
    free();
    try{
      setDefaultCFlags();
      setDirectory();
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return createHull(velements,entries);
  }

  bool ConvexHull::initialize(const vector<xvector<double> >& vcoords,bool has_stoich_coords,bool formation_enthalpy_hull,bool add_artificial_unaries) {
    free();
    try{
      setDefaultCFlags();
      setDirectory();
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return createHull(vcoords,has_stoich_coords,formation_enthalpy_hull,add_artificial_unaries);
  }

  bool ConvexHull::initialize(const vector<ChullPoint>& vpoints,bool formation_enthalpy_hull,bool add_artificial_unaries) {
    free();
    try{
      setDefaultCFlags();
      setDirectory();
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return createHull(vpoints,formation_enthalpy_hull,add_artificial_unaries);
  }

  bool ConvexHull::initialize(const vector<ChullPoint>& vpoints,const vector<string>& velements,bool formation_enthalpy_hull,bool add_artificial_unaries) {
    free();
    try{
      setDefaultCFlags();
      setDirectory();
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return createHull(vpoints,velements,formation_enthalpy_hull,add_artificial_unaries);
  }

  bool ConvexHull::initialize(const aurostd::xoption& vpflow,ostream& oss) {
    xStream::initialize(oss);
    return initialize(vpflow);
  }

  bool ConvexHull::initialize(const aurostd::xoption& vpflow,const string& alloy,ostream& oss) {
    xStream::initialize(oss);
    return initialize(vpflow,alloy);
  }

  bool ConvexHull::initialize(const aurostd::xoption& vpflow,const vector<string>& velements,ostream& oss) {
    xStream::initialize(oss);
    return initialize(vpflow,velements);
  }

  bool ConvexHull::initialize(const aurostd::xoption& vpflow,const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries,ostream& oss) {
    xStream::initialize(oss);
    return initialize(vpflow,velements,entries);
  }

  bool ConvexHull::initialize(const aurostd::xoption& vpflow,const vector<xvector<double> >& vcoords,ostream& oss,bool has_stoich_coords,bool formation_enthalpy_hull,bool add_artificial_unaries) {
    xStream::initialize(oss);
    return initialize(vpflow,vcoords,has_stoich_coords,formation_enthalpy_hull,add_artificial_unaries);
  }

  bool ConvexHull::initialize(const aurostd::xoption& vpflow,const vector<ChullPoint>& vpoints,ostream& oss,bool formation_enthalpy_hull,bool add_artificial_unaries) {
    xStream::initialize(oss);
    return initialize(vpflow,vpoints,formation_enthalpy_hull,add_artificial_unaries);
  }

  bool ConvexHull::initialize(const aurostd::xoption& vpflow,const vector<ChullPoint>& vpoints,const vector<string>& velements,ostream& oss,bool formation_enthalpy_hull,bool add_artificial_unaries) {
    xStream::initialize(oss);
    return initialize(vpflow,vpoints,velements,formation_enthalpy_hull,add_artificial_unaries);
  }

  bool ConvexHull::initialize(const aurostd::xoption& vpflow,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(vpflow);
  }

  bool ConvexHull::initialize(const aurostd::xoption& vpflow,const string& alloy,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(vpflow,alloy);
  }

  bool ConvexHull::initialize(const aurostd::xoption& vpflow,const vector<string>& velements,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(vpflow,velements);
  }

  bool ConvexHull::initialize(const aurostd::xoption& vpflow,const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(vpflow,velements,entries);
  }

  bool ConvexHull::initialize(const aurostd::xoption& vpflow,const vector<xvector<double> >& vcoords,ofstream& FileMESSAGE,ostream& oss,bool has_stoich_coords,bool formation_enthalpy_hull,bool add_artificial_unaries) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(vpflow,vcoords,has_stoich_coords,formation_enthalpy_hull,add_artificial_unaries);
  }

  bool ConvexHull::initialize(const aurostd::xoption& vpflow,const vector<ChullPoint>& vpoints,ofstream& FileMESSAGE,ostream& oss,bool formation_enthalpy_hull,bool add_artificial_unaries) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(vpflow,vpoints,formation_enthalpy_hull,add_artificial_unaries);
  }

  bool ConvexHull::initialize(const aurostd::xoption& vpflow,const vector<ChullPoint>& vpoints,const vector<string>& velements,ofstream& FileMESSAGE,ostream& oss,bool formation_enthalpy_hull,bool add_artificial_unaries) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(vpflow,vpoints,velements,formation_enthalpy_hull,add_artificial_unaries);
  }

  bool ConvexHull::initialize(const aurostd::xoption& vpflow) {
    free();
    try{
      setCFlags(vpflow);
      setDirectory();
      m_initialized=false;  //no points
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }

  bool ConvexHull::initialize(const aurostd::xoption& vpflow,const string& alloy) {
    free();
    try{
      setCFlags(vpflow);
      setDirectory();
      m_initialized=createHull(alloy);
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }

  bool ConvexHull::initialize(const aurostd::xoption& vpflow,const vector<string>& velements) {
    free();
    try{
      setCFlags(vpflow);
      setDirectory();
      m_initialized=createHull(velements);
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }

  bool ConvexHull::initialize(const aurostd::xoption& vpflow,const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries) {
    free();
    try{
      setCFlags(vpflow);
      setDirectory();
      m_initialized=createHull(velements,entries);
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }

  bool ConvexHull::initialize(const aurostd::xoption& vpflow,const vector<xvector<double> >& vcoords,bool has_stoich_coords,bool formation_enthalpy_hull,bool add_artificial_unaries) {
    free();
    try{
      setCFlags(vpflow);
      setDirectory();
      m_initialized=createHull(vcoords,has_stoich_coords,formation_enthalpy_hull,add_artificial_unaries);
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }

  bool ConvexHull::initialize(const aurostd::xoption& vpflow,const vector<ChullPoint>& vpoints,bool formation_enthalpy_hull,bool add_artificial_unaries) {
    free();
    try{
      setCFlags(vpflow);
      setDirectory();
      m_initialized=createHull(vpoints,formation_enthalpy_hull,add_artificial_unaries);
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }

  bool ConvexHull::initialize(const aurostd::xoption& vpflow,const vector<ChullPoint>& vpoints,const vector<string>& velements,bool formation_enthalpy_hull,bool add_artificial_unaries) {
    free();
    try{
      setCFlags(vpflow);
      setDirectory();
      m_initialized=createHull(vpoints,velements,formation_enthalpy_hull,add_artificial_unaries);
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }

  void ConvexHull::initializePoints(const string& alloy){
    loadPoints(alloy);
    structurePoints();
  }

  void ConvexHull::initializePoints(const vector<string>& velements){
    loadPoints(velements);
    structurePoints();
  }

  void ConvexHull::initializePoints(const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries){
    loadPoints(velements,entries);
    structurePoints();
  }

  void ConvexHull::initializePoints(const vector<xvector<double> >& vcoords,bool has_stoich_coords,bool formation_enthalpy_hull,bool add_artificial_unaries){
    loadPoints(vcoords,has_stoich_coords,formation_enthalpy_hull,add_artificial_unaries);
    structurePoints();
  }

  void ConvexHull::initializePoints(const vector<ChullPoint>& vpoints,bool formation_enthalpy_hull,bool add_artificial_unaries){
    loadPoints(vpoints,formation_enthalpy_hull,add_artificial_unaries);
    structurePoints();
  }

  void ConvexHull::initializePoints(const vector<ChullPoint>& vpoints,const vector<string>& velements,bool formation_enthalpy_hull,bool add_artificial_unaries){
    loadPoints(vpoints,velements,formation_enthalpy_hull,add_artificial_unaries);
    structurePoints();
  }

  uint ConvexHull::getDim() const {return m_dim;}
  uint ConvexHull::getEntriesCount(bool only_within_half_hull) const {
    uint i_point=AUROSTD_MAX_UINT,count=0;
    for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
      if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup["+aurostd::utype2string(i_coord_group)+"] is not initialized");}
      for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
        i_point=m_coord_groups[i_coord_group].m_points[i];
        if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
        if(m_points[i_point].m_is_artificial){continue;}
        if(only_within_half_hull&&m_half_hull){
          if(m_points[i_point].isWithinHalfHull(m_lower_hull)){count++;}
        }else{count++;}
      }
    }
    return count;
  }

  uint ConvexHull::getEntriesCount(uint i_nary,bool only_within_half_hull) const {
    if(!m_has_stoich_coords){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Non-stoich coordinates");}
    uint count=0;
    if(i_nary>m_naries.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    for(uint i_alloy=0,fl_size_i_alloy=m_naries[i_nary].m_alloys.size();i_alloy<fl_size_i_alloy;i_alloy++){count+=getEntriesCount(i_nary,i_alloy,only_within_half_hull);}
    return count;
  }

  uint ConvexHull::getEntriesCount(uint i_nary,uint i_alloy,bool only_within_half_hull) const {
    if(!m_has_stoich_coords){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Non-stoich coordinates");}
    if(i_nary>m_naries.size()-1 || i_alloy>m_naries[i_nary].m_alloys.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}
    uint i_point=AUROSTD_MAX_UINT,count=0,i_coord_group=AUROSTD_MAX_UINT;
    for(uint i=0,fl_size_i=m_naries[i_nary].m_alloys[i_alloy].m_coord_groups.size();i<fl_size_i;i++){
      i_coord_group=m_naries[i_nary].m_alloys[i_alloy].m_coord_groups[i];
      if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup["+aurostd::utype2string(i_coord_group)+"] is not initialized");}
      for(uint j=0,fl_size_j=m_coord_groups[i_coord_group].m_points.size();j<fl_size_j;j++){
        i_point=m_coord_groups[i_coord_group].m_points[j];
        if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
        if(m_points[i_point].m_is_artificial){continue;}
        if(only_within_half_hull&&m_half_hull){
          if(m_points[i_point].isWithinHalfHull(m_lower_hull)){count++;}
        }else{count++;}
      }
    }
    return count;
  }

  vector<vector<uint> > ConvexHull::getHullSizes(bool only_within_half_hull) const {
    vector<vector<uint> > counts;
    for(uint i_nary1=0,fl_size_i_nary1=m_naries.size();i_nary1<fl_size_i_nary1;i_nary1++){
      counts.push_back(vector<uint>(0));
      for(uint i_alloy1=0,fl_size_i_alloy1=m_naries[i_nary1].m_alloys.size();i_alloy1<fl_size_i_alloy1;i_alloy1++){
        counts[i_nary1].push_back(0);
      }
    }
    for(uint i_nary1=0,fl_size_i_nary1=m_naries.size();i_nary1<fl_size_i_nary1;i_nary1++){
      for(uint i_alloy1=0,fl_size_i_alloy1=m_naries[i_nary1].m_alloys.size();i_alloy1<fl_size_i_alloy1;i_alloy1++){
        //m_naries[i_nary1].m_alloys[i_alloy1].m_elements_present is relevant hull dimensions
        const xvector<int>& elements_present_hull=m_naries[i_nary1].m_alloys[i_alloy1].m_elements_present;
        for(uint i_nary2=0;i_nary2<=i_nary1;i_nary2++){ //up to i_nary1
          for(uint i_alloy2=0,fl_size_i_alloy2=m_naries[i_nary2].m_alloys.size();i_alloy2<fl_size_i_alloy2;i_alloy2++){
            //space==m_naries[i_nary1].m_alloys[i_alloy1].m_elements_present, subspace==m_naries[i_nary2].m_alloys[i_alloy2].m_elements_present
            if(subspaceBelongs(elements_present_hull,m_naries[i_nary2].m_alloys[i_alloy2].m_elements_present)){
              counts[i_nary1][i_alloy1]+=getEntriesCount(i_nary2,i_alloy2,only_within_half_hull);
            }
          }
        }
      }
    }
    return counts;
  }

  uint ConvexHull::getGStateCount() const {
    uint i_point=AUROSTD_MAX_UINT,count=0;
    for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
      for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
        i_point=m_coord_groups[i_coord_group].m_points[i];
        if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
        if(m_points[i_point].m_is_g_state){count++;}
      }
    }
    return count;
  }

  uint ConvexHull::getGStateCount(uint i_nary) const {
    if(!m_has_stoich_coords){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Non-stoich coordinates");}
    uint i_point=AUROSTD_MAX_UINT,i_coord_group=AUROSTD_MAX_UINT,count=0;
    for(uint i_alloy=0,fl_size_i_alloy=m_naries[i_nary].m_alloys.size();i_alloy<fl_size_i_alloy;i_alloy++){
      for(uint i=0,fl_size_i=m_naries[i_nary].m_alloys[i_alloy].m_coord_groups.size();i<fl_size_i;i++){
        i_coord_group=m_naries[i_nary].m_alloys[i_alloy].m_coord_groups[i];
        if(!m_coord_groups[i_coord_group].m_is_on_hull){continue;}  //if it's not on the hull, definitely won't contain gstates
        for(uint j=0,fl_size_j=m_coord_groups[i_coord_group].m_points.size();j<fl_size_j;j++){
          i_point=m_coord_groups[i_coord_group].m_points[j];
          if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
          if(m_points[i_point].m_is_g_state){count++;}
        }
      }
    }
    return count;
  }

  vector<uint> ConvexHull::getHullPoints(bool sort_stoich_ascending) const { //pure hull-members, not equivalent ones too
    vector<uint> hull_points;
    uint i_point=AUROSTD_MAX_UINT;
    for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
      if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}
      for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
        i_point=m_coord_groups[i_coord_group].m_points[i];
        if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
        if(m_points[i_point].m_is_on_hull){hull_points.push_back(i_point);}
      }
    }
    std::sort(hull_points.begin(),hull_points.end(),sortCHullPoints(m_points,sort_stoich_ascending,true));
    return hull_points;
  }

  vector<uint> ConvexHull::getGStates(bool include_unaries,bool sort_stoich_ascending) const { //pure g_states, not equivalent ones too
    vector<uint> g_states;
    uint i_point=AUROSTD_MAX_UINT;
    for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
      if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}
      if(m_coord_groups[i_coord_group].m_is_on_hull){
        if(!isViablePoint(m_coord_groups[i_coord_group].m_hull_member)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No hull member set for m_coord_group["+aurostd::utype2string(i_coord_group)+"]");}
        i_point=artificialMap(m_coord_groups[i_coord_group].m_hull_member);
        if(!isViableGState(i_point)){continue;}  //could be legitimately missing g-state unary
        if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
        if(!m_points[i_point].m_is_g_state){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Mis-identified ground-state structure");}
        if(!(!include_unaries && m_points[i_point].isUnary())){g_states.push_back(i_point);}
      }
    }
    std::sort(g_states.begin(),g_states.end(),sortCHullPoints(m_points,sort_stoich_ascending,true));
    return g_states;
  }

  uint ConvexHull::getUnaryGState(uint i_alloy) const {
    uint i_nary=0;
    if(i_alloy>m_naries[i_nary].m_alloys.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within alloys");}
    uint g_state=AUROSTD_MAX_UINT;
    if(m_naries[i_nary].m_alloys[i_alloy].m_coord_groups.size()!=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unexpected count of coordgroups for unaries, should only be 1");}
    uint i_coord_group=m_naries[i_nary].m_alloys[i_alloy].m_coord_groups[0];
    uint i_point=AUROSTD_MAX_UINT;
    for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
      i_point=m_coord_groups[i_coord_group].m_points[i];
      if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
      if(m_points[i_point].m_is_g_state){
        g_state=i_point;
        return g_state;
      }
    }
    return g_state;
  }

  bool ConvexHull::isViablePoint(uint i_point) const {return i_point<m_points.size();}
  bool ConvexHull::isViableGState(uint g_state) const {return isViablePoint(g_state) && m_points[g_state].m_has_entry;}

  bool ConvexHull::findPoint(const string& auid,uint& i_point) const {
    uint _i_point;
    //do NOT go through m_points, this may include some duplicates we previously excluded (now not duplicates since we are removing points: AlFe hull)
    for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
      for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
        _i_point=m_coord_groups[i_coord_group].m_points[i];
        if(!m_points[_i_point].m_has_entry){continue;}
        if(m_points[_i_point].m_entry.auid==auid){i_point=_i_point;return true;}
      }
    }
    return false;
  }

  bool ConvexHull::findPoint(const xvector<double>& coords,uint& i_point) const {
    uint _i_point;
    //do NOT go through m_points, this may include some duplicates we previously excluded (now not duplicates since we are removing points: AlFe hull)
    for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
      for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
        _i_point=m_coord_groups[i_coord_group].m_points[i];
        if(aurostd::identical(m_points[_i_point].m_coords,coords,ZERO_TOL)){i_point=_i_point;return true;}
      }
    }
    return false;
  }

  bool ConvexHull::getNariesIndex(uint i_point,uint& i_nary,uint& i_alloy,uint& i_coord_group,bool redo) const {
    if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    return getNariesIndex(m_points[i_point],i_nary,i_alloy,i_coord_group,redo);
  }

  bool ConvexHull::getNariesIndex(const ChullPoint& point,uint& i_nary,uint& i_alloy,uint& i_coord_group,bool redo) const {
    if(!point.m_has_stoich_coords){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point does not have stoich coordinates");}
    return (getAlloyIndex(point,i_nary,i_alloy,redo) && getCoordGroupIndex(point,i_coord_group,redo));
  }

  bool ConvexHull::getCoordGroupIndex(uint i_point,uint& i_coord_group,bool redo) const {
    if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    return getCoordGroupIndex(m_points[i_point],i_coord_group,redo);
  }

  bool ConvexHull::getCoordGroupIndex(const ChullPoint& point,uint& i_coord_group,bool redo) const {
    if(!point.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized chullpoint");}
    if(!redo){
      i_coord_group=point.m_i_coord_group;
      if(i_coord_group<m_coord_groups.size()){return true;}
    }
    return getCoordGroupIndex(point.getStoichiometricCoords(),i_coord_group);
  }

  bool ConvexHull::getCoordGroupIndex(const xvector<double>& r_coords,uint& i_coord_group) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " r_coords=" << r_coords << endl;
      cerr << __AFLOW_FUNC__ << " m_coord_groups.size()=" << m_coord_groups.size() << endl;
      for(uint i=0,fl_size_i=m_coord_groups.size();i<fl_size_i;i++){
        cerr << __AFLOW_FUNC__ << " m_coord_groups[i=" << i << "].m_coords=" << m_coord_groups[i].m_coords << endl;
      }
    }
    for(uint i=0,fl_size_i=m_coord_groups.size();i<fl_size_i;i++){
      if(!m_coord_groups[i].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}
      if(coordsIdentical(m_coord_groups[i].m_coords,r_coords)){
        i_coord_group=i;
        return true;
      }
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_coord_groups[i=" << i << "].m_coords=" << m_coord_groups[i].m_coords << " != " << "r_coords=" << r_coords << endl;}
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " no coord_group index found for r_coords=" << r_coords << endl;}
    return false;
  }

  bool ConvexHull::getAlloyIndex(const ChullPoint& point,uint& i_nary,uint& i_alloy,bool redo) const {
    if(!point.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized chullpoint");}
    if(!point.m_has_stoich_coords){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point does not have stoich coordinates");}
    if(!redo){
      i_nary=point.m_i_nary;
      i_alloy=point.m_i_alloy;
      if(i_nary<m_naries.size() && i_alloy<m_naries[i_nary].m_alloys.size()){return true;}
    }
    return getAlloyIndex(point.m_elements_present,i_nary,i_alloy);
  }

  bool ConvexHull::getAlloyIndex(const CoordGroup& cg,uint& i_nary,uint& i_alloy,bool redo) const {
    if(!cg.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}
    if(!cg.m_has_stoich_coords){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup does not have stoich coordinates");}
    if(!cg.m_points.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup has no points");}
    if(!redo){
      i_nary=cg.m_i_nary;
      i_alloy=cg.m_i_alloy;
      if(i_nary<m_naries.size() && i_alloy<m_naries[i_nary].m_alloys.size()){return true;}
    }
    return getAlloyIndex(cg.getElementsPresent(),i_nary,i_alloy);
  }

  bool ConvexHull::getAlloyIndex(const xvector<int>& elements_present,uint& i_nary,uint& i_alloy) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}
    i_nary=sum(elements_present)-1;
    if(i_nary>m_naries.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " elements_present=" << elements_present << endl;}
    for(uint i=0,fl_size_i=m_naries[i_nary].m_alloys.size();i<fl_size_i;i++){
      if(!m_naries[i_nary].m_alloys[i].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_naries[i_nary=" << i_nary << "].m_alloys[i=" << i << "].m_elements_present=" << m_naries[i_nary].m_alloys[i].m_elements_present << endl;}
      if(m_naries[i_nary].m_alloys[i].m_elements_present==elements_present){
        i_alloy=i;
        return true;
      }
    }
    return false;
  }

  uint ConvexHull::artificialMap(uint i_point) const {
    if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_points");}
    if(!m_points[i_point].m_is_artificial){return i_point;}
    //we found an artificial point, but was it supposed to be here?
    if(!m_add_artificial_unaries){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid artificial point, not sure how to handle mapping");}

    const ChullPoint& art_point=m_points[i_point];
    uint i_coord_group=AUROSTD_MAX_UINT;
    if(!getCoordGroupIndex(art_point,i_coord_group)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup index not set");}
    if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}
    if(m_coord_groups[i_coord_group].m_points.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup["+aurostd::utype2string(i_coord_group)+"] has no points");}

    uint i_point_new=AUROSTD_MAX_UINT;
    for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
      i_point_new=m_coord_groups[i_coord_group].m_points[i];
      if(!m_points[i_point_new].m_is_artificial){return i_point_new;}
    }

    //if we get here, there are no viable replacements, simply return artificial point
    if(m_coord_groups[i_coord_group].m_points.size()==1){if(m_coord_groups[i_coord_group].m_points[0]==i_point){return i_point;}}

    //really bad if we get here
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot determine artificial point mapping");
  }

  bool ConvexHull::write(filetype ftype) const {
    bool written=false;
    try{
      if(ftype==chull_apool_ft){writeAPool();written=true;}
      else if(ftype==json_ft||ftype==txt_ft){writeText(ftype);written=true;}
      else if(ftype==latex_ft){writeLatex();written=true;}
      else if(ftype==chull_web_ft){writeWebApp();written=true;}
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return written;
  }

  void ConvexHull::setDefaultCFlags() {pflow::defaultLoadEntriesFlags(m_cflags, *p_FileMESSAGE, *p_oss, std::string("A"), false, true);}
  void ConvexHull::setCFlags(const aurostd::xoption& vpflow) {m_cflags=vpflow;}
  void ConvexHull::setDirectory() {m_aflags.Directory=getPath(m_cflags,*p_FileMESSAGE,*p_oss);}
  //MOVED TO xStream
  //void ConvexHull::setOFStream(ofstream& FileMESSAGE){p_FileMESSAGE=&FileMESSAGE;}
  //void ConvexHull::setOSS(ostream& oss) {p_oss=&oss;}

  bool ConvexHull::createHull(const string& alloy) {
    try{
      initializePoints(alloy);
      checkStructurePoints(); //some nice checks that everything checks out
      calculate();
      m_initialized=true;
      //hull must be initialized for these analyses
      thermodynamicsPostProcessing(); // will return if not m_thermo_hull
    }
    catch(aurostd::xerror& err){
      pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
      clear();
    }
    return m_initialized;
  }

  bool ConvexHull::createHull(const vector<string>& velements) {
    try{
      initializePoints(velements);
      checkStructurePoints(); //some nice checks that everything checks out
      calculate();
      m_initialized=true;
      //hull must be initialized for these analyses
      thermodynamicsPostProcessing(); // will return if not m_thermo_hull
    }
    catch(aurostd::xerror& err){
      pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
      clear();
    }
    return m_initialized;
  }

  bool ConvexHull::createHull(const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries) {
    try{
      initializePoints(velements,entries);
      checkStructurePoints(); //some nice checks that everything checks out
      calculate();
      m_initialized=true;
      //hull must be initialized for these analyses
      thermodynamicsPostProcessing(); // will return if not m_thermo_hull
    }
    catch(aurostd::xerror& err){
      pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
      clear();
    }
    return m_initialized;
  }

  bool ConvexHull::createHull(const vector<xvector<double> >& vcoords,bool has_stoich_coords,bool formation_energy_hull,bool add_artificial_unaries) {
    try{
      initializePoints(vcoords,has_stoich_coords,formation_energy_hull,add_artificial_unaries);
      checkStructurePoints(); //some nice checks that everything checks out
      calculate();
      m_initialized=true;
      //hull must be initialized for these analyses
      thermodynamicsPostProcessing(); // will return if not m_thermo_hull
    }
    catch(aurostd::xerror& err){
      pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
      clear();
    }
    return m_initialized;
  }

  bool ConvexHull::createHull(const vector<ChullPoint>& vpoints,bool formation_energy_hull,bool add_artificial_unaries) {
    try{
      initializePoints(vpoints,formation_energy_hull,add_artificial_unaries);
      checkStructurePoints(); //some nice checks that everything checks out
      calculate();
      m_initialized=true;
      //hull must be initialized for these analyses
      thermodynamicsPostProcessing(); // will return if not m_thermo_hull
    }
    catch(aurostd::xerror& err){
      pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
      clear();
    }
    return m_initialized;
  }

  bool ConvexHull::createHull(const vector<ChullPoint>& vpoints,const vector<string>& velements,bool formation_energy_hull,bool add_artificial_unaries) {
    try{
      initializePoints(vpoints,velements,formation_energy_hull,add_artificial_unaries);
      checkStructurePoints(); //some nice checks that everything checks out
      calculate();
      m_initialized=true;
      //hull must be initialized for these analyses
      thermodynamicsPostProcessing(); // will return if not m_thermo_hull
    }
    catch(aurostd::xerror& err){
      pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
      clear();
    }
    return m_initialized;
  }

  bool ConvexHull::entryValid(const aflowlib::_aflowlib_entry& entry,bool ignore_bad_database) const {
    string reason;
    return entryValid(entry,reason,ignore_bad_database);
  }

  bool ConvexHull::entryValid(const aflowlib::_aflowlib_entry& entry,string& reason,bool ignore_bad_database) const {
    char LOGGER_TYPE=_LOGGER_OPTION_;
    return entryValid(entry,reason,LOGGER_TYPE,ignore_bad_database);
  }

  bool ConvexHull::entryValid(const aflowlib::_aflowlib_entry& entry,string& reason,char& LOGGER_TYPE,bool ignore_bad_database) const {
    reason="";
    LOGGER_TYPE=_LOGGER_OPTION_;
    //tests of stupidity
    if(entry.vspecies.size()!=entry.vcomposition.size()){
      if(entry.prototype.find("POCC")==string::npos){ //POCC entries have no composition
        //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Bad entry ("+entry.auid+") - vspecies.size!=vcomposition.size()"); //let's not break the code for one bad entry
        reason="Entry (auid="+entry.auid+") is ill-defined: vspecies.size()!=vcomposition.size()";
        reason+=" (please report on AFLOW Forum: aflow.org/forum)";
        LOGGER_TYPE=_LOGGER_WARNING_;
        return false;
      }
    }
    bool found=false;
    for(uint j=0,fl_size_j=entry.vspecies.size();j<fl_size_j;j++){
      found=false;
      for(uint i=0,fl_size_i=m_velements.size();i<fl_size_i && !found;i++){if(m_velements[i]==entry.vspecies[j]){found=true;}}
      if(!found){
        reason="Entry (auid="+entry.auid+") contains "+entry.vspecies[j]+" and does not belong to hull";
        reason+=" (please report on AFLOW Forum: aflow.org/forum)";
        LOGGER_TYPE=_LOGGER_WARNING_;
        return false;
      }
    }
    if(!entry.ldau_TLUJ.empty()){
      reason="calculated with +U parameters";
      return false;
    }
    //mimic SC's approach from aflow_apennsy_gndstate.cpp
    //if(aurostd::substring2bool(entry.entry,"LDAU")){  //should be redundant, but less exclusive than looking for ldau_TLUJ
    //  reason="calculated with +U parameters";
    //  return false;
    //}
    //improved with aflowrc dft_type banning below
    //this helps with BSm, which was run mostly with PAW_GGA
    //if(entry.dft_type!="PAW_PBE" && !((ADD_SC_BSm_EXCEPTION)&&(entry.dft_type=="PAW_GGA" && aurostd::substring2bool(entry.entry,"B_hSm_3")))){ //SC exception
    //  reason="not calculated with PAW_PBE (dft_type=="+entry.dft_type+")";
    //  return false;
    //}
    //filters must be wide to narrow
    if(!m_allow_all_formation_energies){
      if(!aurostd::WithinList(m_allowed_dft_types,entry.dft_type)){
        reason="not calculated with allowed DFT-type: dft_type=="+entry.dft_type+" (allowed="+aurostd::joinWDelimiter(m_allowed_dft_types,",")+")";
        return false;
      }
    }
    if(aurostd::substring2bool(entry.entry,"NUPDOWN")){
      reason="calculated with manual NUPDOWN, thus E-fermi is NOT adjusted for spin-up"; //as explained in ovasp
      return false;
    }
    if(H_f_atom(entry)==AUROSTD_NAN || entry.entropic_temperature==AUROSTD_NAN){  //entry.enthalpy_formation_atom
      reason="enthalpy_formation_atom/entropic_temperature not calculated";
      return false;
    }
    if(ignore_bad_database && entry.ignoreBadDatabase(reason)){return false;}
    //otherwise return true
    return true;
  }

  void ConvexHull::addArtificialUnaries(uint dim){
    //points are really dim+1 dimensional (hidden dimension)
    //really, dim specifies within s_coords unless dim == s_coords.size()
    //then it's the last coord
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(!m_has_stoich_coords){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Non-stoich coordinates");}
    if(dim>m_dim-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid dimension requested");}
    xvector<double> dummy(m_dim);
    if(dim!=m_dim-1){dummy[dim+dummy.lrows]=1;}
    m_points.push_back(ChullPoint(dummy,*p_FileMESSAGE,*p_oss,m_has_stoich_coords,m_formation_energy_hull,true));
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " new artificial point=" << m_points.back().m_coords << endl;}
  }

  void ConvexHull::addArtificialUnaries(){for(uint i=0;i<m_dim;i++){addArtificialUnaries(i);}}

  bool ConvexHull::entryUnique(const vector<uint>& unique_entries,const aflowlib::_aflowlib_entry& entry,string& canonical_auid) const {
    //points have already been created, determined to be unique
    //hack, go backwards, as the way entries are ordered, duplicates occur near each other
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    canonical_auid="";
    for(uint fl_size_i=unique_entries.size(),i=fl_size_i-1;i<fl_size_i;i--){
      const ChullPoint& point=m_points[unique_entries[i]];
      if(point.entryIdentical(entry)){
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " entry[auid=" << point.m_entry.auid << ",compound=" << point.m_entry.compound << ",prototype=" << point.m_entry.prototype << "] == ";
          cerr << "entry[auid=" << entry.auid << ",compound=" << entry.compound << ",prototype=" << entry.prototype << "]" << endl;
        }
        canonical_auid=point.m_entry.auid;
        return false;
      }
    }
    return true;
  }

  void ConvexHull::loadPoints(const string& alloy) {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " initializing alloy, compound=" << alloy << endl;}
    vector<string> velements = aurostd::getElements(alloy, pp_string, *p_FileMESSAGE, true, true, false, *p_oss); //clean and sort, do not keep_pp  //CO20190712
    //[CO20190712 - OBSOLETE]vector<string> velements = pflow::getAlphabeticVectorString(alloy, *p_FileMESSAGE, *p_oss);
    return loadPoints(velements);
  }

  void ConvexHull::loadPoints(const vector<string>& _velements) {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);

    vector<string> velements; 
    for(uint i=0,fl_size_i=_velements.size();i<fl_size_i;i++){velements.push_back(_velements[i]);}
    std::sort(velements.begin(),velements.end()); //safe

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " initializing velements, velements=" << aurostd::joinWDelimiter(velements,",") << endl;}
    
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Loading raw entries from the database", m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    vector<vector<vector<aflowlib::_aflowlib_entry> > > entries;
    if(0){pflow::loadEntries(m_cflags,velements,entries,*p_FileMESSAGE,*p_oss);}
    else{
      //HE new method
      { // ensures that all resources are released as soon as possible
        aflowlib::EntryLoader el;
        el.m_out_debug = LDEBUG; // leave this on for testing
        el.loadAlloy(_velements); // the alloy will be loaded recursive by default, (elements will be cleaned and sorted)
        el.getEntriesThreeLayer(entries);
      }
      uint i_nary=0,i_alloy=0;
      bool sortLIB2=true;
      if(sortLIB2){
        uint i=0;
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " reshuffling unaries to optimize duplicate removal" << endl;}
        //[CO20221111 - VERY SLOW]for(i_alloy=0;i_alloy<entries[i_nary].size();i_alloy++){
        //[CO20221111 - VERY SLOW]  std::stable_sort(entries[i_nary][i_alloy].begin(),entries[i_nary][i_alloy].end(),sortLIB2Entries(velements,*p_FileMESSAGE,*p_oss));
        //[CO20221111 - VERY SLOW]}
        vector<_aflowlib_entry_LIB2sorting> vaes;
        uint i_entry=0;
        for(i_alloy=0;i_alloy<entries[i_nary].size();i_alloy++){
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " i_alloy=" << i_alloy << endl;}
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << " INITIAL ORDER=" << endl;
            for(i=0;i<entries[i_nary][i_alloy].size()-1;i++){
              cerr << __AFLOW_FUNC__ << " i_entry.auid=" << entries[i_nary][i_alloy][i].auid << " i_entry.aurl=" << entries[i_nary][i_alloy][i].aurl << endl;
            }
            cerr << endl;
          }
          vaes.clear();
          for(i_entry=0;i_entry<entries[i_nary][i_alloy].size();i_entry++){
            vaes.emplace_back(_aflowlib_entry_LIB2sorting(entries[i_nary][i_alloy][i_entry],i_entry,velements,*p_FileMESSAGE,*p_oss));
          }
          std::stable_sort(vaes.begin(),vaes.end());
          vector<uint> order_new;
          for(i_entry=0;i_entry<vaes.size();i_entry++){order_new.push_back(vaes[i_entry].m_index);}
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " order_new=" << aurostd::joinWDelimiter(order_new,",") << endl;}
          aurostd::reorder(entries[i_nary][i_alloy],order_new,1); //use mode==1
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << " FINAL ORDER=" << endl;
            for(i=0;i<entries[i_nary][i_alloy].size()-1;i++){
              cerr << __AFLOW_FUNC__ << " i_entry.auid=" << entries[i_nary][i_alloy][i].auid << " i_entry.aurl=" << entries[i_nary][i_alloy][i].aurl << endl;
            }
            cerr << endl;
          }
        }
        //[CO20221111 - VERY SLOW]//need to reshuffle unaries in LIB2
        //[CO20221111 - VERY SLOW]uint i_nary=0,i_alloy=0,i=0,j=0,k=0;
        //[CO20221111 - VERY SLOW]aflowlib::_aflowlib_entry i_entry,j_entry;
        //[CO20221111 - VERY SLOW]vector<string> i_vspecies,j_vspecies;
        //[CO20221111 - VERY SLOW]bool i_found_element=0,j_found_element;
        //[CO20221111 - VERY SLOW]for(i_alloy=0;i_alloy<entries[i_nary].size();i_alloy++){
        //[CO20221111 - VERY SLOW]  for(i=0;i<entries[i_nary][i_alloy].size()-1;i++){
        //[CO20221111 - VERY SLOW]    for(j=i+1;j<entries[i_nary][i_alloy].size();j++){
        //[CO20221111 - VERY SLOW]      const aflowlib::_aflowlib_entry i_entry=entries[i_nary][i_alloy][i];
        //[CO20221111 - VERY SLOW]      const aflowlib::_aflowlib_entry j_entry=entries[i_nary][i_alloy][j];
        //[CO20221111 - VERY SLOW]      if(!i_entry.catalog.empty()&&i_entry.catalog!="LIB2"){continue;}  //only a problem with unaries coming from LIB2
        //[CO20221111 - VERY SLOW]      if(!j_entry.catalog.empty()&&j_entry.catalog!="LIB2"){continue;}  //only a problem with unaries coming from LIB2
        //[CO20221111 - VERY SLOW]      if(i_entry.prototype==j_entry.prototype){
        //[CO20221111 - VERY SLOW]        i_vspecies=i_entry.getSpeciesAURL(*p_FileMESSAGE,*p_oss); //in this case, species from AUID is NOT the real species
        //[CO20221111 - VERY SLOW]        j_vspecies=j_entry.getSpeciesAURL(*p_FileMESSAGE,*p_oss); //in this case, species from AUID is NOT the real species
        //[CO20221111 - VERY SLOW]        if(LDEBUG){
        //[CO20221111 - VERY SLOW]          cerr << __AFLOW_FUNC__ << " i_vspecies=" << aurostd::joinWDelimiter(i_vspecies,",") << endl;
        //[CO20221111 - VERY SLOW]          cerr << __AFLOW_FUNC__ << " j_vspecies=" << aurostd::joinWDelimiter(j_vspecies,",") << endl;
        //[CO20221111 - VERY SLOW]        }
        //[CO20221111 - VERY SLOW]        i_found_element=true;for(k=0;k<i_vspecies.size()&&i_found_element;k++){if(!aurostd::WithinList(velements,i_vspecies[k])){i_found_element=false;}}
        //[CO20221111 - VERY SLOW]        j_found_element=true;for(k=0;k<j_vspecies.size()&&j_found_element;k++){if(!aurostd::WithinList(velements,j_vspecies[k])){j_found_element=false;}}
        //[CO20221111 - VERY SLOW]        if(LDEBUG){cerr << __AFLOW_FUNC__ << " i_found_element=" << i_found_element << " j_found_element=" << j_found_element << endl;}
        //[CO20221111 - VERY SLOW]        bool swap=false;
        //[CO20221111 - VERY SLOW]        if(i_found_element==true && j_found_element==true){
        //[CO20221111 - VERY SLOW]          for(k=0;k<i_vspecies.size();k++){if(j_vspecies[k]<i_vspecies[k]){swap=true;}}
        //[CO20221111 - VERY SLOW]        }
        //[CO20221111 - VERY SLOW]        else if(i_found_element==false && j_found_element==true){swap=true;}
        //[CO20221111 - VERY SLOW]        if(swap){
        //[CO20221111 - VERY SLOW]          cerr << __AFLOW_FUNC__ << " i_vspecies=" << aurostd::joinWDelimiter(i_vspecies,",") << endl;
        //[CO20221111 - VERY SLOW]          cerr << __AFLOW_FUNC__ << " j_vspecies=" << aurostd::joinWDelimiter(j_vspecies,",") << endl;
        //[CO20221111 - VERY SLOW]          cerr << __AFLOW_FUNC__ << " i_found_element=" << i_found_element << " j_found_element=" << j_found_element << endl;
        //[CO20221111 - VERY SLOW]          cerr << __AFLOW_FUNC__ << " i_entry.auid=" << entries[i_nary][i_alloy][i].auid << " i_entry.aurl=" << entries[i_nary][i_alloy][i].aurl << endl;
        //[CO20221111 - VERY SLOW]          cerr << __AFLOW_FUNC__ << " j_entry.auid=" << entries[i_nary][i_alloy][j].auid << " j_entry.aurl=" << entries[i_nary][i_alloy][j].aurl << endl;
        //[CO20221111 - VERY SLOW]          std::iter_swap(entries[i_nary][i_alloy].begin()+i,entries[i_nary][i_alloy].begin()+j);
        //[CO20221111 - VERY SLOW]          cerr << __AFLOW_FUNC__ << " i_entry.auid=" << entries[i_nary][i_alloy][i].auid << " i_entry.aurl=" << entries[i_nary][i_alloy][i].aurl << " SWAPPED" << endl;
        //[CO20221111 - VERY SLOW]          cerr << __AFLOW_FUNC__ << " j_entry.auid=" << entries[i_nary][i_alloy][j].auid << " j_entry.aurl=" << entries[i_nary][i_alloy][j].aurl << " SWAPPED" << endl;
        //[CO20221111 - VERY SLOW]        }
        //[CO20221111 - VERY SLOW]      }
        //[CO20221111 - VERY SLOW]    }
        //[CO20221111 - VERY SLOW]  }
        //[CO20221111 - VERY SLOW]}
      }
    }
    long long int count=0;
    uint i=0;
    for(i=0;i<entries.size();i++){for(uint j=0;j<entries[i].size();j++){count+=entries[i][j].size();}}
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Loaded "+aurostd::utype2string(count)+" raw entries from the database", m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_NOTICE_);
    return loadPoints(velements,entries);
  }

  void ConvexHull::loadPoints(const vector<string>& _velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries) {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    stringstream message;

    vector<string> velements;
    velements.clear(); for(uint i=0,fl_size_i=_velements.size();i<fl_size_i;i++){velements.push_back(_velements[i]);}
    std::sort(velements.begin(),velements.end()); //safe

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " initializing velements WITH entries, velements=" << aurostd::joinWDelimiter(velements,",") << endl;}

    m_has_stoich_coords=true;
    m_formation_energy_hull=!m_cflags.flag("CHULL::ENTROPIC_TEMPERATURE");    //energy vs. entropic_temperature hull
    m_half_hull=true;
    m_lower_hull=m_formation_energy_hull; //energy/entropic_temperature lower/upper hull
    m_add_artificial_unaries=true;

    vector<ChullPoint> points;
    ChullPoint cp;

    for(uint i=0,fl_size_i=entries.size();i<fl_size_i;i++){
      for(uint j=0,fl_size_j=entries[i].size();j<fl_size_j;j++){
        for(uint k=0,fl_size_k=entries[i][j].size();k<fl_size_k;k++){
          cp.initialize(velements,entries[i][j][k],*p_FileMESSAGE,*p_oss,m_formation_energy_hull);
          points.push_back(cp);
          //save icsd entries
          if(aurostd::substring2bool(entries[i][j][k].prototype,"_ICSD_")){m_icsd_entries.push_back(points.size()-1);}
        }
      }
    }
    if(!points.size()){
      message << "No entries loaded";
      //simply always die here, we cannot grab dimensionality of hull without ANY points
      if(0&&m_cflags.flag("FORCE")){pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);}
      else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);}
    }
    return loadPoints(points,velements,m_formation_energy_hull,m_add_artificial_unaries);
  }

  void ConvexHull::loadPoints(const vector<xvector<double> >& vcoords,bool has_stoich_coords,bool formation_energy_hull,bool add_artificial_unaries) {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    vector<ChullPoint> points;
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " initializing vcoords (has_stoich_coords==" << has_stoich_coords << "), count=" << vcoords.size() << endl;}
    ChullPoint cp;
    for(uint i=0,fl_size_i=vcoords.size();i<fl_size_i;i++){
      cp.initialize(vcoords[i],*p_FileMESSAGE,*p_oss,has_stoich_coords,formation_energy_hull);
      points.push_back(cp);
    }
    return loadPoints(points,formation_energy_hull,add_artificial_unaries);
  }

  void ConvexHull::loadPoints(const vector<ChullPoint>& vpoints,bool formation_energy_hull,bool add_artificial_unaries) {
    vector<string> velements;
    return loadPoints(vpoints,velements,formation_energy_hull,add_artificial_unaries);
  }

  void ConvexHull::loadPoints(const vector<ChullPoint>& vpoints,const vector<string>& velements,bool formation_energy_hull,bool add_artificial_unaries) {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    stringstream message;

    if(velements.size()){
      m_velements.clear(); for(uint i=0,fl_size_i=velements.size();i<fl_size_i;i++){m_velements.push_back(velements[i]);}
      std::sort(m_velements.begin(),m_velements.end()); //safe
    }

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " initializing chullpoints, count=" << vpoints.size() << endl;}
    for(uint i=0,fl_size_i=vpoints.size();i<fl_size_i;i++){m_points.push_back(vpoints[i]);};
    if(!m_points.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No points loaded, no way to determine dimensionality of hull");}

    //flag defaults
    m_formation_energy_hull=formation_energy_hull;    //energy vs. entropic_temperature hull
    if(m_formation_energy_hull){m_half_hull=m_lower_hull=true;} //default
    m_half_hull=(m_half_hull && m_cflags.flag("CHULL::FULL_HULL")==false); //energy/entropic_temperature lower/upper hull //override with flag from m_cflags  //HE20210510 - added CHULL::FULL_HULL
    m_lower_hull=(m_formation_energy_hull && m_cflags.flag("CHULL::FULL_HULL")==false); //energy/entropic_temperature lower/upper hull //override with flag from m_cflags  //HE20210510 - added CHULL::FULL_HULL

    m_add_artificial_unaries=add_artificial_unaries;

    //detect for coord types mixture!
    m_has_stoich_coords=( m_points[0].m_has_stoich_coords || m_has_stoich_coords );
    m_dim=m_points[0].m_coords.rows;
    if(m_dim<2){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"1D hulls are trivial");} //this MUST be true: chull cannot find facets for 1D hulls
    //test of stupidity
    for(uint i=0,fl_size_i=m_points.size();i<fl_size_i;i++){
      if(!m_points[i].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized point");}
      if(m_points[i].getDim()!=m_dim){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Dimension mismatch among points");}
      if(m_points[i].m_has_stoich_coords!=m_has_stoich_coords){
        message << "Mismatch among coord types (stoich vs. non-stoich coords), assuming non-stoich coords";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        m_has_stoich_coords=false;
        break;
      }
    }

    if(m_half_hull && m_has_stoich_coords && m_add_artificial_unaries){m_thermo_hull=true;}
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " thermo_hull=" << m_thermo_hull << endl;
      cerr << __AFLOW_FUNC__ << " has_stoich_coords=" << m_has_stoich_coords << endl;
      cerr << __AFLOW_FUNC__ << " half_hull=" << m_half_hull << endl;
      cerr << __AFLOW_FUNC__ << " lower_hull=" << m_lower_hull << endl;
      cerr << __AFLOW_FUNC__ << " add_artificial_unaries=" << m_add_artificial_unaries << endl;
    }

    //ensures proper construction of hull
    //duplicates DO NOT AFFECT performance/accuracy of algorithm
    //we ignore these points after hull construction anyway
    if(m_add_artificial_unaries){addArtificialUnaries();}

    //get s_coords
    if(m_has_stoich_coords){for(uint i=0,fl_size_i=m_points.size();i<fl_size_i;i++){m_points[i].setStoichCoords();}} //repetita iuvant

    //if(0){  //do NOT resort points, keep in same order as user input
    //  if(LDEBUG){cerr << __AFLOW_FUNC__ << " resorting all points (including artificial points) by coord/stoich (descending) and energy (ascending)" << endl;}
    //  std::sort(m_points.begin(),m_points.end());
    //}

    //DO NOT ADD/SUBTRACT/CHANGE-ORDER OF M_POINTS BEYOND THIS FUNCTION
  }

  void ConvexHull::calculateOutlierThreshold(const vector<double>& _energies,double& upper_threshold,double& lower_threshold) {
    xvector<double> energies=aurostd::vector2xvector<double>(_energies);
    return calculateOutlierThreshold(energies,upper_threshold,lower_threshold);
  }

  void ConvexHull::calculateOutlierThreshold(const xvector<double>& energies,double& upper_threshold,double& lower_threshold) {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    stringstream message;

    upper_threshold=AUROSTD_MAX_DOUBLE;   //effectively NOT a threshold
    lower_threshold=-AUROSTD_MAX_DOUBLE;  //effectively NOT a threshold

    //nice solution here! but only works for odd counts
    //http://en.cppreference.com/w/cpp/algorithm/nth_element
    bool iqr_method=true; //unfortunately, MAD is normal distribution dependent, NOT our case here

    uint iqr_count_threshold=4; //3 results in degenerate quartile indices
    if((uint)energies.rows<iqr_count_threshold){ //ALWAYS not enough points to do statistics (need 3 quartiles)
      message << "Not enough degrees of freedom for outlier detection analysis per interquartile-range (count=" << energies.rows << " < " << iqr_count_threshold << ")";
      if(m_cflags.flag("FAKE_HULL")){aurostd::StringstreamClean(message);}  //don't want to see these errors, they are expected
      else if(m_cflags.flag("CHULL::STRICT_OUTLIER_ANALYSIS")&&(!m_cflags.flag("FORCE"))){
        message << " (results may not be reliable). Terminating hull analysis.";throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_RANGE_);
      }else{
        message << ", skipping outlier analysis.";pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
      }
      return;
    }

    if(LDEBUG){xvector<double> temp_eng = energies; sort(temp_eng); cerr << "lastCoords(): " << temp_eng << endl;}
    double q1,q2,q3;
    //so we sort full anyway to be completely robust, should be easy considering how we sorted before
    aurostd::getQuartiles(energies,q1,q2,q3);  //we sort in here

    double lower_anchor,upper_anchor,range,multiplier;
    if(iqr_method){  //classical iqr measure of outlier
      lower_anchor=q1;
      upper_anchor=q3;
      range=q3-q1;      //interquartile range, iqr
      //multiplier=3.25;  //default=1.5, but we need to be more conservative from trials
    }else{  //absolute deviation around the median (MAD)
      lower_anchor=q2;
      upper_anchor=q2;
      range=aurostd::getMAD(energies,q2); //better iqr IF normal distribution, otherwise we need to know type of distribution (quartiles)
      //multiplier=3.25;                    //doi=10.1080/14640749108400962; 3 (very conservative), 2.5 (moderately conservative), 2 (poorly conservative)
    }
    multiplier=DEFAULT_CHULL_OUTLIER_MULTIPLIER;

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " lower_anchor=" << lower_anchor << ", median=" << q2 << ", upper_anchor=" << upper_anchor << ", range=" << range << endl;
    }

    upper_threshold=upper_anchor+(multiplier*range);
    lower_threshold=lower_anchor-(multiplier*range);
  }

  vector<uint> ConvexHull::calculateOutliers(const vector<uint>& points_to_consider) {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    //it is very important that we do not define outliers using std/mean, as these are
    //very sensitive to outliers
    //instead, use median
    //see discussion here:  doi=10.1016/j.jesp.2013.03.013
    vector<uint> outliers;

    //bool keep_outliers=m_cflags.flag("CHULL::INCLUDE_OUTLIERS");
    bool show_warnings=true;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " starting" << endl;}

    //get vector of last coords (we want to find outliers in this dimension)
    vector<double> _energies;  //vector and not xvector because we need push_back()
    uint i_point=AUROSTD_MAX_UINT;
    for(uint i=0,fl_size_i=points_to_consider.size();i<fl_size_i;i++){
      i_point=points_to_consider[i];
      if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
      _energies.push_back(m_points[i_point].getLastCoord());
    }

    double upper_threshold,lower_threshold;
    calculateOutlierThreshold(_energies,upper_threshold,lower_threshold);

    stringstream message;
    char LOGGER_TYPE=_LOGGER_OPTION_;
    if(show_warnings){LOGGER_TYPE=_LOGGER_WARNING_;}  //show warning if we do not remove!
    if(!(m_half_hull && !m_lower_hull)){  //look at lower range
      const double& threshold=lower_threshold;
      for(uint i=0,fl_size_i=points_to_consider.size();i<fl_size_i;i++){
        i_point=points_to_consider[i];
        if(m_points[i_point].getLastCoord()<threshold){
          outliers.push_back(i_point);
          message << "Identified (lower) outlier auid=" << m_points[i_point].m_entry.auid << ", ";
          message << "aurl=" << m_points[i_point].m_entry.aurl << ", ";
          message << "lastCoord()=" << m_points[i_point].getLastCoord() << " (<threshold=" << threshold << ")"; 
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, LOGGER_TYPE);
        }
      }
    }
    if(!(m_half_hull && m_lower_hull)){  //look at upper range
      const double& threshold=upper_threshold;
      for(uint i=0,fl_size_i=points_to_consider.size();i<fl_size_i;i++){
        i_point=points_to_consider[i];
        if(m_points[i_point].getLastCoord()>threshold){
          outliers.push_back(i_point);
          message << "Identified (upper) outlier auid=" << m_points[i_point].m_entry.auid << ", ";
          message << "aurl=" << m_points[i_point].m_entry.aurl << ", ";
          message << "lastCoord()=" << m_points[i_point].getLastCoord() << " (>threshold=" << threshold << ")"; 
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, LOGGER_TYPE);
        }
      }
    }
    //if we request outliers, lets get them, we can neglect them later
    //if(keep_outliers){
    //  message << "NOT removing outliers";
    //  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_OPTION_);
    //  outliers.clear();
    //}
    ////remove outliers
    //std::sort(outliers.rbegin(),outliers.rend()); //descending
    //for(uint i=0,fl_size_i=outliers.size();i<fl_size_i;i++){
    //  message << "Removing outlier auid=" << m_points[outliers[i]].m_entry.auid;
    //  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_OPTION_);
    //  m_points.erase(m_points.begin()+outliers[i]);
    //}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " done" << endl;}
    return outliers;
  }

  vector<uint> ConvexHull::getOutliers() {
    vector<uint> points_to_consider;
    uint i_point=AUROSTD_MAX_UINT;
    for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
      if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}
      for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
        i_point=m_coord_groups[i_coord_group].m_points[i];
        if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
        points_to_consider.push_back(i_point);
      }
    }
    return calculateOutliers(points_to_consider);
  }

  vector<uint> ConvexHull::getOutliers(const xvector<int>& elements_present) {
    stringstream message;
    vector<uint> points_to_consider;
    uint i_point=AUROSTD_MAX_UINT;
    char LOGGER_TYPE=_LOGGER_OPTION_;
    bool silent=0;
    bool see_neglect=m_cflags.flag("CHULL::SEE_NEGLECT");

    if(m_half_hull){  //we only care about points above/below hull
      for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
        if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}
        if(m_coord_groups[i_coord_group].getElementsPresent()!=elements_present){continue;}
        for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
          i_point=m_coord_groups[i_coord_group].m_points[i];
          if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
          if(m_points[i_point].isWithinHalfHull(m_lower_hull)){points_to_consider.push_back(i_point);} //!std::signbit(m_points[i_point].getLastCoord())) //positive
          else{
            silent=(!see_neglect && LOGGER_TYPE==_LOGGER_OPTION_);
            message << "Neglecting ";
            if(m_points[i_point].m_has_entry){message << "[auid=" << m_points[i_point].m_entry.auid << ",aurl=" << m_points[i_point].m_entry.aurl << "] ";}
            else{message << "[i_point=" << i_point << "] ";}
            message << "from outlier analysis: entry not within " << (m_lower_hull?"lower":"upper") << " half hull";
            pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, LOGGER_TYPE, silent);
          }
        }
      }
    }else{
      for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
        if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}
        if(m_coord_groups[i_coord_group].getElementsPresent()!=elements_present){continue;}
        for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
          i_point=m_coord_groups[i_coord_group].m_points[i];
          if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
          points_to_consider.push_back(i_point);
        }
      }
    }

    bool check_count=(m_half_hull&&sum(elements_present)==2);
    if(check_count){
      uint binaries_half_hull_threshold=DEFAULT_CHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES;
      if(points_to_consider.size()<binaries_half_hull_threshold){
        message << "Not enough degrees of freedom for outlier detection analysis per user defined threshold (count=" << points_to_consider.size() << " < " << binaries_half_hull_threshold << ")";
        if(m_cflags.flag("FAKE_HULL")){aurostd::StringstreamClean(message);}  //don't want to see these errors, they are expected
        else if(m_cflags.flag("CHULL::STRICT_OUTLIER_ANALYSIS")&&(!m_cflags.flag("FORCE"))){
          message << " (results may not be reliable). Terminating hull analysis.";throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_RANGE_);
        }else{
          message << ", skipping outlier analysis.";pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        }
        vector<uint> outliers;
        return outliers;
      }
    }

    return calculateOutliers(points_to_consider);
  }

  vector<uint> ConvexHull::findArtificialPoints(uint i_coord_group){
    if(i_coord_group>m_coord_groups.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within coordgroups");}
    if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}

    uint i_point=AUROSTD_MAX_UINT;
    vector<uint> artificial_points;
    for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
      i_point=m_coord_groups[i_coord_group].m_points[i];
      if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
      if(m_points[i_point].m_is_artificial){artificial_points.push_back(i_point);}
    }

    return artificial_points;
  }

  uint ConvexHull::findArtificialUnary(uint i_coord_group){
    if(!m_has_stoich_coords){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No unaries to be found (coordinates are not stoichiometric)");}
    if(i_coord_group>m_coord_groups.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within coordgroups");}
    if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}
    if(!m_coord_groups[i_coord_group].m_points.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup["+aurostd::utype2string(i_coord_group)+"] is empty");}
    if(!m_points[m_coord_groups[i_coord_group].m_points[0]].isUnary()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup["+aurostd::utype2string(i_coord_group)+"] is not unary");}

    vector<uint> artificial_points=findArtificialPoints(i_coord_group);
    if(artificial_points.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Missing artificial points");}
    if(artificial_points.size()!=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Too many artificial points: "+aurostd::utype2string(artificial_points.size()));}

    return artificial_points[0];
  }

  void ConvexHull::organizeHullPoints(uint i_coord_group) {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(i_coord_group>m_coord_groups.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within coordgroups");}
    if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " finding min/max energy point(s) for coord_group[" << i_coord_group << "]" << endl;}

    m_coord_groups[i_coord_group].m_candidate_hull_points.clear();
    uint p_size=m_coord_groups[i_coord_group].m_points.size();
    if(p_size==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup has no points");}

    //the sole purpose of an artificial point is to sit on the hull, ignore all others
    uint i_point=AUROSTD_MAX_UINT;
    if(m_coord_groups[i_coord_group].m_has_artificial_unary){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking for artificial point in coord_group[" << i_coord_group << "]" << endl;}
      i_point=findArtificialUnary(i_coord_group);
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " found artificial point in coord_group[" << i_coord_group << "]: ";
        cerr << "point[" << i_point << "]=" << m_points[i_point].m_coords << endl;
      }
      m_coord_groups[i_coord_group].m_candidate_hull_points.push_back(i_point);
      return;
    }

    if(m_half_hull){
      i_point=m_coord_groups[i_coord_group].m_points[0];        //lowest point for lower_hull, highest point for upper_hull
      if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
      if(m_lower_hull){
        if(aurostd::lessEqualZero(m_points[i_point].getLastCoord())) //std::signbit(m_points[i_point].getLastCoord()))
        { //CO20200106 - patching for auto-indenting
          m_coord_groups[i_coord_group].m_candidate_hull_points.push_back(i_point);
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << " lower half hull point found: ";
            cerr << "point[" << i_point << "]=" << m_points[i_point].m_coords << endl;
          }
        }
      }else{
        if(aurostd::greaterEqualZero(m_points[i_point].getLastCoord()))  //!std::signbit(m_points[i_point].getLastCoord()))
        { //CO20200106 - patching for auto-indenting
          m_coord_groups[i_coord_group].m_candidate_hull_points.push_back(i_point);
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << " upper half hull point found: ";
            cerr << "point[" << i_point << "]=" << m_points[i_point].m_coords << endl;
          }
        }
      }
      return; //below, we add other extreme in energy, which we don't care about for half hulls
    }

    i_point=m_coord_groups[i_coord_group].m_points[0];        //lowest energy point
    if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
    m_coord_groups[i_coord_group].m_candidate_hull_points.push_back(i_point);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " hull point found: ";
      cerr << "point[" << i_point << "]=" << m_points[i_point].m_coords << endl;
    }
    if(p_size>1){ //also grab other extreme
      i_point=m_coord_groups[i_coord_group].m_points[p_size-1];
      if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
      m_coord_groups[i_coord_group].m_candidate_hull_points.push_back(i_point);
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " hull point found: ";
        cerr << "point[" << i_point << "]=" << m_points[i_point].m_coords << endl;
      }
    }
  }

  void ConvexHull::organizeHullPoints() {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(!m_coord_groups.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Points have not been structured correctly");}
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " filtering points definitely NOT on the hull by energy" << endl;
      cerr << __AFLOW_FUNC__ << " only looking for min/max energy points in all coord_groups" << endl;
    }
    for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){organizeHullPoints(i_coord_group);}
  }

  void ConvexHull::initializeNaries() {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    stringstream message;
    //clear
    for(uint i_nary=0,fl_size_i_nary=m_naries.size();i_nary<fl_size_i_nary;i_nary++){
      for(uint i_alloy=0,fl_size_i_alloy=m_naries[i_nary].m_alloys.size();i_alloy<fl_size_i_alloy;i_alloy++){
        m_naries[i_nary].m_alloys[i_alloy].clear();
      }
      m_naries[i_nary].clear();
    }
    m_naries.clear();

    //initialize with combinations of dim
    Nary nary;
    Alloy alloy;
    aurostd::xcombos xc;
    xvector<int> elements_present;
    for(uint i_nary=0;i_nary<m_dim;i_nary++){
      nary.initialize(i_nary+1);
      m_naries.push_back(nary);
      xc.reset(m_dim,i_nary+1,'C');
      while(xc.increment()){
        elements_present=aurostd::vector2xvector<int>(xc.getCombo());
        alloy.initialize(elements_present);
        m_naries[i_nary].m_alloys.push_back(alloy);
      }
    }
    if(m_naries.size()==0){message << "m_naries.size()=0, xcombos may be broken";throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);}
    if(!m_naries[0].m_alloys.size()){message << "m_naries[0].m_alloys.size()==0, xcombos may be broken";throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);}
    //SORT NOW! do not sort later as we populate m_points with i_nary and i_alloy
    for(uint i_nary=0,fl_size_i_nary=m_naries.size();i_nary<fl_size_i_nary;i_nary++){
      std::sort(m_naries[i_nary].m_alloys.rbegin(),m_naries[i_nary].m_alloys.rend()); //descending order
    }
    std::sort(m_naries.begin(),m_naries.end());
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_naries.size()=" << m_naries.size() << endl;}
  }

  void ConvexHull::structurePoints() {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    stringstream message;
    m_coord_groups.clear();
    m_naries.clear();

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " starting" << endl;}

    //NOTE, m_points is UNTOUCHED (unsorted, fully populated) input, i.e., identical input of user
    //HOWEVER, m_coord_groups[].m_points only contain indices to points that we will consider for the hull calculation, i.e., does not contain
    //any points to be neglected (requested or otherwise)
    //if you want input of the user, use m_points
    //if you want points for the hull, go through m_coord_groups
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " gathering points to neglect" << endl;}

    bool fhsc_requested=m_cflags.flag("CHULL::CALCULATE_FAKE_HULL_STABILITY_CRITERION"); //CO20210315
    bool perform_structure_comparison=(1&&(!m_cflags.flag("CHULL::SKIP_STRUCTURE_COMPARISON"))); //(1&&!(m_cflags.flag("CHULL::SKIP_STRUCTURE_COMPARISON")||(!m_cflags.flag("CHULL::MULTI_OUTPUT")&&m_cflags.flag("CHULL::LATEX_DOC")&&m_cflags.flag("CHULL::IMAGE_ONLY"))));
    uint i_point_sc=AUROSTD_MAX_UINT;
    xvector<double> r_coords_sc_input;
    string auid_sc=m_cflags.getattachedscheme("CHULL::CALCULATE_FAKE_HULL_STABILITY_CRITERION");
    if(auid_sc.empty()){fhsc_requested=false;}
    if(fhsc_requested){
      for(uint i=0,fl_size_i=m_points.size();i<fl_size_i&&i_point_sc==AUROSTD_MAX_UINT;i++){
        const ChullPoint& point=m_points[i];
        const aflowlib::_aflowlib_entry& entry=m_points[i].m_entry;
        if(!point.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized point");}
        // start remove points
        if(m_points[i].m_has_entry && entry.auid==auid_sc){
          i_point_sc=i;
          r_coords_sc_input=point.getStoichiometricCoords();
        }
      }
    }
    if(i_point_sc==AUROSTD_MAX_UINT){fhsc_requested=false;}

    bool remove_requested=m_cflags.flag("CHULL::NEGLECT");
    bool see_neglect=m_cflags.flag("CHULL::SEE_NEGLECT");
    bool remove_submodular=true;  //remove AEL-AGL, APL, etc.
    bool fhn1eg_requested=m_cflags.flag("CHULL::CALCULATE_FAKE_HULL_N+1_ENTHALPY_GAIN"); //SK20200327
    bool remove_invalid=true;
    bool remove_duplicate_entries=true;        //we remove duplicate entries from the database, but in general, keep input of user constant
    bool remove_extreme=m_cflags.flag("CHULL::REMOVE_EXTREMA");
    bool perform_outliers_analysis=DEFAULT_CHULL_PERFORM_OUTLIER_ANALYSIS;
    bool keep_outliers=(!perform_outliers_analysis || m_cflags.flag("CHULL::INCLUDE_OUTLIERS"));
    bool remove_outliers=!keep_outliers;
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " remove_outliers=" << remove_outliers << endl;}

    vector<string> points_neglect;
    double extrema_val = 0.0;
    if(remove_requested){aurostd::string2tokens(m_cflags.getattachedscheme("CHULL::NEGLECT"),points_neglect,",");}
    if(points_neglect.size()==0){remove_requested=false;}
    if(remove_extreme){
      extrema_val=aurostd::string2utype<double>(m_cflags.getattachedscheme("CHULL::REMOVE_EXTREMA"));
      if(m_formation_energy_hull){
        if(aurostd::greaterEqualZero(extrema_val)){
          message << "Ignoring remove extreme points flag -- you provided a number >= 0. H_f convex hull sits below 0";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
          remove_extreme=false;
        }
      }else{
        if(aurostd::lessEqualZero(extrema_val)){
          message << "Ignoring remove extreme points flag -- you provided a number <= 0. T_S convex hull sits above 0";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
          remove_extreme=false;
        }
      }
    }

    vector<string> removing_messages;
    if(remove_requested){removing_messages.push_back("undesired");}
    if(remove_invalid){removing_messages.push_back("erroneous");}
    if(remove_duplicate_entries){removing_messages.push_back("duplicate");}
    if(remove_extreme){removing_messages.push_back("extreme");}
    if(remove_outliers){removing_messages.push_back("outlier");}
    if(fhn1eg_requested){removing_messages.push_back(aurostd::utype2string(m_velements.size())+"D (N+1 enthalpy gain)");}  //SK20200327
    if(fhsc_requested){removing_messages.push_back("auid="+auid_sc+" and equivalent (stability criterion)");}  //SK20200327
    if(removing_messages.size()){
      message << "Filtering out " << aurostd::joinWDelimiter(removing_messages,"/") << " entries";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }

    //quick add in of PAW_GGA if not already included in allowed list
    if(m_cflags.flag("CHULL::INCLUDE_PAW_GGA")&&(!aurostd::WithinList(m_allowed_dft_types,"PAW_GGA"))){m_allowed_dft_types.push_back("PAW_GGA");}

    //ignore bad database
    bool ignore_bad_database=(DEFAULT_CHULL_IGNORE_KNOWN_ILL_CONVERGED && !m_cflags.flag("CHULL::INCLUDE_ILL_CONVERGED"));
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " ignore_bad_database=" << ignore_bad_database << endl;}

    //organize into coordgroups
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " organizing into coordgroups" << endl;}
    vector<uint> unique_entries;
    string invalid_reason,canonical_auid;
    char LOGGER_TYPE=_LOGGER_OPTION_;
    bool silent=false;
    uint i_coord_group_sort;  //so it doesn't conflict with i_coord_group in for-loops
    CoordGroup cg;
    xvector<double> r_coords;
    for(uint i=0,fl_size_i=m_points.size();i<fl_size_i;i++){
      const ChullPoint& point=m_points[i];
      const aflowlib::_aflowlib_entry& entry=m_points[i].m_entry;
      if(!point.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized point");}
      // start remove points
      if(m_points[i].m_has_entry){
        if(remove_invalid && !entryValid(entry,invalid_reason,LOGGER_TYPE,ignore_bad_database)){
          if(!invalid_reason.empty()){
            silent=(!see_neglect && LOGGER_TYPE==_LOGGER_OPTION_);
            message << "Neglecting [auid=" << entry.auid << ",aurl=" << entry.aurl << "]: " << invalid_reason;
            pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, LOGGER_TYPE, silent);
          }
          continue;
        }
        if(remove_submodular){
          if(entry.aurl.find("ARUN.AEL_")!=string::npos ||
              entry.aurl.find("ARUN.AGL_")!=string::npos ||
              entry.aurl.find("ARUN.APL_")!=string::npos ||
              entry.aurl.find("ARUN.QHA_")!=string::npos ||
              FALSE){
            silent=true;  //no need to see
            message << "Neglecting [auid=" << entry.auid << ",aurl=" << entry.aurl << "]: sub-module load";
            pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, LOGGER_TYPE, silent);
            continue;
          }
        }
        if(remove_duplicate_entries && !entryUnique(unique_entries,entry,canonical_auid)){
          silent=(!see_neglect);
          message << "Neglecting [auid=" << entry.auid << ",aurl=" << entry.aurl << "]: duplicate database entry (see " << canonical_auid << ")";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_OPTION_, silent);
          continue;
        }
        unique_entries.push_back(i);
        if(aurostd::WithinList(points_neglect,entry.auid)){
          message << "Neglecting [auid=" << entry.auid << ",aurl=" << entry.aurl << "]: as requested";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_OPTION_);
          continue;
        }
        if(fhn1eg_requested) {  //SK20200330
          if(point.getDim()==(point.m_i_nary+1)) {
            message << "Neglecting [auid=" << entry.auid << ",aurl=" << entry.aurl << "] to calculate N+1 enthalpy gain";
            pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_OPTION_);
            continue;
          }
        }
        if(fhsc_requested){ //CO20210315
          if(coordsIdentical(r_coords_sc_input,point.getStoichiometricCoords()) && phasesEquivalent(i_point_sc,i,perform_structure_comparison)){
            message << "Neglecting [auid=" << entry.auid << ",aurl=" << entry.aurl << "] to calculate stability criterion hull";
            pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_OPTION_);
            continue;
          }
        }
        if(remove_extreme){
          if(m_formation_energy_hull){
            if(chull::H_f_atom(entry, _m_) < extrema_val){
              message << "Neglecting [auid=" << entry.auid << ",aurl=" << entry.aurl << "]: flagged as extreme with H_f = " << chull::H_f_atom(entry, _m_) << " (meV/atom)";
              pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_OPTION_);
              continue;
            }
          }else{
            if(chull::T_S(entry) > extrema_val){
              message << "Neglecting [auid=" << entry.auid << ",aurl=" << entry.aurl << "]: flagged as extreme with T_S = " << chull::T_S(entry) << " (K)";
              pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_OPTION_);
              continue;
            }
          }
        }
      }
      //// end remove points
      //organize into m_naries
      //perhaps think about organizing by c_coords instead in the future (int compare vs. double compare), would be more exact
      //this would require adding a priority check for composition, and falling back to stoichiometry if composition is not available (artificial)
      //not really a priority, we don't explore compositions that are so close (to within 1e-8)
      //comparison of stoichiometry is more general, but less reliable than composition (to within tol)
      //nevermind, vcomposition is a double anyway (in anticipation for POCC), so floating point comparisons are inevitable
      r_coords=point.getStoichiometricCoords();
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " point[" << i << "]: m_coords=" << point.m_coords << ", r_coords=" << r_coords << ", compound=\"" << point.m_entry.compound << "\", dim=" << point.m_i_nary+1 << endl;}
      if(!getCoordGroupIndex(r_coords,i_coord_group_sort)){
        cg.initialize(r_coords,point.m_has_stoich_coords);
        m_coord_groups.push_back(cg);
        i_coord_group_sort=m_coord_groups.size()-1;
      }
      m_coord_groups[i_coord_group_sort].m_points.push_back(i);
      if(m_coord_groups[i_coord_group_sort].m_has_stoich_coords && !point.m_has_stoich_coords){
        message << "Mismatch among coord types (stoich vs. non-stoich coords), assuming non-stoich coords";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        m_coord_groups[i_coord_group_sort].m_has_stoich_coords=false;
      }
      if(point.isUnary()&&point.m_is_artificial){m_coord_groups[i_coord_group_sort].m_has_artificial_unary=true;}
    }

    //remove outliers before sort (BY MATCHING COORDS)
    //we do this analysis EVEN if we keep outliers in the end, simply print out useful warnings for user
    uint i_point=AUROSTD_MAX_UINT;
    vector<uint> outliers;

    if(perform_outliers_analysis){
      if(m_has_stoich_coords){
        //proceed through each set of binaries, find all sets of outliers and append together
        xvector<int> elements_present;
        vector<uint> _outliers;
        aurostd::xcombos xc(m_dim,2,'C'); //binaries ONLY for now
        while(xc.increment()){
          elements_present=aurostd::vector2xvector<int>(xc.getCombo());
          _outliers=getOutliers(elements_present);
          outliers.insert(outliers.end(),_outliers.begin(),_outliers.end());
        }
      }else{outliers=getOutliers();}
    }

    if(keep_outliers){
      message << "NOT removing outliers";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_OPTION_);
      outliers.clear();
    }

    bool found_outlier=false;
    vector<uint> points_to_remove;
    uint valid_count=0;
    for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
      points_to_remove.clear();
      for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
        i_point=m_coord_groups[i_coord_group].m_points[i];
        if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
        found_outlier=false;
        for(uint j=0,fl_size_j=outliers.size();j<fl_size_j&&!found_outlier;j++){if(i_point==outliers[j]){found_outlier=true;}}
        if(found_outlier){points_to_remove.push_back(i);} //not i_point, so I can remove this index
        else{if(!m_points[i_point].m_is_artificial){valid_count++;}}
      }
      if(points_to_remove.size()){
        std::sort(points_to_remove.rbegin(),points_to_remove.rend()); //descending
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " before outlier removal, count = " << m_coord_groups[i_coord_group].m_points.size() << endl;}
        for(uint i=0,fl_size_i=points_to_remove.size();i<fl_size_i;i++){
          i_point=m_coord_groups[i_coord_group].m_points[points_to_remove[i]];
          if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
          message << "Removing outlier ";
          if(m_points[i_point].m_has_entry){message << "auid=" << m_points[i_point].m_entry.auid;}
          else{message << "m_coords=" << m_points[i_point].m_coords;}
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_OPTION_);
          m_coord_groups[i_coord_group].m_points.erase(m_coord_groups[i_coord_group].m_points.begin()+points_to_remove[i]);
        }
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " after outlier removal, m_coord_groups[" << i_coord_group << "].m_points.size() = " << m_coord_groups[i_coord_group].m_points.size() << endl;}
      }
    }
    message << "Employing " << valid_count << " total entries for " << pflow::arity_string(m_dim,false,false) << " convex hull analysis";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_NOTICE_);

    //remove empty m_coord_groups
    vector<uint> empty_coord_groups;
    for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
      if(m_coord_groups[i_coord_group].m_points.size()==0){empty_coord_groups.push_back(i_coord_group);}
    }
    std::sort(empty_coord_groups.rbegin(),empty_coord_groups.rend()); //descending
    for(uint i=0,fl_size_i=empty_coord_groups.size();i<fl_size_i;i++){m_coord_groups.erase(m_coord_groups.begin()+empty_coord_groups[i]);}

    //sort
    m_sort_energy_ascending=!(m_half_hull==true && m_lower_hull==false); //upper half hull should sort DESCENDING (ground-state first)
    for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
      std::sort(m_coord_groups[i_coord_group].m_points.begin(),m_coord_groups[i_coord_group].m_points.end(),sortWithinCoordGroup(m_points,m_sort_energy_ascending));  //ascending order
    }
    std::sort(m_coord_groups.rbegin(),m_coord_groups.rend()); //descending order for alphabetic print out later

    //assign coord group indices to points, useful later
    for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
      //set ref state
      if(m_coord_groups[i_coord_group].m_points.size()){
        m_coord_groups[i_coord_group].m_ref_state=artificialMap(m_coord_groups[i_coord_group].m_points[0]);
      }
      for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
        i_point=m_coord_groups[i_coord_group].m_points[i];
        if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
        m_points[i_point].m_i_coord_group=i_coord_group;
      }
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " done organizing into coordgroups" << endl;}

    organizeHullPoints();

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " coord_groups structure:" << endl;
      for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
        cerr << __AFLOW_FUNC__ << " coord_group[" << i_coord_group << "] coords=" << m_coord_groups[i_coord_group].m_coords << endl;
        for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
          cerr << __AFLOW_FUNC__ << " point[" << i << "] compound=" << m_points[m_coord_groups[i_coord_group].m_points[i]].m_entry.compound << " ";
          cerr << "coords=" << m_points[m_coord_groups[i_coord_group].m_points[i]].m_coords << " ";
          cerr << "has_stoich_coords=" << m_points[m_coord_groups[i_coord_group].m_points[i]].m_has_stoich_coords << " ";
          cerr << "formation_energy_coord=" << m_points[m_coord_groups[i_coord_group].m_points[i]].m_formation_energy_coord << endl;
        }
      }
    }

    if(m_has_stoich_coords){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " stoich_coords found, also sorting into n-aries and alloys" << endl;}
      message << "Stoichiometric coordinates detected, structuring entries by arity";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

      initializeNaries(); //create empty bins first, we can do this combinatorially

      if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_coord_groups.size()=" << m_coord_groups.size() << endl;}

      uint i_nary,i_alloy;
      for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_coord_groups[i_coord_group=" << i_coord_group << "].m_points.size()=" << m_coord_groups[i_coord_group].m_points.size() << endl;}
        if(!m_coord_groups[i_coord_group].m_points.size()){continue;}
        //we already filled bins, if we cannot find alloy system, then it's a bust
        if(!getAlloyIndex(m_coord_groups[i_coord_group],i_nary,i_alloy,true)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot get alloys index");}  //redo as we never done before
        m_coord_groups[i_coord_group].m_i_nary=i_nary;
        m_coord_groups[i_coord_group].m_i_alloy=i_alloy;
        m_naries[i_nary].m_alloys[i_alloy].m_coord_groups.push_back(i_coord_group);
        //set i_nary and i_alloy to points too
        for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
          i_point=m_coord_groups[i_coord_group].m_points[i];
          if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
          m_points[i_point].m_i_nary=i_nary;
          m_points[i_point].m_i_alloy=i_alloy;
        }
      }

      bool only_within_half_hull=false; //m_half_hull; //these counts should reflect total entries
      vector<vector<uint> > hull_sizes=getHullSizes(only_within_half_hull);
      if(only_within_half_hull){
        message << "Half hull detected, reducing entry count to those within the relevant hemisphere";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      }
      for(uint i_nary=0,fl_size_i_nary=m_naries.size();i_nary<fl_size_i_nary;i_nary++){
        for(uint i_alloy=0,fl_size_i_alloy=m_naries[i_nary].m_alloys.size();i_alloy<fl_size_i_alloy;i_alloy++){
          message << "Entries structure: employing " << getEntriesCount(i_nary,i_alloy,only_within_half_hull);
          if(m_velements.size()){message << " " << aurostd::joinWDelimiter(alloyToElements(i_nary,i_alloy),"-");}
          message << " [i_nary=" << i_nary << ",i_alloy=" << i_alloy << "]";
          message << " entries, ";
          message << hull_sizes[i_nary][i_alloy] << " entries total for ";
          message << pflow::arity_string(i_nary+1,false,false) << " convex hull analysis";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
        }
      }

      if(LDEBUG){
        for(uint i_nary=0,fl_size_i_nary=m_naries.size();i_nary<fl_size_i_nary;i_nary++){
          cerr << __AFLOW_FUNC__ << " " << pflow::arity_string(i_nary+1,false,false) << ":" << endl;
          for(uint i_alloy=0,fl_size_i_alloy=m_naries[i_nary].m_alloys.size();i_alloy<fl_size_i_alloy;i_alloy++){
            cerr << __AFLOW_FUNC__ << " alloy[" << i_alloy << "]: elements_present=" << m_naries[i_nary].m_alloys[i_alloy].m_elements_present << endl;
            for(uint i=0,fl_size_i=m_naries[i_nary].m_alloys[i_alloy].m_coord_groups.size();i<fl_size_i;i++){
              i_coord_group_sort=m_naries[i_nary].m_alloys[i_alloy].m_coord_groups[i];
              cerr << __AFLOW_FUNC__ << " coord_group[" << i_coord_group_sort << "]=" << m_coord_groups[i_coord_group_sort].m_coords << endl;
              for(uint j=0,fl_size_j=m_coord_groups[i_coord_group_sort].m_points.size();j<fl_size_j;j++){
                cerr << __AFLOW_FUNC__ << " compound=";
                cerr << (m_points[m_coord_groups[i_coord_group_sort].m_points[j]].m_is_artificial ? 
                    string("ARTIFICIAL") : 
                    m_points[m_coord_groups[i_coord_group_sort].m_points[j]].m_entry.compound) << " ";
                cerr << "coords=" << m_points[m_coord_groups[i_coord_group_sort].m_points[j]].m_coords << endl;
              }
            }
          }
        }
      }
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " done" << endl;}
  }

  vector<string> ConvexHull::alloyToElements(const ChullPoint& point) const {return alloyToElements(point.m_i_nary,point.m_i_alloy);}
  vector<string> ConvexHull::alloyToElements(uint i_nary,uint i_alloy) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}
    const xvector<int>& elements_present=getElementsPresent(i_nary,i_alloy);
    if((uint)elements_present.rows!=m_velements.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Dimension mismatch between elements present and composition");}
    vector<string> vout;
    for(int i=elements_present.lrows;i<=elements_present.urows;i++){
      if(elements_present[i]==1){vout.push_back(m_velements[i-elements_present.lrows]);}
    }
    return vout;
  }

  void ConvexHull::checkStructurePoints() {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    stringstream message;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " starting" << endl;}

    bool binary_statistics_check=m_naries.size()>1; //unary hull //true;

    if(m_has_stoich_coords){
      //UNARIES - START
      //tests of stupidity
      if(m_naries.size()!=m_dim){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Missing n-aries");} //these are populated by default if done correctly
      uint i_nary=0;
      if(m_naries[i_nary].m_alloys.size()!=m_dim){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Missing unary alloys");} //these are populated by default if done correctly

      //check for ground-states (non-artificial points)
      bool found_real = false, found_artificial = false;
      uint i_coord_group = 0,i_point = 0,i_point_real = 0;
      for(uint i_alloy=0,fl_size_i_alloy=m_naries[i_nary].m_alloys.size();i_alloy<fl_size_i_alloy;i_alloy++){
        if(m_naries[i_nary].m_alloys[i_alloy].m_coord_groups.size()!=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unexpected count of coordgroups for unaries, should only be 1");}
        found_real=found_artificial=false;
        i_coord_group=m_naries[i_nary].m_alloys[i_alloy].m_coord_groups[0];
        for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i&&!(found_real&&found_artificial);i++){
          i_point=m_coord_groups[i_coord_group].m_points[i];
          if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
          if(m_points[i_point].m_is_artificial){found_artificial=true;}
          else{
            if(!found_real){i_point_real=i_point;}  //grab first i_point that's real
            found_real=true;
          }
        }
        if(m_add_artificial_unaries&&!found_artificial){
          message << "Missing artificial points for";
          if(i_alloy<m_velements.size()){message << " " << m_velements[i_alloy];}
          message << " [i_nary=" << i_nary << ",i_alloy=" << i_alloy << "]";
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message.str()); message.str("");
        }
        if(found_real){
          if(m_thermo_hull){
            if(abs(m_points[i_point_real].getLastCoord())>=ENERGY_TOL){
              message << "Very skewed ground-state end point for";
              if(i_alloy<m_velements.size()){message << " " << m_velements[i_alloy];}
              message << " [i_nary=" << i_nary << ",i_alloy=" << i_alloy << "]";
              message << " (auid=" << m_points[i_point_real].m_entry.auid << ")";
              message << ": abs(" << m_points[i_point_real].getLastCoord() << ")>=" << ENERGY_TOL << " [eV]";
              message << " (please report on AFLOW Forum: aflow.org/forum)";
              if(m_cflags.flag("FAKE_HULL")){aurostd::StringstreamClean(message);}  //don't want to see these errors, they are expected
              else if(m_cflags.flag("FORCE")||m_cflags.flag("CHULL::INCLUDE_SKEWED_HULLS")){pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);}
              else{
                message << ". Override with --force (results may not be reliable).";
                throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
              }
            }
          }
        }else{
          if(!m_cflags.flag("FAKE_HULL")){
            message << "No ground-state available for";
            if(i_alloy<m_velements.size()){message << " " << m_velements[i_alloy];}
            message << " [i_nary=" << i_nary << ",i_alloy=" << i_alloy << "]";
            pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
          }
        }
      }
      //UNARIES - STOP

      //BINARIES - START
      uint count_threshold_binaries_total=DEFAULT_CHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES;
      if(binary_statistics_check){
        i_nary=1;
        uint count;
        for(uint i_alloy=0,fl_size_i_alloy=m_naries[i_nary].m_alloys.size();i_alloy<fl_size_i_alloy;i_alloy++){
          count=getEntriesCount(i_nary,i_alloy,false);  //check total first
          if(count<count_threshold_binaries_total){
            stringstream hull;
            //safe in case m_velements was not set
            if(m_velements.size()){hull << aurostd::joinWDelimiter(alloyToElements(i_nary,i_alloy),"-") << " ";}
            hull << "[i_nary=" << i_nary << ",i_alloy=" << i_alloy << "]";
            message << pflow::arity_string(i_nary+1,true,false) <<  " hull " << hull.str() << " is unreliable (total_entry_count=" << count << " < " << count_threshold_binaries_total << ")";
            if(m_cflags.flag("FAKE_HULL")){aurostd::StringstreamClean(message);}  //don't want to see these errors, they are expected
            else if(m_cflags.flag("FORCE")||m_cflags.flag("CHULL::INCLUDE_UNRELIABLE_HULLS")){pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);}
            else{
              message << ". Override with --force (results may not be reliable).";
              throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
            }
          }
        }
      }
      //BINARIES - STOP
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " done" << endl;}
  }

  void ConvexHull::addPointToFacet(ChullFacet& facet,uint i_point) {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    facet.addVertex(m_points[i_point],i_point);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " facet.m_vertices.size()=" << facet.m_vertices.size() << endl;
      for(uint i=0,fl_size_i=facet.m_vertices.size();i<fl_size_i;i++){
        cerr << facet.m_vertices[i].ch_index << "  h_coords=" << facet.m_vertices[i].ch_point.h_coords << endl;
      }
    }
  }

  void ConvexHull::initializeFacet(ChullFacet& facet,bool check_validity) {facet.initialize(h_reference,h_dim,check_validity);}

  uint ConvexHull::getExtremePoint(uint dim) {
    vector<FacetPoint> points_to_avoid;
    return getExtremePoint(dim,points_to_avoid);
  }

  uint ConvexHull::getExtremePoint(uint dim,const vector<FacetPoint>& points_to_avoid) {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    stringstream message;
    if(!h_points.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No starting points provided");}
    uint i_point=AUROSTD_MAX_UINT;
    double extreme=-1.0;  //since it's the abs() we check below, -1 is FINE
    bool avoid=false;
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " h_points:" << endl;
      for(uint i=0,fl_size_i=h_points.size();i<fl_size_i;i++){
        cerr << h_points[i] << "  h_coords=" << m_points[h_points[i]].h_coords << endl;
      }
      cerr << __AFLOW_FUNC__ << " points_to_avoid:" << endl;
      for(uint j=0,fl_size_j=points_to_avoid.size();j<fl_size_j;j++){
        cerr << points_to_avoid[j].ch_index << "  h_coords=" << m_points[points_to_avoid[j].ch_index].h_coords << endl;
      }
    }
    //HE20220412 START
    //collect directions already spanned by `points_to_avoid`
    vector<xvector<double>> directions;
    for (size_t avoid_i=1; avoid_i<points_to_avoid.size(); avoid_i++){
      // using the first extreme point as our starting point for all directions
      directions.push_back(m_points[points_to_avoid[avoid_i].ch_index].h_coords - m_points[points_to_avoid[0].ch_index].h_coords);
    }
    //HE20220412 END
    int h_coords_index=0;
    xvector<double> new_direction;
    for(uint i=0,fl_size_i=h_points.size();i<fl_size_i;i++){
      h_coords_index=(int)dim+m_points[h_points[i]].h_coords.lrows;
      if(h_coords_index>m_points[h_points[i]].h_coords.urows){
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid coordinate index");
      }
      avoid=false;
      for(uint j=0,fl_size_j=points_to_avoid.size();j<fl_size_j&&!avoid;j++){
        if(h_points[i]==points_to_avoid[j].ch_index){
          avoid=true;
          break;
        }
      } //HE20220412 added break to end loop when point is found
      if(avoid){continue;}
      //ensure that considered points add a new dimension to the initial facet (not collinear)
      //without this check, the first facet may not be spanning a dim-1 space
      //resulting in a wrong normal vector that trips up the search algorithm for the next point to add
      //which can lead to an endless loop in ConvexHull::calculateFacets()
      //this check became necessary as the creation of atomic environments uses ConvexHull to create hulls from arbitrary points
      //HE20220412 START
      if (!directions.empty()) {
        new_direction = m_points[h_points[i]].h_coords - m_points[points_to_avoid[0].ch_index].h_coords;
        for (vector<xvector<double>>::const_iterator direction = directions.begin(); direction != directions.end(); direction++) {
          if (aurostd::isCollinear(*direction, new_direction, ZERO_TOL))  {
            avoid=true;
            break;
          }
        }
        if(avoid){continue;}
      }
      //HE20220412 END
      if(abs(m_points[h_points[i]].h_coords[h_coords_index])>extreme){
        i_point=h_points[i];
        extreme=abs(m_points[h_points[i]].h_coords[h_coords_index]);
      }
    }
    if(i_point==AUROSTD_MAX_UINT){
      message << "No point found, points_to_avoid is too restrictive (points_to_avoid.size()==" << points_to_avoid.size() << ")";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message.str()); message.str("");
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " i_point=" << i_point << endl;}
    return i_point;
  }

  void ConvexHull::setCentroid() {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    vector<xvector<double> > points;
    for(uint i=0,fl_size_i=h_points.size();i<fl_size_i;i++){points.push_back(m_points[h_points[i]].h_coords);}
    h_centroid=h_reference=aurostd::getCentroid(points);  //fix h_reference later
    //if(1){
    //  //h_centroid[0]=0.5;h_centroid[1]=0.5;h_centroid[2]=-0.01;
    //  for(int i=h_centroid.lrows;i<=h_centroid.urows-1;i++){h_centroid[i]=0.5;}
    //  h_centroid[h_centroid.urows]=-0.01;
    //}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " centroid: " << h_centroid << endl;}
  }

  vector<FacetPoint> ConvexHull::getInitialExtremePoints() {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " start" << endl;}
    ChullFacet facet(*p_FileMESSAGE,*p_oss);  //CO20180305
    string error;
    //get first h_dim points from just the extremes, and plug into facet
    //BAD CASE: fewer than h_dim points, this is NOT hull-able
    if(h_points.size()<h_dim){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Not enough points to build an initial facet");}
    for(uint i=0;i<h_dim;i++){addPointToFacet(facet,getExtremePoint(i,facet.m_vertices));}
    if(!facet.isValid(error)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"An initial facet cannot be created per initial extreme points, invalid points input: "+error);}
    initializeFacet(facet,false); //already validated
    //now get furthest point from this facet, initialize outside set first
    FacetPoint f_point;
    for(uint i=0,fl_size_i=h_points.size();i<fl_size_i;i++){
      if(facet.isPointOnFacet(h_points[i])){continue;}
      f_point.initialize(m_points[h_points[i]],h_points[i]);
      facet.f_outside_set.push_back(f_point);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " adding to facet[" << 0 << "] outside set: " << h_points[i] << endl;}
    }
    if(facet.f_outside_set.size()){ //in special case where there are only 2 points in 2D, we simply return facet and not simplex
      facet.setFurthestPoint();
      facet.addVertex(facet.f_furthest_point); //already have a vector<uint> in m_points, simply append and return
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " building initial simplex with extreme points:" << endl;
      for(uint i=0,fl_size_i=facet.m_vertices.size();i<fl_size_i;i++){
        cerr << facet.m_vertices[i].ch_point.h_coords << endl;
      }
    }
    return facet.m_vertices;
  }

  void ConvexHull::setNeighbors() {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);

    //print ridges so far
    if(LDEBUG){
      vector<uint> indices;
      for(uint i=0,fl_size_i=h_facets.size();i<fl_size_i;i++){
        indices=h_facets[i].getCHIndices();
        cerr << __AFLOW_FUNC__ << " facet[" << i << "] has indices: ";
        for(uint j=0,fl_size_j=indices.size();j<fl_size_j;j++){
          cerr << indices[j];
          if(j!=indices.size()-1){cerr << ",";}
        }
        cerr << endl;
      }
    }

    for(uint i_facet=0,fl_size_i_facet=h_facets.size();i_facet<fl_size_i_facet;i_facet++){h_facets[i_facet].f_neighbors.clear();}
    for(uint i=0,fl_size_i=h_facets.size();i<fl_size_i;i++){
      for(uint j=i+1,fl_size_j=h_facets.size();j<fl_size_j;j++){
        if(h_facets[i].shareRidge(h_facets[j])){
          h_facets[i].f_neighbors.push_back(j);
          h_facets[j].f_neighbors.push_back(i);
        }
      }
      //can only perform this check if we have more than one h_facet, otherwise NO neighbors
      if(h_facets.size()>1&&h_facets[i].f_neighbors.size()!=h_dim){
        stringstream message;
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " neighbors for facet[i_facet=" << i << ",coords=";
          for(uint j=0,fl_size_j=h_facets[i].m_vertices.size();j<fl_size_j;j++){cerr << h_facets[i].m_vertices[j].ch_point.h_coords << (j!=h_facets[i].m_vertices.size()-1?", ":"");}
          cerr << "]: ";
          for(uint j=0,fl_size_j=h_facets[i].f_neighbors.size();j<fl_size_j;j++){
            cerr << h_facets[i].f_neighbors[j] << " ";
          }
          cerr << endl;
          cerr << __AFLOW_FUNC__ << " is_vertical=" << h_facets[i].m_is_vertical << endl;
        }
        message << "Neighbor count (" << h_facets[i].f_neighbors.size() << ") and facet dimension (" << h_dim << ") mismatch";
        if(0&&m_cflags.flag("FORCE")){pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);}  //h_facets[i].m_is_vertical
        else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);}
      }
    }
    for(uint i=0,fl_size_i=h_facets.size();i<fl_size_i;i++){std::sort(h_facets[i].f_neighbors.begin(),h_facets[i].f_neighbors.end());}

    if(LDEBUG){
      for(uint i=0,fl_size_i=h_facets.size();i<fl_size_i;i++){
        cerr << __AFLOW_FUNC__ << " neighbors for facet[i_facet=" << i << ",coords=";
        for(uint j=0,fl_size_j=h_facets[i].m_vertices.size();j<fl_size_j;j++){cerr << h_facets[i].m_vertices[j].ch_point.h_coords << (j!=h_facets[i].m_vertices.size()-1?", ":"");}
        cerr << "]: ";
        for(uint j=0,fl_size_j=h_facets[i].f_neighbors.size();j<fl_size_j;j++){
          cerr << h_facets[i].f_neighbors[j] << " ";
        }
        cerr << endl;
      }
    }
  }

  void ConvexHull::createInitializeSimplex() {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    //clear
    h_facets.clear();
    setCentroid();

    vector<FacetPoint> initial_extreme_points=getInitialExtremePoints();
    //fix h_reference now!
    //NOTE: originally we used the hull centroid to determine alignment of facet normals
    //however, this does not work in the case of 3D hulls with no 3D points (only 2D hull points)
    //instead, we align with the centroid of the initial_extreme_points, as this point SHOULD be guaranteed to be in the hull
    //the hull only gets bigger by accepting more points
    //so we set h_reference first to h_centroid, since we don't care about the alignment of the normal inside getInitialExtremePoints()
    //then we fix here and keep for all facets
    vector<xvector<double> > points;
    for(uint i=0,fl_size_i=initial_extreme_points.size();i<fl_size_i;i++){points.push_back(initial_extreme_points[i].ch_point.h_coords);}
    h_reference=aurostd::getCentroid(points);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " new h_reference=" << h_reference << endl;}

    //new first facet of endpoints

    if(initial_extreme_points.size()==h_dim+1){if(LDEBUG){cerr << __AFLOW_FUNC__ << " building initial simplex" << endl;}}
    else if(initial_extreme_points.size()==h_dim){if(LDEBUG){cerr << __AFLOW_FUNC__ << " not enough points to build a simplex, will settle for a facet instead" << endl;}}
    else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Not enough points to construct initial simplex/facet");}

    //get dim combos of points
    string error;
    aurostd::xcombos xc(initial_extreme_points.size(),h_dim,'C');
    vector<int> indices;
    while(xc.increment()){
      //remember, indices are NOT the indices themselves, but a place holder as to whether to include the number at the index
      //e.g., indices==1,1,0; therefore, include indices 0,1
      indices=xc.getCombo();
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " indices=";
        for(uint i=0,fl_size_i=indices.size();i<fl_size_i;i++){cerr << indices[i] << (i!=indices.size()-1?",":"");}
        cerr << endl;
      }
      h_facets.push_back(ChullFacet(*p_FileMESSAGE,*p_oss));  //CO20180305
      for(uint i=0,fl_size_i=indices.size();i<fl_size_i;i++){
        if(indices[i]==1){h_facets.back().addVertex(initial_extreme_points[i]);}
      }
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " initial facet[" << h_facets.size() << "] new point[" << h_facets.back().m_vertices.size() << "] ";
        cerr << "coords=" << m_points[h_facets.back().m_vertices.back().ch_index].h_coords << endl;
      }
      if(!h_facets.back().isValid(error)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"An initial facet cannot be created, invalid points input: "+error);}
      initializeFacet(h_facets.back(),false); //already validated
    }
    if(initial_extreme_points.size()==h_dim+1){if(LDEBUG){cerr << __AFLOW_FUNC__ << " initial simplex built" << endl;}}
    else if(initial_extreme_points.size()==h_dim){if(LDEBUG){cerr << __AFLOW_FUNC__ << " initial (single) facet built" << endl;}}
    else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Not enough points to construct initial simplex/facet");}

    //if net distance, then it must be outside hull
    bool associated;
    FacetPoint f_point;
    for(uint i=0,fl_size_i=h_points.size();i<fl_size_i;i++){
      associated=false;
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " associating point[" << h_points[i] << "] to a facet's outside set" << endl;}
      f_point.initialize(m_points[h_points[i]],h_points[i]);
      for(uint i_facet=0,fl_size_i_facet=h_facets.size();i_facet<fl_size_i_facet && !associated;i_facet++){
        if(h_facets[i_facet].isPointOnFacet(f_point)){associated=true;} //skip the obvious
        if(!associated && h_facets[i_facet].isPointOutside(f_point)){
          h_facets[i_facet].f_outside_set.push_back(f_point);
          associated=true;
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " associating point[" << h_points[i] << "] with facet[" << i_facet << "].f_outside_set" << endl;}
        }
      }
      if(!associated){
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " NOT associating point[" << h_points[i] << "] with any facet outside set" << endl;}
      }
    }
    //get furthest point
    for(uint i_facet=0,fl_size_i_facet=h_facets.size();i_facet<fl_size_i_facet;i_facet++){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " facet[" << i_facet << "] has " << h_facets[i_facet].f_outside_set.size() << " outside points" << endl;}
      h_facets[i_facet].setFurthestPoint();
    }
    setNeighbors();
  }

  void ConvexHull::setVisibleFacets(uint i_facet){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(i_facet>h_facets.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index with h_facets");}
    //clean
    for(uint i=0,fl_size_i=h_facets.size();i<fl_size_i;i++){h_facets[i].f_visited=false;}
    h_visible_facets.clear();

    //initialize
    h_visible_facets.push_back(i_facet);
    h_facets[i_facet].f_visited=true;

    uint i_visible,i_neigh;
    //[CO20190409 - size() must be recalculated for EVERY loop]for(uint i=0,fl_size_i=h_visible_facets.size();i<fl_size_i;i++)
    for(uint i=0;i<h_visible_facets.size();i++)
    { //CO20200106 - patching for auto-indenting
      i_visible=h_visible_facets[i];
      for(uint j=0,fl_size_j=h_facets[i_visible].f_neighbors.size();j<fl_size_j;j++){
        i_neigh=h_facets[i_visible].f_neighbors[j];
        if(h_facets[i_neigh].f_visited){continue;}
        h_facets[i_neigh].f_visited=true;
        if(h_facets[i_neigh].isPointOutside(h_facets[i_facet].f_furthest_point)){h_visible_facets.push_back(i_neigh);}
      }
      //there's no need to set i=0 for this loop (unlike the facets loop in calculateFacets()
      //we iterate through neighbors, visiting only once
      //the neighbors of a facet do not change here
    }

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " visible facets for facet[" << i_facet << "]: ";
      for(uint i=0,fl_size_i=h_visible_facets.size();i<fl_size_i;i++){cerr << h_visible_facets[i] << " ";}
      cerr << endl;
    }
  }

  void ConvexHull::setHorizonRidges(){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(!h_visible_facets.size()){return;}
    h_horizon_ridges.clear();
    bool match;
    for(uint i1=0,fl_size_i1=h_visible_facets.size();i1<fl_size_i1;i1++){
      for(uint j1=0,fl_size_j1=h_facets[h_visible_facets[i1]].m_ridges.size();j1<fl_size_j1;j1++){
        match=false;
        for(uint i2=0,fl_size_i2=h_visible_facets.size();i2<fl_size_i2&&!match;i2++){
          if(i1==i2){continue;}
          for(uint j2=0,fl_size_j2=h_facets[h_visible_facets[i2]].m_ridges.size();j2<fl_size_j2&&!match;j2++){
            if(h_facets[h_visible_facets[i1]].m_ridges[j1].getCHIndices()==h_facets[h_visible_facets[i2]].m_ridges[j2].getCHIndices()){
              match=true;
            }
          }
        }
        if(!match){h_horizon_ridges.push_back(h_facets[h_visible_facets[i1]].m_ridges[j1]);}
      }
    }

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " horizon ridges: ";
      for(uint i=0,fl_size_i=h_horizon_ridges.size();i<fl_size_i;i++){
        for(uint j=0,fl_size_j=h_horizon_ridges[i].m_vertices.size();j<fl_size_j;j++){
          cerr << h_horizon_ridges[i].m_vertices[j].ch_index << (j==h_horizon_ridges[i].m_vertices.size()-1?"":" ");
        }
        cerr << (i==h_horizon_ridges.size()-1?"":", ");
      }
      cerr << endl;
    }
  }

  uint ConvexHull::createNewFacets(FacetPoint furthest_point){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(furthest_point.ch_index>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index with furthest_point");}
    uint old_facet_count=h_facets.size();
    string error; //dummy so we don't recreate every time
    for(uint i=0,fl_size_i=h_horizon_ridges.size();i<fl_size_i;i++){
      //check obvious
      if(h_horizon_ridges[i].isPointOnFacet(furthest_point)){
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " skipping new facet: duplicate point" << endl;}
        continue;
      }
      h_horizon_ridges[i].addVertex(furthest_point);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " furthest point for new facet: " << furthest_point.ch_index << " " << m_points[furthest_point.ch_index].h_coords << endl;}
      if(!h_horizon_ridges[i].isValid(error)){  //corner case, happens when creating new facets on edges or corners of hull
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " skipping new facet: " << error << endl;}
        continue;
      }
      initializeFacet(h_horizon_ridges[i],false); //already validated
      h_facets.push_back(h_horizon_ridges[i]);
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " new facet: ";
        for(uint j=0,fl_size_j=h_facets.back().m_vertices.size();j<fl_size_j;j++){cerr << m_points[h_facets.back().m_vertices[j].ch_index].h_coords << " | ";}
        cerr << endl;
      }
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " " << h_facets.size()-old_facet_count << " new facets" << endl;}
    return h_facets.size()-old_facet_count;
    //do neighbors at the end!
  }

  void ConvexHull::updateOutsideSet(uint new_facet_count){  //they are at the end of the list
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    //we're deleting visible planes soon, so we need to reassess/reassign outside points
    uint i_visible;
    bool associated;
    for(uint i=0,fl_size_i=h_visible_facets.size();i<fl_size_i;i++){
      i_visible=h_visible_facets[i];
      for(uint j=0,fl_size_j=h_facets[i_visible].f_outside_set.size();j<fl_size_j;j++){
        FacetPoint& f_point=h_facets[i_visible].f_outside_set[j];
        associated=false;
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " re-associating point[" << f_point.ch_index << "] with a new facet" << endl;}
        for(uint fl_size_i_facet=h_facets.size(),i_facet=fl_size_i_facet-new_facet_count;i_facet<fl_size_i_facet&&!associated;i_facet++){
          if(h_facets[i_facet].isPointOnFacet(f_point)){associated=true;} //skip the obvious
          if(!associated && h_facets[i_facet].isPointOutside(f_point)){
            h_facets[i_facet].f_outside_set.push_back(f_point);
            associated=true;
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " associating point[" << f_point.ch_index << "] with facet[" << i_facet << "].f_outside_set" << endl;}
          }
        }
      }
    }
    for(uint fl_size_i_facet=h_facets.size(),i_facet=fl_size_i_facet-new_facet_count;i_facet<fl_size_i_facet;i_facet++){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " NEW facet[" << i_facet << "] has " << h_facets[i_facet].f_outside_set.size() << " outside points" << endl;}
      h_facets[i_facet].setFurthestPoint();
    }
  }

  void ConvexHull::deleteVisibleFacets() {
    //sort in descending order because all indices above the point of deletion change
    std::sort(h_visible_facets.rbegin(), h_visible_facets.rend());  //descending  //, std::greater<uint>());
    for(uint i=0,fl_size_i=h_visible_facets.size();i<fl_size_i;i++){h_facets.erase(h_facets.begin()+h_visible_facets[i]);}
  }

  void ConvexHull::removeDuplicateHullPoints() {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    std::sort(h_points.begin(),h_points.end());h_points.erase( std::unique( h_points.begin(), h_points.end() ), h_points.end() );  //first remove duplicate indices

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " starting" << endl;}
    vector<uint> indices_2_remove;
    for(uint i=0,fl_size_i=h_points.size();i<fl_size_i;i++){
      for(uint j=i+1,fl_size_j=h_points.size();j<fl_size_j;j++){
        if(identical(m_points[h_points[i]].h_coords,m_points[h_points[j]].h_coords,ZERO_TOL)){
          indices_2_remove.push_back(j);
        }
      }
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " removing points: ";
      for(uint i=0,fl_size_i=indices_2_remove.size();i<fl_size_i;i++){cerr << indices_2_remove[i] << " ";}
      cerr << endl;
    }
    std::sort(indices_2_remove.rbegin(),indices_2_remove.rend()); //descending
    for(uint i=0,fl_size_i=indices_2_remove.size();i<fl_size_i;i++){h_points.erase(h_points.begin()+indices_2_remove[i]);}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " done" << endl;}
  }

  void ConvexHull::calculateFacets() {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    //only copy those possibly on the hull (does not include endpoints)

    h_facets.clear();
    removeDuplicateHullPoints();
    if(!h_points.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No points found on hull");}
    h_dim=m_points[h_points[0]].h_coords.rows;

    for(uint i=1,fl_size_i=h_points.size();i<fl_size_i;i++){
      if((uint)m_points[h_points[i]].h_coords.rows!=h_dim){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid initialization of hull coordinates");}
    }

    if(LDEBUG){
      for(uint i=0,fl_size_i=h_points.size();i<fl_size_i;i++){
        cerr << __AFLOW_FUNC__ << " m_coords=" << m_points[h_points[i]].m_coords << "; h_coords=" << m_points[h_points[i]].h_coords << endl; //<< "; relevant " << m_points[h_points[i]].isRelevantPoint(elements_present) << endl;
      }
    }

    //main loop START - see algorithm here: 10.1145/235815.235821
    createInitializeSimplex();

    uint new_facet_count;
    //[CO20190409 - size() must be recalculated for EVERY loop]for(uint i_facet=0,fl_size_i_facet=h_facets.size();i_facet<fl_size_i_facet;i_facet++)
    for(uint i_facet=0;i_facet<h_facets.size();i_facet++)
    { //CO20200106 - patching for auto-indenting
      if(!h_facets[i_facet].f_outside_set.size()){continue;}  //the only way we exit this loop is if we continue for ALL facets
      setVisibleFacets(i_facet);
      setHorizonRidges();
      new_facet_count=createNewFacets(h_facets[i_facet].f_furthest_point);
      if(new_facet_count){
        updateOutsideSet(new_facet_count);  //redistribute outside set among new facets (visible planes will go soon)
        for(uint fl_size_i_facet=h_facets.size(),i_facet=fl_size_i_facet-new_facet_count;i_facet<fl_size_i_facet;i_facet++){h_facets[i_facet].setFurthestPoint();}  //get furthest points for new facets
        deleteVisibleFacets();  //delete visible facets
        setNeighbors();
      }
      i_facet=-1; //restart full loop!
    }
    //main loop END - see algorithm here: 10.1145/235815.235821

    //resort points in facets with knowledge of normal
    //if normal is pointed down, sort in descending order
    bool sort_stoich_ascending;
    for(uint i=0,fl_size_i=h_facets.size();i<fl_size_i;i++){
      sort_stoich_ascending=!h_facets[i].m_in_lower_hemisphere;
      std::sort(h_facets[i].m_vertices.begin(),h_facets[i].m_vertices.end(),
          sortThermoPoints(sort_stoich_ascending,m_sort_energy_ascending));
    }
    std::sort(h_facets.begin(),h_facets.end(),sortFacetsByPoints(m_points));  //auto sort

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " hull points:" << endl;
      for(uint i_facet=0,fl_size_i_facet=h_facets.size();i_facet<fl_size_i_facet;i_facet++){
        cerr << __AFLOW_FUNC__ << " facet " << i_facet << ": ";
        for(uint i=0,fl_size_i=h_facets[i_facet].m_vertices.size();i<fl_size_i;i++){
          cerr << m_points[h_facets[i_facet].m_vertices[i].ch_index].h_coords << " - auid " << m_points[h_facets[i_facet].m_vertices[i].ch_index].m_entry.auid << " | ";
        }
        cerr << "normal " << h_facets[i_facet].m_normal << " | angles " << aurostd::getGeneralAngles(h_facets[i_facet].m_normal,ZERO_TOL) << endl;//<< " | x " << cos(h_facets[i_facet].m_angle) << " | y " << sin(h_facets[i_facet].m_angle)  << endl;
      }
    }
  }

  const xvector<int>& ConvexHull::getElementsPresent(uint i_nary,uint i_alloy) const {
    if(i_nary>m_naries.size()-1 || i_alloy>m_naries[i_nary].m_alloys.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}
    return m_naries[i_nary].m_alloys[i_alloy].m_elements_present;
  }
  const xvector<int>& ConvexHull::getElementsPresent(uint i_point) const {
    if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,"ConvexHull::getElementsPresent():","Invalid index within points");}
    return m_points[i_point].m_elements_present;
  }
  xvector<int> ConvexHull::getElementsPresent(const vector<uint>& vcpoints) const {
    //get union of elements_present, just the sum
    if(vcpoints.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No points provided");}
    xvector<int> elements_present=getElementsPresent(vcpoints[0]);
    if(vcpoints.size()==1){return elements_present;}
    for(uint i=1,fl_size_i=vcpoints.size();i<fl_size_i;i++) {
      uint i_point=vcpoints[i];
      const xvector<int>& _elements_present=getElementsPresent(i_point);
      elements_present+=_elements_present;
    }
    return elements_present;
  }

  void ConvexHull::setElementsPresent(uint i_nary,uint i_alloy){m_elements_present=getElementsPresent(i_nary,i_alloy);}

  void ConvexHull::addRelevantUnariesToHullCalculation(uint i_nary,uint i_alloy) {
    if(i_nary>m_naries.size()-1 || i_alloy>m_naries[i_nary].m_alloys.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}
    return addRelevantUnariesToHullCalculation(m_naries[i_nary].m_alloys[i_alloy].m_elements_present);
  }

  void ConvexHull::addRelevantUnariesToHullCalculation(xvector<int>& elements_present) {
    uint i_nary=0;
    if((uint)elements_present.rows!=m_naries[i_nary].m_alloys.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unary alloy is missing from m_naries");}
    uint i_coord_group,i_point;
    for(int i_alloy=elements_present.lrows;i_alloy<=elements_present.urows;i_alloy++){
      if(elements_present[i_alloy]==1){
        if(m_naries[i_nary].m_alloys[i_alloy-elements_present.lrows].m_coord_groups.size()!=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unexpected count of coordgroups for unaries, should only be 1");}
        i_coord_group=m_naries[i_nary].m_alloys[i_alloy-elements_present.lrows].m_coord_groups[0];
        if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}
        for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_candidate_hull_points.size();i<fl_size_i;i++){
          i_point=m_coord_groups[i_coord_group].m_candidate_hull_points[i];
          if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
          addPointToHullCalculation(i_point,elements_present);
        }
      }
    }
  }

  void ConvexHull::addLowerDimensionPointsToHullCalculation(uint i_nary_max){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(i_nary_max>m_naries.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " starting" << endl;}
    //grab from m_naries[i_nary].m_alloys[i_alloy].m_facets
    //we already added unaries, so there will be duplicates, but we safely remove in calculateFacets()
    //don't worry about this yet
    uint i_point=AUROSTD_MAX_UINT,i_facet=AUROSTD_MAX_UINT;
    for(uint i_nary=1;i_nary<i_nary_max;i_nary++){
      if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
      for(uint i_alloy=0,fl_size_i_alloy=m_naries[i_nary].m_alloys.size();i_alloy<fl_size_i_alloy;i_alloy++){
        if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}
        if(m_naries[i_nary].m_alloys[i_alloy].belongs2Hull(m_elements_present)){
          for(uint i=0,fl_size_i=m_naries[i_nary].m_alloys[i_alloy].m_facets.size();i<fl_size_i;i++){
            i_facet=m_naries[i_nary].m_alloys[i_alloy].m_facets[i];
            if(i_facet>m_facets.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_facets");}
            const ChullFacet& facet=m_facets[i_facet];
            for(uint j=0,fl_size_j=facet.m_vertices.size();j<fl_size_j;j++){
              i_point=facet.m_vertices[j].ch_index;
              addPointToHullCalculation(i_point,m_elements_present); //this is hull specific, and set with setElementsPresent()
            }
          }
        }
      }
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " done" << endl;}
  }

  void ConvexHull::addPointToHullCalculation(uint i_point){
    if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    h_points.push_back(i_point);
    m_points[i_point].setHullCoords();
  }

  void ConvexHull::addPointToHullCalculation(uint i_point,xvector<int>& elements_present){
    if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    h_points.push_back(i_point);
    m_points[i_point].setHullCoords(elements_present);
  }

  void ConvexHull::preparePointsForHullCalculation(uint i_nary,uint i_alloy) {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(i_nary>m_naries.size()-1 || i_alloy>m_naries[i_nary].m_alloys.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " starting" << endl;}
    h_points.clear();
    addRelevantUnariesToHullCalculation(i_nary,i_alloy);
    addLowerDimensionPointsToHullCalculation(i_nary); //don't worry about adding unary duplicates, we remove them robustly in calculateFacets()
    uint i_coord_group,i_point;
    for(uint i=0,fl_size_i=m_naries[i_nary].m_alloys[i_alloy].m_coord_groups.size();i<fl_size_i;i++){
      i_coord_group=m_naries[i_nary].m_alloys[i_alloy].m_coord_groups[i];
      if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}
      for(uint j=0,fl_size_j=m_coord_groups[i_coord_group].m_candidate_hull_points.size();j<fl_size_j;j++){
        i_point=m_coord_groups[i_coord_group].m_candidate_hull_points[j];
        if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
        addPointToHullCalculation(i_point,m_naries[i_nary].m_alloys[i_alloy].m_elements_present);
      }
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " done" << endl;}
  }

  void ConvexHull::preparePointsForHullCalculation() {
    h_points.clear();
    //we already checked extremes in last (energy) direction, just add these to the hull
    //don't include those points in between
    uint i_point=AUROSTD_MAX_UINT;
    for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
      if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}
      for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_candidate_hull_points.size();i<fl_size_i;i++){
        i_point=m_coord_groups[i_coord_group].m_candidate_hull_points[i];
        if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
        addPointToHullCalculation(i_point);
      }
    }
  }

  uint ConvexHull::getNearestFacetVertically(const vector<uint>& i_facets,const ChullPoint& point) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(!point.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized point");}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_coords=" << point.m_coords << ", h_coords=" << point.h_coords << endl;}
    return getNearestFacetVertically(i_facets,point.h_coords);
  }

  uint ConvexHull::getNearestFacetVertically(const vector<uint>& i_facets,const xvector<double>& point) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    //determines nearness by vertical distance
    uint i_facet = 0, i_facet_min = 0;
    uint i_facet_artificial=-1; //really large uint
    double vdist,dist=AUROSTD_MAX_DOUBLE;
    for(uint i=0,fl_size_i=i_facets.size();i<fl_size_i;i++){
      i_facet=i_facets[i];
      if(i_facet>m_facets.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_facets");}
      const ChullFacet& facet=m_facets[i_facet];
      if(!facet.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Facet not initialized");}
      if(facet.m_is_hypercollinear){continue;}
      if(facet.m_is_vertical){continue;}
      if(facet.m_is_artificial){
        i_facet_artificial=i_facet;
        continue;
      }
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking at facet[" << i_facet << "]" << endl;}
      vdist=abs(facet.getSignedVerticalDistanceToZero(point));  //abs, //(repetita iuvant)
      if(vdist<dist){
        i_facet_min=i_facet;
        dist=vdist;
      }
    }
    //safety, return artificial facet if no other facets available
    if(dist==AUROSTD_MAX_DOUBLE){
      if(i_facet_artificial>m_facets.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No nearest facet found");}
      return i_facet_artificial;
    }
    return i_facet_min;
  }

  const vector<uint>& ConvexHull::getRelevantFacets(uint i_nary,uint i_alloy) const {
    if(i_nary>m_naries.size()-1 || i_alloy>m_naries[i_nary].m_alloys.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}

    if(i_nary==0){  //unaries are special case
      //try m_facets, otherwise start from top and go down
      if(m_i_facets.size()){return m_i_facets;}
      for(uint fl_size_i_nary=m_naries.size(),i_nary=fl_size_i_nary-1;i_nary<fl_size_i_nary;i_nary--){
        for(uint fl_size_i_alloy=m_naries[i_nary].m_alloys.size(),i_alloy=fl_size_i_alloy-1;i_alloy<=fl_size_i_alloy;i_alloy--){
          if(m_naries[i_nary].m_alloys[i_alloy].m_facets.size()){return m_naries[i_nary].m_alloys[i_alloy].m_facets;}
        }
      }
    }else{
      if(m_naries[i_nary].m_alloys[i_alloy].m_facets.size()){return m_naries[i_nary].m_alloys[i_alloy].m_facets;}
    }
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Facets not calculated");
  }

  void ConvexHull::setHullMembers() {return setHullMembers(m_i_facets);}
  void ConvexHull::setHullMembers(uint i_nary,uint i_alloy) {
    if(i_nary>m_naries.size()-1 || i_alloy>m_naries[i_nary].m_alloys.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}
    const vector<uint>& i_facets=getRelevantFacets(i_nary,i_alloy);//m_naries[i_nary].m_alloys[i_alloy].m_facets;
    return setHullMembers(i_facets);
  }

  void ConvexHull::setHullMembers(const vector<uint>& i_facets) {
    if(!i_facets.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Hull has yet to be calculated");}
    uint i_point=AUROSTD_MAX_UINT,g_state=AUROSTD_MAX_UINT,i_coord_group=AUROSTD_MAX_UINT,i_facet=AUROSTD_MAX_UINT;
    for(uint i=0,fl_size_i=i_facets.size();i<fl_size_i;i++){
      i_facet=i_facets[i];
      if(i_facet>m_facets.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_facets");}
      for(uint j=0,fl_size_j=m_facets[i_facet].m_vertices.size();j<fl_size_j;j++){
        i_point=m_facets[i_facet].m_vertices[j].ch_index;
        if(!getCoordGroupIndex(i_point,i_coord_group)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within coordgroups");}
        //hull member == yes
        m_points[i_point].m_is_on_hull=true;
        m_coord_groups[i_coord_group].m_is_on_hull=true;
        if(m_thermo_hull){m_coord_groups[i_coord_group].m_hull_member=i_point;}  //very safe, only one hull-member per coordgroup
        //g-state == if not artificial point
        g_state=artificialMap(i_point);
        if(!m_points[g_state].m_is_artificial){m_points[g_state].m_is_g_state=true;}
      }
    }
  }

  void ConvexHull::setNearestFacet(uint i_nary,uint i_alloy,uint i_coord_group){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(i_nary>m_naries.size()-1 || i_alloy>m_naries[i_nary].m_alloys.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}
    if(i_coord_group>m_coord_groups.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within coordgroups");}
    if(m_coord_groups[i_coord_group].m_is_on_hull){return;}
    if(m_coord_groups[i_coord_group].m_points.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup["+aurostd::utype2string(i_coord_group)+"] has no points");}
    if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}

    uint i_point=m_coord_groups[i_coord_group].m_points[0];
    if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
    ChullPoint& point=m_points[i_point];
    //it's possible that h_coords has not been set (we only do it to points that make it to hull calc)
    //so set again (repetita iuvant)
    xvector<int>& elements_present=m_naries[i_nary].m_alloys[i_alloy].m_elements_present;
    point.setHullCoords(elements_present);
    const vector<uint>& i_facets=getRelevantFacets(i_nary,i_alloy);//m_naries[i_nary].m_alloys[i_alloy].m_facets;
    uint i_facet=getNearestFacetVertically(i_facets,point);
    m_coord_groups[i_coord_group].m_nearest_facet=i_facet;
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " nearest_facet[i_nary=" << i_nary << ",i_alloy=" << i_alloy << ",i_coord_group=" << i_coord_group;
      cerr << "]=" << m_coord_groups[i_coord_group].m_nearest_facet << endl;
    }
  }

  double ConvexHull::getSignedVerticalDistanceWithinCoordGroup(uint i_coord_group,uint i_point) const {
    if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    return getSignedVerticalDistanceWithinCoordGroup(i_coord_group,m_points[i_point]);
  }

  double ConvexHull::getSignedVerticalDistanceWithinCoordGroup(uint i_coord_group,const ChullPoint& point) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}
    if(m_coord_groups[i_coord_group].m_points.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup["+aurostd::utype2string(i_coord_group)+"] has no points");}
    if(!point.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized point");}
    xvector<double> r_coords=point.getStoichiometricCoords();
    if(!identical(m_coord_groups[i_coord_group].m_coords,r_coords,ZERO_TOL)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point does not belong to coordgroup");}

    uint ref_point=artificialMap(m_coord_groups[i_coord_group].m_points[0]);  //don't need to set distance for artificial point, re-scale to this point
    if(!isViablePoint(ref_point)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup contains invalid point");}
    double ref_dist=0.0;
    //if(m_coord_groups[i_coord_group].m_is_on_hull){ //just a check //bad check, doesn't apply to unaries
    //if(!m_points[ref_point].m_is_on_hull){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup is on the hull, but hull point not found");}  //not the case for unaries
    //}else{  //[CO20200106 - close bracket for indenting]}
    if(!m_coord_groups[i_coord_group].m_is_on_hull){
      uint i_facet=m_coord_groups[i_coord_group].m_nearest_facet;
      if(i_facet>m_facets.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_facets");}
      ref_dist=m_facets[i_facet].getSignedVerticalDistance(m_points[ref_point]);  //signed distance
      //DO NOT CHECK SIGN OF VERTICAL DISTANCE, WE USE THIS FUNCTION FOR STABILITY CRITERION
      //if(m_half_hull){  //tests of stupidity
      //  //do not use signbit, add tol to zero (precision)
      //  //sign of distance:
      //  //independent of lower/upper hull:  above hull is negative, below hull is positive
      //  bool should_be_positive=!m_lower_hull;
      //  bool correct_sign_vertical_distance=chull::correctSignVerticalDistance(ref_dist,should_be_positive);
      //  if( m_lower_hull && !correct_sign_vertical_distance){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"(lower half hull) found point BELOW hull (entry="+m_points[ref_point].m_entry.auid+",dist2Hull="+aurostd::utype2string(ref_dist,4)+")");}
      //  if(!m_lower_hull && !correct_sign_vertical_distance){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"(upper half hull) found point ABOVE hull (entry="+m_points[ref_point].m_entry.auid+",dist2Hull="+aurostd::utype2string(ref_dist,4)+")");}
      //}
    }
    double delta=(point.getLastCoord()-m_points[ref_point].getLastCoord());
    double dist=ref_dist-delta;
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " point.m_entry.auid=" << point.m_entry.auid << endl;
      cerr << __AFLOW_FUNC__ << " m_points[ref_point].m_entry.auid=" << m_points[ref_point].m_entry.auid << endl;
      cerr << __AFLOW_FUNC__ << " point.m_entry.prototype=" << point.m_entry.prototype << endl;
      cerr << __AFLOW_FUNC__ << " m_points[ref_point].m_entry.prototype=" << m_points[ref_point].m_entry.prototype << endl;
      cerr << __AFLOW_FUNC__ << " point.m_coords=" << point.m_coords << endl;
      cerr << __AFLOW_FUNC__ << " m_points[ref_point].m_coords=" << m_points[ref_point].m_coords << endl;
      cerr << __AFLOW_FUNC__ << " ref_dist=" << ref_dist << endl;
      cerr << __AFLOW_FUNC__ << " delta=" << delta << endl;
      cerr << __AFLOW_FUNC__ << " dist=" << dist << endl;
    }
    return dist;
  }

  double ConvexHull::getDistanceToHull(uint i_point,bool redo,bool get_signed_distance) const {
    if(i_point>(m_points.size()-1)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");} //HE20210629 make comparison more precise
    return getDistanceToHull(m_points[i_point],redo,get_signed_distance);
  }

  double ConvexHull::getDistanceToHull(const ChullPoint& point,bool redo,bool get_signed_distance) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(!point.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized point");}
    //const vector<uint>& i_facets=m_i_facets;
    double dist;
    if(!redo){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " attempting to find coordgroup for fast distance calculation" << endl;}
      uint i_coord_group=AUROSTD_MAX_UINT;
      //only if the distances have already been found for all points, and we can identify the coordgroup (point may legitimately be outside)
      //the redo in the getCoordGroupIndex() is VERY safe, so let's keep it
      if(getCoordGroupIndex(point,i_coord_group,true)){
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " found coordgroup[" << i_coord_group << "] (coords=" << m_coord_groups[i_coord_group].m_coords << ")" << endl;}
        dist=getSignedVerticalDistanceWithinCoordGroup(i_coord_group,point);
        if(get_signed_distance==false){dist=abs(dist);} //CO20180828 //CO20190808
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " dist=" << dist << endl;}
        return dist;
      }
      //  uint ref_point=m_coord_groups[i_coord_group].m_ref_state;
      //  if(LDEBUG){cerr << __AFLOW_FUNC__ << " preliminarily setting ref_point=" << ref_point << endl;}
      //  double dist=AUROSTD_NAN;
      //  if(isViablePoint(ref_point)){dist=m_points[ref_point].getDist2Hull(_std_);}
      //  else{
      //    if(m_coord_groups[i_coord_group].m_is_on_hull){
      //      ref_point=m_coord_groups[i_coord_group].m_points[0];
      //      if(isViablePoint(ref_point) && m_points[ref_point].m_is_on_hull){dist=0.0;} //some safety checks, but ultimately set to 0
      //    }
      //  }
      //  if(LDEBUG){
      //    cerr << __AFLOW_FUNC__ << " m_points[ref_point=" << ref_point << "].m_coords=";
      //    for(int i=m_points[ref_point].m_coords.lrows;i<=m_points[ref_point].m_coords.urows;i++){cerr << m_points[ref_point].m_coords[i] << " ";}
      //    cerr << endl;
      //    cerr << __AFLOW_FUNC__ << " dist of ref_point = " << dist << endl;
      //  }
      //  if(!m_points[ref_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized point");}
      //  if(dist<AUROSTD_NAN){ //safety
      //    double delta=(m_points[ref_point].getLastCoord()-point.getLastCoord());
      //    if(LDEBUG){
      //      cerr << __AFLOW_FUNC__ << " delta = " << delta << endl;
      //      cerr << __AFLOW_FUNC__ << " dist+delta = " << dist+delta << endl;
      //      cerr << __AFLOW_FUNC__ << " dist-delta = " << dist-delta << endl;
      //      cerr << __AFLOW_FUNC__ << " delta-dist = " << delta-dist << endl;
      //    }
      //    if(LDEBUG){cerr << __AFLOW_FUNC__ << " getSignedVerticalDistanceWithinCoordGroup()=" << getSignedVerticalDistanceWithinCoordGroup(i_coord_group,point) << endl;}
      //    return delta-dist;
      //  }
    }
    vector<uint> i_facets=m_i_facets;
    //if we have stoich_coords, we need to find relevant facets
    if(m_has_stoich_coords){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " stoich_coords detected, trying to find relevant facets" << endl;}
      uint i_nary,i_alloy;
      if(!getAlloyIndex(point,i_nary,i_alloy)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot get alloys index");}
      if(i_nary==0){  //unaries are special cases, they RARELY work with normal getSignedVerticalDistance()
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " found unary, need to use coordgroup to find distance (rarely works with normal getSignedVerticalDistance()" << endl;}
        uint i_coord_group=AUROSTD_MAX_UINT;
        if(!getCoordGroupIndex(point,i_coord_group,true)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot find unary coordgroup");}
        dist=getSignedVerticalDistanceWithinCoordGroup(i_coord_group,point);
        if(get_signed_distance==false){dist=abs(dist);} //CO20180828 //CO20190808
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " dist=" << dist << endl;}
        return dist;
      }
      i_facets=getRelevantFacets(i_nary,i_alloy);
    }else{
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " searching through ALL facets for closest vertical distance" << endl;}
    }
    //otherwise, all facets are relevant (same dimension)
    uint i_facet=getNearestFacetVertically(i_facets,point);
    dist=m_facets[i_facet].getSignedVerticalDistance(point);
    if(get_signed_distance==false){dist=abs(dist);} //CO20180828 //CO20190808
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " dist=" << dist << endl;}
    return dist;
  }

  vector<double> ConvexHull::getDistancesToHull(const vector<string>& vauid,bool redo) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    stringstream message;
    vector<double> vdist2hull;
    uint i_point=AUROSTD_MAX_UINT;
    for(uint i=0,fl_size_i=vauid.size();i<fl_size_i;i++){
      const string& auid=vauid[i];
      if(auid.empty()){message << "Empty auid found";throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);}
      if(!findPoint(auid,i_point)){message << "Specified auid not found on hull (auid=" << auid << ")";throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " grabbing distance to hull for i_point=" << i_point << endl;}
      vdist2hull.push_back(getDistanceToHull(i_point,redo));
    }
    return vdist2hull;
  }

  bool ConvexHull::hullDistanceExtractionRequired() const {
    if(thermoPropertiesExtractionRequired()){return true;}  //dependency
    if(!m_cflags.flag("CHULL::MULTI_OUTPUT")&&(m_cflags.flag("CHULL::LATEX_DOC")||m_cflags.flag("CHULL::PNG_IMAGE"))){
      aurostd::xoption lflags=resolvePlotLabelSettings();
      bool meta_labels=lflags.flag("CHULL::PLOT_META_LABELS");
      bool labels_off_hull=lflags.flag("CHULL::PLOT_LABELS_OFF_HULL");
      if(meta_labels&&labels_off_hull){return true;}  //need distances to hull here
      bool no_doc=m_cflags.flag("CHULL::NO_DOC");
      if(no_doc){return false;}
      bool image_only=m_cflags.flag("CHULL::IMAGE_ONLY");
      if(image_only){return false;}
    }
    return true;
  }

  bool ConvexHull::thermoPropertiesExtractionRequired() const {
    if(thermoPostProcessingExtractionRequired()){return true;}  //dependency
    if(!m_cflags.flag("CHULL::MULTI_OUTPUT")&&(m_cflags.flag("CHULL::LATEX_DOC")||m_cflags.flag("CHULL::PNG_IMAGE"))){
      aurostd::xoption lflags=resolvePlotLabelSettings();
      bool icsd_labels=lflags.flag("CHULL::PLOT_ICSD_LABELS");
      if(icsd_labels){return true;}
      bool no_doc=m_cflags.flag("CHULL::NO_DOC");
      if(no_doc){return false;}
      bool image_only=m_cflags.flag("CHULL::IMAGE_ONLY");
      if(image_only){return false;}
    }
    return true;
  }

  bool ConvexHull::thermoPostProcessingExtractionRequired() const {
    if(!m_cflags.flag("CHULL::MULTI_OUTPUT")&&(m_cflags.flag("CHULL::LATEX_DOC")||m_cflags.flag("CHULL::PNG_IMAGE"))){
      bool no_doc=m_cflags.flag("CHULL::NO_DOC");
      if(no_doc){return false;}
      bool image_only=m_cflags.flag("CHULL::IMAGE_ONLY");
      if(image_only){return false;}
    }
    return true;
  }

  void ConvexHull::setDistancesToHull(uint i_nary,uint i_alloy) {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    stringstream message;
    if(!m_has_stoich_coords){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Non-stoich coordinates");}
    if(i_nary>m_naries.size()-1 || i_alloy>m_naries[i_nary].m_alloys.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}

    //SKIP_HULL_DISTANCE_CALCULATION is NOT recommended
    //the algorithm is already heavily optimized so it's quite fast
    bool perform_hull_distance_calculation=(1 &&
        hullDistanceExtractionRequired() &&
        (!m_cflags.flag("CHULL::SKIP_HULL_DISTANCE_CALCULATION")) &&
        TRUE); //(1&&!(m_cflags.flag("CHULL::SKIP_HULL_DISTANCE_CALCULATION")||(!m_cflags.flag("CHULL::MULTI_OUTPUT")&&m_cflags.flag("CHULL::LATEX_DOC")&&m_cflags.flag("CHULL::IMAGE_ONLY"))));
    if(!perform_hull_distance_calculation){return;}

    message << "Gathering hull distance data for";
    if(m_velements.size()){message << " " << aurostd::joinWDelimiter(alloyToElements(i_nary,i_alloy),"-");}
    message << " [i_nary=" << i_nary << ",i_alloy=" << i_alloy << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    uint i_coord_group=AUROSTD_MAX_UINT;
    for(uint i=0,fl_size_i=m_naries[i_nary].m_alloys[i_alloy].m_coord_groups.size();i<fl_size_i;i++){
      i_coord_group=m_naries[i_nary].m_alloys[i_alloy].m_coord_groups[i];
      if(!m_coord_groups[i_coord_group].m_points.size()){continue;}
      if(!m_coord_groups[i_coord_group].m_initialized){continue;}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking at i_coord_group=" << i_coord_group << endl;}
      if(!m_coord_groups[i_coord_group].m_is_on_hull){setNearestFacet(i_nary,i_alloy,i_coord_group);} //only off-hull need nearest facets
      setDistancesToHull(i_nary,i_alloy,i_coord_group);
    }
  }

  void ConvexHull::setDistancesToHull(uint i_nary,uint i_alloy,uint i_coord_group) {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(!m_has_stoich_coords){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Non-stoich coordinates");}
    if(i_nary>m_naries.size()-1 || i_alloy>m_naries[i_nary].m_alloys.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}
    if(i_coord_group>m_coord_groups.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within coordgroups");}
    if(m_coord_groups[i_coord_group].m_nearest_facet>m_facets.size()-1){setNearestFacet(i_nary,i_alloy,i_coord_group);}
    if(m_coord_groups[i_coord_group].m_points.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup["+aurostd::utype2string(i_coord_group)+"] has no points");}
    if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}

    uint i_point=m_coord_groups[i_coord_group].m_points[0];
    if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
    i_point=artificialMap(i_point); //don't need to set distance for artificial point, re-scale to this point
    ChullPoint& point=m_points[i_point];
    double dist_2_hull=0.0;
    bool check_dist_2_hull=false;
    if(!m_coord_groups[i_coord_group].m_is_on_hull){ //off-hull, this is the reason it works for unaries, unaries coordgroup are declared ON-HULL, so dist==0 automatically
      uint i_facet=m_coord_groups[i_coord_group].m_nearest_facet;
      dist_2_hull=m_facets[i_facet].getSignedVerticalDistance(point);  //do this "expensive" calculation only once, others are simply z-scaled
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " dist[coord=" << point.h_coords << "]=" << dist_2_hull << endl;
        if(check_dist_2_hull){
          cerr << __AFLOW_FUNC__ << " comparing fast dist calculation = "  << dist_2_hull;
          cerr << " vs. slow dist calculation =" << getDistanceToHull(point) << endl;
        }
      }
      if(m_half_hull){  //tests of stupidity
        //do not use signbit, add tol to zero (precision)
        //sign of distance:
        //independent of lower/upper hull:  above hull is negative, below hull is positive
        bool should_be_positive=!m_lower_hull;
        bool correct_sign_vertical_distance=chull::correctSignVerticalDistance(dist_2_hull,should_be_positive);
        if( m_lower_hull && !correct_sign_vertical_distance){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"(lower half hull) found point BELOW hull (entry="+point.m_entry.auid+",dist2Hull="+aurostd::utype2string(dist_2_hull,4)+")");}
        if(!m_lower_hull && !correct_sign_vertical_distance){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"(upper half hull) found point ABOVE hull (entry="+point.m_entry.auid+",dist2Hull="+aurostd::utype2string(dist_2_hull,4)+")");}
        //[OBSOLETE CO20180828]if(m_lower_hull){
        //[OBSOLETE CO20180828]  if(aurostd::notNegative(dist_2_hull,true,ZERO_TOL)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"(lower half hull) found point BELOW hull (entry="+point.m_entry.auid+",dist2Hull="+aurostd::utype2string(dist_2_hull,4)+")");}
        //[OBSOLETE CO20180828]}else{
        //[OBSOLETE CO20180828]  if(aurostd::notPositive(dist_2_hull,true,ZERO_TOL)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"(upper half hull) found point ABOVE hull (entry="+point.m_entry.auid+",dist2Hull="+aurostd::utype2string(dist_2_hull,4)+")");}
        //[OBSOLETE CO20180828]}
      }
    }
    //KEEP THE SIGN CONVENTION AS IS, it is correct without manipulation
    //by convention, distances inside hull are NEGATIVE, following the reaction equation (spontaneous decomposition)
    //which we just checked is the case, so force negative
    //dist_2_hull=-abs(dist_2_hull); //force positive, the value is what we care about now
    //[OBSOLETE CO20180828]point.m_dist_2_hull=m_coord_groups[i_coord_group].m_nearest_distance=dist_2_hull;
    point.m_dist_2_hull=m_coord_groups[i_coord_group].m_nearest_distance=abs(dist_2_hull);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " |dist|=" << point.m_dist_2_hull << endl;}
    uint j_point;
    for(uint j=0,fl_size_j=m_coord_groups[i_coord_group].m_points.size();j<fl_size_j;j++){
      j_point=m_coord_groups[i_coord_group].m_points[j];
      j_point=artificialMap(j_point); //don't need to set distance for artificial point, re-scale to this point
      m_points[j_point].m_dist_2_hull=point.m_dist_2_hull;
      //[OBSOLETE CO20180828]if(i_point!=j_point){m_points[j_point].m_dist_2_hull+=(point.getLastCoord()-m_points[j_point].getLastCoord());}  //simply ADD difference
      if(i_point!=j_point){m_points[j_point].m_dist_2_hull+=abs(point.getLastCoord() - m_points[j_point].getLastCoord());}  //simply ADD difference
      //it's possible that h_coords has not been set (we only do it to points that make it to hull calc)
      //so set again (repetita iuvant)
      if(LDEBUG){
        xvector<int>& elements_present=m_naries[i_nary].m_alloys[i_alloy].m_elements_present;
        m_points[j_point].setHullCoords(elements_present);
        cerr << __AFLOW_FUNC__ << " dist[coord=" << m_points[j_point].h_coords << "]=" << m_points[j_point].m_dist_2_hull << endl;
      }
    }
  }

  vector<uint> ConvexHull::extractDecompositionPhases(const ChullFacet& facet) const {
    uint i_point=AUROSTD_MAX_UINT;
    vector<uint> decomp_phases;
    for(uint i=0,fl_size_i=facet.m_vertices.size();i<fl_size_i;i++){
      i_point=facet.m_vertices[i].ch_index;
      decomp_phases.push_back(i_point);
    }
    if(decomp_phases.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No decomposition phases found");}
    std::sort(decomp_phases.begin(),
        decomp_phases.end(),
        sortCHullPoints(m_points,false,true)); //stoich sorting should be descending, energy sorting default okay here
    return decomp_phases;
  }

  vector<uint> ConvexHull::getDecompositionPhases(uint i_point) const {
    if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    return getDecompositionPhases(m_points[i_point]);
  }

  vector<uint> ConvexHull::getDecompositionPhases(const ChullPoint& point) const {
    if(!point.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized point");}
    if(point.m_is_on_hull){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No decomposition coefficients for hull members");}
    uint i_facet=AUROSTD_MAX_UINT;
    //polymorph possibility, other reductions in dimensionality are prescribed with the coefficients
    uint i_coord_group=AUROSTD_MAX_UINT;
    bool found_coord_group=getCoordGroupIndex(point,i_coord_group);
    if(found_coord_group){  //composition has already been considered by hull, might be g-state
      if(m_coord_groups[i_coord_group].m_is_on_hull){
        vector<uint> dcomp_phases;
        dcomp_phases.push_back(m_coord_groups[i_coord_group].m_ref_state);
        return dcomp_phases;
      }
      i_facet=m_coord_groups[i_coord_group].m_nearest_facet;
    }
    vector<uint> i_facets=m_i_facets;  //default, most broad  //make a copy so we can reset if m_has_stoich_coords, don't worry, only vector<uint> (pointers essentially)
    if(m_has_stoich_coords){
      uint i_nary,i_alloy;
      if(!getAlloyIndex(point,i_nary,i_alloy)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot get alloys index");}
      i_facets=getRelevantFacets(i_nary,i_alloy);
    }
    if(!found_coord_group||(i_facet>m_facets.size()-1)){i_facet=getNearestFacetVertically(i_facets,point);}
    return extractDecompositionPhases(m_facets[i_facet]);
  }

  void ConvexHull::setDecompositionPhases(uint i_nary,uint i_alloy,uint i_coord_group){
    if(i_nary>m_naries.size()-1 || i_alloy>m_naries[i_nary].m_alloys.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}
    if(i_coord_group>m_coord_groups.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within coordgroups");}
    if(m_coord_groups[i_coord_group].m_nearest_facet>m_facets.size()-1){setNearestFacet(i_nary,i_alloy,i_coord_group);}

    if(m_coord_groups[i_coord_group].m_is_on_hull){
      if(!isViablePoint(m_coord_groups[i_coord_group].m_hull_member)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No hull member set for m_coord_group["+aurostd::utype2string(i_coord_group)+"]",_RUNTIME_ERROR_);}
      vector<uint> dcomp_phases;dcomp_phases.push_back(m_coord_groups[i_coord_group].m_hull_member);
      m_coord_groups[i_coord_group].m_decomp_phases=dcomp_phases;
      return;
    }

    uint i_facet=m_coord_groups[i_coord_group].m_nearest_facet;
    ChullFacet& facet=m_facets[i_facet];
    m_coord_groups[i_coord_group].m_decomp_phases=extractDecompositionPhases(facet);
  }

  xvector<double> ConvexHull::getDecompositionCoefficients(uint i_point,vector_reduction_type vred) const {
    if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    return getDecompositionCoefficients(m_points[i_point],vred);
  }

  xvector<double> ConvexHull::getDecompositionCoefficients(const ChullPoint& point,vector_reduction_type vred) const {
    if(!point.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized point");}
    //[returns self]if(point.m_is_on_hull){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No decomposition coefficients for hull members");}
    uint i_coord_group=AUROSTD_MAX_UINT;
    bool found_coord_group=getCoordGroupIndex(point,i_coord_group);
    if(found_coord_group){  //composition has already been considered by hull, might be g-state
      if(m_coord_groups[i_coord_group].m_is_on_hull){
        xvector<double> dcomp_coefs(2);
        dcomp_coefs[dcomp_coefs.lrows]=dcomp_coefs[dcomp_coefs.lrows+1]=1.0;
        return dcomp_coefs;
      }
      if(m_coord_groups[i_coord_group].m_decomp_phases.size()){
        return getDecompositionCoefficients(point,m_coord_groups[i_coord_group].m_decomp_phases,vred);
      }
    }
    return getDecompositionCoefficients(point,getDecompositionPhases(point),vred);
  }

  xvector<double> ConvexHull::getDecompositionCoefficients(uint i_point,const vector<uint>& decomp_phases,vector_reduction_type vred) const {
    if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    return getDecompositionCoefficients(m_points[i_point],decomp_phases,vred);
  }

  xvector<double> ConvexHull::getDecompositionCoefficients(const ChullPoint& point,const vector<uint>& decomp_phases,vector_reduction_type vred) const {
    //do m_coords_group first (REDUCED)
    if(!point.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized point");}
    if(!decomp_phases.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No decomposition phases found");}
    if(decomp_phases.size()==1){  //probably polymorph
      uint i_point_other=decomp_phases[0];
      if(!isViablePoint(i_point_other)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Index to decomposition phase not found");}
      const ChullPoint& point_other=m_points[i_point_other];
      uint i_coord_group=AUROSTD_MAX_UINT,i_coord_group_other=AUROSTD_MAX_UINT;
      bool found_coord_group=getCoordGroupIndex(point,i_coord_group);
      bool found_coord_group_other=getCoordGroupIndex(point_other,i_coord_group_other);
      if(found_coord_group&&found_coord_group_other&&i_coord_group==i_coord_group_other){ //definitely polymorph
        xvector<double> dcomp_coefs(2);
        dcomp_coefs[dcomp_coefs.lrows]=dcomp_coefs[dcomp_coefs.lrows+1]=1.0;
        return dcomp_coefs;
      }
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Not enough decomposition phases provided");
    }
    vector<xvector<double> > lhs,rhs;
    //lhs
    lhs.push_back(point.getTruncatedReducedCoords(m_elements_present,vred));
    //rhs
    for(uint i=0,fl_size_i=decomp_phases.size();i<fl_size_i;i++){
      if(decomp_phases[i]>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
      if(!m_points[decomp_phases[i]].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized point");}
      rhs.push_back(m_points[decomp_phases[i]].getTruncatedReducedCoords(m_elements_present,vred));
    }
    xvector<double> coef=balanceChemicalEquation(lhs,rhs,true,ZERO_TOL);  //normalize

    return coef;
  }

  void ConvexHull::setDecompositionCoefficients(uint i_nary,uint i_alloy,uint i_coord_group){
    if(i_nary>m_naries.size()-1 || i_alloy>m_naries[i_nary].m_alloys.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}
    if(i_coord_group>m_coord_groups.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within coordgroups");}
    if(m_coord_groups[i_coord_group].m_decomp_phases.size()==0){setDecompositionPhases(i_nary,i_alloy,i_coord_group);}

    if(m_coord_groups[i_coord_group].m_is_on_hull){
      xvector<double> dcomp_coefs(2);
      dcomp_coefs[dcomp_coefs.lrows]=dcomp_coefs[dcomp_coefs.lrows+1]=1.0;
      m_coord_groups[i_coord_group].m_decomp_coefs=dcomp_coefs;
      return;
    }

    //we get different coefficients between stoich and composition
    //ALWAYS use composition (even POCC, simply won't reduce), and do NOT mix stoich + composition
    uint i_point=m_coord_groups[i_coord_group].m_points[0];
    if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
    vector<uint>& decomp_phases=m_coord_groups[i_coord_group].m_decomp_phases;
    m_coord_groups[i_coord_group].m_decomp_coefs=getDecompositionCoefficients(i_point,decomp_phases,frac_vrt);

    //[OBSOLETE - reduce by frac_vrt always! so use coord_group values]//now do individual coefficients //[OBSOLETE - reduce by frac always], NO REDUCE
    //[OBSOLETE - reduce by frac_vrt always! so use coord_group values]for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
    //[OBSOLETE - reduce by frac_vrt always! so use coord_group values]  i_point=m_coord_groups[i_coord_group].m_points[i];
    //[OBSOLETE - reduce by frac_vrt always! so use coord_group values]  if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
    //[OBSOLETE - reduce by frac_vrt always! so use coord_group values]  m_points[i_point].m_decomp_coefs=getDecompositionCoefficients(i_point,decomp_phases,frac_vrt);
    //[OBSOLETE - reduce by frac_vrt always! so use coord_group values]}
  }
  void ConvexHull::setOffHullProperties(uint i_nary,uint i_alloy){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    stringstream message;
    if(!m_has_stoich_coords){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Non-stoich coordinates");}
    if(i_nary>m_naries.size()-1 || i_alloy>m_naries[i_nary].m_alloys.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}

    message << "Gathering decomposition reaction data for";
    if(m_velements.size()){message << " " << aurostd::joinWDelimiter(alloyToElements(i_nary,i_alloy),"-");}
    message << " [i_nary=" << i_nary << ",i_alloy=" << i_alloy << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    uint i_coord_group=AUROSTD_MAX_UINT;
    for(uint i=0,fl_size_i=m_naries[i_nary].m_alloys[i_alloy].m_coord_groups.size();i<fl_size_i;i++){
      i_coord_group=m_naries[i_nary].m_alloys[i_alloy].m_coord_groups[i];
      if(!m_coord_groups[i_coord_group].m_points.size()){continue;}
      if(!m_coord_groups[i_coord_group].m_initialized){continue;}
      if(m_coord_groups[i_coord_group].m_is_on_hull){continue;}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking at i_coord_group=" << i_coord_group << endl;}
      //setNearestFacet(i_nary,i_alloy,i_coord_group);    //we now do separately and previously (setDistancesToHull())
      //setDistancesToHull(i_nary,i_alloy,i_coord_group); //we now do separately and previously (setDistancesToHull())
      setDecompositionPhases(i_nary,i_alloy,i_coord_group);
      setDecompositionCoefficients(i_nary,i_alloy,i_coord_group);
      setEquivalentStates(i_nary,i_alloy,i_coord_group);
    }
  }

  vector<uint> ConvexHull::getAdjacentFacets(uint hull_member,bool ignore_hypercollinear,bool ignore_vertical,bool ignore_artificial) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(hull_member>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No hull member has been identified");}
    vector<uint> adjacent_i_facets;
    if(!m_points[hull_member].m_is_on_hull){return adjacent_i_facets;}

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " starting" << endl;}

    uint i_nary=m_points[hull_member].m_i_nary;
    uint i_alloy=m_points[hull_member].m_i_alloy;
    if(i_nary>m_naries.size()-1 || i_alloy>m_naries[i_nary].m_alloys.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}
    const vector<uint>& i_facets=getRelevantFacets(i_nary,i_alloy);//m_naries[i_nary].m_alloys[i_alloy].m_facets;
    uint i_facet;
    for(uint i=0,fl_size_i=i_facets.size();i<fl_size_i;i++){
      i_facet=i_facets[i];
      if(i_facet>m_facets.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_facets");}
      if(ignore_hypercollinear && m_facets[i_facet].m_is_hypercollinear){continue;}  //really (d-1) facets
      if(ignore_vertical && m_facets[i_facet].m_is_vertical){continue;}           //really (d-1) facets
      if(ignore_artificial && m_facets[i_facet].m_is_artificial){continue;}       //contains all unaries
      if(m_facets[i_facet].isPointOnFacet(hull_member)){adjacent_i_facets.push_back(i_facet);}
    }
    if(adjacent_i_facets.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No adjacent facets found");}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " stop" << endl;}
    return adjacent_i_facets;
  }

  vector<vector<uint> > ConvexHull::getEquilibriumPhases(uint hull_member) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(hull_member>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    const ChullPoint& point=m_points[hull_member];
    if(!point.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized point");}
    if(!point.m_is_on_hull){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No equilibrium phases for non hull members");}
    vector<vector<uint> > equilibrium_phases;

    vector<uint> adjacent_i_facets=getAdjacentFacets(hull_member,true,true,true);
    //make full copy of relevant facets so we can sort accordingly, ignore vertical and artificial facets
    vector<ChullFacet> adjacent_facets;
    for(uint i=0,fl_size_i=adjacent_i_facets.size();i<fl_size_i;i++){adjacent_facets.push_back(m_facets[adjacent_i_facets[i]]);}

    //sort
    for(uint i_facet=0,fl_size_i_facet=adjacent_facets.size();i_facet<fl_size_i_facet;i_facet++){
      std::sort(adjacent_facets[i_facet].m_vertices.begin(),adjacent_facets[i_facet].m_vertices.end(),
          sortThermoPoints(false,true)); //stoich sorting should be descending, energy sorting default okay here
    }
    std::sort(adjacent_facets.begin(),adjacent_facets.end(),
        sortFacetsByPoints(m_points,false,false,false,false));  //not auto sort, sort everything descending

    uint i_point=AUROSTD_MAX_UINT;
    for(uint i_facet=0,fl_size_i_facet=adjacent_facets.size();i_facet<fl_size_i_facet;i_facet++){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " found new relevant facet[" << i_facet << "] with points: ";}
      equilibrium_phases.push_back(vector<uint>(0));
      for(uint i=0,fl_size_i=adjacent_facets[i_facet].m_vertices.size();i<fl_size_i;i++){
        i_point=adjacent_facets[i_facet].m_vertices[i].ch_index;
        equilibrium_phases.back().push_back(i_point);
        if(LDEBUG){cerr << m_points[i_point].h_coords << " [" << i_point << "]" << (i!=adjacent_facets[i_facet].m_vertices.size()-1?", ":"");}
      }
      if(LDEBUG){cerr << endl;}
    }
    return equilibrium_phases;
  }

  void ConvexHull::setEquilibriumPhases(uint i_nary,uint i_alloy,uint i_coord_group){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(i_nary>m_naries.size()-1 || i_alloy>m_naries[i_nary].m_alloys.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}
    if(i_coord_group>m_coord_groups.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within coordgroups");}
    if(!m_coord_groups[i_coord_group].m_is_on_hull){return;}
    if(m_coord_groups[i_coord_group].m_points.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup["+aurostd::utype2string(i_coord_group)+"] has no points");}
    if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}

    //get hull-member and check that it's correct
    uint hull_member=m_coord_groups[i_coord_group].m_hull_member;
    if(hull_member>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No hull member has been identified");}
    if(!m_points[hull_member].m_is_on_hull){return;}
    if(m_points[hull_member].isUnary()){return;} //others are in equilibrium with it, must be context of mixture!

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking at hull-member[" << hull_member << "]=" << m_points[hull_member].h_coords << endl;}

    m_coord_groups[i_coord_group].m_equilibrium_phases=getEquilibriumPhases(hull_member);

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " stop" << endl;}
  }

  bool ConvexHull::phasesEquivalent(uint i_point1,uint i_point2,bool perform_structure_comparison) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(i_point1>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    if(i_point2>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    if(!m_points[i_point1].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point1)+"] is not initialized");}
    if(!m_points[i_point2].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point2)+"] is not initialized");}

    if(i_point1==i_point2){return true;} //it is equivalent to self, so we only do checks if necessary

    //return false early for points without an entry
    if(!m_points[i_point1].m_has_entry){return false;}
    if(!m_points[i_point2].m_has_entry){return false;}

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " comparing [auid=" << m_points[i_point1].m_entry.auid << "] and [auid=" << m_points[i_point2].m_entry.auid << "]" << endl;}

    //strict==false, there might be entries in the database without a formation_enthalpy or sg calculated
    //ultimately, we want to compare structures
    if(energiesDiffer(i_point1,i_point2,false)){return false;} //first filter by those with wildly different energies, not strict
    if(spacegroupsDiffer(i_point1,i_point2,false)){return false;} //first filter by those with wildly different spacegroups, not strict
    if(!perform_structure_comparison){return false;}  //only return true if you can compare structures
    if(!structuresEquivalent(i_point1,i_point2)){return false;}

    return true;
  }

  //strict === strictly differ
  //if we don't know, because of AUROSTD_NAN's or NOSG's, we may still return false anyway to 
  //continue on to more strict determination later
  bool ConvexHull::energiesDiffer(uint i_point1,uint i_point2,bool strict) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(i_point1>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    if(i_point2>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}

    if(i_point1==i_point2){return false;}
    if(!m_points[i_point1].m_has_entry){return strict;}  //if strict, then return that they do differ (true)
    if(!m_points[i_point2].m_has_entry){return strict;}  //if strict, then return that they do differ (true)

    const double& energy1=m_points[i_point1].getFormationEnthalpy();
    const double& energy2=m_points[i_point2].getFormationEnthalpy();
    if(LDEBUG){cerr << __AFLOW_FUNC__ << energy1 << " vs. " << energy2 << endl;}
    if(energy1>=AUROSTD_NAN || energy2>=AUROSTD_NAN){return strict;} //if strict, then return that they do differ (true)

    bool differs=aurostd::isdifferent(energy1,energy2,ENERGY_TOL);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " energies do " << (differs?"":"NOT ") << "differ" << endl;}
    return differs;
  }

  //strict === strictly differ
  //if we don't know, because of AUROSTD_NAN's or NOSG's, we may still return false anyway to 
  //continue on to more strict determination later
  bool ConvexHull::spacegroupsDiffer(uint i_point1,uint i_point2,bool strict) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(i_point1>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    if(i_point2>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}

    if(i_point1==i_point2){return false;}
    if(!m_points[i_point1].m_has_entry){return strict;}  //if strict, then return that they do differ (true)
    if(!m_points[i_point2].m_has_entry){return strict;}  //if strict, then return that they do differ (true)

    const string& sg1=m_points[i_point1].getSG();
    const string& sg2=m_points[i_point2].getSG();
    if(LDEBUG){cerr << __AFLOW_FUNC__ << sg1 << " vs. " << sg2 << endl;}
    if(sg1==NOSG || sg2==NOSG){return strict;}  //if strict, then return that they do differ (true)

    bool differs=(sg1!=sg2);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " spacegroups do " << (differs?"":"NOT ") << "differ" << endl;}
    return differs;
  }

  bool ConvexHull::structuresEquivalent(uint i_point1,uint i_point2) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(i_point1>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    if(i_point2>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}

    if(i_point1==i_point2){return true;}
    //we will not be so stringent on throwing errors here, if we cannot find/load structure, simply return false
    if(!m_points[i_point1].m_has_entry){return false;}
    if(!m_points[i_point2].m_has_entry){return false;}

    xstructure a,b;
    if(!m_points[i_point1].getMostRelaxedXstructure(a)){return false;}
    if(!m_points[i_point2].getMostRelaxedXstructure(b)){return false;}
    //if(!m_points[i_point1].loadXstructures(true)){return false;}
    //if(!m_points[i_point2].loadXstructures(true)){return false;}

    //const xstructure& a=m_points[i_point1].m_entry.vstr[0];
    //const xstructure& b=m_points[i_point2].m_entry.vstr[0];
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " loaded both structures" << endl;
      cerr << __AFLOW_FUNC__ << " structure 1" << endl;
      cerr << a;
      cerr << __AFLOW_FUNC__ << " structure 2" << endl;
      cerr << b;
    }
    bool are_equivalent=(compare::structuresMatch(a,b,true,false,false)||compare::structuresMatch(b,a,true,false,false)); //match species and use fast match, but not scale volume, two structures with different volumes (pressures) are different! //DX20180123 - added fast_match = true //DX20190318 - not fast_match but optimized_match=false //CO20221112 - bit OR short circuit, there are some structures where compare_materials is input-specific
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " structures are " << (are_equivalent?"":"NOT ") << "equivalent" << endl;}
    return are_equivalent;
  }

  vector<uint> ConvexHull::getEquivalentEntries(uint i_point) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    const ChullPoint& point=m_points[i_point];
    if(!point.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized point");}
    uint i_coord_group=AUROSTD_MAX_UINT;
    if(!getCoordGroupIndex(point,i_coord_group)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within coordgroups");}
    if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}
    bool perform_structure_comparison=(1&&(!m_cflags.flag("CHULL::SKIP_STRUCTURE_COMPARISON"))); //(1&&!(m_cflags.flag("CHULL::SKIP_STRUCTURE_COMPARISON")||(!m_cflags.flag("CHULL::MULTI_OUTPUT")&&m_cflags.flag("CHULL::LATEX_DOC")&&m_cflags.flag("CHULL::IMAGE_ONLY"))));

    vector<uint> equivalent_entries;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " starting" << endl;}
    uint j_point=AUROSTD_MAX_UINT;
    for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
      j_point=m_coord_groups[i_coord_group].m_points[i];
      if(!phasesEquivalent(i_point,j_point,perform_structure_comparison)){continue;}
      equivalent_entries.push_back(j_point);
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " g-state[" << i_point << "]=" << m_points[i_point].h_coords;
      cerr << " equivalent structures=";
      for(uint i=0,fl_size_i=equivalent_entries.size();i<fl_size_i;i++){cerr << equivalent_entries[i] << (i!=equivalent_entries.size()-1?", ":"");}
      cerr << endl;
    }
    return equivalent_entries;
  }

  bool ConvexHull::isICSD(uint i_point,uint& i_point_icsd) const {
    //this function is important because ICSD_XX is a prototype and may not actually correspond to compound A1B1_ICSD_XX
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    i_point_icsd=AUROSTD_MAX_UINT;
    if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    string point_icsd_number=getICSDNumber(i_point,true);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking for icsd_number=" << point_icsd_number << endl;}
    if(point_icsd_number.empty()){return false;}  //no possible match
    if(m_icsd_entries.size()==0){return false;} //no possible match
    string entry_icsd_number="";
    for(uint i=0,fl_size_i=m_icsd_entries.size();i<fl_size_i;i++){
      if(i_point==m_icsd_entries[i]){i_point_icsd=m_icsd_entries[i];return true;} //instant match
      entry_icsd_number=getICSDNumber(m_icsd_entries[i],true);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " comparing with entry_icsd_number=" << entry_icsd_number << endl;}
      if(entry_icsd_number.empty()){continue;}  //this is weird/bad, but let's not exit the whole program
      //we need to get the STRING number (may have 0's in the front) and species/stoichiometry equivalent
      //this can be done by looking at elements_present, since they map to species
      if(point_icsd_number==entry_icsd_number && m_points[i_point].m_elements_present==m_points[m_icsd_entries[i]].m_elements_present){
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " found match: point[i_point=" << i_point << "].prototype=" << m_points[i_point].m_entry.prototype;
          cerr << " vs. point[i_point=" << m_icsd_entries[i] << "].prototype=" << m_points[m_icsd_entries[i]].m_entry.prototype << endl;
        }
        i_point_icsd=m_icsd_entries[i];
        return true;
      }
    }
    return false;
  }

  void ConvexHull::setEquivalentGStates(uint i_nary,uint i_alloy,uint i_coord_group){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    m_coord_groups[i_coord_group].m_equivalent_g_states.clear();
    m_coord_groups[i_coord_group].m_calculated_equivalent_g_states=false;
    if(i_nary>m_naries.size()-1 || i_alloy>m_naries[i_nary].m_alloys.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}
    if(i_coord_group>m_coord_groups.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within coordgroups");}
    if(!m_coord_groups[i_coord_group].m_is_on_hull){return;}
    if(m_coord_groups[i_coord_group].m_points.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup["+aurostd::utype2string(i_coord_group)+"] has no points");}
    if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}

    //we need to do a structure comparison, so get artificial map (did already before)
    uint g_state=m_coord_groups[i_coord_group].m_ref_state;
    if(g_state>m_points.size()-1){return;}
    if(!m_points[g_state].m_has_entry){return;} //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No entry (structure) found"); //only point in coordgroup

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking at g-state[" << g_state << "]=" << m_points[g_state].h_coords << endl;}
    m_coord_groups[i_coord_group].m_equivalent_g_states=getEquivalentEntries(g_state);
    if(!m_coord_groups[i_coord_group].m_equivalent_g_states.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No equivalent states found (not even self)");}
    m_coord_groups[i_coord_group].m_calculated_equivalent_g_states=true;
    uint i_point=AUROSTD_MAX_UINT,i_point_icsd=AUROSTD_MAX_UINT;
    bool found_icsd=false,_found_icsd=false;  //whether icsd exists among g_states
    uint i_canonical_icsd=AUROSTD_MAX_UINT;       //index of canonical icsd
    uint canonical_icsd_number=AUROSTD_MAX_UINT;  //lowest number
    uint icsd_number=AUROSTD_MAX_UINT;
    for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_equivalent_g_states.size();i<fl_size_i;i++){
      i_point=m_coord_groups[i_coord_group].m_equivalent_g_states[i];
      if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
      m_points[i_point].m_is_equivalent_g_state=true; //g_state and equivalent_g_state should stay separate
      m_points[i_point].m_calculated_equivalent_entries=true;
      m_points[i_point].m_equivalent_entries=m_coord_groups[i_coord_group].m_equivalent_g_states;
      //get icsd information
      //we were strict in creating m_icsd_entries (only _ICSD_), but we need to check LIBs for ICSD_
      if(aurostd::substring2bool(m_points[i_point].m_entry.prototype,"ICSD_")){  //remember, _ICSD_ means icsd catalog, but ICSD_ means LIBs
        _found_icsd=isICSD(i_point,i_point_icsd);
        found_icsd=(found_icsd||_found_icsd); //short circuit bit OR
        if(_found_icsd){  //find lowest icsd number
          icsd_number=aurostd::string2utype<uint>(getICSDNumber(i_point,true));
          if(icsd_number<canonical_icsd_number){
            i_canonical_icsd=i_point_icsd;
            canonical_icsd_number=icsd_number;
          }
        }
        //aurostd::string2tokens(m_points[i_point].m_entry.prototype,tokens,"_");
        //if(tokens.size()==3){ //skip everything else
        //  icsd_number=aurostd::string2utype<uint>(tokens[2]);
        //  if(icsd_number<canonical_icsd_number){
        //    i_canonical_icsd=i_point;
        //    canonical_icsd_number=icsd_number;
        //  }
        //}
      }
    }

    m_coord_groups[i_coord_group].m_found_icsd_g_state=found_icsd;
    m_coord_groups[i_coord_group].m_i_icsd_g_state=i_canonical_icsd;
    
    for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_equivalent_g_states.size();i<fl_size_i;i++){
      i_point=m_coord_groups[i_coord_group].m_equivalent_g_states[i];
      if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
      m_points[i_point].m_found_icsd=found_icsd;
      m_points[i_point].m_i_icsd=i_canonical_icsd;
    }

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " stop" << endl;}
  }
  
  void ConvexHull::setEquivalentStates(uint i_nary,uint i_alloy,uint i_coord_group){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(i_nary>m_naries.size()-1 || i_alloy>m_naries[i_nary].m_alloys.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}
    if(i_coord_group>m_coord_groups.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within coordgroups");}
    if(m_coord_groups[i_coord_group].m_points.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup["+aurostd::utype2string(i_coord_group)+"] has no points");}
    if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}
    stringstream message;

    //clear all
    if(m_coord_groups[i_coord_group].m_is_on_hull){
      m_coord_groups[i_coord_group].m_equivalent_g_states.clear();
      m_coord_groups[i_coord_group].m_calculated_equivalent_g_states=false;
      m_coord_groups[i_coord_group].m_found_icsd_g_state=false;
      m_coord_groups[i_coord_group].m_i_icsd_g_state=AUROSTD_MAX_UINT;
    }
    uint i_point=AUROSTD_MAX_UINT;
    for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
      i_point=m_coord_groups[i_coord_group].m_points[i];
      if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
      if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
      m_points[i_point].m_equivalent_entries.clear();
      m_points[i_point].m_calculated_equivalent_entries=false;
      m_points[i_point].m_found_icsd=false;
      m_points[i_point].m_i_icsd=AUROSTD_MAX_UINT;
    }
     
    //get vpoints2check
    vector<uint> vpoints2check;
    if(m_coord_groups[i_coord_group].m_is_on_hull){
      //we need to do a structure comparison, so get artificial map (did already before)
      uint g_state=m_coord_groups[i_coord_group].m_ref_state;
      vpoints2check.push_back(g_state);
    }
    double dist=0.0;
    if(m_cflags.flag("CHULL::DIST_FOR_EQUIVALENCE_ANALYSIS")){dist=aurostd::string2utype<double>(m_cflags.getattachedscheme("CHULL::DIST_FOR_EQUIVALENCE_ANALYSIS"));}
    bool ND_only=m_cflags.flag("CHULL::DIST_FOR_EQUIVALENCE_ANALYSIS_ND_ONLY");
    if((!(ND_only && i_nary!=(m_dim-1))) && aurostd::nonZeroWithinTol(dist,ZERO_TOL)){
      message << "Performing equivalence analysis for points within " << aurostd::utype2string(dist) << " (meV/atom) for the";
      if(m_velements.size()){message << " " << aurostd::joinWDelimiter(alloyToElements(i_nary,i_alloy),"-") << " hull";}
      message << " [i_nary=" << i_nary << ",i_alloy=" << i_alloy << ",i_coord_group=" << i_coord_group << "]";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
        i_point=m_coord_groups[i_coord_group].m_points[i];
        if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
        if(aurostd::WithinList(vpoints2check,i_point)){continue;}
        if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
        if(m_points[i_point].getDist2Hull(_m_)<=dist){vpoints2check.push_back(i_point);}
      }
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " vpoints2check=" << aurostd::joinWDelimiter(vpoints2check,",") << endl;}

    if(vpoints2check.empty()){return;}
    
    uint j_point=AUROSTD_MAX_UINT,j_point_icsd=AUROSTD_MAX_UINT;
    //
    bool found_icsd=false,_found_icsd=false;  //whether icsd exists among equivalent entries
    uint i_canonical_icsd=AUROSTD_MAX_UINT;       //index of canonical icsd
    uint canonical_icsd_number=AUROSTD_MAX_UINT;  //lowest number
    uint icsd_number=AUROSTD_MAX_UINT;
    
    for(uint i=0,fl_size_i=vpoints2check.size();i<fl_size_i;i++){
      i_point=vpoints2check[i];
      if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
      ChullPoint& point_i=m_points[i_point];
      if(!point_i.m_has_entry){continue;} //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No entry (structure) found"); //only point in coordgroup
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking at i_point[" << i_point << "]=" << m_points[i_point].h_coords << endl;}
      if(point_i.m_calculated_equivalent_entries==true){continue;}  //speed up
      point_i.m_equivalent_entries=getEquivalentEntries(i_point);
      if(!point_i.m_equivalent_entries.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No equivalent states found (not even self)");}
      point_i.m_calculated_equivalent_entries=true;
      if(m_coord_groups[i_coord_group].m_is_on_hull && i_point==m_coord_groups[i_coord_group].m_ref_state){ //is on hull
        m_coord_groups[i_coord_group].m_equivalent_g_states=point_i.m_equivalent_entries;
        m_coord_groups[i_coord_group].m_calculated_equivalent_g_states=true;
      }
      j_point_icsd=AUROSTD_MAX_UINT;
      found_icsd=false,_found_icsd=false;  //whether icsd exists among equivalent entries
      i_canonical_icsd=AUROSTD_MAX_UINT;       //index of canonical icsd
      canonical_icsd_number=AUROSTD_MAX_UINT;  //lowest number
      icsd_number=AUROSTD_MAX_UINT;
      for(uint j=0,fl_size_j=point_i.m_equivalent_entries.size();j<fl_size_j;j++){
        j_point=point_i.m_equivalent_entries[j];
        if(j_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
        if(!m_points[j_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(j_point)+"] is not initialized");}
        ChullPoint& point_j=m_points[j_point];
        //[CO20221112 - should not be necessary, might mess up consistency of data, just let it redo]if(point_j.m_calculated_equivalent_entries==true){continue;}  //speed up
        if(i_point!=j_point){
          point_j.m_equivalent_entries=point_i.m_equivalent_entries;
          point_j.m_calculated_equivalent_entries=true;
        }
        if(m_coord_groups[i_coord_group].m_is_on_hull && i_point==m_coord_groups[i_coord_group].m_ref_state){ //is on hull
          point_j.m_is_equivalent_g_state=true; //g_state and equivalent_g_state should stay separate
        }
        //get icsd information
        //we were strict in creating m_icsd_entries (only _ICSD_), but we need to check LIBs for ICSD_
        if(aurostd::substring2bool(point_j.m_entry.prototype,"ICSD_")){  //remember, _ICSD_ means icsd catalog, but ICSD_ means LIBs
          _found_icsd=isICSD(j_point,j_point_icsd);
          found_icsd=(found_icsd||_found_icsd); //short circuit bit OR
          if(_found_icsd){  //find lowest icsd number
            icsd_number=aurostd::string2utype<uint>(getICSDNumber(j_point,true));
            if(icsd_number<canonical_icsd_number){
              i_canonical_icsd=j_point_icsd;
              canonical_icsd_number=icsd_number;
            }
          }
        }
      }
      for(uint j=0,fl_size_j=point_i.m_equivalent_entries.size();j<fl_size_j;j++){
        j_point=point_i.m_equivalent_entries[j];
        if(!m_points[j_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(j_point)+"] is not initialized");}
        ChullPoint& point_j=m_points[j_point];
        point_j.m_found_icsd=found_icsd;
        point_j.m_i_icsd=i_canonical_icsd;
      }
      if(m_coord_groups[i_coord_group].m_is_on_hull && i_point==m_coord_groups[i_coord_group].m_ref_state){ //is on hull  //if(m_coord_groups[i_coord_group].m_is_on_hull){
        m_coord_groups[i_coord_group].m_found_icsd_g_state=found_icsd;
        m_coord_groups[i_coord_group].m_i_icsd_g_state=i_canonical_icsd;
      }
    }

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " stop" << endl;}
  }

  vector<uint> ConvexHull::getSymEquivalentGStates(uint g_state) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(g_state>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    const ChullPoint& point=m_points[g_state];
    if(!point.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized point");}
    if(!point.m_is_g_state){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No equivalent g-states for non ground-state structures");}
    uint i_coord_group=AUROSTD_MAX_UINT;
    if(!getCoordGroupIndex(point,i_coord_group)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within coordgroups");}

    vector<uint> sym_equivalent_g_states;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " starting" << endl;}
    uint i_point=AUROSTD_MAX_UINT;
    for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
      i_point=m_coord_groups[i_coord_group].m_points[i];
      if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
      if(g_state!=i_point){if(spacegroupsDiffer(g_state,i_point,true)){continue;}}
      sym_equivalent_g_states.push_back(i_point);
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " g-state[" << g_state << "]=" << m_points[g_state].h_coords;
      cerr << " symmetrically equivalent structures=";
      for(uint i=0,fl_size_i=sym_equivalent_g_states.size();i<fl_size_i;i++){cerr << sym_equivalent_g_states[i] << (i!=sym_equivalent_g_states.size()-1?", ":"");}
      cerr << endl;
    }
    return sym_equivalent_g_states;
  }

  void ConvexHull::setSymEquivalentGStates(uint i_nary,uint i_alloy,uint i_coord_group){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    m_coord_groups[i_coord_group].m_sym_equivalent_g_states.clear();
    if(i_nary>m_naries.size()-1 || i_alloy>m_naries[i_nary].m_alloys.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}
    if(i_coord_group>m_coord_groups.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within coordgroups");}
    if(!m_coord_groups[i_coord_group].m_is_on_hull){return;}
    if(m_coord_groups[i_coord_group].m_points.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup["+aurostd::utype2string(i_coord_group)+"] has no points");}
    if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}

    //we need to do a structure comparison, so get artificial map (did already before)
    uint g_state=m_coord_groups[i_coord_group].m_ref_state;
    if(g_state>m_points.size()-1){return;}
    if(!m_points[g_state].m_has_entry){return;} //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No entry (structure) found"); //only point in coordgroup

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking at g-state[" << g_state << "]=" << m_points[g_state].h_coords << endl;}
    m_coord_groups[i_coord_group].m_sym_equivalent_g_states=getSymEquivalentGStates(g_state);
    uint i_point=AUROSTD_MAX_UINT;
    for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_sym_equivalent_g_states.size();i<fl_size_i;i++){
      i_point=m_coord_groups[i_coord_group].m_sym_equivalent_g_states[i];
      if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
      m_points[i_point].m_is_sym_equivalent_g_state=true; //g_state and sym_equivalent_g_state should stay separate
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " stop" << endl;}
  }

  void ConvexHull::setOnHullProperties(uint i_nary,uint i_alloy){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    stringstream message;
    if(!m_has_stoich_coords){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Non-stoich coordinates");}
    if(i_nary>m_naries.size()-1 || i_alloy>m_naries[i_nary].m_alloys.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}

    message << "Gathering equilibrium phases and determining equivalent ground-states for";
    if(m_velements.size()){message << " " << aurostd::joinWDelimiter(alloyToElements(i_nary,i_alloy),"-");}
    message << " [i_nary=" << i_nary << ",i_alloy=" << i_alloy << "]";
    message << ", please be patient";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    uint i_coord_group=AUROSTD_MAX_UINT;
    for(uint i=0,fl_size_i=m_naries[i_nary].m_alloys[i_alloy].m_coord_groups.size();i<fl_size_i;i++){
      i_coord_group=m_naries[i_nary].m_alloys[i_alloy].m_coord_groups[i];
      if(!m_coord_groups[i_coord_group].m_points.size()){continue;}
      if(!m_coord_groups[i_coord_group].m_initialized){continue;}
      if(!m_coord_groups[i_coord_group].m_is_on_hull){continue;}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking at i_coord_group=" << i_coord_group << endl;}
      //setDistancesToHull(i_nary,i_alloy,i_coord_group); //we now do separately and previously (setDistancesToHull())
      //very important that you do not simply go through all facet points and find equilibrium points
      //this will overwrite binary information with ternary information
      //proceed safely with i_coord_group's
      setDecompositionPhases(i_nary,i_alloy,i_coord_group);
      setDecompositionCoefficients(i_nary,i_alloy,i_coord_group);
      setEquilibriumPhases(i_nary,i_alloy,i_coord_group);
      setSymEquivalentGStates(i_nary,i_alloy,i_coord_group);
      setEquivalentStates(i_nary,i_alloy,i_coord_group);
    }
  }

  void ConvexHull::storeHullData(uint i_nary,uint i_alloy){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    stringstream message;
    if(!h_facets.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Hull has yet to be calculated");}
    if(i_nary>m_naries.size()-1 || i_alloy>m_naries[i_nary].m_alloys.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}

    m_naries[i_nary].m_alloys[i_alloy].m_facets.clear();
    for(uint i=0,fl_size_i=h_facets.size();i<fl_size_i;i++){
      m_facets.push_back(h_facets[i]);
      m_naries[i_nary].m_alloys[i_alloy].m_facets.push_back(m_facets.size()-1);
    }
    setHullMembers(i_nary,i_alloy);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " hull points for i_nary=" << i_nary << ",i_alloy=" << i_alloy << endl;
      for(uint i=0,fl_size_i=m_points.size();i<fl_size_i;i++){
        if(m_points[i].m_is_on_hull){
          cerr << __AFLOW_FUNC__ << " hull_point[" << i << "].m_coords=" << m_points[i].m_coords << endl;
        }
      }
    }
    message << "Hull properties stored for";
    if(m_velements.size()){message << " " << aurostd::joinWDelimiter(alloyToElements(i_nary,i_alloy),"-");}
    message << " [i_nary=" << i_nary << ",i_alloy=" << i_alloy << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  void ConvexHull::storeHullData() {
    for(uint i=0,fl_size_i=h_facets.size();i<fl_size_i;i++){
      m_facets.push_back(h_facets[i]);
      m_i_facets.push_back(m_facets.size()-1);
    }
    setHullMembers();
  }

  void ConvexHull::extractThermodynamicProperties(uint i_nary,uint i_alloy){
    stringstream message;
    if(!h_facets.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Hull has yet to be calculated");}
    if(!m_has_stoich_coords){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Non-stoich coordinates");}
    if(i_nary>m_naries.size()-1 || i_alloy>m_naries[i_nary].m_alloys.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within m_naries");}
    if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary");}
    if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized alloy");}

    bool perform_thermo_properties_extraction=(1 && 
        thermoPropertiesExtractionRequired() && 
        (!m_cflags.flag("CHULL::SKIP_HULL_DISTANCE_CALCULATION")) &&
        (!m_cflags.flag("CHULL::SKIP_THERMO_PROPERTIES_EXTRACTION")) &&
        TRUE); //(1&&!(m_cflags.flag("CHULL::SKIP_THERMO_PROPERTIES_EXTRACTION")||(!m_cflags.flag("CHULL::MULTI_OUTPUT")&&m_cflags.flag("CHULL::LATEX_DOC")&&m_cflags.flag("CHULL::IMAGE_ONLY"))));
    if(!perform_thermo_properties_extraction){return;}
    if(!m_thermo_hull){
      message << "Cannot extract thermodynamic properties, thermodynamic hull NOT detected";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
      return;
    }

    message << "Thermodynamic hull detected, gathering on/off hull properties for";
    if(m_velements.size()){message << " " << aurostd::joinWDelimiter(alloyToElements(i_nary,i_alloy),"-");}
    message << " [i_nary=" << i_nary << ",i_alloy=" << i_alloy << "]";
    if(i_nary==0){message << " (last)";}
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    setOffHullProperties(i_nary,i_alloy);
    setOnHullProperties(i_nary,i_alloy);

    message << "Thermodynamic properties calculated for";
    if(m_velements.size()){message << " " << aurostd::joinWDelimiter(alloyToElements(i_nary,i_alloy),"-");}
    message << " [i_nary=" << i_nary << ",i_alloy=" << i_alloy << "]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  void ConvexHull::thermodynamicsPostProcessing(){
    //hull must be initialized for these analyses
    stringstream message;
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Hull not initialized");}
    if(!m_thermo_hull){return;}

    //do not perform thermodynamics postprocessing without having done extractThermodynamicProperties
    //we need equivalent g_states (setOnHullProperties())
    bool perform_thermo_post_processing=(1 &&
        thermoPostProcessingExtractionRequired() &&
        (!m_cflags.flag("CHULL::SKIP_HULL_DISTANCE_CALCULATION")) &&
        (!m_cflags.flag("CHULL::SKIP_THERMO_PROPERTIES_EXTRACTION")) &&
        (!m_cflags.flag("CHULL::SKIP_THERMO_POST_PROCESSING")) &&
        TRUE); //(1&&!(m_cflags.flag("CHULL::SKIP_THERMO_PROPERTIES_EXTRACTION")||(!m_cflags.flag("CHULL::MULTI_OUTPUT")&&m_cflags.flag("CHULL::LATEX_DOC")&&m_cflags.flag("CHULL::IMAGE_ONLY"))));
    if(!perform_thermo_post_processing){return;}

    message << "Performing thermodynamics post-processing";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    bool perform_stability_criterion=(1&&(!m_cflags.flag("CHULL::SKIP_STABILITY_CRITERION_ANALYSIS"))); //(1&&!(m_cflags.flag("CHULL::SKIP_STABILITY_CRITERION_ANALYSIS")||(!m_cflags.flag("CHULL::MULTI_OUTPUT")&&m_cflags.flag("CHULL::LATEX_DOC")&&m_cflags.flag("CHULL::IMAGE_ONLY"))));
    if(perform_stability_criterion){setStabilityCriterion();}

    bool perform_n_plus_1_enthalpy_gain=(1&&(!m_cflags.flag("CHULL::SKIP_N+1_ENTHALPY_GAIN_ANALYSIS"))); //(1&&!(m_cflags.flag("CHULL::SKIP_N+1_ENTHALPY_GAIN_ANALYSIS")||(!m_cflags.flag("CHULL::MULTI_OUTPUT")&&m_cflags.flag("CHULL::LATEX_DOC")&&m_cflags.flag("CHULL::IMAGE_ONLY"))));
    if(perform_n_plus_1_enthalpy_gain){setNPlus1EnthalpyGain();}

    return;
  }

  void ConvexHull::calculate(){
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    stringstream message;

    for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
      if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}
      for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
        m_points[m_coord_groups[i_coord_group].m_points[i]].m_is_on_hull=false; //completely refresh
      }
      m_coord_groups[i_coord_group].m_is_on_hull=false; //completely refresh
    }

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " starting" << endl;}
    //we first run through alloy hulls IF stoich_coords, grabbing hull_members
    if(m_has_stoich_coords){
      if(!m_naries.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Points have yet to be structured");}
      bool calc_highest_hull_only=m_cflags.flag("CHULL::CALCULATE_HIGHEST_DIMENSION_ONLY");
      uint i_nary_start=1;  //start at binaries
      if(calc_highest_hull_only){
        i_nary_start=m_naries.size()-1;
        message << "Calculating the highest dimensional hull ONLY (stoichiometric coordinates)";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      }else{
        message << "Calculating the hull(s) in increasing dimensionality (stoichiometric coordinates)";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      }
      for(uint i_nary=i_nary_start,fl_size_i_nary=m_naries.size();i_nary<fl_size_i_nary;i_nary++){
        for(uint i_alloy=0,fl_size_i_alloy=m_naries[i_nary].m_alloys.size();i_alloy<fl_size_i_alloy;i_alloy++){
          message << "Calculating " << pflow::arity_string(i_nary+1,false,false) << " hull for";
          if(m_velements.size()){message << " " << aurostd::joinWDelimiter(alloyToElements(i_nary,i_alloy),"-");}
          message << " [i_nary=" << i_nary << ",i_alloy=" << i_alloy << "]";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
          cleanHull();
          setElementsPresent(i_nary,i_alloy); //m_has_stoich_coords only
          if(calc_highest_hull_only){
            preparePointsForHullCalculation();  //inject all points
          }else{
            preparePointsForHullCalculation(i_nary,i_alloy);  //will have unary duplicates, but don't worry, we remove in calculateFacets()
          }
          calculateFacets();
          message << pflow::arity_string(i_nary+1,true,false) << " hull calculated for";
          if(m_velements.size()){message << " " << aurostd::joinWDelimiter(alloyToElements(i_nary,i_alloy),"-");}
          message << " [i_nary=" << i_nary << ",i_alloy=" << i_alloy << "]";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
          storeHullData(i_nary,i_alloy);
          setDistancesToHull(i_nary,i_alloy); //do this always
          if(m_thermo_hull){extractThermodynamicProperties(i_nary,i_alloy);}
        }
      }
      //set unary properties last
      //sort of a hack, but actually an intelligent solution
      //unary hull-members will be identified in higher dimensions
      //therefore, we simply need to calculate their properties
      for(uint i_alloy=0,fl_size_i_alloy=m_naries[0].m_alloys.size();i_alloy<fl_size_i_alloy;i_alloy++){setDistancesToHull(0,i_alloy);} //do this always
      if(m_thermo_hull){  //very safe, not sure how these algorithms perform outside of this domain, already know we have m_has_stoich_coords
        //message << "Calculating thermodynamic properties for " << pflow::arity_string(1,false,true) << " (last)";
        //pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
        for(uint i_alloy=0,fl_size_i_alloy=m_naries[0].m_alloys.size();i_alloy<fl_size_i_alloy;i_alloy++){extractThermodynamicProperties(0,i_alloy);} //knows to skip facet storage //setOnHullProperties(0,i_alloy);
      }
      //print counts of g_states by dimension
      uint gstate_count=0,total_gstate_count=0;
      for(uint i_nary=0,fl_size_i_nary=m_naries.size();i_nary<fl_size_i_nary;i_nary++){
        gstate_count=getGStateCount(i_nary);
        message << "Found " << gstate_count << " " << pflow::arity_string(i_nary+1,false,false) << " ";
        if(m_thermo_hull){message << "ground-state phases";}
        else{message << "points on the hull";}
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
        total_gstate_count+=gstate_count;
      }
      message << "Found " << total_gstate_count << " ";
      if(m_thermo_hull){message << "ground-state phases";}
      else{message << "points on the hull";}
      message << " total";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }else{
      message << "Entering default convex hull calculation (full-dimensional)";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      cleanHull();
      preparePointsForHullCalculation();
      calculateFacets();
      message << "Hull calculated";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      //print count of gstates (hull points)
      for(uint i_nary=0,fl_size_i_nary=m_naries.size();i_nary<fl_size_i_nary;i_nary++){
        message << "Found " << getGStateCount() << " points on the hull total";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      }
    }
    storeHullData();  //that way, the largest dim facets get stored in m_facets!

    //resort within coord_groups to expose hull-members!
    for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
      std::sort(m_coord_groups[i_coord_group].m_points.begin(),m_coord_groups[i_coord_group].m_points.end(),sortWithinCoordGroup(m_points,m_sort_energy_ascending));  //ascending order
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " done" << endl;}
  }

  double ConvexHull::getStabilityCriterion(const string& cauid) const {
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Hull not initialized");}
    vector<string> vcauid; vcauid.push_back(cauid);
    vector<double> vsc=getStabilityCriterion(vcauid);
    return vsc[0];
  }

  double ConvexHull::getStabilityCriterion(uint cpoint) const {
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Hull not initialized");}
    if(cpoint>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    vector<uint> vcpoints; vcpoints.push_back(cpoint);
    vector<double> vsc=getStabilityCriterion(vcpoints);
    return vsc[0];
  }

  vector<double> ConvexHull::getStabilityCriterion(const vector<string>& vcauid) const {
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Hull not initialized");}
    stringstream message;

    uint i_point=AUROSTD_MAX_UINT;
    vector<uint> vcpoints; //follows vcauid
    for(uint i=0,fl_size_i=vcauid.size();i<fl_size_i;i++){
      const string& auid=vcauid[i];
      if(auid.empty()){message << "Empty auid found";throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);}
      if(!findPoint(auid,i_point)){message << "Specified auid not found on hull (auid=" << auid << ")";throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);}
      vcpoints.push_back(i_point);
    }
    return getStabilityCriterion(vcpoints);
  }

  vector<double> ConvexHull::getStabilityCriterion(const vector<uint>& vcpoints) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Hull not initialized");}
    stringstream message;

    //set fake_hull
    xvector<int> elements_present_points=getElementsPresent(vcpoints);
    xvector<int> elements_present_hull=elements_present_points;
    if(sum(elements_present_hull)==1){
      //pick ANY other index to be 1 as well
      if(elements_present_hull[elements_present_hull.lrows]==0){elements_present_hull[elements_present_hull.lrows]=1;}
      else{
        if(elements_present_hull.rows<2){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"elements_present.rows<2",_INPUT_ILLEGAL_);} //we cannot have unary hulls as input (facets issue)
        if(elements_present_hull[elements_present_hull.lrows+1]==1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"elements_present mismatch for unary hull",_INDEX_MISMATCH_);}
        elements_present_hull[elements_present_hull.lrows+1]=1;
      }
    }
    if(sum(elements_present_hull)==1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Could not fix sum(elements_present_hull)==1",_RUNTIME_ERROR_);}
    ConvexHull fake_hull;
    getFakeHull(vcpoints,elements_present_hull,elements_present_points,fake_hull);

    //since i_nary and i_alloy don't change, the getDistanceToHull function should work fine (getRelevantFacets())
    vector<double> vscriterion;
    ChullPoint point;
    uint i_point=AUROSTD_MAX_UINT;
    for(uint i=0,fl_size_i=vcpoints.size();i<fl_size_i;i++) {
      i_point=vcpoints[i];
      //already checked validity of these points
      point=m_points[i_point];
      point.cleanPointForHullTransfer();  //clean now
      point.reduceCoords(elements_present_hull);  //reduce to minimum hull necessary to calculate scriterion (for unaries this will be a binary)
      point.setHullCoords(); //set to most general coords (m_coords), this reflects relevantFacets()
      //we need to clean this point of any properties relevant to REAL hull (not to be confused with fake hull)
      //clean first, then set h_coords
      //in the case that it is not a unary, leverage i_alloy of REAL hull first, then clean
      //if(point.isUnary()){
      //  point.setHullCoords(); //set to most general coords (m_coords), this reflects relevantFacets()
      //}else{
      //  if(!getAlloyIndex(point,i_nary,i_alloy)){message << "Alloy index not set (auid=" << point.m_entry.auid << ")";throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);}
      //  xvector<int> elements_present=m_naries[i_nary].m_alloys[i_alloy].m_elements_present;
      //  if(LDEBUG){cerr << __AFLOW_FUNC__ << " elements_present=" << elements_present << endl;}
      //  point.setHullCoords(elements_present);  //just to be sure
      //}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " m.coords=" << point.m_coords << endl;}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " h.coords=" << point.h_coords << endl;}
      vscriterion.push_back(fake_hull.getDistanceToHull(point,false));
    }
    return vscriterion;
  }

  void ConvexHull::getFakeHull(const vector<uint>& vcpoints,ConvexHull& fake_hull) const {
    if(!m_points.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No points loaded, no way to determine dimensionality of hull");}
    xvector<int> elements_present=aurostd::ones_xv<int>(m_points[0].m_coords.lrows,m_points[0].m_coords.lrows+(int)getDim()-1);
    return getFakeHull(vcpoints,elements_present,elements_present,fake_hull);
  }
  void ConvexHull::getFakeHull(const vector<uint>& vcpoints,const xvector<int>& elements_present,ConvexHull& fake_hull) const {
    return getFakeHull(vcpoints,elements_present,elements_present,fake_hull);
  }
  void ConvexHull::getFakeHull(const vector<uint>& vcpoints,const xvector<int>& elements_present_hull,const xvector<int>& elements_present_points,ConvexHull& fake_hull) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Hull not initialized");}
    stringstream message;

    //get full set of points to neglect
    uint i_point=AUROSTD_MAX_UINT,i_coord_group=AUROSTD_MAX_UINT,g_state=AUROSTD_MAX_UINT;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " initial points_to_neglect.size()=" << vcpoints.size() << endl;}

    //vector<ChullPoint> points_to_neglect; //not sure why it needs to be vector of chullpoints, simply hold indices
    vector<uint> points_to_neglect;
    vector<uint> eq_gstates;
    for(uint i=0,fl_size_i=vcpoints.size();i<fl_size_i;i++){
      i_point=vcpoints[i];
      if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
      const ChullPoint& point=m_points[i_point];  //this point may not be on the hull, it may be an equivalent structure, but coordgroup is
      if(!point.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized chullpoint");}
      if(!getCoordGroupIndex(point,i_coord_group)){message << "Coordgroup index not set (input[" << i << "])";throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);}
      if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}
      if(!m_coord_groups[i_coord_group].m_points.size()){message << "No points found within coordgroup (input[" << i << "])";throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);}
      //if(!m_coord_groups[i_coord_group].m_is_on_hull){message << "Coordgroup not on the hull (input[" << i << "])";throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);} //obsolete, now we don't assume coordgroup is on hull, we can remove ANY point
      points_to_neglect.push_back(i_point);
      if(m_coord_groups[i_coord_group].m_is_on_hull){
        g_state=m_coord_groups[i_coord_group].m_ref_state;  //we lose i_point here, so make sure to add it before this if-statement, otherwise we get SUPER degeneracy
        if(!m_points[g_state].m_is_g_state){message << "Coordgroup reference is not a ground-state structure (input[" << i << "])";throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);}
        //const vector<uint>& eq_gstates=m_coord_groups[i_coord_group].m_equivalent_g_states;
        if(m_coord_groups[i_coord_group].m_calculated_equivalent_g_states){eq_gstates=m_coord_groups[i_coord_group].m_equivalent_g_states;}
        else{eq_gstates=getEquivalentEntries(g_state);}
        if(!eq_gstates.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No equivalent states found (not even self)");}
        //for(uint j=0,fl_size_j=eq_gstates.size();j<fl_size_j;j++){points_to_neglect.push_back(m_points[eq_gstates[j]]);}
        for(uint j=0,fl_size_j=eq_gstates.size();j<fl_size_j;j++){points_to_neglect.push_back(eq_gstates[j]);}
      }
    }

    std::sort(points_to_neglect.begin(),points_to_neglect.end());points_to_neglect.erase( std::unique( points_to_neglect.begin(), points_to_neglect.end() ), points_to_neglect.end() );  //first remove duplicate indices

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " points_to_neglect.size()=" << points_to_neglect.size() << endl;
      cerr << __AFLOW_FUNC__ << " points_to_neglect: ";
      for(uint i=0,fl_size_i=points_to_neglect.size();i<fl_size_i;i++){cerr << points_to_neglect[i] << " ";} 
      cerr << endl;
    }

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " elements_present_hull=" << elements_present_hull << endl;
      cerr << __AFLOW_FUNC__ << " elements_present_points=" << elements_present_points << endl;
    }
    if(sum(elements_present_points)>sum(elements_present_hull)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"elements_present_hull must belong to bigger subspace than elements_present_points",_INPUT_ILLEGAL_);}

    //create new set of points
    //bool found=false;
    vector<ChullPoint> new_points;
    //do NOT go through m_points, this may include some duplicates we previously excluded (now not duplicates since we are removing points: AlFe hull)
    for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_coord_groups[i_coord_group=" << i_coord_group << "].getElementsPresent()=" << m_coord_groups[i_coord_group].getElementsPresent() << endl;}
      if(!subspaceBelongs(elements_present_points,m_coord_groups[i_coord_group].getElementsPresent())){
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " Ignoring irrelevant subspace of m_coord_groups[i_coord_group=" << i_coord_group << "] (m_coord_groups[i_coord_group=" << i_coord_group << "].getElementsPresent()=" << m_coord_groups[i_coord_group].getElementsPresent() << ")" << endl;}
        continue;
      }
      for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
        i_point=m_coord_groups[i_coord_group].m_points[i];
        if(aurostd::WithinList(points_to_neglect,i_point)){
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " Removing point from fake hull (i_point=" << i_point << ",auid=" << m_points[i_point].m_entry.auid << ")" << endl;}
          continue;
        }
        const ChullPoint& point=m_points[i_point];
        if(point.m_is_artificial){continue;}  //they will be added again
        //found=false;
        //for(uint j=0,fl_size_j=points_to_neglect.size();j<fl_size_j&&!found;j++){
        //  if(point.m_entry.auid==points_to_neglect[j].m_entry.auid){
        //    if(LDEBUG){cerr << __AFLOW_FUNC__ << " Removing (equivalent) ground-state (auid=" << points_to_neglect[j].m_entry.auid << ")" << endl;}
        //    found=true;
        //  }
        //}
        //if(found){continue;}
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " including point[i_point=" << i_point << ",auid=" << m_points[i_point].m_entry.auid << "].compound=" << point.m_entry.compound << endl;}
        new_points.push_back(point);
        new_points.back().reduceCoords(elements_present_hull);
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " point[" << i_point << "].m_coords=" << new_points.back().m_coords << endl;}
        new_points.back().cleanPointForHullTransfer();
      }
    }
    if(!new_points.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No points found for pseudo convex hull");}

    //patch m_velements to match reduceCoords(elements_present_hull)
    vector<string> velements;
    if(m_velements.size()>0){
      vector<uint> relevant_indices=getRelevantIndices(elements_present_hull);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " relevant_indices=" << aurostd::joinWDelimiter(relevant_indices,",") << endl;}
      for(uint i=0,fl_size_i=relevant_indices.size();i<fl_size_i;i++){velements.push_back(m_velements[relevant_indices[i]-1]);} //-1 because xvector->vecto
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " velements=" << aurostd::joinWDelimiter(velements,",") << endl;}
    }

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " Creating new (pseudo) hull without relevant g-states" << endl;}
    aurostd::xoption cflags=m_cflags;
    //[do NOT do, we need hull distance calculation]cflags.flag("CHULL::SKIP_HULL_DISTANCE_CALCULATION",true);  //no need to get the distances for ALL points, we will do for the one after hull creation
    cflags.flag("CHULL::SKIP_THERMO_PROPERTIES_EXTRACTION",true); //thermo properties NOT needed, just need hull
    //[killed off by SKIP_THERMO_PROPERTIES_EXTRACTION]cflags.flag("CHULL::SKIP_STRUCTURE_COMPARISON",true);       //structure comparison is not needed for 1 hull point distance calculation
    //[killed off by SKIP_THERMO_PROPERTIES_EXTRACTION]cflags.flag("CHULL::SKIP_STABILITY_CRITERION_ANALYSIS",true); //this would be circular
    //[killed off by SKIP_THERMO_PROPERTIES_EXTRACTION]cflags.flag("CHULL::SKIP_N+1_ENTHALPY_GAIN_ANALYSIS",true); //this would be circular
    cflags.flag("FAKE_HULL",true); //need to avoid "Very skewed ground-state..." and "Unreliable hull" issues, we're removing points, so these may (and likely will) come up
    cflags.flag("FORCE",true); //need to avoid outlier issues, we're removing points, so these may (and likely will) come up
    //let's skip all this extra output
    ofstream devnull("/dev/null");  //NULL
    //https://stackoverflow.com/questions/366955/obtain-a-stdostream-either-from-stdcout-or-stdofstreamfile`
    bool see_sub_output=false;//true;
    std::streambuf* buf_fh;
    if(see_sub_output){buf_fh=std::cout.rdbuf();} //to cout
    else{buf_fh=devnull.rdbuf();} //to devnull
    //ostream& oss_empty=cout;if(!see_sub_output){oss_empty.setstate(std::ios_base::badbit);}  //like NULL
    //https://stackoverflow.com/questions/7818371/printing-to-nowhere-with-ostream
    //stringstream oss_empty; //also an option, might be better than setting a badbit to cout
    //ostream oss_empty(0); //setting badbit to entirely new instance of object
    ostream oss_fh(buf_fh);
    fake_hull.initialize(cflags,new_points,velements,devnull,oss_fh,m_half_hull,m_add_artificial_unaries);
    //oss_empty.clear();  //clear badbit, as cout is GLOBAL
    if(!fake_hull.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Could not create pseudo convex hull");}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " New (pseudo) hull created" << endl;}
  }

  double ConvexHull::getNPlus1EnthalpyGain(const string& cauid) const {
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Hull not initialized");}
    vector<string> vcauid; vcauid.push_back(cauid);
    vector<double> vnp1eg=getNPlus1EnthalpyGain(vcauid);
    return vnp1eg[0];
  }
  double ConvexHull::getNPlus1EnthalpyGain(const string& cauid,ConvexHull& fake_hull,bool hull_set) const {
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Hull not initialized");}
    vector<string> vcauid; vcauid.push_back(cauid);
    vector<double> vnp1eg=getNPlus1EnthalpyGain(vcauid,fake_hull,hull_set);
    return vnp1eg[0];
  }
  vector<double> ConvexHull::getNPlus1EnthalpyGain(const vector<string>& vcauid) const {ConvexHull fake_hull;return getNPlus1EnthalpyGain(vcauid,fake_hull,false);}
  vector<double> ConvexHull::getNPlus1EnthalpyGain(const vector<string>& vcauid,ConvexHull& fake_hull,bool hull_set) const {
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Hull not initialized");}
    stringstream message;

    uint i_point=AUROSTD_MAX_UINT;
    vector<uint> vcpoints; //follows vcauid
    for(uint i=0,fl_size_i=vcauid.size();i<fl_size_i;i++){
      const string& auid=vcauid[i];
      if(auid.empty()){message << "Empty auid found";throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);}
      if(!findPoint(auid,i_point)){message << "Specified auid not found on hull (auid=" << auid << ")";throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);}
      vcpoints.push_back(i_point);
    }
    return getNPlus1EnthalpyGain(vcpoints,fake_hull,hull_set);
  }
  vector<double> ConvexHull::getNPlus1EnthalpyGain(const vector<uint>& vcpoints) const {ConvexHull fake_hull;return getNPlus1EnthalpyGain(vcpoints,fake_hull,false);}
  vector<double> ConvexHull::getNPlus1EnthalpyGain(const vector<uint>& vcpoints,ConvexHull& fake_hull,bool hull_set) const {
    vector<double> vnp1eg;
    for(uint i=0,fl_size_i=vcpoints.size();i<fl_size_i;i++){vnp1eg.push_back(getNPlus1EnthalpyGain(vcpoints[i],fake_hull,hull_set));}
    return vnp1eg;
  }
  double ConvexHull::getNPlus1EnthalpyGain(uint i_point) const {ConvexHull fake_hull;return getNPlus1EnthalpyGain(i_point,fake_hull,false);}
  double ConvexHull::getNPlus1EnthalpyGain(uint i_point,ConvexHull& fake_hull,bool hull_set) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    stringstream message;

    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Hull not initialized");}
    if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    ChullPoint point=m_points[i_point];
    if(point.m_i_nary==0){return AUROSTD_MAX_DOUBLE;} //N+1(unary)=formation_enthalpy (distance from 0 point) - TECHNICALLY TRUE (allowed in polytope theory), but our definition for N+1 enthalpy gain really starts for binaries (Delta H[{N|1,...,N-1}]), can NOT have 1|0
    bool force_calc_binaries=false; //set TRUE for debugging
    if(point.m_i_nary==0 || (point.m_i_nary==1 && force_calc_binaries==false) ){ //N+1(unary)=formation_enthalpy (distance from 0 point), N+1(binary)=formation_enthalpy (distance from 0 tieline)
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " found " << (point.m_i_nary==0 ? "unary" : "binary") << ", returning abs(H_f_atom)" << endl;} //enthalpy_formation_atom
      if(point.m_has_entry&&H_f_atom(point.m_entry)!=AUROSTD_NAN){return abs(H_f_atom(point.m_entry));} //point.m_entry.enthalpy_formation_atom
      return AUROSTD_MAX_DOUBLE; //0; //return null
    }

    //set fake_hull
    xvector<int> elements_present_points=point.m_elements_present;
    xvector<int> elements_present_hull=elements_present_points;
    if(sum(elements_present_hull)==1){
      //pick ANY other index to be 1 as well
      if(elements_present_hull[elements_present_hull.lrows]==0){elements_present_hull[elements_present_hull.lrows]=1;}
      else{
        if(elements_present_hull.rows<2){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"elements_present.rows<2",_INPUT_ILLEGAL_);} //we cannot have unary hulls as input (facets issue)
        if(elements_present_hull[elements_present_hull.lrows+1]==1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"elements_present mismatch for unary hull",_INDEX_MISMATCH_);}
        elements_present_hull[elements_present_hull.lrows+1]=1;
      }
    }
    if(sum(elements_present_hull)==1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Could not fix sum(elements_present_hull)==1",_RUNTIME_ERROR_);}

    //get full set of points to neglect
    //neglect ALL points with equal or higher dimensionality
    if(hull_set==false){
      uint j_nary=AUROSTD_MAX_UINT,j_alloy=AUROSTD_MAX_UINT;
      if(!getAlloyIndex(point,j_nary,j_alloy)){message << "Alloy index not set (auid=" << point.m_entry.auid << ")";throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);}
      uint i_coord_group=AUROSTD_MAX_UINT;
      vector<uint> vcpoints;
      for(uint fl_size_i_nary=m_naries.size(),i_nary=(fl_size_i_nary-1);i_nary<fl_size_i_nary&&i_nary>=j_nary;i_nary--){ //go backwards!
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " removing i_nary=" << i_nary << endl;}
        for(uint i_alloy=0,fl_size_i_alloy=m_naries[i_nary].m_alloys.size();i_alloy<fl_size_i_alloy;i_alloy++){
          for(uint i=0,fl_size_i=m_naries[i_nary].m_alloys[i_alloy].m_coord_groups.size();i<fl_size_i;i++){
            i_coord_group=m_naries[i_nary].m_alloys[i_alloy].m_coord_groups[i];
            for(uint j=0,fl_size_j=m_coord_groups[i_coord_group].m_points.size();j<fl_size_j;j++){
              vcpoints.push_back(m_coord_groups[i_coord_group].m_points[j]);
            }
          }
        }
      }

      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " points_to_neglect.size()=" << vcpoints.size() << endl;
        std::sort(vcpoints.begin(),vcpoints.end());
        cerr << __AFLOW_FUNC__ << " points_to_neglect: ";
        for(uint i=0,fl_size_i=vcpoints.size();i<fl_size_i;i++){cerr << vcpoints[i] << " ";} 
        cerr << endl;
      }

      getFakeHull(vcpoints,elements_present_hull,elements_present_points,fake_hull);
    }

    //since j_nary and j_alloy don't change, the getDistanceToHull function should work fine (getRelevantFacets())
    //xvector<int> elements_present=m_naries[j_nary].m_alloys[j_alloy].m_elements_present;
    //if(LDEBUG){cerr << __AFLOW_FUNC__ << " elements_present=" << elements_present << endl;}
    //point.setHullCoords(elements_present);  //just to be sure
    point.cleanPointForHullTransfer();      //clean now
    point.reduceCoords(elements_present_hull);  //reduce to minimum hull necessary to calculate scriterion (for unaries this will be a binary)
    point.setHullCoords(); //set to most general coords (m_coords), this reflects relevantFacets()
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " m.coords=" << point.m_coords << endl;}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " h.coords=" << point.h_coords << endl;}

    return fake_hull.getDistanceToHull(point,false);
  }

  void ConvexHull::setStabilityCriterion() {
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Hull not initialized");}
    stringstream message;

    message << "Determining stability criteria for ground-state structures, please be patient";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    vector<uint> g_states=getGStates(true);
    uint i_point=AUROSTD_MAX_UINT,i_coord_group=AUROSTD_MAX_UINT;
    double scriterion;
    vector<uint> eq_gstates;
    for(uint i=0,fl_size_i=g_states.size();i<fl_size_i;i++){
      i_point=g_states[i];
      if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
      if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking at point[i_point=" << i_point << "].m_coords=" << m_points[i_point].m_coords << endl;}
      scriterion=getStabilityCriterion(i_point);
      //CONSISTENCY CHECKS?
      //[OBSOLETE - this is a DISTANCE so it is always positive]if(m_half_hull){
      //[OBSOLETE - this is a DISTANCE so it is always positive]i  //do not use signbit, add tol to zero
      //[OBSOLETE - this is a DISTANCE so it is always positive]i  //sign of distance:
      //[OBSOLETE - this is a DISTANCE so it is always positive]i  //independent of lower/upper hull:  above hull is negative, below hull is positive
      //[OBSOLETE - this is a DISTANCE so it is always positive]i  if(m_lower_hull){
      //[OBSOLETE - this is a DISTANCE so it is always positive]    if(aurostd::notPositive(scriterion,true,ZERO_TOL)) //std::signbit(scriterion))
      //[OBSOLETE - this is a DISTANCE so it is always positive]    { //CO20200106 - patching for auto-indenting
      //[OBSOLETE - this is a DISTANCE so it is always positive]      message << "(lower half hull) found ground-state structure INSIDE hull, suggesting it was not really a ground-state";
      //[OBSOLETE - this is a DISTANCE so it is always positive]      message << " (i_point=" << i_point;
      //[OBSOLETE - this is a DISTANCE so it is always positive]      message << (m_points[i_point].m_entry.auid.empty()?"":",auid="+m_points[i_point].m_entry.auid);
      //[OBSOLETE - this is a DISTANCE so it is always positive]      message << ",scriterion=" << scriterion << ")";
      //[OBSOLETE - this is a DISTANCE so it is always positive]      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
      //[OBSOLETE - this is a DISTANCE so it is always positive]    }
      //[OBSOLETE - this is a DISTANCE so it is always positive]  }else{
      //[OBSOLETE - this is a DISTANCE so it is always positive]    if(aurostd::notNegative(scriterion,true,ZERO_TOL)) //!std::signbit(scriterion))
      //[OBSOLETE - this is a DISTANCE so it is always positive]    { //CO20200106 - patching for auto-indenting
      //[OBSOLETE - this is a DISTANCE so it is always positive]      message << "(upper half hull) found ground-state structure INSIDE hull, suggesting it was not really a ground-state";
      //[OBSOLETE - this is a DISTANCE so it is always positive]      message << " (i_point=" << i_point;
      //[OBSOLETE - this is a DISTANCE so it is always positive]      message << (m_points[i_point].m_entry.auid.empty()?"":",auid="+m_points[i_point].m_entry.auid);
      //[OBSOLETE - this is a DISTANCE so it is always positive]      message << ",scriterion=" << scriterion << ")";
      //[OBSOLETE - this is a DISTANCE so it is always positive]i      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
      //[OBSOLETE - this is a DISTANCE so it is always positive]i    }
      //[OBSOLETE - this is a DISTANCE so it is always positive]i  }
      //[OBSOLETE - this is a DISTANCE so it is always positive]i}
      //KEEP THE SIGN CONVENTION AS IS, it is correct without manipulation
      //set sign convention, negative for outside hull
      //scriterion=abs(scriterion);
      message << "Ground-state ";
      message << (!m_points[i_point].m_entry.compound.empty()?m_points[i_point].m_entry.compound+" ":"");
      message << "(i_point=" << i_point;
      message << (m_points[i_point].m_entry.auid.empty()?"":",auid="+m_points[i_point].m_entry.auid);
      message << ") shows stability criterion = ";
      message << chull::convertUnits(scriterion,(m_formation_energy_hull?_m_:_std_)) << " ";
      message << (m_formation_energy_hull?"meV/atom":"K");
      //this check should STILL work for unaries, as it's the difference of distances, and
      //the actual hull points should ALWAYS be lower than pseudo hull points, hence a positive scriterion
      //by convention, stability criterion are NEGATIVE as they are outside the hull
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      m_points[i_point].m_stability_criterion=scriterion;
      //set equivalent ones
      if(!getCoordGroupIndex(i_point,i_coord_group)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup index not set");}
      if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}
      if(!m_coord_groups[i_coord_group].m_is_on_hull){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup is not on the hull");}
      m_coord_groups[i_coord_group].m_stability_criterion=m_points[i_point].m_stability_criterion;
      //const vector<uint>& eq_gstates=m_coord_groups[i_coord_group].m_equivalent_g_states;
      if(m_coord_groups[i_coord_group].m_calculated_equivalent_g_states){eq_gstates=m_coord_groups[i_coord_group].m_equivalent_g_states;}
      else{eq_gstates=getEquivalentEntries(i_point);}
      if(!eq_gstates.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No equivalent states found (not even self)");}
      for(uint j=0,fl_size_j=eq_gstates.size();j<fl_size_j;j++){
        if(eq_gstates[j]>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
        if(!m_points[eq_gstates[j]].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(eq_gstates[j])+"] is not initialized");}
        scriterion=m_points[i_point].m_stability_criterion;
        scriterion+=(m_points[i_point].getLastCoord()-m_points[eq_gstates[j]].getLastCoord());
        //if(m_formation_energy_hull){scriterion+=abs(H_f_atom(m_points[i_point])-H_f_atom(m_points[eq_gstates[j]]));}
        //else{scriterion+=abs(T_S(m_points[i_point])-T_S(m_points[eq_gstates[j]]));}
        //unaries could be really close to 0
        //don't set if it becomes positive
        if(std::signbit(m_points[i_point].m_stability_criterion)!=std::signbit(scriterion)){scriterion=0;}  //define these to be 0, AlFe hull look at Fe
        m_points[eq_gstates[j]].m_stability_criterion=scriterion;
      }
    }
  }

  void ConvexHull::setNPlus1EnthalpyGain(uint i_point) {ConvexHull fake_hull;return setNPlus1EnthalpyGain(i_point,fake_hull,false);}
  void ConvexHull::setNPlus1EnthalpyGain(uint i_point,ConvexHull& fake_hull,bool hull_set) {
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Hull not initialized");}
    stringstream message;

    if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
    //[CO20200404 - not useful, we still want values for binaries, this would skip over value assignment]if(test_binary_np1==false&&m_points[i_point].m_i_nary<2){continue;} //unaries would be trivial, and binaries is simply H_f
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking at point[i_point=" << i_point << "].m_coords=" << m_points[i_point].m_coords << ", compound=" << m_points[i_point].m_entry.compound << endl;}
    double np1egain=getNPlus1EnthalpyGain(i_point,fake_hull,hull_set);
    //CONSISTENCY CHECKS?
    //yes - for unaries/binaries N+1 == formation enthalpy
    //[MATHEMATICALLY ALLOWED but our definition of N+1 energy does not allow for unaries, so don't check]if(m_points[i_point].m_i_nary==0 || m_points[i_point].m_i_nary==1)
    if(m_points[i_point].m_i_nary==1) {
      if(m_points[i_point].m_has_entry&&H_f_atom(m_points[i_point].m_entry)!=AUROSTD_NAN){  //m_points[i_point].m_entry.enthalpy_formation_atom
        if(!aurostd::identical(abs(np1egain),abs(H_f_atom(m_points[i_point].m_entry)),ZERO_MEV_TOL)){ //m_points[i_point].m_entry.enthalpy_formation_atom
          message << "abs(np1egain) != abs(H_f_atom(m_points[i_point].m_entry)) [ " << abs(np1egain) << " != " << abs(H_f_atom(m_points[i_point].m_entry)) << " ], please check"; //m_points[i_point].m_entry.enthalpy_formation_atom
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_); //throw an error because we are not in debug mode (see force_calc_binaries)
          //pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        }
      }
    }
    //[OBSOLETE - this is a DISTANCE so it is always positive]if(m_points[i_point].m_i_nary>0){ //skip unaries for this check
    //[OBSOLETE - this is a DISTANCE so it is always positive]  if(m_half_hull){
    //[OBSOLETE - this is a DISTANCE so it is always positive]    //do not use signbit, add tol to zero
    //[OBSOLETE - this is a DISTANCE so it is always positive]    //sign of distance:
    //[OBSOLETE - this is a DISTANCE so it is always positive]    //independent of lower/upper hull:  above hull is negative, below hull is positive
    //[OBSOLETE - this is a DISTANCE so it is always positive]    if(m_lower_hull){
    //[OBSOLETE - this is a DISTANCE so it is always positive]      if(aurostd::notPositive(np1egain,true,ZERO_TOL)) //std::signbit(np1egain))
    //[OBSOLETE - this is a DISTANCE so it is always positive]      { //CO20200106 - patching for auto-indenting
    //[OBSOLETE - this is a DISTANCE so it is always positive]        message << "(lower half hull) found ground-state structure INSIDE hull, suggesting it was not really a ground-state";
    //[OBSOLETE - this is a DISTANCE so it is always positive]        message << " (i_point=" << i_point;
    //[OBSOLETE - this is a DISTANCE so it is always positive]        message << (m_points[i_point].m_entry.auid.empty()?"":",auid="+m_points[i_point].m_entry.auid);
    //[OBSOLETE - this is a DISTANCE so it is always positive]        message << ",np1egain=" << np1egain << ")";
    //[OBSOLETE - this is a DISTANCE so it is always positive]        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
    //[OBSOLETE - this is a DISTANCE so it is always positive]      }
    //[OBSOLETE - this is a DISTANCE so it is always positive]    }else{
    //[OBSOLETE - this is a DISTANCE so it is always positive]      if(aurostd::notNegative(np1egain,true,ZERO_TOL)) //!std::signbit(np1egain))
    //[OBSOLETE - this is a DISTANCE so it is always positive]      { //CO20200106 - patching for auto-indenting
    //[OBSOLETE - this is a DISTANCE so it is always positive]        message << "(upper half hull) found ground-state structure INSIDE hull, suggesting it was not really a ground-state";
    //[OBSOLETE - this is a DISTANCE so it is always positive]        message << " (i_point=" << i_point;
    //[OBSOLETE - this is a DISTANCE so it is always positive]        message << (m_points[i_point].m_entry.auid.empty()?"":",auid="+m_points[i_point].m_entry.auid);
    //[OBSOLETE - this is a DISTANCE so it is always positive]        message << ",np1egain=" << np1egain << ")";
    //[OBSOLETE - this is a DISTANCE so it is always positive]        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
    //[OBSOLETE - this is a DISTANCE so it is always positive]      }
    //[OBSOLETE - this is a DISTANCE so it is always positive]    }
    //[OBSOLETE - this is a DISTANCE so it is always positive]  }
    //[OBSOLETE - this is a DISTANCE so it is always positive]i}
    //KEEP THE SIGN CONVENTION AS IS, it is correct without manipulation
    //set sign convention, negative for outside hull
    //np1egain=abs(np1egain);
    if(m_points[i_point].m_i_nary>0){ //don't print for unaries
      message << "Ground-state ";
      message << (!m_points[i_point].m_entry.compound.empty()?m_points[i_point].m_entry.compound:"")+" ";
      message << "(i_point=" << i_point;
      message << (m_points[i_point].m_entry.auid.empty()?"":",auid="+m_points[i_point].m_entry.auid);
      message << ") shows N+1 enthalpy gain = ";
      message << chull::convertUnits(np1egain,(m_formation_energy_hull?_m_:_std_)) << " ";
      message << (m_formation_energy_hull?"meV/atom":"K");
    }
    //this check should STILL work for unaries, as it's the difference of distances, and
    //the actual hull points should ALWAYS be lower than pseudo hull points, hence a positive np1egain
    //by convention, enthalpy gains are NEGATIVE as they are outside the hull
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    m_points[i_point].m_n_plus_1_enthalpy_gain=np1egain;
    //set equivalent ones
    uint i_coord_group=AUROSTD_MAX_UINT;
    if(!getCoordGroupIndex(i_point,i_coord_group)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup index not set");}
    if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized coordgroup");}
    if(!m_coord_groups[i_coord_group].m_is_on_hull){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup is not on the hull");}
    m_coord_groups[i_coord_group].m_n_plus_1_enthalpy_gain=m_points[i_point].m_n_plus_1_enthalpy_gain;
    vector<uint> eq_gstates;
    //const vector<uint>& eq_gstates=m_coord_groups[i_coord_group].m_equivalent_g_states;
    if(m_coord_groups[i_coord_group].m_calculated_equivalent_g_states){eq_gstates=m_coord_groups[i_coord_group].m_equivalent_g_states;}
    else{eq_gstates=getEquivalentEntries(i_point);}
    if(!eq_gstates.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No equivalent states found (not even self)");}
    for(uint j=0,fl_size_j=eq_gstates.size();j<fl_size_j;j++){
      if(eq_gstates[j]>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
      if(!m_points[eq_gstates[j]].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(eq_gstates[j])+"] is not initialized");}
      np1egain=m_points[i_point].m_n_plus_1_enthalpy_gain;
      np1egain+=(m_points[i_point].getLastCoord()-m_points[eq_gstates[j]].getLastCoord());
      //if(m_formation_energy_hull){np1egain+=abs(H_f_atom(m_points[i_point])-H_f_atom(m_points[eq_gstates[j]]));}
      //else{np1egain+=abs(T_S(m_points[i_point])-T_S(m_points[eq_gstates[j]]));}
      //unaries could be really close to 0
      //don't set if it becomes positive
      if(std::signbit(m_points[i_point].m_n_plus_1_enthalpy_gain)!=std::signbit(np1egain)){np1egain=0;}  //define these to be 0, AlFe hull look at Fe
      m_points[eq_gstates[j]].m_n_plus_1_enthalpy_gain=np1egain;
    }
  }
  void ConvexHull::setNPlus1EnthalpyGain() {
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Hull not initialized");}
    stringstream message;

    message << "Determining N+1 enthalpy gain for ground-state structures, please be patient";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    vector<uint> g_states=getGStates(true);

    //CO20200804 - group g_states by elements_present, they will all need the same fake_hull
    vector<xvector<int> > v_elements_present;
    vector<vector<uint> > v_g_states;

    uint i=0,j=0,fl_size_i=0,fl_size_j=0;
    uint i_point=AUROSTD_MAX_UINT;
    bool found=false;
    for(i=0,fl_size_i=g_states.size();i<fl_size_i;i++){
      i_point=g_states[i];
      if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
      if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
      //[CO20200404 - not useful, we still want values for binaries, this would skip over value assignment]if(test_binary_np1==false&&m_points[i_point].m_i_nary<2){continue;} //unaries would be trivial, and binaries is simply H_f
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking at point[i_point=" << i_point << "].m_coords=" << m_points[i_point].m_coords << ", compound=" << m_points[i_point].m_entry.compound << endl;}

      found=false;
      for(j=0,fl_size_j=v_elements_present.size();j<fl_size_j&&found==false;j++){
        if(m_points[i_point].m_elements_present==v_elements_present[j]){
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << " adding point[i_point=" << i_point << ",auid=" << m_points[i_point].m_entry.auid << "].compound=\"" << m_points[i_point].m_entry.compound << "\" ";
            cerr << "with elements_present=" << m_points[i_point].m_elements_present << " to v_elements_present[j=" << j << "]=" << v_elements_present[j] << endl;
          }
          v_g_states[j].push_back(g_states[i]);
          found=true;
        }
      }
      if(!found){
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " adding point[i_point=" << i_point << ",auid=" << m_points[i_point].m_entry.auid << "].compound=\"" << m_points[i_point].m_entry.compound << "\" ";
          cerr << "with elements_present=" << m_points[i_point].m_elements_present << " to set" << endl;
        }
        v_elements_present.push_back(m_points[i_point].m_elements_present);
        v_g_states.push_back(vector<uint>(0));
        v_g_states.back().push_back(g_states[i]);
      }
    }
    if(v_elements_present.size()!=v_g_states.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"v_elements_present.size()!=v_g_states.size()",_INDEX_MISMATCH_);}

    ConvexHull fake_hull;
    //[CO20200404 - not useful, we still want values for binaries, this would skip over value assignment]bool test_binary_np1=true;  //turn me off when you're sure np1 is working
    for(i=0,fl_size_i=v_g_states.size();i<fl_size_i;i++){
      /////////////////////////////////////////////////////////////
      if(v_g_states[i].empty()){continue;}
      i_point=v_g_states[i][0]; //g_states[i];
      setNPlus1EnthalpyGain(i_point,fake_hull,false);
      for(j=1,fl_size_j=v_g_states[i].size();j<fl_size_j;j++){
        i_point=v_g_states[i][j]; //g_states[i];
        setNPlus1EnthalpyGain(i_point,fake_hull,true);
      }
    }
  }

  void ConvexHull::cleanHull() {
    h_dim=0;
    m_elements_present.clear();
    h_points.clear();
    h_centroid.clear();
    h_reference.clear();
    h_facets.clear();
    h_visible_facets.clear();
    h_horizon_ridges.clear();
  }

  string ConvexHull::prettyPrintCompound(const ChullPoint& point,vector_reduction_type vred,bool exclude1,filetype ftype) const {  // overload
    if(!point.m_has_entry){
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No entry found");
    }
    return prettyPrintCompound(point.m_entry,vred,exclude1,ftype);
  }

  string ConvexHull::prettyPrintCompound(const aflowlib::_aflowlib_entry& entry,vector_reduction_type vred,bool exclude1,filetype ftype) const {  // overload
    if(entry.vspecies.size()!=entry.vcomposition.size()) {
      if(entry.prototype.find("POCC")!=string::npos){ //POCC entries have no composition
        return pflow::prettyPrintCompound(entry.vspecies,entry.vstoichiometry,no_vrt,exclude1,ftype);  //pass in stoichiometry and do not reduce
      }else{
        stringstream message;
        message << "Entry (auid=" << entry.auid << ") is ill-defined: vspecies.size()!=vcomposition.size()";
        message << " (please report on AFLOW Forum: aflow.org/forum)";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        return entry.compound;
      }
    }
    return pflow::prettyPrintCompound(entry.vspecies,entry.vcomposition,vred,exclude1,ftype);  //ME20190628
  }

  //[ME20190628 - moved to pflow_funcs.cpp] string ConvexHull::prettyPrintCompound(const vector<string>& vspecies,const vector<double>& vcomposition,vector_reduction_type vred,bool exclude1,filetype ftype) const {  // overload
  //[ME20190628 - moved to pflow_funcs.cpp]   return prettyPrintCompound(vspecies,aurostd::vector2xvector<double>(vcomposition),vred,exclude1,ftype);
  //[ME20190628 - moved to pflow_funcs.cpp] }

  //[ME20190628 - moved to pflow_funcs.cpp] string ConvexHull::prettyPrintCompound(const vector<string>& vspecies,const xvector<double>& vcomposition,vector_reduction_type vred,bool exclude1,filetype ftype) const {  // main function
  //[ME20190628 - moved to pflow_funcs.cpp]   // creates compound_label for LaTeX and text docs, like adding $_{}$
  //[ME20190628 - moved to pflow_funcs.cpp]   // 2-D, we usually want vred=gcd_vrt true for convex points, and no_vrt elsewhere
  //[ME20190628 - moved to pflow_funcs.cpp]   uint precision=COEF_PRECISION;
  //[ME20190628 - moved to pflow_funcs.cpp]   stringstream output;output.precision(precision);
  //[ME20190628 - moved to pflow_funcs.cpp]   if(vspecies.size()!=(uint)vcomposition.rows) {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"vspecies.size() != vcomposition.rows");}
  //[ME20190628 - moved to pflow_funcs.cpp]   // special case, unary
  //[ME20190628 - moved to pflow_funcs.cpp]   if(vspecies.size() == 1) {
  //[ME20190628 - moved to pflow_funcs.cpp]     output << vspecies[0];
  //[ME20190628 - moved to pflow_funcs.cpp]     if(!exclude1) {output << (vred==gcd_vrt?1:vcomposition[vcomposition.lrows]);}
  //[ME20190628 - moved to pflow_funcs.cpp]     return output.str();
  //[ME20190628 - moved to pflow_funcs.cpp]   }
  //[ME20190628 - moved to pflow_funcs.cpp]   xvector<double> comp=vcomposition;
  //[ME20190628 - moved to pflow_funcs.cpp]   if(vred==gcd_vrt){comp=aurostd::reduceByGCD(comp,ZERO_TOL);}
  //[ME20190628 - moved to pflow_funcs.cpp]   else if(vred==frac_vrt){comp=aurostd::normalizeSumToOne(comp,ZERO_TOL);}
  //[ME20190628 - moved to pflow_funcs.cpp]   else if(vred==no_vrt){;}
  //[ME20190628 - moved to pflow_funcs.cpp]   else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown reduce mode",_INPUT_UNKNOWN_);}
  //[ME20190628 - moved to pflow_funcs.cpp]   if(aurostd::zeroWithinTol(aurostd::sum(comp),ZERO_TOL)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Empty composition");}
  //[ME20190628 - moved to pflow_funcs.cpp]   for(uint i=0,fl_size_i=vspecies.size();i<fl_size_i;i++) {
  //[ME20190628 - moved to pflow_funcs.cpp]     output << vspecies[i];
  //[ME20190628 - moved to pflow_funcs.cpp]     if(!(exclude1 && aurostd::identical(comp[i+comp.lrows],1.0,ZERO_TOL))) {
  //[ME20190628 - moved to pflow_funcs.cpp]       if(ftype==latex_ft) {output << "$_{";
  //[ME20190628 - moved to pflow_funcs.cpp]       }else if(ftype==gnuplot_ft){output<< "_{";}
  //[ME20190628 - moved to pflow_funcs.cpp]       output << comp[i+comp.lrows];
  //[ME20190628 - moved to pflow_funcs.cpp]       if(ftype==latex_ft) {output << "}$";}
  //[ME20190628 - moved to pflow_funcs.cpp]       else if(ftype==gnuplot_ft){output<< "}";}
  //[ME20190628 - moved to pflow_funcs.cpp]     }
  //[ME20190628 - moved to pflow_funcs.cpp]   }
  //[ME20190628 - moved to pflow_funcs.cpp]   return output.str();
  //[ME20190628 - moved to pflow_funcs.cpp] }

  string ConvexHull::getICSDNumber(uint i_point,bool remove_suffix) const {
    if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    return getICSDNumber(m_points[i_point],remove_suffix);
  }

  string ConvexHull::getICSDNumber(const ChullPoint& point,bool remove_suffix) const {
    if(!point.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized point");}
    if(!point.m_has_entry){return "";}
    return getICSDNumber(point.m_entry,remove_suffix);

  }
  string ConvexHull::getICSDNumber(const aflowlib::_aflowlib_entry& entry,bool remove_suffix) const {
    const string& proto=entry.prototype;
    if(!aurostd::substring2bool(proto,"ICSD_")){return "";}
    vector<string> tokens;
    aurostd::string2tokens(proto,tokens,"_");
    string icsd_number;
    if(tokens.size()==2){icsd_number=tokens[1];} //aurostd::substring2bool(proto,"ICSD_")
    else if(tokens.size()==3&&aurostd::substring2bool(proto,"_ICSD_")){icsd_number=tokens[2];}
    if(icsd_number.empty()){return icsd_number;}
    if(remove_suffix){
      aurostd::string2tokens(icsd_number,tokens,".");if(tokens.size()){icsd_number=tokens[0];}
      aurostd::string2tokens(icsd_number,tokens,":");if(tokens.size()){icsd_number=tokens[0];}
    }
    return icsd_number;
  }

  string ConvexHull::prettyPrintPrototype(const ChullPoint& point, bool double_back_slash,bool icsd_label_skim) const {  // overload
    if(!point.m_has_entry){
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No entry found");
    }
    return prettyPrintPrototype(point.m_entry,double_back_slash,icsd_label_skim);
  }

  string ConvexHull::prettyPrintPrototype(const aflowlib::_aflowlib_entry& entry, bool double_back_slash,bool icsd_label_skim) const {  // main function
    // creates prototype_label for LaTeX ONLY, no use for this function otherwise
    // escapes funny characters
    if(entry.prototype.empty()) {
      stringstream message;
      message << "Entry (auid=" << entry.auid << ") is ill-defined: empty prototype";
      message << " (please report on AFLOW Forum: aflow.org/forum)";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
      return entry.prototype;
    }
    bool fix_icsd_labels=true;
    string proto=entry.prototype;
    if(fix_icsd_labels&&aurostd::substring2bool(proto,"ICSD_")){
      vector<string> tokens;
      string string_to_fix=getICSDNumber(entry,icsd_label_skim==true);  //remove suffix if skim label
      //aurostd::string2tokens(proto,tokens,"_");
      //if(tokens.size()==2){string_to_fix=tokens[1];} //aurostd::substring2bool(proto,"ICSD_")
      //else if(tokens.size()==3&&aurostd::substring2bool(proto,"_ICSD_")){string_to_fix=tokens[2];}
      if(!string_to_fix.empty()){
        return (icsd_label_skim?"":"ICSD~")+(double_back_slash?string("\\\\"):string("\\"))+"#"+aurostd::fixStringLatex(string_to_fix,double_back_slash,false);
      }
      //else, just leave alone
    }
    return aurostd::fixStringLatex(proto,double_back_slash,false);
  }

  //[CO20190419 - moved to aurostd_main.cpp]string ConvexHull::fixStringLatex(const string& input, bool double_back_slash,bool symmetry_string) const {
  //[CO20190419 - moved to aurostd_main.cpp]  // deals with special characters for LaTeX, like some characters in prototype
  //[CO20190419 - moved to aurostd_main.cpp]  // see http://tex.stackexchange.com/questions/34580/escape-character-in-latex
  //[CO20190419 - moved to aurostd_main.cpp]  // double_back_slash was needed SOMETIMES for gnuplot output, as one backslash
  //[CO20190419 - moved to aurostd_main.cpp]  // went away when writing to file, and  -- OBSOLETE NOW
  //[CO20190419 - moved to aurostd_main.cpp]  string output;
  //[CO20190419 - moved to aurostd_main.cpp]  vector<char> problem_characters;
  //[CO20190419 - moved to aurostd_main.cpp]  problem_characters.push_back('&');
  //[CO20190419 - moved to aurostd_main.cpp]  problem_characters.push_back('%');
  //[CO20190419 - moved to aurostd_main.cpp]  problem_characters.push_back('$');
  //[CO20190419 - moved to aurostd_main.cpp]  problem_characters.push_back('#');
  //[CO20190419 - moved to aurostd_main.cpp]  if(!symmetry_string) {
  //[CO20190419 - moved to aurostd_main.cpp]    problem_characters.push_back('_');
  //[CO20190419 - moved to aurostd_main.cpp]    problem_characters.push_back('{');
  //[CO20190419 - moved to aurostd_main.cpp]    problem_characters.push_back('}');
  //[CO20190419 - moved to aurostd_main.cpp]  }
  //[CO20190419 - moved to aurostd_main.cpp]  problem_characters.push_back('~');  // different fix
  //[CO20190419 - moved to aurostd_main.cpp]  problem_characters.push_back('^');  // different fix
  //[CO20190419 - moved to aurostd_main.cpp]  string solution_string;
  //[CO20190419 - moved to aurostd_main.cpp]  solution_string = "\\\\";  // has to be string, \\ char does not work
  //[CO20190419 - moved to aurostd_main.cpp]  bool found_escaped_char;
  //[CO20190419 - moved to aurostd_main.cpp]  bool found_hyphen_symmetry = false;
  //[CO20190419 - moved to aurostd_main.cpp]  bool solved_hyphen_symmetry = false;
  //[CO20190419 - moved to aurostd_main.cpp]  for(uint i=0;i<input.length();i++) {
  //[CO20190419 - moved to aurostd_main.cpp]    // we first enter this loop because symmetry_string and input[i]=='-'
  //[CO20190419 - moved to aurostd_main.cpp]    // second enter loop because symmetry_string and found_hyphen_symmetry
  //[CO20190419 - moved to aurostd_main.cpp]    if(symmetry_string && (input[i] == '-' || found_hyphen_symmetry)) {
  //[CO20190419 - moved to aurostd_main.cpp]      if(!found_hyphen_symmetry) {
  //[CO20190419 - moved to aurostd_main.cpp]        // first enter loop, come here
  //[CO20190419 - moved to aurostd_main.cpp]        found_hyphen_symmetry = true;
  //[CO20190419 - moved to aurostd_main.cpp]        output.append((double_back_slash?string("\\"):string(""))+string("\\overline{"));
  //[CO20190419 - moved to aurostd_main.cpp]        // very important, we don't want to add hyphen, just replace
  //[CO20190419 - moved to aurostd_main.cpp]        // with overline, so continue
  //[CO20190419 - moved to aurostd_main.cpp]        continue;
  //[CO20190419 - moved to aurostd_main.cpp]      }else{
  //[CO20190419 - moved to aurostd_main.cpp]        // second enter loop, do nothing but turn this flag on
  //[CO20190419 - moved to aurostd_main.cpp]        // allow us to add input[i]
  //[CO20190419 - moved to aurostd_main.cpp]        found_hyphen_symmetry = false;
  //[CO20190419 - moved to aurostd_main.cpp]        solved_hyphen_symmetry = true;
  //[CO20190419 - moved to aurostd_main.cpp]      }
  //[CO20190419 - moved to aurostd_main.cpp]    }else{
  //[CO20190419 - moved to aurostd_main.cpp]      if(symmetry_string && solved_hyphen_symmetry) {
  //[CO20190419 - moved to aurostd_main.cpp]        // last step of symmetry_string fix, but we have to do this in part of
  //[CO20190419 - moved to aurostd_main.cpp]        // the loop to allow for next character to be identified as problem
  //[CO20190419 - moved to aurostd_main.cpp]        // character as well
  //[CO20190419 - moved to aurostd_main.cpp]        output.append(1, '}');
  //[CO20190419 - moved to aurostd_main.cpp]        solved_hyphen_symmetry = false;
  //[CO20190419 - moved to aurostd_main.cpp]      }
  //[CO20190419 - moved to aurostd_main.cpp]      // go through all problem characters
  //[CO20190419 - moved to aurostd_main.cpp]      for(uint j=0,fl_size_j=problem_characters.size();j<fl_size_j;j++) {
  //[CO20190419 - moved to aurostd_main.cpp]        if(input[i] == problem_characters[j]) {
  //[CO20190419 - moved to aurostd_main.cpp]          if(double_back_slash) {
  //[CO20190419 - moved to aurostd_main.cpp]            // if we find one, but it has double backslash, leave alone
  //[CO20190419 - moved to aurostd_main.cpp]            // doesn't matter what it is, if it has double backslash it's good
  //[CO20190419 - moved to aurostd_main.cpp]            // if we find one, but it only has single backslash, add one
  //[CO20190419 - moved to aurostd_main.cpp]            if(i && i - 1 && input[i - 1] == '\\' && input[i - 2] == '\\') {break;}
  //[CO20190419 - moved to aurostd_main.cpp]            else if(i && input[i - 1] == '\\') {
  //[CO20190419 - moved to aurostd_main.cpp]              output.append(1, '\\');  // just add one
  //[CO20190419 - moved to aurostd_main.cpp]              break;
  //[CO20190419 - moved to aurostd_main.cpp]            }
  //[CO20190419 - moved to aurostd_main.cpp]            // if we find one, give two backslashes
  //[CO20190419 - moved to aurostd_main.cpp]            output.append("\\\\");
  //[CO20190419 - moved to aurostd_main.cpp]            break;
  //[CO20190419 - moved to aurostd_main.cpp]          }else{
  //[CO20190419 - moved to aurostd_main.cpp]            // if we find one, but it has single backslash, leave alone
  //[CO20190419 - moved to aurostd_main.cpp]            // doesn't matter what it is, if it has single backslash it's good
  //[CO20190419 - moved to aurostd_main.cpp]            // if we find one, give single backslash
  //[CO20190419 - moved to aurostd_main.cpp]            if(i && input[i - 1] == '\\') {break;}  
  //[CO20190419 - moved to aurostd_main.cpp]            output.append(1, '\\');
  //[CO20190419 - moved to aurostd_main.cpp]            break;
  //[CO20190419 - moved to aurostd_main.cpp]          }
  //[CO20190419 - moved to aurostd_main.cpp]        }
  //[CO20190419 - moved to aurostd_main.cpp]      }
  //[CO20190419 - moved to aurostd_main.cpp]      // we also have to add {} for these characters
  //[CO20190419 - moved to aurostd_main.cpp]      if(input[i] == '~' || input[i] == '^') {output.append("{}");}
  //[CO20190419 - moved to aurostd_main.cpp]      found_escaped_char = false;
  //[CO20190419 - moved to aurostd_main.cpp]      if(input[i] == '\\') {
  //[CO20190419 - moved to aurostd_main.cpp]        for(uint j=0,fl_size_j=problem_characters.size();j<fl_size_j;j++) {
  //[CO20190419 - moved to aurostd_main.cpp]          // the only way this works if it's serving as an escape for a character
  //[CO20190419 - moved to aurostd_main.cpp]          // don't worry about double backslash here, we get to that when we find
  //[CO20190419 - moved to aurostd_main.cpp]          // the actual character
  //[CO20190419 - moved to aurostd_main.cpp]          if(i != (input.length() - 1) && input[i+1] == problem_characters[j]) {
  //[CO20190419 - moved to aurostd_main.cpp]            found_escaped_char = true;
  //[CO20190419 - moved to aurostd_main.cpp]            break;  // doesn't matter what it is, if it has backslash it's good
  //[CO20190419 - moved to aurostd_main.cpp]          }
  //[CO20190419 - moved to aurostd_main.cpp]        }
  //[CO20190419 - moved to aurostd_main.cpp]        // this is a problem, no way around it--we cannot output single backslash
  //[CO20190419 - moved to aurostd_main.cpp]        if(!found_escaped_char) {
  //[CO20190419 - moved to aurostd_main.cpp]          stringstream message;
  //[CO20190419 - moved to aurostd_main.cpp]          message << "Extraneous backslash found in \"" << input << "\" which may cause problems for LaTeX/gnuplot";
  //[CO20190419 - moved to aurostd_main.cpp]          pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_WARNING_);
  //[CO20190419 - moved to aurostd_main.cpp]          return input;
  //[CO20190419 - moved to aurostd_main.cpp]        }
  //[CO20190419 - moved to aurostd_main.cpp]      }
  //[CO20190419 - moved to aurostd_main.cpp]    }
  //[CO20190419 - moved to aurostd_main.cpp]    // add in character from input
  //[CO20190419 - moved to aurostd_main.cpp]    output.append(1, input[i]);
  //[CO20190419 - moved to aurostd_main.cpp]  }
  //[CO20190419 - moved to aurostd_main.cpp]  return output;
  //[CO20190419 - moved to aurostd_main.cpp]}

  string ConvexHull::getPlotHeaderPDF(char function_mode,const string& column_header,bool display_color_gradient) const {
    // produces addplot latex string
    stringstream message;
    stringstream addplot_output_ss;addplot_output_ss.str("");

    uint dimension=getDim();
    bool reverse_axes=DEFAULT_CHULL_LATEX_REVERSE_AXIS;

    // http://tex.stackexchange.com/questions/59070/pgfplots-remove-darker-borders-on-marks
    // this string WILL cause warnings when compiling, but it's correct, no way to
    // fix it with our current compiler
    stringstream tmp_dark_border_points_command;
    tmp_dark_border_points_command << "scatter/use mapped color={draw=black,fill=mapped color,solid}";

    string dark_border_points_command = tmp_dark_border_points_command.str();

    // first do a check that the function_mode is not out of scope
    if(dimension == 2) {
      if(function_mode != ADDPLOT_MODE_HULL_POINTS &&
          function_mode != ADDPLOT_MODE_OFF_HULL_POINTS &&
          function_mode != ADDPLOT_MODE_HULL_FACETS) {
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Input dimension mode (2D) and function mode mismatch");
      }
      if(m_velements.size() != 2) {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Input dimension mode (2D) and elements.size() mismatch");}
    }else if(dimension == 3) {
      if(function_mode != ADDPLOT_MODE_HULL_POINTS &&
          function_mode != ADDPLOT_MODE_OFF_HULL_POINTS &&
          function_mode != ADDPLOT_MODE_HULL_FACETS &&
          function_mode != ADDPLOT_MODE_HULL_FACETSDROP_SHADOWS &&
          function_mode != ADDPLOT_MODE_HEAT_MAPS) {
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Input dimension mode (3D) and function mode mismatch");
      }
      if(m_velements.size() != 3) {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Input dimension mode (3D) and elements.size() mismatch");}
    }else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Input dimension mode out of scope");}

    //also check column_header
    bool found=false;
    for(uint i=0,fl_size_i=m_velements.size();i<fl_size_i;i++){found=(found||(column_header==m_velements[i]));}
    found=(found||(column_header=="H_f_meVatom"));
    found=(found||(column_header=="H_f_kJmol"));
    found=(found||(column_header=="T_S"));
    found=(found||(column_header=="Dist2hull"));
    if(!found){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Input column_header mode out of scope ("+column_header+")");}

    // main addplot function
    if(dimension == 2) {addplot_output_ss << "\\addplot+[" << endl;}
    else{  // dimension==3
      if(function_mode == ADDPLOT_MODE_HEAT_MAPS) {addplot_output_ss << "\\addplot3[" << endl;}
      else{addplot_output_ss << "\\addplot3+[" << endl;}
    }

    if(function_mode == ADDPLOT_MODE_HULL_FACETS ||
        function_mode == ADDPLOT_MODE_HULL_FACETSDROP_SHADOWS) {
      addplot_output_ss << "mark=none," << endl;
      if(function_mode == ADDPLOT_MODE_HULL_FACETSDROP_SHADOWS) {
        addplot_output_ss << "color=white," << endl;
        addplot_output_ss << "solid," << endl;
        addplot_output_ss << "line width=2.0pt," << endl;
      }else{  // function_mode==ADDPLOT_MODE_HULL_FACETS
        addplot_output_ss << "color=black," << endl;
        addplot_output_ss << "solid," << endl;
        if(dimension == 3) {addplot_output_ss << "line width=2pt," << endl;}
      }
      addplot_output_ss << "] table ";
      if(dimension == 2) {
        if(reverse_axes){addplot_output_ss << "[x=" << m_velements[0];}
        else{addplot_output_ss << "[x=" << m_velements[1];}
        addplot_output_ss << ",y=" << column_header;
        addplot_output_ss << "]";
      }else{  // dimension==3
        addplot_output_ss << "[x=" << m_velements[1];
        addplot_output_ss << ",y=" << m_velements[0];
        addplot_output_ss << ",z=" << m_velements[2] << "]";  // flipping x and y
      }
    }else{  // function_mode==ADDPLOT_MODE_HULL_POINTS||function_mode==ADDPLOT_MODE_OFF_HULL_POINTS||function_mode==ADDPLOT_MODE_HEAT_MAPS
      if(function_mode == ADDPLOT_MODE_HEAT_MAPS) {
        addplot_output_ss << "patch," << endl;
        addplot_output_ss << "patch type=triangle," << endl;
        addplot_output_ss << "shader=interp," << endl;
      }else{
        addplot_output_ss << "only marks," << endl;  // IMPORTANT, no lines
        if(function_mode == ADDPLOT_MODE_HULL_POINTS) {
          addplot_output_ss << "mark=*," << endl;
          if(dimension == 2) {addplot_output_ss << "mark size=4," << endl;}
          else{  // dimension==3
            addplot_output_ss << "mark size=5," << endl;
            addplot_output_ss << "line width=2pt," << endl;
          }
        }else{  // function_mode==ADDPLOT_MODE_OFF_HULL_POINTS
          addplot_output_ss << "mark=x," << endl;
          addplot_output_ss << "mark options={scale=2,line width=2,solid";
        }
      }
      if(display_color_gradient || function_mode == ADDPLOT_MODE_HEAT_MAPS) {  
        // this OR statement doesn't practically matter, since one does not work without the other, 
        // but I keep it here so that the output of the string is complete
        if(function_mode == ADDPLOT_MODE_OFF_HULL_POINTS) {addplot_output_ss << "}," << endl;}
        addplot_output_ss << "point meta=\\thisrow{" << column_header << "}," << endl;  // uses point meta as color data
        addplot_output_ss << "nodes near coords*={}," << endl;  // no labels, but we need this for colors
        // GOT IT! // http://tex.stackexchange.com/questions/59070/pgfplots-remove-darker-borders-on-marks
        if(function_mode == ADDPLOT_MODE_HULL_POINTS) {addplot_output_ss << dark_border_points_command << "," << endl;}
        addplot_output_ss << "visualization depends on={\\thisrow{" << column_header << "} \\as \\" << column_header << "}," << endl;  // defines visualization dependency
      }else{
        if(function_mode == ADDPLOT_MODE_HULL_POINTS) {addplot_output_ss << "mark options={draw=black,fill=blue,solid}," << endl;}
        else{addplot_output_ss << ",draw=red}," << endl;}// function_mode==ADDPLOT_MODE_OFF_HULL_POINTS
      }
      addplot_output_ss << "] table ";
      if(dimension == 2) {
        if(reverse_axes){addplot_output_ss << "[x=" << m_velements[0];}
        else{addplot_output_ss << "[x=" << m_velements[1];}
        addplot_output_ss << ",y=" << column_header;
        addplot_output_ss << "]";
      }else{  // dimension==3
        addplot_output_ss << "[x=" << m_velements[1];
        addplot_output_ss << ",y=" << m_velements[0];
        addplot_output_ss << ",z=" << m_velements[2] << "]";  // flipping x and y
      }
    }
    addplot_output_ss << "{" << endl;
    for(uint i=0,fl_size_i=m_velements.size();i<fl_size_i;i++){addplot_output_ss << aurostd::PaddedPOST(m_velements[i], 30);}
    addplot_output_ss << aurostd::PaddedPOST("H_f_meVatom", 30);
    addplot_output_ss << aurostd::PaddedPOST("H_f_kJmol", 30);
    addplot_output_ss << aurostd::PaddedPOST("T_S", 30);
    if(PRINT_DIST2HULL_COL_TEX){
      addplot_output_ss << aurostd::PaddedPOST("Dist2hull", 30);
    }
    addplot_output_ss << endl;

    return addplot_output_ss.str();
  }

  string ConvexHull::getPlotPointContentPDF(const ChullPoint& point,bool zero_end_point,bool zero_dist_2_hull) const {  //true,false
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " starting" << endl;}

    // initializations
    stringstream addplot_output_ss;
    // no precision
    stringstream num_ss;

    // explicit dimensions
    for(int i=point.s_coords.lrows;i<=point.s_coords.urows;i++) {
      addplot_output_ss << aurostd::PaddedPOST(aurostd::utype2string(point.s_coords[i],FULL_PRECISION), 30);
    }
    // enthalpy of formation, row 4
    // fix for unaries, set to 0
    if((zero_end_point && point.isUnary()) || !point.m_has_entry) {  // IMPORTANT, these must go through 0
      // no need for precision for next few columns, leave it same way as received from AFLOW
      addplot_output_ss << aurostd::PaddedPOST(0, 30); //aurostd::PaddedPOST(aurostd::utype2string(0.0,CHULL_PRECISION,false,ROUNDOFF_TOL,FIXED_STREAM), 30);
      //H_f but in units of kJ/mol, row 5
      addplot_output_ss << aurostd::PaddedPOST(0, 30); //aurostd::PaddedPOST(aurostd::utype2string(0.0,CHULL_PRECISION,false,ROUNDOFF_TOL,FIXED_STREAM), 30);
      // entropic temperature, row 6
      addplot_output_ss << aurostd::PaddedPOST(0, 30); //aurostd::PaddedPOST(aurostd::utype2string(0.0,CHULL_PRECISION,false,ROUNDOFF_TOL,FIXED_STREAM), 30);
    }else{
      // no need for precision for next few columns, leave it same way as received from AFLOW
      num_ss << chull::H_f_atom(point, _m_);
      addplot_output_ss << aurostd::PaddedPOST(num_ss.str(), 30);
      num_ss.str("");
      //H_f but in units of kJ/mol, row 5
      num_ss << meVatom2kJmol * chull::H_f_atom(point, _m_);
      addplot_output_ss << aurostd::PaddedPOST(num_ss.str(), 30);
      num_ss.str("");
      // entropic temperature, row 6
      num_ss << chull::T_S(point);
      addplot_output_ss << aurostd::PaddedPOST(num_ss.str(), 30);
      num_ss.str("");
    }
    if(PRINT_DIST2HULL_COL_TEX){
      // dist_2_hull, row 7
      if(zero_dist_2_hull) {addplot_output_ss << aurostd::PaddedPOST(0, 30);} //aurostd::PaddedPOST(aurostd::utype2string(0.0,CHULL_PRECISION,false,ROUNDOFF_TOL,FIXED_STREAM), 30);
      else{addplot_output_ss << aurostd::PaddedPOST(aurostd::utype2string(point.getDist2Hull(_m_),CHULL_PRECISION), 30);}
    }
    // end line
    addplot_output_ss << endl;
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " done" << endl;}
    return addplot_output_ss.str();
  }

  string ConvexHull::getNodeCoordPosition(const ChullPoint& point) const {
    if(!point.m_has_entry){
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No entry found");
    }
    return getNodeCoordPosition(point.m_entry,point.m_coords);
  }

  string ConvexHull::getNodeCoordPosition(const aflowlib::_aflowlib_entry& entry,const xvector<double>& coord) const {
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    // produces node latex output
    stringstream message;
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " starting" << endl;}

    double sum_C;
    stringstream output;
    // no precision
    stringstream num_ss;

    uint dimension=getDim();
    bool reverse_axes=DEFAULT_CHULL_LATEX_REVERSE_AXIS;

    // first do a check that the function_mode is not out of scope
    if(dimension == 2) {
      if(coord.rows != 2) {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Input dimension mode (2D) and coordinate size mismatch");}
    }else if(dimension == 3) {
      if(coord.rows != 3) {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Input dimension mode (3D) and coordinate size mismatch");}
    }else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Input dimension mode out of scope");}

    output << "axis cs:";  // define axis as relevant "coordinate system"

    if(dimension == 2) {
      // don't need for loop, just 1 coord
      if(!reverse_axes) {output << aurostd::utype2string(1.0 - coord[coord.lrows],FULL_PRECISION);}
      else{output << aurostd::utype2string(coord[coord.lrows],FULL_PRECISION);}
      output << ",";
      if(m_formation_energy_hull) {
        num_ss << chull::H_f_atom(entry, _m_);
        output << num_ss.str();
        num_ss.str("");
      }else{
        num_ss << chull::T_S(entry);
        output << num_ss.str();
        num_ss.str("");
      }
    }else{  // dimension==3
      sum_C = 0.0;
      for(int j=coord.urows-1;j>=coord.lrows;j--) {
        output << aurostd::utype2string(coord[j],FULL_PRECISION);
        output << ",";
        sum_C += coord[j];
      }
      output << aurostd::utype2string(1.0-sum_C,FULL_PRECISION);
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " done" << endl;}
    return output.str();
  }

  string ConvexHull::nodeCreator(stringstream& option, stringstream& position, stringstream& content) const {
    string _option = option.str();
    string _position = position.str();
    string _content = content.str();
    option.str("");
    position.str("");
    content.str("");
    return nodeCreator(_option, _position, _content);
  }

  string ConvexHull::nodeCreator(const string& option, const string& position, const string& content) const {
    // produces node latex output
    stringstream output;
    output << "\\node ";
    if(!option.empty()) {output << "[" << option << "] ";}
    if(!position.empty()) {output << "at (" << position << ") ";}
    if(!content.empty()) {output << "{" << content << "};";}
    output << endl;
    return output.str();
  }

  bool ConvexHull::unwantedFacetLine(uint vi,uint vj,bool check_border) const {  //bool check_border = true;
    vector<vector<uint> > facet_lines;
    return unwantedFacetLine(vi, vj, facet_lines, check_border);
  }

  bool ConvexHull::unwantedFacetLine(uint vi,uint vj,vector<vector<uint> >& facet_lines,bool check_border) const {  //bool check_border = true;
    // checks if the facet created by this combination of chullPoints is necessary
    // for 3D hull
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " starting" << endl;}
    if(vi>m_points.size()-1 || vj>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
    if(check_border) {
      const ChullPoint& ci=m_points[vi];
      const ChullPoint& cj=m_points[vj];
      const xvector<double>& coordAR = ci.getStoichiometricCoords();
      const xvector<double>& coordBR = cj.getStoichiometricCoords();
      bool endpointA = ci.isUnary();
      bool endpointB = cj.isUnary();
      const aflowlib::_aflowlib_entry& entryA = ci.m_entry;
      const aflowlib::_aflowlib_entry& entryB = cj.m_entry;
      // unary to unary
      if(endpointA && endpointB) {
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " unary to unary" << endl;}
        return true;
      }
      // unary to binary, only if dot product of coords is not 0 (unary to crossing binary line)
      if(endpointA && entryB.vspecies.size() == 2 && scalar_product(coordAR, coordBR) >= ZERO_TOL) {
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " unary to binary" << endl;}
        return true;
      }
      // binary to unary, only if dot product of coords is not 0 (unary to crossing binary line)
      if(entryA.vspecies.size() == 2 && endpointB && scalar_product(coordAR, coordBR) >= ZERO_TOL) {
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " binary to unary" << endl;}
        return true;
      }
      // binary to binary, if same species
      if(entryA.vspecies.size() == 2 && entryA.vspecies == entryB.vspecies) {
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " binary to binary" << endl;}
        return true;
      }
    }
    // check that this combination is unique,
    // only pairs inside facet_lines!!!
    // chullPoints, so we're only interested in compositional components of
    // xvector, not the energy
    for(uint i=0,fl_size_i=facet_lines.size();i<fl_size_i;i++) {
      if(facet_lines[i][0]==vi && facet_lines[i][1]==vj){
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " found match! ";
          cerr << facet_lines[i][0] << "==" << vi;
          cerr << ", ";
          cerr << facet_lines[i][1] << "==" << vj;
          cerr << endl;
        }
        return true;
      }
      if(facet_lines[i][0]==vj && facet_lines[i][1]==vi){
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " found match! ";
          cerr << facet_lines[i][0] << "==" << vj;
          cerr << ", ";
          cerr << facet_lines[i][1] << "==" << vi;
          cerr << endl;
        }
        return true;
      }
    }
    facet_lines.push_back(vector<uint>(0));
    facet_lines.back().push_back(vi);
    facet_lines.back().push_back(vj);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " saving facet lines: " << vi << "," << vj << endl;}
    return false;
  }

  string ConvexHull::getPointsPropertyHeaderList(filetype ftype) const {
    bool m_formation_energy_hull=!m_cflags.flag("CHULL::ENTROPIC_TEMPERATURE");    //energy vs. entropic_temperature hull
    bool compounds_column_report=false;
    if(ftype==latex_ft){
      compounds_column_report=DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN; //only grab if necessary, not an inexpensive string search
    }

    string headers="";
    if(!(ftype==latex_ft && !compounds_column_report)){headers+=(!headers.empty()?string(","):string(""))+"compound";}
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"compound_reduced";}
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"compound_reduced_latex";}
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"compound_fractional";}
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"compound_fractional_latex";}
    headers+=(!headers.empty()?string(","):string(""))+"prototype";
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"prototype_latex";}
    headers+=(!headers.empty()?string(","):string(""))+"auid";
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"aurl";}
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"url_entry_page";}
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"nspecies";}
    headers+=(!headers.empty()?string(","):string(""))+"space_group_orig";
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"space_group_orig_latex";}
    headers+=(!headers.empty()?string(","):string(""))+"space_group_relax";
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"space_group_relax_latex";}
    headers+=(!headers.empty()?string(","):string(""))+"spin_atom";
    headers+=(!headers.empty()?string(","):string(""))+"enthalpy_formation_atom";
    headers+=(!headers.empty()?string(","):string(""))+"entropic_temperature";
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"ground_state";}
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"auid_equivalent_structures";}
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"ground_state_icsd";}
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"auid_icsd";}
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"prototype_icsd";}
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"compound_phases_equilibrium";}
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"auid_phases_equilibrium";}
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"compound_phases_decomposition";}
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"auid_phases_decomposition";}
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"coefficient_phases_decomposition";}
    if(m_formation_energy_hull){headers+=(!headers.empty()?string(","):string(""))+"distance_hull_enthalpy_formation_atom";}
    else{headers+=(!headers.empty()?string(","):string(""))+"distance_hull_entropic_temperature";}
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"stability_criterion";}
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"stability_criterion_relative";}
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"N+1_enthalpy_gain";}
    if(ftype==txt_ft || ftype==json_ft){headers+=(!headers.empty()?string(","):string(""))+"entropy_stabilization_coefficient";}

    return headers;
  }

  string ConvexHull::getDelta(bool helvetica_font) const {return (helvetica_font?"\\Updelta":"\\Delta");}
  string ConvexHull::getStabilityCriterionSymbol(bool helvetica_font) const {
    return getDelta(helvetica_font)+" H_{\\mathrm{sc}}";
    //return "\\delta_{\\mathrm{sc}}";
  }

  string ConvexHull::getSnapshotTableHeader(string headers,bool designate_HEADER) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    vector<string> vheaders,vlabels;
    vector<uint> vpaddings;
    aurostd::string2tokens(headers,vheaders,",");
    bool helvetica_font=DEFAULT_CHULL_LATEX_HELVETICA_FONT;

    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]if(0){  //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]  vector<string> valignments_headertable_string;
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]  for(uint i=0,fl_size_i=vheaders.size();i<fl_size_i;i++){
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]    if(vheaders[i]=="compound"){vlabels.push_back("compound"); vpaddings.push_back(80); valignments_headertable_string.push_back("X[3,c,m]");}
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]    else if(vheaders[i]=="prototype"){vlabels.push_back("prototype"); vpaddings.push_back(80); valignments_headertable_string.push_back("X[3,c,m]");}
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]    else if(vheaders[i]=="auid"){vlabels.push_back("auid"); vpaddings.push_back(80); valignments_headertable_string.push_back("X[3,c,m]");}
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]    else if(vheaders[i]=="space_group_orig"){vlabels.push_back("original space group"); vpaddings.push_back(30); valignments_headertable_string.push_back("X[2,c,m]");} //SG original
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]    else if(vheaders[i]=="space_group_relax"){vlabels.push_back("relaxed space group"); vpaddings.push_back(30); valignments_headertable_string.push_back("X[2,c,m]");} //SG relax
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]    else if(vheaders[i]=="spin_atom"){vlabels.push_back("spin $\\left(\\mu_{\\mathrm{B}}\\mathrm{/atom}\\right)$"); vpaddings.push_back(30); valignments_headertable_string.push_back("X[2,c,m]");}
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]    else if(vheaders[i]=="enthalpy_formation_atom"){vlabels.push_back("$H_{\\mathrm{f}}$ (meV/atom)"); vpaddings.push_back(30); valignments_headertable_string.push_back("X[2,c,m]");}
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]    else if(vheaders[i]=="entropic_temperature"){vlabels.push_back("$T_{\\mathrm{S}}$ (K)"); vpaddings.push_back(30); valignments_headertable_string.push_back("X[2,c,m]");}
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]    else if(vheaders[i]=="distance_hull_enthalpy_formation_atom"){vlabels.push_back("$"+getDelta(helvetica_font)+" H_{\\mathrm{hull}}$ (meV/atom)"); vpaddings.push_back(30); valignments_headertable_string.push_back("X[2,c,m]");}
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]    else if(vheaders[i]=="distance_hull_entropic_temperature"){vlabels.push_back("$"+getDelta(helvetica_font)+" T_{\\mathrm{S}}$ (K)"); vpaddings.push_back(30); valignments_headertable_string.push_back("X[2,c,m]");}
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]    else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown property");}
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]  }
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]}

    vector<uint> valignments_headertable_uint;
    uint total_alignment=0;
    for(uint i=0,fl_size_i=vheaders.size();i<fl_size_i;i++){
      if(vheaders[i]=="compound"){vlabels.push_back("compound"); vpaddings.push_back(80); valignments_headertable_uint.push_back(3); total_alignment+=valignments_headertable_uint.back();}
      else if(vheaders[i]=="prototype"){vlabels.push_back("prototype"); vpaddings.push_back(80); valignments_headertable_uint.push_back(3); total_alignment+=valignments_headertable_uint.back();}
      else if(vheaders[i]=="auid"){vlabels.push_back("auid"); vpaddings.push_back(80); valignments_headertable_uint.push_back(3); total_alignment+=valignments_headertable_uint.back();}
      else if(vheaders[i]=="space_group_orig"){vlabels.push_back("original space group"); vpaddings.push_back(30); valignments_headertable_uint.push_back(2); total_alignment+=valignments_headertable_uint.back();} //SG original
      else if(vheaders[i]=="space_group_relax"){vlabels.push_back("relaxed space group"); vpaddings.push_back(30); valignments_headertable_uint.push_back(2); total_alignment+=valignments_headertable_uint.back();} //SG relax
      else if(vheaders[i]=="spin_atom"){vlabels.push_back("spin $\\left(\\mu_{\\mathrm{B}}\\mathrm{/atom}\\right)$"); vpaddings.push_back(30); valignments_headertable_uint.push_back(2); total_alignment+=valignments_headertable_uint.back();}
      else if(vheaders[i]=="enthalpy_formation_atom"){vlabels.push_back("$H_{\\mathrm{f}}$ (meV/atom)"); vpaddings.push_back(30); valignments_headertable_uint.push_back(2); total_alignment+=valignments_headertable_uint.back();}
      else if(vheaders[i]=="entropic_temperature"){vlabels.push_back("$T_{\\mathrm{S}}$ (K)"); vpaddings.push_back(30); valignments_headertable_uint.push_back(2); total_alignment+=valignments_headertable_uint.back();}
      else if(vheaders[i]=="distance_hull_enthalpy_formation_atom"){vlabels.push_back("$"+getDelta(helvetica_font)+" H_{\\mathrm{hull}}$ (meV/atom)"); vpaddings.push_back(30); valignments_headertable_uint.push_back(2); total_alignment+=valignments_headertable_uint.back();}
      else if(vheaders[i]=="distance_hull_entropic_temperature"){vlabels.push_back("$"+getDelta(helvetica_font)+" T_{\\mathrm{S}}$ (K)"); vpaddings.push_back(30); valignments_headertable_uint.push_back(2); total_alignment+=valignments_headertable_uint.back();}
      else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown property");}
    }

    vector<string> vlabels_styled;
    for(uint i=0,fl_size_i=vlabels.size();i<fl_size_i;i++){vlabels_styled.push_back(aurostd::PaddedPOST("{\\footnotesize "+vlabels[i]+"}", vpaddings[i]));}  //used to be small before auid

    double page_width=(LATEX_WIDTH_LETTER_STD-(LATEX_LEFT_MARGIN_LETTER_STD+LATEX_RIGHT_MARGIN_LETTER_STD));
    double penalty_tabcolsep=2*LATEX_TABCOLSEP_STD*LATEX_PT2INCH; //there are two colsep for reach column, then convert to inches
    double penalty_arrayrulewidth=(valignments_headertable_uint.size()+1)*LATEX_ARRAYRULEWIDTH_STD*LATEX_PT2INCH/(valignments_headertable_uint.size()); //there are columns+1 separators, then convert to inches, then convert to per column penalty

    vector<string> valignments_headertable_string;
    double page_width_fraction=0.0;
    for(uint i=0,fl_size_i=valignments_headertable_uint.size();i<fl_size_i;i++){
      page_width_fraction=page_width*(double)valignments_headertable_uint[i]/(double)total_alignment;
      valignments_headertable_string.push_back("X{"+aurostd::utype2string( page_width_fraction - penalty_tabcolsep - penalty_arrayrulewidth )+"in}"); 
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " " << valignments_headertable_string.back() << endl;}
    }

    stringstream output;
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]output << "\\begin{tabu}{|" << aurostd::joinWDelimiter(valignments_headertable_string,"|") << "|}" << (designate_HEADER?" \%HEADER":"") << endl;
    output << "\\begin{tabular}{|" << aurostd::joinWDelimiter(valignments_headertable_string,"|") << "|}" << (designate_HEADER?" \%HEADER":"") << endl;
    output << "\\toprule" << (designate_HEADER?" \%HEADER":"") << endl;
    output << aurostd::PaddedPOST("\\rowcolor{white}", 30) << " " << aurostd::joinWDelimiter(vlabels_styled," & ") << " \\\\" << (designate_HEADER?" \%HEADER":"") << endl;
    output << "\\midrule" << (designate_HEADER?" \%HEADER":"") << endl;
    output << "\\end{tabular}" << (designate_HEADER?" \%HEADER":"") << endl;
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]output << "\\end{tabu}" << (designate_HEADER?" \%HEADER":"") << endl;

    return output.str();
  }

  bool ConvexHull::addInternalHyperlinks(bool internal_links_graph2report,bool internal_links_withinreport) const {
    int links_setting=DEFAULT_CHULL_LATEX_LINKS;
    bool kill_all_links=(links_setting==0); //no links whatsoever
    bool no_internal_links=(links_setting==2);       //no jumping
    bool image_only=m_cflags.flag("CHULL::IMAGE_ONLY");
    if(image_only){kill_all_links=true;}
    if(kill_all_links){no_internal_links=true;}      //no jumping
    if((!internal_links_graph2report)&&(!internal_links_withinreport)){no_internal_links=true;}
    return !no_internal_links;
  }

  bool ConvexHull::addExternalHyperlinks() const {
    int links_setting=DEFAULT_CHULL_LATEX_LINKS;
    bool kill_all_links=(links_setting==0); //no links whatsoever
    bool no_external_links=(links_setting==3);  //no weblinks
    bool image_only=m_cflags.flag("CHULL::IMAGE_ONLY");
    if(image_only){kill_all_links=true;}
    if(kill_all_links){no_external_links=true;} //no weblinks
    return !no_external_links;
  }

  double ConvexHull::getRoundToValue(double point_range) const {
    int order_of_mag = log10(point_range);
    if(order_of_mag<1){order_of_mag=1;} //we need to work with integers since units are meV/atom, anything less is really just noise
    double round_to_value=pow(10, order_of_mag-1) * 2.5;  //round to nearest 2.5, 25, 250
    return round_to_value;
  }

  double ConvexHull::getYTickDistance(double y_range,int approx_num_ticks,double round_to_value) const {
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " y_range=" << y_range << endl;
      cerr << __AFLOW_FUNC__ << " approx_num_ticks=" << approx_num_ticks << endl;
      cerr << __AFLOW_FUNC__ << " round_to_value=" << round_to_value << endl;
    }
    double delta,delta_pow;
    delta=delta_pow=y_range/(double)approx_num_ticks;
    int exponent=0;
    while(delta_pow<round_to_value){exponent+=1;delta_pow=delta*pow(10,exponent);}
    delta=delta_pow;
    int y_tick_distance_int;
    y_tick_distance_int=aurostd::roundDouble(delta,(int)round_to_value,false);  //divisor is ~ # of ticks  //CO20190724 - explicit double->int conversion (floor is fine)
    if(y_tick_distance_int==0){ //aflow_BSm_hull.pdf, we need some NONZERO y_tick_distance_int
      y_tick_distance_int=aurostd::roundDouble(delta,(int)round_to_value,true);  //divisor is ~ # of ticks //CO20190724 - explicit double->int conversion (floor is fine)
    }
    double y_tick_distance=y_tick_distance_int*pow(10,-exponent);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " y_tick_distance=" << y_tick_distance << endl;
    }
    return y_tick_distance;
  }

  vector<string> ConvexHull::grabAcceptableLatexColors(bool replace_pranab_standard,bool allow_dvips_colors,uint count) const {
    vector<string> banned_colors;
    return grabAcceptableLatexColors(banned_colors,replace_pranab_standard,allow_dvips_colors,count);
  }
  vector<string> ConvexHull::grabAcceptableLatexColors(const string& banned_colors_str,bool replace_pranab_standard,bool allow_dvips_colors,uint count) const {
    vector<string> banned_colors;
    aurostd::string2tokens(banned_colors_str,banned_colors,",");
    return grabAcceptableLatexColors(banned_colors,replace_pranab_standard,allow_dvips_colors,count);
  }
  vector<string> ConvexHull::grabAcceptableLatexColors(const vector<string>& banned_colors,bool replace_pranab_standard,bool allow_dvips_colors,uint count) const {
    vector<string> _latex_colors,latex_colors;
    string color,lower_color,colors=LATEX_DEFAULT_COLORS;
    uint loop=0;
    while(true){
      aurostd::string2tokens(colors,_latex_colors,",");
      for(uint i=0,fl_size_i=_latex_colors.size();i<fl_size_i;i++){
        color=_latex_colors[i];
        lower_color=aurostd::tolower(color);
        if(aurostd::WithinList(banned_colors,lower_color)){continue;}
        if(replace_pranab_standard&&(lower_color=="green"||lower_color=="red")){color="pranab_"+lower_color;}  //now in pranab standard!
        latex_colors.push_back(color);
        if(latex_colors.size()==count){return latex_colors;}
      }
      if(allow_dvips_colors&&(loop++%2==0)){colors=LATEX_DVIPS_COLORS;}
      else{colors=LATEX_DEFAULT_COLORS;}
    }
    return latex_colors;
  }

  aurostd::xoption ConvexHull::resolvePlotLabelSettings() const {
    stringstream message;

    bool labels_off_hull=DEFAULT_CHULL_LATEX_LABELS_OFF_HULL;
    bool meta_labels=DEFAULT_CHULL_LATEX_META_LABELS;
    bool compound_labels=true;  //default is to show compound labels only, for binaries this will be ground-state dependent
    bool prototype_labels=(getDim()==2 && labels_off_hull); //default is to show prototype labels for binaries, with compounds labels on ground-states for reference
    bool icsd_labels=false;
    bool no_labels=false;
    string label_options=DEFAULT_CHULL_LATEX_LABEL_NAME;
    bool plot_labels=(!label_options.empty());           //overarching flag, mostly for getting options
    if(plot_labels){
      compound_labels=prototype_labels=icsd_labels=no_labels=false; //kill defaults
      vector<string> vlabelstring;
      aurostd::string2tokens(label_options,vlabelstring, ",");
      for(uint i=0,fl_size_i=vlabelstring.size();i<fl_size_i;i++) {
        if(vlabelstring[i][0] == 'B' || vlabelstring[i][0] == 'b') {  // both
          compound_labels=true;
          prototype_labels=true;
          no_labels=false;
          break;
        }else if(vlabelstring[i][0] == 'I' || vlabelstring[i][0] == 'i') {  // ICSD
          icsd_labels=true;
          no_labels=false;
          if(!compound_labels&&!prototype_labels){compound_labels=true;}  //default for icsd
        }else if(vlabelstring[i][0] == 'N' || vlabelstring[i][0] == 'n' || vlabelstring[i][0] == 'O' || vlabelstring[i][0] == 'o') {  // none,off
          compound_labels=false;
          prototype_labels=false;
          icsd_labels=false;
          meta_labels=false;
          no_labels=true;
          break;
        }else if(vlabelstring[i][0] == 'C' || vlabelstring[i][0] == 'c') {  // compound
          compound_labels=true;
          no_labels=false;
        }else if(vlabelstring[i][0] == 'P' || vlabelstring[i][0] == 'p') {  // prototype
          prototype_labels=true;
          no_labels=false;
        }else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Incorrect input for plot labels \""+vlabelstring[i]+"\"");}
      }
    }
    if(no_labels){
      if(labels_off_hull){
        //[verbose once in writeLatex()]message << "LABEL_NAME set to NONE but LABELS_OFF_HULL requested (fix .aflow.rc), toggling LABELS_OFF_HULL off"
        //[verbose once in writeLatex()]pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        labels_off_hull=false;
      }
      if(meta_labels){
        //[verbose once in writeLatex()]message << "LABEL_NAME set to NONE but META_LABELS requested (fix .aflow.rc), toggling META_LABELS off"
        //[verbose once in writeLatex()]pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        meta_labels=false;
      }
    }
    aurostd::xoption lflags;
    lflags.flag("CHULL::PLOT_LABELS",plot_labels);
    lflags.flag("CHULL::PLOT_LABELS_OFF_HULL",labels_off_hull);
    lflags.flag("CHULL::PLOT_META_LABELS",meta_labels);
    lflags.flag("CHULL::PLOT_COMPOUND_LABELS",compound_labels);
    lflags.flag("CHULL::PLOT_PROTOTYPE_LABELS",prototype_labels);
    lflags.flag("CHULL::PLOT_ICSD_LABELS",icsd_labels);
    lflags.flag("CHULL::PLOT_NO_LABELS",no_labels);
    return lflags;
  }

  void ConvexHull::writeLatex() const {
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    stringstream message;
    if(!aurostd::IsCommandAvailable("pdflatex")) {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"\"pdflatex\" needs to be in your path");}
    if(m_cflags.flag("CHULL::PNG_IMAGE")){
      if(!aurostd::IsCommandAvailable("convert")) {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"\"convert\" needs to be in your path");}
    }
    message << "Starting LaTeX PDF generator";
    pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

    //////////////////////////////////////////////////////////////////////////////
    // START Getting hull attributes
    //////////////////////////////////////////////////////////////////////////////

    uint dimension=getDim();
    vector<uint> hull_points=getHullPoints(false);
    bool draw_all_facet_lines=false;  //only draw uniques, keep FALSE
    bool include_equilibrium_phases=true; //keep true

    //USER INPUTS
    bool doc_only=m_cflags.flag("CHULL::DOC_ONLY");
    bool no_doc=m_cflags.flag("CHULL::NO_DOC");
    bool image_only=m_cflags.flag("CHULL::IMAGE_ONLY");
    bool light_contrast=m_cflags.flag("CHULL::LIGHT_CONTRAST");
    bool large_font=m_cflags.flag("CHULL::LARGE_FONT");
    bool keep_tex=m_cflags.flag("CHULL::KEEP_TEX");
    bool show_latex_output=m_cflags.flag("CHULL::LATEX_OUTPUT");
    bool latex_interactive_mode=m_cflags.flag("CHULL::LATEX_INTERACTIVE");

    //FROM AFLOWRC
    bool print_aflow_logo_full=_AFLOW_CHULL_PRINT_LOGO_1;  //otherwise print aflow_text
    bool print_aflow_text_logo=false;
    if(print_aflow_logo_full&&print_aflow_text_logo){print_aflow_text_logo=false;}
    bool print_logo_2=_AFLOW_CHULL_PRINT_LOGO_2;
    bool print_aflow_webaddress_logo=true;
    if(print_logo_2&&print_aflow_webaddress_logo){print_aflow_webaddress_logo=false;}
    int banner_setting=DEFAULT_CHULL_LATEX_BANNER;
    bool no_banner=(banner_setting==0);
    bool small_banner=(banner_setting==2||image_only);
    if(no_banner&&small_banner){small_banner=false;}
    if(small_banner&&(image_only||no_doc)){print_aflow_logo_full=print_logo_2=false;}
    bool compounds_column_report=DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN;
    bool stoich_header_report=DEFAULT_CHULL_LATEX_STOICH_HEADER;
    bool plot_unaries=DEFAULT_CHULL_LATEX_PLOT_UNARIES;
    string filter_scheme=DEFAULT_CHULL_LATEX_FILTER_SCHEME;
    bool filter_by_z=(!filter_scheme.empty() && (aurostd::toupper(filter_scheme)[0]=='Z' || aurostd::toupper(filter_scheme)[0]=='E'));
    bool filter_by_distance=(!filter_scheme.empty() && aurostd::toupper(filter_scheme)[0]=='D');
    double filter_cutoff=DEFAULT_CHULL_LATEX_FILTER_VALUE;

    int plot_off_hull_setting=DEFAULT_CHULL_LATEX_PLOT_OFF_HULL;     //does not include unstable
    bool plot_off_hull;
    if(plot_off_hull_setting==-1){plot_off_hull=(getDim()==2?true:false);}
    else{plot_off_hull=(plot_off_hull_setting==0?false:true);}
    if(filter_by_z||filter_by_distance){plot_off_hull=true;}

    bool plot_unstable=DEFAULT_CHULL_LATEX_PLOT_UNSTABLE;
    bool reverse_axes=DEFAULT_CHULL_LATEX_REVERSE_AXIS;

    bool display_color_gradient=DEFAULT_CHULL_LATEX_COLOR_GRADIENT;
    bool include_color_bar=DEFAULT_CHULL_LATEX_COLOR_BAR;
    bool show_heat_map=DEFAULT_CHULL_LATEX_HEAT_MAP;
    if(!display_color_gradient){
      include_color_bar=false;
      show_heat_map=false;
    }

    bool hull_drop_shadow=DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW;  //gus paper

    string ternary_label_color_setting=DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR;
    string ternary_label_color=(ternary_label_color_setting.empty()?"white":ternary_label_color_setting);
    if(hull_drop_shadow&&ternary_label_color=="white"){ternary_label_color="yellow";} //white doesn't work here
    if(!show_heat_map&&(ternary_label_color=="white"||ternary_label_color=="yellow")){ternary_label_color="black";} //white/yellow doesn't work here

    string color_map_setting=DEFAULT_CHULL_LATEX_COLOR_MAP;
    string color_map=(color_map_setting.empty()?"rgb(0pt)=(0,0,1); rgb(63pt)=(1,0.644,0)":color_map_setting);

    //labels options
    aurostd::xoption lflags=resolvePlotLabelSettings();
    bool plot_labels=lflags.flag("CHULL::PLOT_LABELS");
    bool labels_off_hull=lflags.flag("CHULL::PLOT_LABELS_OFF_HULL");
    bool meta_labels=lflags.flag("CHULL::PLOT_META_LABELS");
    bool compound_labels=lflags.flag("CHULL::PLOT_COMPOUND_LABELS");
    bool prototype_labels=lflags.flag("CHULL::PLOT_PROTOTYPE_LABELS");
    bool icsd_labels=lflags.flag("CHULL::PLOT_ICSD_LABELS");
    bool no_labels=lflags.flag("CHULL::PLOT_NO_LABELS");
    //[MOVED to resolvePlotLabelSettings()]bool labels_off_hull=DEFAULT_CHULL_LATEX_LABELS_OFF_HULL;
    //[MOVED to resolvePlotLabelSettings()]bool meta_labels=DEFAULT_CHULL_LATEX_META_LABELS;
    //[MOVED to resolvePlotLabelSettings()]bool compound_labels=true;  //default is to show compound labels only, for binaries this will be ground-state dependent
    //[MOVED to resolvePlotLabelSettings()]bool prototype_labels=(getDim()==2 && labels_off_hull); //default is to show prototype labels for binaries, with compounds labels on ground-states for reference
    //[MOVED to resolvePlotLabelSettings()]bool icsd_labels=false;
    //[MOVED to resolvePlotLabelSettings()]bool no_labels=false;
    //[MOVED to resolvePlotLabelSettings()]string label_options=DEFAULT_CHULL_LATEX_LABEL_NAME;
    //[MOVED to resolvePlotLabelSettings()]bool plot_labels=(!label_options.empty());           //overarching flag, mostly for getting options
    //[MOVED to resolvePlotLabelSettings()]if(plot_labels){
    //[MOVED to resolvePlotLabelSettings()]  compound_labels=prototype_labels=icsd_labels=no_labels=false; //kill defaults
    //[MOVED to resolvePlotLabelSettings()]  vector<string> vlabelstring;
    //[MOVED to resolvePlotLabelSettings()]  aurostd::string2tokens(label_options,vlabelstring, ",");
    //[MOVED to resolvePlotLabelSettings()]  for(uint i=0,fl_size_i=vlabelstring.size();i<fl_size_i;i++) {
    //[MOVED to resolvePlotLabelSettings()]    if(vlabelstring[i][0] == 'B' || vlabelstring[i][0] == 'b') {  // both
    //[MOVED to resolvePlotLabelSettings()]      compound_labels=true;
    //[MOVED to resolvePlotLabelSettings()]      prototype_labels=true;
    //[MOVED to resolvePlotLabelSettings()]      icsd_labels=false;
    //[MOVED to resolvePlotLabelSettings()]      no_labels=false;
    //[MOVED to resolvePlotLabelSettings()]      break;
    //[MOVED to resolvePlotLabelSettings()]    }else if(vlabelstring[i][0] == 'I' || vlabelstring[i][0] == 'i') {  // ICSD
    //[MOVED to resolvePlotLabelSettings()]      compound_labels=false;
    //[MOVED to resolvePlotLabelSettings()]      prototype_labels=false;
    //[MOVED to resolvePlotLabelSettings()]      icsd_labels=true;
    //[MOVED to resolvePlotLabelSettings()]      no_labels=false;
    //[MOVED to resolvePlotLabelSettings()]      break;
    //[MOVED to resolvePlotLabelSettings()]    }else if(vlabelstring[i][0] == 'N' || vlabelstring[i][0] == 'n' || vlabelstring[i][0] == 'O' || vlabelstring[i][0] == 'o') {  // none,off
    //[MOVED to resolvePlotLabelSettings()]      compound_labels=false;
    //[MOVED to resolvePlotLabelSettings()]      prototype_labels=false;
    //[MOVED to resolvePlotLabelSettings()]      icsd_labels=false;
    //[MOVED to resolvePlotLabelSettings()]      meta_labels=false;
    //[MOVED to resolvePlotLabelSettings()]      no_labels=true;
    //[MOVED to resolvePlotLabelSettings()]      break;
    //[MOVED to resolvePlotLabelSettings()]    }else if(vlabelstring[i][0] == 'C' || vlabelstring[i][0] == 'c') {  // compound
    //[MOVED to resolvePlotLabelSettings()]      compound_labels=true;
    //[MOVED to resolvePlotLabelSettings()]      icsd_labels=false;
    //[MOVED to resolvePlotLabelSettings()]      no_labels=false;
    //[MOVED to resolvePlotLabelSettings()]    }else if(vlabelstring[i][0] == 'P' || vlabelstring[i][0] == 'p') {  // prototype
    //[MOVED to resolvePlotLabelSettings()]      prototype_labels=true;
    //[MOVED to resolvePlotLabelSettings()]      icsd_labels=false;
    //[MOVED to resolvePlotLabelSettings()]      no_labels=false;
    //[MOVED to resolvePlotLabelSettings()]    }else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Incorrect input for plot labels \""+vlabelstring[i]+"\"");}
    //[MOVED to resolvePlotLabelSettings()]  }
    //[MOVED to resolvePlotLabelSettings()]}
    if(no_labels){
      if(labels_off_hull){
        message << "LABEL_NAME set to NONE but LABELS_OFF_HULL requested (fix .aflow.rc), toggling LABELS_OFF_HULL off";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        //[turned off already in resolvePlotLabelSettings()]labels_off_hull=false;
      }
      if(meta_labels){
        message << "LABEL_NAME set to NONE but META_LABELS requested (fix .aflow.rc), toggling META_LABELS off";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        //[turned off already in resolvePlotLabelSettings()]meta_labels=false;
      }
    }
    if(labels_off_hull){plot_off_hull=true;}

    int plot_reduced_composition_setting=DEFAULT_CHULL_LATEX_PLOT_REDUCED_COMPOSITION;
    bool plot_reduced_composition;
    if(plot_reduced_composition_setting==-1){plot_reduced_composition=(!(getDim()==2&&compound_labels&&plot_labels));}
    else{plot_reduced_composition=(plot_reduced_composition_setting==0?false:true);}
    if(no_labels){
      if(plot_reduced_composition){
        message << "LABEL_NAME set to NONE but PLOT_REDUCED_COMPOSITION requested (fix .aflow.rc), toggling PLOT_REDUCED_COMPOSITION off";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        plot_reduced_composition=false;
      }
    }
    bool helvetica_font=DEFAULT_CHULL_LATEX_HELVETICA_FONT;
    string font_size=DEFAULT_CHULL_LATEX_FONT_SIZE;
    bool rotate_labels=DEFAULT_CHULL_LATEX_ROTATE_LABELS;
    if(no_labels){
      if(rotate_labels){
        message << "LABEL_NAME set to NONE but ROTATE_LABELS requested (fix .aflow.rc), toggling ROTATE_LABELS off";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        rotate_labels=false;
      }
    }

    //new default option here - -1 means NO BOLD except ternaries
    //bool bold_labels=DEFAULT_CHULL_LATEX_BOLD_LABELS;
    int bold_labels_setting=DEFAULT_CHULL_LATEX_BOLD_LABELS;
    bool bold_labels=(bold_labels_setting==1);
    bool bold_labels_ternaries=(bold_labels || bold_labels_setting==-1);

    //override with user inputs
    if(light_contrast){color_map="rgb(0pt)=(0.035,0.270,0.809); rgb(63pt)=(1,0.644,0)";}
    if(font_size.empty()){
      if(large_font) {
        if(helvetica_font) {font_size="huge";}
        else{font_size="Large";}
      }else{
        if(image_only){font_size="LARGE";}
        else{font_size="large";} //safely, I can do large
      }
    }

    double opacity_watermark=0.3;
    string watermark_x_shift="2.5in";

    //////////////////////////////////////////////////////////////////////////////
    // START Getting hull attributes
    //////////////////////////////////////////////////////////////////////////////

    if(dimension>3) {
      doc_only=true;
      message << "CHULL::DOC_ONLY set to TRUE (dimension>3)";
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_);
    }

    // initializing stringstreams to use
    stringstream doc_header_TEX_ss, _doc_header_TEX_ss;               // latex header commands
    stringstream tikzpic_settings_TEX_ss,_tikzpic_settings_TEX_ss;    // settings in tikzpic
    stringstream tikzpic_TEX_ss;                         // tikzpicture commands
    stringstream convex_hull_facets_TEX_ss;              // never empty
    stringstream convex_hull_vertices_TEX_ss;            // never empty
    stringstream convex_hull_facets_drop_shadow_TEX_ss;  // never empty
    stringstream common_settings_TEX_ss;                    // use for common settings to reduce redundancy
    stringstream pseudo_preliminary_axes_TEX_ss;               // unary tikzpicture label commands
    stringstream heat_map_TEX_ss;                   // heatmap stuff
    stringstream points_data_ss;                    // points commands
    stringstream labels_data_ss;
    stringstream report_data_ss, _report_data_ss;
    stringstream equilibrium_phases_TEX_ss;
    stringstream equilibrium_phases_header_TEX_ss;
    stringstream reaction_chem_eq_TEX_ss;
    stringstream main_TEX_ss;
    stringstream node_option_ss;
    stringstream node_position_ss;
    stringstream node_content_ss;
    stringstream misc_ss;
    // no precision
    stringstream num_ss;

    // initializing some strings
    string main_TEX_file = "", main_PDF_file = "", main_PNG_file = "";
    string main_file = "", main_output_file = "";
    string input = "",input_hyphened = "";

    // creating name of output file
    input=aurostd::joinWDelimiter(m_velements,"");
    input_hyphened=aurostd::joinWDelimiter(m_velements,"-");
    main_file="aflow_"+input+"_hull";
    if(image_only) {main_TEX_file = main_file + "_IMAGEONLY.tex";}
    else{main_TEX_file = main_file + ".tex";}
    main_PDF_file = main_file + ".pdf";
    main_PNG_file = main_file + ".png";
    string aflow_logo_full_file = "aflow_logo_full.pdf";
    string aflow_logo_skinny_file = "aflow_logo_skinny.pdf";
    string logo_file_2 = "logo2.png";
    string DEFAULT_PLOT_COLUMN_HEADER = "";
    if(m_formation_energy_hull){DEFAULT_PLOT_COLUMN_HEADER="H_f_meVatom";}
    else{DEFAULT_PLOT_COLUMN_HEADER="T_S";}

    // other initialization
    uint plot_points_count = 0;                // points to put on ternary plot
    uint plot_points_count_no_end_points = 0;  // mostly for count purposes
    vector<uint> chull_points;
    vector<vector<uint> > facet_lines, facet_lines_dropshadow;
    double min_point = 0.0, max_point = 0.0, point_range = 0.0;  // saves min/max energy value to determine if we can
    // have a colorbar
    double z_filter_cutoff = 0.0, dist_filter_cutoff = 0.0;
    string plot_command = "";
    string output_name = "";
    string misc = "";
    vector<string> files_2_move; //, sg_tokens;
    stringstream command;
    uint num_horizontal_planes = 0;  // to determine whether or not we should have heatmaps
    //uint count_hull_entries=getEntriesCount(m_half_hull);
    uint count_total_entries=getEntriesCount(false);
    uint i_point=AUROSTD_MAX_UINT,i_coord_group=AUROSTD_MAX_UINT;
    bool added_header;

    string general_image_font_size="\\Large";
    string label_image_font_size="\\LARGE";
    if(image_only){
      //https://tex.stackexchange.com/questions/24599/what-point-pt-font-size-are-large-etc
      general_image_font_size="\\Huge";
      label_image_font_size="\\fontsize{28}{32}\\selectfont";
    }

    string MARGIN_PICTURE="includeheadfoot,headheight=16pt,margin=0.5in"; //do not make margin any smaller, will screw up fancyheaders elsewhere
    string MARGIN_REPORT="includeheadfoot,headheight="+(print_aflow_logo_full?string("70"):string("50"))+
      "pt,headsep=0.1in,top=0.5in,bottom=0.75in,left="+aurostd::utype2string(LATEX_LEFT_MARGIN_LETTER_STD)+"in,right="+aurostd::utype2string(LATEX_RIGHT_MARGIN_LETTER_STD)+"in,footskip=0.5in";

    string headers=getPointsPropertyHeaderList(latex_ft);
    vector<string> vheaders;
    aurostd::string2tokens(headers,vheaders,",");

    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]//get alignments
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]vector<string> valignments;
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]for(uint i=0,fl_size_i=vheaders.size();i<fl_size_i;i++){
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]  if(vheaders[i]=="compound"){valignments.push_back("X[3,l,m]");}
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]  else if(vheaders[i]=="prototype"){valignments.push_back("X[3,l,m]");}
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]  else if(vheaders[i]=="auid"){valignments.push_back("X[3,l,m]");}
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]  else if(vheaders[i]=="space_group_orig"){valignments.push_back("X[2,l,m]");}
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]  else if(vheaders[i]=="space_group_relax"){valignments.push_back("X[2,l,m]");}
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]  else if(vheaders[i]=="spin_atom"){valignments.push_back("X[2,r,m]");}
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]  else if(vheaders[i]=="enthalpy_formation_atom"){valignments.push_back("X[2,r,m]");}
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]  else if(vheaders[i]=="entropic_temperature"){valignments.push_back("X[2,r,m]");}
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]  else if(vheaders[i]=="distance_hull_enthalpy_formation_atom"){valignments.push_back("X[2,r,m]");}
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]  else if(vheaders[i]=="distance_hull_entropic_temperature"){valignments.push_back("X[2,r,m]");}
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]  else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown property");}
    //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]}

    //get alignments
    vector<uint> valignments_entrytable_uint;
    uint total_alignment_entrytable=0;
    for(uint i=0,fl_size_i=vheaders.size();i<fl_size_i;i++){
      if(vheaders[i]=="compound"){valignments_entrytable_uint.push_back(3); total_alignment_entrytable+=valignments_entrytable_uint.back();}
      else if(vheaders[i]=="prototype"){valignments_entrytable_uint.push_back(3); total_alignment_entrytable+=valignments_entrytable_uint.back();}
      else if(vheaders[i]=="auid"){valignments_entrytable_uint.push_back(3); total_alignment_entrytable+=valignments_entrytable_uint.back();}
      else if(vheaders[i]=="space_group_orig"){valignments_entrytable_uint.push_back(2); total_alignment_entrytable+=valignments_entrytable_uint.back();}
      else if(vheaders[i]=="space_group_relax"){valignments_entrytable_uint.push_back(2); total_alignment_entrytable+=valignments_entrytable_uint.back();}
      else if(vheaders[i]=="spin_atom"){valignments_entrytable_uint.push_back(2); total_alignment_entrytable+=valignments_entrytable_uint.back();}
      else if(vheaders[i]=="enthalpy_formation_atom"){valignments_entrytable_uint.push_back(2); total_alignment_entrytable+=valignments_entrytable_uint.back();}
      else if(vheaders[i]=="entropic_temperature"){valignments_entrytable_uint.push_back(2); total_alignment_entrytable+=valignments_entrytable_uint.back();}
      else if(vheaders[i]=="distance_hull_enthalpy_formation_atom"){valignments_entrytable_uint.push_back(2); total_alignment_entrytable+=valignments_entrytable_uint.back();}
      else if(vheaders[i]=="distance_hull_entropic_temperature"){valignments_entrytable_uint.push_back(2); total_alignment_entrytable+=valignments_entrytable_uint.back();}
      else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown property");}
    }

    double page_width=(LATEX_WIDTH_LETTER_STD-(LATEX_LEFT_MARGIN_LETTER_STD+LATEX_RIGHT_MARGIN_LETTER_STD));
    double penalty_tabcolsep=2*LATEX_TABCOLSEP_STD*LATEX_PT2INCH; //there are two colsep for reach column, then convert to inches
    double penalty_arrayrulewidth=(valignments_entrytable_uint.size()+1)*LATEX_ARRAYRULEWIDTH_STD*LATEX_PT2INCH/(valignments_entrytable_uint.size()); //there are columns+1 separators, then convert to inches, then convert to per column penalty

    vector<string> valignments_entrytable_string;
    double page_width_fraction=0.0;
    for(uint i=0,fl_size_i=vheaders.size();i<fl_size_i;i++){
      page_width_fraction=page_width*(double)valignments_entrytable_uint[i]/(double)total_alignment_entrytable;
      if(vheaders[i]=="compound"){valignments_entrytable_string.push_back("L{"+aurostd::utype2string( page_width_fraction - penalty_tabcolsep - penalty_arrayrulewidth )+"in}");}
      else if(vheaders[i]=="prototype"){valignments_entrytable_string.push_back("L{"+aurostd::utype2string( page_width_fraction - penalty_tabcolsep - penalty_arrayrulewidth )+"in}");}
      else if(vheaders[i]=="auid"){valignments_entrytable_string.push_back("L{"+aurostd::utype2string( page_width_fraction - penalty_tabcolsep - penalty_arrayrulewidth )+"in}");}
      else if(vheaders[i]=="space_group_orig"){valignments_entrytable_string.push_back("L{"+aurostd::utype2string( page_width_fraction - penalty_tabcolsep - penalty_arrayrulewidth )+"in}");}
      else if(vheaders[i]=="space_group_relax"){valignments_entrytable_string.push_back("L{"+aurostd::utype2string( page_width_fraction - penalty_tabcolsep - penalty_arrayrulewidth )+"in}");}
      else if(vheaders[i]=="spin_atom"){valignments_entrytable_string.push_back("R{"+aurostd::utype2string( page_width_fraction - penalty_tabcolsep - penalty_arrayrulewidth )+"in}");}
      else if(vheaders[i]=="enthalpy_formation_atom"){valignments_entrytable_string.push_back("R{"+aurostd::utype2string( page_width_fraction - penalty_tabcolsep - penalty_arrayrulewidth )+"in}");}
      else if(vheaders[i]=="entropic_temperature"){valignments_entrytable_string.push_back("R{"+aurostd::utype2string( page_width_fraction - penalty_tabcolsep - penalty_arrayrulewidth )+"in}");}
      else if(vheaders[i]=="distance_hull_enthalpy_formation_atom"){valignments_entrytable_string.push_back("R{"+aurostd::utype2string( page_width_fraction - penalty_tabcolsep - penalty_arrayrulewidth )+"in}");}
      else if(vheaders[i]=="distance_hull_entropic_temperature"){valignments_entrytable_string.push_back("R{"+aurostd::utype2string( page_width_fraction - penalty_tabcolsep - penalty_arrayrulewidth )+"in}");}
      else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown property");}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " " << valignments_entrytable_string.back() << endl;}
    }

    //for compound name / ground state thermo properties (2 columns)
    uint n_cols_compoundname_thermoprops=2;
    vector<uint> valignments_compoundname_thermopropstable_uint;
    uint total_alignment_compoundname_thermoprops=0;
    for(uint i=0;i<n_cols_compoundname_thermoprops;i++){
      if(i==0){valignments_compoundname_thermopropstable_uint.push_back(3); total_alignment_compoundname_thermoprops+=valignments_compoundname_thermopropstable_uint.back();}
      if(i==1){valignments_compoundname_thermopropstable_uint.push_back(5); total_alignment_compoundname_thermoprops+=valignments_compoundname_thermopropstable_uint.back();}
    }

    vector<string> valignments_compoundname_thermopropstable_string;
    page_width_fraction=0.0;
    for(uint i=0,fl_size_i=valignments_compoundname_thermopropstable_uint.size();i<fl_size_i;i++){
      page_width_fraction=page_width*(double)valignments_compoundname_thermopropstable_uint[i]/(double)total_alignment_compoundname_thermoprops;
      if(i==0){valignments_compoundname_thermopropstable_string.push_back("L{"+aurostd::utype2string( page_width_fraction - penalty_tabcolsep - penalty_arrayrulewidth )+"in}");}
      if(i==1){valignments_compoundname_thermopropstable_string.push_back("R{"+aurostd::utype2string( page_width_fraction - penalty_tabcolsep - penalty_arrayrulewidth )+"in}");}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " " << valignments_compoundname_thermopropstable_string.back() << endl;}
    }

    //for 1 column tables (left align)
    uint n_cols_onecol_left=1;
    vector<uint> valignments_onecol_lefttable_uint;
    uint total_alignment_onecol_left=0;
    for(uint i=0;i<n_cols_onecol_left;i++){
      if(i==0){valignments_onecol_lefttable_uint.push_back(1); total_alignment_onecol_left+=valignments_onecol_lefttable_uint.back();}
    }

    vector<string> valignments_onecol_lefttable_string;
    page_width_fraction=0.0;
    for(uint i=0,fl_size_i=valignments_onecol_lefttable_uint.size();i<fl_size_i;i++){
      page_width_fraction=page_width*(double)valignments_onecol_lefttable_uint[i]/(double)total_alignment_onecol_left;
      if(i==0){valignments_onecol_lefttable_string.push_back("L{"+aurostd::utype2string( page_width_fraction - penalty_tabcolsep - penalty_arrayrulewidth )+"in}");}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " " << valignments_onecol_lefttable_string.back() << endl;}
    }

    //for equilibrium/decomposition reactions
    uint n_cols_reaction=2;
    vector<uint> valignments_reactiontable_uint;
    uint total_alignment_reaction=0;
    for(uint i=0;i<n_cols_reaction;i++){
      if(i==0){valignments_reactiontable_uint.push_back(1); total_alignment_reaction+=valignments_reactiontable_uint.back();}
      if(i==1){valignments_reactiontable_uint.push_back(3); total_alignment_reaction+=valignments_reactiontable_uint.back();}
    }

    vector<string> valignments_reactiontable_string;
    page_width_fraction=0.0;
    for(uint i=0,fl_size_i=valignments_reactiontable_uint.size();i<fl_size_i;i++){
      page_width_fraction=page_width*(double)valignments_reactiontable_uint[i]/(double)total_alignment_reaction;
      if(i==0){valignments_reactiontable_string.push_back("L{"+aurostd::utype2string( page_width_fraction - penalty_tabcolsep - penalty_arrayrulewidth )+"in}");}
      if(i==1){valignments_reactiontable_string.push_back("R{"+aurostd::utype2string( page_width_fraction - penalty_tabcolsep - penalty_arrayrulewidth )+"in}");}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " " << valignments_reactiontable_string.back() << endl;}
    }

    //ymin needs to change with iso_max plot
    bool plot_iso_max_latent_heat=m_cflags.flag("CHULL::PLOT_ISO_MAX_LATENT_HEAT");
    if(plot_iso_max_latent_heat){
      if(!m_formation_energy_hull){
        message << "CHULL::PLOT_ISO_MAX_LATENT_HEAT set to FALSE (limited to formation energy hulls only)";
        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_);
        plot_iso_max_latent_heat=false;
      }
      if(dimension!=2){
        message << "CHULL::PLOT_ISO_MAX_LATENT_HEAT set to FALSE (dimension!=2)";
        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_);
        plot_iso_max_latent_heat=false;
      }
    }
    double iso_max_latent_heat_min=AUROSTD_MAX_DOUBLE; //absolute min
    vector<uint> g_states=getGStates(false);
    if(plot_iso_max_latent_heat){
      double iso_max;
      for(uint i=0,fl_size_i=g_states.size();i<fl_size_i;i++){
        i_point=g_states[i];
        iso_max=isoMaxLatentHeat(m_points[i_point],0.5,_m_);
        if(iso_max<iso_max_latent_heat_min){iso_max_latent_heat_min=iso_max;}
      }
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " iso-max-latent-lines min=" << iso_max_latent_heat_min << endl;}
    }

    bool unary_g_state_unstable=false;  //the only way this is true is if we plot the unary g_state and it's unstable
    vector<uint> plot_points;

    if(!doc_only && (dimension == 2 || dimension == 3)) {
      ////////////////////////////////////////////////////////////////////////////
      // START Z and DIST filter
      ////////////////////////////////////////////////////////////////////////////

      // we need to handle filters first, so we can test to see if a colorbar is
      // possible
      if(filter_by_z || !plot_unstable) {
        z_filter_cutoff = 0.0;  //default
        if(filter_by_z) {
          z_filter_cutoff = filter_cutoff;
          //let's automatically override plot_unstable if z_filter_cutoff is set appropriately
          bool plot_unstable_old=plot_unstable;
          if(m_formation_energy_hull){if(z_filter_cutoff > 0.0){plot_unstable=true;}}
          else{if(z_filter_cutoff < 0.0){plot_unstable=true;}}
          if(plot_unstable_old!=plot_unstable){
            message << "CHULL::PLOT_UNSTABLE set to TRUE, z_filter_cutoff=" << z_filter_cutoff;
            pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_);
          }
          //[OBSOLETE CO20180227]if(!plot_unstable) {
          //[OBSOLETE CO20180227]  if(m_formation_energy_hull){if(z_filter_cutoff > 0.0){z_filter_cutoff = 0.0;}}
          //[OBSOLETE CO20180227]  else{if(z_filter_cutoff < 0.0){z_filter_cutoff = 0.0;}}
          //[OBSOLETE CO20180227]}
          //[OBSOLETE CO20180227]else{z_filter_cutoff = 0.0;}
        }
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " z_filter_cutoff=" << z_filter_cutoff << endl;}
      }
      if(filter_by_distance) {
        dist_filter_cutoff = filter_cutoff; 
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " DIST dist_filter_cutoff=" << dist_filter_cutoff << endl;}
        if(aurostd::zeroWithinTol(dist_filter_cutoff,ZERO_TOL) && plot_off_hull) {
          plot_off_hull=false;
          message << "CHULL::OFF_HULL set to FALSE, filter_by_distance=0.0";
          pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_);
        }
      }

      ////////////////////////////////////////////////////////////////////////////
      // END Z and DIST filter
      ////////////////////////////////////////////////////////////////////////////

      ////////////////////////////////////////////////////////////////////////////
      // Determine if we need to plot unstable (no stable points) (ROUND 1)
      ////////////////////////////////////////////////////////////////////////////

      plot_points_count = 0;
      plot_points_count_no_end_points = 0;
      min_point = max_point = 0.0;
      plot_points.clear();
      unary_g_state_unstable=false;
      //uint i_point;
      // no way to avoid this, we need to figure out if we can have a colorbar before
      // we go through the rest of the loops
      // so let's avoid double filtering, do it once and store to plot_points

      bool point_added;
      for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
        point_added=false;
        for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i&&!(dimension==3&&point_added);i++){  //we only add one point per coordgroup for 3D
          i_point=m_coord_groups[i_coord_group].m_points[i];
          if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
          const ChullPoint& point=m_points[i_point];
          if(!point.m_has_entry){continue;}
          if(!plot_unaries && point.isUnary()){continue;}
          if(!plot_off_hull && !point.m_is_g_state){continue;}  //if point.isUnary() and it gets here, only pass g-state
          if(!point.m_is_g_state){
            if(filter_by_z || !plot_unstable) {
              if(m_formation_energy_hull) {
                if(chull::H_f_atom(point, _m_) > z_filter_cutoff) {
                  //[CO20181226 - print later]if(filter_by_z) {
                  //[CO20181226 - print later]  message << "Excluding entry " << point.m_entry.auid;
                  //[CO20181226 - print later]  message << " with H_f = " << chull::H_f_atom(point, _m_);
                  //[CO20181226 - print later]  message << " (meV/atom) from plot";
                  //[CO20181226 - print later]  pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_,!filter_by_z);  // too much output to screen
                  //[CO20181226 - print later]}
                  continue;
                }
              }else{
                if(chull::T_S(point) < z_filter_cutoff) {
                  //[CO20181226 - print later]if(filter_by_z) {
                  //[CO20181226 - print later]  message << "Excluding entry " << point.m_entry.auid;
                  //[CO20181226 - print later]  message << " with T_S = " << chull::T_S(point);
                  //[CO20181226 - print later]  message << " (K) from plot";
                  //[CO20181226 - print later]  pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_,!filter_by_z);  // too much output to screen
                  //[CO20181226 - print later]}
                  continue;
                }
              }
            }
            if(filter_by_distance) {
              if(m_formation_energy_hull) {
                if(point.getDist2Hull(_m_) > dist_filter_cutoff) {
                  //[CO20181226 - print later]message << "Excluding entry " << point.m_entry.auid;
                  //[CO20181226 - print later]message << " with distance_hull_enthalpy_formation_atom = " << aurostd::utype2string(point.getDist2Hull(_m_),CHULL_PRECISION);
                  //[CO20181226 - print later]message << " (meV/atom) from plot";
                  //[CO20181226 - print later]pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_);  // too much output to screen
                  continue;
                }
              }else{
                if(point.getDist2Hull(_std_) < dist_filter_cutoff) {
                  //[CO20181226 - print later]message << "Excluding entry " << point.m_entry.auid;
                  //[CO20181226 - print later]message << " with distance_hull_entropic_temperature = " << aurostd::utype2string(point.getDist2Hull(_std_),CHULL_PRECISION);
                  //[CO20181226 - print later]message << " (K) from plot";
                  //[CO20181226 - print later]pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_);  // too much output to screen
                  continue;
                }
              }
            }
          }
          if(point.getLastCoord()<min_point){min_point=point.getLastCoord();}
          if(point.getLastCoord()>max_point){max_point=point.getLastCoord();}
          //[CO20181226 - set later]if(point.isUnary()){
          //[CO20181226 - set later]  if(m_formation_energy_hull) {if(chull::H_f_atom(point, _m_) > z_filter_cutoff) {unary_g_state_unstable=true;}}
          //[CO20181226 - set later]  else{if(chull::T_S(point) < z_filter_cutoff) {unary_g_state_unstable=true;}}
          //[CO20181226 - set later]}else{plot_points_count_no_end_points++;}
          plot_points_count++;
          //[CO20181226 - set later]plot_points.push_back(i_point);
          point_added=true;
        }
      }
      if(!plot_points_count){
        filter_by_z=true;filter_by_distance=false;
        plot_unstable=true;
        if(m_formation_energy_hull) {z_filter_cutoff=1000.0;} //1000 meV
        else{z_filter_cutoff=1.0/(((double)KBOLTZEV)*(0.5*log(0.5)+(1.0-0.5)*log(1.0-0.5)));}  //1.0 eV / (kB * 2 * 0.5 log(0.5))
        message << "CHULL::FILTER_SCHEME set to Z-axis (no stable points found)" << endl;
        message << "CHULL::FILTER_VALUE set to " << z_filter_cutoff << " (no stable found)" << endl;
        message << "CHULL::PLOT_UNSTABLE set to TRUE (no stable found)" << endl;
        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_WARNING_);
      }

      ////////////////////////////////////////////////////////////////////////////
      // START Gathering points for plotting (ROUND 2)
      ////////////////////////////////////////////////////////////////////////////

      plot_points_count = 0;
      plot_points_count_no_end_points = 0;
      min_point = max_point = 0.0;
      plot_points.clear();
      unary_g_state_unstable=false;
      //uint i_point;
      // no way to avoid this, we need to figure out if we can have a colorbar before
      // we go through the rest of the loops
      // so let's avoid double filtering, do it once and store to plot_points

      for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
        point_added=false;
        for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i&&!(dimension==3&&point_added);i++){  //we only add one point per coordgroup for 3D
          i_point=m_coord_groups[i_coord_group].m_points[i];
          if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
          const ChullPoint& point=m_points[i_point];
          if(!point.m_has_entry){continue;}
          if(!plot_unaries && point.isUnary()){continue;}
          if(!plot_off_hull && !point.m_is_g_state){continue;}  //if point.isUnary() and it gets here, only pass g-state
          if(!point.m_is_g_state){
            if(filter_by_z || !plot_unstable) {
              if(m_formation_energy_hull) {
                if(chull::H_f_atom(point, _m_) > z_filter_cutoff) {
                  if(filter_by_z) {
                    message << "Excluding entry " << point.m_entry.auid;
                    message << " with H_f = " << chull::H_f_atom(point, _m_);
                    message << " (meV/atom) from plot";
                    pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_,!filter_by_z);  // too much output to screen
                  }
                  continue;
                }
              }else{
                if(chull::T_S(point) < z_filter_cutoff) {
                  if(filter_by_z) {
                    message << "Excluding entry " << point.m_entry.auid;
                    message << " with T_S = " << chull::T_S(point);
                    message << " (K) from plot";
                    pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_,!filter_by_z);  // too much output to screen
                  }
                  continue;
                }
              }
            }
            if(filter_by_distance) {
              if(m_formation_energy_hull) {
                if(point.getDist2Hull(_m_) > dist_filter_cutoff) {
                  message << "Excluding entry " << point.m_entry.auid;
                  message << " with distance_hull_enthalpy_formation_atom = " << aurostd::utype2string(point.getDist2Hull(_m_),CHULL_PRECISION);
                  message << " (meV/atom) from plot";
                  pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_);  // too much output to screen
                  continue;
                }
              }else{
                if(point.getDist2Hull(_std_) < dist_filter_cutoff) {
                  message << "Excluding entry " << point.m_entry.auid;
                  message << " with distance_hull_entropic_temperature = " << aurostd::utype2string(point.getDist2Hull(_std_),CHULL_PRECISION);
                  message << " (K) from plot";
                  pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_);  // too much output to screen
                  continue;
                }
              }
            }
          }
          if(point.getLastCoord()<min_point){min_point=point.getLastCoord();}
          if(point.getLastCoord()>max_point){max_point=point.getLastCoord();}
          if(point.isUnary()){
            if(m_formation_energy_hull) {if(chull::H_f_atom(point, _m_) > z_filter_cutoff) {unary_g_state_unstable=true;}}
            else{if(chull::T_S(point) < z_filter_cutoff) {unary_g_state_unstable=true;}}
          }else{plot_points_count_no_end_points++;}
          plot_points_count++;
          plot_points.push_back(i_point);
          point_added=true;
        }
      }

      //set to doc_only if a worthless image
      //[CO20181226 - OBSOLETE]if(!doc_only && (dimension == 2 || dimension == 3))
      if(!plot_points_count)
      { //CO20200106 - patching for auto-indenting
        doc_only=true;
        message << "CHULL::DOC_ONLY set to TRUE (no plot points found)";
        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_);
      }

      ////////////////////////////////////////////////////////////////////////////
      // END Gathering points for plotting
      ////////////////////////////////////////////////////////////////////////////
    }

    //links options
    int links_setting=DEFAULT_CHULL_LATEX_LINKS;
    bool kill_all_links=(links_setting==0);       //no links whatsoever
    bool external_links=addExternalHyperlinks();      //weblinks
    bool internal_links_graph2report=(no_doc==false&&doc_only==false); //if doc_only==true, turn off links between graph and report
    bool internal_links_withinreport=(no_doc==false); //if no_doc==true, turn off links within report
    bool internal_links=addInternalHyperlinks(internal_links_graph2report,internal_links_withinreport);       //no jumping
    if(!internal_links){internal_links_graph2report=internal_links_withinreport=false;}

    //////////////////////////////////////////////////////////////////////////////
    // START Document header
    //////////////////////////////////////////////////////////////////////////////

    //doc_header_TEX_ss << "\\documentclass[12pt]{article}" << " \%HEADER" << endl; //non-standard
    doc_header_TEX_ss << "\\documentclass[10pt]{article}" << " \%HEADER" << endl;
    doc_header_TEX_ss << "\\usepackage[utf8x]{inputenc}" << " \%HEADER" << endl;
    doc_header_TEX_ss << "\\usepackage[table,dvipsnames]{xcolor}" << " \%HEADER" << endl;
    if(helvetica_font) {
      //doc_header_TEX_ss << "\\usepackage[scaled]{helvet}" << " \%HEADER" << endl;  //CO20180301
      doc_header_TEX_ss << "\\renewcommand\\familydefault{\\sfdefault}" << " \%HEADER" << endl;
      doc_header_TEX_ss << "\\usepackage{sansmath}" << " \%HEADER" << endl;
      doc_header_TEX_ss << "\\sansmath" << " \%HEADER" << endl; //enable sans-serif math for rest of document
      doc_header_TEX_ss << "\\usepackage{upgreek}" << " \%HEADER" << endl; //to fix \\Delta for H_f
    }
    doc_header_TEX_ss << "\\usepackage[T1]{fontenc}" << " \%HEADER" << endl;  // accents https://tex.stackexchange.com/questions/664/why-should-i-use-usepackaget1fontenc
    doc_header_TEX_ss << "\\usepackage{anyfontsize}" << " \%HEADER" << endl;  // arbitrary font sizes
    doc_header_TEX_ss << "\\usepackage{lmodern}" << " \%HEADER" << endl;  // high quality fonts
    doc_header_TEX_ss << "\\usepackage[yyyymmdd,hhmmss]{datetime}" << " \%HEADER" << endl;
    doc_header_TEX_ss << "\\renewcommand{\\dateseparator}{-}" << " \%HEADER" << endl;
    doc_header_TEX_ss << "\\usepackage{geometry} \%showframe" << " \%HEADER" << endl;
    if(!image_only){
      doc_header_TEX_ss << "\\usepackage{fancyhdr} \%column header" << " \%HEADER" << endl;
      doc_header_TEX_ss << "\\usepackage[absolute]{textpos} %showboxes" << " \%HEADER" << endl;
    }
    if(!no_doc && !image_only) {
      doc_header_TEX_ss << "\\usepackage{longtable}" << " \%HEADER" << endl;
      doc_header_TEX_ss << "\\vbadness=10000" << " \%HEADER" << endl; //patches vbadness warnings (hline issue) with longtable - https://tex.stackexchange.com/questions/71096/underfull-vbox-in-each-longtable-can-it-be-fixed-or-must-it-be-ignored
      //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]doc_header_TEX_ss << "\\usepackage{tabu}" << " \%HEADER" << endl;
      doc_header_TEX_ss << "\\usepackage{booktabs} \%midrule and toprule" << " \%HEADER" << endl;
      //http://tex.stackexchange.com/questions/167948/package-rerunfilecheck-warning-file-out-has-changed
      doc_header_TEX_ss << "\\usepackage{multirow} \%column header" << " \%HEADER" << endl;
    }

    if(1||!kill_all_links) { //need for paper reference (CHULL_JOURNAL_LATEX)
      doc_header_TEX_ss << "\\usepackage{hyperref} \\hypersetup{colorlinks=true,citecolor=blue,linkcolor=blue,urlcolor=blue}" << " \%HEADER" << endl;
      doc_header_TEX_ss << "\\usepackage{bookmark} \%hyperref without .out" << " \%HEADER" << endl;
    }

    doc_header_TEX_ss << "\\usepackage{array} \%m{}" << " \%HEADER" << endl;
    doc_header_TEX_ss << "\\usepackage{setspace} \%setstretch linespacing" << " \%HEADER" << endl;

    _doc_header_TEX_ss << "\\newcolumntype{L}[1]{>{\\raggedright\\arraybackslash}p{#1} }" << " \%HEADER" << endl;
    _doc_header_TEX_ss << "\\newcolumntype{C}[1]{>{\\centering  \\arraybackslash}p{#1} }" << " \%HEADER" << endl;
    _doc_header_TEX_ss << "\\newcolumntype{R}[1]{>{\\raggedleft \\arraybackslash}p{#1} }" << " \%HEADER" << endl;
    _doc_header_TEX_ss << "\\newcolumntype{X}[1]{>{\\setstretch{0.5} \\centering  \\arraybackslash}m{#1} } \%setstretch changes line spacing" << " \%HEADER" << endl;
    _doc_header_TEX_ss << "\\newcommand{\\newAlloy}[2]{%" << " \%HEADER" << endl;
    _doc_header_TEX_ss << "\\setAlloy{#1}" << " \%HEADER" << endl;
    _doc_header_TEX_ss << "\\setCountTotalEntries{#2}" << " \%HEADER" << endl;
    _doc_header_TEX_ss << "}" << " \%HEADER" << endl;

    if(!image_only){
      //to set count variable, we need to create a "section" where it is defined
      //alloys are de facto subsections, section usually "binaries" or "Co-type", etc.
      //https://tex.stackexchange.com/questions/186582/defining-a-dynamic-variable
      //_doc_header_TEX_ss << "\\newcommand{\\newAlloyFormatted}[3]{%" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\newcommand{\\newAlloyFormatted}[2]{%" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\clearpage" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\phantomsection" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\addcontentsline{toc}{subsection}{#1}" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\newAlloy{#1}{#2}" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "}" << " \%HEADER" << endl;
      //_doc_header_TEX_ss << "\\setAlloy{#1}" << " \%HEADER" << endl;
      //_doc_header_TEX_ss << "\\setCountTotalEntries{#2}" << " \%HEADER" << endl;
      //_doc_header_TEX_ss << "\\setCountHullEntries{#2}" << " \%HEADER" << endl;
      //_doc_header_TEX_ss << "\\setCountTotalEntries{#3}" << " \%HEADER" << endl;
    }

    //alloy needs to be dynamic as it sets in the header/footer!
    _doc_header_TEX_ss << "\\newcommand{\\alloy}{}"<< " \%HEADER" << endl;
    _doc_header_TEX_ss << "\\newcommand{\\setAlloy}[1]{\\renewcommand{\\alloy}{#1}}" << " \%HEADER" << endl;
    _doc_header_TEX_ss << "\\newcommand{\\printAlloy}{\\alloy}" << " \%HEADER" << endl;
    //count variable needs to be dynamic as it sets in the header/footer!
    //_doc_header_TEX_ss << "\\newcommand{\\countHullEntries}{}"<< " \%HEADER" << endl;
    //_doc_header_TEX_ss << "\\newcommand{\\setCountHullEntries}[1]{\\renewcommand{\\countHullEntries}{#1}}" << " \%HEADER" << endl;
    //_doc_header_TEX_ss << "\\newcommand{\\printCountHullEntries}{\\countHullEntries}" << " \%HEADER" << endl;
    //total count variable needs to be dynamic as it sets in the header/footer!
    _doc_header_TEX_ss << "\\newcommand{\\countTotalEntries}{}"<< " \%HEADER" << endl;
    _doc_header_TEX_ss << "\\newcommand{\\setCountTotalEntries}[1]{\\renewcommand{\\countTotalEntries}{#1}}" << " \%HEADER" << endl;
    _doc_header_TEX_ss << "\\newcommand{\\printCountTotalEntries}{\\countTotalEntries}" << " \%HEADER" << endl;

    // things to add later, depending on we include hull image on first page
    _doc_header_TEX_ss << "\\begin{document}" << " \%HEADER" << endl;

    // BEAUTIFUL SOLUTION to use links with labels
    if(!plot_unaries && internal_links_graph2report) {
      _doc_header_TEX_ss << "\\newcommand{\\hyperrefTitle}[2]{\\hyperref[#1]{#2}} \%we need this because axis doesn't like []" << " \%HEADER" << endl;
    }
    if(!doc_only && (dimension == 2 || dimension == 3)) {
      _doc_header_TEX_ss << "\\definecolor{pranab_green}{rgb}{0.31,0.53,0.10}" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\definecolor{pranab_red}{rgb}{0.85,0.23,0.11}" << " \%HEADER" << endl;
    }

    //moving from bottom so header stuff stays together
    if(!image_only){
      // fancypagestyle1 for convex hull illustration
      _doc_header_TEX_ss << "\\setlength{\\TPHorizModule}{\\paperwidth}\\setlength{\\TPVertModule}{\\paperheight}" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\fancypagestyle{hullImage}{" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\fancyhf{}" << " \%HEADER" << endl;
      if(!(small_banner || no_banner)){
        _doc_header_TEX_ss << "\\fancyhead[C]{" << " \%HEADER" << endl;
        if(print_aflow_logo_full||print_aflow_text_logo){
          _doc_header_TEX_ss << "\\begin{textblock}{0.1}[0.5,1.0](0.09,0.95){\\rotatebox{90}{";
          if(!print_aflow_logo_full){_doc_header_TEX_ss << "\\shortstack[l]{";}
          if(external_links) {_doc_header_TEX_ss << "\\href{http://aflow.org}";}
          if(print_aflow_logo_full) {_doc_header_TEX_ss << "{\\includegraphics[scale=0.2]{" << aflow_logo_full_file << "}}";}
          else{
            _doc_header_TEX_ss << "{\\LARGE AFLOW V" << string(AFLOW_VERSION) << "}";
            _doc_header_TEX_ss << "\\\\";
            _doc_header_TEX_ss << "{\\large C. Oses and S. Curtarolo}";
            _doc_header_TEX_ss << "}";
          }
          _doc_header_TEX_ss << "}}\\end{textblock}" << " \%HEADER" << endl;
        }

        if(print_logo_2||print_aflow_webaddress_logo){
          _doc_header_TEX_ss << "\\begin{textblock}{0.1}[0.5,0.0](0.09,0.05){\\rotatebox{90}{";
          if(print_logo_2){
            if(external_links) {_doc_header_TEX_ss << "\\href{http://www.nomad-coe.eu/}";}
            _doc_header_TEX_ss << "{\\includegraphics[scale=0.1]{" << logo_file_2 << "}}";
          }else{  //print_aflow_webaddress_logo
            if(external_links) {_doc_header_TEX_ss << "\\href{http://aflow.org}";}
            _doc_header_TEX_ss << "{" << label_image_font_size << "{\\fontfamily{\\sfdefault}\\selectfont{aflow.org}}}";  //\\fontfamily{\\phv}\\selectfont
          }
          _doc_header_TEX_ss << "}}\\end{textblock}" << " \%HEADER" << endl;
        }
        _doc_header_TEX_ss << "}" << " \%HEADER" << endl;

        _doc_header_TEX_ss << "\\fancyfoot[C]{" << " \%HEADER" << endl;
        _doc_header_TEX_ss << "\\begin{textblock}{0.05}[0.5,1.0](0.925,0.95){\\rotatebox{90}{" << general_image_font_size << " " << string("\\printCountTotalEntries\\") << " entries}}\\end{textblock}" << " \%HEADER" << endl; //string("\\printCountHullEntries") //hull count
        _doc_header_TEX_ss << "\\begin{textblock}{0.05}[0.5,0.0](0.925,0.05){\\rotatebox{90}{" << general_image_font_size << " \\today}}\\end{textblock}" << " \%HEADER" << endl; //\\today~\\currenttime
        _doc_header_TEX_ss << "}" << " \%HEADER" << endl;
      }
      _doc_header_TEX_ss << "\\renewcommand{\\headrulewidth}{0pt}     \% size of header line" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\renewcommand{\\footrulewidth}{0pt}     \% size of header line" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "}" << " \%HEADER" << endl;
    }
    if(!no_doc && !image_only) {
      _doc_header_TEX_ss << "\\newsavebox{\\snapshotTableHeaderOne}" << " \%this IS dynamic (countTotalEntries), so define command to update it (must be in right environment)" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\newcommand{\\updateSnapshotTableHeaderOne}{" << " \%HEADER" << endl;
      //LRBOXES store as text, so whatever inherits it doesn't know it's a XX (tabu in this case, which conflicts with fancyheader)
      _doc_header_TEX_ss << "\\begin{lrbox}{\\snapshotTableHeaderOne}" << " \%store as text \%HEADER" << endl;
      //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]_doc_header_TEX_ss << "\\begin{tabu}{X[1,l,m]" << (print_logo_2?string("X[1,c,m]"):string("")) << "X[1,r,m]}" << " \%HEADER" << endl;
      double tableheaderone_colsize=(LATEX_WIDTH_LETTER_STD-(LATEX_LEFT_MARGIN_LETTER_STD+LATEX_RIGHT_MARGIN_LETTER_STD))/(print_logo_2?3.0:2.0);
      _doc_header_TEX_ss << "\\noindent\\begin{tabular}{";
      _doc_header_TEX_ss << "@{}L{" << tableheaderone_colsize << "in}";
      if(print_logo_2){_doc_header_TEX_ss << "@{}C{" << tableheaderone_colsize << "in}";}
      _doc_header_TEX_ss << "@{}R{" << tableheaderone_colsize << "in}@{}}" << " \%HEADER" << endl;
      if(!print_aflow_logo_full){_doc_header_TEX_ss << "\\shortstack[l]{";}
      if(external_links) {_doc_header_TEX_ss << "\\href{" << AFLOW_WEB << "}";}
      if(print_aflow_logo_full) {_doc_header_TEX_ss << "{\\raisebox{-.5\\height}{\\includegraphics[scale=0.25]{" << aflow_logo_full_file << "}}}";}
      else{
        _doc_header_TEX_ss << "{\\LARGE AFLOW V" << string(AFLOW_VERSION) << "}";
        _doc_header_TEX_ss << "\\\\";
        _doc_header_TEX_ss << "{\\large C. Oses and S. Curtarolo}";
        _doc_header_TEX_ss << "}";
      }
      _doc_header_TEX_ss << " & ";
      if(print_logo_2){
        if(external_links) {_doc_header_TEX_ss << "\\href{" << NOMAD_WEB << "}";}
        _doc_header_TEX_ss << "{\\raisebox{-.5\\height}{\\includegraphics[scale=0.125]{" << logo_file_2 << "}}}";
        _doc_header_TEX_ss << " & ";
      }
      //include small redundant header - START
      _doc_header_TEX_ss << "\\begin{tabular}{@{}r@{}}";
      _doc_header_TEX_ss << "\\large " << string("\\printAlloy") << "\\ summary"; // (V" << string(AFLOW_VERSION) << ")"; //Materials Snapshot";
      _doc_header_TEX_ss << "\\\\";
      _doc_header_TEX_ss << "\\large " << string("\\printCountTotalEntries\\") << " entries";
      _doc_header_TEX_ss << "\\\\";
      _doc_header_TEX_ss << "\\large \\today";  //\\today~\\currenttime
      _doc_header_TEX_ss << "\\end{tabular}" << " \%HEADER" << endl;
      //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]_doc_header_TEX_ss << "\\end{tabu}" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\end{tabular}" << " \%HEADER" << endl;
      //include small redundant header - END
      _doc_header_TEX_ss << "\\end{lrbox}" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "}" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\newsavebox{\\snapshotTableHeaderTwo}" << " \%this is NOT dynamic, so simply define once with right geometry and keep forever" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\newgeometry{" << MARGIN_REPORT << "}" << " \%HEADER" << endl;  //force geometry of the table early, so columns are properly sized in headers //new geometry introduces \\newpage automatically
      //LRBOXES store as text, so whatever inherits it doesn't know it's a XX (tabu in this case, which conflicts with fancyheader)
      _doc_header_TEX_ss << "\\begin{lrbox}{\\snapshotTableHeaderTwo}" << " \%store as text \%HEADER" << endl;
      _doc_header_TEX_ss << getSnapshotTableHeader(headers,true);
      _doc_header_TEX_ss << "\\end{lrbox}" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\restoregeometry" << " \%HEADER" << endl;
      // fancypagestyle2 for first page of report
      _doc_header_TEX_ss << "\\fancypagestyle{reportPage1}{" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\fancyhf{}" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\fancyhead[C]{\\usebox{\\snapshotTableHeaderOne}}" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\fancyfoot[C]{\\thepage}" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\renewcommand{\\headrulewidth}{0pt}     \% size of header line" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\renewcommand{\\footrulewidth}{0pt}     \% size of header line" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "}" << " \%HEADER" << endl;
      // fancypagestyle3 for beyond first page
      _doc_header_TEX_ss << "\\fancypagestyle{reportPage2}{" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\fancyhf{}" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\fancyhead[C]{\\usebox{\\snapshotTableHeaderTwo}}" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\fancyfoot[C]{\\thepage}" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\renewcommand{\\headrulewidth}{0pt}     \% size of header line" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "\\renewcommand{\\footrulewidth}{0pt}     \% size of header line" << " \%HEADER" << endl;
      _doc_header_TEX_ss << "}" << " \%HEADER" << endl;
    }
    if(!image_only){_doc_header_TEX_ss << "\\newAlloyFormatted{" << input_hyphened << "}{" << count_total_entries << "}" << endl;} //{" << count_total_entries << "}" << endl;

    //////////////////////////////////////////////////////////////////////////////
    // END Document header (for now)
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // START Image on first page
    //////////////////////////////////////////////////////////////////////////////

    if(!doc_only && (dimension == 2 || dimension == 3)) {
      ////////////////////////////////////////////////////////////////////////////
      // START More header commands specific to tikz pic
      ////////////////////////////////////////////////////////////////////////////

      doc_header_TEX_ss << "\\usepackage{pgfplots}" << " \%HEADER" << endl;
      doc_header_TEX_ss << "\\usepackage{pdflscape}" << " \%HEADER" << endl;
      //[MOVED UP]doc_header_TEX_ss << "\\usepackage{array} \%m" << " \%HEADER" << endl;  %brackets
      //[MOVED UP]doc_header_TEX_ss << "\\usepackage{setspace} \%setstretch linespacing" << " \%HEADER" << endl;
      doc_header_TEX_ss << "\\pgfplotsset{compat=1.10}" << " \%HEADER" << endl;
      doc_header_TEX_ss << "\\usepgfplotslibrary{ternary,units}" << " \%HEADER" << endl;
      doc_header_TEX_ss << "\\usetikzlibrary{decorations.pathmorphing,pgfplots.units,backgrounds}" << " \%HEADER" << endl;
      doc_header_TEX_ss << "\\usepackage{tikz}" << " \%HEADER" << endl;
      doc_header_TEX_ss << "\\usetikzlibrary{positioning}" << " \%HEADER" << endl;

      if(dimension == 3) {
        doc_header_TEX_ss << "\\usetikzlibrary{pgfplots.ternary}" << " \%HEADER" << endl;
      }
      if(image_only) {
        doc_header_TEX_ss << "\\pgfrealjobname{CHull}" << " \%HEADER" << endl;  // dummy name
      }
      doc_header_TEX_ss << "\\pgfdeclarelayer{background}" << " \%HEADER" << endl;
      doc_header_TEX_ss << "\\pgfdeclarelayer{foreground}" << " \%HEADER" << endl;
      doc_header_TEX_ss << "\\pgfsetlayers{background,main,foreground}" << " \%HEADER" << endl;
      doc_header_TEX_ss << "\\AtBeginEnvironment{tikzpicture}{\\tracinglostchars=0\\relax}" << " \%HEADER" << endl; //CO20221111 - avoids printing dumb warnings

      //[MOVED UP]doc_header_TEX_ss << "\\newcolumntype{L}[1]{>{\\raggedright\\arraybackslash}p{#1} }" << " \%HEADER" << endl;
      //[MOVED UP]doc_header_TEX_ss << "\\newcolumntype{C}[1]{>{\\centering  \\arraybackslash}p{#1} }" << " \%HEADER" << endl;
      //[MOVED UP]doc_header_TEX_ss << "\\newcolumntype{R}[1]{>{\\raggedleft \\arraybackslash}p{#1} }" << " \%HEADER" << endl;
      //[MOVED UP]doc_header_TEX_ss << "\\newcolumntype{X}[1]{>{\\setstretch{0.5} \\centering  \\arraybackslash}m{#1} } \%setstretch changes line spacing" << " \%HEADER" << endl;

      doc_header_TEX_ss << _doc_header_TEX_ss.str();
      _doc_header_TEX_ss.str("");  // don't repeat

      doc_header_TEX_ss << "\\newgeometry{" << MARGIN_PICTURE << "}" << endl; //new geometry introduces \\newpage automatically
      if(image_only) {doc_header_TEX_ss << "\\beginpgfgraphicnamed{aflow_" << input << "_hull}" << " \%HEADER" << endl;}
      if(!image_only){
        doc_header_TEX_ss << "\\thispagestyle{hullImage}" << endl;  //empty
        doc_header_TEX_ss << "\\begin{landscape}" << endl;
        doc_header_TEX_ss << "\\centering{" << endl;
        doc_header_TEX_ss << "\\topskip0pt" << endl;
        doc_header_TEX_ss << "\\vspace*{\\fill}" << endl;
      }

      ////////////////////////////////////////////////////////////////////////////
      // END More header commands specific to tikz pic
      ////////////////////////////////////////////////////////////////////////////

      ////////////////////////////////////////////////////////////////////////////
      // START Determining if points range enough to have color/heat map
      // functionality
      ////////////////////////////////////////////////////////////////////////////

      int approx_num_ticks=5;
      double extra_padding_multiplier=(rotate_labels?0.225:0.1);

      if(m_formation_energy_hull) {
        min_point = chull::convertUnits(min_point, _m_);
        if(plot_iso_max_latent_heat){ //min_point might change with iso_max
          double old_min_point=min_point;
          min_point=min(min_point,iso_max_latent_heat_min);
          if(abs(old_min_point-min_point)>abs(max_point - old_min_point)/(double)approx_num_ticks){ //is it bigger than ~tic distance?
            extra_padding_multiplier/=2.0; //don't need so much MORE padding, T_S label is flat, and we just made sure T_S is well below rotated labels
          }
        }
        max_point = chull::convertUnits(max_point, _m_);
      }

      //do after everything
      if(image_only){extra_padding_multiplier*=1.5;}

      point_range = max_point - min_point;
      //bool set_max_tick_spacing_2D=false;//true; //max space between ticks=XX (in yticklabel style),  //CO20180227  //this hardly works, we need to set it exactly, COMPLETELY REMOVED
      //int max_tick_spacing=1;
      bool set_y_range_exactly=true;
      double y_tick_distance=0.25;
      double y_tick_distance_right=y_tick_distance; ///10.0;
      double ymin=(m_formation_energy_hull?-ZERO_RANGE_TOL:0);//-ZERO_RANGE_TOL;
      double ymax=(m_formation_energy_hull?0:ZERO_RANGE_TOL);//ZERO_RANGE_TOL;
      double y_range=abs(ymax-ymin);
      double round_to_value=getRoundToValue(point_range);
      double extra_padding=(point_range)*extra_padding_multiplier; //to avoid label clashing with axis, 20% should be enough
      //int extra_padding=aurostd::roundDouble((point_range)*0.2, (int)round_to_value, false);  //round down to keep tight, unless you get 0, need more padding  //CO20190724 - explicit double->int conversion (floor is fine)
      //if(extra_padding==0){extra_padding=aurostd::roundDouble((point_range)*0.2, (int)round_to_value, true);}  //CO20190724 - explicit double->int conversion (floor is fine)

      //adding conversion between meV/atom to kJ/mol
      bool showkJmolaxis = m_formation_energy_hull;  //default
      double kJmol_threshold = 1.0/((double)meVatom2kJmol); //should give about 10, 1 kJ/mol -> 10 meV/atom
      if(0 && showkJmolaxis && y_range<kJmol_threshold){showkJmolaxis=false;}  //we need integer y tick distance, don't show if it's not at least as big as ~10meV/atom, we should now be able to handle fractional y_dist

      if(showkJmolaxis){y_tick_distance_right/=10.0;} //default value, store AFTER setting showkJmolaxis

      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " range" << endl;
        cerr << __AFLOW_FUNC__ << " max:     " << max_point << endl;
        cerr << __AFLOW_FUNC__ << " min:     " << min_point << endl;
        cerr << __AFLOW_FUNC__ << " point_range: " << point_range << endl;
      }

      // range functionality
      // http://tex.stackexchange.com/questions/69248/set-ticklabels-for-colorbar

      if(dimension == 2) {
        include_color_bar=false;
        message << "CHULL::COLOR_BAR set to FALSE (dimension==2)";
        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_);

        show_heat_map=false;
        message << "CHULL::HEAT_MAP set to FALSE (dimension==2)";
        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_);
      }
      // determine whether we can have a colorbar
      // recently added: (plot_points_count==1&&!show_heat_map)
      // we need difference data, so if there's only one "colored" data, no way to establish a color spectrum
      if( !plot_points_count || (dimension==3 && !show_heat_map && plot_points_count<2) || (abs(point_range) < ZERO_TOL)) {
        include_color_bar=false;
        show_heat_map=false;
        display_color_gradient=false;
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " plot points count = " << plot_points_count << endl;
          cerr << __AFLOW_FUNC__ << " plot points count no end points = " << plot_points_count_no_end_points << endl;
        }
        if(!plot_points_count) {
          message << "CHULL::COLOR_BAR set to FALSE, no entries found";
          pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_);
          message << "CHULL::COLOR_GRADIENT set to FALSE, no entries found";
          pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_);
          message << "CHULL::HEAT_MAP set to FALSE, no entries found";
          pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_);
        }
        if(dimension==3 && !show_heat_map && plot_points_count<2) {
          message << "CHULL::COLOR_BAR set to FALSE, not enough non-unary entries found";
          pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_);
          message << "CHULL::COLOR_GRADIENT set to FALSE, not enough non-unary entries found";
          pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_);
          message << "CHULL::HEAT_MAP set to FALSE, not enough non-unary entries found";
          pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_);
        }
        if(aurostd::zeroWithinTol(point_range,ZERO_TOL)) {
          message << "CHULL::COLOR_BAR set to FALSE, hull has no depth";
          pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_);
          message << "CHULL::COLOR_GRADIENT set to FALSE, hull has no depth";
          pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_);
          message << "CHULL::HEAT_MAP set to FALSE, hull has no depth";
          pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_OPTION_);
        }
      }

      ////////////////////////////////////////////////////////////////////////////
      // END Determining if points range enough to have color/heat map
      // functionality
      ////////////////////////////////////////////////////////////////////////////

      // doc_header_TEX_ss ends here, we need to get data first, see
      // pseudo_preliminary_axes_TEX_ss (hrefs)

      ////////////////////////////////////////////////////////////////////////////
      // START Tikz pic settings
      ////////////////////////////////////////////////////////////////////////////

      //_tikzpic_settings_TEX_ss comes before tikzpic_TEX_ss
      _tikzpic_settings_TEX_ss << "\\begin{tikzpicture}[font=" << general_image_font_size;
      _tikzpic_settings_TEX_ss << "]";
      _tikzpic_settings_TEX_ss << endl;
      // end of _tikzpic_settings_TEX_ss, we need to get data first, see
      // pseudo_preliminary_axes_TEX_ss (hrefs)
      // common settings between two axes
      if((!plot_points_count) || (abs(point_range) < ZERO_RANGE_TOL)) {
        //simply do not calculate ymin, ymax (ints)
        //
        //ymin=-ZERO_RANGE_TOL;
        //ymax=ZERO_RANGE_TOL;
      }else{
        // easier to figure out ranges in meV, stays meV
        // between lowest label and bottom line, increase if need more
        // we might have hull endpoints higher than 0
        if(!m_formation_energy_hull && !(plot_unaries && unary_g_state_unstable) && (!plot_unstable || abs(min_point) < ZERO_TOL)) {ymin=0;}
        else{
          ymin=min_point;
          if(dimension==2){ //only do rounding for dim==2, dim==3 we need EXACT values or colorbar gets screwed up
            if(m_formation_energy_hull){ymin-=extra_padding;}  //CO20180227 - we need to avoid labels here
            ymin=aurostd::roundDouble(ymin,(int)round_to_value,false);  //we want min to be MORE NEGATIVE, so round DOWN //CO20190724 - explicit double->int conversion (floor is fine)
          }
        }
        // between highest label and top line, increase if need more
        // we might have hull endpoints higher than 0
        if(m_formation_energy_hull && !(plot_unaries && unary_g_state_unstable) && (!plot_unstable || abs(max_point) < ZERO_TOL)) {ymax=0;}
        else{
          ymax=max_point;
          if(dimension==2){ //only do rounding for dim==2, dim==3 we need EXACT values or colorbar gets screwed up
            if(!m_formation_energy_hull){ymax+=extra_padding;}  //CO20180227 - we need to avoid labels here
            ymax=aurostd::roundDouble(ymax,(int)round_to_value,true);  //we want max to be MORE POSITIVE, so round UP  //CO20190724 - explicit double->int conversion (floor is fine)
          }
        }
        y_range=abs(ymax-ymin);
        //max_tick_spacing=int(y_range/4);
        y_tick_distance=y_tick_distance_right=getYTickDistance(y_range,approx_num_ticks,round_to_value);

        //need new y_tick_distance if showkJmolaxis
        if(showkJmolaxis){
          double round_to_value_right=getRoundToValue(meVatom2kJmol*point_range);
          y_tick_distance_right=getYTickDistance(meVatom2kJmol*y_range,approx_num_ticks,round_to_value_right);
        }

        if(0||LDEBUG){
          cerr << __AFLOW_FUNC__ << " point_range                                        = " << point_range << endl;
          cerr << __AFLOW_FUNC__ << " min_point                                          = " << min_point << endl;
          cerr << __AFLOW_FUNC__ << " max_point                                          = " << max_point << endl;
          cerr << __AFLOW_FUNC__ << " round_to_value                                     = " << round_to_value << endl;
          //cerr << __AFLOW_FUNC__ << " aurostd::roundDouble(min_point,(int)round_to_value,false)        = " << aurostd::roundDouble(min_point, (int)round_to_value, false) << endl; //CO20190724 - explicit double->int conversion (floor is fine)
          //cerr << __AFLOW_FUNC__ << " aurostd::roundDouble(max_point,(int)round_to_value,true)         = " << aurostd::roundDouble(max_point, (int)round_to_value, true) << endl;  //CO20190724 - explicit double->int conversion (floor is fine)
          cerr << __AFLOW_FUNC__ << " extra_padding                                      = " << extra_padding << endl;
          cerr << __AFLOW_FUNC__ << " ymin                                               = " << ymin << endl;
          cerr << __AFLOW_FUNC__ << " ymax                                               = " << ymax << endl;
          //cerr << __AFLOW_FUNC__ << " max_tick_spacing                                   = " << max_tick_spacing << endl;
          cerr << __AFLOW_FUNC__ << " y_tick_distance                                    = " << y_tick_distance << endl;
          cerr << __AFLOW_FUNC__ << " y_tick_distance_right                              = " << y_tick_distance_right << endl;
        }
        //need to convert to kJ/mol for right axis
        //if(dimension==2){
        //  common_settings_TEX_ss << "ymin=" << ymin << "," << endl;
        //  common_settings_TEX_ss << "ymax=" << ymax << "," << endl;
        //}
      }
      common_settings_TEX_ss << "clip=false," << endl;  // for labels

      double image_width=19;
      if(dimension==3){image_width=17;}
      common_settings_TEX_ss << "width=" << image_width << "cm," << endl;

      //hsize is not well defined, especially for image only
      //double image_height_proportion=0.65;
      //if(dimension==3){image_height_proportion=1.0;}
      //string image_height_str=aurostd::utype2string(image_height_proportion)+"\\hsize";  //better to define in terms of page, since this is what we are confined to
      //better to use cm
      double image_width_proportion=0.8;                //approx 15cm given 19cm width
      if(dimension==3){image_width_proportion=1.2;}     //approx 20cm given 17cm width
      common_settings_TEX_ss << "height=" << image_width_proportion*image_width << "cm," << endl;

      if(dimension==3){common_settings_TEX_ss << "axis line style={line width=3pt}," << endl;}
      //if(dimension == 2) {common_settings_TEX_ss << 20;}
      //else{  // dimension==3
      //  //CO change me if we see otherwise
      //  //if(include_color_bar) {
      //  //  if(image_only) {common_settings_TEX_ss << 21;}
      //  //  else{common_settings_TEX_ss << 17;}  // only time we change size
      //  //}else{common_settings_TEX_ss << 21;}
      //}
      //common_settings_TEX_ss << "cm," << endl;
      //if(dimension == 2) {
      //  if(image_only) {common_settings_TEX_ss << "height=0.8\\hsize," << endl;}
      //  else{common_settings_TEX_ss << "height=0.55\\hsize," << endl;}
      //}else{
      //  common_settings_TEX_ss << "height=\\hsize," << endl;
      //  common_settings_TEX_ss << "axis line style={line width=3pt}," << endl;
      //}
      common_settings_TEX_ss << "grid=none," << endl;
      common_settings_TEX_ss << "axis on top," << endl;  // prevents heat map overlap issues

      ////////////////////////////////////////////////////////////////////////////
      // START Axis
      ////////////////////////////////////////////////////////////////////////////

      double unary_distance_axis_3D=0.25;
      if(image_only){unary_distance_axis_3D=0.35;} //CO FIX

      uint g_state;
      if(dimension == 2) {
        uint top_axis_point, bottom_axis_point;
        if(!reverse_axes) {
          top_axis_point = 0;
          bottom_axis_point = 1;
        }else{
          top_axis_point = 1;
          bottom_axis_point = 0;
        }
        // first axis - kJ/atom, top and right axes, x label only for padding
        if(showkJmolaxis) //||!image_only)   //padding purposes
        { //CO20200106 - patching for auto-indenting
          //only needed if we want to pad for pseudo-centering (by eye), it is perfectly centered, but logos are different shapes
          //no need to pad for eye, waste of time
          if(0){ 
            pseudo_preliminary_axes_TEX_ss << "\\begin{axis}[" << endl;
            //BETTER CENTERING without labels
            //if(!image_only){
            //  pseudo_preliminary_axes_TEX_ss << "xlabel={";
            //  g_state=getUnaryGState(top_axis_point);
            //  //opacity==0
            //  //if(!plot_unaries && isViableGState(g_state)) {
            //  //  if(internal_links_graph2report) {
            //  //    pseudo_preliminary_axes_TEX_ss << "\\hyperrefTitle{" << input << "_" << m_points[g_state].m_entry.auid << "}";
            //  //  }else if(no_doc && external_links) {
            //  //    pseudo_preliminary_axes_TEX_ss << "\\href{" << ENTRY_PAGE_URL_PREFIX << m_points[g_state].m_entry.auid << "}";
            //  //  }
            //  //}
            //  pseudo_preliminary_axes_TEX_ss << "{" << m_velements[top_axis_point] << "}";
            //  pseudo_preliminary_axes_TEX_ss << "}," << endl;
            //  pseudo_preliminary_axes_TEX_ss << "xlabel style={font=" << label_image_font_size << ",at={(xticklabel cs:" << top_axis_point << ")}";
            //  pseudo_preliminary_axes_TEX_ss << ",opacity=0}, \%we need labels for proper padding, but avoid bolder labels with opacity==0" << endl;
            //}
            pseudo_preliminary_axes_TEX_ss << "ylabel={{(kJ/mol)}}," << endl;
            pseudo_preliminary_axes_TEX_ss << "ylabel style={font=" << label_image_font_size; // << ",rotate=180";
            pseudo_preliminary_axes_TEX_ss << ",opacity=0}, \%we need labels for proper padding, but avoid bolder labels with opacity==0" << endl;
            pseudo_preliminary_axes_TEX_ss << "axis x line*=top," << endl;
            pseudo_preliminary_axes_TEX_ss << "axis y line*=right," << endl;
            pseudo_preliminary_axes_TEX_ss << "tick pos=right," << endl;
            pseudo_preliminary_axes_TEX_ss << "xtick={1,0.8,0.6,0.4,0.2,0}," << endl;
            pseudo_preliminary_axes_TEX_ss << "xticklabels={}," << endl;
            pseudo_preliminary_axes_TEX_ss << "xticklabel shift=6pt," << endl;
            pseudo_preliminary_axes_TEX_ss << "yticklabel shift=6pt," << endl;
            pseudo_preliminary_axes_TEX_ss << "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5";
            pseudo_preliminary_axes_TEX_ss << ",opacity=0}, \%we need labels for proper padding, but avoid bolder labels with opacity==0" << endl;
            if(set_y_range_exactly){pseudo_preliminary_axes_TEX_ss << "ytick distance={" << y_tick_distance_right << "}," << endl;}
            pseudo_preliminary_axes_TEX_ss << "scaled y ticks=false," << endl;
            pseudo_preliminary_axes_TEX_ss << "xmin=0," << endl;
            pseudo_preliminary_axes_TEX_ss << "xmax=1," << endl;
            //removed ymin/ymax from common settings
            pseudo_preliminary_axes_TEX_ss << "ymin=" << (showkJmolaxis?meVatom2kJmol:1.0) * ymin << "," << endl;
            pseudo_preliminary_axes_TEX_ss << "ymax=" << (showkJmolaxis?meVatom2kJmol:1.0) * ymax << "," << endl;
            // insert common settings
            pseudo_preliminary_axes_TEX_ss << common_settings_TEX_ss.str();
            pseudo_preliminary_axes_TEX_ss << "]" << endl;
            pseudo_preliminary_axes_TEX_ss << "\\end{axis}" << endl;
          }

          // second axis - kJ/atom, top and right axes, x label only for padding
          pseudo_preliminary_axes_TEX_ss << "\\begin{axis}[" << endl;
          //BETTER CENTERING without labels
          //if(!image_only){
          //  pseudo_preliminary_axes_TEX_ss << "xlabel={";
          //  g_state=getUnaryGState(bottom_axis_point);
          //  //opacity==0
          //  //if(!plot_unaries && isViableGState(g_state)) {
          //  //  if(internal_links_graph2report) {
          //  //    pseudo_preliminary_axes_TEX_ss << "\\hyperrefTitle{" << input << "_" << m_points[g_state].m_entry.auid << "}";
          //  //  }else if(no_doc && external_links) {
          //  //    pseudo_preliminary_axes_TEX_ss << "\\href{" << ENTRY_PAGE_URL_PREFIX << m_points[g_state].m_entry.auid << "}";
          //  //  }
          //  //}
          //  pseudo_preliminary_axes_TEX_ss << "{" << m_velements[bottom_axis_point] << "}";
          //  pseudo_preliminary_axes_TEX_ss << "}," << endl;
          //  pseudo_preliminary_axes_TEX_ss << "xlabel style={font=" << label_image_font_size << ",at={(xticklabel cs:" << bottom_axis_point << ")}";
          //  pseudo_preliminary_axes_TEX_ss << ",opacity=0}, \%we need labels for proper padding, but avoid bolder labels with opacity==0" << endl;
          //}
          pseudo_preliminary_axes_TEX_ss << "ylabel={{(kJ/mol)}}," << endl;
          pseudo_preliminary_axes_TEX_ss << "ylabel style={font=" << label_image_font_size; // << ",rotate=180";
          pseudo_preliminary_axes_TEX_ss << "}," << endl;
          pseudo_preliminary_axes_TEX_ss << "axis x line*=top," << endl;
          pseudo_preliminary_axes_TEX_ss << "axis y line*=right," << endl;
          pseudo_preliminary_axes_TEX_ss << "tick pos=right," << endl;
          pseudo_preliminary_axes_TEX_ss << "xtick={1,0.8,0.6,0.4,0.2,0}," << endl;
          pseudo_preliminary_axes_TEX_ss << "xticklabels={}," << endl;
          pseudo_preliminary_axes_TEX_ss << "xticklabel shift=6pt," << endl;
          pseudo_preliminary_axes_TEX_ss << "yticklabel shift=6pt," << endl;
          pseudo_preliminary_axes_TEX_ss << "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5";
          pseudo_preliminary_axes_TEX_ss << "}," << endl;
          if(set_y_range_exactly){pseudo_preliminary_axes_TEX_ss << "ytick distance={" << y_tick_distance_right << "}," << endl;}
          pseudo_preliminary_axes_TEX_ss << "scaled y ticks=false," << endl;
          pseudo_preliminary_axes_TEX_ss << "xmin=0," << endl;
          pseudo_preliminary_axes_TEX_ss << "xmax=1," << endl;
          //removed ymin/ymax from common settings
          pseudo_preliminary_axes_TEX_ss << "ymin=" << (showkJmolaxis?meVatom2kJmol:1.0) * ymin << "," << endl;
          pseudo_preliminary_axes_TEX_ss << "ymax=" << (showkJmolaxis?meVatom2kJmol:1.0) * ymax << "," << endl;
          // insert common settings
          pseudo_preliminary_axes_TEX_ss << common_settings_TEX_ss.str();
          pseudo_preliminary_axes_TEX_ss << "]" << endl;
          pseudo_preliminary_axes_TEX_ss << "\\end{axis}" << endl;
        }

        //third axis - meV/atom, left x label only, make other stuff color white so it's not seen (exists for positioning purposes)
        if(0){
          //this axis is NO longer needed, we use xticklabel exclusive now
          pseudo_preliminary_axes_TEX_ss << "\\begin{axis}[" << endl;
          pseudo_preliminary_axes_TEX_ss << "xlabel={";
          g_state=getUnaryGState(top_axis_point);
          if(!plot_unaries && isViableGState(g_state)) {
            if(internal_links_graph2report) {
              pseudo_preliminary_axes_TEX_ss << "\\hyperrefTitle{" << input << "_" << m_points[g_state].m_entry.auid << "}";
            }else if(no_doc && external_links) {
              pseudo_preliminary_axes_TEX_ss << "\\href{" << ENTRY_PAGE_URL_PREFIX << m_points[g_state].m_entry.auid << "}";
            }
          }
          pseudo_preliminary_axes_TEX_ss << "{" << m_velements[top_axis_point] << "}";
          pseudo_preliminary_axes_TEX_ss << "}," << endl;
          pseudo_preliminary_axes_TEX_ss << "xlabel style={font=" << label_image_font_size << ",at={(axis cs:" << top_axis_point << "," << ymin << ")}}," << endl; //{(xticklabel cs:" << top_axis_point << ")}," << endl;  //CO20200106 - patching for auto-indenting
          //add "max space between ticks=60" to change colorbar tick density
          if(showkJmolaxis){
            pseudo_preliminary_axes_TEX_ss << "axis x line*=bottom," << endl;
            pseudo_preliminary_axes_TEX_ss << "axis y line*=left," << endl;
            pseudo_preliminary_axes_TEX_ss << "tick pos=left," << endl;
          }
          pseudo_preliminary_axes_TEX_ss << "xtick={1,0.8,0.6,0.4,0.2,0}," << endl;
          pseudo_preliminary_axes_TEX_ss << "xticklabel style={opacity=0},  \%we need labels for proper padding, but avoid bolder labels with opacity==0" << endl;
          pseudo_preliminary_axes_TEX_ss << "xticklabel shift=6pt," << endl;
          pseudo_preliminary_axes_TEX_ss << "yticklabel shift=6pt," << endl;
          pseudo_preliminary_axes_TEX_ss << "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5";
          pseudo_preliminary_axes_TEX_ss << ",opacity=0}, \%we need labels for proper padding, but avoid bolder labels with opacity==0" << endl;
          if(set_y_range_exactly){pseudo_preliminary_axes_TEX_ss << "ytick distance={" << y_tick_distance << "}," << endl;}
          pseudo_preliminary_axes_TEX_ss << "scaled y ticks=false," << endl;
          pseudo_preliminary_axes_TEX_ss << "xmin=0," << endl;
          pseudo_preliminary_axes_TEX_ss << "xmax=1," << endl;
          //removed ymin/ymax from common settings
          pseudo_preliminary_axes_TEX_ss << "ymin=" << ymin << "," << endl;
          pseudo_preliminary_axes_TEX_ss << "ymax=" << ymax << "," << endl;
          // insert common settings
          pseudo_preliminary_axes_TEX_ss << common_settings_TEX_ss.str();
          pseudo_preliminary_axes_TEX_ss << "]" << endl;
          pseudo_preliminary_axes_TEX_ss << "\\end{axis}" << endl;
        }

        // real axis - meV/atom, right x label and all points data
        pseudo_preliminary_axes_TEX_ss << "\\begin{axis}[" << endl;
        //xlabel is NO longer needed, we use xticklabel exclusive now
        //pseudo_preliminary_axes_TEX_ss << "xlabel={";
        //g_state=getUnaryGState(bottom_axis_point);
        //if(!plot_unaries && isViableGState(g_state)) {
        //  if(internal_links_graph2report) {
        //    pseudo_preliminary_axes_TEX_ss << "\\hyperrefTitle{" << input << "_" << m_points[g_state].m_entry.auid << "}";
        //  }else if(no_doc && external_links) {
        //    pseudo_preliminary_axes_TEX_ss << "\\href{" << ENTRY_PAGE_URL_PREFIX << m_points[g_state].m_entry.auid << "}";
        //  }
        //}
        //pseudo_preliminary_axes_TEX_ss << "{" << m_velements[bottom_axis_point] << "}";
        //pseudo_preliminary_axes_TEX_ss << "}," << endl;
        //pseudo_preliminary_axes_TEX_ss << "xlabel style={font=" << label_image_font_size << ",at={(axis cs:" << bottom_axis_point << "," << ymin << ")}}," << endl; //{(xticklabel cs:" << bottom_axis_point << ")}," << endl;
        pseudo_preliminary_axes_TEX_ss << "ylabel={{" << (m_formation_energy_hull?"$H_{\\mathrm{f}}$ (meV/atom)":"$T_{\\mathrm{S}}$ (K)") << "}}," << endl;
        //pseudo_preliminary_axes_TEX_ss << "ylabel={{" << (m_formation_energy_hull?"formation enthalpy (meV/atom)":"entropic temperature (K)") << "}}," << endl;
        //if(m_formation_energy_hull) {pseudo_preliminary_axes_TEX_ss << "{formation enthalpy (meV/atom)}";}
        //else{pseudo_preliminary_axes_TEX_ss << "{entropic temperature (K)}";}
        //pseudo_preliminary_axes_TEX_ss << "}," << endl;
        pseudo_preliminary_axes_TEX_ss << "ylabel style={font=" << label_image_font_size << "}," << endl;
        //add "max space between ticks=60" to change colorbar tick density
        if(showkJmolaxis){
          pseudo_preliminary_axes_TEX_ss << "axis x line*=bottom," << endl;
          pseudo_preliminary_axes_TEX_ss << "axis y line*=left," << endl;
          pseudo_preliminary_axes_TEX_ss << "tick pos=left," << endl;
        }
        pseudo_preliminary_axes_TEX_ss << "xtick={1,0.8,0.6,0.4,0.2,0}," << endl;
        pseudo_preliminary_axes_TEX_ss << "xticklabels={";

        g_state=getUnaryGState(bottom_axis_point);
        if(!plot_unaries && isViableGState(g_state)) {
          if(internal_links_graph2report) {
            pseudo_preliminary_axes_TEX_ss << "\\hyperrefTitle{" << input << "_" << m_points[g_state].m_entry.auid << "}";
          }else if(no_doc && external_links) {
            pseudo_preliminary_axes_TEX_ss << "\\href{" << ENTRY_PAGE_URL_PREFIX << m_points[g_state].m_entry.auid << "}";
          }
        }
        pseudo_preliminary_axes_TEX_ss << "{" << label_image_font_size << " " << m_velements[bottom_axis_point] << "}";

        pseudo_preliminary_axes_TEX_ss << ",0.8,0.6,0.4,0.2,";

        g_state=getUnaryGState(top_axis_point);
        if(!plot_unaries && isViableGState(g_state)) {
          if(internal_links_graph2report) {
            pseudo_preliminary_axes_TEX_ss << "\\hyperrefTitle{" << input << "_" << m_points[g_state].m_entry.auid << "}";
          }else if(no_doc && external_links) {
            pseudo_preliminary_axes_TEX_ss << "\\href{" << ENTRY_PAGE_URL_PREFIX << m_points[g_state].m_entry.auid << "}";
          }
        }
        pseudo_preliminary_axes_TEX_ss << "{" << label_image_font_size << " " << m_velements[top_axis_point] << "}";

        pseudo_preliminary_axes_TEX_ss << "}," << endl;
        pseudo_preliminary_axes_TEX_ss << "xticklabel shift=6pt," << endl;
        pseudo_preliminary_axes_TEX_ss << "yticklabel shift=6pt," << endl;
        pseudo_preliminary_axes_TEX_ss << "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5";
        pseudo_preliminary_axes_TEX_ss << "}," << endl;
        if(set_y_range_exactly){pseudo_preliminary_axes_TEX_ss << "ytick distance={" << y_tick_distance << "}," << endl;}
        pseudo_preliminary_axes_TEX_ss << "scaled y ticks=false," << endl;
      }else{  // dimension==3
        //first axis - kJ/mol, only for right color bar
        if(showkJmolaxis){
          if(display_color_gradient) {
            if(include_color_bar) {
              pseudo_preliminary_axes_TEX_ss << "\\begin{ternaryaxis}[" << endl;
              // Ylabel
              pseudo_preliminary_axes_TEX_ss << "ylabel={";
              g_state=getUnaryGState(0);
              if(!plot_unaries && isViableGState(g_state)) {
                if(internal_links_graph2report) {
                  pseudo_preliminary_axes_TEX_ss << "\\hyperrefTitle{" << input << "_" << m_points[g_state].m_entry.auid << "}";
                }else if(no_doc && external_links) {
                  pseudo_preliminary_axes_TEX_ss << "\\href{" << ENTRY_PAGE_URL_PREFIX << m_points[g_state].m_entry.auid << "}";
                }
              }
              pseudo_preliminary_axes_TEX_ss << "{" << m_velements[0] << "}";
              pseudo_preliminary_axes_TEX_ss << "}," << endl;
              pseudo_preliminary_axes_TEX_ss << "ylabel style={font=" << label_image_font_size << ",at={(axis cs:0,1,0)},anchor=north east,below=" << unary_distance_axis_3D << "cm,left=" << unary_distance_axis_3D << "cm,opacity=0}, \%we need labels for proper padding, but avoid bolder labels with opacity==0" << endl;
              // Xlabel
              pseudo_preliminary_axes_TEX_ss << "xlabel={";
              g_state=getUnaryGState(1);
              if(!plot_unaries && isViableGState(g_state)) {
                if(internal_links_graph2report) {
                  pseudo_preliminary_axes_TEX_ss << "\\hyperrefTitle{" << input << "_" << m_points[g_state].m_entry.auid << "}";
                }else if(no_doc && external_links) {
                  pseudo_preliminary_axes_TEX_ss << "\\href{" << ENTRY_PAGE_URL_PREFIX << m_points[g_state].m_entry.auid << "}";
                }
              }
              pseudo_preliminary_axes_TEX_ss << "{" << m_velements[1] << "}";
              pseudo_preliminary_axes_TEX_ss << "}," << endl;
              pseudo_preliminary_axes_TEX_ss << "xlabel style={font=" << label_image_font_size << ",at={(axis cs:1,0,0)},anchor=south,above=" << sqrt(2.0*unary_distance_axis_3D*unary_distance_axis_3D) << "cm,opacity=0}, \%we need labels for proper padding, but avoid bolder labels with opacity==0" << endl;
              // Zlabel
              pseudo_preliminary_axes_TEX_ss << "zlabel={";
              g_state=getUnaryGState(2);
              if(!plot_unaries && isViableGState(g_state)) {
                if(internal_links_graph2report) {
                  pseudo_preliminary_axes_TEX_ss << "\\hyperrefTitle{" << input << "_" << m_points[g_state].m_entry.auid << "}";
                }else if(no_doc && external_links) {
                  pseudo_preliminary_axes_TEX_ss << "\\href{" << ENTRY_PAGE_URL_PREFIX << m_points[g_state].m_entry.auid << "}";
                }
              }
              pseudo_preliminary_axes_TEX_ss << "{" << m_velements[2] << "}";
              //if(!plot_unaries && isViableGState(g_state)) {
              //  if(internal_links_graph2report) {pseudo_preliminary_axes_TEX_ss << "}";}
              //  else if(no_doc && external_links) {pseudo_preliminary_axes_TEX_ss << "}";}
              //}
              pseudo_preliminary_axes_TEX_ss << "}," << endl;
              pseudo_preliminary_axes_TEX_ss << "zlabel style={font=" << label_image_font_size << ",at={(axis cs:0,0,1)},anchor=north west,below=" << unary_distance_axis_3D << "cm,right=" << unary_distance_axis_3D << "cm,opacity=0},  \%we need labels for proper padding, but avoid bolder labels with opacity==0" << endl;
              pseudo_preliminary_axes_TEX_ss << "xmin=0," << endl;
              pseudo_preliminary_axes_TEX_ss << "xmax=1," << endl;
              pseudo_preliminary_axes_TEX_ss << "ymin=0," << endl;
              pseudo_preliminary_axes_TEX_ss << "ymax=1," << endl;
              pseudo_preliminary_axes_TEX_ss << "zmin=0," << endl;
              pseudo_preliminary_axes_TEX_ss << "zmax=1," << endl;
              pseudo_preliminary_axes_TEX_ss << "ticks=none," << endl;
              //removed ymin/ymax from common settings
              pseudo_preliminary_axes_TEX_ss << common_settings_TEX_ss.str();
              // start colorbar
              // http://tex.stackexchange.com/questions/73025/how-to-change-label-and-ticks-of-a-pgfplots-colorbar
              pseudo_preliminary_axes_TEX_ss << "colorbar," << endl;
              pseudo_preliminary_axes_TEX_ss << "colorbar style={" << endl;
              pseudo_preliminary_axes_TEX_ss << "shift={(1.75cm,0cm)}," << endl; //-0cm
              pseudo_preliminary_axes_TEX_ss << "ylabel={{(kJ/mol)}}," << endl;
              pseudo_preliminary_axes_TEX_ss << "ylabel style={font=" << label_image_font_size << "}," << endl;
              pseudo_preliminary_axes_TEX_ss << "axis y line*=left," << endl;
              // ytick={18,20,25,...,45}
              pseudo_preliminary_axes_TEX_ss << "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5}," << endl;
              pseudo_preliminary_axes_TEX_ss << "yticklabel shift=6pt," << endl;
              double round_to_value_right=getRoundToValue(meVatom2kJmol*point_range);
              double y_tick_distance_right=getYTickDistance(meVatom2kJmol*y_range,approx_num_ticks,round_to_value_right);
              if(set_y_range_exactly){pseudo_preliminary_axes_TEX_ss << "ytick distance={" << y_tick_distance_right << "}," << endl;}
              pseudo_preliminary_axes_TEX_ss << "ytick align=outside," << endl;
              pseudo_preliminary_axes_TEX_ss << "scaled y ticks=false," << endl;
              pseudo_preliminary_axes_TEX_ss << "ymin=" << meVatom2kJmol * ymin << "," << endl;
              pseudo_preliminary_axes_TEX_ss << "ymax=" << meVatom2kJmol * ymax << "," << endl;
              pseudo_preliminary_axes_TEX_ss << "}," << endl;
              // end colorbar
              pseudo_preliminary_axes_TEX_ss << "colormap={mymap}{" << color_map << "}," << endl;  // https://www.sharelatex.com/learn/Pgfplots_package
              // pseudo_preliminary_axes_TEX_ss << "colorbar, colormap/jet," << endl;
              pseudo_preliminary_axes_TEX_ss << "]" << endl;
              // add hullPoints and heat map data
              // there might only be one hullPoint, which breaks the colorbar code
              //heat map data first
              num_horizontal_planes=0;
              added_header=false;
              if(show_heat_map) {
                for(uint i=0,fl_size_i=m_i_facets.size();i<fl_size_i;i++) {
                  const ChullFacet& facet=m_facets[m_i_facets[i]];
                  if(facet.m_is_artificial){continue;}
                  if(facet.m_is_vertical) {
                    if(LDEBUG){cerr << __AFLOW_FUNC__ << " NOT plane: " << abs(facet.m_normal[facet.m_normal.urows]) << endl;}
                  }else{
                    if(LDEBUG){cerr << __AFLOW_FUNC__ << " plane: " << abs(facet.m_normal[facet.m_normal.urows]) << endl;}
                    num_horizontal_planes++;
                  }
                  chull_points = facet.getCHIndices();
                  added_header=false;
                  for(uint j=0,fl_size_j=chull_points.size();j<fl_size_j;j++) {
                    const ChullPoint& point = m_points[chull_points[j]];
                    // convex hull facets color
                    if(j == 0) {
                      pseudo_preliminary_axes_TEX_ss << getPlotHeaderPDF(ADDPLOT_MODE_HEAT_MAPS,DEFAULT_PLOT_COLUMN_HEADER,display_color_gradient);
                      added_header=true;
                    }
                    pseudo_preliminary_axes_TEX_ss << getPlotPointContentPDF(point, true, true);
                  }
                  if(added_header){pseudo_preliminary_axes_TEX_ss << "};" << endl;}
                }
              }
              //now hullPoints
              added_header=false;
              for(uint i=0,fl_size_i=hull_points.size();i<fl_size_i;i++) {
                const ChullPoint& point = m_points[hull_points[i]];
                if(!point.m_has_entry) {continue;}
                if(!plot_unaries && point.isUnary()) {continue;}
                if(!added_header){
                  pseudo_preliminary_axes_TEX_ss << getPlotHeaderPDF(ADDPLOT_MODE_HULL_POINTS,"H_f_kJmol",display_color_gradient);
                  added_header=true;
                }
                pseudo_preliminary_axes_TEX_ss << getPlotPointContentPDF(point, false, true);
              }
              if(plot_unaries){
                for(uint i=0;i<dimension;i++) {
                  g_state=getUnaryGState(i);
                  if(g_state>m_points.size()-1){continue;}
                  const ChullPoint& point = m_points[g_state];
                  if(!point.m_has_entry) {continue;}
                  if(!added_header){
                    pseudo_preliminary_axes_TEX_ss << getPlotHeaderPDF(ADDPLOT_MODE_HULL_POINTS,"H_f_kJmol",display_color_gradient);
                    added_header=true;
                  }
                  pseudo_preliminary_axes_TEX_ss << getPlotPointContentPDF(point, false, true);
                }
              }
              if(added_header){pseudo_preliminary_axes_TEX_ss << "};" << endl;}
              if(plot_off_hull){
                added_header=false;
                for(uint i=0,fl_size_i=plot_points.size();i<fl_size_i;i++){
                  const ChullPoint& point=m_points[plot_points[i]];
                  if(point.m_is_g_state) {continue;}  // we already plotted these points
                  if(!added_header){
                    pseudo_preliminary_axes_TEX_ss << getPlotHeaderPDF(ADDPLOT_MODE_OFF_HULL_POINTS,DEFAULT_PLOT_COLUMN_HEADER,display_color_gradient);
                    added_header=true;
                  }
                  pseudo_preliminary_axes_TEX_ss << getPlotPointContentPDF(point, false, false);
                }
                if(added_header){pseudo_preliminary_axes_TEX_ss << "};" << endl;}
              }
              //add entry count, date, and logo here (all opacity = 0)
              if(1){
                node_option_ss << "opacity=" << 0 << ",shift={(" << watermark_x_shift << ",0cm)},anchor=north"; // << ",anchor=north west";
                //if(include_banner){node_position_ss << "axis cs:0.75,0.125,0.125";} //more complicated than necessary
                //else{
                node_position_ss << "axis cs:1,0,0";
                //}
                node_content_ss << "\\includegraphics[scale=0.25]{" << aflow_logo_skinny_file << "}";
                pseudo_preliminary_axes_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);
              }
              if(small_banner)  //image_only)
              { //CO20200106 - patching for auto-indenting
                pseudo_preliminary_axes_TEX_ss << "\\begin{pgfonlayer}{background}" << endl;  // for overlap
                node_option_ss << "opacity=0.0,anchor=south,shift={(-" << watermark_x_shift << "," << sqrt(2.0*unary_distance_axis_3D*unary_distance_axis_3D) << "cm)}";
                node_position_ss << "axis cs:1,0,0";
                node_content_ss << general_image_font_size << " " << count_total_entries << " entries"; // doesn't work with externalizing images printCountTotalEntries
                pseudo_preliminary_axes_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);

                node_option_ss << "opacity=0.0,anchor=south,shift={(" << watermark_x_shift << "," << sqrt(2.0*unary_distance_axis_3D*unary_distance_axis_3D) << "cm)}";
                node_position_ss << "axis cs:1,0,0";
                node_content_ss << general_image_font_size << " \\today";
                pseudo_preliminary_axes_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);
                pseudo_preliminary_axes_TEX_ss << "\\end{pgfonlayer}{background}" << endl;
              }
              pseudo_preliminary_axes_TEX_ss << "\\end{ternaryaxis}" << endl;
            }
          }
        }

        //real axis - meV/atom
        pseudo_preliminary_axes_TEX_ss << "\\begin{ternaryaxis}[" << endl;
        // Ylabel
        pseudo_preliminary_axes_TEX_ss << "ylabel={";
        g_state=getUnaryGState(0);
        if(!plot_unaries && isViableGState(g_state)) {
          if(internal_links_graph2report) {
            pseudo_preliminary_axes_TEX_ss << "\\hyperrefTitle{" << input << "_" << m_points[g_state].m_entry.auid << "}";
          }else if(no_doc && external_links) {
            pseudo_preliminary_axes_TEX_ss << "\\href{" << ENTRY_PAGE_URL_PREFIX << m_points[g_state].m_entry.auid << "}";
          }
        }
        pseudo_preliminary_axes_TEX_ss << "{" << m_velements[0] << "}";
        pseudo_preliminary_axes_TEX_ss << "}," << endl;
        pseudo_preliminary_axes_TEX_ss << "ylabel style={font=" << label_image_font_size << ",at={(axis cs:0,1,0)},anchor=north east,below=" << unary_distance_axis_3D << "cm,left=" << unary_distance_axis_3D << "cm}," << endl;
        // Xlabel
        pseudo_preliminary_axes_TEX_ss << "xlabel={";
        g_state=getUnaryGState(1);
        if(!plot_unaries && isViableGState(g_state)) {
          if(internal_links_graph2report) {
            pseudo_preliminary_axes_TEX_ss << "\\hyperrefTitle{" << input << "_" << m_points[g_state].m_entry.auid << "}";
          }else if(no_doc && external_links) {
            pseudo_preliminary_axes_TEX_ss << "\\href{" << ENTRY_PAGE_URL_PREFIX << m_points[g_state].m_entry.auid << "}";
          }
        }
        pseudo_preliminary_axes_TEX_ss << "{" << m_velements[1] << "}";
        pseudo_preliminary_axes_TEX_ss << "}," << endl;
        pseudo_preliminary_axes_TEX_ss << "xlabel style={font=" << label_image_font_size << ",at={(axis cs:1,0,0)},anchor=south,above=" << sqrt(2.0*unary_distance_axis_3D*unary_distance_axis_3D) << "cm}," << endl;
        // Zlabel
        pseudo_preliminary_axes_TEX_ss << "zlabel={";
        g_state=getUnaryGState(2);
        if(!plot_unaries && isViableGState(g_state)) {
          if(internal_links_graph2report) {
            pseudo_preliminary_axes_TEX_ss << "\\hyperrefTitle{" << input << "_" << m_points[g_state].m_entry.auid << "}";
          }else if(no_doc && external_links) {
            pseudo_preliminary_axes_TEX_ss << "\\href{" << ENTRY_PAGE_URL_PREFIX << m_points[g_state].m_entry.auid << "}";
          }
        }
        pseudo_preliminary_axes_TEX_ss << "{" << m_velements[2] << "}";
        pseudo_preliminary_axes_TEX_ss << "}," << endl;
        pseudo_preliminary_axes_TEX_ss << "zlabel style={font=" << label_image_font_size << ",at={(axis cs:0,0,1)},anchor=north west,below=" << unary_distance_axis_3D << "cm,right=" << unary_distance_axis_3D << "cm}," << endl;
      }
      // tikzpic_settings_TEX starts here
      tikzpic_settings_TEX_ss << "xmin=0," << endl;
      tikzpic_settings_TEX_ss << "xmax=1," << endl;
      if(dimension == 2) {
        tikzpic_settings_TEX_ss << "ymin=" << ymin << "," << endl;
        tikzpic_settings_TEX_ss << "ymax=" << ymax << "," << endl;
      }else{  // dimension==3
        tikzpic_settings_TEX_ss << "ymin=0," << endl;
        tikzpic_settings_TEX_ss << "ymax=1," << endl;
        tikzpic_settings_TEX_ss << "zmin=0," << endl;
        tikzpic_settings_TEX_ss << "zmax=1," << endl;
        tikzpic_settings_TEX_ss << "ticks=none," << endl;
      }
      //removed ymin/ymax from common settings
      tikzpic_settings_TEX_ss << common_settings_TEX_ss.str();

      ////////////////////////////////////////////////////////////////////////////
      // END Tikz pic settings
      ////////////////////////////////////////////////////////////////////////////

      ////////////////////////////////////////////////////////////////////////////
      // START Colorbar settings
      ////////////////////////////////////////////////////////////////////////////

      // won't enter if dimension==2, colorbar automatically turned off
      // start colorbar
      // http://tex.stackexchange.com/questions/73025/how-to-change-label-and-ticks-of-a-pgfplots-colorbar
      if(display_color_gradient) {
        if(include_color_bar) {
          tikzpic_settings_TEX_ss << "colorbar," << endl;
          tikzpic_settings_TEX_ss << "colorbar style={" << endl;
          tikzpic_settings_TEX_ss << "shift={(1.75cm,0cm)}," << endl;  //-0cm
          tikzpic_settings_TEX_ss << "ylabel={{" << (m_formation_energy_hull?"formation enthalpy (meV/atom)":"entropic temperature (K)") << "}}," << endl;
          tikzpic_settings_TEX_ss << "ylabel style={font=" << label_image_font_size; // << ",rotate=180";
          tikzpic_settings_TEX_ss << "}," << endl;
          if(showkJmolaxis){tikzpic_settings_TEX_ss << "axis y line*=right," << endl;}
          // ytick={18,20,25,...,45}
          tikzpic_settings_TEX_ss << "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5}," << endl;
          tikzpic_settings_TEX_ss << "yticklabel shift=6pt," << endl;
          if(set_y_range_exactly){tikzpic_settings_TEX_ss << "ytick distance={" << y_tick_distance << "}," << endl;}
          tikzpic_settings_TEX_ss << "ytick align=outside," << endl;
          tikzpic_settings_TEX_ss << "scaled y ticks=false," << endl;
          tikzpic_settings_TEX_ss << "ymin=" << ymin << "," << endl;
          tikzpic_settings_TEX_ss << "ymax=" << ymax << "," << endl;
          tikzpic_settings_TEX_ss << "}," << endl;
        }
        // end colorbar
        tikzpic_settings_TEX_ss << "colormap={mymap}{" << color_map << "}," << endl;  // https://www.sharelatex.com/learn/Pgfplots_package
        // tikzpic_settings_TEX_ss << "colorbar, colormap/jet," << endl;
      }
      tikzpic_settings_TEX_ss << "]" << endl;

      ////////////////////////////////////////////////////////////////////////////
      // END Colorbar settings
      ////////////////////////////////////////////////////////////////////////////

      ////////////////////////////////////////////////////////////////////////////
      // START Plotting everything on the hull (facet lines, colors, etc.)
      ////////////////////////////////////////////////////////////////////////////

      // for 2d, we go through points2d (in order by stoichiometry), otherwise, we
      // go through planes (order doesn't matter)
      if(dimension == 2) {
        // we use vertices vs. hull_points for facets because vertices
        // contain artificial end points for plotting

        added_header=false;
        for(uint i=0,fl_size_i=hull_points.size();i<fl_size_i;i++) {
          if(!added_header){
            tikzpic_TEX_ss << getPlotHeaderPDF(ADDPLOT_MODE_HULL_FACETS,DEFAULT_PLOT_COLUMN_HEADER,display_color_gradient);
            added_header=true;
          }
          tikzpic_TEX_ss << getPlotPointContentPDF(m_points[hull_points[i]], true, true);
        }
        if(added_header){tikzpic_TEX_ss << "};" << endl;}

        // we plot hullpoints for points on hull ONLY (not lines)
        added_header=false;
        for(uint i=0,fl_size_i=hull_points.size();i<fl_size_i;i++) {
          const ChullPoint point=m_points[hull_points[i]];
          if(!point.m_has_entry) {continue;}
          if(!plot_unaries && point.isUnary()) {continue;}
          if(!added_header){
            tikzpic_TEX_ss << getPlotHeaderPDF(ADDPLOT_MODE_HULL_POINTS,DEFAULT_PLOT_COLUMN_HEADER,display_color_gradient);
            added_header=true;
          }
          tikzpic_TEX_ss << getPlotPointContentPDF(point, false, true);
        }
        if(plot_unaries){
          for(uint i=0;i<dimension;i++) {
            g_state=getUnaryGState(i);
            if(g_state>m_points.size()-1){continue;}
            const ChullPoint& point = m_points[g_state];
            if(!point.m_has_entry) {continue;}
            if(!added_header){
              tikzpic_TEX_ss << getPlotHeaderPDF(ADDPLOT_MODE_HULL_POINTS,DEFAULT_PLOT_COLUMN_HEADER,display_color_gradient);
              added_header=true;
            }
            tikzpic_TEX_ss << getPlotPointContentPDF(point, false, true);
          }
        }
        if(added_header){tikzpic_TEX_ss << "};" << endl;}
      }else{  // dimension==3
        //////////////////////////////////////////////////////////////////////////
        // START 3D Heat maps and facet lines
        //////////////////////////////////////////////////////////////////////////

        // this option allows for overlapping facet lines (may or may not be desired)
        if(draw_all_facet_lines) {
          num_horizontal_planes=0;
          for(uint i=0,fl_size_i=m_i_facets.size();i<fl_size_i;i++) {
            const ChullFacet& facet=m_facets[m_i_facets[i]];
            if(facet.m_is_artificial){continue;}
            if(facet.m_is_vertical) {
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " NOT plane: " << abs(facet.m_normal[facet.m_normal.urows]) << endl;}
            }else{
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " plane: " << abs(facet.m_normal[facet.m_normal.urows]) << endl;}
              num_horizontal_planes++;
            }
            chull_points = facet.getCHIndices();
            chull_points.push_back(chull_points[0]);  // that way we get full facet
            added_header=false;
            for(uint j=0,fl_size_j=chull_points.size();j<fl_size_j;j++) {
              const ChullPoint& point = m_points[chull_points[j]];
              // keep first coord
              if(j == 0) {
                // convex hull facets color
                heat_map_TEX_ss << getPlotHeaderPDF(ADDPLOT_MODE_HEAT_MAPS,DEFAULT_PLOT_COLUMN_HEADER,display_color_gradient);
                // convex hull facets (lines)
                convex_hull_facets_TEX_ss << getPlotHeaderPDF(ADDPLOT_MODE_HULL_FACETS,DEFAULT_PLOT_COLUMN_HEADER,display_color_gradient);
                // thick white line for contrast
                convex_hull_facets_drop_shadow_TEX_ss << getPlotHeaderPDF(ADDPLOT_MODE_HULL_FACETSDROP_SHADOWS,DEFAULT_PLOT_COLUMN_HEADER,display_color_gradient);
                added_header=true;
              }
              convex_hull_facets_TEX_ss << getPlotPointContentPDF(point, true, true);
              convex_hull_facets_drop_shadow_TEX_ss << getPlotPointContentPDF(point, true, true);
              if(j!=chull_points.size()-1) {  // heatmap doesn't need first coord repeated
                heat_map_TEX_ss << getPlotPointContentPDF(point, true, true);
              }
            }
            // end line for connecting facet lines
            if(added_header){
              convex_hull_facets_TEX_ss << "};" << endl;
              convex_hull_facets_drop_shadow_TEX_ss << "};" << endl;
              heat_map_TEX_ss << "};" << endl;
            }
          }
        }else{
          // no choice but to go through planes twice, once for
          // heatmap, the other for facet lines/drop shadows

          // heatmap
          num_horizontal_planes=0;
          if(show_heat_map) {
            for(uint i=0,fl_size_i=m_i_facets.size();i<fl_size_i;i++) {
              const ChullFacet& facet=m_facets[m_i_facets[i]];
              if(facet.m_is_artificial){continue;}
              if(facet.m_is_vertical) {
                if(LDEBUG){cerr << __AFLOW_FUNC__ << " NOT plane: " << abs(facet.m_normal[facet.m_normal.urows]) << endl;}
              }else{
                if(LDEBUG){cerr << __AFLOW_FUNC__ << " plane: " << abs(facet.m_normal[facet.m_normal.urows]) << endl;}
                num_horizontal_planes++;
              }
              chull_points = facet.getCHIndices();
              added_header=false;
              for(uint j=0,fl_size_j=chull_points.size();j<fl_size_j;j++) {
                const ChullPoint& point = m_points[chull_points[j]];
                // convex hull facets color
                if(j == 0) {
                  heat_map_TEX_ss << getPlotHeaderPDF(ADDPLOT_MODE_HEAT_MAPS,DEFAULT_PLOT_COLUMN_HEADER,display_color_gradient);
                  added_header=true;
                }
                heat_map_TEX_ss << getPlotPointContentPDF(point, true, true);
              }
              if(added_header){heat_map_TEX_ss << "};" << endl;}
            }
          }
          // save time by not saving WHOLE chullPoints, just compositional part of
          // xvector for facet lines
          //vector<vector<xvector<double> > > facet_lines, facet_lines_dropshadow;
          //fix UNWANTEDFACETLINE to check indices, not actual coords
          num_horizontal_planes=0;
          for(uint ii=0,fl_size_ii=facet_lines.size();ii<fl_size_ii;ii++){facet_lines[ii].clear();} facet_lines.clear();
          for(uint ii=0,fl_size_ii=facet_lines_dropshadow.size();ii<fl_size_ii;ii++){facet_lines_dropshadow[ii].clear();} facet_lines_dropshadow.clear();
          for(uint i=0,fl_size_i=m_i_facets.size();i<fl_size_i;i++) {
            const ChullFacet& facet=m_facets[m_i_facets[i]];
            if(facet.m_is_artificial){continue;}
            if(facet.m_is_vertical) {
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " NOT plane: " << abs(facet.m_normal[facet.m_normal.urows]) << endl;}
            }else{
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " plane: " << abs(facet.m_normal[facet.m_normal.urows]) << endl;}
              num_horizontal_planes++;
            }
            chull_points = facet.getCHIndices();
            chull_points.push_back(chull_points[0]);  // that way we get full facet
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking for all unwanted facets" << endl;}
            for(uint l=0,fl_size_l=chull_points.size();l<fl_size_l - 1;l++) {
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking at point l=" << l << " and l=" << l + 1 << endl;}
              if(!unwantedFacetLine(chull_points[l], chull_points[l+1], facet_lines, true)) {
                if(LDEBUG){cerr << __AFLOW_FUNC__ << " plotting point l=" << l << " and l=" << l + 1 << endl;}
                // convex hull facets (lines)
                convex_hull_facets_TEX_ss << getPlotHeaderPDF(ADDPLOT_MODE_HULL_FACETS,DEFAULT_PLOT_COLUMN_HEADER,display_color_gradient);
                for(uint j=0;j<2;j++) {
                  const ChullPoint& point = m_points[chull_points[l+j]];
                  convex_hull_facets_TEX_ss << getPlotPointContentPDF(point, true, true);
                }
                // end line for connecting facet lines
                convex_hull_facets_TEX_ss << "};" << endl;
              }else{
                if(LDEBUG){
                  cerr << __AFLOW_FUNC__ << " skipping l=" << l << " and l=" << l + 1 << endl;
                  cerr << __AFLOW_FUNC__ << " l=" << l << " is " << m_points[chull_points[l]].m_coords << endl;
                  cerr << __AFLOW_FUNC__ << " l=" << l+1 << " is " << m_points[chull_points[l+1]].m_coords << endl;
                }
              }
              if(hull_drop_shadow) {
                if(!unwantedFacetLine(chull_points[l], chull_points[l+1], facet_lines_dropshadow, false)){
                  if(LDEBUG){cerr << __AFLOW_FUNC__ << " plotting drop shadow l=" << l << " and l=" << l + 1 << endl;}
                  // thick white line for contrast
                  convex_hull_facets_drop_shadow_TEX_ss << getPlotHeaderPDF(ADDPLOT_MODE_HULL_FACETSDROP_SHADOWS,DEFAULT_PLOT_COLUMN_HEADER,display_color_gradient);
                  for(uint j=0;j<2;j++) {
                    const ChullPoint& point = m_points[chull_points[l+j]];
                    convex_hull_facets_drop_shadow_TEX_ss << getPlotPointContentPDF(point, true, true);
                  }
                  // end line for connecting facet lines
                  convex_hull_facets_drop_shadow_TEX_ss << "};" << endl;
                }else{
                  if(LDEBUG){
                    cerr << __AFLOW_FUNC__ << " skipping drop shadow l=" << l << " and l=" << l + 1 << endl;
                    cerr << __AFLOW_FUNC__ << " l=" << l << " is " << m_points[chull_points[l]].m_coords << endl;
                    cerr << __AFLOW_FUNC__ << " l=" << l+1 << " is " << m_points[chull_points[l+1]].m_coords << endl;
                  }
                }
              }
            }
          }
        }

        //////////////////////////////////////////////////////////////////////////
        // END 3D Heat maps and facet lines
        //////////////////////////////////////////////////////////////////////////

        // hullPoints only, so they don't repeat
        added_header=false;
        for(uint i=0,fl_size_i=hull_points.size();i<fl_size_i;i++) {
          const ChullPoint& point = m_points[hull_points[i]];
          if(!point.m_has_entry) {continue;}
          if(!plot_unaries && point.isUnary()) {continue;}
          if(!added_header){
            convex_hull_vertices_TEX_ss << getPlotHeaderPDF(ADDPLOT_MODE_HULL_POINTS,DEFAULT_PLOT_COLUMN_HEADER,display_color_gradient);
            added_header=true;
          }
          convex_hull_vertices_TEX_ss << getPlotPointContentPDF(point, false, true);
        }
        if(plot_unaries){
          for(uint i=0;i<dimension;i++) {
            g_state=getUnaryGState(i);
            if(g_state>m_points.size()-1){continue;}
            const ChullPoint& point = m_points[g_state];
            if(!point.m_has_entry) {continue;}
            if(!added_header){
              convex_hull_vertices_TEX_ss << getPlotHeaderPDF(ADDPLOT_MODE_HULL_POINTS,DEFAULT_PLOT_COLUMN_HEADER,display_color_gradient);
              added_header=true;
            }
            convex_hull_vertices_TEX_ss << getPlotPointContentPDF(point, false, true);
          }
        }
        if(added_header){convex_hull_vertices_TEX_ss << "};" << endl;}

        // big merge
        // white lines first, then black
        if(hull_drop_shadow) {tikzpic_TEX_ss << convex_hull_facets_drop_shadow_TEX_ss.str();}
        tikzpic_TEX_ss << convex_hull_facets_TEX_ss.str();
        tikzpic_TEX_ss << convex_hull_vertices_TEX_ss.str();
        convex_hull_facets_TEX_ss.str("");
        convex_hull_vertices_TEX_ss.str("");
        convex_hull_facets_drop_shadow_TEX_ss.str("");
      }

      ////////////////////////////////////////////////////////////////////////////
      // END Plotting everything on the hull (facet lines, colors, etc.)
      ////////////////////////////////////////////////////////////////////////////

      ////////////////////////////////////////////////////////////////////////////
      // START Plotting off hull points
      ////////////////////////////////////////////////////////////////////////////

      if(LDEBUG){cerr << __AFLOW_FUNC__ << " start plotting off hull points" << endl;}
      if(plot_off_hull || !no_labels || internal_links_graph2report || no_doc) {
        string compound_label;
        bool found_icsd_label=false;
        for(uint i=0,fl_size_i=plot_points.size();i<fl_size_i;i++){
          const ChullPoint& point=m_points[plot_points[i]];
          const aflowlib::_aflowlib_entry& entry = point.m_entry;
          // get coords for table
          if(plot_off_hull) {
            if(!point.m_is_g_state) {  // we already plotted these points
              points_data_ss << getPlotPointContentPDF(point, false, false);
            }
          }
          if(!plot_off_hull && !point.m_is_g_state) {continue;}

          //////////////////////////////////////////////////////////////////////
          // START Creating clickable links at points
          //////////////////////////////////////////////////////////////////////

          if(internal_links_graph2report || (no_doc && external_links)) {
            node_option_ss << "opacity=0.0";// get node option
            node_position_ss << getNodeCoordPosition(point); // get node position
            // get node content
            node_content_ss << "\\tiny";
            node_content_ss << "{";
            if(no_doc && external_links) {node_content_ss << "\\href{" << ENTRY_PAGE_URL_PREFIX << entry.auid << "}";}
            else{
              //uint i_coord_group;
              if(getCoordGroupIndex(point,i_coord_group)){
                uint ref_state=m_coord_groups[i_coord_group].m_ref_state;
                if(isViablePoint(ref_state)){
                  node_content_ss << "\\hyperref[" << input << "_" << m_points[ref_state].m_entry.auid << "]";  // for hyperref to right part of document
                }
              }
            }
            node_content_ss << "{O}";
            node_content_ss << "}";
            labels_data_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss); // create node
          }

          //////////////////////////////////////////////////////////////////////
          // END Creating clickable links at points
          //////////////////////////////////////////////////////////////////////

          // now we create labels
          // do not create another node for unaries unless metalabels or
          // protolables
          if(!no_labels && !(!labels_off_hull && !point.m_is_g_state) && !(point.isUnary() && !(meta_labels || prototype_labels))) {
            ////////////////////////////////////////////////////////////////////
            // START Meta-labels (lots of info)
            ////////////////////////////////////////////////////////////////////

            if(meta_labels) {
              uint precision_tmp=0;
              double tmp_roundoff_tol=5.0*pow(10,-((int)precision_tmp)-1);
              // no node option

              // get node position
              node_position_ss << getNodeCoordPosition(point);

              // get node content
              if(dimension == 3 && entry.vspecies.size() == 3) {node_content_ss << "\\color{" << ternary_label_color << "}";}

              node_content_ss << "\\tiny{";   //these are sort of "special" font settings meant to reduce clutter, no need to add to aflowrc
              node_content_ss << "\\shortstack{";
              compound_label = prettyPrintCompound(point,(plot_reduced_composition?gcd_vrt:no_vrt),true,latex_ft);  // don't print prototype if not
              node_content_ss << compound_label;
              // equal to compound
              output_name=prettyPrintPrototype(point,false,false);  // only one backslash needed
              if(compound_label != output_name) {
                node_content_ss << "::";
                node_content_ss << output_name;
              }
              // shortstack newline
              node_content_ss << "\\\\";
              // enthalpy of formation, row 4
              // no need for precision for next few columns, leave it same way
              // as received from AFLOW
              node_content_ss << "$H_{\\mathrm{f}}$=" << aurostd::utype2string(chull::H_f_atom(point,_m_),precision_tmp,true,tmp_roundoff_tol,FIXED_STREAM) << " meV/atom";
              //num_ss << chull::H_f_atom(point,_m_);
              //node_content_ss << "$H_{\\mathrm{f}}$=" << num_ss.str() << " meV/atom";
              //num_ss.str("");
              // shortstack newline
              node_content_ss << "\\\\";
              // entropic temperature, row 5
              node_content_ss << "$T_{\\mathrm{S}}$=" << aurostd::utype2string(chull::T_S(point),precision_tmp,true,tmp_roundoff_tol,FIXED_STREAM) << " K";
              //num_ss << chull::T_S(point);
              //node_content_ss << "$T_{\\mathrm{S}}$=" << num_ss.str() << " K";
              //num_ss.str("");
              // dist_2_hull, row 6
              if(!point.isGState()) {
                // shortstack newline
                node_content_ss << "\\\\";
                if(m_formation_energy_hull) {
                  node_content_ss << "$" << getDelta(helvetica_font) << " H_{\\mathrm{hull}}$=" << aurostd::utype2string(point.getDist2Hull(_m_),precision_tmp,true,tmp_roundoff_tol,FIXED_STREAM) << " meV/atom";  //CHULL_PRECISION
                }else{
                  node_content_ss << "$" << getDelta(helvetica_font) << " T_{\\mathrm{S}}$=" << aurostd::utype2string(point.getDist2Hull(_std_),precision_tmp,true,tmp_roundoff_tol,FIXED_STREAM) << " K"; //CHULL_PRECISION
                }
              }
              node_content_ss << "}";
              node_content_ss << "}";
              // create node
              labels_data_ss << nodeCreator(node_option_ss, node_position_ss,node_content_ss);

              //////////////////////////////////////////////////////////////////
              // END Meta-labels (lots of info)
              //////////////////////////////////////////////////////////////////

            }else{  // no metadata
              // special case
              // create two nodes, one with compound name below ground-state
              // node, another with prototype above each node

              //////////////////////////////////////////////////////////////////
              // START Special case labels (2D with prototypes)
              //////////////////////////////////////////////////////////////////

              if(dimension == 2 && !plot_labels && prototype_labels && rotate_labels) {
                if(point.m_is_g_state && !point.isUnary()) {
                  // get node option
                  node_option_ss << "rotate=90,anchor=" << (m_formation_energy_hull?"east":"west");
                  // get node position
                  node_position_ss << getNodeCoordPosition(point);
                  // get node content
                  node_content_ss << "\\" << font_size << "{";  //leave g-states as default font
                  if(bold_labels/* && !helvetica_font*/) {node_content_ss << "\\textbf{";}
                  if(!m_formation_energy_hull){node_content_ss << "~~";}   //pre - 1
                  node_content_ss << prettyPrintCompound(point,(plot_reduced_composition?gcd_vrt:no_vrt),true,latex_ft);
                  if(m_formation_energy_hull){node_content_ss << "~~";}  //post - 2: only works if we put two?
                  if(bold_labels/* && !helvetica_font*/) {node_content_ss << "}";}
                  node_content_ss << "}";
                  // create node
                  labels_data_ss << nodeCreator(node_option_ss, node_position_ss,node_content_ss);
                }
                // prototype label
                // get node option
                if(point.isUnary()) {
                  if(entry.vspecies[0] == m_velements[0]) {node_option_ss << "anchor=east";}
                  else{node_option_ss << "anchor=west";}
                }else{node_option_ss << "anchor=" << (m_formation_energy_hull?"south,above=0.1cm":"north,below=0.1cm");}
                // get node position
                node_position_ss << getNodeCoordPosition(point);
                // get node content
                node_content_ss << "\\scriptsize{"; //these are sort of "special" font settings meant to reduce clutter, no need to add to aflowrc
                output_name=prettyPrintPrototype(point,false,false);  // only one backslash needed
                if(point.isUnary() && entry.vspecies[0] == m_velements[1]) {node_content_ss << "~~";}
                node_content_ss << output_name;
                if(point.isUnary() && entry.vspecies[0] == m_velements[0]) {node_content_ss << "~~";}
                // enclose brackets
                node_content_ss << "}";
                // create node
                labels_data_ss << nodeCreator(node_option_ss, node_position_ss,node_content_ss);

                ////////////////////////////////////////////////////////////////
                // END Special case labels (2D with prototypes)
                ////////////////////////////////////////////////////////////////

              }else{  // normal case

                ////////////////////////////////////////////////////////////////
                // START Normal labels
                ////////////////////////////////////////////////////////////////

                // get node option
                if(rotate_labels) {
                  if(dimension == 2) {
                    if(!point.isUnary()) {  //deal with unaries later
                      if(point.m_is_g_state) {node_option_ss << "rotate=90,anchor=" << (m_formation_energy_hull?"east":"west");}
                      else{node_option_ss << "anchor=" << (m_formation_energy_hull?"south,above=0.1cm":"north,below=0.1cm");}
                    }
                  }else{  // dimension==3
                    // ordered to optimize speed (binaries are bulk, then ternaries, then unaries)
                    if(entry.vspecies.size() == 2) {
                      // / side of triangle
                      if(entry.vspecies[0] == m_velements[0] && entry.vspecies[1] == m_velements[1]) {
                        node_option_ss << "rotate=-30,anchor=east";
                        // \ side of triangle
                      }else if(entry.vspecies[0] == m_velements[1] && entry.vspecies[1] == m_velements[2]) {
                        node_option_ss << "rotate=30,anchor=west";
                        // _ side of triangle
                      }else if(entry.vspecies[0] == m_velements[0] && entry.vspecies[1] == m_velements[2]) {
                        node_option_ss << "rotate=90,anchor=east";
                      }
                    }else if(entry.vspecies.size() == 3) {
                      if(!show_heat_map) {node_option_ss << "anchor=south,above=0.1cm";}
                    }else if(entry.vspecies.size() == 1) {
                      // bottom two edges of triangle
                      if(entry.vspecies[0] == m_velements[0] || entry.vspecies[0] == m_velements[2]) {
                        node_option_ss << "anchor=north,below=0.1cm";
                        // top of triangle
                      }else{node_option_ss << "anchor=south,above=0.1cm";}
                    }
                  }
                }
                if(dimension == 2) {  //deal with unaries separately
                  if(point.isUnary()) {
                    if(entry.vspecies[0] == m_velements[0]) {node_option_ss << "anchor=east";}
                    else{node_option_ss << "anchor=west";}
                  }
                }
                // get node position
                node_position_ss << getNodeCoordPosition(point);
                // get node content
                if(dimension == 3 && entry.vspecies.size() == 3) {node_content_ss << "\\color{" << ternary_label_color << "}";}
                if(compound_labels && prototype_labels) {node_content_ss << "\\scriptsize{";} //these are sort of "special" font settings meant to reduce clutter, no need to add to aflowrc
                else{node_content_ss << "\\" << font_size << "{";}  //already gets wrapped with brackets in nodeCreator(), just leave alone so you mitigate spaces
                if(bold_labels || (entry.vspecies.size() == 3 && bold_labels_ternaries)/* && !helvetica_font*/) {node_content_ss << "\\textbf{";}

                // add space for rotation
                if(rotate_labels) {
                  if(dimension == 2) {
                    if(!point.isUnary()){ //handle unaries later
                      if(!m_formation_energy_hull && point.m_is_g_state){node_content_ss << "~";}
                    }
                  }else{  // dimension==3
                    if(entry.vspecies.size() == 2) {
                      // \ side of triangle
                      if(entry.vspecies[0] == m_velements[1] && entry.vspecies[1] == m_velements[2]) {
                        node_content_ss << "~";
                      }
                    }
                  }
                }
                if(dimension == 2) {
                  if(point.isUnary()){  //unaries
                    if(entry.vspecies[0] == m_velements[1]) {node_content_ss << "~";}
                  }
                }

                found_icsd_label=false;
                if(icsd_labels){  //icsd_labels mutually exclusive with compound/prototype
                  //uint i_coord_group;
                  if(getCoordGroupIndex(point,i_coord_group)){
                    uint canonical_icsd=m_coord_groups[i_coord_group].m_i_icsd_g_state;
                    if(isViablePoint(canonical_icsd)){
                      compound_label = prettyPrintPrototype(m_points[canonical_icsd],false,true);  // only one backslash needed  //icsd_skim_label=true
                      if(!compound_label.empty()){
                        node_content_ss << compound_label;
                        found_icsd_label=true;
                      }
                    }
                  }
                }
                if(!(icsd_labels&&found_icsd_label)){
                  if(compound_labels) {
                    compound_label = prettyPrintCompound(point,(plot_reduced_composition?gcd_vrt:no_vrt),true,latex_ft);  // don't print prototype if
                    node_content_ss << compound_label;
                    // equal to compound
                  }
                  if(prototype_labels) {
                    output_name=prettyPrintPrototype(point,false,false);  // only one backslash needed
                    if(compound_labels) {
                      if(compound_label != output_name) {
                        node_content_ss << "::";
                        node_content_ss << output_name;
                      }
                    }else{node_content_ss << output_name;}
                  }
                }

                if(rotate_labels) {
                  if(dimension == 2) {
                    if(!point.isUnary()){ //handle unaries later
                      if(m_formation_energy_hull && point.m_is_g_state){node_content_ss << "~~";} // only works if we put two?
                    }
                  }else{  // dimension==3
                    if(entry.vspecies.size() == 2) {
                      // / side of triangle
                      if(entry.vspecies[0] == m_velements[0] && entry.vspecies[1] == m_velements[1]) {
                        node_content_ss << "~~";  // only works if we put two?
                        // _ side of triangle
                      }else if(entry.vspecies[0] == m_velements[0] && entry.vspecies[1] == m_velements[2]) {
                        node_content_ss << "~~";  // only works if we put two?
                      }
                    }
                  }
                }
                if(dimension == 2) {
                  if(point.isUnary()){ //unaries
                    if(entry.vspecies[0] == m_velements[0]) {node_content_ss << "~~";}
                  }
                }

                // enclose brackets
                if(bold_labels || (entry.vspecies.size() == 3 && bold_labels_ternaries)/* && !helvetica_font*/) {node_content_ss << "}";}
                node_content_ss << "}";
                // create node
                labels_data_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);

                ////////////////////////////////////////////////////////////////
                // END Normal labels
                ////////////////////////////////////////////////////////////////
              }
            }
          }
        }
        // merge labels data
        if(plot_off_hull && !points_data_ss.str().empty()) {
          tikzpic_TEX_ss << getPlotHeaderPDF(ADDPLOT_MODE_OFF_HULL_POINTS,DEFAULT_PLOT_COLUMN_HEADER,display_color_gradient); // header for points off hull
          tikzpic_TEX_ss << points_data_ss.str();
          tikzpic_TEX_ss << "};" << endl;
        }
        //add isomax data
        if(plot_iso_max_latent_heat){
          string color;
          uint g_point;
          vector<string> latex_colors=grabAcceptableLatexColors(LATEX_COLORS_TO_AVOID,true,true,g_states.size());
          if(latex_colors.size()!=g_states.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Issue grabbing latex colors");}
          for(uint gi=0,fl_size_gi=g_states.size();gi<fl_size_gi;gi++){
            color=latex_colors[gi];
            //add curves
            tikzpic_TEX_ss << "\\addplot[domain=0.001:0.999,smooth,variable=\\x," << color << "] {";
            g_point=g_states[gi];
            tikzpic_TEX_ss << ((double)KBOLTZEV) << "*" << 1.0e3 << "*" << T_S(m_points[g_point]) << "*(\\x*ln(\\x)+(1-\\x)*ln(1-\\x))};" << endl;
            //add labels
            labels_data_ss << "\\node [fill=white] at (axis cs:0.5," << isoMaxLatentHeat(m_points[g_point],0.5,_m_);
            labels_data_ss << ") {\\" << font_size << "{$T_{\\mathrm{S}}=" << int(T_S(m_points[g_point])+0.5) << "~\\mathrm{K}$}};" << endl;
          }
        }
        if(!labels_data_ss.str().empty()) {
          tikzpic_TEX_ss << "\\begin{pgfonlayer}{foreground}" << endl;  // for overlap
          tikzpic_TEX_ss << labels_data_ss.str();
          tikzpic_TEX_ss << "\\end{pgfonlayer}{foreground}" << endl;
        }
      }
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " start plotting logos" << endl;}
      if(1){
        //add node with logo inside first
        tikzpic_TEX_ss << "\\begin{pgfonlayer}{background}" << endl;  // for overlap
        //tikzpic_TEX_ss << "\\begin{scope}[overlay,on background layer]" << endl;
        if(dimension==2){
          ////top left
          //node_option_ss << "opacity=" << opacity_watermark;
          //node_position_ss << "axis cs:0.1," << ymin*0.05;
          //node_content_ss << "\\includegraphics[scale=0.15]{" << aflow_logo_skinny_file << "}";
          //tikzpic_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);
          ////top right
          //node_option_ss << "opacity=" << opacity_watermark;
          //node_position_ss << "axis cs:0.9," << ymin*0.05;
          //node_content_ss << "\\includegraphics[scale=0.15]{" << aflow_logo_skinny_file << "}";
          //tikzpic_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);
          ////bottom left
          //node_option_ss << "opacity=" << opacity_watermark;
          //node_position_ss << "axis cs:0.1," << ymin*0.95;
          //node_content_ss << "\\includegraphics[scale=0.15]{" << aflow_logo_skinny_file << "}";
          //tikzpic_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);
          //bottom right
          //node_option_ss << "opacity=" << opacity_watermark;
          //node_position_ss << "axis cs:0.9," << ymin*0.95;
          //node_content_ss << "\\includegraphics[scale=0.15]{" << aflow_logo_skinny_file << "}";
          //tikzpic_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);
          ////center
          //node_option_ss << "opacity=" << opacity_watermark;
          //node_position_ss << "axis cs:0.5," << ymin*0.5;
          //node_content_ss << "\\includegraphics[scale=0.15]{" << aflow_logo_skinny_file << "}";
          //tikzpic_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);
          //node_position_ss << "axis cs:0.125," << ymin*0.95;
          node_option_ss << "opacity=" << opacity_watermark << ",anchor=south west";
          node_position_ss << "axis cs:0," << ymin;
          node_content_ss << "\\includegraphics[scale=0.25]{" << aflow_logo_skinny_file << "}";
          tikzpic_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);
        }else{  //dimension==3
          //// / side
          //node_option_ss << "opacity=" << opacity_watermark << ",rotate=60,anchor=south";
          //node_position_ss << "axis cs:0.5,0.5,0";
          //node_content_ss << "\\includegraphics[scale=0.15]{" << aflow_logo_skinny_file << "}";
          //tikzpic_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);
          //// \ side
          //node_option_ss << "opacity=" << opacity_watermark << ",rotate=-60,anchor=south";
          //node_position_ss << "axis cs:0.5,0,0.5";
          //node_content_ss << "\\includegraphics[scale=0.15]{" << aflow_logo_skinny_file << "}";
          //tikzpic_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);
          //// _ side
          //node_option_ss << "opacity=" << opacity_watermark << ",anchor=north";
          //node_position_ss << "axis cs:0,0.5,0.5";
          //node_content_ss << "\\includegraphics[scale=0.15]{" << aflow_logo_skinny_file << "}";
          //tikzpic_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);
          //node_position_ss << "axis cs:0,0.5,0.5";
          //node_option_ss << ",shift={(0in,-1in)}";
          if(!small_banner) //&& include_color_bar)
          { //CO20200106 - patching for auto-indenting
            node_option_ss << "opacity=" << opacity_watermark << ",shift={(" << watermark_x_shift << ",0cm)},anchor=north"; // << ",anchor=north west";
            //if(include_color_bar){node_position_ss << "axis cs:0.75,0.125,0.125";}  //more complicated than necessary
            //else{
            node_position_ss << "axis cs:1,0,0";
            //}
          }else{
            node_option_ss << "opacity=" << opacity_watermark << ",shift={(-" << watermark_x_shift << ",0cm)},anchor=north"; // << ",anchor=north west";
            node_position_ss << "axis cs:1,0,0";
          }
          node_content_ss << "\\includegraphics[scale=0.25]{" << aflow_logo_skinny_file << "}";
          tikzpic_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);
        }
        if(small_banner) //image_only) //note: this must also be added to pseudo_preliminary_axes_TEX_ss for proper padding (color bar), ONLY AFFECTS 3D
        { //CO20200106 - patching for auto-indenting
          if(dimension==2){
            node_option_ss << "anchor=south west";
            node_position_ss << "axis cs:0," << ymax;
            node_content_ss << general_image_font_size << " " << count_total_entries << " entries"; // doesn't work with externalizing images printCountTotalEntries
            tikzpic_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);

            node_option_ss << "anchor=south east";
            node_position_ss << "axis cs:1," << ymax;
            node_content_ss << general_image_font_size << " \\today";
            tikzpic_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);
          }else{  //dimension==3
            //node_option_ss << "shift={(0in,5in)}";
            //node_position_ss << "axis cs:0,1,0";
            node_option_ss << "anchor=south,shift={(-" << watermark_x_shift << "," << sqrt(2.0*unary_distance_axis_3D*unary_distance_axis_3D) << "cm)}";
            node_position_ss << "axis cs:1,0,0";
            node_content_ss << general_image_font_size << " " << count_total_entries << " entries"; // doesn't work with externalizing images printCountTotalEntries
            tikzpic_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);

            //node_option_ss << "shift={(0in,5in)}";
            //node_position_ss << "axis cs:0,0,1";
            node_option_ss << "anchor=south,shift={(" << watermark_x_shift << "," << sqrt(2.0*unary_distance_axis_3D*unary_distance_axis_3D) << "cm)}";
            node_position_ss << "axis cs:1,0,0";
            node_content_ss << general_image_font_size << " \\today";
            tikzpic_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);
          }
        }
        tikzpic_TEX_ss << "\\end{pgfonlayer}{background}" << endl;
        //tikzpic_TEX_ss << "\\end{scope}" << endl;
      }

      ////////////////////////////////////////////////////////////////////////////
      // END Plotting off hull points
      ////////////////////////////////////////////////////////////////////////////

      if(dimension == 2) {tikzpic_TEX_ss << "\\end{axis}" << endl;}
      else{tikzpic_TEX_ss << "\\end{ternaryaxis}" << endl;}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " ending axis" << endl;}

      ////////////////////////////////////////////////////////////////////////////
      // END Axis
      ////////////////////////////////////////////////////////////////////////////

      //this section sucks to maintain, I fix with fancyheaders
      //I only keep here as reference to some potentially useful code
      //WARNING: the banner is EASILY the most complicated part of the code, as it's very
      //options dependent.
      //modify at your own peril 
      //if(0&&!no_banner) { //gone!
      //  if(dimension == 2) {
      //    ////////////////////////////////////////////////////////////////////////
      //    // START Logo on top left
      //    ////////////////////////////////////////////////////////////////////////
      //
      //    if(image_only){tikzpic_TEX_ss << "\\begin{scope}[remember picture,on background layer]" << endl;}
      //    else{tikzpic_TEX_ss << "\\begin{scope}[remember picture,overlay,on background layer]" << endl;}
      //    
      //    if(!small_banner) {
      //      // get node option
      //      if(print_aflow_logo_full) {
      //        if(image_only) {
      //          if(print_logo_2){node_option_ss << "shift={(2.75cm,0.75cm)}";}  //DONE
      //          else{node_option_ss << "shift={(2.75cm,0.5cm)}";} //DONE
      //        }else{node_option_ss << "shift={(2.25cm,0.5cm)}";} //DONE
      //      }else{
      //        if(image_only){node_option_ss << "shift={(4.0cm,0.75cm)}";} //DONE
      //        else{node_option_ss << "shift={(3.25cm,0.75cm)}";} //DONE
      //			}
      //      // get node position
      //      node_position_ss << "current bounding box.north west";
      //      // get node content
      //      if(print_aflow_logo_full) {
      //        if(external_links) {node_content_ss << "\\href{" << AFLOW_WEB << "}";}
      //        node_content_ss << "{\\includegraphics[scale=0.2]{" << aflow_logo_full_file << "}}";
      //      }else{
      //        node_content_ss << "\\shortstack[l]{";
      //        if(external_links) {node_content_ss << "\\href{" << AFLOW_WEB << "}";}
      //        node_content_ss << "{\\LARGE AFLOW V" << string(AFLOW_VERSION) << "}";
      //        node_content_ss << "\\\\";
      //        node_content_ss << "{\\large C. Oses and S. Curtarolo}";
      //        node_content_ss << "}";
      //      }
      //      // create node
      //      tikzpic_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);
      //    }
      //
      //    ////////////////////////////////////////////////////////////////////////
      //    // END Logo on top left
      //    ////////////////////////////////////////////////////////////////////////
      //
      //    ////////////////////////////////////////////////////////////////////////
      //    // START AFLOW.org on top right
      //    // nomad logo too
      //    ////////////////////////////////////////////////////////////////////////
      //
      //    // get node option
      //    if(image_only){
      //      if(!small_banner&&print_logo_2){
      //        if(print_aflow_logo_full){node_option_ss << "shift={(-2cm,-0.5cm)}";} //1
      //        else{node_option_ss << "shift={(-2cm,-0.75cm)}";}   //2
      //      }else{
      //        if(print_aflow_logo_full){
      //          if(!small_banner) {node_option_ss << "shift={(-2cm,-1.0cm)}";} //3
      //          else{node_option_ss << "shift={(-2cm,0.75cm)}";} //4
      //        }else{node_option_ss << "shift={(-2cm,-0.75cm)}";}  //5
      //      }
      //    }else{node_option_ss << "shift={(-1.25cm,0.75cm)}";}  //6
      //    // get node position
      //    node_position_ss << "current bounding box.north east";
      //    // get node content
      //    if(!small_banner&&print_logo_2){
      //      if(external_links) {node_content_ss << "\\href{" << NOMAD_WEB << "}";}
      //      node_content_ss << "{\\includegraphics[scale=0.1]{" << logo_file_2 << "}}";
      //    }else{
      //      node_content_ss << "\\large";
      //      node_content_ss << "{\\fontfamily{phv}\\selectfont";
      //      if(external_links) {node_content_ss << "\\href{" << AFLOW_WEB << "}";}
      //      node_content_ss << "{" << AFLOWLIB_MATERIALS_SERVER << "}";
      //      node_content_ss << "}";
      //    }
      //    // create node
      //    tikzpic_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);
      //
      //    ////////////////////////////////////////////////////////////////////////
      //    // END AFLOW.org on top right
      //    ////////////////////////////////////////////////////////////////////////
      //
      //    if(!small_banner) {
      //      //////////////////////////////////////////////////////////////////////
      //      // START Count on bottom left
      //      //////////////////////////////////////////////////////////////////////
      //
      //      // get node option
      //      if(image_only) {node_option_ss << "shift={(2.5cm,-1.0cm)}";}  //
      //      else{node_option_ss << "shift={(1.5cm,-0.75cm)}";}  //
      //      // get node position
      //      node_position_ss << "current bounding box.south west";
      //      // get node content
      //      node_content_ss << "count=" << count_entries;
      //      // create node
      //      tikzpic_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);
      //
      //      //////////////////////////////////////////////////////////////////////
      //      // END Count on bottom left
      //      //////////////////////////////////////////////////////////////////////
      //
      //      //////////////////////////////////////////////////////////////////////
      //      // START Date on bottom right
      //      //////////////////////////////////////////////////////////////////////
      //
      //      // get node option
      //      if(image_only) {node_option_ss << "shift={(-4.25cm,0.375cm)}";} //
      //      else{node_option_ss << "shift={(-2.5cm,-0.75cm)}";} //
      //      // get node position
      //      node_position_ss << "current bounding box.south east";
      //      // get node content
      //      node_content_ss << "\\today~\\currenttime";
      //      // create node
      //      tikzpic_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);
      //
      //      //////////////////////////////////////////////////////////////////////
      //      // END Date on bottom right
      //      //////////////////////////////////////////////////////////////////////
      //    }
      //    tikzpic_TEX_ss << "\\end{scope}" << endl;
      //  
      //  }else{  // dimension==3
      //
      //    ////////////////////////////////////////////////////////////////////////
      //    // START Full banner on left (logo,count,date)
      //    ////////////////////////////////////////////////////////////////////////
      //
      //    tikzpic_TEX_ss << "\\begin{scope}[remember picture,overlay,on background layer]" << endl;
      //    if(!small_banner) {
      //      // get node option
      //      node_content_ss << "\\shortstack[l]{";
      //      if(include_color_bar) {
      //        if(print_aflow_logo_full) {
      //          if(image_only){node_option_ss << "shift={(3.75cm,-2.5cm)}";}
      //          else{node_option_ss << "shift={(3.25cm,-1.0cm)}";}
      //        }else{
      //          if(image_only){node_option_ss << "shift={(4cm,-2.25cm)}";}
      //          else{node_option_ss << "shift={(3.25cm,-1.25cm)}";}
      //        }
      //      }else{
      //        if(print_aflow_logo_full) {
      //          if(image_only){node_option_ss << "shift={(3.25cm,-2.25cm)}";}
      //          else{node_option_ss << "shift={(2.25cm,-1.5cm)}";}
      //        }else{
      //          if(image_only){node_option_ss << "shift={(4cm,-2.25cm)}";}
      //          else{node_option_ss << "shift={(2cm,-1.5cm)}";}
      //        }
      //      }
      //      // get node position
      //      node_position_ss << "current bounding box.north west";
      //      // get node content
      //      if(print_aflow_logo_full) {
      //        if(external_links) {node_content_ss << "\\href{" << AFLOW_WEB << "}";}
      //        node_content_ss << "{\\includegraphics[scale=0.25]{" << aflow_logo_full_file << "}}";
      //      }else{
      //        if(external_links) {node_content_ss << "\\href{" << AFLOW_WEB << "}";}
      //        node_content_ss << "{\\LARGE AFLOW V" << string(AFLOW_VERSION) << "}";
      //        node_content_ss << "\\\\";
      //        node_content_ss << "{\\large C. Oses and S. Curtarolo}";
      //      }
      //      // get node content
      //      node_content_ss << "\\\\";
      //      node_content_ss << "\\\\";
      //      node_content_ss << "\\large count=" << count_entries;
      //      node_content_ss << "\\\\";
      //      node_content_ss << "\\large \\today~\\currenttime";
      //      node_content_ss << "}";
      //      // create node
      //      tikzpic_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);
      //    }
      //    ////////////////////////////////////////////////////////////////////////
      //    // END Full banner on top left (logo,count,date)
      //    ////////////////////////////////////////////////////////////////////////
      //
      //    ////////////////////////////////////////////////////////////////////////
      //    // START AFLOW.org on top right
      //    // nomad logo too
      //    ////////////////////////////////////////////////////////////////////////
      //
      //    // get node option
      //    if(!small_banner&&print_logo_2){
      //      if(include_color_bar) {
      //        if(image_only){node_option_ss << "shift={(-2cm,-0.75cm)}";}
      //        else{node_option_ss << "shift={(-1.5cm,0.0cm)}";} //
      //      }else{
      //        if(image_only) {node_option_ss << "shift={(-2.5cm,-1.25cm)}";}
      //        else{node_option_ss << "shift={(-0.25cm,-0.5cm)}";}
      //      }
      //    }else{
      //      if(include_color_bar) {
      //        if(image_only){node_option_ss << "shift={(-2cm,-0.75cm)}";}
      //        else{node_option_ss << "shift={(-1.5cm,0.0cm)}";} //
      //      }else{
      //        if(image_only) {node_option_ss << "shift={(-2.5cm,-1.25cm)}";}
      //        else{node_option_ss << "shift={(-0.25cm,-0.5cm)}";}
      //      }
      //    }
      //    // get node position
      //    node_position_ss << "current bounding box.north east";
      //    // get node content
      //    //node_content_ss << "\\shortstack[l]{";
      //    node_content_ss << "{\\large";
      //    node_content_ss << "{\\fontfamily{phv}\\selectfont";
      //    if(external_links) {node_content_ss << "\\href{" << AFLOW_WEB << "}";}
      //    node_content_ss << "{" << AFLOWLIB_MATERIALS_SERVER << "}";
      //    node_content_ss << "}";
      //    node_content_ss << "}";
      //    //node_content_ss << "}";
      //    // create node
      //    tikzpic_TEX_ss << nodeCreator(node_option_ss, node_position_ss, node_content_ss);
      //
      //    ////////////////////////////////////////////////////////////////////////
      //    // END AFLOW.org on top right
      //    ////////////////////////////////////////////////////////////////////////
      //
      //    tikzpic_TEX_ss << "\\end{scope}" << endl;
      //  }
      //}
      tikzpic_TEX_ss << "\\end{tikzpicture}" << endl;
      if(image_only) {tikzpic_TEX_ss << "\\endpgfgraphicnamed" << endl;}
      else{
        tikzpic_TEX_ss << "\\vspace*{\\fill}" << endl;
        tikzpic_TEX_ss << "}" << endl;
        tikzpic_TEX_ss << "\\end{landscape}" << endl;
        //if(!no_doc) {tikzpic_TEX_ss << "\\restoregeometry" << endl;}
      }
      tikzpic_TEX_ss << "\\restoregeometry" << endl;
    }else{
      // make sure to add this in
      doc_header_TEX_ss << "\\usepackage{graphicx}  \%HEADER" << endl;
      // contains begin{document}, which needs to go in sooner
      doc_header_TEX_ss << _doc_header_TEX_ss.str();
      _doc_header_TEX_ss.str("");
    }

    // BIG MERGE
    main_TEX_ss << doc_header_TEX_ss.str();
    doc_header_TEX_ss.str("");

    if(!doc_only && (dimension == 2 || dimension == 3)) {
      main_TEX_ss << _tikzpic_settings_TEX_ss.str();
      _tikzpic_settings_TEX_ss.str("");
      main_TEX_ss << pseudo_preliminary_axes_TEX_ss.str();
      pseudo_preliminary_axes_TEX_ss.str("");
      main_TEX_ss << tikzpic_settings_TEX_ss.str();
      tikzpic_settings_TEX_ss.str("");
      if(show_heat_map && num_horizontal_planes > 1) {
        main_TEX_ss << heat_map_TEX_ss.str();
        heat_map_TEX_ss.str("");
      }
      main_TEX_ss << tikzpic_TEX_ss.str();
      tikzpic_TEX_ss.str("");
    }

    //////////////////////////////////////////////////////////////////////////////
    // END Image on first page
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // START Report
    //////////////////////////////////////////////////////////////////////////////

    if(!no_doc && !image_only) {
      //if(!doc_only) {main_TEX_ss << "\\newpage" << endl;}
      ////////////////////////////////////////////////////////////////////////////
      // START Column header setup
      ////////////////////////////////////////////////////////////////////////////

      main_TEX_ss << "\\setlength{\\tabcolsep}{" << LATEX_TABCOLSEP_STD << "pt}" << endl;
      main_TEX_ss << "\\setlength\\arrayrulewidth{" << LATEX_ARRAYRULEWIDTH_STD << "pt}" << endl;

      main_TEX_ss << "\\newgeometry{" << MARGIN_REPORT << "}" << endl; //new geometry introduces \\newpage automatically
      main_TEX_ss << "\\updateSnapshotTableHeaderOne" << endl;

      main_TEX_ss << "\\thispagestyle{reportPage1}" << endl;
      main_TEX_ss << "\\pagestyle{reportPage2}" << endl;

      //citation
      main_TEX_ss << " \%CITATION" << endl;
      main_TEX_ss << "\\noindent " << CHULL_CITE << "" << " \%CITATION" << endl;
      main_TEX_ss << " \%CITATION" << endl;
      main_TEX_ss << "\\noindent\\begin{quote}" << " \%CITATION" << endl;
      main_TEX_ss << CHULL_AUTHORS << ", ";
      main_TEX_ss << "\\textit{" << CHULL_TITLE << "}, ";
      main_TEX_ss << CHULL_JOURNAL_LATEX << "." << " \%CITATION" << endl;
      main_TEX_ss << "\\end{quote}" << " \%CITATION" << endl;

      ////////////////////////////////////////////////////////////////////////////
      // END Column header setup
      ////////////////////////////////////////////////////////////////////////////

      uint counter;
      bool putColumnHeader = true;  // only put it for the first reported stoich
      vector<vector<uint> > equilibrium_phases;
      vector<vector<ChullPoint> > equilibrium_phases_CP;
      vector<ChullPoint> dummyDCP;
      vector<string> equilibrium_phases_vs, _equilibrium_phases_vs;
      vector<uint> decomposition_phases;
      //xvector<double> decomposition_coefficients;
      vector<string> decompositionCoefPhase_vs;  // combined coef * phase
      string misc;
      uint i_phase;
      bool added_nary_tag=false;  //for unary, binary, ternary, etc.
      bool add_vspace=true; //ONLY once, first TERNARY is too close to column headers
      string pdftable_font_sizes;
      vector<string> chpoint_properties;
      bool print_scriterion=true,print_np1=true;
      stringstream scriterion_data_ss,np1_data_ss;
      //uint num_cols_scriterion=2;
      //uint num_cols_np1=3;
      uint precision_tmp=0;
      double tmp_roundoff_tol;
      if(compounds_column_report){pdftable_font_sizes="\\fontsize{4}{6}\\selectfont";}
      else{pdftable_font_sizes="\\fontsize{5}{7}\\selectfont";}

      ////////////////////////////////////////////////////////////////////////////
      // START Stoichiometry group loop
      ////////////////////////////////////////////////////////////////////////////

      for(uint fl_size_i_nary=m_naries.size(),i_nary=(fl_size_i_nary-1);i_nary<fl_size_i_nary;i_nary--){ //go backwards!
        added_nary_tag=false;
        for(uint i_alloy=0,fl_size_i_alloy=m_naries[i_nary].m_alloys.size();i_alloy<fl_size_i_alloy;i_alloy++){
          for(uint i=0,fl_size_i=m_naries[i_nary].m_alloys[i_alloy].m_coord_groups.size();i<fl_size_i;i++){
            i_coord_group=m_naries[i_nary].m_alloys[i_alloy].m_coord_groups[i];
            counter = 0;

            //////////////////////////////////////////////////////////////////////////
            // START Stoichiometry group points loop
            //////////////////////////////////////////////////////////////////////////

            for(uint j=0,fl_size_j=m_coord_groups[i_coord_group].m_points.size();j<fl_size_j;j++) {
              ////////////////////////////////////////////////////////////////////////
              // START Entries filter
              ////////////////////////////////////////////////////////////////////////

              const ChullPoint& point = m_points[m_coord_groups[i_coord_group].m_points[j]];
              if(!point.m_has_entry) {continue;}
              const aflowlib::_aflowlib_entry& entry = point.m_entry;

              ////////////////////////////////////////////////////////////////////////
              // END Entries filter
              ////////////////////////////////////////////////////////////////////////

              if(!counter) {
                //////////////////////////////////////////////////////////////////////
                // START Concentration label
                //////////////////////////////////////////////////////////////////////

                if(stoich_header_report) {
                  for(int j=point.m_elements_present.lrows;j<=point.m_elements_present.urows;j++) {
                    if(point.m_elements_present[j]==1){
                      _report_data_ss << m_velements[j] << "$_{" << aurostd::utype2string(point.s_coords[j],COEF_PRECISION) << "}$";
                    }
                  }
                }else{_report_data_ss << prettyPrintCompound(entry,gcd_vrt,true,latex_ft);}

                //////////////////////////////////////////////////////////////////////
                // START Concentration label
                //////////////////////////////////////////////////////////////////////

                //////////////////////////////////////////////////////////////////////
                // START Gathering info about equilibrium phases / decomposition
                // reaction
                //////////////////////////////////////////////////////////////////////

                if(m_coord_groups[i_coord_group].m_is_on_hull) {

                  const vector<vector<uint> >& equilibrium_phases=m_coord_groups[i_coord_group].m_equilibrium_phases;
                  //nice fast way to toggle equivalent_phases on/off
                  //if(include_equilibrium_phases) {equilibrium_phases = m_coord_groups[i_coord_group].m_equilibrium_phases;}
                  //else{equilibrium_phases.clear();}

                  if(include_equilibrium_phases && !equilibrium_phases.empty()) {
                    for(uint k=0,fl_size_k=equilibrium_phases.size();k<fl_size_k;k++) {
                      for(uint l=0,fl_size_l=equilibrium_phases[k].size();l<fl_size_l;l++) {
                        i_phase=artificialMap(equilibrium_phases[k][l]);
                        const ChullPoint& eq_phase=m_points[i_phase];
                        if(!eq_phase.m_has_entry) {
                          // we need to adjust for missing unaries
                          // look for coords
                          output_name=aurostd::joinWDelimiter(alloyToElements(eq_phase),"");  //unary, so "" delimiter doesn't play a role
                        }else{
                          const aflowlib::_aflowlib_entry& equation_entry = eq_phase.m_entry;

                          output_name=prettyPrintCompound(equation_entry,gcd_vrt,true,latex_ft);
                          // do not hyperlink current point (pointless)
                          if(!aurostd::identical(point.getStoichiometricCoords(), eq_phase.getStoichiometricCoords(), ZERO_TOL)) {
                            if(internal_links_withinreport) {
                              misc_ss << "\\hyperref[" << input << "_" << equation_entry.auid << "]{";
                              misc_ss << output_name;
                              misc_ss << "}";
                              output_name = misc_ss.str();
                              misc_ss.str("");
                            }
                          }
                        }
                        _equilibrium_phases_vs.push_back(output_name);
                      }
                      if(m_coord_groups[i_coord_group].getDim() > 6) {equilibrium_phases_vs.push_back(aurostd::joinWDelimiter(_equilibrium_phases_vs, " -- "));}  // that way, we don't run off the line
                      else{equilibrium_phases_vs.push_back(aurostd::joinWDelimiter(_equilibrium_phases_vs, "--"));}
                      _equilibrium_phases_vs.clear();
                    }
                    equilibrium_phases_TEX_ss << aurostd::joinWDelimiter(equilibrium_phases_vs, ", ", " and ", ", and ");
                    equilibrium_phases_vs.clear();
                    equilibrium_phases_CP.clear();
                    //get header
                    //misc = prettyPrintCompound(entry,gcd_vrt,true,latex_ft);
                    equilibrium_phases_header_TEX_ss << "vertex of facets: (" << m_coord_groups[i_coord_group].getDim() << "-phase~equilibria)";  //patched ugly looking discontinuity
                    //equilibrium_phases_header_TEX_ss << " with " << misc;
                    //equilibrium_phases_header_TEX_ss << ":";
                  }
                }else{
                  // decomposition equation
                  const vector<uint>& decomposition_phases = m_coord_groups[i_coord_group].m_decomp_phases;
                  const xvector<double>& decomposition_coefficients = m_coord_groups[i_coord_group].m_decomp_coefs;
                  if(!decomposition_phases.empty() && scalar_product(decomposition_coefficients, decomposition_coefficients) >= ZERO_TOL) {
                    if(decomposition_phases.size() != (uint)decomposition_coefficients.rows-1) {
                      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Size of decomposition phases != size of decomposition coefficients for coordgroup["+aurostd::utype2string(i)+"]");
                    }
                    // write out decomposition equation
                    reaction_chem_eq_TEX_ss << prettyPrintCompound(entry,frac_vrt,true,latex_ft);
                    reaction_chem_eq_TEX_ss << " $\\to$ ";
                    for(uint k=0,fl_size_k=decomposition_phases.size();k<fl_size_k;k++) {
                      if(abs(decomposition_coefficients[decomposition_coefficients.lrows+k+1]) < ZERO_TOL) {continue;}
                      i_phase=artificialMap(decomposition_phases[k]);
                      const ChullPoint& dc_phase=m_points[i_phase];
                      if(!dc_phase.m_has_entry) {
                        // we need to adjust for missing unaries
                        // look for coords
                        output_name=aurostd::joinWDelimiter(alloyToElements(dc_phase),"");  //unary, so "" delimiter doesn't play a role
                      }else{
                        const aflowlib::_aflowlib_entry& equation_entry = dc_phase.m_entry;
                        output_name=prettyPrintCompound(equation_entry,frac_vrt,true,latex_ft);
                        if(internal_links_withinreport) {
                          misc_ss << "\\hyperref[" << input << "_" << equation_entry.auid << "]{";
                          misc_ss << output_name;
                          misc_ss << "}";
                          output_name = misc_ss.str();
                          misc_ss.str("");
                        }
                      }
                      misc_ss << aurostd::utype2string(abs(decomposition_coefficients[decomposition_coefficients.lrows+k+1]),COEF_PRECISION) << "~";
                      misc_ss << output_name;
                      decompositionCoefPhase_vs.push_back(misc_ss.str());
                      misc_ss.str("");
                    }
                    reaction_chem_eq_TEX_ss << aurostd::joinWDelimiter(decompositionCoefPhase_vs, " + ");
                    decompositionCoefPhase_vs.clear();
                  }
                }

                //////////////////////////////////////////////////////////////////////
                // END Gathering info about equilibrium phases / decomposition
                // reaction
                //////////////////////////////////////////////////////////////////////
              }

              ////////////////////////////////////////////////////////////////////////
              // START Row properties set up
              ////////////////////////////////////////////////////////////////////////

              //sg_tokens=entry.vsg;
              //if(!entry.sg.empty()) {
              //  aurostd::string2tokens(entry.sg, sg_tokens, ",");
              //  if(sg_tokens.size() != 3) {
              //    continue;
              //  }
              //}
              if(point.isGState()) {report_data_ss << aurostd::PaddedPOST("\\rowcolor{green!85!blue} ", 30);} //red!25
              else if(point.m_is_sym_equivalent_g_state) {report_data_ss << aurostd::PaddedPOST("\\rowcolor{orange!85} ", 30);}  // odd should be white
              else if(counter % 2) {report_data_ss << aurostd::PaddedPOST("\\rowcolor{white} ", 30);}  // odd should be white
              else{report_data_ss << aurostd::PaddedPOST("\\rowcolor{gray!25} ", 30);}
              chpoint_properties.clear();
              for(uint i=0,fl_size_i=vheaders.size();i<fl_size_i;i++){chpoint_properties.push_back(grabCHPointProperty(point,vheaders[i],latex_ft));}
              report_data_ss << aurostd::joinWDelimiter(chpoint_properties," & ");
              report_data_ss << " \\\\" << endl;
              counter++;
            }

            //////////////////////////////////////////////////////////////////////////
            // END Entry properties output
            //////////////////////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////////////////////
            // START Entry table creation
            //////////////////////////////////////////////////////////////////////////

            if(!report_data_ss.str().empty()) {
              if(putColumnHeader){
                main_TEX_ss << "\\noindent{" << pdftable_font_sizes << endl;
                main_TEX_ss << getSnapshotTableHeader(headers,false);
                main_TEX_ss << "}" << endl;
              }
              if(!added_nary_tag){
                if(add_vspace){
                  main_TEX_ss << endl;
                  main_TEX_ss << "\\vspace{5pt}" << endl;
                }
                main_TEX_ss << "\\centering{" << general_image_font_size << "{" << pflow::arity_string(i_nary+1,false,true) << "}}" << endl; //lowercase per SC's request
                main_TEX_ss << endl;
                main_TEX_ss << "\\vspace{-15pt}" << endl;
                added_nary_tag=true;
              }
              //main_TEX_ss << "{" << pdftable_font_sizes << endl;
              ////[CO20190226 - TABU IS BROKEN IN TeX Live 2019]main_TEX_ss << "\\begin{longtabu}{|" << aurostd::joinWDelimiter(valignments,"|") << "|}" << endl;
              //main_TEX_ss << "\\begin{longtable}{|" << aurostd::joinWDelimiter(valignments_entrytable_string,"|") << "|}" << endl;
              //putColumnHeader = false;
              //if(internal_links && isViablePoint(m_coord_groups[i_coord_group].m_ref_state)) {
              //  uint ref_state=m_coord_groups[i_coord_group].m_ref_state;
              //  main_TEX_ss << "\\multicolumn{" << vheaders.size() << "}{l}{\\phantomsection\\label{"+input + "_" + m_points[ref_state].m_entry.auid + "}} \\\\[0.1cm]" << endl;
              //}else{main_TEX_ss << "\\multicolumn{" << vheaders.size() << "}{l}{} \\\\[0.1cm]" << endl;}  //padding and spacing preservation
              print_scriterion=false;
              print_np1=false;
              scriterion_data_ss.str("");
              np1_data_ss.str("");
              if(m_coord_groups[i_coord_group].m_is_on_hull) {
                _report_data_ss << " " << "(ground-state)"; // if ground-state
                if(m_coord_groups[i_coord_group].m_stability_criterion<AUROSTD_NAN){
                  print_scriterion=true;
                  precision_tmp=0;
                  tmp_roundoff_tol=5.0*pow(10,-((int)precision_tmp)-1);
                  scriterion_data_ss << "$" << getStabilityCriterionSymbol(helvetica_font) << "="; //this delta is okay, should be italicized
                  scriterion_data_ss << aurostd::utype2string(convertUnits(m_coord_groups[i_coord_group].m_stability_criterion,(m_formation_energy_hull?_m_:_std_)),precision_tmp,true,tmp_roundoff_tol,FIXED_STREAM);
                  scriterion_data_ss << "$~" << (m_formation_energy_hull?string("meV/atom"):string("K"));
                }
                if(m_coord_groups[i_coord_group].m_i_nary>=1 && m_coord_groups[i_coord_group].m_n_plus_1_enthalpy_gain<AUROSTD_NAN){ //print binaries anyway... //print only for ternaries and up, binaries is trivial formation enthalpy (save space)
                  print_np1=true;
                  precision_tmp=0;
                  tmp_roundoff_tol=5.0*pow(10,-((int)precision_tmp)-1);
                  if(0){
                    np1_data_ss << "$" << getDelta(helvetica_font) << " H[N|\\{1,\\cdots,N-1\\}]="; //this delta is okay, should be italicized
                  }else{
                    np1_data_ss << "$" << getDelta(helvetica_font) << " H[" << m_coord_groups[i_coord_group].m_i_nary+1 << "|"; //this delta is okay, should be italicized
                    if(m_coord_groups[i_coord_group].m_i_nary==0){np1_data_ss << "0";}  //not allowed anymore, but we can print the value
                    else if(m_coord_groups[i_coord_group].m_i_nary==1){np1_data_ss << "1";}
                    else if(m_coord_groups[i_coord_group].m_i_nary==2){np1_data_ss << "\\{1,2\\}";}
                    else if(m_coord_groups[i_coord_group].m_i_nary==3){np1_data_ss << "\\{1,2,3\\}";}
                    else if(m_coord_groups[i_coord_group].m_i_nary>3){np1_data_ss << "\\{1,\\cdots," << m_coord_groups[i_coord_group].m_i_nary << "\\}";}
                    np1_data_ss << "]=";
                  }
                  np1_data_ss << aurostd::utype2string(convertUnits(m_coord_groups[i_coord_group].m_n_plus_1_enthalpy_gain,(m_formation_energy_hull?_m_:_std_)),precision_tmp,true,tmp_roundoff_tol,FIXED_STREAM);
                  np1_data_ss << "$~" << (m_formation_energy_hull?string("meV/atom"):string("K"));
                }
              }else{
                _report_data_ss << " " << "(unstable)"; // if above hull	  
              }
              // compound name
              main_TEX_ss << "\\begin{longtable}{" << aurostd::joinWDelimiter( (m_coord_groups[i_coord_group].m_is_on_hull?valignments_compoundname_thermopropstable_string:valignments_onecol_lefttable_string) ,"") << "}" << endl;
              if(internal_links && isViablePoint(m_coord_groups[i_coord_group].m_ref_state)) {
                uint ref_state=m_coord_groups[i_coord_group].m_ref_state;
                main_TEX_ss << "\\multicolumn{" << (m_coord_groups[i_coord_group].m_is_on_hull?n_cols_compoundname_thermoprops:n_cols_onecol_left) << "}{l}{\\phantomsection\\label{"+input + "_" + m_points[ref_state].m_entry.auid + "}} \\\\[0.1cm]" << endl;
              }else{main_TEX_ss << "\\multicolumn{" << (m_coord_groups[i_coord_group].m_is_on_hull?n_cols_compoundname_thermoprops:n_cols_onecol_left) << "}{l}{} \\\\[0.1cm]" << endl;}  //padding and spacing preservation
              //main_TEX_ss << "\\multicolumn{" << vheaders.size()-((print_scriterion?num_cols_scriterion:0)+(print_np1?num_cols_np1:0)) << "}{l}{";
              main_TEX_ss << "\\cellcolor{white}\\normalsize{" + _report_data_ss.str() << "}";
              //main_TEX_ss << "}";
              if(print_scriterion || print_np1){
                main_TEX_ss << " & ";
                //main_TEX_ss << "\\multicolumn{" << (print_scriterion?num_cols_scriterion:0)+(print_np1?num_cols_np1:0) << "}{r}{";
                main_TEX_ss << "\\cellcolor{white}\\normalsize{" + (print_scriterion?scriterion_data_ss.str():"") + (print_scriterion&&print_np1?", ":"") + (print_np1?np1_data_ss.str():"") << "}";
                //main_TEX_ss << "}";
              }
              main_TEX_ss << " \\\\[0.05cm]" << endl;  //0.1cm
              main_TEX_ss << "\\end{longtable}" << endl;
              main_TEX_ss << "\\vspace{-23pt}" << endl;

              //main table
              main_TEX_ss << "{" << pdftable_font_sizes << endl;
              //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]main_TEX_ss << "\\begin{longtabu}{|" << aurostd::joinWDelimiter(valignments,"|") << "|}" << endl;
              main_TEX_ss << "\\begin{longtable}{|" << aurostd::joinWDelimiter(valignments_entrytable_string,"|") << "|}" << endl;
              putColumnHeader = false;
              //if(internal_links && isViablePoint(m_coord_groups[i_coord_group].m_ref_state)) {
              //  uint ref_state=m_coord_groups[i_coord_group].m_ref_state;
              //  main_TEX_ss << "\\multicolumn{" << vheaders.size() << "}{l}{\\phantomsection\\label{"+input + "_" + m_points[ref_state].m_entry.auid + "}} \\\\[0.1cm]" << endl;
              //}else{main_TEX_ss << "\\multicolumn{" << vheaders.size() << "}{l}{} \\\\[0.1cm]" << endl;}  //padding and spacing preservation
              main_TEX_ss << "\\hline" << endl;
              main_TEX_ss << report_data_ss.str();
              main_TEX_ss << "\\hline" << endl;
              main_TEX_ss << "\\end{longtable}" << endl;
              //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]main_TEX_ss << "\\end{longtabu}" << endl;
              main_TEX_ss << "}" << endl;

              ////////////////////////////////////////////////////////////////////////
              // END Entry table creation
              ////////////////////////////////////////////////////////////////////////

              ////////////////////////////////////////////////////////////////////////
              // START Equilibrium phases / decomposition reaction table
              ////////////////////////////////////////////////////////////////////////

              if(!(equilibrium_phases_header_TEX_ss.str().empty() && equilibrium_phases_TEX_ss.str().empty() && reaction_chem_eq_TEX_ss.str().empty())) {
                main_TEX_ss << "\\vspace{-20pt}" << endl;
                //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]main_TEX_ss << "\\begin{longtabu}{X[1,l]X[3,r]}" << endl;
                main_TEX_ss << "\\begin{longtable}{" << aurostd::joinWDelimiter(valignments_reactiontable_string,"") << "}" << endl;
                if(m_coord_groups[i_coord_group].m_is_on_hull) {
                  if(include_equilibrium_phases && !(equilibrium_phases_header_TEX_ss.str().empty() && equilibrium_phases_TEX_ss.str().empty())) {
                    main_TEX_ss << equilibrium_phases_header_TEX_ss.str();
                    main_TEX_ss << " & ";
                    main_TEX_ss << equilibrium_phases_TEX_ss.str() << endl;
                  }
                }else{
                  if(!reaction_chem_eq_TEX_ss.str().empty()) {
                    main_TEX_ss << "decomposition reaction: & ";
                    main_TEX_ss << reaction_chem_eq_TEX_ss.str() << endl;
                  }
                }
                main_TEX_ss << "\\end{longtable}" << endl;
                //[CO20190226 - TABU IS BROKEN IN TeX Live 2019]main_TEX_ss << "\\end{longtabu}" << endl;
              }

              ////////////////////////////////////////////////////////////////////////
              // END Equilibrium phases / decomposition reaction table
              ////////////////////////////////////////////////////////////////////////

              equilibrium_phases_header_TEX_ss.str("");
              equilibrium_phases_TEX_ss.str("");
              reaction_chem_eq_TEX_ss.str("");
              if(!((i_alloy==m_naries[i_nary].m_alloys.size()-1)&&(i==m_naries[i_nary].m_alloys[i_alloy].m_coord_groups.size()-1))) {
                main_TEX_ss << "\\vspace{-20pt}" << endl;
              }
              //if(i!=m_coord_groups.size()-1) {main_TEX_ss << "\\vspace{-20pt}" << endl;}
            }
            scriterion_data_ss.str("");
            np1_data_ss.str("");
            _report_data_ss.str("");
            report_data_ss.str("");
            reaction_chem_eq_TEX_ss.str("");
          }

          ////////////////////////////////////////////////////////////////////////////
          // END Stoichiometry group points loop
          ////////////////////////////////////////////////////////////////////////////
        }
      }
      main_TEX_ss << "\\restoregeometry" << endl;
    }
    //////////////////////////////////////////////////////////////////////////////
    // END Stoichiometry group loop
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // END Report
    //////////////////////////////////////////////////////////////////////////////

    main_TEX_ss << "\\end{document}" << " \%FOOTER" << endl;

    //////////////////////////////////////////////////////////////////////////////
    // START Create tmp directory for compilation of .tex document
    //////////////////////////////////////////////////////////////////////////////

    string PWD = getPath();
    string path = getPath(m_cflags, *p_FileMESSAGE, *p_oss);
    string destination = path+main_PDF_file;
    main_output_file=main_PDF_file;
    if(m_cflags.flag("CHULL::PNG_IMAGE")&&!m_cflags.flag("CHULL::LATEX_DOC")){destination=path+main_PNG_file;main_output_file=main_PNG_file;}
    string LATEX_dir = aurostd::TmpDirectoryCreate("chullLATEX");
    //[CO20221027 - suppressing warning]if(0) chdir(LATEX_dir.c_str());
    aurostd::stringstream2file(main_TEX_ss, LATEX_dir+"/"+main_TEX_file);
    if(!aurostd::FileExist(LATEX_dir+"/"+main_TEX_file)) {
      //[CO20221027 - suppressing warning]if(0) chdir(PWD.c_str());
#ifndef _AFLOW_TEMP_PRESERVE_
      aurostd::RemoveDirectory(LATEX_dir);
#endif
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Could not write "+main_TEX_file+" to "+LATEX_dir);
    }
    //watermark
    if(!doc_only){
      aurostd::EmbData::save_to_file("AFLOW_LOGO_SKINNY.pdf", LATEX_dir+"/"+aflow_logo_skinny_file);
      if(!aurostd::FileExist(LATEX_dir+"/"+aflow_logo_skinny_file)) {
        //[CO20221027 - suppressing warning]if(0) chdir(PWD.c_str());
#ifndef _AFLOW_TEMP_PRESERVE_
        aurostd::RemoveDirectory(LATEX_dir);
#endif
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Could not write "+aflow_logo_skinny_file+" to "+LATEX_dir);
      }
    }
    if(print_aflow_logo_full) {
      aurostd::EmbData::save_to_file("AFLOW_LOGO_FULL.pdf", LATEX_dir+"/"+aflow_logo_full_file);
      if(!aurostd::FileExist(LATEX_dir+"/"+aflow_logo_full_file)) {
        //[CO20221027 - suppressing warning]if(0) chdir(PWD.c_str());
#ifndef _AFLOW_TEMP_PRESERVE_
        aurostd::RemoveDirectory(LATEX_dir);
#endif
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Could not write "+aflow_logo_full_file+" to "+LATEX_dir);
      }
    }
    if(print_logo_2) {
      aurostd::EmbData::save_to_file("NOMAD_LOGO.png", LATEX_dir+"/"+logo_file_2);
      if(!aurostd::FileExist(LATEX_dir+"/"+logo_file_2)) {
        //[CO20221027 - suppressing warning]if(0) chdir(PWD.c_str());
#ifndef _AFLOW_TEMP_PRESERVE_
        aurostd::RemoveDirectory(LATEX_dir);
#endif
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Could not write "+logo_file_2+" to "+LATEX_dir);
      }
    }
    command.str("");
    stringstream clean_command;
    uint num_compile;
    command << "cd " << LATEX_dir << " && ";
    command << XHOST.command("pdflatex") << " ";
    clean_command << XHOST.command("pdflatex") << " ";
    if(!show_latex_output) {command << "-interaction=nonstopmode -halt-on-error ";}  //-interaction=batchmode
    if(image_only) {
      command << "--jobname=" << main_file << " " << main_TEX_file << " ";
      clean_command << "--jobname=" << main_file << " " << main_TEX_file << " ";
      num_compile = 1;
    }else{
      command << main_TEX_file << " ";
      clean_command << main_TEX_file << " ";

      //if(doc_only && no internal_links) { //[CO20200106 - close bracket for indenting]}
      if(!internal_links) {num_compile = 1;}
      else{num_compile = 2;}
    }
    if(!show_latex_output) {command << "1>/dev/null ";}
    message << "Attempting to compile " << main_TEX_file;  //CO20180220 //the .tex file";
    pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

    if(show_latex_output) {
      if(latex_interactive_mode) {for (uint i = 0; i < num_compile; i++) {aurostd::execute(command.str());}} // will not save output, allows you to interact with LaTEX
      else{for (uint i = 0; i < num_compile; i++) {*p_oss << aurostd::execute2string(command.str()) << endl;}} // saves output
    }else{for (uint i = 0; i < num_compile; i++) {aurostd::execute(command.str());}} // no output to save
    if(!aurostd::FileExist(LATEX_dir+"/"+main_PDF_file)) {
      message << main_PDF_file << " was not created successfully, likely a LaTeX issue";
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_ERROR_);
      files_2_move.clear(); //only move these files
      files_2_move.push_back(LATEX_dir+"/"+main_TEX_file);
      if(!doc_only){files_2_move.push_back(LATEX_dir+"/"+aflow_logo_skinny_file);}
      if(print_aflow_logo_full){files_2_move.push_back(LATEX_dir+"/"+aflow_logo_full_file);}
      if(print_logo_2){files_2_move.push_back(LATEX_dir+"/"+logo_file_2);}

      message << "Moving " << aurostd::joinWDelimiter(files_2_move,", "," and ",", and ") << " to " << path; //CO20180220 - current directory";
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
      message << "Try running \"" << aurostd::RemoveWhiteSpacesFromTheBack(clean_command.str()) << "\"";
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

      aurostd::file2directory(files_2_move, path);
      //[CO20221027 - suppressing warning]if(0) chdir(PWD.c_str());
#ifndef _AFLOW_TEMP_PRESERVE_
      aurostd::RemoveDirectory(LATEX_dir);
#endif
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Issues compiling .tex file");
    }
    if(m_cflags.flag("CHULL::PNG_IMAGE")){
      command.str("");clean_command.str("");
      uint default_resolution=DEFAULT_CHULL_PNG_RESOLUTION,resolution=default_resolution;
      string resolution_input=m_cflags.getattachedscheme("CHULL::PNG_RESOLUTION");
      if(!resolution_input.empty()){
        if(!aurostd::isfloat(resolution_input)){
          message << "PNG_RESOLUTION input is not a number, defaulting to " << default_resolution;
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        }else{resolution=aurostd::string2utype<int>(resolution_input);}
        if(resolution==0){
          message << "PNG_RESOLUTION input is 0, defaulting to " << default_resolution;
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
          resolution=default_resolution;
        }
      }
      message << "Attempting to convert " << main_PDF_file << " to " << main_PNG_file;  //CO20180220 //the .tex file";
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
      command << "cd " << LATEX_dir << " && ";
      command << XHOST.command("convert") << " -density " << resolution << " " << main_PDF_file << " " << main_PNG_file;
      clean_command << XHOST.command("convert") << " -density " << resolution << " " << main_PDF_file << " " << main_PNG_file;
      command << " 1>/dev/null 2>&1";
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " executing: " << clean_command.str() << endl;}
      string convert_output=aurostd::execute2string(command.str());
      if(!aurostd::RemoveWhiteSpaces(convert_output).empty()){
        message << main_PNG_file << " was not created successfully, likely a convert issue";
        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_ERROR_);
        files_2_move.clear(); //only move these files
        files_2_move.push_back(LATEX_dir+"/"+main_PDF_file);

        message << "Moving " << aurostd::joinWDelimiter(files_2_move,", "," and ",", and ") << " to " << path; //CO20180220 - current directory";
        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
        message << "Try running \"" << aurostd::RemoveWhiteSpacesFromTheBack(clean_command.str()) << "\"";
        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

        aurostd::file2directory(files_2_move, path);
        //[CO20221027 - suppressing warning]if(0) chdir(PWD.c_str());
#ifndef _AFLOW_TEMP_PRESERVE_
        aurostd::RemoveDirectory(LATEX_dir);
#endif
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Issues converting PDF to PNG file");
      }
    }
    if(m_cflags.flag("CHULL::LATEX_DOC")){files_2_move.push_back(LATEX_dir+"/"+main_PDF_file);}
    if(m_cflags.flag("CHULL::PNG_IMAGE")){files_2_move.push_back(LATEX_dir+"/"+main_PNG_file);}
    if(keep_tex) {
      files_2_move.push_back(LATEX_dir+"/"+main_TEX_file);
      if(!doc_only){files_2_move.push_back(LATEX_dir+"/"+aflow_logo_skinny_file);}
      if(print_aflow_logo_full){files_2_move.push_back(LATEX_dir+"/"+aflow_logo_full_file);} //files_2_move.push_back(LATEX_dir+"/"+aflow_logo_skinny_file);
      if(print_logo_2){files_2_move.push_back(LATEX_dir+"/"+logo_file_2);}
      message << "Moving " << aurostd::joinWDelimiter(files_2_move,", "," and ",", and ") << " to " << path; //CO20180220 - current directory";
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
    }
    if(!aurostd::file2directory(files_2_move, path)) {
      //[CO20221027 - suppressing warning]if(0) chdir(PWD.c_str());
#ifndef _AFLOW_TEMP_PRESERVE_
      aurostd::RemoveDirectory(LATEX_dir);
#endif
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unable to move files out of temporary compilation directory");
    }
    //[CO20221027 - suppressing warning]if(0) chdir(PWD.c_str());
#ifndef _AFLOW_TEMP_PRESERVE_
    aurostd::RemoveDirectory(LATEX_dir);
#endif
    if(!aurostd::FileExist(destination)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unable to write "+main_output_file+" to "+path);}
    message << main_output_file << " was created successfully, see destination=" << path;
    pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_COMPLETE_);

    //////////////////////////////////////////////////////////////////////////////
    // END Create tmp directory for compilation of .tex document
    //////////////////////////////////////////////////////////////////////////////
  }

  string ConvexHull::getPlainTextHeader() const {
    stringstream main_text_ss;main_text_ss.str("");
    stringstream misc_ss;misc_ss.str("");
    const uint width_page_txt = 130;   // helps with formatting of output
    //misc_ss << aurostd::joinWDelimiter(m_velements,"") << " summary"; //Materials Snapshot
    //main_text_ss << aurostd::PaddedCENTER(misc_ss.str(), width_page_txt+2); misc_ss.str("");  //CENTER is off from PRE/POST by 2 spaces (see where i starts)
    //main_text_ss << endl;
    if(0){
      misc_ss << aurostd::joinWDelimiter(m_velements,"-") << " summary"; // (V" << string(AFLOW_VERSION) << ")";
      main_text_ss << aurostd::PaddedPOST(misc_ss.str(), width_page_txt / 2); misc_ss.str("");
      misc_ss << AFLOW_WEB;
      main_text_ss << aurostd::PaddedPRE(misc_ss.str(), width_page_txt / 2); misc_ss.str("");
      main_text_ss << endl;

      //misc_ss << "count=" << getEntriesCount(m_half_hull) << "/" << getEntriesCount(false);
      misc_ss << "count_entries_total=" << getEntriesCount(false);
      main_text_ss << aurostd::PaddedPOST(misc_ss.str(), width_page_txt / 2); misc_ss.str("");
      main_text_ss << aurostd::PaddedPRE(aurostd::get_datetime_formatted("-",false), width_page_txt / 2); //CO20180819 - we don't need minute/second precision, date is fine
      main_text_ss << endl;
    }

    main_text_ss << aurostd::joinWDelimiter(m_velements,"-") << " summary" << " (" << AFLOWLIB_MATERIALS_SERVER << ")" << endl;
    main_text_ss << aurostd::get_datetime_formatted("-",false) << endl;

    //cite as
    main_text_ss << endl;
    main_text_ss << CHULL_CITE << endl;
    main_text_ss << "    " << CHULL_AUTHORS << "," << endl;
    main_text_ss << "    " << CHULL_TITLE << "," << endl;
    main_text_ss << "    " << CHULL_JOURNAL_TXT << "." << endl;
    main_text_ss << endl;

    main_text_ss << endl;
    main_text_ss << "HULL DATA" << endl;
    //gstate_counts
    if(m_has_stoich_coords){
      uint gstate_count=0,total_gstate_count=0;
      for(uint i_nary=0,fl_size_i_nary=m_naries.size();i_nary<fl_size_i_nary;i_nary++){
        gstate_count=getGStateCount(i_nary);total_gstate_count+=gstate_count;
        main_text_ss << "Count ground-state " << pflow::arity_string(i_nary+1,false,true) << " = " << gstate_count << endl;
      }
      main_text_ss << "Count ground-state total = " << total_gstate_count;
      main_text_ss << endl;
    }
    main_text_ss << "Count entries total = " << getEntriesCount(false) << endl;
    main_text_ss << endl;

    return main_text_ss.str();
  }

  string ConvexHull::getJSONHeader() const {
    vector<string> vout;
    stringstream misc_ss;misc_ss.str("");

    misc_ss << "\"alloy\":\"" << aurostd::joinWDelimiter(m_velements,"") << "\"";
    vout.push_back(misc_ss.str()); misc_ss.str("");

    misc_ss << "\"aflow_version\":\"aflow" << string(AFLOW_VERSION) << "\"";
    vout.push_back(misc_ss.str()); misc_ss.str("");

    misc_ss << "\"aflow_website\":\"" << AFLOW_WEB << "\"";
    vout.push_back(misc_ss.str()); misc_ss.str("");

    misc_ss << "\"date\":\"" << aurostd::get_datetime_formatted("-",false) << "\""; //CO20180819 - we don't need minute/second precision, date is fine //previously datetime
    vout.push_back(misc_ss.str()); misc_ss.str("");

    //misc_ss << "\"count_hull\":\"" << getEntriesCount(m_half_hull) << "\"";
    //vout.push_back(misc_ss.str()); misc_ss.str("");
    //misc_ss << "\"count_total\":\"" << getEntriesCount(false) << "\"";

    if(m_has_stoich_coords){
      uint gstate_count=0,total_gstate_count=0;
      for(uint i_nary=0,fl_size_i_nary=m_naries.size();i_nary<fl_size_i_nary;i_nary++){
        gstate_count=getGStateCount(i_nary);total_gstate_count+=gstate_count;
        misc_ss << "\"count_gstate_" << pflow::arity_string(i_nary+1,false,false) << "\":" << gstate_count;
        vout.push_back(misc_ss.str()); misc_ss.str("");
      }
      misc_ss << "\"count_gstate_total\":" << total_gstate_count;
      vout.push_back(misc_ss.str()); misc_ss.str("");
    }

    misc_ss << "\"count_entries_total\":" << getEntriesCount(false);
    vout.push_back(misc_ss.str()); misc_ss.str("");

    misc_ss << "\"publication\":\"";
    misc_ss << CHULL_AUTHORS << ", ";
    misc_ss << CHULL_TITLE << ", ";
    misc_ss << CHULL_JOURNAL_TXT;
    misc_ss << "\"";
    vout.push_back(misc_ss.str()); misc_ss.str("");

    return aurostd::wrapString(aurostd::joinWDelimiter(vout,","),"{","}");
  }

  string ConvexHull::grabCHPointProperty(const ChullPoint& point,const string& property,filetype ftype) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " start" << endl;}
    uint precision=COEF_PRECISION,precision_tmp=COEF_PRECISION;
    double tmp_roundoff_tol=5.0*pow(10,-((int)precision)-1);
    double d_tmp=0.0;
    string value="";
    string equilibrium_phases_delimiter="-";
    string string_wrapper="";
    string list_prefix="";
    string list_suffix="";
    string null_value="-";
    uint padding=30;  //latex only
    if(ftype==json_ft){
      equilibrium_phases_delimiter=",";
      string_wrapper="\"";
      list_prefix="[";
      list_suffix="]";
      null_value="null";
    }
    string string_wrapper_math_mode="$";  //latex only
    bool compounds_column_report=false;
    bool external_links=false;
    if(ftype==latex_ft){
      null_value="N/A";
      compounds_column_report=DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN; //only grab if necessary, not an inexpensive string search
      external_links=addExternalHyperlinks(); //only grab if necessary, not an inexpensive string search
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " starting property=" << property << endl;}
    bool latex_property=false;
    bool math_mode=false;
    //[OBSOLETE]bool reduced=false;
    if(property=="compound"||property=="compound_latex"||property=="compound_reduced"||property=="compound_reduced_latex"||property=="compound_fractional"||property=="compound_fractional_latex"){
      math_mode=false;  //do not wrap these, they are automatically wrapped inside prettyPrint
      latex_property=(property=="compound_latex"||property=="compound_reduced_latex"||property=="compound_fractional_latex");
      vector_reduction_type vred=no_vrt;
      if(property=="compound_reduced"||property=="compound_reduced_latex"){vred=gcd_vrt;}
      if(property=="compound_fractional"||property=="compound_fractional_latex"){vred=frac_vrt;}
      //[OBSOLETE]reduced=(property=="compound_reduced"||property=="compound_reduced_latex");
      value=prettyPrintCompound(point,vred,(ftype==latex_ft||latex_property||vred==gcd_vrt||vred==frac_vrt),(latex_property?latex_ft:ftype));
      if((ftype==latex_ft||latex_property)&&math_mode&&!value.empty()){value=aurostd::wrapString(value,string_wrapper_math_mode);}
      value=aurostd::wrapString(value,string_wrapper);
    }
    else if(property=="prototype"||property=="prototype_latex"){
      latex_property=(property=="prototype_latex");
      if(ftype==latex_ft||latex_property){value=prettyPrintPrototype(point,ftype==json_ft,false);}  //JSON needs double backslash
      else{value=point.m_entry.prototype;}
      value=aurostd::wrapString(value,string_wrapper);
      padding=100;  //latex only
    }
    else if(property=="auid"){
      if(ftype==latex_ft){value="\\texttt{"+aurostd::wrapString(point.m_entry.auid,string_wrapper)+"}";}
      else{value=aurostd::wrapString(point.m_entry.auid,string_wrapper);}
    }
    else if(property=="aurl"){value=aurostd::wrapString(point.m_entry.aurl,string_wrapper);}
    else if(property=="url_entry_page"){value=aurostd::wrapString(ENTRY_PAGE_URL_PREFIX+point.m_entry.auid,string_wrapper);}
    else if(property=="space_group_orig"||property=="space_group_orig_latex"){
      math_mode=true;
      latex_property=(property=="space_group_orig_latex");
      string sg=point.getVSG().front();
      if(ftype==latex_ft||latex_property){sg=aurostd::fixStringLatex(sg,ftype==json_ft,true);}
      value=sg;
      if((ftype==latex_ft||latex_property)&&math_mode&&!value.empty()){value=aurostd::wrapString(value,string_wrapper_math_mode);}
      value=aurostd::wrapString(value,string_wrapper);
    }
    else if(property=="space_group_relax"||property=="space_group_relax_latex"){
      math_mode=true;
      latex_property=(property=="space_group_relax_latex");
      string sg=point.getVSG().back();
      if(ftype==latex_ft||latex_property){sg=aurostd::fixStringLatex(sg,ftype==json_ft,true);}
      value=sg;
      if((ftype==latex_ft||latex_property)&&math_mode&&!value.empty()){value=aurostd::wrapString(value,string_wrapper_math_mode);}
      value=aurostd::wrapString(value,string_wrapper);
    }
    else if(property=="spin_atom"){
      d_tmp=point.m_entry.spin_atom;
      if(d_tmp!=AUROSTD_NAN){
        precision_tmp=precision;
        if(ftype==latex_ft){precision_tmp=2;tmp_roundoff_tol=5.0*pow(10,-((int)precision_tmp)-1);}
        value=aurostd::utype2string(d_tmp,precision_tmp,true,tmp_roundoff_tol,FIXED_STREAM);
      }
    }
    else if(property=="enthalpy_formation_atom"){
      precision_tmp=precision;
      if(ftype==latex_ft){precision_tmp=0;tmp_roundoff_tol=5.0*pow(10,-((int)precision_tmp)-1);}
      value=aurostd::utype2string(H_f_atom(point,_m_),precision_tmp,true,tmp_roundoff_tol,FIXED_STREAM);
    }
    else if(property=="entropic_temperature"){
      precision_tmp=precision;
      if(ftype==latex_ft){precision_tmp=0;tmp_roundoff_tol=5.0*pow(10,-((int)precision_tmp)-1);}
      value=aurostd::utype2string(point.m_entry.entropic_temperature,precision_tmp,true,tmp_roundoff_tol,FIXED_STREAM);
    }
    else if(property=="ground_state"){value=(point.isGState()?"true":"false");}
    else if(property=="auid_equivalent_structures"){
      if(!(ftype==txt_ft || ftype==json_ft)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No latex rule defined for "+property);}
      if(!point.m_equivalent_entries.empty()){
        const vector<uint>& equivalent_g_states=point.m_equivalent_entries;
        vector<string> auids;
        uint i_point=AUROSTD_MAX_UINT;
        for(uint i=0,fl_size_i=equivalent_g_states.size();i<fl_size_i;i++){
          i_point=artificialMap(equivalent_g_states[i]);
          if(m_points[i_point].m_has_entry){
            auids.push_back(aurostd::wrapString(m_points[i_point].m_entry.auid,string_wrapper));
          }else{auids.push_back(null_value);}
        }
        value=aurostd::wrapString(aurostd::joinWDelimiter(auids,","),list_prefix,list_suffix);
      }
      //[CO20221112 - OBSOLETE]if(point.isGState()){
      //[CO20221112 - OBSOLETE]  //need to grab from coord_group
      //[CO20221112 - OBSOLETE]  uint i_coord_group=AUROSTD_MAX_UINT;
      //[CO20221112 - OBSOLETE]  if(!getCoordGroupIndex(point,i_coord_group)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup index not set");}
      //[CO20221112 - OBSOLETE]  if(m_coord_groups[i_coord_group].m_is_on_hull){
      //[CO20221112 - OBSOLETE]    vector<string> auids;
      //[CO20221112 - OBSOLETE]    uint i_point=AUROSTD_MAX_UINT;
      //[CO20221112 - OBSOLETE]    const vector<uint>& equivalent_g_states=m_coord_groups[i_coord_group].m_equivalent_g_states;
      //[CO20221112 - OBSOLETE]    if(equivalent_g_states.size()){
      //[CO20221112 - OBSOLETE]      for(uint i=0,fl_size_i=equivalent_g_states.size();i<fl_size_i;i++){
      //[CO20221112 - OBSOLETE]        i_point=artificialMap(equivalent_g_states[i]);
      //[CO20221112 - OBSOLETE]        if(m_points[i_point].m_has_entry){
      //[CO20221112 - OBSOLETE]          auids.push_back(aurostd::wrapString(m_points[i_point].m_entry.auid,string_wrapper));
      //[CO20221112 - OBSOLETE]        }else{auids.push_back(null_value);}
      //[CO20221112 - OBSOLETE]      }
      //[CO20221112 - OBSOLETE]      value=aurostd::wrapString(aurostd::joinWDelimiter(auids,","),list_prefix,list_suffix);
      //[CO20221112 - OBSOLETE]    }
      //[CO20221112 - OBSOLETE]  }
      //[CO20221112 - OBSOLETE]}
    }
    else if(property=="ground_state_icsd"){
      bool icsd_g_state=false;
      if(point.isGState()){
        //need to grab from coord_group
        uint i_coord_group=AUROSTD_MAX_UINT;
        if(!getCoordGroupIndex(point,i_coord_group)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup index not set");}
        if(m_coord_groups[i_coord_group].m_is_on_hull){
          icsd_g_state=m_coord_groups[i_coord_group].m_found_icsd_g_state;
        }
      }
      value=(icsd_g_state?"true":"false");
    }
    else if(property=="auid_icsd"){
      if(point.m_found_icsd){
        uint i_point=point.m_i_icsd;
        if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
        value=aurostd::wrapString(m_points[i_point].m_entry.auid,string_wrapper);
      }
      //[CO20221112 - OBSOLETE]if(point.isGState()){
      //[CO20221112 - OBSOLETE]  //need to grab from coord_group
      //[CO20221112 - OBSOLETE]  uint i_coord_group=AUROSTD_MAX_UINT;
      //[CO20221112 - OBSOLETE]  if(!getCoordGroupIndex(point,i_coord_group)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup index not set");}
      //[CO20221112 - OBSOLETE]  if(m_coord_groups[i_coord_group].m_is_on_hull){
      //[CO20221112 - OBSOLETE]    if(m_coord_groups[i_coord_group].m_found_icsd_g_state){
      //[CO20221112 - OBSOLETE]      uint i_point=m_coord_groups[i_coord_group].m_i_icsd_g_state;
      //[CO20221112 - OBSOLETE]      if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
      //[CO20221112 - OBSOLETE]      if(m_points[i_point].m_has_entry){value=aurostd::wrapString(m_points[i_point].m_entry.auid,string_wrapper);}
      //[CO20221112 - OBSOLETE]    }
      //[CO20221112 - OBSOLETE]  }
      //[CO20221112 - OBSOLETE]}
    }
    else if(property=="prototype_icsd"){
      if(point.m_found_icsd){
        uint i_point=point.m_i_icsd;
        if(i_point>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within points");}
        value=aurostd::wrapString(m_points[i_point].m_entry.prototype,string_wrapper);
      }
    }
    else if(property=="compound_phases_equilibrium"){
      if(!(ftype==txt_ft || ftype==json_ft)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No latex rule defined for "+property);}
      if(point.isGState()&&!point.isUnary()){ //unaries are always gstates, but do NOT have any mixture context
        //need to grab from coord_group
        uint i_coord_group=AUROSTD_MAX_UINT;
        if(!getCoordGroupIndex(point,i_coord_group)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup index not set");}
        if(!m_cflags.flag("CHULL::CALCULATE_HIGHEST_DIMENSION_ONLY")){
          if(m_coord_groups[i_coord_group].m_equilibrium_phases.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Equilibrium phases not set");}
        }
        vector<string> compounds,_compounds;
        uint i_point=AUROSTD_MAX_UINT;
        const vector<vector<uint> >& equilibrium_phases=m_coord_groups[i_coord_group].m_equilibrium_phases;
        if(equilibrium_phases.size()){
          for(uint i=0,fl_size_i=equilibrium_phases.size();i<fl_size_i;i++){
            if(equilibrium_phases[i].size()==0){compounds.push_back(null_value);continue;}
            _compounds.clear();
            for(uint j=0,fl_size_j=equilibrium_phases[i].size();j<fl_size_j;j++){
              i_point=artificialMap(equilibrium_phases[i][j]);
              if(m_points[i_point].m_has_entry){
                _compounds.push_back(aurostd::wrapString(m_points[i_point].m_entry.compound,string_wrapper));
              }else{_compounds.push_back(null_value);}
            }
            compounds.push_back(aurostd::wrapString(aurostd::joinWDelimiter(_compounds,equilibrium_phases_delimiter),list_prefix,list_suffix));
          }
          value=aurostd::wrapString(aurostd::joinWDelimiter(compounds,","),list_prefix,list_suffix);
        }
      }
    }
    else if(property=="auid_phases_equilibrium"){
      if(!(ftype==txt_ft || ftype==json_ft)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No latex rule defined for "+property);}
      if(point.isGState()&&!point.isUnary()){
        //need to grab from coord_group
        uint i_coord_group=AUROSTD_MAX_UINT;
        if(!getCoordGroupIndex(point,i_coord_group)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup index not set");}
        if(!m_cflags.flag("CHULL::CALCULATE_HIGHEST_DIMENSION_ONLY")){
          if(m_coord_groups[i_coord_group].m_equilibrium_phases.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Equilibrium phases not set");}
        }
        vector<string> auids,_auids;
        uint i_point=AUROSTD_MAX_UINT;
        const vector<vector<uint> >& equilibrium_phases=m_coord_groups[i_coord_group].m_equilibrium_phases;
        if(equilibrium_phases.size()){
          for(uint i=0,fl_size_i=equilibrium_phases.size();i<fl_size_i;i++){
            if(equilibrium_phases[i].size()==0){auids.push_back(null_value);continue;}
            _auids.clear();
            for(uint j=0,fl_size_j=equilibrium_phases[i].size();j<fl_size_j;j++){
              i_point=artificialMap(equilibrium_phases[i][j]);
              if(m_points[i_point].m_has_entry){
                _auids.push_back(aurostd::wrapString(m_points[i_point].m_entry.auid,string_wrapper));
              }else{_auids.push_back(null_value);}
            }
            auids.push_back(aurostd::wrapString(aurostd::joinWDelimiter(_auids,equilibrium_phases_delimiter),list_prefix,list_suffix));
          }
          value=aurostd::wrapString(aurostd::joinWDelimiter(auids,","),list_prefix,list_suffix);
        }
      }
    }
    else if(property=="compound_phases_decomposition"){
      if(!(ftype==txt_ft || ftype==json_ft)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No latex rule defined for "+property);}
      if(GET_DECOMPOSITION_POLYMORPHS||!point.isGState()){
        //need to grab from coord_group
        uint i_coord_group=AUROSTD_MAX_UINT;
        if(!getCoordGroupIndex(point,i_coord_group)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup index not set");}
        if(GET_DECOMPOSITION_POLYMORPHS||!m_coord_groups[i_coord_group].m_is_on_hull){
          if(m_coord_groups[i_coord_group].m_decomp_phases.size()){
            vector<string> compounds;
            uint i_point=AUROSTD_MAX_UINT;
            for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_decomp_phases.size();i<fl_size_i;i++){
              i_point=artificialMap(m_coord_groups[i_coord_group].m_decomp_phases[i]);
              const xvector<double>& decomposition_coefficients=m_coord_groups[i_coord_group].m_decomp_coefs;
              //[OBSOLETE - reduce by frac_vrt always! so use coord_group values]const xvector<double>& decomposition_coefficients=point.m_decomp_coefs;
              if((decomposition_coefficients.lrows+i+1<=(uint)decomposition_coefficients.urows)&&(aurostd::nonZeroWithinTol(decomposition_coefficients[decomposition_coefficients.lrows+i+1],ZERO_TOL))){
                if(m_points[i_point].m_has_entry){
                  compounds.push_back(aurostd::wrapString(m_points[i_point].m_entry.compound,string_wrapper));
                }else{compounds.push_back(null_value);} //already did artificialMap()
              } //else{compounds.push_back(null_value);}
            }
            value=aurostd::wrapString(aurostd::joinWDelimiter(compounds,","),list_prefix,list_suffix);
          }
        }
      }
    }
    else if(property=="auid_phases_decomposition"){
      if(!(ftype==txt_ft || ftype==json_ft)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No latex rule defined for "+property);}
      if(GET_DECOMPOSITION_POLYMORPHS||!point.isGState()){
        //need to grab from coord_group
        uint i_coord_group=AUROSTD_MAX_UINT;
        if(!getCoordGroupIndex(point,i_coord_group)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup index not set");}
        if(GET_DECOMPOSITION_POLYMORPHS||!m_coord_groups[i_coord_group].m_is_on_hull){
          if(m_coord_groups[i_coord_group].m_decomp_phases.size()){
            vector<string> auids;
            uint i_point=AUROSTD_MAX_UINT;
            for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_decomp_phases.size();i<fl_size_i;i++){
              i_point=artificialMap(m_coord_groups[i_coord_group].m_decomp_phases[i]);
              const xvector<double>& decomposition_coefficients=m_coord_groups[i_coord_group].m_decomp_coefs;
              //[OBSOLETE - reduce by frac_vrt always! so use coord_group values]const xvector<double>& decomposition_coefficients=point.m_decomp_coefs;
              if((decomposition_coefficients.lrows+i+1<=(uint)decomposition_coefficients.urows)&&(aurostd::nonZeroWithinTol(decomposition_coefficients[decomposition_coefficients.lrows+i+1],ZERO_TOL))){
                if(m_points[i_point].m_has_entry){
                  auids.push_back(aurostd::wrapString(m_points[i_point].m_entry.auid,string_wrapper));
                }else{auids.push_back(null_value);} //already did artificialMap()
              } //else{auids.push_back(null_value);}
            }
            value=aurostd::wrapString(aurostd::joinWDelimiter(auids,","),list_prefix,list_suffix);
          }
        }
      }
    }
    else if(property=="coefficient_phases_decomposition"){
      if(!(ftype==txt_ft || ftype==json_ft)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No latex rule defined for "+property);}
      if(GET_DECOMPOSITION_POLYMORPHS||!point.isGState()){
        //need to grab from coord_group
        uint i_coord_group=AUROSTD_MAX_UINT;
        if(!getCoordGroupIndex(point,i_coord_group)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup index not set");}
        if(GET_DECOMPOSITION_POLYMORPHS||!m_coord_groups[i_coord_group].m_is_on_hull){
          if(m_coord_groups[i_coord_group].m_decomp_phases.size()){
            vector<double> nonzero_coefs;
            //[OBSOLETE - reduce by frac_vrt always! so use coord_group values]for(int i=point.m_decomp_coefs.lrows;i<=point.m_decomp_coefs.urows;i++){
            //[OBSOLETE - reduce by frac_vrt always! so use coord_group values]  if(aurostd::nonZeroWithinTol(point.m_decomp_coefs[i],ZERO_TOL)){nonzero_coefs.push_back(point.m_decomp_coefs[i]);}
            //[OBSOLETE - reduce by frac_vrt always! so use coord_group values]}
            for(int i=m_coord_groups[i_coord_group].m_decomp_coefs.lrows;i<=m_coord_groups[i_coord_group].m_decomp_coefs.urows;i++){
              if(aurostd::nonZeroWithinTol(m_coord_groups[i_coord_group].m_decomp_coefs[i],ZERO_TOL)){nonzero_coefs.push_back(m_coord_groups[i_coord_group].m_decomp_coefs[i]);}
            }
            precision_tmp=precision;
            tmp_roundoff_tol=5.0*pow(10,-((int)precision_tmp)-1);
            value=aurostd::wrapString(aurostd::joinWDelimiter(aurostd::vecDouble2vecString(nonzero_coefs,precision_tmp,true,tmp_roundoff_tol,FIXED_STREAM),","),list_prefix,list_suffix);
          }
        }
      }
    }
    else if(property=="nspecies"){
      uint i_nary=point.m_i_nary;
      if(i_nary!=AUROSTD_MAX_UINT){i_nary+=1;}  //m_i_nary is just index, increment by count to get nspecies
      else{i_nary=0;}  //set to something reasonable
      if(i_nary==0){    //potential other ways of defining nspecies
        if(point.m_has_entry){i_nary=point.m_entry.nspecies;}
      }
      value=aurostd::utype2string(i_nary);  //no precision needed here, simple uint
    }
    else if(property=="distance_hull_enthalpy_formation_atom"){
      d_tmp=point.getDist2Hull(_m_);
      if(d_tmp!=AUROSTD_MAX_DOUBLE){
        precision_tmp=precision;
        if(ftype==latex_ft){precision_tmp=0;tmp_roundoff_tol=5.0*pow(10,-((int)precision_tmp)-1);}
        value=aurostd::utype2string(d_tmp,precision_tmp,true,tmp_roundoff_tol,FIXED_STREAM);
      }
    }
    else if(property=="distance_hull_entropic_temperature"){
      d_tmp=point.getDist2Hull(_std_);
      if(d_tmp!=AUROSTD_MAX_DOUBLE){
        precision_tmp=precision;
        if(ftype==latex_ft){precision_tmp=0;tmp_roundoff_tol=5.0*pow(10,-((int)precision_tmp)-1);}
        value=aurostd::utype2string(d_tmp,precision_tmp,true,tmp_roundoff_tol,FIXED_STREAM); //will never be _m_ units
      }
    }
    else if(property=="stability_criterion"){
      if(m_cflags.flag("CHULL::SKIP_STABILITY_CRITERION_ANALYSIS")||point.getStabilityCriterion(_m_)>=AUROSTD_NAN){value=null_value;}
      else{
        if(point.isGState()){
          d_tmp=point.getStabilityCriterion(_m_);
          if(d_tmp!=AUROSTD_MAX_DOUBLE){
            precision_tmp=precision;
            if(ftype==latex_ft){precision_tmp=0;tmp_roundoff_tol=5.0*pow(10,-((int)precision_tmp)-1);}
            value=aurostd::utype2string(d_tmp,precision_tmp,true,tmp_roundoff_tol,FIXED_STREAM); //can be _m_ units, but smart enough to switch if T_S
          }
        }
      }
    }
    else if(property=="stability_criterion_relative"){
      if(m_cflags.flag("CHULL::SKIP_STABILITY_CRITERION_ANALYSIS")||point.getRelativeStabilityCriterion()>=AUROSTD_NAN){value=null_value;}
      else{
        if(point.isGState()){
          d_tmp=point.getRelativeStabilityCriterion();
          if(d_tmp!=AUROSTD_MAX_DOUBLE){
            precision_tmp=precision;
            if(ftype==latex_ft){precision_tmp=0;tmp_roundoff_tol=5.0*pow(10,-((int)precision_tmp)-1);}
            value=aurostd::utype2string(d_tmp,precision_tmp,true,tmp_roundoff_tol,FIXED_STREAM); //delivers as decimal, show as percentage  //CO20180409 - not showing as fraction anymore, not necessarily out of 100%
          }
        }
      }
    }
    else if(property=="N+1_enthalpy_gain"){
      if(m_cflags.flag("CHULL::SKIP_N+1_ENTHALPY_GAIN_ANALYSIS")||point.getNPlus1EnthalpyGain(_m_)>=AUROSTD_NAN){value=null_value;}
      else{
        if(point.isGState()){
          d_tmp=point.getNPlus1EnthalpyGain(_m_);
          if(d_tmp!=AUROSTD_MAX_DOUBLE){
            precision_tmp=precision;
            if(ftype==latex_ft){precision_tmp=0;tmp_roundoff_tol=5.0*pow(10,-((int)precision_tmp)-1);}
            value=aurostd::utype2string(d_tmp,precision_tmp,true,tmp_roundoff_tol,FIXED_STREAM); //delivers as decimal, show as percentage  //CO20180409 - not showing as fraction anymore, not necessarily out of 100%
          }
        }
      }
    }
    else if(property=="entropy_stabilization_coefficient"){
      d_tmp=point.getEntropyStabilizationCoefficient();
      if(d_tmp!=AUROSTD_MAX_DOUBLE){
        precision_tmp=precision;
        if(ftype==latex_ft){precision_tmp=0;tmp_roundoff_tol=5.0*pow(10,-((int)precision_tmp)-1);}
        value=aurostd::utype2string(d_tmp,precision_tmp,true,tmp_roundoff_tol,FIXED_STREAM);
      }
    }
    else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown property");}
    if(value.empty()){value=null_value;}
    if(ftype==latex_ft){
      if(external_links && (property=="compound"||(!compounds_column_report&&property=="prototype"))){
        value="\\href{"+ENTRY_PAGE_URL_PREFIX+point.m_entry.auid+"}{"+value+"}";
      }
      value=aurostd::PaddedPOST(value,padding);
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " done" << endl;}
    return value;
  }

  string ConvexHull::grabCHFacetProperty(const ChullFacet& facet,const string& property,filetype ftype) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    uint precision=COEF_PRECISION;
    double roundoff_tol=5.0*pow(10,-((int)precision)-1);
    string value="";
    string vector_delimiter=";";
    string string_wrapper="";
    string list_prefix="";
    string list_suffix="";
    string null_value="-";
    if(ftype==json_ft){
      vector_delimiter=",";
      string_wrapper="\"";
      list_prefix="[";
      list_suffix="]";
      null_value="null";
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " starting property=" << property << endl;}
    if(property=="position_vertices"){
      vector<string> vstr;
      for(uint i=0,fl_size_i=facet.m_vertices.size();i<fl_size_i;i++){
        vstr.push_back(aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(facet.m_vertices[i].ch_point.h_coords,precision,true,roundoff_tol,FIXED_STREAM),vector_delimiter));
      }
      value=aurostd::wrapString(aurostd::joinWDelimiter(aurostd::wrapVecEntries(vstr,"[","]"),","),list_prefix,list_suffix);
    }
    else if(property=="compound_vertices"){
      vector<string> compounds;
      for(uint i=0,fl_size_i=facet.m_vertices.size();i<fl_size_i;i++){
        if(facet.m_vertices[i].ch_point.m_has_entry){
          if(!m_points[facet.m_vertices[i].ch_index].m_entry.compound.empty()){
            compounds.push_back(aurostd::wrapString(m_points[facet.m_vertices[i].ch_index].m_entry.compound,string_wrapper));
          }else{compounds.push_back(null_value);}
        }
        else if(facet.m_vertices[i].ch_point.m_is_artificial){
          uint i_nary=facet.m_vertices[i].ch_point.m_i_nary;
          uint i_alloy=facet.m_vertices[i].ch_point.m_i_alloy;
          string elements=aurostd::joinWDelimiter(alloyToElements(i_nary,i_alloy),"");  //unary, so "" delimiter doesn't play a role
          compounds.push_back(aurostd::wrapString("artificial:"+elements,string_wrapper));
        }else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown print setting");}
      }
      value=aurostd::wrapString(aurostd::joinWDelimiter(compounds,","),list_prefix,list_suffix);
    }
    else if(property=="auid_vertices"){
      vector<string> auids;
      for(uint i=0,fl_size_i=facet.m_vertices.size();i<fl_size_i;i++){
        if(facet.m_vertices[i].ch_point.m_has_entry){
          if(!m_points[facet.m_vertices[i].ch_index].m_entry.auid.empty()){
            auids.push_back(aurostd::wrapString(m_points[facet.m_vertices[i].ch_index].m_entry.auid,string_wrapper));
          }else{auids.push_back(null_value);}
        }
        else if(facet.m_vertices[i].ch_point.m_is_artificial){
          uint i_nary=facet.m_vertices[i].ch_point.m_i_nary;
          uint i_alloy=facet.m_vertices[i].ch_point.m_i_alloy;
          string elements=aurostd::joinWDelimiter(alloyToElements(i_nary,i_alloy),"");  //unary, so "" delimiter doesn't play a role
          auids.push_back(aurostd::wrapString("artificial:"+elements,string_wrapper));
        }else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown print setting");}
      }
      value=aurostd::wrapString(aurostd::joinWDelimiter(auids,","),list_prefix,list_suffix);
    }
    else if(property=="content"){value=aurostd::utype2string(facet.m_content,precision,true,roundoff_tol,FIXED_STREAM);}
    else if(property=="normal"){value="["+aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(facet.m_normal,precision,true,roundoff_tol,FIXED_STREAM),vector_delimiter)+"]";}
    else if(property=="offset"){value=aurostd::utype2string(facet.m_offset,precision,true,roundoff_tol,FIXED_STREAM);}
    else if(property=="centroid"){value="["+aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(facet.m_facet_centroid,precision,true,roundoff_tol,FIXED_STREAM),vector_delimiter)+"]";}
    else if(property=="is_hypercollinear"){value=(facet.m_is_hypercollinear?"true":"false");}
    else if(property=="is_vertical"){value=(facet.m_is_vertical?"true":"false");}
    else if(property=="is_artificial"){value=(facet.m_is_artificial?"true":"false");}
    else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown property");}
    if(value.empty()){value=null_value;}
    return value;
  }

  vector<vector<string> > ConvexHull::getPointsData(const string& properties_str,vector<string>& headers,filetype ftype) const {
    stringstream message;

    vector<string> vproperties;
    uint i_point=AUROSTD_MAX_UINT;
    string value;

    vector<vector<string> > ventries;
    aurostd::string2tokens(properties_str,vproperties,",");
    for(uint i_coord_group=0,fl_size_i_coord_group=m_coord_groups.size();i_coord_group<fl_size_i_coord_group;i_coord_group++){
      if(!m_coord_groups[i_coord_group].m_points.size()){continue;}
      if(!m_coord_groups[i_coord_group].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Coordgroup["+aurostd::utype2string(i_coord_group)+"] is not initialized");}
      for(uint i=0,fl_size_i=m_coord_groups[i_coord_group].m_points.size();i<fl_size_i;i++){
        i_point=m_coord_groups[i_coord_group].m_points[i];
        if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
        const ChullPoint& point=m_points[i_point];
        if(!point.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
        if(!point.m_has_entry){continue;}
        ventries.push_back(vector<string>(0));
        for(uint j=0,fl_size_j=vproperties.size();j<fl_size_j;j++){
          value=grabCHPointProperty(point,vproperties[j],ftype);
          ventries.back().push_back(value);
        }
      }
    }

    string header;
    headers.clear();
    if(ftype==json_ft){
      //for json response, we want TRUE keywords as they would appear in aflowlib.out
      for(uint i=0,fl_size_i=vproperties.size();i<fl_size_i;i++){
        header=vproperties[i];
        //add any processing here
        headers.push_back(header);
      }
    }else{
      //for text response, we decorate nicely by replacing _atom with _[/atom]
      for(uint i=0,fl_size_i=vproperties.size();i<fl_size_i;i++){
        header=vproperties[i];
        //while this is "pretty", it's NON-STANDARD
        //define in the paper, and leave it be
        //if(vproperties[i]=="url_entry_page"){header="entry_page_url";}
        //if(vproperties[i]=="enthalpy_formation_atom"){header="formation_enthalpy";}
        //if(vproperties[i]=="spin_atom"){header="spin";}
        //if(vproperties[i]=="compound_phases_equilibrium"){header="equilibrium_phases";}
        //if(vproperties[i]=="auid_phases_equilibrium"){header="equilibrium_phases_auids";}
        //if(vproperties[i]=="compound_phases_decomposition"){header="decomposition_phases";}
        //if(vproperties[i]=="auid_phases_decomposition"){header="decomposition_phases_auids";}
        //if(vproperties[i]=="coefficient_phases_decomposition"){header="decomposition_coefficients";}
        //if(vproperties[i]=="distance_hull_enthalpy_formation_atom"){header="formation_enthalpy_difference";}
        //add any processing here before upper
        header=aurostd::toupper(header);
        //units are nice, but again, NON-STANDARD
        //define in the ppaer, and leave it be
        //add any processing here after upper
        //if(vproperties[i]=="enthalpy_formation_atom"){header+="_[meV/atom]";}
        //if(vproperties[i]=="entropic_temperature"){header+="_[K]";}
        //if(vproperties[i]=="spin_atom"){header+="_[mu_B/atom]";}
        //if(vproperties[i]=="distance_hull_enthalpy_formation_atom"){header+="_[meV/atom]";}
        //if(vproperties[i]=="distance_hull_entropic_temperature"){header+="_[K]";}
        headers.push_back(header);
      }
    }
    return ventries;
  }

  vector<vector<vector<vector<string> > > > ConvexHull::getFacetsData(const string& facet_properties_str,vector<string>& headers,filetype ftype) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    stringstream message;

    if(m_naries.size()<2){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No naries (larger than 1) found");}

    vector<string> vproperties;
    string value;
    vector<vector<vector<vector<string> > > > ventries;

    aurostd::string2tokens(facet_properties_str,vproperties,",");
    uint i_facet;
    for(uint i_nary=1,fl_size_i_nary=m_naries.size();i_nary<fl_size_i_nary;i_nary++){
      if(!m_naries[i_nary].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary["+aurostd::utype2string(i_nary)+"]");}
      if(m_naries[i_nary].m_alloys.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No alloys found in nary["+aurostd::utype2string(i_nary)+"]");}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking at i_nary=" << i_nary << endl;}
      ventries.push_back(vector<vector<vector<string> > >(0));
      for(uint i_alloy=0,fl_size_i_alloy=m_naries[i_nary].m_alloys.size();i_alloy<fl_size_i_alloy;i_alloy++){
        if(!m_naries[i_nary].m_alloys[i_alloy].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized nary[i_nary="+aurostd::utype2string(i_nary)+",i_alloy="+aurostd::utype2string(i_alloy)+"]");}
        if(!m_cflags.flag("CHULL::CALCULATE_HIGHEST_DIMENSION_ONLY")){
          if(m_naries[i_nary].m_alloys[i_alloy].m_facets.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No facets found in nary[i_nary="+aurostd::utype2string(i_nary)+",i_alloy="+aurostd::utype2string(i_alloy)+"]");}
        }
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking at i_nary=" << i_nary << ",i_alloy=" << i_alloy << endl;}
        ventries.back().push_back(vector<vector<string> >(0));
        for(uint i=0,fl_size_i=m_naries[i_nary].m_alloys[i_alloy].m_facets.size();i<fl_size_i;i++){
          i_facet=m_naries[i_nary].m_alloys[i_alloy].m_facets[i];
          const ChullFacet& facet=m_facets[i_facet];
          if(!facet.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Facet[i_nary="+aurostd::utype2string(i_nary)+",i_alloy="+aurostd::utype2string(i_alloy)+",i_facet="+aurostd::utype2string(i_facet)+"] is not initialized");}
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking at i_nary=" << i_nary << ",i_alloy=" << i_alloy << ",i_facet=" << i_facet << endl;}
          ventries.back().back().push_back(vector<string>(0));
          for(uint j=0,fl_size_j=vproperties.size();j<fl_size_j;j++){
            value=grabCHFacetProperty(facet,vproperties[j],ftype);
            ventries.back().back().back().push_back(value);
          }
        }
      }
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " created ventries" << endl;}

    string header;
    headers.clear();
    if(ftype==json_ft){
      //for json response, we want TRUE keywords as they would appear in aflowlib.out
      for(uint i=0,fl_size_i=vproperties.size();i<fl_size_i;i++){
        header=vproperties[i];
        //add any processing here
        headers.push_back(header);
      }
    }else{
      //for text response, we decorate nicely by replacing _atom with _[/atom]
      for(uint i=0,fl_size_i=vproperties.size();i<fl_size_i;i++){
        header=vproperties[i];
        //add any processing here before upper
        header=aurostd::toupper(header);
        //add any processing here after upper
        headers.push_back(header);
      }
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " created headers" << endl;}
    return ventries;
  }

  void ConvexHull::getPlainTextColumnSizes(const vector<string>& headers,const vector<vector<string> >& ventries,vector<uint>& sizes) const {
    for(uint i=0,fl_size_i=ventries.size();i<fl_size_i;i++){
      if(headers.size()!=ventries[i].size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Dimension mismatch between headers and ventries");}
    }
    if(sizes.size()==0){  //initialize
      sizes.resize(headers.size());
      for(uint i=0,fl_size_i=headers.size();i<fl_size_i;i++){sizes[i]=headers[i].size();}
    }
    if(sizes.size()!=headers.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Dimension mismatch between sizes and headers");}
    uint p_size;
    for(uint i=0,fl_size_i=ventries.size();i<fl_size_i;i++){
      for(uint j=0,fl_size_j=ventries[i].size();j<fl_size_j;j++){
        p_size=ventries[i][j].size();
        sizes[j]=std::max(sizes[j],p_size);
      }
    }
  }

  void ConvexHull::getPlainTextColumnSizesPoints(const vector<string>& headers,const vector<vector<string> >& ventries,vector<uint>& sizes) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    stringstream message;
    sizes.clear();
    getPlainTextColumnSizes(headers,ventries,sizes);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " sizes=" << sizes.size() << endl;}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " determined sizes" << endl;}
  }

  void ConvexHull::getPlainTextColumnSizesFacets(const vector<string>& headers,const vector<vector<vector<vector<string> > > >& ventries,vector<uint>& sizes) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    stringstream message;
    sizes.clear();
    for(uint i=0,fl_size_i=ventries.size();i<fl_size_i;i++){
      for(uint j=0,fl_size_j=ventries[i].size();j<fl_size_j;j++){
        getPlainTextColumnSizes(headers,ventries[i][j],sizes);
      }
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " sizes=" << sizes.size() << endl;}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " determined sizes" << endl;}
  }

  string ConvexHull::getPlainTextTable(const vector<string>& headers,const vector<vector<string> >& ventries,const vector<uint>& sizes) const {
    stringstream out_ss;
    for(uint i=0,fl_size_i=headers.size();i<fl_size_i;i++){out_ss << aurostd::PaddedPOST(headers[i],sizes[i]+5);}
    out_ss << endl;
    for(uint i=0,fl_size_i=ventries.size();i<fl_size_i;i++){
      for(uint j=0,fl_size_j=ventries[i].size();j<fl_size_j;j++){
        out_ss << aurostd::PaddedPOST(ventries[i][j],sizes[j]+5);
      }
      out_ss << endl;
    }
    return out_ss.str();
  }

  string ConvexHull::getJSONTable(const vector<string>& headers,const vector<vector<string> >& ventries) const {
    stringstream misc_ss;misc_ss.str("");
    vector<string> _vout;_vout.clear();
    vector<string> vout;vout.clear();
    for(uint i=0,fl_size_i=ventries.size();i<fl_size_i;i++){
      for(uint j=0,fl_size_j=ventries[i].size();j<fl_size_j;j++){
        misc_ss << aurostd::wrapString(headers[j],"\"") << ":" << ventries[i][j];
        _vout.push_back(misc_ss.str()); misc_ss.str("");
      }
      vout.push_back(aurostd::wrapString(aurostd::joinWDelimiter(_vout,","),"{","}")); _vout.clear();
    }
    return aurostd::wrapString(aurostd::joinWDelimiter(vout,","),"[","]");
  }

  void ConvexHull::writeText(filetype ftype) const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    stringstream message;

    if(!(ftype==txt_ft || ftype==json_ft)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown mode");}

    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Hull not initialized");}

    if(ftype==txt_ft){message << "Starting plain text generator";}
    else{message << "Starting JSON generator";}
    pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

    string properties_str_points=getPointsPropertyHeaderList(ftype);

    bool terse_output=false;
    string facet_properties_str="position_vertices";
    if(!terse_output){facet_properties_str+=",compound_vertices";}
    facet_properties_str+=",auid_vertices";
    facet_properties_str+=",content";
    facet_properties_str+=",normal";
    facet_properties_str+=",offset";
    facet_properties_str+=",centroid";
    facet_properties_str+=",is_hypercollinear";
    facet_properties_str+=",is_vertical";
    facet_properties_str+=",is_artificial";

    stringstream main_text_ss;main_text_ss.str("");
    vector<string> vout;        //json only
    vector<uint> column_sizes;  //plain text only

    //POINTS DATA
    vector<string> headers_points;
    vector<vector<string> > ventries_points=getPointsData(properties_str_points,headers_points,ftype);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " got points data" << endl;}

    //FACETS DATA
    vector<string> headers_facets;
    //first layer=nary, second=alloy, third=facet, fourth=properties
    vector<vector<vector<vector<string> > > > ventries_facets=getFacetsData(facet_properties_str,headers_facets,ftype);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " got facets data" << endl;}

    //HEADER
    if(ftype==json_ft){vout.push_back("\"hull_data\":"+getJSONHeader());}
    else{main_text_ss << getPlainTextHeader();}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " created doc header" << endl;}

    //points data
    if(ftype==json_ft){vout.push_back("\"points_data\":"+getJSONTable(headers_points,ventries_points));}
    else{
      getPlainTextColumnSizesPoints(headers_points,ventries_points,column_sizes);
      main_text_ss << endl;
      main_text_ss << "POINTS DATA" << endl;
      main_text_ss << getPlainTextTable(headers_points,ventries_points,column_sizes);
      main_text_ss << endl;
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " added points data" << endl;}

    //facets data
    if(ftype==json_ft){
      stringstream misc_ss;misc_ss.str("");
      vector<string> _vout;
      for(uint i=0,fl_size_i=ventries_facets.size();i<fl_size_i;i++){
        for(uint j=0,fl_size_j=ventries_facets[i].size();j<fl_size_j;j++){
          misc_ss << "\"" << i+2 << "-nary:";
          if(m_velements.size()){misc_ss << aurostd::joinWDelimiter(alloyToElements(i+1,j),"-");}
          else{misc_ss << j+1;}
          misc_ss << "\":" << getJSONTable(headers_facets,ventries_facets[i][j]);
          _vout.push_back(misc_ss.str()); misc_ss.str("");
        }
      }
      vout.push_back("\"facets_data\":"+aurostd::wrapString(aurostd::joinWDelimiter(_vout,","),"{","}"));
    }else{
      getPlainTextColumnSizesFacets(headers_facets,ventries_facets,column_sizes);
      main_text_ss << endl;
      for(uint i=0,fl_size_i=ventries_facets.size();i<fl_size_i;i++){
        for(uint j=0,fl_size_j=ventries_facets[i].size();j<fl_size_j;j++){
          main_text_ss << "FACETS DATA " << i+2 << "-nary ";
          if(m_velements.size()){main_text_ss << "(" << aurostd::joinWDelimiter(alloyToElements(i+1,j),"-") << ")";}
          else{main_text_ss << "(" << j+1 << ")";}
          main_text_ss << endl;
          main_text_ss << getPlainTextTable(headers_facets,ventries_facets[i][j],column_sizes);
          if(!(i==ventries_facets.size()-1 && j==ventries_facets[i].size()-1)){main_text_ss << endl;}
        }
        if(i!=ventries_facets.size()-1){main_text_ss << endl;}
      }
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " added facets data" << endl;}

    if(ftype==json_ft){main_text_ss << aurostd::wrapString(aurostd::joinWDelimiter(vout,","),"{","}");}

    string file_name="aflow_"+aurostd::joinWDelimiter(m_velements,"")+"_hull";
    if(ftype==json_ft){file_name+=".json";}
    else{file_name+=".txt";}

    if(m_cflags.flag("CHULL::SCREEN_ONLY")){
      *p_oss << main_text_ss.str();
      return;
    }

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " now writting to file" << endl;}
    string path = getPath(m_cflags, *p_FileMESSAGE, *p_oss);
    string destination = path + file_name;
    aurostd::stringstream2file(main_text_ss,destination);
    if(!aurostd::FileExist(destination)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unable to write "+file_name+" to "+path);}
    message << file_name << " was created successfully, see destination=" << path;
    pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_COMPLETE_);
  }

  void ConvexHull::writeWebApp() const {
    bool LDEBUG=(FALSE || _DEBUG_CHULL_ || XHOST.DEBUG);
    stringstream message;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Starting web-specific JSONifier", m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    // initializing stringstreams to use
    stringstream main_JSON_ss;
    stringstream id_data_JSON_ss;
    stringstream species_data_JSON_ss;
    vector<string> points_data_JSON_vs;
    vector<string> distances_data_JSON_vs;
    stringstream distances_data_JSON_ss;
    vector<string> stoich_comp_JSON_vs;
    stringstream stoich_comp_JSON_ss;
    vector<string> stoich_data_JSON_vs;
    stringstream stoich_data_JSON_ss;
    vector<string> stoich_points_data_JSON_vs;
    vector<string> vertices_data_JSON_vs;
    stringstream vertices_data_JSON_ss;
    vector<string> hull_points_data_JSON_vs;
    stringstream hull_points_data_JSON_ss;
    vector<string> planes_data_JSON_vs;
    vector<string> decomp_data_JSON_vs;
    stringstream decomp_data_JSON_ss;
    stringstream data_helper_ss;
    stringstream num_ss;  // no precision, for properties calculated by AFLOW

    // initializing some strings
    string main_JSON_file;
    string input;//,input_hyphened;

    // creating name of output file
    input=aurostd::joinWDelimiter(m_velements,"");
    //input_hyphened=aurostd::joinWDelimiter(m_velements,"-");
    main_JSON_file="aflow_"+input; //SK20200406
    //[SK20200325 - OBSOLETE]main_JSON_file="aflow_"+input+"_hull_web.json"; //WSCHMITT20190620
    //SK20200331 START
    bool fhsc_requested=m_cflags.flag("CHULL::CALCULATE_FAKE_HULL_STABILITY_CRITERION");  //only neglect feature via web
    bool fhn1eg_requested=m_cflags.flag("CHULL::CALCULATE_FAKE_HULL_N+1_ENTHALPY_GAIN");
    vector<string> sc_point;
    string delimiter="";
    // naming stability criterion files
    if(fhsc_requested) {
      aurostd::string2tokens(m_cflags.getattachedscheme("CHULL::CALCULATE_FAKE_HULL_STABILITY_CRITERION"),sc_point,",");
      std::sort(sc_point.begin(),sc_point.end()); //CO20200404 - this sort is NOT necessary, as web only removes 1 point a time, but this is SAFE
      delimiter = "_sc_";
      // limiting to the characters after "aflow:" because ":" is a reserved character for php query calls, also shortening queries
      main_JSON_file=main_JSON_file + delimiter + sc_point[0].substr(6); // restricting to single auid to limit file name growth, the substr(6) removes 'aflow:' from string which would cause issues for the filename/web
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " main_JSON_file=" << main_JSON_file << endl;}
    }
    // naming n+1 enthalpy gain files
    if (fhn1eg_requested) {
      delimiter = "_n1eg";  //CO20200404 - n+1 is the same for all points on a hull
      main_JSON_file=main_JSON_file + delimiter;
    }
    //SK20200331 END
    main_JSON_file=main_JSON_file+"_hull_web.json"; //SK20200327
    species_data_JSON_ss << aurostd::joinWDelimiter(aurostd::wrapVecEntries(m_velements,"\""),",");
    //for (uint i = 0; i < m_velements.size(); i++) {
    //  main_JSON_file.append(m_velements[i]);
    //  species_data_JSON_ss << "\"" << m_velements[i] << "\"";
    //  if(i != m_velements.size() - 1) {
    //    species_data_JSON_ss << ",";
    //  }
    //}
    //input = main_JSON_file;
    //main_JSON_file.append("_hull.json");

    // other initializations
    bool stoich_groups_set;
    xvector<double> coord;
    xvector<double> normal;

    //////////////////////////////////////////////////////////////////////////////
    // START Stoichiometry group loop
    //////////////////////////////////////////////////////////////////////////////

    uint i_point=AUROSTD_MAX_UINT,i_phase=AUROSTD_MAX_UINT;
    for(uint i=0,fl_size_i=m_coord_groups.size();i<fl_size_i;i++) {
      stoich_groups_set = false;

      ////////////////////////////////////////////////////////////////////////////
      // START Stoichiometry group points loop
      ////////////////////////////////////////////////////////////////////////////

      for(uint j=0,fl_size_j=m_coord_groups[i].m_points.size();j<fl_size_j;j++) {
        //////////////////////////////////////////////////////////////////////////
        // START Entries filter
        //////////////////////////////////////////////////////////////////////////
        i_point=m_coord_groups[i].m_points[j];
        if(!m_points[i_point].m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Point["+aurostd::utype2string(i_point)+"] is not initialized");}
        const ChullPoint& point = m_points[i_point];
        const aflowlib::_aflowlib_entry& entry = point.m_entry;

        //////////////////////////////////////////////////////////////////////////
        // END Entries filter
        //////////////////////////////////////////////////////////////////////////

        if(!stoich_groups_set) {
          ////////////////////////////////////////////////////////////////////////
          // START Stoichiometry group properties output
          ////////////////////////////////////////////////////////////////////////

          stoich_data_JSON_ss << "{";
          // get stoich id and stoichiometries at the same time
          // id_data_JSON_ss will take care of id
          for(int j=point.m_elements_present.lrows;j<=point.m_elements_present.urows;j++) {
            if(point.m_elements_present[j]==1){
              stoich_comp_JSON_ss << "{";
              stoich_comp_JSON_ss << "\"element\":\"" << m_velements[j-point.m_elements_present.lrows] << "\",";
              id_data_JSON_ss << m_velements[j-point.m_elements_present.lrows] << aurostd::utype2string(point.s_coords[j],CHULL_PRECISION); // make id out of nonzero components
              stoich_comp_JSON_ss << "\"stoichiometry\":" << aurostd::utype2string(point.s_coords[j],CHULL_PRECISION) << "}";
              stoich_comp_JSON_vs.push_back(stoich_comp_JSON_ss.str());
              stoich_comp_JSON_ss.str("");
            }
          }
          stoich_data_JSON_ss << "\"stoichiometries\":[" << aurostd::joinWDelimiter(stoich_comp_JSON_vs, ',') << "],";
          stoich_data_JSON_ss << "\"id\":\"" << id_data_JSON_ss.str() << "\",";
          id_data_JSON_ss.str("");
          // get is_hull
          stoich_data_JSON_ss << "\"isGroundState\":" << (m_coord_groups[i].m_is_on_hull?"true":"false") << ",";
          // get endPoint
          stoich_data_JSON_ss << "\"endPoint\":" << (m_coord_groups[i].m_i_nary==0?"true":"false") << ",";
          // get decomposition information for all phases
          const vector<uint>& decomposition_phases = m_coord_groups[i].m_decomp_phases;
          const xvector<double>& decomposition_coefficients = m_coord_groups[i].m_decomp_coefs;

          ////////////////////////////////////////////////////////////////////////
          // START Gathering info about decomposition reaction (equilibrium phase
          // stuff handled by javascript with FACES)
          ////////////////////////////////////////////////////////////////////////

          if(!decomposition_phases.empty() && scalar_product(decomposition_coefficients, decomposition_coefficients) >= ZERO_TOL) {
            if(decomposition_phases.size() != (uint)decomposition_coefficients.rows-1) {
              throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Size of decomposition phases != size of decomposition coefficients for coordgroup["+aurostd::utype2string(i)+"]");
            }
            for(uint k=0,fl_size_k=decomposition_phases.size();k<fl_size_k;k++) {
              if(abs(decomposition_coefficients[decomposition_coefficients.lrows+k+1]) < ZERO_TOL) {continue;}
              decomp_data_JSON_ss << "{";
              i_phase=artificialMap(decomposition_phases[k]);
              const ChullPoint& dc_phase=m_points[i_phase];
              if(!dc_phase.m_has_entry) {
                // we need to adjust for missing unaries
                // look for coords
                decomp_data_JSON_ss << "\"entry\":\"" << AFLOW_HULL_ENDPOINT_STRING << ":" << aurostd::joinWDelimiter(alloyToElements(dc_phase),"") << "\","; //unary, so "" delimiter doesn't play a role
              }else{
                const aflowlib::_aflowlib_entry& equation_entry = dc_phase.m_entry;
                decomp_data_JSON_ss << "\"entry\":\"" << equation_entry.auid << "\",";
              }
              decomp_data_JSON_ss << "\"coefficient\":" << aurostd::utype2string(decomposition_coefficients[decomposition_coefficients.lrows+k+1],CHULL_PRECISION);
              decomp_data_JSON_ss << "}";
              decomp_data_JSON_vs.push_back(decomp_data_JSON_ss.str());
              decomp_data_JSON_ss.str("");
            }
          }
          // print even if empty
          stoich_data_JSON_ss << "\"decomposition_phases\":[" << aurostd::joinWDelimiter(decomp_data_JSON_vs, ',') << "],";
          decomp_data_JSON_vs.clear();
          // filling in "points" is last
          stoich_groups_set = true;

          ////////////////////////////////////////////////////////////////////////
          // END Gathering info about decomposition reaction
          ////////////////////////////////////////////////////////////////////////

          ////////////////////////////////////////////////////////////////////////
          // END Stoichiometry group properties output
          ////////////////////////////////////////////////////////////////////////
        }

        //////////////////////////////////////////////////////////////////////////
        // START Entry properties output
        //////////////////////////////////////////////////////////////////////////

        if(point.m_has_entry) {
          data_helper_ss << "\"" << entry.auid << "\"";
          stoich_points_data_JSON_vs.push_back(data_helper_ss.str());
          data_helper_ss.str("");
          // remove artificial points from distances
          distances_data_JSON_ss << "{";
          distances_data_JSON_ss << "\"auid\":";
        }
        // get points, distances, vertices
        if(point.m_is_g_state) {
          hull_points_data_JSON_ss << "{";
          hull_points_data_JSON_ss << "\"auid\":";  // why not entry? or all auid? why two different names?
        }
        if(!point.m_has_entry) {
          data_helper_ss << "\"" << AFLOW_HULL_ENDPOINT_STRING << ":" << aurostd::joinWDelimiter(alloyToElements(point),"") << "\"";  //unary, so "" delimiter doesn't play a role
          if(point.m_is_g_state) {hull_points_data_JSON_ss << data_helper_ss.str();}
          data_helper_ss.str(""); //WSCHMITT20190731 - patching quaternary hull writer issues
        }else{
          data_helper_ss << "\"" << entry.auid << "\"";
          points_data_JSON_vs.push_back(data_helper_ss.str());
          distances_data_JSON_ss << data_helper_ss.str();
          if(point.m_is_g_state) {hull_points_data_JSON_ss << data_helper_ss.str();}
          data_helper_ss.str("");
          // wrap up distances data
          distances_data_JSON_ss << ",";
          distances_data_JSON_ss << "\"distanceToHull\":" << aurostd::utype2string(point.getDist2Hull(_std_),CHULL_PRECISION);
          // ADDED BY EGOSS
          //
          // Changed distances_data to be points data. I believe 
          // points_data_JSON may be removed and replaced with this.
          //
          distances_data_JSON_ss << ",";
          distances_data_JSON_ss << "\"compound\": \""<< entry.compound << "\"";
          distances_data_JSON_ss << ",";
          distances_data_JSON_ss << "\"composition\":[";
          // explicit dimensions
          const xvector<double>& coord = point.s_coords;
          distances_data_JSON_ss << aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(coord,CHULL_PRECISION,false),",");
          //[CO20190423 - OBSOLETE]for(int k=coord.lrows;k<=coord.urows;k++) {
          //[CO20190423 - OBSOLETE] distances_data_JSON_ss << aurostd::utype2string(coord(k),CHULL_PRECISION);  // is 3 digits okay? I normally do 15
          //[CO20190423 - OBSOLETE] if(k!=coord.urows){distances_data_JSON_ss << ",";}
          //[CO20190423 - OBSOLETE]}
          distances_data_JSON_ss << "],";
          // fix for unaries, set to 0
          if(!point.m_has_entry) {  // these are only hull_members, so they only happen to
            // endpoints
            // enthalpy of formation, row 4
            // no need for precision for next few columns, leave it same way as
            // received from AFLOW
            distances_data_JSON_ss << "\"enthalpyFormationAtom\":" << aurostd::utype2string(0.0,CHULL_PRECISION);
            distances_data_JSON_ss << ",";
            // entropic temperature, row 5
            distances_data_JSON_ss << "\"entropicTemperature\":" << aurostd::utype2string(0.0,CHULL_PRECISION);
            //WSCHMITT20190620 - adding in np1 and stab criterion to webApp - START
            distances_data_JSON_ss << ",";
            distances_data_JSON_ss << "\"nPlus1EnthalpyGain\":" << aurostd::utype2string(0.0, CHULL_PRECISION);
            distances_data_JSON_ss << ",";
            distances_data_JSON_ss << "\"stabilityCriterion\":" << aurostd::utype2string(0.0, CHULL_PRECISION);
            //WSCHMITT20190620 - adding in np1 and stab criterion to webApp - STOP
          }else{
            // enthalpy of formation, row 4
            // no need for precision for next few columns, leave it same way as
            // received from AFLOW
            num_ss << chull::H_f_atom(entry, _std_);
            distances_data_JSON_ss << "\"enthalpyFormationAtom\":" << num_ss.str();
            distances_data_JSON_ss << ",";
            num_ss.str("");
            // entropic temperature, row 5
            num_ss << chull::T_S(entry);
            distances_data_JSON_ss << "\"entropicTemperature\":" << num_ss.str();
            num_ss.str("");

            //WSCHMITT20190620 - adding in np1 and stab criterion to webApp, "N/A" for non ground-states - START
            distances_data_JSON_ss << ",";

            num_ss << ConvexHull::grabCHPointProperty(point,"N+1_enthalpy_gain",json_ft);
            distances_data_JSON_ss << "\"nPlus1EnthalpyGain\":" << num_ss.str();
            distances_data_JSON_ss << ",";
            num_ss.str("");

            num_ss << ConvexHull::grabCHPointProperty(point,"stability_criterion",json_ft);
            distances_data_JSON_ss << "\"stabilityCriterion\":" << num_ss.str();
            num_ss.str("");
            //WSCHMITT20190620 - adding in np1 and stab criterion to webApp, "N/A" for non ground-states - STOP

            //SK20200825 start
            //adding decompositionAuids for webapp
            distances_data_JSON_ss << ",";
            num_ss.str("");

            num_ss << ConvexHull::grabCHPointProperty(point,"auid_phases_decomposition",json_ft);
            distances_data_JSON_ss << "\"decompositionAuids\":" << num_ss.str();
            num_ss.str("");
            //SK20200825 end
          }
          distances_data_JSON_ss << "}";
          distances_data_JSON_vs.push_back(distances_data_JSON_ss.str());
          distances_data_JSON_ss.str("");
        }
        // distances
        // vertices
        if(point.m_is_g_state) {
          hull_points_data_JSON_ss << ",";
          hull_points_data_JSON_ss << "\"compound\": \""<< entry.compound << "\"";
          hull_points_data_JSON_ss << ",";
          hull_points_data_JSON_ss << "\"composition\":[";
          // explicit dimensions
          const xvector<double>& coord = point.s_coords;
          hull_points_data_JSON_ss << aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(coord,CHULL_PRECISION,false),",");
          //[CO20190423 - OBSOLETE]for(int k=coord.lrows;k<=coord.urows;k++) {
          //[CO20190423 - OBSOLETE]  hull_points_data_JSON_ss << aurostd::utype2string(coord(k),CHULL_PRECISION);  // is 3 digits okay? I normally do 15
          //[CO20190423 - OBSOLETE]  if(k!=coord.urows){hull_points_data_JSON_ss << ",";}
          //[CO20190423 - OBSOLETE]}
          // implicit dimension
          hull_points_data_JSON_ss << "],";
          // fix for unaries, set to 0
          if(!point.m_has_entry) {  // these are only hull_members, so they only happen to
            // endpoints
            // enthalpy of formation, row 4
            // no need for precision for next few columns, leave it same way as
            // received from AFLOW
            hull_points_data_JSON_ss << "\"enthalpyFormationAtom\":" << aurostd::utype2string(0.0,CHULL_PRECISION);
            hull_points_data_JSON_ss << ",";
            // entropic temperature, row 5
            hull_points_data_JSON_ss << "\"entropicTemperature\":" << aurostd::utype2string(0.0,CHULL_PRECISION);
            //WSCHMITT20190620 - adding in np1 and stab criterion to webApp - START
            hull_points_data_JSON_ss << ",";
            hull_points_data_JSON_ss << "\"nPlus1EnthalpyGain\":" << aurostd::utype2string(0.0, CHULL_PRECISION);
            hull_points_data_JSON_ss << ",";
            hull_points_data_JSON_ss << "\"stabilityCriterion\":" << aurostd::utype2string(0.0, CHULL_PRECISION);
            //WSCHMITT20190620 - adding in np1 and stab criterion to webApp - STOP
          }else{
            // enthalpy of formation, row 4
            // no need for precision for next few columns, leave it same way as
            // received from AFLOW
            num_ss << chull::H_f_atom(entry, _std_);
            hull_points_data_JSON_ss << "\"enthalpyFormationAtom\":" << num_ss.str();
            hull_points_data_JSON_ss << ",";
            num_ss.str("");
            // entropic temperature, row 5
            num_ss << chull::T_S(entry);
            hull_points_data_JSON_ss << "\"entropicTemperature\":" << num_ss.str();
            num_ss.str("");

            //WSCHMITT20190620 - adding in np1 and stab criterion to webApp - START
            hull_points_data_JSON_ss << ",";

            num_ss << ConvexHull::grabCHPointProperty(point,"N+1_enthalpy_gain",json_ft);
            hull_points_data_JSON_ss << "\"nPlus1EnthalpyGain\":" << num_ss.str();
            hull_points_data_JSON_ss << ",";
            num_ss.str("");

            num_ss << ConvexHull::grabCHPointProperty(point,"stability_criterion",json_ft);
            hull_points_data_JSON_ss << "\"stabilityCriterion\":" << num_ss.str();
            num_ss.str("");
            //WSCHMITT20190620 - adding in np1 and stab criterion to webApp - STOP
          }
          hull_points_data_JSON_ss << "}";
          hull_points_data_JSON_vs.push_back(hull_points_data_JSON_ss.str());
          hull_points_data_JSON_ss.str("");
        }

        //////////////////////////////////////////////////////////////////////////
        // END Entry properties output
        //////////////////////////////////////////////////////////////////////////
      }
      if(stoich_groups_set) {
        stoich_data_JSON_ss << "\"points\":[" << aurostd::joinWDelimiter(stoich_points_data_JSON_vs, ',') << "]";
        stoich_data_JSON_ss << "}";
        stoich_data_JSON_vs.push_back(stoich_data_JSON_ss.str());
        stoich_points_data_JSON_vs.clear();
        stoich_data_JSON_ss.str("");
      }

      ////////////////////////////////////////////////////////////////////////////
      // END Stoichiometry group points loop
      ////////////////////////////////////////////////////////////////////////////
    }

    //////////////////////////////////////////////////////////////////////////////
    // END Stoichiometry group loop
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // START Grabbing planes data
    //////////////////////////////////////////////////////////////////////////////
    // index specific, so we can't loop through stoichGroups

    vector<uint> vertices=getHullPoints(false); //moved up for faces
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " vertices=" << aurostd::joinWDelimiter(vertices,",") << endl;}

    vector<uint> v_ch_indices,v_vertices_indices;
    bool found=false;
    for(uint i=0,fl_size_i=m_facets.size();i<fl_size_i;i++) {
      //const xvector<double> normal = m_facets[i].m_normal;
      if(!(/*m_facets[i].m_is_hypercollinear||*/m_facets[i].m_is_vertical||m_facets[i].m_is_artificial)&&m_facets[i].m_dim==m_dim){  //keep is_hypercollinear, might create small gaps in visualization otherwise //CO20190423 - EG needs three-dimensional only for 3D
        v_ch_indices=m_facets[i].getCHIndices();
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " v_ch_indices=" << aurostd::joinWDelimiter(v_ch_indices,",") << endl;}
        v_vertices_indices.clear();
        for(uint j=0,fl_size_j=v_ch_indices.size();j<fl_size_j;j++){
          found=false;
          for(uint k=0,fl_size_k=vertices.size();k<fl_size_k&&!found;k++){
            if(v_ch_indices[j]==vertices[k]){
              v_vertices_indices.push_back(k);
              found=true;
            }
          }
          if(!found){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"ChullPoint["+aurostd::utype2string(v_ch_indices[j])+"] could not be found in getHullPoints()",_INDEX_ILLEGAL_);}
        }
        //data_helper_ss << "[" << aurostd::joinWDelimiter(m_facets[i].getCHIndices(), ',') << "]";
        data_helper_ss << "[" << aurostd::joinWDelimiter(v_vertices_indices, ',') << "]"; //CO20190423 - EG needs references to vertices, not points
        planes_data_JSON_vs.push_back(data_helper_ss.str());
        data_helper_ss.str("");
      }
    }
    //////////////////////////////////////////////////////////////////////////////
    // END Grabbing planes data
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // START Grabbing vertex data (for 3-D visualization of hull)
    //////////////////////////////////////////////////////////////////////////////
    // index specific, so we can't loop through stoichGroups

    //[CO20190423 - moving up for faces]vector<uint> vertices=getHullPoints(false);
    for(uint i=0,fl_size_i=vertices.size();i<fl_size_i;i++) {
      vertices_data_JSON_ss << "{";
      vertices_data_JSON_ss << "\"auid\":";  // why not entry? or all auid? why two different names?

      const ChullPoint& point = m_points[vertices[i]];
      const aflowlib::_aflowlib_entry& entry = point.m_entry;
      const xvector<double>& coord = point.s_coords;

      if(!point.m_has_entry) {
        vertices_data_JSON_ss << "\"" << AFLOW_HULL_ENDPOINT_STRING << ":" << aurostd::joinWDelimiter(alloyToElements(point),"") << "\""; //unary, so "" delimiter doesn't play a role
      }else{vertices_data_JSON_ss << "\"" << entry.auid << "\"";}

      vertices_data_JSON_ss << ",";
      vertices_data_JSON_ss << "\"compound\": \"" << entry.compound << "\","; 
      vertices_data_JSON_ss << "\"composition\":[";
      // explicit dimensions
      vertices_data_JSON_ss << aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(coord,CHULL_PRECISION,false),",");
      //[CO20190423 - OBSOLETE]for(int k=coord.lrows;k<=coord.urows;k++) {
      //[CO20190423 - OBSOLETE] vertices_data_JSON_ss << aurostd::utype2string(coord(k),CHULL_PRECISION); // is 3 digits okay? I normally do 15
      //[CO20190423 - OBSOLETE] if(k!=coord.rows-1){vertices_data_JSON_ss << ",";}
      //[CO20190423 - OBSOLETE]}
      vertices_data_JSON_ss << "],";
      // fix for unaries, set to 0
      if(!point.m_has_entry) {  // these are only hull_members, so they only happen to
        // endpoints
        // enthalpy of formation, row 4
        // no need for precision for next few columns, leave it same way as
        // received from AFLOW
        vertices_data_JSON_ss << "\"enthalpyFormationAtom\":" << aurostd::utype2string(0.0,CHULL_PRECISION);
        vertices_data_JSON_ss << ",";
        // entropic temperature, row 5
        vertices_data_JSON_ss << "\"entropicTemperature\":" << aurostd::utype2string(0.0,CHULL_PRECISION);
        //WSCHMITT20190620 - adding in np1 and stab criterion to webApp - START
        vertices_data_JSON_ss << ",";
        vertices_data_JSON_ss << "\"nPlus1EnthalpyGain\":" << aurostd::utype2string(0.0, CHULL_PRECISION);
        vertices_data_JSON_ss << ",";
        vertices_data_JSON_ss << "\"stabilityCriterion\":" << aurostd::utype2string(0.0, CHULL_PRECISION);
        //WSCHMITT20190620 - adding in np1 and stab criterion to webApp - STOP
      }else{
        // enthalpy of formation, row 4
        // no need for precision for next few columns, leave it same way as
        // received from AFLOW
        num_ss << chull::H_f_atom(entry, _std_);
        vertices_data_JSON_ss << "\"enthalpyFormationAtom\":" << num_ss.str();
        vertices_data_JSON_ss << ",";
        num_ss.str("");
        // entropic temperature, row 5
        num_ss << chull::T_S(entry);
        vertices_data_JSON_ss << "\"entropicTemperature\":" << num_ss.str();
        num_ss.str("");

        //WSCHMITT20190620 - adding in np1 and stab criterion to webApp, "N/A" for non ground-states - START
        vertices_data_JSON_ss << ",";

        num_ss << ConvexHull::grabCHPointProperty(point,"N+1_enthalpy_gain",json_ft);
        vertices_data_JSON_ss << "\"nPlus1EnthalpyGain\":" << num_ss.str();
        vertices_data_JSON_ss << ",";
        num_ss.str("");

        num_ss << ConvexHull::grabCHPointProperty(point,"stability_criterion",json_ft);
        vertices_data_JSON_ss << "\"stabilityCriterion\":" << num_ss.str();
        num_ss.str("");
        //WSCHMITT20190620 - adding in np1 and stab criterion to webApp, "N/A" for non ground-states - STOP
      }
      vertices_data_JSON_ss << "}";
      vertices_data_JSON_vs.push_back(vertices_data_JSON_ss.str());
      vertices_data_JSON_ss.str("");
    }
    //////////////////////////////////////////////////////////////////////////////
    // END Grabbing vertex data
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // START Output amalgamation and writing
    //////////////////////////////////////////////////////////////////////////////

    main_JSON_ss << "{";
    main_JSON_ss << "\"name\":\"" << input << "\","; //input_hyphened
    main_JSON_ss << "\"species\":[" << species_data_JSON_ss.str() << "],";
    species_data_JSON_ss.str("");
    // CHANGE BY EGOSS 
    // replaced points_data_JSON_vs with distance_data_JSON_vs
    //
    // OLD
    //main_JSON_ss << "\"points\":[" << aurostd::joinWDelimiter(points_data_JSON_vs, ',')
    //             << "],";
    //points_data_JSON_vs.clear();
    main_JSON_ss << "\"points\":[" << aurostd::joinWDelimiter(distances_data_JSON_vs, ',') << "],";
    distances_data_JSON_vs.clear();
    main_JSON_ss << "\"groundStates\":[" << aurostd::joinWDelimiter(hull_points_data_JSON_vs, ',') << "],";
    hull_points_data_JSON_vs.clear();
    main_JSON_ss << "\"vertices\":[" << aurostd::joinWDelimiter(vertices_data_JSON_vs, ',') << "],";
    vertices_data_JSON_vs.clear();
    // CHANGE BY EGOSS
    // distanceToHull is removed. Now distances are contained within the entries 
    // of points
    //
    //main_JSON_ss << "\"distanceToHull\":["
    //             << aurostd::joinWDelimiter(distances_data_JSON_vs, ',') << "],";
    //distances_data_JSON_vs.clear();

    // CHANGE BY EGOSS
    // At this time I do not require stoichiometryGroups for the visualization. 
    // It should be hidden behind a flag for now.

    if(false) {
      main_JSON_ss << "\"stoichiometryGroups\":[" << aurostd::joinWDelimiter(stoich_data_JSON_vs, ',') << "],";
      stoich_data_JSON_vs.clear();
    }

    main_JSON_ss << "\"faces\":[" << aurostd::joinWDelimiter(planes_data_JSON_vs, ',') << "]";
    planes_data_JSON_vs.clear();  // no comma
    main_JSON_ss << "}";

    if(m_cflags.flag("CHULL::SCREEN_ONLY")){
      *p_oss << main_JSON_ss.str();
      return;
    }

    string path = getPath(m_cflags, *p_FileMESSAGE, *p_oss);
    string destination = path + main_JSON_file;
    aurostd::stringstream2file(main_JSON_ss, destination);
    if(!aurostd::FileExist(destination)) {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Could not write "+main_JSON_file+" to "+path);}

    //////////////////////////////////////////////////////////////////////////////
    // END Output amalgamation and writing
    //////////////////////////////////////////////////////////////////////////////

    message << main_JSON_file << " was created successfully, see destination=" << path;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
  }

  void ConvexHull::writeAPool() const {
  }

  bool ConvexHull::sortWithinCoordGroup::operator() (uint i,uint j) {
    //ascending order
    if( (i>m_points.size()-1) || (j>m_points.size()-1) ) {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index within CoordGroup");}  //safety
    const ChullPoint& ci=m_points[i];
    const ChullPoint& cj=m_points[j];
    if(ci.isGState()!=cj.isGState()){return ci.isGState()>cj.isGState();}
    //always sort by last coord first
    bool energy_sort=(m_sort_energy_ascending ? (ci.getLastCoord()<cj.getLastCoord()) : (ci.getLastCoord()>cj.getLastCoord()));
    if(ci.m_has_entry && cj.m_has_entry){
      //if entry, then we also sort by proto, compound, and then aurl to make final print out pretty
      return ( energy_sort ||
          ((ci.getLastCoord() == cj.getLastCoord()) && 
           (ci.m_entry.prototype<cj.m_entry.prototype)) ||
          ((ci.getLastCoord() == cj.getLastCoord()) &&
           (ci.m_entry.prototype == cj.m_entry.prototype) && 
           (ci.m_entry.compound<cj.m_entry.compound)) ||
          ((ci.getLastCoord() == cj.getLastCoord()) &&
           (ci.m_entry.prototype == cj.m_entry.prototype) && 
           (ci.m_entry.compound == cj.m_entry.compound) && 
           (ci.m_entry.aurl<cj.m_entry.aurl)) );  //aurl is guaranteed to be unique (more so than auid)
    }else{return energy_sort;}
  }

  bool ConvexHull::sortCHullPoints::operator() (uint i,uint j) const {
    if(i>m_points.size()-1 || j>m_points.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid index for m_points");}
    const ChullPoint& ci=m_points[i];
    const ChullPoint& cj=m_points[j];
    sortThermoPoints stp(m_sort_stoich_ascending,m_sort_energy_ascending);
    return stp.operator()(ci,cj);
  }

  bool ConvexHull::sortFacetsByPoints::operator() (const ChullFacet& fi,const ChullFacet& fj) const {
    if(fi.m_vertices.size()!=fj.m_vertices.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Dimension mismatch among facets");} //return m_ascending_order ? fi.m_vertices.size()<rj.m_vertices.size() : fi.m_vertices.size()>rj.m_vertices.size();
    if(fi.m_normal.rows!=fj.m_normal.rows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Dimension mismatch between facet normals");}
    if(!(fi.m_initialized&&fj.m_initialized)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Uninitialized facet");}       //ensure we have inward pointing normal and angle

    if(m_auto_sort_energy){
      if(fi.m_in_lower_hemisphere!=fj.m_in_lower_hemisphere){return fi.m_in_lower_hemisphere>fj.m_in_lower_hemisphere;} //upper hemisphere last
      if(fi.m_is_hypercollinear!=fj.m_is_hypercollinear){return fi.m_is_hypercollinear>fj.m_is_hypercollinear;} //hypercollinearity last within hemisphere
      if(fi.m_is_vertical!=fj.m_is_vertical){return fi.m_is_vertical>fj.m_is_vertical;} //vertical above hypercollinearity
    }

    for(uint i=0,fl_size_i=fi.m_vertices.size();i<fl_size_i;i++){
      if(fi.m_vertices[i].ch_index!=fj.m_vertices[i].ch_index){
        const ChullPoint& ci=m_points[fi.m_vertices[i].ch_index];
        const ChullPoint& cj=m_points[fj.m_vertices[i].ch_index];
        bool sort_stoich_ascending=m_sort_stoich_ascending;
        bool sort_energy_ascending=m_sort_energy_ascending;
        if(m_auto_sort_stoich){sort_stoich_ascending=!fi.m_in_lower_hemisphere;} //left to right in lower hemisphere
        if(m_auto_sort_energy){sort_energy_ascending=fi.m_in_lower_hemisphere;}  //bottom to top in lower hemisphere
        sortThermoPoints stp(sort_stoich_ascending,sort_energy_ascending);
        return stp.operator()(ci,cj);
      }
    }

    return false;
  }
} // namespace chull

#endif  // _AFLOW_CHULL_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *           Aflow COREY OSES - Duke University 2013-2021                  *
// *                                                                         *
// ***************************************************************************
