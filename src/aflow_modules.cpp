// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *            Aflow MARCO ESTERS - Duke University 2018-2021               *
// *                                                                         *
// ***************************************************************************

// This file provides serves as an AflowIn reader for AFLOW modules (e.g. APL, AEL, AGL).
// It loads all parameters from the aflow.in file (or the defaults if no file
// is present).

#include "aflow.h"

#define _ASTROPT_APL_OLD_ string("[AFLOW_PHONONS]")
#define _ASTROPT_APL_ string("[AFLOW_APL]")
#define _ASTROPT_QHA_ string("[AFLOW_QHA]")
#define _ASTROPT_AAPL_ string("[AFLOW_AAPL]")
#define _ASTROPT_AEL_ string("[AFLOW_AEL]")
#define _ASTROPT_AGL_ string("[AFLOW_AGL]")

#define FLAG_PRECISION 7  //CO20181226

#define DEBUG_MODULES false

using aurostd::utype2string;

namespace KBIN {

  //setModules//////////////////////////////////////////////////////////////////
  // Set all module flags to their default values. This is used when no
  // aflow.in file is availabe.
  void setModules(_xvasp& xvasp) {
    _xinput xinput(xvasp);
    setModules(xinput);
    xvasp = xinput.xvasp;
  }

  // General case
  void setModules(_xinput& xinput) {
    _moduleOptions module_opts;
    module_opts.aplflags = loadDefaultsAPL();
    module_opts.aaplflags = loadDefaultsAAPL();
    module_opts.aelflags = loadDefaultsAEL();
    module_opts.aglflags = loadDefaultsAGL();  
    module_opts.qhaflags = loadDefaultsQHA(); //AS20200302
    // The readParameters functions are necessary to set xvasp
    string placeholder = "";  // acts as pseudo-aflow.in
    readParametersAPL(placeholder, module_opts, xinput);
    readParametersAAPL(placeholder, module_opts, xinput);
    readParametersQHA(placeholder, module_opts, xinput);
  }

  //readModulesfromAflowIn//////////////////////////////////////////////////////
  // Reads all module flags from an AflowIn file.
  void readModulesFromAflowIn(const string& AflowIn,
      _kflags& kflags, _xvasp& xvasp) {
    _xinput xinput(xvasp);
    readModulesFromAflowIn(AflowIn, kflags, xinput);
    xvasp = xinput.xvasp;
  }

  // General case
  void readModulesFromAflowIn(const string& AflowIn,
      _kflags& kflags, _xinput& xinput) {
    _moduleOptions module_opts;
    module_opts.aplflags = loadDefaultsAPL();
    module_opts.aaplflags = loadDefaultsAAPL();
    module_opts.aelflags = loadDefaultsAEL();
    module_opts.aglflags = loadDefaultsAGL();  
    module_opts.qhaflags = loadDefaultsQHA();
    readParametersAPL(AflowIn, module_opts, xinput);
    readParametersAAPL(AflowIn, module_opts, xinput);
    readParametersQHA(AflowIn, module_opts, xinput);
    kflags.KBIN_MODULE_OPTIONS = module_opts;
  }                 

  // APL-related functions -----------------------------------------------------

  //loadDefaultsAPL/////////////////////////////////////////////////////////////
  // Sets all APL flags to their default values.
  vector<aurostd::xoption> loadDefaultsAPL() {
    bool LDEBUG = (FALSE || XHOST.DEBUG || DEBUG_MODULES);
    vector<aurostd::xoption> aplflags;
    aurostd::xoption opt;
    opt.keyword="RELAX"; opt.option = DEFAULT_APL_RELAX; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();
    opt.keyword="RELAX_COMMENSURATE"; opt.option = DEFAULT_APL_RELAX_COMMENSURATE; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();  //ME20200427
    opt.keyword="HIBERNATE"; opt.option = DEFAULT_APL_HIBERNATE; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();
    opt.keyword="ENGINE"; opt.xscheme = DEFAULT_APL_ENGINE; aplflags.push_back(opt); opt.clear();
    opt.keyword="SUPERCELL"; opt.xscheme = ""; aplflags.push_back(opt); opt.clear();
    opt.keyword="MINATOMS"; opt.xscheme = utype2string<int>(DEFAULT_APL_MINATOMS); aplflags.push_back(opt); opt.clear();
    opt.keyword="MINATOMS_UNIFORM"; opt.xscheme = utype2string<int>(DEFAULT_APL_MINATOMS); aplflags.push_back(opt); opt.clear();
    opt.keyword="MINSHELL"; opt.xscheme = utype2string<int>(DEFAULT_APL_MINSHELL); aplflags.push_back(opt); opt.clear();
    opt.keyword="POLAR"; opt.option = DEFAULT_APL_POLAR; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();
    opt.keyword="DMAG"; opt.xscheme = utype2string<double>(DEFAULT_APL_DMAG, FLAG_PRECISION); aplflags.push_back(opt); opt.clear();
    opt.keyword="DXYZONLY"; opt.option = DEFAULT_APL_DXYZONLY; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();
    opt.keyword="DSYMMETRIZE"; opt.option = DEFAULT_APL_DSYMMETRIZE; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();
    opt.keyword="DINEQUIV_ONLY"; opt.option = DEFAULT_APL_DINEQUIV_ONLY; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear(); //CO20190131
    //[ME20181226 - now a default in .aflow.rc]// Special case: DPM can be true, false, or empty
    opt.keyword="DPM"; opt.xscheme = DEFAULT_APL_DPM; opt.option = (opt.xscheme=="ON"?true:false); aplflags.push_back(opt); opt.clear();  //CO20181226
    //[ME20181226 - now a default in .aflow.rc]// Special case: k-points options can be empty
    opt.keyword="KPPRA"; opt.xscheme = utype2string<int>(DEFAULT_PHONONS_KPPRA); aplflags.push_back(opt); opt.clear(); //CO20181226 //ME20190112
    opt.keyword="KSCHEME"; opt.xscheme = DEFAULT_PHONONS_KSCHEME; aplflags.push_back(opt); opt.clear();  //ME20190109 - KPPRA can be taken from STATIC, but KSCHEME should default to G
    opt.keyword="KPOINTS"; aplflags.push_back(opt); opt.clear();
    opt.keyword="KPOINTS_GRID"; aplflags.push_back(opt); opt.clear();  //ME20200427
    opt.keyword="KPOINTS_SHIFT"; aplflags.push_back(opt); opt.clear();  //ME20200427
    opt.keyword="PREC"; opt.xscheme = DEFAULT_APL_PREC; aplflags.push_back(opt); opt.clear();
    opt.keyword="ZEROSTATE"; opt.option = DEFAULT_APL_ZEROSTATE; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();
    opt.keyword="FREQFORMAT"; opt.xscheme = DEFAULT_APL_FREQFORMAT; aplflags.push_back(opt); opt.clear();
    opt.keyword="DC"; opt.option = DEFAULT_APL_DC; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();
    opt.keyword="DCPOINTS"; opt.xscheme = utype2string<int>(DEFAULT_APL_DCPOINTS); aplflags.push_back(opt); opt.clear();
    opt.keyword="DCPATH"; opt.xscheme = DEFAULT_APL_DCPATH; aplflags.push_back(opt); opt.clear();
    opt.keyword="DCINITCOORDSFRAC"; opt.xscheme = ""; aplflags.push_back(opt); opt.clear();
    opt.keyword="DCINITCOORDSCART"; opt.xscheme = ""; aplflags.push_back(opt); opt.clear();
    opt.keyword="DCINITCOORDSLABELS"; opt.xscheme = ""; aplflags.push_back(opt); opt.clear();
    opt.keyword="DCUSERPATH"; opt.xscheme = ""; aplflags.push_back(opt); opt.clear();
    opt.keyword="DOS"; opt.option = DEFAULT_APL_DOS; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();
    opt.keyword="DOSMETHOD"; opt.xscheme = DEFAULT_APL_DOSMETHOD; aplflags.push_back(opt); opt.clear();
    opt.keyword="DOSMESH"; opt.xscheme = DEFAULT_APL_DOSMESH; aplflags.push_back(opt); opt.clear();
    opt.keyword="DOSSMEAR"; opt.xscheme = utype2string<double>(DEFAULT_APL_DOSSMEAR, FLAG_PRECISION); aplflags.push_back(opt); opt.clear();
    opt.keyword="DOSPOINTS"; opt.xscheme = utype2string<int>(DEFAULT_APL_DOSPOINTS); aplflags.push_back(opt); opt.clear();
    opt.keyword="DOS_PROJECT"; opt.option = DEFAULT_APL_DOS_PROJECT; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();  //ME20200213
    opt.keyword="DOSPROJECTIONS_CART"; opt.xscheme = ""; aplflags.push_back(opt); opt.clear();
    opt.keyword="DOSPROJECTIONS_FRAC"; opt.xscheme = ""; aplflags.push_back(opt); opt.clear();
    opt.keyword="TP"; opt.option = DEFAULT_APL_TP; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();
    opt.keyword="DISPLACEMENTS"; opt.option = DEFAULT_APL_DISPLACEMENTS; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();
    opt.keyword="TPT"; opt.xscheme = DEFAULT_APL_TPT; aplflags.push_back(opt); opt.clear();
    opt.keyword="GROUP_VELOCITY"; opt.option = DEFAULT_APL_GVEL; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();
    //opt.keyword="ZEROSTATE_CHGCAR"; opt.option = DEFAULT_APL_ZEROSTATE_CHGCAR; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();  // OBSOLETE ME20200507
    if (LDEBUG) {
      for (uint i = 0; i < aplflags.size(); i++) {
        std::cerr << __AFLOW_FUNC__ << " key: " << aplflags[i].keyword << ", xscheme: " << aplflags[i].xscheme << ", option: " << aplflags[i].option << std::endl;
      }
    }
    return aplflags;
  }

  //writeFlagAPL///////////////////////////////////////////////////////////////
  // Determines whether flag should be written to aflow.in
  //CO20181226
  bool writeFlagAPL(const string& key,const xoption& xopt){  //ME20190113
    // return true;  OBSOLETE ME20190113
    if (xopt.isentry) {return true;}  //ME20190116 - Do not remove user entries
    if(key=="RELAX"){return true;}
    if(key=="HIBERNATE"){if((xopt.option == AFLOWRC_DEFAULT_APL_HIBERNATE) && (xopt.option == DEFAULT_APL_HIBERNATE)) {return false;}}  //ME20190113
    if(key=="ENGINE"){return true;}
    if(key=="SUPERCELL"){return true;}
    if(key=="MINATOMS"){return true;}
    if(key=="MINSHELL"){return true;}
    if(key=="POLAR"){return true;}
    if(key=="DMAG"){return true;}
    if(key=="DXYZONLY"){if((xopt.option == AFLOWRC_DEFAULT_APL_DXYZONLY) && (xopt.option == DEFAULT_APL_DXYZONLY)) {return false;}}  //ME20190113
    if(key=="DSYMMETRIZE"){if((xopt.option == AFLOWRC_DEFAULT_APL_DSYMMETRIZE) && (xopt.option == DEFAULT_APL_DSYMMETRIZE)) {return false;}}  //ME20190113
    if(key=="DINEQUIV_ONLY"){if((xopt.option == AFLOWRC_DEFAULT_APL_DINEQUIV_ONLY) && (xopt.option == DEFAULT_APL_DINEQUIV_ONLY)) {return false;}}  //ME20190113
    if(key=="DPM"){if(AFLOWRC_DEFAULT_APL_DPM==xopt.xscheme && DEFAULT_APL_DPM==xopt.xscheme){return false;}}  //ME20190113
    if(key=="PREC"){if(AFLOWRC_DEFAULT_APL_PREC==xopt.xscheme && DEFAULT_APL_PREC==xopt.xscheme){return false;}}  //ME20190113
    if(key=="ZEROSTATE"){return true;}
    if(key=="FREQFORMAT"){if(AFLOWRC_DEFAULT_APL_FREQFORMAT==xopt.xscheme && DEFAULT_APL_FREQFORMAT==xopt.xscheme){return false;}}  //ME20190113
    if(key=="DC"){return true;}
    if(key=="DCPOINTS"){if(utype2string<int>(AFLOWRC_DEFAULT_APL_DCPOINTS)==xopt.xscheme && utype2string<int>(DEFAULT_APL_DCPOINTS)==xopt.xscheme){return false;}}  //ME20190113
    if(key=="DCPATH"){if(AFLOWRC_DEFAULT_APL_DCPATH==xopt.xscheme && DEFAULT_APL_DCPATH==xopt.xscheme){return false;}}  //ME20190113
    if(key=="DOS"){return true;}
    if(key=="DOSMETHOD"){if(AFLOWRC_DEFAULT_APL_DOSMETHOD==xopt.xscheme && DEFAULT_APL_DOSMETHOD==xopt.xscheme){return false;}}  //ME20190113
    if(key=="DOSMESH"){return true;}  //ME20190113 - should always write
    if(key=="DOSSMEAR"){if(utype2string<double>(AFLOWRC_DEFAULT_APL_DOSSMEAR, FLAG_PRECISION)==xopt.xscheme && utype2string<double>(DEFAULT_APL_DOSSMEAR, FLAG_PRECISION)==xopt.xscheme){return false;}}  //ME20190113
    if(key=="DOSPOINTS"){if(utype2string<int>(AFLOWRC_DEFAULT_APL_DOSPOINTS)==xopt.xscheme && utype2string<int>(DEFAULT_APL_DOSPOINTS)==xopt.xscheme){return false;}}  //ME20190113
    if(key=="TP"){return true;}
    if(key=="TPT"){return true;}
    return true;
  }

  //readParametersAPL///////////////////////////////////////////////////////////
  // Reads APL flags from an aflow.in file.
  void readParametersAPL(const string& AflowIn,
      _moduleOptions& module_opts, _xinput& xinput) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || DEBUG_MODULES);
    string key, entry, xvaspflag;
    vector<bool> supercell_method(4, false);
    for (uint i = 0; i < module_opts.aplflags.size(); i++) {
      key = module_opts.aplflags[i].keyword;
      entry = _ASTROPT_APL_ + key + "=|" + _ASTROPT_APL_OLD_ + key + "="; //CO20181226
      entry += "|" + _ASTROPT_AAPL_ + key + "=|" + _ASTROPT_QHA_ + key + "="; //CO20181226
      if (key == "DMAG") { // for backwards compatibility
        entry += "|" + _ASTROPT_APL_ + "DISMAG"+ "=|" + _ASTROPT_APL_OLD_ + "DISMAG" + "="; //CO20181226
        entry += "|" + _ASTROPT_AAPL_ + "DISMAG" + "=|" + _ASTROPT_QHA_ + "DISMAG" + "="; //CO20181226
        entry += "|" + _ASTROPT_APL_ + "TDMAG"+ "=|" + _ASTROPT_APL_OLD_ + "TDMAG" + "="; //CO20181226
        entry += "|" + _ASTROPT_AAPL_ + "TDMAG" + "="; //CO20181226
        entry += "|" + _ASTROPT_APL_ + "TDISMAG"+ "=|" + _ASTROPT_APL_OLD_ + "TDISMAG" + "="; //CO20181226
        entry += "|" + _ASTROPT_AAPL_ + "TDISMAG" + "="; //CO20181226
      }
      else if (key == "DSYMMETRIZE") {  // for backwards compatibility
        entry += "|" + _ASTROPT_APL_ + "SYMMETRIZE=";  //CO20181226
        entry += "|" + _ASTROPT_APL_ + "SYM=";  //CO20181226
      }
      else if (key == "DINEQUIV_ONLY") {  // for backwards compatibility
        entry += "|" + _ASTROPT_APL_ + "INEQUIVONLY=";  //CO20181226
      }
      else if (key == "MINATOMS_UNIFORM") {  // for backwards compatibility
        entry += "|" + _ASTROPT_APL_ + "MINATOMS_UNIF=";        //CO20181226
        entry += "|" + _ASTROPT_APL_ + "MINATOMS_UNF=";         //CO20181226
        entry += "|" + _ASTROPT_APL_ + "MINATOMS_CUBE=";        //CO20181226
        entry += "|" + _ASTROPT_APL_ + "MINATOMS_RESTRICTED=";  //CO20181226
      }


      module_opts.aplflags[i].options2entry(AflowIn, entry, module_opts.aplflags[i].option, module_opts.aplflags[i].xscheme);

      // options2entry sets the keyword to _ASTROPT_ + key + "=", so reset
      module_opts.aplflags[i].keyword = key;

      //[ME20181226 - now a default in .aflow.rc]// Special case: set DPM to AUTO if not an entry
      //[ME20181226 - now a default in .aflow.rc]if ((key == "DPM") && (!module_opts.aplflags[i].isentry)) module_opts.aplflags[i].xscheme = "AUTO";

      // Write xvasp
      if(xinput.AFLOW_MODE_VASP) {
        xvaspflag = "AFLOWIN_FLAG::APL_" + key;
        //[CO20181226 - need to revise]if(writeFlagAPL(key,module_opts.aplflags[i].xscheme)){xinput.xvasp.aplopts.flag(xvaspflag, TRUE);}  //CO20181226
        xinput.xvasp.aplopts.flag(xvaspflag, TRUE);  //CO20181226
        xinput.xvasp.aplopts.push_attached(xvaspflag, module_opts.aplflags[i].xscheme); //this should become pop/push or changeattachedscheme (eventually)
      }

      // Supercell options
      if ((key == "SUPERCELL") && (module_opts.aplflags[i].isentry)) {
        supercell_method[0] = true;
        continue;
      }
      if ((key == "MINATOMS") && (module_opts.aplflags[i].isentry)) {
        supercell_method[1] = true;
        continue;
      }
      if ((key == "MINATOMS_UNIFORM") && (module_opts.aplflags[i].isentry)) {
        supercell_method[2] = true;
        continue;
      }
      // MAXSHELL will be skipped because it's not documented at all
      if ((key == "MINSHELL") && (module_opts.aplflags[i].isentry)) {
        supercell_method[3] = true;
        continue;
      }
    }

    // Was a supercell entry found? If not, switch to MINATOMS
    bool supercell_found = false;
    for (uint i = 0; i < supercell_method.size() && !supercell_found; i++) {
      if (supercell_method[i]) supercell_found = true;
    }
    if (!supercell_found) supercell_method[1] = true;

    // Unset contradicting/unnecessary xvasp flags
    if (xinput.AFLOW_MODE_VASP) {
      // Engine-related parameters
      if (xinput.xvasp.aplopts.getattachedscheme("AFLOWIN_FLAG::APL_ENGINE") == "LR") {
        // Unset DM parameters - do not unset DMAG because AAPL may need it
        xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_DPM", false);
        xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_DSYMMETRIZE", false); //CO20190131
        xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_DINEQUIV_ONLY", false); //CO20190131
        xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_DXYZONLY", false);
        xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_ZEROSTATE", false);
      } else {
        // DPM
        //[ME20181226 - now a default in .aflow.rc]if (xinput.xvasp.aplopts.getattachedscheme("AFLOWIN_FLAG::APL_DPM") == "AUTO") {
        //[ME20181226 - now a default in .aflow.rc]  xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_DPM", false);
        //[ME20181226 - now a default in .aflow.rc]}
      }

      // Supercell
      xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_SUPERCELL", supercell_method[0]);
      xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_MINATOMS", supercell_method[1]);
      xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_MINATOMS_UNIFORM", supercell_method[2]);
      xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_MINSHELL", supercell_method[3]);

      // Phonon dispersion
      if (xinput.xvasp.aplopts.getattachedscheme("AFLOWIN_FLAG::APL_DCPATH") == "LATTICE") {
        xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_DCINITCOORDSFRAC", false);
        xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_DCINITCOORDSCART", false);
        xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_DCINITCOORDSLABELS", false);
        xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_DCUSERPATH", false);
      }

      // DOS
      if (xinput.xvasp.aplopts.getattachedscheme("AFLOWIN_FLAG::APL_DOSMETHOD") == "LT") {
        xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_DOSSMEAR", false);
      }
    }
    if (LDEBUG) {
      for (uint i = 0; i < module_opts.aplflags.size(); i++) {
        std::cerr << __AFLOW_FUNC__ << "  " << module_opts.aplflags[i].keyword << " = " << module_opts.aplflags[i].xscheme << std::endl;
      }
    }
  }

  // AAPL-related functions ----------------------------------------------------

  //loadDefaultsAAPL////////////////////////////////////////////////////////////
  // Sets all AAPL flags to their default values.
  vector<aurostd::xoption> loadDefaultsAAPL() {
    bool LDEBUG = (FALSE || XHOST.DEBUG || DEBUG_MODULES);
    vector<aurostd::xoption> aaplflags;
    aurostd::xoption opt;
    opt.keyword="BTE"; opt.xscheme = DEFAULT_AAPL_BTE; aaplflags.push_back(opt); opt.clear();
    //opt.keyword="FOURTH_ORDER"; opt.option = DEFAULT_AAPL_FOURTH_ORDER; opt.xscheme = (opt.option?"ON":"OFF"); aaplflags.push_back(opt); opt.clear();  //ME20220420 - not working, do not use
    opt.keyword="CUT_RAD"; opt.xscheme = DEFAULT_AAPL_CUT_RAD; aaplflags.push_back(opt); opt.clear();
    opt.keyword="CUT_SHELL"; opt.xscheme = DEFAULT_AAPL_CUT_SHELL; aaplflags.push_back(opt); opt.clear();
    opt.keyword="THERMALGRID"; opt.xscheme = DEFAULT_AAPL_THERMALGRID; aaplflags.push_back(opt); opt.clear();
    opt.keyword="TCT"; opt.xscheme = DEFAULT_AAPL_TCT; aaplflags.push_back(opt); opt.clear();
    opt.keyword="SUMRULE"; opt.xscheme = utype2string<double>(DEFAULT_AAPL_SUMRULE, FLAG_PRECISION); aaplflags.push_back(opt); opt.clear();
    opt.keyword="SUMRULE_MAX_ITER"; opt.xscheme = utype2string<int>(DEFAULT_AAPL_SUMRULE_MAX_ITER); aaplflags.push_back(opt); opt.clear();
    opt.keyword="MIXING_COEFFICIENT"; opt.xscheme = utype2string<double>(DEFAULT_AAPL_MIXING_COEFFICIENT, FLAG_PRECISION); aaplflags.push_back(opt); opt.clear();
    opt.keyword="ISOTOPE"; opt.option = DEFAULT_AAPL_ISOTOPE; opt.xscheme = (opt.option?"ON":"OFF"); aaplflags.push_back(opt); opt.clear();
    opt.keyword="BOUNDARY"; opt.option = DEFAULT_AAPL_BOUNDARY; opt.xscheme = (opt.option?"ON":"OFF"); aaplflags.push_back(opt); opt.clear();
    opt.keyword="CUMULATIVEK"; opt.option = DEFAULT_AAPL_CUMULATIVEK; opt.xscheme = (opt.option?"ON":"OFF"); aaplflags.push_back(opt); opt.clear();
    opt.keyword="NANO_SIZE"; opt.xscheme = utype2string<double>(DEFAULT_AAPL_NANO_SIZE, FLAG_PRECISION); aaplflags.push_back(opt); opt.clear();
    opt.keyword="KPPRA_AAPL"; opt.xscheme = utype2string<int>(-1); aaplflags.push_back(opt);
    if (LDEBUG) {
      for (uint i = 0; i < aaplflags.size(); i++) {
        std::cerr << __AFLOW_FUNC__ << " key: " << aaplflags[i].keyword << ", xscheme: " << aaplflags[i].xscheme << ", option: " << aaplflags[i].option << std::endl;
      }
    }
    return aaplflags;
  }

  //writeFlagAAPL///////////////////////////////////////////////////////////////
  // Determines whether flag should be written to aflow.in
  //CO20181226
  bool writeFlagAAPL(const string& key,const xoption& xopt){
    // return true; OBSOLETE ME20190113
    if (xopt.isentry) {return true;}  //ME20190116 - Do not remove user entries
    if(key=="BTE"){return true;}
    if(key=="FOURTH_ORDER"){return true;}  //ME20190113 - should always write to be explicit
    if(key=="CUT_RAD"){return true;}
    if(key=="CUT_SHELL"){return true;}
    if(key=="THERMALGRID"){return true;}
    if(key=="TCT"){return true;}
    if(key=="SUMRULE"){if(utype2string<double>(AFLOWRC_DEFAULT_AAPL_SUMRULE, FLAG_PRECISION)==xopt.xscheme && utype2string<double>(DEFAULT_AAPL_SUMRULE, FLAG_PRECISION)==xopt.xscheme){return false;}}  //ME20190113
    if(key=="SUMRULE_MAX_ITER"){if(utype2string<int>(AFLOWRC_DEFAULT_AAPL_SUMRULE_MAX_ITER)==xopt.xscheme && utype2string<int>(DEFAULT_AAPL_SUMRULE_MAX_ITER)==xopt.xscheme){return false;}}  //ME20190113
    if(key=="MIXING_COEFFICIENT"){if(utype2string<double>(AFLOWRC_DEFAULT_AAPL_MIXING_COEFFICIENT, FLAG_PRECISION)==xopt.xscheme && utype2string<double>(DEFAULT_AAPL_MIXING_COEFFICIENT, FLAG_PRECISION)==xopt.xscheme){return false;}}  //ME20190113
    if(key=="ISOTOPE"){return true;}
    if(key=="BOUNDARY"){return true;}
    if(key=="CUMULATIVEK"){if((xopt.option == AFLOWRC_DEFAULT_AAPL_CUMULATIVEK) && (xopt.option == DEFAULT_AAPL_CUMULATIVEK)) {return false;}}  //ME20190113
    if(key=="NANO_SIZE"){return true;}
    return true;
  }

  //readParametersAAPL//////////////////////////////////////////////////////////
  // Reads AAPL flags from an aflow.in file.
  void readParametersAAPL(const string& AflowIn,
      _moduleOptions& module_opts, _xinput& xinput) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || DEBUG_MODULES);
    string key, entry, xvaspflag;
    for (uint i = 0; i < module_opts.aaplflags.size(); i++) {
      key = module_opts.aaplflags[i].keyword;
      entry = _ASTROPT_AAPL_ + key + "=|" + _ASTROPT_APL_OLD_ + key + "=";  //CO20181226
      module_opts.aaplflags[i].options2entry(AflowIn, entry, module_opts.aaplflags[i].option, module_opts.aaplflags[i].xscheme);
      // options2entry sets the keyword to _ASTROPT_ + key + "=", so reset
      module_opts.aaplflags[i].keyword = key;
      if (xinput.AFLOW_MODE_VASP) {
        xvaspflag = "AFLOWIN_FLAG::AAPL_" + key;
        //[CO20181226 - need to revise]if(writeFlagAAPL(key,module_opts.aaplflags[i].xscheme)){xinput.xvasp.aaplopts.flag(xvaspflag, TRUE);}  //CO20181226
        xinput.xvasp.aaplopts.flag(xvaspflag, TRUE);  //CO20181226
        xinput.xvasp.aaplopts.push_attached(xvaspflag, module_opts.aaplflags[i].xscheme); //this should become pop/push or changeattachedscheme (eventually)
      }
      // Special rules for certain keywords
      //ME20190408 START
      // If KPPRA_AAPL is not set, use APL KPPRA
      if (key == "KPPRA_AAPL" && module_opts.aaplflags[i].content_int < 1) {
        xinput.xvasp.aaplopts.flag("AFLOWIN_FLAG::AAPL_KPPRA_AAPL", false);
        continue;
      }
      //ME20190408 END
    }
    if (LDEBUG) {
      for (uint i = 0; i < module_opts.aaplflags.size(); i++) {
        std::cerr << __AFLOW_FUNC__ << "  " << module_opts.aaplflags[i].keyword << " = " << module_opts.aaplflags[i].xscheme << std::endl;
      }
    }
  }


  // AEL-related functions -----------------------------------------------------

  //loadDefaultsAEL/////////////////////////////////////////////////////////////
  // Sets all AEL flags to their default values.
  vector<aurostd::xoption> loadDefaultsAEL() {
    bool LDEBUG = (FALSE || XHOST.DEBUG || DEBUG_MODULES);
    vector<aurostd::xoption> aelflags;
    aurostd::xoption opt;
    opt.keyword="STRAIN_SYMMETRY"; opt.option = DEFAULT_AEL_STRAIN_SYMMETRY; opt.xscheme = (opt.option?"ON":"OFF"); aelflags.push_back(opt); opt.clear();
    opt.keyword="NNORMAL_STRAINS"; opt.xscheme = utype2string<int>(DEFAULT_AEL_NNORMAL_STRAINS); aelflags.push_back(opt); opt.clear();
    opt.keyword="NSHEAR_STRAINS"; opt.xscheme = utype2string<int>(DEFAULT_AEL_NNORMAL_STRAINS); aelflags.push_back(opt); opt.clear();
    opt.keyword="NORMAL_STRAIN_STEP"; opt.xscheme = utype2string<double>(DEFAULT_AEL_NORMAL_STRAIN_STEP, FLAG_PRECISION); aelflags.push_back(opt); opt.clear();
    opt.keyword="SHEAR_STRAIN_STEP"; opt.xscheme = utype2string<double>(DEFAULT_AEL_SHEAR_STRAIN_STEP, FLAG_PRECISION); aelflags.push_back(opt); opt.clear();
    opt.keyword="ORIGIN_STRAIN_CALC"; opt.option = DEFAULT_AEL_ORIGIN_STRAIN_CALC; opt.xscheme = (opt.option?"ON":"OFF"); aelflags.push_back(opt); opt.clear();
    opt.keyword="ORIGIN_STRAIN_FIT"; opt.option = DEFAULT_AEL_ORIGIN_STRAIN_FIT; opt.xscheme = (opt.option?"ON":"OFF"); aelflags.push_back(opt); opt.clear();
    opt.keyword="RELAXED_STRUCT_FIT"; opt.option = DEFAULT_AEL_RELAXED_STRUCT_FIT; opt.xscheme = (opt.option?"ON":"OFF"); aelflags.push_back(opt); opt.clear();
    opt.keyword="NEG_STRAINS"; opt.option = DEFAULT_AEL_NEG_STRAINS; opt.xscheme = (opt.option?"ON":"OFF"); aelflags.push_back(opt); opt.clear();
    opt.keyword="NIND_STRAIN_DIRS"; opt.xscheme = utype2string<int>(DEFAULT_AEL_NIND_STRAIN_DIRS); aelflags.push_back(opt); opt.clear();
    opt.keyword="VASPSYM"; opt.option = DEFAULT_AEL_VASPSYM; opt.xscheme = (opt.option?"ON":"OFF"); aelflags.push_back(opt); opt.clear();
    opt.keyword="PRECACC_ALGONORM"; opt.option = DEFAULT_AEL_PRECACC_ALGONORM; opt.xscheme = (opt.option?"ON":"OFF"); aelflags.push_back(opt); opt.clear();
    opt.keyword="VASPRUNXML_STRESS"; opt.option = DEFAULT_AEL_VASPRUNXML_STRESS; opt.xscheme = (opt.option?"ON":"OFF"); aelflags.push_back(opt); opt.clear();
    opt.keyword="AUTOSKIP_FAILED_ARUNS"; opt.option = DEFAULT_AEL_AUTOSKIP_FAILED_ARUNS; opt.xscheme = (opt.option?"ON":"OFF"); aelflags.push_back(opt); opt.clear();
    opt.keyword="SKIP_ARUNS_MAX"; opt.xscheme = utype2string<int>(DEFAULT_AEL_SKIP_ARUNS_MAX); aelflags.push_back(opt); opt.clear();
    opt.keyword="CHECK_ELASTIC_SYMMETRY"; opt.option = DEFAULT_AEL_CHECK_ELASTIC_SYMMETRY; opt.xscheme = (opt.option?"ON":"OFF"); aelflags.push_back(opt); opt.clear();
    opt.keyword="SYMMETRIZE"; opt.option = DEFAULT_AEL_SYMMETRIZE; opt.xscheme = (opt.option?"ON":"OFF"); aelflags.push_back(opt); opt.clear();
    opt.keyword="WRITE_FULL_RESULTS"; opt.option = DEFAULT_AEL_WRITE_FULL_RESULTS; opt.xscheme = (opt.option?"ON":"OFF"); aelflags.push_back(opt); opt.clear();
    opt.keyword="DIRNAME_ARUN"; opt.option = DEFAULT_AEL_DIRNAME_ARUN; opt.xscheme = (opt.option?"ON":"OFF"); aelflags.push_back(opt); opt.clear();    
    if (LDEBUG) {
      for (uint i = 0; i < aelflags.size(); i++) {
        std::cerr << __AFLOW_FUNC__ << " key: " << aelflags[i].keyword << ", xscheme: " << aelflags[i].xscheme << ", option: " << aelflags[i].option << std::endl;
      }
    }
    return aelflags;
  }

  //writeFlagAEL///////////////////////////////////////////////////////////////
  // Determines whether flag should be written to aflow.in
  bool writeFlagAEL(const string& key,const xoption& xopt){ 
    if (xopt.isentry) {return true;}  
    if(key=="STRAIN_SYMMETRY"){return true;}
    if(key=="NNORMAL_STRAINS"){return true;}  
    if(key=="NSHEAR_STRAINS"){return true;}
    if(key=="NORMAL_STRAIN_STEP"){return true;}
    if(key=="SHEAR_STRAIN_STEP"){return true;}
    if(key=="ORIGIN_STRAIN_CALC"){return true;}
    if(key=="ORIGIN_STRAIN_FIT"){return true;}
    if(key=="RELAXED_STRUCT_FIT"){return true;}
    if(key=="NEG_STRAINS"){return true;}
    if(key=="NIND_STRAIN_DIRS"){return true;}
    if(key=="VASPSYM"){return true;}
    if(key=="PRECACC_ALGONORM"){return true;}
    if(key=="VASPRUNXML_STRESS"){return true;}
    if(key=="AUTOSKIP_FAILED_ARUNS"){return true;}
    if(key=="SKIP_ARUNS_MAX"){return true;}
    if(key=="CHECK_ELASTIC_SYMMETRY"){return true;}
    if(key=="SYMMETRIZE"){if((xopt.option == AFLOWRC_DEFAULT_AEL_SYMMETRIZE) && (xopt.option == DEFAULT_AEL_SYMMETRIZE)) {return false;}}
    if(key=="WRITE_FULL_RESULTS"){return true;}
    if(key=="DIRNAME_ARUN"){return true;}
    return true;
  }

  // AGL-related functions -----------------------------------------------------

  //loadDefaultsAGL/////////////////////////////////////////////////////////////
  // Sets all AGL flags to their default values.
  vector<aurostd::xoption> loadDefaultsAGL() {
    bool LDEBUG = (FALSE || XHOST.DEBUG || DEBUG_MODULES);
    vector<aurostd::xoption> aglflags;
    aurostd::xoption opt;
    opt.keyword="AEL_POISSON_RATIO"; opt.option = DEFAULT_AGL_AEL_POISSON_RATIO; opt.xscheme = (opt.option?"ON":"OFF"); aglflags.push_back(opt); opt.clear();
    opt.keyword="NSTRUCTURES"; opt.xscheme = utype2string<int>(DEFAULT_AGL_NSTRUCTURES); aglflags.push_back(opt); opt.clear();
    opt.keyword="STRAIN_STEP"; opt.xscheme = utype2string<double>(DEFAULT_AGL_STRAIN_STEP, FLAG_PRECISION); aglflags.push_back(opt); opt.clear();
    opt.keyword="AUTOSKIP_FAILED_ARUNS"; opt.option = DEFAULT_AGL_AUTOSKIP_FAILED_ARUNS; opt.xscheme = (opt.option?"ON":"OFF"); aglflags.push_back(opt); opt.clear();
    opt.keyword="SKIP_ARUNS_MAX"; opt.xscheme = utype2string<int>(DEFAULT_AGL_SKIP_ARUNS_MAX); aglflags.push_back(opt); opt.clear();
    opt.keyword="NTEMPERATURE"; opt.xscheme = utype2string<int>(DEFAULT_AGL_NTEMPERATURE); aglflags.push_back(opt); opt.clear();
    opt.keyword="STEMPERATURE"; opt.xscheme = utype2string<double>(DEFAULT_AGL_STEMPERATURE, FLAG_PRECISION); aglflags.push_back(opt); opt.clear();
    opt.keyword="NPRESSURE"; opt.xscheme = utype2string<int>(DEFAULT_AGL_NPRESSURE); aglflags.push_back(opt); opt.clear();
    opt.keyword="SPRESSURE"; opt.xscheme = utype2string<double>(DEFAULT_AGL_SPRESSURE, FLAG_PRECISION); aglflags.push_back(opt); opt.clear();
    opt.keyword="POISSON_RATIO"; opt.xscheme = utype2string<double>(DEFAULT_AGL_POISSON_RATIO, FLAG_PRECISION); aglflags.push_back(opt); opt.clear();
    opt.keyword="IEOS"; opt.xscheme = utype2string<int>(DEFAULT_AGL_IEOS); aglflags.push_back(opt); opt.clear();
    opt.keyword="IDEBYE"; opt.xscheme = utype2string<int>(DEFAULT_AGL_IDEBYE); aglflags.push_back(opt); opt.clear();
    opt.keyword="FIT_TYPE"; opt.xscheme = utype2string<int>(DEFAULT_AGL_FIT_TYPE); aglflags.push_back(opt); opt.clear();
    opt.keyword="CHECK_EV_CONCAVITY"; opt.option = DEFAULT_AGL_CHECK_EV_CONCAVITY; opt.xscheme = (opt.option?"ON":"OFF"); aglflags.push_back(opt); opt.clear();
    opt.keyword="CHECK_EV_MIN"; opt.option = DEFAULT_AGL_CHECK_EV_MIN; opt.xscheme = (opt.option?"ON":"OFF"); aglflags.push_back(opt); opt.clear();
    opt.keyword="HUGONIOT_CALC"; opt.option = DEFAULT_AGL_HUGONIOT_CALC; opt.xscheme = (opt.option?"ON":"OFF"); aglflags.push_back(opt); opt.clear();
    opt.keyword="HUGONIOT_EXTRAPOLATE"; opt.option = DEFAULT_AGL_HUGONIOT_EXTRAPOLATE; opt.xscheme = (opt.option?"ON":"OFF"); aglflags.push_back(opt); opt.clear();
    opt.keyword="RUN_ALL_PRESSURE_TEMPERATURE"; opt.option = DEFAULT_AGL_RUN_ALL_PRESSURE_TEMPERATURE; opt.xscheme = (opt.option?"ON":"OFF"); aglflags.push_back(opt); opt.clear();
    opt.keyword="WRITE_FULL_RESULTS"; opt.option = DEFAULT_AGL_WRITE_FULL_RESULTS; opt.xscheme = (opt.option?"ON":"OFF"); aglflags.push_back(opt); opt.clear();
    opt.keyword="DIRNAME_ARUN"; opt.option = DEFAULT_AGL_DIRNAME_ARUN; opt.xscheme = (opt.option?"ON":"OFF"); aglflags.push_back(opt); opt.clear();
    opt.keyword="WRITE_GIBBS_INPUT"; opt.option = DEFAULT_AGL_WRITE_GIBBS_INPUT; opt.xscheme = (opt.option?"ON":"OFF"); aglflags.push_back(opt); opt.clear();
    opt.keyword="PLOT_RESULTS"; opt.option = DEFAULT_AGL_PLOT_RESULTS; opt.xscheme = (opt.option?"ON":"OFF"); aglflags.push_back(opt); opt.clear();

    if (LDEBUG) {
      for (uint i = 0; i < aglflags.size(); i++) {
        std::cerr << __AFLOW_FUNC__ << " key: " << aglflags[i].keyword << ", xscheme: " << aglflags[i].xscheme << ", option: " << aglflags[i].option << std::endl;
      }
    }
    return aglflags;
  }

  //writeFlagAGL///////////////////////////////////////////////////////////////
  // Determines whether flag should be written to aflow.in
  bool writeFlagAGL(const string& key,const xoption& xopt){ 
    if (xopt.isentry) {return true;}  
    if(key=="AEL_POISSON_RATIO"){return true;}
    if(key=="NSTRUCTURES"){return true;}  
    if(key=="STRAIN_STEP"){return true;}
    if(key=="AUTOSKIP_FAILED_ARUNS"){return true;}
    if(key=="SKIP_ARUNS_MAX"){return true;}
    if(key=="NTEMPERATURE"){return true;}
    if(key=="STEMPERATURE"){return true;}
    if(key=="NPRESSURE"){return true;}
    if(key=="SPRESSURE"){return true;}
    if(key=="IEOS"){if(utype2string<int>(AFLOWRC_DEFAULT_AGL_IEOS)==xopt.xscheme && utype2string<int>(DEFAULT_AGL_FIT_TYPE)==xopt.xscheme){return false;}}
    if(key=="IDEBYE"){if(utype2string<int>(AFLOWRC_DEFAULT_AGL_IDEBYE)==xopt.xscheme && utype2string<int>(DEFAULT_AGL_IDEBYE)==xopt.xscheme){return false;}}
    if(key=="FIT_TYPE"){if(utype2string<int>(AFLOWRC_DEFAULT_AGL_FIT_TYPE)==xopt.xscheme && utype2string<int>(DEFAULT_AGL_FIT_TYPE)==xopt.xscheme){return false;}}
    if(key=="CHECK_EV_CONCAVITY"){if((xopt.option == AFLOWRC_DEFAULT_AGL_CHECK_EV_CONCAVITY) && (xopt.option == DEFAULT_AGL_CHECK_EV_CONCAVITY)) {return false;}}
    if(key=="CHECK_EV_MIN"){if((xopt.option == AFLOWRC_DEFAULT_AGL_CHECK_EV_MIN) && (xopt.option == DEFAULT_AGL_CHECK_EV_MIN)) {return false;}}
    if(key=="HUGONIOT_CALC"){return true;}
    if(key=="HUGONIOT_EXTRAPOLATE"){if((xopt.option == AFLOWRC_DEFAULT_AGL_HUGONIOT_EXTRAPOLATE) && (xopt.option == DEFAULT_AGL_HUGONIOT_EXTRAPOLATE)) {return false;}}
    if(key=="RUN_ALL_PRESSURE_TEMPERATURE"){return true;}
    if(key=="WRITE_FULL_RESULTS"){return true;}
    if(key=="DIRNAME_ARUN"){return true;}
    if(key=="WRITE_GIBBS_INPUT"){if((xopt.option == AFLOWRC_DEFAULT_AGL_WRITE_GIBBS_INPUT) && (xopt.option == DEFAULT_AGL_WRITE_GIBBS_INPUT)) {return false;}}
    if(key=="PLOT_RESULTS"){return true;}

    return true;
  }

  // QHA-related functions -----------------------------------------------------
  //AS20200302
  //loadDefaultsQHA/////////////////////////////////////////////////////////////
  // Sets all QHA flags to their default values.
  vector<aurostd::xoption> loadDefaultsQHA() {
    bool LDEBUG = (FALSE || XHOST.DEBUG || DEBUG_MODULES);
    vector<aurostd::xoption> qhaflags;
    aurostd::xoption opt;
    opt.keyword="MODE"; opt.xscheme = DEFAULT_QHA_MODE; qhaflags.push_back(opt); opt.clear();
    opt.keyword="EOS"; opt.option = DEFAULT_QHA_EOS; opt.xscheme = (opt.option?"ON":"OFF"); qhaflags.push_back(opt); opt.clear();
    opt.keyword="EOS_DISTORTION_RANGE"; opt.xscheme = DEFAULT_QHA_EOS_DISTORTION_RANGE; qhaflags.push_back(opt); opt.clear();
    opt.keyword="EOS_MODEL"; opt.xscheme = DEFAULT_QHA_EOS_MODEL; qhaflags.push_back(opt); opt.clear();//AS20200818
    opt.keyword="GP_DISTORTION"; opt.xscheme = utype2string<double>(DEFAULT_QHA_GP_DISTORTION); qhaflags.push_back(opt); opt.clear();
    opt.keyword="TAYLOR_EXPANSION_ORDER"; opt.xscheme = utype2string<double>(DEFAULT_QHA_TAYLOR_EXPANSION_ORDER); qhaflags.push_back(opt); opt.clear();//AS20200602
    opt.keyword="INCLUDE_ELEC_CONTRIB"; opt.option = DEFAULT_QHA_INCLUDE_ELEC_CONTRIB; opt.xscheme = (opt.option?"ON":"OFF"); qhaflags.push_back(opt); opt.clear();
    //AS20200528
    opt.keyword="SOMMERFELD_EXPANSION"; opt.option = DEFAULT_QHA_SOMMERFELD_EXPANSION; opt.xscheme = (opt.option?"ON":"OFF"); qhaflags.push_back(opt); opt.clear();
    opt.keyword="PDIS_T"; opt.xscheme = DEFAULT_QHA_PDIS_T; qhaflags.push_back(opt); opt.clear();
    //AS20200508 BEGIN
    opt.keyword="GP_FINITE_DIFF"; opt.option = DEFAULT_QHA_GP_FINITE_DIFF; opt.xscheme = (opt.option?"ON":"OFF"); qhaflags.push_back(opt); opt.clear();
    opt.keyword="IGNORE_IMAGINARY"; opt.option = DEFAULT_QHA_IGNORE_IMAGINARY; opt.xscheme = (opt.option?"ON":"OFF"); qhaflags.push_back(opt); opt.clear();
    //AS20200508 END
    opt.keyword="RELAX_IONS_CELL"; opt.option = DEFAULT_QHA_RELAX_IONS_CELL; opt.xscheme = (opt.option?"ON":"OFF"); qhaflags.push_back(opt); opt.clear();//AS20201123

    if (LDEBUG) {
      for (uint i = 0; i < qhaflags.size(); i++) {
        std::cerr << __AFLOW_FUNC__ << " key: " << qhaflags[i].keyword << ", xscheme: " << qhaflags[i].xscheme << ", option: " << qhaflags[i].option << std::endl;
      }
    }
    return qhaflags;
  }

  //writeFlagQHA////////////////////////////////////////////////////////////////
  // Determines whether flag should be written to aflow.in
  bool writeFlagQHA(const string& key,const xoption& xopt){
    if (xopt.isentry) {return true;}
    if(key=="MODE"){return true;}
    if(key=="EOS"){return true;}
    if(key=="EOS_DISTORTION_RANGE"){return true;}
    if(key=="EOS_MODEL"){return true;}//AS20200818
    if(key=="GP_DISTORTION"){return true;}
    if(key=="TAYLOR_EXPANSION_ORDER"){return false;}//AS20200602
    if(key=="INCLUDE_ELEC_CONTRIB"){return true;}
    if(key=="SOMMERFELD_EXPANSION"){return true;}//AS20200528
    if(key=="PDIS_T"){return false;}
    //AS20200508 BEGIN
    if(key=="GP_FINITE_DIFF"){return true;}
    if(key=="IGNORE_IMAGINARY"){return false;}
    if(key=="RELAX_IONS_CELL"){return false;}
    //AS20200508 END

    return true;
  }

  //readParametersQHA///////////////////////////////////////////////////////////
  // Reads QHA flags from an aflow.in file.
  void readParametersQHA(const string& AflowIn,
      _moduleOptions& module_opts, _xinput& xinput) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || DEBUG_MODULES);
    string key="", entry="", xvaspflag="";
    for (uint i = 0; i < module_opts.qhaflags.size(); i++) {
      key = module_opts.qhaflags[i].keyword;
      entry = _ASTROPT_QHA_ + key + "=";
      module_opts.qhaflags[i].options2entry(AflowIn, entry, module_opts.qhaflags[i].option, module_opts.qhaflags[i].xscheme);

      module_opts.qhaflags[i].keyword = key;

      // Write xvasp
      if(xinput.AFLOW_MODE_VASP) {
        xvaspflag = "AFLOWIN_FLAG::QHA_" + key;
        xinput.xvasp.qhaopts.flag(xvaspflag, TRUE);
        xinput.xvasp.qhaopts.push_attached(xvaspflag, module_opts.qhaflags[i].xscheme); //this should become pop/push or changeattachedscheme (eventually)
      }
    }

    if (LDEBUG) {
      for (uint i = 0; i < module_opts.qhaflags.size(); i++) {
        std::cerr << __AFLOW_FUNC__ << "  " << module_opts.qhaflags[i].keyword << " = " << module_opts.qhaflags[i].xscheme << std::endl;
      }
    }
  }

}  // namespace KBIN

//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *            Aflow MARCO ESTERS - Duke University 2018-2021               *
// *                                                                         *
//****************************************************************************
